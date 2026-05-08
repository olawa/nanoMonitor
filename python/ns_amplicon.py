# Filename: ns_amplicon.py
# Created: 2025-11-21 21:00 CET

import edlib
import numpy as np
from collections import defaultdict, Counter
from ns_core import reverse_complement, get_streamer, detect_adapter_position, detect_internal_adapters
import ns_core
from ns_resources import ONT_ADAPTERS
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
import gc
import time

# --- UTILS ---


def get_biological_ends(read, end_length=150, cached_adapter_pos=None):
    """
    Extracts the biological start and end sequences of the insert.
    For mapped reads (BAM), uses alignment soft-clipping/start.
    For unmapped reads (FASTQ/BAM), uses adapter detection or cached position.
    """
    seq = read.query_sequence
    if not seq: return None, None, None, None
    
    # 1. Mapped Reads: Use alignment information
    if not read.is_unmapped:
        s = read.query_alignment_start
        e = read.query_alignment_end
        # Ensure we don't go out of bounds
        if s + end_length > len(seq): start_seq = seq[s:]
        else: start_seq = seq[s : s + end_length]
        
        if e - end_length < 0: end_seq = seq[:e]
        else: end_seq = seq[e - end_length : e]
        
        return start_seq, end_seq, s, e

    # 2. Unmapped Reads: Use Adapter Detection
    read_len = len(seq)
    
    # Use cached position if available (Fast Mode)
    if cached_adapter_pos is not None:
        start_idx = min(cached_adapter_pos, read_len)
        end_idx = read_len
    else:
        # Slow Mode: Align adapters for every read (only if detection failed and no fallback yet)
        # This usually only happens for the first few batches before fallback triggers
        # or if adapters are inconsistent.
        # We use the first adapter pair for detection.
        adapter_5p = ONT_ADAPTERS[0]
        
        # Align 5'
        res_5p = edlib.align(adapter_5p, seq[:300], mode="HW", task="locations")
        if res_5p["editDistance"] != -1 and res_5p["locations"]:
            start_idx = res_5p["locations"][0][1] + 1
        else:
            start_idx = 0 # Fallback to 0 if not found
        end_idx = read_len
            
    
    roi_start = seq[start_idx : min(start_idx + end_length, read_len)]
    roi_end = seq[max(0, end_idx - end_length) : end_idx]
    
    return roi_start, roi_end, start_idx, end_idx

def annotate_amplicon(chrom, start, end, gene_models):
    if not gene_models: return None
    if isinstance(gene_models, list):
        for r in gene_models:
            if r["chrom"] == chrom and max(start, r["start"]) < min(end, r["end"]):
                return r["name"]
    elif hasattr(gene_models, 'items'):  # Check if it's a dict-like object
        for name, d in gene_models.items():
            if d["chrom"] == chrom and max(start, d["start"]) < min(end, d["end"]):
                return name 
    return None

def cluster_kmers(candidates, dist_threshold=3):
    """
    Cluster similar k-mers to avoid duplicates from shifts/mismatches.
    candidates: list of (kmer_seq, count), sorted by count desc.
    Returns: clustered list of (kmer_seq, count).
    """
    if not candidates: return []
    
    clusters = [] # List of [center_seq, total_count]
    
    for seq, count in candidates:
        matched = False
        for i, (center, c_count) in enumerate(clusters):
            # Fast length check first
            if abs(len(seq) - len(center)) > dist_threshold:
                continue
                
            # Align
            res = edlib.align(seq, center, mode="NW", task="distance")
            if res["editDistance"] <= dist_threshold:
                clusters[i][1] += count
                matched = True
                break
        
        if not matched:
            clusters.append([seq, count])
            
    # Return as tuples, sorted by count
    clusters.sort(key=lambda x: x[1], reverse=True)
    return [tuple(c) for c in clusters]

def parse_amplicon_name(name):
    """
    Parse amplicon name to extract gene information.
    
    Examples:
        "chr17:41196312-41277500(BRCA1_ex2-5)" -> 
            {chrom: "chr17", start: 41196312, end: 41277500, 
             gene_name: "BRCA1_ex2-5", genes: ["BRCA1"]}
        
        "chr17:41196312-41277500(BRCA1_ex2-5,TP53_ex1)" ->
            {chrom: "chr17", start: 41196312, end: 41277500,
             gene_name: "BRCA1_ex2-5,TP53_ex1", genes: ["BRCA1", "TP53"]}
        
        "chr17:41196312-41277500" ->
            {chrom: "chr17", start: 41196312, end: 41277500,
             gene_name: None, genes: []}
        
        "PRIMER1...PRIMER2" ->
            {chrom: None, start: None, end: None,
             gene_name: None, genes: []}
    
    

    
    Returns:
        dict with keys: chrom, start, end, gene_name, genes
    """
    result = {
        "chrom": None,
        "start": None,
        "end": None,
        "gene_name": None,
        "genes": []
    }
    
    # Check if this is a sequence-based name (contains "...")
    if "..." in name:
        return result
    
    # Check if this is a coordinate-based name
    if ":" not in name or "-" not in name:
        return result
    
    try:
        # Split into coordinate part and gene part
        if "(" in name and ")" in name:
            coord_part = name.split("(")[0]
            gene_part = name.split("(")[1].rstrip(")")
            result["gene_name"] = gene_part
            
            # Extract individual gene names (before underscore or exon info)
            gene_entries = gene_part.split(",")
            for entry in gene_entries:
                # Extract base gene name (before _ex or other suffixes)
                base_gene = entry.split("_")[0] if "_" in entry else entry
                if base_gene and base_gene not in result["genes"]:
                    result["genes"].append(base_gene)
        else:
            coord_part = name
        
        # Parse coordinates - handle formats like:
        # "chr17:41196312-41277500" or "PrimerA-PrimerB:chr1:17322355-17327066"
        # Split by ":" and take the last two parts (chrom and interval)
        parts = coord_part.split(":")
        if len(parts) >= 2:
            # Last part should be "start-end", second-to-last should be chrom
            interval = parts[-1]
            chrom = parts[-2]
            
            if "-" in interval:
                start_str, end_str = interval.split("-")
                result["chrom"] = chrom
                result["start"] = int(start_str)
                result["end"] = int(end_str)
        
    except (ValueError, IndexError) as e:
        # If parsing fails, return empty result
        pass
    
    return result

# --- ALGORITHMS ---
def discover_primers_from_batch(reads, start_length=19, end_length=30, cached_adapter_pos=None):
    """
    Counts k-mer pairs (start, end) from reads to identify amplicons directly.
    For mapped reads, counts (ref_name, ref_start, ref_end) tuples.
    Returns: (pair_counts, valid_reads)
    """
    pair_counts = Counter()
    valid_reads = 0
    
    for read in reads:
        # Unified discovery: Use sequence k-mers for both Mapped and Unmapped reads
        # get_biological_ends handles the coordinate extraction for mapped reads (using alignment)
        # and adapter detection for unmapped reads.
        start_seq, end_seq, _, _ = get_biological_ends(read, end_length, cached_adapter_pos)
        
        if start_seq and end_seq and len(start_seq) >= start_length and len(end_seq) >= start_length:
            # Use first 30bp (end_length) as the signature
            s_sig = start_seq[:end_length]
            e_sig = end_seq[-end_length:]
            
            pair = (s_sig, e_sig)
            pair_counts[pair] += 1
            valid_reads += 1
            
    return pair_counts, valid_reads

def discover_primers_from_batch_parallel(reads, start_length=19, end_length=30, cached_adapter_pos=None):
    return discover_primers_from_batch(reads, start_length, end_length, cached_adapter_pos)

def init_discovery_stats():
    return defaultdict(lambda: {"count": 0, "total_accuracy": 0.0, "lengths": []})  # "total_qs": 0.0 - TODO: Add QS tracking

def match_kmer_pairs_batch(reads, top_pairs, end_length=150, max_edit_dist=3, cached_adapter_pos=None):
    """
    Match reads against identified (start, end) pairs.
    top_pairs: list of (start_seq, end_seq) tuples.
    Returns: (matches_found, local_stats, read_amplicon_map)
    """
    local_stats = defaultdict(lambda: {"count": 0, "total_accuracy": 0.0, "lengths": []})  # "total_qs": 0.0 - TODO: Add QS tracking
    read_amplicon_map = {}
    matches_found = 0
    
    # Pre-compute reverse complements if needed? 
    # No, top_pairs are extracted directly from reads, so they match the read orientation.
    
    for read in reads:
        roi_start, roi_end, _, _ = get_biological_ends(read, end_length, cached_adapter_pos)
        if not roi_start or not roi_end: continue 
        
        read_len = read.query_length
        try: acc = (1 - read.get_tag("de")) * 100.0
        except KeyError: acc = 0.0
        
        best_pair = None
        min_dist = float('inf')
        
        # Check against top pairs
        # Optimization: Check exact match first?
        # Or just use edlib for robustness.
        
        # Determine if we are matching positions or sequences
        if top_pairs and len(top_pairs[0]) == 3 and isinstance(top_pairs[0][1], int):
            # Position Matching (Mapped)
            if read.is_unmapped: continue
            
            r_chrom = read.reference_name
            r_start = read.reference_start
            r_end = read.reference_end
            
            for p_chrom, p_start, p_end in top_pairs:
                if r_chrom == p_chrom:
                    # Check proximity (e.g. 50bp)
                    if abs(r_start - p_start) < 50 and abs(r_end - p_end) < 50:
                        best_pair = (p_chrom, p_start, p_end)
                        # Name resolution (Gene Name or Coords)
                        # We need to use the same name as in discovery!
                        # In discovery we used candidates list to send names to UI.
                        # But here we need to key the stats.
                        # Let's use the coord string as key, or lookup?
                        # For simplicity, use coord string.
                        break # Greedy match
        else:
            # Sequence Matching (Unmapped)
            for p_start, p_end in top_pairs:
                # We match the signature (30bp) against the ROI (150bp)
                # Or just match signature against signature?
                # If we used 30bp for discovery, p_start is 30bp.
                # roi_start is 150bp.
                
                res_s = edlib.align(p_start, roi_start, mode="HW", task="distance", k=max_edit_dist)
                if res_s["editDistance"] == -1: continue
                
                res_e = edlib.align(p_end, roi_end, mode="HW", task="distance", k=max_edit_dist)
                if res_e["editDistance"] == -1: continue
                
                total_dist = res_s["editDistance"] + res_e["editDistance"]
                if total_dist < min_dist:
                    min_dist = total_dist
                    best_pair = (p_start, p_end)
        
        if best_pair:
            matches_found += 1
            # Use the pair signature as the name
            if len(best_pair) == 3 and isinstance(best_pair[1], int):
                name = f"{best_pair[0]}:{best_pair[1]}-{best_pair[2]}"
                # Try to resolve gene name if possible?
                # Ideally we pass the gene map to this function.
                # But for now, let's stick to coords.
                # Wait, if we used Gene Name in discovery, we should use it here too.
                # But we don't have the gene map here easily.
                # Let's rely on finalize_discovery_stats to rename?
                # Or just use coords.
            else:
                name = f"{best_pair[0]}...{best_pair[1]}"
                
            # TODO: Calculate QS
            # qs = 0.0
            # quals = read.query_qualities
            # if quals: qs = float(np.mean(quals))

            s = local_stats[name]
            s["count"] += 1
            s["total_accuracy"] += acc
            # s["total_qs"] += qs  # TODO: Add QS tracking
            s["lengths"].append(read_len)
            
            # Track read mapping
            read_amplicon_map[read.query_name] = name
            
    return matches_found, local_stats, read_amplicon_map

def cluster_discovered_amplicons(amplicon_stats, read_amplicon_map, similarity_threshold=0.90, max_edit_dist=5):
    """
    Cluster similar amplicons by comparing their primer sequences.
    Merges variants that are highly similar (>90% identity).
    
    Args:
        amplicon_stats: Dict of {amplicon_name: stats}
        read_amplicon_map: Dict of {read_id: amplicon_name}
        similarity_threshold: Minimum similarity to merge (default 0.90)
        max_edit_dist: Maximum edit distance for edlib alignment (default 5)
    
    Returns:
        tuple: (clustered_stats, updated_read_amplicon_map)
    """
    # Sort by count (most abundant first)
    sorted_amps = sorted(amplicon_stats.items(), key=lambda x: x[1]['count'], reverse=True)
    
    clustered_stats = {}
    amplicon_mapping = {}  # Maps old name -> new representative name
    used = set()
    
    for amp_name, amp_data in sorted_amps:
        if amp_name in used:
            continue
        
        # Extract primer sequences from amplicon name
        # Format: "PRIMER1...PRIMER2" or "chr:start-end"
        if "..." in amp_name:
            # Sequence-based name
            parts = amp_name.split("...")
            if len(parts) == 2:
                fwd_primer = parts[0]
                rev_primer = parts[1]
            else:
                # Can't cluster this one
                clustered_stats[amp_name] = amp_data
                amplicon_mapping[amp_name] = amp_name
                used.add(amp_name)
                continue
        else:
            # Position-based name (chr:start-end) - can't cluster by sequence
            clustered_stats[amp_name] = amp_data
            amplicon_mapping[amp_name] = amp_name
            used.add(amp_name)
            continue
        
        # Find similar amplicons
        cluster_members = [amp_name]
        
        for other_name, other_data in sorted_amps:
            if other_name in used or other_name == amp_name:
                continue
            
            if "..." not in other_name:
                continue
            
            other_parts = other_name.split("...")
            if len(other_parts) != 2:
                continue
            
            other_fwd = other_parts[0]
            other_rev = other_parts[1]
            
            # Compare primers using edlib
            fwd_match = edlib.align(fwd_primer, other_fwd, mode="NW", task="distance", k=max_edit_dist)
            rev_match = edlib.align(rev_primer, other_rev, mode="NW", task="distance", k=max_edit_dist)
            
            if fwd_match["editDistance"] == -1 or rev_match["editDistance"] == -1:
                continue
            
            # Calculate similarity
            fwd_len = max(len(fwd_primer), len(other_fwd))
            rev_len = max(len(rev_primer), len(other_rev))
            
            fwd_sim = 1.0 - (fwd_match["editDistance"] / fwd_len) if fwd_len > 0 else 0
            rev_sim = 1.0 - (rev_match["editDistance"] / rev_len) if rev_len > 0 else 0
            
            avg_sim = (fwd_sim + rev_sim) / 2.0
            
            if avg_sim >= similarity_threshold:
                cluster_members.append(other_name)
                amplicon_mapping[other_name] = amp_name
                used.add(other_name)
        
        # Merge cluster stats
        merged_data = {
            "count": sum(amplicon_stats[m]["count"] for m in cluster_members),
            "total_accuracy": sum(amplicon_stats[m]["total_accuracy"] for m in cluster_members),
            "lengths": []
        }
        
        for member in cluster_members:
            merged_data["lengths"].extend(amplicon_stats[member]["lengths"])
        
        clustered_stats[amp_name] = merged_data
        amplicon_mapping[amp_name] = amp_name
        used.add(amp_name)
    
    # Update read_amplicon_map with new clustered names
    updated_map = {}
    for read_id, old_amp_name in read_amplicon_map.items():
        # Strip |CONCAT suffix if present
        clean_name = old_amp_name.replace("|CONCAT", "")
        concat_suffix = "|CONCAT" if "|CONCAT" in old_amp_name else ""
        
        new_amp_name = amplicon_mapping.get(clean_name, clean_name)
        updated_map[read_id] = new_amp_name + concat_suffix
    
    return clustered_stats, updated_map

def finalize_discovery_stats(amplicon_stats, gene_models=None, known_names=None):
    """
    Finalizes stats: calculates mean accuracy, converts lists to arrays.
    Also maps coord-based names to Gene Names if gene_models provided.
    known_names: Optional dict mapping 'chrom:start-end' -> 'GeneName' to avoid re-lookup.
    """
    final_stats = {}
    
    for name, data in amplicon_stats.items():
        count = data["count"]
        if count == 0: continue
        
        mean_acc = data["total_accuracy"] / count
        lengths = np.array(data["lengths"], dtype=np.int32)
        
        # Use known mapping if available, otherwise keep original name
        # (Discovery mode resolves names at the end and passes them here)
        final_name = name
        if known_names and name in known_names:
            final_name = known_names[name]
        
        # Parse gene information from the final name
        gene_info = parse_amplicon_name(final_name)
            
        final_stats[final_name] = {
            "count": count,
            "average_accuracy": mean_acc,
            # "average_qs": data["total_qs"] / count if count > 0 else 0.0,  # TODO: Add QS tracking
            "median_length": np.median(lengths),
            "stdev_length": np.std(lengths),
            "raw_lengths": lengths,
            "positions": [], # We don't store per-read positions anymore
            # Gene metadata for programmatic access
            "chrom": gene_info["chrom"],
            "start": gene_info["start"],
            "end": gene_info["end"],
            "region": f"{gene_info['chrom']}:{gene_info['start']}-{gene_info['end']}" if gene_info["chrom"] and gene_info["start"] is not None and gene_info["end"] is not None else None,
            "gene_name": gene_info["gene_name"],
            "genes": gene_info["genes"]
        }
        
    return final_stats


def identify_primers_combinatorial(read, primer_list, gene_models, max_edit_dist=3, end_length=100, cached_adapter_pos=None, primer_tolerance=0):
    """
    Identify primers for a single read using combinatorial search (Start/End vs All Primers).
    Returns: (amplicon_name, length, accuracy, read_id) or None
    """
    roi_start_seq, roi_end_seq, start_pos, end_pos = get_biological_ends(read, end_length=end_length, cached_adapter_pos=cached_adapter_pos)
    if roi_start_seq is None or roi_end_seq is None:
        return None
        
    best_start = None
    best_end = None
    min_dist_start = max_edit_dist + 1
    min_dist_end = max_edit_dist + 1
    
    # Iterate all primers to find best match at Start and End
    # We check both Primer and Primer_RC against the read ends
    for name, seq in primer_list:
        # Check Start
        res_s = edlib.align(seq, roi_start_seq, mode="HW", task="distance", k=max_edit_dist)
        if res_s['editDistance'] != -1 and res_s['editDistance'] < min_dist_start:
            min_dist_start = res_s['editDistance']
            best_start = name

        # Check End (Expect Reverse Complement of primer if it's the other end of the amplicon)
        # Note: If the read is Fwd strand, end should match RevPrimer_RC. 
        # But since we scan ALL primers, we just check RC of everything against the end.
        seq_rc = reverse_complement(seq)
        res_e = edlib.align(seq_rc, roi_end_seq, mode="HW", task="distance", k=max_edit_dist)
        if res_e['editDistance'] != -1 and res_e['editDistance'] < min_dist_end:
            min_dist_end = res_e['editDistance']
            best_end = name
            
    if best_start and best_end:
        # Found a pair!
        p1 = best_start
        p2 = best_end
        
        # Sort to ensure P1-P2 is same as P2-P1
        pair_names = sorted([p1, p2])
        base_name = f"{pair_names[0]}-{pair_names[1]}"
        
        # Annotate with Gene Info if mapped
        gene_info = ""
        coord_info = ""
        if not read.is_unmapped:
             # Use alignment coordinates
             chrom = read.reference_name
             # Use mapped start/end (approximate amplicon bounds)
             start = read.reference_start
             end = read.reference_end
             
             coord_info = f":{chrom}:{start}-{end}"
             
             if gene_models:
                 # Annotate
                 g_name = annotate_amplicon(chrom, start, end, gene_models)
                 if g_name:
                     gene_info = f"({g_name})"
                 
        amplicon_name = f"{base_name}{coord_info}{gene_info}"
        
        length = read.query_length
        try:
            acc = (1 - read.get_tag("de")) * 100.0
        except KeyError:
            acc = 0.0
            
        return (amplicon_name, length, acc, read.query_name)

    return None

def process_amplicon_batch(reads, primer_list, gene_models, amplicon_stats, max_edit_dist=3, end_length=100, primer_tolerance=0):
    """
    Process a batch of reads sequentially (for non-parallel mode).
    """
    for read in reads:
        result = identify_primers_combinatorial(read, primer_list, gene_models, max_edit_dist, end_length, primer_tolerance=primer_tolerance)
        if result:
            name, length, acc, _ = result  # , qs - TODO: Add QS tracking
            amplicon_stats[name]["count"] += 1
            amplicon_stats[name]["total_accuracy"] += acc
            # amplicon_stats[name]["total_qs"] += qs  # TODO: Add QS tracking
            amplicon_stats[name]["lengths"].append(length)

def process_amplicon_batch_parallel(reads, primer_list, gene_models, max_edit_dist=4, end_length=120, primer_tolerance=0):
    """
    Process a batch of reads in parallel (returns results, no shared state).
    """
    results = []
    for read in reads:
        result = identify_primers_combinatorial(read, primer_list, gene_models, max_edit_dist, end_length, primer_tolerance=primer_tolerance)
        if result:
            results.append(result)
    return results

def _run_parallel_primer_mode(streamer, threads, stop_check_cb, partial_cb, progress_cb,
                               adapter_kmers, kmer_size, primer_list, gene_models, 
                               stats_lock, amplicon_stats, read_amplicon_map, primer_tolerance=0): # Added primer_tolerance # Added primer_tolerance
    """Helper for parallel primer-based analysis."""
    internal_adapter_count = 0
    last_update_time = time.time()
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        
        for batch, meta_batch in streamer.stream_batches(progress_cb):
            if stop_check_cb and stop_check_cb():
                print("DEBUG: Stop requested. Aborting analysis.")
                break
            
            # Live updates
            current_time = time.time()
            if partial_cb:
                payload = {
                    "metadata": meta_batch,
                    "read_amplicon_map": dict(read_amplicon_map)
                }
                
                if current_time - last_update_time > 50.0:
                    with stats_lock:
                        formatted_amps = {}
                        for name, d in amplicon_stats.items():
                            c = d["count"]
                            if c > 0:
                                formatted_amps[name] = {
                                    "count": c, 
                                    "average_accuracy": d["total_accuracy"]/c,
                                    "median_length": 0, # Placeholder for speed
                                    "stdev_length": 0,
                                    "raw_lengths": [] # Don't send full list for live updates
                                }
                                # Optional: Calculate median if needed for UI
                                if d["lengths"]:
                                    formatted_amps[name]["median_length"] = d["lengths"][-1] # Use last length as proxy or calc median
                                    # actually np.median is fine
                                    # formatted_amps[name]["median_length"] = np.median(d["lengths"])
                        payload["amplicons"] = formatted_amps
                    last_update_time = current_time
                
                partial_cb(payload)
            
            # Detect concatemers in this batch (fast k-mer check)
            batch_concatemer_count, concatemer_ids = detect_internal_adapters(batch, adapter_kmers, kmer_size)
            internal_adapter_count += batch_concatemer_count
            
            # Mark concatemer reads in the mapping (use special marker)
            for read_id in concatemer_ids:
                read_amplicon_map[read_id] = read_amplicon_map.get(read_id, "Unknown") + "|CONCAT"
            
            # Submit batch to worker pool
            future = executor.submit(process_amplicon_batch_parallel, batch, primer_list, gene_models, 3, 150, primer_tolerance)
            futures.append(future)
            
            # Process completed futures (limit queue size)
            if len(futures) > threads * 2:
                from concurrent.futures import wait, FIRST_COMPLETED
                done, _ = wait(futures, return_when=FIRST_COMPLETED)
                for f in done:
                    results = f.result()
                    # Aggregate results (thread-safe)
                    with stats_lock:
                        for name, length, acc, read_id in results:
                            amplicon_stats[name]["count"] += 1
                            amplicon_stats[name]["total_accuracy"] += acc
                            amplicon_stats[name]["lengths"].append(length)
                            # Track read-to-amplicon mapping
                            read_amplicon_map[read_id] = name
                    futures.remove(f)
        
        # Process remaining futures
        for f in as_completed(futures):
            results = f.result()
            with stats_lock:
                for name, length, acc, read_id in results:
                    amplicon_stats[name]["count"] += 1
                    amplicon_stats[name]["total_accuracy"] += acc
                    amplicon_stats[name]["lengths"].append(length)
                    # Track read-to-amplicon mapping
                    read_amplicon_map[read_id] = name
                    
    return internal_adapter_count

def _run_parallel_discovery_mode(streamer, threads, stop_check_cb, partial_cb, progress_cb,
                                 adapter_kmers, kmer_size, gene_models,
                                 stats_lock, read_amplicon_map, debug_log, result_payload):
    """Helper for parallel discovery mode analysis."""
    internal_adapter_count = 0
    total_valid_reads_for_discovery = 0
    global_kmer_counts = Counter()
    discovery_stats = init_discovery_stats()
    
    DISCOVERY_LIMIT = 50000
    buffer = []
    discovery_done = False
    top_kmers = []
    cached_adapter_pos = None
    total_matches = 0
    last_update_time = time.time()
    adapter_detection_attempted = False
    coord_to_name_map = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        
        for batch, meta_batch in streamer.stream_batches(progress_cb):
            if stop_check_cb and stop_check_cb():
                print("DEBUG: Stop requested. Aborting analysis.")
                break
            
            # Live updates
            current_time = time.time()
            if partial_cb:
                payload = {"metadata": meta_batch}
                if current_time - last_update_time > 50.0:
                    # Send amplicon stats update
                    with stats_lock:
                        snapshot = {k: v.copy() for k, v in discovery_stats.items()}
                        for k in snapshot:
                            snapshot[k]["lengths"] = list(snapshot[k]["lengths"])
                    
                    final_snapshot = finalize_discovery_stats(snapshot, gene_models)
                    payload["amplicons"] = final_snapshot
                    last_update_time = current_time
                
                partial_cb(payload)
            
            # Detect concatemers
            batch_concatemer_count, concatemer_ids = detect_internal_adapters(batch, adapter_kmers, kmer_size)
            internal_adapter_count += batch_concatemer_count
            
            for read_id in concatemer_ids:
                read_amplicon_map[read_id] = read_amplicon_map.get(read_id, "Unknown") + "|CONCAT"
            
            if not discovery_done:
                # Phase 1: Buffering & Discovery
                # For mapped reads (BAM), we skip adapter detection and assume reads are ready or use alignment.
                # User explicitly requested to skip adapter detection for mapped reads.
                cached_adapter_pos = 0 
                
                buffer.extend(batch)
                
                if len(buffer) >= DISCOVERY_LIMIT:
                    debug_log.append(f"Buffered {len(buffer)} reads. Running primer discovery...")
                    
                    # Run discovery on buffer (parallel)
                    chunk_size = max(1, len(buffer) // threads)
                    chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
                    
                    kmer_futures = [executor.submit(discover_primers_from_batch_parallel, chunk, 19, 30, cached_adapter_pos) for chunk in chunks]
                    
                    for f in as_completed(kmer_futures):
                        p_counts, valid_n = f.result()
                        global_kmer_counts.update(p_counts)
                        total_valid_reads_for_discovery += valid_n
                    
                    # Identify Top Pairs / Positions
                    MIN_COUNT = max(10, 0.001 * len(buffer))
                    sorted_pairs = global_kmer_counts.most_common(200)
                    
                    top_kmers = []
                    candidates = []
                    
                    # Check if position-based (mapped) or sequence-based (unmapped)
                    if sorted_pairs and len(sorted_pairs[0][0]) == 3 and isinstance(sorted_pairs[0][0][1], int):
                        # Position-Based Clustering
                        by_chrom = defaultdict(list)
                        for pair, count in sorted_pairs:
                            if count < 5: continue
                            chrom, start, end = pair
                            by_chrom[chrom].append({'start': start, 'end': end, 'count': count, 'pair': pair})
                        
                        final_clusters = []
                        for chrom, items in by_chrom.items():
                            items.sort(key=lambda x: x['count'], reverse=True)
                            active_clusters = []
                            for item in items:
                                merged = False
                                for cluster in active_clusters:
                                    if abs(item['start'] - cluster['start']) < 50 and abs(item['end'] - cluster['end']) < 50:
                                        cluster['count'] += item['count']
                                        merged = True
                                        break
                                if not merged:
                                    active_clusters.append(item)
                            
                            for c in active_clusters:
                                if c['count'] >= MIN_COUNT:
                                    final_clusters.append(c)
                        
                        final_clusters.sort(key=lambda x: x['count'], reverse=True)
                        
                        for c in final_clusters[:50]: # Limit to top 50
                            chrom, start, end = c['pair']
                            coord_key = f"{chrom}:{start}-{end}"
                            gene_name = coord_key
                            
                            # Gene Mapping
                            if gene_models:
                                try:
                                    # Robust Chromosome Lookup
                                    target_chrom = None
                                    norm_chrom = chrom.replace("chr", "")
                                    
                                    if chrom in gene_models:
                                        target_chrom = chrom
                                    elif norm_chrom in gene_models:
                                        target_chrom = norm_chrom
                                    elif f"chr{norm_chrom}" in gene_models:
                                        target_chrom = f"chr{norm_chrom}"
                                    
                                    if target_chrom:
                                        overlaps = gene_models[target_chrom].overlap(start, end)
                                        
                                        best_gene = None
                                        found_genes = set()
                                        found_exons = defaultdict(set)
                                        
                                        for interval in overlaps:
                                            data = interval.data
                                            gname = data.get("name")
                                            ftype = data.get("type")
                                            if gname:
                                                found_genes.add(gname)
                                                best_gene = gname
                                                if ftype == "exon":
                                                    enum = data.get("exon")
                                                    if enum and enum != "?":
                                                        found_exons[gname].add(enum)
                                        
                                        if found_genes:
                                            gene_strs = []
                                            for g in sorted(found_genes):
                                                if g in found_exons and found_exons[g]:
                                                    try:
                                                        exs = sorted([int(e) for e in found_exons[g]])
                                                        if len(exs) > 1: ex_str = f"ex{min(exs)}-{max(exs)}"
                                                        else: ex_str = f"ex{exs[0]}"
                                                        gene_strs.append(f"{g}_{ex_str}")
                                                    except:
                                                        gene_strs.append(f"{g}_ex{','.join(sorted(found_exons[g]))}")
                                                else:
                                                    gene_strs.append(g)
                                            gene_name = f"{chrom}:{start}-{end}({','.join(gene_strs)})"
                                        elif best_gene:
                                            gene_name = best_gene
                                except Exception as e:
                                    print(f"DEBUG: Gene resolution error: {e}")
                                    pass
                            
                            # Store mapping for consistency
                            if gene_name != coord_key:
                                coord_to_name_map[coord_key] = gene_name

                            top_kmers.append(c['pair'])
                            candidates.append((gene_name, c['count']))
                            
                    else:
                        # Sequence Based
                        for pair, count in sorted_pairs:
                            if count >= MIN_COUNT:
                                top_kmers.append(pair)
                                candidates.append((f"{pair[0]}...{pair[1]}", count))
                    
                    debug_log.append(f"Identified {len(top_kmers)} candidate amplicon pairs from buffer.")
                    result_payload["suggested_primers"] = candidates
                    
                    # Process Buffer with identified pairs
                    debug_log.append("Processing buffered reads...")
                    match_futures = [executor.submit(match_kmer_pairs_batch, chunk, top_kmers, 150, 3, cached_adapter_pos) for chunk in chunks]
                    
                    for f in as_completed(match_futures):
                        m_found, local_s, local_map = f.result()
                        with stats_lock:
                            for name, d in local_s.items():
                                s = discovery_stats[name]
                                s["count"] += d["count"]
                                s["total_accuracy"] += d["total_accuracy"]
                                s["lengths"].extend(d["lengths"])
                            read_amplicon_map.update(local_map)
                    
                    # Cluster similar amplicons
                    debug_log.append("Clustering similar amplicons...")
                    discovery_stats, read_amplicon_map = cluster_discovered_amplicons(
                        discovery_stats, read_amplicon_map, similarity_threshold=0.90
                    )
                    debug_log.append(f"Clustering complete. Final amplicon count: {len(discovery_stats)}")
                    
                    if partial_cb:
                        partial_cb({
                            "amplicons": discovery_stats,
                            "read_amplicon_map": dict(read_amplicon_map)
                        })
                        
                    buffer = []
                    global_kmer_counts.clear()
                    gc.collect()
                    discovery_done = True
                    debug_log.append("Discovery complete. Switching to streaming analysis.")
            else:
                # Phase 2: Streaming Analysis
                future = executor.submit(match_kmer_pairs_batch, batch, top_kmers, 150, 3, cached_adapter_pos)
                futures.append(future)
                
                if len(futures) > threads * 2:
                    from concurrent.futures import wait, FIRST_COMPLETED
                    done, _ = wait(futures, return_when=FIRST_COMPLETED)
                    for f in done:
                        m_found, local_s, local_map = f.result()
                        total_matches += m_found
                        with stats_lock:
                            for name, d in local_s.items():
                                s = discovery_stats[name]
                                s["count"] += d["count"]
                                s["total_accuracy"] += d["total_accuracy"]
                                s["lengths"].extend(d["lengths"])
                            read_amplicon_map.update(local_map)
                        futures.remove(f)

        # Process remaining futures
        for f in as_completed(futures):
            m_found, local_s, local_map = f.result()
            total_matches += m_found
            with stats_lock:
                for name, d in local_s.items():
                    s = discovery_stats[name]
                    s["count"] += d["count"]
                    s["total_accuracy"] += d["total_accuracy"]
                    s["lengths"].extend(d["lengths"])
                read_amplicon_map.update(local_map)
        
        # If discovery_done was never reached
        if not discovery_done and buffer:
            debug_log.append(f"Total reads ({len(buffer)}) less than DISCOVERY_LIMIT. Running primer discovery on all available reads.")
            chunk_size = max(1, len(buffer) // threads)
            chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
            
            kmer_futures = [executor.submit(discover_primers_from_batch_parallel, chunk, 19, 30, cached_adapter_pos) for chunk in chunks]
            for f in as_completed(kmer_futures):
                k_counts, valid_n = f.result()
                global_kmer_counts.update(k_counts)
                total_valid_reads_for_discovery += valid_n
            
            ADAPTER_THRESHOLD = 0.25 * len(buffer)
            sorted_kmers = global_kmer_counts.most_common(200)
            candidates = []
            adapters = []
            for kmer, count in sorted_kmers:
                if count > ADAPTER_THRESHOLD: adapters.append((kmer, count))
                else: candidates.append((kmer, count))
            
            # Handle mapped vs unmapped for top_kmers
            top_kmers = []
            if candidates and len(candidates[0][0]) == 3 and isinstance(candidates[0][0][1], int):
                 # Mapped reads logic (already handled above)
                 pass 
            else:
                 # Unmapped reads: Cluster similar k-mers
                 candidates = cluster_kmers(candidates, dist_threshold=3)
                 
                 for kmer, count in candidates[:50]:
                     top_kmers.append(kmer)

            result_payload["adapters_found"] = adapters
            result_payload["suggested_primers"] = candidates[:50]

            match_futures = [executor.submit(match_kmer_pairs_batch, chunk, top_kmers, 150, 3, cached_adapter_pos) for chunk in chunks]
            for f in as_completed(match_futures):
                batch_matches, local_stats, local_map = f.result()
                total_matches += batch_matches
                with stats_lock:
                    for name, data in local_stats.items():
                        discovery_stats[name]["count"] += data["count"]
                        discovery_stats[name]["total_accuracy"] += data["total_accuracy"]
                        discovery_stats[name]["lengths"].extend(data["lengths"])
                    read_amplicon_map.update(local_map)
            buffer = []
            global_kmer_counts.clear()
            gc.collect()
        
        debug_log.append(f"Pairing Logic: Scanned {streamer.total_reads} reads, found {total_matches} valid primer pairs.")
    
    # Update read_amplicon_map with resolved gene names
    # This ensures that per-read data matches the final table names for filtering
    if coord_to_name_map:
        # We need to be careful about concurrency if this is running while threads are active?
        # No, threads are joined by now (executor context exit).
        for read_id, name in read_amplicon_map.items():
            # Handle potential concatemer suffix
            is_concat = False
            base_name = name
            if "|CONCAT" in name:
                is_concat = True
                base_name = name.replace("|CONCAT", "")
            
            if base_name in coord_to_name_map:
                new_name = coord_to_name_map[base_name]
                if is_concat:
                    new_name += "|CONCAT"
                read_amplicon_map[read_id] = new_name
        
    return internal_adapter_count, discovery_stats, total_valid_reads_for_discovery, global_kmer_counts, coord_to_name_map


def _run_sequential_mode(streamer, partial_cb, progress_cb, is_discovery,
                         debug_log, global_kmer_counts, total_valid_reads_for_discovery,
                         primer_list, gene_models, amplicon_stats):
    """Helper for sequential analysis."""
    cached_adapter_pos = None
    first_batch = True
    last_update_time = time.time()
    
    for batch, meta_batch in streamer.stream_batches(progress_cb):
        current_time = time.time()
        if partial_cb: 
            payload = {"metadata": meta_batch}
            
            if not is_discovery and current_time - last_update_time > 50.0:
                formatted_amps = {}
                for name, d in amplicon_stats.items():
                    c = d["count"]
                    if c > 0:
                        formatted_amps[name] = {
                            "count": c, 
                            "average_accuracy": d["total_accuracy"]/c,
                            "median_length": 0,
                            "stdev_length": 0,
                            "raw_lengths": []
                        }
                payload["amplicons"] = formatted_amps
                last_update_time = current_time
            
            partial_cb(payload)
        if is_discovery:
            # Detect adapter position from first batch
            if first_batch:
                debug_log.append("Detecting adapter position from first batch (up to 10000 reads)...")
                cached_adapter_pos, n_detected = detect_adapter_position(batch, max_sample=10000)
                if cached_adapter_pos:
                    debug_log.append(f"Adapter position detected: {cached_adapter_pos}bp (from {n_detected} reads, will use for fast trimming)")
                else:
                    debug_log.append("No adapters detected (using alignment for all reads)")
                first_batch = False
            
            k_counts, valid_n = discover_primers_from_batch(batch, cached_adapter_pos=cached_adapter_pos)
            global_kmer_counts.update(k_counts)
            total_valid_reads_for_discovery += valid_n
        else:
            process_amplicon_batch(batch, primer_list, gene_models, amplicon_stats)
            
    return total_valid_reads_for_discovery




def run_analysis(bam_path, primers_dict, gene_models, filter_settings, threads, progress_cb, partial_cb, stop_check_cb=None, qc_only=False, primer_tolerance=0):
    
    print(f"DEBUG: Starting analysis for {bam_path}")
    streamer = ns_core.get_streamer(bam_path, filter_settings, threads=threads)
    is_discovery = primers_dict is None
    amplicon_stats = defaultdict(lambda: {"count": 0, "total_accuracy": 0.0, "lengths": [], "positions": []})  # "total_qs": 0.0 - TODO: Add QS tracking
    stats_lock = Lock()  # Thread-safe stats updates
    active_primer_pairs = []
    global_kmer_counts = Counter()
    internal_adapter_count = 0 
    
    # Debug tracking
    debug_log = []
    read_amplicon_map = {}  # Track {read_id: amplicon_name} for QS tracking
    
    # Initialize result_payload early to avoid UnboundLocalError
    result_payload = {
        "summary": {}, "raw_accuracies": [], "read_lengths": [],
        "amplicons": {}, "adapters_found": [], "suggested_primers": [], "debug_log": debug_log,
        "read_amplicon_map": read_amplicon_map  # Include mapping for downstream QS tracking
    }
    total_valid_reads_for_discovery = 0
    last_update_time = time.time()
    
    if not is_discovery:
        # Just create a flat list of primers to check combinatorially
        primer_list = []
        for name, seq in primers_dict.items():
            if seq:
                # Store (PrimerID, Sequence)
                # Strip suffixes for display name cleanliness if preferred, 
                # but user gave IDs, so we use them.
                # Actually user said "names = sorted... replace _FWD" etc.
                # But for combinatorial, we want the exact primer ID to report P1-P2.
                primer_list.append((name, seq.upper()))
    
    # Precompute adapter k-mers for fast concatemer detection (all adapters + reverse complements)
    adapter_kmers = set()
    kmer_size = 15
    for adapter in ONT_ADAPTERS:  # Use ALL adapters for maximum sensitivity
        # Add forward k-mers
        for i in range(len(adapter) - kmer_size + 1):
            adapter_kmers.add(adapter[i:i+kmer_size])
        
        # Add reverse complement k-mers
        rc_adapter = reverse_complement(adapter)
        for i in range(len(rc_adapter) - kmer_size + 1):
            adapter_kmers.add(rc_adapter[i:i+kmer_size])

    # Parallel processing for non-discovery mode
    if qc_only:
        # QC Only Mode: Just stream and report metadata
        debug_log.append("QC Only Mode: Skipping primer detection/analysis.")
        last_update_time = time.time()
        for batch, meta_batch in streamer.stream_batches(progress_cb):
            if stop_check_cb and stop_check_cb():
                break
            
            # Live updates
            current_time = time.time()
            if partial_cb:
                payload = {"metadata": meta_batch}
                # No amplicons to report
                partial_cb(payload)
                
            # Detect concatemers (optional for QC? Let's keep it as it's fast and useful for QC)
            batch_concatemer_count, concatemer_ids = detect_internal_adapters(batch, adapter_kmers, kmer_size)
            internal_adapter_count += batch_concatemer_count
            
    elif not is_discovery and threads > 1:
        internal_adapter_count = _run_parallel_primer_mode(
            streamer, threads, stop_check_cb, partial_cb, progress_cb,
            adapter_kmers, kmer_size, primer_list, gene_models,
            stats_lock, amplicon_stats, read_amplicon_map, primer_tolerance
        )
        
        # --- FALLBACK CHECK ---
        total_assigned = sum(d["count"] for d in amplicon_stats.values())
        total_reads = streamer.total_reads
        
        # Threshold: If < 5% reads assigned and we processed a decent amount, try discovery
        if total_reads > 500 and (total_assigned / total_reads) < 0.05:
             print(f"DEBUG: Low primer assignment ({total_assigned}/{total_reads}). Falling back to Discovery Mode.")
             debug_log.append(f"WARNING: Only {total_assigned} reads assigned to primers. Falling back to K-mer Discovery.")
             
             # Reset for Discovery
             is_discovery = True
             streamer.reset() # Reset stream to start
             amplicon_stats.clear()
             read_amplicon_map.clear()
             
             # Run Discovery Mode
             internal_adapter_count, discovery_stats, total_valid_reads_for_discovery, global_kmer_counts, coord_to_name_map = _run_parallel_discovery_mode(
                streamer, threads, stop_check_cb, partial_cb, progress_cb,
                adapter_kmers, kmer_size, gene_models,
                stats_lock, read_amplicon_map, debug_log, result_payload
            )

    elif is_discovery and threads > 1:
        internal_adapter_count, discovery_stats, total_valid_reads_for_discovery, global_kmer_counts, coord_to_name_map = _run_parallel_discovery_mode(
            streamer, threads, stop_check_cb, partial_cb, progress_cb,
            adapter_kmers, kmer_size, gene_models,
            stats_lock, read_amplicon_map, debug_log, result_payload
        )
            
    else:
        # Sequential processing (single thread)
        # Note: Sequential mode doesn't return map yet, so we use empty dict
        total_valid_reads_for_discovery = _run_sequential_mode(
            streamer, partial_cb, progress_cb, is_discovery,
            debug_log, global_kmer_counts, total_valid_reads_for_discovery,
            primer_list, gene_models, amplicon_stats
        )
        coord_to_name_map = {}

    stats = streamer.get_summary()
    raw_accuracies = stats.pop("raw_accuracies", [])
    read_lengths = stats.pop("read_lengths", [])
    stats["internal_adapter_count"] = internal_adapter_count
    
    # --- DEBUG INFO ---
    debug_log.append(f"BAM Analysis Finished. Total reads: {stats['total_reads_processed']}")
    if is_discovery:
        debug_log.append(f"Discovery Mode: Analyzed {total_valid_reads_for_discovery} reads with valid ends.")
        debug_log.append(f"Total unique K-mers found: {len(global_kmer_counts)}")

    # Update result_payload with final stats
    result_payload.update({
        "summary": stats, "raw_accuracies": raw_accuracies, "read_lengths": read_lengths,
        "debug_log": debug_log
    })
    
    if is_discovery:
        # Finalize discovery stats
        final_stats = finalize_discovery_stats(discovery_stats, gene_models, known_names=coord_to_name_map)
        result_payload["amplicons"] = final_stats
        debug_log.append(f"Discovered {len(final_stats)} unique amplicon pairs.")
        result_payload["debug_log"] = debug_log # Update log
        
    else:
        formatted_amps = {}
        for name, d in amplicon_stats.items():
            c = d["count"]
            if c > 0:
                lens = np.array(d["lengths"])
                
                # Parse region from amplicon name (format: "PrimerA-PrimerB:chr:start-end" or "PrimerA-PrimerB:chr:start-end(Gene)")
                gene_info = parse_amplicon_name(name)
                region = None
                if gene_info["chrom"] and gene_info["start"] is not None and gene_info["end"] is not None:
                    region = f"{gene_info['chrom']}:{gene_info['start']}-{gene_info['end']}"
                
                formatted_amps[name] = {
                    "count": c, "average_accuracy": d["total_accuracy"]/c,
                    "median_length": np.median(lens), "stdev_length": np.std(lens), "raw_lengths": d["lengths"],
                    "chrom": gene_info["chrom"],
                    "start": gene_info["start"],
                    "end": gene_info["end"],
                    "region": region,
                    "gene_name": gene_info["gene_name"],
                    "genes": gene_info["genes"]
                }
        result_payload["amplicons"] = formatted_amps

    return result_payload