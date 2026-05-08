# Filename: rna_logic.py
# Created: 2025-11-20 11:00 CET

import pysam
import time
from collections import defaultdict
import os

def load_gene_bed(filepath):
    """
    Loads a simple 4-column BED file: chrom, start, end, gene_name
    Returns a list of gene tuples.
    """
    genes = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        genes.append({
                            "chrom": parts[0],
                            "start": int(parts[1]),
                            "end": int(parts[2]),
                            "name": parts[3]
                        })
    except Exception as e:
        print(f"Error loading BED file: {e}")
        return []
    return genes

def run_rna_analysis(bam_file_path, genes_list, filter_obj, progress_callback, results_callback, partial_callback=None):
    """
    Counts reads mapping to specific genes.
    Now supports partial updates and accuracy collection.
    """
    gene_counts = defaultdict(int)
    raw_accuracies = [] # NEW: Collect accuracy for plot
    total_reads = 0
    counted_reads = 0
    start_time = time.time()
    
    CHUNK_SIZE = 200000 # Updated to 200k
    
    try:
        with pysam.AlignmentFile(bam_file_path, "rb", threads=4) as bamfile:
            
            genes_by_chrom = defaultdict(list)
            for g in genes_list:
                genes_by_chrom[g["chrom"]].append(g)

            read_buffer = []

            def process_buffer(reads_to_process):
                nonlocal counted_reads
                filtered_reads = filter_obj.process_chunk(reads_to_process)
                
                for read in filtered_reads:
                    # Collect accuracy
                    try:
                        de_tag = read.get_tag("de")
                        raw_accuracies.append((1 - de_tag) * 100.0)
                    except KeyError:
                        raw_accuracies.append(0.0)

                    if read.reference_name in genes_by_chrom:
                        ref_start = read.reference_start
                        ref_end = read.reference_end
                        
                        for gene in genes_by_chrom[read.reference_name]:
                            if max(ref_start, gene["start"]) < min(ref_end, gene["end"]):
                                gene_counts[gene["name"]] += 1
                                counted_reads += 1
            
            for read in bamfile:
                total_reads += 1
                
                read_buffer.append(read)
                
                if len(read_buffer) >= CHUNK_SIZE:
                    process_buffer(read_buffer)
                    read_buffer = []
                    
                    # Emit Progress & Partial Data
                    progress_callback(total_reads)
                    if partial_callback:
                        partial_callback({
                            "raw_accuracies": raw_accuracies,
                            "genes": dict(gene_counts)
                        })
            
            if read_buffer:
                process_buffer(read_buffer)

    except Exception as e:
        results_callback(None, f"Error in RNA analysis: {e}")
        return

    final_results = {
        "summary": {
            "total_reads": total_reads,
            "counted_reads": counted_reads,
            "time_s": time.time() - start_time
        },
        "genes": dict(gene_counts),
        "raw_accuracies": raw_accuracies # Return final accuracies
    }
    
    results_callback(final_results, None)