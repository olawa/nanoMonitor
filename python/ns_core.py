# Filename: ns_core.py
# Created: 2025-11-21 19:00 CET

import pysam
import time
import os
import numpy as np

# --- UTILITIES ---
def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

def detect_adapter_position(reads_sample, max_sample=10000):
    """
    Detect the most common adapter end position from a sample of reads.
    Returns the median position, or None if no adapters found.
    """
    from ns_resources import ONT_ADAPTERS
    import edlib
    
    positions = []
    for i, read in enumerate(reads_sample):
        if i >= max_sample:
            break
        if read.is_unmapped:
            seq = read.query_sequence
            if not seq:
                continue
            for adapter in ONT_ADAPTERS:
                prefix = seq[:100]
                res = edlib.align(adapter, prefix, mode="HW", task="locations", k=3)
                if res["editDistance"] > -1 and res["locations"]:
                    end_pos = res["locations"][-1][1] + 1
                    positions.append(end_pos)
                    break  # Found adapter, move to next read
    
    if positions:
        # Return median position (more robust than mean)
        return int(np.median(positions)), len(positions)
    return None, 0

def detect_internal_adapters(reads_batch, adapter_kmers=None, kmer_size=15, min_internal_dist=500):
    """
    Fast detection of internal adapters (concatemers) using k-mer matching.
    
    Args:
        reads_batch: List of reads to check
        adapter_kmers: Set of k-mers from adapter sequences (precomputed for speed)
        kmer_size: Size of k-mers to use (default 15 for speed/accuracy balance)
        min_internal_dist: Minimum distance from read ends to consider "internal" (default 500bp)
    
    Returns:
        tuple: (count, list of read IDs with internal adapters)
    """
    from ns_resources import ONT_ADAPTERS
    
    if adapter_kmers is None:
        # Precompute adapter k-mers from ALL ONT_ADAPTERS
        adapter_kmers = set()
        for adapter in ONT_ADAPTERS:  # Use ALL adapters for maximum sensitivity
            # Add forward k-mers
            for i in range(len(adapter) - kmer_size + 1):
                adapter_kmers.add(adapter[i:i+kmer_size])
            
            # Add reverse complement k-mers
            rc_adapter = reverse_complement(adapter)
            for i in range(len(rc_adapter) - kmer_size + 1):
                adapter_kmers.add(rc_adapter[i:i+kmer_size])
    
    internal_count = 0
    concatemer_read_ids = []
    
    for read in reads_batch:
        seq = read.query_sequence
        if not seq or len(seq) < min_internal_dist * 2:
            continue
        
        # Only check middle section (exclude ends)
        middle_seq = seq[min_internal_dist:-min_internal_dist]
        
        # Fast k-mer scan
        for i in range(len(middle_seq) - kmer_size + 1):
            kmer = middle_seq[i:i+kmer_size]
            if kmer in adapter_kmers:
                internal_count += 1
                concatemer_read_ids.append(read.query_name)
                break  # Count each read only once
    
    return internal_count, concatemer_read_ids

def calculate_pore_stats(all_meta):
    """
    Calculate pore statistics (gap times) from read metadata.
    Returns a dictionary of stats or None if calculation fails.
    """
    try:
        import pandas as pd
        df = pd.DataFrame(all_meta)
        if 'ch' in df.columns and 'mx' in df.columns and 'st' in df.columns:
            # Ensure numeric types for sorting
            df['ch'] = pd.to_numeric(df['ch'], errors='coerce').fillna(-1).astype(int)
            df['mx'] = pd.to_numeric(df['mx'], errors='coerce').fillna(-1).astype(int)
            df['st'] = pd.to_numeric(df['st'], errors='coerce')
            
            # Filter valid tags
            df = df[(df['ch'] != -1) & (df['st'] > 0)].copy()
            if not df.empty:
                df = df.sort_values(['ch', 'mx', 'st'])
                # Calculate gaps
                df['gap'] = df.groupby(['ch', 'mx'])['st'].diff()
                
                # User feedback: 400bps is approx, so gaps can be negative.
                # We care about LARGE gaps.
                # Filter out NaNs (first read in each group) and invalid gaps
                valid_gaps = df[(df['gap'].notna()) & (df['gap'] < 3600)].copy()
                
                if not valid_gaps.empty:
                    # For stats, maybe we treat negative gaps as 0? 
                    # Or just report raw mean (which might be slightly lower due to negatives)
                    # Let's report raw mean/median but also count "Long Gaps" (> 60s)
                    
                    mean_gap = valid_gaps['gap'].mean()
                    median_gap = valid_gaps['gap'].median()
                    
                    # Count long gaps
                    long_gaps = valid_gaps[valid_gaps['gap'] > 60].shape[0]
                    total_gaps = valid_gaps.shape[0]
                    long_gap_pct = (long_gaps / total_gaps * 100) if total_gaps > 0 else 0
                    
                    # Per channel stats
                    ch_stats = valid_gaps.groupby(['ch', 'mx'])['gap'].agg(['mean', 'median', lambda x: (x > 60).sum()]).reset_index()
                    ch_stats.columns = ['ch', 'mx', 'mean', 'median', 'long_gaps']
                    
                    return {
                        'global_mean_gap': float(mean_gap),
                        'global_median_gap': float(median_gap),
                        'long_gap_pct': float(long_gap_pct),
                        'channel_stats': ch_stats.to_dict('records')
                    }
    except ImportError:
        pass # No pandas, skip analysis
    except Exception as e:
        print(f"Pore analysis failed: {e}")
    return None



# --- FILTERING ---
class ReadFilter:
    """Handles filtering of BAM reads based on tags and properties."""
    def __init__(self, min_qs=0, min_length=0, allow_unmapped=False, duplex_only=False):
        self.min_qs = min_qs
        self.min_length = min_length
        self.allow_unmapped = allow_unmapped
        self.duplex_only = duplex_only

    def passes(self, read):
        """Checks thresholds. Returns True if read should be kept."""
        if read.is_unmapped and not self.allow_unmapped:
            return False
        if self.duplex_only:
            # Check for dx:i:1
            try:
                if read.get_tag("dx") != 1: return False
            except KeyError:
                return False
        if read.query_length < self.min_length:
            return False
        if self.min_qs > 0:
            try:
                qs = read.get_tag("qs")
                if qs < self.min_qs:
                    return False
            except KeyError:
                pass 
        return True

# --- BAM/FASTQ ENGINE ---

class FastqRead:
    """
    Duck-typing wrapper for FASTQ records to mimic pysam.AlignedSegment.
    Parses tags from the comment line (e.g., qs:f:24.7 dx:i:1).
    """
    def __init__(self, record):
        self.query_name = record.name
        self.query_sequence = record.sequence
    def __init__(self, name, sequence, quality, comment):
        self.query_name = name
        self.query_sequence = sequence
        self.query_length = len(sequence)
        self.is_unmapped = True # FASTQ reads are always unmapped
        self.is_reverse = False
        self.reference_name = None
        self.reference_start = None
        self.reference_end = None
        self.mapping_quality = 0
        self._quality = quality
        self.tags = {}
        
        # Mapping defaults (unmapped)
        self.reference_id = -1
        self.reference_start = -1
        self.is_reverse = False
        
        # Parse comment tags
        if comment:
            # User report: tags can be tab separated. split() handles both space and tab.
            parts = comment.split()
            for p in parts:
                if ':' in p:
                    # Expecting tag:type:value
                    fields = p.split(':', 2) # Limit split to 2 to handle colons in value (e.g. timestamps)
                    if len(fields) >= 3:
                        tag = fields[0]
                        val_type = fields[1]
                        val = fields[2]
                        
                        # Type conversion
                        if val_type == 'i':
                            try: self.tags[tag] = int(val)
                            except: pass
                        elif val_type == 'f':
                            try: self.tags[tag] = float(val)
                            except: pass
                        else:
                            self.tags[tag] = val      # Convert 'st' to timestamp if present
        if 'st' in self.tags:
            try:
                # Remove Z if present (though usually +00:00)
                ts_str = self.tags['st']
                # Basic ISO parsing
                from datetime import datetime
                dt = datetime.fromisoformat(ts_str)
                self.tags['st_ts'] = dt.timestamp()
            except Exception:
                pass

    def get_tag(self, tag):
        if tag in self.tags:
            return self.tags[tag]
            
        # Fallback/Emulation
        if tag == "qs":
            # If not in tags, calculate from quality string
            if not self._quality: return 0
            q_scores = np.frombuffer(self._quality.encode(), dtype=np.int8) - 33
            # Mean Error Probability: -10 * log10( mean( 10^(-Q/10) ) )
            p_scores = 10.0 ** (-q_scores / 10.0)
            mean_p = np.mean(p_scores)
            return -10.0 * np.log10(mean_p)
        if tag == "de":
            # Estimate error from qs if available
            qs = self.get_tag("qs")
            return 10 ** (-qs / 10.0)
            
        raise KeyError(tag)
        
    @property
    def query_alignment_start(self): return 0
    @property
    def query_alignment_end(self): return self.query_length
    
    def has_tag(self, tag):
        # Only return True for tags that exist in self.tags or can be emulated
        if tag in self.tags:
            return True
        # Emulated tags that get_tag can handle
        if tag in ["qs", "de"]:
            return True
        # st_ts is derived from st
        if tag == "st_ts" and "st" in self.tags:
            return True
        return False


class BaseStreamer:
    """Base class for file streamers."""
    def __init__(self, filepath, filter_settings=None, threads=8, chunk_size=2000, allow_unmapped=False):
        self.filepath = filepath
        self.threads = threads
        self.chunk_size = chunk_size
        
        filters = filter_settings if filter_settings else {}
        self.read_filter = ReadFilter(
            min_qs=filters.get("min_qs", 0),
            min_length=filters.get("min_len", 0),
            allow_unmapped=allow_unmapped,
            duplex_only=filters.get("duplex_only", False)
        )
        
        self.total_reads = 0
        self.mapped_reads = 0
        self.start_time = time.time()
        self.raw_accuracies = [] 
        self.read_lengths = []   
        self.sv_links = []       
        self.chrom_lengths = {}  

    def get_summary(self):
        return {
            "total_reads_processed": self.total_reads,
            "mapped_passed_reads": self.mapped_reads,
            "processing_time_s": time.time() - self.start_time,
            "raw_accuracies": self.raw_accuracies,
            "read_lengths": self.read_lengths,
            "sv_links": self.sv_links,
            "chrom_lengths": self.chrom_lengths
        }
        
    def reset(self):
        """Reset counters for re-analysis."""
        self.total_reads = 0
        self.mapped_reads = 0
        self.raw_accuracies = [] 
        self.read_lengths = []   
        self.sv_links = []       
        self.chrom_lengths = {}
        self.start_time = time.time()
        
    def _process_read(self, read, read_len, acc, qs, extract_sequences=False):
        # Common processing for stats
        self.total_reads += 1
        self.raw_accuracies.append(acc)
        self.read_lengths.append(read_len)
        
        dx = 0
        try: dx = read.get_tag("dx")
        except KeyError: pass
        
        ch = read.get_tag("ch") if read.has_tag("ch") else -1
        mx = read.get_tag("mx") if read.has_tag("mx") else -1
        
        st = 0.0
        if read.has_tag("st_ts"):
            st = read.get_tag("st_ts")
        elif read.has_tag("st"):
            # BAM usually has 'st' as start time (samples or timestamp)
            # We use it for sorting, so any monotonic value works.
            st_val = read.get_tag("st")
            if isinstance(st_val, str):
                try:
                    from datetime import datetime
                    dt = datetime.fromisoformat(st_val.replace("Z", "+00:00"))
                    st = dt.timestamp()
                except:
                    pass
            else:
                try:
                    st = float(st_val)
                except:
                    pass
            
        # Include sequence ends for duplex discovery (optimized memory)
        # Store first 500bp (Head) and last 500bp (Tail)
        head_seq = ""
        tail_seq = ""
        
        if extract_sequences:
            full_seq = read.query_sequence
            seq_len = len(full_seq)
            limit = 500
            
            if seq_len <= limit:
                head_seq = full_seq
                tail_seq = full_seq
            else:
                head_seq = full_seq[:limit]
                tail_seq = full_seq[-limit:]
            
        # Mapping info
        rid = getattr(read, 'reference_id', -1)
        pos = getattr(read, 'reference_start', -1)
        rev = getattr(read, 'is_reverse', False)
        
        return {
            'acc': acc, 'qs': qs, 'len': read_len, 'dx': dx, 'ch': ch, 'mx': mx, 'st': st, 
            'head': head_seq, 'tail': tail_seq, 'id': read.query_name,
            'rid': rid, 'pos': pos, 'rev': rev
        }


class BamStreamer(BaseStreamer):
    """
    Engine to read BAM files.
    """
    def stream_batches(self, progress_callback=None, extract_sequences=False):
        read_buffer = []
        meta_buffer = []
        
        try:
            with pysam.AlignmentFile(self.filepath, "rb", check_sq=False, threads=self.threads) as bamfile:
                # Check for SQ lines to determine if we should allow unmapped reads by default
                if not bamfile.header.get("SQ"):
                    self.read_filter.allow_unmapped = True
                    
                self.chrom_lengths = dict(zip(bamfile.references, bamfile.lengths))
                
                for read in bamfile:
                    # Extract Metadata
                    read_len = read.query_length
                    acc = 0.0
                    qs = 0.0
                    
                    # Accuracy from 'de' tag
                    try: 
                        acc = (1 - read.get_tag("de")) * 100.0
                    except KeyError: 
                        pass
                    # Fallback to mean Q-score if 'de' tag is missing
                        #if read.query_qualities:
                            # Calculate Q-score first
                         ##   qs_val = sum(read.query_qualities) / len(read.query_qualities)
#acc = (1 - 10**(-qs_val/10.0)) * 100.0
                            
                    # Q-Score from 'qs' tag
                    try:
                        qs = float(read.get_tag("qs"))
                    except KeyError:
                        # Fallback to calculation
                        if read.query_qualities:
                            # Mean Error Probability
                            q_scores = np.array(read.query_qualities, dtype=np.float64)
                            p_scores = 10.0 ** (-q_scores / 10.0)
                            mean_p = np.mean(p_scores)
                            qs = -10.0 * np.log10(mean_p)
                    
                    meta = self._process_read(read, read_len, acc, qs, extract_sequences)
                    meta_buffer.append(meta)
                    
                    if progress_callback and self.total_reads % self.chunk_size == 0:
                        progress_callback(self.total_reads)

                    if not self.read_filter.passes(read): continue

                    self.mapped_reads += 1
                    
                    # SA Tag Parsing
                    if not read.is_unmapped and read.has_tag("SA"):
                        sa_tag = read.get_tag("SA")
                        splits = sa_tag.split(';')
                        for split in splits:
                            if not split: continue
                            parts = split.split(',')
                            if len(parts) >= 2:
                                dest_chrom = parts[0]
                                try:
                                    dest_pos = int(parts[1])
                                    if dest_chrom in self.chrom_lengths:
                                        self.sv_links.append(
                                            ((read.reference_name, read.reference_start), 
                                             (dest_chrom, dest_pos),
                                             read.query_name)
                                        )
                                except ValueError: pass

                    read_buffer.append(read)

                    if len(read_buffer) >= self.chunk_size:
                        yield read_buffer, meta_buffer
                        read_buffer = []
                        meta_buffer = []

                if read_buffer or meta_buffer:
                    yield read_buffer, meta_buffer

        except Exception as e:
            print(f"BAM Stream Error: {e}")
            raise e


class FastqStreamer(BaseStreamer):
    """
    Engine to read FASTQ files (plain or gz).
    """
    def __init__(self, filepath, filter_settings=None, threads=8, chunk_size=2000):
        super().__init__(filepath, filter_settings, threads, chunk_size, allow_unmapped=True)

    def stream_batches(self, progress_callback=None, extract_sequences=False):
        read_buffer = []
        meta_buffer = []
        
        try:
            # pysam.FastxFile handles gzip automatically
            with pysam.FastxFile(self.filepath) as fh:
                for entry in fh:
                    read = FastqRead(entry.name, entry.sequence, entry.quality, entry.comment)
                    
                    # Extract Metadata
                    read_len = read.query_length
                    qs = read.get_tag("qs") # Use calculated Q-score
                    
                    # Calculate Accuracy from Q-Score
                    # acc = (1 - 10^(-qs/10)) * 100
                    acc = (1 - 10**(-qs/10.0)) * 100.0
                    
                    meta = self._process_read(read, read_len, acc, qs, extract_sequences)
                    meta_buffer.append(meta)
                    
                    if progress_callback and self.total_reads % self.chunk_size == 0:
                        progress_callback(self.total_reads)

                    if not self.read_filter.passes(read): continue

                    self.mapped_reads += 1 # Count as passed reads
                    
                    read_buffer.append(read)

                    if len(read_buffer) >= self.chunk_size:
                        yield read_buffer, meta_buffer
                        read_buffer = []
                        meta_buffer = []

                if read_buffer or meta_buffer:
                    yield read_buffer, meta_buffer

        except Exception as e:
            print(f"FASTQ Stream Error: {e}")
            raise e

def get_streamer(filepath, *args, **kwargs):
    """Factory to get appropriate streamer."""
    if filepath.endswith(".bam") or filepath.endswith(".cram"):
        return BamStreamer(filepath, *args, **kwargs)
    else:
        return FastqStreamer(filepath, *args, **kwargs)