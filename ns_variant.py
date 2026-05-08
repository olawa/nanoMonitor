# Filename: ns_variant.py
# Created: 2025-11-20 15:00 CET

import os
import pysam
import subprocess
import tempfile
import concurrent.futures
from PyQt6.QtCore import QThread, pyqtSignal


def _proc_variant_task(bam_path, task, min_af):
    """Top-level helper for ProcessPoolExecutor to ensure pickleability."""
    name, chrom, start, end, target_ids = task
    caller = SimplePileupCaller(bam_path, target_ids)
    return name, caller.call_variants(chrom, start, end, min_af)


class SimplePileupCaller:
    """
    Calculates allele frequencies using pysam pileup.
    Fast and sufficient for amplicons.
    """
    def __init__(self, bam_path, target_read_ids=None):
        self.bam_path = bam_path
        self.target_read_ids = target_read_ids

    def call_variants(self, chrom, start, end, min_af=0.02):
        """
        Runs pileup on region and returns variants above min_af.
        """
        variants = []
        try:
            with pysam.AlignmentFile(self.bam_path, "rb") as sam:
                # Get reference sequence for the region if possible
                # But for now, we'll just use the majority base as 'REF' if we don't have a fasta
                # Actually, sam.pileup gives us base counts.
                
                for pileupcolumn in sam.pileup(chrom, start, end, truncate=True, stepper="samtools"):
                    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                    total = 0
                    
                    for pileupread in pileupcolumn.pileups:
                        if self.target_read_ids and pileupread.alignment.query_name not in self.target_read_ids:
                            continue
                            
                        if not pileupread.is_del and not pileupread.is_refskip:
                            base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                            if base in counts:
                                counts[base] += 1
                                total += 1
                    
                    if total < 10: continue # Minimum depth
                    
                    # Find potential variants (non-majority bases)
                    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
                    ref_base, ref_count = sorted_counts[0]
                    
                    for alt_base, alt_count in sorted_counts[1:]:
                        if alt_count == 0: continue
                        af = alt_count / total
                        if af >= min_af:
                            variants.append({
                                'chrom': chrom,
                                'pos': pileupcolumn.pos + 1,
                                'ref': ref_base,
                                'alt': alt_base,
                                'af': af,
                                'depth': total
                            })
        except Exception as e:
            print(f"Pileup Error: {e}")
            
        return variants


class VariantWorker(QThread):
    """Thread to run the variant calling process without freezing GUI."""
    finished = pyqtSignal(bool, str, list) # success, message, variants_list
    
    def __init__(self, bam_path, chrom, start, end, target_read_ids=None, min_af=0.02):
        super().__init__()
        self.bam_path = bam_path
        self.chrom = chrom
        self.v_start = start
        self.v_end = end
        self.target_read_ids = target_read_ids
        self.min_af = min_af

    def run(self):
        try:
            caller = SimplePileupCaller(self.bam_path, self.target_read_ids)
            variants = caller.call_variants(self.chrom, self.v_start, self.v_end, self.min_af)
            
            if not variants:
                self.finished.emit(True, "No variants found above threshold.", [])
                return
                
            msg = f"Found {len(variants)} variants in {self.chrom}:{self.v_start}-{self.v_end}"
            self.finished.emit(True, msg, variants)
                
        except Exception as e:
            self.finished.emit(False, str(e), [])

class BatchVariantWorker(QThread):
    """Thread to run variant calling for multiple amplicons sequentially."""
    finished = pyqtSignal(dict) # {name: [variants]}
    progress = pyqtSignal(int, int) # current, total
    
    def __init__(self, bam_path, tasks, min_af=0.02):
        """
        tasks: list of (name, chrom, start, end, target_ids)
        """
        super().__init__()
        self.bam_path = bam_path
        self.tasks = tasks
        self.min_af = min_af

    def run(self):
        results = {}
        total = len(self.tasks)
        
        # Parallelize using ProcessPoolExecutor for >100% CPU (bypassing GIL)
        # max_workers=8 is good for modern Macs
        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
            # Map tasks to executor
            futures = [
                executor.submit(_proc_variant_task, self.bam_path, task, self.min_af) 
                for task in self.tasks
            ]
            
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                try:
                    name, variants = future.result()
                    results[name] = variants
                except Exception as e:
                    print(f"Batch Process Error: {e}")
                    # We might not know which 'name' failed here easily without result
                
                completed += 1
                self.progress.emit(completed, total)
        
        self.finished.emit(results)