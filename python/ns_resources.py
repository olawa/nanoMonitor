# Filename: ns_resources.py
# Created: 2025-11-21 16:00 CET

import pysam
import os
import csv
from collections import defaultdict
try:
    from intervaltree import IntervalTree
except ImportError:
    IntervalTree = None

# --- CONSTANTS ---
ONT_ADAPTERS = [
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT",
    "AGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT",
    "AGGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACCT",
    "AGGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAATAGGGAAACACGATAGAATCCGAACAGCACCT",
    "AGGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCT",
    "AGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACCT",
    "AGGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACCT",
    "AGGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACCT",
    "AGGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACCT",
    "AGGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACCT",
    "AGGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACCT",
    "AGGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACCT",
    "AGGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAGAACGACTTCCATACTCGTGTGACAGCACCT",
    "AGGTGCTGTCACACGAGTATGGAAGTCGTTCTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAACGAGTCTCTTGGGACCCATAGACAGCACCT",
    "AGGTGCTGTCTATGGGTCCCAAGAGACTCGTTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAAGGTCTACCTCGCTAACACCACTGCAGCACCT",
    "AGGTGCTGCAGTGGTGTTAGCGAGGTAGACCTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACGTCAACTGACAGTGGTTCGTACTCAGCACCT",
    "AGGTGCTGAGTACGAACCACTGTCAGTTGACGTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAACCCTCCAGGAAAGTACCTCTGATCAGCACCT",
    "AGGTGCTGATCAGAGGTACTTTCCTGGAGGGTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACCAAACCCAACAACCTAGATAGGCCAGCACCT",
    "AGGTGCTGGCCTATCTAGGTTGTTGGGTTTGGTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGTTCCTCGTGCAGTGTCAAGAGATCAGCACCT",
    "AGGTGCTGATCTCTTGACACTGCACGAGGAACTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAATTGCGTCCTGTTACGAGAACTCATCAGCACCT",
    "AGGTGCTGATGAGTTCTCGTAACAGGACGCAATTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGAGCCTCTCATTGTCCGTTCTCTACAGCACCT",
    "AGGTGCTGTAGAGAACGGACAATGAGAGGCTCTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAACCACTGCCATGTATCAAAGTACGCAGCACCT",
    "AGGTGCTGCGTACTTTGATACATGGCAGTGGTTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACTTACTACCCAGTGAACCTCCTCGCAGCACCT",
    "AGGTGCTGCGAGGAGGTTCACTGGGTAGTAAGTTAACCTTAGCAATACGTAACTGAACGAAGT",
    "AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGCATAGTTCTGCATGATGGGTTAGCAGCACCT",
    "AGGTGCTGCTAACCCATCATGCAGAACTATGCTTAACCTTAGCAATACGTAACTGAACGAAGT",
]

import gzip

class ResourceManager:
    """
    Central handler for external genomic resources (FASTA, GTF, BED, Primers).
    """
    def __init__(self):
        self.fasta = None
        self.genes = {}      # Dictionary of gene models from GTF (Legacy)
        self.gene_trees = defaultdict(IntervalTree) # IntervalTrees by chromosome
        self.primers = {}    # Dictionary of primers
        self.bed_regions = [] # List of simple regions from BED

    # --- FASTA HANDLING ---
    def load_fasta(self, filepath):
        """Loads an indexed FASTA file using pysam."""
        try:
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"FASTA file not found: {filepath}")
            
            # pysam.FastaFile requires an indexed fasta (.fai)
            self.fasta = pysam.FastaFile(filepath)
            return True, f"FASTA loaded: {len(self.fasta.references)} contigs."
        except Exception as e:
            return False, f"Error loading FASTA: {e}"

    def get_sequence(self, chrom, start, end):
        """Fetches sequence from loaded FASTA."""
        if self.fasta:
            try:
                return self.fasta.fetch(chrom, start, end).upper()
            except ValueError:
                return None
        return None

    # --- GTF/GFF HANDLING (Gene Models) ---
    def parse_gtf_line(self, line):
        """Helper to parse a single GTF line."""
        if line.startswith("#"): return None
        parts = line.strip().split('\t')
        if len(parts) < 9: return None
        
        chrom = parts[0]
        norm_chrom = chrom.replace("chr", "")
        feature_type = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        attributes = parts[8]
        
        attr_dict = {}
        for attr in attributes.split(';'):
            if not attr.strip(): continue
            try:
                if ' ' in attr.strip():
                    key, val = attr.strip().split(' ', 1)
                    attr_dict[key] = val.replace('"', '')
                elif '=' in attr.strip():
                    key, val = attr.strip().split('=', 1)
                    attr_dict[key] = val
            except ValueError:
                pass
        
        gene_name = attr_dict.get("gene_name") or attr_dict.get("gene_id") or attr_dict.get("Name")
        transcript_id = attr_dict.get("transcript_id", "unknown")
        exon_number = attr_dict.get("exon_number", "?")
        
        if not gene_name: return None
        
        data = {
            "name": gene_name,
            "type": feature_type,
            "id": transcript_id,
            "exon": exon_number
        }
        return norm_chrom, start, end, data

    def load_gtf(self, filepath):
        """
        Parses a GTF file to extract Exon/CDS structures for gene models.
        If a .tbi index exists, returns an IndexedGeneModels object (lazy loading).
        Otherwise, builds IntervalTrees in memory (legacy).
        """
        if IntervalTree is None:
            return False, "intervaltree library not installed. Please run: pip install intervaltree"

        # Check for index
        if filepath.endswith('.gz') and os.path.exists(filepath + '.tbi'):
            try:
                self.gene_trees = IndexedGeneModels(filepath, self)
                return True, f"GTF loaded (Indexed): {filepath}"
            except Exception as e:
                print(f"Failed to load indexed GTF, falling back to full load: {e}")

        try:
            count = 0
            self.genes = defaultdict(lambda: {"transcripts": defaultdict(list)})
            self.gene_trees = defaultdict(IntervalTree)
            
            # Handle gzip
            if filepath.endswith('.gz'):
                f = gzip.open(filepath, 'rt') # rt = read text
            else:
                f = open(filepath, 'r')
                
            with f:
                for line in f:
                    res = self.parse_gtf_line(line)
                    if res:
                        norm_chrom, start, end, data = res
                        self.gene_trees[norm_chrom].addi(start, end + 1, data)
                        count += 1
                        
            return True, f"GTF loaded: {count} features parsed into IntervalTrees."
        except Exception as e:
            return False, f"Error loading GTF: {e}"

    # --- PRIMER HANDLING ---
    def load_primers(self, filepath):
        """Loads primers from a TSV file (name\tsequence)."""
        self.primers = {}
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            # Format: Name [TAB] Sequence
                            name = parts[0].strip()
                            seq = parts[1].strip().upper()
                            self.primers[name] = seq
            return self.primers 
        except Exception as e:
            print(f"Error loading primers: {e}")
            return None

    def load_known_mutations(self, path):
        """
        Loads known pathogenic mutations from TSV.
        """
        known = {}
        try:
            with open(path, 'r') as f:
                header = next(f, None)
                for line in f:
                    if not line.strip(): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        chrom = parts[0]
                        try:
                            start = int(parts[1])
                            bas_parts = parts[3].split('>')
                            if len(bas_parts) == 2:
                                ref = bas_parts[0].strip()
                                alt = bas_parts[1].strip()
                                aa_change = parts[4].strip()
                                known[(chrom, start, ref, alt)] = aa_change
                        except ValueError:
                            continue
            return known
        except Exception as e:
            print(f"Error loading clinical mutations: {e}")
            return {}

    def load_common_snps(self, path):
        """
        Loads common SNPs from a BED file.
        """
        snps = {}
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip(): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        chrom = parts[0]
                        # Standard BED is 0-based. start(0-based) to end(1-based).
                        # Example: chr17 7674108 7674109 -> variant at 7674109 (1-based)
                        start_0 = int(parts[1])
                        pos_1 = start_0 + 1
                        rs_id = parts[3]
                        snps[(chrom, pos_1)] = rs_id
            return snps
        except Exception as e:
            print(f"Error loading SNPs: {e}")
            return {}

    # --- BED HANDLING ---
    def load_simple_bed(self, filepath):
        """Loads a simple BED file for ROI (chrom, start, end, name)."""
        regions = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            chrom = parts[0]
                            start = int(parts[1]) 
                            end = int(parts[2])
                            name = parts[3] if len(parts) > 3 else f"Region_{len(regions)+1}"
                            regions.append({"chrom": chrom, "start": start, "end": end, "name": name})
            self.bed_regions = regions
            return regions 
        except Exception as e:
            print(f"Error loading BED file: {e}")
            return []

# --- HELPERS FOR INDEXED LOADING ---

class IntervalObject:
    """Mimics the object returned by IntervalTree.overlap"""
    def __init__(self, start, end, data):
        self.begin = start
        self.end = end
        self.data = data

class IndexedChromosome:
    """Mimics an IntervalTree for a specific chromosome, using Tabix."""
    def __init__(self, tabix, chrom, parser):
        self.tabix = tabix
        self.chrom = chrom
        self.parser = parser

    def __getstate__(self):
        # Don't pickle the tabix object
        state = self.__dict__.copy()
        if 'tabix' in state:
            del state['tabix']
        # We also need the filepath to reconstruct it, but IndexedGeneModels should handle this
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Tabix will be restored by the parent IndexedGeneModels if needed, 
        # or we assume it's just a proxy. 
        # Actually, IndexedChromosome is usually created on the fly by IndexedGeneModels.__getitem__
        self.tabix = None 

    def overlap(self, start, end):
        results = []
        try:
            iterator = self.tabix.fetch(self.chrom, start, end)
            for line in iterator:
                res = self.parser.parse_gtf_line(line)
                if res:
                    _, f_start, f_end, data = res
                    results.append(IntervalObject(f_start, f_end + 1, data))
        except (ValueError, OSError):
            pass 
        return results

class IndexedGeneModels:
    """Mimics the defaultdict(IntervalTree) structure."""
    def __init__(self, filepath, parser):
        self.filepath = filepath
        self.tabix = pysam.TabixFile(filepath)
        self.parser = parser
        self.contigs = set(self.tabix.contigs)

    def __getstate__(self):
        state = self.__dict__.copy()
        if 'tabix' in state:
            del state['tabix']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if os.path.exists(self.filepath):
            self.tabix = pysam.TabixFile(self.filepath)
        else:
            self.tabix = None
            print(f"Warning: Could not restore TabixFile, path not found: {self.filepath}")
        
    def __contains__(self, key):
        if key in self.contigs: return True
        if f"chr{key}" in self.contigs: return True
        return False
        
    def __getitem__(self, key):
        if not self.tabix:
             if os.path.exists(self.filepath):
                 self.tabix = pysam.TabixFile(self.filepath)
             else:
                 return None

        if key in self.contigs:
            return IndexedChromosome(self.tabix, key, self.parser)
        elif f"chr{key}" in self.contigs:
            return IndexedChromosome(self.tabix, f"chr{key}", self.parser)
        else:
            return IndexedChromosome(self.tabix, key, self.parser)

# Global instance for shared access if needed
resource_manager = ResourceManager()