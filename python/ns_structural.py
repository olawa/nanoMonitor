# Filename: ns_structural.py
# Created: 2025-11-21 13:00 CET

import numpy as np

class GenomeLinearizer:
    """
    Helper class to map (Chrom, Pos) -> Global Linear Position (0..GenomeSize).
    Essential for plotting whole-genome matrices.
    """
    def __init__(self, chrom_lengths):
        self.chrom_lengths = chrom_lengths
        self.chrom_order = self._sort_chromosomes(chrom_lengths.keys())
        self.offsets = {}
        self.total_length = 0
        
        current_offset = 0
        for chrom in self.chrom_order:
            self.offsets[chrom] = current_offset
            current_offset += chrom_lengths[chrom]
        
        self.total_length = current_offset

    def _sort_chromosomes(self, chroms):
        """Sorts chromosomes naturally (1, 2, ... 10, ... X, Y)."""
        def sort_key(c):
            c = c.replace('chr', '')
            if c.isdigit(): return int(c)
            if c == 'X': return 100
            if c == 'Y': return 101
            if c == 'M' or c == 'MT': return 102
            return 200
        return sorted([c for c in chroms if '_' not in c and 'EBV' not in c], key=sort_key)

    def to_global(self, chrom, pos):
        """Converts (chr, pos) to global X coordinate."""
        if chrom not in self.offsets: return None
        return self.offsets[chrom] + pos

    def from_global(self, global_pos):
        """Converts global X coordinate back to (chr, pos)."""
        for chrom in reversed(self.chrom_order):
            if global_pos >= self.offsets[chrom]:
                return chrom, int(global_pos - self.offsets[chrom])
        return None, None

# Known Fusions Database (Example for demo)
KNOWN_FUSIONS = {
    "BCR-ABL1": [("chr22", 23522000), ("chr9", 133738000)],
    "PML-RARA": [("chr15", 74329000), ("chr17", 38492000)],
}

def prepare_matrix_data(sv_links, linearizer):
    """
    Converts raw SV links ((chrA, posA), (chrB, posB)) into plotting coordinates.
    Filters out small local indels (<100kb) to clean up the diagonal.
    """
    x_coords = []
    y_coords = []
    colors = []
    
    for link in sv_links:
        # Handle both 2-element (legacy) and 3-element (with read_id) tuples
        if len(link) == 3:
            (chrA, posA), (chrB, posB), _ = link
        else:
            (chrA, posA), (chrB, posB) = link
        gx = linearizer.to_global(chrA, posA)
        gy = linearizer.to_global(chrB, posB)
        
        if gx is None or gy is None: continue
        
        # Filter local noise (e.g., small deletions) on the diagonal
        if chrA == chrB and abs(posA - posB) < 100000:
            continue
            
        # Ensure x < y for symmetric upper-triangle plotting
        if gx > gy:
            gx, gy = gy, gx
            
        x_coords.append(gx)
        y_coords.append(gy)
        
        # Color logic: Inter-chromosomal = Red, Intra-chromosomal = Blue
        if chrA != chrB:
            colors.append("#E91E63") # Pink/Red for Translocations
        else:
            colors.append("#2196F3") # Blue for large SVs
            
    return x_coords, y_coords, colors