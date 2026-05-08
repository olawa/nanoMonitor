# Filename: ns_visualizer.py
# Created: 2025-11-20 14:00 CET

import pysam
import heapq
from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QPushButton, 
                             QLabel, QLineEdit, QGraphicsView, QGraphicsScene, 
                             QGraphicsRectItem, QWidget, QMessageBox, QComboBox,
                             QCheckBox, QSpinBox)
from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QColor, QPen, QBrush, QPainter, QPixmap, QImageReader
import subprocess
import tempfile

class LayoutEngine:
    """
    Calculates the vertical stacking (lanes) for reads to avoid overlap.
    This is the part that replaces the slow Python logic in bamsnap.
    """
    @staticmethod
    def calculate_lanes_python(reads):
        """
        Optimized Python implementation using a Min-Heap.
        Complexity: O(N log L) where N=reads, L=lanes.
        reads: list of (start, end, read_obj) sorted by start.
        Returns: list of (read_obj, lane_index)
        """
        # Sort by start position
        reads.sort(key=lambda x: x[0])
        
        # Heap stores: (end_position_of_lane, lane_index)
        lane_heap = [] 
        
        results = []
        
        # Track Next Available Lane Index
        next_new_lane = 0
        
        for start, end, read in reads:
            lane_index = -1
            
            # Check if we can fit in an existing lane (where lane_end < current_start)
            if lane_heap and lane_heap[0][0] < start:
                _, lane_index = heapq.heappop(lane_heap)
            else:
                # Create new lane
                lane_index = next_new_lane
                next_new_lane += 1
            
            # Add to results
            results.append((read, lane_index))
            
            # Push back to heap with new end position (plus minimal gap)
            heapq.heappush(lane_heap, (end + 5, lane_index)) # 5bp visual gap
            
        return results, next_new_lane

class ReadViewer(QGraphicsView):
    def __init__(self):
        super().__init__()
        self.scene = QGraphicsScene()
        self.setScene(self.scene)
        self.setRenderHint(QPainter.RenderHint.Antialiasing, False) # False = Faster for rectangles
        self.setBackgroundBrush(QColor("#FFFFFF"))
        self.scale(1, 1) # Initial scale

    def draw_reads(self, placed_reads, total_lanes, region_start, region_end, gene_data=None):
        self.scene.clear()
        
        # Settings
        lane_height = 10
        lane_gap = 2
        
        # Draw gene models at the top first
        gene_height = 0
        if gene_data:
            gene_height = self.draw_gene_models(gene_data, region_start, region_end, y_offset=0)
            gene_height += 20  # Add gap between genes and reads
        
        # Draw Reads below genes
        pen = QPen(Qt.PenStyle.NoPen) # No border is faster
        
        for read, lane in placed_reads:
            # Calculate geometry (offset by gene track height)
            x = read.reference_start
            w = read.reference_length
            y = gene_height + (lane * (lane_height + lane_gap))
            
            # Color based on strand
            color = QColor("#4CAF50") if not read.is_reverse else QColor("#2196F3")
            if read.mapping_quality < 10:
                color.setAlpha(100) # Transparent for low MAPQ
            
            rect = QGraphicsRectItem(x, y, w, lane_height)
            rect.setBrush(QBrush(color))
            rect.setPen(pen)
            
            # Tooltip
            rect.setToolTip(f"{read.query_name}\\nLen: {read.query_length}\\nMAPQ: {read.mapping_quality}")
            
            self.scene.addItem(rect)
        
        # Calculate total height
        reads_height = total_lanes * (lane_height + lane_gap)
        total_height = gene_height + reads_height + 30
            
        # Fit view
        # Ensure the scene rect covers the gene area (y=0 to total_height)
        # and the genomic region (region_start to region_end)
        self.setSceneRect(region_start, 0, region_end - region_start, total_height)
        self.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
    
    def draw_gene_models(self, gene_data, region_start, region_end, y_offset):
        """
        Draw gene models (exons and transcripts) from GTF data.
        gene_data: list of intervals from GTF query
        y_offset: vertical position to start drawing genes
        """
        if not gene_data:
            return 0
        
        # Group by gene name
        from collections import defaultdict
        genes = defaultdict(list)
        for interval in gene_data:
            data = interval.data
            gene_name = data.get("name", "Unknown")
            feature_type = data.get("type", "").lower() # Case-insensitive
            
            # Only draw exons for now
            if feature_type == "exon":
                genes[gene_name].append({
                    "start": interval.begin,
                    "end": interval.end
                })
        
        # Draw each gene
        gene_height = 12 # Thicker exons
        gene_gap = 10
        current_y = y_offset
        
        print(f"DEBUG: Drawing {len(genes)} genes.")
        
        for gene_name, exons in genes.items():
            if not exons:
                continue
            
            # Merge overlapping exons to create a canonical model
            exons.sort(key=lambda x: x["start"])
            merged_exons = []
            if exons:
                curr = exons[0]
                for next_ex in exons[1:]:
                    if next_ex["start"] < curr["end"]: # Overlap
                        curr["end"] = max(curr["end"], next_ex["end"])
                    else:
                        merged_exons.append(curr)
                        curr = next_ex
                merged_exons.append(curr)
            
            # Draw gene line (intron line)
            gene_start = min(e["start"] for e in merged_exons)
            gene_end = max(e["end"] for e in merged_exons)
            
            # Draw thick black line for gene span
            line_pen = QPen(QColor("#000000")) # Black
            line_pen.setWidth(2)
            line = self.scene.addLine(gene_start, current_y + gene_height/2, 
                                      gene_end, current_y + gene_height/2, line_pen)
            
            # Draw exons as thick rectangles
            exon_pen = QPen(QColor("#000000")) # Black border
            exon_pen.setWidth(1)
            exon_brush = QBrush(QColor("#FF9800"))  # Orange for exons
            
            for exon in merged_exons:
                rect = QGraphicsRectItem(exon["start"], current_y, 
                                        exon["end"] - exon["start"], gene_height)
                rect.setBrush(exon_brush)
                rect.setPen(exon_pen)
                rect.setToolTip(f"{gene_name}\\nExon: {exon['start']}-{exon['end']}")
                self.scene.addItem(rect)
            
            # Add gene label
            from PyQt6.QtWidgets import QGraphicsTextItem
            label = QGraphicsTextItem(gene_name)
            label.setPos(gene_start, current_y - 20) # Above gene
            label.setDefaultTextColor(QColor("#000000"))
            from PyQt6.QtGui import QFont
            font = QFont()
            font.setPointSize(10) # Larger font
            font.setBold(True)
            label.setFont(font)
            self.scene.addItem(label)
            
            current_y += gene_height + gene_gap + 15 # Extra space for label
        
        return current_y - y_offset  # Return height used


class BamSnapDialog(QDialog):
    def __init__(self, bam_path, parent=None, gene_models=None, gene_file=None, target_read_ids=None):
        super().__init__(parent)
        self.bam_path = bam_path
        self.gene_models = gene_models
        self.gene_file = gene_file
        self.target_read_ids = target_read_ids
        self.setWindowTitle(f"NanoSnap: {os.path.basename(bam_path)}")
        self.resize(1200, 800)
        
        layout = QVBoxLayout(self)
        
        # Controls
        ctrl_layout = QHBoxLayout()
        self.input_region = QLineEdit()
        self.input_region.setPlaceholderText("chr:start-end (e.g., chr1:10000-20000)")
        
        btn_go = QPushButton("Snap!")
        btn_go.clicked.connect(self.run_rsnap)
        
        btn_zoom_in = QPushButton("+")
        btn_zoom_in.clicked.connect(lambda: self.viewer.scale(1.2, 1.2))
        
        btn_zoom_out = QPushButton("-")
        btn_zoom_out.clicked.connect(lambda: self.viewer.scale(0.8, 0.8))

        ctrl_layout.addWidget(QLabel("Region:"))
        ctrl_layout.addWidget(self.input_region)
        
        # Max Reads Control
        self.chk_limit = QCheckBox("Limit Reads:")
        self.chk_limit.setChecked(True)
        self.chk_limit.stateChanged.connect(lambda s: self.spin_limit.setEnabled(self.chk_limit.isChecked()))
        
        self.spin_limit = QSpinBox()
        self.spin_limit.setRange(1, 100000)
        self.spin_limit.setValue(100)
        self.spin_limit.setToolTip("Maximum number of reads to show (-m)")
        
        ctrl_layout.addWidget(self.chk_limit)
        ctrl_layout.addWidget(self.spin_limit)
        
        # Filter Toggle
        self.chk_filter = QCheckBox("Filter by Selection")
        self.chk_filter.setChecked(True)
        self.chk_filter.setToolTip("Only show reads matching selected amplicon and filters")
        if not self.target_read_ids:
            self.chk_filter.setChecked(False)
            self.chk_filter.setEnabled(False)
            
        ctrl_layout.addWidget(self.chk_filter)
        
        # Rsnap Options
        self.chk_squash = QCheckBox("Squash")
        self.chk_density = QCheckBox("Density")
        
        self.spin_show_ins = QSpinBox()
        self.spin_show_ins.setRange(0, 10000)
        self.spin_show_ins.setValue(0)
        self.spin_show_ins.setToolTip("Minimum insertion length to show (--show-ins)")
        self.spin_show_ins.setSuffix(" bp")
        
        ctrl_layout.addWidget(self.chk_squash)
        ctrl_layout.addWidget(self.chk_density)
        ctrl_layout.addWidget(QLabel("Show Ins:"))
        ctrl_layout.addWidget(self.spin_show_ins)
        
        ctrl_layout.addWidget(btn_go)
        ctrl_layout.addWidget(btn_zoom_in)
        ctrl_layout.addWidget(btn_zoom_out)
        
        layout.addLayout(ctrl_layout)
        
        # Viewer
        self.viewer = ReadViewer()
        self.viewer.setDragMode(QGraphicsView.DragMode.ScrollHandDrag) # Enable panning
        layout.addWidget(self.viewer)
        
        self.snap_count = 0
        
    def run_rsnap(self):
        region_str = self.input_region.text().strip()
        if not region_str: return
        
        # Fix image allocation limit
        QImageReader.setAllocationLimit(0)
        
        temp_bam_path = None
        
        try:
            # Parse region (robust)
            parts = region_str.split(":")
            if len(parts) >= 2:
                # Handle standard chr:start-end or Prefix:chr:start-end
                # We assume the LAST two parts are always chrom and interval if len >= 2
                # e.g. "chr1", "100-200" OR "PrimerA-PrimerB", "chr1", "100-200"
                
                # Careful: The interval part must contain "-"
                if "-" in parts[-1]:
                    chrom = parts[-2]
                    interval = parts[-1].split('(')[0] # Strip (Gene) if present manually
                    start_s, end_s = interval.split("-")
                    start = int(start_s.replace(",",""))
                    end = int(end_s.replace(",",""))
                else:
                    return
            else:
                 return

            # Add padding (10%)
            p_start = max(1, start - int((end-start)*0.1))
            p_end = end + int((end-start)*0.1)      
            padded_region = f"{chrom}:{p_start}-{p_end}"
            
            self.snap_count += 1
            bam_base = os.path.basename(self.bam_path)
            output_file = f"{bam_base}.snap{self.snap_count}.png"
            
            bam_to_use = self.bam_path
            
            # Filter reads if target_read_ids provided AND filter checked
            if self.target_read_ids and self.chk_filter.isChecked():
                print(f"Filtering {len(self.target_read_ids)} reads for visualization...")
                
                # Create temp BAM
                tf = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
                temp_bam_path = tf.name
                tf.close()
                
                count_written = 0
                with pysam.AlignmentFile(self.bam_path, "rb") as infile:
                    with pysam.AlignmentFile(temp_bam_path, "wb", template=infile) as outfile:
                        for read in infile.fetch(chrom, start, end):
                            if read.query_name in self.target_read_ids:
                                outfile.write(read)
                                count_written += 1
                
                print(f"Written {count_written} reads to temp BAM.")
                
                if count_written == 0:
                    QMessageBox.warning(self, "No Reads", "No reads matched the filter criteria in this region.")
                    os.remove(temp_bam_path)
                    return
                
                # Index temp BAM
                pysam.index(temp_bam_path)
                bam_to_use = temp_bam_path
            
            # Construct command
            # Construct command
            # Use --viewer instead of -o
            cmd = ["rsnap", "--viewer", "-b", bam_to_use, "-p", padded_region]
            
            if self.gene_file:
                cmd.extend(["-g", self.gene_file])
                
            if self.chk_limit.isChecked():
                cmd.extend(["-m", str(self.spin_limit.value())])
            
            if self.chk_squash.isChecked():
                cmd.append("--squash")
            if self.chk_density.isChecked():
                cmd.append("--density")
            if self.spin_show_ins.value() > 0:
                cmd.extend(["--show-ins", str(self.spin_show_ins.value())])
                
            print(f"Running rsnap: {' '.join(cmd)}")
            
            # Run blocking so we can clean up temp file afterwards
            subprocess.run(cmd, check=True)
            
            # Viewer handles display, so we don't load anything back into self.viewer
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Error", f"Failed to run rsnap: {str(e)}")
        finally:
            # Cleanup temp files
            if temp_bam_path and os.path.exists(temp_bam_path):
                try:
                    os.remove(temp_bam_path)
                    if os.path.exists(temp_bam_path + ".bai"):
                        os.remove(temp_bam_path + ".bai")
                except Exception:
                    pass

    def load_region(self):
        self.run_rsnap()

import os