# Filename: ns_plotting.py
# Created: 2025-11-21 19:00 CET

from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QDialogButtonBox, QMessageBox, 
                             QFileDialog, QLabel, QHBoxLayout, QCheckBox)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.stats import gaussian_kde
import ns_structural

# --- ACCURACY PLOT DIALOG ---
class AccuracyPlotDialog(QDialog):
    def __init__(self, raw_accuracies, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Read Accuracy Density (KDE)")
        self.resize(800, 600)
        
        layout = QVBoxLayout(self)
        self.canvas = FigureCanvasQTAgg(Figure(figsize=(8, 5)))
        layout.addWidget(self.canvas)
        
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)
        
        self.plot(raw_accuracies)

    def plot(self, accuracies):
        fig = self.canvas.figure
        ax = fig.add_subplot(111)
        ax.clear()
        
        data = np.array([a for a in accuracies if a >= 95 and a <= 100.0])
        
        if len(data) < 5:
            ax.text(0.5, 0.5, "Insufficient data (>95%)", ha='center')
        else:
            kde = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 200)
            ax.plot(x, kde(x), color='#2196F3', lw=2)
            ax.fill_between(x, kde(x), color='#BBDEFB', alpha=0.5)
            
            mean = np.mean(data)
            ax.axvline(mean, color='red', ls='--', label=f'Mean: {mean:.2f}%')
            ax.set_xlabel("Accuracy (%)")
            ax.legend()
            
        self.canvas.draw()

# --- FUSION MATRIX DIALOG (Moved from Main) ---
class FusionMatrixDialog(QDialog):
    def __init__(self, sv_links, chrom_lengths, parent=None, filtered_read_ids=None):
        super().__init__(parent)
        self.setWindowTitle("Genome-Wide Contact Matrix (Structural Variants)")
        self.resize(900, 900)
        self.sv_links = sv_links
        self.filtered_read_ids = filtered_read_ids
        self.linearizer = ns_structural.GenomeLinearizer(chrom_lengths)
        self.parent_gui = parent # Reference to main window for callbacks
        
        layout = QVBoxLayout(self)
        
        # Controls
        ctrl_layout = QHBoxLayout()
        ctrl_layout.addWidget(QLabel("Click points to inspect in Snap View. Red = Translocation, Blue = Intra-chromosomal SV."))
        
        self.chk_filter = QCheckBox("Filter by Selection")
        self.chk_filter.setToolTip("Only show SVs from reads matching current filters (QS, Len, Amplicon)")
        self.chk_filter.setChecked(True if self.filtered_read_ids else False)
        self.chk_filter.setEnabled(True if self.filtered_read_ids else False)
        self.chk_filter.stateChanged.connect(self.draw_matrix)
        
        ctrl_layout.addStretch()
        ctrl_layout.addWidget(self.chk_filter)
        layout.addLayout(ctrl_layout)
        
        self.figure = Figure(figsize=(10, 10))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self) 
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        
        self.draw_matrix()
        self.canvas.mpl_connect('pick_event', self.on_pick)

    def draw_matrix(self):
        ax = self.figure.add_subplot(111)
        ax.clear()
        
        # Filter links
        links_to_plot = self.sv_links
        if self.chk_filter.isChecked() and self.filtered_read_ids:
            # Check if sv_links has read_ids (3-element tuple)
            # If so, filter. If not (legacy), show all? Or warn?
            # We assume updated ns_core.py
            filtered = []
            for link in self.sv_links:
                if len(link) == 3:
                    if link[2] in self.filtered_read_ids:
                        filtered.append(link)
                else:
                    filtered.append(link) # Keep legacy links (no ID to filter)
            links_to_plot = filtered
            
        xs, ys, colors = ns_structural.prepare_matrix_data(links_to_plot, self.linearizer)
        
        if not xs:
            ax.text(0.5, 0.5, "No split reads (SVs) found.", ha='center', transform=ax.transAxes)
            self.canvas.draw()
            return

        # Draw Grid
        ticks = []
        labels = []
        curr = 0
        for chrom in self.linearizer.chrom_order:
            length = self.linearizer.chrom_lengths[chrom]
            ax.axvline(curr, color='#ddd', linewidth=0.5)
            ax.axhline(curr, color='#ddd', linewidth=0.5)
            ticks.append(curr + length/2)
            labels.append(chrom.replace('chr', ''))
            curr += length
            
        ax.set_xlim(0, self.linearizer.total_length)
        ax.set_ylim(0, self.linearizer.total_length)
        
        # Scatter
        self.scatter = ax.scatter(xs, ys, c=colors, s=15, alpha=0.6, edgecolors='none', picker=5)
        
        # Annotations (Known Fusions)
        for name, coords in ns_structural.KNOWN_FUSIONS.items():
            gx = self.linearizer.to_global(coords[0][0], coords[0][1])
            gy = self.linearizer.to_global(coords[1][0], coords[1][1])
            if gx and gy:
                if gx > gy: gx, gy = gy, gx 
                
                # Use Arrow instead of Star covering the point
                # Point to the intersection
                ax.annotate(name, xy=(gx, gy), xytext=(gx + 50000000, gy - 50000000),
                            arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                            fontsize=9, fontweight='bold', color='black')
                
                # Optional: Add a small, transparent star in background
                ax.plot(gx, gy, 'y*', markersize=10, alpha=0.3, zorder=-1)

        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_title(f"Structural Variant Map ({len(xs)} links)", fontsize=14)
        
        self.figure.tight_layout()
        self.canvas.draw()

    def on_pick(self, event):
        ind = event.ind[0] 
        x_data, y_data = self.scatter.get_offsets().T
        global_x = x_data[ind]
        global_y = y_data[ind]
        
        chr1, pos1 = self.linearizer.from_global(global_x)
        chr2, pos2 = self.linearizer.from_global(global_y)
        
        if chr1 and chr2:
            msg = f"Jump to breakpoint?\n{chr1}:{pos1} <--> {chr2}:{pos2}"
            reply = QMessageBox.question(self, "Inspect SV", msg, QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            
            if reply == QMessageBox.StandardButton.Yes:
                region_str = f"{chr2}:{pos2-500}-{pos2+500}"
                # Callback via parent
                if hasattr(self.parent_gui, "open_snap_view_region"):
                     self.parent_gui.open_snap_view_region(region_str)

# --- PDF REPORT ---
def generate_pdf_report(stats, mode, total_reads, parent_widget=None):
    """Generates PDF report."""
    path, _ = QFileDialog.getSaveFileName(parent_widget, "Save Report", "nanostream_report.pdf", "PDF Files (*.pdf)")
    if not path: return

    try:
        with PdfPages(path) as pdf:
            # Title Page
            fig = Figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.axis('off')
            ax.text(0.5, 0.8, "NanoStream Report", ha='center', size=20, weight='bold')
            ax.text(0.5, 0.7, f"Mode: {mode}", ha='center', size=14)
            ax.text(0.5, 0.6, f"Total Reads: {total_reads:,}", ha='center', size=12)
            pdf.savefig(fig)
            
            # Histograms (Amplicon only)
            if mode == "Amplicon":
                sorted_items = sorted(stats.items(), key=lambda x: x[1]['count'], reverse=True)
                for name, data in sorted_items[:20]:
                    if "raw_lengths" in data and len(data["raw_lengths"]) > 10:
                        fig = Figure(figsize=(8,4))
                        ax = fig.add_subplot(111)
                        ax.hist(data["raw_lengths"], bins=50, color='skyblue', edgecolor='black')
                        ax.set_title(f"{name} (n={data['count']})")
                        pdf.savefig(fig)
                        
        QMessageBox.information(parent_widget, "Saved", f"Report saved to {path}")
    except Exception as e:
        QMessageBox.critical(parent_widget, "Error", str(e))