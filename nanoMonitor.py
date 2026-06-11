# Filename: ns_main.py
# Created: 2025-11-21 23:00 CET

import sys
import os
import re
import json
import argparse  
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton, 
                             QFileDialog, QTextEdit, QProgressBar, QLabel, QHBoxLayout, 
                             QGroupBox, QMessageBox, QTabWidget, QTableWidget, 
                             QTableWidgetItem, QHeaderView, QComboBox, QSpinBox, 
                             QStackedWidget, QFormLayout, QDialog, QSplitter, QCheckBox,
                             QListWidget, QListWidgetItem, QLineEdit, QDialogButtonBox,
                             QGridLayout)
from PyQt6.QtCore import Qt, QTimer
import subprocess
import tempfile
import pysam
import threading

# --- PLOTLY IMPORTS ---
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
import plotly.io as pio

import numpy as np
from scipy.stats import gaussian_kde

# Import modules
sys.path.append(os.path.join(os.path.dirname(__file__), "python"))
import ns_core
import ns_workers
import ns_plotting
import ns_visualizer 
import ns_variant
import ns_resources
import ns_structural

# --- HTML TEMPLATE ---
PLOTLY_HTML_TEMPLATE = """
<!DOCTYPE html>
<html><head><script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>html, body { margin: 0; padding: 0; height: 100%; overflow: hidden; }</style></head>
<body><div id="plot" style="width:100%; height:100%;"></div>
<script>
    Plotly.newPlot('plot', [], {}, {responsive: true, displayModeBar: false});
    function updatePlot(data, layout) {
        layout.uirevision = 'true'; 
        Plotly.react('plot', data, layout);
    }
</script></body></html>
"""


DEFAULT_PRIMERS = {
    "BCR-ABL1majorF": "tgaccaactcgtgtgtgaaactc",
    "BCR-ABL1minorF": "accgcatgttccgggacaaaag",
    "BCR-ABL1_R": "tccacttcgtctgagatactggatt",
    "TP53_3kb_A1_F": "GAAGGCTGGTGGCTTCATAG",
    "TP53_3kb_A1_R": "CCAGGGACGAGTGTGGATAC",
    "TP53_D2_F": "AGGCAGAGTTAGGGGGATTG",
    "TP53_D2_R": "CACTCTCAAAGAGGCCAAGG",
    "TP53_I4_F": "AGACAGGTCTGAAGCCTGGA",
    "TP53_I4_R": "GTGGCTGCTCTTCTCTGTCC",
    "FLT3_TKD_3kb_1_F": "GGTGCTTTCACGTTGGTTTT",
    "FLT3_TKD_3kb_1_R": "CGTGCTTCATGCTTGGACTA",
    "TP53_F2_F": "GCCAGGCATTGAAGTCTCAT",
    "TP53_F2_R": "AGGGGATGTTTTGTCAGTGC",
    "TP53_J3_F": "GCTATGATGTTCCTTAGATTAGGTG",
    "TP53_J3_R": "GTTTCTTTGCTGCCGTCTTC",
    "IDH1_1_F": "TGGCCAGGATGAAAGGATAG",
    "IDH1_1_R": "GACAGTGGTCTGGGCAATTT",
    "TP53_3kb_B1_F": "ATCCTGCCACTTTCTGATGG",
    "TP53_3kb_B1_R": "GTCGCATGCACATGTAGTCC",
    "NPM1_5_F": "TTGCTGTTCCATTTGACTGC",
    "NPM1_5_R": "CACCTGACCTCAACCTGGAT",
    "TP53_G4_F": "ACTCGTGAGGCTGCTAGAGG",
    "TP53_G4_R": "GCAGGATTCCTCCAAAATGA",
    "TP53_C10_F": "GAAGGCAGGATGAGAATGGA",
    "TP53_C10_R": "GCAGAGTTAGGGGGATTGCT",
    "TP53_H4_F": "TCACTTCCACGACTGACAGC",
    "TP53_H4_R": "TTCATCTCCCCAGACTCCAC",
    "TP53_E1_F": "TCTACTCCCAACCACCCTTG",
    "TP53_E1_R": "CAGCCATTCTTTTCCTGCTC",
    "IDH2_1_F": "CCCTGGCTTATCCAATCAGA",
    "IDH2_1_R": "CTGCCAACCTCTCTCCAAAG",
    "FLT3_ITD_3kb_1_F": "GGACGAGGATGGAATCAAGA",
    "FLT3_ITD_3kb_1_R": "CTTATTTGCCCTCAGCTTGC"
}

# --- REMOTE FILE DIALOG ---
class RemoteFileDialog(QDialog):
    def __init__(self, server_address, secret=None, parent=None, mode="File"):
        super().__init__(parent)
        self.setWindowTitle(f"Remote {mode} Browser - {server_address}")
        self.resize(600, 400)
        self.server_address = server_address
        self.secret = secret
        self.mode = mode
        self.current_path = "."
        
        import ns_workers
        self.client = ns_workers.RemoteClient(server_address, secret)
        
        self.setup_ui()
        # Delay load to allow dialog to show first
        QTimer.singleShot(100, lambda: self.load_path("."))
        
    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        # Top Bar
        top = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.returnPressed.connect(self.on_path_entered)
        top.addWidget(QLabel("Path:"))
        top.addWidget(self.path_edit)
        
        btn_up = QPushButton("Up")
        btn_up.clicked.connect(self.go_up)
        top.addWidget(btn_up)
        
        btn_go = QPushButton("Go")
        btn_go.clicked.connect(self.on_path_entered)
        top.addWidget(btn_go)
        
        layout.addLayout(top)
        
        # List
        self.list_widget = QListWidget()
        self.list_widget.itemDoubleClicked.connect(self.on_item_double_clicked)
        layout.addWidget(self.list_widget)
        
        # Buttons
        btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Open | QDialogButtonBox.StandardButton.Cancel)
        btns.accepted.connect(self.accept_selection)
        btns.rejected.connect(self.reject)
        
        if self.mode == "Directory":
            btn_choose = QPushButton("Select This Folder")
            btn_choose.clicked.connect(self.accept_current_folder)
            btns.addButton(btn_choose, QDialogButtonBox.ButtonRole.ActionRole)
            
        layout.addWidget(btns)
        
        self.selected_file = None
        
    def load_path(self, path):
        try:
            data = self.client.list_dir(path)
            self.current_path = data.get('current_path', path)
            self.path_edit.setText(self.current_path)
            
            self.list_widget.clear()
            items = data.get('items', [])
            
            for item in items:
                name = item['name']
                icon = "📁" if item['type'] == 'dir' else "📄"
                display = f"{icon} {name}"
                if item['type'] == 'file':
                     size_mb = item['size'] / (1024*1024)
                     display += f" ({size_mb:.2f} MB)"
                
                li = QListWidgetItem(display)
                li.setData(Qt.ItemDataRole.UserRole, item)
                self.list_widget.addItem(li)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to list directory: {e}")

    def on_path_entered(self):
        self.load_path(self.path_edit.text())
        
    def go_up(self):
        parent = os.path.dirname(self.current_path)
        self.load_path(parent)
        
    def on_item_double_clicked(self, item):
        data = item.data(Qt.ItemDataRole.UserRole)
        if data['type'] == 'dir':
            self.load_path(data['path'])
        else:
            self.selected_file = data['path']
            self.accept()
            
    def accept_selection(self):
        item = self.list_widget.currentItem()
        if item:
            data = item.data(Qt.ItemDataRole.UserRole)
            if data['type'] == 'file':
                if self.mode == "File":
                    self.selected_file = data['path']
                    self.accept()
            else:
                self.load_path(data['path'])
                
    def accept_current_folder(self):
        if self.mode == "Directory":
            self.selected_file = self.current_path
            self.accept()

# --- MAIN WINDOW ---
class MainWindow(QMainWindow):
    def __init__(self, args=None):
        super().__init__()
        self.setWindowTitle("NanoStream Monitor v9.1 (FASTQ Support)")
        self.setGeometry(100, 100, 1600, 1000)
        
        # State
        self.mode = "Amplicon"
        self.resources = ns_resources.ResourceManager()
        self.primer_dict = DEFAULT_PRIMERS.copy()
        self.gene_list = None
        self.gene_models = None
        self.monitor_dir = None
        self.file_queue = []
        self.global_stats = {}
        self.session_reads_processed = 0 
        self.current_file_reads_processed = 0 
        # Memory-efficient storage: Use numpy arrays instead of list of dicts
        self.read_qs = np.array([], dtype=np.float32)
        self.read_acc = np.array([], dtype=np.float32)
        self.read_len = np.array([], dtype=np.int32)
        self.read_dx = np.array([], dtype=np.int8)
        self.read_ids = []  # List of read IDs (parallel to above arrays)
        self.read_amplicons = []  # List of amplicon names (parallel to above arrays)
        self.read_concatemers = []  # List of booleans indicating concatemer status (parallel to above arrays)
        self.global_stats = {} # {name: {count, acc, ...}}
        self.current_file_stats = {} # {name: {count, acc, ...}} (Live)
        self.file_stats = {} # {filename: {amplicon_stats}} - Per-file storage
        self.selected_file = "All Files" # Current view selection
        self.file_plot_data = {} # {filename: {qs, acc, len, dx}} - Per-file plot data
        self.sv_links = [] 
        self.chrom_lengths = {} 
        self.current_bam_path = None 
        self.is_monitoring = False
        self.is_processing = False
        self.threads = 8
        self.discovery_mode = False
        self.last_plot_data_len = 0 
        self.max_reads = 1000000 # Default max reads
        self.secret = None
        self.server_address = None
        self.gene_file_path = None
        self.primer_file_path = None
        self.reference_path = None
        self.amplicon_variants = {} # Cache for automated variant results
        self.common_snps = {} # {(chrom, pos): rsID}
        self.amplicon_variants = {} # Cache for automated variant results
        self.common_snps = {} # {(chrom, pos): rsID}
        self.known_mutations = {} # {(chrom, pos, ref, alt): "Name"}
        self.all_seen_files = {} # {filename: full_path} - For "Scan All" feature
        self.pending_file_batches = [] # List of (bam_path, batch_dict) for sequential processing
        
        # Performance: Debounce UI updates
        self.refresh_timer = QTimer()
        self.refresh_timer.setSingleShot(True)
        self.refresh_timer.timeout.connect(lambda: self.refresh_table())
        
        # Auto-load resources
        if os.path.exists("dbSNP_variants.bed"):
            self.common_snps = self.resources.load_common_snps("dbSNP_variants.bed")
        if os.path.exists("clinical_mutations.tsv"):
            self.known_mutations = self.resources.load_known_mutations("clinical_mutations.tsv")
            print(f"DEBUG: Loaded {len(self.known_mutations)} clinical mutations")
            if self.known_mutations:
                sample_keys = list(self.known_mutations.keys())[:3]
                for k in sample_keys:
                    print(f"  -> {k}: {self.known_mutations[k]}")
            
        # Multi-Barcode Support
        self.barcode_stats = {} # {barcode: {amplicon: stats}}
        self.detected_barcodes = set()
        self.current_barcode = "All Barcodes"

        # Automated Variant Calling State
        self.batch_variant_worker = None
        self.variant_queue = {} # {name: region}
        self.scanned_amplicons = set()
        self.loaded_stats_only = False
        self.reset_current_file_buffers()
        
        self.setup_ui()
        self.table.cellDoubleClicked.connect(self.on_table_double_click)
        
        self.plot_timer = QTimer(self)
        self.plot_timer.setInterval(1000) 
        self.plot_timer.timeout.connect(self.update_plots_if_needed)
        self.plot_timer.start()

        if args:
            self.handle_cli_args(args)

    def extract_barcode(self, filepath):
        """Extract barcode from file path."""
        # Check for typical patterns: /barcodeXX/ or /BCXX/ or filename_barcodeXX.fastq
        # Regex for folder or filename
        
        # 1. Check parent folder
        folder = os.path.basename(os.path.dirname(filepath))
        m = re.search(r'(barcode\d+|BC\d+)', folder, re.IGNORECASE)
        if m: return m.group(1).lower()
        
        # 2. Check filename
        filename = os.path.basename(filepath)
        m = re.search(r'(barcode\d+|BC\d+)', filename, re.IGNORECASE)
        if m: return m.group(1).lower()
        
        return "unknown"
        
    def on_barcode_selected(self, bc):
        self.current_barcode = bc
        self.refresh_table()


    def setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        # Base styling to make the dense control layout easier to scan.
        self.setStyleSheet("""
            QMainWindow { background-color: #E3F2FD; }
            QGroupBox {
                font-weight: 600;
                border: 1px solid #c7d8e6;
                border-radius: 6px;
                margin-top: 8px;
                padding-top: 8px;
                background: #edf5fb;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 8px;
                padding: 0 4px;
            }
            QPushButton { min-height: 24px; padding: 3px 8px; }
            QSpinBox, QComboBox, QLineEdit { min-height: 24px; }
            QTabBar::tab { min-height: 24px; padding: 6px 10px; }
            QSplitter::handle { background: #d7e5f0; }
        """)
        main_layout = QVBoxLayout(central)
        
        # --- TOP SPLITTER ---
        top_splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # 1. LEFT PANEL (Controls)
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)
        left_widget.setMaximumWidth(340)

        mode_group = QGroupBox("Session")
        mode_layout = QHBoxLayout()
        self.combo_mode = QComboBox()
        self.combo_mode.addItems(["Amplicon", "RNA-Seq", "DNA"])
        self.combo_mode.currentTextChanged.connect(self.change_mode)
        mode_layout.addWidget(QLabel("Mode:"))
        mode_layout.addWidget(self.combo_mode, 1)
        mode_group.setLayout(mode_layout)
        left_layout.addWidget(mode_group)

        perf_grp = QGroupBox("Performance")
        perf_lay = QGridLayout()
        self.chk_qc_only = QCheckBox("QC Only")
        self.chk_qc_only.setToolTip("Skip primer detection for faster QC")
        self.s_memory_cap = QSpinBox()
        self.s_memory_cap.setRange(100000, 10000000)
        self.s_memory_cap.setValue(1000000)
        self.s_memory_cap.setSingleStep(100000)
        self.s_memory_cap.setSuffix(" reads")
        self.s_memory_cap.setToolTip("Max reads retained in memory during monitoring")
        self.s_memory_cap.valueChanged.connect(self.update_max_reads)
        self.s_tolerance = QSpinBox()
        self.s_tolerance.setRange(0, 500)
        self.s_tolerance.setValue(50)
        self.s_tolerance.setSuffix(" bp")
        self.s_tolerance.setToolTip("Tolerance for fuzzy coordinate matching")
        perf_lay.addWidget(self.chk_qc_only, 0, 0, 1, 2)
        perf_lay.addWidget(QLabel("Memory cap"), 1, 0)
        perf_lay.addWidget(self.s_memory_cap, 1, 1)
        perf_lay.addWidget(QLabel("Primer tol"), 2, 0)
        perf_lay.addWidget(self.s_tolerance, 2, 1)
        perf_grp.setLayout(perf_lay)
        left_layout.addWidget(perf_grp)

        monitor_group = QGroupBox("Resources and Input")
        monitor_layout = QVBoxLayout()
        self.stack = QStackedWidget()

        p_amp = QWidget(); l_amp = QGridLayout(p_amp)
        self.b_load_p = QPushButton("Load Primers"); self.b_load_p.clicked.connect(self.load_primers)
        self.b_load_gtf_amp = QPushButton("Load GTF"); self.b_load_gtf_amp.clicked.connect(self.load_gtf)
        self.l_primers = QLabel("No Primers")
        self.b_load_p.setMaximumWidth(76)
        self.b_load_gtf_amp.setMaximumWidth(76)
        l_amp.addWidget(QLabel("Primers"), 0, 0)
        l_amp.addWidget(self.l_primers, 0, 1)
        l_amp.addWidget(self.b_load_p, 0, 2)
        l_amp.addWidget(QLabel("GTF"), 1, 0)
        l_amp.addWidget(QLabel("Optional"), 1, 1)
        l_amp.addWidget(self.b_load_gtf_amp, 1, 2)

        p_rna = QWidget(); l_rna = QGridLayout(p_rna)
        self.b_load_bed = QPushButton("Load BED"); self.b_load_bed.clicked.connect(self.load_bed)
        self.b_load_gtf_rna = QPushButton("Load GTF"); self.b_load_gtf_rna.clicked.connect(self.load_gtf)
        self.l_genes = QLabel("No Genes")
        self.b_load_bed.setMaximumWidth(76)
        self.b_load_gtf_rna.setMaximumWidth(76)
        l_rna.addWidget(QLabel("BED"), 0, 0)
        l_rna.addWidget(self.l_genes, 0, 1)
        l_rna.addWidget(self.b_load_bed, 0, 2)
        l_rna.addWidget(QLabel("GTF"), 1, 0)
        l_rna.addWidget(QLabel("Optional"), 1, 1)
        l_rna.addWidget(self.b_load_gtf_rna, 1, 2)

        p_dna = QWidget(); l_dna = QVBoxLayout(p_dna)
        self.l_dna_info = QLabel("<b>Plain DNA Mode (QC Only):</b><br>Quality & length control for general sequencing runs.<br>No reference primers or gene BED/GTF files needed.")
        self.l_dna_info.setWordWrap(True)
        self.l_dna_info.setStyleSheet("color: #455a64; font-style: italic;")
        l_dna.addWidget(self.l_dna_info)
        l_dna.addStretch()

        self.stack.addWidget(p_amp)
        self.stack.addWidget(p_rna)
        self.stack.addWidget(p_dna)
        monitor_layout.addWidget(self.stack)

        # Input Selection Buttons Redesign
        input_btn_layout = QHBoxLayout()
        self.b_load_file = QPushButton("📂 Load File")
        self.b_load_file.setStyleSheet("background-color: #BBDEFB; font-weight: bold; min-height: 28px; border-radius: 4px;")
        self.b_load_file.clicked.connect(self.select_single_file)
        
        self.b_monitor_folder = QPushButton("👁️ Monitor Folder")
        self.b_monitor_folder.setStyleSheet("background-color: #C8E6C9; font-weight: bold; min-height: 28px; border-radius: 4px;")
        self.b_monitor_folder.clicked.connect(self.select_monitor_dir)
        
        input_btn_layout.addWidget(self.b_load_file)
        input_btn_layout.addWidget(self.b_monitor_folder)
        monitor_layout.addLayout(input_btn_layout)

        # State Indicator Checkbox (Read-only/Disabled for clear status tracking)
        self.chk_monitor_dir = QCheckBox("Directory Monitoring Active")
        self.chk_monitor_dir.setEnabled(False)
        self.chk_monitor_dir.setStyleSheet("color: #37474F; font-weight: bold;")
        self.chk_monitor_dir.toggled.connect(self.on_monitor_dir_toggled)
        monitor_layout.addWidget(self.chk_monitor_dir)

        self.lbl_dir = QLabel("None")
        self.lbl_dir.setStyleSheet("color: #466;")
        self.b_toggle = QPushButton("Start"); self.b_toggle.setCheckable(True); self.b_toggle.setEnabled(False)
        self.b_toggle.setStyleSheet("font-weight: 700; min-height: 34px;")
        self.b_toggle.clicked.connect(self.toggle_monitor)
        self.b_toggle.setVisible(False)
        monitor_layout.addWidget(QLabel("Selected source"))
        monitor_layout.addWidget(self.lbl_dir)
        self.chk_auto_variant = QCheckBox("Auto-Scan Variants")
        self.chk_auto_variant.setToolTip("Automatically run variant calling on new amplicons during monitoring")
        self.chk_auto_variant.setChecked(False)
        self.chk_primer_analysis = QCheckBox("Primer analysis")
        self.chk_primer_analysis.setChecked(False)
        self.chk_primer_analysis.toggled.connect(self.on_primer_analysis_toggled)
        monitor_layout.addWidget(self.chk_primer_analysis)
        monitor_layout.addWidget(self.chk_auto_variant)
        monitor_layout.addWidget(self.b_toggle)
        monitor_group.setLayout(monitor_layout)
        left_layout.addWidget(monitor_group)

        view_group = QGroupBox("View and Status")
        view_layout = QGridLayout()
        self.progress = QProgressBar(); self.progress.setRange(0, 1000000); self.progress.setValue(0); self.progress.setFormat("0 reads")
        self.l_status = QLabel("Idle")
        self.l_filtered_count = QLabel("Total: 0 | Filtered: 0")
        self.l_filtered_count.setStyleSheet("color: #555; font-weight: bold;")
        self.combo_file_selector = QComboBox()
        self.combo_file_selector.addItem("All Files")
        self.combo_file_selector.currentTextChanged.connect(self.on_file_selected)
        self.combo_file_selector.setStyleSheet("font-weight: 600;")
        self.combo_barcode = QComboBox()
        self.combo_barcode.addItem("All Barcodes")
        self.combo_barcode.currentTextChanged.connect(self.on_barcode_selected)
        self.combo_barcode.setStyleSheet("font-weight: 600; min-width: 100px;")
        self.b_clear = QPushButton("Clear")
        self.b_clear.setStyleSheet("background-color: #E0E0E0; font-weight: bold; min-height: 28px;")
        self.b_clear.clicked.connect(self.on_clear_clicked)
        self.b_snap = QPushButton("Snap")
        self.b_snap.clicked.connect(self.open_snap_view)
        self.b_snap.setEnabled(False)
        self.b_snap.setStyleSheet("background-color: #E1BEE7; font-weight: bold;")
        self.b_variant = QPushButton("Variants"); self.b_variant.clicked.connect(self.run_variant_calling); self.b_variant.setEnabled(False)
        self.b_scan_all = QPushButton("Scan All"); self.b_scan_all.clicked.connect(self.scan_all_files); self.b_scan_all.setEnabled(False)
        self.b_scan_all.setToolTip("Run variant calling on ALL monitoring files")
        self.b_matrix = QPushButton("Matrix"); self.b_matrix.clicked.connect(self.open_fusion_matrix); self.b_matrix.setEnabled(False)
        view_layout.addWidget(self.l_status, 0, 0, 1, 2)
        view_layout.addWidget(self.l_filtered_count, 1, 0, 1, 2)
        view_layout.addWidget(self.progress, 2, 0, 1, 2)
        view_layout.addWidget(QLabel("File"), 3, 0)
        view_layout.addWidget(self.combo_file_selector, 3, 1)
        view_layout.addWidget(QLabel("Barcode"), 4, 0)
        view_layout.addWidget(self.combo_barcode, 4, 1)
        view_layout.addWidget(self.b_snap, 5, 0)
        view_layout.addWidget(self.b_clear, 5, 1)
        view_layout.addWidget(self.b_variant, 6, 0)
        view_layout.addWidget(self.b_scan_all, 6, 1)
        view_layout.addWidget(self.b_matrix, 7, 0, 1, 2)
        view_group.setLayout(view_layout)
        left_layout.addWidget(view_group)
        left_layout.addStretch()
        
        top_splitter.addWidget(left_widget)
        
        # 2. RIGHT PANEL (Table/Log)
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(8)

        f_grp = QGroupBox("Interactive Filters")
        f_lay = QGridLayout()
        self.s_qs = QSpinBox(); self.s_qs.setValue(10); self.s_qs.setRange(0, 60)
        self.s_len = QSpinBox(); self.s_len.setRange(0, 50000)
        self.s_max_len = QSpinBox()
        self.s_max_len.setRange(0, 500000)
        self.s_max_len.setValue(0)
        self.s_max_len.setSpecialValueText("None")
        self.s_max_len.setToolTip("Maximum read length limit (0 = none)")
        
        self.s_filter_max_reads = QSpinBox()
        self.s_filter_max_reads.setRange(0, 1000000)
        self.s_filter_max_reads.setValue(0)
        self.s_filter_max_reads.setSpecialValueText("All")
        self.s_filter_max_reads.setToolTip("Limit selected/exported reads per amplicon (0 = all)")
        self.s_qs.valueChanged.connect(self.force_update_plots) 
        self.s_len.valueChanged.connect(self.force_update_plots)
        self.s_max_len.valueChanged.connect(self.force_update_plots)
        self.chk_duplex = QCheckBox("Duplex Only")
        self.chk_duplex.stateChanged.connect(self.force_update_plots)
        
        self.chk_use_rust = QCheckBox("Use Rust (nanostream)")
        self.chk_use_rust.setToolTip("Use fast Rust backend for primer matching (requires nanostream binary)")
        self.chk_use_rust.setChecked(True)  # Default to Rust if available
        
        self.b_run_duplex = QPushButton("Run Duplex Discovery")
        self.b_run_duplex.setToolTip("Re-scan file with current filters to find duplex pairs")
        self.b_recalc = QPushButton("Recalculate Table")
        self.b_recalc.setToolTip("Update table with current filters (QS, Len, Max Len)")
        self.b_recalc.clicked.connect(self.recalculate_table)
        
        self.b_auto_variant = QPushButton("Auto-Variant")
        self.b_auto_variant.setToolTip("Run variant calling for top amplicons")
        self.b_auto_variant.clicked.connect(self.run_batch_variant_calling)
 
        self.b_rsnap_var = QPushButton("rsnap Variant")
        self.b_rsnap_var.setToolTip("Run rsnap variant caller on selected amplicon")
        self.b_rsnap_var.clicked.connect(self.run_rsnap_variants)
        
        self.b_run_duplex.clicked.connect(self.run_duplex_discovery)
 
        f_lay.addWidget(QLabel("Min QS"), 0, 0)
        f_lay.addWidget(self.s_qs, 0, 1)
        f_lay.addWidget(QLabel("Min Len"), 0, 2)
        f_lay.addWidget(self.s_len, 0, 3)
        f_lay.addWidget(QLabel("Max Len"), 0, 4)
        f_lay.addWidget(self.s_max_len, 0, 5)
        f_lay.addWidget(QLabel("Max Reads"), 0, 6)
        f_lay.addWidget(self.s_filter_max_reads, 0, 7)
        f_lay.addWidget(self.chk_duplex, 1, 0, 1, 2)
        f_lay.addWidget(self.chk_use_rust, 1, 2, 1, 2)
        f_lay.addWidget(self.b_recalc, 1, 4, 1, 2)
        f_grp.setLayout(f_lay)
        right_layout.addWidget(f_grp)

        analysis_grp = QGroupBox("Advanced Actions")
        analysis_lay = QHBoxLayout()
        analysis_lay.addWidget(self.b_auto_variant)
        analysis_lay.addWidget(self.b_rsnap_var)
        analysis_lay.addWidget(self.b_run_duplex)
        analysis_lay.addStretch()
        analysis_grp.setLayout(analysis_lay)
        right_layout.addWidget(analysis_grp)
        
        self.tabs = QTabWidget()
        self.table = QTableWidget()
        
        # Enable multi-row selection
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.MultiSelection)
        self.table.itemSelectionChanged.connect(self.update_selection_count)
        
        # Configure column resize modes
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch) # Name gets most space
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents) # Count
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents) # Med Len
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents) # SD Len
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents) # Avg QS
        
        # Add selection counter and export button below table
        selection_widget = QWidget()
        selection_layout = QHBoxLayout(selection_widget)
        selection_layout.setContentsMargins(5, 2, 5, 2)
        self.l_selected_count = QLabel("Selected: 0 amplicons, 0 reads")
        self.b_export = QPushButton("Export Selected")
        self.b_export.clicked.connect(self.export_selected_reads)
        self.b_export.setEnabled(False)
        self.b_export.setStyleSheet("height: 30px; font-weight: bold;")
        selection_layout.addWidget(self.l_selected_count)
        selection_layout.addStretch()
        selection_layout.addWidget(self.b_export)
        
        # Add table and selection controls to tab
        table_container = QWidget()
        table_container_layout = QVBoxLayout(table_container)
        table_container_layout.setContentsMargins(0, 0, 0, 0)
        table_container_layout.addWidget(self.table)
        table_container_layout.addWidget(selection_widget)
        
        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        self.tabs.addTab(table_container, "Results")
        self.tabs.addTab(self.log_view, "Log")
        # self.tabs.setMaximumHeight(350) # Removed to allow expansion 
        right_layout.addWidget(self.tabs)
        top_splitter.addWidget(right_widget)
        top_splitter.setStretchFactor(0, 1); top_splitter.setStretchFactor(1, 2)
        
        # --- MAIN SPLITTER (Top vs Bottom) ---
        main_splitter = QSplitter(Qt.Orientation.Vertical)
        main_splitter.addWidget(top_splitter)
        
        # 3. BOTTOM SECTION (Plotly Charts Side-by-Side)
        bottom_widget = QWidget()
        bottom_layout = QHBoxLayout(bottom_widget)
        bottom_layout.setContentsMargins(0, 5, 0, 0)
        
        # --- PLOTS ---
        self.web_acc = QWebEngineView()
        self.web_acc.setHtml(PLOTLY_HTML_TEMPLATE)
        
        self.web_qs = QWebEngineView()
        self.web_qs.setHtml(PLOTLY_HTML_TEMPLATE)
        
        self.web_hist = QWebEngineView()
        self.web_hist.setHtml(PLOTLY_HTML_TEMPLATE)
        
        # Wrap histogram in a layout with bases/counts toggle
        hist_container = QWidget()
        hist_vlay = QVBoxLayout(hist_container)
        hist_vlay.setContentsMargins(0, 0, 0, 0)
        hist_vlay.setSpacing(2)
        
        hist_header = QHBoxLayout()
        self.chk_hist_bases = QCheckBox("Show bases per bar")
        self.chk_hist_bases.setToolTip("Toggle between read counts and total bases per length bin")
        self.chk_hist_bases.stateChanged.connect(self.force_update_plots)
        self.chk_hist_bases.setStyleSheet("font-weight: bold; color: #311B92;")
        
        hist_header.addWidget(QLabel("<b>Fragment Length Distribution</b>"))
        hist_header.addStretch()
        hist_header.addWidget(self.chk_hist_bases)
        
        hist_vlay.addLayout(hist_header)
        hist_vlay.addWidget(self.web_hist)
        
        # Top Row: Acc, QS, Hist (All Horizontal)
        plot_splitter = QSplitter(Qt.Orientation.Horizontal)
        plot_splitter.addWidget(self.web_acc)
        plot_splitter.addWidget(self.web_qs)
        plot_splitter.addWidget(hist_container)
        
        # Sizes: Acc=0.75, QS=1, Hist=2 (Approx)
        plot_splitter.setSizes([225, 300, 675])
        
        bottom_layout.addWidget(plot_splitter)
        
        main_splitter.addWidget(bottom_widget)
        main_splitter.setSizes([500, 400])
        main_layout.addWidget(main_splitter)

        self.watcher_thread = None
        self.worker_thread = None
        self.variant_thread = None 

    def reset_current_file_buffers(self):
        self.current_file_qs = np.array([], dtype=np.float32)
        self.current_file_acc = np.array([], dtype=np.float32)
        self.current_file_len = np.array([], dtype=np.int32)
        self.current_file_dx = np.array([], dtype=np.int8)
        self.current_file_amplicons = []
        self.current_file_ids = []

    # --- PLOT UPDATES ---
    def update_plots_if_needed(self):
        if hasattr(self, 'update_clear_button_state'):
            self.update_clear_button_state()
        curr_len = len(self.read_qs)
        if curr_len > self.last_plot_data_len:
            self.update_plots()
            self.last_plot_data_len = curr_len
            
    def force_update_plots(self):
        self.update_plots()

    def update_plots(self):
        qs, acc, lengths = self.get_filtered_data()
        self.update_acc_plot_js(acc)
        self.update_qs_plot_js(qs)
        self.update_hist_plot_js(lengths)

    def get_filtered_data(self):
        # Get data based on selected file
        if self.selected_file == "All Files":
            # Use global arrays (all files combined)
            qs_data = self.read_qs
            acc_data = self.read_acc
            len_data = self.read_len
            dx_data = self.read_dx
        else:
            # Use selected file's data
            file_data = self.file_plot_data.get(self.selected_file, {})
            qs_data = file_data.get("qs", np.array([]))
            acc_data = file_data.get("acc", np.array([]))
            len_data = file_data.get("len", np.array([]))
            dx_data = file_data.get("dx", np.array([]))
        
        if len(qs_data) == 0: 
            self.l_filtered_count.setText("Total: 0 | Filtered: 0")
            return np.array([]), np.array([]), np.array([])
            
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        max_len = self.s_max_len.value()
        duplex_only = self.chk_duplex.isChecked()
        
        # Create filter mask
        mask = (qs_data >= min_qs) & (len_data >= min_len)
        if max_len > 0:
            mask = mask & (len_data <= max_len)
        if duplex_only:
            mask = mask & (dx_data == 1)
        
        filtered_qs = qs_data[mask]
        filtered_acc = acc_data[mask]
        filtered_len = len_data[mask]
        
        total = len(qs_data)
        filtered_count = len(filtered_qs)
        
        self.l_filtered_count.setText(f"Total: {total:,} | Filtered: {filtered_count:,}")
        
        return filtered_qs, filtered_acc, filtered_len

    def update_acc_plot_js(self, data):
        # Handle Clear (Empty Data)
        if len(data) == 0:
            self.web_acc.page().runJavaScript("Plotly.newPlot('plot', [], {}, {responsive: true, displayModeBar: false});")
            return

        if len(data) < 10: return
        
        # Accuracy Plot (Always 95-100 or 0-100)
        x_range = [95, 100]
        x_title = "Accuracy (%)"
        visible_data = data[(data >= 0) & (data <= 100)]
        if len(visible_data) < 10: return

        try:
            # Check for singular data (all values same)
            if np.std(visible_data) < 1e-6:
                # Singular data: Plot a single vertical line
                mode_val = float(np.mean(visible_data))
                x_vals = [mode_val, mode_val]
                y_vals = [0, 1] # Arbitrary height
                trace = {
                    "x": x_vals, "y": y_vals,
                    "mode": 'lines',
                    "line": {"color": '#4CAF50', "width": 4},
                    "name": 'Density'
                }
            else:
                kde = gaussian_kde(visible_data)
                x_vals = np.linspace(visible_data.min(), visible_data.max(), 200)
                y_vals = kde(x_vals)
                peak_idx = np.argmax(y_vals)
                mode_val = float(x_vals[peak_idx])
                max_y = float(y_vals.max()) if hasattr(y_vals, 'max') else 1.0
                
                trace = {
                    "x": x_vals.tolist(), "y": y_vals.tolist(),
                    "fill": 'tozeroy', "mode": 'lines',
                    "line": {"color": '#4CAF50', "width": 2},
                    "name": 'Density'
                }
                
            max_y_val = float(y_vals.max()) if hasattr(y_vals, 'max') else 1.0
            if isinstance(y_vals, list): max_y_val = 1.0 # Fallback for singular case list

            layout = {
                "margin": {"l": 40, "r": 20, "t": 30, "b": 30},
                "xaxis": {"title": x_title, "range": x_range},
                "yaxis": {"title": "Density"},
                "showlegend": False,
                "shapes": [{"type": "line", "x0": mode_val, "x1": mode_val, "y0": 0, "y1": max_y_val, "line": {"color": "red", "width": 2, "dash": "dash"}}],
                "annotations": [{"x": mode_val, "y": max_y_val, "text": f"Mode: {mode_val:.2f}", "showarrow": True, "arrowhead": 1, "ax": 0, "ay": -40, "font": {"color": "red"}}]
            }
            self.web_acc.page().runJavaScript(f"updatePlot([ {json.dumps(trace)} ], {json.dumps(layout)});")
        except Exception as e:
            print(f"Acc Plot Error: {e}")
            pass

    def update_qs_plot_js(self, data):
        # Handle Clear
        if len(data) == 0:
            self.web_qs.page().runJavaScript("Plotly.newPlot('plot', [], {}, {responsive: true, displayModeBar: false});")
            return

        if len(data) < 10: return
        
        # Q-Score Plot (0-60)
        x_range = [0, 60]
        x_title = "Q-Score (Phred)"
        visible_data = data[(data >= 0) & (data <= 60)]
        if len(visible_data) < 10: return

        try:
            # Check for singular data
            if np.std(visible_data) < 1e-6:
                mode_val = float(np.mean(visible_data))
                x_vals = [mode_val, mode_val]
                y_vals = [0, 1]
                trace = {
                    "x": x_vals, "y": y_vals,
                    "mode": 'lines',
                    "line": {"color": '#2196F3', "width": 4},
                    "name": 'Density'
                }
            else:
                kde = gaussian_kde(visible_data)
                x_vals = np.linspace(visible_data.min(), visible_data.max(), 200)
                y_vals = kde(x_vals)
                peak_idx = np.argmax(y_vals)
                mode_val = float(x_vals[peak_idx])
                
                trace = {
                    "x": x_vals.tolist(), "y": y_vals.tolist(),
                    "fill": 'tozeroy', "mode": 'lines',
                    "line": {"color": '#2196F3', "width": 2},
                    "name": 'Density'
                }
                
            max_y_val = float(y_vals.max()) if hasattr(y_vals, 'max') else 1.0
            if isinstance(y_vals, list): max_y_val = 1.0

            layout = {
                "margin": {"l": 40, "r": 20, "t": 30, "b": 30},
                "xaxis": {"title": x_title, "range": x_range},
                "yaxis": {"title": "Density"},
                "showlegend": False,
                "shapes": [{"type": "line", "x0": mode_val, "x1": mode_val, "y0": 0, "y1": max_y_val, "line": {"color": "red", "width": 2, "dash": "dash"}}],
                "annotations": [{"x": mode_val, "y": max_y_val, "text": f"Mode: {mode_val:.2f}", "showarrow": True, "arrowhead": 1, "ax": 0, "ay": -40, "font": {"color": "red"}}]
            }
            self.web_qs.page().runJavaScript(f"updatePlot([ {json.dumps(trace)} ], {json.dumps(layout)});")
        except Exception as e:
            print(f"QS Plot Error: {e}")
            pass

    def update_hist_plot_js(self, data):
        # Handle Clear
        if len(data) == 0:
            self.web_hist.page().runJavaScript("Plotly.newPlot('plot', [], {}, {responsive: true, displayModeBar: false});")
            return

        visible = data[data < 50000]
        if len(visible) == 0:
            visible = data
        if len(visible) == 0: return
        median_len = np.median(visible)
        
        # Determine max length for histogram range and bins
        filter_max_len = self.s_max_len.value() if hasattr(self, 's_max_len') else 0
        if filter_max_len > 0:
            max_len = filter_max_len
        elif len(data) > 0:
            max_len = int(np.max(data))
        else:
            max_len = 6000
            
        max_len = max(max_len, 1000)
        
        # Dynamic bin size: targeting around 60 bins
        bin_size = max(10, int(max_len / 60))
        # Round to clean intervals
        if bin_size > 1000:
            bin_size = (bin_size // 1000) * 1000
        elif bin_size > 500:
            bin_size = (bin_size // 500) * 500
        elif bin_size > 100:
            bin_size = (bin_size // 100) * 100
        elif bin_size > 50:
            bin_size = (bin_size // 50) * 50
        elif bin_size > 10:
            bin_size = (bin_size // 10) * 10
        bin_size = max(bin_size, 10)
        
        # Check toggle state for Show Bases
        show_bases = self.chk_hist_bases.isChecked() if hasattr(self, 'chk_hist_bases') else False
        
        if show_bases:
            trace = {
                "x": data.tolist(), 
                "y": data.tolist(),
                "histfunc": "sum",
                "type": "histogram", 
                "xbins": {"start": 0, "end": max_len, "size": bin_size},
                "marker": {"color": '#673AB7', "opacity": 0.75}, 
                "name": 'Bases'
            }
            yaxis_title = "Total Bases (bp)"
        else:
            trace = {
                "x": data.tolist(), 
                "type": "histogram", 
                "xbins": {"start": 0, "end": max_len, "size": bin_size},
                "marker": {"color": '#673AB7', "opacity": 0.75}, 
                "name": 'Reads'
            }
            yaxis_title = "Count"
        
        shapes = [{"type": "line", "x0": median_len, "x1": median_len, "y0": 0, "y1": 1, "yref": "paper", "line": {"color": "orange", "width": 2, "dash": "dash"}}]
        annotations = [{"x": median_len, "y": 1, "yref": "paper", "text": f"Med: {int(median_len)}", "showarrow": False, "yshift": 10, "font": {"color": "orange"}}]

        if self.mode == "Amplicon" and self.global_stats:
            sorted_amps = sorted(self.global_stats.items(), key=lambda x: x[1]['count'], reverse=True)
            for i, (name, stats) in enumerate(sorted_amps[:10]):
                med_len = stats.get("median_len", 0)
                if med_len > 0:
                    short_name = name.split(' ')[0][:10] + "..." if len(name) > 15 else name
                    y_pos = 0.9 - (i * 0.05)
                    shapes.append({
                        "type": "line", "x0": med_len, "x1": med_len,
                        "y0": 0, "y1": y_pos, "yref": "paper",
                        "line": {"color": "gray", "width": 1, "dash": "dot"}
                    })
                    annotations.append({
                        "x": med_len, "y": y_pos, "yref": "paper",
                        "text": short_name, "xanchor": "left", "showarrow": False,
                        "font": {"size": 10, "color": "gray"}
                    })

        layout = {
            "margin": {"l": 50, "r": 20, "t": 30, "b": 40},
            "xaxis": {"title": "Length (bp)", "range": [0, max_len]}, 
            "yaxis": {"title": yaxis_title},
            "bargap": 0.1,
            "hovermode": "x",
            "shapes": shapes,
            "annotations": annotations
        }
        traces = [trace]
        
        # Highlight selected amplicons
        if hasattr(self, 'selected_amplicons') and self.selected_amplicons:
            selected_indices = []
            for i, amp in enumerate(self.read_amplicons):
                if amp in self.selected_amplicons:
                    selected_indices.append(i)
            
            if selected_indices:
                selected_lengths = self.read_len[selected_indices]
                
                trace_selected = {
                    "x": selected_lengths.tolist(),
                    "type": "histogram",
                    "xbins": {"start": 0, "end": max_len, "size": bin_size},
                    "marker": {"color": 'red', "opacity": 0.6},
                    "name": 'Selected'
                }
                traces.append(trace_selected)
        
        # Overlay concatemers as black bars
        if hasattr(self, 'read_concatemers') and any(self.read_concatemers):
            concatemer_indices = [i for i, is_concat in enumerate(self.read_concatemers) if is_concat]
            
            if concatemer_indices:
                concatemer_lengths = self.read_len[concatemer_indices]
                
                trace_concatemers = {
                    "x": concatemer_lengths.tolist(),
                    "type": "histogram",
                    "xbins": {"start": 0, "end": max_len, "size": bin_size},
                    "marker": {"color": 'black', "opacity": 0.7},
                    "name": 'Concatemers'
                }
                traces.append(trace_concatemers)

        layout = {
            "margin": {"l": 50, "r": 20, "t": 30, "b": 40},
            "xaxis": {"title": "Length (bp)", "range": [0, max_len]}, 
            "yaxis": {"title": yaxis_title},
            "bargap": 0.1,
            "hovermode": "x",
            "barmode": "overlay", # Overlay histograms
            "shapes": shapes,
            "annotations": annotations
        }
        self.web_hist.page().runJavaScript(f"updatePlot({json.dumps(traces)}, {json.dumps(layout)});")

    # --- DATA HANDLER ---
    def update_live_data(self, partial):
        if "metadata" in partial:
            # Extract only numeric fields
            new_qs = np.array([r['qs'] for r in partial["metadata"]], dtype=np.float32)
            new_acc = np.array([r['acc'] for r in partial["metadata"]], dtype=np.float32)
            new_len = np.array([r['len'] for r in partial["metadata"]], dtype=np.int32)
            new_dx = np.array([r.get('dx', 0) for r in partial["metadata"]], dtype=np.int8)
            new_ids = [r['id'] for r in partial["metadata"]]  # Extract read IDs
            
            # Append to global arrays (all files combined)
            self.read_qs = np.concatenate([self.read_qs, new_qs])
            self.read_acc = np.concatenate([self.read_acc, new_acc])
            self.read_len = np.concatenate([self.read_len, new_len])
            self.read_dx = np.concatenate([self.read_dx, new_dx])
            self.read_ids.extend(new_ids)
            
            # Initialize amplicon assignments as "Unknown" for now
            self.read_amplicons.extend(["Unknown"] * len(new_ids))
            
            # Initialize concatemer flags as False for now
            self.read_concatemers.extend([False] * len(new_ids))
            
            # Also append to current file arrays
            self.current_file_qs = np.concatenate([self.current_file_qs, new_qs])
            self.current_file_acc = np.concatenate([self.current_file_acc, new_acc])
            self.current_file_len = np.concatenate([self.current_file_len, new_len])
            self.current_file_dx = np.concatenate([self.current_file_dx, new_dx])
            
            if hasattr(self, 'current_file_ids'):
                self.current_file_ids.extend(new_ids)
            if hasattr(self, 'current_file_amplicons'):
                self.current_file_amplicons.extend(["Unknown"] * len(new_ids))
                
            # Rolling Window (Memory Optimization)
            # Check if we exceeded max_reads (only during live directory monitoring)
            current_total = len(self.read_qs)
            if self.is_monitoring and current_total > self.max_reads:
                excess = current_total - self.max_reads
                # Slice arrays to keep last N
                self.read_qs = self.read_qs[excess:]
                self.read_acc = self.read_acc[excess:]
                self.read_len = self.read_len[excess:]
                self.read_dx = self.read_dx[excess:]
                self.read_ids = self.read_ids[excess:]
                self.read_amplicons = self.read_amplicons[excess:]
                self.read_concatemers = self.read_concatemers[excess:]
                
                # Also trim current file arrays if needed (though usually we want to keep current file complete?)
                # If we are monitoring a directory, "current file" might be small, but "global" grows.
                # If we are in single file mode, current == global.
                # Let's apply limit to global arrays primarily for visualization.
                # But we should probably also trim current file arrays to avoid memory leak if single file is huge.
                
                curr_file_total = len(self.current_file_qs)
                if curr_file_total > self.max_reads:
                    excess_file = curr_file_total - self.max_reads
                    self.current_file_qs = self.current_file_qs[excess_file:]
                    self.current_file_acc = self.current_file_acc[excess_file:]
                    self.current_file_len = self.current_file_len[excess_file:]
                    self.current_file_dx = self.current_file_dx[excess_file:]
                    if hasattr(self, 'current_file_ids'):
                        self.current_file_ids = self.current_file_ids[excess_file:]
                    if hasattr(self, 'current_file_amplicons'):
                        self.current_file_amplicons = self.current_file_amplicons[excess_file:]
        
        # Process amplicon assignments
        if "read_amplicon_map" in partial:
            self.update_amplicon_assignments(partial["read_amplicon_map"])
            
        # Handle Live Amplicon Stats
        if "amplicons" in partial:
            # print(f"DEBUG: Received {len(partial['amplicons'])} amplicons from worker.")
            for name, data in partial["amplicons"].items():
                if name not in self.global_stats: 
                    self.global_stats[name] = {"count":0, "acc":0, "median_len":0, "stdev_len":0, "raw_lengths":[]}
                
                # Update counts (Note: data contains TOTAL counts from worker snapshot)
                # Wait, worker sends snapshot of TOTAL stats? 
                # ns_amplicon.py: snapshot = {k: v.copy() ...}
                # Yes, it sends the ACCUMULATED stats.
                # So we should REPLACE global_stats for this name, not add?
                # But wait, if we are processing multiple files (Directory Monitor), global_stats should accumulate across files.
                # But within one file, the worker sends accumulated stats for THAT file.
                # So we need to be careful.
                # If we just replace, we lose previous files' data.
                # If we add, we double count.
                
                # Solution:
                # The worker sends stats for the CURRENT file.
                # We need to track stats for the current file separately and merge into global?
                # Or, since we process files sequentially:
                # When a file starts, we could have a "current_file_stats".
                # When it finishes, we merge into "global_stats".
                # But for live view, we want to see global + current.
                
                # Actually, let's look at how on_results handles it.
                # on_results adds: self.global_stats[name]["count"] += data["count"]
                # This implies on_results receives the FINAL delta or total?
                # ns_amplicon.py returns final_stats at the end.
                
                # If partial sends TOTAL stats for the current file so far:
                # We should use a separate dict for "current_file_stats".
                # And refresh_table should display global + current.
                
                # Let's add self.current_file_stats = {}
                # And reset it in process_next.
                
                self.current_file_stats[name] = data
            
            self.refresh_table()
    
    def update_amplicon_assignments(self, read_amplicon_map):
        """Update amplicon assignments for reads based on mapping from ns_amplicon."""
        if not read_amplicon_map:
            return
        # Create a reverse lookup for efficient matching
        # Build index of read_ids to array positions
        read_id_index = {rid: i for i, rid in enumerate(self.read_ids)}
        
        # Also build index for current file if available
        current_id_index = {}
        if hasattr(self, 'current_file_ids'):
            current_id_index = {rid: i for i, rid in enumerate(self.current_file_ids)}
        
        # Update amplicon assignments
        count_updated = 0
        for read_id, amplicon_name in read_amplicon_map.items():
            # Check if this read is a concatemer (marked with |CONCAT suffix)
            is_concatemer = "|CONCAT" in amplicon_name
            clean_amplicon_name = amplicon_name.replace("|CONCAT", "") if is_concatemer else amplicon_name
            
            if read_id in read_id_index:
                idx = read_id_index[read_id]
                self.read_amplicons[idx] = clean_amplicon_name
                self.read_concatemers[idx] = is_concatemer
                count_updated += 1
            
            if read_id in current_id_index and hasattr(self, 'current_file_amplicons'):
                c_idx = current_id_index[read_id]
                self.current_file_amplicons[c_idx] = clean_amplicon_name
    
    def calculate_amplicon_stats_from_metadata(self, amplicons=None, qs=None, acc=None, lengths=None):
        """Calculate amplicon statistics from metadata arrays, including QS."""
        if amplicons is None: amplicons = self.read_amplicons
        if qs is None: qs = self.read_qs
        if acc is None: acc = self.read_acc
        if lengths is None: lengths = self.read_len
        
        amplicon_stats = {}
        
        # Group data by amplicon
        for i in range(len(amplicons)):
            amp_name = amplicons[i]
            if amp_name == "Unknown":
                continue
            
            # Ensure index is within bounds for all arrays
            if i >= len(qs) or i >= len(acc) or i >= len(lengths):
                continue
            
            if amp_name not in amplicon_stats:
                amplicon_stats[amp_name] = {
                    "count": 0,
                    "qs_values": [],
                    "acc_values": [],
                    "lengths": []
                }
            
            amplicon_stats[amp_name]["count"] += 1
            amplicon_stats[amp_name]["qs_values"].append(float(qs[i]))
            amplicon_stats[amp_name]["acc_values"].append(float(acc[i]))
            amplicon_stats[amp_name]["lengths"].append(int(lengths[i]))
        
        # Calculate averages and add to formatted stats
        final_stats = {}
        for amp_name, data in amplicon_stats.items():
            chrom, start, end = self.parse_genomic_region(amp_name)
            region_str = f"{chrom}:{start}-{end}" if chrom else None
            
            final_stats[amp_name] = {
                "count": data["count"],
                "average_qs": np.mean(data["qs_values"]),
                "average_accuracy": np.mean(data["acc_values"]),
                "median_length": np.median(data["lengths"]),
                "stdev_length": np.std(data["lengths"]),
                "raw_lengths": data["lengths"],
                "region": region_str
            }
        
        return final_stats

    # --- CLI & ARGS ---
    def handle_cli_args(self, args):
        self.log("Initializing from CLI arguments...")
        if args.primers:
            p_path = os.path.expanduser(args.primers)
            self.primer_file_path = os.path.abspath(p_path)
            self.primer_dict = self.resources.load_primers(self.primer_file_path)
            if self.primer_dict:
                self.l_primers.setText("Loaded (CLI)")
                self.combo_mode.setCurrentText("Amplicon"); self.mode = "Amplicon"
        if args.secret:
            self.secret = args.secret
            
        if args.server:
            self.server_address = args.server
            self.statusBar().showMessage(f"Mode: Remote Server ({self.server_address})")
            
        if args.ref:
            ref_path = os.path.expanduser(args.ref)
            if os.path.exists(ref_path):
                self.reference_path = os.path.abspath(ref_path)
                self.log(f"Reference genome set: {self.reference_path}")
            else:
                self.log(f"Warning: Reference file not found: {ref_path}")
        if args.genes:
            gene_path = os.path.expanduser(args.genes)
            gene_path = os.path.abspath(gene_path)
            if gene_path.lower().endswith(('.gtf', '.gff', '.gz')):
                success, msg = self.resources.load_gtf(gene_path)
                if success:
                    self.gene_models = self.resources.gene_trees # Use Trees!
                    self.gene_file_path = gene_path
                    self.l_genes.setText("GTF Loaded (CLI)")
                    self.log(f"DEBUG: GTF loaded successfully from {gene_path}")
                    
                    # Only set mode if not already set by primers
                    if not self.primer_dict:
                        self.combo_mode.setCurrentText("Amplicon"); self.mode = "Amplicon"
                else:
                    self.log(f"Error loading GTF: {msg}")
            else:
                self.gene_list = self.resources.load_simple_bed(gene_path)
                if self.gene_list:
                    self.gene_models = self.gene_list # Ensure gene_models is set!
                    self.l_genes.setText("BED Loaded (CLI)")
                    
                    if not self.primer_dict:
                        self.combo_mode.setCurrentText("RNA-Seq"); self.mode = "RNA-Seq"
        
        input_path = args.input or args.bam
        if input_path:
            input_path = os.path.expanduser(input_path)
            if os.path.isfile(input_path):
                input_path = os.path.abspath(input_path)
                self.log(f"CLI: File {input_path}")
                self.chk_monitor_dir.setChecked(False)
                self.file_queue = [input_path]; self.current_bam_path = input_path;
                self.update_selected_source_label()
                self.check_ready()
                self.b_snap.setEnabled(True); self.b_variant.setEnabled(True)
                if self.mode == "Amplicon" or (self.mode == "RNA-Seq" and self.gene_list): 
                    self.process_next()
                elif self.mode == "Amplicon" and not self.primer_dict:
                     self.discovery_mode = True
                     self.log("CLI: Starting in Discovery Mode.")
                     self.process_next()
                else:
                     self.log(f"CLI Warning: Missing references.")
            elif os.path.isdir(input_path):
                self.log(f"CLI: Directory {input_path}")
                self.chk_monitor_dir.setChecked(True)
                self.monitor_dir = input_path
                self.update_selected_source_label()
                self.check_ready()
                if self.b_toggle.isEnabled(): self.b_toggle.setChecked(True); self.toggle_monitor(True)

    def select_single_file(self):
        if self.is_processing or self.is_monitoring:
            QMessageBox.warning(self, "Busy", "Cannot run single file while monitoring or processing queue.")
            return
        self.is_processing = False 
        self.file_queue = []
        self.discovery_mode = (self.mode == "Amplicon" and self.primer_dict is None)
        self.chk_monitor_dir.setChecked(False)
        
        if self.server_address:
            # Remote File Selection
            dlg = RemoteFileDialog(self.server_address, secret=self.secret, parent=self)
            if dlg.exec():
                p = dlg.selected_file
            else:
                p = None
        else:
            # Local File Selection
            p, _ = QFileDialog.getOpenFileName(self, "Open Read File", "", "Read Files (*.bam *.fastq *.fastq.gz *.fq *.fq.gz);;BAM Files (*.bam);;FASTQ Files (*.fastq *.fastq.gz *.fq *.fq.gz)")
            
        if p:
            self.current_bam_path = p
            self.monitor_dir = None
            self.setWindowTitle(f"NanoStream Monitor - {os.path.basename(p)}")
            self.update_selected_source_label()
            self.check_ready()
            self.start_selected_file_analysis(include_primer_analysis=False, reset_metadata=True)

    def select_monitor_dir(self):
        self.chk_monitor_dir.setChecked(True)
        if self.server_address:
             dlg = RemoteFileDialog(self.server_address, secret=self.secret, parent=self, mode="Directory")
             if dlg.exec():
                 directory = dlg.selected_file
             else:
                 return
        else:
             directory = QFileDialog.getExistingDirectory(self, "Monitor Directory")
        if directory:
            self.monitor_dir = directory
            self.current_bam_path = None
            self.update_selected_source_label()
            self.check_ready()

    def select_input_source(self):
        if self.chk_monitor_dir.isChecked():
            self.select_monitor_dir()
        else:
            self.select_single_file()

    def on_monitor_dir_toggled(self, checked):
        if checked:
            self.current_bam_path = None
            self.file_queue = []
        else:
            self.monitor_dir = None
        self.b_toggle.setText("Start")
        self.update_selected_source_label()
        self.check_ready()

    def on_primer_analysis_toggled(self, checked):
        if self.chk_monitor_dir.isChecked():
            return
        if not self.current_bam_path or self.is_processing:
            return
        if checked and self.mode == "RNA-Seq" and self.gene_list is None and self.gene_models is None:
            QMessageBox.warning(self, "Missing", "Load BED or GTF before RNA analysis.")
            self.chk_primer_analysis.setChecked(False)
            return
        if checked:
            self.log(f"Running primer analysis for {os.path.basename(self.current_bam_path)}...")
            self.start_selected_file_analysis(include_primer_analysis=True, reset_metadata=False)
        else:
            self.global_stats = {}
            self.current_file_stats = {}
            self.table.setRowCount(0)
            self.refresh_table()
            self.log("Primer analysis disabled. Showing statistics only.")

    def update_selected_source_label(self):
        if self.chk_monitor_dir.isChecked():
            source = self.monitor_dir
        else:
            source = self.current_bam_path
        self.lbl_dir.setText(os.path.basename(source) if source else "None")

    def start_selected_file_analysis(self, include_primer_analysis=False, reset_metadata=True):
        if not self.current_bam_path:
            return
        if self.is_processing:
            return

        if reset_metadata:
            self.clear_session()
            self.reset_current_file_buffers()
            self.file_queue = [self.current_bam_path]
        else:
            self.file_queue = [self.current_bam_path]
            self.current_file_stats = {}

        self.is_monitoring = False
        self.is_processing = True
        self.loaded_stats_only = not include_primer_analysis
        self.l_status.setText(
            f"Analyzing {os.path.basename(self.current_bam_path)}..."
            if include_primer_analysis else
            f"Loading statistics: {os.path.basename(self.current_bam_path)}..."
        )
        self.progress.setValue(0)
        self.b_toggle.setChecked(False)

        genes = self.gene_models if self.gene_models else self.gene_list
        config = {
            "primers": self.primer_dict,
            "genes": genes,
            "qc_only": not include_primer_analysis,
        }
        filters = {
            "min_qs": self.s_qs.value(),
            "min_len": 300,
            "max_len": self.s_max_len.value(),
            "duplex_only": self.chk_duplex.isChecked(),
        }

        if self.server_address:
            self.worker_thread = ns_workers.RemoteAnalysisWorker(
                self.server_address,
                self.current_bam_path,
                self.mode,
                config,
                filters,
                threads=self.threads,
                secret=self.secret,
            )
        else:
            use_rust = (
                include_primer_analysis and
                self.chk_use_rust.isChecked() and
                self.primer_file_path and
                self.mode == "Amplicon"
            )
            if use_rust:
                self.worker_thread = ns_workers.NanostreamWorker(
                    self.current_bam_path, self.primer_file_path, threads=self.threads,
                    primer_tolerance=self.s_tolerance.value()
                )
            else:
                self.worker_thread = ns_workers.AnalysisWorker(
                    self.current_bam_path,
                    self.mode,
                    config,
                    filters,
                    threads=self.threads,
                    primer_tolerance=self.s_tolerance.value(),
                    collect_metadata=reset_metadata,
                )

        self.worker_thread.progress.connect(self.update_prog)
        self.worker_thread.results.connect(self.on_results)
        self.worker_thread.partial_results.connect(self.update_live_data)
        self.worker_thread.finished.connect(self.on_finished)
        self.worker_thread.error.connect(self.on_error)
        self.worker_thread.start()
    def change_mode(self, m):
        self.mode = m
        self.stack.setCurrentIndex(0 if m == "Amplicon" else 1 if m == "RNA-Seq" else 2)
        is_amplicon = (m == "Amplicon")
        self.chk_primer_analysis.setEnabled(is_amplicon)
        if not is_amplicon and self.chk_primer_analysis.isChecked():
            self.chk_primer_analysis.setChecked(False)
        self.check_ready()
    def check_ready_with_warning(self):
        if self.mode == "Amplicon" and self.primer_dict is None:
             # QMessageBox.information(self, "Discovery", "Running in Primer Discovery Mode.")
             return True
        if self.mode == "RNA-Seq" and self.gene_list is None and self.gene_models is None:
             QMessageBox.warning(self, "Missing", "Load BED or GTF."); return False
        return True
    def check_ready(self):
        ready = False
        has_source = bool(self.monitor_dir) if self.chk_monitor_dir.isChecked() else bool(self.current_bam_path)
        if has_source:
            if self.mode in ["Amplicon", "DNA"]:
                ready = True
            elif self.mode == "RNA-Seq":
                ready = (self.gene_list is not None or self.gene_models is not None)
        self.b_toggle.setEnabled(ready)
        self.b_toggle.setText("Start")
        self.b_toggle.setVisible(self.chk_monitor_dir.isChecked())
    # def toggle_monitor(self, active): Removed duplicate
    #    if active:
    #        if not self.check_ready_with_warning(): self.b_toggle.setChecked(False); return
    #        self.is_monitoring = True; self.b_toggle.setText("Stop"); self.start_watcher()
    #    else:
    #        self.is_monitoring = False; self.b_toggle.setText("Start"); self.stop_watcher()
    def start_watcher(self):
        if self.server_address:
             self.watcher = ns_workers.RemoteDirectoryWatcher(self.server_address, self.monitor_dir, secret=self.secret)
        else:
             self.watcher = ns_workers.DirectoryWatcher(self.monitor_dir)
        self.watcher.new_files_found.connect(self.add_files)
        self.watcher.start()
    def stop_watcher(self):
        if hasattr(self, 'watcher'): self.watcher.stop()
        if self.worker_thread: self.worker_thread.terminate(); self.is_processing=False; self.l_status.setText("Stopped")
    
    def stop_processing(self):
        """Stop current processing."""
        if self.worker_thread:
            self.worker_thread.terminate()
            self.is_processing = False
            self.l_status.setText("Stopped by User")
            self.log("Processing stopped by user.")
            
    def add_files(self, files):
        new_bcs = False
        for f in files:
            if f not in self.file_queue: self.file_queue.append(f)
            
            # Track full path for "Scan All" feature
            fname = os.path.basename(f)
            self.all_seen_files[fname] = f
            
            # Detect Barcode
            bc = self.extract_barcode(f)
            if bc and bc != "unknown" and bc not in self.detected_barcodes:
                self.detected_barcodes.add(bc)
                self.combo_barcode.addItem(bc)
                new_bcs = True
        
        # Enable Scan All if we have files
        if self.all_seen_files:
            self.b_scan_all.setEnabled(True)
            
        if self.is_monitoring: self.process_next()
    
    def on_file_selected(self, filename):
        """Handle file selection from dropdown."""
        self.selected_file = filename
        self.refresh_table()
        self.update_plots()  # Force plot update
    
    def process_next(self):
        if self.is_processing or not self.file_queue: return
        f = self.file_queue[0]
        
        # Detect barcode for current file
        self.current_file_barcode = self.extract_barcode(f)
        
        # Don't reset global arrays - accumulate across files for "All Files" view
        # Only reset current file arrays
        self.reset_current_file_buffers()
        
        self.is_processing = True; self.l_status.setText(f"Processing: {os.path.basename(f)} ({self.current_file_barcode})")
        genes = self.gene_models if self.gene_models else self.gene_list
        config = {
            "primers": self.primer_dict, 
            "genes": genes, 
            "qc_only": self.chk_qc_only.isChecked()
        }
        filters = {"min_qs": 0, "min_len": 300, "max_len": self.s_max_len.value(), "duplex_only": self.chk_duplex.isChecked()} 
        if self.server_address:
            from ns_workers import RemoteAnalysisWorker
            print(f"Starting Remote Analysis on {self.server_address}...")
            self.worker_thread = RemoteAnalysisWorker(
                self.server_address,
                f,
                self.mode,
                config,
                filters,
                threads=self.threads,
                secret=self.secret
            )
        else:
            # Check if we should use nanostream (Rust backend)
            use_rust = self.chk_use_rust.isChecked() and self.primer_file_path and self.mode == "Amplicon"
            
            if use_rust:
                self.log(f"Using nanostream (Rust) for {os.path.basename(f)}...")
                self.worker_thread = ns_workers.NanostreamWorker(
                    f, self.primer_file_path, threads=self.threads,
                    primer_tolerance=self.s_tolerance.value()
                )
            else:
                self.worker_thread = ns_workers.AnalysisWorker(f, self.mode, config, filters, threads=self.threads)
        self.worker_thread.progress.connect(self.update_prog)
        self.worker_thread.results.connect(self.on_results)
        self.worker_thread.partial_results.connect(self.update_live_data)
        self.worker_thread.finished.connect(self.on_finished)
        self.worker_thread.error.connect(self.on_error)
        self.worker_thread.start()

    def update_max_reads(self, val):
        self.max_reads = val

    def toggle_monitor(self, active=None):
        should_start = active if active is not None else not self.is_monitoring
        
        if should_start:
            if self.is_monitoring: return
            # START
            if not self.chk_monitor_dir.isChecked():
                # Single File Mode
                if not self.current_bam_path:
                    QMessageBox.warning(self, "Error", "Please select a valid file.")
                    return
                if not self.server_address and not os.path.exists(self.current_bam_path):
                    QMessageBox.warning(self, "Error", "Please select a valid file.")
                    return
                
                self.is_monitoring = True
                self.b_toggle.setText("Stop")
                self.l_status.setText(f"Analyzing {os.path.basename(self.current_bam_path)}...")
                self.progress.setValue(0)
                self.clear_session() # Clear previous data
                
                # Start Worker
                filters = {
                    "min_qs": self.s_qs.value(),
                    "min_len": 300, # Enforce 300bp min length for analysis
                    "max_len": self.s_max_len.value(),
                    "duplex_only": self.chk_duplex.isChecked()
                }
                
                config = {
                    "primers": self.primer_dict,
                    "genes": self.gene_list,
                    "qc_only": self.chk_qc_only.isChecked()
                }
                
                
                if self.server_address:
                    self.worker_thread = ns_workers.RemoteAnalysisWorker(
                        self.server_address,
                        self.current_bam_path,
                        self.mode,
                        config,
                        filters,
                        threads=self.threads,
                        secret=self.secret
                    )
                else:
                    self.worker_thread = ns_workers.AnalysisWorker(
                        self.current_bam_path, self.mode, config, filters, threads=self.threads
                    )
                
                self.worker_thread.progress.connect(self.update_prog)
                self.worker_thread.partial_results.connect(self.update_live_data)
                self.worker_thread.results.connect(self.on_results)
                self.worker_thread.error.connect(self.on_error)
                self.worker_thread.finished.connect(self.on_finished)
                self.worker_thread.start()
                
            else:
                # Directory Monitor Mode
                if not self.server_address:
                    if not self.monitor_dir or not os.path.exists(self.monitor_dir):
                        QMessageBox.warning(self, "Error", "Please select a valid directory.")
                        return
                elif not self.monitor_dir:
                     # Check if monitor_dir is set (even if remote)
                     QMessageBox.warning(self, "Error", "Please select a directory to monitor.")
                     return
                
                self.is_monitoring = True
                self.b_toggle.setText("Stop")
                self.l_status.setText(f"Monitoring {os.path.basename(self.monitor_dir)}...")
                
                # Start Watcher
                if self.server_address:
                     self.watcher_thread = ns_workers.RemoteDirectoryWatcher(self.server_address, self.monitor_dir, secret=self.secret)
                else:
                     self.watcher_thread = ns_workers.DirectoryWatcher(self.monitor_dir)
                     
                self.watcher_thread.new_files_found.connect(self.add_files)
                self.watcher_thread.start()
                
        else:
            # STOP
            self.is_monitoring = False
            self.b_toggle.setText("Start")
            self.l_status.setText("Stopped")
            
            if self.watcher_thread:
                self.watcher_thread.stop()
                self.watcher_thread = None
            
            if self.worker_thread:
                self.worker_thread.stop()
                self.worker_thread = None

    def on_finished(self):
        if self.file_queue:
            f = self.file_queue.pop(0)
            self.log(f"Finished {os.path.basename(f)}")
            self.current_bam_path = f 
        else:
            self.log("Finished processing (Queue empty).")

        self.session_reads_processed += self.current_file_reads_processed
        self.current_file_reads_processed = 0
        self.l_status.setText("None")
        if self.current_bam_path:
             self.b_snap.setEnabled(True); self.b_variant.setEnabled(True)
        self.is_processing = False  # Reset processing flag
        if self.is_monitoring: 
            self.process_next()
        elif self.file_queue:
            self.process_next()

    def load_primers(self):
        p, _ = QFileDialog.getOpenFileName(self, "Primers", "", "TSV (*.txt *.tsv)")
        if p: 
            self.primer_file_path = p
            self.primer_dict = self.resources.load_primers(p)
            self.l_primers.setText("Loaded")
            self.check_ready()
    def load_bed(self):
        p, _ = QFileDialog.getOpenFileName(self, "BED", "", "BED (*.bed)")
        if p: self.gene_list = self.resources.load_simple_bed(p); self.gene_models=self.gene_list; self.l_genes.setText("Loaded"); self.check_ready()
    def load_gtf(self):
        p, _ = QFileDialog.getOpenFileName(self, "GTF", "", "GTF/GFF (*.gtf *.gff *.gtf.gz *.gff.gz)")
        if p:
            success, msg = self.resources.load_gtf(p)
            if success:
                self.gene_models = self.resources.gene_trees # Use Trees!
                self.gene_file_path = p
                if self.mode == "Amplicon": QMessageBox.information(self, "Loaded", f"Loaded GTF: {os.path.basename(p)}")
                else: self.l_genes.setText(f"GTF Loaded")
                self.check_ready()
            else: QMessageBox.warning(self, "Error", msg)
    def load_genes(self):
        p, _ = QFileDialog.getOpenFileName(self, "Load Genes (BED/GTF)", "", "Gene Files (*.bed *.gtf *.gff *.gtf.gz *.gff.gz);;BED (*.bed);;GTF/GFF (*.gtf *.gff *.gtf.gz *.gff.gz)")
        if p:
            if p.lower().endswith(('.gtf', '.gff', '.gz')):
                success, msg = self.resources.load_gtf(p)
                if success:
                    self.gene_models = self.resources.gene_trees # Use Trees!
                    self.gene_file_path = p
                    self.l_genes.setText(f"GTF Loaded")
                    self.check_ready()
                else: QMessageBox.warning(self, "Error", msg)
            else:
                self.gene_list = self.resources.load_simple_bed(p)
                if self.gene_list:
                    self.gene_models = self.gene_list
                    self.l_genes.setText("BED Loaded")
                    self.check_ready()
    def update_prog(self, count):
        self.current_file_reads_processed = count
        total = self.session_reads_processed + count
        val = total % 1000000
        self.progress.setValue(val)
        self.progress.setFormat(f"Total: {total:,}")
    def on_error(self, err_msg):
        self.log(f"Error: {err_msg}")
        self.is_processing = False
        self.l_status.setText("Error")
        if self.file_queue: self.file_queue.pop(0) # Remove failed file
        if self.is_monitoring: self.process_next()

    def on_results(self, res):
        self.update_live_data(res)
        self.update_plots() 
        
        # Store per-file stats and plot data
        if self.current_bam_path:
            filename = os.path.basename(self.current_bam_path)
            if self.mode == "Amplicon":
                # Calculate stats for this file including QS
                has_assigned_amplicons = (
                    hasattr(self, 'current_file_amplicons') and
                    any(a and a != "Unknown" for a in self.current_file_amplicons)
                )
                if has_assigned_amplicons:
                    file_stats = self.calculate_amplicon_stats_from_metadata(
                        amplicons=self.current_file_amplicons,
                        qs=self.current_file_qs,
                        acc=self.current_file_acc,
                        lengths=self.current_file_len
                    )
                    self.file_stats[filename] = file_stats
                elif "amplicons" in res:
                    # Fallback to result from worker (no QS)
                    self.file_stats[filename] = res.get("amplicons", {})
                
                # Add to dropdown if not already there
                if self.combo_file_selector.findText(filename) == -1:
                    self.combo_file_selector.addItem(filename)
            
            # Store plot data for this file
            self.file_plot_data[filename] = {
                "qs": self.current_file_qs.copy(),
                "acc": self.current_file_acc.copy(),
                "len": self.current_file_len.copy(),
                "dx": self.current_file_dx.copy(),
                "ids": list(self.current_file_ids),
                "amplicons": list(self.current_file_amplicons),
            }
        
        if "summary" in res:
            summary = res["summary"]
            if "sv_links" in summary:
                new_links = summary["sv_links"]
                if new_links: self.sv_links.extend(new_links); self.b_matrix.setEnabled(True); self.log(f"Added {len(new_links)} SV links.")
            if "chrom_lengths" in summary: self.chrom_lengths.update(summary["chrom_lengths"])
        if self.mode == "Amplicon":
            if "internal_adapter_count" in res.get("summary", {}): self.log(f"Concatemers: {res['summary']['internal_adapter_count']:,}")
            
            bc = getattr(self, "current_file_barcode", "unknown")
            if bc != "unknown" and bc not in self.barcode_stats:
                self.barcode_stats[bc] = {}
                
            for name, data in res.get("amplicons", {}).items():
                # 1. Update Global Stats
                if name not in self.global_stats: 
                    self.global_stats[name] = {
                        "count":0, "average_accuracy":0, "average_qs":0, 
                        "median_len":0, "stdev_len":0, "raw_lengths":[], 
                        "region": data.get("region")
                    }
                
                g = self.global_stats[name]
                g["count"] += data["count"]
                g["average_accuracy"] = data["average_accuracy"] 
                g["average_qs"] = data.get("average_qs", 0)
                if "raw_lengths" in data: g["raw_lengths"].extend(data["raw_lengths"])
                
                # 2. Update Barcode Stats
                if bc != "unknown":
                    if name not in self.barcode_stats[bc]:
                        self.barcode_stats[bc][name] = {
                            "count":0, "average_accuracy":0, "average_qs":0, 
                            "median_len":0, "stdev_len":0, "raw_lengths":[], 
                            "region": data.get("region")
                        }
                    b = self.barcode_stats[bc][name]
                    b["count"] += data["count"]
                    b["average_accuracy"] = data["average_accuracy"]
                    b["average_qs"] = data.get("average_qs", 0)
                    if "raw_lengths" in data: b["raw_lengths"].extend(data["raw_lengths"])
            
            # Clear live stats for this file as they are now in global
            self.current_file_stats = {}
            
            # --- Auto-Trigger Variant Calling ---
            # Gate: Only run if Auto-Scan is enabled
            if self.chk_auto_variant.isChecked():
                # Identify amplicons with regions that haven't been scanned
                amplicons_data = res.get("amplicons", {})
                if amplicons_data:
                    # Determine filename safely
                    filename = "unknown"
                    if self.current_bam_path:
                        filename = os.path.basename(self.current_bam_path)
                    elif "filename" in res:
                        filename = res["filename"]
                
                
                # print(f"DEBUG Auto-Variant: mode={self.mode}, amplicons_count={len(amplicons_data)} in {filename}")
                
                new_candidates = {}
                for name, data in amplicons_data.items():
                    region = data.get("region") if isinstance(data, dict) else None
                    # Fallback for Nanostream which provides raw coords
                    if not region and isinstance(data, dict):
                        c = data.get("chrom")
                        s = data.get("start")
                        e = data.get("end")
                        if c and s is not None and e is not None:
                            region = f"{c}:{s}-{e}"
                    
                    # Faceted Key: (filename, amplicon_name) ensures we call variants for each unique sample/file
                    scan_key = (filename, name)
                    if region and scan_key not in self.scanned_amplicons and name not in self.variant_queue:
                        new_candidates[name] = region
                
                if new_candidates:
                    # print(f"DEBUG: Found {len(new_candidates)} new variant candidates")
                    self.variant_queue.update(new_candidates)
                    self.process_variant_queue()
            else:
                pass
                # print(f"DEBUG: No new variant candidates. Amplicons={len(amplicons_data)}, Scanned={len(self.scanned_amplicons)}")

            
        elif self.mode == "RNA-Seq":
            for k, v in res.get("genes", {}).items(): self.global_stats[k] = self.global_stats.get(k, 0) + v
            
        if "pore_stats" in res:
            stats = res["pore_stats"]
            self.log(f"\n--- Pore Performance (Gap Analysis) ---")
            self.log(f"Global Mean Gap: {stats['global_mean_gap']:.4f} s")
            self.log(f"Global Median Gap: {stats['global_median_gap']:.4f} s")
            
            pct = stats.get('long_gap_pct', 0)
            if pct > 5.0: # Alert threshold
                self.log(f"⚠️ ALERT: {pct:.2f}% of gaps are > 60s!")
            else:
                self.log(f"Long Gaps (>60s): {pct:.2f}%")
                
            self.log(f"Channel Stats (Top 10 by Long Gaps):")
            # Sort by long_gaps desc
            ch_stats = stats.get('channel_stats', [])
            ch_stats.sort(key=lambda x: x.get('long_gaps', 0), reverse=True)
            
            for s in ch_stats[:10]:
                self.log(f"  Ch {s['ch']} Mx {s['mx']}: Mean={s['mean']:.2f}s, Med={s['median']:.2f}s, >60s={s.get('long_gaps',0)}")
            if len(ch_stats) > 10: self.log(f"  ... and {len(ch_stats)-10} more channels.")
            
            if 'potential_duplex' in stats:
                self.log(f"\nPotential Duplex Pairs Found: {stats['potential_duplex']}")
                if stats['potential_duplex'] > 0:
                    self.log(f"  -> Saved list to 'duplex_candidates.txt'")
            
        self.refresh_table()
    def refresh_table(self, custom_stats=None):
        self.table.setRowCount(0)
        if self.mode == "Amplicon":
            self.update_table_headers()
            
            # Determine which stats to display
            if custom_stats:
                display_stats = custom_stats
            elif self.selected_file == "All Files":
                # Check Barcode Selection
                if self.current_barcode != "All Barcodes":
                    display_stats = self.barcode_stats.get(self.current_barcode, {})
                else:
                    # Performance: Use pre-aggregated global stats
                    display_stats = self.global_stats
            else:
                # Show stats for selected file only
                display_stats = self.file_stats.get(self.selected_file, {})
            
            self.table.setRowCount(len(display_stats))
            sorted_items = sorted(display_stats.items(), key=lambda x: x[1]['count'], reverse=True)
            for r, (name, d) in enumerate(sorted_items):
                item = QTableWidgetItem(str(name))
                if d.get("region"):
                    item.setToolTip(f"Genomic Position: {d['region']}")
                self.table.setItem(r, 0, item)
                
                self.table.setItem(r, 1, QTableWidgetItem(f"{d['count']:,}"))
                med = d.get("median_len", d.get("median_length", 0))
                sd = d.get("stdev_len", d.get("stdev_length", 0))
                raw_lengths = d.get("raw_lengths", [])
                if raw_lengths is not None and len(raw_lengths) > 0:
                    med = np.median(raw_lengths)
                    sd = np.std(raw_lengths)
                
                self.table.setItem(r, 2, QTableWidgetItem(f"{int(med)}"))
                self.table.setItem(r, 3, QTableWidgetItem(f"{sd:.1f}"))
                
                # Add QS column
                avg_qs = d.get("average_qs", d.get("avg_qs", 0))
                self.table.setItem(r, 4, QTableWidgetItem(f"{avg_qs:.1f}"))
                
                # Add Variants column (Automated)
                var_info = ""
                full_tooltip = ""
                if name in self.amplicon_variants:
                    vars_found = self.amplicon_variants[name]
                    if vars_found:
                        # Summarize top 2 variants: "A123T (0.45), ..."
                        # Summarize top 2 variants (Prioritize Novel)
                        sorted_vars = sorted(vars_found, key=lambda x: x['af'], reverse=True)
                        
                        # Separate novel vs common
                        novel_vars = []
                        common_vars = []
                        clinical_vars = []
                        
                        for v in sorted_vars:
                            # Check clinical first
                            k_key = (v['chrom'], v['pos'], v['ref'], v['alt'])
                            if k_key in self.known_mutations:
                                clinical_vars.append(v)
                            elif (v['chrom'], v['pos']) in self.common_snps:
                                common_vars.append(v)
                            else:
                                novel_vars.append(v)
                        
                        # Prioritize showing clinical -> novel -> common
                        display_list = clinical_vars + novel_vars
                        if not display_list:
                            display_list = common_vars
                        
                        v_strs = []
                        tooltip_lines = []
                        
                        if clinical_vars:
                            tooltip_lines.append(f"<b><font color='red'>CLINICAL VARIANTS ({len(clinical_vars)}):</font></b>")
                            for v in clinical_vars:
                                k_name = self.known_mutations.get((v['chrom'], v['pos'], v['ref'], v['alt']), "Unknown")
                                tooltip_lines.append(f"<font color='red'>{v['ref']}{v['pos']}{v['alt']} ({k_name}) AF={v['af']:.3f}</font>")
                        
                        if novel_vars:
                            tooltip_lines.append(f"<b>Novel Variants ({len(novel_vars)}):</b>")
                            for v in novel_vars[:5]:
                                tooltip_lines.append(f"{v['ref']}{v['pos']}{v['alt']} AF={v['af']:.3f}")
                        
                        if common_vars:
                            tooltip_lines.append(f"<b>dbSNP Variants ({len(common_vars)}):</b>")
                            for v in common_vars[:5]:
                                rsid = self.common_snps.get((v['chrom'], v['pos']), "dbSNP")
                                tooltip_lines.append(f"{v['ref']}{v['pos']}{v['alt']} ({rsid}) AF={v['af']:.3f}")

                        # Text Summary (First 2 interesting ones)
                        for v in display_list[:2]:
                             v_strs.append(f"{v['ref']}{v['pos']}{v['alt']} ({v['af']:.3f})")
                        
                        var_info = ", ".join(v_strs)
                        remaining = len(vars_found) - 2
                        if remaining > 0:
                            var_info += f" (+{remaining})"
                            
                        # If we showed novel only, but there are dbSNP ones, maybe mention it
                        if novel_vars and common_vars:
                            var_info += f" [{len(common_vars)} dbSNP]"
                            
                        full_tooltip = "<br>".join(tooltip_lines)
                    else:
                        var_info = "None"
                        full_tooltip = "No variants > 2%"
                
                var_item = QTableWidgetItem(var_info)
                if full_tooltip:
                    var_item.setToolTip(full_tooltip)
                self.table.setItem(r, 5, var_item)

        elif self.mode == "DNA":
            self.update_table_headers()
            
            rows_to_show = []
            min_qs = self.s_qs.value()
            min_len = self.s_len.value()
            max_len = self.s_max_len.value()
            duplex_only = self.chk_duplex.isChecked()

            def get_stats_for_arrays(qs_arr, len_arr, dx_arr):
                if len(qs_arr) == 0:
                    return 0, 0, 0.0, 0.0
                mask = (qs_arr >= min_qs) & (len_arr >= min_len)
                if max_len > 0:
                    mask = mask & (len_arr <= max_len)
                if duplex_only:
                    mask = mask & (dx_arr == 1)
                
                f_qs = qs_arr[mask]
                f_len = len_arr[mask]
                count = len(f_len)
                if count > 0:
                    return count, int(np.median(f_len)), float(np.std(f_len)), float(np.mean(f_qs))
                return 0, 0, 0.0, 0.0

            if self.selected_file == "All Files":
                # Whole Run Row
                count, med, sd, avg_q = get_stats_for_arrays(self.read_qs, self.read_len, self.read_dx)
                rows_to_show.append(("Whole Run", count, med, sd, avg_q))
                
                # Per File Rows
                for filename, file_data in sorted(self.file_plot_data.items()):
                    f_qs = file_data.get("qs", np.array([]))
                    f_len = file_data.get("len", np.array([]))
                    f_dx = file_data.get("dx", np.array([]))
                    count, med, sd, avg_q = get_stats_for_arrays(f_qs, f_len, f_dx)
                    rows_to_show.append((filename, count, med, sd, avg_q))
            else:
                # Selected File Row Only
                file_data = self.file_plot_data.get(self.selected_file, {})
                f_qs = file_data.get("qs", np.array([]))
                f_len = file_data.get("len", np.array([]))
                f_dx = file_data.get("dx", np.array([]))
                count, med, sd, avg_q = get_stats_for_arrays(f_qs, f_len, f_dx)
                rows_to_show.append((self.selected_file, count, med, sd, avg_q))

            self.table.setRowCount(len(rows_to_show))
            for r, (label, count, med, sd, avg_q) in enumerate(rows_to_show):
                self.table.setItem(r, 0, QTableWidgetItem(str(label)))
                self.table.setItem(r, 1, QTableWidgetItem(f"{count:,}"))
                self.table.setItem(r, 2, QTableWidgetItem(f"{med}"))
                self.table.setItem(r, 3, QTableWidgetItem(f"{sd:.1f}"))
                self.table.setItem(r, 4, QTableWidgetItem(f"{avg_q:.1f}"))

        else:
            self.update_table_headers()
            sorted_items = sorted(self.global_stats.items(), key=lambda x: x[1], reverse=True)
            self.table.setRowCount(len(sorted_items))
            for r, (k, v) in enumerate(sorted_items):
                self.table.setItem(r, 0, QTableWidgetItem(str(k)))
                self.table.setItem(r, 1, QTableWidgetItem(f"{v:,}"))
    def update_table_headers(self):
        """Reset table headers based on mode."""
        if self.mode == "Amplicon":
            if self.discovery_mode:
                self.table.setColumnCount(4)
                self.table.setHorizontalHeaderLabels(["Amplicon (Candidate)", "Count", "Med. Length", "Avg. QS"])
            else:
                self.table.setColumnCount(6)
                self.table.setHorizontalHeaderLabels(["Amplicon Name", "Count", "Med. Length", "SD Length", "Avg. QS", "Variants (>2%)"])
            
            # Widen first column (Amplicon Name)
            self.table.setColumnWidth(0, 450)
            
            # Allow numerical columns to shrink, name to stretch/stay wide
            header = self.table.horizontalHeader()
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
            for i in range(1, self.table.columnCount()):
                header.setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)
            
        elif self.mode == "RNA-Seq":
            self.table.setColumnCount(2)
            self.table.setHorizontalHeaderLabels(["Gene", "Reads"])
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        elif self.mode == "DNA":
            self.table.setColumnCount(5)
            self.table.setHorizontalHeaderLabels(["Sample / Barcode", "Count", "Med. Length", "SD Length", "Avg. QS"])
            self.table.setColumnWidth(0, 300)
            header = self.table.horizontalHeader()
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
            for i in range(1, 5):
                header.setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)

    def update_selection_count(self):
        """Update the selection counter when table selection changes."""
        selected_rows = self.table.selectionModel().selectedRows()
        total_reads = 0
        amplicon_names = []
        
        for row_index in selected_rows:
            row = row_index.row()
            # Get amplicon name
            name_item = self.table.item(row, 0)
            if name_item:
                amplicon_names.append(name_item.text())
            
            # Get count
            count_item = self.table.item(row, 1)
            if count_item:
                count_str = count_item.text().replace(",", "")
                try:
                    total_reads += int(count_str)
                except ValueError:
                    pass
        
        self.selected_amplicons = amplicon_names
        self.l_selected_count.setText(f"Selected: {len(amplicon_names)} amplicons, {total_reads:,} reads")
        self.b_export.setEnabled(len(amplicon_names) > 0)
        
        # Trigger plot update to show highlighting
        self.update_plots()
        
        # Store selected amplicons for export
        self.selected_amplicons = amplicon_names
    
    def export_selected_reads(self):
        """Export reads for selected amplicons to a new file."""
        if not hasattr(self, 'selected_amplicons') or not self.selected_amplicons:
            QMessageBox.warning(self, "No Selection", "Please select amplicons to export.")
            return
        
        if not self.current_bam_path:
            QMessageBox.warning(self, "No Source File", "No BAM/FASTQ file has been processed yet.")
            return
        
        # Ask for output file
        default_name = f"exported_{len(self.selected_amplicons)}_amplicons.bam"
        output_path, _ = QFileDialog.getSaveFileName(
            self, "Export Reads", default_name,
            "BAM Files (*.bam);;FASTQ Files (*.fastq *.fastq.gz)"
        )
        
        if not output_path:
            return
        
        # Get filter settings
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        max_len = self.s_max_len.value()
        duplex_only = self.chk_duplex.isChecked()
        
        # Show progress dialog
        progress_dialog = QMessageBox(self)
        progress_dialog.setWindowTitle("Exporting Reads")
        progress_dialog.setText(f"Exporting {len(self.selected_amplicons)} amplicons...\nThis may take a while.")
        progress_dialog.setStandardButtons(QMessageBox.StandardButton.NoButton)
        progress_dialog.show()
        
        # Start export in background
        from ns_workers import ExportWorker
        self.export_worker = ExportWorker(
            self.current_bam_path,
            output_path,
            self.selected_amplicons,
            self.primer_dict,
            min_qs,
            min_len,
            duplex_only,
            max_len=max_len
        )
        self.export_worker.finished.connect(lambda count: self.on_export_finished(count, progress_dialog))
        self.export_worker.error.connect(lambda err: self.on_export_error(err, progress_dialog))
        self.export_worker.start()
    
    def on_export_finished(self, count, dialog):
        """Handle export completion."""
        dialog.close()
        QMessageBox.information(
            self, "Export Complete",
            f"Successfully exported {count:,} reads to the output file."
        )
        self.log(f"Exported {count:,} reads for {len(self.selected_amplicons)} amplicons.")
    
    def on_export_error(self, error_msg, dialog):
        """Handle export error."""
        dialog.close()
        QMessageBox.critical(self, "Export Error", f"Failed to export reads:\n{error_msg}")
        self.log(f"Export error: {error_msg}")

    def make_pdf(self): 
        # Ask for filtering options
        dialog = QDialog(self)
        dialog.setWindowTitle("Snap Options")
        layout = QVBoxLayout()
        
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        
        chk_filter = QCheckBox(f"Apply current filters (QS >= {min_qs}, Len >= {min_len})")
        chk_filter.setChecked(False)
        layout.addWidget(chk_filter)
        
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        dialog.setLayout(layout)
        
        if dialog.exec() == QDialog.DialogCode.Accepted:
            if chk_filter.isChecked():
                # Apply filters
                mask = (self.read_qs >= min_qs) & (self.read_len >= min_len)
                
                # Filter arrays
                f_qs = self.read_qs[mask]
                f_acc = self.read_acc[mask]
                f_len = self.read_len[mask]
                
                # Handle read_amplicons (list)
                f_amplicons = np.array(self.read_amplicons)[mask]
                
                # Calculate stats for filtered data
                stats = self.calculate_amplicon_stats_from_metadata(f_amplicons, f_qs, f_acc, f_len)
                total_reads = len(f_qs)
                
                ns_plotting.generate_pdf_report(stats, self.mode, total_reads, self)
            else:
                # Use all data (recalculate to ensure consistency)
                stats = self.calculate_amplicon_stats_from_metadata()
                total_reads = len(self.read_qs)
                ns_plotting.generate_pdf_report(stats, self.mode, total_reads, self)
    def recalculate_table(self):
        """Recalculate table stats based on current filters."""
        if self.mode == "DNA":
            self.log("Recalculating DNA statistics with current filters...")
            self.refresh_table()
            return
            
        if len(self.read_amplicons) == 0: return
        
        self.log("Recalculating table with current filters...")
        
        # 1. Get Mask
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        max_len = self.s_max_len.value()
        duplex_only = self.chk_duplex.isChecked()
        
        qs_data = self.read_qs
        len_data = self.read_len
        dx_data = self.read_dx
        
        # Handle case where arrays might be shorter than list (race condition?)
        n = min(len(qs_data), len(len_data), len(self.read_amplicons))
        
        mask = (qs_data[:n] >= min_qs) & (len_data[:n] >= min_len)
        if max_len > 0:
            mask = mask & (len_data[:n] <= max_len)
        if duplex_only:
            mask = mask & (dx_data[:n] == 1)
            
        # 2. Filter Data
        indices = np.where(mask)[0]
        filtered_amps = [self.read_amplicons[i] for i in indices]
        filtered_lens = len_data[indices]
        filtered_qs = qs_data[indices]
        
        # 3. Aggregate
        from collections import defaultdict
        amp_stats = defaultdict(lambda: {"count": 0, "lengths": [], "qs_sum": 0.0})
        
        for i, amp in enumerate(filtered_amps):
            if not amp or amp == "Unknown": continue # Skip empty or Unknown amplicons
            amp_stats[amp]["count"] += 1
            amp_stats[amp]["lengths"].append(filtered_lens[i])
            amp_stats[amp]["qs_sum"] += filtered_qs[i]
            
        # 4. Format for update_table
        final_stats = {"amplicons": {}}
        for amp, data in amp_stats.items():
            count = data["count"]
            lengths = np.array(data["lengths"])
            mean_qs = data["qs_sum"] / count if count > 0 else 0
            
            final_stats["amplicons"][amp] = {
                "count": count,
                "median_length": np.median(lengths) if count > 0 else 0,
                "stdev_length": np.std(lengths) if count > 0 else 0,
                "average_qs": mean_qs,
                "raw_lengths": lengths # Optional, for PDF
            }
            
        # 5. Update Table
        self.refresh_table(custom_stats=final_stats["amplicons"])
        self.l_filtered_count.setText(f"Total: {len(self.read_ids):,} | Filtered: {len(indices):,}")
        self.log(f"Recalculation complete. Showing {len(indices)} reads.")

    def open_snap_view(self): 
        if not self.current_bam_path: 
            QMessageBox.warning(self, "No BAM", "Please load a BAM file first.")
            return
        
        # Get selected rows
        selected_rows = self.table.selectionModel().selectedRows()
        
        if len(selected_rows) == 1:
            # Single amplicon mode: Create temp BAM
            row = selected_rows[0].row()
            full_name = self.table.item(row, 0).text()
            
            # Use pre-parsed region if possible
            display_stats = {}
            if self.selected_file == "All Files":
                display_stats = self.calculate_amplicon_stats_from_metadata()
            else:
                display_stats = self.file_stats.get(self.selected_file, {})
                
            coords = display_stats.get(full_name, {}).get("region")
            
            if not coords:
                # Fallback to robust parsing
                chrom, start, end = self.parse_genomic_region(full_name)
                if chrom:
                    coords = f"{chrom}:{start}-{end}"
            
            if not coords:
                QMessageBox.warning(self, "Invalid Region", f"Could not parse region from: {full_name}")
                return
            
            # Collect read IDs matching selection and current filters
            min_qs = self.s_qs.value()
            min_len = self.s_len.value()
            
            target_ids = set()
            
            # Determine which arrays to use (All vs Selective)
            if self.selected_file == "All Files":
                ids_data = self.read_ids
                amp_data = self.read_amplicons
                qs_data = self.read_qs
                len_data = self.read_len
            else:
                f_data = self.file_plot_data.get(self.selected_file, {})
                ids_data = f_data.get("ids", [])
                amp_data = f_data.get("amplicons", [])
                qs_data = f_data.get("qs", np.array([]))
                len_data = f_data.get("len", np.array([]))
            
            # Find matching reads
            for i in range(len(ids_data)):
                if i < len(amp_data) and i < len(qs_data) and i < len(len_data):
                    if amp_data[i] == full_name:
                        if qs_data[i] >= min_qs and len_data[i] >= min_len:
                            target_ids.add(ids_data[i])
            
            if not target_ids:
                QMessageBox.warning(self, "No Reads", f"No reads found matching filters for: {full_name}")
                return
                
            self.launch_rsnap(self.current_bam_path, region=coords, target_read_ids=target_ids)
            
        else:
            # Multiple or no selection: Use original BAM
            region = None
            if len(selected_rows) > 1:
                # Use region of first selected row
                row = selected_rows[0].row()
                full_name = self.table.item(row, 0).text()
                
                # Try to get structured region
                display_stats = {}
                if self.selected_file == "All Files":
                    display_stats = self.calculate_amplicon_stats_from_metadata()
                else:
                    display_stats = self.file_stats.get(self.selected_file, {})
                
                region = display_stats.get(full_name, {}).get("region")
                
                if not region:
                    # Fallback to parsing
                    chrom, start, end = self.parse_genomic_region(full_name)
                    if chrom:
                        region = f"{chrom}:{start}-{end}"
                
            self.launch_rsnap(self.current_bam_path, region=region)

    def launch_rsnap(self, bam_path, region=None, target_read_ids=None):
        """Launches rsnap in a separate process."""
        if not bam_path: return
        
        temp_files = []
        bam_to_use = bam_path
        
        # 1. Selection-based filtering (One amplicon selected)
        if target_read_ids:
            self.log(f"Filtering {len(target_read_ids)} reads for specialized SNAP...")
            try:
                # Use a local temp directory to avoid /var/folders quarantine issues on macOS
                local_tmp = os.path.expanduser("~/nanoStream_tmp")
                os.makedirs(local_tmp, exist_ok=True)
                
                tf = tempfile.NamedTemporaryFile(suffix=".bam", delete=False, dir=local_tmp)
                tmp_bam = tf.name
                tf.close()
                temp_files.append(tmp_bam)
                
                with pysam.AlignmentFile(bam_path, "rb") as infile:
                    with pysam.AlignmentFile(tmp_bam, "wb", template=infile) as outfile:
                        # Optimization: Fetch only from region if available
                        if region and ":" in region:
                            try:
                                r_chrom, r_int = region.split(':')
                                r_start, r_end = map(int, r_int.split('-'))
                                it = infile.fetch(r_chrom, r_start, r_end)
                            except: it = infile.fetch()
                        else:
                            it = infile.fetch()
                            
                        count = 0
                        for read in it:
                            if read.query_name in target_read_ids:
                                outfile.write(read)
                                count += 1
                
                if count > 0:
                    pysam.index(tmp_bam)
                    bam_to_use = tmp_bam
                    
                    # CRITICAL: Strip quarantine attributes (macOS) to prevent 'killed' signal
                    if sys.platform == 'darwin':
                        try:
                            subprocess.run(["xattr", "-c", tmp_bam], check=False)
                            subprocess.run(["xattr", "-c", tmp_bam + ".bai"], check=False)
                        except: pass
                        
                    self.log(f"Created temp BAM with {count} reads.")
                else:
                    self.log("Warning: No reads matched target IDs. Using original BAM.")
            except Exception as e:
                self.log(f"Error creating temp BAM: {e}")
        
        # 2. Construct Command
        cmd = ["rsnap", "--viewer", "-b", bam_to_use]
        
        if region and ":" in region:
            # Add padding
            try:
                chrom, interval = region.split(':')
                s_e = interval.split('-')
                start = int(s_e[0].replace(",",""))
                end = int(s_e[1].replace(",",""))
                pad = int((end-start)*0.1)
                padded = f"{chrom}:{max(1, start-pad)}-{end+pad}"
                cmd.extend(["-p", padded])
            except:
                cmd.extend(["-p", region])
        
        if self.gene_file_path:
            cmd.extend(["-g", self.gene_file_path])
            
        if self.reference_path:
            # Pass the custom reference genome to rsnap
            cmd.extend(["-r", self.reference_path])
        
        self.log(f"Launching rsnap: {' '.join(cmd)}")
        
        # 3. Async Spawn
        try:
            proc = subprocess.Popen(cmd)
            
            # Cleanup thread
            if temp_files:
                def cleanup():
                    proc.wait()
                    for f in temp_files:
                        try:
                            os.remove(f); os.remove(f + ".bai")
                        except: pass
                threading.Thread(target=cleanup, daemon=True).start()
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to launch rsnap: {e}")

    def run_variant_calling(self):
        """Runs fast variant calling on selected amplicon."""
        bam_path = self.current_bam_path
        if not bam_path and self.selected_file and os.path.exists(self.selected_file):
            bam_path = self.selected_file
            
        if not bam_path:
            QMessageBox.warning(self, "No BAM", "Please load a BAM file first.")
            return
            
        selected_rows = self.table.selectionModel().selectedRows()
        if not selected_rows:
            QMessageBox.warning(self, "Select Amplicon", "Please select an amplicon row first.")
            return

        row = selected_rows[0].row()
        full_name = self.table.item(row, 0).text()
        chrom, start, end = self.parse_genomic_region(full_name)
        if not chrom:
            QMessageBox.warning(self, "Invalid Region", f"Could not parse region: {full_name}")
            return
            
        # Collect filtered read IDs (reuse logic)
        target_ids = self.get_target_ids_for_amplicon(full_name)
        if not target_ids:
            QMessageBox.warning(self, "No Reads", "No reads matching current filters.")
            return

        self.log(f"Starting Variant Calling for {full_name}...")
        self.b_variant.setEnabled(False)

        region = f"{chrom}:{start}-{end}"
        from ns_workers import RsnapVariantWorker
        self.var_worker = RsnapVariantWorker(bam_path, region, 0.02)
        self.var_worker.finished.connect(self.on_variant_results)
        self.var_worker.start()

    def get_target_ids_for_amplicon(self, full_name):
        """Returns set of read IDs for an amplicon passing current filters (with subsampling)."""
        all_matches = []
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        max_reads = self.s_filter_max_reads.value()
        
        if self.selected_file == "All Files":
            ids_d, amp_d, qs_d, len_d = self.read_ids, self.read_amplicons, self.read_qs, self.read_len
        else:
            f_data = self.file_plot_data.get(self.selected_file, {})
            ids_d, amp_d, qs_d, len_d = f_data.get("ids", []), f_data.get("amplicons", []), f_data.get("qs", []), f_data.get("len", [])
            
        for i in range(len(ids_d)):
            if i < len(amp_d) and i < len(qs_d) and i < len(len_d):
                if amp_d[i] == full_name and qs_d[i] >= min_qs and len_d[i] >= min_len:
                    all_matches.append(ids_d[i])
        
        # Subsampling
        if max_reads > 0 and len(all_matches) > max_reads:
            # For speed, just take the first N. Alternative: random.sample(all_matches, max_reads)
            return set(all_matches[:max_reads])
            
        return set(all_matches)

    def run_batch_variant_calling(self):
        """Automated variant calling for top amplicons."""
        if not self.current_bam_path: return
        
        # Get top 20 amplicons currently in display
        display_stats = {}
        if self.selected_file == "All Files":
            display_stats = self.calculate_amplicon_stats_from_metadata()
        else:
            display_stats = self.file_stats.get(self.selected_file, {})
            
        sorted_amps = sorted(display_stats.items(), key=lambda x: x[1]['count'], reverse=True)[:20]
        
        tasks = []
        for name, _ in sorted_amps:
            chrom, start, end = self.parse_genomic_region(name)
            if chrom:
                t_ids = self.get_target_ids_for_amplicon(name)
                if t_ids:
                    tasks.append((name, chrom, start, end, t_ids))
        
        if not tasks: return
        
        self.log(f"Starting Batch Variant Calling for {len(tasks)} amplicons...")
        self.b_auto_variant.setEnabled(False)
        self.b_auto_variant.setText(f"Wait (0/{len(tasks)})...")
        
        self.batch_var_worker = ns_variant.BatchVariantWorker(self.current_bam_path, tasks, min_af=0.02)
        self.batch_var_worker.progress.connect(lambda cur, tot: self.b_auto_variant.setText(f"Wait ({cur}/{tot})..."))
        self.batch_var_worker.finished.connect(self.on_batch_variant_results)
        self.batch_var_worker.start()

    def on_batch_variant_results(self, results):
        self.b_auto_variant.setEnabled(True)
        self.b_auto_variant.setText("Auto-Variant")
        self.log(f"Batch Variant Calling complete ({len(results)} amplicons).")
        
        # Merge results into cache
        self.amplicon_variants.update(results)
        self.refresh_table()

    def parse_genomic_region(self, full_name):
        """Robustly parses chrom, start, end from complex amplicon names."""
        name = full_name.split('(')[0] if '(' in full_name else full_name
        parts = name.split(':')
        if len(parts) >= 2:
            chrom = parts[-2]
            interval = parts[-1]
            if '-' in interval:
                try:
                    s_str, e_str = interval.split('-')
                    start = int(s_str.replace(",",""))
                    end = int(e_str.replace(",",""))
                    return chrom, start, end
                except: pass
        return None, None, None

    def on_variant_results(self, success, msg, variants):
        self.b_variant.setEnabled(True)
        self.log(msg)
        
        if not success:
            QMessageBox.critical(self, "Variant Error", msg)
            return
            
        if not variants:
            QMessageBox.information(self, "No Variants", "No variants detected above threshold.")
            return
            
        self.show_variant_dialog(variants, msg)

    def show_variant_dialog(self, variants, title_msg="Variant Results"):
        """Displays interactive variant table dialog."""
        dlg = QDialog(self)
        dlg.setWindowTitle("Variant Results")
        dlg.resize(700, 500)
        vbox = QVBoxLayout(dlg)
        
        lbl = QLabel(f"<b>{title_msg}</b><br/>Select a row and click 'Snap' to view in rsnap (+/- 250bp window).")
        vbox.addWidget(lbl)
        
        table = QTableWidget()
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels(["Pos", "Ref", "Alt", "AF", "Depth", "ID"])
        table.setRowCount(len(variants))
        table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        
        # Sort variants by position
        variants = sorted(variants, key=lambda x: x['pos'])
        
        for r, v in enumerate(variants):
            rs_id = self.common_snps.get((v['chrom'], v['pos']), ".")
            k_key = (v['chrom'], v['pos'], v['ref'], v['alt'])
            clin_name = self.known_mutations.get(k_key, None)
            
            table.setItem(r, 0, QTableWidgetItem(f"{v['chrom']}:{v['pos']}"))
            table.setItem(r, 1, QTableWidgetItem(v['ref']))
            table.setItem(r, 2, QTableWidgetItem(v['alt']))
            table.setItem(r, 3, QTableWidgetItem(f"{v['af']:.3f}"))
            table.setItem(r, 4, QTableWidgetItem(str(v['depth'])))
            
            id_text = rs_id
            bg_color = None
            
            if clin_name:
                id_text = f"{clin_name} (Clinical)"
                bg_color = Qt.GlobalColor.red
            elif rs_id != ".":
                bg_color = Qt.GlobalColor.lightGray
            
            id_item = QTableWidgetItem(id_text)
            if bg_color:
                id_item.setBackground(bg_color)
                if bg_color == Qt.GlobalColor.red:
                    id_item.setForeground(Qt.GlobalColor.white)
            table.setItem(r, 5, id_item)
        
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        vbox.addWidget(table)
        
        bbox = QHBoxLayout()
        btn_snap = QPushButton("Snap to Variant")
        btn_snap.clicked.connect(lambda: self.snap_to_selected_variant(table, variants))
        bbox.addWidget(btn_snap)
        
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(dlg.accept)
        bbox.addWidget(btn_close)
        vbox.addLayout(bbox)
        
        dlg.exec()

    def on_table_double_click(self, row, col):
        """Handle double clicks on table."""
        # Col 5 is Variants
        if col == 5:
            full_name = self.table.item(row, 0).text()
            if full_name in self.amplicon_variants:
                vars_found = self.amplicon_variants[full_name]
                if vars_found:
                    self.show_variant_dialog(vars_found, f"Variants for {full_name}")
                else:
                    QMessageBox.information(self, "No Variants", "No variants found for this amplicon.")

    def snap_to_selected_variant(self, table, variants):
        """Launches rsnap centered on selected variant row."""
        rows = table.selectionModel().selectedRows()
        if not rows:
            QMessageBox.warning(table, "Select Row", "Please select a variant row first.")
            return
            
        idx = rows[0].row()
        v = variants[idx]
        
        chrom = v['chrom']
        pos = v['pos']
        window = 250
        region = f"{chrom}:{max(1, pos-window)}-{pos+window}"
        
        self.log(f"Snapping to variant at {region}...")
        self.launch_rsnap(self.current_bam_path, region=region)

    def run_rsnap_variants(self):
        """Runs rsnap variant calling on selected amplicon."""
        if not self.current_bam_path: return
        
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            QMessageBox.warning(self, "Select Amplicon", "Please select an amplicon row first.")
            return
            
        full_name = self.table.item(rows[0].row(), 0).text()
        
        # Parse region from name "PrimerA-PrimerB:chrom:start-end"
        region = ""
        if ":" in full_name and "-" in full_name:
             try:
                parts = full_name.split(':')
                chrom = parts[-2]
                coords = parts[-1].split('(')[0] # Handle "start-end(Gene)" if present
                region = f"{chrom}:{coords}"
             except:
                 region = ""
        
        if not region:
            QMessageBox.warning(self, "Error", f"Could not parse region from: {full_name}")
            return
            
        self.log(f"Running rsnap variant calling on {region}...")
        self.b_rsnap_var.setEnabled(False)
        self.statusBar().showMessage("Running rsnap...")
        
    def run_rsnap_variants(self):
        """Runs rsnap variant calling on selected amplicon."""
        if not self.current_bam_path: return
        
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            QMessageBox.warning(self, "Select Amplicon", "Please select an amplicon row first.")
            return
            
        full_name = self.table.item(rows[0].row(), 0).text()
        
        # Parse region from name "PrimerA-PrimerB:chrom:start-end"
        region = ""
        if ":" in full_name and "-" in full_name:
             try:
                parts = full_name.split(':')
                chrom = parts[-2]
                coords = parts[-1].split('(')[0] # Handle "start-end(Gene)" if present
                region = f"{chrom}:{coords}"
             except:
                 region = ""
        
        if not region:
            QMessageBox.warning(self, "Error", f"Could not parse region from: {full_name}")
            return
            
        self.log(f"Running rsnap variant calling on {region}...")
        self.b_rsnap_var.setEnabled(False)
        self.statusBar().showMessage("Running rsnap...")
        
        from ns_workers import RsnapVariantWorker
        self.rsnap_worker = RsnapVariantWorker(self.current_bam_path, region, 0.015)
        self.rsnap_worker.finished.connect(self.on_rsnap_variant_results)
        self.rsnap_worker.start()
        
    def on_rsnap_variant_results(self, success, msg, variants):
        self.b_rsnap_var.setEnabled(True)
        self.statusBar().showMessage(msg)
        self.log(msg)
        
        if success:
            if not variants:
                QMessageBox.information(self, "No Variants", "rsnap found no variants above threshold.")
            else:
                self.show_variant_dialog(variants, f"rsnap Results: {len(variants)} found")
        else:
             QMessageBox.critical(self, "rsnap Error", msg)

    def scan_all_files(self):
        """Manually triggers variant calling for ALL files and ALL amplicons."""
        count_queued = 0
        min_qs = self.s_qs.value()
        
        print(f"DEBUG: Scan All triggered. Tracking {len(self.all_seen_files)} files. Min QS: {min_qs}")
        self.log(f"Starting Scan All (Filter QS >= {min_qs})...")
        
        self.pending_file_batches = []
        
        # Iterate over all files we have stats for
        for filename, full_path in self.all_seen_files.items():
            # Do we have amplicons for this file?
            if filename in self.file_stats:
                amplicons_data = self.file_stats[filename]
                
                # Identify candidates
                candidates = {}
                for name, data in amplicons_data.items():
                    # QV Filter
                    avg_qs = data.get("average_qs", data.get("avg_qs", 0))
                    if avg_qs < min_qs:
                        continue
                    
                    region = data.get("region") if isinstance(data, dict) else None
                    
                    # Faceted Key Check (ensure we don't re-queue duplicates if already scanned)
                    scan_key = (filename, name)
                    if region and scan_key not in self.scanned_amplicons:
                         candidates[name] = region
                
                if candidates:
                    self.pending_file_batches.append((full_path, candidates))
                    count_queued += len(candidates)
        
        if self.pending_file_batches:
            self.log(f"Queued {len(self.pending_file_batches)} files with {count_queued} amplicons.")
            self.process_variant_queue() # Trigger processing
        else:
            self.log("No new amplicons found to scan (check QS filter).")

    def open_fusion_matrix(self):
        if not self.sv_links: return
        
        # Calculate filtered read IDs if filter is active
        filtered_ids = None
        # Check if we have active filters (QS, Len, Selection)
        # We can reuse the logic from open_snap_view_region or get_filtered_data
        # But get_filtered_data returns arrays, not IDs.
        # However, self.read_ids corresponds to self.read_qs etc.
        
        # Let's get the mask from get_filtered_data logic
        # But get_filtered_data uses self.selected_file.
        # Matrix view usually shows global data? Or current file?
        # sv_links are accumulated globally in self.sv_links?
        # Wait, self.sv_links is reset in clear_session.
        # But where is it populated?
        # It seems it's NOT populated in update_live_data.
        # It must be populated in on_results or on_finished?
        # I need to check on_results/on_finished/on_file_completed.
        # Assuming sv_links contains all links.
        
        # For now, let's just pass the set of IDs that pass the CURRENT filters on the CURRENT file (or all files).
        # If "All Files" selected:
        qs_data = self.read_qs
        len_data = self.read_len
        ids_data = self.read_ids
        dx_data = self.read_dx
        amp_data = self.read_amplicons
        
        if self.selected_file != "All Files":
             # If specific file selected, maybe we should only show SVs from that file?
             # But sv_links doesn't store filename.
             # So we can only filter by Read ID.
             # If we pass a set of allowed Read IDs, FusionMatrixDialog can filter.
             # So we should construct the set of "Visible Read IDs".
             
             # Get data for selected file
             file_data = self.file_plot_data.get(self.selected_file, {})
             # This is complex because file_plot_data stores arrays but maybe not IDs?
             # Let's use the global arrays and filter by what's visible.
             pass

        # Construct mask
        min_qs = self.s_qs.value()
        min_len = self.s_len.value()
        duplex_only = self.chk_duplex.isChecked()
        
        mask = (qs_data >= min_qs) & (len_data >= min_len)
        if duplex_only:
            mask = mask & (dx_data == 1)
            
        # Also filter by selected rows in table (Amplicons)
        selected_rows = self.table.selectionModel().selectedRows()
        selected_amplicons = set()
        for row in selected_rows:
            selected_amplicons.add(self.table.item(row.row(), 0).text())
            
        # If amplicons selected, further refine mask
        if selected_amplicons:
            # This requires matching amplicon names.
            # self.read_amplicons is a list.
            # Convert to numpy for fast comparison? Or list comp.
            # mask is numpy bool array.
            # Let's iterate.
            
            # Optimization: Get indices where mask is True
            candidate_indices = np.where(mask)[0]
            final_ids = set()
            
            for i in candidate_indices:
                if i < len(amp_data):
                    amp = amp_data[i]
                    # Check if amp matches selected (handling gene names etc if needed)
                    # Exact match for now as table uses same names
                    if amp in selected_amplicons:
                        if i < len(ids_data):
                            final_ids.add(ids_data[i])
        else:
            # No amplicon selection -> All reads passing QS/Len
            final_ids = set()
            candidate_indices = np.where(mask)[0]
            for i in candidate_indices:
                if i < len(ids_data):
                    final_ids.add(ids_data[i])
        
        ns_plotting.FusionMatrixDialog(self.sv_links, self.chrom_lengths, self, filtered_read_ids=final_ids).exec()
    def open_snap_view_region(self, region_str, target_read_ids=None):
        if self.current_bam_path:
             self.launch_rsnap(self.current_bam_path, region=region_str, target_read_ids=target_read_ids)
    def save_session_data(self): pass 
    def load_session_data(self): pass 
    def run_duplex_discovery(self):
        if not self.current_bam_path:
            QMessageBox.warning(self, "No File", "Please load a file first.")
            return
            
        self.log("\nStarting Duplex Discovery (On-Demand)...")
        self.b_run_duplex.setEnabled(False)
        
        filters = {
            "min_qs": self.s_qs.value(),
            "min_len": self.s_len.value()
        }
        
        self.duplex_worker = ns_workers.DuplexWorker(self.current_bam_path, filters)
        self.duplex_worker.progress.connect(lambda p: self.l_status.setText(f"Scanning: {p} reads"))
        self.duplex_worker.results.connect(self.on_duplex_results)
        self.duplex_worker.error.connect(self.on_error)
        self.duplex_worker.finished.connect(lambda: self.b_run_duplex.setEnabled(True))
        self.duplex_worker.start()
        
    # --- BATCH VARIANT CALLING ---
    def process_variant_queue(self):
        """Checks queue and starts batch worker if idle."""
        # 1. Manage Meta-Queue (Pending Batches)
        # If main queue is empty, try to pop from pending file batches
        bam_path_override = None
        
        if not self.variant_queue and hasattr(self, 'pending_file_batches') and self.pending_file_batches:
             # Pop next batch
             next_bam, next_batch = self.pending_file_batches.pop(0)
             self.variant_queue = next_batch
             bam_path_override = next_bam
             self.log(f"Processing batch for {os.path.basename(next_bam)} ({len(next_batch)} amplicons remaining: {len(self.pending_file_batches)} files)")
        
        if not self.variant_queue: return
        
        # If worker exists and is running, wait
        if self.batch_variant_worker and self.batch_variant_worker.isRunning():
            return
            
        # Prepare batch
        # Take a snapshot of current queue
        batch = self.variant_queue.copy()
        self.variant_queue = {} # Clear queue (will be processed now)
        
        
        self.log(f"Auto-scanning variants for {len(batch)} amplicons...")
        
        # Determine BAM to use (Current file? Or rely on specific file?)
        # For now, we use the current file being processed or the selected file?
        # Crucial Issue: Batch worker needs a BAM. on_results comes from a specific file.
        # But 'batch' might contain amplicons from mixed files if we queue indiscriminately.
        # SIMPLIFICATION: We only run on the current file being processed in `process_next`.
        # However, `on_results` emits partials.
        
        # Fix: We need the BAM path. `self.current_bam_path` might be useful, 
        # but `process_next` updates it? No, `process_next` populates `f`.
        # Let's use `self.file_queue[0]`? No, relying on `process_next` logic.
        # Best approach: Pass BAM path in `on_results` or similar.
        # For now, let's use `self.current_bam_path` if available.
        # Wait, `self.current_bam_path` isn't reliably set in `process_next` in the snippets I saw.
        
        # Let's grab it from the worker thread if possible, or just use `self.selected_file` if it's a file?
        # Actually, `ns_workers` passes `res`.
        # Let's assume we are processing `self.worker_thread.bam_path` if it exists.
        
        bam_path = None
        if hasattr(self, 'worker_thread'):
            if hasattr(self.worker_thread, 'bam_file'):
                bam_path = self.worker_thread.bam_file
            elif hasattr(self.worker_thread, 'bam_path'):
                bam_path = self.worker_thread.bam_path
        
        if not bam_path and self.file_queue:
             bam_path = self.file_queue[0]
             
        if not bam_path:
            # Try to grab from selected file if it's a real path
            if self.selected_file:
                 bam_path = self.selected_file
        
        # Override if provided by pending batch logic
        if bam_path_override:
            bam_path = bam_path_override
                 
        if not bam_path:
            print(f"DEBUG: process_variant_queue - No BAM found. Worker: {hasattr(self, 'worker_thread')}, Queue: {len(self.file_queue)}")
            return

        # In server mode, the file exists on the server, not locally.
        if not self.server_address and not os.path.exists(bam_path):
            self.log(f"Skipping variant scan: BAM not found locally.")
            return

        from ns_workers import BatchRsnapWorker
        self.batch_variant_worker = BatchRsnapWorker(
            bam_path, 
            batch, 
            0.015, # Default AF
            self.reference_path,
            server_address=self.server_address,
            secret=self.secret
        )
        
        # Update scanned set with (bam, name) keys NOW that we have a bam_path
        # Fix: Use basename to match on_results logic
        bam_filename = os.path.basename(bam_path)
        for name in batch.keys():
            self.scanned_amplicons.add((bam_filename, name))
            
        self.batch_variant_worker.partial_result.connect(self.on_batch_variant_result)
        self.batch_variant_worker.finished.connect(self.on_batch_variant_finished)
        self.batch_variant_worker.start()
        
    def on_batch_variant_result(self, name, variants):
        """Handle per-amplicon result from batch worker, aggregating results."""
        if name not in self.amplicon_variants:
            self.amplicon_variants[name] = []
        
        # Merge new variants with existing ones based on unique mutation key
        # Key: (chrom, pos, ref, alt)
        existing_map = {}
        for idx, v in enumerate(self.amplicon_variants[name]):
            key = (v.get('chrom'), v.get('pos'), v.get('ref'), v.get('alt'))
            existing_map[key] = idx
            
        for v in variants:
            key = (v.get('chrom'), v.get('pos'), v.get('ref'), v.get('alt'))
            if key in existing_map:
                # Variant exists. Update if this finding has higher AF?
                # Or just keep the first one? Let's keep the one with higher AF.
                idx = existing_map[key]
                old_v = self.amplicon_variants[name][idx]
                if v.get('af', 0) > old_v.get('af', 0):
                    self.amplicon_variants[name][idx] = v
            else:
                # New variant
                self.amplicon_variants[name].append(v)
                # Update map in case duplicates in the same batch
                existing_map[key] = len(self.amplicon_variants[name]) - 1

        # Trigger table update (debounced)
        self.refresh_timer.start(500) # Wait 500ms of quiet before refreshing
        
    def on_batch_variant_finished(self, success, msg):
        if success:
            self.log("Batch variant scan complete.")
        else:
            self.log(f"Batch variant scan error: {msg}")
            
        # Check queue again
        self.process_variant_queue()

    def clear_session(self):
        """Resets all session data and UI."""
        self.stop_processing()
        if self.batch_variant_worker:
            self.batch_variant_worker.stop()
            self.batch_variant_worker.wait()
        
        # Reset Arrays
        self.read_qs = np.array([], dtype=np.float32)
        self.read_acc = np.array([], dtype=np.float32)
        self.read_len = np.array([], dtype=np.int32)
        self.read_dx = np.array([], dtype=np.int8)
        self.read_ids = []
        self.read_amplicons = []
        self.read_concatemers = []
        self.reset_current_file_buffers()
        self.session_reads_processed = 0
        self.global_stats = {}
        self.barcode_stats = {}
        self.file_stats = {}
        self.file_plot_data = {}
        self.selected_file = "All Files"
        self.detected_barcodes = set()
        self.combo_barcode.clear()
        self.combo_barcode.addItem("All Barcodes")
        self.combo_file_selector.clear()
        self.combo_file_selector.addItem("All Files")
        self.sv_links = []
        self.current_file_stats = {} # Reset live stats
        
        # Reset Automated Variant State
        self.variant_queue = {}
        self.scanned_amplicons = set()
        self.amplicon_variants = {}
        self.all_seen_files = {} # Reset specific file tracking
        
        # Reset UI
        self.table.setRowCount(0)
        self.l_status.setText("Idle")
        self.l_filtered_count.setText("Total: 0 | Filtered: 0")
        self.progress.setValue(0)
        self.log("Session cleared.")
        
        # Clear Plots (by sending empty data)
        self.update_acc_plot_js(np.array([]))
        self.update_qs_plot_js(np.array([]))
        self.update_hist_plot_js(np.array([]))
        
        self.b_snap.setEnabled(False)
        self.b_variant.setEnabled(False)
        self.b_matrix.setEnabled(False)
        
    def on_duplex_results(self, msg):
        self.log(f"Duplex Results: {msg}")
        self.l_status.setText("Duplex Done")
        QMessageBox.information(self, "Duplex Discovery", msg)

    def log(self, msg): self.log_view.append(msg)

    def update_clear_button_state(self):
        if self.is_processing or self.is_monitoring:
            self.b_clear.setText("🛑 Stop")
            self.b_clear.setStyleSheet("background-color: #d32f2f; color: white; font-weight: bold; min-height: 28px;")
        else:
            self.b_clear.setText("Clear")
            self.b_clear.setStyleSheet("background-color: #E0E0E0; font-weight: bold; min-height: 28px;")

    def on_clear_clicked(self):
        if self.is_processing or self.is_monitoring:
            if self.is_monitoring:
                self.stop_watcher()
                self.is_monitoring = False
            self.stop_processing()
            self.l_status.setText("Stopped")
            self.log("Analysis stopped by user.")
            self.update_clear_button_state()
        else:
            self.clear_session()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?")
    parser.add_argument("--bam")
    parser.add_argument("--primers")
    parser.add_argument("--genes", default="~/data/gencode.v46.annotation.sorted.gtf.gz")
    parser.add_argument("--ref", default="~/data/homo_sapiens.fasta")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--dir", help="Start in directory mode")
    parser.add_argument("--server", help="Address of NanoStream server (e.g. tcp://127.0.0.1:5555)") # New arg
    parser.add_argument("--secret", help="Authentication Secret")
    args = parser.parse_args()
    app = QApplication(sys.argv)
    w = MainWindow(args)
    w.show()
    sys.exit(app.exec())
