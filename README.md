# nanoStream

`nanoStream` is a high-performance, real-time monitoring and analysis suite for Oxford Nanopore Technologies (ONT) sequencing data. It provides a comprehensive graphical interface (PyQt6) to visualize sequencing metrics, detect structural variants, and monitor targeted amplicons as data is generated.

![nanoStream Overview](snapshot.png)

## Key Features

*   **Real-time Monitoring**: Stream data directly from BAM or FASTQ (including `.gz`) files as they are written by the sequencer.
*   **Interactive DNA Metrics**: Live-updating, high-fidelity plots for Read Accuracy, Quality Scores (Phred), and Read Length distributions.
*   **Structural Variant (SV) Discovery**: Interactive genome-wide contact matrix for real-time identification of translocations and fusions.
*   **Targeted Amplicon Analysis**: Optimized monitoring for specific gene regions with integrated primer and barcode management.
*   **Duplex Read Detection**: Specialized workflows to monitor and identify duplex read pairs and rates.
*   **Remote Scalability**: A robust client-server architecture powered by ZeroMQ, enabling remote heavy-lifting analysis with local visualization.
*   **Genomic Context**: Full support for reference genomes (FASTA), gene models (GTF/Tabix), and BED-based regions of interest.

## Requirements

- **Python**: 3.9+
- **Core Dependencies**:
    - `PyQt6`, `PyQt6-WebEngine`: UI and interactive plot rendering.
    - `pysam`: Genomic data access (BAM, FASTQ, Tabix).
    - `numpy`, `scipy`: Numerical processing and KDE density estimation.
    - `plotly`: Interactive web-based visualizations.
    - `matplotlib`: Static plotting and SV matrix rendering.
    - `pyzmq`: High-performance networking for remote analysis.
    - `intervaltree`: Genomic interval management.
- **Rust**: Required for building the `ns_rust` highly-optimized lane management module.

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/[your-repo]/nanoStream.git
cd nanoStream
```

### 2. Install Python dependencies
```bash
pip install PyQt6 PyQt6-WebEngine pysam numpy scipy plotly matplotlib pyzmq intervaltree
```

### 3. Build the Rust Optimization Module
`nanoStream` uses a Rust extension for performance-critical visualization tasks.
```bash
# Install maturin to build the extension
pip install maturin
maturin develop --release
```

## Usage

### Local Monitoring
To monitor a local BAM file or a directory being written to:
```bash
# Single file monitoring
python nanoMonitor.py --bam /path/to/sequencing_data.bam --primers primers.tsv

# Directory monitoring
python nanoMonitor.py --dir /path/to/output_fastq_dir/ --threads 12
```

### Remote Analysis Mode
`nanoStream` can run the analysis on a powerful remote server while you monitor the progress on your local workstation.

1.  **Start the Server** (on the analysis machine):
    ```bash
    python ns_server.py --rep-port 5555 --pub-port 5556 --secret my_secure_token
    ```

2.  **Start the Monitor** (on your local machine):
    ```bash
    python nanoMonitor.py --server tcp://server_ip:5555 --secret my_secure_token
    ```

### CLI Arguments

| Argument | Description |
| :--- | :--- |
| `input` | Optional positional argument for input file. |
| `--bam` | Path to the input BAM file. |
| `--primers` | Path to a TSV file containing primer names and sequences. |
| `--genes` | Path to a GTF or BED file for gene models. |
| `--threads` | Number of processing threads (default: 8). |
| `--dir` | Path to a directory for continuous monitoring. |
| `--server` | Address of a remote NanoStream server (e.g., `tcp://10.0.0.1:5555`). |
| `--secret` | Authentication secret token for the remote server. |

## Module Architecture

- `nanoMonitor.py`: Main application entry point and PyQt6 GUI.
- `ns_server.py`: Remote analysis server implementation.
- `ns_core.py`: Core streaming engines for BAM and FASTQ data.
- `ns_workers.py`: Multi-threaded worker management for local and remote analysis.
- `ns_structural.py`: Structural variant and contact matrix logic.
- `ns_visualizer.py`: High-performance BAM/FASTQ region viewer (Snap View).
- `ns_rust.rs`: Rust implementation for optimized alignment lane calculation.

---
© 2025 Genomics Suite - Optimized for Nanopore Data Monitoring.
