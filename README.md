# nanoMonitor

`nanoMonitor` is a real-time monitoring and analysis suite for Oxford Nanopore Technologies (ONT) sequencing data. It provides a PyQt6 graphical interface for sequencing metrics, targeted amplicon monitoring, structural-variant exploration, and local/remote analysis workflows.

![nanoMonitor Overview](snapshot.png)

## Key features

- **Real-time monitoring** of BAM or FASTQ/FASTQ.GZ files as they are written by a sequencer.
- **Interactive read metrics** for accuracy, quality score, read length, barcode and amplicon summaries.
- **Targeted amplicon analysis** with primer and barcode support.
- **Structural variant exploration** using whole-genome contact-style plots and region drill-down.
- **Remote analysis mode** using ZeroMQ so heavy analysis can run on another machine while the GUI runs locally.
- **Genomic context** using FASTA, GTF/Tabix and BED-style resources.

## Requirements

- Python `>=3.10,<3.14`
- A project-local virtual environment is recommended.
- Optional: Rust + `maturin` if building the optional Rust acceleration module.

Python dependencies are declared in:

- `pyproject.toml` for `uv`
- `requirements.txt` for standard `pip`/`venv`

Core dependencies include `PyQt6`, `PyQt6-WebEngine`, `pysam`, `numpy`, `scipy`, `plotly`, `matplotlib`, `pyzmq`, `intervaltree`, `edlib`, and `pandas`.

## Installation

See [INSTALL.md](INSTALL.md) for detailed setup notes.

Recommended quick start with `uv`:

```bash
git clone git@github.com:olawa/nanoMonitor.git
cd nanoMonitor

uv python install 3.12
uv venv --python 3.12
source .venv/bin/activate
uv sync
```

Start the GUI:

```bash
python nanoMonitor.py
```

Alternative with standard `venv`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python nanoMonitor.py
```

## Usage

### Local monitoring

Single BAM file:

```bash
python nanoMonitor.py --bam /path/to/sequencing_data.bam --primers primers.tsv
```

Directory monitoring:

```bash
python nanoMonitor.py --dir /path/to/output_fastq_dir --threads 12
```

### Remote analysis mode

Start the server on the analysis machine:

```bash
python ns_server.py --rep-port 5555 --pub-port 5556 --secret my_secure_token
```

Start the monitor locally:

```bash
python nanoMonitor.py --server tcp://server_ip:5555 --secret my_secure_token
```

## CLI arguments

| Argument | Description |
| :--- | :--- |
| `input` | Optional positional input file. |
| `--bam` | Path to an input BAM file. |
| `--primers` | Path to a TSV file containing primer names and sequences. |
| `--genes` | Path to a GTF or BED file for gene models. |
| `--threads` | Number of processing threads, default `8`. |
| `--dir` | Directory for continuous monitoring. |
| `--server` | Remote server address, for example `tcp://10.0.0.1:5555`. |
| `--secret` | Authentication secret for remote server mode. |

## Module architecture

- `nanoMonitor.py`: main application entry point and PyQt6 GUI.
- `python/ns_core.py`: BAM/FASTQ streaming and core metrics.
- `python/ns_workers.py`: PyQt worker threads and remote ZeroMQ client logic.
- `python/ns_amplicon.py`: targeted amplicon analysis.
- `python/ns_rna.py`: RNA/gene counting support.
- `python/ns_structural.py`: structural variant/contact matrix logic.
- `python/ns_visualizer.py`: BAM/FASTQ region viewer.
- `python/ns_resources.py`: FASTA/GTF/BED/primer resource loading.

---

© 2025 Genomics Suite - Optimized for Nanopore Data Monitoring.
