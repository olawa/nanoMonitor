# Installation

This project is a Python/PyQt6 desktop application with optional Rust acceleration.
Use a project-local virtual environment; do **not** install the Python dependencies into Homebrew's global Python.

## Recommended: uv

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

Run against a BAM file:

```bash
python nanoMonitor.py --bam /path/to/sequencing_data.bam --primers primers.tsv
```

Monitor a directory containing BAM/FASTQ files:

```bash
python nanoMonitor.py --dir /path/to/output_fastq_dir --threads 12
```

## Alternative: standard venv + pip

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python nanoMonitor.py
```

## Optional Rust acceleration

The README mentions an optional Rust extension for optimized lane management. If the Rust extension is present in the checkout, install the development dependency group and build it with maturin:

```bash
uv sync --group dev
maturin develop --release
```

If there is no Rust extension source in the checkout, skip this step; the Python fallback remains available.

## macOS notes

Homebrew Python is externally managed, so this will fail and should be avoided:

```bash
pip3 install PyQt6
```

Use `.venv` instead. After activation, `python` and `pip` refer to the virtual environment:

```bash
source .venv/bin/activate
python --version
python -m pip --version
```

For GUI use on macOS, run the application from a normal terminal session with access to the desktop environment.

## Troubleshooting

### `ModuleNotFoundError: No module named 'PyQt6'`

Activate the virtual environment and install dependencies:

```bash
source .venv/bin/activate
uv sync
```

### `ModuleNotFoundError: No module named 'PyQt6.QtWebEngineWidgets'`

Install or resync dependencies; `PyQt6-WebEngine` is required separately from `PyQt6`:

```bash
uv sync
```

or:

```bash
python -m pip install PyQt6-WebEngine
```

### `ModuleNotFoundError: No module named 'edlib'` or `pandas`

These are included in `pyproject.toml` and `requirements.txt`. Resync the environment:

```bash
uv sync
```
