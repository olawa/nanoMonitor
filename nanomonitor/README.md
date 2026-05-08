# nanomonitor (Rust GUI)

Native desktop monitor for nanopore workflows, built with `egui`/`eframe`.

## Run
```bash
cargo run -p nanomonitor
```

## Run with CLI paths
```bash
cargo run -p nanomonitor -- \
  --mode amplicon \
  --input /path/to/sample.bam \
  --primers /path/to/primers.tsv \
  --reference /path/to/ref.fa \
  --gtf /path/to/genes.gtf
```

For monitor-directory startup:
```bash
cargo run -p nanomonitor -- --monitor-dir /path/to/run_folder --start
```

## Build release
```bash
cargo build -p nanomonitor --release
```

## Scope (current prototype)
- Replicated layout of the existing NanoStream monitor:
  - left control rail
  - top filter strip
  - results/log panel
  - bottom analytics plots
- Includes three modes in the UI:
  - Amplicon
  - RNA-Seq
  - WGS + CNV (initial placeholder plot)
- Keeps `nanoparse` as separate CLI executable (command builder integrated in UI).
- Amplicon mode now runs `nanoparse` from the GUI and refreshes:
  - result table
  - length histogram
  - Q-score density
  - accuracy density
- File handling:
  - GUI browse for input BAM/FASTQ file and monitor directory
  - GUI browse for primers/reference/GTF-BED files
  - CLI flags for input/monitor-dir/reference/gtf/primers/nanoparse-bin
- Monitor-directory mode now runs continuously:
  - queues existing BAM/FASTQ files on start
  - watches for newly created files and processes them sequentially
  - tracks processed and failed files in the Run Control panel

Note: current Rust `nanoparse` integration processes BAM input. FASTQ paths can be selected and monitored, but are currently skipped with a log message until FASTQ analysis is added to `nanoparse`.

## Cross-platform
- Target desktop OS: Windows, macOS, Linux.
- Uses `egui`/`eframe` with OpenGL (`glow`) backend to keep runtime dependencies low.
- Remote headless analysis is planned via a separate agent process (see `REBUILD_PLAN.md`).
