# Rust Architecture Review

## Current Shape

The Rust workspace is organized as six crates:

- `nanoseq-core`: shared low-level nanopore utilities such as FASTQ path detection, header parsing, quality helpers, pore ranges, and time windows.
- `nanoparse-core`: analysis logic for statistics, enrichment, amplicon matching, primer parsing, QV helpers, output, and pore idle-time statistics.
- `nanoparse`: CLI wrapper around `nanoparse-core`.
- `nanofilter-core`: transformation/filtering logic for BAM/FASTQ filtering, barcode split/discovery, UMI detection, clustering, consensus, and reports.
- `nanofilter`: CLI wrapper around `nanofilter-core`.
- `nanomonitor`: native GUI that already links `nanoparse-core`, `nanofilter-core`, and `nanoseq-core` directly.

This is mostly a healthy shape: core libraries contain reusable behavior, thin binaries provide CLI entry points, and the GUI can call library code without shelling out.

## Main Issue

The user-facing Rust surface is fragmented. A viewer or CLI user has to know whether a command lives in `nanoparse` or `nanofilter`, even though both operate on the same BAM/FASTQ sequencing data and share metadata concepts through `nanoseq-core`.

That fragmentation is starting to show up in small usability and maintenance issues:

- `pore-stats` was naturally discussed as possibly belonging to either `nanoparse` or `nanofilter`.
- time slicing exists in `nanofilter`, while timing QC now exists in `nanoparse`.
- `nanomonitor` still has UI language and command previews centered on `nanoparse`, even though it already depends on both core libraries.
- the old Python monitor shells out to `nanoparse`, creating another binary discovery path that will need updating if the primary tool changes.

## Recommended Target

Create one user-facing Rust binary, tentatively named `nanostream`, while keeping the core crates separate for now.

Recommended CLI layout:

```text
nanostream stats <input>
nanostream amplicons <input> --primers primers.tsv
nanostream enrichment <input> --bed targets.bed
nanostream pore-stats [input] --sequencing-summary sequencing_summary.txt
nanostream filter <input> --output filtered.bam
nanostream extract <input> --channel-range 1-512 --output subset.bam
nanostream split <input> --barcodes barcodes.tsv --output-dir out
nanostream discover <input>
nanostream make-pe <input>
nanostream umi <input> --output-dir out
nanostream monitor
```

This gives the viewer and CLI users one executable, one help tree, and one stable integration point.

## What To Keep Separate

Do not merge `nanoparse-core` and `nanofilter-core` immediately. They have different responsibilities:

- `nanoparse-core` answers analysis questions: what is in the reads?
- `nanofilter-core` produces transformed read outputs: which reads should be kept, split, or converted?
- `nanoseq-core` should continue to absorb shared primitives that both need.

The better near-term move is a unified binary crate that depends on both cores. Core merging can be revisited later if shared code keeps growing.

## Refactor Plan

1. Add a new binary crate `nanostream-cli` or rename one existing binary crate to own the combined CLI.
2. Move the `nanoparse/src/main.rs` and `nanofilter/src/main.rs` command definitions into modules under the combined binary.
3. Preserve existing subcommands and flags first, so current scripts keep working with minimal changes.
4. Add compatibility binaries or aliases for `nanoparse` and `nanofilter` during a transition period.
5. Update `nanomonitor` command naming from `nanoparse` to the unified binary where CLI preview or binary discovery is still used.
6. Update the Python `NanoparseWorker` discovery path to prefer the unified binary, while accepting `nanoparse` as a fallback.
7. Move shared sequencing-summary parsing, BAM/FASTQ timing metadata, and file format helpers toward `nanoseq-core` as they stabilize.

## Cleanup Opportunities

- Move sequencing summary parsing from `nanoparse-core::pore_stats` into `nanoseq-core` if it will also support filtering, slicing, or QC elsewhere.
- Move BAM/FASTQ input classification into a shared enum in `nanoseq-core` instead of repeatedly checking extensions.
- Remove `clap` from core crates where possible. CLI parsing should live in binary crates; core crates should expose plain config structs.
- Consider standardizing on one BAM stack. `nanoparse-core` uses `rust-htslib`; `nanofilter-core` uses `noodles`. Keeping both may be fine, but it should be a deliberate choice because it increases dependency and maintenance surface.
- Split very large files over time: `nanomonitor/src/app.rs`, `nanoparse-core/src/matcher.rs`, `nanoparse-core/src/pore_stats.rs`, and several `nanofilter-core` modules are large enough that feature work will get harder.

## Suggested First Implementation Step

The lowest-risk first step is to add a new `nanostream` binary crate that re-exports the current command surfaces:

- `nanostream parse stats|amplicons|enrichment|pore-stats` or direct top-level aliases.
- `nanostream filter`, `nanostream split`, `nanostream umi`, and the other current `nanofilter` commands.

Once that builds and tests pass, update GUI/Python integrations to use `nanostream` by default. Only after that should we decide whether to retire or keep the old binary names.

## Follow-up Investigation

`~/dev/methylartist/methylartist-rs/crates/genomics-bam-core` is useful, but it is not a BAM backend today. It intentionally avoids choosing between `rust-htslib`, `noodles`, or another reader. It currently provides shared scan orchestration pieces: regions, windows, ordered chunk execution, thread options, and runtime summaries.

That means it could help standardize future region/window/chunk processing, but it does not by itself solve the "one BAM path" goal. For that, we still need to choose or build a concrete BAM reader abstraction.

Current BAM state:

- `nanoparse-core` uses `rust-htslib`.
- `nanofilter-core` uses `noodles`.
- `nanostream` currently calls both cores, so both stacks are still present.
- BAM filtering/export in `nanofilter-core` now has threaded reader/writer entry points for the unified CLI.

Recommended BAM direction:

- Prefer `rust-htslib` for the shared read-analysis path if practical, because it already supports threaded BAM IO in the current parser/statistics code and tends to be straightforward for tag-heavy ONT BAMs.
- Keep `noodles` where it is already working for write-heavy filtering until we have tests around equivalent `htslib` output.
- Use `genomics-bam-core` later for chunk/window orchestration once we have a shared record adapter.

Current FASTQ state:

- FASTQ filtering is still mostly sequential.
- FASTQ split and UMI paths already use worker threads.
- `pigz` is available on this machine at `/opt/homebrew/bin/pigz`; using it as an optional external gzip backend is a good next optimization for `.fastq.gz` read/write paths.
