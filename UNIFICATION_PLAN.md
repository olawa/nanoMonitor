# NanoStream / Nanofilter Unification Plan

## Goal

Build one Rust-backed processing stack that serves:

- `nanofilter`-style CLI filtering, demultiplexing, conversion, pore/time splitting
- `nanoparse`-style CLI analysis, amplicon calling, BED enrichment
- `nanomonitor` GUI using the exact same library code paths as the CLI

The GUI must not reimplement filtering or analysis logic. It should call shared library functions directly.

## Non-Goals

- Do not force an immediate rewrite of all BAM handling onto one backend.
- Do not merge all logic into one large crate.
- Do not allow GUI behavior to drift from CLI behavior.

## Canonical Quality Definition

This becomes the shared definition of per-read QV everywhere:

1. If BAM tag `QS` is present, use it.
2. Otherwise compute:

```python
-10 * log10(mean(10 ** (-q / 10) for q in quals))
```

This matches `nanofilter` and should replace the arithmetic mean currently used in `nanoparse`.

## Current State

### `nanoparse`

Primary strengths:

- Primer loading and k-mer indexing
- Amplicon matching
- BED enrichment and pore statistics
- JSON output for GUI consumption

Current files:

- `nanoparse/src/main.rs`
- `nanoparse/src/stats.rs`
- `nanoparse/src/matcher.rs`
- `nanoparse/src/primers.rs`
- `nanoparse/src/enrichment.rs`

Current issues relevant to unification:

- Duplicated FASTQ parsing logic in `stats.rs` and `matcher.rs`
- Wrong QV definition for filtering/stats
- `run_stats` loads all BAM records into memory

### `nanofilter`

Primary strengths:

- Streaming FASTQ/BAM filtering
- BAM to FASTQ and FASTQ to BAM conversion
- Barcode split and barcode discovery
- Nanopore header parsing
- Better QV semantics

Relevant source files:

- `/Users/olwal516/dev/nanofilter/src/main.rs`
- `/Users/olwal516/dev/nanofilter/src/filter.rs`
- `/Users/olwal516/dev/nanofilter/src/header.rs`
- `/Users/olwal516/dev/nanofilter/src/barcode.rs`
- `/Users/olwal516/dev/nanofilter/src/bam.rs`
- `/Users/olwal516/dev/nanofilter/src/fastq.rs`

Current issues relevant to unification:

- Mixed domain scope inside one binary
- BAM backend differs from `nanoparse`
- GUI does not consume this logic directly

## Target Workspace Layout

Keep a Rust workspace with separate crates:

1. `nanoseq-core`
2. `nanofilter-core`
3. `nanoparse-core`
4. `nanofilter` CLI
5. `nanoparse` CLI
6. `nanomonitor` GUI

### Proposed responsibilities

#### `nanoseq-core`

Shared low-level and cross-cutting utilities.

Suggested modules:

- `quality`
- `format`
- `header`
- `tags`
- `filters`
- `split`
- `intervals`

Suggested API areas:

- file format detection
- shared QV calculation
- `QS`, `cm`, `dx`, and time-tag accessors
- nanopore FASTQ header parsing
- common filter structs
- shared pore/time range parsing
- BED/interval helpers

This crate should be backend-neutral where practical.

#### `nanofilter-core`

Streaming file transformation and demux.

Suggested modules:

- `filter_fastq`
- `filter_bam`
- `convert`
- `barcode`
- `discover`
- `split_by_barcode`
- `split_by_pore`
- `split_by_time`

This remains allowed to use `noodles` internally.

#### `nanoparse-core`

Mapped-read and primer-based analysis.

Suggested modules:

- `primers`
- `amplicons`
- `stats`
- `enrichment`

This remains allowed to use `rust-htslib` internally.

#### `nanomonitor`

GUI only.

Responsibilities:

- choose files/directories
- configure filters and analysis options
- execute library operations
- present tables, plots, logs, and export actions

No duplicate filtering or parsing logic should live here.

## Shared Data Model

Define shared config types in `nanoseq-core`.

### `ReadFilter`

```rust
pub struct ReadFilter {
    pub min_qv: Option<f64>,
    pub min_len: Option<usize>,
    pub max_len: Option<usize>,
    pub duplex_only: bool,
    pub pore_range: Option<PoreRange>,
    pub time_window: Option<TimeWindow>,
    pub max_reads: Option<usize>,
}
```

### `PoreRange`

```rust
pub struct PoreRange {
    pub min: i64,
    pub max: i64,
}
```

### `TimeWindow`

```rust
pub struct TimeWindow {
    pub start: chrono::DateTime<chrono::FixedOffset>,
    pub end: chrono::DateTime<chrono::FixedOffset>,
}
```

### `SplitMode`

```rust
pub enum SplitMode {
    Barcode,
    PoreBuckets(Vec<PoreRange>),
    TimeBuckets(Vec<TimeWindow>),
}
```

## Time Handling

Support both BAM and FASTQ.

Priority order:

1. BAM tag with explicit start/read time if present
2. FASTQ header fields such as `start_time`
3. Fail with a clear error if time splitting is requested but no time metadata exists

This must be implemented once in shared code, not independently in CLI and GUI.

## Immediate Extraction Candidates

These should move first because they are duplicated or should be canonical.

### Move into `nanoseq-core`

From `nanofilter`:

- `calculate_phred_avg` from `/Users/olwal516/dev/nanofilter/src/filter.rs`
- `parse_nanopore_header` and `NanoporeMetadata` from `/Users/olwal516/dev/nanofilter/src/header.rs`
- `reverse_complement` byte-oriented helper from `/Users/olwal516/dev/nanofilter/src/barcode.rs`
- tag helpers currently embedded in `/Users/olwal516/dev/nanofilter/src/bam.rs`

From `nanoparse`:

- FASTQ path detection duplicated in `nanoparse/src/stats.rs` and `nanoparse/src/matcher.rs`
- BED interval parsing and merging from `nanoparse/src/enrichment.rs`

### Keep in `nanoparse-core`

- `nanoparse/src/primers.rs`
- amplicon matching logic from `nanoparse/src/matcher.rs`
- analysis JSON result structs

### Keep in `nanofilter-core`

- barcode matching and discovery from `/Users/olwal516/dev/nanofilter/src/barcode.rs`
- BAM/FASTQ split pipelines
- unaligned BAM writer
- synthetic paired-end generation

## Recommended Sequence

### Phase 1: Standardize behavior

1. Create `nanoseq-core`.
2. Add shared QV calculation.
3. Add shared FASTQ/header/tag utilities.
4. Refactor `nanoparse` stats and FASTQ filtering to use shared QV logic.
5. Leave all CLI behavior otherwise unchanged.

Acceptance criteria:

- `nanofilter` and `nanoparse` return the same per-read QV for the same sequence/qualities
- BAM `QS` is preferred consistently
- no GUI changes yet

### Phase 2: Shared filtering primitives

1. Add shared `ReadFilter`, `PoreRange`, and time parsing.
2. Refactor `nanofilter` filtering to consume shared filter structs.
3. Add `cm` pore filtering to `nanofilter-core` pipelines.
4. Add time-window filtering to `nanofilter-core`.

Acceptance criteria:

- CLI can filter and split by pore and time using shared config types
- filter semantics are identical across BAM and FASTQ where metadata exists

### Phase 3: Workspace split

1. Rename current `nanoparse` crate internals to `nanoparse-core`.
2. Create thin `nanoparse` CLI wrapper.
3. Import `nanofilter` into this workspace as:
   - `nanofilter-core`
   - `nanofilter` CLI
4. Keep backend-specific code inside each core crate.

Acceptance criteria:

- all binaries build from one workspace
- no cross-project copy-paste remains for shared utilities

### Phase 4: GUI unification

1. Extend `nanomonitor` to call:
   - `nanofilter-core` for filtering/splitting/conversion
   - `nanoparse-core` for stats/amplicons/enrichment
2. Add GUI panels for:
   - QV/length filtering
   - barcode split
   - split by pore
   - split by time
   - amplicon analysis
   - enrichment over BED
3. Make GUI actions pure wrappers over core library calls.

Acceptance criteria:

- same input and options produce the same outputs from GUI and CLI
- no GUI-only logic for filtering/demux/analysis

## Backend Strategy

Do not unify on `rust-htslib` or `noodles` immediately.

Short-term rule:

- `nanofilter-core` keeps `noodles`
- `nanoparse-core` keeps `rust-htslib`
- `nanoseq-core` stays backend-neutral where possible

Revisit backend unification only if:

- performance profiling shows a concrete bottleneck
- duplicated backend adapters become hard to maintain
- one backend blocks a needed feature

## Risks

### Risk 1: Semantic drift during migration

Mitigation:

- freeze QV and filter semantics first
- add fixture-based tests with the same expected outputs in both binaries

### Risk 2: FASTQ and BAM time metadata differ

Mitigation:

- define a single normalized time extraction API
- explicitly report unsupported cases

### Risk 3: GUI feature growth outruns library cleanup

Mitigation:

- do not add GUI actions until the corresponding library API exists

## First Implementation Slice

This is the first slice worth doing immediately:

1. Add `nanoseq-core` to this workspace.
2. Implement:
   - shared QV calculation
   - shared `QS` accessor
   - shared FASTQ path detection
   - shared nanopore header parser
3. Refactor `nanoparse/src/stats.rs` and `nanoparse/src/matcher.rs` to use shared QV logic.
4. Add regression tests proving the QV calculation is identical to `nanofilter`.

This gives immediate correctness wins without forcing a backend rewrite.

## Suggested Near-Term Backlog

1. Phase 1: `nanoseq-core` with quality/header/tag utilities
2. Fix `nanoparse` stats to stream rather than collect all BAM reads
3. Import `nanofilter` as workspace crates
4. Add shared pore/time split model
5. Expose filter/split operations in `nanomonitor`
