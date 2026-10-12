# BCR-ABL1 mode for nanoMonitor — plan

Status: plan only, nothing implemented. Written 2026-10-12 after a full review of the workspace
(see "Prerequisites") and a read of nagini's BCR-ABL path (`~/dev/cgu_pipelines/nagini`).

## Decisions (confirmed by Ola)

1. GUI is Rust/egui (`crates/nanomonitor`). The PyQt app is not extended.
2. nanoMonitor does the whole front half: MinKNOW FASTQ -> filter (length, qs) -> **splice
   alignment** -> aligned BAM per barcode, which is handed to nagini as `--input-type aligned_bam`.
3. nanoMonitor runs on the sequencing computer, next to MinKNOW (local filesystem, no remote layer).
4. Reference build is hg19, as nagini's `bcr_abl` assay (`config.yaml:40-48`).

## Why the handoff is safe (what nagini does with an aligned BAM)

Verified in `stage_aligned_bam.smk` and `nagini_routing.py`:

- `aligned_bam` input is **not realigned and not read-filtered**: it only runs
  `samtools addreplacerg -m overwrite_all`, `sort`, `index`. The assay forces `read_filter` on
  only for the FASTQ/unmapped-BAM paths.
- So the BAM we hand over must already be exactly what nagini would have produced itself:
  - filtered with the assay's settings (`read_filter`: `min_qv 20` against dorado's `qs`,
    `length_range 500-3000`) using the same `nanostream filter` code;
  - aligned with `minimap2 -a -x splice -y -Y --MD` against the same hg19 FASTA nagini uses
    (`/data/references/hg19/hg19.fa` on the workstation profile). A prebuilt index silently keeps
    its own `-k/-w`, so a `.mmi` must be built with `-x splice` or we align against the FASTA.
- minimap2 aligns each read independently, so chunked/incremental alignment gives the same
  alignments as one pass over the whole FASTQ. This is what makes a golden test possible (below).
- nagini rewrites `@RG` itself, so we do not need to match its read-group IDs, but we must
  keep `MM/ML`-free standard tags and `qs`/`ch` comments (`-y` copies FASTQ header comments).
- Provenance: the engine writes `provenance.json` (nanostream git commit, minimap2 version +
  args, reference path + `.fai` checksum, filter parameters, input file list). nagini-side
  verification of that file is a follow-up for the nagini thread, not a blocker.

## Architecture

```
MinKNOW fastq_pass/barcodeXX/*.fastq.gz   (already demultiplexed by MinKNOW/dorado)
   |  watcher: file is "ready" when size+mtime are stable for N s (or a newer chunk exists);
   |           ledger of processed files + byte counts persisted -> resumable, never re-read
   v
filter   (nanofilter-core, same call path as `nanostream filter`; qs from header)
   v
align    (minimap2 -a -x splice -y -Y --MD, see "Alignment process" below)
   v
segments  <outdir>/barcodeXX/seg_000123.bam   (append-only, unsorted, one per batch)
   |                       |
   |                       +-- live stats (Rust): fusion isoform, Branford-site depth/VAF
   v
snapshot  every ~10 min: merge+sort+index segments -> snapshots/barcodeXX.bam(.bai)  (tmp + rename)
   v
prelim analysis (nagini's scripts or Rust, see Phase 2) -> state.json, prelim/*.tsv|json
   v
GUI reads only <outdir>  (state.json, prelim/, logs)             [Start full analysis]
                                                                      |
   final: merge all segments -> alignment/barcodeXX.bam, write barcode map, run
   start_nagini.sh --input-type aligned_bam --assay bcr_abl --library-type cdna \
       --reference-build hg19 --barcode-map ... --analysis-dir ...
   v
GUI shows results/reports/bcr_abl/*, summary xlsx, and prelim vs final comparison
```

- The engine is a separate headless process (`nanostream bcr-abl watch`), started and observed by
  the GUI but not owned by it. A GUI crash must not stop alignment, and the same engine can be run
  over ssh. GUI and engine communicate through the output directory (`state.json`, ledger, logs).
- New library crate (working name `nanobcr-core`) holds: watcher, ledger, align runner, segment
  store, live counters, state schema. `nanostream` gets the `bcr-abl` subcommand, `nanomonitor`
  gets the tab. Core crates stay free of clap and of `build_global()` rayon calls.

### Alignment process (the one real technical risk)

Whole-genome hg19 splice indexing takes minutes and several GB, and minimap2 on stdin waits for a
full mini-batch (`-K`, 500 Mb default) before emitting anything. Options, decide by benchmark in
Phase 1 (gate: new reads appear in a segment within ~30 s of being written by MinKNOW, steady
state keeps up with the run's output rate):

- **A (preferred):** one long-lived minimap2 per run, FASTQ fed over stdin, small `-K`
  (a few Mb), SAM on stdout parsed/written by the engine. Index loaded once.
- **B:** a prebuilt `-x splice` `.mmi`, minimap2 started once per batch (index load from page
  cache is seconds).
- **C (rejected for now):** align against chr9+chr22 only. Faster, but off-target reads would
  map differently from nagini's whole-genome alignment and break equivalence.

`minimap2` and `samtools` must be on PATH on the sequencing computer; the engine checks versions
and records them in provenance. No in-process aligner in this plan.

## Phases

### Phase 0 — prerequisites (before any new code)

- Split the working tree into three commits, order A -> B -> C (each builds):
  A pe refactor (`pe.rs`, `nanofilter-core/lib.rs`, `nanostream/lib.rs`, two `fastq.rs` hunks);
  B `--top-n` (`matcher.rs`, `nanoparse/main.rs`, `nanostream/{cli,main}.rs`, `nanomonitor/app.rs`);
  C qs/MinKNOW headers (`header.rs`, `filter.rs`, `bam.rs`, three `fastq.rs` hunks, the
  `qs:f:` parse + tests in `matcher.rs`). `matcher.rs` needs hand-split patches. Do this only
  after the nagini thread has stopped editing the tree.
- Fix before relying on the filter path: `QUAL=*` overflow (`bam.rs:617`, `pe.rs:246`);
  basecall model truncation without barcode in `RG` (`header.rs:72`).
- Not blockers for this feature but tracked: primer-matching bugs in `matcher.rs`
  (`semi_global_align`, `is_fwd_rev_pair`, tid-vs-name, supplementary records) — they affect
  nagini's `filter_bcr_abl_amplicons.py` today and should be reported to the nagini thread.
- GUI bugs that must be fixed before the tab is built on them (all in `app.rs`):
  Progress overwriting `self.data` (~1104), `create_worker_sender_if_none` replacing the
  analysis channel (~613), per-callback full-state clone (matcher.rs ~1367/1687), and
  snapshot/variant regions built from amplicon names.
- Decide what happens to `ns_server.py` (`pickle.loads` before auth, binds `tcp://*`): bind to
  localhost + JSON, or delete. Independent of this feature, but do not leave it.

### Phase 1 — engine: watch, filter, align, segments (no GUI)

- Watcher + ledger + readiness rule; handles `fastq_pass` and `fastq_fail` (default both, as
  nagini does; configurable), `.fastq` and `.fastq.gz`, and a growing last file.
- Filter + align + segment writer; `provenance.json`; `state.json` (per barcode: reads seen /
  passed filter / aligned / primary-mapped, bases, last file, last-update time).
- `--replay <finished run dir> --speed N` feeds existing files in at MinKNOW's cadence. This is
  how everything below is tested without a sequencer.
- Barcode map: same CSV as nagini (`barcode, sample, type`), plus MinKNOW `sample_sheet*.csv`
  auto-detection like `start_nagini.sh`.
- **Golden test (the acceptance test for the whole handoff):** take a finished run
  (e.g. `bcr_abl_exp7_9`, a few barcodes), run it (a) through nagini's own FASTQ path and
  (b) through replay -> engine -> `aligned_bam`; require identical primary alignments
  (compare `samtools view` flagstat + per-read primary position/CIGAR hash) and identical
  `fusion_summary.tsv`, `branford_mutations.tsv`, `phasing.tsv`. Also run nagini `--dry-run`
  against the engine's output to prove `prepare_nagini_inputs.py` accepts it (open question:
  whether `aligned_bam` + `--barcode-map` over a directory maps barcodes to samples; if not,
  the engine writes `samples.tsv`/`units.tsv` itself or nagini gets a small change — raise with
  the nagini thread).
- Benchmark gate for the alignment-process choice above.

#### Phase 1 status (2026-10-12): implemented in `crates/nanobcr-core`, CLI `nanostream bcr-abl watch|snapshot|replay`

Done and verified:
- Watcher (mtime + size stability, truncated gzip retried, `fastq_fail` optional), JSONL ledger
  (resume after restart, torn last line tolerated), filter via the same
  `fastq_record_passes` predicate as `nanostream filter`, one `minimap2 -a -x splice -y -Y --MD`
  run per cycle with in-process per-barcode SAM split and counting, segments compacted
  behind an atomic manifest, `state.json`, `provenance.json`, `replay`.
- **Golden test on real data** (31,110 reads sampled from the real hg19 BCR-ABL barcode24 BAM,
  written as a MinKNOW-style run, replayed live in 8 half-written chunks): the engine's BAM
  (4 cycles, aligned against a `-x splice` `.mmi`) has **byte-identical records, every field and tag,
  to `nanostream filter` + one `minimap2 -x splice` run on the FASTA** — 5,039 records, same md5.
  The synthetic-genome version of this comparison runs in `cargo test`.
- **Benchmark gate (hg19, 84 contigs, 3.1 GB FASTA): passed, with a memory cost.** Index build
  35 s / 23 GB RSS (once); the splice `.mmi` is 11.6 GB; each cycle pays ~4 s index load + ~3 s
  teardown with 14–18 GB peak RSS; alignment of unfiltered long reads (median 4.4 kb) ran ~150 reads/s
  on 8 threads, so the filtered 0.5–3 kb amplicon reads are far cheaper. Per-batch invocation (option B)
  is therefore fine and the persistent-process variant (A) is not needed. **The sequencing
  computer needs >= ~24 GB free RAM** while a cycle runs; a smaller reference would cut this, but
  restricting to chr9/chr22 would break equivalence with nagini (see option C).
- nagini accepts the handoff: `prepare_nagini_inputs.py --input-type aligned_bam --barcode-map map.csv`
  over a directory of `barcodeNN.bam` (+`.bai`) produced `units.tsv`/`samples.tsv` with
  `input_type=aligned_bam`, `assay=bcr_abl`, `library_type=cdna`. So the Phase 4 handoff layout is
  `<handoff>/barcodeNN.bam(.bai)` + a barcode map CSV (`barcode,sample,type,assay,library_type`).
  (A full `snakemake` dry run was not possible here: the workstation hg19 paths do not exist on this Mac.)

Known gaps / follow-ups:
- Not yet run against a live MinKNOW run; only replay. Verify the stability rule against real file
  rotation, and that MinKNOW's own `final_summary_*.txt` location matches `run_finished`.
- A file modified after it was processed is not re-read (a warning is not yet logged).
- A `.mmi` as `--reference` leaves no `.fai` fingerprint in `provenance.json`; pass the FASTA.
- samtools truncates multi-word FASTQ comments to the first tag (`qs:f:` survives, `ch:i:`/`st:Z:` do
  not); nagini's own path behaves the same, so this is equivalent, but those tags are not in the BAM.
- Segment/compaction cost grows with run length (each compaction merges into the whole BAM); revisit
  if multi-day runs make snapshots slow.

### Phase 2 — preliminary analysis

- Snapshot job every ~10 min (configurable) per barcode.
- Two sources, kept clearly labelled in the UI:
  - **Live (Rust, per batch, cheap):** per-read fusion isoform (e13a2/e14a2/e1a2/…) from the
    primary + SA alignments against `bcr_abl1_exons_hg19.bed`; depth and raw alt fraction at
    the 73 Branford sites from CIGAR/MD; reads spanning BCR-ABL junction; target depth. Marked
    "indicative".
  - **Authoritative prelim (nagini's scripts on the snapshot):** `bcr_abl_fusion_counter.py`
    -> `bcr_abl_branford_caller.py` -> `bcr_abl_phaser.py`, run as a subprocess in a
    configured Python env, outputs under `prelim/`. `bcr_abl_exon_coverage.py` is too slow for
    repeated runs: replace with a Rust/htslib per-exon depth, or run it only on the final BAM.
- Read all resources from nagini (`workflow/resources/bcr_abl/*`, thresholds in `config.yaml:423-430`)
  via a configured nagini path; do **not** copy them into nanoMonitor. Retire nanoMonitor's own
  hg38 `clinical_mutations.tsv` and hg19/hg38-mixed `KNOWN_FUSIONS` for this mode.
- Keep history (one row per snapshot) so the GUI can show VAF / isoform fraction vs reads.
- Regression test: live Rust numbers vs script outputs on the golden dataset, within stated
  tolerances; differences documented rather than hidden.
- "Enough data" indicator per sample: Branford-site min depth and fusion `PRIMARY` status
  (>=5 reads) from nagini's thresholds; MinKNOW `final_summary_*.txt` marks run end.

#### Phase 2 status (2026-10-12): implemented in `nanobcr-core::prelim`; `watch --nagini DIR` and `bcr-abl analyze`

Decision: **the preliminary analysis runs nagini's own scripts; no Rust re-implementation.** On every new
snapshot a background worker (never blocking an alignment cycle; queued snapshots of one barcode are
coalesced to the newest) runs `bcr_abl_fusion_counter.py`, `bcr_abl_branford_caller.py` and
`bcr_abl_phaser.py` with the arguments, resources and thresholds from nagini's own `bcr_abl.smk` and
`config.yaml`, on a private hard-linked copy of the snapshot. It writes `<out>/<barcode>/prelim/`
(the full script outputs, `summary.json` for the GUI, and one line per snapshot in `history.jsonl`),
and `state.json` carries `prelim_generation/status/unix/error/summary` per barcode. A failing stage
is reported (`status: partial`) without losing the others. `summary.json` also has a `sufficiency`
block: a PRIMARY fusion call and >= 95 % of Branford sites at nagini's `min_depth`; its inputs are recorded
so the heuristic can be judged. `bcr_abl_exon_coverage.py` is not run live (slow; final pipeline only).

The Rust per-read live counters planned for this phase were dropped: on the real barcode24 BAM
the fusion counter takes 1.7 s and the phaser 7 s even at 637k records; only the Branford caller
scales with depth (0.5 s at 2.7k reads, 92 s at ~46,000x). With a snapshot every 2-5 minutes the
scripts keep up, and a second implementation would only add a way for the live and final numbers to
disagree. The engine's own per-cycle counters (reads, filter outcome, mapped/supplementary) stay live.

Verified on the real data, replayed live with a snapshot every 20 s: generation 1 (329 reads) already
calls e1a2 (220 fusion reads); generation 2 (1,869 reads) shows Y253H 2.40 % and T315I 2.36 %; the final
generation (2,742 reads) gives Y253H 2.46 %, T315I 2.35 %, CIS phasing 90 %, identical to running the scripts
by hand on the final BAM. Tests use a stub nagini (canned script outputs, recorded argv) so the plumbing,
nagini-config passthrough, partial failure and history are covered without pysam or clinical data.

Open: the Branford step can take minutes on a full-depth run; the worker keeps up only if the snapshot
interval is not much shorter than that. The GUI should show "analysing generation N" from `prelim_status`.

### Phase 3 — egui BCR-ABL mode

- Run setup: MinKNOW run dir (or replay dir), barcode map, reference FASTA, nagini path,
  output dir; preflight checks (minimap2/samtools versions, FASTA + `.fai`, nagini script env,
  disk space) shown before start.
- Live view: per-barcode table (reads, passed filter, aligned, last file age), isoform fractions,
  KD mutation list with VAF/depth/status and TKI relevance, phasing pairs, per-exon depth, a
  time-series of the above, engine log tab. Staleness indicators (engine heartbeat age).
- Add `egui_extras` image loaders so `fusion_plot.png` renders. Reuse existing density plots
  for filtered read length/qs.
- Engine control: start/stop/attach to a running engine via `state.json` + a pid/lock file.

### Phase 4 — start full pipeline and show results

- "Run nagini" button: final snapshot -> merge segments -> `alignment/barcodeXX.bam` + `.bai` ->
  assemble the `start_nagini.sh` command (shown to the user before launch) -> spawn in its own
  process group, stdout/stderr to `nagini_run.log`.
- Progress from `.snakemake/log/*.snakemake.log` ("N of M steps (x%) done"), per-rule logs for
  failures; Stop = SIGINT to the process group, then `--unlock` offer. Success = exit 0 **and**
  `summary_report.xlsx` present (profiles use `keep-going`, so failures arrive at the end).
- Results view: `results/reports/bcr_abl/{s}_{t}.fusion_summary.json`, `branford_mutations.tsv`,
  `phasing.json`, `exon_coverage.tsv`, `fusion_plot.png`; link/open the xlsx workbooks and
  MultiQC html; table comparing prelim vs final per sample, flagging any change in isoform or
  mutation calls.

## Not in scope / deferred

- In-process alignment (no minimap2 subprocess), remote/ZeroMQ monitoring, hg38 BCR-ABL.
- Changing nagini (except the possible `aligned_bam` + barcode-map and `filter_bcr_abl_amplicons.py`
  config-path items, which go to the nagini thread).
- Fixing the unrelated findings (UMI strand handling, consensus orientation, barcode panics,
  sequencing-summary relative time window). Tracked in the review; not part of this feature.

## Open risks

- minimap2 batch latency/index load (benchmark gate in Phase 1).
- Partially written / rotated MinKNOW files (stable-size rule + ledger; verify against a real
  live run, not only replay).
- nagini `aligned_bam` + barcode map behaviour (verify in Phase 1 with `--dry-run`).
- Memory growth in live counters over a multi-day run (bounded histograms, no per-read vectors).
- Nagini's Python/containers on the sequencing computer (the prelim scripts run in a configured
  env; if unavailable, the Rust live view still works and the UI says so).
