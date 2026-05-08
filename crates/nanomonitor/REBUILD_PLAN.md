# nanoMonitor Rust Rebuild Plan

## Goals
- Replace Python/PyQt monitor UI with a native Rust GUI (`egui`/`eframe`).
- Keep `nanoparse` as an independent CLI executable for local and remote pipelines.
- Support Windows, macOS, and Linux with minimal runtime dependencies.
- Preserve current functional workflows (Amplicon, RNA-Seq) and add WGS + CNV view.
- Keep remote execution first-class for headless nodes (secondary in phase 1, hardened in phase 2+).

## Proposed Architecture
- `nanomonitor` (this crate): desktop GUI only.
- `nanoparse` (existing crate): analysis CLI with expandable subcommands.
- `nanomonitor-agent` (new future crate): headless daemon on remote compute nodes.
  - Accepts authenticated requests.
  - Runs `nanoparse` jobs.
  - Streams structured progress and partial results.

## UI Direction (replicate + improve)
- Keep the current mental model:
  - Left rail for resources/run controls/remote/session actions.
  - Top filter strip with quick toggles and run actions.
  - Center results table + log tab.
  - Bottom charts (accuracy, Q-score, length).
- Improvements:
  - Reduce control clutter with grouped cards and progressive disclosure.
  - Keep critical controls always visible (mode, start/stop, filters, remote status).
  - Move advanced actions into compact menus (variant tuning, export options, debug controls).
  - Mode-specific panels:
    - Amplicon: current table and length/QS summaries.
    - RNA-Seq: target coverage/gene expression panel.
    - WGS: CNV panel with chromosome navigation and segment calls.

## Execution Plan

### Phase 0 - Foundation (current)
- Workspace includes both `nanoparse` and `nanomonitor`.
- `nanomonitor` egui shell mirrors existing screen layout.
- `nanoparse` command builder integrated in GUI state.

### Phase 1 - Functional Parity (local)
- Implement file load and monitor directory workflows.
- Parse and render real `nanoparse` outputs (JSON schema).
- Replace placeholder charts with real data pipelines.
- Implement robust table interactions (sorting, filtering, multi-select, export).

### Phase 2 - Remote Workflows
- Define stable JSON protocol (no pickle).
- Build `nanomonitor-agent` with auth token + job lifecycle.
- Add reconnect/resume behavior and node capability checks.
- Stream partial metrics so remote mode feels identical to local.

### Phase 3 - WGS + CNV
- Add WGS mode data contracts in `nanoparse`.
- CNV plots:
  - binned log2 ratio
  - segmentation overlays
  - gene/region tooltips
- Support chromosome jumping and copy-number threshold presets.

### Phase 4 - Hardening
- End-to-end integration tests (local + remote).
- Packaging:
  - Windows MSI or zip
  - macOS app bundle
  - Linux AppImage/tarball
- Telemetry/logging and crash-safe recovery.

## Recommended Technical Improvements
- Shared typed models:
  - Put request/response structs in a small shared crate (`nanomonitor-protocol`).
- Process management:
  - Use async job supervisor for subprocess lifecycle, timeouts, cancellation, restart.
- Data flow:
  - Keep per-mode immutable snapshots + diff updates to avoid UI lockups.
- Security:
  - Token auth + optional mTLS for cluster installs.
  - No code-loading serializers; use JSON or MessagePack.
- Performance:
  - Incremental downsampling for large plots.
  - Background aggregation for >10M read sessions.

## Immediate Next Steps
1. Define `nanoparse` output schema v1 (Amplicon + RNA-Seq + WGS placeholder).
2. Implement local process runner in `nanomonitor` and bind real results table.
3. Replace placeholder plots with computed distributions from live results.
4. Freeze remote protocol draft before building `nanomonitor-agent`.
