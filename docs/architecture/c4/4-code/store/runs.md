# Runs: code

C4 level 4 (code) of the component **Runs** (`store.runs`) in the container **nyx-store** (crate `nyx_store`). [Component view of the container](../../3-components/store.md) · [All code pages](../README.md)

> Append-only run records in runs.jsonl, the digest derived from them and the filter language of the runs table.

**Decisions:** [ADR 64](../../../adr/0064-identity-hashes-one-persisted-the-rest-are-caches.md)  
**Size:** 1 file, 400 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-store/src/runs.rs` | 400 | Run records: one append-only line per run in `runs.jsonl`, the digest derived from them, and the small filter language of the runs table. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| container `workspaceStore` (Workspace data) |  (declared: Appends runs.jsonl) |  |
| [Data model and card](./model.md) | model | Digest, Range |
| [Utilities](./util.md) | util | io_err, run_suffix |

## Used by

| Component | Modules | Items |
|---|---|---|
| [engine: Store link](../engine/storeLink.md) | runs | Artifact, RunRecord, RunStatus |
| [Workflow folder](./flow.md) | runs | RunPatch, RunRecord |
| [Search](./search.md) | runs | RunRecord, parse_quantity, read |
| [Workspace](./workspace.md) | runs | RunRecord, RunStatus |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class RunStatus {
      <<enumeration>>
      Ok
      Failed
      Cancelled
  }
  class Artifact {
      <<struct>>
      +kind String
      +path String
      +bytes u64
  }
  class RunRecord {
      <<struct>>
      +id String
      +t String
      +sec f64
      +status RunStatus
      +spec String
      +femcore String
      +version u32
      +params BTreeMap
      +more_6_fields
  }
  class RunPatch {
      <<struct>>
      +star Option~bool~
      +tags Option
      +note Option~String~
  }
  class Filter {
      <<struct>>
      +parse(text) Result
      +is_empty() bool
      +matches(r) bool
  }
  class fns_runs_rs {
      <<module runs.rs>>
      +append(path, r) Result
      +append_patch(path, id, patch) Result
      +read(path) tuple
      +parse(text) tuple
      +digest(runs) Digest
      +next_id(runs) String
      +parse_quantity(s) Option~f64~
      +to_csv(runs) String
  }
  RunRecord --> RunStatus
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `RunStatus` | `runs.rs:15` | variants: `Ok`, `Failed`, `Cancelled` |
| struct | `Artifact` | `runs.rs:23` | fields: `kind`, `path`, `bytes` |
| struct | `RunRecord` | `runs.rs:32` | fields: `id`, `t`, `sec`, `status`, `spec`, `femcore`, `version`, `params`, `metrics`, `tags`, `star`, `note`, ... |
| struct | `RunPatch` | `runs.rs:58` | fields: `star`, `tags`, `note` |
| fn | `append` | `runs.rs:82` | `pub fn append(path: &Path, r: &RunRecord) -> Result<(), String>` |
| fn | `append_patch` | `runs.rs:86` | `pub fn append_patch(path: &Path, id: &str, patch: RunPatch) -> Result<(), String>` |
| fn | `read` | `runs.rs:91` | `pub fn read(path: &Path) -> (Vec<RunRecord>, usize)` —  The runs in file order with patches applied, and the number of lines that could not be read (a torn last line, say). |
| fn | `parse` | `runs.rs:96` | `pub fn parse(text: &str) -> (Vec<RunRecord>, usize)` |
| fn | `digest` | `runs.rs:122` | `pub fn digest(runs: &[RunRecord]) -> Digest` —  Counts and ranges of `runs` (failed and cancelled runs count but add no ranges). |
| fn | `next_id` | `runs.rs:154` | `pub fn next_id(runs: &[RunRecord]) -> String` —  The next run id for `runs`: `r0007_7hb` (number, then three characters of entropy so two writers cannot collide). |
| fn | `parse_quantity` | `runs.rs:162` | `pub fn parse_quantity(s: &str) -> Option<f64>` —  A number with an optional unit suffix, in SI: `450 nm`, `1.55um`, `193 THz`, `2e-7`. |
| struct | `Filter` | `runs.rs:205` |  A parsed filter such as `params.w between 450 nm and 550 nm and metrics.neff_re > 2.4 and star`. |
| fn | `to_csv` | `runs.rs:314` | `pub fn to_csv(runs: &[&RunRecord]) -> String` —  CSV of runs: id, time, status, version, then every parameter and metric column (the union, sorted). |
| fn | `Filter::parse` | `runs.rs:211` | `pub fn parse(text: &str) -> Result<Filter, String>` —  Fields: `params.NAME`, `metrics.NAME`, `version`, `sec`; `status = ok\|failed\|cancelled`, `tag = NAME`, `star` and `note` alone. |
| fn | `Filter::is_empty` | `runs.rs:275` | `pub fn is_empty(&self) -> bool` |
| fn | `Filter::matches` | `runs.rs:279` | `pub fn matches(&self, r: &RunRecord) -> bool` |
