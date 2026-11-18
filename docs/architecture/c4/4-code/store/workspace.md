# Workspace: code

C4 level 4 (code) of the component **Workspace** (`store.workspace`) in the container **nyx-store** (crate `nyx_store`). [Component view of the container](../../3-components/store.md) · [All code pages](../README.md)

> The folder behind the dashboard: settings, the list of workflows (manifests only) and the sections of the dashboard.

**Decisions:** [ADR 14](../../../adr/0014-workflow-store-json-manifests-and-nyx-archives.md)  
**Size:** 2 files, 306 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-store/src/lib.rs` | 22 | nyx-store: where nyx keeps its work. |
| `crates/nyx-store/src/workspace.rs` | 284 | The workspace: the folder behind the dashboard, its settings (`workspace.json`), the list of workflows (manifests only) and the dashboard retention rules. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Workflow folder](./flow.md) | flow | WorkflowDir |
| [Data model and card](./model.md) | model |  |
| [Runs](./runs.md) | runs | RunRecord, RunStatus |
| [Utilities](./util.md) | util |  |

## Used by

| Component | Modules | Items |
|---|---|---|
| [engine: Registry and facade](../engine/registry.md) | crate root |  |
| [Search](./search.md) | workspace | Entry |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace workspace_rs {
    class Settings {
        <<struct>>
        +hide_after_days u32
        +dashboard_cap usize
        +keep_heavy Option~usize~
        +failed_heavy_days Option~u32~
        +share_default ShareLevel
        +strip_paths bool
        +run_shared_scripts String
        +run_at_once usize
        +more_11_fields
    }
    class Entry {
        <<struct>>
        +dir PathBuf
        +wf Workflow
    }
    class Sections {
        <<struct>>
        +running Vec~usize~
        +recent Vec~usize~
        +earlier Vec~usize~
    }
    class Workspace {
        <<struct>>
        +root PathBuf
        +settings Settings
        +open(root) Result
        +save_settings() Result
        +list() Vec~Entry~
        +find(id) Option~Entry~
        +sections(entries, running, now) Sections
        +housekeeping(w, now) Result
    }
    class fns_workspace_rs {
        <<module workspace.rs>>
        +default_root() PathBuf
    }
  }
  Workspace --> Settings
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `Settings` | `workspace.rs:17` | fields: `hide_after_days`, `dashboard_cap`, `keep_heavy`, `failed_heavy_days`, `share_default`, `strip_paths`, `run_shared_scripts`, `run_at_once`, `length_unit`, `frequency_unit`, `theme`, `assistant`, ... |
| struct | `Entry` | `workspace.rs:80` | fields: `dir`, `wf` —  One workflow on disk as the dashboard sees it. |
| struct | `Sections` | `workspace.rs:86` | fields: `running`, `recent`, `earlier` |
| struct | `Workspace` | `workspace.rs:93` | fields: `root`, `settings` |
| fn | `default_root` | `workspace.rs:99` | `pub fn default_root() -> PathBuf` —  `Documents\nyx` of the user (the current folder when there is no home). |
| fn | `Workspace::open` | `workspace.rs:105` | `pub fn open(root: &Path) -> Result<Workspace, String>` —  Opens (and creates) the workspace folder; a missing or unreadable `workspace.json` gives the default settings. |
| fn | `Workspace::save_settings` | `workspace.rs:111` | `pub fn save_settings(&self) -> Result<(), String>` |
| fn | `Workspace::list` | `workspace.rs:117` | `pub fn list(&self) -> Vec<Entry>` —  Every workflow folder with a readable manifest, newest id first. Only `workflow.json` is read. |
| fn | `Workspace::find` | `workspace.rs:132` | `pub fn find(&self, id: &str) -> Option<Entry>` |
| fn | `Workspace::sections` | `workspace.rs:141` | `pub fn sections(&self, entries: &[Entry], running: &[String], now: i64) -> Sections` —  Dashboard sections for `entries` at time `now`: |
| fn | `Workspace::housekeeping` | `workspace.rs:167` | `pub fn housekeeping(&self, w: &mut WorkflowDir, now: i64) -> Result<u64, String>` —  Removes heavy artifacts (everything under `runs/<id>/` except the log) according to the settings. Returns the bytes freed. |

Trait implementations in this component: `Settings: Default`
