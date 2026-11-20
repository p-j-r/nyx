# Scripts and batch: code

C4 level 4 (code) of the component **Scripts and batch** (`app.headless`) in the container **app** (crate `nyx`). [Component view of the container](../../3-components/app.md) · [All code pages](../README.md)

> Runs without a person at the window: scripted runs with screenshots (--script), batch runs for clusters, command lists (nyx check, nyx run) and the command reference (nyx api), and the optional user Python script.

**Decisions:** [ADR 23](../../../adr/0023-app-is-a-thin-shell-over-one-project.md), [ADR 89](../../../adr/0089-command-lists-are-a-front-end-checked-planned-and-run.md), [ADR 93](../../../adr/0093-scripts-are-python-text-run-by-monty-in-a-worker-process.md)  
**Size:** 5 files, 733 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/app/src/apicmd.rs` | 71 |  |
| `crates/app/src/batch.rs` | 58 |  |
| `crates/app/src/pyscript.rs` | 55 |  |
| `crates/app/src/runfile.rs` | 231 |  |
| `crates/app/src/script.rs` | 318 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `engine` (nyx-engine)](../../3-components/engine.md) | apiref, commandlist, crate root, engine | Done, Engine, Problem, Scalar |
| [container `script` (nyx-script)](../../3-components/script.md) | crate root, host, stubs | Limits, Options, check, nyx_pyi, run |
| [Shell and dashboard](./shell.md) | examples, shell | Page, Part, Sheet, Start, Sub, list |
| [Look and widgets](./style.md) | cvd | Kind, simulate_rgba |
| [Window](./window.md) | app, state | App, Msg |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Window](./window.md) | pyscript | spawn |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace apicmd_rs {
    class fns_apicmd_rs {
        <<module apicmd.rs>>
        +run(args) Result
        +Exit type
    }
  }
  namespace batch_rs {
    class fns_batch_rs {
        <<module batch.rs>>
        +run(args) Result
    }
  }
  namespace pyscript_rs {
    class fns_pyscript_rs {
        <<module pyscript.rs>>
        +spawn(script, project_json, csv, tx, ctx)
    }
  }
  namespace runfile_rs {
    class fns_runfile_rs {
        <<module runfile.rs>>
        +check(args) Result
        +run(args) Result
        +Exit type
    }
  }
  namespace script_rs {
    class App {
        <<impl block>>
        +script_step(ctx)
    }
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| type | `Exit` | `apicmd.rs:12` |  |
| fn | `run` | `apicmd.rs:18` | `pub fn run(args: &[String]) -> Result<(), Exit>` |
| fn | `run` | `batch.rs:16` | `pub fn run(args: &[String]) -> Result<(), String>` |
| fn | `spawn` | `pyscript.rs:11` | `pub fn spawn(script: PathBuf, project_json: String, csv: Option<String>, tx: Sender<(u64, Msg)>, ctx: egui::Context)` |
| type | `Exit` | `runfile.rs:25` |  |
| fn | `check` | `runfile.rs:148` | `pub fn check(args: &[String]) -> Result<(), Exit>` —  `nyx check FILE`. |
| fn | `run` | `runfile.rs:179` | `pub fn run(args: &[String]) -> Result<(), Exit>` —  `nyx run FILE`. |
| fn | `App::script_step` | `script.rs:18` | `pub fn script_step(&mut self, ctx: &egui::Context)` |
