# Script host: code

C4 level 4 (code) of the component **Script host** (`script.host`) in the container **nyx-script** (crate `nyx_script`). [Component view of the container](../../3-components/script.md) · [All code pages](../README.md)

> Runs a script (Python syntax, `*.nys`) through a worker process and answers each `nyx.<command>(...)` it makes with Engine::call; the limits, the clock and random policy, the path guard turned on for the run, the script kept with the runs it makes; `check` parses and type-checks without running. JSON to and from script values.

**Decisions:** [ADR 93](../../../adr/0093-scripts-are-python-text-run-by-monty-in-a-worker-process.md), [ADR 57](../../../adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md)  
**Size:** 3 files, 217 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-script/src/convert.rs` | 65 |  |
| `crates/nyx-script/src/host.rs` | 41 |  |
| `crates/nyx-script/src/lib.rs` | 111 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `engine` (nyx-engine)](../../3-components/engine.md) | engine | Engine |
| [Script kernel](./kernel.md) | kernel | Kernel, Step |
| [Worker program](./worker.md) |  (declared: Frames over a pipe: the script pauses at each command and the host answers) |  |

## Used by

| Component | Modules | Items |
|---|---|---|
| [app: Scripts and batch](../app/headless.md) | crate root, host | Limits, Options, check, run |
| [app: Window](../app/window.md) | crate root, host | Limits, Options, check, find_worker, mcp_dir, run |
| [Script kernel](./kernel.md) | convert, crate root | Limits, Options, Outcome, args_json, find_worker, from_object, to_object |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace convert_rs {
    class fns_convert_rs {
        <<module convert.rs>>
        +to_object(v) MontyObject
        +from_object(o) Result
        +args_json(args, skip) Result
    }
  }
  namespace host_rs {
    class fns_host_rs {
        <<module host.rs>>
        +run(engine, name, src, limits, opts) Result
        +check(engine, name, src, limits, opts) Result
    }
  }
  namespace lib_rs {
    class Limits {
        <<struct>>
        +seconds f64
        +memory_bytes usize
        +depth usize
        +calls usize
        +output_bytes usize
    }
    class Options {
        <<struct>>
        +root PathBuf
        +worker Option~PathBuf~
        +new(root) Options
    }
    class Outcome {
        <<struct>>
        +value Value
        +printed Vec~String~
        +calls u64
    }
    class fns_lib_rs {
        <<module lib.rs>>
        +find_worker() Option~PathBuf~
        +mcp_dir() PathBuf
    }
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `to_object` | `convert.rs:10` | `pub fn to_object(v: &Value) -> MontyObject` —  JSON to a script value. |
| fn | `from_object` | `convert.rs:25` | `pub fn from_object(o: ObjectRef) -> Result<Value, String>` —  A script value to JSON; a value JSON cannot hold (a function, a set, a dict key that is not text, an integer past 64 bits) is an error. |
| fn | `args_json` | `convert.rs:50` | `pub fn args_json(args: &CallArgs, skip: usize) -> Result<Value, String>` —  The arguments of one command call as the engine's JSON object: the keyword arguments, or one positional dict. `skip` leading positional |
| fn | `run` | `host.rs:23` | `pub fn run(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<Outcome, String>` —  Runs a script. Its value is its last expression. The engine's path guard is on for the run and restored after it. A person runs their own |
| fn | `check` | `host.rs:35` | `pub fn check(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<(), String>` —  Parses and type-checks a script against the commands of `engine` and runs nothing: a misspelled command or argument, a missing argument, a |
| struct | `Limits` | `lib.rs:33` | fields: `seconds`, `memory_bytes`, `depth`, `calls`, `output_bytes` —  What the host allows a script to use. The defaults are strict; a run may raise them (S4 decision). |
| struct | `Options` | `lib.rs:54` | fields: `root`, `worker` —  How a script is run. |
| struct | `Outcome` | `lib.rs:69` | fields: `value`, `printed`, `calls` —  What a finished script gave back. |
| fn | `find_worker` | `lib.rs:81` | `pub fn find_worker() -> Option<PathBuf>` —  The script worker program: named by `NYX_SCRIPT_WORKER`, else `nyx-script-worker` next to the running program |
| fn | `mcp_dir` | `lib.rs:104` | `pub fn mcp_dir() -> PathBuf` —  The folder the MCP server keeps its workspace in by default: `Documents/nyx_runs/mcp` of the person (the temporary folder when there is no |
| fn | `Options::new` | `lib.rs:62` | `pub fn new(root: impl Into<PathBuf>) -> Options` |

Trait implementations in this component: `Limits: Default`
