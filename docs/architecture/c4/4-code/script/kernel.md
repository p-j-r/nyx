# Script kernel: code

C4 level 4 (code) of the component **Script kernel** (`script.kernel`) in the container **nyx-script** (crate `nyx_script`). [Component view of the container](../../3-components/script.md) · [All code pages](../README.md)

> One Monty session in a worker process that outlives a script: `run`, `wait`, `cancel` and `check`; variables stay from one script to the next; a command that waits for a job (`solve_run`, `convergence_check`, `job_wait`, `queue_wait`) is started and polled by the host so no call blocks for ever and a script can be left `running`; `nyx.aio` futures; the destructive commands refused unless allowed; a dead worker replaced from a dump.

**Decisions:** [ADR 94](../../../adr/0094-the-mcp-server-is-four-script-tools-over-one-kernel.md), [ADR 93](../../../adr/0093-scripts-are-python-text-run-by-monty-in-a-worker-process.md)  
**Size:** 1 file, 660 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-script/src/kernel.rs` | 660 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `engine` (nyx-engine)](../../3-components/engine.md) | engine | Engine |
| [Script host](./host.md) | convert, crate root | Limits, Options, Outcome, args_json, find_worker, from_object, to_object |
| [Type stub](./stubs.md) | stubs |  |
| [Tracing of script functions](./tracing.md) | trace |  |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Script host](./host.md) | kernel | Kernel, Step |
| [Type stub](./stubs.md) | kernel | POLLABLE |
| [Tracing of script functions](./tracing.md) | kernel | POLLABLE |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Progress {
      <<struct>>
      +printed Vec~String~
      +calls u64
      +waiting String
  }
  class Step {
      <<enumeration>>
      Done
      Running
  }
  class Kernel {
      <<struct>>
      +new(engine, limits, opts) Result
      +set_script_name(name)
      +waiting() Option~String~
      +waiting_job() Option~u64~
      +run(engine, name, src, params, wait, allow_destructive) Result
      +wait(engine, wait) Result
      +cancel(engine) Result
      +check(engine, name, src) Result
  }
  class fns_kernel_rs {
      <<module kernel.rs>>
      +is_destructive(name) bool
      ~explain(e, limits) String
      +POLLABLE const
  }
  Step --> Progress
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `POLLABLE` | `kernel.rs:34` | `&[&str]` —  The commands that wait for a job. Each is served by starting the job and polling it, and each has an `nyx.aio` twin. |
| fn | `is_destructive` | `kernel.rs:40` | `pub fn is_destructive(name: &str) -> bool` —  Commands that discard the project or a part of it. |
| struct | `Progress` | `kernel.rs:46` | fields: `printed`, `calls`, `waiting` —  What a script that is not finished is waiting for. |
| enum | `Step` | `kernel.rs:56` | variants: `Done`, `Running` —  The answer to `run` and `wait`. |
| fn | `explain` | `kernel.rs:157` | `pub(crate) fn explain(e: PoolError, limits: &Limits) -> String` —  What the pool's error means to the person who wrote the script. |
| struct | `Kernel` | `kernel.rs:239` |  A Monty session that outlives a script. |
| fn | `Kernel::new` | `kernel.rs:272` | `pub fn new(engine: &Engine, limits: Limits, opts: Options) -> Result<Kernel, String>` —  A kernel for `engine`'s catalog. No process is started until the first script. |
| fn | `Kernel::set_script_name` | `kernel.rs:305` | `pub fn set_script_name(&mut self, name: &str)` —  The name the static check gives the script in its messages (`file.nys:2:1`). It is the name of the session, so it is set before the first |
| fn | `Kernel::waiting` | `kernel.rs:310` | `pub fn waiting(&self) -> Option<String>` —  What the script that has not finished is waiting for. |
| fn | `Kernel::waiting_job` | `kernel.rs:318` | `pub fn waiting_job(&self) -> Option<u64>` —  The engine job the script that has not finished is waiting for (none when it waits for futures): a host asks the engine about it. |
| fn | `Kernel::run` | `kernel.rs:333` | `pub fn run(&mut self, engine: &mut Engine, name: &str, src: &str, params: Option<&Value>, wait: Option<Duration>, allow_destructive: bool) -> Result<Step, String>` —  Feeds `src` to the session. It returns when the script is finished, or when `wait` is spent and it is waiting for a job. |
| fn | `Kernel::wait` | `kernel.rs:356` | `pub fn wait(&mut self, engine: &mut Engine, wait: Option<Duration>) -> Result<Step, String>` —  Continues the script that is waiting; at the end of `wait` it may still be. |
| fn | `Kernel::cancel` | `kernel.rs:372` | `pub fn cancel(&mut self, engine: &mut Engine) -> Result<(), String>` —  Abandons the script that is waiting: its job is cancelled and the session starts again from the last script that finished. |
| fn | `Kernel::check` | `kernel.rs:389` | `pub fn check(&mut self, engine: &mut Engine, name: &str, src: &str) -> Result<(), String>` —  Parses and type-checks `src` against the commands and the variables of the session and runs nothing that matters: the session it is |

Trait implementations in this component: `Kernel: Drop`
