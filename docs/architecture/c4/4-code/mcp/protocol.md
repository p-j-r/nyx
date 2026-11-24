# Protocol loop: code

C4 level 4 (code) of the component **Protocol loop** (`mcp.protocol`) in the container **nyx-mcp** (crate `nyx_mcp`). [Component view of the container](../../3-components/mcp.md) · [All code pages](../README.md)

> Newline-delimited JSON-RPC 2.0 over stdio, no async runtime and no SDK; initialize, ping, tools/list and tools/call. `--dir` (workflows and saved scripts), `--root` (the folder a script's file arguments may name); `--job` is the engine's solve worker.

**Decisions:** [ADR 15](../../../adr/0015-ai-interface-separate-workspace-and-mcp-server.md), [ADR 57](../../../adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md), [ADR 94](../../../adr/0094-the-mcp-server-is-four-script-tools-over-one-kernel.md)  
**Size:** 1 file, 813 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-mcp/src/main.rs` | 813 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `engine` (nyx-engine)](../../3-components/engine.md) |  (declared: Server instructions text) |  |

## Used by

Nothing else in the workspace names this component (entry points are reached by the program or by tests).

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Server {
      <<struct>>
      -new(dir, root, queue, max_wait, rx) Server
      -initialize(params) Value
      -handle(msg) Option~Value~
  }
  class fns_main_rs {
      <<module main.rs>>
      -home() Option~PathBuf~
      -default_dir() PathBuf
      -default_root(dir_given, dir) PathBuf
      -session_id() String
      -serve(dir, root, queue, max_wait)
      -usage()
      -take_opt(args, flag) Option~PathBuf~
      -main()
      -NAME const
      -VERSION const
      -PROTOCOLS const
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `NAME` | `main.rs:29` | `&str` |
| const | `VERSION` | `main.rs:30` | `&str` |
| const | `PROTOCOLS` | `main.rs:32` | `&[&str]` —  Protocol revisions this server speaks, newest first. It answers `initialize` with the client's if listed, else the newest. |
| fn | `home` | `main.rs:47` | `fn home() -> Option<PathBuf>` |
| fn | `default_dir` | `main.rs:51` | `fn default_dir() -> PathBuf` |
| fn | `default_root` | `main.rs:57` | `fn default_root(dir_given: bool, dir: &Path) -> PathBuf` —  The folder the file arguments of a script's commands must lie in (the engine's path guard, ADR 93): `--root`, else the workspace when |
| fn | `session_id` | `main.rs:109` | `fn session_id() -> String` —  A short id for this process: the session a result belongs to. |
| struct | `Server` | `main.rs:266` |  |
| fn | `serve` | `main.rs:735` | `fn serve(dir: PathBuf, root: PathBuf, queue: PathBuf, max_wait: f64)` |
| fn | `usage` | `main.rs:757` | `fn usage()` |
| fn | `take_opt` | `main.rs:765` | `fn take_opt(args: &mut Vec<String>, flag: &str) -> Option<PathBuf>` —  Takes `--flag VALUE` out of the arguments. |
| fn | `main` | `main.rs:776` | `fn main()` |
| fn | `Server::new` | `main.rs:288` | `fn new(dir: PathBuf, root: PathBuf, queue: PathBuf, max_wait: f64, rx: mpsc::Receiver<Incoming>) -> Server` |
| fn | `Server::initialize` | `main.rs:394` | `fn initialize(&mut self, params: &Value) -> Value` |
| fn | `Server::handle` | `main.rs:699` | `fn handle(&mut self, msg: &Value) -> Option<Value>` —  One incoming message; `None` means no reply (notifications). |
