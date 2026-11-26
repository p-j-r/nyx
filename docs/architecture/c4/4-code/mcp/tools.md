# Script tools: code

C4 level 4 (code) of the component **Script tools** (`mcp.tools`) in the container **nyx-mcp** (crate `nyx_mcp`). [Component view of the container](../../3-components/mcp.md) · [All code pages](../README.md)

> Four tools instead of one per command: `script_run` (a script on the kernel; `params`, `path`, `save_as`, `wait_s`, `allow_destructive`), `script_wait` (continue or cancel a script that is waiting for a job), `script_check` (the static check) and `api` (the guide, the command list, one command). One answer, as text, with the project header.

**Decisions:** [ADR 94](../../../adr/0094-the-mcp-server-is-four-script-tools-over-one-kernel.md), [ADR 57](../../../adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md)  
**Size:** 1 file, 1077 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-mcp/src/lib.rs` | 1077 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `engine` (nyx-engine)](../../3-components/engine.md) |  (declared: The one engine the scripts drive, and its project header) |  |
| [container `script` (nyx-script)](../../3-components/script.md) |  (declared: api: the command list and one command's documentation) |  |

## Used by

Nothing else in the workspace names this component (entry points are reached by the program or by tests).

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class fns_lib_rs {
      <<module lib.rs>>
      +DEFAULT_WAIT_S const
      +MAX_WAIT_S const
      +MAX_VALUE_BYTES const
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `DEFAULT_WAIT_S` | `lib.rs:34` | `f64` —  How long `script_run` and `script_wait` give a script before answering `running`, and the most they may be asked to. |
| const | `MAX_WAIT_S` | `lib.rs:37` | `f64` —  The longest wait unless `--max-wait` says otherwise. Hosts end a tool call after 60 to 180 s (Antigravity: 180 s); an answer `running` is cheap, |
| const | `MAX_VALUE_BYTES` | `lib.rs:45` | `usize` —  A script's value is cut past this many bytes of JSON: the point of a script is that fields stay in the engine. |
