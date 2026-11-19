# 93. Scripts are Python-syntax text run by Monty in a worker process, behind an engine path guard

Date: 2026-11-18

## Status

**Amended by ADR 94** (provisional number, productivity plan S7): the MCP server offers four script tools, `script_run`, `script_wait`, `script_check` and `api`, instead of one tool per engine command, after Paul's answer to Q9 (2026-11-19); the engine commands themselves are unchanged.

Accepted as step S4 of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3b), after Paul's choice (2026-11-18). **Number 93: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.** Amends ADR 57 (the front-ends over `Engine::call` are now three: the window, the MCP server, scripts) and uses ADR 86 (work that can take the process down runs in a child process). Python stays out in the sense of ADR 13: no CPython is embedded or needed; Monty is a Python-syntax interpreter written in Rust, and Paul accepted it on that reading.

## Context

An agent or a person who wants a workflow with loops, decisions and sweeps had two choices: one tool call per command (ADR 15, the cost S1 to S3 reduced), or a JSON command list (ADR 89), which has no loops and no decisions. The plan wanted a third thin front-end: a text file in an embedded, pure-Rust language that calls engine commands. The S4 spike bound three candidates to the whole catalog and ran the same sample in each (`docs/PRODUCTIVITY_STATUS.md`, S4 spike). The numbers were the same in all three, and the interpreter's cost (1 to 25 microseconds per call) is not a factor against solves of seconds. What differed:

* **Starlark** turns on `serde_json`'s `arbitrary_precision` for every crate of the build, which breaks nine tests (three untagged enums, one of them part of the file format); it cannot loop with `while`, has no `try`, no `%.2f`, and a Python writer's first top-level `for` is an error.
* **Rhai** is clean to build and bounds memory in-process, but its syntax is new to every reader and every model.
* **Monty** takes real Python syntax (nothing to unlearn; an agent writes it natively), has no side effect on the build, and gives CPython-style errors with the line. It needs a process boundary for memory: a script that doubles a list in a loop aborts the process the interpreter runs in.

Two gates of the plan were not met by any candidate in the spike: a script must not read or write a file outside the workspace, and a script that allocates without bound must be stopped.

## Decision

* **Language: Monty 1.6.0** (crates `monty-pool`, `monty-types`, `monty-proto`, `monty-alloc`, MIT, pinned with `=`). A script is a text file `*.nys` (the extension is ours: `*.py` is git-ignored in this repository on purpose, CLAUDE.md). Its value is its last expression; what it prints and returns is all that goes back. Engine commands are `nyx.<command>(name=value, ...)`, one method per entry of `Engine::catalog()` (nothing is listed by hand, ADR 57), plus `nyx.call(name, **args)`; an engine error is raised as `RuntimeError("<command>: <message>")`, which the script can catch.
* **The interpreter runs in a worker process, never in the engine's** (`nyx-script-worker`, a program of the `nyx-script` package: Monty's own worker mode, vendored as one 100-line file so that building nyx builds its worker and nothing is installed beside it). The engine stays in the host process: each command a script makes pauses the script, the host answers with `Engine::call`, and the script goes on (Monty's native model). A script that allocates without bound, overflows the stack or hangs ends one worker and the host goes on (`tests/script.rs`: doubling a list and a string, `[0] * 10**12`, `while True`). One worker per run, none kept.
* **Limits are strict by default and raised per run** (`Limits`): 30 s of script time (**the time the host spends answering a command is not counted**, so a 40 s solve is not script time), 256 MiB for the worker, call depth 100, 20,000 commands, 1 MiB of output. `nyx run FILE.nys --seconds N --memory MB --calls N` raises them. The pool also kills a worker that is silent past the time limit plus ten seconds.
* **No ambient access** (ADR 11, 13). The sandbox has no file, network, environment or process access; `open`, `pathlib` and `os` raise `PermissionError` with the path. **The clock is refused** (`RuntimeError`), **`random` starts from the fixed seed 0** (CPython's `random.seed(0)` sequence), `sleep` returns at once: the same script and inputs give the same numbers.
* **The path guard is in the engine** (`nyx_engine::path_guard`, `Engine::set_path_root`), not in the script host: with a root set, every argument named `path`, `folder` or `workspace` of any command is rewritten to an absolute path inside the root (a relative path is taken from the root) or the command fails before it runs; `..`, an absolute path elsewhere, a sibling folder with the same prefix and a link or junction out of the root are refused. `PATH_ARGS` names the guarded arguments and a test over the catalog fails when a command takes a path-like argument not in it (it found `store_save.workspace` while this was written). The host turns the guard on for a run, with the script's workspace (else the script's folder) as root, and restores it after.
* **A static check runs before every script and as `nyx check FILE.nys`.** Monty type-checks the script against a stub generated from the catalog (`nyx_script::stubs::nyx_pyi`, also `nyx api stubs`): a misspelled command or argument, a missing required argument, an unknown module and a syntax error fail in milliseconds with the line named and nothing ran. **Names and required-ness are checked, not value types**: the catalog says `number` for a field that also takes an expression text (`"w/2"`), and a check that refuses what the engine accepts is worse than none. The same stub, saved as `nyx.pyi` next to a script, gives an editor the commands and their documentation.
* **A run keeps the script that made it**: while a script runs the engine holds its text, and each run it records (`solve_run`, a sweep point) gets the file `script-<hash>.nys` in its folder (blake3, 32 hex digits, as the store hashes everything) with an artifact of kind `script`. No field of the run record changes, so no file format and no exemplar moves.
* **`nyx check FILE.nys` and `nyx run FILE.nys`** extend S1's commands to scripts (a file ending in `.nys`); `run` checks first.

## Consequences

* The worker is a second program to ship: `nyx-script-worker` next to `nyx` (35 MB, most of it the bundled type checker). The program finds it in its own folder or in `NYX_SCRIPT_WORKER`; without it a script fails with that sentence.
* Building the workspace now also builds Monty and its type checker (`ty`), about two minutes cold. `Cargo.lock` pins `salsa` to 0.28.2, the version Monty's own lock file pins (0.28.5, which cargo picks otherwise, does not compile with `ruff_python_ast` 0.0.14). The `monty-*` crates must move together.
* Monty is nine months old; its Python subset is large but not complete (no third-party modules, `assert` messages are pytest-style, a script that ends in an assignment returns nothing). Errors say what is missing (`ModuleNotFoundError`, a type-check message). Its sandbox is young: the path guard and the worker's isolation are ours and are tested here, and the interpreter's own isolation is relied on, not audited.
* A script's start-up costs a process: about 0.14 s measured (spawn and type check; 0.85 s the first time after a build, while the 35 MB worker is read from disk). A script is for a workflow, not for a loop of cheap commands; per-element work is a missing engine command (ADR 90).

## Not in this decision

* **A plan (dry run) of a script.** A script's commands depend on its results (a convergence loop, `if res["stale"]`), so a dry run that answers commands with placeholders cannot follow it. `nyx check FILE.nys` is the static stage only. What the script would mesh and solve can be shown from a run (the queue's estimates, ADR 86), not before it. A way to plan a script is an open design point.
* Submitting a script to the queue as one request with an owner (plan Q3), the MCP script mode (S7), a script browser in the dashboard (`scripts/` is a folder convention: the workspace's scripts live there, nothing lists them yet), N2 formulas and N3 tracing of script functions (S5, S6).
* A pool of kept workers: one worker per run costs about 0.14 s and is the simplest thing that isolates.

## Code

`crates/nyx-script/src/{lib,host,convert,stubs}.rs`, `crates/nyx-script/src/bin/nyx-script-worker.rs`, `crates/nyx-engine/src/path_guard.rs` (and `Engine::{set_path_root, set_script}`, `commands.rs` `call`, `store_link.rs`, `sweep.rs`), `crates/app/src/{runfile,apicmd}.rs`. Tests: `crates/nyx-script/tests/script.rs`, `crates/nyx-script/tests/contract.rs`, `crates/nyx-engine/tests/path_guard.rs` and the unit tests of `path_guard.rs`. Samples: `docs/scripting-samples/strip_sweep.nys`.

## Contract

`script_api.v1` (`docs/architecture/contracts/script_api.v1.json`): what a script can rely on: the extension, the object and its call forms, the default limits, the guarded argument names, the clock and random policy, the script artifact. The commands themselves are `engine_commands`.
