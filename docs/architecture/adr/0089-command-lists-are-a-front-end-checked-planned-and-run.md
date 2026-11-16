# 89. Command lists are a front-end: checked, planned and run

Date: 2026-11-16

## Status

Accepted as a step of the productivity plan (S1, `docs/NYX2_PRODUCTIVITY_PLAN.md` section 3b). **The number 89 is the one the plan reserves and is provisional**: take the next free ADR number at merge and renumber this file, the index row and the code comments. Amends ADR 57 (front-ends are thin over `Engine::call`) with a third front-end; the scripting language of step S4 will extend this record, not replace it. **Built on** the queue (ADR 87) and on P0's command lists and the power-throttling opt-out (ADR 11 note, `process.rs`): the S1 bookmark is a merge of those lines.

## Context

An agent that wants a workflow (open, change a width, solve, read the index, five times) makes one tool call per command and finds a wrong argument at the call that uses it, perhaps an hour into a run. A person who wants the same workflow again has nothing to re-run. P0 wrote the workflows as JSON lines of `{cmd, args}` for timing; the same files are the smallest useful script.

## Decision

* **A command list is a file of JSON lines, `{"cmd": <command of Engine::catalog()>, "args": {...}}`**, sent to `Engine::call` in order. It is not a second model format and not a language: no loops, no values carried from one line to the next (the replay keys `job_arg`, `after_stage`, `delay_s`, `expect_error` and `note` exist for the cancel workflows and the timing tables). Blank lines and `#` lines are skipped. A language comes with step S4.
* **Three stages, the first two before anything runs** (`nyx_engine::commandlist`): *parse* (every line is an object with a `cmd` and known keys); *static* (the command is in the catalog, with a did-you-mean when it is not, and `args` fit its JSON Schema: type, enum, `oneOf`, required, no unknown argument; a keyword the checker does not know is not checked, so a valid command is never refused); *plan*, a dry run on a scratch copy of the engine in which every command that changes the project is applied for real and the commands that start work are only counted. The plan finds what only the project can tell (a solve with no physics at that point, a feature that does not exist, a problem that does not check, a result read before anything solves it) and says how many meshes and solve jobs the list would make, with the meshes the stage cache would reuse (ADR 87). It meshes nothing, solves nothing, writes nothing and starts no job. A bad line is reported with its line number and command and nothing runs; every problem of the static stage is reported at once.
* **Run** sends the lines to an engine in order, through the same commands the window and the MCP server use, so a solve waits in the queue like any other. It stops at the first failing command and at a replay key that cannot be honoured (a cancel that found no running job, an `expect_error` line that succeeded), so a step that did not happen cannot pass.
* **Front ends.** `nyx check FILE` prints the plan (`--json` for an agent), exit 2 on a problem; `nyx run FILE` checks and then runs (`--no-check` to skip the check, `--workspace DIR` to record every solve as a run of a workflow there, `--json` for one JSON line per command); exit 1 when a command fails. The `workflow` example replays through the same `run`.
* **The reference is generated from the same catalog** (`nyx_engine::apiref`, S2): `docs/api/` holds one page per command (title, access, description, a table of the arguments flattened from the JSON Schema with their types, required flags, defaults and the unit or meaning the schema states, and an example line made from the schema that passes `nyx check`) and an index; `nyx api list` and `nyx api show COMMAND` print the same text. A test regenerates the pages and fails when `docs/api/` differs, so the reference cannot disagree with a command, and another requires every argument to have a description. Rewrite after a deliberate change with `NYX_WRITE_API=1 cargo test --release -p nyx-engine --test api_docs`.
* **Where the logic lives.** In `nyx-engine` (`commandlist.rs`), with the window and MCP server free to use it later; the app crate only parses arguments and prints (`runfile.rs`).

## Not in this decision

* **A script language, values between lines, loops, the run record holding the script text** (steps S4 to S6). A list says what to do, not what to do next.
* **Cost estimates in the plan** (plan Q2): the plan counts meshes and solves and does not say how long they take.
* **Planning commands that read results** (`result_info`, `analyse`, ...): they are checked for a result to read, not run.
* **`open_results` in a plan**: it replaces the project with a file's, so the lines after it are checked without the checks that need to know what was solved.

## Consequences

* The P0 workflows are the first lists and are checked by a test (every list checks clean and the check starts nothing).
* `Engine::sweep_prepare` is the check of a sweep without submitting it, shared by `sweep_run` and the plan.
* No contract changes: the commands are the catalog's.

## Code

`crates/nyx-engine/src/{commandlist,apiref,sweep}.rs`, `crates/app/src/{runfile,apicmd,main}.rs`, `docs/api/`, `crates/nyx-engine/examples/workflow.rs`; tests `crates/nyx-engine/tests/{commandlist,api_docs}.rs`, `crates/app/tests/cli.rs`.
