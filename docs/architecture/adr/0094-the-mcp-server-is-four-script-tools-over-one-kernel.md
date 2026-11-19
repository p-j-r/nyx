# 94. The MCP server is four script tools over one long-lived kernel

Date: 2026-11-19

## Status

Accepted as step S7 of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3b), after Paul's answer to Q9 on 2026-11-19: **script tools only; the one-tool-per-command catalog is deleted, not kept behind a flag.** **Number 94: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.** Amends ADR 15 (the MCP tool list), ADR 57 (the MCP server's tools are no longer the catalog) and ADR 93 (the script host becomes a kernel that outlives one script).

## Context

ADR 57 made every engine command an MCP tool: 99 tools, about 50 KB of definitions, one call per command, and every result sent twice (as text and as `structuredContent`). S1 to S3 cut the calls an agent needs, S4 gave scripts, S5 and S6 gave them formulas. An agent that has to make five calls to change one width, resend a whole model description to do it, and read the full project JSON to find one feature is paying for the interface and not for the physics. Anthropic's "Code execution with MCP" and Cloudflare's "Code Mode" make the same point: the agent writes code against a typed API and only what the code prints comes back.

Monty (ADR 93) is more than "a script runs". Measured on 1.6.0 with a stand-in host (2026-11-19):

* **A paused script is a coroutine the host steps.** It stays paused at a host call for as long as the host likes: a 20 s gap passed with a 5 s feed limit and a 15 s request timeout. Host time is not script time.
* **State persists between feeds.** Variables and functions of one feed are known to the next, and the static check knows them too: `x = 41` then `x + 1` passes, a misspelled name or argument in a later feed is still refused with its line. An input (`params`) declared in the stub type-checks.
* **A session can be dumped and restored**: 1.4 KB idle, 2.5 KB paused at a host call, restored in a new worker with variables, functions and the host module binding intact. The dump holds the interpreter only, not the engine.
* **Futures work.** A host answers a call to an async function with `ResumeValue::Future` and settles it on `TurnEvent::ResolveFutures`; top-level `await`, `asyncio.gather` (also over a list comprehension, and from an `async def` of the script), start now and await later, and an exception raised at the `await` all work. `gather(return_exceptions=True)` and `asyncio.create_task` do not exist in 1.6.0. Three 0.5 s timer jobs gathered took 0.51 s against 1.53 s awaited one by one. **That gain is not carried over to solves**: the engine's heavy lane runs one solve at a time (ADR 87), so futures give pipelining (post-process each point as it finishes) and overlap with light work, not faster heavy solves.
* Monty has no screenshots.

## Decision

* **Four tools, and no others.**
  * `script_run {source | path, params, wait_s, allow_destructive, save_as}`: runs a script and answers `done`, `running` or `error`.
  * `script_wait {wait_s, cancel}`: continues a script that answered `running`, or abandons it.
  * `script_check {source | path}`: the static check of ADR 93 (parse, type-check against the stub, the traced-function guard of ADR 90 N3), nothing runs.
  * `api {topic | command}`: the API text. Without arguments: the instructions and the list of commands with their titles. A guide card (`geometry`, `optical-modes`, `troubleshooting`) or one command with its full description and arguments.
  The 99 per-command tools, `get_guide` as a tool, `nyx-mcp --list` and `--call`, the per-tool annotations and the per-result header are removed. The engine commands stay: the window and the scripts use them (`engine_commands` is unchanged). `--job` stays: it is the engine's solve worker (ADR 86).
* **One kernel.** The server creates a kernel on the first script call (so `initialize` and `tools/list` need no worker): one runtime, one worker pool, one Monty session that lives until the server stops. Variables, functions and imports of one `script_run` are there in the next. One script runs at a time; `script_run` while another is `running` is an error that names it.
* **A call never blocks for ever.** `solve_run`, `convergence_check`, `job_wait` and `queue_wait` are served by the host as: start the job without waiting, then poll `job_wait` with `timeout_s` 0 until it finished, the call's own timeout passed (`job_wait`, `queue_wait`), or the time given by `wait_s` (default 25 s) is spent. In the last case the script stays paused at the call and `script_run` answers `running` with what the script printed so far; `script_wait` continues the same poll. The value the script gets is the value the blocking call would have given. A job that cannot start (the engine says why) is an error, not an endless `running`. The host's own polling is not a script command: it does not count against `Limits.calls` and is not in the run record.
* **`params`.** A dict given to `script_run` is the script's input named `params` (declared in the stub, so it type-checks); a hand-written script reads `params.get("w", 500e-9)`. A saved script is run again with other values and nothing is resent.
* **Saved scripts** are `<workspace>/scripts/<name>.nys`. `path` is the name (letters, digits, `_ - .`, no folder); `save_as` stores the script just run. Nothing else in the file system is reachable (ADR 93's guard is unchanged).
* **Destructive commands need a flag.** `project_new`, `project_load`, `project_open`, `store_open` and every `remove_*` raise `RuntimeError` in a script unless `script_run` was called with `allow_destructive: true`, which the host's approval dialog shows with the call. The list lives in `nyx-script` (the MCP server used to carry it). `nyx run` allows them: a person runs their own file.
* **Results are text, once.** One JSON object: `status`, `value`, `printed` (capped at `Limits.output_bytes`), `calls`, the script's hash and the project header `nyx {session, revision, project_hash}` once per answer instead of once per command. No `structuredContent`.
* **A crashed worker costs the call, not the session.** After every finished script the host keeps an in-memory dump; when the worker dies (memory limit, stack) the next call restores the last dump into a new worker, so variables survive. Dumps are not written to disk: a restart of the server starts a new kernel (the interpreter state is not the engine's, and a dump is not authenticated).
* **`script_run`, `nyx run` and `nyx check` share one driver** (`nyx_script::kernel`); `host::run` and `host::check` are thin wrappers over it.
* **Futures (`nyx.aio`).** Decided after the poll loop above exists: the async twins `await nyx.aio.solve_run(...)`, `nyx.aio.convergence_check`, `nyx.aio.job_wait` reuse it. A future left pending when a feed ends, aborts or fails has its engine job cancelled.

## Consequences

* An agent works in one to three calls where it made five to forty, and fields and tables stay in the engine; only what the script prints comes back.
* The `mcp_tools` contract is a major version: four tools. The 50 KB `mcp_tools.v2.json` snapshot is replaced; the full catalog remains pinned by `engine_commands`, which is what the window and scripts use. `script_api` gains `params` and the kernel semantics (additive).
* `nyx-mcp` depends on `nyx-script` (layering table), and `nyx-script-worker` must sit next to `nyx-mcp` (`NYX_SCRIPT_WORKER`). `cargo test -p nyx-mcp` alone does not build the worker; the tests that run a script say so.
* **Known limits.** The traced-function guard (ADR 90 N3) reads the source of one feed: a function defined in an earlier `script_run` and passed as a formula in a later one is not seen by it. The dump does not hold the engine: a restarted server has to reopen the project. The window and the server are separate processes with their own engines (ADR 57): a project the agent saves is reopened in the window, there is no live link. `gather(return_exceptions=True)` and `create_task` are missing in Monty 1.6.0.
* Not built: scripts saved by the window ("export this model as a script"), a plan (dry run) mode, dumps on disk.

## Export as script (S7b)

A person who builds a model in the window, with no script, can hand it to an agent as a script: **File > Export as script...** writes `.nys` text that rebuilds the project, and the engine command `project_script` (read-only, so `nyx.project_script()` in a script and an MCP agent have it too) returns the same text. The agent then runs it with `script_run` and changes a value with `params`, and the person's model is the starting point of a script and not a JSON blob.

* **The commands, then the text.** The engine builds a list of `(command, arguments)` and renders it as Python; it does not format text first. The list is, in this order: `project_new(dim, title)`; `set_param` per parameter in stored order, written `value=params.get("w", 5e-07)`, so `script_run path=... params={"w": 3e-7}` is the variant with nothing edited; `add_<kind>(name=..., <args>)` per feature (`geometry::describe` gives the arguments in the form the add command takes), then `suppress_feature` and `tag_feature`; `set_mesh` with the fields that differ from the default; `add_material` per project material; `set_problem` per physics with the problem as a dict; `set_active`.
* **Verified before it is offered.** The list is replayed on a scratch engine and the result must equal the project. Equal means the project JSON with every feature's `origin` set to `Code` on both sides (the window's features are `Gui`, `add_*` makes `Code`; origin changes no result), and, when no overlay has pins or pinned regions, with the feature ids renumbered (a window model whose features were removed has id gaps; replay renumbers them). Three forms, first that verifies: `commands` (the list above); `geometry_json` (the same, with the geometry as one `set_geometry`: exact ids, for models whose pins refer to them); `project_json` (one `project_load` of the whole project: always exact, not readable). The answer says which, and so does the header comment of the script. The tests require `commands` for every example project, so a fallback cannot hide a bug.
* **Rendering.** `true`, `false` and `null` become `True`, `False` and `None`; numbers keep the shortest text that reads back to the same value; strings are written as JSON strings (valid Python); a dict or list longer than a line is laid out over lines. The text holds no file path, so it does not depend on the path guard's root, and it starts with `project_new`: over MCP it needs `script_run allow_destructive: true`, which the header comment says.
* **The window** pushes its buffered edits first (as Save does), asks the engine, and writes the file (default folder the MCP server's `scripts` folder, `Documents/nyx_runs/mcp/scripts`, so an agent runs it with `script_run path=<name>`; `nyx_script::mcp_dir` is the one place that path is made). Contract: `engine_commands` 1.19 (additive: `project_script`).
* **Not exported:** results, meshes and the run history (they are outputs), the window's view settings and the title bar state. A script never solves; the person or the agent adds `nyx.solve_run()`.

## Code

`nyx-script::{kernel,host,stubs,prelude}`, `nyx-mcp::main`, `nyx-engine::{guide,cards,scriptout}`, `app::app` (File menu); contracts `mcp_tools` 3.0, `script_api` 1.2, `engine_commands` 1.19.
