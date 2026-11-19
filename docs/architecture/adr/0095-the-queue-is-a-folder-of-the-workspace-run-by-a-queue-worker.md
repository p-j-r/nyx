# 95. The queue is a folder of the workspace, run by a queue worker that outlives the window

Date: 2026-11-19

## Status

Accepted as step Q3 of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3a), with Paul's "yes" of 2026-11-19. **The number 95 is provisional** (renumber when the lines merge). Builds on ADR 86 (solves in worker processes, memory admission) and ADR 87 (lanes, priorities, owners); amends ADR 57 (the front-ends share one queue, not only the commands) and ADR 86 (the queue worker is a second kind of worker).

## Context

ADR 87 gave each `Engine` a queue; ADR 86 put every solve in a child process. The window and the MCP server are still two programs with an engine each, so a person's Compute and an agent's solve run side by side on one machine, with their memory added up, and closing the window ends what it was waiting for. Plan Q3 asks for one queue that every front-end reads and steers, that survives a front-end, and that survives the death of a worker, with no daemon to install and no network.

## Decision

* **The queue is the folder `queue/` of a workspace.** An item is a folder `queue/items/<id>/` that a front-end writes whole (a temporary folder renamed into place): `request.json` (owner, priority, label, targets, the hash of the frozen project, the memory estimate the front-end made), `project.nyx.json` and `job.json` (exactly what a worker of ADR 86 is given). The id sorts by submission time. A front-end writes nothing else but **markers**: an empty file `cancel`, an empty file `hold` (removed to release), a file `priority`. **No file has two writers**: the queue worker alone writes `state.json` (an item with none is queued), `order.json` (the ids in the order they will start, which is what a front-end shows as the place in the queue), the results (`results.h5`, `stderr.txt`) and the lease.
* **The queue worker is `PROGRAM --worker WORKSPACE`**, a process of its own (`nyx`, `nyx-mcp` and `nyx-worker` all answer it). It takes the lease, ranks the items with the job runner's rules (priority, lowered one level for each minute waited, then the owners in turns, then submission order), runs the head in a solver child exactly as ADR 86 does (`worker::run_child`, the one place a worker is started; `tests/worker.rs`'s bit comparison covers it), admits it by memory (ADR 86's `admit`, here and nowhere else, so nothing is admitted twice), appends the run to `solve_history.jsonl` (so the estimates are learned whoever was watching), and ends after 30 s with nothing to do (`NYX_QUEUE_IDLE_SECS`). It does **not** end with the program that started it, unlike a job's worker: it is started detached (no console, no standard input, outside the job object of the starter), so closing the window or the terminal leaves it running. The solver children it starts keep ADR 86's rule: they end with it.
* **The lease** is `queue/runner.lock`, taken with `create_new` and kept alive by touching its modification time every second; one not touched for six seconds is dead. A process id alone is not trusted (Windows reuses them). Two programs that find a dead lease both try; one `create_new` wins and each reads the file back to see whose token is in it. A queue worker that finds the previous one's item `running` marks it `failed: the queue worker that was running this item stopped` (its solver child went with it).
* **A front-end that waits is a job in the new `shared` lane.** `solve_run`, `queue_submit` and `sweep_run` take `shared` (the window and the MCP server set it as the default with `Engine::set_shared_queue(folder, true)`; a solve that cannot be isolated runs in the process as before). The job meshes (one mesh at a time across such jobs, so fifty sweep points do not start fifty meshers), writes the item, makes sure a queue worker runs (and starts another if the lease goes stale while it waits), shows the item's place and then the queue worker's stage as its own progress, turns a cancel of the job into the `cancel` marker, reads the results file as it does a worker's, and removes the item. The lane has no limit that matters, so a job that is only queued remotely holds no core and no memory of the process and blocks neither the engine's heavy lane nor a convergence check. A cancel of a running item kills its solver at once (the queue worker polls the marker every 100 ms).
* **`queue_list` shows the whole queue**: the engine's own jobs and the items other programs wrote, marked `foreign`, each with a `shared_id`. `queue_get`, `queue_cancel`, `queue_hold`, `queue_release` and `queue_set_priority` take `shared_id` instead of `id` for those; `queue_pause` and `queue_resume` also write the queue's `pause` marker. `queue_submit {detached: true}` writes an item and returns at once with nothing waiting for it; the answer names `results`, the file `open_results` reads once the item is `done`.
* **The command line** is `nyx queue submit|list|cancel|hold|release|priority|pause|resume|clear WORKSPACE ...` (also `nyx-worker queue ...`).
* **Which folder.** The window uses its workspace (`NYX_WORKSPACE`, else `Documents/nyx`); the MCP server uses `--dir` when given, and otherwise the window's workspace (its own default, `Documents/nyx_runs/mcp`, keeps the scripts and nothing else of the queue). The cost history is the queue folder's file for any engine that has one set.
* **Ordering never changes a number** (ADR 11): an item is a frozen project; nothing carries from one item to the next, which is also the rule that confines the warm starts of the sweep plan (P5b) to one sweep request.

## Not in this decision

* **No run record for an item that finishes while nothing watches it.** The results wait in the item (`results.h5`, removed after 24 hours if nobody collects them); `open_results` reads them. A front-end that is open records the run as before.
* **Convergence checks, mesh builds and `open_results` stay in the process** and are not coordinated across programs; two programs can still run one of these each beside the queue worker's solve.
* **A hold, a priority change or a hold-at-submit for a job that is still meshing** (before its item exists) is refused; `hold: true` at submit works.
* **The order of the meshing stage between jobs of one engine is first come, not by priority** (the gate is a lock).
* **Position in the queue** is the queue worker's ranking at its last look (a quarter of a second old).
* **No backfill** (one heavy item at a time, as in ADR 86), no permissions per owner, no Q4 (workers on other machines).

## Consequences

* Closing the window leaves its queued and running items to the queue worker; reopening it lists them.
* A window and an MCP server on one workspace run one solve at a time, in order, instead of two at once.
* A queue worker process may remain for up to 30 s after the last item; it holds the workspace folder open, so deleting a workspace right after a run can fail once.
* `engine_commands` gains arguments and answer fields (additive, provisional 1.17); `mcp_tools` is unchanged (the four tools reach the commands through the script).

## Code

`crates/nyx-engine/src/{shared_queue,worker,queue,engine,jobs,estimate,commands}.rs`, `crates/nyx-engine/src/bin/nyx-worker.rs`, `crates/app/src/{main,app}.rs`, `crates/nyx-mcp/src/main.rs`; tests `crates/nyx-engine/tests/shared_queue.rs`.
