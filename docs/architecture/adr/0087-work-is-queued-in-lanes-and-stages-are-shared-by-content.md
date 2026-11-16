# 87. Work is queued in lanes; stages are shared by content

Date: 2026-11-16

## Status

Accepted as a step of the productivity plan (Q1, `docs/NYX2_PRODUCTIVITY_PLAN.md` section 3a). **The number 87 is provisional**: take the next free ADR number at merge and renumber this file, the index row and the code comments. **Depends on ADR 84 (P1) and ADR 85 (P2)**: this bookmark is built on `productivity-2-jobs` and cannot merge without both. Builds on ADR 57 (the front-ends are thin), ADR 64 (cache keys, never persisted) and ADR 84 (the problem key).

## Context

The runner started every job on a thread of its own at once. The engine refused nothing, so a host could start three solves and they ran together on the one rayon pool with their memory added up; the window kept its own "one at a time" rule beside it. A request that came after another of the same project meshed again, and a convergence check re-solved the mesh the user had just solved.

## Decision

* **Jobs run in lanes.** The *heavy* lane (mesh, solve, convergence check, opening a results file: they hold every core and the memory) runs one job at a time. The *light* lane (writing a results file) runs beside it, two at a time. The lane is fixed by the kind of job; nothing guesses a size (estimates are plan step Q2).
* **A job submitted while its lane is busy waits in the queue, and the runner starts the next one itself** when a job ends, from the finishing job's own thread. No front-end has to poll for the queue to move; `wait` works on a queued job. The state `queued` is new, `JobInfo` carries `lane`, `owner`, `priority`, `held`, `position` and `queued_seconds`.
* **Order:** priority (1 a person's request, 2 an agent's, 3 sweeps and batches; the window's engine is the owner `gui` and defaults to 1, others to 2), with ageing (a waiting job competes one level higher after a minute, never above 1); then the owner served longest ago (round robin, so a 50-point sweep of an agent delays a person by at most one item); then submission order. The order never enters a number: a job carries a frozen copy of the project (it always did) and reads no state another job left (a test submits the same three solves in three orders and compares the numbers bit for bit).
* **Cancel, hold, supersede, pause.** A queued job that is cancelled never runs and finishes `Cancelled` at once; a request cancelled before it ran leaves no run record. `hold` keeps a queued job from starting until `release`; `supersede` (a key of the kind and targets) cancels a queued job of the same owner and key when a newer one arrives, which is what a new Compute of the same target does; `pause` stops anything new from starting.
* **Stages are shared by content.** The engine keeps the last two meshes by mesh key and the last eight solve summaries by `(problem key, chain)`. A solve or mesh request that finds its mesh there does not build it ("mesh: the same mesh was built by an earlier request" in its log); a convergence check whose coarse solve (the physics alone on the present mesh) is there uses its numbers and solves only the finer mesh (`coarse_reused: true`, and a log line). The cache key includes the **chain** (the targets of the job that solved it), because the problem key alone is not enough: a physics solved alone and the same physics solved with the physics it imports from share a key and not their numbers (the limit recorded in ADR 84). A cancelled request keeps the mesh it built if it had built it. The summaries are the numbers a check compares, not the fields.
* **A failure stops the chain after it, not before it.** `run_chain_partial` returns the physics solved before the one that failed; the engine stores them (they are valid results of their own problems), records the run as failed and returns the error. A cancel discards the chain as before.
* **Sweeps are a queue of points, and a point is an edit list.** `sweep_run` takes `points` (each a list of `{ cmd, args }`: any command whose access is `edit` except those that replace the project) or `series` (one command repeated with a field set to each value), applies them to copies of the project and checks them all before anything is submitted (a bad point fails the request, naming it, and nothing runs), then submits one solve per point to the heavy lane at priority 3. What varies is the case's to say: a geometry parameter, the wavelength in the study, a material, a mesh setting, several at once. Each point is a job with its own frozen project, uses the stage caches (points that share a mesh build it once), is recorded in the run table as its own run (`StoreLink::record_point` freezes the point's project in the run and does *not* commit it as the experiment's project) and never touches the engine's project, mesh or results: the table (`sweep_table`, optionally waiting `wait_s`) holds each point's numbers, unknowns, cells and run, and a result's fields are dropped after the run is written. `sweep_cancel` cancels the points that have not finished; a point cancelled before it ran leaves no record. A point's numbers equal those of the same edits made by hand (a test compares them bit for bit) and do not depend on the order of the points. The azimuthal sweep of capability R9 becomes a thin command over this at merge.
* **Commands.** `queue_submit` (solve, mesh or convergence), `queue_list`, `queue_get`, `queue_cancel`, `queue_hold`, `queue_release`, `queue_set_priority`, `queue_pause`, `queue_resume`, and `queue_wait` with a **required** `timeout_s`, which returns the job's entry with `waiting: true` at the timeout instead of blocking for ever. `job_wait` takes the same optional `timeout_s`. An entry says what was asked (`kind`, `targets`) and whether the project has changed since (`project_changed`).

## Not in this decision

* **Estimates and admission by memory** (plan Q2, with P4): the heavy lane's limit is one job, not a memory budget.
* **A queue shared by the window and an MCP server** (Q3): each process has its own engine and queue; two processes still run their heavy jobs together.
* **Stage graphs with refcounted nodes.** Requests share a mesh and a coarse solve through the content cache; a request is still one job, not a graph of stages that fail and cancel separately. The partial result of a failed chain is the one case that crosses that line.
* **Fair share by an owner's permissions** (an AI host may only cancel its own requests) and moving an item up or down by drag; the priority number is the order control.

## Consequences

* `start_solve` and the other starters return at once with an id as before, but the job may be queued; `job_status` says so. Code that assumed two jobs run at the same time (a test that holds one job while waiting for the next) needs a gate; none in the tree did.
* Two solves that a host starts with `wait: false` now run one after the other, in order of priority.

## Code

`crates/nyx-engine/src/{jobs,queue,sweep,engine,solve,convergence,commands,store_link}.rs`; tests `crates/nyx-engine/tests/{queue,sweep}.rs`.

## Contract

`engine_commands` 1.8 (provisional): the ten `queue_*` commands, `job_wait` takes `timeout_s`, the replies of `job_status` and `job_list` gain the queue fields, `convergence_check` gains `coarse_reused`. `mcp_tools` 2.4 (provisional): the same tools. `engine_commands` 1.9 and `mcp_tools` 2.5 (provisional): `sweep_run`, `sweep_table`, `sweep_cancel`.
