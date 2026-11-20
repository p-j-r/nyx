# Using nyx from an AI host: what costs time and what to set

For the people who connect nyx to an AI host (Antigravity, Claude Code, Claude Desktop, others) and for the agents that run there. It comes from the benchmark session of 2026-11-19 (two silicon-on-sapphire waveguides, `LATENCY_AND_WORKFLOW_ANALYSIS.md` in the benchmark workspace): the solves that gave the final numbers took 3.8 s and 4.2 s; the sessions took 35 and 13 minutes. Over 90 % of the time was waiting and model turns, not solving. ADR 97 records what nyx changed because of it.

## What nyx does since L1 (ADR 97)

- **A call never waits longer than 90 s** (`--max-wait SECONDS` changes it, 1 to 900). `wait_s` above the limit is cut to it and the answer says so. A host whose tool timeout is shorter than 90 s must be given `--max-wait` below its timeout, in the `args` of its MCP configuration. Antigravity ends a call after 180 s.
- **A waiting call can be ended.** While `script_run` or `script_wait` waits for a solve, the server reads what the host sends: a newer `script_run` or `script_wait` ends the wait at once (the first call is answered `running`, with a note), `notifications/cancelled` for the waiting call stops the wait and sends no answer, `ping` and `tools/list` are answered during the wait. The solve itself goes on in its worker process and the next `script_wait` picks it up. Before this a call the host had given up on kept the server blocked for its whole wait, and the follow-up calls queued behind it timed out one after the other (12 minutes in the benchmark).
- **Progress.** A call that carries `_meta.progressToken` gets `notifications/progress` every 5 s while it waits; hosts that restart their timeout on progress then never time out.
- **`running` answers say where the job is:** `job` has its stage, seconds, fraction, queue place and the estimate of what such a solve took on this machine (absent when there is none).
- **`mesh_build` answers `cost`** (cells, expected seconds and memory from earlier solves, a warning for a large mesh). `convergence_check` takes `dry_run` (the cells of both meshes, nothing solved), `max_cells` and `max_seconds` (refused before any solve).

## What an agent should do

1. **One `script_run` per task.** Build the geometry, `mesh_build`, `solve_run`, the reductions (`result_info`, `analyse`, `confinement`, `probe`) and `print` in one script; save it with `save_as` and run it again with `params`. In the benchmark 14 separate calls paid 5 to 18 s of model time each, which was more than all the solves.
2. **Read the guides in one batch.** `api` (instructions), `api topic=geometry`, `api topic=optical-modes` and the commands you need are independent: ask for them in one turn, not one per turn (the benchmark used 7 turns, 2.5 minutes, before the first question).
3. **Look at `cost` before the solve.** A mesh of 22 714 cells for a 2D mode took minutes to solve and `convergence_check` at 0.6 would have made it about 170 000 cells. A 2D mode needs 2 000 to 4 000 cells. Run `convergence_check` with `dry_run` first, or with `max_cells`.
4. **Keep `wait_s` under the host's timeout** (at most 90 is enforced) and answer `running` with `script_wait`, not with a new `script_run`.
5. **Do not scan when you can aim.** In a rib guide solve `n_modes` 1 at each `n_target` instead of eight modes five times.
6. **Do not quote noise.** A lossless guide on a tilted crystal has an imaginary part near 1e-10 and a Q near 1e10: `convergence_check` marks those rows `noise`.

## What to choose in the host

- A **single model that executes** is faster for this work than a planner and an executor in turn: the benchmark spent 4.5 minutes streaming plans and 2 minutes on two hand-backs, in the session with the planner.
- A host that lets you **raise the tool timeout** should have it raised only as far as `--max-wait` plus a margin: nyx answers `running` and the host loses nothing by asking again.
- The server needs `nyx-script-worker` next to `nyx-mcp`, and its folder (`--dir`) is where workflows and saved scripts go.
