# 97. The MCP server reads while it waits, and says what a job costs

Date: 2026-11-20

## Status

Accepted as the latency work L1 to L3 (Paul, 2026-11-20: "do cap and the reader; yes fix the convergences; yes estimated times in all of them are useful; check the wildcard issue"). Amends ADR 94 (the waits of `script_run` and `script_wait`) and ADR 91 (what a convergence check counts). The number was free after ADR 96.

## Context

A benchmark session on 2026-11-19 (two silicon-on-sapphire waveguides driven from Antigravity through `nyx-mcp`) needed 3.8 s and 4.2 s of solving for the numbers it reported and took 35 and 13 minutes. What nyx contributed:

1. `script_run` accepted `wait_s` up to 900 and the server handled one line of stdin at a time. The host ends a call after 180 s; the server was still inside its 600 s wait, so the three calls the agent sent next queued behind it and each timed out after 180 s: 12 minutes.
2. `convergence_check` compares every number by its relative difference. A lossless guide on a tilted crystal runs the complex solver and has an imaginary part of 1e-10 that changes by 100 % between two meshes, and a quality factor of 1e10 made of it: `converged` stayed false and the agent went back to its planner twice.
3. A mesh of 22 714 cells for a 2D mode (a high `points_per_wavelength` with a small `adapt_scale`) was found out by solving it; a convergence check at 0.6 would have made it about 170 000.
4. A material given for `Air*` made `mesh_build` with `physics_sizing` fail ("the sizing names region 'Air*', which the model does not have"): the solver expanded the star, the mesher's wave sizing did not.

The rest (the model turns, the two-model planning, fourteen calls where one script would do) belongs to the host and the agent: `docs/HOST_SETUP.md` and the guide cards say what to do.

## Decision

* **L1, the reader.** `nyx-mcp` reads stdin on a thread of its own. The engine and the script kernel stay on the main thread: a wait is made of slices of 200 ms (`Kernel::wait` with a short wait) and between two slices the server looks at what arrived. A newer `script_run` or `script_wait` ends the wait of the call before it (answered `running`, with a note; the job goes on, the next `script_wait` picks it up); `notifications/cancelled` for the waiting call stops the wait and sends no answer; `ping` and `tools/list` are answered during a wait; a request cancelled before its turn is dropped. Every other message waits for its turn, as before. The printed lines of all slices are kept and a command is counted once.
* **L1, the limit.** A call waits at most `--max-wait` seconds (default 90, 1 to 900). A longer `wait_s` is cut and the answer says `wait_s: {asked, used}`; it is not an error, because an error costs the agent a turn. The tool schema carries the maximum. A call with `_meta.progressToken` gets `notifications/progress` every 5 s.
* **L1, what a `running` answer says.** The job the script waits for: stage, seconds in it, fraction, queue place and the estimate (`Kernel::waiting_job`, `Engine::queue_entry`).
* **L2, noise in a convergence check.** A mode whose `eigenvalue[i].im` is below 1e-9 times its real part on both meshes (a quality factor above 5e8 is not told from rounding) has its `im` row and its quality-factor rows (`Q<n>`, `Q <n> = Re/(2|Im|) of neff`, `f im <n> (Hz)`) marked `noise: true`; they are listed, passed when there is a tolerance and not counted in `worst_relative` or `converged`. The answer lists `noise_modes`. A mode with loss is never noise: a cavity with Q 1e6 has `im/re` 5e-7. Chosen over comparing `im` against the complex magnitude, which would hide a loss figure that does not converge.
* **L3, estimates.** `mesh_build` answers `cost` (cells, dimension, the estimate of a solve of the default chain from the solves this machine has run: seconds, memory, runs, extrapolated; `null` without a basis, no number is invented; a warning above 15 000 cells in 2D and 100 000 in 3D). `convergence_check` takes `dry_run` (both meshes are built, nothing is solved: the cells of each, their cost, the factor between them), `max_cells` and `max_seconds` (refused before any solve, after the finer mesh is built; no refusal on time without an estimate). The finer mesh is built once and solved on, not meshed again.
* **L3, the star.** `Sizing::expand_wildcards`: a name ending in `*` stands for the region and its pieces `name.*` in the label map of the model; a name given in full keeps its own wavelength. Applied where the mesher reads the sizing (`mesh_for_sizing`, `wave_sizes_by_region`).
* **Cards.** `nyx-instructions` (one script per task, `wait_s`, `cost`, `dry_run`), `troubleshooting` (what makes a mesh too large, the noise rows, a loop variable), `optical-modes` (aim `n_target` in a rib guide, the energy fraction in a tensor crystal). The card limit of 3600 bytes holds.

## Not in this decision

* **No estimate from a model of the mesher.** `cost.solve` is null until this machine has solved a few meshes of the same chain and dimension; the history does not yet know the complex solver from the real one, so a fit across both is wide (the answer carries the number of runs and says when it extrapolates).
* **Concurrent tool calls.** A host that sends `api` or `script_check` while a script waits gets them after the wait (at most the limit); only `script_run` and `script_wait` end a wait.
* **A solve is not made faster.** The wait is shorter for the host; the solve takes what it takes.
* **Connecting a host such as Antigravity** (`nyx-connect` has no entry for it) is a separate decision.
