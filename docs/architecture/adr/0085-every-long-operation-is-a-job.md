# 85. Every long operation is a job that can be watched and cancelled

Date: 2026-11-16

## Status

Accepted as a step of the productivity plan (P2, `docs/NYX2_PRODUCTIVITY_PLAN.md`). **Number 85: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.** **Depends on ADR 84 (P1)**: this bookmark is built on `productivity-1-stale` and cannot merge without it. Builds on ADR 57 (the front-ends are thin) and ADR 61 (convergence is an engine command).

## Context

Measured in P0 (`docs/PRODUCTIVITY_STATUS.md`): the window froze for 0.35 to 1.3 s after every recorded run while it wrote `results.h5`, and for 0.1 to 0.3 s in `open_results` (read, a blocking re-mesh, compare); a mesh job could not be cancelled and logged no stage; a mesh convergence check ran on a bare thread with no cancel, no progress and no place in the job list, beside whatever else the window started.

## Decision

* **A mesh is a cancellable job with stages.** `start_mesh` installs the job's `Progress` and `FemcoreMesher::mesh` opens a stage; the adaptive drivers (2D and 3D, each iteration and each split cascade) and the 3D lattice mesher (between its phases) check the cancel flag. Checks are on the driving thread, not inside parallel passes (worker threads do not see the progress), so a cancel lands at the next pass boundary. A cancelled mesh is `Cancelled`, not an invalid model.
* **A convergence check is a job** (`Engine::start_convergence`, `convergence_check {wait: false}`). Its two solves run on the job's thread under its `Progress`, as two stages, through the same `run_chain` as a solve; it is listed by `job_list`, reports its stage and is cancelled by `job_cancel`. The synchronous command (default `wait: true`) is unchanged for hosts.
* **Opening a results file is a job** (`start_open_results`, `open_results {wait: false}`): reading, meshing again and the comparison with the stored mesh run on the job; the engine is not touched until the job is applied, then the file's project replaces the project and the stored results become the latest results (stored under the keys of ADR 84). A failed or cancelled open changes nothing, with no restore step. The window makes its draft workflow when the job ends. (The plan proposed showing the stored results against the stored mesh at once and re-meshing as the job; the whole open as a job removes the freeze without the viewer having to draw from a mesh no model exists for.)
* **Writing a results file is a job** (`start_export_results`, `export_results {wait: false}`) that writes a snapshot of the project, mesh and results of the moment of the request, to a temporary name that is renamed when the file is complete. The window starts it after each recorded run and goes on; a failure shows when the job ends.
* **Honest cancel.** `JobInfo` carries `cancel_requested` and `stage_seconds` (job_status and job_list), and the toolbar shows the seconds in the stage the job waits for; a cancel lasts as long as the stage it is in cannot be interrupted (a factorisation), and nothing estimates it.
* **A job started for a project that was replaced is dropped.** The engine counts the projects it has held (`set_project`); a mesh or solve that ends for an earlier one returns an error and changes nothing (it would otherwise store its result into the new project, ADR 84's rule for an empty slot).
* **Gate:** a test reads the window's sources and fails when they call `mesh_now`, `solve_now`, `engine.wait`, `open_results`, `export_results` or `convergence_check` (the headless `batch.rs` may wait).

## Not in this decision

* **A queue.** The window keeps its one running job and its one waiting slot; running several jobs in order, shared between the window and a host, is the scheduler of plan step Q1 (ADR [87]).
* **Cancel inside a factorisation or an assembly pass**: not possible in process; the way to stop one at once is a child process (plan step P4).
* **Waiting for a pending results file when the window closes**: the temporary name means a closed window leaves no file that looks complete, but the file of the last run can be missing.

## Consequences

* The window can no longer be frozen by the engine on its own thread; the first 3D freeze measured in P0 is gone.
* `job_status` and `job_list` gain two fields; `convergence_check`, `export_results` and `open_results` gain `wait`; `solve_run`-style replies are unchanged.
* Edits made in the window while a results file is being opened are lost when the file's project replaces it (the old window could not be edited during the open).

## Code

`crates/nyx-engine/src/{engine,jobs,convergence,results_file,meshing,commands}.rs`, `crates/femcore/src/adapt/{drive,drive3}.rs`, `crates/femcore/src/meshing/mesher3d.rs`, `crates/app/src/{app,shell,tabs}.rs`; tests `crates/nyx-engine/tests/jobs.rs` and the window tests in `app.rs`.

## Contract

`engine_commands` 1.10: `convergence_check`, `export_results`, `open_results` take `wait`; replies of `job_status` and `job_list` gain `cancel_requested` and `stage_seconds`. `mcp_tools` 2.5: the same three tools take `wait`.
