# 96. A sweep may start each point from the modes of the one before

Date: 2026-11-19

## Status

Accepted as step P5b of the productivity plan (`docs/NYX2_PRODUCTIVITY_PLAN.md`, section 3, "Reuse between solves"), with Paul's "yes build p5b" of 2026-11-19 after the measurement in `docs/PRODUCTIVITY_STATUS.md`. **Number 96: settled at the merge of the productivity line with main (2026-11-19), the next free number after the capability line.** Builds on ADR 87 and ADR 95 (the queue never lets one item change another's numbers) and ADR 91 (validated is a state).

## Context

A sweep that keeps its mesh solves the same eigenproblem many times with a small change. Each solve starts its shift-and-invert eigensolver from a pseudo-random block, and the iterations that follow take most of the time that is not the factorisation. Measured with a throw-away probe (status file, "P5b measured"): the previous point's wanted modes as the start block cut the eigensolve by about 40 % where the modes move, and a solve by **12 % to 20 %**; a sweep that changes the mesh (a width, a wavelength that sizes the mesh) has no start block to reuse and gains nothing. Plan P5a (bit-identical reuse of the symbolic factorisation) was measured at 1 to 14 % and skipped.

## Decision

* **Opt-in: `sweep_run {reuse: "warm"}`.** Without it nothing changes, bit for bit: no context is installed and the eigensolvers take the path they always took (`tests/worker.rs`, the goldens and `sweep.rs` run unchanged). There is no warm start outside a sweep; a request that is not a point of a warm sweep has no link, so **nothing ever carries between queue items** (ADR 95's rule).
* **The numerical core** (`femcore::warm`): a `Warm` context installed on the solving thread, like `progress`. The two shift-and-invert eigensolvers (`driver::eigs_near`, the Krylov-Schur path of the 3D cavity and the complex pencils, and `arnoldi_real::eigs_general`, the 2D mode solver) ask for a `Ticket`; the i-th eigensolve of a solve takes the i-th entry left by the previous solve if it is the same kind of solve on the same number of unknowns for the same number of modes, and starts cold otherwise (never an error: a start block can make a solve faster or slower, not wrong). `krylov_schur` takes a start of several columns (`given.min(block)`; one column, or none, is the block it always had). **All `nev` columns are warm**: replacing even one with a random column put the run back to cold speed in the measurement.
* **The chain.** Point i's link names the cell of point i-1 and the job that solves it. A point **waits** for its predecessor only in the `shared` lane (ADR 95), where it holds no core; in the heavy lane the lane runs points in order, so an unfinished predecessor is one that cannot finish (held, cancelled): the point is cold and says why (waiting would hold the lane for ever; tested). A point's cell is filled when its job ends, however it ends (`CellGuard`), so a point never waits for one that failed or never ran. Only a point that finished its whole chain hands its modes on. A warm sweep through the shared queue submits one point at a time (each waits for the one before), so the queue shows one point of it at a time.
* **The vectors travel as bytes** (`femcore::warm::to_bytes`): to a worker process in the job folder (`warm_in.bin` in, `warm_out.bin` and `warm_calls.json` out), and so also to the queue worker, whose item folder is a job folder. The three ways a solve runs give one answer (tested: the numbers are equal to the bit between in-process and worker, and equal through the queue worker).
* **The record.** A point's row and its run record (`reuse`, optional in `nyx_store::RunRecord`) say: `warm`, `from_point`, the restarts or cycles of every eigensolve, and for a warm point `exploratory: true` with the note that it agrees with a cold solve to the solver's tolerance, not bit for bit. A point that started cold says `why_cold` (the first point, a predecessor that had not finished, one that left no modes, nothing of the same size to start from).
* **A warm solve is never a cached coarse solve.** The engine does not keep a warm point's summary in the stage cache, so a mesh convergence check of the same problem solves its coarse mesh cold.

## Not in this decision (stated)

* **A warm result is not marked `exploratory` when it is opened from its run** (`open_results`): the results file does not carry `reuse`; the run record and the table row do. `validated` of ADR 91 applies to the engine's own results, which a sweep never replaces.
* **A start block of the old modes can miss a mode that enters the wanted window and is orthogonal to them.** In the measurement rounding noise let every entering mode be found (in as many restarts as a cold solve), but it is not guaranteed; that is why a warm point is not confirmed until solved cold. No second random block is added (it gave back the whole gain in the measurement).
* **Only the shift-and-invert eigensolvers are covered**: not LOBPCG, the direct solves of the static and transient studies, or the cavity solvers that do not use `eigs_near`. The factorisation is not reused (a new matrix at every point).
* **The window's sweep section sweeps geometry parameters, which change the mesh**, so it has no warm option; the command, scripts and the MCP server have.
* **Parallel warm points are not possible** (each needs the one before); a warm sweep is a chain.
* **No mode tracking** (R11) and no shifts taken from the previous eigenvalues: the shift stays what the study says, so a point's result does not depend on where it started except through the iteration.

## Consequences

* A warm sweep of a case that keeps its mesh solves 12 % to 20 % faster per point after the first; one that does not keep its mesh runs as a cold sweep and says why in every row.
* `engine_commands` gains the `reuse` argument of `sweep_run` and a `reuse` object in the rows of its table (additive, provisional 1.18). The run record gains an optional `reuse`.

## Code

`crates/femcore/src/{warm,linalg/driver,linalg/arnoldi_real,linalg/eps/krylov_schur}.rs`, `crates/nyx-engine/src/{reuse,engine,sweep,worker,queue,store_link,commands}.rs`, `crates/nyx-store/src/runs.rs`; tests `crates/femcore/tests/warm.rs`, `crates/nyx-engine/tests/warm_sweep.rs`.
