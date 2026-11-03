# 9. Every kernel is three stages behind one entry point

Date: 2026-11-03

## Status

Accepted (commits `fa6eff2`, `106e7f8`; table in `docs/KERNEL_PORT_STATUS.md`).

## Context

Kernels were single functions from mesh and spec to result. Swapping a linear-algebra backend (ADR 8), reusing an assembled system across a sweep, or exposing per-cell fields all required cutting them open.

## Decision

* `solve(mesh, spec, state)` / `solve_with(.., progress)` is the **one entry point for all physics**: it checks the element order and numbers (`validate`), dispatches by `Physics` (and by mesh dimension for EM waves), refuses a result containing NaN or infinity, and appends the plain-language stage log to `SolveResult::log`. A cancelled solve returns `Err("cancelled")`.
* Every kernel is **assemble -> solve -> fields**, and `solve_*` is their composition:
  1. *assemble*: matrices on the free dofs, constraints, loads, the material data later stages need;
  2. *solve*: the raw result (complex eigenpairs, a vector, transient states); the backend is chosen here (`spec.solver`, or the eigensolver the stage calls);
  3. *fields*: plot fields (node-averaged, as before) and `Solution::cell_fields`, the same quantities per cell and local node, not averaged, with wafer names and complex values where the quantity is complex.
* The node average of a cell field is the plot field of the same name (checked on every golden case).
* Splitting changed no number: the golden files are unchanged.

## Consequences

* A `femcore::solvers` adapter replaces one function per kernel.
* Stage function names are public API (`assemble_modes`, `solve_mode_pairs`, `mode_fields_cells`, ...); they are listed in contract `kernel_stages` because renaming one silently breaks callers outside the repository.
* The `Physics::Emw` dispatch by mesh dimension is a quiet rule (R10).
* `progress` is a thread-local, so stage functions called outside `solve_with` run unobserved and uncancellable (R19).

## Code

`crates/femcore/src/solve.rs`, `validate.rs`, `progress.rs`, `fields.rs`, and the kernels `dfs.rs`, `ems.rs`, `ems_pml.rs`, `ems_axisym.rs`, `ems3d.rs`, `sms.rs`, `sms_pml.rs`, `sms_floquet.rs`, `sms_lrb.rs`, `sms_kz.rs`, `sms_axisym.rs`, `sbs.rs`, `sbs_ring.rs`.

## Contract

`contracts/kernel_stages`: `crates/femcore/tests/stages.rs` (composition equals the one-call solver) plus `contracts.rs::kernel_stage_api_is_public`, which names every *public* stage function so removing one fails to compile. The stages of plain `sms`, `sms_lrb` and `ems3d` are `pub(crate)` (R26), so only 7 of 10 kernel families are reachable from outside the crate.
