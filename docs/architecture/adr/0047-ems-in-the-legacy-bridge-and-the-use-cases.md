# 47. `ems` in the legacy bridge, and the use cases take their solver as an argument

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-ems`). Applies ADR 41, item 6, to `ems`; amends ADR 22 (use cases live in `femcore`) for how they solve.

## Context

The window and the MCP session still build the old `ProblemSpec` for EM waves and read the old `SolveResult`, including its `analysis` (the Stokes parameters and the Poincare sphere). The old problem has no EM study kind: a 2D mesh meant waveguide modes, a 3D mesh meant a cavity. The use cases `usecases::soi` and `usecases::bragg` solved by calling `femcore::solve::solve`, which cannot reach `nyx-ems`; they are used by the window (the example projects), the MCP session (`slab_neff`, wire parameters), tests and examples.

## Decision

1. **`nyx_engine::legacy::ems`** converts the old problem to `nyx_ems::Problem` (`problem_from_old`): a 2D mesh with a `Modes` study is `Study::Modes`, a 3D mesh with an `Eigenfrequency` study is `Study::Cavity`; any other combination answers the text the old dispatch gave ("mode analysis needs a Modes study", "EM waves on a 3D mesh use the Eigenfrequency study"). Materials (scalars, crystals, PML) and Floquet pairs are copied. It runs `nyx_ems::run_analysed`, converts the result with `result_back` and puts the `Analysis` into `SolveResult::analysis`. `legacy::solve_with` routes `Physics::Emw` there; `femcore::solve` answers an error.
2. **Use cases get a solver.** `usecases::Solver` is a function from a mesh and an old spec to an old result. `SoiWire::solve`, `BraggGrating::{solve, band_edges, band_edges_of, bragg_frequency_estimate, matched_period}`, `fundamental_neff` and `group_index` take one. The window, the tests and the examples pass `nyx_engine::legacy::solve`. The geometry, mesh and spec parts of the use cases are unchanged and still in `femcore`, so the `nyx-ems` golden tests can use `SoiWire` without the bridge. The use cases and the bridge are deleted together in phase 7 (the plan retires `usecases/`).
3. **Where the mixed tests went.** `usecases.rs` (it needs the old result with its analysis) is `crates/nyx-engine/tests/usecases.rs`; the studies and benchmarks `soi_modes`, `bragg_gratings`, `cavity_bench` and `complex_bench` are `crates/nyx-engine/examples/` (`tools/baseline.sh` runs `cavity_bench` from there). Every other ems test is in `crates/nyx-ems/tests/` and runs the typed problem directly.

## Consequences

* The window behaves as before: same results, same analysis tab, same log text.
* `femcore` has no EM solver left; `femcore::solve` is down to `sms`.
* Two result shapes still exist until phase 7, as for `dfs`.

## Code

`crates/nyx-engine/src/legacy/ems.rs`, `crates/femcore/src/usecases/*`, `crates/nyx-engine/tests/usecases.rs`, `crates/nyx-engine/examples/*`.

## Contract

`kernel_stages` 1.3 (the ems stages are checked in `crates/nyx-ems/tests/contracts.rs`).
