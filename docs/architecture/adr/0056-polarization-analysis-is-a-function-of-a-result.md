# 56. The polarization analysis is a function of a result (`nyx-analysis`)

Date: 2026-11-05

## Status

Accepted (restructure phase 5, bookmark `restructure-5-analysis`). Applies the plan's phase 5 and resolves ROUGH_SPOTS R28; amends ADR 21 (the analysis lived in `femcore`) and ADR 45 and 47 (the mode kernels computed the analysis and the window got it through the result).

## Context

`femcore::{analysis, polarization}` computed the Stokes parameters, Poincare points, mode classes and degenerate groups of optical modes. The `ems` kernels called it at the end of a mode study, wrote one text scalar `polarization N` per mode and the table into the log, and handed the whole `Analysis` out next to the result for the window (`Solved`, `run_analysed`, `SolveResult::analysis` of the old window result). So a result was not just data: a front-end that was not the window, or a stored result, had no way to ask for the analysis, and the analysis could not be changed without touching a physics. The analysis also named its axes after the local frame of the mode solvers (x lateral, y normal, z propagation), so `Ey` meant two things (R28).

## Decision

1. **`nyx-analysis` holds the analysis.** `polarization` (Stokes, Poincare, classes, drawings) and `modes` (`ModeAnalysis`, the table, CSV, groups) moved from `femcore` with their tests; the numbers are unchanged. `polarization(mesh, result)` is the typed entry point; `PolarizationAnalysis` is the same behind the `Analysis` trait (plain JSON). The old struct `Analysis` is `ModeAnalysis` now, so the name is free for the trait. `nyx-analysis` depends on `nyx-physics` and `femcore` (mesh geometry and the quadrature rule); the layering table says so.
2. **A result carries what the analysis needs, not the analysis.** A mode study of `ems` puts the complex electric field of every mode in its cell fields (`Ex`, `Ey`, `Ez` in wafer axes; the 2D kernels already did, the cavity too) and states its **frame** in the scalar `field frame` (`nyx_physics::FieldFrame`): which wafer axes are the lateral, the normal and the propagation direction, the handedness of that triple, and the domains that are absorbing layers and are left out. The text is readable (`lateral x, normal z, propagation y, left-handed, outside domains 2`). A result without a frame is not analysed (`applies_to` is false); the analysis never guesses. The cavity result also has the complex frequency as `Solution::eigenvalue` (Hz), as the field's documentation always said it could.
3. **Sampling.** A cross-section is sampled at the 24 points of `femcore::optical::fine_rule` in every cell outside the absorbing domains, with the field interpolated from the cell nodes (P1 or P2). The mode fields are piecewise quadratic, so this is exact; the reproduced Stokes lines of the five mode goldens agree with the old ones (tests below). A cavity is sampled at the four vertices of every cell, each standing for a quarter of the volume, as before.
4. **The kernels stop computing it.** `Solved`, `run_analysed`, `with_analysis` and the `polarization N` scalars are gone from `nyx-ems`; `nyx_ems::run` returns the result. The table is no longer appended to the log of the new result. The unused public `axisym::mode_samples` and the `optical_set_*` helpers that only fed the analysis were deleted.
5. **The window.** `nyx_engine::legacy::solve_with_analysis` returns the old-shape result and the `ModeAnalysis` of an EM-wave study; `ResultView::analysis` holds it (the old `SolveResult::analysis` is gone). The bridge still appends the table to the old result's log, and sets `Solution::te_fraction` (one number, for the use cases, the run table and nothing else) until phase 7 deletes the old result type. The MCP session reads the analysis from the same call. The engine re-exports the crate as `nyx_engine::analysis` so the window needs no new dependency.
6. **R28.** The samples are `lat`, `nor`, `prop`: what the component is, not an axis number. The documentation maps them to the wafer axes x, z and y; the labels (`quasi-TE (x)`, `quasi-TM (z)`, the `Ey` column for the longitudinal share) were already wafer names and are unchanged, so the Analysis tab reads the same.
7. **Sign convention of the longitudinal field (a fix found on the way).** The real mode solver flips the transverse field so that its largest component is positive but did not flip `E_y` with it, so the relative sign of `E_y` between two modes was arbitrary. The old analysis evaluated the true field and did not see it; an analysis that reads the result does, in the overlap and the resolved pair of a degenerate group. `E_y` is now flipped together with the transverse field, in the plot field and in the cell field. This changes the sign of the `Ey` probe of 6 modes in three goldens (below) and nothing else.

## Consequences

* The five mode goldens lose their `polarization N` scalars and gain `field frame`. The polarization lines they used to hold are kept unchanged in `crates/nyx-ems/tests/golden/polarization.json` and the analysis is checked against them (class as text, rounded numbers to 1.5e-3, 5e-3 for the silicon-wire files at the G2 tolerances). `GOLDEN_CHANGES.md` records both changes.
* Any front-end can run the analysis on a stored result and its mesh; a new physics with complex vector fields can state a frame and get the same analysis.
* `femcore` has no analysis code left; `Solution::te_fraction` and `legacy::solve_with_analysis` go with the old result type in phase 7.
* The analysis table is not in the log of a new-shape result any more; the window's log is unchanged.

## Code

`crates/nyx-analysis/src/{lib,polarization,modes,samples}.rs`, `crates/nyx-physics/src/frame.rs`, `crates/nyx-ems/src/{modes,modes_complex,cavity,optical,lib}.rs`, `crates/nyx-engine/src/legacy/{mod,ems}.rs`, `crates/app/src/{worker,state,app,tabs}.rs`, `crates/nyx-session/src/tools.rs`.

## Contract

`physics_ems_result` 1.1 (scalar `field frame` instead of `polarization N`; the cavity solutions have an eigenvalue). The frame text is documented in `nyx_physics::frame`. No exemplar of its own yet: `analysis_api` is versioned when a second analysis exists.
