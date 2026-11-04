# 50. Wavelength sizing in 3D: one element size per region, asked from the neutral `Sizing`

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). Carries ADR 32 step 6 (the neutral `femcore::meshing::Sizing`) into 3D; lifts the "wavelength sizing is not extended to 3D" item of ADR 40. Paul: adaptive meshing, wavelength sizing and structured PML are essential in 3D.

## Context

A physics reports the shortest wavelength it wants resolved in every region (`Sizing { order, headroom, wavelengths: name -> metres }`); the mesher turns it into element sizes with its points-per-wavelength calibration (`physics_size::cell_sizes_for`: 12 per wavelength at order 1, 6 at order 2, relaxed by `pml_relax` in absorbing pieces). In 2D `mesh_for_problem` did this for a `ProblemSpec`, with a global size (simple mesher) or a metric (adaptive mesher). In 3D it answered an error, because it called the 2D mesher; and `physics_size::metrics` is 2D only.

## Decision

* **`meshing::mesh_for_sizing(list, plan, sizing, sz)`** takes the neutral `Sizing` and works in either dimension. Region names are resolved with the label map of the model (a name the model does not have is an error, never silence); absorbing pieces are known from the model (`Evaluated::pml_regions`). In 2D it is `mesh_for_problem` without a `ProblemSpec` (both now share `mesh_waves_2d`, which changes no number). `mesh_for_problem` answers a clear error for a 3D model and points to `mesh_for_sizing`.
* **A size per region, not one global size.** In 3D the lattice mesher already refines a region that asks for a finer element than the global one (Rivara bisection, ADR 40). The wave size of a region is therefore handed to `GeomDomain3` as one more request on the region (`Sizes3::region_h`, finer only: the plan's `h` stays the ceiling, as in 2D), next to the `Size` rules and the layers of an absorbing piece. A cavity with a glass ball (wavelength 0.6, `h` 0.1) in air (wavelength 2.4, `h` 0.4) gets 0.1 in the ball and 0.4 in the air; 2D could only use the smallest size everywhere with the simple mesher.
* **The lattice is made at the coarsest size any meshed region needs** (never above the plan's `h`), and the finer regions are refined from it. A study whose waves are short everywhere is meshed at the right size at once instead of at a coarse size and bisected afterwards.
* **Adaptive mesher.** The start mesh is made from the model's own sizes (as in 2D); the metric is the intersection of the geometry metric (`adapt3d::vertex_sizes3`) and the wave size, both times `Adaptive::scale`, so `scale` stays one knob (ADR 51).
* Not done: curvature or corner sizing of the waves (calibration stays on the mesher side, `MISSING_FEATURES.md` E6); anisotropic sizes (the metric is isotropic, as in 2D).

## Consequences

* The engine (`FemcoreMesher::mesh`) can drop its "physics sizing is unsupported" answer for 3D models once it calls `mesh_for_sizing` or `mesh_request` (ADR 54); not done here (outside `meshing` and `adapt`).
* Measured (`tests/mesh3d_features.rs`): the ball above has 20 times the cells of the unsized mesh, the mean longest edge in the ball is below 0.13 for a wanted 0.1 (bisection ratio 1.3), the air stays at 0.2 or more; volumes are as exact as before.

## Code

`crates/femcore/src/meshing/{mod,geomdrive3}.rs` (`mesh_for_sizing`, `wave_sizes_by_region`, `Sizes3`, `GeomDomain3::new_with`).

## Contract

`mesh_request_and_report_match_their_exemplars_and_schemas` (the request carries the `Sizing`), behaviour in `crates/femcore/tests/mesh3d_features.rs` (`wavelength_sizing_*`, `absorbing_pieces_are_relaxed_by_pml_relax`).
