# 83. Size rules on faces and `max_dev` act in the 3D lattice mesher

Date: 2026-11-15

## Status

Implemented on `capability-r12-curved-elements`.

## Context

`GeomDomain3` kept only the `Size` rules on regions, and `max_dev` was read only by the adaptive path, so `add_size` on a face and `max_dev` were accepted and silently ignored in 3D. Separately, the `set_mesh` command rejected `max_dev`, `structured_pml` and `points_per_wavelength` as unknown fields although its schema lists them and the window sends them (found while doing this work, fixed here).

## Decision

* Per face, `h_face = min(rule h, 2 sqrt(2 r max_dev))`, the second term for a face with a radius of curvature `r` (`facet_radius`, the formula of `adapt3d`); `max_dev` comes from a rule on the face or from the request (`GeomDomain3::with_max_dev`). With no `max_dev` the lattice path has no chord limit, so existing meshes do not change (the adaptive path keeps its `h / 50` default).
* `Domain3::size_at(p)` returns the smallest `h_face + 0.3 d` over faces at distance `d` from `p` (growth 0.3 per length: a ratio of 1.3 between neighbours); step 8 of `mesh_3d` uses `min(region size, size_at)` and runs when any face asks for less than the global size (`asks_size_at`). New boundary vertices are placed on the exact face by `project_edge`.
* A `Size` rule on crease edges or vertices has no meaning in 3D, `min_layers` has no effect, and a `max_deviation` on a region cannot apply: `MeshReport.warnings` names each instead of silence.
* A face size that reaches a periodic face is refused with a message (the pair would be refined differently); congruent periodic faces are R30 of `capability-r7-periodic-cells`.
* `set_mesh` accepts `max_dev`, `structured_pml`, `points_per_wavelength` and `curved`.

## Consequences

Meshes without a face size or `max_dev` are unchanged (tested). `mesh_report` 1.1 gains `warnings` and `curved`.

## Code

`crates/femcore/src/meshing/geomdrive3.rs`, `mesher3d.rs`, `request.rs`; `crates/nyx-engine/src/commands.rs`.
