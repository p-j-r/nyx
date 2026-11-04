# 51. Adaptive meshing of 3D geom models: `adapt3` with the exact faces as the oracle

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). Lifts "the adaptive mesher and error estimators are not extended to 3D" of ADR 40. Uses ADR 11 (determinism), ADR 18 (no curved elements), ADR 40 (lattice mesher), ADR 50 (wave sizes).

## Context

`Mesher::Adaptive` meant the 2D metric adaptation. In 3D the pieces existed apart: the start mesh (ADR 40), `adapt::drive3::adapt3` (split, collapse, flip, edge removal, smoothing on tetrahedra, with surfaces, creases and corners classified from the start mesh itself) and an oracle trait `Geometry3` that only knew triangles of the start mesh or a hand-given sphere projector. Run on a curved model, `adapt3` shrank the body (every new surface vertex lies on a chord) and it never carried the interface facets of ADR 40 through.

## Decision

`Mesher::Adaptive` on a 3D model is: lattice start mesh, a vertex metric, `adapt3_with` (`meshing::adapt3d`).

* **The oracle** (`GeomOracle3`, `FaceTable`). A surface vertex slides over its patch and a new one lands on the nearest *exact* geom face among the faces of the patch's signature (`Complex::facet_project`: plane, sphere, torus, cylinder wall); a boundary patch takes the faces of its boundary number, an interface patch the faces between exactly its two domains. Creases and corners keep the polylines of the start mesh: straight creases (boxes, layer stacks, ribs) are exact, curved ones as good as the start mesh, like the lattice mesher (ADR 40).
* **Creases are by label, not by angle.** `TopoOptions3::crease_angle_deg` is 120 for geom models: faces with different signatures always meet in a crease, but a coarse facetted sphere (dihedral angles above 40 degrees) is not a crease that pins its vertices to the chords. A same-label fold stays safe because the oracle projects onto the nearest exact face of that label.
* **Interface facets.** `Work3::to_mesh` returns outer boundary facets only and `Work3::from_mesh` rejects interior facets. The interface facets of the lattice mesh are removed before and written back after, labelled by the nearest geom face between the two domains (`FaceTable::interface_label`, the rule of `GeomDomain3::boundary_label`), oriented out of the lower-numbered domain.
* **Sizes** (`vertex_sizes3`, `a_priori_metrics3`, `geometry_floor3`): region sizes (global, `Size` rules, layers of absorbing pieces), the chord limit `2 sqrt(2 r max_dev)` on spheres, and the wave sizes of ADR 50. `geometry_floor3` is the part the geometry insists on, the floor of an estimator metric (ADR 53). Size rules on edges are not meaningful for the lattice and stay ignored.
* **Options for geom models** (`Adaptive` maps to `AdaptOptions3`): at most 12 iterations, stop at 2 % of the tetrahedra changed per iteration (a curved surface never lets splits and collapses cancel, so the defaults of 30 and 0.5 % are spent on churn), and `sweep_frac` 0.3 % (see below).
* **Fixes in `adapt3` found on the way** (all in `adapt::{ops3,topo3,drive3}`; the default behaviour of the 3D loops is unchanged except where it was wrong):
  * `edge_removal_pass` could turn a triangle of the ring inside out (the ends of the edge still lay on its two sides), which overlaps the new tetrahedra; the triangles must now all face the way the ring does. It could also create a diagonal that exists elsewhere in the mesh, and a face with four tetrahedra. Edges and faces made earlier in the same pass are tracked (the snapshot does not know them), also in `flip_pass` for the 3-2 and 2-3 flips.
  * `Work3::check` (run after every pass with `FEMMESH_CHECK=1`) now also finds *overlap*: the two tetrahedra of an interior face must lie on opposite sides of it. Positive volumes and at most two tetrahedra per face do not exclude a folded mesh (the symptom was a total volume of 65.95 for a box of 64).
  * Dead tetrahedra were never removed from the working arrays, so every sweep walked an array several times the size of the mesh (`Work3::compact_tets` at the start of every loop of a sweep). `AdaptOptions3::sweep_frac` (default 0: the old rule, repeat until a sweep changes nothing) lets a sweep loop stop when it changes at most that fraction of the tetrahedra; every sweep rebuilds the incidence of the whole mesh.
* Labels pass through unchanged: domains, boundary facets and (rewritten) interface facets keep the numbers of the model's `LabelMap` (contract `mesh_geom3d_adaptive_and_structured_meshes_keep_the_label_map`).

## Consequences

* A sphere of radius 0.9 in a box refined from `h` 0.6 to 0.3: every vertex of the sphere is on the sphere to 1e-6 (the run without the oracle shrinks the ball to 64 % of its volume), the box faces are exact to 1e-14, the total volume is 64 to 1e-14, the mesh is identical for 1, 3 and 8 threads.
* A strongly graded case (wave size 0.15 in a ball inside a box with the start mesh at 0.8, 17 000 tetrahedra) takes about 30 s on this machine (183 s before the fixes and the sweep rule); the passes are still serial (`MISSING_FEATURES.md` T4) and rebuild the incidence per sweep, which is the next thing to make cheaper. The worst quality at the end is about 1e-3 on the sphere (surface flips on curved patches, T1, are not done).
* Not done: curved creases from the exact edges (the polylines of the start mesh), surface flips on curved patches, threaded passes, mixed sweep axes (geom).

## Code

`crates/femcore/src/meshing/adapt3d.rs`, dispatch in `meshing/mod.rs` (`mesh_3d_plan`); fixes in `crates/femcore/src/adapt/{ops3,topo3,drive3}.rs`; the loop that uses it with an estimator in ADR 53.

## Contract

`mesh_geom3d_adaptive_and_structured_meshes_keep_the_label_map` (`tests/contracts.rs`); behaviour in `tests/mesh3d_features.rs` (`the_adaptive_mesher_*`, `adaptation_is_deterministic_*`, `meshes_are_identical_for_any_thread_count`).
