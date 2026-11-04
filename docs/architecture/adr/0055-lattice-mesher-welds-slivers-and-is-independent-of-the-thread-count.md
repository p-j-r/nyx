# 55. The lattice mesher welds rounding-distance vertices and does not depend on the thread count

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). Amends ADR 40, which listed "a test that the mesh is identical for any `NYX_THREADS` (ADR 11)" as still to add. Found while meshing layer stacks and absorbing windows (ADR 52).

## Context

1. **Slivers.** A surface of the model that lies on a lattice plane (a layer interface at `z = 0`, the faces of an absorbing window) makes the lattice vertices on it count as warped onto the surface, and the lattice edges next to them are also cut, so a cut point can be a few ulp from a warped vertex. The tetrahedra between the two have volume 1e-16 to 1e-30. A box of air over glass at `h` 0.5 had 140 of 3300; their boundary triangles are zero-area, which turned the prisms of the structured layers inside out (`build_pml`: "inverted, scaled Jacobian 0"). The same cells would reach any solver.
2. **Thread count.** Two things in `mesh_3d` depended on how the lattice was divided into blocks, that is on the thread count: the cells were appended block by block as (whole tetrahedra, cut tetrahedra), so the order of the cells changed with the number of blocks, and equal cut ids made by two blocks with different last bits were deduplicated by an unstable sort of the id alone.

## Decision

* After the cells are made and before the facets are derived, every cell edge shorter than `1e-8` of the smallest lattice spacing welds its ends (union-find, representative = lowest id, flags merged); cells with a repeated vertex are dropped. Nothing of ordinary size is touched. With it the same box has a smallest cell volume of 5e-3 (1e-16 before), the structured layers build, and the sliver count of the benchmark models is zero.
* All whole tetrahedra come first, then all cut ones, each in lattice order, independent of the blocks; cut points of equal id are sorted by id and then by position before the deduplication.
* Test: `meshes_are_identical_for_any_thread_count` runs the probe `probe_mesh_hashes` in child processes with `NYX_THREADS` 1, 3 and 8 and compares the hash of points, cells, facets and labels for the simple, adaptive, structured and adaptive-structured plans of a ball in a box with an absorbing window. It failed (the cell order differed from index 19112 on) before the second fix.

## Consequences

Meshes of models without plane-aligned surfaces (the sphere cavity of the benchmark, the torus) are unchanged cell for cell; meshes with them lose their slivers and have a different but equivalent cell numbering. Nothing in the goldens uses a 3D geom mesh.

## Code

`crates/femcore/src/meshing/mesher3d.rs` (steps 4, 5 and 5b).

## Contract

`tests/mesh3d_features.rs` (`meshes_are_identical_for_any_thread_count`, the volume checks of every structured test).
