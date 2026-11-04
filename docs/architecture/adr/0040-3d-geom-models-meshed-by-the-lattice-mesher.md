# 40. A 3D geom model is meshed by the lattice mesher through `Domain3`

Date: 2026-11-04

## Status

Accepted (restructure phase 2, bookmark `restructure-2b-mesh3d-adapter`). Reverses the "test fixture only" half of R05 (ADR 31, ADR 32) for `mesher3d` and `refine3d`, on paul's request ("sure lets have a 3d mesh adapter for these", 2026-11-04). Uses ADR 18 (no curved elements) and ADR 32 (1-based labels, `CellMesh` as the native mesh).

## Context

`geom` evaluates 3D models (prismatic layer stacks, curved bodies inside a prismatic host, polyhedral models) with exact faces, but `meshing::geomdrive::mesh_geom_2d` refused every model that is not 2D, so no 3D geom model could be meshed. The isosurface-stuffing mesher (`mesher3d`) and its Rivara refinement (`refine3d`) were parked in `tests/common/solid3d` because they worked from their own `Solid` enum (box, sphere, cylinder), whose signed distance and face numbering duplicated what geom knows, and nothing in the library called them.

## Decision

* `meshing::mesher3d` and `meshing::refine3d` are library code again. The mesher knows no primitive: it asks the trait `Domain3` for the region at a point, the box of each region, the labels, and a projection of a new boundary vertex onto the exact face. `Solid::sdf` and `face_at` are gone.
* The surfaces that are fitted are the boundaries of the regions. A point is inside or outside a region (`Domain3::region_at`); where a lattice edge changes side the mesher bisects on that, so cut points are exact to bisection precision on a sphere, a torus or a plane alike. Where two surfaces coincide (an interface is the boundary of both of its regions, and the boundary of the solid can coincide with the one of a region) the points on it are classed from the two ends of the bisection bracket: every surface whose side changes across it counts as "on".
* The geom adapter is `meshing::geomdrive3::GeomDomain3`: `Evaluated::classify` for the region, `Complex::facet_project` for the exact closest point of a face, `LabelMap` for the integers. `mesh_geom_3d` and the quick problem `sphere_cavity` are its entry points. `meshing::mesh_evaluated_report` dispatches on `list.dim`, so `mesh_model` meshes 3D models with `Mesher::Simple` and inline absorbing layers (the adaptive mesher and the structured PML are 2D and are refused with a message).
* The mesher returns the native `mesh::cells::CellMesh<3>` (all tetrahedra, labelled cells and facets); `MeshedGeom.mesh` is `SimplexMesh::from_cells` of it, the solver view. `refine3d` works on `CellMesh<3>`.
* Labels (ADR 32): cells take the domain number of their region (overlay order from 1, pins honoured); boundary facets are the faces of cells that have the void on the other side and the faces between cells of different regions, and take the label of the geom face they lie on (the nearest face with exactly those two regions on its sides). Names and selectors (`MeshedGeom::boundary("air|ball")`, `select("plane(z = 0)")`) therefore mean the same as in 2D.
* Element sizes: the global `h`, lowered per region by a `Size` rule on regions and by the layers of an absorbing piece; a finer region is refined by Rivara bisection and its new boundary vertices are projected onto the exact face.
* A hole or void host contributes no cells, so a PEC cavity is a sphere inside a hole box (`sphere_cavity`); geom needed no change for that.

## Not done (scope)

* **No change in `geom`.** The prismatic restrictions stay as geom states them: curved bodies must sit strictly inside a prismatic host (touching or overlapping curved bodies are refused), a sphere or torus alone has no host (`stack.rs:203`), and the cylinders and extrusions of one model share a sweep axis (mixed axes need curved-surface intersections, `scene.rs:255`). The adapter meshes whatever geom evaluates. General 3D placement of curved solids is a geom feature (an arrangement of curved surfaces), not a meshing one; if paul wants it, it is a separate geom task and this adapter needs no change for it.
* **Creases are not fitted.** The old fixture fitted sharp edges and corners of boxes and cylinders that do not lie on a lattice plane. The adapter drops that: the lattice is aligned with the box of the model, so the creases of a box host and of lattice-aligned layers are exact; a slanted crease is cut by one face and warped by the other, with an error of order `alpha h`. Geom has the data for it (`Edge3`); it is a later feature.
* The adaptive mesher, error estimators and `mesh_for_problem` (wavelength sizing) are not extended to 3D. The 3D structured PML (`pml3`, `pml_box`) stays a separate generator.
* The new mesher is block-parallel with fixed blocks and deterministic cut ids, but a test that the mesh is identical for any `NYX_THREADS` (ADR 11) and a benchmark are still to add.

## Consequences

* The fixture `tests/common/solid3d` and `meshing/README.md` ("moved") are deleted; `emw3d` and `cavity_bench` build the sphere from geom.
* The old fixture's own 15 tests are not brought back; `tests/mesh3d.rs` tests the adapter, and `emw3d::spherical_cavity_lowest_mode_on_a_stuffed_mesh` checks the sphere end to end.
* Accuracy at the sizes tested: a sphere (radius 0.7, h 0.2) loses 1.9 % of its volume to the inscribed polyhedron; a torus (R 1, r 0.3, h 0.1) has volume error 0.9 % and area error 0.2 % against `2 pi^2 R r^2` and `4 pi^2 R r`; the PEC sphere cavity (R 0.05, h 0.012, order 2) stays within 1 % of the closed-form lowest mode.

## Code

`crates/femcore/src/meshing/{mesher3d,refine3d,geomdrive3}.rs`, dispatch in `meshing/mod.rs`.

## Contract

`mesh_geom3d_numbers_cells_and_facets_with_the_model_label_map` in `crates/femcore/tests/contracts.rs` (contract `mesh` 1.2, additive); behaviour in `crates/femcore/tests/mesh3d.rs` and `emw3d.rs`.
