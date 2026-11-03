# 6. One meshing interface: MeshPlan, with a simple and an adaptive mesher behind it

Date: 2026-11-02

## Status

Accepted (commit `a94dad3`, "Unified meshing interface").

## Context

Two ways of meshing a geom model existed: the constrained Delaunay mesher the app had always used (`geomdrive::mesh_geom_2d`) and the new metric adaptation of the mesher port. Callers (app, use cases, the MCP session) should not care which one made the cells, and names, selectors and boundary conditions must mean the same on either mesh.

## Decision

`meshing::mesh_model(list, &MeshPlan)` is the one way to mesh a geom model. `MeshPlan { opts, mesher }` chooses `Mesher::Simple` (chords of each geom edge, CDT, Ruppert refinement, lattice points in finer regions) or `Mesher::Adaptive(Adaptive { scale, grade, max_cells })` (the same start, then metric adaptation toward a size field built from the model: region sizes, `Size` rules, absorbing layers, edge curvature; graded). The result is a `MeshedGeom` with names, labels and exact geometry.

* Cell and facet labels are carried through unchanged by adaptation.
* Vertices that move or are created on a boundary or interface land on the **exact** geom edge (`GeomOracle2`), not on the chord polygon; vertices of the geometry never move.
* Adapting to an ungraded metric with sharp features is unreliable, so `grade` defaults to 0.4 and is always applied.
* `adapt_to_metrics` is the hook for an error estimator; a `MeshCache` gives warm starts and reuse across parameter sweeps (`154b64e`, tied to `geom::diff`).

## Consequences

* The app has a switch (`GeometrySpec::adapt_scale`) and a checkbox for physics sizing (`GeometrySpec::physics_sizing`, ADR 17).
* The 3D route (`physics_size::adapt_tets`, `pml3`, `estimate3`) has no app entry because there is no 3D geom mesh to adapt (extruded models do produce tets; the hook is not written).
* Mesh equality across the two meshers is only by names and exact geometry, not by cell count.

## Code

`crates/femcore/src/meshing.rs`, `meshing_pml.rs`, `geomdrive.rs`, `mesher.rs`, `meshers/adapt/{geom,field,drive,topo,ops}.rs`; `crates/app/src/state.rs`. History: `a94dad3`, `3fbb24d`, `154b64e`, `e69fefd`.

## Contract

`contracts/mesh` (labels survive adaptation) and `contracts/geom_model`.
