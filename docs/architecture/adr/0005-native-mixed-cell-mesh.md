# 5. Meshes stay native: one cell type per cell, simplex split only as a fallback

Date: 2026-11-02

## Status

Accepted as a constraint on the mesher side (paul, 2026-11-02, quoted in `docs/ports/meshers/FEATURE_DEPENDENCIES.md`). The solver side that consumes native cells is being written elsewhere and does not exist in this repository.

## Context

The structured PML layers are quadrilaterals (2D) and prisms and hexahedra (3D). Cutting them into simplices doubles the cell count (2D) or multiplies it by 3 to 6 (3D) and gives different finite-element spaces than the native cells. The solvers are being upgraded, so nothing should assume simplex-only or low order.

## Decision

* Meshes are native `CellMesh` values with a **type per cell** (Tri, Quad, Tet, Prism, Hex). `pml3::PmlBox::cells` and `meshing_pml::structured_pml_native` hold the uncut cells in nyx's tensor-product node order (quad `0,1,3,2` around the cell).
* The simplex split (`SimplexMesh::from_cells_split`) is a **fallback** that lets today's simplex-only solvers run; it is conforming because it depends only on global vertex ids.
* The mesher keeps its own dimension-generic types (`meshers::cellmesh::CellMesh<D>`, `meshers::mesh::SimplexMesh<D>`) as working types; `meshers::interop` translates at the border, including node order, and facets become the elements of one dimension lower with 0-based numbers.
* The mesher assumes no maximum element order and the Nedelec families; it states its node-order conventions in `meshers/cellmesh.rs` (Tri/Quad counter-clockwise, Tet right-handed, Hex bottom quad then top, Prism base then top with the handedness opposite to the VTK wedge).

## Consequences

* Two `CellMesh` types coexist (R01) and the label bases differ between them (R09). One of them has to go when the solver interface is fixed; the open questions are on the solver side (cell types and node orders it wants, whether boundary quads and triangles are both expected).
* The 3D PML box (`pml3`) and `meshers::adapt::pml` generate the same structure twice (R04).
* Nothing consumes native cells today; they are checked only by counts, orientation and split equivalence.

## Code

`crates/femcore/src/cellmesh.rs`, `meshers/cellmesh.rs`, `meshers/interop.rs`, `meshers/adapt/pml.rs`, `pml3.rs`, `meshing_pml.rs`. History: `b9a135e`, `154b64e`, `79936d6`.

## Contract

`contracts/mesh`: pins the node-order permutations of `interop` and the label bases.
