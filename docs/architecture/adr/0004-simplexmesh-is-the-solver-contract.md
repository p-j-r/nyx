# 4. SimplexMesh is the one mesh the physics kernels read

Date: 2026-10-28

## Status

Accepted for today; **superseded in direction by ADR 5** (native mixed cells), which is not implemented on the solver side yet.

## Context

The kernels (`dfs`, `sms`, `ems`, `sbs`) are finite-element codes on P1 and P2 Lagrange elements and Nedelec elements of order 1 and 2 on triangles (2D) and tetrahedra (3D). They need numbering of unknowns (vertices, then edge midpoints for P2), per-cell domain numbers and per-facet boundary numbers.

## Decision

`femcore::mesh::SimplexMesh` is the single mesh type the solvers take: runtime `dim` (2 or 3), flat `points`, `cells` (`dim + 1` vertices each), 1-based `cell_labels`, `facets` (edges in 2D, triangles in 3D) with 1-based `facet_labels`. Edges and the per-cell edge table are built on construction (`SimplexMesh::new`); P2 degrees of freedom number vertices first, then edge midpoints. Meshes with quadrilaterals, prisms or hexahedra reach the solvers through `from_cells`/`from_cells_split`, which cut those cells into simplices without adding vertices (cone from the smallest vertex id, quad diagonal through the smallest id, so neighbouring cells agree).

`solve(mesh, spec, state)` takes `&SimplexMesh`; `validate::check_inputs` refuses what it cannot solve (order other than 1 or 2, NaN geometry).

## Consequences

* Every mesher in the repository ends in a `SimplexMesh` (`geomdrive`, `meshing`, `mesher3d`, `meshers::interop::simplex_from_adapt`, the PML builders through the split).
* Solution vectors are indexed by dof in this numbering; a result is only meaningful with the mesh that produced it. `Solution` stores no mesh.
* Quads, prisms and hexes are never solved natively; the split changes the discrete space (Q1 versus P1, Nedelec on quads versus on triangles), so native-versus-split agreement can only be to discretisation error (`ports/meshers/STATUS.md`, M2).
* `SimplexMesh` and the mesher's `meshers::mesh::SimplexMesh<D>` are different types (rough spot R01).

## Code

`crates/femcore/src/mesh.rs`, `cellmesh.rs`, `solve.rs`, `validate.rs`, `view.rs`, `fem.rs`, `meshers/interop.rs`. History: `6de7458` (initial), `7a9a85e`.

## Contract

`contracts/mesh` (`crates/femcore/tests/contracts.rs`, tests `mesh_*`): label bases, facet counts, orientation, split volume conservation, round trip through `interop`.
