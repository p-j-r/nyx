# 52. Structured absorbing layers in 3D for geom models

Date: 2026-11-05

## Status

Accepted (bookmark `restructure-4m-mesh3d-features`). Lifts "the structured PML is 2D" of ADR 40. Uses ADR 32 step 5 (one PML generator, `pml_box::build_pml`, which keeps the physical vertex and cell numbers), ADR 50, ADR 51.

## Context

`PmlMode::Structured` rebuilt the absorbing layers of a 2D model as quadrilaterals around the physical window (`pml_geom::structured_pml`). In 3D `pml3::add_pml_3d` existed for a box-shaped tetrahedral mesh, with no link to a geom model: no domain numbers of the pieces, no boundary labels of the outer wall, no mesher that makes the physical window box-shaped.

## Decision

`PmlMode::Structured` on a 3D model (`meshing::pml_geom3`):

1. **The physical window is meshed alone.** The absorbing pieces of the model are left out of the lattice (`Sizes3::skip`; `GeomDomain3` treats a face between a meshed and a skipped region like a face to the void). The lattice mesher aligns its lattice with the box of the regions it meshes (ADR 40), so the window planes, edges and corners are lattice planes and the physical boundary lies on the window box exactly, which `build_pml` demands. A mesh made the ordinary way (pieces inline) works too: `structured_pml` drops the cells of the pieces; the window planes are surfaces of the model, cut exactly even off the lattice, so the boundary is box-shaped (test `layers_can_replace_the_inline_pieces_of_a_mesh_made_the_ordinary_way`).
2. Optionally the physical mesh is adapted (ADR 51) before the layers are built (the flow of the mesher chart, item P3: adapt the physical mesh, regenerate the layers). Box faces, edges and corners survive adaptation exactly.
3. **`build_pml`** sweeps prisms on the faces, hexahedra on the edges and corners, with the number of layers and the thickness of the `Pml` feature (face pieces; a piece that ends inside the layer, or two pieces of one side that disagree, is an error with a message).
4. **Labels are geom's.** Layer cells take the domain number of the absorbing piece they fill (`air.pml_xmax`, edge and corner pieces included; found by stretch direction and parent material); the outer wall takes the number of the geom face it lies on (nearest outer face of the piece); facets of the physical window that face a layer are interior now and are gone.
5. **Two meshes come out.** `structured_pml_3d_native` returns the uncut `CellMesh<3>` (tetrahedra, prisms, hexahedra; triangular and quadrilateral boundary facets; `CellMesh::check` passes: boundary facets only). `structured_pml_3d` cuts it for the simplex solvers (`SimplexMesh::from_cells_split`), orients every tetrahedron to positive volume (the lattice mesher does; the cut does not), and writes the faces between cells of different domains as facets with the number of the geom face between those two domains, as the inline mesh has them (ADR 40).
* **Physical numbering is kept**: vertices and cells of the physical mesh come first, in the same ids, order, connectivity and labels, in the native and in the cut mesh (contract test, 3D counterpart of the `build_pml` test).
* Restriction: the physical window must be box-shaped (all its boundary on the box faces); a physical region that is not (its boundary meets the box faces and something else) is refused by `build_pml` with its message. A model whose absorbing pieces belong to no meshed region (a sphere cavity in a hole host, whose host is the only region with pieces) has no layers to build.

## Consequences

* Volumes of every domain agree with the inline mesh to 1e-9 of the box; labels equal; conformity (boundary facets = faces of the void plus faces between domains) holds.
* The physics side needs nothing new: `MeshedGeom::pml_materials` already keys the stretch by domain number in 3D.
* The adaptation of a mesh *with* layers is not supported (the layers are rebuilt from the physical mesh each time, as in 2D).

## Code

`crates/femcore/src/meshing/{pml_geom3,pml_geom,geomdrive3,mod}.rs`.

## Contract

`mesh_geom3d_structured_pml_keeps_the_physical_vertex_and_cell_numbers`, `mesh_geom3d_adaptive_and_structured_meshes_keep_the_label_map` (`tests/contracts.rs`); `tests/mesh3d_features.rs` (`structured_layers_*`).
