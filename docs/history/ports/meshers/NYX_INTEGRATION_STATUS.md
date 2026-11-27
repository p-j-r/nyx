# Mesher work on top of `meshing-unification` (jj bookmark `mesher-work`)

What the solvers actually consume (read from the code, 2026-11): the physics solvers take `SimplexMesh` (P1/P2, Nédélec orders 1 and 2 on
triangles and tetrahedra); PML is per-domain-label stretching from each domain's bounding box (`emw_pml::cell_stretches`). Quadrilaterals,
prisms and hexahedra reach them through `SimplexMesh::from_cells_split`. Solvers for native mixed cells, curved nodes and higher orders are
being written elsewhere; the mesh side below keeps the native `CellMesh` available (`pml3::PmlBox::cells`,
`meshing_pml::structured_pml_native`) and assumes nothing about their interface.

## Added

| Where | What | Checked against |
|---|---|---|
| `physics_size` | element size per domain from the wavelength of the study (EM modes / eigenfrequencies, solid eigenfrequencies with the slowest anisotropic plane wave, Brillouin), points per wavelength by element order (12 / 6), PML relaxation, eigenfrequency headroom | PEC guide TE10 index at orders 1 and 2 (7e-5, 5e-7), index-3 core gets ~2.7x finer edges |
| `meshing::mesh_for_problem` | 2D geom pipeline: physics metric intersected with the geometry metric | as above |
| `estimate::Estimator::{combine, e_field}` | E components as one vector field, one metric | dielectric core in PEC box: equal error per cell to intersected components and to `\|E\|` |
| `meshing_pml::structured_pml_native` | 2D structured layers as quads, uncut | counts, orientation, split equals `structured_pml` |
| `meshers::adapt::{topo3, geom3, field3, ops3, drive3}` | 3D adaptation: classification (interior / surface / crease / corner), facet-and-crease geometry oracle with sphere and SDF projectors, split, collapse (interior, surface, crease), 3-2 / 2-3 flips, edge removal (4-7 tets), flat-boundary flips, length-equalising smoothing | volume, boundary area, corners and interface exact; box 6144 tets 100 % in band; graded 5092 tets 98 % in band, worst tet 0.59 |
| `physics_size::adapt_tets` | physics-sized, graded 3D adaptation of a tet mesh | PEC box cavity: order 2, 1601 tets, 9e-5 on 4 modes |
| `pml3` | 3D structured PML (prisms on faces, hexahedra on edges and corners) for tet blocks, `Pml` per layer domain, materials copied to layer labels | real stretch 2 over a layer equals the longer cavity (1.8e-4); six unit-stretch layers equal the bigger cavity (2.3e-4); mesh convergence of damped modes (ignored study) |
| app `GeometrySpec::physics_sizing` | opt-in physics sizing for the problem a mesh is built for (batch and interactive rebuild); no GUI control yet | unit test |

## Not done / limits

* Lossy PML is checked for convergence with the mesh, not for reflection (eigen solvers give no driven response).
* 3D adaptation is serial; closest-point queries are linear in patch size; flat-boundary flips only (curved surfaces and interfaces keep their
  triangulation); edge removal stops at 7 tets; curved-surface collapse is limited to ~25 degrees of normal turn. No Floquet-periodic or
  mirror-constrained 3D adaptation (frozen labels exist, untested in 3D).
* No 3D geom pipeline exists, so `adapt_tets` starts from a tet mesh (block, extrusion, imported).
* Points per wavelength are starting values tested on smooth modes only.
* Not started: Maxwell residual estimator, goal-oriented adaptivity, field transfer, curved elements (waiting on the solver interface).
* `FEMMESH_CHECK=1` validates the 3D working mesh after every pass; `FEMMESH_VERBOSE=1` prints iterations.
