# 32. One native mesh; `SimplexMesh` becomes the solver view of it

Date: 2026-11-04

## Status

Accepted (phase 2). Steps 1 and 2 are done (bookmark `restructure-2-mesh`); steps 3 to 5 are in progress; the removal of `SimplexMesh` is a phase 4 deliverable (below). Carries out ADR 5; resolves R01, R09 and R04 when done. Answers from paul (2026-11-04): 1-based labels everywhere, `interop` and the nyx `CellMesh`/`ElementType` deleted, one PML generator.

## Context

Four mesh types exist today: `mesh::SimplexMesh` (the solver contract: P1/P2 numbering, edges, 1-based labels), `mesh::cellmesh::CellMesh` (nyx's list-of-element-types form, boundary numbers 0-based), and the ported `mesh::native::SimplexMesh<D>` and `mesh::cells::CellMesh<D>` (dimension-generic, used by adaptation and the structured PML generator). `mesh::interop` converts between them and shifts boundary numbers. Callers of the bridge today: `adapt::{estimate, estimate3}`, `meshing::{mod, physics_size, pml, pml3}`, tests `adapt3d`, `contracts`, `meshing`, `split`. 0-based boundary numbers survive in `cellmesh`, in the `.fmsh`/view reader and in `meshing::pml`.

## Decision (target)

* `mesh::cells::CellMesh<D>` is the one mesh: points, per-cell-type connectivity with a domain number per cell, boundary elements with a boundary number, **both numbered from 1**; nodes per cell are defined by `(cell type, order)` (curved-node storage stays out, ADR 18).
* `mesh::native::SimplexMesh<D>` disappears: the adaptation works on a `CellMesh` whose cell types are all simplices.
* `mesh::SimplexMesh` stays as the solver's view (`SolverMesh` view): `SimplexMesh::from_cells(&CellMesh<3|2>, SplitMode)` splits non-simplex cells, builds edges and the unknown numbering. Kernels keep using it unchanged until they take cells directly.
* `mesh::cellmesh` (nyx), `ElementType` and `mesh::interop` are deleted.
* `meshing::pml` (2D structured PML for the nyx mesh) is rebuilt on `meshing::pml_box` (one generator); the PML output is a `CellMesh`.
* `Sizing` (wavelength per region, order, headroom) is defined in `femcore::meshing`; `physics_size` takes it instead of `ProblemSpec`; `nyx_physics::Sizing` re-exports it.

## Steps (each its own commit, suite at baseline after each)

1. Labels: 1-based boundary numbers in `pml_geom.rs`, the view reader and `from_cells` (delete the shift in `interop`); rewrite the `mesh_*` contract on purpose (ADR 19, `VERSIONS.md`).
2. `pml.rs`/`pml3.rs` return `cells::CellMesh`; `SimplexMesh::from_cells_split*` takes it; rebuild the hex/prism/quad cases of `tests/split.rs` on it.
3. **Amended (paul and the coordinator, 2026-11-04):** the adaptation keeps a private working mesh. `mesh::native` moves to `adapt::work` and `mesh::interop` is dissolved (the two conversion functions live in `adapt::work`, `pml_domains` in `meshing::pml_box`), so the public mesh vocabulary is `mesh::cells::CellMesh` plus the `SimplexMesh` solver view; adaptation is not ported onto `CellMesh` (cost below).
4. `mesh::cellmesh` and `mesh::interop` are gone (step 2 and 3); `adapt::work::SimplexMesh<D>` stays internal to `adapt`.
5. One PML generator; then `Sizing`; then rewrite the remaining `mesh_*` contracts and schemas and re-measure the benchmarks (best of three, idle machine).

## Consequences

Risk is in step 3 (adaptation speed and determinism; ADR 11) and in the numbering change of step 1 (anything that selects boundaries by number: `geomdrive` label maps, the app's selectors). Both are covered by the existing goldens and `tests/{meshers,adapt3d,meshing,physics_meshing,split}.rs`.

## Code

`crates/femcore/src/{mesh,meshing,adapt}/`.

## Contract

`mesh_*` in `crates/femcore/tests/contracts.rs` (rewritten at step 1 and 5).

## Amendment 2026-11-04: what stays simplex, and why (audit for paul's question "is the simplex still useful?")

* **Solver kernels need simplices today.** `dfs` (P1/P2 on triangles and tetrahedra), `sms` (P1/P2 triangles and tets), `ems` (Nedelec order 1 and 2 on triangles with P1/P2 phi), `ems3d` (Nedelec order 1 and 2 on tetrahedra, 6 and 20 dofs), `ems_pml`, `sms_pml`, periodic pairs, fields, analysis and `h5io` all read `mesh::SimplexMesh`, with its edge list and unknown numbering (28 files of `femcore/src` outside the mesh trees, plus `app/src/{state,worker}.rs`, `nyx-session/src/tools.rs`, `nyx-physics/src/lib.rs`). The kernels are being upgraded and extended by other threads; none reads quads, prisms or hexes, so structured PML layers are cut into simplices (`SimplexMesh::from_cells_split`).
* **Adaptation is intrinsically simplex.** Edge split, collapse, flip, smoothing, 3D edge removal and face flips, Rivara bisection and the Hessian recovery (`adapt/ops.rs` 421 lines, `ops3.rs` 863, `topo.rs` 434, `topo3.rs` 467, `estimate*.rs` 631) only exist for triangles and tetrahedra and run on their own topology structures (`Work2`, `Work3`). `metric.rs` (1009 lines) is dimension-generic algebra and does not touch a mesh. The mesh type appears only at entry and exit of those structures (about 21 `SimplexMesh` references and 35 field accesses in `topo*.rs`, `work.rs`, `field*.rs`, `estimate*.rs`, `drive*.rs`).
* **Mixed cells never enter adaptation.** The pipeline is: adapt the physical simplex mesh, then build the structured layers (`meshing::pml_box::build_pml`, which keeps the physical vertex ids), then split for the solvers. Porting the entry and exit of adaptation onto `CellMesh<D>` restricted to simplices would be a day of mechanical work (about 60 lines of conversion code) with no gain, so the private working type stays; a later `Work2::from_cells(&CellMesh<2>)` is possible without touching the algorithms.

## Phase 4 deliverable: remove `SimplexMesh`

When the physics crates exist (phase 4), each physics reads the native `CellMesh` and numbers its own unknowns (edges, faces, nodes for its element family). Then `mesh::SimplexMesh`, `from_cells_split` and the split report go away; quadrilaterals, prisms and hexahedra are solved natively, which also removes the unknowns that the cut adds in PML layers. Until then `SimplexMesh` is the solver view of the one mesh.

## Step 5 (done): one PML generator

`meshing::pml_box::build_pml` (2D quadrilateral layers, 3D prisms and hexahedra) is the only generator of structured layers. `meshing::pml_geom` (`structured_pml`, `structured_pml_native`) and `meshing::pml3` (`add_pml_3d`) are adapters: they turn a geom model's `Pml` feature, or a box-shaped tetrahedral mesh, into a `PmlSpec`, call `build_pml`, and cut the result for the simplex solvers. The inline mode (geom's `Pml` pieces meshed as ordinary regions, `PmlMode::Inline`) is a property of the geom model, not a second generator. Guarantee, now a contract test (`mesh_pml_generator_keeps_the_physical_vertex_and_cell_numbers`): the output mesh starts with the physical vertices and cells unchanged in id, order, connectivity and label, and the layers are appended; this is what the mesher feature chart item P3 (adapt the physical mesh, regenerate the layers, map fields back by id) needs, and it holds in 2D and 3D.

## Step 6 (done): the neutral `Sizing`

`femcore::meshing::Sizing { order, headroom, wavelengths: region name -> shortest wavelength }` is defined by the mesher (`meshing/sizing.rs`) and re-exported as `nyx_physics::Sizing`; the serde shape is unchanged, so the `physics_descriptors` contract did not move. `physics_size` computes element sizes from a `Sizing` (`cell_sizes_for`, `vertex_sizes_for`, `metrics_for`, `adapt_tets_for`); a size is `wavelength / headroom / points-per-wavelength`, relaxed in PML domains. Names are resolved to domain labels by the caller (`Sizing::resolve`; `sizing::label_names` for a mesh with numbers only). The functions that take a `ProblemSpec` (`cell_sizes`, `vertex_sizes`, `metrics`, `metrics_d`, `adapt_tets`, `mesh_for_problem`) remain as the adapter for the old kernels through `spec_sizing`; they are deleted in phase 4 when each physics crate answers `Physics::sizing`. One float detail: the headroom of an eigenfrequency study used to multiply the frequency before the wavelength was formed; it now divides the wavelength, which can differ in the last bit.
