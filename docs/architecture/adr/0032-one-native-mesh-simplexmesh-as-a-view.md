# 32. One native mesh; `SimplexMesh` becomes the solver view of it

Date: 2026-11-04

## Status

Proposed (phase 2, plan written, **not implemented**). Carries out ADR 5; resolves R01, R09 and R04 when done. Answers from paul (2026-11-04): 1-based labels everywhere, `interop` and the nyx `CellMesh`/`ElementType` deleted, one PML generator.

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

1. Labels: 1-based boundary numbers in `pml.rs`, the view reader and `from_cells` (delete the shift in `interop`); rewrite the `mesh_*` contract on purpose (ADR 19, `VERSIONS.md`).
2. `pml.rs`/`pml3.rs` return `cells::CellMesh`; `SimplexMesh::from_cells_split*` takes it; rebuild the hex/prism/quad cases of `tests/split.rs` on it.
3. Adaptation entry points take `CellMesh` (simplex only); the `simplex_to_adapt`/`simplex_from_adapt` call sites in `estimate*`, `meshing`, `physics_size` and tests go away.
4. Delete `mesh::cellmesh`, `mesh::interop`, `mesh::native::SimplexMesh`.
5. One PML generator; then `Sizing`; then rewrite the remaining `mesh_*` contracts and schemas and re-measure the benchmarks (best of three, idle machine).

## Consequences

Risk is in step 3 (adaptation speed and determinism; ADR 11) and in the numbering change of step 1 (anything that selects boundaries by number: `geomdrive` label maps, the app's selectors). Both are covered by the existing goldens and `tests/{meshers,adapt3d,meshing,physics_meshing,split}.rs`.

## Code

`crates/femcore/src/{mesh,meshing,adapt}/`.

## Contract

`mesh_*` in `crates/femcore/tests/contracts.rs` (rewritten at step 1 and 5).
