# 41. A physics crate owns its typed problem; the window reaches it through a legacy bridge until phase 7

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-physics`). Carries out ADR 25 for the kernels: `dfs` is the first physics moved (this ADR was written with it); `sms` and `ems` follow the same pattern. Amends ADR 32 (what "remove `SimplexMesh`" means) and ADR 9 (the stages now live in the physics crates).

## Context

ADR 25 created the `Physics` trait and the envelope, and left the kernels in `femcore`. They read `femcore::problem::ProblemSpec` (one struct for every physics, with enums for the studies and conditions), return `femcore::problem::SolveResult` (flat field names, a `disp` array, a `neff` pair, an `analysis`), and are called by `femcore::solve::solve` and, through it, by the window (`app`), the MCP session (`nyx-session`), about forty tests and some examples. The plan's definition of done for phase 4 is that the golden tests give the same numbers, that `solve()` is replaced by `registry.get(id).solve(...)`, and that nothing outside the physics crates names a physics. The window and the session are rewired only in phase 7 (ADR 37).

## Decision

1. **Typed problem per crate.** Each physics crate has `problem.rs` with its own typed `Problem` (study, conditions on mesh numbers, materials, options) and `decode(spec, ctx) -> Result<Problem, Error>` that reads the envelope: parameters through `nyx_physics::decode::Params` (descriptor defaults), region names and selectors through `ctx.resolver` (`materials_by_label`, `resolve_conditions`). The kernels take the typed problem; code that already has mesh numbers (tests, the window bridge) builds it directly. A condition that does not belong to the study is refused in `decode`.
2. **Shared per-domain data stays in femcore.** `femcore::props` holds `Materials` (scalars per domain, anisotropic tensors with a crystal orientation, PML per domain), `Pml`, `PmlDir`, `TensorKind`, `Periodic` and `LinearSolver`; the old `problem` module re-exports them. They are what phase 6 (material files) replaces. The envelope carries scalar properties only; anisotropic tensors reach a kernel through the typed problem until the material files exist. `fem` keeps the helpers that more than one physics uses (`check_boundary_labels`, `cell_values_at_nodes`, assembly, constraints, `solve_spd`).
3. **Three stages, public, in the crate** (ADR 9), plus `run(mesh, &Problem, &Arc<Progress>) -> Result<SolveResult, Error>`: finite-number checks of the mesh and the problem, the progress stage, the kernel, the finite check of the result (`nyx_physics::check_result_finite`), the stage log appended. `Physics::solve` is `validate_spec`, `validate`, `decode`, `run`.
4. **Results are the new `nyx_physics::SolveResult`**: serialisable, a vector field is one `Field` with `components` (a vector `E` with components `x`, `z`), per-cell fields are `CellField`. `Solution::flat_fields()` and `Solution::flat(name)` give the flat names the golden files and the window use (`|E|`, `Ex`, `Ez`, in that order). A coupling export (`dfs.T`) is the last solution's field `T`; the engine hands it to the next physics through `Ctx::inputs`, so `SolveState` is not part of the new path.
5. **Mesh (done).** Kernels name the mesh only as `nyx_physics::Mesh` (what the trait takes: `femcore::mesh::Mesh`, the native `CellMesh` in 2D or 3D, in `mesh/exchange.rs`) and `nyx_physics::Discretisation` (what a solve assembles on, from `nyx_physics::discretise(mesh)`; check a mesh with `nyx_physics::check_mesh`). *What "remove `SimplexMesh`" means* (ADR 32, phase 4 deliverable; reading agreed with the coordinator): the type has left the physics interface; `Mesh` is the one native cell mesh and `discretise` builds, for each solve, the numbering of unknowns (edges, P2 dofs) from it, **still cutting non-simplex cells into simplices exactly as `SimplexMesh::from_cells_split` does**, so every golden keeps its numbers (`Mesh::from_solver_mesh` and `discretise` round-trip a solver mesh with identical vertices, cells, facets, labels and edges; a unit test pins it). The struct itself stays in `femcore` as the numbering of unknowns (its name, `Discretisation`, is the alias kernels use) because the mesher and the window produce and hold it until phase 7; renaming it is mechanical then. Quadrilaterals, prisms and hexahedra solved natively would change the unknown count of PML layers and with it the goldens; that is a later, golden-changing step with its own ADR.
6. **Window bridge.** `nyx_engine::legacy` (`solve`, `solve_with`) takes the old `ProblemSpec` and `SolveState` and returns the old `SolveResult`: one file per physics converts the old problem to the typed one (same mesh numbers, no names), runs `run`, and converts the result back with `result_back`. A physics whose kernels are still in `femcore` is run by `femcore::solve`. `app` and `nyx-session` call the bridge (a temporary layering edge `nyx-session -> nyx-engine`). The bridge, `femcore::problem` and `femcore::solve` are deleted together in phase 7.
7. **Tests live in the crate they test.** Goldens move with their kernel to `crates/<crate>/tests/golden/` unchanged; the harness is `nyx_physics::testkit` (`golden`, `record`, `cell_fields_reproduce_the_plot_fields`), the mesh fixtures are `femcore::testkit` (`tests/common` re-exports it). A test that needs two physics (heat feeding thermal expansion) lives in `nyx-engine/tests` and uses the bridge. `femcore` keeps only tests of what stays in `femcore` (meshing, adaptation, linear algebra, materials).
8. **Adaptation loops are generic over the result**: `femcore::adapt::PlotFields` (plot fields per solution), implemented for the old result and, in `nyx-physics`, for the new one; `solve_adaptively` and `solve_adaptively3` take `R: PlotFields`.
9. **Benchmarks that drive several physics live in `nyx-engine/examples`** (`bench`, `stage_log`); `tools/baseline.sh` runs `bench` from there. The other four timings stay in `femcore` until their physics move.

## The pattern for the next physics (checklist)

1. Typed `Problem` and `decode` in `problem.rs`; descriptors refined (a physics-level `solver` option if it has one) and the exemplar `physics_<id>.v1.json` rewritten with a `VERSIONS.md` line.
2. Kernels copied, then edited only for types: `&Discretisation`, the typed problem, `Field::scalar/vector`, `cell_fields_of`. No change of formula or order of operations.
3. `run`, `Physics::solve`, `Physics::validate` (orders, finite mesh).
4. `nyx_engine::legacy/<id>.rs`: old to typed problem, result back (`disp` and `neff` for solid and optical results), and the arm in `legacy::solve_with`.
5. Tests: moved from the femcore test files with the same numbers and tolerances; goldens moved and passing without being rewritten; a test that the envelope and the typed problem give the same bits; a refusal test; stage-API contract test.
6. Delete the old module, its arm in `femcore::solve`, the `CODE_TO_ADR.md` rows; state in the status file which tests moved where.

## Consequences

* Positive: a physics is a crate with its own problem type; the window keeps working; every moved golden still passes unchanged; the engine's coupled chain and the window use the same kernel code.
* Negative: for a while two result shapes exist (the bridge converts); `femcore::problem::ProblemSpec` and `femcore::solve` survive until phase 7; the envelope cannot yet express anisotropic tensors or PML (phase 6 and the ems/sms descriptors).
* Open: the old adapters (`physics_size::spec_sizing`, `mesh_for_problem`) go when the window is rewired.

## Code

`crates/nyx-dfs/src/*`, `crates/nyx-physics/src/{decode,testkit}.rs` and the result helpers, `crates/femcore/src/{props,testkit}.rs`, `crates/nyx-engine/src/legacy/*`.

## Contract

`physics_dfs.v1.json` 1.1 (the `solver` option); `physics_api` (unchanged); the stage-API contract test in `crates/nyx-dfs/tests/contracts.rs`.
