# Code to ADR index

Every source, test, example, data and doc file of the repository maps to the ADRs that explain why it is the way it is. Files are grouped by module, because one ADR per file would repeat itself: the unit of a decision is a module, and each ADR lists its files under "Code".

* First match wins (top to bottom). `*` matches across `/`.
* "none" means that no decision of its own is recorded; the file is covered by the module's neighbours, or it is a benchmark or study whose result is in the git history.
* `docs/architecture/check_code_to_adr.sh` fails when a file matches no row, so a new file cannot be added without deciding which ADR owns it.
* ADR titles are in `ADR_INDEX.md`.

<!-- map:begin -->
| Glob | ADRs | Note |
|---|---|---|
| `Cargo.toml` | 13, 20 | workspace, profile choices |
| `Cargo.lock` | 13 | |
| `rust-toolchain.toml` | 13 | pinned compiler |
| `CLAUDE.md` | 1, 10, 12, 13, 20 | conventions |
| `README.md` | 13 | |
| `.gitignore` | 13, 20 | `*.py` ignored on purpose |
| `docs/architecture/adr/*` | 19 | the ADRs themselves |
| `docs/architecture/contracts/*` | 19 | exemplars |
| `docs/architecture/*` | 19, 20 | model, rough spots, indexes |
| `docs/RESTRUCTURE_STATUS.md` | 20, 19 | handoff for the restructure plan |
| `docs/FRAME_UNIFICATION_STATUS.md` | 1 | |
| `docs/KERNEL_PORT_STATUS.md` | 1, 9, 10, 12 | |
| `docs/GEOMETRY_ARCHITECTURE.md` | 2, 3, 18 | historical |
| `docs/GEOM_HANDOFF*.md` | 2, 3, 18 | historical |
| `docs/geom_handoff.zip` | 2, 3 | |
| `docs/GEOMETRY_GUI_PLAN.md` | 2, 23 | |
| `docs/WORKFLOWS_DESIGN.md` | 14, 23 | |
| `docs/MCP_ASSISTANT_PLAN.md` | 15 | |
| `docs/SWEEP_CACHING_PROPOSAL.md` | 6 | proposal; historical |
| `docs/HOW_IT_WORKS.md` | 4, 9 | |
| `docs/chip_report.md` | 22 | |
| `docs/photonic_use_cases.md` | 21, 22 | |
| `docs/brand/*` | 23 | |
| `docs/mockups/*` | 14, 23 | |
| `docs/ports/README.md` | 5, 7, 8 | |
| `docs/ports/NOTICES.md` | 7, 8 | licences |
| `docs/ports/meshers/*` | 5, 7, 16, 17, 18 | |
| `docs/ports/solvers/*` | 8, 11 | |
| `crates/geom/src/model.rs` | 2 | feature list |
| `crates/geom/src/expr.rs` | 2, 3 | duplicate of femcore::expr (R06) |
| `crates/geom/src/print.rs` | 2 | "Copy as Rust" |
| `crates/geom/src/diff.rs` | 2, 6 | sweep change detection, mesh cache |
| `crates/geom/src/eval.rs` | 2, 6 | tiers, cache |
| `crates/geom/src/error.rs` | 3 | |
| `crates/geom/src/labels.rs` | 2 | LabelMap |
| `crates/geom/src/select.rs` | 2 | selectors |
| `crates/geom/src/sim.rs` | 2, 7 | PML, size rules, periodic pairs as data |
| `crates/geom/src/stack.rs` | 3, 1 | layer stack |
| `crates/geom/src/topo3.rs` | 3 | |
| `crates/geom/src/query.rs` | 3, 18 | exact queries, closest point |
| `crates/geom/src/curve.rs` | 3, 18 | |
| `crates/geom/src/tess.rs` | 3 | display only |
| `crates/geom/src/*` | 3 | arrangement, complex, planar, prim, scene, math, lib |
| `crates/geom/tests/contracts.rs` | 19, 3 | |
| `crates/geom/tests/goldens.rs` | 10, 2 | |
| `crates/geom/tests/wafer.rs` | 1 | |
| `crates/geom/tests/sweep_axes.rs` | 1 | |
| `crates/geom/tests/selectors.rs` | 2 | |
| `crates/geom/tests/*` | 3 | |
| `crates/geom/golden/*` | 10, 2 | |
| `crates/geom/examples/*` | none | benchmarks and studies |
| `crates/geom/*.md` | 2, 3 | |
| `crates/geom/LICENSE` | 3 | |
| `crates/geom/rustfmt.toml` | 13 | only formatted crate (R23) |
| `crates/geom/Cargo.toml` | 3, 13 | serde only |
| `crates/femcore/src/solvers/*` | 8, 11 | ported solver library |
| `crates/femcore/src/mesh/mod.rs` | 4 | `SimplexMesh`, the solver contract (was `mesh.rs`) |
| `crates/femcore/src/mesh/view.rs` | 4 | draw-only mesh |
| `crates/femcore/src/mesh/cellmesh.rs` | 5, 4 | the nyx `CellMesh`, to be deleted in phase 2 |
| `crates/femcore/src/mesh/native.rs` | 5, 4 | dimension-generic native mesh (was `meshers/mesh.rs`) |
| `crates/femcore/src/mesh/cells.rs` | 5 | per-cell-type mesh (was `meshers/cellmesh.rs`) |
| `crates/femcore/src/mesh/interop.rs` | 5, 4 | border between the two mesh families, to be deleted in phase 2 |
| `crates/femcore/src/mesh/io/*` | 16 | `.fmsh` |
| `crates/femcore/src/mesh/{quality,fixtures}.rs` | 7 | quality statistics, test meshes |
| `crates/femcore/src/meshing/README.md` | 4, 7 | where mesher3d and refine3d went (R05) |
| `crates/femcore/src/meshing/mod.rs` | 6, 7, 18 | `MeshPlan`, `mesh_model` (was `meshing.rs`) |
| `crates/femcore/src/meshing/geomdrive.rs` | 2, 6 | |
| `crates/femcore/src/meshing/mesher.rs` | 6, 2 | 2D triangulation core |
| `crates/femcore/tests/common/solid3d/*` | 4, 7 | R05: stuffed solid meshes (sphere cavity), test fixture only, was `meshing/mesher3d.rs`, `refine3d.rs` |
| `crates/femcore/src/meshing/pml3.rs` | 5, 7 | |
| `crates/femcore/src/meshing/pml.rs` | 5, 7 | structured PML for the nyx mesh (was `meshing_pml.rs`) |
| `crates/femcore/src/meshing/pml_box.rs` | 7, 5 | structured PML generator (was `meshers/adapt/pml.rs`) |
| `crates/femcore/src/meshing/physics_size.rs` | 17 | |
| `crates/femcore/src/meshing/validate/*` | 7 | |
| `crates/femcore/src/adapt/par.rs` | 11 | R03 |
| `crates/femcore/src/adapt/estimate*.rs` | 7, 17 | |
| `crates/femcore/src/adapt/*` | 7 | metric adaptation, timer, mod |
| `crates/femcore/src/problem.rs` | 1, 9, 19 | |
| `crates/femcore/src/material.rs` | 1 | |
| `crates/femcore/src/expr.rs` | 2 | duplicate of geom::expr (R06) |
| `crates/femcore/src/validate.rs` | 9 | |
| `crates/femcore/src/progress.rs` | 9, 11 | |
| `crates/femcore/src/par.rs` | 11 | |
| `crates/femcore/src/solve.rs` | 9, 1 | |
| `crates/femcore/src/section.rs` | 1 | |
| `crates/femcore/src/dfs.rs` | 9, 12 | |
| `crates/femcore/src/sms*.rs` | 1, 9, 12 | |
| `crates/femcore/src/periodic.rs` | 12, 9 | Floquet pairs; belongs to no single physics |
| `crates/femcore/src/ems*.rs` | 1, 9, 12 | |
| `crates/femcore/src/optical.rs` | 21, 9 | optical mode sets and field sampling (split out of `sbs.rs` in phase 1) |
| `crates/femcore/src/usecases/*` | 22, 2 | |
| `crates/femcore/src/analysis.rs` | 21 | |
| `crates/femcore/src/polarization.rs` | 21 | |
| `crates/femcore/src/fields.rs` | 9 | cell fields |
| `crates/femcore/src/fem.rs` | 4, 8 | element basics, original stack |
| `crates/femcore/src/sparse.rs` | 8, 11 | original stack (R02) |
| `crates/femcore/src/csparse.rs` | 8 | original stack |
| `crates/femcore/src/linalg.rs` | 8 | original stack |
| `crates/femcore/src/amg.rs` | 8 | original stack |
| `crates/femcore/src/eigen.rs` | 8 | original stack |
| `crates/femcore/src/lobpcg.rs` | 8 | original stack |
| `crates/femcore/src/xmlio.rs` | 13, 19 | |
| `crates/femcore/src/h5io.rs` | 13, 16 | |
| `crates/femcore/src/lib.rs` | 12 | module map |
| `crates/femcore/tests/contracts.rs` | 19, 1, 4, 5, 9, 16 | |
| `crates/femcore/tests/golden_frames*` | 10, 1 | |
| `crates/femcore/tests/stages.rs` | 9 | |
| `crates/femcore/tests/frame*.rs` | 1 | |
| `crates/femcore/tests/wafer_names.rs` | 1 | |
| `crates/femcore/tests/ring_frame.rs` | 1 | |
| `crates/femcore/tests/determinism.rs` | 11 | |
| `crates/femcore/tests/solvers/*` | 8, 11 | |
| `crates/femcore/tests/solvers.rs` | 8 | |
| `crates/femcore/tests/meshers.rs` | 7, 5 | |
| `crates/femcore/tests/adapt3d.rs` | 7 | |
| `crates/femcore/tests/meshing.rs` | 6 | |
| `crates/femcore/tests/physics_meshing.rs` | 17, 7 | |
| `crates/femcore/tests/geom_meshes.rs` | 2, 6 | |
| `crates/femcore/tests/split.rs` | 4, 5 | |
| `crates/femcore/tests/finite_guard.rs` | 9 | |
| `crates/femcore/tests/cancel_and_log.rs` | 9 | |
| `crates/femcore/tests/usecases.rs` | 22 | |
| `crates/femcore/tests/photonics.rs` | 21, 22 | |
| `crates/femcore/tests/bragg_stack.rs` | 22 | |
| `crates/femcore/tests/emw2d.rs` | 1, 12 | old prefix (R13) |
| `crates/femcore/tests/emw3d.rs` | 1, 12 | old prefix (R13) |
| `crates/femcore/tests/pml.rs` | 7 | |
| `crates/femcore/tests/*` | 9, 12 | physics tests: scalar, solid, solid_kz, solid_pml, anisotropic, aniso_optics, amg_solve, common |
| `crates/femcore/examples/fs_*.rs` | 8 | ported CLI tools |
| `crates/femcore/examples/mesh_*.rs` | 7, 16 | mesher studies |
| `crates/femcore/examples/chip_*.rs` | 22 | |
| `crates/femcore/examples/{soi_modes,bragg_gratings,rib_sweep,rib_shapes,wg800}.rs` | 22 | |
| `crates/femcore/examples/*` | none | benchmarks and studies; the numbers are in the history |
| `crates/femcore/Cargo.toml` | 13, 8, 11 | features `hdf5-metno`, `phdf5`, `portable-kernels` |
| `crates/nyx-store/*` | 14 | |
| `crates/app/src/connect.rs` | 15 | duplicate of crates/nyx-connect (R07) |
| `crates/app/src/{shell,wfview,sheets,dashboard}.rs` | 14, 23 | |
| `crates/app/src/*` | 23 | |
| `crates/app/*` | 23 | |
| `crates/nyx-mcp/tests/contracts.rs` | 19, 15 | |
| `crates/nyx-mcp/*` | 15 | moved from `ai/` into the root workspace (ADR 26) |
| `crates/nyx-session/*` | 15 | dissolved in phase 7 |
| `crates/nyx-connect/*` | 15 | |
| `tools/*` | 20 | baseline and helper scripts |
| `.cargo/*` | 26, 27 | shared target directory |
| `crates/nyx-physics/*` | 25, 19 | the Physics trait, envelope, descriptors, registry |
| `crates/nyx-dfs/*` | 25, 12 | skeleton: descriptors; kernels move in phase 4 |
| `crates/nyx-sms/*` | 25, 12 | skeleton: descriptors; kernels move in phase 4 |
| `crates/nyx-ems/*` | 25, 12 | skeleton: descriptors; kernels move in phase 4 |
| `crates/nyx-analysis/*` | 25, 21 | skeleton: the Analysis trait; content moves in phase 5 |
| `crates/nyx-engine/tests/layering.rs` | 26 | the layering table |
| `crates/nyx-engine/*` | 25, 23, 26 | skeleton: registry from cargo features; Project and commands in phase 7 |
<!-- map:end -->

## Gaps this index shows

* Files with no ADR of their own: benchmarks and studies in `examples/`. They are measurements, not decisions.
* Modules whose ADR is retrospective and thin on rationale because the history gives only the subject line: `fem.rs`, `sparse.rs`, `csparse.rs`, `linalg.rs`, `amg.rs`, `eigen.rs`, `lobpcg.rs` (ADR 8 covers them as the "original stack"), `mesher3d.rs`, `refine3d.rs` (ADR 4, 7), `dfs.rs`, `sms*.rs`, `ems*.rs`, `sbs*.rs` as physics (ADR 9, 12 are about structure and naming, not about the choice of formulation; `docs/HOW_IT_WORKS.md` is the explanation of the formulations).
* Decisions with no ADR because they are not recorded anywhere I could read: choice of `faer`, `spade`, `eframe`; choice of P1/P2 and Nedelec 1/2 as the element families (ADR 4 only states what exists); the Brillouin gain formulation. Ask paul if they should be written down.
