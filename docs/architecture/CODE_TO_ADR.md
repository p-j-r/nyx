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
| `docs/UX_STATUS.md` | 20, 88 | handoff for the UX review and plan |
| `docs/ux/*` | 88 | UX review, action plan, colour tokens, contrast script and screenshots of the window |
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
| `crates/femcore/src/linalg/*` | 33, 8, 11 | ported solver library; `driver.rs` is the door the kernels use, `csr.rs` the one matrix type, `amg_elastic/` the elasticity multigrid, `arnoldi_real.rs` the real waveguide eigensolver |
| `crates/femcore/src/mesh/mod.rs` | 4 | `SimplexMesh`, the solver contract (was `mesh.rs`) |
| `crates/femcore/src/mesh/view.rs` | 4 | draw-only mesh |
| `crates/femcore/src/mesh/cells.rs` | 5 | per-cell-type mesh (was `meshers/cellmesh.rs`) |
| `crates/femcore/src/mesh/exchange.rs` | 32, 41 | the mesh a physics is given (native cell mesh, 2D or 3D) and its numbering of unknowns |
| `crates/femcore/src/mesh/io/*` | 16 | `.fmsh` |
| `crates/femcore/src/mesh/{quality,fixtures}.rs` | 7 | quality statistics, test meshes |
| `crates/femcore/src/meshing/mod.rs` | 6, 7, 18, 50, 51, 52 | `MeshPlan`, `mesh_model` (was `meshing.rs`) |
| `crates/femcore/src/meshing/geomdrive.rs` | 2, 6 | |
| `crates/femcore/src/meshing/mesher.rs` | 6, 2 | 2D triangulation core |
| `crates/femcore/src/meshing/mesher3d.rs` | 40, 55 | lattice (isosurface stuffing) tetrahedral mesher on the trait `Domain3`; was the R05 fixture |
| `crates/femcore/src/meshing/refine3d.rs` | 40 | Rivara longest-edge refinement of a `CellMesh<3>` |
| `crates/femcore/src/meshing/geomdrive3.rs` | 40, 2, 50, 52 | adapter: a 3D geom model as `Domain3` (classify, facet_project, LabelMap), `mesh_geom_3d`, `sphere_cavity` |
| `crates/femcore/src/meshing/pml3.rs` | 5, 7 | |
| `crates/femcore/src/meshing/pml_geom3.rs` | 52, 5 | structured absorbing layers for a meshed 3D geom model: the physical window, `build_pml`, geom labels, the tetrahedral cut |
| `crates/femcore/src/meshing/adapt3d.rs` | 51, 40 | adaptive meshing of a 3D geom model: the exact faces as the oracle of `adapt3`, interface facets, vertex sizes |
| `crates/femcore/src/meshing/request.rs` | 54 | the serialisable mesh request and report, their JSON Schemas |
| `crates/femcore/src/meshing/pml_geom.rs` | 5, 7, 52 | structured PML for a meshed geom model, an adapter over `pml_box` (was `meshing_pml.rs`, then `meshing/pml.rs`) |
| `crates/femcore/src/meshing/pml_box.rs` | 7, 5 | structured PML generator (was `meshers/adapt/pml.rs`) |
| `crates/femcore/src/meshing/physics_size.rs` | 17 | |
| `crates/femcore/src/meshing/sizing.rs` | 17, 25 | the neutral Sizing a physics reports to the mesher (re-exported by nyx-physics) |
| `crates/femcore/src/meshing/validate/*` | 7 | |
| `crates/femcore/src/par.rs` | 11 | R03; the fixed-chunk reductions of the adaptation live here too (was `adapt/par.rs`) |
| `crates/femcore/src/adapt/estimate*.rs` | 7, 17, 53 | |
| `crates/femcore/src/adapt/work.rs` | 5, 4, 7 | the adaptation's private simplex working mesh and its conversion from and to the solver mesh (was `mesh/native.rs`, `mesh/interop.rs`) |
| `crates/femcore/src/adapt/maxwell.rs` | 53 | residual error estimator for edge elements and the equidistributing metric |
| `crates/femcore/src/adapt/{ops3,topo3,drive3}.rs` | 7, 51 | 3D adaptation operators, working mesh, driver; overlap check, edge-removal and flip guards, compaction (ADR 51) |
| `crates/femcore/src/adapt/*` | 7 | metric adaptation, timer, mod |
| `crates/femcore/src/material/*` | 1, 29, 30, 34, 35 | spec, units, elastic, tensor, optics, schema (phase 6) |
| `materials/*` | 29, 30, 34, 35 | the shipped material library (TOML) |
| `docs/MATERIALS.md` | 29, 30, 34, 35 | |
| `crates/femcore/src/expr.rs` | 2 | duplicate of geom::expr (R06) |
| `crates/femcore/src/validate.rs` | 9 | |
| `crates/femcore/src/progress.rs` | 9, 11 | |
| `crates/femcore/src/section.rs` | 1 | |
| `crates/femcore/src/props.rs` | 41, 25 | shared per-domain data (materials, PML, periodic pairs, linear solver); replaced by material files in phase 6 |
| `crates/femcore/src/testkit/*` | 41, 10, 60 | mesh fixtures of the tests of every physics crate; `soi.rs`, the cross-section of the SOI wire |
| `examples/*` | 60 | the example projects, written by `nyx-engine/tests/make_examples.rs` |
| `crates/femcore/src/periodic.rs` | 12, 9 | Floquet pairs; belongs to no single physics |
| `crates/femcore/src/edge.rs` | 45, 9 | the 2D edge-element basis shared by the ems kernels and the field sampler (split out of `ems.rs`) |
| `crates/femcore/src/pml.rs` | 45, 7 | coordinate stretching shared by the ems and sms kernels (split out of `ems_pml.rs`) |
| `crates/femcore/src/optical.rs` | 21, 9, 45 | optical mode sets, field sampling and the complex material tensors (split out of `sbs.rs` in phase 1; the solvers that fill a set are in `nyx-ems`) |
| `crates/femcore/src/fields.rs` | 9 | cell fields |
| `crates/femcore/src/fem.rs` | 4, 8 | element basics, original stack |
| `crates/femcore/src/lib.rs` | 12 | module map |
| `crates/femcore/tests/linalg_contract.rs` | 33, 19 | driver contract and closed-form checks |
| `crates/femcore/tests/contracts.rs` | 19, 1, 4, 5, 9, 16, 51, 52, 53, 54 | |
| `crates/femcore/tests/material*.rs` | 29, 30, 34, 35, 19 | behaviour and contract tests of the material format |
| `crates/femcore/tests/golden_frames*` | 10, 1 | |
| `crates/femcore/tests/frame*.rs` | 1 | |
| `crates/femcore/tests/ring_frame.rs` | 1 | |
| `crates/femcore/tests/linalg/*` | 8, 11 | |
| `crates/femcore/tests/linalg/main.rs` | 8 | |
| `crates/femcore/tests/meshers.rs` | 7, 5 | |
| `crates/femcore/tests/adapt3d.rs` | 7 | |
| `crates/femcore/tests/mesh3d_features.rs` | 50, 51, 52, 55 | wavelength sizing, adaptive meshing, structured layers and thread-count independence for 3D geom models |
| `crates/nyx-ems/tests/maxwell_estimator.rs` | 53 | the Maxwell estimator and loop on a PEC cavity and on a sphere from geom |
| `crates/femcore/tests/meshing.rs` | 6 | |
| `crates/femcore/tests/geom_meshes.rs` | 2, 6 | |
| `crates/femcore/tests/*` | 9, 12 | physics tests: solid, solid_pml, anisotropic, common |
| `crates/femcore/examples/fs_*.rs` | 8 | ported CLI tools |
| `crates/femcore/examples/mesh_*.rs` | 7, 16 | mesher studies |
| `crates/femcore/examples/*` | none | benchmarks and studies; the numbers are in the history |
| `crates/femcore/Cargo.toml` | 13, 8, 11 | features `hdf5-metno`, `phdf5`, `portable-kernels` |
| `crates/nyx-store/*` | 14, 39 | JSON specs and the card function (39) |
| `crates/app/src/{tokens,theme,kit,menu,regions,activity,sweepui,cvd}.rs` | 88 | colour tokens as data, the looks built from them, the widgets, the one table of actions (palette and native menu), the colour slots of regions (solid, hatched, grey), the busy chip, queue popover and stale band from the engine queue, the Sweep section of the Build page, the colour-vision simulation of the scripted shots |
| `crates/app/tests/colour_literals.rs` | 88 | the ratchet: colour literals in the window code only go down |
| `crates/app/src/{shell,wfview,sheets,dashboard}.rs` | 14, 23 | |
| `crates/app/src/*` | 23 | |
| `crates/app/*` | 23 | |
| `crates/nyx-mcp/tests/contracts.rs` | 19, 15 | |
| `crates/nyx-mcp/*` | 15 | moved from `ai/` into the root workspace (ADR 26) |
| `crates/nyx-connect/*` | 15, 28, 19 | one copy of the Connect logic (R07 resolved)|
| `tools/*` | 20 | baseline and helper scripts |
| `.cargo/*` | 26, 27 | shared target directory |
| `crates/nyx-physics/src/decode.rs` | 41, 25 | typed parameters and name resolution shared by every physics |
| `crates/nyx-physics/src/testkit.rs` | 41, 10 | the golden-file harness |
| `crates/nyx-physics/src/frame.rs` | 56, 25 | the frame a result states for the analyses of its vector fields (axes, handedness, absorbing domains) |
| `crates/nyx-physics/*` | 25, 19 | the Physics trait, envelope, descriptors, registry |
| `crates/nyx-dfs/tests/golden/*` | 10, 41 | the dfs goldens, moved unchanged from femcore |
| `crates/nyx-dfs/*` | 41, 25, 12, 9 | diffusion: typed problem, electrostatics and heat kernels (moved from femcore::dfs), descriptors |
| `crates/nyx-sms/tests/golden/*` | 10, 42 | the sms goldens (and `ring_model_axes`), moved unchanged from femcore |
| `crates/nyx-sms/src/{axisym,cross_section}.rs` | 43, 42 | bodies of revolution on their half-plane (library entry, not in the envelope) |
| `crates/nyx-sms/*` | 42, 41, 25, 12, 9 | solid mechanics: typed problem, static, eigenfrequency, absorbing-layer and Floquet kernels (moved from femcore::sms*), descriptors |
| `crates/nyx-ems/tests/golden/*` | 10, 45 | the ems goldens, moved unchanged from femcore (the two SOI wire files at the G2 tolerances) |
| `crates/nyx-ems/tests/axisym.rs` | 46, 45 | the axisymmetric kernel against closed forms (gate G1) |
| `crates/nyx-ems/tests/{emw2d,emw3d,pml,aniso_optics,bragg_stack,photonics,physics_meshing,meshing,moved}.rs` | 45, 1, 12 | the ems tests moved from femcore |
| `crates/nyx-ems/*` | 45, 25, 12, 9, 46 | electromagnetic waves: typed problem, mode, complex-mode, axisymmetric and cavity kernels (moved from femcore::ems*), descriptors |
| `crates/nyx-analysis/*` | 56, 25, 21 | the polarization analysis (Stokes, Poincare, classes, groups) as a function of a result; moved from `femcore::{analysis, polarization}` in phase 5 |
| `crates/nyx-engine/tests/layering.rs` | 26 | the layering table |
| `crates/nyx-engine/tests/jobs.rs` | 85 | every long operation is a job: cancel tests per kind, the dropped job of a replaced project, the window never blocks |
| `crates/nyx-engine/tests/queue.rs` | 87 | the queue: lanes, priority, round robin, ageing, supersede, hold, shared stages, order never changes a number |
| `crates/nyx-engine/tests/commandlist.rs` | 89 | command lists: static check against the catalog schemas, the dry-run plan, a bad line named and nothing run, the P0 workflows from their files |
| `crates/nyx-engine/tests/api_docs.rs` | 89, 57 | the API reference: `docs/api/` is what the catalog generates, every command has a page and an index row, every example passes `nyx check`, every argument is described |
| `docs/api/*` | 89, 57 | the API reference, one page per engine command, generated from `Engine::catalog()` (never edited by hand) |
| `crates/app/tests/cli.rs` | 89 | `nyx check` and `nyx run`: the plan, exit codes, a bad list runs nothing |
| `crates/nyx-engine/tests/sweep.rs` | 87 | sweeps: a point is an edit list, one recorded run per point, the numbers of a hand-made solve, the engine's own project and results untouched, a bad point stops the request |
| `crates/nyx-engine/tests/worker.rs` | 86 | a worker gives the numbers of the same solve in this process to the bit (five examples, a coupled chain), its stage log is the job log, a cancel kills it at once, what cannot be isolated says so, sweep points in workers |
| `crates/nyx-engine/tests/admission.rs` | 86 | estimates from worker runs, no estimate without two sizes, a solve waits for memory and says why, refusal unless forced, a cancel ends the wait, the history survives in a file |
| `crates/nyx-engine/tests/worker_crash.rs` | 86 | a worker that dies fails its job only and the next job runs |
| `crates/nyx-engine/tests/stale.rs` | 84, 64 | staleness per physics, late jobs, the `stale` fields of the replies |
| `crates/nyx-engine/tests/format.rs` | 63 | the `format` key at project, feature list and problem level; old and foreign files refused |
| `crates/femcore/tests/expr_grammar.rs` | 65 | the two expression grammars (`geom::expr`, `femcore::expr`) accept and refuse the same corpus (R06) |
| `crates/nyx-engine/tests/golden/*` | 10, 42 | `solid_static_thermal_from_heat`, moved unchanged from femcore |
| `crates/nyx-engine/examples/*` | 41, 9 | benchmarks and logs that drive every physics |
| `crates/nyx-engine/*` | 37, 38, 25, 26 | `Project`, commands, jobs, store wiring and geometry descriptors (phase 7, early); the registry from cargo features |
<!-- map:end -->

## Gaps this index shows

* Files with no ADR of their own: benchmarks and studies in `examples/`. They are measurements, not decisions.
* Modules whose ADR is retrospective and thin on rationale because the history gives only the subject line: `fem.rs`, `mesher3d.rs`, `refine3d.rs` (ADR 4, 7), `sms*.rs`, `ems*.rs`, `sbs*.rs` as physics (ADR 9, 12 are about structure and naming, not about the choice of formulation; `docs/HOW_IT_WORKS.md` is the explanation of the formulations).
* Decisions with no ADR because they are not recorded anywhere I could read: choice of `faer`, `spade`, `eframe`; choice of P1/P2 and Nedelec 1/2 as the element families (ADR 4 only states what exists); the Brillouin gain formulation. Ask paul if they should be written down.
