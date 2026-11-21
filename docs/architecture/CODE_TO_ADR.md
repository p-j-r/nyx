# Code to ADR index

Every source, test, example, data and doc file of the repository maps to the ADRs that explain why it is the way it is. Files are grouped by module, because one ADR per file would repeat itself: the unit of a decision is a module, and each ADR lists its files under "Code".

* First match wins (top to bottom). `*` matches across `/`.
* "none" means that no decision of its own is recorded; the file is covered by the module's neighbours, or it is a benchmark or study whose result is in the git history.
* `docs/architecture/check_code_to_adr.sh` fails when a file matches no row, so a new file cannot be added without deciding which ADR owns it.
* ADR titles are in `ADR_INDEX.md`.

<!-- map:begin -->
| Glob | ADRs | Note |
|---|---|---|
| `Cargo.toml` | 13, 20, 107 | workspace, profile choices, pure-Rust `hdf5-pure` + `zlib-rs` |
| `Cargo.lock` | 13, 107 | |
| `rust-toolchain.toml` | 13 | pinned compiler |
| `CLAUDE.md` | 1, 10, 12, 13, 20 | conventions |
| `AGENTS.md` | 1, 10, 12, 13, 20 | agent-neutral repository guide |
| `README.md` | 13 | |
| `LICENSE` | 100 | the AGPL-3.0 text |
| `NOTICE` | 100 | the dual licence and the parts that keep their own |
| `.gitignore` | 13, 20 | `*.py` ignored on purpose |
| `docs/architecture/adr/*` | 19 | the ADRs themselves |
| `docs/architecture/contracts/*` | 19 | exemplars |
| `docs/architecture/*` | 19, 20 | model, rough spots, indexes |
| `docs/RESTRUCTURE_STATUS.md` | 20, 19 | handoff for the restructure plan |
| `docs/UX_STATUS.md` | 20, 88 | handoff for the UX review and plan |
| `docs/HOST_SETUP.md` | 97, 98, 94 | what costs time when an AI host drives nyx, and what to set in the host and do in the agent |
| `docs/LATENCY_STATUS.md` | 97, 20 | handoff for the latency work L1 to L3 |
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
| `docs/PLOTTING_PLAN.md` | 99 | design note of the mode figures (the script, the results layout, the stages) |
| `docs/GDSFACTORY_INTEGRATION.md` | 57, 94 | design note: a gdsfactory plugin that drives the engine; not built |
| `docs/MARKET_NOTE.md` | - | market note; no architecture decision |
| `docs/results_format_measure.py.txt` | 62 | the measurement script of the results format study (`.py.txt`: `*.py` is git-ignored) |
| `docs/RESULTS_FORMAT_STUDY.md` | 62, 99, 107 | study of the results file layout and pure-Rust `hdf5-pure` decision |
| `docs/figures/example_modes_gallery.png` | 99 | example output of the `figure` command |
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
| `crates/geom/src/poly_batch.rs` | 102, 2 | batch polygon cleanup and shared-corner-pinned RDP simplification |
| `crates/geom/src/*` | 3 | arrangement, complex, planar, prim, scene, synth, math, lib |
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
| `crates/femcore/src/mesh/mod.rs` | 4, 82 | `SimplexMesh`, the solver contract (was `mesh.rs`) |
| `crates/femcore/src/mesh/view.rs` | 4 | draw-only mesh |
| `crates/femcore/src/mesh/cells.rs` | 5, 82 | per-cell-type mesh (was `meshers/cellmesh.rs`) |
| `crates/femcore/src/mesh/exchange.rs` | 32, 41, 82 | the mesh a physics is given (native cell mesh, 2D or 3D) and its numbering of unknowns |
| `crates/femcore/src/mesh/io/*` | 16 | `.fmsh` |
| `crates/femcore/src/mesh/{quality,fixtures}.rs` | 7 | quality statistics, test meshes |
| `crates/femcore/src/meshing/mod.rs` | 6, 7, 18, 50, 51, 52, 82 | `MeshPlan`, `mesh_model` (was `meshing.rs`) |
| `crates/femcore/src/meshing/geomdrive.rs` | 2, 6, 82 | |
| `crates/femcore/src/meshing/mesher.rs` | 6, 2, 66 | 2D triangulation core; `extrude_y` with `Layer::fill` |
| `crates/femcore/src/meshing/mesher3d.rs` | 40, 55, 82, 83 | lattice (isosurface stuffing) tetrahedral mesher on the trait `Domain3`; was the R05 fixture |
| `crates/femcore/src/meshing/refine3d.rs` | 40 | Rivara longest-edge refinement of a `CellMesh<3>` |
| `crates/femcore/src/meshing/curve3.rs` | 82 | curved edges: mid-edge geometry nodes from the exact surface, Jacobian validity, `CurveReport` |
| `crates/femcore/src/meshing/geomdrive3.rs` | 40, 2, 50, 52, 82, 83 | adapter: a 3D geom model as `Domain3` (classify, facet_project, LabelMap), `mesh_geom_3d`, `sphere_cavity` |
| `crates/femcore/src/meshing/pml3.rs` | 5, 7 | |
| `crates/femcore/src/meshing/pml_geom3.rs` | 52, 5 | structured absorbing layers for a meshed 3D geom model: the physical window, `build_pml`, geom labels, the tetrahedral cut |
| `crates/femcore/src/meshing/adapt3d.rs` | 51, 40 | adaptive meshing of a 3D geom model: the exact faces as the oracle of `adapt3`, interface facets, vertex sizes |
| `crates/femcore/src/meshing/request.rs` | 54, 82, 83 | the serialisable mesh request and report, their JSON Schemas |
| `crates/femcore/src/meshing/pml_geom.rs` | 5, 7, 52 | structured PML for a meshed geom model, an adapter over `pml_box` (was `meshing_pml.rs`, then `meshing/pml.rs`) |
| `crates/femcore/src/meshing/pml_box.rs` | 7, 5 | structured PML generator (was `meshers/adapt/pml.rs`) |
| `crates/femcore/src/meshing/physics_size.rs` | 17 | |
| `crates/femcore/src/meshing/sizing.rs` | 17, 25 | the neutral Sizing a physics reports to the mesher (re-exported by nyx-physics) |
| `crates/femcore/src/meshing/validate/*` | 7 | |
| `crates/femcore/src/meshing/sheet3d.rs` | 102 | conformal 3D sheet-stack mesher (`Op::Sheet`, `<sheet>:edge`, Dompierre prism splitting, single-facet interior boundaries) |
| `crates/femcore/src/par.rs` | 11 | R03; the fixed-chunk reductions of the adaptation live here too (was `adapt/par.rs`) |
| `crates/femcore/src/adapt/estimate*.rs` | 7, 17, 53 | |
| `crates/femcore/src/adapt/work.rs` | 5, 4, 7 | the adaptation's private simplex working mesh and its conversion from and to the solver mesh (was `mesh/native.rs`, `mesh/interop.rs`) |
| `crates/femcore/src/adapt/maxwell.rs` | 53 | residual error estimator for edge elements and the equidistributing metric |
| `crates/femcore/src/adapt/{ops3,topo3,drive3}.rs` | 7, 51 | 3D adaptation operators, working mesh, driver; overlap check, edge-removal and flip guards, compaction (ADR 51) |
| `crates/femcore/src/adapt/*` | 7 | metric adaptation, timer, mod |
| `crates/femcore/src/material/*` | 1, 29, 30, 34, 35, 101, 104 | spec, units, elastic, tensor, optics, schema (phase 6), Zimmermann BCS conductivity (`bcs.rs`) |
| `materials/*` | 29, 30, 34, 35, 101, 104 | the shipped material library (TOML, including superconductors) |
| `docs/MATERIALS.md` | 29, 30, 34, 35, 101 | |
| `crates/femcore/src/expr.rs` | 2 | duplicate of geom::expr (R06) |
| `crates/femcore/src/fieldexpr.rs` | 90 | the typed expression tree over complex field values and its evaluator (N1; N2 added comparisons, `where` and more functions) |
| `crates/femcore/src/formula.rs` | 90 | the text formula parser that produces that tree: grammar, vectors, names resolved by the caller (N2) |
| `crates/femcore/src/validate.rs` | 9 | |
| `crates/femcore/src/progress.rs` | 9, 11 | |
| `crates/femcore/src/warm.rs` | 96 | warm starts of the eigensolvers: a context on the solving thread, the previous solve's modes as the start block, bytes for a worker process |
| `crates/femcore/tests/warm.rs` | 96 | same modes in fewer iterations; no context or an empty one is the same bits; another size is a cold start |
| `crates/nyx-engine/src/reuse.rs` | 96 | the chain of a warm sweep: the link, the cell a point leaves its modes in, the wait in the shared lane, the reuse record |
| `crates/nyx-engine/tests/warm_sweep.rs` | 96 | a warm sweep agrees with a cold one, later points start warm, a worker and the queue worker give one answer, a predecessor that cannot finish does not hold the lane |
| `crates/nyx-engine/tests/latency.rs` | 97 | a starred material sizes a mesh with absorbing layers, `mesh_build` answers `cost`, a convergence dry run and its limits, the noise rows of a lossless guide on a tilted crystal |
| `crates/femcore/src/section.rs` | 1 | |
| `crates/femcore/src/props.rs` | 41, 25 | shared per-domain data (materials, PML, periodic pairs, linear solver); replaced by material files in phase 6 |
| `crates/femcore/src/testkit/*` | 41, 10, 60 | mesh fixtures of the tests of every physics crate; `soi.rs`, the cross-section of the SOI wire |
| `examples/*` | 60 | the example projects, written by `nyx-engine/tests/make_examples.rs` |
| `crates/femcore/src/periodic.rs` | 12, 9 | Floquet pairs; belongs to no single physics |
| `crates/femcore/src/edge.rs` | 45, 9 | the 2D edge-element basis shared by the ems kernels and the field sampler (split out of `ems.rs`) |
| `crates/femcore/src/pml.rs` | 45, 7, 69 | coordinate stretching shared by the ems and sms kernels (split out of `ems_pml.rs`; `radial_coordinate` is the cylindrical stretch both axisymmetric kernels use) |
| `crates/femcore/src/optical.rs` | 21, 9, 45 | optical mode sets, field sampling and the complex material tensors (split out of `sbs.rs` in phase 1; the solvers that fill a set are in `nyx-ems`) |
| `crates/femcore/src/fields.rs` | 9 | cell fields |
| `crates/femcore/src/fem.rs` | 4, 8, 82 | element basics, original stack |
| `crates/femcore/src/lib.rs` | 12 | module map |
| `crates/femcore/tests/linalg_contract.rs` | 33, 19 | driver contract and closed-form checks |
| `crates/femcore/tests/contracts.rs` | 19, 1, 4, 5, 9, 16, 51, 52, 53, 54 | |
| `crates/femcore/tests/material*.rs` | 29, 30, 34, 35, 19, 101, 104 | behaviour, BCS, and contract tests of the material format |
| `crates/femcore/tests/sheet_mesh3d.rs` | 102 | tests of conformal 3D sheet-stack meshing and `<sheet>:edge` tags |
| `crates/femcore/tests/golden_frames*` | 10, 1 | |
| `crates/femcore/tests/frame*.rs` | 1 | |
| `crates/femcore/tests/ring_frame.rs` | 1 | |
| `crates/femcore/tests/linalg/*` | 8, 11 | |
| `crates/femcore/tests/linalg/main.rs` | 8 | |
| `crates/femcore/tests/meshers.rs` | 7, 5 | |
| `crates/femcore/tests/adapt3d.rs` | 7 | |
| `crates/femcore/tests/mesh3d_features.rs` | 50, 51, 52, 55 | wavelength sizing, adaptive meshing, structured layers and thread-count independence for 3D geom models |
| `crates/nyx-ems/tests/maxwell_estimator.rs` | 53 | the Maxwell estimator and loop on a PEC cavity and on a sphere from geom |
| `crates/femcore/tests/meshing.rs` | 6, 66 | |
| `crates/femcore/tests/curved.rs` | 82 | curved tetrahedra: the map, the quadrature, storage, sphere volume, thread independence |
| `crates/femcore/tests/face_size.rs` | 83 | face sizes and `max_dev` in the lattice mesher |
| `crates/nyx-ems/tests/curved.rs` | 82 | the cavity kernel on curved tetrahedra against closed forms and the PEC sphere cavity |
| `crates/nyx-engine/tests/curved.rs` | 82 | a periodic dielectric sphere with `curved` on and off |
| `crates/femcore/tests/geom_meshes.rs` | 2, 6 | |
| `crates/femcore/tests/*` | 9, 12 | physics tests: solid, solid_pml, anisotropic, common |
| `crates/femcore/examples/fs_*.rs` | 8 | ported CLI tools |
| `crates/femcore/examples/mesh_*.rs` | 7, 16 | mesher studies |
| `crates/femcore/examples/*` | none | benchmarks and studies; the numbers are in the history |
| `crates/femcore/Cargo.toml` | 13, 8, 11 | feature `portable-kernels` |
| `crates/nyx-storage/*` | 107, 62, 99 | pure-Rust HDF5 storage engine (`hdf5-pure` + `zlib-rs`), packed 5D Shuffle+LZF layout, streaming reader/writer, XDMF, sweep checkpoint |
| `crates/nyx-store/*` | 14, 39 | JSON specs and the card function (39) |
| `crates/app/src/{tokens,theme,kit,menu,regions,activity,sweepui,cvd}.rs` | 88 | colour tokens as data, the looks built from them, the widgets, the one table of actions (palette and native menu), the colour slots of regions (solid, hatched, grey), the busy chip, queue popover and stale band from the engine queue, the Sweep section of the Build page, the colour-vision simulation of the scripted shots |
| `crates/app/tests/colour_literals.rs` | 88 | the ratchet: colour literals in the window code only go down |
| `crates/app/src/{feedback,spans}.rs` | 92 | the problem and size lines of an edit, the polarization thread, and the timing spans of the frame log |
| `crates/app/src/{shell,wfview,sheets,dashboard}.rs` | 14, 23 | |
| `crates/app/src/*` | 23 | |
| `crates/app/*` | 23 | |
| `crates/nyx-mcp/tests/{contracts,stdio}.rs` | 94, 97, 19, 15 | the `mcp_tools.v3` snapshot (four tools), the handshake, and the gate of S7: the P0 workflows in script mode against the per-command baseline |
| `crates/nyx-mcp/*` | 94, 97, 15 | the MCP server: four script tools (`script_run`, `script_wait`, `script_check`, `api`) over one script kernel, no tool per command (S7; ADR 15 amended); moved from `ai/` into the root workspace (ADR 26) |
| `crates/nyx-connect/*` | 15, 28, 19, 98 | one copy of the Connect logic (R07 resolved)|
| `tools/*` | 20 | baseline and helper scripts |
| `.cargo/*` | 26, 27, 107 | shared target directory (no `HDF5_DIR` C dependency) |
| `crates/nyx-physics/src/decode.rs` | 41, 25 | typed parameters and name resolution shared by every physics |
| `crates/nyx-physics/src/testkit.rs` | 41, 10 | the golden-file harness |
| `crates/nyx-physics/src/frame.rs` | 56, 25 | the frame a result states for the analyses of its vector fields (axes, handedness, absorbing domains) |
| `crates/nyx-physics/src/cqed.rs` | 103, 104, 105, 106 | typed cQED payloads (`CapacitanceMatrixData`, `InductanceMatrixData`, `TlsMeta`, `CircuitRom`, `DrivenPortData`) |
| `crates/nyx-physics/*` | 25, 19 | the Physics trait, envelope, descriptors, registry |
| `crates/nyx-dfs/tests/golden/*` | 10, 41 | the dfs goldens, moved unchanged from femcore |
| `crates/nyx-dfs/src/magneto.rs` | 103 | gauge-stabilized London magnetostatics (`cross_section`, `planar_stream`, `mixed_3d`) and inductance matrix extraction |
| `crates/nyx-dfs/tests/capacitance_qubit.rs` | 103 | multi-terminal capacitance matrix & magnetostatic inductance benchmarks (Test 1) |
| `crates/nyx-dfs/*` | 41, 25, 12, 9, 103 | diffusion: typed problem, multi-RHS electrostatics, magnetostatics, and heat kernels, descriptors |
| `crates/nyx-sms/tests/golden/*` | 10, 42 | the sms goldens (and `ring_model_axes`), moved unchanged from femcore |
| `crates/nyx-sms/src/{axisym,cross_section,waveguide}.rs` | 43, 42, 68, 69, 80, 101 | bodies of revolution on their half-plane and 2.5D waveguide cross-sections (`exp(-j q y)`): the library entries (`solve_solid_axisym`, `solve_solid_waveguide`), the `axisymmetric_eigenfrequency` and `modes` studies, axis conditions, cylindrical and planar layers, the scalars that identify a mode (numeric frequency parts, energy shares per component and domain) |
| `crates/nyx-sms/tests/axisym_mode_scalars.rs` | 80 | the numeric frequency parts and the energy shares against closed forms whose family is known |
| `crates/nyx-sms/tests/{pml_floquet,waveguide_modes}.rs` | 101 | Floquet-periodic absorbing layers (`P^H K P` congruence, inactive mass skip) and native 2.5D elastic waveguide modes (`study = "modes"`, closed-form SH cutoff, 3D swept Floquet agreement, PML complex frequencies and shares) |
| `crates/nyx-engine/src/convergence.rs` | 61, 80, 97, 101 | mesh convergence on a scratch engine; solutions compared by index, paired by nearest eigenvalue, or one mode followed; lossless imaginary floor and quality-factor noise rows |
| `crates/nyx-engine/src/solve.rs` | 61, 80 | the solve chain in coupling order; puts the region names into the scalars a kernel names by domain number |
| `crates/nyx-sms/*` | 42, 41, 25, 12, 9, 68, 101 | solid mechanics: typed problem, static, eigenfrequency, axisymmetric eigenfrequency, waveguide modes, absorbing-layer and Floquet kernels (moved from femcore::sms*), descriptors |
| `crates/nyx-ems/tests/golden/*` | 10, 45 | the ems goldens, moved unchanged from femcore (the two SOI wire files at the G2 tolerances) |
| `crates/nyx-ems/tests/axisym.rs` | 46, 45 | the axisymmetric kernel against closed forms (gate G1) |
| `crates/nyx-ems/tests/{emw2d,emw3d,pml,aniso_optics,bragg_stack,photonics,physics_meshing,meshing,moved}.rs` | 45, 1, 12, 67 | the ems tests moved from femcore |
| `crates/nyx-ems/src/{ibc,rom}.rs` | 104 | superconducting thin-film IBC, rank-1 `lumped_jj`, physical $H$/$K$ normalization, passive `CircuitRom`, and `TlsMeta` |
| `crates/nyx-ems/tests/kicq_cpw_bcs.rs` | 104 | 3D CPW BCS kinetic inductance shift and lumped JJ equipartition tests (Test 3) |
| `crates/nyx-ems/src/{driven,ports}.rs` | 106 | driven frequency-domain solver, `lumped_port`, 2D boundary-mode `waveport`, and Hermite Tangential Galerkin PROM |
| `crates/nyx-ems/tests/driven_prom.rs` | 106 | 3D driven PROM, waveport, and BBQ admittance verification (Test 6) |
| `crates/nyx-ems/*` | 45, 25, 12, 9, 46, 67, 104, 106 | electromagnetic waves: typed problem, mode, complex-mode, axisymmetric, cavity, and driven kernels |
| `crates/nyx-analysis/src/capacitance.rs` | 103 | `capacitance_matrix` and `inductance_matrix` analyses of `dfs` results |
| `crates/nyx-analysis/src/{epr,hamiltonian,tls}.rs` | 105 | `epr` KICQ analysis, `hamiltonian` (`normal_fock` / `craig_bampton`), and `tls_participation` |
| `crates/nyx-analysis/tests/{tls_benchmark,park2026_transmon_benchmark,fluxonium_benchmark}.rs` | 105 | quantitative cQED acceptance benchmarks (Tests 2, 4, 5) |
| `crates/nyx-analysis/*` | 56, 25, 21, 99, 103, 105 | analyses of results (polarization, capacitance/inductance, EPR, circuit Hamiltonian, TLS, mode figures) |
| `crates/nyx-engine/tests/layering.rs` | 26, 107 | the layering table |
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
| `crates/nyx-engine/tests/preview.rs` | 92 | the problem of an edit, the trial mesh and its estimate, the polarization task: none of them changes the engine |
| `crates/nyx-engine/src/preview.rs` | 92 | feedback while editing: check an edit, size and cost from a trial mesh, the polarization as a function that owns its inputs |
| `crates/nyx-engine/tests/validated.rs` | 91, 101 | the three states of a result, a stale result never checked, the verdict kept by problem key, check not affordable from the memory estimate, the verdict in the run record, lossless tilted-crystal waveguide convergence with coupled-beta rounding marked as noise |
| `crates/nyx-engine/tests/reduce.rs` | 90 | the reductions against closed forms and against loops on the raw arrays, to the bit, for any thread count, on a worker's result and on real physics |
| `crates/nyx-engine/tests/stale.rs` | 84, 64 | staleness per physics, late jobs, the `stale` fields of the replies |
| `crates/nyx-engine/tests/format.rs` | 63 | the `format` key at project, feature list and problem level; old and foreign files refused |
| `crates/femcore/tests/expr_grammar.rs` | 65 | the two expression grammars (`geom::expr`, `femcore::expr`) accept and refuse the same corpus (R06) |
| `crates/nyx-engine/tests/golden/*` | 10, 42 | `solid_static_thermal_from_heat`, moved unchanged from femcore |
| `crates/nyx-engine/tests/swept_cell.rs` | 66 | a periodic swept cell with filled layers: closed form and refusal |
| `crates/nyx-engine/examples/*` | 41, 9 | benchmarks and logs that drive every physics |
| `crates/nyx-engine/src/validated.rs` | 91 | the validated state (exploratory, checked, check not affordable), the checks kept by problem key and the verdict written into the run record |
| `crates/nyx-engine/src/reduce.rs` | 90 | field handles and the reductions `integrate`, `extremum`, `probe`, `line_cut`, `overlap`, `confinement`: the quadrature rule, the order of the sums, the cell locator |
| `crates/nyx-engine/src/scriptout.rs` | 94 | export as script (`project_script`): the project as a command list replayed on a scratch engine to verify it, then rendered as Python; three forms (`commands`, `geometry_json`, `project_json`) |
| `crates/nyx-engine/src/shared_queue.rs` | 95 | the queue in the workspace folder: items and markers, the lease, the queue worker (`--worker`), `run_shared` for the job that waits for an item, the `nyx queue` commands |
| `crates/nyx-engine/tests/shared_queue.rs` | 95 | two front-ends on one queue, a submitter that exits, a solver killed, a queue worker killed, cancel and steering from another program |
| `crates/nyx-engine/tests/cqed_validated.rs` | 102, 103, 105, 91 | end-to-end superconducting CPW `Engine` workflow reaching `validated: checked` and running `capacitance_matrix` & `tls_participation` |
| `crates/nyx-engine/*` | 37, 38, 25, 26, 62, 99, 102, 103, 104, 105, 107 | `Project`, commands, jobs, store wiring, pure-Rust HDF5 results (`h5io.rs`, `xdmf.rs`, `results_file.rs`) and geometry descriptors; the registry from cargo features |
| `crates/nyx-script/src/{lib,host,convert,stubs}.rs` | 93, 94, 57 | scripts: the host (run and check a script through a worker process), JSON to and from script values, the `nyx` stub generated from the catalog |
| `crates/nyx-script/src/{trace.rs,prelude.nys}` | 90, 93 | tracing of script functions (N3): the Python prelude the worker is fed, the formula functions, the `nyx` wrapper that turns a function given as a formula into its text, and the static guard against a branch on a value of the point |
| `crates/nyx-script/src/kernel.rs` | 94, 93 | the kernel: one Monty session that outlives a script; commands that wait for a job served by polling so no call blocks for ever; `nyx.aio` futures; the refusal of the destructive commands; a dead worker replaced from a dump; `host::run` and `check` are one script on a kernel of their own |
| `crates/nyx-script/src/bin/nyx-script-worker.rs` | 93, 86 | the script worker program: Monty's worker mode, the only place the interpreter runs |
| `crates/nyx-script/{Cargo.toml,tests/script.rs,tests/trace.rs,tests/kernel.rs,tests/export.rs,tests/contract.rs}` | 94, 93, 90 | the nyx-script package, its tests (limits, the path guard, determinism, the static check, same numbers as plain commands) and the `script_api.v1` exemplar test |
| `crates/nyx-engine/src/path_guard.rs` | 93 | the path guard: file and folder arguments of commands kept inside one folder |
| `crates/nyx-engine/tests/path_guard.rs` | 93 | the path guard: every path-like argument of the catalog is guarded, no command writes or reads outside the root |
| `docs/scripting-samples/*` | 93 | a sample script (the strip sweep) |
| `docs/architecture/contracts/script_api.v1.json` | 93, 19 | what a script can rely on (the binding's rules; the commands are `engine_commands`) |
<!-- map:end -->

## Gaps this index shows

* Files with no ADR of their own: benchmarks and studies in `examples/`. They are measurements, not decisions.
* Modules whose ADR is retrospective and thin on rationale because the history gives only the subject line: `fem.rs`, `mesher3d.rs`, `refine3d.rs` (ADR 4, 7), `sms*.rs`, `ems*.rs`, `sbs*.rs` as physics (ADR 9, 12 are about structure and naming, not about the choice of formulation; `docs/HOW_IT_WORKS.md` is the explanation of the formulations).
* Decisions with no ADR because they are not recorded anywhere I could read: choice of `faer`, `spade`, `eframe`; choice of P1/P2 and Nedelec 1/2 as the element families (ADR 4 only states what exists); the Brillouin gain formulation. Ask paul if they should be written down.
