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
| `CONTRIBUTING.md` | 13, 19, 100 | contributor workflow, VVUQ checks, and inbound contribution licence |
| `SECURITY.md` | 93, 100, 131 | vulnerability reporting policy and security boundaries |
| `CITATION.cff` | 13, 100 | Citation File Format metadata for academic attribution |
| `.gitignore` | 13, 20 | `*.py` ignored on purpose |
| `docs/architecture/adr/*` | 19 | the ADRs themselves |
| `docs/architecture/contracts/*` | 19 | exemplars |
| `docs/architecture/*` | 19, 20 | model, rough spots, indexes |
| `docs/history/architecture/*` | 19, 20 | archived restructure plan, proposals and decision rounds |
| `docs/history/README.md` | 20 | what the archive is |
| `docs/README.md` | 20 | index of the docs folder |
| `docs/history/RESTRUCTURE_STATUS.md` | 20, 19 | handoff for the restructure plan |
| `docs/history/UX_STATUS.md` | 20, 88 | handoff for the UX review and plan |
| `docs/HOST_SETUP.md` | 97, 98, 94 | what costs time when an AI host drives nyx, and what to set in the host and do in the agent |
| `docs/history/LATENCY_STATUS.md` | 97, 20 | handoff for the latency work L1 to L3 |
| `docs/ux/*` | 88, 123 | colour tokens (`tokens.json`, compared with the app by a test) |
| `docs/history/ux/*` | 88, 123 | archived UX review and plans |
| `docs/vvuq/*` | 10, 19, 61, 91, 112, 120, 131 | ASME VVUQ & Rust formal verification plan and EGM-refined plan |
| `benchmarks/validation/*` | 133 | ASME V&V 20 versioned validation datasets (Park 2026, Palace/Q3D/SQuADDS, `scqubits`/TLS, photonics, MEMS/piezo, NAFEMS T1–T4) |
| `docs/history/FRAME_UNIFICATION_STATUS.md` | 1 | |
| `docs/history/KERNEL_PORT_STATUS.md` | 1, 9, 10, 12 | |
| `docs/history/GEOMETRY_ARCHITECTURE.md` | 2, 3, 18 | historical |
| `docs/history/GEOM_HANDOFF*.md` | 2, 3, 18 | historical |
| `docs/history/GEOMETRY_GUI_PLAN.md` | 2, 23 | |
| `docs/history/WORKFLOWS_DESIGN.md` | 14, 23 | |
| `docs/history/MCP_ASSISTANT_PLAN.md` | 15 | |
| `docs/history/SWEEP_CACHING_PROPOSAL.md` | 6 | proposal; historical |
| `docs/history/PML_FEATURE_REQUEST.md` | 128, 5, 7, 52, 69, 82, 101, 120 | COMSOL-inspired structured/curvilinear PML & open-boundary specification (ADR 128) |
| `docs/history/PML_IMPLEMENTATION_PLAN.md` | 128, 5, 7, 52, 69, 82, 101, 120 | implementation plan for native mixed-cell curvilinear PML/IE & physics scaling (ADR 128) |
| `docs/HOW_IT_WORKS.md` | 4, 9 | |
| `docs/history/PLOTTING_PLAN.md` | 99 | design note of the mode figures (the script, the results layout, the stages) |
| `docs/history/GDSFACTORY_INTEGRATION.md` | 57, 94 | design note: a gdsfactory plugin that drives the engine; not built |
| `docs/history/RESULTS_FORMAT_STUDY.md` | 62, 99, 107 | study of the results file layout and pure-Rust `hdf5-pure` decision |
| `docs/figures/example_modes_gallery.png` | 99 | example output of the `figure` command |
| `docs/history/chip_report.md` | 22 | |
| `docs/history/photonic_use_cases.md` | 21, 22 | |
| `docs/SQUADDS_BENCHMARK_REPORT.md` | 102, 103, 104, 105, 106, 114, 115, 116, 117, 118, 119, 120, 121 | SQuADDS cQED statistical parity, AF-1..AF-6 ablation, and cryogenic hardware validation report |
| `docs/history/SQUADDS_BUG_FEATURE_REQUESTS.md` | 11, 102, 103, 104, 105, 106, 114, 115, 116, 117, 118, 120, 121, 125, 126, 127 | Post-SQuADDS benchmark bug-fix and solver/meshing feature requests (FR-25..FR-33) covering Core Bugs 1-5 and Limitations L1-L4 |
| `docs/brand/*` | 23 | |
| `docs/ports/README.md` | 5, 7, 8 | |
| `docs/ports/NOTICES.md` | 7, 8 | licences |
| `docs/ports/meshers/*` | 5, 7, 16, 17, 18 | |
| `docs/history/ports/meshers/*` | 5, 7, 16, 17, 18 | archived staging notes of the femmesh port |
| `docs/ports/solvers/*` | 8, 11 | |
| `docs/history/ports/solvers/*` | 8, 11 | archived staging notes of the femsolve port |
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
| `crates/geom/src/topo3.rs` | 3, 132 | shell 2-manifold edge pairing & closed directed loop checks |
| `crates/geom/src/query.rs` | 3, 18 | exact queries, closest point |
| `crates/geom/src/curve.rs` | 3, 18 | |
| `crates/geom/src/tess.rs` | 3 | display only |
| `crates/geom/src/bvh.rs` | 116, 3, 132 | flat 2D AABB BVH (`AabbBvh2`, `ComplexBvh`) with depth > 64 spill stack accelerating arrangement, `Complex::locate`, and stack face lookup |
| `crates/geom/src/poly_batch.rs` | 102, 2, 116, 121 | batch polygon cleanup, spatial-hash shared-corner-pinned RDP simplification, and 1-DOF endpoint-anchored circular arc fitting |
| `crates/geom/src/gds.rs` | 114, 121, 3, 131 | pure-Rust binary GDSII stream parser, IBM float codec, hierarchy flattener with `GdsImportLimits`/ROI clipping, arc fitting, and unioned island/pin extractor |
| `crates/geom/src/*` | 3, 116, 121, 132 | arrangement (`Verts::cluster`, exact conic ray crossing), complex (`check3` quadric flux), planar, prim, scene, synth, math, lib |
| `crates/geom/tests/contracts.rs` | 19, 3 | |
| `crates/geom/tests/goldens.rs` | 10, 2 | |
| `crates/geom/tests/wafer.rs` | 1 | |
| `crates/geom/tests/sweep_axes.rs` | 1 | |
| `crates/geom/tests/selectors.rs` | 2 | |
| `crates/geom/tests/gds.rs` | 114, 3, 131 | binary GDSII parser, IBM hex float, hierarchy, `PATH`, keyhole cutline, `TEXT` pin, and `GdsImportLimits`/ROI DoS tests |
| `crates/geom/tests/arc_fitting.rs` | 121, 102, 114 | 1-DOF circular arc fitting (`Polygons::arc_tol` & `GdsImportOptions::arc_tol`) and `Curve::Circle` preservation tests |
| `crates/geom/tests/thick_sheet_stack.rs` | 121, 102 | conformal finite-thickness metal (`thickness`) and substrate trench (`etch_depth`) 3D stack tests |
| `crates/geom/tests/{fuzz_gds,bvh}.rs` | 131, 132, 114, 116, 3 | Phase 1 `bolero` GDSII fuzz target and `proptest`/deep-skewed BVH brute-force equivalence suite |
| `crates/geom/tests/*` | 3, 131, 132, 133 | geometry kernel unit/integration suites and Tier A `RAYON_NUM_THREADS` 0-ULP determinism tests |
| `crates/geom/golden/*` | 10, 2 | |
| `crates/geom/examples/*` | none | benchmarks and studies |
| `crates/geom/*.md` | 2, 3 | |
| `crates/geom/LICENSE` | 3 | |
| `crates/geom/rustfmt.toml` | 13 | only formatted crate (R23) |
| `crates/geom/Cargo.toml` | 3, 13 | serde only |
| `crates/femcore/src/linalg/equilibrate.rs` | 129, 106, 33 | 8-pass symmetric Ruiz equilibration (`RCsr`, `CCsr`, affine term terms) shared by `nyx-sms` and `nyx-ems` |
| `crates/femcore/src/linalg/prom.rs` | 129, 106, 33, 133 | multi-port Hermite + isometric-QR reduced-order model (`PROM`) sweep builder shared by `nyx-ems` and `nyx-sms` (Phase 3 QR/solve residual diagnostics) |
| `crates/femcore/src/linalg/*` | 33, 8, 11, 117, 125, 127, 129, 131, 132, 133 | ported solver library; `driver.rs` is the door the kernels use, `csr.rs` the one matrix type, `amg_elastic/` the elasticity multigrid, `arnoldi_real.rs` the real waveguide eigensolver, `direct.rs` pivot/inertia/backward-residual diagnostics, perturbed-shift unshift, & R31 recovery |
| `crates/femcore/src/mesh/mod.rs` | 4, 82, 131 | `SimplexMesh`, the solver contract (was `mesh.rs`; `SimplexMesh::check()`) |
| `crates/femcore/src/mesh/view.rs` | 4 | draw-only mesh |
| `crates/femcore/src/mesh/cells.rs` | 5, 82 | per-cell-type mesh (was `meshers/cellmesh.rs`) |
| `crates/femcore/src/mesh/exchange.rs` | 32, 41, 82, 128 | the mesh a physics is given (native cell mesh, 2D or 3D) and its numbering of unknowns |
| `crates/femcore/src/mesh/mixed.rs` | 128, 4, 5, 32 | native mixed-cell (`Tri`/`Quad`, `Tet`/`Prism`/`Hex`) solver discretisation, $D_3$/$D_4$ face canonicalization, and zero-pyramid conformity check |
| `crates/femcore/src/mesh/io/*` | 16 | `.fmsh` |
| `crates/femcore/src/mesh/{quality,fixtures}.rs` | 7 | quality statistics, test meshes |
| `crates/femcore/src/meshing/mod.rs` | 6, 7, 18, 50, 51, 52, 82 | `MeshPlan`, `mesh_model` (was `meshing.rs`) |
| `crates/femcore/src/meshing/geomdrive.rs` | 2, 6, 82, 131, 132 | |
| `crates/femcore/src/meshing/mesher.rs` | 6, 2, 66, 110, 126, 131, 132 | 2D triangulation core; `extrude_y` with `Layer::fill`, slanted periodic wall splitting (R29), `extrude_y_with_remap`, `sweep_model_3d`, and cross-gap dielectric edge splitting |
| `crates/femcore/src/meshing/mesher3d.rs` | 40, 55, 82, 83, 110, 131 | lattice (isosurface stuffing) tetrahedral mesher on `Domain3`, with parent->child periodic axis conformance |
| `crates/femcore/src/meshing/refine3d.rs` | 40, 110 | Rivara longest-edge refinement of a `CellMesh<3>`, with synchronous parent->child periodic twin bisection |
| `crates/femcore/src/meshing/curve3.rs` | 82, 110, 131, 132 | curved edges: mid-edge geometry nodes, 20-coeff Bernstein Jacobian certification (non-QP) vs interior `q8` (QP), periodic mid-edge sync |
| `crates/femcore/src/meshing/geomdrive3.rs` | 40, 2, 50, 52, 82, 83, 110, 132 | adapter: a 3D geom model as `Domain3` (classify, facet_project, curved crease projection R33, LabelMap, `PeriodicAxis3`), `mesh_geom_3d` |
| `crates/femcore/src/meshing/pml3.rs` | 5, 7 | |
| `crates/femcore/src/meshing/pml_geom3.rs` | 52, 5, 131 | structured absorbing layers for a meshed 3D geom model: the physical window, `build_pml`, geom labels, the tetrahedral cut |
| `crates/femcore/src/meshing/adapt3d.rs` | 51, 40 | adaptive meshing of a 3D geom model: the exact faces as the oracle of `adapt3`, interface facets, vertex sizes |
| `crates/femcore/src/meshing/request.rs` | 54, 82, 83, 114, 116, 120, 131 | the serialisable mesh request and report (`curved`, `gap_cells`, `z_coarsen`, `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, `max_steiner_multiplier`, `strict_mesh`, `singular_edges`, `quarter_point_xi`, `corner_grading`), their JSON Schemas |
| `crates/femcore/src/meshing/pml_geom.rs` | 5, 7, 52, 131 | structured PML for a meshed geom model, an adapter over `pml_box` (was `meshing_pml.rs`, then `meshing/pml.rs`) |
| `crates/femcore/src/meshing/pml_box.rs` | 7, 5, 128 | structured PML generator (was `meshers/adapt/pml.rs`) |
| `crates/femcore/src/meshing/{pml_cyl,pml_sphere}.rs` | 128, 52, 82 | structured 2D/3D cylindrical (`Quad`, barrel/endcap `Prism`, rim `Hex`) and spherical (`Quad` at $r=0$, concentric `Prism`) mesh sweepers with lockstep curved shells |
| `crates/femcore/src/meshing/physics_size.rs` | 17 | |
| `crates/femcore/src/meshing/sizing.rs` | 17, 25 | the neutral Sizing a physics reports to the mesher (re-exported by nyx-physics) |
| `crates/femcore/src/meshing/validate/*` | 7 | |
| `crates/femcore/src/meshing/sheet3d.rs` | 102, 114, 116, 120, 121, 126, 131, 132 | conformal 3D sheet-stack mesher (`Op::Sheet`, `<sheet>:edge`, PSLG guard, `gap_cells` medial-axis seeding, `z_coarsen` `max_cells` enforcement, thin slabs, curved/QP enrichment, periodic Dompierre sync) |
| `crates/femcore/src/meshing/sheet3d_amr.rs` | 120, 102, 114 | `OverlayOracle2` 2D-overlay and `full_3d` sheet-preserving adaptive mesh refinement helpers |
| `crates/femcore/src/meshing/sheet3d_curve.rs` | 121, 82, 102, 132 | exact 2D `Curve::Circle` / `Curve::Ellipse` mid-edge node projection and periodic twin sync on extruded `sheet3d` meshes |
| `crates/femcore/src/meshing/singularity.rs` | 120, 115, 82, 132 | Barsoum quarter-point P2 rim shifting (with curved-rim clamping, periodic twin sync) and Morrison–Lewis radical corner grading |
| `crates/femcore/src/par.rs` | 11, 133 | R03; unified fixed-chunk deterministic reductions and `chunk_map` folds across femcore and physics crates |
| `crates/femcore/src/adapt/electro.rs` | 120, 7, 53 | goal-oriented a posteriori electrostatic flux-jump estimator, 2D column projection, isotropic corner balls, and 3D metrics |
| `crates/femcore/src/adapt/estimate*.rs` | 7, 17, 53 | |
| `crates/femcore/src/adapt/work.rs` | 5, 4, 7 | the adaptation's private simplex working mesh and its conversion from and to the solver mesh (was `mesh/native.rs`, `mesh/interop.rs`) |
| `crates/femcore/src/adapt/maxwell.rs` | 53, 133 | analytic P1/P2/CurvedTet residual error estimator for edge elements and the equidistributing metric |
| `crates/femcore/src/adapt/{ops3,topo3,drive3}.rs` | 7, 51, 131 | 3D adaptation operators, working mesh, driver; overlap check, edge-removal and flip guards, compaction (ADR 51) |
| `crates/femcore/src/adapt/*` | 7 | metric adaptation, timer, mod |
| `crates/femcore/src/material/*` | 1, 29, 30, 34, 35, 101, 104, 108, 109, 133 | spec, units, elastic, tensor, optics, schema (`anyOf` / closed `temp_curve` & `model`), T-curves (`state_t.rs` `eval_checked` / `eval_strict`), `Source` uncertainty flags, Zimmermann BCS (`bcs.rs`), coupled tensors (`coupled.rs`), database importers (`import.rs`) |
| `crates/femcore/tests/material_import.rs` | 108, 109 | pure-Rust importer & catalog regression tests |
| `crates/femcore/tests/material_bcs.rs` | 104, 109 | Mattis–Bardeen BCS surface impedance and sub-gap mK loss floor tests |
| `materials/*` | 29, 30, 34, 35, 101, 104, 108, 109 | the shipped multi-source material library (TOML, `format = 2`) |
| `docs/MATERIALS.md` | 29, 30, 34, 35, 101, 108, 109 | material library guide, `format = 2` sections, variant catalog, and database import commands |
| `crates/nyx-engine/src/material_fetch.rs` | 108, 109, 131 | engine catalog, remote search, offline database conversion, and allowlisted/SSRF-hardened HTTPS material fetch |
| `crates/app/src/matcmd.rs` | 108, 109 | `nyx material` CLI subcommands (`list`, `show`, `eval`, `bcs`, `search`, `import`, `fetch`) |
| `crates/femcore/src/expr.rs` | 2 | duplicate of geom::expr (R06) |
| `crates/femcore/src/fieldexpr.rs` | 90 | the typed expression tree over complex field values and its evaluator (N1; N2 added comparisons, `where` and more functions) |
| `crates/femcore/src/formula.rs` | 90 | the text formula parser that produces that tree: grammar, vectors, names resolved by the caller (N2) |
| `crates/femcore/src/validate.rs` | 9 | |
| `crates/femcore/src/vvuq.rs` | 132, 133, 19 | shared VVUQ Layer 2 & Layer 3 zero-sized `SourceField`/`UniformSource` MMS trait, `assert_convergence_rate`, and `roache_gci` |
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
| `crates/femcore/src/farfield.rs` | 120, 4 | asymptotic 2D/3D Robin and 3D symmetric Bayliss–Turkel Laplace–Beltrami far-field boundary assemblers |
| `crates/femcore/src/rt.rs` | 120, 4 | $\mathrm{RT}_0/\mathrm{RT}_1$ $H(\mathrm{div})$ elements, sheet-split vertex-patch flux equilibration, and Prager–Synge Loewner bounds |
| `crates/femcore/src/fem.rs` | 4, 8, 82, 120 | element basics, original stack |
| `crates/femcore/src/{fem_mixed,edge_mixed}.rs` | 128, 4, 9 | native tensor-product Lagrange ($Q_1$, serendipity $Q_2'$, full $Q_2$, curved) and Nédélec ($N_1$, $N_2$) elements plus algebraic de Rham `AMS`/`AmsP2` operators |
| `crates/femcore/src/lib.rs` | 12 | module map |
| `crates/femcore/tests/linalg_contract.rs` | 33, 19 | driver contract and closed-form checks |
| `crates/femcore/tests/contracts.rs` | 19, 1, 4, 5, 9, 16, 51, 52, 53, 54 | |
| `crates/femcore/tests/material*.rs` | 29, 30, 34, 35, 19, 101, 104, 108, 109 | behaviour, BCS, and contract tests of the material format |
| `crates/femcore/tests/sheet_mesh3d.rs` | 102 | tests of conformal 3D sheet-stack meshing and `<sheet>:edge` tags |
| `crates/femcore/tests/sheet3d_qcap.rs` | 114, 102 | PSLG curved-hole regression, coplanar `gap_cells`, crack-free pure-tet `z_coarsen`, and thread determinism tests |
| `crates/femcore/tests/electro_amr.rs` | 120, 7 | electrostatic flux-jump estimator, corner-ball 2D metrics, `OverlayOracle2`, `full_3d` sheet AMR, and thread bit-determinism tests |
| `crates/femcore/tests/singularity_qp.rs` | 120, 115, 132 | Barsoum quarter-point shifting, curved-rim clamping, periodic twin sync, Morrison–Lewis corner grading, and `0/0`-free singular gradient tests |
| `crates/femcore/tests/rt_equilibration.rs` | 120, 4 | $\mathrm{RT}_0/\mathrm{RT}_1$ exact polynomial recovery, sheet half-star equilibration, and Prager–Synge Loewner bound tests |
| `crates/femcore/tests/sheet3d_curved_arcs.rs` | 121, 82, 102, 132 | curved `sheet3d` arc projection, 216-point collapsed quadrature, periodic twin sync, and concentric annular capacitance tests |
| `crates/femcore/tests/periodic_mesh3d.rs` | 110, 40 | 3D parent->child periodic vertex, edge, and `Ned3` triangle congruence across lattice & Rivara refinement |
| `crates/femcore/tests/golden_frames*` | 10, 1 | |
| `crates/femcore/tests/frame*.rs` | 1 | |
| `crates/femcore/tests/ring_frame.rs` | 1 | |
| `crates/femcore/tests/linalg/fr16.rs` | 117, 33 | supernodal `LBL^H`/`Sym LBL^T` direct solver benchmarks/tests and $P_2$ Nédélec `AmsP2` FGMRES/JD tests |
| `crates/femcore/tests/linalg/*` | 8, 11, 117, 132 | |
| `crates/femcore/tests/linalg/main.rs` | 8 | |
| `crates/femcore/tests/meshers.rs` | 7, 5 | |
| `crates/femcore/tests/adapt3d.rs` | 7 | |
| `crates/femcore/tests/mesh3d_features.rs` | 50, 51, 52, 55 | wavelength sizing, adaptive meshing, structured layers and thread-count independence for 3D geom models |
| `crates/nyx-ems/tests/maxwell_estimator.rs` | 53 | the Maxwell estimator and loop on a PEC cavity and on a sphere from geom |
| `crates/femcore/tests/meshing.rs` | 6, 66, 132 | |
| `crates/femcore/tests/curved.rs` | 82, 132 | curved tetrahedra: the map, the quadrature, storage, sphere volume, periodic twin creases, thread independence |
| `crates/femcore/tests/face_size.rs` | 83 | face sizes and `max_dev` in the lattice mesher |
| `crates/nyx-ems/tests/curved.rs` | 82, 132 | the cavity kernel on curved tetrahedra against closed forms, $O(h^4)$ PEC sphere cavity, and R34 curved PML post-solve integrals |
| `crates/nyx-engine/tests/curved.rs` | 82 | a periodic dielectric sphere with `curved` on and off |
| `crates/femcore/tests/geom_meshes.rs` | 2, 6 | |
| `crates/femcore/tests/pml_curvilinear.rs` | 128, 52, 82 | verification of all 7 canonical PML topologies, zero-pyramid conformity, lockstep curved shells, ray-invariance guards, and thread determinism |
| `crates/femcore/tests/{portable_kernels,vvuq_props}.rs` | 131, 11, 33 | `portable-kernels` `< 1e-13` differential parity and `proptest` mesh/Csr/formula invariants |
| `crates/femcore/tests/*` | 9, 12, 133 | physics tests: solid, solid_pml, anisotropic, common, `material_range_and_uncertainty` |
| `crates/femcore/examples/fs_*.rs` | 8 | ported CLI tools |
| `crates/femcore/examples/mesh_*.rs` | 7, 16 | mesher studies |
| `crates/femcore/examples/*` | none | benchmarks and studies; the numbers are in the history |
| `crates/femcore/Cargo.toml` | 13, 8, 11 | feature `portable-kernels` |
| `crates/nyx-storage/*` | 107, 62, 99, 131, 133 | pure-Rust HDF5 storage engine (`hdf5-pure` + `zlib-rs`), packed 5D Shuffle+LZF layout, streaming reader/writer with `max_solutions` & byte budgets, XDMF, sweep checkpoint (`SkippedSweepPoint`) |
| `crates/nyx-store/*` | 14, 39, 130, 131, 133 | JSON specs, settings (`scale_mode`, `zoom`), `.nys`/`.xdmf` archive rules, `store_run` `1.1` (`BuildProvenance`, `RuntimeProvenance`, `RunDiagnosticsSummary`, optional `Convergence` GCI, `material_hashes`), and the card function (39) |
| `crates/app/src/{tokens,theme,kit,menu,regions,activity,sweepui,cvd}.rs` | 88, 123, 124, 130, 133 | colour tokens as data, the looks built from them (`weak_text_color`), the widgets, the one table of actions (palette and native menu), the colour slots of regions (solid, hatched, grey), the busy chip, queue popover and stale band from the engine queue, the Sweep section of the Build page, the colour-vision simulation of the scripted shots |
| `crates/app/src/console.rs` | 123, 122, 131 | the console of the command-line modes on Windows: the program is of the GUI subsystem, the command-line modes join the console of the terminal that started it unless the caller redirected the handles (`// SAFETY:` audited) |
| `crates/app/tests/colour_literals.rs` | 88 | the ratchet: colour literals in the window code only go down |
| `crates/app/src/{feedback,spans}.rs` | 92, 130 | the problem and size lines of an edit (`Feedback::is_busy`), the polarization thread, and the timing spans of the frame log |
| `crates/app/src/{shell,wfview,sheets,dashboard}.rs` | 14, 23, 124, 130 | the workflow view: the one bar, the five tabs with the viewport always on, high-DPI scaling/persistence (`NYX_SCALE`, `scale_mode`), WSL workspace hint, responsive cards, non-overlapping toolbars |
| `crates/app/src/*` | 23, 130, 133 | window startup/WSLg backend selection, repaint bursts, opaque clear, pointer/zoom handling, mode plot integer ticks, Study validation/AMR controls, Results error-budget/diagnostic/tier UI, shell `--script` driver |
| `crates/app/*` | 23, 130 | |
| `crates/nyx-mcp/tests/{contracts,stdio}.rs` | 94, 97, 19, 15 | the `mcp_tools.v3` snapshot (four tools), the handshake, and the gate of S7: the P0 workflows in script mode against the per-command baseline |
| `crates/nyx-mcp/tests/{reset_session,multi_session}.rs` | 113, 94 | `reset_session` REPL global reset and isolated concurrent multi-session (`session`) MCP tests |
| `crates/nyx-mcp/*` | 94, 97, 15, 113, 131 | the MCP server: four script tools (`script_run`, `script_wait`, `script_check`, `api`) over isolated `McpSession` kernels with bounded framing/sessions |
| `crates/nyx-connect/*` | 15, 28, 19, 98 | one copy of the Connect logic (R07 resolved)|
| `.github/workflows/*` | 122, 131, 20, 133 | 3-tier VVUQ CI workflow (`vvuq-ci.yml`) and cross-OS Windows/Linux/macOS binary release packaging (`release.yml` via `nyx-pack`) |
| `.github/*` | 19, 100 | repository ownership (`CODEOWNERS`) and PR checklist (`pull_request_template.md`) |
| `tools/vvuq_ci.*` | 131, 20, 133 | local Linux/WSL and Windows PowerShell VVUQ 3-tier CI playbook runners |
| `tools/vvuq_telemetry.py` | 133, 20 | best-of-N quiet-runner wall-time, peak RSS (`ru_maxrss`/`K32GetProcessMemoryInfo`), and ASME V&V 20 telemetry ledger |
| `tools/*` | 20 | baseline and helper scripts |
| `.cargo/*` | 26, 27, 107 | shared target directory (no `HDF5_DIR` C dependency) |
| `crates/nyx-physics/src/decode.rs` | 41, 25, 109, 113, 133 | typed parameters, inline tensor decoding (`PropValue`), name resolution, and active-study `audit_active_materials` shared by every physics |
| `crates/nyx-physics/src/testkit.rs` | 41, 10, 110 | the golden-file harness |
| `crates/nyx-physics/src/frame.rs` | 56, 25 | the frame a result states for the analyses of its vector fields (axes, handedness, absorbing domains) |
| `crates/nyx-physics/src/cqed.rs` | 103, 104, 105, 106, 118, 120, 121 | typed cQED payloads (`CapacitanceMatrixData` bounds/loss/AMR history, `InductanceMatrixData`, `TlsMeta`, `CircuitRom`, `DrivenPortData`, `ElectroAmrData`) |
| `crates/nyx-physics/tests/inline_tensor.rs` | 113, 25, 129 | inline diagonal, 6-component, 3x3 symmetric, and `3x6` piezoelectric tensor (`PropValue` / `PropsTensor` / `Oriented`) tests |
| `crates/nyx-physics/*` | 25, 19, 113, 118, 120, 121, 129, 131, 133 | the Physics trait, envelope, descriptors, registry, `check_mesh`, `SolveResult` diagnostics/error budgets (`2.1`), and `check_result_finite` |
| `crates/nyx-dfs/tests/golden/*` | 10, 41, 132 | the dfs goldens, moved unchanged from femcore, plus Phase 2 `es3d_curved_coax` and `mag3d_mixed_loop` |
| `crates/nyx-dfs/src/ac_loss.rs` | 121, 103 | quasi-electrostatic AC dielectric conductance $G_{ij}(\omega)$, $Q_{ij}$ matrix, and positive complex Schur $\Delta C^{(2)} > 0$ |
| `crates/nyx-dfs/src/magneto.rs` | 103, 132 | gauge-stabilized London magnetostatics (`cross_section`, `planar_stream`, `mixed_3d`), `SourceField` MMS, and inductance matrix extraction |
| `crates/nyx-dfs/src/planar_bem.rs` | 115, 103, 121 | 2D layered-dielectric static & complex-loss boundary-element capacitance solver (`planar_bem`, image series + Wilton panels) |
| `crates/nyx-dfs/tests/ac_loss_electro.rs` | 121, 103 | $13\,\mathrm{GHz}$ silicon loss participation, DC $\sigma$ conductance, and strong-loss positive $\Delta C^{(2)} > 0$ tests |
| `crates/nyx-dfs/tests/farfield_electro.rs` | 120, 103 | isolated sphere ($4\pi\varepsilon_0 R$) and circular disk ($8\varepsilon_0 R$) single-terminal `far_field` capacitance tests |
| `crates/nyx-dfs/tests/capacitance_qubit.rs` | 103 | multi-terminal capacitance matrix & magnetostatic inductance benchmarks (Test 1) |
| `crates/nyx-dfs/tests/qcap_solvers.rs` | 115, 103, 132 | vertical-line AMG Block-CG bit-determinism, quarter-point P2 rim edges, and `planar_bem` SPD/Richardson verification |
| `crates/nyx-dfs/tests/open_boundaries_ie.rs` | 128, 120, 41 | 3-tier open-boundary verification (`heat_stationary` & `magneto` curl-Robin `far_field`, native `Quad`/`Prism`/`Hex` Rational & Exponential Infinite Elements, `heat_transient` guard) |
| `crates/nyx-dfs/tests/finite_guard.rs` | 131, 41 | non-finite parameter, material, mesh, and BC guard tests for `nyx-dfs` |
| `crates/nyx-dfs/*` | 41, 25, 12, 9, 103, 115, 118, 120, 121, 131, 132, 133 | diffusion: typed problem, multi-RHS electrostatics, Prager–Synge bounds, AC loss, planar BEM, magnetostatics, heat kernels, Layer 2 ASME 6-pillar & R32 verification, Phase 3 `diag.rs` & solver diagnostics/error budgets |
| `crates/nyx-sms/tests/golden/*` | 10, 42, 132 | the sms goldens (and `ring_model_axes`), moved unchanged from femcore, plus Phase 2 `solid_3d_lrb_viscoelastic` |
| `crates/nyx-sms/src/{piezo_mat,piezo_asm,piezo,piezo_driven}.rs` | 129, 42, 106, 128, 132, 133 | coupled piezoelectric `(u, phi)` constitutive condensation, simplex/mixed assembly, equilibrated direct / Schur Krylov-Schur / QEP / Schur-AMG Right-FGMRES static/eigen solvers, and driven multi-port `Y`/`Z`/`S` PROM sweep |
| `crates/nyx-sms/src/{axisym,cross_section,waveguide}.rs` | 43, 42, 68, 69, 80, 101, 129, 132, 133 | bodies of revolution on their half-plane and 2.5D waveguide cross-sections (`exp(-j q y)`): library entries, `axisymmetric_eigenfrequency` and `modes` studies, piezoelectric 4-block extensions, exact `m=±1` MPC (`R6a`), quadrature `s_theta` (`R6d`), cylindrical/planar PMLs, and mode scalars |
| `crates/nyx-sms/tests/axisym_mode_scalars.rs` | 80 | the numeric frequency parts and the energy shares against closed forms whose family is known |
| `crates/nyx-sms/tests/{pml_floquet,waveguide_modes}.rs` | 101 | Floquet-periodic absorbing layers (`P^H K P` congruence, inactive mass skip) and native 2.5D elastic waveguide modes (`study = "modes"`, closed-form SH cutoff, 3D swept Floquet agreement, PML complex frequencies and shares) |
| `crates/nyx-sms/tests/mixed_pml_and_ie.rs` | 128, 42, 69 | native mixed-cell elastodynamic PML, Christoffel backward-wave stabilization, signed-$Q$ amplifying-mode detection, outer clamping, and static Kelvin IE tests |
| `crates/nyx-sms/tests/{piezo,piezo_driven,piezo_waveguide,piezo_axisym}.rs` | 129 | coupled piezoelectric static/eigen/Schur-AMG (`piezo.rs`), multi-port driven admittance PROM (`piezo_driven.rs`), 2.5D waveguide SAW/Lamb (`piezo_waveguide.rs`), and axisymmetric `6mm` (`piezo_axisym.rs`) test suites |
| `crates/nyx-engine/src/convergence.rs` | 61, 80, 97, 101, 112, 119, 120, 133 | mesh convergence on a scratch engine: 3-grid bidirectional `factor`, `scale_sizes`, mode-signature `track`, `ModeScalarMeta`, Roache GCI & `r_eff < 1.3` fallback (`R11a`/`R11b`), sub-pF relative tolerance, and log-space SPD Richardson extrapolation |
| `crates/nyx-engine/src/solve.rs` | 61, 80, 133 | the solve chain in coupling order; puts the region names into the scalars a kernel names by domain number; auto-AMR dispatch for `ems::cavity` and `sms` |
| `crates/nyx-sms/*` | 42, 41, 25, 12, 9, 68, 101, 109, 129, 132, 133 | solid mechanics: typed problem, static, eigenfrequency, axisymmetric eigenfrequency, waveguide modes, driven multi-port piezo, absorbing-layer, viscoelastic damping (`r`-weighted `R6b`), Floquet kernels, Layer 2 6-pillar & `R32` verification, Phase 3 `diag.rs` & solver diagnostics/error budgets |
| `crates/nyx-ems/tests/golden/*` | 10, 45, 132 | the ems goldens, moved unchanged from femcore, plus Phase 2 R22 `.fmsh` solver-only oracles and 3D cavity/Bloch/PROM snapshots |
| `crates/nyx-ems/tests/axisym.rs` | 46, 45 | the axisymmetric kernel against closed forms (gate G1) |
| `crates/nyx-ems/tests/{emw2d,emw3d,pml,aniso_optics,bragg_stack,photonics,physics_meshing,meshing,moved}.rs` | 45, 1, 12, 67, 132 | the ems tests moved from femcore (plus un-ignored `soi_uniform_convergence_table`) |
| `crates/nyx-ems/src/bloch.rs` | 110, 45, 9 | 3-stage fixed-frequency 3D Floquet-Bloch waveguide and leaky grating solver (`Study::Bloch`) |
| `crates/nyx-ems/tests/{bloch3d,group_index_and_scalars}.rs` | 110, 45 | 3D Bloch modes, Brillouin-zone shifts, dispersive Poynting `v_g`/`n_g`, mode shares, and midplane symmetry tests |
| `crates/nyx-ems/src/{ibc,rom}.rs` | 104, 109, 125, 127, 132, 133 | superconducting thin-film IBC (`direct`/`perturbative`/`auto`, curved Jacobians), per-terminal `rank1` `lumped_jj`, physical $H$/$K$ normalization, passive `CircuitRom`, and `TlsMeta` |
| `crates/nyx-ems/tests/kicq_cpw_bcs.rs` | 104 | 3D CPW BCS kinetic inductance shift and lumped JJ equipartition tests (Test 3) |
| `crates/nyx-ems/src/{driven,ports}.rs` | 106, 125, 127, 128, 132, 133 | driven frequency-domain solver (`auto`/`direct`/`lblt`/`blr`/`ams_fgmres`), per-terminal `lumped_port`, curved 2D boundary-mode `waveport`, curvilinear PML power `P_pml`, and Hermite Tangential Galerkin PROM |
| `crates/nyx-ems/tests/driven_prom.rs` | 106, 127 | 3D driven PROM, waveport, BBQ admittance, and `FR-33` `blr`/`ams_fgmres` symbolic-reuse verification (Test 6) |
| `crates/nyx-ems/src/mixed_ops.rs` | 128, 45, 9, 132, 133 | native mixed-cell Maxwell assemblers, Gauss-point curvilinear PML pullbacks, mixed ports/IBC/ROM, and `Ams`/`AmsP2` curl-curl drivers |
| `crates/nyx-ems/tests/pml_mixed_curvilinear.rs` | 128, 45, 106 | native mixed-cell EM verification ($r=0$ spherical Bessel modes, $C_{4v}$ doublets, 2D Cartesian/Cylindrical reflection, `Ams`/`AmsP2`, driven `P_pml` PROM, Option B DOF & `LBL^H` speedup) |
| `crates/nyx-ems/tests/finite_guard.rs` | 131, 45 | non-finite parameter, material, mesh, and BC guard tests for `nyx-ems` |
| `crates/nyx-ems/*` | 45, 25, 12, 9, 46, 67, 104, 106, 109, 110, 125, 127, 128, 131, 132, 133 | electromagnetic waves: typed problem, mode, complex-mode (`emw_pml_x` `+4` Krylov/`layer_fraction`), axisymmetric (`R6d`), cavity (`SourceField` MMS, curved R34 integrals), Bloch, driven, Layer 2 6-pillar & `R32` verification, Phase 3 `diag.rs` & solver diagnostics/error budgets |
| `crates/nyx-analysis/src/capacitance.rs` | 103, 115, 118, 120, 121, 125, 126 | `capacitance_matrix` (Schur elimination, `Xmon`/`TransmonPocket` reduction, 3D `metal_thickness_m` outer-die guard, Prager–Synge Loewner brackets, AC conductance/$Q$, 3-level Stieltjes–Padé AMR extrapolation) and `inductance_matrix` |
| `crates/nyx-analysis/src/view3d.rs` | 111, 99 | un-averaged 8-subtet P2 3D slices, cutaways, marching-tets isosurfaces, mesh views, and Bloch phase tiling |
| `crates/nyx-analysis/src/network.rs` | 118, 105, 106, 133 | unified `NetworkAnalysis` (`id: "network"`): companion QEP, ROM `Z/Y/S`, Schur BBQ, Taubin $S_{21}$ circle fit, Touchstone v1.1, §5.5 convergence/residual diagnostics |
| `crates/nyx-analysis/src/{epr,hamiltonian,tls}.rs` | 105, 118, 121, 125, 127, 133 | `epr` KICQ analysis, `hamiltonian` (`normal_fock` / `craig_bampton`, `f_cav_Hz`/`g_matrix_Hz`/`kappa_cav_Hz`/`t1_purcell_s`, `sweep_port`, DC-SQUID), and BVH/TOML/hybrid quasi-electrostatic `tls_participation` with §5.5 diagnostics |
| `crates/nyx-analysis/tests/wave1_cqed_analyses.rs` | 118, 105 | verification suite for `NetworkAnalysis`, `hamiltonian` SQUID/Purcell/CB, 3D `I_edge` thickness capacitance, and BVH/TOML TLS |
| `crates/nyx-analysis/tests/capacitance_bounds.rs` | 120, 103, 115 | Prager–Synge Loewner bound propagation (`[E_C_low_Hz, E_C_up_Hz]`) and 3D transmon midpoint accuracy gain tests |
| `crates/nyx-analysis/tests/tls_thick_sheet.rs` | 121, 105 | conformal thick-metal and substrate-trench 3D TLS sidewall classification without knife-edge double-counting |
| `crates/nyx-analysis/tests/cqed_capacitance_reduction.rs` | 115, 103 | analytical verification of Schur elimination, Eqs. 17-18 transmon reduction, and resonator `g_Hz` / `chi_Hz` |
| `crates/nyx-analysis/tests/figure3d.rs` | 111, 99 | 3D figure kinds, dielectric jump preservation, Bloch tiling, latency, and `NYX_THREADS` bit-exact PNG tests |
| `crates/nyx-analysis/tests/{tls_benchmark,park2026_transmon_benchmark,fluxonium_benchmark}.rs` | 105 | quantitative cQED acceptance benchmarks (Tests 2, 4, 5) |
| `crates/nyx-analysis/tests/squadds_smoke.rs` | 103, 105, 118, 120, 125 | PR CI smoke suite for TransmonCross Mode G2 + `far_field` + `I_edge` (`bug1`), DC-SQUID `Dissipator` sweep, and `MUNINN` 3D TLS |
| `crates/nyx-analysis/src/piezo_resonator.rs` | 129, 118 | piezoelectric resonator mBVD 6-param fit, Schur-reduced multi-port admittance, phase-slope `Q_s`/`Q_p`, and IEEE `k_t^2`/`k_eff^2` |
| `crates/nyx-analysis/*` | 56, 25, 21, 99, 103, 105, 111, 115, 127, 129, 131, 133 | analyses of results (polarization, capacitance/inductance, EPR, circuit Hamiltonian, TLS, piezo_resonator, zero-solve `uq_cqed` with ASME V&V 20 interval separation, 2D/3D mode figures, §5.5 diagnostics) |
| `crates/nyx-engine/tests/layering.rs` | 26, 107 | the layering table |
| `crates/nyx-engine/tests/jobs.rs` | 85 | every long operation is a job: cancel tests per kind, the dropped job of a replaced project, the window never blocks |
| `crates/nyx-engine/tests/queue.rs` | 87 | the queue: lanes, priority, round robin, ageing, supersede, hold, shared stages, order never changes a number |
| `crates/nyx-engine/tests/commandlist.rs` | 89 | command lists: static check against the catalog schemas, the dry-run plan, a bad line named and nothing run, the P0 workflows from their files |
| `crates/nyx-engine/tests/api_docs.rs` | 89, 57 | the API reference: `docs/api/` is what the catalog generates, every command has a page and an index row, every example passes `nyx check`, every argument is described |
| `docs/api/*` | 89, 57 | the API reference, one page per engine command, generated from `Engine::catalog()` (never edited by hand) |
| `crates/app/tests/cli.rs` | 89 | `nyx check` and `nyx run`: the plan, exit codes, a bad list runs nothing |
| `crates/nyx-engine/tests/sweep.rs` | 87, 133 | sweeps: a point is an edit list, one recorded run per point, the numbers of a hand-made solve, the engine's own project and results untouched, a bad point stops the request |
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
| `crates/nyx-engine/examples/*` | 41, 9, 133 | quantitative benchmarks and logs that drive every physics (`examples.rs` checks closed forms, passivity, convergence) |
| `crates/nyx-engine/src/validated.rs` | 91, 112, 133 | the validated state (`exploratory`, `checked`, `check_not_affordable`), `state` + `tier` (`MeshConverged` / `BoundCertified`), the checks kept by problem key and the verdict written into `/provenance/validated_json` and the run record |
| `crates/nyx-engine/src/reduce.rs` | 90, 112, 133 | field handles and the reductions `integrate`, `extremum` (with P2 edge-midpoint bound), `probe`, `line_cut`, `overlap`, `confinement` (with quadrature warnings), including single-pass bit-exact batch evaluation |
| `crates/nyx-engine/src/amr.rs` | 120, 61, 85, 133 | multi-pass goal-oriented electrostatic, `ems::cavity`, and `sms` AMR loops (`sheet_2d_z`, `full_3d`, and 2D), `retain_amr_data`, and `AdaptPassRecord` tracking |
| `crates/nyx-engine/src/gds.rs` | 114, 120, 121, 93 | `Engine::import_gds` (`arc_tol`, `lateral_pad`) and `Engine::add_chip_stack` 3D stack builder behind `path_guard` |
| `crates/nyx-engine/cards/cqed.md` | 113, 94 | embedded `cqed` guide card for `get_guide` and `nyx-mcp`'s `api` tool |
| `crates/nyx-engine/src/scriptout.rs` | 94, 114, 116, 120 | export as script (`project_script`): the project as a command list replayed on a scratch engine to verify it, then rendered as Python; three forms (`commands`, `geometry_json`, `project_json`) |
| `crates/nyx-engine/src/shared_queue.rs` | 95, 133 | the queue in the workspace folder: items and markers, configurable `LEASE_DEAD` (`NYX_QUEUE_LEASE_DEAD_SECS`), the queue worker (`--worker`), `run_shared`, and the `nyx queue` commands |
| `crates/nyx-engine/tests/shared_queue.rs` | 95, 133 | two front-ends on one queue, configured lease timeout recovery, submitter/solver/worker kill, cancel and steering tests |
| `crates/nyx-engine/tests/cqed_validated.rs` | 102, 103, 105, 91 | end-to-end superconducting CPW `Engine` workflow reaching `validated: checked` and running `capacitance_matrix` & `tls_participation` |
| `crates/nyx-engine/tests/{convergence_sizes,convergence_extrapolate}.rs` | 112, 61, 80, 91 | `scale_sizes`, mode-signature `track`, bidirectional `factor`, and log-space SPD Richardson extrapolation tests |
| `crates/nyx-engine/tests/figure3d.rs` | 111, 99 | `Engine::figure` 3D kinds, auto-resolved `periodic_shift` / `periodic_phase`, and latency tests |
| `crates/nyx-engine/tests/import_gds.rs` | 114, 102, 120 | synthetic GDSII transmon `Engine::import_gds`, `add_chip_stack`, `sheet3d` area conservation, `<sheet>:edge`, and `path_guard` tests |
| `crates/nyx-engine/tests/periodic_photonics_e2e.rs` | 110, 111, 112, 113 | end-to-end periodic photonics (`PFR-1..PFR-7`) integration tests through `Engine::call` |
| `crates/nyx-engine/tests/cqed_qcap_e2e.rs` | 112, 113, 114, 115 | end-to-end superconducting qubit capacitance (`QFR-1..QFR-7`) integration tests through `Engine::call` |
| `crates/nyx-engine/tests/squadds_fr_e2e.rs` | 116, 117, 118, 119 | end-to-end SQuADDS (`FR-01..FR-24`) integration tests through `Engine::call` |
| `crates/nyx-engine/tests/accuracy_features_e2e.rs` | 120, 121 | end-to-end Phase A sub-$0.1\%$ electrostatic accuracy features (`AF-1..AF-6`) and Cross-Feature Compatibility Rules 1–5 tests |
| `crates/nyx-engine/tests/squadds_engine_smoke.rs` | 102, 104, 105, 106, 116, 118, 125, 126, 11 | PR CI `Engine::call` smoke suite for `CLT` QEP vs Taubin driven fit, `QubitCavity` M1/M2/M3 equivalence, `Hole`-post airbridges, `bug2`–`bug5`, and thread determinism |
| `crates/nyx-engine/tests/fr25_33_e2e.rs` | 125, 126, 127 | end-to-end post-SQuADDS (`FR-25..FR-33`) integration tests through `Engine::call` (`ibc_mode: "perturbative"`, hybrid quasi-electrostatic TLS, and `ems::driven` `blr`/`ams_fgmres` + memory calibration) |
| `crates/nyx-engine/tests/piezo_driven_network.rs` | 129, 106, 118 | Mason 1D BAW driven piezo FEM -> `NetworkAnalysis` and `piezo_resonator` mBVD round-trip tests |
| `crates/nyx-engine/*` | 37, 38, 25, 26, 62, 99, 102, 103, 104, 105, 107, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120, 121, 127, 129, 130, 131, 132, 133 | `Project`, `build.rs` provenance metadata, commands, jobs, `material_warnings.rs`, `qmc.rs`, `provenance.rs` (material library BLAKE3 drift check), `vvuq_report.rs` (`vvuq_export`), store wiring, pure-Rust HDF5 results (`h5io.rs`, `xdmf.rs`, `results_file.rs` v4), Phase 3–5 tests (`vvuq_phase3.rs`, `vvuq_phase4_*.rs`, `vvuq_phase5_*.rs`, `validation_*.rs`, `common/mod.rs`) and geometry descriptors |
| `crates/nyx-script/src/{lib,host,convert,stubs}.rs` | 93, 94, 57, 113 | scripts: the host (run and check a script through a worker process, BOM stripping), JSON to and from script values, the `nyx` stub generated from the catalog |
| `crates/nyx-script/src/{trace.rs,prelude.nys}` | 90, 93, 113 | tracing of script functions (N3, including dict recursion): the Python prelude, formula functions, `nyx` wrapper, and static guard |
| `crates/nyx-script/src/kernel.rs` | 94, 93, 113, 131 | the kernel: Monty session with per-turn suspension budget, proactive dump recycle, `reset`, `nyx.aio` futures, timeout/output clamping, network refusal, and worker recovery |
| `crates/nyx-script/src/bin/nyx-script-worker.rs` | 93, 86 | the script worker program: Monty's worker mode, the only place the interpreter runs |
| `crates/nyx-script/tests/resilience.rs` | 113, 93, 94 | per-turn suspension budget rollover, fatal turn recovery, `reset`, UTF-8 BOM, and dict-traced lambda tests |
| `crates/nyx-script/tests/accuracy_features_script.rs` | 120, 121, 93 | Monty `.nys` script static check and execution of `AF-1..AF-6` electrostatic accuracy features |
| `crates/nyx-script/{Cargo.toml,tests/script.rs,tests/trace.rs,tests/kernel.rs,tests/export.rs,tests/contract.rs,tests/fuzz_script_trace.rs}` | 94, 93, 90, 131 | the nyx-script package, its tests (limits, the path guard, determinism, the static check, same numbers as plain commands), `bolero` trace/formula/timeout fuzz target, and the `script_api.v1` exemplar test |
| `crates/nyx-engine/src/path_guard.rs` | 93 | the path guard: file and folder arguments of commands kept inside one folder |
| `crates/nyx-engine/tests/path_guard.rs` | 93 | the path guard: every path-like argument of the catalog is guarded, no command writes or reads outside the root |
| `docs/scripting-samples/*` | 93 | a sample script (the strip sweep) |
| `docs/architecture/contracts/script_api.v1.json` | 93, 19, 131 | what a script can rely on (the binding's rules; the commands are `engine_commands`) |
| `docs/architecture/contracts/install_manifest.v1.*` | 122, 19 | versioned installer receipt (`format = 1`) and JSON Schema |
| `docs/architecture/contracts/install_report.v1.*` | 122, 19 | machine-readable `nyx-setup --json` report and JSON Schema |
| `crates/nyx-install/*` | 122, 28, 88, 95, 100, 107, 131 | cross-platform pure-Rust installer (`nyx-setup`, `nyx-uninstall`), `nyx-pack` packager, OS desktop/PATH/file-assoc/ARP/Nyx.app drivers, 7-phase atomic journal, optional `glow` GUI wizard, and contract/hermetic/E2E tests |
| `website/*` | 134, 13, 60, 88, 89, 108, 133 | static GitHub Pages project hub (Zola `0.23.6`), standalone `site-prep` extractor & confidentiality gate, generated tokens/data, and Swiss Lab Monograph templates |
| `tools/gen_site_figures.*` | 134, 99, 111 | public website mode/3D cutaway plate figure generator |
<!-- map:end -->

## Gaps this index shows

* Files with no ADR of their own: benchmarks and studies in `examples/`. They are measurements, not decisions.
* Modules whose ADR is retrospective and thin on rationale because the history gives only the subject line: `fem.rs`, `mesher3d.rs`, `refine3d.rs` (ADR 4, 7), `sms*.rs`, `ems*.rs`, `sbs*.rs` as physics (ADR 9, 12 are about structure and naming, not about the choice of formulation; `docs/HOW_IT_WORKS.md` is the explanation of the formulations).
* Decisions with no ADR because they are not recorded anywhere I could read: choice of `faer`, `spade`, `eframe`; choice of P1/P2 and Nedelec 1/2 as the element families (ADR 4 only states what exists); the Brillouin gain formulation. Ask paul if they should be written down.
