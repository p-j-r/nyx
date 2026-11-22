# Architecture decision records

Nygard format (`Date`, `Status`, `Context`, `Decision`, `Consequences`), plus `Code` (files, commits) and `Contract` (what pins it). Written retrospectively on 2026-11-03 from `git log`, the status and handoff documents in `docs/`, `docs/ports/` and the code; the date of each is the date of the decision where the history gives one. An ADR is never edited to say something different: a later ADR supersedes it.

| # | Decision | Module |
|---|---|---|
| 1 | [One public wafer frame, propagation along +y](adr/0001-one-public-wafer-frame.md) | `section`, kernels, geom |
| 2 | [Geometry is a feature list addressed by names and selectors](adr/0002-geometry-as-named-feature-list.md) | `geom`, `geomdrive` |
| 3 | [geom is a lean, separate, analytic kernel](adr/0003-geom-is-a-lean-separate-kernel.md) | `geom` |
| 4 | [SimplexMesh is the one mesh the kernels read](adr/0004-simplexmesh-is-the-solver-contract.md) | `mesh`, `solve` |
| 5 | [Meshes stay native: one cell type per cell](adr/0005-native-mixed-cell-mesh.md) | `cellmesh`, `meshers`, `pml3` |
| 6 | [One meshing interface (MeshPlan)](adr/0006-one-meshing-interface.md) | `meshing` |
| 7 | [Metric adaptation and structured PML](adr/0007-metric-adaptation-and-structured-pml.md) | `meshers`, `estimate` |
| 8 | [Ported solver library merged, not wired](adr/0008-ported-solver-library-staged-not-wired.md) | `linalg` |
| 9 | [Every kernel is three stages](adr/0009-kernels-are-three-stages.md) | kernels, `solve` |
| 10 | [Golden files are the numerical oracle](adr/0010-golden-files-are-the-numerical-oracle.md) | tests |
| 11 | [Determinism by fixed chunks](adr/0011-determinism-by-fixed-chunks.md) | `par`, `linalg::par` |
| 12 | [Three-letter physics codes](adr/0012-three-letter-physics-codes.md) | `dfs sms ems sbs` |
| 13 | [Pure Rust, no Python, HDF5 optional](adr/0013-pure-rust-no-python.md) | build |
| 14 | [Workflow store: JSON manifests, `.nyx` archives](adr/0014-workflow-store-json-manifests-and-nyx-archives.md) | `nyx-store` |
| 15 | [AI interface: separate workspace, stdio MCP](adr/0015-ai-interface-separate-workspace-and-mcp-server.md) | `ai/` |
| 16 | [`.fmsh` lossless mesh format](adr/0016-fmsh-lossless-mesh-and-result-format.md) | `meshers::io` |
| 17 | [Physics-driven mesh sizing](adr/0017-physics-driven-mesh-sizing.md) | `physics_size` |
| 18 | [Curved elements deferred](adr/0018-curved-elements-deferred.md) | mesher, geom |
| 19 | [Versioned interface contracts](adr/0019-versioned-interface-contracts.md) | contracts |
| 20 | [Version control: jj, branch per task](adr/0020-version-control-jj-branch-per-task-status-files.md) | process |
| 21 | [Polarization from complex fields](adr/0021-polarization-analysis-from-complex-fields.md) | `analysis` |
| 22 | [Use cases and chips live in femcore](adr/0022-use-cases-and-chips-live-in-femcore.md) | `usecases`, `chips` |
| 23 | [The app is a shell over one Project](adr/0023-app-is-a-thin-shell-over-one-project.md) | `app` |
| 24 | [Brillouin and the ring chips removed](adr/0024-brillouin-and-chips-removed.md) | `sbs`, `chips` (gone) |
| 25 | [A physics is a crate behind the `Physics` trait](adr/0025-physics-trait-registry-and-problem-envelope.md) | `nyx-physics`, `nyx-dfs/sms/ems`, `nyx-engine` |
| 26 | [One workspace and a layering test](adr/0026-one-workspace-and-a-layering-test.md) | `Cargo.toml`, layering |
| 27 | [`nyx2`: second repository, a jj workspace per session](adr/0027-nyx2-repository-and-jj-workspace-per-working-session.md) | process |
| 28 | [`nyx-connect` is one crate with versioned reports](adr/0028-nyx-connect-is-one-crate-with-a-versioned-report.md) | `nyx-connect`, app Connect menu |
| 29 | [Materials are TOML files with explicit units and a checked stiffness](adr/0029-material-files-are-toml-with-explicit-units.md) | `material`, `materials/` |
| 30 | [Optical properties are dispersion models evaluated at the study wavelength](adr/0030-optical-dispersion-models-evaluated-at-the-study-wavelength.md) | `material::optics` |
| 31 | [femcore holds meshes, meshing and adaptation in three module trees](adr/0031-mesh-meshing-adapt-module-layout.md) | `femcore::{mesh, meshing, adapt}` |
| 32 | [One native mesh; `SimplexMesh` as the solver view](adr/0032-one-native-mesh-simplexmesh-as-a-view.md) | `femcore::mesh` (proposed) |
| 33 | [One linear algebra: the ported library behind one driver](adr/0033-one-linear-algebra-driver-facade.md) | `linalg::driver`, phase 3 |
| 34 | [A material may extend another and replaces whole tables](adr/0034-a-material-may-extend-another-and-replaces-whole-tables.md) | `material::spec` (`extends`) |
| 35 | [Provenance, oscillator models and the checked library](adr/0035-material-provenance-oscillator-models-and-the-checked-library.md) | `material::{source,optics}`, `materials/` |
| 37 | [The engine core is built before the kernels, with two seams](adr/0037-engine-core-before-the-kernels.md) | `nyx-engine` |
| 38 | [Geometry operations are described once](adr/0038-geometry-operations-are-described-once.md) | `nyx-engine::geomdesc` |
| 39 | [The store takes spec text and a card function](adr/0039-the-store-takes-spec-text-and-a-card-function.md) | `nyx-store`, `nyx-engine::store_link` |
| 40 | [3D geom models meshed by the lattice mesher (`Domain3`)](adr/0040-3d-geom-models-meshed-by-the-lattice-mesher.md) | `femcore::meshing::{mesher3d, refine3d, geomdrive3}` |
| 41 | [A physics crate owns its typed problem; the window reaches it through a legacy bridge](adr/0041-physics-crate-pattern-and-the-legacy-bridge.md) | `nyx-dfs`, `nyx-physics::decode`, `nyx-engine::legacy` |
| 42 | [Solid mechanics moves into `nyx-sms`](adr/0042-solid-mechanics-moves-into-nyx-sms.md) | `nyx-sms`, `nyx-engine::legacy::sms` |
| 43 | [Bodies of revolution stay a library entry of `nyx-sms`](adr/0043-axisymmetric-solid-kernel-stays-a-library-entry.md) | `nyx-sms::{axisym, cross_section}` |
| 44 | [What the sms move leaves in femcore](adr/0044-what-the-sms-move-leaves-in-femcore.md) | `femcore::{periodic, ems_pml::cell_stretches, progress}` |
| 45 | [The EM kernels move into `nyx-ems`; what they share stays in `femcore`](adr/0045-ems-moves-into-nyx-ems.md) | `nyx-ems`, `femcore::{edge, pml, optical}` |
| 46 | [The axisymmetric EM kernel is checked against closed forms and gets a study](adr/0046-axisymmetric-ems-checked-and-given-a-study.md) | `nyx-ems::axisym` |
| 47 | [`ems` in the legacy bridge; the use cases take their solver as an argument](adr/0047-ems-in-the-legacy-bridge-and-the-use-cases.md) | `nyx-engine::legacy::ems`, `femcore::usecases` |
| 50 | [Wavelength sizing in 3D: one element size per region](adr/0050-wavelength-sizing-in-3d.md) | `femcore::meshing::mesh_for_sizing` |
| 51 | [Adaptive meshing of 3D geom models: `adapt3` with the exact faces as the oracle](adr/0051-adaptive-meshing-of-3d-geom-models.md) | `femcore::meshing::adapt3d`, `adapt::{ops3,topo3,drive3}` |
| 52 | [Structured absorbing layers in 3D for geom models](adr/0052-structured-absorbing-layers-in-3d.md) | `femcore::meshing::pml_geom3` |
| 53 | [A Maxwell residual estimator, and an adaptation loop that takes any estimator](adr/0053-maxwell-residual-estimator-and-the-generic-adaptation-loop.md) | `femcore::adapt::{maxwell,estimate3}` |
| 54 | [A mesh is asked for and reported in data](adr/0054-mesh-request-and-report-as-data.md) | `femcore::meshing::request` |
| 55 | [The lattice mesher welds slivers and does not depend on the thread count](adr/0055-lattice-mesher-welds-slivers-and-is-independent-of-the-thread-count.md) | `femcore::meshing::mesher3d` |
| 56 | [The polarization analysis is a function of a result](adr/0056-polarization-analysis-is-a-function-of-a-result.md) | `nyx-analysis`, `nyx_physics::frame` |
| 57 | [The window and the MCP server are thin front-ends: every tool is an engine command](adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md) | `nyx-mcp`, `nyx-engine`, `app` |
| 58 | [The engine resolves named materials from the library and the project](adr/0058-the-engine-resolves-materials-from-the-library-and-the-project.md) | `nyx-physics::decode`, `nyx-engine::materials` |
| 59 | [The store knows nothing of a project](adr/0059-the-store-knows-nothing-of-a-project.md) | `nyx-store` |
| 60 | [The examples are project files; the use cases are gone](adr/0060-examples-are-project-files-and-use-cases-are-gone.md) | `examples/`, `app/examples.rs` |
| 61 | [Mesh convergence is an engine command](adr/0061-convergence-is-an-engine-command.md) | `nyx-engine::convergence` |
| 62 | [Results and meshes are stored in HDF5; JSON is for small things](adr/0062-results-are-stored-in-hdf5.md) | `nyx-engine::h5io` |
| 63 | [One `format` key, JSON for small files, HDF5 for bulk data](adr/0063-the-format-key-and-no-json-for-bulk-data.md) | `nyx-engine::project`, `geom::model`, `nyx-physics::spec`, `nyx-store` |
| 64 | [Identity: one persisted hash, the rest are cache keys](adr/0064-identity-hashes-one-persisted-the-rest-are-caches.md) | `nyx-store::util`, `nyx-engine::engine`, `geom::math` |
| 65 | [Closing the restructure: what the final structure replaces](adr/0065-closing-the-restructure-what-the-structure-replaced.md) | `docs/architecture`, `Cargo.toml` |
| 66 | [A swept section can be a periodic cell: sorted vertices and `fill` layers](adr/0066-sweep-layers-can-be-filled-and-periodic-cells-are-sorted.md) | `femcore::meshing::mesher`, `nyx-engine::meshing`, `nyx-engine::project` |
| 67 | [`pmc`: a magnetic wall for `ems` cavities and Floquet cells](adr/0067-a-magnetic-wall-boundary-for-ems-cavities.md) | `nyx-ems::problem`, `nyx-ems::cavity` |
| 68 | [The axisymmetric solid gets a study](adr/0068-the-axisymmetric-solid-gets-a-study.md) | `nyx-sms::{problem,lib,cross_section}` |
| 69 | [Axis conditions and the cylindrical stretch of the axisymmetric solid kernel](adr/0069-axis-conditions-and-cylindrical-stretch-of-the-axisymmetric-solid.md) | `nyx-sms::axisym`, `femcore::pml` |
| 80 | [Modes carry what identifies them, and a convergence check follows a mode](adr/0080-modes-carry-what-identifies-them-and-a-convergence-check-follows-a-mode.md) | `nyx-sms::{cross_section,axisym}`, `nyx-engine::{convergence,solve}` |
| 82 | [Curved (isoparametric, quadratic) tetrahedra for the `ems` cavity study; amends ADR 18](adr/0082-curved-order-2-geometry.md) | `femcore::fem`, `femcore::meshing::curve3`, `nyx-ems::cavity` |
| 83 | [Size rules on faces and `max_dev` act in the 3D lattice mesher](adr/0083-face-size-rules-and-max-dev-in-the-lattice-mesher.md) | `femcore::meshing::geomdrive3`, `mesher3d` |
| 84 | [Results are stale per physics: the key is the problem key](adr/0084-results-are-stale-per-physics-the-key-is-the-problem-key.md) | `nyx-engine::engine`, `app::app` |
| 85 | [Every long operation is a job that can be watched and cancelled](adr/0085-every-long-operation-is-a-job.md) | `nyx-engine::{engine,jobs,convergence,results_file}`, `app::app` |
| 86 | [Solves run in a worker process; their memory is estimated from what workers measured](adr/0086-solves-run-in-a-worker-process.md) | `nyx-engine::{worker,estimate,engine,jobs,process}`, `app::main`, `nyx-mcp::main` |
| 87 | [Work is queued in lanes; stages are shared by content](adr/0087-work-is-queued-in-lanes-and-stages-are-shared-by-content.md) | `nyx-engine::{jobs,queue,sweep,engine,solve,convergence}` |
| 88 | [Colour has roles, and the window follows the system appearance (accepted)](adr/0088-colour-has-roles-and-the-window-follows-the-system.md) | `app::{tokens,theme,kit,menu,colormap,shell}` |
| 89 | [Command lists are a front-end: checked, planned and run](adr/0089-command-lists-are-a-front-end-checked-planned-and-run.md) | `nyx-engine::commandlist`, `app::runfile` |
| 90 | [Whole-array operations on results are engine commands over a typed expression tree evaluated in Rust](adr/0090-whole-array-operations-are-engine-commands-over-a-typed-tree.md) | `nyx-engine::{reduce,engine}`, `femcore::fieldexpr`, `app::tabs` |
| 91 | [Validated is a state: a result is checked only by a convergence check on its own problem key](adr/0091-validated-is-a-state-checked-against-the-problem-key.md) | `nyx-engine::{validated,convergence,engine}`, `nyx-store::runs`, `app::{activity,kit}` |
| 92 | [Feedback while editing is a value, never a job](adr/0092-feedback-while-editing-is-a-value-never-a-job.md) | `nyx-engine::preview`, `app::{feedback,spans}` |
| 93 | [Scripts are Python-syntax text run by Monty in a worker process, behind an engine path guard](adr/0093-scripts-are-python-text-run-by-monty-in-a-worker-process.md) | `nyx-script`, `nyx-engine::path_guard`, `app::{runfile,apicmd}` |
| 94 | [The MCP server is four script tools over one long-lived kernel](adr/0094-the-mcp-server-is-four-script-tools-over-one-kernel.md) | `nyx-script::{kernel,host,stubs,prelude}`, `nyx-mcp::main`, `nyx-engine::{guide,cards}` |
| 95 | [The queue is a folder of the workspace, run by a queue worker that outlives the window](adr/0095-the-queue-is-a-folder-of-the-workspace-run-by-a-queue-worker.md) | `nyx-engine::{shared_queue,worker,queue,engine}`, `nyx-mcp::main`, `app::main` |
| 96 | [A sweep may start each point from the modes of the one before](adr/0096-a-sweep-may-start-each-point-from-the-modes-of-the-one-before.md) | `femcore::warm`, `femcore::linalg::{driver,arnoldi_real,eps}`, `nyx-engine::{reuse,sweep,engine,worker}` |
| 97 | [The MCP server reads while it waits, and says what a job costs](adr/0097-the-mcp-server-reads-while-it-waits-and-says-what-a-job-costs.md) | `nyx-mcp`, `nyx-script::kernel`, `nyx-engine::{convergence,preview}`, `femcore::meshing::sizing` |
| 98 | [Connect knows Antigravity and gives an AI the instructions](adr/0098-connect-knows-antigravity-and-gives-an-ai-the-instructions.md) | `nyx-connect::{hosts,server}`, `nyx-connect` command, the app's Connect menu |
| 99 | [Figures of modes, and one node order in the results file](adr/0099-figures-of-modes-and-one-node-order-in-the-results-file.md) | `nyx-analysis::{section,modeview,canvas,raster,figure}`, `nyx-engine::{figure,h5io,xdmf}` |
| 100 | [nyx is dual licensed: AGPL-3.0 or a commercial licence](adr/0100-nyx-is-dual-licensed-agpl-3-or-commercial.md) | `LICENSE`, `NOTICE`, the `license` of the crates, the About sheet |
| 101 | [Solid waveguide modes, Floquet-periodic absorbing layers, and lossless convergence floor](adr/0101-sms-waveguide-modes-and-floquet-pml.md) | `nyx-sms::{pml,waveguide,cross_section,problem,lib,solid}`, `nyx-engine::convergence`, `materials/` |
| 102 | [Conformal interior 2D sheets (`Sheet`) and batch polygons (`Polygons`)](adr/0102-conformal-interior-sheets-and-batch-polygons.md) | `geom::poly_batch`, `Op::{Polygons,Sheet}`, `femcore::meshing::sheet3d`, `nyx-engine::{geomdesc,reduce}` |
| 103 | [Multi-RHS capacitance matrices and gauge-stabilized magnetostatics](adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md) | `nyx-dfs::{electro,magneto}`, `nyx-analysis::capacitance`, `nyx-physics::cqed` |
| 104 | [BCS superconducting surface impedance, thin-film IBC, and rank-1 lumped Josephson ports](adr/0104-bcs-surface-impedance-and-lumped-josephson-ports.md) | `femcore::material::bcs`, `materials/*.toml`, `nyx-ems::{ibc,rom,cavity}`, `nyx-engine::materials` |
| 105 | [KICQ-EPR, hierarchical Craig–Bampton quantization, and edge-regularized TLS participation](adr/0105-kicq-epr-craig-bampton-quantization-and-tls-participation.md) | `nyx-analysis::{epr,hamiltonian,tls}`, `nyx-engine::analyse` |
| 106 | [Driven frequency-domain solver, waveports, and Hermite Tangential Galerkin PROM](adr/0106-driven-frequency-domain-waveports-and-hermite-prom.md) | `nyx-ems::{driven,ports}`, `nyx-physics::cqed::DrivenPortData` |
| 107 | [Pure-Rust HDF5 storage (`hdf5-pure` + `zlib-rs`) and packed 5D Shuffle+LZF layout](adr/0107-pure-rust-hdf5-storage-and-packed-lzf-layout.md) | `nyx-storage`, `nyx-engine::h5io`, `Cargo.toml` |
| 108 | [Multi-source material library and on-demand database importers](adr/0108-multi-source-material-library-and-database-importers.md) | `materials/*`, `femcore::material::import`, `nyx-engine::material_fetch`, `nyx material` CLI |
| 109 | [Material `format = 2`: multi-physics constitutive models, Mattis–Bardeen BCS surface impedance, and coupled-field tensors](adr/0109-material-format-v2-multiphysics-bcs-coupled-tensors.md) | `femcore::material::{bcs,coupled,state_t,spec,schema,reader,optics,import}`, `nyx-ems::ibc`, `nyx-sms::solid`, `nyx-physics::decode`, `nyx-engine::{materials,material_fetch}`, `app::{matcmd,tabs}` |
| 110 | [Periodic 3D lattice conformance, per-layer region `sweep`, and fixed-frequency 3D Bloch modes with dispersive Poynting group index](adr/0110-periodic-3d-lattice-conformance-and-bloch-poynting-modes.md) | `femcore::linalg::driver`, `femcore::meshing::{mesher3d,refine3d,geomdrive3,mesher,curve3}`, `nyx-ems::{bloch,modes,modes_complex,cavity,problem,lib}`, `nyx-engine::{project,meshing,preview}` |
| 111 | [Un-averaged P2 3D field figures (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`) and zero-copy Bloch tiling](adr/0111-3d-field-figures-and-bloch-unit-cell-tiling.md) | `nyx-analysis::{view3d,figure,raster}`, `nyx-engine::{figure,commands}` |
| 112 | [Single-pass bit-exact batch reductions, mode-signature convergence tracking, and SPD Richardson extrapolation](adr/0112-single-pass-batch-reductions-and-convergence-tracking-extrapolation.md) | `nyx-engine::{reduce,convergence,validated,engine,commands,commandlist}` |
| 113 | [Script per-turn suspension resilience, UTF-8 BOM tolerance, inline tensor properties, and multi-session MCP](adr/0113-script-budget-resilience-multi-session-mcp-and-inline-tensors.md) | `nyx-script::{kernel,host,prelude,trace}`, `nyx-mcp::main`, `nyx-physics::{spec,validate,decode,schema}`, `femcore::material::spec`, `nyx-engine::{guide,cards,project,commandlist}` |
| 114 | [Pure-Rust binary GDSII `import_gds`, `sheet3d` PSLG curved-hole guard, coplanar `gap_cells`, and pure-tet `z_coarsen`](adr/0114-gdsii-sheet-import-coplanar-gap-sizing-and-z-coarsening.md) | `geom::gds`, `femcore::meshing::{sheet3d,geomdrive,request}`, `nyx-engine::{gds,project,meshing,scriptout,commands}` |
| 115 | [Vertical-line P1-iso-P2 Block-CG AMG, quarter-point sheet edges, layered `planar_bem`, and cQED Schur/qubit reduction](adr/0115-line-amg-block-cg-quarter-point-sheets-planar-bem-and-cqed-reduction.md) | `femcore::{fem,linalg}`, `nyx-dfs::{planar_bem,electro,problem,lib}`, `nyx-analysis::capacitance`, `nyx-engine::commands` |
| 116 | [Geometry BVH acceleration, polymorphic selectors, stepped airbridges, sheet3d Z-slab/max_cells guards, Ruppert diagnostics, and curved sheet projection](adr/0116-geom-bvh-meshing-guards-curved-sheets-and-airbridges.md) | `geom::{bvh,complex,arrange,stack,poly_batch,select,scene,model,planar}`, `femcore::meshing::{sheet3d,mesher,request,geomdrive,adapt3d,pml_geom3,curve3}`, `nyx-engine::{project,meshing,commands,scriptout,convergence}`, `app::app` |
| 117 | [Indefinite supernodal `LBL^H`/`Sym LBL^T` direct solver and $P_2$ Nédélec Hiptmair–Xu AMS](adr/0117-indefinite-lblt-direct-solver-and-p2-nedelec-ams.md) | `femcore::linalg::{direct,driver,ams,maxwell,sweep,blr,mixed,sparse_lu}` |
| 118 | [cQED open Robin BC, cavity lumped ports, driven PROM/waveports, unified `NetworkAnalysis`, DC-SQUID sweeps, 3D metal-thickness capacitance, and TOML TLS](adr/0118-cqed-cavity-driven-open-robin-network-bbq-and-squid.md) | `nyx-physics::{cqed,decode}`, `nyx-dfs::{electro,problem,lib}`, `nyx-ems::{driven,ports,ibc,cavity,rom,problem,lib}`, `nyx-analysis::{network,epr,hamiltonian,capacitance,tls}`, `nyx-engine::{analyse,commands,cards}` |
| 119 | [Field-less worker HDF5 I/O, per-mode read budgets, memory admission timeouts, and sub-picofarad linear convergence checks](adr/0119-fieldless-worker-io-admission-timeout-and-sub-pf-convergence.md) | `nyx-storage`, `nyx-engine::{h5io,xdmf,store_link,estimate,jobs,worker,queue,shared_queue,engine,sweep,commands,convergence}`, `app::runfile` |
