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
| 84 | [Results are stale per physics: the key is the problem key](adr/0084-results-are-stale-per-physics-the-key-is-the-problem-key.md) | `nyx-engine::engine`, `app::app` |
