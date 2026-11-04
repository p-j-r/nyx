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
| 8 | [Ported solver library merged, not wired](adr/0008-ported-solver-library-staged-not-wired.md) | `solvers` |
| 9 | [Every kernel is three stages](adr/0009-kernels-are-three-stages.md) | kernels, `solve` |
| 10 | [Golden files are the numerical oracle](adr/0010-golden-files-are-the-numerical-oracle.md) | tests |
| 11 | [Determinism by fixed chunks](adr/0011-determinism-by-fixed-chunks.md) | `par`, `solvers::par` |
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
| 37 | [The engine core is built before the kernels, with two seams](adr/0037-engine-core-before-the-kernels.md) | `nyx-engine` |
| 38 | [Geometry operations are described once](adr/0038-geometry-operations-are-described-once.md) | `nyx-engine::geomdesc` |
| 39 | [The store takes spec text and a card function](adr/0039-the-store-takes-spec-text-and-a-card-function.md) | `nyx-store`, `nyx-engine::store_link` |
