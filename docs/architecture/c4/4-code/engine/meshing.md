# Meshing seam: code

C4 level 4 (code) of the component **Meshing seam** (`engine.meshing`) in the container **nyx-engine** (crate `nyx_engine`). [Component view of the container](../../3-components/engine.md) · [All code pages](../README.md)

> MeshRequest to MeshedModel (mesh, names, label resolver); the only caller of the femcore mesher; absorbing layers across the two halves of a project; trial meshes and edit checks for the feedback while editing (ADR 92).

**Decisions:** [ADR 54](../../../adr/0054-mesh-request-and-report-as-data.md), [ADR 50](../../../adr/0050-wavelength-sizing-in-3d.md), [ADR 92](../../../adr/0092-feedback-while-editing-is-a-value-never-a-job.md)  
**Size:** 4 files, 1616 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-engine/src/amr.rs` | 335 |  |
| `crates/nyx-engine/src/meshing.rs` | 559 |  |
| `crates/nyx-engine/src/pml.rs` | 533 |  |
| `crates/nyx-engine/src/preview.rs` | 189 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `analysis` (nyx-analysis)](../../3-components/analysis.md) | modes, polarization | ModeAnalysis |
| [container `femcore` (femcore)](../../3-components/femcore.md) | adapt::electro, material, material::spec, mesh, mesh::cells, meshing, meshing::geomdrive, meshing::geomdrive3 +8 | Adaptive, CANCELLED, CellType, ElectroMetric2dOptions, ElectroMetric3dOptions, Layer, Library, M6, MeshReport, MeshRequest +31 |
| [container `geom` (geom)](../../3-components/geom.md) | eval, labels, model, select, sim | EntityKind, EntitySet, Evaluated, FeatureList, LabelMap, PmlCoordKind, PmlGrowMode, PmlRegion, PmlRole, Selector |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | cqed, crate root, error, physics, registry, result, spec | AdaptPassRecord, BcSpec, CAPACITANCE_KEY, CapacitanceMatrixData, ElectroAmrData, Error, LabelResolver, MaterialAssign, Mesh, ProblemSpec +5 |
| [Commands](./commands.md) | engine | Engine |
| [Job runner and queue](./jobs.md) | estimate, solve | Estimate, History, chain_key, chain_key_for_project, solve_one_single |
| [Project](./project.md) | crate root, project | MeshSettings, PmlModeRequest, Project, SweepLayer |
| [Registry and facade](./registry.md) | crate root | physics_for |

## Used by

| Component | Modules | Items |
|---|---|---|
| [app: Window](../app/window.md) | meshing, preview | MeshedModel, SizePreview |
| [Analyses](./analyse.md) | meshing | MeshedModel |
| [Commands](./commands.md) | amr, meshing | FemcoreMesher, MeshedModel, Mesher, project_uses_electro_amr |
| [Job runner and queue](./jobs.md) | amr, meshing, preview | FemcoreMesher, MeshedModel, Mesher, cost_of_project, project_uses_electro_amr, solve_dfs_with_amr |
| [Registry and facade](./registry.md) | meshing | FemcoreMesher, MeshedModel, Mesher |
| [Results files](./resultsIO.md) | amr, meshing | MeshedModel, Mesher, project_has_electro_amr |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace amr_rs {
    class fns_amr_rs {
        <<module amr.rs>>
        +project_uses_electro_amr(project, targets) bool
        +project_has_electro_amr(project) bool
        ~solve_dfs_with_amr(registry, project, initial_model, library, progress, done) Result
    }
  }
  namespace meshing_rs {
    class MeshedModel {
        <<struct>>
        +from_exchange(mesh, eval, labels, swept) Result
        +new(mesh, eval, labels, swept) Result
        ~with_report(report) MeshedModel
        ~with_amr_adapted() MeshedModel
        +is_amr_adapted() bool
        +mesh() Arc~Mesh~
        +solver_mesh() Arc~SimplexMesh~
        +evaluated() Arc~Evaluated~
        +more_7_methods()
    }
    class Mesher {
        <<interface>>
        +mesh(geometry, settings, sizing) Result
    }
    class FemcoreMesher {
        <<struct>>
    }
    class fns_meshing_rs {
        <<module meshing.rs>>
        ~build_mesh_request(s, sizing) Result
    }
  }
  namespace pml_rs {
    class Engine {
        <<impl block>>
        +add_pml_conditions(physics, strength) Result
        +add_pml_conditions_with(physics, strength, profile, auto_scale, target_reflection) Result
        +pml_diagnostics() Result
        +default_targets() Vec~String~
        ~cost_json(targets, dim, cells) Value
        +check_edits(geometry, mesh) Result
        +size_preview(geometry, mesh) Result
        +polarization_task(physics) Option
    }
  }
  namespace preview_rs {
    class SizePreview {
        <<struct>>
        +dim u8
        +cells u64
        +mesh_seconds f64
        +chain String
        +estimate Option~Estimate~
    }
    class SizeTask {
        <<struct>>
        +run() Result
    }
    class fns_preview_rs {
        <<module preview.rs>>
        ~cost_of_chain(history, chain, dim, cells) Value
        ~cost_of(history, targets, dim, cells) Value
        ~cost_of_project(history, project, targets, dim, cells) Value
    }
  }
  FemcoreMesher ..|> Mesher
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `project_uses_electro_amr` | `amr.rs:66` | `pub fn project_uses_electro_amr(project: &Project, targets: &[String]) -> bool` —  True when `project` and `targets` request multi-pass electrostatic AMR (`adapt_passes > 1`) on an unswept mesh. |
| fn | `project_has_electro_amr` | `amr.rs:74` | `pub fn project_has_electro_amr(project: &Project) -> bool` —  True when `project` carries a `dfs` `electrostatic` study configured with `adapt_passes > 1`. |
| fn | `solve_dfs_with_amr` | `amr.rs:132` | `pub(crate) fn solve_dfs_with_amr( registry: &Registry, project: &Project, initial_model: Arc<MeshedModel>, library: &Arc<femcore::material::Library>, progress: &Arc<Progress>, done: &[(String, SolveResult)], ) -> Result<(Arc<MeshedModel>, SolveResult), Error>` —  Runs the multi-pass goal-oriented AMR loop for `dfs` `electrostatic` and returns `(final_model, final_result)`. |
| struct | `MeshedModel` | `meshing.rs:17` |  |
| trait | `Mesher` | `meshing.rs:271` | methods: `mesh` —  Builds the mesh of a model. `sizing` is what the physics of the project want resolved (`MeshSettings::physics_sizing`); the engine |
| struct | `FemcoreMesher` | `meshing.rs:291` |  The mesher of `femcore`, driven by a `MeshRequest` (ADR 54): 2D and 3D models, simple or adaptive, inline or structured absorbing |
| fn | `build_mesh_request` | `meshing.rs:327` | `pub(crate) fn build_mesh_request( s: &MeshSettings, sizing: Option<&Sizing>, ) -> Result<femcore::meshing::request::MeshRequest, Error>` —  Validates `MeshSettings` and builds the corresponding [`femcore::meshing::request::MeshRequest`]. |
| fn | `MeshedModel::from_exchange` | `meshing.rs:36` | `pub fn from_exchange( mesh: Mesh, eval: Arc<Evaluated>, labels: LabelMap, swept: bool, ) -> Result<MeshedModel, Error>` —  From an uncut mixed-cell exchange [`Mesh`] (`Quad`, `Prism`, `Hex` preserved) without eager simplex splitting. |
| fn | `MeshedModel::new` | `meshing.rs:95` | `pub fn new( mesh: femcore::mesh::SimplexMesh, eval: Arc<Evaluated>, labels: LabelMap, swept: bool, ) -> Result<MeshedModel, Error>` —  From the mesh the mesher made (the solver mesh of femcore) to the exchange mesh the physics are given. |
| fn | `MeshedModel::with_report` | `meshing.rs:108` | `pub(crate) fn with_report(mut self, report: femcore::meshing::request::MeshReport) -> MeshedModel` —  The same, with the report the mesher made of it. |
| fn | `MeshedModel::with_amr_adapted` | `meshing.rs:114` | `pub(crate) fn with_amr_adapted(mut self) -> MeshedModel` —  Marks this model as produced by multi-pass AMR (`AF-1`). |
| fn | `MeshedModel::is_amr_adapted` | `meshing.rs:120` | `pub fn is_amr_adapted(&self) -> bool` —  True when this model was adapted during multi-pass AMR rather than built solely from `MeshSettings`. |
| fn | `MeshedModel::mesh` | `meshing.rs:124` | `pub fn mesh(&self) -> &Arc<Mesh>` |
| fn | `MeshedModel::solver_mesh` | `meshing.rs:130` | `pub fn solver_mesh(&self) -> Arc<femcore::mesh::SimplexMesh>` —  The mesh in the numbering the plots read (vertices, edges, P2 dofs): the exchange mesh cut into simplices, as every solve sees it. |
| fn | `MeshedModel::evaluated` | `meshing.rs:137` | `pub fn evaluated(&self) -> &Arc<Evaluated>` —  The geometry the mesh was made from, evaluated: what a front-end needs to name an entity it points at (selector suggestions). |
| fn | `MeshedModel::labels` | `meshing.rs:142` | `pub fn labels(&self) -> &LabelMap` —  The names of the model: region and boundary numbers. |
| fn | `MeshedModel::swept_ends` | `meshing.rs:147` | `pub fn swept_ends(&self) -> Option<(i32, i32)>` —  The boundary labels of the two end faces of a model swept along y. |
| fn | `MeshedModel::report` | `meshing.rs:151` | `pub fn report(&self) -> Option<&femcore::meshing::request::MeshReport>` |
| fn | `MeshedModel::regions` | `meshing.rs:156` | `pub fn regions(&self) -> Vec<(i32, String)>` —  Region labels with their names. |
| fn | `MeshedModel::info` | `meshing.rs:161` | `pub fn info(&self) -> Value` —  What a person or a host needs to choose names: the regions, the boundary names and the size of the mesh. |
| fn | `MeshedModel::suggest` | `meshing.rs:528` | `pub fn suggest(&self, labels: &[i32]) -> Vec<String>` —  Selectors, simplest first, that pick exactly the boundaries with these numbers; empty when no selector does. What the window offers |
| fn | `MeshedModel::boundary_names` | `meshing.rs:555` | `pub fn boundary_names(&self) -> Vec<(i32, String)>` —  The names of the boundaries with their numbers (a boundary name may stand for several numbers). |
| fn | `Engine::add_pml_conditions` | `pml.rs:135` | `pub fn add_pml_conditions(&mut self, physics: &str, strength: f64) -> Result<Value, Error>` —  Adds a `pml` (or `infinite_element` for static studies) condition for every absorbing piece of the geometry to `physics` |
| fn | `Engine::add_pml_conditions_with` | `pml.rs:142` | `pub fn add_pml_conditions_with( &mut self, physics: &str, strength: Option<f64>, profile: Option<&str>, auto_scale: Option<bool>, target_reflection: Option<f64>, ) -> Result<Value, Error>` —  Full absorbing-condition generator supporting explicit or auto-calibrated `strength`, `profile` (`"auto"`, `"kappa_polynomial"`, |
| fn | `Engine::pml_diagnostics` | `pml.rs:369` | `pub fn pml_diagnostics(&mut self) -> Result<Value, Error>` —  Pre-solve and post-mesh absorbing-layer diagnostics: uncut element counts, minimum PML scaled Jacobian, |
| fn | `Engine::default_targets` | `pml.rs:511` | `pub fn default_targets(&self) -> Vec<String>` —  The physics solved when none is named: the active one and, before it, those of the project that export what it imports (a heat |
| struct | `SizePreview` | `preview.rs:31` | fields: `dim`, `cells`, `mesh_seconds`, `chain`, `estimate` —  The size and cost of a model that is being edited. |
| fn | `cost_of_chain` | `preview.rs:51` | `pub(crate) fn cost_of_chain(history: &History, chain: String, dim: u8, cells: u64) -> serde_json::Value` —  What a solve of `targets` on a mesh of `cells` cells costs, from the solves `history` has seen (the answer of `mesh_build` and `convergence_check` |
| fn | `cost_of` | `preview.rs:64` | `pub(crate) fn cost_of(history: &History, targets: &[String], dim: u8, cells: u64) -> serde_json::Value` |
| fn | `cost_of_project` | `preview.rs:68` | `pub(crate) fn cost_of_project( history: &History, project: &Project, targets: &[String], dim: u8, cells: u64, ) -> serde_json::Value` |
| struct | `SizeTask` | `preview.rs:88` |  A trial mesh waiting to run. It owns copies of everything it needs and is `Send`, so the caller chooses the thread. |
| fn | `SizeTask::run` | `preview.rs:101` | `pub fn run(self) -> Result<SizePreview, Error>` —  Meshes the copy of the edits and looks up the estimate. Blocks for as long as the mesher takes (milliseconds for the models this is |
| fn | `Engine::cost_json` | `preview.rs:117` | `pub(crate) fn cost_json(&self, targets: &[String], dim: u8, cells: u64) -> serde_json::Value` —  What a solve of `targets` on a mesh of `cells` cells costs, as the answer of `mesh_build` and `convergence_check` shows it: the estimate from the |
| fn | `Engine::check_edits` | `preview.rs:123` | `pub fn check_edits(&self, geometry: &FeatureList, mesh: &MeshSettings) -> Result<(), Error>` —  The first problem of this project with `geometry` and `mesh` in place of its own, in plain words, or `Ok` when it passes |
| fn | `Engine::size_preview` | `preview.rs:133` | `pub fn size_preview(&self, geometry: FeatureList, mesh: MeshSettings) -> Result<SizeTask, Error>` —  A trial mesh of the project with `geometry` and `mesh` applied, as a task that can run on any thread. Refused (with the reason, in plain |
| fn | `Engine::polarization_task` | `preview.rs:161` | `pub fn polarization_task( &self, physics: &str, ) -> Option<(Arc<SolveResult>, impl FnOnce() -> Result<ModeAnalysis, Error> + Send + 'static)>` —  The polarization analysis of the latest result of `physics` as a function that does not borrow the engine, so it can run beside the |

Trait implementations in this component: `MeshedModel: LabelResolver`, `FemcoreMesher: Mesher`
