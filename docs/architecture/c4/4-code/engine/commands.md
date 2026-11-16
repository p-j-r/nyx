# Commands: code

C4 level 4 (code) of the component **Commands** (`engine.commands`) in the container **nyx-engine** (crate `nyx_engine`). [Component view of the container](../../3-components/engine.md) · [All code pages](../README.md)

> The Engine and its catalog: every action by name with JSON arguments, JSON Schema and JSON reply; the dispatcher; the guide cards a host reads on demand. The window and the MCP server use the same entries.

**Decisions:** [ADR 57](../../../adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md), [ADR 38](../../../adr/0038-geometry-operations-are-described-once.md), [ADR 15](../../../adr/0015-ai-interface-separate-workspace-and-mcp-server.md)  
**Size:** 3 files, 1410 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-engine/src/commands.rs` | 696 | The commands: every action a front-end can take, by name, with JSON arguments, a JSON schema and a JSON reply. The window calls the typed methods of `Engine` directly or these commands; the MCP server turns `catalog()` into its tool list... |
| `crates/nyx-engine/src/engine.rs` | 681 | The engine: one project, the mesh built for it, the latest result of each physics, the jobs that produce them and the binding to a workflow folder. The window and the MCP server hold an `Engine` and nothing else; every action is a method... |
| `crates/nyx-engine/src/guide.rs` | 33 | Cards: short texts a host reads on demand through the `get_guide` command, embedded in the binary. They describe the commands, not a physics: how to build a geometry, assign materials, solve and read a result. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material::spec, mesh, progress | Library, SimplexMesh, install, stage |
| [container `geom` (geom)](../../3-components/geom.md) | model | FeatureList, Op |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | contract, crate root, error, registry, result, spec, validate | BcSpec, Error, MaterialAssign, ProblemSpec, Registry, Sizing, SolveResult, StudySpec, canonical_json, validate_spec |
| [container `store` (nyx-store)](../../3-components/store.md) | util | hash_bytes |
| [Geometry operations](./geometry.md) | geomdesc, geometry |  |
| [Job runner and queue](./jobs.md) | jobs, queue, solve, sweep | JobId, JobMeta, JobOutput, JobRunner, Lane, Opened, QueueOpts, SolveSummary, Stages, Sweep +2 |
| [Materials](./materials.md) | materials | default_dirs, load_library |
| [Meshing seam](./meshing.md) | meshing | FemcoreMesher, MeshedModel, Mesher |
| [Project](./project.md) | project | MeshSettings, Project |
| [Registry and facade](./registry.md) | crate root | physics_for, registry |
| [Results files](./resultsIO.md) | results_file | read_and_mesh, write_snapshot |
| [Store link](./storeLink.md) | store_link | RunOutcome, StoreLink |

## Used by

| Component | Modules | Items |
|---|---|---|
| [app: Scripts and batch](../app/headless.md) | engine | Engine |
| [app: Shell and dashboard](../app/shell.md) | engine | Engine |
| [app: Window](../app/window.md) | engine | Engine |
| [Analyses](./analyse.md) | engine | Engine |
| [Command lists](./commandlist.md) | commands, engine | Access, Engine |
| [Job runner and queue](./jobs.md) | commands, engine | Access, Engine, Purpose |
| [Materials](./materials.md) | engine | Engine |
| [Meshing seam](./meshing.md) | engine | Engine |
| [Registry and facade](./registry.md) | commands, engine | Access, CommandDesc, Engine |
| [Results files](./resultsIO.md) | engine | Engine |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace commands_rs {
    class Access {
        <<enumeration>>
        Read
        Edit
        Expensive
    }
    class CommandDesc {
        <<struct>>
        +name String
        +title String
        +description String
        +access Access
        +schema Value
    }
    class Engine {
        <<impl block>>
        +registry Registry
        +mesher Arc
        +project Project
        +model Option
        +results BTreeMap
        +jobs JobRunner
        +pending BTreeMap
        +meta BTreeMap
        +more_10_fields
        +catalog() Vec~CommandDesc~
        +call(name, args) R
        +new() Engine
        +with(registry, mesher) Engine
        +registry() Registry
        +set_mesher(mesher)
        +project() Project
        +jobs() JobRunner
        +more_27_methods()
    }
  }
  namespace engine_rs {
    class Stored {
        <<struct>>
        +result Arc~SolveResult~
        +key String
        +job JobId
        +at Instant
    }
    class Purpose {
        <<enumeration>>
        Current
        Point
    }
    class Pending {
        <<enumeration>>
        Mesh
        Solve
        Point
        Convergence
        Open
        Export
    }
  }
  namespace guide_rs {
    class fns_guide_rs {
        <<module guide.rs>>
        +card(topic) Option
        +TOPICS const
        +SERVER_INSTRUCTIONS const
    }
  }
  CommandDesc --> Access
  Engine --> Stored
  Engine --> Pending
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `Access` | `commands.rs:17` | variants: `Read`, `Edit`, `Expensive` —  How much a command changes or costs; a host may ask the person before an expensive or editing command. |
| struct | `CommandDesc` | `commands.rs:27` | fields: `name`, `title`, `description`, `access`, `schema` |
| fn | `Engine::catalog` | `commands.rs:105` | `pub fn catalog(&self) -> Vec<CommandDesc>` —  Every command of this build, with its schema. The geometry `add_*` commands come from `geomdesc`; the choices of physics from |
| fn | `Engine::call` | `commands.rs:229` | `pub fn call(&mut self, name: &str, args: &Value) -> R` —  Runs the command `name` with JSON arguments and returns its JSON reply. |
| struct | `Stored` | `engine.rs:22` | fields: `result`, `key`, `job`, `at` |
| enum | `Purpose` | `engine.rs:34` | variants: `Current`, `Point` —  What a solve job is for: the engine's own project, or one point of a sweep. |
| enum | `Pending` | `engine.rs:39` | variants: `Mesh`, `Solve`, `Point`, `Convergence`, `Open`, `Export` |
| struct | `Engine` | `engine.rs:61` | fields: `registry`, `mesher`, `project`, `model`, `results`, `jobs`, `pending`, `meta`, `stages`, `owner`, `sweeps`, `next_sweep`, ... |
| fn | `Engine::new` | `engine.rs:96` | `pub fn new() -> Engine` —  An engine with the physics compiled into this build and an empty 2D project. |
| fn | `Engine::with` | `engine.rs:104` | `pub fn with(registry: Registry, mesher: Arc<dyn Mesher>) -> Engine` |
| fn | `Engine::registry` | `engine.rs:127` | `pub fn registry(&self) -> &Registry` |
| fn | `Engine::set_mesher` | `engine.rs:133` | `pub fn set_mesher(&mut self, mesher: Arc<dyn Mesher>)` —  Replaces the mesher (a front-end with its own mesher, or a test that counts or holds its calls). The mesh the engine holds is kept; the |
| fn | `Engine::project` | `engine.rs:137` | `pub fn project(&self) -> &Project` |
| fn | `Engine::jobs` | `engine.rs:141` | `pub fn jobs(&self) -> &JobRunner` |
| fn | `Engine::set_workspace` | `engine.rs:147` | `pub fn set_workspace(&mut self, root: PathBuf, auto_store: bool)` —  Names the workspace folder workflows are stored in. With `auto_store` the project is bound to a new workflow there at the first solve and |
| fn | `Engine::sizing` | `engine.rs:154` | `pub fn sizing(&self) -> Option<Sizing>` —  What the physics of the project want the mesh to resolve, merged: the shortest wavelength of each region, the highest order and the |
| fn | `Engine::mesh_key` | `engine.rs:174` | `pub fn mesh_key(&self) -> String` —  The hash of the geometry, the mesh settings and the sizing: what a mesh depends on. |
| fn | `Engine::stamp` | `engine.rs:181` | `pub fn stamp(&self) -> String` —  The hash of the whole project. Staleness no longer uses it (`problem_key` does); it is kept for tools that want an identity of the whole |
| fn | `Engine::problem_key` | `engine.rs:189` | `pub fn problem_key(&self, physics: &str) -> Option<String>` —  What a result of `physics` depends on, as a cache key: the problem of that physics, the materials the project carries, the mesh key and the |
| fn | `Engine::problem_keys` | `engine.rs:195` | `pub fn problem_keys(&self) -> BTreeMap<String, String>` —  The problem key of every physics of the project, with the mesh key computed once. |
| fn | `Engine::is_stale` | `engine.rs:233` | `pub fn is_stale(&self, physics: &str) -> bool` —  True when the result stored for `physics` was not made from the problem the project has now (or the project has no such problem any more). |
| fn | `Engine::result_age` | `engine.rs:241` | `pub fn result_age(&self, physics: &str) -> Option<f64>` —  Seconds since the result of `physics` was stored. |
| fn | `Engine::set_project` | `engine.rs:246` | `pub fn set_project(&mut self, project: Project)` —  Replaces the project (a new one, or one just opened). The mesh and results are dropped and the store binding is released. |
| fn | `Engine::model` | `engine.rs:256` | `pub fn model(&self) -> Option<&Arc<MeshedModel>>` —  The mesh if one has been built for the current geometry and settings. |
| fn | `Engine::result` | `engine.rs:261` | `pub fn result(&self, physics: &str) -> Option<(&Arc<SolveResult>, bool)>` —  The latest result of a physics and whether what it depends on has changed since (its problem key, not the whole project). |
| fn | `Engine::start_mesh` | `engine.rs:269` | `pub fn start_mesh(&mut self) -> Result<JobId, Error>` —  Starts meshing the project on a job (queued when the heavy lane is busy). |
| fn | `Engine::start_mesh_with` | `engine.rs:274` | `pub fn start_mesh_with(&mut self, opts: &QueueOpts) -> Result<JobId, Error>` —  `start_mesh` with queue options. |
| fn | `Engine::start_solve` | `engine.rs:302` | `pub fn start_solve(&mut self, targets: Option<Vec<String>>) -> Result<JobId, Error>` —  Starts solving `targets` (default: the active physics, else all) on a job. The mesh is built inside the job when none is current. |
| fn | `Engine::start_solve_with` | `engine.rs:307` | `pub fn start_solve_with(&mut self, targets: Option<Vec<String>>, opts: &QueueOpts) -> Result<JobId, Error>` —  `start_solve` with queue options. The job waits in the queue while the heavy lane is busy; it carries a frozen copy of the project. |
| fn | `Engine::submit_solve` | `engine.rs:313` | `pub(crate) fn submit_solve( &mut self, src: Option<&Engine>, targets: Option<Vec<String>>, opts: &QueueOpts, purpose: Purpose, ) -> Result<JobId, Error>` —  The one place a solve is submitted. `src` is the engine whose project is solved (a point of a sweep is a scratch engine with the point's |
| fn | `Engine::start_convergence` | `engine.rs:409` | `pub fn start_convergence(&mut self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<JobId, Error>` —  Starts a mesh convergence check (two solves, the second on a finer mesh) on a job. Both solves run on the job's thread under its progress, |
| fn | `Engine::start_convergence_with` | `engine.rs:415` | `pub fn start_convergence_with( &mut self, physics: &str, factor: f64, tolerance: Option<f64>, opts: &QueueOpts, ) -> Result<JobId, Error>` —  `start_convergence` with queue options. When the solve of `physics` alone on the present mesh was done already by an earlier request, the |
| fn | `Engine::start_open_results` | `engine.rs:445` | `pub fn start_open_results(&mut self, path: &str) -> Result<JobId, Error>` —  Starts reading a results file, meshing its project again and checking the mesh against the stored one, on a job. The engine is not touched |
| fn | `Engine::start_export_results` | `engine.rs:460` | `pub fn start_export_results(&mut self, path: &str, parallel: bool) -> Result<JobId, Error>` —  Starts writing the results the engine holds now (the project, the mesh and the latest result of every physics) into an HDF5 file on a job. |
| fn | `Engine::poll` | `engine.rs:479` | `pub fn poll(&mut self) -> Vec<(JobId, Result<Value, Error>)>` —  Takes the finished jobs and applies them: the mesh and results become current, a solve is recorded in the run table when the project |
| fn | `Engine::wait` | `engine.rs:491` | `pub fn wait(&mut self, id: JobId) -> Result<Value, Error>` —  Waits for a job and applies it. |
| fn | `Engine::apply` | `engine.rs:506` | `pub(crate) fn apply(&mut self, id: JobId, outcome: Result<JobOutput, Error>) -> Result<Value, Error>` |
| fn | `Engine::mesh_now` | `engine.rs:633` | `pub fn mesh_now(&mut self) -> Result<Value, Error>` —  Mesh now (start and wait). |
| fn | `Engine::solve_now` | `engine.rs:639` | `pub fn solve_now(&mut self, targets: Option<Vec<String>>) -> Result<Value, Error>` —  Solve now (start and wait). |
| fn | `Engine::save_to_workspace` | `engine.rs:647` | `pub fn save_to_workspace(&mut self, workspace: Option<PathBuf>) -> Result<Value, Error>` —  Binds the project to a new workflow folder in `workspace` (or commits it, when already bound). |
| fn | `Engine::open_workflow` | `engine.rs:666` | `pub fn open_workflow(&mut self, folder: &std::path::Path, experiment: Option<&str>) -> Result<Value, Error>` —  Opens a workflow folder and makes its experiment the project. |
| const | `TOPICS` | `guide.rs:4` | `&[&str]` |
| fn | `card` | `guide.rs:6` | `pub fn card(topic: &str) -> Option<&'static str>` |
| const | `SERVER_INSTRUCTIONS` | `guide.rs:17` | `&str` —  The text a server sends with `initialize`: just enough to find the guide. |

Trait implementations in this component: `Engine: Default`
