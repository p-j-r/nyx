# Commands: code

C4 level 4 (code) of the component **Commands** (`engine.commands`) in the container **nyx-engine** (crate `nyx_engine`). [Component view of the container](../../3-components/engine.md) · [All code pages](../README.md)

> The Engine and its catalog: every action by name with JSON arguments, JSON Schema and JSON reply; the dispatcher; the guide cards a host reads on demand; the reductions of results (field handles, integrals, extrema, probes, overlaps). The window and the MCP server use the same entries.

**Decisions:** [ADR 57](../../../adr/0057-front-ends-are-thin-the-mcp-tools-are-engine-commands.md), [ADR 38](../../../adr/0038-geometry-operations-are-described-once.md), [ADR 15](../../../adr/0015-ai-interface-separate-workspace-and-mcp-server.md), [ADR 90](../../../adr/0090-whole-array-operations-are-engine-commands-over-a-typed-tree.md), [ADR 91](../../../adr/0091-validated-is-a-state-checked-against-the-problem-key.md), [ADR 93](../../../adr/0093-scripts-are-python-text-run-by-monty-in-a-worker-process.md), [ADR 94](../../../adr/0094-the-mcp-server-is-four-script-tools-over-one-kernel.md), [ADR 102](../../../adr/0102-conformal-interior-sheets-and-batch-polygons.md)  
**Size:** 7 files, 7318 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-engine/src/commands.rs` | 1285 |  |
| `crates/nyx-engine/src/engine.rs` | 1427 |  |
| `crates/nyx-engine/src/guide.rs` | 36 |  |
| `crates/nyx-engine/src/path_guard.rs` | 159 |  |
| `crates/nyx-engine/src/reduce.rs` | 3578 |  |
| `crates/nyx-engine/src/scriptout.rs` | 346 |  |
| `crates/nyx-engine/src/validated.rs` | 487 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem, fieldexpr, formula, linalg::csr::assembly, material, material::spec, mesh, par +3 | Anisotropy, C, Fun, Library, MAX_POWER, Named, Names, Op, RefElement, SimplexMesh +12 |
| [container `geom` (geom)](../../3-components/geom.md) | model | FeatureList, Op |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | contract, cqed, crate root, decode, error, physics, registry, result +2 | BcSpec, CapacitanceMatrixData, Ctx, DiagnosticSeverity, DiagnosticWarning, Error, LabelResolver, MaterialAssign, NumericalErrorBudget, ProblemSpec +8 |
| [container `store` (nyx-store)](../../3-components/store.md) | runs, util | Convergence, RunPatch, hash_bytes |
| [Geometry operations](./geometry.md) | geomdesc, geometry |  |
| [Job runner and queue](./jobs.md) | convergence, estimate, jobs, process, queue, reuse, shared_queue, solve +2 | Estimate, History, JobId, JobMeta, JobOutput, JobRunner, Lane, Limits, Matching, MemoryProbe +34 |
| [Materials](./materials.md) | materials | default_dirs, load_library, study_wavelength |
| [Meshing seam](./meshing.md) | amr, meshing | FemcoreMesher, MeshedModel, Mesher, project_uses_amr |
| [Project](./project.md) | crate root, project | Auto, MeshSettings, Project |
| [Registry and facade](./registry.md) | crate root | physics_for, registry |
| [Results files](./resultsIO.md) | results_file | read_and_mesh, write_snapshot |
| [Store link](./storeLink.md) | store_link | RunOutcome, StoreLink |

## Used by

| Component | Modules | Items |
|---|---|---|
| [app: Scripts and batch](../app/headless.md) | engine | Engine |
| [app: Inspector tabs](../app/inspector.md) | engine | Engine |
| [app: Shell and dashboard](../app/shell.md) | engine, validated | Engine, Validated, Validation |
| [app: Window](../app/window.md) | engine, validated | Engine, Validated |
| [Analyses](./analyse.md) | engine | Engine |
| [API reference](./apiref.md) | commands, engine | Access, CommandDesc, Engine |
| [Command lists](./commandlist.md) | commands, engine | Access, Engine |
| [Geometry operations](./geometry.md) | engine, path_guard | Engine, confine |
| [Job runner and queue](./jobs.md) | commands, engine, validated | Access, Check, Engine, Pending, Purpose, Stored |
| [Materials](./materials.md) | engine | Engine |
| [Meshing seam](./meshing.md) | engine | Engine |
| [Registry and facade](./registry.md) | commands, engine | Access, CommandDesc, Engine |
| [Results files](./resultsIO.md) | engine, validated | Check, Engine |
| [script: Script host](../script/host.md) | engine | Engine |
| [script: Script kernel](../script/kernel.md) | engine | Engine |
| [script: Type stub](../script/stubs.md) | commands, engine | CommandDesc, Engine |
| [script: Tracing of script functions](../script/tracing.md) | engine | Engine |

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
        +more_22_fields
        +catalog() Vec~CommandDesc~
        +call(name, args) R
        +new() Engine
        +with(registry, mesher) Engine
        +set_script(script)
        +set_path_root(root) Result
        +path_root() Option
        +set_history_file(path)
        +more_63_methods()
    }
  }
  namespace engine_rs {
    class Stored {
        <<struct>>
        +result Arc~SolveResult~
        +key String
        +job JobId
        +at Instant
        +model Arc~MeshedModel~
        +run Option~String~
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
        QmcSweep
        Convergence
        Open
        Export
    }
    class fns_engine_rs {
        <<module engine.rs>>
        ~hash_library(lib) String
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
  namespace path_guard_rs {
    class fns_path_guard_rs {
        <<module path_guard.rs>>
        ~absolute_root(root) Result
        +confine(root, requested) Result
        ~confine_args(root, command, args) Result
        +PATH_ARGS const
    }
  }
  namespace reduce_rs {
    class What {
        <<enumeration>>
        Field
        Coord
        Material
        Region
        Normal
    }
    class Spec {
        <<struct>>
        +physics String
        +solution usize
        +name String
        +what What
    }
    class Defaults {
        <<struct>>
    }
    class Prep {
        <<struct>>
    }
    class Namer {
        <<impl block>>
    }
    class fns_reduce_rs {
        <<module reduce.rs>>
        +BLOCK const
    }
  }
  namespace scriptout_rs {
    class Exported {
        <<struct>>
        +script String
        +form str
        +commands usize
    }
    class fns_scriptout_rs {
        <<module scriptout.rs>>
        +export(p) Result
    }
  }
  namespace validated_rs {
    class Check {
        <<struct>>
        +factor f64
        +tolerance f64
        +matching String
        +passed bool
        +worst_relative f64
        +worst_quantity String
        +quantities usize
        +cells_coarse Option~u64~
        +more_9_fields
        +observable() String
        ~from_answer(v) Option~Check~
        +to_json() Value
        +to_record() Convergence
    }
    class Validation {
        <<enumeration>>
        Exploratory
        Checked
        NotAffordable
        +as_str() str
    }
    class ValidationTier {
        <<enumeration>>
        BoundCertified
        GciVerified
        TwoGridChecked
        ResidualChecked
        Exploratory
        NotAffordable
        NonConvergent
        +as_str() str
        +weakness() u8
    }
    class Validated {
        <<struct>>
        +state Validation
        +tier ValidationTier
        +reason String
        +check Option~Check~
        +to_json() Value
    }
    class fns_validated_rs {
        <<module validated.rs>>
        +PLANNED_FACTOR const
    }
  }
  CommandDesc --> Access
  Engine --> Stored
  Engine --> Pending
  Spec --> What
  Validated --> Validation
  Validated --> ValidationTier
  Validated --> Check
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `Access` | `commands.rs:20` | variants: `Read`, `Edit`, `Expensive` —  How much a command changes or costs; a host may ask the person before an expensive or editing command. |
| struct | `CommandDesc` | `commands.rs:30` | fields: `name`, `title`, `description`, `access`, `schema` |
| fn | `Engine::catalog` | `commands.rs:131` | `pub fn catalog(&self) -> Vec<CommandDesc>` —  Every command of this build, with its schema. The geometry `add_*` commands come from `geomdesc`; the choices of physics from |
| fn | `Engine::call` | `commands.rs:487` | `pub fn call(&mut self, name: &str, args: &Value) -> R` —  Runs the command `name` with JSON arguments and returns its JSON reply. |
| struct | `Stored` | `engine.rs:24` | fields: `result`, `key`, `job`, `at`, `model`, `run` |
| enum | `Purpose` | `engine.rs:40` | variants: `Current`, `Point` —  What a solve job is for: the engine's own project, or one point of a sweep. |
| enum | `Pending` | `engine.rs:45` | variants: `Mesh`, `Solve`, `Point`, `QmcSweep`, `Convergence`, `Open`, `Export` |
| struct | `Engine` | `engine.rs:75` | fields: `registry`, `mesher`, `project`, `model`, `results`, `jobs`, `pending`, `meta`, `stages`, `owner`, `sweeps`, `next_sweep`, ... |
| fn | `hash_library` | `engine.rs:144` | `pub(crate) fn hash_library(lib: &femcore::material::Library) -> String` |
| fn | `Engine::new` | `engine.rs:156` | `pub fn new() -> Engine` —  An engine with the physics compiled into this build and an empty 2D project. |
| fn | `Engine::with` | `engine.rs:165` | `pub fn with(registry: Registry, mesher: Arc<dyn Mesher>) -> Engine` |
| fn | `Engine::set_script` | `engine.rs:206` | `pub fn set_script(&mut self, script: Option<(String, String)>)` —  Names the script that drives this engine (`name`, `text`): every run recorded while it is set keeps the text, with its hash in the file name, |
| fn | `Engine::set_path_root` | `engine.rs:213` | `pub fn set_path_root(&mut self, root: Option<&std::path::Path>) -> Result<(), String>` —  Keeps every file and folder argument of the commands inside `root` (S4, ADR [89]): a relative path is taken from `root`, a path that |
| fn | `Engine::path_root` | `engine.rs:222` | `pub fn path_root(&self) -> Option<&std::path::Path>` —  The folder the path guard keeps commands in, when it is on. |
| fn | `Engine::set_history_file` | `engine.rs:228` | `pub fn set_history_file(&mut self, path: &std::path::Path)` —  Keeps what solves cost in `path` (JSON lines, appended to) and learns from what it holds: the memory and time estimates of later solves |
| fn | `Engine::history` | `engine.rs:233` | `pub fn history(&self) -> &Arc<crate::estimate::History>` —  What the engine has learned about the cost of solves. |
| fn | `Engine::set_memory_probe` | `engine.rs:238` | `pub fn set_memory_probe(&mut self, probe: crate::estimate::MemoryProbe)` —  Replaces how the engine asks the machine for `(total, available)` memory in bytes (a test that plays a machine that is full). |
| fn | `Engine::set_memory_wait_timeout` | `engine.rs:244` | `pub fn set_memory_wait_timeout(&mut self, timeout_s: Option<f64>)` —  Sets the engine's caller-level default `memory_wait_timeout_s` (`None` waits indefinitely unless `NYX_MEMORY_WAIT_TIMEOUT_S` or an |
| fn | `Engine::memory_wait_timeout` | `engine.rs:249` | `pub fn memory_wait_timeout(&self) -> Option<f64>` —  Returns the engine's caller-level default `memory_wait_timeout_s`. |
| fn | `Engine::set_worker` | `engine.rs:257` | `pub fn set_worker(&mut self, worker: PathBuf, by_default: bool)` —  Names the worker program (`nyx-worker`, built from this crate) that solves run in when they are isolated: a solve in a child process |
| fn | `Engine::set_shared_queue` | `engine.rs:268` | `pub fn set_shared_queue(&mut self, folder: &std::path::Path, by_default: bool)` —  Names the folder whose `queue` this engine shares with every other program that names it (the window, the MCP server, `nyx queue`), and |
| fn | `Engine::set_self_as_worker` | `engine.rs:277` | `pub fn set_self_as_worker(&mut self, by_default: bool) -> bool` —  Makes the running program its own worker: a front-end (`nyx`, `nyx-mcp`) that answers `--job FOLDER` with `worker::serve` calls this, so |
| fn | `Engine::base_library` | `engine.rs:287` | `pub(crate) fn base_library(&self) -> &Arc<femcore::material::Library>` |
| fn | `Engine::set_base_library_arc` | `engine.rs:291` | `pub(crate) fn set_base_library_arc(&mut self, library: Arc<femcore::material::Library>)` |
| fn | `Engine::registry` | `engine.rs:296` | `pub fn registry(&self) -> &Registry` |
| fn | `Engine::set_mesher` | `engine.rs:302` | `pub fn set_mesher(&mut self, mesher: Arc<dyn Mesher>)` —  Replaces the mesher (a front-end with its own mesher, or a test that counts or holds its calls). The mesh the engine holds is kept; the |
| fn | `Engine::project` | `engine.rs:307` | `pub fn project(&self) -> &Project` |
| fn | `Engine::jobs` | `engine.rs:311` | `pub fn jobs(&self) -> &JobRunner` |
| fn | `Engine::set_workspace` | `engine.rs:317` | `pub fn set_workspace(&mut self, root: PathBuf, auto_store: bool)` —  Names the workspace folder workflows are stored in. With `auto_store` the project is bound to a new workflow there at the first solve and |
| fn | `Engine::sizing` | `engine.rs:324` | `pub fn sizing(&self) -> Option<Sizing>` —  What the physics of the project want the mesh to resolve, merged: the shortest wavelength of each region, the highest order and the |
| fn | `Engine::mesh_key` | `engine.rs:344` | `pub fn mesh_key(&self) -> String` —  The hash of the geometry, the mesh settings and the sizing: what a mesh depends on. |
| fn | `Engine::stamp` | `engine.rs:351` | `pub fn stamp(&self) -> String` —  The hash of the whole project. Staleness no longer uses it (`problem_key` does); it is kept for tools that want an identity of the whole |
| fn | `Engine::problem_key` | `engine.rs:359` | `pub fn problem_key(&self, physics: &str) -> Option<String>` —  What a result of `physics` depends on, as a cache key: the problem of that physics, the materials the project carries, the mesh key and the |
| fn | `Engine::problem_keys` | `engine.rs:365` | `pub fn problem_keys(&self) -> BTreeMap<String, String>` —  The problem key of every physics of the project, with the mesh key computed once. |
| fn | `Engine::is_stale` | `engine.rs:413` | `pub fn is_stale(&self, physics: &str) -> bool` —  True when the result stored for `physics` was not made from the problem the project has now (or the project has no such problem any more). |
| fn | `Engine::result_age` | `engine.rs:421` | `pub fn result_age(&self, physics: &str) -> Option<f64>` —  Seconds since the result of `physics` was stored. |
| fn | `Engine::set_project` | `engine.rs:426` | `pub fn set_project(&mut self, project: Project)` —  Replaces the project (a new one, or one just opened). The mesh, results, checks and store binding are cleared. |
| fn | `Engine::model` | `engine.rs:437` | `pub fn model(&self) -> Option<&Arc<MeshedModel>>` —  The mesh if one has been built for the current geometry and settings. |
| fn | `Engine::result` | `engine.rs:442` | `pub fn result(&self, physics: &str) -> Option<(&Arc<SolveResult>, bool)>` —  The latest result of a physics and whether what it depends on has changed since (its problem key, not the whole project). |
| fn | `Engine::start_mesh` | `engine.rs:450` | `pub fn start_mesh(&mut self) -> Result<JobId, Error>` —  Starts meshing the project on a job (queued when the heavy lane is busy). |
| fn | `Engine::start_mesh_with` | `engine.rs:455` | `pub fn start_mesh_with(&mut self, opts: &QueueOpts) -> Result<JobId, Error>` —  `start_mesh` with queue options. |
| fn | `Engine::start_solve` | `engine.rs:483` | `pub fn start_solve(&mut self, targets: Option<Vec<String>>) -> Result<JobId, Error>` —  Starts solving `targets` (default: the active physics, else all) on a job. The mesh is built inside the job when none is current. |
| fn | `Engine::start_solve_with` | `engine.rs:488` | `pub fn start_solve_with(&mut self, targets: Option<Vec<String>>, opts: &QueueOpts) -> Result<JobId, Error>` —  `start_solve` with queue options. The job waits in the queue while the heavy lane is busy; it carries a frozen copy of the project. |
| fn | `Engine::submit_solve` | `engine.rs:494` | `pub(crate) fn submit_solve( &mut self, src: Option<&Engine>, targets: Option<Vec<String>>, opts: &QueueOpts, purpose: Purpose, ) -> Result<JobId, Error>` —  The one place a solve is submitted. `src` is the engine whose project is solved (a point of a sweep is a scratch engine with the point's |
| fn | `Engine::queue_detached` | `engine.rs:865` | `pub fn queue_detached(&mut self, targets: Option<Vec<String>>, opts: &QueueOpts) -> Result<Value, Error>` —  Writes a solve into the queue of the workspace and leaves it there: nothing of this engine waits for it, so closing the front-end (or the |
| fn | `Engine::queue_detached` | `engine.rs:921` | `pub fn queue_detached(&mut self, _targets: Option<Vec<String>>, _opts: &QueueOpts) -> Result<Value, Error>` |
| fn | `Engine::start_convergence` | `engine.rs:929` | `pub fn start_convergence(&mut self, physics: &str, factor: f64, tolerance: Option<f64>) -> Result<JobId, Error>` —  Starts a mesh convergence check (two solves, the second on a finer mesh) on a job. Both solves run on the job's thread under its progress, |
| fn | `Engine::start_convergence_with` | `engine.rs:935` | `pub fn start_convergence_with( &mut self, physics: &str, factor: f64, tolerance: Option<f64>, opts: &QueueOpts, ) -> Result<JobId, Error>` —  `start_convergence` with queue options. When the solve of `physics` alone on the present mesh was done already by an earlier request, the |
| fn | `Engine::start_convergence_matched` | `engine.rs:947` | `pub fn start_convergence_matched( &mut self, physics: &str, factor: f64, tolerance: Option<f64>, matching: &crate::convergence::Matching, opts: &QueueOpts, ) -> Result<JobId, Error>` —  `start_convergence_with` with the modes of the two meshes paired as `matching` says (ADR 80); anything but the default pairing solves the |
| fn | `Engine::start_convergence_limited` | `engine.rs:967` | `pub fn start_convergence_limited( &mut self, physics: &str, factor: f64, tolerance: Option<f64>, matching: &crate::convergence::Matching, opts: &QueueOpts, limits: crate::convergence::Limits, ) -> Result<JobId, Error>` —  `start_convergence_matched` that refuses, before any solve, a check whose finer mesh has more cells than `limits.max_cells` or is estimated |
| fn | `Engine::start_open_results` | `engine.rs:1005` | `pub fn start_open_results(&mut self, path: &str) -> Result<JobId, Error>` —  Starts reading a results file, meshing its project again and checking the mesh against the stored one, on a job. The engine is not touched |
| fn | `Engine::start_export_results` | `engine.rs:1020` | `pub fn start_export_results(&mut self, path: &str, parallel: bool) -> Result<JobId, Error>` —  Starts writing the results the engine holds now (the project, the mesh and the latest result of every physics) into an HDF5 file on a job. |
| fn | `Engine::poll` | `engine.rs:1039` | `pub fn poll(&mut self) -> Vec<(JobId, Result<Value, Error>)>` —  Takes the finished jobs and applies them: the mesh and results become current, a solve is recorded in the run table when the project |
| fn | `Engine::job_blocked` | `engine.rs:1052` | `pub fn job_blocked(&self, id: JobId) -> Option<&'static str>` —  Why a queued job cannot start however long it is waited for (it is held, or the queue is paused); `None` for a job that is running or |
| fn | `Engine::wait` | `engine.rs:1073` | `pub fn wait(&mut self, id: JobId) -> Result<Value, Error>` —  Waits for a job and applies it. |
| fn | `Engine::apply` | `engine.rs:1088` | `pub(crate) fn apply(&mut self, id: JobId, outcome: Result<JobOutput, Error>) -> Result<Value, Error>` |
| fn | `Engine::learn` | `engine.rs:1296` | `pub(crate) fn learn( &self, project: &Project, results: &[SolveResult], cells: u64, dim: u8, worker: &Option<crate::jobs::WorkerRun>, )` —  Learns what a finished worker run cost: its peak memory against the size of the problem, for the estimates of later solves. Only runs |
| fn | `Engine::mesh_now` | `engine.rs:1379` | `pub fn mesh_now(&mut self) -> Result<Value, Error>` —  Mesh now (start and wait). |
| fn | `Engine::solve_now` | `engine.rs:1385` | `pub fn solve_now(&mut self, targets: Option<Vec<String>>) -> Result<Value, Error>` —  Solve now (start and wait). |
| fn | `Engine::save_to_workspace` | `engine.rs:1393` | `pub fn save_to_workspace(&mut self, workspace: Option<PathBuf>) -> Result<Value, Error>` —  Binds the project to a new workflow folder in `workspace` (or commits it, when already bound). |
| fn | `Engine::open_workflow` | `engine.rs:1412` | `pub fn open_workflow(&mut self, folder: &std::path::Path, experiment: Option<&str>) -> Result<Value, Error>` —  Opens a workflow folder and makes its experiment the project. |
| const | `TOPICS` | `guide.rs:6` | `&[&str]` |
| fn | `card` | `guide.rs:8` | `pub fn card(topic: &str) -> Option<&'static str>` |
| const | `SERVER_INSTRUCTIONS` | `guide.rs:20` | `&str` —  The text a server sends with `initialize`: just enough to find the guide (the `api` tool) and to know that the work is scripts. |
| const | `PATH_ARGS` | `path_guard.rs:17` | `&[&str]` —  The names of the arguments that hold a file or a folder, in every command of the catalog. `tests/path_guard.rs` fails when a command takes |
| fn | `absolute_root` | `path_guard.rs:20` | `pub(crate) fn absolute_root(root: &Path) -> Result<PathBuf, String>` —  A root made absolute, as the guard keeps it. |
| fn | `confine` | `path_guard.rs:29` | `pub fn confine(root: &Path, requested: &str) -> Result<PathBuf, String>` —  `requested` as an absolute path inside `root`, or the reason it is not allowed. |
| fn | `confine_args` | `path_guard.rs:72` | `pub(crate) fn confine_args(root: &Path, command: &str, args: &mut Map<String, Value>) -> Result<(), Error>` —  Rewrites the path arguments of one command call, or refuses it. |
| const | `BLOCK` | `reduce.rs:52` | `usize` —  Consecutive cells summed together before the blocks are added in order. |
| enum | `What` | `reduce.rs:60` | variants: `Field`, `Coord`, `Material`, `Region`, `Normal` —  What a leaf reads (N2): a field, a coordinate, a material value or the indicator of a region. All of them live on the mesh of one result. |
| struct | `Spec` | `reduce.rs:74` | fields: `physics`, `solution`, `name`, `what` —  One leaf of an expression: a named field of one solution of one physics (or, with `what`, a coordinate, a material value, a region). |
| struct | `Defaults` | `reduce.rs:83` |  Where the arguments say a bare name comes from. |
| struct | `Prep` | `reduce.rs:131` |  The leaves of one request resolved against the stored results, and the mesh they live on. |
| fn | `Engine::prepare` | `reduce.rs:860` | `pub(crate) fn prepare(&self, specs: &[Spec]) -> Result<Prep, Error>` —  Resolves the leaves against the stored results. |
| fn | `Engine::expressions` | `reduce.rs:971` | `pub(crate) fn expressions(&self, args: &Value, exprs: &[&Value]) -> Result<(Vec<Tree>, Prep), Error>` —  Parses expressions that share one leaf table (each a handle, a name or a number), resolved against the stored results. |
| fn | `Engine::expressions_for_solutions` | `reduce.rs:977` | `pub(crate) fn expressions_for_solutions( &self, d: &Defaults, solutions: &[usize], batch_solutions: bool, exprs: &[&Value], ) -> Result<(Vec<Tree>, Prep), Error>` —  Parses `exprs` across `solutions` (`sol_idx` outer, `expr_idx` inner) against one deduplicated `Prep`. |
| fn | `Engine::result_regions` | `reduce.rs:1105` | `pub fn result_regions(&self, physics: &str) -> Vec<String>` —  The names of the regions of the mesh a result was solved on (what the Results tab offers to reduce over). |
| fn | `Engine::cmd_field` | `reduce.rs:1111` | `pub(crate) fn cmd_field(&self, args: &Value) -> R` —  `field`: a handle to a field of a solution, with what it holds. A vector name (`E` when `Ex`, `Ez` exist) answers one handle per component. |
| fn | `Engine::cmd_field_op` | `reduce.rs:1172` | `pub(crate) fn cmd_field_op(&self, args: &Value) -> R` —  `field_op`: arithmetic on handles. The operation is checked against the results (every leaf must resolve) and the handle comes back. |
| fn | `Engine::cmd_integrate` | `reduce.rs:2871` | `pub(crate) fn cmd_integrate(&self, args: &Value) -> R` |
| fn | `Engine::cmd_extremum` | `reduce.rs:3133` | `pub(crate) fn cmd_extremum(&self, args: &Value) -> R` |
| fn | `Engine::cmd_probe` | `reduce.rs:3226` | `pub(crate) fn cmd_probe(&self, args: &Value) -> R` |
| fn | `Engine::cmd_line_cut` | `reduce.rs:3244` | `pub(crate) fn cmd_line_cut(&self, args: &Value) -> R` |
| fn | `Engine::cmd_overlap` | `reduce.rs:3281` | `pub(crate) fn cmd_overlap(&self, args: &Value) -> R` |
| fn | `Engine::cmd_confinement` | `reduce.rs:3315` | `pub(crate) fn cmd_confinement(&self, args: &Value) -> R` |
| struct | `Exported` | `scriptout.rs:27` | fields: `script`, `form`, `commands` —  The script and how it was made. |
| fn | `export` | `scriptout.rs:35` | `pub fn export(p: &Project) -> Result<Exported, Error>` —  The project as a script that rebuilds it (it does not solve). |
| const | `PLANNED_FACTOR` | `validated.rs:23` | `f64` —  The element-size ratio a check is planned with when none was made (the step the command's own help recommends). |
| struct | `Check` | `validated.rs:31` | fields: `factor`, `tolerance`, `matching`, `passed`, `worst_relative`, `worst_quantity`, `quantities`, `cells_coarse`, `cells_refined`, `seconds`, `grids`, `cells_medium`, ... —  What a mesh-convergence check with a tolerance found, for one problem key. |
| enum | `Validation` | `validated.rs:157` | variants: `Exploratory`, `Checked`, `NotAffordable` —  The three states, in the words every front-end shows. |
| enum | `ValidationTier` | `validated.rs:179` | variants: `BoundCertified`, `GciVerified`, `TwoGridChecked`, `ResidualChecked`, `Exploratory`, `NotAffordable`, `NonConvergent` —  Granular validation tier ordered strictly from strongest certainty to non-convergent. |
| struct | `Validated` | `validated.rs:218` | fields: `state`, `tier`, `reason`, `check` —  The state of one result with its validation tier, reason, and the check it rests on or the last one that did not pass. |
| fn | `Check::observable` | `validated.rs:60` | `pub fn observable(&self) -> String` —  What was compared and how, in words (the `observable` of the run record). |
| fn | `Check::from_answer` | `validated.rs:75` | `pub(crate) fn from_answer(v: &Value) -> Option<Check>` —  The part of a convergence answer that is a verdict; `None` for a check that had no tolerance (it says nothing is passed or failed). |
| fn | `Check::to_json` | `validated.rs:118` | `pub fn to_json(&self) -> Value` |
| fn | `Check::to_record` | `validated.rs:144` | `pub fn to_record(&self) -> nyx_store::runs::Convergence` |
| fn | `Validation::as_str` | `validated.rs:167` | `pub fn as_str(self) -> &'static str` |
| fn | `ValidationTier::as_str` | `validated.rs:190` | `pub fn as_str(self) -> &'static str` |
| fn | `ValidationTier::weakness` | `validated.rs:204` | `pub fn weakness(self) -> u8` —  Ordering of certainty: `"bound_certified"` / `"gci_verified"` &lt; `"two_grid_checked"` &lt; |
| fn | `Validated::to_json` | `validated.rs:226` | `pub fn to_json(&self) -> Value` |
| fn | `Engine::validated` | `validated.rs:280` | `pub fn validated(&self, physics: &str) -> Validated` —  The validated state of the latest result of `physics`. A physics without a result is `exploratory` with that as the reason. |
| fn | `Engine::validated_json` | `validated.rs:434` | `pub fn validated_json(&self, physics: &str) -> Value` —  `validated` as the JSON every result answer carries. |
| fn | `Engine::validated_of` | `validated.rs:439` | `pub(crate) fn validated_of(&self, physics: &[String]) -> Value` —  The weakest state of several results (a reduction over more than one physics), with one row per physics. |
| fn | `Engine::note_check` | `validated.rs:462` | `pub(crate) fn note_check(&mut self, physics: &str, key: &str, answer: &Value) -> Option<(Check, Option<String>)>` —  Keeps what a finished convergence check found, under the problem key it was started for, and writes it into the run that solved that |
| fn | `Engine::checked_answer` | `validated.rs:477` | `pub(crate) fn checked_answer(&mut self, physics: &str, key: &str, mut v: Value) -> Result<Value, Error>` —  The answer of a finished check with the state it leaves the result in. |

Trait implementations in this component: `Engine: Default`, `Namer: Names`
