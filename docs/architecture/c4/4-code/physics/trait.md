# Physics trait and registry: code

C4 level 4 (code) of the component **Physics trait and registry** (`physics.trait`) in the container **nyx-physics** (crate `nyx_physics`). [Component view of the container](../../3-components/physics.md) · [All code pages](../README.md)

> The Physics trait (everything the rest of nyx may ask of a physics), the Registry of the physics a build contains, and the plain-data error every method returns.

**Decisions:** [ADR 25](../../../adr/0025-physics-trait-registry-and-problem-envelope.md)  
**Size:** 4 files, 343 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-physics/src/error.rs` | 66 |  |
| `crates/nyx-physics/src/lib.rs` | 81 |  |
| `crates/nyx-physics/src/physics.rs` | 103 |  |
| `crates/nyx-physics/src/registry.rs` | 93 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material::spec, mesh, mesh::exchange, mesh::mixed, progress | FacetType, Library, Mesh, MixedDiscretisation, Progress, SimplexMesh |
| [cQED result payloads](./cqed.md) | cqed | AdaptPassRecord, CAPACITANCE_KEY, CIRCUIT_ROM_KEY, CapacitanceMatrixData, CircuitRom, DRIVEN_PORT_KEY, DrivenPortData, ELECTRO_AMR_KEY, ElectroAmrData, INDUCTANCE_KEY +7 |
| [Descriptors and schemas](./descriptors.md) | crate root, desc | BcDesc, FieldDesc, FieldLink, ParamDesc, PropDesc, Sizing, StudyDesc |
| [Problem envelope](./envelope.md) | decode, spec, validate | BcSpec, FORMAT, MaterialAssign, ProblemSpec, PropValue, StudySpec, canonical_json, check_params, check_result_finite, tls_interface_from_material_name +2 |
| [Result and frame](./result.md) | frame, result | CellField, DiagnosticSeverity, DiagnosticWarning, FRAME_KEY, Field, FieldFrame, MixedToSimplexProjector, NumericalErrorBudget, Scalar, Solution +5 |

## Used by

| Component | Modules | Items |
|---|---|---|
| [analysis: Analysis API](../analysis/api.md) | crate root, error | Error, Mesh |
| [analysis: Capacitance and inductance matrices](../analysis/capacitance.md) | crate root, error | Error, Mesh, discretise |
| [analysis: cQED quantization and TLS participation](../analysis/cqed.md) | crate root, error | Discretisation, Error, Mesh, discretise |
| [analysis: Mode figures](../analysis/figures.md) | crate root, error | Discretisation, Error, Mesh, discretise |
| [analysis: Piezo resonator analysis](../analysis/piezo.md) | crate root, error | Error, Mesh |
| [analysis: Field samples](../analysis/samples.md) | crate root, error | Error, Mesh, discretise |
| [dfs: Physics adapter](../dfs/adapter.md) | crate root |  |
| [dfs: Shared helpers](../dfs/common.md) | crate root | Discretisation, Mesh |
| [dfs: Electrostatics kernel](../dfs/electro.md) | crate root | Discretisation, Mesh, discretise, discretise_mixed |
| [dfs: Heat kernel](../dfs/heat.md) | crate root | Discretisation, Mesh, discretise, discretise_mixed |
| [dfs: Magnetostatics kernel](../dfs/magneto.md) | crate root | Discretisation, Mesh, discretise, discretise_mixed |
| [dfs: Typed problem](../dfs/problem.md) | error, physics | Ctx, Error, Physics |
| [ems: Physics adapter](../ems/adapter.md) | crate root | Mesh |
| [ems: Axisymmetric kernel](../ems/axisym.md) | crate root | Discretisation, Mesh |
| [ems: 3D cavity kernel](../ems/cavity.md) | crate root | Discretisation, Mesh |
| [ems: Driven frequency-domain and PROM](../ems/driven.md) | crate root | Discretisation, Mesh |
| [ems: Thin-film IBC and Josephson ROM](../ems/ibc.md) | crate root | Discretisation |
| [ems: Native mixed-cell EM](../ems/mixed.md) | crate root | Discretisation |
| [ems: Mode kernel](../ems/modes.md) | crate root | Discretisation, Mesh |
| [ems: Complex-mode kernel](../ems/modescx.md) | crate root | Discretisation, Mesh |
| [ems: Optical mode sets](../ems/optical.md) | crate root | Discretisation |
| [ems: Typed problem](../ems/problem.md) | crate root, error, physics | Ctx, Discretisation, Error, MixedDiscretisation, Physics |
| [engine: Analyses](../engine/analyse.md) | error, physics | Error, LabelResolver |
| [engine: Command lists](../engine/commandlist.md) | error | Error |
| [engine: Commands](../engine/commands.md) | crate root, error, physics, registry | Ctx, Error, LabelResolver, Registry, Sizing |
| [engine: Geometry operations](../engine/geometry.md) | error | Error |
| [engine: Job runner and queue](../engine/jobs.md) | error, physics, registry | Ctx, Error, Physics, Registry |
| [engine: Materials](../engine/materials.md) | error, physics | Error, Physics |
| [engine: Meshing seam](../engine/meshing.md) | crate root, error, physics, registry | Error, LabelResolver, Mesh, Registry, Sizing |
| [engine: Project](../engine/project.md) | error, registry | Error, Registry |
| [engine: Registry and facade](../engine/registry.md) | crate root, error, physics, registry | Error, Physics, Registry |
| [engine: Results files](../engine/resultsIO.md) | error, registry | Error, Registry |
| [engine: Store link](../engine/storeLink.md) | error | Error |
| [Descriptors and schemas](./descriptors.md) | physics | Physics |
| [Problem envelope](./envelope.md) | error, physics | Ctx, Error, LabelResolver, Physics |
| [Result and frame](./result.md) | crate root, error | Error, Mesh, discretise |
| [Test support](./testsupport.md) | crate root, physics | Discretisation, Physics |
| [sms: Physics adapter](../sms/adapter.md) | crate root | Mesh |
| [sms: Axisymmetric kernel](../sms/axisym.md) | crate root | Discretisation |
| [sms: Cross-section study kernel](../sms/crosssection.md) | crate root | Discretisation, Mesh |
| [sms: Floquet kernel](../sms/floquet.md) | crate root | Discretisation, Mesh |
| [sms: Low-reflecting boundary](../sms/lrb.md) | crate root | Discretisation, Mesh |
| [sms: Piezoelectric kernel](../sms/piezo.md) | crate root | Discretisation, Mesh |
| [sms: Absorbing-layer kernel](../sms/pml.md) | crate root | Discretisation, Mesh |
| [sms: Typed problem](../sms/problem.md) | error, physics | Ctx, Error, Physics |
| [sms: Solid kernel](../sms/solid.md) | crate root | Discretisation, Mesh, discretise_mixed |
| [sms: Waveguide kernel](../sms/waveguide.md) | crate root | Discretisation |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace error_rs {
    class Error {
        <<enumeration>>
        Invalid
        Unsupported
        Solver
        NonFinite
        Cancelled
        NotYetMoved
        +code() str
        +from_kernel(text) Error
        +invalid(text) Error
    }
  }
  namespace lib_rs {
    class fns_lib_rs {
        <<module lib.rs>>
        +discretise(mesh) Result
        +discretise_mixed(mesh) Result
        +check_mesh(mesh) Result
        +Discretisation type
    }
  }
  namespace physics_rs {
    class LabelResolver {
        <<interface>>
        +domains(region) Result
        +boundaries(selector) Result
        +sheet_thickness_m() BTreeMap
        +sheet_etch_depth_m() BTreeMap
    }
    class Ctx {
        <<struct>>
        +resolver LabelResolver
        +progress Arc~Progress~
        +library Option
        +wavelength Option~f64~
        +frequency Option~f64~
        +temperature Option~f64~
        +inputs BTreeMap
        +new(resolver) Ctx
        +check() Result
    }
    class Physics {
        <<interface>>
        +id() str
        +name() str
        +studies() Vec~StudyDesc~
        +bc_kinds() Vec~BcDesc~
        +properties() Vec~PropDesc~
        +params() Vec~ParamDesc~
        +fields(study) Vec~FieldDesc~
        +exports() Vec~FieldLink~
        +more_4_methods()
    }
  }
  namespace registry_rs {
    class Registry {
        <<struct>>
        +new() Registry
        +register(p) Result
        +get(id) Option
        +ids() Vec
        +len() usize
        +is_empty() bool
        +solve_order(ids) Result
    }
    class fns_registry_rs {
        <<module registry.rs>>
        +is_valid_code(id) bool
    }
  }
  Ctx --> LabelResolver
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `Error` | `error.rs:8` | variants: `Invalid`, `Unsupported`, `Solver`, `NonFinite`, `Cancelled`, `NotYetMoved` |
| fn | `Error::code` | `error.rs:25` | `pub fn code(&self) -> &'static str` —  A stable short name for the kind, for logs, JSON and MCP error codes. |
| fn | `Error::from_kernel` | `error.rs:37` | `pub fn from_kernel(text: String) -> Error` —  Maps the `String` errors of the kernels in `femcore` (the text `"cancelled"` is the cancel signal there). |
| fn | `Error::invalid` | `error.rs:47` | `pub fn invalid(text: impl Into<String>) -> Error` |
| type | `Discretisation` | `lib.rs:49` |  The numbering of unknowns (edges, P2 dofs) a kernel assembles on. Each solve builds its own from the `Mesh` with `discretise`. The type |
| fn | `discretise` | `lib.rs:52` | `pub fn discretise(mesh: &Mesh) -> Result<Discretisation, Error>` —  The `Discretisation` a physics assembles on, made from the `Mesh` it was given (non-simplex cells are cut into simplices). |
| fn | `discretise_mixed` | `lib.rs:57` | `pub fn discretise_mixed(mesh: &Mesh) -> Result<MixedDiscretisation, Error>` —  Native mixed-cell discretisation made directly from `mesh` without splitting `Quad`, `Prism`, or `Hex` cells. |
| fn | `check_mesh` | `lib.rs:63` | `pub fn check_mesh(mesh: &Mesh) -> Result<(), Error>` —  Refuses a mesh with a NaN or infinite coordinate, naming the vertex and the axis, or an invalid |
| trait | `LabelResolver` | `physics.rs:21` | methods: `domains`, `boundaries`, `sheet_thickness_m`, `sheet_etch_depth_m` —  Turns the names of a problem (`geom` region names and boundary selectors) into the integer labels of the mesh being solved. |
| struct | `Ctx` | `physics.rs:38` | fields: `resolver`, `progress`, `library`, `wavelength`, `frequency`, `temperature`, `inputs` —  What a solve may use besides the mesh and the problem: the label resolver, progress and cancellation, and the fields that other |
| trait | `Physics` | `physics.rs:76` | methods: `id`, `name`, `studies`, `bc_kinds`, `properties`, `params`, `fields`, `exports`, `imports`, `validate`, `sizing`, `solve` |
| fn | `Ctx::new` | `physics.rs:54` | `pub fn new(resolver: &'a dyn LabelResolver) -> Ctx<'a>` |
| fn | `Ctx::check` | `physics.rs:67` | `pub fn check(&self) -> Result<(), Error>` —  `Err(Cancelled)` when the run was cancelled; kernels call it between stages. |
| struct | `Registry` | `registry.rs:12` |  |
| fn | `is_valid_code` | `registry.rs:17` | `pub fn is_valid_code(id: &str) -> bool` —  A code is exactly three lowercase ASCII letters (ADR 12). |
| fn | `Registry::new` | `registry.rs:22` | `pub fn new() -> Registry` |
| fn | `Registry::register` | `registry.rs:27` | `pub fn register(&mut self, p: Arc<dyn Physics>) -> Result<(), Error>` —  Adds a physics. Refuses a code that is not three lowercase letters and a code that is already taken. |
| fn | `Registry::get` | `registry.rs:39` | `pub fn get(&self, id: &str) -> Option<&dyn Physics>` |
| fn | `Registry::ids` | `registry.rs:44` | `pub fn ids(&self) -> Vec<&'static str>` —  The codes in alphabetical order. |
| fn | `Registry::len` | `registry.rs:48` | `pub fn len(&self) -> usize` |
| fn | `Registry::is_empty` | `registry.rs:52` | `pub fn is_empty(&self) -> bool` |
| fn | `Registry::solve_order` | `registry.rs:58` | `pub fn solve_order(&self, ids: &[&str]) -> Result<Vec<String>, Error>` —  The order in which `ids` must be solved so that every import finds its export first. An import with no matching export among `ids` is |

Trait implementations in this component: `Error: Display`, `Error: Error`
