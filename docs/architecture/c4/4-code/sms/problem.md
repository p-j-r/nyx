# Typed problem: code

C4 level 4 (code) of the component **Typed problem** (`sms.problem`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> The typed sms problem the kernels read, decoded from the envelope.

**Decisions:** [ADR 42](../../../adr/0042-solid-mechanics-moves-into-nyx-sms.md), [ADR 68](../../../adr/0068-the-axisymmetric-solid-gets-a-study.md)  
**Size:** 1 file, 375 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/problem.rs` | 375 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | props, validate | LinearSolver, Materials, Periodic, Pml, PmlDir, non_finite_token |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | decode, error, physics, spec | Ctx, Error, Params, Physics, ProblemSpec, materials_by_label, resolve_conditions |
| [Physics adapter](./adapter.md) | crate root | Sms |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | problem | Bc, BcKind, Plane, Problem, Study, ThermalExpansion, ThermalSource, decode |
| [Axisymmetric study kernel](./crosssection.md) | problem |  |
| [Floquet kernel](./floquet.md) | problem |  |
| [Low-reflecting boundary](./lrb.md) | problem | Problem |
| [Absorbing-layer kernel](./pml.md) | problem |  |
| [Solid kernel](./solid.md) | problem |  |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Plane {
      <<enumeration>>
      Strain
      Stress
  }
  class BcKind {
      <<enumeration>>
      Fixed
      Load
      LowReflecting
      +name() str
  }
  class Bc {
      <<struct>>
      +labels Vec~i32~
      +kind BcKind
  }
  class ThermalSource {
      <<enumeration>>
      Uniform
      Imported
  }
  class ThermalExpansion {
      <<struct>>
      +source ThermalSource
      +t_ref f64
  }
  class Study {
      <<enumeration>>
      Static
      Eigenfrequency
      AxisymmetricEigenfrequency
      +kind() str
      +is_eigen() bool
  }
  class Problem {
      <<struct>>
      +order usize
      +study Study
      +materials Materials
      +bcs Vec~Bc~
      +plane Plane
      +thickness f64
      +thermal Option~ThermalExpansion~
      +solver LinearSolver
      +more_2_fields
      +new(study, order) Problem
      +check_finite() Result
      +check_wafer(dim) Result
      +check_axisymmetric() Result
  }
  class fns_problem_rs {
      <<module problem.rs>>
      +decode(spec, ctx) Result
  }
  Bc --> BcKind
  ThermalExpansion --> ThermalSource
  Problem --> Study
  Problem --> Bc
  Problem --> Plane
  Problem --> ThermalExpansion
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `Plane` | `problem.rs:14` | variants: `Strain`, `Stress` |
| enum | `BcKind` | `problem.rs:20` | variants: `Fixed`, `Load`, `LowReflecting` |
| struct | `Bc` | `problem.rs:41` | fields: `labels`, `kind` —  A condition on boundary numbers of the mesh. |
| enum | `ThermalSource` | `problem.rs:47` | variants: `Uniform`, `Imported` |
| struct | `ThermalExpansion` | `problem.rs:54` | fields: `source`, `t_ref` |
| enum | `Study` | `problem.rs:60` | variants: `Static`, `Eigenfrequency`, `AxisymmetricEigenfrequency` |
| struct | `Problem` | `problem.rs:92` | fields: `order`, `study`, `materials`, `bcs`, `plane`, `thickness`, `thermal`, `solver`, `inactive`, `periodic` |
| fn | `decode` | `problem.rs:245` | `pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error>` —  Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and |
| fn | `BcKind::name` | `problem.rs:30` | `pub fn name(&self) -> &'static str` |
| fn | `Study::kind` | `problem.rs:77` | `pub fn kind(&self) -> &'static str` —  The `study.kind` of the envelope. |
| fn | `Study::is_eigen` | `problem.rs:86` | `pub fn is_eigen(&self) -> bool` —  Whether the study finds eigenfrequencies (`Eigenfrequency` or `AxisymmetricEigenfrequency`). |
| fn | `Problem::new` | `problem.rs:111` | `pub fn new(study: Study, order: usize) -> Problem` |
| fn | `Problem::check_finite` | `problem.rs:128` | `pub fn check_finite(&self) -> Result<(), String>` —  Every number the user can type must be finite: NaN or infinity does not stop a factorisation or a Krylov iteration, it flows into |
| fn | `Problem::check_wafer` | `problem.rs:183` | `pub fn check_wafer(&self, dim: usize) -> Result<(), String>` —  Refuses what a wafer-frame 2D problem cannot say: there is no y axis in the plane of a cross-section. |
| fn | `Problem::check_axisymmetric` | `problem.rs:199` | `pub fn check_axisymmetric(&self) -> Result<(), Error>` —  The rules of the axisymmetric study that do not need the mesh (ADR 68): no boundary load, no Floquet pair, no low-reflecting boundary |
