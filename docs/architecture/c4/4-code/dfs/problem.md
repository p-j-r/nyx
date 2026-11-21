# Typed problem: code

C4 level 4 (code) of the component **Typed problem** (`dfs.problem`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> The typed dfs problem the kernels read, decoded from the ProblemSpec envelope.

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md), [ADR 103](../../../adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md)  
**Size:** 1 file, 314 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/problem.rs` | 314 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material, material::bcs, props | FilmModel, LinearSolver, Materials, bcs_surface_impedance |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | decode, error, physics, spec | Ctx, Error, Params, Physics, ProblemSpec, domain_names_by_label, materials_by_label, resolve_conditions |
| [Physics adapter](./adapter.md) | crate root | Dfs |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | problem | Bc, BcKind, Problem, Study, decode |
| [Electrostatics kernel](./electro.md) | problem | BcKind, Problem |
| [Heat kernel](./heat.md) | problem | BcKind, Problem, Study |
| [Magnetostatics kernel](./magneto.md) | problem | BcKind, Problem, Study |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class BcKind {
      <<enumeration>>
      Potential
      Ground
      SurfaceCharge
      Terminal
      Temperature
      Flux
      Convection
      LondonSheet
      CurrentTerminal
      FluxLoop
      +name() str
      +is_electrostatic() bool
      +is_heat() bool
      +is_magnetostatic() bool
  }
  class Bc {
      <<struct>>
      +labels Vec~i32~
      +kind BcKind
  }
  class Study {
      <<enumeration>>
      Electrostatic
      HeatStationary
      HeatTransient
      Magnetostatic
      +kind() str
  }
  class Problem {
      <<struct>>
      +order usize
      +study Study
      +materials Materials
      +domain_names BTreeMap
      +bcs Vec~Bc~
      +solver LinearSolver
      +new(study, order) Problem
      +check_finite() Result
  }
  class fns_problem_rs {
      <<module problem.rs>>
      +decode(spec, ctx) Result
  }
  Bc --> BcKind
  Problem --> Study
  Problem --> Bc
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `BcKind` | `problem.rs:16` | variants: `Potential`, `Ground`, `SurfaceCharge`, `Terminal`, `Temperature`, `Flux`, `Convection`, `LondonSheet`, `CurrentTerminal`, `FluxLoop` |
| struct | `Bc` | `problem.rs:81` | fields: `labels`, `kind` —  A condition on boundary numbers of the mesh. |
| enum | `Study` | `problem.rs:87` | variants: `Electrostatic`, `HeatStationary`, `HeatTransient`, `Magnetostatic` |
| struct | `Problem` | `problem.rs:115` | fields: `order`, `study`, `materials`, `domain_names`, `bcs`, `solver` |
| fn | `decode` | `problem.rs:192` | `pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error>` —  Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and |
| fn | `BcKind::name` | `problem.rs:48` | `pub fn name(&self) -> &'static str` |
| fn | `BcKind::is_electrostatic` | `problem.rs:63` | `pub fn is_electrostatic(&self) -> bool` |
| fn | `BcKind::is_heat` | `problem.rs:67` | `pub fn is_heat(&self) -> bool` |
| fn | `BcKind::is_magnetostatic` | `problem.rs:71` | `pub fn is_magnetostatic(&self) -> bool` |
| fn | `Study::kind` | `problem.rs:104` | `pub fn kind(&self) -> &'static str` —  The `study.kind` of the envelope. |
| fn | `Problem::new` | `problem.rs:126` | `pub fn new(study: Study, order: usize) -> Problem` |
| fn | `Problem::check_finite` | `problem.rs:138` | `pub fn check_finite(&self) -> Result<(), String>` —  Every number the user can type must be finite: NaN or infinity does not stop a factorisation, it flows into every field. |
