# Typed problem: code

C4 level 4 (code) of the component **Typed problem** (`dfs.problem`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> The typed dfs problem the kernels read, decoded from the ProblemSpec envelope.

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md)  
**Size:** 1 file, 197 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/problem.rs` | 197 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | props | LinearSolver, Materials |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | decode, error, physics, spec | Ctx, Error, Params, Physics, ProblemSpec, materials_by_label, resolve_conditions |
| [Physics adapter](./adapter.md) | crate root | Dfs |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | problem | Bc, BcKind, Problem, Study, decode |
| [Electrostatics kernel](./electro.md) | problem | BcKind, Problem |
| [Heat kernel](./heat.md) | problem | BcKind, Problem, Study |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class BcKind {
      <<enumeration>>
      Potential
      Ground
      SurfaceCharge
      Temperature
      Flux
      Convection
      +name() str
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
      +kind() str
  }
  class Problem {
      <<struct>>
      +order usize
      +study Study
      +materials Materials
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
| enum | `BcKind` | `problem.rs:14` | variants: `Potential`, `Ground`, `SurfaceCharge`, `Temperature`, `Flux`, `Convection` |
| struct | `Bc` | `problem.rs:47` | fields: `labels`, `kind` —  A condition on boundary numbers of the mesh. |
| enum | `Study` | `problem.rs:53` | variants: `Electrostatic`, `HeatStationary`, `HeatTransient` |
| struct | `Problem` | `problem.rs:76` | fields: `order`, `study`, `materials`, `bcs`, `solver` |
| fn | `decode` | `problem.rs:140` | `pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error>` —  Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and |
| fn | `BcKind::name` | `problem.rs:29` | `pub fn name(&self) -> &'static str` |
| fn | `Study::kind` | `problem.rs:66` | `pub fn kind(&self) -> &'static str` —  The `study.kind` of the envelope. |
| fn | `Problem::new` | `problem.rs:86` | `pub fn new(study: Study, order: usize) -> Problem` |
| fn | `Problem::check_finite` | `problem.rs:91` | `pub fn check_finite(&self) -> Result<(), String>` —  Every number the user can type must be finite: NaN or infinity does not stop a factorisation, it flows into every field. |
