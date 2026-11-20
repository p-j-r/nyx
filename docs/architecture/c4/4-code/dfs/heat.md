# Heat kernel: code

C4 level 4 (code) of the component **Heat kernel** (`dfs.heat`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> Heat transfer on P1/P2 simplices, stationary and transient.

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md), [ADR 9](../../../adr/0009-kernels-are-three-stages.md)  
**Size:** 1 file, 275 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/heat.rs` | 275 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | expr, fem, fields, linalg::csr::assembly, linalg::driver, material, props, section | AmgSystem, CellFields, DirectSolver, LinearSolver, M3, RCsr, TensorKind, amg_wanted, component_names, dense_gen_eigh +2 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Field, Scalar, Solution, SolveResult, cell_fields_of |
| [Shared helpers](./common.md) | common | cell_tensor_gradient, components, minmax, plane_tensors |
| [Typed problem](./problem.md) | problem | BcKind, Problem, Study |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | heat | HeatStates, HeatSystem, assemble_heat, heat_fields, solve_heat, solve_heat_states |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class HeatSystem {
      <<struct>>
      +k Csr
      +f Vec~f64~
      +ktensor Vec~M3~
      +fixed Vec~usize~
      +vals Vec~f64~
      +has_convection bool
  }
  class HeatStates {
      <<struct>>
      +times Vec~f64~
      +temps Vec
      +log Vec~String~
  }
  class fns_heat_rs {
      <<module heat.rs>>
      +assemble_heat(mesh, spec) Result
      +heat_fields(mesh, t, ktensor, order) tuple
      +solve_heat_states(mesh, spec, sys) Result
      +solve_heat(mesh, spec) Result
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `HeatSystem` | `heat.rs:20` | fields: `k`, `f`, `ktensor`, `fixed`, `vals`, `has_convection` —  Stage 1 of heat transfer: conductivity matrix, load, Dirichlet dofs and the conductivity tensors (for the flux). |
| fn | `assemble_heat` | `heat.rs:29` | `pub fn assemble_heat(mesh: &Discretisation, spec: &Problem) -> Result<HeatSystem, String>` |
| fn | `heat_fields` | `heat.rs:72` | `pub fn heat_fields(mesh: &Discretisation, t: &[f64], ktensor: &[M3], order: usize) -> (Vec<Field>, CellFields)` —  Stage 3 of heat transfer: `T` and `q` (flat names `\|q\|`, `qx`, `qz` come from `Solution::flat_fields`) and the same per cell without |
| struct | `HeatStates` | `heat.rs:86` | fields: `times`, `temps`, `log` —  Stage 2 of heat transfer: the temperature at every dof, for the stationary study (one state) or at each output |
| fn | `solve_heat_states` | `heat.rs:93` | `pub fn solve_heat_states(mesh: &Discretisation, spec: &Problem, sys: &HeatSystem) -> Result<HeatStates, String>` |
| fn | `solve_heat` | `heat.rs:239` | `pub fn solve_heat(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String>` |
