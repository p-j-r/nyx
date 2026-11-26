# Heat kernel: code

C4 level 4 (code) of the component **Heat kernel** (`dfs.heat`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> Heat transfer on P1/P2 simplices, stationary and transient (`SourceField` MMS, BDF2 time-varying forcing & temporal error estimate, R32 curved guard).

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md), [ADR 9](../../../adr/0009-kernels-are-three-stages.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md), [ADR 133](../../../adr/0133-vvuq-error-budgets-gci-provenance-and-uq.md)  
**Size:** 1 file, 1283 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/heat.rs` | 1283 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | expr, fem, fem_mixed, fields, linalg::csr::assembly, linalg::driver, material, mesh::mixed +3 | AmgSystem, CellFields, DirectSolveCheck, DirectSolver, LagrangeSpace, LinearSolver, M3, MixedDiscretisation, MixedLagrangeDofMap, RCsr +10 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult, cell_fields_of, discretise, discretise_mixed |
| [Shared helpers](./common.md) | common, diag | DfsDiag, amg_solve_free_checked, assemble_boundary_load_field, assemble_boundary_load_mixed, assemble_boundary_mass_mixed, assemble_farfield_robin_mixed, assemble_ie_outer_robin_mixed, assemble_laplace_ie_mixed, assemble_source_field, assemble_source_mixed +16 |
| [Typed problem](./problem.md) | problem | BcKind, Problem, Study |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | heat | HeatStates, HeatSystem, assemble_heat, heat_fields, solve_heat, solve_heat_mixed, solve_heat_states |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class TransientHeatOptions {
      <<struct>>
      +bdf2_substeps usize
      +force_bdf2 bool
  }
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
      +assemble_heat_with_source(mesh, spec, source) Result
      +heat_fields(mesh, t, ktensor, order) tuple
      +solve_heat_states(mesh, spec, sys) Result
      ~solve_heat_states_checked(mesh, spec, sys) Result
      +solve_heat_states_with_source(mesh, spec, sys, source, opts) Result
      ~solve_heat_states_with_source_checked(mesh, spec, sys, source, opts) Result
      +solve_heat_with_source(mesh, spec, source) Result
      +solve_heat_bdf2_with_source(mesh, spec, source, opts) Result
      +solve_heat(mesh, spec) Result
      +more_1_functions()
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `TransientHeatOptions` | `heat.rs:63` | fields: `bdf2_substeps`, `force_bdf2` —  Options controlling transient BDF2 time integration in [`solve_heat_states_with_source`] and |
| struct | `HeatSystem` | `heat.rs:77` | fields: `k`, `f`, `ktensor`, `fixed`, `vals`, `has_convection` —  Stage 1 of heat transfer: conductivity matrix, load, Dirichlet dofs and the conductivity tensors (for the flux). |
| fn | `assemble_heat` | `heat.rs:86` | `pub fn assemble_heat(mesh: &Discretisation, spec: &Problem) -> Result<HeatSystem, String>` |
| fn | `assemble_heat_with_source` | `heat.rs:94` | `pub fn assemble_heat_with_source<const D: usize, S: SourceField<D>>( mesh: &Discretisation, spec: &Problem, source: &S, ) -> Result<HeatSystem, String>` |
| fn | `heat_fields` | `heat.rs:284` | `pub fn heat_fields(mesh: &Discretisation, t: &[f64], ktensor: &[M3], order: usize) -> (Vec<Field>, CellFields)` —  Stage 3 of heat transfer: `T` and `q` (flat names `\|q\|`, `qx`, `qz` come from `Solution::flat_fields`) and the same per cell without |
| struct | `HeatStates` | `heat.rs:298` | fields: `times`, `temps`, `log` —  Stage 2 of heat transfer: the temperature at every dof, for the stationary study (one state) or at each output |
| fn | `solve_heat_states` | `heat.rs:330` | `pub fn solve_heat_states(mesh: &Discretisation, spec: &Problem, sys: &HeatSystem) -> Result<HeatStates, String>` |
| fn | `solve_heat_states_checked` | `heat.rs:334` | `pub(crate) fn solve_heat_states_checked( mesh: &Discretisation, spec: &Problem, sys: &HeatSystem, ) -> Result<(HeatStates, DfsDiag), String>` |
| fn | `solve_heat_states_with_source` | `heat.rs:346` | `pub fn solve_heat_states_with_source<const D: usize, S: SourceField<D>>( mesh: &Discretisation, spec: &Problem, sys: &HeatSystem, source: &S, opts: Option<TransientHeatOptions>, ) -> Result<HeatStates, String>` |
| fn | `solve_heat_states_with_source_checked` | `heat.rs:356` | `pub(crate) fn solve_heat_states_with_source_checked<const D: usize, S: SourceField<D>>( mesh: &Discretisation, spec: &Problem, sys: &HeatSystem, source: &S, opts: Option<TransientHeatOptions>, ) -> Result<(HeatStates, DfsDiag), String>` |
| fn | `solve_heat_with_source` | `heat.rs:852` | `pub fn solve_heat_with_source<const D: usize, S: SourceField<D>>( mesh: &Discretisation, spec: &Problem, source: &S, ) -> Result<SolveResult, String>` |
| fn | `solve_heat_bdf2_with_source` | `heat.rs:860` | `pub fn solve_heat_bdf2_with_source<const D: usize, S: SourceField<D>>( mesh: &Discretisation, spec: &Problem, source: &S, opts: TransientHeatOptions, ) -> Result<SolveResult, String>` |
| fn | `solve_heat` | `heat.rs:906` | `pub fn solve_heat(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String>` |
| fn | `solve_heat_mixed` | `heat.rs:971` | `pub fn solve_heat_mixed(mesh: &Mesh, spec: &Problem) -> Result<SolveResult, String>` —  Solves stationary or transient heat transfer natively on a mixed mesh (`Quad`, `Prism`, `Hex`) |

Trait implementations in this component: `TransientHeatOptions: Default`
