# Physics adapter: code

C4 level 4 (code) of the component **Physics adapter** (`dfs.adapter`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> The Physics implementation of dfs: descriptors, validation, R32 curved guards, solve entry.

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md), [ADR 25](../../../adr/0025-physics-trait-registry-and-problem-envelope.md), [ADR 103](../../../adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md), [ADR 133](../../../adr/0133-vvuq-error-budgets-gci-provenance-and-uq.md)  
**Size:** 1 file, 653 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/lib.rs` | 653 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material::spec, progress | Library, Progress |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root |  |
| [Electrostatics kernel](./electro.md) | electro, planar_bem | EPS0, ElectroAssembly, LayeredStack, assemble_electrostatics, electrostatics_fields, layered_greens_function, solve_electrostatics, solve_electrostatics_mixed, solve_planar_bem, solve_potential |
| [Heat kernel](./heat.md) | heat | HeatStates, HeatSystem, assemble_heat, heat_fields, solve_heat, solve_heat_mixed, solve_heat_states |
| [Magnetostatics kernel](./magneto.md) | magneto | MU0, solve_magnetostatics, solve_magnetostatics_mixed |
| [Typed problem](./problem.md) | problem | Bc, BcKind, IeCoordKind, IeProfile, InfiniteElementSpec, Problem, Study, decode |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Typed problem](./problem.md) | crate root | Dfs |
| [engine: Registry and facade](../engine/registry.md) | crate root | Dfs |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Dfs {
      <<struct>>
  }
  class fns_lib_rs {
      <<module lib.rs>>
      +run(mesh, problem, progress) Result
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `Dfs` | `lib.rs:34` |  |
| fn | `run` | `lib.rs:536` | `pub fn run(mesh: &Mesh, problem: &Problem, progress: &Arc<Progress>) -> Result<SolveResult, Error>` —  Solves a typed problem: checks the numbers, runs the kernel of the study under `progress` (stage log, cancellation) and appends the |

Trait implementations in this component: `Dfs: Physics`
