# Physics adapter: code

C4 level 4 (code) of the component **Physics adapter** (`sms.adapter`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> The Physics implementation of sms with five studies (static, eigenfrequency, axisymmetric eigenfrequency, waveguide modes, driven frequency sweep): descriptors, validation, solve entry.

**Decisions:** [ADR 42](../../../adr/0042-solid-mechanics-moves-into-nyx-sms.md), [ADR 25](../../../adr/0025-physics-trait-registry-and-problem-envelope.md), [ADR 68](../../../adr/0068-the-axisymmetric-solid-gets-a-study.md), [ADR 101](../../../adr/0101-sms-waveguide-modes-and-floquet-pml.md), [ADR 129](../../../adr/0129-coupled-piezoelectric-fem-in-nyx-sms.md)  
**Size:** 1 file, 670 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/lib.rs` | 670 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material, material::coupled, material::spec, meshing::physics_size, progress | EPSILON_0, Library, Progress, elastic_wavelength, isotropic_stiffness, piezo_e_global, piezo_wave_speed_range, stiffness_global, tensor_global |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, decode | Mesh, Params, materials_by_region |
| [Cross-section study kernel](./crosssection.md) | cross_section | CrossModes, solve_axisymmetric, solve_axisymmetric_mixed, solve_solid_axisym, solve_solid_waveguide, solve_waveguide_modes, solve_waveguide_modes_mixed |
| [Typed problem](./problem.md) | problem | Bc, BcKind, ConductorRole, IeProfile, InfiniteElement, PiezoMode, Plane, Problem, Study, ThermalExpansion +2 |
| [Solid kernel](./solid.md) | solid | Setup, solve_solid, solve_solid_mixed |

## Used by

| Component | Modules | Items |
|---|---|---|
| [engine: Registry and facade](../engine/registry.md) | crate root | Sms |
| [Typed problem](./problem.md) | crate root | Sms |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Sms {
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
| struct | `Sms` | `lib.rs:38` |  |
| fn | `run` | `lib.rs:595` | `pub fn run(mesh: &Mesh, problem: &Problem, progress: &Arc<Progress>) -> Result<SolveResult, Error>` —  Solves a typed problem: checks the numbers, assembles natively on `MixedDiscretisation` whenever |

Trait implementations in this component: `Sms: Physics`
