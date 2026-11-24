# Analysis API: code

C4 level 4 (code) of the component **Analysis API** (`analysis.api`) in the container **nyx-analysis** (crate `nyx_analysis`). [Component view of the container](../../3-components/analysis.md) · [All code pages](../README.md)

> The entry points an engine calls: an analysis is a function of a SolveResult and its mesh.

**Decisions:** [ADR 56](../../../adr/0056-polarization-analysis-is-a-function-of-a-result.md), [ADR 103](../../../adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md), [ADR 105](../../../adr/0105-kicq-epr-craig-bampton-quantization-and-tls-participation.md)  
**Size:** 1 file, 127 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-analysis/src/lib.rs` | 127 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, error, result | Error, Mesh, SolveResult |
| [Capacitance and inductance matrices](./capacitance.md) | capacitance | CapacitanceMatrixAnalysis, InductanceMatrixAnalysis, capacitance_matrix, capacitance_matrix_applies, capacitance_matrix_with_mesh_and_params, capacitance_matrix_with_params, inductance_matrix, inductance_matrix_applies |
| [cQED quantization and TLS participation](./cqed.md) | tls | TlsParticipationAnalysis, TlsReport, tls_participation, tls_participation_applies, tls_participation_with_params |
| [Mode report](./modes.md) | modes | AnalysisOptions, GroupReport, ModeAnalysis, ModeInput, ModeRow, ResolveAxis, Resolved, analyze_modes |
| [Field samples](./samples.md) | samples | PolarizationInput, polarization_input |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Capacitance and inductance matrices](./capacitance.md) | crate root | Analysis |
| [cQED quantization and TLS participation](./cqed.md) | crate root | Analysis |
| [engine: Analyses](../engine/analyse.md) | crate root | Analysis, PolarizationAnalysis, polarization_applies |
| [engine: Registry and facade](../engine/registry.md) | crate root |  |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class Analysis {
      <<interface>>
      +id() str
      +name() str
      +applies_to(result) bool
      +run(mesh, result) Result
      +run_with_params(mesh, result, _params) Result
  }
  class PolarizationAnalysis {
      <<struct>>
  }
  class fns_lib_rs {
      <<module lib.rs>>
      +polarization(mesh, result) Result
      +polarization_applies(result) bool
  }
  PolarizationAnalysis ..|> Analysis
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| trait | `Analysis` | `lib.rs:55` | methods: `id`, `name`, `applies_to`, `run`, `run_with_params` —  One analysis: a name, a test of whether it applies to a result, and the run, which returns plain JSON data. |
| fn | `polarization` | `lib.rs:69` | `pub fn polarization(mesh: &Mesh, result: &SolveResult) -> Result<ModeAnalysis, Error>` —  The polarization of the modes of `result`, solved on `mesh`. |
| fn | `polarization_applies` | `lib.rs:75` | `pub fn polarization_applies(result: &SolveResult) -> bool` —  Whether `polarization` can analyse `result`. |
| struct | `PolarizationAnalysis` | `lib.rs:80` |  The polarization analysis behind the `Analysis` trait. |

Trait implementations in this component: `PolarizationAnalysis: Analysis`
