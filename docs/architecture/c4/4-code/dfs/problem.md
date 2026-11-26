# Typed problem: code

C4 level 4 (code) of the component **Typed problem** (`dfs.problem`) in the container **nyx-dfs** (crate `nyx_dfs`). [Component view of the container](../../3-components/dfs.md) · [All code pages](../README.md)

> The typed dfs problem the kernels read, decoded from the ProblemSpec envelope.

**Decisions:** [ADR 41](../../../adr/0041-physics-crate-pattern-and-the-legacy-bridge.md), [ADR 103](../../../adr/0103-multi-rhs-capacitance-and-gauge-stabilized-magnetostatics.md)  
**Size:** 1 file, 897 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-dfs/src/problem.rs` | 897 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | material, material::bcs, props, validate | FilmModel, LinearSolver, Materials, PmlDir, bcs_surface_impedance, non_finite_token |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | cqed, decode, error, physics, spec | Ctx, Error, Params, Physics, ProblemSpec, TlsInterfaces, domain_names_by_label, materials_by_label, resolve_conditions, tls_interfaces_from_spec |
| [Physics adapter](./adapter.md) | crate root | Dfs |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | problem | Bc, BcKind, IeCoordKind, IeProfile, InfiniteElementSpec, Problem, Study, decode |
| [Shared helpers](./common.md) | problem | BcKind, IeCoordKind, IeProfile, InfiniteElementSpec, Problem |
| [Electrostatics kernel](./electro.md) | problem | BcKind, Problem |
| [Heat kernel](./heat.md) | problem | BcKind, Problem, Study |
| [Magnetostatics kernel](./magneto.md) | problem | BcKind, Problem, Study |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class IeCoordKind {
      <<enumeration>>
      Auto
      Cartesian
      Cylindrical
      Spherical
      +parse(s) Result
  }
  class IeProfile {
      <<enumeration>>
      Rational
      Exponential
      +parse(s) Result
  }
  class InfiniteElementSpec {
      <<struct>>
      +coord IeCoordKind
      +profile IeProfile
      +x PmlDir
      +y PmlDir
      +z PmlDir
      +center Option
      +axis slice
      +r_in Option~f64~
      +more_5_fields
  }
  class BcKind {
      <<enumeration>>
      Potential
      Ground
      SurfaceCharge
      Terminal
      Open
      FarField
      InfiniteElement
      Temperature
      Flux
      Convection
      LondonSheet
      CurrentTerminal
      more_1_variants
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
      +tol f64
      +amg_line_smooth bool
      +more_15_fields
      +new(study, order) Problem
      +with_tls_interfaces(tls_interfaces) Self
      +with_flux_bounds(flux_bounds) Self
      +with_freq(freq) Self
      +with_farfield_t_inf(t_inf) Self
      +has_infinite_element() bool
      +with_adapt_passes(adapt_passes) Self
      +with_adapt_estimator(adapt_estimator) Self
      +more_6_methods()
  }
  class fns_problem_rs {
      <<module problem.rs>>
      +decode(spec, ctx) Result
  }
  InfiniteElementSpec --> IeCoordKind
  InfiniteElementSpec --> IeProfile
  BcKind --> InfiniteElementSpec
  Bc --> BcKind
  Problem --> Study
  Problem --> Bc
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `IeCoordKind` | `problem.rs:19` | variants: `Auto`, `Cartesian`, `Cylindrical`, `Spherical` —  Coordinate topology of a domain-applied static Infinite Element (`BcKind::InfiniteElement`). |
| enum | `IeProfile` | `problem.rs:43` | variants: `Rational`, `Exponential` —  Radial/normal coordinate-stretching profile inside a static Infinite Element domain. |
| struct | `InfiniteElementSpec` | `problem.rs:63` | fields: `coord`, `profile`, `x`, `y`, `z`, `center`, `axis`, `r_in`, `thickness`, `pole_dist`, `decay`, `bounds_lo`, ... —  Configuration of a domain-applied static Infinite Element (`BcKind::InfiniteElement`). |
| enum | `BcKind` | `problem.rs:109` | variants: `Potential`, `Ground`, `SurfaceCharge`, `Terminal`, `Open`, `FarField`, `InfiniteElement`, `Temperature`, `Flux`, `Convection`, `LondonSheet`, `CurrentTerminal`, ... |
| struct | `Bc` | `problem.rs:210` | fields: `labels`, `kind` —  A condition on boundary numbers of the mesh. |
| enum | `Study` | `problem.rs:216` | variants: `Electrostatic`, `HeatStationary`, `HeatTransient`, `Magnetostatic` |
| struct | `Problem` | `problem.rs:244` | fields: `order`, `study`, `materials`, `domain_names`, `bcs`, `solver`, `tol`, `amg_line_smooth`, `formulation`, `singular_edges`, `richardson`, `tls_interfaces`, ... |
| fn | `decode` | `problem.rs:547` | `pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error>` —  Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and |
| fn | `IeCoordKind::parse` | `problem.rs:28` | `pub fn parse(s: &str) -> Result<Self, String>` |
| fn | `IeProfile::parse` | `problem.rs:52` | `pub fn parse(s: &str) -> Result<Self, String>` |
| fn | `BcKind::name` | `problem.rs:153` | `pub fn name(&self) -> &'static str` |
| fn | `BcKind::is_electrostatic` | `problem.rs:171` | `pub fn is_electrostatic(&self) -> bool` |
| fn | `BcKind::is_heat` | `problem.rs:184` | `pub fn is_heat(&self) -> bool` |
| fn | `BcKind::is_magnetostatic` | `problem.rs:195` | `pub fn is_magnetostatic(&self) -> bool` |
| fn | `Study::kind` | `problem.rs:233` | `pub fn kind(&self) -> &'static str` —  The `study.kind` of the envelope. |
| fn | `Problem::new` | `problem.rs:289` | `pub fn new(study: Study, order: usize) -> Problem` |
| fn | `Problem::with_tls_interfaces` | `problem.rs:317` | `pub fn with_tls_interfaces(mut self, tls_interfaces: Option<TlsInterfaces>) -> Self` |
| fn | `Problem::with_flux_bounds` | `problem.rs:322` | `pub fn with_flux_bounds(mut self, flux_bounds: impl Into<String>) -> Self` |
| fn | `Problem::with_freq` | `problem.rs:327` | `pub fn with_freq(mut self, freq: f64) -> Self` |
| fn | `Problem::with_farfield_t_inf` | `problem.rs:332` | `pub fn with_farfield_t_inf(mut self, t_inf: f64) -> Self` |
| fn | `Problem::has_infinite_element` | `problem.rs:338` | `pub fn has_infinite_element(&self) -> bool` —  True if any condition is an `InfiniteElement` or `materials.pml` has an active domain. |
| fn | `Problem::with_adapt_passes` | `problem.rs:342` | `pub fn with_adapt_passes(mut self, adapt_passes: usize) -> Self` |
| fn | `Problem::with_adapt_estimator` | `problem.rs:347` | `pub fn with_adapt_estimator(mut self, adapt_estimator: impl Into<String>) -> Self` |
| fn | `Problem::with_formulation` | `problem.rs:351` | `pub fn with_formulation(mut self, formulation: impl Into<String>) -> Self` |
| fn | `Problem::with_singular_edges` | `problem.rs:356` | `pub fn with_singular_edges(mut self, singular_edges: impl Into<String>) -> Self` |
| fn | `Problem::with_richardson` | `problem.rs:361` | `pub fn with_richardson(mut self, richardson: bool) -> Self` |
| fn | `Problem::with_amg_line_smooth` | `problem.rs:366` | `pub fn with_amg_line_smooth(mut self, amg_line_smooth: bool) -> Self` |
| fn | `Problem::with_tol` | `problem.rs:371` | `pub fn with_tol(mut self, tol: f64) -> Self` |
| fn | `Problem::check_finite` | `problem.rs:377` | `pub fn check_finite(&self) -> Result<(), String>` —  Every number the user can type must be finite: NaN or infinity does not stop a factorisation, it flows into every field. |

Trait implementations in this component: `InfiniteElementSpec: Default`
