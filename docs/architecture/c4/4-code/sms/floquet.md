# Floquet kernel: code

C4 level 4 (code) of the component **Floquet kernel** (`sms.floquet`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals).

**Decisions:** [ADR 42](../../../adr/0042-solid-mechanics-moves-into-nyx-sms.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md)  
**Size:** 1 file, 389 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/floquet.rs` | 389 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem_mixed, linalg::csr, linalg::csr::assembly, linalg::driver, mesh::mixed, periodic, props | C, CCsr, FacetType, LagrangeSpace, MixedDiscretisation, MixedLagrangeDofMap, Periodic, QUAD_EDGES, RCsr, Reduction +7 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult, cell_fields_of |
| [Physics adapter](./adapter.md) | diag | SmsDiag, rigal_gaches_complex_residuals |
| [Typed problem](./problem.md) | problem |  |
| [Solid kernel](./solid.md) | solid | Setup, complex_mode_fields, complex_mode_fields_mixed |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Cross-section study kernel](./crosssection.md) | floquet | mixed_nodal_relations |
| [Piezoelectric kernel](./piezo.md) | floquet | mixed_nodal_relations |
| [Absorbing-layer kernel](./pml.md) | floquet | mixed_nodal_relations |
| [Solid kernel](./solid.md) | floquet | solve_floquet, solve_floquet_mixed |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class FloquetAssembly {
      <<struct>>
      +order usize
      +n usize
      +nred usize
      +n_pairs usize
      +red Reduction
      +kc CCsr
      +mc CCsr
      +n_modes usize
      +more_1_fields
  }
  class FloquetSolve {
      <<struct>>
      +modes Vec
      +vectors Vec
      +residuals Vec~f64~
      +restarts usize
      +dim usize
      +converged bool
  }
  class fns_floquet_rs {
      <<module floquet.rs>>
      ~mixed_nodal_relations(mixed, dof_map, ncomp, p) Result
      ~solve_floquet(mesh, spec, st, k, m, n_modes, f_search, t0) Result
      ~solve_floquet_mixed(raw_mesh, mixed, dof_map, spec, fixed, k, m, n_modes, f_search, t0) Result
      +assemble_floquet(mesh, spec, st, k, m, n_modes, f_search) Result
      +solve_floquet_modes(asm) Result
      +floquet_fields(mesh, asm, solved) Result
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `mixed_nodal_relations` | `floquet.rs:26` | `pub(crate) fn mixed_nodal_relations( mixed: &MixedDiscretisation, dof_map: &MixedLagrangeDofMap, ncomp: usize, p: &Periodic, ) -> Result<Vec<(usize, usize, C)>, String>` —  Builds Floquet-periodic DOF relations `(dest_dof, src_dof, exp(-i phase))` directly on a |
| fn | `solve_floquet` | `floquet.rs:171` | `pub(crate) fn solve_floquet( mesh: &Discretisation, spec: &Problem, st: &Setup, k: &Csr, m: &Csr, n_modes: usize, f_search: f64, t0: Instant, ) -> Result<SolveResult, String>` |
| fn | `solve_floquet_mixed` | `floquet.rs:188` | `pub(crate) fn solve_floquet_mixed( raw_mesh: &Mesh, mixed: &MixedDiscretisation, dof_map: &MixedLagrangeDofMap, spec: &Problem, fixed: &[usize], k: &Csr, m: &Csr, n_modes: usize, f_search: f64, t0: Instant, ) -> Result<SolveResult, String>` |
| struct | `FloquetAssembly` | `floquet.rs:263` | fields: `order`, `n`, `nred`, `n_pairs`, `red`, `kc`, `mc`, `n_modes`, `f_search` —  Stage 1 of the Floquet-periodic solid kernel: the stiffness and mass reduced to the independent dofs of the |
| struct | `FloquetSolve` | `floquet.rs:276` | fields: `modes`, `vectors`, `residuals`, `restarts`, `dim`, `converged` —  Stage 2 output: frequency and index of each converged mode (by frequency), the reduced vectors. |
| fn | `assemble_floquet` | `floquet.rs:285` | `pub fn assemble_floquet( mesh: &Discretisation, spec: &Problem, st: &Setup, k: &Csr, m: &Csr, n_modes: usize, f_search: f64, ) -> Result<FloquetAssembly, String>` |
| fn | `solve_floquet_modes` | `floquet.rs:315` | `pub fn solve_floquet_modes(asm: &FloquetAssembly) -> Result<FloquetSolve, String>` —  Stage 2: the complex Hermitian eigenpairs of the reduced problem nearest the search frequency. |
| fn | `floquet_fields` | `floquet.rs:342` | `pub fn floquet_fields( mesh: &Discretisation, asm: &FloquetAssembly, solved: &FloquetSolve, ) -> Result<SolveResult, String>` —  Stage 3: plot fields and per-cell complex displacement of every mode. |
