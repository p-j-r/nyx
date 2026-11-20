# Low-reflecting boundary: code

C4 level 4 (code) of the component **Low-reflecting boundary** (`sms.lrb`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> Low-reflecting boundaries for solid mechanics.

**Decisions:** [ADR 42](../../../adr/0042-solid-mechanics-moves-into-nyx-sms.md)  
**Size:** 1 file, 298 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/lrb.rs` | 298 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem, linalg::csr, linalg::csr::assembly, linalg::driver, material | C, CCoo, CCsr, Coo, M6, RCsr, c, eigs_complex |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Scalar, Solution, SolveResult, cell_fields_of |
| [Typed problem](./problem.md) | problem | Problem |
| [Solid kernel](./solid.md) | solid | Setup, complex_mode_fields |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Axisymmetric study kernel](./crosssection.md) | lrb | damping_matrix, quadratic_modes |
| [Solid kernel](./solid.md) | lrb | solve_lrb |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class LrbAssembly {
      <<struct>>
      +order usize
      +d usize
      +n usize
      +kf CCsr
      +mf CCsr
      +cf Csr
      +free Vec~usize~
      +n_modes usize
      +more_1_fields
  }
  class LrbSolve {
      <<struct>>
      +omegas Vec~C~
      +vectors Vec
      +converged bool
      +dim usize
  }
  class fns_lrb_rs {
      <<module lrb.rs>>
      +voigt_speeds(c6, rho) tuple
      +damping_matrix(mesh, c6, rho, labels, order, ncomp) Csr
      +quadratic_modes(kf, mf, cf, omega_s, nev) Result
      +solve_lrb(mesh, spec, st, k, m, rho_cells, labels, n_modes, f_search, t0) Result
      +assemble_lrb(mesh, spec, st, k, m, rho_cells, labels, n_modes, f_search) Result
      +solve_lrb_modes(asm) Result
      +lrb_fields(mesh, asm, solved) Result
      +QuadraticModes type
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `voigt_speeds` | `lrb.rs:26` | `pub fn voigt_speeds(c6: &M6, rho: f64) -> (f64, f64)` —  Pressure and shear speeds of an anisotropic stiffness (Voigt averages) and density. |
| fn | `damping_matrix` | `lrb.rs:36` | `pub fn damping_matrix( mesh: &Discretisation, c6: &[M6], rho: &[f64], labels: &[i32], order: usize, ncomp: usize, ) -> Csr` —  The real symmetric damping matrix `C` over `ncomp` component blocks of the boundaries `labels`. |
| type | `QuadraticModes` | `lrb.rs:122` |  Eigenvalues, eigenvectors, whether all converged, and the Krylov dimension used. |
| fn | `quadratic_modes` | `lrb.rs:126` | `pub fn quadratic_modes(kf: &CCsr, mf: &CCsr, cf: &Csr, omega_s: f64, nev: usize) -> Result<QuadraticModes, String>` —  Eigenpairs `omega` (complex, rad/s) near `omega_s` of `(K + i omega C - omega^2 M) u = 0` over the free |
| fn | `solve_lrb` | `lrb.rs:179` | `pub fn solve_lrb( mesh: &Discretisation, spec: &Problem, st: &crate::solid::Setup, k: &Csr, m: &Csr, rho_cells: &[f64], labels: &[i32], n_modes: usize, f_search: f64, t0: std::time::Instant, ) -> Result<SolveResult, String>` —  Eigenfrequencies with low-reflecting boundaries for the real stiffness and mass matrices of |
| struct | `LrbAssembly` | `lrb.rs:199` | fields: `order`, `d`, `n`, `kf`, `mf`, `cf`, `free`, `n_modes`, `f_search` —  Stage 1 of the solid kernel with low-reflecting boundaries: stiffness, mass and the boundary damping on the free dofs |
| struct | `LrbSolve` | `lrb.rs:212` | fields: `omegas`, `vectors`, `converged`, `dim` —  Stage 2 output: complex angular frequencies and the free-dof vectors of the converged modes. |
| fn | `assemble_lrb` | `lrb.rs:220` | `pub fn assemble_lrb( mesh: &Discretisation, spec: &Problem, st: &crate::solid::Setup, k: &Csr, m: &Csr, rho_cells: &[f64], labels: &[i32], n_modes: usize, f_search: f64, ) -> Result<LrbAssembly, String>` |
| fn | `solve_lrb_modes` | `lrb.rs:253` | `pub fn solve_lrb_modes(asm: &LrbAssembly) -> Result<LrbSolve, String>` —  Stage 2: the quadratic eigenpairs nearest the search frequency. |
| fn | `lrb_fields` | `lrb.rs:263` | `pub fn lrb_fields(mesh: &Discretisation, asm: &LrbAssembly, solved: &LrbSolve) -> Result<SolveResult, String>` —  Stage 3: plot fields and per-cell complex displacement of every mode, the frequencies and Q. |
