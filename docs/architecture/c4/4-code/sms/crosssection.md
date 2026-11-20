# Axisymmetric study kernel: code

C4 level 4 (code) of the component **Axisymmetric study kernel** (`sms.crosssection`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> Axisymmetric study kernel: constraints and axis conditions, absorbing layers, complex eigensolve and fields on the (r, z) half plane (the shared cross-section stages); every mode carries numeric frequency parts, a layer fraction and its energy shares per component and domain.

**Decisions:** [ADR 43](../../../adr/0043-axisymmetric-solid-kernel-stays-a-library-entry.md), [ADR 68](../../../adr/0068-the-axisymmetric-solid-gets-a-study.md), [ADR 69](../../../adr/0069-axis-conditions-and-cylindrical-stretch-of-the-axisymmetric-solid.md), [ADR 80](../../../adr/0080-modes-carry-what-identifies-them-and-a-convergence-check-follows-a-mode.md)  
**Size:** 1 file, 504 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/cross_section.rs` | 504 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem, fields, linalg::csr, linalg::csr::assembly, linalg::driver, material, props, section | C, CCsr, CellFields, M6, PmlDir, RCsr, c, check_boundary_labels, eigs_complex, plane_stiffness |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Field, Scalar, Solution, SolveResult, cell_fields_of |
| [Axisymmetric kernel](./axisym.md) | axisym | assemble_mass_axisym, assemble_stiffness_axisym, average_about_axis, axis_constraints, cylindrical_stretches, mode_energy_axisym |
| [Low-reflecting boundary](./lrb.md) | lrb | damping_matrix, quadratic_modes |
| [Typed problem](./problem.md) | problem |  |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | cross_section | CrossModes, solve_axisymmetric, solve_solid_axisym |
| [Solid kernel](./solid.md) | cross_section | solve_axisymmetric |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class CrossModes {
      <<struct>>
      +m_az f64
      +freqs Vec~C~
      +vectors Vec
      +c6 Vec~M6~
      +rho Vec~f64~
      +converged bool
  }
  class CrossSection {
      <<enumeration>>
      Axisymmetric
  }
  class CrossAssembly {
      <<struct>>
      +kind CrossSection
      +order usize
      +n usize
      +n1 usize
      +kf CCsr
      +mf CCsr
      +cf Option~Csr~
      +free Vec~usize~
      +more_8_fields
  }
  class CrossSolve {
      <<struct>>
      +modes CrossModes
      +dim usize
  }
  class fns_cross_section_rs {
      <<module cross_section.rs>>
      +solve_solid_axisym(mesh, spec, m) Result
      +solve_axisymmetric(mesh, problem) Result
      +solve_cross_section(mesh, spec, kind) Result
      +assemble_cross_section(mesh, spec, kind) Result
      +solve_cross_section_modes(asm) Result
      +cross_section_fields(mesh, asm, solved) Result
  }
  CrossAssembly --> CrossSection
  CrossSolve --> CrossModes
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `CrossModes` | `cross_section.rs:49` | fields: `m_az`, `freqs`, `vectors`, `c6`, `rho`, `converged` —  The complex eigenvectors behind a cross-section solid result: displacement over all `3 * n_dofs(order)` dofs (component blocks), the complex |
| enum | `CrossSection` | `cross_section.rs:61` | variants: `Axisymmetric` —  What the third displacement component and the third dimension are: the azimuth of an axisymmetric body with azimuthal order `m` |
| fn | `solve_solid_axisym` | `cross_section.rs:67` | `pub fn solve_solid_axisym(mesh: &Discretisation, spec: &Problem, m: f64) -> Result<(SolveResult, CrossModes), String>` —  Elastic modes of a body of revolution with azimuthal dependence `exp(-i m theta)` on the (r, z) half |
| fn | `solve_axisymmetric` | `cross_section.rs:74` | `pub fn solve_axisymmetric(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String>` —  The axisymmetric eigenfrequency study (ADR 68): the azimuthal order is taken from the study, the result carries `physics` and `study`. |
| fn | `solve_cross_section` | `cross_section.rs:100` | `pub fn solve_cross_section( mesh: &Discretisation, spec: &Problem, kind: CrossSection, ) -> Result<(SolveResult, CrossModes), String>` |
| struct | `CrossAssembly` | `cross_section.rs:115` | fields: `kind`, `order`, `n`, `n1`, `kf`, `mf`, `cf`, `free`, `fixed`, `c6`, `rho`, `n_modes`, ... —  Stage 1 of the cross-section solid kernels: stiffness, mass and (low-reflecting) damping on the free dofs, the |
| fn | `assemble_cross_section` | `cross_section.rs:139` | `pub fn assemble_cross_section( mesh: &Discretisation, spec: &Problem, kind: CrossSection, ) -> Result<CrossAssembly, String>` |
| struct | `CrossSolve` | `cross_section.rs:276` | fields: `modes`, `dim` —  Stage 2 output: the mode shapes over all dofs and what the eigensolver says about them. |
| fn | `solve_cross_section_modes` | `cross_section.rs:284` | `pub fn solve_cross_section_modes(asm: &CrossAssembly) -> Result<CrossSolve, String>` —  Stage 2: the complex eigenpairs (`K - w^2 M` or the quadratic problem with the low-reflecting damping) nearest the |
| fn | `cross_section_fields` | `cross_section.rs:357` | `pub fn cross_section_fields( mesh: &Discretisation, asm: &CrossAssembly, solved: &CrossSolve, ) -> Result<SolveResult, String>` —  Stage 3: the plot fields `\|u\|`, `Re ux` (radial), `Re uz` (axial), `Re uy` (azimuthal), `\|uy\|` of every mode (largest component |
