# Axisymmetric kernel: code

C4 level 4 (code) of the component **Axisymmetric kernel** (`sms.axisym`) in the container **nyx-sms** (crate `nyx_sms`). [Component view of the container](../../3-components/sms.md) · [All code pages](../README.md)

> Elastic modes of a body of revolution on its half-plane cross-section.

**Decisions:** [ADR 43](../../../adr/0043-axisymmetric-solid-kernel-stays-a-library-entry.md)  
**Size:** 1 file, 232 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-sms/src/axisym.rs` | 232 | Elastic modes of a body of revolution (a ring, a disk, a wheel) on its half-plane cross-section: the displacement is `u(r, z) exp(-i m theta)` with azimuthal order `m`, and the mesh has `x = r`, `y = z`. The three components are stored a... |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem, linalg::csr, linalg::csr::assembly, material | C, CCoo, CCsr, M6, c |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root | Discretisation |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Cross-section kernel](./crosssection.md) | axisym | assemble_mass_axisym, assemble_stiffness_axisym, average_about_axis |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class fns_axisym_rs {
      <<module axisym.rs>>
      +average_about_axis(c6, axis, samples) M6
      +assemble_stiffness_axisym(mesh, c6, order, m, stretch) Result
      +assemble_mass_axisym(mesh, rho, order, stretch) CCsr
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| fn | `average_about_axis` | `axisym.rs:24` | `pub fn average_about_axis(c6: &M6, axis: usize, samples: usize) -> M6` —  The stiffness (solver Voigt order xx, yy, zz, xy, yz, xz in the frame of the cross-section) |
| fn | `assemble_stiffness_axisym` | `axisym.rs:81` | `pub fn assemble_stiffness_axisym( mesh: &Discretisation, c6: &[M6], order: usize, m: f64, stretch: Option<&[[C; 3]]>, ) -> Result<CCsr, String>` |
| fn | `assemble_mass_axisym` | `axisym.rs:189` | `pub fn assemble_mass_axisym(mesh: &Discretisation, rho: &[f64], order: usize, stretch: Option<&[[C; 3]]>) -> CCsr` —  Mass matrix `rho r` over `3 * n_dofs` unknowns, with the PML stretch determinant as a complex weight. |
