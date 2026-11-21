# Driven frequency-domain and PROM: code

C4 level 4 (code) of the component **Driven frequency-domain and PROM** (`ems.driven`) in the container **nyx-ems** (crate `nyx_ems`). [Component view of the container](../../3-components/ems.md) · [All code pages](../README.md)

> Driven 3D frequency-domain curl-curl solver with lumped/waveport excitations and Hermite Tangential Galerkin PROM sweep acceleration.

**Decisions:** [ADR 106](../../../adr/0106-driven-frequency-domain-waveports-and-hermite-prom.md)  
**Size:** 2 files, 2102 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-ems/src/driven.rs` | 1270 |  |
| `crates/nyx-ems/src/ports.rs` | 832 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | edge, fem, fields, linalg::csr, linalg::csr::assembly, linalg::driver, mesh, periodic +4 | C, CCsr, CellFields, Coo, CurvedTet, DirectSolver, Geometry, Match, RCsr, Reduction +16 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | cqed, crate root, frame, result | Discretisation, DrivenPortData, Field, FieldFrame, Scalar, Solution, SolveResult, cell_fields_of |
| [3D cavity kernel](./cavity.md) | cavity | Ned3, assemble_maxwell, basis3 |
| [Thin-film IBC and Josephson ROM](./ibc.md) | ibc, rom | AssembledJjPort, MU0, SheetCellKind, assemble_ibc_and_jj_3d, build_face_to_cell_3d, build_tls_meta, reduce_real_csr, reduce_real_vec, sorted3 |
| [Mode kernel](./modes.md) | crate root, modes | C0, EpsCell, assemble_modes_with_bcs, solve_mode_pairs |
| [Typed problem](./problem.md) | problem | LumpedPortBc, Problem, Study, SurfaceImpedanceBc, WavePortBc |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | driven, ports | AssembledLumpedPort, AssembledWavePort, DrivenAssembly, assemble_emw_driven, assemble_lumped_ports_3d, assemble_waveports_3d, solve_driven |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace driven_rs {
    class DrivenPortKind {
        <<enumeration>>
        Lumped
        Wave
    }
    class DrivenPortEntry {
        <<struct>>
        +name String
        +kind DrivenPortKind
        +b_hat Vec~f64~
        +excitation_coeffs(omega, waveports) tuple
    }
    class DrivenAssembly {
        <<struct>>
        +order usize
        +n usize
        +nred usize
        +ned Ned3
        +red Reduction
        +f_min f64
        +f_max f64
        +n_points usize
        +more_10_fields
        +term_coeffs(omega) Vec
        +assemble_matrix_at(omega) CCsr
        +apply_matrix_deriv(omega, order_deriv, x) Vec~C~
    }
    class fns_driven_rs {
        <<module driven.rs>>
        +assemble_emw_driven(mesh, problem) Result
        +solve_driven(mesh, problem) Result
    }
  }
  namespace ports_rs {
    class AssembledLumpedPort {
        <<struct>>
        +bc LumpedPortBc
        +d_hat slice
        +w_hat slice
        +n_hat slice
        +length_m f64
        +width_m f64
        +triangles Vec
        +b_full Vec~f64~
        +more_2_fields
    }
    class AssembledWavePort {
        <<struct>>
        +bc WavePortBc
        +n_out slice
        +u_hat slice
        +v_hat slice
        +triangles Vec
        +b_full Vec~f64~
        +b_red Vec~f64~
        +w_full RCsr
        +more_2_fields
        +beta_and_derivs(omega) tuple
    }
    class fns_ports_rs {
        <<module ports.rs>>
        +eval_triangle_tangential_basis_ordered(mesh, ned, geo, face_to_cell, tri_sorted) Resultusizef64f64
        +assemble_lumped_ports_3d(mesh, ned, pec_labels, lumped_bcs) Result
        +assemble_waveports_3d(mesh, ned, pec_labels, pmc_labels, ibc_bcs, waveport_bcs, cell_eps_3d, cell_mur_3d, f_min, f_max) Result
    }
  }
  DrivenPortEntry --> DrivenPortKind
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| enum | `DrivenPortKind` | `driven.rs:110` | variants: `Lumped`, `Wave` —  Kind of driven port in the unified port list `0..P`. |
| struct | `DrivenPortEntry` | `driven.rs:117` | fields: `name`, `kind`, `b_hat` —  Unified metadata for one driven port $p \in \{0..P-1\}$. |
| struct | `DrivenAssembly` | `driven.rs:152` | fields: `order`, `n`, `nred`, `ned`, `red`, `f_min`, `f_max`, `n_points`, `prom_tol`, `d_ruiz`, `terms_hat`, `ports`, ... —  Stage 1 of the 3D driven frequency-domain solver: Nedelec numbering, PEC/periodic reduction, |
| fn | `assemble_emw_driven` | `driven.rs:226` | `pub fn assemble_emw_driven(mesh: &Discretisation, problem: &Problem) -> Result<DrivenAssembly, String>` —  Assembles Stage 1 (`DrivenAssembly`) for a 3D `Study::Driven` problem. |
| fn | `solve_driven` | `driven.rs:1024` | `pub fn solve_driven(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String>` —  Solves a 3D `Study::Driven` frequency sweep via Hermite Tangential Galerkin PROM (`prom_tol = Some(tol)`) |
| fn | `DrivenPortEntry::excitation_coeffs` | `driven.rs:128` | `pub fn excitation_coeffs(&self, omega: f64, waveports: &[AssembledWavePort]) -> (f64, C, C, C)` —  Returns `(z0_p(omega), gamma_p(omega), d(gamma_p)/d(omega), d^2(gamma_p)/d(omega)^2)` at angular frequency `omega` (rad/s), |
| fn | `DrivenAssembly::term_coeffs` | `driven.rs:177` | `pub fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)>` —  Evaluates the complex scalar coefficients $(\alpha_m(\omega), \alpha_m'(\omega), \alpha_m''(\omega))$ for each term $\hat{\mathbf{A}}_m$. |
| fn | `DrivenAssembly::assemble_matrix_at` | `driven.rs:196` | `pub fn assemble_matrix_at(&self, omega: f64) -> CCsr` —  Assembles the full equilibrated matrix $\hat{\mathbf{A}}(\omega) = \sum_m \alpha_m(\omega) \hat{\mathbf{A}}_m$. |
| fn | `DrivenAssembly::apply_matrix_deriv` | `driven.rs:208` | `pub fn apply_matrix_deriv(&self, omega: f64, order_deriv: usize, x: &[C]) -> Vec<C>` —  Computes $\mathbf{y} = \hat{\mathbf{A}}'(\omega) \mathbf{x}$ or $\hat{\mathbf{A}}''(\omega) \mathbf{x}$. |
| struct | `AssembledLumpedPort` | `ports.rs:55` | fields: `bc`, `d_hat`, `w_hat`, `n_hat`, `length_m`, `width_m`, `triangles`, `b_full`, `b_red`, `c_port_full` —  Assembled operators and geometric metadata for one `lumped_port` terminal port. |
| struct | `AssembledWavePort` | `ports.rs:73` | fields: `bc`, `n_out`, `u_hat`, `v_hat`, `triangles`, `b_full`, `b_red`, `w_full`, `w_red`, `beta2_poly` —  Assembled operators, 2D modal dispersion, and geometric metadata for one `waveport`. |
| fn | `eval_triangle_tangential_basis_ordered` | `ports.rs:118` | `pub fn eval_triangle_tangential_basis_ordered( mesh: &Discretisation, ned: &Ned3, geo: &femcore::fem::Geometry, face_to_cell: &HashMap<[u32; 3], (usize, [usize; 3])>, tri_sorted: [u32; 3], ) -> Result<(usize, f64, [f64` —  Evaluates the 3D tangential Nedelec basis vectors `psi_{s,\|\|} = w_s - (w_s . n) n` at the 6-point Dunavant |
| fn | `assemble_lumped_ports_3d` | `ports.rs:197` | `pub fn assemble_lumped_ports_3d( mesh: &Discretisation, ned: &Ned3, pec_labels: &[i32], lumped_bcs: &[LumpedPortBc], ) -> Result<Vec<AssembledLumpedPort>, String>` —  Assembles all 3D `lumped_port` boundary conditions following the single-facet rule. |
| fn | `assemble_waveports_3d` | `ports.rs:490` | `pub fn assemble_waveports_3d( mesh: &Discretisation, ned: &Ned3, pec_labels: &[i32], pmc_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], waveport_bcs: &[WavePortBc], cell_eps_3d: &[[[f64; 3]; 3]], cell_mur_3d: &[f64], f_min: f64, f_max: f64, ) -> Result<Vec<AssembledWavePort>, String>` —  Assembles all 3D `waveport` boundary conditions by extracting the 2D boundary mesh on each waveport plane |
| fn | `AssembledWavePort::beta_and_derivs` | `ports.rs:94` | `pub fn beta_and_derivs(&self, omega: f64) -> (f64, f64, f64, f64)` —  Evaluates `(beta_p(omega), d(beta_p)/d(omega), d^2(beta_p)/d(omega)^2, z0_p(omega))` at angular frequency `omega` (rad/s). |
