# Thin-film IBC and Josephson ROM: code

C4 level 4 (code) of the component **Thin-film IBC and Josephson ROM** (`ems.ibc`) in the container **nyx-ems** (crate `nyx_ems`). [Component view of the container](../../3-components/ems.md) · [All code pages](../README.md)

> Superconducting thin-film surface impedance boundary condition (Leontovich/two-fluid, curved Jacobian integrals), rank-1 lumped Josephson port contribution, passive CircuitRom, and TLS metadata.

**Decisions:** [ADR 104](../../../adr/0104-bcs-surface-impedance-and-lumped-josephson-ports.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md)  
**Size:** 2 files, 3222 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-ems/src/ibc.rs` | 1758 |  |
| `crates/nyx-ems/src/rom.rs` | 1464 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | edge, edge_mixed, fem, fem_mixed, linalg::amg_elastic, linalg::csr, linalg::csr::assembly, linalg::dense +6 | C, CCoo, CCsr, CellType, Coo, CurvedTet, Dense, DirectSolver, FacetType, Geometry +36 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | cqed, crate root | CircuitRom, Discretisation, ModeEnergyRow, RomPort, TlsMeta |
| [3D cavity kernel](./cavity.md) | cavity | CURVED_QUAD_DEGREE, EmwAssembly, Ned3, basis3_keyed, gradient_full, safe_curved_jacobian |
| [Native mixed-cell EM](./mixed.md) | mixed_ops | gradient_mixed_full, scalar_nullspace_count_mixed |
| [Mode kernel](./modes.md) | crate root | C0 |
| [Typed problem](./problem.md) | problem | LumpedJjBc, Problem, SurfaceImpedanceBc, Symmetry, symmetry_facet_flags_mixed |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | rom | rom_complex_eigenfrequencies, rom_eigenfrequencies |
| [Axisymmetric kernel](./axisym.md) | ibc | add_ibc_2d_axisym, add_ibc_2d_axisym_mixed, pec_edges_and_nodes_2d |
| [3D cavity kernel](./cavity.md) | ibc, rom | AssembledJjPort, EPS0, MU0, PerturbativeIbcEvaluator, PerturbativeIbcTri, SheetCellKind, assemble_ibc_and_jj_3d, assemble_ibc_and_jj_3d_cavity +10 |
| [Driven frequency-domain and PROM](./driven.md) | crate root, ibc, rom | AssembledJjPort, MU0, SheetCellKind, assemble_ibc_and_jj_3d, build_face_to_cell_3d, build_tls_meta, reduce_real_csr, reduce_real_vec +1 |
| [Native mixed-cell EM](./mixed.md) | ibc | AssembledJjPort, ETA0, IbcJjAssembly3D, MU0, PerturbativeIbcTri, SheetCellKind, ibc_parallel_params, jj_inv_l_tot +1 |
| [Mode kernel](./modes.md) | ibc | add_ibc_2d_waveguide, add_ibc_2d_waveguide_mixed, pec_edges_and_nodes_2d, pec_edges_and_nodes_2d_mixed |
| [Complex-mode kernel](./modescx.md) | ibc | EPS0, MU0, add_ibc_2d_waveguide, add_ibc_2d_waveguide_mixed, pec_edges_and_nodes_2d, pec_edges_and_nodes_2d_mixed |
| [Typed problem](./problem.md) | ibc | resolve_jj_l_and_ej |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace ibc_rs {
    class AssembledJjPort {
        <<struct>>
        +bc LumpedJjBc
        +d_hat slice
        +w_hat slice
        +n_hat slice
        +length_m f64
        +width_m f64
        +u_min f64
        +u_max f64
        +more_6_fields
    }
    class SheetCellKind {
        <<enumeration>>
        Ibc
        Jj
        Pec
    }
    class PerturbativeIbcTri {
        <<struct>>
        +tri slice
        +l_k_sq f64
        +r_s f64
        +f_ref f64
    }
    class IbcJjAssembly3D {
        <<struct>>
        +pec_tris Vec
        +ibc_tris Vec
        +pert_ibc_tris Vec~PerturbativeIbcTri~
        +direct_ibc_labels Vec~i32~
        +jj_tris Vec
        +k_ibc_full RCsr
        +c_ibc_full RCsr
        +kc_ibc_full CCsr
        +more_5_fields
    }
    class PerturbativeIbcEvaluator {
        <<struct>>
        +new(mesh, ned, geo, mur_phys, pert_tris) Result
        +new_mixed(mesh, ned, mur_phys, pert_tris) Result
        +eval_mode(x_full) tuple
        +eval_rom_grams(q_full_cols) tuple
    }
    class fns_ibc_rs {
        <<module ibc.rs>>
        +jj_flux_scale() f64
        +resolve_jj_l_and_ej(l_j, e_j_hz) tuple
        +jj_inv_l_tot(l_j, e_l_hz) f64
        +ibc_parallel_params(r_s, l_k_sq, omega_ref) tuple
        +sorted3(a, b, c) slice
        +build_face_to_cells_3d(mesh) HashMapu32
        +assemble_ibc_and_jj_3d(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f0) Result
        +assemble_ibc_and_jj_3d_cavity(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f0) Result
        +reduce_real_csr(red, pmat, a) RCsr
        +reduce_real_vec(red, b_full) Vec~f64~
        +more_9_functions()
        +MU0 const
        +EPS0 const
        +ETA0 const
        +PHI0 const
        +more_1_items
    }
  }
  namespace rom_rs {
    class fns_rom_rs {
        <<module rom.rs>>
        +compute_mode_energy_and_scale(f, x_red, rot, asm) tuple
        +build_circuit_rom(mesh, problem, asm, solved_vectors, freq_indices, modes, pert_eval) Result
        +build_circuit_rom_mixed(mesh, ned, problem, asm, solved_vectors, freq_indices, modes, pert_eval) Result
        +rom_eigenfrequencies(rom) Result
        +rom_complex_eigenfrequencies(rom) Result
        +build_tls_meta(mesh, problem) TlsMeta
    }
  }
  IbcJjAssembly3D --> PerturbativeIbcTri
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `MU0` | `ibc.rs:25` | `f64` —  Vacuum magnetic permeability `mu_0` in H/m. |
| const | `EPS0` | `ibc.rs:27` | `f64` —  Vacuum electric permittivity `eps_0 = 1 / (mu_0 c_0^2)` in F/m, consistent with `MU0` and `C0` to machine precision. |
| const | `ETA0` | `ibc.rs:29` | `f64` —  Vacuum wave impedance `eta_0 = mu_0 c_0` in Ohm. |
| const | `PHI0` | `ibc.rs:31` | `f64` —  Superconducting magnetic flux quantum `Phi_0 = h / (2e)` in Wb. |
| const | `PLANCK_H` | `ibc.rs:33` | `f64` —  Planck constant `h` in J*s. |
| fn | `jj_flux_scale` | `ibc.rs:39` | `pub fn jj_flux_scale() -> f64` —  Conversion scale `(Phi_0 / (2 pi))^2 / h` relating `L_J` (H) and `E_J / h` (Hz): `L_J * E_J_Hz = jj_flux_scale()`. |
| fn | `resolve_jj_l_and_ej` | `ibc.rs:46` | `pub fn resolve_jj_l_and_ej(l_j: f64, e_j_hz: f64) -> (f64, f64)` —  Resolves `(L_j, E_j_Hz)` when one is given and the other is `<= 0.0`. |
| fn | `jj_inv_l_tot` | `ibc.rs:59` | `pub fn jj_inv_l_tot(l_j: f64, e_l_hz: f64) -> f64` —  Total inverse linear inductance `1 / L_tot = 1 / L_J + 1 / L_p` (1/H) of a lumped JJ port. |
| fn | `ibc_parallel_params` | `ibc.rs:69` | `pub fn ibc_parallel_params(r_s: f64, l_k_sq: f64, omega_ref: f64) -> (f64, f64)` —  Series-to-parallel conversion of `Z_s(omega_ref) = R_s + i omega_ref L_k_sq`: |
| fn | `sorted3` | `ibc.rs:84` | `pub fn sorted3(a: u32, b: u32, c: u32) -> [u32` |
| struct | `AssembledJjPort` | `ibc.rs:118` | fields: `bc`, `d_hat`, `w_hat`, `n_hat`, `length_m`, `width_m`, `u_min`, `u_max`, `inv_l_tot`, `triangles`, `b_full`, `b_red`, ... —  Assembled operators and geometric metadata for one `lumped_jj` terminal port. |
| enum | `SheetCellKind` | `ibc.rs:141` | variants: `Ibc`, `Jj`, `Pec` —  Sheet surface-current evaluator rule attached to a 3D cell face. |
| struct | `PerturbativeIbcTri` | `ibc.rs:149` | fields: `tri`, `l_k_sq`, `r_s`, `f_ref` —  One boundary triangle handled by perturbative surface-current IBC (`ibc_mode = "perturbative"` or `"auto"`, FR-30). |
| struct | `IbcJjAssembly3D` | `ibc.rs:157` | fields: `pec_tris`, `ibc_tris`, `pert_ibc_tris`, `direct_ibc_labels`, `jj_tris`, `k_ibc_full`, `c_ibc_full`, `kc_ibc_full`, `k_jj_full`, `m_jj_full`, `jj_ports`, `cell_sheet`, ... —  Combined 3D assembly of `pec`, `surface_impedance`, and `lumped_jj` boundary conditions. |
| fn | `build_face_to_cells_3d` | `ibc.rs:190` | `pub fn build_face_to_cells_3d(mesh: &Discretisation) -> HashMap<[u32` —  Maps each sorted 3D triangle `[v0, v1, v2]` to all adjacent tetrahedra `(ci, lf)` (1 for exterior faces, |
| struct | `PerturbativeIbcEvaluator` | `ibc.rs:209` |  Precomputed quadrature-point surface-current basis on perturbative `surface_impedance` triangles (FR-30). |
| fn | `assemble_ibc_and_jj_3d` | `ibc.rs:703` | `pub fn assemble_ibc_and_jj_3d( mesh: &Discretisation, ned: &Ned3, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], jj_bcs: &[LumpedJjBc], f0: f64, ) -> Result<IbcJjAssembly3D, String>` —  Deduplicates and assembles all 3D `pec`, `surface_impedance`, and `lumped_jj` boundary conditions |
| fn | `assemble_ibc_and_jj_3d_cavity` | `ibc.rs:716` | `pub fn assemble_ibc_and_jj_3d_cavity( mesh: &Discretisation, ned: &Ned3, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], jj_bcs: &[LumpedJjBc], f0: f64, ) -> Result<IbcJjAssembly3D, String>` —  Deduplicates and assembles all 3D `pec`, `surface_impedance`, and `lumped_jj` boundary conditions |
| fn | `reduce_real_csr` | `ibc.rs:1086` | `pub fn reduce_real_csr(red: &femcore::periodic::Reduction, pmat: &CCsr, a: &RCsr) -> RCsr` —  Real congruence transform `P^H A P` for a real matrix `A` under a periodic/Dirichlet reduction `red`. |
| fn | `reduce_real_vec` | `ibc.rs:1094` | `pub fn reduce_real_vec(red: &femcore::periodic::Reduction, b_full: &[f64]) -> Vec<f64>` —  Restricts a full real load vector `b_full` (`n`) to the reduced space `b_red = P^T b_full` (`nred`). |
| fn | `ruiz_equilibrate` | `ibc.rs:1110` | `pub fn ruiz_equilibrate(kc: &CCsr, mc: &CCsr, g_red: &CCsr, sigma: f64) -> (Vec<f64>, CCsr, CCsr, CCsr, CCsr)` —  Computes iterative symmetric Ruiz diagonal scaling factors `d` on `\|K\| + sigma \|M\|`, |
| fn | `ruiz_unscale_and_rayleigh` | `ibc.rs:1209` | `pub fn ruiz_unscale_and_rayleigh(d: &[f64], y: &[C], kc: &CCsr, mc: &CCsr) -> (Vec<C>, C)` —  Unscales a Ruiz-equilibrated eigenvector `x = D y` and computes its exact Rayleigh-quotient |
| fn | `ruiz_unscale_and_rayleigh_damped` | `ibc.rs:1225` | `pub fn ruiz_unscale_and_rayleigh_damped( d: &[f64], y: &[C], kc_no_damp: &CCsr, c_damp: &RCsr, mc: &CCsr, ) -> (Vec<C>, C)` —  Unscales a Ruiz-equilibrated eigenvector `x = D y` and evaluates the mode-exact damped Rayleigh quotient |
| fn | `pec_edges_and_nodes_2d` | `ibc.rs:1288` | `pub fn pec_edges_and_nodes_2d( mesh: &Discretisation, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], pmc_labels: &[i32], ) -> (Vec<usize>, Vec<usize>)` —  Collects PEC boundary edge indices `be` and vertex indices `bn` for a 2D mesh, honoring explicit `pec` |
| fn | `add_ibc_2d_waveguide` | `ibc.rs:1319` | `pub fn add_ibc_2d_waveguide( mesh: &Discretisation, order: usize, freq: f64, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], a_coo: &mut CCoo, b_coo: &mut CCoo, ) -> Result<(), String>` —  Assembles 2D `surface_impedance` contributions into `(a_coo, b_coo)` for `modes` / `modes_complex`. |
| fn | `add_ibc_2d_axisym` | `ibc.rs:1417` | `pub fn add_ibc_2d_axisym( mesh: &Discretisation, order: usize, f_search: f64, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], a_coo: &mut CCoo, ) -> Result<(), String>` —  Assembles 2D `surface_impedance` contributions into `a_coo` for `axisymmetric_modes` (`(r, z)` half-plane). |
| fn | `pec_edges_and_nodes_2d_mixed` | `ibc.rs:1518` | `pub fn pec_edges_and_nodes_2d_mixed( mesh: &MixedDiscretisation, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], pmc_labels: &[i32], symmetry: Symmetry, ) -> Result<(Vec<usize>, Vec<usize>), String>` |
| fn | `add_ibc_2d_waveguide_mixed` | `ibc.rs:1565` | `pub fn add_ibc_2d_waveguide_mixed( mesh: &MixedDiscretisation, ned: &NedMixed, phi_map: &MixedLagrangeDofMap, freq: f64, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], a_coo: &mut CCoo, b_coo: &mut CCoo, ) -> Result<(), String>` |
| fn | `add_ibc_2d_axisym_mixed` | `ibc.rs:1662` | `pub fn add_ibc_2d_axisym_mixed( mesh: &MixedDiscretisation, ned: &NedMixed, phi_map: &MixedLagrangeDofMap, f_search: f64, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], a_coo: &mut CCoo, ) -> Result<(), String>` |
| fn | `PerturbativeIbcEvaluator::new` | `ibc.rs:215` | `pub fn new( mesh: &Discretisation, ned: &Ned3, geo: &femcore::fem::Geometry, mur_phys: &[f64], pert_tris: &[PerturbativeIbcTri], ) -> Result<Self, String>` |
| fn | `PerturbativeIbcEvaluator::new_mixed` | `ibc.rs:317` | `pub fn new_mixed( mesh: &MixedDiscretisation, ned: &NedMixed, mur_phys: &[f64], pert_tris: &[PerturbativeIbcTri], ) -> Result<Self, String>` —  Mixed-cell (`Tet + Prism + Hex`) constructor for perturbative `surface_impedance` triangles. |
| fn | `PerturbativeIbcEvaluator::eval_mode` | `ibc.rs:452` | `pub fn eval_mode(&self, x_full: &[C]) -> (f64, f64)` —  Evaluates `(sum_T (L_k_sq / mu_0) int_T \|K_tilde\|^2 dS, sum_T R_s int_T \|K_tilde\|^2 dS)` |
| fn | `PerturbativeIbcEvaluator::eval_rom_grams` | `ibc.rs:473` | `pub fn eval_rom_grams(&self, q_full_cols: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>)` —  Evaluates the symmetric `r x r` Gram matrices `(dk_kin, g_r_mat)` on full-DOF real basis columns `q_full_cols`: |
| fn | `compute_mode_energy_and_scale` | `rom.rs:80` | `pub fn compute_mode_energy_and_scale(f: C, x_red: &[C], rot: C, asm: &EmwAssembly) -> (f64, ModeEnergyRow, f64)` —  Computes the physical scaling factor `s_E > 0` such that `U_elec,m = (1/4) eps_0 (s_E x_m)^H M_tot (s_E x_m) = 1 J`, |
| fn | `build_circuit_rom` | `rom.rs:562` | `pub fn build_circuit_rom( mesh: &Discretisation, problem: &Problem, asm: &EmwAssembly, solved_vectors: &[Vec<C>], freq_indices: &[(C, usize)], modes: Vec<ModeEnergyRow>, pert_eval: Option<&PerturbativeIbcEvaluator>, ) -> Result<CircuitRom, String>` —  Assembles the passive Craig-Bampton reduced-order model `CircuitRom` (`Q = orth_{M_0}([Q_static, Q_dyn])`). |
| fn | `build_circuit_rom_mixed` | `rom.rs:576` | `pub fn build_circuit_rom_mixed( mesh: &MixedDiscretisation, ned: &NedMixed, problem: &Problem, asm: &EmwAssembly, solved_vectors: &[Vec<C>], freq_indices: &[(C, usize)], modes: Vec<ModeEnergyRow>, pert_eval: Option<&PerturbativeIbcEvaluator>, ) -> Result<CircuitRom, String>` —  Assembles `CircuitRom` on a native 3D `MixedDiscretisation`. |
| fn | `rom_eigenfrequencies` | `rom.rs:1275` | `pub fn rom_eigenfrequencies(rom: &CircuitRom) -> Result<Vec<f64>, String>` —  Evaluates the linearized resonant eigenfrequencies (Hz, ascending) of a `CircuitRom` by solving the |
| fn | `rom_complex_eigenfrequencies` | `rom.rs:1313` | `pub fn rom_complex_eigenfrequencies(rom: &CircuitRom) -> Result<Vec<C>, String>` —  Evaluates the complex damped resonant eigenfrequencies (Hz, positive real part, ascending by `Re(f)`) |
| fn | `build_tls_meta` | `rom.rs:1384` | `pub fn build_tls_meta(mesh: &Discretisation, problem: &Problem) -> TlsMeta` —  Builds `TlsMeta` for dielectric TLS surface participation analysis on a 3D `ems` result. |
