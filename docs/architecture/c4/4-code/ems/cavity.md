# 3D cavity kernel: code

C4 level 4 (code) of the component **3D cavity kernel** (`ems.cavity`) in the container **nyx-ems** (crate `nyx_ems`). [Component view of the container](../../3-components/ems.md) · [All code pages](../README.md)

> Electromagnetic eigenmodes of 3D cavities and resonators, with superconducting thin-film IBC, lumped Josephson junctions, and physical H/K export.

**Decisions:** [ADR 45](../../../adr/0045-ems-moves-into-nyx-ems.md), [ADR 104](../../../adr/0104-bcs-surface-impedance-and-lumped-josephson-ports.md)  
**Size:** 2 files, 2145 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-ems/src/bloch.rs` | 630 |  |
| `crates/nyx-ems/src/cavity.rs` | 1515 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | fem, fields, linalg::amg_elastic, linalg::arnoldi_real, linalg::csr, linalg::csr::assembly, linalg::driver, par +5 | C, CCoo, CCsr, CellFields, Coo, CurvedTet, DirectSolver, EigOpts, Match, PmlDir +22 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, frame, result | Discretisation, Field, FieldFrame, Scalar, Solution, SolveResult, cell_fields_of |
| [Driven frequency-domain and PROM](./driven.md) | ports | AssembledLumpedPort, assemble_lumped_ports_3d, validate_port_contacts_3d |
| [Thin-film IBC and Josephson ROM](./ibc.md) | ibc, rom | AssembledJjPort, EPS0, MU0, SheetCellKind, assemble_ibc_and_jj_3d, build_circuit_rom, build_tls_meta, compute_mode_energy_and_scale, reduce_real_csr, reduce_real_vec +3 |
| [Mode kernel](./modes.md) | crate root | C0 |
| [Optical mode sets](./optical.md) | optical | absorbing_domains |
| [Typed problem](./problem.md) | problem | Problem, Study, symmetry_facet_flags |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | bloch, cavity | BlochAssembly, BlochSolve, EmwAssembly, EmwSolve, assemble_bloch, assemble_emw_eigen, bloch_fields_cells, emw_eigen_fields +4 |
| [Driven frequency-domain and PROM](./driven.md) | cavity | Ned3, assemble_maxwell, basis3 |
| [Thin-film IBC and Josephson ROM](./ibc.md) | cavity | EmwAssembly, Ned3, basis3, gradient_full |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace bloch_rs {
    class BlochAssembly {
        <<struct>>
        +order usize
        +n usize
        +nred usize
        +ned Ned3
        +red Reduction
        +a0_red CCsr
        +a1_red CCsr
        +a2_red CCsr
        +more_14_fields
    }
    class BlochSolve {
        <<struct>>
        +modes Vec
        +vectors Vec
        +residuals Vec~f64~
        +dim usize
        +converged bool
    }
    class fns_bloch_rs {
        <<module bloch.rs>>
        +solve_bloch(mesh, problem) Result
        +assemble_bloch(mesh, problem) Result
        +solve_bloch_eigen(asm) Result
        +bloch_fields_cells(mesh, problem, asm, solved) Result
    }
  }
  namespace cavity_rs {
    class Ned3 {
        <<struct>>
        +order usize
        +n_edges usize
        +n_faces usize
        +faces Vec
        +cell_dofs Vec
        +vkey Vec~u32~
        +new(mesh, order, vkey) Ned3
        +n_dofs() usize
        +keys(mesh, ci) slice
        +face_dofs(v) slice
        +nb() usize
        +boundary_dofs(mesh, v) Vec~usize~
    }
    class EmwAssembly {
        <<struct>>
        +order usize
        +n usize
        +nred usize
        +ned Ned3
        +red Reduction
        +kc CCsr
        +kc_no_damp CCsr
        +mc CCsr
        +more_26_fields
    }
    class EmwSolve {
        <<struct>>
        +freqs Vec
        +vectors Vec
        +dim usize
        +converged bool
    }
    class Mode3dMetrics {
        <<struct>>
        +vg f64
        +comp_shares slice
        +domain_shares Vec
        +layer_fraction f64
    }
    class fns_cavity_rs {
        <<module cavity.rs>>
        +basis3(order, lam, g, gid, w, curl)
        ~assemble_complex(nc, n, f) CCsr
        +assemble_maxwell(mesh, ned, eps, mu_inv) tuple
        ~gradient_full(mesh, ned) Csr
        ~vertex_keys(mesh, matches) Vec~u32~
        ~ned_relations(mesh, ned, m, phase) Result
        +solve_cavity(mesh, problem) Result
        +assemble_emw_eigen(mesh, problem) Result
        +solve_emw_modes(asm) Result
        +emw_eigen_fields(mesh, problem, asm, solved) Result
        +more_1_functions()
        ~C3 type
        +CURVED_QUAD_DEGREE const
    }
  }
  BlochAssembly --> Ned3
  EmwAssembly --> Ned3
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `BlochAssembly` | `bloch.rs:38` | fields: `order`, `n`, `nred`, `ned`, `red`, `a0_red`, `a1_red`, `a2_red`, `freq`, `k0`, `lambda0`, `period_y`, ... —  Stage 1 of the 3D fixed-frequency Bloch kernel: Nedelec numbering, reduced quadratic pencil `(A0, A1, A2)` in `nu = k_y / k_0`, |
| struct | `BlochSolve` | `bloch.rs:64` | fields: `modes`, `vectors`, `residuals`, `dim`, `converged` —  Stage 2 output: converged complex Bloch indices `(n_B, vec_index)` sorted by distance to `sigma`. |
| fn | `solve_bloch` | `bloch.rs:73` | `pub fn solve_bloch(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String>` —  Solves a 3D `Study::Bloch` problem through its three stages (`assemble_bloch`, `solve_bloch_eigen`, `bloch_fields_cells`). |
| fn | `assemble_bloch` | `bloch.rs:201` | `pub fn assemble_bloch(mesh: &Discretisation, problem: &Problem) -> Result<BlochAssembly, String>` —  Stage 1: validates `y`-periodicity and transverse PMLs, assembles `(A0, A1, A2)` on `Ned3`, and reduces them via `P^H (.) P` |
| fn | `solve_bloch_eigen` | `bloch.rs:396` | `pub fn solve_bloch_eigen(asm: &BlochAssembly) -> Result<BlochSolve, String>` —  Stage 2: solves `(A0 + nu A1 + nu^2 A2) u = 0` near `sigma`, filters true QEP relative residuals `< 1e-5`, |
| fn | `bloch_fields_cells` | `bloch.rs:423` | `pub fn bloch_fields_cells( mesh: &Discretisation, problem: &Problem, asm: &BlochAssembly, solved: &BlochSolve, ) -> Result<SolveResult, String>` —  Stage 3: reconstructs physical `E = e_tilde exp(-j k_y (y - y_min))` and |
| struct | `Ned3` | `cavity.rs:49` | fields: `order`, `n_edges`, `n_faces`, `faces`, `cell_dofs`, `vkey` —  Global numbering of the degrees of freedom. |
| fn | `basis3` | `cavity.rs:177` | `pub fn basis3( order: usize, lam: &[f64; 4], g: &[[f64; 3]; 4], gid: &[u32; 4], w: &mut Vec<[f64; 3]>, curl: &mut Vec<[f64; 3]>, )` —  Basis functions and their curls at the barycentric point `lam` of a tetrahedron with vertex |
| type | `C3` | `cavity.rs:217` |  |
| fn | `assemble_complex` | `cavity.rs:219` | `pub(crate) fn assemble_complex<F>(nc: usize, n: usize, f: F) -> CCsr` |
| const | `CURVED_QUAD_DEGREE` | `cavity.rs:242` | `usize` —  Degree of the quadrature rule on cells with a curved edge (ADR 82): the integrand is rational in the barycentric coordinates there. |
| fn | `assemble_maxwell` | `cavity.rs:286` | `pub fn assemble_maxwell(mesh: &Discretisation, ned: &Ned3, eps: &[C3], mu_inv: &[[C; 3]]) -> (CCsr, CCsr)` —  The curl-curl matrix `K` (weights `1/mu` per axis) and the mass matrix `M` (permittivity tensor) |
| fn | `gradient_full` | `cavity.rs:337` | `pub(crate) fn gradient_full(mesh: &Discretisation, ned: &Ned3) -> Csr` —  The discrete gradient over all nodes: one column per vertex (the Whitney edge coefficients |
| fn | `vertex_keys` | `cavity.rs:367` | `pub(crate) fn vertex_keys(mesh: &Discretisation, matches: &[femcore::periodic::Match]) -> Vec<u32>` —  Equivalence classes of the vertices of periodic boundaries: the key of a vertex is the smallest |
| fn | `ned_relations` | `cavity.rs:391` | `pub(crate) fn ned_relations( mesh: &Discretisation, ned: &Ned3, m: &femcore::periodic::Match, phase: f64, ) -> Result<Vec<(usize, usize, C)>, String>` —  Relations between the dofs of the two sides of a periodic pair (all one to one because the basis |
| fn | `solve_cavity` | `cavity.rs:433` | `pub fn solve_cavity(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String>` —  Eigenfrequencies of the cavity for a `Cavity` study: the three stages in order. The result carries the complex field of every mode and |
| struct | `EmwAssembly` | `cavity.rs:444` | fields: `order`, `n`, `nred`, `ned`, `red`, `kc`, `kc_no_damp`, `mc`, `g_red`, `g_h`, `f0`, `sigma`, ... —  Stage 1 of the 3D EM eigenfrequency kernel: Nedelec numbering, `K` and `M` (PML stretch, anisotropy, IBC, JJ) reduced by the PEC |
| struct | `EmwSolve` | `cavity.rs:482` | fields: `freqs`, `vectors`, `dim`, `converged` —  Stage 2 output: complex frequencies (Hz) with the index of their reduced vector, sorted by real part. |
| fn | `assemble_emw_eigen` | `cavity.rs:489` | `pub fn assemble_emw_eigen(mesh: &Discretisation, problem: &Problem) -> Result<EmwAssembly, String>` |
| fn | `solve_emw_modes` | `cavity.rs:962` | `pub fn solve_emw_modes(asm: &EmwAssembly) -> Result<EmwSolve, String>` —  Stage 2: the eigenpairs `K x = k0^2 M x` nearest the shift with the gradients projected out (real symmetric route with |
| fn | `emw_eigen_fields` | `cavity.rs:1132` | `pub fn emw_eigen_fields( mesh: &Discretisation, problem: &Problem, asm: &EmwAssembly, solved: &EmwSolve, ) -> Result<SolveResult, String>` —  Stage 3: plot fields and per-cell complex E and H (plus sheet K when sheets are present) of every mode at the quadratic nodes, |
| struct | `Mode3dMetrics` | `cavity.rs:1385` | fields: `vg`, `comp_shares`, `domain_shares`, `layer_fraction` —  3D mode energy shares, PML layer fraction, and (when `compute_vg` is true) single-solve Poynting `v_g` |
| fn | `compute_3d_mode_metrics` | `cavity.rs:1393` | `pub(crate) fn compute_3d_mode_metrics( mesh: &Discretisation, ned: &Ned3, geo: &Geometry, x: &[C], ky_envelope: C, omega: f64, pml_cells: &[bool], eps_phys_re: &[[[f64; 3]; 3]], d_omega_eps: &[[[f64; 3]; 3]], mur_phys: &[f64], compute_vg: bool, ) -> Mode3dMetrics` |
| fn | `Ned3::new` | `cavity.rs:71` | `pub fn new(mesh: &Discretisation, order: usize, vkey: Vec<u32>) -> Ned3` |
| fn | `Ned3::n_dofs` | `cavity.rs:112` | `pub fn n_dofs(&self) -> usize` |
| fn | `Ned3::keys` | `cavity.rs:121` | `pub fn keys(&self, mesh: &Discretisation, ci: usize) -> [u32` —  Ordering keys of the vertices of cell `ci`. |
| fn | `Ned3::face_dofs` | `cavity.rs:127` | `pub fn face_dofs(&self, v: [u32; 3]) -> [usize` —  Dofs of the functions of face `v` (sorted triple of vertex numbers). |
| fn | `Ned3::nb` | `cavity.rs:132` | `pub fn nb(&self) -> usize` |
| fn | `Ned3::boundary_dofs` | `cavity.rs:142` | `pub fn boundary_dofs(&self, mesh: &Discretisation, v: [u32; 3]) -> Vec<usize>` —  Dofs of the tangential trace on a boundary face with vertices `v`: the face functions and the |
