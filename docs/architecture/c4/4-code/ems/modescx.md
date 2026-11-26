# Complex-mode kernel: code

C4 level 4 (code) of the component **Complex-mode kernel** (`ems.modescx`) in the container **nyx-ems** (crate `nyx_ems`). [Component view of the container](../../3-components/ems.md) · [All code pages](../README.md)

> Complex-valued mode analysis: lossy media and perfectly matched layers (`+4` Krylov + `layer_fraction` filter, R32 check).

**Decisions:** [ADR 45](../../../adr/0045-ems-moves-into-nyx-ems.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md)  
**Size:** 1 file, 1559 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-ems/src/modes_complex.rs` | 1559 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | edge, edge_mixed, fem, fem_mixed, fields, linalg::csr, linalg::csr::assembly, linalg::driver +5 | C, CCoo, CCsr, CellFields, EpsC, IsoparametricCell, LagrangeSpace, MixedDiscretisation, MixedLagrangeDofMap, MuC +23 |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | crate root, result | Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult, cell_fields_of |
| [Thin-film IBC and Josephson ROM](./ibc.md) | ibc | EPS0, MU0, add_ibc_2d_waveguide, add_ibc_2d_waveguide_mixed, pec_edges_and_nodes_2d, pec_edges_and_nodes_2d_mixed |
| [Native mixed-cell EM](./mixed.md) | mixed_ops | C3, EmPmlTable, MixedCellMaterials, complex_eps_c3, iso_mu_inv_c3 |
| [Mode kernel](./modes.md) | crate root | C0 |
| [Optical mode sets](./optical.md) | optical | mode_frame |
| [Typed problem](./problem.md) | problem | Problem, Study, SurfaceImpedanceBc, Symmetry, pec_edges_and_nodes_2d_with_symmetry |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Physics adapter](./adapter.md) | modes_complex | ModeAssemblyC, ModesC, assemble_modes_complex, assemble_modes_complex_mixed, assemble_modes_complex_with_bcs, mode_fields_complex, mode_fields_complex_cells, solve_mode_pairs_complex +3 |
| [Mode kernel](./modes.md) | modes_complex | compute_waveguide_2d_metrics, compute_waveguide_2d_metrics_mixed, solve_waveguide_complex, solve_waveguide_complex_mixed |
| [Optical mode sets](./optical.md) | modes_complex | ModesC, solve_modes_complex |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class ModesC {
      <<struct>>
      +neff Vec~C~
      +et Vec
      +phi Vec
      +converged bool
  }
  class ModeAssemblyC {
      <<struct>>
      +a CCsr
      +b CCsr
      +a1 Option~CCsr~
      +free Vec~usize~
      +n_tot usize
      +net usize
      +k0 f64
  }
  class WaveguideModeMetrics {
      <<struct>>
      +vg f64
      +comp_shares slice
      +domain_shares Vec
      +layer_fraction f64
      +push_scalars(mode_num, scalars)
  }
  class fns_modes_complex_rs {
      <<module modes_complex.rs>>
      +solve_modes_complex(mesh, freq, eps, mu, n_target, k, order) Result
      +assemble_modes_complex(mesh, freq, eps, mu, order) Result
      +assemble_modes_complex_with_bcs(mesh, freq, eps, mu, order, pec_labels, ibc_bcs, pmc_labels) Result
      ~assemble_modes_complex_with_bcs_and_symmetry(mesh, freq, eps, mu, order, pec_labels, ibc_bcs, pmc_labels, symmetry) Result
      +solve_mode_pairs_complex(asm, n_target, k) Result
      +mode_fields_complex(mesh, modes, j, freq, order) Vec~Field~
      +mode_fields_complex_cells(mesh, modes, j, freq, order) tuple
      +solve_waveguide_complex(mesh, problem) Result
      ~compute_waveguide_2d_metrics(mesh, problem, freq, neff, et, phi) Result
      +solve_waveguide_complex_mixed(mesh, mixed, problem) Result
      +more_2_functions()
  }
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `ModesC` | `modes_complex.rs:59` | fields: `neff`, `et`, `phi`, `converged` |
| fn | `solve_modes_complex` | `modes_complex.rs:67` | `pub fn solve_modes_complex( mesh: &Discretisation, freq: f64, eps: &[EpsC], mu: &[MuC], n_target: C, k: usize, order: usize, ) -> Result<ModesC, String>` —  Modes with complex effective index nearest `n_target` (complex). |
| struct | `ModeAssemblyC` | `modes_complex.rs:82` | fields: `a`, `b`, `a1`, `free`, `n_tot`, `net`, `k0` —  Stage 1 of the complex mode solver: `A x = lambda B x` (`lambda = -beta^2`) on the free dofs, plus the part of `A` |
| fn | `assemble_modes_complex` | `modes_complex.rs:93` | `pub fn assemble_modes_complex( mesh: &Discretisation, freq: f64, eps: &[EpsC], mu: &[MuC], order: usize, ) -> Result<ModeAssemblyC, String>` |
| fn | `assemble_modes_complex_with_bcs` | `modes_complex.rs:103` | `pub fn assemble_modes_complex_with_bcs( mesh: &Discretisation, freq: f64, eps: &[EpsC], mu: &[MuC], order: usize, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], pmc_labels: &[i32], ) -> Result<ModeAssemblyC, String>` |
| fn | `assemble_modes_complex_with_bcs_and_symmetry` | `modes_complex.rs:126` | `pub(crate) fn assemble_modes_complex_with_bcs_and_symmetry( mesh: &Discretisation, freq: f64, eps: &[EpsC], mu: &[MuC], order: usize, pec_labels: &[i32], ibc_bcs: &[SurfaceImpedanceBc], pmc_labels: &[i32], symmetry: Symmetry, ) -> Result<ModeAssemblyC, String>` |
| fn | `solve_mode_pairs_complex` | `modes_complex.rs:298` | `pub fn solve_mode_pairs_complex(asm: &ModeAssemblyC, n_target: C, k: usize) -> Result<ModesC, String>` —  Stage 2: the `k` eigenpairs nearest the complex effective index `n_target`, as effective indices, edge coefficients |
| fn | `mode_fields_complex` | `modes_complex.rs:422` | `pub fn mode_fields_complex(mesh: &Discretisation, modes: &ModesC, j: usize, freq: f64, order: usize) -> Vec<Field>` —  Nodal magnitude fields of a complex mode, scaled to max \|E\| = 1, under wafer names: `\|E\|, \|Exz\|, \|Ey\|, Re Ex, |
| fn | `mode_fields_complex_cells` | `modes_complex.rs:428` | `pub fn mode_fields_complex_cells( mesh: &Discretisation, modes: &ModesC, j: usize, freq: f64, order: usize, ) -> (Vec<Field>, CellFields)` —  Stage 3: `mode_fields_complex` and, next to them, the complex components `Ex`, `Ez`, `Ey` at the local dofs of every |
| fn | `solve_waveguide_complex` | `modes_complex.rs:508` | `pub fn solve_waveguide_complex(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String>` —  The modes of a `Modes` study in the complex solver: absorbing layers, lossy or complex materials, a crystal tilted against the guide. |
| struct | `WaveguideModeMetrics` | `modes_complex.rs:639` | fields: `vg`, `comp_shares`, `domain_shares`, `layer_fraction` —  Poynting group velocity / group index and non-PML electric energy shares for a 2D waveguide mode. |
| fn | `compute_waveguide_2d_metrics` | `modes_complex.rs:664` | `pub(crate) fn compute_waveguide_2d_metrics( mesh: &Discretisation, problem: &Problem, freq: f64, neff: &[C], et: &[Vec<C>], phi: &[Vec<C>], ) -> Result<Vec<WaveguideModeMetrics>, String>` —  Evaluates `P_y`, dispersive `U_EM`, `v_g`, `n_g`, component shares (`Ex`, `Ey`, `Ez`), domain shares, |
| fn | `solve_waveguide_complex_mixed` | `modes_complex.rs:1060` | `pub fn solve_waveguide_complex_mixed( mesh: &Mesh, mixed: &MixedDiscretisation, problem: &Problem, ) -> Result<SolveResult, String>` —  Solves a 2D complex `Modes` study natively on `Tri + Quad` (`MixedDiscretisation`). |
| fn | `assemble_modes_complex_mixed` | `modes_complex.rs:1170` | `pub fn assemble_modes_complex_mixed( mesh: &MixedDiscretisation, problem: &Problem, freq: f64, ) -> Result<(ModeAssemblyC, NedMixed, MixedLagrangeDofMap), String>` |
| fn | `compute_waveguide_2d_metrics_mixed` | `modes_complex.rs:1412` | `pub(crate) fn compute_waveguide_2d_metrics_mixed( mesh: &MixedDiscretisation, ned: &NedMixed, phi_map: &MixedLagrangeDofMap, problem: &Problem, freq: f64, neff: &[C], et: &[Vec<C>], phi: &[Vec<C>], ) -> Result<Vec<WaveguideModeMetrics>, String>` |
| fn | `WaveguideModeMetrics::push_scalars` | `modes_complex.rs:647` | `pub fn push_scalars(&self, mode_num: usize, scalars: &mut Vec<(String, Scalar)>)` |
