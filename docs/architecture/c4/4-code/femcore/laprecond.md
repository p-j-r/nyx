# Preconditioners: code

C4 level 4 (code) of the component **Preconditioners** (`femcore.laprecond`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> Smoothed-aggregation multigrid (scalar and elastic), auxiliary-space Maxwell and divergence solvers, Schwarz, field-split and incomplete factorisations.

**Decisions:** [ADR 33](../../../adr/0033-one-linear-algebra-driver-facade.md), [ADR 8](../../../adr/0008-ported-solver-library-staged-not-wired.md)  
**Size:** 8 files, 4019 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/linalg/ads.rs` | 196 | **Auxiliary-space divergence solver** (ADS, Hiptmair–Xu 2007; Kolev–Vassilevski 2012 for the practical variant) for lowest-order `H(div)` (Raviart–Thomas, face element) matrices `A = Dᵀ W D + β M_f` (port item P1, the `H(div)` c... |
| `crates/femcore/src/linalg/amg.rs` | 327 | Smoothed-aggregation algebraic multigrid for real and **complex symmetric** matrices (port items "AMG for complex symmetric matrices" and the multigrid half of the shifted-Laplacian preconditioner, P2). |
| `crates/femcore/src/linalg/amg_elastic/mod.rs` | 1377 | Smoothed-aggregation algebraic multigrid and a preconditioned conjugate gradient solver for large symmetric positive definite systems (3D elasticity and diffusion with millions of unknowns, where a sparse direct factorisation needs too m... |
| `crates/femcore/src/linalg/amg_elastic/single.rs` | 68 | The single-precision matrix copy used by the multigrid cycles, and a dot product. |
| `crates/femcore/src/linalg/ams.rs` | 972 | **Auxiliary-space Maxwell solver** (AMS, Hiptmair–Xu) for lowest-order `H(curl)` matrices (port item P1), for the curl-curl-plus-mass operator `A = K + βM` (also complex symmetric). |
| `crates/femcore/src/linalg/fieldsplit.rs` | 280 | Field-split and Schur-complement preconditioners (port item P5; the role of PETSc `PCFIELDSPLIT`) for coupled problems: saddle points, piezoelectric `(u, φ)`, thermo-mechanical or opto-mechanical couplings. |
| `crates/femcore/src/linalg/ilu.rs` | 445 | Incomplete factorisations as preconditioners (port item P3): **ILU(k)**, **ILUT** and **incomplete Cholesky** (IC(k); Hermitian `LLᴴ` or complex-symmetric `LLᵀ`), for real and complex matrices. Sequential and deterministic; natural o... |
| `crates/femcore/src/linalg/schwarz.rs` | 354 | Overlapping **additive and restricted additive Schwarz** preconditioners (port item P4; the role of PETSc `PCASM`) with exact local solves by sparse factorisation and an optional Nicolaides coarse space (one constant per subdomain), for ... |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Matrices and vectors](./lacore.md) | linalg::csr, linalg::csr::assembly, linalg::faer_glue, linalg::linop, linalg::par, linalg::scalar, linalg::vecops | C64, Csr, FaerScalar, Preconditioner, RCsr, RawCols, Scalar, chunk_map, for_chunks_mut |
| [Direct solvers](./ladirect.md) | linalg::direct | Definiteness, DirectSolve, ShiftedFactorization, ShiftedSolver, SolverError |
| [Linear algebra driver](./ladriver.md) | linalg::driver | DirectSolver |
| [Krylov solvers](./laksp.md) | linalg::ksp, linalg::ksp::gmres | Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances |
| [Runtime services](./runtime.md) | par, progress | check, lap_since, note, stage, threads |

## Used by

| Component | Modules | Items |
|---|---|---|
| [dfs: Electrostatics kernel](../dfs/electro.md) | linalg::amg_elastic | AmgOptions |
| [ems: 3D cavity kernel](../ems/cavity.md) | linalg::amg_elastic | spmm |
| [ems: Driven frequency-domain and PROM](../ems/driven.md) | linalg::amg_elastic, linalg::ams | AmsFgmresOptions, AmsFgmresSolver, AmsP2Topology, assemble_shifted_curl_curl, spmm |
| [ems: Thin-film IBC and Josephson ROM](../ems/ibc.md) | linalg::amg_elastic | spmm |
| [ems: Native mixed-cell EM](../ems/mixed.md) | linalg::ams | Ams, AmsOptions |
| [Finite-element basics](./fem.md) | linalg::amg_elastic, linalg::ams | Amg, AmgOptions, AmsP2Topology, block_pcg, pcg |
| [Matrices and vectors](./lacore.md) | linalg::ams | Ams, AmsFgmresOptions, AmsFgmresSolver, AmsFgmresStats, AmsOptions, AmsP2, AmsP2Options, AmsP2Topology +3 |
| [Linear algebra driver](./ladriver.md) | linalg::ams | AmsP2, AmsP2Options, AmsP2Topology |
| [sms: Piezoelectric kernel](../sms/piezo.md) | linalg::amg_elastic | AmgOptions |
| [sms: Solid kernel](../sms/solid.md) | linalg::amg_elastic | Amg, AmgOptions, pcg |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace linalg_ads_rs {
    class AdsOptions {
        <<struct>>
        +amg AmgOptions
        +ams AmsOptions
        +smoother_sweeps usize
        +edge_shift f64
    }
    class Ads {
        <<struct>>
        +new(a, curl, gradient, edge_vectors, edge_mass, pi, opts) Result
    }
  }
  namespace linalg_amg_rs {
    class AmgOptions {
        <<struct>>
        +theta f64
        +sweeps usize
        +chebyshev usize
        +max_coarse usize
        +max_levels usize
        +gamma usize
        +smooth_p usize
        +block bool
        +more_5_fields
        +with_env() Self
    }
    class AmgStats {
        <<struct>>
        +level_sizes Vec~usize~
        +level_nnz Vec~usize~
        +operator_complexity f64
    }
    class Amg {
        <<struct>>
        +new(a, opts) Result
        +stats() AmgStats
        +new(a, nullspace, bs, opts) Result
        +new_p2(a, nullspace, bs, opts, lin) Result
        +n_levels() usize
        +sizes() Vec~usize~
        +line_smooth_active() bool
        +n_line_columns() usize
        +more_3_methods()
    }
    class fns_linalg_amg_rs {
        <<module linalg/amg.rs>>
        ~gauss_seidel(a, inv_diag, b, x, forward)
        +shifted_laplacian_operator(k, mass, k2, beta) Csr~C64~
    }
  }
  namespace linalg_amg_elastic_mod_rs {
    class PcgResult {
        <<struct>>
        +x Vec~f64~
        +iterations usize
        +relative_residual f64
    }
    class BlockPcgResult {
        <<struct>>
        +x Vec~f64~
        +iterations usize
        +relative_residual f64
        +residuals Vec~f64~
    }
    class fns_linalg_amg_elastic_mod_rs {
        <<module linalg/amg_elastic/mod.rs>>
        +spmm(a, b) Csr
        +block_pcg(a, b, k, amg, tol, max_iter) BlockPcgResult
        +pcg(a, b, amg, tol, max_iter) PcgResult
    }
  }
  namespace linalg_amg_elastic_single_rs {
    class Csr32 {
        <<struct>>
        +nrows usize
        +indptr Vec~usize~
        +indices Vec~u32~
        +data Vec~f32~
        +from_csr(a) Csr32
        +nnz() usize
        +matvec_into(x, y)
        +matvec(x) Vec~f64~
    }
    class fns_linalg_amg_elastic_single_rs {
        <<module linalg/amg_elastic/single.rs>>
        +dot(a, b) f64
    }
  }
  namespace linalg_ams_rs {
    class AmsOptions {
        <<struct>>
        +amg AmgOptions
        +smoother_sweeps usize
    }
    class Ams {
        <<struct>>
        +new(a, g, edge_vectors, opts) Result
        +nnz_storage() usize
    }
    class AmsP2Topology {
        <<struct>>
        +p1_to_p2 Csr~f64~
        +gradient_p1 Csr~f64~
        +gradient_p2 Option
        +edge_vectors Vec
        +blocks Vec
    }
    class AmsP2Options {
        <<struct>>
        +coarse_ams AmsOptions
        +amg_g2 AmgOptions
        +smoother_sweeps usize
        +coarse_sweeps usize
    }
    class AmsP2 {
        <<struct>>
        +new(a, topo, opts) Result
        +dim() usize
        +nnz_storage() usize
    }
    class AmsFgmresOptions {
        <<struct>>
        +rtol f64
        +atol f64
        +max_it usize
        +restart usize
        +flexible bool
        +ams AmsOptions
        +ams_p2 AmsP2Options
    }
    class AmsFgmresStats {
        <<struct>>
        +total_iterations usize
        +max_relative_residual f64
        +all_converged bool
        +preconditioner_nnz usize
    }
    class AmsFgmresSolver {
        <<struct>>
        +new_p1(a, pc_matrix, g, edge_vectors, opts) Result
        +new_p2(a, pc_matrix, topo, opts) Result
        +update_matrix(a, pc_matrix) Result
        +dim() usize
        +preconditioner_nnz() usize
        +options() AmsFgmresOptions
        +solve_one(rhs) tuple
        +solve_many(rhs_cols) tuple
        +more_2_methods()
    }
    class AmsPrecondKind {
        <<impl block>>
    }
    class fns_linalg_ams_rs {
        <<module linalg/ams.rs>>
        +assemble_shifted_curl_curl(k, c_damp, m, k0) tuple
        +solve_ams_fgmres_many(a, pc_matrix, g, edge_vectors, rhs_cols, opts) Result
        +solve_ams_p2_fgmres_many(a, pc_matrix, topo, rhs_cols, opts) Result
    }
  }
  namespace linalg_fieldsplit_rs {
    class SplitKind {
        <<enumeration>>
        Additive
        Multiplicative
        SymmetricMultiplicative
    }
    class FieldSplit {
        <<struct>>
        +new(a, fields, make, kind) Result
    }
    class SchurFact {
        <<enumeration>>
        Diag
        Lower
        Upper
        Full
    }
    class SchurPre {
        <<enumeration>>
        A11
        SelfP
        User
    }
    class SchurFieldSplit {
        <<struct>>
        +new(a, i0, i1, fact, pre, negate_schur, make_a00, make_s) Result
    }
    class fns_linalg_fieldsplit_rs {
        <<module linalg/fieldsplit.rs>>
        +extract_block(a, rows, cols) Csr~T~
        +direct_maker(definiteness) FnCsrTResultBoxdynPreconditionerTSolverError
        +MakePc type
    }
  }
  namespace linalg_ilu_rs {
    class Ilu {
        <<struct>>
        +perturbed_pivots usize
        +ilu_k(a, k) Self
        +ilut(a, droptol, lfil) Self
        +nnz() usize
        +solve(r, z)
    }
    class Ic {
        <<struct>>
        +shift f64
        +attempts usize
        +ic_k(a, k, hermitian) Self
        +nnz() usize
        +solve(r, z)
    }
  }
  namespace linalg_schwarz_rs {
    class SchwarzKind {
        <<enumeration>>
        Basic
        Restricted
    }
    class SchwarzOptions {
        <<struct>>
        +overlap usize
        +kind SchwarzKind
        +coarse bool
        +workers usize
        +definiteness Definiteness
    }
    class Schwarz {
        <<struct>>
        +new(a, part, opts) Result
        +subdomains() usize
        +sizes() Vec~usize~
    }
    class fns_linalg_schwarz_rs {
        <<module linalg/schwarz.rs>>
        +contiguous_partition(n, parts) Vec~usize~
        +graph_partition(a, parts) Vec~usize~
    }
  }
  AdsOptions --> AmgOptions
  AdsOptions --> AmsOptions
  AmsOptions --> AmgOptions
  AmsP2Options --> AmsOptions
  AmsP2Options --> AmgOptions
  AmsFgmresOptions --> AmsOptions
  AmsFgmresOptions --> AmsP2Options
  SchwarzOptions --> SchwarzKind
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `AdsOptions` | `linalg/ads.rs:35` | fields: `amg`, `ams`, `smoother_sweeps`, `edge_shift` |
| struct | `Ads` | `linalg/ads.rs:63` |  |
| fn | `Ads::new` | `linalg/ads.rs:82` | `pub fn new( a: &Csr<T>, curl: &Csr<f64>, gradient: &Csr<f64>, edge_vectors: &[[f64; 3]], edge_mass: &Csr<f64>, pi: &[Csr<f64>; 3], opts: &AdsOptions, ) -> Result<Self, SolverError>` —  * `a`: the face matrix; |
| struct | `AmgOptions` | `linalg/amg.rs:29` | fields: `theta`, `max_coarse`, `max_levels`, `smooth_prolongator`, `sweeps`, `hermitian` |
| struct | `AmgStats` | `linalg/amg.rs:66` | fields: `level_sizes`, `level_nnz`, `operator_complexity` —  Summary of the hierarchy, for logs and tests. |
| struct | `Amg` | `linalg/amg.rs:73` |  |
| fn | `gauss_seidel` | `linalg/amg.rs:172` | `pub(crate) fn gauss_seidel<T: Scalar>(a: &Csr<T>, inv_diag: &[T], b: &[T], x: &mut [T], forward: bool)` —  One Gauss–Seidel sweep (`forward` or backward) for `A x = b`, in place on `x`. |
| fn | `shifted_laplacian_operator` | `linalg/amg.rs:314` | `pub fn shifted_laplacian_operator<T: Scalar>( k: &Csr<T>, mass: Option<&Csr<T>>, k2: f64, beta: crate::linalg::scalar::C64, ) -> Csr<crate::linalg::scalar::C64>` —  Builds the shifted-Laplacian operator `K − (β₁ + iβ₂)·k²·M` (Erlangga et al.): the matrix that |
| fn | `Amg::new` | `linalg/amg.rs:81` | `pub fn new(a: &Csr<T>, opts: AmgOptions) -> Result<Self, crate::linalg::direct::SolverError>` |
| fn | `Amg::stats` | `linalg/amg.rs:125` | `pub fn stats(&self) -> AmgStats` |
| struct | `AmgOptions` | `linalg/amg_elastic/mod.rs:19` | fields: `theta`, `sweeps`, `chebyshev`, `max_coarse`, `max_levels`, `gamma`, `smooth_p`, `block`, `no_p1_coarse`, `smooth_lin`, `line_smooth`, `p2_columns`, ... |
| struct | `Amg` | `linalg/amg_elastic/mod.rs:340` |  |
| fn | `spmm` | `linalg/amg_elastic/mod.rs:352` | `pub fn spmm(a: &Csr, b: &Csr) -> Csr` —  Sparse matrix product, on all cores by blocks of rows of `a` for large matrices. |
| struct | `PcgResult` | `linalg/amg_elastic/mod.rs:1166` | fields: `x`, `iterations`, `relative_residual` |
| struct | `BlockPcgResult` | `linalg/amg_elastic/mod.rs:1172` | fields: `x`, `iterations`, `relative_residual`, `residuals` |
| fn | `block_pcg` | `linalg/amg_elastic/mod.rs:1211` | `pub fn block_pcg(a: &Csr, b: &[f64], k: usize, amg: &Amg, tol: f64, max_iter: usize) -> BlockPcgResult` —  Simultaneous multi-RHS Block-CG preconditioned by one AMG V-cycle per iteration across all `k` |
| fn | `pcg` | `linalg/amg_elastic/mod.rs:1327` | `pub fn pcg(a: &Csr, b: &[f64], amg: &Amg, tol: f64, max_iter: usize) -> PcgResult` —  Conjugate gradients preconditioned by one AMG V-cycle. |
| fn | `AmgOptions::with_env` | `linalg/amg_elastic/mod.rs:70` | `pub fn with_env(mut self) -> Self` —  Overrides from the environment (NYX_AMG_THETA, NYX_AMG_SWEEPS), for experiments. |
| fn | `Amg::new` | `linalg/amg_elastic/mod.rs:788` | `pub fn new(a: &Csr, nullspace: Vec<Vec<f64>>, bs: usize, opts: &AmgOptions) -> Result<Amg, String>` —  `nullspace`: the near-nullspace vectors in dof order (node-major, block size `bs`). |
| fn | `Amg::new_p2` | `linalg/amg_elastic/mod.rs:796` | `pub fn new_p2(a: &Csr, nullspace: Vec<Vec<f64>>, bs: usize, opts: &AmgOptions, lin: &Csr) -> Result<Amg, String>` —  Like `new` for a quadratic (P2) discretisation: the first coarse space is the linear (P1) |
| fn | `Amg::n_levels` | `linalg/amg_elastic/mod.rs:921` | `pub fn n_levels(&self) -> usize` |
| fn | `Amg::sizes` | `linalg/amg_elastic/mod.rs:926` | `pub fn sizes(&self) -> Vec<usize>` —  Sizes of the matrices on each level, finest first. |
| fn | `Amg::line_smooth_active` | `linalg/amg_elastic/mod.rs:933` | `pub fn line_smooth_active(&self) -> bool` —  Whether vertical-line Thomas smoothing is active on at least one level of the hierarchy. |
| fn | `Amg::n_line_columns` | `linalg/amg_elastic/mod.rs:938` | `pub fn n_line_columns(&self) -> usize` —  Number of vertical columns smoothed on the primary P1/P2 levels. |
| fn | `Amg::apply` | `linalg/amg_elastic/mod.rs:1147` | `pub fn apply(&self, b: &[f64]) -> Vec<f64>` —  One V-cycle applied to `b` from a zero guess. |
| fn | `Amg::apply_block_into` | `linalg/amg_elastic/mod.rs:1154` | `pub fn apply_block_into(&self, b: &[f64], x: &mut [f64], k: usize)` —  One V-cycle applied simultaneously to `k` right-hand sides stored row-major (`[n x k]`). |
| fn | `Amg::apply_block` | `linalg/amg_elastic/mod.rs:1159` | `pub fn apply_block(&self, b: &[f64], k: usize) -> Vec<f64>` —  One V-cycle applied simultaneously to `k` right-hand sides stored row-major (`[n x k]`). |
| struct | `Csr32` | `linalg/amg_elastic/single.rs:18` | fields: `nrows`, `indptr`, `indices`, `data` —  Compressed rows with 32-bit column indices and single-precision values, for matrix-vector |
| fn | `dot` | `linalg/amg_elastic/single.rs:65` | `pub fn dot(a: &[f64], b: &[f64]) -> f64` —  Dot product. |
| fn | `Csr32::from_csr` | `linalg/amg_elastic/single.rs:27` | `pub fn from_csr(a: &RCsr) -> Csr32` |
| fn | `Csr32::nnz` | `linalg/amg_elastic/single.rs:30` | `pub fn nnz(&self) -> usize` |
| fn | `Csr32::matvec_into` | `linalg/amg_elastic/single.rs:41` | `pub fn matvec_into(&self, x: &[f64], y: &mut [f64])` |
| fn | `Csr32::matvec` | `linalg/amg_elastic/single.rs:55` | `pub fn matvec(&self, x: &[f64]) -> Vec<f64>` |
| struct | `AmsOptions` | `linalg/ams.rs:37` | fields: `amg`, `smoother_sweeps` |
| struct | `Ams` | `linalg/ams.rs:58` |  |
| struct | `AmsP2Topology` | `linalg/ams.rs:189` | fields: `p1_to_p2`, `gradient_p1`, `gradient_p2`, `edge_vectors`, `blocks` —  Topological operators connecting a second-order ($P_2$) Nédélec `H(curl)` space to its |
| struct | `AmsP2Options` | `linalg/ams.rs:210` | fields: `coarse_ams`, `amg_g2`, `smoother_sweeps`, `coarse_sweeps` |
| struct | `AmsP2` | `linalg/ams.rs:415` |  **Auxiliary-space Maxwell solver for second-order ($P_2$) Nédélec `H(curl)` spaces** (p-multigrid |
| struct | `AmsFgmresOptions` | `linalg/ams.rs:567` | fields: `rtol`, `atol`, `max_it`, `restart`, `flexible`, `ams`, `ams_p2` —  Options for [`AmsFgmresSolver`] (deterministic right-preconditioned FGMRES/GMRES driven by |
| struct | `AmsFgmresStats` | `linalg/ams.rs:600` | fields: `total_iterations`, `max_relative_residual`, `all_converged`, `preconditioner_nnz` —  Summary statistics of one [`AmsFgmresSolver::solve_one`] or [`AmsFgmresSolver::solve_many`] call. |
| fn | `assemble_shifted_curl_curl` | `linalg/ams.rs:647` | `pub fn assemble_shifted_curl_curl( k: &Csr<C64>, c_damp: Option<&Csr<C64>>, m: &Csr<C64>, k0: f64, ) -> (Csr<C64>, Csr<C64>)` —  Assembles the complex-symmetric driven curl-curl operator and its shifted-Laplacian auxiliary |
| struct | `AmsFgmresSolver` | `linalg/ams.rs:668` |  Deterministic preconditioned FGMRES/GMRES solver driven by [`Ams`] ($P_1$) or [`AmsP2`] ($P_2$) |
| fn | `solve_ams_fgmres_many` | `linalg/ams.rs:853` | `pub fn solve_ams_fgmres_many<T: FaerScalar>( a: &Csr<T>, pc_matrix: Option<&Csr<T>>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], rhs_cols: &[Vec<T>], opts: &AmsFgmresOptions, ) -> Result<(Vec<Vec<T>>, AmsFgmresStats), SolverError>` —  Convenience helper solving `A X = B` for multiple RHS columns using $P_1$ [`Ams`] + deterministic FGMRES. |
| fn | `solve_ams_p2_fgmres_many` | `linalg/ams.rs:866` | `pub fn solve_ams_p2_fgmres_many<T: FaerScalar>( a: &Csr<T>, pc_matrix: Option<&Csr<T>>, topo: &AmsP2Topology, rhs_cols: &[Vec<T>], opts: &AmsFgmresOptions, ) -> Result<(Vec<Vec<T>>, AmsFgmresStats), SolverError>` —  Convenience helper solving `A X = B` for multiple RHS columns using $P_2$ [`AmsP2`] + deterministic FGMRES. |
| fn | `Ams::new` | `linalg/ams.rs:78` | `pub fn new(a: &Csr<T>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: &AmsOptions) -> Result<Self, SolverError>` —  `a`: the edge matrix (SPD, or complex symmetric); `g`: discrete gradient (edges × interior |
| fn | `Ams::nnz_storage` | `linalg/ams.rs:137` | `pub fn nnz_storage(&self) -> usize` —  Approximate total sparse matrix non-zeros stored across this AMS hierarchy (for memory comparisons). |
| fn | `AmsP2::new` | `linalg/ams.rs:431` | `pub fn new(a: &Csr<T>, topo: &AmsP2Topology, opts: &AmsP2Options) -> Result<Self, SolverError>` —  Builds the $P_2$ AMS preconditioner for `a` (`ndofs_p2 × ndofs_p2`, SPD or complex symmetric |
| fn | `AmsP2::dim` | `linalg/ams.rs:491` | `pub fn dim(&self) -> usize` |
| fn | `AmsP2::nnz_storage` | `linalg/ams.rs:497` | `pub fn nnz_storage(&self) -> usize` —  Approximate total sparse matrix non-zeros stored across this $P_2$ AMS preconditioner |
| fn | `AmsFgmresSolver::new_p1` | `linalg/ams.rs:680` | `pub fn new_p1( a: &Csr<T>, pc_matrix: Option<&Csr<T>>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: AmsFgmresOptions, ) -> Result<Self, SolverError>` —  Builds a $P_1$ `Ams`-preconditioned FGMRES/GMRES solver for `a`. |
| fn | `AmsFgmresSolver::new_p2` | `linalg/ams.rs:708` | `pub fn new_p2( a: &Csr<T>, pc_matrix: Option<&Csr<T>>, topo: &AmsP2Topology, opts: AmsFgmresOptions, ) -> Result<Self, SolverError>` —  Builds a $P_2$ `AmsP2`-preconditioned FGMRES/GMRES solver for `a`. |
| fn | `AmsFgmresSolver::update_matrix` | `linalg/ams.rs:733` | `pub fn update_matrix(&mut self, a: &Csr<T>, pc_matrix: Option<&Csr<T>>) -> Result<(), SolverError>` —  Updates the system matrix `a` (and optional auxiliary preconditioner matrix `pc_matrix`) at a |
| fn | `AmsFgmresSolver::dim` | `linalg/ams.rs:752` | `pub fn dim(&self) -> usize` |
| fn | `AmsFgmresSolver::preconditioner_nnz` | `linalg/ams.rs:756` | `pub fn preconditioner_nnz(&self) -> usize` |
| fn | `AmsFgmresSolver::options` | `linalg/ams.rs:760` | `pub fn options(&self) -> &AmsFgmresOptions` |
| fn | `AmsFgmresSolver::solve_one` | `linalg/ams.rs:765` | `pub fn solve_one(&self, rhs: &[T]) -> (Vec<T>, AmsFgmresStats)` —  Solves `A x = rhs` for a single right-hand side using deterministic FGMRES/GMRES. |
| fn | `AmsFgmresSolver::solve_many` | `linalg/ams.rs:797` | `pub fn solve_many(&self, rhs_cols: &[Vec<T>]) -> (Vec<Vec<T>>, AmsFgmresStats)` —  Solves `A X = B` across multiple right-hand-side columns using deterministic FGMRES/GMRES. |
| fn | `AmsFgmresSolver::new_shifted_complex_p1` | `linalg/ams.rs:824` | `pub fn new_shifted_complex_p1( k: &Csr<C64>, c_damp: Option<&Csr<C64>>, m: &Csr<C64>, k0: f64, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: AmsFgmresOptions, ) -> Result<Self, SolverError>` —  Builds a $P_1$ `Ams<C64>` FGMRES/GMRES solver for `(K + i k_0 C - k_0^2 M) X = B` |
| fn | `AmsFgmresSolver::new_shifted_complex_p2` | `linalg/ams.rs:839` | `pub fn new_shifted_complex_p2( k: &Csr<C64>, c_damp: Option<&Csr<C64>>, m: &Csr<C64>, k0: f64, topo: &AmsP2Topology, opts: AmsFgmresOptions, ) -> Result<Self, SolverError>` —  Builds a $P_2$ `AmsP2<C64>` FGMRES/GMRES solver for `(K + i k_0 C - k_0^2 M) X = B` |
| type | `MakePc` | `linalg/fieldsplit.rs:26` | `<'a, T>` —  Builds a sub-preconditioner from a diagonal block. |
| fn | `extract_block` | `linalg/fieldsplit.rs:29` | `pub fn extract_block<T: Scalar>(a: &Csr<T>, rows: &[usize], cols: &[usize]) -> Csr<T>` —  Extracts `A[rows, cols]` (rows/cols are global indices; the result is renumbered locally). |
| enum | `SplitKind` | `linalg/fieldsplit.rs:48` | variants: `Additive`, `Multiplicative`, `SymmetricMultiplicative` |
| struct | `FieldSplit` | `linalg/fieldsplit.rs:59` |  |
| enum | `SchurFact` | `linalg/fieldsplit.rs:157` | variants: `Diag`, `Lower`, `Upper`, `Full` |
| enum | `SchurPre` | `linalg/fieldsplit.rs:164` | variants: `A11`, `SelfP`, `User` |
| struct | `SchurFieldSplit` | `linalg/fieldsplit.rs:173` |  |
| fn | `direct_maker` | `linalg/fieldsplit.rs:275` | `pub fn direct_maker<T: crate::linalg::faer_glue::FaerScalar>( definiteness: crate::linalg::direct::Definiteness, ) -> impl Fn(&Csr<T>) -> Result<Box<dyn Preconditioner<T>>, SolverError>` —  Convenience builder for the common case: an exact (sparse direct) solver for a block. |
| fn | `FieldSplit::new` | `linalg/fieldsplit.rs:68` | `pub fn new(a: &Csr<T>, fields: Vec<Vec<usize>>, make: MakePc<'_, T>, kind: SplitKind) -> Result<Self, SolverError>` |
| fn | `SchurFieldSplit::new` | `linalg/fieldsplit.rs:188` | `pub fn new( a: &Csr<T>, i0: Vec<usize>, i1: Vec<usize>, fact: SchurFact, pre: SchurPre<T>, negate_schur: bool, make_a00: MakePc<'_, T>, make_s: MakePc<'_, T>, ) -> Result<Self, SolverError>` —  `i0`, `i1`: the two fields (disjoint, covering all unknowns). `make_a00` / `make_s` build the solvers for |
| struct | `Ilu` | `linalg/ilu.rs:27` | fields: `perturbed_pivots` —  Row-wise sparse triangular storage used by both factorisations: strictly-lower part, diagonal, strictly-upper part. |
| struct | `Ic` | `linalg/ilu.rs:307` | fields: `shift`, `attempts` —  Incomplete Cholesky `A ≈ L Lᴴ` (Hermitian) or `L Lᵀ` (complex symmetric / real symmetric). |
| fn | `Ilu::ilu_k` | `linalg/ilu.rs:103` | `pub fn ilu_k(a: &Csr<T>, k: u32) -> Self` —  ILU(`k`): `k = 0` is ILU(0). |
| fn | `Ilu::ilut` | `linalg/ilu.rs:145` | `pub fn ilut(a: &Csr<T>, droptol: f64, lfil: usize) -> Self` —  ILUT with drop tolerance `droptol` (relative to the 2-norm of the row of `A`) and at most `lfil` kept |
| fn | `Ilu::nnz` | `linalg/ilu.rs:269` | `pub fn nnz(&self) -> usize` |
| fn | `Ilu::solve` | `linalg/ilu.rs:274` | `pub fn solve(&self, r: &[T], z: &mut [T])` —  `z ← (LU)⁻¹ r`. |
| fn | `Ic::ic_k` | `linalg/ilu.rs:321` | `pub fn ic_k(a: &Csr<T>, k: u32, hermitian: bool) -> Self` —  IC(`k`) of the lower triangle of `a` (assumed Hermitian if `hermitian`, else symmetric). |
| fn | `Ic::nnz` | `linalg/ilu.rs:412` | `pub fn nnz(&self) -> usize` |
| fn | `Ic::solve` | `linalg/ilu.rs:417` | `pub fn solve(&self, r: &[T], z: &mut [T])` —  `z ← (L Lᴴ)⁻¹ r` (or `(L Lᵀ)⁻¹ r`). |
| enum | `SchwarzKind` | `linalg/schwarz.rs:33` | variants: `Basic`, `Restricted` |
| struct | `SchwarzOptions` | `linalg/schwarz.rs:39` | fields: `overlap`, `kind`, `coarse`, `workers`, `definiteness` |
| fn | `contiguous_partition` | `linalg/schwarz.rs:63` | `pub fn contiguous_partition(n: usize, parts: usize) -> Vec<usize>` —  `part[i] = i·parts/n`: contiguous index ranges of (almost) equal size. |
| fn | `graph_partition` | `linalg/schwarz.rs:69` | `pub fn graph_partition<T: Scalar>(a: &Csr<T>, parts: usize) -> Vec<usize>` —  Breadth-first region growing on the (symmetrised) graph of `a`: seeds are spread by repeated |
| struct | `Schwarz` | `linalg/schwarz.rs:168` |  |
| fn | `Schwarz::new` | `linalg/schwarz.rs:180` | `pub fn new(a: &Csr<T>, part: &[usize], opts: &SchwarzOptions) -> Result<Self, SolverError>` —  `part[i]` assigns dof `i` to a subdomain (`0..nparts`). |
| fn | `Schwarz::subdomains` | `linalg/schwarz.rs:289` | `pub fn subdomains(&self) -> usize` |
| fn | `Schwarz::sizes` | `linalg/schwarz.rs:294` | `pub fn sizes(&self) -> Vec<usize>` —  Sizes of the (overlapping) subdomains. |

Trait implementations in this component: `AdsOptions: Default`, `Ads: Preconditioner`, `AmgOptions: Default`, `Amg: Preconditioner`, `AmgOptions: Default`, `AmsOptions: Default`, `Ams: Preconditioner`, `AmsP2Options: Default`, `AmsP2: Preconditioner`, `AmsFgmresOptions: Default`, `AmsPrecondKind: Preconditioner`, `FieldSplit: Preconditioner`, `SchurFieldSplit: Preconditioner`, `Ilu: Preconditioner`, `Ic: Preconditioner`, `SchwarzOptions: Default`, `Schwarz: Preconditioner`
