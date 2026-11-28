//! **Auxiliary-space Maxwell solver** (AMS, Hiptmair–Xu) for lowest-order `H(curl)` matrices
//! (port item P1), for the curl-curl-plus-mass operator `A = K + βM` (also complex symmetric).
//!
//! The null space of `K` (the gradients) defeats ordinary multigrid. AMS removes it with two
//! auxiliary nodal spaces, using the discrete gradient `G` and the vertex coordinates:
//!
//! * `G`: edges × nodes; `K G = 0`. The correction `G · AMG(Gᵀ A G)⁻¹ · Gᵀ` handles the gradient
//!   (kernel) part.
//! * `Π = [Π_x Π_y Π_z]`: the Nédélec interpolation of nodal vector fields,
//!   `(Π_c)_{e,j} = ½ |G_{e,j}| · (t_e)_c` with `t_e = Σ_j G_{e,j} x_j` the edge vector. The correction
//!   `Σ_c Π_c · AMG(Π_cᵀ A Π_c)⁻¹ · Π_cᵀ` handles the smooth (non-gradient) part. As in the simplest
//!   hypre variant the vector Laplacian is approximated by its three scalar diagonal blocks.
//!
//! The cycle is the symmetric multiplicative one (`01210` in hypre's notation): forward
//! Gauss–Seidel, `G` correction, `Π` correction (the three components additive), `G` correction,
//! backward Gauss–Seidel. The
//! resulting preconditioner is symmetric (`Mᵀ = M`), so it drives CG (SPD `A`) and COCG (complex
//! symmetric `A`) alike.
//!
//! Limits: lowest-order elements only; `G` must be the *interior* gradient consistent with the
//! boundary conditions of `A` (same constrained edges, constrained nodes removed); no
//! rotated/anisotropic coefficient handling; Gauss–Seidel is sequential.
//!
//! Written from the Hiptmair–Xu paper (2007) and Kolev–Vassilevski (2009), not from hypre's code.

use crate::linalg::amg::{
    Amg, AmgF32, AmgOptions, CsrF32, SmootherKind, chebyshev_smooth, chebyshev_smooth_f32,
    estimate_dinv_a_spectral_radius, f32_safe_scale, gauss_seidel, gauss_seidel_f32,
};
use crate::linalg::csr::Csr;
use crate::linalg::device::CpuVec;
use crate::linalg::direct::SolverError;
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::ksp::{Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
use crate::linalg::linop::Preconditioner;
use crate::linalg::par;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{self, CHUNK};
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct AmsOptions {
    /// Options of the AMG hierarchies on the auxiliary nodal operators.
    pub amg: AmgOptions,
    /// Gauss–Seidel sweeps of the edge smoother (each direction, when `smoother == GaussSeidel`).
    pub smoother_sweeps: usize,
    /// Fine-edge smoother selection (default [`SmootherKind::GaussSeidel`]).
    pub smoother: SmootherKind,
}

impl Default for AmsOptions {
    fn default() -> Self {
        Self { amg: AmgOptions::default(), smoother_sweeps: 1, smoother: SmootherKind::GaussSeidel }
    }
}

struct Scratch<T> {
    x: Vec<T>,
    res: Vec<T>,
    tmp_edge: Vec<T>,
    tmp_node: Vec<T>,
    corr_node: Vec<T>,
}

pub struct Ams<T: FaerScalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    lmax: f64,
    smoother: SmootherKind,
    g: Csr<T>,
    gt: Csr<T>,
    amg_g: Amg<T>,
    pi: [Csr<T>; 3],
    pit: [Csr<T>; 3],
    amg_pi: [Amg<T>; 3],
    sweeps: usize,
    scratch: RefCell<Scratch<T>>,
}

impl<T: FaerScalar> Ams<T> {
    /// `a`: the edge matrix (SPD, or complex symmetric); `g`: discrete gradient (edges × interior
    /// nodes); `edge_vectors[e]`: the vector `x_b − x_a` from the first to the second end node of
    /// edge `e` (the orientation `g`'s row uses: `−1` at `a`, `+1` at `b`).
    ///
    /// The interpolation `Π_c` takes `½ (t_e)_c` at every *interior* end node of `e`; an end node
    /// removed by the boundary condition carries the value zero, so it contributes nothing.
    pub fn new(a: &Csr<T>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: &AmsOptions) -> Result<Self, SolverError> {
        let ne = a.nrows();
        assert_eq!(a.ncols(), ne);
        assert_eq!(g.nrows(), ne);
        assert_eq!(edge_vectors.len(), ne);
        let nn = g.ncols();
        let gt_real = g.transpose();
        let gm: Csr<T> = g.map(T::from_f64);
        let gtm: Csr<T> = gt_real.map(T::from_f64);
        let a_g = gtm.matmul(&a.matmul(&gm));
        let amg_g = Amg::new(&a_g, opts.amg.clone())?;
        let mut pis = Vec::new();
        let mut pits = Vec::new();
        let mut amgs = Vec::new();
        for c in 0..3 {
            let mut trip = Vec::with_capacity(g.nnz());
            for (e, t) in edge_vectors.iter().enumerate() {
                for &j in g.row(e).0 {
                    trip.push((e, j as usize, 0.5 * t[c]));
                }
            }
            let pic: Csr<f64> = Csr::from_triplets(ne, nn, &trip).unwrap();
            let pm: Csr<T> = pic.map(T::from_f64);
            let pmt: Csr<T> = pic.transpose().map(T::from_f64);
            let a_pi = pmt.matmul(&a.matmul(&pm));
            amgs.push(Amg::new(&a_pi, opts.amg.clone())?);
            pis.push(pm);
            pits.push(pmt);
        }
        let inv_diag: Vec<T> =
            a.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();
        let lmax = if matches!(opts.smoother, SmootherKind::Chebyshev { .. }) {
            estimate_dinv_a_spectral_radius(a, &inv_diag)
        } else {
            1.0
        };
        let z = |n: usize| vec![T::ZERO; n];
        let to3 = |mut v: Vec<Csr<T>>| -> [Csr<T>; 3] {
            let c = v.pop().unwrap();
            let b = v.pop().unwrap();
            let a = v.pop().unwrap();
            [a, b, c]
        };
        let to3a = |mut v: Vec<Amg<T>>| -> [Amg<T>; 3] {
            let c = v.pop().unwrap();
            let b = v.pop().unwrap();
            let a = v.pop().unwrap();
            [a, b, c]
        };
        Ok(Self {
            a: a.clone(),
            inv_diag,
            lmax,
            smoother: opts.smoother,
            g: gm,
            gt: gtm,
            amg_g,
            pi: to3(pis),
            pit: to3(pits),
            amg_pi: to3a(amgs),
            sweeps: opts.smoother_sweeps,
            scratch: RefCell::new(Scratch { x: z(ne), res: z(ne), tmp_edge: z(ne), tmp_node: z(nn), corr_node: z(nn) }),
        })
    }

    /// Approximate total sparse matrix non-zeros stored across this AMS hierarchy (for memory comparisons).
    pub fn nnz_storage(&self) -> usize {
        let amg_nnz = |amg: &Amg<T>| amg.stats().level_nnz.iter().sum::<usize>();
        let pi_nnz: usize = (0..3)
            .map(|c| self.pi[c].nnz() + self.pit[c].nnz() + amg_nnz(&self.amg_pi[c]))
            .sum();
        self.a.nnz() + self.g.nnz() + self.gt.nnz() + amg_nnz(&self.amg_g) + pi_nnz
    }

    /// `x += Σ_i Q_i · B_i⁻¹ · Q_iᵀ (b − A x)`: the corrections are *additive* (all computed from the
    /// same residual), which for the three `Π_c` is the block-diagonal vector-Laplacian solve and keeps
    /// the whole cycle symmetric.
    fn correct(&self, parts: &[(&Csr<T>, &Csr<T>, &Amg<T>)], b: &[T], sc: &mut Scratch<T>) {
        self.a.spmv(&sc.x, &mut sc.res);
        for (r, bb) in sc.res.iter_mut().zip(b) {
            *r = *bb - *r;
        }
        for &(q, qt, amg) in parts {
            qt.spmv(&sc.res, &mut sc.tmp_node);
            amg.apply(&sc.tmp_node, &mut sc.corr_node);
            q.spmv(&sc.corr_node, &mut sc.tmp_edge);
            for (x, t) in sc.x.iter_mut().zip(&sc.tmp_edge) {
                *x += *t;
            }
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for Ams<T> {
    fn apply(&self, b: &[T], z: &mut [T]) {
        let mut guard = self.scratch.borrow_mut();
        let sc = &mut *guard;
        sc.x.fill(T::ZERO);
        match self.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.sweeps {
                    gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, true);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth(
                    &self.a,
                    &self.inv_diag,
                    self.lmax,
                    degree.max(1),
                    b,
                    &mut sc.x,
                    &mut sc.res,
                    &mut sc.tmp_edge,
                    true,
                );
            }
        }
        let grad = [(&self.g, &self.gt, &self.amg_g)];
        let pis = [0, 1, 2].map(|c| (&self.pi[c], &self.pit[c], &self.amg_pi[c]));
        self.correct(&grad, b, sc);
        self.correct(&pis, b, sc);
        self.correct(&grad, b, sc);
        match self.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.sweeps {
                    gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, false);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth(
                    &self.a,
                    &self.inv_diag,
                    self.lmax,
                    degree.max(1),
                    b,
                    &mut sc.x,
                    &mut sc.res,
                    &mut sc.tmp_edge,
                    false,
                );
            }
        }
        z.copy_from_slice(&sc.x);
    }
}

impl<T: FaerScalar> crate::linalg::device::DevicePreconditioner<CpuVec<T>> for Ams<T> {
    #[inline]
    fn apply(&self, r: &CpuVec<T>, z: &mut CpuVec<T>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}

// ------------------------------------------------------------------------------------------------------------ Mixed-Precision AmsF32

struct ScratchAmsF32 {
    rhs: Vec<f32>,
    x: Vec<f32>,
    res: Vec<f32>,
    tmp_edge: Vec<f32>,
    tmp_node: Vec<f32>,
    corr_node: Vec<f32>,
}

/// Mixed-precision single-precision (`f32`) Hiptmair–Xu AMS preconditioner driving outer `f64`
/// PCG / FGMRES solves (`FR-43a`).
///
/// Both the fine-edge smoother (symmetric Gauss–Seidel or parallel Chebyshev) and all four
/// auxiliary nodal AMG hierarchies (`AmgF32`) execute in `f32`.
pub struct AmsF32 {
    a: CsrF32,
    inv_diag: Vec<f32>,
    lmax: f32,
    smoother: SmootherKind,
    sweeps: usize,
    g: CsrF32,
    gt: CsrF32,
    amg_g: AmgF32,
    pi: [CsrF32; 3],
    pit: [CsrF32; 3],
    amg_pi: [AmgF32; 3],
    scratch: RefCell<ScratchAmsF32>,
}

impl AmsF32 {
    pub fn new(a: &Csr<f64>, g: &Csr<f64>, edge_vectors: &[[f64; 3]], opts: &AmsOptions) -> Result<Self, SolverError> {
        let ne = a.nrows();
        if a.ncols() != ne {
            return Err(SolverError::InvalidArgument("AmsF32 matrix must be square"));
        }
        if g.nrows() != ne || edge_vectors.len() != ne {
            return Err(SolverError::InvalidArgument("AmsF32 gradient and edge_vectors must match matrix dimension"));
        }
        if a.values().iter().any(|v| !v.is_finite()) || g.values().iter().any(|v| !v.is_finite()) {
            return Err(SolverError::InvalidArgument("AmsF32 input contains non-finite entries"));
        }
        let nn = g.ncols();
        let gt_f64 = g.transpose();
        let a_g = gt_f64.matmul(&a.matmul(g));
        let amg_g = AmgF32::new(&a_g, opts.amg.clone())?;

        let mut pis = Vec::with_capacity(3);
        let mut pits = Vec::with_capacity(3);
        let mut amgs = Vec::with_capacity(3);
        for c in 0..3 {
            let mut trip = Vec::with_capacity(g.nnz());
            for (e, t) in edge_vectors.iter().enumerate() {
                for &j in g.row(e).0 {
                    trip.push((e, j as usize, 0.5 * t[c]));
                }
            }
            let pic = Csr::from_triplets(ne, nn, &trip)
                .map_err(|_| SolverError::InvalidArgument("failed to assemble AmsF32 Pi operator"))?;
            let pict = pic.transpose();
            let a_pi = pict.matmul(&a.matmul(&pic));
            amgs.push(AmgF32::new(&a_pi, opts.amg.clone())?);
            pis.push(CsrF32::from_csr_f64(&pic));
            pits.push(CsrF32::from_csr_f64(&pict));
        }
        let inv_diag_f64: Vec<f64> =
            a.diagonal().into_iter().map(|d| if d == 0.0 { 1.0 } else { 1.0 / d }).collect();
        let lmax = estimate_dinv_a_spectral_radius(a, &inv_diag_f64) as f32;
        let inv_diag: Vec<f32> = inv_diag_f64.into_iter().map(|d| d as f32).collect();

        let [pi0, pi1, pi2]: [CsrF32; 3] = pis
            .try_into()
            .map_err(|_| SolverError::InvalidArgument("internal Pi count mismatch"))?;
        let [pit0, pit1, pit2]: [CsrF32; 3] = pits
            .try_into()
            .map_err(|_| SolverError::InvalidArgument("internal Pit count mismatch"))?;
        let [amg0, amg1, amg2]: [AmgF32; 3] = amgs
            .try_into()
            .map_err(|_| SolverError::InvalidArgument("internal AmgF32 count mismatch"))?;

        Ok(Self {
            a: CsrF32::from_csr_f64(a),
            inv_diag,
            lmax,
            smoother: opts.smoother,
            sweeps: opts.smoother_sweeps,
            g: CsrF32::from_csr_f64(g),
            gt: CsrF32::from_csr_f64(&gt_f64),
            amg_g,
            pi: [pi0, pi1, pi2],
            pit: [pit0, pit1, pit2],
            amg_pi: [amg0, amg1, amg2],
            scratch: RefCell::new(ScratchAmsF32 {
                rhs: vec![0.0; ne],
                x: vec![0.0; ne],
                res: vec![0.0; ne],
                tmp_edge: vec![0.0; ne],
                tmp_node: vec![0.0; nn],
                corr_node: vec![0.0; nn],
            }),
        })
    }

    fn correct_f32(&self, parts: &[(&CsrF32, &CsrF32, &AmgF32)], sc: &mut ScratchAmsF32) {
        let ne = self.a.nrows;
        let par = par::worth_parallel(ne);
        self.a.spmv(&sc.x, &mut sc.res);
        par::for_chunks_mut(&mut sc.res, CHUNK, par, |c, rc| {
            let s = c * CHUNK;
            for (i, ri) in rc.iter_mut().enumerate() {
                *ri = sc.rhs[s + i] - *ri;
            }
        });
        for &(q, qt, amg) in parts {
            qt.spmv(&sc.res, &mut sc.tmp_node);
            amg.apply_f32(&sc.tmp_node, &mut sc.corr_node);
            q.spmv(&sc.corr_node, &mut sc.tmp_edge);
            par::for_chunks_mut(&mut sc.x, CHUNK, par, |c, xc| {
                let s = c * CHUNK;
                for (i, xi) in xc.iter_mut().enumerate() {
                    *xi += sc.tmp_edge[s + i];
                }
            });
        }
    }
}

impl Preconditioner<f64> for AmsF32 {
    fn apply(&self, b: &[f64], z: &mut [f64]) {
        let (to_s, from_s) = f32_safe_scale(b);
        let mut guard = self.scratch.borrow_mut();
        let sc = &mut *guard;
        for (dst, &src) in sc.rhs.iter_mut().zip(b) {
            *dst = (src * to_s) as f32;
        }
        sc.x.fill(0.0);
        match self.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.sweeps {
                    gauss_seidel_f32(&self.a, &self.inv_diag, &sc.rhs, &mut sc.x, true);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth_f32(
                    &self.a,
                    &self.inv_diag,
                    self.lmax,
                    degree.max(1),
                    &sc.rhs,
                    &mut sc.x,
                    &mut sc.res,
                    &mut sc.tmp_edge,
                    true,
                );
            }
        }
        let grad = [(&self.g, &self.gt, &self.amg_g)];
        let pis = [0, 1, 2].map(|c| (&self.pi[c], &self.pit[c], &self.amg_pi[c]));
        self.correct_f32(&grad, sc);
        self.correct_f32(&pis, sc);
        self.correct_f32(&grad, sc);
        match self.smoother {
            SmootherKind::GaussSeidel => {
                for _ in 0..self.sweeps {
                    gauss_seidel_f32(&self.a, &self.inv_diag, &sc.rhs, &mut sc.x, false);
                }
            }
            SmootherKind::Chebyshev { degree } => {
                chebyshev_smooth_f32(
                    &self.a,
                    &self.inv_diag,
                    self.lmax,
                    degree.max(1),
                    &sc.rhs,
                    &mut sc.x,
                    &mut sc.res,
                    &mut sc.tmp_edge,
                    false,
                );
            }
        }
        for (dst, &src) in z.iter_mut().zip(&sc.x) {
            *dst = (src as f64) * from_s;
        }
    }
}

impl crate::linalg::device::DevicePreconditioner<CpuVec<f64>> for AmsF32 {
    #[inline]
    fn apply(&self, r: &CpuVec<f64>, z: &mut CpuVec<f64>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}

// ------------------------------------------------------------------------------------------------------------ P2 Nédélec AMS

/// Topological operators connecting a second-order ($P_2$) Nédélec `H(curl)` space to its
/// lowest-order ($P_1$) subspace and scalar `H^1` gradient spaces, for [`AmsP2`].
#[derive(Debug, Clone)]
pub struct AmsP2Topology {
    /// Prolongation `P_{1→2}` (`ndofs_p2 × n_edges_p1`) embedding lowest-order $P_1$ edge DOFs
    /// into the $P_2$ space (for hierarchical $P_2$ Nédélec bases, `P_{2*e, e} = 1`).
    pub p1_to_p2: Csr<f64>,
    /// Lowest-order discrete gradient `G_{P1}` (`n_edges_p1 × n_nodes_p1`) on the $P_1$ subspace.
    pub gradient_p1: Csr<f64>,
    /// Optional full $P_2$ discrete gradient `G_{P2}` (`ndofs_p2 × n_scalar_p2`) mapping the $P_2$
    /// scalar `H^1` space (interior vertices + interior edge bubbles) into `ker(curl) ⊂ V_{P2}`.
    /// Strongly recommended for curl-curl-dominated or shifted problems so higher-order edge-bubble
    /// gradients are eliminated in `O(1)` iterations.
    pub gradient_p2: Option<Csr<f64>>,
    /// Edge vectors `x_b − x_a` for every interior $P_1$ edge (`n_edges_p1`), matching the row
    /// orientation of `gradient_p1`.
    pub edge_vectors: Vec<[f64; 3]>,
    /// Optional disjoint local DOF blocks (e.g. 2-DOF edge pairs `[2*e, 2*e+1]` and 2-DOF face
    /// pairs `[2*ne+2*f, 2*ne+2*f+1]`) for symmetric block Gauss–Seidel smoothing. When empty,
    /// scalar symmetric Gauss–Seidel is used.
    pub blocks: Vec<Vec<usize>>,
}

#[derive(Debug, Clone)]
pub struct AmsP2Options {
    /// Options of the coarse lowest-order $P_1$ [`Ams`] preconditioner.
    pub coarse_ams: AmsOptions,
    /// Options of the AMG hierarchy on the full $P_2$ gradient space `G_{P2}ᵀ A G_{P2}`.
    pub amg_g2: AmgOptions,
    /// Fine-level smoother sweeps (each direction: forward before, backward after).
    pub smoother_sweeps: usize,
    /// Number of coarse $P_1$ AMS V-cycles per $P_2$ cycle (default 1).
    pub coarse_sweeps: usize,
}

impl Default for AmsP2Options {
    fn default() -> Self {
        Self {
            coarse_ams: AmsOptions::default(),
            amg_g2: AmgOptions::default(),
            smoother_sweeps: 1,
            coarse_sweeps: 1,
        }
    }
}

/// Dense inverse of one small diagonal block `A_{bb}` for symmetric block Gauss–Seidel.
struct BlockInv<T: Scalar> {
    dofs: Vec<usize>,
    /// Row-major `m × m` inverse of `A_{bb}`.
    inv: Vec<T>,
}

fn invert_small_dense<T: Scalar>(a: &[T], m: usize) -> Option<Vec<T>> {
    if m == 0 {
        return Some(Vec::new());
    }
    if m == 1 {
        let d = a[0];
        if d.modulus() <= 1e-300 {
            return None;
        }
        return Some(vec![T::ONE / d]);
    }
    if m == 2 {
        let (a00, a01, a10, a11) = (a[0], a[1], a[2], a[3]);
        let det = a00 * a11 - a01 * a10;
        let scale = (a00.modulus() * a11.modulus() + a01.modulus() * a10.modulus()).max(1e-300);
        if det.modulus() <= 1e-14 * scale {
            return None;
        }
        let inv_det = T::ONE / det;
        return Some(vec![a11 * inv_det, -(a01 * inv_det), -(a10 * inv_det), a00 * inv_det]);
    }
    // General small m×m Gauss–Jordan with partial pivoting.
    let mut aug = vec![T::ZERO; m * 2 * m];
    for i in 0..m {
        for j in 0..m {
            aug[i * 2 * m + j] = a[i * m + j];
        }
        aug[i * 2 * m + m + i] = T::ONE;
    }
    for col in 0..m {
        let mut piv_row = col;
        let mut piv_abs = aug[col * 2 * m + col].modulus();
        for r in col + 1..m {
            let v = aug[r * 2 * m + col].modulus();
            if v > piv_abs {
                piv_abs = v;
                piv_row = r;
            }
        }
        if piv_abs <= 1e-300 || !piv_abs.is_finite() {
            return None;
        }
        if piv_row != col {
            for j in 0..2 * m {
                aug.swap(col * 2 * m + j, piv_row * 2 * m + j);
            }
        }
        let inv_p = T::ONE / aug[col * 2 * m + col];
        for j in col..2 * m {
            aug[col * 2 * m + j] *= inv_p;
        }
        for r in 0..m {
            if r == col {
                continue;
            }
            let f = aug[r * 2 * m + col];
            if f == T::ZERO {
                continue;
            }
            for j in col..2 * m {
                let v = aug[col * 2 * m + j];
                aug[r * 2 * m + j] -= f * v;
            }
        }
    }
    let mut inv = vec![T::ZERO; m * m];
    for i in 0..m {
        for j in 0..m {
            inv[i * m + j] = aug[i * 2 * m + m + j];
        }
    }
    Some(inv)
}

fn build_block_smoother<T: Scalar>(a: &Csr<T>, blocks: &[Vec<usize>]) -> Result<Option<Vec<BlockInv<T>>>, SolverError> {
    if blocks.is_empty() {
        return Ok(None);
    }
    let n = a.nrows();
    let mut covered = vec![false; n];
    let mut out = Vec::with_capacity(blocks.len());
    for blk in blocks {
        if blk.is_empty() {
            continue;
        }
        let m = blk.len();
        for &d in blk {
            if d >= n {
                return Err(SolverError::InvalidArgument("block smoother DOF index out of range"));
            }
            covered[d] = true;
        }
        let mut dense = vec![T::ZERO; m * m];
        for (p, &row_dof) in blk.iter().enumerate() {
            let (cols, vals) = a.row(row_dof);
            for (&c, &v) in cols.iter().zip(vals) {
                let c = c as usize;
                for (q, &col_dof) in blk.iter().enumerate() {
                    if c == col_dof {
                        dense[p * m + q] += v;
                    }
                }
            }
        }
        let inv = match invert_small_dense(&dense, m) {
            Some(inv) => inv,
            None => {
                // Singular local block fallback: diagonal inverse
                let mut d_inv = vec![T::ZERO; m * m];
                for p in 0..m {
                    let diag = dense[p * m + p];
                    d_inv[p * m + p] = if diag == T::ZERO { T::ONE } else { T::ONE / diag };
                }
                d_inv
            }
        };
        out.push(BlockInv { dofs: blk.clone(), inv });
    }
    // Cover any remaining DOFs with 1×1 blocks so every unknown is smoothed.
    for (d, &is_cov) in covered.iter().enumerate() {
        if !is_cov {
            let (cols, vals) = a.row(d);
            let diag = cols.iter().zip(vals).filter(|&(&c, _)| c as usize == d).fold(T::ZERO, |acc, (_, &v)| acc + v);
            let inv = if diag == T::ZERO { T::ONE } else { T::ONE / diag };
            out.push(BlockInv { dofs: vec![d], inv: vec![inv] });
        }
    }
    Ok(Some(out))
}

fn block_gauss_seidel<T: Scalar>(a: &Csr<T>, blocks: &[BlockInv<T>], b: &[T], x: &mut [T], forward: bool) {
    let mut r_buf = Vec::new();
    let mut sweep = |blk: &BlockInv<T>| {
        let m = blk.dofs.len();
        r_buf.clear();
        for &row_dof in &blk.dofs {
            let (cols, vals) = a.row(row_dof);
            let mut s = b[row_dof];
            for (&j, &v) in cols.iter().zip(vals) {
                s -= v * x[j as usize];
            }
            r_buf.push(s);
        }
        for p in 0..m {
            let mut dx = T::ZERO;
            let row_off = p * m;
            for q in 0..m {
                dx += blk.inv[row_off + q] * r_buf[q];
            }
            x[blk.dofs[p]] += dx;
        }
    };
    if forward {
        blocks.iter().for_each(&mut sweep);
    } else {
        blocks.iter().rev().for_each(&mut sweep);
    }
}

struct ScratchP2<T> {
    x: Vec<T>,
    res: Vec<T>,
    tmp_p2: Vec<T>,
    tmp_p1: Vec<T>,
    corr_p1: Vec<T>,
    tmp_g2: Vec<T>,
    corr_g2: Vec<T>,
}

/// **Auxiliary-space Maxwell solver for second-order ($P_2$) Nédélec `H(curl)` spaces** (p-multigrid
/// `P_1 → P_2` coarse-space [`Ams`] + full $P_2$ scalar gradient AMG + symmetric block/point
/// Gauss–Seidel smoother).
///
/// The cycle is strictly symmetric (`0 − G₂ − P₁(AMS) − G₂ − 0`), so `AmsP2` drives PCG (SPD),
/// COCG (complex symmetric), FGMRES/GMRES (indefinite shifted curl-curl, PML, driven solves), and
/// Jacobi–Davidson interior eigensolvers with `O(N)` memory.
pub struct AmsP2<T: FaerScalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    blocks: Option<Vec<BlockInv<T>>>,
    p1_to_p2: Csr<T>,
    p2_to_p1: Csr<T>,
    coarse_ams: Ams<T>,
    g2: Option<(Csr<T>, Csr<T>, Amg<T>)>,
    sweeps: usize,
    coarse_sweeps: usize,
    scratch: RefCell<ScratchP2<T>>,
}

impl<T: FaerScalar> AmsP2<T> {
    /// Builds the $P_2$ AMS preconditioner for `a` (`ndofs_p2 × ndofs_p2`, SPD or complex symmetric
    /// auxiliary operator such as `K + βM` or shifted-Laplacian `K − (σ + iγ)M`).
    pub fn new(a: &Csr<T>, topo: &AmsP2Topology, opts: &AmsP2Options) -> Result<Self, SolverError> {
        let n = a.nrows();
        if a.ncols() != n {
            return Err(SolverError::InvalidArgument("AmsP2 matrix must be square"));
        }
        if topo.p1_to_p2.nrows() != n {
            return Err(SolverError::InvalidArgument("p1_to_p2 row count must match P2 matrix dimension"));
        }
        let ne = topo.p1_to_p2.ncols();
        if topo.gradient_p1.nrows() != ne || topo.edge_vectors.len() != ne {
            return Err(SolverError::InvalidArgument("gradient_p1 and edge_vectors must match P1 edge count"));
        }

        let p1_to_p2: Csr<T> = topo.p1_to_p2.map(T::from_f64);
        let p2_to_p1: Csr<T> = topo.p1_to_p2.transpose().map(T::from_f64);
        let a_p1 = p2_to_p1.matmul(&a.matmul(&p1_to_p2));
        let coarse_ams = Ams::new(&a_p1, &topo.gradient_p1, &topo.edge_vectors, &opts.coarse_ams)?;

        let (g2, ng2) = match &topo.gradient_p2 {
            Some(g) => {
                if g.nrows() != n {
                    return Err(SolverError::InvalidArgument("gradient_p2 row count must match P2 matrix dimension"));
                }
                let ng = g.ncols();
                let gm: Csr<T> = g.map(T::from_f64);
                let gtm: Csr<T> = g.transpose().map(T::from_f64);
                let a_g2 = gtm.matmul(&a.matmul(&gm));
                let amg_g2 = Amg::new(&a_g2, opts.amg_g2.clone())?;
                (Some((gm, gtm, amg_g2)), ng)
            }
            None => (None, 0),
        };

        let blocks = build_block_smoother(a, &topo.blocks)?;
        let inv_diag: Vec<T> =
            a.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();

        let z = |k: usize| vec![T::ZERO; k];
        Ok(Self {
            a: a.clone(),
            inv_diag,
            blocks,
            p1_to_p2,
            p2_to_p1,
            coarse_ams,
            g2,
            sweeps: opts.smoother_sweeps,
            coarse_sweeps: opts.coarse_sweeps.max(1),
            scratch: RefCell::new(ScratchP2 {
                x: z(n),
                res: z(n),
                tmp_p2: z(n),
                tmp_p1: z(ne),
                corr_p1: z(ne),
                tmp_g2: z(ng2),
                corr_g2: z(ng2),
            }),
        })
    }

    pub fn dim(&self) -> usize {
        self.a.nrows()
    }

    /// Approximate total sparse matrix non-zeros stored across this $P_2$ AMS preconditioner
    /// (operator `A`, prolongation/restriction, full $P_2$ gradient AMG, and coarse $P_1$ AMS).
    pub fn nnz_storage(&self) -> usize {
        let g2_nnz = self.g2.as_ref().map_or(0, |(g, gt, amg)| {
            g.nnz() + gt.nnz() + amg.stats().level_nnz.iter().sum::<usize>()
        });
        let blk_entries = self.blocks.as_ref().map_or(0, |bs| bs.iter().map(|b| b.inv.len()).sum());
        self.a.nnz() + self.p1_to_p2.nnz() + self.p2_to_p1.nnz() + self.coarse_ams.nnz_storage() + g2_nnz + blk_entries
    }

    fn smooth(&self, b: &[T], x: &mut [T], forward: bool) {
        if let Some(blocks) = &self.blocks {
            block_gauss_seidel(&self.a, blocks, b, x, forward);
        } else {
            gauss_seidel(&self.a, &self.inv_diag, b, x, forward);
        }
    }

    fn compute_residual(&self, b: &[T], sc: &mut ScratchP2<T>) {
        self.a.spmv(&sc.x, &mut sc.res);
        for (r, bb) in sc.res.iter_mut().zip(b) {
            *r = *bb - *r;
        }
    }

    fn correct_g2(&self, b: &[T], sc: &mut ScratchP2<T>) {
        let Some((g, gt, amg)) = &self.g2 else { return };
        self.compute_residual(b, sc);
        gt.spmv(&sc.res, &mut sc.tmp_g2);
        amg.apply(&sc.tmp_g2, &mut sc.corr_g2);
        g.spmv(&sc.corr_g2, &mut sc.tmp_p2);
        for (xi, ti) in sc.x.iter_mut().zip(&sc.tmp_p2) {
            *xi += *ti;
        }
    }

    fn correct_p1(&self, b: &[T], sc: &mut ScratchP2<T>) {
        for _ in 0..self.coarse_sweeps {
            self.compute_residual(b, sc);
            self.p2_to_p1.spmv(&sc.res, &mut sc.tmp_p1);
            self.coarse_ams.apply(&sc.tmp_p1, &mut sc.corr_p1);
            self.p1_to_p2.spmv(&sc.corr_p1, &mut sc.tmp_p2);
            for (xi, ti) in sc.x.iter_mut().zip(&sc.tmp_p2) {
                *xi += *ti;
            }
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for AmsP2<T> {
    fn apply(&self, b: &[T], z: &mut [T]) {
        let mut guard = self.scratch.borrow_mut();
        let sc = &mut *guard;
        sc.x.fill(T::ZERO);
        for _ in 0..self.sweeps {
            self.smooth(b, &mut sc.x, true);
        }
        self.correct_g2(b, sc);
        self.correct_p1(b, sc);
        self.correct_g2(b, sc);
        for _ in 0..self.sweeps {
            self.smooth(b, &mut sc.x, false);
        }
        z.copy_from_slice(&sc.x);
    }
}

impl<T: FaerScalar> crate::linalg::device::DevicePreconditioner<CpuVec<T>> for AmsP2<T> {
    #[inline]
    fn apply(&self, r: &CpuVec<T>, z: &mut CpuVec<T>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}

// ------------------------------------------------------------------------------------------------------------ Deterministic AMS + FGMRES/GMRES solver (FR-33)

/// Options for [`AmsFgmresSolver`] (deterministic right-preconditioned FGMRES/GMRES driven by
/// [`Ams`] for $P_1$ or [`AmsP2`] for $P_2$).
#[derive(Debug, Clone)]
pub struct AmsFgmresOptions {
    /// Target relative residual `‖A x − b‖₂ / ‖b‖₂` (default `1e-9`).
    pub rtol: f64,
    /// Absolute residual stopping tolerance (default `1e-14`).
    pub atol: f64,
    /// Maximum Krylov iterations per RHS column (default `500`).
    pub max_it: usize,
    /// Restart dimension `m` of FGMRES/GMRES(m) (default `50`).
    pub restart: usize,
    /// Use Flexible GMRES (`true`, default) instead of standard right-preconditioned GMRES.
    pub flexible: bool,
    /// Options of the $P_1$ Hiptmair–Xu [`Ams`] preconditioner.
    pub ams: AmsOptions,
    /// Options of the $P_2$ Hiptmair–Xu [`AmsP2`] preconditioner.
    pub ams_p2: AmsP2Options,
}

impl Default for AmsFgmresOptions {
    fn default() -> Self {
        Self {
            rtol: 1e-9,
            atol: 1e-14,
            max_it: 500,
            restart: 50,
            flexible: true,
            ams: AmsOptions::default(),
            ams_p2: AmsP2Options::default(),
        }
    }
}

/// Summary statistics of one [`AmsFgmresSolver::solve_one`] or [`AmsFgmresSolver::solve_many`] call.
#[derive(Debug, Clone, Default)]
pub struct AmsFgmresStats {
    /// Total Krylov iterations summed across all solved RHS columns.
    pub total_iterations: usize,
    /// Worst independently recomputed relative residual `max_j ‖A x_j − b_j‖₂ / ‖b_j‖₂`.
    pub max_relative_residual: f64,
    /// `true` iff every RHS column met the convergence tolerance.
    pub all_converged: bool,
    /// Approximate total sparse non-zeros stored by the `Ams` / `AmsP2` preconditioner hierarchy.
    pub preconditioner_nnz: usize,
}

enum AmsTopologyKind {
    P1 { g: Csr<f64>, edge_vectors: Vec<[f64; 3]> },
    P2(AmsP2Topology),
}

enum AmsPrecondKind<T: FaerScalar> {
    P1(Box<Ams<T>>),
    P2(Box<AmsP2<T>>),
}

impl<T: FaerScalar> AmsPrecondKind<T> {
    fn nnz_storage(&self) -> usize {
        match self {
            Self::P1(a) => a.nnz_storage(),
            Self::P2(a) => a.nnz_storage(),
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for AmsPrecondKind<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        match self {
            Self::P1(a) => a.apply(r, z),
            Self::P2(a) => a.apply(r, z),
        }
    }
}

/// Assembles the complex-symmetric driven curl-curl operator and its shifted-Laplacian auxiliary
/// preconditioner matrix:
/// - `A = K + i k_0 C - k_0^2 M`
/// - `M_pc = K + i k_0 C + k_0^2 M`
///
/// Using `+k_0^2 M` in `M_pc` keeps the real part of the auxiliary nodal operators (`G^T M_pc G`
/// and `Π_c^T M_pc Π_c`) elliptic/positive definite across the frequency band so `Amg<C64>` stays
/// non-singular while clustering the gradient kernel of `M_pc^{-1} A` tightly near `-1`.
pub fn assemble_shifted_curl_curl(
    k: &Csr<C64>,
    c_damp: Option<&Csr<C64>>,
    m: &Csr<C64>,
    k0: f64,
) -> (Csr<C64>, Csr<C64>) {
    let k0_sq = k0 * k0;
    let kc = match c_damp {
        Some(c) if k0 != 0.0 => Csr::axpby(C64::ONE, k, C64::new(0.0, k0), c),
        _ => k.clone(),
    };
    let a_sys = Csr::axpby(C64::ONE, &kc, C64::new(-k0_sq, 0.0), m);
    let a_pc = Csr::axpby(C64::ONE, &kc, C64::new(k0_sq, 0.0), m);
    (a_sys, a_pc)
}

/// Deterministic preconditioned FGMRES/GMRES solver driven by [`Ams`] ($P_1$) or [`AmsP2`] ($P_2$)
/// for multi-RHS real or complex-symmetric curl-curl systems `(K + i k_0 C - k_0^2 M) X = B` (`FR-33`).
///
/// All inner products, norms, basis updates, and SpMVs use `femcore::linalg::vecops` / `basis`
/// fixed-chunk reductions (`CHUNK = 2048`), giving bit-identical results across `NYX_THREADS` (`ADR 11`).
pub struct AmsFgmresSolver<T: FaerScalar> {
    a: Csr<T>,
    topo: AmsTopologyKind,
    pc: AmsPrecondKind<T>,
    opts: AmsFgmresOptions,
}

impl<T: FaerScalar> AmsFgmresSolver<T> {
    /// Builds a $P_1$ `Ams`-preconditioned FGMRES/GMRES solver for `a`.
    ///
    /// Pass `pc_matrix = Some(m_pc)` (e.g. `K + i k_0 C + k_0^2 M`) when `a` is indefinite,
    /// or `None` to build [`Ams`] directly on `a`.
    pub fn new_p1(
        a: &Csr<T>,
        pc_matrix: Option<&Csr<T>>,
        g: &Csr<f64>,
        edge_vectors: &[[f64; 3]],
        opts: AmsFgmresOptions,
    ) -> Result<Self, SolverError> {
        let n = a.nrows();
        if n == 0 || a.ncols() != n {
            return Err(SolverError::InvalidArgument("AmsFgmresSolver matrix must be square and non-empty"));
        }
        if g.nrows() != n || edge_vectors.len() != n {
            return Err(SolverError::InvalidArgument("gradient and edge_vectors must match matrix dimension"));
        }
        let pc_mat = pc_matrix.unwrap_or(a);
        if pc_mat.nrows() != n || pc_mat.ncols() != n {
            return Err(SolverError::InvalidArgument("pc_matrix must match matrix dimension"));
        }
        let ams = Ams::new(pc_mat, g, edge_vectors, &opts.ams)?;
        Ok(Self {
            a: a.clone(),
            topo: AmsTopologyKind::P1 { g: g.clone(), edge_vectors: edge_vectors.to_vec() },
            pc: AmsPrecondKind::P1(Box::new(ams)),
            opts,
        })
    }

    /// Builds a $P_2$ `AmsP2`-preconditioned FGMRES/GMRES solver for `a`.
    pub fn new_p2(
        a: &Csr<T>,
        pc_matrix: Option<&Csr<T>>,
        topo: &AmsP2Topology,
        opts: AmsFgmresOptions,
    ) -> Result<Self, SolverError> {
        let n = a.nrows();
        if n == 0 || a.ncols() != n {
            return Err(SolverError::InvalidArgument("AmsFgmresSolver matrix must be square and non-empty"));
        }
        let pc_mat = pc_matrix.unwrap_or(a);
        if pc_mat.nrows() != n || pc_mat.ncols() != n {
            return Err(SolverError::InvalidArgument("pc_matrix must match matrix dimension"));
        }
        let ams = AmsP2::new(pc_mat, topo, &opts.ams_p2)?;
        Ok(Self {
            a: a.clone(),
            topo: AmsTopologyKind::P2(topo.clone()),
            pc: AmsPrecondKind::P2(Box::new(ams)),
            opts,
        })
    }

    /// Updates the system matrix `a` (and optional auxiliary preconditioner matrix `pc_matrix`) at a
    /// new frequency while reusing the stored discrete gradient / prolongation topology.
    pub fn update_matrix(&mut self, a: &Csr<T>, pc_matrix: Option<&Csr<T>>) -> Result<(), SolverError> {
        let n = self.a.nrows();
        if a.nrows() != n || a.ncols() != n {
            return Err(SolverError::InvalidArgument("updated matrix dimension does not match solver dimension"));
        }
        let pc_mat = pc_matrix.unwrap_or(a);
        if pc_mat.nrows() != n || pc_mat.ncols() != n {
            return Err(SolverError::InvalidArgument("pc_matrix must match solver dimension"));
        }
        self.pc = match &self.topo {
            AmsTopologyKind::P1 { g, edge_vectors } => {
                AmsPrecondKind::P1(Box::new(Ams::new(pc_mat, g, edge_vectors, &self.opts.ams)?))
            }
            AmsTopologyKind::P2(topo) => AmsPrecondKind::P2(Box::new(AmsP2::new(pc_mat, topo, &self.opts.ams_p2)?)),
        };
        self.a = a.clone();
        Ok(())
    }

    pub fn dim(&self) -> usize {
        self.a.nrows()
    }

    pub fn preconditioner_nnz(&self) -> usize {
        self.pc.nnz_storage()
    }

    pub fn options(&self) -> &AmsFgmresOptions {
        &self.opts
    }

    /// Solves `A x = rhs` for a single right-hand side using deterministic FGMRES/GMRES.
    pub fn solve_one(&self, rhs: &[T]) -> (Vec<T>, AmsFgmresStats) {
        let n = self.a.nrows();
        assert_eq!(rhs.len(), n, "rhs length must match matrix dimension");
        let gmres_opts = GmresOptions {
            tol: Tolerances { rtol: self.opts.rtol, atol: self.opts.atol, max_it: self.opts.max_it, ..Default::default() },
            restart: self.opts.restart.max(2),
            side: PcSide::Right,
            orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
            flexible: self.opts.flexible,
            ..Default::default()
        };
        let mut gmres = Gmres::<T>::new(gmres_opts);
        let mut x = vec![T::ZERO; n];
        let info = gmres.solve(&self.a, &self.pc, rhs, &mut x);
        let mut r = vec![T::ZERO; n];
        self.a.spmv(&x, &mut r);
        vecops::aypx(-T::ONE, rhs, &mut r);
        let b_norm = vecops::norm2(rhs);
        let rel = if b_norm == 0.0 { vecops::norm2(&r) } else { vecops::norm2(&r) / b_norm };
        let converged = info.reason.is_converged() && rel <= (self.opts.rtol * 10.0).max(self.opts.atol);
        (
            x,
            AmsFgmresStats {
                total_iterations: info.iterations,
                max_relative_residual: rel,
                all_converged: converged,
                preconditioner_nnz: self.pc.nnz_storage(),
            },
        )
    }

    /// Solves `A X = B` across multiple right-hand-side columns using deterministic FGMRES/GMRES.
    pub fn solve_many(&self, rhs_cols: &[Vec<T>]) -> (Vec<Vec<T>>, AmsFgmresStats) {
        let mut xs = Vec::with_capacity(rhs_cols.len());
        let mut total_iterations = 0usize;
        let mut max_relative_residual = 0.0f64;
        let mut all_converged = true;
        for b in rhs_cols {
            let (x, st) = self.solve_one(b);
            total_iterations += st.total_iterations;
            max_relative_residual = max_relative_residual.max(st.max_relative_residual);
            all_converged &= st.all_converged;
            xs.push(x);
        }
        (
            xs,
            AmsFgmresStats {
                total_iterations,
                max_relative_residual,
                all_converged,
                preconditioner_nnz: self.pc.nnz_storage(),
            },
        )
    }
}

impl AmsFgmresSolver<C64> {
    /// Builds a $P_1$ `Ams<C64>` FGMRES/GMRES solver for `(K + i k_0 C - k_0^2 M) X = B`
    /// preconditioned by the shifted-Laplacian auxiliary matrix `K + i k_0 C + k_0^2 M`.
    pub fn new_shifted_complex_p1(
        k: &Csr<C64>,
        c_damp: Option<&Csr<C64>>,
        m: &Csr<C64>,
        k0: f64,
        g: &Csr<f64>,
        edge_vectors: &[[f64; 3]],
        opts: AmsFgmresOptions,
    ) -> Result<Self, SolverError> {
        let (a_sys, a_pc) = assemble_shifted_curl_curl(k, c_damp, m, k0);
        Self::new_p1(&a_sys, Some(&a_pc), g, edge_vectors, opts)
    }

    /// Builds a $P_2$ `AmsP2<C64>` FGMRES/GMRES solver for `(K + i k_0 C - k_0^2 M) X = B`
    /// preconditioned by the shifted-Laplacian auxiliary matrix `K + i k_0 C + k_0^2 M`.
    pub fn new_shifted_complex_p2(
        k: &Csr<C64>,
        c_damp: Option<&Csr<C64>>,
        m: &Csr<C64>,
        k0: f64,
        topo: &AmsP2Topology,
        opts: AmsFgmresOptions,
    ) -> Result<Self, SolverError> {
        let (a_sys, a_pc) = assemble_shifted_curl_curl(k, c_damp, m, k0);
        Self::new_p2(&a_sys, Some(&a_pc), topo, opts)
    }
}

/// Convenience helper solving `A X = B` for multiple RHS columns using $P_1$ [`Ams`] + deterministic FGMRES.
pub fn solve_ams_fgmres_many<T: FaerScalar>(
    a: &Csr<T>,
    pc_matrix: Option<&Csr<T>>,
    g: &Csr<f64>,
    edge_vectors: &[[f64; 3]],
    rhs_cols: &[Vec<T>],
    opts: &AmsFgmresOptions,
) -> Result<(Vec<Vec<T>>, AmsFgmresStats), SolverError> {
    let solver = AmsFgmresSolver::new_p1(a, pc_matrix, g, edge_vectors, opts.clone())?;
    Ok(solver.solve_many(rhs_cols))
}

/// Convenience helper solving `A X = B` for multiple RHS columns using $P_2$ [`AmsP2`] + deterministic FGMRES.
pub fn solve_ams_p2_fgmres_many<T: FaerScalar>(
    a: &Csr<T>,
    pc_matrix: Option<&Csr<T>>,
    topo: &AmsP2Topology,
    rhs_cols: &[Vec<T>],
    opts: &AmsFgmresOptions,
) -> Result<(Vec<Vec<T>>, AmsFgmresStats), SolverError> {
    let solver = AmsFgmresSolver::new_p2(a, pc_matrix, topo, opts.clone())?;
    Ok(solver.solve_many(rhs_cols))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::maxwell::{nedelec_cube, nedelec_p2_cube};
    use crate::linalg::rng::random_vector;

    #[test]
    fn ams_fgmres_complex_symmetric_p1_multi_rhs_and_deterministic() {
        let cube = nedelec_cube(4);
        let k = cube.stiffness.to_complex();
        let m = cube.mass.to_complex();
        // Complex-symmetric driven curl-curl system: K + i k0 C - k0^2 M (indefinite, A = A^T != A^H)
        let c_damp = Csr::axpby(C64::new(0.4, 0.0), &m, C64::ZERO, &m);
        let k0 = 3.2f64; // k0^2 = 10.24 < first resonance 2*pi^2 ~= 19.74, plus lossy damping
        let opts = AmsFgmresOptions { rtol: 1e-9, max_it: 250, restart: 40, ..Default::default() };

        let mut solver = AmsFgmresSolver::new_shifted_complex_p1(
            &k,
            Some(&c_damp),
            &m,
            k0,
            &cube.gradient,
            &cube.edge_vectors,
            opts.clone(),
        )
        .unwrap();
        let n = solver.dim();
        assert!(solver.preconditioner_nnz() > 0);
        let rhs_cols = vec![
            random_vector::<C64>(n, 11),
            random_vector::<C64>(n, 22),
            random_vector::<C64>(n, 33),
        ];

        let run_in = |nt: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(nt).build().unwrap();
            pool.install(|| {
                let s = AmsFgmresSolver::new_shifted_complex_p1(
                    &k,
                    Some(&c_damp),
                    &m,
                    k0,
                    &cube.gradient,
                    &cube.edge_vectors,
                    opts.clone(),
                )
                .unwrap();
                s.solve_many(&rhs_cols)
            })
        };

        let (sols_1, st_1) = run_in(1);
        assert!(st_1.all_converged, "P1 AMS+FGMRES did not converge: {st_1:?}");
        assert!(st_1.max_relative_residual < 1e-8, "residual {:e}", st_1.max_relative_residual);

        let (sols_8, st_8) = run_in(8);
        assert_eq!(st_1.total_iterations, st_8.total_iterations);
        for (c1, c8) in sols_1.iter().zip(&sols_8) {
            for (a, b) in c1.iter().zip(c8) {
                assert_eq!(a.re.to_bits(), b.re.to_bits(), "re bit mismatch across threads");
                assert_eq!(a.im.to_bits(), b.im.to_bits(), "im bit mismatch across threads");
            }
        }

        // Frequency update at a higher frequency (k0^2 > first resonance = strongly indefinite!)
        let (a2, pc2) = assemble_shifted_curl_curl(&k, Some(&c_damp), &m, 4.8);
        solver.update_matrix(&a2, Some(&pc2)).unwrap();
        let (_sols2, st2) = solver.solve_many(&rhs_cols[..2]);
        assert!(st2.all_converged, "P1 AMS+FGMRES indefinite frequency did not converge: {st2:?}");
        assert!(st2.max_relative_residual < 1e-8);
    }

    #[test]
    fn ams_p2_fgmres_complex_symmetric_multi_rhs() {
        let cube = nedelec_p2_cube(3);
        let k = cube.stiffness.to_complex();
        let m = cube.mass.to_complex();
        let c_damp = Csr::axpby(C64::new(0.4, 0.0), &m, C64::ZERO, &m);
        let topo = AmsP2Topology {
            p1_to_p2: cube.p1_to_p2.clone(),
            gradient_p1: cube.gradient_p1.clone(),
            gradient_p2: Some(cube.gradient_p2.clone()),
            edge_vectors: cube.edge_vectors.clone(),
            blocks: cube.blocks.clone(),
        };
        let opts = AmsFgmresOptions { rtol: 1e-9, max_it: 250, restart: 40, ..Default::default() };
        let solver = AmsFgmresSolver::new_shifted_complex_p2(&k, Some(&c_damp), &m, 3.5, &topo, opts).unwrap();
        let n = solver.dim();
        let rhs_cols = vec![random_vector::<C64>(n, 44), random_vector::<C64>(n, 55)];
        let (_sols, st) = solver.solve_many(&rhs_cols);
        assert!(st.all_converged, "P2 AMS+FGMRES did not converge: {st:?}");
        assert!(st.max_relative_residual < 1e-8, "P2 residual {:e}", st.max_relative_residual);
        assert!(st.preconditioner_nnz > 0);
    }
}
