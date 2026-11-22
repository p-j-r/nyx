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

use crate::linalg::amg::{Amg, AmgOptions, gauss_seidel};
use crate::linalg::csr::Csr;
use crate::linalg::direct::SolverError;
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::Preconditioner;
use crate::linalg::scalar::Scalar;
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct AmsOptions {
    /// Options of the AMG hierarchies on the auxiliary nodal operators.
    pub amg: AmgOptions,
    /// Gauss–Seidel sweeps of the edge smoother (each direction).
    pub smoother_sweeps: usize,
}

impl Default for AmsOptions {
    fn default() -> Self {
        Self { amg: AmgOptions::default(), smoother_sweeps: 1 }
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
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, true);
        }
        let grad = [(&self.g, &self.gt, &self.amg_g)];
        let pis = [0, 1, 2].map(|c| (&self.pi[c], &self.pit[c], &self.amg_pi[c]));
        self.correct(&grad, b, sc);
        self.correct(&pis, b, sc);
        self.correct(&grad, b, sc);
        for _ in 0..self.sweeps {
            gauss_seidel(&self.a, &self.inv_diag, b, &mut sc.x, false);
        }
        z.copy_from_slice(&sc.x);
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
