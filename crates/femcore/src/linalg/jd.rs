//! Jacobi–Davidson interior eigensolver (port item E2; SLEPc's `EPSJD`/`EPSGD` are the nearest
//! references): eigenpairs near a *target* without factorising `A − τB`.
//!
//! * Harmonic Rayleigh–Ritz extraction (default) for interior targets, or standard extraction;
//!   restarts with the best Ritz vectors, locking with deflation.
//! * The correction equation `Q_L (A−θB) Q_R t = −r`, with `Q_R = I − X(BX)ᴴ` and
//!   `Q_L = I − BX·Xᴴ` projecting out the converged vectors `X` and the current one, is solved
//!   approximately by a few steps of left-preconditioned GMRES, with the preconditioner
//!   `K⁻¹ = (I − M⁻¹BX (BXᴴM⁻¹BX)⁻¹ BXᴴ) M⁻¹` built from a user-supplied approximation `M⁻¹` of
//!   `(A−τB)⁻¹` (an incomplete factorisation, AMG on a complex-shifted operator, a low-precision
//!   factorisation, …). **The quality of `M⁻¹` decides everything**; with the identity JD
//!   degrades to a slow subspace iteration.
//! * Supported problems (everything is complex arithmetic; promote real problems with
//!   [`Csr::to_complex`]):
//!   - [`JdOptions::hermitian`]: `A` Hermitian and `B` Hermitian positive definite (or absent):
//!     `B`-orthonormal basis, eigenvectors `B`-orthonormal.
//!   - standard non-Hermitian (`B = None`, e.g. complex symmetric PML): Euclidean basis and *Schur
//!     vector* deflation, the eigenvectors are recovered from the triangular factor.
//!     Non-Hermitian *pencils* (`B ≠ I`, `A` non-Hermitian) are not supported here: use the
//!     contour method or shift-and-invert Krylov–Schur.

use crate::linalg::basis::Basis;
use crate::linalg::csr::Csr;
use crate::linalg::dense::{Dense, DenseError, solve_small, triangular_eigenvector};
use crate::linalg::faer_glue::generalized_eig_c64;
use crate::linalg::ksp::{Gmres, GmresOptions, Tolerances};
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::rng::random_vector;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{dotc, norm2, norm2_sq, scale_real};
use std::cell::RefCell;
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct JdOptions {
    pub nev: usize,
    /// The eigenvalues closest to this point are wanted.
    pub target: C64,
    /// Convergence: `‖Au − θBu‖ ≤ tol·|θ|` (absolute if `|θ|` is tiny).
    pub tol: f64,
    /// Search-space dimension kept after a restart / allowed at most.
    pub m_min: usize,
    pub m_max: usize,
    /// Outer iteration limit.
    pub max_it: usize,
    /// GMRES steps for each correction equation and the relative accuracy asked of them.
    pub inner_its: usize,
    pub inner_tol: f64,
    /// Harmonic extraction (interior targets); `false` = standard Rayleigh–Ritz.
    pub harmonic: bool,
    /// See the module docs.
    pub hermitian: bool,
    pub seed: u64,
}

impl JdOptions {
    pub fn new(nev: usize, target: C64, hermitian: bool) -> Self {
        Self {
            nev,
            target,
            tol: 1e-8,
            m_min: 8,
            m_max: 24,
            max_it: 400,
            inner_its: 12,
            inner_tol: 0.1,
            harmonic: hermitian,
            hermitian,
            seed: 0x1D_5EED,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum JdError {
    Dense(DenseError),
    InvalidArgument(&'static str),
    Unsupported(&'static str),
}

impl std::fmt::Display for JdError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            JdError::Dense(e) => write!(f, "projected problem failed: {e}"),
            JdError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
            JdError::Unsupported(m) => write!(f, "unsupported: {m}"),
        }
    }
}

impl std::error::Error for JdError {}

impl From<DenseError> for JdError {
    fn from(e: DenseError) -> Self {
        JdError::Dense(e)
    }
}

#[derive(Debug, Clone)]
pub struct JdHistoryEntry {
    pub iteration: usize,
    pub nconv: usize,
    pub dim: usize,
    /// Relative residual of the pair currently being refined.
    pub residual: f64,
}

#[derive(Debug)]
pub struct JdResult {
    pub converged: bool,
    /// Eigenvalues in order of convergence.
    pub values: Vec<C64>,
    /// Eigenvectors (columns): `B`-orthonormal (Hermitian) or unit (non-Hermitian).
    pub vectors: Basis<C64>,
    /// Independently recomputed relative residuals `‖Ax − λBx‖/(|λ|‖x‖)`.
    pub residuals: Vec<f64>,
    pub iterations: usize,
    pub matvecs: usize,
    pub pc_applies: usize,
    pub inner_iterations: usize,
    pub history: Vec<JdHistoryEntry>,
    pub seconds: f64,
}

/// `Q_L (A − θB) Q_R` on the complement of the locked and current vectors.
struct ProjectedOp<'a> {
    a: &'a Csr<C64>,
    b: Option<&'a Csr<C64>>,
    theta: C64,
    x: &'a Basis<C64>,
    bx: &'a Basis<C64>,
    k: usize,
    scratch: RefCell<(Vec<C64>, Vec<C64>, Vec<C64>)>,
}

impl LinearOperator<C64> for ProjectedOp<'_> {
    fn dim(&self) -> usize {
        self.a.nrows()
    }
    fn apply(&self, x: &[C64], y: &mut [C64]) {
        let n = x.len();
        let mut sc = self.scratch.borrow_mut();
        let (t, bt, h) = &mut *sc;
        t.copy_from_slice(x);
        // t = Q_R x = x − X (BXᴴ x)
        self.bx.view().dots_block(self.k, x, 1, h, None);
        self.x.view().sub_combination(self.k, h, t);
        // y = (A − θB) t
        self.a.spmv(t, y);
        match self.b {
            Some(b) => b.spmv(t, bt),
            None => bt.copy_from_slice(t),
        }
        for i in 0..n {
            y[i] -= self.theta * bt[i];
        }
        // y = Q_L y = y − BX (Xᴴ y)
        self.x.view().dots_block(self.k, y, 1, h, None);
        self.bx.view().sub_combination(self.k, h, y);
    }
}

/// `K⁻¹ r = M⁻¹r − M⁻¹BX (BXᴴ M⁻¹BX)⁻¹ BXᴴ M⁻¹r`: the preconditioner of the correction equation.
struct ProjectedPc<'a> {
    pc: &'a dyn Preconditioner<C64>,
    bx: &'a Basis<C64>,
    y: Basis<C64>,
    s: Dense<C64>,
    k: usize,
    count: RefCell<usize>,
}

impl Preconditioner<C64> for ProjectedPc<'_> {
    fn apply(&self, r: &[C64], z: &mut [C64]) {
        self.pc.apply(r, z);
        *self.count.borrow_mut() += 1;
        let mut h = vec![C64::ZERO; self.k];
        self.bx.view().dots_block(self.k, z, 1, &mut h, None);
        if let Some(c) = solve_small(&self.s, &h) {
            self.y.view().sub_combination(self.k, &c, z);
        }
    }
}

fn gram_update(g: &mut Dense<C64>, left: &Basis<C64>, right: &Basis<C64>, j: usize) {
    // g[i,j] = left_iᴴ right_j (i ≤ j) and g[j,i] = left_jᴴ right_i (i < j)
    let mut col = vec![C64::ZERO; j + 1];
    left.view().dots_block(j + 1, right.col(j), 1, &mut col, None);
    for (i, v) in col.iter().enumerate() {
        g[(i, j)] = *v;
    }
    if j > 0 {
        let mut row = vec![C64::ZERO; j];
        right.view().dots_block(j, left.col(j), 1, &mut row, None);
        for (i, v) in row.iter().enumerate() {
            g[(j, i)] = v.conj();
        }
    }
}

fn gram_full(g: &mut Dense<C64>, left: &Basis<C64>, right: &Basis<C64>, m: usize) {
    for j in 0..m {
        gram_update(g, left, right, j);
    }
}

pub fn jacobi_davidson(
    a: &Csr<C64>,
    b: Option<&Csr<C64>>,
    pc: &dyn Preconditioner<C64>,
    start: Option<&[C64]>,
    opts: &JdOptions,
) -> Result<JdResult, JdError> {
    let t_start = Instant::now();
    let n = a.nrows();
    if a.ncols() != n || b.is_some_and(|b| b.nrows() != n || b.ncols() != n) {
        return Err(JdError::InvalidArgument("A and B must be square of the same size"));
    }
    if b.is_some() && !opts.hermitian {
        return Err(JdError::Unsupported("non-Hermitian pencils: use the contour or Krylov–Schur solvers"));
    }
    if opts.nev == 0 || opts.nev >= n || opts.m_min < 2 || opts.m_max <= opts.m_min + 1 || opts.m_max + 1 >= n {
        return Err(JdError::InvalidArgument("need 1 ≤ nev < n and 2 ≤ m_min < m_max − 1 < n − 2"));
    }
    if start.is_some_and(|s| s.len() != n) {
        return Err(JdError::InvalidArgument("start vector has wrong length"));
    }
    let weighted = b.is_some(); // B inner product (Hermitian pencil)
    let tau = opts.target;
    let (mmax, nev) = (opts.m_max, opts.nev);

    let mut v = Basis::<C64>::new(n, mmax + 1);
    let mut av = Basis::<C64>::new(n, mmax + 1);
    let mut bv = Basis::<C64>::new(n, mmax + 1);
    let mut xq = Basis::<C64>::new(n, nev + 1); // converged vectors; slot `k` holds the current Ritz vector
    let mut bxq = Basis::<C64>::new(n, nev + 1);
    let mut tmat = Dense::<C64>::zeros(nev, nev); // Schur factor (non-Hermitian)
    let mut lam: Vec<C64> = Vec::new();
    let (mut ha, mut hb) = (Dense::<C64>::zeros(mmax, mmax), Dense::<C64>::zeros(mmax, mmax));
    let (mut gaa, mut gab, mut gbb) =
        (Dense::<C64>::zeros(mmax, mmax), Dense::<C64>::zeros(mmax, mmax), Dense::<C64>::zeros(mmax, mmax));

    let (mut matvecs, mut pc_applies, mut inner_total) = (0usize, 0usize, 0usize);
    let mut history = Vec::new();

    let bmul = |x: &[C64], y: &mut [C64]| match b {
        Some(b) => b.spmv(x, y),
        None => y.copy_from_slice(x),
    };

    // Orthonormalises `t` (and its B-image `bt`) against the locked vectors and the basis, twice.
    let orth = |t: &mut [C64],
                bt: &mut [C64],
                v: &Basis<C64>,
                bv: &Basis<C64>,
                m: usize,
                xq: &Basis<C64>,
                bxq: &Basis<C64>,
                k: usize|
     -> f64 {
        let mut h = vec![C64::ZERO; mmax + nev + 1];
        for _ in 0..2 {
            if k > 0 {
                xq.view().dots_block(k, bt, 1, &mut h, None);
                xq.view().sub_combination(k, &h, t);
                bxq.view().sub_combination(k, &h, bt);
            }
            if m > 0 {
                v.view().dots_block(m, bt, 1, &mut h, None);
                v.view().sub_combination(m, &h, t);
                bv.view().sub_combination(m, &h, bt);
            }
        }
        dotc(t, bt).re.max(0.0).sqrt()
    };

    // Starting vector.
    let mut m = 0usize;
    {
        let mut t = match start {
            Some(s) => s.to_vec(),
            None => random_vector::<C64>(n, opts.seed),
        };
        let mut bt = vec![C64::ZERO; n];
        bmul(&t, &mut bt);
        let nrm = orth(&mut t, &mut bt, &v, &bv, 0, &xq, &bxq, 0);
        if nrm.is_nan() || nrm <= 0.0 {
            return Err(JdError::InvalidArgument("start vector must be non-zero"));
        }
        scale_real(1.0 / nrm, &mut t);
        scale_real(1.0 / nrm, &mut bt);
        v.col_mut(0).copy_from_slice(&t);
        bv.col_mut(0).copy_from_slice(&bt);
        let mut w = vec![C64::ZERO; n];
        a.spmv(&t, &mut w);
        matvecs += 1;
        av.col_mut(0).copy_from_slice(&w);
        m = m.max(1);
        gram_update(&mut ha, &v, &av, 0);
        gram_update(&mut hb, &v, &bv, 0);
        if opts.harmonic {
            gram_update(&mut gaa, &av, &av, 0);
            gram_update(&mut gab, &av, &bv, 0);
            gram_update(&mut gbb, &bv, &bv, 0);
        }
    }

    // Replaces the (empty) search space by one random vector, orthogonal to the locked set.
    macro_rules! reseed {
        ($salt:expr) => {{
            let mut t = random_vector::<C64>(n, opts.seed.wrapping_add(1_000_003 * $salt as u64));
            let mut bt = vec![C64::ZERO; n];
            bmul(&t, &mut bt);
            let nrm = orth(&mut t, &mut bt, &v, &bv, 0, &xq, &bxq, lam.len());
            scale_real(1.0 / nrm, &mut t);
            scale_real(1.0 / nrm, &mut bt);
            v.col_mut(0).copy_from_slice(&t);
            bv.col_mut(0).copy_from_slice(&bt);
            let mut w = vec![C64::ZERO; n];
            a.spmv(&t, &mut w);
            matvecs += 1;
            if !opts.hermitian && !lam.is_empty() {
                let mut h = vec![C64::ZERO; lam.len()];
                xq.view().dots_block(lam.len(), &w, 1, &mut h, None);
                xq.view().sub_combination(lam.len(), &h, &mut w);
            }
            av.col_mut(0).copy_from_slice(&w);
            m = 1;
        }};
    }

    let mut inner = Gmres::<C64>::new(GmresOptions {
        tol: Tolerances { rtol: opts.inner_tol, max_it: opts.inner_its, ..Default::default() },
        restart: opts.inner_its,
        ..Default::default()
    });
    let mut converged_all = false;
    let mut it = 0usize;
    while it < opts.max_it {
        it += 1;
        let k = lam.len();

        // ---- Extraction on the m-dimensional search space.
        let (mu, cvecs): (Vec<C64>, Dense<C64>) = {
            let (da, db) = if opts.harmonic {
                // G c = μ F c with G = WᴴW, F = WᴴBV, W = (A − τB)V.
                let tc = tau;
                let g = Dense::from_fn(m, m, |i, j| {
                    gaa[(i, j)] - tc * gab[(i, j)] - tc.conj() * gab[(j, i)].conj() + tc.norm_sqr() * gbb[(i, j)]
                });
                let f = Dense::from_fn(m, m, |i, j| gab[(i, j)] - tc.conj() * gbb[(i, j)]);
                (g, f)
            } else {
                (Dense::from_fn(m, m, |i, j| ha[(i, j)]), Dense::from_fn(m, m, |i, j| hb[(i, j)]))
            };
            generalized_eig_c64(&da, &db)?
        };
        // Order by closeness to the target (harmonic: smallest |μ|; standard: |θ − τ|).
        let mut order: Vec<usize> = (0..m).filter(|&i| mu[i].re.is_finite()).collect();
        order.sort_by(|&i, &j| {
            let key = |x: usize| if opts.harmonic { mu[x].norm() } else { (mu[x] - tau).norm() };
            key(i).total_cmp(&key(j))
        });
        if order.is_empty() {
            return Err(JdError::Dense(DenseError::NoConvergence));
        }

        // ---- Best Ritz pair.
        let best = order[0];
        let c = cvecs.col(best);
        let (mut u, mut au, mut bu) = (vec![C64::ZERO; n], vec![C64::ZERO; n], vec![C64::ZERO; n]);
        v.view().add_combination(m, c, &mut u);
        av.view().add_combination(m, c, &mut au);
        bv.view().add_combination(m, c, &mut bu);
        let unrm = if weighted { dotc(&u, &bu).re.max(0.0).sqrt() } else { norm2(&u) };
        for vec in [&mut u, &mut au, &mut bu] {
            scale_real(1.0 / unrm, vec);
        }
        let theta = dotc(&u, &au) / dotc(&u, &bu);
        let mut r = au.clone();
        for i in 0..n {
            r[i] -= theta * bu[i];
        }
        let rnorm = norm2(&r);
        let rel = rnorm / theta.norm().max(1e-300);
        history.push(JdHistoryEntry { iteration: it, nconv: k, dim: m, residual: rel });

        if rel <= opts.tol {
            // ---- Lock: u joins the converged set X (and, non-Hermitian, the Schur factor).
            xq.col_mut(k).copy_from_slice(&u);
            bxq.col_mut(k).copy_from_slice(&bu);
            if !opts.hermitian {
                // T[0..=k, k] = Qᴴ (A u) (the deflated Schur column); diagonal = θ.
                let mut col = vec![C64::ZERO; k + 1];
                xq.view().dots_block(k + 1, &au_unprojected(&au, &xq, k, a, &u), 1, &mut col, None);
                for (i, v) in col.iter().enumerate().take(k) {
                    tmat[(i, k)] = *v;
                }
                tmat[(k, k)] = theta;
            }
            if !opts.hermitian {
                // Schur deflation: the stored columns of A·V must live on the complement of u too.
                for i in 0..m {
                    let h = dotc(&u, av.col(i));
                    for (x, y) in av.col_mut(i).iter_mut().zip(&u) {
                        *x -= h * *y;
                    }
                }
            }
            lam.push(theta);
            if lam.len() >= nev {
                converged_all = true;
                break;
            }
            // Restart from the remaining Ritz vectors (the locked direction is dropped).
            let keep: Vec<usize> =
                order.iter().copied().skip(1).take(opts.m_min.min(m.saturating_sub(1)).max(1)).collect();
            m = restart(
                &mut v,
                &mut av,
                &mut bv,
                m,
                &keep,
                &cvecs,
                &Locked { x: &xq, bx: &bxq, lam: &lam, hermitian: opts.hermitian },
            );
            if m == 0 {
                // Everything kept was dependent on the locked set: start over from a random vector.
                reseed!(it);
            }
            rebuild(&mut ha, &mut hb, &mut gaa, &mut gab, &mut gbb, &v, &av, &bv, m, opts.harmonic);
            continue;
        }

        // ---- Restart if the search space is full.
        if m >= mmax {
            let keep: Vec<usize> = order.iter().copied().take(opts.m_min).collect();
            m = restart(
                &mut v,
                &mut av,
                &mut bv,
                m,
                &keep,
                &cvecs,
                &Locked { x: &xq, bx: &bxq, lam: &lam, hermitian: opts.hermitian },
            );
            if m == 0 {
                // Everything kept was dependent on the locked set: start over from a random vector.
                reseed!(it);
            }
            rebuild(&mut ha, &mut hb, &mut gaa, &mut gab, &mut gbb, &v, &av, &bv, m, opts.harmonic);
            // After restarting the previous best vector is span-contained; fall through to expand.
        }

        // ---- Correction equation: Q_L (A − θB) Q_R t = −r, preconditioned.
        xq.col_mut(k).copy_from_slice(&u);
        bxq.col_mut(k).copy_from_slice(&bu);
        let kk = k + 1;
        let mut yb = Basis::<C64>::new(n, kk);
        for j in 0..kk {
            let mut z = vec![C64::ZERO; n];
            pc.apply(bxq.col(j), &mut z);
            pc_applies += 1;
            yb.col_mut(j).copy_from_slice(&z);
        }
        let mut s = Dense::<C64>::zeros(kk, kk);
        {
            let mut col = vec![C64::ZERO; kk];
            for j in 0..kk {
                bxq.view().dots_block(kk, yb.col(j), 1, &mut col, None);
                for (i, v) in col.iter().enumerate() {
                    s[(i, j)] = *v;
                }
            }
        }
        let op = ProjectedOp {
            a,
            b,
            theta,
            x: &xq,
            bx: &bxq,
            k: kk,
            scratch: RefCell::new((vec![C64::ZERO; n], vec![C64::ZERO; n], vec![C64::ZERO; mmax + nev + 2])),
        };
        let kpc = ProjectedPc { pc, bx: &bxq, y: yb, s, k: kk, count: RefCell::new(0) };
        let rhs: Vec<C64> = r.iter().map(|x| -*x).collect();
        let mut t = vec![C64::ZERO; n];
        let info = inner.solve(&op, &kpc, &rhs, &mut t);
        inner_total += info.iterations;
        matvecs += info.matvecs;
        pc_applies += *kpc.count.borrow();

        // ---- Expand the search space with t.
        let mut bt = vec![C64::ZERO; n];
        bmul(&t, &mut bt);
        let mut nrm = orth(&mut t, &mut bt, &v, &bv, m, &xq, &bxq, k);
        if nrm.is_nan() || nrm <= 1e-10 * norm2(&rhs).max(1e-300) || !nrm.is_finite() {
            // Stagnation: continue with a random direction.
            t = random_vector::<C64>(n, opts.seed.wrapping_add(it as u64));
            bmul(&t, &mut bt);
            nrm = orth(&mut t, &mut bt, &v, &bv, m, &xq, &bxq, k);
        }
        scale_real(1.0 / nrm, &mut t);
        scale_real(1.0 / nrm, &mut bt);
        v.col_mut(m).copy_from_slice(&t);
        bv.col_mut(m).copy_from_slice(&bt);
        let mut w = vec![C64::ZERO; n];
        a.spmv(&t, &mut w);
        matvecs += 1;
        if !opts.hermitian && k > 0 {
            // Schur deflation: A is applied on the complement of the locked vectors.
            let mut h = vec![C64::ZERO; k];
            xq.view().dots_block(k, &w, 1, &mut h, None);
            xq.view().sub_combination(k, &h, &mut w);
        }
        av.col_mut(m).copy_from_slice(&w);
        gram_update(&mut ha, &v, &av, m);
        gram_update(&mut hb, &v, &bv, m);
        if opts.harmonic {
            gram_update(&mut gaa, &av, &av, m);
            gram_update(&mut gab, &av, &bv, m);
            gram_update(&mut gbb, &bv, &bv, m);
        }
        m += 1;
    }

    // ---- Results.
    let k = lam.len();
    let mut vectors = Basis::<C64>::new(n, k.max(1));
    let mut values = lam.clone();
    if opts.hermitian {
        for j in 0..k {
            vectors.col_mut(j).copy_from_slice(xq.col(j));
        }
    } else {
        // Schur vectors → eigenvectors through the triangular factor.
        let tt = Dense::from_fn(k, k, |i, j| if i <= j { tmat[(i, j)] } else { C64::ZERO });
        let mut y = Dense::<C64>::zeros(k, k);
        for p in 0..k {
            for (i, c) in triangular_eigenvector(&tt, k, p).iter().enumerate() {
                y[(i, p)] = *c;
            }
        }
        for p in 0..k {
            let mut x = vec![C64::ZERO; n];
            xq.view().add_combination(k, y.col(p), &mut x);
            let nx = norm2(&x);
            scale_real(1.0 / nx, &mut x);
            vectors.col_mut(p).copy_from_slice(&x);
        }
        values = (0..k).map(|p| tt[(p, p)]).collect();
    }
    let mut residuals = Vec::with_capacity(k);
    let (mut ax, mut bx) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    for (j, &l) in values.iter().enumerate() {
        let x = vectors.col(j);
        a.spmv(x, &mut ax);
        bmul(x, &mut bx);
        let r2: f64 = ax.iter().zip(&bx).map(|(p, q)| (*p - l * *q).norm_sqr()).sum();
        residuals.push(r2.sqrt() / (l.norm().max(1e-300) * norm2(x)));
    }
    let _ = norm2_sq(&ax);
    // Safety net: the locked set was accepted on the *recurrence* residual (`A V` is updated, not
    // recomputed), so confirm it with the true residuals computed above.
    let verified = residuals.iter().all(|&x| x <= 100.0 * opts.tol);
    Ok(JdResult {
        converged: converged_all && verified,
        values,
        vectors,
        residuals,
        iterations: it,
        matvecs,
        pc_applies,
        inner_iterations: inner_total,
        history,
        seconds: t_start.elapsed().as_secs_f64(),
    })
}

/// `A u` without the Schur projection, for the triangular column (recomputed to keep `av`'s
/// deflated form separate from the true product).
fn au_unprojected(_au: &[C64], _xq: &Basis<C64>, _k: usize, a: &Csr<C64>, u: &[C64]) -> Vec<C64> {
    let mut w = vec![C64::ZERO; u.len()];
    a.spmv(u, &mut w);
    w
}

/// The converged vectors (and `B` times them) and their values.
struct Locked<'a> {
    x: &'a Basis<C64>,
    bx: &'a Basis<C64>,
    lam: &'a [C64],
    hermitian: bool,
}

/// Replaces the basis by the Ritz vectors `keep` (columns of `cvecs`), re-orthonormalised
/// (`B`-inner product when `weighted`), applying the same transformation to `A V` and `B V`.
fn restart(
    v: &mut Basis<C64>,
    av: &mut Basis<C64>,
    bv: &mut Basis<C64>,
    m: usize,
    keep: &[usize],
    cvecs: &Dense<C64>,
    locked: &Locked,
) -> usize {
    let mk = keep.len();
    let ld = m;
    let q: Vec<C64> = keep.iter().flat_map(|&i| cvecs.col(i).to_vec()).collect();
    v.mult_in_place(0, m, mk, &q, ld);
    av.mult_in_place(0, m, mk, &q, ld);
    bv.mult_in_place(0, m, mk, &q, ld);
    // Modified Gram–Schmidt on the new columns with the inner product ⟨x,y⟩ = xᴴ(B y) = xᴴ·bv_y.
    let mut out = 0usize; // columns kept so far; nearly dependent ones are dropped
    for j in 0..mk {
        let before = dotc(v.col(j), bv.col(j)).re.max(0.0).sqrt();
        for _ in 0..2 {
            // Keep the basis B-orthogonal to the locked vectors: harmonic (and non-normal) Ritz
            // vectors are not orthogonal to the pair that was just locked. The image of the
            // removed component is λ B x (Hermitian) or zero (Schur-deflated `A V`).
            for (i, &l) in locked.lam.iter().enumerate() {
                let h = dotc(locked.x.col(i), bv.col(j));
                for (a, b) in v.col_mut(j).iter_mut().zip(locked.x.col(i)) {
                    *a -= h * *b;
                }
                for (a, b) in bv.col_mut(j).iter_mut().zip(locked.bx.col(i)) {
                    *a -= h * *b;
                }
                if locked.hermitian {
                    for (a, b) in av.col_mut(j).iter_mut().zip(locked.bx.col(i)) {
                        *a -= h * l * *b;
                    }
                }
            }
            for i in 0..out {
                let h = dotc(v.col(i), bv.col(j));
                let (vi, vj) = (v.col(i).to_vec(), v.col_mut(j));
                for (a, b) in vj.iter_mut().zip(&vi) {
                    *a -= h * *b;
                }
                let (ai, aj) = (av.col(i).to_vec(), av.col_mut(j));
                for (a, b) in aj.iter_mut().zip(&ai) {
                    *a -= h * *b;
                }
                let (bi, bj) = (bv.col(i).to_vec(), bv.col_mut(j));
                for (a, b) in bj.iter_mut().zip(&bi) {
                    *a -= h * *b;
                }
            }
        }
        let nrm = dotc(v.col(j), bv.col(j)).re.max(0.0).sqrt();
        if !nrm.is_finite() || nrm <= 1e-10 * before.max(1e-300) {
            continue; // linearly dependent on the columns already kept
        }
        if out != j {
            for basis in [&mut *v, &mut *av, &mut *bv] {
                let src = basis.col(j).to_vec();
                basis.col_mut(out).copy_from_slice(&src);
            }
        }
        for basis in [&mut *v, &mut *av, &mut *bv] {
            scale_real(1.0 / nrm, basis.col_mut(out));
        }
        out += 1;
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn rebuild(
    ha: &mut Dense<C64>,
    hb: &mut Dense<C64>,
    gaa: &mut Dense<C64>,
    gab: &mut Dense<C64>,
    gbb: &mut Dense<C64>,
    v: &Basis<C64>,
    av: &Basis<C64>,
    bv: &Basis<C64>,
    m: usize,
    harmonic: bool,
) {
    gram_full(ha, v, av, m);
    gram_full(hb, v, bv, m);
    if harmonic {
        gram_full(gaa, av, av, m);
        gram_full(gab, av, bv, m);
        gram_full(gbb, bv, bv, m);
    }
}
