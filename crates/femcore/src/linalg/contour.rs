//! Contour-integral eigensolver (port item E3): *all* eigenvalues of `A x = λ B x` inside an
//! elliptical region of the complex plane, by block Sakurai–Sugiura moments with Rayleigh–Ritz
//! extraction (the CIRR variant; SLEPc's `EPSCISS` is the reference).
//!
//! For quadrature nodes `z_j` on the boundary, the shifted systems `(A − z_j B) Y_j = B V` with a
//! block `V` of `L` random vectors are solved independently (each node is one numeric
//! factorisation on a reused symbolic analysis, with multi-right-hand-side solves; nodes can run
//! on several threads). The moments `S_k = Σ_j w_j ((z_j−c)/ρ)^k Y_j` span, up to quadrature
//! error, the eigenspace of the eigenvalues inside the contour; a rank-revealing SVD of
//! `[S_0 … S_{M−1}]` gives an orthonormal basis `Q`, and the small pencil `(QᴴAQ, QᴴBQ)` yields the
//! eigenpairs, which are kept if they lie inside the region and have a small *explicitly
//! recomputed* residual.
//!
//! Two safety nets guard against missed or spurious modes:
//! * [`ContourResult::estimated_count`]: the stochastic trace estimate `tr(P) ≈ Re tr(Vᴴ S_0)/L`
//!   of the number of eigenvalues inside the contour (compare it with the number returned);
//! * [`ContourResult::saturated`]: the moment matrix has full numerical rank `M·L`, so there may be
//!   more eigenvalues than can be captured — increase `block`/`moments` or split the region.
//!
//! Complex arithmetic is used throughout (the nodes are complex), so real problems are promoted
//! to complex and the eigenvectors are returned as complex vectors. Exploiting conjugate symmetry
//! of real problems (half the nodes) is not done.

use crate::linalg::basis::Basis;
use crate::linalg::csr::Csr;
use crate::linalg::dense::{Dense, DenseError};
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError, Threads};
use crate::linalg::faer_glue::{FaerScalar, generalized_eig_c64, thin_svd_c64};
use crate::linalg::rng::SplitMix64;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::norm2;
use std::time::Instant;

/// The region `{ z : ((Re z − c_re)/a)² + ((Im z − c_im)/b)² ≤ 1 }`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    pub center: C64,
    /// Semi-axis along the real direction.
    pub a: f64,
    /// Semi-axis along the imaginary direction.
    pub b: f64,
}

impl Ellipse {
    pub fn circle(center: C64, radius: f64) -> Self {
        Self { center, a: radius, b: radius }
    }

    /// An interval `[lo, hi]` of the real axis, widened into an ellipse with imaginary semi-axis
    /// `aspect · (hi − lo)/2` (a circle for `aspect = 1`).
    pub fn real_interval(lo: f64, hi: f64, aspect: f64) -> Self {
        let a = 0.5 * (hi - lo);
        Self { center: C64::new(0.5 * (lo + hi), 0.0), a, b: aspect * a }
    }

    pub fn contains(&self, z: C64) -> bool {
        let d = z - self.center;
        (d.re / self.a).powi(2) + (d.im / self.b).powi(2) <= 1.0
    }

    /// Node `j` of `n` (trapezoid rule on the boundary, offset by half a step) and its weight
    /// `z'(θ)/(i·n)`, such that `(1/2πi)∮ f dz ≈ Σ_j w_j f(z_j)`.
    pub fn node(&self, j: usize, n: usize) -> (C64, C64) {
        let th = 2.0 * std::f64::consts::PI * (j as f64 + 0.5) / n as f64;
        let z = self.center + C64::new(self.a * th.cos(), self.b * th.sin());
        let dz = C64::new(-self.a * th.sin(), self.b * th.cos());
        (z, dz / C64::new(0.0, n as f64))
    }
}

#[derive(Debug, Clone)]
pub struct ContourOptions {
    pub region: Ellipse,
    /// Quadrature nodes (default 24). More nodes sharpen the filter at the boundary.
    pub nodes: usize,
    /// Random vectors `L` per node (default 8): must be at least the largest multiplicity of an
    /// eigenvalue inside the region, and `moments·block` at least their total number.
    pub block: usize,
    /// Moments `M` (default 4).
    pub moments: usize,
    /// Singular values below `rank_tol · σ₁` are dropped (default 1e-12).
    pub rank_tol: f64,
    /// Pairs with relative residual `‖Ax−λBx‖/(|λ|‖x‖)` above this are rejected (default 1e-8).
    pub residual_tol: f64,
    pub seed: u64,
    /// Nodes factorised and solved concurrently, each with its own factorisation (default 1:
    /// one factorisation at a time, using the shared pool inside faer). Results do not depend on
    /// this number: moments are folded in node order.
    pub workers: usize,
}

impl ContourOptions {
    pub fn new(region: Ellipse) -> Self {
        Self {
            region,
            nodes: 24,
            block: 8,
            moments: 4,
            rank_tol: 1e-12,
            residual_tol: 1e-8,
            seed: 0xC0DE_5EED,
            workers: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ContourError {
    Solver(SolverError),
    Dense(DenseError),
    InvalidArgument(&'static str),
}

impl std::fmt::Display for ContourError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ContourError::Solver(e) => write!(
                f,
                "factorisation failed at a quadrature node: {e} (a node on an eigenvalue? shift the region slightly)"
            ),
            ContourError::Dense(e) => write!(f, "projected problem failed: {e}"),
            ContourError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
        }
    }
}

impl std::error::Error for ContourError {}

impl From<SolverError> for ContourError {
    fn from(e: SolverError) -> Self {
        ContourError::Solver(e)
    }
}
impl From<DenseError> for ContourError {
    fn from(e: DenseError) -> Self {
        ContourError::Dense(e)
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ContourTimings {
    pub factor: f64,
    pub solve: f64,
    pub moments: f64,
    pub svd: f64,
    pub project: f64,
    pub total: f64,
}

#[derive(Debug)]
pub struct ContourResult {
    /// Accepted eigenvalues inside the region (unordered across the region; sorted by real part,
    /// then imaginary part).
    pub values: Vec<C64>,
    /// Unit-norm eigenvectors (columns), in the order of `values`.
    pub vectors: Basis<C64>,
    /// Independently recomputed relative residuals.
    pub residuals: Vec<f64>,
    /// Stochastic (Hutchinson) estimate of the number of eigenvalues inside the contour, and its
    /// standard error from the spread over the `block` random vectors. The error is typically a
    /// few units for `block` around 10: use it to *flag* a shortfall (found ≪ estimate − 3·error),
    /// not to count.
    pub estimated_count: f64,
    pub count_stderr: f64,
    /// Numerical rank of the moment matrix (≤ `moments · block`). It is an *upper bound* on the
    /// number of eigenvalues inside: eigenvalues just outside the contour leak into the moments
    /// at the quadrature-error level (about 1e-6 for 24–48 nodes), above `rank_tol`; those extra
    /// Ritz vectors are discarded by the region and residual tests.
    pub rank: usize,
    /// The moment matrix has full rank `moments · block`: more eigenvalues may exist than were
    /// captured (it can also fire when the leakage alone fills the subspace; then enlarge it).
    pub saturated: bool,
    /// Ritz values inside the region rejected for a large residual (a symptom of too few nodes or
    /// eigenvalues just outside the contour).
    pub rejected: usize,
    pub timings: ContourTimings,
}

/// Solution of one node: the stacked `(zB − A)⁻¹ B V`, factor time, solve time.
type NodeSolve = Result<(Vec<C64>, f64, f64), SolverError>;

/// All eigenvalues of `A x = λ B x` (`B = None`: identity) inside `opts.region`.
pub fn contour_eigs<T: FaerScalar>(
    a: &Csr<T>,
    b: Option<&Csr<T>>,
    opts: &ContourOptions,
) -> Result<ContourResult, ContourError> {
    let t_start = Instant::now();
    let n = a.nrows();
    if opts.nodes < 4 || opts.block == 0 || opts.moments == 0 || opts.workers == 0 {
        return Err(ContourError::InvalidArgument("need at least 4 nodes, one vector, one moment and one worker"));
    }
    if !(opts.region.a > 0.0 && opts.region.b > 0.0) {
        return Err(ContourError::InvalidArgument("ellipse semi-axes must be positive"));
    }
    let (l, mm) = (opts.block, opts.moments);
    if l * mm > n {
        return Err(ContourError::InvalidArgument("moments · block exceeds the problem size; use a dense solver"));
    }
    let ac = a.map(|v| v.to_c64());
    let bc = b.map(|m| m.map(|v| v.to_c64()));
    let mut tm = ContourTimings::default();

    // Random sign block V, and B·V.
    let mut rng = SplitMix64::new(opts.seed);
    let v: Vec<C64> = (0..n * l).map(|_| C64::new(if rng.next_u64() & 1 == 0 { 1.0 } else { -1.0 }, 0.0)).collect();
    let mut bv = v.clone();
    if let Some(bm) = &bc {
        for c in 0..l {
            bm.spmv(&v[c * n..(c + 1) * n], &mut bv[c * n..(c + 1) * n]);
        }
    }

    // One factorisation object per worker (own numeric storage; symbolic analysis once each).
    let workers = opts.workers.min(opts.nodes);
    let mut facs: Vec<ShiftedFactorization<C64>> = Vec::with_capacity(workers);
    for _ in 0..workers {
        let mut f = ShiftedFactorization::new(&ac, bc.as_ref(), Definiteness::General)?;
        if workers > 1 {
            f.threads.factor = Threads::Sequential; // do not stack the workers on top of a faer pool
        }
        facs.push(f);
    }

    let rho = opts.region.a.max(opts.region.b);
    let mut s = vec![C64::ZERO; n * l * mm];
    let mut start = 0;
    while start < opts.nodes {
        let chunk: Vec<usize> = (start..(start + workers).min(opts.nodes)).collect();
        let t0 = Instant::now();
        let solve_node = |fac: &mut ShiftedFactorization<C64>, j: usize| -> Result<(Vec<C64>, f64, f64), SolverError> {
            let (z, _) = opts.region.node(j, opts.nodes);
            let t = Instant::now();
            fac.set_shift(z)?;
            let tf = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let mut y = bv.clone();
            fac.solve_block_in_place(&mut y, l);
            // (zB − A)⁻¹ B V = −(A − zB)⁻¹ B V
            for e in &mut y {
                *e = -*e;
            }
            Ok((y, tf, t.elapsed().as_secs_f64()))
        };
        let results: Vec<NodeSolve> = if workers == 1 {
            vec![solve_node(&mut facs[0], chunk[0])]
        } else {
            std::thread::scope(|sc| {
                let solve_node = &solve_node;
                let handles: Vec<_> = facs
                    .iter_mut()
                    .zip(chunk.iter().copied())
                    .map(|(f, j)| sc.spawn(move || solve_node(f, j)))
                    .collect();
                handles.into_iter().map(|h| h.join().expect("worker thread")).collect()
            })
        };
        let wall = t0.elapsed().as_secs_f64();
        // Fold in node order so that the result does not depend on the number of workers.
        let t1 = Instant::now();
        for (&j, res) in chunk.iter().zip(results) {
            let (y, tf, ts) = res?;
            tm.factor += tf / workers as f64;
            tm.solve += ts / workers as f64;
            let (z, w) = opts.region.node(j, opts.nodes);
            let t = (z - opts.region.center).scale(1.0 / rho);
            let mut coef = w;
            for k in 0..mm {
                let sk = &mut s[k * n * l..(k + 1) * n * l];
                for (si, yi) in sk.iter_mut().zip(&y) {
                    *si += coef * *yi;
                }
                coef *= t;
            }
        }
        tm.moments += t1.elapsed().as_secs_f64();
        let _ = wall;
        start += workers;
    }

    // Eigenvalue count estimate: tr(P) ≈ Re tr(Vᴴ S_0)/L  (S_0 is the first block of `s`).
    let per_vector: Vec<f64> =
        (0..l).map(|c| (0..n).map(|i| (v[c * n + i].conj() * s[c * n + i]).re).sum::<f64>()).collect();
    let est = per_vector.iter().sum::<f64>() / l as f64;
    let count_stderr = if l > 1 {
        (per_vector.iter().map(|x| (x - est).powi(2)).sum::<f64>() / (l as f64 - 1.0) / l as f64).sqrt()
    } else {
        f64::INFINITY
    };

    // Rank-revealing basis of the moment space.
    let t = Instant::now();
    let cols = l * mm;
    let (sigma, u) = thin_svd_c64(&s, n, cols)?;
    tm.svd = t.elapsed().as_secs_f64();
    let r = sigma.iter().take_while(|&&x| x > opts.rank_tol * sigma[0]).count();
    if r == 0 || sigma[0] == 0.0 {
        tm.total = t_start.elapsed().as_secs_f64();
        return Ok(ContourResult {
            values: Vec::new(),
            vectors: Basis::new(n, 1),
            residuals: Vec::new(),
            estimated_count: est,
            count_stderr,
            rank: 0,
            saturated: false,
            rejected: 0,
            timings: tm,
        });
    }
    let mut q = Basis::<C64>::new(n, r);
    for j in 0..r {
        q.col_mut(j).copy_from_slice(&u[j * n..(j + 1) * n]);
    }

    // Rayleigh–Ritz on the pencil.
    let t = Instant::now();
    let (mut aq, mut bq) = (vec![C64::ZERO; n * r], vec![C64::ZERO; n * r]);
    for j in 0..r {
        ac.spmv(q.col(j), &mut aq[j * n..(j + 1) * n]);
        match &bc {
            Some(bm) => bm.spmv(q.col(j), &mut bq[j * n..(j + 1) * n]),
            None => bq[j * n..(j + 1) * n].copy_from_slice(q.col(j)),
        }
    }
    let (mut ha, mut hb) = (vec![C64::ZERO; r * r], vec![C64::ZERO; r * r]);
    q.view().dots_block(r, &aq, r, &mut ha, None);
    q.view().dots_block(r, &bq, r, &mut hb, None);
    let da = Dense::from_fn(r, r, |i, j| ha[j * r + i]);
    let db = Dense::from_fn(r, r, |i, j| hb[j * r + i]);
    let (lam, vecs) = generalized_eig_c64(&da, &db)?;

    let mut accepted: Vec<(C64, Vec<C64>, f64)> = Vec::new();
    let mut rejected = 0;
    let (mut ax, mut bx) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    for (i, &lm) in lam.iter().enumerate() {
        if !lm.re.is_finite() || !opts.region.contains(lm) {
            continue;
        }
        let mut x = vec![C64::ZERO; n];
        q.view().add_combination(r, vecs.col(i), &mut x);
        let nx = norm2(&x);
        if nx == 0.0 {
            continue;
        }
        crate::linalg::vecops::scale_real(1.0 / nx, &mut x);
        ac.spmv(&x, &mut ax);
        match &bc {
            Some(bm) => bm.spmv(&x, &mut bx),
            None => bx.copy_from_slice(&x),
        }
        let res = ax.iter().zip(&bx).map(|(p, q2)| (*p - lm * *q2).norm_sqr()).sum::<f64>().sqrt()
            / lm.norm().max(f64::MIN_POSITIVE);
        if res <= opts.residual_tol {
            accepted.push((lm, x, res));
        } else {
            rejected += 1;
        }
    }
    accepted.sort_by(|x, y| (x.0.re, x.0.im).partial_cmp(&(y.0.re, y.0.im)).unwrap());
    let mut vectors = Basis::<C64>::new(n, accepted.len().max(1));
    for (i, (_, x, _)) in accepted.iter().enumerate() {
        vectors.col_mut(i).copy_from_slice(x);
    }
    tm.project = t.elapsed().as_secs_f64();
    tm.total = t_start.elapsed().as_secs_f64();
    Ok(ContourResult {
        values: accepted.iter().map(|x| x.0).collect(),
        residuals: accepted.iter().map(|x| x.2).collect(),
        vectors,
        estimated_count: est,
        count_stderr,
        rank: r,
        saturated: r == cols,
        rejected,
        timings: tm,
    })
}
