//! Nonlinear eigenvalue problems `T(λ) x = 0` in the split form `T(λ) = Σ_j f_j(λ) A_j` (port item
//! N2): dispersive materials (Drude/Lorentz permittivity), frequency-dependent PML stretching,
//! delay and radiation boundary conditions, and the polynomial case.
//!
//! Two methods, both reusing the sparse direct layer:
//!
//! * [`residual_inverse_iteration`] (Neumaier 1985): *one* eigenpair near a shift `σ` from a single
//!   factorisation of `T(σ)`: `x ← x − T(σ)⁻¹ T(λ) x` with `λ` from the scalar equation
//!   `xᴴ T(σ)⁻¹ T(λ) x = 0`. Linear convergence with rate proportional to the distance of `σ` from the
//!   eigenvalue; this is the cheap first step the port list names. Use it to polish a pair.
//! * [`beyn_contour`] (Beyn 2012; block Hankel form of Asakura–Sakurai et al.): *all* eigenvalues
//!   inside an ellipse from quadrature of `T(z)⁻¹ V` — one factorisation per node (values are replaced
//!   on one symbolic analysis), moments `Ŵᴴ S_p`, block Hankel matrices, a rank-revealing SVD and a small
//!   dense eigenproblem. The moment SVD gives the number of eigenvalues inside; a full-rank flag warns
//!   when `block·moments` is too small. Eigenpairs are filtered by region and by an independently
//!   recomputed residual.
//!
//! NLEIGS (rational Krylov on a rational interpolant, the robust method for many eigenvalues of large problems)
//! lives in [`crate::linalg::nleigs`]. **Not done:** nonlinear Rayleigh–Ritz refinement, see `docs/STATUS.md`.

use crate::linalg::basis::Basis;
use crate::linalg::contour::Ellipse;
use crate::linalg::csr::Csr;
use crate::linalg::dense::{Dense, DenseError};
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::linalg::faer_glue::{generalized_eig_c64, thin_svd_uv_c64};
use crate::linalg::rng::random_vector;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{dotc, norm2};

/// One term `f(λ)·A` of the split form.
pub struct NepTerm {
    pub matrix: Csr<C64>,
    pub f: Box<dyn Fn(C64) -> C64>,
}

/// `T(λ) = Σ f_j(λ) A_j`; all matrices must have the same dimensions.
pub struct NonlinearEigenProblem {
    pub terms: Vec<NepTerm>,
    n: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum NepError {
    Solver(SolverError),
    Dense(DenseError),
    InvalidArgument(&'static str),
}

impl std::fmt::Display for NepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NepError::Solver(e) => write!(f, "factorisation failed: {e}"),
            NepError::Dense(e) => write!(f, "dense step failed: {e}"),
            NepError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
        }
    }
}

impl std::error::Error for NepError {}

impl From<SolverError> for NepError {
    fn from(e: SolverError) -> Self {
        NepError::Solver(e)
    }
}
impl From<DenseError> for NepError {
    fn from(e: DenseError) -> Self {
        NepError::Dense(e)
    }
}

impl NonlinearEigenProblem {
    pub fn new(terms: Vec<NepTerm>) -> Result<Self, NepError> {
        let first = terms.first().ok_or(NepError::InvalidArgument("at least one term is required"))?;
        let n = first.matrix.nrows();
        if terms.iter().any(|t| t.matrix.nrows() != n || t.matrix.ncols() != n) {
            return Err(NepError::InvalidArgument("all terms must be square of the same size"));
        }
        Ok(Self { terms, n })
    }

    /// The polynomial problem `Σ_j λ^j A_j`.
    pub fn polynomial(coefficients: Vec<Csr<C64>>) -> Result<Self, NepError> {
        let terms = coefficients
            .into_iter()
            .enumerate()
            .map(|(j, matrix)| NepTerm { matrix, f: Box::new(move |z: C64| z.powi(j as i32)) })
            .collect();
        Self::new(terms)
    }

    pub fn dim(&self) -> usize {
        self.n
    }

    /// `T(z)` on the union pattern of the terms (explicit zeros are kept, so the pattern does not
    /// depend on `z`).
    pub fn assemble(&self, z: C64) -> Csr<C64> {
        let mut acc = Csr::axpby(C64::ZERO, &self.terms[0].matrix, (self.terms[0].f)(z), &self.terms[0].matrix);
        for t in &self.terms[1..] {
            acc = Csr::axpby(C64::ONE, &acc, (t.f)(z), &t.matrix);
        }
        acc
    }

    /// `T(λ) x` without assembling.
    pub fn apply(&self, lambda: C64, x: &[C64], y: &mut [C64]) {
        y.fill(C64::ZERO);
        let mut tmp = vec![C64::ZERO; self.n];
        for t in &self.terms {
            t.matrix.spmv(x, &mut tmp);
            let f = (t.f)(lambda);
            for (yi, ti) in y.iter_mut().zip(&tmp) {
                *yi += f * *ti;
            }
        }
    }

    /// `‖T(λ)x‖ / Σ_j |f_j(λ)|‖A_j x‖`: a backward-error-type relative residual in `[0, 1]`.
    pub fn relative_residual(&self, lambda: C64, x: &[C64]) -> f64 {
        let mut y = vec![C64::ZERO; self.n];
        self.apply(lambda, x, &mut y);
        let mut tmp = vec![C64::ZERO; self.n];
        let mut scale = 0.0;
        for t in &self.terms {
            t.matrix.spmv(x, &mut tmp);
            scale += (t.f)(lambda).norm() * norm2(&tmp);
        }
        norm2(&y) / scale.max(1e-300)
    }
}

// ---------------------------------------------------------------------------------------------
// Residual inverse iteration
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct RiiOptions {
    /// Relative residual (see [`NonlinearEigenProblem::relative_residual`]) to stop at.
    pub tol: f64,
    pub max_it: usize,
    /// Starting eigenvalue estimate (default: the shift).
    pub lambda0: Option<C64>,
    pub seed: u64,
}

impl Default for RiiOptions {
    fn default() -> Self {
        Self { tol: 1e-12, max_it: 100, lambda0: None, seed: 0x4E50 }
    }
}

#[derive(Debug, Clone)]
pub struct RiiResult {
    pub value: C64,
    pub vector: Vec<C64>,
    pub residual: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// One eigenpair near `sigma` by residual inverse iteration with a *fixed* factorisation of `T(σ)`.
pub fn residual_inverse_iteration(
    p: &NonlinearEigenProblem,
    sigma: C64,
    start: Option<&[C64]>,
    opts: &RiiOptions,
) -> Result<RiiResult, NepError> {
    let n = p.dim();
    let t_sigma = p.assemble(sigma);
    let mut fac = ShiftedFactorization::new(&t_sigma, None, Definiteness::General)?;
    fac.set_shift(C64::ZERO)?;
    let mut x = match start {
        Some(s) if s.len() == n => s.to_vec(),
        Some(_) => return Err(NepError::InvalidArgument("start vector has the wrong length")),
        None => random_vector::<C64>(n, opts.seed),
    };
    let nx = norm2(&x);
    if nx == 0.0 {
        return Err(NepError::InvalidArgument("start vector must be non-zero"));
    }
    x.iter_mut().for_each(|v| *v = v.scale(1.0 / nx));
    let mut lambda = opts.lambda0.unwrap_or(sigma);
    let mut tmp = vec![C64::ZERO; n];
    let mut iterations = 0;
    // The best iterate is returned: a shift sitting (numerically) on the eigenvalue makes T(σ) singular
    // and the iteration can lose accuracy, which must never make a good start worse.
    let mut best = (lambda, x.clone(), p.relative_residual(lambda, &x));
    if best.2 <= opts.tol {
        return Ok(RiiResult { value: best.0, vector: best.1, residual: best.2, iterations: 0, converged: true });
    }
    for it in 1..=opts.max_it {
        iterations = it;
        // w_j = T(σ)⁻¹ A_j x ; then g(λ) = Σ f_j(λ)(xᴴ w_j) = xᴴ T(σ)⁻¹ T(λ) x
        let mut w: Vec<Vec<C64>> = Vec::with_capacity(p.terms.len());
        let mut coef = Vec::with_capacity(p.terms.len());
        for t in &p.terms {
            t.matrix.spmv(&x, &mut tmp);
            fac.solve_block_in_place(&mut tmp, 1);
            coef.push(dotc(&x, &tmp));
            w.push(tmp.clone());
        }
        let g = |z: C64| -> C64 { p.terms.iter().zip(&coef).map(|(t, c)| (t.f)(z) * *c).sum() };
        // secant iteration for g(λ) = 0 from λ and a nearby point
        let (mut a, mut b) = (lambda, lambda + C64::new(1e-6, 1e-6) * (1.0 + lambda.norm()));
        let (mut ga, mut gb) = (g(a), g(b));
        for _ in 0..60 {
            if gb == ga {
                break;
            }
            let c = b - gb * (b - a) / (gb - ga);
            a = b;
            ga = gb;
            b = c;
            gb = g(b);
            if (b - a).norm() <= 1e-15 * b.norm().max(1.0) {
                break;
            }
        }
        lambda = b;
        // x ← x − T(σ)⁻¹ T(λ) x
        let mut next = x.clone();
        for (wj, t) in w.iter().zip(&p.terms) {
            let f = (t.f)(lambda);
            for (nv, wv) in next.iter_mut().zip(wj) {
                *nv -= f * *wv;
            }
        }
        let nn = norm2(&next);
        if nn == 0.0 || !nn.is_finite() {
            return Err(NepError::Dense(DenseError::NonFinite));
        }
        next.iter_mut().for_each(|v| *v = v.scale(1.0 / nn));
        x = next;
        let residual = p.relative_residual(lambda, &x);
        if residual < best.2 {
            best = (lambda, x.clone(), residual);
        }
        if residual <= opts.tol {
            break;
        }
    }
    Ok(RiiResult { value: best.0, vector: best.1, residual: best.2, iterations, converged: best.2 <= opts.tol })
}

// ---------------------------------------------------------------------------------------------
// Beyn's contour method (block Hankel form)
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct BeynOptions {
    pub region: Ellipse,
    /// Quadrature nodes (default 32).
    pub nodes: usize,
    /// Probing vectors `ℓ` (default 8): at least the largest multiplicity inside the region.
    pub block: usize,
    /// Moments `K` (default 3): `block·moments` should exceed the number of eigenvalues inside.
    pub moments: usize,
    /// Singular values of the Hankel matrix below `rank_tol·σ₁` are dropped (default 1e-10).
    pub rank_tol: f64,
    /// Pairs with relative residual above this are rejected (default 1e-6; the Hankel extraction is
    /// less accurate than Rayleigh–Ritz: polish accepted pairs with [`residual_inverse_iteration`]).
    pub residual_tol: f64,
    pub seed: u64,
}

impl BeynOptions {
    pub fn new(region: Ellipse) -> Self {
        Self { region, nodes: 32, block: 8, moments: 3, rank_tol: 1e-10, residual_tol: 1e-6, seed: 0xBE17 }
    }
}

#[derive(Debug, Clone)]
pub struct BeynResult {
    pub values: Vec<C64>,
    pub vectors: Basis<C64>,
    pub residuals: Vec<f64>,
    /// Numerical rank of the Hankel matrix: the number of eigenvalues inside (with multiplicity), up
    /// to quadrature leakage.
    pub rank: usize,
    /// The Hankel matrix has full rank `block·moments`: there may be more eigenvalues inside.
    pub saturated: bool,
    /// Candidates inside the region rejected for a large residual.
    pub rejected: usize,
}

pub fn beyn_contour(p: &NonlinearEigenProblem, opts: &BeynOptions) -> Result<BeynResult, NepError> {
    let n = p.dim();
    let (l, k) = (opts.block, opts.moments);
    if opts.nodes < 4 || l == 0 || k == 0 || k * l > n {
        return Err(NepError::InvalidArgument("need nodes ≥ 4, block ≥ 1, moments ≥ 1 and block·moments ≤ n"));
    }
    let rho = opts.region.a.max(opts.region.b);
    let center = opts.region.center;
    let v_hat: Vec<C64> = random_vector::<C64>(n * l, opts.seed);
    let w_hat: Vec<C64> = random_vector::<C64>(n * l, opts.seed ^ 0x9E37_79B9);
    // S_p = Σ_j w_j ζ_j^p T(z_j)⁻¹ V̂ for p = 0..2K−1 ; we keep S_0..S_{K−1} (for the vectors) and the
    // small projections M_p = Ŵᴴ S_p for all 2K moments.
    let nm = 2 * k;
    let mut s_keep: Vec<Vec<C64>> = vec![vec![C64::ZERO; n * l]; k];
    let mut m_small: Vec<Dense<C64>> = (0..nm).map(|_| Dense::zeros(l, l)).collect();
    let mut fac: Option<ShiftedFactorization<C64>> = None;
    for j in 0..opts.nodes {
        let (z, w) = opts.region.node(j, opts.nodes);
        let zeta = (z - center) / C64::new(rho, 0.0);
        let tz = p.assemble(z);
        let f = match fac.as_mut() {
            Some(f) => {
                f.replace_matrices(&tz, None)?;
                f
            }
            None => fac.insert(ShiftedFactorization::new(&tz, None, Definiteness::General)?),
        };
        f.set_shift(C64::ZERO)?;
        let mut y = v_hat.clone();
        f.solve_block_in_place(&mut y, l);
        let mut zp = C64::ONE;
        for pidx in 0..nm {
            let coef = w * zp;
            if pidx < k {
                for (s, yv) in s_keep[pidx].iter_mut().zip(&y) {
                    *s += coef * *yv;
                }
            }
            for a in 0..l {
                for b in 0..l {
                    let d = dotc(&w_hat[a * n..(a + 1) * n], &y[b * n..(b + 1) * n]);
                    m_small[pidx][(a, b)] += coef * d;
                }
            }
            zp *= zeta;
        }
    }
    // Block Hankel matrices H0 = [M_{a+b}], H1 = [M_{a+b+1}] of size Kℓ.
    let dim = k * l;
    let hank = |shift: usize| {
        Dense::from_fn(dim, dim, |i, j| {
            let (bi, bj) = (i / l, j / l);
            m_small[bi + bj + shift][(i % l, j % l)]
        })
    };
    let (h0, h1) = (hank(0), hank(1));
    let h0_cm: Vec<C64> = (0..dim * dim).map(|idx| h0[(idx % dim, idx / dim)]).collect();
    let (sv, u, v) = thin_svd_uv_c64(&h0_cm, dim, dim)?;
    let rank = sv.iter().take_while(|&&s| s > opts.rank_tol * sv[0]).count();
    let saturated = rank == dim;
    if rank == 0 {
        return Ok(BeynResult {
            values: vec![],
            vectors: Basis::new(n, 1),
            residuals: vec![],
            rank: 0,
            saturated: false,
            rejected: 0,
        });
    }
    // B = U_rᴴ H1 V_r Σ_r⁻¹
    let ucol = |c: usize| &u[c * dim..(c + 1) * dim];
    let vcol = |c: usize| &v[c * dim..(c + 1) * dim];
    let h1v: Vec<Vec<C64>> =
        (0..rank).map(|c| (0..dim).map(|i| (0..dim).map(|j| h1[(i, j)] * vcol(c)[j]).sum()).collect()).collect(); // H1 v_c
    let b = Dense::from_fn(rank, rank, |a, c| dotc(ucol(a), &h1v[c]).scale(1.0 / sv[c]));
    let (zetas, y) = generalized_eig_c64(&b, &Dense::identity(rank))?;
    // eigenvectors x_i = [S_0 … S_{K−1}] V_r Σ_r⁻¹ s_i
    let mut cand: Vec<(C64, Vec<C64>)> = Vec::new();
    for (i, zeta) in zetas.iter().enumerate() {
        if !zeta.re.is_finite() {
            continue;
        }
        let lam = center + zeta.scale(rho);
        // coefficient vector in the Kℓ space: V_r Σ⁻¹ s_i
        let mut coef = vec![C64::ZERO; dim];
        for c in 0..rank {
            let sc = y[(c, i)].scale(1.0 / sv[c]);
            for (r, cv) in coef.iter_mut().enumerate() {
                *cv += vcol(c)[r] * sc;
            }
        }
        let mut x = vec![C64::ZERO; n];
        for (bi, sk) in s_keep.iter().enumerate() {
            for a in 0..l {
                let c = coef[bi * l + a];
                for (xv, sv_) in x.iter_mut().zip(&sk[a * n..(a + 1) * n]) {
                    *xv += c * *sv_;
                }
            }
        }
        let nx = norm2(&x);
        if nx > 0.0 {
            x.iter_mut().for_each(|v| *v = v.scale(1.0 / nx));
            cand.push((lam, x));
        }
    }
    let mut values = Vec::new();
    let mut residuals = Vec::new();
    let mut kept: Vec<Vec<C64>> = Vec::new();
    let mut rejected = 0;
    for (lam, x) in cand {
        if !opts.region.contains(lam) {
            continue;
        }
        let r = p.relative_residual(lam, &x);
        if r <= opts.residual_tol {
            values.push(lam);
            residuals.push(r);
            kept.push(x);
        } else {
            rejected += 1;
        }
    }
    let mut vectors = Basis::<C64>::new(n, kept.len().max(1));
    for (j, x) in kept.iter().enumerate() {
        vectors.col_mut(j).copy_from_slice(x);
    }
    Ok(BeynResult { values, vectors, residuals, rank, saturated, rejected })
}
