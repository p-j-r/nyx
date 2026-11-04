//! **NLEIGS** (Güttel, Van Beeumen, Meerbergen & Michiels 2014): a rational Krylov method for the nonlinear
//! eigenproblem `T(λ) x = 0`, `T(λ) = Σ_t f_t(λ) A_t`, in a region `Σ` of the complex plane (port item N2, the
//! robust method after [`crate::linalg::nep`]'s Beyn contour integral and residual inverse iteration).
//!
//! # Idea
//! 1. **Rational interpolation.** Each scalar function `f_t` is interpolated on `Σ` by a rational function in a
//!    *Newton basis* with nodes `σ_k` in `Σ` and poles `ξ_k` near the singularities of the `f_t` (branch points,
//!    poles): `b_0 = 1`, `b_{k+1}(z) = b_k(z)(z − σ_k) / (β_{k+1}(1 − z/ξ_{k+1}))`, `f_t ≈ Σ_{k=0}^{d} c_{t,k} b_k`. Nodes
//!    and poles are chosen by the **Leja–Bagby** greedy rule on discretisations of `Σ` and of the singularity set
//!    `Ξ` (so the interpolant converges geometrically on `Σ` even for `√z`-type functions), the last pole is
//!    `∞`. `T` becomes `T̃(λ) = Σ_k b_k(λ) D_k`, `D_k = Σ_t c_{t,k} A_t`.
//! 2. **Linearisation.** `T̃(λ)x = 0` is equivalent to the linear pencil `(λ𝐀 − 𝐁) y = 0` of size `d·n` with
//!    `y = [b_0 x; …; b_{d−1} x]` (derived in the code; first block row `Σ D_k y_k + …`, then one row per basis
//!    recurrence). The nonlinear eigenvalues in `Σ` are eigenvalues of this pencil.
//! 3. **Shift-and-invert Krylov–Schur** on the pencil. The operator `(𝐁 − s𝐀)⁻¹𝐀` is applied *without forming
//!    the big matrices*: block substitution reduces the solve to **one sparse factorisation of `T̃(s)`** and
//!    a few products with the `A_t`. The existing [`krylov_schur`] is reused unchanged.
//! 4. The eigenvalues in `Σ` are returned with the true residual `‖T(λ)x‖` of the *nonlinear* problem
//!    ([`NonlinearEigenProblem::relative_residual`]) and, by default, polished by residual inverse iteration.
//!
//! # Compared with [`beyn_contour`](crate::linalg::nep::beyn_contour)
//! NLEIGS needs one factorisation, not one per contour node, and finds eigenvalues by Krylov convergence, so it
//! scales to many eigenvalues of large problems; the price is `d` times longer Krylov vectors (`d·n`, no TOAR-style
//! compression here: memory `ncv · d · n` — the compact representation of the paper is not implemented) and the need
//! to describe the singularities of the `f_t`.
//!
//! # Usage notes
//! * `region` is the ellipse `Σ`; give the singularities of the `f_t` that lie *outside* `Σ` in `singularities`
//!   (discretise a branch cut, e.g. with [`ray_samples`]); leave it empty for entire functions (then the rational
//!   interpolant degenerates to a Newton polynomial interpolant).
//! * The shift is the centre of `Σ` unless given. Eigenvalues converge fastest near the shift.
//! * The reported `interpolation_error` is the maximum relative deviation of the interpolants from the `f_t` on a
//!   set of test points of `Σ`; the achievable eigenvalue accuracy is of that order (polishing removes it).

use crate::linalg::contour::Ellipse;
use crate::linalg::csr::Csr;
use crate::linalg::dense::Dense;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::eps::{EpsOptions, Which, krylov_schur};
use crate::linalg::linop::LinearOperator;
use crate::linalg::nep::{NepError, NonlinearEigenProblem, RiiOptions, residual_inverse_iteration};
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::axpy;

/// Points `from + direction·r_k` with `r_k` growing geometrically from `r_min` to `r_max` (and `from` itself):
/// a discretisation of a branch cut `from + direction·[0, ∞)`. `direction` need not be normalised.
pub fn ray_samples(from: C64, direction: C64, r_min: f64, r_max: f64, count: usize) -> Vec<C64> {
    let dir = direction / direction.norm();
    let mut v = vec![from];
    for k in 0..count {
        let t = if count > 1 { k as f64 / (count - 1) as f64 } else { 0.0 };
        v.push(from + dir.scale(r_min * (r_max / r_min).powf(t)));
    }
    v
}

#[derive(Debug, Clone)]
pub struct NleigsOptions {
    /// The region in which eigenvalues are wanted.
    pub region: Ellipse,
    /// Discretisation of the singularity set of the `f_t` (outside `region`); empty = entire functions.
    pub singularities: Vec<C64>,
    /// Number of eigenpairs requested from the Krylov solver (those closest to the shift); the result lists the
    /// converged ones that lie in the region.
    pub nev: usize,
    /// Shift of the spectral transformation (default: the centre of the region).
    pub shift: Option<C64>,
    /// Target accuracy of the rational interpolation; the degree is increased until it is met.
    pub interpolation_tol: f64,
    /// Largest interpolation degree `d` tried.
    pub max_degree: usize,
    /// Krylov–Schur tolerance (on the linearised problem).
    pub tol: f64,
    pub ncv: Option<usize>,
    pub max_restarts: usize,
    /// Boundary points of `region` used to discretise `Σ` for the Leja–Bagby choice (default 400).
    pub boundary_samples: usize,
    /// Polish accepted eigenpairs by residual inverse iteration on the nonlinear problem (one factorisation each).
    pub polish: bool,
    pub seed: u64,
}

impl NleigsOptions {
    pub fn new(region: Ellipse, nev: usize) -> Self {
        Self {
            region,
            singularities: Vec::new(),
            nev,
            shift: None,
            interpolation_tol: 1e-11,
            max_degree: 80,
            tol: 1e-10,
            ncv: None,
            max_restarts: 300,
            boundary_samples: 400,
            polish: true,
            seed: 0x4E4C,
        }
    }
}

#[derive(Debug, Clone)]
pub struct NleigsResult {
    /// Eigenvalues inside the region, ordered by distance to the shift.
    pub values: Vec<C64>,
    /// Unit-norm eigenvectors of the nonlinear problem.
    pub vectors: Vec<Vec<C64>>,
    /// `‖T(λ)x‖ / Σ_t |f_t(λ)|‖A_t x‖` (after polishing, if enabled).
    pub residuals: Vec<f64>,
    /// Degree `d` of the rational interpolant.
    pub degree: usize,
    pub interpolation_error: f64,
    /// Converged Krylov pairs, including those outside the region.
    pub converged_total: usize,
    pub matvecs: usize,
    pub restarts: usize,
}

/// The rational Newton interpolants of the scalar functions of a [`NonlinearEigenProblem`].
#[derive(Debug, Clone)]
pub struct RationalInterpolant {
    /// All nodes, poles and arguments are measured from this point (the centre of the region), so that a singularity
    /// at the origin of the original variable is representable (the basis divides by `1 − z/ξ`).
    pub center: C64,
    /// Nodes `σ_0 … σ_d` (relative to `center`).
    pub sigma: Vec<C64>,
    /// `1/ξ_k` for `k = 1..=d` (index 0 unused, relative to `center`); zero for a pole at infinity.
    pub inv_xi: Vec<C64>,
    /// Scalings `β_k` for `k = 1..=d` (index 0 unused).
    pub beta: Vec<f64>,
    /// `coef[t][k]`: the `k`-th rational divided difference of `f_t`.
    pub coef: Vec<Vec<C64>>,
}

impl RationalInterpolant {
    pub fn degree(&self) -> usize {
        self.sigma.len() - 1
    }

    /// `b_0(z), …, b_d(z)`.
    pub fn basis(&self, z: C64) -> Vec<C64> {
        let z = z - self.center;
        let d = self.degree();
        let mut b = Vec::with_capacity(d + 1);
        b.push(C64::ONE);
        for k in 0..d {
            let den = (C64::ONE - z * self.inv_xi[k + 1]).scale(self.beta[k + 1]);
            let prev = b[k];
            b.push(prev * (z - self.sigma[k]) / den);
        }
        b
    }

    /// The interpolant of term `t` at `z`.
    pub fn eval(&self, t: usize, z: C64) -> C64 {
        self.basis(z).iter().zip(&self.coef[t]).fold(C64::ZERO, |s, (b, c)| s + *b * *c)
    }
}

/// Leja–Bagby nodes/poles and divided differences for degree `d`.
fn leja_bagby(
    center: C64,
    sigma_hat: &[C64],
    xi_hat: &[C64],
    d: usize,
    funcs: &[&dyn Fn(C64) -> C64],
) -> Option<RationalInterpolant> {
    let sigma_hat: Vec<C64> = sigma_hat.iter().map(|z| *z - center).collect();
    let xi_hat: Vec<C64> = xi_hat.iter().map(|z| *z - center).collect();
    let mut b_s = vec![C64::ONE; sigma_hat.len()];
    let mut b_x = vec![C64::ONE; xi_hat.len()];
    let argmax = |v: &[C64]| {
        v.iter().enumerate().fold((0, -1.0f64), |m, (i, x)| if x.norm() > m.1 { (i, x.norm()) } else { m }).0
    };
    let mut sigma = vec![sigma_hat[0]];
    let mut inv_xi = vec![C64::ZERO];
    let mut beta = vec![1.0];
    for j in 0..d {
        let inv = if j + 1 == d || xi_hat.is_empty() {
            C64::ZERO
        } else {
            let k = b_x
                .iter()
                .enumerate()
                .fold((0, f64::INFINITY), |m, (i, x)| if x.norm() < m.1 { (i, x.norm()) } else { m })
                .0;
            C64::ONE / xi_hat[k]
        };
        let ratio = |z: C64| (z - sigma[j]) / (C64::ONE - z * inv);
        let r_s: Vec<C64> = sigma_hat.iter().zip(&b_s).map(|(z, b)| *b * ratio(*z)).collect();
        let bt = r_s.iter().map(|x| x.norm()).fold(0.0, f64::max);
        if !(bt.is_finite() && bt > 0.0) {
            return None;
        }
        b_s = r_s.into_iter().map(|x| x.scale(1.0 / bt)).collect();
        for (z, b) in xi_hat.iter().zip(b_x.iter_mut()) {
            let den = C64::ONE - *z * inv;
            *b = if den.norm() < 1e-300 {
                C64::new(f64::INFINITY, 0.0)
            } else {
                (*b * (*z - sigma[j]) / den).scale(1.0 / bt)
            };
        }
        inv_xi.push(inv);
        beta.push(bt);
        sigma.push(sigma_hat[argmax(&b_s)]);
    }
    let mut it = RationalInterpolant { center, sigma, inv_xi, beta, coef: Vec::new() };
    // lower-triangular system  f(σ_i) = Σ_{k ≤ i} c_k b_k(σ_i)
    let basis_at_nodes: Vec<Vec<C64>> = it.sigma.iter().map(|s| it.basis(*s + it.center)).collect();
    for f in funcs {
        let mut c = vec![C64::ZERO; d + 1];
        for i in 0..=d {
            let mut r = f(it.sigma[i] + it.center);
            for k in 0..i {
                r -= c[k] * basis_at_nodes[i][k];
            }
            let diag = basis_at_nodes[i][i];
            if diag.norm() == 0.0 {
                return None;
            }
            c[i] = r / diag;
        }
        if !c.iter().all(|x| x.is_finite()) {
            return None;
        }
        it.coef.push(c);
    }
    Some(it)
}

/// Builds the interpolant, raising the degree until the test error is below `tol`; returns it with that error.
pub fn build_interpolant(
    p: &NonlinearEigenProblem,
    region: &Ellipse,
    singularities: &[C64],
    tol: f64,
    max_degree: usize,
    boundary_samples: usize,
) -> Result<(RationalInterpolant, f64), NepError> {
    let m = boundary_samples.max(32);
    let sigma_hat: Vec<C64> = (0..m)
        .map(|j| {
            let th = 2.0 * std::f64::consts::PI * j as f64 / m as f64;
            region.center + C64::new(region.a * th.cos(), region.b * th.sin())
        })
        .collect();
    // test points: boundary points between the nodes of the discretisation, and interior points
    let mut tests: Vec<C64> = (0..m)
        .map(|j| {
            let th = 2.0 * std::f64::consts::PI * (j as f64 + 0.5) / m as f64;
            region.center + C64::new(region.a * th.cos(), region.b * th.sin())
        })
        .collect();
    for (r, k) in [(0.0, 1usize), (0.5, 16), (0.85, 24)] {
        for q in 0..k {
            let th = 2.0 * std::f64::consts::PI * (q as f64 + 0.3) / k as f64;
            tests.push(region.center + C64::new(r * region.a * th.cos(), r * region.b * th.sin()));
        }
    }
    let funcs: Vec<&dyn Fn(C64) -> C64> = p.terms.iter().map(|t| &*t.f as &dyn Fn(C64) -> C64).collect();
    let scale: Vec<f64> =
        funcs.iter().map(|f| tests.iter().map(|z| f(*z).norm()).fold(0.0, f64::max).max(1e-300)).collect();
    let mut degrees: Vec<usize> =
        [4usize, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128].into_iter().filter(|&d| d < max_degree).collect();
    degrees.push(max_degree.max(3));
    let mut best: Option<(RationalInterpolant, f64)> = None;
    for d in degrees {
        let Some(it) = leja_bagby(region.center, &sigma_hat, singularities, d, &funcs) else { continue };
        let err = funcs
            .iter()
            .enumerate()
            .map(|(t, f)| tests.iter().map(|z| (it.eval(t, *z) - f(*z)).norm()).fold(0.0, f64::max) / scale[t])
            .fold(0.0, f64::max);
        let better = best.as_ref().is_none_or(|b| err < b.1);
        if better {
            best = Some((it, err));
        }
        if err <= tol {
            break;
        }
    }
    best.ok_or(NepError::InvalidArgument("rational interpolation failed (degenerate nodes)"))
}

/// The operator `(𝐁 − s𝐀)⁻¹ 𝐀` of the NLEIGS linearisation, applied through one factorisation of `T̃(s)`.
pub struct NleigsOperator<'a> {
    p: &'a NonlinearEigenProblem,
    it: RationalInterpolant,
    /// The shift relative to the interpolant's centre.
    s: C64,
    n: usize,
    d: usize,
    /// `b_j(s)`.
    bs: Vec<C64>,
    fac: ShiftedFactorization<C64>,
}

impl<'a> NleigsOperator<'a> {
    pub fn new(p: &'a NonlinearEigenProblem, it: RationalInterpolant, s: C64) -> Result<Self, NepError> {
        let d = it.degree();
        if d < 2 {
            return Err(NepError::InvalidArgument("interpolation degree must be at least 2"));
        }
        let bs = it.basis(s);
        // T̃(s) = Σ_t (Σ_k b_k(s) c_{t,k}) A_t
        let mut acc: Option<Csr<C64>> = None;
        for (t, term) in p.terms.iter().enumerate() {
            let ft = bs.iter().zip(&it.coef[t]).fold(C64::ZERO, |a, (b, c)| a + *b * *c);
            acc = Some(match acc {
                None => Csr::axpby(C64::ZERO, &term.matrix, ft, &term.matrix),
                Some(a) => Csr::axpby(C64::ONE, &a, ft, &term.matrix),
            });
        }
        let ts = acc.expect("at least one term");
        let mut fac = ShiftedFactorization::new(&ts, None, Definiteness::General)?;
        fac.set_shift(C64::ZERO)?;
        let s = s - it.center;
        Ok(Self { p, it, s, n: p.dim(), d, bs, fac })
    }

    pub fn interpolant(&self) -> &RationalInterpolant {
        &self.it
    }

    /// The dense pencil `(𝐀, 𝐁)` of size `d·n` — for tests and small problems only.
    pub fn dense_pencil(&self) -> (Dense<C64>, Dense<C64>) {
        let (n, d) = (self.n, self.d);
        let it = &self.it;
        let dk = |k: usize| -> Dense<C64> {
            let mut m = Dense::<C64>::zeros(n, n);
            for (t, term) in self.p.terms.iter().enumerate() {
                let dense = term.matrix.to_dense();
                for i in 0..n {
                    for j in 0..n {
                        m[(i, j)] += it.coef[t][k] * dense[i * n + j];
                    }
                }
            }
            m
        };
        let mut a = Dense::<C64>::zeros(d * n, d * n);
        let mut b = Dense::<C64>::zeros(d * n, d * n);
        let put = |m: &mut Dense<C64>, bi: usize, bj: usize, blk: &Dense<C64>, f: C64| {
            for i in 0..n {
                for j in 0..n {
                    m[(bi * n + i, bj * n + j)] += f * blk[(i, j)];
                }
            }
        };
        let eye = Dense::<C64>::identity(n);
        let dd = dk(d);
        put(&mut a, 0, d - 1, &dd, C64::new(1.0 / it.beta[d], 0.0));
        for j in 0..d - 1 {
            put(&mut b, 0, j, &dk(j), -C64::ONE);
        }
        put(&mut b, 0, d - 1, &dk(d - 1), -C64::ONE);
        put(&mut b, 0, d - 1, &dd, it.sigma[d - 1].scale(1.0 / it.beta[d]));
        for j in 0..d - 1 {
            put(&mut a, j + 1, j, &eye, C64::ONE);
            put(&mut a, j + 1, j + 1, &eye, it.inv_xi[j + 1].scale(it.beta[j + 1]));
            put(&mut b, j + 1, j, &eye, it.sigma[j]);
            put(&mut b, j + 1, j + 1, &eye, C64::new(it.beta[j + 1], 0.0));
        }
        (a, b)
    }
}

impl LinearOperator<C64> for NleigsOperator<'_> {
    fn dim(&self) -> usize {
        self.d * self.n
    }

    fn apply(&self, v: &[C64], y: &mut [C64]) {
        let (n, d, s) = (self.n, self.d, self.s);
        let it = &self.it;
        let nt = self.p.terms.len();
        fn blk_of(x: &[C64], j: usize, n: usize) -> &[C64] {
            &x[j * n..(j + 1) * n]
        }
        // w = 𝐀 v
        let mut w = vec![C64::ZERO; d * n];
        {
            // w_0 = (1/β_d) Σ_t c_{t,d} A_t v_{d−1}
            let mut tmp = vec![C64::ZERO; n];
            for t in 0..nt {
                self.p.terms[t].matrix.spmv(blk_of(v, d - 1, n), &mut tmp);
                axpy(it.coef[t][d].scale(1.0 / it.beta[d]), &tmp, &mut w[..n]);
            }
        }
        for j in 0..d - 1 {
            let f = it.inv_xi[j + 1].scale(it.beta[j + 1]);
            let wj1 = &mut w[(j + 1) * n..(j + 2) * n];
            wj1.copy_from_slice(blk_of(v, j, n));
            axpy(f, blk_of(v, j + 1, n), wj1);
        }
        // c_0 = 0,  c_{j+1} = (w_{j+1} − (σ_j − s) c_j) / (β_{j+1}(1 − s/ξ_{j+1}))
        let mut c = vec![C64::ZERO; d * n];
        for j in 0..d - 1 {
            let den = (C64::ONE - s * it.inv_xi[j + 1]).scale(it.beta[j + 1]);
            let (prev, next) = c.split_at_mut((j + 1) * n);
            let cj = &prev[j * n..];
            let cn = &mut next[..n];
            let g = (it.sigma[j] - s) / den;
            for i in 0..n {
                cn[i] = w[(j + 1) * n + i] / den - g * cj[i];
            }
        }
        // rhs = −w_0 − Σ_t A_t u_t,  u_t = Σ_{j≤d−2} c_{t,j} c_j + (c_{t,d−1} + c_{t,d}(s−σ_{d−1})/β_d) c_{d−1}
        let mut rhs: Vec<C64> = w[..n].iter().map(|x| -*x).collect();
        let mut u = vec![C64::ZERO; n];
        let mut tmp = vec![C64::ZERO; n];
        for t in 0..nt {
            u.fill(C64::ZERO);
            for j in 1..d - 1 {
                axpy(it.coef[t][j], blk_of(&c, j, n), &mut u);
            }
            let last = it.coef[t][d - 1] + it.coef[t][d] * (s - it.sigma[d - 1]).scale(1.0 / it.beta[d]);
            axpy(last, blk_of(&c, d - 1, n), &mut u);
            self.p.terms[t].matrix.spmv(&u, &mut tmp);
            axpy(-C64::ONE, &tmp, &mut rhs);
        }
        self.fac.solve_block_in_place(&mut rhs, 1);
        // z_j = b_j(s) z_0 + c_j
        for j in 0..d {
            let yj = &mut y[j * n..(j + 1) * n];
            for i in 0..n {
                yj[i] = self.bs[j] * rhs[i] + c[j * n + i];
            }
        }
    }
}

/// Eigenvalues of `p` inside `opts.region` (see the module documentation).
pub fn nleigs(p: &NonlinearEigenProblem, opts: &NleigsOptions) -> Result<NleigsResult, NepError> {
    let (it, interp_err) = build_interpolant(
        p,
        &opts.region,
        &opts.singularities,
        opts.interpolation_tol,
        opts.max_degree,
        opts.boundary_samples,
    )?;
    let d = it.degree();
    let s = opts.shift.unwrap_or(opts.region.center);
    let op = NleigsOperator::new(p, it, s)?;
    let n = p.dim();
    let big = d * n;
    let nev = opts.nev.min(big.saturating_sub(1)).max(1);
    let eps = EpsOptions {
        tol: opts.tol,
        ncv: opts.ncv,
        max_it: Some(opts.max_restarts),
        which: Which::LargestMagnitude,
        seed: opts.seed,
        ..EpsOptions::general(nev)
    };
    let res = krylov_schur::<C64, _>(&op, None, None, &eps)
        .map_err(|_| NepError::InvalidArgument("Krylov solver rejected its input"))?;
    let mut found: Vec<(C64, Vec<C64>, f64)> = Vec::new();
    for (i, theta) in res.values.iter().enumerate().take(res.nconv.min(res.values.len())) {
        if theta.norm() == 0.0 {
            continue;
        }
        let lambda = s + C64::ONE / *theta;
        if !opts.region.contains(lambda) {
            continue;
        }
        let y = res.vector(i);
        let mut x = y[..n].to_vec();
        let nx = crate::linalg::vecops::norm2(&x);
        if nx == 0.0 {
            continue;
        }
        x.iter_mut().for_each(|v| *v = v.scale(1.0 / nx));
        let mut r = p.relative_residual(lambda, &x);
        let (mut lam, mut xv) = (lambda, x);
        if opts.polish {
            let ropts = RiiOptions { tol: 1e-13, max_it: 30, lambda0: Some(lam), ..Default::default() };
            let polished = residual_inverse_iteration(p, lam, Some(&xv), &ropts).ok().filter(|pol| {
                pol.residual < r
                    && opts.region.contains(pol.value)
                    && (pol.value - lam).norm() < 1e-3 * (1.0 + lam.norm())
            });
            if let Some(pol) = polished {
                lam = pol.value;
                xv = pol.vector;
                r = pol.residual;
            }
        }
        found.push((lam, xv, r));
    }
    found.sort_by(|a, b| (a.0 - s).norm().total_cmp(&(b.0 - s).norm()));
    Ok(NleigsResult {
        values: found.iter().map(|f| f.0).collect(),
        residuals: found.iter().map(|f| f.2).collect(),
        vectors: found.into_iter().map(|f| f.1).collect(),
        degree: d,
        interpolation_error: interp_err,
        converged_total: res.nconv,
        matvecs: res.matvecs,
        restarts: res.iterations,
    })
}
