//! Optimisation (port item O1; the role of PETSc `TAO`): **L-BFGS**, **nonlinear conjugate
//! gradients** (Polak–Ribière⁺) with a strong-Wolfe line search for smooth unconstrained problems,
//! and the **method of moving asymptotes** (MMA, Svanberg 1987) for box-constrained problems with
//! inequality constraints, as used in topology optimisation. Gradients are the caller's job:
//! [`adjoint_gradient`] shows the adjoint pattern for objectives that depend on the solution of a
//! linear system `K(p) u = f(p)`, and [`crate::linalg::sens`] gives eigenvalue derivatives for eigenvalue
//! objectives.
//!
//! MMA here solves the convex separable subproblem through its *dual*: `x_j(λ)` is closed form, the dual
//! variables `λ_i ∈ [0, c]` (the upper bound `c` is the exact-penalty form of Svanberg's artificial
//! variables, which keeps every subproblem feasible) are found by coordinate-wise bisection on the constraint
//! values, which is exact for one constraint and converges for a few. Large constraint counts would need
//! the primal–dual interior point method of the reference implementation (not done).

use crate::linalg::csr::Csr;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    GradientTolerance,
    FunctionTolerance,
    StepTolerance,
    MaxIterations,
    LineSearchFailed,
    KktTolerance,
}

impl StopReason {
    pub fn is_converged(self) -> bool {
        matches!(
            self,
            StopReason::GradientTolerance
                | StopReason::FunctionTolerance
                | StopReason::StepTolerance
                | StopReason::KktTolerance
        )
    }
}

#[derive(Debug, Clone)]
pub struct MinResult {
    pub x: Vec<f64>,
    pub f: f64,
    pub grad_norm: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub reason: StopReason,
    /// Objective after each iteration (index 0 = start).
    pub history: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct MinOptions {
    /// Stop when `‖g‖ ≤ gtol·max(1, |f|)`... (relative to `max(1, ‖x‖)`: `‖g‖ ≤ gtol·max(1, ‖x‖)`).
    pub gtol: f64,
    /// Stop when `|f_k − f_{k−1}| ≤ ftol·max(1, |f_k|)` for three consecutive iterations.
    pub ftol: f64,
    pub max_it: usize,
    /// L-BFGS memory.
    pub memory: usize,
    pub c1: f64,
    pub c2: f64,
}

impl Default for MinOptions {
    fn default() -> Self {
        Self { gtol: 1e-8, ftol: 0.0, max_it: 1000, memory: 10, c1: 1e-4, c2: 0.9 }
    }
}

type Fg<'a> = &'a mut dyn FnMut(&[f64], &mut [f64]) -> f64;
/// Line-search evaluation `α ↦ (f, f')` along the search direction (also fills the trial point and gradient).
type LineEval<'a> = &'a mut dyn FnMut(f64, &mut [f64], &mut [f64]) -> (f64, f64);

/// Strong-Wolfe line search (Nocedal & Wright, Algorithms 3.5/3.6) along `d` from `x`; returns the step
/// `α`, `f(x+αd)` and fills `x_new`, `g_new`.
#[allow(clippy::too_many_arguments)]
fn wolfe_search(
    fg: Fg,
    x: &[f64],
    d: &[f64],
    f0: f64,
    g0: &[f64],
    alpha_init: f64,
    c1: f64,
    c2: f64,
    x_new: &mut [f64],
    g_new: &mut [f64],
    evals: &mut usize,
) -> Option<(f64, f64)> {
    let dg0 = dot(g0, d);
    if dg0 >= 0.0 {
        return None;
    }
    let mut eval = |a: f64, xn: &mut [f64], gn: &mut [f64]| -> (f64, f64) {
        for ((o, &xi), &di) in xn.iter_mut().zip(x).zip(d) {
            *o = xi + a * di;
        }
        let f = fg(xn, gn);
        *evals += 1;
        (f, dot(gn, d))
    };
    let round = 1e-14 * f0.abs().max(1.0);
    let (mut a_prev, mut f_prev, mut dg_prev) = (0.0, f0, dg0);
    let mut a = alpha_init;
    for i in 0..40 {
        let (f, dg) = eval(a, x_new, g_new);
        if !f.is_finite() {
            a = 0.5 * (a_prev + a);
            continue;
        }
        // the sufficient-decrease test is relaxed by a rounding allowance, otherwise the search fails once the
        // expected decrease is below the resolution of `f` (flat minima, tight gradient tolerances)
        if f > f0 + c1 * a * dg0 + round || (i > 0 && f >= f_prev + round) {
            return zoom(&mut eval, a_prev, f_prev, dg_prev, a, f, dg, f0, dg0, c1, c2, round, x_new, g_new);
        }
        if dg.abs() <= -c2 * dg0 {
            return Some((a, f));
        }
        if dg >= 0.0 {
            return zoom(&mut eval, a, f, dg, a_prev, f_prev, dg_prev, f0, dg0, c1, c2, round, x_new, g_new);
        }
        a_prev = a;
        f_prev = f;
        dg_prev = dg;
        a *= 2.0;
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn zoom(
    eval: LineEval,
    mut lo: f64,
    mut flo: f64,
    mut dlo: f64,
    mut hi: f64,
    mut fhi: f64,
    mut dhi: f64,
    f0: f64,
    dg0: f64,
    c1: f64,
    c2: f64,
    round: f64,
    x_new: &mut [f64],
    g_new: &mut [f64],
) -> Option<(f64, f64)> {
    for _ in 0..40 {
        // cubic interpolation between lo and hi (safeguarded to the inner 80 %)
        let (d1, d2);
        let dd = hi - lo;
        d1 = dlo + dhi - 3.0 * (flo - fhi) / (lo - hi);
        let disc = d1 * d1 - dlo * dhi;
        let mut a = if disc >= 0.0 {
            d2 = disc.sqrt() * (hi - lo).signum();
            hi - dd * (dhi + d2 - d1) / (dhi - dlo + 2.0 * d2)
        } else {
            0.5 * (lo + hi)
        };
        let (mn, mx) = (lo.min(hi), lo.max(hi));
        if !a.is_finite() || a <= mn + 0.1 * (mx - mn) || a >= mx - 0.1 * (mx - mn) {
            a = 0.5 * (lo + hi);
        }
        let (f, dg) = eval(a, x_new, g_new);
        if f > f0 + c1 * a * dg0 + round || f >= flo + round {
            hi = a;
            fhi = f;
            dhi = dg;
        } else {
            if dg.abs() <= -c2 * dg0 {
                return Some((a, f));
            }
            if dg * (hi - lo) >= 0.0 {
                hi = lo;
                fhi = flo;
                dhi = dlo;
            }
            lo = a;
            flo = f;
            dlo = dg;
        }
        if (hi - lo).abs() < 1e-14 * lo.abs().max(1.0) {
            break;
        }
    }
    None
}

/// L-BFGS. `fg(x, g)` writes the gradient and returns the objective.
pub fn lbfgs(fg: Fg, x0: &[f64], opts: &MinOptions) -> MinResult {
    minimise(fg, x0, opts, Method::Lbfgs)
}

/// Nonlinear conjugate gradients, Polak–Ribière⁺ with restarts (use `opts.c2 = 0.1`).
pub fn nonlinear_cg(fg: Fg, x0: &[f64], opts: &MinOptions) -> MinResult {
    minimise(fg, x0, &MinOptions { c2: opts.c2.min(0.4), ..opts.clone() }, Method::Cg)
}

#[derive(Clone, Copy, PartialEq)]
enum Method {
    Lbfgs,
    Cg,
}

fn minimise(fg: Fg, x0: &[f64], opts: &MinOptions, method: Method) -> MinResult {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut g = vec![0.0; n];
    let mut f = fg(&x, &mut g);
    let mut evals = 1usize;
    let mut history = vec![f];
    let (mut xn, mut gn) = (vec![0.0; n], vec![0.0; n]);
    let mut d = vec![0.0; n];
    let mut s_hist: Vec<Vec<f64>> = Vec::new();
    let mut y_hist: Vec<Vec<f64>> = Vec::new();
    let mut small_f = 0;
    let (mut last_alpha, mut last_dg) = (f64::NAN, f64::NAN);
    let mut reason = StopReason::MaxIterations;
    let mut it = 0;
    while it < opts.max_it {
        if norm(&g) <= opts.gtol * norm(&x).max(1.0) {
            reason = StopReason::GradientTolerance;
            break;
        }
        // direction
        match method {
            Method::Lbfgs => {
                let mut q = g.clone();
                let m = s_hist.len();
                let mut alpha = vec![0.0; m];
                for i in (0..m).rev() {
                    let rho = 1.0 / dot(&y_hist[i], &s_hist[i]);
                    alpha[i] = rho * dot(&s_hist[i], &q);
                    for (qj, yj) in q.iter_mut().zip(&y_hist[i]) {
                        *qj -= alpha[i] * yj;
                    }
                }
                let gamma = if m > 0 {
                    dot(&s_hist[m - 1], &y_hist[m - 1]) / dot(&y_hist[m - 1], &y_hist[m - 1])
                } else {
                    1.0 / norm(&g).max(1.0)
                };
                for qj in q.iter_mut() {
                    *qj *= gamma;
                }
                for i in 0..m {
                    let rho = 1.0 / dot(&y_hist[i], &s_hist[i]);
                    let beta = rho * dot(&y_hist[i], &q);
                    for (qj, sj) in q.iter_mut().zip(&s_hist[i]) {
                        *qj += (alpha[i] - beta) * sj;
                    }
                }
                for (dj, qj) in d.iter_mut().zip(&q) {
                    *dj = -qj;
                }
            }
            Method::Cg => {
                if s_hist.is_empty() {
                    for (dj, gj) in d.iter_mut().zip(&g) {
                        *dj = -gj;
                    }
                } else {
                    // β = max(0, gnewᵀ(gnew − gold)/‖gold‖²); y_hist holds gnew − gold, s_hist the previous direction
                    let y = &y_hist[0];
                    let gold2 = dot(&g, &g) - 2.0 * dot(&g, y) + dot(y, y); // ‖gold‖² = ‖g − y‖²
                    let beta = (dot(&g, y) / gold2.max(1e-300)).max(0.0);
                    for ((dj, gj), pj) in d.iter_mut().zip(&g).zip(&s_hist[0]) {
                        *dj = -gj + beta * pj;
                    }
                    if dot(&d, &g) >= 0.0 {
                        for (dj, gj) in d.iter_mut().zip(&g) {
                            *dj = -gj;
                        }
                    }
                }
            }
        }
        let dg_now = dot(&g, &d);
        let alpha0 = if it == 0 && (method == Method::Cg || s_hist.is_empty()) {
            (1.0 / norm(&d).max(1e-300)).min(1.0)
        } else if method == Method::Cg && last_alpha.is_finite() && dg_now < 0.0 {
            // Nocedal–Wright (3.60): scale the previous step by the ratio of directional derivatives
            (last_alpha * last_dg / dg_now).clamp(1e-3, 1e3)
        } else {
            1.0
        };
        let Some((alpha, fnew)) =
            wolfe_search(&mut *fg, &x, &d, f, &g, alpha0, opts.c1, opts.c2, &mut xn, &mut gn, &mut evals)
        else {
            // a steepest-descent restart if a quasi-Newton/CG direction failed
            if !s_hist.is_empty() {
                s_hist.clear();
                y_hist.clear();
                continue;
            }
            reason = StopReason::LineSearchFailed;
            break;
        };
        last_alpha = alpha;
        last_dg = dg_now;
        let s: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&g).map(|(a, b)| a - b).collect();
        match method {
            Method::Lbfgs => {
                if dot(&s, &y) > 1e-12 * norm(&s) * norm(&y) {
                    s_hist.push(s.clone());
                    y_hist.push(y);
                    if s_hist.len() > opts.memory {
                        s_hist.remove(0);
                        y_hist.remove(0);
                    }
                }
            }
            Method::Cg => {
                s_hist = vec![d.iter().map(|v| v * alpha).collect()];
                s_hist[0] = d.clone();
                y_hist = vec![y];
            }
        }
        let df = (f - fnew).abs();
        x.copy_from_slice(&xn);
        g.copy_from_slice(&gn);
        f = fnew;
        history.push(f);
        it += 1;
        if opts.ftol > 0.0 && df <= opts.ftol * f.abs().max(1.0) {
            small_f += 1;
            if small_f >= 3 {
                reason = StopReason::FunctionTolerance;
                break;
            }
        } else {
            small_f = 0;
        }
    }
    MinResult { x, f, grad_norm: norm(&g), iterations: it, evaluations: evals, reason, history }
}

// ---------------------------------------------------------------------------------------------
// Adjoint gradient
// ---------------------------------------------------------------------------------------------

/// Adjoint gradient of `J(p) = gᵀ u` with `K(p) u = f` for real `K`:
/// `dJ/dp_k = (∂f/∂p_k)ᵀ w − wᵀ (∂K/∂p_k) u` with `K w = g` (symmetric `K`; pass `Kᵀ` otherwise).
/// `dk_du(k)` returns `(∂K/∂p_k) u` and `df(k)` returns `∂f/∂p_k` (or `None` if `f` does not depend on `p`).
/// One factorisation is shared by the forward and adjoint solves.
pub fn adjoint_gradient(
    k_mat: &Csr<f64>,
    f: &[f64],
    g: &[f64],
    n_params: usize,
    mut dk_u: impl FnMut(usize, &[f64]) -> Vec<f64>,
    mut df: impl FnMut(usize) -> Option<Vec<f64>>,
    definiteness: Definiteness,
) -> Result<(f64, Vec<f64>, Vec<f64>), SolverError> {
    let mut fac = ShiftedFactorization::new(k_mat, None, definiteness)?;
    fac.set_shift(0.0)?;
    let mut u = f.to_vec();
    fac.solve_block_in_place(&mut u, 1);
    let mut w = g.to_vec();
    fac.solve_block_in_place(&mut w, 1);
    let mut grad = vec![0.0; n_params];
    for (k, gk) in grad.iter_mut().enumerate() {
        let ku = dk_u(k, &u);
        *gk = -dot(&w, &ku);
        if let Some(dfk) = df(k) {
            *gk += dot(&dfk, &w);
        }
    }
    Ok((dot(g, &u), grad, u))
}

// ---------------------------------------------------------------------------------------------
// MMA
// ---------------------------------------------------------------------------------------------

/// Function values and gradients for MMA: objective `f0`, constraints `c_i(x) ≤ 0`.
#[derive(Debug, Clone)]
pub struct MmaEval {
    pub f0: f64,
    pub g0: Vec<f64>,
    pub c: Vec<f64>,
    pub gc: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct MmaOptions {
    pub xmin: Vec<f64>,
    pub xmax: Vec<f64>,
    /// Move limit as a fraction of `xmax − xmin` per iteration.
    pub move_limit: f64,
    pub max_it: usize,
    /// Stop when the largest design change is below this.
    pub xtol: f64,
    /// Upper bound of the dual variables (the penalty of violating a constraint).
    pub penalty: f64,
    pub asymptote_init: f64,
    pub asymptote_increase: f64,
    pub asymptote_decrease: f64,
}

impl MmaOptions {
    pub fn new(xmin: Vec<f64>, xmax: Vec<f64>) -> Self {
        Self {
            xmin,
            xmax,
            move_limit: 0.2,
            max_it: 200,
            xtol: 1e-6,
            penalty: 1e3,
            asymptote_init: 0.5,
            asymptote_increase: 1.2,
            asymptote_decrease: 0.7,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MmaResult {
    pub x: Vec<f64>,
    pub f0: f64,
    pub constraints: Vec<f64>,
    pub iterations: usize,
    pub history: Vec<f64>,
    pub reason: StopReason,
}

struct Subproblem<'a> {
    xl: &'a [f64],
    low: &'a [f64],
    upp: &'a [f64],
    alpha: &'a [f64],
    beta: &'a [f64],
    p: Vec<Vec<f64>>, // [i][j], i = 0 objective
    q: Vec<Vec<f64>>,
    r: Vec<f64>,
}

impl Subproblem<'_> {
    #[allow(clippy::needless_range_loop)]
    fn x_of(&self, lam: &[f64], x: &mut [f64]) {
        for j in 0..x.len() {
            let mut pj = self.p[0][j];
            let mut qj = self.q[0][j];
            for (i, l) in lam.iter().enumerate() {
                pj += l * self.p[i + 1][j];
                qj += l * self.q[i + 1][j];
            }
            let (sp, sq) = (pj.max(1e-30).sqrt(), qj.max(1e-30).sqrt());
            x[j] = ((sp * self.low[j] + sq * self.upp[j]) / (sp + sq)).clamp(self.alpha[j], self.beta[j]);
        }
    }
    /// Value of the approximated constraint `i + 1` at `x`.
    fn constraint(&self, i: usize, x: &[f64]) -> f64 {
        self.r[i + 1]
            + (0..x.len())
                .map(|j| self.p[i + 1][j] / (self.upp[j] - x[j]) + self.q[i + 1][j] / (x[j] - self.low[j]))
                .sum::<f64>()
    }
}

/// Method of moving asymptotes. `eval(x)` returns the objective, constraints and gradients.
pub fn mma(eval: &mut dyn FnMut(&[f64]) -> MmaEval, x0: &[f64], opts: &MmaOptions) -> MmaResult {
    let n = x0.len();
    let range: Vec<f64> = (0..n).map(|j| (opts.xmax[j] - opts.xmin[j]).max(1e-12)).collect();
    let mut x = x0.to_vec();
    let (mut x1, mut x2) = (vec![0.0; n], vec![0.0; n]);
    let (mut low, mut upp) = (vec![0.0; n], vec![0.0; n]);
    let mut history = Vec::new();
    let mut reason = StopReason::MaxIterations;
    let mut last = eval(&x);
    let m = last.c.len();
    let mut iterations = 0;
    for k in 1..=opts.max_it {
        history.push(last.f0);
        // asymptotes
        for j in 0..n {
            if k <= 2 {
                low[j] = x[j] - opts.asymptote_init * range[j];
                upp[j] = x[j] + opts.asymptote_init * range[j];
            } else {
                let zz = (x[j] - x1[j]) * (x1[j] - x2[j]);
                let factor = if zz > 0.0 {
                    opts.asymptote_increase
                } else if zz < 0.0 {
                    opts.asymptote_decrease
                } else {
                    1.0
                };
                low[j] = x[j] - factor * (x1[j] - low[j]);
                upp[j] = x[j] + factor * (upp[j] - x1[j]);
                low[j] = low[j].clamp(x[j] - 10.0 * range[j], x[j] - 0.01 * range[j]);
                upp[j] = upp[j].clamp(x[j] + 0.01 * range[j], x[j] + 10.0 * range[j]);
            }
        }
        let alpha: Vec<f64> = (0..n)
            .map(|j| opts.xmin[j].max(low[j] + 0.1 * (x[j] - low[j])).max(x[j] - opts.move_limit * range[j]))
            .collect();
        let beta: Vec<f64> = (0..n)
            .map(|j| opts.xmax[j].min(upp[j] - 0.1 * (upp[j] - x[j])).min(x[j] + opts.move_limit * range[j]))
            .collect();
        // convex separable approximations of f0 and the constraints
        let (mut p, mut q, mut r) = (vec![vec![0.0; n]; m + 1], vec![vec![0.0; n]; m + 1], vec![0.0; m + 1]);
        for i in 0..=m {
            let (g, fval) = if i == 0 { (&last.g0, last.f0) } else { (&last.gc[i - 1], last.c[i - 1]) };
            let raa = 1e-5;
            let mut rr = fval;
            for j in 0..n {
                let (up, dn) = (g[j].max(0.0), (-g[j]).max(0.0));
                p[i][j] = (upp[j] - x[j]).powi(2) * (1.001 * up + 0.001 * dn + raa / range[j]);
                q[i][j] = (x[j] - low[j]).powi(2) * (0.001 * up + 1.001 * dn + raa / range[j]);
                rr -= p[i][j] / (upp[j] - x[j]) + q[i][j] / (x[j] - low[j]);
            }
            r[i] = rr;
        }
        let sp = Subproblem { xl: &x, low: &low, upp: &upp, alpha: &alpha, beta: &beta, p, q, r };
        let _ = sp.xl;
        // dual: coordinate-wise bisection on λ_i ∈ [0, penalty] so that the approximated constraint is ≤ 0
        let mut lam = vec![0.0; m];
        let mut xs = x.clone();
        for _sweep in 0..50 {
            let mut changed = 0.0f64;
            for i in 0..m {
                let eval_c = |l: f64, lam: &mut Vec<f64>, xs: &mut Vec<f64>| -> f64 {
                    lam[i] = l;
                    sp.x_of(lam, xs);
                    sp.constraint(i, xs)
                };
                let old = lam[i];
                if eval_c(0.0, &mut lam, &mut xs) <= 0.0 {
                    lam[i] = 0.0;
                } else if eval_c(opts.penalty, &mut lam, &mut xs) >= 0.0 {
                    lam[i] = opts.penalty;
                } else {
                    let (mut lo, mut hi) = (0.0, opts.penalty);
                    for _ in 0..100 {
                        let mid = 0.5 * (lo + hi);
                        if eval_c(mid, &mut lam, &mut xs) > 0.0 {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                        if hi - lo <= 1e-14 * hi.max(1e-10) {
                            break;
                        }
                    }
                    lam[i] = 0.5 * (lo + hi);
                }
                changed = changed.max((lam[i] - old).abs());
            }
            if m <= 1 || changed < 1e-12 {
                break;
            }
        }
        sp.x_of(&lam, &mut xs);
        x2.copy_from_slice(&x1);
        x1.copy_from_slice(&x);
        let change = xs.iter().zip(&x).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        x.copy_from_slice(&xs);
        iterations = k;
        last = eval(&x);
        if change <= opts.xtol {
            reason = StopReason::StepTolerance;
            break;
        }
    }
    history.push(last.f0);
    MmaResult { x, f0: last.f0, constraints: last.c, iterations, history, reason }
}
