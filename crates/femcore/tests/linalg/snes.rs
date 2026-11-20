// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Nonlinear solvers on manufactured problems: `F(u) = L u + u³ − g` with `g` chosen so that a given
//! `u*` is the exact solution; the arctan system on which full Newton steps diverge; an analytic complex
//! cubic; and Bratu-type 2-D problems for Newton against JFNK with an AMG preconditioner.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::amg::{Amg, AmgOptions};
use femcore::linalg::csr::{laplacian_1d, laplacian_2d};
use femcore::linalg::direct::Definiteness;
use femcore::linalg::linop::Preconditioner;
use femcore::linalg::snes::{
    FixedPointOptions, JfnkOptions, LineSearch, NewtonOptions, Reason, anderson, jfnk, newton, picard,
};
use femcore::linalg::{C64, Csr};

fn sines(n: usize) -> Vec<f64> {
    (1..=n)
        .map(|j| {
            let t = j as f64 / (n as f64 + 1.0);
            (std::f64::consts::PI * t).sin() + 0.3 * (3.0 * std::f64::consts::PI * t).sin()
        })
        .collect()
}

/// `L u + u³ − g` and its Jacobian, with `g = L u* + u*³`.
struct Cubic {
    l: Csr<f64>,
    g: Vec<f64>,
    n: usize,
}

impl Cubic {
    fn new(n: usize, scale: f64) -> (Self, Vec<f64>) {
        let l = laplacian_1d(n).map(|v| v * scale);
        let ustar: Vec<f64> = sines(n).iter().map(|v| 3.0 * v).collect();
        let mut g = vec![0.0; n];
        l.spmv(&ustar, &mut g);
        for (gi, u) in g.iter_mut().zip(&ustar) {
            *gi += u * u * u;
        }
        (Self { l, g, n }, ustar)
    }
    fn f(&self, u: &[f64], out: &mut [f64]) {
        self.l.spmv(u, out);
        for i in 0..self.n {
            out[i] += u[i] * u[i] * u[i] - self.g[i];
        }
    }
    fn jac(&self, u: &[f64]) -> Csr<f64> {
        let d: Vec<_> = (0..self.n).map(|i| (i, i, 3.0 * u[i] * u[i])).collect();
        Csr::axpby(1.0, &self.l, 1.0, &Csr::from_triplets(self.n, self.n, &d).unwrap())
    }
}

fn max_err(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

#[test]
fn newton_converges_quadratically_to_the_manufactured_solution() {
    let (p, ustar) = Cubic::new(400, 1.0);
    let r = newton(
        &mut |u, o| p.f(u, o),
        &mut |u| p.jac(u),
        &vec![0.0; 400],
        &NewtonOptions { rtol: 1e-12, ..Default::default() },
    )
    .unwrap();
    assert!(r.converged(), "{:?}", r.reason);
    assert!(max_err(&r.x, &ustar) < 1e-9, "{}", max_err(&r.x, &ustar));
    assert!(r.iterations <= 8, "{} iterations", r.iterations);
    // quadratic convergence at the end: r_{k+1} ≲ C r_k² once r_k is small
    let h = &r.residual_norms;
    let k = h.len() - 3;
    assert!(h[k + 1] < 1e-2 * h[k].max(1e-30) || h[k + 1] < 1e-10, "{h:?}");
    assert_eq!(r.symbolic_analyses, 1, "values are replaced on one symbolic analysis");
}

#[test]
fn line_search_rescues_a_diverging_newton_iteration() {
    // F_i(x) = atan(x_i): full Newton steps from |x₀| ≳ 1.39 diverge.
    let n = 10;
    let mut f = |x: &[f64], o: &mut [f64]| x.iter().zip(o.iter_mut()).for_each(|(a, b)| *b = a.atan());
    let mut jac = |x: &[f64]| {
        let t: Vec<_> = x.iter().enumerate().map(|(i, a)| (i, i, 1.0 / (1.0 + a * a))).collect();
        Csr::from_triplets(n, n, &t).unwrap()
    };
    let x0 = vec![2.0; n];
    let plain = newton(
        &mut f,
        &mut jac,
        &x0,
        &NewtonOptions { line_search: LineSearch::None, max_it: 20, ..Default::default() },
    )
    .unwrap();
    assert!(!plain.converged(), "full steps should diverge, got {:?}", plain.reason);
    assert_ne!(plain.reason, Reason::ConvergedStep);
    let bt = newton(&mut f, &mut jac, &x0, &NewtonOptions::default()).unwrap();
    assert!(bt.converged(), "{:?}", bt.reason);
    assert!(bt.x.iter().all(|v| v.abs() < 1e-6));
    assert!(bt.backtracks > 0);
}

#[test]
fn jacobian_lag_trades_factorisations_for_iterations() {
    let (p, ustar) = Cubic::new(300, 1.0);
    let run = |lag: usize| {
        newton(
            &mut |u, o| p.f(u, o),
            &mut |u| p.jac(u),
            &vec![0.0; 300],
            &NewtonOptions { jacobian_lag: lag, ..Default::default() },
        )
        .unwrap()
    };
    let (every, lagged) = (run(1), run(3));
    assert!(every.converged() && lagged.converged());
    assert!(max_err(&lagged.x, &ustar) < 1e-6);
    assert!(
        lagged.factorizations < every.factorizations + lagged.iterations
            && lagged.jacobian_evals < every.jacobian_evals + 1
    );
    assert!(lagged.iterations >= every.iterations);
    eprintln!(
        "Newton its {} / {} factorisations; lag 3: {} / {}",
        every.iterations, every.factorizations, lagged.iterations, lagged.factorizations
    );
}

#[test]
fn jfnk_matches_newton_with_and_without_a_preconditioner() {
    let (p, ustar) = Cubic::new(200, 1.0);
    let x0: Vec<f64> = ustar.iter().map(|v| 0.6 * v).collect(); // the preconditioner is built here
    let reference =
        newton(&mut |u, o| p.f(u, o), &mut |u| p.jac(u), &x0, &NewtonOptions { rtol: 1e-12, ..Default::default() })
            .unwrap();
    let opts = JfnkOptions { rtol: 1e-10, ..Default::default() };
    let plain = jfnk(&mut |u, o| p.f(u, o), &x0, None, &opts);
    assert!(plain.converged(), "{:?}", plain.reason);
    assert!(max_err(&plain.x, &reference.x) < 1e-7);
    // LU of the Jacobian at the starting guess as preconditioner
    let j0 = p.jac(&x0);
    let mut fac = femcore::linalg::direct::ShiftedFactorization::new(&j0, None, Definiteness::General).unwrap();
    use femcore::linalg::direct::ShiftedSolver;
    fac.set_shift(0.0).unwrap();
    struct LuPc<'a>(&'a femcore::linalg::direct::ShiftedFactorization<f64>);
    impl Preconditioner<f64> for LuPc<'_> {
        fn apply(&self, r: &[f64], z: &mut [f64]) {
            z.copy_from_slice(r);
            self.0.solve_block_in_place(z, 1);
        }
    }
    let pc = jfnk(&mut |u, o| p.f(u, o), &x0, Some(&LuPc(&fac)), &opts);
    assert!(pc.converged(), "{:?}", pc.reason);
    assert!(max_err(&pc.x, &reference.x) < 1e-7);
    eprintln!(
        "JFNK: {} Newton / {} Krylov / {} F evaluations without PC; {} / {} / {} with LU(J(x0))",
        plain.iterations,
        plain.linear_iterations,
        plain.function_evals,
        pc.iterations,
        pc.linear_iterations,
        pc.function_evals
    );
    assert!(
        pc.linear_iterations * 2 < plain.linear_iterations,
        "{} vs {}",
        pc.linear_iterations,
        plain.linear_iterations
    );
}

#[test]
fn anderson_accelerates_picard() {
    // G(x) = (L + 4I)⁻¹ (g − x³ + 4x): a contraction (3x² < 4) whose rate is close to 1 for the smooth modes
    let (p, ustar) = Cubic::new(30, 1.0);
    let ustar: Vec<f64> = ustar.iter().map(|v| 0.7 * v).collect();
    let mut gvec = vec![0.0; 30];
    p.l.spmv(&ustar, &mut gvec);
    for (gi, u) in gvec.iter_mut().zip(&ustar) {
        *gi += u * u * u;
    }
    let p = Cubic { g: gvec, ..p };
    let shift = Csr::axpby(1.0, &p.l, 4.0, &Csr::identity(30));
    let mut fac = femcore::linalg::direct::ShiftedFactorization::new(&shift, None, Definiteness::PositiveDefinite).unwrap();
    use femcore::linalg::direct::ShiftedSolver;
    fac.set_shift(0.0).unwrap();
    let mut g = |x: &[f64], out: &mut [f64]| {
        for i in 0..30 {
            out[i] = p.g[i] - x[i] * x[i] * x[i] + 4.0 * x[i];
        }
        fac.solve_block_in_place(out, 1);
    };
    let x0 = vec![0.0; 30];
    let opts = FixedPointOptions { rtol: 1e-10, max_it: 3000, ..Default::default() };
    let pic = picard(&mut g, &x0, &opts);
    let and = anderson(&mut g, &x0, &FixedPointOptions { depth: 6, ..opts.clone() });
    eprintln!(
        "Picard {} iterations ({:?}), Anderson(6) {} ({:?})",
        pic.iterations, pic.reason, and.iterations, and.reason
    );
    assert!(and.converged(), "{:?}", and.reason);
    assert!(max_err(&and.x, &ustar) < 1e-6);
    // Picard is so slow here that it may not reach the tolerance within the cap
    assert!(!pic.converged() || and.iterations * 3 < pic.iterations, "{} vs {}", and.iterations, pic.iterations);
}

#[test]
fn complex_analytic_cubic() {
    // F(u) = (L − k²) u + (1 + 0.5i) u³ − g,  u* complex.
    let n = 150;
    let k2 = 0.4;
    let a: Csr<C64> = Csr::axpby(C64::ONE, &laplacian_1d(n).to_complex(), C64::new(-k2, 0.0), &Csr::identity(n));
    let kappa = C64::new(1.0, 0.5);
    let ustar: Vec<C64> = sines(n).iter().enumerate().map(|(i, v)| C64::new(*v, 0.2 * (i as f64 / n as f64))).collect();
    let mut g = vec![C64::ZERO; n];
    a.spmv(&ustar, &mut g);
    for i in 0..n {
        g[i] += kappa * ustar[i] * ustar[i] * ustar[i];
    }
    let mut f = |u: &[C64], o: &mut [C64]| {
        a.spmv(u, o);
        for i in 0..n {
            o[i] += kappa * u[i] * u[i] * u[i] - g[i];
        }
    };
    let mut jac = |u: &[C64]| {
        let d: Vec<_> = (0..n).map(|i| (i, i, kappa * C64::new(3.0, 0.0) * u[i] * u[i])).collect();
        Csr::axpby(C64::ONE, &a, C64::ONE, &Csr::from_triplets(n, n, &d).unwrap())
    };
    let x0: Vec<C64> = ustar.iter().map(|v| *v * 0.5 + C64::new(0.1, 0.0)).collect();
    let r = newton(&mut f, &mut jac, &x0, &NewtonOptions { rtol: 1e-12, ..Default::default() }).unwrap();
    assert!(r.converged(), "{:?}", r.reason);
    let err = r.x.iter().zip(&ustar).map(|(a, b)| (*a - *b).norm()).fold(0.0, f64::max);
    assert!(err < 1e-9, "{err}");
}

/// Bratu-type problem `−Δu − μ e^u = 0` on a grid (Dirichlet), unit spacing scaled by `h⁻²`.
struct Bratu {
    l: Csr<f64>,
    mu: f64,
    n: usize,
}

impl Bratu {
    fn new(nx: usize, mu: f64) -> Self {
        let h2 = 1.0 / ((nx + 1) as f64).powi(2);
        Self { l: laplacian_2d(nx, nx).map(|v| v / h2), mu, n: nx * nx }
    }
    fn f(&self, u: &[f64], o: &mut [f64]) {
        self.l.spmv(u, o);
        for i in 0..self.n {
            o[i] -= self.mu * u[i].exp();
        }
    }
    fn jac(&self, u: &[f64]) -> Csr<f64> {
        let d: Vec<_> = (0..self.n).map(|i| (i, i, -self.mu * u[i].exp())).collect();
        Csr::axpby(1.0, &self.l, 1.0, &Csr::from_triplets(self.n, self.n, &d).unwrap())
    }
}

#[test]
fn bratu_newton_and_jfnk_with_amg_agree() {
    let p = Bratu::new(80, 6.0); // below the turning point μ* ≈ 6.8
    let x0 = vec![0.0; p.n];
    let nr = newton(
        &mut |u, o| p.f(u, o),
        &mut |u| p.jac(u),
        &x0,
        &NewtonOptions { rtol: 1e-10, definiteness: Definiteness::PositiveDefinite, ..Default::default() },
    )
    .unwrap();
    assert!(nr.converged(), "{:?}", nr.reason);
    // AMG of the *linear* part −Δ as a fixed preconditioner (the exponential term is mild)
    let amg = Amg::new(&p.l, AmgOptions::default()).unwrap();
    struct P<'a>(&'a Amg<f64>);
    impl Preconditioner<f64> for P<'_> {
        fn apply(&self, r: &[f64], z: &mut [f64]) {
            self.0.apply(r, z)
        }
    }
    let jr = jfnk(&mut |u, o| p.f(u, o), &x0, Some(&P(&amg)), &JfnkOptions { rtol: 1e-10, ..Default::default() });
    assert!(jr.converged(), "{:?} after {} its", jr.reason, jr.iterations);
    assert!(
        max_err(&jr.x, &nr.x) < 1e-6 * nr.x.iter().fold(0.0f64, |m, v| m.max(v.abs())),
        "{}",
        max_err(&jr.x, &nr.x)
    );
    eprintln!(
        "Bratu 80²: Newton {} its / {} factorisations; JFNK+AMG {} its / {} Krylov / {} F",
        nr.iterations, nr.factorizations, jr.iterations, jr.linear_iterations, jr.function_evals
    );
}

#[test]
fn failures_are_reported_not_hidden() {
    let n = 5;
    let jac = |_: &[f64]| Csr::<f64>::identity(n);
    // F returns NaN
    let r = newton(&mut |_, o| o.fill(f64::NAN), &mut |x| jac(x), &vec![1.0; n], &NewtonOptions::default()).unwrap();
    assert_eq!(r.reason, Reason::DivergedFunctionNorm);
    // iteration cap
    let mut slow = |x: &[f64], o: &mut [f64]| x.iter().zip(o.iter_mut()).for_each(|(a, b)| *b = a.powi(3));
    let mut jslow = |x: &[f64]| {
        let t: Vec<_> = x.iter().enumerate().map(|(i, a)| (i, i, 3.0 * a * a)).collect();
        Csr::from_triplets(n, n, &t).unwrap()
    };
    let r =
        newton(&mut slow, &mut jslow, &vec![1.0; n], &NewtonOptions { max_it: 3, rtol: 1e-12, ..Default::default() })
            .unwrap();
    assert_eq!(r.reason, Reason::DivergedMaxIt);
    // singular Jacobian
    let r = newton(
        &mut |x, o| o.copy_from_slice(x),
        &mut |_| Csr::<f64>::from_triplets(n, n, &[]).unwrap(),
        &vec![1.0; n],
        &NewtonOptions::default(),
    );
    assert!(r.is_err() || !r.unwrap().converged());
    assert!(norm(&[3.0, 4.0]) == 5.0);
}
