//! Tier C checks for the time steppers: closed-form solutions on tridiagonal systems whose
//! eigenvectors are sines, convergence orders at fixed step, error-versus-tolerance for the adaptive
//! methods, and the algebraic properties (energy conservation, high-frequency damping).
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::laplacian_1d;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::ts::bdf::{BdfOptions, bdf_solve};
use femcore::linalg::ts::newmark::{NewmarkParams, newmark_solve};
use femcore::linalg::ts::rosenbrock::{RosenbrockOptions, rosenbrock_solve};
use femcore::linalg::ts::theta::theta_solve;
use femcore::linalg::{C64, Csr};
use std::f64::consts::PI;

fn sine(n: usize, k: usize) -> Vec<f64> {
    (1..=n).map(|j| (k as f64 * PI * j as f64 / (n as f64 + 1.0)).sin()).collect()
}

fn max_err(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

/// Heat problem `y' + L y = 0`, `y0 = s_a + s_b` (a slow and a stiff mode): exact solution.
fn heat_exact(n: usize, (ka, kb): (usize, usize), t: f64) -> Vec<f64> {
    let ev = laplacian_1d_eigenvalues(n);
    let (sa, sb) = (sine(n, ka), sine(n, kb));
    (0..n).map(|i| (-ev[ka - 1] * t).exp() * sa[i] + (-ev[kb - 1] * t).exp() * sb[i]).collect()
}

#[test]
fn bdf_error_decreases_with_the_tolerance_and_reuses_factorisations() {
    let n = 200;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let modes = (40, n);
    let y0 = heat_exact(n, modes, 0.0);
    let t_eval: Vec<f64> = (1..=10).map(|i| 0.5 * i as f64).collect();
    let mut errs = Vec::new();
    for rtol in [1e-3, 1e-5, 1e-7] {
        let sol =
            bdf_solve(&m, &k, &mut |_, f| f.fill(0.0), 0.0, 5.0, &y0, &t_eval, &BdfOptions::new(rtol, 1e-10)).unwrap();
        let e = sol.t.iter().zip(&sol.y).map(|(&t, y)| max_err(y, &heat_exact(n, modes, t))).fold(0.0, f64::max);
        eprintln!(
            "rtol {rtol:e}: error {e:.2e}, steps {}, rejected {}, factorisations {}, by order {:?}",
            sol.stats.steps, sol.stats.rejected, sol.stats.factorizations, sol.stats.steps_per_order
        );
        assert_eq!(sol.t.len(), t_eval.len());
        assert!(sol.stats.symbolic_analyses <= 2, "symbolic analyses {}", sol.stats.symbolic_analyses);
        assert!(sol.stats.factorizations * 3 <= sol.stats.steps + 3 * sol.stats.rejected + 10, "{:?}", sol.stats);
        errs.push(e);
    }
    assert!(errs[0] > errs[1] && errs[1] > errs[2], "{errs:?}");
    assert!(errs[0] < 50.0 * 1e-3 && errs[1] < 50.0 * 1e-5 && errs[2] < 50.0 * 1e-7, "{errs:?}");
}

#[test]
fn bdf_forced_response_matches_the_scalar_ode() {
    // y' + L y = sin(ωt) s_a :  c' + λ c = sin ωt, c(0) = 0.
    let n = 100;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let (ka, omega) = (7usize, 3.0);
    let lam = laplacian_1d_eigenvalues(n)[ka - 1];
    let sa = sine(n, ka);
    let mut force = |t: f64, f: &mut [f64]| {
        for (fi, s) in f.iter_mut().zip(&sa) {
            *fi = (omega * t).sin() * s;
        }
    };
    let t_end = 6.0;
    let sol =
        bdf_solve(&m, &k, &mut force, 0.0, t_end, &vec![0.0; n], &[t_end], &BdfOptions::new(1e-7, 1e-10)).unwrap();
    let c = (lam * (omega * t_end).sin() - omega * (omega * t_end).cos() + omega * (-lam * t_end).exp())
        / (lam * lam + omega * omega);
    let want: Vec<f64> = sa.iter().map(|s| c * s).collect();
    let e = max_err(&sol.y[0], &want);
    assert!(e < 1e-5, "error {e}");
}

#[test]
fn bdf_complex_symmetric_matrix() {
    // K = L + iγ I: exact decay e^{-(λ+iγ)t} of the sine mode.
    let n = 120;
    let gamma = 0.7;
    let k: Csr<C64> = Csr::axpby(C64::ONE, &laplacian_1d(n).to_complex(), C64::new(0.0, gamma), &Csr::identity(n));
    let m: Csr<C64> = Csr::identity(n);
    let ka = 20usize;
    let lam = laplacian_1d_eigenvalues(n)[ka - 1];
    let y0: Vec<C64> = sine(n, ka).iter().map(|&s| C64::new(s, 0.0)).collect();
    let t_end = 3.0;
    let sol =
        bdf_solve(&m, &k, &mut |_, f| f.fill(C64::ZERO), 0.0, t_end, &y0, &[t_end], &BdfOptions::new(1e-7, 1e-10))
            .unwrap();
    let decay = C64::new(-lam * t_end, -gamma * t_end);
    let factor = C64::new(decay.re.exp() * decay.im.cos(), decay.re.exp() * decay.im.sin());
    let e = y0.iter().zip(&sol.y[0]).map(|(a, b)| (*a * factor - *b).norm()).fold(0.0, f64::max);
    assert!(e < 1e-5, "error {e}");
}

#[test]
fn bdf_fixed_step_order_of_convergence() {
    let n = 100;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let modes = (5, 12);
    let y0 = heat_exact(n, modes, 0.0);
    let t_end = 1.0;
    let err_at = |steps: usize, max_order: usize| {
        let opts = BdfOptions { fixed_step: Some(t_end / steps as f64), max_order, ..BdfOptions::new(1e-6, 1e-9) };
        let sol = bdf_solve(&m, &k, &mut |_, f| f.fill(0.0), 0.0, t_end, &y0, &[t_end], &opts).unwrap();
        max_err(&sol.y[0], &heat_exact(n, modes, t_end))
    };
    // BDF1 and BDF2 (the start-up steps of order 1 limit the higher orders to 2 for this fixed-step
    // ramp; the adaptive test above exercises the high orders)
    for (order, lo, hi) in [(1usize, 1.7, 2.3), (2, 3.3, 4.7)] {
        let (e1, e2) = (err_at(40, order), err_at(80, order));
        let ratio = e1 / e2;
        eprintln!("BDF{order}: {e1:.3e} -> {e2:.3e}, ratio {ratio:.2}");
        assert!(ratio > lo && ratio < hi, "order {order}: error ratio {ratio}");
    }
}

#[test]
fn theta_method_orders() {
    let n = 100;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let modes = (5, 12);
    let y0 = heat_exact(n, modes, 0.0);
    let err_at = |steps: usize, theta: f64| {
        let sol =
            theta_solve(&m, &k, &mut |_, f| f.fill(0.0), 0.0, 1.0, steps, &y0, theta, steps, Definiteness::General)
                .unwrap();
        max_err(sol.y.last().unwrap(), &heat_exact(n, modes, 1.0))
    };
    let be = err_at(40, 1.0) / err_at(80, 1.0);
    let cn = err_at(40, 0.5) / err_at(80, 0.5);
    assert!(be > 1.7 && be < 2.3, "backward Euler ratio {be}");
    assert!(cn > 3.4 && cn < 4.6, "Crank–Nicolson ratio {cn}");
}

fn wave_exact(n: usize, ka: usize, t: f64) -> (Vec<f64>, Vec<f64>) {
    let w = laplacian_1d_eigenvalues(n)[ka - 1].sqrt();
    let s = sine(n, ka);
    (s.iter().map(|x| x * (w * t).cos()).collect(), s.iter().map(|x| -w * x * (w * t).sin()).collect())
}

#[test]
fn newmark_conserves_energy_and_is_second_order() {
    let n = 150;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let ka = 10usize;
    let (u0, v0) = wave_exact(n, ka, 0.0);
    let energy = |u: &[f64], v: &[f64]| {
        let mut ku = vec![0.0; n];
        k.spmv(u, &mut ku);
        0.5 * v.iter().map(|x| x * x).sum::<f64>() + 0.5 * u.iter().zip(&ku).map(|(a, b)| a * b).sum::<f64>()
    };
    let e0 = energy(&u0, &v0);
    let run = |steps: usize, p: NewmarkParams| {
        newmark_solve(&m, None, &k, &mut |_, f| f.fill(0.0), 0.0, 20.0, steps, &u0, &v0, p, 1, Definiteness::General)
            .unwrap()
    };
    let sol = run(400, NewmarkParams::average_acceleration());
    for (u, v) in sol.u.iter().zip(&sol.v) {
        assert!((energy(u, v) - e0).abs() < 1e-11 * e0, "energy drift {}", (energy(u, v) - e0) / e0);
    }
    assert_eq!(sol.stats.factorizations, 2); // S and M, once each
    // second order for both families
    for (name, p) in [
        ("average acceleration", NewmarkParams::average_acceleration()),
        ("generalised-α ρ∞ = 0.8", NewmarkParams::generalized_alpha(0.8)),
    ] {
        let err = |steps: usize| {
            let s = run(steps, p);
            let (uw, _) = wave_exact(n, ka, 20.0);
            max_err(s.u.last().unwrap(), &uw)
        };
        let ratio = err(200) / err(400);
        eprintln!("{name}: ratio {ratio:.2}");
        assert!(ratio > 3.4 && ratio < 4.6, "{name}: error ratio {ratio}");
    }
}

#[test]
fn generalised_alpha_damps_the_highest_modes_newmark_does_not() {
    let n = 100;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let u0 = sine(n, n); // highest mode, ω ≈ 2
    let v0 = vec![0.0; n];
    let norm_end = |p: NewmarkParams| {
        let s = newmark_solve(
            &m,
            None,
            &k,
            &mut |_, f| f.fill(0.0),
            0.0,
            100.0,
            20,
            &u0,
            &v0,
            p,
            20,
            Definiteness::General,
        )
        .unwrap(); // h = 5, ωh ≈ 10
        let (u, v) = (s.u.last().unwrap(), s.v.last().unwrap());
        (u.iter().map(|x| x * x).sum::<f64>() + v.iter().map(|x| x * x).sum::<f64>()).sqrt()
    };
    let start = norm(&u0);
    assert!(norm_end(NewmarkParams::average_acceleration()) > 0.2 * start);
    assert!(
        norm_end(NewmarkParams::generalized_alpha(0.3)) < 1e-4 * start,
        "{}",
        norm_end(NewmarkParams::generalized_alpha(0.3))
    );
}

#[test]
fn damped_second_order_system_decays_like_the_modal_solution() {
    // M = I, C = c I, K = L: mode k: u'' + c u' + λ u = 0 with a known underdamped solution.
    let n = 80;
    let (m, k) = (Csr::<f64>::identity(n), laplacian_1d(n));
    let cmat = Csr::axpby(0.2, &Csr::identity(n), 0.0, &Csr::identity(n));
    let ka = 6usize;
    let lam = laplacian_1d_eigenvalues(n)[ka - 1];
    let s = sine(n, ka);
    let (alpha, wd) = (0.1, (lam - 0.01f64).sqrt());
    let t_end = 15.0;
    let sol = newmark_solve(
        &m,
        Some(&cmat),
        &k,
        &mut |_, f| f.fill(0.0),
        0.0,
        t_end,
        600,
        &s,
        &vec![0.0; n],
        NewmarkParams::average_acceleration(),
        600,
        Definiteness::General,
    )
    .unwrap();
    let c = (-alpha * t_end).exp() * ((wd * t_end).cos() + alpha / wd * (wd * t_end).sin());
    let want: Vec<f64> = s.iter().map(|x| c * x).collect();
    assert!(max_err(sol.u.last().unwrap(), &want) < 2e-4);
}

#[test]
fn rosenbrock_logistic_fixed_step_order_two_and_adaptive_tolerance() {
    let n = 5;
    let y0: Vec<f64> = (0..n).map(|i| 0.05 + 0.1 * i as f64).collect();
    let m = Csr::<f64>::identity(n);
    let exact = |t: f64| -> Vec<f64> { y0.iter().map(|&y| 1.0 / (1.0 + (1.0 / y - 1.0) * (-t).exp())).collect() };
    let run = |opts: &RosenbrockOptions| {
        let mut f = |y: &[f64], out: &mut [f64]| {
            for (o, &v) in out.iter_mut().zip(y) {
                *o = v * (1.0 - v);
            }
        };
        let mut jac = |y: &[f64]| {
            let t: Vec<_> = y.iter().enumerate().map(|(i, &v)| (i, i, 1.0 - 2.0 * v)).collect();
            Csr::from_triplets(n, n, &t).unwrap()
        };
        rosenbrock_solve(&m, &mut f, &mut jac, 0.0, 6.0, &y0, &[6.0], opts).unwrap()
    };
    let fixed = |h: f64| {
        let o = RosenbrockOptions { fixed_step: true, ..RosenbrockOptions::new(1e-6, 1e-9, h) };
        max_err(&run(&o).y[0], &exact(6.0))
    };
    let ratio = fixed(0.1) / fixed(0.05);
    assert!(ratio > 3.4 && ratio < 4.6, "ROS2 fixed-step ratio {ratio}");
    let mut errs = Vec::new();
    for rtol in [1e-3, 1e-5, 1e-7] {
        let sol = run(&RosenbrockOptions::new(rtol, 1e-10, 0.01));
        errs.push(max_err(&sol.y[0], &exact(6.0)));
    }
    eprintln!("ROS2 adaptive errors {errs:?}");
    assert!(errs[0] > errs[1] && errs[1] > errs[2], "{errs:?}");
    // Per-step error control gives a global error ∝ tol^(p/(p+1)) = tol^(2/3) for this second-order
    // method: a factor ≈ 22 per two decades, not 100.
    assert!(errs[0] / errs[1] > 10.0 && errs[1] / errs[2] > 10.0, "{errs:?}");
    assert!(errs[2] < 1e-4 && errs[0] < 1e-1, "{errs:?}");
}

#[test]
fn rosenbrock_stiff_nonlinear_heat_agrees_with_a_tight_reference() {
    // M y' = −κ L y − y³ + 1 (autonomous), n = 100, stiff (κ L has eigenvalues up to ≈ 4κ).
    let n = 100;
    let kappa = 2000.0;
    let l = laplacian_1d(n);
    let m = Csr::<f64>::identity(n);
    let y0 = vec![0.0; n];
    let run = |rtol: f64| {
        let mut f = |y: &[f64], out: &mut [f64]| {
            l.spmv(y, out);
            for (o, &v) in out.iter_mut().zip(y) {
                *o = -kappa * *o - v * v * v + 1.0;
            }
        };
        let mut jac = |y: &[f64]| {
            let d: Vec<_> = y.iter().enumerate().map(|(i, &v)| (i, i, -3.0 * v * v)).collect();
            Csr::axpby(-kappa, &l, 1.0, &Csr::from_triplets(n, n, &d).unwrap())
        };
        rosenbrock_solve(&m, &mut f, &mut jac, 0.0, 1.0, &y0, &[1.0], &RosenbrockOptions::new(rtol, 1e-10, 1e-4))
            .unwrap()
    };
    let reference = run(1e-9);
    let sol = run(1e-5);
    let e = max_err(&sol.y[0], &reference.y[0]);
    eprintln!("stiff nonlinear heat: error {e:.2e}, steps {} (reference {})", sol.stats.steps, reference.stats.steps);
    assert!(e < 1e-3, "error {e}");
    // L-stability: far fewer steps than the explicit stability limit h ≈ 2/(4κ) = 2.5e-4 → 4000 steps
    assert!(sol.stats.steps < 400, "{} steps", sol.stats.steps);
}
