//! LOBPCG: closed-form spectra (real, generalised, complex Hermitian), constraints, highest
//! eigenvalues, soft locking and the stall test.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::amg::{Amg, AmgOptions};
use femcore::solvers::csr::{laplacian_1d, laplacian_2d, mass_1d};
use femcore::solvers::linop::Preconditioner;
use femcore::solvers::lobpcg::{LobpcgOptions, LobpcgReason, LobpcgWhich, lobpcg};
use femcore::solvers::{C64, Csr, Scalar};
use std::f64::consts::PI;

fn fe_spectrum(n: usize) -> Vec<f64> {
    (1..=n)
        .map(|k| {
            let c = (k as f64 * PI / (n as f64 + 1.0)).cos();
            6.0 * (1.0 - c) / (2.0 + c)
        })
        .collect()
}

fn sine(n: usize, k: usize) -> Vec<f64> {
    (1..=n).map(|j| (k as f64 * PI * j as f64 / (n as f64 + 1.0)).sin()).collect()
}

struct AmgPc(Amg<f64>);
impl Preconditioner<f64> for AmgPc {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.0.apply(r, z)
    }
}

#[test]
fn lowest_modes_of_the_2d_laplacian_with_an_amg_preconditioner() {
    let (nx, ny) = (40, 40);
    let a = laplacian_2d(nx, ny);
    let exact = laplacian_2d_eigenvalues(nx, ny);
    let pc = AmgPc(Amg::new(&a, AmgOptions::default()).unwrap());
    let r = lobpcg(&a, None, Some(&pc), None, &[], &LobpcgOptions { tol: 1e-9, ..LobpcgOptions::new(8) }).unwrap();
    assert_eq!(r.reason, LobpcgReason::Converged, "{:?} after {} its", r.reason, r.iterations);
    for (g, w) in r.values.iter().zip(&exact) {
        assert!((g - w).abs() < 1e-9, "{g} vs {w}");
    }
    assert!(r.residuals.iter().all(|&e| e < 1e-8));
    // B-orthonormal (here Euclidean) vectors
    for i in 0..8 {
        for j in 0..=i {
            let d: f64 = r.vectors[i].iter().zip(&r.vectors[j]).map(|(p, q)| p * q).sum();
            assert!((d - if i == j { 1.0 } else { 0.0 }).abs() < 1e-8);
        }
    }
    // the preconditioner matters
    let plain =
        lobpcg(&a, None, None, None, &[], &LobpcgOptions { tol: 1e-9, max_it: 2000, ..LobpcgOptions::new(8) }).unwrap();
    eprintln!("2-D Laplacian: {} iterations with AMG, {} without ({:?})", r.iterations, plain.iterations, plain.reason);
    assert!(r.iterations * 3 < plain.iterations || plain.reason != LobpcgReason::Converged);
}

#[test]
fn generalised_finite_element_problem_and_complex_hermitian() {
    let n = 400;
    let a = laplacian_1d(n);
    let b = mass_1d(n, 1.0);
    let r = lobpcg(&a, Some(&b), None, None, &[], &LobpcgOptions { tol: 1e-8, max_it: 3000, ..LobpcgOptions::new(5) })
        .unwrap();
    assert_eq!(r.reason, LobpcgReason::Converged, "{:?} after {}", r.reason, r.iterations);
    for (g, w) in r.values.iter().zip(fe_spectrum(n)) {
        assert!((g - w).abs() < 1e-8 * w.max(1.0), "{g} vs {w}");
    }
    // complex Hermitian chain with phase-twisted hoppings: unitarily equivalent to the Laplacian
    let m = 300;
    let mut t = Vec::new();
    for i in 0..m {
        t.push((i, i, C64::new(2.0, 0.0)));
        if i + 1 < m {
            let ph = C64::new(0.0, 0.3 * (1 + i % 3) as f64).exp();
            t.push((i, i + 1, ph.scale(-1.0)));
            t.push((i + 1, i, ph.conj().scale(-1.0)));
        }
    }
    let h = Csr::from_triplets(m, m, &t).unwrap();
    assert!(h.is_hermitian(1e-14));
    let rc =
        lobpcg(&h, None, None, None, &[], &LobpcgOptions { tol: 1e-9, max_it: 4000, ..LobpcgOptions::new(4) }).unwrap();
    assert_eq!(rc.reason, LobpcgReason::Converged, "{:?} after {}", rc.reason, rc.iterations);
    for (g, w) in rc.values.iter().zip(laplacian_1d_eigenvalues(m)) {
        assert!((g - w).abs() < 1e-8, "{g} vs {w}");
    }
}

#[test]
fn constraints_deflate_known_modes() {
    let n = 300;
    let a = laplacian_1d(n);
    let known: Vec<Vec<f64>> = (1..=3).map(|k| sine(n, k)).collect();
    let r = lobpcg(&a, None, None, None, &known, &LobpcgOptions { tol: 1e-9, max_it: 3000, ..LobpcgOptions::new(4) })
        .unwrap();
    assert_eq!(r.reason, LobpcgReason::Converged, "{:?} after {}", r.reason, r.iterations);
    let ev = laplacian_1d_eigenvalues(n);
    for (g, w) in r.values.iter().zip(&ev[3..7]) {
        assert!((g - w).abs() < 1e-8, "{g} vs {w}");
    }
    // orthogonal to the constraints
    for v in &r.vectors {
        for y in &known {
            let d: f64 = v.iter().zip(y).map(|(p, q)| p * q).sum::<f64>() / norm(y);
            assert!(d.abs() < 1e-7, "{d}");
        }
    }
}

#[test]
fn highest_eigenvalues_soft_locking_and_the_stall_test() {
    let n = 200;
    let a = laplacian_1d(n);
    let ev = laplacian_1d_eigenvalues(n);
    let hi = lobpcg(
        &a,
        None,
        None,
        None,
        &[],
        &LobpcgOptions { which: LobpcgWhich::Highest, tol: 1e-8, max_it: 3000, ..LobpcgOptions::new(3) },
    )
    .unwrap();
    assert_eq!(hi.reason, LobpcgReason::Converged, "{:?}", hi.reason);
    for (g, w) in hi.values.iter().zip(ev.iter().rev()) {
        assert!((g - w).abs() < 1e-7, "{g} vs {w}");
    }
    // soft locking: the number of active columns shrinks before the iteration ends
    let r = lobpcg(
        &a,
        None,
        None,
        None,
        &[],
        &LobpcgOptions { tol: 1e-8, max_it: 4000, block: 10, ..LobpcgOptions::new(8) },
    )
    .unwrap();
    assert_eq!(r.reason, LobpcgReason::Converged);
    let min_active = *r.active_history.iter().min().unwrap();
    eprintln!("active columns: start {}, minimum {min_active}, {} iterations", r.active_history[0], r.iterations);
    assert!(min_active < r.active_history[0]);
    // stall: an unreachable tolerance ends with Stalled (or max_it), not an endless loop
    let st = lobpcg(
        &a,
        None,
        None,
        None,
        &[],
        &LobpcgOptions { tol: 1e-17, max_it: 3000, stall_window: 15, ..LobpcgOptions::new(3) },
    )
    .unwrap();
    assert_ne!(st.reason, LobpcgReason::Converged);
    assert!(st.iterations < 3000 || st.reason == LobpcgReason::MaxIterations);
    // invalid input
    assert!(lobpcg(&a, None, None, None, &[], &LobpcgOptions::new(0)).is_err());
    let small = laplacian_1d(6);
    assert!(lobpcg(&small, None, None, None, &[], &LobpcgOptions::new(3)).is_err());
    let _ = Scalar::modulus(1.0f64);
}
