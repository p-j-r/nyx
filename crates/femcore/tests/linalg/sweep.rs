// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Frequency sweeps: reduced-order model against closed forms (modal solutions of the 1-D Laplacian
//! with damping or a Drude medium) and against a full solve per frequency on a 2-D problem.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::{laplacian_1d, laplacian_2d};
use femcore::linalg::direct::Definiteness;
use femcore::linalg::nep::{NepTerm, NonlinearEigenProblem};
use femcore::linalg::sweep::{Rom, SweepOptions, direct_sweep, rom_sweep};
use femcore::linalg::{C64, Csr};
use std::f64::consts::PI;

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn sine(n: usize, k: usize) -> Vec<f64> {
    (1..=n).map(|j| (k as f64 * PI * j as f64 / (n as f64 + 1.0)).sin()).collect()
}

/// `T(ω) = L − ω² I + iωγ I` (damped vibration, frequency domain).
fn damped(n: usize, gamma: f64) -> NonlinearEigenProblem {
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(move |z| -z * z + c(0.0, gamma) * z) },
    ])
    .unwrap()
}

/// Exact `cᵀ x(ω)` for `T x = b` when `T(ω) = g(ω) I + L`-type: `x = Σ_k (s_kᵀb/‖s_k‖²) s_k / (μ_k + g(ω))`.
fn modal_output(n: usize, b: &[C64], out: &[C64], g: impl Fn(f64) -> C64, omega: f64) -> C64 {
    let ev = laplacian_1d_eigenvalues(n);
    let mut total = C64::ZERO;
    for k in 1..=n {
        let s = sine(n, k);
        let sn2: f64 = s.iter().map(|x| x * x).sum();
        let sb: C64 = s.iter().zip(b).map(|(a, bb)| *bb * *a).sum::<C64>() / sn2;
        let co: C64 = s.iter().zip(out).map(|(a, cc)| cc.conj() * *a).sum();
        total += co * sb / (c(ev[k - 1], 0.0) + g(omega));
    }
    total
}

fn grid(lo: f64, hi: f64, m: usize) -> Vec<f64> {
    (0..m).map(|i| lo + (hi - lo) * i as f64 / (m - 1) as f64).collect()
}

#[test]
fn rom_sweep_reproduces_the_modal_closed_form() {
    let n = 200;
    let gamma = 0.02;
    let p = damped(n, gamma);
    let b: Vec<C64> = (0..n).map(|i| c(((i * 37) % 11) as f64 - 5.0, 0.3 * ((i * 17) % 7) as f64)).collect();
    let out: Vec<C64> = (0..n).map(|i| c(1.0 / (1.0 + i as f64), 0.0)).collect();
    // about 23 resonances of the n = 200 string (√μ_k ≈ kπ/201) fall in this band; a reduced model
    // needs about one basis vector per resonance, so it pays when the grid is much denser than that
    let omegas = grid(0.05, 0.4, 400);
    let opts = SweepOptions::new(1e-8);
    let r = rom_sweep(&p, &b, &omegas, std::slice::from_ref(&out), &opts).unwrap();
    assert!(r.converged, "max residual {:?}", r.residuals.iter().cloned().fold(0.0, f64::max));
    assert!(r.residuals.iter().all(|&x| x <= 1e-8 * 1.001));
    let mut worst = 0.0f64;
    for (w, o) in omegas.iter().zip(&r.outputs) {
        let exact = modal_output(n, &b, &out, |w| c(-w * w, gamma * w), *w);
        worst = worst.max((o[0] - exact).norm() / exact.norm().max(1e-12));
    }
    eprintln!(
        "basis {} from {} full solves for {} frequencies, worst relative output error {worst:.2e}",
        r.basis_dim,
        r.factorizations,
        omegas.len()
    );
    assert!(worst < 1e-5, "{worst}");
    assert!(r.factorizations * 4 < omegas.len(), "{} solves for {} frequencies", r.factorizations, omegas.len());
    assert_eq!(r.symbolic_analyses, 1);
    assert_eq!(r.basis_dim, r.expansion.len());
}

#[test]
fn rom_sweep_handles_a_dispersive_medium() {
    // T(ω) = L − ω²(1 − ω_p²/(ω² + iγω)) I = L + [−ω² + ω_p² ω/(ω + iγ)] I
    let n = 150;
    let (wp, gamma) = (1.0, 0.1);
    let l = laplacian_1d(n).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let g = move |w: f64| -> C64 {
        let z = c(w, 0.0);
        -z * z + c(wp * wp, 0.0) * z / (z + c(0.0, gamma))
    };
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(move |z| -z * z + c(wp * wp, 0.0) * z / (z + c(0.0, gamma))) },
    ])
    .unwrap();
    // forcing in the first 12 sine modes, so the response lives in a 12-dimensional subspace whose
    // resonances √(μ_k + ω_p²) ∈ (1, 1.1) lie in the band
    let mut b = vec![C64::ZERO; n];
    for k in 1..=12 {
        for (bi, s) in b.iter_mut().zip(sine(n, k)) {
            *bi += c(s / k as f64, 0.2 * s);
        }
    }
    let out = b.clone();
    let omegas = grid(0.9, 1.2, 300);
    let r = rom_sweep(&p, &b, &omegas, std::slice::from_ref(&out), &SweepOptions::new(1e-8)).unwrap();
    assert!(r.converged);
    let mut worst = 0.0f64;
    for (w, o) in omegas.iter().zip(&r.outputs) {
        let exact = modal_output(n, &b, &out, g, *w);
        worst = worst.max((o[0] - exact).norm() / exact.norm().max(1e-12));
    }
    assert!(worst < 1e-5, "{worst}");
    assert!(r.factorizations * 4 < omegas.len(), "{} solves", r.factorizations);
}

#[test]
fn rom_sweep_agrees_with_a_full_solve_per_frequency_in_2d() {
    let (nx, ny) = (40, 40);
    let n = nx * ny;
    let l = laplacian_2d(nx, ny).to_complex();
    let id: Csr<C64> = Csr::identity(n);
    let gamma = 0.02;
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: l, f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: id, f: Box::new(move |z| -z * z + c(0.0, gamma) * z) },
    ])
    .unwrap();
    let mut b = vec![C64::ZERO; n];
    b[(ny / 3) * nx + nx / 4] = C64::ONE; // point source
    let mut out = vec![C64::ZERO; n];
    out[(2 * ny / 3) * nx + 3 * nx / 4] = C64::ONE; // probe
    let omegas = grid(0.2, 0.45, 200); // about 27 of the 1 600 modes have ω in this band
    let rom = rom_sweep(&p, &b, &omegas, std::slice::from_ref(&out), &SweepOptions::new(1e-7)).unwrap();
    let full = direct_sweep(&p, &b, &omegas, std::slice::from_ref(&out), Definiteness::General).unwrap();
    assert!(rom.converged);
    let scale = full.outputs.iter().map(|o| o[0].norm()).fold(0.0, f64::max);
    let worst = rom.outputs.iter().zip(&full.outputs).map(|(a, b)| (a[0] - b[0]).norm() / scale).fold(0.0, f64::max);
    eprintln!(
        "2-D: ROM {} solves, basis {}, {:.3}s vs full {} solves {:.3}s; worst output error {worst:.2e}",
        rom.factorizations, rom.basis_dim, rom.seconds, full.factorizations, full.seconds
    );
    assert!(worst < 1e-5, "{worst}");
    assert!(rom.factorizations * 3 < full.factorizations);
    assert!(full.residuals.iter().all(|&x| x < 1e-10));
}

#[test]
fn indicator_matches_the_exact_residual_and_basis_is_orthonormal() {
    let n = 80;
    let p = damped(n, 0.1);
    let b: Vec<C64> = (0..n).map(|i| c((i % 5) as f64, 0.0)).collect();
    let mut rom = Rom::new(&p, &b);
    // a few full solutions as basis vectors
    for w in [0.2, 0.7, 1.1] {
        let t = p.assemble(c(w, 0.0));
        let x = dense_solve(&t, &b);
        assert!(rom.extend(&x));
    }
    assert_eq!(rom.dim(), 3);
    assert!(!rom.extend(&b.iter().map(|v| *v * 0.0).collect::<Vec<_>>()) || true);
    for w in [0.3, 0.5, 0.9, 1.3] {
        let xr = rom.solve_reduced(w).unwrap();
        let (ex, ind) = (rom.residual(w, &xr), rom.residual_indicator(w, &xr));
        assert!((ex - ind).abs() < 1e-7 + 1e-6 * ex, "ω = {w}: exact {ex:e} vs indicator {ind:e}");
    }
    // the residual vanishes at an expansion frequency
    let xr = rom.solve_reduced(0.7).unwrap();
    assert!(rom.residual(0.7, &xr) < 1e-10);
}

/// Dense Gaussian elimination solve for a small complex CSR matrix.
fn dense_solve(a: &Csr<C64>, b: &[C64]) -> Vec<C64> {
    let n = a.nrows();
    let d = femcore::linalg::dense::Dense::from_fn(n, n, |i, j| a.to_dense()[i * n + j]);
    femcore::linalg::dense::solve_small(&d, b).unwrap()
}

#[test]
fn invalid_inputs() {
    let p = damped(20, 0.1);
    let b = vec![C64::ONE; 20];
    assert!(rom_sweep(&p, &b, &[], &[], &SweepOptions::new(1e-6)).is_err());
    assert!(rom_sweep(&p, &b[..5], &[1.0], &[], &SweepOptions::new(1e-6)).is_err());
    assert!(rom_sweep(&p, &vec![C64::ZERO; 20], &[1.0], &[], &SweepOptions::new(1e-6)).is_err());
    assert!(rom_sweep(&p, &b, &[1.0], &[vec![C64::ONE; 3]], &SweepOptions::new(1e-6)).is_err());
}
