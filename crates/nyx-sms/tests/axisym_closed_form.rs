// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The axisymmetric solid kernel (`sms::axisym`, through `solve_solid_axisym`) against closed forms (gate G1: paul asked that its
//! implementation be checked; the electromagnetic one is checked in `femcore/tests/axisym_closed_form.rs`). Printed relative errors are
//! in the test output (`-- --nocapture`).
mod common;
use nyx_sms::*;

fn ring(m: f64, n_modes: usize, f_search: f64) -> Vec<f64> {
    // radius 1 m, section 0.02 x 0.02 m (mesh x = r, y = z), isotropic silicon-like
    let (rr, w) = (1.0, 0.02);
    let mesh =
        common::parts_mesh(&[common::part(rr - w / 2.0, -w / 2.0, rr + w / 2.0, w / 2.0, 1, Some(w / 4.0))], w / 4.0);
    let mut sp = Problem::new(Study::Static, 2);
    sp.order = 2;
    sp.materials.set(1, "E", 169e9);
    sp.materials.set(1, "nu", 0.28);
    sp.materials.set(1, "rho", 2329.0);
    sp.study = Study::Eigenfrequency { n_modes, f_search };
    let (_, modes) = nyx_sms::solve_solid_axisym(&mesh, &sp, m).unwrap();
    let mut f: Vec<f64> = modes.freqs.iter().map(|z| z.re).collect();
    f.sort_by(|a, b| a.total_cmp(b));
    f
}

fn nearest(f: &[f64], target: f64) -> f64 {
    f.iter().copied().min_by(|a, b| (a - target).abs().total_cmp(&(b - target).abs())).unwrap()
}

/// Thin ring: breathing mode `f = sqrt(E / rho) / (2 pi R)` (corrections of order (w/R)^2), in-plane flexural modes
/// `f_m = m (m^2 - 1) / sqrt(m^2 + 1) * (w / sqrt 12) / (2 pi R^2) * sqrt(E / rho)` (Euler-Bernoulli, error of order (m w / R)^2).
#[test]
fn a_thin_ring_has_the_breathing_and_flexural_frequencies_of_ring_theory() {
    let (e, rho, rr, w) = (169e9_f64, 2329.0_f64, 1.0_f64, 0.02_f64);
    let c = (e / rho).sqrt();
    let breathing = c / (2.0 * std::f64::consts::PI * rr);
    let f0 = ring(0.0, 4, breathing);
    let got = nearest(&f0, breathing);
    println!(
        "ring breathing (m = 0): {got:.4} Hz against {breathing:.4} Hz, relative error {:.2e}",
        (got - breathing).abs() / breathing
    );
    assert!((got - breathing).abs() < 2e-3 * breathing, "{got} against {breathing}");
    for m in [2.0_f64, 3.0] {
        let flex =
            m * (m * m - 1.0) / (m * m + 1.0).sqrt() * (w / 12f64.sqrt()) / (2.0 * std::f64::consts::PI * rr * rr) * c;
        let f = ring(m, 4, flex);
        let got = nearest(&f, flex);
        println!(
            "ring flexure m = {m}: {got:.4} Hz against {flex:.4} Hz, relative error {:.2e}",
            (got - flex).abs() / flex
        );
        assert!((got - flex).abs() < 2e-2 * flex, "m = {m}: {got} against {flex}; modes {f:?}");
    }
}
