// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What identifies a mode of the axisymmetric study (ADR 80): the numeric scalars `f re N (Hz)` and `f im N (Hz)`, and the energy shares
//! (`int rho |u|^2 r` outside the layers) per displacement component and per domain. Checked on closed forms whose family is known
//! (torsion: azimuthal only; radial plane strain: radial only), on a two-domain body and on a body with an absorbing layer.
mod common;
use femcore::props::*;
use nyx_physics::{Discretisation, Scalar, SolveResult};
use nyx_sms::*;

const PI: f64 = std::f64::consts::PI;

/// Bessel function of the first kind by the periodic trapezoid rule on `J_n(x) = (1 / 2 pi) int cos(n t - x sin t) dt`.
fn bessel(n: i32, x: f64) -> f64 {
    let steps = 400;
    (0..steps).map(|k| 2.0 * PI * k as f64 / steps as f64).map(|t| (n as f64 * t - x * t.sin()).cos()).sum::<f64>()
        / steps as f64
}

fn problem(n_modes: usize, f_search: f64, domains: &[i32]) -> Problem {
    let mut sp = Problem::new(Study::Static, 2);
    sp.order = 2;
    for &d in domains {
        sp.materials.set(d, "E", 169e9);
        sp.materials.set(d, "nu", 0.28);
        sp.materials.set(d, "rho", 2329.0);
    }
    sp.study = Study::Eigenfrequency { n_modes, f_search };
    sp
}

fn num(r: &SolveResult, name: &str) -> f64 {
    match r.scalars.iter().find(|(n, _)| n == name) {
        Some((_, Scalar::Num(x))) => *x,
        other => panic!("scalar '{name}' is not a number: {other:?}"),
    }
}

fn has(r: &SolveResult, name: &str) -> bool {
    r.scalars.iter().any(|(n, _)| n == name)
}

/// The mode nearest `target`, one-based.
fn nearest(r: &SolveResult, target: f64) -> usize {
    (1..=r.solutions.len())
        .min_by(|a, b| {
            (num(r, &format!("f re {a} (Hz)")) - target)
                .abs()
                .total_cmp(&(num(r, &format!("f re {b} (Hz)")) - target).abs())
        })
        .unwrap()
}

fn shares(r: &SolveResult, j: usize) -> [f64; 3] {
    [num(r, &format!("radial share {j}")), num(r, &format!("axial share {j}")), num(r, &format!("azimuthal share {j}"))]
}

fn end_labels(mesh: &Discretisation, len: f64) -> Vec<i32> {
    let mut labels = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        if fv.iter().all(|&v| (mesh.point(v as usize)[1].abs() - len / 2.0).abs() < 1e-9 * len) {
            labels.insert(mesh.facet_labels[f]);
        }
    }
    labels.into_iter().collect()
}

/// T1 and T2 (lossless): every mode has numeric `f re` and `f im` equal to its eigenvalue, shares that sum to one, and the torsional mode of
/// a solid cylinder (`f = 5.13562 c_s / (2 pi R)`) has azimuthal share 1.
#[test]
fn a_lossless_mode_has_numeric_frequency_parts_and_a_torsional_mode_is_azimuthal() {
    let cs = (169e9_f64 / (2.0 * 1.28) / 2329.0).sqrt();
    let want = 5.13562 * cs / (2.0 * PI);
    let mesh = common::parts_mesh(&[common::part(0.0, -0.5, 1.0, 0.5, 1, Some(0.1))], 0.1);
    let (r, _) = solve_solid_axisym(&mesh, &problem(6, want, &[1]), 0.0).unwrap();
    for (j, s) in r.solutions.iter().enumerate() {
        let [re, im] = s.eigenvalue.unwrap();
        assert_eq!(num(&r, &format!("f re {} (Hz)", j + 1)), re);
        assert_eq!(num(&r, &format!("f im {} (Hz)", j + 1)), im);
        let sh = shares(&r, j + 1);
        assert!((sh.iter().sum::<f64>() - 1.0).abs() < 1e-12, "{sh:?}");
        assert!((num(&r, &format!("share in domain 1 {}", j + 1)) - 1.0).abs() < 1e-12);
    }
    let k = nearest(&r, want);
    let sh = shares(&r, k);
    println!("torsion: f = {} Hz against {want}, shares {sh:?}", num(&r, &format!("f re {k} (Hz)")));
    assert!((num(&r, &format!("f re {k} (Hz)")) - want).abs() < 5e-4 * want);
    assert!((sh[2] - 1.0).abs() < 1e-12 && sh[0] < 1e-12 && sh[1] < 1e-12, "{sh:?}");
}

/// T2: the radial plane-strain mode of a cylinder held in `z` at its ends (`k R = 2.0288` for `nu = 0.28`) is radial.
#[test]
fn a_radial_plane_strain_mode_is_radial() {
    let (e, nu, rho) = (169e9_f64, 0.28_f64, 2329.0_f64);
    let cp = (e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) / rho).sqrt();
    // first root of x J0(x) / J1(x) = (1 - 2 nu) / (1 - nu) (the closed form of `axisym_axis.rs`, T2 there), by bisection on (1.5, 2.4)
    let g = |x: f64| x * bessel(0, x) / bessel(1, x) - (1.0 - 2.0 * nu) / (1.0 - nu);
    let (mut lo, mut hi) = (1.5_f64, 2.4_f64);
    assert!(g(lo) * g(hi) < 0.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if g(lo) * g(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let want = 0.5 * (lo + hi) * cp / (2.0 * PI);
    let mesh = common::parts_mesh(&[common::part(0.0, -0.25, 1.0, 0.25, 1, Some(0.1))], 0.1);
    let mut sp = problem(6, want, &[1]);
    sp.bcs.push(Bc {
        labels: end_labels(&mesh, 0.5),
        kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] },
    });
    let (r, _) = solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
    let k = nearest(&r, want);
    let sh = shares(&r, k);
    println!("radial: f = {} Hz against {want}, shares {sh:?}", num(&r, &format!("f re {k} (Hz)")));
    assert!((num(&r, &format!("f re {k} (Hz)")) - want).abs() < 1e-3 * want);
    assert!(sh[2] < 1e-12, "a radial mode at m = 0 has no azimuthal part: {sh:?}");
    assert!(sh[0] > 1.0 - 1e-9, "{sh:?}");
}

/// T2: the shares of a body of two domains sum to one over the domains, and each domain holds some of the energy.
#[test]
fn the_domain_shares_of_a_two_domain_body_sum_to_one() {
    let parts = [common::part(0.0, 0.0, 1.0, 0.5, 1, Some(0.1)), common::part(1.0, 0.0, 2.0, 0.5, 2, Some(0.1))];
    let mesh = common::parts_mesh(&parts, 0.1);
    let (r, _) = solve_solid_axisym(&mesh, &problem(6, 0.0, &[1, 2]), 0.0).unwrap();
    for j in 1..=r.solutions.len() {
        let (a, b) = (num(&r, &format!("share in domain 1 {j}")), num(&r, &format!("share in domain 2 {j}")));
        assert!((a + b - 1.0).abs() < 1e-12, "{a} + {b}");
        assert!(a > 0.0 && b > 0.0);
    }
}

/// T1 (lossy): with an absorbing layer `f re` and `f im` are numbers next to the text `fN (Hz)`, `f im` is the imaginary part of the
/// eigenvalue, the shares are taken outside the layer (no share for the layer's domain) and still sum to one.
#[test]
fn a_lossy_mode_has_numeric_frequency_parts_and_shares_outside_the_layer() {
    let parts = [common::part(4.0, 0.0, 6.0, 0.4, 1, Some(0.1)), common::part(6.0, 0.0, 7.0, 0.4, 2, Some(0.1))];
    let mesh = common::parts_mesh(&parts, 0.1);
    let mut sp = problem(4, 0.0, &[1, 2]);
    sp.materials.pml.insert(
        2,
        Pml {
            x: PmlDir::Plus,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 3.0,
            order: 2.0,
            re: 1.0,
            ..Default::default()
        },
    );
    let (r, _) = solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
    let mut lossy = 0;
    for (j, s) in r.solutions.iter().enumerate() {
        let [re, im] = s.eigenvalue.unwrap();
        assert_eq!(num(&r, &format!("f re {} (Hz)", j + 1)), re);
        assert_eq!(num(&r, &format!("f im {} (Hz)", j + 1)), im);
        if im != 0.0 {
            lossy += 1;
            assert!(matches!(
                r.scalars.iter().find(|(n, _)| *n == format!("f{} (Hz)", j + 1)),
                Some((_, Scalar::Text(_)))
            ));
        }
        assert!((shares(&r, j + 1).iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((num(&r, &format!("share in domain 1 {}", j + 1)) - 1.0).abs() < 1e-12);
        assert!(!has(&r, &format!("share in domain 2 {}", j + 1)), "the layer has no share");
    }
    assert!(lossy > 0, "no mode got an imaginary part from the layer");
}
