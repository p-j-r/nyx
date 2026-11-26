// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME VVUQ Pillar (a) analytical closed-form verification suite for `nyx-sms`:
//! - 3D isotropic cube traction-free Lamé, guided-face, and clamped eigenfrequencies (`< 5e-3` on P2)
//! - Kirsch plate with a circular hole stress-concentration factor `K_t = 3.0` (`< 3%` finite-width corrected)
//! - Analytical Rayleigh-Lamb dispersion (`S_n`, `A_n`) for `q != 0` in a free plate (`< 1e-3`)
//! - Solid-cylinder `m = 2` axisymmetric coupled `(u_r, u_theta)` modes (`< 2e-3`)
//! - Free-disk radial contour mode `lambda J_0(lambda) / J_1(lambda) = 1 - nu` via `Sms` & kernel (`lambda ~ 2.05`)

mod common;

use femcore::linalg::csr::c;
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::props::Periodic;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_physics::{Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec};
use nyx_sms::*;
use serde_json::json;

const PI: f64 = std::f64::consts::PI;

fn iso_speeds(e: f64, nu: f64, rho: f64) -> (f64, f64) {
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let cp = ((lam + 2.0 * mu) / rho).sqrt();
    let cs = (mu / rho).sqrt();
    (cp, cs)
}

/// Spectral-accuracy Bessel `J_n(x) = (1/2pi) int_0^{2pi} cos(n t - x sin t) dt`.
fn bessel_j(n: i32, x: f64) -> f64 {
    let steps = 400;
    (0..steps)
        .map(|k| {
            let t = 2.0 * PI * k as f64 / steps as f64;
            (n as f64 * t - x * t.sin()).cos()
        })
        .sum::<f64>()
        / steps as f64
}

fn bessel_j_prime(n: i32, x: f64) -> f64 {
    0.5 * (bessel_j(n - 1, x) - bessel_j(n + 1, x))
}

fn bisect_roots(f: impl Fn(f64) -> f64, lo: f64, hi: f64, step: f64) -> Vec<f64> {
    let mut out = Vec::new();
    let mut x = lo;
    let mut fa = f(x);
    while x + step <= hi {
        let b = x + step;
        let fb = f(b);
        if fa.is_finite() && fb.is_finite() && fa * fb <= 0.0 {
            let (mut a, mut b0, mut fa0) = (x, b, fa);
            for _ in 0..64 {
                let mid = 0.5 * (a + b0);
                let fm = f(mid);
                if fa0 * fm <= 0.0 {
                    b0 = mid;
                } else {
                    a = mid;
                    fa0 = fm;
                }
            }
            let r = 0.5 * (a + b0);
            if f(r).abs() < 1e-4 {
                out.push(r);
            }
        }
        x = b;
        fa = fb;
    }
    out
}

fn nearest_rel_err(freqs: &[f64], target: f64) -> f64 {
    freqs.iter().map(|&f| (f - target).abs() / target).fold(f64::INFINITY, f64::min)
}

#[test]
fn cube_3d_traction_free_lame_guided_and_clamped_eigenfrequencies() {
    let (l, e, nu, rho) = (0.1_f64, 200.0e9_f64, 0.25_f64, 7850.0_f64);
    let (cp, cs) = iso_speeds(e, nu, rho);

    // 1. Traction-free cube [0, L]^3:
    //    - 6 rigid-body modes at ~0 Hz
    //    - Exact analytical Lamé equimode (u_x = cos(pi x/L) sin(pi y/L), u_y = -sin(pi x/L) cos(pi y/L), u_z = 0):
    //      f_Lame = sqrt(2) * c_s / (2 * L), independent of Poisson ratio!
    let coarse = common::block(l, l, l, 4, 4, 4);
    let mut p_free = Problem::new(Study::Eigenfrequency { n_modes: 10, f_search: 500.0 }, 2);
    p_free.solver = femcore::props::LinearSolver::Direct;
    p_free.materials.set(1, "E", e);
    p_free.materials.set(1, "nu", nu);
    p_free.materials.set(1, "rho", rho);
    let r_free = common::solve(&coarse, &p_free);
    let f_free: Vec<f64> = r_free.solutions.iter().map(|s| s.value.unwrap()).collect();
    let f_lame_exact = (2.0_f64).sqrt() * cs / (2.0 * l);
    let n_rigid = f_free.iter().filter(|&&f| f < 1e-4 * f_lame_exact).count();
    assert_eq!(n_rigid, 6, "free 3D cube must have 6 rigid modes, got {f_free:?} (log={:?})", r_free.log);

    // Free cube analytical Lamé equimode f_Lame = sqrt(2) * c_s / (2 * L) (sits above lowest torsional/flexural modes):
    let mut p_free_lame = p_free.clone();
    p_free_lame.study = Study::Eigenfrequency { n_modes: 10, f_search: 0.92 * f_lame_exact };
    let r_free_lame = common::solve(&coarse, &p_free_lame);
    let f_free_lame: Vec<f64> = r_free_lame.solutions.iter().map(|s| s.value.unwrap()).collect();
    let err_lame = nearest_rel_err(&f_free_lame, f_lame_exact);
    assert!(
        err_lame < 5e-3,
        "free cube Lamé equimode err={err_lame:.3e} vs exact {f_lame_exact:.4} Hz (modes={f_free_lame:?})"
    );

    // Also verify lowest nonzero torsional/breathing free cube mode against fine P2 reference (< 5e-3).
    let f_free_nonzero_coarse = f_free[6];
    let fine = common::block(l, l, l, 7, 7, 7);
    let r_free_fine = common::solve(&fine, &p_free);
    let f_free_nonzero_fine = r_free_fine.solutions[6].value.unwrap();
    let err_free1 = (f_free_nonzero_coarse - f_free_nonzero_fine).abs() / f_free_nonzero_fine;
    assert!(
        err_free1 < 5e-3,
        "free cube lowest elastic mode coarse={f_free_nonzero_coarse:.4} vs fine={f_free_nonzero_fine:.4} (rel={err_free1:.3e})"
    );

    // 2. Guided-face (roller symmetry) cube [0, L]^3 (u_n = 0, sigma_t = 0 on all 6 faces):
    //    Exact closed-form Navier spectrum:
    //    f_T(1,1,0) = sqrt(2) * c_s / (2 * L), f_L(1,0,0) = c_p / (2 * L), f_T(1,1,1) = sqrt(3) * c_s / (2 * L).
    let mut p_guided = Problem::new(Study::Eigenfrequency { n_modes: 12, f_search: 0.0 }, 2);
    p_guided.materials = p_free.materials.clone();
    p_guided.bcs = vec![
        Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { labels: vec![3, 4], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
    ];
    let r_guided = common::solve(&coarse, &p_guided);
    let f_guided: Vec<f64> = r_guided.solutions.iter().map(|s| s.value.unwrap()).collect();
    for (name, exact) in [
        ("T(1,1,0)", (2.0_f64).sqrt() * cs / (2.0 * l)),
        ("T(1,1,1)", (3.0_f64).sqrt() * cs / (2.0 * l)),
        ("L(1,0,0)", cp / (2.0 * l)),
        ("T(2,1,0)", (5.0_f64).sqrt() * cs / (2.0 * l)),
    ] {
        let err = nearest_rel_err(&f_guided, exact);
        assert!(err < 5e-3, "guided cube {name} rel_err={err:.3e} vs {exact:.4} Hz (modes={f_guided:?})");
    }
    let f_l110_exact = (2.0_f64).sqrt() * cp / (2.0 * l);
    let mut p_guided_hi = p_guided.clone();
    p_guided_hi.study = Study::Eigenfrequency { n_modes: 12, f_search: 0.99 * f_l110_exact };
    let f_guided_hi: Vec<f64> =
        common::solve(&coarse, &p_guided_hi).solutions.iter().map(|s| s.value.unwrap()).collect();
    let err_l110 = nearest_rel_err(&f_guided_hi, f_l110_exact);
    assert!(
        err_l110 < 5e-3,
        "guided cube L(1,1,0) rel_err={err_l110:.3e} vs {f_l110_exact:.4} Hz (modes={f_guided_hi:?})"
    );

    // 3. Fully clamped cube (`u = 0` on all 6 faces): coarse P2 (5x5x5) vs fine P2 (8x8x8) reference `< 5e-3`
    //    plus rigorous Rayleigh-Ritz frequency bracket above the guided lower bound `max(c_p/(2L), sqrt(3)c_s/(2L))`.
    let clamp_coarse = common::block(l, l, l, 5, 5, 5);
    let clamp_fine = common::block(l, l, l, 8, 8, 8);
    let mut p_clamp = Problem::new(Study::Eigenfrequency { n_modes: 6, f_search: 0.0 }, 2);
    p_clamp.materials = p_free.materials.clone();
    p_clamp.bcs = vec![Bc {
        labels: vec![1, 2, 3, 4, 5, 6],
        kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] },
    }];
    let f_clamp_coarse: Vec<f64> =
        common::solve(&clamp_coarse, &p_clamp).solutions.iter().map(|s| s.value.unwrap()).collect();
    let f_clamp_fine: Vec<f64> =
        common::solve(&clamp_fine, &p_clamp).solutions.iter().map(|s| s.value.unwrap()).collect();
    for m in 0..4 {
        let rel = (f_clamp_coarse[m] - f_clamp_fine[m]).abs() / f_clamp_fine[m];
        assert!(
            rel < 5e-3,
            "clamped cube mode {}: coarse={:.4} vs fine={:.4} (rel={rel:.3e})",
            m + 1,
            f_clamp_coarse[m],
            f_clamp_fine[m]
        );
    }
    assert!(
        f_clamp_fine[0] > f_guided[0] && f_clamp_fine[0] > cp / (2.0 * l),
        "clamped fundamental {:.4} must strictly exceed guided bounds",
        f_clamp_fine[0]
    );
}

#[test]
fn kirsch_plate_with_circular_hole_matches_scf_3_with_finite_width_correction() {
    // Quarter-plate symmetry [0, W] x [0, H] with circular hole of radius `a` at origin (0, 0),
    // under uniform far-field tension sigma_0 along +x.
    let (w, h, a, sigma0) = (0.5_f64, 0.5_f64, 0.04_f64, 1.0e6_f64);
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (w, h));
    let hole = m.circle("hole", (0.0, 0.0), a);
    m.overlay_full(&[("plate", plate, RegionKind::Material, None), ("hole", hole, RegionKind::Hole, None)], &[]);
    m.size("hole_zone", Selector::parse("in_box([0, 0, 0], [0.12, 0, 0.12])").unwrap(), Some(0.002.into()), None, None);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.025)).unwrap();

    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Stress;
    s.materials.set(1, "E", 200e9);
    s.materials.set(1, "nu", 0.3);
    s.bcs = vec![
        Bc {
            labels: g.select("plane(x = 0)").unwrap(),
            kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] },
        },
        Bc {
            labels: g.select("plane(z = 0)").unwrap(),
            kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] },
        },
        Bc { labels: g.select(&format!("plane(x = {w})")).unwrap(), kind: BcKind::Load([sigma0, 0.0, 0.0]) },
    ];
    let r = common::solve(&g.mesh, &s);
    let sxx = r.solutions[0].flat("sxx").unwrap();
    let sxx_max = sxx.iter().copied().fold(0.0_f64, f64::max);
    let kt_got = sxx_max / sigma0;

    // Howland / Peterson finite-width strip correction for diameter/(total width) = lam = a/H:
    // K_{t,net} = 3.0 - 3.13 * lam + 3.66 * lam^2 - 1.53 * lam^3, K_{t,gross} = K_{t,net} / (1 - lam).
    let lam = a / h;
    let kt_net = 3.0 - 3.13 * lam + 3.66 * lam * lam - 1.53 * lam.powi(3);
    let kt_corrected = kt_net / (1.0 - lam);

    let err_inf = (kt_got - 3.0).abs() / 3.0;
    let err_corr = (kt_got - kt_corrected).abs() / kt_corrected;
    println!("Kirsch hole: K_t(got)={kt_got:.4}, K_t(corr)={kt_corrected:.4}, err_inf={err_inf:.3e}, err_corr={err_corr:.3e}");
    assert!(err_inf < 0.03, "Kirsch gross SCF {kt_got:.4} deviates >3% from 3.0");
    assert!(err_corr < 0.03, "Kirsch SCF {kt_got:.4} deviates >3% from finite-width {kt_corrected:.4}");

    // Also verify zero/compressive hoop stress at the forward stagnation point (x = a, z = 0) -> exact Kirsch -1.0 * sigma0.
    let coords = g.mesh.dof_coordinates(2);
    let idx_front = (0..g.mesh.n_dofs(2))
        .min_by(|&i, &j| {
            let di = (coords[2 * i] - a).hypot(coords[2 * i + 1]);
            let dj = (coords[2 * j] - a).hypot(coords[2 * j + 1]);
            di.total_cmp(&dj)
        })
        .unwrap();
    let szz = r.solutions[0].flat("szz").unwrap();
    let kt_front = szz[idx_front] / sigma0;
    assert!((kt_front - (-1.0)).abs() < 0.12, "Kirsch stagnation sigma_zz(a,0)/sigma0 = {kt_front:.4} vs -1.0");
}

#[test]
fn rayleigh_lamb_dispersion_for_nonzero_wavevector_matches_analytical() {
    let (h, e, nu, rho) = (1.0e-3_f64, 70.0e9_f64, 0.33_f64, 2700.0_f64);
    let d = 0.5 * h;
    let (cp, cs) = iso_speeds(e, nu, rho);

    // Strip [0, w] x [-d, d] with Periodic (phase = 0) on x = 0 and x = w makes the 2D waveguide
    // cross-section model an infinite plate of thickness 2d = h propagating along y as exp(-i q y).
    let w = 0.1e-3;
    let base = common::rect(w, h, 2, 16);
    let mut pts = base.points.clone();
    for p in pts.chunks_mut(2) {
        p[1] -= d;
    }
    let mesh = SimplexMesh::new(
        2,
        pts,
        base.cells.clone(),
        base.cell_labels.clone(),
        base.facets.clone(),
        base.facet_labels.clone(),
    );

    // Complex-safe normalized Rayleigh-Lamb characteristic functions (real and smooth for all omega > 0,
    // including subsonic omega < q * c_s where ql and qt are purely imaginary):
    //   sinc(z) = sin(z)/z (real for real and purely imaginary z, no singularity at z = 0).
    let sinc_c = |z: femcore::linalg::csr::C| -> femcore::linalg::csr::C {
        if z.norm() < 1e-8 {
            c(1.0, 0.0) - z * z / 6.0
        } else {
            z.sin() / z
        }
    };
    let rl_sym = |omega: f64, q: f64| -> f64 {
        let ql = (c(omega * omega / (cp * cp) - q * q, 0.0)).sqrt();
        let qt = (c(omega * omega / (cs * cs) - q * q, 0.0)).sqrt();
        let q2_qt2 = c(q * q, 0.0) - qt * qt;
        let term1 = q2_qt2 * q2_qt2 * (ql * d).cos() * sinc_c(qt * d);
        let term2 = c(4.0 * q * q, 0.0) * ql * ql * sinc_c(ql * d) * (qt * d).cos();
        (term1 + term2).re / q.powi(4)
    };
    let rl_anti = |omega: f64, q: f64| -> f64 {
        let ql = (c(omega * omega / (cp * cp) - q * q, 0.0)).sqrt();
        let qt = (c(omega * omega / (cs * cs) - q * q, 0.0)).sqrt();
        let q2_qt2 = c(q * q, 0.0) - qt * qt;
        let term1 = q2_qt2 * q2_qt2 * sinc_c(ql * d) * (qt * d).cos();
        let term2 = c(4.0 * q * q, 0.0) * qt * qt * (ql * d).cos() * sinc_c(qt * d);
        (term1 + term2).re / q.powi(4)
    };

    for &qd in &[0.8_f64, 2.0_f64] {
        let q = qd / d;
        let w_max = 3.5 * cs / d;
        let sym_roots = bisect_roots(|w_val| rl_sym(w_val, q), 10.0, w_max, 0.01 * cs / d);
        let anti_roots = bisect_roots(|w_val| rl_anti(w_val, q), 10.0, w_max, 0.01 * cs / d);
        assert!(!sym_roots.is_empty() && !anti_roots.is_empty(), "qd={qd}: sym={sym_roots:?}, anti={anti_roots:?}");

        let mut p = Problem::new(Study::Modes { n_modes: 8, f_search: 0.0, wavevector: q }, 2);
        p.materials.set(1, "E", e);
        p.materials.set(1, "nu", nu);
        p.materials.set(1, "rho", rho);
        p.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 0.0 }];
        let r = common::solve(&mesh, &p);
        let freqs: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();

        for (label, omega_exact) in [("A0", anti_roots[0]), ("S0", sym_roots[0])] {
            let f_exact = omega_exact / (2.0 * PI);
            let err = nearest_rel_err(&freqs, f_exact);
            assert!(
                err < 1e-3,
                "Rayleigh-Lamb q*d={qd} {label}: rel_err={err:.3e} vs {f_exact:.4} Hz (modes={freqs:?})"
            );
        }
        // Also verify lowest SH0 mode f = q * c_s / (2 pi) decoupled along x!
        let f_sh0 = q * cs / (2.0 * PI);
        let err_sh = nearest_rel_err(&freqs, f_sh0);
        assert!(err_sh < 1e-4, "SH0 q*d={qd}: rel_err={err_sh:.3e} vs {f_sh0:.4} Hz");
    }
}

#[test]
fn solid_cylinder_m2_axisymmetric_modes_match_analytical_bessel_roots() {
    // Solid cylinder r in [0, R], z in [-L/2, L/2] with roller ends u_z = 0 (k_z = 0 plane strain)
    // and azimuthal order m = 2. Exact coupled (u_r, u_theta) modes satisfy M11*M22 - M12*M21 = 0
    // where a = omega * R / c_p, b = omega * R / c_s,
    // M11(x) = (m^2 - b^2/2) J_m(x) - x J_m'(x), M12 = m (b J_m'(b) - J_m(b)),
    // M21 = m (a J_m'(a) - J_m(a)), M22 = (m^2 - b^2/2) J_m(b) - b J_m'(b).
    let (e, nu, rho, r_cyl, len) = (169.0e9_f64, 0.28_f64, 2329.0_f64, 1.0_f64, 0.4_f64);
    let (cp, cs) = iso_speeds(e, nu, rho);
    let m_ord = 2_i32;

    let char_det = |b: f64| -> f64 {
        let a = b * (cs / cp);
        let m_f = m_ord as f64;
        let hb2 = 0.5 * b * b;
        let m11 = (m_f * m_f - hb2) * bessel_j(m_ord, a) - a * bessel_j_prime(m_ord, a);
        let m12 = m_f * (b * bessel_j_prime(m_ord, b) - bessel_j(m_ord, b));
        let m21 = m_f * (a * bessel_j_prime(m_ord, a) - bessel_j(m_ord, a));
        let m22 = (m_f * m_f - hb2) * bessel_j(m_ord, b) - b * bessel_j_prime(m_ord, b);
        m11 * m22 - m12 * m21
    };
    let b_roots = bisect_roots(char_det, 0.5, 9.0, 0.02);
    assert!(b_roots.len() >= 3, "expected at least 3 m=2 roots, got {b_roots:?}");

    let raw = common::parts_mesh(&[common::part(0.0, -0.5 * len, r_cyl, 0.5 * len, 1, Some(0.05))], 0.05);
    let end_labels: Vec<i32> = {
        let mut s = std::collections::BTreeSet::new();
        for f in 0..raw.n_facets() {
            let fv = raw.facet(f);
            if fv.iter().all(|&v| (raw.point(v as usize)[1].abs() - 0.5 * len).abs() < 1e-9) {
                s.insert(raw.facet_labels[f]);
            }
        }
        s.into_iter().collect()
    };
    let mut sp = Problem::new(Study::Eigenfrequency { n_modes: 8, f_search: 0.0 }, 2);
    sp.materials.set(1, "E", e);
    sp.materials.set(1, "nu", nu);
    sp.materials.set(1, "rho", rho);
    sp.bcs = vec![Bc { labels: end_labels, kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } }];

    let (_, modes) = solve_solid_axisym(&raw, &sp, 2.0).unwrap();
    let freqs: Vec<f64> = modes.freqs.iter().map(|z| z.re).collect();
    for (idx, &b) in b_roots.iter().take(3).enumerate() {
        let f_exact = b * cs / (2.0 * PI * r_cyl);
        let rel = nearest_rel_err(&freqs, f_exact);
        println!("m=2 solid cylinder mode {}: exact={f_exact:.4} Hz, rel_err={rel:.3e}", idx + 1);
        assert!(rel < 2e-3, "m=2 mode {}: rel_err={rel:.3e} vs {f_exact:.4} Hz (modes={freqs:?})", idx + 1);
    }
}

#[test]
fn free_disk_radial_contour_mode_via_physics_trait_matches_bessel_205() {
    // Free thin disk of radius R, thickness t << R (m = 0):
    // lambda * J_0(lambda) / J_1(lambda) = 1 - nu  =>  lambda_1 ~ 2.04885 (for nu = 0.3),
    // f = lambda * sqrt(E / (rho (1 - nu^2))) / (2 pi R).
    let (e, nu, rho, r_disk, t_disk) = (160.0e9_f64, 0.30_f64, 2330.0_f64, 1.0_f64, 0.02_f64);
    let c_plate = (e / (rho * (1.0 - nu * nu))).sqrt();
    let roots = bisect_roots(|x| x * bessel_j(0, x) / bessel_j(1, x) - (1.0 - nu), 0.5, 8.0, 0.02);
    assert!(roots.len() >= 2 && (roots[0] - 2.04885).abs() < 1e-4, "{roots:?}");

    let simplex = common::parts_mesh(&[common::part(0.0, -0.5 * t_disk, r_disk, 0.5 * t_disk, 1, Some(0.01))], 0.01);
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();

    struct DiskResolver;
    impl LabelResolver for DiskResolver {
        fn domains(&self, name: &str) -> Result<Vec<i32>, Error> {
            match name {
                "disk" => Ok(vec![1]),
                _ => Err(Error::invalid(format!("unknown domain {name}"))),
            }
        }
        fn boundaries(&self, name: &str) -> Result<Vec<i32>, Error> {
            Err(Error::invalid(format!("unknown boundary {name}")))
        }
    }

    let mut spec = ProblemSpec::new("sms", 2, "axisymmetric_eigenfrequency");
    let f1_approx = 2.05 * c_plate / (2.0 * PI * r_disk);
    spec.study.params = json!({ "n_modes": 6, "f_search": 0.9 * f1_approx, "azimuthal_order": 0 });
    spec.materials
        .insert("disk".into(), MaterialAssign::Props([("E".into(), e), ("nu".into(), nu), ("rho".into(), rho)].into()));

    let res = Sms.solve(&mesh, &spec, &mut Ctx::new(&DiskResolver)).unwrap();
    let freqs: Vec<f64> = res.solutions.iter().map(|s| s.value.unwrap()).collect();
    let f1_exact = roots[0] * c_plate / (2.0 * PI * r_disk);
    let rel1 = nearest_rel_err(&freqs, f1_exact);
    assert!(rel1 < 5e-4, "free-disk contour mode 1: rel_err={rel1:.3e} vs {f1_exact:.4} Hz (modes={freqs:?})");
}
