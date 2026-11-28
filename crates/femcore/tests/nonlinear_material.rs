// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification of `[nonlinear]` TOML v2 material definitions (`Chi2`, `Chi3`, `QpmSpec`, ADR 143):
//! - Neumann `PointGroup` centrosymmetry & crystal-class selection-rule enforcement,
//! - Bunge Z-X-Z Euler rotation of `LiNbO3` (`3m`) contracted $d_{iJ}$ tensor (including $120^\circ$
//!   3-fold invariance and X-cut / Z-cut axis permutations),
//! - Periodic quasi-phase-matching (`QpmSpec::poling_sign_at`) sign switching and duty cycle,
//! - Vectorial Maker–Terhune $\chi^{(3)}$ (`A`, `B`, $2/3$ circular-to-linear ratio, cubic `m-3m` anisotropy),
//! - Shared time-domain nonlinear polarization & exponential carrier-lifetime evaluators,
//! - Curated `[nonlinear]` sections in `lithium_niobate`, `silicon_nitride`, and `silicon`.

use femcore::material::*;
use num_complex::Complex64;
use std::path::Path;

fn load_repo_library() -> Library {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    Library::load_dir(&dir).unwrap_or_else(|errs| {
        panic!(
            "failed to load materials library:\n{}",
            errs.iter().map(ToString::to_string).collect::<Vec<_>>().join("\n")
        )
    })
}

#[test]
fn centrosymmetric_point_groups_reject_chi2_unless_explicitly_broken() {
    for pg in ["m-3m", "isotropic", "-3m", "6/mmm", "mmm", "4/mmm", "2/m", "-1", "432"] {
        let toml_bad = format!(
            r#"
format = 2
name = "centro_test"
point_group = "{pg}"
source = {{ cite = "APPROXIMATE test" }}
[nonlinear.chi2]
unit = "pm/V"
d33 = 10.0
source = {{ cite = "APPROXIMATE test" }}
"#
        );
        let err = MaterialSpec::from_toml(&toml_bad)
            .expect_err(&format!("point_group {pg} must reject bulk chi2"))
            .to_string();
        assert!(
            err.contains("forbids bulk second-order susceptibility chi2"),
            "unexpected error for {pg}: {err}"
        );
    }

    // Electric-field-induced / strain-poled broken centrosymmetry opt-in succeeds cleanly.
    let toml_efish = r#"
format = 2
name = "poled_sin"
point_group = "isotropic"
source = { cite = "APPROXIMATE optically poled Si3N4 test" }
[nonlinear]
broken_centrosymmetry = true
[nonlinear.chi2]
unit = "pm/V"
d33 = 0.3
d31 = 0.1
d15 = 0.1
source = { cite = "APPROXIMATE photo-galvanic poled Si3N4" }
"#;
    let spec = MaterialSpec::from_toml(toml_efish).expect("broken_centrosymmetry opt-in must succeed");
    let c2 = spec.nonlinear.as_ref().unwrap().chi2.as_ref().unwrap();
    assert!(c2.broken_centrosymmetry);
    assert!((c2.d[2][2] - 0.3e-12).abs() < 1e-24);
}

#[test]
fn linbo3_3m_expansion_threefold_invariance_and_euler_zxz_rotations() {
    let lib = load_repo_library();
    let ln = lib.get("lithium_niobate").expect("lithium_niobate in library");
    let nl = ln.nonlinear.as_ref().expect("lithium_niobate must carry [nonlinear]");
    let c2 = nl.chi2.as_ref().expect("lithium_niobate must carry [nonlinear.chi2]");

    // Verify 3m symmetry expansion: d33=-25.2 pm/V, d31=d32=-4.6 pm/V, d15=d24=-4.6 pm/V, d22=2.1 pm/V, d21=d16=-2.1 pm/V.
    assert!((c2.d[2][2] - (-25.2e-12)).abs() < 1e-24);
    assert!((c2.d[2][0] - (-4.6e-12)).abs() < 1e-24);
    assert!((c2.d[2][1] - (-4.6e-12)).abs() < 1e-24);
    assert!((c2.d[0][4] - (-4.6e-12)).abs() < 1e-24);
    assert!((c2.d[1][3] - (-4.6e-12)).abs() < 1e-24);
    assert!((c2.d[1][1] - 2.1e-12).abs() < 1e-24);
    assert!((c2.d[1][0] - (-2.1e-12)).abs() < 1e-24);
    assert!((c2.d[0][5] - (-2.1e-12)).abs() < 1e-24);

    // 1. Trigonal 3m MUST be invariant under +/-120 degree rotation about crystal Z ([120, 0, 0] and [240, 0, 0]),
    // which specifically tests that contracted optical d_iJ has NO spurious 1/2 shear factor (d16 = -d22, NOT -2 d22)!
    let d0 = c2.d_global([0.0, 0.0, 0.0]);
    for phi_z in [120.0, -120.0, 240.0] {
        let d_rot = c2.d_global([phi_z, 0.0, 0.0]);
        for i in 0..3 {
            for j in 0..6 {
                assert!(
                    (d_rot[i][j] - d0[i][j]).abs() < 1e-24,
                    "3m 120-deg invariance failed at d{}{ } for phi={phi_z}: {} vs {}",
                    i + 1,
                    j + 1,
                    d_rot[i][j],
                    d0[i][j]
                );
            }
        }
    }
    // Conversely, a 60-degree rotation flips the mirror-symmetric d22 component: d22 -> -d22 while keeping d33, d31, d15 unchanged!
    let d60 = c2.d_global([60.0, 0.0, 0.0]);
    assert!((d60[1][1] + d0[1][1]).abs() < 1e-24);
    assert!((d60[2][2] - d0[2][2]).abs() < 1e-24);

    // 2. X-cut TFLN with optic Z-axis rotated into global lateral +X axis so TE modes (Ex) access d33:
    // Bunge Z-X-Z euler = [90.0, 90.0, 0.0]:
    // R = Rz(90) * Rx(90):
    // col 0 (cryst X) = [0, 1, 0], col 1 (cryst Y) = [0, 0, 1], col 2 (cryst Z) = [1, 0, 0]!
    // So crystal Z maps to +X_global => global d_11 (d_xxx) MUST equal crystal d_33!
    let d_xcut = c2.d_global([90.0, 90.0, 0.0]);
    assert!(
        (d_xcut[0][0] - d0[2][2]).abs() < 1e-23,
        "X-cut (cryst Z -> +X_glob) must map d33 -> d11: got {:e} vs {:e}",
        d_xcut[0][0],
        d0[2][2]
    );
    // And crystal X->+Y, Y->+Z, Z->+X means P_x = 2 eps0 (d33 Ex^2 + d31 Ey^2 + d32 Ez^2 + ...), so d_12 == d31 and d_13 == d32!
    assert!((d_xcut[0][1] - d0[2][0]).abs() < 1e-23);
    assert!((d_xcut[0][2] - d0[2][1]).abs() < 1e-23);

    // 3. Round-trip TOML serialization preserves d_iJ, T_range, and chi3 exactly.
    let roundtrip = MaterialSpec::from_toml(&ln.to_toml()).expect("LiNbO3 roundtrip TOML");
    assert_eq!(roundtrip.nonlinear, ln.nonlinear);
}

#[test]
fn qpm_periodic_poling_sign_duty_cycle_and_arbitrary_axis() {
    let toml_qpm = r#"
format = 2
name = "ppln_test"
point_group = "3m"
source = { cite = "APPROXIMATE PPLN test" }
[nonlinear.chi2]
unit = "pm/V"
d33 = -25.2
source = { cite = "APPROXIMATE PPLN test" }
[nonlinear.qpm]
sign = 1.0
Lambda_qpm = "4.0 um"
axis = [0.0, 1.0, 0.0]
duty_cycle = 0.5
phase_offset = "0.0 um"
"#;
    let spec = MaterialSpec::from_toml(toml_qpm).expect("parse PPLN");
    let qpm = spec.nonlinear.as_ref().unwrap().qpm.as_ref().unwrap();
    assert!((qpm.period.unwrap() - 4.0e-6).abs() < 1e-18);

    // Period = 4 um along +y, duty_cycle = 0.5 => [0, 2) um -> +1, [2, 4) um -> -1, periodic for negative y too!
    assert_eq!(qpm.poling_sign_at([0.0, 0.5e-6, 0.0]), 1.0);
    assert_eq!(qpm.poling_sign_at([0.0, 1.99e-6, 0.0]), 1.0);
    assert_eq!(qpm.poling_sign_at([0.0, 2.01e-6, 0.0]), -1.0);
    assert_eq!(qpm.poling_sign_at([0.0, 3.99e-6, 0.0]), -1.0);
    assert_eq!(qpm.poling_sign_at([0.0, 4.5e-6, 0.0]), 1.0);
    assert_eq!(qpm.poling_sign_at([0.0, -0.5e-6, 0.0]), -1.0);
    assert_eq!(qpm.poling_sign_at([0.0, -3.5e-6, 0.0]), 1.0);

    // Non-50% duty cycle (D = 0.25) and inverted base sign (-1)
    let qpm2 = QpmSpec {
        sign: -1.0,
        period: Some(10.0e-6),
        axis: [3.0 / 5.0, 4.0 / 5.0, 0.0],
        duty_cycle: 0.25,
        phase_offset: 0.0,
    };
    qpm2.validate().unwrap();
    // Along unit vector (0.6, 0.8, 0): s = 1.0 um < 2.5 um -> -1.0; s = 5.0 um >= 2.5 um -> +1.0
    assert_eq!(qpm2.poling_sign_at([0.6e-6, 0.8e-6, 0.0]), -1.0);
    assert_eq!(qpm2.poling_sign_at([3.0e-6, 4.0e-6, 0.0]), 1.0);
}

#[test]
fn chi3_maker_terhune_circular_polarization_and_time_domain_evaluators() {
    let lib = load_repo_library();
    let sin = lib.get("silicon_nitride").unwrap();
    let c3_sin = sin.nonlinear.as_ref().unwrap().chi3.as_ref().unwrap();
    assert!((c3_sin.n_2.unwrap() - 2.4e-19).abs() < 1e-30);

    let n0 = 2.0;
    let e0 = 1.0e7; // 10 MV/m peak phasor

    // Linear polarization E = [e0, 0, 0]: P^(3)_x = eps0 * (A + 0.5 B) * e0^3 = eps0 * c0 * eps0 * n0^2 * n2 * e0^3
    let p_lin = eval_chi3_envelope_time(
        c3_sin,
        n0,
        [0.0; 3],
        [Complex64::new(e0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
    );
    let expected_lin = EPSILON_0 * (C0_NL * EPSILON_0 * n0 * n0 * c3_sin.n_2.unwrap()) * e0.powi(3);
    assert!((p_lin[0].re - expected_lin).abs() / expected_lin < 1e-13);
    assert!(p_lin[0].im.abs() < 1e-30);

    // Circular polarization E = [e0, i e0, 0] / sqrt(2) (same total intensity I!):
    // E . E = 0 => the 0.5 * B * (E . E) E^* term VANISHES identically!
    // Under electronic Kleinman symmetry (A = B), circular polarization experiences 2/3 of the linear Kerr shift!
    let s = e0 / 2.0_f64.sqrt();
    let e_circ = [Complex64::new(s, 0.0), Complex64::new(0.0, s), Complex64::new(0.0, 0.0)];
    let p_circ = eval_chi3_envelope_time(c3_sin, n0, [0.0; 3], e_circ);
    let ratio = (p_circ[0].re / s) / (p_lin[0].re / e0);
    assert!((ratio - 2.0 / 3.0).abs() < 1e-13, "Maker-Terhune circular/linear Kerr ratio must be 2/3, got {ratio}");

    // Instantaneous time-domain evaluator matches cycle peak (4/3 of monochromatic 3/4 factor since cos^3(wt) = 3/4 cos(wt) + 1/4 cos(3wt)):
    let p_td = eval_chi3_kerr_time(c3_sin, n0, [0.0; 3], [e0, 0.0, 0.0]);
    assert!(((0.75 * p_td[0]) - p_lin[0].re).abs() / expected_lin < 1e-13);

    // Instantaneous chi2 tangent matches central finite difference of eval_chi2_polarization_time to < 1e-10:
    let ln = lib.get("lithium_niobate").unwrap();
    let d_glob = ln.nonlinear.as_ref().unwrap().chi2.as_ref().unwrap().d_global([25.0, 40.0, 15.0]);
    let e_pt = [1.2e6, -0.8e6, 2.3e6];
    let jac = eval_chi2_tangent_time(&d_glob, -1.0, e_pt);
    let h = 100.0; // 100 V/m step (quadratic P^(2) has zero 3rd derivative, avoiding float cancellation on e_pt ~ 2e6)
    for m in 0..3 {
        let mut ep = e_pt;
        let mut em = e_pt;
        ep[m] += h;
        em[m] -= h;
        let pp = eval_chi2_polarization_time(&d_glob, -1.0, ep);
        let pm = eval_chi2_polarization_time(&d_glob, -1.0, em);
        for i in 0..3 {
            let fd = (pp[i] - pm[i]) / (2.0 * h * EPSILON_0);
            assert!((jac[i][m] - fd).abs() < 1e-10 * jac[i][m].abs().max(1e-12), "jac[{i}][{m}]: {} vs {fd}", jac[i][m]);
        }
    }

    // Silicon TPA + FCA/FCD steady-state carrier density matches long-time limit of drude_carrier_step:
    let si = lib.get("silicon").unwrap();
    let c3_si = si.nonlinear.as_ref().unwrap().chi3.as_ref().unwrap();
    let omega = 2.0 * std::f64::consts::PI * C0_NL / 1.55e-6;
    let n_si = 3.4757;
    let nc_steady = c3_si.steady_carrier_density(
        omega,
        n_si,
        &[Complex64::new(e0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
    );
    assert!(nc_steady > 0.0);
    let intensity = 0.5 * C0_NL * EPSILON_0 * n_si * e0 * e0;
    let mut nc = 0.0;
    let dt = 0.1e-9; // 0.1 ns step (tau_c = 1.0 ns)
    for _ in 0..300 {
        nc = drude_carrier_step(
            nc,
            dt,
            omega,
            c3_si.beta_tpa.unwrap(),
            c3_si.tau_c.unwrap(),
            intensity,
        );
    }
    assert!((nc - nc_steady).abs() / nc_steady < 1e-11, "carrier step steady limit: {nc} vs {nc_steady}");
}
