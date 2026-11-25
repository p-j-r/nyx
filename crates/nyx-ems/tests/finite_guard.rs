// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! NaN or infinite input must be refused with a message that names the culprit, and a solve whose
//! result contains non-finite numbers must not be returned as a success.

mod common;

use common::{emw, solve};
use femcore::material::{Anisotropy, Tensor};
use femcore::props::{Periodic, Pml, PmlCoord, PmlDir, PmlProfile, PmlStabilization};
use nyx_ems::{LumpedJjBc, LumpedPortBc, Problem, Study, SurfaceImpedanceBc, WavePortBc};
use nyx_physics::{check_result_finite, Field, Solution, SolveResult, TlsInterfaceSpec, TlsInterfaces};

const NON_FINITES: [f64; 3] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

fn valid_emw_problem() -> Problem {
    let mut p = emw();
    p.materials.default.insert("eps_r".into(), 1.0);
    p.materials.default.insert("mu_r".into(), 1.0);
    p
}

#[test]
fn a_valid_problem_still_solves() {
    let m = common::rect(0.02286, 0.01016, 6, 3);
    assert!(solve(&m, &valid_emw_problem()).is_ok());

    // Quoted strings containing "inf" / "NaN" must not trigger a false positive in check_finite.
    let mut quoted = valid_emw_problem();
    quoted.surface_impedance.push(SurfaceImpedanceBc {
        labels: vec![1],
        l_k_sq: 1e-12,
        r_s: 0.0,
        material: "infinity_NaN_nb".into(),
        thickness: 100e-9,
        temperature: 0.02,
        f_ref: 5e9,
        film_model: "thin_limit".into(),
        lossless: true,
        ibc_mode: "direct".into(),
    });
    quoted.lumped_ports.push(LumpedPortBc {
        labels: vec![2],
        name: "infinite_nan_port".into(),
        z0: 50.0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: false,
        mode: "single".into(),
    });
    quoted.lumped_jj.push(LumpedJjBc {
        labels: vec![3],
        name: "nan_jj_inf".into(),
        l_j: 10e-9,
        e_j_hz: 16e9,
        c_j: 0.0,
        e_l_hz: 0.0,
        phi_ext_rad: 0.0,
        n_g: 0.0,
        squid_asymmetry: None,
        direction: [0.0, 0.0, 1.0],
        formulation: "rank1".into(),
    });
    quoted.waveports.push(WavePortBc {
        labels: vec![4],
        name: "inf_waveport".into(),
        mode: 1,
        n_eff_guess: Some(1.0),
        z0_ohm: Some(50.0),
    });
    assert!(quoted.check_finite().is_ok());
}

#[test]
fn nan_and_infinite_materials_are_refused_by_name() {
    let m = common::rect(0.02286, 0.01016, 6, 3);
    for key in ["eps_r", "eps_r_im", "mu_r", "sigma"] {
        for bad in NON_FINITES {
            let mut p = valid_emw_problem();
            p.materials.default.insert(key.into(), bad);
            let e = solve(&m, &p).unwrap_err();
            assert!(e.contains(&format!("'{key}'")), "expected '{key}' in error: {e}");

            let mut pd = valid_emw_problem();
            pd.materials.set(1, key, bad);
            let ed = solve(&m, &pd).unwrap_err();
            assert!(ed.contains("domain 1") && ed.contains(&format!("'{key}'")), "{ed}");
        }
    }

    for bad in NON_FINITES {
        let mut p_aniso = valid_emw_problem();
        p_aniso.materials.aniso.insert(
            1,
            Anisotropy { permittivity: Some(Tensor::Diagonal { x: bad, y: 1.0, z: 1.0 }), ..Default::default() },
        );
        assert!(solve(&m, &p_aniso).is_err());

        let mut p_im = valid_emw_problem();
        p_im.materials.aniso.insert(
            1,
            Anisotropy { permittivity_im: Some(Tensor::Uniaxial { perp: bad, par: 0.0 }), ..Default::default() },
        );
        assert!(solve(&m, &p_im).is_err());

        let mut p_euler = valid_emw_problem();
        p_euler.materials.aniso.insert(1, Anisotropy { euler_deg: [bad, 0.0, 0.0], ..Default::default() });
        assert!(solve(&m, &p_euler).is_err());

        for pml in [
            Pml { x: PmlDir::Plus, strength: bad, ..Default::default() },
            Pml {
                coord: PmlCoord::Spherical { center: [bad, 0.0, 0.0], r_in: 0.01, thickness: 0.005 },
                ..Default::default()
            },
            Pml { profile: PmlProfile::Rational { curvature: bad }, ..Default::default() },
            Pml { stabilization: PmlStabilization::Mpml { ratio: bad }, ..Default::default() },
            Pml { bounds: Some(([0.0; 3], [bad, 1.0, 1.0])), ..Default::default() },
        ] {
            let mut p_pml = valid_emw_problem();
            p_pml.materials.pml.insert(1, pml);
            assert!(solve(&m, &p_pml).is_err());
        }

        let mut p_per = valid_emw_problem();
        p_per.periodic.push(Periodic { source: vec![1], dest: vec![4], phase: bad });
        assert!(solve(&m, &p_per).is_err());

        let mut p_disp = valid_emw_problem();
        p_disp.d_omega_eps.insert(1, [[bad, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(solve(&m, &p_disp).is_err());

        for spec in [
            TlsInterfaceSpec { thickness_m: bad, eps_r: 10.0, tan_delta: 1e-3 },
            TlsInterfaceSpec { thickness_m: 2e-9, eps_r: bad, tan_delta: 1e-3 },
            TlsInterfaceSpec { thickness_m: 2e-9, eps_r: 10.0, tan_delta: bad },
        ] {
            let mut p_tls = valid_emw_problem();
            p_tls.tls_interfaces = Some(TlsInterfaces { ma: None, ms: Some(spec), sa: None });
            assert!(solve(&m, &p_tls).is_err());
        }
    }
}

#[test]
fn nan_geometry_is_refused() {
    for bad in NON_FINITES {
        let mut m = common::rect(0.02286, 0.01016, 6, 3);
        m.points[3] = bad;
        let e = solve(&m, &valid_emw_problem()).unwrap_err();
        assert!(e.contains("non-finite coordinate"), "{e}");
    }
}

#[test]
fn nan_and_infinite_study_parameters_are_refused() {
    let m = common::rect(0.02286, 0.01016, 6, 3);
    for bad in NON_FINITES {
        let studies = [
            Study::Modes { freq: bad, n_target: 1.0, evanescent: false, n_modes: 2 },
            Study::Modes { freq: 10e9, n_target: bad, evanescent: false, n_modes: 2 },
            Study::Bloch { freq: bad, n_target: 1.0, n_modes: 2, brillouin_zone: 0 },
            Study::Bloch { freq: 10e9, n_target: bad, n_modes: 2, brillouin_zone: 0 },
            Study::AxisymmetricModes { f_search: bad, azimuthal_order: 0, n_modes: 2 },
            Study::Cavity { n_modes: 2, f_search: bad },
            Study::Driven {
                f_min: bad,
                f_max: 10e9,
                n_points: 5,
                freqs_hz: None,
                prom_tol: None,
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            Study::Driven {
                f_min: 5e9,
                f_max: bad,
                n_points: 5,
                freqs_hz: None,
                prom_tol: None,
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            Study::Driven {
                f_min: 5e9,
                f_max: 10e9,
                n_points: 2,
                freqs_hz: Some(vec![5e9, bad]),
                prom_tol: None,
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            Study::Driven {
                f_min: 5e9,
                f_max: 10e9,
                n_points: 5,
                freqs_hz: None,
                prom_tol: Some(bad),
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            Study::Driven {
                f_min: 5e9,
                f_max: 10e9,
                n_points: 5,
                freqs_hz: None,
                prom_tol: None,
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![bad],
                solver: "auto".into(),
                blr_tol: 0.0,
            },
            Study::Driven {
                f_min: 5e9,
                f_max: 10e9,
                n_points: 5,
                freqs_hz: None,
                prom_tol: None,
                max_hdm_solves: 12,
                strict_prom_tol: false,
                seed_freqs_hz: vec![],
                solver: "auto".into(),
                blr_tol: bad,
            },
        ];

        for study in studies {
            let mut p = valid_emw_problem();
            p.study = study;
            assert!(solve(&m, &p).is_err());
        }
    }
}

#[test]
fn nan_and_infinite_boundary_conditions_are_refused() {
    let m = common::rect(0.02286, 0.01016, 6, 3);
    for bad in NON_FINITES {
        let ibc_cases = [
            SurfaceImpedanceBc {
                labels: vec![1],
                l_k_sq: bad,
                r_s: 0.0,
                material: String::new(),
                thickness: 100e-9,
                temperature: 0.02,
                f_ref: 5e9,
                film_model: "thin_limit".into(),
                lossless: true,
                ibc_mode: "direct".into(),
            },
            SurfaceImpedanceBc {
                labels: vec![1],
                l_k_sq: 1e-12,
                r_s: bad,
                material: String::new(),
                thickness: 100e-9,
                temperature: 0.02,
                f_ref: 5e9,
                film_model: "thin_limit".into(),
                lossless: false,
                ibc_mode: "direct".into(),
            },
            SurfaceImpedanceBc {
                labels: vec![1],
                l_k_sq: 1e-12,
                r_s: 0.0,
                material: String::new(),
                thickness: bad,
                temperature: 0.02,
                f_ref: 5e9,
                film_model: "thin_limit".into(),
                lossless: true,
                ibc_mode: "direct".into(),
            },
            SurfaceImpedanceBc {
                labels: vec![1],
                l_k_sq: 1e-12,
                r_s: 0.0,
                material: String::new(),
                thickness: 100e-9,
                temperature: bad,
                f_ref: 5e9,
                film_model: "thin_limit".into(),
                lossless: true,
                ibc_mode: "direct".into(),
            },
            SurfaceImpedanceBc {
                labels: vec![1],
                l_k_sq: 1e-12,
                r_s: 0.0,
                material: String::new(),
                thickness: 100e-9,
                temperature: 0.02,
                f_ref: bad,
                film_model: "thin_limit".into(),
                lossless: true,
                ibc_mode: "direct".into(),
            },
        ];
        for ibc in ibc_cases {
            let mut p = valid_emw_problem();
            p.surface_impedance.push(ibc);
            assert!(solve(&m, &p).is_err());
        }

        let jj_base = LumpedJjBc {
            labels: vec![2],
            name: "jj1".into(),
            l_j: 10e-9,
            e_j_hz: 16e9,
            c_j: 0.0,
            e_l_hz: 0.0,
            phi_ext_rad: 0.0,
            n_g: 0.0,
            squid_asymmetry: None,
            direction: [0.0, 0.0, 1.0],
            formulation: "rank1".into(),
        };
        let jj_cases = [
            LumpedJjBc { l_j: bad, ..jj_base.clone() },
            LumpedJjBc { e_j_hz: bad, ..jj_base.clone() },
            LumpedJjBc { c_j: bad, ..jj_base.clone() },
            LumpedJjBc { e_l_hz: bad, ..jj_base.clone() },
            LumpedJjBc { phi_ext_rad: bad, ..jj_base.clone() },
            LumpedJjBc { n_g: bad, ..jj_base.clone() },
            LumpedJjBc { squid_asymmetry: Some(bad), ..jj_base.clone() },
            LumpedJjBc { direction: [bad, 0.0, 0.0], ..jj_base.clone() },
        ];
        for jj in jj_cases {
            let mut p = valid_emw_problem();
            p.lumped_jj.push(jj);
            assert!(solve(&m, &p).is_err());
        }

        for lp in [
            LumpedPortBc {
                labels: vec![2],
                name: "p1".into(),
                z0: bad,
                direction: [0.0, 0.0, 1.0],
                complex_damping: false,
                mode: "single".into(),
            },
            LumpedPortBc {
                labels: vec![2],
                name: "p1".into(),
                z0: 50.0,
                direction: [0.0, bad, 0.0],
                complex_damping: false,
                mode: "single".into(),
            },
        ] {
            let mut p = valid_emw_problem();
            p.lumped_ports.push(lp);
            assert!(solve(&m, &p).is_err());
        }

        for wp in [
            WavePortBc { labels: vec![3], name: "wp1".into(), mode: 1, n_eff_guess: Some(bad), z0_ohm: None },
            WavePortBc { labels: vec![3], name: "wp1".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(bad) },
        ] {
            let mut p = valid_emw_problem();
            p.waveports.push(wp);
            assert!(solve(&m, &p).is_err());
        }
    }
}

#[test]
fn a_non_finite_result_is_an_error() {
    let mut r = SolveResult::default();
    r.solutions.push(Solution {
        label: "mode 1".into(),
        fields: vec![Field::scalar("Ex", "V/m", vec![0.0, f64::NAN, 1.0])],
        ..Default::default()
    });
    let e = check_result_finite(&r).unwrap_err().to_string();
    assert!(e.contains("'Ex'") && e.contains("1 NaN"), "{e}");
    r.solutions[0].fields[0].values[1] = 0.5;
    assert!(check_result_finite(&r).is_ok());
}
