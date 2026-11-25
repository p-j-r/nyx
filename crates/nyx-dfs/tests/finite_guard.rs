// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! NaN or infinite input must be refused with a message that names the culprit, and a solve whose
//! result contains non-finite numbers must not be returned as a success.

mod common;

use common::try_solve;
use femcore::material::{Anisotropy, Tensor};
use femcore::mesh::SimplexMesh;
use femcore::props::{Pml, PmlCoord, PmlDir, PmlProfile, PmlStabilization};
use nyx_dfs::{Bc, BcKind, InfiniteElementSpec, Problem, Study};
use nyx_physics::{check_result_finite, Field, Solution, SolveResult, TlsInterfaceSpec, TlsInterfaces};

const NON_FINITES: [f64; 3] = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];

fn electro_spec(mesh: &SimplexMesh) -> Problem {
    let mut p = Problem::new(Study::Electrostatic, 1);
    p.materials.default.insert("eps_r".into(), 1.0);
    p.materials.default.insert("rho_v".into(), 0.0);
    let labels = mesh.boundary_labels();
    p.bcs = vec![
        Bc { labels: vec![labels[0]], kind: BcKind::Ground },
        Bc { labels: vec![*labels.last().unwrap()], kind: BcKind::Potential(1.0) },
    ];
    p
}

#[test]
fn a_valid_problem_still_solves() {
    let m = common::rect(1.0, 0.5, 4, 2);
    assert!(try_solve(&m, &electro_spec(&m)).is_ok());

    // Quoted strings containing "inf" / "NaN" must not trigger a false positive in check_finite.
    let mut quoted = electro_spec(&m);
    quoted.bcs.push(Bc { labels: vec![m.boundary_labels()[1]], kind: BcKind::Terminal("infinity_NaN_pad".into()) });
    assert!(quoted.check_finite().is_ok());
}

#[test]
fn nan_and_infinite_materials_are_refused_by_name() {
    let m = common::rect(1.0, 0.5, 4, 2);
    for key in ["eps_r", "eps_r_im", "rho_v", "k", "rho", "Cp", "Q", "mu_r"] {
        for bad in NON_FINITES {
            let mut p = electro_spec(&m);
            p.materials.default.insert(key.into(), bad);
            let e = try_solve(&m, &p).unwrap_err();
            assert!(e.contains(&format!("'{key}'")), "expected '{key}' in error: {e}");

            let mut pd = electro_spec(&m);
            pd.materials.set(1, key, bad);
            let ed = try_solve(&m, &pd).unwrap_err();
            assert!(ed.contains("domain 1") && ed.contains(&format!("'{key}'")), "{ed}");
        }
    }

    for bad in NON_FINITES {
        let mut p_aniso = electro_spec(&m);
        p_aniso.materials.aniso.insert(
            1,
            Anisotropy { permittivity: Some(Tensor::Diagonal { x: bad, y: 1.0, z: 1.0 }), ..Default::default() },
        );
        assert!(try_solve(&m, &p_aniso).is_err());

        let mut p_euler = electro_spec(&m);
        p_euler.materials.aniso.insert(1, Anisotropy { euler_deg: [0.0, bad, 0.0], ..Default::default() });
        assert!(try_solve(&m, &p_euler).is_err());

        let mut p_pml = electro_spec(&m);
        p_pml.materials.pml.insert(
            1,
            Pml {
                x: PmlDir::Plus,
                strength: bad,
                coord: PmlCoord::Cartesian,
                profile: PmlProfile::KappaPolynomial { kappa_max: 1.0 },
                stabilization: PmlStabilization::Auto,
                ..Default::default()
            },
        );
        assert!(try_solve(&m, &p_pml).is_err());

        for spec in [
            TlsInterfaceSpec { thickness_m: bad, eps_r: 10.0, tan_delta: 1e-3 },
            TlsInterfaceSpec { thickness_m: 2e-9, eps_r: bad, tan_delta: 1e-3 },
            TlsInterfaceSpec { thickness_m: 2e-9, eps_r: 10.0, tan_delta: bad },
        ] {
            let mut p_tls = electro_spec(&m);
            p_tls.tls_interfaces = Some(TlsInterfaces { ma: Some(spec), ms: None, sa: None });
            assert!(try_solve(&m, &p_tls).is_err());
        }
    }
}

#[test]
fn nan_geometry_is_refused() {
    for bad in NON_FINITES {
        let mut m = common::rect(1.0, 0.5, 4, 2);
        m.points[3] = bad;
        let e = try_solve(&m, &electro_spec(&m)).unwrap_err();
        assert!(e.contains("non-finite coordinate"), "{e}");
    }
}

#[test]
fn nan_and_infinite_study_and_problem_parameters_are_refused() {
    let m = common::rect(1.0, 0.5, 4, 2);
    for bad in NON_FINITES {
        let mut p_tol = electro_spec(&m);
        p_tol.tol = bad;
        let e = try_solve(&m, &p_tol).unwrap_err();
        assert!(e.contains("'tol'"), "{e}");

        // Check freq, adapt_rel_tol, adapt_growth on both Electrostatic and non-Electrostatic studies.
        for study in [Study::Electrostatic, Study::HeatStationary] {
            let mut p_freq = Problem::new(study.clone(), 1);
            p_freq.freq = bad;
            assert!(p_freq.check_finite().is_err());

            let mut p_rel = Problem::new(study.clone(), 1);
            p_rel.adapt_rel_tol = bad;
            assert!(p_rel.check_finite().is_err());

            let mut p_gr = Problem::new(study, 1);
            p_gr.adapt_growth = bad;
            assert!(p_gr.check_finite().is_err());
        }

        let mut p_tinf = electro_spec(&m);
        p_tinf.farfield_t_inf = bad;
        assert!(try_solve(&m, &p_tinf).is_err());

        let mut p_thk = electro_spec(&m);
        p_thk.sheet_thickness_m.insert(1, bad);
        assert!(try_solve(&m, &p_thk).is_err());

        let mut p_etch = electro_spec(&m);
        p_etch.sheet_etch_depth_m.insert(1, bad);
        assert!(try_solve(&m, &p_etch).is_err());

        let mut p_heat = electro_spec(&m);
        p_heat.study = Study::HeatTransient { t_end: bad, n_out: 4, t0: "300".into() };
        assert!(p_heat.check_finite().is_err());
    }
}

#[test]
fn nan_and_infinite_boundary_conditions_are_refused() {
    let m = common::rect(1.0, 0.5, 4, 2);
    let lbl = *m.boundary_labels().last().unwrap();

    for bad in NON_FINITES {
        let electro_bcs = [
            BcKind::Potential(bad),
            BcKind::SurfaceCharge(bad),
            BcKind::Open { origin: Some([bad, 0.0, 0.0]) },
            BcKind::FarField { order: 1, decay_power: bad, center: None },
            BcKind::FarField { order: 1, decay_power: 1.0, center: Some([0.0, bad, 0.0]) },
            BcKind::InfiniteElement(InfiniteElementSpec { decay: bad, ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { axis: [0.0, 0.0, bad], ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { center: Some([bad, 0.0, 0.0]), ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { r_in: Some(bad), ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { thickness: Some(bad), ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { pole_dist: Some(bad), ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { bounds_lo: Some([bad, 0.0, 0.0]), ..Default::default() }),
            BcKind::InfiniteElement(InfiniteElementSpec { bounds_hi: Some([0.0, 0.0, bad]), ..Default::default() }),
        ];
        for kind in electro_bcs {
            let mut p = electro_spec(&m);
            p.bcs[1] = Bc { labels: vec![lbl], kind };
            let e = try_solve(&m, &p).unwrap_err();
            assert!(e.contains("NaN or infinite"), "{e}");
        }

        let other_bcs = [
            BcKind::Temperature(bad),
            BcKind::Flux(bad),
            BcKind::Convection { h: bad, text: 300.0 },
            BcKind::Convection { h: 10.0, text: bad },
            BcKind::LondonSheet { l_k_sq: bad },
            BcKind::CurrentTerminal { name: "src".into(), current: bad, direction: [1.0, 0.0, 0.0] },
            BcKind::CurrentTerminal { name: "src".into(), current: 1.0, direction: [0.0, bad, 0.0] },
            BcKind::FluxLoop { name: "loop1".into(), flux_wb: bad },
        ];
        for kind in other_bcs {
            let mut p = electro_spec(&m);
            p.bcs[1] = Bc { labels: vec![lbl], kind };
            let e = p.check_finite().unwrap_err();
            assert!(e.contains("NaN or infinite"), "{e}");
        }
    }
}

#[test]
fn a_non_finite_result_is_an_error() {
    let mut r = SolveResult::default();
    r.solutions.push(Solution {
        label: "electrostatic".into(),
        fields: vec![Field::scalar("V", "V", vec![0.0, f64::NAN, 1.0])],
        ..Default::default()
    });
    let e = check_result_finite(&r).unwrap_err().to_string();
    assert!(e.contains("'V'") && e.contains("1 NaN"), "{e}");
    r.solutions[0].fields[0].values[1] = 0.5;
    assert!(check_result_finite(&r).is_ok());
}
