// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end M2 round-trip checks connecting `nyx-sms` `Study::Driven` -> `nyx-analysis::NetworkAnalysis`
//! plus `nyx-engine` `study_frequency` and `chain_key_for_project` calibration keys.

use femcore::material::{coupled::EPSILON_0, Anisotropy, Stiffness, Tensor};
use femcore::props::Materials;
use femcore::testkit::block;
use nyx_analysis::{network_applies, network_with_params, piezo_resonator_applies, piezo_resonator_with_params};
use nyx_engine::analyse::analyses;
use nyx_engine::estimate::chain_key_for_project;
use nyx_engine::materials::study_frequency;
use nyx_engine::project::Project;
use nyx_physics::{Mesh, ProblemSpec};
use nyx_sms::problem::{Bc, BcKind, Problem, Study};
use nyx_sms::Sms;
use serde_json::json;
use std::f64::consts::PI;

#[test]
fn sms_driven_piezo_round_trips_through_network_analysis_and_engine_keys() {
    let simplex = block(1.0e-6, 1.0e-6, 1.0e-6, 2, 2, 4);
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();

    let mut mats = Materials::default();
    mats.default.insert("rho".into(), 3260.0);
    mats.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: 389e9, c44: 125e9 }),
            permittivity: Some(Tensor::Diagonal { x: 9.21, y: 9.21, z: 10.12 }),
            piezo_e: Some([
                [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
                [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
                [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0],
            ]),
            ..Default::default()
        },
    );

    let mut prob = Problem::new(
        Study::Driven {
            f_min: 2.3e9,
            f_max: 2.9e9,
            n_points: 15,
            freqs_hz: None,
            prom_tol: Some(1e-5),
            max_hdm_solves: 10,
            strict_prom_tol: false,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        1,
    );
    prob.materials = mats;
    prob.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } },
    ];

    let res = nyx_sms::run(&mesh, &prob, &femcore::progress::Progress::new()).unwrap();
    assert!(network_applies(&res), "NetworkAnalysis must apply to sms driven SolveResult");
    let net = network_with_params(&res, &json!({ "touchstone_format": "RI" })).unwrap();
    assert_eq!(net.ports, vec!["top"]);
    assert_eq!(net.freqs_hz.len(), 15);
    assert!(net.touchstone_snp.as_ref().unwrap().contains("# Hz S RI R 50"));

    // Engine study_frequency + chain_key_for_project checks
    let mut project = Project::new("PiezoDriven", 3);
    let mut driven_spec = ProblemSpec::new("sms", 2, "driven");
    driven_spec.study.params = json!({ "f_min": 2.0e9, "f_max": 3.0e9, "solver": "blr" });
    assert_eq!(study_frequency(&Sms, &driven_spec), Some(2.5e9));

    project.physics.push(driven_spec);
    assert_eq!(chain_key_for_project(&project, &["sms".into()]), "sms:driven:blr:p2");

    // Static piezo with electrical BCs uses sms:piezo:p2; pure elastic keeps "sms"
    let mut piezo_static = ProblemSpec::new("sms", 2, "static");
    piezo_static.bcs.push(nyx_physics::BcSpec { on: "bot".into(), kind: "ground".into(), params: json!({}) });
    project.physics[0] = piezo_static;
    assert_eq!(chain_key_for_project(&project, &["sms".into()]), "sms:piezo:p2");

    project.physics[0] = ProblemSpec::new("sms", 1, "static");
    assert_eq!(chain_key_for_project(&project, &["sms".into()]), "sms");
}

#[test]
fn sms_driven_baw_recovers_mbvd_and_analytical_kt2_via_piezo_resonator() {
    let (lx, ly, d) = (0.4e-6, 0.4e-6, 2.0e-6);
    let area = lx * ly;
    let simplex = block(lx, ly, d, 1, 1, 16);
    let mesh = Mesh::from_solver_mesh(&simplex).unwrap();

    let (c33_e, e33, eps33_s, rho) = (389e9_f64, 1.55_f64, 10.12_f64, 3260.0_f64);
    let mut visc = [[0.0; 6]; 6];
    for row in &mut visc {
        for v in row.iter_mut() {
            *v = 0.0;
        }
    }
    for (i, row) in visc.iter_mut().enumerate() {
        row[i] = 2.0e-4;
    }
    let mut mats = Materials::default();
    mats.default.insert("rho".into(), rho);
    mats.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Hexagonal { c11: 410e9, c12: 149e9, c13: 99e9, c33: c33_e, c44: 125e9 }),
            permittivity: Some(Tensor::Diagonal { x: 9.21, y: 9.21, z: eps33_s }),
            permittivity_im: Some(Tensor::Diagonal { x: 9.21 * 1e-3, y: 9.21 * 1e-3, z: eps33_s * 1e-3 }),
            piezo_e: Some([
                [0.0, 0.0, 0.0, 0.0, -0.48, 0.0],
                [0.0, 0.0, 0.0, -0.48, 0.0, 0.0],
                [-0.58, -0.58, e33, 0.0, 0.0, 0.0],
            ]),
            viscosity: Some(visc),
            ..Default::default()
        },
    );

    let mut prob = Problem::new(
        Study::Driven {
            f_min: 2.65e9,
            f_max: 2.90e9,
            n_points: 121,
            freqs_hz: None,
            prom_tol: Some(1e-7),
            max_hdm_solves: 14,
            strict_prom_tol: true,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-8,
        },
        2,
    );
    prob.materials = mats;
    prob.bcs = vec![
        Bc { labels: vec![1, 6], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { labels: vec![2, 5], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Ground },
        Bc { labels: vec![4], kind: BcKind::Terminal { name: "top".into(), z0_ohm: 50.0 } },
    ];

    let res = nyx_sms::run(&mesh, &prob, &femcore::progress::Progress::new()).unwrap();
    assert!(piezo_resonator_applies(&res));
    assert!(analyses().iter().any(|a| a.id() == "piezo_resonator"));

    let rep = piezo_resonator_with_params(&res, &json!({ "fix_r_s_ohm": 0.0 })).unwrap();
    let c33_d = c33_e + e33 * e33 / (eps33_s * EPSILON_0);
    let fp_exact = (c33_d / rho).sqrt() / (2.0 * d);
    let kt2_exact = (e33 * e33) / (c33_d * eps33_s * EPSILON_0);
    let c0_exact = eps33_s * EPSILON_0 * area / d;

    assert!((rep.f_p_hz / fp_exact - 1.0).abs() < 2e-3, "f_p mismatch: {} vs {fp_exact}", rep.f_p_hz);
    assert!((rep.k_t2_ieee / kt2_exact - 1.0).abs() < 2e-2, "k_t^2 mismatch: {} vs {kt2_exact}", rep.k_t2_ieee);
    assert!((rep.c0_f / c0_exact - 1.0).abs() < 2e-2, "C0 mismatch: {} vs {c0_exact}", rep.c0_f);
    assert!(rep.fit_residual < 5e-3, "fit_residual {}", rep.fit_residual);
    let _ = PI;
}
