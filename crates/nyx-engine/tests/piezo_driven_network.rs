// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end M2 round-trip checks connecting `nyx-sms` `Study::Driven` -> `nyx-analysis::NetworkAnalysis`
//! plus `nyx-engine` `study_frequency` and `chain_key_for_project` calibration keys.

use femcore::material::{Anisotropy, Stiffness, Tensor};
use femcore::props::Materials;
use femcore::testkit::block;
use nyx_analysis::{network_applies, network_with_params};
use nyx_engine::estimate::chain_key_for_project;
use nyx_engine::materials::study_frequency;
use nyx_engine::project::Project;
use nyx_physics::{Mesh, ProblemSpec};
use nyx_sms::problem::{Bc, BcKind, Problem, Study};
use nyx_sms::Sms;
use serde_json::json;

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
