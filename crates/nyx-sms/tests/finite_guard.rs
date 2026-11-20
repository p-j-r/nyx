// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! NaN or infinite input must be refused with a message that names the culprit, and a solve whose
//! result contains non-finite numbers must not be returned as a success.
mod common;
use common::try_solve;
use femcore::mesh::SimplexMesh;
use nyx_physics::{check_result_finite, Field, Solution, SolveResult};
use nyx_sms::*;

fn solid_spec(mesh: &SimplexMesh) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.order = 1;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    let labels = mesh.boundary_labels();
    s.bcs = vec![
        Bc { labels: vec![labels[0]], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![*labels.last().unwrap()], kind: BcKind::Load([0.0, -1e5, 0.0]) },
    ];
    s
}

#[test]
fn a_valid_problem_still_solves() {
    let m = common::rect(1.0, 0.2, 10, 3);
    assert!(try_solve(&m, &solid_spec(&m)).is_ok());
}

#[test]
fn nan_and_infinite_materials_are_refused_by_name() {
    let m = common::rect(1.0, 0.2, 10, 3);
    for (key, bad) in [("E", f64::NAN), ("nu", f64::INFINITY), ("rho", f64::NEG_INFINITY)] {
        let mut s = solid_spec(&m);
        s.materials.default.insert(key.into(), bad);
        let e = try_solve(&m, &s).unwrap_err();
        assert!(e.contains(&format!("'{key}'")), "{e}");
    }
    let mut s = solid_spec(&m);
    s.materials.set(1, "E", f64::NAN);
    let e = try_solve(&m, &s).unwrap_err();
    assert!(e.contains("domain 1"), "{e}");
}

#[test]
fn nan_geometry_is_refused() {
    let mut m = common::rect(1.0, 0.2, 10, 3);
    m.points[3] = f64::NAN;
    let e = try_solve(&m, &solid_spec(&m)).unwrap_err();
    assert!(e.contains("non-finite coordinate"), "{e}");
}

#[test]
fn nan_boundary_values_and_study_parameters_are_refused() {
    let m = common::rect(1.0, 0.2, 10, 3);
    let mut s = solid_spec(&m);
    s.bcs[1].kind = BcKind::Load([0.0, f64::NAN, 0.0]);
    let e = try_solve(&m, &s).unwrap_err();
    assert!(e.contains("boundary condition 2"), "{e}");
    let mut s = solid_spec(&m);
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: f64::NAN };
    let e = try_solve(&m, &s).unwrap_err();
    assert!(e.contains("search frequency"), "{e}");
    let mut s = solid_spec(&m);
    s.thickness = f64::INFINITY;
    assert!(try_solve(&m, &s).is_err());
}

#[test]
fn a_non_finite_result_is_an_error() {
    let mut r = SolveResult::default();
    r.solutions.push(Solution {
        label: "static".into(),
        fields: vec![Field::scalar("ux", "m", vec![0.0, f64::NAN, 1.0])],
        ..Default::default()
    });
    let e = check_result_finite(&r).unwrap_err().to_string();
    assert!(e.contains("'ux'") && e.contains("1 NaN"), "{e}");
    r.solutions[0].fields[0].values[1] = 0.5;
    assert!(check_result_finite(&r).is_ok());
}
