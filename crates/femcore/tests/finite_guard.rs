//! NaN or infinite input must be refused with a message that names the culprit, and a solve whose
//! result contains non-finite numbers must not be returned as a success.
use femcore::mesh::SimplexMesh;
use femcore::problem::*;
use femcore::solve::solve;

fn solid_spec(mesh: &SimplexMesh) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.order = 1;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    let labels = mesh.boundary_labels();
    s.bcs = vec![Bc { on: None, labels: vec![labels[0]], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![*labels.last().unwrap()], kind: BcKind::Load([0.0, -1e5, 0.0]) }];
    s
}

#[test]
fn a_valid_problem_still_solves() {
    let m = SimplexMesh::rectangle(1.0, 0.2, 10, 3);
    assert!(solve(&m, &solid_spec(&m), &mut SolveState::default()).is_ok());
}

#[test]
fn nan_and_infinite_materials_are_refused_by_name() {
    let m = SimplexMesh::rectangle(1.0, 0.2, 10, 3);
    for (key, bad) in [("E", f64::NAN), ("nu", f64::INFINITY), ("rho", f64::NEG_INFINITY)] {
        let mut s = solid_spec(&m);
        s.materials.default.insert(key.into(), bad);
        let e = solve(&m, &s, &mut SolveState::default()).unwrap_err();
        assert!(e.contains(&format!("'{key}'")), "{e}");
    }
    let mut s = solid_spec(&m);
    s.materials.set(1, "E", f64::NAN);
    let e = solve(&m, &s, &mut SolveState::default()).unwrap_err();
    assert!(e.contains("domain 1"), "{e}");
}

#[test]
fn nan_geometry_is_refused() {
    let m = SimplexMesh::rectangle(f64::NAN, 0.2, 10, 3);
    let e = solve(&m, &solid_spec(&m), &mut SolveState::default()).unwrap_err();
    assert!(e.contains("non-finite coordinate"), "{e}");
}

#[test]
fn nan_boundary_values_and_study_parameters_are_refused() {
    let m = SimplexMesh::rectangle(1.0, 0.2, 10, 3);
    let mut s = solid_spec(&m);
    s.bcs[1].kind = BcKind::Load([0.0, f64::NAN, 0.0]);
    let e = solve(&m, &s, &mut SolveState::default()).unwrap_err();
    assert!(e.contains("boundary condition 2"), "{e}");
    let mut s = solid_spec(&m);
    s.study = Study::Eigenfrequency { n_modes: 4, f_search: f64::NAN };
    let e = solve(&m, &s, &mut SolveState::default()).unwrap_err();
    assert!(e.contains("search frequency"), "{e}");
    let mut s = solid_spec(&m);
    s.thickness = f64::INFINITY;
    assert!(solve(&m, &s, &mut SolveState::default()).is_err());
}

#[test]
fn electrostatics_and_heat_refuse_nan_materials_too() {
    let m = SimplexMesh::rectangle(1.0, 1.0, 6, 6);
    let mut es = ProblemSpec::new(Physics::Electrostatics);
    es.materials.default.insert("eps_r".into(), f64::NAN);
    es.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    assert!(solve(&m, &es, &mut SolveState::default()).unwrap_err().contains("'eps_r'"));
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.materials.default.insert("k".into(), f64::NAN);
    heat.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(1.0) }];
    assert!(solve(&m, &heat, &mut SolveState::default()).unwrap_err().contains("'k'"));
}

#[test]
fn a_non_finite_result_is_an_error() {
    let mut r = SolveResult::default();
    r.solutions.push(Solution { label: "static".into(), fields: vec![Field { name: "ux".into(), unit: "m".into(), values: vec![0.0, f64::NAN, 1.0] }], ..Default::default() });
    let e = femcore::validate::check_result(&r).unwrap_err();
    assert!(e.contains("'ux'") && e.contains("1 NaN"), "{e}");
    r.solutions[0].fields[0].values[1] = 0.5;
    assert!(femcore::validate::check_result(&r).is_ok());
}
