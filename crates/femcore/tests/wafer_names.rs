//! Every field a ported 2D kernel returns carries a wafer name straight from the kernel (no renaming step): the in-plane components of a
//! cross-section are `x` and `z`, the out-of-plane one is `y`.
mod common;
use femcore::material::*;
use femcore::problem::*;
use std::collections::BTreeSet;

fn names(r: &SolveResult) -> BTreeSet<String> {
    let mut all = BTreeSet::new();
    for s in &r.solutions {
        let these: BTreeSet<String> = s.fields.iter().map(|f| f.name.clone()).collect();
        if all.is_empty() {
            all = these;
        } else {
            assert_eq!(all, these, "every solution of a study has the same fields");
        }
    }
    all
}

fn set(v: &[&str]) -> BTreeSet<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn scalar_kernels_name_their_components_x_and_z() {
    let mesh = common::rect(1.0, 0.5, 6, 3);
    let mut es = ProblemSpec::new(Physics::Electrostatics);
    es.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    assert_eq!(names(&femcore::solve::solve(&mesh, &es, &mut SolveState::default()).unwrap()), set(&["V", "|E|", "Ex", "Ez"]));
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.materials.default.insert("k".into(), 1.0);
    heat.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Temperature(0.0) }];
    assert_eq!(names(&femcore::solve::solve(&mesh, &heat, &mut SolveState::default()).unwrap()), set(&["T", "|q|", "qx", "qz"]));
}

#[test]
fn optical_mode_solvers_name_their_components_in_wafer_axes() {
    let mesh = common::rect(0.03, 0.015, 12, 6);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.aniso.insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }), ..Default::default() });
    s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 2 };
    // the real solver
    let real = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    assert_eq!(names(&real), set(&["|E|", "|Exz|", "|Ey|", "Ex", "Ey", "Ez"]));
    // the complex solver (a tilted crystal)
    s.materials.aniso.get_mut(&1).unwrap().euler_deg = [0.0, 30.0, 0.0];
    let complex = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    assert_eq!(names(&complex), set(&["|E|", "|Exz|", "|Ey|", "Re Ex", "Re Ey", "Re Ez"]));
}

#[test]
fn an_initial_temperature_on_a_cross_section_is_written_in_x_and_z() {
    // 1 x 0.2 section (x lateral, z vertical), T = 0 on the bottom and the top: T0 = sin(pi z / 0.2) decays as exp(-k (pi / 0.2)^2 t)
    let mesh = common::rect(1.0, 0.2, 20, 4);
    let mut s = ProblemSpec::new(Physics::Heat);
    s.materials.default.insert("k".into(), 2.0);
    s.materials.default.insert("rho".into(), 1.0);
    s.materials.default.insert("Cp".into(), 1.0);
    s.bcs = vec![Bc { on: None, labels: vec![2, 3], kind: BcKind::Temperature(0.0) }];
    s.study = Study::Transient { t_end: 0.001, n_out: 2, t0: "sin(pi*z/0.2)".into() };
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let t = &r.solutions[2].field("T").unwrap().values;
    let xz = mesh.dof_coordinates(2);
    let decay = (-2.0 * (std::f64::consts::PI / 0.2).powi(2) * 0.001).exp();
    let exact: Vec<f64> = (0..t.len()).map(|i| (std::f64::consts::PI * xz[2 * i + 1] / 0.2).sin() * decay).collect();
    assert!(common::maxdiff(t, &exact) < 5e-3, "{}", common::maxdiff(t, &exact));
}
