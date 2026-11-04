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
