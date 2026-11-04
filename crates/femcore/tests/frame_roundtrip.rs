//! What a wafer-frame problem on a 2D cross-section says, solved through `solve::solve`, against physics that has a direction: every
//! physics and every result field it makes is named in the wafer axes (an unnamed field is a failure; `solve` also asserts it in debug
//! builds), the fields that carry a direction point the right way; the component that points out of the plane is wafer y.
//! The mesh of `common::rect` has the columns (x, z); its boundaries are 1 x = 0, 2 z = 0 (bottom), 3 z = H (top), 4 x = L.
mod common;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;

fn names(r: &SolveResult) -> Vec<String> {
    let mut v: Vec<String> = r.solutions[0].fields.iter().map(|f| f.name.clone()).collect();
    v.sort();
    v
}

fn sorted(v: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = v.iter().map(|s| s.to_string()).collect();
    v.sort();
    v
}

fn mean(r: &SolveResult, field: &str) -> f64 {
    let v = &r.solutions[0].field(field).unwrap_or_else(|| panic!("no field {field} in {:?}", names(r))).values;
    v.iter().sum::<f64>() / v.len() as f64
}

fn peak(r: &SolveResult, sol: usize, field: &str) -> f64 {
    r.solutions[sol].field(field).unwrap_or_else(|| panic!("no field {field}")).values.iter().fold(0.0, |m, v| f64::max(m, v.abs()))
}

fn solve(mesh: &SimplexMesh, spec: &ProblemSpec) -> SolveResult {
    femcore::solve::solve(mesh, spec, &mut SolveState::default()).unwrap()
}

#[test]
fn a_plate_pressed_from_above_moves_down_in_z_and_is_stressed_in_z_and_y() {
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut s = ProblemSpec::new(Physics::Solid);
    s.materials.default.insert("E".into(), 70e9);
    s.materials.default.insert("nu".into(), 0.3);
    // fixed on the bottom, a downward traction on the top
    s.bcs = vec![
        Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { on: None, labels: vec![3], kind: BcKind::Load([0.0, 0.0, -1e6]) },
    ];
    let r = solve(&mesh, &s);
    assert_eq!(names(&r), sorted(&["sxx", "sxz", "syy", "szz", "ux", "uz", "von Mises", "|u|"]));
    assert!(mean(&r, "uz") < 0.0, "the plate is pushed down: {}", mean(&r, "uz"));
    // plane strain: sigma_zz = -1 MPa, and the out-of-plane stress is sigma_yy = nu (sigma_xx + sigma_zz), about -0.3 MPa (sigma_xx is not
    // zero near the clamped bottom)
    assert!((mean(&r, "szz") / -1e6 - 1.0).abs() < 0.02, "{}", mean(&r, "szz"));
    assert!(mean(&r, "syy") < -0.2e6 && mean(&r, "syy") > -0.45e6, "{}", mean(&r, "syy"));
}

#[test]
fn the_lowest_mode_of_a_pec_guide_has_its_field_along_the_short_side_z() {
    // 30 x 15 mm, air, 12 GHz: TE10 has E along the short side, which in the wafer frame is z
    let mesh = common::rect(0.03, 0.015, 24, 12);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.materials.default.insert("eps_r".into(), 1.0);
    s.study = Study::Modes { freq: 12e9, n_target: 0.9, evanescent: false, n_modes: 1 };
    let r = solve(&mesh, &s);
    assert_eq!(names(&r), sorted(&["Ex", "Ey", "Ez", "|E|", "|Ey|", "|Exz|"]));
    assert!(peak(&r, 0, "Ez") > 20.0 * peak(&r, 0, "Ex"), "Ez {} Ex {}", peak(&r, 0, "Ez"), peak(&r, 0, "Ex"));
}
