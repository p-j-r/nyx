//! The input checks of solid mechanics moved to `nyx-sms/tests/finite_guard.rs`; this keeps the result check of the old result shape (`validate::check_result`,
//! used by the window path for the physics that have not moved yet).
mod common;
use femcore::problem::*;

#[test]
fn a_non_finite_result_is_an_error() {
    let mut r = SolveResult::default();
    r.solutions.push(Solution { label: "static".into(), fields: vec![Field { name: "ux".into(), unit: "m".into(), values: vec![0.0, f64::NAN, 1.0] }], ..Default::default() });
    let e = femcore::validate::check_result(&r).unwrap_err();
    assert!(e.contains("'ux'") && e.contains("1 NaN"), "{e}");
    r.solutions[0].fields[0].values[1] = 0.5;
    assert!(femcore::validate::check_result(&r).is_ok());
}
