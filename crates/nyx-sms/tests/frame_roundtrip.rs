// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! What a wafer-frame problem on a 2D cross-section says, solved through `nyx_sms::run`, against physics that has a direction: every result
//! field is named in the wafer axes (flat names `ux`, `uz`, `|u|`), the fields that carry a direction point the right way; the component that
//! points out of the plane is wafer y. The mesh of `common::rect` has the columns (x, z); its boundaries are 1 x = 0, 2 z = 0 (bottom),
//! 3 z = H (top), 4 x = L. (The electromagnetic part stays in femcore.)
mod common;
use nyx_physics::SolveResult;
use nyx_sms::*;

fn names(r: &SolveResult) -> Vec<String> {
    let mut v: Vec<String> = r.solutions[0].flat_fields().into_iter().map(|(n, _, _)| n).collect();
    v.sort();
    v
}

fn sorted(v: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = v.iter().map(|s| s.to_string()).collect();
    v.sort();
    v
}

fn mean(r: &SolveResult, field: &str) -> f64 {
    let v = r.solutions[0].flat(field).unwrap_or_else(|| panic!("no field {field} in {:?}", names(r)));
    v.iter().sum::<f64>() / v.len() as f64
}

#[test]
fn a_plate_pressed_from_above_moves_down_in_z_and_is_stressed_in_z_and_y() {
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut s = Problem::new(Study::Static, 2);
    s.materials.default.insert("E".into(), 70e9);
    s.materials.default.insert("nu".into(), 0.3);
    // fixed on the bottom, a downward traction on the top
    s.bcs = vec![
        Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } },
        Bc { labels: vec![3], kind: BcKind::Load([0.0, 0.0, -1e6]) },
    ];
    let r = common::solve(&mesh, &s);
    assert_eq!(names(&r), sorted(&["sxx", "sxz", "syy", "szz", "ux", "uz", "von Mises", "|u|"]));
    assert!(mean(&r, "uz") < 0.0, "the plate is pushed down: {}", mean(&r, "uz"));
    // plane strain: sigma_zz = -1 MPa, and the out-of-plane stress is sigma_yy = nu (sigma_xx + sigma_zz), about -0.3 MPa (sigma_xx is not
    // zero near the clamped bottom)
    assert!((mean(&r, "szz") / -1e6 - 1.0).abs() < 0.02, "{}", mean(&r, "szz"));
    assert!(mean(&r, "syy") < -0.2e6 && mean(&r, "syy") > -0.45e6, "{}", mean(&r, "syy"));
}
