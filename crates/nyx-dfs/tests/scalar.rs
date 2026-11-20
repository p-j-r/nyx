// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electrostatics and heat transfer against closed forms.
mod common;
use common::*;
use nyx_dfs::{Bc, BcKind, Problem, Study};

fn es_spec(eps: f64) -> Problem {
    let mut s = Problem::new(Study::Electrostatic, 2);
    s.materials.default.insert("eps_r".into(), eps);
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![4], kind: BcKind::Ground }];
    s
}

#[test]
fn es2d_parallel_plates_have_a_linear_potential() {
    // 1 x 0.5 rectangle, V = 1 on boundary 1 (x = 0), ground on 4 (x = 1), eps_r 4.2: V = 1 - x, which P2 reproduces exactly
    let mesh = common::rect(1.0, 0.5, 6, 3);
    let r = nyx_dfs::solve_electrostatics(&mesh, &es_spec(4.2)).unwrap();
    let v = &r.solutions[0].field("V").unwrap().values;
    let x = mesh.dof_coordinates(2);
    let exact: Vec<f64> = (0..v.len()).map(|i| 1.0 - x[2 * i]).collect();
    assert!(maxdiff(v, &exact) < 1e-12, "{}", maxdiff(v, &exact));
    // C = eps0 eps_r (0.5 / 1) per metre
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (nyx_dfs::EPS0 * 4.2) - 0.5).abs() < 1e-6);
}

#[test]
fn es3d_parallel_plates_have_a_linear_potential() {
    // V = 1 on the face x = 0 (boundary 1), ground on x = lx (boundary 4 of the block is z = lz, so use 6)
    let mesh = common::block(1.0, 0.4, 0.3, 4, 2, 2);
    let mut s = es_spec(4.2);
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![6], kind: BcKind::Ground }];
    let r = nyx_dfs::solve_electrostatics(&mesh, &s).unwrap();
    let v = &r.solutions[0].field("V").unwrap().values;
    let x = mesh.dof_coordinates(2);
    let exact: Vec<f64> = (0..v.len()).map(|i| 1.0 - x[3 * i]).collect();
    assert!(maxdiff(v, &exact) < 1e-9, "{}", maxdiff(v, &exact));
}

#[test]
fn es_series_capacitor_with_two_dielectrics() {
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut s = es_spec(1.0);
    s.materials.set(1, "eps_r", 2.0);
    s.materials.set(2, "eps_r", 8.0);
    let r = nyx_dfs::solve_electrostatics(&mesh, &s).unwrap();
    let exact = nyx_dfs::EPS0 / (0.5 / 2.0 + 0.5 / 8.0);
    assert!((r.scalar("capacitance").unwrap() / exact - 1.0).abs() < 1e-9);
}

fn heat_spec(k: f64) -> Problem {
    let mut s = Problem::new(Study::HeatStationary, 2);
    s.materials.default.insert("k".into(), k);
    s
}

#[test]
fn heat_convection_gives_the_1d_profile() {
    // 1 x 0.2 plate, T = 100 at x = 0, convection (h, text) at x = 1, the long sides insulated:
    // T = T0 - (T0 - text) h x / (k + h L), linear, so P2 reproduces it exactly
    let mesh = common::rect(1.0, 0.2, 8, 2);
    let (k, h, text, t0) = (2.0, 5.0, 20.0, 100.0);
    let mut s = heat_spec(k);
    s.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(t0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h, text } },
    ];
    let r = nyx_dfs::solve_heat(&mesh, &s).unwrap();
    let t = &r.solutions[0].field("T").unwrap().values;
    let x = mesh.dof_coordinates(2);
    let exact: Vec<f64> = (0..t.len()).map(|i| t0 - (t0 - text) * h * x[2 * i] / (k + h)).collect();
    assert!(maxdiff(t, &exact) < 1e-9, "{}", maxdiff(t, &exact));
}

#[test]
fn heat_source_gives_the_parabolic_profile() {
    // -k T'' = Q, T = 0 at x = 0, insulated at x = L: T = Q (L x - x^2 / 2) / k, quadratic, so P2 is exact
    let mesh = common::rect(1.0, 0.2, 8, 2);
    let (k, q) = (2.0, 3.0);
    let mut s = heat_spec(k);
    s.materials.default.insert("Q".into(), q);
    s.bcs = vec![Bc { labels: vec![1], kind: BcKind::Temperature(0.0) }];
    let r = nyx_dfs::solve_heat(&mesh, &s).unwrap();
    let t = &r.solutions[0].field("T").unwrap().values;
    let x = mesh.dof_coordinates(2);
    let exact: Vec<f64> = (0..t.len()).map(|i| q * (x[2 * i] - 0.5 * x[2 * i] * x[2 * i]) / k).collect();
    assert!(maxdiff(t, &exact) < 1e-9, "{}", maxdiff(t, &exact));
}

#[test]
fn transient_heat_follows_the_decaying_sine() {
    // 1 x 0.2, k 2, rho cp 1, T0 = sin(pi x), T = 0 on boundaries 1 and 4, 21 output times 0..0.2:
    // T = sin(pi x) exp(-k pi^2 t); the dense semi-discrete solution is exact in time, so only the mesh remains
    let mesh = common::rect(1.0, 0.2, 20, 2);
    let mut s = heat_spec(2.0);
    s.materials.default.insert("rho".into(), 1.0);
    s.materials.default.insert("Cp".into(), 1.0);
    s.bcs = vec![Bc { labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    s.study = Study::HeatTransient { t_end: 0.2, n_out: 20, t0: "sin(pi*x)".into() };
    let r = nyx_dfs::solve_heat(&mesh, &s).unwrap();
    assert_eq!(r.solutions.len(), 21);
    assert!(!r.log[0].starts_with("BDF2"), "{:?}", r.log);
    let x = mesh.dof_coordinates(2);
    for (j, sol) in r.solutions.iter().enumerate() {
        let t = &sol.field("T").unwrap().values;
        let decay = (-2.0 * std::f64::consts::PI.powi(2) * 0.2 * j as f64 / 20.0).exp();
        let exact: Vec<f64> = (0..t.len()).map(|i| (std::f64::consts::PI * x[2 * i]).sin() * decay).collect();
        let e = maxdiff(t, &exact);
        assert!(e < 1e-4, "time index {j}: {e}");
    }
}

#[test]
fn transient_bdf2_agrees_with_decay() {
    // big enough to leave the dense path: analytic decay of the first mode
    let mesh = common::rect(1.0, 0.2, 60, 14);
    let mut s = heat_spec(2.0);
    s.materials.default.insert("rho".into(), 1.0);
    s.materials.default.insert("Cp".into(), 1.0);
    s.bcs = vec![Bc { labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    s.study = Study::HeatTransient { t_end: 0.2, n_out: 10, t0: "sin(pi*x)".into() };
    let r = nyx_dfs::solve_heat(&mesh, &s).unwrap();
    assert!(r.log[0].starts_with("BDF2"), "{:?}", r.log);
    let x = mesh.dof_coordinates(2);
    let t = &r.solutions[10].field("T").unwrap().values;
    let decay = (-2.0 * std::f64::consts::PI.powi(2) * 0.2).exp();
    let exact: Vec<f64> = (0..t.len()).map(|i| (std::f64::consts::PI * x[2 * i]).sin() * decay).collect();
    assert!(maxdiff(t, &exact) < 2e-3, "{}", maxdiff(t, &exact));
}

#[test]
fn transient_convection_cools_to_ambient() {
    let mesh = common::rect(1.0, 1.0, 6, 6);
    let mut s = heat_spec(100.0);
    s.materials.default.insert("rho".into(), 1000.0);
    s.materials.default.insert("Cp".into(), 500.0);
    s.bcs = vec![Bc { labels: vec![1, 2, 3, 4], kind: BcKind::Convection { h: 2.0, text: 300.0 } }];
    s.study = Study::HeatTransient { t_end: 1e5, n_out: 4, t0: "400".into() };
    let r = nyx_dfs::solve_heat(&mesh, &s).unwrap();
    let maxes: Vec<f64> =
        r.solutions.iter().map(|s| s.field("T").unwrap().values.iter().cloned().fold(f64::MIN, f64::max)).collect();
    assert!((maxes[0] - 400.0).abs() < 1e-9);
    assert!(maxes.windows(2).all(|w| w[1] < w[0]), "{maxes:?}");
    // small Biot number (0.02): lumped estimate, tau = rho cp A / (h P) = 1000*500*1/(2*4) = 62500 s
    let lumped = 300.0 + 100.0 * (-1e5f64 / 62500.0).exp();
    assert!((maxes[4] - lumped).abs() < 1.0, "{} vs {lumped}", maxes[4]);
}

#[test]
fn errors_are_reported() {
    let mesh = common::rect(1.0, 1.0, 4, 4);
    assert!(nyx_dfs::solve_electrostatics(&mesh, &Problem::new(Study::Electrostatic, 2)).is_err());
    let mut s = es_spec(1.0);
    s.bcs[0].labels = vec![9];
    assert!(nyx_dfs::solve_electrostatics(&mesh, &s).unwrap_err().contains("do not exist"));
    assert!(nyx_dfs::solve_heat(&mesh, &heat_spec(1.0)).is_err());
    let mut h = heat_spec(1.0);
    h.bcs = vec![Bc { labels: vec![1], kind: BcKind::Temperature(1.0) }];
    h.study = Study::HeatTransient { t_end: 1.0, n_out: 2, t0: "foo(x)".into() };
    assert!(nyx_dfs::solve_heat(&mesh, &h).unwrap_err().contains("initial temperature"));
}
