// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Frozen reference output of the dfs kernels (ADR 10): the golden files in `tests/golden/` were taken through the pre-restructure `solve::solve`
//! (femcore tests/golden_frames) and moved here with the kernels. The kernels must reproduce them to round-off (tolerance 1e-9).
//! `NYX_WRITE_GOLDEN=<prefix> cargo test --release -p nyx-dfs --test golden` rewrites them; do that only for a kernel known to be right.
mod common;
use femcore::material::*;
use nyx_dfs::{Bc, BcKind, Problem, Study};
use nyx_physics::testkit::{cell_fields_reproduce_the_plot_fields, golden};
use nyx_physics::{Discretisation, SolveResult};

fn dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// Solves and checks that the cell fields of the fields stage reproduce the plot fields.
fn checked(mesh: &Discretisation, p: &Problem) -> SolveResult {
    let r = common::solve(mesh, p);
    cell_fields_reproduce_the_plot_fields(mesh, &r);
    r
}

#[test]
fn golden_electrostatics_anisotropic_plates() {
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    // a generic orientation, so that every entry of the tensor (xz, xy, zy) is non-zero
    spec.materials.aniso.insert(
        1,
        Anisotropy {
            permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }),
            euler_deg: [30.0, 40.0, 10.0],
            ..Default::default()
        },
    );
    spec.bcs = vec![Bc { labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { labels: vec![3], kind: BcKind::Ground }];
    golden(&dir(), "es_aniso_plates", &checked(&mesh, &spec), 1e-9);
}

#[test]
fn golden_electrostatics_series_with_charge() {
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.materials.set(2, "rho_v", 1e-9);
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Potential(1.0) },
        Bc { labels: vec![4], kind: BcKind::Ground },
        Bc { labels: vec![3], kind: BcKind::SurfaceCharge(1e-10) },
    ];
    golden(&dir(), "es_series_charge", &checked(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_anisotropic_stationary() {
    let mesh = common::rect(1.0, 0.4, 10, 4);
    let mut spec = Problem::new(Study::HeatStationary, 2);
    spec.materials.aniso.insert(
        1,
        Anisotropy {
            conductivity: Some(Tensor::Diagonal { x: 30.0, y: 5.0, z: 2.0 }),
            euler_deg: [10.0, 30.0, 50.0],
            ..Default::default()
        },
    );
    spec.materials.default.insert("Q".into(), 1e4);
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h: 20.0, text: 280.0 } },
        Bc { labels: vec![3], kind: BcKind::Flux(500.0) },
    ];
    golden(&dir(), "heat_aniso_stationary", &checked(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_transient_dense() {
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let mut spec = Problem::new(Study::HeatStationary, 2);
    spec.materials.default.insert("k".into(), 2.0);
    spec.materials.default.insert("rho".into(), 1.0);
    spec.materials.default.insert("Cp".into(), 1.0);
    spec.bcs = vec![Bc { labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    spec.study = Study::HeatTransient { t_end: 0.2, n_out: 4, t0: "sin(pi*x)".into() };
    golden(&dir(), "heat_transient", &checked(&mesh, &spec), 1e-9);
}
