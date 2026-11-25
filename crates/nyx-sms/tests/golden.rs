// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Frozen reference output of the sms kernels (ADR 10): the golden files in `tests/golden/` were taken through the pre-restructure
//! `solve::solve` (femcore tests/golden_frames) and moved here with the kernels, unchanged. The kernels must reproduce them (tolerances as
//! before, 1e-8 for the static cases and 1e-7 for the complex-frequency ones).
//! `NYX_WRITE_GOLDEN=solid cargo test --release -p nyx-sms --test golden` rewrites them; do that only for a kernel known to be right.
mod common;
use femcore::material::*;
use femcore::props::*;
use nyx_physics::testkit::{cell_fields_reproduce_the_plot_fields, golden};
use nyx_physics::{Discretisation, SolveResult};
use nyx_sms::*;

fn dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// Solves and checks that the cell fields of the fields stage reproduce the plot fields.
fn solve(mesh: &Discretisation, p: &Problem) -> SolveResult {
    let r = common::solve(mesh, p);
    cell_fields_reproduce_the_plot_fields(mesh, &r);
    r
}

fn sapphire_euler() -> [f64; 3] {
    [30.0, 40.0, 10.0]
}

fn sapphire() -> Anisotropy {
    let st = Stiffness::Trigonal {
        c11: 497.35e9,
        c12: 163.97e9,
        c13: 112.20e9,
        c14: -23.58e9,
        c33: 499.11e9,
        c44: 147.39e9,
    };
    Anisotropy {
        stiffness: Some(st),
        expansion: Some(Tensor::Diagonal { x: 5.0e-6, y: 5.0e-6, z: 6.0e-6 }),
        euler_deg: sapphire_euler(),
        ..Default::default()
    }
}

fn clamped_left(spec: &mut Problem) {
    spec.bcs.push(Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
}

#[test]
fn golden_solid_static_tilted_sapphire_plane_strain() {
    // a generic crystal orientation turns every entry of the stiffness; loads along x and z, a sliding wall at the bottom
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Strain;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    clamped_left(&mut s);
    s.bcs.push(Bc { labels: vec![4], kind: BcKind::Load([2.0e6, 0.0, -1.5e6]) });
    s.bcs.push(Bc { labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0, 0.0, 1.0e-7] } });
    golden(&dir(), "solid_static_sapphire_strain", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_static_plane_stress_with_thermal_expansion() {
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut s = Problem::new(Study::Static, 2);
    s.plane = Plane::Stress;
    s.materials.aniso.insert(1, sapphire());
    s.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
    clamped_left(&mut s);
    golden(&dir(), "solid_static_sapphire_stress_thermal", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_eigenfrequencies_tilted_sapphire() {
    let mesh = common::rect(1.0e-3, 0.2e-3, 10, 3);
    let mut s = Problem::new(Study::Eigenfrequency { n_modes: 5, f_search: 0.0 }, 2);
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    clamped_left(&mut s);
    golden(&dir(), "solid_eigen_sapphire", &solve(&mesh, &s), 1e-8);
}

#[test]
fn golden_solid_eigenfrequencies_floquet_cell() {
    // a phononic cell: the sides x = 0 and x = lx are periodic with a phase
    let mesh = common::rect(1.0e-3, 0.3e-3, 8, 3);
    let mut s = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 1.0e6 }, 2);
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 1.2 }];
    golden(&dir(), "solid_eigen_floquet", &solve(&mesh, &s), 1e-7);
}

#[test]
fn golden_solid_eigenfrequencies_with_a_pml() {
    // a clamped isotropic bar whose right end is an absorbing layer along x
    let mesh = common::rect(1.0, 0.1, 30, 3);
    let mut s = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    s.plane = Plane::Strain;
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        s.materials.default.insert(k.into(), v);
    }
    clamped_left(&mut s);
    s.materials.pml.insert(
        1,
        Pml {
            x: PmlDir::Plus,
            strength: 0.1,
            order: 0.0,
            clamp_outer: false,
            stabilization: PmlStabilization::Off,
            ..Pml::default()
        },
    );
    golden(&dir(), "solid_eigen_pml_x", &solve(&mesh, &s), 1e-7);
}

#[test]
fn golden_solid_eigenfrequencies_with_a_pml_along_z() {
    // a plate on a substrate: the layer is the lower part of a tall section and absorbs along the surface normal z
    let mesh = common::rect(0.2, 1.0, 3, 30);
    let mut s = Problem::new(Study::Eigenfrequency { n_modes: 4, f_search: 0.0 }, 2);
    s.plane = Plane::Strain;
    s.materials.aniso.insert(1, sapphire());
    s.materials.default.insert("rho".into(), 3980.0);
    s.bcs.push(Bc { labels: vec![3], kind: BcKind::Fixed { fixed: [true, true, true], values: [0.0; 3] } });
    s.materials.pml.insert(
        1,
        Pml {
            z: PmlDir::Minus,
            strength: 0.1,
            order: 0.0,
            clamp_outer: false,
            stabilization: PmlStabilization::Off,
            ..Pml::default()
        },
    );
    golden(&dir(), "solid_eigen_pml_z", &solve(&mesh, &s), 1e-7);
}
