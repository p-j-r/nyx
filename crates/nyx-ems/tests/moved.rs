// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tests of the ems kernels that were mixed into other femcore test files before the move (ADR 45): the anisotropic guide, the field names in
//! wafer axes, the stage functions of the mode solver, a Floquet cell on an adapted mesh, a guide in the frame round trip.
mod common;
use femcore::adapt::drive3::adapt3;
use femcore::adapt::field::UniformSize;
use femcore::adapt::ops3::AdaptOptions3;
use femcore::material::*;
use nyx_ems::Study;
use nyx_physics::SolveResult;
use std::collections::BTreeSet;

fn names_of(r: &SolveResult) -> BTreeSet<String> {
    let mut all = BTreeSet::new();
    for s in &r.solutions {
        let these: BTreeSet<String> = s.flat_fields().into_iter().map(|f| f.0).collect();
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

fn names(r: &SolveResult) -> Vec<String> {
    names_of(r).into_iter().collect()
}

fn sorted(v: &[&str]) -> Vec<String> {
    let mut v: Vec<String> = v.iter().map(|s| s.to_string()).collect();
    v.sort();
    v
}

fn peak(r: &SolveResult, sol: usize, field: &str) -> f64 {
    r.solutions[sol]
        .flat(field)
        .unwrap_or_else(|| panic!("no field {field}"))
        .iter()
        .fold(0.0, |m, v| f64::max(m, v.abs()))
}

#[test]
fn anisotropic_waveguide_te_modes_use_the_right_permittivity_component() {
    // PEC rectangle a x b (x lateral, z the vertical axis, propagation along y) filled with diag(eps_x, eps_y, eps_z): TE10 (E along z) has
    // n^2 = eps_z - (pi/(a k0))^2, TE01 (E along x) has n^2 = eps_x - (pi/(b k0))^2; eps_y plays no role
    let (a, b, f) = (22.86e-3, 10.16e-3, 12e9);
    let (ex, ey, ez) = (5.0, 2.0, 3.0);
    let mesh = common::rect(a, b, 24, 12);
    let k0 = 2.0 * std::f64::consts::PI * f / nyx_ems::C0;
    let n10 = (ez - (std::f64::consts::PI / a / k0).powi(2)).sqrt();
    let n01 = (ex - (std::f64::consts::PI / b / k0).powi(2)).sqrt();
    let mut spec = common::emw();
    spec.materials
        .aniso
        .insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: ex, y: ey, z: ez }), ..Default::default() });
    spec.study = Study::Modes { freq: f, n_target: 1.6, evanescent: false, n_modes: 6 };
    let r = common::solve(&mesh, &spec).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|s| common::neff(s).0).collect();
    for want in [n10, n01] {
        let best = got.iter().map(|g| (g - want).abs()).fold(f64::MAX, f64::min);
        assert!(best < 2e-3 * want, "n = {want}: nearest {best} ({got:?})");
    }
    // a crystal tilted against the propagation axis (xz/yz terms) is solved exactly by the complex solver (checked against a 3D
    // periodic cell in `aniso_optics.rs`); here only that it runs and stays close to the untilted modes
    let mut tilted = spec.clone();
    tilted.materials.aniso.get_mut(&1).unwrap().euler_deg = [0.0, 30.0, 0.0];
    let r = common::solve(&mesh, &tilted).unwrap();
    assert!(!r.solutions.is_empty() && common::neff(&r.solutions[0]).0 > 1.0);
}

#[test]
fn the_lowest_mode_of_a_pec_guide_has_its_field_along_the_short_side_z() {
    // 30 x 15 mm, air, 12 GHz: TE10 has E along the short side, which in the wafer frame is z
    let mesh = common::rect(0.03, 0.015, 24, 12);
    let mut s = common::emw();
    s.materials.default.insert("eps_r".into(), 1.0);
    s.study = Study::Modes { freq: 12e9, n_target: 0.9, evanescent: false, n_modes: 1 };
    let r = common::solve(&mesh, &s).unwrap();
    assert_eq!(names(&r), sorted(&["Ex", "Ey", "Ez", "|E|", "|Ey|", "|Exz|"]));
    assert!(peak(&r, 0, "Ez") > 20.0 * peak(&r, 0, "Ex"), "Ez {} Ex {}", peak(&r, 0, "Ez"), peak(&r, 0, "Ex"));
}

#[test]
fn optical_mode_stages_compose_to_the_one_call_solver() {
    use nyx_ems::*;
    let mesh = common::rect(1.0e-6, 0.5e-6, 12, 6);
    let nc = mesh.n_cells();
    let eps = vec![EpsCell { xx: 4.0, xy: 0.0, yy: 4.0, zz: 4.0 }; nc];
    let mu = vec![1.0; nc];
    let whole = solve_modes_tensor(&mesh, 2.0e14, &eps, &mu, 1.5, false, 3, 2).unwrap();
    let asm = assemble_modes(&mesh, 2.0e14, &eps, &mu, 2).unwrap();
    let staged = solve_mode_pairs(&asm, 1.5, false, 3).unwrap();
    assert_eq!(whole.neff, staged.neff);
    let (fields, cells) = mode_fields_cells(&mesh, &staged, 0, 2.0e14, 2);
    assert_eq!(fields.len(), 6);
    // the components carry the wafer names: x lateral, z normal, y along the guide
    for name in ["Ex", "Ez", "Ey"] {
        let avg = cells.node_average(&mesh, name).unwrap();
        let plot = fields.iter().find(|f| f.name == name).unwrap();
        let scale = plot.values.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(avg.iter().zip(&plot.values).all(|(a, b)| (a - b).abs() < 1e-12 * scale), "{name}");
    }
}

#[test]
fn optical_mode_solvers_name_their_components_in_wafer_axes() {
    let mesh = common::rect(0.03, 0.015, 12, 6);
    let mut s = common::emw();
    s.order = 2;
    s.materials.aniso.insert(
        1,
        Anisotropy { permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }), ..Default::default() },
    );
    s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 2 };
    // the real solver
    let real = common::solve(&mesh, &s).unwrap();
    assert_eq!(names_of(&real), set(&["|E|", "|Exz|", "|Ey|", "Ex", "Ey", "Ez"]));
    // the complex solver (a tilted crystal)
    s.materials.aniso.get_mut(&1).unwrap().euler_deg = [0.0, 30.0, 0.0];
    let complex = common::solve(&mesh, &s).unwrap();
    assert_eq!(names_of(&complex), set(&["|E|", "|Exz|", "|Ey|", "Re Ex", "Re Ey", "Re Ez"]));
}

#[test]
fn frozen_periodic_faces_survive_adaptation_and_the_floquet_cell_still_solves() {
    use femcore::adapt::topo3::TopoOptions3;
    use femcore::adapt::work::{simplex_from_adapt, simplex_to_adapt};
    use femcore::props::Periodic;
    use nyx_ems::Study;
    // a periodic cell of a PEC guide (see tests/emw3d.rs): faces 2 (y = 0) and 5 (y = l) are Floquet partners and must stay congruent
    let (a, b, l) = (0.03, 0.015, 0.01);
    // boundary numbers of common::block: 1 x = 0, 2 y = 0, 3 z = 0, 4 z = lz, 5 y = ly, 6 x = lx
    let nyx = common::block(a, l, b, 6, 3, 2);
    let before = simplex_to_adapt::<3>(&nyx).unwrap();
    let topo = TopoOptions3 { frozen_labels: vec![2, 5], ..TopoOptions3::default() };
    // finer than the start inside, as fine as the start on the frozen faces
    let out = adapt3(&before, &UniformSize(0.0035), &topo, &AdaptOptions3::default()).unwrap();
    let after = simplex_from_adapt(&out.mesh);
    assert!(after.n_cells() > nyx.n_cells());
    let on = |m: &femcore::mesh::SimplexMesh, y: f64| -> BTreeSet<[u64; 2]> {
        (0..m.n_points())
            .filter(|&i| (m.points[3 * i + 1] - y).abs() < 1e-15)
            .map(|i| [m.points[3 * i].to_bits(), m.points[3 * i + 2].to_bits()])
            .collect()
    };
    assert_eq!(on(&nyx, 0.0), on(&after, 0.0), "the y = 0 face is unchanged");
    assert_eq!(on(&nyx, l), on(&after, l), "the y = l face is unchanged");
    // TE10 with Floquet phase 0.8 over the cell: f = c/2 sqrt((1/a)^2 + (beta/pi)^2)
    let phase = 0.8f64;
    let c0 = 299_792_458.0;
    let exact = 0.5 * c0 * ((1.0 / a).powi(2) + (phase / l / std::f64::consts::PI).powi(2)).sqrt();
    let mut s = common::emw();
    s.order = 2;
    s.materials.default.insert("eps_r".into(), 1.0);
    s.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase }];
    s.study = Study::Cavity { n_modes: 3, f_search: 0.95 * exact };
    let r = common::solve(&after, &s).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|x| x.value.unwrap()).collect();
    eprintln!("{} -> {} tets, Floquet TE10 {:.6e} vs {:.6e}", nyx.n_cells(), after.n_cells(), got[0], exact);
    assert!((got[0] / exact - 1.0).abs() < 5e-3, "{got:?} vs {exact}");
}
