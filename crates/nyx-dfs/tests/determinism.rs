// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tier A (`<= 60k` DOFs) subprocess `0-ULP` (`to_bits()`) determinism matrix for `nyx-dfs`:
//! verifies `NYX_THREADS=1` vs `NYX_THREADS=4` bit-exact equality across electrostatic (3D AMG-CG,
//! 3D direct, and 2D layered `planar_bem`), heat (3D stationary AMG-CG/direct and transient BDF2),
//! and magnetostatic (`cross_section` and `planar_stream`) solves.

mod common;

use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::props::LinearSolver;
use geom::model::{Model, RegionKind};
use nyx_dfs::planar_bem::{solve_planar_bem_with_stack, LayeredStack};
use nyx_dfs::{solve_electrostatics, solve_heat, solve_magnetostatics, Bc, BcKind, Problem, Study};
use nyx_physics::{Scalar, SolveResult};
use std::f64::consts::PI;

fn push_solve_bits(out: &mut Vec<u64>, r: &SolveResult) {
    assert!(r.dofs <= 60_000, "Tier A problem must stay <= 60k DOFs, got {}", r.dofs);
    out.push(r.dofs as u64);
    out.push(r.solutions.len() as u64);
    for sol in &r.solutions {
        out.push(sol.fields.len() as u64);
        for f in &sol.fields {
            out.push(f.values.len() as u64);
            out.extend(f.values.iter().map(|v| v.to_bits()));
        }
        out.push(sol.cell_fields.len() as u64);
        for cf in &sol.cell_fields {
            out.extend(cf.re.iter().map(|v| v.to_bits()));
            out.extend(cf.im.iter().map(|v| v.to_bits()));
        }
    }
    for (_, s) in &r.scalars {
        if let Scalar::Num(v) = s {
            out.push(v.to_bits());
        }
    }
    if let Some(ref d) = r.diagnostics {
        out.push(d.max_algebraic_residual.to_bits());
    }
}

fn electro_3d_bits(solver: LinearSolver) -> Vec<u64> {
    // 10 x 8 x 6 hexes = 2880+ tets (> 2000 cells, triggering parallel assembly and AMG/direct paths)
    let m = common::block(1.0, 0.8, 0.5, 10, 8, 6);
    let mut p = Problem::new(Study::Electrostatic, 2);
    p.solver = solver;
    p.tol = 1e-11;
    p.materials.default.insert("eps_r".into(), 3.9);
    p.bcs = vec![Bc { labels: vec![3], kind: BcKind::Ground }, Bc { labels: vec![4], kind: BcKind::Potential(2.5) }];
    let r = common::solve(&m, &p);
    let mut out = Vec::new();
    push_solve_bits(&mut out, &r);
    out
}

fn planar_bem_bits() -> Vec<u64> {
    // >= 48 panels to engage parallel BEM matrix assembly in `planar_bem.rs`
    let mut m = Model::new_2d();
    let s1 = m.rect("s1", (-120e-6, -80e-6), (-20e-6, 80e-6));
    let s2 = m.rect("s2", (20e-6, -80e-6), (120e-6, 80e-6));
    m.overlay(&[("s1", s1), ("s2", s2)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(22e-6)).unwrap();
    assert!(g.mesh.n_cells() >= 48, "expected >= 48 triangles for parallel BEM");

    let (mut l1, mut l2) = (Vec::new(), Vec::new());
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        let xm = 0.5 * (g.mesh.point(fv[0] as usize)[0] + g.mesh.point(fv[1] as usize)[0]);
        let lab = g.mesh.facet_labels[f];
        if xm < 0.0 && !l1.contains(&lab) {
            l1.push(lab);
        } else if xm > 0.0 && !l2.contains(&lab) {
            l2.push(lab);
        }
    }
    let mut p = Problem::new(Study::Electrostatic, 1).with_formulation("planar_bem");
    p.materials.default.insert("eps_r".into(), 6.225);
    p.bcs = vec![
        Bc { labels: l1, kind: BcKind::Terminal("pad_1".into()) },
        Bc { labels: l2, kind: BcKind::Terminal("pad_2".into()) },
    ];
    let stack = LayeredStack { eps_top: 1.0, eps_sub: 11.45, eps_bot: 1.0, t_sub: 300e-6, grounded_backplane: true };
    let r = solve_planar_bem_with_stack(&g.mesh, &p, Some(stack)).unwrap();
    let mut out = Vec::new();
    push_solve_bits(&mut out, &r);
    out
}

fn heat_stationary_3d_bits(solver: LinearSolver) -> Vec<u64> {
    let m = common::block(1.0, 0.6, 0.8, 10, 6, 8);
    let mut p = Problem::new(Study::HeatStationary, 2);
    p.solver = solver;
    p.tol = 1e-11;
    p.materials.default.insert("k".into(), 3.2);
    p.materials.default.insert("Q".into(), 18.0);
    p.bcs = vec![
        Bc { labels: vec![3], kind: BcKind::Temperature(300.0) },
        Bc { labels: vec![4], kind: BcKind::Convection { h: 6.0, text: 293.0 } },
    ];
    let r = solve_heat(&m, &p).unwrap();
    let mut out = Vec::new();
    push_solve_bits(&mut out, &r);
    out
}

fn heat_transient_bits(solver: LinearSolver) -> Vec<u64> {
    let m = common::block(1.0, 0.4, 0.6, 10, 4, 6);
    let mut p =
        Problem::new(Study::HeatTransient { t_end: 0.4, n_out: 4, t0: "sin(pi * x) * sin(pi * z / 0.6)".into() }, 2);
    p.solver = solver;
    p.tol = 1e-11;
    p.materials.default.insert("k".into(), 1.8);
    p.materials.default.insert("rho".into(), 1.5);
    p.materials.default.insert("Cp".into(), 0.9);
    p.bcs = vec![Bc { labels: vec![1, 3, 4, 6], kind: BcKind::Temperature(0.0) }];
    let r = solve_heat(&m, &p).unwrap();
    let mut out = Vec::new();
    push_solve_bits(&mut out, &r);
    out
}

fn magneto_bits() -> Vec<u64> {
    let mut out = Vec::new();

    // 1. 2D coaxial cable (`cross_section`)
    let (a, b) = (0.2_f64, 0.8_f64);
    let mut m_coax = Model::new_2d();
    let outer = m_coax.circle("diel", (0.0, 0.0), b);
    let inner = m_coax.circle("core", (0.0, 0.0), a);
    m_coax.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let g_coax = mesh_geom_2d(evaluate(&m_coax.list).unwrap(), &GeomMeshOpts::new(0.04)).unwrap();
    let r_mid = 0.5 * (a + b);
    let (mut in_labs, mut out_labs) = (Vec::new(), Vec::new());
    for f in 0..g_coax.mesh.n_facets() {
        let fv = g_coax.mesh.facet(f);
        let p0 = g_coax.mesh.point(fv[0] as usize);
        let p1 = g_coax.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g_coax.mesh.facet_labels[f];
        if rm < r_mid && !in_labs.contains(&lab) {
            in_labs.push(lab);
        } else if rm >= r_mid && !out_labs.contains(&lab) {
            out_labs.push(lab);
        }
    }
    let mut p_coax = Problem::new(Study::Magnetostatic { formulation: "cross_section".into() }, 2);
    p_coax.bcs = vec![
        Bc {
            labels: in_labs,
            kind: BcKind::CurrentTerminal { name: "center".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc { labels: out_labs, kind: BcKind::Ground },
    ];
    let r_coax = solve_magnetostatics(&g_coax.mesh, &p_coax).unwrap();
    push_solve_bits(&mut out, &r_coax);

    // 2. 2D superconducting ring (`planar_stream`)
    let (r_in, r_out) = (90e-6_f64, 110e-6_f64);
    let r_loop = 0.5 * (r_in + r_out);
    let mut m_loop = Model::new_2d();
    let outer_c = m_loop.circle("washer", (0.0, 0.0), r_out);
    let inner_c = m_loop.circle("hole", (0.0, 0.0), r_in);
    m_loop.overlay_full(
        &[("washer", outer_c, RegionKind::Material, None), ("hole", inner_c, RegionKind::Hole, None)],
        &[],
    );
    let g_loop = mesh_geom_2d(evaluate(&m_loop.list).unwrap(), &GeomMeshOpts::new(4.0e-6)).unwrap();
    let (mut hole_labs, mut outer_labs) = (Vec::new(), Vec::new());
    for f in 0..g_loop.mesh.n_facets() {
        let fv = g_loop.mesh.facet(f);
        let p0 = g_loop.mesh.point(fv[0] as usize);
        let p1 = g_loop.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g_loop.mesh.facet_labels[f];
        if rm < r_loop && !hole_labs.contains(&lab) {
            hole_labs.push(lab);
        } else if rm >= r_loop && !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }
    let mut p_ps = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_ps.bcs = vec![
        Bc { labels: outer_labs, kind: BcKind::Ground },
        Bc { labels: hole_labs, kind: BcKind::FluxLoop { name: "loop".into(), flux_wb: 1.0 } },
    ];
    let r_ps = solve_magnetostatics(&g_loop.mesh, &p_ps).unwrap();
    push_solve_bits(&mut out, &r_ps);

    // 3. Direct multi-terminal 3D electrostatic check via solve_electrostatics too
    let m3 = SimplexMesh::new(
        3,
        common::block(0.5, 0.5, 0.5, 6, 6, 6).points,
        common::block(0.5, 0.5, 0.5, 6, 6, 6).cells,
        common::block(0.5, 0.5, 0.5, 6, 6, 6).cell_labels,
        common::block(0.5, 0.5, 0.5, 6, 6, 6).facets,
        common::block(0.5, 0.5, 0.5, 6, 6, 6).facet_labels,
    );
    let mut p_e = Problem::new(Study::Electrostatic, 2);
    p_e.materials.default.insert("eps_r".into(), 2.0 + 0.0 * PI);
    p_e.bcs = vec![Bc { labels: vec![3], kind: BcKind::Ground }, Bc { labels: vec![4], kind: BcKind::Potential(1.0) }];
    push_solve_bits(&mut out, &solve_electrostatics(&m3, &p_e).unwrap());

    out
}

#[test]
fn dfs_electrostatic_heat_and_magnetostatic_are_bit_exact_across_nyx_threads() {
    if let Ok(out_path) = std::env::var("NYX_DFS_DET_OUT") {
        let mut all = Vec::new();
        all.extend_from_slice(&electro_3d_bits(LinearSolver::Iterative));
        all.extend_from_slice(&electro_3d_bits(LinearSolver::Direct));
        all.extend_from_slice(&planar_bem_bits());
        all.extend_from_slice(&heat_stationary_3d_bits(LinearSolver::Iterative));
        all.extend_from_slice(&heat_stationary_3d_bits(LinearSolver::Direct));
        all.extend_from_slice(&heat_transient_bits(LinearSolver::Iterative));
        all.extend_from_slice(&heat_transient_bits(LinearSolver::Direct));
        all.extend_from_slice(&magneto_bits());
        let bytes: Vec<u8> = all.iter().flat_map(|w| w.to_le_bytes()).collect();
        std::fs::write(&out_path, bytes).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-dfs-det-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args([
                "--exact",
                "dfs_electrostatic_heat_and_magnetostatic_are_bit_exact_across_nyx_threads",
                "--nocapture",
            ])
            .env("NYX_THREADS", nt)
            .env_remove("NYX_DETERMINISTIC_FACTOR")
            .env("NYX_DFS_DET_OUT", &out_path)
            .output()
            .expect("spawn dfs determinism subprocess");
        assert!(
            out.status.success(),
            "child NYX_THREADS={nt} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        bytes
    };

    let b1 = run_child("1");
    let b4 = run_child("4");
    assert!(!b1.is_empty());
    assert_eq!(b1, b4, "nyx-dfs Tier A outputs must be 0-ULP bit-identical for NYX_THREADS=1 and 4");
}
