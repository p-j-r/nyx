// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Poisson on a disk meshed by the unified meshing interface: the adaptive mesher and the solve-estimate-adapt loop improve the electrostatic
//! solution. (Moved from the femcore meshing tests with the dfs kernel.)
use femcore::adapt::estimate::{solve_adaptively, LoopSettings};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_dfs::{Bc, BcKind, Problem, Study};

fn disk(r: f64) -> Model {
    let mut m = Model::new_2d();
    let d = m.circle("disk", (0.0, 0.0), r);
    m.overlay(&[("disk", d)]);
    m
}

#[test]
fn poisson_on_the_adapted_disk_is_at_least_as_accurate_as_on_the_simple_one() {
    // -div(eps0 grad V) = eps0 on a disk of radius R, V = 0 on the rim: V max = R^2/4
    let r = 0.5;
    let m = disk(r);
    let run = |plan: &MeshPlan| {
        let g = mesh_model(&m.list, plan).unwrap();
        let mut s = Problem::new(Study::Electrostatic, 2);
        s.materials.default.insert("rho_v".into(), nyx_dfs::EPS0);
        s.bcs = vec![Bc { labels: g.select("outer").unwrap(), kind: BcKind::Ground }];
        let res = nyx_dfs::solve_electrostatics(&g.mesh, &s).unwrap();
        ((res.scalar("V max").unwrap() / (r * r / 4.0) - 1.0).abs(), g.mesh.n_cells())
    };
    let (e_simple, n_simple) = run(&MeshPlan::simple(0.15));
    let (e_adapt, n_adapt) = run(&MeshPlan::adaptive(0.15, 0.5));
    let (e_fine, n_fine) = run(&MeshPlan::simple(0.15 * 0.5));
    eprintln!("poisson disk: simple h=0.15 {e_simple:.2e} ({n_simple} cells); adaptive 0.5x {e_adapt:.2e} ({n_adapt}); simple h=0.075 {e_fine:.2e} ({n_fine})");
    assert!(n_adapt > n_simple);
    assert!(e_adapt < e_simple, "adapted {e_adapt} ({n_adapt} cells) vs simple {e_simple} ({n_simple} cells)");
}

#[test]
fn the_loop_improves_poisson_on_a_disk_and_reports_its_cycles() {
    let r = 0.5;
    let m = disk(r);
    let mut st = LoopSettings::new(0.25);
    st.max_cycles = 4;
    st.rel_tol = 1e-6;
    let mut solve = |g: &femcore::meshing::geomdrive::MeshedGeom| {
        let mut s = Problem::new(Study::Electrostatic, 2);
        s.materials.default.insert("rho_v".into(), nyx_dfs::EPS0);
        s.bcs = vec![Bc { labels: g.select("outer").unwrap(), kind: BcKind::Ground }];
        let res = nyx_dfs::solve_electrostatics(&g.mesh, &s)?;
        let v = res.scalar("V max").ok_or("no V max")?;
        Ok((res, v))
    };
    let out = solve_adaptively(&m.list, &st, &mut solve).unwrap();
    let exact = r * r / 4.0;
    for c in &out.cycles {
        eprintln!(
            "  cells {:6}  Vmax err {:.2e}  change {:?}  mesh {:.2}s solve {:.2}s",
            c.cells,
            (c.quantity / exact - 1.0).abs(),
            c.change,
            c.mesh_seconds,
            c.solve_seconds
        );
    }
    assert!(out.cycles.len() >= 2);
    let first = (out.cycles[0].quantity / exact - 1.0).abs();
    let last = (out.cycles.last().unwrap().quantity / exact - 1.0).abs();
    assert!(last < first, "{last} vs {first}");
    assert!(out.cycles.last().unwrap().cells > out.cycles[0].cells);
}
