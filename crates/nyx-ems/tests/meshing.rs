// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The SOI wire (`femcore::usecases::soi`) through the simple and the adaptive mesher and the adaptive loop, solved by the ems mode solver:
//! moved from the femcore meshing tests when the kernels moved into this crate (ADR 45). The ignored tests are studies (`--ignored --nocapture`).
mod common;
use femcore::adapt::estimate::{solve_adaptively, Estimator, LoopSettings};
use femcore::meshing::{mesh_model, MeshPlan};
use std::collections::BTreeSet;

fn labels(mesh: &femcore::mesh::SimplexMesh) -> (BTreeSet<i32>, BTreeSet<i32>) {
    (mesh.cell_labels.iter().copied().collect(), mesh.facet_labels.iter().copied().collect())
}

/// The real use case: the quasi-TE index of a 500 x 300 nm SOI wire on the simple mesher, on its a priori graded adaptation, and against a
/// reference on a mesh with every size halved. Prints the table; asserts only what must hold (same labels, close to the reference).
#[test]
fn soi_wire_modes_on_the_simple_and_the_adaptive_mesh() {
    use femcore::testkit::soi::SoiWire;
    let wire = SoiWire::default();
    let te = |plan: &MeshPlan, w: &SoiWire| {
        let g = mesh_model(&w.model().list, plan).unwrap();
        let spec = common::wire_problem(w, &g);
        let res = common::solve(&g.mesh, &spec).unwrap();
        let neff = res.solutions.iter().filter_map(|s| s.eigenvalue.map(|n| n[0])).fold(0.0, f64::max);
        (neff, g.mesh.n_cells(), labels(&g.mesh))
    };
    let (n_simple, c_simple, l_simple) = te(&MeshPlan::simple(wire.h_far), &wire);
    let (n_ad1, c_ad1, l_ad1) = te(&MeshPlan::adaptive(wire.h_far, 1.0), &wire);
    let (n_ad2, c_ad2, l_ad2) = te(&MeshPlan::adaptive(wire.h_far, 0.6), &wire);
    let fine =
        SoiWire { h_core: wire.h_core * 0.4, h_near: wire.h_near * 0.4, h_far: wire.h_far * 0.4, ..wire.clone() };
    let (n_ref, c_ref, _) = te(&MeshPlan::simple(fine.h_far), &fine);
    eprintln!("SOI quasi-TE n_eff (reference {n_ref:.6}, {c_ref} cells):");
    eprintln!("  simple            {n_simple:.6}  err {:.1e}  {c_simple} cells", (n_simple - n_ref).abs());
    eprintln!("  adaptive x1.0     {n_ad1:.6}  err {:.1e}  {c_ad1} cells", (n_ad1 - n_ref).abs());
    eprintln!("  adaptive x0.6     {n_ad2:.6}  err {:.1e}  {c_ad2} cells", (n_ad2 - n_ref).abs());
    assert_eq!(l_simple, l_ad1);
    assert_eq!(l_simple, l_ad2);
    assert!((n_ad2 - n_ref).abs() < 2e-3);
}

#[test]
fn the_loop_on_the_soi_wire_converges_n_eff() {
    use femcore::testkit::soi::SoiWire;
    let wire = SoiWire::default();
    let mut st = LoopSettings::new(wire.h_far);
    st.estimator = Estimator { fields: vec!["|E|".into()], ..Estimator::default() };
    st.max_cycles = 4;
    st.rel_tol = 2e-6;
    let mut solve = |g: &femcore::meshing::geomdrive::MeshedGeom| {
        let res = common::solve(&g.mesh, &common::wire_problem(&wire, g))?;
        let n = res.solutions.iter().filter_map(|s| s.eigenvalue.map(|n| n[0])).fold(0.0, f64::max);
        Ok((res, n))
    };
    let out = solve_adaptively(&wire.model().list, &st, &mut solve).unwrap();
    // reference: every size times 0.22 on the simple mesher (16k cells; the study below shows it still drifts by about 1e-5)
    let fine =
        SoiWire { h_core: wire.h_core * 0.22, h_near: wire.h_near * 0.22, h_far: wire.h_far * 0.22, ..wire.clone() };
    let gref = mesh_model(&fine.model().list, &MeshPlan::simple(fine.h_far)).unwrap();
    let rr = common::solve(&gref.mesh, &common::wire_problem(&fine, &gref)).unwrap();
    let nref = rr.solutions.iter().filter_map(|s| s.eigenvalue.map(|n| n[0])).fold(0.0, f64::max);
    eprintln!("SOI loop, reference n_eff {nref:.6} on {} cells", gref.mesh.n_cells());
    for c in &out.cycles {
        eprintln!(
            "  cells {:6}  n_eff {:.6}  err {:.1e}  change {:?}  mesh {:.2}s solve {:.2}s",
            c.cells,
            c.quantity,
            (c.quantity - nref).abs(),
            c.change,
            c.mesh_seconds,
            c.solve_seconds
        );
    }
    let last = out.cycles.last().unwrap();
    assert!((last.quantity - nref).abs() < 2e-5, "{} vs {nref}", last.quantity);
}

/// Study, not a check (run with `--ignored --nocapture`): n_eff of the SOI wire on uniformly scaled simple meshes, to see what the
/// adaptive loop is converging to.
#[test]
#[ignore]
fn soi_uniform_convergence_table() {
    use femcore::testkit::soi::SoiWire;
    let wire = SoiWire::default();
    for f in [1.0, 0.6, 0.4, 0.3, 0.22, 0.16] {
        let w = SoiWire { h_core: wire.h_core * f, h_near: wire.h_near * f, h_far: wire.h_far * f, ..wire.clone() };
        let g = mesh_model(&w.model().list, &MeshPlan::simple(w.h_far)).unwrap();
        let t = std::time::Instant::now();
        let r = common::solve(&g.mesh, &common::wire_problem(&w, &g)).unwrap();
        let n = r.solutions.iter().filter_map(|s| s.eigenvalue.map(|n| n[0])).fold(0.0, f64::max);
        eprintln!(
            "uniform x{f:<5} cells {:6}  n_eff {n:.7}  solve {:.1}s",
            g.mesh.n_cells(),
            t.elapsed().as_secs_f64()
        );
    }
}

/// Study (run with `--ignored --nocapture`): the adaptive loop at a new SOI width, cold and warm-started from the neighbouring width.
#[test]
#[ignore]
fn warm_start_study() {
    use femcore::testkit::soi::SoiWire;
    let run = |w: f64, warm: Option<femcore::adapt::estimate::WarmStart>| {
        let wire = SoiWire { width: w, ..SoiWire::default() };
        let mut st = LoopSettings::new(wire.h_far);
        st.estimator = Estimator { fields: vec!["|E|".into()], ..Estimator::default() };
        st.max_cycles = 5;
        st.rel_tol = 3e-6;
        st.warm = warm;
        let mut solves = 0;
        let mut solve = |g: &femcore::meshing::geomdrive::MeshedGeom| {
            solves += 1;
            let res = common::solve(&g.mesh, &common::wire_problem(&wire, g))?;
            let n = res.solutions.iter().filter_map(|s| s.eigenvalue.map(|n| n[0])).fold(0.0, f64::max);
            Ok((res, n))
        };
        let t = std::time::Instant::now();
        let out = solve_adaptively(&wire.model().list, &st, &mut solve).unwrap();
        (out, solves, t.elapsed().as_secs_f64())
    };
    let (a, sa, ta) = run(450e-9, None);
    eprintln!(
        "w=450 cold: {sa} solves, {:.1}s, n_eff {:.7}, {} cells",
        ta,
        a.cycles.last().unwrap().quantity,
        a.geom.mesh.n_cells()
    );
    let (b, sb, tb) = run(470e-9, None);
    eprintln!(
        "w=470 cold: {sb} solves, {:.1}s, n_eff {:.7}, {} cells",
        tb,
        b.cycles.last().unwrap().quantity,
        b.geom.mesh.n_cells()
    );
    let (c, sc, tc) = run(470e-9, Some(a.warm.clone()));
    for cy in &c.cycles {
        eprintln!("    warm cycle: cells {} n_eff {:.7} change {:?}", cy.cells, cy.quantity, cy.change);
    }
    eprintln!(
        "w=470 warm: {sc} solves, {:.1}s, n_eff {:.7}, {} cells",
        tc,
        c.cycles.last().unwrap().quantity,
        c.geom.mesh.n_cells()
    );
}
