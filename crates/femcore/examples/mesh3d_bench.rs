// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Timing of the 3D meshing paths that exist before and after ADR 50 to 55, with fixed arguments so that a control run of the parent commit
//! can be compared case by case: the lattice mesher on a sphere in a box and on a layer stack, the sphere cavity of `cavity_bench`, and plain
//! `adapt3` on a block. Prints one line per case: `name seconds cells`. Run three times and take the best of each (`tools/bench_mesh3d.sh`).
//! Does not use anything the parent commit lacks, so the same file builds in both.
use femcore::adapt::drive3::adapt3;
use femcore::adapt::metric::isotropic;
use femcore::adapt::ops3::AdaptOptions3;
use femcore::adapt::topo3::TopoOptions3;
use femcore::adapt::work::simplex_to_adapt;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use std::time::Instant;

fn time<T>(name: &str, f: impl FnOnce() -> (T, usize)) {
    let t = Instant::now();
    let (_, cells) = f();
    println!("{name:28} {:9.3} s {cells:9} cells", t.elapsed().as_secs_f64());
}

fn main() {
    // lattice mesher: a sphere in a box
    let mut ball = Model::new_3d();
    let air = ball.box_("air", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let b = ball.sphere("ball", (2.0, 2.0, 2.0), 0.9);
    ball.overlay(&[("air", air), ("ball", b)]);
    time("lattice sphere in box 0.1", || {
        let g = mesh_model(&ball.list, &MeshPlan::simple(0.1)).unwrap();
        (0, g.mesh.n_cells())
    });
    // lattice mesher: a layer stack with planes on the lattice (the case with slivers before ADR 55)
    let mut stack = Model::new_3d();
    let a = stack.box_("air", (-2.0, -2.0, -2.0), (2.0, 2.0, 2.0));
    let g = stack.box_("glass", (-2.0, -2.0, -2.0), (2.0, 2.0, 0.0));
    stack.overlay(&[("air", a), ("glass", g)]);
    time("lattice layer stack 0.05", || {
        let g = mesh_model(&stack.list, &MeshPlan::simple(0.05)).unwrap();
        (0, g.mesh.n_cells())
    });
    // the cavity of `cavity_bench`
    time("sphere cavity 0.05 / 0.004", || {
        let g = femcore::meshing::geomdrive3::sphere_cavity(0.05, 0.004).unwrap();
        (0, g.mesh.n_cells())
    });
    // the one PML generator on a tetrahedral block (it exists before ADR 52 and gains from the removal of a quadratic copy)
    let block = femcore::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.025, 40).unwrap().mesh;
    time("add_pml_3d on a 40^3 block", || {
        let spec = femcore::meshing::pml_box::PmlSpec::<3> { outer_label: Some(1), ..femcore::meshing::pml_box::PmlSpec::<3>::uniform(0.1, 3) };
        let pb = femcore::meshing::pml3::add_pml_3d(&block, &spec, 3.0, 2.0, 1.0).unwrap();
        (0, pb.mesh.n_cells())
    });
    // plain 3D adaptation of a block towards a graded size
    let start = femcore::testkit::block(1.0, 1.0, 1.0, 6, 6, 6);
    let old = simplex_to_adapt::<3>(&start).unwrap();
    time("adapt3 block graded", || {
        let field = |x: &[f64; 3]| isotropic::<3>(0.04 + 0.2 * x[0]);
        let out = adapt3(&old, &field, &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
        println!("  adapt3: {} iterations, converged {}, qmin {:.3}", out.report.iters.len(), out.report.converged, out.report.iters.last().map_or(0.0, |s| s.q_min));
        (0, out.mesh.n_cells())
    });
}
