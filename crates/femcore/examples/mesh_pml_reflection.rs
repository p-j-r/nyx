//! M2 experiment: PML truncation error vs number of layers, structured quad layers (Q1) against
//! the same layers split into triangles (P1), on identical physical meshes. Quad runs are compared
//! with a quad reference and split runs with a split reference.
//!
//! Error = relative L2 difference over the physical vertices from a reference solution on the same
//! interface mesh with a layer twice as thick and many more layers, so that the discretisation
//! error of the physical domain cancels and only truncation / layer discretisation remain.
//!
//! Usage: `cargo run --release --example pml_reflection -- [wavelengths=3] [points_per_wavelength=12]`

use femcore::meshers::adapt::pml::{PmlSpec, build_pml};
use femcore::meshers::cellmesh::CellMesh;
use femcore::meshers::fixtures::box_tri;
use femcore::meshers::validate::cplx::C64;
use femcore::meshers::validate::helmholtz::{Profile, Solution, gaussian_source, rel_l2, solve};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nw: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(3.0);
    let ppw: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(12);
    let k = 2.0 * std::f64::consts::PI; // wavelength 1
    let n = (nw * ppw as f64).round() as usize;
    let phys = box_tri([[0.0, nw], [0.0, nw]], n, n);
    let src = gaussian_source([nw / 2.0, nw / 2.0], 0.15);
    println!(
        "# PML truncation error vs layers (wavelength 1, domain {nw}x{nw}, {ppw} pts/wavelength, {} physical vertices)",
        phys.n_points()
    );

    let run_with = |thick: f64, layers: usize, grading: f64, split: bool, r0: f64| -> (Vec<C64>, Solution, f64) {
        let mut spec = PmlSpec::<2>::uniform(thick, layers);
        spec.grading = grading;
        let pm = build_pml(&phys, &spec).unwrap();
        let mesh = if split { CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap()) } else { pm.mesh.clone() };
        let t = Instant::now();
        let s = solve(&mesh, Some((&pm, Profile { r0, ..Default::default() })), k, &src).unwrap();
        let dt = t.elapsed().as_secs_f64();
        (s.u.clone(), s, dt)
    };

    let run = |thick: f64, layers: usize, grading: f64, split: bool| run_with(thick, layers, grading, split, 1e-8);
    let (reference, rs, rt) = run(2.0, 48, 1.0, false);
    println!(
        "reference: thickness 2, 48 layers, {} dof, band {}, {:.2}s, residual {:.1e}",
        rs.n_dof, rs.bandwidth, rt, rs.residual
    );
    // Each discretisation is also compared with a reference of its own kind, so that the layer
    // discretisation of the reference cannot favour one of them.
    let (reference_tri, _, _) = run(2.0, 48, 1.0, true);
    println!(
        "quad-reference vs tri-reference (physical vertices): {:.2e}\n",
        rel_l2(&reference_tri, &reference, phys.n_points())
    );
    println!("| thickness | layers | quad dof | quad err | quad s | tri err | tri s |");
    println!("|---|---|---|---|---|---|---|");
    for thick in [0.5, 1.0] {
        for layers in [1usize, 2, 3, 4, 6, 8, 12] {
            let (uq, sq, tq) = run(thick, layers, 1.0, false);
            let (ut, _st, tt) = run(thick, layers, 1.0, true);
            println!(
                "| {thick} | {layers} | {} | {:.2e} | {:.2} | {:.2e} | {:.2} |",
                sq.n_dof,
                rel_l2(&uq, &reference, phys.n_points()),
                tq,
                rel_l2(&ut, &reference_tri, phys.n_points()),
                tt
            );
        }
    }

    println!("\n## Tuned: best error over design reflection R0 in 1e-1..1e-8 and grading in {{1, 4}} (thickness 1)\n");
    println!("| layers | dof | quad best err | quad (R0, grading) | tri best err | tri (R0, grading) |");
    println!("|---|---|---|---|---|---|");
    for layers in [2usize, 3, 4, 6, 8, 12] {
        let mut best = [(f64::INFINITY, 0.0, 0.0); 2];
        let mut dof = 0;
        for (mi, split) in [false, true].into_iter().enumerate() {
            for r0 in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8] {
                for grading in [1.0, 4.0] {
                    let (u, s, _) = run_with(1.0, layers, grading, split, r0);
                    dof = s.n_dof;
                    let e = rel_l2(&u, if split { &reference_tri } else { &reference }, phys.n_points());
                    if e < best[mi].0 {
                        best[mi] = (e, r0, grading);
                    }
                }
            }
        }
        println!(
            "| {layers} | {dof} | {:.2e} | ({:.0e}, {}) | {:.2e} | ({:.0e}, {}) |",
            best[0].0, best[0].1, best[0].2, best[1].0, best[1].1, best[1].2
        );
    }
}
