//! M7 experiment: 3D PML truncation error, native prism/hex layers (P1 x linear / trilinear)
//! against the same layers split into tetrahedra (P1). Each variant is compared with a reference
//! of its own kind (thick layer, many layers) on the same interface mesh, so the physical
//! discretisation error cancels.
//!
//! Usage: `cargo run --release --example pml_3d -- [phys_cells=6] [wavelengths=1]`

use femcore::meshing::pml_box::{PmlSpec, build_pml};
use femcore::mesh::cells::CellMesh;
use femcore::mesh::fixtures::box_tet;
use femcore::meshing::validate::cplx::C64;
use femcore::meshing::validate::helmholtz::{Profile, gaussian_source, rel_l2, solve};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(6);
    let nw: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1.0);
    let k = 2.0 * std::f64::consts::PI;
    let phys = box_tet([[0.0, nw]; 3], [n; 3]);
    let src = gaussian_source([nw / 2.0; 3], 0.12 * nw.max(1.0));
    println!(
        "# 3D PML truncation error (wavelength 1, cube {nw}, {n}^3 physical cells, {} physical vertices)\n",
        phys.n_points()
    );
    let run = |thick: f64, layers: usize, split: bool| -> (Vec<C64>, usize, usize, f64) {
        let pm = build_pml(&phys, &PmlSpec::<3>::uniform(thick, layers)).unwrap();
        let mesh = if split { CellMesh::from_simplex(&pm.mesh.split_to_simplices().unwrap()) } else { pm.mesh.clone() };
        let t = Instant::now();
        let s = solve(&mesh, Some((&pm, Profile::default())), k, &src).unwrap();
        assert!(s.residual < 1e-8, "residual {}", s.residual);
        (s.u, s.n_dof, s.bandwidth, t.elapsed().as_secs_f64())
    };
    let (rq, nd, bw, t) = run(0.6, 6, false);
    println!("reference (native): thickness 0.6, 6 layers, {nd} dof, band {bw}, {t:.1}s");
    let (rt, nd, bw, t) = run(0.6, 6, true);
    println!("reference (split): {nd} dof, band {bw}, {t:.1}s");
    println!("native-reference vs split-reference (physical vertices): {:.2e}\n", rel_l2(&rt, &rq, phys.n_points()));
    println!("| thickness | layers | dof | native err | native s | split err | split s |");
    println!("|---|---|---|---|---|---|---|");
    for thick in [0.3, 0.5] {
        for layers in [1usize, 2, 3, 5] {
            let (uq, nd, _, tq) = run(thick, layers, false);
            let (ut, _, _, tt) = run(thick, layers, true);
            println!(
                "| {thick} | {layers} | {nd} | {:.2e} | {tq:.1} | {:.2e} | {tt:.1} |",
                rel_l2(&uq, &rq, phys.n_points()),
                rel_l2(&ut, &rt, phys.n_points())
            );
        }
    }
}
