//! Timing of the complex shift-invert stage on a rectangular guide with an absorbing layer (`ems_pml`), at a chosen element size:
//!
//!   cargo run --release -p femcore --example complex_bench -- <h in mm> [modes]
//!
//! Prints the number of unknowns, the solve time and the effective indices. Used to compare the sequential and the column-parallel
//! solves of the new driver (`NYX_PAR_SOLVE_MIN_N=0` or `=5000`) with the old stack (build this file at the commit `restructure-foundation`).
use femcore::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::problem::*;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let h: f64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(2.5) * 1e-3;
    let n_modes: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let (lx, ly, xc) = (0.03, 0.015, 0.02);
    let mut m = geom::model::Model::new_2d();
    let ids: Vec<(&str, geom::model::FeatureId)> = vec![("a", m.rect("a", (0.0, 0.0), (xc, ly))), ("b", m.rect("b", (xc, 0.0), (lx, ly)))];
    m.overlay(&ids);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h)).unwrap();
    let (l, b, t, r) = (g.select("plane(x = 0)").unwrap(), g.select("plane(z = 0)").unwrap(), g.select(&format!("plane(z = {ly:e})")).unwrap(), g.select(&format!("plane(x = {lx:e})")).unwrap());
    let mut mesh = g.mesh;
    let old = mesh.facet_labels.clone();
    for (f, lab) in mesh.facet_labels.iter_mut().enumerate() {
        *lab = if l.contains(&old[f]) { 1 } else if b.contains(&old[f]) { 2 } else if t.contains(&old[f]) { 3 } else if r.contains(&old[f]) { 4 } else { old[f] };
    }
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.set(1, "eps_r", 4.0);
    s.materials.set(2, "eps_r", 1.0);
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Pml::default() });
    s.study = Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes };
    let n = femcore::ems::n_unknowns(&mesh, 2);
    let mut best = f64::INFINITY;
    let mut neffs = String::new();
    for _ in 0..3 {
        let t0 = Instant::now();
        let res = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
        best = best.min(t0.elapsed().as_secs_f64());
        neffs = res.solutions.iter().map(|x| format!("{:.6}", x.neff.map_or(0.0, |z| z.0))).collect::<Vec<_>>().join(" ");
    }
    println!("h {:.2} mm: {} cells, about {} unknowns, {} modes, best of 3: {:.3} s, neff {}", h * 1e3, mesh.n_cells(), n, n_modes, best, neffs);
}
