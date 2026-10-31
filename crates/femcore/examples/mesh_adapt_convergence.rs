//! M4 study: P1 interpolation error against number of elements for metric-adapted 2D meshes
//! (Hessian-driven, with and without gradation) versus uniform meshes of the same size.
//!
//! Error norms are recomputed here by quadrature on the output mesh, independently of the adaptation code.
//! Usage: `cargo run --release --example adapt_convergence -- [max_target_elements=128000]`

use femcore::meshers::adapt::drive::{AdaptOptions, adapt2};
use femcore::meshers::adapt::field::BackgroundMetric2;
use femcore::meshers::adapt::metric::*;
use femcore::meshers::fixtures::{Field, Gaussian, ObliqueLayer, box_tri, rect_tri};
use femcore::meshers::mesh::SimplexMesh;
use femcore::meshers::quality::{edge_stats, mesh_quality};
use std::time::Instant;

/// 7-point degree-5 rule on the reference triangle (barycentric points and weights summing to 1).
fn rule() -> Vec<([f64; 3], f64)> {
    let a = (6.0 - 15f64.sqrt()) / 21.0;
    let b = (6.0 + 15f64.sqrt()) / 21.0;
    let wa = (155.0 - 15f64.sqrt()) / 1200.0;
    let wb = (155.0 + 15f64.sqrt()) / 1200.0;
    vec![
        ([1.0 / 3.0; 3], 9.0 / 40.0),
        ([1.0 - 2.0 * a, a, a], wa),
        ([a, 1.0 - 2.0 * a, a], wa),
        ([a, a, 1.0 - 2.0 * a], wa),
        ([1.0 - 2.0 * b, b, b], wb),
        ([b, 1.0 - 2.0 * b, b], wb),
        ([b, b, 1.0 - 2.0 * b], wb),
    ]
}

/// (L2, max) error of the P1 interpolant.
fn errors(m: &SimplexMesh<2>, f: &dyn Fn(&[f64; 2]) -> f64) -> (f64, f64) {
    let vals: Vec<f64> = m.points.iter().map(f).collect();
    let r = rule();
    let (mut l2, mut mx) = (0.0, 0.0f64);
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let p: Vec<[f64; 2]> = v.iter().map(|&i| m.points[i as usize]).collect();
        let area = m.measure(c);
        for (l, w) in &r {
            let x =
                [l[0] * p[0][0] + l[1] * p[1][0] + l[2] * p[2][0], l[0] * p[0][1] + l[1] * p[1][1] + l[2] * p[2][1]];
            let li = l[0] * vals[v[0] as usize] + l[1] * vals[v[1] as usize] + l[2] * vals[v[2] as usize];
            let e = li - f(&x);
            l2 += w * area * e * e;
            mx = mx.max(e.abs());
        }
    }
    (l2.sqrt(), mx)
}

fn run(name: &str, field: &dyn Field<2>, hmin: f64, bg_n: usize, max_target: f64) {
    println!("\n### {name}\n");
    println!(
        "| target | variant | elements | adapt s | in-band | q min | q mean | L2 error | max error | uniform L2 | uniform max | L2 gain | max gain | converged |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|");
    let lin = |n: usize| (0..=n).map(|i| i as f64 / n as f64).collect::<Vec<_>>();
    let bg = rect_tri(&lin(bg_n), &lin(bg_n));
    let hess: Vec<Sym<2>> = bg.points.iter().map(|p| field.hessian(p)).collect();
    let vol = bg.vertex_volumes();
    let edges = bg.edges();
    let f = |x: &[f64; 2]| field.value(x);
    for target in [500.0, 2000.0, 8000.0, 32000.0, 128000.0].into_iter().filter(|&t| t <= max_target) {
        for (variant, p, grade_g) in
            [("p=inf", f64::INFINITY, None), ("p=2", 2.0, None), ("p=2 graded", 2.0, Some(0.4))]
        {
            let opt = LpOptions {
                p,
                target_complexity: complexity_for_elements(target, 2),
                hmin,
                hmax: 0.5,
                ..Default::default()
            };
            let mut lp = lp_metrics(&hess, &vol, &opt);
            if let Some(g) = grade_g {
                grade(&mut lp.metrics, &bg.points, &edges, g, 50);
            }
            let field_m = BackgroundMetric2::new(&bg, &lp.metrics);
            let k0 = ((target / 40.0).sqrt() as usize).max(4);
            let start = box_tri([[0.0, 1.0], [0.0, 1.0]], k0, k0);
            let t = Instant::now();
            // runaway guard: six times the element count the (graded) field itself asks for
            let expected = elements_for_complexity(complexity(&lp.metrics, &vol), 2);
            let ao = AdaptOptions { max_cells: Some((6.0 * expected) as usize), ..Default::default() };
            let a = adapt2(&start, &field_m, &ao).unwrap();
            let secs = t.elapsed().as_secs_f64();
            assert!(
                a.metrics.iter().flatten().flatten().all(|x| x.is_finite()),
                "non-finite metric in {name} {variant} {target}"
            );
            let n = a.mesh.n_cells();
            let k = ((n as f64 / 2.0).sqrt().round() as usize).max(1);
            let uni = box_tri([[0.0, 1.0], [0.0, 1.0]], k, k);
            let (l2, mx) = errors(&a.mesh, &f);
            let (ul2, umx) = errors(&uni, &f);
            eprintln!("  [{target} {variant}] {} iterations | {}", a.report.iters.len(), a.report.phases.summary());
            let q = mesh_quality(&a.mesh, Some(&a.metrics));
            let e = edge_stats(&a.mesh, Some(&a.metrics));
            println!(
                "| {target} | {variant} | {n} | {secs:.2} | {:.3} | {:.3} | {:.3} | {l2:.2e} | {mx:.2e} | {ul2:.2e} | {umx:.2e} | {:.1}x | {:.1}x | {} |",
                e.in_band,
                q.min,
                q.mean,
                ul2 / l2,
                umx / mx,
                if a.report.converged {
                    "yes"
                } else if a.report.capped {
                    "CAPPED"
                } else {
                    "no"
                }
            );
        }
    }
}

fn main() {
    let max_target: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(128000.0);
    println!("# P1 interpolation error vs element count (adapted vs uniform with the same number of triangles)");
    run(
        "anisotropic Gaussian, widths (0.12, 0.03)",
        &Gaussian::<2> { centre: [0.5, 0.5], width: [0.12, 0.03] },
        1e-4,
        400,
        max_target,
    );
    run(
        "oblique layer tanh((x + 0.5 y - 0.8)/0.01)",
        &ObliqueLayer::<2> { normal: [1.0, 0.5], offset: 0.8 / 1.25f64.sqrt(), delta: 0.01 / 1.25f64.sqrt() },
        1e-4,
        400,
        max_target,
    );
    run(
        "thin Gaussian ridge, widths (0.3, 0.005)",
        &Gaussian::<2> { centre: [0.5, 0.4], width: [0.3, 0.005] },
        1e-4,
        400,
        max_target,
    );
}
