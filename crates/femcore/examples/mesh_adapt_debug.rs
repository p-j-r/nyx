// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Single adaptation run with a per-iteration log, for diagnosing the cases of `adapt_convergence`.
//! Usage: `cargo run --release --example adapt_debug -- <p|inf> <target> <gauss|layer|ridge> [grade_g]`

use femcore::adapt::drive::{AdaptOptions, adapt2};
use femcore::adapt::field::BackgroundMetric2;
use femcore::adapt::metric::*;
use femcore::mesh::fixtures::{Field, Gaussian, ObliqueLayer, box_tri, rect_tri};
use femcore::mesh::quality::{edge_stats, mesh_quality};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let p = if a[1] == "inf" { f64::INFINITY } else { a[1].parse().unwrap() };
    let target: f64 = a[2].parse().unwrap();
    let g: Option<f64> = a.get(4).and_then(|s| s.parse().ok());
    let field: Box<dyn Field<2>> = match a[3].as_str() {
        "gauss" => Box::new(Gaussian::<2> { centre: [0.5, 0.5], width: [0.12, 0.03] }),
        "layer" => Box::new(ObliqueLayer::<2> {
            normal: [1.0, 0.5],
            offset: 0.8 / 1.25f64.sqrt(),
            delta: 0.01 / 1.25f64.sqrt(),
        }),
        _ => Box::new(Gaussian::<2> { centre: [0.5, 0.4], width: [0.3, 0.005] }),
    };
    let lin = |n: usize| (0..=n).map(|i| i as f64 / n as f64).collect::<Vec<_>>();
    let bg = rect_tri(&lin(400), &lin(400));
    let hess: Vec<Sym<2>> = bg.points.iter().map(|q| field.hessian(q)).collect();
    let opt = LpOptions {
        p,
        target_complexity: complexity_for_elements(target, 2),
        hmin: 1e-4,
        hmax: 0.5,
        ..Default::default()
    };
    let mut lp = lp_metrics(&hess, &bg.vertex_volumes(), &opt);
    println!("scale {:.3e}, complexity {:.3}, saturated {}", lp.scale, lp.complexity, lp.saturated);
    if let Some(g) = g {
        grade(&mut lp.metrics, &bg.points, &bg.edges(), g, 50);
    }
    let fm = BackgroundMetric2::new(&bg, &lp.metrics);
    let k0 = ((target / 40.0).sqrt() as usize).max(4);
    let t0 = std::time::Instant::now();
    let r = adapt2(
        &box_tri([[0.0, 1.0], [0.0, 1.0]], k0, k0),
        &fm,
        &AdaptOptions { max_cells: Some(40_000), ..Default::default() },
    )
    .unwrap();
    println!("converged {} capped {}", r.report.converged, r.report.capped);
    println!("adapt time {:.2}s | {}", t0.elapsed().as_secs_f64(), r.report.phases.summary());
    for it in &r.report.iters {
        println!(
            "iter {:2}: {:6} tris +{} -{} flips {} moves {} in-band {:.3} qmin {:.3}",
            it.iter, it.n_triangles, it.splits, it.collapses, it.flips, it.moves, it.in_band, it.q_min
        );
    }
    let q = mesh_quality(&r.mesh, Some(&r.metrics));
    let e = edge_stats(&r.mesh, Some(&r.metrics));
    println!(
        "final {} tris, q min {:.3e}, edge {:.3}..{:.3}, in-band {:.3}, hist {:?}",
        r.mesh.n_cells(),
        q.min,
        e.min,
        e.max,
        e.in_band,
        e.histogram
    );
    // worst triangle
    let mut worst = (f64::INFINITY, 0);
    for c in 0..r.mesh.n_cells() {
        let v = r.mesh.cell(c);
        let pts: Vec<[f64; 2]> = v.iter().map(|&i| r.mesh.points[i as usize]).collect();
        let ms = [r.metrics[v[0] as usize], r.metrics[v[1] as usize], r.metrics[v[2] as usize]];
        let qq = femcore::mesh::quality::simplex_quality::<2>(&pts, Some(&ms));
        if qq < worst.0 {
            worst = (qq, c);
        }
    }
    let v = r.mesh.cell(worst.1);
    println!("worst tri {} q {:.3e}:", worst.1, worst.0);
    for &i in v {
        println!("   vertex {} at {:?}, sizes {:?}", i, r.mesh.points[i as usize], sizes(&r.metrics[i as usize]).0);
    }
}
