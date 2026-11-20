// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Replays the driver loop on one case and reports the first operator pass that creates a
//! triangle of metric quality below a threshold. Usage: `adapt_trace <p|inf> <target> [grade_g]`

use femcore::adapt::drive::AdaptOptions;
use femcore::adapt::field::BackgroundMetric2;
use femcore::adapt::geom::PolylineGeometry;
use femcore::adapt::metric::*;
use femcore::adapt::ops::*;
use femcore::adapt::topo::Work2;
use femcore::mesh::fixtures::{Field, Gaussian, box_tri, rect_tri};

fn minq(w: &Work2) -> (f64, usize) {
    let mut m = (f64::INFINITY, 0);
    for (t, tr) in w.tris.iter().enumerate() {
        if !w.talive[t] {
            continue;
        }
        let p = |i: usize| w.pts[tr[i] as usize];
        let q =
            tri_quality([p(0), p(1), p(2)], [&w.met[tr[0] as usize], &w.met[tr[1] as usize], &w.met[tr[2] as usize]]);
        if q < m.0 {
            m = (q, t);
        }
    }
    m
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let p = if a[1] == "inf" { f64::INFINITY } else { a[1].parse().unwrap() };
    let target: f64 = a[2].parse().unwrap();
    let g: Option<f64> = a.get(3).and_then(|s| s.parse().ok());
    let field = Gaussian::<2> { centre: [0.5, 0.5], width: [0.12, 0.03] };
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
    if let Some(g) = g {
        grade(&mut lp.metrics, &bg.points, &bg.edges(), g, 50);
    }
    let fm = BackgroundMetric2::new(&bg, &lp.metrics);
    let k0 = ((target / 40.0).sqrt() as usize).max(4);
    let ao = AdaptOptions::default();
    let mut w = Work2::from_mesh(&box_tri([[0.0, 1.0], [0.0, 1.0]], k0, k0), &fm, &ao.topo).unwrap();
    let geom = PolylineGeometry::from_work(&w);
    let mut prev = minq(&w).0;
    let mut step = 0;
    let mut log = |name: &str, w: &Work2, n: usize| {
        step += 1;
        let (q, t) = minq(w);
        if q < 1e-3 && prev >= 1e-3 || (step <= 12) {
            let tr = w.tris[t];
            println!("step {step:3} {name:8} n={n:6}  minq {q:.3e} (was {prev:.3e}) tri {t} verts {:?}", tr);
        }
        prev = q;
    };
    for it in 0..4 {
        for _ in 0..30 {
            let n = split_pass(&mut w, &fm, &geom, &ao, it < 2);
            log("split", &w, n);
            if n == 0 {
                break;
            }
        }
        for _ in 0..10 {
            let n = collapse_pass(&mut w, &ao);
            log("collapse", &w, n);
            if n == 0 {
                break;
            }
        }
        for _ in 0..10 {
            let n = flip_pass(&mut w);
            log("flip", &w, n);
            if n == 0 {
                break;
            }
        }
        for _ in 0..3 {
            let n = smooth_pass(&mut w, &fm, &geom, &ao);
            log("smooth", &w, n);
        }
        println!("-- iteration {it} done");
    }
}
