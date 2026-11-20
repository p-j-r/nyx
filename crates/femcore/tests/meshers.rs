// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! End-to-end tests of the 2D adaptation driver against analytic metrics. All quality, length
//! and conservation figures are recomputed here from the output coordinates, independently of
//! the code under test.

use femcore::adapt::drive::{AdaptOptions, adapt2, adapt2_with};
use femcore::adapt::field::{BackgroundMetric2, MetricField};
use femcore::adapt::geom::{PolylineGeometry, circle_projector};
use femcore::adapt::metric::*;
use femcore::mesh::fixtures::{Field, Gaussian, box_tri, rect_tri, side_label};
use femcore::adapt::work::SimplexMesh;
use femcore::mesh::quality::{edge_stats, mesh_quality};

fn iso(h: f64) -> impl Fn(&[f64; 2]) -> Sym<2> + Sync {
    move |_| isotropic::<2>(h)
}

fn side_length(m: &SimplexMesh<2>, label: i32) -> f64 {
    (0..m.n_facets())
        .filter(|&f| m.facet_label[f] == label)
        .map(|f| {
            let v = m.facet(f);
            let (a, b) = (m.points[v[0] as usize], m.points[v[1] as usize]);
            (a[0] - b[0]).hypot(a[1] - b[1])
        })
        .sum()
}

fn label_area(m: &SimplexMesh<2>, label: i32) -> f64 {
    (0..m.n_cells()).filter(|&c| m.cell_label[c] == label).map(|c| m.measure(c)).sum()
}

fn report(name: &str, a: &femcore::adapt::drive::Adapted) {
    let q = mesh_quality(&a.mesh, Some(&a.metrics));
    let e = edge_stats(&a.mesh, Some(&a.metrics));
    println!(
        "{name}: {} vertices {} triangles, {} iterations, in-band {:.3}, q min {:.3} mean {:.3}, edge length {:.2}..{:.2} | {}",
        a.mesh.n_points(),
        a.mesh.n_cells(),
        a.report.iters.len(),
        e.in_band,
        q.min,
        q.mean,
        e.min,
        e.max,
        a.report.phases.summary()
    );
}

#[test]
fn uniform_refinement_of_a_coarse_mesh() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 4, 4);
    let h = 0.05;
    let a = adapt2(&m, &iso(h), &AdaptOptions::default()).unwrap();
    report("uniform h=0.05", &a);
    let expected = 1.0 / (unit_simplex_measure(2) * h * h);
    let n = a.mesh.n_cells() as f64;
    assert!(n > 0.8 * expected && n < 1.25 * expected, "{n} vs {expected}");
    let e = edge_stats(&a.mesh, Some(&a.metrics));
    assert!(e.in_band > 0.85, "in-band {}", e.in_band);
    assert!(mesh_quality(&a.mesh, Some(&a.metrics)).min > 0.3);
    assert!((a.mesh.total_measure() - 1.0).abs() < 1e-12);
    for ax in 0..2 {
        for s in [false, true] {
            assert!((side_length(&a.mesh, side_label(ax, s)) - 1.0).abs() < 1e-12);
        }
    }
    // the four corners are still vertices and nothing lies outside the box
    for c in [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]] {
        assert!(a.mesh.points.contains(&c));
    }
    assert!(a.mesh.points.iter().all(|p| (0.0..=1.0).contains(&p[0]) && (0.0..=1.0).contains(&p[1])));
}

#[test]
fn coarsening_of_a_fine_mesh() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 40, 40);
    let h = 0.1;
    let a = adapt2(&m, &iso(h), &AdaptOptions::default()).unwrap();
    report("coarsen to h=0.1", &a);
    let expected = 1.0 / (unit_simplex_measure(2) * h * h);
    let n = a.mesh.n_cells() as f64;
    assert!(n > 0.7 * expected && n < 1.5 * expected, "{n} vs {expected}");
    assert!(edge_stats(&a.mesh, Some(&a.metrics)).in_band > 0.8);
    assert!((a.mesh.total_measure() - 1.0).abs() < 1e-12);
}

#[test]
fn graded_isotropic_size() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5);
    let f = |x: &[f64; 2]| isotropic::<2>(0.01 + 0.09 * x[0]);
    let a = adapt2(&m, &f, &AdaptOptions::default()).unwrap();
    report("graded", &a);
    assert!(edge_stats(&a.mesh, Some(&a.metrics)).in_band > 0.8);
    assert!(mesh_quality(&a.mesh, Some(&a.metrics)).min > 0.25);
    // many more vertices near x = 0 than near x = 1
    let left = a.mesh.points.iter().filter(|p| p[0] < 0.2).count();
    let right = a.mesh.points.iter().filter(|p| p[0] > 0.8).count();
    assert!(left > 5 * right, "{left} vs {right}");
}

#[test]
fn anisotropic_boundary_layer() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 8, 8);
    // size 0.1 along x everywhere; along y it grows from 0.004 at y = 0 to 0.1
    let f = |x: &[f64; 2]| from_sizes(&[0.1, 0.004 + 0.096 * x[1]], &identity::<2>());
    let a = adapt2(&m, &f, &AdaptOptions::default()).unwrap();
    report("boundary layer", &a);
    let e = edge_stats(&a.mesh, Some(&a.metrics));
    assert!(e.in_band > 0.75, "in-band {}", e.in_band);
    assert!(mesh_quality(&a.mesh, Some(&a.metrics)).min > 0.2);
    // complexity predicts the element count: integral of sqrt(det M) / unit area
    let ny = 4000;
    let integral: f64 = (0..ny)
        .map(|j| {
            let y = (j as f64 + 0.5) / ny as f64;
            1.0 / (0.1 * (0.004 + 0.096 * y)) / ny as f64
        })
        .sum();
    let expected = integral / unit_simplex_measure(2);
    let n = a.mesh.n_cells() as f64;
    assert!(n > 0.7 * expected && n < 1.4 * expected, "{n} vs {expected}");
    // elements near the wall are flat: their x extent is many times their y extent
    let flat = (0..a.mesh.n_cells())
        .filter(|&c| {
            let p: Vec<[f64; 2]> = a.mesh.cell(c).iter().map(|&v| a.mesh.points[v as usize]).collect();
            let ext = |k: usize| {
                p.iter().map(|q| q[k]).fold(f64::NEG_INFINITY, f64::max)
                    - p.iter().map(|q| q[k]).fold(f64::INFINITY, f64::min)
            };
            p.iter().all(|q| q[1] < 0.03) && ext(0) > 3.0 * ext(1)
        })
        .count();
    assert!(flat > 30, "flat triangles near the wall: {flat}");
}

#[test]
fn material_interface_and_labels_are_preserved() {
    let mut m = box_tri([[0.0, 2.0], [0.0, 1.0]], 6, 3);
    for c in 0..m.n_cells() {
        let cx: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][0]).sum::<f64>() / 3.0;
        m.cell_label[c] = if cx < 1.0 { 1 } else { 2 };
    }
    let f = |x: &[f64; 2]| isotropic::<2>(0.04 + 0.1 * (x[0] - 1.0).abs());
    let a = adapt2(&m, &f, &AdaptOptions::default()).unwrap();
    report("interface", &a);
    assert!((label_area(&a.mesh, 1) - 1.0).abs() < 1e-12);
    assert!((label_area(&a.mesh, 2) - 1.0).abs() < 1e-12);
    // every edge between a label-1 and a label-2 triangle lies on x = 1
    let fc = a.mesh.facet_cells();
    let mut n_if = 0;
    for (k, cs) in &fc {
        if cs.len() == 2 && a.mesh.cell_label[cs[0] as usize] != a.mesh.cell_label[cs[1] as usize] {
            n_if += 1;
            for &v in k.iter().take(2) {
                assert!((a.mesh.points[v as usize][0] - 1.0).abs() < 1e-12);
            }
        }
    }
    assert!(n_if >= 6);
    // and the whole line x = 1 is covered by interface edges
    let covered: f64 = fc
        .iter()
        .filter(|(_, cs)| cs.len() == 2 && a.mesh.cell_label[cs[0] as usize] != a.mesh.cell_label[cs[1] as usize])
        .map(|(k, _)| (a.mesh.points[k[0] as usize][1] - a.mesh.points[k[1] as usize][1]).abs())
        .sum();
    assert!((covered - 1.0).abs() < 1e-12);
}

#[test]
fn curved_boundary_converges_to_the_circle() {
    // disk of radius 1 meshed as a coarse polygon by mapping a square mesh (concentric squares)
    let n = 16;
    let mut m = SimplexMesh::<2>::default();
    let mut id = std::collections::HashMap::new();
    for j in 0..=n {
        for i in 0..=n {
            let (u, v) = (2.0 * i as f64 / n as f64 - 1.0, 2.0 * j as f64 / n as f64 - 1.0);
            let r = u.abs().max(v.abs());
            let p = if r == 0.0 {
                [0.0, 0.0]
            } else {
                let s = r / u.hypot(v);
                [u * s, v * s]
            };
            id.insert((i, j), m.points.len() as u32);
            m.points.push(p);
        }
    }
    // keep the cells in the square index space, they now fill the disk
    for j in 0..n {
        for i in 0..n {
            let q = [id[&(i, j)], id[&(i + 1, j)], id[&(i + 1, j + 1)], id[&(i, j + 1)]];
            let [t0, t1] = femcore::mesh::cells::split_quad(&q);
            m.push_cell(&t0, 1);
            m.push_cell(&t1, 1);
        }
    }
    m.orient_positive();
    let keys = m.boundary_keys();
    for k in keys {
        m.push_facet(&[k[0], k[1]], 9);
    }
    // facets must be oriented consistently only for output; check() uses keys
    m.check().unwrap();
    let area0 = m.total_measure();
    let f = iso(0.08);
    let opt = AdaptOptions::default();
    let plain = adapt2(&m, &f, &opt).unwrap();
    let curved = adapt2_with(&m, &f, &opt, |w| {
        PolylineGeometry::from_work(w).with_projector(9, circle_projector([0.0, 0.0], 1.0))
    })
    .unwrap();
    report("disk polyline", &plain);
    report("disk circle", &curved);
    let a_circle = curved.mesh.total_measure();
    println!(
        "areas: input polygon {area0:.5}, polyline-adapted {:.5}, circle-adapted {a_circle:.5}, pi {:.5}",
        plain.mesh.total_measure(),
        std::f64::consts::PI
    );
    // boundary vertices of the circle-adapted mesh lie on the circle
    let bk = curved.mesh.boundary_keys();
    for k in &bk {
        for &v in k.iter().take(2) {
            let p = curved.mesh.points[v as usize];
            assert!((p[0].hypot(p[1]) - 1.0).abs() < 1e-12, "boundary vertex off the circle: {p:?}");
        }
    }
    // the inscribed polygon underestimates the area; the circle projection closes the gap
    assert!((std::f64::consts::PI - a_circle).abs() < (std::f64::consts::PI - area0).abs());
    assert!((std::f64::consts::PI - a_circle).abs() < 0.01);
    assert!(mesh_quality(&curved.mesh, Some(&curved.metrics)).min > 0.2);
}

#[test]
fn deterministic() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5);
    let f = |x: &[f64; 2]| isotropic::<2>(0.03 + 0.1 * x[1]);
    let a = adapt2(&m, &f, &AdaptOptions::default()).unwrap();
    let b = adapt2(&m, &f, &AdaptOptions::default()).unwrap();
    assert_eq!(a.mesh.cells, b.mesh.cells);
    assert_eq!(a.mesh.cell_label, b.mesh.cell_label);
    for (p, q) in a.mesh.points.iter().zip(&b.mesh.points) {
        assert_eq!(p[0].to_bits(), q[0].to_bits());
        assert_eq!(p[1].to_bits(), q[1].to_bits());
    }
}

/// P1 interpolation error of `f` on a mesh, in the max norm over a fine sampling of every triangle.
fn interp_error_max(m: &SimplexMesh<2>, f: &dyn Fn(&[f64; 2]) -> f64) -> f64 {
    let vals: Vec<f64> = m.points.iter().map(f).collect();
    let mut e: f64 = 0.0;
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        let p: Vec<[f64; 2]> = v.iter().map(|&i| m.points[i as usize]).collect();
        for &(a, b) in &[
            (1.0 / 3.0, 1.0 / 3.0),
            (0.5, 0.5),
            (0.5, 0.0),
            (0.0, 0.5),
            (0.2, 0.2),
            (0.6, 0.2),
            (0.2, 0.6),
            (0.05, 0.9),
        ] {
            let c0 = 1.0 - a - b;
            let x = [c0 * p[0][0] + a * p[1][0] + b * p[2][0], c0 * p[0][1] + a * p[1][1] + b * p[2][1]];
            let li = c0 * vals[v[0] as usize] + a * vals[v[1] as usize] + b * vals[v[2] as usize];
            e = e.max((li - f(&x)).abs());
        }
    }
    e
}

#[test]
fn hessian_driven_adaptation_beats_uniform_meshes_of_equal_size() {
    // anisotropic Gaussian bump: narrow in y, wider in x
    let g = Gaussian::<2> { centre: [0.5, 0.5], width: [0.12, 0.03] };
    let f = |x: &[f64; 2]| g.value(x);
    // background mesh and per-vertex Hessians -> metric with a target element count
    let bg = rect_tri(
        &(0..=150).map(|i| i as f64 / 150.0).collect::<Vec<_>>(),
        &(0..=150).map(|i| i as f64 / 150.0).collect::<Vec<_>>(),
    );
    let hess: Vec<Sym<2>> = bg.points.iter().map(|p| g.hessian(p)).collect();
    let vol = bg.vertex_volumes();
    let target_elements = 1500.0;
    let opt = LpOptions {
        p: f64::INFINITY,
        target_complexity: complexity_for_elements(target_elements, 2),
        hmin: 1e-3,
        hmax: 0.5,
        ..Default::default()
    };
    let lp = lp_metrics(&hess, &vol, &opt);
    let field = BackgroundMetric2::new(&bg, &lp.metrics);
    let start = box_tri([[0.0, 1.0], [0.0, 1.0]], 6, 6);
    let a = adapt2(&start, &field, &AdaptOptions::default()).unwrap();
    report("hessian-driven", &a);
    for it in &a.report.iters {
        println!(
            "  iter {}: {} tris, +{} split -{} collapse {} flip {} move, in-band {:.3}, qmin {:.3}",
            it.iter, it.n_triangles, it.splits, it.collapses, it.flips, it.moves, it.in_band, it.q_min
        );
    }
    let n_adapted = a.mesh.n_cells();
    // uniform mesh with (about) the same number of triangles
    let k = ((n_adapted as f64 / 2.0).sqrt()).round() as usize;
    let uni = box_tri([[0.0, 1.0], [0.0, 1.0]], k, k);
    let e_adapt = interp_error_max(&a.mesh, &f);
    let e_uni = interp_error_max(&uni, &f);
    println!(
        "elements: adapted {n_adapted}, uniform {}; max interpolation error: adapted {e_adapt:.3e}, uniform {e_uni:.3e} ({:.1}x)",
        uni.n_cells(),
        e_uni / e_adapt
    );
    assert!(e_adapt * 5.0 < e_uni, "adapted {e_adapt}, uniform {e_uni}");
    assert!(n_adapted as f64 > 0.5 * target_elements && (n_adapted as f64) < 2.0 * target_elements);
    let _ = <BackgroundMetric2 as MetricField<2>>::at;
}

#[test]
fn frozen_boundary_is_untouched_and_others_adapt() {
    // the frozen side x = 0 is already meshed at the requested size (25 edges of 0.04)
    let xs: Vec<f64> = (0..=5).map(|i| i as f64 / 5.0).collect();
    let ys: Vec<f64> = (0..=25).map(|i| i as f64 / 25.0).collect();
    let m = rect_tri(&xs, &ys);
    let frozen = side_label(0, false); // x = 0
    let opt = AdaptOptions {
        topo: femcore::adapt::topo::TopoOptions { frozen_labels: vec![frozen], ..Default::default() },
        ..Default::default()
    };
    let a = adapt2(&m, &iso(0.04), &opt).unwrap();
    report("frozen x=0", &a);
    // the vertices on x = 0 are exactly the original ones, in the same positions
    let mut before: Vec<[u64; 2]> =
        m.points.iter().filter(|p| p[0] == 0.0).map(|p| [p[0].to_bits(), p[1].to_bits()]).collect();
    let mut after: Vec<[u64; 2]> =
        a.mesh.points.iter().filter(|p| p[0] == 0.0).map(|p| [p[0].to_bits(), p[1].to_bits()]).collect();
    before.sort();
    after.sort();
    assert_eq!(before, after);
    // other sides were refined
    let on_right = a.mesh.points.iter().filter(|p| p[0] == 1.0).count();
    assert!(on_right > 20, "{on_right}");
    assert_eq!(a.report.unresolved_frozen, 0);
    assert!((a.mesh.total_measure() - 1.0).abs() < 1e-12);
    assert!(mesh_quality(&a.mesh, Some(&a.metrics)).min > 0.3);
}

#[test]
fn unresolved_frozen_edges_terminate_and_are_reported() {
    // a frozen side far coarser than the metric asks for: cannot be fixed, must not run away
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5);
    let opt = AdaptOptions {
        topo: femcore::adapt::topo::TopoOptions { frozen_labels: vec![side_label(0, false)], ..Default::default() },
        ..Default::default()
    };
    let a = adapt2(&m, &iso(0.04), &opt).unwrap();
    report("frozen coarse", &a);
    assert_eq!(a.report.unresolved_frozen, 5);
    assert!(a.mesh.n_cells() < 5000, "{}", a.mesh.n_cells());
    assert!((a.mesh.total_measure() - 1.0).abs() < 1e-12);
}

#[test]
fn l_shaped_domain_keeps_its_reentrant_corner() {
    // [0,2]^2 minus [1,2]x[1,2], built from a uniform square mesh
    let sq = box_tri([[0.0, 2.0], [0.0, 2.0]], 6, 6);
    let mut m = SimplexMesh::<2> { points: sq.points.clone(), ..Default::default() };
    for c in 0..sq.n_cells() {
        let (mut cx, mut cy) = (0.0, 0.0);
        for &v in sq.cell(c) {
            cx += sq.points[v as usize][0] / 3.0;
            cy += sq.points[v as usize][1] / 3.0;
        }
        if !(cx > 1.0 && cy > 1.0) {
            m.push_cell(sq.cell(c), 1);
        }
    }
    for k in m.boundary_keys() {
        m.push_facet(&[k[0], k[1]], 5);
    }
    m.check().unwrap();
    assert!((m.total_measure() - 3.0).abs() < 1e-12);
    let a = adapt2(&m, &iso(0.06), &AdaptOptions::default()).unwrap();
    report("L-shape", &a);
    assert!((a.mesh.total_measure() - 3.0).abs() < 1e-12);
    for corner in [[1.0, 1.0], [0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 2.0], [0.0, 2.0]] {
        assert!(a.mesh.points.contains(&corner), "corner {corner:?} lost");
    }
    // no vertex inside the removed quadrant
    assert!(a.mesh.points.iter().all(|p| !(p[0] > 1.0 + 1e-12 && p[1] > 1.0 + 1e-12)));
    assert!(mesh_quality(&a.mesh, Some(&a.metrics)).min > 0.3);
}
