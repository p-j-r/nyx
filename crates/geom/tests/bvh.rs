// SPDX-License-Identifier: BSD-2-Clause
// Copyright (C) 2026 p-j-r
//! Property-based tests (`proptest`) verifying `AabbBvh2` and `ComplexBvh` against brute-force linear scans.

use geom::bvh::{AabbBvh2, Bbox2, ComplexBvh};
use geom::eval::Tier;
use geom::expr::Expr;
use geom::math::V2;
use geom::model::Model;
use proptest::prelude::*;

fn arb_coord() -> impl Strategy<Value = f64> {
    (-50.0_f64..50.0).prop_map(|x| (x * 4.0).round() * 0.25)
}

fn arb_bbox2() -> impl Strategy<Value = Bbox2> {
    prop_oneof![
        1 => Just(Bbox2::EMPTY),
        6 => (arb_coord(), arb_coord(), 0.0_f64..20.0, 0.0_f64..20.0).prop_map(|(x, y, dx, dy)| {
            Bbox2::new([x, y], [x + (dx * 4.0).round() * 0.25, y + (dy * 4.0).round() * 0.25])
        }),
    ]
}

fn arb_query() -> impl Strategy<Value = (V2, Bbox2, f64)> {
    ([arb_coord(), arb_coord()], arb_bbox2(), prop_oneof![Just(0.0_f64), 0.0_f64..1.0])
}

fn collect_bvh_sorted(f: impl FnOnce(&mut dyn FnMut(u32))) -> Vec<u32> {
    let mut out = Vec::new();
    let mut push = |id: u32| out.push(id);
    f(&mut push);
    let raw_len = out.len();
    out.sort_unstable();
    out.dedup();
    assert_eq!(out.len(), raw_len, "BVH query visited an item more than once");
    out
}

fn assert_bvh_matches_brute_force(boxes: &[Bbox2], bvh: &AabbBvh2, p: V2, rect: &Bbox2, pad: f64) {
    let expected_pt: Vec<u32> = (0..boxes.len() as u32)
        .filter(|&i| !boxes[i as usize].is_empty() && boxes[i as usize].contains_point(p, pad))
        .collect();
    let got_pt = collect_bvh_sorted(|cb| bvh.query_point(p, pad, cb));
    assert_eq!(got_pt, expected_pt, "query_point mismatch for p={p:?} pad={pad}");

    let expected_box: Vec<u32> = (0..boxes.len() as u32)
        .filter(|&i| !boxes[i as usize].is_empty() && !rect.is_empty() && boxes[i as usize].overlaps(rect, pad))
        .collect();
    let got_box = collect_bvh_sorted(|cb| bvh.query_box(rect, pad, cb));
    assert_eq!(got_box, expected_box, "query_box mismatch for rect={rect:?} pad={pad}");

    let expected_ray: Vec<u32> = (0..boxes.len() as u32)
        .filter(|&i| !boxes[i as usize].is_empty() && boxes[i as usize].overlaps_ray_x(p, pad))
        .collect();
    let got_ray = collect_bvh_sorted(|cb| bvh.query_ray_x(p, pad, cb));
    assert_eq!(got_ray, expected_ray, "query_ray_x mismatch for p={p:?} pad={pad}");
}

proptest! {
    #[test]
    fn aabb_bvh2_matches_brute_force_queries(
        boxes in prop::collection::vec(arb_bbox2(), 0..96),
        queries in prop::collection::vec(arb_query(), 1..24),
    ) {
        let bvh = AabbBvh2::build(boxes.clone());
        let expected_empty = boxes.iter().all(Bbox2::is_empty);
        prop_assert_eq!(bvh.is_empty(), expected_empty);

        for (p, rect, pad) in queries {
            assert_bvh_matches_brute_force(&boxes, &bvh, p, &rect, pad);
        }
    }

    #[test]
    fn complex_bvh_locate_and_loop_edge_bvhs_match_brute_force(
        n_ngon in 32usize..44,
        cx_off in -2.0_f64..2.0,
        cy_off in -2.0_f64..2.0,
        r in 0.4_f64..1.5,
        samples in prop::collection::vec(((-4.0_f64..4.0), (-4.0_f64..4.0)), 1..24),
    ) {
        let mut m = Model::new_2d();
        let outer = m.rect("outer", (-5.0, -5.0), (5.0, 5.0));
        let pts: Vec<[Expr; 2]> = (0..n_ngon)
            .map(|k| {
                let th = std::f64::consts::TAU * (k as f64) / (n_ngon as f64);
                [Expr::from(cx_off + r * th.cos()), Expr::from(cy_off + r * th.sin())]
            })
            .collect();
        let hole_poly = m.polygons("ngon", vec![pts], vec![]);
        let frame = m.difference("frame", outer, &[hole_poly]);
        m.overlay(&[("frame", frame), ("hole", hole_poly)]);

        let ev = m.eval(Tier::Complex).expect("valid 2D model with >=32-edge polygon");
        let cx = ev.complex();
        cx.check().expect("complex valid");

        // Loop with >=32 edges must populate its dedicated loop_edge_bvh.
        prop_assert!(
            cx.bvh.loop_edge_bvhs.iter().any(Option::is_some),
            "expected at least one loop_edge_bvh for {n_ngon}-gon"
        );

        let mut brute_cx = cx.clone();
        brute_cx.bvh = ComplexBvh::default();

        let edge_boxes = &cx.bvh.edges.item_boxes;
        let loop_boxes = &cx.bvh.loops.item_boxes;
        let region_boxes = &cx.bvh.regions.item_boxes;
        let pad = cx.tol.max(1e-12);

        for (ux, uy) in samples {
            let p = [ux, uy];
            // Locate must agree with unaccelerated brute-force linear region scan.
            prop_assert_eq!(cx.locate(p), brute_cx.locate(p));

            let query_rect = Bbox2::new([ux - 0.5, uy - 0.5], [ux + 0.5, uy + 0.5]);
            assert_bvh_matches_brute_force(region_boxes, &cx.bvh.regions, p, &query_rect, pad);
            assert_bvh_matches_brute_force(loop_boxes, &cx.bvh.loops, p, &query_rect, pad);
            assert_bvh_matches_brute_force(edge_boxes, &cx.bvh.edges, p, &query_rect, pad);

            for (li, maybe_bvh) in cx.bvh.loop_edge_bvhs.iter().enumerate() {
                if let Some(lbvh) = maybe_bvh {
                    let lboxes = &lbvh.item_boxes;
                    prop_assert_eq!(lboxes.len(), cx.loops[li].edges.len());
                    assert_bvh_matches_brute_force(lboxes, lbvh, p, &query_rect, pad);
                }
            }
        }
    }
}
