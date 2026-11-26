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

#[test]
fn aabb_bvh2_empty_and_singleton() {
    let empty = AabbBvh2::build(vec![]);
    assert!(empty.is_empty());
    assert_bvh_matches_brute_force(&[], &empty, [0.0, 0.0], &Bbox2::new([-1.0, -1.0], [1.0, 1.0]), 0.1);

    let all_empty = vec![Bbox2::EMPTY; 4];
    let bvh_empty = AabbBvh2::build(all_empty.clone());
    assert!(bvh_empty.is_empty());
    assert_bvh_matches_brute_force(&all_empty, &bvh_empty, [0.0, 0.0], &Bbox2::new([-1.0, -1.0], [1.0, 1.0]), 0.1);

    let single = vec![Bbox2::new([1.0, 2.0], [3.0, 4.0])];
    let bvh_single = AabbBvh2::build(single.clone());
    assert!(!bvh_single.is_empty());
    for pt in [[2.0, 3.0], [0.95, 2.0], [5.0, 5.0]] {
        assert_bvh_matches_brute_force(&single, &bvh_single, pt, &Bbox2::new([0.5, 1.5], [2.5, 3.5]), 0.1);
    }
}

#[test]
fn aabb_bvh2_coincident_boxes_match_brute_force() {
    let boxes = vec![Bbox2::new([-0.5, -0.5], [0.5, 0.5]); 200];
    let bvh = AabbBvh2::build(boxes.clone());
    assert!(!bvh.is_empty());
    for pt in [[0.0, 0.0], [0.5, 0.5], [-1.0, 0.0], [1.0, 0.0]] {
        assert_bvh_matches_brute_force(&boxes, &bvh, pt, &Bbox2::new([-0.25, -0.25], [0.25, 0.25]), 1e-6);
    }
}

#[test]
fn aabb_bvh2_deep_skewed_over_64_spills_stack_safely() {
    use geom::bvh::BvhNode2;
    let depth = 80usize;
    let item_boxes: Vec<Bbox2> = (0..=depth).map(|i| Bbox2::new([-(i as f64) - 1.0, -1.0], [1.0, 1.0])).collect();
    let items: Vec<u32> = (0..=depth as u32).collect();
    // Build a deep skewed internal chain (depth = 80 > inline stack capacity 64) where every internal
    // node's bounding box contains the origin [0.0, 0.0] so all queries must descend the full stack.
    let mut nodes: Vec<BvhNode2> = Vec::with_capacity(2 * depth + 1);
    for i in 0..depth {
        nodes.push(BvhNode2 {
            bbox: Bbox2::new([-(depth as f64) - 1.0, -1.0], [1.0, 1.0]),
            left_or_start: (i + 1) as u32,
            right: (depth + 1 + i) as u32,
            count: 0,
        });
    }
    // Deepest left leaf at index `depth`
    nodes.push(BvhNode2 { bbox: item_boxes[0], left_or_start: 0, right: 0, count: 1 });
    // Right leaf siblings at indices `depth + 1 ..= 2 * depth`
    for i in 0..depth {
        nodes.push(BvhNode2 { bbox: item_boxes[i + 1], left_or_start: (i + 1) as u32, right: 0, count: 1 });
    }
    let bvh = AabbBvh2 { nodes, items, item_boxes: item_boxes.clone() };
    for pt in [[0.0, 0.0], [-40.0, 0.0], [-90.0, 0.0], [0.0, 5.0]] {
        assert_bvh_matches_brute_force(&item_boxes, &bvh, pt, &Bbox2::new([-50.0, -0.5], [0.5, 0.5]), 1e-6);
    }
}

#[test]
fn aabb_bvh2_10k_primitives_brute_force_parity() {
    let mut seed = 0x5EED_BEEF_CAFE_u64;
    let mut next_f = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / ((1u64 << 53) as f64)
    };
    let mut boxes = Vec::with_capacity(10_000);
    for i in 0..10_000 {
        if i % 97 == 0 {
            boxes.push(Bbox2::EMPTY);
            continue;
        }
        let x0 = next_f() * 200.0 - 100.0;
        let y0 = next_f() * 200.0 - 100.0;
        let w = next_f() * 8.0;
        let h = next_f() * 8.0;
        boxes.push(Bbox2::new([x0, y0], [x0 + w, y0 + h]));
    }
    let bvh = AabbBvh2::build(boxes.clone());
    for _ in 0..32 {
        let px = next_f() * 220.0 - 110.0;
        let py = next_f() * 220.0 - 110.0;
        let pad = if next_f() < 0.25 { 0.0 } else { 0.25 };
        let rect = Bbox2::new([px - 4.0, py - 3.0], [px + 4.0, py + 3.0]);
        assert_bvh_matches_brute_force(&boxes, &bvh, [px, py], &rect, pad);
    }
}
