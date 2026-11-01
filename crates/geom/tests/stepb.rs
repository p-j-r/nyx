//! Step B: normalisation, cache behaviour, queries, tessellation, parameter-sweep support (signature, stable ids, fractions, diff),
//! and the printed-code round trip.

mod common;
use geom::diff::{diff, ChangeKind};
use geom::eval::{Evaluator, Tier};
use geom::math::*;
use geom::model::*;
use geom::query::DiscretizeOpts;
use geom::sim::EntityRef;
use std::f64::consts::PI;

fn strip(w: f64, extra: bool) -> FeatureList {
    let mut m = Model::new_2d();
    m.param("w", w).param("t", 300e-9);
    let clad = m.rect("clad", ("-1.5e-6", "-2e-6"), ("1.5e-6", "t+1e-6"));
    let bx = m.rect("box", ("-1.5e-6", "-2e-6"), ("1.5e-6", 0.0));
    let core = m.rect("core", ("-w/2", 0.0), ("w/2", "t"));
    let hole = m.circle("hole", (0.9e-6, 0.0), 0.2e-6);
    let mut regs = vec![("clad", clad), ("BOX", bx), ("Si", core), ("hole", hole)];
    if extra {
        let r = m.rect("extra", ("-1.4e-6", "0.5e-6"), ("-1.2e-6", "0.7e-6"));
        regs.push(("extra", r));
    }
    m.overlay(&regs);
    m.list
}

#[test]
fn normalisation_removes_redundancy_without_changing_the_answer() {
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (2.0, 1.0));
    let b = m.circle("b", (0.5, 0.5), 0.2);
    let far = m.circle("far", (10.0, 10.0), 0.2);
    let u1 = m.union("u1", &[a, b]);
    let u2 = m.union("u2", &[u1, a]); // nested and duplicated
    let aa = m.intersect("aa", &[u2, u2]); // A cap A
    let d = m.difference("d", aa, &[far]); // far cannot touch: dropped, so it adds no primitive
    m.overlay(&[("r", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.scene().prims.len(), 2);
    assert!((ev.volume(0) - 2.0).abs() < 1e-12);
    // disjoint intersection is empty
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    let b = m.rect("b", (3.0, 0.0), (4.0, 1.0));
    let i = m.intersect("i", &[a, b]);
    m.overlay(&[("a", a), ("i", i)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.regions().len(), 2);
    assert!(ev.warnings.iter().any(|w| w.contains("'i'")), "{:?}", ev.warnings);
}

#[test]
fn a_parameter_sweep_recomputes_only_what_depends_on_it() {
    // 40 separate rectangles: moving one must rebuild one primitive and the pairs that involve it
    let mut m = Model::new_2d();
    m.param("x", 0.0);
    let mut regs = vec![];
    for i in 0..40 {
        let x0 = if i == 7 { "x".to_string() } else { format!("{}", 3 * i) };
        let id = m.rect(
            &format!("r{i}"),
            (x0.as_str(), 0.0),
            (format!("{}+2", if i == 7 { "x".to_string() } else { (3 * i).to_string() }).as_str(), 1.0),
        );
        regs.push((format!("R{i}"), id));
    }
    let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
    m.overlay(&regs);
    let mut list = m.list.clone();
    list.params[0].value = Expr::Num(21.0);
    let mut ev = Evaluator::new();
    ev.eval(&list, Tier::Complex).unwrap();
    let s0 = ev.stats();
    assert_eq!(s0.prims_built, 40);
    let list2 = list.with(&[("x", 21.5)]);
    ev.eval(&list2, Tier::Complex).unwrap();
    let s1 = ev.stats();
    assert_eq!(s1.prims_built - s0.prims_built, 1, "{s1:?}");
    assert!(s1.pairs_cached > s0.pairs_cached);
    let again = ev.stats();
    ev.eval(&list2, Tier::Complex).unwrap();
    assert_eq!(ev.stats().prims_built, again.prims_built);
    assert_eq!(ev.stats().pairs_computed, again.pairs_computed);
}

use geom::expr::Expr;

#[test]
fn deviation_driven_discretisation_hits_its_target() {
    let mut m = Model::new_2d();
    let c = m.circle("c", (0.0, 0.0), 1.0);
    let e = m.ellipse("e", (5.0, 0.0), (2.0, 0.5), 25.0);
    m.overlay(&[("c", c), ("e", e)]);
    let ev = m.eval(Tier::Complex).unwrap();
    for (i, ed) in ev.edges().iter().enumerate() {
        for d in [1e-1, 1e-2, 1e-4] {
            let pts = ev.discretize_edge(i as u32, DiscretizeOpts { max_dev: d, max_len: 1e9 });
            assert_eq!(pts.first().unwrap().0, ed.t0);
            assert_eq!(pts.last().unwrap().0, ed.t1);
            let worst = pts.windows(2).map(|w| ev.chord_deviation(i as u32, w[0].0, w[1].0)).fold(0.0, f64::max);
            assert!(worst <= d * 1.0000001, "{} d={d}: {worst}", ed.name);
            // the straight facets really are within d of the exact edge
            let poly: Vec<V3> = pts.iter().map(|p| p.1).collect();
            assert!(ev.facet_deviation(i as u32, &poly) <= d * 1.0000001 + 1e-15);
            // and the sizes are not wasteful: the count is within 2x of the minimum for circles
            if ed.name == "c.edge" {
                let n = (PI / (1.0 - d).acos()).ceil() as usize;
                assert!(pts.len() - 1 <= 2 * n.max(4), "{} chords for d={d}", pts.len() - 1);
            }
        }
        // size limit
        let pts = ev.discretize_edge(i as u32, DiscretizeOpts { max_dev: 1.0, max_len: 0.3 });
        assert!(pts.windows(2).all(|w| dist2([w[0].1[0], w[0].1[2]], [w[1].1[0], w[1].1[2]]) <= 0.3 + 1e-12));
    }
    assert!((ev.min_radius_of_curvature(0) - 1.0).abs() < 1e-14);
    let ell = ev.edges().iter().position(|e| e.name == "e.edge").unwrap() as u32;
    assert!((ev.min_radius_of_curvature(ell) - 0.25 / 2.0).abs() < 1e-14);
}

#[test]
fn midpoints_are_on_the_curve_and_cross_the_seam() {
    let mut m = Model::new_2d();
    let c = m.circle("c", (1.0, 2.0), 3.0);
    m.overlay(&[("c", c)]);
    let ev = m.eval(Tier::Complex).unwrap();
    let ed = &ev.edges()[0];
    // the chord across the seam t = 2 pi - 0.2 .. 2 pi + 0.2 gives the point at angle 0
    let (t, p) = ev.point_at_mid(0, ed.t1 - 0.2, ed.t1 + 0.2);
    assert!((t - ed.t1).abs() < 1e-12 && dist2([p[0], p[2]], [4.0, 2.0]) < 1e-12);
    let a = ev.edge_eval(0, ed.t1 - 0.2);
    let b = ev.edge_eval(0, ed.t1 + 0.2);
    let q = ev.project_midpoint(0, a, b);
    assert!(dist2([q[0], q[2]], [4.0, 2.0]) < 1e-12);
    // curvature and tangent
    assert!((ev.edge_curvature(0, 1.0) - 1.0 / 3.0).abs() < 1e-14);
    let tg = ev.edge_tangent(0, 0.0);
    assert!((tg[2] - 1.0).abs() < 1e-14);
}

#[test]
fn tessellation_covers_the_regions_and_carries_ids() {
    let g = common::plate_with_hole();
    let ev = g.model.eval(Tier::Display).unwrap();
    let mesh = ev.tessellate(1e-4);
    let area: f64 = mesh
        .tris
        .iter()
        .map(|t| {
            0.5 * ((mesh.points[t[1] as usize][0] - mesh.points[t[0] as usize][0])
                * (mesh.points[t[2] as usize][2] - mesh.points[t[0] as usize][2])
                - (mesh.points[t[1] as usize][2] - mesh.points[t[0] as usize][2])
                    * (mesh.points[t[2] as usize][0] - mesh.points[t[0] as usize][0]))
        })
        .sum();
    assert!((area - ev.volume(0)).abs() < 1e-3 * ev.volume(0), "{area} vs {}", ev.volume(0));
    assert!(mesh.tri_face.iter().all(|&f| f == 0));
    assert_eq!(mesh.segs.len(), mesh.seg_edge.len());
    // every chord on the circle edge is within max_dev of the circle
    for (s, &e) in mesh.segs.iter().zip(&mesh.seg_edge) {
        let (a, b) = (mesh.points[s[0] as usize], mesh.points[s[1] as usize]);
        let m = [0.5 * (a[0] + b[0]), 0.0, 0.5 * (a[2] + b[2])];
        let q = ev.edge_closest(e, m).0;
        assert!(dist2([m[0], m[2]], [q[0], q[2]]) <= 1.0001e-4);
    }
    // a three-region model with nested holes triangulates to the right areas
    let g = common::concentric_discs();
    let ev = g.model.eval(Tier::Display).unwrap();
    let mesh = ev.tessellate(1e-5);
    let mut by = [0.0; 3];
    for (t, &f) in mesh.tris.iter().zip(&mesh.tri_face) {
        let p = |i: u32| mesh.points[i as usize];
        by[f as usize] += 0.5
            * ((p(t[1])[0] - p(t[0])[0]) * (p(t[2])[2] - p(t[0])[2])
                - (p(t[1])[2] - p(t[0])[2]) * (p(t[2])[0] - p(t[0])[0]));
    }
    for (r, a) in ev.regions().iter().zip(by) {
        assert!((a - r.measure).abs() < 1e-4 * r.measure, "{}: {a} vs {}", r.name, r.measure);
    }
}

#[test]
fn signature_ids_and_names_across_a_sweep() {
    let mut e = Evaluator::new();
    let base = e.eval(&strip(500e-9, false), Tier::Complex).unwrap();
    let sig = base.topology_signature();
    for w in [300e-9, 450e-9, 700e-9, 1.2e-6, 1.3e-6] {
        let ev = e.eval(&strip(w, false), Tier::Complex).unwrap();
        assert_eq!(ev.topology_signature(), sig, "w = {w}");
        assert_eq!(ev.regions().len(), base.regions().len());
        let names = |x: &geom::eval::Evaluated| x.edges().iter().map(|e| e.name.clone()).collect::<Vec<_>>();
        assert_eq!(names(&ev), names(&base));
        for (a, b) in ev.edges().iter().zip(base.edges()) {
            assert_eq!((a.v0, a.v1, a.left, a.right, &a.prov), (b.v0, b.v1, b.left, b.right, &b.prov));
        }
        assert_eq!(
            ev.vertices().iter().map(|v| v.name.clone()).collect::<Vec<_>>(),
            base.vertices().iter().map(|v| v.name.clone()).collect::<Vec<_>>()
        );
        assert_eq!(diff(&base, &ev).kind, if w == 500e-9 { ChangeKind::None } else { ChangeKind::Morph });
    }
    // a width so large that the core reaches the hole: the topology changes
    let wide = e.eval(&strip(1.6e-6, false), Tier::Complex);
    // the core edges now cross the hole circle: the structure changes
    let w = wide.unwrap();
    assert_ne!(w.topology_signature(), sig);
    assert_eq!(diff(&base, &w).kind, ChangeKind::Remesh);
    // an unrelated feature changes the signature but not the names of the existing edges
    let more = e.eval(&strip(500e-9, true), Tier::Complex).unwrap();
    assert_ne!(more.topology_signature(), sig);
    for ed in base.edges() {
        assert!(more.edges().iter().any(|x| x.name == ed.name), "{}", ed.name);
    }
    let d = diff(&base, &more);
    assert_eq!(d.kind, ChangeKind::Remesh);
    assert!(d.changed_features.len() >= 2); // the new rectangle and the overlay
                                            // removing a feature (suppressing the hole) is a Remesh as well
    let mut l = strip(500e-9, false);
    l.features.iter_mut().find(|f| f.name == "hole").unwrap().suppressed = true;
    let nh = e.eval(&l, Tier::Complex).unwrap();
    assert_eq!(diff(&base, &nh).kind, ChangeKind::Remesh);
}

#[test]
fn fractions_round_trip_and_move_continuously() {
    let mut e = Evaluator::new();
    let a = e.eval(&strip(500e-9, false), Tier::Complex).unwrap();
    let l = a.bounds().diagonal();
    for (i, ed) in a.edges().iter().enumerate() {
        for k in 0..=10 {
            let f = geom::query::Fractions { u: k as f64 / 10.0, v: 0.0 };
            let p = a.position_at(EntityRef::Edge(i as u32), f).unwrap();
            let g = a.fractions_of(EntityRef::Edge(i as u32), p).unwrap();
            let q = a.position_at(EntityRef::Edge(i as u32), g).unwrap();
            assert!(dist2([p[0], p[1]], [q[0], q[1]]) <= 1e-12 * l, "{} at {}", ed.name, f.u);
        }
    }
    // continuity: move the width a little; points keep their place to the order of the step (including the circle seam)
    let step = 1e-9;
    let b = e.eval(&strip(500e-9 + step, false), Tier::Complex).unwrap();
    for i in 0..a.edges().len() {
        for k in 0..=20 {
            let f = geom::query::Fractions { u: k as f64 / 20.0, v: 0.0 };
            let p = a.position_at(EntityRef::Edge(i as u32), f).unwrap();
            let g = b.fractions_of(EntityRef::Edge(i as u32), p).unwrap();
            let q = b.position_at(EntityRef::Edge(i as u32), g).unwrap();
            assert!(dist2([p[0], p[1]], [q[0], q[1]]) <= 2.0 * step + 1e-12 * l, "edge {i} f={}", f.u);
            let q2 = b.position_at(EntityRef::Edge(i as u32), f).unwrap();
            assert!(dist2([p[0], p[1]], [q2[0], q2[1]]) <= 2.0 * step + 1e-12 * l, "edge {i} same fraction jumped");
        }
    }
    // a circle whose radius changes: same fraction stays near the same angle (no seam jump)
    let mk = |r: f64| {
        let mut m = Model::new_2d();
        let c = m.circle("c", (0.0, 0.0), r);
        m.overlay(&[("c", c)]);
        m.eval(Tier::Complex).unwrap()
    };
    let (c1, c2) = (mk(1.0), mk(1.001));
    for k in 0..=40 {
        let f = geom::query::Fractions { u: k as f64 / 40.0, v: 0.0 };
        let (p, q) = (c1.position_at(EntityRef::Edge(0), f).unwrap(), c2.position_at(EntityRef::Edge(0), f).unwrap());
        assert!(dist2([p[0], p[1]], [q[0], q[1]]) <= 1.1e-3);
    }
}

#[test]
fn to_rust_round_trips_as_data() {
    for g in common::all_models() {
        let code = g.model.list.to_rust();
        let back = common::interp::build(&code);
        assert_eq!(back.list, g.model.list, "{}\n{code}", g.name);
        assert_eq!(code, back.list.to_rust());
    }
    let mut m = Model::new_2d();
    m.param("n", 3.0);
    let r = m.rect("a b", ("-1e-6", "0"), (1.0, "n*2"));
    let ps = m.polygon_named(
        "tri",
        &[(Expr::from(0.0), Expr::from(0.0)), (Expr::from(1.0), Expr::from(0.0)), (Expr::from(0.0), Expr::from("n"))],
        &["x", "y", "z"],
    );
    let h = m.circle("h", (0.0, 0.0), 0.1);
    let hp = m.linear_pattern("hp", h, (0.5, 0.0, 0.0), 3.0);
    let d = m.difference("d", r, &[hp]);
    m.overlay_full(
        &[("one", d, RegionKind::Material, Some(5)), ("tri", ps, RegionKind::Hole, None)],
        &[("name(a b.left)", 9)],
    );
    m.suppress(h, false);
    m.tag_feature(r, "frame");
    assert_eq!(common::interp::build(&m.list.to_rust()).list, m.list);
}

#[test]
fn evaluation_is_deterministic() {
    for g in common::all_2d() {
        let (a, b) = (g.model.eval(Tier::Complex).unwrap(), g.model.eval(Tier::Complex).unwrap());
        assert_eq!(a.topology_signature(), b.topology_signature());
        for (x, y) in a.regions().iter().zip(b.regions()) {
            assert_eq!(x.measure.to_bits(), y.measure.to_bits());
        }
        for (x, y) in a.edges().iter().zip(b.edges()) {
            assert_eq!((x.t0.to_bits(), x.t1.to_bits(), &x.name), (y.t0.to_bits(), y.t1.to_bits(), &y.name));
        }
    }
}
