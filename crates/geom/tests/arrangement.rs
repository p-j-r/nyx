//! 2D arrangement: closed forms, refusal paths, classification against the boundary loops, transforms and patterns.

use geom::eval::Tier;
use geom::math::*;
use geom::model::*;
use std::f64::consts::PI;

fn areas(m: &Model) -> Vec<(String, f64)> {
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    ev.regions().iter().map(|r| (r.name.clone(), r.measure)).collect()
}
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}
fn lens(r1: f64, r2: f64, d: f64) -> f64 {
    let a = ((d * d + r1 * r1 - r2 * r2) / (2.0 * d * r1)).acos();
    let b = ((d * d + r2 * r2 - r1 * r1) / (2.0 * d * r2)).acos();
    r1 * r1 * (a - a.sin() * a.cos()) + r2 * r2 * (b - b.sin() * b.cos())
}

#[test]
fn overlapping_discs_give_the_lens_area() {
    let mut m = Model::new_2d();
    let a = m.circle("a", (0.0, 0.0), 1.0);
    let b = m.circle("b", (1.2, 0.3), 0.8);
    let i = m.intersect("ab", &[a, b]);
    m.overlay(&[("lens", i)]);
    let d = (1.2f64).hypot(0.3);
    assert!(close(areas(&m)[0].1, lens(1.0, 0.8, d)));
    // painter order instead: a, then b on top
    let mut m = Model::new_2d();
    let a = m.circle("a", (0.0, 0.0), 1.0);
    let b = m.circle("b", (1.2, 0.3), 0.8);
    m.overlay(&[("a", a), ("b", b)]);
    let ar = areas(&m);
    assert!(close(ar[0].1, PI - lens(1.0, 0.8, d)) && close(ar[1].1, PI * 0.64));
}

#[test]
fn ellipse_rounded_rect_and_triangle() {
    let mut m = Model::new_2d();
    let e = m.ellipse("e", (3.0, -1.0), (1.5, 0.6), 30.0);
    let r = m.rounded_rect("r", (0.0, 0.0), (4.0, 2.0), 0.5);
    let t = m.polygon("t", &[(10.0, 0.0), (11.0, 0.0), (10.5, 2.0)]);
    m.overlay(&[("e", e), ("r", r), ("t", t)]);
    let a = areas(&m);
    assert!(close(a[0].1, PI * 1.5 * 0.6) && close(a[1].1, 8.0 - (4.0 - PI) * 0.25) && close(a[2].1, 1.0));
}

#[test]
fn clockwise_polygons_and_mirrors_have_the_same_area() {
    let mut m = Model::new_2d();
    let cw = m.polygon("cw", &[(0.0, 0.0), (0.0, 1.0), (2.0, 0.0)]);
    let mir = m.push("mir", Op::Transform { input: cw, steps: vec![TStep::Mirror { axis: Axis::X, at: 5.0.into() }] });
    let rot = m.push(
        "rot",
        Op::Transform {
            input: cw,
            steps: vec![
                TStep::Rotate { axis: Axis::Y, deg: 37.0.into(), about: [3.0.into(), 0.0.into(), 4.0.into()] },
                TStep::Translate([20.0.into(), 0.0.into(), 0.0.into()]),
            ],
        },
    );
    m.overlay(&[("cw", cw), ("mir", mir), ("rot", rot)]);
    for (_, a) in areas(&m) {
        assert!(close(a, 1.0));
    }
}

#[test]
fn pattern_of_holes() {
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (10.0, 4.0));
    let h = m.circle("hole", (1.0, 2.0), 0.3);
    let hs = m.linear_pattern("holes", h, (1.0, 0.0, 0.0), 9);
    let d = m.difference("pd", plate, &[hs]);
    m.overlay(&[("plate", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(close(ev.regions()[0].measure, 40.0 - 9.0 * PI * 0.09));
    assert!(ev.edges().iter().any(|e| e.name == "hole[3].edge"));
    assert_eq!(ev.regions()[0].loops.len(), 10);
}

#[test]
fn coincident_and_t_junction_edges_conform() {
    // two rectangles share part of an edge, a third touches the middle of that edge (T junction)
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (2.0, 1.0));
    let b = m.rect("b", (2.0, 0.5), (3.0, 1.5));
    let c = m.rect("c", (0.5, 1.0), (1.5, 2.0));
    m.overlay(&[("a", a), ("b", b), ("c", c)]);
    let ev = m.eval(Tier::Complex).unwrap();
    let cx = ev.complex();
    cx.check().unwrap();
    assert!(
        close(ev.regions()[0].measure, 2.0) && close(ev.regions()[1].measure, 1.0) && close(ev.regions()[2].measure, 1.0)
    );
    assert!(close(cx.interface_length(0, 1), 0.5) && close(cx.interface_length(0, 2), 1.0));
}

fn err_of(m: &Model) -> String {
    m.eval(Tier::Complex).err().expect("expected an error").0
}

#[test]
fn grazing_and_ambiguous_input_is_refused_with_names() {
    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (2.0, 1.0));
    let c = m.circle("c2", (1.7, 0.5), 0.3); // touches x = 2 tangentially
    m.overlay(&[("plate", plate), ("c2", c)]);
    let e = err_of(&m);
    assert!(e.contains("plate.right") && e.contains("c2.edge") && e.contains("tangentially"), "{e}");

    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    let b = m.rect("b", (1.0 + 5e-8, 0.0), (2.0, 1.0)); // a gap of 1e-13: neither coincident nor apart
    m.overlay(&[("a", a), ("b", b)]);
    assert!(err_of(&m).contains("close"), "{}", err_of(&m));

    let mut m = Model::new_2d();
    let bow = m.polygon("bow", &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.5, -1.0)]);
    m.overlay(&[("bow", bow)]);
    let e = err_of(&m);
    assert!(e.contains("bow") && e.contains("crosses itself"), "{e}");

    let mut m = Model::new_2d();
    let r = m.rect("r", (0.0, 0.0), (0.0, 1.0));
    m.overlay(&[("r", r)]);
    assert!(err_of(&m).contains("width"));

    let mut m = Model::new_2d();
    m.param("w", "1/0");
    let r = m.rect("r", (0.0, 0.0), ("w", 1.0));
    m.overlay(&[("r", r)]);
    assert!(err_of(&m).contains("parameter 'w'"));
}

#[test]
fn classify_agrees_with_the_boundary_loops() {
    let g = {
        let mut m = Model::new_2d();
        let clad = m.rect("clad", (-3.0, -3.0), (3.0, 3.0));
        let ring_o = m.circle("ro", (0.0, 0.0), 2.0);
        let ring_i = m.circle("ri", (0.0, 0.0), 1.0);
        let ring = m.difference("ring", ring_o, &[ring_i]);
        let e = m.ellipse("e", (0.5, 0.2), (1.8, 0.5), 20.0);
        let t = m.polygon("t", &[(-2.5, -2.5), (-1.0, -2.5), (-2.0, -1.0)]);
        m.overlay(&[("clad", clad), ("ring", ring), ("e", e), ("t", t)]);
        m
    };
    let ev = g.eval(Tier::Complex).unwrap();
    let cx = ev.complex();
    cx.check().unwrap();
    let (mut s, mut n) = (12345u64, 0);
    for _ in 0..4000 {
        let mut r = || {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let p = [r() * 7.0 - 3.5, 0.0, r() * 7.0 - 3.5];
        if ev.sdf_bound(p, 1.0) < 1e-3 {
            continue;
        }
        n += 1;
        assert_eq!(ev.classify(p), cx.locate([p[0], p[2]]), "at {p:?}");
    }
    assert!(n > 3000);
}

#[test]
fn small_features_at_micrometre_offsets() {
    // nanometre structure 5 um away from the origin: tolerance is relative to the model size
    let mut m = Model::new_2d();
    let big = m.rect("big", (5e-6, 5e-6), (5.002e-6, 5.001e-6));
    let core = m.rect("core", (5.0005e-6, 5.0002e-6), (5.0015e-6, 5.0008e-6));
    m.overlay(&[("big", big), ("core", core)]);
    let a = areas(&m);
    assert!(close(a[1].1, 1e-9 * 6e-10) && close(a[0].1, 2e-9 * 1e-9 - 6e-19));
}

#[test]
fn vertices_and_oracle() {
    let mut m = Model::new_2d();
    let d = m.circle("d", (1.0, 1.0), 0.5);
    let r = m.rect("r", (0.0, 0.0), (3.0, 3.0));
    m.overlay(&[("r", r), ("d", d)]);
    let ev = m.eval(Tier::Complex).unwrap();
    let o = ev.oracle();
    let e = ev.edges().iter().position(|e| e.name == "d|r").unwrap() as u32;
    // the world frame of a 2D model is the x-z plane
    let q = o.project([3.0, 0.0, 1.0], geom::sim::EntityRef::Edge(e)).unwrap();
    assert!((q[0] - 1.5).abs() < 1e-14 && q[1] == 0.0 && (q[2] - 1.0).abs() < 1e-14);
    assert!(dist2([q[0], q[2]], [1.5, 1.0]) < 1e-14);
}

#[test]
fn cut_and_pml_are_data_and_geometry() {
    let mut m = Model::new_2d();
    let r = m.rect("r", (-1.0, 0.0), (1.0, 1.0));
    m.overlay(&[("r", r)]);
    m.cut("sym", (1.0, 0.0, 0.0), 0.0, Keep::Above);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(close(ev.regions()[0].measure, 1.0));
    assert!(ev.edges().iter().any(|e| e.name == "sym.cut"));
    let tagged = ev.resolve(&geom::select::Selector::Tag("symmetry".into())).unwrap();
    assert_eq!(tagged.set.ids.len(), 1);
}
