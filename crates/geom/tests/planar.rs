//! Polyhedral models: rotated and tapered solids and slanted cuts, through the plane-by-plane arrangement.

use geom::eval::Tier;
use geom::model::*;
use geom::select::Selector;
use geom::sim::EntityRef;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-10 * a.abs().max(b.abs())
}

fn rot(axis: Axis, deg: f64, about: [f64; 3]) -> TStep {
    TStep::Rotate { axis, deg: deg.into(), about: [about[0].into(), about[1].into(), about[2].into()] }
}

fn counts(ev: &geom::eval::Evaluated) -> (usize, usize, usize) {
    (ev.faces().len(), ev.edges3().len(), ev.vertices().len())
}

#[test]
fn slanted_cut_of_a_box() {
    let mut m = Model::new_3d();
    let b = m.box_("block", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("block", b)]);
    m.cut("slant", (1.0, 0.0, 1.0), 1.0, Keep::Below);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 0.5));
    // a triangular prism: 5 faces, 9 edges, 6 vertices (Euler: 6 - 9 + 5 = 2)
    assert_eq!(counts(&ev), (5, 9, 6));
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert!(
        names.contains(&"slant.cut") && names.contains(&"block.xmin") && names.contains(&"block.zmin"),
        "{names:?}"
    );
    let cut = ev.resolve(&Selector::parse("name(slant.cut)").unwrap()).unwrap().set.ids[0];
    assert!(close(ev.faces()[cut as usize].area, 2f64.sqrt()));
    assert_eq!(ev.resolve(&Selector::parse("tag(symmetry)").unwrap()).unwrap().set.ids, vec![cut]);
    // the slanted face is selected by its normal; the oracle projects onto it exactly
    assert_eq!(ev.resolve(&Selector::parse("normal([1, 0, 1], 1)").unwrap()).unwrap().set.ids, vec![cut]);
    let q = ev.oracle().project([1.0, 0.5, 1.0], EntityRef::Face(cut)).unwrap();
    assert!((q[0] - 0.5).abs() < 1e-12 && (q[1] - 0.5).abs() < 1e-12 && (q[2] - 0.5).abs() < 1e-12, "{q:?}");
    assert_eq!(ev.classify([0.1, 0.5, 0.1]), Some(0));
    assert_eq!(ev.classify([0.9, 0.5, 0.9]), None);
    // volume is right for any offset: the cut at x + z = 0.5 keeps a smaller prism
    let mut m2 = m.clone();
    for f in &mut m2.list.features {
        if let Op::Cut { offset, .. } = &mut f.op {
            *offset = 0.5.into();
        }
    }
    assert!(close(m2.eval(Tier::Complex).unwrap().volume(0), 0.125));
}

#[test]
fn rotated_boxes_and_booleans() {
    // a cube turned about x keeps its volume and its topology
    let mut m = Model::new_3d();
    let b = m.box_("cube", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let t = m.transform("turned", b, vec![rot(Axis::X, 30.0, [0.0, 0.5, 0.5])]);
    m.overlay(&[("cube", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 1.0));
    assert_eq!(counts(&ev), (6, 12, 8));
    assert!(ev.edges3().iter().all(|e| close(e.length, 1.0)));
    // a small turned cube inside a big box: two regions, an A|B interface of six faces
    let mut m = Model::new_3d();
    let big = m.box_("big", (0.0, 0.0, 0.0), (4.0, 4.0, 4.0));
    let c = m.box_("small", (1.5, 1.5, 1.5), (2.5, 2.5, 2.5));
    let t = m.transform("tilt", c, vec![rot(Axis::X, 45.0, [2.0, 2.0, 2.0]), rot(Axis::Y, 20.0, [2.0, 2.0, 2.0])]);
    m.overlay(&[("big", big), ("small", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 63.0) && close(ev.volume(1), 1.0));
    assert_eq!(ev.faces().iter().filter(|f| f.name == "big|small").count(), 6);
    assert_eq!(counts(&ev), (12, 24, 16));
    // a turned cube straddling a wall: the overlap is exactly half of it (it is centred on the plane x = 2)
    let mut m = Model::new_3d();
    let a = m.box_("wall", (0.0, 0.0, 0.0), (2.0, 4.0, 4.0));
    let c = m.box_("blade", (1.5, 1.5, 1.5), (2.5, 2.5, 2.5));
    let t = m.transform("tilt", c, vec![rot(Axis::X, 45.0, [2.0, 2.0, 2.0])]);
    m.overlay(&[("wall", a), ("blade", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(1), 1.0) && close(ev.volume(0), 32.0 - 0.5));
    // fractions round trip on a face of a tilted solid
    let l = ev.bounds().diagonal();
    for (i, _) in ev.faces().iter().enumerate() {
        let p = ev.position_at(EntityRef::Face(i as u32), geom::query::Fractions { u: 0.4, v: 0.6 }).unwrap();
        let g = ev.fractions_of(EntityRef::Face(i as u32), p).unwrap();
        let q = ev.position_at(EntityRef::Face(i as u32), g).unwrap();
        let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
        assert!(d <= 1e-12 * l, "face {i}: {d}");
    }
    // tessellation covers every face exactly
    let mesh = ev.tessellate(1e-3);
    let mut by = vec![0.0; ev.faces().len()];
    for (t, &f) in mesh.tris.iter().zip(&mesh.tri_face) {
        let (a, b, c) = (mesh.points[t[0] as usize], mesh.points[t[1] as usize], mesh.points[t[2] as usize]);
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        by[f as usize] += 0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
    }
    for (f, a) in ev.faces().iter().zip(by) {
        assert!((a - f.area).abs() <= 1e-9 * f.area, "{}: {a} vs {}", f.name, f.area);
    }
}

#[test]
fn tapered_and_concave_extrusions() {
    // a 2 x 2 square tapered by 10 degrees over height 1: a frustum
    let mut m = Model::new_3d();
    let p = m.rect("outline", (0.0, 0.0), (2.0, 2.0));
    let e = m.extrude_tapered("frustum", p, 0.0, 1.0, 10.0);
    m.overlay(&[("frustum", e)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let top = 2.0 - 2.0 * 10f64.to_radians().tan();
    assert!(close(ev.volume(0), (4.0 + top * top + 2.0 * top) / 3.0));
    assert_eq!(counts(&ev), (6, 12, 8));
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert!(names.contains(&"frustum.side:left") && names.contains(&"frustum.ymax") && names.contains(&"frustum.ymin"));
    assert_eq!(ev.resolve(&Selector::parse("kind(plane)").unwrap()).unwrap().set.ids.len(), 6);
    // an L-shaped profile, turned about x: eight faces; the triangulation diagonals leave no trace
    let mut m = Model::new_3d();
    let l = m.polygon_named(
        "ell",
        &[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)],
        &["a", "b", "c", "d", "e", "f"],
    );
    let ex = m.extrude("beam", l, 0.0, 2.0);
    let t = m.transform("turned", ex, vec![rot(Axis::X, 15.0, [0.0, 0.0, 0.0])]);
    m.overlay(&[("beam", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 5.0 * 2.0));
    let mut names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), 8, "{names:?}");
    assert!(names.iter().all(|n| !n.contains("diagonal")));
    // pattern of tilted bars: volumes add up
    let mut m = Model::new_3d();
    let bar = m.box_("bar", (0.0, 0.0, 0.0), (0.2, 1.0, 0.2));
    let tilted = m.transform("tilted", bar, vec![rot(Axis::X, 10.0, [0.0, 0.0, 0.0])]);
    let row = m.linear_pattern("row", tilted, (1.0, 0.0, 0.0), 4.0);
    m.overlay(&[("bars", row)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert!(close(ev.volume(0), 4.0 * 0.04));
    // faces on the same plane with the same pair of regions are one face (several components): 4 shared planes + 2 x 4 end faces
    assert_eq!(ev.faces().len(), 12);
}

#[test]
fn what_the_polyhedral_path_refuses() {
    let msg = |f: &dyn Fn(&mut Model) -> Vec<(&'static str, FeatureId)>| {
        let mut m = Model::new_3d();
        let r = f(&mut m);
        m.overlay(&r);
        m.eval(Tier::Complex).err().expect("expected a refusal").0
    };
    let e = msg(&|m| {
        let c = m.cylinder("pin", (0.0, 0.0, 0.0), 0.5, 1.0);
        let t = m.transform("tilt", c, vec![rot(Axis::X, 30.0, [0.0; 3])]);
        vec![("a", t)]
    });
    assert!(e.contains("Cylinder 'pin'") && e.contains("rotated, tapered or slanted"), "{e}");
    let e = msg(&|m| {
        let r = m.rounded_rect("rr", (0.0, 0.0), (2.0, 2.0), 0.5);
        let x = m.extrude_tapered("t", r, 0.0, 1.0, 5.0);
        vec![("a", x)]
    });
    assert!(e.contains("RoundedRect 'rr'"), "{e}");
    let e = msg(&|m| {
        let r = m.polygon("ell", &[(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (1.0, 1.0), (1.0, 3.0), (0.0, 3.0)]);
        let x = m.extrude_tapered("t", r, 0.0, 1.0, 5.0);
        vec![("a", x)]
    });
    assert!(e.contains("convex profile"), "{e}");
    let e = msg(&|m| {
        let r = m.rect("sq", (0.0, 0.0), (1.0, 1.0));
        let x = m.extrude_tapered("t", r, 0.0, 10.0, 30.0); // closes before the top
        vec![("a", x)]
    });
    assert!(e.contains("taper"), "{e}");
    // two solids that share a face within the band but are not coincident
    let e = msg(&|m| {
        let a = m.box_("a", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let b = m.box_("b", (1.0 + 5e-8, 0.0, 0.0), (2.0, 1.0, 1.0));
        let t = m.transform("tb", b, vec![rot(Axis::X, 0.0, [0.0; 3])]);
        vec![("a", a), ("b", t)]
    });
    assert!(!e.is_empty());
}

#[test]
fn pml_around_a_window_in_a_polyhedral_model() {
    // a tilted block inside an air box; the window is the box minus a margin, with absorbing layers on all sides
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, -2.0, -2.0), (3.0, 3.0, 3.0));
    let blk = m.box_("blk", (0.2, 0.2, 0.2), (0.8, 0.8, 0.8));
    let t = m.transform("tilt", blk, vec![rot(Axis::X, 30.0, [0.5, 0.5, 0.5])]);
    let win = m.box_("window", (-1.0, -1.0, -1.0), (2.0, 2.0, 2.0));
    m.overlay(&[("air", air), ("block", t)]);
    m.pml("pml", win, &[], 1.0, 4.0, 2.0);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    assert_eq!(ev.pml_regions().len(), 26);
    assert!(close(ev.regions().iter().map(|r| r.measure).sum::<f64>(), 125.0));
    let corner = ev.pml_regions().iter().find(|p| p.sides == [1, 1, 1]).unwrap();
    assert_eq!(ev.regions()[corner.region as usize].name, "air.pml_xmax_ymax_zmax");
    assert!(close(ev.volume(corner.region), 1.0));
    // the tilted block lies in the physical part: its volume and its interface with the air are untouched
    let b = ev.regions().iter().position(|r| r.name == "block").unwrap() as u32;
    assert!(close(ev.volume(b), 0.6f64.powi(3)));
}

#[test]
fn polyhedral_path_refuses_sheets_with_clear_error_instead_of_dropping_them() {
    let mut m = Model::new_3d();
    let blk = m.box_("blk", (-50e-6, -50e-6, -20e-6), (50e-6, 50e-6, 20e-6));
    let ramp = m.box_("ramp", (-10e-6, -5e-6, 0.0), (10e-6, 5e-6, 3e-6));
    let slanted = m.transform("slanted_ramp", ramp, vec![rot(Axis::Y, 25.0, [0.0, 0.0, 0.0])]);
    let pad = m.rect("pad_rect", (-20e-6, -20e-6), (20e-6, 20e-6));
    m.sheet("cpw_pad", pad, 0.0);
    m.overlay(&[("blk", blk), ("ramp", slanted)]);

    assert!(geom::planar::needs_planar(&m.list));
    let err = m.eval(Tier::Complex).err().expect("slanted solid + Op::Sheet must error").0;
    assert!(err.contains("Sheet 'cpw_pad'") && err.contains("stepped_airbridge"), "unexpected error message: {err}");
}
