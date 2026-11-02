//! The sweep axis is an implementation detail: a profile swept along y and the same profile swept along z (y and z exchanged
//! in the world frame) must answer every public query alike, up to that exchange.
//! (Internally each sweep axis has its own sweep frame, see `scene::axis_perm`; this pins that the frame never shows.)

use geom::eval::{Evaluated, Tier};
use geom::model::*;
use geom::select::Selector;
use geom::sim::EntityRef;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1.0)
}

/// y <-> z
fn sw(v: [f64; 3]) -> [f64; 3] {
    [v[0], v[2], v[1]]
}

fn close3(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| close(a[i], b[i]))
}

/// A rib (a polygon-free rectangle profile) on a slab inside a box, swept along `axis`; the box and every world coordinate are
/// written in the frame of that axis.
fn model(axis: Axis) -> Model {
    let mut m = Model::new_3d();
    let (lo, hi) = match axis {
        Axis::Y => ((-2.0, 0.0, -1.0), (2.0, 10.0, 1.0)),
        _ => ((-2.0, -1.0, 0.0), (2.0, 1.0, 10.0)),
    };
    let air = m.box_("air", lo, hi);
    let p = m.rect("outline", (-0.25, 0.0), (0.25, 0.22));
    let q = m.rect("slab", (-1.5, -0.3), (1.5, 0.0));
    let wg = m.extrude_along("wg", axis, p, 1.0, 9.0);
    let sl = m.extrude_along("sl", axis, q, 0.5, 9.5);
    m.overlay(&[("air", air), ("sl", sl), ("wg", wg)]);
    m
}

/// `air.ymin` <-> `air.zmin` etc.: the roles of a box are named after world axes
fn swap_name(n: &str) -> String {
    let mut out = n.to_string();
    for (a, b) in [("ymin", "#1"), ("ymax", "#2"), ("zmin", "#3"), ("zmax", "#4")] {
        out = out.replace(a, b);
    }
    for (a, b) in [("#1", "zmin"), ("#2", "zmax"), ("#3", "ymin"), ("#4", "ymax")] {
        out = out.replace(a, b);
    }
    out
}

fn count(ev: &Evaluated, sel: &str) -> usize {
    ev.resolve(&Selector::parse(sel).unwrap()).unwrap().set.ids.len()
}

#[test]
fn a_y_sweep_and_a_z_sweep_answer_alike() {
    let ey = model(Axis::Y).eval(Tier::Complex).unwrap();
    let ez = model(Axis::Z).eval(Tier::Complex).unwrap();
    ey.complex().check().unwrap();
    ez.complex().check().unwrap();

    // bounds
    let (by, bz) = (ey.bounds(), ez.bounds());
    assert!(close3(by.lo, sw(bz.lo)) && close3(by.hi, sw(bz.hi)), "{by:?} {bz:?}");

    // regions: names, measures, boxes
    assert_eq!(ey.regions().len(), ez.regions().len());
    for (a, b) in ey.regions().iter().zip(ez.regions()) {
        assert_eq!(a.name, b.name);
        assert!(close(a.measure, b.measure), "{}: {} {}", a.name, a.measure, b.measure);
        assert!(close3(a.bbox.lo, sw(b.bbox.lo)) && close3(a.bbox.hi, sw(b.bbox.hi)), "{}", a.name);
    }
    for r in 0..ey.regions().len() as u32 {
        assert!(close(ey.volume(r), ez.volume(r)));
    }

    // faces: the same set (names with the roles of a box exchanged), measures, boxes, normals
    assert_eq!(ey.faces().len(), ez.faces().len());
    let mut fy: Vec<_> = ey.faces().iter().enumerate().map(|(i, f)| (f.name.clone(), f.area, f.bbox, i as u32)).collect();
    let mut fz: Vec<_> = ez.faces().iter().enumerate().map(|(i, f)| (swap_name(&f.name), f.area, sw_box(f.bbox), i as u32)).collect();
    let key = |t: &(String, f64, geom::math::Aabb, u32)| (t.0.clone(), (t.2.lo.map(|x| (x * 1e6).round() as i64), t.2.hi.map(|x| (x * 1e6).round() as i64)));
    fy.sort_by_key(key);
    fz.sort_by_key(key);
    for (a, b) in fy.iter().zip(&fz) {
        assert_eq!(a.0, b.0);
        assert!(close(a.1, b.1), "{}: area {} {}", a.0, a.1, b.1);
        assert!(close3(a.2.lo, b.2.lo) && close3(a.2.hi, b.2.hi), "{}: box", a.0);
        let (na, nb) = (ey.facet_normals(a.3), ez.facet_normals(b.3));
        assert_eq!(na.len(), nb.len(), "{}", a.0);
        let probe = [0.01, 5.0, 0.3];
        assert!(close3(ey.complex().facet_project(a.3, probe), sw(ez.complex().facet_project(b.3, sw(probe)))), "{}: projection", a.0);
        assert!(na.iter().zip(&nb).all(|(p, q)| close3(*p, sw(*q))), "{}: normals", a.0);
        assert_eq!(ey.complex().facet_kind(a.3), ez.complex().facet_kind(b.3));
        assert!(close(ey.complex().facet_measure(a.3), ez.complex().facet_measure(b.3)));
    }

    // crease edges and vertices
    assert_eq!(ey.edges3().len(), ez.edges3().len());
    assert_eq!(ey.vertices().len(), ez.vertices().len());
    let sorted = |ev: &Evaluated, swap: bool| {
        let mut v: Vec<[i64; 3]> = ev
            .vertices()
            .iter()
            .map(|v| {
                let p = if swap { sw(v.p) } else { v.p };
                p.map(|x| (x * 1e6).round() as i64)
            })
            .collect();
        v.sort();
        v
    };
    assert_eq!(sorted(&ey, false), sorted(&ez, true));

    // classification at world points
    for p in [[0.0, 5.0, 0.1], [0.0, 0.5, 0.1], [0.0, 5.0, -0.15], [1.2, 5.0, -0.15], [1.8, 5.0, -0.5], [0.0, 0.2, -0.1], [0.0, 9.8, 0.1], [0.0, 5.0, 0.5]] {
        assert_eq!(ey.classify(p), ez.classify(sw(p)), "classify {p:?}");
    }

    // selectors: counts of what is picked by world planes and extremes
    for (sy, sz) in [
        ("plane(z = 0)", "plane(y = 0)"),
        ("plane(z = 0.22)", "plane(y = 0.22)"),
        ("plane(y = 9)", "plane(z = 9)"),
        ("extreme(+z)", "extreme(+y)"),
        ("extreme(-z)", "extreme(-y)"),
        ("extreme(+y)", "extreme(+z)"),
        ("extreme(+x)", "extreme(+x)"),
        ("normal(+z)", "normal(+y)"),
        ("normal(+y)", "normal(+z)"),
        ("outer", "outer"),
        ("outer & extreme(+z)", "outer & extreme(+y)"),
        ("boundary(region(wg))", "boundary(region(wg))"),
        ("in_box([-0.3, 0.9, -0.1], [0.3, 9.1, 0.3])", "in_box([-0.3, -0.1, 0.9], [0.3, 0.3, 9.1])"),
    ] {
        assert_eq!(count(&ey, sy), count(&ez, sz), "{sy} against {sz}");
    }

    // the oracle: projection onto a face and onto a crease edge
    for i in 0..ey.faces().len() as u32 {
        let name = &ey.faces()[i as usize].name;
        let j = ez.faces().iter().position(|f| swap_name(&f.name) == *name && close(f.area, ey.faces()[i as usize].area) && close3(sw(f.bbox.lo), ey.faces()[i as usize].bbox.lo) && close3(sw(f.bbox.hi), ey.faces()[i as usize].bbox.hi)).unwrap() as u32;
        let p = [0.3, 4.0, 0.4];
        let (a, b) = (ey.oracle().project(p, EntityRef::Face(i)).unwrap(), ez.oracle().project(sw(p), EntityRef::Face(j)).unwrap());
        assert!(close3(a, sw(b)), "{name}: {a:?} {b:?}");
    }

    // the label map names the same boundaries
    let (ly, lz) = (ey.label_map().unwrap(), ez.label_map().unwrap());
    assert_eq!(ly.domains.len(), lz.domains.len());
    assert_eq!(ly.boundary("air|wg").len(), lz.boundary("air|wg").len());

    // the display mesh covers the same extent
    let (ty, tz) = (ey.tessellate(1e-3), ez.tessellate(1e-3));
    assert_eq!(ty.tris.len(), tz.tris.len());
}

fn sw_box(b: geom::math::Aabb) -> geom::math::Aabb {
    geom::math::Aabb { lo: sw(b.lo), hi: sw(b.hi) }
}
