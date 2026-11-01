//! The wafer frame is native: x lateral, y along the propagation direction, z the surface normal. A 2D model is the x-z
//! cross-section (entities at y = 0); a 3D waveguide is that profile swept along y. No conversions anywhere.

use geom::eval::Tier;
use geom::model::*;
use geom::select::Selector;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.0)
}

/// Substrate with a rib, as a 2D cross-section: points are (x, z).
fn section() -> Model {
    let mut m = Model::new_2d();
    let sub = m.rect("sub", (-2.0, -1.0), (2.0, 0.0));
    let rib = m.rect("rib", (-0.25, 0.0), (0.25, 0.22));
    let air = m.rect("air", (-2.0, 0.0), (2.0, 1.0));
    m.overlay(&[("air", air), ("sub", sub), ("rib", rib)]);
    m
}

#[test]
fn a_2d_model_lies_in_the_x_z_plane() {
    let ev = section().eval(Tier::Complex).unwrap();
    let b = ev.bounds();
    assert_eq!((b.lo, b.hi), ([-2.0, 0.0, -1.0], [2.0, 0.0, 1.0]));
    for v in ev.vertices() {
        assert_eq!(v.p[1], 0.0);
    }
    // the top of the rib is the edge at the largest z; on an interface the normal points from the first region (air) to the second,
    // on the outer boundary it points out of the model
    let top = ev.resolve(&Selector::parse("boundary(region(rib)) & extreme(+z)").unwrap()).unwrap().set;
    assert_eq!(top.ids.len(), 1);
    let bb = ev.complex().facet_bbox(top.ids[0]);
    assert!(close(bb.lo[2], 0.22) && close(bb.hi[2], 0.22) && close(bb.lo[0], -0.25) && close(bb.hi[0], 0.25));
    let n = ev.complex().facet_normals(top.ids[0]);
    assert!(n.iter().all(|v| close(v[2], -1.0) && v[0].abs() < 1e-14 && v[1] == 0.0));
    let lid = ev.resolve(&Selector::parse("outer & extreme(+z)").unwrap()).unwrap().set;
    assert!(ev.complex().facet_normals(lid.ids[0]).iter().all(|v| close(v[2], 1.0) && v[1] == 0.0));
    assert_eq!(ev.resolve(&Selector::parse("outer & extreme(+z)").unwrap()).unwrap().set.ids.len(), 1);
    assert_eq!(ev.resolve(&Selector::parse("outer & extreme(-z)").unwrap()).unwrap().set.ids.len(), 1);
    assert_eq!(ev.resolve(&Selector::parse("plane(z = 0)").unwrap()).unwrap().set.ids.len(), 3);
    // classification takes world points
    assert_eq!(ev.classify([0.0, 0.0, 0.1]), Some(2));
    assert_eq!(ev.classify([1.0, 0.0, 0.1]), Some(0));
    assert_eq!(ev.classify([1.0, 0.0, -0.5]), Some(1));
    // the oracle projects world points onto the exact edge
    let e = top.ids[0];
    let q = ev.oracle().project([0.1, 0.0, 0.5], geom::sim::EntityRef::Edge(e)).unwrap();
    assert!(close(q[0], 0.1) && q[1] == 0.0 && close(q[2], 0.22));
}

#[test]
fn pml_sides_of_a_2d_model_are_x_and_z() {
    let mut m = section();
    let win = m.rect("window", (-1.5, -0.5), (1.5, 0.8));
    m.pml("pml", win, &[], 0.5, 4.0, 2.0);
    let ev = m.eval(Tier::Complex).unwrap();
    let names: Vec<&str> = ev.regions().iter().map(|r| r.name.as_str()).collect();
    assert!(names.contains(&"air.pml_zmax") && names.contains(&"sub.pml_xmin_zmin"), "{names:?}");
    assert!(names.iter().all(|n| !n.contains("ymin") && !n.contains("ymax")), "{names:?}");
    let p = ev.pml_regions().iter().find(|p| ev.regions()[p.region as usize].name == "sub.pml_xmin_zmin").unwrap();
    assert_eq!(p.sides, [-1, 0, -1]);
    assert_eq!((p.window.lo, p.window.hi), ([-1.5, 0.0, -0.5], [1.5, 0.0, 0.8]));
    // a y side does not exist in a 2D model
    let mut m = section();
    let win = m.rect("window", (-1.5, -0.5), (1.5, 0.8));
    m.pml("pml", win, &[Side::YMax], 0.5, 4.0, 2.0);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("ymax"));
}

/// A 3D waveguide: the (x, z) profile swept along y.
fn waveguide() -> Model {
    let mut m = Model::new_3d();
    let air = m.box_("air", (-2.0, 0.0, -1.0), (2.0, 10.0, 1.0));
    let p = m.rect("outline", (-0.25, 0.0), (0.25, 0.22));
    let wg = m.extrude("wg", p, 1.0, 9.0);
    m.overlay(&[("air", air), ("wg", wg)]);
    m
}

#[test]
fn a_3d_extrusion_sweeps_the_cross_section_along_y() {
    let ev = waveguide().eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let r = ev.regions().iter().find(|r| r.name == "wg").unwrap();
    assert!(close(r.measure, 0.5 * 0.22 * 8.0));
    assert_eq!((r.bbox.lo, r.bbox.hi), ([-0.25, 1.0, 0.0], [0.25, 9.0, 0.22]));
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    for n in ["air|wg", "air.ymin", "air.ymax", "air.zmax", "air.zmin"] {
        assert!(names.contains(&n), "{n} in {names:?}");
    }
    // the profile edges keep their roles: the top of the rib faces +z, the ends face along y
    let top = ev.resolve(&Selector::parse("from(outline, side:top)").unwrap()).unwrap().set;
    assert_eq!(top.ids.len(), 1);
    let f = top.ids[0];
    assert!(ev.facet_normals(f).iter().all(|n| close(n[2], -1.0))); // air -> wg
    let lid = ev.resolve(&Selector::parse("name(air.zmax)").unwrap()).unwrap().set;
    assert!(ev.facet_normals(lid.ids[0]).iter().all(|n| close(n[2], 1.0)));
    assert!(close(ev.faces()[f as usize].area, 0.5 * 8.0));
    let end = ev.resolve(&Selector::parse("name(air|wg) & plane(y = 9)").unwrap()).unwrap().set;
    assert_eq!(end.ids.len(), 1);
    assert!(ev.facet_normals(end.ids[0]).iter().all(|n| close(n[1].abs(), 1.0)));
    assert!(close(ev.faces()[end.ids[0] as usize].area, 0.5 * 0.22));
    assert_eq!(ev.classify([0.0, 5.0, 0.1]), Some(1));
    assert_eq!(ev.classify([0.0, 0.5, 0.1]), Some(0));
    // the stack axis is y: the layer stack's levels are y positions
    let ls = ev.layer_stack().unwrap();
    assert_eq!(ls.axis, Axis::Y);
}

#[test]
fn the_cross_section_of_a_3d_model_is_its_2d_model() {
    // slicing the swept model at any y gives the 2D section: same areas
    let ev3 = waveguide().eval(Tier::Complex).unwrap();
    let mut m = Model::new_2d();
    let air = m.rect("air", (-2.0, -1.0), (2.0, 1.0));
    let p = m.rect("rib", (-0.25, 0.0), (0.25, 0.22));
    m.overlay(&[("air", air), ("wg", p)]);
    let ev2 = m.eval(Tier::Complex).unwrap();
    let a3 = ev3.regions().iter().find(|r| r.name == "wg").unwrap().measure;
    let a2 = ev2.regions().iter().find(|r| r.name == "wg").unwrap().measure;
    assert!(close(a3, a2 * 8.0));
}

#[test]
fn rotations_are_right_handed_about_the_named_axis() {
    // +90 degrees about y turns +z into +x (and +x into -z), in 2D and in 3D alike
    let rot = |axis| TStep::Rotate { axis, deg: 90.0.into(), about: [0.0.into(), 0.0.into(), 0.0.into()] };
    let mut m = Model::new_2d();
    let r = m.rect("r", (1.0, -0.1), (2.0, 0.1)); // along +x
    let t = m.transform("t", r, vec![rot(Axis::Y)]);
    m.overlay(&[("r", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.classify([0.0, 0.0, -1.5]), Some(0));
    assert_eq!(ev.classify([1.5, 0.0, 0.0]), None);
    let mut m = Model::new_3d();
    let r = m.box_("r", (1.0, 0.0, -0.1), (2.0, 1.0, 0.1));
    let t = m.transform("t", r, vec![rot(Axis::Y)]);
    m.overlay(&[("r", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.classify([0.0, 0.5, -1.5]), Some(0));
    assert_eq!(ev.classify([1.5, 0.5, 0.0]), None);
    // a profile extruded along y and turned about y goes through the prismatic path and agrees
    let mut m = Model::new_3d();
    let p = m.rect("p", (1.0, -0.1), (2.0, 0.1));
    let e = m.extrude("e", p, 0.0, 1.0);
    let t = m.transform("t", e, vec![rot(Axis::Y)]);
    m.overlay(&[("r", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.classify([0.0, 0.5, -1.5]), Some(0));
    // about another axis the model is polyhedral; the same rotation convention applies
    let mut m = Model::new_3d();
    let r = m.box_("r", (1.0, 0.0, -0.1), (2.0, 1.0, 0.1));
    let t = m.transform("t", r, vec![rot(Axis::Z)]); // +x -> +y
    m.overlay(&[("r", t)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.classify([1.5, 0.5, 0.0]), None); // where the box was before the turn
    assert_eq!(ev.classify([-0.5, 1.5, 0.0]), Some(0));
}

#[test]
fn tapered_extrusions_run_along_y_and_name_their_caps() {
    let mut m = Model::new_3d();
    let p = m.rect("outline", (-1.0, -1.0), (1.0, 1.0));
    let e = m.extrude_tapered("taper", p, 0.0, 1.0, 10.0);
    m.overlay(&[("taper", e)]);
    let ev = m.eval(Tier::Complex).unwrap();
    ev.complex().check().unwrap();
    let top = 2.0 - 2.0 * 10f64.to_radians().tan();
    assert!(close(ev.volume(0), (4.0 + top * top + 2.0 * top) / 3.0));
    let b = ev.bounds();
    assert_eq!((b.lo[1], b.hi[1]), (0.0, 1.0)); // swept along y
    assert_eq!((b.lo[0], b.hi[0], b.lo[2], b.hi[2]), (-1.0, 1.0, -1.0, 1.0)); // the wide end is at y = 0
    let names: Vec<&str> = ev.faces().iter().map(|f| f.name.as_str()).collect();
    assert!(names.contains(&"taper.ymin") && names.contains(&"taper.ymax"), "{names:?}");
}

#[test]
fn the_model_serialises_with_the_new_field_names_and_reads_the_old_ones() {
    let m = waveguide();
    let j = serde_json::to_string(&m.list).unwrap();
    assert!(j.contains("\"from\"") && j.contains("\"to\"") && j.contains("\"Y\""), "{j}");
    let old = j.replace("\"from\"", "\"z0\"").replace("\"to\"", "\"z1\"");
    let back: FeatureList = serde_json::from_str(&old).unwrap();
    assert_eq!(back, m.list);
    // the Rust form of a model uses the plain builder
    assert!(m.list.to_rust().contains("m.extrude(\"wg\""));
}
