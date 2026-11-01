//! The refusal catalogue: everything the kernel declines to decide or build is an `Err` whose message names the features involved
//! and says what to change. One case per refusal; nothing here may panic.

use geom::eval::Tier;
use geom::model::*;
use geom::select::Selector;

fn refused(m: &Model) -> String {
    m.eval(Tier::Complex).err().expect("expected a refusal").0
}

fn two_d(build: impl FnOnce(&mut Model) -> Vec<(&'static str, FeatureId)>) -> String {
    let mut m = Model::new_2d();
    let r = build(&mut m);
    m.overlay(&r);
    refused(&m)
}

fn three_d(build: impl FnOnce(&mut Model) -> Vec<(&'static str, FeatureId)>) -> String {
    let mut m = Model::new_3d();
    let r = build(&mut m);
    m.overlay(&r);
    refused(&m)
}

#[test]
fn undecidable_geometry_2d() {
    // tangent contact of two shapes
    let e = two_d(|m| {
        let a = m.rect("plate", (0.0, 0.0), (2.0, 1.0));
        let b = m.circle("disc", (1.7, 0.5), 0.3);
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("plate.right") && e.contains("disc.edge") && e.contains("tangentially"), "{e}");
    // nearly, but not exactly, coincident edges
    let e = two_d(|m| {
        let a = m.rect("left", (0.0, 0.0), (1.0, 1.0));
        let b = m.rect("right", (1.0 + 5e-8, 0.0), (2.0, 1.0));
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("close") && (e.contains("left") || e.contains("right")), "{e}");
    // a corner that nearly touches an edge
    let e = two_d(|m| {
        let a = m.rect("base", (0.0, 0.0), (2.0, 1.0));
        let b = m.rect("tab", (0.5, 1.0 + 5e-8), (1.0, 2.0));
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("tab") && e.contains("base"), "{e}");
    // a profile that crosses itself
    let e = two_d(|m| {
        let a = m.polygon("bow", &[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.5, -1.0)]);
        vec![("a", a)]
    });
    assert!(e.contains("bow") && e.contains("crosses itself"), "{e}");
    // a sliver thinner than the smallest distinguishable feature
    let e = two_d(|m| {
        let a = m.rect("big", (0.0, 0.0), (10.0, 10.0));
        let b = m.rect("hair", (1.0, 1.0), (1.0 + 1e-10, 9.0));
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("hair") && (e.contains("thick") || e.contains("shorter")), "{e}");
}

#[test]
fn invalid_input_names_the_input() {
    let e = two_d(|m| {
        let a = m.rect("r", (0.0, 0.0), (0.0, 1.0));
        vec![("a", a)]
    });
    assert!(e.contains("Rect 'r'") && e.contains("width"), "{e}");
    let e = two_d(|m| {
        let a = m.circle("c", (0.0, 0.0), f64::NAN);
        vec![("a", a)]
    });
    assert!(e.contains("radius of feature 'c'"), "{e}");
    let e = two_d(|m| {
        m.param("w", "1/0");
        let a = m.rect("r", (0.0, 0.0), ("w", 1.0));
        vec![("a", a)]
    });
    assert!(e.contains("parameter 'w'"), "{e}");
    let e = two_d(|m| {
        let a = m.rect("r", (0.0, 0.0), ("nope + 1", 1.0));
        vec![("a", a)]
    });
    assert!(e.contains("unknown name 'nope'"), "{e}");
    let e = two_d(|m| {
        let a = m.circle("h", (0.0, 0.0), 0.1);
        let p = m.linear_pattern("p", a, (1.0, 0.0, 0.0), 2.5);
        vec![("a", p)]
    });
    assert!(e.contains("Pattern 'p'") && e.contains("whole number"), "{e}");
    let e = two_d(|m| {
        let a = m.rounded_rect("rr", (0.0, 0.0), (1.0, 1.0), 0.8);
        vec![("a", a)]
    });
    assert!(e.contains("RoundedRect 'rr'"), "{e}");
    // list structure
    let mut m = Model::new_2d();
    let a = m.rect("same", (0.0, 0.0), (1.0, 1.0));
    m.rect("same", (2.0, 0.0), (3.0, 1.0));
    m.overlay(&[("a", a)]);
    assert!(refused(&m).contains("two features are named 'same'"));
    let mut m = Model::new_2d();
    m.rect("bad.name", (0.0, 0.0), (1.0, 1.0));
    assert!(refused(&m).contains("reserved"));
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    m.overlay(&[("a", a)]);
    m.push("second", Op::Overlay { regions: vec![], pins: vec![], naming: NameScheme::Feature });
    assert!(refused(&m).contains("only one Overlay"));
    let m = Model::new_2d();
    assert!(refused(&m).contains("no Overlay"));
    let mut l = FeatureList::new(2);
    l.features.push(Feature {
        id: 1,
        name: "u".into(),
        op: Op::Union { ids: vec![2] },
        suppressed: false,
        origin: Origin::Code,
        tags: vec![],
    });
    assert!(l.validate().unwrap_err().0.contains("does not come before"));
}

#[test]
fn unsupported_3d_content_is_refused_with_the_remedy() {
    let e = three_d(|m| {
        let a = m.rect("flat", (0.0, 0.0), (1.0, 1.0));
        vec![("a", a)]
    });
    assert!(e.contains("must be the input of an Extrude"), "{e}");
    let mut m = Model::new_2d();
    let b = m.box_("b", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    m.overlay(&[("b", b)]);
    assert!(refused(&m).contains("3D solid"));
    // a slanted cut of a curved prism is not polyhedral
    let mut m = Model::new_3d();
    let a = m.cylinder("pin", (0.0, 0.0, 0.0), 0.5, 1.0);
    m.overlay(&[("a", a)]);
    m.cut("slant", (1.0, 0.0, 1.0), 0.5, Keep::Below);
    assert!(refused(&m).contains("Cylinder 'pin'") && refused(&m).contains("slanted"));
    // cylinders along different axes in one model
    let e = three_d(|m| {
        let a = m.cylinder_along("along_x", Axis::X, (0.0, 0.0, 0.0), 0.5, 1.0);
        let b = m.cylinder("along_z", (3.0, 0.0, 0.0), 0.5, 1.0);
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("'along_x' runs along x") && e.contains("'along_z' along z"), "{e}");
    let e = three_d(|m| {
        let a = m.cylinder("pin", (0.0, 0.0, 0.0), 0.5, 1.0);
        let t = m.push(
            "tilt",
            Op::Transform {
                input: a,
                steps: vec![TStep::Rotate {
                    axis: Axis::X,
                    deg: 10.0.into(),
                    about: [0.0.into(), 0.0.into(), 0.0.into()],
                }],
            },
        );
        vec![("a", t)]
    });
    assert!(e.contains("Cylinder 'pin'") && e.contains("rotated, tapered or slanted"), "{e}");
    // z levels that are almost, but not exactly, equal
    let e = three_d(|m| {
        let a = m.box_("low", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
        let b = m.box_("high", (0.0, 0.0, 1.0 + 3e-8), (1.0, 1.0, 2.0));
        vec![("a", a), ("b", b)]
    });
    assert!(e.contains("z levels") && e.contains("low") && e.contains("high"), "{e}");
}

#[test]
fn curved_bodies_that_touch_something_are_refused() {
    let e = three_d(|m| {
        let b = m.box_("box", (0.0, 0.0, 0.0), (2.0, 2.0, 2.0));
        let s = m.sphere("ball", (1.0, 1.0, 1.9), 0.5);
        vec![("a", b), ("b", s)]
    });
    assert!(e.contains("sphere 'ball'") && e.contains("box.zmax"), "{e}");
    let e = three_d(|m| {
        let b = m.box_("box", (0.0, 0.0, 0.0), (6.0, 6.0, 6.0));
        let s1 = m.sphere("s1", (2.0, 3.0, 3.0), 1.0);
        let s2 = m.sphere("s2", (3.0, 3.0, 3.0), 1.0);
        vec![("a", b), ("b", s1), ("c", s2)]
    });
    assert!(e.contains("'s1'") && e.contains("'s2'") && e.contains("intersect"), "{e}");
    let e = three_d(|m| {
        let b = m.box_("box", (0.0, 0.0, 0.0), (6.0, 6.0, 6.0));
        let t = m.torus("ring", (3.0, 3.0, 3.0), 1.0, 0.3);
        let s = m.sphere("ball", (3.0, 3.0, 3.0), 0.5); // the torus' bounding box encloses the ball
        vec![("a", b), ("b", t), ("c", s)]
    });
    assert!(e.contains("torus 'ring'") && e.contains("sphere 'ball'"), "{e}");
    let e = three_d(|m| {
        let s = m.sphere("alone", (0.0, 0.0, 0.0), 1.0);
        vec![("a", s)]
    });
    assert!(e.contains("prismatic host"), "{e}");
    let e = three_d(|m| {
        let t = m.torus("fat", (0.0, 0.0, 0.0), 1.0, 1.5);
        vec![("a", t)]
    });
    assert!(e.contains("Torus 'fat'"), "{e}");
}

#[test]
fn labels_selectors_and_periodic_pairs() {
    // two regions pinned to one label
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    let b = m.rect("b", (2.0, 0.0), (3.0, 1.0));
    m.overlay_defs(
        vec![
            RegionDef { name: "A".into(), input: a, kind: RegionKind::Material, pinned_id: Some(1) },
            RegionDef { name: "B".into(), input: b, kind: RegionKind::Material, pinned_id: Some(1) },
        ],
        vec![],
    );
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev.label_map().unwrap_err().0.contains("both pinned to the domain label 1"));
    // a pin that matches nothing
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    m.overlay_full(&[("A", a, RegionKind::Material, None)], &[("name(nothing.here)", 5)]);
    assert!(m.eval(Tier::Complex).unwrap().label_map().unwrap_err().0.contains("matches nothing"));
    // mixed entity kinds, unknown terms and an id inside a pin
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 1.0));
    m.overlay(&[("A", a)]);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev.resolve(&Selector::parse("region(A) | outer").unwrap()).unwrap_err().0.contains("cannot combine"));
    assert!(Selector::parse("frobnicate(a)").unwrap_err().0.contains("unknown selector term"));
    // empty selectors are warnings, not errors
    assert!(!ev.resolve(&Selector::parse("name(none)").unwrap()).unwrap().warnings.is_empty());
    // a periodic pair that is not a translation of itself
    let mut m = Model::new_2d();
    let a = m.rect("cell", (0.0, 0.0), (1.0, 2.0));
    m.overlay(&[("cell", a)]);
    m.periodic("per", Selector::Name("cell.left".into()), Selector::Name("cell.right".into()), (1.0, 0.5, 0.0));
    assert!(refused(&m).contains("Periodic 'per'"));
}
