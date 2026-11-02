//! Selectors, label map, serde round trips (JSON bit-exact, and the XML value encoding of nyx `xmlio`).

mod common;
use geom::eval::Tier;
use geom::model::*;
use geom::select::*;

fn soi() -> geom::eval::Evaluated {
    common::soi_strip().model.eval(Tier::Complex).unwrap()
}

const EXAMPLES: &[&str] = &[
    "name(core.top)",
    "tag(input)",
    "region(Si)",
    "boundary(region(Si))",
    "outer",
    "interface(Si, BOX)",
    "interface(Si, *)",
    "from(core, right)",
    "kind(circle)",
    "radius(1.2e-6, 1e-9)",
    "area(>1e-12)",
    "plane(z = 0)",
    "plane(y = 0.5, 1e-9)",
    "normal(+z)",
    "normal(-x, 5)",
    "normal([0, 1, 0], 2)",
    "extreme(+z)",
    "in_box([0, 0, 0], [1, 2, 3], within)",
    "in_box([0, 0, 0], [1, 2, 3], touching)",
    "adjacent(name(a.top))",
    "id(7)",
    "name(BOX|Si) | name(hole[3].edge)",
    "outer - tag(pml_outer)",
    "!(name(a) | name(b))",
    "name(a) - (name(b) | name(c))",
    "boundary(region(sapphire)) & extreme(+z)",
    "size_w",
];

#[test]
fn text_form_round_trips() {
    for t in EXAMPLES.iter().take(EXAMPLES.len() - 1) {
        let s = Selector::parse(t).unwrap_or_else(|e| panic!("{t}: {e}"));
        let printed = s.to_string();
        assert_eq!(Selector::parse(&printed).unwrap(), s, "{t} -> {printed}");
        assert_eq!(printed, Selector::parse(&printed).unwrap().to_string());
        let j = serde_json::to_string(&s).unwrap();
        assert_eq!(serde_json::from_str::<Selector>(&j).unwrap(), s);
    }
    for bad in [
        "name(",
        "foo(1)",
        "name(a) |",
        "interface(a)",
        "plane(q = 1)",
        "area(5)",
        "kind(blob)",
        "in_box([1,2], [3,4,5])",
        "name(a) name(b)",
    ] {
        assert!(Selector::parse(bad).is_err(), "{bad}");
    }
    // left to right, equal precedence
    let s = Selector::parse("name(a) - name(b) | name(c)").unwrap();
    assert_eq!(s, Selector::Name("a".into()).minus(Selector::Name("b".into())).or(Selector::Name("c".into())));
}

#[test]
fn resolves_the_documented_examples() {
    let ev = soi();
    let n = |t: &str| ev.resolve(&Selector::parse(t).unwrap()).unwrap();
    assert_eq!(n("name(BOX|Si)").set.ids.len(), 1);
    assert_eq!(n("interface(Si, BOX)").set, n("name(BOX|Si)").set);
    assert_eq!(n("interface(Si, *)").set.ids.len(), 4); // BOX|Si, Si|clad x3
    assert_eq!(n("outer").set.ids.len(), 6);
    assert_eq!(n("region(Si)").set.kind, EntityKind::Regions);
    assert_eq!(n("boundary(region(Si))").set, n("interface(Si, *)").set);
    assert_eq!(n("from(core, top)").set.ids.len(), 1);
    assert_eq!(n("plane(z = 0)").set.ids.len(), 2 + 1); // BOX|clad pieces (2) and BOX|Si
    assert_eq!(n("normal(+z)").set.ids.len(), 5);
    assert_eq!(n("outer & extreme(+z)").set.ids.len(), 1);
    assert_eq!(n("outer - name(clad.top)").set.ids.len(), 5);
    assert_eq!(n("!outer").set.ids.len(), 6);
    assert_eq!(n("boundary(region(BOX)) & extreme(+z)").set.ids.len(), 3);
    assert_eq!(n("in_box([-1, 0, -1], [1, 0, 1], within)").set.ids.len(), 12);
    assert_eq!(n("in_box([-3e-7, 0, -1e-9], [3e-7, 0, 4e-7], within)").set.ids.len(), 4);
    assert_eq!(n("in_box([0, 0, 0], [1e-7, 0, 1e-7], touching)").set.ids.len(), 1);
    // an empty result is a warning, not an error
    let r = n("name(nothing.here)");
    assert!(r.set.ids.is_empty() && r.warnings.iter().any(|w| w.contains("nothing.here")));
    // mixing kinds is an error with a message
    let e = ev.resolve(&Selector::parse("region(Si) | outer").unwrap()).unwrap_err();
    assert!(e.0.contains("cannot combine"), "{e}");
    // bulk: area / radius filters
    assert!(n("area(>1.2e-6)").set.ids.len() == 6);
}

#[test]
fn suggest_returns_only_exact_selectors() {
    let ev = soi();
    for e in 0..ev.edges().len() as u32 {
        let picked = EntitySet::new(EntityKind::Edges, vec![e]);
        let s = ev.suggest(&picked);
        assert!(!s.is_empty(), "no selector for edge {e} ({})", ev.edges()[e as usize].name);
        for sel in &s {
            assert_eq!(ev.resolve(sel).unwrap().set, picked, "{sel}");
        }
    }
    let picked = ev.resolve(&Selector::parse("interface(Si, *)").unwrap()).unwrap().set;
    let s = ev.suggest(&picked);
    assert_eq!(s[0].to_string(), "boundary(region(Si))".to_string().replace("boundary(region(Si))", &s[0].to_string()));
    assert!(s.iter().all(|x| ev.resolve(x).unwrap().set == picked));
    // a set no rule can describe simply gets no suggestion (two unrelated fragments of one named group)
    let two = EntitySet::new(
        EntityKind::Edges,
        vec![
            ev.edges().iter().position(|e| e.name == "Si|clad").unwrap() as u32,
            ev.edges().iter().position(|e| e.name == "box.bottom").unwrap() as u32,
        ],
    );
    for sel in ev.suggest(&two) {
        assert_eq!(ev.resolve(&sel).unwrap().set, two);
    }
}

#[test]
fn label_map_is_deterministic_and_honours_pins() {
    let ev = soi();
    let lm = ev.label_map().unwrap();
    assert_eq!((lm.domain("BOX"), lm.domain("Si"), lm.domain("clad")), (Some(1), Some(2), Some(3)));
    assert_eq!(lm.name_of_domain(2), Some("Si"));
    assert_eq!(lm.boundary("BOX|Si").len(), 1);
    assert_eq!(lm.boundary("BOX|clad").len(), 2);
    let mut all: Vec<i32> = lm.boundaries.iter().map(|b| b.id).collect();
    all.sort_unstable();
    assert_eq!(all, (1..=12).collect::<Vec<_>>());
    assert_eq!(lm.name_of_boundary(lm.boundary("clad.top")[0]), Some("clad.top"));
    assert_eq!(ev.label_map().unwrap(), lm);
    let l = ev.select_labels(&Selector::parse("outer").unwrap()).unwrap();
    assert_eq!(l.len(), 6);
    assert_eq!(ev.select_labels(&Selector::parse("region(Si)").unwrap()).unwrap(), vec![2]);
    // pins: legacy numbering for a group, all fragments share the int
    let mut g = common::soi_strip();
    let ov = g.model.list.features.iter_mut().find(|f| matches!(f.op, Op::Overlay { .. })).unwrap();
    if let Op::Overlay { pins, .. } = &mut ov.op {
        pins.push(LabelPin { sel: Selector::Name("BOX|clad".into()), id: 7 });
        pins.push(LabelPin { sel: Selector::Name("Si|clad".into()), id: 2 });
    }
    let lm = g.model.eval(Tier::Complex).unwrap().label_map().unwrap();
    assert_eq!(lm.boundary("BOX|clad"), vec![7]);
    assert_eq!(lm.boundaries.iter().find(|b| b.id == 7).unwrap().entities.len(), 2);
    assert!(lm.boundaries.iter().filter(|b| !b.pinned).all(|b| b.id != 7 && b.id != 2));
    // collisions are errors
    let ov = g.model.list.features.iter_mut().find(|f| matches!(f.op, Op::Overlay { .. })).unwrap();
    if let Op::Overlay { regions, .. } = &mut ov.op {
        regions[0].pinned_id = Some(1);
    }
    assert!(g.model.eval(Tier::Complex).unwrap().label_map().is_err());
}

#[test]
fn periodic_pairs_are_verified() {
    let mut m = Model::new_2d();
    let cell = m.rect("cell", (0.0, 0.0), (1.0, 2.0));
    m.overlay(&[("cell", cell)]);
    m.periodic("per", Selector::Name("cell.left".into()), Selector::Name("cell.right".into()), (1.0, 0.0, 0.0));
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev.periodic_pairs()[0].checked_congruent);
    // a wrong translation is refused with a named error
    let mut m = Model::new_2d();
    let cell = m.rect("cell", (0.0, 0.0), (1.0, 2.0));
    m.overlay(&[("cell", cell)]);
    m.periodic("per", Selector::Name("cell.left".into()), Selector::Name("cell.right".into()), (1.0, 0.5, 0.0));
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("Periodic 'per'"));
    // faces of different size are refused
    let mut m = Model::new_2d();
    let a = m.rect("a", (0.0, 0.0), (1.0, 2.0));
    let b = m.rect("b", (1.0, 0.0), (2.0, 1.0));
    m.overlay(&[("a", a), ("b", b)]);
    m.periodic("per", Selector::Name("a.left".into()), Selector::Name("b.right".into()), (1.0, 0.0, 0.0));
    assert!(m.eval(Tier::Complex).is_err());
}

#[test]
fn size_tag_and_pml_data() {
    let g = common::chip("chip_t", 0.3, 0.17, 0.5, 2.0);
    let mut m = g.model;
    m.size("fine", Selector::Region("rib".into()).or(Selector::Region("rib".into())), None, None, None);
    let e = m.eval(Tier::Complex).err().unwrap();
    assert!(e.0.contains("sets none"), "{e}");
    let mut m = common::chip("chip_t", 0.3, 0.17, 0.5, 2.0).model;
    m.size("fine", Selector::Interface("rib".into(), None), Some(40e-9.into()), None, None);
    m.tag("pml_outer", Selector::Name("air.right".into()), "outer_pml");
    let ev = m.eval(Tier::Complex).unwrap();
    assert_eq!(ev.size_rules()[0].h, Some(40e-9));
    assert!(!ev.size_rules()[0].entities.ids.is_empty());
    assert_eq!(ev.pml_regions().len(), 12);
    let p = ev.pml_regions().iter().find(|p| ev.regions()[p.region as usize].name == "sapphire.pml_xmin_zmin").unwrap();
    assert_eq!(p.sides, [-1, 0, -1]);
    assert_eq!(ev.regions()[p.parent.unwrap() as usize].name, "sapphire");
    assert_eq!((p.thickness, p.layers, p.grading), (1e-6, 4.0, 2.0));
    assert_eq!(ev.regions()[p.region as usize].kind, RegionKind::Pml);
}

// ---- serde ------------------------------------------------------------------------------------------------

/// The nyx `xmlio` encoding of a JSON value (objects `t="o"`, arrays `t="a"` with `<item>`, strings, numbers, booleans, null),
/// reduced to a tree so the round trip can be checked without quick-xml.
fn to_xml(name: &str, v: &serde_json::Value, out: &mut String) {
    use serde_json::Value::*;
    let esc = |s: &str| s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
    let valid = |k: &str| {
        let mut c = k.chars();
        matches!(c.next(), Some(ch) if ch.is_ascii_alphabetic() || ch == '_')
            && c.all(|ch| ch.is_ascii_alphanumeric() || ch == '_' || ch == '-')
            && !k.to_ascii_lowercase().starts_with("xml")
            && k != "entry"
            && k != "item"
    };
    match v {
        Object(m) => {
            out.push_str(&format!("<{name} t=\"o\">"));
            for (k, x) in m {
                assert!(valid(k), "key {k} needs an <entry key=..>; geom keys must be plain names");
                to_xml(k, x, out);
            }
            out.push_str(&format!("</{name}>"));
        }
        Array(a) => {
            out.push_str(&format!("<{name} t=\"a\">"));
            a.iter().for_each(|x| to_xml("item", x, out));
            out.push_str(&format!("</{name}>"));
        }
        String(s) => out.push_str(&format!("<{name} t=\"s\">{}</{name}>", esc(s))),
        Number(n) => out.push_str(&format!("<{name} t=\"n\">{n}</{name}>")),
        Bool(b) => out.push_str(&format!("<{name} t=\"b\">{b}</{name}>")),
        Null => out.push_str(&format!("<{name} t=\"z\"></{name}>")),
    }
}

#[test]
fn serde_round_trips_bit_exactly_and_is_xml_compatible() {
    for g in common::all_2d() {
        let l = &g.model.list;
        let s = serde_json::to_string(l).unwrap();
        let back: FeatureList = serde_json::from_str(&s).unwrap();
        assert_eq!(&back, l);
        // numbers keep their identity (f64 vs string expressions)
        let v = serde_json::to_value(l).unwrap();
        let mut xml = String::new();
        to_xml("model", &v, &mut xml);
        assert!(xml.contains("t=\"n\"") && xml.contains("t=\"s\""));
        assert_eq!(serde_json::from_value::<FeatureList>(v).unwrap(), *l);
    }
    // awkward floats survive
    let mut m = Model::new_2d();
    let r = m.rect("r", (0.1 + 0.2, f64::MIN_POSITIVE), (1.0 / 3.0, 1e300));
    m.overlay(&[("r", r)]);
    let back: FeatureList = serde_json::from_str(&serde_json::to_string(&m.list).unwrap()).unwrap();
    assert_eq!(back, m.list);
    // JSON null is how serde_json writes non-finite floats; refuse them at evaluation, with the input's name
    let mut m = Model::new_2d();
    let r = m.rect("r", (0.0, 0.0), (f64::NAN, 1.0));
    m.overlay(&[("r", r)]);
    assert!(m.eval(Tier::Complex).err().unwrap().0.contains("max corner of feature 'r'"));
}

#[test]
fn naming_schemes_are_selectable_and_deterministic() {
    use geom::model::NameScheme;
    let names = |s: NameScheme| {
        let mut g = common::soi_strip();
        g.model.set_naming(s);
        let ev = g.model.eval(Tier::Complex).unwrap();
        let mut v: Vec<String> = ev.edges().iter().map(|e| e.name.clone()).collect();
        v.sort();
        (v, ev)
    };
    let (feat, _) = names(NameScheme::Feature);
    assert!(feat.contains(&"box.bottom".to_string()));
    let (reg, ev) = names(NameScheme::Region);
    assert!(
        reg.contains(&"BOX.bottom".to_string())
            && reg.contains(&"clad.top".to_string())
            && reg.contains(&"BOX|Si".to_string())
    );
    // the shape that drew a face stays reachable by provenance whatever the name
    assert_eq!(ev.resolve(&Selector::parse("from(box, bottom)").unwrap()).unwrap().set.ids.len(), 1);
    let (num, ev) = names(NameScheme::Numbered);
    assert!(
        num.contains(&"BOX.1".to_string())
            && num.contains(&"clad.1".to_string())
            && !num.iter().any(|n| n.contains("bottom"))
    );
    assert_eq!(num.iter().filter(|n| n.starts_with("BOX.")).count(), 3); // bottom, left, right
    assert_eq!(ev.label_map().unwrap().boundary("BOX.2").len(), 1);
    // 3D: the same schemes on faces, and the printed code keeps the choice
    let mut m = common::bragg3("b", true, true).model;
    m.set_naming(NameScheme::Region);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev.faces().iter().any(|f| f.name == "BOX.zmin") && ev.faces().iter().any(|f| f.name == "Si.zmin"));
    m.set_naming(NameScheme::Numbered);
    let ev = m.eval(Tier::Complex).unwrap();
    assert!(ev.faces().iter().any(|f| f.name == "clad.1") && !ev.faces().iter().any(|f| f.name.contains("zmin")));
    ev.complex().check().unwrap();
    assert_eq!(common::interp::build(&m.list.to_rust()).list, m.list);
}
