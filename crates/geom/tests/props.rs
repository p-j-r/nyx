//! Property tests: random feature lists never panic and never fail a self-check; union area is associative and commutative;
//! `A cap A = A`; unrelated features do not rename existing edges.

use geom::eval::Tier;
use geom::expr::Expr;
use geom::model::*;
use proptest::prelude::*;

#[derive(Clone, Debug)]
enum Shape {
    Rect(f64, f64, f64, f64),
    Circle(f64, f64, f64),
    Ellipse(f64, f64, f64, f64, f64),
    Poly(Vec<(f64, f64)>),
    Rounded(f64, f64, f64, f64, f64),
}

fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![(-5i32..5).prop_map(|v| v as f64), -5.0f64..5.0, Just(0.0), Just(1e-9), Just(-3.5)]
}

/// Coordinates without features near the tolerance: the tolerance is relative to the model size, so a far-away addition
/// legitimately changes what counts as "the same point" for input that is within a few `1e-9 L` of degenerate.
fn clean() -> impl Strategy<Value = f64> {
    prop_oneof![(-5i32..5).prop_map(|v| v as f64), (-50i32..50).prop_map(|v| v as f64 / 10.0)]
}

fn shape() -> impl Strategy<Value = Shape> {
    prop_oneof![
        (coord(), coord(), coord(), coord()).prop_map(|(a, b, c, d)| Shape::Rect(a, b, c, d)),
        (coord(), coord(), -1.0f64..3.0).prop_map(|(a, b, r)| Shape::Circle(a, b, r)),
        (coord(), coord(), 0.1f64..3.0, 0.1f64..3.0, 0.0f64..180.0)
            .prop_map(|(a, b, c, d, e)| Shape::Ellipse(a, b, c, d, e)),
        proptest::collection::vec((coord(), coord()), 0..6).prop_map(Shape::Poly),
        (coord(), coord(), 0.5f64..4.0, 0.5f64..4.0, 0.0f64..2.0).prop_map(|(a, b, c, d, r)| Shape::Rounded(
            a,
            b,
            a + c,
            b + d,
            r
        )),
    ]
}

fn add(m: &mut Model, name: &str, s: &Shape) -> FeatureId {
    match s {
        Shape::Rect(a, b, c, d) => m.rect(name, (*a, *b), (*c, *d)),
        Shape::Circle(a, b, r) => m.circle(name, (*a, *b), *r),
        Shape::Ellipse(a, b, c, d, e) => m.ellipse(name, (*a, *b), (*c, *d), *e),
        Shape::Poly(p) => {
            let p: Vec<(Expr, Expr)> = p.iter().map(|(x, y)| (Expr::from(*x), Expr::from(*y))).collect();
            m.polygon(name, &p)
        }
        Shape::Rounded(a, b, c, d, r) => m.rounded_rect(name, (*a, *b), (*c, *d), *r),
    }
}

fn ok_or_refusal(r: Result<geom::eval::Evaluated, geom::GeomError>) -> Option<geom::eval::Evaluated> {
    match r {
        Ok(e) => Some(e),
        Err(e) => {
            assert!(!e.0.contains("self-check") && !e.0.contains("internal error"), "{e}");
            None
        }
    }
}

fn clean_shape() -> impl Strategy<Value = Shape> {
    prop_oneof![
        (clean(), clean(), clean(), clean()).prop_map(|(a, b, c, d)| Shape::Rect(a, b, c, d)),
        (clean(), clean(), 0.1f64..3.0).prop_map(|(a, b, r)| Shape::Circle(a, b, r)),
        proptest::collection::vec((clean(), clean()), 3..6).prop_map(Shape::Poly),
    ]
}

fn solid() -> impl Strategy<Value = (f64, f64, f64, f64, f64, f64, bool)> {
    (coord(), coord(), coord(), 0.5f64..3.0, 0.5f64..3.0, 0.2f64..2.0, any::<bool>())
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 200, failure_persistence: None, ..ProptestConfig::default() })]

    #[test]
    fn random_3d_models_never_panic_and_stay_closed(solids in proptest::collection::vec(solid(), 1..5), op in 0u8..3) {
        let mut m = Model::new_3d();
        let ids: Vec<FeatureId> = solids
            .iter()
            .enumerate()
            .map(|(i, &(x, y, z, a, b, h, cyl))| {
                if cyl { m.cylinder(&format!("s{i}"), (x, y, z), a.min(b) * 0.5, h) } else { m.box_(&format!("s{i}"), (x, y, z), (x + a, y + b, z + h)) }
            })
            .collect();
        let mut regs: Vec<(String, FeatureId)> = ids.iter().enumerate().map(|(i, &id)| (format!("R{i}"), id)).collect();
        if ids.len() > 1 {
            let id = match op {
                0 => m.union("u", &ids[..2]),
                1 => m.difference("d", ids[0], &ids[1..2]),
                _ => m.intersect("i", &ids[..2]),
            };
            regs.push(("X".into(), id));
        }
        let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
        m.overlay(&regs);
        if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
            ev.complex().check().unwrap();
            prop_assert!(ev.layer_stack().is_some());
            for r in ev.regions() { prop_assert!(r.measure >= 0.0); }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 300, failure_persistence: None, ..ProptestConfig::default() })]

    #[test]
    fn random_feature_lists_never_panic(shapes in proptest::collection::vec(shape(), 1..6), ops in proptest::collection::vec(0u8..5, 0..4)) {
        let mut m = Model::new_2d();
        let ids: Vec<FeatureId> = shapes.iter().enumerate().map(|(i, s)| add(&mut m, &format!("s{i}"), s)).collect();
        let mut regs: Vec<(String, FeatureId)> = ids.iter().enumerate().map(|(i, &id)| (format!("R{i}"), id)).collect();
        for (k, op) in ops.iter().enumerate() {
            let (a, b) = (ids[k % ids.len()], ids[(k + 1) % ids.len()]);
            let id = match op {
                0 => m.union(&format!("u{k}"), &[a, b]),
                1 => m.difference(&format!("d{k}"), a, &[b]),
                2 => m.intersect(&format!("i{k}"), &[a, b]),
                3 => m.translate(&format!("t{k}"), a, (0.5, -0.25, 0.0)),
                _ => m.linear_pattern(&format!("p{k}"), a, (1.0, 0.0, 0.0), 2.0),
            };
            regs.push((format!("X{k}"), id));
        }
        let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
        m.overlay(&regs);
        if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
            ev.complex().check().unwrap();
            for r in ev.regions() { prop_assert!(r.measure >= 0.0); }
            let lm = ev.label_map().unwrap();
            prop_assert_eq!(lm.boundaries.iter().map(|b| b.entities.len()).sum::<usize>(), ev.edges().len());
        }
    }

    #[test]
    fn union_area_is_associative_and_commutative(a in shape(), b in shape(), c in shape()) {
        let area = |order: [usize; 3], nested: bool| -> Option<f64> {
            let sh = [&a, &b, &c];
            let mut m = Model::new_2d();
            let ids: Vec<FeatureId> = order.iter().enumerate().map(|(i, &o)| add(&mut m, &format!("s{i}"), sh[o])).collect();
            let u = if nested {
                let ab = m.union("ab", &[ids[0], ids[1]]);
                m.union("abc", &[ab, ids[2]])
            } else {
                let bc = m.union("bc", &[ids[1], ids[2]]);
                m.union("abc", &[ids[0], bc])
            };
            m.overlay(&[("u", u)]);
            let ev = ok_or_refusal(m.eval(Tier::Complex))?;
            Some(ev.regions().iter().map(|r| r.measure).sum())
        };
        let base = area([0, 1, 2], true);
        for (o, n) in [([0, 1, 2], false), ([2, 0, 1], true), ([1, 2, 0], false), ([2, 1, 0], true)] {
            if let (Some(x), Some(y)) = (base, area(o, n)) {
                prop_assert!((x - y).abs() <= 1e-9 * x.abs().max(y.abs()).max(1e-12), "{x} vs {y}");
            }
        }
    }

    #[test]
    fn a_cap_a_is_a(a in shape()) {
        let mut m = Model::new_2d();
        let id = add(&mut m, "a", &a);
        let i = m.intersect("aa", &[id, id]);
        m.overlay(&[("i", i)]);
        let mut m2 = Model::new_2d();
        let id2 = add(&mut m2, "a", &a);
        m2.overlay(&[("i", id2)]);
        match (ok_or_refusal(m.eval(Tier::Complex)), ok_or_refusal(m2.eval(Tier::Complex))) {
            (Some(x), Some(y)) => prop_assert_eq!(x.regions()[0].measure.to_bits(), y.regions()[0].measure.to_bits()),
            (None, None) => {}
            _ => prop_assert!(false, "one side was refused and the other not"),
        }
    }

    #[test]
    fn an_unrelated_feature_keeps_existing_names_and_ids(a in clean_shape(), b in clean_shape()) {
        let build = |extra: bool| {
            let mut m = Model::new_2d();
            let x = add(&mut m, "a", &a);
            let y = add(&mut m, "b", &b);
            let mut regs = vec![("A", x), ("B", y)];
            if extra {
                let far = m.rect("far", (7.0, 7.0), (8.0, 8.5));
                regs.push(("FAR", far));
            }
            m.overlay(&regs);
            ok_or_refusal(m.eval(Tier::Complex))
        };
        if let (Some(base), Some(more)) = (build(false), build(true)) {
            // the far rectangle's edges sort after everything else (last feature), so every existing id keeps its name
            for (i, e) in base.edges().iter().enumerate() {
                prop_assert_eq!(&more.edges()[i].name, &e.name);
            }
            prop_assert_eq!(more.regions().len(), base.regions().len() + 1);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 150, failure_persistence: None, ..ProptestConfig::default() })]

    #[test]
    fn random_polyhedral_models_never_panic_and_stay_closed(
        boxes in proptest::collection::vec((coord(), coord(), coord(), 0.5f64..2.5, 0.5f64..2.5, 0.5f64..2.5, 0.0f64..90.0, 0u8..3), 1..4),
        op in 0u8..3,
        cut in proptest::option::of((-1i32..2, -1i32..2, -1i32..2, -2.0f64..2.0)),
    ) {
        let mut m = Model::new_3d();
        let ids: Vec<FeatureId> = boxes
            .iter()
            .enumerate()
            .map(|(i, &(x, y, z, a, b, c, deg, ax))| {
                let bx = m.box_(&format!("b{i}"), (x, y, z), (x + a, y + b, z + c));
                let axis = [Axis::X, Axis::Y, Axis::Z][ax as usize];
                m.transform(&format!("t{i}"), bx, vec![TStep::Rotate { axis, deg: deg.into(), about: [x.into(), y.into(), z.into()] }])
            })
            .collect();
        let mut regs: Vec<(String, FeatureId)> = ids.iter().enumerate().map(|(i, &id)| (format!("R{i}"), id)).collect();
        if ids.len() > 1 {
            let id = match op {
                0 => m.union("u", &ids[..2]),
                1 => m.difference("d", ids[0], &ids[1..2]),
                _ => m.intersect("i", &ids[..2]),
            };
            regs.push(("X".into(), id));
        }
        let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
        m.overlay(&regs);
        if let Some((a, b, c, off)) = cut {
            if (a, b, c) != (0, 0, 0) {
                m.cut("sym", (a as f64, b as f64, c as f64), off, Keep::Below);
            }
        }
        if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
            ev.complex().check().unwrap();
            for r in ev.regions() { prop_assert!(r.measure >= -1e-9); }
        }
    }
}
