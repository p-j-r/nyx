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
                prop_assert!((x - y).abs() <= 1e-7 * x.abs().max(y.abs()).max(1e-12), "{x} vs {y}");
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

#[test]
fn hairpin_near_collinear_polygon_is_cleanly_refused_or_valid() {
    let pts = [(1e-9, -1.0), (0.0, 2.6790266390872066), (0.0, 0.0), (-0.10259072690960207, 4.950451230554963)];
    let mut m = Model::new_2d();
    let epts: Vec<(Expr, Expr)> = pts.iter().map(|&(x, y)| (Expr::from(x), Expr::from(y))).collect();
    let p = m.polygon("s0", &epts);
    m.overlay(&[("R0", p)]);
    if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
        ev.complex().check().unwrap();
    }
}

#[test]
fn near_tolerance_union_area_is_associative_regression() {
    let a = Shape::Ellipse(4.0, 0.0, 0.1, 0.1, 0.0);
    let b = Shape::Poly(vec![(0.0, 0.0), (0.0, -1.0), (-3.5, 0.0)]);
    let c = Shape::Poly(vec![(-3.5, 1e-9), (0.0, -1.0), (-2.742422836991286, 0.0)]);
    let area = |order: [usize; 3], nested: bool| -> Option<f64> {
        let sh = [&a, &b, &c];
        let mut m = Model::new_2d();
        let ids: Vec<FeatureId> =
            order.iter().enumerate().map(|(i, &o)| add(&mut m, &format!("s{i}"), sh[o])).collect();
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
            assert!((x - y).abs() <= 1e-7 * x.abs().max(y.abs()).max(1e-12), "{x} vs {y}");
        }
    }
}

#[test]
fn pinched_csg_edge_between_rotated_boxes_passes_check3() {
    let boxes = [
        (-1.5090423952238703, 1e-9, 0.0, 0.5, 2.349053304024697, 0.5, 85.67248748751473, 2u8),
        (-3.5, 1e-9, 0.0, 0.5, 0.5, 0.5, 79.62470893893675, 1u8),
        (-3.5, 1e-9, 0.0, 0.5, 0.5, 0.5, 31.484207691471557, 1u8),
    ];
    let mut m = Model::new_3d();
    let ids: Vec<FeatureId> = boxes
        .iter()
        .enumerate()
        .map(|(i, &(x, y, z, a, b, c, deg, ax))| {
            let bx = m.box_(&format!("b{i}"), (x, y, z), (x + a, y + b, z + c));
            let axis = [Axis::X, Axis::Y, Axis::Z][ax as usize];
            m.transform(
                &format!("t{i}"),
                bx,
                vec![TStep::Rotate { axis, deg: deg.into(), about: [x.into(), y.into(), z.into()] }],
            )
        })
        .collect();
    let mut regs: Vec<(String, FeatureId)> = ids.iter().enumerate().map(|(i, &id)| (format!("R{i}"), id)).collect();
    let id = m.intersect("i", &ids[..2]);
    regs.push(("X".into(), id));
    let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
    m.overlay(&regs);
    if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
        ev.complex().check().unwrap();
        for r in ev.regions() {
            assert!(r.measure >= 0.0);
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
            for r in ev.regions() { prop_assert!(r.measure >= 0.0); }
        }
    }
}

mod vvuq_expansion {
    use geom::error::GeomError;
    use geom::eval::Tier;
    use geom::expr::Expr;
    use geom::gds::{encode_ibm_f64, parse_gds, GdsAref, GdsBoundary, GdsCell, GdsLibrary, GdsPath, GdsSref, GdsText};
    use geom::math::V3;
    use geom::model::{Axis, Keep, Model, TStep};
    use geom::sim::EntityRef;
    use proptest::prelude::*;

    fn ok_or_refusal<T>(r: Result<T, GeomError>) -> Option<T> {
        match r {
            Ok(x) => Some(x),
            Err(e) => {
                assert!(!e.0.is_empty());
                None
            }
        }
    }

    fn dist3(a: V3, b: V3) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    fn push_rec(out: &mut Vec<u8>, tag: u16, payload: &[u8]) {
        let len = 4 + payload.len() + (payload.len() % 2);
        out.extend_from_slice(&(len as u16).to_be_bytes());
        out.extend_from_slice(&tag.to_be_bytes());
        out.extend_from_slice(payload);
        if payload.len() % 2 == 1 {
            out.push(0);
        }
    }

    fn push_i16(out: &mut Vec<u8>, tag: u16, vals: &[i16]) {
        let mut p = Vec::with_capacity(vals.len() * 2);
        for &v in vals {
            p.extend_from_slice(&v.to_be_bytes());
        }
        push_rec(out, tag, &p);
    }

    fn push_i32(out: &mut Vec<u8>, tag: u16, vals: &[i32]) {
        let mut p = Vec::with_capacity(vals.len() * 4);
        for &v in vals {
            p.extend_from_slice(&v.to_be_bytes());
        }
        push_rec(out, tag, &p);
    }

    fn push_str(out: &mut Vec<u8>, tag: u16, s: &str) {
        push_rec(out, tag, s.as_bytes());
    }

    fn push_f64(out: &mut Vec<u8>, tag: u16, vals: &[f64]) {
        let mut p = Vec::with_capacity(vals.len() * 8);
        for &v in vals {
            p.extend_from_slice(&encode_ibm_f64(v));
        }
        push_rec(out, tag, &p);
    }

    fn encode_gds_library(lib: &GdsLibrary) -> Vec<u8> {
        let mut out = Vec::new();
        push_i16(&mut out, 0x0002, &[lib.version]); // HEADER
        push_i16(&mut out, 0x0102, &[2026, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0]); // BGNLIB
        push_str(&mut out, 0x0206, &lib.name); // LIBNAME
        push_f64(&mut out, 0x0305, &[lib.user_unit, lib.dbu_m]); // UNITS

        let ts = [2026i16, 1, 1, 0, 0, 0, 2026, 1, 1, 0, 0, 0];
        for cell in &lib.cells {
            push_i16(&mut out, 0x0502, &ts); // BGNSTR
            push_str(&mut out, 0x0606, &cell.name); // STRNAME
            for b in &cell.boundaries {
                push_rec(&mut out, 0x0800, &[]);
                push_i16(&mut out, 0x0D02, &[b.layer]);
                push_i16(&mut out, 0x0E02, &[b.datatype]);
                let xy: Vec<i32> = b.xy.iter().flat_map(|p| [p[0], p[1]]).collect();
                push_i32(&mut out, 0x1003, &xy);
                push_rec(&mut out, 0x1100, &[]);
            }
            for p in &cell.paths {
                push_rec(&mut out, 0x0900, &[]);
                push_i16(&mut out, 0x0D02, &[p.layer]);
                push_i16(&mut out, 0x0E02, &[p.datatype]);
                push_i16(&mut out, 0x2102, &[p.pathtype]);
                push_i32(&mut out, 0x0F03, &[p.width]);
                if p.pathtype == 4 {
                    push_i32(&mut out, 0x3003, &[p.bgnextn]);
                    push_i32(&mut out, 0x3103, &[p.endextn]);
                }
                let xy: Vec<i32> = p.xy.iter().flat_map(|pt| [pt[0], pt[1]]).collect();
                push_i32(&mut out, 0x1003, &xy);
                push_rec(&mut out, 0x1100, &[]);
            }
            for s in &cell.srefs {
                push_rec(&mut out, 0x0A00, &[]);
                push_str(&mut out, 0x1206, &s.sname);
                push_i16(&mut out, 0x1A01, &[s.strans as i16]);
                push_f64(&mut out, 0x1B05, &[s.mag]);
                push_f64(&mut out, 0x1C05, &[s.angle_deg]);
                push_i32(&mut out, 0x1003, &[s.xy[0], s.xy[1]]);
                push_rec(&mut out, 0x1100, &[]);
            }
            for a in &cell.arefs {
                push_rec(&mut out, 0x0B00, &[]);
                push_str(&mut out, 0x1206, &a.sname);
                push_i16(&mut out, 0x1A01, &[a.strans as i16]);
                push_f64(&mut out, 0x1B05, &[a.mag]);
                push_f64(&mut out, 0x1C05, &[a.angle_deg]);
                push_i16(&mut out, 0x1302, &[a.cols, a.rows]);
                let xy: Vec<i32> = a.xy.iter().flat_map(|pt| [pt[0], pt[1]]).collect();
                push_i32(&mut out, 0x1003, &xy);
                push_rec(&mut out, 0x1100, &[]);
            }
            for t in &cell.texts {
                push_rec(&mut out, 0x0C00, &[]);
                push_i16(&mut out, 0x0D02, &[t.layer]);
                push_i16(&mut out, 0x1602, &[t.texttype]);
                push_i16(&mut out, 0x1A01, &[t.strans as i16]);
                push_f64(&mut out, 0x1B05, &[t.mag]);
                push_f64(&mut out, 0x1C05, &[t.angle_deg]);
                push_i32(&mut out, 0x1003, &[t.xy[0], t.xy[1]]);
                push_str(&mut out, 0x1906, &t.string);
                push_rec(&mut out, 0x1100, &[]);
            }
            push_rec(&mut out, 0x0700, &[]);
        }
        push_rec(&mut out, 0x0400, &[]);
        out
    }

    type PolyBoxSpec = (i8, i8, i8, i8, i8, i8, u8, i16);

    fn build_random_3d_polyhedral_model(specs: &[PolyBoxSpec], op: u8, cut: Option<(i8, i8, i8, f64)>) -> Model {
        let mut m = Model::new_3d();
        let ids: Vec<_> = specs
            .iter()
            .enumerate()
            .map(|(i, &(x, y, z, a, b, c, ax, deg))| {
                let (x, y, z, a, b, c) = (x as f64, y as f64, z as f64, a as f64, b as f64, c as f64);
                let bx = m.box_(&format!("b{i}"), (x, y, z), (x + a, y + b, z + c));
                let axis = [Axis::X, Axis::Y, Axis::Z][ax as usize];
                m.transform(
                    &format!("t{i}"),
                    bx,
                    vec![TStep::Rotate {
                        axis,
                        deg: Expr::from(i32::from(deg)),
                        about: [Expr::from(x), Expr::from(y), Expr::from(z)],
                    }],
                )
            })
            .collect();
        let mut regs: Vec<(String, _)> = ids.iter().enumerate().map(|(i, &id)| (format!("R{i}"), id)).collect();
        if ids.len() > 1 {
            let id = match op {
                0 => m.union("u", &ids[..2]),
                1 => m.difference("d", ids[0], &ids[1..2]),
                _ => m.intersect("i", &ids[..2]),
            };
            regs.push(("X".into(), id));
        }
        let reg_refs: Vec<(&str, _)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
        m.overlay(&reg_refs);
        if let Some((a, b, c, off)) = cut {
            if (a, b, c) != (0, 0, 0) {
                m.cut("sym", (a as f64, b as f64, c as f64), off, Keep::Below);
            }
        }
        m
    }

    proptest! {
        #[test]
        fn op_sheet_thickness_and_etch_depth_never_panic(
            w0 in 20.0e-6_f64..80.0e-6,
            w1 in 20.0e-6_f64..80.0e-6,
            t0 in 0.0_f64..15.0e-6,
            e0 in 0.0_f64..15.0e-6,
            t1 in 0.0_f64..15.0e-6,
            e1 in 0.0_f64..15.0e-6,
            same_etch in any::<bool>(),
        ) {
            let mut m = Model::new_3d();
            let sub = m.box_("Si", (-200e-6, -200e-6, -100e-6), (200e-6, 200e-6, 0.0));
            let vac = m.box_("Vac", (-200e-6, -200e-6, 0.0), (200e-6, 200e-6, 100e-6));
            let p0 = m.rect("p0", (-20e-6 - w0, -40e-6), (-10e-6, 40e-6));
            let p1 = m.rect("p1", (10e-6, -40e-6), (20e-6 + w1, 40e-6));
            m.sheet_thick("c0", p0, 0.0, t0, e0);
            m.sheet_thick("c1", p1, 0.0, t1, if same_etch { e0 } else { e1 });
            m.overlay(&[("Si", sub), ("Vac", vac)]);
            if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
                ev.complex().check().unwrap();
                for r in ev.regions() {
                    prop_assert!(r.measure >= 0.0, "negative measure on {}: {}", r.name, r.measure);
                }
            }
        }

        #[test]
        fn op_polygons_with_arc_tol_gt_zero_preserves_check_and_area(
            n in 16usize..40,
            r in 1.0_f64..5.0,
            ry_ratio in 0.7_f64..1.3,
            simplify_tol in 0.0_f64..0.02,
            arc_tol in 1e-4_f64..0.05,
        ) {
            let pts_raw: Vec<[f64; 2]> = (0..n)
                .map(|k| {
                    let th = std::f64::consts::TAU * (k as f64) / (n as f64);
                    [r * th.cos(), r * ry_ratio * th.sin()]
                })
                .collect();
            let poly_area = {
                let mut sum = 0.0;
                for i in 0..pts_raw.len() {
                    let a = pts_raw[i];
                    let b = pts_raw[(i + 1) % pts_raw.len()];
                    sum += a[0] * b[1] - b[0] * a[1];
                }
                0.5 * sum.abs()
            };
            let pts: Vec<[Expr; 2]> = pts_raw.into_iter().map(|p| [Expr::from(p[0]), Expr::from(p[1])]).collect();

            let mut m = Model::new_2d();
            let id = m.polygons_with_arcs("arc_poly", vec![pts], vec![], simplify_tol, 0.0, arc_tol);
            m.overlay(&[("R", id)]);
            if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
                ev.complex().check().unwrap();
                prop_assert_eq!(ev.regions().len(), 1);
                let meas = ev.regions()[0].measure;
                prop_assert!(meas > 0.0);
                let true_area = std::f64::consts::PI * r * (r * ry_ratio);
                let tol = (arc_tol + simplify_tol) * r * std::f64::consts::TAU * 2.0 + 1e-6;
                let lo = poly_area.min(true_area) - tol;
                let hi = poly_area.max(true_area) + tol;
                prop_assert!(
                    meas >= lo && meas <= hi,
                    "meas={meas}, poly={poly_area}, ellipse={true_area}, bounds=[{lo}, {hi}]"
                );
            }
        }

        #[test]
        fn synthetic_gds_parse_roundtrip(
            version in 100i16..=600,
            layer in 0i16..=63,
            datatype in 0i16..=15,
            pathtype in prop_oneof![Just(0i16), Just(1), Just(2), Just(4)],
            width in prop_oneof![-500i32..=-10, 10i32..=1000],
            bgnextn in -200i32..=200,
            endextn in -200i32..=200,
            cols in 1i16..=8,
            rows in 1i16..=8,
            mag in 0.25_f64..4.0,
            angle_deg in -180.0_f64..180.0,
            dx in 100i32..=5000,
            dy in 100i32..=5000,
        ) {
            let lib = GdsLibrary {
                version,
                name: "SYNTH".into(),
                user_unit: 1e-3,
                dbu_m: 1e-9,
                cells: vec![
                    GdsCell {
                        name: "LEAF".into(),
                        boundaries: vec![GdsBoundary {
                            layer,
                            datatype,
                            xy: vec![[0, 0], [dx, 0], [dx, dy], [0, dy], [0, 0]],
                        }],
                        paths: vec![GdsPath {
                            layer,
                            datatype,
                            pathtype,
                            width,
                            bgnextn: if pathtype == 4 { bgnextn } else { 0 },
                            endextn: if pathtype == 4 { endextn } else { 0 },
                            xy: vec![[0, 0], [dx, 0], [dx, dy]],
                        }],
                        srefs: vec![],
                        arefs: vec![],
                        texts: vec![GdsText {
                            layer,
                            texttype: datatype,
                            string: "pinA".into(),
                            xy: [dx / 2, dy / 2],
                            strans: 0,
                            mag,
                            angle_deg,
                        }],
                    },
                    GdsCell {
                        name: "TOP".into(),
                        boundaries: vec![],
                        paths: vec![],
                        srefs: vec![GdsSref {
                            sname: "LEAF".into(),
                            strans: 0,
                            mag,
                            angle_deg,
                            xy: [100, -100],
                        }],
                        arefs: vec![GdsAref {
                            sname: "LEAF".into(),
                            strans: 0,
                            mag,
                            angle_deg,
                            cols,
                            rows,
                            xy: [[0, 0], [dx * 2 * i32::from(cols), 0], [0, dy * 2 * i32::from(rows)]],
                        }],
                        texts: vec![],
                    },
                ],
            };

            let bytes = encode_gds_library(&lib);
            let parsed = parse_gds(&bytes).expect("synthetic GDS bytes must parse cleanly");
            prop_assert_eq!(parsed.version, lib.version);
            prop_assert_eq!(&parsed.name, &lib.name);
            prop_assert!((parsed.user_unit - lib.user_unit).abs() <= 1e-12 * lib.user_unit);
            prop_assert!((parsed.dbu_m - lib.dbu_m).abs() <= 1e-12 * lib.dbu_m);
            prop_assert_eq!(parsed.cells.len(), 2);

            let p_leaf = &parsed.cells[0];
            prop_assert_eq!(&p_leaf.boundaries, &lib.cells[0].boundaries);
            prop_assert_eq!(&p_leaf.paths, &lib.cells[0].paths);
            prop_assert_eq!(p_leaf.texts.len(), 1);
            prop_assert_eq!(&p_leaf.texts[0].string, &lib.cells[0].texts[0].string);
            prop_assert!((p_leaf.texts[0].mag - mag).abs() <= 1e-12);

            let p_top = &parsed.cells[1];
            prop_assert_eq!(p_top.srefs.len(), 1);
            prop_assert_eq!(p_top.arefs.len(), 1);
            prop_assert_eq!(p_top.arefs[0].cols, cols);
            prop_assert_eq!(p_top.arefs[0].rows, rows);
            prop_assert!((p_top.srefs[0].mag - mag).abs() <= 1e-12);
            prop_assert!((p_top.srefs[0].angle_deg - angle_deg).abs() <= 1e-11);
        }

        #[test]
        fn facet_project_is_idempotent(
            specs in prop::collection::vec((-3i8..=3, -3i8..=3, -3i8..=3, 2i8..=6, 2i8..=6, 2i8..=6, 0u8..3, -90i16..=90), 1..=3),
            op in 0u8..3,
            queries in prop::collection::vec((-8.0_f64..8.0, -8.0_f64..8.0, -8.0_f64..8.0), 1..=10),
        ) {
            let m = build_random_3d_polyhedral_model(&specs, op, None);
            if let Some(ev) = ok_or_refusal(m.eval(Tier::Complex)) {
                let cx = ev.complex();
                let oracle = ev.oracle();
                for fi in 0..cx.faces.len().min(12) {
                    for &p in &queries {
                        let pt = [p.0, p.1, p.2];
                        let q1 = cx.facet_project(fi as u32, pt);
                        let q2 = cx.facet_project(fi as u32, q1);
                        let err = dist3(q1, q2);
                        prop_assert!(err <= 1e-12, "facet_project not idempotent: face={fi}, err={err}");
                        prop_assert_eq!(oracle.project(pt, EntityRef::Face(fi as u32)), Some(q1));
                    }
                }
            }
        }
    }

    #[cfg(feature = "parallel")]
    proptest! {
        #[test]
        fn parallel_eval_is_bit_identical_between_1_and_4_threads(
            specs in prop::collection::vec((-3i8..=3, -3i8..=3, -3i8..=3, 2i8..=6, 2i8..=6, 2i8..=6, 0u8..3, -90i16..=90), 1..=4),
            op in 0u8..3,
            cut in prop::option::of((-2i8..=2, -2i8..=2, -2i8..=2, -2.0_f64..2.0)),
        ) {
            let pool1 = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
            let pool4 = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
            let m1 = build_random_3d_polyhedral_model(&specs, op, cut);
            let m4 = build_random_3d_polyhedral_model(&specs, op, cut);

            let r1 = pool1.install(|| m1.eval(Tier::Complex));
            let r4 = pool4.install(|| m4.eval(Tier::Complex));
            match (r1, r4) {
                (Ok(ev1), Ok(ev4)) => {
                    let reg1 = ev1.regions();
                    let reg4 = ev4.regions();
                    prop_assert_eq!(reg1.len(), reg4.len());
                    for (a, b) in reg1.iter().zip(reg4.iter()) {
                        prop_assert_eq!(&a.name, &b.name);
                        prop_assert_eq!(a.measure.to_bits(), b.measure.to_bits());
                    }
                    let c1 = ev1.complex();
                    let c4 = ev4.complex();
                    prop_assert_eq!(c1.vertices.len(), c4.vertices.len());
                    prop_assert_eq!(c1.faces.len(), c4.faces.len());
                    for (v1, v4) in c1.vertices.iter().zip(c4.vertices.iter()) {
                        for k in 0..3 {
                            prop_assert_eq!(v1.p[k].to_bits(), v4.p[k].to_bits());
                        }
                    }
                }
                (Err(e1), Err(e4)) => prop_assert_eq!(e1.0, e4.0),
                _ => prop_assert!(false, "1-thread vs 4-thread eval diverged between Ok and Err"),
            }
        }
    }
}
