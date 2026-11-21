//! `FeatureList::to_rust`: prints the builder calls that recreate a list (the "Copy as Rust" of the GUI). One statement per feature,
//! `let fN = m.method(...)` where `N` is the feature id; operations without a builder shorthand fall back to `m.push(name, Op::...)`
//! in `Debug` form. Ids are assumed to be 1..n in order (as the builder assigns them).

use crate::expr::Expr;
use crate::model::*;

fn e(x: &Expr) -> String {
    match x {
        Expr::Num(v) => format!("{v:?}"),
        Expr::Text(s) => format!("{s:?}"),
    }
}
fn ex(x: &Expr) -> String {
    format!("Expr::from({})", e(x))
}
fn t2(v: &[Expr; 2]) -> String {
    format!("({}, {})", e(&v[0]), e(&v[1]))
}
fn t3(v: &[Expr; 3]) -> String {
    format!("({}, {}, {})", e(&v[0]), e(&v[1]), e(&v[2]))
}
fn ids(v: &[FeatureId]) -> String {
    v.iter().map(|i| format!("f{i}")).collect::<Vec<_>>().join(", ")
}
fn opt(x: &Option<Expr>) -> String {
    match x {
        Some(v) => format!("Some({})", ex(v)),
        None => "None".into(),
    }
}
fn strs(v: &[String]) -> String {
    v.iter().map(|s| format!("{s:?}")).collect::<Vec<_>>().join(", ")
}

pub fn to_rust(l: &FeatureList) -> String {
    let mut s = String::new();
    s.push_str(&format!("let mut m = Model::new_{}d();\n", l.dim));
    for p in &l.params {
        s.push_str(&format!("m.param({:?}, {});\n", p.name, e(&p.value)));
    }
    for f in &l.features {
        let n = format!("{:?}", f.name);
        let id = f.id;
        let naming = if let Op::Overlay { naming, .. } = &f.op { *naming } else { NameScheme::Feature };
        let call = match &f.op {
            Op::Rect { min, max } => format!("m.rect({n}, {}, {})", t2(min), t2(max)),
            Op::RoundedRect { min, max, r } => format!("m.rounded_rect({n}, {}, {}, {})", t2(min), t2(max), e(r)),
            Op::Circle { center, r } => format!("m.circle({n}, {}, {})", t2(center), e(r)),
            Op::Revolve { profile, axis: Axis::Z, center } => format!("m.revolve({n}, f{profile}, {})", t3(center)),
            Op::Revolve { profile, axis, center } => {
                format!("m.revolve_along({n}, Axis::{axis:?}, f{profile}, {})", t3(center))
            }
            Op::Cone { base, axis: Axis::Z, r0, r1, h } => {
                format!("m.cone({n}, {}, {}, {}, {})", t3(base), e(r0), e(r1), e(h))
            }
            Op::Cone { base, axis, r0, r1, h } => {
                format!("m.cone_along({n}, Axis::{axis:?}, {}, {}, {}, {})", t3(base), e(r0), e(r1), e(h))
            }
            Op::Sphere { center, r } => format!("m.sphere({n}, {}, {})", t3(center), e(r)),
            Op::Torus { center, major, minor } => format!("m.torus({n}, {}, {}, {})", t3(center), e(major), e(minor)),
            Op::Box { min, max } => format!("m.box_({n}, {}, {})", t3(min), t3(max)),
            Op::Cylinder { base, axis: Axis::Z, r, h } => format!("m.cylinder({n}, {}, {}, {})", t3(base), e(r), e(h)),
            Op::Extrude { profile, axis: Axis::Y, from: z0, to: z1, taper_deg: Some(t) } => {
                format!("m.extrude_tapered({n}, f{profile}, {}, {}, {})", e(z0), e(z1), e(t))
            }
            Op::Cylinder { base, axis, r, h } => {
                format!("m.cylinder_along({n}, Axis::{axis:?}, {}, {}, {})", t3(base), e(r), e(h))
            }
            Op::Extrude { profile, axis, from: z0, to: z1, taper_deg: None } if *axis != Axis::Y => {
                format!("m.extrude_along({n}, Axis::{axis:?}, f{profile}, {}, {})", e(z0), e(z1))
            }
            Op::Extrude { profile, axis: Axis::Y, from: z0, to: z1, taper_deg: None } => {
                format!("m.extrude({n}, f{profile}, {}, {})", e(z0), e(z1))
            }
            Op::Sheet { input, z, thickness, etch_depth } => {
                if *thickness == Expr::Num(0.0) && *etch_depth == Expr::Num(0.0) {
                    format!("m.sheet({n}, f{input}, {})", e(z))
                } else {
                    format!("m.sheet_thick({n}, f{input}, {}, {}, {})", e(z), e(thickness), e(etch_depth))
                }
            }
            Op::Ellipse { center, radii, angle_deg } => {
                format!("m.ellipse({n}, {}, {}, {})", t2(center), t2(radii), e(angle_deg))
            }
            Op::Polygon { pts, names } => {
                let p: Vec<String> = pts.iter().map(|q| format!("({}, {})", ex(&q[0]), ex(&q[1]))).collect();
                if names.is_empty() {
                    format!("m.polygon({n}, &[{}])", p.join(", "))
                } else {
                    format!("m.polygon_named({n}, &[{}], &[{}])", p.join(", "), strs(names))
                }
            }
            Op::Transform { input, steps } => match steps.as_slice() {
                [TStep::Translate(v)] => format!("m.translate({n}, f{input}, {})", t3(v)),
                _ => format!("m.push({n}, {})", dbg(&f.op)),
            },
            Op::Pattern { input, kind: PatternKind::Linear { step }, n: cnt } => {
                format!("m.linear_pattern({n}, f{input}, {}, {})", t3(step), e(cnt))
            }
            Op::Union { ids: v } => format!("m.union({n}, &[{}])", ids(v)),
            Op::Intersect { ids: v } => format!("m.intersect({n}, &[{}])", ids(v)),
            Op::Difference { a, b } => format!("m.difference({n}, f{a}, &[{}])", ids(b)),
            Op::Overlay { regions, pins, .. } if f.name == "overlay" => {
                if pins.is_empty() && regions.iter().all(|r| r.kind == RegionKind::Material && r.pinned_id.is_none()) {
                    let r: Vec<String> = regions.iter().map(|r| format!("({:?}, f{})", r.name, r.input)).collect();
                    format!("m.overlay(&[{}])", r.join(", "))
                } else {
                    let r: Vec<String> = regions
                        .iter()
                        .map(|r| {
                            format!(
                                "({:?}, f{}, RegionKind::{:?}, {})",
                                r.name,
                                r.input,
                                r.kind,
                                r.pinned_id.map_or("None".to_string(), |p| format!("Some({p})"))
                            )
                        })
                        .collect();
                    let p: Vec<String> = pins.iter().map(|p| format!("({:?}, {})", p.sel.to_string(), p.id)).collect();
                    format!("m.overlay_full(&[{}], &[{}])", r.join(", "), p.join(", "))
                }
            }
            Op::Pml { window, sides, thickness, layers, grading, regions } if regions.is_empty() => {
                let sd: Vec<String> = sides.iter().map(|s| format!("Side::{s:?}")).collect();
                format!("m.pml({n}, f{window}, &[{}], {}, {}, {})", sd.join(", "), e(thickness), e(layers), e(grading))
            }
            Op::Cut { normal, offset, keep } if normal.iter().all(|x| matches!(x, Expr::Num(_))) => {
                let c = |i: usize| e(&normal[i]);
                format!("m.cut({n}, ({}, {}, {}), {}, Keep::{keep:?})", c(0), c(1), c(2), e(offset))
            }
            Op::Periodic { sel_a, sel_b, translation } => format!(
                "m.periodic({n}, Selector::parse({:?})?, Selector::parse({:?})?, {})",
                sel_a.to_string(),
                sel_b.to_string(),
                t3(translation)
            ),
            Op::Size { sel, h, max_deviation, min_layers } => format!(
                "m.size({n}, Selector::parse({:?})?, {}, {}, {})",
                sel.to_string(),
                opt(h),
                opt(max_deviation),
                opt(min_layers)
            ),
            Op::Tag { sel, name } => format!("m.tag({n}, Selector::parse({:?})?, {name:?})", sel.to_string()),
            op => format!("m.push({n}, {})", dbg(op)),
        };
        s.push_str(&format!("let f{id} = {call};\n"));
        if naming != NameScheme::Feature {
            s.push_str(&format!("m.set_naming(NameScheme::{naming:?});\n"));
        }
        if f.suppressed {
            s.push_str(&format!("m.suppress(f{id}, true);\n"));
        }
        for t in &f.tags {
            s.push_str(&format!("m.tag_feature(f{id}, {t:?});\n"));
        }
    }
    s
}

/// `Debug` form of an operation (valid with `use geom::{expr::Expr::*, model::*, select::*}` in scope).
fn dbg(op: &Op) -> String {
    format!("{op:?}")
}
