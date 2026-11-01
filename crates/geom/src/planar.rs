//! General planar-faced solids: convex polyhedra (boxes in any orientation, tapered extrusions, half-space cuts at any angle)
//! combined by Booleans. Used when a model is not a stack of z-prisms (a rotation about x or y, a tapered extrusion, a slanted `Cut`).
//!
//! The arrangement reuses the 2D one plane by plane. For every plane that carries a face, each solid contributes its cross-section with
//! that plane (the face itself if it lies on the plane, otherwise the convex hull of the crossings); the 2D arrangement of those sections,
//! with the composite class (region above, region below), gives the cells of the plane, and the cells whose two sides differ are
//! the faces. Edges are the cell boundaries merged by 3D line (so T-junctions between planes are resolved), vertices their end points.
//! Region volumes follow from the divergence theorem; the closure and edge incidence checks apply as for stacks.

use crate::arrange::arrange_with;
use crate::complex::*;
use crate::error::{GeomError, GeomResult};
use crate::expr::Env;
use crate::math::*;
use crate::model::*;
use crate::prim::*;
use crate::scene::{Cache, Node, RegionX, Scene};
use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;

/// True when the model needs the polyhedral path.
pub fn needs_planar(list: &FeatureList) -> bool {
    let sweep = crate::scene::stack_axis(list).unwrap_or(2);
    list.dim == 3
        && list.features.iter().filter(|f| !f.suppressed).any(|f| match &f.op {
            Op::Transform { steps, .. } => {
                steps.iter().any(|s| matches!(s, TStep::Rotate { axis, .. } if axis.index() != sweep))
            }
            Op::Pattern { kind: PatternKind::Polar { axis, .. }, .. } => axis.index() != sweep,
            Op::Extrude { taper_deg, .. } => taper_deg.is_some(),
            Op::Cut { normal, .. } => {
                let nz = |e: &crate::expr::Expr| !matches!(e, crate::expr::Expr::Num(v) if *v == 0.0);
                nz(&normal[2]) as u8 + nz(&normal[0]) as u8 + nz(&normal[1]) as u8 > 1
            }
            _ => false,
        })
}

struct PB<'a> {
    list: &'a FeatureList,
    env: &'a Env,
    prims: Vec<Arc<Prim>>,
    memo: HashMap<(FeatureId, [u64; 12], String), usize>,
    index: HashMap<FeatureId, usize>,
    depth: usize,
}

fn rings_of_box(lo: V3, hi: V3) -> Vec<(Vec<V3>, String)> {
    let c = |x: f64, y: f64, z: f64| [x, y, z];
    let (x0, y0, z0, x1, y1, z1) = (lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
    vec![
        (vec![c(x0, y0, z0), c(x0, y1, z0), c(x0, y1, z1), c(x0, y0, z1)], "xmin".into()),
        (vec![c(x1, y0, z0), c(x1, y1, z0), c(x1, y1, z1), c(x1, y0, z1)], "xmax".into()),
        (vec![c(x0, y0, z0), c(x1, y0, z0), c(x1, y0, z1), c(x0, y0, z1)], "ymin".into()),
        (vec![c(x0, y1, z0), c(x1, y1, z0), c(x1, y1, z1), c(x0, y1, z1)], "ymax".into()),
        (vec![c(x0, y0, z0), c(x1, y0, z0), c(x1, y1, z0), c(x0, y1, z0)], "zmin".into()),
        (vec![c(x0, y0, z1), c(x1, y0, z1), c(x1, y1, z1), c(x0, y1, z1)], "zmax".into()),
    ]
}

/// Side walls, bottom and top of a convex counter-clockwise polygon swept from `z0` to `z1`, the top offset inwards by `off`.
fn prism_rings(
    pts: &[V2],
    roles: &[String],
    z0: f64,
    z1: f64,
    off: f64,
    who: &str,
) -> GeomResult<Vec<(Vec<V3>, String)>> {
    let n = pts.len();
    let top: Vec<V2> = if off == 0.0 {
        pts.to_vec()
    } else {
        // move every edge inwards by `off` and intersect neighbouring offset lines
        let lines: Vec<(V2, V2)> = (0..n)
            .map(|i| {
                let (a, b) = (pts[i], pts[(i + 1) % n]);
                let d = sub2(b, a);
                let nrm = mul2(perp(d), 1.0 / norm2(d));
                (add2(a, mul2(nrm, off)), d)
            })
            .collect();
        let mut out = vec![];
        for i in 0..n {
            let (p, d) = lines[(i + n - 1) % n];
            let (q, e) = lines[i];
            let den = cross2(e, d);
            if den.abs() < 1e-12 * norm2(d) * norm2(e) {
                return Err(GeomError(format!(
                    "{who}: two neighbouring edges are parallel; a tapered extrusion needs corners"
                )));
            }
            out.push(add2(q, mul2(e, cross2(sub2(p, q), d) / den)));
        }
        let area = polygon_area(&out);
        let flipped = (0..n).any(|i| dot2(sub2(out[(i + 1) % n], out[i]), sub2(pts[(i + 1) % n], pts[i])) <= 0.0);
        if area <= 1e-9 * polygon_area(pts) || flipped {
            return Err(GeomError(format!(
                "{who}: the taper closes the profile before reaching the top; use a smaller angle or height"
            )));
        }
        out
    };
    let v3 = |p: V2, z: f64| [p[0], p[1], z];
    let mut rings = vec![
        (pts.iter().map(|p| v3(*p, z0)).collect(), "bottom".to_string()),
        (top.iter().map(|p| v3(*p, z1)).collect(), "top".to_string()),
    ];
    for i in 0..n {
        let j = (i + 1) % n;
        rings
            .push((vec![v3(pts[i], z0), v3(pts[j], z0), v3(top[j], z1), v3(top[i], z1)], format!("side:{}", roles[i])));
    }
    Ok(rings)
}

/// Rings built in the profile frame `(u, v, w)` of a sweep along `axis`, in world coordinates, with the caps named after the axis.
fn world_rings(rings: Vec<(Vec<V3>, String)>, axis: Axis) -> Vec<(Vec<V3>, String)> {
    let perm = crate::scene::axis_perm(axis.index());
    rings
        .into_iter()
        .map(|(r, role)| {
            let role = match role.as_str() {
                "bottom" => crate::scene::cap_role(axis.index(), false),
                "top" => crate::scene::cap_role(axis.index(), true),
                _ => role,
            };
            let r = r
                .into_iter()
                .map(|p| {
                    let mut o = [0.0; 3];
                    for i in 0..3 {
                        o[perm[i]] = p[i];
                    }
                    o
                })
                .collect();
            (r, role)
        })
        .collect()
}

impl PB<'_> {
    fn ev(&self, e: &crate::expr::Expr, f: &Feature, what: &str) -> GeomResult<f64> {
        self.env.eval(e, &format!("{what} of feature '{}'", f.name))
    }
    fn ev3(&self, e: &[crate::expr::Expr; 3], f: &Feature, what: &str) -> GeomResult<V3> {
        Ok([self.ev(&e[0], f, what)?, self.ev(&e[1], f, what)?, self.ev(&e[2], f, what)?])
    }

    fn leaf(&mut self, f: &Feature, rings: Vec<(Vec<V3>, String)>, xf: &Aff3, suffix: &str, tag: &str) -> Node {
        let key = (f.id, xf.key(), format!("{suffix}{tag}"));
        if let Some(&i) = self.memo.get(&key) {
            return Node::Prim(i);
        }
        let rings = rings.into_iter().map(|(r, role)| (r.into_iter().map(|p| xf.apply(p)).collect(), role)).collect();
        self.prims.push(Arc::new(Prim::new_poly(f.id, format!("{}{suffix}", f.name), Poly3::new(rings), false)));
        self.memo.insert(key, self.prims.len() - 1);
        Node::Prim(self.prims.len() - 1)
    }

    fn step(&self, s: &TStep, f: &Feature) -> GeomResult<Aff3> {
        let axis = |a: Axis| match a {
            Axis::X => [1.0, 0.0, 0.0],
            Axis::Y => [0.0, 1.0, 0.0],
            Axis::Z => [0.0, 0.0, 1.0],
        };
        Ok(match s {
            TStep::Translate(v) => Aff3::translate(self.ev3(v, f, "translation")?),
            TStep::Rotate { axis: a, deg, about } => Aff3::rotation(
                axis(*a),
                self.ev(deg, f, "rotation angle")?.to_radians(),
                self.ev3(about, f, "rotation centre")?,
            ),
            TStep::Mirror { axis: a, at } => Aff3::mirror(a.index(), self.ev(at, f, "mirror position")?),
            TStep::Scale { factor, about } => {
                let k = self.ev(factor, f, "scale factor")?;
                if k <= 0.0 {
                    return Err(GeomError(format!("Transform '{}' has scale factor {k}; it must be positive", f.name)));
                }
                Aff3::scaling(k, self.ev3(about, f, "scale centre")?)
            }
        })
    }

    fn inst(&mut self, id: FeatureId, xf: &Aff3, suffix: &str) -> GeomResult<Node> {
        let Some(&fi) = self.index.get(&id) else { return Err(GeomError(format!("unknown feature id {id}"))) };
        let list = self.list;
        let f = &list.features[fi];
        if f.suppressed {
            return Ok(Node::Empty);
        }
        self.depth += 1;
        if self.depth > 256 {
            return Err(GeomError(format!("feature '{}' is nested more than 256 levels deep", f.name)));
        }
        let r = self.inst_op(f, xf, suffix);
        self.depth -= 1;
        r
    }

    /// The counter-clockwise polygon of a `Rect` or `Polygon` profile with its edge roles.
    fn profile(&self, id: FeatureId) -> GeomResult<(Vec<V2>, Vec<String>, String)> {
        let Some(&fi) = self.index.get(&id) else { return Err(GeomError(format!("unknown feature id {id}"))) };
        let f = &self.list.features[fi];
        let ev2 = |e: &[crate::expr::Expr; 2], w: &str| -> GeomResult<V2> {
            Ok([self.ev(&e[0], f, w)?, self.ev(&e[1], f, w)?])
        };
        match &f.op {
            Op::Rect { min, max } => {
                let (a, b) = (ev2(min, "min corner")?, ev2(max, "max corner")?);
                let (lo, hi) = ([a[0].min(b[0]), a[1].min(b[1])], [a[0].max(b[0]), a[1].max(b[1])]);
                if hi[0] <= lo[0] || hi[1] <= lo[1] {
                    return Err(GeomError(format!("Rect '{}' has zero or negative width or height", f.name)));
                }
                Ok((vec![lo, [hi[0], lo[1]], hi, [lo[0], hi[1]]], ["bottom", "right", "top", "left"].map(String::from).to_vec(), f.name.clone()))
            }
            Op::Polygon { pts, names } => {
                if pts.len() < 3 {
                    return Err(GeomError(format!("Polygon '{}' needs at least 3 points", f.name)));
                }
                let p: Vec<V2> = pts.iter().map(|q| ev2(q, "point")).collect::<GeomResult<_>>()?;
                let nm: Vec<String> = if names.is_empty() { (0..p.len()).map(|i| format!("e{i}")).collect() } else { names.clone() };
                if nm.len() != p.len() || polygon_area(&p).abs() <= 1e-12 {
                    return Err(GeomError(format!("Polygon '{}' has a wrong number of edge names or zero area", f.name)));
                }
                let (p, nm) = ccw_polygon(&p, &nm);
                Ok((p, nm, f.name.clone()))
            }
            op => Err(GeomError(format!(
                "{} '{}' cannot be extruded together with rotated, tapered or slanted geometry: only Rect and Polygon profiles are polyhedral",
                op.kind(),
                f.name
            ))),
        }
    }

    fn inst_op(&mut self, f: &Feature, xf: &Aff3, suffix: &str) -> GeomResult<Node> {
        match &f.op {
            Op::Box { min, max } => {
                let (a, b) = (self.ev3(min, f, "min corner")?, self.ev3(max, f, "max corner")?);
                let (lo, hi) = ([a[0].min(b[0]), a[1].min(b[1]), a[2].min(b[2])], [a[0].max(b[0]), a[1].max(b[1]), a[2].max(b[2])]);
                if (0..3).any(|k| hi[k] <= lo[k]) {
                    return Err(GeomError(format!("Box '{}' has zero or negative extent", f.name)));
                }
                Ok(self.leaf(f, rings_of_box(lo, hi), xf, suffix, ""))
            }
            Op::Extrude { profile, axis, from: z0, to: z1, taper_deg } => {
                let (a, b) = (self.ev(z0, f, "z0")?, self.ev(z1, f, "z1")?);
                if (b - a).abs() <= 0.0 {
                    return Err(GeomError(format!("Extrude '{}' has zero height", f.name)));
                }
                let (z0, z1) = (a.min(b), a.max(b));
                let (pts, roles, _) = self.profile(*profile)?;
                let taper = match taper_deg {
                    Some(t) => self.ev(t, f, "taper angle")?,
                    None => 0.0,
                };
                let convex = (0..pts.len()).all(|i| cross2(sub2(pts[(i + 1) % pts.len()], pts[i]), sub2(pts[(i + 2) % pts.len()], pts[(i + 1) % pts.len()])) >= -1e-14);
                if convex {
                    let off = (z1 - z0) * taper.to_radians().tan();
                    let rings = world_rings(prism_rings(&pts, &roles, z0, z1, off, &format!("Extrude '{}'", f.name))?, *axis);
                    return Ok(self.leaf(f, rings, xf, suffix, ""));
                }
                if taper != 0.0 {
                    return Err(GeomError(format!("Extrude '{}': a tapered extrusion needs a convex profile", f.name)));
                }
                // concave: a union of convex prisms over a triangulation; the diagonals are interior and vanish in the arrangement
                let p3: Vec<V3> = pts.iter().map(|p| [p[0], p[1], 0.0]).collect();
                let tris = crate::tess::triangulate(&p3, (0..pts.len() as u32).collect(), vec![]);
                let mut parts = vec![];
                for (k, t) in tris.iter().enumerate() {
                    let idx = [t[0] as usize, t[1] as usize, t[2] as usize];
                    let tp: Vec<V2> = idx.iter().map(|&i| pts[i]).collect();
                    let tr: Vec<String> = (0..3)
                        .map(|e| {
                            let (i, j) = (idx[e], idx[(e + 1) % 3]);
                            if (i + 1) % pts.len() == j {
                                roles[i].clone()
                            } else {
                                "diagonal".to_string()
                            }
                        })
                        .collect();
                    let rings = world_rings(prism_rings(&tp, &tr, z0, z1, 0.0, &format!("Extrude '{}'", f.name))?, *axis);
                    parts.push(self.leaf(f, rings, xf, suffix, &format!("#{k}")));
                }
                Ok(Node::Union(parts))
            }
            Op::Transform { input, steps } => {
                let mut t = Aff3::ID;
                for s in steps {
                    t = self.step(s, f)?.after(&t);
                }
                self.inst(*input, &xf.after(&t), suffix)
            }
            Op::Pattern { input, kind, n } => {
                let nn = self.ev(n, f, "count")?;
                if nn.fract() != 0.0 || !(1.0..=100_000.0).contains(&nn) {
                    return Err(GeomError(format!("Pattern '{}' count {nn} must be a whole number between 1 and 100000", f.name)));
                }
                let nn = nn as usize;
                let mut out = vec![];
                for i in 0..nn {
                    let t = match kind {
                        PatternKind::Linear { step } => Aff3::translate(mul3(self.ev3(step, f, "step")?, i as f64)),
                        PatternKind::Polar { axis, center, total_deg } => {
                            let k = match axis {
                                Axis::X => [1.0, 0.0, 0.0],
                                Axis::Y => [0.0, 1.0, 0.0],
                                Axis::Z => [0.0, 0.0, 1.0],
                            };
                            Aff3::rotation(k, (self.ev(total_deg, f, "angle")? * i as f64 / nn as f64).to_radians(), self.ev3(center, f, "centre")?)
                        }
                    };
                    out.push(self.inst(*input, &xf.after(&t), &format!("{suffix}[{i}]"))?);
                }
                Ok(Node::Union(out))
            }
            Op::Union { ids } => {
                let v = ids.iter().map(|i| self.inst(*i, xf, suffix)).collect::<GeomResult<Vec<_>>>()?;
                Ok(Node::Union(v.into_iter().filter(|n| *n != Node::Empty).collect()))
            }
            Op::Intersect { ids } => {
                let v: Vec<Node> = ids.iter().map(|i| self.inst(*i, xf, suffix)).collect::<GeomResult<Vec<_>>>()?.into_iter().filter(|n| *n != Node::Empty).collect();
                Ok(if v.is_empty() { Node::Empty } else { Node::Inter(v) })
            }
            Op::Difference { a, b } => {
                let a = self.inst(*a, xf, suffix)?;
                if a == Node::Empty {
                    return Ok(a);
                }
                let b = b.iter().map(|i| self.inst(*i, xf, suffix)).collect::<GeomResult<Vec<_>>>()?;
                Ok(Node::Diff(Box::new(a), b.into_iter().filter(|n| *n != Node::Empty).collect()))
            }
            Op::Rect { .. } | Op::Polygon { .. } => Err(GeomError(format!("{} '{}' is a 2D profile: in a 3D model it must be the input of an Extrude", f.op.kind(), f.name))),
            op => Err(GeomError(format!(
                "{} '{}' cannot be combined with rotated, tapered or slanted geometry in this build (curved surfaces need the prismatic or the body path)",
                op.kind(),
                f.name
            ))),
        }
    }
}

/// Scene of a polyhedral model.
pub fn build_scene(list: &FeatureList, cache: &mut Cache) -> GeomResult<Scene> {
    let _ = cache;
    list.validate()?;
    let mut env = Env::default();
    for p in &list.params {
        let v = env.eval(&p.value, &format!("parameter '{}'", p.name))?;
        env.vars.push((p.name.clone(), v));
    }
    let overlay = list
        .features
        .iter()
        .find(|f| matches!(f.op, Op::Overlay { .. }) && !f.suppressed)
        .ok_or_else(|| GeomError("the model has no Overlay feature; add one that lists the regions".into()))?;
    let Op::Overlay { regions, pins, naming } = &overlay.op else { unreachable!() };
    let index = list.features.iter().enumerate().map(|(i, f)| (f.id, i)).collect();
    let mut b = PB { list, env: &env, prims: vec![], memo: HashMap::new(), index, depth: 0 };
    let mut rx: Vec<RegionX> = vec![];
    for r in regions {
        let expr = b.inst(r.input, &Aff3::ID, "")?;
        rx.push(RegionX {
            name: r.name.clone(),
            kind: r.kind,
            feature: r.input,
            pinned_id: r.pinned_id,
            expr,
            parent: None,
            pml: None,
        });
    }
    let model_box = |prims: &[Arc<Prim>]| prims.iter().filter(|p| !p.aux).fold(Aabb::EMPTY, |a, p| a.union(&p.bbox3()));
    let diag = model_box(&b.prims).diagonal();
    let after = list.features.iter().position(|f| f.id == overlay.id).unwrap_or(0);
    let mut auto_tags = vec![];
    for f in list.features[after..].iter().filter(|f| !f.suppressed) {
        match &f.op {
            Op::Cut { normal, offset, keep } => {
                let n = b.ev3(normal, f, "normal")?;
                let nl = norm3(n);
                if nl < 1e-12 {
                    return Err(GeomError(format!("Cut '{}' needs a nonzero normal", f.name)));
                }
                let (n, off) = (mul3(n, 1.0 / nl), b.ev(offset, f, "offset")? / nl);
                let sgn = if *keep == Keep::Above { 1.0 } else { -1.0 };
                let a = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                let u = unit3(cross3(n, a));
                let v = cross3(n, u);
                let big = 8.0 * diag + 8.0 * off.abs() + 1.0;
                let q = mul3(n, off);
                let corner = |x: f64, y: f64, z: f64| {
                    add3(add3(q, add3(mul3(u, x * big), mul3(v, y * big))), mul3(n, sgn * z * big))
                };
                let c = |i: usize, j: usize, k: usize| corner([-1.0, 1.0][i], [-1.0, 1.0][j], k as f64);
                let rings = vec![
                    (vec![c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0)], "cut".to_string()),
                    (vec![c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)], "far0".to_string()),
                    (vec![c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(0, 0, 1)], "far1".to_string()),
                    (vec![c(0, 1, 0), c(1, 1, 0), c(1, 1, 1), c(0, 1, 1)], "far2".to_string()),
                    (vec![c(0, 0, 0), c(0, 1, 0), c(0, 1, 1), c(0, 0, 1)], "far3".to_string()),
                    (vec![c(1, 0, 0), c(1, 1, 0), c(1, 1, 1), c(1, 0, 1)], "far4".to_string()),
                ];
                b.prims.push(Arc::new(Prim::new_poly(f.id, f.name.clone(), Poly3::new(rings), true)));
                let h = Node::Prim(b.prims.len() - 1);
                for r in &mut rx {
                    if r.expr != Node::Empty {
                        r.expr = Node::Inter(vec![std::mem::replace(&mut r.expr, Node::Empty), h.clone()]);
                    }
                }
                auto_tags.push((f.id, "cut".to_string(), "symmetry".to_string()));
            }
            Op::Pml { window, sides, thickness, layers, grading, regions: only } => {
                let (t, nl, gr) =
                    (b.ev(thickness, f, "thickness")?, b.ev(layers, f, "layers")?, b.ev(grading, f, "grading")?);
                if t <= 0.0 || nl < 1.0 || gr <= 0.0 {
                    return Err(GeomError(format!(
                        "Pml '{}' needs thickness > 0, layers >= 1 and grading > 0",
                        f.name
                    )));
                }
                let mark = b.prims.len();
                let wnode = b.inst(*window, &Aff3::ID, "")?;
                let win = wnode.bbox(&b.prims);
                b.prims.truncate(mark);
                b.memo.retain(|_, v| *v < mark);
                if win.is_empty() {
                    return Err(GeomError(format!("Pml '{}' has an empty window", f.name)));
                }
                let all = [Side::XMin, Side::XMax, Side::YMin, Side::YMax, Side::ZMin, Side::ZMax];
                let sel: Vec<Side> = if sides.is_empty() { all.to_vec() } else { sides.clone() };
                let mut combos: Vec<[i8; 3]> = vec![];
                for code in 0..27 {
                    let (mut c, mut k) = ([0i8; 3], code);
                    for ci in c.iter_mut() {
                        *ci = (k % 3) as i8 - 1;
                        k /= 3;
                    }
                    if c != [0, 0, 0]
                        && (0..3)
                            .all(|a| c[a] == 0 || sel.contains(&if c[a] < 0 { all[2 * a] } else { all[2 * a + 1] }))
                    {
                        combos.push(c);
                    }
                }
                combos.sort_by_key(|c| {
                    let n = c.iter().filter(|&&x| x != 0).count();
                    let single = (0..3).find(|&a| c[a] != 0).map_or(0, |a| 2 * a + (c[a] > 0) as usize);
                    (n, if n == 1 { single } else { 0 }, (0..3).rev().map(|a| c[a] + 1).collect::<Vec<_>>())
                });
                let mut pieces: Vec<(Node, [i8; 3], String)> = vec![];
                for c in &combos {
                    let rng = |k: usize| match c[k] {
                        -1 => (win.lo[k] - t, win.lo[k]),
                        0 => (win.lo[k], win.hi[k]),
                        _ => (win.hi[k], win.hi[k] + t),
                    };
                    let piece = (0..3)
                        .filter(|&a| c[a] != 0)
                        .map(|a| format!("{}{}", ["x", "y", "z"][a], if c[a] < 0 { "min" } else { "max" }))
                        .collect::<Vec<_>>()
                        .join("_");
                    let (x, y, z) = (rng(0), rng(1), rng(2));
                    b.prims.push(Arc::new(Prim::new_poly(
                        f.id,
                        format!("{}[{piece}]", f.name),
                        Poly3::new(rings_of_box([x.0, y.0, z.0], [x.1, y.1, z.1])),
                        true,
                    )));
                    pieces.push((Node::Prim(b.prims.len() - 1), *c, piece));
                }
                let pad_bb: Vec<Aabb> = pieces.iter().map(|(n, _, _)| n.bbox(&b.prims)).collect();
                let mut out: Vec<RegionX> = vec![];
                for r in rx.iter() {
                    let candidate = r.pml.is_none()
                        && r.parent.is_none()
                        && matches!(r.kind, RegionKind::Material)
                        && (only.is_empty() || only.contains(&r.name));
                    if !candidate || r.expr == Node::Empty {
                        out.push(r.clone());
                        continue;
                    }
                    let eb = r.expr.bbox(&b.prims);
                    let mut phys = r.clone();
                    phys.expr = Node::Diff(Box::new(r.expr.clone()), pieces.iter().map(|p| p.0.clone()).collect());
                    out.push(phys);
                    for ((node, c, piece), bb) in pieces.iter().zip(&pad_bb) {
                        if eb.overlaps(bb) {
                            out.push(RegionX {
                                name: format!("{}.pml_{piece}", r.name),
                                kind: RegionKind::Pml,
                                feature: r.feature,
                                pinned_id: None,
                                expr: Node::Inter(vec![r.expr.clone(), node.clone()]),
                                parent: Some(r.name.clone()),
                                pml: Some(crate::scene::PmlInfo {
                                    feature: f.id,
                                    sides: *c,
                                    window: win,
                                    thickness: t,
                                    layers: nl,
                                    grading: gr,
                                }),
                            });
                        }
                    }
                }
                rx = out;
            }
            _ => {}
        }
    }
    let mut prims = b.prims;
    for r in &mut rx {
        r.expr = std::mem::replace(&mut r.expr, Node::Empty).simplify(&prims);
    }
    let mut leaves = vec![];
    rx.iter().for_each(|r| r.expr.leaves(&mut leaves));
    let mut map = vec![usize::MAX; prims.len()];
    let mut kept = vec![];
    for (i, p) in prims.drain(..).enumerate() {
        if leaves.contains(&i) {
            map[i] = kept.len();
            kept.push(p);
        }
    }
    rx.iter_mut().for_each(|r| r.expr.remap(&map));
    let bounds = model_box(&kept);
    Ok(Scene {
        dim: 3,
        env,
        prims: kept,
        regions: rx,
        pins: pins.clone(),
        bounds,
        warnings: vec![],
        naming: *naming,
        axis: 2,
        auto_tags,
    })
}

fn map_node(n: &Node, map: &[Option<usize>]) -> Node {
    match n {
        Node::Prim(i) => map[*i].map_or(Node::Empty, Node::Prim),
        Node::Union(v) => Node::Union(v.iter().map(|x| map_node(x, map)).collect()),
        Node::Inter(v) => Node::Inter(v.iter().map(|x| map_node(x, map)).collect()),
        Node::Diff(a, b) => Node::Diff(Box::new(map_node(a, map)), b.iter().map(|x| map_node(x, map)).collect()),
        Node::Empty => Node::Empty,
    }
}

fn hull(mut p: Vec<V2>) -> Vec<V2> {
    p.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    p.dedup_by(|a, b| a == b);
    if p.len() < 3 {
        return p;
    }
    let mut h: Vec<V2> = vec![];
    for pass in 0..2 {
        let start = h.len();
        let it: Vec<V2> = if pass == 0 { p.clone() } else { p.iter().rev().copied().collect() };
        for q in it {
            while h.len() >= start + 2 && cross2(sub2(h[h.len() - 1], h[h.len() - 2]), sub2(q, h[h.len() - 2])) <= 0.0 {
                h.pop();
            }
            h.push(q);
        }
        h.pop();
    }
    h
}

struct Plane {
    n: V3,
    d: f64,
    o: V3,
    u: V3,
    v: V3,
}

/// The region complex of a polyhedral scene.
pub fn build(scene: &Scene, cache: &mut Cache) -> GeomResult<Complex> {
    let diag = scene.bounds.diagonal();
    if diag.is_nan() || diag <= 0.0 || !diag.is_finite() {
        return Err(GeomError("the model has no extent: no region has any volume".into()));
    }
    let tol = 1e-9 * diag;
    let band = crate::arrange::GRAZE_K * tol;
    let c0 = scene.bounds.center();
    let prims = &scene.prims;
    // planes of all faces
    let mut planes: Vec<Plane> = vec![];
    let mut owners: Vec<Vec<(usize, usize, f64)>> = vec![];
    for (i, p) in prims.iter().enumerate() {
        for (j, f) in p.poly.as_ref().unwrap().faces.iter().enumerate() {
            if f.ring.iter().any(|q| (dot3(f.n, *q) - f.d).abs() > 8.0 * tol) {
                return Err(GeomError(format!("a face of '{}' is not planar within the tolerance", p.name)));
            }
            let flip = (0..3).find(|&k| f.n[k].abs() > 1e-9).is_some_and(|k| f.n[k] < 0.0);
            let (n, d, sgn) = if flip { (mul3(f.n, -1.0), -f.d, -1.0) } else { (f.n, f.d, 1.0) };
            let pid =
                planes.iter().position(|q| dot3(q.n, n) > 1.0 - 1e-10 && (q.d - d).abs() <= tol).unwrap_or_else(|| {
                    let a = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    let u = unit3(cross3(n, a));
                    planes.push(Plane { n, d, o: sub3(c0, mul3(n, dot3(n, c0) - d)), u, v: cross3(n, u) });
                    owners.push(vec![]);
                    planes.len() - 1
                });
            owners[pid].push((i, j, sgn));
        }
    }
    let nsr = scene.regions.len();
    let mut fin: Vec<Option<u32>> = vec![None; nsr];
    // generated (Pml) regions that end up with no face anywhere are dropped below
    let mut regions: Vec<Region> = vec![];
    for (i, r) in scene.regions.iter().enumerate() {
        if !matches!(r.kind, RegionKind::Hole | RegionKind::Void) {
            fin[i] = Some(regions.len() as u32);
            regions.push(Region {
                name: r.name.clone(),
                kind: r.kind,
                feature: r.feature,
                pinned_id: r.pinned_id,
                scene_index: i,
                measure: 0.0,
                bbox: Aabb::EMPTY,
                loops: vec![],
                parent: None,
            });
        }
    }
    let leaves: Vec<Vec<usize>> = scene
        .regions
        .iter()
        .map(|r| {
            let mut v = vec![];
            r.expr.leaves(&mut v);
            v
        })
        .collect();
    let h = diag / (2.0 * 2f64.sqrt());
    let plane_bounds = Aabb::new([-h, -h, 0.0], [h, h, 0.0]);

    let mut faces: Vec<Face> = vec![];
    let mut plane_data: Vec<PlaneData> = vec![];
    let mut cell_face: Vec<HashMap<u32, u32>> = vec![]; // per plane: 2D region id -> face id (filled after sorting)
    struct Cand {
        plane: usize,
        cell: u32,
        face: Face,
    }
    let mut cands: Vec<Cand> = vec![];
    for (pi, pl) in planes.iter().enumerate() {
        let to2 = |p: V3| -> V2 {
            let d = sub3(p, pl.o);
            [dot3(d, pl.u), dot3(d, pl.v)]
        };
        let mut pmap: Vec<Option<usize>> = vec![None; prims.len()];
        let mut sprims: Vec<Arc<Prim>> = vec![];
        let mut interior_side: Vec<f64> = vec![]; // 0 for non-owners
        let mut section: Vec<Vec<V2>> = vec![];
        for (i, p) in prims.iter().enumerate() {
            let poly = p.poly.as_ref().unwrap();
            let own = owners[pi].iter().find(|o| o.0 == i);
            let pts: Vec<V2> = if let Some(&(_, j, _)) = own {
                poly.faces[j].ring.iter().map(|q| to2(*q)).collect()
            } else {
                let mut cr: Vec<V2> = vec![];
                for f in &poly.faces {
                    for k in 0..f.ring.len() {
                        let (a, b) = (f.ring[k], f.ring[(k + 1) % f.ring.len()]);
                        let (da, db) = (dot3(pl.n, a) - pl.d, dot3(pl.n, b) - pl.d);
                        if da.abs() <= tol {
                            cr.push(to2(a));
                        } else if (da < 0.0) != (db < 0.0) && db.abs() > tol {
                            cr.push(to2(add3(a, mul3(sub3(b, a), da / (da - db)))));
                        }
                    }
                }
                // merge crossings closer than the undecidable band (a plane through a corner crosses its edges at nearly one point)
                let mut uniq: Vec<V2> = vec![];
                for q in cr {
                    if !uniq.iter().any(|u| dist2(*u, q) <= band) {
                        uniq.push(q);
                    }
                }
                hull(uniq)
            };
            let mut pts = pts;
            if pts.len() < 3 || polygon_area(&pts).abs() < 1e-3 * band * band {
                continue;
            }
            if polygon_area(&pts) < 0.0 {
                pts.reverse();
            }
            let names: Vec<String> = (0..pts.len()).map(|k| format!("e{k}")).collect();
            let prim = Prim::new(p.feature, p.name.clone(), LShape::Poly { pts: pts.clone(), names }, Xf2::ID, p.aux)?;
            pmap[i] = Some(sprims.len());
            sprims.push(Arc::new(prim));
            interior_side.push(own.map_or(0.0, |o| -o.2));
            section.push(pts);
        }
        if sprims.is_empty() {
            plane_data.push(PlaneData { n: pl.n, d: pl.d, o: pl.o, u: pl.u, v: pl.v, cx: empty_complex(tol) });
            cell_face.push(HashMap::new());
            continue;
        }
        let exprs: Vec<Node> = scene.regions.iter().map(|r| map_node(&r.expr, &pmap).simplify(&sprims)).collect();
        let painter = |bit: &dyn Fn(usize) -> bool| -> Option<usize> {
            for (i, r) in scene.regions.iter().enumerate().rev() {
                if exprs[i].eval(bit) {
                    return if matches!(r.kind, RegionKind::Hole | RegionKind::Void) { None } else { Some(i) };
                }
            }
            None
        };
        let classes: RefCell<Vec<(Option<usize>, Option<usize>)>> = RefCell::new(vec![]);
        let classify = |bit: &dyn Fn(usize) -> bool| -> Option<usize> {
            let above = painter(&|j| bit(j) && interior_side[j] >= 0.0);
            let below = painter(&|j| bit(j) && interior_side[j] <= 0.0);
            if above.is_none() && below.is_none() {
                return None;
            }
            let mut c = classes.borrow_mut();
            Some(c.iter().position(|t| *t == (above, below)).unwrap_or_else(|| {
                c.push((above, below));
                c.len() - 1
            }))
        };
        let synth = || -> Vec<RegionX> {
            (0..classes.borrow().len())
                .map(|n| RegionX {
                    name: format!("c{n}"),
                    kind: RegionKind::Material,
                    feature: 0,
                    pinned_id: None,
                    expr: Node::Empty,
                    parent: None,
                    pml: None,
                })
                .collect()
        };
        let ps = Scene {
            dim: 2,
            env: scene.env.clone(),
            prims: sprims.clone(),
            regions: vec![],
            pins: vec![],
            bounds: plane_bounds,
            warnings: vec![],
            naming: NameScheme::Feature,
            axis: 2,
            auto_tags: vec![],
        };
        let cx2 = arrange_with(&ps, cache, &classify, &synth).map_err(|e| {
            GeomError(format!("plane n = ({:.3}, {:.3}, {:.3}), d = {:.6e}: {}", pl.n[0], pl.n[1], pl.n[2], pl.d, e.0))
        })?;
        let tuples = classes.into_inner();
        let lift = |q: V2| add3(pl.o, add3(mul3(pl.u, q[0]), mul3(pl.v, q[1])));
        for (r, reg) in cx2.regions.iter().enumerate() {
            let (above, below) = tuples[reg.scene_index];
            if above == below || reg.loops.is_empty() || reg.measure <= 0.0 {
                continue;
            }
            let pt = crate::stack::interior_point(&cx2, r as u32).unwrap_or([0.0, 0.0]);
            let (first, second) = (below.and_then(|s| fin[s]), above.and_then(|s| fin[s]));
            let name_prov = match (first, second) {
                (Some(x), Some(y)) => {
                    let (mut a, mut b) = (regions[x as usize].name.clone(), regions[y as usize].name.clone());
                    if a > b {
                        std::mem::swap(&mut a, &mut b);
                    }
                    (format!("{a}|{b}"), vec![])
                }
                _ => {
                    let (rs, side) = if second.is_some() { (above.unwrap(), 1.0) } else { (below.unwrap(), -1.0) };
                    let mut best: Option<(u8, u8, u32, String, String)> = None;
                    for &(i, j, sgn) in &owners[pi] {
                        let Some(k) = pmap[i] else { continue };
                        let s = &section[k];
                        let inside =
                            (0..s.len()).all(|e| cross2(sub2(s[(e + 1) % s.len()], s[e]), sub2(pt, s[e])) >= -band);
                        if !inside {
                            continue;
                        }
                        let p = &prims[i];
                        let cand = (
                            u8::from(!leaves[rs].contains(&i)),
                            u8::from(-sgn != side),
                            p.feature,
                            p.name.clone(),
                            p.poly.as_ref().unwrap().faces[j].role.clone(),
                        );
                        if best.as_ref().is_none_or(|b| cand < *b) {
                            best = Some(cand);
                        }
                    }
                    let rname = &scene.regions[rs].name;
                    match best {
                        Some((_, _, feature, inst, role)) => {
                            let nm = match scene.naming {
                                NameScheme::Feature => format!("{inst}.{role}"),
                                _ => format!("{rname}.{role}"),
                            };
                            (nm, vec![Prov { feature, instance: inst, role }])
                        }
                        None => (format!("{rname}.face"), vec![]),
                    }
                }
            };
            let mut bbox = Aabb::EMPTY;
            for &l in &reg.loops {
                for &(e, _) in &cx2.loops[l as usize].edges {
                    let ed = &cx2.edges[e as usize];
                    let c = &cx2.curves[ed.curve as usize];
                    bbox.grow(lift(c.eval(ed.t0)));
                    bbox.grow(lift(c.eval(ed.t1)));
                }
            }
            cands.push(Cand {
                plane: pi,
                cell: r as u32,
                face: Face {
                    name: name_prov.0,
                    geom: FaceGeom::Planar { plane: pi as u32, cell: r as u32 },
                    first,
                    second,
                    prov: name_prov.1,
                    area: reg.measure,
                    bbox,
                },
            });
        }
        plane_data.push(PlaneData { n: pl.n, d: pl.d, o: pl.o, u: pl.u, v: pl.v, cx: cx2 });
        cell_face.push(HashMap::new());
    }
    // deterministic face order: first origin, then plane and cell
    cands.sort_by_key(|c| {
        let k = c.face.prov.iter().map(|p| (p.feature, p.instance.clone(), p.role.clone())).min().unwrap_or((
            u32::MAX,
            c.face.name.clone(),
            String::new(),
        ));
        (k, c.plane, c.cell)
    });
    for (fi, c) in cands.iter().enumerate() {
        cell_face[c.plane].insert(c.cell, fi as u32);
    }
    faces.extend(cands.into_iter().map(|c| c.face));
    if scene.naming == NameScheme::Numbered {
        let mut count = vec![0usize; regions.len()];
        for f in faces.iter_mut().filter(|f| f.first.is_none() || f.second.is_none()) {
            let r = f.first.or(f.second).unwrap() as usize;
            count[r] += 1;
            f.name = format!("{}.{}", regions[r].name, count[r]);
        }
    }
    // volumes by the divergence theorem, bounding boxes from the faces
    for f in &faces {
        let FaceGeom::Planar { plane, .. } = f.geom else { continue };
        let pd = &plane_data[plane as usize];
        let off = pd.d - dot3(pd.n, c0);
        if let Some(r) = f.first {
            regions[r as usize].measure += off * f.area / 3.0;
        }
        if let Some(r) = f.second {
            regions[r as usize].measure -= off * f.area / 3.0;
        }
        for r in [f.first, f.second].into_iter().flatten() {
            regions[r as usize].bbox = regions[r as usize].bbox.union(&f.bbox);
        }
    }
    // edges: cell boundaries merged by 3D line
    struct Seg {
        a: f64,
        b: f64,
        faces: Vec<u32>,
    }
    struct Line {
        o: V3,
        d: V3,
        segs: Vec<Seg>,
    }
    let mut lines: Vec<Line> = vec![];
    for (pi, pd) in plane_data.iter().enumerate() {
        let lift = |q: V2| add3(pd.o, add3(mul3(pd.u, q[0]), mul3(pd.v, q[1])));
        for ed in &pd.cx.edges {
            let fs: Vec<u32> =
                [ed.left, ed.right].into_iter().flatten().filter_map(|c| cell_face[pi].get(&c).copied()).collect();
            if fs.is_empty() {
                continue;
            }
            let c = &pd.cx.curves[ed.curve as usize];
            let (p, q) = (lift(c.eval(ed.t0)), lift(c.eval(ed.t1)));
            let d = unit3(sub3(q, p));
            let li = lines
                .iter()
                .position(|l| {
                    norm3(cross3(l.d, sub3(p, l.o))) <= 2.0 * tol && norm3(cross3(l.d, sub3(q, l.o))) <= 2.0 * tol
                })
                .unwrap_or_else(|| {
                    let flip = (0..3).find(|&k| d[k].abs() > 1e-9).is_some_and(|k| d[k] < 0.0);
                    lines.push(Line { o: p, d: if flip { mul3(d, -1.0) } else { d }, segs: vec![] });
                    lines.len() - 1
                });
            let (sa, sb) = (dot3(lines[li].d, sub3(p, lines[li].o)), dot3(lines[li].d, sub3(q, lines[li].o)));
            lines[li].segs.push(Seg { a: sa.min(sb), b: sa.max(sb), faces: fs });
        }
    }
    let mut edges3: Vec<(u32, usize, f64, Edge3)> = vec![];
    for (li, l) in lines.iter().enumerate() {
        let mut bp: Vec<f64> = l.segs.iter().flat_map(|s| [s.a, s.b]).collect();
        bp.sort_by(f64::total_cmp);
        bp.dedup_by(|a, b| (*a - *b).abs() <= 2.0 * tol);
        let mut pieces: Vec<(f64, f64, Vec<u32>)> = vec![];
        for w in bp.windows(2) {
            let mut fs: Vec<u32> = l
                .segs
                .iter()
                .filter(|s| s.a <= w[0] + 2.0 * tol && s.b >= w[1] - 2.0 * tol)
                .flat_map(|s| s.faces.iter().copied())
                .collect();
            fs.sort_unstable();
            fs.dedup();
            if fs.len() < 2 {
                continue;
            }
            match pieces.last_mut() {
                Some(last) if (last.1 - w[0]).abs() <= 2.0 * tol && last.2 == fs => last.1 = w[1],
                _ => pieces.push((w[0], w[1], fs)),
            }
        }
        for (a, b, fs) in pieces {
            let (pa, pb) = (add3(l.o, mul3(l.d, a)), add3(l.o, mul3(l.d, b)));
            let mut names: Vec<&str> = fs.iter().map(|&f| faces[f as usize].name.as_str()).collect();
            names.sort_unstable();
            let mut bb = Aabb::EMPTY;
            bb.grow(pa);
            bb.grow(pb);
            edges3.push((
                fs[0],
                li,
                a,
                Edge3 {
                    name: names.join("/"),
                    geom: EdgeGeom::Segment { a: pa, b: pb },
                    v0: 0,
                    v1: 0,
                    faces: fs.clone(),
                    length: b - a,
                    bbox: bb,
                },
            ));
        }
    }
    edges3.sort_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)).then(x.2.total_cmp(&y.2)));
    let mut vertices: Vec<Vertex> = vec![];
    let mut edges: Vec<Edge3> = vec![];
    for (_, _, _, mut e) in edges3 {
        let EdgeGeom::Segment { a, b } = e.geom else { unreachable!() };
        let mut ids = [0u32; 2];
        for (k, p) in [a, b].into_iter().enumerate() {
            ids[k] = match vertices.iter().position(|v| norm3(sub3(v.p, p)) <= 2.0 * tol) {
                Some(i) => i as u32,
                None => {
                    vertices.push(Vertex { p, name: String::new(), edges: vec![] });
                    (vertices.len() - 1) as u32
                }
            };
        }
        (e.v0, e.v1) = (ids[0], ids[1]);
        edges.push(e);
    }
    for (ei, e) in edges.iter().enumerate() {
        for v in [e.v0, e.v1] {
            vertices[v as usize].edges.push(ei as u32);
        }
    }
    for v in &mut vertices {
        let mut names: Vec<&str> = v
            .edges
            .iter()
            .flat_map(|&e| edges[e as usize].faces.iter().map(|&f| faces[f as usize].name.as_str()))
            .collect();
        names.sort_unstable();
        names.dedup();
        v.name = names.join("/");
    }
    // drop generated (Pml) regions that have no face at all (painted over or outside every solid)
    let mut used = vec![false; regions.len()];
    for f in &faces {
        for r in [f.first, f.second].into_iter().flatten() {
            used[r as usize] = true;
        }
    }
    let keep: Vec<bool> =
        regions.iter().enumerate().map(|(i, r)| used[i] || scene.regions[r.scene_index].pml.is_none()).collect();
    if keep.iter().any(|k| !k) {
        let mut newid = vec![None; regions.len()];
        let mut kept = vec![];
        for (i, r) in regions.drain(..).enumerate() {
            if keep[i] {
                newid[i] = Some(kept.len() as u32);
                kept.push(r);
            }
        }
        regions = kept;
        for f in &mut faces {
            f.first = f.first.and_then(|r| newid[r as usize]);
            f.second = f.second.and_then(|r| newid[r as usize]);
        }
        for r in &mut regions {
            r.parent = r.parent.and_then(|p| newid[p as usize]);
        }
        for x in fin.iter_mut() {
            *x = x.and_then(|r| newid[r as usize]);
        }
    }
    Ok(Complex {
        dim: 3,
        tol,
        origin: [c0[0], c0[1]],
        curves: vec![],
        regions,
        edges: vec![],
        vertices,
        loops: vec![],
        scene_to_region: fin,
        faces,
        edges3: edges,
        stack: None,
        planar: plane_data,
        axis: 2,
    })
}

fn empty_complex(tol: f64) -> Complex {
    Complex {
        dim: 2,
        tol,
        origin: [0.0; 2],
        curves: vec![],
        regions: vec![],
        edges: vec![],
        vertices: vec![],
        loops: vec![],
        scene_to_region: vec![],
        faces: vec![],
        edges3: vec![],
        stack: None,
        planar: vec![],
        axis: 2,
    }
}
