//! Prismatic 3D models as a *layer stack*: the model is cut at every z level where a solid starts or ends; each layer is a 2D
//! arrangement of the solids active in it, and the *overlay* is the common refinement of all layers (every layer boundary is an
//! edge of it, each overlay cell maps to one region per layer). 3D faces are derived from these: lateral faces are overlay
//! edges swept over runs of layers (merged while the sides, name and edge stay the same), caps are overlay cells at a level where
//! the region below differs from the region above. Names follow the 2D rules: `A|B` for interfaces, `<feature>.<role>` for boundaries.
//!
//! Exactness: volumes, face areas and projections are analytic (a lateral patch is a 2D curve times an interval; a cap is a
//! 2D cell), and the closure of every region (sum of outward area vectors = 0) is a self-check.

use crate::arrange::{arrange, arrange_with};
use crate::complex::*;
use crate::error::{GeomError, GeomResult};
use crate::math::*;
use crate::model::{Axis, NameScheme, RegionKind};
use crate::scene::{Cache, Node, RegionX, Scene};
use std::cell::RefCell;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Layer {
    pub z0: f64,
    pub z1: f64,
    /// the 2D partition of this layer (its own region ids; see `region_ids`)
    pub partition: Complex,
    /// final region id of every region of `partition`
    pub region_ids: Vec<Option<u32>>,
    /// final region id of every overlay cell in this layer (`None`: void)
    pub region_of: Vec<Option<u32>>,
}

#[derive(Clone, Debug)]
pub struct LayerStack {
    pub axis: Axis,
    pub layers: Vec<Layer>,
    /// the common refinement of all layer partitions
    pub(crate) overlay: Complex,
}

impl LayerStack {
    /// Number of cells of the common refinement (a count: the overlay itself lives in the internal sweep frame).
    pub fn overlay_region_count(&self) -> usize {
        self.overlay.regions.len()
    }
    pub fn overlay_edge_count(&self) -> usize {
        self.overlay.edges.len()
    }

    /// Z levels where layers start or end, ascending.
    pub fn levels(&self) -> Vec<f64> {
        let mut v: Vec<f64> = self.layers.iter().flat_map(|l| [l.z0, l.z1]).collect();
        v.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.0));
        v
    }
}

fn substitute(n: &Node, active: &[bool]) -> Node {
    match n {
        Node::Prim(i) if !active[*i] => Node::Empty,
        Node::Union(v) => Node::Union(v.iter().map(|x| substitute(x, active)).collect()),
        Node::Inter(v) => Node::Inter(v.iter().map(|x| substitute(x, active)).collect()),
        Node::Diff(a, b) => {
            Node::Diff(Box::new(substitute(a, active)), b.iter().map(|x| substitute(x, active)).collect())
        }
        o => o.clone(),
    }
}

fn painter(regions: &[RegionX], exprs: &[Node], bit: &dyn Fn(usize) -> bool) -> Option<usize> {
    for (i, r) in regions.iter().enumerate().rev() {
        if exprs[i].eval(bit) {
            return if matches!(r.kind, RegionKind::Hole | RegionKind::Void) { None } else { Some(i) };
        }
    }
    None
}

/// A point strictly inside an overlay cell (halfway along an inward ray from the middle of its first boundary edge).
pub(crate) fn interior_point(cx: &Complex, region: u32) -> Option<V2> {
    let r = &cx.regions[region as usize];
    let l = &cx.loops[*r.loops.iter().find(|&&l| cx.loops[l as usize].area > 0.0)? as usize];
    let &(e, fwd) = l.edges.first()?;
    let ed = &cx.edges[e as usize];
    let t = 0.5 * (ed.t0 + ed.t1);
    let c = &cx.curves[ed.curve as usize];
    let (m, tg) = (c.eval(t), c.tangent(t));
    let tg = if fwd { tg } else { mul2(tg, -1.0) };
    let n = perp(tg); // the region is on the left of the travel direction
    let mut best = f64::INFINITY;
    for &lp in &r.loops {
        for &(e2, _) in &cx.loops[lp as usize].edges {
            let e2d = &cx.edges[e2 as usize];
            let c2 = &cx.curves[e2d.curve as usize];
            let k = if matches!(c2, crate::curve::Curve::Line { .. }) { 1 } else { 128 };
            let mut prev = c2.eval(e2d.t0);
            for i in 1..=k {
                let q = c2.eval(e2d.t0 + (e2d.t1 - e2d.t0) * i as f64 / k as f64);
                // intersect the ray m + s n (s > 0) with the chord prev-q
                let d = sub2(q, prev);
                let den = cross2(n, d);
                if den.abs() > 1e-300 {
                    let w = sub2(prev, m);
                    let (s, u) = (cross2(w, d) / den, cross2(w, n) / den);
                    if s > 1e-12 * (1.0 + best.min(1.0)) && (0.0..=1.0).contains(&u) && s < best {
                        best = s;
                    }
                }
                prev = q;
            }
        }
    }
    if best.is_finite() {
        Some(add2(m, mul2(n, 0.5 * best)))
    } else {
        None
    }
}

pub fn build(scene: &Scene, cache: &mut Cache) -> GeomResult<Complex> {
    let diag = scene.bounds.diagonal();
    if diag.is_nan() || diag <= 0.0 || !diag.is_finite() {
        return Err(GeomError("the model has no extent: no region has any volume".into()));
    }
    let tol = 1e-9 * diag;
    let band = crate::arrange::GRAZE_K * tol;
    let prims = &scene.prims;
    let np = prims.len();

    // z levels
    let mut zs: Vec<(f64, usize)> = prims
        .iter()
        .enumerate()
        .filter(|(_, p)| !p.is_curved_solid())
        .flat_map(|(i, p)| [(p.z.0, i), (p.z.1, i)])
        .collect();
    zs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut levels: Vec<f64> = vec![];
    for (k, &(z, i)) in zs.iter().enumerate() {
        match levels.last() {
            Some(&l) if z - l <= tol => {}
            Some(&l) if z - l <= band => {
                let other = zs[..k].iter().rev().find(|x| (x.0 - l).abs() <= tol).map_or(0, |x| x.1);
                return Err(GeomError(format!(
                    "the z levels of '{}' and '{}' differ by {:.3e} m, close to but not within the tolerance {:.1e} m: make them equal or move them apart",
                    prims[other].name, prims[i].name, z - l, tol
                )));
            }
            _ => levels.push(z),
        }
    }
    let snap = |z: f64| -> f64 { levels.iter().copied().find(|l| (z - l).abs() <= tol).unwrap_or(z) };

    // layers: the 2D partition of the solids active between two levels
    struct Raw {
        z0: f64,
        z1: f64,
        exprs: Vec<Node>,
        cx: Complex,
    }
    let mut raws: Vec<Raw> = vec![];
    for w in levels.windows(2) {
        let (z0, z1) = (w[0], w[1]);
        let active: Vec<bool> =
            prims.iter().map(|p| !p.is_curved_solid() && snap(p.z.0) <= z0 && snap(p.z.1) >= z1).collect();
        let exprs: Vec<Node> = scene.regions.iter().map(|r| substitute(&r.expr, &active).simplify(prims)).collect();
        if exprs.iter().all(|e| *e == Node::Empty) {
            continue;
        }
        // compact the primitives of this layer
        let mut used = vec![false; np];
        let mut lv = vec![];
        exprs.iter().for_each(|e| e.leaves(&mut lv));
        lv.iter().for_each(|&i| used[i] = true);
        let mut map = vec![usize::MAX; np];
        let mut kept: Vec<Arc<crate::prim::Prim>> = vec![];
        for i in 0..np {
            if used[i] {
                map[i] = kept.len();
                kept.push(prims[i].clone());
            }
        }
        let mut regions = scene.regions.clone();
        for (r, e) in regions.iter_mut().zip(&exprs) {
            r.expr = e.clone();
            r.expr.remap(&map);
        }
        let ls = Scene {
            dim: 2,
            env: scene.env.clone(),
            prims: kept,
            regions,
            pins: vec![],
            bounds: scene.bounds,
            warnings: vec![],
            // numbering is applied to the 3D faces at the end; the layers keep role names so fragments can be matched and merged
            naming: if scene.naming == NameScheme::Numbered { NameScheme::Region } else { scene.naming },
            axis: 2,
            auto_tags: vec![],
        };
        let cx = arrange(&ls, cache).map_err(|e| GeomError(format!("layer z = {z0:.6e}..{z1:.6e}: {}", e.0)))?;
        raws.push(Raw { z0, z1, exprs, cx });
    }
    if raws.is_empty() {
        return Err(GeomError("the model has no prismatic host: curved bodies (sphere, torus) must sit inside a box, cylinder or extruded solid in this build".into()));
    }

    // final regions: every declared region, and generated (Pml) ones that exist in some layer
    let nsr = scene.regions.len();
    let mut present = vec![false; nsr];
    for (i, r) in scene.regions.iter().enumerate() {
        if matches!(r.kind, RegionKind::Hole | RegionKind::Void) {
            continue;
        }
        present[i] = r.pml.is_none() || raws.iter().any(|l| l.cx.regions.iter().any(|x| x.scene_index == i));
    }
    let mut fin: Vec<Option<u32>> = vec![None; nsr];
    let mut regions: Vec<Region> = vec![];
    for (i, r) in scene.regions.iter().enumerate() {
        if present[i] {
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
    for i in 0..regions.len() {
        if let Some(pn) = &scene.regions[regions[i].scene_index].parent {
            regions[i].parent = regions.iter().position(|r| &r.name == pn).map(|p| p as u32);
        }
    }
    for l in &raws {
        for r in &l.cx.regions {
            if let Some(g) = fin[r.scene_index] {
                let g = &mut regions[g as usize];
                g.measure += r.measure * (l.z1 - l.z0);
                let mut b = r.bbox;
                b.lo[2] = l.z0;
                b.hi[2] = l.z1;
                g.bbox = g.bbox.union(&b);
            }
        }
    }

    // overlay: classes are tuples of the per-layer regions
    let classes: RefCell<Vec<Vec<Option<u32>>>> = RefCell::new(vec![]);
    let classify = |bit: &dyn Fn(usize) -> bool| -> Option<usize> {
        let tuple: Vec<Option<u32>> =
            raws.iter().map(|l| painter(&scene.regions, &l.exprs, bit).and_then(|s| fin[s])).collect();
        if tuple.iter().all(|t| t.is_none()) {
            return None;
        }
        let mut c = classes.borrow_mut();
        Some(c.iter().position(|t| *t == tuple).unwrap_or_else(|| {
            c.push(tuple);
            c.len() - 1
        }))
    };
    let synth = || -> Vec<RegionX> {
        (0..classes.borrow().len())
            .map(|n| RegionX {
                name: format!("cell{n}"),
                kind: RegionKind::Material,
                feature: 0,
                pinned_id: None,
                expr: Node::Empty,
                parent: None,
                pml: None,
            })
            .collect()
    };
    let ov_scene = Scene {
        dim: 2,
        env: scene.env.clone(),
        prims: prims.clone(),
        regions: vec![],
        pins: vec![],
        bounds: scene.bounds,
        warnings: vec![],
        naming: NameScheme::Feature,
        axis: 2,
        auto_tags: vec![],
    };
    let mut overlay = arrange_with(&ov_scene, cache, &classify, &synth)?;
    let tuples = classes.into_inner();
    let ncell = overlay.regions.len();
    let cell_tuple: Vec<&Vec<Option<u32>>> = overlay.regions.iter().map(|r| &tuples[r.scene_index]).collect();
    for r in &mut overlay.regions {
        r.name = format!("cell{}", r.scene_index);
    }

    // per-layer tables
    let layers: Vec<Layer> = raws
        .iter()
        .enumerate()
        .map(|(k, l)| Layer {
            z0: l.z0,
            z1: l.z1,
            partition: l.cx.clone(),
            region_ids: l.cx.regions.iter().map(|r| fin[r.scene_index]).collect(),
            region_of: (0..ncell).map(|c| cell_tuple[c][k]).collect(),
        })
        .collect();

    // lateral faces: overlay edges against the edges of each layer (which carry the names)
    struct Lat {
        edge: u32,
        k0: usize,
        k1: usize,
        name: String,
        first: Option<u32>,
        second: Option<u32>,
        prov: Vec<Prov>,
    }
    let mut lats: Vec<Lat> = vec![];
    for (ei, oe) in overlay.edges.iter().enumerate() {
        let c = &overlay.curves[oe.curve as usize];
        let mid = c.eval(0.5 * (oe.t0 + oe.t1));
        let tg = c.tangent(0.5 * (oe.t0 + oe.t1));
        for (k, l) in layers.iter().enumerate() {
            let mut best: Option<(f64, u32)> = None;
            for (li, le) in l.partition.edges.iter().enumerate() {
                let bb = le.bbox.pad(8.0 * tol);
                if mid[0] < bb.lo[0] || mid[0] > bb.hi[0] || mid[1] < bb.lo[1] || mid[1] > bb.hi[1] {
                    continue;
                }
                let d = dist2(l.partition.project_edge(li as u32, mid), mid);
                if d <= 8.0 * tol && best.is_none_or(|b| d < b.0) {
                    best = Some((d, li as u32));
                }
            }
            let Some((_, li)) = best else { continue };
            let le = &l.partition.edges[li as usize];
            let lc = &l.partition.curves[le.curve as usize];
            let ltg = lc.tangent(lc.closest(mid).0);
            let (a, b) =
                (le.left.and_then(|r| l.region_ids[r as usize]), le.right.and_then(|r| l.region_ids[r as usize]));
            let (first, second) = if dot2(ltg, tg) >= 0.0 { (a, b) } else { (b, a) };
            if let Some(prev) = lats.last_mut() {
                if prev.edge == ei as u32
                    && prev.k1 + 1 == k
                    && prev.first == first
                    && prev.second == second
                    && prev.name == le.name
                    && (layers[k - 1].z1 - l.z0).abs() <= tol
                {
                    prev.k1 = k;
                    continue;
                }
            }
            lats.push(Lat {
                edge: ei as u32,
                k0: k,
                k1: k,
                name: le.name.clone(),
                first,
                second,
                prov: le.prov.clone(),
            });
        }
    }
    let mut faces: Vec<Face> = vec![];
    for f in lats {
        let oe = &overlay.edges[f.edge as usize];
        let (z0, z1) = (layers[f.k0].z0, layers[f.k1].z1);
        let mut bbox = oe.bbox;
        bbox.lo[2] = z0;
        bbox.hi[2] = z1;
        faces.push(Face {
            name: f.name,
            geom: FaceGeom::Lateral { edge: f.edge, z0, z1 },
            first: f.first,
            second: f.second,
            prov: f.prov,
            area: oe.length * (z1 - z0),
            bbox,
        });
    }

    // caps: at every level, cells whose region below differs from the region above
    let mut cell_pt: Vec<Option<V2>> = vec![];
    for c in 0..ncell {
        cell_pt.push(interior_point(&overlay, c as u32));
    }
    let nl = layers.len();
    for k in 0..=nl {
        let below = if k > 0 { Some(k - 1) } else { None };
        let above = if k < nl { Some(k) } else { None };
        // a level between two layers that do not touch is two levels
        let touching = matches!((below, above), (Some(b), Some(a)) if (layers[b].z1 - layers[a].z0).abs() <= tol);
        let mut jobs: Vec<(Option<usize>, Option<usize>, f64)> = vec![];
        match (below, above) {
            (Some(b), Some(a)) if !touching => {
                jobs.push((Some(b), None, layers[b].z1));
                jobs.push((None, Some(a), layers[a].z0));
            }
            (Some(b), Some(a)) => jobs.push((Some(b), Some(a), layers[b].z1.max(layers[a].z0))),
            (Some(b), None) => jobs.push((Some(b), None, layers[b].z1)),
            (None, Some(a)) => jobs.push((None, Some(a), layers[a].z0)),
            (None, None) => {}
        }
        for (b, a, z) in jobs {
            #[allow(clippy::needless_range_loop)]
            for c in 0..ncell {
                let (rb, ra) = (b.and_then(|b| layers[b].region_of[c]), a.and_then(|a| layers[a].region_of[c]));
                if rb == ra {
                    continue;
                }
                let name_prov = match (rb, ra) {
                    (Some(x), Some(y)) => {
                        let (mut p, mut q) = (regions[x as usize].name.clone(), regions[y as usize].name.clone());
                        if p > q {
                            std::mem::swap(&mut p, &mut q);
                        }
                        (format!("{p}|{q}"), vec![])
                    }
                    _ => {
                        // the solid that ends (top cap) or starts (bottom cap) here and contains the cell
                        let (layer, top) = if rb.is_some() { (b.unwrap(), true) } else { (a.unwrap(), false) };
                        let region = rb.or(ra).unwrap();
                        let sidx = regions[region as usize].scene_index;
                        let mut lv = vec![];
                        scene.regions[sidx].expr.leaves(&mut lv);
                        let pt = cell_pt[c].unwrap_or([0.0, 0.0]);
                        let (lz0, lz1) = (layers[layer].z0, layers[layer].z1);
                        let mut best: Option<(u8, u32, String)> = None;
                        for (i, p) in prims.iter().enumerate() {
                            let ends = if top { (p.z.1 - lz1).abs() <= tol } else { (p.z.0 - lz0).abs() <= tol };
                            if !ends || p.z.0 > lz0 + tol || p.z.1 < lz1 - tol || p.sdf_near(pt, band) >= 0.0 {
                                continue;
                            }
                            let cand = (
                                u8::from(!lv.contains(&i)),
                                p.feature,
                                format!("{}.{}", p.name, p.caps[usize::from(top)]),
                            );
                            if best.as_ref().is_none_or(|x| cand < *x) {
                                best = Some(cand);
                            }
                        }
                        match best {
                            Some((_, feature, nm)) => {
                                let (inst, role) = nm.rsplit_once('.').unwrap();
                                let shown = if scene.naming == NameScheme::Feature {
                                    nm.clone()
                                } else {
                                    format!("{}.{role}", regions[region as usize].name)
                                };
                                (shown, vec![Prov { feature, instance: inst.to_string(), role: role.to_string() }])
                            }
                            None => (
                                format!(
                                    "{}.{}",
                                    regions[region as usize].name,
                                    crate::scene::cap_role(scene.axis, top)
                                ),
                                vec![],
                            ),
                        }
                    }
                };
                let cell = &overlay.regions[c];
                let mut bbox = cell.bbox;
                bbox.lo[2] = z;
                bbox.hi[2] = z;
                faces.push(Face {
                    name: name_prov.0,
                    geom: FaceGeom::Cap { cell: c as u32, z },
                    first: rb,
                    second: ra,
                    prov: name_prov.1,
                    area: cell.measure,
                    bbox,
                });
            }
        }
    }

    // deterministic order: first origin, then lateral before caps, then edge / cell and height
    let fkey = |f: &Face| -> (u32, String, String, u8, u32, i64) {
        let (fe, ins, ro) = f.prov.iter().map(|p| (p.feature, p.instance.clone(), p.role.clone())).min().unwrap_or((
            u32::MAX,
            f.name.clone(),
            String::new(),
        ));
        match f.geom {
            FaceGeom::Lateral { edge, z0, .. } => (fe, ins, ro, 0, edge, (z0 / tol).round() as i64),
            FaceGeom::Cap { cell, z } => (fe, ins, ro, 1, cell, (z / tol).round() as i64),
            FaceGeom::Quadric(_) | FaceGeom::Planar { .. } => (fe, ins, ro, 2, 0, 0),
        }
    };
    let mut keyed: Vec<_> = faces.into_iter().map(|f| (fkey(&f), f)).collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0));
    let mut faces: Vec<Face> = keyed.into_iter().map(|x| x.1).collect();
    if scene.naming == NameScheme::Numbered {
        let mut count = vec![0usize; regions.len()];
        for f in faces.iter_mut().filter(|f| f.first.is_none() || f.second.is_none()) {
            let r = f.first.or(f.second).unwrap() as usize;
            count[r] += 1;
            f.name = format!("{}.{}", regions[r].name, count[r]);
        }
    }

    let stack = LayerStack { axis: [Axis::X, Axis::Y, Axis::Z][scene.axis], layers, overlay: overlay.clone() };
    let cx = Complex {
        dim: 3,
        tol,
        origin: [scene.bounds.center()[0], scene.bounds.center()[1]],
        curves: overlay.curves.clone(),
        regions,
        edges: vec![],
        vertices: vec![],
        loops: vec![],
        scene_to_region: fin,
        faces,
        edges3: vec![],
        planar: vec![],
        axis: 2,
        stack: Some(Box::new(stack)),
    };
    let mut cx = cx;
    let circles = add_bodies(scene, &mut cx, band)?;
    let (edges3, vertices) = crate::topo3::build(&cx);
    cx.edges3 = edges3;
    cx.vertices = vertices;
    for (c, r, fs) in circles {
        // a circle edge has one seam vertex at angle 0
        let seam = add3(c, [r, 0.0, 0.0]);
        cx.vertices.push(Vertex { p: seam, name: String::new(), edges: vec![cx.edges3.len() as u32] });
        let v = (cx.vertices.len() - 1) as u32;
        let mut names: Vec<String> = fs.iter().map(|&f| cx.faces[f as usize].name.clone()).collect();
        names.sort();
        let mut bb = Aabb::EMPTY;
        bb.grow(sub3(c, [r, r, 0.0]));
        bb.grow(add3(c, [r, r, 0.0]));
        cx.vertices[v as usize].name = names.join("/");
        cx.edges3.push(Edge3 {
            name: names.join("/"),
            geom: EdgeGeom::Circle { c, r },
            v0: v,
            v1: v,
            faces: fs,
            length: TAU * r,
            bbox: bb,
        });
    }
    cx.set_world_frame(scene.axis);
    Ok(cx)
}

/// An analytic surface that bounds faces (what an isosurface-stuffing mesher intersects its lattice with).
#[derive(Clone, Debug, PartialEq)]
pub enum Surface {
    /// `n . p = d`, `n` a unit normal
    Plane {
        n: V3,
        d: f64,
    },
    /// the axis through `c` along `axis` (unit), radius `r`
    Cylinder {
        c: V3,
        axis: V3,
        r: f64,
    },
    Sphere {
        c: V3,
        r: f64,
    },
    Torus {
        c: V3,
        major: f64,
        minor: f64,
    },
}

/// One surface with the faces on it and the extent of those faces, so a mesher keeps a *local* surface list per lattice cell.
#[derive(Clone, Debug)]
pub struct StuffSurface {
    pub surface: Surface,
    pub face_aabb: Aabb,
    pub faces: Vec<u32>,
}

impl crate::eval::Evaluated {
    /// Surfaces (planes and cylinders; coincident faces share one) with the face extents, in a deterministic order.
    pub fn surfaces_for_stuffing(&self) -> Vec<StuffSurface> {
        let cx = self.complex();
        let Some(st) = &cx.stack else { return vec![] };
        let tol = cx.tol;
        let mut out: Vec<StuffSurface> = vec![];
        for (i, f) in cx.faces.iter().enumerate() {
            let surf = match f.geom {
                FaceGeom::Lateral { edge, .. } => {
                    match &st.overlay.curves[st.overlay.edges[edge as usize].curve as usize] {
                        crate::curve::Curve::Line { o, d } => {
                            let n = [d[1], -d[0], 0.0];
                            Surface::Plane { n: cx.wv(n), d: n[0] * o[0] + n[1] * o[1] }
                        }
                        crate::curve::Curve::Circle { c, r } => {
                            Surface::Cylinder { c: cx.wv([c[0], c[1], 0.0]), axis: cx.wv([0.0, 0.0, 1.0]), r: *r }
                        }
                        _ => continue,
                    }
                }
                FaceGeom::Cap { z, .. } => Surface::Plane { n: cx.wv([0.0, 0.0, 1.0]), d: z },
                FaceGeom::Planar { plane, .. } => {
                    let pd = &cx.planar[plane as usize];
                    Surface::Plane { n: pd.n, d: pd.d }
                }
                FaceGeom::Quadric(crate::prim::Body::CylSide { c, r, .. }) => {
                    Surface::Cylinder { c: cx.wv([c[0], c[1], 0.0]), axis: cx.wv([0.0, 0.0, 1.0]), r }
                }
                FaceGeom::Quadric(crate::prim::Body::Annulus { c, up, .. }) => Surface::Plane {
                    n: cx.wv([0.0, 0.0, if up { 1.0 } else { -1.0 }]),
                    d: c[2] * if up { 1.0 } else { -1.0 },
                },
                FaceGeom::Quadric(crate::prim::Body::Disk { c, up, .. }) => Surface::Plane {
                    n: cx.wv([0.0, 0.0, if up { 1.0 } else { -1.0 }]),
                    d: c[2] * if up { 1.0 } else { -1.0 },
                },
                FaceGeom::Quadric(crate::prim::Body::Sphere { c, r }) => Surface::Sphere { c: cx.wv(c), r },
                FaceGeom::Quadric(crate::prim::Body::Cone { .. } | crate::prim::Body::ConeSide { .. }) => continue,
                FaceGeom::Quadric(crate::prim::Body::Torus { c, big, small }) => {
                    Surface::Torus { c, major: big, minor: small }
                }
            };
            let same = |a: &Surface| match (a, &surf) {
                (Surface::Plane { n: n1, d: d1 }, Surface::Plane { n: n2, d: d2 }) => {
                    norm3(sub3(*n1, *n2)) < 1e-12 && (d1 - d2).abs() <= tol
                }
                (Surface::Cylinder { c: c1, r: r1, .. }, Surface::Cylinder { c: c2, r: r2, .. }) => {
                    norm3(sub3(*c1, *c2)) <= tol && (r1 - r2).abs() <= tol
                }
                _ => false,
            };
            match out.iter_mut().find(|s| same(&s.surface)) {
                Some(s) => {
                    s.face_aabb = s.face_aabb.union(&f.bbox);
                    s.faces.push(i as u32);
                }
                None => out.push(StuffSurface { surface: surf, face_aabb: f.bbox, faces: vec![i as u32] }),
            }
        }
        out
    }
}

/// Curved bodies (sphere, torus) are added to the host stack as closed faces without edges, provided every body is clear of every
/// other face and body (disjoint or nested spheres): the distance is exact for spheres, a bounding-box test for tori. Anything that
/// touches or crosses is refused with both names; the intersection curves of general quadrics are not implemented.
/// Circle edges of the bodies added: centre, radius and the faces that meet there.
type BodyCircles = Vec<(V3, f64, Vec<u32>)>;

/// Faces `(role, geometry)` of a curved solid and its rim circles `(centre, radius, indices of the two faces)`.
type SolidParts = (Vec<(String, crate::prim::Body)>, Vec<(V3, f64, [usize; 2])>);

/// A curved solid of the scene: a single-surface body or a solid of revolution.
enum Solid<'a> {
    B(crate::prim::Body),
    R(&'a crate::prim::Revolved),
}

impl Solid<'_> {
    fn kind(&self) -> &'static str {
        match self {
            Solid::B(b) => b.kind(),
            Solid::R(_) => "revolve",
        }
    }
    fn bbox(&self) -> Aabb {
        match self {
            Solid::B(b) => b.bbox(),
            Solid::R(r) => r.bbox(),
        }
    }
    fn volume(&self) -> f64 {
        match self {
            Solid::B(b) => b.volume(),
            Solid::R(r) => r.volume(),
        }
    }
    fn surface_point(&self) -> V3 {
        use crate::prim::Body;
        match self {
            Solid::B(Body::Sphere { c, r }) => [c[0] + r, c[1], c[2]],
            Solid::B(Body::Torus { c, big, small }) => [c[0] + big + small, c[1], c[2]],
            Solid::B(Body::Cone { c, r0, r1, h }) => [c[0] + 0.5 * (r0 + r1), c[1], c[2] + 0.5 * h],
            Solid::B(_) => unreachable!(),
            Solid::R(r) => r.surface_point(),
        }
    }
    /// The faces (role, geometry) and the circle edges as pairs of face indices into that list with centre and radius.
    fn parts(&self) -> SolidParts {
        use crate::prim::Body;
        match self {
            Solid::B(Body::Cone { c, r0, r1, h }) => {
                let (c, r0, r1, h) = (*c, *r0, *r1, *h);
                let mut v = vec![("side".to_string(), Body::ConeSide { c, r0, r1, h, flip: false })];
                let mut circles = vec![];
                if r0 > 0.0 {
                    v.push(("bottom".to_string(), Body::Disk { c, r: r0, up: false }));
                    circles.push((c, r0, [0, v.len() - 1]));
                }
                if r1 > 0.0 {
                    let top = [c[0], c[1], c[2] + h];
                    v.push(("top".to_string(), Body::Disk { c: top, r: r1, up: true }));
                    circles.push((top, r1, [0, v.len() - 1]));
                }
                (v, circles)
            }
            Solid::B(b) => (vec![("surface".to_string(), *b)], vec![]),
            Solid::R(r) => {
                let faces = r.faces();
                let n = r.pts.len();
                let index_of = |seg: usize| faces.iter().position(|f| f.seg == seg);
                let mut circles = vec![];
                for i in 0..n {
                    // the corner pts[i] joins segment i-1 and segment i
                    if r.pts[i][0] > 0.0 {
                        if let (Some(a), Some(b)) = (index_of((i + n - 1) % n), index_of(i)) {
                            circles.push(([r.c[0], r.c[1], r.c[2] + r.pts[i][1]], r.pts[i][0], [a, b]));
                        }
                    }
                }
                (faces.into_iter().map(|f| (f.role, f.body)).collect(), circles)
            }
        }
    }
}

fn add_bodies(scene: &Scene, cx: &mut Complex, band: f64) -> GeomResult<BodyCircles> {
    use crate::prim::Body;
    let prims = &scene.prims;
    let solids: Vec<(usize, Solid)> = prims
        .iter()
        .enumerate()
        .filter_map(|(i, p)| match (&p.body, &p.rev) {
            (Some(b), _) => Some((i, Solid::B(*b))),
            (_, Some(r)) => Some((i, Solid::R(r))),
            _ => None,
        })
        .collect();
    let mut circles: BodyCircles = vec![];
    if solids.is_empty() {
        return Ok(circles);
    }
    let nhost = cx.faces.len();
    for (a, (ia, sa)) in solids.iter().enumerate() {
        let name = &prims[*ia].name;
        for fi in 0..nhost {
            let fname = &cx.faces[fi].name;
            let clear = match sa {
                Solid::B(Body::Sphere { c, r }) => norm3(sub3(cx.facet_project(fi as u32, *c), *c)) > r + band,
                _ => !sa.bbox().pad(band).overlaps(&cx.faces[fi].bbox),
            };
            if !clear {
                return Err(GeomError(format!(
                    "{} '{name}' touches, crosses or encloses the face '{fname}'; a curved body must lie clear of every other face (its intersection curves are not implemented): move it by more than {band:.1e} m",
                    sa.kind()
                )));
            }
        }
        for (ib, sb) in &solids[a + 1..] {
            let other = &prims[*ib].name;
            let ok = match (sa, sb) {
                (Solid::B(Body::Sphere { c: c1, r: r1 }), Solid::B(Body::Sphere { c: c2, r: r2 })) => {
                    let d = norm3(sub3(*c1, *c2));
                    d > r1 + r2 + band || d + r1.min(*r2) < r1.max(*r2) - band
                }
                _ => !sa.bbox().pad(band).overlaps(&sb.bbox()),
            };
            if !ok {
                return Err(GeomError(format!(
                    "{} '{name}' and {} '{other}' intersect or touch; only disjoint bodies, or one sphere strictly inside another, are supported",
                    sa.kind(),
                    sb.kind()
                )));
            }
        }
    }
    for (ia, sa) in &solids {
        let p0 = sa.surface_point();
        let region_at = |inside: bool| -> Option<u32> {
            let bits: Vec<bool> =
                prims.iter().enumerate().map(|(i, q)| if i == *ia { inside } else { q.sdf3(p0, band) < 0.0 }).collect();
            scene.region_of(&|i| bits[i]).and_then(|s| cx.scene_to_region[s])
        };
        let (rin, rout) = (region_at(true), region_at(false));
        if rin == rout {
            continue;
        }
        let v = sa.volume();
        if let Some(r) = rin {
            cx.regions[r as usize].measure += v;
            cx.regions[r as usize].bbox = cx.regions[r as usize].bbox.union(&sa.bbox());
        }
        if let Some(r) = rout {
            cx.regions[r as usize].measure -= v;
        }
        let p = &prims[*ia];
        let interface = match (rin, rout) {
            (Some(x), Some(y)) => {
                let (mut a, mut b) = (cx.regions[x as usize].name.clone(), cx.regions[y as usize].name.clone());
                if a > b {
                    std::mem::swap(&mut a, &mut b);
                }
                Some(format!("{a}|{b}"))
            }
            _ => None,
        };
        let (mut parts, rims) = sa.parts();
        for (role, _) in &mut parts {
            // the end disks of a cone are the caps of its sweep axis
            match role.as_str() {
                "bottom" => *role = crate::scene::cap_role(scene.axis, false),
                "top" => *role = crate::scene::cap_role(scene.axis, true),
                _ => {}
            }
        }
        let first_id = cx.faces.len() as u32;
        for (role, part) in &parts {
            let name = interface.clone().unwrap_or_else(|| match scene.naming {
                crate::model::NameScheme::Feature => format!("{}.{role}", p.name),
                _ => format!("{}.{role}", cx.regions[rin.or(rout).unwrap() as usize].name),
            });
            cx.faces.push(Face {
                name,
                geom: FaceGeom::Quadric(*part),
                first: rin,
                second: rout,
                prov: vec![Prov { feature: p.feature, instance: p.name.clone(), role: role.clone() }],
                area: part.area(),
                bbox: part.bbox(),
            });
        }
        for (c, r, [i, j]) in rims {
            circles.push((c, r, vec![first_id + i as u32, first_id + j as u32]));
        }
    }
    if scene.naming == NameScheme::Numbered {
        let mut count = vec![0usize; cx.regions.len()];
        for f in cx.faces.iter_mut().filter(|f| f.first.is_none() || f.second.is_none()) {
            let r = f.first.or(f.second).unwrap() as usize;
            count[r] += 1;
            f.name = format!("{}.{}", cx.regions[r].name, count[r]);
        }
    }
    Ok(circles)
}
