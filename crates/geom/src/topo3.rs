//! Edges and vertices of a 3D (layer stack) model, derived from its face fragments so that the topology is closed and checkable.
//!
//! Rules (all on fragments, matching how faces are defined):
//! * a **vertical edge** runs along an overlay vertex through a run of layers when at least two distinct lateral faces meet there
//!   (two at a corner, three where three regions meet); consecutive layers with the same set of faces are one edge;
//! * a **horizontal edge** is an overlay edge at a level where at least two distinct faces meet it (the lateral faces below and
//!   above, and the caps of the cells on its two sides);
//! * a **vertex** is an overlay vertex at a level, created where an edge ends.
//!
//! Names: an edge is named by the sorted names of the faces that meet at it, joined by `/` (`box.xmin/box.ymin`); a vertex by the sorted
//! unique names of the faces around it. Ids are assigned in a deterministic order (first face, kind, overlay id, level).

use crate::complex::*;
use crate::math::*;
use std::collections::BTreeMap;

pub fn build(cx: &Complex) -> (Vec<Edge3>, Vec<Vertex>) {
    let st = cx.stack.as_ref().expect("a layer stack");
    let (ov, layers, tol) = (&st.overlay, &st.layers, cx.tol);
    // levels (deduplicated within tol) and the layer that ends / starts at each
    let mut levels: Vec<f64> = layers.iter().flat_map(|l| [l.z0, l.z1]).collect();
    levels.sort_by(f64::total_cmp);
    levels.dedup_by(|a, b| (*a - *b).abs() <= tol);
    let lev = |z: f64| levels.iter().position(|l| (l - z).abs() <= tol).unwrap_or(0);
    let below: Vec<Option<usize>> =
        levels.iter().map(|&z| layers.iter().position(|l| (l.z1 - z).abs() <= tol)).collect();
    let above: Vec<Option<usize>> =
        levels.iter().map(|&z| layers.iter().position(|l| (l.z0 - z).abs() <= tol)).collect();
    // lateral face of (overlay edge, layer) and cap of (cell, level)
    let mut lat: BTreeMap<(u32, usize), u32> = BTreeMap::new();
    let mut cap: BTreeMap<(u32, usize), u32> = BTreeMap::new();
    for (fi, f) in cx.faces.iter().enumerate() {
        match f.geom {
            FaceGeom::Lateral { edge, z0, z1 } => {
                for (k, l) in layers.iter().enumerate() {
                    if l.z0 >= z0 - tol && l.z1 <= z1 + tol {
                        lat.insert((edge, k), fi as u32);
                    }
                }
            }
            FaceGeom::Cap { cell, z } => {
                cap.insert((cell, lev(z)), fi as u32);
            }
            FaceGeom::Quadric(_) | FaceGeom::Planar { .. } => {}
        }
    }
    struct Raw {
        key: (u32, u8, u32, usize),
        geom: EdgeGeom,
        ends: [(u32, usize); 2],
        faces: Vec<u32>,
    }
    let mut raws: Vec<Raw> = vec![];
    // vertical edges: per overlay vertex, runs of layers with the same face set
    for (v, vx) in ov.vertices.iter().enumerate() {
        let mut run: Option<(usize, usize, Vec<u32>)> = None;
        let flush = |run: &mut Option<(usize, usize, Vec<u32>)>, raws: &mut Vec<Raw>| {
            if let Some((k0, k1, faces)) = run.take() {
                let (z0, z1) = (layers[k0].z0, layers[k1].z1);
                raws.push(Raw {
                    key: (faces[0], 1, v as u32, lev(z0)),
                    geom: EdgeGeom::Vertical { vertex: v as u32, z0, z1 },
                    ends: [(v as u32, lev(z0)), (v as u32, lev(z1))],
                    faces,
                });
            }
        };
        for k in 0..layers.len() {
            let mut f: Vec<u32> = vx.edges.iter().filter_map(|&e| lat.get(&(e, k)).copied()).collect();
            f.sort_unstable();
            f.dedup();
            let touching = k > 0 && (layers[k - 1].z1 - layers[k].z0).abs() <= tol;
            match run.as_mut() {
                Some(r) if touching && r.1 + 1 == k && r.2 == f => r.1 = k,
                _ => {
                    flush(&mut run, &mut raws);
                    if f.len() >= 2 {
                        run = Some((k, k, f));
                    }
                }
            }
        }
        flush(&mut run, &mut raws);
    }
    // horizontal edges: per overlay edge and level
    for (e, ed) in ov.edges.iter().enumerate() {
        for (li, &z) in levels.iter().enumerate() {
            let mut f: Vec<u32> = vec![];
            for k in [below[li], above[li]].into_iter().flatten() {
                f.extend(lat.get(&(e as u32, k)));
            }
            for c in [ed.left, ed.right].into_iter().flatten() {
                f.extend(cap.get(&(c, li)));
            }
            f.sort_unstable();
            f.dedup();
            if f.len() >= 2 {
                raws.push(Raw {
                    key: (f[0], 0, e as u32, li),
                    geom: EdgeGeom::Horizontal { edge: e as u32, z },
                    ends: [(ed.v0, li), (ed.v1, li)],
                    faces: f,
                });
            }
        }
    }
    raws.sort_by_key(|a| a.key);
    // vertices in order of first appearance
    let mut vid: BTreeMap<(u32, usize), u32> = BTreeMap::new();
    let mut vertices: Vec<Vertex> = vec![];
    let mut edges: Vec<Edge3> = vec![];
    for r in raws {
        let mut ids = [0u32; 2];
        for (i, &(v, li)) in r.ends.iter().enumerate() {
            ids[i] = *vid.entry((v, li)).or_insert_with(|| {
                let p = ov.vertices[v as usize].p;
                vertices.push(Vertex { p: [p[0], p[1], levels[li]], name: String::new(), edges: vec![] });
                (vertices.len() - 1) as u32
            });
        }
        let mut names: Vec<&str> = r.faces.iter().map(|&f| cx.faces[f as usize].name.as_str()).collect();
        names.sort_unstable();
        let (length, bbox) = match r.geom {
            EdgeGeom::Vertical { vertex, z0, z1 } => {
                let p = ov.vertices[vertex as usize].p;
                (z1 - z0, Aabb::new([p[0], p[1], z0], [p[0], p[1], z1]))
            }
            EdgeGeom::Segment { .. } | EdgeGeom::Circle { .. } => unreachable!(),
            EdgeGeom::Horizontal { edge, z } => {
                let e = &ov.edges[edge as usize];
                let mut b = e.bbox;
                b.lo[2] = z;
                b.hi[2] = z;
                (e.length, b)
            }
        };
        edges.push(Edge3 { name: names.join("/"), geom: r.geom, v0: ids[0], v1: ids[1], faces: r.faces, length, bbox });
    }
    for (ei, e) in edges.iter().enumerate() {
        for v in [e.v0, e.v1] {
            if vertices[v as usize].edges.last() != Some(&(ei as u32)) {
                vertices[v as usize].edges.push(ei as u32);
            }
        }
    }
    for v in &mut vertices {
        let mut names: Vec<&str> = v
            .edges
            .iter()
            .flat_map(|&e| edges[e as usize].faces.iter().map(|&f| cx.faces[f as usize].name.as_str()))
            .collect();
        names.sort_unstable();
        names.dedup();
        v.name = names.join("/");
    }
    (edges, vertices)
}

impl Complex {
    /// Point of a 3D edge at fraction `u` of its parameter range.
    pub fn edge3_point(&self, e: u32, u: f64) -> V3 {
        self.wv(self.edge3_point_internal(e, u))
    }
    fn edge3_point_internal(&self, e: u32, u: f64) -> V3 {
        let u = u.clamp(0.0, 1.0);
        if let EdgeGeom::Segment { a, b } = self.edges3[e as usize].geom {
            return add3(a, mul3(sub3(b, a), u));
        }
        if let EdgeGeom::Circle { c, r } = self.edges3[e as usize].geom {
            return add3(c, [r * (u * TAU).cos(), r * (u * TAU).sin(), 0.0]);
        }
        let ov = &self.stack.as_ref().expect("a layer stack").overlay;
        match self.edges3[e as usize].geom {
            EdgeGeom::Segment { .. } | EdgeGeom::Circle { .. } => unreachable!(),
            EdgeGeom::Vertical { vertex, z0, z1 } => {
                let p = ov.vertices[vertex as usize].p;
                [p[0], p[1], z0 + u * (z1 - z0)]
            }
            EdgeGeom::Horizontal { edge, z } => {
                let ed = &ov.edges[edge as usize];
                let q = ov.curves[ed.curve as usize].eval(ed.t0 + u * (ed.t1 - ed.t0));
                [q[0], q[1], z]
            }
        }
    }

    /// Closest point of a 3D edge to `p`, with its fraction.
    pub fn edge3_closest(&self, e: u32, p: V3) -> (V3, f64) {
        let (q, u) = self.edge3_closest_internal(e, self.iv(p));
        (self.wv(q), u)
    }
    fn edge3_closest_internal(&self, e: u32, p: V3) -> (V3, f64) {
        if let EdgeGeom::Circle { c, r } = self.edges3[e as usize].geom {
            let d = sub3(p, c);
            let u = (d[1].atan2(d[0]) / TAU).rem_euclid(1.0);
            return (add3(c, [r * (u * TAU).cos(), r * (u * TAU).sin(), 0.0]), u);
        }
        if let EdgeGeom::Segment { a, b } = self.edges3[e as usize].geom {
            let d = sub3(b, a);
            let u = (dot3(sub3(p, a), d) / dot3(d, d)).clamp(0.0, 1.0);
            return (add3(a, mul3(d, u)), u);
        }
        let ov = &self.stack.as_ref().expect("a layer stack").overlay;
        match self.edges3[e as usize].geom {
            EdgeGeom::Segment { .. } | EdgeGeom::Circle { .. } => unreachable!(),
            EdgeGeom::Vertical { vertex, z0, z1 } => {
                let q = ov.vertices[vertex as usize].p;
                let z = p[2].clamp(z0, z1);
                ([q[0], q[1], z], if z1 > z0 { (z - z0) / (z1 - z0) } else { 0.0 })
            }
            EdgeGeom::Horizontal { edge, z } => {
                let ed = &ov.edges[edge as usize];
                let (t, q) = ov.edge_closest(edge, [p[0], p[1]]);
                ([q[0], q[1], z], ((t - ed.t0) / (ed.t1 - ed.t0)).clamp(0.0, 1.0))
            }
        }
    }

    /// Topology checks of a 3D complex: every edge is shared by at least two distinct faces, and the edges around every face add up to
    /// its boundary (lateral: top and bottom plus the two verticals; cap: its cell's loops).
    pub(crate) fn check_topology3(&self) -> Result<(), String> {
        let scale = self.regions.iter().map(|r| r.bbox.diagonal()).fold(0.0, f64::max).max(1e-300);
        if self.stack.is_none() {
            let mut sum = vec![0.0f64; self.faces.len()];
            for e in &self.edges3 {
                if e.faces.len() < 2 {
                    return Err(format!("edge '{}' is used by a single face", e.name));
                }
                e.faces.iter().for_each(|&f| sum[f as usize] += e.length);
            }
            for (i, f) in self.faces.iter().enumerate() {
                let FaceGeom::Planar { plane, cell } = f.geom else { continue };
                let cx = &self.planar[plane as usize].cx;
                let want: f64 = cx.regions[cell as usize]
                    .loops
                    .iter()
                    .flat_map(|&l| cx.loops[l as usize].edges.iter())
                    .map(|&(e, _)| cx.edges[e as usize].length)
                    .sum();
                if (sum[i] - want).abs() > 32.0 * self.tol + 1e-10 * scale {
                    return Err(format!(
                        "the edges around face '{}' have length {} but its boundary has {want}",
                        f.name, sum[i]
                    ));
                }
            }
            return Ok(());
        }
        let st = self.stack.as_ref().ok_or("no layer stack")?;
        let ov = &st.overlay;
        let scale = self.regions.iter().map(|r| r.bbox.diagonal()).fold(0.0, f64::max).max(1e-300);
        let mut sum = vec![0.0f64; self.faces.len()];
        for e in &self.edges3 {
            if e.faces.len() < 2 {
                return Err(format!("edge '{}' is used by a single face", e.name));
            }
            for &f in &e.faces {
                sum[f as usize] += e.length;
            }
        }
        for (i, f) in self.faces.iter().enumerate() {
            let want = match f.geom {
                FaceGeom::Lateral { edge, z0, z1 } => {
                    let e = &ov.edges[edge as usize];
                    2.0 * e.length + if e.v0 == e.v1 { 0.0 } else { 2.0 * (z1 - z0) }
                }
                FaceGeom::Quadric(b) => b.perimeter(),
                FaceGeom::Planar { .. } => 0.0,
                FaceGeom::Cap { cell, .. } => ov.regions[cell as usize]
                    .loops
                    .iter()
                    .flat_map(|&l| ov.loops[l as usize].edges.iter())
                    .map(|&(e, _)| ov.edges[e as usize].length)
                    .sum(),
            };
            if (sum[i] - want).abs() > 32.0 * self.tol + 1e-10 * scale {
                return Err(format!(
                    "the edges around face '{}' have length {} but its boundary has {want}",
                    f.name, sum[i]
                ));
            }
        }
        Ok(())
    }
}
