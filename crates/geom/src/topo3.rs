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
            if f.len() >= 2
                || (f.len() == 1
                    && cx.faces[f[0] as usize].is_sheet()
                    && matches!(cx.faces[f[0] as usize].geom, FaceGeom::Cap { .. }))
            {
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

    /// Topology checks of a 3D complex: every edge is shared by at least two distinct faces (except free rims of open sheets),
    /// the edges around every face form closed directed boundary loops matching its exact perimeter, and every non-sheet
    /// shell around an edge pairs half-edges in opposite directions (2-manifold per shell).
    pub(crate) fn check_topology3(&self) -> Result<(), String> {
        let scale = self.regions.iter().map(|r| r.bbox.diagonal()).fold(0.0, f64::max).max(1e-300);
        let tol_s = 32.0 * self.tol + 1e-10 * scale;
        // Map (face_idx) -> list of directed half-edges (edge3_idx, fwd_on_face)
        let mut face_half_edges: Vec<Vec<(u32, bool)>> = vec![vec![]; self.faces.len()];
        let mut sum = vec![0.0f64; self.faces.len()];
        // Map (edge3_idx, face_idx) -> orientation on that face
        let mut edge_face_dir: Vec<Vec<(u32, i8)>> = vec![vec![]; self.edges3.len()];

        for (ei, e) in self.edges3.iter().enumerate() {
            if e.length.is_nan() || e.length <= 0.0 {
                return Err(format!("edge '{}' has invalid length {}", e.name, e.length));
            }
            if (e.v0 as usize) >= self.vertices.len() || (e.v1 as usize) >= self.vertices.len() {
                return Err(format!("edge '{}' references out-of-range vertex ({}, {})", e.name, e.v0, e.v1));
            }
            if norm3(sub3(self.vertices[e.v0 as usize].p, self.edge3_point(ei as u32, 0.0))) > tol_s
                || norm3(sub3(self.vertices[e.v1 as usize].p, self.edge3_point(ei as u32, 1.0))) > tol_s
            {
                return Err(format!("edge '{}' endpoints disagree with its vertex coordinates", e.name));
            }
            if matches!(e.geom, EdgeGeom::Circle { .. }) && e.v0 != e.v1 {
                return Err(format!("circle edge '{}' must start and end at its seam vertex", e.name));
            }
            if e.faces.windows(2).any(|w| w[0] == w[1]) || {
                let mut sorted = e.faces.clone();
                sorted.sort_unstable();
                sorted.windows(2).any(|w| w[0] == w[1])
            } {
                return Err(format!("edge '{}' lists a face more than once", e.name));
            }
            if self.stack.is_none() {
                if e.faces.len() < 2 {
                    return Err(format!("edge '{}' is used by a single face", e.name));
                }
            } else if e.faces.len() < 2 && !(e.faces.len() == 1 && self.faces[e.faces[0] as usize].is_sheet()) {
                return Err(format!("edge '{}' is used by a single face", e.name));
            }
            for &f in &e.faces {
                sum[f as usize] += e.length;
            }
        }

        if self.stack.is_none() {
            let mut planar_segs: Vec<Vec<(V3, V3, V3, f64)>> = vec![vec![]; self.faces.len()];
            for (fi, f) in self.faces.iter().enumerate() {
                let FaceGeom::Planar { plane, cell } = f.geom else {
                    return Err(format!("face '{}' is non-planar without a layer stack", f.name));
                };
                let pd = &self.planar[plane as usize];
                let lift = |q: V2| add3(pd.o, add3(mul3(pd.u, q[0]), mul3(pd.v, q[1])));
                for &li in &pd.cx.regions[cell as usize].loops {
                    for &(e2, fwd2) in &pd.cx.loops[li as usize].edges {
                        let ed = &pd.cx.edges[e2 as usize];
                        let c = &pd.cx.curves[ed.curve as usize];
                        let (q0, q1) =
                            if fwd2 { (c.eval(ed.t0), c.eval(ed.t1)) } else { (c.eval(ed.t1), c.eval(ed.t0)) };
                        let (p0, p1) = (lift(q0), lift(q1));
                        let d = sub3(p1, p0);
                        let len = norm3(d);
                        if len > 0.0 {
                            planar_segs[fi].push((p0, mul3(d, 1.0 / len), p1, len));
                        }
                    }
                }
            }
            for (ei, e) in self.edges3.iter().enumerate() {
                let EdgeGeom::Segment { a, b } = e.geom else {
                    return Err(format!("edge '{}' is not a segment in a polyhedral complex", e.name));
                };
                let de = unit3(sub3(b, a));
                let mid = mul3(add3(a, b), 0.5);
                for &fi in &e.faces {
                    let mut matched: Option<bool> = None;
                    for &(p0, dir2, _, len2) in &planar_segs[fi as usize] {
                        let rel = sub3(mid, p0);
                        let t = dot3(rel, dir2);
                        if t < -tol_s || t > len2 + tol_s {
                            continue;
                        }
                        let perp = norm3(sub3(rel, mul3(dir2, t)));
                        let par = dot3(de, dir2);
                        if perp <= tol_s && par.abs() > 1.0 - 1e-5 {
                            if matched.is_some() {
                                return Err(format!(
                                    "edge '{}' matches multiple boundary edges on face '{}'",
                                    e.name, self.faces[fi as usize].name
                                ));
                            }
                            matched = Some(par > 0.0);
                        }
                    }
                    let Some(fwd) = matched else {
                        return Err(format!(
                            "edge '{}' does not lie on a boundary edge of face '{}'",
                            e.name, self.faces[fi as usize].name
                        ));
                    };
                    face_half_edges[fi as usize].push((ei as u32, fwd));
                    edge_face_dir[ei].push((fi, if fwd { 1 } else { -1 }));
                }
            }
        } else {
            let st = self.stack.as_ref().ok_or("no layer stack")?;
            let ov = &st.overlay;
            let mut cap_dir: Vec<BTreeMap<u32, bool>> = vec![BTreeMap::new(); self.faces.len()];
            for (fi, f) in self.faces.iter().enumerate() {
                if let FaceGeom::Cap { cell, .. } = f.geom {
                    for &li in &ov.regions[cell as usize].loops {
                        for &(e2, fwd2) in &ov.loops[li as usize].edges {
                            cap_dir[fi].insert(e2, fwd2);
                        }
                    }
                }
            }
            for (ei, e) in self.edges3.iter().enumerate() {
                for &fi in &e.faces {
                    let f = &self.faces[fi as usize];
                    let fwd = match (&f.geom, &e.geom) {
                        (FaceGeom::Lateral { edge: fe, z0, z1 }, EdgeGeom::Horizontal { edge: ee, z })
                            if *ee == *fe =>
                        {
                            if (*z - *z0).abs() <= tol_s {
                                true
                            } else if (*z - *z1).abs() <= tol_s {
                                false
                            } else {
                                return Err(format!(
                                    "horizontal edge '{}' does not lie on z0/z1 of lateral face '{}'",
                                    e.name, f.name
                                ));
                            }
                        }
                        (FaceGeom::Lateral { edge: fe, z0, z1 }, EdgeGeom::Vertical { vertex, z0: vz0, z1: vz1 }) => {
                            let ov_e = &ov.edges[*fe as usize];
                            if ov_e.v0 == ov_e.v1 || *vz0 < *z0 - tol_s || *vz1 > *z1 + tol_s {
                                return Err(format!(
                                    "vertical edge '{}' is out of bounds for lateral face '{}'",
                                    e.name, f.name
                                ));
                            }
                            if *vertex == ov_e.v1 {
                                true
                            } else if *vertex == ov_e.v0 {
                                false
                            } else {
                                return Err(format!(
                                    "vertical edge '{}' vertex does not bound lateral face '{}'",
                                    e.name, f.name
                                ));
                            }
                        }
                        (FaceGeom::Cap { z, .. }, EdgeGeom::Horizontal { edge: ee, z: ez })
                            if (*ez - *z).abs() <= tol_s =>
                        {
                            let Some(&fwd2) = cap_dir[fi as usize].get(ee) else {
                                return Err(format!(
                                    "horizontal edge '{}' is not part of cap face '{}'",
                                    e.name, f.name
                                ));
                            };
                            fwd2
                        }
                        (FaceGeom::Quadric(b), EdgeGeom::Circle { c: ec, r: er }) => match *b {
                            crate::prim::Body::Disk { c, r, up } => {
                                if norm3(sub3(*ec, c)) > tol_s || (*er - r).abs() > tol_s {
                                    return Err(format!(
                                        "circle edge '{}' does not match disk face '{}'",
                                        e.name, f.name
                                    ));
                                }
                                up
                            }
                            crate::prim::Body::Annulus { c, r_in, r_out, up } => {
                                if norm3(sub3(*ec, c)) > tol_s {
                                    return Err(format!(
                                        "circle edge '{}' center does not match annulus face '{}'",
                                        e.name, f.name
                                    ));
                                }
                                if (*er - r_out).abs() <= tol_s {
                                    up
                                } else if (*er - r_in).abs() <= tol_s {
                                    !up
                                } else {
                                    return Err(format!(
                                        "circle edge '{}' radius {er} does not match annulus face '{}'",
                                        e.name, f.name
                                    ));
                                }
                            }
                            crate::prim::Body::CylSide { c, r, h, inward } => {
                                if (ec[0] - c[0]).abs() > tol_s
                                    || (ec[1] - c[1]).abs() > tol_s
                                    || (*er - r).abs() > tol_s
                                {
                                    return Err(format!(
                                        "circle edge '{}' does not match cylinder side '{}'",
                                        e.name, f.name
                                    ));
                                }
                                if (ec[2] - c[2]).abs() <= tol_s {
                                    !inward
                                } else if (ec[2] - (c[2] + h)).abs() <= tol_s {
                                    inward
                                } else {
                                    return Err(format!(
                                        "circle edge '{}' z = {} is not on cylinder side '{}'",
                                        e.name, ec[2], f.name
                                    ));
                                }
                            }
                            crate::prim::Body::ConeSide { c, r0, r1, h, flip } => {
                                if (ec[0] - c[0]).abs() > tol_s || (ec[1] - c[1]).abs() > tol_s {
                                    return Err(format!(
                                        "circle edge '{}' axis does not match cone side '{}'",
                                        e.name, f.name
                                    ));
                                }
                                if (ec[2] - c[2]).abs() <= tol_s && (*er - r0).abs() <= tol_s {
                                    !flip
                                } else if (ec[2] - (c[2] + h)).abs() <= tol_s && (*er - r1).abs() <= tol_s {
                                    flip
                                } else {
                                    return Err(format!(
                                        "circle edge '{}' does not match cone side '{}'",
                                        e.name, f.name
                                    ));
                                }
                            }
                            _ => {
                                return Err(format!(
                                    "quadric face '{}' cannot be bounded by edge '{}'",
                                    f.name, e.name
                                ));
                            }
                        },
                        _ => {
                            return Err(format!("edge '{}' has incompatible geometry for face '{}'", e.name, f.name));
                        }
                    };
                    face_half_edges[fi as usize].push((ei as u32, fwd));
                    edge_face_dir[ei].push((fi, if fwd { 1 } else { -1 }));
                }
            }
        }

        // Verify exact face perimeter length and closed directed boundary loops per face.
        for (i, f) in self.faces.iter().enumerate() {
            let want = match f.geom {
                FaceGeom::Planar { plane, cell } => {
                    let cx = &self.planar[plane as usize].cx;
                    cx.regions[cell as usize]
                        .loops
                        .iter()
                        .flat_map(|&l| cx.loops[l as usize].edges.iter())
                        .map(|&(e, _)| cx.edges[e as usize].length)
                        .sum()
                }
                FaceGeom::Lateral { edge, z0, z1 } => {
                    let ov = &self.stack.as_ref().unwrap().overlay;
                    let e = &ov.edges[edge as usize];
                    2.0 * e.length + if e.v0 == e.v1 { 0.0 } else { 2.0 * (z1 - z0) }
                }
                FaceGeom::Quadric(b) => b.perimeter(),
                FaceGeom::Cap { cell, .. } => {
                    let ov = &self.stack.as_ref().unwrap().overlay;
                    ov.regions[cell as usize]
                        .loops
                        .iter()
                        .flat_map(|&l| ov.loops[l as usize].edges.iter())
                        .map(|&(e, _)| ov.edges[e as usize].length)
                        .sum()
                }
            };
            if (sum[i] - want).abs() > tol_s {
                return Err(format!(
                    "the edges around face '{}' have length {} but its boundary has {want}",
                    f.name, sum[i]
                ));
            }
            let mut deg: BTreeMap<u32, i32> = BTreeMap::new();
            for &(ei, fwd) in &face_half_edges[i] {
                let e = &self.edges3[ei as usize];
                let (vs, ve) = if fwd { (e.v0, e.v1) } else { (e.v1, e.v0) };
                *deg.entry(vs).or_insert(0) += 1;
                *deg.entry(ve).or_insert(0) -= 1;
            }
            if let Some((&v, &bal)) = deg.iter().find(|&(_, &b)| b != 0) {
                return Err(format!(
                    "face '{}' boundary edge loop is not closed at vertex '{}' (net degree {bal})",
                    f.name, self.vertices[v as usize].name
                ));
            }
        }

        // Verify 2-manifold shell edge pairing: per shell (Option<u32>), every undirected edge is used by
        // balanced opposite-orientation half-edge pairs (`2k` half-edges summing to `0`; `k = 1` normally, `k > 1` at
        // pinched CSG edges where two closed lobes meet along a shared 1D edge). Zero-thickness interior sheets with
        // `first == second` are two-sided internal baffles embedded inside one region and do not bound a 3D shell.
        for (ei, dirs) in edge_face_dir.iter().enumerate() {
            let mut per_shell: BTreeMap<Option<u32>, Vec<i8>> = BTreeMap::new();
            for &(fi, dir_on_face) in dirs {
                let f = &self.faces[fi as usize];
                if f.first == f.second {
                    continue;
                }
                per_shell.entry(f.first).or_default().push(dir_on_face);
                per_shell.entry(f.second).or_default().push(-dir_on_face);
            }
            for (shell, hs) in per_shell {
                let net: i32 = hs.iter().map(|&d| i32::from(d)).sum();
                if hs.len() < 2 || !hs.len().is_multiple_of(2) || net != 0 {
                    let sname = shell.map_or_else(|| "void".to_string(), |r| self.regions[r as usize].name.clone());
                    return Err(format!(
                        "edge '{}' is non-manifold or misoriented on shell '{sname}': half-edge orientations {hs:?}",
                        self.edges3[ei].name
                    ));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::eval::Tier;
    use crate::model::{Axis, Model};
    use crate::prim::Body;

    #[test]
    fn analytic_quadric_area_vector_and_topology3_accept_valid_and_reject_corrupt_shells() {
        let mut m = Model::new_3d();
        let air = m.box_("air", (-5.0, -5.0, -2.0), (5.0, 5.0, 6.0));
        let frustum = m.cone("frustum", (0.0, 0.0, 0.0), 1.5, 0.5, 2.0);
        let prof = m.rect("tube_prof", (2.0, 3.0), (3.0, 5.0));
        let tube = m.revolve("tube", prof, (0.0, 0.0, 0.0));
        let cone_z = m.cone("horn", (-3.0, 3.0, 0.0), 0.6, 0.2, 1.5);
        let ball = m.sphere("ball", (-3.0, -3.0, 2.0), 0.7);
        m.overlay(&[("air", air), ("frustum", frustum), ("tube", tube), ("horn", cone_z), ("ball", ball)]);
        let ev = m.eval(Tier::Complex).unwrap();
        let cx = ev.complex().clone();
        cx.check().unwrap();

        let mut mx = Model::new_3d();
        let host_x = mx.box_("host_x", (-4.0, -4.0, -4.0), (4.0, 4.0, 4.0));
        let horn_x = mx.cone_along("horn_x", Axis::X, (-2.0, 0.0, 0.0), 0.8, 0.2, 3.0);
        mx.overlay(&[("host_x", host_x), ("horn_x", horn_x)]);
        mx.eval(Tier::Complex).unwrap().complex().check().unwrap();

        // 1. Corrupting Disk `up` flag breaks analytic quadric flux closure.
        let mut flipped_disk = cx.clone();
        for f in &mut flipped_disk.faces {
            if let crate::complex::FaceGeom::Quadric(Body::Disk { up, .. }) = &mut f.geom {
                *up = !*up;
                break;
            }
        }
        assert!(flipped_disk.check().is_err(), "flipped quadric disk normal must fail check3");

        // 2. Corrupting ConeSide `flip` flag breaks analytic quadric flux closure AND shell edge orientation.
        let mut flipped_cone = cx.clone();
        for f in &mut flipped_cone.faces {
            if let crate::complex::FaceGeom::Quadric(Body::ConeSide { flip, .. }) = &mut f.geom {
                *flip = !*flip;
                break;
            }
        }
        assert!(flipped_cone.check().is_err(), "flipped cone side must fail check3");

        // 3. Swapping `first` / `second` on Sphere keeps area_vector == 0 but flips divergence volume flux.
        let mut inside_out_sphere = cx.clone();
        for f in &mut inside_out_sphere.faces {
            if matches!(f.geom, crate::complex::FaceGeom::Quadric(Body::Sphere { .. })) {
                std::mem::swap(&mut f.first, &mut f.second);
                break;
            }
        }
        assert!(inside_out_sphere.check().is_err(), "inside-out sphere must fail volume flux check3");

        // 4. Corrupting Edge3 endpoint (`v0`/`v1` mismatch) breaks closed directed face edge loop check.
        let mut broken_loop = cx.clone();
        let ei = broken_loop
            .edges3
            .iter()
            .position(|e| matches!(e.geom, crate::complex::EdgeGeom::Vertical { .. }))
            .unwrap();
        let other_v = broken_loop
            .edges3
            .iter()
            .find_map(|e| (e.v0 != broken_loop.edges3[ei].v0 && e.v0 != broken_loop.edges3[ei].v1).then_some(e.v0))
            .unwrap();
        broken_loop.edges3[ei].v1 = other_v;
        assert!(broken_loop.check().is_err(), "open/broken edge endpoint must fail check_topology3");

        // 5. Dropping a face reference from an Edge3 breaks 2-manifold shell pairing.
        let mut non_manifold = cx;
        let ei_h = non_manifold
            .edges3
            .iter()
            .position(|e| matches!(e.geom, crate::complex::EdgeGeom::Horizontal { .. }) && e.faces.len() >= 2)
            .unwrap();
        non_manifold.edges3[ei_h].faces.pop();
        assert!(non_manifold.check().is_err(), "edge missing half-edge partner must fail check_topology3");
    }
}
