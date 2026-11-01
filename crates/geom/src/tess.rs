//! Display tessellation (viewport and picking only, never used for meshing): every edge is discretised once by deviation, region
//! loops become polygons that are triangulated by hole bridging and ear clipping, and each triangle / segment carries the id of its
//! region / edge so a pick maps back to the model.

use crate::complex::{Complex, FaceGeom};
use crate::eval::Evaluated;
use crate::math::*;
use crate::query::DiscretizeOpts;

#[derive(Clone, Debug, Default)]
pub struct DisplayMesh {
    pub points: Vec<V3>,
    pub tris: Vec<[u32; 3]>,
    /// region of each triangle
    pub tri_face: Vec<u32>,
    pub segs: Vec<[u32; 2]>,
    /// edge of each segment
    pub seg_edge: Vec<u32>,
}

fn area2(p: &[V3], a: u32, b: u32, c: u32) -> f64 {
    let (a, b, c) = (p[a as usize], p[b as usize], p[c as usize]);
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn in_tri(p: &[V3], t: [u32; 3], q: u32) -> bool {
    let pq = p[q as usize];
    if t.iter().any(|&v| p[v as usize][0] == pq[0] && p[v as usize][1] == pq[1]) {
        return false;
    }
    let (a, b, c) = (area2(p, t[0], t[1], q), area2(p, t[1], t[2], q), area2(p, t[2], t[0], q));
    a >= 0.0 && b >= 0.0 && c >= 0.0
}

/// Triangulates a counter-clockwise outer loop with clockwise holes (indices into `p`).
pub(crate) fn triangulate(p: &[V3], outer: Vec<u32>, mut holes: Vec<Vec<u32>>) -> Vec<[u32; 3]> {
    let mut poly = outer;
    let xmax = |h: &Vec<u32>| h.iter().map(|&i| p[i as usize][0]).fold(f64::NEG_INFINITY, f64::max);
    holes.sort_by(|a, b| xmax(b).total_cmp(&xmax(a)));
    for h in holes {
        let m = (0..h.len()).max_by(|&a, &b| p[h[a] as usize][0].total_cmp(&p[h[b] as usize][0])).unwrap();
        let mp = p[h[m] as usize];
        // nearest polygon edge hit by the ray to +x
        let (mut best, mut bx) = (None, f64::INFINITY);
        for i in 0..poly.len() {
            let (a, b) = (p[poly[i] as usize], p[poly[(i + 1) % poly.len()] as usize]);
            if (a[1] <= mp[1]) != (b[1] <= mp[1]) {
                let x = a[0] + (mp[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if x >= mp[0] && x < bx {
                    bx = x;
                    best = Some(if a[0] > b[0] { i } else { (i + 1) % poly.len() });
                }
            }
        }
        let Some(mut ip) = best else { continue };
        // a reflex vertex inside the triangle (M, hit, P) would block the bridge: take the one closest to the ray
        let hit = [bx, mp[1], 0.0];
        let mut pts = p.to_vec();
        let (hi, pi) = (pts.len() as u32, poly[ip]);
        pts.push(hit);
        let tri = [h[m], hi, pi];
        let mut bestc = f64::NEG_INFINITY;
        for (j, &v) in poly.iter().enumerate() {
            if v != pi
                && in_tri(
                    &pts,
                    if area2(&pts, tri[0], tri[1], tri[2]) >= 0.0 { tri } else { [tri[0], tri[2], tri[1]] },
                    v,
                )
            {
                let d = sub3(p[v as usize], mp);
                let c = d[0] / (d[0].hypot(d[1]).max(1e-300));
                if c > bestc {
                    bestc = c;
                    ip = j;
                }
            }
        }
        let mut merged: Vec<u32> = poly[..=ip].to_vec();
        merged.extend((0..=h.len()).map(|k| h[(m + k) % h.len()]));
        merged.extend_from_slice(&poly[ip..]);
        poly = merged;
    }
    // ear clipping
    let mut out = vec![];
    let mut guard = 0;
    while poly.len() > 3 && guard < 10 * poly.len() * poly.len() + 100 {
        guard += 1;
        let n = poly.len();
        let mut clipped = false;
        for i in 0..n {
            let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
            if area2(p, a, b, c) <= 0.0 {
                continue;
            }
            if poly.iter().any(|&q| q != a && q != b && q != c && in_tri(p, [a, b, c], q)) {
                continue;
            }
            out.push([a, b, c]);
            poly.remove(i);
            clipped = true;
            break;
        }
        if !clipped {
            // degenerate (touching bridge): drop the flattest vertex
            let i = (0..n)
                .min_by(|&x, &y| {
                    area2(p, poly[(x + n - 1) % n], poly[x], poly[(x + 1) % n])
                        .abs()
                        .total_cmp(&area2(p, poly[(y + n - 1) % n], poly[y], poly[(y + 1) % n]).abs())
                })
                .unwrap();
            poly.remove(i);
        }
    }
    if poly.len() == 3 && area2(p, poly[0], poly[1], poly[2]) > 0.0 {
        out.push([poly[0], poly[1], poly[2]]);
    }
    out
}

impl Complex {
    /// 2D: triangles carrying region ids and polyline segments carrying edge ids.
    pub fn tessellate2d(&self, max_dev: f64, big: f64) -> DisplayMesh {
        let cx = self;
        let mut m = DisplayMesh::default();
        let mut vid: Vec<Option<u32>> = vec![None; cx.vertices.len()];
        let mut pl: Vec<Vec<u32>> = vec![];
        for (e, ed) in cx.edges.iter().enumerate() {
            let pts = cx.discretize_edge(e as u32, DiscretizeOpts { max_dev, max_len: big });
            let mut ids = vec![];
            let last = pts.len() - 1;
            for (k, (_, p)) in pts.iter().enumerate() {
                let v = if k == 0 {
                    Some(ed.v0)
                } else if k == last {
                    Some(ed.v1)
                } else {
                    None
                };
                let id = match v {
                    Some(v) => *vid[v as usize].get_or_insert_with(|| {
                        m.points.push(cx.iv(cx.vertices[v as usize].p));
                        (m.points.len() - 1) as u32
                    }),
                    None => {
                        m.points.push(*p);
                        (m.points.len() - 1) as u32
                    }
                };
                ids.push(id);
            }
            for w in ids.windows(2) {
                m.segs.push([w[0], w[1]]);
                m.seg_edge.push(e as u32);
            }
            pl.push(ids);
        }
        let ring = |l: &crate::complex::Loop| -> Vec<u32> {
            let mut v = vec![];
            for &(e, fwd) in &l.edges {
                let p = &pl[e as usize];
                if fwd {
                    v.extend_from_slice(&p[..p.len() - 1]);
                } else {
                    v.extend(p[1..].iter().rev());
                }
            }
            v
        };
        for (ri, r) in cx.regions.iter().enumerate() {
            for &l in &r.loops {
                let lp = &cx.loops[l as usize];
                if lp.area <= 0.0 {
                    continue;
                }
                let holes: Vec<Vec<u32>> = r
                    .loops
                    .iter()
                    .filter(|&&h| cx.loops[h as usize].outer == Some(l))
                    .map(|&h| ring(&cx.loops[h as usize]))
                    .collect();
                for t in triangulate(&m.points, ring(lp), holes) {
                    m.tris.push(t);
                    m.tri_face.push(ri as u32);
                }
            }
        }
        m
    }
}

impl Evaluated {
    /// Triangles carrying region ids (2D) or face ids (3D) and polyline segments carrying edge ids (2D) or face ids (3D); chords
    /// deviate from the exact curves by at most `max_dev`.
    pub fn tessellate(&self, max_dev: f64) -> DisplayMesh {
        let cx = self.complex();
        let big = self.bounds().diagonal().max(1.0) * 4.0;
        if cx.dim == 2 {
            let mut m = cx.tessellate2d(max_dev, big);
            m.points.iter_mut().for_each(|p| *p = cx.wv(*p));
            return m;
        }
        let empty = Complex::default();
        let ov = cx.stack.as_ref().map_or(&empty, |s| &s.overlay);
        let flat = if cx.stack.is_some() { ov.tessellate2d(max_dev, big) } else { DisplayMesh::default() };
        let mut m = DisplayMesh::default();
        for (fi, f) in cx.faces.iter().enumerate() {
            let fi = fi as u32;
            match f.geom {
                FaceGeom::Lateral { edge, z0, z1 } => {
                    let pts = ov.discretize_edge(edge, DiscretizeOpts { max_dev, max_len: big });
                    let base = m.points.len() as u32;
                    for (_, p) in &pts {
                        m.points.push(cx.wv([p[0], p[1], z0]));
                        m.points.push(cx.wv([p[0], p[1], z1]));
                    }
                    for k in 0..pts.len() as u32 - 1 {
                        let (a, b, c, d) = (base + 2 * k, base + 2 * k + 1, base + 2 * k + 2, base + 2 * k + 3);
                        m.tris.push([a, c, b]);
                        m.tris.push([b, c, d]);
                        m.tri_face.extend([fi, fi]);
                        m.segs.extend([[a, c], [b, d]]);
                        m.seg_edge.extend([fi, fi]);
                    }
                }
                FaceGeom::Planar { plane, cell } => {
                    let pd = &cx.planar[plane as usize];
                    let flat = pd.cx.tessellate2d(max_dev, big);
                    let base = m.points.len() as u32;
                    m.points.extend(flat.points.iter().map(|p| cx.plane_lift(plane, [p[0], p[1]])));
                    for (t, &c) in flat.tris.iter().zip(&flat.tri_face) {
                        if c == cell {
                            m.tris.push([base + t[0], base + t[1], base + t[2]]);
                            m.tri_face.push(fi);
                        }
                    }
                }
                FaceGeom::Quadric(b) => {
                    let n = ((std::f64::consts::TAU
                        / (2.0 * (1.0 - max_dev / b.radius_of_curvature()).clamp(-1.0, 1.0).acos().max(1e-3)))
                    .ceil() as usize)
                        .clamp(12, 256);
                    let base = m.points.len() as u32;
                    for iu in 0..n {
                        for iv in 0..=n {
                            m.points.push(cx.wv(b.point_at(iu as f64 / n as f64, iv as f64 / n as f64)));
                        }
                    }
                    let w = n as u32 + 1;
                    for iu in 0..n as u32 {
                        for iv in 0..n as u32 {
                            let (a, bb, c, d) = (
                                base + iu * w + iv,
                                base + iu * w + iv + 1,
                                base + ((iu + 1) % n as u32) * w + iv,
                                base + ((iu + 1) % n as u32) * w + iv + 1,
                            );
                            m.tris.push([a, c, bb]);
                            m.tris.push([bb, c, d]);
                            m.tri_face.extend([fi, fi]);
                        }
                    }
                }
                FaceGeom::Cap { cell, z } => {
                    let base = m.points.len() as u32;
                    m.points.extend(flat.points.iter().map(|p| cx.wv([p[0], p[1], z])));
                    for (t, &c) in flat.tris.iter().zip(&flat.tri_face) {
                        if c == cell {
                            m.tris.push([base + t[0], base + t[1], base + t[2]]);
                            m.tri_face.push(fi);
                        }
                    }
                }
            }
        }
        m
    }
}
