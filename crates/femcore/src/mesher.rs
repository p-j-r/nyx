//! A small mesher: 2D shapes (rectangles, circles, polygons, painted in order so later shapes
//! override earlier ones) become a constrained Delaunay triangulation with Ruppert refinement
//! (via the `spade` crate), and 2D meshes can be extruded into layered 3D tetrahedral meshes.
//!
//! Domain numbers come from the shapes; boundary numbers are those of the shape edges (a
//! rectangle numbers its edges 1 x=min, 2 y=min, 3 y=max, 4 x=max; a circle is one
//! boundary; each polygon edge is its own boundary) continued across shapes in order.

use crate::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use spade::{AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation};
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Shape {
    Rect { x0: f64, y0: f64, x1: f64, y1: f64 },
    Circle { cx: f64, cy: f64, r: f64 },
    Polygon(Vec<[f64; 2]>),
}

impl Shape {
    pub fn name(&self) -> &'static str {
        match self {
            Shape::Rect { .. } => "Rectangle",
            Shape::Circle { .. } => "Circle",
            Shape::Polygon(_) => "Polygon",
        }
    }

    pub(crate) fn contains(&self, p: [f64; 2]) -> bool {
        match self {
            Shape::Rect { x0, y0, x1, y1 } => p[0] >= x0.min(*x1) && p[0] <= x0.max(*x1) && p[1] >= y0.min(*y1) && p[1] <= y0.max(*y1),
            Shape::Circle { cx, cy, r } => (p[0] - cx).powi(2) + (p[1] - cy).powi(2) <= r * r,
            Shape::Polygon(v) => {
                let mut inside = false;
                let n = v.len();
                for i in 0..n {
                    let (a, b) = (v[i], v[(i + 1) % n]);
                    if (a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0] {
                        inside = !inside;
                    }
                }
                inside
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Part {
    pub shape: Shape,
    /// domain number given to the cells of this shape (0 or negative: a hole, no cells)
    pub domain: i32,
    /// element size inside and along this shape (default: the global size)
    pub h: Option<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Spec2d {
    pub parts: Vec<Part>,
    /// global target element size
    pub h: f64,
    pub min_angle_deg: f64,
}

impl Default for Spec2d {
    fn default() -> Self {
        Spec2d { parts: vec![], h: 0.1, min_angle_deg: 25.0 }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Layer {
    pub thickness: f64,
    pub n: usize,
    /// added to the domain numbers of this layer's cells
    pub domain_offset: i32,
}

/// One closed polyline edge chain: sample points and the boundary label of every segment.
struct Chain {
    pts: Vec<[f64; 2]>,
    label: i32,
    part: usize,
}

pub(crate) fn dist_point_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
    let l2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if l2 > 0.0 { ((ap[0] * ab[0] + ap[1] * ab[1]) / l2).clamp(0.0, 1.0) } else { 0.0 };
    ((a[0] + t * ab[0] - p[0]).powi(2) + (a[1] + t * ab[1] - p[1]).powi(2)).sqrt()
}

fn sample(a: [f64; 2], b: [f64; 2], h: f64) -> Vec<[f64; 2]> {
    let len = ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
    let n = (len / h).ceil().max(1.0) as usize;
    (0..=n).map(|i| [a[0] + (b[0] - a[0]) * i as f64 / n as f64, a[1] + (b[1] - a[1]) * i as f64 / n as f64]).collect()
}

/// Makes the polylines of all shapes mutually conforming: vertices that coincide up to rounding are
/// merged, and every vertex that lies on another segment (a T-junction, a collinear overlap) or
/// every crossing of two segments becomes a vertex of both. Without this, overlapping shapes
/// leave vertices a rounding error apart, and the triangulation contains slivers of zero area.
fn conform_chains(chains: &mut [Chain], scale: f64) {
    let tol = 1e-9 * scale;
    // canonical vertices on a hash grid
    let vcell = tol * 8.0;
    let mut verts: HashMap<(i64, i64), Vec<[f64; 2]>> = HashMap::new();
    let key = |p: [f64; 2]| ((p[0] / vcell).floor() as i64, (p[1] / vcell).floor() as i64);
    let canon = |p: [f64; 2], verts: &mut HashMap<(i64, i64), Vec<[f64; 2]>>| -> [f64; 2] {
        let (kx, ky) = key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                if let Some(list) = verts.get(&(kx + dx, ky + dy)) {
                    for q in list {
                        if (q[0] - p[0]).abs() <= tol && (q[1] - p[1]).abs() <= tol {
                            return *q;
                        }
                    }
                }
            }
        }
        verts.entry((kx, ky)).or_default().push(p);
        p
    };
    for ch in chains.iter_mut() {
        for p in ch.pts.iter_mut() {
            *p = canon(*p, &mut verts);
        }
        ch.pts.dedup();
    }
    struct Seg {
        chain: usize,
        idx: usize,
        a: [f64; 2],
        b: [f64; 2],
    }
    let mut segs: Vec<Seg> = vec![];
    for (ci, ch) in chains.iter().enumerate() {
        for (i, w) in ch.pts.windows(2).enumerate() {
            segs.push(Seg { chain: ci, idx: i, a: w[0], b: w[1] });
        }
    }
    // broad phase: segments by the grid cells of their bounding boxes
    let gs = (scale / 256.0).max(tol * 16.0);
    let gk = |x: f64| (x / gs).floor() as i64;
    let mut grid: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (i, sg) in segs.iter().enumerate() {
        let (x0, x1) = (sg.a[0].min(sg.b[0]) - tol, sg.a[0].max(sg.b[0]) + tol);
        let (y0, y1) = (sg.a[1].min(sg.b[1]) - tol, sg.a[1].max(sg.b[1]) + tol);
        for cx in gk(x0)..=gk(x1) {
            for cy in gk(y0)..=gk(y1) {
                grid.entry((cx, cy)).or_default().push(i);
            }
        }
    }
    let mut extras: HashMap<usize, Vec<[f64; 2]>> = HashMap::new();
    let mut seen: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
    let add = |si: usize, p: [f64; 2], extras: &mut HashMap<usize, Vec<[f64; 2]>>| {
        let sg = &segs[si];
        if p == sg.a || p == sg.b {
            return;
        }
        extras.entry(si).or_default().push(p);
    };
    for list in grid.values() {
        for x in 0..list.len() {
            for y in x + 1..list.len() {
                let (i, j) = (list[x].min(list[y]), list[x].max(list[y]));
                if !seen.insert((i, j)) {
                    continue;
                }
                let (si, sj) = (&segs[i], &segs[j]);
                // endpoints of one lying on the other (T-junctions, collinear overlaps)
                for (e, host) in [(sj.a, i), (sj.b, i), (si.a, j), (si.b, j)] {
                    let h = &segs[host];
                    if dist_point_segment(e, h.a, h.b) < tol {
                        let c = canon(e, &mut verts);
                        add(host, c, &mut extras);
                    }
                }
                // proper crossings
                let (r, t) = ([si.b[0] - si.a[0], si.b[1] - si.a[1]], [sj.b[0] - sj.a[0], sj.b[1] - sj.a[1]]);
                let d = r[0] * t[1] - r[1] * t[0];
                let (lr, lt) = ((r[0] * r[0] + r[1] * r[1]).sqrt(), (t[0] * t[0] + t[1] * t[1]).sqrt());
                if d.abs() > 1e-12 * lr * lt {
                    let w = [sj.a[0] - si.a[0], sj.a[1] - si.a[1]];
                    let u = (w[0] * t[1] - w[1] * t[0]) / d;
                    let v = (w[0] * r[1] - w[1] * r[0]) / d;
                    let (eu, ev) = (tol / lr, tol / lt);
                    if u > eu && u < 1.0 - eu && v > ev && v < 1.0 - ev {
                        let p = canon([si.a[0] + u * r[0], si.a[1] + u * r[1]], &mut verts);
                        add(i, p, &mut extras);
                        add(j, p, &mut extras);
                    }
                }
            }
        }
    }
    if extras.is_empty() {
        return;
    }
    // rebuild the chains with the extra vertices inserted in order along their segments
    let mut per_chain: Vec<Vec<[f64; 2]>> = vec![vec![]; chains.len()];
    let mut next_seg = 0;
    for (ci, ch) in chains.iter().enumerate() {
        let out = &mut per_chain[ci];
        for (i, w) in ch.pts.windows(2).enumerate() {
            debug_assert_eq!((segs[next_seg].chain, segs[next_seg].idx), (ci, i));
            out.push(w[0]);
            if let Some(ex) = extras.get(&next_seg) {
                let mut ex = ex.clone();
                let dir = [w[1][0] - w[0][0], w[1][1] - w[0][1]];
                ex.sort_by(|p, q| {
                    let tp = (p[0] - w[0][0]) * dir[0] + (p[1] - w[0][1]) * dir[1];
                    let tq = (q[0] - w[0][0]) * dir[0] + (q[1] - w[0][1]) * dir[1];
                    tp.partial_cmp(&tq).unwrap()
                });
                ex.dedup();
                out.extend(ex);
            }
            next_seg += 1;
        }
        if let Some(l) = ch.pts.last() {
            out.push(*l);
        }
    }
    for (ch, pts) in chains.iter_mut().zip(per_chain) {
        ch.pts = pts;
    }
}

/// Meshes the shapes of `spec` with triangles.
pub fn mesh_2d(spec: &Spec2d) -> Result<SimplexMesh, String> {
    if spec.parts.is_empty() {
        return Err("add at least one shape".into());
    }
    if spec.h <= 0.0 || !spec.h.is_finite() {
        return Err("the element size must be positive".into());
    }
    // chains with their boundary labels
    let mut chains: Vec<Chain> = vec![];
    let mut next_label = 1;
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for (pi, part) in spec.parts.iter().enumerate() {
        let h = part.h.unwrap_or(spec.h).min(spec.h.max(part.h.unwrap_or(spec.h)));
        if h <= 0.0 {
            return Err("element sizes must be positive".into());
        }
        let corners: Vec<[f64; 2]> = match &part.shape {
            Shape::Rect { x0, y0, x1, y1 } => {
                if x0 == x1 || y0 == y1 {
                    return Err("a rectangle needs a non-zero width and height".into());
                }
                let (xa, xb, ya, yb) = (x0.min(*x1), x0.max(*x1), y0.min(*y1), y0.max(*y1));
                // edges: left, bottom, top, right (left, bottom, top, right)
                for (a, b) in [([xa, ya], [xa, yb]), ([xa, ya], [xb, ya]), ([xa, yb], [xb, yb]), ([xb, ya], [xb, yb])] {
                    chains.push(Chain { pts: sample(a, b, h), label: next_label, part: pi });
                    next_label += 1;
                }
                vec![[xa, ya], [xb, yb]]
            }
            Shape::Circle { cx, cy, r } => {
                if *r <= 0.0 {
                    return Err("a circle needs a positive radius".into());
                }
                let n = ((2.0 * std::f64::consts::PI * r / h).ceil() as usize).max(12);
                let pts: Vec<[f64; 2]> = (0..=n).map(|i| {
                    let t = 2.0 * std::f64::consts::PI * (i % n) as f64 / n as f64;
                    [cx + r * t.cos(), cy + r * t.sin()]
                }).collect();
                chains.push(Chain { pts, label: next_label, part: pi });
                next_label += 1;
                vec![[cx - r, cy - r], [cx + r, cy + r]]
            }
            Shape::Polygon(v) => {
                if v.len() < 3 {
                    return Err("a polygon needs at least three vertices".into());
                }
                for i in 0..v.len() {
                    chains.push(Chain { pts: sample(v[i], v[(i + 1) % v.len()], h), label: next_label, part: pi });
                    next_label += 1;
                }
                v.iter().fold(vec![[f64::INFINITY; 2], [f64::NEG_INFINITY; 2]], |mut b, p| {
                    for k in 0..2 {
                        b[0][k] = b[0][k].min(p[k]);
                        b[1][k] = b[1][k].max(p[k]);
                    }
                    b
                })
            }
        };
        for k in 0..2 {
            lo[k] = lo[k].min(corners[0][k]);
            hi[k] = hi[k].max(corners[1][k]);
        }
    }
    let scale = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    conform_chains(&mut chains, scale);
    if std::env::var_os("NYX_MESH_DEBUG").is_some() {
        for ch in &chains {
            for w in ch.pts.windows(2) {
                if (w[0][0] - 3.5e-6).abs() < 2e-7 && (w[0][1] - 9e-6).abs() < 1e-9 {
                    eprintln!("chain part {} label {}: {:?} -> {:?}", ch.part, ch.label, w[0], w[1]);
                }
            }
        }
    }
    // triangulate in coordinates centred on the shapes and scaled to about one, so rounding errors are
    // relative to the size of the model whatever its units or position (nanometre structures
    // sit at 1e-6 offsets and would otherwise get slivers)
    let c0 = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
    let nrm = |p: [f64; 2]| [(p[0] - c0[0]) / scale, (p[1] - c0[1]) / scale];
    let den = |p: [f64; 2]| [p[0] * scale + c0[0], p[1] * scale + c0[1]];
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> = ConstrainedDelaunayTriangulation::new();
    for ch in &chains {
        let mut prev = None;
        for p in &ch.pts {
            let pn = nrm(*p);
            let v = cdt.insert(Point2::new(pn[0], pn[1])).map_err(|e| format!("meshing failed: {e:?}"))?;
            if let Some(q) = prev {
                if q != v {
                    cdt.add_constraint_and_split(q, v, |p| p);
                }
            }
            prev = Some(v);
        }
    }
    // interior lattice points where a shape asks for a finer size than the global one. A point is
    // inserted only where this shape is the one that paints the cell (later shapes override), so
    // overlapping shapes do not interleave two lattices, and never close to a boundary sample or to
    // an earlier lattice point (which would give needle triangles)
    let h_fine_min = spec.parts.iter().filter_map(|p| p.h).fold(spec.h, f64::min).max(1e-300);
    let reach = 0.5 * spec.parts.iter().filter_map(|p| p.h).fold(0.0, f64::max).max(spec.h);
    let gcell = reach.max(1e-9 * scale);
    let gk = |x: f64| (x / gcell).floor() as i64;
    let lcell = (0.6 * h_fine_min).max(1e-9 * scale);
    let lk = |x: f64| (x / lcell).floor() as i64;
    type Segment = ([f64; 2], [f64; 2]);
    let mut seg_grid: HashMap<(i64, i64), Vec<Segment>> = HashMap::new();
    for ch in &chains {
        for w in ch.pts.windows(2) {
            let (x0, x1) = (w[0][0].min(w[1][0]) - reach, w[0][0].max(w[1][0]) + reach);
            let (y0, y1) = (w[0][1].min(w[1][1]) - reach, w[0][1].max(w[1][1]) + reach);
            for cx in gk(x0)..=gk(x1) {
                for cy in gk(y0)..=gk(y1) {
                    seg_grid.entry((cx, cy)).or_default().push((w[0], w[1]));
                }
            }
        }
    }
    let mut lattice_pts: HashMap<(i64, i64), Vec<[f64; 2]>> = HashMap::new();
    for (pi, part) in spec.parts.iter().enumerate() {
        if part.domain <= 0 {
            continue;
        }
        if let Some(h) = part.h.filter(|h| *h < spec.h * 0.999) {
            let (blo, bhi) = match &part.shape {
                Shape::Rect { x0, y0, x1, y1 } => ([x0.min(*x1), y0.min(*y1)], [x0.max(*x1), y0.max(*y1)]),
                Shape::Circle { cx, cy, r } => ([cx - r, cy - r], [cx + r, cy + r]),
                Shape::Polygon(v) => (
                    [v.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min), v.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min)],
                    [v.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max), v.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max)],
                ),
            };
            let dy = h * 3f64.sqrt() / 2.0;
            let mut row = 0;
            let mut y = blo[1] + dy;
            while y < bhi[1] {
                let mut x = blo[0] + if row % 2 == 0 { 0.5 * h } else { h };
                while x < bhi[0] {
                    let owner = spec.parts.iter().rposition(|p| p.shape.contains([x, y]));
                    if owner == Some(pi) {
                        let (kx, ky) = (gk(x), gk(y));
                        let (lx, ly) = (lk(x), lk(y));
                        let span = (0.6 * h / lcell).ceil() as i64 + 1;
                        let near_boundary = seg_grid.get(&(kx, ky)).is_some_and(|l| l.iter().any(|s| dist_point_segment([x, y], s.0, s.1) < 0.4 * h));
                        let mut near_point = false;
                        'o: for dx in -span..=span {
                            for dyy in -span..=span {
                                if let Some(l) = lattice_pts.get(&(lx + dx, ly + dyy)) {
                                    if l.iter().any(|q| (q[0] - x).hypot(q[1] - y) < 0.6 * h) {
                                        near_point = true;
                                        break 'o;
                                    }
                                }
                            }
                        }
                        if !near_boundary && !near_point {
                            let pn = nrm([x, y]);
                            let _ = cdt.insert(Point2::new(pn[0], pn[1]));
                            lattice_pts.entry((lx, ly)).or_default().push([x, y]);
                        }
                    }
                    x += h;
                }
                y += dy;
                row += 1;
            }
        }
    }
    let params = RefinementParameters::<f64>::new()
        // no exclude_outer_faces: it decides inside/outside by constraint parity, which nested or
        // overlapping shapes break; unwanted cells are dropped by the shape test below instead
        .keep_constraint_edges()
        .with_angle_limit(AngleLimit::from_deg(spec.min_angle_deg.clamp(0.0, 33.0)))
        .with_max_allowed_area(0.5 * (spec.h / scale) * (spec.h / scale));
    cdt.refine(params);

    // cells: faces whose centroid is inside some shape (the last containing shape wins)
    let nvert = cdt.num_vertices();
    let mut cells: Vec<u32> = vec![];
    let mut cell_labels: Vec<i32> = vec![];
    for f in cdt.inner_faces() {
        let [a, b, c] = f.vertices();
        let pos = [a.position(), b.position(), c.position()];
        // a face whose vertices are collinear to rounding is not a cell: Ruppert refinement can split a
        // boundary edge at a point a rounding error off the line, leaving such a sliver outside the shapes
        let area2 = ((pos[1].x - pos[0].x) * (pos[2].y - pos[0].y) - (pos[2].x - pos[0].x) * (pos[1].y - pos[0].y)).abs();
        let longest2 = (0..3).map(|i| (pos[i].x - pos[(i + 1) % 3].x).powi(2) + (pos[i].y - pos[(i + 1) % 3].y).powi(2)).fold(0.0, f64::max);
        if area2 < 1e-8 * longest2 {
            continue;
        }
        let cen = den([(pos[0].x + pos[1].x + pos[2].x) / 3.0, (pos[0].y + pos[1].y + pos[2].y) / 3.0]);
        if let Some(part) = spec.parts.iter().rev().find(|p| p.shape.contains(cen)).filter(|p| p.domain > 0) {
            cells.extend_from_slice(&[a.fix().index() as u32, b.fix().index() as u32, c.fix().index() as u32]);
            cell_labels.push(part.domain);
        }
    }
    if cells.is_empty() {
        return Err("the shapes enclose no area".into());
    }
    // compact vertex numbering
    let mut used = vec![u32::MAX; nvert];
    let mut points: Vec<f64> = vec![];
    for v in cells.iter_mut() {
        if used[*v as usize] == u32::MAX {
            let p = cdt.vertex(spade::handles::FixedVertexHandle::from_index(*v as usize)).position();
            used[*v as usize] = (points.len() / 2) as u32;
            let pd = den([p.x, p.y]);
            points.extend_from_slice(&pd);
        }
        *v = used[*v as usize];
    }
    // boundary facets: edges used by exactly one kept cell
    let mut edge_count: HashMap<(u32, u32), u32> = HashMap::new();
    for t in cells.chunks(3) {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *edge_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut boundary: Vec<(u32, u32)> = edge_count.iter().filter(|(_, &n)| n == 1).map(|(e, _)| *e).collect();
    boundary.sort_unstable();
    let (mut facets, mut facet_labels) = (vec![], vec![]);
    let tol = 1e-6 * scale;
    for (a, b) in boundary {
        let pa = [points[2 * a as usize], points[2 * a as usize + 1]];
        let pb = [points[2 * b as usize], points[2 * b as usize + 1]];
        let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0];
        // the chain segment closest to the midpoint; the last shape wins ties
        let mut best: Option<(f64, usize, i32)> = None;
        for ch in &chains {
            for s in ch.pts.windows(2) {
                let d = dist_point_segment(mid, s[0], s[1]);
                if d < tol * 10.0 + 1e-12 && best.is_none_or(|(bd, bp, _)| d < bd - tol || (d <= bd + tol && ch.part >= bp)) {
                    best = Some((d, ch.part, ch.label));
                }
            }
        }
        let lab = best.map(|b| b.2).ok_or_else(|| format!("a boundary edge near ({:.4e}, {:.4e}) belongs to no shape edge", mid[0], mid[1]))?;
        facets.extend_from_slice(&[a, b]);
        facet_labels.push(lab);
    }
    let nc = cell_labels.len();
    debug_assert_eq!(nc * 3, cells.len());
    Ok(SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels))
}

/// Extrudes a 2D mesh into tetrahedra along z through the given layers (each prism is split into
/// three tetrahedra consistently with the vertex numbers, so neighbours always match). Side
/// boundaries keep the 2D boundary numbers; the bottom is `nb + 1` and the top `nb + 2`, with
/// `nb` the largest 2D boundary number.
pub fn extrude(m2: &SimplexMesh, layers: &[Layer], z0: f64) -> Result<SimplexMesh, String> {
    if m2.dim != 2 {
        return Err("only 2D meshes can be extruded".into());
    }
    if layers.is_empty() || layers.iter().any(|l| l.thickness <= 0.0 || l.n == 0) {
        return Err("every layer needs a positive thickness and at least one element layer".into());
    }
    let n2 = m2.n_points();
    let mut z = vec![z0];
    let mut layer_of_slab = vec![];
    for (li, l) in layers.iter().enumerate() {
        for _ in 0..l.n {
            z.push(z.last().unwrap() + l.thickness / l.n as f64);
            layer_of_slab.push(li);
        }
    }
    let nslab = z.len() - 1;
    let mut points = Vec::with_capacity((nslab + 1) * n2 * 3);
    for zk in &z {
        for p in m2.points.chunks(2) {
            points.extend_from_slice(&[p[0], p[1], *zk]);
        }
    }
    let id = |v: u32, k: usize| (k * n2) as u32 + v;
    let (mut cells, mut cell_labels) = (vec![], vec![]);
    for c in 0..m2.n_cells() {
        let mut v = m2.cell(c).to_vec();
        v.sort_unstable();
        for k in 0..nslab {
            let (b0, b1, b2) = (id(v[0], k), id(v[1], k), id(v[2], k));
            let (t0, t1, t2) = (id(v[0], k + 1), id(v[1], k + 1), id(v[2], k + 1));
            for tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                cells.extend_from_slice(&tet);
                cell_labels.push(m2.cell_labels[c] + layers[layer_of_slab[k]].domain_offset);
            }
        }
    }
    let nb = m2.facet_labels.iter().copied().max().unwrap_or(0);
    let (mut facets, mut facet_labels) = (vec![], vec![]);
    for c in 0..m2.n_cells() {
        let v = m2.cell(c);
        facets.extend_from_slice(&[id(v[0], 0), id(v[1], 0), id(v[2], 0)]);
        facet_labels.push(nb + 1);
        facets.extend_from_slice(&[id(v[0], nslab), id(v[1], nslab), id(v[2], nslab)]);
        facet_labels.push(nb + 2);
    }
    for f in 0..m2.n_facets() {
        let e = m2.facet(f);
        let (a, b) = (e[0].min(e[1]), e[0].max(e[1]));
        for k in 0..nslab {
            facets.extend_from_slice(&[id(a, k), id(b, k), id(b, k + 1)]);
            facets.extend_from_slice(&[id(a, k), id(b, k + 1), id(a, k + 1)]);
            facet_labels.push(m2.facet_labels[f]);
            facet_labels.push(m2.facet_labels[f]);
        }
    }
    Ok(SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels))
}
