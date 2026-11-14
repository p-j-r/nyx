//! The triangulation core shared by the geom path (`geomdrive`): chains of boundary points become a constrained Delaunay triangulation with
//! Ruppert refinement (via the `spade` crate), cells labelled by a domain function; and the sweep of a section along y into layered
//! tetrahedra (`extrude_y`).

use crate::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use spade::{AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation};
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Layer {
    pub thickness: f64,
    pub n: usize,
    /// added to the domain numbers of this layer's cells
    pub domain_offset: i32,
    /// when set, every cell of the layer takes this domain number whatever its 2D region (a layer of cladding above and below a patterned
    /// slab: the pattern's regions are all the cladding there); the offset is then not added
    #[serde(default)]
    pub fill: Option<i32>,
}

/// One closed polyline edge chain: sample points and the boundary label of every segment.
pub(crate) struct Chain {
    pub pts: Vec<[f64; 2]>,
    pub label: i32,
    /// ties between chains on the same segment go to the larger `part`
    pub part: usize,
}

pub(crate) fn dist_point_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let (ab, ap) = ([b[0] - a[0], b[1] - a[1]], [p[0] - a[0], p[1] - a[1]]);
    let l2 = ab[0] * ab[0] + ab[1] * ab[1];
    let t = if l2 > 0.0 { ((ap[0] * ab[0] + ap[1] * ab[1]) / l2).clamp(0.0, 1.0) } else { 0.0 };
    ((a[0] + t * ab[0] - p[0]).powi(2) + (a[1] + t * ab[1] - p[1]).powi(2)).sqrt()
}

/// One lattice fill: a region that wants elements of size `h` (finer than the global size) inside the box `lo..hi`; `owner` says whether a point is
/// painted by this region.
pub(crate) struct Fill<'a> {
    pub lo: [f64; 2],
    pub hi: [f64; 2],
    pub h: f64,
    pub owner: Box<dyn Fn([f64; 2]) -> bool + 'a>,
}

pub(crate) struct MeshParams {
    /// global element size (the area limit of the refinement)
    pub h: f64,
    /// the smallest size asked for anywhere (spacing of the lattice test) and the largest (reach of the boundary test)
    pub h_fine_min: f64,
    pub h_reach: f64,
    pub min_angle_deg: f64,
    /// bounding box of the chains
    pub lo: [f64; 2],
    pub hi: [f64; 2],
}

/// The shared tail of the meshers: constrained Delaunay triangulation of the chains (conforming, labelled), lattice points, Ruppert refinement,
/// cells labelled by `domain_at` (None: no cell) and boundary facets labelled by the chain they lie on.
pub(crate) fn mesh_chains(
    chains: &[Chain],
    fills: &[Fill],
    mp: &MeshParams,
    domain_at: &dyn Fn([f64; 2]) -> Option<i32>,
) -> Result<SimplexMesh, String> {
    let (lo, hi, h_global) = (mp.lo, mp.hi, mp.h);
    let scale = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    // triangulate in coordinates centred on the shapes and scaled to about one, so rounding errors are
    // relative to the size of the model whatever its units or position (nanometre structures
    // sit at 1e-6 offsets and would otherwise get slivers)
    let c0 = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
    let nrm = |p: [f64; 2]| [(p[0] - c0[0]) / scale, (p[1] - c0[1]) / scale];
    let den = |p: [f64; 2]| [p[0] * scale + c0[0], p[1] * scale + c0[1]];
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> = ConstrainedDelaunayTriangulation::new();
    for ch in chains {
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
    let h_fine_min = mp.h_fine_min.max(1e-300);
    let reach = 0.5 * mp.h_reach;
    let gcell = reach.max(1e-9 * scale);
    let gk = |x: f64| (x / gcell).floor() as i64;
    let lcell = (0.6 * h_fine_min).max(1e-9 * scale);
    let lk = |x: f64| (x / lcell).floor() as i64;
    type Segment = ([f64; 2], [f64; 2]);
    let mut seg_grid: HashMap<(i64, i64), Vec<Segment>> = HashMap::new();
    for ch in chains {
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
    for fill in fills {
        let (h, blo, bhi) = (fill.h, fill.lo, fill.hi);
        let dy = h * 3f64.sqrt() / 2.0;
        let mut row = 0;
        let mut y = blo[1] + dy;
        while y < bhi[1] {
            let mut x = blo[0] + if row % 2 == 0 { 0.5 * h } else { h };
            while x < bhi[0] {
                if (fill.owner)([x, y]) {
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
    let params = RefinementParameters::<f64>::new()
        // no exclude_outer_faces: it decides inside/outside by constraint parity, which nested or
        // overlapping shapes break; unwanted cells are dropped by the shape test below instead
        .keep_constraint_edges()
        .with_angle_limit(AngleLimit::from_deg(mp.min_angle_deg.clamp(0.0, 33.0)))
        .with_max_allowed_area(0.5 * (h_global / scale) * (h_global / scale));
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
        if let Some(domain) = domain_at(cen) {
            cells.extend_from_slice(&[a.fix().index() as u32, b.fix().index() as u32, c.fix().index() as u32]);
            cell_labels.push(domain);
        }
    }
    if cells.is_empty() {
        return Err("the geometry encloses no area".into());
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
        for ch in chains {
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

/// Sweeps the x-z cross-section of the wafer frame (a 2D mesh with the columns `(x, z)`) along the propagation axis y through the given
/// layers: a point `(x, z)` becomes `(x, y, z)`. Each prism is split into three tetrahedra consistently with the vertex numbers, so
/// neighbours always match. The side walls keep the 2D boundary numbers, the end at `y0` is `nb + 1` and the far end `nb + 2`, with `nb` the
/// largest 2D boundary number. The domain offsets of the layers are added to the cells' domain numbers.
pub fn extrude_y(m2: &SimplexMesh, layers: &[Layer], y0: f64) -> Result<SimplexMesh, String> {
    extrude_with(m2, layers, y0, |p, t| [p[0], t, p[1]])
}

fn extrude_with(m2: &SimplexMesh, layers: &[Layer], z0: f64, place: impl Fn([f64; 2], f64) -> [f64; 3]) -> Result<SimplexMesh, String> {
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
            points.extend_from_slice(&place([p[0], p[1]], *zk));
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
                let l = &layers[layer_of_slab[k]];
                cell_labels.push(l.fill.unwrap_or(m2.cell_labels[c] + l.domain_offset));
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
