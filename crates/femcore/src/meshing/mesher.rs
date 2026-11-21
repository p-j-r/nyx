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

/// A sweep [`Layer`] paired with an optional domain-label remap `(source_domain -> target_domain)` for that layer.
pub type RemappedLayer = (Layer, std::collections::BTreeMap<i32, i32>);

/// Sweeps the x-z cross-section of the wafer frame (a 2D mesh with the columns `(x, z)`) along the propagation axis y through the given
/// layers: a point `(x, z)` becomes `(x, y, z)`. Each prism is split into three tetrahedra consistently with the vertex numbers, so
/// neighbours always match. The side walls keep the 2D boundary numbers, the end at `y0` is `nb + 1` and the far end `nb + 2`, with `nb` the
/// largest 2D boundary number. The domain offsets of the layers are added to the cells' domain numbers.
pub fn extrude_y(m2: &SimplexMesh, layers: &[Layer], y0: f64) -> Result<SimplexMesh, String> {
    let mapped: Vec<RemappedLayer> =
        layers.iter().cloned().map(|l| (l, std::collections::BTreeMap::new())).collect();
    extrude_y_with_remap(m2, &mapped, y0)
}

/// Like [`extrude_y`], with an optional per-layer domain map `(source_domain -> target_domain)` applied when `l.fill` is `None`.
pub fn extrude_y_with_remap(
    m2: &SimplexMesh,
    layers: &[RemappedLayer],
    y0: f64,
) -> Result<SimplexMesh, String> {
    extrude_with(m2, layers, y0, |p, t| [p[0], t, p[1]])
}

fn extrude_with(
    m2: &SimplexMesh,
    layers: &[RemappedLayer],
    z0: f64,
    place: impl Fn([f64; 2], f64) -> [f64; 3],
) -> Result<SimplexMesh, String> {
    if m2.dim != 2 {
        return Err("only 2D meshes can be extruded".into());
    }
    if layers.is_empty() || layers.iter().any(|(l, _)| l.thickness <= 0.0 || l.n == 0) {
        return Err("every layer needs a positive thickness and at least one element layer".into());
    }
    let n2 = m2.n_points();
    let mut z = vec![z0];
    let mut layer_of_slab = vec![];
    for (li, (l, _)) in layers.iter().enumerate() {
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
            let (l, remap) = &layers[layer_of_slab[k]];
            let base = m2.cell_labels[c];
            let dom = l.fill.or_else(|| remap.get(&base).copied()).unwrap_or(base + l.domain_offset);
            for tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                cells.extend_from_slice(&tet);
                cell_labels.push(dom);
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

/// Sweeps a 3D `geom` model along the propagation axis `y` through `layers`: extracts the combined 2D `(x, z)` cross-section
/// arrangement across all layer midpoints, triangulates it with local `Size` rules and `sort_nodes`, extrudes conformal Dompierre
/// prisms from the model's `ymin`, assigns cell domain labels from `l.fill` / `l.regions` / `Evaluated::classify`, and tags outer,
/// cap, and material-interface facets via `GeomDomain3::boundary_label` so all 3D selectors and periodic pairs work directly.
pub fn sweep_model_3d(
    ev: &geom::eval::Evaluated,
    labels: &geom::labels::LabelMap,
    h_global: f64,
    max_dev: Option<f64>,
    min_angle_deg: f64,
    sizes: &crate::meshing::geomdrive3::Sizes3,
    layers: &[RemappedLayer],
) -> Result<SimplexMesh, String> {
    use crate::meshing::geomdrive::sort_nodes;
    use crate::meshing::geomdrive3::GeomDomain3;
    use crate::meshing::mesher3d::Domain3;
    use geom::query::DiscretizeOpts;

    if ev.list.dim != 3 {
        return Err("sweep_model_3d needs a 3D model".into());
    }
    if !(h_global > 0.0 && h_global.is_finite()) {
        return Err("the element size must be positive".into());
    }
    if layers.is_empty() || layers.iter().any(|(l, _)| l.thickness <= 0.0 || l.n == 0) {
        return Err("every layer needs a positive thickness and at least one element layer".into());
    }
    let dom = GeomDomain3::new_with(ev, labels, h_global, sizes)?.with_max_dev(max_dev);
    let nr = ev.regions().len();
    let y0 = (0..nr)
        .filter(|&r| dom.is_meshed(r as u32))
        .map(|r| ev.regions()[r].bbox.lo[1])
        .fold(f64::INFINITY, f64::min);
    let y0 = if y0.is_finite() { y0 } else { ev.bounds().lo[1] };

    let mut y_planes = vec![y0];
    let mut layer_of_slab = Vec::new();
    let mut sample_ys = Vec::new();
    for (li, (l, _)) in layers.iter().enumerate() {
        let ys = *y_planes.last().unwrap();
        sample_ys.push(ys + 0.5 * l.thickness);
        for _ in 0..l.n {
            let yn = *y_planes.last().unwrap() + l.thickness / l.n as f64;
            sample_ys.push(0.5 * (*y_planes.last().unwrap() + yn));
            y_planes.push(yn);
            layer_of_slab.push(li);
        }
    }
    for r in 0..nr {
        if dom.is_meshed(r as u32) {
            let b = ev.regions()[r].bbox;
            if b.hi[1] > b.lo[1] {
                sample_ys.push(0.5 * (b.lo[1] + b.hi[1]));
            }
        }
    }
    let tol = ev.tol().max(1e-12);
    sample_ys.sort_by(|a, b| a.total_cmp(b));
    sample_ys.dedup_by(|a, b| (*a - *b).abs() <= tol);

    let md = max_dev.unwrap_or(h_global / 50.0);
    let has_any_fill = layers.iter().any(|(l, _)| l.fill.is_some());

    // Build combined 2D (x, z) cross-section mesh
    let m2_raw = if ev.complex().axis == 1 && ev.complex().stack.is_some() {
        let st = ev.complex().stack.as_ref().unwrap();
        let ov = &st.overlay;
        let ncell = ov.regions.len();
        let nedge = ov.edges.len();
        let mut cell_h = vec![h_global; ncell];
        let mut cell_active = vec![false; ncell];
        for c in 0..ncell {
            for l in &st.layers {
                if let Some(r) = l.region_of[c] {
                    if dom.is_meshed(r) {
                        cell_active[c] = true;
                        cell_h[c] = cell_h[c].min(dom.region_sizes()[r as usize]);
                    } else if has_any_fill {
                        cell_active[c] = true;
                    }
                }
            }
        }
        let mut chains = Vec::with_capacity(nedge);
        let mut edge_h = vec![h_global; nedge];
        for (ei, ed) in ov.edges.iter().enumerate() {
            let eh = [ed.left, ed.right]
                .into_iter()
                .flatten()
                .map(|c| cell_h[c as usize])
                .fold(h_global, f64::min);
            edge_h[ei] = eh;
            let mut pts: Vec<[f64; 2]> = ov
                .discretize_edge(ei as u32, DiscretizeOpts { max_dev: md, max_len: eh })
                .into_iter()
                .map(|(_, p)| [p[0], p[1]])
                .collect();
            let last = pts.len() - 1;
            pts[0] = [ov.vertices[ed.v0 as usize].p[0], ov.vertices[ed.v0 as usize].p[1]];
            pts[last] = [ov.vertices[ed.v1 as usize].p[0], ov.vertices[ed.v1 as usize].p[1]];
            chains.push(Chain { pts, label: (ei + 1) as i32, part: ei });
        }
        let mut ov_bb = geom::math::Aabb::EMPTY;
        for r in &ov.regions {
            ov_bb = ov_bb.union(&r.bbox);
        }
        let fills: Vec<Fill> = (0..ncell)
            .filter(|&c| cell_active[c] && cell_h[c] < h_global * 0.999)
            .map(|c| {
                let bb = ov.regions[c].bbox;
                Fill {
                    lo: [bb.lo[0], bb.lo[1]],
                    hi: [bb.hi[0], bb.hi[1]],
                    h: cell_h[c],
                    owner: Box::new(move |p| ov.locate(p) == Some(c as u32)),
                }
            })
            .collect();
        let h_fine_min = cell_h.iter().copied().chain(edge_h).fold(h_global, f64::min).max(1e-12 * h_global);
        mesh_chains(
            &chains,
            &fills,
            &MeshParams {
                h: h_global,
                h_fine_min,
                h_reach: h_global,
                min_angle_deg,
                lo: [ov_bb.lo[0], ov_bb.lo[1]],
                hi: [ov_bb.hi[0], ov_bb.hi[1]],
            },
            &|cen| {
                let c = ov.locate(cen)?;
                cell_active[c as usize].then_some(c as i32 + 1)
            },
        )?
    } else {
        let cx = ev.complex();
        let mut xs = Vec::new();
        let mut zs = Vec::new();
        for r in 0..nr {
            if dom.is_meshed(r as u32) || has_any_fill {
                let b = ev.regions()[r].bbox;
                if !b.is_empty() {
                    xs.extend_from_slice(&[b.lo[0], b.hi[0]]);
                    zs.extend_from_slice(&[b.lo[2], b.hi[2]]);
                }
            }
        }
        for fi in 0..cx.facet_count() as u32 {
            let b = cx.facet_bbox(fi);
            if !b.is_empty() {
                xs.extend_from_slice(&[b.lo[0], b.hi[0]]);
                zs.extend_from_slice(&[b.lo[2], b.hi[2]]);
            }
        }
        for v in ev.vertices() {
            xs.push(v.p[0]);
            zs.push(v.p[2]);
        }
        xs.retain(|x| x.is_finite());
        zs.retain(|z| z.is_finite());
        xs.sort_by(|a, b| a.total_cmp(b));
        xs.dedup_by(|a, b| (*a - *b).abs() <= tol);
        zs.sort_by(|a, b| a.total_cmp(b));
        zs.dedup_by(|a, b| (*a - *b).abs() <= tol);
        if xs.len() < 2 || zs.len() < 2 {
            return Err("the 3D model has no (x, z) cross-section extent".into());
        }
        let (nx, nz) = (xs.len() - 1, zs.len() - 1);
        let mut sig: Vec<Vec<Vec<Option<u32>>>> = vec![vec![Vec::new(); nz]; nx];
        let mut active = vec![vec![false; nz]; nx];
        let mut cell_h = vec![vec![h_global; nz]; nx];
        for ix in 0..nx {
            let cx0 = 0.5 * (xs[ix] + xs[ix + 1]);
            for iz in 0..nz {
                let cz0 = 0.5 * (zs[iz] + zs[iz + 1]);
                let s: Vec<Option<u32>> = sample_ys.iter().map(|&ym| dom.region_at([cx0, ym, cz0])).collect();
                let mut act = s.iter().any(|r| r.is_some());
                if !act && has_any_fill && sample_ys.iter().any(|&ym| ev.classify([cx0, ym, cz0]).is_some()) {
                    act = true;
                }
                let mut h = h_global;
                for &r in s.iter().flatten() {
                    h = h.min(dom.region_sizes()[r as usize]);
                }
                if dom.asks_size_at() {
                    for &ym in &sample_ys {
                        if let Some(sh) = dom.size_at([cx0, ym, cz0]) {
                            h = h.min(sh);
                        }
                    }
                }
                sig[ix][iz] = s;
                active[ix][iz] = act;
                cell_h[ix][iz] = h;
            }
        }
        let subdiv = |p0: [f64; 2], p1: [f64; 2], h_seg: f64| -> Vec<[f64; 2]> {
            let len = (p1[0] - p0[0]).hypot(p1[1] - p0[1]);
            let n = (len / h_seg.max(1e-15)).ceil().max(1.0) as usize;
            (0..=n)
                .map(|i| {
                    if i == 0 {
                        p0
                    } else if i == n {
                        p1
                    } else {
                        let t = i as f64 / n as f64;
                        [p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])]
                    }
                })
                .collect()
        };
        let mut chains = Vec::new();
        for iz in 0..=nz {
            for ix in 0..nx {
                let below = (iz > 0 && active[ix][iz - 1]).then(|| &sig[ix][iz - 1]);
                let above = (iz < nz && active[ix][iz]).then(|| &sig[ix][iz]);
                if below != above {
                    let mut h_seg = h_global;
                    if iz > 0 && active[ix][iz - 1] {
                        h_seg = h_seg.min(cell_h[ix][iz - 1]);
                    }
                    if iz < nz && active[ix][iz] {
                        h_seg = h_seg.min(cell_h[ix][iz]);
                    }
                    let id = chains.len();
                    chains.push(Chain {
                        pts: subdiv([xs[ix], zs[iz]], [xs[ix + 1], zs[iz]], h_seg),
                        label: (id + 1) as i32,
                        part: id,
                    });
                }
            }
        }
        for ix in 0..=nx {
            for iz in 0..nz {
                let left = (ix > 0 && active[ix - 1][iz]).then(|| &sig[ix - 1][iz]);
                let right = (ix < nx && active[ix][iz]).then(|| &sig[ix][iz]);
                if left != right {
                    let mut h_seg = h_global;
                    if ix > 0 && active[ix - 1][iz] {
                        h_seg = h_seg.min(cell_h[ix - 1][iz]);
                    }
                    if ix < nx && active[ix][iz] {
                        h_seg = h_seg.min(cell_h[ix][iz]);
                    }
                    let id = chains.len();
                    chains.push(Chain {
                        pts: subdiv([xs[ix], zs[iz]], [xs[ix], zs[iz + 1]], h_seg),
                        label: (id + 1) as i32,
                        part: id,
                    });
                }
            }
        }
        let mut fills = Vec::new();
        let mut h_fine_min = h_global;
        for ix in 0..nx {
            for iz in 0..nz {
                if active[ix][iz] && cell_h[ix][iz] < h_global * 0.999 {
                    let (x0, x1, z0, z1) = (xs[ix], xs[ix + 1], zs[iz], zs[iz + 1]);
                    let ch = cell_h[ix][iz];
                    h_fine_min = h_fine_min.min(ch);
                    fills.push(Fill {
                        lo: [x0, z0],
                        hi: [x1, z1],
                        h: ch,
                        owner: Box::new(move |p| p[0] > x0 && p[0] < x1 && p[1] > z0 && p[1] < z1),
                    });
                }
            }
        }
        let locate_grid = |p: [f64; 2]| -> Option<i32> {
            let ix = xs.partition_point(|&x| x <= p[0]).checked_sub(1)?;
            let iz = zs.partition_point(|&z| z <= p[1]).checked_sub(1)?;
            if ix < nx && iz < nz && active[ix][iz] {
                Some((iz * nx + ix + 1) as i32)
            } else {
                None
            }
        };
        mesh_chains(
            &chains,
            &fills,
            &MeshParams {
                h: h_global,
                h_fine_min: h_fine_min.max(1e-12 * h_global),
                h_reach: h_global,
                min_angle_deg,
                lo: [xs[0], zs[0]],
                hi: [*xs.last().unwrap(), *zs.last().unwrap()],
            },
            &locate_grid,
        )?
    };

    let m2 = sort_nodes(&m2_raw);
    let n2 = m2.n_points();
    let ntri2 = m2.n_cells();
    let nslab = y_planes.len() - 1;
    let region_of_dom: HashMap<i32, u32> = labels.domains.iter().map(|d| (d.id, d.region)).collect();

    // Determine active prisms (ti, k) -> (dom_id, eff_region_id)
    let mut prism_info: Vec<Option<(i32, u32)>> = vec![None; nslab * ntri2];
    let mut used_v = vec![u32::MAX; (nslab + 1) * n2];
    for k in 0..nslab {
        let ym = 0.5 * (y_planes[k] + y_planes[k + 1]);
        let (l, remap) = &layers[layer_of_slab[k]];
        for ti in 0..ntri2 {
            let tv = m2.cell(ti);
            let cx0 = (m2.points[2 * tv[0] as usize] + m2.points[2 * tv[1] as usize] + m2.points[2 * tv[2] as usize]) / 3.0;
            let cz0 = (m2.points[2 * tv[0] as usize + 1] + m2.points[2 * tv[1] as usize + 1] + m2.points[2 * tv[2] as usize + 1]) / 3.0;
            let r_raw = dom.region_at([cx0, ym, cz0]);
            let info = if let Some(fdom) = l.fill {
                if r_raw.is_some() || ev.classify([cx0, ym, cz0]).is_some() {
                    let reff = region_of_dom.get(&fdom).copied().or(r_raw).unwrap_or(0);
                    Some((fdom, reff))
                } else {
                    None
                }
            } else if let Some(r) = r_raw {
                let base = dom.region_label(r);
                let did = remap.get(&base).copied().unwrap_or(base + l.domain_offset);
                let reff = region_of_dom.get(&did).copied().unwrap_or(r);
                Some((did, reff))
            } else {
                None
            };
            if info.is_some() {
                prism_info[k * ntri2 + ti] = info;
                for &v in tv {
                    used_v[k * n2 + v as usize] = 0;
                    used_v[(k + 1) * n2 + v as usize] = 0;
                }
            }
        }
    }

    let mut points = Vec::new();
    let mut n_pts3 = 0u32;
    for (k, &yk) in y_planes.iter().enumerate() {
        for v2 in 0..n2 {
            let slot = &mut used_v[k * n2 + v2];
            if *slot == 0 {
                *slot = n_pts3;
                n_pts3 += 1;
                points.extend_from_slice(&[m2.points[2 * v2], yk, m2.points[2 * v2 + 1]]);
            }
        }
    }
    let vid = |v2: u32, k: usize| -> u32 { used_v[k * n2 + v2 as usize] };
    let pt3 = |v3: u32| -> [f64; 3] {
        let i = 3 * v3 as usize;
        [points[i], points[i + 1], points[i + 2]]
    };
    let det3 = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| -> f64 {
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let sub3 = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] };

    let mut cells: Vec<u32> = Vec::new();
    let mut cell_labels: Vec<i32> = Vec::new();
    let mut cell_regions: Vec<u32> = Vec::new();
    for k in 0..nslab {
        for ti in 0..ntri2 {
            let Some((did, reff)) = prism_info[k * ntri2 + ti] else { continue };
            let mut v = [m2.cell(ti)[0], m2.cell(ti)[1], m2.cell(ti)[2]];
            v.sort_unstable();
            let (b0, b1, b2) = (vid(v[0], k), vid(v[1], k), vid(v[2], k));
            let (t0, t1, t2) = (vid(v[0], k + 1), vid(v[1], k + 1), vid(v[2], k + 1));
            for mut tet in [[b0, b1, b2, t2], [b0, b1, t1, t2], [b0, t0, t1, t2]] {
                let (p0, p1, p2, p3) = (pt3(tet[0]), pt3(tet[1]), pt3(tet[2]), pt3(tet[3]));
                if det3(sub3(p1, p0), sub3(p2, p0), sub3(p3, p0)) < 0.0 {
                    tet.swap(1, 2);
                }
                cells.extend_from_slice(&tet);
                cell_labels.push(did);
                cell_regions.push(reff);
            }
        }
    }
    if cell_labels.is_empty() {
        return Err("the swept 3D model encloses no volume".into());
    }

    // Extract outer and material-interface facets and tag them via GeomDomain3::boundary_label
    let nc = cell_labels.len();
    let mut keyed: Vec<([u32; 3], u32, u8)> = Vec::with_capacity(nc * 4);
    for c in 0..nc {
        let t = &cells[4 * c..4 * c + 4];
        for skip in 0..4u8 {
            let mut f = [0u32; 3];
            let mut m = 0;
            for (q, &v) in t.iter().enumerate() {
                if q != skip as usize {
                    f[m] = v;
                    m += 1;
                }
            }
            f.sort_unstable();
            keyed.push((f, c as u32, skip));
        }
    }
    keyed.sort_unstable();

    let outward = |&(_, c, skip): &([u32; 3], u32, u8)| -> [u32; 3] {
        let t = &cells[4 * c as usize..4 * c as usize + 4];
        let mut f = [0u32; 3];
        let mut m = 0;
        for (q, &v) in t.iter().enumerate() {
            if q != skip as usize {
                f[m] = v;
                m += 1;
            }
        }
        let (a, b, cc, d) = (pt3(f[0]), pt3(f[1]), pt3(f[2]), pt3(t[skip as usize]));
        if det3(sub3(b, a), sub3(cc, a), sub3(d, a)) > 0.0 {
            f.swap(1, 2);
        }
        f
    };

    let mut facets: Vec<u32> = Vec::new();
    let mut facet_labels: Vec<i32> = Vec::new();
    let mut i = 0;
    while i < keyed.len() {
        let mut j = i + 1;
        while j < keyed.len() && keyed[j].0 == keyed[i].0 {
            j += 1;
        }
        if j - i == 1 {
            let tri = outward(&keyed[i]);
            let r1 = cell_regions[keyed[i].1 as usize];
            let (p0, p1, p2) = (pt3(tri[0]), pt3(tri[1]), pt3(tri[2]));
            let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
            facets.extend_from_slice(&tri);
            facet_labels.push(dom.boundary_label(cen, Some(r1), None));
        } else if j - i == 2 {
            let (c1, c2) = (keyed[i].1 as usize, keyed[i + 1].1 as usize);
            if cell_labels[c1] != cell_labels[c2] {
                let first = if cell_labels[c1] < cell_labels[c2] { i } else { i + 1 };
                let tri = outward(&keyed[first]);
                let (r1, r2) = (cell_regions[c1], cell_regions[c2]);
                let (p0, p1, p2) = (pt3(tri[0]), pt3(tri[1]), pt3(tri[2]));
                let cen = [(p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0, (p0[2] + p1[2] + p2[2]) / 3.0];
                let sides = if r1 != r2 { (Some(r1.min(r2)), Some(r1.max(r2))) } else { (Some(r1), None) };
                facets.extend_from_slice(&tri);
                facet_labels.push(dom.boundary_label(cen, sides.0, sides.1));
            }
        }
        i = j;
    }

    Ok(SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels))
}
