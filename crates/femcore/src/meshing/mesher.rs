//! The triangulation core shared by the geom path (`geomdrive`): chains of boundary points become a constrained Delaunay triangulation with
//! Ruppert refinement (via the `spade` crate), cells labelled by a domain function; and the sweep of a section along y into layered
//! tetrahedra (`extrude_y`).

use crate::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use spade::{AngleLimit, ConstrainedDelaunayTriangulation, Point2, RefinementParameters, Triangulation};
use std::collections::{HashMap, HashSet};

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
    /// the smallest size asked for anywhere (base spacing of the octave lattice grid)
    pub h_fine_min: f64,
    pub min_angle_deg: f64,
    /// bounding box of the chains
    pub lo: [f64; 2],
    pub hi: [f64; 2],
    /// cap on additional Steiner vertices as a multiple of initial vertices (`None`: `(20 * n_init).max(100_000)`)
    pub max_steiner_multiplier: Option<usize>,
    /// fail with `Err` if Ruppert refinement stops before completion or leaves unhealed interior slivers
    pub strict_mesh: bool,
    /// optional upper cap on 2D cell / additional vertex budget
    pub max_cells: Option<usize>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct MeshChainsReport {
    #[cfg_attr(not(test), allow(dead_code))]
    pub n_init: usize,
    #[cfg_attr(not(test), allow(dead_code))]
    pub added_vertices: usize,
    #[cfg_attr(not(test), allow(dead_code))]
    pub refinement_complete: bool,
    #[cfg_attr(not(test), allow(dead_code))]
    pub healed_caps: usize,
    #[cfg_attr(not(test), allow(dead_code))]
    pub dropped_slivers: usize,
    pub warnings: Vec<String>,
}

#[derive(Clone, Copy)]
struct SegEntry {
    a: [f64; 2],
    b: [f64; 2],
    part: usize,
    label: i32,
    orig_idx: usize,
    bbox: [f64; 4],
}

enum SegBvhKind {
    Leaf { start: usize, end: usize },
    Branch { left: usize, right: usize },
}

struct SegBvhNode {
    bbox: [f64; 4],
    kind: SegBvhKind,
}

/// Flat `Vec`-backed 2D AABB BVH over chain segments for $O(\log N)$ boundary proximity and facet labelling.
struct SegBvh2 {
    segs: Vec<SegEntry>,
    nodes: Vec<SegBvhNode>,
}

impl SegBvh2 {
    fn build(chains: &[Chain]) -> Self {
        let mut segs = Vec::new();
        for ch in chains {
            for w in ch.pts.windows(2) {
                let orig_idx = segs.len();
                let bbox = [w[0][0].min(w[1][0]), w[0][1].min(w[1][1]), w[0][0].max(w[1][0]), w[0][1].max(w[1][1])];
                segs.push(SegEntry { a: w[0], b: w[1], part: ch.part, label: ch.label, orig_idx, bbox });
            }
        }
        let mut nodes = Vec::with_capacity(segs.len().saturating_mul(2));
        if !segs.is_empty() {
            let n = segs.len();
            Self::build_rec(&mut segs, 0, n, &mut nodes);
        }
        Self { segs, nodes }
    }

    fn build_rec(segs: &mut [SegEntry], start: usize, end: usize, nodes: &mut Vec<SegBvhNode>) -> usize {
        let mut bbox = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
        let mut cbox = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
        for s in &segs[start..end] {
            bbox[0] = bbox[0].min(s.bbox[0]);
            bbox[1] = bbox[1].min(s.bbox[1]);
            bbox[2] = bbox[2].max(s.bbox[2]);
            bbox[3] = bbox[3].max(s.bbox[3]);
            let cx = 0.5 * (s.bbox[0] + s.bbox[2]);
            let cy = 0.5 * (s.bbox[1] + s.bbox[3]);
            cbox[0] = cbox[0].min(cx);
            cbox[1] = cbox[1].min(cy);
            cbox[2] = cbox[2].max(cx);
            cbox[3] = cbox[3].max(cy);
        }
        let idx = nodes.len();
        nodes.push(SegBvhNode { bbox, kind: SegBvhKind::Leaf { start, end } });
        if end - start <= 8 {
            return idx;
        }
        let axis = usize::from(cbox[3] - cbox[1] > cbox[2] - cbox[0]);
        let mid = start + (end - start) / 2;
        segs[start..end].select_nth_unstable_by(mid - start, |x, y| {
            let cx = 0.5 * (x.bbox[axis] + x.bbox[axis + 2]);
            let cy = 0.5 * (y.bbox[axis] + y.bbox[axis + 2]);
            cx.total_cmp(&cy).then_with(|| x.orig_idx.cmp(&y.orig_idx))
        });
        let left = Self::build_rec(segs, start, mid, nodes);
        let right = Self::build_rec(segs, mid, end, nodes);
        nodes[idx].kind = SegBvhKind::Branch { left, right };
        idx
    }

    #[inline]
    fn box_dist_sq(b: &[f64; 4], p: [f64; 2]) -> f64 {
        let dx = (b[0] - p[0]).max(p[0] - b[2]).max(0.0);
        let dy = (b[1] - p[1]).max(p[1] - b[3]).max(0.0);
        dx * dx + dy * dy
    }

    fn any_within(&self, p: [f64; 2], r: f64) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        let r2 = r * r;
        let mut stack = [0usize; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp]];
            if Self::box_dist_sq(&node.bbox, p) >= r2 {
                continue;
            }
            match node.kind {
                SegBvhKind::Leaf { start, end } => {
                    for s in &self.segs[start..end] {
                        if Self::box_dist_sq(&s.bbox, p) < r2 && dist_point_segment(p, s.a, s.b) < r {
                            return true;
                        }
                    }
                }
                SegBvhKind::Branch { left, right } => {
                    stack[sp] = right;
                    stack[sp + 1] = left;
                    sp += 2;
                }
            }
        }
        false
    }

    fn candidates_within(&self, p: [f64; 2], r: f64, out: &mut Vec<usize>) {
        out.clear();
        if self.nodes.is_empty() {
            return;
        }
        let r2 = r * r;
        let mut stack = [0usize; 64];
        let mut sp = 1usize;
        while sp > 0 {
            sp -= 1;
            let node = &self.nodes[stack[sp]];
            if Self::box_dist_sq(&node.bbox, p) >= r2 {
                continue;
            }
            match node.kind {
                SegBvhKind::Leaf { start, end } => {
                    for (i, s) in self.segs[start..end].iter().enumerate() {
                        if Self::box_dist_sq(&s.bbox, p) < r2 {
                            out.push(start + i);
                        }
                    }
                }
                SegBvhKind::Branch { left, right } => {
                    stack[sp] = right;
                    stack[sp + 1] = left;
                    sp += 2;
                }
            }
        }
        out.sort_unstable_by_key(|&idx| self.segs[idx].orig_idx);
    }
}

struct OctaveLevel {
    cell_size: f64,
    pts: HashMap<(i64, i64), Vec<[f64; 2]>>,
    has_finer: HashSet<(i64, i64)>,
}

/// Hierarchical octave-bucketed hash grid (`l_k = 2^k * 0.6 * h_fine_min`) bounding proximity probes to $\le 5 \times 5$
/// per active octave while preserving exact distance rejection `< 0.6 * h` across arbitrary size ratios.
struct OctaveFillGrid {
    l0: f64,
    h_fine_min: f64,
    levels: Vec<OctaveLevel>,
}

impl OctaveFillGrid {
    fn new(l0: f64, h_fine_min: f64, fills: &[Fill]) -> Self {
        let mut k_max = 0usize;
        for f in fills {
            let k = Self::octave_of_static(f.h, h_fine_min);
            k_max = k_max.max(k);
        }
        let levels = (0..=k_max)
            .map(|m| OctaveLevel {
                cell_size: l0 * ((1u64 << m) as f64),
                pts: HashMap::new(),
                has_finer: HashSet::new(),
            })
            .collect();
        Self { l0, h_fine_min, levels }
    }

    #[inline]
    fn octave_of_static(h: f64, h_fine_min: f64) -> usize {
        let ratio = (h / (h_fine_min * (1.0 + 1e-11))).max(1.0);
        (ratio.log2().floor() as usize).min(30)
    }

    #[inline]
    fn octave_of(&self, h: f64) -> usize {
        Self::octave_of_static(h, self.h_fine_min).min(self.levels.len().saturating_sub(1))
    }

    #[inline]
    fn base_cell(&self, x: f64) -> i64 {
        (x / self.l0).floor() as i64
    }

    fn insert(&mut self, p: [f64; 2], h: f64) {
        let k = self.octave_of(h);
        let cx0 = self.base_cell(p[0]);
        let cy0 = self.base_cell(p[1]);
        self.levels[k].pts.entry((cx0 >> k, cy0 >> k)).or_default().push(p);
        for m in (k + 1)..self.levels.len() {
            self.levels[m].has_finer.insert((cx0 >> m, cy0 >> m));
        }
    }

    fn near_any(&self, p: [f64; 2], r: f64, h: f64) -> bool {
        let k = self.octave_of(h);
        let r_pad = r * (1.0 + 1e-12);
        let bx_lo = self.base_cell(p[0] - r_pad);
        let bx_hi = self.base_cell(p[0] + r_pad);
        let by_lo = self.base_cell(p[1] - r_pad);
        let by_hi = self.base_cell(p[1] + r_pad);

        // 1. Active octaves `m >= k`: window `(cx_lo..=cx_hi) x (cy_lo..=cy_hi)` is guaranteed <= 5x5 cells
        for m in k..self.levels.len() {
            let lvl = &self.levels[m];
            if lvl.pts.is_empty() {
                continue;
            }
            let (cx_lo, cx_hi) = (bx_lo >> m, bx_hi >> m);
            let (cy_lo, cy_hi) = (by_lo >> m, by_hi >> m);
            debug_assert!(cx_hi - cx_lo <= 4 && cy_hi - cy_lo <= 4);
            for cx in cx_lo..=cx_hi {
                for cy in cy_lo..=cy_hi {
                    if let Some(bucket) = lvl.pts.get(&(cx, cy)) {
                        if bucket.iter().any(|q| (q[0] - p[0]).hypot(q[1] - p[1]) < r) {
                            return true;
                        }
                    }
                }
            }
        }

        // 2. Finer octaves `m < k`: prune empty / distant / wholly-contained coarse cells via `has_finer`
        if k > 0 && !self.levels[k].has_finer.is_empty() {
            let mut stack: Vec<(usize, i64, i64)> = Vec::with_capacity(32);
            let (cx_lo, cx_hi) = (bx_lo >> k, bx_hi >> k);
            let (cy_lo, cy_hi) = (by_lo >> k, by_hi >> k);
            for cx in cx_lo..=cx_hi {
                for cy in cy_lo..=cy_hi {
                    if self.levels[k].has_finer.contains(&(cx, cy)) {
                        stack.push((k, cx, cy));
                    }
                }
            }
            while let Some((m, cx, cy)) = stack.pop() {
                let lm = self.levels[m].cell_size;
                let (x0, x1) = (cx as f64 * lm, (cx + 1) as f64 * lm);
                let (y0, y1) = (cy as f64 * lm, (cy + 1) as f64 * lm);
                let pad = 1e-11 * x0.abs().max(x1.abs()).max(y0.abs()).max(y1.abs()).max(lm);
                let dx_min = (x0 - pad - p[0]).max(p[0] - (x1 + pad)).max(0.0);
                let dy_min = (y0 - pad - p[1]).max(p[1] - (y1 + pad)).max(0.0);
                if dx_min.hypot(dy_min) >= r {
                    continue;
                }
                let dx_max = (p[0] - (x0 - pad)).abs().max((p[0] - (x1 + pad)).abs());
                let dy_max = (p[1] - (y0 - pad)).abs().max((p[1] - (y1 + pad)).abs());
                if dx_max.hypot(dy_max) < r {
                    return true;
                }
                let child_m = m - 1;
                let child_lvl = &self.levels[child_m];
                for dcx in 0..=1 {
                    for dcy in 0..=1 {
                        let (ccx, ccy) = ((cx << 1) + dcx, (cy << 1) + dcy);
                        if let Some(bucket) = child_lvl.pts.get(&(ccx, ccy)) {
                            if bucket.iter().any(|q| (q[0] - p[0]).hypot(q[1] - p[1]) < r) {
                                return true;
                            }
                        }
                        if child_m > 0 && child_lvl.has_finer.contains(&(ccx, ccy)) {
                            stack.push((child_m, ccx, ccy));
                        }
                    }
                }
            }
        }
        false
    }
}

/// The shared tail of the meshers: constrained Delaunay triangulation of the chains (conforming, labelled), lattice points, Ruppert refinement,
/// cells labelled by `domain_at` (None: no cell) and boundary facets labelled by the chain they lie on.
pub(crate) fn mesh_chains(
    chains: &[Chain],
    fills: &[Fill],
    mp: &MeshParams,
    domain_at: &dyn Fn([f64; 2]) -> Option<i32>,
) -> Result<(SimplexMesh, MeshChainsReport), String> {
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
    let lcell = (0.6 * h_fine_min).max(1e-9 * scale);
    let seg_bvh = SegBvh2::build(chains);
    if !fills.is_empty() {
        let mut octave_grid = OctaveFillGrid::new(lcell, h_fine_min, fills);
        for fill in fills {
            let (h, blo, bhi) = (fill.h, fill.lo, fill.hi);
            let dy = h * 3f64.sqrt() / 2.0;
            let mut row = 0;
            let mut y = blo[1] + dy;
            while y < bhi[1] {
                let mut x = blo[0] + if row % 2 == 0 { 0.5 * h } else { h };
                while x < bhi[0] {
                    if (fill.owner)([x, y]) {
                        let near_boundary = seg_bvh.any_within([x, y], 0.4 * h);
                        if !near_boundary && !octave_grid.near_any([x, y], 0.6 * h, h) {
                            let pn = nrm([x, y]);
                            let _ = cdt.insert(Point2::new(pn[0], pn[1]));
                            octave_grid.insert([x, y], h);
                        }
                    }
                    x += h;
                }
                y += dy;
                row += 1;
            }
        }
    }
    let n_init = cdt.num_vertices();
    let max_add = match mp.max_steiner_multiplier {
        Some(mult) => n_init.saturating_mul(mult),
        None => n_init.saturating_mul(20).max(100_000),
    };
    let max_add = match mp.max_cells {
        Some(mc) => max_add.min(mc),
        None => max_add,
    };
    let params = RefinementParameters::<f64>::new()
        // no exclude_outer_faces: it decides inside/outside by constraint parity, which nested or
        // overlapping shapes break; unwanted cells are dropped by the shape test below instead
        .keep_constraint_edges()
        .with_angle_limit(AngleLimit::from_deg(mp.min_angle_deg.clamp(0.0, 33.0)))
        .with_max_allowed_area(0.5 * (h_global / scale) * (h_global / scale))
        .with_min_required_area(1e-14)
        .with_max_additional_vertices(max_add);
    let ref_res = cdt.refine(params);
    let added_vertices = cdt.num_vertices().saturating_sub(n_init);

    // cells: faces whose centroid is inside some shape (the last containing shape wins)
    let nvert = cdt.num_vertices();
    let weld_tol = 1e-9_f64;
    let weld_tol2 = weld_tol * weld_tol;
    let vpos: Vec<[f64; 2]> = (0..nvert)
        .map(|i| {
            let p = cdt.vertex(spade::handles::FixedVertexHandle::from_index(i)).position();
            [p.x, p.y]
        })
        .collect();
    let mut order: Vec<usize> = (0..nvert).collect();
    order.sort_by(|&a, &b| vpos[a][0].total_cmp(&vpos[b][0]).then_with(|| a.cmp(&b)));
    let mut canon: Vec<usize> = (0..nvert).collect();
    for ii in 0..nvert {
        let i = order[ii];
        let mut jj = ii;
        while jj > 0 {
            jj -= 1;
            let j = order[jj];
            let dx = vpos[i][0] - vpos[j][0];
            if dx > weld_tol {
                break;
            }
            let dy = vpos[i][1] - vpos[j][1];
            if dx * dx + dy * dy <= weld_tol2 {
                let (lo, hi) = (canon[i].min(canon[j]), canon[i].max(canon[j]));
                canon[i] = lo;
                canon[j] = lo;
                canon[hi] = lo;
            }
        }
    }
    for i in 0..nvert {
        let mut r = i;
        while canon[r] != r {
            r = canon[r];
        }
        canon[i] = r;
    }

    let mut cells: Vec<u32> = vec![];
    let mut cell_labels: Vec<i32> = vec![];
    let mut caps: Vec<(u32, u32, u32)> = vec![];
    for f in cdt.inner_faces() {
        let [a, b, c] = f.vertices();
        let (ia, ib, ic) = (canon[a.fix().index()], canon[b.fix().index()], canon[c.fix().index()]);
        if ia == ib || ib == ic || ia == ic {
            continue;
        }
        let pos = [vpos[ia], vpos[ib], vpos[ic]];
        // a face whose vertices are collinear to rounding is not a cell: Ruppert refinement can split a
        // boundary edge at a point a rounding error off the line, leaving such a sliver outside the shapes
        let area2 = ((pos[1][0] - pos[0][0]) * (pos[2][1] - pos[0][1]) - (pos[2][0] - pos[0][0]) * (pos[1][1] - pos[0][1])).abs();
        let l2 = [
            (pos[0][0] - pos[1][0]).powi(2) + (pos[0][1] - pos[1][1]).powi(2),
            (pos[1][0] - pos[2][0]).powi(2) + (pos[1][1] - pos[2][1]).powi(2),
            (pos[2][0] - pos[0][0]).powi(2) + (pos[2][1] - pos[0][1]).powi(2),
        ];
        let longest2 = l2[0].max(l2[1]).max(l2[2]);
        if area2 < 1e-8 * longest2 {
            let verts = [ia as u32, ib as u32, ic as u32];
            let k = if l2[0] >= l2[1] && l2[0] >= l2[2] {
                0
            } else if l2[1] >= l2[2] {
                1
            } else {
                2
            };
            caps.push((verts[k], verts[(k + 1) % 3], verts[(k + 2) % 3]));
            continue;
        }
        let cen = den([(pos[0][0] + pos[1][0] + pos[2][0]) / 3.0, (pos[0][1] + pos[1][1] + pos[2][1]) / 3.0]);
        if let Some(domain) = domain_at(cen) {
            cells.extend_from_slice(&[ia as u32, ib as u32, ic as u32]);
            cell_labels.push(domain);
        }
    }
    let mut healed_caps = 0usize;
    let mut dropped_slivers = 0usize;
    if !caps.is_empty() {
        let mut cap_healed = vec![false; caps.len()];
        // An interior collinear cap `(u, v, m)` (stopped by `min_required_area`) leaves `(u, m)`, `(m, v)` and `(u, v)`
        // mismatched across `u-v`; split the single kept triangle across `(u, v)` at `m` so no interior crack remains.
        for _pass in 0..2 {
            let mut e_cnt: HashMap<(u32, u32), u32> = HashMap::new();
            for t in cells.chunks(3) {
                for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    *e_cnt.entry((a.min(b), a.max(b))).or_insert(0) += 1;
                }
            }
            let mut any = false;
            for (ci_cap, &(u, v, m)) in caps.iter().enumerate() {
                if cap_healed[ci_cap] {
                    continue;
                }
                let e_uv = (u.min(v), u.max(v));
                let e_um = (u.min(m), u.max(m));
                let e_mv = (m.min(v), m.max(v));
                if e_cnt.get(&e_uv).copied() == Some(1)
                    && e_cnt.get(&e_um).copied() == Some(1)
                    && e_cnt.get(&e_mv).copied() == Some(1)
                {
                    if let Some(ci) = cells.chunks(3).position(|t| {
                        [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])].iter().any(|&(a, b)| (a.min(b), a.max(b)) == e_uv)
                    }) {
                        let t = [cells[3 * ci], cells[3 * ci + 1], cells[3 * ci + 2]];
                        let (d, p, q) = if (t[0].min(t[1]), t[0].max(t[1])) == e_uv {
                            (t[2], t[0], t[1])
                        } else if (t[1].min(t[2]), t[1].max(t[2])) == e_uv {
                            (t[0], t[1], t[2])
                        } else {
                            (t[1], t[2], t[0])
                        };
                        let dom = cell_labels[ci];
                        cells[3 * ci..3 * ci + 3].copy_from_slice(&[d, p, m]);
                        cells.extend_from_slice(&[d, m, q]);
                        cell_labels.push(dom);
                        *e_cnt.entry(e_uv).or_insert(0) = 0;
                        *e_cnt.entry(e_um).or_insert(0) += 1;
                        *e_cnt.entry(e_mv).or_insert(0) += 1;
                        *e_cnt.entry((d.min(m), d.max(m))).or_insert(0) += 2;
                        cap_healed[ci_cap] = true;
                        healed_caps += 1;
                        any = true;
                    }
                }
            }
            if !any {
                break;
            }
        }
        if cap_healed.iter().any(|&h| !h) {
            let mut post_e_cnt: HashMap<(u32, u32), u32> = HashMap::new();
            for t in cells.chunks(3) {
                for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    *post_e_cnt.entry((a.min(b), a.max(b))).or_insert(0) += 1;
                }
            }
            for (ci_cap, &(u, v, m)) in caps.iter().enumerate() {
                if cap_healed[ci_cap] {
                    continue;
                }
                let e_uv = (u.min(v), u.max(v));
                let e_um = (u.min(m), u.max(m));
                let e_mv = (m.min(v), m.max(v));
                let c_uv = post_e_cnt.get(&e_uv).copied().unwrap_or(0);
                let c_um = post_e_cnt.get(&e_um).copied().unwrap_or(0);
                let c_mv = post_e_cnt.get(&e_mv).copied().unwrap_or(0);
                if c_uv > 0 {
                    dropped_slivers += 1;
                } else if (c_um > 0 || c_mv > 0) && (c_um < 2 || c_mv < 2) {
                    let (pu, pv) = (den(vpos[u as usize]), den(vpos[v as usize]));
                    let (dx, dy) = (pv[0] - pu[0], pv[1] - pu[1]);
                    let len = dx.hypot(dy);
                    if len > 0.0 {
                        let (nx, ny) = (-dy / len, dx / len);
                        let eps = 1e-5 * len;
                        let (mx, my) = (0.5 * (pu[0] + pv[0]), 0.5 * (pu[1] + pv[1]));
                        if domain_at([mx + eps * nx, my + eps * ny]).is_some()
                            && domain_at([mx - eps * nx, my - eps * ny]).is_some()
                        {
                            dropped_slivers += 1;
                        }
                    }
                }
            }
        }
    }
    crate::progress::note(format!(
        "2D Ruppert refinement: n_init={n_init}, added_vertices={added_vertices}, complete={}, healed_caps={healed_caps}, dropped_slivers={dropped_slivers}",
        ref_res.refinement_complete
    ));
    let mut warnings = Vec::new();
    if !ref_res.refinement_complete {
        warnings.push(format!(
            "2D Ruppert refinement stopped before completion after adding {added_vertices} Steiner vertices (n_init={n_init}, max_additional_vertices={max_add}); increase max_steiner_multiplier, lower min_angle_deg, or increase h"
        ));
    }
    if dropped_slivers > 0 {
        warnings.push(format!(
            "2D triangulation dropped {dropped_slivers} unhealed collinear sliver triangle(s) (healed_caps={healed_caps})"
        ));
    }
    if mp.strict_mesh && (!ref_res.refinement_complete || dropped_slivers > 0) {
        return Err(format!(
            "strict_mesh check failed: refinement_complete={}, dropped_slivers={dropped_slivers} (n_init={n_init}, added_vertices={added_vertices}, max_additional_vertices={max_add})",
            ref_res.refinement_complete
        ));
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
    let search_r = tol * 10.0 + 1e-12;
    let mut cand_idxs: Vec<usize> = Vec::new();
    for (a, b) in boundary {
        let pa = [points[2 * a as usize], points[2 * a as usize + 1]];
        let pb = [points[2 * b as usize], points[2 * b as usize + 1]];
        let mid = [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0];
        // the chain segment closest to the midpoint; the last shape wins ties
        seg_bvh.candidates_within(mid, search_r, &mut cand_idxs);
        let mut best: Option<(f64, usize, i32)> = None;
        for &si in &cand_idxs {
            let s = &seg_bvh.segs[si];
            let d = dist_point_segment(mid, s.a, s.b);
            if d < search_r && best.is_none_or(|(bd, bp, _)| d < bd - tol || (d <= bd + tol && s.part >= bp)) {
                best = Some((d, s.part, s.label));
            }
        }
        let lab = best.map(|b| b.2).ok_or_else(|| format!("a boundary edge near ({:.4e}, {:.4e}) belongs to no shape edge", mid[0], mid[1]))?;
        facets.extend_from_slice(&[a, b]);
        facet_labels.push(lab);
    }
    let nc = cell_labels.len();
    debug_assert_eq!(nc * 3, cells.len());
    let report = MeshChainsReport {
        n_init,
        added_vertices,
        refinement_complete: ref_res.refinement_complete,
        healed_caps,
        dropped_slivers,
        warnings,
    };
    Ok((SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels), report))
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
                min_angle_deg,
                lo: [ov_bb.lo[0], ov_bb.lo[1]],
                hi: [ov_bb.hi[0], ov_bb.hi[1]],
                max_steiner_multiplier: None,
                strict_mesh: false,
                max_cells: None,
            },
            &|cen| {
                let c = ov.locate(cen)?;
                cell_active[c as usize].then_some(c as i32 + 1)
            },
        )?
        .0
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
                min_angle_deg,
                lo: [xs[0], zs[0]],
                hi: [*xs.last().unwrap(), *zs.last().unwrap()],
                max_steiner_multiplier: None,
                strict_mesh: false,
                max_cells: None,
            },
            &locate_grid,
        )?
        .0
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn octave_fill_grid_1000_to_1_size_ratio_under_250ms_and_deterministic() {
        let h_global = 0.25;
        let h_fine = 0.00025; // exact 1000:1 size ratio
        let h_mid = 0.01;

        let subdiv = |a: [f64; 2], b: [f64; 2], n: usize| -> Vec<[f64; 2]> {
            (0..=n)
                .map(|i| {
                    let t = i as f64 / n as f64;
                    [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]
                })
                .collect()
        };
        let mut chains = Vec::new();
        // Outer 1.0 x 1.0 box subdivided into 2,000 boundary segments (500 per side)
        for (i, (a, b)) in [
            ([0.0, 0.0], [1.0, 0.0]),
            ([1.0, 0.0], [1.0, 1.0]),
            ([1.0, 1.0], [0.0, 1.0]),
            ([0.0, 1.0], [0.0, 0.0]),
        ]
        .into_iter()
        .enumerate()
        {
            chains.push(Chain { pts: subdiv(a, b, 500), label: (i + 1) as i32, part: i });
        }
        // Fine inner box [0.499, 0.501]^2
        let (f0, f1) = (0.499, 0.501);
        for (i, (a, b)) in [
            ([f0, f0], [f1, f0]),
            ([f1, f0], [f1, f1]),
            ([f1, f1], [f0, f1]),
            ([f0, f1], [f0, f0]),
        ]
        .into_iter()
        .enumerate()
        {
            chains.push(Chain { pts: subdiv(a, b, 8), label: 10, part: 10 + i });
        }

        let make = || {
            let fills = vec![
                Fill {
                    lo: [f0, f0],
                    hi: [f1, f1],
                    h: h_fine,
                    owner: Box::new(move |p| p[0] > f0 && p[0] < f1 && p[1] > f0 && p[1] < f1),
                },
                Fill {
                    lo: [0.45, 0.45],
                    hi: [0.55, 0.55],
                    h: h_mid,
                    owner: Box::new(move |p| {
                        p[0] > 0.45
                            && p[0] < 0.55
                            && p[1] > 0.45
                            && p[1] < 0.55
                            && !(p[0] >= f0 && p[0] <= f1 && p[1] >= f0 && p[1] <= f1)
                    }),
                },
            ];
            let mp = MeshParams {
                h: h_global,
                h_fine_min: h_fine,
                min_angle_deg: 25.0,
                lo: [0.0, 0.0],
                hi: [1.0, 1.0],
                max_steiner_multiplier: None,
                strict_mesh: true,
                max_cells: None,
            };
            mesh_chains(&chains, &fills, &mp, &|p| {
                if p[0] < 0.0 || p[0] > 1.0 || p[1] < 0.0 || p[1] > 1.0 {
                    None
                } else if p[0] > f0 && p[0] < f1 && p[1] > f0 && p[1] < f1 {
                    Some(2)
                } else {
                    Some(1)
                }
            })
            .unwrap()
        };

        let t0 = std::time::Instant::now();
        let (m1, rep1) = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap().install(make);
        let dt_ms = t0.elapsed().as_secs_f64() * 1000.0;
        assert!(
            dt_ms < 250.0,
            "1000:1 octave fill grid + 2000-segment boundary took {dt_ms:.2} ms (limit 250 ms)"
        );
        assert!(rep1.n_init > 2000 && rep1.added_vertices > 0);
        assert!(rep1.refinement_complete, "Ruppert must complete on 1000:1 fill grid");
        assert_eq!(rep1.dropped_slivers, 0, "healed_caps={}, dropped={}", rep1.healed_caps, rep1.dropped_slivers);
        assert!(rep1.warnings.is_empty());

        let (m8, _) = rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap().install(make);
        assert_eq!(
            m1.points.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            m8.points.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
        assert_eq!(m1.cells, m8.cells);
        assert_eq!(m1.cell_labels, m8.cell_labels);
        assert_eq!(m1.facets, m8.facets);
        assert_eq!(m1.facet_labels, m8.facet_labels);
    }
}
