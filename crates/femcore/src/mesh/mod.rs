//! Simplex meshes (triangles in 2D, tetrahedra in 3D) with their boundary facets, the numbering of P1/P2 unknowns, and
//! structured generators with numbered boundaries.
//!
//! A `SimplexMesh` stores the vertices, the cells with a domain number each (1-based), the boundary facets (edges in 2D,
//! triangles in 3D) with a boundary number each (1-based), and all distinct edges. Unknown numbering: P1 unknowns are the
//! vertices; P2 unknowns are the vertices followed by one per edge (its midpoint), in the sorted order of the edges.
//! `rectangle` numbers its boundaries 1 x=0, 2 y=0, 3 y=ly, 4 x=lx; `block` 1 x=0, 2 y=0, 3 z=0, 4 z=lz, 5 y=ly,
//! 6 x=lx. Meshes given as lists of cells may contain quadrilaterals, hexahedra and prisms; `from_cells_split` cuts
//! them into simplices without new vertices (`SplitMode::Minimal`).

pub mod cells;
pub mod exchange;
pub mod fixtures;

pub mod io;

pub mod quality;
pub mod view;

pub use exchange::Mesh;
use crate::mesh::cells::CellMesh as NativeCells;
use std::collections::BTreeSet;

/// How hexahedra, prisms and boundary quadrilaterals are cut into simplices in 3D.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SplitMode {
    /// Cut along the face diagonal that starts at the lowest-numbered vertex of each quadrilateral
    /// face, and cone each cell from its lowest-numbered vertex (Dompierre et al.): 6 tetrahedra per
    /// hexahedron, 3 per prism, no new vertices. Neighbouring cells always agree on shared faces
    /// because the rule depends only on the global vertex numbers. This keeps the dof count of the
    /// original mesh and gives much better conditioned systems.
    #[default]
    Minimal,
    /// A vertex at the centre of every quadrilateral face and of every cell: 24 tetrahedra per
    /// hexahedron, 14 per prism. Symmetric, but about five times as many unknowns and poorer elements.
    Star,
}

/// Two triangles of a quadrilateral given in cyclic order, cut along the diagonal from its
/// lowest-numbered vertex.
fn quad_tris(cyc: [u32; 4]) -> [[u32; 3]; 2] {
    let m = (0..4).min_by_key(|&i| cyc[i]).unwrap();
    let c = |k: usize| cyc[(m + k) % 4];
    [[c(0), c(1), c(2)], [c(0), c(2), c(3)]]
}

/// What `from_cells` had to do to make a mesh usable by the simplex solvers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SplitReport {
    pub quads: usize,
    pub hexes: usize,
    pub prisms: usize,
    pub boundary_quads: usize,
    pub cells_after: usize,
}

impl SplitReport {
    pub fn is_empty(&self) -> bool {
        self.quads + self.hexes + self.prisms + self.boundary_quads == 0
    }
    pub fn describe(&self) -> String {
        let mut parts = vec![];
        if self.quads > 0 {
            parts.push(format!("{} quadrilaterals into triangles", self.quads));
        }
        if self.hexes > 0 {
            parts.push(format!("{} hexahedra into tetrahedra", self.hexes));
        }
        if self.prisms > 0 {
            parts.push(format!("{} prisms into tetrahedra", self.prisms));
        }
        format!("split {} ({} simplex cells in all)", parts.join(", "), self.cells_after)
    }
}

/// The mesh type the physics crates take (`nyx_physics::Mesh`). Today it is the simplex mesh; restructure phase 2 makes it the one native
/// per-cell-type mesh (ADR 5). Code outside femcore names the mesh through this alias so that change is one line.

#[derive(Clone, Debug)]
pub struct SimplexMesh {
    pub dim: usize,
    /// flat coordinates, `dim` per vertex
    pub points: Vec<f64>,
    /// flat vertex indices, `dim + 1` per cell
    pub cells: Vec<u32>,
    /// domain number of each cell (1-based)
    pub cell_labels: Vec<i32>,
    /// flat vertex indices of boundary facets (edges in 2D, triangles in 3D), `dim` each
    pub facets: Vec<u32>,
    /// boundary number of each facet (1-based)
    pub facet_labels: Vec<i32>,
    /// all distinct edges as sorted vertex pairs, in sorted order
    edges: Vec<(u32, u32)>,
    /// per cell, the edge index of each local vertex pair (i<j in lexicographic order)
    cell_edges: Vec<u32>,
    /// 3D: the curved edges (ADR 82), as (edge index, position of the mid-edge geometry node), sorted by edge index; an edge that is not
    /// listed is straight. Empty for every mesh made by `new`.
    edge_nodes: Vec<(u32, [f64; 3])>,
}

/// The distinct edges of a mesh as sorted vertex pairs, and for every cell the edge index of each
/// local vertex pair (i < j in lexicographic order). Built by bucketing the pairs by their lower
/// vertex and sorting the small buckets in parallel, so large meshes take seconds, not minutes.
fn build_edges(dim: usize, cells: &[u32], nv: usize) -> (Vec<(u32, u32)>, Vec<u32>) {
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicU32, Ordering::Relaxed};
    let nloc = dim + 1;
    let nc = cells.len() / nloc;
    let npair = nloc * (nloc - 1) / 2;
    let nt = if nc < 20_000 { 1 } else { crate::par::threads() };
    let cchunk = nc.div_ceil(nt).max(1);
    let pairs_of = |c: usize, f: &mut dyn FnMut(u32, u32)| {
        let v = &cells[c * nloc..(c + 1) * nloc];
        for i in 0..nloc {
            for j in i + 1..nloc {
                f(v[i].min(v[j]), v[i].max(v[j]));
            }
        }
    };
    // 1. count the pairs per lower vertex
    let counts: Vec<AtomicU32> = (0..nv + 1).map(|_| AtomicU32::new(0)).collect();
    (0..nt).into_par_iter().for_each(|t| {
        for c in t * cchunk..((t + 1) * cchunk).min(nc) {
            pairs_of(c, &mut |a, _| {
                counts[a as usize + 1].fetch_add(1, Relaxed);
            });
        }
    });
    let mut start = vec![0usize; nv + 1];
    for v in 0..nv {
        start[v + 1] = start[v] + counts[v + 1].load(Relaxed) as usize;
    }
    // 2. scatter the upper vertices into their buckets
    let cursor: Vec<AtomicU32> = (0..nv).map(|_| AtomicU32::new(0)).collect();
    let buf: Vec<AtomicU32> = (0..start[nv]).map(|_| AtomicU32::new(0)).collect();
    (0..nt).into_par_iter().for_each(|t| {
        for c in t * cchunk..((t + 1) * cchunk).min(nc) {
            pairs_of(c, &mut |a, b| {
                let pos = start[a as usize] + cursor[a as usize].fetch_add(1, Relaxed) as usize;
                buf[pos].store(b, Relaxed);
            });
        }
    });
    let mut buf: Vec<u32> = buf.into_iter().map(|x| x.into_inner()).collect();
    // 3. sort and deduplicate every bucket, by blocks of vertices
    let vchunk = nv.div_ceil(nt).max(1);
    let mut slices: Vec<(usize, &mut [u32])> = vec![];
    {
        let mut rest: &mut [u32] = &mut buf;
        let mut v0 = 0;
        while v0 < nv {
            let v1 = (v0 + vchunk).min(nv);
            let (head, tail) = rest.split_at_mut(start[v1] - start[v0]);
            slices.push((v0, head));
            rest = tail;
            v0 = v1;
        }
    }
    type EdgeBlock = (Vec<(u32, u32)>, Vec<usize>);
    let blocks: Vec<EdgeBlock> = slices
        .into_par_iter()
        .map(|(v0, slice)| {
            let v1 = (v0 + vchunk).min(nv);
            let base = start[v0];
            let mut edges = Vec::new();
            let mut unique = Vec::with_capacity(v1 - v0);
            for v in v0..v1 {
                let bucket = &mut slice[start[v] - base..start[v + 1] - base];
                bucket.sort_unstable();
                let before = edges.len();
                let mut last = u32::MAX;
                for &w in bucket.iter() {
                    if w != last {
                        edges.push((v as u32, w));
                        last = w;
                    }
                }
                unique.push(edges.len() - before);
            }
            (edges, unique)
        })
        .collect();
    let total: usize = blocks.iter().map(|b| b.0.len()).sum();
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(total);
    let mut estart = vec![0usize; nv + 1];
    let mut v = 0;
    for (e, unique) in blocks {
        edges.extend(e);
        for u in unique {
            estart[v + 1] = estart[v] + u;
            v += 1;
        }
    }
    // 4. the edge index of every cell pair
    let mut cell_edges = vec![0u32; nc * npair];
    cell_edges.par_chunks_mut(cchunk * npair).enumerate().for_each(|(t, out)| {
        for (lc, o) in out.chunks_mut(npair).enumerate() {
            let c = t * cchunk + lc;
            let v = &cells[c * nloc..(c + 1) * nloc];
            let mut k = 0;
            for i in 0..nloc {
                for j in i + 1..nloc {
                    let (a, b) = (v[i].min(v[j]), v[i].max(v[j]));
                    let range = &edges[estart[a as usize]..estart[a as usize + 1]];
                    let pos = range.binary_search_by_key(&b, |e| e.1).expect("edge present");
                    o[k] = (estart[a as usize] + pos) as u32;
                    k += 1;
                }
            }
        }
    });
    (edges, cell_edges)
}

impl SimplexMesh {
    pub fn new(dim: usize, points: Vec<f64>, cells: Vec<u32>, cell_labels: Vec<i32>, facets: Vec<u32>, facet_labels: Vec<i32>) -> Self {
        assert!(dim == 2 || dim == 3);
        let nv = points.len() / dim;
        let (edges, cell_edges) = build_edges(dim, &cells, nv);
        SimplexMesh { dim, points, cells, cell_labels, facets, facet_labels, edges, cell_edges, edge_nodes: Vec::new() }
    }

    /// Sets the curved edges of a 3D mesh: `(vertex pair, node position)`; pairs that are not an edge of the mesh are ignored (the
    /// connectivity is what decides), and the list is kept sorted by edge index.
    pub fn set_edge_nodes(&mut self, nodes: &[([u32; 2], [f64; 3])]) {
        if self.dim != 3 {
            self.edge_nodes.clear();
            return;
        }
        let mut v: Vec<(u32, [f64; 3])> = nodes
            .iter()
            .filter_map(|&([a, b], p)| self.edges.binary_search(&(a.min(b), a.max(b))).ok().map(|e| (e as u32, p)))
            .collect();
        v.sort_by_key(|e| e.0);
        v.dedup_by_key(|e| e.0);
        self.edge_nodes = v;
    }
    /// The curved edges as `(sorted vertex pair, node position)`, sorted by pair (the exchange form of `CellMesh::edge_nodes`).
    pub fn edge_node_list(&self) -> Vec<([u32; 2], [f64; 3])> {
        let mut v: Vec<([u32; 2], [f64; 3])> = self.edge_nodes.iter().map(|&(e, p)| ([self.edges[e as usize].0, self.edges[e as usize].1], p)).collect();
        v.sort_by_key(|e| e.0);
        v
    }
    /// Position of the mid-edge geometry node of edge `e` (index into `edges()`); `None` when the edge is straight.
    pub fn edge_node(&self, e: usize) -> Option<[f64; 3]> {
        self.edge_nodes.binary_search_by_key(&(e as u32), |x| x.0).ok().map(|i| self.edge_nodes[i].1)
    }
    /// Whether the mesh has any curved edge.
    pub fn has_curved(&self) -> bool {
        !self.edge_nodes.is_empty()
    }
    /// Whether one of the edges of cell `c` is curved (3D only).
    pub fn cell_is_curved(&self, c: usize) -> bool {
        !self.edge_nodes.is_empty() && self.dim == 3 && self.cell_edge_ids(c).iter().any(|&e| self.edge_node(e as usize).is_some())
    }
    /// The ten geometry nodes of tetrahedron `c`: its four vertices, then the nodes of its six edges in the order of the local vertex
    /// pairs (0,1), (0,2), (0,3), (1,2), (1,3), (2,3); a straight edge gives its midpoint.
    pub fn cell_geometry_nodes(&self, c: usize) -> [[f64; 3]; 10] {
        let v = self.cell(c);
        let p = |i: usize| -> [f64; 3] { [self.points[v[i] as usize * 3], self.points[v[i] as usize * 3 + 1], self.points[v[i] as usize * 3 + 2]] };
        let mut out = [[0.0; 3]; 10];
        for (i, o) in out.iter_mut().enumerate().take(4) {
            *o = p(i);
        }
        let ids = self.cell_edge_ids(c);
        let mut k = 4;
        for i in 0..4 {
            for j in i + 1..4 {
                out[k] = self.edge_node(ids[k - 4] as usize).unwrap_or_else(|| {
                    let (a, b) = (p(i), p(j));
                    [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0]
                });
                k += 1;
            }
        }
        out
    }
    pub fn n_points(&self) -> usize {
        self.points.len() / self.dim
    }
    pub fn n_cells(&self) -> usize {
        self.cells.len() / (self.dim + 1)
    }
    pub fn n_facets(&self) -> usize {
        self.facets.len() / self.dim
    }
    /// Which boundary facets lie on the outside of the mesh (used by exactly one cell). Meshes
    /// from files and from the shape mesher also carry facets on interfaces between domains, which
    /// are not walls.
    pub fn exterior_facets(&self) -> Vec<bool> {
        use std::collections::HashMap;
        let nloc = self.dim + 1;
        let key = |v: &[u32]| -> [u32; 3] {
            let mut k = [u32::MAX; 3];
            for (i, &x) in v.iter().enumerate() {
                k[i] = x;
            }
            k[..v.len()].sort_unstable();
            k
        };
        let mut count: HashMap<[u32; 3], u32> = HashMap::new();
        for c in 0..self.n_cells() {
            let v = self.cell(c);
            for skip in 0..nloc {
                let f: Vec<u32> = (0..nloc).filter(|&i| i != skip).map(|i| v[i]).collect();
                *count.entry(key(&f)).or_insert(0) += 1;
            }
        }
        (0..self.n_facets()).map(|f| count.get(&key(self.facet(f))).copied().unwrap_or(0) == 1).collect()
    }

    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }
    pub fn edges(&self) -> &[(u32, u32)] {
        &self.edges
    }
    pub fn point(&self, i: usize) -> &[f64] {
        &self.points[i * self.dim..(i + 1) * self.dim]
    }
    pub fn cell(&self, c: usize) -> &[u32] {
        &self.cells[c * (self.dim + 1)..(c + 1) * (self.dim + 1)]
    }
    pub fn facet(&self, f: usize) -> &[u32] {
        &self.facets[f * self.dim..(f + 1) * self.dim]
    }
    /// Local pair edges of cell `c` (index into `edges()`), pairs (0,1),(0,2),...
    pub fn cell_edge_ids(&self, c: usize) -> &[u32] {
        let npair = self.dim * (self.dim + 1) / 2;
        &self.cell_edges[c * npair..(c + 1) * npair]
    }
    pub fn edge_index(&self, a: u32, b: u32) -> usize {
        self.edges.binary_search(&(a.min(b), a.max(b))).expect("edge in mesh")
    }
    pub fn n_dofs(&self, order: usize) -> usize {
        self.n_points() + if order == 2 { self.n_edges() } else { 0 }
    }
    pub fn nloc(&self, order: usize) -> usize {
        let n = self.dim + 1;
        if order == 1 { n } else { n + n * (n - 1) / 2 }
    }
    /// Dof indices of cell `c` (P1: vertices; P2: vertices then edge midpoints in the order of
    /// local vertex pairs).
    pub fn cell_dofs(&self, c: usize, order: usize, out: &mut Vec<usize>) {
        out.clear();
        out.extend(self.cell(c).iter().map(|&v| v as usize));
        if order == 2 {
            let nv = self.n_points();
            out.extend(self.cell_edge_ids(c).iter().map(|&e| nv + e as usize));
        }
    }
    /// All cell dofs, flat with stride `nloc(order)`.
    pub fn all_cell_dofs(&self, order: usize) -> Vec<usize> {
        let nloc = self.nloc(order);
        let mut v = Vec::with_capacity(self.n_cells() * nloc);
        let mut tmp = Vec::new();
        for c in 0..self.n_cells() {
            self.cell_dofs(c, order, &mut tmp);
            v.extend_from_slice(&tmp);
        }
        v
    }
    /// Dofs of boundary facet `f` (same local order as for cells).
    pub fn facet_dofs(&self, f: usize, order: usize, out: &mut Vec<usize>) {
        out.clear();
        let fv = self.facet(f);
        out.extend(fv.iter().map(|&v| v as usize));
        if order == 2 {
            let nv = self.n_points();
            for i in 0..fv.len() {
                for j in i + 1..fv.len() {
                    out.push(nv + self.edge_index(fv[i], fv[j]));
                }
            }
        }
    }
    /// Coordinates of every dof (vertices, then edge midpoints).
    pub fn dof_coordinates(&self, order: usize) -> Vec<f64> {
        let d = self.dim;
        let mut x = self.points.clone();
        if order == 2 {
            for &(a, b) in &self.edges {
                for k in 0..d {
                    x.push(0.5 * (self.points[a as usize * d + k] + self.points[b as usize * d + k]));
                }
            }
        }
        x
    }
    pub fn domain_labels(&self) -> Vec<i32> {
        self.cell_labels.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
    }
    pub fn boundary_labels(&self) -> Vec<i32> {
        self.facet_labels.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
    }
    pub fn bbox(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in self.points.chunks(self.dim) {
            for k in 0..self.dim {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        for k in self.dim..3 {
            lo[k] = 0.0;
            hi[k] = 0.0;
        }
        (lo, hi)
    }

    /// Relabels cells whose centroid lies in a box and boundary facets whose midpoint lies in a
    /// box. A box is `[min..., max...]` with `dim` values each. Later entries win.
    pub fn relabel(&mut self, domains: &[(Vec<f64>, i32)], boundaries: &[(Vec<f64>, i32)]) -> Result<(), String> {
        let d = self.dim;
        let inside = |p: &[f64], b: &[f64]| (0..d).all(|k| p[k] >= b[k] - 1e-12 && p[k] <= b[d + k] + 1e-12);
        for (b, _) in domains.iter().chain(boundaries.iter()) {
            if b.len() != 2 * d {
                return Err(format!("a box needs {} numbers (min and max per axis) for a {d}D mesh", 2 * d));
            }
        }
        for (b, lab) in domains {
            for c in 0..self.n_cells() {
                let cen = self.centroid(self.cell(c));
                if inside(&cen, b) {
                    self.cell_labels[c] = *lab;
                }
            }
        }
        for (b, lab) in boundaries {
            for f in 0..self.n_facets() {
                let cen = self.centroid(self.facet(f));
                if inside(&cen, b) {
                    self.facet_labels[f] = *lab;
                }
            }
        }
        Ok(())
    }

    fn centroid(&self, v: &[u32]) -> Vec<f64> {
        let d = self.dim;
        let mut c = vec![0.0; d];
        for &i in v {
            for k in 0..d {
                c[k] += self.points[i as usize * d + k];
            }
        }
        c.iter_mut().for_each(|x| *x /= v.len() as f64);
        c
    }

    /// Mesh from a mixed-cell mesh (`mesh::cells::CellMesh`). Domain and boundary numbers are 1-based in both and stay as they are.
    ///
    /// Cells that are not simplices are split into simplices so the solvers can use them:
    /// quadrilaterals into two triangles (along the shorter diagonal), and in 3D hexahedra and
    /// prisms into tetrahedra (`SplitMode`). The original vertices keep their indices; new vertices are appended.
    /// Straight-sided cells are reproduced exactly. Cells come out grouped by shape in the order of first appearance
    /// (cells, then boundary facets), each group in the order of the input.
    pub fn from_cells<const D: usize>(nm: &NativeCells<D>) -> Result<Self, String> {
        Self::from_cells_split(nm).map(|(m, _)| m)
    }

    pub fn from_cells_split<const D: usize>(nm: &NativeCells<D>) -> Result<(Self, SplitReport), String> {
        Self::from_cells_split_with(nm, SplitMode::default())
    }

    pub fn from_cells_split_with<const D: usize>(nm: &NativeCells<D>, mode: SplitMode) -> Result<(Self, SplitReport), String> {
        use crate::mesh::cells::CellType;
        let d = D;
        if d != 2 && d != 3 {
            return Err(format!("mesh dimension {d} is not supported"));
        }
        let mut report = SplitReport::default();
        let mut points: Vec<f64> = nm.points.iter().flat_map(|p| p.iter().copied()).collect();
        let (mut cells, mut cell_labels) = (vec![], vec![]);
        let (mut facets, mut facet_labels) = (vec![], vec![]);
        let dist2 = |pts: &[f64], a: u32, b: u32| -> f64 { (0..d).map(|k| (pts[a as usize * d + k] - pts[b as usize * d + k]).powi(2)).sum() };
        let cell_kind = |c: usize| match nm.types[c] {
            CellType::Tri => "tri",
            CellType::Quad => "quad",
            CellType::Tet => "tet",
            CellType::Prism => "prism",
            CellType::Hex => "hex",
        };
        let facet_kind = |f: usize| match nm.facet(f).len() {
            2 => "edg",
            3 => "tri",
            _ => "quad",
        };
        // shapes in the order of first appearance, cells first
        let mut kinds: Vec<&str> = vec![];
        for k in (0..nm.n_cells()).map(cell_kind).chain((0..nm.n_facets()).map(facet_kind)) {
            if !kinds.contains(&k) {
                kinds.push(k);
            }
        }
        // the quadrilateral and hexahedron vertices in tensor-product order (i + 2 j + 4 k), which the cuts below use
        let tensor_quad = |v: &[u32]| [v[0], v[1], v[3], v[2]];
        let tensor_hex = |v: &[u32]| [v[0], v[1], v[3], v[2], v[4], v[5], v[7], v[6]];
        let mut centers: std::collections::HashMap<[u32; 4], u32> = std::collections::HashMap::new();
        let mut center_of = |points: &mut Vec<f64>, q: [u32; 4]| -> u32 {
            let mut key = q;
            key.sort_unstable();
            if let Some(&c) = centers.get(&key) {
                return c;
            }
            let id = (points.len() / 3) as u32;
            let mut c = [0.0; 3];
            for &v in &q {
                for k in 0..3 {
                    c[k] += points[v as usize * 3 + k] / 4.0;
                }
            }
            points.extend_from_slice(&c);
            centers.insert(key, id);
            id
        };
        // faces of the reference cells in tensor-product vertex order [p00, p10, p01, p11]
        const HEX_FACES: [[usize; 4]; 6] = [[0, 2, 4, 6], [1, 3, 5, 7], [0, 1, 4, 5], [2, 3, 6, 7], [0, 1, 2, 3], [4, 5, 6, 7]];
        const PRISM_QUADS: [[usize; 4]; 3] = [[0, 1, 3, 4], [1, 2, 4, 5], [0, 2, 3, 5]];
        for kind in kinds {
            match (d, kind) {
                (2, "tri") | (3, "tet") => {
                    for c in (0..nm.n_cells()).filter(|&c| cell_kind(c) == kind) {
                        cells.extend_from_slice(nm.cell(c));
                        cell_labels.push(nm.label[c]);
                    }
                }
                (2, "quad") => {
                    for c in (0..nm.n_cells()).filter(|&c| cell_kind(c) == kind) {
                        let [a, b, cc, e] = tensor_quad(nm.cell(c));
                        let tris = if dist2(&points, a, e) <= dist2(&points, b, cc) { [[a, b, e], [a, e, cc]] } else { [[a, b, cc], [b, e, cc]] };
                        for tr in tris {
                            cells.extend_from_slice(&tr);
                            cell_labels.push(nm.label[c]);
                        }
                        report.quads += 1;
                    }
                }
                (2, "edg") | (3, "tri") => {
                    for f in (0..nm.n_facets()).filter(|&f| facet_kind(f) == kind) {
                        facets.extend_from_slice(nm.facet(f));
                        facet_labels.push(nm.facet_label[f]);
                    }
                }
                (3, "quad") => {
                    for f in (0..nm.n_facets()).filter(|&f| facet_kind(f) == kind) {
                        let q = tensor_quad(nm.facet(f));
                        let cyc = [q[0], q[1], q[3], q[2]];
                        if mode == SplitMode::Minimal {
                            for tr in quad_tris(cyc) {
                                facets.extend_from_slice(&tr);
                                facet_labels.push(nm.facet_label[f]);
                            }
                        } else {
                            let c = center_of(&mut points, q);
                            for k in 0..4 {
                                facets.extend_from_slice(&[c, cyc[k], cyc[(k + 1) % 4]]);
                                facet_labels.push(nm.facet_label[f]);
                            }
                        }
                        report.boundary_quads += 1;
                    }
                }
                (3, "hex") | (3, "prism") => {
                    let hex = kind == "hex";
                    for c in (0..nm.n_cells()).filter(|&c| cell_kind(c) == kind) {
                        let rv: Vec<u32> = if hex { tensor_hex(nm.cell(c)).to_vec() } else { nm.cell(c).to_vec() };
                        let r = &rv[..];
                        let label = nm.label[c];
                        if mode == SplitMode::Minimal {
                            // cone from the lowest-numbered vertex over the faces that do not contain it
                            let v0 = (0..r.len()).min_by_key(|&k| r[k]).unwrap();
                            let mut far: Vec<[u32; 3]> = vec![];
                            let quad_faces: Vec<[usize; 4]> = if hex { HEX_FACES.to_vec() } else { PRISM_QUADS.to_vec() };
                            for f in &quad_faces {
                                if !f.contains(&v0) {
                                    far.extend(quad_tris([r[f[0]], r[f[1]], r[f[3]], r[f[2]]]));
                                }
                            }
                            if !hex {
                                for tri in [[0usize, 1, 2], [3, 4, 5]] {
                                    if !tri.contains(&v0) {
                                        far.push([r[tri[0]], r[tri[1]], r[tri[2]]]);
                                    }
                                }
                            }
                            for tr in far {
                                cells.extend_from_slice(&[r[v0], tr[0], tr[1], tr[2]]);
                                cell_labels.push(label);
                            }
                            if hex {
                                report.hexes += 1;
                            } else {
                                report.prisms += 1;
                            }
                            continue;
                        }
                        let nv = r.len();
                        let mut cc = [0.0; 3];
                        for &v in r {
                            for k in 0..3 {
                                cc[k] += points[v as usize * 3 + k] / nv as f64;
                            }
                        }
                        let cid = (points.len() / 3) as u32;
                        points.extend_from_slice(&cc);
                        let mut tris: Vec<[u32; 3]> = vec![];
                        let quad_faces: Vec<[usize; 4]> = if hex { HEX_FACES.to_vec() } else { PRISM_QUADS.to_vec() };
                        for f in &quad_faces {
                            let q = [r[f[0]], r[f[1]], r[f[2]], r[f[3]]];
                            let c = center_of(&mut points, q);
                            let cyc = [q[0], q[1], q[3], q[2]];
                            for k in 0..4 {
                                tris.push([c, cyc[k], cyc[(k + 1) % 4]]);
                            }
                        }
                        if !hex {
                            tris.push([r[0], r[1], r[2]]);
                            tris.push([r[3], r[4], r[5]]);
                        }
                        for tr in tris {
                            cells.extend_from_slice(&[cid, tr[0], tr[1], tr[2]]);
                            cell_labels.push(label);
                        }
                        if hex {
                            report.hexes += 1;
                        } else {
                            report.prisms += 1;
                        }
                    }
                }
                (_, other) => return Err(format!("unsupported element shape '{other}' in a {d}D mesh")),
            }
        }
        if cells.is_empty() || facets.is_empty() {
            return Err("the mesh has no cells or no boundary elements".into());
        }
        report.cells_after = cells.len() / (d + 1);
        let mut mesh = SimplexMesh::new(d, points, cells, cell_labels, facets, facet_labels);
        // curved edges (ADR 82) survive only where the cells are the cells that were given: a mesh that had to be split keeps none
        if d == 3 && !nm.edge_nodes.is_empty() && report.is_empty() {
            let nodes: Vec<([u32; 2], [f64; 3])> = nm.edge_nodes.iter().map(|e| (e.0, [e.1[0], e.1[1], e.1[2]])).collect();
            mesh.set_edge_nodes(&nodes);
        }
        Ok((mesh, report))
    }
}

