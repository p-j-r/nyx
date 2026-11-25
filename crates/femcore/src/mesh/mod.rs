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
pub mod mixed;

pub mod io;

pub mod quality;
pub mod view;

pub use exchange::Mesh;
pub use mixed::{FacetType, MixedDiscretisation};
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

/// Structural, geometric, or edge-endpoint defect detected by [`SimplexMesh::check`].
#[derive(Clone, Debug, PartialEq)]
pub enum MeshIntegrityError {
    InvalidDimension(usize),
    LengthMismatch(&'static str),
    NonFiniteCoordinate { point: usize, axis: usize },
    VertexOutOfBounds { context: &'static str, index: usize, vertex: u32, n_points: usize },
    DuplicateVertexInSimplex { context: &'static str, index: usize, vertex: u32 },
    NonPositiveMeasure { cell: usize, det: f64 },
    NonManifoldFacet { vertices: Vec<u32>, incident_cells: usize },
    InconsistentFacetTag { vertices: Vec<u32>, detail: String },
    FacetNotInAnyCell { facet: usize, vertices: Vec<u32> },
    InvalidCurvedEdge { edge_index: u32, detail: String },
    InvalidQpEdge { edge_index: u32, detail: String },
}

impl std::fmt::Display for MeshIntegrityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDimension(d) => write!(f, "unsupported mesh dimension {d}"),
            Self::LengthMismatch(msg) => write!(f, "mesh array length mismatch: {msg}"),
            Self::NonFiniteCoordinate { point, axis } => {
                write!(f, "non-finite coordinate at point {point}, axis {axis}")
            }
            Self::VertexOutOfBounds { context, index, vertex, n_points } => {
                write!(f, "{context} {index}: vertex {vertex} out of bounds (n_points = {n_points})")
            }
            Self::DuplicateVertexInSimplex { context, index, vertex } => {
                write!(f, "{context} {index}: duplicate vertex {vertex}")
            }
            Self::NonPositiveMeasure { cell, det } => {
                write!(f, "cell {cell}: non-positive orientation determinant {det:e}")
            }
            Self::NonManifoldFacet { vertices, incident_cells } => {
                write!(f, "non-manifold facet {vertices:?} shared by {incident_cells} cells")
            }
            Self::InconsistentFacetTag { vertices, detail } => {
                write!(f, "inconsistent facet tag on {vertices:?}: {detail}")
            }
            Self::FacetNotInAnyCell { facet, vertices } => {
                write!(f, "facet {facet} ({vertices:?}) does not belong to any cell")
            }
            Self::InvalidCurvedEdge { edge_index, detail } => {
                write!(f, "invalid curved edge {edge_index}: {detail}")
            }
            Self::InvalidQpEdge { edge_index, detail } => {
                write!(f, "invalid quarter-point edge {edge_index}: {detail}")
            }
        }
    }
}

impl std::error::Error for MeshIntegrityError {}

impl From<MeshIntegrityError> for String {
    fn from(e: MeshIntegrityError) -> String {
        e.to_string()
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
    /// 3D: sorted edge indices carrying Barsoum quarter-point singularity nodes (`AF-2`); empty by default.
    qp_edges: Vec<u32>,
}

/// The distinct edges of a mesh as sorted vertex pairs, and for every cell the edge index of each
/// local vertex pair (i < j in lexicographic order). Built by bucketing the pairs by their lower
/// vertex and sorting the small buckets in parallel, so large meshes take seconds, not minutes.
fn build_edges(dim: usize, cells: &[u32], nv: usize) -> (Vec<(u32, u32)>, Vec<u32>) {
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicU32, Ordering::Relaxed};
    let nloc = dim + 1;
    let npair = nloc * (nloc - 1) / 2;
    if nv == 0 || !cells.len().is_multiple_of(nloc) || cells.iter().any(|&v| (v as usize) >= nv) {
        return (Vec::new(), Vec::new());
    }
    let nc = cells.len() / nloc;
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
        SimplexMesh { dim, points, cells, cell_labels, facets, facet_labels, edges, cell_edges, edge_nodes: Vec::new(), qp_edges: Vec::new() }
    }

    /// Verifies mesh coordinate finiteness, vertex index bounds, distinct simplex vertices, strictly
    /// positive cell orientation/measure (`det > 0.0`), manifold facet incidence (`<= 2` incident
    /// cells), consistent sheet/boundary tags, facet vertex membership, and curved/QP edge endpoints.
    pub fn check(&self) -> Result<(), MeshIntegrityError> {
        use std::collections::HashMap;
        if self.dim != 2 && self.dim != 3 {
            return Err(MeshIntegrityError::InvalidDimension(self.dim));
        }
        if !self.points.len().is_multiple_of(self.dim) {
            return Err(MeshIntegrityError::LengthMismatch("points length is not a multiple of dim"));
        }
        let n_verts = self.points.len() / self.dim;
        let nloc = self.dim + 1;
        if !self.cells.len().is_multiple_of(nloc) {
            return Err(MeshIntegrityError::LengthMismatch("cells length is not a multiple of dim + 1"));
        }
        let n_cells = self.cells.len() / nloc;
        if self.cell_labels.len() != n_cells {
            return Err(MeshIntegrityError::LengthMismatch("cell_labels length does not match cell count"));
        }
        if !self.facets.len().is_multiple_of(self.dim) {
            return Err(MeshIntegrityError::LengthMismatch("facets length is not a multiple of dim"));
        }
        let n_facets = self.facets.len() / self.dim;
        if self.facet_labels.len() != n_facets {
            return Err(MeshIntegrityError::LengthMismatch("facet_labels length does not match facet count"));
        }
        let npair = self.dim * nloc / 2;
        for (p_idx, pt) in self.points.chunks_exact(self.dim).enumerate() {
            for (axis, &coord) in pt.iter().enumerate() {
                if !coord.is_finite() {
                    return Err(MeshIntegrityError::NonFiniteCoordinate { point: p_idx, axis });
                }
            }
        }
        let mut facet_cells: HashMap<[u32; 3], Vec<(usize, i32)>> =
            HashMap::with_capacity(n_cells.saturating_mul(nloc) / 2);
        let facet_key = |vs: &[u32]| -> [u32; 3] {
            if vs.len() == 2 {
                [vs[0].min(vs[1]), vs[0].max(vs[1]), u32::MAX]
            } else {
                let mut k = [vs[0], vs[1], vs[2]];
                k.sort_unstable();
                k
            }
        };
        let key_verts = |k: [u32; 3]| -> Vec<u32> {
            if self.dim == 2 { vec![k[0], k[1]] } else { k.to_vec() }
        };
        for c in 0..n_cells {
            let vs = &self.cells[c * nloc..(c + 1) * nloc];
            for (i, &v) in vs.iter().enumerate() {
                if (v as usize) >= n_verts {
                    return Err(MeshIntegrityError::VertexOutOfBounds {
                        context: "cell",
                        index: c,
                        vertex: v,
                        n_points: n_verts,
                    });
                }
                if vs[..i].contains(&v) {
                    return Err(MeshIntegrityError::DuplicateVertexInSimplex {
                        context: "cell",
                        index: c,
                        vertex: v,
                    });
                }
            }
            let det = if self.dim == 2 {
                let p = |v: u32| -> [f64; 2] {
                    let i = v as usize * 2;
                    [self.points[i], self.points[i + 1]]
                };
                let (a, b, cc) = (p(vs[0]), p(vs[1]), p(vs[2]));
                (b[0] - a[0]) * (cc[1] - a[1]) - (cc[0] - a[0]) * (b[1] - a[1])
            } else {
                let p = |v: u32| -> [f64; 3] {
                    let i = v as usize * 3;
                    [self.points[i], self.points[i + 1], self.points[i + 2]]
                };
                let (p0, p1, p2, p3) = (p(vs[0]), p(vs[1]), p(vs[2]), p(vs[3]));
                let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                let w = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
                u[0] * (v[1] * w[2] - v[2] * w[1])
                    - u[1] * (v[0] * w[2] - v[2] * w[0])
                    + u[2] * (v[0] * w[1] - v[1] * w[0])
            };
            if !(det.is_finite() && det > 0.0) {
                return Err(MeshIntegrityError::NonPositiveMeasure { cell: c, det });
            }
            let lab = self.cell_labels[c];
            let mut fbuf = [0u32; 3];
            for skip in 0..nloc {
                let mut idx = 0;
                for (i, &x) in vs.iter().enumerate() {
                    if i != skip {
                        fbuf[idx] = x;
                        idx += 1;
                    }
                }
                let key = facet_key(&fbuf[..self.dim]);
                let entry = facet_cells.entry(key).or_default();
                entry.push((c, lab));
                if entry.len() > 2 {
                    return Err(MeshIntegrityError::NonManifoldFacet {
                        vertices: key_verts(key),
                        incident_cells: entry.len(),
                    });
                }
            }
        }
        // Track duplicate geometric facets and ensure interior sheet labels never clash with exterior/interface roles.
        let mut seen_facets: HashMap<[u32; 3], i32> = HashMap::with_capacity(n_facets);
        let mut tag_roles: HashMap<i32, (bool, bool)> = HashMap::new();
        for f in 0..n_facets {
            let fv = &self.facets[f * self.dim..(f + 1) * self.dim];
            for (i, &v) in fv.iter().enumerate() {
                if (v as usize) >= n_verts {
                    return Err(MeshIntegrityError::VertexOutOfBounds {
                        context: "facet",
                        index: f,
                        vertex: v,
                        n_points: n_verts,
                    });
                }
                if fv[..i].contains(&v) {
                    return Err(MeshIntegrityError::DuplicateVertexInSimplex {
                        context: "facet",
                        index: f,
                        vertex: v,
                    });
                }
            }
            let key = facet_key(fv);
            let lab = self.facet_labels[f];
            let Some(cs) = facet_cells.get(&key) else {
                return Err(MeshIntegrityError::FacetNotInAnyCell {
                    facet: f,
                    vertices: key_verts(key),
                });
            };
            if let Some(&prev_lab) = seen_facets.get(&key) {
                return Err(MeshIntegrityError::InconsistentFacetTag {
                    vertices: key_verts(key),
                    detail: format!("duplicate facet entry (labels {prev_lab} and {lab})"),
                });
            }
            seen_facets.insert(key, lab);
            let role = tag_roles.entry(lab).or_insert((false, false));
            if cs.len() == 2 && cs[0].1 == cs[1].1 {
                role.0 = true;
            } else {
                role.1 = true;
            }
            if role.0 && role.1 {
                return Err(MeshIntegrityError::InconsistentFacetTag {
                    vertices: key_verts(key),
                    detail: format!("label {lab} mixes interior same-domain sheet and non-sheet facets"),
                });
            }
        }
        if self.cell_edges.len() != n_cells * npair {
            return Err(MeshIntegrityError::LengthMismatch("cell_edges length does not match cells"));
        }
        for (i, &(u, v)) in self.edges.iter().enumerate() {
            if u >= v || (v as usize) >= n_verts {
                return Err(MeshIntegrityError::InvalidCurvedEdge {
                    edge_index: i as u32,
                    detail: format!("mesh edge ({u}, {v}) invalid for {n_verts} points"),
                });
            }
        }
        for &eid in &self.cell_edges {
            if (eid as usize) >= self.edges.len() {
                return Err(MeshIntegrityError::InvalidCurvedEdge {
                    edge_index: eid,
                    detail: format!("cell_edges reference {eid} out of range (n_edges = {})", self.edges.len()),
                });
            }
        }
        let mut prev_edge: Option<u32> = None;
        for &(eid, pos) in &self.edge_nodes {
            if let Some(prev) = prev_edge {
                if eid <= prev {
                    return Err(MeshIntegrityError::InvalidCurvedEdge {
                        edge_index: eid,
                        detail: "edge_nodes are not strictly sorted by edge index".into(),
                    });
                }
            }
            prev_edge = Some(eid);
            let Some(&(u, v)) = self.edges.get(eid as usize) else {
                return Err(MeshIntegrityError::InvalidCurvedEdge {
                    edge_index: eid,
                    detail: format!("edge index {eid} out of bounds (n_edges = {})", self.edges.len()),
                });
            };
            if u >= v || (u as usize) >= n_verts || (v as usize) >= n_verts {
                return Err(MeshIntegrityError::InvalidCurvedEdge {
                    edge_index: eid,
                    detail: format!("invalid endpoints ({u}, {v}) for {n_verts} points"),
                });
            }
            if !pos.iter().all(|x| x.is_finite()) {
                return Err(MeshIntegrityError::InvalidCurvedEdge {
                    edge_index: eid,
                    detail: "curved edge midpoint has non-finite coordinate".into(),
                });
            }
        }
        let mut prev_qp: Option<u32> = None;
        for &eid in &self.qp_edges {
            if let Some(prev) = prev_qp {
                if eid <= prev {
                    return Err(MeshIntegrityError::InvalidQpEdge {
                        edge_index: eid,
                        detail: "qp_edges are not strictly sorted by edge index".into(),
                    });
                }
            }
            prev_qp = Some(eid);
            let Some(&(u, v)) = self.edges.get(eid as usize) else {
                return Err(MeshIntegrityError::InvalidQpEdge {
                    edge_index: eid,
                    detail: format!("edge index {eid} out of bounds (n_edges = {})", self.edges.len()),
                });
            };
            if u >= v || (u as usize) >= n_verts || (v as usize) >= n_verts {
                return Err(MeshIntegrityError::InvalidQpEdge {
                    edge_index: eid,
                    detail: format!("invalid endpoints ({u}, {v}) for {n_verts} points"),
                });
            }
            if self.edge_nodes.binary_search_by_key(&eid, |x| x.0).is_err() {
                return Err(MeshIntegrityError::InvalidQpEdge {
                    edge_index: eid,
                    detail: "quarter-point edge has no matching edge_nodes entry".into(),
                });
            }
        }
        Ok(())
    }

    /// Sets the curved edges of a mesh: `(vertex pair, node position)`; pairs that are not an edge of the mesh are ignored (the
    /// connectivity is what decides), and the list is kept sorted by edge index.
    pub fn set_edge_nodes(&mut self, nodes: &[([u32; 2], [f64; 3])]) {
        if self.dim != 2 && self.dim != 3 {
            self.edge_nodes.clear();
            self.qp_edges.clear();
            return;
        }
        let mut v: Vec<(u32, [f64; 3])> = nodes
            .iter()
            .filter_map(|&([a, b], p)| self.edges.binary_search(&(a.min(b), a.max(b))).ok().map(|e| (e as u32, p)))
            .collect();
        v.sort_by_key(|e| e.0);
        v.dedup_by_key(|e| e.0);
        self.edge_nodes = v;
        self.qp_edges.retain(|e| self.edge_nodes.binary_search_by_key(e, |x| x.0).is_ok());
    }

    /// Marks which curved edges carry Barsoum quarter-point singularity shifts (`AF-2`), sorted by edge index.
    pub fn set_qp_edges(&mut self, edges: &[[u32; 2]]) {
        if self.dim != 2 && self.dim != 3 {
            self.qp_edges.clear();
            return;
        }
        let mut v: Vec<u32> = edges
            .iter()
            .filter_map(|&[a, b]| {
                let e = self.edges.binary_search(&(a.min(b), a.max(b))).ok()? as u32;
                self.edge_nodes.binary_search_by_key(&e, |x| x.0).ok().map(|_| e)
            })
            .collect();
        v.sort_unstable();
        v.dedup();
        self.qp_edges = v;
    }

    /// Whether edge `e` (index into `edges()`) carries a quarter-point singularity node.
    pub fn edge_is_qp(&self, e: usize) -> bool {
        self.qp_edges.binary_search(&(e as u32)).is_ok()
    }

    /// Whether the mesh carries any quarter-point shifted edge.
    pub fn has_qp(&self) -> bool {
        !self.qp_edges.is_empty()
    }

    /// Whether at least one edge of cell `c` is a quarter-point shifted edge (`AF-2`).
    pub fn cell_is_qp(&self, c: usize) -> bool {
        !self.qp_edges.is_empty() && self.dim == 3 && self.cell_edge_ids(c).iter().any(|&e| self.edge_is_qp(e as usize))
    }

    /// The quarter-point edges as sorted vertex pairs `[a, b]`, sorted by pair (`CellMesh::qp_edges`).
    pub fn qp_edge_list(&self) -> Vec<[u32; 2]> {
        let mut v: Vec<[u32; 2]> = self.qp_edges.iter().map(|&e| [self.edges[e as usize].0, self.edges[e as usize].1]).collect();
        v.sort_unstable();
        v
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
    /// Facets carrying the given 1-based boundary label from `LabelMap`.
    ///
    /// In a 3D sheet-stack mesh (`meshing::sheet3d`), each interior 2D conductive sheet triangle (`Op::Sheet`)
    /// appears **exactly once** in `facets` / `facet_labels` with its 1-based boundary label, shared conformally
    /// by the lower and upper tetrahedra.
    pub fn sheet_facets(&self, label: i32) -> Vec<&[u32]> {
        (0..self.n_facets()).filter(|&f| self.facet_labels[f] == label).map(|f| self.facet(f)).collect()
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
        let n_pts = nm.n_points();
        for c in 0..nm.n_cells() {
            let vs = nm.cell(c);
            if vs.len() != nm.types[c].nv() {
                return Err(format!("cell {c}: wrong node count {}", vs.len()));
            }
            if let Some(&v) = vs.iter().find(|&&v| (v as usize) >= n_pts) {
                return Err(format!("cell {c}: vertex {v} out of range (n_points = {n_pts})"));
            }
        }
        for f in 0..nm.n_facets() {
            let fv = nm.facet(f);
            if let Some(&v) = fv.iter().find(|&&v| (v as usize) >= n_pts) {
                return Err(format!("facet {f}: vertex {v} out of range (n_points = {n_pts})"));
            }
        }
        let mut report = SplitReport::default();
        let mut points: Vec<f64> = nm.points.iter().flat_map(|p| p.iter().copied()).collect();
        let (mut cells, mut cell_labels) = (vec![], vec![]);
        let (mut facets, mut facet_labels) = (vec![], vec![]);
        let orient_tet = |pts: &[f64], mut tet: [u32; 4]| -> [u32; 4] {
            let p = |v: u32| [pts[v as usize * 3], pts[v as usize * 3 + 1], pts[v as usize * 3 + 2]];
            let (p0, p1, p2, p3) = (p(tet[0]), p(tet[1]), p(tet[2]), p(tet[3]));
            let (u, v, w) = (
                [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]],
                [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]],
                [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]],
            );
            let det = u[0] * (v[1] * w[2] - v[2] * w[1])
                - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]);
            if det < 0.0 {
                tet.swap(2, 3);
            }
            tet
        };
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
                        let vs = nm.cell(c);
                        for mut tr in quad_tris([vs[0], vs[1], vs[2], vs[3]]) {
                            let (ax, ay) = (points[tr[0] as usize * 2], points[tr[0] as usize * 2 + 1]);
                            let (bx, by) = (points[tr[1] as usize * 2], points[tr[1] as usize * 2 + 1]);
                            let (cx, cy) = (points[tr[2] as usize * 2], points[tr[2] as usize * 2 + 1]);
                            if (bx - ax) * (cy - ay) - (cx - ax) * (by - ay) < 0.0 {
                                tr.swap(1, 2);
                            }
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
                                cells.extend_from_slice(&orient_tet(&points, [r[v0], tr[0], tr[1], tr[2]]));
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
                            cells.extend_from_slice(&orient_tet(&points, [cid, tr[0], tr[1], tr[2]]));
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
        // Retain all curved edges and quarter-point edges whose endpoints remain an edge in `mesh` (even when outer PML cells were split).
        if !nm.edge_nodes.is_empty() {
            let nodes: Vec<([u32; 2], [f64; 3])> = nm
                .edge_nodes
                .iter()
                .map(|&(pr, p)| {
                    let mut p3 = [0.0; 3];
                    p3[..d].copy_from_slice(&p);
                    (pr, p3)
                })
                .collect();
            mesh.set_edge_nodes(&nodes);
            if !nm.qp_edges.is_empty() {
                mesh.set_qp_edges(&nm.qp_edges);
            }
        }
        Ok((mesh, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_2d() -> SimplexMesh {
        SimplexMesh::new(
            2,
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![0, 1, 2, 1, 3, 2],
            vec![1, 1],
            vec![0, 1, 1, 3, 3, 2, 2, 0],
            vec![1, 2, 3, 4],
        )
    }

    fn valid_3d() -> SimplexMesh {
        SimplexMesh::new(
            3,
            vec![
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
                1.0, 1.0, 1.0,
            ],
            vec![0, 1, 2, 3, 1, 2, 3, 4],
            vec![1, 1],
            vec![0, 1, 2, 0, 1, 3],
            vec![1, 2],
        )
    }

    #[test]
    fn simplex_mesh_check_accepts_valid_and_rejects_defects() {
        valid_2d().check().unwrap();
        valid_3d().check().unwrap();

        let mut bad_coord = valid_2d();
        bad_coord.points[3] = f64::NAN;
        assert!(matches!(bad_coord.check(), Err(MeshIntegrityError::NonFiniteCoordinate { .. })));

        let oob = SimplexMesh::new(2, vec![0.0; 6], vec![0, 1, 9], vec![1], vec![], vec![]);
        assert!(matches!(oob.check(), Err(MeshIntegrityError::VertexOutOfBounds { .. })));

        let dup_v = SimplexMesh::new(2, vec![0.0; 6], vec![0, 1, 1], vec![1], vec![], vec![]);
        assert!(matches!(dup_v.check(), Err(MeshIntegrityError::DuplicateVertexInSimplex { .. })));

        let mut inverted_2d = valid_2d();
        inverted_2d.cells[1] = 2;
        inverted_2d.cells[2] = 1;
        assert!(matches!(inverted_2d.check(), Err(MeshIntegrityError::NonPositiveMeasure { cell: 0, .. })));

        let mut inverted_3d = valid_3d();
        inverted_3d.cells.swap(2, 3);
        assert!(matches!(inverted_3d.check(), Err(MeshIntegrityError::NonPositiveMeasure { cell: 0, .. })));

        let non_manifold = SimplexMesh::new(
            2,
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.5, 0.5],
            vec![0, 1, 2, 0, 3, 1, 0, 1, 4],
            vec![1, 1, 1],
            vec![],
            vec![],
        );
        assert!(matches!(non_manifold.check(), Err(MeshIntegrityError::NonManifoldFacet { .. })));

        let mut orphan_facet = valid_2d();
        orphan_facet.facets.extend_from_slice(&[0, 3]);
        orphan_facet.facet_labels.push(9);
        assert!(matches!(orphan_facet.check(), Err(MeshIntegrityError::FacetNotInAnyCell { .. })));

        // Sheet consistency: same label cannot be both an interior same-domain sheet and an exterior wall.
        let mut mixed_sheet = valid_2d();
        mixed_sheet.facets.extend_from_slice(&[1, 2]);
        mixed_sheet.facet_labels.push(1);
        assert!(matches!(mixed_sheet.check(), Err(MeshIntegrityError::InconsistentFacetTag { .. })));

        // Non-finite curved edge midpoint must fail check().
        let mut bad_curve = valid_3d();
        bad_curve.set_edge_nodes(&[([0, 1], [0.5, f64::INFINITY, 0.0])]);
        assert!(matches!(bad_curve.check(), Err(MeshIntegrityError::InvalidCurvedEdge { .. })));
    }
}

