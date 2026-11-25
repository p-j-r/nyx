//! Mixed-cell meshes: triangles and quads in 2D; tetrahedra, prisms and hexahedra in 3D.
//!
//! Stand-in for nyx's `CellMesh`. Node orderings (all positively oriented):
//!
//! * `Tri`, `Quad`: counter-clockwise.
//! * `Tet`: `(v1 - v0, v2 - v0, v3 - v0)` right-handed.
//! * `Hex`: nodes 0..3 are the bottom quad counter-clockwise seen from the top, nodes 4..7 the
//!   top quad above them.
//! * `Prism`: nodes 0..2 are the base triangle, nodes 3..5 the top triangle above them, with
//!   `(0,1,2,3)` right-handed (the opposite handedness of the VTK wedge).
//!
//! [`CellMesh::split_to_simplices`] is the fallback that lets today's simplex-only solvers run on
//! a mixed mesh. It is conforming because the diagonal of every quad face is a function of the
//! face's global vertex ids alone (the diagonal through its smallest id), and every hex or prism is
//! split as a cone from its smallest-id vertex.

use crate::adapt::work::{FacetKey, SimplexMesh, facet_key, simplex_signed_measure};
use std::collections::HashMap;

/// Cell shapes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CellType {
    Tri,
    Quad,
    Tet,
    Prism,
    Hex,
}

const TRI_FACES: &[&[u8]] = &[&[0, 1], &[1, 2], &[2, 0]];
const QUAD_FACES: &[&[u8]] = &[&[0, 1], &[1, 2], &[2, 3], &[3, 0]];
const TET_FACES: &[&[u8]] = &[&[1, 2, 3], &[0, 3, 2], &[0, 1, 3], &[0, 2, 1]];
const PRISM_FACES: &[&[u8]] = &[&[0, 2, 1], &[3, 4, 5], &[0, 1, 4, 3], &[1, 2, 5, 4], &[2, 0, 3, 5]];
const HEX_FACES: &[&[u8]] = &[&[0, 3, 2, 1], &[4, 5, 6, 7], &[0, 1, 5, 4], &[1, 2, 6, 5], &[2, 3, 7, 6], &[3, 0, 4, 7]];

impl CellType {
    /// Number of nodes.
    pub fn nv(self) -> usize {
        match self {
            CellType::Tri => 3,
            CellType::Quad | CellType::Tet => 4,
            CellType::Prism => 6,
            CellType::Hex => 8,
        }
    }

    /// Spatial dimension of the shape.
    pub fn dim(self) -> usize {
        match self {
            CellType::Tri | CellType::Quad => 2,
            _ => 3,
        }
    }

    /// Local node lists of the faces (edges for 2D shapes).
    pub fn faces(self) -> &'static [&'static [u8]] {
        match self {
            CellType::Tri => TRI_FACES,
            CellType::Quad => QUAD_FACES,
            CellType::Tet => TET_FACES,
            CellType::Prism => PRISM_FACES,
            CellType::Hex => HEX_FACES,
        }
    }

    /// Code used by the binary format.
    pub fn code(self) -> u8 {
        self as u8
    }

    /// Inverse of [`CellType::code`].
    pub fn from_code(c: u8) -> Option<CellType> {
        Some(match c {
            0 => CellType::Tri,
            1 => CellType::Quad,
            2 => CellType::Tet,
            3 => CellType::Prism,
            4 => CellType::Hex,
            _ => return None,
        })
    }
}

/// Mixed-cell mesh in `D` dimensions.
#[derive(Clone, Debug, Default)]
pub struct CellMesh<const D: usize> {
    /// Vertex coordinates.
    pub points: Vec<[f64; D]>,
    /// Shape of every cell.
    pub types: Vec<CellType>,
    /// `offsets[c]..offsets[c + 1]` is the slice of `conn` for cell `c`.
    pub offsets: Vec<u32>,
    /// Flat connectivity.
    pub conn: Vec<u32>,
    /// Domain label of every cell.
    pub label: Vec<i32>,
    /// Boundary-facet offsets into `facet_conn`.
    pub facet_offsets: Vec<u32>,
    /// Flat boundary-facet connectivity (2 nodes in 2D; 3 or 4 in 3D).
    pub facet_conn: Vec<u32>,
    /// Label of every boundary facet.
    pub facet_label: Vec<i32>,
    /// Curved edges (ADR 82): for an edge, given as its sorted vertex pair, the position of its mid-edge geometry node, so that the
    /// cells around it are quadratic (isoparametric) there. Sorted by pair. Empty by default; an edge that is not listed is straight
    /// (its node is the midpoint), so an empty list is the straight mesh of ADR 18.
    pub edge_nodes: Vec<([u32; 2], [f64; D])>,
    /// Optional curved quad-face center nodes (for 9-node `Quad`, 18-node `Prism`, 27-node `Hex`), keyed by [`canonical_quad_face`]
    /// and sorted by key. Empty by default; unlisted quad faces use transfinite / bilinear face centers.
    pub quad_face_nodes: Vec<([u32; 4], [f64; D])>,
    /// Barsoum quarter-point singularity edges (`AF-2`), given as sorted vertex pairs `[a, b]` (`a < b`), sorted by pair.
    pub qp_edges: Vec<[u32; 2]>,
}

/// Canonical cyclic key of a quadrilateral `[q0, q1, q2, q3]`: starts at the smallest vertex id
/// and steps in the cyclic direction of the smaller adjacent vertex id. All 8 dihedral (`D_4`)
/// traversal orders of the same untwisted quad map to the identical `[u32; 4]` key.
pub fn canonical_quad_face(q: [u32; 4]) -> [u32; 4] {
    let m = (0..4).min_by_key(|&i| q[i]).unwrap();
    let next = q[(m + 1) % 4];
    let prev = q[(m + 3) % 4];
    if next <= prev {
        [q[m], next, q[(m + 2) % 4], prev]
    } else {
        [q[m], prev, q[(m + 2) % 4], next]
    }
}

impl<const D: usize> CellMesh<D> {
    /// Empty mesh with the offset arrays primed.
    pub fn new() -> Self {
        CellMesh { offsets: vec![0], facet_offsets: vec![0], ..Default::default() }
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.types.len()
    }

    /// Number of vertices.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// Number of labelled boundary facets.
    pub fn n_facets(&self) -> usize {
        self.facet_label.len()
    }

    /// Node ids of cell `c`.
    pub fn cell(&self, c: usize) -> &[u32] {
        &self.conn[self.offsets[c] as usize..self.offsets[c + 1] as usize]
    }

    /// Node ids of boundary facet `f`.
    pub fn facet(&self, f: usize) -> &[u32] {
        &self.facet_conn[self.facet_offsets[f] as usize..self.facet_offsets[f + 1] as usize]
    }

    /// Facets carrying the given 1-based boundary label from `LabelMap`.
    ///
    /// In a 3D sheet-stack mesh (`meshing::sheet3d`), each interior 2D conductive sheet triangle (`Op::Sheet`)
    /// appears **exactly once** in `facet_conn` / `facet_label` (and `SimplexMesh::facets` / `facet_labels`)
    /// with its 1-based boundary label, shared conformally by the lower and upper tetrahedra.
    pub fn sheet_facets(&self, label: i32) -> Vec<&[u32]> {
        (0..self.n_facets()).filter(|&f| self.facet_label[f] == label).map(|f| self.facet(f)).collect()
    }

    /// Append a cell.
    pub fn push_cell(&mut self, t: CellType, vs: &[u32], label: i32) {
        assert_eq!(t.dim(), D, "cell dimension does not match mesh dimension");
        assert_eq!(vs.len(), t.nv());
        if self.offsets.is_empty() {
            self.offsets.push(0);
        }
        self.types.push(t);
        self.conn.extend_from_slice(vs);
        self.offsets.push(self.conn.len() as u32);
        self.label.push(label);
    }

    /// Append a boundary facet.
    pub fn push_facet(&mut self, vs: &[u32], label: i32) {
        if self.facet_offsets.is_empty() {
            self.facet_offsets.push(0);
        }
        self.facet_conn.extend_from_slice(vs);
        self.facet_offsets.push(self.facet_conn.len() as u32);
        self.facet_label.push(label);
    }

    /// Sets the curved edges: pairs are put in ascending vertex order, the list is sorted by pair, and when an edge is given twice
    /// the last position wins.
    pub fn set_edge_nodes(&mut self, nodes: impl IntoIterator<Item = ([u32; 2], [f64; D])>) {
        let mut v: Vec<([u32; 2], [f64; D])> = nodes.into_iter().map(|([a, b], p)| ([a.min(b), a.max(b)], p)).collect();
        v.sort_by_key(|e| e.0);
        let mut out: Vec<([u32; 2], [f64; D])> = Vec::with_capacity(v.len());
        for e in v {
            match out.last_mut() {
                Some(l) if l.0 == e.0 => *l = e,
                _ => out.push(e),
            }
        }
        self.edge_nodes = out;
    }

    /// Position of the mid-edge node of the edge `a`-`b`, `None` when the edge is straight.
    pub fn edge_node(&self, a: u32, b: u32) -> Option<[f64; D]> {
        let key = [a.min(b), a.max(b)];
        self.edge_nodes.binary_search_by_key(&key, |e| e.0).ok().map(|i| self.edge_nodes[i].1)
    }

    /// Sets the curved quad-face center nodes: cyclic quadruples are canonicalised with [`canonical_quad_face`], sorted, and when a face
    /// is given twice the last position wins.
    pub fn set_quad_face_nodes(&mut self, nodes: impl IntoIterator<Item = ([u32; 4], [f64; D])>) {
        let mut v: Vec<([u32; 4], [f64; D])> = nodes.into_iter().map(|(q, p)| (canonical_quad_face(q), p)).collect();
        v.sort_by_key(|e| e.0);
        let mut out: Vec<([u32; 4], [f64; D])> = Vec::with_capacity(v.len());
        for e in v {
            match out.last_mut() {
                Some(l) if l.0 == e.0 => *l = e,
                _ => out.push(e),
            }
        }
        self.quad_face_nodes = out;
    }

    /// Position of the curved center node of the quad `q` (given in cyclic order), `None` when not listed.
    pub fn quad_face_node(&self, q: [u32; 4]) -> Option<[f64; D]> {
        let key = canonical_quad_face(q);
        self.quad_face_nodes.binary_search_by_key(&key, |e| e.0).ok().map(|i| self.quad_face_nodes[i].1)
    }

    /// Drops the curved edges and quad-face centers that do not belong to any cell (after the connectivity changed). Never adds one.
    pub fn retain_existing_edge_nodes(&mut self) {
        if self.edge_nodes.is_empty() && self.qp_edges.is_empty() && self.quad_face_nodes.is_empty() {
            return;
        }
        let mut exist: std::collections::HashSet<[u32; 2]> = std::collections::HashSet::new();
        let mut exist_quads: std::collections::HashSet<[u32; 4]> = std::collections::HashSet::new();
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
            let pairs: &[(usize, usize)] = match self.types[c] {
                CellType::Tri => &[(0, 1), (0, 2), (1, 2)],
                CellType::Quad => &[(0, 1), (1, 2), (2, 3), (3, 0)],
                CellType::Tet => &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
                CellType::Prism => &[(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)],
                CellType::Hex => &[(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
            };
            for &(i, j) in pairs {
                exist.insert([vs[i].min(vs[j]), vs[i].max(vs[j])]);
            }
            if !self.quad_face_nodes.is_empty() {
                if self.types[c] == CellType::Quad {
                    exist_quads.insert(canonical_quad_face([vs[0], vs[1], vs[2], vs[3]]));
                } else {
                    for f in self.types[c].faces() {
                        if f.len() == 4 {
                            exist_quads.insert(canonical_quad_face([
                                vs[f[0] as usize],
                                vs[f[1] as usize],
                                vs[f[2] as usize],
                                vs[f[3] as usize],
                            ]));
                        }
                    }
                }
            }
        }
        self.edge_nodes.retain(|e| exist.contains(&e.0));
        self.qp_edges.retain(|e| exist.contains(e));
        self.quad_face_nodes.retain(|q| exist_quads.contains(&q.0));
    }

    /// Number of cells of each shape.
    pub fn counts(&self) -> HashMap<CellType, usize> {
        let mut m = HashMap::new();
        for &t in &self.types {
            *m.entry(t).or_insert(0) += 1;
        }
        m
    }

    /// Keys and vertex counts of the faces of cell `c`.
    pub fn cell_face_keys(&self, c: usize) -> Vec<FacetKey> {
        let vs = self.cell(c);
        let mut out = Vec::with_capacity(6);
        for f in self.types[c].faces() {
            let mut tmp = [0u32; 4];
            for (i, &l) in f.iter().enumerate() {
                tmp[i] = vs[l as usize];
            }
            out.push(facet_key(&tmp[..f.len()]));
        }
        out
    }

    /// Face key to the cells that contain it.
    pub fn face_cells(&self) -> HashMap<FacetKey, Vec<u32>> {
        let mut m: HashMap<FacetKey, Vec<u32>> = HashMap::new();
        for c in 0..self.n_cells() {
            for k in self.cell_face_keys(c) {
                m.entry(k).or_default().push(c as u32);
            }
        }
        m
    }

    /// Faces that belong to exactly one cell, sorted.
    pub fn boundary_keys(&self) -> Vec<FacetKey> {
        let mut v: Vec<FacetKey> =
            self.face_cells().into_iter().filter(|(_, cs)| cs.len() == 1).map(|(k, _)| k).collect();
        v.sort_unstable();
        v
    }

    /// Boundary faces (belonging to exactly one cell) with their nodes in the cyclic order of the
    /// owning cell's face table, in cell order (deterministic).
    pub fn boundary_faces_ordered(&self) -> Vec<Vec<u32>> {
        let fc = self.face_cells();
        let mut out = Vec::new();
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
            for f in self.types[c].faces() {
                let nodes: Vec<u32> = f.iter().map(|&l| vs[l as usize]).collect();
                if fc[&facet_key(&nodes)].len() == 1 {
                    out.push(nodes);
                }
            }
        }
        out
    }

    /// Positive measure of cell `c` (area or volume); negative if the cell is inverted.
    pub fn measure(&self, c: usize) -> f64 {
        let vs = self.cell(c);
        let p = |i: usize| self.points[vs[i] as usize];
        match self.types[c] {
            CellType::Tri => simplex_signed_measure::<D>(&[p(0), p(1), p(2)]),
            CellType::Quad => {
                simplex_signed_measure::<D>(&[p(0), p(1), p(2)]) + simplex_signed_measure::<D>(&[p(0), p(2), p(3)])
            }
            CellType::Tet => simplex_signed_measure::<D>(&[p(0), p(1), p(2), p(3)]),
            CellType::Prism | CellType::Hex => {
                let mut tets = Vec::new();
                split_cell3(self.types[c], vs, &self.points, &mut tets);
                let pt = |v: u32| self.points[v as usize];
                tets.iter().map(|t| simplex_signed_measure::<D>(&[pt(t[0]), pt(t[1]), pt(t[2]), pt(t[3])])).sum()
            }
        }
    }

    /// Sum of all cell measures.
    pub fn total_measure(&self) -> f64 {
        (0..self.n_cells()).map(|c| self.measure(c)).sum()
    }

    /// Validity check: positive measure of every cell (and of every sub-simplex of a split cell),
    /// every face shared by at most two cells, labelled facets lie on the boundary.
    pub fn check(&self) -> Result<(), String> {
        let np = self.n_points() as u32;
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
            if vs.len() != self.types[c].nv() {
                return Err(format!("cell {c}: wrong node count"));
            }
            for (i, &v) in vs.iter().enumerate() {
                if v >= np {
                    return Err(format!("cell {c}: vertex {v} out of range"));
                }
                if vs[..i].contains(&v) {
                    return Err(format!("cell {c}: repeated vertex {v}"));
                }
            }
            let m = self.measure(c);
            if !(m > 0.0) {
                return Err(format!("cell {c} ({:?}): non-positive measure {m:e}", self.types[c]));
            }
            // the measure of a split prism / hex ignores inversion, the corner Jacobians do not
            if !matches!(self.types[c], CellType::Tri | CellType::Tet) {
                let sj = crate::mesh::quality::cell_scaled_jacobian_min(self, c);
                if !(sj > 0.0) {
                    return Err(format!("cell {c} ({:?}): inverted (min scaled Jacobian {sj:e})", self.types[c]));
                }
            }
        }
        let fc = self.face_cells();
        for (k, cs) in &fc {
            if cs.len() > 2 {
                return Err(format!("face {k:?} shared by {} cells", cs.len()));
            }
        }
        for f in 0..self.n_facets() {
            match fc.get(&facet_key(self.facet(f))) {
                Some(cs) if cs.len() == 1 => {}
                Some(_) => return Err(format!("labelled facet {f} is interior")),
                None => return Err(format!("labelled facet {f} not in any cell")),
            }
        }
        Ok(())
    }

    /// Convert a simplex mesh without changing anything.
    pub fn from_simplex(m: &SimplexMesh<D>) -> Self {
        let t = if D == 2 { CellType::Tri } else { CellType::Tet };
        let mut out = CellMesh::<D>::new();
        out.points = m.points.clone();
        for c in 0..m.n_cells() {
            out.push_cell(t, m.cell(c), m.cell_label[c]);
        }
        for f in 0..m.n_facets() {
            out.push_facet(m.facet(f), m.facet_label[f]);
        }
        out
    }

    /// Whether every cell is a triangle (2D) or tetrahedron (3D) and every facet has `D` vertices.
    pub fn is_pure_simplex(&self) -> bool {
        let t = if D == 2 { CellType::Tri } else { CellType::Tet };
        self.types.iter().all(|&x| x == t) && (0..self.n_facets()).all(|f| self.facet(f).len() == D)
    }

    /// `Some` simplex mesh when every cell is a triangle / tetrahedron and every facet an
    /// edge / triangle, else `None` (use [`CellMesh::split_to_simplices`] to convert).
    pub fn as_simplex(&self) -> Option<SimplexMesh<D>> {
        if !self.is_pure_simplex() {
            return None;
        }
        Some(SimplexMesh {
            points: self.points.clone(),
            cells: self.conn.clone(),
            cell_label: self.label.clone(),
            facets: self.facet_conn.clone(),
            facet_label: self.facet_label.clone(),
        })
    }

    /// Split into a pure simplex mesh (triangles / tetrahedra), preserving cell and facet labels.
    pub fn split_to_simplices(&self) -> Result<SimplexMesh<D>, String> {
        let mut out = SimplexMesh::<D> { points: self.points.clone(), ..Default::default() };
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
            let lab = self.label[c];
            match self.types[c] {
                CellType::Tri | CellType::Tet => out.push_cell(vs, lab),
                CellType::Quad => {
                    let [t0, t1] = split_quad(&[vs[0], vs[1], vs[2], vs[3]]);
                    out.push_cell(&t0, lab);
                    out.push_cell(&t1, lab);
                }
                CellType::Prism | CellType::Hex => {
                    let mut tets = Vec::new();
                    split_cell3(self.types[c], vs, &self.points, &mut tets);
                    for t in tets {
                        out.push_cell(&t, lab);
                    }
                }
            }
        }
        for f in 0..self.n_facets() {
            let vs = self.facet(f);
            let lab = self.facet_label[f];
            match (D, vs.len()) {
                (2, 2) | (3, 3) => out.push_facet(vs, lab),
                (3, 4) => {
                    let [t0, t1] = split_quad(&[vs[0], vs[1], vs[2], vs[3]]);
                    out.push_facet(&t0, lab);
                    out.push_facet(&t1, lab);
                }
                _ => return Err(format!("facet {f}: unsupported node count {}", vs.len())),
            }
        }
        Ok(out)
    }
}

/// Split a quad `[a, b, c, d]` (cyclic order) into two triangles along the diagonal through its
/// smallest vertex id. The result depends only on the vertex ids, so two cells sharing the quad
/// agree on the diagonal.
pub fn split_quad(q: &[u32; 4]) -> [[u32; 3]; 2] {
    let m = (0..4).min_by_key(|&i| q[i]).unwrap();
    let at = |i: usize| q[(m + i) % 4];
    [[at(0), at(1), at(2)], [at(0), at(2), at(3)]]
}

/// Split a 3D cell (prism or hex) into tetrahedra as a cone from its smallest-id vertex over the
/// faces that do not contain it. Tets are returned positively oriented when the point geometry
/// allows it.
fn split_cell3<const D: usize>(t: CellType, vs: &[u32], pts: &[[f64; D]], out: &mut Vec<[u32; 4]>) {
    let apex_local = (0..vs.len()).min_by_key(|&i| vs[i]).unwrap() as u8;
    let apex = vs[apex_local as usize];
    for f in t.faces() {
        if f.contains(&apex_local) {
            continue;
        }
        let g: Vec<u32> = f.iter().map(|&l| vs[l as usize]).collect();
        let tris: Vec<[u32; 3]> =
            if g.len() == 3 { vec![[g[0], g[1], g[2]]] } else { split_quad(&[g[0], g[1], g[2], g[3]]).to_vec() };
        for tr in tris {
            let mut tet = [apex, tr[0], tr[1], tr[2]];
            let p = |v: u32| pts[v as usize];
            if simplex_signed_measure::<D>(&[p(tet[0]), p(tet[1]), p(tet[2]), p(tet[3])]) < 0.0 {
                tet.swap(2, 3);
            }
            out.push(tet);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_hex_points() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    }

    #[test]
    fn hex_splits_into_six_tets_with_full_volume_for_any_vertex_numbering() {
        // permute global ids so every vertex takes a turn as the minimum
        for rot in 0..8u32 {
            let mut m = CellMesh::<3>::new();
            let vs: Vec<u32> = (0..8).map(|i| (i + rot) % 8).collect();
            // renumber: local node i has global id vs[i], so relabel points accordingly
            let mut pts = vec![[0.0; 3]; 8];
            for i in 0..8 {
                pts[vs[i] as usize] = unit_hex_points()[i];
            }
            m.points = pts;
            m.push_cell(CellType::Hex, &vs, 1);
            let s = m.split_to_simplices().unwrap();
            assert_eq!(s.n_cells(), 6);
            s.check().unwrap();
            assert!((s.total_measure() - 1.0).abs() < 1e-14, "rot {rot}");
            assert_eq!(s.boundary_keys().len(), 12);
        }
    }

    #[test]
    fn prism_splits_into_three_tets() {
        let mut m = CellMesh::<3>::new();
        m.points =
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]];
        m.push_cell(CellType::Prism, &[0, 1, 2, 3, 4, 5], 1);
        m.check().unwrap();
        let s = m.split_to_simplices().unwrap();
        assert_eq!(s.n_cells(), 3);
        s.check().unwrap();
        assert!((s.total_measure() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn two_hexes_split_conformingly_whatever_their_orientation() {
        // two unit cubes side by side sharing the x = 1 face; the second is numbered so the shared
        // face is traversed in the opposite sense and with different minimum ids.
        let mut pts = Vec::new();
        for k in 0..2 {
            for j in 0..2 {
                for i in 0..3 {
                    pts.push([i as f64, j as f64, k as f64]);
                }
            }
        }
        let id = |i: usize, j: usize, k: usize| (k * 6 + j * 3 + i) as u32;
        let mut m = CellMesh::<3>::new();
        m.points = pts;
        let h = |i0: usize| {
            [
                id(i0, 0, 0),
                id(i0 + 1, 0, 0),
                id(i0 + 1, 1, 0),
                id(i0, 1, 0),
                id(i0, 0, 1),
                id(i0 + 1, 0, 1),
                id(i0 + 1, 1, 1),
                id(i0, 1, 1),
            ]
        };
        m.push_cell(CellType::Hex, &h(0), 1);
        m.push_cell(CellType::Hex, &h(1), 1);
        m.check().unwrap();
        let s = m.split_to_simplices().unwrap();
        s.check().unwrap();
        assert_eq!(s.n_cells(), 12);
        assert!((s.total_measure() - 2.0).abs() < 1e-14);
        // outer surface of a 2x1x1 box: 2*(2*1 + 2*1 + 1*1) = 10 quads -> 20 triangles
        assert_eq!(s.boundary_keys().len(), 20);
    }

    #[test]
    fn prism_and_hex_share_a_quad_face_conformingly() {
        // unit cube (hex) with a triangular prism attached on its x = 1 face.
        let mut m = CellMesh::<3>::new();
        m.points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
        ];
        m.push_cell(CellType::Hex, &[0, 1, 2, 3, 4, 5, 6, 7], 1);
        // prism with triangle base (1, 8, 2) at z = 0 and top (5, 9, 6) at z = 1
        m.push_cell(CellType::Prism, &[1, 8, 2, 5, 9, 6], 2);
        m.check().unwrap();
        let s = m.split_to_simplices().unwrap();
        s.check().unwrap();
        assert!((s.total_measure() - 1.5).abs() < 1e-14);
    }

    #[test]
    fn quads_split_conformingly() {
        let mut m = CellMesh::<2>::new();
        m.points = vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 1.0]];
        m.push_cell(CellType::Quad, &[0, 1, 4, 3], 1);
        m.push_cell(CellType::Quad, &[1, 2, 5, 4], 1);
        m.push_facet(&[0, 1], 9);
        m.check().unwrap();
        let s = m.split_to_simplices().unwrap();
        s.check().unwrap();
        assert_eq!(s.n_cells(), 4);
        assert!((s.total_measure() - 2.0).abs() < 1e-15);
    }
}
