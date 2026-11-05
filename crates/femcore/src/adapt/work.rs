//! Simplex meshes (triangles for `D = 2`, tetrahedra for `D = 3`).
//!
//! This is the stand-in for nyx's `SimplexMesh`: flat connectivity, an integer domain label per
//! cell, and labelled boundary facets. The adaptation code only relies on the API in this file, so
//! porting into `femcore` means replacing this type by the real one behind the same accessors.
//!
//! Facets are edges in 2D and triangles in 3D. A *facet key* is the sorted vertex tuple, padded
//! with `u32::MAX`, and identifies a facet independently of orientation.

use std::collections::HashMap;

/// Sorted, padded vertex tuple identifying a facet (edge in 2D, triangle or quad in 3D).
pub type FacetKey = [u32; 4];

/// Build the orientation-independent key of a facet with 2 to 4 vertices.
pub fn facet_key(vs: &[u32]) -> FacetKey {
    debug_assert!((2..=4).contains(&vs.len()));
    let mut k = [u32::MAX; 4];
    k[..vs.len()].copy_from_slice(vs);
    k[..vs.len()].sort_unstable();
    k
}

/// Conforming simplex mesh in `D` dimensions.
#[derive(Clone, Debug, Default)]
pub struct SimplexMesh<const D: usize> {
    /// Vertex coordinates.
    pub points: Vec<[f64; D]>,
    /// Flat cell connectivity, `D + 1` vertex ids per cell.
    pub cells: Vec<u32>,
    /// Domain label of every cell.
    pub cell_label: Vec<i32>,
    /// Flat boundary-facet connectivity, `D` vertex ids per facet.
    pub facets: Vec<u32>,
    /// Label of every boundary facet.
    pub facet_label: Vec<i32>,
}

impl<const D: usize> SimplexMesh<D> {
    /// Vertices per cell.
    pub const NV: usize = D + 1;

    /// Number of vertices.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// Number of cells.
    pub fn n_cells(&self) -> usize {
        self.cell_label.len()
    }

    /// Number of labelled boundary facets.
    pub fn n_facets(&self) -> usize {
        self.facet_label.len()
    }

    /// Vertex ids of cell `c`.
    pub fn cell(&self, c: usize) -> &[u32] {
        &self.cells[c * (D + 1)..(c + 1) * (D + 1)]
    }

    /// Vertex ids of boundary facet `f`.
    pub fn facet(&self, f: usize) -> &[u32] {
        &self.facets[f * D..(f + 1) * D]
    }

    /// Append a cell.
    pub fn push_cell(&mut self, vs: &[u32], label: i32) {
        assert_eq!(vs.len(), D + 1);
        self.cells.extend_from_slice(vs);
        self.cell_label.push(label);
    }

    /// Append a boundary facet.
    pub fn push_facet(&mut self, vs: &[u32], label: i32) {
        assert_eq!(vs.len(), D);
        self.facets.extend_from_slice(vs);
        self.facet_label.push(label);
    }

    /// Signed measure (area in 2D, volume in 3D) of a simplex given by its vertices.
    pub fn signed_measure(p: &[[f64; D]]) -> f64 {
        simplex_signed_measure::<D>(p)
    }

    /// Signed measure of cell `c`.
    pub fn measure(&self, c: usize) -> f64 {
        let vs = self.cell(c);
        let mut p = [[0.0; D]; 4];
        for (i, &v) in vs.iter().enumerate() {
            p[i] = self.points[v as usize];
        }
        simplex_signed_measure::<D>(&p[..D + 1])
    }

    /// Sum of all (signed) cell measures.
    pub fn total_measure(&self) -> f64 {
        (0..self.n_cells()).map(|c| self.measure(c)).sum()
    }

    /// Swap two vertices of every negatively oriented cell so that all measures are positive.
    /// Returns the number of cells that were flipped.
    pub fn orient_positive(&mut self) -> usize {
        let mut n = 0;
        for c in 0..self.n_cells() {
            if self.measure(c) < 0.0 {
                self.cells.swap(c * (D + 1), c * (D + 1) + 1);
                n += 1;
            }
        }
        n
    }

    /// Axis-aligned bounding box `[min, max]` per axis.
    pub fn bbox(&self) -> [[f64; 2]; D] {
        bbox_of(&self.points)
    }

    /// All unique edges as sorted pairs, in sorted order (deterministic).
    pub fn edges(&self) -> Vec<[u32; 2]> {
        let mut e = Vec::with_capacity(self.n_cells() * D * (D + 1) / 2);
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
            for i in 0..=D {
                for j in i + 1..=D {
                    let (a, b) = (vs[i], vs[j]);
                    e.push([a.min(b), a.max(b)]);
                }
            }
        }
        e.sort_unstable();
        e.dedup();
        e
    }

    /// Facets of cell `c` (the face opposite vertex `i` is returned at index `i`).
    pub fn cell_facets(&self, c: usize) -> [[u32; 3]; 4] {
        let vs = self.cell(c);
        let mut out = [[u32::MAX; 3]; 4];
        for i in 0..=D {
            let mut k = 0;
            for (j, &v) in vs.iter().enumerate() {
                if j != i {
                    out[i][k] = v;
                    k += 1;
                }
            }
        }
        out
    }

    /// Facet-to-cell incidence: key of every facet and the cells that contain it.
    pub fn facet_cells(&self) -> HashMap<FacetKey, Vec<u32>> {
        let mut m: HashMap<FacetKey, Vec<u32>> = HashMap::with_capacity(self.n_cells() * 2);
        for c in 0..self.n_cells() {
            let fs = self.cell_facets(c);
            for f in fs.iter().take(D + 1) {
                m.entry(facet_key(&f[..D])).or_default().push(c as u32);
            }
        }
        m
    }

    /// Facets that belong to exactly one cell, sorted by key (deterministic).
    pub fn boundary_keys(&self) -> Vec<FacetKey> {
        let mut v: Vec<FacetKey> =
            self.facet_cells().into_iter().filter(|(_, cs)| cs.len() == 1).map(|(k, _)| k).collect();
        v.sort_unstable();
        v
    }

    /// Lumped vertex volumes: each cell gives `measure / (D + 1)` to each of its vertices.
    /// Accumulated in cell order, so the result is deterministic.
    pub fn vertex_volumes(&self) -> Vec<f64> {
        let mut vol = vec![0.0; self.n_points()];
        for c in 0..self.n_cells() {
            let m = self.measure(c).abs() / (D as f64 + 1.0);
            for &v in self.cell(c) {
                vol[v as usize] += m;
            }
        }
        vol
    }

    /// Structural and geometric validity check (Tier A properties):
    /// positive measure, every facet shared by one or two cells, every labelled boundary facet is a
    /// true boundary facet, no vertex ids out of range, and no duplicated vertices in a cell.
    pub fn check(&self) -> Result<(), String> {
        let np = self.n_points() as u32;
        if self.cells.len() != self.n_cells() * (D + 1) {
            return Err("cell array length does not match labels".into());
        }
        if self.facets.len() != self.n_facets() * D {
            return Err("facet array length does not match labels".into());
        }
        for c in 0..self.n_cells() {
            let vs = self.cell(c);
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
                return Err(format!("cell {c}: non-positive measure {m:e}"));
            }
        }
        let fc = self.facet_cells();
        for (k, cs) in &fc {
            if cs.len() > 2 {
                return Err(format!("facet {k:?} shared by {} cells", cs.len()));
            }
        }
        for f in 0..self.n_facets() {
            let key = facet_key(self.facet(f));
            match fc.get(&key) {
                Some(cs) if cs.len() == 1 => {}
                Some(_) => return Err(format!("labelled facet {f} is interior")),
                None => return Err(format!("labelled facet {f} not in any cell")),
            }
        }
        Ok(())
    }
}

/// Signed area (2D) or volume (3D) of the simplex with the given `D + 1` vertices.
pub fn simplex_signed_measure<const D: usize>(p: &[[f64; D]]) -> f64 {
    match D {
        2 => {
            let (a, b, c) = (p[0], p[1], p[2]);
            0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        }
        3 => {
            let (a, b, c, d) = (p[0], p[1], p[2], p[3]);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
            (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]))
                / 6.0
        }
        _ => panic!("only D = 2 and D = 3 are supported"),
    }
}

/// Bounding box of a point set.
pub fn bbox_of<const D: usize>(pts: &[[f64; D]]) -> [[f64; 2]; D] {
    let mut b = [[f64::INFINITY, f64::NEG_INFINITY]; D];
    for p in pts {
        for a in 0..D {
            b[a][0] = b[a][0].min(p[a]);
            b[a][1] = b[a][1].max(p[a]);
        }
    }
    b
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_square() -> SimplexMesh<2> {
        let mut m =
            SimplexMesh::<2> { points: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]], ..Default::default() };
        m.push_cell(&[0, 1, 2], 1);
        m.push_cell(&[0, 2, 3], 1);
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0)] {
            m.push_facet(&[a, b], 7);
        }
        m
    }

    #[test]
    fn square_is_valid() {
        let m = unit_square();
        m.check().unwrap();
        assert!((m.total_measure() - 1.0).abs() < 1e-15);
        assert_eq!(m.edges().len(), 5);
        assert_eq!(m.boundary_keys().len(), 4);
        let vv = m.vertex_volumes();
        assert!((vv.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn detects_inverted_and_orients() {
        let mut m = unit_square();
        m.cells.swap(0, 1);
        assert!(m.check().is_err());
        assert_eq!(m.orient_positive(), 1);
        m.check().unwrap();
    }

    #[test]
    fn detects_non_manifold_facet() {
        let mut m = unit_square();
        m.points.push([2.0, 0.5]);
        m.push_cell(&[0, 2, 4], 1); // shares edge (0,2) with both others
        assert!(m.check().is_err());
    }

    #[test]
    fn tet_volume() {
        let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!((simplex_signed_measure::<3>(&p) - 1.0 / 6.0).abs() < 1e-16);
    }
}

// ---- conversion from and to the solver's mesh (was mesh::interop) ----
use crate::mesh::SimplexMesh as NyxMesh;

/// A nyx simplex mesh of dimension `D` as the mesher's mesh.
pub fn simplex_to_adapt<const D: usize>(m: &NyxMesh) -> Result<SimplexMesh<D>, String> {
    if m.dim != D {
        return Err(format!("the mesh is {}D, not {D}D", m.dim));
    }
    let points = (0..m.n_points()).map(|i| std::array::from_fn(|a| m.points[i * D + a])).collect();
    let mut out = SimplexMesh::<D> {
        points,
        cells: m.cells.clone(),
        cell_label: m.cell_labels.clone(),
        facets: m.facets.clone(),
        facet_label: m.facet_labels.clone(),
    };
    // meshes read from files may be clockwise; the working mesh needs positive cells
    out.orient_positive();
    Ok(out)
}

/// The mesher's simplex mesh as a nyx mesh (edges and numbering of unknowns are rebuilt).
pub fn simplex_from_adapt<const D: usize>(m: &SimplexMesh<D>) -> NyxMesh {
    let points = m.points.iter().flat_map(|p| p.iter().copied()).collect();
    NyxMesh::new(D, points, m.cells.clone(), m.cell_label.clone(), m.facets.clone(), m.facet_label.clone())
}


#[cfg(test)]
mod bridge_tests {
    use super::*;
    use crate::mesh::cells::CellType;
    use crate::meshing::pml_box::pml_domains;
    use crate::props::PmlDir;
    use crate::adapt::drive::{AdaptOptions, adapt2};
    use crate::adapt::field::UniformSize;
    use crate::meshing::pml_box::{PmlSpec, build_pml};

    fn area(m: &NyxMesh) -> f64 {
        (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                let p = |i: usize| (m.points[v[i] as usize * 2], m.points[v[i] as usize * 2 + 1]);
                let (a, b, c) = (p(0), p(1), p(2));
                0.5 * ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs()
            })
            .sum()
    }

    fn volume(m: &NyxMesh) -> f64 {
        (0..m.n_cells())
            .map(|c| {
                let v = m.cell(c);
                let p = |i: usize| -> [f64; 3] { std::array::from_fn(|k| m.points[v[i] as usize * 3 + k]) };
                let (a, b, c, d) = (p(0), p(1), p(2), p(3));
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let x = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
                let det = u[0] * (w[1] * x[2] - w[2] * x[1]) - u[1] * (w[0] * x[2] - w[2] * x[0]) + u[2] * (w[0] * x[1] - w[1] * x[0]);
                det.abs() / 6.0
            })
            .sum()
    }

    #[test]
    fn a_nyx_mesh_survives_the_round_trip_and_the_adaptation() {
        let nyx = crate::meshing::geomdrive::plate(2.0, 1.0, 0.2).unwrap().mesh;
        let a = simplex_to_adapt::<2>(&nyx).unwrap();
        a.check().unwrap();
        assert!((a.total_measure().abs() - 2.0).abs() < 1e-12);
        let back = simplex_from_adapt(&a);
        assert_eq!(back.cells, nyx.cells);
        assert_eq!(back.facet_labels, nyx.facet_labels);
        assert_eq!(back.n_edges(), nyx.n_edges());
        // adapt to half the size: the area of the domain and the boundary numbers stay
        let r = adapt2(&a, &UniformSize(0.1), &AdaptOptions::default()).unwrap();
        let fine = simplex_from_adapt(&r.mesh);
        assert!(fine.n_cells() > 2 * nyx.n_cells());
        assert!((area(&fine) - 2.0).abs() < 1e-9, "area {}", area(&fine));
        let mut labels = fine.facet_labels.clone();
        labels.sort_unstable();
        labels.dedup();
        assert_eq!(labels, vec![1, 2, 3, 4]);
        assert!(simplex_to_adapt::<3>(&nyx).is_err());
    }

    #[test]
    fn a_2d_pml_mesh_splits_into_a_valid_nyx_mesh_with_matching_pml_domains() {
        let nyx = crate::meshing::geomdrive::plate(1.0, 1.0, 0.125).unwrap().mesh;
        let phys = simplex_to_adapt::<2>(&nyx).unwrap();
        let spec = PmlSpec { outer_label: Some(1), ..PmlSpec::<2>::uniform(0.25, 3) };
        let pml = build_pml(&phys, &spec).unwrap();
        let cells = &pml.mesh;
        assert!(cells.types.contains(&CellType::Quad) && cells.types.contains(&CellType::Tri) && cells.n_facets() > 0);
        let (mesh, report) = NyxMesh::from_cells_split(cells).unwrap();
        // four sides of 8 edges x 3 layers, four corners of 3 x 3
        assert_eq!(report.quads, 4 * 8 * 3 + 4 * 3 * 3);
        // the full box is 1.5 x 1.5 and every cell has positive area
        assert!((area(&mesh) - 2.25).abs() < 1e-9, "area {}", area(&mesh));
        let doms = pml_domains(&pml, 3.0, 2.0, 1.0);
        // 4 faces + 4 corners, one material
        assert_eq!(doms.len(), 8);
        let right = doms.values().filter(|p| p.x == PmlDir::Plus && p.z == PmlDir::Off).count();
        let corner = doms.values().filter(|p| p.x != PmlDir::Off && p.z != PmlDir::Off).count();
        assert_eq!((right, corner), (1, 4));
        // every domain number used by the cells is either physical (1) or a PML domain
        for l in mesh.domain_labels() {
            assert!(l == 1 || doms.contains_key(&l), "label {l}");
        }
    }

    #[test]
    fn a_3d_pml_mesh_splits_into_tetrahedra_of_the_right_volume() {
        let nyx = crate::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.34, 3).unwrap().mesh;
        let phys = simplex_to_adapt::<3>(&nyx).unwrap();
        let spec = PmlSpec { outer_label: Some(1), ..PmlSpec::<3>::uniform(0.5, 1) };
        let pml = build_pml(&phys, &spec).unwrap();
        let cells = &pml.mesh;
        assert!(cells.types.contains(&CellType::Hex) && cells.types.contains(&CellType::Prism));
        let (mesh, report) = NyxMesh::from_cells_split(cells).unwrap();
        assert!(report.hexes > 0 && report.prisms > 0);
        assert!((volume(&mesh) - 8.0).abs() < 1e-9, "volume {}", volume(&mesh));
        let doms = pml_domains(&pml, 3.0, 2.0, 1.0);
        assert_eq!(doms.len(), 26);
    }
}
