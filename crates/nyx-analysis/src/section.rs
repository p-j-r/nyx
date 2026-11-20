// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The cross-section of a result as a mesh of triangles with the nodes the cell fields live on: the whole mesh for a 2D result, the face of the
//! 3D cells that lies on a plane of constant propagation coordinate for a 3D periodic slab (a mode of a Bloch study along y). Plain functions of
//! the mesh and the order of the result; no physics is named.
//!
//! The cell node order is that of nyx (`cell_order` of a results file): P1 `v0 v1 v2`, P2 `v0 v1 v2 e01 e02 e12`, and for a tetrahedron
//! `v0 v1 v2 v3 e01 e02 e03 e12 e13 e23`. A face `(i, j, k)` of a tetrahedron, `i < j < k`, becomes the triangle `vi vj vk eij eik ejk`, which is
//! again in the order of a P2 triangle.

use femcore::mesh::SimplexMesh;

/// Local edge pairs of a tetrahedron in the order of its P2 edge nodes.
const TET_PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

/// What the triangles say of one edge: the label of the first, how many hold it, whether two labels met on it, its midpoint node.
type EdgeUse = (i32, u32, bool, Option<u32>);

#[derive(Clone, Debug)]
pub struct Section {
    /// 1 or 2: the order the cell fields were solved with
    pub order: usize,
    /// nodes per triangle: 3 (P1) or 6 (P2)
    pub nt: usize,
    /// node coordinates (lateral x, normal z) in the wafer frame
    pub xz: Vec<[f64; 2]>,
    /// the nodes of every triangle, `nt` per triangle, in the cell node order above
    pub tris: Vec<u32>,
    /// the cell label (domain number) of every triangle
    pub labels: Vec<i32>,
    /// for every triangle node, where the cell fields keep its value: `cell * nloc + local position`
    pub src: Vec<u32>,
    /// the plane the section lies on (`None` for a 2D mesh)
    pub at: Option<f64>,
}

impl Section {
    pub fn n_tris(&self) -> usize {
        self.labels.len()
    }

    pub fn n_nodes(&self) -> usize {
        self.xz.len()
    }

    /// The nodes of triangle `t`.
    pub fn tri(&self, t: usize) -> &[u32] {
        &self.tris[t * self.nt..(t + 1) * self.nt]
    }

    /// The section of `mesh` for a result of `order`. A 3D mesh is cut on the plane y = `at` (default: its smallest y); the cells that have a
    /// face on that plane give the triangles. A section in the interior of the cells would need interpolation and is not offered.
    pub fn new(mesh: &SimplexMesh, order: usize, at: Option<f64>) -> Result<Section, String> {
        let order = if order == 2 { 2 } else { 1 };
        let dim = mesh.dim;
        if dim != 2 && dim != 3 {
            return Err(format!("a section needs a 2D or 3D mesh, this one is {dim}D"));
        }
        let nloc = mesh.nloc(order);
        let coords = mesh.dof_coordinates(order);
        let nt = if order == 2 { 6 } else { 3 };
        let mut dofs = Vec::new();
        let mut node_of = vec![u32::MAX; mesh.n_dofs(order)];
        let mut sec = Section { order, nt, xz: vec![], tris: vec![], labels: vec![], src: vec![], at: None };
        let mut push = |sec: &mut Section, cell: usize, dofs: &[usize], pos: &[usize]| {
            for &p in pos {
                let d = dofs[p];
                if node_of[d] == u32::MAX {
                    node_of[d] = sec.xz.len() as u32;
                    let c = &coords[d * dim..(d + 1) * dim];
                    sec.xz.push(if dim == 2 { [c[0], c[1]] } else { [c[0], c[2]] });
                }
                sec.tris.push(node_of[d]);
                sec.src.push((cell * nloc + p) as u32);
            }
            sec.labels.push(mesh.cell_labels[cell]);
        };
        if dim == 2 {
            let pos: Vec<usize> = (0..nloc).collect();
            for c in 0..mesh.n_cells() {
                mesh.cell_dofs(c, order, &mut dofs);
                push(&mut sec, c, &dofs, &pos);
            }
        } else {
            let ys = (0..mesh.n_points()).map(|v| mesh.points[v * 3 + 1]);
            let (ymin, ymax) = ys.fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), y| (a.min(y), b.max(y)));
            let y0 = at.unwrap_or(ymin);
            let tol = 1e-9 * (ymax - ymin).max(f64::MIN_POSITIVE) + 1e-18;
            sec.at = Some(y0);
            for c in 0..mesh.n_cells() {
                let verts = mesh.cell(c);
                let on: Vec<usize> =
                    (0..4).filter(|&i| (mesh.points[verts[i] as usize * 3 + 1] - y0).abs() <= tol).collect();
                if on.len() != 3 {
                    continue;
                }
                let e = |i: usize, j: usize| 4 + TET_PAIRS.iter().position(|&p| p == (i.min(j), i.max(j))).unwrap();
                let (i, j, k) = (on[0], on[1], on[2]);
                let pos_all = [i, j, k, e(i, j), e(i, k), e(j, k)];
                mesh.cell_dofs(c, order, &mut dofs);
                push(&mut sec, c, &dofs, &pos_all[..nt]);
            }
            if sec.n_tris() == 0 {
                return Err(format!("no cell has a face on the plane y = {y0}; give `at` as the y of a mesh plane"));
            }
        }
        Ok(sec)
    }

    /// The triangles to draw: every P2 triangle split in four linear ones (`v0 e01 e02`, `e01 v1 e12`, `e02 e12 v2`, `e01 e12 e02`), which
    /// shows a quadratic field exactly at its nodes; a P1 triangle as it is.
    pub fn linear_tris(&self) -> Vec<[u32; 3]> {
        let mut out = Vec::with_capacity(self.n_tris() * if self.nt == 6 { 4 } else { 1 });
        for t in 0..self.n_tris() {
            let n = self.tri(t);
            if self.nt == 6 {
                out.push([n[0], n[3], n[4]]);
                out.push([n[3], n[1], n[5]]);
                out.push([n[4], n[5], n[2]]);
                out.push([n[3], n[5], n[4]]);
            } else {
                out.push([n[0], n[1], n[2]]);
            }
        }
        out
    }

    /// The label of every linear triangle of `linear_tris`.
    pub fn linear_labels(&self) -> Vec<i32> {
        let k = if self.nt == 6 { 4 } else { 1 };
        self.labels.iter().flat_map(|&l| std::iter::repeat_n(l, k)).collect()
    }

    /// The outline of the regions: the edges (as segments of two nodes, through the edge midpoint of a P2 triangle) that lie on the border of the
    /// section or between two cell labels.
    pub fn outline(&self) -> Vec<[u32; 2]> {
        use std::collections::HashMap;
        // key: the two vertex nodes of an edge; value: the label of the first triangle, how many triangles hold the edge, whether two labels
        // met on it, the midpoint node
        let mut edges: HashMap<(u32, u32), EdgeUse> = HashMap::new();
        for t in 0..self.n_tris() {
            let n = self.tri(t);
            let l = self.labels[t];
            let mid = |a: usize| if self.nt == 6 { Some(n[a]) } else { None };
            for (a, b, m) in [(0usize, 1usize, mid(3)), (0, 2, mid(4)), (1, 2, mid(5))] {
                let key = (n[a].min(n[b]), n[a].max(n[b]));
                let ent = edges.entry(key).or_insert((l, 0, false, m));
                ent.1 += 1;
                if ent.0 != l {
                    ent.2 = true;
                }
            }
        }
        let mut keys: Vec<_> = edges.into_iter().filter(|(_, v)| v.1 == 1 || v.2).collect();
        keys.sort_by_key(|(k, _)| *k);
        let mut out = vec![];
        for ((a, b), v) in keys {
            match v.3 {
                Some(m) => {
                    out.push([a, m]);
                    out.push([m, b]);
                }
                None => out.push([a, b]),
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit square cut in two triangles.
    fn square() -> SimplexMesh {
        SimplexMesh::new(
            2,
            vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0],
            vec![0, 1, 2, 0, 2, 3],
            vec![1, 2],
            vec![],
            vec![],
        )
    }

    #[test]
    fn a_2d_mesh_is_its_own_section_and_the_outline_separates_labels() {
        let s = Section::new(&square(), 2, None).unwrap();
        assert_eq!((s.n_tris(), s.nt, s.n_nodes()), (2, 6, 9));
        assert_eq!(s.linear_tris().len(), 8);
        // four border edges and the shared diagonal that separates the two labels: 5 edges, 10 segments through the midpoints
        assert_eq!(s.outline().len(), 10);
        let t = s.tri(0);
        let (a, b, m) = (s.xz[t[0] as usize], s.xz[t[1] as usize], s.xz[t[3] as usize]);
        assert!((m[0] - 0.5 * (a[0] + b[0])).abs() < 1e-15 && (m[1] - 0.5 * (a[1] + b[1])).abs() < 1e-15);
    }

    #[test]
    fn a_3d_slab_is_cut_on_its_first_plane() {
        // one tetrahedron with the face (0, 1, 2) in the plane y = 0 and the apex above it
        let pts = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let m = SimplexMesh::new(3, pts, vec![0, 1, 2, 3], vec![1], vec![], vec![]);
        let s = Section::new(&m, 2, None).unwrap();
        assert_eq!((s.n_tris(), s.n_nodes()), (1, 6));
        assert_eq!(s.at, Some(0.0));
        let xs = &s.xz;
        assert!(xs.contains(&[0.5, 0.0]) && xs.contains(&[0.0, 0.5]) && xs.contains(&[0.5, 0.5]));
    }
}
