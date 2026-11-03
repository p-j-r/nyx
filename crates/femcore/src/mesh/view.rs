//! Draw-only mesh for display: works for simplex meshes and for meshes the solvers cannot use
//! (quadrilaterals, hexahedra, prisms...), of which only the surface is kept.

use crate::mesh::SimplexMesh;
use crate::cellmesh::CellMesh;
use std::collections::BTreeSet;

#[derive(Clone, Debug, Default)]
pub struct ViewMesh {
    pub dim: usize,
    pub points: Vec<[f32; 3]>,
    /// surface triangles: the cells in 2D, the boundary facets in 3D
    pub tris: Vec<[u32; 3]>,
    /// domain (2D) or boundary (3D) number of each surface triangle
    pub tri_label: Vec<i32>,
    /// boundary segments (2D only) with their boundary numbers
    pub segs: Vec<[u32; 2]>,
    pub seg_label: Vec<i32>,
    pub domain_labels: Vec<i32>,
    pub boundary_labels: Vec<i32>,
    pub n_cells: usize,
    pub bbox: ([f32; 3], [f32; 3]),
    /// element type names present in the source mesh
    pub types: Vec<String>,
    /// why the solvers cannot use this mesh (None for a simplex mesh)
    pub view_only: Option<String>,
}

fn finish(mut v: ViewMesh) -> ViewMesh {
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for p in &v.points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    if v.points.is_empty() {
        lo = [0.0; 3];
        hi = [1.0; 3];
    }
    v.bbox = (lo, hi);
    v
}

impl ViewMesh {
    pub fn from_simplex(m: &SimplexMesh) -> ViewMesh {
        let d = m.dim;
        let points: Vec<[f32; 3]> = m
            .points
            .chunks(d)
            .map(|p| {
                let mut q = [0f32; 3];
                for k in 0..d {
                    q[k] = p[k] as f32;
                }
                q
            })
            .collect();
        let (tris, tri_label, segs, seg_label) = if d == 2 {
            (
                m.cells.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(),
                m.cell_labels.clone(),
                m.facets.chunks(2).map(|c| [c[0], c[1]]).collect(),
                m.facet_labels.clone(),
            )
        } else {
            (m.facets.chunks(3).map(|c| [c[0], c[1], c[2]]).collect(), m.facet_labels.clone(), vec![], vec![])
        };
        finish(ViewMesh {
            dim: d,
            points,
            tris,
            tri_label,
            segs,
            seg_label,
            domain_labels: m.domain_labels(),
            boundary_labels: m.boundary_labels(),
            n_cells: m.n_cells(),
            bbox: ([0.0; 3], [0.0; 3]),
            types: vec![if d == 2 { "tri".into() } else { "tet".into() }],
            view_only: None,
        })
    }

    pub fn from_cells(nm: &CellMesh) -> ViewMesh {
        let d = nm.sdim;
        let np = nm.n_points();
        let points: Vec<[f32; 3]> = (0..np)
            .map(|i| {
                let mut q = [0f32; 3];
                for k in 0..d.min(3) {
                    q[k] = nm.points[i * d + k] as f32;
                }
                q
            })
            .collect();
        let mut tris = vec![];
        let mut tri_label = vec![];
        let mut segs = vec![];
        let mut seg_label = vec![];
        // domains are 1-based in the file, boundaries 0-based
        let shift = if d == 2 { 0 } else { 1 };
        for t in &nm.types {
            if t.name == "tri" || t.name == "quad" {
                for i in 0..t.len() {
                    let r = t.row(i);
                    let lab = t.dom.get(i).copied().unwrap_or(0) + shift;
                    tris.push([r[0], r[1], r[2]]);
                    tri_label.push(lab);
                    if t.name == "quad" {
                        tris.push([r[0], r[2], r[3]]);
                        tri_label.push(lab);
                    }
                }
            } else if t.name == "edg" && d == 2 {
                for i in 0..t.len() {
                    let r = t.row(i);
                    segs.push([r[0], r[1]]);
                    seg_label.push(t.dom.get(i).copied().unwrap_or(0) + 1);
                }
            }
        }
        let cell_names: &[&str] = if d == 2 { &["tri", "quad"] } else { &["tet", "hex", "prism", "pyr"] };
        let n_cells = nm.types.iter().filter(|t| cell_names.contains(&t.name.as_str())).map(|t| t.len()).sum();
        let types: Vec<String> = nm.types.iter().map(|t| t.name.clone()).collect();
        let domain_labels: Vec<i32> = if d == 2 {
            tri_label.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
        } else {
            nm.types
                .iter()
                .filter(|t| cell_names.contains(&t.name.as_str()))
                .flat_map(|t| t.dom.iter().copied())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect()
        };
        let boundary_labels: Vec<i32> = if d == 2 {
            seg_label.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
        } else {
            tri_label.iter().copied().collect::<BTreeSet<_>>().into_iter().collect()
        };
        let known: &[&str] = if d == 2 { &["vtx", "edg", "tri", "quad"] } else { &["vtx", "edg", "tri", "quad", "tet", "hex", "prism"] };
        let odd: Vec<&str> = types.iter().map(|s| s.as_str()).filter(|n| !known.contains(n)).collect();
        let view_only = if odd.is_empty() { None } else { Some(format!("This mesh has {} elements, which are not supported; it is shown for inspection only.", odd.join(", "))) };
        finish(ViewMesh {
            dim: d,
            points,
            tris,
            tri_label,
            segs,
            seg_label,
            domain_labels,
            boundary_labels,
            n_cells,
            bbox: ([0.0; 3], [0.0; 3]),
            types,
            view_only,
        })
    }
}
