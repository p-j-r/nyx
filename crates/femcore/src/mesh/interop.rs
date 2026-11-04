//! Bridges between the nyx mesh types (`crate::mesh::SimplexMesh`, `mesh::cells::CellMesh`) and the dimension-generic
//! types the ported mesher works on (`mesh::native::SimplexMesh<D>`, `mesh::cells::CellMesh<D>`).
//!
//! The mesher was written against stand-ins for the nyx types (see `docs/ports/meshers/HANDOFF.md`); rather than rewrite
//! ten thousand lines of tested adaptation code onto flat arrays, the stand-ins stay as the mesher's working types and
//! these functions translate at the boundary:
//!
//! * `simplex_to_adapt` / `simplex_from_adapt`: triangles and tetrahedra with their domain and boundary numbers. The
//!   numbers are the same 1-based labels on both sides.
//! * `pml_domains`: the nyx `Pml` description (`ProblemSpec::optical_pml`) of every PML domain label a `PmlMesh` uses. nyx
//!   stretches along the bounding box of a domain's cells, which is exactly how the mesher cuts its regions, so the
//!   direction of the stretch is the sign of the region code.

use crate::mesh::SimplexMesh as NyxMesh;
use crate::meshing::pml_box::PmlMesh;
use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::native::SimplexMesh;
use crate::problem::{Pml, PmlDir};
use std::collections::BTreeMap;

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

/// The `Pml` of every domain label of a PML mesh, keyed by label, to merge into `ProblemSpec::optical_pml`.
/// `strength`, `order` and `re` are the stretch parameters of `Pml`.
pub fn pml_domains<const D: usize>(pml: &PmlMesh<D>, strength: f64, order: f64, re: f64) -> BTreeMap<i32, Pml> {
    let dir = |c: i8| match c {
        0 => PmlDir::Off,
        c if c > 0 => PmlDir::Plus,
        _ => PmlDir::Minus,
    };
    let mut out = BTreeMap::new();
    for r in &pml.regions {
        let axis = |a: usize| if a < D { dir(r.code[a]) } else { PmlDir::Off };
        for &(_, label) in &r.labels {
            // a 2D mesh is the x-z cross-section: its second axis is z
            let (y, z) = if D == 2 { (PmlDir::Off, axis(1)) } else { (axis(1), axis(2)) };
            out.insert(label, Pml { x: axis(0), y, z, strength, order, re });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
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
