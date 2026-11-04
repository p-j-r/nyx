//! The mesh a physics is given (`Physics::solve`): the one native cell mesh (`cells::CellMesh`) in 2D or in 3D, as plain data. It replaces
//! `SimplexMesh` as the exchange type of the physics interface (ADR 32 phase 4 deliverable, ADR 41). Each solve turns it into the numbering of
//! unknowns it assembles on with `Mesh::discretise` (edges, P2 dofs; non-simplex cells are cut into simplices without new vertices, as
//! `SimplexMesh::from_cells_split` always did, so the numbers of every solution are unchanged).

use crate::mesh::cells::{CellMesh, CellType};
use crate::mesh::SimplexMesh;

/// A mesh in 2D (a cross-section, the x-z plane of the wafer frame) or in 3D.
#[derive(Clone, Debug)]
pub enum Mesh {
    D2(CellMesh<2>),
    D3(CellMesh<3>),
}

fn cells_of<const D: usize>(m: &SimplexMesh) -> CellMesh<D> {
    let mut out = CellMesh::<D>::new();
    out.points = m
        .points
        .chunks(D)
        .map(|p| {
            let mut q = [0.0; D];
            q.copy_from_slice(p);
            q
        })
        .collect();
    let t = if D == 2 { CellType::Tri } else { CellType::Tet };
    for c in 0..m.n_cells() {
        out.push_cell(t, m.cell(c), m.cell_labels[c]);
    }
    for f in 0..m.n_facets() {
        out.push_facet(m.facet(f), m.facet_labels[f]);
    }
    out
}

impl Mesh {
    /// The exchange form of a solver mesh: same vertices, cells, facets and numbers, in the same order.
    pub fn from_solver_mesh(m: &SimplexMesh) -> Result<Mesh, String> {
        match m.dim {
            2 => Ok(Mesh::D2(cells_of::<2>(m))),
            3 => Ok(Mesh::D3(cells_of::<3>(m))),
            d => Err(format!("mesh dimension {d} is not supported")),
        }
    }

    /// A simplex mesh from flat arrays (as `SimplexMesh::new`: coordinates `dim` per vertex, `dim + 1` vertices per cell, `dim` per facet, one
    /// number per cell and per facet). For small meshes written out in tests; panics when `dim` is not 2 or 3.
    pub fn new(dim: usize, points: Vec<f64>, cells: Vec<u32>, cell_labels: Vec<i32>, facets: Vec<u32>, facet_labels: Vec<i32>) -> Mesh {
        Mesh::from_solver_mesh(&SimplexMesh::new(dim, points, cells, cell_labels, facets, facet_labels)).expect("dimension 2 or 3")
    }

    /// 2 or 3.
    pub fn dim(&self) -> usize {
        match self {
            Mesh::D2(_) => 2,
            Mesh::D3(_) => 3,
        }
    }

    pub fn n_cells(&self) -> usize {
        match self {
            Mesh::D2(m) => m.n_cells(),
            Mesh::D3(m) => m.n_cells(),
        }
    }

    pub fn n_points(&self) -> usize {
        match self {
            Mesh::D2(m) => m.n_points(),
            Mesh::D3(m) => m.n_points(),
        }
    }

    /// The first vertex and axis whose coordinate is NaN or infinite.
    pub fn first_non_finite(&self) -> Option<(usize, usize)> {
        fn find<const D: usize>(m: &CellMesh<D>) -> Option<(usize, usize)> {
            m.points.iter().enumerate().find_map(|(i, p)| p.iter().position(|v| !v.is_finite()).map(|k| (i, k)))
        }
        match self {
            Mesh::D2(m) => find(m),
            Mesh::D3(m) => find(m),
        }
    }

    /// The numbering of unknowns a kernel assembles on. Cells that are not simplices are cut (`SimplexMesh::from_cells_split`).
    pub fn discretise(&self) -> Result<SimplexMesh, String> {
        match self {
            Mesh::D2(m) => SimplexMesh::from_cells(m),
            Mesh::D3(m) => SimplexMesh::from_cells(m),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_solver_mesh_survives_the_round_trip_through_the_exchange_mesh() {
        for m in [crate::testkit::rect(1.0, 0.5, 5, 3), crate::testkit::block(1.0, 0.5, 0.4, 3, 2, 2)] {
            let back = Mesh::from_solver_mesh(&m).unwrap().discretise().unwrap();
            assert_eq!((back.dim, &back.points, &back.cells), (m.dim, &m.points, &m.cells));
            assert_eq!((&back.cell_labels, &back.facets, &back.facet_labels), (&m.cell_labels, &m.facets, &m.facet_labels));
            assert_eq!(back.edges(), m.edges());
        }
    }
}
