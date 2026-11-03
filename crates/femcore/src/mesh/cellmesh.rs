//! A mesh as a list of cells of mixed type (segments, triangles, quadrilaterals, tetrahedra, hexahedra, prisms), the
//! form in which meshes from outside programs and the structured test generators arrive. `SimplexMesh::from_cells`
//! turns it into the triangles or tetrahedra the solvers use, cutting other cells without adding vertices.

#[derive(Clone, Debug, Default)]
pub struct ElementType {
    pub name: String,
    /// vertices per element
    pub nep: usize,
    /// flat vertex indices, `nep` per element
    pub elem: Vec<u32>,
    /// geometric entity index of each element (0-based for boundaries, 1-based for domains:
    /// entity 0 of the highest dimension is the exterior)
    pub dom: Vec<i32>,
}

impl ElementType {
    pub fn len(&self) -> usize {
        self.elem.len().checked_div(self.nep).unwrap_or(0)
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn row(&self, i: usize) -> &[u32] {
        &self.elem[i * self.nep..(i + 1) * self.nep]
    }
}

#[derive(Clone, Debug, Default)]
pub struct CellMesh {
    pub sdim: usize,
    /// flat coordinates, `sdim` per point
    pub points: Vec<f64>,
    pub types: Vec<ElementType>,
}

impl CellMesh {
    pub fn n_points(&self) -> usize {
        self.points.len().checked_div(self.sdim).unwrap_or(0)
    }
    pub fn ty(&self, name: &str) -> Option<&ElementType> {
        self.types.iter().find(|t| t.name == name)
    }
}
