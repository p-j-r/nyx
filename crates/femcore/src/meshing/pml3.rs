//! Structured absorbing layers around a tetrahedral mesh in a box: prisms on the faces, hexahedra on the edges and corners, the labels
//! and `Pml` entries the solvers read, and the simplex mesh they run on today.
//!
//! `add_pml_3d` takes any conforming tetrahedral mesh whose boundary lies on the faces of its bounding box (a block, or such a block
//! adapted by `physics_size::adapt_tets`, which keeps the box faces, edges and corners exactly) and sweeps the layers outward from it
//! (`meshing::pml_box::build_pml`). The result has
//!
//! * `layout.mesh`: the mixed mesh as a `mesh::cells::CellMesh` (tetrahedra, prisms, hexahedra, boundary triangles and quadrilaterals), for
//!   solvers that take those cells natively;
//! * `mesh`: the same mesh cut into tetrahedra without new vertices (`SimplexMesh::from_cells_split`), what the solvers of today run on;
//! * `pml`: the `Pml` stretch of every layer domain, keyed by its label, one domain per region (a face has one stretched direction, an
//!   edge two, a corner three) because `femcore::pml::cell_stretches` takes the bounding box of a domain's cells.
//!
//! A material that crosses the interface continues into a PML copy of itself; `apply_materials` copies the material data of the physical
//! domains to the layer labels and registers the `Pml` entries, so a `ProblemSpec` for the physical mesh becomes one for the full mesh.

use crate::mesh::{SimplexMesh, SplitReport};
use crate::meshing::pml_box::{build_pml, PmlMesh, PmlSpec};
use crate::adapt::work::simplex_to_adapt;
use crate::meshing::pml_box::pml_domains;
use crate::props::{Materials, Pml};
use std::collections::BTreeMap;

pub struct PmlBox {
    /// the same mesh cut into tetrahedra
    pub mesh: SimplexMesh,
    pub split: SplitReport,
    /// stretch of every layer domain, by label
    pub pml: BTreeMap<i32, Pml>,
    /// the generator's description: regions, labels, interface box
    pub layout: PmlMesh<3>,
}

/// Adds structured layers to the tetrahedral mesh `phys`. `strength`, `order` and `re` are the stretch `s = re - j strength t^order` of
/// every layer (see `problem::Pml`); `spec` has the thickness, number of layers and label scheme per side. The label scheme must not
/// collide with the labels of `phys` (the default starts at 1000).
pub fn add_pml_3d(phys: &SimplexMesh, spec: &PmlSpec<3>, strength: f64, order: f64, re: f64) -> Result<PmlBox, String> {
    if phys.dim != 3 {
        return Err("add_pml_3d needs a tetrahedral mesh".into());
    }
    let a = simplex_to_adapt::<3>(phys)?;
    let used: Vec<i32> = phys.cell_labels.clone();
    if used.iter().any(|&l| l >= spec.label_base) {
        return Err(format!("a physical domain label is not below the PML label base {}", spec.label_base));
    }
    let layout = build_pml(&a, spec)?;
    let (mesh, split) = SimplexMesh::from_cells_split(&layout.mesh)?;
    let pml = pml_domains(&layout, strength, order, re);
    Ok(PmlBox { mesh, split, pml, layout })
}

/// Gives every layer domain the material data of the physical domain it continues and its `Pml` entry. Domains of `materials` that are not
/// listed by label (only `default` values) need nothing.
pub fn apply_materials(pb: &PmlBox, materials: &mut Materials) {
    for r in &pb.layout.regions {
        for &(material, label) in &r.labels {
            if let Some(m) = materials.domains.get(&material).cloned() {
                materials.domains.insert(label, m);
            }
            if let Some(a) = materials.aniso.get(&material).cloned() {
                materials.aniso.insert(label, a);
            }
            if let Some(p) = pb.pml.get(&label) {
                materials.pml.insert(label, p.clone());
            }
        }
    }
}
