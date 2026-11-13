//! The wafer frame: x-y plane of the wafer, z the surface normal, a 2D cross-section is the x-z plane.
mod common;
use femcore::props::*;

#[test]
fn the_absorbing_pieces_of_a_geom_model_come_out_in_wafer_names() {
    use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
    use geom::model::Model;
    let mut m = Model::new_2d();
    let air = m.rect("air", (-2.0, -1.0), (2.0, 1.0));
    let win = m.rect("window", (-1.5, -0.5), (1.5, 0.8));
    m.overlay(&[("air", air)]);
    m.pml("pml", win, &[], 0.5, 4.0, 2.0);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.25)).unwrap();
    let pml = g.pml_materials(3.0);
    let corner = pml[&g.domain("air.pml_xmin_zmin").unwrap()].clone();
    assert_eq!((corner.x, corner.y, corner.z), (PmlDir::Minus, PmlDir::Off, PmlDir::Minus));
    let top = pml[&g.domain("air.pml_zmax").unwrap()].clone();
    assert_eq!((top.x, top.y, top.z), (PmlDir::Off, PmlDir::Off, PmlDir::Plus));
}

