//! The wafer frame: x-y plane of the wafer, z the surface normal, a 2D cross-section is the x-z plane.
mod common;
use common::*;
use femcore::material::*;
use femcore::problem::*;

#[test]
fn a_permittivity_tensor_in_the_wafer_frame_acts_on_the_vertical_axis_of_the_x_z_plane() {
    // plates on the bottom (z = 0) and top (z = H) of an x-z section of length L: C = eps0 eps_zz L / H per metre
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.aniso.insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }), ..Default::default() });
    spec.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    let r = femcore::solve::solve(&mesh, &spec.clone(), &mut SolveState::default()).unwrap();
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (femcore::dfs::EPS0 * 5.0 * 1.0 / 0.5) - 1.0).abs() < 1e-9, "{c}");
    // the field of a uniform vertical field is called Ez
    assert!(r.solutions[0].field("Ez").is_some() && r.solutions[0].field("Ey").is_none());
    // plates on the sides see eps_xx
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    let r = femcore::solve::solve(&mesh, &spec, &mut SolveState::default()).unwrap();
    assert!((r.scalar("capacitance").unwrap() / (femcore::dfs::EPS0 * 7.0 * 0.5 / 1.0) - 1.0).abs() < 1e-9);
}

#[test]
fn a_loaded_cubic_plate_in_the_x_z_plane_feels_a_turn_of_the_crystal_about_y_but_not_about_z_in_the_same_way() {
    // wafer description: fixed at x = 0 (x and z components), a load along -z on the right end; the crystal is turned by 20 degrees about y
    // (a rotation in the plane of the section) or by 20 degrees about the normal z (a rotation that mixes the in-plane axis x with the
    // axis y out of the section): the two are different crystals and give different deflections
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let rot_y = |t: f64| -> M3 {
        let (s, c) = t.to_radians().sin_cos();
        [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]
    };
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.plane = Plane::Stress;
    spec.materials.set(1, "E", 1e9);
    spec.materials.set(1, "nu", 0.3);
    spec.materials.aniso.insert(1, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: euler_from_rotation(&rot_y(20.0)), ..Default::default() });
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } }, Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e6]) }];
    let tip = |spec: &ProblemSpec| -> f64 {
        let r = femcore::solve::solve(&mesh, spec, &mut SolveState::default()).unwrap();
        r.solutions[0].field("uz").unwrap().values.iter().fold(0.0f64, |a, b| a.max(b.abs()))
    };
    let about_y = tip(&spec);
    spec.materials.aniso.get_mut(&1).unwrap().euler_deg = [20.0, 0.0, 0.0];
    let about_z = tip(&spec);
    spec.materials.aniso.get_mut(&1).unwrap().euler_deg = [0.0; 3];
    let straight = tip(&spec);
    assert!(about_y > 0.0 && (about_y / straight - 1.0).abs() > 1e-3, "{about_y} against {straight}");
    // the turn about z is not the turn about y
    assert!((about_z / about_y - 1.0).abs() > 1e-4, "{about_z} against {about_y}");
}

#[test]
fn the_chip_crystals_are_oriented_in_the_wafer_frame() {
    use femcore::chips::*;
    let b = 57.6f64.to_radians();
    // psi = 0: the sapphire c axis lies in the y-z plane (y = propagation, up to sign), 57.6 degrees from the surface normal z
    let (r_si, r_sap) = crystal_rotations(0.0);
    let c_axis = [r_sap[0][2], r_sap[1][2], r_sap[2][2]];
    assert!(c_axis[0].abs() < 1e-12 && (c_axis[1].abs() - b.sin()).abs() < 1e-12 && (c_axis[2] - b.cos()).abs() < 1e-12, "{c_axis:?}");
    let a_axis = [r_sap[0][0], r_sap[1][0], r_sap[2][0]];
    assert!((a_axis[0] - 1.0).abs() < 1e-12, "{a_axis:?}");
    // silicon [100] along the propagation direction y, [010] along x, [001] along z
    assert!((r_si[1][0].abs() - 1.0).abs() < 1e-12 && (r_si[0][1] - 1.0).abs() < 1e-12 && (r_si[2][2] - 1.0).abs() < 1e-12);
    // psi = 90: the projection of the c axis lies along x
    let (_, r90) = crystal_rotations(90.0);
    assert!((r90[0][2] - b.sin()).abs() < 1e-12 && r90[1][2].abs() < 1e-12);
    // every rotation is proper
    for psi in [0.0, 30.0, 45.0, 90.0] {
        let (s, p) = crystal_rotations(psi);
        for r in [s, p] {
            let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
            assert!((det - 1.0).abs() < 1e-12);
        }
    }
}

#[test]
fn a_layer_along_y_does_not_exist_on_a_2d_cross_section() {
    // layers of a 2D cross-section are written along x and z; there is no y axis in the plane
    let mesh = common::rect(1.0, 0.5, 4, 2);
    let mut bad = ProblemSpec::new(Physics::Solid);
    bad.study = Study::Eigenfrequency { n_modes: 2, f_search: 0.0 };
    bad.materials.pml.insert(1, Pml { y: PmlDir::Plus, ..Default::default() });
    let e = femcore::solve::solve(&mesh, &bad, &mut SolveState::default()).err().unwrap();
    assert!(e.contains("absorbs along y") && e.contains("use z"), "{e}");
    // the Brillouin optical layers follow the same rule
    let mut b = ProblemSpec::new(Physics::Brillouin);
    b.study = Study::Eigenfrequency { n_modes: 2, f_search: 1e9 };
    b.brillouin.as_mut().unwrap().optical_pml.insert(1, Pml { y: PmlDir::Plus, ..Default::default() });
    let e = femcore::solve::solve(&mesh, &b, &mut SolveState::default()).err().unwrap();
    assert!(e.contains("absorbs along y") && e.contains("use z"), "{e}");
}

#[test]
fn the_absorbing_pieces_of_a_geom_model_come_out_in_wafer_names() {
    use femcore::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
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

