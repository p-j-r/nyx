//! The wafer frame: x-y plane of the wafer, z the surface normal, a 2D cross-section is the x-z plane.
mod common;
use common::*;
use femcore::material::*;
use femcore::mesh::SimplexMesh;
use femcore::problem::*;

#[test]
fn a_permittivity_tensor_in_the_wafer_frame_acts_on_the_vertical_axis_of_the_x_z_plane() {
    // plates on the bottom (z = 0) and top (z = H) of an x-z section of length L: C = eps0 eps_zz L / H per metre
    let mesh = SimplexMesh::rectangle(1.0, 0.5, 8, 4);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.aniso.insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }), ..Default::default() });
    spec.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    let r = femcore::solve::solve(&mesh, &spec.clone(), &mut SolveState::default()).unwrap();
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (femcore::scalar::EPS0 * 5.0 * 1.0 / 0.5) - 1.0).abs() < 1e-9, "{c}");
    // the field of a uniform vertical field is called Ez
    assert!(r.solutions[0].field("Ez").is_some() && r.solutions[0].field("Ey").is_none());
    // plates on the sides see eps_xx
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }];
    let r = femcore::solve::solve(&mesh, &spec, &mut SolveState::default()).unwrap();
    assert!((r.scalar("capacitance").unwrap() / (femcore::scalar::EPS0 * 7.0 * 0.5 / 1.0) - 1.0).abs() < 1e-9);
}

#[test]
fn the_stiffness_of_a_tilted_crystal_in_the_model_frame_is_the_wafer_tensor_turned_by_the_frame_change() {
    // sapphire with its c axis tilted by 57.6 degrees from the normal z towards y; the crystal axes in the wafer frame are
    // a = x, c x a, c; the model axes are (x, z, -y). The 6x6 matrix is compared with a direct fourth-rank transform.
    let b = 57.6f64.to_radians();
    let r_lab: M3 = [[1.0, 0.0, 0.0], [0.0, b.cos(), b.sin()], [0.0, -b.sin(), b.cos()]];
    let eul = euler_from_rotation(&r_lab);
    assert!({
        let r = rotation(eul);
        (0..3).all(|i| (0..3).all(|j| (r[i][j] - r_lab[i][j]).abs() < 1e-12))
    });
    let st = Stiffness::Trigonal { c11: 497.35e9, c12: 163.97e9, c13: 112.20e9, c14: -23.58e9, c33: 499.11e9, c44: 147.39e9 };
    let c_lab = stiffness_global(&st, eul); // solver Voigt order xx yy zz xy yz xz, wafer axes
    let s: M3 = femcore::frame_adapter::WAFER_TO_MODEL;
    // rank 4 from the 6x6
    let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, bb)) in vp.iter().enumerate() {
        for (j, &(c, d)) in vp.iter().enumerate() {
            for (p, q) in [(a, bb), (bb, a)] {
                for (u, v) in [(c, d), (d, c)] {
                    t[p][q][u][v] = c_lab[i][j];
                }
            }
        }
    }
    let at = |i: usize, j: usize, k: usize, l: usize| -> f64 {
        let mut acc = 0.0;
        for a in 0..3 {
            for bb in 0..3 {
                for c in 0..3 {
                    for d in 0..3 {
                        acc += s[i][a] * s[j][bb] * s[k][c] * s[l][d] * t[a][bb][c][d];
                    }
                }
            }
        }
        acc
    };
    let c_model = stiffness_global(&st, femcore::frame_adapter::euler_to_model_frame(eul));
    for (i, &(a, bb)) in vp.iter().enumerate() {
        for (j, &(c, d)) in vp.iter().enumerate() {
            let want = at(a, bb, c, d);
            assert!((c_model[i][j] - want).abs() < 1e-6 * 500e9, "{i}{j}: {} against {want}", c_model[i][j]);
        }
    }
}

#[test]
fn a_loaded_cubic_plate_in_the_x_z_plane_deflects_like_the_same_plate_in_the_model_frame() {
    // wafer description: fixed at x = 0 (x and z components), a load along -z on the right end; the crystal is turned by 20 degrees about y
    let mesh = SimplexMesh::rectangle(1.0, 0.2, 10, 2);
    let rot_y = |t: f64| -> M3 {
        let (s, c) = t.to_radians().sin_cos();
        [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]
    };
    let eul = euler_from_rotation(&rot_y(20.0));
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.plane = Plane::Stress;
    spec.materials.set(1, "E", 1e9);
    spec.materials.set(1, "nu", 0.3);
    spec.materials.aniso.insert(1, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }), euler_deg: eul, ..Default::default() });
    let mut w = spec.clone();
    w.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } }, Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e6]) }];
    let rw = femcore::solve::solve(&mesh, &w, &mut SolveState::default()).unwrap();
    // the model-frame description: the same rotation expressed in (x, z)
    let mut m = spec;
    m.materials.aniso.get_mut(&1).unwrap().euler_deg = femcore::frame_adapter::euler_to_model_frame(eul);
    m.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, true, false], values: [0.0; 3] } }, Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, -1e6, 0.0]) }];
    let rm = femcore::solid::solve_solid(&mesh, &m, &SolveState::default()).unwrap();
    let (uz_w, uy_m) = (&rw.solutions[0].field("uz").unwrap().values, &rm.solutions[0].field("uy").unwrap().values);
    let scale = uy_m.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    assert!(scale > 0.0 && maxdiff(uz_w, uy_m) < 1e-9 * scale, "{} of {scale}", maxdiff(uz_w, uy_m));
    // turning the crystal about the normal z instead changes the answer (a different rotation)
    let mut other = m.clone();
    other.materials.aniso.get_mut(&1).unwrap().euler_deg = [20.0, 0.0, 0.0];
    other.bcs = w.bcs.clone();
    let ro = femcore::solve::solve(&mesh, &other, &mut SolveState::default()).unwrap();
    let uz_o = &ro.solutions[0].field("uz").unwrap().values;
    assert!(maxdiff(uz_w, uz_o) > 1e-6 * scale);
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
fn layers_of_a_2d_cross_section_are_written_along_x_and_z_and_a_positive_k_runs_along_plus_y() {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Minus, strength: 3.0, order: 2.0, re: 1.0 });
    s.out_of_plane = Some(5.0);
    let m = femcore::frame_adapter::to_model_frame(&s, 2);
    // the solvers' second axis is z, and their out-of-plane axis is -y
    let p = &m.materials.pml[&2];
    assert_eq!((p.x, p.y, p.z), (PmlDir::Plus, PmlDir::Minus, PmlDir::Off));
    assert_eq!(m.out_of_plane, Some(-5.0));
    // a y layer does not exist on a 2D cross-section
    let mesh = SimplexMesh::rectangle(1.0, 0.5, 4, 2);
    let mut bad = ProblemSpec::new(Physics::Solid);
    bad.study = Study::Eigenfrequency { n_modes: 2, f_search: 0.0 };
    bad.materials.pml.insert(1, Pml { y: PmlDir::Plus, ..Default::default() });
    let e = femcore::solve::solve(&mesh, &bad, &mut SolveState::default()).err().unwrap();
    assert!(e.contains("absorbs along y") && e.contains("use z"), "{e}");
    // the Brillouin optical layers follow the same rule
    let mut b = ProblemSpec::new(Physics::Brillouin);
    b.brillouin.as_mut().unwrap().optical_pml.insert(3, Pml { x: PmlDir::Off, y: PmlDir::Off, z: PmlDir::Plus, ..Default::default() });
    let m = femcore::frame_adapter::to_model_frame(&b, 2);
    let p = &m.brillouin.as_ref().unwrap().optical_pml[&3];
    assert_eq!((p.y, p.z), (PmlDir::Plus, PmlDir::Off));
    // a 3D problem is untouched
    let m3 = femcore::frame_adapter::to_model_frame(&s, 3);
    assert_eq!(m3.out_of_plane, Some(5.0));
    assert_eq!(m3.materials.pml[&2].z, PmlDir::Minus);
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

#[test]
fn a_problem_saved_in_the_legacy_model_axes_is_refused_and_the_wafer_key_is_not_written() {
    let mut s = ProblemSpec::new(Physics::Electrostatics);
    s.materials.default.insert("eps_r".into(), 1.0);
    s.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    let json = serde_json::to_value(&s).unwrap();
    assert!(json["materials"].get("frame").is_none(), "there is no frame to save any more");
    // files of the brief window with an Axes toggle say what they are
    let mut old = json.clone();
    old["materials"]["frame"] = serde_json::json!("Model");
    let legacy: ProblemSpec = serde_json::from_value(old).unwrap();
    let mesh = SimplexMesh::rectangle(1.0, 0.5, 4, 2);
    let e = femcore::solve::solve(&mesh, &legacy, &mut SolveState::default()).err().unwrap();
    assert!(e.contains("legacy model axes"), "{e}");
    // and the old default, "Wafer", is simply what every problem is
    let mut wafer = json;
    wafer["materials"]["frame"] = serde_json::json!("Wafer");
    let ok: ProblemSpec = serde_json::from_value(wafer).unwrap();
    assert!(femcore::solve::solve(&mesh, &ok, &mut SolveState::default()).is_ok());
}
