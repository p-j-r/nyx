//! The wafer frame: x-y plane of the wafer, z the surface normal, a 2D cross-section is the x-z plane (the solid mechanics part; the
//! electromagnetic and geometry parts stay in femcore).
mod common;
use common::try_solve;
use femcore::material::*;
use femcore::props::*;
use nyx_sms::*;

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
    let mut spec = Problem::new(Study::Static, 2);
    spec.plane = Plane::Stress;
    spec.materials.set(1, "E", 1e9);
    spec.materials.set(1, "nu", 0.3);
    spec.materials.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Cubic { c11: 165.7e9, c12: 63.9e9, c44: 79.6e9 }),
            euler_deg: euler_from_rotation(&rot_y(20.0)),
            ..Default::default()
        },
    );
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] } },
        Bc { labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e6]) },
    ];
    let tip = |spec: &Problem| -> f64 {
        let r = common::solve(&mesh, spec);
        r.solutions[0].flat("uz").unwrap().iter().fold(0.0f64, |a, b| a.max(b.abs()))
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
fn a_layer_along_y_does_not_exist_on_a_2d_cross_section() {
    // layers of a 2D cross-section are written along x and z; there is no y axis in the plane
    let mesh = common::rect(1.0, 0.5, 4, 2);
    let mut bad = Problem::new(Study::Eigenfrequency { n_modes: 2, f_search: 0.0 }, 2);
    bad.materials.pml.insert(1, Pml { y: PmlDir::Plus, ..Default::default() });
    let e = try_solve(&mesh, &bad).err().unwrap();
    assert!(e.contains("absorbs along y") && e.contains("use z"), "{e}");
}
