//! The ring (axisymmetric) path in the wafer frame.
//!
//! A ring resonator is solved on the (r, z) half-plane: mesh x = r, mesh y = z, and the azimuthal direction theta is the third
//! displacement/field component. In the wafer frame the ring's local axes are (x, y, z) = (r, theta, z): x outward, y along the ring,
//! z the surface normal. `solve` turns a 2D wafer-frame problem into the solvers' model frame (`ProblemSpec::to_model_frame`) and
//! renames the result fields (`rename_wafer_fields`); both were written for the straight-waveguide kernels. These tests pin down what
//! the ring path needs from them, without any reference value: only invariances and a naming check.
//!
//! 1. A crystal ring is averaged over rotations about the ring axis, which must be the wafer normal z. So the elastic modes cannot
//!    depend on the azimuth of the crystal: a c axis turned about z by any angle (sapphire is trigonal, so only a lab rotation about z is a symmetry, not
//!    "tilting towards x instead of y", which also turns the crystal about its own c axis) gives the same
//!    frequencies, while a different tilt angle does not (the test can fail).
//! 2. The same crystal described in the legacy model axes (2D mesh = x-y plane, the ring axis is y) gives the same frequencies as its
//!    wafer description.
//! 3. A torsional mode (m = 0, only u_theta) is called `uy` in the wafer frame, never `ux` or `uz`.
//!
//! Not covered (open, needs a decision): the sign of the u_theta component. The kernel's third component is +theta (azimuth, exp(-i m theta)),
//! while `rename_wafer_fields` maps it to `-uy` (the model frame's out-of-plane axis is -y). With y = +theta (right-handed wafer axes) that
//! sign is wrong for rings; eigenvector signs are arbitrary, so an eigenmode cannot show it.
use femcore::material::*;
use femcore::mesher::*;
use femcore::problem::*;

type R3 = [[f64; 3]; 3];

fn mul(a: &R3, b: &R3) -> R3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}

/// Columns of the result are the crystal axes in the lab frame; the c axis is the third column.
fn tilt_about_x(deg: f64) -> R3 {
    let (s, c) = deg.to_radians().sin_cos();
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]] // c axis (0, s, c): tilted from z towards +y
}
fn tilt_about_y(deg: f64) -> R3 {
    let (s, c) = deg.to_radians().sin_cos();
    [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]] // c axis (s, 0, c): tilted from z towards +x
}
fn turn_about_z(deg: f64) -> R3 {
    let (s, c) = deg.to_radians().sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// A ring of mean radius 1 m with a 0.02 m x 0.02 m cross-section (mesh x = r, y = z).
fn ring_mesh() -> femcore::mesh::SimplexMesh {
    let (rr, a) = (1.0, 0.02);
    let parts = vec![Part { shape: Shape::Rect { x0: rr - a / 2.0, y0: -a / 2.0, x1: rr + a / 2.0, y1: a / 2.0 }, domain: 1, h: Some(a / 3.0) }];
    mesh_2d(&Spec2d { parts, h: a / 3.0, min_angle_deg: 25.0 }).unwrap()
}

/// Sapphire ring with its crystal axes given by the rotation `r`, in the frame `frame`.
fn sapphire_ring(r: &R3, frame: Frame, m: f64) -> ProblemSpec {
    let (stiffness, _, _, _, rho) = sapphire();
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = 2;
    sp.materials.frame = frame;
    sp.materials.default.insert("rho".into(), rho);
    sp.materials.aniso.insert(1, Anisotropy { stiffness: Some(stiffness), euler_deg: euler_from_rotation(r), ..Default::default() });
    sp.study = Study::Eigenfrequency { n_modes: 6, f_search: 40.0 * (1.0 + m) };
    sp
}

/// Lowest elastic frequencies of the ring at azimuthal order `m`, the way `solve` hands a wafer-frame problem to the ring kernels.
fn ring_frequencies(spec: &ProblemSpec, m: f64) -> Vec<f64> {
    let inner = spec.to_model_frame(2);
    let (_, modes) = femcore::solid_kz::solve_solid_axisym(&ring_mesh(), &inner, m).unwrap();
    let mut f: Vec<f64> = modes.freqs.iter().map(|z| z.re).collect();
    f.sort_by(|a, b| a.total_cmp(b));
    f
}

fn assert_same(a: &[f64], b: &[f64], what: &str) {
    assert_eq!(a.len(), b.len(), "{what}: {a:?} against {b:?}");
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-6 * y.abs().max(1.0), "{what}: {a:?} against {b:?}");
    }
}

#[test]
fn the_elastic_modes_of_a_crystal_ring_do_not_depend_on_the_azimuth_of_its_c_axis() {
    let m = 2.0;
    let b = 57.6;
    let reference = ring_frequencies(&sapphire_ring(&tilt_about_x(b), Frame::Wafer, m), m);
    assert!(reference.len() >= 4, "{reference:?}");
    // a rotation of the whole crystal about the ring axis z (a rotation applied in the lab frame, on the left) cannot matter
    for deg in [37.0, 90.0, -90.0, 180.0] {
        let f = ring_frequencies(&sapphire_ring(&mul(&turn_about_z(deg), &tilt_about_x(b)), Frame::Wafer, m), m);
        assert_same(&f, &reference, &format!("crystal turned {deg} degrees about z"));
    }
    // the check can fail: a different angle between the c axis and the normal is a different crystal ring
    let other = ring_frequencies(&sapphire_ring(&tilt_about_x(20.0), Frame::Wafer, m), m);
    let worst = other.iter().zip(&reference).map(|(a, b)| (a - b).abs() / b.abs()).fold(0.0, f64::max);
    assert!(worst > 1e-3, "a 20 degree tilt gave the same modes as 57.6 degrees: {other:?} against {reference:?}");
}

#[test]
fn a_ring_described_in_the_wafer_frame_and_in_the_model_axes_has_the_same_modes() {
    let m = 2.0;
    let b = 57.6;
    // wafer axes: c axis tilted from the normal z towards x (the radial direction)
    let wafer = ring_frequencies(&sapphire_ring(&tilt_about_y(b), Frame::Wafer, m), m);
    // model axes of the (r, z) half-plane: x = r, y = z (the ring axis), the c axis (sin b, cos b, 0) is tilted from y towards x
    let to_model_y: R3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]; // sends the crystal c axis to model y
    let (s, c) = b.to_radians().sin_cos();
    let about_model_z: R3 = [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]; // turns model y towards model x
    let r_model = mul(&about_model_z, &to_model_y);
    let c_axis = [r_model[0][2], r_model[1][2], r_model[2][2]];
    assert!((c_axis[0] - s).abs() < 1e-12 && (c_axis[1] - c).abs() < 1e-12 && c_axis[2].abs() < 1e-12, "{c_axis:?}");
    let model = ring_frequencies(&sapphire_ring(&r_model, Frame::Model, m), m);
    assert_same(&model, &wafer, "model-axes description");
}

#[test]
fn a_torsional_ring_mode_is_the_wafer_y_component() {
    // an isotropic ring at m = 0: besides the breathing and the in-plane/out-of-plane modes there are torsional modes that move only along theta
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.order = 2;
    sp.materials.default.insert("E".into(), 70e9);
    sp.materials.default.insert("nu".into(), 0.33);
    sp.materials.default.insert("rho".into(), 2700.0);
    sp.study = Study::Eigenfrequency { n_modes: 8, f_search: 500.0 };
    let (mut res, _) = femcore::solid_kz::solve_solid_axisym(&ring_mesh(), &sp.to_model_frame(2), 0.0).unwrap();
    let energy = |s: &femcore::problem::Solution, name: &str| s.field(name).map_or(0.0, |f| f.values.iter().map(|v| v * v).sum::<f64>());
    // before the renaming the kernel's components are (u_r, u_z, u_theta) = (ux, uy, uz)
    let torsional: Vec<usize> = (0..res.solutions.len())
        .filter(|&i| {
            let s = &res.solutions[i];
            let (r, z, t) = (energy(s, "ux"), energy(s, "uy"), energy(s, "uz"));
            t > 0.0 && t > 100.0 * (r + z)
        })
        .collect();
    assert!(!torsional.is_empty(), "no torsional mode among the {} modes found", res.solutions.len());
    femcore::solve::rename_wafer_fields(&mut res);
    for &i in &torsional {
        let s = &res.solutions[i];
        let (x, y, z) = (energy(s, "ux"), energy(s, "uy"), energy(s, "uz"));
        assert!(y > 100.0 * (x + z), "torsion should be the wafer y component: ux {x:.3e}, uy {y:.3e}, uz {z:.3e}");
    }
}
