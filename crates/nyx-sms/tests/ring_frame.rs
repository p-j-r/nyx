// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The ring (axisymmetric) kernels in the wafer frame.
//!
//! A ring resonator is solved on the (r, z) half-plane: mesh x = r, mesh z = z, and the azimuthal direction theta is the third
//! displacement/field component. In the wafer frame the ring's local axes are (x, y, z) = (r, theta, z): x outward, y along the ring,
//! z the surface normal (right-handed), so the kernels read a wafer-frame problem as it is, with the components (u_r, u_z, u_theta) =
//! (ux, uz, uy). These tests pin down what the ring needs:
//!
//! 1. A crystal ring is averaged over rotations about the ring axis, which must be the wafer normal z. So the elastic modes cannot
//!    depend on the azimuth of the crystal: a c axis turned about z by any angle (sapphire is trigonal, so only a lab rotation about z is a symmetry, not
//!    "tilting towards x instead of y", which also turns the crystal about its own c axis) gives the same
//!    frequencies, while a different tilt angle does not (the test can fail).
//! 2. The frequencies of a tilted crystal ring equal the numbers frozen from the kernels when the model axes (r, z, theta) were still a
//!    public frame.
//! 3. A torsional mode (m = 0, only u_theta) is called `uy` in the wafer frame, never `ux` or `uz`.
//! 4. The direction of increasing theta: an eigenvector has an arbitrary sign, so a mode cannot show it, but a rigid translation can. Along a
//!    fixed Cartesian X it is u_r = cos(theta), u_theta = -sin(theta), whose exp(-i theta) part is (1, -i): the stiffness of the m = 1 kernel
//!    must annihilate (u_r, u_z, u_theta) = (1, 0, -i) and not (1, 0, +i). So the azimuthal dependence is exp(-i m theta) with theta along wafer +y.
//! 5. The optical ring modes: with the time convention e^{j omega t} of the rest of the library, exp(-i m theta) with m > 0 carries power
//!    along +theta, i.e. +y (and m < 0 along -y).
//! 6. The crystal data of a ring are averaged about the ring axis, which makes the sign convention of the third axis irrelevant for
//!    orientations.
mod common;
use femcore::material::*;
use nyx_sms::*;

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
    common::parts_mesh(&[common::part(rr - a / 2.0, -a / 2.0, rr + a / 2.0, a / 2.0, 1, Some(a / 3.0))], a / 3.0)
}

/// Sapphire ring with its crystal axes given by the rotation `r`, in the frame `frame`.
fn sapphire_ring(r: &R3, m: f64) -> Problem {
    let (stiffness, _, _, _, rho) = sapphire();
    let mut sp = Problem::new(Study::Static, 2);
    sp.order = 2;
    sp.materials.default.insert("rho".into(), rho);
    sp.materials
        .aniso
        .insert(1, Anisotropy { stiffness: Some(stiffness), euler_deg: euler_from_rotation(r), ..Default::default() });
    sp.study = Study::Eigenfrequency { n_modes: 6, f_search: 40.0 * (1.0 + m) };
    sp
}

/// Lowest elastic frequencies of the ring at azimuthal order `m` of a wafer-frame problem.
fn ring_frequencies(spec: &Problem, m: f64) -> Vec<f64> {
    let (_, modes) = nyx_sms::solve_solid_axisym(&ring_mesh(), spec, m).unwrap();
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
    let reference = ring_frequencies(&sapphire_ring(&tilt_about_x(b), m), m);
    assert!(reference.len() >= 4, "{reference:?}");
    // a rotation of the whole crystal about the ring axis z (a rotation applied in the lab frame, on the left) cannot matter
    for deg in [37.0, 90.0, -90.0, 180.0] {
        let f = ring_frequencies(&sapphire_ring(&mul(&turn_about_z(deg), &tilt_about_x(b)), m), m);
        assert_same(&f, &reference, &format!("crystal turned {deg} degrees about z"));
    }
    // the check can fail: a different angle between the c axis and the normal is a different crystal ring
    let other = ring_frequencies(&sapphire_ring(&tilt_about_x(20.0), m), m);
    let worst = other.iter().zip(&reference).map(|(a, b)| (a - b).abs() / b.abs()).fold(0.0, f64::max);
    assert!(worst > 1e-3, "a 20 degree tilt gave the same modes as 57.6 degrees: {other:?} against {reference:?}");
}

#[test]
fn a_ring_in_the_wafer_frame_has_the_frozen_frequencies_of_the_model_axes() {
    let m = 2.0;
    let b = 57.6;
    // wafer axes: c axis tilted from the normal z towards x (the radial direction)
    let wafer = ring_frequencies(&sapphire_ring(&tilt_about_y(b), m), m);
    // frozen when the same ring was described in the model axes of the (r, z) half-plane (NYX_WRITE_GOLDEN=1 rewrites the file)
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/ring_model_axes.json");
    if std::env::var("NYX_WRITE_GOLDEN").is_ok() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let text = format!("{{\"description\": \"sapphire ring, c axis 57.6 deg from z towards r, model axes (r, z, theta), azimuthal order 2: lowest elastic frequencies (Hz)\", \"m\": 2, \"frequencies\": {:?}}}\n", wafer);
        std::fs::write(&path, text).unwrap();
    }
    let golden: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&path).expect("tests/golden/ring_model_axes.json")).unwrap();
    let frozen: Vec<f64> = golden["frequencies"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    assert_same(&wafer, &frozen, "frozen frequencies");
}

#[test]
fn a_torsional_ring_mode_is_the_wafer_y_component() {
    // an isotropic ring at m = 0: besides the breathing and the in-plane/out-of-plane modes there are torsional modes that move only along theta
    let mut sp = Problem::new(Study::Static, 2);
    sp.order = 2;
    sp.materials.default.insert("E".into(), 70e9);
    sp.materials.default.insert("nu".into(), 0.33);
    sp.materials.default.insert("rho".into(), 2700.0);
    sp.study = Study::Eigenfrequency { n_modes: 8, f_search: 500.0 };
    let (res, _) = nyx_sms::solve_solid_axisym(&ring_mesh(), &sp, 0.0).unwrap();
    let energy = |s: &nyx_physics::Solution, name: &str| {
        s.field(name).map_or(0.0, |f| f.values.iter().map(|v| v * v).sum::<f64>())
    };
    // the components (u_r, u_z, u_theta) are named Re ux, Re uz, Re uy
    let torsional: Vec<usize> = (0..res.solutions.len())
        .filter(|&i| {
            let s = &res.solutions[i];
            let (x, y, z) = (energy(s, "Re ux"), energy(s, "Re uy"), energy(s, "Re uz"));
            y > 0.0 && y > 100.0 * (x + z)
        })
        .collect();
    assert!(!torsional.is_empty(), "no torsional mode among the {} modes found", res.solutions.len());
}

#[test]
fn the_ring_kernel_has_exp_minus_i_m_theta_with_theta_along_wafer_plus_y() {
    use femcore::linalg::csr::c;
    let mesh = ring_mesh();
    let c6 = vec![isotropic_stiffness(70e9, 0.33); mesh.n_cells()];
    let k = nyx_sms::axisym::assemble_stiffness_axisym(&mesh, &c6, 2, 1.0, None).unwrap();
    let n1 = mesh.n_dofs(2);
    // components (u_r, u_z, u_theta) as blocks of n1 dofs; uniform fields
    let uniform = |s: f64| {
        let mut u = vec![c(0.0, 0.0); 3 * n1];
        for i in 0..n1 {
            u[i] = c(1.0, 0.0);
            u[2 * n1 + i] = c(0.0, s);
        }
        u
    };
    let size = |v: Vec<femcore::linalg::csr::C>| v.iter().map(|z| z.norm()).fold(0.0, f64::max);
    let translation = size(k.matvec(&uniform(-1.0)));
    let wrong = size(k.matvec(&uniform(1.0)));
    assert!(wrong > 1e3, "the check cannot fail: |K u(+i)| = {wrong}");
    assert!(translation < 1e-6 * wrong, "a rigid translation along X costs {translation} against {wrong}");
}

/// Fourth-rank transform of a Voigt matrix (xx yy zz xy yz xz) by the 3x3 matrix `s` (a rotation or a reflection).
fn transform_c6(c6: &M6, s: &R3) -> M6 {
    let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            for (p, q) in [(a, b), (b, a)] {
                for (u, v) in [(cc, d), (d, cc)] {
                    t[p][q][u][v] = c6[i][j];
                }
            }
        }
    }
    let mut out = [[0.0; 6]; 6];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            let mut acc = 0.0;
            for p in 0..3 {
                for q in 0..3 {
                    for u in 0..3 {
                        for v in 0..3 {
                            acc += s[a][p] * s[b][q] * s[cc][u] * s[d][v] * t[p][q][u][v];
                        }
                    }
                }
            }
            out[i][j] = acc;
        }
    }
    out
}

#[test]
fn averaging_about_the_ring_axis_makes_the_sign_of_the_third_axis_irrelevant() {
    // sapphire with its c axis tilted from z towards the ring direction y, in the wafer frame
    let (stiffness, _, _, _, _) = sapphire();
    let wafer = stiffness_global(&stiffness, euler_from_rotation(&tilt_about_x(57.6)));
    // kernel axes (r, z, theta) = wafer (x, z, +y) for a ring, (x, z, -y) for a straight section: reflections / rotations as 3x3 maps
    let ring: R3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]];
    let section: R3 = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]];
    let (cr, cs) = (transform_c6(&wafer, &ring), transform_c6(&wafer, &section));
    let diff = |a: &M6, b: &M6| (0..6).flat_map(|i| (0..6).map(move |j| (a[i][j] - b[i][j]).abs())).fold(0.0, f64::max);
    // before the averaging the two conventions are different crystals (the check can fail) ...
    assert!(diff(&cr, &cs) > 1e9, "{}", diff(&cr, &cs));
    // ... after it, the same ring
    let (ar, a_s) = (nyx_sms::axisym::average_about_axis(&cr, 1, 16), nyx_sms::axisym::average_about_axis(&cs, 1, 16));
    assert!(diff(&ar, &a_s) < 1e-6 * 500e9, "{}", diff(&ar, &a_s));
}
