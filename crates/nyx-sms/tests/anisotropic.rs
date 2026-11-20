// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Anisotropic materials: crystal stiffness, tensor conductivity/permittivity/expansion.
mod common;
use common::*;
use femcore::material::*;
use nyx_sms::*;

fn fixed(f: [bool; 3]) -> BcKind {
    BcKind::Fixed { fixed: f, values: [0.0; 3] }
}

fn solid_spec(plane: Plane, study: Study) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.plane = plane;
    s.study = study;
    s.materials.default.insert("rho".into(), 2329.0);
    s
}

#[test]
fn cubic_with_isotropic_constants_reproduces_the_isotropic_solution() {
    let mesh = common::block(1.0, 0.2, 0.2, 6, 2, 2);
    let (e, nu) = (70e9, 0.33);
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let bcs = vec![
        Bc { labels: vec![1], kind: fixed([true; 3]) },
        Bc { labels: vec![6], kind: BcKind::Load([1e6, 2e5, -3e5]) },
    ];
    let mut iso = solid_spec(Plane::Strain, Study::Static);
    iso.materials.default.insert("E".into(), e);
    iso.materials.default.insert("nu".into(), nu);
    iso.bcs = bcs.clone();
    let mut cub = solid_spec(Plane::Strain, Study::Static);
    cub.bcs = bcs;
    cub.materials.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Cubic { c11: lam + 2.0 * mu, c12: lam, c44: mu }),
            euler_deg: [20.0, 40.0, 60.0],
            ..Default::default()
        },
    );
    let a = solve_solid(&mesh, &iso).unwrap();
    let b = solve_solid(&mesh, &cub).unwrap();
    for name in ["ux", "uy", "uz", "von Mises", "sxy", "syz"] {
        let (x, y) = (&a.solutions[0].flat(name).unwrap(), &b.solutions[0].flat(name).unwrap());
        let sc = x.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(maxdiff(x, y) < 1e-9 * sc, "{name}: {}", maxdiff(x, y));
    }
}

#[test]
fn cubic_plate_stiffness_depends_on_orientation_like_the_compliance_formula() {
    // plane stress bar under uniaxial traction along x; the cubic crystal is rotated about z.
    // 1/E_x(theta) = S11 - 2 (S11 - S12 - S44/2) sin^2 cos^2, from the cubic compliances
    let (c11, c12, c44) = (165.7e9, 63.9e9, 79.6e9);
    let den = (c11 - c12) * (c11 + 2.0 * c12);
    let (s11, s12, s44) = ((c11 + c12) / den, -c12 / den, 1.0 / c44);
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let sigma = 1e6;
    for theta in [0.0f64, 45.0] {
        let mut spec = solid_spec(Plane::Stress, Study::Static);
        spec.materials.aniso.insert(
            1,
            Anisotropy {
                stiffness: Some(Stiffness::Cubic { c11, c12, c44 }),
                euler_deg: [theta, 0.0, 0.0],
                ..Default::default()
            },
        );
        spec.bcs = vec![
            Bc { labels: vec![1], kind: fixed([true, false, false]) },
            Bc { labels: vec![2], kind: fixed([false, true, false]) },
            Bc { labels: vec![4], kind: BcKind::Load([sigma, 0.0, 0.0]) },
        ];
        let r = solve_solid(&mesh, &spec).unwrap();
        let ux = &r.solutions[0].flat("ux").unwrap();
        let x = mesh.dof_coordinates(2);
        let t = theta.to_radians();
        let inv_e = s11 - 2.0 * (s11 - s12 - 0.5 * s44) * (t.sin() * t.cos()).powi(2);
        let exact: Vec<f64> = (0..ux.len()).map(|i| sigma * inv_e * x[2 * i]).collect();
        assert!(
            maxdiff(ux, &exact) < 1e-9 * exact.iter().cloned().fold(0.0, f64::max),
            "theta {theta}: {}",
            maxdiff(ux, &exact)
        );
    }
    // the two orientations really differ (Si is about 30 % stiffer along <110>)
    let e0 = 1.0 / s11;
    let e45 = 1.0 / (s11 - 0.5 * (s11 - s12 - 0.5 * s44));
    assert!(e45 / e0 > 1.2);
}

#[test]
fn sapphire_modes_follow_the_crystal_symmetry() {
    let mesh = common::block(1.0, 0.3, 0.3, 8, 3, 3);
    let (st, ..) = sapphire();
    let freqs = |euler: [f64; 3]| {
        let mut spec = solid_spec(Plane::Strain, Study::Eigenfrequency { n_modes: 5, f_search: 0.0 });
        spec.materials.default.insert("rho".into(), 3980.0);
        spec.materials
            .aniso
            .insert(1, Anisotropy { stiffness: Some(st.clone()), euler_deg: euler, ..Default::default() });
        spec.bcs = vec![Bc { labels: vec![1], kind: fixed([true; 3]) }];
        let r = solve_solid(&mesh, &spec).unwrap();
        r.solutions.iter().map(|s| s.value.unwrap()).collect::<Vec<f64>>()
    };
    let f0 = freqs([0.0, 0.0, 0.0]);
    let f120 = freqs([120.0, 0.0, 0.0]);
    for (a, b) in f0.iter().zip(&f120) {
        assert!((a / b - 1.0).abs() < 1e-9, "{f0:?} vs {f120:?}");
    }
    // the c axis turned along the bar (x): the Young.s modulus along the bar changes, so the bending frequencies change
    let f90 = freqs([90.0, 90.0, 0.0]);
    assert!(f0.iter().zip(&f90).any(|(a, b)| (a / b - 1.0).abs() > 1e-3), "{f0:?} vs {f90:?}");
}

#[test]
fn anisotropic_thermal_expansion_of_a_free_plate_gives_no_stress() {
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let (ax, ay) = (8e-6, 2e-6);
    let dt = 100.0;
    for plane in [Plane::Stress, Plane::Strain] {
        let mut spec = solid_spec(plane, Study::Static);
        spec.materials.default.insert("E".into(), 100e9);
        spec.materials.default.insert("nu".into(), 0.25);
        spec.materials
            .aniso
            .insert(1, Anisotropy { expansion: Some(Tensor::Diagonal { x: ax, y: ay, z: ay }), ..Default::default() });
        spec.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(293.15 + dt), t_ref: 293.15 });
        spec.bcs = vec![
            Bc { labels: vec![1], kind: fixed([true, false, false]) },
            Bc { labels: vec![2], kind: fixed([false, false, true]) },
        ];
        let r = solve_solid(&mesh, &spec).unwrap();
        let s = &r.solutions[0];
        let x = mesh.dof_coordinates(2);
        let ux = &s.flat("ux").unwrap();
        let uz = &s.flat("uz").unwrap();
        // plane stress: free expansion alpha dT; plane strain: eps_yy = 0 (y is out of the section) adds nu alpha_y dT in plane
        let (az, nu) = (ay, 0.25); // az: the expansion along the axis out of the section (y)
        let (kx, ky) = if plane == Plane::Stress { (ax, ay) } else { (ax + nu * az, ay + nu * az) };
        let ex: Vec<f64> = (0..ux.len()).map(|i| kx * dt * x[2 * i]).collect();
        let ey: Vec<f64> = (0..ux.len()).map(|i| ky * dt * x[2 * i + 1]).collect();
        assert!(maxdiff(ux, &ex) < 1e-9 * kx * dt, "{plane:?} ux {}", maxdiff(ux, &ex));
        assert!(maxdiff(uz, &ey) < 1e-9 * ky * dt, "{plane:?} uz {}", maxdiff(uz, &ey));
        for name in ["sxx", "szz", "sxz"] {
            let v = s.flat(name).unwrap().iter().fold(0.0f64, |m, x| m.max(x.abs()));
            assert!(v < 1e-3, "{plane:?} {name} should vanish, got {v} Pa");
        }
        let syy = s.flat("syy").unwrap().iter().cloned().fold(0.0f64, |m, x| if x.abs() > m.abs() { x } else { m });
        let want = if plane == Plane::Stress { 0.0 } else { -100e9 * az * dt };
        assert!((syy - want).abs() < 1e-6 * want.abs().max(1.0), "{plane:?} syy {syy} vs {want}");
    }
}

#[test]
fn isotropic_plane_strain_free_expansion_is_stress_free() {
    // the physical check of the thermal stress recovery: uniform heating of an unconstrained plate
    let mesh = common::rect(1.0, 0.2, 10, 3);
    let mut spec = solid_spec(Plane::Strain, Study::Static);
    for (k, v) in [("E", 70e9), ("nu", 0.33), ("alpha", 23e-6)] {
        spec.materials.default.insert(k.into(), v);
    }
    spec.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
    spec.bcs = vec![
        Bc { labels: vec![1], kind: fixed([true, false, false]) },
        Bc { labels: vec![2], kind: fixed([false, false, true]) },
    ];
    let r = solve_solid(&mesh, &spec).unwrap();
    // in-plane stress vanishes; the constrained out-of-plane direction (y) carries -E alpha dT
    for name in ["sxx", "szz", "sxz"] {
        let v = r.solutions[0].flat(name).unwrap().iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(v < 1e-3, "{name}: {v} Pa");
    }
    let syy = r.solutions[0].flat("syy").unwrap()[0];
    assert!((syy + 70e9 * 23e-6 * 100.0).abs() < 1e-3, "{syy}");
}

#[test]
fn nonpositive_definite_stiffness_is_rejected() {
    let mesh = common::rect(1.0, 0.2, 4, 2);
    let mut spec = solid_spec(Plane::Strain, Study::Static);
    spec.materials.aniso.insert(
        1,
        Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 10e9, c12: 60e9, c44: 5e9 }), ..Default::default() },
    );
    spec.bcs = vec![Bc { labels: vec![1], kind: fixed([true; 3]) }];
    assert!(solve_solid(&mesh, &spec).unwrap_err().contains("positive definite"));
}
