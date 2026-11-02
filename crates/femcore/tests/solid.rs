//! Solid mechanics against beam theory and closed forms.
mod common;
use common::*;
use femcore::problem::*;
use femcore::solid::solve_solid;

fn steel() -> Materials {
    let mut m = Materials::default();
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        m.default.insert(k.into(), v);
    }
    m
}

fn cantilever(plane: Plane, study: Study) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.materials = steel();
    s.plane = plane;
    s.study = study;
    s.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    s
}

fn planes() -> [(Plane, &'static str); 2] {
    [(Plane::Stress, "stress"), (Plane::Strain, "strain")]
}

/// Plane modulus: E for plane stress, E / (1 - nu^2) for plane strain.
fn plane_modulus(plane: Plane, e: f64, nu: f64) -> f64 {
    match plane {
        Plane::Stress => e,
        Plane::Strain => e / (1.0 - nu * nu),
    }
}

#[test]
fn cantilever_tip_deflection_matches_beam_theory() {
    // cantilever 1 x 0.1, E 200e9, nu 0.3, edge load -1e5 N/m^2 on boundary 4 (F = 1e4 N per metre of depth):
    // Euler-Bernoulli F L^3 / (3 E' I) plus the shear deflection F L / (kappa G A), kappa = 5/6
    let mesh = common::rect(1.0, 0.1, 40, 4);
    for (plane, name) in planes() {
        let mut spec = cantilever(plane, Study::Static);
        spec.bcs.push(Bc { on: None, labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e5]) });
        let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
        let uz = &r.solutions[0].field("uz").unwrap().values;
        let tip = uz.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let (f, l, h, e, nu) = (1e4f64, 1.0f64, 0.1f64, 200e9f64, 0.3f64);
        let ep = plane_modulus(plane, e, nu);
        let g = e / (2.0 * (1.0 + nu));
        let exact = f * l.powi(3) / (3.0 * ep * h.powi(3) / 12.0) + f * l / (5.0 / 6.0 * g * h);
        assert!((tip / exact - 1.0).abs() < 0.04, "{name}: tip {tip:e} against beam theory {exact:e}");
        // the largest von Mises stress is the bending stress at the clamp, 6 F L / h^2, up to the clamp's stress concentration
        let vm = r.solutions[0].field("von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
        let bending = 6.0 * f * l / (h * h);
        assert!(vm > 0.85 * bending && vm < 1.5 * bending, "{name}: von Mises {vm:e} against {bending:e}");
    }
}

#[test]
fn cantilever_eigenfrequencies_match_beam_theory() {
    // Euler-Bernoulli: f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E' I / (rho A)), beta_n L = 1.8751, 4.6941
    let mesh = common::rect(1.0, 0.1, 40, 4);
    for (plane, name) in planes() {
        let r = solve_solid(&mesh, &cantilever(plane, Study::Eigenfrequency { n_modes: 6, f_search: 0.0 }), &SolveState::default()).unwrap();
        let mine: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
        let (l, h, e, nu, rho) = (1.0f64, 0.1f64, 200e9, 0.3, 7850.0);
        let ep = plane_modulus(plane, e, nu);
        let unit = (ep * h * h / 12.0 / rho).sqrt() / (2.0 * std::f64::consts::PI * l * l);
        for (i, bl) in [1.8751f64, 4.6941].iter().enumerate() {
            let exact = bl * bl * unit;
            // shear deformation and rotary inertia lower the second mode of this stocky beam (h / L = 0.1) by about 5 %
            let tol = if i == 0 { 0.03 } else { 0.08 };
            assert!((mine[i] / exact - 1.0).abs() < tol, "{name} mode {i}: {} against {exact}", mine[i]);
        }
    }
}

#[test]
fn stress_at_triangle_centroids_is_uniaxial() {
    // roller supports (ux = 0 on the left, uy = 0 at the bottom) and 1 MPa on the right: sigma_xx = 1 MPa in every
    // triangle, plane strain adds sigma_zz = nu sigma_xx
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let nu = 0.3;
    for (plane, name) in planes() {
        let mut spec = ProblemSpec::new(Physics::Solid);
        spec.plane = plane;
        spec.materials.default.insert("E".into(), 1e9);
        spec.materials.default.insert("nu".into(), nu);
        spec.bcs = vec![
            Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
            Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
            Bc { on: None, labels: vec![4], kind: BcKind::Load([1e6, 0.0, 0.0]) },
        ];
        let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
        let s = &r.solutions[0];
        let n1 = mesh.n_dofs(2);
        let mut u = vec![0.0; 2 * n1];
        u[..n1].copy_from_slice(&s.field("ux").unwrap().values);
        u[n1..].copy_from_slice(&s.field("uz").unwrap().values);
        let nc = mesh.n_cells();
        let c6 = vec![femcore::material::isotropic_stiffness(1e9, nu); nc];
        let cs = femcore::solid::centroid_stress(&mesh, &u, &c6, plane, 2);
        // our row: [sxx, sxz, szx, szz, syy (out of the plane), vm]
        let szz = if plane == Plane::Strain { nu * 1e6 } else { 0.0 };
        let vm = (0.5 * ((1e6f64 - 0.0).powi(2) + (0.0 - szz).powi(2) + (szz - 1e6).powi(2))).sqrt();
        for (c, row) in cs.iter().enumerate() {
            for (k, want) in [1e6, 0.0, 0.0, 0.0, szz, vm].iter().enumerate() {
                assert!((row[k] - want).abs() < 1e-6 * 1e6, "{name} cell {c} column {k}: {} against {want}", row[k]);
            }
        }
    }
}

#[test]
fn thermal_expansion_of_a_free_bar() {
    // bar 1 x 0.1, E 70e9, nu 0.33, alpha 23e-6, uniform 393.15 K, reference 293.15 K, roller supports only (free to
    // expand): u = a (x, y) with a = alpha dT (plane stress) or (1 + nu) alpha dT (plane strain, sigma_zz = -E alpha dT)
    let mesh = common::rect(1.0, 0.1, 10, 2);
    let (e, nu, alpha, dt) = (70e9, 0.33, 23e-6, 100.0);
    for (plane, name) in planes() {
        let mut spec = cantilever(plane, Study::Static);
        for (k, v) in [("E", e), ("nu", nu), ("alpha", alpha)] {
            spec.materials.default.insert(k.into(), v);
        }
        spec.bcs = vec![
            Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
            Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        ];
        spec.thermal = Some(ThermalExpansion { source: ThermalSource::Uniform(393.15), t_ref: 293.15 });
        let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
        let s = &r.solutions[0];
        let a = if plane == Plane::Strain { (1.0 + nu) * alpha * dt } else { alpha * dt };
        let x = mesh.dof_coordinates(2);
        let (ux, uz) = (&s.field("ux").unwrap().values, &s.field("uz").unwrap().values);
        let ex: Vec<f64> = (0..ux.len()).map(|i| a * x[2 * i]).collect();
        let ey: Vec<f64> = (0..ux.len()).map(|i| a * x[2 * i + 1]).collect();
        assert!(maxdiff(ux, &ex) < 1e-9 * a, "{name} ux {}", maxdiff(ux, &ex));
        assert!(maxdiff(uz, &ey) < 1e-9 * a, "{name} uz {}", maxdiff(uz, &ey));
        let vm = s.field("von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
        let want = if plane == Plane::Strain { e * alpha * dt } else { 0.0 };
        assert!((vm - want).abs() < 1e-6 * e * alpha * dt, "{name} von Mises {vm} against {want}");
    }
}

#[test]
fn heat_result_drives_thermal_expansion() {
    let mesh = common::rect(1.0, 0.1, 10, 2);
    let mut heat = ProblemSpec::new(Physics::Heat);
    heat.materials.default.insert("k".into(), 1.0);
    heat.bcs = vec![Bc { on: None, labels: vec![1, 4], kind: BcKind::Temperature(393.15) }];
    let mut state = SolveState::default();
    femcore::scalar::solve_heat(&mesh, &heat, &mut state).unwrap();
    let mk = |source| {
        let mut spec = cantilever(Plane::Strain, Study::Static);
        for (k, v) in [("E", 70e9), ("nu", 0.33), ("alpha", 23e-6)] {
            spec.materials.default.insert(k.into(), v);
        }
        spec.thermal = Some(ThermalExpansion { source, t_ref: 293.15 });
        spec
    };
    let a = solve_solid(&mesh, &mk(ThermalSource::FromHeat), &state).unwrap();
    let b = solve_solid(&mesh, &mk(ThermalSource::Uniform(393.15)), &state).unwrap();
    let (ua, ub) = (&a.solutions[0].field("ux").unwrap().values, &b.solutions[0].field("ux").unwrap().values);
    assert!(maxdiff(ua, ub) < 1e-12);
    assert!(solve_solid(&mesh, &mk(ThermalSource::FromHeat), &SolveState::default()).is_err());
}
#[test]
fn block_3d_matches_beam_theory() {
    // block 1 x 0.1 x 0.1, fixed at x = 0 (boundary 1), load (0, 0, -1e5) N/m^2 on x = 1 (boundary 6): F = 1e3 N
    let mesh = common::block(1.0, 0.1, 0.1, 12, 2, 2);
    let mut spec = cantilever(Plane::Strain, Study::Static);
    spec.bcs.push(Bc { on: None, labels: vec![6], kind: BcKind::Load([0.0, 0.0, -1e5]) });
    let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
    let uz = &r.solutions[0].field("uz").unwrap().values;
    let tip = uz.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let (f, l, h, e, nu) = (1e3f64, 1.0f64, 0.1f64, 200e9f64, 0.3f64);
    let g = e / (2.0 * (1.0 + nu));
    let exact = f * l.powi(3) / (3.0 * e * h.powi(4) / 12.0) + f * l / (5.0 / 6.0 * g * h * h);
    assert!((tip / exact - 1.0).abs() < 0.08, "tip {tip:e} against beam theory {exact:e}");
    // two degenerate bending modes (y and z) at the Euler-Bernoulli frequency
    let spec = cantilever(Plane::Strain, Study::Eigenfrequency { n_modes: 6, f_search: 0.0 });
    let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
    let f1 = (e * h * h / 12.0 / 7850.0f64).sqrt() * 1.8751f64.powi(2) / (2.0 * std::f64::consts::PI * l * l);
    for sol in r.solutions.iter().take(2) {
        assert!((sol.value.unwrap() / f1 - 1.0).abs() < 0.06, "{} against {f1}", sol.value.unwrap());
    }
}

#[test]
fn roller_constraint_gives_uniaxial_strain_field() {
    // ux = 0 on the left face only, uy = 0 on the bottom only, nu = 0, load 1 MPa on the right:
    // uniaxial stress, ux = sigma x / E exactly (the P2 field is exactly linear)
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.plane = Plane::Stress;
    spec.materials.default.insert("E".into(), 1e9);
    spec.materials.default.insert("nu".into(), 0.0);
    spec.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { on: None, labels: vec![4], kind: BcKind::Load([1e6, 0.0, 0.0]) },
    ];
    let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
    let ux = &r.solutions[0].field("ux").unwrap().values;
    let x = mesh.dof_coordinates(2);
    let exact: Vec<f64> = (0..ux.len()).map(|i| 1e6 * x[2 * i] / 1e9).collect();
    assert!(maxdiff(ux, &exact) < 1e-12);
    let vm = r.solutions[0].field("von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
    assert!((vm / 1e6 - 1.0).abs() < 1e-9);
}

#[test]
fn prescribed_displacement_and_two_materials() {
    // two layers in series with a prescribed end displacement: stress equal in both layers
    let mesh = common::rect_split_x(1.0, 0.2, 0.5, 0.1);
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.plane = Plane::Stress;
    spec.materials.set(1, "E", 100e9);
    spec.materials.set(2, "E", 400e9);
    spec.materials.default.insert("nu".into(), 0.0);
    spec.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true, false, false], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } },
        Bc { on: None, labels: vec![4], kind: BcKind::Fixed { fixed: [true, false, false], values: [1e-3, 0.0, 0.0] } },
    ];
    let r = solve_solid(&mesh, &spec, &SolveState::default()).unwrap();
    // strain split: sigma = E1 e1 = E2 e2, e1 * 0.5 + e2 * 0.5 = 1e-3
    let e1 = 1e-3 / (0.5 + 0.5 * 100.0 / 400.0);
    let sigma = 100e9 * e1;
    let vm = r.solutions[0].field("von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
    assert!((vm / sigma - 1.0).abs() < 1e-6, "{vm} vs {sigma}");
}
