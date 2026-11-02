//! Elastic modes of a cross-section with an out-of-plane wave number.
mod common;
use femcore::problem::*;

fn spec(kz: f64, f_search: f64, n: usize) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Solid);
    s.materials.default.insert("E".into(), 70e9);
    s.materials.default.insert("nu".into(), 0.3);
    s.materials.default.insert("rho".into(), 2700.0);
    s.study = Study::Eigenfrequency { n_modes: n, f_search };
    s.out_of_plane = Some(kz);
    s
}

#[test]
fn clamped_plate_shear_modes_at_zero_wave_number() {
    // kz = 0 decouples u_z: clamped on all sides, f = (1/2) sqrt(mu/rho) sqrt(1/a^2 + 1/b^2)
    let (a, b) = (0.02, 0.01);
    let mesh = common::rect(a, b, 24, 12);
    let mut s = spec(0.0, 0.0, 4);
    s.bcs = vec![Bc { on: None, labels: vec![1, 2, 3, 4], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let mu = 70e9 / (2.0 * 1.3);
    let want = 0.5 * (mu / 2700.0f64).sqrt() * (1.0 / (a * a) + 1.0 / (b * b)).sqrt();
    let found = r.solutions.iter().map(|x| x.value.unwrap()).collect::<Vec<_>>();
    assert!(found.iter().any(|f| (f - want).abs() < 2e-3 * want), "expected {want}, found {found:?}");
}

#[test]
fn free_rod_extensional_and_flexural_branches() {
    // a free square rod at small kz a: extensional wave omega = sqrt(E/rho) kz (Young's rod), bending
    // omega = sqrt(E I / (rho A)) kz^2 with I/A = a^2 / 12 (Euler-Bernoulli)
    let a = 0.01;
    let mesh = common::rect(a, a, 16, 16);
    let kz = 20.0; // kz a = 0.2: corrections are of order a few percent for bending
    let c_ext = (70e9f64 / 2700.0).sqrt();
    let f_ext = c_ext * kz / (2.0 * std::f64::consts::PI);
    let f_flex = (a / 12f64.sqrt()) * c_ext * kz * kz / (2.0 * std::f64::consts::PI);
    let s = spec(kz, f_ext, 8);
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let freqs: Vec<f64> = r.solutions.iter().map(|x| x.value.unwrap()).collect();
    assert!(freqs.iter().any(|f| (f - f_ext).abs() < 0.03 * f_ext), "extensional {f_ext}: {freqs:?}");
    let s = spec(kz, 0.7 * f_flex, 6);
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let freqs: Vec<f64> = r.solutions.iter().map(|x| x.value.unwrap()).collect();
    // two degenerate flexural modes (x and y bending)
    let close = freqs.iter().filter(|f| (**f - f_flex).abs() < 0.06 * f_flex).count();
    assert!(close >= 2, "flexural {f_flex}: {freqs:?}");
}

#[test]
fn pml_in_the_cross_section_gives_the_analytic_complex_frequency() {
    // axial modes of a bar clamped at x = 0 with a complex stretch s = 1 - 0.1 j: f_n = (2n-1) c / (4 s L),
    // Q = 5; with kz = 0 the u_z (shear) modes appear as well, so pick the pair that matches
    let (e, rho, l) = (200e9f64, 7850.0f64, 1.0f64);
    let mesh = common::rect(l, 0.1, 80, 2);
    let mut s = ProblemSpec::new(Physics::Solid);
    for (k, v) in [("E", e), ("nu", 0.0), ("rho", rho)] {
        s.materials.default.insert(k.into(), v);
    }
    s.bcs = vec![
        Bc { on: None, labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { on: None, labels: vec![2, 3], kind: BcKind::Fixed { fixed: [false, true, false], values: [0.0; 3] } },
    ];
    s.materials.pml.insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.1, order: 0.0, re: 1.0 });
    let cs = (e / rho).sqrt();
    let want = cs / (4.0 * l) / 1.01; // first axial mode, Re f = f0 / |s|^2 * Re(1/s...) -> f0 (1 + 0.1j) / 1.01
    s.study = Study::Eigenfrequency { n_modes: 6, f_search: want };
    s.out_of_plane = Some(0.0);
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let hit = r.solutions.iter().position(|x| (x.value.unwrap() / want - 1.0).abs() < 2e-3);
    let hit = hit.unwrap_or_else(|| panic!("no mode near {want}: {:?}", r.solutions.iter().map(|x| x.label.clone()).collect::<Vec<_>>()));
    let q = r.scalar(&format!("Q{}", hit + 1)).unwrap();
    assert!((q / 5.0 - 1.0).abs() < 3e-3, "Q = {q}");
}

#[test]
fn a_floquet_periodic_3d_cell_reproduces_the_cross_section_solver() {
    // a free square rod: a 3D cell of length l along y, periodic in y with phase k l, against the 2D solver at kz
    use femcore::mesher::{extrude_y, Layer};
    let a = 0.01;
    let m2 = common::rect(a, a, 6, 6);
    let l = 0.01;
    let kz = 20.0;
    let mesh3 = extrude_y(&m2, &[Layer { thickness: l, n: 2, domain_offset: 0 }], 0.0).unwrap();
    let c_ext = (70e9f64 / 2700.0).sqrt();
    let f_ext = c_ext * kz / (2.0 * std::f64::consts::PI);
    let mut s3 = spec(kz, f_ext, 6);
    s3.out_of_plane = None;
    s3.periodic = vec![Periodic { source: vec![5], dest: vec![6], phase: kz * l }];
    let r3 = femcore::solve::solve(&mesh3, &s3, &mut SolveState::default()).unwrap();
    let s2 = spec(kz, f_ext, 6);
    let r2 = femcore::solve::solve(&m2, &s2, &mut SolveState::default()).unwrap();
    let f3: Vec<f64> = r3.solutions.iter().map(|s| s.value.unwrap()).collect();
    let f2: Vec<f64> = r2.solutions.iter().map(|s| s.value.unwrap()).collect();
    assert!(f3.iter().any(|f| (f / f_ext - 1.0).abs() < 0.04), "3D {f3:?} vs {f_ext}");
    // the modes of the 3D cell near the extensional wave appear in the 2D list too
    for f in f3.iter().filter(|f| (**f / f_ext - 1.0).abs() < 0.04) {
        assert!(f2.iter().any(|g| (f / g - 1.0).abs() < 0.03), "3D {f} not in 2D {f2:?}");
    }
}
