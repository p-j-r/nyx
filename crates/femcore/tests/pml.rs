//! Perfectly matched layers and complex materials in the mode solver.
mod common;
use femcore::csparse::c;
use femcore::ems::{self, C0};
use femcore::ems_pml::*;
use femcore::problem::*;

fn k0(f: f64) -> f64 {
    2.0 * std::f64::consts::PI * f / C0
}

fn spec(f: f64, target: f64, n: usize) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Emw);
    s.study = Study::Modes { freq: f, n_target: target, evanescent: false, n_modes: n };
    s
}

fn closest(r: &SolveResult, re: f64) -> (f64, f64) {
    r.solutions.iter().map(|s| s.neff.unwrap()).min_by(|a, b| (a.0 - re).abs().partial_cmp(&(b.0 - re).abs()).unwrap()).unwrap()
}

#[test]
fn complex_solver_without_pml_reproduces_the_real_solver() {
    let (a, b, f) = (22.86e-3, 10.16e-3, 20e9);
    let mesh = common::rect(a, b, 20, 10);
    let nc = mesh.n_cells();
    let eps = vec![EpsC { xx: c(1.0, 0.0), xy: c(0.0, 0.0), yy: c(1.0, 0.0), zz: c(1.0, 0.0), xz: c(0.0, 0.0), yz: c(0.0, 0.0) }; nc];
    let mu = vec![MuC { x: c(1.0, 0.0), y: c(1.0, 0.0), z: c(1.0, 0.0) }; nc];
    let cm = solve_modes_complex(&mesh, f, &eps, &mu, c(0.9, 0.0), 4, 2).unwrap();
    let rm = ems::solve_modes(&mesh, f, &vec![1.0; nc], &vec![1.0; nc], 0.9, false, 4, 2).unwrap();
    for j in 0..4 {
        assert!((cm.neff[j].re - rm.neff[j].0).abs() < 1e-9 && cm.neff[j].im.abs() < 1e-9, "mode {j}: {:?} vs {:?}", cm.neff[j], rm.neff[j]);
    }
}

#[test]
fn real_stretch_equals_a_wider_guide() {
    // stretching a guide of width a by s = 1.5 in x acts like a guide of width 1.5 a
    let (a, b, f) = (15e-3, 10e-3, 15e9);
    let mesh = common::parts_mesh(&[common::part(0.0, 0.0, a, b, 2, None)], (a / 24.0).min(b / 12.0));
    let mut s = spec(f, 0.9, 4);
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.0, order: 0.0, re: 1.5 });
    s.materials.set(1, "eps_r", 1.0);
    s.materials.set(2, "eps_r", 1.0);
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let want = (1.0 - (std::f64::consts::PI / (1.5 * a) / k0(f)).powi(2)).sqrt();
    let (re, im) = closest(&r, want);
    assert!((re - want).abs() < 1e-3 && im.abs() < 1e-9, "{re} {im} vs {want}");
}

#[test]
fn complex_constant_stretch_gives_the_analytic_complex_index() {
    // s = 1 - 0.3 j along x: beta^2 = k0^2 - (pi / (a s))^2
    let (a, b, f) = (15e-3, 10e-3, 15e9);
    let mesh = common::parts_mesh(&[common::part(0.0, 0.0, a, b, 2, None)], (a / 24.0).min(b / 12.0));
    let mut s = spec(f, 0.9, 4);
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.3, order: 0.0, re: 1.0 });
    let r = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    let sx = c(1.0, -0.3);
    let kc = c(std::f64::consts::PI / a, 0.0) / sx;
    let neff = (c(k0(f) * k0(f), 0.0) - kc * kc).sqrt() / k0(f);
    let best = r.solutions.iter().map(|s| {
        let z = s.neff.unwrap();
        ((z.0 - neff.re).powi(2) + (z.1 - neff.im).powi(2)).sqrt()
    }).fold(f64::MAX, f64::min);
    assert!(best < 2e-3, "analytic neff {neff:?}; nearest found {best} away");
}

#[test]
fn pml_terminated_slab_keeps_the_guided_mode_and_absorbs_the_rest() {
    // a dielectric core (eps 4) in an air-filled PEC box 40 x 10 mm; PML layers on both x ends
    let (lx, ly, f) = (40e-3, 10e-3, 12e9);
    let mesh = common::parts_mesh(
        &[common::part(0.0, 0.0, lx, ly, 1, None), common::part(15e-3, 0.0, 25e-3, ly, 2, None), common::part(0.0, 0.0, 8e-3, ly, 3, None), common::part(32e-3, 0.0, lx, ly, 4, None)],
        (lx / 80.0).min(ly / 8.0),
    );
    let mut pml = spec(f, 1.9, 6);
    pml.materials.default.insert("eps_r".into(), 1.0);
    pml.materials.set(2, "eps_r", 4.0);
    pml.materials.pml.insert(3, Pml { x: PmlDir::Minus, y: PmlDir::Off, z: PmlDir::Off, strength: 4.0, order: 2.0, re: 1.0 });
    pml.materials.pml.insert(4, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 4.0, order: 2.0, re: 1.0 });
    let with = femcore::solve::solve(&mesh, &pml, &mut SolveState::default()).unwrap();
    // reference: the same box without PML (PEC walls), real solver
    let mut plain = spec(f, 1.9, 6);
    plain.materials.default.insert("eps_r".into(), 1.0);
    plain.materials.set(2, "eps_r", 4.0);
    let without = femcore::solve::solve(&mesh, &plain, &mut SolveState::default()).unwrap();
    let guided = without.solutions[0].neff.unwrap().0;
    let (re, im) = closest(&with, guided);
    assert!((re - guided).abs() < 2e-3, "guided mode: {re} with PML vs {guided} without");
    assert!(im.abs() < 2e-3, "a guided mode should hardly be damped, Im = {im}");
    // a mode that lives in the cladding (neff near 1) is absorbed: much larger imaginary part
    let lossy = with.solutions.iter().map(|s| s.neff.unwrap().1.abs()).fold(0.0, f64::max);
    assert!(lossy > 10.0 * im.abs().max(1e-6) || with.solutions.len() < 6, "no absorbed mode found among {:?}", with.solutions.iter().map(|s| s.neff).collect::<Vec<_>>());
}

#[test]
fn pml_with_off_axes_changes_nothing() {
    let mesh = common::rect(22.86e-3, 10.16e-3, 16, 8);
    let mut s = spec(20e9, 0.9, 3);
    s.materials.pml.insert(1, Pml::default()); // both axes Off
    let a = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    s.materials.pml.clear();
    let b = femcore::solve::solve(&mesh, &s, &mut SolveState::default()).unwrap();
    for (x, y) in a.solutions.iter().zip(&b.solutions) {
        assert_eq!(x.neff, y.neff);
    }
}
