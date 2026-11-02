//! Waveguide mode analysis against closed forms.
mod common;
use femcore::emw::*;
use femcore::problem::*;

fn emw_spec(eps: f64, freq: f64, target: f64, evanescent: bool, n: usize) -> ProblemSpec {
    let mut s = ProblemSpec::new(Physics::Emw);
    s.materials.default.insert("eps_r".into(), eps);
    s.study = Study::Modes { freq, n_target: target, evanescent, n_modes: n };
    s
}

fn real_neff(r: &SolveResult) -> Vec<f64> {
    r.solutions.iter().map(|s| s.neff.unwrap().0).collect()
}

const WR90: (f64, f64) = (22.86e-3, 10.16e-3);

/// neff^2 = eps_r - kc^2 / k0^2 of the (m, n) mode of a PEC rectangle filled with eps_r (TE_mn and TM_mn share it).
fn neff_squared(eps: f64, freq: f64, m: usize, n: usize) -> f64 {
    let k0 = 2.0 * std::f64::consts::PI * freq / C0;
    let kc2 = (m as f64 * std::f64::consts::PI / WR90.0).powi(2) + (n as f64 * std::f64::consts::PI / WR90.1).powi(2);
    eps - kc2 / (k0 * k0)
}

#[test]
fn filled_wr90_propagating_modes_match_the_closed_form() {
    // WR-90 section 22.86 x 10.16 mm filled with eps_r 2.1, 12 GHz: TE10, TE20, TE01 and the degenerate TE11/TM11
    let m = common::rect(WR90.0, WR90.1, 24, 12);
    for (mm, nn) in [(1, 0), (2, 0), (0, 1), (1, 1)] {
        let exact = neff_squared(2.1, 12e9, mm, nn).sqrt();
        let r = solve_emw(&m, &emw_spec(2.1, 12e9, exact * 1.001, false, 3)).unwrap();
        let best = real_neff(&r).iter().map(|n| (n - exact).abs()).fold(f64::MAX, f64::min);
        assert!(best < 3e-4, "mode ({mm},{nn}) with neff {exact}: nearest found {best} away");
    }
}

#[test]
fn filled_wr90_evanescent_modes_match_the_closed_form() {
    // below cutoff neff is imaginary, neff = -i sqrt(kc^2 / k0^2 - eps_r)
    let m = common::rect(WR90.0, WR90.1, 24, 12);
    let exact = neff_squared(2.1, 12e9, 3, 0).abs().sqrt();
    let r = solve_emw(&m, &emw_spec(2.1, 12e9, exact, true, 4)).unwrap();
    let best = r.solutions.iter().map(|s| (s.neff.unwrap().1.abs() - exact).abs()).fold(f64::MAX, f64::min);
    assert!(best < 3e-3, "TE30 attenuation {exact}: nearest found {best} away");
}

#[test]
fn ez_field_of_the_tm11_mode_is_sin_sin() {
    // Ez of the TM11 mode is sin(pi x / a) sin(pi y / b); the eigenvector's scale is arbitrary. TE11 has no Ez, so
    // any vector of the degenerate pair is either the TM11 field or zero
    let (a, b) = WR90;
    let m = common::rect(a, b, 24, 12);
    let eps = vec![2.1; m.n_cells()];
    let mur = vec![1.0; m.n_cells()];
    let exact = neff_squared(2.1, 12e9, 1, 1).sqrt();
    let modes = solve_modes(&m, 12e9, &eps, &mur, exact * 1.001, false, 2, 2).unwrap();
    let norm = |v: &Vec<f64>| v.iter().map(|x| x * x).sum::<f64>();
    let j = (0..modes.phi.len()).max_by(|&p, &q| norm(&modes.phi[p]).partial_cmp(&norm(&modes.phi[q])).unwrap()).unwrap();
    let x = m.dof_coordinates(2);
    let shape: Vec<f64> = x.chunks(2).map(|p| (std::f64::consts::PI * p[0] / a).sin() * (std::f64::consts::PI * p[1] / b).sin()).collect();
    let p = &modes.phi[j];
    let c = p.iter().zip(&shape).map(|(u, v)| u * v).sum::<f64>() / shape.iter().map(|v| v * v).sum::<f64>();
    let num: f64 = p.iter().zip(&shape).map(|(u, v)| (u - c * v).powi(2)).sum::<f64>().sqrt();
    let den: f64 = p.iter().map(|u| u * u).sum::<f64>().sqrt();
    assert!(num / den < 5e-3, "{}", num / den);
}

#[test]
fn empty_wr90_closed_form() {
    // empty WR-90 at 20 GHz: neff = sqrt(1 - (kc/k0)^2), kc^2 = (m pi/a)^2 + (n pi/b)^2
    let (a, b) = (22.86e-3, 10.16e-3);
    let m = common::rect(a, b, 24, 12);
    let r = solve_emw(&m, &emw_spec(1.0, 20e9, 0.95, false, 4)).unwrap();
    let k0 = 2.0 * std::f64::consts::PI * 20e9 / C0;
    let mut exact: Vec<f64> = vec![];
    for (mm, nn) in [(1, 0), (2, 0), (0, 1), (1, 1), (3, 0)] {
        let kc2 = (mm as f64 * std::f64::consts::PI / a).powi(2) + (nn as f64 * std::f64::consts::PI / b).powi(2);
        let v = 1.0 - kc2 / (k0 * k0);
        if v > 0.0 {
            exact.push(v.sqrt());
        }
    }
    let got = real_neff(&r);
    for e in exact.iter().take(3) {
        let best = got.iter().map(|g| (g - e).abs()).fold(f64::MAX, f64::min);
        assert!(best < 2e-4, "mode with neff {e}: nearest found {best} away ({got:?})");
    }
}

#[test]
fn slab_loaded_guide_matches_transcendental_root() {
    // WR-90-like guide half... eps_r 4 in a layer of width d next to the wall, the rest air: the
    // lowest-order LSE mode with Ey = 0 satisfies k1 cot(k1 d) + k2 cot(k2 (a - d)) = 0
    let (a, b, d, f) = (22.86e-3, 10.16e-3, 11.43e-3, 8e9);
    let m = common::parts_mesh(&[common::part(0.0, 0.0, a, b, 1, None), common::part(0.0, 0.0, d, b, 2, None)], (a / 36.0).min(b / 8.0));
    let mut s = ProblemSpec::new(Physics::Emw);
    s.materials.set(1, "eps_r", 1.0);
    s.materials.set(2, "eps_r", 4.0);
    s.study = Study::Modes { freq: f, n_target: 1.7, evanescent: false, n_modes: 3 };
    let r = solve_emw(&m, &s).unwrap();
    let k0 = 2.0 * std::f64::consts::PI * f / C0;
    // the TE-to-y mode with a dielectric slab at x in [0, d]: solve for beta by bisection on
    // tan(k1 d) / k1 = -tan(k2 (a-d)) / k2 with k1^2 = 4 k0^2 - b^2, k2^2 = k0^2 - b^2
    let g = |neff: f64| {
        let (k1, k2) = ((4.0 * k0 * k0 - (neff * k0).powi(2)).sqrt(), (k0 * k0 - (neff * k0).powi(2)).abs().sqrt());
        if neff * neff < 1.0 {
            // air region oscillatory
            (k1 * d).sin() / k1 * (k2 * (a - d)).cos() + (k2 * (a - d)).sin() / k2 * (k1 * d).cos()
        } else {
            // air region evanescent: sin -> sinh
            (k1 * d).sin() / k1 * (k2 * (a - d)).cosh() + (k2 * (a - d)).sinh() / k2 * (k1 * d).cos()
        }
    };
    // scan for sign changes of g between neff 1.0 and 2.0 and compare with found modes
    let mut roots = vec![];
    let mut prev = (1.0001, g(1.0001));
    let mut x = 1.0001;
    while x < 1.999 {
        x += 0.0005;
        let v = g(x);
        if prev.1 * v < 0.0 {
            let (mut lo, mut hi) = (prev.0, x);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if g(lo) * g(mid) <= 0.0 { hi = mid } else { lo = mid }
            }
            roots.push(0.5 * (lo + hi));
        }
        prev = (x, v);
    }
    assert!(!roots.is_empty());
    let got = real_neff(&r);
    let best = got.iter().map(|g| (g - roots[roots.len() - 1]).abs()).fold(f64::MAX, f64::min);
    assert!(best < 2e-3, "highest root {} not found among {got:?} (roots {roots:?})", roots[roots.len() - 1]);
}

#[test]
fn order_one_converges_toward_order_two() {
    let (a, b) = (22.86e-3, 10.16e-3);
    let m = common::rect(a, b, 24, 12);
    let k0 = 2.0 * std::f64::consts::PI * 20e9 / C0;
    let exact = (1.0 - (std::f64::consts::PI / a / k0).powi(2)).sqrt();
    let mut e = vec![];
    for order in [1usize, 2] {
        let mut s = emw_spec(1.0, 20e9, 0.945, false, 2);
        s.order = order;
        let r = solve_emw(&m, &s).unwrap();
        e.push(real_neff(&r).iter().map(|n| (n - exact).abs()).fold(f64::MAX, f64::min));
    }
    assert!(e[1] < e[0] * 0.2, "order 1 error {} vs order 2 error {}", e[0], e[1]);
}

#[test]
fn mode_fields_are_normalised_and_finite() {
    let m = common::rect(WR90.0, WR90.1, 12, 6);
    let r = solve_emw(&m, &emw_spec(2.1, 12e9, 0.54, false, 3)).unwrap();
    for s in &r.solutions {
        let e = s.field("|E|").unwrap();
        assert!((e.values.iter().cloned().fold(0.0, f64::max) - 1.0).abs() < 1e-12);
        assert!(s.fields.iter().all(|f| f.values.iter().all(|v| v.is_finite())));
        assert_eq!(e.values.len(), m.n_dofs(2));
    }
}
