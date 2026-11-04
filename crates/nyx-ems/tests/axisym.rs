//! The axisymmetric kernel (bodies of revolution on the (r, z) half plane) against closed forms (gate G1: paul asked that its implementation be
//! checked; ADR 46) and the direction of the azimuthal wave.
mod common;
use nyx_ems::{solve_axisym_modes, AxiMaterial};

const C0: f64 = 299_792_458.0;

/// A PEC cylinder of radius `R` (a thin conducting rod of radius `a << R` on the axis keeps the mesh off the axis) and height `h`: its
/// `TM_{m,1,0}` modes have `E_z = J_m(x_{m,1} r / R) exp(-j m theta)`, uniform in `z`, and the frequency `f = c x_{m,1} / (2 pi R)` with
/// `x_{m,1}` the first zero of the Bessel function `J_m` (`x_{1,1} = 3.8317`, `x_{2,1} = 5.1356`, `x_{3,1} = 6.3802`). Unlike the coaxial
/// TEM modes this exercises the azimuthal terms of the curl (`m / r`) and the `E_theta` unknown at `m` of 1, 2 and 3. The rod perturbs the
/// frequency by about `(a/R)^2` (4e-4 here).
#[test]
fn a_pec_cylinder_cavity_has_the_bessel_zero_frequencies_at_nonzero_azimuthal_order() {
    let (rr, a, h) = (0.10, 0.002, 0.12);
    let mesh = common::parts_mesh(&[common::part(a, 0.0, rr, h, 1, Some(0.015))], 0.015);
    let mat = AxiMaterial::from_indices(&mesh, &Default::default());
    for (m, x) in [(1.0, 3.8317060), (2.0, 5.1356223), (3.0, 6.3801619)] {
        let exact = C0 * x / (2.0 * std::f64::consts::PI * rr);
        let modes = solve_axisym_modes(&mesh, &mat, m, exact * 1.002, 1, 2).unwrap();
        let rel = (modes.freqs[0] - exact).abs() / exact;
        println!("PEC cylinder TM_{m}10: {:.6e} Hz against {exact:.6e} Hz, relative error {rel:.2e}", modes.freqs[0]);
        assert!(rel < 2e-3, "m = {m}: {} against {exact}", modes.freqs[0]);
    }
}

/// A coaxial cavity with perfectly conducting walls, `a < r < b`, `0 < z < h`, vacuum: its TEM modes (E_r ~ sin(p pi z / h) / r, H_theta ~
/// cos(p pi z / h) / r, azimuthal order m = 0) have `f = p c / (2 h)` whatever `a` and `b` are.
#[test]
fn coaxial_cavity_tem_modes_have_the_closed_form_frequencies() {
    let (a, b, h) = (0.10, 0.20, 0.15);
    let mesh = common::parts_mesh(&[common::part(a, 0.0, b, h, 1, Some(0.02))], 0.02);
    let mat = nyx_ems::AxiMaterial::from_indices(&mesh, &Default::default());
    for p in [1.0, 2.0] {
        let exact = p * C0 / (2.0 * h);
        let modes = nyx_ems::solve_axisym_modes(&mesh, &mat, 0.0, exact * 1.001, 1, 2).unwrap();
        let rel = (modes.freqs[0] - exact).abs() / exact;
        println!("coax TEM p={p}: {:.6e} Hz against {exact:.6e} Hz, relative error {rel:.2e}", modes.freqs[0]);
        assert!(rel < 2e-4, "p = {p}: {} against {exact}", modes.freqs[0]);
    }
    // the same cavity at azimuthal order 1 has no TEM mode: nothing of it may sit at the TEM frequency
    let m1 = nyx_ems::solve_axisym_modes(&mesh, &mat, 1.0, C0 / (2.0 * h), 1, 2).unwrap();
    assert!((m1.freqs[0] - C0 / (2.0 * h)).abs() > 1e-3 * C0 / (2.0 * h), "{}", m1.freqs[0]);
}

/// The m = 1 error of the cylinder test is the rod, not the discretisation: halving the rod radius takes it to about a quarter.
#[test]
fn the_remaining_error_of_the_cylinder_mode_at_m_1_is_the_rod_on_the_axis() {
    let (rr, h, x) = (0.10, 0.12, 3.8317060);
    let exact = C0 * x / (2.0 * std::f64::consts::PI * rr);
    let err = |a: f64| {
        let mesh = common::parts_mesh(&[common::part(a, 0.0, rr, h, 1, Some(0.015))], 0.015);
        let mat = AxiMaterial::from_indices(&mesh, &Default::default());
        let f = solve_axisym_modes(&mesh, &mat, 1.0, exact * 1.002, 1, 2).unwrap().freqs[0];
        (f - exact).abs() / exact
    };
    let (e1, e2) = (err(0.002), err(0.001));
    println!("rod 2 mm: {e1:.2e}, rod 1 mm: {e2:.2e}");
    assert!(e2 < 0.45 * e1, "{e1:e} then {e2:e}");
}

#[test]
fn a_whispering_gallery_mode_with_positive_azimuthal_order_carries_power_along_plus_y() {
    // the silicon wire in silica on a ring of 20 mm radius: exp(-i m theta) with the time convention e^{j omega t} is a wave along +theta
    let um = 1e-6;
    let radius = 20e-3;
    let rect = |x0: f64, z0: f64, x1: f64, z1: f64, d: i32, h: Option<f64>| {
        common::part(radius + x0 * um, z0 * um, radius + x1 * um, z1 * um, d, h.map(|v| v * um))
    };
    let mesh = common::parts_mesh(
        &[rect(-1.2, -1.0, 1.2, 1.0, 1, None), rect(-0.225, -0.11, 0.225, 0.11, 2, Some(0.04))],
        0.3 * um,
    );
    let mut idx = std::collections::BTreeMap::new();
    idx.insert(1, 1.444);
    idx.insert(2, 3.476);
    let mat = nyx_ems::AxiMaterial::from_indices(&mesh, &idx);
    let f0 = 299_792_458.0 / 1.55e-6;
    let m = (2.35 * 2.0 * std::f64::consts::PI * f0 / 299_792_458.0 * radius).round();
    let (plus, minus) = (
        nyx_ems::solve_axisym_modes(&mesh, &mat, m, f0, 1, 2).unwrap(),
        nyx_ems::solve_axisym_modes(&mesh, &mat, -m, f0, 1, 2).unwrap(),
    );
    let (pp, pm) =
        (nyx_ems::power_along_theta(&mesh, &mat, &plus, 0), nyx_ems::power_along_theta(&mesh, &mat, &minus, 0));
    // the same ring both ways: the frequencies agree, the powers are opposite
    assert!(
        (plus.freqs[0] - minus.freqs[0]).abs() < 1e-6 * plus.freqs[0],
        "{} against {}",
        plus.freqs[0],
        minus.freqs[0]
    );
    assert!(pp > 0.0 && pm < 0.0 && (pp + pm).abs() < 1e-3 * pp.abs(), "P(m) = {pp:e}, P(-m) = {pm:e}");
}
