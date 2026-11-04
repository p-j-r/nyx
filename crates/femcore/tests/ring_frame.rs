//! The optical ring modes in the wafer frame: with the time convention e^{j omega t} of the rest of the library, exp(-i m theta) with m > 0
//! carries power along +theta, i.e. +y (and m < 0 along -y). The elastic ring tests (frame of the model axes, frozen frequencies, torsion,
//! the sign of theta, the average about the ring axis) moved to `nyx-sms/tests/ring_frame.rs` with the solid kernel.
mod common;

#[test]
fn a_whispering_gallery_mode_with_positive_azimuthal_order_carries_power_along_plus_y() {
    // the silicon wire in silica on a ring of 20 mm radius: exp(-i m theta) with the time convention e^{j omega t} is a wave along +theta
    let um = 1e-6;
    let radius = 20e-3;
    let rect = |x0: f64, z0: f64, x1: f64, z1: f64, d: i32, h: Option<f64>| common::part(radius + x0 * um, z0 * um, radius + x1 * um, z1 * um, d, h.map(|v| v * um));
    let mesh = common::parts_mesh(&[rect(-1.2, -1.0, 1.2, 1.0, 1, None), rect(-0.225, -0.11, 0.225, 0.11, 2, Some(0.04))], 0.3 * um);
    let mut idx = std::collections::BTreeMap::new();
    idx.insert(1, 1.444);
    idx.insert(2, 3.476);
    let mat = femcore::ems_axisym::AxiMaterial::from_indices(&mesh, &idx);
    let f0 = 299_792_458.0 / 1.55e-6;
    let m = (2.35 * 2.0 * std::f64::consts::PI * f0 / 299_792_458.0 * radius).round();
    let (plus, minus) = (femcore::ems_axisym::solve_axisym_modes(&mesh, &mat, m, f0, 1, 2).unwrap(), femcore::ems_axisym::solve_axisym_modes(&mesh, &mat, -m, f0, 1, 2).unwrap());
    let (pp, pm) = (femcore::ems_axisym::power_along_theta(&mesh, &mat, &plus, 0), femcore::ems_axisym::power_along_theta(&mesh, &mat, &minus, 0));
    // the same ring both ways: the frequencies agree, the powers are opposite
    assert!((plus.freqs[0] - minus.freqs[0]).abs() < 1e-6 * plus.freqs[0], "{} against {}", plus.freqs[0], minus.freqs[0]);
    assert!(pp > 0.0 && pm < 0.0 && (pp + pm).abs() < 1e-3 * pp.abs(), "P(m) = {pp:e}, P(-m) = {pm:e}");
}
