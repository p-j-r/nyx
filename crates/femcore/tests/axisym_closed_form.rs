//! The axisymmetric electromagnetic kernel (`ems_axisym`) against a closed form (gate G1: paul asked that its implementation be checked; the solid
//! one moved to `nyx-sms/tests/axisym_closed_form.rs`). Printed relative errors are in the test output (`-- --nocapture`).
mod common;

const C0: f64 = 299_792_458.0;

/// A coaxial cavity with perfectly conducting walls, `a < r < b`, `0 < z < h`, vacuum: its TEM modes (E_r ~ sin(p pi z / h) / r, H_theta ~
/// cos(p pi z / h) / r, azimuthal order m = 0) have `f = p c / (2 h)` whatever `a` and `b` are.
#[test]
fn coaxial_cavity_tem_modes_have_the_closed_form_frequencies() {
    let (a, b, h) = (0.10, 0.20, 0.15);
    let mesh = common::parts_mesh(&[common::part(a, 0.0, b, h, 1, Some(0.02))], 0.02);
    let mat = femcore::ems_axisym::AxiMaterial::from_indices(&mesh, &Default::default());
    for p in [1.0, 2.0] {
        let exact = p * C0 / (2.0 * h);
        let modes = femcore::ems_axisym::solve_axisym_modes(&mesh, &mat, 0.0, exact * 1.001, 1, 2).unwrap();
        let rel = (modes.freqs[0] - exact).abs() / exact;
        println!("coax TEM p={p}: {:.6e} Hz against {exact:.6e} Hz, relative error {rel:.2e}", modes.freqs[0]);
        assert!(rel < 2e-4, "p = {p}: {} against {exact}", modes.freqs[0]);
    }
    // the same cavity at azimuthal order 1 has no TEM mode: nothing of it may sit at the TEM frequency
    let m1 = femcore::ems_axisym::solve_axisym_modes(&mesh, &mat, 1.0, C0 / (2.0 * h), 1, 2).unwrap();
    assert!((m1.freqs[0] - C0 / (2.0 * h)).abs() > 1e-3 * C0 / (2.0 * h), "{}", m1.freqs[0]);
}

