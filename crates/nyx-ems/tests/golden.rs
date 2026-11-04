//! Frozen reference output of the ems mode kernels (ADR 10): the golden files in `tests/golden/` were taken through the pre-restructure
//! `solve::solve` (femcore tests/golden_frames) and moved here with the kernels; the kernels reproduce them (tolerance 1e-7 and 1e-8). The
//! two silicon-wire files are the exception: they never reproduced on this machine (R22), and gate G2 (option B: loose tolerances, no
//! reference platform) sets their tolerances to what the physics justifies: about 5e-5 relative on n_eff and 1e-3 of the maximum on the fields.
//! `NYX_WRITE_GOLDEN=<prefix> cargo test --release -p nyx-ems --test golden` rewrites them; do that only for a kernel known to be right.
mod common;
use femcore::material::*;
use femcore::props::{Pml, PmlDir};
use femcore::usecases::soi::SoiWire;
use nyx_ems::{Problem, Study};
use nyx_physics::testkit::{cell_fields_reproduce_the_plot_fields, golden};
use nyx_physics::{Discretisation, SolveResult};
use serde_json::Value;

fn dir() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

/// Solves and checks that the cell fields of the fields stage reproduce the plot fields.
fn checked(mesh: &Discretisation, p: &Problem) -> SolveResult {
    let r = common::solve(mesh, p).unwrap();
    cell_fields_reproduce_the_plot_fields(mesh, &r);
    r
}

/// The G2 tolerances for the silicon-wire files: the effective indices to 5e-5 relative, then the rest of the record (fields against their own
/// maximum, scalars) to 1e-3. The polarization lines are texts of rounded numbers; only the class of the mode is compared.
fn golden_g2(name: &str, r: &SolveResult) {
    let want: Value =
        serde_json::from_str(&std::fs::read_to_string(dir().join(format!("{name}.json"))).unwrap()).unwrap();
    let sols = want["solutions"].as_array().unwrap();
    assert_eq!(sols.len(), r.solutions.len(), "{name}: number of modes");
    for (i, (w, s)) in sols.iter().zip(&r.solutions).enumerate() {
        let (wn, gn) = (w["neff"].as_array().unwrap(), common::neff(s));
        let (wr, wi) = (wn[0].as_f64().unwrap(), wn[1].as_f64().unwrap());
        assert!(((wr - gn.0) / wr).abs() < 5e-5, "{name}: n_eff of mode {} is {} against {wr}", i + 1, gn.0);
        assert!((wi - gn.1).abs() < 5e-5 * wr, "{name}: Im n_eff of mode {} is {} against {wi}", i + 1, gn.1);
    }
    // the polarization lines: compare the class, which opens the text, and leave the rounded numbers to the tolerances above and below
    let mut r = r.clone();
    for (k, v) in r.scalars.iter_mut() {
        if let (true, nyx_physics::Scalar::Text(t)) = (k.starts_with("polarization"), v) {
            let w = want["scalars"][k.as_str()].as_str().unwrap().to_string();
            let class = |t: &str| t.split("  TE ").next().unwrap().to_string();
            assert_eq!(class(t), class(&w), "{name}: {k}");
            *t = w;
        }
    }
    golden(&dir(), name, &r, 1e-3);
}

#[test]
fn golden_soi_wire_modes() {
    let wire = SoiWire { n_modes: 2, ..Default::default() };
    let g = wire.mesh().unwrap();
    let r = common::solve(&g.mesh, &common::from_old(&wire.spec(&g))).unwrap();
    golden_g2("emw_soi_wire", &r);
}

#[test]
fn golden_soi_wire_with_sloped_wall_hybrid_modes() {
    let wire = SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), n_modes: 2, ..Default::default() };
    let g = wire.mesh().unwrap();
    let r = common::solve(&g.mesh, &common::from_old(&wire.spec(&g))).unwrap();
    golden_g2("emw_soi_hybrid", &r);
}

fn crystal(tilt_deg: f64) -> Anisotropy {
    Anisotropy {
        permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }),
        euler_deg: [0.0, tilt_deg, 0.0],
        ..Default::default()
    }
}

#[test]
fn golden_tilted_crystal_in_a_pec_guide() {
    // the coupled (xz/yz) path of the complex solver, then the uncoupled real path
    for (name, tilt) in [("emw_pec_tilted", 40.0), ("emw_pec_straight", 0.0)] {
        let mesh = common::rect(0.03, 0.015, 12, 6);
        let mut s = common::emw();
        s.materials.aniso.insert(1, crystal(tilt));
        s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 3 };
        golden(&dir(), name, &checked(&mesh, &s), 1e-7);
    }
}

#[test]
fn golden_pec_guide_with_a_pml() {
    let mesh = common::rect_split_x(0.03, 0.015, 0.02, 0.0025);
    let mut s = common::emw();
    s.materials.set(1, "eps_r", 4.0);
    s.materials.set(2, "eps_r", 1.0);
    // absorbing layer on the right (domain 2) along x
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Pml::default() });
    s.study = Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 3 };
    golden(&dir(), "emw_pml_x", &checked(&mesh, &s), 1e-7);
}
