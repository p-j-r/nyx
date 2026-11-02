//! Optical modes of a crystal whose axes are tilted against the propagation axis: the permittivity has xz/yz terms,
//! the 2D mode problem becomes quadratic in the propagation constant and is solved exactly. Checked against the 3D
//! eigenmode solver on a periodic cell of the same PEC guide (a different formulation sharing only the tensor). Both are written in the
//! wafer axes: the guide runs along y, the crystal axis is tilted in the y-z plane against the propagation direction, the 2D problem goes
//! through `solve::solve` (so the frame adapter turns the crystal into the kernel's axes) and the 3D cell is periodic along y.
//! The last test drives the 2D kernel directly in its own axes (`emw::solve_emw`).
mod common;
use femcore::material::*;
use femcore::problem::*;

const C0: f64 = 299_792_458.0;

fn crystal(tilt_deg: f64) -> Anisotropy {
    Anisotropy { permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }), euler_deg: [0.0, tilt_deg, 0.0], ..Default::default() }
}

#[test]
fn a_tilted_uniaxial_crystal_in_a_pec_guide_has_the_same_modes_in_2d_and_in_a_3d_periodic_cell() {
    let (a, b, l, phase) = (0.03, 0.015, 0.01, 1.5);
    let beta = phase / l;
    let tilt = 40.0;
    // 3D cell
    let mesh3 = common::block(a, l, b, 6, 2, 3);
    let mut s3 = ProblemSpec::new(Physics::Emw);
    s3.order = 2;
    s3.materials.aniso.insert(1, crystal(tilt));
    s3.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase }];
    let f_est = C0 / (2.0 * std::f64::consts::PI * 4.8f64.sqrt()) * ((std::f64::consts::PI / a).powi(2) + beta * beta).sqrt();
    s3.study = Study::Eigenfrequency { n_modes: 4, f_search: 0.9 * f_est };
    let r3 = femcore::solve::solve(&mesh3, &s3, &mut SolveState::default()).unwrap();
    let f3: Vec<f64> = r3.solutions.iter().map(|x| x.value.unwrap()).collect();
    println!("3D frequencies (GHz): {:?}", f3.iter().map(|f| f / 1e9).collect::<Vec<_>>());
    // 2D at every 3D frequency: some mode must have beta
    let mesh2 = common::rect(a, b, 24, 12);
    let mut matched = 0;
    for &f in &f3 {
        let mut s2 = ProblemSpec::new(Physics::Emw);
        s2.order = 2;
        s2.materials.aniso.insert(1, crystal(tilt));
        let k0 = 2.0 * std::f64::consts::PI * f / C0;
        s2.study = Study::Modes { freq: f, n_target: beta / k0, evanescent: false, n_modes: 3 };
        let r2 = femcore::solve::solve(&mesh2, &s2, &mut SolveState::default()).unwrap();
        let betas: Vec<f64> = r2.solutions.iter().map(|s| s.neff.unwrap().0 * k0).collect();
        println!("f = {:.4} GHz: 2D betas {betas:.2?} (cell: {beta:.2})", f / 1e9);
        if betas.iter().any(|x| (x / beta - 1.0).abs() < 4e-3) {
            matched += 1;
        }
    }
    assert!(matched >= 2, "2D and 3D agree on {matched} of {} modes", f3.len());
}

#[test]
fn the_coupled_solver_reduces_to_the_uncoupled_one_without_tilt() {
    // a crystal with its axis along z (no xz/yz terms) goes through the real solver; forcing the complex path with a PML-free
    // complex solver must give the same effective indices
    let mesh = common::rect(0.03, 0.015, 12, 6);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.aniso.insert(1, crystal(0.0));
    s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 3 };
    let real = femcore::emw::solve_emw(&mesh, &s).unwrap();
    let (eps, mu) = femcore::emw_pml::cell_materials(&mesh, &s.materials).unwrap();
    let cm = femcore::emw_pml::solve_modes_complex(&mesh, 12e9, &eps, &mu, femcore::csparse::c(1.6, 0.0), 3, 2).unwrap();
    for (r, z) in real.solutions.iter().zip(&cm.neff) {
        assert!((r.neff.unwrap().0 - z.re).abs() < 1e-8 && z.im.abs() < 1e-8, "{} against {}", r.neff.unwrap().0, z.re);
    }
}
