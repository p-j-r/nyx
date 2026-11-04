//! ROM sweep vs one full solve per frequency: 2-D damped Helmholtz-like problem, 120 × 120 unknowns.
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::nep::{NepTerm, NonlinearEigenProblem};
use femcore::linalg::sweep::{SweepOptions, direct_sweep, rom_sweep};
use femcore::linalg::{C64, Csr};

fn main() {
    let (nx, ny) = (120usize, 120usize);
    let n = nx * ny;
    let gamma = 0.02;
    let p = NonlinearEigenProblem::new(vec![
        NepTerm { matrix: laplacian_2d(nx, ny).to_complex(), f: Box::new(|_| C64::ONE) },
        NepTerm { matrix: Csr::identity(n), f: Box::new(move |z| -z * z + C64::new(0.0, gamma) * z) },
    ])
    .unwrap();
    let mut b = vec![C64::ZERO; n];
    b[(ny / 3) * nx + nx / 4] = C64::ONE;
    let mut out = vec![C64::ZERO; n];
    out[(2 * ny / 3) * nx + 3 * nx / 4] = C64::ONE;
    let omegas: Vec<f64> = (0..300).map(|i| 0.1 + 0.25 * i as f64 / 299.0).collect();
    let rom = rom_sweep(&p, &b, &omegas, std::slice::from_ref(&out), &SweepOptions::new(1e-7)).unwrap();
    let full = direct_sweep(&p, &b, &omegas, std::slice::from_ref(&out), Definiteness::General).unwrap();
    let scale = full.outputs.iter().map(|o| o[0].norm()).fold(0.0, f64::max);
    let err = rom.outputs.iter().zip(&full.outputs).map(|(a, b)| (a[0] - b[0]).norm() / scale).fold(0.0, f64::max);
    println!(
        "n = {n}, {} frequencies: ROM {} factorisations, basis {}, {:.2} s | direct {} factorisations, {:.2} s | speed-up {:.1}x, worst output error {err:.1e}",
        omegas.len(),
        rom.factorizations,
        rom.basis_dim,
        rom.seconds,
        full.factorizations,
        full.seconds,
        full.seconds / rom.seconds
    );
}
