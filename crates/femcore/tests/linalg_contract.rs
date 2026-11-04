//! Contract of the linear-algebra driver (docs/architecture/adr/0031): the shape of the eigensolver options, the names the kernels call,
//! and closed-form answers that do not depend on frozen floating-point values. `NYX_WRITE_CONTRACT=1` rewrites the exemplar
//! (only with a line in `contracts/VERSIONS.md`).
use femcore::progress::{self, Progress};
use femcore::linalg::csr::{laplacian_1d, Csr};
use femcore::linalg::direct::Definiteness;
use femcore::linalg::driver::{dense_eigh, dense_gen_eigh, eigs_near, eigs_real_general_near, DirectSolver, EigenOptions};
use femcore::linalg::problem::Symmetry;
use femcore::linalg::scalar::C64;
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts").join(name)
}

#[test]
fn eigen_options_have_the_documented_shape_and_defaults() {
    let got = json!({ "default": EigenOptions::default(), "custom": EigenOptions { tol: 1e-12, block: 4, ncv: 40, par_solve_min_n: 5000 } });
    let path = contract_file("linalg_options.v1.json");
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::write(&path, format!("{}\n", serde_json::to_string_pretty(&got).unwrap())).unwrap();
        return;
    }
    let want: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    assert_eq!(want, got, "linalg_options changed: this is an interface change (ADR 19)");
    let back: EigenOptions = serde_json::from_value(got["custom"].clone()).unwrap();
    assert_eq!(back, EigenOptions { tol: 1e-12, block: 4, ncv: 40, par_solve_min_n: 5000 });
}

#[test]
fn the_driver_surface_kernels_call_is_public() {
    let _ = (dense_eigh, dense_gen_eigh, eigs_real_general_near);
    let _ = (DirectSolver::<f64>::new, DirectSolver::<C64>::new, eigs_near::<f64>, eigs_near::<C64>);
}

/// `K u = 1` for the 1D Laplacian (2, -1): `u_i = i (n + 1 - i) / 2`.
#[test]
fn a_direct_solve_reproduces_the_closed_form_poisson_solution() {
    let n = 200;
    let k = laplacian_1d(n);
    let s = DirectSolver::new(&k, true).unwrap();
    let u = s.solve(&vec![1.0; n]);
    for (i, v) in u.iter().enumerate() {
        let exact = (i + 1) as f64 * (n - i) as f64 / 2.0;
        assert!((v - exact).abs() <= 1e-9 * exact.max(1.0), "{i}: {v} vs {exact}");
    }
    // many right-hand sides give the same as one at a time
    let many = s.solve_many(&[vec![1.0; n], vec![2.0; n]]);
    assert_eq!(many[0], u);
    assert_eq!(many[1], u.iter().map(|x| 2.0 * x).collect::<Vec<_>>());
}

/// Eigenvalues of the 1D Laplacian: `2 - 2 cos(k pi / (n + 1))`.
#[test]
fn shift_invert_finds_the_closed_form_laplacian_eigenvalues() {
    let n = 300;
    let (a, b) = (laplacian_1d(n), Csr::<f64>::identity(n));
    let p = eigs_near(&a, &b, -0.01, 4, Symmetry::Hermitian, Definiteness::PositiveDefinite, &EigenOptions::default(), None).unwrap();
    assert!(p.converged);
    for (k, lam) in p.values.iter().enumerate() {
        let exact = 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
        assert!((lam.re - exact).abs() <= 1e-9 * exact, "{k}: {} vs {exact}", lam.re);
        assert!(p.residuals[k] * lam.re <= 1e-8, "residual {}", p.residuals[k]);
    }
    // a repeated run gives the same bits
    let q = eigs_near(&a, &b, -0.01, 4, Symmetry::Hermitian, Definiteness::PositiveDefinite, &EigenOptions::default(), None).unwrap();
    assert_eq!(p.values, q.values);
}

/// A symmetric tridiagonal (2, -1) matrix scaled by a diagonal similarity `D^-1 T D` (D from 1 to 2): not symmetric, same real eigenvalues
/// `2 - 2 cos(k pi / (n + 1))`, mildly non-normal.
#[test]
fn a_real_non_symmetric_pencil_gives_its_real_eigenvalues_with_real_vectors() {
    let n = 120;
    let d: Vec<f64> = (0..n).map(|i| 1.0 + i as f64 / n as f64).collect();
    let mut t = vec![];
    for i in 0..n {
        t.push((i, i, 2.0));
        if i > 0 {
            t.push((i, i - 1, -d[i - 1] / d[i]));
        }
        if i + 1 < n {
            t.push((i, i + 1, -d[i + 1] / d[i]));
        }
    }
    let a = Csr::from_triplets(n, n, &t).unwrap();
    let p = eigs_real_general_near(&a, &Csr::identity(n), 0.0, 3, &EigenOptions::default()).unwrap();
    assert_eq!(p.values.len(), 3);
    for (k, lam) in p.values.iter().enumerate() {
        let exact = 2.0 - 2.0 * ((k + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
        assert!((lam.re - exact).abs() <= 1e-7 * exact, "{k}: {} vs {exact}", lam.re);
        assert!(p.vectors[k].iter().all(|x| x.is_finite()));
    }
}

#[test]
fn dense_generalised_eigenproblem_matches_its_closed_form() {
    // K = diag(1,4,9), M = diag(1,1,4): eigenvalues 1, 4, 9/4, ascending 1, 2.25, 4
    let k = [1.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 9.0];
    let m = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 4.0];
    let (vals, vecs) = dense_gen_eigh(&k, &m, 3).unwrap();
    for (v, w) in vals.iter().zip([1.0, 2.25, 4.0]) {
        assert!((v - w).abs() < 1e-12);
    }
    let mnorm: f64 = vecs[1].iter().zip([1.0, 1.0, 4.0]).map(|(x, d)| x * x * d).sum();
    assert!((mnorm - 1.0).abs() < 1e-12);
    let (w2, _) = dense_eigh(&[2.0, 1.0, 1.0, 2.0], 2).unwrap();
    assert!((w2[0] - 1.0).abs() < 1e-12 && (w2[1] - 3.0).abs() < 1e-12);
}

#[test]
fn a_cancel_during_the_eigensolve_returns_the_cancelled_string() {
    let n = 400;
    let (a, b) = (laplacian_1d(n), Csr::<f64>::identity(n));
    let p = Progress::new();
    let _g = progress::install(&p);
    let pc = p.clone();
    // the projector runs after the first operator application and raises the flag; the next application must stop the run
    let project = move |_x: &mut Vec<f64>| pc.cancel();
    let r = eigs_near(&a, &b, -0.01, 3, Symmetry::Hermitian, Definiteness::PositiveDefinite, &EigenOptions::default(), Some(&project));
    assert_eq!(r.err().as_deref(), Some(progress::CANCELLED));
}

/// A real matrix whose eigenvalues nearest the shift are a complex-conjugate pair (the case of a leaky mode of a real mode pencil): the real
/// solver must report it (`complex_found`) rather than drop it, so that the caller solves in complex arithmetic. Eigenvalues: `1 +- 0.5 i` and
/// the real numbers 12, 13, ..., 69.
fn complex_pair_matrix(n: usize) -> Csr<f64> {
    let mut t = vec![(0, 0, 1.0), (0, 1, -0.5), (1, 0, 0.5), (1, 1, 1.0)];
    for k in 2..n {
        t.push((k, k, 10.0 + k as f64));
    }
    Csr::from_triplets(n, n, &t).unwrap()
}

#[test]
fn a_converged_complex_pair_of_a_real_pencil_is_reported_and_never_dropped() {
    use femcore::linalg::driver::{eigs_general, EigOpts};
    let n = 70;
    let (a, b) = (complex_pair_matrix(n), Csr::<f64>::identity(n));
    // shift at the pair: the two nearest eigenvalues are 1 +- 0.5i
    let r = eigs_general(&a, &b, 1.0, 2, &EigOpts::default()).unwrap();
    assert!(r.complex_found, "the pair 1 +- 0.5i is the nearest and must be reported");
    // shift among the real eigenvalues: the pair is far (rank well beyond the wanted ones) and is not what the caller asked for
    let r = eigs_general(&a, &b, 20.2, 2, &EigOpts::default()).unwrap();
    assert!(!r.complex_found);
    assert_eq!(r.values.len(), 2);
    assert!((r.values[0] - 20.0).abs() < 1e-8 && (r.values[1] - 21.0).abs() < 1e-8, "{:?}", r.values);
}
