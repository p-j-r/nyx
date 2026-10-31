//! Incomplete factorisations: exactness when no fill is dropped, iteration-count ordering of the
//! variants, complex symmetric and Hermitian incomplete Cholesky, and the diagonal-shift fallback.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::solvers::csr::{laplacian_1d, laplacian_2d};
use femcore::solvers::ilu::{Ic, Ilu};
use femcore::solvers::ksp::{Cg, CgKind, CgOptions, Gmres, GmresOptions, Tolerances};
use femcore::solvers::linop::{Identity, Preconditioner};
use femcore::solvers::pc::Jacobi;
use femcore::solvers::rng::random_vector;
use femcore::solvers::{C64, Csr, Scalar};

fn convdiff(nx: usize, bx: f64, by: f64) -> Csr<f64> {
    // 5-point convection–diffusion on an nx × nx grid, h = 1/(nx+1), centred convection
    let h = 1.0 / (nx as f64 + 1.0);
    let mut t = Vec::new();
    for j in 0..nx {
        for i in 0..nx {
            let p = j * nx + i;
            t.push((p, p, 4.0));
            if i > 0 {
                t.push((p, p - 1, -1.0 - 0.5 * h * bx));
            }
            if i + 1 < nx {
                t.push((p, p + 1, -1.0 + 0.5 * h * bx));
            }
            if j > 0 {
                t.push((p, p - nx, -1.0 - 0.5 * h * by));
            }
            if j + 1 < nx {
                t.push((p, p + nx, -1.0 + 0.5 * h * by));
            }
        }
    }
    Csr::from_triplets(nx * nx, nx * nx, &t).unwrap()
}

fn apply_a<T: Scalar>(a: &Csr<T>, x: &[T]) -> Vec<T> {
    let mut y = vec![T::ZERO; a.nrows()];
    a.spmv(x, &mut y);
    y
}

fn cg_its<M: Preconditioner<f64>>(a: &Csr<f64>, m: &M) -> usize {
    let n = a.nrows();
    let b = random_vector::<f64>(n, 4);
    let mut x = vec![0.0; n];
    let info =
        Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, max_it: 5000, ..Default::default() }, ..Default::default() })
            .solve(a, m, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?}", info.reason);
    info.iterations
}

fn gmres_its<M: Preconditioner<f64>>(a: &Csr<f64>, m: &M) -> usize {
    let n = a.nrows();
    let b = random_vector::<f64>(n, 5);
    let mut x = vec![0.0; n];
    let info = Gmres::new(GmresOptions {
        tol: Tolerances { rtol: 1e-8, max_it: 3000, ..Default::default() },
        restart: 50,
        ..Default::default()
    })
    .solve(a, m, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?}", info.reason);
    info.iterations
}

#[test]
fn ilu0_of_a_tridiagonal_matrix_is_the_exact_lu() {
    let a = laplacian_1d(50);
    let f = Ilu::ilu_k(&a, 0);
    let b = random_vector::<f64>(50, 1);
    let mut z = vec![0.0; 50];
    f.solve(&b, &mut z);
    let r = apply_a(&a, &z);
    assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-12));
    assert_eq!(f.perturbed_pivots, 0);
}

#[test]
fn ilu_with_enough_levels_is_exact_on_a_2d_grid() {
    let a = convdiff(8, 20.0, -10.0);
    let n = a.nrows();
    let f = Ilu::ilu_k(&a, 1000);
    let b = random_vector::<f64>(n, 1);
    let mut z = vec![0.0; n];
    f.solve(&b, &mut z);
    let r = apply_a(&a, &z);
    let err = r.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    assert!(err < 1e-12, "{err}");
    // fill grows with k
    let nz: Vec<usize> = [0, 1, 2, 3].iter().map(|&k| Ilu::ilu_k(&a, k).nnz()).collect();
    assert!(nz.windows(2).all(|w| w[0] < w[1]), "{nz:?}");
    assert_eq!(nz[0], a.nnz());
}

#[test]
fn ic_iteration_counts_improve_with_fill() {
    let a = laplacian_2d(60, 60);
    let jac = cg_its(&a, &Jacobi::new(&a));
    let ic0 = cg_its(&a, &Ic::ic_k(&a, 0, true));
    let ic1 = cg_its(&a, &Ic::ic_k(&a, 1, true));
    let ic3 = cg_its(&a, &Ic::ic_k(&a, 3, true));
    eprintln!("CG its: Jacobi {jac}, IC(0) {ic0}, IC(1) {ic1}, IC(3) {ic3}");
    assert!(ic0 < jac && ic1 <= ic0 && ic3 < ic1, "{jac} {ic0} {ic1} {ic3}");
    assert_eq!(Ic::ic_k(&a, 0, true).shift, 0.0);
    let _ = Identity;
}

#[test]
fn ilu_and_ilut_precondition_gmres_on_convection_diffusion() {
    let a = convdiff(60, 40.0, 25.0);
    let none = gmres_its(&a, &Identity);
    let i0 = gmres_its(&a, &Ilu::ilu_k(&a, 0));
    let i2 = gmres_its(&a, &Ilu::ilu_k(&a, 2));
    let t = Ilu::ilut(&a, 1e-3, 20);
    let it = gmres_its(&a, &t);
    eprintln!(
        "GMRES(50) its: none {none}, ILU(0) {i0}, ILU(2) {i2}, ILUT(1e-3, 20) {it} (nnz {} vs A {})",
        t.nnz(),
        a.nnz()
    );
    assert!(i0 < none && i2 < i0 && it < i0, "{none} {i0} {i2} {it}");
    assert!(t.nnz() <= a.nnz() * 5);
}

#[test]
fn ilut_with_zero_droptol_and_big_fill_is_exact() {
    let a = convdiff(9, 30.0, 10.0);
    let n = a.nrows();
    let f = Ilu::ilut(&a, 0.0, n);
    let b = random_vector::<f64>(n, 8);
    let mut z = vec![0.0; n];
    f.solve(&b, &mut z);
    let r = apply_a(&a, &z);
    let err = r.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    assert!(err < 1e-11, "{err}");
}

#[test]
fn complex_symmetric_ic_drives_cocg_and_hermitian_ic_drives_cg() {
    // complex symmetric: L + iγ I − k² I with small damping (indefinite real part would break IC; use k² < λ_min)
    let n = 30 * 30;
    let l = laplacian_2d(30, 30).to_complex();
    let sym = Csr::axpby(C64::ONE, &l, C64::new(0.01, 0.3), &Csr::identity(n));
    let b = random_vector::<C64>(n, 3);
    let mut x = vec![C64::ZERO; n];
    let ic = Ic::ic_k(&sym, 1, false);
    let info = Cg::new(CgOptions {
        kind: CgKind::Symmetric,
        tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() },
        ..Default::default()
    })
    .solve(&sym, &ic, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    let r = apply_a(&sym, &x);
    let err = r.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&b);
    assert!(err < 1e-5, "{err}");
    let plain = {
        let mut x2 = vec![C64::ZERO; n];
        Cg::new(CgOptions {
            kind: CgKind::Symmetric,
            tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() },
            ..Default::default()
        })
        .solve(&sym, &Identity, &b, &mut x2)
        .iterations
    };
    assert!(info.iterations * 2 < plain, "{} vs {plain}", info.iterations);
    // Hermitian: gauge-transformed Laplacian
    let herm_pairs = [(0usize, 1usize)];
    let _ = herm_pairs;
    let mut t = Vec::new();
    for p in 0..n {
        t.push((p, p, C64::new(4.0, 0.0)));
        if p % 30 != 0 {
            t.push((p, p - 1, C64::new(0.0, 0.5).exp().scale(-1.0)));
        }
        if p % 30 != 29 {
            t.push((p, p + 1, C64::new(0.0, -0.5).exp().scale(-1.0)));
        }
        if p >= 30 {
            t.push((p, p - 30, C64::new(-1.0, 0.0)));
        }
        if p + 30 < n {
            t.push((p, p + 30, C64::new(-1.0, 0.0)));
        }
    }
    let h = Csr::from_triplets(n, n, &t).unwrap();
    assert!(h.is_hermitian(1e-14));
    let ich = Ic::ic_k(&h, 0, true);
    let mut xh = vec![C64::ZERO; n];
    let info =
        Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() }, ..Default::default() })
            .solve(&h, &ich, &b, &mut xh);
    assert!(info.reason.is_converged(), "{:?}", info.reason);
}

#[test]
fn ic_falls_back_to_a_diagonal_shift_when_a_pivot_breaks_down() {
    // a symmetric matrix with a (barely) negative direction under IC(0): strong off-diagonals
    let n = 20;
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, 1.0));
        for j in 0..n {
            if i != j && (i as isize - j as isize).abs() <= 2 {
                t.push((i, j, 0.7));
            }
        }
    }
    let a = Csr::from_triplets(n, n, &t).unwrap();
    let ic = Ic::ic_k(&a, 0, true);
    assert!(ic.shift > 0.0 && ic.attempts > 1, "shift {} after {} attempts", ic.shift, ic.attempts);
    // still a usable (finite) preconditioner
    let b = random_vector::<f64>(n, 2);
    let mut z = vec![0.0; n];
    ic.solve(&b, &mut z);
    assert!(z.iter().all(|v| v.is_finite()));
}
