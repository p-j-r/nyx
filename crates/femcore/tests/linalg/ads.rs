// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::ads::{Ads, AdsOptions};
use femcore::linalg::ksp::{Cg, CgKind, CgOptions, Tolerances};
use femcore::linalg::linop::{Identity, Preconditioner};
use femcore::linalg::maxwell::{FaceElementCube, nedelec_cube, raviart_thomas_cube};
use femcore::linalg::pc::Jacobi;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr};

/// `Dᵀ W D + β M_f`.
fn face_matrix(f: &FaceElementCube, beta: f64) -> Csr<f64> {
    Csr::axpby(1.0, &f.div_stiffness, beta, &f.mass)
}

fn build_ads(n: usize, beta: f64) -> (Csr<f64>, Ads<f64>) {
    let e = nedelec_cube(n);
    let f = raviart_thomas_cube(n);
    let a = face_matrix(&f, beta);
    let ads = Ads::new(&a, &f.curl, &e.gradient, &e.edge_vectors, &e.mass, &f.pi, &AdsOptions::default()).unwrap();
    (a, ads)
}

fn cg_iterations<M: Preconditioner<f64>>(a: &Csr<f64>, m: &M) -> (usize, bool, f64) {
    let n = a.nrows();
    let b = random_vector::<f64>(n, 11);
    let mut x = vec![0.0; n];
    let info =
        Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, max_it: 3000, ..Default::default() }, ..Default::default() })
            .solve(a, m, &b, &mut x);
    let mut ax = vec![0.0; n];
    a.spmv(&x, &mut ax);
    let res = ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    (info.iterations, info.reason.is_converged(), res)
}

#[test]
fn pcg_with_ads_is_nearly_independent_of_the_mesh_while_jacobi_degrades() {
    let (mut ads_its, mut jac_its) = (Vec::new(), Vec::new());
    for n in [4usize, 8, 12] {
        let (a, ads) = build_ads(n, 1.0);
        let (it, ok, res) = cg_iterations(&a, &ads);
        assert!(ok && res < 1e-6, "ADS-PCG failed for n = {n} (residual {res})");
        ads_its.push(it);
        let (it, ok, _) = cg_iterations(&a, &Jacobi::new(&a));
        assert!(ok);
        jac_its.push(it);
    }
    eprintln!("ADS {ads_its:?}   Jacobi {jac_its:?}");
    assert!(*ads_its.last().unwrap() <= 40, "{ads_its:?}");
    assert!(
        *ads_its.last().unwrap() as f64 <= 1.8 * ads_its[0] as f64 + 6.0,
        "iterations grow with the mesh: {ads_its:?}"
    );
    assert!(*jac_its.last().unwrap() >= 3 * ads_its.last().unwrap(), "Jacobi {jac_its:?} vs ADS {ads_its:?}");
}

#[test]
fn the_curl_correction_is_what_removes_the_null_space() {
    // β = 1e-4: the grad-div term dominates, the curls are an almost exact kernel
    let (a, ads) = build_ads(8, 1e-4);
    let (it_ads, ok, _) = cg_iterations(&a, &ads);
    assert!(ok, "ADS-PCG not converged");
    let (it_jac, _, _) = cg_iterations(&a, &Jacobi::new(&a));
    let (it_none, _, _) = cg_iterations(&a, &Identity);
    eprintln!("β = 1e-4: ADS {it_ads}, Jacobi {it_jac}, none {it_none}");
    assert!(it_ads <= 50 && it_jac >= 4 * it_ads && it_none >= 4 * it_ads);
}

#[test]
fn the_preconditioner_is_symmetric_and_positive() {
    let (a, ads) = build_ads(6, 0.5);
    let n = a.nrows();
    let (x, y) = (random_vector::<f64>(n, 1), random_vector::<f64>(n, 2));
    let (mut mx, mut my) = (vec![0.0; n], vec![0.0; n]);
    ads.apply(&x, &mut mx);
    ads.apply(&y, &mut my);
    let dot = |p: &[f64], q: &[f64]| p.iter().zip(q).map(|(a, b)| a * b).sum::<f64>();
    let (s1, s2) = (dot(&y, &mx), dot(&x, &my));
    assert!((s1 - s2).abs() < 1e-10 * s1.abs().max(1.0), "{s1} vs {s2}");
    assert!(dot(&x, &mx) > 0.0 && dot(&y, &my) > 0.0);
}

#[test]
fn complex_symmetric_face_matrix_with_cocg() {
    let e = nedelec_cube(8);
    let f = raviart_thomas_cube(8);
    let a: Csr<C64> = Csr::axpby(C64::ONE, &f.div_stiffness.to_complex(), C64::new(1.0, 0.4), &f.mass.to_complex());
    let ads = Ads::new(&a, &f.curl, &e.gradient, &e.edge_vectors, &e.mass, &f.pi, &AdsOptions::default()).unwrap();
    let n = a.nrows();
    let b = random_vector::<C64>(n, 3);
    let mut x = vec![C64::ZERO; n];
    let info = Cg::new(CgOptions {
        tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() },
        kind: CgKind::Symmetric,
        ..Default::default()
    })
    .solve(&a, &ads, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(info.iterations <= 60, "{}", info.iterations);
    let mut ax = vec![C64::ZERO; n];
    a.spmv(&x, &mut ax);
    let res = ax.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm(&b);
    assert!(res < 1e-6, "{res}");
}
