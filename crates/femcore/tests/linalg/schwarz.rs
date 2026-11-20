// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Additive and restricted Schwarz: exactness for one subdomain, effect of overlap, two-level
//! scalability, RAS with GMRES, complex symmetric COCG and independence from the worker count.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::direct::Definiteness;
use femcore::linalg::ksp::{Cg, CgKind, CgOptions, Gmres, GmresOptions, Tolerances};
use femcore::linalg::linop::Preconditioner;
use femcore::linalg::rng::random_vector;
use femcore::linalg::schwarz::{Schwarz, SchwarzKind, SchwarzOptions, contiguous_partition, graph_partition};
use femcore::linalg::{C64, Csr};

fn convdiff(nx: usize, bx: f64, by: f64) -> Csr<f64> {
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

fn cg_its<M: Preconditioner<f64>>(a: &Csr<f64>, m: &M) -> usize {
    let n = a.nrows();
    let b = random_vector::<f64>(n, 9);
    let mut x = vec![0.0; n];
    let info =
        Cg::new(CgOptions { tol: Tolerances { rtol: 1e-8, max_it: 2000, ..Default::default() }, ..Default::default() })
            .solve(a, m, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?}", info.reason);
    info.iterations
}

#[test]
fn one_subdomain_is_a_direct_solver() {
    let a = convdiff(10, 30.0, 10.0);
    let n = a.nrows();
    let s = Schwarz::new(
        &a,
        &vec![0; n],
        &SchwarzOptions { overlap: 0, kind: SchwarzKind::Restricted, ..Default::default() },
    )
    .unwrap();
    let b = random_vector::<f64>(n, 1);
    let mut z = vec![0.0; n];
    s.apply(&b, &mut z);
    let mut r = vec![0.0; n];
    a.spmv(&z, &mut r);
    assert!(r.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-11));
    assert_eq!(s.subdomains(), 1);
}

#[test]
fn overlap_reduces_iterations_and_sizes_grow() {
    let a = laplacian_2d(48, 48);
    let part = contiguous_partition(a.nrows(), 8);
    let mut its = Vec::new();
    let mut sizes = Vec::new();
    for ov in [0usize, 1, 2, 4] {
        let s = Schwarz::new(
            &a,
            &part,
            &SchwarzOptions {
                overlap: ov,
                kind: SchwarkBasic(),
                definiteness: Definiteness::PositiveDefinite,
                ..Default::default()
            },
        )
        .unwrap();
        sizes.push(s.sizes().iter().sum::<usize>());
        its.push(cg_its(&a, &s));
    }
    eprintln!("ASM(8 contiguous blocks) CG iterations by overlap 0/1/2/4: {its:?}; total subdomain sizes {sizes:?}");
    assert!(its.windows(2).all(|w| w[1] <= w[0]), "{its:?}");
    assert!(its[3] * 10 < its[0] * 7, "{its:?}");
    assert!(sizes.windows(2).all(|w| w[1] > w[0]));
}

#[allow(non_snake_case)]
fn SchwarkBasic() -> SchwarzKind {
    SchwarzKind::Basic
}

#[test]
fn coarse_space_removes_the_dependence_on_the_number_of_subdomains() {
    let a = laplacian_2d(64, 64);
    let n = a.nrows();
    let mut one = Vec::new();
    let mut two = Vec::new();
    for p in [4usize, 16, 64] {
        let part = graph_partition(&a, p);
        assert_eq!(part.iter().copied().max().unwrap() + 1, p);
        let base = SchwarzOptions {
            overlap: 2,
            kind: SchwarzKind::Basic,
            definiteness: Definiteness::PositiveDefinite,
            ..Default::default()
        };
        one.push(cg_its(&a, &Schwarz::new(&a, &part, &base).unwrap()));
        two.push(cg_its(&a, &Schwarz::new(&a, &part, &SchwarzOptions { coarse: true, ..base }).unwrap()));
    }
    eprintln!("CG iterations for 4/16/64 subdomains: one level {one:?}, two level {two:?} (n = {n})");
    assert!(one[2] > 2 * one[0], "one-level iterations should grow with the number of subdomains: {one:?}");
    assert!(two[2] * 10 <= two[0] * 18, "two-level iterations should stay about constant: {two:?}");
    assert!(two[2] < one[2]);
}

#[test]
fn restricted_schwarz_with_gmres_on_convection_diffusion() {
    let a = convdiff(60, 40.0, 25.0);
    let n = a.nrows();
    let part = graph_partition(&a, 9);
    let b = random_vector::<f64>(n, 5);
    let its = |kind: SchwarzKind, ov: usize| {
        let s = Schwarz::new(&a, &part, &SchwarzOptions { overlap: ov, kind, ..Default::default() }).unwrap();
        let mut x = vec![0.0; n];
        let info = Gmres::new(GmresOptions {
            tol: Tolerances { rtol: 1e-8, max_it: 2000, ..Default::default() },
            restart: 50,
            ..Default::default()
        })
        .solve(&a, &s, &b, &mut x);
        assert!(info.reason.is_converged(), "{:?}", info.reason);
        let mut ax = vec![0.0; n];
        a.spmv(&x, &mut ax);
        assert!(ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() < 1e-6 * norm(&b));
        info.iterations
    };
    let (ras1, asm1, ras3) =
        (its(SchwarzKind::Restricted, 1), its(SchwarzKind::Basic, 1), its(SchwarzKind::Restricted, 3));
    eprintln!("GMRES(50) its: RAS(1) {ras1}, ASM(1) {asm1}, RAS(3) {ras3}");
    assert!(ras1 <= asm1, "RAS should not need more iterations than ASM: {ras1} vs {asm1}");
    assert!(ras3 <= ras1);
}

#[test]
fn complex_symmetric_matrix_with_cocg() {
    let n = 40 * 40;
    let l = laplacian_2d(40, 40).to_complex();
    let a = Csr::axpby(C64::ONE, &l, C64::new(-0.05, 0.4), &Csr::identity(n));
    let part = graph_partition(&a, 6);
    let s = Schwarz::new(
        &a,
        &part,
        &SchwarzOptions { overlap: 2, kind: SchwarzKind::Basic, coarse: true, ..Default::default() },
    )
    .unwrap();
    let b = random_vector::<C64>(n, 3);
    let mut x = vec![C64::ZERO; n];
    let info = Cg::new(CgOptions {
        kind: CgKind::Symmetric,
        tol: Tolerances { rtol: 1e-8, max_it: 500, ..Default::default() },
        ..Default::default()
    })
    .solve(&a, &s, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    let mut ax = vec![C64::ZERO; n];
    a.spmv(&x, &mut ax);
    assert!(ax.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() < 1e-5 * norm(&b));
}

#[test]
fn result_does_not_depend_on_the_number_of_workers() {
    let a = laplacian_2d(30, 30);
    let n = a.nrows();
    let part = graph_partition(&a, 7);
    let r = random_vector::<f64>(n, 2);
    let mut out = Vec::new();
    for w in [1usize, 2, 3, 7] {
        let s = Schwarz::new(&a, &part, &SchwarzOptions { workers: w, overlap: 2, coarse: true, ..Default::default() })
            .unwrap();
        let mut z = vec![0.0; n];
        s.apply(&r, &mut z);
        out.push(z);
    }
    for z in &out[1..] {
        assert!(z.iter().zip(&out[0]).all(|(a, b)| a.to_bits() == b.to_bits()), "worker count changed the result");
    }
}

#[test]
fn partitions_and_invalid_inputs() {
    let a = laplacian_2d(10, 10);
    let c = contiguous_partition(100, 4);
    assert_eq!(c.iter().filter(|&&p| p == 0).count(), 25);
    let g = graph_partition(&a, 4);
    assert!((0..4).all(|p| g.iter().filter(|&&q| q == p).count() > 10));
    assert!(Schwarz::new(&a, &[0usize; 5], &SchwarzOptions::default()).is_err());
}
