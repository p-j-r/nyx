// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Field-split and Schur-complement preconditioners on saddle-point and coupled problems.
#![allow(clippy::needless_range_loop)]
use crate::common;

use common::*;
use femcore::linalg::csr::laplacian_2d;
use femcore::linalg::direct::{Definiteness, DirectSolve};
use femcore::linalg::fieldsplit::{FieldSplit, SchurFact, SchurFieldSplit, SchurPre, SplitKind, direct_maker, extract_block};
use femcore::linalg::ksp::{Gmres, GmresOptions, Minres, MinresOptions, Tolerances};
use femcore::linalg::linop::{Identity, Preconditioner};
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Csr, Scalar};

/// `[[K, Bᵀ], [B, −εM]]`: K = 2-D Laplacian (n0 = nx²), B a sparse full-row-rank constraint matrix.
fn saddle(nx: usize, n1: usize, eps: f64) -> (Csr<f64>, Vec<usize>, Vec<usize>) {
    let n0 = nx * nx;
    assert!(3 * n1 <= n0);
    let k = laplacian_2d(nx, nx);
    let mut t = Vec::new();
    for i in 0..n0 {
        let (cs, vs) = k.row(i);
        for (&c, &v) in cs.iter().zip(vs) {
            t.push((i, c as usize, v));
        }
    }
    for r in 0..n1 {
        for (c, w) in [(3 * r, 1.0), (3 * r + 1, -1.0), (3 * r + 2, 0.5)] {
            t.push((n0 + r, c, w));
            t.push((c, n0 + r, w));
        }
        if eps != 0.0 {
            t.push((n0 + r, n0 + r, -eps));
        }
    }
    (Csr::from_triplets(n0 + n1, n0 + n1, &t).unwrap(), (0..n0).collect(), (n0..n0 + n1).collect())
}

fn gmres_its(a: &Csr<f64>, m: &dyn Preconditioner<f64>, rtol: f64) -> usize {
    let b = random_vector::<f64>(a.nrows(), 2);
    let mut x = vec![0.0; a.nrows()];
    let info = Gmres::new(GmresOptions {
        tol: Tolerances { rtol, max_it: 400, ..Default::default() },
        restart: 400,
        ..Default::default()
    })
    .solve(a, m, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {} its", info.reason, info.iterations);
    let mut ax = vec![0.0; a.nrows()];
    a.spmv(&x, &mut ax);
    let res = ax.iter().zip(&b).map(|(p, q)| (p - q).powi(2)).sum::<f64>().sqrt() / norm(&b);
    // left preconditioning monitors the preconditioned residual: the true one may be larger by ‖M‖
    assert!(res < 1e-4, "true relative residual {res:e} after {} iterations", info.iterations);
    info.iterations
}

/// The exact Schur complement `A11 − A10 A00⁻¹ A01` as a (dense-pattern) CSR matrix.
fn exact_schur(a: &Csr<f64>, i0: &[usize], i1: &[usize]) -> Csr<f64> {
    let a00 = extract_block(a, i0, i0);
    let (a01, a10, a11) = (extract_block(a, i0, i1), extract_block(a, i1, i0), extract_block(a, i1, i1));
    let solve = DirectSolve::new(&a00, Definiteness::PositiveDefinite).unwrap();
    let n1 = i1.len();
    let mut t = Vec::new();
    for j in 0..n1 {
        let mut col = vec![0.0; i0.len()];
        for r in 0..i0.len() {
            let (cs, vs) = a01.row(r);
            for (&c, &v) in cs.iter().zip(vs) {
                if c as usize == j {
                    col[r] = v;
                }
            }
        }
        let mut y = vec![0.0; i0.len()];
        solve.apply(&col, &mut y);
        let mut ay = vec![0.0; n1];
        a10.spmv(&y, &mut ay);
        for i in 0..n1 {
            let aij = {
                let (cs, vs) = a11.row(i);
                cs.iter().zip(vs).find(|(c, _)| **c as usize == j).map_or(0.0, |(_, v)| *v)
            };
            t.push((i, j, aij - ay[i]));
        }
    }
    Csr::from_triplets(n1, n1, &t).unwrap()
}

#[test]
fn exact_blocks_give_the_textbook_iteration_counts() {
    let (a, i0, i1) = saddle(10, 30, 0.0);
    let s = exact_schur(&a, &i0, &i1);
    let direct_spd = direct_maker::<f64>(Definiteness::General);
    let counts: Vec<(SchurFact, usize)> = [SchurFact::Diag, SchurFact::Lower, SchurFact::Upper, SchurFact::Full]
        .into_iter()
        .map(|f| {
            let pc = SchurFieldSplit::new(
                &a,
                i0.clone(),
                i1.clone(),
                f,
                SchurPre::User(s.clone()),
                false,
                &direct_spd,
                &direct_spd,
            )
            .unwrap();
            (f, gmres_its(&a, &pc, 1e-10))
        })
        .collect();
    eprintln!("GMRES iterations with exact blocks: {counts:?}");
    assert!(counts[0].1 <= 3, "{counts:?}");
    assert!(counts[1].1 <= 2 && counts[2].1 <= 2, "{counts:?}");
    assert!(counts[3].1 <= 1, "{counts:?}");
}

#[test]
fn selfp_schur_with_exact_a00_is_a_good_preconditioner_for_a_saddle_point_problem() {
    let (a, i0, i1) = saddle(14, 60, 0.0);
    let direct = direct_maker::<f64>(Definiteness::General);
    let none = {
        let b = random_vector::<f64>(a.nrows(), 2);
        let mut x = vec![0.0; a.nrows()];
        Gmres::new(GmresOptions {
            tol: Tolerances { rtol: 1e-8, max_it: 400, ..Default::default() },
            restart: 400,
            ..Default::default()
        })
        .solve(&a, &Identity, &b, &mut x)
        .iterations
    };
    let mut its = Vec::new();
    for f in [SchurFact::Diag, SchurFact::Upper, SchurFact::Full] {
        let pc = SchurFieldSplit::new(&a, i0.clone(), i1.clone(), f, SchurPre::SelfP, false, &direct, &direct).unwrap();
        its.push(gmres_its(&a, &pc, 1e-8));
    }
    eprintln!("selfp Schur: no PC {none}, diag/upper/full {its:?}");
    assert!(its[2] <= its[1] && its[1] <= its[0] + 1, "{its:?}");
    assert!(its[0] * 3 < none || none >= 400, "{its:?} vs {none}");
}

#[test]
fn minres_with_the_positive_definite_block_diagonal_schur_preconditioner() {
    let (a, i0, i1) = saddle(12, 40, 0.0);
    let s = exact_schur(&a, &i0, &i1);
    let direct = direct_maker::<f64>(Definiteness::General);
    // diag(A00⁻¹, −S⁻¹) is SPD for A11 = 0
    let pc = SchurFieldSplit::new(&a, i0, i1, SchurFact::Diag, SchurPre::User(s), true, &direct, &direct).unwrap();
    let b = random_vector::<f64>(a.nrows(), 4);
    let mut x = vec![0.0; a.nrows()];
    let info = Minres::new(MinresOptions {
        tol: Tolerances { rtol: 1e-9, max_it: 100, ..Default::default() },
        ..Default::default()
    })
    .solve(&a, &pc, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    assert!(info.iterations <= 4, "{}", info.iterations); // 3 distinct eigenvalues
}

#[test]
fn multiplicative_field_split_beats_additive_on_a_coupled_three_field_problem() {
    // three coupled fields on the same 2-D grid: A = I₃ ⊗ (L + I) + C ⊗ I with a non-symmetric coupling C (‖C‖ ≈ 0.9)
    let nx = 20;
    let l = laplacian_2d(nx, nx);
    let n0 = nx * nx;
    let c = [[0.0, 0.8, 0.0], [0.5, 0.0, 0.6], [0.0, 0.7, 0.0]];
    let mut t = Vec::new();
    for f in 0..3 {
        for i in 0..n0 {
            let (cs, vs) = l.row(i);
            for (&cc, &v) in cs.iter().zip(vs) {
                // diagonal blocks L + I keep the coupled matrix nonsingular for the strong coupling C below
                t.push((f * n0 + i, f * n0 + cc as usize, if cc as usize == i { v + 1.0 } else { v }));
            }
            for g in 0..3 {
                if c[f][g] != 0.0 {
                    t.push((f * n0 + i, g * n0 + i, c[f][g]));
                }
            }
        }
    }
    let a = Csr::from_triplets(3 * n0, 3 * n0, &t).unwrap();
    let fields: Vec<Vec<usize>> = (0..3).map(|f| (f * n0..(f + 1) * n0).collect()).collect();
    let direct = direct_maker::<f64>(Definiteness::General);
    let its = |kind| {
        let pc = FieldSplit::new(&a, fields.clone(), &direct, kind).unwrap();
        gmres_its(&a, &pc, 1e-8)
    };
    let (add, mult, sym) =
        (its(SplitKind::Additive), its(SplitKind::Multiplicative), its(SplitKind::SymmetricMultiplicative));
    eprintln!(
        "3-field problem, GMRES iterations: additive {add}, multiplicative {mult}, symmetric multiplicative {sym}"
    );
    assert!(mult < add && sym <= mult + 1, "{add} {mult} {sym}");
    // invalid splits are rejected
    assert!(FieldSplit::new(&a, vec![fields[0].clone(), fields[1].clone()], &direct, SplitKind::Additive).is_err());
    assert!(
        FieldSplit::new(
            &a,
            vec![fields[0].clone(), fields[0].clone(), fields[2].clone()],
            &direct,
            SplitKind::Additive
        )
        .is_err()
    );
}

#[test]
fn complex_symmetric_piezo_like_coupling() {
    // [[K − ω²M + iωD, G], [Gᵀ, −Kφ]] with complex entries: Schur "selfp" + exact blocks, GMRES
    let nx = 12;
    let n0 = nx * nx;
    let n1 = 40;
    let k = laplacian_2d(nx, nx).to_complex();
    let shift = Csr::axpby(C64::ONE, &k, C64::new(-0.3, 0.2), &Csr::identity(n0));
    let mut t = Vec::new();
    for i in 0..n0 {
        let (cs, vs) = shift.row(i);
        for (&c, &v) in cs.iter().zip(vs) {
            t.push((i, c as usize, v));
        }
    }
    for r in 0..n1 {
        for (c, w) in [(3 * r, C64::new(1.0, 0.1)), (3 * r + 1, C64::new(-0.7, 0.0)), (3 * r + 2, C64::new(0.3, -0.2))]
        {
            t.push((n0 + r, c, w));
            t.push((c, n0 + r, w));
        }
        t.push((n0 + r, n0 + r, C64::new(-1.0, 0.0)));
    }
    let a = Csr::from_triplets(n0 + n1, n0 + n1, &t).unwrap();
    let (i0, i1): (Vec<usize>, Vec<usize>) = ((0..n0).collect(), (n0..n0 + n1).collect());
    let direct = direct_maker::<C64>(Definiteness::General);
    let pc = SchurFieldSplit::new(&a, i0, i1, SchurFact::Upper, SchurPre::SelfP, false, &direct, &direct).unwrap();
    let b = random_vector::<C64>(n0 + n1, 8);
    let mut x = vec![C64::ZERO; n0 + n1];
    let info = Gmres::new(GmresOptions {
        tol: Tolerances { rtol: 1e-9, max_it: 300, ..Default::default() },
        restart: 300,
        ..Default::default()
    })
    .solve(&a, &pc, &b, &mut x);
    assert!(info.reason.is_converged(), "{:?} after {}", info.reason, info.iterations);
    let mut ax = vec![C64::ZERO; n0 + n1];
    a.spmv(&x, &mut ax);
    assert!(ax.iter().zip(&b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() < 1e-6 * norm(&b));
    let _ = Scalar::modulus(1.0f64);
}
