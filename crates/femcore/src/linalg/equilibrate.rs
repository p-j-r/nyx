// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Symmetric Ruiz diagonal equilibration (`D A D`) for multi-scale and multi-physics
//! sparse matrices (real `RCsr` and complex `CCsr`), preserving symmetry (`A = A^T` and `A = A^H`).

use crate::linalg::csr::{CCsr, Csr, RCsr};
use crate::linalg::par;
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::CHUNK;

/// Default number of symmetric Ruiz scaling sweeps.
pub const DEFAULT_RUIZ_PASSES: usize = 8;

/// A weighted sparse matrix term contributing to the joint infinity-norm envelope in
/// [`ruiz_symmetric_terms`], allowing real and complex affine terms to be equilibrated together.
#[derive(Clone, Copy)]
pub enum RuizTerm<'a> {
    Real(&'a RCsr, f64),
    Complex(&'a CCsr, f64),
}

/// Computes the symmetric Ruiz diagonal scaling vector `d` (`n`) over `passes` sweeps for
/// a list of weighted sparse matrix terms `sum_k w_k A_k`, clamped to `[1e-14, 1e14]`.
///
/// The scaled matrices `D A_k D` (with `D = diag(d)`) preserve structural and numerical symmetry
/// (`A = A^T` or `A = A^H`) while reducing row/column dynamic range across disparate physical blocks.
pub fn ruiz_symmetric_terms(n: usize, terms: &[RuizTerm<'_>], passes: usize) -> Vec<f64> {
    let mut d = vec![1.0_f64; n];
    if n == 0 || terms.is_empty() || passes == 0 {
        return d;
    }
    let total_nnz: usize = terms
        .iter()
        .map(|t| match *t {
            RuizTerm::Real(m, w) if w != 0.0 => m.nnz(),
            RuizTerm::Complex(m, w) if w != 0.0 => m.nnz(),
            _ => 0,
        })
        .sum();
    if total_nnz == 0 {
        return d;
    }
    let parallel = par::worth_parallel(total_nnz);
    let mut rmax = vec![0.0_f64; n];

    for _ in 0..passes {
        par::for_chunks_mut(&mut rmax, CHUNK, parallel, |chunk_idx, rc| {
            let r0 = chunk_idx * CHUNK;
            for (off, slot) in rc.iter_mut().enumerate() {
                let i = r0 + off;
                let di = d[i];
                let mut m_val = 0.0_f64;
                for term in terms {
                    match *term {
                        RuizTerm::Real(mat, w) if mat.nnz() > 0 && w != 0.0 => {
                            let w_abs = w.abs();
                            let (cols, vals) = mat.row(i);
                            for (&j, &v) in cols.iter().zip(vals) {
                                m_val = m_val.max(di * w_abs * v.abs() * d[j as usize]);
                            }
                        }
                        RuizTerm::Complex(mat, w) if mat.nnz() > 0 && w != 0.0 => {
                            let w_abs = w.abs();
                            let (cols, vals) = mat.row(i);
                            for (&j, &v) in cols.iter().zip(vals) {
                                m_val = m_val.max(di * w_abs * v.norm() * d[j as usize]);
                            }
                        }
                        _ => {}
                    }
                }
                *slot = m_val;
            }
        });
        for i in 0..n {
            if rmax[i] > 1e-300 {
                d[i] = (d[i] / rmax[i].sqrt()).clamp(1e-14, 1e14);
            }
        }
    }
    d
}

/// Computes `DEFAULT_RUIZ_PASSES` symmetric Ruiz diagonal scaling weights for a single CSR matrix.
pub fn ruiz_symmetric<T: Scalar>(a: &Csr<T>) -> Vec<f64> {
    ruiz_symmetric_weighted(&[(a, 1.0)], a.nrows(), DEFAULT_RUIZ_PASSES)
}

/// Computes symmetric Ruiz diagonal scaling weights for weighted homogeneous `Csr<T>` terms.
pub fn ruiz_symmetric_weighted<T: Scalar>(terms: &[(&Csr<T>, f64)], n: usize, passes: usize) -> Vec<f64> {
    let mut d = vec![1.0_f64; n];
    if n == 0 || terms.is_empty() || passes == 0 {
        return d;
    }
    let total_nnz: usize = terms.iter().filter(|(_, w)| *w != 0.0).map(|(m, _)| m.nnz()).sum();
    if total_nnz == 0 {
        return d;
    }
    let parallel = par::worth_parallel(total_nnz);
    let mut rmax = vec![0.0_f64; n];
    for _ in 0..passes {
        par::for_chunks_mut(&mut rmax, CHUNK, parallel, |chunk_idx, rc| {
            let r0 = chunk_idx * CHUNK;
            for (off, slot) in rc.iter_mut().enumerate() {
                let i = r0 + off;
                let di = d[i];
                let mut m_val = 0.0_f64;
                for &(mat, w) in terms {
                    if mat.nnz() == 0 || w == 0.0 {
                        continue;
                    }
                    let w_abs = w.abs();
                    let (cols, vals) = mat.row(i);
                    for (&j, &v) in cols.iter().zip(vals) {
                        m_val = m_val.max(di * w_abs * v.modulus() * d[j as usize]);
                    }
                }
                *slot = m_val;
            }
        });
        for i in 0..n {
            if rmax[i] > 1e-300 {
                d[i] = (d[i] / rmax[i].sqrt()).clamp(1e-14, 1e14);
            }
        }
    }
    d
}

/// Forms `D A D` symmetrically on the identical sparsity pattern of `a`.
pub fn scale_symmetric<T: Scalar>(a: &Csr<T>, d: &[f64]) -> Csr<T> {
    let n = a.nrows();
    assert_eq!(a.ncols(), n, "scale_symmetric requires a square matrix");
    assert_eq!(d.len(), n, "scale_symmetric diagonal length mismatch");
    if a.nnz() == 0 {
        return Csr::zeros(n, n);
    }
    let mut vals = a.values().to_vec();
    for i in 0..n {
        let di = d[i];
        let range = a.row_range(i);
        let cols = &a.indices()[range.clone()];
        let row_vals = &mut vals[range];
        for (v, &col) in row_vals.iter_mut().zip(cols) {
            *v = v.scale(di * d[col as usize]);
        }
    }
    Csr::try_new(n, n, a.indptr().to_vec(), a.indices().to_vec(), vals).expect("identical CSR pattern")
}

/// Scales a vector component-wise: `out[i] = v[i] * d[i]`.
pub fn scale_vec<T: Scalar>(v: &[T], d: &[f64]) -> Vec<T> {
    assert_eq!(v.len(), d.len());
    v.iter().zip(d).map(|(&vi, &di)| vi.scale(di)).collect()
}

/// Unscales a vector component-wise: `out[i] = v[i] / d[i]`.
pub fn unscale_vec<T: Scalar>(v: &[T], d: &[f64]) -> Vec<T> {
    assert_eq!(v.len(), d.len());
    v.iter().zip(d).map(|(&vi, &di)| vi.scale(1.0 / di)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::csr::c;

    #[test]
    fn ruiz_equilibrates_22_order_saddle_point_and_preserves_symmetry() {
        // 2x2 block system with 1e11 vs 1e-11 contrast
        let a = Csr::from_triplets(
            2,
            2,
            &[(0, 0, 1.0e11), (0, 1, 1.5), (1, 0, 1.5), (1, 1, -1.0e-11)],
        )
        .unwrap();
        let d = ruiz_symmetric(&a);
        let a_hat = scale_symmetric(&a, &d);
        assert!(a_hat.is_symmetric(1e-14));
        for i in 0..2 {
            let rmax = a_hat.row(i).1.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
            assert!(rmax > 0.85 && rmax < 1.05, "row {i} max after 8 passes = {rmax}");
        }
        let d16 = ruiz_symmetric_weighted(&[(&a, 1.0)], 2, 16);
        let a_hat16 = scale_symmetric(&a, &d16);
        for i in 0..2 {
            let rmax = a_hat16.row(i).1.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
            assert!((rmax - 1.0).abs() < 1e-3, "row {i} max after 16 passes = {rmax}");
        }
        // Complex variant
        let ac = CCsr::from_triplets(
            2,
            2,
            &[
                (0, 0, c(1.0e11, 2.0e8)),
                (0, 1, c(1.5, 0.0)),
                (1, 0, c(1.5, 0.0)),
                (1, 1, c(-1.0e-11, 1.0e-14)),
            ],
        )
        .unwrap();
        let dc = ruiz_symmetric_terms(2, &[RuizTerm::Complex(&ac, 1.0)], DEFAULT_RUIZ_PASSES);
        let ac_hat = scale_symmetric(&ac, &dc);
        assert!(ac_hat.is_symmetric(1e-14));
        let r0 = ac_hat.row(0).1.iter().fold(0.0_f64, |acc, v| acc.max(v.norm()));
        let r1 = ac_hat.row(1).1.iter().fold(0.0_f64, |acc, v| acc.max(v.norm()));
        assert!(r0 > 0.85 && r0 < 1.05 && r1 > 0.85 && r1 < 1.05, "r0={r0}, r1={r1}");
    }
}
