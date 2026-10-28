//! Locally optimal block preconditioned conjugate gradients (Knyazev 2001) for the lowest
//! eigenpairs of `K x = lambda M x` with `K` symmetric positive (semi)definite and `M` symmetric
//! positive definite.
//!
//! Unlike shift-and-invert Krylov methods it needs no factorisation and no accurate linear solves:
//! every iteration applies the preconditioner (one multigrid cycle for `K`) once per vector of the
//! block, so it is the method for problems whose factorisation does not fit in memory. Blocks are
//! stored row-major (all columns of a row together) so the matrix-block products and the dense
//! Gram matrices run on all cores in blocks of rows. The block is orthonormalised in the `M`
//! inner product by the SVQB method (eigendecomposition of the Gram matrix, dropping directions
//! that have collapsed), which is robust when residuals become small.

use crate::eigen::EigPairs;
use crate::linalg::eigh;
use crate::sparse::Csr;
use rayon::prelude::*;

fn threads() -> usize {
    crate::par::threads()
}

/// Blocks of work below this many multiply-adds run on the calling thread (a pool task costs ~10 microseconds).
const PAR_MIN_WORK: usize = 60_000;

/// Rows per piece of a Gram matrix sum (see `gram`).
const GRAM_ROWS: usize = 1024;

/// `Y = A X` for a row-major block `X` of `s` columns.
fn spmm_block(a: &Csr, x: &[f64], s: usize) -> Vec<f64> {
    let n = a.nrows;
    let mut y = vec![0.0; n * s];
    let nt = if a.nnz() * s < PAR_MIN_WORK { 1 } else { threads() };
    let chunk = n.div_ceil(nt).max(1);
    y.par_chunks_mut(chunk * s).enumerate().for_each(|(ci, ys)| {
        for (i, row) in ys.chunks_mut(s).enumerate() {
            let r = ci * chunk + i;
            for k in a.indptr[r]..a.indptr[r + 1] {
                let (c, v) = (a.indices[k], a.data[k]);
                let xr = &x[c * s..(c + 1) * s];
                for j in 0..s {
                    row[j] += v * xr[j];
                }
            }
        }
    });
    y
}

/// `G = A^T B` (`sa x sb`, row-major) for row-major blocks with `n` rows.
fn gram(a: &[f64], sa: usize, b: &[f64], sb: usize, n: usize) -> Vec<f64> {
    // The rows are cut into fixed pieces of GRAM_ROWS rows whatever the number of threads and the piece sums are
    // added in order, so the result does not depend on how many threads run it (a marginal LOBPCG iteration with a
    // weak preconditioner converged or not depending on the rounding of these sums).
    let chunk = GRAM_ROWS;
    let piece = |t: usize| {
        let mut g = vec![0.0; sa * sb];
        for r in t * chunk..((t + 1) * chunk).min(n) {
            let (ar, br) = (&a[r * sa..(r + 1) * sa], &b[r * sb..(r + 1) * sb]);
            for i in 0..sa {
                let ai = ar[i];
                for j in 0..sb {
                    g[i * sb + j] += ai * br[j];
                }
            }
        }
        g
    };
    let pieces = n.div_ceil(chunk);
    let parts: Vec<Vec<f64>> = if n * sa * sb < PAR_MIN_WORK { (0..pieces).map(piece).collect() } else { (0..pieces).into_par_iter().map(piece).collect() };
    let mut g = vec![0.0; sa * sb];
    for p in parts {
        g.iter_mut().zip(&p).for_each(|(x, y)| *x += y);
    }
    g
}

/// `Y = X C` with `X` (`n x s`) and `C` (`s x c`), row-major.
fn combine(x: &[f64], s: usize, c: &[f64], cc: usize, n: usize) -> Vec<f64> {
    let mut y = vec![0.0; n * cc];
    let nt = if n * s * cc < PAR_MIN_WORK { 1 } else { threads() };
    let chunk = n.div_ceil(nt).max(1);
    y.par_chunks_mut(chunk * cc).enumerate().for_each(|(ci, ys)| {
        for (i, row) in ys.chunks_mut(cc).enumerate() {
            let r = ci * chunk + i;
            let xr = &x[r * s..(r + 1) * s];
            for k in 0..s {
                let xv = xr[k];
                if xv != 0.0 {
                    let crow = &c[k * cc..(k + 1) * cc];
                    for j in 0..cc {
                        row[j] += xv * crow[j];
                    }
                }
            }
        }
    });
    y
}

/// Columns of two row-major blocks side by side.
fn hcat(blocks: &[(&[f64], usize)], n: usize) -> Vec<f64> {
    let s: usize = blocks.iter().map(|b| b.1).sum();
    let mut out = vec![0.0; n * s];
    for r in 0..n {
        let mut off = 0;
        for (b, w) in blocks {
            out[r * s + off..r * s + off + w].copy_from_slice(&b[r * w..(r + 1) * w]);
            off += w;
        }
    }
    out
}

/// Scales the columns of `v`, `kv` and `mv` (blocks `n x s`, with `mv = M v`) so that every column
/// has unit `M`-norm; zero columns stay zero. Keeps the Gram matrix well scaled when residuals
/// are many orders of magnitude smaller than the vectors themselves.
fn normalize_columns(v: &mut [f64], kv: &mut [f64], mv: &mut [f64], s: usize, n: usize) {
    let g = gram(v, s, mv, s, n);
    for j in 0..s {
        let nrm = g[j * s + j].max(0.0).sqrt();
        if nrm > 0.0 {
            let inv = 1.0 / nrm;
            for r in 0..n {
                v[r * s + j] *= inv;
                kv[r * s + j] *= inv;
                mv[r * s + j] *= inv;
            }
        }
    }
}

fn column(x: &[f64], s: usize, j: usize, n: usize) -> Vec<f64> {
    (0..n).map(|r| x[r * s + j]).collect()
}

/// Symmetric dense matrix (row-major, `n x n`) eigendecomposition: ascending eigenvalues and the
/// eigenvectors as the columns of a row-major matrix.
fn eig_dense(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    let mut sym = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            sym[i * n + j] = 0.5 * (a[i * n + j] + a[j * n + i]);
        }
    }
    let (vals, vecs) = eigh(&sym, n)?;
    let mut u = vec![0.0; n * n];
    for (j, v) in vecs.iter().enumerate() {
        for i in 0..n {
            u[i * n + j] = v[i];
        }
    }
    Ok((vals, u))
}

/// Start of the error text returned when the stall guard gives up (see `lobpcg_guarded`).
pub const STALLED: &str = "LOBPCG stalled";

/// Iterations over which the stall guard compares the best residual with the best one before.
const STALL_WINDOW: usize = 10;

/// `Some(best recent residual)` when the history `worst` of the worst residual per iteration has stopped improving: from
/// iteration `2 * STALL_WINDOW` on, the best value of the last `STALL_WINDOW` iterations is not below 90 % of the best value
/// of the `STALL_WINDOW` before.
fn stalled(worst: &[f64]) -> Option<f64> {
    let n = worst.len();
    if n < 2 * STALL_WINDOW + 1 {
        return None;
    }
    let best = |s: &[f64]| s.iter().cloned().fold(f64::INFINITY, f64::min);
    let recent = best(&worst[n - STALL_WINDOW..]);
    let before = best(&worst[n - 2 * STALL_WINDOW..n - STALL_WINDOW]);
    (recent > 0.9 * before).then_some(recent)
}

/// The `nev` lowest eigenpairs of `K x = lambda M x`. `precond(r)` applies an approximation of
/// `K^-1` (for example one multigrid cycle) to a vector. `extra` additional vectors are carried
/// to speed up convergence of the wanted ones. Eigenvectors are `M`-orthonormal.
pub fn lobpcg(k: &Csr, m: &Csr, nev: usize, extra: usize, precond: &dyn Fn(&[f64]) -> Vec<f64>, tol: f64, max_iter: usize) -> Result<EigPairs, String> {
    lobpcg_guarded(k, m, nev, extra, precond, tol, max_iter, false)
}

/// Like `lobpcg`; with `stall_guard` the iteration gives up with an error starting with `STALLED` as soon as
/// the worst residual of the wanted vectors has stopped decreasing (the best residual of the last
/// `STALL_WINDOW` iterations is not even 10 % below the best one of the `STALL_WINDOW` before, from iteration
/// `2 * STALL_WINDOW` on). Thin or stretched parts can make the multigrid preconditioner good enough for the
/// static problem but not for the eigenvectors: the residual then sits at 1e-3..1e-2 for hundreds of iterations
/// (a 70,000-unknown slender beam burnt 130-240 s before the caller fell back to the factorisation, which takes 11 s).
#[allow(clippy::too_many_arguments)]
pub fn lobpcg_guarded(k: &Csr, m: &Csr, nev: usize, extra: usize, precond: &dyn Fn(&[f64]) -> Vec<f64>, tol: f64, max_iter: usize, stall_guard: bool) -> Result<EigPairs, String> {
    let n = k.nrows;
    if nev == 0 || nev + extra + 1 >= n / 3 {
        return Err("number of eigenvalues must be small compared with the number of unknowns".into());
    }
    let _stage = crate::progress::stage(format!("LOBPCG eigensolver with multigrid preconditioner ({nev} modes, {n} unknowns)"));
    let p = nev + extra;
    // random start, M-orthonormalised, then Rayleigh-Ritz
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        ((seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    };
    let x0: Vec<f64> = (0..n * p).map(|_| rnd()).collect();
    let mx0 = spmm_block(m, &x0, p);
    let g0 = gram(&x0, p, &mx0, p, n);
    let (d0, u0) = eig_dense(&g0, p)?;
    let dmax = d0.iter().cloned().fold(0.0, f64::max);
    let mut q0 = vec![0.0; p * p];
    for j in 0..p {
        let sc = if d0[j] > 1e-12 * dmax { 1.0 / d0[j].sqrt() } else { 0.0 };
        for i in 0..p {
            q0[i * p + j] = u0[i * p + j] * sc;
        }
    }
    let mut x = combine(&x0, p, &q0, p, n);
    let mut kx = spmm_block(k, &x, p);
    let mut mxb = spmm_block(m, &x, p);
    {
        let a = gram(&x, p, &kx, p, n);
        let (th, y) = eig_dense(&a, p)?;
        x = combine(&x, p, &y, p, n);
        kx = combine(&kx, p, &y, p, n);
        mxb = combine(&mxb, p, &y, p, n);
        let _ = th;
    }
    let mut lam: Vec<f64> = vec![0.0; p];
    // search directions
    let mut pdir: Option<(Vec<f64>, Vec<f64>, Vec<f64>)> = None;
    let mut converged = false;
    let mut iters = 0;
    let mut worst: Vec<f64> = vec![];
    for it in 0..max_iter {
        iters = it + 1;
        // Rayleigh quotients and residuals
        let xkx = gram(&x, p, &kx, p, n);
        for j in 0..p {
            lam[j] = xkx[j * p + j];
        }
        let mut r = vec![0.0; n * p];
        for row in 0..n {
            for j in 0..p {
                r[row * p + j] = kx[row * p + j] - lam[j] * mxb[row * p + j];
            }
        }
        let rr = gram(&r, p, &r, p, n);
        let kk = gram(&kx, p, &kx, p, n);
        let rel: Vec<f64> = (0..p).map(|j| (rr[j * p + j] / kk[j * p + j].max(1e-300)).sqrt()).collect();
        if std::env::var_os("NYX_LOBPCG_TRACE").is_some() && it % 10 == 0 {
            eprintln!("  lobpcg it {it}: lam0 {:.6e} rel {:?}", lam[0], &rel[..nev.min(3)]);
        }
        crate::progress::check()?;
        if it % 10 == 0 {
            let worst_now = (0..nev).map(|j| rel[j]).fold(0.0, f64::max);
            crate::progress::note(format!("iteration {it}: worst residual of the {nev} wanted modes {worst_now:.1e} (goal {tol:.0e})"));
            // how far the residual has come from 1 down to the goal, on a log scale
            crate::progress::fraction(worst_now.max(1e-300).ln() / tol.ln());
        }
        if (0..nev).all(|j| rel[j] < tol) {
            converged = true;
            break;
        }
        if stall_guard {
            worst.push((0..nev).map(|j| rel[j]).fold(0.0, f64::max));
            if let Some(recent) = stalled(&worst) {
                return Err(format!("{STALLED} after {} iterations (residual {recent:.1e})", it + 1));
            }
        }
        // preconditioned residuals (converged vectors are not updated)
        let mut w = vec![0.0; n * p];
        for j in 0..p {
            if rel[j] < 0.1 * tol {
                continue;
            }
            let t = precond(&column(&r, p, j, n));
            for row in 0..n {
                w[row * p + j] = t[row];
            }
        }
        let mut kw = spmm_block(k, &w, p);
        let mut mw = spmm_block(m, &w, p);
        normalize_columns(&mut w, &mut kw, &mut mw, p, n);
        let (s, ks, ms, sw) = match &pdir {
            Some((pp, kp, mp)) => (hcat(&[(&x, p), (&w, p), (pp, p)], n), hcat(&[(&kx, p), (&kw, p), (kp, p)], n), hcat(&[(&mxb, p), (&mw, p), (mp, p)], n), 3 * p),
            None => (hcat(&[(&x, p), (&w, p)], n), hcat(&[(&kx, p), (&kw, p)], n), hcat(&[(&mxb, p), (&mw, p)], n), 2 * p),
        };
        let g = gram(&s, sw, &ms, sw, n);
        let a = gram(&s, sw, &ks, sw, n);
        // SVQB: M-orthonormal basis of the span, dropping collapsed directions
        let (dg, ug) = eig_dense(&g, sw)?;
        let dmax = dg.iter().cloned().fold(0.0, f64::max);
        let keep: Vec<usize> = (0..sw).filter(|&j| dg[j] > 1e-10 * dmax).collect();
        let rk = keep.len();
        if rk < p {
            return Err("the block of vectors collapsed (is the mass matrix positive definite?)".into());
        }
        let mut qm = vec![0.0; sw * rk];
        for (c, &j) in keep.iter().enumerate() {
            let sc = 1.0 / dg[j].sqrt();
            for i in 0..sw {
                qm[i * rk + c] = ug[i * sw + j] * sc;
            }
        }
        // projected matrix Q^T A Q
        let mut aq = vec![0.0; sw * rk];
        for i in 0..sw {
            for c in 0..rk {
                let mut acc = 0.0;
                for l in 0..sw {
                    acc += a[i * sw + l] * qm[l * rk + c];
                }
                aq[i * rk + c] = acc;
            }
        }
        let mut ap = vec![0.0; rk * rk];
        for i in 0..rk {
            for c in 0..rk {
                let mut acc = 0.0;
                for l in 0..sw {
                    acc += qm[l * rk + i] * aq[l * rk + c];
                }
                ap[i * rk + c] = acc;
            }
        }
        let (_th, y) = eig_dense(&ap, rk)?;
        // coefficients of the p lowest Ritz vectors in the original basis: C = Q Y[:, :p]
        let mut cmat = vec![0.0; sw * p];
        for i in 0..sw {
            for j in 0..p {
                let mut acc = 0.0;
                for l in 0..rk {
                    acc += qm[i * rk + l] * y[l * rk + j];
                }
                cmat[i * p + j] = acc;
            }
        }
        let xn = combine(&s, sw, &cmat, p, n);
        let kxn = combine(&ks, sw, &cmat, p, n);
        let mxn = combine(&ms, sw, &cmat, p, n);
        // new directions: the part of the update that does not come from the old X
        let mut cp = cmat.clone();
        for i in 0..p {
            for j in 0..p {
                cp[i * p + j] = 0.0;
            }
        }
        let (mut pn, mut kpn, mut mpn) = (combine(&s, sw, &cp, p, n), combine(&ks, sw, &cp, p, n), combine(&ms, sw, &cp, p, n));
        normalize_columns(&mut pn, &mut kpn, &mut mpn, p, n);
        pdir = Some((pn, kpn, mpn));
        x = xn;
        kx = kxn;
        mxb = mxn;
        if (it + 1) % 8 == 0 {
            // the recurrences for K x and M x drift by rounding: recompute them now and then
            kx = spmm_block(k, &x, p);
            mxb = spmm_block(m, &x, p);
        }
    }
    let xkx = gram(&x, p, &kx, p, n);
    let mut out = EigPairs { values: vec![], vectors: vec![], dim: p, converged };
    for j in 0..nev {
        out.values.push(xkx[j * p + j]);
        out.vectors.push(column(&x, p, j, n));
    }
    let _ = iters;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Coo;

    #[test]
    fn the_stall_guard_tells_flat_histories_from_converging_ones() {
        // converging by a factor 0.8 per iteration: never stalled
        let good: Vec<f64> = (0..80).map(|i| 0.8f64.powi(i)).collect();
        for n in 1..=good.len() {
            assert!(stalled(&good[..n]).is_none(), "{n}");
        }
        // fast at first, then flat around 5e-3 (what a slender 3D beam does): detected shortly after it flattens
        let mut bad: Vec<f64> = (0..12).map(|i| 0.5f64.powi(i)).collect();
        bad.extend((0..60).map(|i| 4e-4 + 1e-4 * ((i * 7 % 5) as f64) / 5.0));
        let first = (1..=bad.len()).find(|&n| stalled(&bad[..n]).is_some()).expect("flat history is detected");
        assert!(first <= 45, "detected only after {first} iterations");
        // too short a history is never judged
        assert!(stalled(&[1.0; 20]).is_none());
    }

    /// 1D Laplacian with Dirichlet ends and identity mass: eigenvalues 2 - 2 cos(k pi / (n + 1)).
    #[test]
    fn finds_the_lowest_eigenvalues_of_a_chain() {
        let n = 400;
        let mut kc = Coo::new(n, n);
        let mut mc = Coo::new(n, n);
        for i in 0..n {
            kc.push(i, i, 2.0);
            if i + 1 < n {
                kc.push(i, i + 1, -1.0);
                kc.push(i + 1, i, -1.0);
            }
            mc.push(i, i, 1.0);
        }
        let (k, m) = (kc.to_csr(), mc.to_csr());
        // a deliberately weak preconditioner: the inverse diagonal
        let pre = |r: &[f64]| r.iter().map(|x| x / 2.0).collect::<Vec<f64>>();
        let res = lobpcg(&k, &m, 4, 4, &pre, 1e-6, 4000).unwrap();
        assert!(res.converged);
        for (j, &lam) in res.values.iter().enumerate() {
            let exact = 2.0 - 2.0 * ((j + 1) as f64 * std::f64::consts::PI / (n + 1) as f64).cos();
            assert!((lam - exact).abs() < 1e-8 * exact, "{j}: {lam} vs {exact}");
        }
    }
}
