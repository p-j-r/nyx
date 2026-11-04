//! Real block Arnoldi with shift-and-invert for the real, non-symmetric pencil of the 2D waveguide modes (`ems`): the eigenvalues of
//! `A x = lambda B x` nearest a real shift, real vectors, complex Ritz values skipped. The library's Krylov-Schur solves a real
//! non-symmetric pencil only in complex arithmetic, where the real eigenvalues of this non-normal pencil pick up imaginary parts of
//! 1e-8 to 1e-5 relative (finding of phase 3); here the arithmetic stays real, so real eigenvalues stay real. Ported unchanged from
//! the old `eigen.rs` onto the one matrix type and `driver::DirectSolver`.

use crate::linalg::csr::RCsr as Csr;
use crate::linalg::driver::DirectSolver;
use rayon::prelude::*;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}


#[derive(Clone, Debug)]
pub struct EigOpts {
    /// relative residual `|OP x - nu x| / |nu x|` required of every returned pair
    pub tol: f64,
    /// cap on the Krylov basis size (memory is about 3 n dim doubles)
    pub max_dim: usize,
    /// block size; 0 picks `min(nev, 8)`
    pub block: usize,
}

impl Default for EigOpts {
    fn default() -> Self {
        EigOpts { tol: 1e-10, max_dim: 0, block: 0 }
    }
}

#[derive(Clone, Debug)]
pub struct EigPairs {
    /// eigenvalues of the original problem, closest to the shift first
    pub values: Vec<f64>,
    /// matching eigenvectors (unit norm: M-norm for the symmetric solver, 2-norm otherwise)
    pub vectors: Vec<Vec<f64>>,
    /// size of the Krylov basis that was used
    pub dim: usize,
    /// whether every requested pair met the tolerance
    pub converged: bool,
    /// a *converged complex pair* (conjugate Ritz values: a leaky or complex-propagation-constant mode) is among the wanted eigenvalues.
    /// The real solver then stops and returns what it has; the caller must solve the pencil in complex arithmetic. Never set by the
    /// symmetric solvers.
    pub complex_found: bool,
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }
}

fn max_dim_for(n: usize, nev: usize, opts: &EigOpts) -> usize {
    let cap = if opts.max_dim > 0 { opts.max_dim } else { (60_000_000 / n.max(1)).clamp(4 * nev + 24, 600) };
    cap.min(n)
}

fn par_threads(work: usize) -> usize {
    // a task on the shared pool costs about 10 microseconds, a few thousand multiply-adds
    if work < 100_000 { 1 } else { crate::par::threads() }
}

/// `<q_k, w>` for every vector `q_k`, one thread per group of vectors.
fn par_dots(vs: &[Vec<f64>], w: &[f64]) -> Vec<f64> {
    let m = vs.len();
    let mut out = vec![0.0; m];
    let nt = par_threads(m * w.len());
    let chunk = m.div_ceil(nt).max(1);
    out.par_chunks_mut(chunk).enumerate().for_each(|(ci, o)| {
        for (k, x) in o.iter_mut().enumerate() {
            *x = dot(&vs[ci * chunk + k], w);
        }
    });
    out
}

/// `sum_k y_k q_k` with the rows split over threads.
fn par_combine(vs: &[Vec<f64>], y: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    let nt = par_threads(vs.len() * n);
    let chunk = n.div_ceil(nt).max(1);
    out.par_chunks_mut(chunk).enumerate().for_each(|(ci, o)| {
        let r0 = ci * chunk;
        let len = o.len();
        for (q, &yk) in vs.iter().zip(y) {
            if yk != 0.0 {
                for (x, a) in o.iter_mut().zip(&q[r0..r0 + len]) {
                    *x += yk * a;
                }
            }
        }
    });
    out
}

/// `w -= sum_k c_k q_k`, rows split over threads.
fn par_subtract(w: &mut [f64], vs: &[Vec<f64>], coef: &[f64]) {
    let nt = par_threads(vs.len() * w.len());
    let chunk = w.len().div_ceil(nt).max(1);
    w.par_chunks_mut(chunk).enumerate().for_each(|(ci, o)| {
        let r0 = ci * chunk;
        let len = o.len();
        for (q, &ck) in vs.iter().zip(coef) {
            for (x, a) in o.iter_mut().zip(&q[r0..r0 + len]) {
                *x -= ck * a;
            }
        }
    });
}

/// General real pencil `A x = lambda B x` with real eigenvalues wanted (A - sigma B must be
/// nonsingular), nearest to the real shift `sigma`. Eigenpairs whose Ritz value is complex are
/// skipped. Used for the (indefinite) mode-analysis pencil.
pub fn eigs_general(a: &Csr, b: &Csr, sigma: f64, nev: usize, opts: &EigOpts) -> Result<EigPairs, String> {
    let n = a.nrows();
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    let shifted = a.add_scaled(b, -sigma);
    let solver = DirectSolver::new(&shifted, false)?;
    let p = if opts.block > 0 { opts.block } else { nev.min(8) };
    let max_dim = max_dim_for(n, nev, opts);

    let mut v: Vec<Vec<f64>> = Vec::new();
    let mut ov: Vec<Vec<f64>> = Vec::new();
    let mut h: Vec<Vec<f64>> = Vec::new(); // h[i][j] = v_i . (OP v_j)
    let mut rng = Rng(0xD1B54A32D192ED03);
    let mut cand: Vec<Vec<f64>> = (0..p).map(|_| (0..n).map(|_| rng.next()).collect()).collect();
    loop {
        crate::progress::check()?;
        let before = v.len();
        for mut w in cand.drain(..) {
            let w0 = dot(&w, &w).sqrt();
            for _pass in 0..2 {
                let coef = par_dots(&v, &w);
                par_subtract(&mut w, &v, &coef);
            }
            let nrm = dot(&w, &w).sqrt();
            if nrm > 1e-10 * w0 && nrm > 0.0 && v.len() < max_dim {
                w.iter_mut().for_each(|x| *x /= nrm);
                v.push(w);
            }
        }
        let e = v.len();
        if e == before {
            break;
        }
        let rhs: Vec<Vec<f64>> = (before..e).map(|j| b.matvec(&v[j])).collect();
        ov.extend(solver.solve_many(&rhs));
        for row in h.iter_mut() {
            row.resize(e, 0.0);
        }
        h.resize(e, vec![0.0; e]);
        for j in before..e {
            let col = par_dots(&v[..e], &ov[j]);
            for i in 0..e {
                h[i][j] = col[i];
            }
        }
        for i in before..e {
            let row = par_dots(&ov[..before], &v[i]);
            h[i][..before].copy_from_slice(&row[..before]);
        }
        let flat: Vec<f64> = (0..e).flat_map(|i| (0..e).map(move |j| (i, j))).map(|(i, j)| h[i][j]).collect();
        let ritz = eig_general(&flat, e)?;
        let mut idx: Vec<usize> = (0..e).collect();
        idx.sort_by(|&a, &b| {
            let (za, zb) = (&ritz[a].value, &ritz[b].value);
            (zb.0.hypot(zb.1)).partial_cmp(&za.0.hypot(za.1)).unwrap()
        });
        let mut out = EigPairs { values: vec![], vectors: vec![], dim: e, converged: false, complex_found: false };
        let mut all_ok = true;
        for (considered, &k) in idx.iter().enumerate() {
            if out.values.len() >= nev || considered >= 3 * nev + 6 {
                break;
            }
            let (re, im) = ritz[k].value;
            // Real or complex: the Ritz values of a real matrix come from its real Schur form, so a real value has an imaginary part of
            // exactly zero and a complex one comes in a conjugate pair; the cutoff only guards against a pair of nearly equal real values
            // that the projection merged, whose imaginary part is below the accuracy `tol` asks of the value itself.
            if im.abs() > 10.0 * opts.tol * re.hypot(im) {
                if considered < nev + 2 {
                    // a genuine complex pair is a wanted mode only if it has converged; an unconverged one is a Ritz artefact, skipped
                    let (yr, yi): (Vec<f64>, Vec<f64>) = ritz[k].vector.iter().map(|z| (z.0, z.1)).unzip();
                    let (xr, xi) = (par_combine(&v[..e], &yr, n), par_combine(&v[..e], &yi, n));
                    let (our, oui) = (par_combine(&ov[..e], &yr, n), par_combine(&ov[..e], &yi, n));
                    let (mut rr, mut ri) = (0.0, 0.0);
                    for t in 0..n {
                        let a = our[t] - (re * xr[t] - im * xi[t]);
                        let b = oui[t] - (re * xi[t] + im * xr[t]);
                        rr += a * a;
                        ri += b * b;
                    }
                    let xn = (dot(&xr, &xr) + dot(&xi, &xi)).sqrt();
                    if (rr + ri).sqrt() / (re.hypot(im) * xn).max(1e-300) <= opts.tol {
                        out.complex_found = true;
                        out.converged = false;
                        return Ok(out);
                    }
                }
                continue;
            }
            // rotate the eigenvector to be real
            let big = ritz[k].vector.iter().cloned().fold((0.0f64, 0.0f64), |m, z| if z.0.hypot(z.1) > m.0.hypot(m.1) { z } else { m });
            let mag = big.0.hypot(big.1);
            let (cr, ci) = (big.0 / mag, -big.1 / mag);
            let y: Vec<f64> = ritz[k].vector.iter().map(|z| z.0 * cr - z.1 * ci).collect();
            let mut x = par_combine(&v[..e], &y, n);
            let mut r = par_combine(&ov[..e], &y, n);
            for t in 0..n {
                r[t] -= re * x[t];
            }
            let xn = dot(&x, &x).sqrt();
            let rel = dot(&r, &r).sqrt() / (re.abs() * xn).max(1e-300);
            if rel > opts.tol {
                all_ok = false;
            }
            x.iter_mut().for_each(|t| *t /= xn);
            out.values.push(sigma + 1.0 / re);
            out.vectors.push(x);
        }
        all_ok &= out.values.len() == nev;
        if all_ok || e + p > max_dim {
            out.converged = all_ok;
            return Ok(out);
        }
        cand = (before..e).map(|j| ov[j].clone()).collect();
        if cand.is_empty() {
            break;
        }
    }
    Err("the Krylov space was exhausted before enough eigenpairs were found".into())
}




/// One eigenpair of a general real matrix.
struct GeneralEig {
    value: (f64, f64),
    vector: Vec<(f64, f64)>,
}

/// Eigenpairs of a general real matrix (row-major), for the projected matrix of the Arnoldi process.
fn eig_general(a: &[f64], n: usize) -> Result<Vec<GeneralEig>, String> {
    let m = faer::Mat::from_fn(n, n, |i, j| a[i * n + j]);
    let e = m.eigen().map_err(|e| format!("eigendecomposition failed: {e:?}"))?;
    let (u, s) = (e.U(), e.S());
    Ok((0..n).map(|j| {
        let z = s.column_vector()[j];
        GeneralEig { value: (z.re, z.im), vector: (0..n).map(|i| (u[(i, j)].re, u[(i, j)].im)).collect() }
    }).collect())
}
