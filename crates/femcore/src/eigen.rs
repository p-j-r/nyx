//! Sparse eigenvalue solvers: block shift-and-invert Krylov methods.
//!
//! Both solvers find the eigenvalues of `A x = lambda B x` closest to a real shift `sigma`
//! by building a Krylov basis of `OP = (A - sigma B)^-1 B` from a *block* of starting vectors
//! (so repeated eigenvalues, common in symmetric structures, are found) and extracting Ritz
//! pairs of the projected matrix. Residuals of the Ritz pairs are evaluated explicitly against
//! the operator, so a returned pair is converged to `tol` regardless of the basis size.

use crate::linalg::{eig_general, eigh, SparseSolver};
use crate::sparse::{dot, Csr};
use rayon::prelude::*;

/// The solve of `(A - sigma M) x = b` supplied by the caller of `eigs_sym_with`.
pub type SolveFn<'a> = &'a dyn Fn(&[f64]) -> Result<Vec<f64>, String>;

/// A projection applied to start vectors and operator results (removes a known null space).
pub type Projector<'a, T> = &'a dyn Fn(&mut Vec<T>);

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

/// Symmetric generalised problem `A x = lambda M x` (A symmetric, M symmetric positive
/// definite), eigenvalues nearest to `sigma`.
pub fn eigs_sym(a: &Csr, m: &Csr, sigma: f64, nev: usize, opts: &EigOpts) -> Result<EigPairs, String> {
    let n = a.nrows;
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    let shifted = a.add_scaled(m, -sigma);
    let solver = SparseSolver::new(&shifted, sigma <= 0.0)?;
    eigs_sym_with(a, m, sigma, nev, opts, &|b: &[f64]| Ok(solver.solve(b)))
}

/// `eigs_sym` with the solves of `(A - sigma M) x = b` supplied by the caller, for instance by
/// multigrid-preconditioned conjugate gradients when `A - sigma M` is positive definite and too
/// large to factorise. The solves must be accurate to about `opts.tol` relative.
pub fn eigs_sym_with(a: &Csr, m: &Csr, sigma: f64, nev: usize, opts: &EigOpts, solve: SolveFn) -> Result<EigPairs, String> {
    eigs_sym_projected(a, m, sigma, nev, opts, solve, None)
}

/// `eigs_sym_with` that also removes, from every start vector and every operator application, the
/// component along a known null space of `A` (for example the discrete gradients of the curl-curl
/// operator, whose eigenvalue 0 would otherwise crowd out the physical modes). `project` must be an
/// `M`-orthogonal projection onto the complement of that null space.
pub fn eigs_sym_projected(a: &Csr, m: &Csr, sigma: f64, nev: usize, opts: &EigOpts, solve: SolveFn, project: Option<Projector<f64>>) -> Result<EigPairs, String> {
    let _ = a;
    let n = a.nrows;
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    let p = if opts.block > 0 { opts.block } else { nev.min(8) };
    let max_dim = max_dim_for(n, nev, opts);
    let _stage = crate::progress::stage(format!("shift-and-invert Krylov eigensolver ({nev} modes, {n} unknowns)"));

    let mnorm2 = |x: &[f64], mx: &[f64]| dot(x, mx);
    let mut v: Vec<Vec<f64>> = Vec::new(); // M-orthonormal basis
    let mut mv: Vec<Vec<f64>> = Vec::new(); // M v
    let mut ov: Vec<Vec<f64>> = Vec::new(); // OP v
    let mut h: Vec<Vec<f64>> = Vec::new(); // h[i][j] = (M v_i) . (OP v_j)

    // start block
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut cand: Vec<Vec<f64>> = (0..p).map(|_| (0..n).map(|_| rng.next()).collect()).collect();
    if let Some(pr) = project {
        cand.iter_mut().for_each(pr);
    }
    let mut s = 0usize;
    loop {
        crate::progress::check()?;
        // M-orthonormalise `cand` against the basis and among themselves
        let before = v.len();
        for mut w in cand.drain(..) {
            let w0 = mnorm2(&w, &m.matvec(&w)).sqrt();
            let mut mw = m.matvec(&w);
            for _pass in 0..2 {
                for i in 0..v.len() {
                    let c = dot(&mv[i], &w);
                    if c != 0.0 {
                        for (x, y) in w.iter_mut().zip(&v[i]) {
                            *x -= c * y;
                        }
                    }
                }
                mw = m.matvec(&w);
            }
            let nrm = dot(&w, &mw).max(0.0).sqrt();
            if nrm > 1e-10 * w0 && nrm > 0.0 && v.len() < max_dim {
                w.iter_mut().for_each(|x| *x /= nrm);
                mw.iter_mut().for_each(|x| *x /= nrm);
                v.push(w);
                mv.push(mw);
            }
        }
        let e = v.len();
        if e == before {
            break; // Krylov space exhausted
        }
        // apply OP to the new columns and extend the projected matrix
        for j in before..e {
            let mut x = solve(&mv[j])?;
            if let Some(pr) = project {
                pr(&mut x);
            }
            ov.push(x);
        }
        for row in h.iter_mut() {
            row.resize(e, 0.0);
        }
        h.resize(e, vec![0.0; e]);
        for j in before..e {
            for i in 0..e {
                h[i][j] = dot(&mv[i], &ov[j]);
            }
        }
        for i in before..e {
            for j in 0..before {
                h[i][j] = dot(&mv[i], &ov[j]);
            }
        }
        s = e;
        // Ritz pairs of the projected matrix
        let flat: Vec<f64> = (0..e).flat_map(|i| (0..e).map(move |j| (i, j))).map(|(i, j)| 0.5 * (h[i][j] + h[j][i])).collect();
        let (nu, y) = eigh(&flat, e)?;
        let mut idx: Vec<usize> = (0..e).collect();
        idx.sort_by(|&a, &b| nu[b].abs().partial_cmp(&nu[a].abs()).unwrap());
        let want = &idx[..nev.min(e)];
        let mut all_ok = want.len() == nev;
        let mut worst_rel = 0.0f64;
        let mut out = EigPairs { values: vec![], vectors: vec![], dim: e, converged: false };
        for &k in want {
            let mut x = vec![0.0; n];
            let mut r = vec![0.0; n];
            for j in 0..e {
                let c = y[k][j];
                if c != 0.0 {
                    for t in 0..n {
                        x[t] += c * v[j][t];
                        r[t] += c * ov[j][t];
                    }
                }
            }
            for t in 0..n {
                r[t] -= nu[k] * x[t];
            }
            let rn = dot(&r, &m.matvec(&r)).max(0.0).sqrt();
            let xn = dot(&x, &m.matvec(&x)).max(0.0).sqrt();
            let rel = rn / (nu[k].abs() * xn).max(1e-300);
            if rel > opts.tol {
                all_ok = false;
            }
            worst_rel = worst_rel.max(rel);
            x.iter_mut().for_each(|t| *t /= xn);
            out.values.push(sigma + 1.0 / nu[k]);
            out.vectors.push(x);
        }
        if all_ok || e + p > max_dim {
            out.converged = all_ok;
            crate::progress::note(format!("Krylov space of {e} vectors, worst residual {worst_rel:.1e}{}", if all_ok { "" } else { " (not converged)" }));
            return Ok(out);
        }
        crate::progress::note(format!("Krylov space {e} vectors, worst residual {worst_rel:.1e}"));
        // next block: the operator applied to the newest basis vectors
        cand = (before..e).map(|j| ov[j].clone()).collect();
        if cand.is_empty() {
            break;
        }
    }
    let _ = s;
    Err("the Krylov space was exhausted before enough eigenpairs were found".into())
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
    let n = a.nrows;
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    let shifted = a.add_scaled(b, -sigma);
    let solver = SparseSolver::new(&shifted, false)?;
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
        let mut out = EigPairs { values: vec![], vectors: vec![], dim: e, converged: false };
        let mut all_ok = true;
        for (considered, &k) in idx.iter().enumerate() {
            if out.values.len() >= nev || considered >= 3 * nev + 6 {
                break;
            }
            let (re, im) = ritz[k].value;
            if im.abs() > 1e-8 * re.abs().max(1e-300) {
                continue; // complex Ritz value: not a physical (real) eigenvalue
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::Coo;

    /// 5-point Laplacian on an n x n grid of interior nodes (Dirichlet), eigenvalues
    /// 4 - 2cos(i pi/(n+1)) - 2cos(j pi/(n+1)); every i != j value is double.
    fn grid(n: usize) -> Csr {
        let id = |i: usize, j: usize| i * n + j;
        let mut c = Coo::new(n * n, n * n);
        for i in 0..n {
            for j in 0..n {
                c.push(id(i, j), id(i, j), 4.0);
                if i > 0 {
                    c.push(id(i, j), id(i - 1, j), -1.0);
                }
                if i + 1 < n {
                    c.push(id(i, j), id(i + 1, j), -1.0);
                }
                if j > 0 {
                    c.push(id(i, j), id(i, j - 1), -1.0);
                }
                if j + 1 < n {
                    c.push(id(i, j), id(i, j + 1), -1.0);
                }
            }
        }
        c.to_csr()
    }

    fn exact(n: usize) -> Vec<f64> {
        let pi = std::f64::consts::PI;
        let mut v = vec![];
        for i in 1..=n {
            for j in 1..=n {
                v.push(4.0 - 2.0 * (i as f64 * pi / (n + 1) as f64).cos() - 2.0 * (j as f64 * pi / (n + 1) as f64).cos());
            }
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v
    }

    #[test]
    fn symmetric_finds_double_eigenvalues() {
        let n = 20;
        let a = grid(n);
        let m = Csr::identity(n * n);
        let r = eigs_sym(&a, &m, -0.1, 7, &EigOpts::default()).unwrap();
        assert!(r.converged);
        let mut got = r.values.clone();
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let ex = exact(n);
        for k in 0..7 {
            assert!((got[k] - ex[k]).abs() < 1e-10, "{} {} {}", k, got[k], ex[k]);
        }
        // eigenvectors satisfy the equation
        for (lam, x) in r.values.iter().zip(&r.vectors) {
            let ax = a.matvec(x);
            let res: f64 = ax.iter().zip(x).map(|(p, q)| (p - lam * q).powi(2)).sum::<f64>().sqrt();
            assert!(res < 1e-8);
        }
    }

    #[test]
    fn symmetric_interior_shift() {
        let n = 16;
        let a = grid(n);
        let m = Csr::identity(n * n);
        let ex = exact(n);
        let target = ex[100] + 1e-3;
        let r = eigs_sym(&a, &m, target, 3, &EigOpts::default()).unwrap();
        let mut got = r.values.clone();
        got.sort_by(|p, q| p.partial_cmp(q).unwrap());
        let mut want: Vec<f64> = ex.to_vec();
        want.sort_by(|p, q| (p - target).abs().partial_cmp(&(q - target).abs()).unwrap());
        let mut want3 = want[..3].to_vec();
        want3.sort_by(|p, q| p.partial_cmp(q).unwrap());
        for k in 0..3 {
            assert!((got[k] - want3[k]).abs() < 1e-9);
        }
    }

    #[test]
    fn general_indefinite_pencil() {
        // A = grid Laplacian, B = diag(+1/-1 alternating): real eigenvalues of A x = lambda B x are
        // not known in closed form; check residuals and symmetry of the problem instead
        let n = 12;
        let a = grid(n);
        let d: Vec<f64> = (0..n * n).map(|i| if i % 7 < 4 { 1.0 } else { -0.5 }).collect();
        let b = Csr::from_diag(&d);
        let r = eigs_general(&a, &b, 1.3, 4, &EigOpts::default()).unwrap();
        assert!(r.converged, "dim {}", r.dim);
        for (lam, x) in r.values.iter().zip(&r.vectors) {
            let ax = a.matvec(x);
            let bx = b.matvec(x);
            let res: f64 = ax.iter().zip(&bx).map(|(p, q)| (p - lam * q).powi(2)).sum::<f64>().sqrt();
            assert!(res < 1e-7, "{res}");
        }
    }
}
