//! Complex sparse matrices and a complex block shift-and-invert Arnoldi eigensolver, for problems with absorbing layers
//! (perfectly matched layers) or loss, whose matrices are complex symmetric.
//!
//! `CCoo`/`CCsr` mirror `Coo`/`Csr` for `c64` entries. `CSolver` is a sparse complex LU (faer, on the factorisation pool).
//! `eigs_complex` finds the eigenpairs of `A x = lambda B x` nearest a complex shift with full orthogonalisation (done in
//! parallel); convergence is judged by the exact residual of the Ritz pairs, obtained from the part of the newest block
//! that is not yet in the basis. Every step checks for cancellation and writes a note to the stage log.

use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::Lu;
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::{c64, Mat};
use rayon::prelude::*;

pub type C = c64;

pub fn c(re: f64, im: f64) -> C {
    C::new(re, im)
}

#[derive(Clone, Debug, Default)]
pub struct CCoo {
    pub nrows: usize,
    pub ncols: usize,
    pub rows: Vec<u32>,
    pub cols: Vec<u32>,
    pub vals: Vec<C>,
}

impl CCoo {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        CCoo { nrows, ncols, ..Default::default() }
    }
    pub fn push_block(&mut self, rows: &[usize], cols: &[usize], local: &[C]) {
        for (a, &r) in rows.iter().enumerate() {
            for (b, &cc) in cols.iter().enumerate() {
                let v = local[a * cols.len() + b];
                if v.re != 0.0 || v.im != 0.0 {
                    self.rows.push(r as u32);
                    self.cols.push(cc as u32);
                    self.vals.push(v);
                }
            }
        }
    }
    pub fn to_csr(&self) -> CCsr {
        let n = self.nrows;
        let mut count = vec![0usize; n + 1];
        for &r in &self.rows {
            count[r as usize + 1] += 1;
        }
        for i in 0..n {
            count[i + 1] += count[i];
        }
        let mut order = vec![0usize; self.rows.len()];
        let mut next = count.clone();
        for k in 0..self.rows.len() {
            let r = self.rows[k] as usize;
            order[next[r]] = k;
            next[r] += 1;
        }
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        let mut row: Vec<(u32, C)> = vec![];
        for r in 0..n {
            row.clear();
            for &k in &order[count[r]..count[r + 1]] {
                row.push((self.cols[k], self.vals[k]));
            }
            row.sort_unstable_by_key(|e| e.0);
            let mut i = 0;
            while i < row.len() {
                let col = row[i].0;
                let mut s = c(0.0, 0.0);
                while i < row.len() && row[i].0 == col {
                    s += row[i].1;
                    i += 1;
                }
                indices.push(col as usize);
                data.push(s);
            }
            indptr.push(indices.len());
        }
        CCsr { nrows: n, ncols: self.ncols, indptr, indices, data }
    }
}

#[derive(Clone, Debug)]
pub struct CCsr {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<C>,
}

impl CCsr {
    pub fn matvec(&self, x: &[C]) -> Vec<C> {
        (0..self.nrows)
            .map(|r| {
                let mut s = c(0.0, 0.0);
                for k in self.indptr[r]..self.indptr[r + 1] {
                    s += self.data[k] * x[self.indices[k]];
                }
                s
            })
            .collect()
    }
    pub fn transpose(&self) -> CCsr {
        let mut count = vec![0usize; self.ncols + 1];
        for &cc in &self.indices {
            count[cc + 1] += 1;
        }
        for i in 0..self.ncols {
            count[i + 1] += count[i];
        }
        let mut next = count.clone();
        let mut indices = vec![0usize; self.data.len()];
        let mut data = vec![c(0.0, 0.0); self.data.len()];
        for r in 0..self.nrows {
            for k in self.indptr[r]..self.indptr[r + 1] {
                let cc = self.indices[k];
                indices[next[cc]] = r;
                data[next[cc]] = self.data[k];
                next[cc] += 1;
            }
        }
        CCsr { nrows: self.ncols, ncols: self.nrows, indptr: count, indices, data }
    }
    /// `self + alpha * other`.
    pub fn add_scaled(&self, other: &CCsr, alpha: C) -> CCsr {
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        for r in 0..self.nrows {
            let (mut i, ie) = (self.indptr[r], self.indptr[r + 1]);
            let (mut j, je) = (other.indptr[r], other.indptr[r + 1]);
            while i < ie || j < je {
                let ci = if i < ie { self.indices[i] } else { usize::MAX };
                let cj = if j < je { other.indices[j] } else { usize::MAX };
                if ci == cj {
                    indices.push(ci);
                    data.push(self.data[i] + alpha * other.data[j]);
                    i += 1;
                    j += 1;
                } else if ci < cj {
                    indices.push(ci);
                    data.push(self.data[i]);
                    i += 1;
                } else {
                    indices.push(cj);
                    data.push(alpha * other.data[j]);
                    j += 1;
                }
            }
            indptr.push(indices.len());
        }
        CCsr { nrows: self.nrows, ncols: self.ncols, indptr, indices, data }
    }
    pub fn submatrix(&self, rows: &[usize], cols: &[usize]) -> CCsr {
        let mut map = vec![usize::MAX; self.ncols];
        for (j, &cc) in cols.iter().enumerate() {
            map[cc] = j;
        }
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        let mut row: Vec<(usize, C)> = vec![];
        for &r in rows {
            row.clear();
            for k in self.indptr[r]..self.indptr[r + 1] {
                let j = map[self.indices[k]];
                if j != usize::MAX {
                    row.push((j, self.data[k]));
                }
            }
            row.sort_unstable_by_key(|e| e.0);
            for &(j, v) in &row {
                indices.push(j);
                data.push(v);
            }
            indptr.push(indices.len());
        }
        CCsr { nrows: rows.len(), ncols: cols.len(), indptr, indices, data }
    }
}

pub struct CSolver {
    lu: Lu<usize, C>,
    n: usize,
}

impl CSolver {
    pub fn new(a: &CCsr) -> Result<Self, String> {
        let n = a.nrows;
        let t = a.transpose(); // CSR of A^T is CSC of A
        let sym = SymbolicSparseColMatRef::new_checked(n, n, &t.indptr, None, &t.indices);
        let mat = SparseColMatRef::new(sym, &t.data);
        crate::par::on_factor_pool(|| match mat.sp_lu() {
            Ok(lu) => Ok(CSolver { lu, n }),
            Err(e) => Err(format!("complex sparse factorisation failed (singular matrix?): {e:?}")),
        })
    }
    pub fn solve(&self, b: &[C]) -> Vec<C> {
        let mut rhs = Mat::<C>::from_fn(self.n, 1, |i, _| b[i]);
        self.lu.solve_in_place(rhs.as_mut());
        (0..self.n).map(|i| rhs[(i, 0)]).collect()
    }

    /// Several right-hand sides at once.
    pub fn solve_many(&self, bs: &[Vec<C>]) -> Vec<Vec<C>> {
        let nt = if bs.len() > 1 && self.n > 5_000 { threads().min(bs.len()) } else { 1 };
        if nt <= 1 {
            let mut rhs = Mat::<C>::from_fn(self.n, bs.len(), |i, j| bs[j][i]);
            self.lu.solve_in_place(rhs.as_mut());
            return (0..bs.len()).map(|j| (0..self.n).map(|i| rhs[(i, j)]).collect()).collect();
        }
        // the triangular solves of different columns are independent: one thread per group of columns
        let chunk = bs.len().div_ceil(nt);
        let n = self.n;
        let lu = &self.lu;
        let parts: Vec<Vec<Vec<C>>> = bs
            .par_chunks(chunk)
            .map(|cols| {
                let mut rhs = Mat::<C>::from_fn(n, cols.len(), |i, j| cols[j][i]);
                lu.solve_in_place(rhs.as_mut());
                (0..cols.len()).map(|j| (0..n).map(|i| rhs[(i, j)]).collect()).collect::<Vec<Vec<C>>>()
            })
            .collect();
        parts.into_iter().flatten().collect()
    }
}

#[derive(Clone, Debug)]
pub struct CEigPairs {
    pub values: Vec<C>,
    pub vectors: Vec<Vec<C>>,
    /// relative residual of every pair (the last measured)
    pub residuals: Vec<f64>,
    pub dim: usize,
    pub converged: bool,
}

fn cdot(a: &[C], b: &[C]) -> C {
    // conjugate-linear in the first argument
    a.iter().zip(b).fold(c(0.0, 0.0), |s, (x, y)| s + x.conj() * *y)
}

fn norm(a: &[C]) -> f64 {
    a.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt()
}

fn threads() -> usize {
    crate::par::threads()
}

/// `<q_k, w>` (conjugate-linear in `q_k`) for every vector of `vs`, one thread per group of vectors.
fn par_cdots(vs: &[Vec<C>], w: &[C]) -> Vec<C> {
    let m = vs.len();
    let mut out = vec![c(0.0, 0.0); m];
    let nt = if m * w.len() < 400_000 { 1 } else { threads() };
    let chunk = m.div_ceil(nt).max(1);
    out.par_chunks_mut(chunk).enumerate().for_each(|(ci, o)| {
        for (k, x) in o.iter_mut().enumerate() {
            *x = cdot(&vs[ci * chunk + k], w);
        }
    });
    out
}

/// `w -= sum_i <v_i, w> v_i` over `basis`, done twice (classical Gram-Schmidt with
/// reorthogonalisation), with the dot products and the update on all cores.
fn orthogonalise(w: &mut [C], basis: &[Vec<C>]) {
    let (n, m) = (w.len(), basis.len());
    if m == 0 {
        return;
    }
    let nt = if n * m < 200_000 { 1 } else { threads() };
    for _ in 0..2 {
        // coefficients: one dot product per basis vector
        let mut coef = vec![c(0.0, 0.0); m];
        {
            let wr: &[C] = w;
            let chunk = m.div_ceil(nt).max(1);
            coef.par_chunks_mut(chunk).enumerate().for_each(|(ci, out)| {
                for (k, o) in out.iter_mut().enumerate() {
                    *o = cdot(&basis[ci * chunk + k], wr);
                }
            });
        }
        // update by blocks of rows
        let chunk = n.div_ceil(nt).max(1);
        let coef = &coef;
        w.par_chunks_mut(chunk).enumerate().for_each(|(ci, ws)| {
            let r0 = ci * chunk;
            let len = ws.len();
            for (i, q) in basis.iter().enumerate() {
                let s = coef[i];
                for (x, y) in ws.iter_mut().zip(&q[r0..r0 + len]) {
                    *x -= s * *y;
                }
            }
        });
    }
}

/// Eigenpairs of `A x = lambda B x` nearest the complex shift `sigma` (block shift-and-invert
/// Arnoldi with full orthogonalisation; A - sigma B must be nonsingular). Convergence is judged
/// by the exact residual of the Ritz pairs, which for a block Arnoldi process follows from the
/// part of the newest block that is not yet in the basis, so no vectors have to be assembled
/// until the pairs are final.
pub fn eigs_complex(a: &CCsr, b: &CCsr, sigma: C, nev: usize, tol: f64, max_dim: usize) -> Result<CEigPairs, String> {
    eigs_complex_projected(a, b, sigma, nev, tol, max_dim, None)
}

/// `eigs_complex` with a projection that removes a known null space of `A` from every start vector and
/// every operator application (see `eigen::eigs_sym_projected`).
pub fn eigs_complex_projected(a: &CCsr, b: &CCsr, sigma: C, nev: usize, tol: f64, max_dim: usize, project: Option<crate::eigen::Projector<C>>) -> Result<CEigPairs, String> {
    let n = a.nrows;
    if nev == 0 || nev >= n {
        return Err("number of eigenvalues must be between 1 and n-1".into());
    }
    crate::progress::check()?;
    let _stage = crate::progress::stage(format!("complex eigenmodes ({nev} wanted, {n} unknowns)"));
    let t_all = std::time::Instant::now();
    let (mut t_solve, mut t_eig, mut t_orth) = (0.0f64, 0.0f64, 0.0f64);
    let solver = {
        let _f = crate::progress::stage(format!("factorising the complex matrix ({n} unknowns); this cannot be interrupted"));
        CSolver::new(&a.add_scaled(b, -sigma))?
    };
    crate::progress::check()?;
    let p = nev.min(8);
    let max_dim = max_dim.min(n).max(nev + 8);
    let mut v: Vec<Vec<C>> = vec![];
    let mut ov: Vec<Vec<C>> = vec![];
    let mut h: Vec<Vec<C>> = vec![];
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        ((seed >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    };
    let mut cand: Vec<Vec<C>> = (0..p).map(|_| (0..n).map(|_| c(rnd(), rnd())).collect()).collect();
    if let Some(pr) = project {
        cand.iter_mut().for_each(pr);
    }
    // Ritz pairs of the previous step: values nu, coefficient vectors over the basis, and the index
    // where the last block of the basis starts
    type Ritz = (Vec<C>, Vec<Vec<C>>, usize);
    let mut prev: Option<Ritz> = None;
    let mut last_pairs: Option<Ritz> = None;
    loop {
        let t = std::time::Instant::now();
        let before = v.len();
        let mut rems: Vec<Vec<C>> = vec![];
        for mut w in cand.drain(..) {
            let w0 = norm(&w);
            orthogonalise(&mut w, &v[..before]);
            if prev.is_some() {
                rems.push(w.clone());
            }
            let added: Vec<Vec<C>> = v[before..].to_vec();
            orthogonalise(&mut w, &added);
            let nr = norm(&w);
            if nr > 1e-10 * w0 && nr > 0.0 && v.len() < max_dim {
                w.iter_mut().for_each(|x| *x /= nr);
                v.push(w);
            }
        }
        t_orth += t.elapsed().as_secs_f64();
        // residuals of the previous Ritz pairs from the remainders of the last block
        let mut all_ok = false;
        let mut res_now: Vec<f64> = vec![];
        if let Some((nus, ys, lb)) = &prev {
            all_ok = true;
            for (nu, y) in nus.iter().zip(ys) {
                let mut r = vec![c(0.0, 0.0); n];
                for (jj, rem) in rems.iter().enumerate() {
                    let cf = y[lb + jj];
                    for (x, q) in r.iter_mut().zip(rem) {
                        *x += cf * *q;
                    }
                }
                let rel = norm(&r) / nu.norm().max(1e-300);
                res_now.push(rel);
                if rel > tol {
                    all_ok = false;
                }
            }
        }
        let e = v.len();
        // soft stop: the leading modes are converged and the rest of the block is clustered radiation
        // modes that would need many more steps; the caller drops pairs that are not converged
        let soft = prev.is_some() && before >= 8 * nev + 80 && res_now.first().is_some_and(|r| *r < 1e-8);
        if all_ok || soft || e == before || (e + p > max_dim && (prev.is_some() || last_pairs.is_some())) {
            // final pairs: those of the previous step (or the last ones computed)
            let (nus, ys, used) = prev.take().or(last_pairs.take()).ok_or("the Krylov space was exhausted before enough eigenpairs were found")?;
            let _ = used;
            let mut out = CEigPairs { values: vec![], vectors: vec![], residuals: res_now.clone(), dim: before, converged: all_ok };
            for (nu, y) in nus.iter().zip(&ys) {
                let mut x = vec![c(0.0, 0.0); n];
                for (j, cf) in y.iter().enumerate() {
                    for (t, q) in x.iter_mut().zip(&v[j]) {
                        *t += *cf * *q;
                    }
                }
                let xn = norm(&x);
                x.iter_mut().for_each(|t| *t /= xn);
                out.values.push(sigma + c(1.0, 0.0) / *nu);
                out.vectors.push(x);
            }
            crate::progress::note(format!("Krylov space of {before} vectors: orthogonalising {t_orth:.1} s, solves {t_solve:.1} s, small eigenproblems {t_eig:.1} s, total {:.1} s{}", t_all.elapsed().as_secs_f64(), if all_ok { "" } else { " (not all modes converged)" }));
            return Ok(out);
        }
        crate::progress::check()?;
        if !res_now.is_empty() {
            let worst = res_now.iter().cloned().fold(0.0, f64::max);
            crate::progress::note(format!("Krylov space {e} vectors, worst residual of the {} wanted modes {worst:.1e} (goal {tol:.0e})", res_now.len()));
        }
        // apply the operator to the new block
        let t = std::time::Instant::now();
        let new: Vec<Vec<C>> = (before..e).map(|j| b.matvec(&v[j])).collect();
        let mut sols = solver.solve_many(&new);
        if let Some(pr) = project {
            sols.iter_mut().for_each(pr);
        }
        ov.extend(sols);
        t_solve += t.elapsed().as_secs_f64();
        for row in h.iter_mut() {
            row.resize(e, c(0.0, 0.0));
        }
        h.resize(e, vec![c(0.0, 0.0); e]);
        for j in before..e {
            let col = par_cdots(&v[..e], &ov[j]);
            for i in 0..e {
                h[i][j] = col[i];
            }
        }
        for i in before..e {
            let row = par_cdots(&ov[..before], &v[i]);
            for j in 0..before {
                // <v_i, ov_j> = conj(<ov_j, v_i>)
                h[i][j] = row[j].conj();
            }
        }
        // Ritz pairs of the projected matrix (every step for small bases, then every other step)
        let t = std::time::Instant::now();
        if e <= 120 || (e / p) % 2 == 0 {
            let hm = Mat::<C>::from_fn(e, e, |i, j| h[i][j]);
            let ev = hm.eigen().map_err(|er| format!("eigendecomposition failed: {er:?}"))?;
            let (u, s) = (ev.U(), ev.S());
            let mut idx: Vec<usize> = (0..e).collect();
            idx.sort_by(|&x, &y| s.column_vector()[y].norm().partial_cmp(&s.column_vector()[x].norm()).unwrap());
            let (mut nus, mut ys) = (vec![], vec![]);
            for &k in idx.iter().take(nev) {
                nus.push(s.column_vector()[k]);
                ys.push((0..e).map(|j| u[(j, k)]).collect::<Vec<C>>());
            }
            // the residual of a pair in the next step involves the rows of the block just added
            prev = Some((nus.clone(), ys.clone(), before));
            last_pairs = Some((nus, ys, before));
        } else {
            prev = None;
        }
        t_eig += t.elapsed().as_secs_f64();
        cand = (before..e).map(|j| ov[j].clone()).collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_eigenvalues_of_a_shifted_diagonal_plus_coupling() {
        // A = tridiagonal with complex diagonal, B = I; compare with dense eigenvalues
        let n = 40;
        let mut coo = CCoo::new(n, n);
        for i in 0..n {
            coo.push_block(&[i], &[i], &[c(i as f64 * 0.5, 0.1 * (i % 5) as f64)]);
            if i + 1 < n {
                coo.push_block(&[i], &[i + 1], &[c(0.3, 0.05)]);
                coo.push_block(&[i + 1], &[i], &[c(0.3, 0.05)]);
            }
        }
        let a = coo.to_csr();
        let mut id = CCoo::new(n, n);
        for i in 0..n {
            id.push_block(&[i], &[i], &[c(1.0, 0.0)]);
        }
        let b = id.to_csr();
        let r = eigs_complex(&a, &b, c(5.0, 0.0), 3, 1e-10, 40).unwrap();
        assert!(r.converged);
        for (lam, x) in r.values.iter().zip(&r.vectors) {
            let ax = a.matvec(x);
            let res: f64 = ax.iter().zip(x).map(|(p, q)| (*p - *lam * *q).norm_sqr()).sum::<f64>().sqrt();
            assert!(res < 1e-8, "{res}");
        }
    }
}
