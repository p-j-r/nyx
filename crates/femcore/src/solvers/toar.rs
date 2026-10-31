//! **TOAR** for the quadratic eigenproblem `(λ² M + λ C + K) x = 0` (port item N1): Arnoldi on the
//! shift-and-invert operator of the linearisation, with the `2n`-dimensional Krylov basis stored in
//! compressed form.
//!
//! With `z = [λx; x]` the problem is the pencil `A z = λ B z`, `A = [−C −K; I 0]`, `B = [M 0; 0 I]`.
//! The shift-and-invert operator `(A − σB)⁻¹ B` needs **one** factorisation of `P(σ) = σ²M + σC + K`
//! (size `n`, not `2n`): for `v = [v₁; v₂]`
//!
//! ```text
//! w₂ = −P(σ)⁻¹ (M v₁ + (C + σM) v₂),      w₁ = v₂ + σ w₂.
//! ```
//!
//! TOAR (Lu, Su, Bai 2016) observes that the Krylov vectors satisfy
//! `V_j = [Q U₁ⱼ; Q U₂ⱼ]` with one orthonormal `Q ∈ ℂ^{n×q}`, `q ≤ j + 2`, and small coefficient
//! blocks `U₁, U₂`. Each step orthogonalises the two `n`-blocks of `w` against `Q` (extending it by
//! at most two columns), then orthogonalises the *coefficient* vector against the previous columns in
//! the `2q`-dimensional coefficient space: the basis uses `n(j+2)` scalars instead of `2n(j+1)`, and
//! the orthogonalisation of the long vectors is against `q ≈ j` columns, not `2j`.
//!
//! Restarting is **Krylov–Schur in the compressed form**: after `ncv` steps the Schur form `H Z = Z S` of
//! the Hessenberg matrix is reordered so the Ritz values nearest `σ` come first; the `k` kept Schur vectors
//! `V Z_k` have coefficient blocks `U Z_k`, and the `2(k+1)` coefficient columns (with the residual
//! vector) are compressed to an orthonormal `q × r` basis `W`, `Q ← Q W`, which also caps `q` at `2(k+1)`.
//! There is no locking of converged pairs.
//! Complex arithmetic throughout (real problems are promoted).

use crate::solvers::basis::Basis;
use crate::solvers::csr::Csr;
use crate::solvers::dense::{Dense, DenseError, reorder, schur, triangular_eigenvector};
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::{C64, Scalar};
use crate::solvers::vecops::norm2;

/// `λ² M + λ C + K`.
#[derive(Debug, Clone)]
pub struct QuadraticProblem {
    pub m: Csr<C64>,
    pub c: Csr<C64>,
    pub k: Csr<C64>,
}

impl QuadraticProblem {
    pub fn dim(&self) -> usize {
        self.k.nrows()
    }

    /// `‖P(λ)x‖ / (|λ|²‖Mx‖ + |λ|‖Cx‖ + ‖Kx‖)`.
    pub fn relative_residual(&self, lambda: C64, x: &[C64]) -> f64 {
        let n = self.dim();
        let (mut mx, mut cx, mut kx) = (vec![C64::ZERO; n], vec![C64::ZERO; n], vec![C64::ZERO; n]);
        self.m.spmv(x, &mut mx);
        self.c.spmv(x, &mut cx);
        self.k.spmv(x, &mut kx);
        let l2 = lambda * lambda;
        let r: f64 = (0..n).map(|i| (l2 * mx[i] + lambda * cx[i] + kx[i]).norm_sqr()).sum::<f64>().sqrt();
        r / (l2.norm() * norm2(&mx) + lambda.norm() * norm2(&cx) + norm2(&kx)).max(1e-300)
    }

    /// The `2n × 2n` pencil `(A, B)` of the linearisation (for reference solvers and tests).
    pub fn linearisation(&self) -> (Csr<C64>, Csr<C64>) {
        let n = self.dim();
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (blk, mat, sign) in [((0, 0), &self.c, -1.0), ((0, 1), &self.k, -1.0)] {
            for i in 0..n {
                let (cols, vals) = mat.row(i);
                for (&j, &v) in cols.iter().zip(vals) {
                    a.push((blk.0 * n + i, blk.1 * n + j as usize, v.scale(sign)));
                }
            }
        }
        for i in 0..n {
            a.push((n + i, i, C64::ONE));
            b.push((n + i, n + i, C64::ONE));
            let (cols, vals) = self.m.row(i);
            for (&j, &v) in cols.iter().zip(vals) {
                b.push((i, j as usize, v));
            }
        }
        (Csr::from_triplets(2 * n, 2 * n, &a).unwrap(), Csr::from_triplets(2 * n, 2 * n, &b).unwrap())
    }
}

#[derive(Debug, Clone)]
pub struct ToarOptions {
    pub nev: usize,
    /// Arnoldi steps (the Krylov dimension).
    pub ncv: usize,
    pub sigma: C64,
    /// Relative Ritz error estimate `‖r‖ / |θ|` to accept a pair.
    pub tol: f64,
    pub seed: u64,
    /// Krylov–Schur restarts allowed (0 = a single run of `ncv` steps).
    pub max_restarts: usize,
    /// Fraction of `ncv` kept at a restart (at least `nev + 1`).
    pub keep: f64,
}

impl ToarOptions {
    pub fn new(nev: usize, ncv: usize, sigma: C64) -> Self {
        Self { nev, ncv, sigma, tol: 1e-9, seed: 0x70A2, max_restarts: 100, keep: 0.5 }
    }
}

#[derive(Debug, Clone)]
pub struct ToarResult {
    /// The `nev` eigenvalues nearest `σ` among the Ritz values, nearest first.
    pub values: Vec<C64>,
    /// Their eigenvectors `x` (unit 2-norm), one per column.
    pub vectors: Basis<C64>,
    /// Residual estimates (Arnoldi) and independently recomputed relative residuals.
    pub estimates: Vec<f64>,
    pub residuals: Vec<f64>,
    /// All `nev` pairs have an estimate ≤ `tol`.
    pub converged: bool,
    /// Arnoldi steps taken in the last cycle (less than `ncv` after an invariant subspace was found).
    pub steps: usize,
    pub restarts: usize,
    /// Applications of the shift-and-invert operator (one solve with `P(σ)` each).
    pub operator_applications: usize,
    /// Columns of `Q` at the end.
    pub q_columns: usize,
    /// Scalars stored for the Krylov basis: TOAR (`Q`, `U₁`, `U₂`) versus the `2n × (steps+1)` basis
    /// of Arnoldi on the linearisation.
    pub basis_scalars: usize,
    pub basis_scalars_linearised: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ToarError {
    Solver(SolverError),
    Dense(DenseError),
    InvalidArgument(&'static str),
}

impl std::fmt::Display for ToarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ToarError::Solver(e) => write!(f, "factorisation of P(σ) failed: {e}"),
            ToarError::Dense(e) => write!(f, "dense step failed: {e}"),
            ToarError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
        }
    }
}

impl std::error::Error for ToarError {}
impl From<SolverError> for ToarError {
    fn from(e: SolverError) -> Self {
        ToarError::Solver(e)
    }
}
impl From<DenseError> for ToarError {
    fn from(e: DenseError) -> Self {
        ToarError::Dense(e)
    }
}

/// `Qᴴ w` then `w ← w − Q (Qᴴ w)`, twice (CGS2); returns the accumulated coefficients (length `q`).
fn orth_against(q: &Basis<C64>, qcols: usize, w: &mut [C64]) -> Vec<C64> {
    let mut acc = vec![C64::ZERO; qcols];
    let mut h = vec![C64::ZERO; qcols.max(1)];
    for _ in 0..2 {
        if qcols > 0 {
            q.view().dots_block(qcols, w, 1, &mut h, None);
            q.view().sub_combination(qcols, &h, w);
            for (a, hv) in acc.iter_mut().zip(&h) {
                *a += *hv;
            }
        }
    }
    acc
}

/// Reallocates `q` with more room when all `used` columns are taken (`Q` can need up to `n` columns;
/// in exact arithmetic it needs about `ncv + 2`, but a restart may leave a larger rank).
fn grow_if_full(q: &mut Basis<C64>, used: usize) {
    if used >= q.capacity() {
        let n = q.len();
        let mut bigger = Basis::<C64>::new(n, (q.capacity() * 3 / 2 + 2).min(n));
        bigger.cols_mut(0, used).copy_from_slice(q.cols(0, used));
        *q = bigger;
    }
}

/// Orthonormal basis (modified Gram–Schmidt twice, rank-revealing) of the columns of the `q × cols`
/// matrix `cols_data` (each column a `Vec` of length `q`); returns `W` as columns of length `q`.
fn orthonormal_basis(cols_data: &[Vec<C64>], q: usize) -> Vec<Vec<C64>> {
    let mut w: Vec<Vec<C64>> = Vec::new();
    let scale = cols_data.iter().map(|c| norm2(&c[..q])).fold(0.0, f64::max).max(1e-300);
    for col in cols_data {
        let mut v: Vec<C64> = col[..q].to_vec();
        for _ in 0..2 {
            for b in &w {
                let d: C64 = b.iter().zip(&v).map(|(x, y)| x.conj() * *y).sum();
                for (vi, bi) in v.iter_mut().zip(b) {
                    *vi -= d * *bi;
                }
            }
        }
        let nv = norm2(&v);
        if nv > 1e-10 * scale && w.len() < q {
            v.iter_mut().for_each(|x| *x = x.scale(1.0 / nv));
            w.push(v);
        }
    }
    w
}

pub fn toar_shift_invert(p: &QuadraticProblem, opts: &ToarOptions) -> Result<ToarResult, ToarError> {
    let n = p.dim();
    let m_steps = opts.ncv;
    if n == 0 || p.m.nrows() != n || p.c.nrows() != n || opts.nev == 0 || opts.nev > m_steps || m_steps > 2 * n {
        return Err(ToarError::InvalidArgument("need matching M, C, K, 0 < nev ≤ ncv ≤ 2n"));
    }
    let sigma = opts.sigma;
    // P(σ) = σ²M + σC + K
    let p_sigma = Csr::axpby(C64::ONE, &Csr::axpby(sigma * sigma, &p.m, sigma, &p.c), C64::ONE, &p.k);
    let mut fac = ShiftedFactorization::new(&p_sigma, None, Definiteness::General)?;
    fac.set_shift(C64::ZERO)?;
    // C + σM, for the second term of the right-hand side
    let c_sm = Csr::axpby(C64::ONE, &p.c, sigma, &p.m);

    // Q has at most n columns (it is a basis of a subspace of ℂⁿ), whatever the Krylov dimension
    let qmax = (m_steps + 2).min(n);
    let mut q = Basis::<C64>::new(n, qmax);
    let mut qc = 0usize; // columns of Q in use
    // coefficient columns, each of length qmax (zero beyond the current q)
    // (q ≤ 2·ncv + 2 always: at most two new directions per step, and a restart leaves ≤ 2(k+1))
    let ulen = (2 * m_steps + 2).min(n);
    let mut u1: Vec<Vec<C64>> = vec![vec![C64::ZERO; ulen]; m_steps + 1];
    let mut u2: Vec<Vec<C64>> = vec![vec![C64::ZERO; ulen]; m_steps + 1];
    let mut hess = Dense::<C64>::zeros(m_steps + 1, m_steps);

    // extends Q by the (up to two) directions of r1, r2 and returns the coefficients of the two
    // n-blocks in the enlarged Q, given their coefficients against the old Q
    let extend = |q: &mut Basis<C64>,
                  qc: &mut usize,
                  mut r1: Vec<C64>,
                  mut r2: Vec<C64>,
                  mut c1: Vec<C64>,
                  mut c2: Vec<C64>,
                  scale: f64| {
        let thr = 1e-12 * scale.max(1e-300);
        let a1 = norm2(&r1);
        let added1 = a1 > thr && *qc < q.len();
        if added1 {
            grow_if_full(q, *qc);
        }
        if added1 {
            r1.iter_mut().for_each(|v| *v = v.scale(1.0 / a1));
            q.col_mut(*qc).copy_from_slice(&r1);
            *qc += 1;
            // the second block's remainder is orthogonalised against the new column
            let b1: C64 = r1.iter().zip(&r2).map(|(a, b)| a.conj() * *b).sum();
            for (x, y) in r2.iter_mut().zip(&r1) {
                *x -= b1 * *y;
            }
            c1.push(C64::new(a1, 0.0));
            c2.push(b1);
        }
        let a2 = norm2(&r2);
        if a2 > thr && *qc < q.len() {
            grow_if_full(q, *qc);
            r2.iter_mut().for_each(|v| *v = v.scale(1.0 / a2));
            q.col_mut(*qc).copy_from_slice(&r2);
            *qc += 1;
            c1.push(C64::ZERO);
            c2.push(C64::new(a2, 0.0));
        }
        (c1, c2)
    };

    // starting vector
    {
        let s = random_vector::<C64>(2 * n, opts.seed);
        let (s1, s2) = (s[..n].to_vec(), s[n..].to_vec());
        let (c1, c2) = extend(&mut q, &mut qc, s1, s2, Vec::new(), Vec::new(), 1.0);
        let nu = (c1.iter().chain(&c2).map(|v| v.norm_sqr()).sum::<f64>()).sqrt();
        for (i, v) in c1.iter().enumerate() {
            u1[0][i] = v.scale(1.0 / nu);
        }
        for (i, v) in c2.iter().enumerate() {
            u2[0][i] = v.scale(1.0 / nu);
        }
    }

    let (mut x1, mut x2) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    let (mut rhs, mut tmp) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    let mut j_start = 0usize; // first Arnoldi step of the current cycle
    let mut restarts = 0usize;
    let mut applications = 0usize;
    let upper = m_steps.saturating_sub(1).max(1);
    let keep_target = ((m_steps as f64 * opts.keep) as usize).clamp((opts.nev + 1).min(upper), upper);

    loop {
        // ---- Arnoldi steps j_start..m_steps
        let mut steps = j_start;
        for j in j_start..m_steps {
            x1.fill(C64::ZERO);
            x2.fill(C64::ZERO);
            q.view().add_combination(qc, &u1[j][..qc], &mut x1);
            q.view().add_combination(qc, &u2[j][..qc], &mut x2);
            // w2 = −P⁻¹ (M x1 + (C + σM) x2), w1 = x2 + σ w2
            p.m.spmv(&x1, &mut rhs);
            c_sm.spmv(&x2, &mut tmp);
            for (r, t) in rhs.iter_mut().zip(&tmp) {
                *r = -(*r + *t);
            }
            fac.solve_block_in_place(&mut rhs, 1);
            applications += 1;
            let w2 = rhs.clone();
            let w1: Vec<C64> = x2.iter().zip(&w2).map(|(&a, &b)| a + sigma * b).collect();
            let wnorm = (norm2(&w1).powi(2) + norm2(&w2).powi(2)).sqrt();
            let (mut r1, mut r2) = (w1, w2);
            let c1 = orth_against(&q, qc, &mut r1);
            let c2 = orth_against(&q, qc, &mut r2);
            let (c1, c2) = extend(&mut q, &mut qc, r1, r2, c1, c2, wnorm);
            let qn = c1.len();
            // coefficient-space orthogonalisation against columns 0..=j (CGS2)
            let (mut v1, mut v2) = (c1, c2);
            let mut hcol = vec![C64::ZERO; j + 2];
            for _ in 0..2 {
                for k in 0..=j {
                    let mut d = C64::ZERO;
                    for i in 0..qn {
                        d += u1[k][i].conj() * v1[i] + u2[k][i].conj() * v2[i];
                    }
                    for i in 0..qn {
                        v1[i] -= d * u1[k][i];
                        v2[i] -= d * u2[k][i];
                    }
                    hcol[k] += d;
                }
            }
            let nu = (v1.iter().chain(&v2).map(|v| v.norm_sqr()).sum::<f64>()).sqrt();
            for (k, h) in hcol.iter().enumerate().take(j + 1) {
                hess[(k, j)] = *h;
            }
            hess[(j + 1, j)] = C64::new(nu, 0.0);
            steps = j + 1;
            if nu <= 1e-13 * wnorm.max(1e-300) {
                hess[(j + 1, j)] = C64::ZERO;
                break; // invariant subspace: the Ritz pairs are exact
            }
            for i in 0..qn {
                u1[j + 1][i] = v1[i].scale(1.0 / nu);
                u2[j + 1][i] = v2[i].scale(1.0 / nu);
            }
        }

        // ---- Ritz extraction: ordered Schur form of H_m
        let mm = steps;
        let hm = Dense::from_fn(mm, mm, |i, j| hess[(i, j)]);
        let (mut s, mut z) = schur(&hm)?;
        reorder(&mut s, &mut z, mm, |c| -c.norm()); // largest |θ| (nearest σ) first
        let h_last = hess[(mm, mm - 1)];
        // residual coupling row r = h_{m+1,m} Zᵀ-row m−1 ; estimate of pair p: ‖r·y_p‖/‖y_p‖/|θ_p|
        let rrow: Vec<C64> = (0..mm).map(|i| h_last * z[(mm - 1, i)]).collect();
        let nev = opts.nev.min(mm);
        let mut estimates = Vec::with_capacity(nev);
        let mut ys: Vec<Vec<C64>> = Vec::with_capacity(nev);
        for pidx in 0..nev {
            let y = triangular_eigenvector(&s, mm, pidx);
            let ynorm = norm2(&y);
            let val: C64 = (0..=pidx).map(|i| rrow[i] * y[i]).sum();
            estimates.push(val.norm() / ynorm / s[(pidx, pidx)].norm().max(1e-300));
            ys.push(y);
        }
        let done = (nev == opts.nev && estimates.iter().all(|&e| e <= opts.tol)) || h_last.norm() == 0.0;
        if done || restarts >= opts.max_restarts {
            // ---- final extraction
            let mut values = Vec::with_capacity(nev);
            let mut vectors = Basis::<C64>::new(n, nev.max(1));
            for (slot, y) in ys.iter().enumerate() {
                let theta = s[(slot, slot)];
                // Z y combines the Schur vectors; x = Q (U2[:, :m] (Z y))
                let zy: Vec<C64> = (0..mm).map(|i| (0..mm).map(|k| z[(i, k)] * y[k]).sum()).collect();
                let mut coef = vec![C64::ZERO; qc];
                for (k, yk) in zy.iter().enumerate() {
                    for (cv, uv) in coef.iter_mut().zip(&u2[k][..qc]) {
                        *cv += *uv * *yk;
                    }
                }
                let x = vectors.col_mut(slot);
                q.view().add_combination(qc, &coef, x);
                let nx = norm2(x);
                x.iter_mut().for_each(|v| *v = v.scale(1.0 / nx));
                values.push(sigma + C64::ONE / theta);
            }
            let residuals: Vec<f64> =
                values.iter().enumerate().map(|(i, &l)| p.relative_residual(l, vectors.col(i))).collect();
            let converged = values.len() == opts.nev && estimates.iter().all(|&e| e <= opts.tol);
            return Ok(ToarResult {
                values,
                vectors,
                estimates,
                residuals,
                converged,
                steps,
                restarts,
                operator_applications: applications,
                q_columns: qc,
                basis_scalars: n * qc + 2 * qc * (steps + 1),
                basis_scalars_linearised: 2 * n * (steps + 1),
            });
        }

        // ---- Krylov–Schur restart in compressed form
        restarts += 1;
        let k = keep_target.min(mm - 1).max(1);
        // coefficient columns of V Z_k (k of them) and of the residual vector v_{m+1}
        let q_old = qc;
        let mut new1: Vec<Vec<C64>> = Vec::with_capacity(k + 1);
        let mut new2: Vec<Vec<C64>> = Vec::with_capacity(k + 1);
        for i in 0..k {
            let (mut a1, mut a2) = (vec![C64::ZERO; q_old], vec![C64::ZERO; q_old]);
            for l in 0..mm {
                let zl = z[(l, i)];
                for r in 0..q_old {
                    a1[r] += u1[l][r] * zl;
                    a2[r] += u2[l][r] * zl;
                }
            }
            new1.push(a1);
            new2.push(a2);
        }
        new1.push(u1[mm][..q_old].to_vec());
        new2.push(u2[mm][..q_old].to_vec());
        // compress: orthonormal W spanning all 2(k+1) coefficient columns
        let all_cols: Vec<Vec<C64>> = new1.iter().chain(new2.iter()).cloned().collect();
        let w = orthonormal_basis(&all_cols, q_old);
        let r = w.len();
        let wflat: Vec<C64> = w.iter().flat_map(|c| c.iter().copied()).collect();
        q.mult_in_place(0, q_old, r, &wflat, q_old);
        qc = r;
        for col in u1.iter_mut().chain(u2.iter_mut()) {
            col.fill(C64::ZERO);
        }
        for i in 0..=k {
            for (a, wa) in w.iter().enumerate() {
                let (mut d1, mut d2) = (C64::ZERO, C64::ZERO);
                for t in 0..q_old {
                    d1 += wa[t].conj() * new1[i][t];
                    d2 += wa[t].conj() * new2[i][t];
                }
                u1[i][a] = d1;
                u2[i][a] = d2;
            }
        }
        // Hessenberg: T_k on the leading block, coupling row below it
        hess = Dense::zeros(m_steps + 1, m_steps);
        for i in 0..k {
            for j in i..k {
                hess[(i, j)] = s[(i, j)];
            }
            hess[(k, i)] = rrow[i];
        }
        j_start = k;
    }
}
