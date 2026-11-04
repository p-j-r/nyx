//! LOBPCG (Knyazev 2001) for the generalised Hermitian problem `A x = λ B x` (port item E4), real or complex,
//! with the improvements the port list names:
//!
//! * **complex Hermitian** matrices (the same code as real: inner products are conjugated);
//! * **constraints**: the search is kept `B`-orthogonal to given vectors `Y` (deflation of known modes);
//! * **soft locking**: converged Ritz vectors stay in the block (so they keep improving through the Rayleigh–Ritz
//!   step) but their residuals no longer generate search directions, so the cost per iteration shrinks as
//!   the block converges;
//! * a **stall test**: the iteration reports `Stalled` when the worst active residual has not improved by 1 %
//!   over `stall_window` iterations (instead of running to `max_it`);
//! * a robust basis: the new directions `W` (preconditioned residuals) and `P` are `B`-orthonormalised
//!   against the block, the constraints and each other (two Gram–Schmidt passes, dependent columns dropped),
//!   which avoids the ill-conditioned Gram matrices of the textbook variant.
//!
//! `which` selects the lowest or the highest eigenvalues. For interior eigenvalues use Jacobi–Davidson
//! ([`crate::solvers::jd`]) or the contour solver ([`crate::solvers::contour`]).

use crate::solvers::csr::Csr;
use crate::solvers::dense::{Dense, DenseError, hermitian_eig};
use crate::solvers::faer_glue::FaerScalar;
use crate::solvers::linop::Preconditioner;
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::{dotc, norm2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LobpcgWhich {
    Lowest,
    Highest,
}

#[derive(Debug, Clone)]
pub struct LobpcgOptions {
    pub nev: usize,
    /// Block size (≥ `nev`; a few extra vectors speed up convergence of the last wanted ones).
    pub block: usize,
    pub tol: f64,
    pub max_it: usize,
    pub which: LobpcgWhich,
    /// Iterations without 1 % improvement of the worst active residual before reporting a stall.
    pub stall_window: usize,
    /// Relative improvement of the worst active residual that counts as progress (0.01 = 1 %).
    pub stall_gain: f64,
    pub seed: u64,
}

impl LobpcgOptions {
    pub fn new(nev: usize) -> Self {
        Self {
            nev,
            block: nev + nev.div_ceil(4).max(2),
            tol: 1e-8,
            max_it: 500,
            which: LobpcgWhich::Lowest,
            stall_window: 200,
            stall_gain: 0.01,
            seed: 0x10B9,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LobpcgReason {
    Converged,
    MaxIterations,
    Stalled,
}

#[derive(Debug, Clone)]
pub struct LobpcgResult<T: Scalar> {
    pub values: Vec<f64>,
    /// `nev` vectors of length `n`, `B`-orthonormal, stored column after column.
    pub vectors: Vec<Vec<T>>,
    /// Relative residuals `‖Ax − λBx‖ / (‖Ax‖ + |λ|‖Bx‖)`.
    pub residuals: Vec<f64>,
    pub iterations: usize,
    pub reason: LobpcgReason,
    /// Matrix products with `A` (block columns counted individually).
    pub matvecs: usize,
    /// Active (non-locked) block columns per iteration.
    pub active_history: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LobpcgError {
    InvalidArgument(&'static str),
    Dense(DenseError),
}

impl std::fmt::Display for LobpcgError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LobpcgError::InvalidArgument(m) => write!(f, "invalid argument: {m}"),
            LobpcgError::Dense(e) => write!(f, "dense step failed: {e}"),
        }
    }
}

impl std::error::Error for LobpcgError {}
impl From<DenseError> for LobpcgError {
    fn from(e: DenseError) -> Self {
        LobpcgError::Dense(e)
    }
}

/// `B x` (or `x` for `B = None`).
fn apply_b<T: Scalar>(b: Option<&Csr<T>>, x: &[T]) -> Vec<T> {
    match b {
        Some(b) => {
            let mut y = vec![T::ZERO; x.len()];
            b.spmv(x, &mut y);
            y
        }
        None => x.to_vec(),
    }
}

/// `B`-orthonormalises the column `v` (with its `B`-image `bv`) against the sets `against` (each a list of `(vector, B·vector)`)
/// twice; returns `false` if the remainder is negligible.
fn orth_column<T: Scalar>(v: &mut [T], bv: &mut [T], against: &[&Vec<(Vec<T>, Vec<T>)>]) -> bool {
    let n0 = norm2(v);
    if n0 == 0.0 {
        return false;
    }
    for _ in 0..2 {
        for set in against {
            for (u, bu) in set.iter() {
                let h = dotc(bu, v); // uᴴ B v
                for ((vi, bvi), (ui, bui)) in v.iter_mut().zip(bv.iter_mut()).zip(u.iter().zip(bu)) {
                    *vi -= h * *ui;
                    *bvi -= h * *bui;
                }
            }
        }
    }
    let nb = dotc(v, bv).re();
    if nb.is_nan() || nb <= 1e-20 * n0 * n0 {
        return false;
    }
    let s = 1.0 / nb.sqrt();
    v.iter_mut().for_each(|x| *x = x.scale(s));
    bv.iter_mut().for_each(|x| *x = x.scale(s));
    true
}

/// LOBPCG. `x0`: optional starting block (`n`-vectors; random if absent). `constraints`: vectors (not
/// necessarily orthonormal) to stay `B`-orthogonal to. `pc` approximates `(A − σB)⁻¹` for a shift below the wanted part of the
/// spectrum (identity if `None`).
pub fn lobpcg<T: FaerScalar>(
    a: &Csr<T>,
    b: Option<&Csr<T>>,
    pc: Option<&dyn Preconditioner<T>>,
    x0: Option<&[Vec<T>]>,
    constraints: &[Vec<T>],
    opts: &LobpcgOptions,
) -> Result<LobpcgResult<T>, LobpcgError> {
    let n = a.nrows();
    let m = opts.block.max(opts.nev);
    if opts.nev == 0 || 3 * m + constraints.len() >= n || a.ncols() != n {
        return Err(LobpcgError::InvalidArgument("need nev ≥ 1, a square matrix and 3·block + constraints < n"));
    }
    let sign = if opts.which == LobpcgWhich::Lowest { 1.0 } else { -1.0 };
    // constraint basis (B-orthonormal)
    let mut cons: Vec<(Vec<T>, Vec<T>)> = Vec::new();
    for y in constraints {
        let (mut v, mut bv) = (y.clone(), apply_b(b, y));
        if orth_column(&mut v, &mut bv, &[&cons]) {
            cons.push((v, bv));
        }
    }
    // starting block
    let mut x: Vec<(Vec<T>, Vec<T>)> = Vec::new();
    let mut seed = opts.seed;
    for j in 0..m {
        loop {
            let mut v = match x0.and_then(|s| s.get(j)) {
                Some(s) if seed == opts.seed || x.len() > j => s.clone(),
                _ => random_vector::<T>(n, seed),
            };
            seed = seed.wrapping_add(0x9E37);
            let mut bv = apply_b(b, &v);
            if orth_column(&mut v, &mut bv, &[&cons, &x]) {
                x.push((v, bv));
                break;
            }
        }
    }
    let mut ax: Vec<Vec<T>> = x
        .iter()
        .map(|(v, _)| {
            let mut y = vec![T::ZERO; n];
            a.spmv(v, &mut y);
            y
        })
        .collect();
    let mut matvecs = m;
    let ident: Box<dyn Preconditioner<T>>;
    let pcr: &dyn Preconditioner<T> = match pc {
        Some(p) => p,
        None => {
            ident = Box::new(crate::solvers::linop::Identity);
            &*ident
        }
    };
    let mut lambda = vec![0.0f64; m];
    // initial Rayleigh–Ritz in the starting block
    {
        let g = Dense::from_fn(m, m, |i, j| dotc(&x[i].0, &ax[j]));
        let g = hermitise(&g, sign);
        let (w, q) = hermitian_eig(&g)?;
        rotate(&mut x, &mut ax, &q, m);
        for i in 0..m {
            lambda[i] = sign * w[i];
        }
    }
    let mut p_dirs: Vec<(Vec<T>, Vec<T>, Vec<T>)> = Vec::new(); // (P, B·P, A·P) for the active columns
    let mut active_history = Vec::new();
    let (mut best_worst, mut since_improve) = (f64::INFINITY, 0usize);
    let mut iterations = 0;
    let mut reason = LobpcgReason::MaxIterations;
    let mut residual = vec![0.0f64; m];
    for it in 0..opts.max_it {
        // residuals and convergence of the block
        let mut r: Vec<Vec<T>> = Vec::with_capacity(m);
        for j in 0..m {
            let mut rj = ax[j].clone();
            for (ri, bi) in rj.iter_mut().zip(&x[j].1) {
                *ri -= bi.scale(lambda[j]);
            }
            let scale = norm2(&ax[j]) + lambda[j].abs() * norm2(&x[j].1);
            residual[j] = norm2(&rj) / scale.max(1e-300);
            r.push(rj);
        }
        let wanted_ok = residual.iter().take(opts.nev).all(|&e| e <= opts.tol);
        let active: Vec<usize> = (0..m).filter(|&j| residual[j] > opts.tol).collect();
        active_history.push(active.len());
        iterations = it;
        if wanted_ok {
            reason = LobpcgReason::Converged;
            break;
        }
        let worst = active.iter().filter(|&&j| j < opts.nev).map(|&j| residual[j]).fold(0.0, f64::max);
        if worst < (1.0 - opts.stall_gain) * best_worst {
            best_worst = worst;
            since_improve = 0;
        } else {
            since_improve += 1;
            if since_improve >= opts.stall_window {
                reason = LobpcgReason::Stalled;
                break;
            }
        }
        // new directions: W = T R (active), orthonormalised against constraints, X and each other
        let mut w_dirs: Vec<(Vec<T>, Vec<T>, Vec<T>)> = Vec::new();
        let mut w_basis: Vec<(Vec<T>, Vec<T>)> = Vec::new();
        for &j in &active {
            let mut wj = vec![T::ZERO; n];
            pcr.apply(&r[j], &mut wj);
            let mut bw = apply_b(b, &wj);
            if orth_column(&mut wj, &mut bw, &[&cons, &x, &w_basis]) {
                let mut aw = vec![T::ZERO; n];
                a.spmv(&wj, &mut aw);
                matvecs += 1;
                w_basis.push((wj.clone(), bw.clone()));
                w_dirs.push((wj, bw, aw));
            }
        }
        // previous directions, orthonormalised against everything else
        let mut p_new: Vec<(Vec<T>, Vec<T>, Vec<T>)> = Vec::new();
        let mut p_basis: Vec<(Vec<T>, Vec<T>)> = Vec::new();
        for (pv, _, _) in p_dirs.iter() {
            let mut v = pv.clone();
            let mut bv = apply_b(b, &v);
            if orth_column(&mut v, &mut bv, &[&cons, &x, &w_basis, &p_basis]) {
                let mut av = vec![T::ZERO; n];
                a.spmv(&v, &mut av);
                matvecs += 1;
                p_basis.push((v.clone(), bv.clone()));
                p_new.push((v, bv, av));
            }
        }
        // Rayleigh–Ritz on S = [X, W, P]
        let mut s_vecs: Vec<&Vec<T>> = x.iter().map(|e| &e.0).collect();
        let mut as_vecs: Vec<&Vec<T>> = ax.iter().collect();
        for (v, _, av) in w_dirs.iter().chain(p_new.iter()) {
            s_vecs.push(v);
            as_vecs.push(av);
        }
        let d = s_vecs.len();
        let g = Dense::from_fn(d, d, |i, j| dotc(s_vecs[i], as_vecs[j]));
        let g = hermitise(&g, sign);
        let (w, q) = hermitian_eig(&g)?;
        // new X = S Q[:, :m], new P = (W,P part) Q[:, :m]
        let cols: Vec<Vec<T>> = (0..m).map(|c| (0..d).map(|i| q[(i, c)]).collect()).collect();
        let combine = |vecs: &[&Vec<T>], coef: &[T], from: usize| -> Vec<T> {
            let mut out = vec![T::ZERO; n];
            for (k, v) in vecs.iter().enumerate().skip(from) {
                let c = coef[k];
                for (o, vi) in out.iter_mut().zip(v.iter()) {
                    *o += c * *vi;
                }
            }
            out
        };
        let bs_vecs: Vec<&Vec<T>> =
            x.iter().map(|e| &e.1).chain(w_dirs.iter().map(|e| &e.1)).chain(p_new.iter().map(|e| &e.1)).collect();
        let mut new_x = Vec::with_capacity(m);
        let mut new_ax = Vec::with_capacity(m);
        let mut new_p = Vec::new();
        for col in cols.iter().take(m) {
            let xv = combine(&s_vecs, col, 0);
            let bxv = combine(&bs_vecs, col, 0);
            let axv = combine(&as_vecs, col, 0);
            new_x.push((xv, bxv));
            new_ax.push(axv);
            if m < d {
                let pv = combine(&s_vecs, col, m);
                let bpv = combine(&bs_vecs, col, m);
                let apv = combine(&as_vecs, col, m);
                new_p.push((pv, bpv, apv));
            }
        }
        for i in 0..m {
            lambda[i] = sign * w[i];
        }
        x = new_x;
        ax = new_ax;
        // soft locking: only the directions of still-active columns are kept for the next iteration
        let act_after: Vec<usize> = active.clone();
        p_dirs = new_p.into_iter().enumerate().filter(|(j, _)| act_after.contains(j)).map(|(_, p)| p).collect();
        if p_dirs.is_empty() {
            p_dirs = Vec::new();
        }
    }
    // final residuals (recomputed from the final block)
    let nev = opts.nev;
    let mut vectors = Vec::with_capacity(nev);
    let mut residuals = Vec::with_capacity(nev);
    for j in 0..nev {
        let mut rj = ax[j].clone();
        for (ri, bi) in rj.iter_mut().zip(&x[j].1) {
            *ri -= bi.scale(lambda[j]);
        }
        residuals.push(norm2(&rj) / (norm2(&ax[j]) + lambda[j].abs() * norm2(&x[j].1)).max(1e-300));
        vectors.push(x[j].0.clone());
    }
    Ok(LobpcgResult { values: lambda[..nev].to_vec(), vectors, residuals, iterations, reason, matvecs, active_history })
}

/// Hermitian part of `g`, with the sign flipped for "highest" so that the ascending eigenvalue order is the wanted order.
fn hermitise<T: FaerScalar>(g: &Dense<T>, sign: f64) -> Dense<T> {
    let d = g.rows();
    Dense::from_fn(d, d, |i, j| (g[(i, j)] + g[(j, i)].conj()).scale(0.5 * sign))
}

fn rotate<T: Scalar>(x: &mut [(Vec<T>, Vec<T>)], ax: &mut [Vec<T>], q: &Dense<T>, m: usize) {
    let n = x[0].0.len();
    fn mix<T: Scalar>(vecs: &[&Vec<T>], q: &Dense<T>, c: usize, n: usize) -> Vec<T> {
        let mut out = vec![T::ZERO; n];
        for (i, v) in vecs.iter().enumerate() {
            let coef = q[(i, c)];
            for (o, vi) in out.iter_mut().zip(v.iter()) {
                *o += coef * *vi;
            }
        }
        out
    }
    let xs: Vec<&Vec<T>> = x.iter().map(|e| &e.0).collect();
    let bxs: Vec<&Vec<T>> = x.iter().map(|e| &e.1).collect();
    let axs: Vec<&Vec<T>> = ax.iter().collect();
    let new_x: Vec<Vec<T>> = (0..m).map(|c| mix(&xs, q, c, n)).collect();
    let new_bx: Vec<Vec<T>> = (0..m).map(|c| mix(&bxs, q, c, n)).collect();
    let new_ax: Vec<Vec<T>> = (0..m).map(|c| mix(&axs, q, c, n)).collect();
    for c in 0..m {
        x[c] = (new_x[c].clone(), new_bx[c].clone());
        ax[c] = new_ax[c].clone();
    }
}
