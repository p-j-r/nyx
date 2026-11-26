use super::orth::{OrthWork, orthogonalize};
use super::{EigenResult, EpsError, EpsOptions, EpsReason, HistoryEntry, Kind, Timings};
use crate::linalg::basis::Basis;
use crate::linalg::dense::{Dense, hermitian_eig, reorder, schur, triangular_eigenvector};
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::LinearOperator;
use crate::linalg::rng::random_vector;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{norm2, scale_real};
use std::time::Instant;

/// Ritz data of the active window, in "most wanted first" order.
struct Ritz<T: Scalar> {
    theta: Vec<C64>,
    /// Residual norm estimates `‖F·Q·y‖/‖y‖` (`y` = unit vector or triangular eigenvector).
    rho: Vec<f64>,
    /// Orthogonal/unitary transformation (`m × m`, columns in wanted order).
    q: Dense<T>,
    /// Upper triangular Schur factor in wanted order (general kind only).
    s: Option<Dense<C64>>,
}

/// Thick-restart Krylov–Schur with locking for `op·x = λ x`.
///
/// * `weight`: Hermitian positive (semi)definite `B`; when given, the basis is orthonormal in
///   `⟨x,y⟩ = xᴴBy` and `op` must be self-adjoint in it (e.g. `(A−σB)⁻¹B`). Use `None` for the
///   Euclidean inner product.
/// * `start`: first start vector (the other block columns are pseudo-random); `None` uses a
///   fixed pseudo-random vector.
///
/// Follows SLEPc's `EPSKRYLOVSCHUR`: restart keeps `max(1, ⌊(nv−k)·keep⌋)` Ritz vectors beyond
/// the `k` converged ones, convergence is the leading run of Ritz pairs (in wanted order) with
/// `‖r‖ < tol·|θ|`, and defaults for `ncv`/`max_it` are SLEPc's. For `block_size = 1` and the
/// Hermitian kind this reproduces SLEPc's restart and matvec counts on the problems tested (see
/// `docs/STATUS.md`).
///
/// Residual estimates for [`Kind::General`] use the eigenvector of the triangular factor,
/// ignoring the (locked-part) components of the eigenvector; callers that need guaranteed
/// accuracy must recompute `‖Ax − λBx‖` themselves, which [`crate::linalg::problem`] does.
pub fn krylov_schur<T, A>(
    op: &A,
    weight: Option<&dyn LinearOperator<T>>,
    start: Option<&[T]>,
    opts: &EpsOptions,
) -> Result<EigenResult<T>, EpsError>
where
    T: Scalar + FaerScalar,
    A: LinearOperator<T> + ?Sized,
{
    let t_start = Instant::now();
    let n = op.dim();
    let nev = opts.nev;
    let b = opts.block_size;
    if nev == 0 || nev > n {
        return Err(EpsError::InvalidArgument("nev must be in 1..=n"));
    }
    if b == 0 {
        return Err(EpsError::InvalidArgument("block_size must be at least 1"));
    }
    if opts.tol.is_nan() || opts.tol <= 0.0 {
        return Err(EpsError::InvalidArgument("tol must be positive"));
    }
    if !(0.1..=0.9).contains(&opts.keep) {
        return Err(EpsError::InvalidArgument("keep must be in [0.1, 0.9]"));
    }
    if opts.kind == Kind::General && !T::IS_COMPLEX {
        return Err(EpsError::Unsupported("the general (non-Hermitian) kind needs a complex scalar"));
    }
    if weight.is_some_and(|w| w.dim() != n) {
        return Err(EpsError::InvalidArgument("weight matrix dimension differs from the operator"));
    }
    if start.is_some_and(|s| s.is_empty() || s.len() % n != 0) {
        return Err(EpsError::InvalidArgument("start vector has wrong length"));
    }
    // The decomposition needs `ncv + b` basis columns, which must fit in the space.
    let default_ncv = (2 * nev).max(nev + 15).max(nev + 2 * b);
    let ncv = match opts.ncv {
        Some(c) => c,
        None => default_ncv,
    }
    .min(n.saturating_sub(b));
    if ncv < nev + 1 || ncv < b + 1 {
        return Err(EpsError::InvalidArgument(
            "problem too small for this nev/block size (n >= nev + block_size + 1 is needed); use a dense solver",
        ));
    }
    let max_it = opts.max_it.unwrap_or_else(|| 100.max(2 * n / ncv));
    let weighted = weight.is_some();
    let ld = ncv + b;

    let mut v = Basis::<T>::new(n, ncv + b);
    let mut t = Dense::<T>::zeros(ld, ncv);
    let mut lam = vec![C64::ZERO; ncv];
    let mut errs = vec![f64::INFINITY; ncv];
    let mut work = OrthWork::<T>::new(n, ncv + b, b, weighted, opts.seed ^ 0xB0B0_B0B0);
    let mut tm = Timings::default();
    let mut history = Vec::new();
    let mut matvecs = 0usize;

    // Start block: the user's columns first (one vector, or the several columns of a warm start, at most `b`), pseudo-random fill, then
    // B-orthonormalise. A single start vector, or none, gives the block it always did.
    {
        let cols = v.cols_mut(0, b);
        let given = start.map_or(0, |s| (s.len() / n).min(b));
        match start {
            Some(s) => cols[..n * given].copy_from_slice(&s[..n * given]),
            None => cols[..n].copy_from_slice(&random_vector::<T>(n, opts.seed)),
        }
        for c in given.max(1)..b {
            cols[c * n..(c + 1) * n].copy_from_slice(&random_vector::<T>(n, opts.seed.wrapping_add(0x1000 + c as u64)));
        }
        let nrm = norm2(&cols[..n]);
        if nrm.is_nan() || nrm <= 0.0 || !nrm.is_finite() {
            return Err(EpsError::InvalidArgument("start vector must be non-zero and finite"));
        }
        let in_range = opts.start_in_range.unwrap_or(weighted && opts.kind == Kind::Hermitian);
        if in_range {
            // Force the start block into the range of the operator (see `EpsOptions`).
            let mut tmp = vec![T::ZERO; n * b];
            let t0 = Instant::now();
            op.apply_block(&cols[..n * b], &mut tmp, b);
            tm.operator += t0.elapsed().as_secs_f64();
            matvecs += b;
            cols[..n * b].copy_from_slice(&tmp);
            let nr = norm2(&cols[..n]);
            if nr.is_nan() || nr <= 0.0 {
                return Err(EpsError::InvalidArgument("start vector lies in the null space of the operator"));
            }
        }
        let t0 = Instant::now();
        let (head, w) = v.split_cols_mut(0);
        orthogonalize(head, 0, &mut w[..n * b], b, weight, &mut work);
        tm.orthogonalize += t0.elapsed().as_secs_f64();
        if weighted {
            std::mem::swap(&mut work.u, &mut work.bx);
        }
    }

    let (mut nconv, mut kept) = (0usize, 0usize);
    let mut its = 0usize;
    let reason;
    loop {
        its += 1;

        // ---- Expand the Krylov–Schur decomposition from `kept` columns to `nv`.
        let mut j0 = kept;
        while j0 + b <= ncv {
            let k = j0 + b;
            {
                let (head, tail) = v.split_cols_mut(k);
                let y = &mut tail[..n * b];
                let t0 = Instant::now();
                if weighted {
                    op.apply_block_with_bx(head.cols(j0, k), &work.bx, y, b);
                } else {
                    op.apply_block(head.cols(j0, k), y, b);
                }
                tm.operator += t0.elapsed().as_secs_f64();
                matvecs += b;
                let (t1, bt0) = (Instant::now(), work.b_time);
                orthogonalize(head, k, y, b, weight, &mut work);
                let dt = t1.elapsed().as_secs_f64() - (work.b_time - bt0).as_secs_f64();
                tm.orthogonalize += dt;
            }
            for c in 0..b {
                for i in 0..k {
                    t[(i, j0 + c)] = work.h[i + c * k];
                }
                for i in 0..=c {
                    t[(k + i, j0 + c)] = work.r[i + c * b];
                }
            }
            if weighted {
                std::mem::swap(&mut work.u, &mut work.bx);
            }
            j0 += b;
        }
        let nv = j0;

        // ---- Projected problem on the active window [nconv, nv).
        let m = nv - nconv;
        let t0 = Instant::now();
        let ritz = solve_projected(&t, nconv, nv, b, opts)?;
        tm.dense += t0.elapsed().as_secs_f64();

        let mut k = nconv;
        let mut run = true;
        let (mut cycle_errs, mut cycle_ritz) = (Vec::new(), Vec::new());
        for p in 0..m {
            let (theta, rho) = (ritz.theta[p], ritz.rho[p]);
            let err = if rho == 0.0 {
                0.0
            } else if theta.norm() > 0.0 {
                match opts.original_shift {
                    None => rho / theta.norm(),
                    Some(sigma) => {
                        let lam_norm = (sigma + C64::ONE / theta).norm();
                        // Never relax stopping below transformed residual `rho / |theta|` when |lambda| > 1
                        // (e.g. SI-unit FEM pencils where omega^2 or k0^2 >> 1); only tighten when |lambda| < 1.
                        let denom = theta.norm() * lam_norm.min(1.0);
                        if denom > 0.0 { rho / denom } else { f64::MAX }
                    }
                }
            } else {
                f64::MAX
            };
            lam[nconv + p] = theta;
            errs[nconv + p] = err;
            if run && err < opts.tol {
                k += 1;
            } else {
                run = false;
            }
            if opts.record_history && p < nev + 2 {
                cycle_errs.push(err);
                cycle_ritz.push(theta);
            }
        }
        if opts.record_history {
            history.push(HistoryEntry { cycle: its, nconv: k, nv, matvecs, errors: cycle_errs, ritz: cycle_ritz });
        }

        let stop = if k >= nev {
            Some(EpsReason::Converged)
        } else if its >= max_it {
            Some(EpsReason::DivergedIts)
        } else {
            None
        };
        let l = if stop.is_some() || k == nv {
            0
        } else {
            // SLEPc: l = max(1, (nv-k)*keep); never leave less than one block of room.
            1.max(((nv - k) as f64 * opts.keep) as usize).min(ncv - k - b)
        };

        // ---- Restart: V[:, nconv..] ← V[:, nconv..nv]·Q, rebuild the leading part of T.
        let t0 = Instant::now();
        let cols = k - nconv + l;
        let new_kept = nconv + cols;
        // Coupling rows F·Q (b × cols) and, for the general kind, the coupling of the
        // locked part to the new Schur vectors.
        let q_sel = &ritz.q.as_slice()[..m * cols];
        let mut fq = Dense::<T>::zeros(b, cols);
        for c in 0..cols {
            for r in 0..b {
                let mut s = T::ZERO;
                for i in 0..m {
                    s += t[(nv + r, nconv + i)] * q_sel[c * m + i];
                }
                fq[(r, c)] = s;
            }
        }
        let mut locked_coupling = Dense::<T>::zeros(nconv, cols);
        if opts.kind == Kind::General {
            for c in 0..cols {
                for r in 0..nconv {
                    let mut s = T::ZERO;
                    for i in 0..m {
                        s += t[(r, nconv + i)] * q_sel[c * m + i];
                    }
                    locked_coupling[(r, c)] = s;
                }
            }
        }
        if cols > 0 {
            v.mult_in_place(nconv, m, cols, q_sel, m);
        }
        // Clear everything from column `nconv` on, then write the new structure.
        for c in nconv..ncv {
            t.col_mut(c).fill(T::ZERO);
        }
        for c in 0..cols {
            let th = ritz.theta[c];
            match &ritz.s {
                None => t[(nconv + c, nconv + c)] = T::from_f64(th.re),
                Some(s) => {
                    for r in 0..=c {
                        t[(nconv + r, nconv + c)] = T::from_parts(s[(r, c)].re, s[(r, c)].im);
                    }
                    for r in 0..nconv {
                        t[(r, nconv + c)] = locked_coupling[(r, c)];
                    }
                }
            }
            // Residual coupling of the pairs that stay active; locked pairs count as exact.
            if nconv + c >= k {
                for r in 0..b {
                    t[(new_kept + r, nconv + c)] = fq[(r, c)];
                }
            }
        }
        lam[nconv..nconv + cols].copy_from_slice(&ritz.theta[..cols]);
        nconv = k;
        if stop.is_none() {
            v.move_cols(nv, new_kept, b);
            kept = new_kept;
        }
        tm.restart += t0.elapsed().as_secs_f64();

        if let Some(r) = stop {
            reason = r;
            break;
        }
    }

    // ---- Collect results: eigenvectors, ordered from most to least wanted.
    let t0 = Instant::now();
    if opts.kind == Kind::General && nconv > 0 {
        // Schur vectors → eigenvectors: V·Y with Y the eigenvector matrix of the (locked)
        // upper triangular block.
        let tc = Dense::<C64>::from_fn(nconv, nconv, |i, j| if i <= j { t[(i, j)].to_c64() } else { C64::ZERO });
        let mut y = Dense::<T>::zeros(nconv, nconv);
        for p in 0..nconv {
            let col = triangular_eigenvector(&tc, nconv, p);
            for (i, c) in col.iter().enumerate() {
                y[(i, p)] = T::from_parts(c.re, c.im);
            }
        }
        v.mult_in_place(0, nconv, nconv, y.as_slice(), nconv);
        for p in 0..nconv {
            let nrm = norm2(v.col(p));
            scale_real(1.0 / nrm, v.col_mut(p));
        }
    }
    let mut idx: Vec<usize> = (0..nconv).collect();
    idx.sort_by(|&x, &y| opts.which.key(lam[x]).total_cmp(&opts.which.key(lam[y])));
    let mut vectors = Basis::<T>::new(n, nconv.max(1));
    for (dst, &src) in idx.iter().enumerate() {
        vectors.col_mut(dst).copy_from_slice(v.col(src));
    }
    tm.restart += t0.elapsed().as_secs_f64();
    tm.weight = work.b_time.as_secs_f64();
    tm.total = t_start.elapsed().as_secs_f64();
    Ok(EigenResult {
        reason,
        values: idx.iter().map(|&i| lam[i]).collect(),
        errors: idx.iter().map(|&i| errs[i]).collect(),
        vectors,
        nconv,
        ncv,
        iterations: its,
        matvecs,
        b_products: work.b_products,
        breakdowns: work.breakdowns,
        orth_calls: work.calls,
        orth_refinements: work.refinements,
        timings: tm,
        history,
    })
}

/// Solves the projected problem of the active window `[nconv, nv)` and orders the Ritz pairs.
fn solve_projected<T: Scalar + FaerScalar>(
    t: &Dense<T>,
    nconv: usize,
    nv: usize,
    b: usize,
    opts: &EpsOptions,
) -> Result<Ritz<T>, EpsError> {
    let m = nv - nconv;
    let f = |r: usize, i: usize| t[(nv + r, nconv + i)];
    match opts.kind {
        Kind::Hermitian => {
            // Use the computed upper triangle and mirror it (the lower part holds only the
            // restart coupling).
            let tact = Dense::<T>::from_fn(m, m, |i, j| {
                if i < j {
                    t[(nconv + i, nconv + j)]
                } else if i > j {
                    t[(nconv + j, nconv + i)].conj()
                } else {
                    T::from_f64(t[(nconv + i, nconv + i)].re())
                }
            });
            let (w, q) = hermitian_eig(&tact)?;
            let mut order: Vec<usize> = (0..m).collect();
            order.sort_by(|&x, &y| opts.which.key(C64::new(w[x], 0.0)).total_cmp(&opts.which.key(C64::new(w[y], 0.0))));
            let qs = Dense::<T>::from_fn(m, m, |i, p| q[(i, order[p])]);
            let mut rho = vec![0.0; m];
            for p in 0..m {
                let mut s = 0.0;
                for r in 0..b {
                    let mut z = T::ZERO;
                    for i in 0..m {
                        z += f(r, i) * qs[(i, p)];
                    }
                    s += z.modulus_sq();
                }
                rho[p] = s.sqrt();
            }
            Ok(Ritz { theta: order.iter().map(|&i| C64::new(w[i], 0.0)).collect(), rho, q: qs, s: None })
        }
        Kind::General => {
            let tact = Dense::<C64>::from_fn(m, m, |i, j| t[(nconv + i, nconv + j)].to_c64());
            let (mut s, mut z) = schur(&tact)?;
            reorder(&mut s, &mut z, m, |c| opts.which.key(c));
            // F·Z (b × m)
            let fz =
                Dense::<C64>::from_fn(b, m, |r, c| (0..m).fold(C64::ZERO, |a, i| a + f(r, i).to_c64() * z[(i, c)]));
            let mut rho = vec![0.0; m];
            for (p, rho_p) in rho.iter_mut().enumerate() {
                let y = triangular_eigenvector(&s, m, p);
                let ynorm = y.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
                let mut acc = 0.0;
                for r in 0..b {
                    let val = (0..=p).fold(C64::ZERO, |a, i| a + fz[(r, i)] * y[i]);
                    acc += val.norm_sqr();
                }
                *rho_p = acc.sqrt() / ynorm;
            }
            let theta = (0..m).map(|p| s[(p, p)]).collect();
            let q = Dense::<T>::from_fn(m, m, |i, j| T::from_parts(z[(i, j)].re, z[(i, j)].im));
            Ok(Ritz { theta, rho, q, s: Some(s) })
        }
    }
}
