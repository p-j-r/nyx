//! Block Gram–Schmidt for the Krylov bases, with an optional `B`-weighted inner product
//! `⟨x, y⟩ = xᴴ B y` (generalised Hermitian problems).
//!
//! One call orthogonalises a block `W` of `b` new vectors against the first `k` basis columns
//! (classical Gram–Schmidt with a DGKS refinement pass when more than `1 − 1/√2` of the norm
//! cancels, the SLEPc default), then orthonormalises the columns of `W` among themselves
//! (column-wise CGS2), returning the coefficients `H` (`k × b`) and the triangular factor `R`
//! (`b × b`) of the Krylov decomposition `A·X = [V W]·[H; R]`.
//!
//! A column whose residual norm is below `100·ε` times its original norm is a *breakdown*
//! (an invariant subspace was found): it is replaced by a random vector orthogonal to
//! everything and its diagonal entry in `R` is set to exactly zero, so the decomposition stays
//! valid and the iteration continues without a restart from scratch.

use crate::solvers::basis::BasisView;
use crate::solvers::linop::LinearOperator;
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::{axpy, dotc, norm2_sq, scale_real};
use std::time::{Duration, Instant};

/// DGKS threshold on the *squared* norm ratio (`η² = 1/2`).
const ETA_SQ: f64 = 0.5;

pub(crate) struct OrthWork<T: Scalar> {
    n: usize,
    b: usize,
    /// `B·W` for the block just orthogonalised (weighted problems only).
    pub u: Vec<T>,
    /// `B·X` for the block that is about to be expanded (weighted problems only).
    pub bx: Vec<T>,
    /// Coefficients `H`, `k × b` column-major.
    pub h: Vec<T>,
    h2: Vec<T>,
    n0: Vec<f64>,
    n1: Vec<f64>,
    /// Triangular factor `R`, `b × b` column-major.
    pub r: Vec<T>,
    seed: u64,
    pub breakdowns: usize,
    /// Vectors multiplied by `B`.
    pub b_products: usize,
    /// Block orthogonalisations performed / how many needed the refinement pass.
    pub calls: usize,
    pub refinements: usize,
    pub b_time: Duration,
}

impl<T: Scalar> OrthWork<T> {
    pub fn new(n: usize, max_k: usize, b: usize, weighted: bool, seed: u64) -> Self {
        Self {
            n,
            b,
            u: if weighted { vec![T::ZERO; n * b] } else { Vec::new() },
            bx: if weighted { vec![T::ZERO; n * b] } else { Vec::new() },
            h: vec![T::ZERO; max_k * b],
            h2: vec![T::ZERO; max_k * b],
            n0: vec![0.0; b],
            n1: vec![0.0; b],
            r: vec![T::ZERO; b * b],
            seed,
            breakdowns: 0,
            b_products: 0,
            calls: 0,
            refinements: 0,
            b_time: Duration::ZERO,
        }
    }

    fn next_seed(&mut self) -> u64 {
        self.seed = self.seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        self.seed
    }
}

fn apply_b<T: Scalar>(
    bop: &dyn LinearOperator<T>,
    w: &[T],
    u: &mut [T],
    b: usize,
    wk_products: &mut usize,
    wk_time: &mut Duration,
) {
    let t0 = Instant::now();
    bop.apply_block(w, u, b);
    *wk_products += b;
    *wk_time += t0.elapsed();
}

/// Orthogonalises the `b` columns of `w` against `head[:, :k]` and among themselves. On return
/// `wk.h` holds `H` (`k × b`), `wk.r` holds `R`, and, if weighted, `wk.u = B·w`.
pub(crate) fn orthogonalize<T: Scalar>(
    head: BasisView<'_, T>,
    k: usize,
    w: &mut [T],
    b: usize,
    weight: Option<&dyn LinearOperator<T>>,
    wk: &mut OrthWork<T>,
) {
    let n = wk.n;
    assert!(b == wk.b && w.len() >= n * b && head.len() == n);
    let w = &mut w[..n * b];
    wk.calls += 1;

    match weight {
        None => head.dots_block(k, w, b, &mut wk.h, Some(&mut wk.n0)),
        Some(bop) => {
            apply_b(bop, w, &mut wk.u, b, &mut wk.b_products, &mut wk.b_time);
            head.dots_block(k, &wk.u, b, &mut wk.h, None);
            for c in 0..b {
                wk.n0[c] = dotc(&w[c * n..(c + 1) * n], &wk.u[c * n..(c + 1) * n]).re();
            }
        }
    }
    match weight {
        None => head.sub_block_then_dots(k, &wk.h, w, b, &mut wk.h2, &mut wk.n1),
        Some(bop) => {
            head.sub_block(k, &wk.h, w, b);
            apply_b(bop, w, &mut wk.u, b, &mut wk.b_products, &mut wk.b_time);
            head.dots_block(k, &wk.u, b, &mut wk.h2, None);
            for c in 0..b {
                wk.n1[c] = dotc(&w[c * n..(c + 1) * n], &wk.u[c * n..(c + 1) * n]).re();
            }
        }
    }
    if (0..b).any(|c| wk.n1[c] < ETA_SQ * wk.n0[c]) {
        wk.refinements += 1;
        head.sub_block(k, &wk.h2, w, b);
        for (a, c) in wk.h[..k * b].iter_mut().zip(&wk.h2[..k * b]) {
            *a += *c;
        }
        if let Some(bop) = weight {
            apply_b(bop, w, &mut wk.u, b, &mut wk.b_products, &mut wk.b_time);
        }
    }
    intra_block(head, k, w, b, weight, wk);
}

/// Column-wise CGS2 inside the block; fills `wk.r`.
fn intra_block<T: Scalar>(
    head: BasisView<'_, T>,
    k: usize,
    w: &mut [T],
    b: usize,
    weight: Option<&dyn LinearOperator<T>>,
    wk: &mut OrthWork<T>,
) {
    let n = wk.n;
    let weighted = weight.is_some();
    wk.r.fill(T::ZERO);
    // Take `u` out of the work struct so its column slices can be held while `wk` is borrowed
    // mutably by the breakdown path.
    let mut u_all = std::mem::take(&mut wk.u);
    for c in 0..b {
        let scale0 = wk.n0[c].max(0.0).sqrt();
        let (prev, cur) = w.split_at_mut(c * n);
        let cur = &mut cur[..n];
        let (uprev, ucur): (&[T], &mut [T]) = if weighted {
            let (p, q) = u_all.split_at_mut(c * n);
            (&*p, &mut q[..n])
        } else {
            (&[], &mut [])
        };
        if c > 0 {
            for _pass in 0..2 {
                for i in 0..c {
                    let wi = &prev[i * n..(i + 1) * n];
                    let g = if weighted { dotc(wi, ucur) } else { dotc(wi, cur) };
                    axpy(-g, wi, cur);
                    if weighted {
                        axpy(-g, &uprev[i * n..(i + 1) * n], ucur);
                    }
                    wk.r[i + c * b] += g;
                }
            }
        }
        let nsq = if weighted { dotc(cur, ucur).re() } else { norm2_sq(cur) };
        let nrm = nsq.max(0.0).sqrt();
        if nrm.is_nan() || nrm <= 100.0 * f64::EPSILON * scale0 || !nrm.is_finite() {
            // Breakdown: continue with a fresh vector orthogonal to everything so far.
            wk.breakdowns += 1;
            let seed = wk.next_seed();
            replace_column(head, k, prev, c, cur, ucur, weight, seed, wk);
            wk.r[c + c * b] = T::ZERO;
        } else {
            wk.r[c + c * b] = T::from_f64(nrm);
            scale_real(1.0 / nrm, cur);
            if weighted {
                scale_real(1.0 / nrm, ucur);
            }
        }
    }
    wk.u = u_all;
}

/// Fills `x` with a unit vector orthogonal to `head[:, :k]` and to the `c` already-final block
/// columns `prev`; sets `u = B·x` if weighted.
#[allow(clippy::too_many_arguments)]
fn replace_column<T: Scalar>(
    head: BasisView<'_, T>,
    k: usize,
    prev: &[T],
    c: usize,
    x: &mut [T],
    u: &mut [T],
    weight: Option<&dyn LinearOperator<T>>,
    mut seed: u64,
    wk: &mut OrthWork<T>,
) {
    let n = wk.n;
    let mut h = vec![T::ZERO; k.max(1)];
    for attempt in 0..8u64 {
        seed = seed.wrapping_add(attempt);
        x.copy_from_slice(&random_vector::<T>(n, seed));
        let mut before = 0.0;
        for pass in 0..3 {
            if let Some(bop) = weight {
                apply_b(bop, x, u, 1, &mut wk.b_products, &mut wk.b_time);
            }
            // Classical Gram–Schmidt: all coefficients come from the same `ux = B·x`.
            let gs: Vec<T> = {
                let ux: &[T] = if weight.is_some() { &*u } else { &*x };
                if pass == 0 {
                    before = dotc(x, ux).re().max(0.0).sqrt();
                }
                if k > 0 {
                    head.dots(k, ux, &mut h[..k]);
                }
                (0..c).map(|i| dotc(&prev[i * n..(i + 1) * n], ux)).collect()
            };
            if k > 0 {
                head.sub_combination(k, &h[..k], x);
            }
            for (i, g) in gs.iter().enumerate() {
                axpy(-*g, &prev[i * n..(i + 1) * n], x);
            }
        }
        let nsq = if let Some(bop) = weight {
            apply_b(bop, x, u, 1, &mut wk.b_products, &mut wk.b_time);
            dotc(x, u).re()
        } else {
            norm2_sq(x)
        };
        let nrm = nsq.max(0.0).sqrt();
        if nrm > 1e-8 * before && nrm.is_finite() && nrm > 0.0 {
            scale_real(1.0 / nrm, x);
            if weight.is_some() {
                scale_real(1.0 / nrm, u);
            }
            return;
        }
    }
    // The space is exhausted (only possible when the basis already spans everything).
    x.fill(T::ZERO);
    if weight.is_some() {
        u.fill(T::ZERO);
    }
}
