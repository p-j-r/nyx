//! **GCRO-DR** (Parks, de Sturler, Mackey, Johnson & Maiti 2006): restarted GMRES that *recycles* an
//! approximate invariant subspace, for the `A x_i = b_i` sequences that a frequency sweep, a Newton
//! iteration or a time stepper produces.
//!
//! Why: restarted GMRES(m) throws its Krylov space away every `m` steps, and with it the knowledge of
//! the few small (near-zero) eigenvalues that make the convergence stall. GCRO-DR keeps `k` harmonic
//! Ritz vectors `U` (those belonging to the harmonic Ritz values of smallest magnitude) together with
//! `C = A U`, `CᴴC = I`, and
//!
//! * projects every residual onto the complement of `C` (`x += U Cᴴ r`, `r -= C Cᴴ r`), so those
//!   eigen-directions are *solved exactly*, and
//! * runs the Arnoldi process on `(I − CCᴴ)A`, which never sees them again.
//!
//! The recycled space survives between calls of [`GcroDr::solve`], so a second right-hand side (or a
//! slightly different matrix: call [`GcroDr::matrix_changed`] first) starts with the deflation of the
//! first one. Per cycle the algorithm is the one of the paper (Algorithm 2) with these choices:
//!
//! * **Right preconditioning** (`A M⁻¹`): the monitored norm is the true residual `‖b − Ax‖`, and `U`
//!   lives in the preconditioned space, so one `M⁻¹` per cycle is needed for the solution update.
//! * **LS problem**. With `r ⟂ C` the least-squares problem `min ‖[0; β e₁] − Ḡ y‖`, `Ḡ = [D B; 0 H̄]`
//!   decouples: the `D` block (diagonal, `1/‖u_i‖`) can always be made to vanish, so the residual is that of
//!   the Hessenberg problem `min ‖β e₁ − H̄ y_v‖` alone (monitored with Givens rotations as in GMRES) and
//!   `x += M⁻¹(V y_v − U B y_v)`.
//! * **Harmonic Ritz vectors** from `ḠᴴḠ g = θ Ḡᴴ(V̂ᴴŴ) g` with `Ŵ = [U D, V]`, `V̂ = [C, V]`: reduced to a
//!   standard problem by a Cholesky factor of `ḠᴴḠ` and solved with the dense complex Schur routines of
//!   [`crate::solvers::dense`]. For real matrices complex pairs are replaced by their real and imaginary parts.
//! * The residual is recomputed from `b − Ax` at every cycle and projected again, which repairs
//!   rounding drift and any loss of orthogonality to `C`.
//!
//! Not done: the "GCROT"-style truncation and the optional modification that selects Ritz instead of
//! harmonic Ritz vectors (the paper's `GCRO-DR(m,k)` with harmonic vectors is the robust default).

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::solvers::dense::{Dense, schur, triangular_eigenvector};
use crate::solvers::linop::{LinearOperator, Preconditioner};
use crate::solvers::scalar::{C64, Scalar};
use crate::solvers::vecops::{axpy, dotc, norm2, scale_real};

#[derive(Debug, Clone)]
pub struct GcroDrOptions {
    pub tol: Tolerances,
    /// Cycle length `m`: recycled vectors plus Arnoldi steps per cycle (memory `n·(2k + m + 1)`).
    pub m: usize,
    /// Dimension `k < m` of the recycled subspace; `0` gives plain restarted GMRES(m) with right preconditioning.
    pub k: usize,
    pub initial_guess_nonzero: bool,
    pub record_history: bool,
}

impl Default for GcroDrOptions {
    fn default() -> Self {
        Self { tol: Tolerances::default(), m: 30, k: 10, initial_guess_nonzero: false, record_history: false }
    }
}

/// Counters of one solve beyond [`SolveInfo`].
#[derive(Debug, Clone, Default)]
pub struct GcroDrStats {
    /// Restart cycles.
    pub cycles: usize,
    /// Dimension of the recycled space when the solve started.
    pub recycled_at_start: usize,
    /// Dimension of the recycled space when the solve ended.
    pub recycled_at_end: usize,
}

pub struct GcroDr<T: Scalar> {
    pub opts: GcroDrOptions,
    u: Vec<Vec<T>>,
    c: Vec<Vec<T>>,
    stale: bool,
    stats: GcroDrStats,
}

/// Small helper: `Σ_i coef[i]·cols[i]` into `out` (zeroed first).
fn combine<T: Scalar>(cols: &[Vec<T>], coef: &[T], out: &mut [T]) {
    out.fill(T::ZERO);
    for (c, &a) in cols.iter().zip(coef) {
        if a != T::ZERO {
            axpy(a, c, out);
        }
    }
}

/// Plane rotation with real cosine that maps `(a, b)` to `(r, 0)`: `[c s; −s̄ c]`.
fn givens<T: Scalar>(a: T, b: T) -> (f64, T) {
    let (na, nb) = (a.modulus(), b.modulus());
    if nb == 0.0 {
        return (1.0, T::ZERO);
    }
    if na == 0.0 {
        return (0.0, b.conj().scale(1.0 / nb));
    }
    let norm = na.hypot(nb);
    (na / norm, (a.scale(1.0 / na) * b.conj()).scale(1.0 / norm))
}

impl<T: Scalar> GcroDr<T> {
    pub fn new(opts: GcroDrOptions) -> Self {
        assert!(opts.m >= 2, "m must be at least 2");
        assert!(opts.k < opts.m, "the recycled dimension k must be smaller than the cycle length m");
        Self { opts, u: Vec::new(), c: Vec::new(), stale: false, stats: GcroDrStats::default() }
    }

    /// Dimension of the subspace currently kept for the next solve.
    pub fn recycled_dim(&self) -> usize {
        self.u.len()
    }

    pub fn stats(&self) -> &GcroDrStats {
        &self.stats
    }

    /// Tells the solver that the matrix (or the preconditioner) changed since the last solve: `C = A U`
    /// is recomputed (and re-orthonormalised, `U` accordingly) at the start of the next one, `k` matrix
    /// products.
    pub fn matrix_changed(&mut self) {
        self.stale = true;
    }

    /// Diagnostics for tests: `(max |CᴴC − I|, max_i ‖Â u_i − c_i‖ / ‖c_i‖)` with `Â = A M⁻¹`. Both are rounding
    /// level for a healthy recycled space (the second only if the matrix did not change since it was built).
    pub fn recycle_defects<A: LinearOperator<T> + ?Sized, M: Preconditioner<T> + ?Sized>(
        &self,
        a: &A,
        m: &M,
    ) -> (f64, f64) {
        let n = a.dim();
        let (mut z, mut w) = (vec![T::ZERO; n], vec![T::ZERO; n]);
        let mut orth = 0.0f64;
        let mut link = 0.0f64;
        for (i, ci) in self.c.iter().enumerate() {
            for (j, cj) in self.c.iter().enumerate() {
                let d = dotc(ci, cj) - if i == j { T::ONE } else { T::ZERO };
                orth = orth.max(d.modulus());
            }
            m.apply(&self.u[i], &mut z);
            a.apply(&z, &mut w);
            axpy(-T::ONE, ci, &mut w);
            link = link.max(norm2(&w) / norm2(ci).max(f64::MIN_POSITIVE));
        }
        (orth, link)
    }

    /// Forgets the recycled space.
    pub fn reset(&mut self) {
        self.u.clear();
        self.c.clear();
        self.stale = false;
    }

    /// Applies `A M⁻¹` to `v`.
    fn op<A: LinearOperator<T> + ?Sized, M: Preconditioner<T> + ?Sized>(
        a: &A,
        m: &M,
        v: &[T],
        z: &mut [T],
        w: &mut [T],
        info: &mut SolveInfo,
    ) {
        m.apply(v, z);
        a.apply(z, w);
        info.matvecs += 1;
        info.pc_applies += 1;
    }

    /// `C = Â U` with a thin QR (modified Gram–Schmidt, `U ← U R⁻¹`), for a changed matrix.
    fn rebuild_c<A: LinearOperator<T> + ?Sized, M: Preconditioner<T> + ?Sized>(
        &mut self,
        a: &A,
        m: &M,
        info: &mut SolveInfo,
    ) {
        let n = a.dim();
        let (mut z, mut w) = (vec![T::ZERO; n], vec![T::ZERO; n]);
        let mut cs: Vec<Vec<T>> = Vec::new();
        let mut us: Vec<Vec<T>> = Vec::new();
        for u in &self.u {
            Self::op(a, m, u, &mut z, &mut w, info);
            let mut uu = u.clone();
            // orthogonalise against the accepted columns, carrying U along (so that C = Â U stays true)
            for _ in 0..2 {
                for (c, uc) in cs.iter().zip(&us) {
                    let h = dotc(c, &w);
                    axpy(-h, c, &mut w);
                    axpy(-h, uc, &mut uu);
                }
            }
            let nrm = norm2(&w);
            if nrm > 1e-12 * norm2(&uu).max(f64::MIN_POSITIVE) && nrm.is_finite() {
                scale_real(1.0 / nrm, &mut w);
                scale_real(1.0 / nrm, &mut uu);
                cs.push(w.clone());
                us.push(uu);
            }
        }
        self.c = cs;
        self.u = us;
        self.stale = false;
    }

    pub fn solve<A: LinearOperator<T> + ?Sized, M: Preconditioner<T> + ?Sized>(
        &mut self,
        a: &A,
        m: &M,
        b: &[T],
        x: &mut [T],
    ) -> SolveInfo {
        let n = a.dim();
        assert!(b.len() == n && x.len() == n, "dimension mismatch");
        let opts = self.opts.clone();
        let mut info = SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedIts,
            residual_norm: 0.0,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };
        self.stats = GcroDrStats::default();
        if !opts.initial_guess_nonzero {
            x.fill(T::ZERO);
        }
        if self.u.first().is_some_and(|u| u.len() != n) {
            self.reset();
        }
        if self.stale && !self.u.is_empty() {
            self.rebuild_c(a, m, &mut info);
        }
        self.stats.recycled_at_start = self.u.len();

        let mut r = vec![T::ZERO; n];
        let (mut z, mut w, mut tmp) = (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
        if opts.initial_guess_nonzero {
            a.apply(x, &mut r);
            info.matvecs += 1;
            for (ri, bi) in r.iter_mut().zip(b) {
                *ri = *bi - *ri;
            }
        } else {
            r.copy_from_slice(b);
        }
        let mut conv = ConvTest::new(opts.tol);
        let rn0 = norm2(&r);
        conv.set_reference(rn0);
        info.residual_norm = rn0;
        if opts.record_history {
            info.history.push(rn0);
        }
        if !rn0.is_finite() {
            info.reason = ConvergedReason::DivergedNanOrInf;
            return info;
        }
        if let Some(reason) = conv.check(rn0) {
            info.reason = reason;
            self.stats.recycled_at_end = self.u.len();
            return info;
        }

        let mut its = 0usize;
        'cycles: loop {
            self.stats.cycles += 1;
            // ---- project r onto the complement of C ---------------------------------------------
            let kc = self.c.len();
            if kc > 0 {
                for _ in 0..2 {
                    let coef: Vec<T> = self.c.iter().map(|c| dotc(c, &r)).collect();
                    combine(&self.u, &coef, &mut tmp);
                    m.apply(&tmp, &mut z);
                    info.pc_applies += 1;
                    axpy(T::ONE, &z, x);
                    for (c, &h) in self.c.iter().zip(&coef) {
                        axpy(-h, c, &mut r);
                    }
                }
            }
            let beta = norm2(&r);
            info.residual_norm = beta;
            if let Some(reason) = conv.check(beta) {
                info.reason = reason;
                break;
            }
            if its >= opts.tol.max_it {
                info.reason = ConvergedReason::DivergedIts;
                break;
            }
            let msteps = (opts.m - kc).max(1).min(opts.tol.max_it - its);

            // ---- Arnoldi on (I − CCᴴ) A M⁻¹ -------------------------------------------------------
            let mut v: Vec<Vec<T>> = Vec::with_capacity(msteps + 1);
            let mut v0 = r.clone();
            scale_real(1.0 / beta, &mut v0);
            v.push(v0);
            let mut bmat: Vec<Vec<T>> = Vec::new(); // bmat[j][i] = (Cᴴ A V)_{i,j}
            let mut hmat: Vec<Vec<T>> = Vec::new(); // rotated (triangular) columns of H̄, for the solve
            let mut hraw: Vec<Vec<T>> = Vec::new(); // un-rotated columns of H̄ (length j+2), for the recycling
            let mut cs = vec![0.0f64; msteps];
            let mut sn = vec![T::ZERO; msteps];
            let mut rs = vec![T::ZERO; msteps + 1];
            rs[0] = T::from_f64(beta);
            let mut done = 0usize;
            let mut breakdown = false;
            for j in 0..msteps {
                Self::op(a, m, &v[j], &mut z, &mut w, &mut info);
                its += 1;
                info.iterations = its;
                let norm_in = norm2(&w);
                let mut bcol = vec![T::ZERO; kc];
                for _ in 0..2 {
                    for (i, c) in self.c.iter().enumerate() {
                        let h = dotc(c, &w);
                        axpy(-h, c, &mut w);
                        bcol[i] += h;
                    }
                }
                let mut hcol = vec![T::ZERO; j + 2];
                for pass in 0..2 {
                    for (i, vi) in v.iter().enumerate() {
                        let h = dotc(vi, &w);
                        axpy(-h, vi, &mut w);
                        hcol[i] += h;
                    }
                    if pass == 0 && norm2(&w) > std::f64::consts::FRAC_1_SQRT_2 * norm_in {
                        break;
                    }
                }
                let hn = norm2(&w);
                hcol[j + 1] = T::from_f64(hn);
                bmat.push(bcol);
                if hn > 1e-14 * norm_in.max(f64::MIN_POSITIVE) && hn.is_finite() {
                    let mut vn = w.clone();
                    scale_real(1.0 / hn, &mut vn);
                    v.push(vn);
                } else {
                    breakdown = true;
                }
                hraw.push(hcol.clone());
                // apply previous rotations, then form a new one
                for i in 0..j {
                    let (hi, hi1) = (hcol[i], hcol[i + 1]);
                    hcol[i] = hi.scale(cs[i]) + sn[i] * hi1;
                    hcol[i + 1] = hi1.scale(cs[i]) - sn[i].conj() * hi;
                }
                let (c, s) = givens(hcol[j], hcol[j + 1]);
                cs[j] = c;
                sn[j] = s;
                let (a0, a1) = (hcol[j], hcol[j + 1]);
                hcol[j] = a0.scale(c) + s * a1;
                hcol[j + 1] = T::ZERO;
                let r0 = rs[j];
                rs[j] = r0.scale(c);
                rs[j + 1] = -(s.conj() * r0);
                hmat.push(hcol);
                done = j + 1;
                let res = rs[j + 1].modulus();
                if opts.record_history {
                    info.history.push(res);
                }
                if conv.check(res).is_some() || breakdown {
                    break;
                }
            }

            // ---- solution update: y from the rotated triangular system ------------------------------
            let mut y = vec![T::ZERO; done];
            for i in (0..done).rev() {
                let mut s = rs[i];
                for k in i + 1..done {
                    s -= hmat[k][i] * y[k];
                }
                y[i] = s / hmat[i][i];
            }
            // update direction  V y − U (B y)   (in the preconditioned space)
            tmp.fill(T::ZERO);
            for (j, &yj) in y.iter().enumerate() {
                axpy(yj, &v[j], &mut tmp);
            }
            if kc > 0 {
                let mut by = vec![T::ZERO; kc];
                for (j, &yj) in y.iter().enumerate() {
                    for (i, byi) in by.iter_mut().enumerate() {
                        *byi += bmat[j][i] * yj;
                    }
                }
                for (u, &h) in self.u.iter().zip(&by) {
                    axpy(-h, u, &mut tmp);
                }
            }
            m.apply(&tmp, &mut z);
            info.pc_applies += 1;
            axpy(T::ONE, &z, x);

            // ---- recycling: harmonic Ritz vectors of the cycle's Ḡ ----------------------------------
            if opts.k > 0 && !breakdown && done >= 1 {
                self.update_recycle(&v, &bmat, &hraw, done, kc);
            }

            // ---- true residual for the next cycle -----------------------------------------------------
            a.apply(x, &mut w);
            info.matvecs += 1;
            for ((ri, bi), wi) in r.iter_mut().zip(b).zip(&w) {
                *ri = *bi - *wi;
            }
            let rtrue = norm2(&r);
            if !rtrue.is_finite() {
                info.reason = ConvergedReason::DivergedNanOrInf;
                break 'cycles;
            }
            info.residual_norm = rtrue;
            if let Some(reason) = conv.check(rtrue) {
                info.reason = reason;
                break;
            }
            if breakdown {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
        }
        self.stats.recycled_at_end = self.u.len();
        info
    }

    /// Replaces `(U, C)` by the `k` harmonic Ritz vectors of the cycle (see the module documentation).
    /// `v`: Arnoldi vectors `V_1..V_{j+1}`; `bmat[j]`: column `j` of `B = Cᴴ A V`; `hbar[j]`: column `j` of the
    /// *unrotated* `H̄` (length `j+2`); `j`: Arnoldi steps done; `kc`: current recycled dimension.
    fn update_recycle(&mut self, v: &[Vec<T>], bmat: &[Vec<T>], hbar: &[Vec<T>], j: usize, kc: usize) {
        let n = v[0].len();
        let nc = kc + j; // columns of Ḡ
        let nr = kc + j + 1; // rows of Ḡ
        if v.len() < j + 1 {
            return;
        }
        let want = self.opts.k.min(nc);
        if want == 0 {
            return;
        }
        // Û = U D
        let dscale: Vec<f64> = self.u.iter().map(|u| 1.0 / norm2(u).max(f64::MIN_POSITIVE)).collect();
        // Ḡ
        let mut g = Dense::<C64>::zeros(nr, nc);
        for i in 0..kc {
            g[(i, i)] = C64::new(dscale[i], 0.0);
        }
        for jj in 0..j {
            for i in 0..kc {
                g[(i, kc + jj)] = bmat[jj][i].to_c64();
            }
            for i in 0..jj + 2 {
                g[(kc + i, kc + jj)] = hbar[jj][i].to_c64();
            }
        }
        // F = V̂ᴴ Ŵ
        let mut f = Dense::<C64>::zeros(nr, nc);
        for b in 0..kc {
            for a in 0..kc {
                f[(a, b)] = dotc(&self.c[a], &self.u[b]).to_c64().scale(dscale[b]);
            }
            for a in 0..=j {
                f[(kc + a, b)] = dotc(&v[a], &self.u[b]).to_c64().scale(dscale[b]);
            }
        }
        for jj in 0..j {
            f[(kc + jj, kc + jj)] = C64::ONE;
        }
        let Some(p) = harmonic_ritz_coefficients::<T>(&g, &f, want) else { return };
        let knew = p.len();
        if knew == 0 {
            return;
        }
        // M = Ḡ P  (nr × knew), thin QR by MGS
        let mut q: Vec<Vec<T>> = Vec::with_capacity(knew); // columns of Q (length nr)
        let mut rmat = vec![vec![T::ZERO; knew]; knew];
        let mut keep = 0usize;
        for pc in p.iter() {
            let mut mcol = vec![T::ZERO; nr];
            for (cidx, &pv) in pc.iter().enumerate() {
                if pv == T::ZERO {
                    continue;
                }
                for (i, mi) in mcol.iter_mut().enumerate() {
                    *mi += T::from_parts(g[(i, cidx)].re, g[(i, cidx)].im) * pv;
                }
            }
            let norm0 = norm2(&mcol);
            for _ in 0..2 {
                for (qi, qv) in q.iter().enumerate() {
                    let h = dotc(qv, &mcol);
                    axpy(-h, qv, &mut mcol);
                    rmat[qi][keep] += h;
                }
            }
            let nrm = norm2(&mcol);
            if nrm <= 1e-10 * norm0.max(f64::MIN_POSITIVE) || !nrm.is_finite() {
                break; // rank deficient: keep what we have
            }
            rmat[keep][keep] = T::from_f64(nrm);
            scale_real(1.0 / nrm, &mut mcol);
            q.push(mcol);
            keep += 1;
        }
        let knew = keep;
        if knew == 0 {
            return;
        }
        // C_new = V̂ Q,   Ŵ P = Û P_top + V P_bot, then U_new = (Ŵ P) R⁻¹
        let mut c_new: Vec<Vec<T>> = Vec::with_capacity(knew);
        let mut u_new: Vec<Vec<T>> = Vec::with_capacity(knew);
        for jcol in 0..knew {
            let mut cc = vec![T::ZERO; n];
            for (ca, &qa) in self.c.iter().zip(&q[jcol]) {
                axpy(qa, ca, &mut cc);
            }
            for (va, &qa) in v.iter().zip(&q[jcol][kc..]) {
                axpy(qa, va, &mut cc);
            }
            c_new.push(cc);
            let mut wp = vec![T::ZERO; n];
            for a in 0..kc {
                axpy(p[jcol][a].scale(dscale[a]), &self.u[a], &mut wp);
            }
            for a in 0..j {
                axpy(p[jcol][kc + a], &v[a], &mut wp);
            }
            for (i, ui) in u_new.iter().enumerate() {
                axpy(-rmat[i][jcol], ui, &mut wp);
            }
            scale_real(1.0 / rmat[jcol][jcol].re(), &mut wp);
            u_new.push(wp);
        }
        self.u = u_new;
        self.c = c_new;
        self.stale = false;
    }
}

/// The `want` harmonic Ritz vectors (coefficient vectors `P`, each of length `Ḡ.cols()`, already real for real
/// `T`) with smallest harmonic Ritz values `θ`: `ḠᴴḠ g = θ Ḡᴴ F g`. Returns `None` if the dense problem is
/// numerically singular or the Schur iteration fails.
fn harmonic_ritz_coefficients<T: Scalar>(g: &Dense<C64>, f: &Dense<C64>, want: usize) -> Option<Vec<Vec<T>>> {
    let nc = g.cols();
    let gh = g.adjoint();
    let a1 = gh.matmul(g);
    let b1 = gh.matmul(f);
    // Cholesky a1 = L Lᴴ
    let mut l = Dense::<C64>::zeros(nc, nc);
    let scale = a1.frobenius().max(f64::MIN_POSITIVE);
    for j in 0..nc {
        let mut d = a1[(j, j)].re;
        for k in 0..j {
            d -= l[(j, k)].norm_sqr();
        }
        if d <= 1e-14 * scale || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        l[(j, j)] = C64::new(d, 0.0);
        for i in j + 1..nc {
            let mut s = a1[(i, j)];
            for k in 0..j {
                s -= l[(i, k)] * l[(j, k)].conj();
            }
            l[(i, j)] = s / d;
        }
    }
    // Cm = L⁻¹ B1 L⁻ᴴ
    let mut x = b1.clone(); // solve L X = B1 (forward substitution, column by column)
    for col in 0..nc {
        for i in 0..nc {
            let mut s = x[(i, col)];
            for k in 0..i {
                s -= l[(i, k)] * x[(k, col)];
            }
            x[(i, col)] = s / l[(i, i)];
        }
    }
    // Cm = X L⁻ᴴ  ⇔  Cm Lᴴ = X, solve row by row (Lᴴ upper triangular)
    let mut cm = x;
    for row in 0..nc {
        for j in 0..nc {
            let mut s = cm[(row, j)];
            for k in 0..j {
                s -= cm[(row, k)] * l[(j, k)].conj();
            }
            cm[(row, j)] = s / l[(j, j)].conj();
        }
    }
    let (t, z) = schur(&cm).ok()?;
    // eigenpairs μ_i = T_ii, y_i = Z · (triangular eigenvector); θ_i = 1/μ_i
    let mut order: Vec<usize> = (0..nc).collect();
    order.sort_by(|&i, &j| t[(j, j)].norm().total_cmp(&t[(i, i)].norm())); // largest |μ| (smallest |θ|) first
    let mut vecs: Vec<(C64, Vec<C64>)> = Vec::new();
    for &i in order.iter().take(want.max(2) + 2) {
        let yt = triangular_eigenvector(&t, nc, i);
        let y: Vec<C64> = (0..nc).map(|r| (0..nc).fold(C64::ZERO, |s, k| s + z[(r, k)] * yt[k])).collect();
        // g = L⁻ᴴ y (backward substitution with Lᴴ)
        let mut gv = y;
        for r in (0..nc).rev() {
            let mut s = gv[r];
            for k in r + 1..nc {
                s -= l[(k, r)].conj() * gv[k];
            }
            gv[r] = s / l[(r, r)];
        }
        vecs.push((t[(i, i)], gv));
    }
    let mut out: Vec<Vec<T>> = Vec::new();
    if T::IS_COMPLEX {
        for (_, gv) in vecs.iter().take(want) {
            out.push(gv.iter().map(|c| T::from_parts(c.re, c.im)).collect());
        }
        return Some(out);
    }
    // real matrices: a complex vector and its conjugate partner give the same real 2-D subspace
    let mut used: Vec<C64> = Vec::new();
    for (mu, gv) in &vecs {
        if out.len() >= want {
            break;
        }
        if used.iter().any(|u| (u.conj() - *mu).norm() <= 1e-8 * mu.norm().max(1e-300)) {
            continue; // conjugate partner of a pair already taken
        }
        // rotate the phase so that the largest component is real
        let big = gv.iter().max_by(|a, b| a.norm().total_cmp(&b.norm()))?;
        let ph = big.conj() / big.norm();
        let gr: Vec<C64> = gv.iter().map(|c| *c * ph).collect();
        let imag_part = gr.iter().map(|c| c.im.abs()).fold(0.0, f64::max);
        let real_part = gr.iter().map(|c| c.re.abs()).fold(0.0, f64::max);
        out.push(gr.iter().map(|c| T::from_f64(c.re)).collect());
        if imag_part > 1e-8 * real_part && out.len() < want {
            out.push(gr.iter().map(|c| T::from_f64(c.im)).collect());
        }
        used.push(*mu);
    }
    Some(out)
}
