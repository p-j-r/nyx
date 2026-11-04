use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances, ensure_len};
use crate::linalg::basis::Basis;
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops;

/// Which side the preconditioner is applied on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PcSide {
    /// Solve `M⁻¹A x = M⁻¹b`; the monitored norm is `‖M⁻¹r‖` (PETSc default).
    #[default]
    Left,
    /// Solve `A M⁻¹ y = b`, `x = M⁻¹y`; the monitored norm is the true `‖r‖`.
    Right,
}

/// Gram–Schmidt variant used to orthogonalise each new Krylov vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orthogonalization {
    /// Classical Gram–Schmidt, one pass (PETSc default; one fused multi-dot per step).
    #[default]
    ClassicalGramSchmidt,
    /// Classical GS with a second pass when the norm drops by more than `1/√2`.
    ClassicalGramSchmidtIfNeeded,
    ModifiedGramSchmidt,
}

#[derive(Debug, Clone)]
pub struct GmresOptions {
    pub tol: Tolerances,
    pub restart: usize,
    pub side: PcSide,
    pub orthog: Orthogonalization,
    /// Flexible GMRES (FGMRES): right preconditioning with the preconditioned vectors stored, so
    /// the preconditioner may change at every application (inner iterations, adaptive AMG).
    /// Forces `side = Right`; costs one extra basis of memory.
    pub flexible: bool,
    pub initial_guess_nonzero: bool,
    pub record_history: bool,
}

impl Default for GmresOptions {
    fn default() -> Self {
        Self {
            tol: Tolerances::default(),
            restart: 30,
            side: PcSide::Left,
            orthog: Orthogonalization::ClassicalGramSchmidt,
            flexible: false,
            initial_guess_nonzero: false,
            record_history: false,
        }
    }
}

/// Restarted GMRES(m) / FGMRES with complex Givens-rotation least-squares updates, for real or
/// complex (including complex symmetric) systems.
#[derive(Debug, Clone)]
pub struct Gmres<T: Scalar> {
    pub opts: GmresOptions,
    v: Option<Basis<T>>,
    zb: Option<Basis<T>>,
    r: Vec<T>,
    t: Vec<T>,
    u: Vec<T>,
}

const HAPPY_TOL: f64 = 1e-30;
const BREAKDOWN_TOL: f64 = 0.1;

/// Plane rotation `[c s; −s̄ c]` with real `c` that maps `(a, b)` to `(r, 0)`.
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

impl<T: Scalar> Gmres<T> {
    pub fn new(opts: GmresOptions) -> Self {
        Self { opts, v: None, zb: None, r: Vec::new(), t: Vec::new(), u: Vec::new() }
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
        let restart = self.opts.restart.max(1);
        let opts = self.opts.clone();
        let flexible = opts.flexible;
        let left = opts.side == PcSide::Left && !flexible;
        if self.v.as_ref().is_none_or(|v| v.len() != n || v.capacity() != restart + 1) {
            self.v = Some(Basis::new(n, restart + 1));
        }
        if flexible && self.zb.as_ref().is_none_or(|v| v.len() != n || v.capacity() != restart) {
            self.zb = Some(Basis::new(n, restart));
        }
        for v in [&mut self.r, &mut self.t, &mut self.u] {
            ensure_len(v, n);
        }
        let ld = restart + 1;
        let mut h = vec![T::ZERO; ld * restart];
        let mut cs = vec![0.0f64; restart];
        let mut sn = vec![T::ZERO; restart];
        let mut rs = vec![T::ZERO; ld];
        let Self { v, zb, r, t, u, .. } = self;
        let v = v.as_mut().unwrap();

        let mut conv = ConvTest::new(opts.tol);
        let mut info = SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedIts,
            residual_norm: 0.0,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };
        let mut guess_nonzero = opts.initial_guess_nonzero;
        if !guess_nonzero {
            x.fill(T::ZERO);
        }
        let mut first_cycle = true;
        let mut last_rnorm = -1.0f64;
        let mut its = 0usize;

        'outer: loop {
            // Initial residual of this cycle in v[0]: left: M⁻¹(b − A x), else b − A x.
            if guess_nonzero {
                a.apply(x, r);
                info.matvecs += 1;
                vecops::aypx(-T::ONE, b, r);
            } else {
                r.copy_from_slice(b);
            }
            if left {
                m.apply(r, v.col_mut(0));
                info.pc_applies += 1;
            } else {
                v.col_mut(0).copy_from_slice(r);
            }
            let res = vecops::norm2(v.col(0));
            if !res.is_finite() {
                info.reason = ConvergedReason::DivergedNanOrInf;
                info.residual_norm = res;
                break;
            }
            if first_cycle {
                let reference = if guess_nonzero {
                    let s = if left {
                        m.apply(b, t);
                        info.pc_applies += 1;
                        vecops::norm2(t)
                    } else {
                        vecops::norm2(b)
                    };
                    if s == 0.0 { res } else { s }
                } else {
                    res
                };
                conv.set_reference(reference);
                if opts.record_history {
                    info.history.push(res);
                }
            } else if last_rnorm > 0.0
                && (res - last_rnorm).abs() > BREAKDOWN_TOL * conv.rnorm0().max(f64::MIN_POSITIVE)
            {
                info.reason = ConvergedReason::DivergedBreakdown;
                info.residual_norm = res;
                break;
            }
            info.residual_norm = res;
            if res == 0.0 {
                info.reason = ConvergedReason::Atol;
                break;
            }
            if let Some(reason) = conv.check(res) {
                info.reason = reason;
                break;
            }
            vecops::scale_real(1.0 / res, v.col_mut(0));
            rs.fill(T::ZERO);
            rs[0] = T::from_f64(res);

            let mut k = 0usize;
            let mut reason = None;
            let mut res_k = res;
            while k < restart && its < opts.tol.max_it {
                {
                    let (vk, w) = v.col_pair_mut(k, k + 1);
                    if left {
                        a.apply(vk, t);
                        m.apply(t, w);
                    } else if flexible {
                        let z = zb.as_mut().unwrap();
                        m.apply(vk, z.col_mut(k));
                        a.apply(z.col(k), w);
                    } else {
                        m.apply(vk, t);
                        a.apply(t, w);
                    }
                }
                info.matvecs += 1;
                info.pc_applies += 1;
                let hcol = &mut h[k * ld..(k + 1) * ld];
                let tt = orthogonalize(v, k, opts.orthog, hcol);
                if !tt.is_finite() {
                    reason = Some(ConvergedReason::DivergedNanOrInf);
                    break;
                }
                if tt != 0.0 {
                    vecops::scale_real(1.0 / tt, v.col_mut(k + 1));
                }
                let hapbnd = (tt / rs[k].modulus()).min(HAPPY_TOL);
                let happy = tt < hapbnd;
                for j in 0..k {
                    let (c, s) = (cs[j], sn[j]);
                    let (hj, hj1) = (hcol[j], hcol[j + 1]);
                    hcol[j] = hj.scale(c) + s * hj1;
                    hcol[j + 1] = hj1.scale(c) - s.conj() * hj;
                }
                hcol[k + 1] = T::from_f64(tt);
                if !happy {
                    let (c, s) = givens(hcol[k], hcol[k + 1]);
                    if c == 1.0 && hcol[k] == T::ZERO && hcol[k + 1] == T::ZERO {
                        reason = Some(ConvergedReason::DivergedBreakdown);
                        break;
                    }
                    cs[k] = c;
                    sn[k] = s;
                    let rk = rs[k];
                    rs[k + 1] = -(s.conj() * rk);
                    rs[k] = rk.scale(c);
                    hcol[k] = hcol[k].scale(c) + s * hcol[k + 1];
                    res_k = rs[k + 1].modulus();
                } else {
                    res_k = 0.0;
                }
                hcol[k + 1] = T::ZERO;
                k += 1;
                its += 1;
                info.iterations = its;
                info.residual_norm = res_k;
                if opts.record_history {
                    info.history.push(res_k);
                }
                reason = conv.check(res_k);
                if happy && reason.is_none() {
                    reason = Some(ConvergedReason::DivergedBreakdown);
                }
                if reason.is_some() {
                    break;
                }
            }

            // Solve the k×k upper-triangular system H y = rs and update x.
            if k > 0 {
                let mut y = vec![T::ZERO; k];
                for i in (0..k).rev() {
                    let mut s = rs[i];
                    for j in i + 1..k {
                        s -= h[j * ld + i] * y[j];
                    }
                    let d = h[i * ld + i];
                    if d == T::ZERO {
                        reason = Some(ConvergedReason::DivergedBreakdown);
                    } else {
                        y[i] = s / d;
                    }
                }
                if left {
                    v.view().add_combination(k, &y, x);
                } else if flexible {
                    zb.as_ref().unwrap().view().add_combination(k, &y, x);
                } else {
                    t.fill(T::ZERO);
                    v.view().add_combination(k, &y, t);
                    m.apply(t, u);
                    info.pc_applies += 1;
                    vecops::axpy(T::ONE, u, x);
                }
            }
            last_rnorm = res_k;
            if let Some(reason) = reason {
                info.reason = reason;
                break 'outer;
            }
            if its >= opts.tol.max_it {
                info.reason = ConvergedReason::DivergedIts;
                break;
            }
            guess_nonzero = true;
            first_cycle = false;
        }
        info
    }
}

/// Orthogonalise `v[k+1]` against `v[0..=k]`, filling `h[0..=k]`; returns the remaining norm.
fn orthogonalize<T: Scalar>(v: &mut Basis<T>, k: usize, kind: Orthogonalization, h: &mut [T]) -> f64 {
    let (head, tail) = v.split_cols_mut(k + 1);
    let w = &mut tail[..head.len()];
    match kind {
        Orthogonalization::ClassicalGramSchmidt => {
            head.dots(k + 1, w, h);
            head.sub_combination(k + 1, h, w);
            vecops::norm2(w)
        }
        Orthogonalization::ClassicalGramSchmidtIfNeeded => {
            let before = vecops::norm2(w);
            head.dots(k + 1, w, h);
            head.sub_combination(k + 1, h, w);
            let mut after = vecops::norm2(w);
            if after < std::f64::consts::FRAC_1_SQRT_2 * before {
                let mut h2 = vec![T::ZERO; k + 1];
                head.dots(k + 1, w, &mut h2);
                head.sub_combination(k + 1, &h2, w);
                for (a, b) in h.iter_mut().zip(&h2) {
                    *a += *b;
                }
                after = vecops::norm2(w);
            }
            after
        }
        Orthogonalization::ModifiedGramSchmidt => {
            for (j, hj) in h.iter_mut().enumerate().take(k + 1) {
                let c = vecops::dotc(head.col(j), w);
                *hj = c;
                vecops::axpy(-c, head.col(j), w);
            }
            vecops::norm2(w)
        }
    }
}
