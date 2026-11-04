//! BiCGStab (van der Vorst 1992), right-preconditioned: it monitors the *true* residual `‖b − Ax‖`
//! and the preconditioner may be any (also nonsymmetric) approximation of `A⁻¹`. Real and complex.
//! For complex symmetric systems the conjugated shadow-residual products are the standard choice
//! (the same algorithm as the real one for real scalars).

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{axpy, dotc, norm2};

#[derive(Debug, Clone, Default)]
pub struct BicgstabOptions {
    pub tol: Tolerances,
    pub record_history: bool,
}

pub struct Bicgstab<T: Scalar> {
    pub opts: BicgstabOptions,
    _p: std::marker::PhantomData<T>,
}

impl<T: Scalar> Bicgstab<T> {
    pub fn new(opts: BicgstabOptions) -> Self {
        Self { opts, _p: Default::default() }
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
        x.fill(T::ZERO);
        let mut info = SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedIts,
            residual_norm: 0.0,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };
        let mut conv = ConvTest::new(self.opts.tol);
        let mut r = b.to_vec();
        let rhat = r.clone();
        let rn0 = norm2(&r);
        conv.set_reference(rn0);
        info.residual_norm = rn0;
        if self.opts.record_history {
            info.history.push(rn0);
        }
        if let Some(reason) = conv.check(rn0) {
            info.reason = reason;
            return info;
        }
        let (mut rho, mut alpha, mut omega) = (T::ONE, T::ONE, T::ONE);
        let (mut v, mut p) = (vec![T::ZERO; n], vec![T::ZERO; n]);
        let (mut phat, mut s, mut shat, mut t) =
            (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
        for k in 1..=self.opts.tol.max_it {
            info.iterations = k;
            let rho_new = dotc(&rhat, &r);
            if rho_new == T::ZERO || !rho_new.is_finite() {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            let beta = (rho_new / rho) * (alpha / omega);
            for i in 0..n {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
            m.apply(&p, &mut phat);
            info.pc_applies += 1;
            a.apply(&phat, &mut v);
            info.matvecs += 1;
            let den = dotc(&rhat, &v);
            if den == T::ZERO || !den.is_finite() {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            alpha = rho_new / den;
            for i in 0..n {
                s[i] = r[i] - alpha * v[i];
            }
            m.apply(&s, &mut shat);
            info.pc_applies += 1;
            a.apply(&shat, &mut t);
            info.matvecs += 1;
            let tt = dotc(&t, &t);
            omega = if tt == T::ZERO { T::ZERO } else { dotc(&t, &s) / tt };
            axpy(alpha, &phat, x);
            axpy(omega, &shat, x);
            for i in 0..n {
                r[i] = s[i] - omega * t[i];
            }
            rho = rho_new;
            let rn = norm2(&r);
            info.residual_norm = rn;
            if self.opts.record_history {
                info.history.push(rn);
            }
            if let Some(reason) = conv.check(rn) {
                info.reason = reason;
                break;
            }
            if omega == T::ZERO {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
        }
        info
    }
}
