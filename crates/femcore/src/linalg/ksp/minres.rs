//! MINRES (Paige & Saunders 1975) for real symmetric or complex Hermitian *indefinite* systems with an
//! SPD preconditioner, in the formulation of Choi & Saunders. The monitored norm is the
//! preconditioned residual norm `√(rᴴM⁻¹r)` obtained from the recurrence (PETSc's `KSPMINRES` monitors the same
//! quantity), so the iteration counts can be compared with PETSc's.

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{axpy, dotc, norm2};

#[derive(Debug, Clone, Default)]
pub struct MinresOptions {
    pub tol: Tolerances,
    pub record_history: bool,
}

pub struct Minres<T: Scalar> {
    pub opts: MinresOptions,
    _p: std::marker::PhantomData<T>,
}

impl<T: Scalar> Minres<T> {
    pub fn new(opts: MinresOptions) -> Self {
        Self { opts, _p: Default::default() }
    }

    /// Solves `A x = b` from `x = 0`. `A` Hermitian (possibly indefinite), `M` Hermitian positive definite.
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
        let mut r1 = b.to_vec();
        let mut y = vec![T::ZERO; n];
        m.apply(&r1, &mut y);
        info.pc_applies += 1;
        let b1 = dotc(&r1, &y).re();
        if b1 < 0.0 {
            info.reason = ConvergedReason::DivergedIndefinitePc;
            return info;
        }
        let beta1 = b1.sqrt();
        conv.set_reference(beta1);
        info.residual_norm = beta1;
        if self.opts.record_history {
            info.history.push(beta1);
        }
        if let Some(r) = conv.check(beta1) {
            info.reason = r;
            return info;
        }
        let mut r2 = r1.clone();
        let (mut oldb, mut beta, mut dbar, mut epsln, mut phibar) = (0.0f64, beta1, 0.0f64, 0.0f64, beta1);
        let (mut cs, mut sn) = (-1.0f64, 0.0f64);
        let (mut w, mut w1, mut w2) = (vec![T::ZERO; n], vec![T::ZERO; n], vec![T::ZERO; n]);
        let mut v = vec![T::ZERO; n];
        for k in 1..=self.opts.tol.max_it {
            info.iterations = k;
            let s = 1.0 / beta;
            for (vi, yi) in v.iter_mut().zip(&y) {
                *vi = yi.scale(s);
            }
            a.apply(&v, &mut y);
            info.matvecs += 1;
            if k >= 2 {
                axpy(T::from_f64(-beta / oldb), &r1, &mut y);
            }
            let alfa = dotc(&v, &y).re();
            axpy(T::from_f64(-alfa / beta), &r2, &mut y);
            std::mem::swap(&mut r1, &mut r2); // r1 ← r2
            r2.copy_from_slice(&y);
            m.apply(&r2, &mut y);
            info.pc_applies += 1;
            oldb = beta;
            let b2 = dotc(&r2, &y).re();
            if b2 < 0.0 {
                info.reason = ConvergedReason::DivergedIndefinitePc;
                break;
            }
            beta = b2.sqrt();
            // previous rotation, then the new one
            let oldeps = epsln;
            let delta = cs * dbar + sn * alfa;
            let gbar = sn * dbar - cs * alfa;
            epsln = sn * beta;
            dbar = -cs * beta;
            let gamma = (gbar * gbar + beta * beta).sqrt();
            if gamma == 0.0 {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            cs = gbar / gamma;
            sn = beta / gamma;
            let phi = cs * phibar;
            phibar *= sn;
            let denom = 1.0 / gamma;
            std::mem::swap(&mut w1, &mut w2); // w1 ← old w2
            std::mem::swap(&mut w2, &mut w); // w2 ← old w
            for i in 0..n {
                w[i] = (v[i] - w1[i].scale(oldeps) - w2[i].scale(delta)).scale(denom);
            }
            axpy(T::from_f64(phi), &w, x);
            let rnorm = phibar.abs();
            info.residual_norm = rnorm;
            if self.opts.record_history {
                info.history.push(rnorm);
            }
            if let Some(reason) = conv.check(rnorm) {
                info.reason = reason;
                break;
            }
            if beta == 0.0 {
                info.reason = ConvergedReason::Atol;
                break;
            }
        }
        let _ = norm2::<T>;
        info
    }
}
