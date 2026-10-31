//! COCR (Sogabe & Zhang 2007): the conjugate *residual* method with the bilinear form `xᵀy`, for complex
//! symmetric systems. The literature reports smoother residual histories than COCG's on indefinite complex
//! symmetric problems (not verified here: the tests only check convergence and the solution). Left preconditioning with a complex symmetric `M⁻¹`:
//! it is conjugate residuals for `M⁻¹A`, which is symmetric in the `M`-bilinear form, so only `M⁻¹`
//! is needed:
//!
//! ```text
//! z = M⁻¹r, p = z, q = A p, s = M⁻¹q
//! α = (zᵀAz) / (qᵀs),  x += α p,  z −= α s
//! β = (z'ᵀAz') / (zᵀAz),  p = z' + β p,  q = A z' + β q,  s = M⁻¹ q
//! ```
//!
//! One product with `A` and one application of `M⁻¹` per iteration. The monitored norm is `‖z‖₂`
//! (preconditioned residual).

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::solvers::linop::{LinearOperator, Preconditioner};
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::{axpy, dotu, norm2};

#[derive(Debug, Clone, Default)]
pub struct CocrOptions {
    pub tol: Tolerances,
    pub record_history: bool,
}

pub struct Cocr<T: Scalar> {
    pub opts: CocrOptions,
    _p: std::marker::PhantomData<T>,
}

impl<T: Scalar> Cocr<T> {
    pub fn new(opts: CocrOptions) -> Self {
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
        let mut z = vec![T::ZERO; n];
        m.apply(b, &mut z);
        info.pc_applies += 1;
        let zn = norm2(&z);
        conv.set_reference(zn);
        info.residual_norm = zn;
        if self.opts.record_history {
            info.history.push(zn);
        }
        if let Some(r) = conv.check(zn) {
            info.reason = r;
            return info;
        }
        let mut p = z.clone();
        let mut q = vec![T::ZERO; n];
        a.apply(&p, &mut q);
        info.matvecs += 1;
        let mut s = vec![T::ZERO; n];
        m.apply(&q, &mut s);
        info.pc_applies += 1;
        let mut az = q.clone(); // A z
        let mut zaz = dotu(&z, &az);
        for k in 1..=self.opts.tol.max_it {
            info.iterations = k;
            let qs = dotu(&q, &s);
            if qs == T::ZERO || !qs.is_finite() || zaz == T::ZERO {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            let alpha = zaz / qs;
            axpy(alpha, &p, x);
            axpy(-alpha, &s, &mut z);
            let zn = norm2(&z);
            info.residual_norm = zn;
            if self.opts.record_history {
                info.history.push(zn);
            }
            if let Some(reason) = conv.check(zn) {
                info.reason = reason;
                break;
            }
            a.apply(&z, &mut az);
            info.matvecs += 1;
            let zaz_new = dotu(&z, &az);
            let beta = zaz_new / zaz;
            zaz = zaz_new;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
                q[i] = az[i] + beta * q[i];
            }
            m.apply(&q, &mut s);
            info.pc_applies += 1;
        }
        info
    }
}
