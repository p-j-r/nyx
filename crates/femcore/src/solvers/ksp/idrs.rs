//! IDR(s) (Sonneveld & van Gijzen 2008; the biorthogonal variant of van Gijzen & Sonneveld, ACM TOMS 38, 2011):
//! a short-recurrence method for general (also complex) systems that needs only `s + 1` matrix products per `s + 1`
//! dimension reductions and `O(s)` vectors, typically between BiCGStab and GMRES in iteration count and memory.
//! Right-preconditioned; the monitored norm is the true residual. `s = 1` is close to BiCGStab.

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::solvers::linop::{LinearOperator, Preconditioner};
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::Scalar;
use crate::solvers::vecops::{axpy, dotc, norm2};

#[derive(Debug, Clone)]
pub struct IdrsOptions {
    pub tol: Tolerances,
    /// Dimension of the shadow space.
    pub s: usize,
    /// Minimum cosine of the angle between `t` and `r` in the `ω` choice ("maintaining the convergence", 0.7).
    pub angle: f64,
    pub seed: u64,
    pub record_history: bool,
}

impl Default for IdrsOptions {
    fn default() -> Self {
        Self { tol: Tolerances::default(), s: 4, angle: 0.7, seed: 0x1D25, record_history: false }
    }
}

pub struct Idrs<T: Scalar> {
    pub opts: IdrsOptions,
    _p: std::marker::PhantomData<T>,
}

impl<T: Scalar> Idrs<T> {
    pub fn new(opts: IdrsOptions) -> Self {
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
        let s = self.opts.s.max(1).min(n);
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
        // shadow space P: s random vectors, orthonormalised
        let mut p: Vec<Vec<T>> = Vec::with_capacity(s);
        for j in 0..s {
            let mut v = random_vector::<T>(n, self.opts.seed.wrapping_add(31 * j as u64 + 1));
            for q in &p {
                let h = dotc(q, &v);
                axpy(-h, q, &mut v);
            }
            let nv = norm2(&v);
            v.iter_mut().for_each(|e| *e = e.scale(1.0 / nv));
            p.push(v);
        }
        let mut g: Vec<Vec<T>> = vec![vec![T::ZERO; n]; s];
        let mut u: Vec<Vec<T>> = vec![vec![T::ZERO; n]; s];
        let mut mm = vec![vec![T::ZERO; s]; s];
        for (i, row) in mm.iter_mut().enumerate() {
            row[i] = T::ONE;
        }
        let mut omega = T::ONE;
        let mut f: Vec<T> = p.iter().map(|pi| dotc(pi, &r)).collect();
        let (mut v, mut t) = (vec![T::ZERO; n], vec![T::ZERO; n]);
        let mut iter = 0usize;
        'outer: while iter < self.opts.tol.max_it {
            for k in 0..s {
                // solve the lower-triangular system M[k..s, k..s] c = f[k..s]
                let mut c = vec![T::ZERO; s];
                for i in k..s {
                    let mut sum = f[i];
                    for j in k..i {
                        sum -= mm[i][j] * c[j];
                    }
                    if mm[i][i] == T::ZERO {
                        info.reason = ConvergedReason::DivergedBreakdown;
                        break 'outer;
                    }
                    c[i] = sum / mm[i][i];
                }
                v.copy_from_slice(&r);
                for j in k..s {
                    axpy(-c[j], &g[j], &mut v);
                }
                let mut vp = vec![T::ZERO; n];
                m.apply(&v, &mut vp);
                info.pc_applies += 1;
                // u_k = Σ c_j u_j + ω v
                let mut uk = vec![T::ZERO; n];
                for j in k..s {
                    axpy(c[j], &u[j], &mut uk);
                }
                axpy(omega, &vp, &mut uk);
                let mut gk = vec![T::ZERO; n];
                a.apply(&uk, &mut gk);
                info.matvecs += 1;
                // make g_k orthogonal to p_0 .. p_{k-1}
                for i in 0..k {
                    let alpha = dotc(&p[i], &gk) / mm[i][i];
                    axpy(-alpha, &g[i], &mut gk);
                    axpy(-alpha, &u[i], &mut uk);
                }
                for i in k..s {
                    mm[i][k] = dotc(&p[i], &gk);
                }
                if mm[k][k] == T::ZERO {
                    info.reason = ConvergedReason::DivergedBreakdown;
                    break 'outer;
                }
                let beta = f[k] / mm[k][k];
                axpy(-beta, &gk, &mut r);
                axpy(beta, &uk, x);
                g[k] = gk;
                u[k] = uk;
                iter += 1;
                info.iterations = iter;
                let rn = norm2(&r);
                info.residual_norm = rn;
                if self.opts.record_history {
                    info.history.push(rn);
                }
                if let Some(reason) = conv.check(rn) {
                    info.reason = reason;
                    break 'outer;
                }
                for i in k + 1..s {
                    f[i] -= beta * mm[i][k];
                }
            }
            // entering the next G space
            m.apply(&r, &mut v);
            info.pc_applies += 1;
            a.apply(&v, &mut t);
            info.matvecs += 1;
            let (tt, tr) = (norm2(&t), dotc(&t, &r));
            let rn = norm2(&r);
            let rho = if tt * rn > 0.0 { tr.modulus() / (tt * rn) } else { 0.0 };
            let mut om = if tt > 0.0 { tr / T::from_f64(tt * tt) } else { T::ZERO };
            if rho < self.opts.angle && rho > 0.0 {
                om = om.scale(self.opts.angle / rho);
            }
            if om == T::ZERO {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            omega = om;
            axpy(-omega, &t, &mut r);
            axpy(omega, &v, x);
            let rn = norm2(&r);
            info.residual_norm = rn;
            iter += 1;
            info.iterations = iter;
            if self.opts.record_history {
                info.history.push(rn);
            }
            if let Some(reason) = conv.check(rn) {
                info.reason = reason;
                break;
            }
            for i in 0..s {
                f[i] = dotc(&p[i], &r);
            }
        }
        info
    }
}
