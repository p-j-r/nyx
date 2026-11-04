//! QMR (Freund & Nachtigal 1991): the quasi-minimal residual method built on the two-sided
//! (unsymmetric) Lanczos process, for general real or complex systems, plus the
//! **complex-symmetric** variant (Freund 1992) that needs only *one* product with `A` per iteration.
//!
//! Why QMR: BiCG-type methods (BiCGStab is the other member here) can see their residual jump
//! up and down by orders of magnitude; QMR minimises the *quasi*-residual, which removes the
//! peaks and gives a smooth, nearly monotone history. Short recurrences (no stored basis, unlike
//! GMRES), at the price of no optimality guarantee for the true residual and possible breakdown of
//! the Lanczos process (no look-ahead here: breakdown is reported as
//! [`ConvergedReason::DivergedBreakdown`]).
//!
//! The recurrences are those of the *Templates* book (Barrett et al. 1994, Algorithm "Quasi-Minimal
//! Residual") with a left preconditioner `M⁻¹` only (`M₂ = I`), generalised to complex scalars by
//! using the unconjugated bilinear form `xᵀy` for the two Lanczos sequences and moduli where the real
//! algorithm uses absolute values:
//!
//! ```text
//! ṽ₁ = r₀,  y = M⁻¹ṽ₁,  ρ₁ = ‖y‖;   w̃₁ = r₀,  ξ₁ = ‖w̃₁‖          (r₀ = b, x₀ = 0)
//! for i = 1, 2, …
//!     v = ṽ/ρ,  y /= ρ,  w = w̃/ξ,  z = w̃/ξ
//!     δ = zᵀy,   z̃ = M⁻ᵀ z
//!     p = y − (ξδ/ε) p,   q = z̃ − (ρδ/ε) q                       (i > 1)
//!     p̃ = A p,   ε = qᵀp̃,   β = ε/δ
//!     ṽ = p̃ − βv,  y = M⁻¹ṽ,  ρ⁺ = ‖y‖;   w̃ = Aᵀq − βw,  ξ⁺ = ‖w̃‖
//!     θ = ρ⁺/(γ|β|),  γ⁺ = 1/√(1+θ²),  η⁺ = −η ρ γ⁺² / (β γ²)
//!     d = η⁺p + (θγ⁺)² d,   s = η⁺p̃ + (θγ⁺)² s,   x += d,   r −= s
//! ```
//!
//! Two matrix products per iteration (`A` and `Aᵀ`, the *transpose*, not the adjoint). If
//! `A = Aᵀ` and `M = Mᵀ` (complex symmetric: PML, Helmholtz with absorption) the second Lanczos
//! sequence is a multiple of the first: choosing `ξ = ρ` makes `w̃ = ṽ` and `q = p`, so
//! [`Qmr::solve_symmetric`] skips the transposed product and the transposed preconditioner solve:
//! one `A` product and one `M⁻¹` per iteration, the same cost as COCG.
//!
//! The monitored norm is the true residual `‖b − Ax‖₂` (updated recursively, so after many
//! iterations it may drift from the true one; the tests compare both).

use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances};
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{axpy, dotu, norm2, scale_real};

#[derive(Debug, Clone, Default)]
pub struct QmrOptions {
    pub tol: Tolerances,
    pub record_history: bool,
}

pub struct Qmr<T: Scalar> {
    pub opts: QmrOptions,
    _p: std::marker::PhantomData<T>,
}

/// The second Lanczos sequence: a transposed operator and preconditioner, or the shortcut for
/// symmetric `A`, `M`.
enum Left<'a, T: Scalar> {
    Symmetric,
    General { at: &'a dyn LinearOperator<T>, mt: &'a dyn Preconditioner<T> },
}

impl<T: Scalar> Qmr<T> {
    pub fn new(opts: QmrOptions) -> Self {
        Self { opts, _p: Default::default() }
    }

    /// General systems. `at` must apply `Aᵀ` and `mt` must apply `M⁻ᵀ` (transposes, not adjoints; for a
    /// real matrix they coincide). `mt` has to be the exact transpose of `m`: the left Lanczos sequence
    /// is a Krylov sequence of `(M⁻¹A)ᵀ`, and any other operator breaks the biorthogonality.
    pub fn solve<A, At, M, Mt>(&mut self, a: &A, at: &At, m: &M, mt: &Mt, b: &[T], x: &mut [T]) -> SolveInfo
    where
        A: LinearOperator<T> + ?Sized,
        At: LinearOperator<T>,
        M: Preconditioner<T> + ?Sized,
        Mt: Preconditioner<T>,
    {
        self.run(a, m, Left::General { at, mt }, b, x)
    }

    /// Complex-symmetric (or real symmetric) `A = Aᵀ` with a symmetric preconditioner `M⁻¹ = M⁻ᵀ`:
    /// one product with `A` and one preconditioner solve per iteration.
    pub fn solve_symmetric<A, M>(&mut self, a: &A, m: &M, b: &[T], x: &mut [T]) -> SolveInfo
    where
        A: LinearOperator<T> + ?Sized,
        M: Preconditioner<T> + ?Sized,
    {
        self.run(a, m, Left::Symmetric, b, x)
    }

    fn run<A, M>(&mut self, a: &A, m: &M, left: Left<'_, T>, b: &[T], x: &mut [T]) -> SolveInfo
    where
        A: LinearOperator<T> + ?Sized,
        M: Preconditioner<T> + ?Sized,
    {
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
        let symmetric = matches!(left, Left::Symmetric);

        let zero = vec![T::ZERO; n];
        let mut vt = r.clone(); // ṽ
        let mut wt = r.clone(); // w̃
        let mut y = zero.clone();
        m.apply(&vt, &mut y);
        info.pc_applies += 1;
        let mut rho = norm2(&y);
        let mut xi = if symmetric { rho } else { norm2(&wt) };
        let (mut p, mut q, mut pt) = (zero.clone(), zero.clone(), zero.clone());
        let (mut d, mut s) = (zero.clone(), zero.clone());
        let (mut v, mut w, mut zt) = (zero.clone(), zero.clone(), zero.clone());
        let (mut gamma_prev, mut eta_prev, mut theta_prev) = (1.0f64, -T::ONE, 0.0f64);
        let mut eps_prev = T::ONE;

        for i in 1..=self.opts.tol.max_it {
            if rho == 0.0 || xi == 0.0 || !rho.is_finite() || !xi.is_finite() {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            info.iterations = i;
            // normalise both Lanczos vectors (and the preconditioned right one)
            v.copy_from_slice(&vt);
            scale_real(1.0 / rho, &mut v);
            scale_real(1.0 / rho, &mut y);
            let delta;
            match &left {
                Left::Symmetric => {
                    // z = w̃/ξ = v, z̃ = M⁻¹ z = y
                    delta = dotu(&v, &y);
                    zt.copy_from_slice(&y);
                }
                Left::General { mt, .. } => {
                    w.copy_from_slice(&wt);
                    scale_real(1.0 / xi, &mut w);
                    delta = dotu(&w, &y);
                    mt.apply(&w, &mut zt);
                    info.pc_applies += 1;
                }
            }
            if delta == T::ZERO || !delta.is_finite() {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            if i == 1 {
                p.copy_from_slice(&y);
                q.copy_from_slice(&zt);
            } else {
                // p = y − (ξδ/ε) p,  q = z̃ − (ρδ/ε) q
                let cp = delta.scale(xi) / eps_prev;
                let cq = delta.scale(rho) / eps_prev;
                for k in 0..n {
                    p[k] = y[k] - cp * p[k];
                    q[k] = zt[k] - cq * q[k];
                }
            }
            a.apply(&p, &mut pt);
            info.matvecs += 1;
            let eps = if symmetric { dotu(&p, &pt) } else { dotu(&q, &pt) };
            if eps == T::ZERO || !eps.is_finite() {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            let beta = eps / delta;
            if beta == T::ZERO {
                info.reason = ConvergedReason::DivergedBreakdown;
                break;
            }
            // ṽ = p̃ − β v
            vt.copy_from_slice(&pt);
            axpy(-beta, &v, &mut vt);
            m.apply(&vt, &mut y);
            info.pc_applies += 1;
            let rho_next = norm2(&y);
            let xi_next = match &left {
                Left::Symmetric => rho_next,
                Left::General { at, .. } => {
                    // w̃ = Aᵀq − βw
                    at.apply(&q, &mut wt);
                    info.matvecs += 1;
                    axpy(-beta, &w, &mut wt);
                    norm2(&wt)
                }
            };
            let theta = rho_next / (gamma_prev * beta.modulus());
            let gamma = 1.0 / (1.0 + theta * theta).sqrt();
            let eta = -eta_prev * T::from_f64(rho * gamma * gamma) / (beta * T::from_f64(gamma_prev * gamma_prev));
            if i == 1 {
                for k in 0..n {
                    d[k] = eta * p[k];
                    s[k] = eta * pt[k];
                }
            } else {
                let c = (theta_prev * gamma) * (theta_prev * gamma);
                for k in 0..n {
                    d[k] = eta * p[k] + d[k].scale(c);
                    s[k] = eta * pt[k] + s[k].scale(c);
                }
            }
            axpy(T::ONE, &d, x);
            axpy(-T::ONE, &s, &mut r);
            let rn = norm2(&r);
            info.residual_norm = rn;
            if self.opts.record_history {
                info.history.push(rn);
            }
            if let Some(reason) = conv.check(rn) {
                info.reason = reason;
                break;
            }
            rho = rho_next;
            xi = xi_next;
            gamma_prev = gamma;
            eta_prev = eta;
            theta_prev = theta;
            eps_prev = eps;
        }
        info
    }
}
