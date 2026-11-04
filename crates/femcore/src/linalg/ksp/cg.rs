use super::{ConvTest, ConvergedReason, SolveInfo, Tolerances, ensure_len};
use crate::linalg::linop::{LinearOperator, Preconditioner, Scaling};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops;

/// Which inner product CG uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CgKind {
    /// `xᴴy`: Hermitian positive definite systems (identical to the real algorithm for real scalars).
    #[default]
    Hermitian,
    /// `xᵀy` without conjugation: **COCG**, for complex *symmetric* (not Hermitian) systems.
    /// Breaks down when `pᵀAp = 0`, which is possible for indefinite complex symmetric matrices;
    /// use GMRES then. For real scalars it is the same as `Hermitian`.
    Symmetric,
}

/// Which residual norm CG monitors (PETSc `KSPSetNormType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CgNorm {
    /// `‖M⁻¹r‖₂` (PETSc default).
    #[default]
    Preconditioned,
    /// `‖r‖₂`.
    Unpreconditioned,
    /// `√|rᴴM⁻¹r|` (or the bilinear form for [`CgKind::Symmetric`]).
    Natural,
}

#[derive(Debug, Clone, Default)]
pub struct CgOptions {
    pub tol: Tolerances,
    pub kind: CgKind,
    pub norm: CgNorm,
    /// Start from the contents of `x` instead of zero.
    pub initial_guess_nonzero: bool,
    pub record_history: bool,
    /// Reproduce PETSc's breakdown test for [`CgKind::Symmetric`]: PETSc applies the "indefinite
    /// matrix" sign test to `Re(pᵀAp)` even for COCG, where `pᵀAp` is complex and its real part
    /// legitimately changes sign, so PETSc stops with `DIVERGED_INDEFINITE_MAT` on problems
    /// COCG solves fine. Off by default; on only for iteration-count comparisons with PETSc.
    pub petsc_sign_test: bool,
}

/// Preconditioned CG / COCG, reusing its work vectors across solves.
#[derive(Debug, Clone)]
pub struct Cg<T: Scalar> {
    pub opts: CgOptions,
    r: Vec<T>,
    z: Vec<T>,
    p: Vec<T>,
    w: Vec<T>,
}

impl<T: Scalar> Cg<T> {
    pub fn new(opts: CgOptions) -> Self {
        Self { opts, r: Vec::new(), z: Vec::new(), p: Vec::new(), w: Vec::new() }
    }

    /// Solves `A x = b`. `x` is the initial guess if `opts.initial_guess_nonzero`, and is
    /// overwritten with the solution in every case.
    pub fn solve<A: LinearOperator<T> + ?Sized, M: Preconditioner<T> + ?Sized>(
        &mut self,
        a: &A,
        m: &M,
        b: &[T],
        x: &mut [T],
    ) -> SolveInfo {
        let n = a.dim();
        assert!(b.len() == n && x.len() == n, "dimension mismatch");
        for v in [&mut self.r, &mut self.z, &mut self.p, &mut self.w] {
            ensure_len(v, n);
        }
        let Self { opts, r, z, p, w } = self;
        let conj = opts.kind == CgKind::Hermitian;
        let dot = |u: &[T], v: &[T]| if conj { vecops::dotc(u, v) } else { vecops::dotu(u, v) };
        let mut conv = ConvTest::new(opts.tol);
        let mut info = SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedIts,
            residual_norm: 0.0,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };

        // r ← b − A x
        if opts.initial_guess_nonzero {
            a.apply(x, r);
            info.matvecs += 1;
            vecops::aypx(-T::ONE, b, r);
        } else {
            x.fill(T::ZERO);
            r.copy_from_slice(b);
        }

        // Element-wise preconditioners are fused into the vector updates; for the identity `z`
        // aliases `r` and is never materialised.
        let diag: Option<Vec<T>> = match m.scaling() {
            Some(Scaling::Diagonal(d)) => Some(d.to_vec()),
            _ => None,
        };
        let identity = matches!(m.scaling(), Some(Scaling::Identity));

        let monitored = |norm: CgNorm, zr: T, zz: f64, rr: f64| match norm {
            CgNorm::Preconditioned => zz.sqrt(),
            CgNorm::Unpreconditioned => rr.sqrt(),
            CgNorm::Natural => zr.modulus().sqrt(),
        };
        // Initial z, beta = zᴴr (or zᵀr), monitored norm.
        let (mut beta, mut dp) = if identity {
            let rr = vecops::norm2_sq(r);
            (if conj { T::from_f64(rr) } else { vecops::dotu(r, r) }, monitored(opts.norm, T::from_f64(rr), rr, rr))
        } else {
            match &diag {
                Some(d) => {
                    for ((zi, di), ri) in z.iter_mut().zip(d).zip(r.iter()) {
                        *zi = *di * *ri;
                    }
                }
                None => m.apply(r, z),
            }
            info.pc_applies += 1;
            let (zr, zz) = vecops::dot_and_norm_sq(z, r, conj);
            (zr, monitored(opts.norm, zr, zz, vecops::norm2_sq(r)))
        };
        // With a non-zero guess PETSc measures rtol against the preconditioned/natural norm of
        // the right-hand side instead of the initial residual.
        let reference = if opts.initial_guess_nonzero {
            let mut zb = vec![T::ZERO; n];
            m.apply(b, &mut zb);
            info.pc_applies += 1;
            let snorm = match opts.norm {
                CgNorm::Preconditioned => vecops::norm2(&zb),
                CgNorm::Unpreconditioned => vecops::norm2(b),
                CgNorm::Natural => dot(&zb, b).modulus().sqrt(),
            };
            if snorm == 0.0 { dp } else { snorm }
        } else {
            dp
        };
        conv.set_reference(reference);
        info.residual_norm = dp;
        if opts.record_history {
            info.history.push(dp);
        }
        if let Some(reason) = conv.check(dp) {
            info.reason = reason;
            return info;
        }

        let (mut betaold, mut dpi_old) = (T::ONE, T::ZERO);
        let mut i = 0usize;
        loop {
            info.iterations = i + 1;
            if beta == T::ZERO {
                info.reason = ConvergedReason::Atol;
                break;
            }
            if !T::IS_COMPLEX && i > 0 && beta.re() * betaold.re() < 0.0 {
                info.reason = ConvergedReason::DivergedIndefinitePc;
                break;
            }
            let zv: &[T] = if identity { r } else { z };
            if i == 0 {
                p.copy_from_slice(zv);
            } else {
                vecops::aypx(beta / betaold, zv, p);
            }
            let dpi = if conj { a.apply_dotc(p, w) } else { a.apply_dotu(p, w) };
            info.matvecs += 1;
            betaold = beta;
            let sign_flip = i > 0 && (conj || opts.petsc_sign_test) && dpi.re() * dpi_old.re() < 0.0;
            if dpi == T::ZERO || sign_flip {
                info.reason = if conj || opts.petsc_sign_test {
                    ConvergedReason::DivergedIndefiniteMat
                } else {
                    ConvergedReason::DivergedBreakdown
                };
                break;
            }
            if !dpi.is_finite() {
                info.reason = ConvergedReason::DivergedNanOrInf;
                break;
            }
            dpi_old = dpi;
            let alpha = beta / dpi;
            if identity {
                let (zr, zz, rr) = vecops::cg_update(alpha, p, x, w, r, &mut [], None, conj);
                info.pc_applies += 1;
                beta = zr;
                dp = monitored(opts.norm, zr, zz, rr);
            } else if let Some(d) = &diag {
                let (zr, zz, rr) = vecops::cg_update(alpha, p, x, w, r, z, Some(d), conj);
                info.pc_applies += 1;
                beta = zr;
                dp = monitored(opts.norm, zr, zz, rr);
            } else {
                vecops::axpy(alpha, p, x);
                vecops::axpy(-alpha, w, r);
                m.apply(r, z);
                info.pc_applies += 1;
                let (zr, zz) = vecops::dot_and_norm_sq(z, r, conj);
                beta = zr;
                dp = monitored(opts.norm, zr, zz, vecops::norm2_sq(r));
            }
            info.residual_norm = dp;
            if opts.record_history {
                info.history.push(dp);
            }
            if let Some(reason) = conv.check(dp) {
                info.reason = reason;
                break;
            }
            i += 1;
            if i >= opts.tol.max_it {
                info.reason = ConvergedReason::DivergedIts;
                info.iterations = i;
                break;
            }
        }
        info
    }
}
