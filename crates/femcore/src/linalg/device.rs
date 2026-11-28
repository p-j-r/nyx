// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Device-agnostic vector and linear-operator abstractions, deterministic fixed-order tree
//! reductions, sliced ELLPACK (`SellCSigma`), and 64-bit sparse row matrices (`Csr64`).
//!
//! All Krylov (`cg_vec`, `fgmres_vec`) and block-eigensolver (`lobpcg_vec`) iterations keep
//! tall vectors on the target device via [`Vector`], transferring only scalar inner products or
//! small `≤ 3m × 3m` Rayleigh–Ritz matrices to the host. [`CpuVec`] delegates directly to
//! [`crate::linalg::vecops`] and [`crate::linalg::kernels`], preserving 0-ULP bit-exact parity
//! with the existing slice-based CPU solvers across all `NYX_THREADS` counts.

use crate::linalg::basis::BCHUNK;
use crate::linalg::csr::{Csr, CsrError};
use crate::linalg::dense::{Dense, hermitian_eig};
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::ksp::{
    CgKind, CgNorm, CgOptions, ConvTest, ConvergedReason, GmresOptions, Orthogonalization, PcSide, SolveInfo,
};
use crate::linalg::linop::{Identity, LinearOperator, Preconditioner};
use crate::linalg::lobpcg::{LobpcgError, LobpcgOptions, LobpcgReason};
use crate::linalg::par;
use crate::linalg::rng::random_vector;
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::{self, CHUNK};
use std::fmt;
use std::ops::{Deref, DerefMut};

// ---------------------------------------------------------------------------------------
// Vector & Device Operator Traits
// ---------------------------------------------------------------------------------------

/// Device-resident vector over [`Self::Scalar`] supporting Level-1 BLAS operations and
/// deterministic reductions.
pub trait Vector: Clone + Send + Sync {
    type Scalar: Scalar;

    /// Number of scalar elements in the vector.
    fn len(&self) -> usize;

    #[inline]
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Conjugated inner product `Σ conj(self_i) · other_i`.
    fn dotc(&self, other: &Self) -> Self::Scalar;

    /// Unconjugated bilinear form `Σ self_i · other_i` (used by COCG).
    fn dotu(&self, other: &Self) -> Self::Scalar;

    /// Squared Euclidean norm `Σ |self_i|²`.
    #[inline]
    fn norm2_sq(&self) -> f64 {
        self.dotc(self).re()
    }

    /// Euclidean norm `√(Σ |self_i|²)`.
    #[inline]
    fn norm2(&self) -> f64 {
        self.norm2_sq().sqrt()
    }

    /// `self ← self + a · x`.
    fn axpy(&mut self, a: Self::Scalar, x: &Self);

    /// `self ← a · self`.
    fn scale(&mut self, a: Self::Scalar);

    /// `self ← r · self` for a real factor `r`.
    fn scale_real(&mut self, r: f64);

    /// `self ← x + beta · self` (PETSc `VecAYPX`).
    fn aypx(&mut self, beta: Self::Scalar, x: &Self);

    /// `self ← other`.
    fn copy_from(&mut self, other: &Self);

    /// Sets every entry of `self` to zero.
    fn fill_zero(&mut self) {
        self.scale_real(0.0);
    }

    /// Copies the device vector contents to a host `Vec<Self::Scalar>`.
    fn to_host(&self) -> Vec<Self::Scalar>;

    /// Overwrites `self` from `host` (`host.len() == self.len()`).
    #[allow(clippy::wrong_self_convention)]
    fn from_host(&mut self, host: &[Self::Scalar]);

    /// Allocates a zero-initialised vector of the same length on the same device as `self`.
    fn zeros_like(&self) -> Self;

    /// `self ← self + sign · Σ_j coeffs[j] · cols[j]`.
    ///
    /// Overridden by [`CpuVec`] to use the 4-way cache-blocked combination of
    /// [`crate::linalg::basis::BasisView`] for 0-ULP parity with `Gmres`.
    fn add_combination(&mut self, cols: &[&Self], coeffs: &[Self::Scalar], sign: f64) {
        for (&v, &c) in cols.iter().zip(coeffs) {
            self.axpy(c.scale(sign), v);
        }
    }
}

/// Host-memory vector wrapper delegating every operation to [`vecops`] with zero overhead and
/// 0-ULP bit-for-bit parity.
#[derive(Debug, Clone, PartialEq)]
pub struct CpuVec<T: Scalar>(pub Vec<T>);

impl<T: Scalar> CpuVec<T> {
    #[inline]
    pub fn new(data: Vec<T>) -> Self {
        Self(data)
    }

    #[inline]
    pub fn zeros(n: usize) -> Self {
        Self(vec![T::ZERO; n])
    }

    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.0
    }

    #[inline]
    pub fn into_inner(self) -> Vec<T> {
        self.0
    }
}

impl<T: Scalar> From<Vec<T>> for CpuVec<T> {
    #[inline]
    fn from(v: Vec<T>) -> Self {
        Self(v)
    }
}

impl<T: Scalar> Deref for CpuVec<T> {
    type Target = [T];
    #[inline]
    fn deref(&self) -> &[T] {
        &self.0
    }
}

impl<T: Scalar> DerefMut for CpuVec<T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut [T] {
        &mut self.0
    }
}

impl<T: Scalar> Vector for CpuVec<T> {
    type Scalar = T;

    #[inline]
    fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    fn dotc(&self, other: &Self) -> T {
        vecops::dotc(&self.0, &other.0)
    }

    #[inline]
    fn dotu(&self, other: &Self) -> T {
        vecops::dotu(&self.0, &other.0)
    }

    #[inline]
    fn norm2_sq(&self) -> f64 {
        vecops::norm2_sq(&self.0)
    }

    #[inline]
    fn norm2(&self) -> f64 {
        vecops::norm2(&self.0)
    }

    #[inline]
    fn axpy(&mut self, a: T, x: &Self) {
        vecops::axpy(a, &x.0, &mut self.0);
    }

    #[inline]
    fn scale(&mut self, a: T) {
        vecops::scale(a, &mut self.0);
    }

    #[inline]
    fn scale_real(&mut self, r: f64) {
        vecops::scale_real(r, &mut self.0);
    }

    #[inline]
    fn aypx(&mut self, beta: T, x: &Self) {
        vecops::aypx(beta, &x.0, &mut self.0);
    }

    #[inline]
    fn copy_from(&mut self, other: &Self) {
        self.0.copy_from_slice(&other.0);
    }

    #[inline]
    fn fill_zero(&mut self) {
        self.0.fill(T::ZERO);
    }

    #[inline]
    fn to_host(&self) -> Vec<T> {
        self.0.clone()
    }

    #[inline]
    fn from_host(&mut self, host: &[T]) {
        if self.0.len() == host.len() {
            self.0.copy_from_slice(host);
        } else {
            self.0.clear();
            self.0.extend_from_slice(host);
        }
    }

    #[inline]
    fn zeros_like(&self) -> Self {
        Self(vec![T::ZERO; self.0.len()])
    }

    fn add_combination(&mut self, cols: &[&Self], coeffs: &[T], sign: f64) {
        let k = cols.len().min(coeffs.len());
        if k == 0 {
            return;
        }
        let n = self.0.len();
        let nch = n.div_ceil(BCHUNK);
        let yc = par::RawCols::new(&mut self.0, n, 1);
        par::chunk_map(nch, par::worth_parallel(n * k), |ch| {
            let (s, e) = (ch * BCHUNK, ((ch + 1) * BCHUNK).min(n));
            let dst = yc.chunk_mut(0, s, e);
            let len = dst.len();
            let mut j = 0usize;
            while j + 4 <= k {
                let (c0, c1, c2, c3) = (
                    coeffs[j].scale(sign),
                    coeffs[j + 1].scale(sign),
                    coeffs[j + 2].scale(sign),
                    coeffs[j + 3].scale(sign),
                );
                let (v0, v1, v2, v3) = (
                    &cols[j].0[s..s + len],
                    &cols[j + 1].0[s..s + len],
                    &cols[j + 2].0[s..s + len],
                    &cols[j + 3].0[s..s + len],
                );
                for i in 0..len {
                    dst[i] += (c0 * v0[i] + c1 * v1[i]) + (c2 * v2[i] + c3 * v3[i]);
                }
                j += 4;
            }
            while j < k {
                T::k_axpy(coeffs[j].scale(sign), &cols[j].0[s..s + len], dst);
                j += 1;
            }
        });
    }
}

/// Square linear operator `y ← A x` acting on device-resident vectors `V`.
pub trait DeviceLinearOperator<V: Vector> {
    /// Dimension `n` of the square operator.
    fn dim(&self) -> usize;

    /// `y ← A x`.
    fn apply(&self, x: &V, y: &mut V);

    /// `y ← A x`, returning `xᴴ y`.
    fn apply_dotc(&self, x: &V, y: &mut V) -> V::Scalar {
        self.apply(x, y);
        x.dotc(y)
    }

    /// `y ← A x`, returning `xᵀ y` without conjugation.
    fn apply_dotu(&self, x: &V, y: &mut V) -> V::Scalar {
        self.apply(x, y);
        x.dotu(y)
    }
}

impl<V: Vector, O: DeviceLinearOperator<V> + ?Sized> DeviceLinearOperator<V> for &O {
    #[inline]
    fn dim(&self) -> usize {
        (**self).dim()
    }
    #[inline]
    fn apply(&self, x: &V, y: &mut V) {
        (**self).apply(x, y);
    }
    #[inline]
    fn apply_dotc(&self, x: &V, y: &mut V) -> V::Scalar {
        (**self).apply_dotc(x, y)
    }
    #[inline]
    fn apply_dotu(&self, x: &V, y: &mut V) -> V::Scalar {
        (**self).apply_dotu(x, y)
    }
}

/// Adapter wrapping any host [`LinearOperator<T>`] as a [`DeviceLinearOperator<CpuVec<T>>`].
pub struct HostOp<'a, O: ?Sized>(pub &'a O);

impl<T: Scalar, O: LinearOperator<T> + ?Sized> DeviceLinearOperator<CpuVec<T>> for HostOp<'_, O> {
    #[inline]
    fn dim(&self) -> usize {
        self.0.dim()
    }
    #[inline]
    fn apply(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) {
        self.0.apply(&x.0, &mut y.0);
    }
    #[inline]
    fn apply_dotc(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) -> T {
        self.0.apply_dotc(&x.0, &mut y.0)
    }
    #[inline]
    fn apply_dotu(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) -> T {
        self.0.apply_dotu(&x.0, &mut y.0)
    }
}

impl<T: Scalar> DeviceLinearOperator<CpuVec<T>> for Csr<T> {
    #[inline]
    fn dim(&self) -> usize {
        LinearOperator::dim(self)
    }
    #[inline]
    fn apply(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) {
        self.spmv(&x.0, &mut y.0);
    }
    #[inline]
    fn apply_dotc(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) -> T {
        self.spmv_dotc(&x.0, &mut y.0)
    }
}

/// Preconditioner `z ← M⁻¹ r` acting on device-resident vectors `V`.
pub trait DevicePreconditioner<V: Vector> {
    fn apply(&self, r: &V, z: &mut V);
}

impl<V: Vector, M: DevicePreconditioner<V> + ?Sized> DevicePreconditioner<V> for &M {
    #[inline]
    fn apply(&self, r: &V, z: &mut V) {
        (**self).apply(r, z);
    }
}

/// Identity preconditioner `z ← r` for any [`Vector`].
#[derive(Debug, Clone, Copy, Default)]
pub struct DeviceIdentity;

impl<V: Vector> DevicePreconditioner<V> for DeviceIdentity {
    #[inline]
    fn apply(&self, r: &V, z: &mut V) {
        z.copy_from(r);
    }
}

impl<T: Scalar> DevicePreconditioner<CpuVec<T>> for Identity {
    #[inline]
    fn apply(&self, r: &CpuVec<T>, z: &mut CpuVec<T>) {
        z.copy_from(r);
    }
}

/// Adapter wrapping any host [`Preconditioner<T>`] as a [`DevicePreconditioner<CpuVec<T>>`].
pub struct HostPc<'a, M: ?Sized>(pub &'a M);

impl<T: Scalar, M: Preconditioner<T> + ?Sized> DevicePreconditioner<CpuVec<T>> for HostPc<'_, M> {
    #[inline]
    fn apply(&self, r: &CpuVec<T>, z: &mut CpuVec<T>) {
        self.0.apply(&r.0, &mut z.0);
    }
}

// ---------------------------------------------------------------------------------------
// Generic Device-Resident Krylov & Block Eigensolvers
// ---------------------------------------------------------------------------------------

/// Preconditioned Conjugate Gradient / COCG over device-resident [`Vector`]s.
///
/// On [`CpuVec<T>`] this reproduces [`crate::linalg::ksp::Cg::solve`] bit-for-bit (`0-ULP`).
pub fn cg_vec<V, A, M>(a: &A, m: &M, b: &V, x: &mut V, opts: &CgOptions) -> SolveInfo
where
    V: Vector,
    A: DeviceLinearOperator<V> + ?Sized,
    M: DevicePreconditioner<V> + ?Sized,
{
    type S<V> = <V as Vector>::Scalar;
    let n = a.dim();
    if b.len() != n || x.len() != n {
        return SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedBreakdown,
            residual_norm: f64::NAN,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };
    }

    let conj = opts.kind == CgKind::Hermitian;
    let dot = |u: &V, v: &V| if conj { u.dotc(v) } else { u.dotu(v) };
    let monitored = |norm: CgNorm, zr: S<V>, zz: f64, rr: f64| match norm {
        CgNorm::Preconditioned => zz.sqrt(),
        CgNorm::Unpreconditioned => rr.sqrt(),
        CgNorm::Natural => zr.modulus().sqrt(),
    };

    let mut conv = ConvTest::new(opts.tol);
    let mut info = SolveInfo {
        iterations: 0,
        reason: ConvergedReason::DivergedIts,
        residual_norm: 0.0,
        history: Vec::new(),
        matvecs: 0,
        pc_applies: 0,
    };

    let mut r = b.zeros_like();
    let mut z = b.zeros_like();
    let mut p = b.zeros_like();
    let mut w = b.zeros_like();

    if opts.initial_guess_nonzero {
        a.apply(x, &mut r);
        info.matvecs += 1;
        r.aypx(-S::<V>::ONE, b);
    } else {
        x.fill_zero();
        r.copy_from(b);
    }

    m.apply(&r, &mut z);
    info.pc_applies += 1;
    let zr0 = dot(&z, &r);
    let zz0 = z.norm2_sq();
    let rr0 = r.norm2_sq();
    let mut beta = zr0;
    let mut dp = monitored(opts.norm, zr0, zz0, rr0);

    let reference = if opts.initial_guess_nonzero {
        let mut zb = b.zeros_like();
        m.apply(b, &mut zb);
        info.pc_applies += 1;
        let snorm = match opts.norm {
            CgNorm::Preconditioned => zb.norm2(),
            CgNorm::Unpreconditioned => b.norm2(),
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

    let mut betaold = S::<V>::ONE;
    let mut dpi_old = S::<V>::ZERO;
    let mut i = 0usize;
    loop {
        info.iterations = i + 1;
        if beta == S::<V>::ZERO {
            info.reason = ConvergedReason::Atol;
            break;
        }
        if !S::<V>::IS_COMPLEX && i > 0 && beta.re() * betaold.re() < 0.0 {
            info.reason = ConvergedReason::DivergedIndefinitePc;
            break;
        }
        if i == 0 {
            p.copy_from(&z);
        } else {
            p.aypx(beta / betaold, &z);
        }
        let dpi = if conj { a.apply_dotc(&p, &mut w) } else { a.apply_dotu(&p, &mut w) };
        info.matvecs += 1;
        betaold = beta;
        let sign_flip = i > 0 && (conj || opts.petsc_sign_test) && dpi.re() * dpi_old.re() < 0.0;
        if dpi == S::<V>::ZERO || sign_flip {
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
        x.axpy(alpha, &p);
        r.axpy(-alpha, &w);
        m.apply(&r, &mut z);
        info.pc_applies += 1;
        let zr = dot(&z, &r);
        let zz = z.norm2_sq();
        let rr = r.norm2_sq();
        beta = zr;
        dp = monitored(opts.norm, zr, zz, rr);
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

const HAPPY_TOL: f64 = 1e-30;
const BREAKDOWN_TOL: f64 = 0.1;

fn givens_vec<T: Scalar>(a: T, b: T) -> (f64, T) {
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

fn orthogonalize_vec<V: Vector>(v: &[V], w: &mut V, k: usize, kind: Orthogonalization, h: &mut [V::Scalar]) -> f64 {
    let cols: Vec<&V> = v.iter().take(k + 1).collect();
    match kind {
        Orthogonalization::ClassicalGramSchmidt => {
            for j in 0..=k {
                h[j] = v[j].dotc(w);
            }
            w.add_combination(&cols, &h[..=k], -1.0);
            w.norm2()
        }
        Orthogonalization::ClassicalGramSchmidtIfNeeded => {
            let before = w.norm2();
            for j in 0..=k {
                h[j] = v[j].dotc(w);
            }
            w.add_combination(&cols, &h[..=k], -1.0);
            let mut after = w.norm2();
            if after < std::f64::consts::FRAC_1_SQRT_2 * before {
                let mut h2 = vec![V::Scalar::ZERO; k + 1];
                for j in 0..=k {
                    h2[j] = v[j].dotc(w);
                }
                w.add_combination(&cols, &h2, -1.0);
                for (a, b) in h[..=k].iter_mut().zip(&h2) {
                    *a += *b;
                }
                after = w.norm2();
            }
            after
        }
        Orthogonalization::ModifiedGramSchmidt => {
            for (j, hj) in h.iter_mut().enumerate().take(k + 1) {
                let c = v[j].dotc(w);
                *hj = c;
                w.axpy(-c, &v[j]);
            }
            w.norm2()
        }
    }
}

/// Restarted GMRES(m) / Flexible GMRES (FGMRES) over device-resident [`Vector`]s.
///
/// On [`CpuVec<T>`] this reproduces [`crate::linalg::ksp::Gmres::solve`] bit-for-bit (`0-ULP`).
pub fn fgmres_vec<V, A, M>(a: &A, m: &M, b: &V, x: &mut V, opts: &GmresOptions) -> SolveInfo
where
    V: Vector,
    A: DeviceLinearOperator<V> + ?Sized,
    M: DevicePreconditioner<V> + ?Sized,
{
    type S<V> = <V as Vector>::Scalar;
    let n = a.dim();
    if b.len() != n || x.len() != n {
        return SolveInfo {
            iterations: 0,
            reason: ConvergedReason::DivergedBreakdown,
            residual_norm: f64::NAN,
            history: Vec::new(),
            matvecs: 0,
            pc_applies: 0,
        };
    }

    let restart = opts.restart.max(1);
    let flexible = opts.flexible;
    let left = opts.side == PcSide::Left && !flexible;

    let mut v: Vec<V> = (0..restart + 1).map(|_| b.zeros_like()).collect();
    let mut zb: Vec<V> = if flexible { (0..restart).map(|_| b.zeros_like()).collect() } else { Vec::new() };
    let mut r = b.zeros_like();
    let mut t = b.zeros_like();
    let mut u = b.zeros_like();

    let ld = restart + 1;
    let mut h = vec![S::<V>::ZERO; ld * restart];
    let mut cs = vec![0.0f64; restart];
    let mut sn = vec![S::<V>::ZERO; restart];
    let mut rs = vec![S::<V>::ZERO; ld];

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
        x.fill_zero();
    }
    let mut first_cycle = true;
    let mut last_rnorm = -1.0f64;
    let mut its = 0usize;

    'outer: loop {
        if guess_nonzero {
            a.apply(x, &mut r);
            info.matvecs += 1;
            r.aypx(-S::<V>::ONE, b);
        } else {
            r.copy_from(b);
        }
        if left {
            m.apply(&r, &mut v[0]);
            info.pc_applies += 1;
        } else {
            v[0].copy_from(&r);
        }
        let res = v[0].norm2();
        if !res.is_finite() {
            info.reason = ConvergedReason::DivergedNanOrInf;
            info.residual_norm = res;
            break;
        }
        if first_cycle {
            let reference = if guess_nonzero {
                let s = if left {
                    m.apply(b, &mut t);
                    info.pc_applies += 1;
                    t.norm2()
                } else {
                    b.norm2()
                };
                if s == 0.0 { res } else { s }
            } else {
                res
            };
            conv.set_reference(reference);
            if opts.record_history {
                info.history.push(res);
            }
        } else if last_rnorm > 0.0 && (res - last_rnorm).abs() > BREAKDOWN_TOL * conv.rnorm0().max(f64::MIN_POSITIVE) {
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
        v[0].scale_real(1.0 / res);
        rs.fill(S::<V>::ZERO);
        rs[0] = S::<V>::from_f64(res);

        let mut k = 0usize;
        let mut reason = None;
        let mut res_k = res;
        while k < restart && its < opts.tol.max_it {
            {
                let (head, tail) = v.split_at_mut(k + 1);
                let vk = &head[k];
                let w = &mut tail[0];
                if left {
                    a.apply(vk, &mut t);
                    m.apply(&t, w);
                } else if flexible {
                    m.apply(vk, &mut zb[k]);
                    a.apply(&zb[k], w);
                } else {
                    m.apply(vk, &mut t);
                    a.apply(&t, w);
                }
                info.matvecs += 1;
                info.pc_applies += 1;
                let hcol = &mut h[k * ld..(k + 1) * ld];
                let tt = orthogonalize_vec(head, w, k, opts.orthog, hcol);
                if !tt.is_finite() {
                    reason = Some(ConvergedReason::DivergedNanOrInf);
                    break;
                }
                if tt != 0.0 {
                    w.scale_real(1.0 / tt);
                }
                let hapbnd = (tt / rs[k].modulus()).min(HAPPY_TOL);
                let happy = tt < hapbnd;
                for j in 0..k {
                    let (c, s) = (cs[j], sn[j]);
                    let (hj, hj1) = (hcol[j], hcol[j + 1]);
                    hcol[j] = hj.scale(c) + s * hj1;
                    hcol[j + 1] = hj1.scale(c) - s.conj() * hj;
                }
                hcol[k + 1] = S::<V>::from_f64(tt);
                if !happy {
                    let (c, s) = givens_vec(hcol[k], hcol[k + 1]);
                    if c == 1.0 && hcol[k] == S::<V>::ZERO && hcol[k + 1] == S::<V>::ZERO {
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
                hcol[k + 1] = S::<V>::ZERO;
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
        }

        if k > 0 {
            let mut y = vec![S::<V>::ZERO; k];
            for i in (0..k).rev() {
                let mut s = rs[i];
                for j in i + 1..k {
                    s -= h[j * ld + i] * y[j];
                }
                let d = h[i * ld + i];
                if d == S::<V>::ZERO {
                    reason = Some(ConvergedReason::DivergedBreakdown);
                } else {
                    y[i] = s / d;
                }
            }
            if left {
                let cols: Vec<&V> = v.iter().take(k).collect();
                x.add_combination(&cols, &y, 1.0);
            } else if flexible {
                let cols: Vec<&V> = zb.iter().take(k).collect();
                x.add_combination(&cols, &y, 1.0);
            } else {
                t.fill_zero();
                let cols: Vec<&V> = v.iter().take(k).collect();
                t.add_combination(&cols, &y, 1.0);
                m.apply(&t, &mut u);
                info.pc_applies += 1;
                x.axpy(S::<V>::ONE, &u);
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

/// Result of [`lobpcg_vec`], keeping converged Ritz vectors on the target device `V`.
#[derive(Debug, Clone)]
pub struct LobpcgVecResult<V: Vector> {
    pub values: Vec<f64>,
    pub vectors: Vec<V>,
    pub residuals: Vec<f64>,
    pub iterations: usize,
    pub reason: LobpcgReason,
    pub matvecs: usize,
    pub active_history: Vec<usize>,
}

fn apply_b_vec<V, B>(b: Option<&B>, x: &V) -> V
where
    V: Vector,
    B: DeviceLinearOperator<V> + ?Sized,
{
    match b {
        Some(op) => {
            let mut y = x.zeros_like();
            op.apply(x, &mut y);
            y
        }
        None => x.clone(),
    }
}

fn orth_column_vec<V: Vector>(v: &mut V, bv: &mut V, against: &[&Vec<(V, V)>]) -> bool {
    let n0 = v.norm2();
    if n0 == 0.0 {
        return false;
    }
    for _ in 0..2 {
        for set in against {
            for (u, bu) in set.iter() {
                let h = bu.dotc(v);
                v.axpy(-h, u);
                bv.axpy(-h, bu);
            }
        }
    }
    let nb = v.dotc(bv).re();
    if nb.is_nan() || nb <= 1e-20 * n0 * n0 {
        return false;
    }
    let s = 1.0 / nb.sqrt();
    v.scale_real(s);
    bv.scale_real(s);
    true
}

fn hermitise_dense<T: FaerScalar>(g: &Dense<T>, sign: f64) -> Dense<T> {
    let d = g.rows();
    Dense::from_fn(d, d, |i, j| (g[(i, j)] + g[(j, i)].conj()).scale(0.5 * sign))
}

fn rotate_vec<V: Vector>(x: &mut [(V, V)], ax: &mut [V], q: &Dense<V::Scalar>, m: usize) {
    let mix = |vecs: &[&V], c: usize| -> V {
        let mut out = vecs[0].zeros_like();
        for (i, v) in vecs.iter().enumerate() {
            out.axpy(q[(i, c)], v);
        }
        out
    };
    let xs: Vec<&V> = x.iter().map(|e| &e.0).collect();
    let bxs: Vec<&V> = x.iter().map(|e| &e.1).collect();
    let axs: Vec<&V> = ax.iter().collect();
    let new_x: Vec<V> = (0..m).map(|c| mix(&xs, c)).collect();
    let new_bx: Vec<V> = (0..m).map(|c| mix(&bxs, c)).collect();
    let new_ax: Vec<V> = (0..m).map(|c| mix(&axs, c)).collect();
    for c in 0..m {
        x[c] = (new_x[c].clone(), new_bx[c].clone());
        ax[c] = new_ax[c].clone();
    }
}

/// Device-resident LOBPCG for `A x = λ B x`.
///
/// Tall-skinny blocks (`X`, `W`, `P`, `A·X`, `B·X`) remain resident on the device `V`, while the
/// small `≤ 3m × 3m` Rayleigh–Ritz subspace eigensolve runs on the host via `faer`.
pub fn lobpcg_vec<V, A, B, M>(
    a: &A,
    b: Option<&B>,
    pc: Option<&M>,
    proto: &V,
    x0: Option<&[V]>,
    constraints: &[V],
    opts: &LobpcgOptions,
) -> Result<LobpcgVecResult<V>, LobpcgError>
where
    V: Vector,
    V::Scalar: FaerScalar,
    A: DeviceLinearOperator<V> + ?Sized,
    B: DeviceLinearOperator<V> + ?Sized,
    M: DevicePreconditioner<V> + ?Sized,
{
    let n = a.dim();
    if proto.len() != n {
        return Err(LobpcgError::InvalidArgument("prototype vector length does not match operator dimension"));
    }
    if let Some(b_op) = b {
        if b_op.dim() != n {
            return Err(LobpcgError::InvalidArgument("mass operator dimension does not match A"));
        }
    }
    let m = opts.block.max(opts.nev);
    if opts.nev == 0 || 3 * m + constraints.len() >= n {
        return Err(LobpcgError::InvalidArgument("need nev >= 1, a square operator and 3*block + constraints < n"));
    }
    let sign = if opts.which == crate::linalg::lobpcg::LobpcgWhich::Lowest { 1.0 } else { -1.0 };

    let mut cons: Vec<(V, V)> = Vec::new();
    for y in constraints {
        if y.len() != n {
            return Err(LobpcgError::InvalidArgument("constraint vector length mismatch"));
        }
        let mut v = y.clone();
        let mut bv = apply_b_vec(b, y);
        if orth_column_vec(&mut v, &mut bv, &[&cons]) {
            cons.push((v, bv));
        }
    }

    let mut x: Vec<(V, V)> = Vec::new();
    let mut seed = opts.seed;
    for j in 0..m {
        loop {
            let mut v = match x0.and_then(|s| s.get(j)) {
                Some(s) if (seed == opts.seed || x.len() > j) && s.len() == n => s.clone(),
                _ => {
                    let host = random_vector::<V::Scalar>(n, seed);
                    let mut tmp = proto.zeros_like();
                    tmp.from_host(&host);
                    tmp
                }
            };
            seed = seed.wrapping_add(0x9E37);
            let mut bv = apply_b_vec(b, &v);
            if orth_column_vec(&mut v, &mut bv, &[&cons, &x]) {
                x.push((v, bv));
                break;
            }
        }
    }

    let mut ax: Vec<V> = x
        .iter()
        .map(|(v, _)| {
            let mut y = proto.zeros_like();
            a.apply(v, &mut y);
            y
        })
        .collect();
    let mut matvecs = m;
    let mut lambda = vec![0.0f64; m];

    {
        let g = Dense::from_fn(m, m, |i, j| x[i].0.dotc(&ax[j]));
        let g = hermitise_dense(&g, sign);
        let (w, q) = hermitian_eig(&g)?;
        rotate_vec(&mut x, &mut ax, &q, m);
        for i in 0..m {
            lambda[i] = sign * w[i];
        }
    }

    let mut p_dirs: Vec<(V, V, V)> = Vec::new();
    let mut active_history = Vec::new();
    let (mut best_worst, mut since_improve) = (f64::INFINITY, 0usize);
    let mut iterations = 0usize;
    let mut reason = LobpcgReason::MaxIterations;
    let mut residual = vec![0.0f64; m];

    for it in 0..opts.max_it {
        let mut r: Vec<V> = Vec::with_capacity(m);
        for j in 0..m {
            let mut rj = ax[j].clone();
            rj.axpy(V::Scalar::from_f64(-lambda[j]), &x[j].1);
            let scale = ax[j].norm2() + lambda[j].abs() * x[j].1.norm2();
            residual[j] = rj.norm2() / (scale + 1e-30);
            r.push(rj);
        }
        let wanted_ok = residual.iter().take(opts.nev).all(|&e| e <= opts.tol);
        let active: Vec<usize> = (0..m).filter(|&j| residual[j] > opts.tol).collect();
        active_history.push(active.len());
        iterations = it;
        if wanted_ok {
            reason = LobpcgReason::Converged;
            break;
        }
        let worst = active.iter().filter(|&&j| j < opts.nev).map(|&j| residual[j]).fold(0.0, f64::max);
        if worst < (1.0 - opts.stall_gain) * best_worst {
            best_worst = worst;
            since_improve = 0;
        } else {
            since_improve += 1;
            if since_improve >= opts.stall_window {
                reason = LobpcgReason::Stalled;
                break;
            }
        }

        let mut w_dirs: Vec<(V, V, V)> = Vec::new();
        let mut w_basis: Vec<(V, V)> = Vec::new();
        for &j in &active {
            let mut wj = proto.zeros_like();
            match pc {
                Some(p_op) => p_op.apply(&r[j], &mut wj),
                None => wj.copy_from(&r[j]),
            }
            let mut bw = apply_b_vec(b, &wj);
            if orth_column_vec(&mut wj, &mut bw, &[&cons, &x, &w_basis]) {
                let mut aw = proto.zeros_like();
                a.apply(&wj, &mut aw);
                matvecs += 1;
                w_basis.push((wj.clone(), bw.clone()));
                w_dirs.push((wj, bw, aw));
            }
        }

        let mut p_new: Vec<(V, V, V)> = Vec::new();
        let mut p_basis: Vec<(V, V)> = Vec::new();
        for (pv, _, _) in p_dirs.iter() {
            let mut v = pv.clone();
            let mut bv = apply_b_vec(b, &v);
            if orth_column_vec(&mut v, &mut bv, &[&cons, &x, &w_basis, &p_basis]) {
                let mut av = proto.zeros_like();
                a.apply(&v, &mut av);
                matvecs += 1;
                p_basis.push((v.clone(), bv.clone()));
                p_new.push((v, bv, av));
            }
        }

        let mut s_vecs: Vec<&V> = x.iter().map(|e| &e.0).collect();
        let mut as_vecs: Vec<&V> = ax.iter().collect();
        for (v, _, av) in w_dirs.iter().chain(p_new.iter()) {
            s_vecs.push(v);
            as_vecs.push(av);
        }
        let d = s_vecs.len();
        let g = Dense::from_fn(d, d, |i, j| s_vecs[i].dotc(as_vecs[j]));
        let g = hermitise_dense(&g, sign);
        let (w, q) = hermitian_eig(&g)?;

        let cols: Vec<Vec<V::Scalar>> = (0..m).map(|c| (0..d).map(|i| q[(i, c)]).collect()).collect();
        let combine = |vecs: &[&V], coef: &[V::Scalar], from: usize| -> V {
            let mut out = proto.zeros_like();
            for (k, v) in vecs.iter().enumerate().skip(from) {
                out.axpy(coef[k], v);
            }
            out
        };
        let bs_vecs: Vec<&V> =
            x.iter().map(|e| &e.1).chain(w_dirs.iter().map(|e| &e.1)).chain(p_new.iter().map(|e| &e.1)).collect();
        let mut new_x = Vec::with_capacity(m);
        let mut new_ax = Vec::with_capacity(m);
        let mut new_p = Vec::new();
        for col in cols.iter().take(m) {
            let xv = combine(&s_vecs, col, 0);
            let bxv = combine(&bs_vecs, col, 0);
            let axv = combine(&as_vecs, col, 0);
            new_x.push((xv, bxv));
            new_ax.push(axv);
            if m < d {
                let pv = combine(&s_vecs, col, m);
                let bpv = combine(&bs_vecs, col, m);
                let apv = combine(&as_vecs, col, m);
                new_p.push((pv, bpv, apv));
            }
        }
        for i in 0..m {
            lambda[i] = sign * w[i];
        }
        x = new_x;
        ax = new_ax;
        let act_after = active;
        p_dirs = new_p.into_iter().enumerate().filter(|(j, _)| act_after.contains(j)).map(|(_, p)| p).collect();
    }

    let nev = opts.nev;
    let mut vectors = Vec::with_capacity(nev);
    let mut residuals = Vec::with_capacity(nev);
    for j in 0..nev {
        let mut rj = ax[j].clone();
        rj.axpy(V::Scalar::from_f64(-lambda[j]), &x[j].1);
        residuals.push(rj.norm2() / (ax[j].norm2() + lambda[j].abs() * x[j].1.norm2() + 1e-30));
        vectors.push(x[j].0.clone());
    }
    Ok(LobpcgVecResult {
        values: lambda[..nev].to_vec(),
        vectors,
        residuals,
        iterations,
        reason,
        matvecs,
        active_history,
    })
}

// ---------------------------------------------------------------------------------------
// Deterministic Fixed-Order Two-Level Tree Reductions (Zero Floating-Point Atomics)
// ---------------------------------------------------------------------------------------

/// Fixed workgroup block width of the deterministic tree reductions (`BLOCK = 256`).
pub const TREE_REDUCE_BLOCK: usize = 256;

#[inline]
fn reduce_block_256<T: Scalar>(slice: &[T]) -> T {
    let mut buf = [T::ZERO; TREE_REDUCE_BLOCK];
    for (dst, &src) in buf.iter_mut().zip(slice) {
        *dst = src;
    }
    let mut s = TREE_REDUCE_BLOCK / 2;
    while s > 0 {
        for tid in 0..s {
            let v = buf[tid + s];
            buf[tid] += v;
        }
        s >>= 1;
    }
    buf[0]
}

#[inline]
fn reduce_block_256_f32(slice: &[f32]) -> f32 {
    let mut buf = [0.0f32; TREE_REDUCE_BLOCK];
    for (dst, &src) in buf.iter_mut().zip(slice) {
        *dst = src;
    }
    let mut s = TREE_REDUCE_BLOCK / 2;
    while s > 0 {
        for tid in 0..s {
            buf[tid] += buf[tid + s];
        }
        s >>= 1;
    }
    buf[0]
}

/// Deterministic hierarchical binary-tree sum over `BLOCK = 256` tiles with zero floating-point
/// atomics. Identical across all `NYX_THREADS` counts.
pub fn tree_reduce_sum<T: Scalar>(xs: &[T]) -> T {
    if xs.is_empty() {
        return T::ZERO;
    }
    let nblocks = xs.len().div_ceil(TREE_REDUCE_BLOCK);
    let mut partials = par::chunk_map(nblocks, par::worth_parallel(xs.len()), |b| {
        let s = b * TREE_REDUCE_BLOCK;
        let e = (s + TREE_REDUCE_BLOCK).min(xs.len());
        reduce_block_256(&xs[s..e])
    });
    while partials.len() > 1 {
        let nb = partials.len().div_ceil(TREE_REDUCE_BLOCK);
        partials = (0..nb)
            .map(|b| {
                let s = b * TREE_REDUCE_BLOCK;
                let e = (s + TREE_REDUCE_BLOCK).min(partials.len());
                reduce_block_256(&partials[s..e])
            })
            .collect();
    }
    partials[0]
}

/// Deterministic hierarchical binary-tree inner product `Σ conj(x_i) · y_i` (`BLOCK = 256`, zero atomics).
pub fn tree_reduce_dotc<T: Scalar>(x: &[T], y: &[T]) -> T {
    let n = x.len().min(y.len());
    if n == 0 {
        return T::ZERO;
    }
    let nblocks = n.div_ceil(TREE_REDUCE_BLOCK);
    let mut partials = par::chunk_map(nblocks, par::worth_parallel(n), |b| {
        let s = b * TREE_REDUCE_BLOCK;
        let e = (s + TREE_REDUCE_BLOCK).min(n);
        let mut buf = [T::ZERO; TREE_REDUCE_BLOCK];
        for (tid, (&xi, &yi)) in x[s..e].iter().zip(&y[s..e]).enumerate() {
            buf[tid] = xi.conj() * yi;
        }
        let mut stride = TREE_REDUCE_BLOCK / 2;
        while stride > 0 {
            for tid in 0..stride {
                let v = buf[tid + stride];
                buf[tid] += v;
            }
            stride >>= 1;
        }
        buf[0]
    });
    while partials.len() > 1 {
        let nb = partials.len().div_ceil(TREE_REDUCE_BLOCK);
        partials = (0..nb)
            .map(|b| {
                let s = b * TREE_REDUCE_BLOCK;
                let e = (s + TREE_REDUCE_BLOCK).min(partials.len());
                reduce_block_256(&partials[s..e])
            })
            .collect();
    }
    partials[0]
}

/// Single-precision `f32` hierarchical binary-tree sum (`BLOCK = 256`, zero atomics), matching the
/// exact reduction schedule of the WGSL GPU compute shader.
pub fn tree_reduce_sum_f32(xs: &[f32]) -> f32 {
    if xs.is_empty() {
        return 0.0;
    }
    let nblocks = xs.len().div_ceil(TREE_REDUCE_BLOCK);
    let mut partials = par::chunk_map(nblocks, par::worth_parallel(xs.len()), |b| {
        let s = b * TREE_REDUCE_BLOCK;
        let e = (s + TREE_REDUCE_BLOCK).min(xs.len());
        reduce_block_256_f32(&xs[s..e])
    });
    while partials.len() > 1 {
        let nb = partials.len().div_ceil(TREE_REDUCE_BLOCK);
        partials = (0..nb)
            .map(|b| {
                let s = b * TREE_REDUCE_BLOCK;
                let e = (s + TREE_REDUCE_BLOCK).min(partials.len());
                reduce_block_256_f32(&partials[s..e])
            })
            .collect();
    }
    partials[0]
}

/// Single-precision `f32` hierarchical binary-tree dot product (`BLOCK = 256`, zero atomics),
/// matching the WGSL `tree_reduce_dot_l1` + `tree_reduce_sum_level` compute shaders.
pub fn tree_reduce_dot_f32(x: &[f32], y: &[f32]) -> f32 {
    let n = x.len().min(y.len());
    if n == 0 {
        return 0.0;
    }
    let nblocks = n.div_ceil(TREE_REDUCE_BLOCK);
    let mut partials = par::chunk_map(nblocks, par::worth_parallel(n), |b| {
        let s = b * TREE_REDUCE_BLOCK;
        let e = (s + TREE_REDUCE_BLOCK).min(n);
        let mut buf = [0.0f32; TREE_REDUCE_BLOCK];
        for (tid, (&xi, &yi)) in x[s..e].iter().zip(&y[s..e]).enumerate() {
            buf[tid] = xi * yi;
        }
        let mut stride = TREE_REDUCE_BLOCK / 2;
        while stride > 0 {
            for tid in 0..stride {
                buf[tid] += buf[tid + stride];
            }
            stride >>= 1;
        }
        buf[0]
    });
    while partials.len() > 1 {
        let nb = partials.len().div_ceil(TREE_REDUCE_BLOCK);
        partials = (0..nb)
            .map(|b| {
                let s = b * TREE_REDUCE_BLOCK;
                let e = (s + TREE_REDUCE_BLOCK).min(partials.len());
                reduce_block_256_f32(&partials[s..e])
            })
            .collect();
    }
    partials[0]
}

// ---------------------------------------------------------------------------------------
// SELL-C-sigma Sparse Matrix Format (C = 32, sigma = 128)
// ---------------------------------------------------------------------------------------

/// Default warp-aligned slice height `C = 32` for [`SellCSigma`].
pub const DEFAULT_SELL_C: usize = 32;
/// Default local row-sorting window `σ = 128` for [`SellCSigma`].
pub const DEFAULT_SELL_SIGMA: usize = 128;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SellCSigmaError {
    InvalidParams(&'static str),
    DimensionTooLarge,
}

impl fmt::Display for SellCSigmaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParams(m) => write!(f, "invalid SELL-C-sigma parameters: {m}"),
            Self::DimensionTooLarge => write!(f, "matrix dimension or storage exceeds u32 range"),
        }
    }
}

impl std::error::Error for SellCSigmaError {}

/// Sliced ELLPACK with local `σ`-window row-length sorting (`SELL-C-σ`, Kreutzer et al. 2014).
///
/// Within each contiguous window of `sigma` rows (a multiple of `c`), rows are stably sorted in
/// descending order of `nnz` so every `c`-row warp slice has minimal padding, while `perm` stays
/// strictly confined to its `sigma`-window so parallel window chunks write to disjoint slices of `y`.
#[derive(Debug, Clone, PartialEq)]
pub struct SellCSigma<T: Scalar> {
    nrows: usize,
    ncols: usize,
    nnz: usize,
    c: usize,
    sigma: usize,
    slice_ptr: Vec<usize>,
    col_idx: Vec<u32>,
    values: Vec<T>,
    row_lengths: Vec<u32>,
    perm: Vec<u32>,
    inv_perm: Vec<u32>,
}

impl<T: Scalar> SellCSigma<T> {
    /// Converts `csr` to `SELL-C-σ` with the default `C = 32, σ = 128` parameters.
    pub fn from_csr(csr: &Csr<T>) -> Result<Self, SellCSigmaError> {
        Self::from_csr_with(csr, DEFAULT_SELL_C, DEFAULT_SELL_SIGMA)
    }

    /// Converts `csr` to `SELL-C-σ` with explicit slice height `c` and sorting window `sigma`
    /// (`c >= 1`, `sigma >= c`, and `sigma % c == 0`).
    pub fn from_csr_with(csr: &Csr<T>, c: usize, sigma: usize) -> Result<Self, SellCSigmaError> {
        if c == 0 {
            return Err(SellCSigmaError::InvalidParams("slice height C must be >= 1"));
        }
        if sigma < c || !sigma.is_multiple_of(c) {
            return Err(SellCSigmaError::InvalidParams("sigma must be a positive multiple of C"));
        }
        let nrows = csr.nrows();
        let ncols = csr.ncols();
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err(SellCSigmaError::DimensionTooLarge);
        }

        let mut perm = vec![0u32; nrows];
        let mut inv_perm = vec![0u32; nrows];
        for win_start in (0..nrows).step_by(sigma) {
            let win_end = (win_start + sigma).min(nrows);
            let mut order: Vec<(usize, usize)> =
                (win_start..win_end).map(|r| (csr.row(r).0.len(), r)).collect();
            // Stable descending sort by row length; ties keep ascending row index.
            order.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
            for (offset, &(_, orig_r)) in order.iter().enumerate() {
                let sorted_r = win_start + offset;
                perm[sorted_r] = orig_r as u32;
                inv_perm[orig_r] = sorted_r as u32;
            }
        }

        let n_slices = nrows.div_ceil(c);
        let padded_rows = n_slices.checked_mul(c).ok_or(SellCSigmaError::DimensionTooLarge)?;
        let mut slice_ptr = Vec::with_capacity(n_slices + 1);
        slice_ptr.push(0usize);
        let mut row_lengths = vec![0u32; padded_rows];

        for s in 0..n_slices {
            let mut max_len = 0usize;
            for lane in 0..c {
                let sr = s * c + lane;
                if sr < nrows {
                    let rl = csr.row(perm[sr] as usize).0.len();
                    let rl_u32 = u32::try_from(rl).map_err(|_| SellCSigmaError::DimensionTooLarge)?;
                    row_lengths[sr] = rl_u32;
                    if rl > max_len {
                        max_len = rl;
                    }
                }
            }
            let added = max_len.checked_mul(c).ok_or(SellCSigmaError::DimensionTooLarge)?;
            let next = slice_ptr[s].checked_add(added).ok_or(SellCSigmaError::DimensionTooLarge)?;
            slice_ptr.push(next);
        }

        let total_slots = slice_ptr[n_slices];
        let mut col_idx = vec![0u32; total_slots];
        let mut values = vec![T::ZERO; total_slots];
        for s in 0..n_slices {
            let base = slice_ptr[s];
            for lane in 0..c {
                let sr = s * c + lane;
                if sr < nrows {
                    let orig_r = perm[sr] as usize;
                    let (cols, vals) = csr.row(orig_r);
                    for (k, (&col, &val)) in cols.iter().zip(vals).enumerate() {
                        let pos = base + k * c + lane;
                        col_idx[pos] = col;
                        values[pos] = val;
                    }
                }
            }
        }

        Ok(Self {
            nrows,
            ncols,
            nnz: csr.nnz(),
            c,
            sigma,
            slice_ptr,
            col_idx,
            values,
            row_lengths,
            perm,
            inv_perm,
        })
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    pub fn nnz(&self) -> usize {
        self.nnz
    }
    pub fn stored_entries(&self) -> usize {
        self.values.len()
    }
    pub fn c(&self) -> usize {
        self.c
    }
    pub fn sigma(&self) -> usize {
        self.sigma
    }
    pub fn slice_ptr(&self) -> &[usize] {
        &self.slice_ptr
    }
    pub fn col_idx(&self) -> &[u32] {
        &self.col_idx
    }
    pub fn values(&self) -> &[T] {
        &self.values
    }
    pub fn row_lengths(&self) -> &[u32] {
        &self.row_lengths
    }
    pub fn perm(&self) -> &[u32] {
        &self.perm
    }
    pub fn inv_perm(&self) -> &[u32] {
        &self.inv_perm
    }

    /// Ratio of stored slots (including slice padding) to true nonzeros (`>= 1.0`, or `1.0` if empty).
    pub fn fill_overhead(&self) -> f64 {
        if self.nnz == 0 { 1.0 } else { self.values.len() as f64 / self.nnz as f64 }
    }

    /// Deterministic parallel `y ← A x`, bit-for-bit identical to [`Csr::spmv`] across all thread counts.
    pub fn spmv(&self, x: &[T], y: &mut [T]) {
        if x.len() != self.ncols || y.len() != self.nrows || self.nrows == 0 {
            return;
        }
        // Each chunk spans an integer number of `sigma` windows (`perm` never crosses a `sigma` boundary),
        // so every parallel chunk writes strictly into its own `yc` sub-slice without `unsafe` or atomics.
        let win_chunk = CHUNK.div_ceil(self.sigma) * self.sigma;
        let slices_per_win = win_chunk / self.c;
        let par = par::worth_parallel(self.nnz);
        par::for_chunks_mut(y, win_chunk, par, |w_idx, yc| {
            let row_base = w_idx * win_chunk;
            let s_start = w_idx * slices_per_win;
            let s_end = (s_start + yc.len().div_ceil(self.c)).min(self.slice_ptr.len() - 1);
            for s in s_start..s_end {
                let base = self.slice_ptr[s];
                for lane in 0..self.c {
                    let sr = s * self.c + lane;
                    if sr < self.nrows {
                        let rlen = self.row_lengths[sr] as usize;
                        let mut acc = T::ZERO;
                        for k in 0..rlen {
                            let pos = base + k * self.c + lane;
                            acc += self.values[pos] * x[self.col_idx[pos] as usize];
                        }
                        let orig_r = self.perm[sr] as usize;
                        yc[orig_r - row_base] = acc;
                    }
                }
            }
        });
    }
}

impl<T: Scalar> LinearOperator<T> for SellCSigma<T> {
    fn dim(&self) -> usize {
        self.nrows
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y);
    }
}

impl<T: Scalar> DeviceLinearOperator<CpuVec<T>> for SellCSigma<T> {
    fn dim(&self) -> usize {
        self.nrows
    }
    fn apply(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) {
        self.spmv(&x.0, &mut y.0);
    }
}

// ---------------------------------------------------------------------------------------
// Csr64: 64-bit Sparse Row Matrix (`u64` indptr and column indices)
// ---------------------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Csr64Error {
    BadIndptr(&'static str),
    LengthMismatch,
    RowOutOfRange { row: u64, nrows: u64 },
    ColumnOutOfRange { row: u64, col: u64, ncols: u64 },
    AllocationOverflow,
    DimensionMismatch,
}

impl fmt::Display for Csr64Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BadIndptr(m) => write!(f, "invalid u64 indptr: {m}"),
            Self::LengthMismatch => write!(f, "indices/values length does not match u64 indptr"),
            Self::RowOutOfRange { row, nrows } => write!(f, "row {row} out of range (nrows = {nrows})"),
            Self::ColumnOutOfRange { row, col, ncols } => {
                write!(f, "column {col} out of range in row {row} (ncols = {ncols})")
            }
            Self::AllocationOverflow => write!(f, "matrix storage overflows host usize"),
            Self::DimensionMismatch => write!(f, "vector lengths do not match Csr64 dimensions"),
        }
    }
}

impl std::error::Error for Csr64Error {}

/// Compressed sparse row matrix with `u64` row pointers and `u64` column indices, lifting the
/// `u32::MAX` ceiling on both `ncols` and `nnz`.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr64<T: Scalar> {
    nrows: u64,
    ncols: u64,
    indptr: Vec<u64>,
    indices: Vec<u64>,
    values: Vec<T>,
}

impl<T: Scalar> Csr64<T> {
    pub fn try_new(
        nrows: u64,
        ncols: u64,
        indptr: Vec<u64>,
        indices: Vec<u64>,
        values: Vec<T>,
    ) -> Result<Self, Csr64Error> {
        let nrows_usize = usize::try_from(nrows).map_err(|_| Csr64Error::AllocationOverflow)?;
        let expected_indptr = nrows_usize.checked_add(1).ok_or(Csr64Error::AllocationOverflow)?;
        if indptr.len() != expected_indptr {
            return Err(Csr64Error::BadIndptr("length must be nrows + 1"));
        }
        if indptr[0] != 0 {
            return Err(Csr64Error::BadIndptr("must start at 0"));
        }
        if indptr.windows(2).any(|w| w[0] > w[1]) {
            return Err(Csr64Error::BadIndptr("must be non-decreasing"));
        }
        let nnz_u64 = indptr[nrows_usize];
        let nnz = usize::try_from(nnz_u64).map_err(|_| Csr64Error::AllocationOverflow)?;
        if indices.len() != nnz || values.len() != nnz {
            return Err(Csr64Error::LengthMismatch);
        }
        for row in 0..nrows_usize {
            let lo = indptr[row] as usize;
            let hi = indptr[row + 1] as usize;
            for &c in &indices[lo..hi] {
                if c >= ncols {
                    return Err(Csr64Error::ColumnOutOfRange { row: row as u64, col: c, ncols });
                }
            }
        }
        Ok(Self { nrows, ncols, indptr, indices, values })
    }

    /// Builds a `Csr64` from `(row, col, value)` 64-bit triplets; duplicates are summed and
    /// columns are sorted in ascending order within every row.
    ///
    /// Memory is `O(nrows + triplets.len())` regardless of `ncols`, supporting `ncols > u32::MAX`.
    pub fn from_triplets64(nrows: u64, ncols: u64, triplets: &[(u64, u64, T)]) -> Result<Self, Csr64Error> {
        let nrows_usize = usize::try_from(nrows).map_err(|_| Csr64Error::AllocationOverflow)?;
        let count_len = nrows_usize.checked_add(1).ok_or(Csr64Error::AllocationOverflow)?;
        let mut count = vec![0usize; count_len];
        for &(r, c, _) in triplets {
            if r >= nrows {
                return Err(Csr64Error::RowOutOfRange { row: r, nrows });
            }
            if c >= ncols {
                return Err(Csr64Error::ColumnOutOfRange { row: r, col: c, ncols });
            }
            count[r as usize + 1] += 1;
        }
        for i in 0..nrows_usize {
            count[i + 1] = count[i + 1].checked_add(count[i]).ok_or(Csr64Error::AllocationOverflow)?;
        }
        let mut next = count.clone();
        let mut entries = vec![(0u64, T::ZERO); triplets.len()];
        for &(r, c, v) in triplets {
            let slot = next[r as usize];
            entries[slot] = (c, v);
            next[r as usize] = slot + 1;
        }
        let mut indptr = Vec::with_capacity(count_len);
        let mut indices = Vec::with_capacity(triplets.len());
        let mut values = Vec::with_capacity(triplets.len());
        indptr.push(0u64);
        for r in 0..nrows_usize {
            let row = &mut entries[count[r]..count[r + 1]];
            row.sort_by_key(|&(c, _)| c);
            let mut last: Option<u64> = None;
            for &(c, v) in row.iter() {
                if last == Some(c) {
                    if let Some(prev) = values.last_mut() {
                        *prev += v;
                    }
                } else {
                    indices.push(c);
                    values.push(v);
                    last = Some(c);
                }
            }
            indptr.push(indices.len() as u64);
        }
        Self::try_new(nrows, ncols, indptr, indices, values)
    }

    /// Widens a 32-bit [`Csr<T>`] into `Csr64<T>`.
    pub fn from_csr(csr: &Csr<T>) -> Self {
        Self {
            nrows: csr.nrows() as u64,
            ncols: csr.ncols() as u64,
            indptr: csr.indptr().iter().map(|&p| p as u64).collect(),
            indices: csr.indices().iter().map(|&c| c as u64).collect(),
            values: csr.values().to_vec(),
        }
    }

    /// Checked narrowing conversion to 32-bit [`Csr<T>`]; returns [`CsrError::DimensionTooLarge`]
    /// if `nrows`, `ncols`, or `nnz` exceeds `u32::MAX`.
    pub fn try_to_csr(&self) -> Result<Csr<T>, CsrError> {
        if self.exceeds_u32_limits() {
            return Err(CsrError::DimensionTooLarge);
        }
        let nrows = usize::try_from(self.nrows).map_err(|_| CsrError::DimensionTooLarge)?;
        let ncols = usize::try_from(self.ncols).map_err(|_| CsrError::DimensionTooLarge)?;
        let mut indptr = Vec::with_capacity(self.indptr.len());
        for &p in &self.indptr {
            indptr.push(usize::try_from(p).map_err(|_| CsrError::DimensionTooLarge)?);
        }
        let mut indices = Vec::with_capacity(self.indices.len());
        for &c in &self.indices {
            indices.push(u32::try_from(c).map_err(|_| CsrError::DimensionTooLarge)?);
        }
        Csr::try_new(nrows, ncols, indptr, indices, self.values.clone())
    }

    /// Returns `true` if `nrows > u32::MAX`, `ncols > u32::MAX`, or `nnz > u32::MAX`.
    #[inline]
    pub fn exceeds_u32_limits(&self) -> bool {
        let cap = u32::MAX as u64;
        self.nrows > cap || self.ncols > cap || self.nnz() > cap
    }

    pub fn nrows(&self) -> u64 {
        self.nrows
    }
    pub fn ncols(&self) -> u64 {
        self.ncols
    }
    pub fn nnz(&self) -> u64 {
        self.values.len() as u64
    }
    pub fn indptr(&self) -> &[u64] {
        &self.indptr
    }
    pub fn indices(&self) -> &[u64] {
        &self.indices
    }
    pub fn values(&self) -> &[T] {
        &self.values
    }

    pub fn row(&self, i: usize) -> (&[u64], &[T]) {
        let lo = self.indptr[i] as usize;
        let hi = self.indptr[i + 1] as usize;
        (&self.indices[lo..hi], &self.values[lo..hi])
    }

    pub fn map<U: Scalar>(&self, f: impl Fn(T) -> U) -> Csr64<U> {
        Csr64 {
            nrows: self.nrows,
            ncols: self.ncols,
            indptr: self.indptr.clone(),
            indices: self.indices.clone(),
            values: self.values.iter().map(|&v| f(v)).collect(),
        }
    }

    pub fn to_complex(&self) -> Csr64<C64> {
        self.map(|v| v.to_c64())
    }

    /// Checked deterministic parallel `y ← A x`.
    pub fn spmv_checked(&self, x: &[T], y: &mut [T]) -> Result<(), Csr64Error> {
        if x.len() as u64 != self.ncols || y.len() as u64 != self.nrows {
            return Err(Csr64Error::DimensionMismatch);
        }
        let par = par::worth_parallel(self.values.len());
        par::for_chunks_mut(y, CHUNK, par, |c, yc| {
            let start = c * CHUNK;
            let ip = &self.indptr[start..=start + yc.len()];
            for (yi, w) in yc.iter_mut().zip(ip.windows(2)) {
                let (lo, hi) = (w[0] as usize, w[1] as usize);
                let mut s = T::ZERO;
                for (&col, &val) in self.indices[lo..hi].iter().zip(&self.values[lo..hi]) {
                    s += val * x[col as usize];
                }
                *yi = s;
            }
        });
        Ok(())
    }

    /// Deterministic parallel `y ← A x` (no-op if slice lengths do not match `(ncols, nrows)`).
    pub fn spmv(&self, x: &[T], y: &mut [T]) {
        let _ = self.spmv_checked(x, y);
    }
}

impl<T: Scalar> From<&Csr<T>> for Csr64<T> {
    fn from(csr: &Csr<T>) -> Self {
        Self::from_csr(csr)
    }
}

impl<T: Scalar> TryFrom<&Csr64<T>> for Csr<T> {
    type Error = CsrError;
    fn try_from(csr64: &Csr64<T>) -> Result<Self, CsrError> {
        csr64.try_to_csr()
    }
}

impl<T: Scalar> LinearOperator<T> for Csr64<T> {
    fn dim(&self) -> usize {
        self.nrows as usize
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y);
    }
}

impl<T: Scalar> DeviceLinearOperator<CpuVec<T>> for Csr64<T> {
    fn dim(&self) -> usize {
        self.nrows as usize
    }
    fn apply(&self, x: &CpuVec<T>, y: &mut CpuVec<T>) {
        self.spmv(&x.0, &mut y.0);
    }
}
