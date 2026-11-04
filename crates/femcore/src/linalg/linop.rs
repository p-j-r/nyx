//! Operator and preconditioner interfaces.
//!
//! The traits carry the *fusion hooks* that made the `rsolvers` kernels fast (SpMV + dot,
//! preconditioner + two dots, element-wise preconditioner scaling). They have defaults, so a
//! plain `apply` is enough to implement them, but a fused implementation can be added later
//! without changing any solver. Hooks must give the same bits as the unfused default (tested
//! for [`Csr`](crate::linalg::csr::Csr)).
//!
//! Block (multi-right-hand-side) application is a first-class method because both the block
//! Krylov–Schur and `faer`'s sparse triangular solves are much cheaper per vector when several
//! vectors are processed together.

use crate::linalg::scalar::Scalar;
use crate::linalg::vecops;

/// A square linear operator `y = A x`.
pub trait LinearOperator<T: Scalar> {
    /// Dimension `n` of the (square) operator.
    fn dim(&self) -> usize;

    /// `y ← A x`. Both slices have length [`dim`](Self::dim).
    fn apply(&self, x: &[T], y: &mut [T]);

    /// `y[:, c] ← A x[:, c]` for `c < nrhs`; columns are contiguous, each of length `dim()`.
    /// Override when the operator can process several vectors faster than one at a time.
    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        let n = self.dim();
        for c in 0..nrhs {
            self.apply(&x[c * n..(c + 1) * n], &mut y[c * n..(c + 1) * n]);
        }
    }

    /// `y ← A x`, returning `xᴴ y`. Implementations may fuse the reduction into the product.
    fn apply_dotc(&self, x: &[T], y: &mut [T]) -> T {
        self.apply(x, y);
        vecops::dotc(x, y)
    }

    /// `y ← A x`, returning the unconjugated `xᵀ y` (COCG needs the bilinear form).
    fn apply_dotu(&self, x: &[T], y: &mut [T]) -> T {
        self.apply(x, y);
        vecops::dotu(x, y)
    }

    /// Block form of `apply` when the caller already holds `bx = B x` for the weight matrix
    /// `B` the operator was built from (generalised problems: `A_op = (A−σB)⁻¹ B`). The default
    /// ignores `bx`; shift-and-invert operators use it to skip one product with `B` per vector.
    fn apply_block_with_bx(&self, x: &[T], _bx: &[T], y: &mut [T], nrhs: usize) {
        self.apply_block(x, y, nrhs)
    }
}

impl<T: Scalar, O: LinearOperator<T> + ?Sized> LinearOperator<T> for &O {
    fn dim(&self) -> usize {
        (**self).dim()
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        (**self).apply(x, y)
    }
    fn apply_block(&self, x: &[T], y: &mut [T], nrhs: usize) {
        (**self).apply_block(x, y, nrhs)
    }
    fn apply_dotc(&self, x: &[T], y: &mut [T]) -> T {
        (**self).apply_dotc(x, y)
    }
    fn apply_dotu(&self, x: &[T], y: &mut [T]) -> T {
        (**self).apply_dotu(x, y)
    }
    fn apply_block_with_bx(&self, x: &[T], bx: &[T], y: &mut [T], nrhs: usize) {
        (**self).apply_block_with_bx(x, bx, y, nrhs)
    }
}

/// How a preconditioner acts element-wise, so solvers can fuse it into their vector updates.
#[derive(Debug, Clone, Copy)]
pub enum Scaling<'a, T> {
    /// `z = r`.
    Identity,
    /// `z[i] = d[i]·r[i]`.
    Diagonal(&'a [T]),
}

/// A preconditioner application `z ← M⁻¹ r`.
pub trait Preconditioner<T: Scalar> {
    fn apply(&self, r: &[T], z: &mut [T]);

    /// `z ← M⁻¹ r`, returning `(zᴴ r, zᴴ z)`; fused implementations must match the default
    /// bit for bit.
    fn apply_dots(&self, r: &[T], z: &mut [T]) -> (T, T) {
        self.apply(r, z);
        (vecops::dotc(z, r), vecops::dotc(z, z))
    }

    /// Element-wise structure, if any (see [`Scaling`]).
    fn scaling(&self) -> Option<Scaling<'_, T>> {
        None
    }
}

/// `M = I`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Identity;

impl<T: Scalar> Preconditioner<T> for Identity {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
    }
    fn scaling(&self) -> Option<Scaling<'_, T>> {
        Some(Scaling::Identity)
    }
}

/// Wraps a closure `(x, y)` as a [`LinearOperator`] (the analogue of PETSc `MatShell`).
pub struct FnOperator<F> {
    n: usize,
    f: F,
}

impl<F> FnOperator<F> {
    pub fn new(n: usize, f: F) -> Self {
        Self { n, f }
    }
}

impl<T: Scalar, F: Fn(&[T], &mut [T])> LinearOperator<T> for FnOperator<F> {
    fn dim(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        (self.f)(x, y)
    }
}
