//! Preconditioners (PETSc `PC`). More (ILU/IC, Schwarz, AMG) follow in later items.

use crate::solvers::csr::Csr;
use crate::solvers::linop::{Preconditioner, Scaling};
use crate::solvers::scalar::Scalar;

/// Point Jacobi: `M⁻¹ = diag(A)⁻¹` (PETSc `PCJACOBI`). A zero diagonal entry is replaced by 1
/// as PETSc does. For complex symmetric matrices the plain (unconjugated) inverse is used.
#[derive(Debug, Clone)]
pub struct Jacobi<T: Scalar> {
    inv: Vec<T>,
}

impl<T: Scalar> Jacobi<T> {
    pub fn new(a: &Csr<T>) -> Self {
        Self::from_diagonal(a.diagonal())
    }

    /// From the diagonal alone, for operators that are never assembled (see [`crate::solvers::matfree`]).
    pub fn from_diagonal(diag: Vec<T>) -> Self {
        Self { inv: diag.into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect() }
    }
}

impl<T: Scalar> Preconditioner<T> for Jacobi<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        for ((zi, di), ri) in z.iter_mut().zip(&self.inv).zip(r) {
            *zi = *di * *ri;
        }
    }
    fn scaling(&self) -> Option<Scaling<'_, T>> {
        Some(Scaling::Diagonal(&self.inv))
    }
}
