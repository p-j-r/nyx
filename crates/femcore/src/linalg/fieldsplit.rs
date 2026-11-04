//! Field-split and Schur-complement preconditioners (port item P5; the role of PETSc `PCFIELDSPLIT`) for
//! coupled problems: saddle points, piezoelectric `(u, φ)`, thermo-mechanical or opto-mechanical couplings.
//!
//! * [`FieldSplit`]: any number of fields (index sets); one sub-preconditioner per diagonal block
//!   `A_ii` (built by a user closure: an exact [`DirectSolve`], AMG, ILU, ...); combined as block
//!   Jacobi ([`SplitKind::Additive`]), block Gauss–Seidel ([`SplitKind::Multiplicative`]) or its
//!   symmetric form ([`SplitKind::SymmetricMultiplicative`]).
//! * [`SchurFieldSplit`]: two fields `[[A00, A01], [A10, A11]]` with the Schur complement
//!   `S = A11 − A10 A00⁻¹ A01` approximated by `A11`, by the "selfp" matrix `A11 − A10 diag(A00)⁻¹ A01`
//!   (PETSc's `-pc_fieldsplit_schur_precondition selfp`) or by a user matrix, and applied as the
//!   block diagonal, lower, upper or full factorisation. With exact `A00⁻¹` and `S⁻¹` the preconditioned
//!   operator has 3 (diagonal), 2 (triangular) or 1 (full) distinct eigenvalues, so GMRES converges in
//!   that many iterations: the tests use this as an exactness check.
//!
//! All solves are sequential; `SchurFact::Diag` applies `diag(A00⁻¹, S⁻¹)` (set `negate_schur` for the
//! `diag(A00⁻¹, −S⁻¹)` that PETSc uses and that is positive definite for saddle points with `A11 = 0`,
//! as MINRES needs).

use crate::linalg::csr::Csr;
use crate::linalg::direct::{DirectSolve, SolverError};
use crate::linalg::linop::Preconditioner;
use crate::linalg::scalar::Scalar;
use std::cell::RefCell;

/// Builds a sub-preconditioner from a diagonal block.
pub type MakePc<'a, T> = &'a dyn Fn(&Csr<T>) -> Result<Box<dyn Preconditioner<T>>, SolverError>;

/// Extracts `A[rows, cols]` (rows/cols are global indices; the result is renumbered locally).
pub fn extract_block<T: Scalar>(a: &Csr<T>, rows: &[usize], cols: &[usize]) -> Csr<T> {
    let mut local = vec![usize::MAX; a.ncols()];
    for (l, &g) in cols.iter().enumerate() {
        local[g] = l;
    }
    let mut trip = Vec::new();
    for (li, &g) in rows.iter().enumerate() {
        let (cs, vs) = a.row(g);
        for (&c, &v) in cs.iter().zip(vs) {
            let lc = local[c as usize];
            if lc != usize::MAX {
                trip.push((li, lc, v));
            }
        }
    }
    Csr::from_triplets(rows.len(), cols.len(), &trip).expect("block extraction")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitKind {
    Additive,
    Multiplicative,
    SymmetricMultiplicative,
}

struct Field<T: Scalar> {
    idx: Vec<usize>,
    pc: Box<dyn Preconditioner<T>>,
}

pub struct FieldSplit<T: Scalar> {
    n: usize,
    a: Csr<T>,
    fields: Vec<Field<T>>,
    kind: SplitKind,
    scratch: RefCell<Vec<T>>,
}

impl<T: Scalar> FieldSplit<T> {
    pub fn new(a: &Csr<T>, fields: Vec<Vec<usize>>, make: MakePc<'_, T>, kind: SplitKind) -> Result<Self, SolverError> {
        let n = a.nrows();
        let mut seen = vec![false; n];
        for f in &fields {
            for &i in f {
                if i >= n || seen[i] {
                    return Err(SolverError::InvalidArgument(
                        "fields must be disjoint index sets within the matrix dimension",
                    ));
                }
                seen[i] = true;
            }
        }
        if seen.iter().any(|s| !s) {
            return Err(SolverError::InvalidArgument("fields must cover all unknowns"));
        }
        let fields = fields
            .into_iter()
            .map(|idx| {
                let blk = extract_block(a, &idx, &idx);
                Ok(Field { pc: make(&blk)?, idx })
            })
            .collect::<Result<Vec<_>, SolverError>>()?;
        Ok(Self { n, a: a.clone(), fields, kind, scratch: RefCell::new(vec![T::ZERO; n]) })
    }

    fn sweep(&self, order: &[usize], r: &[T], z: &mut [T]) {
        // block Gauss–Seidel in the given order, maintaining the residual r − A z
        let mut res = r.to_vec();
        let mut ax = self.scratch.borrow_mut();
        for &fi in order {
            let f = &self.fields[fi];
            let rf: Vec<T> = f.idx.iter().map(|&i| res[i]).collect();
            let mut dz = vec![T::ZERO; f.idx.len()];
            f.pc.apply(&rf, &mut dz);
            // scatter the correction, update the residual with the corresponding columns of A
            let mut full = vec![T::ZERO; self.n];
            for (&i, &d) in f.idx.iter().zip(&dz) {
                z[i] += d;
                full[i] = d;
            }
            self.a.spmv(&full, &mut ax);
            for (ri, ai) in res.iter_mut().zip(ax.iter()) {
                *ri -= *ai;
            }
        }
    }
}

impl<T: Scalar> Preconditioner<T> for FieldSplit<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.fill(T::ZERO);
        match self.kind {
            SplitKind::Additive => {
                for f in &self.fields {
                    let rf: Vec<T> = f.idx.iter().map(|&i| r[i]).collect();
                    let mut dz = vec![T::ZERO; f.idx.len()];
                    f.pc.apply(&rf, &mut dz);
                    for (&i, &d) in f.idx.iter().zip(&dz) {
                        z[i] += d;
                    }
                }
            }
            SplitKind::Multiplicative => {
                let order: Vec<usize> = (0..self.fields.len()).collect();
                self.sweep(&order, r, z);
            }
            SplitKind::SymmetricMultiplicative => {
                let fwd: Vec<usize> = (0..self.fields.len()).collect();
                let bwd: Vec<usize> = (0..self.fields.len()).rev().collect();
                self.sweep(&fwd, r, z);
                // second sweep on the updated residual
                let mut res = r.to_vec();
                let mut az = vec![T::ZERO; self.n];
                self.a.spmv(z, &mut az);
                for (ri, ai) in res.iter_mut().zip(&az) {
                    *ri -= *ai;
                }
                let mut z2 = vec![T::ZERO; self.n];
                self.sweep(&bwd, &res, &mut z2);
                for (zi, di) in z.iter_mut().zip(&z2) {
                    *zi += *di;
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchurFact {
    Diag,
    Lower,
    Upper,
    Full,
}

pub enum SchurPre<T: Scalar> {
    /// `S ≈ A11`.
    A11,
    /// `S ≈ A11 − A10 diag(A00)⁻¹ A01`.
    SelfP,
    /// A user-supplied approximation (e.g. the exact Schur complement in tests, a mass matrix for Stokes).
    User(Csr<T>),
}

pub struct SchurFieldSplit<T: Scalar> {
    i0: Vec<usize>,
    i1: Vec<usize>,
    a00: Box<dyn Preconditioner<T>>,
    s: Box<dyn Preconditioner<T>>,
    a01: Csr<T>,
    a10: Csr<T>,
    fact: SchurFact,
    negate: bool,
}

impl<T: Scalar> SchurFieldSplit<T> {
    /// `i0`, `i1`: the two fields (disjoint, covering all unknowns). `make_a00` / `make_s` build the solvers for
    /// `A00` and for the Schur approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        a: &Csr<T>,
        i0: Vec<usize>,
        i1: Vec<usize>,
        fact: SchurFact,
        pre: SchurPre<T>,
        negate_schur: bool,
        make_a00: MakePc<'_, T>,
        make_s: MakePc<'_, T>,
    ) -> Result<Self, SolverError> {
        if i0.len() + i1.len() != a.nrows() {
            return Err(SolverError::InvalidArgument("the two fields must cover all unknowns"));
        }
        let a00 = extract_block(a, &i0, &i0);
        let a01 = extract_block(a, &i0, &i1);
        let a10 = extract_block(a, &i1, &i0);
        let a11 = extract_block(a, &i1, &i1);
        let s_mat = match pre {
            SchurPre::A11 => a11,
            SchurPre::SelfP => {
                let dinv: Vec<(usize, usize, T)> = a00
                    .diagonal()
                    .iter()
                    .enumerate()
                    .map(|(i, &d)| (i, i, if d == T::ZERO { T::ONE } else { T::ONE / d }))
                    .collect();
                let dm = Csr::from_triplets(a00.nrows(), a00.nrows(), &dinv).unwrap();
                let prod = a10.matmul(&dm.matmul(&a01));
                Csr::axpby(T::ONE, &a11, -T::ONE, &prod)
            }
            SchurPre::User(m) => m,
        };
        Ok(Self { a00: make_a00(&a00)?, s: make_s(&s_mat)?, a01, a10, i0, i1, fact, negate: negate_schur })
    }
}

impl<T: Scalar> Preconditioner<T> for SchurFieldSplit<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        let (n0, n1) = (self.i0.len(), self.i1.len());
        let r0: Vec<T> = self.i0.iter().map(|&i| r[i]).collect();
        let r1: Vec<T> = self.i1.iter().map(|&i| r[i]).collect();
        let (mut z0, mut z1) = (vec![T::ZERO; n0], vec![T::ZERO; n1]);
        let sign = if self.negate { -T::ONE } else { T::ONE };
        let mut tmp0 = vec![T::ZERO; n0];
        let mut tmp1 = vec![T::ZERO; n1];
        match self.fact {
            SchurFact::Diag => {
                self.a00.apply(&r0, &mut z0);
                self.s.apply(&r1, &mut z1);
                z1.iter_mut().for_each(|v| *v *= sign);
            }
            SchurFact::Lower => {
                self.a00.apply(&r0, &mut z0);
                self.a10.spmv(&z0, &mut tmp1);
                let rhs: Vec<T> = r1.iter().zip(&tmp1).map(|(a, b)| *a - *b).collect();
                self.s.apply(&rhs, &mut z1);
            }
            SchurFact::Upper => {
                self.s.apply(&r1, &mut z1);
                self.a01.spmv(&z1, &mut tmp0);
                let rhs: Vec<T> = r0.iter().zip(&tmp0).map(|(a, b)| *a - *b).collect();
                self.a00.apply(&rhs, &mut z0);
            }
            SchurFact::Full => {
                let mut t0 = vec![T::ZERO; n0];
                self.a00.apply(&r0, &mut t0);
                self.a10.spmv(&t0, &mut tmp1);
                let rhs1: Vec<T> = r1.iter().zip(&tmp1).map(|(a, b)| *a - *b).collect();
                self.s.apply(&rhs1, &mut z1);
                self.a01.spmv(&z1, &mut tmp0);
                let mut corr = vec![T::ZERO; n0];
                self.a00.apply(&tmp0, &mut corr);
                for i in 0..n0 {
                    z0[i] = t0[i] - corr[i];
                }
            }
        }
        for (&i, &v) in self.i0.iter().zip(&z0) {
            z[i] = v;
        }
        for (&i, &v) in self.i1.iter().zip(&z1) {
            z[i] = v;
        }
    }
}

/// Convenience builder for the common case: an exact (sparse direct) solver for a block.
pub fn direct_maker<T: crate::linalg::faer_glue::FaerScalar>(
    definiteness: crate::linalg::direct::Definiteness,
) -> impl Fn(&Csr<T>) -> Result<Box<dyn Preconditioner<T>>, SolverError> {
    move |a| Ok(Box::new(DirectSolve::new(a, definiteness)?) as Box<dyn Preconditioner<T>>)
}
