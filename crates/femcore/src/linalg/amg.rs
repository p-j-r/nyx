//! Smoothed-aggregation algebraic multigrid for real and **complex symmetric** matrices
//! (port items "AMG for complex symmetric matrices" and the multigrid half of the shifted-Laplacian
//! preconditioner, P2).
//!
//! * Strength of connection uses moduli relative to the largest off-diagonal of the row:
//!   `|a_ij| ≥ θ·max_k |a_ik|`.
//! * Aggregation: the standard three-pass greedy scheme; one real near-null-space vector
//!   (the constant) is carried down the hierarchy; the tentative prolongator is normalised per
//!   aggregate and then smoothed by one damped Jacobi step with the *complex* diagonal.
//! * Coarse operators are Galerkin products `R A P` with `R = Pᵀ` (plain transpose), which keeps a
//!   complex symmetric matrix complex symmetric. Set [`AmgOptions::hermitian`] to use `R = Pᴴ`.
//! * The smoother is symmetric Gauss–Seidel applied as *forward before, backward after* the coarse
//!   correction, so the V-cycle is a symmetric (`Mᵀ = M`) preconditioner and can drive COCG.
//! * The coarsest level is solved directly with a sparse LU.
//!
//! This is not the nyx AMG (rigid-body near-null-space, 3×3 node blocks, Chebyshev smoother, f32
//! cycle): those are real/elasticity specific and remain in nyx; this one covers the scalar
//! complex Helmholtz/PML family. Gauss–Seidel sweeps are sequential here (deterministic);
//! a coloured or Chebyshev smoother is the route to threads (see `docs/STATUS.md`).

use crate::linalg::csr::Csr;
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::faer_glue::FaerScalar;
use crate::linalg::linop::Preconditioner;
use crate::linalg::scalar::Scalar;
use std::cell::RefCell;

#[derive(Debug, Clone)]
pub struct AmgOptions {
    /// Strength threshold `θ` (relative to the row's largest off-diagonal modulus).
    pub theta: f64,
    /// Stop coarsening below this many unknowns.
    pub max_coarse: usize,
    pub max_levels: usize,
    /// Smooth the tentative prolongator (smoothed aggregation); `false` = plain aggregation.
    pub smooth_prolongator: bool,
    /// Gauss–Seidel sweeps before and after the coarse correction.
    pub sweeps: usize,
    /// Use `R = Pᴴ` (Hermitian problems) instead of `R = Pᵀ` (complex symmetric).
    pub hermitian: bool,
}

impl Default for AmgOptions {
    fn default() -> Self {
        Self { theta: 0.25, max_coarse: 150, max_levels: 20, smooth_prolongator: true, sweeps: 1, hermitian: false }
    }
}

struct Level<T: Scalar> {
    a: Csr<T>,
    inv_diag: Vec<T>,
    /// Prolongation `P` (fine × coarse) and restriction `R` (coarse × fine); empty on the last level.
    p: Option<Csr<T>>,
    r: Option<Csr<T>>,
}

struct Scratch<T> {
    rhs: Vec<T>,
    x: Vec<T>,
    res: Vec<T>,
    tmp: Vec<T>,
}

/// Summary of the hierarchy, for logs and tests.
#[derive(Debug, Clone)]
pub struct AmgStats {
    pub level_sizes: Vec<usize>,
    pub level_nnz: Vec<usize>,
    /// `Σ nnz(A_l) / nnz(A_0)`.
    pub operator_complexity: f64,
}

pub struct Amg<T: FaerScalar> {
    levels: Vec<Level<T>>,
    coarse: ShiftedFactorization<T>,
    scratch: RefCell<Vec<Scratch<T>>>,
    opts: AmgOptions,
}

impl<T: FaerScalar> Amg<T> {
    pub fn new(a: &Csr<T>, opts: AmgOptions) -> Result<Self, crate::linalg::direct::SolverError> {
        assert_eq!(a.nrows(), a.ncols());
        let mut levels: Vec<Level<T>> = Vec::new();
        let mut cur = a.clone();
        let mut nullvec = vec![1.0f64; a.nrows()];
        loop {
            let n = cur.nrows();
            let inv_diag: Vec<T> =
                cur.diagonal().into_iter().map(|d| if d == T::ZERO { T::ONE } else { T::ONE / d }).collect();
            if n <= opts.max_coarse || levels.len() + 1 >= opts.max_levels {
                levels.push(Level { a: cur, inv_diag, p: None, r: None });
                break;
            }
            let agg = aggregate(&cur, opts.theta);
            let nagg = agg.iter().copied().max().map_or(0, |m| m + 1);
            if nagg == 0 || nagg >= n {
                // No coarsening possible (e.g. diagonal matrix): finish with a direct solve here.
                levels.push(Level { a: cur, inv_diag, p: None, r: None });
                break;
            }
            let (p_tent, coarse_null) = tentative_prolongator(&agg, nagg, &nullvec);
            let p = if opts.smooth_prolongator { smooth(&cur, &inv_diag, &p_tent) } else { p_tent };
            let r = {
                let t = p.transpose();
                if opts.hermitian { t.map(|v| v.conj()) } else { t }
            };
            let coarse = r.matmul(&cur.matmul(&p));
            levels.push(Level { a: cur, inv_diag, p: Some(p), r: Some(r) });
            cur = coarse;
            nullvec = coarse_null;
        }
        let last = &levels.last().unwrap().a;
        let mut coarse = ShiftedFactorization::new(last, None, Definiteness::General)?;
        coarse.set_shift(T::ZERO)?;
        let scratch = levels
            .iter()
            .map(|l| {
                let n = l.a.nrows();
                Scratch { rhs: vec![T::ZERO; n], x: vec![T::ZERO; n], res: vec![T::ZERO; n], tmp: vec![T::ZERO; n] }
            })
            .collect();
        Ok(Self { levels, coarse, scratch: RefCell::new(scratch), opts })
    }

    pub fn stats(&self) -> AmgStats {
        let level_sizes: Vec<usize> = self.levels.iter().map(|l| l.a.nrows()).collect();
        let level_nnz: Vec<usize> = self.levels.iter().map(|l| l.a.nnz()).collect();
        let operator_complexity = level_nnz.iter().sum::<usize>() as f64 / level_nnz[0] as f64;
        AmgStats { level_sizes, level_nnz, operator_complexity }
    }

    fn cycle(&self, lvl: usize, sc: &mut [Scratch<T>]) {
        let (cur, rest) = sc.split_first_mut().unwrap();
        let level = &self.levels[lvl];
        if level.p.is_none() {
            cur.x.copy_from_slice(&cur.rhs);
            self.coarse
                .solve_block_in_place(&mut cur.x, 1)
                .expect("AMG coarse solver is always factored in construction");
            return;
        }
        cur.x.fill(T::ZERO);
        for _ in 0..self.opts.sweeps {
            gauss_seidel(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, true);
        }
        // res = rhs − A x
        level.a.spmv(&cur.x, &mut cur.res);
        for (r, b) in cur.res.iter_mut().zip(&cur.rhs) {
            *r = *b - *r;
        }
        level.r.as_ref().unwrap().spmv(&cur.res, &mut rest[0].rhs);
        self.cycle(lvl + 1, rest);
        level.p.as_ref().unwrap().spmv(&rest[0].x, &mut cur.tmp);
        for (x, t) in cur.x.iter_mut().zip(&cur.tmp) {
            *x += *t;
        }
        for _ in 0..self.opts.sweeps {
            gauss_seidel(&level.a, &level.inv_diag, &cur.rhs, &mut cur.x, false);
        }
    }
}

impl<T: FaerScalar> Preconditioner<T> for Amg<T> {
    /// One V-cycle from a zero initial guess: `z ≈ A⁻¹ r`.
    fn apply(&self, r: &[T], z: &mut [T]) {
        let mut sc = self.scratch.borrow_mut();
        sc[0].rhs.copy_from_slice(r);
        self.cycle(0, &mut sc);
        z.copy_from_slice(&sc[0].x);
    }
}

/// One Gauss–Seidel sweep (`forward` or backward) for `A x = b`, in place on `x`.
pub(crate) fn gauss_seidel<T: Scalar>(a: &Csr<T>, inv_diag: &[T], b: &[T], x: &mut [T], forward: bool) {
    let n = a.nrows();
    let mut sweep = |i: usize| {
        let (cols, vals) = a.row(i);
        let mut s = b[i];
        for (&j, &v) in cols.iter().zip(vals) {
            s -= v * x[j as usize];
        }
        // s = b_i − Σ_j a_ij x_j (diagonal included); add back the diagonal contribution.
        x[i] += s * inv_diag[i];
    };
    if forward {
        (0..n).for_each(&mut sweep);
    } else {
        (0..n).rev().for_each(&mut sweep);
    }
}

/// Standard three-pass aggregation on the strength graph; returns the aggregate of each node.
fn aggregate<T: Scalar>(a: &Csr<T>, theta: f64) -> Vec<usize> {
    let n = a.nrows();
    // Classical strength: |a_ij| >= theta * max_{k != i} |a_ik| (robust to diagonal shifts, which
    // matters for shifted-Laplacian operators whose diagonal carries the complex shift).
    let rowmax: Vec<f64> = (0..n)
        .map(|i| {
            let (cols, vals) = a.row(i);
            cols.iter().zip(vals).filter(|&(&j, _)| j as usize != i).map(|(_, v)| v.modulus()).fold(0.0, f64::max)
        })
        .collect();
    let strong = |i: usize| -> Vec<(usize, f64)> {
        let (cols, vals) = a.row(i);
        cols.iter()
            .zip(vals)
            .filter(|&(&j, v)| j as usize != i && v.modulus() > 0.0 && v.modulus() >= theta * rowmax[i])
            .map(|(&j, v)| (j as usize, v.modulus()))
            .collect()
    };
    let nbrs: Vec<Vec<(usize, f64)>> = (0..n).map(strong).collect();
    const NONE: usize = usize::MAX;
    let mut agg = vec![NONE; n];
    let mut next = 0usize;
    // Distance-2 maximal independent set of roots: no two roots within two strong-graph hops.
    // (Taking a root with only its immediate neighbours gives aggregates of ~5 nodes in 2-D, ~7
    // in 3-D and a coarse-operator complexity that is too high; MIS-2 gives ~3^d.)
    let mut blocked = vec![false; n];
    for i in 0..n {
        if !blocked[i] && !nbrs[i].is_empty() {
            blocked[i] = true;
            agg[i] = next;
            for &(j, _) in &nbrs[i] {
                blocked[j] = true;
                if agg[j] == NONE {
                    agg[j] = next;
                }
                for &(k, _) in &nbrs[j] {
                    blocked[k] = true;
                }
            }
            next += 1;
        }
    }
    // Nodes two hops from a root: join the strongest already-assigned neighbour's aggregate
    // (repeat until nothing changes so that every node connected to a root is reached).
    loop {
        let snapshot = agg.clone();
        let mut changed = false;
        for i in 0..n {
            if agg[i] != NONE {
                continue;
            }
            if let Some(&(j, _)) =
                nbrs[i].iter().filter(|&&(j, _)| snapshot[j] != NONE).max_by(|x, y| x.1.total_cmp(&y.1))
            {
                agg[i] = snapshot[j];
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    // Isolated nodes (no strong connection) become singleton aggregates.
    for a in agg.iter_mut() {
        if *a == NONE {
            *a = next;
            next += 1;
        }
    }
    agg
}

/// `P_tent` with columns normalised per aggregate for the near-null-space vector `b`; also
/// returns the coarse near-null-space vector (the aggregate norms).
fn tentative_prolongator<T: Scalar>(agg: &[usize], nagg: usize, b: &[f64]) -> (Csr<T>, Vec<f64>) {
    let mut nrm = vec![0.0f64; nagg];
    for (i, &g) in agg.iter().enumerate() {
        nrm[g] += b[i] * b[i];
    }
    for v in &mut nrm {
        *v = v.sqrt();
    }
    let n = agg.len();
    let indptr: Vec<usize> = (0..=n).collect();
    let indices: Vec<u32> = agg.iter().map(|&g| g as u32).collect();
    let values: Vec<T> =
        (0..n).map(|i| T::from_f64(if nrm[agg[i]] > 0.0 { b[i] / nrm[agg[i]] } else { 0.0 })).collect();
    (Csr::try_new(n, nagg, indptr, indices, values).unwrap(), nrm)
}

/// `P = (I − ω D⁻¹ A) P_tent` with `ω = 4/(3ρ)` and `ρ` the Gershgorin bound of `D⁻¹A` in moduli.
fn smooth<T: Scalar>(a: &Csr<T>, inv_diag: &[T], p_tent: &Csr<T>) -> Csr<T> {
    let n = a.nrows();
    let rho = (0..n)
        .map(|i| {
            let (_, v) = a.row(i);
            v.iter().map(|x| x.modulus()).sum::<f64>() * inv_diag[i].modulus()
        })
        .fold(0.0, f64::max)
        .max(f64::MIN_POSITIVE);
    let omega = 4.0 / (3.0 * rho);
    let ap = a.matmul(p_tent);
    // Scale rows of AP by −ω D⁻¹ and add P_tent.
    let scaled = Csr::try_new(
        ap.nrows(),
        ap.ncols(),
        ap.indptr().to_vec(),
        ap.indices().to_vec(),
        (0..n)
            .flat_map(|i| {
                let (_, v) = ap.row(i);
                v.iter().map(move |x| *x * inv_diag[i].scale(-omega)).collect::<Vec<_>>()
            })
            .collect(),
    )
    .unwrap();
    Csr::axpby(T::ONE, p_tent, T::ONE, &scaled)
}

/// Builds the shifted-Laplacian operator `K − (β₁ + iβ₂)·k²·M` (Erlangga et al.): the matrix that
/// is *inverted approximately* (one AMG V-cycle) to precondition the Helmholtz/PML operator
/// `K − k²M (+ damping)`. The usual choice is `β = 1 + 0.5i`; larger imaginary parts make the
/// multigrid easier and the preconditioner a worse approximation.
pub fn shifted_laplacian_operator<T: Scalar>(
    k: &Csr<T>,
    mass: Option<&Csr<T>>,
    k2: f64,
    beta: crate::linalg::scalar::C64,
) -> Csr<crate::linalg::scalar::C64> {
    let kc = k.to_complex();
    let shift = beta.scale(-k2);
    match mass {
        Some(m) => Csr::axpby(crate::linalg::scalar::C64::ONE, &kc, shift, &m.to_complex()),
        None => Csr::axpby(crate::linalg::scalar::C64::ONE, &kc, shift, &Csr::identity(k.nrows())),
    }
}
