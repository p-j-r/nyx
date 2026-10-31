//! Eigenvalue sensitivities and mode tracking (port item N3; groundwork for the optimisation of
//! item O1).
//!
//! **Derivatives.** For a simple eigenvalue of `A(p) x = λ B(p) x` with right vector `x` and left
//! vector `y` (`yᴴ A = λ yᴴ B`),
//!
//! ```text
//! dλ/dp = yᴴ (A' − λ B') x / (yᴴ B x).
//! ```
//!
//! For a nonlinear problem `T(λ; p) x = 0` (the split form of [`crate::solvers::nep`], including the quadratic
//! `λ²M + λC + K`): `dλ/dp = − yᴴ (∂T/∂p) x / (yᴴ T'(λ) x)` with `T'(λ) = Σ f_j'(λ) A_j`.
//! The cost is one left-vector solve per eigenpair plus a few sparse products *per parameter*, so the
//! gradient of an eigenvalue objective with many design parameters is cheap (the adjoint property).
//! Left vectors come from [`left_eigenvector`] (inverse iteration on the conjugate transpose, one
//! factorisation) or are free for symmetric structure: `y = x` (Hermitian, `B` Hermitian) and `y = x̄`
//! (complex symmetric, `A` and `B` symmetric).
//!
//! **Tracking.** [`match_modes`] assigns the modes of the new parameter value to the modes of the
//! previous one by maximal `B`-weighted overlap `|x_oldᴴ B x_new|² / (‖x_old‖_B² ‖x_new‖_B²)` with the
//! Hungarian algorithm, and flags assignments whose overlap is low or ambiguous (crossings).

use crate::solvers::csr::Csr;
use crate::solvers::direct::{Definiteness, ShiftedFactorization, ShiftedSolver, SolverError};
use crate::solvers::nep::NonlinearEigenProblem;
use crate::solvers::rng::random_vector;
use crate::solvers::scalar::{C64, Scalar};
use crate::solvers::vecops::{dotc, norm2};

fn conj_transpose(a: &Csr<C64>) -> Csr<C64> {
    a.transpose().map(|v| v.conj())
}

/// Left eigenvector of `A x = λ B x`: `Aᴴ y = λ̄ Bᴴ y`, by inverse iteration with the
/// factorisation of `Aᴴ − σ̄ Bᴴ` at a shift slightly off `λ` (relative `1e-6`), starting from `start` (e.g. the
/// right vector, which is already a good guess for nearly normal problems).
pub fn left_eigenvector(
    a: &Csr<C64>,
    b: Option<&Csr<C64>>,
    lambda: C64,
    start: &[C64],
    iterations: usize,
) -> Result<Vec<C64>, SolverError> {
    let n = a.nrows();
    let ah = conj_transpose(a);
    let bh = b.map(conj_transpose);
    let sigma = lambda + C64::new(1e-6 * lambda.norm().max(1.0), 1e-6 * lambda.norm().max(1.0));
    let mut fac = ShiftedFactorization::new(&ah, bh.as_ref(), Definiteness::General)?;
    fac.set_shift(sigma.conj())?;
    let mut y = if start.len() == n { start.to_vec() } else { random_vector::<C64>(n, 0x1EF7) };
    for _ in 0..iterations.max(1) {
        // y ← (Aᴴ − σ̄Bᴴ)⁻¹ Bᴴ y
        let mut rhs = match &bh {
            Some(bh) => {
                let mut t = vec![C64::ZERO; n];
                bh.spmv(&y, &mut t);
                t
            }
            None => y.clone(),
        };
        fac.solve_block_in_place(&mut rhs, 1);
        let nr = norm2(&rhs);
        y = rhs.into_iter().map(|v| v.scale(1.0 / nr)).collect();
    }
    Ok(y)
}

/// `dλ/dp` for the generalised linear problem, given `dA/dp` and optionally `dB/dp`.
pub fn eigenvalue_derivative(
    lambda: C64,
    x: &[C64],
    y: &[C64],
    b: Option<&Csr<C64>>,
    da: &Csr<C64>,
    db: Option<&Csr<C64>>,
) -> C64 {
    let n = x.len();
    let mut t = vec![C64::ZERO; n];
    da.spmv(x, &mut t);
    if let Some(db) = db {
        let mut u = vec![C64::ZERO; n];
        db.spmv(x, &mut u);
        for (ti, ui) in t.iter_mut().zip(&u) {
            *ti -= lambda * *ui;
        }
    }
    let num = dotc(y, &t);
    let den = match b {
        Some(b) => {
            let mut u = vec![C64::ZERO; n];
            b.spmv(x, &mut u);
            dotc(y, &u)
        }
        None => dotc(y, x),
    };
    num / den
}

/// `dλ/dp = −yᴴ (∂T/∂p) x / (yᴴ T'(λ) x)` for `T(λ;p) = Σ f_j(λ) A_j(p)`; `dt` holds the split form of `∂T/∂p`
/// (its own `f` functions, evaluated at `λ`). `T'(λ)` uses a central difference of the `f_j` of `p`
/// (complex step `h = 1e-6·max(1,|λ|)`), accurate to about 1e-9 relative for smooth `f_j`.
pub fn nonlinear_eigenvalue_derivative(
    p: &NonlinearEigenProblem,
    dt: &NonlinearEigenProblem,
    lambda: C64,
    x: &[C64],
    y: &[C64],
) -> C64 {
    let n = x.len();
    let mut dtx = vec![C64::ZERO; n];
    dt.apply(lambda, x, &mut dtx);
    let h = C64::new(1e-6 * lambda.norm().max(1.0), 0.0);
    let (mut tp, mut tm) = (vec![C64::ZERO; n], vec![C64::ZERO; n]);
    p.apply(lambda + h, x, &mut tp);
    p.apply(lambda - h, x, &mut tm);
    let tpx: Vec<C64> = tp.iter().zip(&tm).map(|(a, b)| (*a - *b) / (h * 2.0)).collect();
    -dotc(y, &dtx) / dotc(y, &tpx)
}

/// Result of matching new modes to previous ones.
#[derive(Debug, Clone)]
pub struct ModeMatch {
    /// `assignment[i]` = index of the new mode matched to old mode `i`.
    pub assignment: Vec<usize>,
    /// Overlap (in `[0, 1]`) of every matched pair.
    pub overlap: Vec<f64>,
    /// Matches with overlap below `min_overlap`, or whose runner-up overlap is within `ambiguity`× of the best.
    pub uncertain: Vec<bool>,
}

/// Hungarian algorithm (Kuhn–Munkres, `O(n³)`) minimising the total cost of a square assignment.
fn hungarian(cost: &[Vec<f64>]) -> Vec<usize> {
    let n = cost.len();
    let (mut u, mut v) = (vec![0.0; n + 1], vec![0.0; n + 1]);
    let (mut p, mut way) = (vec![0usize; n + 1], vec![0usize; n + 1]);
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![f64::INFINITY; n + 1];
        let mut used = vec![false; n + 1];
        loop {
            used[j0] = true;
            let (i0, mut delta, mut j1) = (p[j0], f64::INFINITY, 0);
            for j in 1..=n {
                if !used[j] {
                    let cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
                    if cur < minv[j] {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if minv[j] < delta {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut assignment = vec![0usize; n];
    for j in 1..=n {
        assignment[p[j] - 1] = j - 1;
    }
    assignment
}

/// Matches `new` modes (columns of length `n`) to `old` modes by `B`-overlap (`b = None`: Euclidean).
/// `old` and `new` need the same number of modes.
pub fn match_modes(
    old: &[Vec<C64>],
    new: &[Vec<C64>],
    b: Option<&Csr<C64>>,
    min_overlap: f64,
    ambiguity: f64,
) -> ModeMatch {
    let m = old.len();
    assert_eq!(new.len(), m, "old and new must hold the same number of modes");
    let weigh = |x: &[C64]| -> Vec<C64> {
        match b {
            Some(b) => {
                let mut t = vec![C64::ZERO; x.len()];
                b.spmv(x, &mut t);
                t
            }
            None => x.to_vec(),
        }
    };
    let bnew: Vec<Vec<C64>> = new.iter().map(|x| weigh(x)).collect();
    let nrm_new: Vec<f64> = new.iter().zip(&bnew).map(|(x, bx)| dotc(x, bx).re.max(0.0)).collect();
    let nrm_old: Vec<f64> = old.iter().map(|x| dotc(x, &weigh(x)).re.max(0.0)).collect();
    let mut overlap = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in 0..m {
            let d = dotc(&old[i], &bnew[j]);
            overlap[i][j] = d.norm_sqr() / (nrm_old[i] * nrm_new[j]).max(1e-300);
        }
    }
    let cost: Vec<Vec<f64>> = overlap.iter().map(|r| r.iter().map(|o| 1.0 - o).collect()).collect();
    let assignment = if m == 0 { Vec::new() } else { hungarian(&cost) };
    let mut ov = Vec::with_capacity(m);
    let mut uncertain = Vec::with_capacity(m);
    for i in 0..m {
        let best = overlap[i][assignment[i]];
        let runner = (0..m).filter(|&j| j != assignment[i]).map(|j| overlap[i][j]).fold(0.0, f64::max);
        ov.push(best);
        uncertain.push(best < min_overlap || runner * ambiguity > best);
    }
    ModeMatch { assignment, overlap: ov, uncertain }
}

/// Follows modes along a parameter sweep. Feed the modes (values and vectors) at each parameter step to
/// [`Tracker::step`]; the returned order is stable: entry `i` is always the same physical mode.
pub struct Tracker {
    prev: Vec<Vec<C64>>,
    pub min_overlap: f64,
    pub ambiguity: f64,
    /// Steps on which at least one match was flagged uncertain.
    pub uncertain_steps: Vec<usize>,
    steps: usize,
}

impl Tracker {
    pub fn new(initial: Vec<Vec<C64>>) -> Self {
        Self { prev: initial, min_overlap: 0.5, ambiguity: 2.0, uncertain_steps: Vec::new(), steps: 0 }
    }

    /// Reorders `values`/`vectors` of the new step to follow the tracked modes; returns the match.
    pub fn step(&mut self, values: &mut [C64], vectors: &mut [Vec<C64>], b: Option<&Csr<C64>>) -> ModeMatch {
        self.steps += 1;
        let mm = match_modes(&self.prev, vectors, b, self.min_overlap, self.ambiguity);
        let (v0, x0) = (values.to_vec(), vectors.to_vec());
        for (i, &j) in mm.assignment.iter().enumerate() {
            values[i] = v0[j];
            vectors[i] = x0[j].clone();
        }
        if mm.uncertain.iter().any(|&u| u) {
            self.uncertain_steps.push(self.steps);
        }
        self.prev = vectors.to_vec();
        mm
    }
}
