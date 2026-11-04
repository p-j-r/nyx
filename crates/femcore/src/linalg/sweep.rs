//! Frequency sweeps by **reduced-order modelling** (port item W1): the response `x(ω)` of
//! `T(ω) x = b`, `T(ω) = Σ_j f_j(ω) A_j` (the split form of [`crate::linalg::nep`]: `K − ω²M + iωC`,
//! dispersive materials, frequency-dependent PML), for hundreds of frequencies from a handful of full
//! solves.
//!
//! Method: multi-point Galerkin projection (the reduced-basis / rational-Krylov idea). The full
//! solutions `x(ω_i)` at *expansion frequencies* span a basis `V`; the projected matrices
//! `A_j,r = Vᴴ A_j V` give the reduced model `T_r(ω) = Σ f_j(ω) A_j,r`, an `r × r` dense system per
//! frequency. Because the frequency dependence is affine in the terms, the reduced model is
//! evaluated at any `ω` without touching the large matrices. The expansion frequencies are chosen
//! *greedily*: evaluate the full-space relative residual `‖T(ω) V x_r − b‖/‖b‖` of the reduced
//! solution on the whole sweep grid (one pass of sparse products per frequency, no factorisation), add
//! the full solve at the worst frequency, and repeat until the worst residual is below the
//! tolerance. The residual is an *a-posteriori* check on every requested frequency, so a returned
//! sweep is certified up to the chosen tolerance (modulo the conditioning of `T(ω)`: a small residual
//! bounds the error by `‖T⁻¹‖·tol·‖b‖`, large near resonances).
//!
//! Each full solve is one numeric factorisation on a single symbolic analysis
//! ([`crate::linalg::direct::ShiftedFactorization::replace_matrices`]); [`direct_sweep`] is the baseline with one
//! factorisation per frequency. No Padé/Lanczos derivative matching (value matching only) and no
//! parametric (multi-parameter) models: see `docs/STATUS.md`.

#![allow(clippy::needless_range_loop)] // index loops over several parallel arrays (terms, Gram blocks)

use crate::linalg::csr::Csr;
use crate::linalg::dense::{Dense, solve_small};
use crate::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use crate::linalg::nep::{NepError, NonlinearEigenProblem};
use crate::linalg::scalar::{C64, Scalar};
use crate::linalg::vecops::norm2;

#[derive(Debug, Clone)]
pub struct SweepOptions {
    /// Target for the maximum relative residual over the sweep grid.
    pub tol: f64,
    /// Expansion frequencies before the greedy loop (evenly spread over the grid, at least 2).
    pub initial_points: usize,
    /// Stop (unconverged) after this many full solves.
    pub max_solves: usize,
    pub definiteness: Definiteness,
}

impl SweepOptions {
    pub fn new(tol: f64) -> Self {
        Self { tol, initial_points: 3, max_solves: 200, definiteness: Definiteness::General }
    }
}

#[derive(Debug, Clone)]
pub struct SweepResult {
    pub omegas: Vec<f64>,
    /// `outputs[i][k] = c_kᴴ x(ω_i)` for the requested functionals `c_k`.
    pub outputs: Vec<Vec<C64>>,
    /// Final relative residual of the reduced solution at every `ω_i`.
    pub residuals: Vec<f64>,
    /// Expansion frequencies in the order they were added.
    pub expansion: Vec<f64>,
    pub basis_dim: usize,
    pub converged: bool,
    /// Full-space numeric factorisations / symbolic analyses.
    pub factorizations: usize,
    pub symbolic_analyses: usize,
    /// Sparse residual passes over the grid (cost of the error indicator).
    pub residual_passes: usize,
    pub seconds: f64,
}

/// A reduced model: basis `V`, projected term matrices and the Gram matrices of the images `A_j V`,
/// which give the full-space residual norm of a reduced solution *without touching `n`-vectors*:
/// `‖Σ_j f_j A_j V x − b‖² = Σ_jk f̄_j f_k xᴴ G_jk x − 2 Re Σ_j f_j Σ_l x_l ḡ_jl + ‖b‖²` with
/// `G_jk = (A_jV)ᴴ(A_kV)`, `g_jl = (A_jV_l)ᴴ b`. This costs `terms²·r` dot products per added basis
/// vector and `terms²·r²` flops per frequency; it loses accuracy below about `1e-7` relative (squared
/// cancellation), so [`Rom::residual`] (exact) certifies the final result.
pub struct Rom<'a> {
    p: &'a NonlinearEigenProblem,
    /// Orthonormal basis columns (each of length `n`).
    v: Vec<Vec<C64>>,
    /// `av[j][i] = A_j V_i`.
    av: Vec<Vec<Vec<C64>>>,
    /// `a_r[j] = Vᴴ A_j V` (`r × r`).
    a_r: Vec<Dense<C64>>,
    /// `gram[j][k] = (A_jV)ᴴ(A_kV)`.
    gram: Vec<Vec<Dense<C64>>>,
    /// `gb[j][l] = (A_jV_l)ᴴ b`.
    gb: Vec<Vec<C64>>,
    b_r: Vec<C64>,
    b: &'a [C64],
    bnorm: f64,
}

impl<'a> Rom<'a> {
    pub fn new(p: &'a NonlinearEigenProblem, b: &'a [C64]) -> Self {
        let q = p.terms.len();
        Self {
            p,
            v: Vec::new(),
            av: vec![Vec::new(); q],
            a_r: (0..q).map(|_| Dense::zeros(0, 0)).collect(),
            gram: (0..q).map(|_| (0..q).map(|_| Dense::zeros(0, 0)).collect()).collect(),
            gb: vec![Vec::new(); q],
            b_r: Vec::new(),
            b,
            bnorm: norm2(b),
        }
    }

    pub fn dim(&self) -> usize {
        self.v.len()
    }

    /// Orthonormalises `x` against the basis and appends it (returns `false` if it was already in the span).
    pub fn extend(&mut self, x: &[C64]) -> bool {
        let mut w = x.to_vec();
        let nx = norm2(&w);
        for _ in 0..2 {
            for vi in &self.v {
                let d: C64 = vi.iter().zip(&w).map(|(a, b)| a.conj() * *b).sum();
                for (wi, vv) in w.iter_mut().zip(vi) {
                    *wi -= d * *vv;
                }
            }
        }
        let nw = norm2(&w);
        if nw <= 1e-13 * nx.max(1e-300) {
            return false;
        }
        w.iter_mut().for_each(|c| *c = c.scale(1.0 / nw));
        let (r, n, q) = (self.v.len(), w.len(), self.p.terms.len());
        let dot = |a: &[C64], b: &[C64]| -> C64 { a.iter().zip(b).map(|(x, y)| x.conj() * *y).sum() };
        let aw: Vec<Vec<C64>> = self
            .p
            .terms
            .iter()
            .map(|t| {
                let mut y = vec![C64::ZERO; n];
                t.matrix.spmv(&w, &mut y);
                y
            })
            .collect();
        // projected matrices
        for j in 0..q {
            let mut grown = Dense::<C64>::zeros(r + 1, r + 1);
            for i in 0..r {
                for l in 0..r {
                    grown[(i, l)] = self.a_r[j][(i, l)];
                }
                grown[(i, r)] = dot(&self.v[i], &aw[j]); // V_iᴴ A_j w
                grown[(r, i)] = dot(&w, &self.av[j][i]); // wᴴ A_j V_i
            }
            grown[(r, r)] = dot(&w, &aw[j]);
            self.a_r[j] = grown;
        }
        // Gram matrices of the images
        for j in 0..q {
            for k in 0..q {
                let mut grown = Dense::<C64>::zeros(r + 1, r + 1);
                for i in 0..r {
                    for l in 0..r {
                        grown[(i, l)] = self.gram[j][k][(i, l)];
                    }
                    grown[(i, r)] = dot(&self.av[j][i], &aw[k]);
                    grown[(r, i)] = dot(&aw[j], &self.av[k][i]);
                }
                grown[(r, r)] = dot(&aw[j], &aw[k]);
                self.gram[j][k] = grown;
            }
        }
        for j in 0..q {
            self.gb[j].push(dot(&aw[j], self.b));
        }
        for (j, y) in aw.into_iter().enumerate() {
            self.av[j].push(y);
        }
        self.b_r.push(dot(&w, self.b));
        self.v.push(w);
        true
    }

    fn coefficients(&self, omega: f64) -> Vec<C64> {
        self.p.terms.iter().map(|t| (t.f)(C64::new(omega, 0.0))).collect()
    }

    /// Reduced coordinates `x_r(ω)`; `None` if the reduced matrix is singular at `ω`.
    pub fn solve_reduced(&self, omega: f64) -> Option<Vec<C64>> {
        let r = self.dim();
        let fs = self.coefficients(omega);
        let t = Dense::from_fn(r, r, |i, j| self.a_r.iter().zip(&fs).map(|(a, f)| a[(i, j)] * *f).sum());
        solve_small(&t, &self.b_r)
    }

    /// `V x_r` (the approximate field).
    pub fn field(&self, xr: &[C64]) -> Vec<C64> {
        let n = self.v[0].len();
        let mut x = vec![C64::ZERO; n];
        for (vi, c) in self.v.iter().zip(xr) {
            for (xk, vv) in x.iter_mut().zip(vi) {
                *xk += *c * *vv;
            }
        }
        x
    }

    /// Exact full-space `‖T(ω)x − b‖ / ‖b‖` for the field of the reduced solution.
    pub fn residual(&self, omega: f64, xr: &[C64]) -> f64 {
        let x = self.field(xr);
        let mut y = vec![C64::ZERO; x.len()];
        self.p.apply(C64::new(omega, 0.0), &x, &mut y);
        let r: f64 = y.iter().zip(self.b).map(|(a, b)| (*a - *b).norm_sqr()).sum::<f64>().sqrt();
        r / self.bnorm.max(1e-300)
    }

    /// The same residual from the Gram matrices (cheap, accurate to about `1e-7` relative).
    pub fn residual_indicator(&self, omega: f64, xr: &[C64]) -> f64 {
        let (q, r) = (self.p.terms.len(), self.dim());
        let fs = self.coefficients(omega);
        let mut yy = 0.0;
        for j in 0..q {
            for k in 0..q {
                let g = &self.gram[j][k];
                let mut s = C64::ZERO;
                for i in 0..r {
                    let mut row = C64::ZERO;
                    for l in 0..r {
                        row += g[(i, l)] * xr[l];
                    }
                    s += xr[i].conj() * row;
                }
                yy += (fs[j].conj() * fs[k] * s).re;
            }
        }
        let mut by = 0.0;
        for j in 0..q {
            let s: C64 = (0..r).map(|l| xr[l] * self.gb[j][l].conj()).sum();
            by += (fs[j] * s).re;
        }
        let r2 = yy - 2.0 * by + self.bnorm * self.bnorm;
        r2.max(0.0).sqrt() / self.bnorm.max(1e-300)
    }
}

fn check_inputs(p: &NonlinearEigenProblem, b: &[C64], omegas: &[f64], outputs: &[Vec<C64>]) -> Result<(), NepError> {
    if b.len() != p.dim() || outputs.iter().any(|c| c.len() != p.dim()) {
        return Err(NepError::InvalidArgument("right-hand side and output vectors must have length n"));
    }
    if omegas.is_empty() || norm2(b) == 0.0 {
        return Err(NepError::InvalidArgument("need at least one frequency and a non-zero right-hand side"));
    }
    Ok(())
}

fn functional(c: &[C64], x: &[C64]) -> C64 {
    c.iter().zip(x).map(|(a, b)| a.conj() * *b).sum()
}

/// Full solve at `ω` through a persistent factorisation object.
struct FullSolver<'a> {
    p: &'a NonlinearEigenProblem,
    fac: Option<ShiftedFactorization<C64>>,
    definiteness: Definiteness,
}

impl FullSolver<'_> {
    fn solve(&mut self, omega: f64, b: &[C64]) -> Result<Vec<C64>, NepError> {
        let t: Csr<C64> = self.p.assemble(C64::new(omega, 0.0));
        let fac = match self.fac.as_mut() {
            Some(f) => {
                f.replace_matrices(&t, None)?;
                f
            }
            None => self.fac.insert(ShiftedFactorization::new(&t, None, self.definiteness)?),
        };
        fac.set_shift(C64::ZERO)?;
        let mut x = b.to_vec();
        fac.solve_block_in_place(&mut x, 1);
        Ok(x)
    }
    fn counts(&self) -> (usize, usize) {
        self.fac.as_ref().map_or((0, 0), |f| (f.stats().numeric_count, f.stats().symbolic_count))
    }
}

/// Baseline: one full factorisation and solve per frequency.
pub fn direct_sweep(
    p: &NonlinearEigenProblem,
    b: &[C64],
    omegas: &[f64],
    outputs: &[Vec<C64>],
    definiteness: Definiteness,
) -> Result<SweepResult, NepError> {
    check_inputs(p, b, omegas, outputs)?;
    let t0 = std::time::Instant::now();
    let mut s = FullSolver { p, fac: None, definiteness };
    let mut outs = Vec::with_capacity(omegas.len());
    let mut residuals = Vec::with_capacity(omegas.len());
    for &w in omegas {
        let x = s.solve(w, b)?;
        outs.push(outputs.iter().map(|c| functional(c, &x)).collect());
        let mut y = vec![C64::ZERO; x.len()];
        p.apply(C64::new(w, 0.0), &x, &mut y);
        residuals.push(y.iter().zip(b).map(|(a, bb)| (*a - *bb).norm_sqr()).sum::<f64>().sqrt() / norm2(b));
    }
    let (nf, ns) = s.counts();
    Ok(SweepResult {
        omegas: omegas.to_vec(),
        outputs: outs,
        residuals,
        expansion: omegas.to_vec(),
        basis_dim: 0,
        converged: true,
        factorizations: nf,
        symbolic_analyses: ns,
        residual_passes: 0,
        seconds: t0.elapsed().as_secs_f64(),
    })
}

/// Greedy reduced-order sweep (see the module docs). `omegas` is the grid on which the residual is
/// certified; `outputs` are the functionals `c_k` whose values `c_kᴴ x(ω)` are returned.
pub fn rom_sweep(
    p: &NonlinearEigenProblem,
    b: &[C64],
    omegas: &[f64],
    outputs: &[Vec<C64>],
    opts: &SweepOptions,
) -> Result<SweepResult, NepError> {
    check_inputs(p, b, omegas, outputs)?;
    let t0 = std::time::Instant::now();
    let mut full = FullSolver { p, fac: None, definiteness: opts.definiteness };
    let mut rom = Rom::new(p, b);
    let mut expansion: Vec<f64> = Vec::new();
    let add = |w: f64, rom: &mut Rom, full: &mut FullSolver, expansion: &mut Vec<f64>| -> Result<bool, NepError> {
        let x = full.solve(w, b)?;
        expansion.push(w);
        Ok(rom.extend(&x))
    };
    // initial expansion frequencies spread over the grid
    let m = opts.initial_points.clamp(2, omegas.len().max(2));
    for k in 0..m {
        let idx = if m == 1 { 0 } else { k * (omegas.len() - 1) / (m - 1) };
        if !expansion.contains(&omegas[idx]) {
            add(omegas[idx], &mut rom, &mut full, &mut expansion)?;
        }
    }
    let mut passes = 0usize;
    let mut exact = false; // switches on after the cheap indicator first claims convergence
    let (mut residuals, mut converged);
    loop {
        passes += 1;
        residuals = Vec::with_capacity(omegas.len());
        let mut worst = (0usize, 0.0f64);
        for (i, &w) in omegas.iter().enumerate() {
            let r = match rom.solve_reduced(w) {
                Some(xr) if exact => rom.residual(w, &xr),
                Some(xr) => rom.residual_indicator(w, &xr),
                None => f64::INFINITY,
            };
            if r > worst.1 || (r.is_nan() && worst.1.is_finite()) {
                worst = (i, if r.is_nan() { f64::INFINITY } else { r });
            }
            residuals.push(r);
        }
        converged = worst.1 <= opts.tol;
        if converged && !exact {
            // certify with the exact residuals; continue greedily with them if the indicator was optimistic
            exact = true;
            continue;
        }
        if converged || expansion.len() >= opts.max_solves {
            break;
        }
        let w = omegas[worst.0];
        if expansion.contains(&w) && !exact {
            // the indicator has reached its noise floor (≈ 1e-8): continue with exact residuals
            exact = true;
            continue;
        }
        if expansion.contains(&w) || !add(w, &mut rom, &mut full, &mut expansion)? {
            break; // the worst frequency is already an expansion point (or adds nothing): stagnation
        }
    }
    let mut outs = Vec::with_capacity(omegas.len());
    for &w in omegas {
        let x = rom.solve_reduced(w).map(|xr| rom.field(&xr));
        outs.push(match x {
            Some(x) => outputs.iter().map(|c| functional(c, &x)).collect(),
            None => vec![C64::new(f64::NAN, f64::NAN); outputs.len()],
        });
    }
    let (nf, ns) = full.counts();
    Ok(SweepResult {
        omegas: omegas.to_vec(),
        outputs: outs,
        residuals,
        expansion,
        basis_dim: rom.dim(),
        converged,
        factorizations: nf,
        symbolic_analyses: ns,
        residual_passes: passes,
        seconds: t0.elapsed().as_secs_f64(),
    })
}
