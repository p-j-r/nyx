// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Frequency-domain nonlinear optics solvers (`Study::NonlinearDriven` and `Study::NonlinearModes`, ADR 143):
//!
//! 1. **Single-tone $\chi^{(3)}$ Kerr + TPA + FCA/FCD Harmonic Balance FEM (`process = "chi3"`)**:
//!    - Solves the real $2N_{\text{red}} \times 2N_{\text{red}}$ split $(\operatorname{Re}\hat{\mathbf{y}}, \operatorname{Im}\hat{\mathbf{y}})$
//!      curl-curl system
//!      $$\hat{\mathbf{F}}(\hat{\mathbf{y}}, \alpha) = \hat{\mathbf{A}}(\omega)\hat{\mathbf{y}} - k_0^2 \hat{\mathbf{g}}_{\text{NL}}(\hat{\mathbf{y}}) - \alpha \hat{\mathbf{b}}_{\text{ref}} = \mathbf{0}$$
//!      with full cross-polarization vector Kerr susceptibility
//!      $\mathbf{P}^{(3,\text{Kerr})}(\omega) = \varepsilon_0\big[A(\mathbf{E}\cdot\mathbf{E}^*)\mathbf{E} + \frac{1}{2}B(\mathbf{E}\cdot\mathbf{E})\mathbf{E}^*\big]$
//!      (plus rotated cubic `m-3m` anisotropy), two-photon absorption $\beta_{\text{TPA}}$, and steady-state
//!      TPA-generated free-carrier dispersion/absorption $N_c(\mathbf{E}) = \frac{\tau_c \beta_{\text{TPA}}}{2\hbar\omega} I(\mathbf{E})^2$.
//!    - Exact Wirtinger Fréchet Jacobian $\mathbf{J} = \partial\hat{\mathbf{F}}/\partial(\operatorname{Re}\hat{\mathbf{y}}, \operatorname{Im}\hat{\mathbf{y}})$
//!      reusing a single static CSR symbolic analysis across `femcore::linalg::snes::newton` (or `snes::jfnk`).
//!    - Supports both monotonic `power_ramp` continuation and bordered $(2N_{\text{red}} + 1)$ Keller
//!      **pseudo-arclength continuation** (`continuation = "arclength"`) tracing S-shaped optical bistability
//!      curves cleanly across upper and lower saddle-node turning points.
//!
//! 2. **Two-tone $\chi^{(3)}$ Cross-Phase Modulation (`process = "chi3_xpm"`)**:
//!    - Coupled $(\omega_1, \omega_2)$ Harmonic Balance FEM capturing both co-polarized ($2 n_2 I_1$) and
//!      orthogonally polarized ($\frac{2}{3} n_2 I_1$) Maker–Terhune cross-phase modulation and two-photon cross-absorption.
//!
//! 3. **Multi-Harmonic $\chi^{(2)}$ Coupled-Wave FEM (`process = "shg" | "sfg" | "dfg" | "opa"`)**:
//!    - Degenerate second-harmonic generation ($\omega_2 = 2\omega_1$) and 3-wave mixing ($\omega_3 = \omega_1 + \omega_2$)
//!      with chromatically dispersive `[optical]` evaluation at each harmonic $\lambda_m = c_0 / f_m$,
//!      Gauss-quadrature quasi-phase-matching (`QpmSpec::poling_sign_at`), fast `undepleted` sequential solves,
//!      and exact `depleted` coupled Newton (`snes::newton` / `snes::jfnk`) satisfying discrete Manley–Rowe
//!      photon-flux conservation to solver tolerance.
//!
//! 4. **2D Cross-Section Nonlinear Eigenmodes (`Study::NonlinearModes`)**:
//!    - Full-vectorial waveguide nonlinear parameter $\gamma$ ($\text{W}^{-1}\text{m}^{-1}$), effective nonlinear
//!      mode area $A_{\text{eff}}$ ($\text{m}^2$), TPA figure of merit, and fixed-power bordered Newton continuation
//!      computing power-dependent complex effective index $n_{\text{eff}}(P)$ on $(\mathbf{e}_t, \tilde{\boldsymbol{\phi}}, \beta)$.

use crate::cavity::{basis3_keyed, safe_curved_jacobian};
use crate::driven::{assemble_emw_driven, DrivenAssembly};
use crate::ibc::{EPS0, MU0};
use crate::modes::C0;
use crate::modes_complex::{
    assemble_modes_complex_with_bcs_and_symmetry_and_sheets, mode_fields_complex_cells, solve_mode_pairs_complex,
    ModeAssemblyC, ModesC,
};
use crate::problem::{Problem, Study};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, quad_rule, shape, CurvedTet};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, RCsr, C};
use femcore::linalg::direct::{Definiteness, ShiftedFactorization, ShiftedSolver};
use femcore::linalg::linop::Preconditioner;
use femcore::linalg::snes::{jfnk, newton, JfnkOptions, LineSearch, NewtonOptions};
use femcore::material::{Chi3, QpmSpec, HBAR};
use femcore::optical::{cell_materials, EpsC, MuC};
use nyx_physics::{cell_fields_of, Discretisation, DrivenPortData, Field, Scalar, Solution, SolveResult};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::time::Instant;

/// Allowed choices for `Study::NonlinearDriven` `process`.
pub const NONLINEAR_PROCESS_CHOICES: [&str; 7] = ["auto", "chi3", "chi3_xpm", "shg", "sfg", "dfg", "opa"];

/// Allowed choices for `Study::NonlinearDriven` `continuation`.
pub const NONLINEAR_CONTINUATION_CHOICES: [&str; 2] = ["power_ramp", "arclength"];

/// Allowed choices for `Study::NonlinearDriven` `coupling`.
pub const CHI2_COUPLING_CHOICES: [&str; 2] = ["depleted", "undepleted"];

/// Allowed choices for `Study::NonlinearDriven` `nl_solver`.
pub const NONLINEAR_SOLVER_CHOICES: [&str; 2] = ["newton", "jfnk"];

/// Nonlinear optical process solved by `Study::NonlinearDriven`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NonlinearProcess {
    #[default]
    Auto,
    Chi3,
    Chi3Xpm,
    Shg,
    Sfg,
    Dfg,
    Opa,
}

impl NonlinearProcess {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "auto" | "" => Ok(Self::Auto),
            "chi3" => Ok(Self::Chi3),
            "chi3_xpm" => Ok(Self::Chi3Xpm),
            "shg" => Ok(Self::Shg),
            "sfg" => Ok(Self::Sfg),
            "dfg" => Ok(Self::Dfg),
            "opa" => Ok(Self::Opa),
            other => {
                Err(format!("unknown nonlinear process '{other}' (expected one of {NONLINEAR_PROCESS_CHOICES:?})"))
            }
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Chi3 => "chi3",
            Self::Chi3Xpm => "chi3_xpm",
            Self::Shg => "shg",
            Self::Sfg => "sfg",
            Self::Dfg => "dfg",
            Self::Opa => "opa",
        }
    }
}

/// Continuation strategy for `Study::NonlinearDriven`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NonlinearContinuation {
    #[default]
    PowerRamp,
    Arclength,
}

impl NonlinearContinuation {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "power_ramp" | "" => Ok(Self::PowerRamp),
            "arclength" => Ok(Self::Arclength),
            other => {
                Err(format!("unknown continuation '{other}' (expected one of {NONLINEAR_CONTINUATION_CHOICES:?})"))
            }
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::PowerRamp => "power_ramp",
            Self::Arclength => "arclength",
        }
    }
}

/// Depleted vs undepleted coupling mode for multi-harmonic $\chi^{(2)}$ and XPM solves.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Chi2CouplingMode {
    #[default]
    Depleted,
    Undepleted,
}

impl Chi2CouplingMode {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "depleted" | "" => Ok(Self::Depleted),
            "undepleted" => Ok(Self::Undepleted),
            other => Err(format!("unknown coupling '{other}' (expected one of {CHI2_COUPLING_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Depleted => "depleted",
            Self::Undepleted => "undepleted",
        }
    }
}

/// Nonlinear algebraic solver backend (`snes::newton` vs preconditioned `snes::jfnk`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum NonlinearSolverKind {
    #[default]
    Newton,
    Jfnk,
}

impl NonlinearSolverKind {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "newton" | "" => Ok(Self::Newton),
            "jfnk" => Ok(Self::Jfnk),
            other => Err(format!("unknown nl_solver '{other}' (expected one of {NONLINEAR_SOLVER_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Newton => "newton",
            Self::Jfnk => "jfnk",
        }
    }
}

/// Typed parameters of `Study::NonlinearDriven` (`"nonlinear_driven"`).
#[derive(Clone, Debug, PartialEq)]
pub struct NonlinearDrivenSpec {
    pub process: NonlinearProcess,
    /// Primary pump frequency $f_1 = \omega_1 / (2\pi)$ in Hz.
    pub freq: f64,
    /// Optional second frequency $f_2$ in Hz for `sfg` ($\omega_3 = \omega_1 + \omega_2$), `dfg`/`opa` ($\omega_1, \omega_2$ with pump $\omega_3 = \omega_1 + \omega_2$), or `chi3_xpm`.
    pub freq_2: f64,
    /// Nominal input pump power $P_{\text{in}}$ in Watts (used when `powers_w` is empty and `n_steps == 1`).
    pub power_w: f64,
    /// Optional explicit sequence of input powers in Watts for `power_ramp` sweeps.
    pub powers_w: Vec<f64>,
    /// Minimum sweep power in Watts when `n_steps > 1` and `powers_w` is empty.
    pub power_min_w: f64,
    /// Maximum sweep power in Watts (or target power ceiling for `arclength` continuation).
    pub power_max_w: f64,
    /// Number of power steps for `power_ramp` when `powers_w` is empty (default `1`).
    pub n_steps: usize,
    /// 0-based driven port index receiving the primary pump (default `0`).
    pub drive_port: usize,
    /// Optional per-port incident powers (W) at $\omega_1$ overriding `power_w` on `drive_port` (scaled proportionally during sweeps).
    pub port_powers_w: Vec<f64>,
    /// Optional per-port incident phases (rad) at $\omega_1$.
    pub port_phases_rad: Vec<f64>,
    /// Secondary/signal input power $P_2$ in Watts for `sfg`, `dfg`, `opa`, or `chi3_xpm`.
    pub signal_power_w: f64,
    /// 0-based port index for the secondary/signal wave (default `0`).
    pub signal_port: usize,
    /// Secondary/signal input phase in radians.
    pub signal_phase_rad: f64,
    pub continuation: NonlinearContinuation,
    pub coupling: Chi2CouplingMode,
    pub nl_solver: NonlinearSolverKind,
    pub rtol: f64,
    pub max_it: usize,
    /// Dimensionless initial pseudo-arclength step $\Delta s$ (default `0.08`).
    pub arclength_ds: f64,
    /// Maximum number of pseudo-arclength continuation steps (default `240`).
    pub arclength_max_steps: usize,
}

impl Default for NonlinearDrivenSpec {
    fn default() -> Self {
        Self {
            process: NonlinearProcess::Auto,
            freq: 1.934_144_89e14,
            freq_2: 0.0,
            power_w: 1.0,
            powers_w: Vec::new(),
            power_min_w: 0.0,
            power_max_w: 0.0,
            n_steps: 1,
            drive_port: 0,
            port_powers_w: Vec::new(),
            port_phases_rad: Vec::new(),
            signal_power_w: 0.0,
            signal_port: 0,
            signal_phase_rad: 0.0,
            continuation: NonlinearContinuation::PowerRamp,
            coupling: Chi2CouplingMode::Depleted,
            nl_solver: NonlinearSolverKind::Newton,
            rtol: 1e-10,
            max_it: 40,
            arclength_ds: 0.08,
            arclength_max_steps: 240,
        }
    }
}

impl NonlinearDrivenSpec {
    pub fn is_finite(&self) -> bool {
        self.freq.is_finite()
            && self.freq_2.is_finite()
            && self.power_w.is_finite()
            && self.power_min_w.is_finite()
            && self.power_max_w.is_finite()
            && self.signal_power_w.is_finite()
            && self.signal_phase_rad.is_finite()
            && self.rtol.is_finite()
            && self.arclength_ds.is_finite()
            && self.powers_w.iter().all(|v| v.is_finite())
            && self.port_powers_w.iter().all(|v| v.is_finite())
            && self.port_phases_rad.iter().all(|v| v.is_finite())
    }
}

/// Typed parameters of `Study::NonlinearModes` (`"nonlinear_modes"`).
#[derive(Clone, Debug, PartialEq)]
pub struct NonlinearModesSpec {
    pub freq: f64,
    pub n_target: f64,
    pub evanescent: bool,
    /// 1-based mode index from the linear eigensolve to track across optical power `P` (default `1`).
    pub mode_index: usize,
    /// Single target optical power $P$ in Watts when `powers_w` is empty (default `1.0` W).
    pub power_w: f64,
    /// Optional sequence of guided optical powers $P_k \ge 0$ in Watts for continuation.
    pub powers_w: Vec<f64>,
    pub rtol: f64,
    pub max_it: usize,
}

impl Default for NonlinearModesSpec {
    fn default() -> Self {
        Self {
            freq: 1.934_144_89e14,
            n_target: 1.8,
            evanescent: false,
            mode_index: 1,
            power_w: 1.0,
            powers_w: Vec::new(),
            rtol: 1e-10,
            max_it: 30,
        }
    }
}

impl NonlinearModesSpec {
    pub fn is_finite(&self) -> bool {
        self.freq.is_finite()
            && self.n_target.is_finite()
            && self.power_w.is_finite()
            && self.rtol.is_finite()
            && self.powers_w.iter().all(|v| v.is_finite())
    }
}

// =============================================================================================
// Pre-rotated per-domain nonlinear constitutive data & CSR fast-update pattern
// =============================================================================================

#[derive(Clone, Debug)]
struct DomainChi2Data {
    /// Global-frame rank-3 $d^{\text{glob}}_{ijk}$ ($\text{m/V}$, symmetric in `(j, k)`).
    d_rank3: [[[f64; 3]; 3]; 3],
    qpm: QpmSpec,
}

#[derive(Clone, Debug)]
struct DomainChi3Data {
    chi3: Chi3,
    n0: f64,
    euler_deg: [f64; 3],
    has_cubic: bool,
    a_mt: f64,
    b_mt: f64,
    rank4_glob: [[[[f64; 3]; 3]; 3]; 3],
}

/// Precomputed cell basis evaluations at Gauss quadrature points for fast deterministic assembly.
#[derive(Clone, Debug)]
struct TetQuadPoint3d {
    wq_vol: f64,
    pos: [f64; 3],
    w_basis: Vec<[f64; 3]>,
}

#[derive(Clone, Debug)]
struct ActiveCell3d {
    dom: i32,
    /// Reduced DOF index (`usize::MAX` if fixed PEC) and complex Floquet*Ruiz factor `d_ruiz[r] * ratio[full_dof]` per local basis slot.
    red_map: Vec<(usize, C)>,
    qpts: Vec<TetQuadPoint3d>,
}

/// Static sorted CSR sparsity pattern allowing $O(\text{nnz}\log d)$ zero-allocation numeric assembly
/// with constant `indptr` and `indices` across every Newton and continuation step.
#[derive(Clone, Debug)]
pub struct StaticRealPattern {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
}

impl StaticRealPattern {
    fn from_pairs(nrows: usize, ncols: usize, pairs: impl IntoIterator<Item = (usize, usize)>) -> Self {
        let mut rows: Vec<Vec<usize>> = vec![Vec::new(); nrows];
        for (r, c_idx) in pairs {
            if r < nrows && c_idx < ncols {
                rows[r].push(c_idx);
            }
        }
        let mut indptr = Vec::with_capacity(nrows + 1);
        let mut indices = Vec::new();
        indptr.push(0);
        for row in &mut rows {
            row.sort_unstable();
            row.dedup();
            indices.extend_from_slice(row);
            indptr.push(indices.len());
        }
        Self { nrows, ncols, indptr, indices }
    }

    #[inline]
    pub fn zero_values(&self) -> Vec<f64> {
        vec![0.0; self.indices.len()]
    }

    #[inline]
    fn add_entry(&self, vals: &mut [f64], row: usize, col: usize, v: f64) {
        if v == 0.0 {
            return;
        }
        let start = self.indptr[row];
        let end = self.indptr[row + 1];
        if let Ok(off) = self.indices[start..end].binary_search(&col) {
            vals[start + off] += v;
        }
    }

    /// Adds a complex Wirtinger pair `(j1, j2)` at complex block `(r, s)` into the `2x2` real split slots `(2r..2r+2, 2s..2s+2)`.
    #[inline]
    fn add_wirtinger_2x2(&self, vals: &mut [f64], r: usize, s: usize, j1: C, j2: C) {
        let rr = 2 * r;
        let cc = 2 * s;
        let a00 = (j1 + j2).re;
        let a01 = -(j1 - j2).im;
        let a10 = (j1 + j2).im;
        let a11 = (j1 - j2).re;
        self.add_entry(vals, rr, cc, a00);
        self.add_entry(vals, rr, cc + 1, a01);
        self.add_entry(vals, rr + 1, cc, a10);
        self.add_entry(vals, rr + 1, cc + 1, a11);
    }

    pub fn into_csr(&self, vals: Vec<f64>) -> RCsr {
        RCsr::from_raw(self.nrows, self.ncols, self.indptr.clone(), self.indices.clone(), vals)
    }
}

/// Linear LU preconditioner wrapping a `ShiftedFactorization<f64>` of the linear Jacobian $\mathbf{J}(\mathbf{0})$ for `snes::jfnk`.
struct FactoredRealPc {
    fac: ShiftedFactorization<f64>,
}

impl Preconditioner<f64> for FactoredRealPc {
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.copy_from_slice(x);
        let _ = self.fac.solve_block_in_place(y, 1);
    }
}

#[inline]
fn join_real_to_complex(r: &[f64]) -> Vec<C> {
    let n = r.len() / 2;
    (0..n).map(|i| c(r[2 * i], r[2 * i + 1])).collect()
}

fn clone_problem_at_freq(problem: &Problem, freq_hz: f64) -> Result<Problem, String> {
    let orig_freq = problem.study.reference_freq();
    let mut p = problem.clone();
    p.study = Study::Driven {
        f_min: freq_hz,
        f_max: freq_hz,
        n_points: 1,
        freqs_hz: Some(vec![freq_hz]),
        prom_tol: None,
        max_hdm_solves: 4,
        strict_prom_tol: false,
        seed_freqs_hz: Vec::new(),
        solver: "direct".into(),
        blr_tol: 1e-6,
    };
    // Re-evaluate library [optical] permittivity when cloning to a distinct harmonic frequency
    // (omega_m != orig_freq) so multi-harmonic SHG/SFG/DFG/OPA captures chromatic dispersion n(omega_m) != n(omega_1),
    // while preserving the primary study's decoded permittivity (including any explicit Ctx::wavelength) at orig_freq.
    let lambda_m = C0 / freq_hz.max(1.0);
    let temp_k = problem.fallback_temperature;
    if (freq_hz - orig_freq).abs() > 1e-12 * orig_freq.max(1.0) && lambda_m <= 3e-4 {
        for (&dom, td) in &problem.thermo_optic_domains {
            if td.spec.optical.is_some() {
                let eps_cryst = td.spec.optical_permittivity(lambda_m, temp_k, [0.0; 3]).map_err(|e| e.to_string())?;
                let eps_glob =
                    td.spec.optical_permittivity(lambda_m, temp_k, td.euler_deg).map_err(|e| e.to_string())?;
                let re_cryst = [
                    [eps_cryst[0][0].re, eps_cryst[0][1].re, eps_cryst[0][2].re],
                    [eps_cryst[1][0].re, eps_cryst[1][1].re, eps_cryst[1][2].re],
                    [eps_cryst[2][0].re, eps_cryst[2][1].re, eps_cryst[2][2].re],
                ];
                let re_glob = [
                    [eps_glob[0][0].re, eps_glob[0][1].re, eps_glob[0][2].re],
                    [eps_glob[1][0].re, eps_glob[1][1].re, eps_glob[1][2].re],
                    [eps_glob[2][0].re, eps_glob[2][1].re, eps_glob[2][2].re],
                ];
                let scale = re_glob[0][0].abs().max(re_glob[1][1].abs()).max(re_glob[2][2].abs()).max(1.0);
                let iso = (re_glob[0][0] - re_glob[1][1]).abs()
                    + (re_glob[1][1] - re_glob[2][2]).abs()
                    + re_glob[0][1].abs()
                    + re_glob[0][2].abs()
                    + re_glob[1][2].abs()
                    < 1e-12 * scale;
                if iso {
                    p.materials.set(dom, "eps_r", re_glob[0][0]);
                    if let Some(an) = p.materials.aniso.get_mut(&dom) {
                        an.permittivity = None;
                    }
                } else {
                    p.materials.aniso.entry(dom).or_default().permittivity =
                        Some(femcore::material::Tensor::Symmetric {
                            xx: re_cryst[0][0],
                            yy: re_cryst[1][1],
                            zz: re_cryst[2][2],
                            xy: 0.5 * (re_cryst[0][1] + re_cryst[1][0]),
                            yz: 0.5 * (re_cryst[1][2] + re_cryst[2][1]),
                            xz: 0.5 * (re_cryst[0][2] + re_cryst[2][0]),
                        });
                }
            }
        }
    }
    Ok(p)
}

fn build_active_cells_3d(
    mesh: &Discretisation,
    asm: &DrivenAssembly,
    active_domains: &BTreeSet<i32>,
) -> Vec<ActiveCell3d> {
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    let nb = asm.ned.nb();
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];
    let mut out = Vec::new();

    for ci in 0..mesh.n_cells() {
        if asm.pml_cells.get(ci).copied().unwrap_or(false) {
            continue;
        }
        let dom = mesh.cell_labels[ci];
        if !active_domains.contains(&dom) {
            continue;
        }
        let cdofs = &asm.ned.cell_dofs[ci];
        let red_map: Vec<(usize, C)> = cdofs
            .iter()
            .map(|&gd| {
                let r = asm.red.index_of[gd];
                if r == usize::MAX {
                    (usize::MAX, c(0.0, 0.0))
                } else {
                    (r, asm.red.ratio[gd] * asm.d_ruiz[r])
                }
            })
            .collect();
        if red_map.iter().all(|&(r, _)| r == usize::MAX) {
            continue;
        }
        let gid = asm.ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for k in 0..3 {
                gs[b][k] = geo.g[ci * 12 + b * 3 + k];
            }
        }
        let v_tet = mesh.cell(ci);
        let pts_tet = [
            mesh.point(v_tet[0] as usize),
            mesh.point(v_tet[1] as usize),
            mesh.point(v_tet[2] as usize),
            mesh.point(v_tet[3] as usize),
        ];
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        let mut qpts = Vec::with_capacity(quad.w.len());
        for (q, l) in quad.lam.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            let (vol_q, pos_q) = match &curved {
                Some(ct) => {
                    let gq = safe_curved_jacobian(ct, &lam, &gs);
                    basis3_keyed(asm.order, &lam, &gq, &gid, &mut w, &mut cu);
                    let det_j = ct.jacobian(&lam).1.abs();
                    let vol_c = if det_j.is_finite() && det_j > 0.0 { det_j / 6.0 } else { geo.vol[ci] };
                    (vol_c, ct.position(&lam))
                }
                None => {
                    basis3_keyed(asm.order, &lam, &gs, &gid, &mut w, &mut cu);
                    let mut xq = [0.0; 3];
                    for b in 0..4 {
                        for k in 0..3 {
                            xq[k] += lam[b] * pts_tet[b][k];
                        }
                    }
                    (geo.vol[ci], xq)
                }
            };
            qpts.push(TetQuadPoint3d { wq_vol: quad.w[q] * vol_q, pos: pos_q, w_basis: w[..nb].to_vec() });
        }
        out.push(ActiveCell3d { dom, red_map, qpts });
    }
    out
}

fn extract_domain_nonlinear_maps(
    problem: &Problem,
    asm: &DrivenAssembly,
) -> Result<(BTreeMap<i32, DomainChi2Data>, BTreeMap<i32, DomainChi3Data>), String> {
    let mut dom_n0: BTreeMap<i32, f64> = BTreeMap::new();
    for (ci, eps_c) in asm.eps_cells.iter().enumerate() {
        let dom = problem.materials.domains.keys().copied().find(|&d| d == (ci as i32)).unwrap_or(0);
        let _ = dom;
        let trace_re = ((eps_c[0][0].re + eps_c[1][1].re + eps_c[2][2].re) / 3.0).max(1.0);
        let _ = trace_re;
    }
    let mut chi2_map = BTreeMap::new();
    let mut chi3_map = BTreeMap::new();
    for (&dom, an) in &problem.materials.aniso {
        let Some(ref nl) = an.nonlinear else { continue };
        let n0 = an
            .permittivity
            .as_ref()
            .map(|t| {
                let m = t.matrix();
                ((m[0][0] + m[1][1] + m[2][2]) / 3.0).max(1.0).sqrt()
            })
            .or_else(|| problem.materials.get(dom, "eps_r").map(|e| e.max(1.0).sqrt()))
            .unwrap_or(1.0);
        dom_n0.insert(dom, n0);
        if let Some(ref c2) = nl.chi2 {
            let d_glob = c2.d_global(an.euler_deg);
            let max_d = d_glob.iter().flatten().fold(0.0_f64, |a, &x| a.max(x.abs()));
            if max_d > 0.0 {
                let qpm = nl.qpm.clone().unwrap_or_default();
                qpm.validate()?;
                chi2_map.insert(dom, DomainChi2Data { d_rank3: femcore::material::d_voigt_to_rank3(&d_glob), qpm });
            }
        }
        if let Some(ref c3) = nl.chi3 {
            c3.validate()?;
            if c3.is_active() {
                let (a_mt, b_mt) = c3.isotropic_ab(n0);
                let has_cubic = c3.has_cubic_anisotropy(n0);
                let rank4_glob = c3.rank4_global(n0, an.euler_deg);
                chi3_map.insert(
                    dom,
                    DomainChi3Data { chi3: c3.clone(), n0, euler_deg: an.euler_deg, has_cubic, a_mt, b_mt, rank4_glob },
                );
            }
        }
    }
    Ok((chi2_map, chi3_map))
}

// =============================================================================================
// Pointwise Wirtinger Constitutive Kernels (chi3 + TPA + FCA/FCD and multi-harmonic chi2)
// =============================================================================================

/// Evaluates $\frac{1}{\varepsilon_0}\mathbf{P}_{\text{NL}}(\mathbf{E})$ and its two $3 \times 3$ complex Wirtinger
/// Jacobian blocks $J^{(1)}_{im} = \frac{\partial (P_{\text{NL},i}/\varepsilon_0)}{\partial E_m}$ and
/// $J^{(2)}_{im} = \frac{\partial (P_{\text{NL},i}/\varepsilon_0)}{\partial E_m^*}$ at a single quadrature point.
#[inline]
fn eval_chi3_point_wirtinger(
    d3: &DomainChi3Data,
    omega: f64,
    k0: f64,
    e: &[C; 3],
) -> ([C; 3], [[C; 3]; 3], [[C; 3]; 3], f64) {
    let s = e[0].norm_sqr() + e[1].norm_sqr() + e[2].norm_sqr();
    let mut p_over_eps0 = [c(0.0, 0.0); 3];
    let mut j1 = [[c(0.0, 0.0); 3]; 3];
    let mut j2 = [[c(0.0, 0.0); 3]; 3];

    if !d3.has_cubic {
        let a_mt = d3.a_mt;
        let b_mt = d3.b_mt;
        if a_mt != 0.0 || b_mt != 0.0 {
            let e_dot_e = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
            let ca = a_mt * s;
            let cb = 0.5 * b_mt * e_dot_e;
            for i in 0..3 {
                p_over_eps0[i] += e[i] * ca + cb * e[i].conj();
                j1[i][i] += c(ca, 0.0);
                for m in 0..3 {
                    // d/dE_m [ A (E.E^*) E_i + 0.5 B (E.E) E_i^* ] = A s delta_im + A E_i E_m^* + B E_i^* E_m
                    j1[i][m] += e[i] * e[m].conj() * a_mt + e[i].conj() * e[m] * b_mt;
                    // d/dE_m^* [ A (E.E^*) E_i + 0.5 B (E.E) E_i^* ] = A E_i E_m + 0.5 B (E.E) delta_im
                    j2[i][m] += e[i] * e[m] * a_mt;
                }
                j2[i][i] += cb;
            }
        }
    } else {
        for i in 0..3 {
            let mut acc = c(0.0, 0.0);
            for j in 0..3 {
                for k in 0..3 {
                    let ejk = e[j] * e[k];
                    for l in 0..3 {
                        let chi = d3.rank4_glob[i][j][k][l];
                        if chi != 0.0 {
                            acc += ejk * e[l].conj() * (0.75 * chi);
                            // Symmetry in (j, k) gives factor 2 * 0.75 = 1.5 for d/dE_j
                            j1[i][j] += e[k] * e[l].conj() * (1.5 * chi);
                            j2[i][l] += ejk * (0.75 * chi);
                        }
                    }
                }
            }
            p_over_eps0[i] += acc;
        }
    }

    // Two-photon absorption (TPA) + steady-state free-carrier dispersion/absorption (FCD/FCA):
    // I(s) = 0.5 * c0 * eps0 * n0 * s = kappa_I * s
    // N_c(s) = (tau_c * beta_tpa / (2 * hbar * omega)) * I(s)^2 = c_nc * s^2
    // alpha_tilde(s) * E:
    //   Re(alpha_tilde) = 2 * n0 * k_fcd * N_c(s) = c_fcd * s^2
    //   Im(alpha_tilde) = -(n0 / k0) * (beta_tpa * I(s) + sigma_fca * N_c(s)) = -(c_tpa * s + c_fca * s^2)
    let kappa_i = 0.5 * C0 * EPS0 * d3.n0;
    let beta_tpa = d3.chi3.beta_tpa.unwrap_or(0.0);
    let sigma_fca = d3.chi3.sigma_fca.unwrap_or(0.0);
    let k_fcd = d3.chi3.k_fcd.unwrap_or(0.0);
    let tau_c = d3.chi3.tau_c.unwrap_or(0.0);

    let mut nc_steady = 0.0;
    if beta_tpa != 0.0 || sigma_fca != 0.0 || k_fcd != 0.0 {
        let c_tpa = (d3.n0 / k0.max(1e-30)) * beta_tpa * kappa_i;
        let c_nc = if beta_tpa > 0.0 && tau_c > 0.0 && omega > 0.0 {
            (tau_c * beta_tpa / (2.0 * HBAR * omega)) * kappa_i * kappa_i
        } else {
            0.0
        };
        nc_steady = c_nc * s * s;
        let c_fcd = 2.0 * d3.n0 * k_fcd * c_nc;
        let c_fca = (d3.n0 / k0.max(1e-30)) * sigma_fca * c_nc;

        let alpha_s = c(c_fcd * s * s, -(c_tpa * s + c_fca * s * s));
        let dalpha_ds = c(2.0 * c_fcd * s, -(c_tpa + 2.0 * c_fca * s));
        for i in 0..3 {
            p_over_eps0[i] += alpha_s * e[i];
            j1[i][i] += alpha_s;
            for m in 0..3 {
                j1[i][m] += dalpha_ds * e[i] * e[m].conj();
                j2[i][m] += dalpha_ds * e[i] * e[m];
            }
        }
    }

    (p_over_eps0, j1, j2, nc_steady)
}

// =============================================================================================
// Stage 1 / Stage 2 / Stage 3 of 3D Study::NonlinearDriven
// =============================================================================================

/// Stage 1 of `Study::NonlinearDriven`: harmonic linear `DrivenAssembly` blocks, active nonlinear cell
/// quadrature tables, and static CSR Jacobian sparsity patterns.
pub struct NonlinearDrivenAssembly {
    pub spec: NonlinearDrivenSpec,
    pub resolved_process: NonlinearProcess,
    /// Harmonic frequencies $[f_1]$ (`chi3`), $[f_1, 2 f_1]$ (`shg`), $[f_1, f_2]$ (`chi3_xpm`), or $[f_1, f_2, f_1 + f_2]$ (`sfg`/`dfg`/`opa`).
    pub harmonic_freqs_hz: Vec<f64>,
    pub harmonics: Vec<DrivenAssembly>,
    active_cells: Vec<Vec<ActiveCell3d>>,
    chi2_domains: BTreeMap<i32, DomainChi2Data>,
    chi3_domains: BTreeMap<i32, DomainChi3Data>,
    /// Static $2 N_0 \times 2 N_0$ real CSR pattern for single-tone `chi3` Newton solves.
    pub single_pattern: Option<StaticRealPattern>,
    /// Static $(2 N_0 + 1) \times (2 N_0 + 1)$ bordered real CSR pattern for Keller pseudo-arclength continuation.
    pub bordered_pattern: Option<StaticRealPattern>,
    /// Static coupled real CSR pattern (`4 N` or `6 N`) for multi-harmonic `depleted` Newton solves.
    pub coupled_pattern: Option<StaticRealPattern>,
}

/// Record of one converged operating point along a power ramp or pseudo-arclength continuation curve.
#[derive(Clone, Debug)]
pub struct NonlinearDrivenStep {
    /// Primary input power $P_{\text{in}}$ in Watts.
    pub power_in_w: f64,
    /// Per-harmonic unscaled reduced field vectors `x_red[m]` ($\text{V/m}$ peak phasor).
    pub x_red: Vec<Vec<C>>,
    /// Per-harmonic, per-port incident power $P_{\text{in},p}^{(m)}$ in Watts.
    pub port_p_in_w: Vec<Vec<f64>>,
    /// Per-harmonic, per-port outgoing power $P_{\text{out},p}^{(m)} = |b_{\text{out},p}^{(m)}|^2$ in Watts.
    pub port_p_out_w: Vec<Vec<f64>>,
    /// Per-harmonic, per-port outgoing complex power-wave amplitude $b_{\text{out},p}^{(m)}$ ($\text{W}^{1/2}$).
    pub port_b_out: Vec<Vec<C>>,
    /// Stored electromagnetic electric-mode energy $\frac{1}{2}\int_{\Omega\setminus\text{PML}} \varepsilon_0 \mathbf{E}_1^* \cdot \operatorname{Re}(\boldsymbol{\varepsilon}_r)\mathbf{E}_1\,\mathrm{d}V$ (J) at $\omega_1$.
    pub stored_energy_j: f64,
    /// Peak steady-state free-carrier density $\max_{\mathbf{x}} N_c(\mathbf{x})$ ($\text{m}^{-3}$).
    pub max_carrier_density_m3: f64,
    /// Signed parameter tangent $t_\alpha = \mathrm{d}\alpha / \mathrm{d}s$ along pseudo-arclength continuation (`1.0` for `power_ramp`).
    pub tangent_alpha: f64,
    pub newton_iterations: usize,
    pub residual_norm: f64,
}

/// Stage 2 output of `Study::NonlinearDriven`.
#[derive(Clone, Debug)]
pub struct NonlinearDrivenSolve {
    pub steps: Vec<NonlinearDrivenStep>,
    /// Turning-point input powers $P_{\text{in},*}$ (W) detected along pseudo-arclength continuation where $\dot{\alpha}(s) = 0$.
    pub turning_point_powers_w: Vec<f64>,
    pub total_newton_iterations: usize,
    pub total_factorizations: usize,
    pub symbolic_analyses: usize,
    pub solve_seconds: f64,
}

fn build_chi3_patterns(
    asm: &DrivenAssembly,
    omega: f64,
    cells: &[ActiveCell3d],
    with_bordered: bool,
) -> (StaticRealPattern, Option<StaticRealPattern>) {
    let nred = asm.nred;
    let n_real = 2 * nred;
    let a_lin = asm.assemble_matrix_at(omega);
    let mut complex_pairs: Vec<(usize, usize)> = Vec::new();
    for r in 0..nred {
        complex_pairs.push((r, r));
        let (cols, _) = a_lin.row(r);
        for &col in cols {
            complex_pairs.push((r, col as usize));
        }
    }
    for cell in cells {
        for &(r, _) in &cell.red_map {
            if r == usize::MAX {
                continue;
            }
            for &(s, _) in &cell.red_map {
                if s != usize::MAX {
                    complex_pairs.push((r, s));
                }
            }
        }
    }
    complex_pairs.sort_unstable();
    complex_pairs.dedup();

    let mut real_pairs = Vec::with_capacity(4 * complex_pairs.len());
    for &(r, s) in &complex_pairs {
        for dr in 0..2 {
            for ds in 0..2 {
                real_pairs.push((2 * r + dr, 2 * s + ds));
            }
        }
    }
    let single = StaticRealPattern::from_pairs(n_real, n_real, real_pairs.iter().copied());
    let bordered = if with_bordered {
        let border = n_real;
        for i in 0..n_real {
            real_pairs.push((i, border));
            real_pairs.push((border, i));
        }
        real_pairs.push((border, border));
        Some(StaticRealPattern::from_pairs(n_real + 1, n_real + 1, real_pairs))
    } else {
        None
    };
    (single, bordered)
}

fn build_coupled_pattern(
    harmonics: &[DrivenAssembly],
    omegas: &[f64],
    active_cells: &[Vec<ActiveCell3d>],
) -> StaticRealPattern {
    let m_cnt = harmonics.len();
    let mut offsets = vec![0usize; m_cnt + 1];
    for m in 0..m_cnt {
        offsets[m + 1] = offsets[m] + 2 * harmonics[m].nred;
    }
    let n_tot_real = offsets[m_cnt];
    let mut real_pairs = Vec::new();

    for m in 0..m_cnt {
        let a_lin = harmonics[m].assemble_matrix_at(omegas[m]);
        let base = offsets[m];
        for r in 0..harmonics[m].nred {
            let (cols, _) = a_lin.row(r);
            for &col in cols {
                let s = col as usize;
                for dr in 0..2 {
                    for ds in 0..2 {
                        real_pairs.push((base + 2 * r + dr, base + 2 * s + ds));
                    }
                }
            }
        }
    }
    if !active_cells.is_empty() {
        let n_active = active_cells[0].len();
        for idx in 0..n_active {
            for m1 in 0..m_cnt {
                let base1 = offsets[m1];
                for &(r, _) in &active_cells[m1][idx].red_map {
                    if r == usize::MAX {
                        continue;
                    }
                    for m2 in 0..m_cnt {
                        let base2 = offsets[m2];
                        for &(s, _) in &active_cells[m2][idx].red_map {
                            if s == usize::MAX {
                                continue;
                            }
                            for dr in 0..2 {
                                for ds in 0..2 {
                                    real_pairs.push((base1 + 2 * r + dr, base2 + 2 * s + ds));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    StaticRealPattern::from_pairs(n_tot_real, n_tot_real, real_pairs)
}

/// Assembles Stage 1 (`NonlinearDrivenAssembly`) for a 3D `Study::NonlinearDriven` problem.
pub fn assemble_nonlinear_driven(mesh: &Discretisation, problem: &Problem) -> Result<NonlinearDrivenAssembly, String> {
    if mesh.dim != 3 {
        return Err("the 'nonlinear_driven' study requires a 3D mesh".into());
    }
    let Study::NonlinearDriven(ref spec) = problem.study else {
        return Err("assemble_nonlinear_driven requires Study::NonlinearDriven".into());
    };
    if !(spec.freq > 0.0 && spec.freq.is_finite()) {
        return Err(format!("nonlinear_driven requires positive finite freq (got {})", spec.freq));
    }
    if spec.power_w < 0.0 || spec.power_min_w < 0.0 || spec.power_max_w < 0.0 || spec.signal_power_w < 0.0 {
        return Err("nonlinear_driven optical powers must be non-negative".into());
    }

    // Build harmonic-1 linear DrivenAssembly first to inspect ports and domain nonlinearities.
    let p1 = clone_problem_at_freq(problem, spec.freq)?;
    let asm1 = assemble_emw_driven(mesh, &p1)?;
    if asm1.ports.is_empty() {
        return Err("nonlinear_driven requires at least one 'waveport' or 'lumped_port'".into());
    }
    if spec.drive_port >= asm1.ports.len() {
        return Err(format!(
            "nonlinear_driven drive_port {} is out of range (problem has {} port(s))",
            spec.drive_port,
            asm1.ports.len()
        ));
    }

    let (chi2_domains, chi3_domains) = extract_domain_nonlinear_maps(problem, &asm1)?;
    let resolved_process = match spec.process {
        NonlinearProcess::Auto => {
            if spec.freq_2 > 0.0 && !chi2_domains.is_empty() {
                NonlinearProcess::Sfg
            } else if !chi2_domains.is_empty() && chi3_domains.is_empty() {
                NonlinearProcess::Shg
            } else {
                NonlinearProcess::Chi3
            }
        }
        other => other,
    };

    let harmonic_freqs_hz = match resolved_process {
        NonlinearProcess::Auto | NonlinearProcess::Chi3 => vec![spec.freq],
        NonlinearProcess::Shg => vec![spec.freq, 2.0 * spec.freq],
        NonlinearProcess::Chi3Xpm => {
            if !(spec.freq_2 > 0.0 && spec.freq_2.is_finite()) {
                return Err("nonlinear_driven process 'chi3_xpm' requires positive 'freq_2'".into());
            }
            vec![spec.freq, spec.freq_2]
        }
        NonlinearProcess::Sfg => {
            if !(spec.freq_2 > 0.0 && spec.freq_2.is_finite()) {
                return Err("nonlinear_driven process 'sfg' requires positive 'freq_2'".into());
            }
            vec![spec.freq, spec.freq_2, spec.freq + spec.freq_2]
        }
        NonlinearProcess::Dfg | NonlinearProcess::Opa => {
            // Pump at `freq` = f_3 (highest frequency), signal at `freq_2` = f_1 (< freq), idler at f_2 = freq - freq_2!
            // Or if freq_2 > freq: f_1 = freq, f_2 = freq_2 - freq, f_3 = freq_2.
            if !(spec.freq_2 > 0.0 && spec.freq_2.is_finite()) {
                return Err(format!(
                    "nonlinear_driven process '{}' requires positive signal frequency 'freq_2'",
                    resolved_process.as_str()
                ));
            }
            if (spec.freq - spec.freq_2).abs() <= 1e-9 * spec.freq {
                return Err(
                    "for degenerate 2:1 mixing use process = 'shg' (non-degenerate dfg/opa requires freq != freq_2)"
                        .into(),
                );
            }
            if spec.freq > spec.freq_2 {
                // Ordered (omega_1 = signal, omega_2 = idler = pump - signal, omega_3 = pump) so omega_3 = omega_1 + omega_2!
                vec![spec.freq_2, spec.freq - spec.freq_2, spec.freq]
            } else {
                vec![spec.freq, spec.freq_2 - spec.freq, spec.freq_2]
            }
        }
    };

    let mut harmonics = Vec::with_capacity(harmonic_freqs_hz.len());
    for &fm in &harmonic_freqs_hz {
        if (fm - spec.freq).abs() <= 1e-15 * spec.freq {
            harmonics.push(assemble_emw_driven(mesh, &p1)?);
        } else {
            let pm = clone_problem_at_freq(problem, fm)?;
            harmonics.push(assemble_emw_driven(mesh, &pm)?);
        }
    }
    for h in &harmonics {
        if spec.signal_port >= h.ports.len() {
            return Err(format!(
                "nonlinear_driven signal_port {} is out of range (harmonic has {} port(s))",
                spec.signal_port,
                h.ports.len()
            ));
        }
    }

    // Refine per-domain refractive index n0 from actual harmonic-1 cell permittivity tensors!
    let mut chi3_domains = chi3_domains;
    let mut dom_eps_acc: BTreeMap<i32, (f64, usize)> = BTreeMap::new();
    for ci in 0..mesh.n_cells() {
        if harmonics[0].pml_cells.get(ci).copied().unwrap_or(false) {
            continue;
        }
        let dom = mesh.cell_labels[ci];
        let e = &harmonics[0].eps_cells[ci];
        let tr = ((e[0][0].re + e[1][1].re + e[2][2].re) / 3.0).max(1.0);
        let entry = dom_eps_acc.entry(dom).or_insert((0.0, 0));
        entry.0 += tr;
        entry.1 += 1;
    }
    for (&dom, d3) in &mut chi3_domains {
        if let Some(&(sum_eps, cnt)) = dom_eps_acc.get(&dom) {
            if cnt > 0 {
                let n0 = (sum_eps / cnt as f64).max(1.0).sqrt();
                d3.n0 = n0;
                let (a_mt, b_mt) = d3.chi3.isotropic_ab(n0);
                d3.a_mt = a_mt;
                d3.b_mt = b_mt;
                d3.has_cubic = d3.chi3.has_cubic_anisotropy(n0);
                d3.rank4_glob = d3.chi3.rank4_global(n0, d3.euler_deg);
            }
        }
    }

    let active_dom_set: BTreeSet<i32> = match resolved_process {
        NonlinearProcess::Auto | NonlinearProcess::Chi3 | NonlinearProcess::Chi3Xpm => {
            chi3_domains.keys().copied().collect()
        }
        NonlinearProcess::Shg | NonlinearProcess::Sfg | NonlinearProcess::Dfg | NonlinearProcess::Opa => {
            chi2_domains.keys().copied().collect()
        }
    };

    let active_cells: Vec<Vec<ActiveCell3d>> =
        harmonics.iter().map(|h| build_active_cells_3d(mesh, h, &active_dom_set)).collect();

    let omegas: Vec<f64> = harmonic_freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let (single_pattern, bordered_pattern, coupled_pattern) = if harmonics.len() == 1 {
        let with_bordered = spec.continuation == NonlinearContinuation::Arclength;
        let (s_pat, b_pat) = build_chi3_patterns(&harmonics[0], omegas[0], &active_cells[0], with_bordered);
        (Some(s_pat), b_pat, None)
    } else {
        let c_pat = build_coupled_pattern(&harmonics, &omegas, &active_cells);
        (None, None, Some(c_pat))
    };

    Ok(NonlinearDrivenAssembly {
        spec: spec.clone(),
        resolved_process,
        harmonic_freqs_hz,
        harmonics,
        active_cells,
        chi2_domains,
        chi3_domains,
        single_pattern,
        bordered_pattern,
        coupled_pattern,
    })
}

/// Evaluates the single-tone $\chi^{(3)}$ + TPA + FCA/FCD residual $\hat{\mathbf{F}}(\hat{\mathbf{y}}) \in \mathbb{R}^{2 N_{\text{red}}}$
/// and (optionally) its exact Wirtinger Jacobian values on `asm.single_pattern`, returning peak steady-state carrier density $N_c$ ($\text{m}^{-3}$).
pub fn eval_chi3_residual_and_jacobian(
    asm: &NonlinearDrivenAssembly,
    y_real: &[f64],
    rhs_hat: &[C],
    res_out: &mut [f64],
    mut jac_vals: Option<&mut [f64]>,
) -> f64 {
    let h0 = &asm.harmonics[0];
    let nred = h0.nred;
    let omega = 2.0 * PI * asm.harmonic_freqs_hz[0];
    let k0 = omega / C0;
    let k0_2 = k0 * k0;

    let y_c = join_real_to_complex(y_real);
    let a_lin = h0.assemble_matrix_at(omega);
    let ay = a_lin.matvec(&y_c);
    for r in 0..nred {
        let diff = ay[r] - rhs_hat[r];
        res_out[2 * r] = diff.re;
        res_out[2 * r + 1] = diff.im;
    }

    let pat = asm.single_pattern.as_ref();
    if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
        jv.fill(0.0);
        for r in 0..nred {
            let (cols, vals) = a_lin.row(r);
            for (&col, &val) in cols.iter().zip(vals) {
                p.add_wirtinger_2x2(jv, r, col as usize, val, c(0.0, 0.0));
            }
        }
    }

    let mut max_nc = 0.0_f64;
    for cell in &asm.active_cells[0] {
        let Some(d3) = asm.chi3_domains.get(&cell.dom) else { continue };
        let nb = cell.red_map.len();
        let mut x_loc = vec![c(0.0, 0.0); nb];
        for s in 0..nb {
            let (r, scale) = cell.red_map[s];
            if r != usize::MAX {
                x_loc[s] = scale * y_c[r];
            }
        }
        for qp in &cell.qpts {
            let mut e_pt = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let c_s = x_loc[s];
                if c_s.re != 0.0 || c_s.im != 0.0 {
                    e_pt[0] += c_s * qp.w_basis[s][0];
                    e_pt[1] += c_s * qp.w_basis[s][1];
                    e_pt[2] += c_s * qp.w_basis[s][2];
                }
            }
            let (p_nl, j1_pt, j2_pt, nc_pt) = eval_chi3_point_wirtinger(d3, omega, k0, &e_pt);
            max_nc = max_nc.max(nc_pt);
            let w_k2 = -k0_2 * qp.wq_vol;

            for a in 0..nb {
                let (ra, sc_a) = cell.red_map[a];
                if ra == usize::MAX {
                    continue;
                }
                let wa = qp.w_basis[a];
                let dot_wp = p_nl[0] * wa[0] + p_nl[1] * wa[1] + p_nl[2] * wa[2];
                let contrib = sc_a.conj() * dot_wp * w_k2;
                res_out[2 * ra] += contrib.re;
                res_out[2 * ra + 1] += contrib.im;

                if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
                    let w_j1 = [
                        j1_pt[0][0] * wa[0] + j1_pt[1][0] * wa[1] + j1_pt[2][0] * wa[2],
                        j1_pt[0][1] * wa[0] + j1_pt[1][1] * wa[1] + j1_pt[2][1] * wa[2],
                        j1_pt[0][2] * wa[0] + j1_pt[1][2] * wa[1] + j1_pt[2][2] * wa[2],
                    ];
                    let w_j2 = [
                        j2_pt[0][0] * wa[0] + j2_pt[1][0] * wa[1] + j2_pt[2][0] * wa[2],
                        j2_pt[0][1] * wa[0] + j2_pt[1][1] * wa[1] + j2_pt[2][1] * wa[2],
                        j2_pt[0][2] * wa[0] + j2_pt[1][2] * wa[1] + j2_pt[2][2] * wa[2],
                    ];
                    let sca_conj_wk2 = sc_a.conj() * w_k2;
                    for b in 0..nb {
                        let (rb, sc_b) = cell.red_map[b];
                        if rb == usize::MAX {
                            continue;
                        }
                        let wb = qp.w_basis[b];
                        let m1 = w_j1[0] * wb[0] + w_j1[1] * wb[1] + w_j1[2] * wb[2];
                        let m2 = w_j2[0] * wb[0] + w_j2[1] * wb[1] + w_j2[2] * wb[2];
                        let j1_ab = sca_conj_wk2 * sc_b * m1;
                        let j2_ab = sca_conj_wk2 * sc_b.conj() * m2;
                        p.add_wirtinger_2x2(jv, ra, rb, j1_ab, j2_ab);
                    }
                }
            }
        }
    }
    max_nc
}

/// Builds the reference peak-phasor equilibrated RHS $\hat{\mathbf{b}}_{\text{ref}}$ and port incident power-waves $a_{\text{in},p}$ ($\text{W}^{1/2}$)
/// for a given total scale `power_scale` (so incident port powers equal `base_port_powers * power_scale`).
fn build_harmonic_drive_rhs(
    h: &DrivenAssembly,
    omega: f64,
    base_port_powers_w: &[f64],
    port_phases_rad: &[f64],
    power_scale: f64,
) -> (Vec<C>, Vec<C>, Vec<f64>) {
    let n_ports = h.ports.len();
    let mut rhs_hat = vec![c(0.0, 0.0); h.nred];
    let mut a_in = vec![c(0.0, 0.0); n_ports];
    let mut p_in = vec![0.0; n_ports];
    for p in 0..n_ports {
        let pw = base_port_powers_w.get(p).copied().unwrap_or(0.0).max(0.0) * power_scale.max(0.0);
        let ph = port_phases_rad.get(p).copied().unwrap_or(0.0);
        p_in[p] = pw;
        if pw > 0.0 {
            let phase_c = c(ph.cos(), ph.sin());
            a_in[p] = phase_c * pw.sqrt();
            // Peak phasor has sqrt(2) * RMS amplitude!
            let amp_peak = phase_c * (2.0 * pw).sqrt();
            let (_, gamma_p, _, _) = h.ports[p].excitation_coeffs(omega, &h.waveports, &h.floquet_ports);
            let scale = amp_peak * gamma_p;
            let b_vec = h.ports[p].drive_hat();
            for (ri, &bi) in rhs_hat.iter_mut().zip(&b_vec) {
                *ri += scale * bi;
            }
        }
    }
    (rhs_hat, a_in, p_in)
}

/// Computes outgoing complex power-wave amplitudes $b_{\text{out},p}$ ($\text{W}^{1/2}$) and powers $P_{\text{out},p} = |b_{\text{out},p}|^2$ (W)
/// from an equilibrated solution vector `y_hat` and incident amplitudes `a_in` ($\text{W}^{1/2}$).
fn compute_port_outputs(h: &DrivenAssembly, omega: f64, y_hat: &[C], a_in: &[C]) -> (Vec<C>, Vec<f64>) {
    let n_ports = h.ports.len();
    let mut b_out = Vec::with_capacity(n_ports);
    let mut p_out = Vec::with_capacity(n_ports);
    for p in 0..n_ports {
        let (z0_wave, _, _, _) = h.ports[p].excitation_coeffs(omega, &h.waveports, &h.floquet_ports);
        let v_peak = h.ports[p].recv_dot(y_hat);
        let b_p = v_peak / (2.0 * z0_wave.max(1e-30)).sqrt() - a_in[p];
        p_out.push(b_p.norm_sqr());
        b_out.push(b_p);
    }
    (b_out, p_out)
}

/// Computes time-averaged stored electric-mode energy $U \approx \frac{1}{2}\int_{\Omega\setminus\text{PML}} \varepsilon_0 \mathbf{E}^* \cdot \operatorname{Re}(\boldsymbol{\varepsilon}_r)\mathbf{E}\,\mathrm{d}V$ (Joules)
/// for peak-phasor unscaled reduced field `x_red` ($\text{V/m}$).
fn compute_stored_energy_3d(mesh: &Discretisation, h: &DrivenAssembly, x_red: &[C]) -> f64 {
    let x_full = h.red.expand(x_red);
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    let nb = h.ned.nb();
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];
    let mut u_sum = 0.0_f64;

    for ci in 0..mesh.n_cells() {
        if h.pml_cells.get(ci).copied().unwrap_or(false) {
            continue;
        }
        let cdofs = &h.ned.cell_dofs[ci];
        if cdofs.iter().all(|&gd| x_full[gd].re == 0.0 && x_full[gd].im == 0.0) {
            continue;
        }
        let gid = h.ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for k in 0..3 {
                gs[b][k] = geo.g[ci * 12 + b * 3 + k];
            }
        }
        let eps_c = &h.eps_cells[ci];
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        for (q, l) in quad.lam.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            match &curved {
                Some(ct) => {
                    let gq = safe_curved_jacobian(ct, &lam, &gs);
                    basis3_keyed(h.order, &lam, &gq, &gid, &mut w, &mut cu);
                }
                None => basis3_keyed(h.order, &lam, &gs, &gid, &mut w, &mut cu),
            }
            let mut e_q = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let cf = x_full[cdofs[s]];
                e_q[0] += cf * w[s][0];
                e_q[1] += cf * w[s][1];
                e_q[2] += cf * w[s][2];
            }
            let mut e_eps_e = 0.0;
            for a in 0..3 {
                for b in 0..3 {
                    e_eps_e += eps_c[a][b].re * (e_q[a].conj() * e_q[b]).re;
                }
            }
            // Total stored EM energy of a resonant mode equals 2 * time-averaged electric energy = 2 * (1/4) eps0 int E^* eps E dV = (1/2) eps0 int E^* eps E dV.
            u_sum += 0.5 * EPS0 * e_eps_e * quad.w[q] * geo.vol[ci];
        }
    }
    u_sum * h.symmetry_factor
}

fn build_base_port_powers_w(spec: &NonlinearDrivenSpec, n_ports: usize) -> (Vec<f64>, f64) {
    if !spec.port_powers_w.is_empty() {
        let mut v = vec![0.0; n_ports];
        let mut tot = 0.0;
        for (i, &pw) in spec.port_powers_w.iter().take(n_ports).enumerate() {
            v[i] = pw.max(0.0);
            tot += v[i];
        }
        let p_nom = if tot > 0.0 { tot } else { 1.0 };
        for vi in &mut v {
            *vi /= p_nom;
        }
        (v, if spec.power_w > 0.0 && !spec.powers_w.is_empty() { tot.max(1.0) } else { tot })
    } else {
        let mut v = vec![0.0; n_ports];
        v[spec.drive_port.min(n_ports - 1)] = 1.0;
        (v, spec.power_w)
    }
}

fn solve_chi3_power_ramp(
    mesh: &Discretisation,
    asm: &NonlinearDrivenAssembly,
    t0: Instant,
) -> Result<NonlinearDrivenSolve, String> {
    let spec = &asm.spec;
    let h0 = &asm.harmonics[0];
    let nred = h0.nred;
    let n_real = 2 * nred;
    let omega = 2.0 * PI * asm.harmonic_freqs_hz[0];
    let pat = asm.single_pattern.as_ref().ok_or("missing single_pattern in chi3 solve")?;

    let (unit_port_weights, default_power_w) = build_base_port_powers_w(spec, h0.ports.len());
    let powers: Vec<f64> = if !spec.powers_w.is_empty() {
        spec.powers_w.clone()
    } else if spec.n_steps > 1 {
        let p_hi = if spec.power_max_w > 0.0 { spec.power_max_w } else { default_power_w };
        (0..spec.n_steps)
            .map(|k| spec.power_min_w + (p_hi - spec.power_min_w) * (k as f64) / ((spec.n_steps - 1) as f64))
            .collect()
    } else {
        vec![default_power_w]
    };

    let mut y_prev = vec![0.0_f64; n_real];
    let mut prev_p = 0.0_f64;
    let mut steps = Vec::with_capacity(powers.len());
    let mut tot_iters = 0usize;
    let mut tot_facs = 0usize;
    let mut tot_sym = 0usize;

    // Optional JFNK linear preconditioner factored once at zero field:
    let jfnk_pc = if spec.nl_solver == NonlinearSolverKind::Jfnk {
        let mut dummy_res = vec![0.0; n_real];
        let mut j0_vals = pat.zero_values();
        eval_chi3_residual_and_jacobian(asm, &y_prev, &vec![c(0.0, 0.0); nred], &mut dummy_res, Some(&mut j0_vals));
        let j0_csr = pat.into_csr(j0_vals);
        let mut fac = ShiftedFactorization::new(&j0_csr, None, Definiteness::General).map_err(|e| e.to_string())?;
        fac.set_shift(0.0).map_err(|e| e.to_string())?;
        tot_facs += 1;
        tot_sym += 1;
        Some(FactoredRealPc { fac })
    } else {
        None
    };

    let run_one_power =
        |p_target: f64,
         y_init: &[f64]|
         -> Result<(femcore::linalg::snes::NonlinearResult<f64>, Vec<C>, Vec<C>, Vec<f64>), String> {
            let (rhs_hat, a_in, p_in) =
                build_harmonic_drive_rhs(h0, omega, &unit_port_weights, &spec.port_phases_rad, p_target);
            let rhs_norm = rhs_hat.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            let atol = (1e-14 * rhs_norm).max(1e-30);

            let nl_res = match spec.nl_solver {
                NonlinearSolverKind::Newton => {
                    let mut f_clos = |x: &[f64], out: &mut [f64]| {
                        eval_chi3_residual_and_jacobian(asm, x, &rhs_hat, out, None);
                    };
                    let mut j_clos = |x: &[f64]| -> RCsr {
                        let mut dummy = vec![0.0; n_real];
                        let mut jv = pat.zero_values();
                        eval_chi3_residual_and_jacobian(asm, x, &rhs_hat, &mut dummy, Some(&mut jv));
                        pat.into_csr(jv)
                    };
                    let opts = NewtonOptions {
                        rtol: spec.rtol,
                        atol,
                        stol: 1e-13,
                        max_it: spec.max_it,
                        line_search: LineSearch::default(),
                        jacobian_lag: 1,
                        definiteness: Definiteness::General,
                    };
                    newton(&mut f_clos, &mut j_clos, y_init, &opts).map_err(|e| format!("chi3 Newton failed: {e}"))?
                }
                NonlinearSolverKind::Jfnk => {
                    let mut f_clos = |x: &[f64], out: &mut [f64]| {
                        eval_chi3_residual_and_jacobian(asm, x, &rhs_hat, out, None);
                    };
                    let opts = JfnkOptions {
                        rtol: spec.rtol,
                        atol,
                        max_it: spec.max_it,
                        line_search: LineSearch::default(),
                        ..Default::default()
                    };
                    jfnk(&mut f_clos, y_init, jfnk_pc.as_ref().map(|p| p as &dyn Preconditioner<f64>), &opts)
                }
            };
            Ok((nl_res, rhs_hat, a_in, p_in))
        };

    for &pw in &powers {
        femcore::progress::check()?;
        // Adaptive homotopy sub-stepping when a large jump into stiff FCA/FCD saturation exceeds Newton's single-step basin:
        let mut pending = vec![pw];
        let mut bisections = 0usize;
        let mut final_iters = 0usize;
        let mut final_rnorm = 0.0_f64;
        let mut rhs_hat = Vec::new();
        let mut a_in = Vec::new();
        let mut p_in = Vec::new();

        while let Some(p_try) = pending.pop() {
            femcore::progress::check()?;
            let mut y_init = if prev_p > 0.0 && p_try > 0.0 {
                let ratio = (p_try / prev_p).sqrt();
                y_prev.iter().map(|&v| v * ratio).collect::<Vec<f64>>()
            } else {
                y_prev.clone()
            };
            if p_try == 0.0 {
                y_init.fill(0.0);
            }

            let (mut nl_res, mut rh, mut ai, mut pi) = run_one_power(p_try, &y_init)?;
            if !nl_res.converged() && prev_p > 0.0 && p_try > prev_p {
                // Under strong TPA/FCA clamping, unscaled y_prev is often closer than sqrt(p_try/prev_p) * y_prev:
                let (nl_alt, rh_alt, ai_alt, pi_alt) = run_one_power(p_try, &y_prev)?;
                tot_iters += nl_res.iterations;
                tot_facs += nl_res.factorizations;
                tot_sym += nl_res.symbolic_analyses;
                nl_res = nl_alt;
                rh = rh_alt;
                ai = ai_alt;
                pi = pi_alt;
            }

            if nl_res.converged() {
                tot_iters += nl_res.iterations;
                tot_facs += nl_res.factorizations;
                tot_sym += nl_res.symbolic_analyses;
                y_prev.copy_from_slice(&nl_res.x);
                prev_p = p_try;
                if pending.is_empty() {
                    final_iters = nl_res.iterations;
                    final_rnorm = *nl_res.residual_norms.last().unwrap_or(&0.0);
                    rhs_hat = rh;
                    a_in = ai;
                    p_in = pi;
                }
            } else if bisections < 16 && (p_try - prev_p).abs() > 1e-12 * p_try.abs().max(1e-30) {
                tot_iters += nl_res.iterations;
                tot_facs += nl_res.factorizations;
                tot_sym += nl_res.symbolic_analyses;
                let p_mid = 0.5 * (prev_p + p_try);
                pending.push(p_try);
                pending.push(p_mid);
                bisections += 1;
            } else {
                return Err(format!(
                    "chi3 nonlinear solver did not converge at P_in = {p_try:.4e} W: {:?} (residual = {:.3e})",
                    nl_res.reason,
                    nl_res.residual_norms.last().copied().unwrap_or(f64::NAN)
                ));
            }
        }

        let mut res_chk = vec![0.0; n_real];
        let max_nc = eval_chi3_residual_and_jacobian(asm, &y_prev, &rhs_hat, &mut res_chk, None);
        let y_hat_c = join_real_to_complex(&y_prev);
        let x_red: Vec<C> = y_hat_c.iter().zip(&h0.d_ruiz).map(|(&y, &d)| y * d).collect();
        let (b_out, p_out) = compute_port_outputs(h0, omega, &y_hat_c, &a_in);
        let u_j = compute_stored_energy_3d(mesh, h0, &x_red);

        steps.push(NonlinearDrivenStep {
            power_in_w: pw,
            x_red: vec![x_red],
            port_p_in_w: vec![p_in],
            port_p_out_w: vec![p_out],
            port_b_out: vec![b_out],
            stored_energy_j: u_j,
            max_carrier_density_m3: max_nc,
            tangent_alpha: 1.0,
            newton_iterations: final_iters,
            residual_norm: final_rnorm,
        });
    }

    Ok(NonlinearDrivenSolve {
        steps,
        turning_point_powers_w: Vec::new(),
        total_newton_iterations: tot_iters,
        total_factorizations: tot_facs,
        symbolic_analyses: tot_sym,
        solve_seconds: t0.elapsed().as_secs_f64(),
    })
}

/// Solves single-tone $\chi^{(3)}$ Harmonic Balance FEM along a Keller **pseudo-arclength continuation** path
/// parameterized by $\alpha = \sqrt{P_{\text{in}} / P_{\text{ref}}} \ge 0$, using a nonsingular $(2N_{\text{red}} + 1)$
/// bordered sparse LU factorization to traverse upper and lower optical-bistability turning points cleanly.
fn solve_chi3_arclength(
    mesh: &Discretisation,
    asm: &NonlinearDrivenAssembly,
    t0: Instant,
) -> Result<NonlinearDrivenSolve, String> {
    let spec = &asm.spec;
    let h0 = &asm.harmonics[0];
    let nred = h0.nred;
    let n_real = 2 * nred;
    let border = n_real;
    let omega = 2.0 * PI * asm.harmonic_freqs_hz[0];
    let s_pat = asm.single_pattern.as_ref().ok_or("missing single_pattern")?;
    let b_pat = asm.bordered_pattern.as_ref().ok_or("missing bordered_pattern")?;

    let (unit_port_weights, default_power_w) = build_base_port_powers_w(spec, h0.ports.len());
    let p_max = if spec.power_max_w > 0.0 {
        spec.power_max_w
    } else if !spec.powers_w.is_empty() {
        spec.powers_w.iter().copied().fold(0.0, f64::max)
    } else {
        default_power_w.max(1e-6)
    };
    let p_ref = p_max.max(1e-30);
    let alpha_stop = 1.0_f64; // alpha = sqrt(P / P_ref), so alpha = 1 corresponds to P = p_max!

    let (rhs_ref_hat, a_in_ref, p_in_ref) =
        build_harmonic_drive_rhs(h0, omega, &unit_port_weights, &spec.port_phases_rad, p_ref);
    let mut f_alpha = vec![0.0; n_real]; // dF / d(alpha) = -b_ref_real
    for r in 0..nred {
        f_alpha[2 * r] = -rhs_ref_hat[r].re;
        f_alpha[2 * r + 1] = -rhs_ref_hat[r].im;
    }
    let rhs_ref_norm = f_alpha.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);

    // 1. Compute linear response v0 = J(0)^(-1) * b_ref_real to calibrate balanced Keller field scale y_scale = ||v0||_2:
    let mut y_cur = vec![0.0; n_real];
    let mut alpha_cur = (spec.power_min_w.max(0.0) / p_ref).sqrt();
    let mut res_buf = vec![0.0; n_real];
    let mut j0_single = s_pat.zero_values();
    eval_chi3_residual_and_jacobian(asm, &y_cur, &vec![c(0.0, 0.0); nred], &mut res_buf, Some(&mut j0_single));
    let j0_row_scale = (j0_single.iter().map(|v| v * v).sum::<f64>() / (n_real.max(1) as f64)).sqrt().max(1.0);

    // Equilibrate column `border` by 1 / rhs_ref_norm (unknown = rhs_ref_norm * d_alpha) and row `border` by v0_norm
    // so every block of the bordered matrix is O(1) even inside ultra-high-Q (Q >= 1e5) resonators!
    let inv_rhs_norm = 1.0 / rhs_ref_norm;
    let f_alpha_unit: Vec<f64> = f_alpha.iter().map(|&v| v * inv_rhs_norm).collect();
    let assemble_bordered = |j_single: &[f64], wy_scaled: &[f64], w_alpha_col: f64| -> RCsr {
        let mut b_vals = b_pat.zero_values();
        for r in 0..n_real {
            let s0 = s_pat.indptr[r];
            let s1 = s_pat.indptr[r + 1];
            for p in s0..s1 {
                b_pat.add_entry(&mut b_vals, r, s_pat.indices[p], j_single[p]);
            }
            b_pat.add_entry(&mut b_vals, r, border, f_alpha_unit[r]);
            b_pat.add_entry(&mut b_vals, border, r, wy_scaled[r]);
        }
        b_pat.add_entry(&mut b_vals, border, border, w_alpha_col);
        b_pat.into_csr(b_vals)
    };

    // Factor bordered matrix at (0, 0) with wy = 0, w_alpha_col = 1 => solves J(0) v0 + f_alpha_unit * (rhs_ref_norm * 1) = 0:
    let zero_wy = vec![0.0; n_real];
    let b0_csr = assemble_bordered(&j0_single, &zero_wy, 1.0);
    let mut fac = ShiftedFactorization::new(&b0_csr, None, Definiteness::General).map_err(|e| e.to_string())?;
    fac.set_shift(0.0).map_err(|e| e.to_string())?;
    let mut tang_rhs = vec![0.0; n_real + 1];
    tang_rhs[border] = rhs_ref_norm;
    fac.solve_block_in_place(&mut tang_rhs, 1).map_err(|e| e.to_string())?;
    let v0_norm = tang_rhs[..n_real].iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-30);
    let inv_yscale2 = 1.0 / (v0_norm * v0_norm);

    let normalize_w = |ty: &mut [f64], ta: &mut f64| {
        let ny2: f64 = ty.iter().map(|v| v * v).sum();
        let nrm = (ny2 * inv_yscale2 + (*ta) * (*ta)).sqrt().max(1e-30);
        for v in ty.iter_mut() {
            *v /= nrm;
        }
        *ta /= nrm;
    };

    let mut t_y: Vec<f64> = tang_rhs[..n_real].to_vec();
    let mut t_alpha = tang_rhs[border] * inv_rhs_norm;
    normalize_w(&mut t_y, &mut t_alpha);

    // If power_min_w > 0, converge the starting point at alpha_cur first!
    if alpha_cur > 0.0 {
        for i in 0..n_real {
            y_cur[i] = alpha_cur * tang_rhs[i];
        }
        let rhs_start: Vec<C> = rhs_ref_hat.iter().map(|&z| z * alpha_cur).collect();
        let mut f_clos = |x: &[f64], out: &mut [f64]| {
            eval_chi3_residual_and_jacobian(asm, x, &rhs_start, out, None);
        };
        let mut j_clos = |x: &[f64]| -> RCsr {
            let mut dummy = vec![0.0; n_real];
            let mut jv = s_pat.zero_values();
            eval_chi3_residual_and_jacobian(asm, x, &rhs_start, &mut dummy, Some(&mut jv));
            s_pat.into_csr(jv)
        };
        let start_res = newton(
            &mut f_clos,
            &mut j_clos,
            &y_cur,
            &NewtonOptions { rtol: spec.rtol, atol: 1e-13 * rhs_ref_norm, max_it: spec.max_it, ..Default::default() },
        )
        .map_err(|e| e.to_string())?;
        if !start_res.converged() {
            return Err("arclength start point failed to converge".into());
        }
        y_cur = start_res.x;
    }

    let make_step_record = |y_vec: &[f64], alpha_val: f64, ta: f64, iters: usize, rnorm: f64| -> NonlinearDrivenStep {
        let pw = alpha_val * alpha_val * p_ref;
        let rhs_cur: Vec<C> = rhs_ref_hat.iter().map(|&z| z * alpha_val).collect();
        let a_in_cur: Vec<C> = a_in_ref.iter().map(|&z| z * alpha_val).collect();
        let p_in_cur: Vec<f64> = p_in_ref.iter().map(|&p| p * alpha_val * alpha_val).collect();
        let mut rchk = vec![0.0; n_real];
        let max_nc = eval_chi3_residual_and_jacobian(asm, y_vec, &rhs_cur, &mut rchk, None);
        let y_hat_c = join_real_to_complex(y_vec);
        let x_red: Vec<C> = y_hat_c.iter().zip(&h0.d_ruiz).map(|(&y, &d)| y * d).collect();
        let (b_out, p_out) = compute_port_outputs(h0, omega, &y_hat_c, &a_in_cur);
        let u_j = compute_stored_energy_3d(mesh, h0, &x_red);
        NonlinearDrivenStep {
            power_in_w: pw,
            x_red: vec![x_red],
            port_p_in_w: vec![p_in_cur],
            port_p_out_w: vec![p_out],
            port_b_out: vec![b_out],
            stored_energy_j: u_j,
            max_carrier_density_m3: max_nc,
            tangent_alpha: ta,
            newton_iterations: iters,
            residual_norm: rnorm,
        }
    };

    let mut steps = Vec::new();
    steps.push(make_step_record(&y_cur, alpha_cur, t_alpha, 0, 0.0));

    let mut ds = spec.arclength_ds.clamp(1e-4, 0.5);
    let ds_min = 1e-5;
    let ds_max = (2.5 * spec.arclength_ds).clamp(0.02, 0.35);
    let mut s_cur = 0.0_f64;
    let mut turning_points_w = Vec::new();
    let mut tot_iters = 0usize;

    for _step_idx in 0..spec.arclength_max_steps {
        femcore::progress::check()?;
        let wy: Vec<f64> = t_y.iter().map(|&v| v * inv_yscale2).collect();
        let wy_scaled: Vec<f64> = t_y.iter().map(|&v| v / v0_norm).collect();
        let w_alpha = t_alpha;
        let w_alpha_col = v0_norm * t_alpha * inv_rhs_norm;

        let mut accepted = None;
        for _retry in 0..8 {
            // Predictor:
            let mut y_m: Vec<f64> = y_cur.iter().zip(&t_y).map(|(&y, &ty)| y + ds * ty).collect();
            let mut alpha_m = alpha_cur + ds * t_alpha;
            if alpha_m < 0.0 {
                ds *= 0.5;
                if ds < ds_min {
                    break;
                }
                continue;
            }

            // Bordered Newton corrector on (y, alpha):
            let mut converged = false;
            let mut it_used = 0usize;
            let mut final_res = f64::INFINITY;
            let mut j_single = s_pat.zero_values();

            for it in 0..spec.max_it {
                let rhs_m: Vec<C> = rhs_ref_hat.iter().map(|&z| z * alpha_m).collect();
                eval_chi3_residual_and_jacobian(asm, &y_m, &rhs_m, &mut res_buf, Some(&mut j_single));
                let mut n_constraint = -ds;
                for i in 0..n_real {
                    n_constraint += wy[i] * (y_m[i] - y_cur[i]);
                }
                n_constraint += w_alpha * (alpha_m - alpha_cur);

                let f_norm = res_buf.iter().map(|v| v * v).sum::<f64>().sqrt();
                let y_norm = y_m.iter().map(|v| v * v).sum::<f64>().sqrt();
                let op_scale = rhs_ref_norm + j0_row_scale * y_norm;
                let f_tol = (spec.rtol * rhs_ref_norm).max(5e-13 * op_scale);
                let n_tol = spec.rtol.max(1e-12) * ds.max(1.0);
                final_res = f_norm / op_scale.max(1e-30);
                if f_norm <= f_tol && n_constraint.abs() <= n_tol {
                    converged = true;
                    it_used = it;
                    break;
                }

                let b_csr = assemble_bordered(&j_single, &wy_scaled, w_alpha_col);
                if fac.replace_matrices(&b_csr, None).is_err() || fac.set_shift(0.0).is_err() {
                    break;
                }
                let mut delta = Vec::with_capacity(n_real + 1);
                for &fi in &res_buf {
                    delta.push(-fi);
                }
                delta.push(-v0_norm * n_constraint);
                if fac.solve_block_in_place(&mut delta, 1).is_err() || delta.iter().any(|v| !v.is_finite()) {
                    break;
                }
                let d_alpha = delta[border] * inv_rhs_norm;

                // Damped step if alpha_m + d_alpha would cross below 0:
                let mut damp = 1.0_f64;
                if alpha_m + d_alpha < 0.0 {
                    damp = (0.8 * alpha_m / (-d_alpha).max(1e-30)).clamp(0.05, 1.0);
                }
                for i in 0..n_real {
                    y_m[i] += damp * delta[i];
                }
                alpha_m += damp * d_alpha;
                it_used = it + 1;
            }

            if converged && alpha_m >= 0.0 {
                // Solve for next unit tangent on the converged bordered Jacobian!
                let rhs_m: Vec<C> = rhs_ref_hat.iter().map(|&z| z * alpha_m).collect();
                eval_chi3_residual_and_jacobian(asm, &y_m, &rhs_m, &mut res_buf, Some(&mut j_single));
                let b_csr = assemble_bordered(&j_single, &wy_scaled, w_alpha_col);
                if fac.replace_matrices(&b_csr, None).is_ok() && fac.set_shift(0.0).is_ok() {
                    let mut next_t = vec![0.0; n_real + 1];
                    next_t[border] = v0_norm;
                    if fac.solve_block_in_place(&mut next_t, 1).is_ok() && next_t.iter().all(|v| v.is_finite()) {
                        let mut nt_y = next_t[..n_real].to_vec();
                        let mut nt_a = next_t[border] * inv_rhs_norm;
                        normalize_w(&mut nt_y, &mut nt_a);
                        accepted = Some((y_m, alpha_m, nt_y, nt_a, it_used, final_res));
                        break;
                    }
                }
            }

            ds *= 0.5;
            if ds < ds_min {
                break;
            }
        }

        let Some((y_next, alpha_next, nt_y, nt_a, it_used, final_res)) = accepted else {
            return Err(format!(
                "pseudo-arclength continuation stalled near P_in = {:.4e} W (ds dropped below {ds_min:.1e})",
                alpha_cur * alpha_cur * p_ref
            ));
        };

        tot_iters += it_used;
        let ds_used = ds;

        // Detect saddle-node turning point when t_alpha changes sign between s_cur and s_cur + ds_used!
        if (t_alpha > 0.0 && nt_a <= 0.0) || (t_alpha < 0.0 && nt_a >= 0.0) {
            // Cubic Hermite interpolation of alpha(s) for s = u * ds_used, u in [0, 1]:
            // alpha(u) = (2u^3 - 3u^2 + 1) a0 + (u^3 - 2u^2 + u) (ds * ta0) + (-2u^3 + 3u^2) a1 + (u^3 - u^2) (ds * ta1)
            // Find root of alpha'(u) = 0 in [0, 1] by bisection (monotone derivative near simple fold):
            let a0 = alpha_cur;
            let a1 = alpha_next;
            let m0 = ds_used * t_alpha;
            let m1 = ds_used * nt_a;
            let d_hermite = |u: f64| -> f64 {
                (6.0 * u * u - 6.0 * u) * (a0 - a1) + (3.0 * u * u - 4.0 * u + 1.0) * m0 + (3.0 * u * u - 2.0 * u) * m1
            };
            let eval_hermite = |u: f64| -> f64 {
                let u2 = u * u;
                let u3 = u2 * u;
                (2.0 * u3 - 3.0 * u2 + 1.0) * a0
                    + (u3 - 2.0 * u2 + u) * m0
                    + (-2.0 * u3 + 3.0 * u2) * a1
                    + (u3 - u2) * m1
            };
            let mut lo = 0.0_f64;
            let mut hi = 1.0_f64;
            let d_lo = d_hermite(lo);
            for _ in 0..50 {
                let mid = 0.5 * (lo + hi);
                let dm = d_hermite(mid);
                if (d_lo >= 0.0 && dm >= 0.0) || (d_lo <= 0.0 && dm <= 0.0) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let alpha_star = eval_hermite(0.5 * (lo + hi)).max(0.0);
            turning_points_w.push(alpha_star * alpha_star * p_ref);
        }

        s_cur += ds_used;
        let _ = s_cur;
        y_cur = y_next;
        alpha_cur = alpha_next;
        t_y = nt_y;
        t_alpha = nt_a;

        steps.push(make_step_record(&y_cur, alpha_cur, t_alpha, it_used, final_res));

        // Stop cleanly once we are on the forward upper branch (t_alpha > 0) and have reached P_max (alpha >= 1.0):
        if alpha_cur >= alpha_stop && t_alpha > 0.0 {
            break;
        }

        // Gentle step adaptation: hold step steady near folds (|t_alpha| <= 0.35) and grow mildly on smooth branches:
        if it_used <= 3 && t_alpha.abs() > 0.35 {
            ds = (ds * 1.15).min(ds_max);
        } else if it_used >= 7 {
            ds = (ds * 0.75).max(ds_min);
        }
    }

    let stats = fac.stats();
    Ok(NonlinearDrivenSolve {
        steps,
        turning_point_powers_w: turning_points_w,
        total_newton_iterations: tot_iters,
        total_factorizations: stats.numeric_count,
        symbolic_analyses: stats.symbolic_count,
        solve_seconds: t0.elapsed().as_secs_f64(),
    })
}

// =============================================================================================
// Multi-Harmonic Chi^(2) (SHG / SFG / DFG / OPA) & Two-Tone Chi^(3) XPM Residual + Jacobian
// =============================================================================================

/// Evaluates the coupled real residual and (optionally) exact Wirtinger block Jacobian on `asm.coupled_pattern`
/// for multi-harmonic `shg` ($2$ harmonics), `chi3_xpm` ($2$ harmonics), or `sfg`/`dfg`/`opa` ($3$ harmonics).
pub fn eval_chi2_residual_and_jacobian(
    asm: &NonlinearDrivenAssembly,
    y_real: &[f64],
    rhs_hats: &[Vec<C>],
    include_backaction: bool,
    res_out: &mut [f64],
    mut jac_vals: Option<&mut [f64]>,
) {
    let m_cnt = asm.harmonics.len();
    let mut offsets = vec![0usize; m_cnt + 1];
    let mut c_offsets = vec![0usize; m_cnt + 1];
    for m in 0..m_cnt {
        c_offsets[m + 1] = c_offsets[m] + asm.harmonics[m].nred;
        offsets[m + 1] = offsets[m] + 2 * asm.harmonics[m].nred;
    }

    let omegas: Vec<f64> = asm.harmonic_freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let k0_2: Vec<f64> = omegas.iter().map(|&w| (w / C0).powi(2)).collect();

    let y_c: Vec<Vec<C>> = (0..m_cnt).map(|m| join_real_to_complex(&y_real[offsets[m]..offsets[m + 1]])).collect();

    let pat = asm.coupled_pattern.as_ref();
    if let Some(jv) = jac_vals.as_deref_mut() {
        jv.fill(0.0);
    }

    // Linear curl-curl blocks A_m(omega_m) y_m - rhs_m:
    for m in 0..m_cnt {
        let h = &asm.harmonics[m];
        let a_lin = h.assemble_matrix_at(omegas[m]);
        let ay = a_lin.matvec(&y_c[m]);
        let base_r = offsets[m];
        let base_c = c_offsets[m];
        for r in 0..h.nred {
            let diff = ay[r] - rhs_hats[m][r];
            res_out[base_r + 2 * r] = diff.re;
            res_out[base_r + 2 * r + 1] = diff.im;
        }
        if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
            for r in 0..h.nred {
                let (cols, vals) = a_lin.row(r);
                for (&col, &val) in cols.iter().zip(vals) {
                    p.add_wirtinger_2x2(jv, base_c + r, base_c + col as usize, val, c(0.0, 0.0));
                }
            }
        }
    }

    let n_active = asm.active_cells.first().map_or(0, Vec::len);
    for idx in 0..n_active {
        let dom = asm.active_cells[0][idx].dom;

        if asm.resolved_process == NonlinearProcess::Chi3Xpm {
            let Some(d3) = asm.chi3_domains.get(&dom) else { continue };
            let c1 = &asm.active_cells[0][idx];
            let c2 = &asm.active_cells[1][idx];
            let nb = c1.red_map.len();
            let mut x1 = vec![c(0.0, 0.0); nb];
            let mut x2 = vec![c(0.0, 0.0); nb];
            for s in 0..nb {
                let (r1, sc1) = c1.red_map[s];
                if r1 != usize::MAX {
                    x1[s] = sc1 * y_c[0][r1];
                }
                let (r2, sc2) = c2.red_map[s];
                if r2 != usize::MAX {
                    x2[s] = sc2 * y_c[1][r2];
                }
            }
            let a_mt = d3.a_mt;
            let b_mt = d3.b_mt;
            let kappa_i = 0.5 * C0 * EPS0 * d3.n0;
            let beta_tpa = d3.chi3.beta_tpa.unwrap_or(0.0);

            for q_idx in 0..c1.qpts.len() {
                let qp = &c1.qpts[q_idx];
                let mut e1 = [c(0.0, 0.0); 3];
                let mut e2 = [c(0.0, 0.0); 3];
                for s in 0..nb {
                    for k in 0..3 {
                        e1[k] += x1[s] * qp.w_basis[s][k];
                        e2[k] += x2[s] * qp.w_basis[s][k];
                    }
                }
                // SPM + Maker-Terhune XPM at (omega_1, omega_2):
                // P_NL(omega_2)/eps0 = SPM(E2) + A [(E1.E1^*) E2 + (E1^*.E2) E1] + B (E1.E2) E1^* - i (2 n0 / k0_2) beta_tpa I_1 E_2
                let (mut p1, mut j1_11, j2_11, _) = eval_chi3_point_wirtinger(d3, omegas[0], omegas[0] / C0, &e1);
                let (mut p2, mut j1_22, j2_22, _) = eval_chi3_point_wirtinger(d3, omegas[1], omegas[1] / C0, &e2);

                let e1_n2 = e1[0].norm_sqr() + e1[1].norm_sqr() + e1[2].norm_sqr();
                let e2_n2 = e2[0].norm_sqr() + e2[1].norm_sqr() + e2[2].norm_sqr();
                let e1c_e2 = e1[0].conj() * e2[0] + e1[1].conj() * e2[1] + e1[2].conj() * e2[2];
                let e1_e2 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
                let x_tpa_2 = -2.0 * (d3.n0 / (omegas[1] / C0).max(1e-30)) * beta_tpa * kappa_i;
                let x_tpa_1 = -2.0 * (d3.n0 / (omegas[0] / C0).max(1e-30)) * beta_tpa * kappa_i;

                let mut j1_21 = [[c(0.0, 0.0); 3]; 3];
                let mut j2_21 = [[c(0.0, 0.0); 3]; 3];
                let mut j1_12 = [[c(0.0, 0.0); 3]; 3];
                let mut j2_12 = [[c(0.0, 0.0); 3]; 3];

                for i in 0..3 {
                    let xpm_2 = e2[i] * c(a_mt * e1_n2, x_tpa_2 * e1_n2)
                        + e1[i] * (a_mt * e1c_e2)
                        + e1[i].conj() * (b_mt * e1_e2);
                    p2[i] += xpm_2;
                    j1_22[i][i] += c(a_mt * e1_n2, x_tpa_2 * e1_n2);
                    for m in 0..3 {
                        j1_22[i][m] += e1[i] * e1[m].conj() * a_mt + e1[i].conj() * e1[m] * b_mt;
                        j1_21[i][m] += e2[i] * e1[m].conj() * c(a_mt, x_tpa_2) + e1[i].conj() * e2[m] * b_mt;
                        if i == m {
                            j1_21[i][m] += e1c_e2 * a_mt;
                        }
                        j2_21[i][m] += e2[i] * e1[m] * c(a_mt, x_tpa_2) + e1[i] * e2[m] * a_mt;
                        if i == m {
                            j2_21[i][m] += e1_e2 * b_mt;
                        }
                    }

                    if include_backaction {
                        let xpm_1 = e1[i] * c(a_mt * e2_n2, x_tpa_1 * e2_n2)
                            + e2[i] * (a_mt * e1c_e2.conj())
                            + e2[i].conj() * (b_mt * e1_e2);
                        p1[i] += xpm_1;
                        j1_11[i][i] += c(a_mt * e2_n2, x_tpa_1 * e2_n2);
                        for m in 0..3 {
                            j1_11[i][m] += e2[i] * e2[m].conj() * a_mt + e2[i].conj() * e2[m] * b_mt;
                            j1_12[i][m] += e1[i] * e2[m].conj() * c(a_mt, x_tpa_1) + e2[i].conj() * e1[m] * b_mt;
                            if i == m {
                                j1_12[i][m] += e1c_e2.conj() * a_mt;
                            }
                            j2_12[i][m] += e1[i] * e2[m] * c(a_mt, x_tpa_1) + e2[i] * e1[m] * a_mt;
                            if i == m {
                                j2_12[i][m] += e1_e2 * b_mt;
                            }
                        }
                    }
                }

                let p_arr = [p1, p2];
                let cells_pair = [c1, c2];
                for m_row in 0..2 {
                    let w_k2 = -k0_2[m_row] * qp.wq_vol;
                    let base_r = offsets[m_row];
                    for a in 0..nb {
                        let (ra, sc_a) = cells_pair[m_row].red_map[a];
                        if ra == usize::MAX {
                            continue;
                        }
                        let wa = qp.w_basis[a];
                        let dot_wp = p_arr[m_row][0] * wa[0] + p_arr[m_row][1] * wa[1] + p_arr[m_row][2] * wa[2];
                        let contrib = sc_a.conj() * dot_wp * w_k2;
                        res_out[base_r + 2 * ra] += contrib.re;
                        res_out[base_r + 2 * ra + 1] += contrib.im;
                    }
                }

                if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
                    let blocks = [
                        (0usize, 0usize, &j1_11, &j2_11),
                        (0usize, 1usize, &j1_12, &j2_12),
                        (1usize, 0usize, &j1_21, &j2_21),
                        (1usize, 1usize, &j1_22, &j2_22),
                    ];
                    for (mr, mc, j1_b, j2_b) in blocks {
                        if !include_backaction && mr == 0 && mc == 1 {
                            continue;
                        }
                        let w_k2 = -k0_2[mr] * qp.wq_vol;
                        let bc_r = c_offsets[mr];
                        let bc_c = c_offsets[mc];
                        for a in 0..nb {
                            let (ra, sc_a) = cells_pair[mr].red_map[a];
                            if ra == usize::MAX {
                                continue;
                            }
                            let wa = qp.w_basis[a];
                            let sca_w = sc_a.conj() * w_k2;
                            for b in 0..nb {
                                let (rb, sc_b) = cells_pair[mc].red_map[b];
                                if rb == usize::MAX {
                                    continue;
                                }
                                let wb = qp.w_basis[b];
                                let (mut m1, mut m2) = (c(0.0, 0.0), c(0.0, 0.0));
                                for ii in 0..3 {
                                    for mm in 0..3 {
                                        let wab = wa[ii] * wb[mm];
                                        m1 += j1_b[ii][mm] * wab;
                                        m2 += j2_b[ii][mm] * wab;
                                    }
                                }
                                p.add_wirtinger_2x2(
                                    jv,
                                    bc_r + ra,
                                    bc_c + rb,
                                    sca_w * sc_b * m1,
                                    sca_w * sc_b.conj() * m2,
                                );
                            }
                        }
                    }
                }
            }
            continue;
        }

        // Second-order chi^(2): SHG (m_cnt == 2) or 3-wave SFG/DFG/OPA (m_cnt == 3):
        let Some(d2) = asm.chi2_domains.get(&dom) else { continue };
        let d_r3 = &d2.d_rank3;
        let nb = asm.active_cells[0][idx].red_map.len();

        if m_cnt == 2 {
            // Degenerate SHG (omega_2 = 2 omega_1):
            // P(omega_2)/eps0 = s_q * sum_{j,k} d_{ijk} E_{1,j} E_{1,k}
            // P(omega_1)/eps0 = 2 * s_q * sum_{i,k} d_{ijk} E_{2,i} E_{1,k}^* (when include_backaction)
            let c1 = &asm.active_cells[0][idx];
            let c2 = &asm.active_cells[1][idx];
            let mut x1 = vec![c(0.0, 0.0); nb];
            let mut x2 = vec![c(0.0, 0.0); nb];
            for s in 0..nb {
                let (r1, sc1) = c1.red_map[s];
                if r1 != usize::MAX {
                    x1[s] = sc1 * y_c[0][r1];
                }
                let (r2, sc2) = c2.red_map[s];
                if r2 != usize::MAX {
                    x2[s] = sc2 * y_c[1][r2];
                }
            }
            for q_idx in 0..c1.qpts.len() {
                let qp = &c1.qpts[q_idx];
                let sq = d2.qpm.poling_sign_at(qp.pos);
                let mut e1 = [c(0.0, 0.0); 3];
                let mut e2 = [c(0.0, 0.0); 3];
                for s in 0..nb {
                    for k in 0..3 {
                        e1[k] += x1[s] * qp.w_basis[s][k];
                        e2[k] += x2[s] * qp.w_basis[s][k];
                    }
                }
                let mut p_w2 = [c(0.0, 0.0); 3];
                let mut p_w1 = [c(0.0, 0.0); 3];
                // Local 3x3 Wirtinger derivatives:
                // d(P2_i/eps0)/dE1_m = 2 s_q sum_k d_{imk} E1_k  (holomorphic J1_{2,1})
                // d(P1_j/eps0)/dE2_m = 2 s_q sum_k d_{mjk} E1_k^* (holomorphic J1_{1,2})
                // d(P1_j/eps0)/dE1_m^* = 2 s_q sum_i d_{ijm} E2_i (conjugate J2_{1,1})
                let mut j1_21 = [[c(0.0, 0.0); 3]; 3];
                let mut j1_12 = [[c(0.0, 0.0); 3]; 3];
                let mut j2_11 = [[c(0.0, 0.0); 3]; 3];

                for i in 0..3 {
                    for j in 0..3 {
                        for k in 0..3 {
                            let dijk = sq * d_r3[i][j][k];
                            if dijk == 0.0 {
                                continue;
                            }
                            p_w2[i] += e1[j] * e1[k] * dijk;
                            j1_21[i][j] += e1[k] * (2.0 * dijk);
                            if include_backaction {
                                p_w1[j] += e2[i] * e1[k].conj() * (2.0 * dijk);
                                j1_12[j][i] += e1[k].conj() * (2.0 * dijk);
                                j2_11[j][k] += e2[i] * (2.0 * dijk);
                            }
                        }
                    }
                }

                let wk2_1 = -k0_2[0] * qp.wq_vol;
                let wk2_2 = -k0_2[1] * qp.wq_vol;
                let base_r1 = offsets[0];
                let base_r2 = offsets[1];
                for a in 0..nb {
                    let wa = qp.w_basis[a];
                    if include_backaction {
                        let (ra1, sc_a1) = c1.red_map[a];
                        if ra1 != usize::MAX {
                            let dot1 = p_w1[0] * wa[0] + p_w1[1] * wa[1] + p_w1[2] * wa[2];
                            let c_val = sc_a1.conj() * dot1 * wk2_1;
                            res_out[base_r1 + 2 * ra1] += c_val.re;
                            res_out[base_r1 + 2 * ra1 + 1] += c_val.im;
                        }
                    }
                    let (ra2, sc_a2) = c2.red_map[a];
                    if ra2 != usize::MAX {
                        let dot2 = p_w2[0] * wa[0] + p_w2[1] * wa[1] + p_w2[2] * wa[2];
                        let c_val = sc_a2.conj() * dot2 * wk2_2;
                        res_out[base_r2 + 2 * ra2] += c_val.re;
                        res_out[base_r2 + 2 * ra2 + 1] += c_val.im;
                    }
                }

                if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
                    let bc1 = c_offsets[0];
                    let bc2 = c_offsets[1];
                    for a in 0..nb {
                        let wa = qp.w_basis[a];
                        // Row in omega_2, col in omega_1 (J1_21):
                        let (ra2, sc_a2) = c2.red_map[a];
                        if ra2 != usize::MAX {
                            let sca_w2 = sc_a2.conj() * wk2_2;
                            for b in 0..nb {
                                let (rb1, sc_b1) = c1.red_map[b];
                                if rb1 == usize::MAX {
                                    continue;
                                }
                                let wb = qp.w_basis[b];
                                let mut m1 = c(0.0, 0.0);
                                for ii in 0..3 {
                                    for mm in 0..3 {
                                        m1 += j1_21[ii][mm] * (wa[ii] * wb[mm]);
                                    }
                                }
                                p.add_wirtinger_2x2(jv, bc2 + ra2, bc1 + rb1, sca_w2 * sc_b1 * m1, c(0.0, 0.0));
                            }
                        }
                        // Back-action rows in omega_1:
                        if include_backaction {
                            let (ra1, sc_a1) = c1.red_map[a];
                            if ra1 != usize::MAX {
                                let sca_w1 = sc_a1.conj() * wk2_1;
                                for b in 0..nb {
                                    let wb = qp.w_basis[b];
                                    // Col in omega_1 (J2_11):
                                    let (rb1, sc_b1) = c1.red_map[b];
                                    if rb1 != usize::MAX {
                                        let mut m2 = c(0.0, 0.0);
                                        for ii in 0..3 {
                                            for mm in 0..3 {
                                                m2 += j2_11[ii][mm] * (wa[ii] * wb[mm]);
                                            }
                                        }
                                        p.add_wirtinger_2x2(
                                            jv,
                                            bc1 + ra1,
                                            bc1 + rb1,
                                            c(0.0, 0.0),
                                            sca_w1 * sc_b1.conj() * m2,
                                        );
                                    }
                                    // Col in omega_2 (J1_12):
                                    let (rb2, sc_b2) = c2.red_map[b];
                                    if rb2 != usize::MAX {
                                        let mut m1 = c(0.0, 0.0);
                                        for ii in 0..3 {
                                            for mm in 0..3 {
                                                m1 += j1_12[ii][mm] * (wa[ii] * wb[mm]);
                                            }
                                        }
                                        p.add_wirtinger_2x2(jv, bc1 + ra1, bc2 + rb2, sca_w1 * sc_b2 * m1, c(0.0, 0.0));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else if m_cnt == 3 {
            // Non-degenerate 3-wave mixing (omega_3 = omega_1 + omega_2):
            // P_i(omega_3)/eps0 = 2 s_q sum_{j,k} d_{ijk} E_{1,j} E_{2,k}
            // P_j(omega_1)/eps0 = 2 s_q sum_{i,k} d_{ijk} E_{3,i} E_{2,k}^*
            // P_k(omega_2)/eps0 = 2 s_q sum_{i,j} d_{ijk} E_{3,i} E_{1,j}^*
            let c1 = &asm.active_cells[0][idx];
            let c2 = &asm.active_cells[1][idx];
            let c3 = &asm.active_cells[2][idx];
            let mut x1 = vec![c(0.0, 0.0); nb];
            let mut x2 = vec![c(0.0, 0.0); nb];
            let mut x3 = vec![c(0.0, 0.0); nb];
            for s in 0..nb {
                if c1.red_map[s].0 != usize::MAX {
                    x1[s] = c1.red_map[s].1 * y_c[0][c1.red_map[s].0];
                }
                if c2.red_map[s].0 != usize::MAX {
                    x2[s] = c2.red_map[s].1 * y_c[1][c2.red_map[s].0];
                }
                if c3.red_map[s].0 != usize::MAX {
                    x3[s] = c3.red_map[s].1 * y_c[2][c3.red_map[s].0];
                }
            }
            let undepleted_sfg = !include_backaction && asm.resolved_process == NonlinearProcess::Sfg;
            let undepleted_dfg =
                !include_backaction && matches!(asm.resolved_process, NonlinearProcess::Dfg | NonlinearProcess::Opa);

            for q_idx in 0..c1.qpts.len() {
                let qp = &c1.qpts[q_idx];
                let sq2 = 2.0 * d2.qpm.poling_sign_at(qp.pos);
                let mut e1 = [c(0.0, 0.0); 3];
                let mut e2 = [c(0.0, 0.0); 3];
                let mut e3 = [c(0.0, 0.0); 3];
                for s in 0..nb {
                    for k in 0..3 {
                        e1[k] += x1[s] * qp.w_basis[s][k];
                        e2[k] += x2[s] * qp.w_basis[s][k];
                        e3[k] += x3[s] * qp.w_basis[s][k];
                    }
                }
                let mut p_w = [[c(0.0, 0.0); 3]; 3];
                let mut j1_31 = [[c(0.0, 0.0); 3]; 3];
                let mut j1_32 = [[c(0.0, 0.0); 3]; 3];
                let mut j1_13 = [[c(0.0, 0.0); 3]; 3];
                let mut j2_12 = [[c(0.0, 0.0); 3]; 3];
                let mut j1_23 = [[c(0.0, 0.0); 3]; 3];
                let mut j2_21 = [[c(0.0, 0.0); 3]; 3];

                let active_w3 = !undepleted_dfg; // in undepleted DFG/OPA, pump omega_3 is fixed!
                let active_w12 = !undepleted_sfg; // in undepleted SFG, pumps omega_1 & omega_2 are fixed!

                for i in 0..3 {
                    for j in 0..3 {
                        for k in 0..3 {
                            let d_val = sq2 * d_r3[i][j][k];
                            if d_val == 0.0 {
                                continue;
                            }
                            if active_w3 {
                                p_w[2][i] += e1[j] * e2[k] * d_val;
                                j1_31[i][j] += e2[k] * d_val;
                                j1_32[i][k] += e1[j] * d_val;
                            }
                            if active_w12 {
                                p_w[0][j] += e3[i] * e2[k].conj() * d_val;
                                p_w[1][k] += e3[i] * e1[j].conj() * d_val;
                                j1_13[j][i] += e2[k].conj() * d_val;
                                j2_12[j][k] += e3[i] * d_val;
                                j1_23[k][i] += e1[j].conj() * d_val;
                                j2_21[k][j] += e3[i] * d_val;
                            }
                        }
                    }
                }

                let cells_3 = [c1, c2, c3];
                for m in 0..3 {
                    if (m == 2 && !active_w3) || (m < 2 && !active_w12) {
                        continue;
                    }
                    let wk2 = -k0_2[m] * qp.wq_vol;
                    let base_r = offsets[m];
                    for a in 0..nb {
                        let (ra, sc_a) = cells_3[m].red_map[a];
                        if ra == usize::MAX {
                            continue;
                        }
                        let wa = qp.w_basis[a];
                        let dot_m = p_w[m][0] * wa[0] + p_w[m][1] * wa[1] + p_w[m][2] * wa[2];
                        let cv = sc_a.conj() * dot_m * wk2;
                        res_out[base_r + 2 * ra] += cv.re;
                        res_out[base_r + 2 * ra + 1] += cv.im;
                    }
                }

                if let (Some(jv), Some(p)) = (jac_vals.as_deref_mut(), pat) {
                    let zero3 = [[c(0.0, 0.0); 3]; 3];
                    let couplings: &[(usize, usize, &[[C; 3]; 3], &[[C; 3]; 3], bool)] = &[
                        (2, 0, &j1_31, &zero3, active_w3),
                        (2, 1, &j1_32, &zero3, active_w3),
                        (0, 2, &j1_13, &zero3, active_w12),
                        (0, 1, &zero3, &j2_12, active_w12),
                        (1, 2, &j1_23, &zero3, active_w12),
                        (1, 0, &zero3, &j2_21, active_w12),
                    ];
                    for &(mr, mc, j1_b, j2_b, en) in couplings {
                        if !en {
                            continue;
                        }
                        let wk2 = -k0_2[mr] * qp.wq_vol;
                        let bcr = c_offsets[mr];
                        let bcc = c_offsets[mc];
                        for a in 0..nb {
                            let (ra, sc_a) = cells_3[mr].red_map[a];
                            if ra == usize::MAX {
                                continue;
                            }
                            let wa = qp.w_basis[a];
                            let sca_w = sc_a.conj() * wk2;
                            for b in 0..nb {
                                let (rb, sc_b) = cells_3[mc].red_map[b];
                                if rb == usize::MAX {
                                    continue;
                                }
                                let wb = qp.w_basis[b];
                                let (mut m1, mut m2) = (c(0.0, 0.0), c(0.0, 0.0));
                                for ii in 0..3 {
                                    for mm in 0..3 {
                                        let wab = wa[ii] * wb[mm];
                                        m1 += j1_b[ii][mm] * wab;
                                        m2 += j2_b[ii][mm] * wab;
                                    }
                                }
                                p.add_wirtinger_2x2(
                                    jv,
                                    bcr + ra,
                                    bcc + rb,
                                    sca_w * sc_b * m1,
                                    sca_w * sc_b.conj() * m2,
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

fn solve_multi_harmonic_coupled(
    mesh: &Discretisation,
    asm: &NonlinearDrivenAssembly,
    t0: Instant,
) -> Result<NonlinearDrivenSolve, String> {
    let spec = &asm.spec;
    let m_cnt = asm.harmonics.len();
    let omegas: Vec<f64> = asm.harmonic_freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let pat = asm.coupled_pattern.as_ref().ok_or("missing coupled_pattern")?;

    let mut offsets = vec![0usize; m_cnt + 1];
    for m in 0..m_cnt {
        offsets[m + 1] = offsets[m] + 2 * asm.harmonics[m].nred;
    }
    let n_tot_real = offsets[m_cnt];

    let (unit_pump_weights, default_power_w) = build_base_port_powers_w(spec, asm.harmonics[0].ports.len());
    let powers: Vec<f64> = if !spec.powers_w.is_empty() {
        spec.powers_w.clone()
    } else if spec.n_steps > 1 {
        let p_hi = if spec.power_max_w > 0.0 { spec.power_max_w } else { default_power_w };
        (0..spec.n_steps)
            .map(|k| spec.power_min_w + (p_hi - spec.power_min_w) * (k as f64) / ((spec.n_steps - 1) as f64))
            .collect()
    } else {
        vec![default_power_w]
    };

    let include_backaction = spec.coupling == Chi2CouplingMode::Depleted;
    let mut y_prev = vec![0.0_f64; n_tot_real];
    let mut steps = Vec::with_capacity(powers.len());
    let mut tot_iters = 0usize;
    let mut tot_facs = 0usize;
    let mut tot_sym = 0usize;

    // Optional JFNK block-diagonal linear preconditioner:
    let jfnk_pc = if spec.nl_solver == NonlinearSolverKind::Jfnk {
        let zero_rhs: Vec<Vec<C>> = asm.harmonics.iter().map(|h| vec![c(0.0, 0.0); h.nred]).collect();
        let mut dummy = vec![0.0; n_tot_real];
        let mut j0_vals = pat.zero_values();
        eval_chi2_residual_and_jacobian(asm, &y_prev, &zero_rhs, include_backaction, &mut dummy, Some(&mut j0_vals));
        let j0_csr = pat.into_csr(j0_vals);
        let mut fac = ShiftedFactorization::new(&j0_csr, None, Definiteness::General).map_err(|e| e.to_string())?;
        fac.set_shift(0.0).map_err(|e| e.to_string())?;
        tot_facs += 1;
        tot_sym += 1;
        Some(FactoredRealPc { fac })
    } else {
        None
    };

    for (step_idx, &pw) in powers.iter().enumerate() {
        femcore::progress::check()?;

        // Assign drive powers to each harmonic according to process:
        let mut rhs_hats = Vec::with_capacity(m_cnt);
        let mut a_ins = Vec::with_capacity(m_cnt);
        let mut p_ins = Vec::with_capacity(m_cnt);

        for m in 0..m_cnt {
            let h = &asm.harmonics[m];
            let np = h.ports.len();
            let mut sig_weights = vec![0.0; np];
            sig_weights[spec.signal_port.min(np - 1)] = 1.0;
            let mut sig_phases = vec![0.0; np];
            sig_phases[spec.signal_port.min(np - 1)] = spec.signal_phase_rad;

            let (rh, ai, pi) = match asm.resolved_process {
                NonlinearProcess::Shg => {
                    if m == 0 {
                        build_harmonic_drive_rhs(h, omegas[m], &unit_pump_weights, &spec.port_phases_rad, pw)
                    } else {
                        build_harmonic_drive_rhs(h, omegas[m], &sig_weights, &sig_phases, spec.signal_power_w)
                    }
                }
                NonlinearProcess::Chi3Xpm | NonlinearProcess::Sfg => {
                    if m == 0 {
                        build_harmonic_drive_rhs(h, omegas[m], &unit_pump_weights, &spec.port_phases_rad, pw)
                    } else if m == 1 {
                        build_harmonic_drive_rhs(h, omegas[m], &sig_weights, &sig_phases, spec.signal_power_w)
                    } else {
                        build_harmonic_drive_rhs(h, omegas[m], &vec![0.0; np], &[], 0.0)
                    }
                }
                NonlinearProcess::Dfg | NonlinearProcess::Opa => {
                    // Harmonic 2 is the high-frequency pump omega_3 = omega_1 + omega_2; harmonic 0 is signal omega_1!
                    if m == 2 {
                        build_harmonic_drive_rhs(h, omegas[m], &unit_pump_weights, &spec.port_phases_rad, pw)
                    } else if m == 0 {
                        build_harmonic_drive_rhs(h, omegas[m], &sig_weights, &sig_phases, spec.signal_power_w)
                    } else {
                        build_harmonic_drive_rhs(h, omegas[m], &vec![0.0; np], &[], 0.0)
                    }
                }
                _ => build_harmonic_drive_rhs(h, omegas[m], &unit_pump_weights, &spec.port_phases_rad, pw),
            };
            rhs_hats.push(rh);
            a_ins.push(ai);
            p_ins.push(pi);
        }

        // Fast exact sequential solve for undepleted SHG / SFG (each harmonic is an N x N complex linear solve!):
        if !include_backaction && matches!(asm.resolved_process, NonlinearProcess::Shg | NonlinearProcess::Sfg) {
            let mut y_seq = vec![0.0_f64; n_tot_real];
            for m in 0..(m_cnt - 1) {
                let a_m = asm.harmonics[m].assemble_matrix_at(omegas[m]);
                let mut fac_m =
                    ShiftedFactorization::new(&a_m, None, Definiteness::General).map_err(|e| e.to_string())?;
                fac_m.set_shift(c(0.0, 0.0)).map_err(|e| e.to_string())?;
                let mut ym = rhs_hats[m].clone();
                fac_m.solve_block_in_place(&mut ym, 1).map_err(|e| e.to_string())?;
                tot_facs += 1;
                tot_sym += 1;
                for (i, &zi) in ym.iter().enumerate() {
                    y_seq[offsets[m] + 2 * i] = zi.re;
                    y_seq[offsets[m] + 2 * i + 1] = zi.im;
                }
            }
            let last = m_cnt - 1;
            let mut res_nl_eval = vec![0.0_f64; n_tot_real];
            eval_chi2_residual_and_jacobian(asm, &y_seq, &rhs_hats, false, &mut res_nl_eval, None);
            let n_last = asm.harmonics[last].nred;
            let mut rhs_last: Vec<C> = (0..n_last)
                .map(|i| c(-res_nl_eval[offsets[last] + 2 * i], -res_nl_eval[offsets[last] + 2 * i + 1]))
                .collect();
            let a_last = asm.harmonics[last].assemble_matrix_at(omegas[last]);
            let mut fac_last =
                ShiftedFactorization::new(&a_last, None, Definiteness::General).map_err(|e| e.to_string())?;
            fac_last.set_shift(c(0.0, 0.0)).map_err(|e| e.to_string())?;
            fac_last.solve_block_in_place(&mut rhs_last, 1).map_err(|e| e.to_string())?;
            tot_facs += 1;
            tot_sym += 1;
            for (i, &zi) in rhs_last.iter().enumerate() {
                y_seq[offsets[last] + 2 * i] = zi.re;
                y_seq[offsets[last] + 2 * i + 1] = zi.im;
            }
            y_prev.copy_from_slice(&y_seq);

            let mut x_reds = Vec::with_capacity(m_cnt);
            let mut port_b_outs = Vec::with_capacity(m_cnt);
            let mut port_p_outs = Vec::with_capacity(m_cnt);
            for m in 0..m_cnt {
                let h = &asm.harmonics[m];
                let ym_c = join_real_to_complex(&y_prev[offsets[m]..offsets[m + 1]]);
                let xr: Vec<C> = ym_c.iter().zip(&h.d_ruiz).map(|(&y, &d)| y * d).collect();
                let (bo, po) = compute_port_outputs(h, omegas[m], &ym_c, &a_ins[m]);
                x_reds.push(xr);
                port_b_outs.push(bo);
                port_p_outs.push(po);
            }
            let u_j = compute_stored_energy_3d(mesh, &asm.harmonics[0], &x_reds[0]);
            steps.push(NonlinearDrivenStep {
                power_in_w: pw,
                x_red: x_reds,
                port_p_in_w: p_ins,
                port_p_out_w: port_p_outs,
                port_b_out: port_b_outs,
                stored_energy_j: u_j,
                max_carrier_density_m3: 0.0,
                tangent_alpha: 1.0,
                newton_iterations: 1,
                residual_norm: 0.0,
            });
            tot_iters += 1;
            continue;
        }

        // Fast undepleted seed on step 0 (solves undepleted sequential linear/coupled step in 1-2 iterations as warm start!):
        let mut y_init = y_prev.clone();
        if step_idx == 0 {
            let mut seed_res = vec![0.0; n_tot_real];
            let mut j0_vals = pat.zero_values();
            eval_chi2_residual_and_jacobian(asm, &y_init, &rhs_hats, false, &mut seed_res, Some(&mut j0_vals));
            let j0_csr = pat.into_csr(j0_vals);
            if let Ok(mut fac) = ShiftedFactorization::new(&j0_csr, None, Definiteness::General) {
                if fac.set_shift(0.0).is_ok() {
                    let mut d1: Vec<f64> = seed_res.iter().map(|&v| -v).collect();
                    if fac.solve_block_in_place(&mut d1, 1).is_ok() {
                        y_init = d1;
                        tot_facs += 1;
                        tot_sym += 1;
                    }
                }
            }
        }

        let rhs_tot_norm: f64 = rhs_hats.iter().flat_map(|v| v.iter()).map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        let atol = (1e-14 * rhs_tot_norm).max(1e-30);

        let nl_res = match spec.nl_solver {
            NonlinearSolverKind::Newton => {
                let mut f_clos = |x: &[f64], out: &mut [f64]| {
                    eval_chi2_residual_and_jacobian(asm, x, &rhs_hats, include_backaction, out, None);
                };
                let mut j_clos = |x: &[f64]| -> RCsr {
                    let mut dummy = vec![0.0; n_tot_real];
                    let mut jv = pat.zero_values();
                    eval_chi2_residual_and_jacobian(asm, x, &rhs_hats, include_backaction, &mut dummy, Some(&mut jv));
                    pat.into_csr(jv)
                };
                let opts = NewtonOptions {
                    rtol: spec.rtol,
                    atol,
                    stol: 1e-13,
                    max_it: spec.max_it,
                    line_search: LineSearch::default(),
                    jacobian_lag: 1,
                    definiteness: Definiteness::General,
                };
                newton(&mut f_clos, &mut j_clos, &y_init, &opts)
                    .map_err(|e| format!("coupled nonlinear Newton failed: {e}"))?
            }
            NonlinearSolverKind::Jfnk => {
                let mut f_clos = |x: &[f64], out: &mut [f64]| {
                    eval_chi2_residual_and_jacobian(asm, x, &rhs_hats, include_backaction, out, None);
                };
                let opts = JfnkOptions {
                    rtol: spec.rtol,
                    atol,
                    max_it: spec.max_it,
                    line_search: LineSearch::default(),
                    restart: 40,
                    max_linear_it: 160,
                    ..Default::default()
                };
                jfnk(&mut f_clos, &y_init, jfnk_pc.as_ref().map(|p| p as &dyn Preconditioner<f64>), &opts)
            }
        };

        if !nl_res.converged() {
            return Err(format!(
                "multi-harmonic '{}' nonlinear solve did not converge at P_in = {pw:.4e} W: {:?} (residual = {:.3e})",
                asm.resolved_process.as_str(),
                nl_res.reason,
                nl_res.residual_norms.last().copied().unwrap_or(f64::NAN)
            ));
        }

        tot_iters += nl_res.iterations;
        tot_facs += nl_res.factorizations;
        tot_sym += nl_res.symbolic_analyses;
        y_prev.copy_from_slice(&nl_res.x);

        let mut x_reds = Vec::with_capacity(m_cnt);
        let mut port_b_outs = Vec::with_capacity(m_cnt);
        let mut port_p_outs = Vec::with_capacity(m_cnt);
        for m in 0..m_cnt {
            let h = &asm.harmonics[m];
            let ym_c = join_real_to_complex(&y_prev[offsets[m]..offsets[m + 1]]);
            let xr: Vec<C> = ym_c.iter().zip(&h.d_ruiz).map(|(&y, &d)| y * d).collect();
            let (bo, po) = compute_port_outputs(h, omegas[m], &ym_c, &a_ins[m]);
            x_reds.push(xr);
            port_b_outs.push(bo);
            port_p_outs.push(po);
        }
        let u_j = compute_stored_energy_3d(mesh, &asm.harmonics[0], &x_reds[0]);

        steps.push(NonlinearDrivenStep {
            power_in_w: pw,
            x_red: x_reds,
            port_p_in_w: p_ins,
            port_p_out_w: port_p_outs,
            port_b_out: port_b_outs,
            stored_energy_j: u_j,
            max_carrier_density_m3: 0.0,
            tangent_alpha: 1.0,
            newton_iterations: nl_res.iterations,
            residual_norm: *nl_res.residual_norms.last().unwrap_or(&0.0),
        });
    }

    Ok(NonlinearDrivenSolve {
        steps,
        turning_point_powers_w: Vec::new(),
        total_newton_iterations: tot_iters,
        total_factorizations: tot_facs,
        symbolic_analyses: tot_sym,
        solve_seconds: t0.elapsed().as_secs_f64(),
    })
}

/// Stage 2 of `Study::NonlinearDriven`: solves the single-tone or multi-harmonic nonlinear FEM system.
pub fn solve_nonlinear_driven_system(
    mesh: &Discretisation,
    asm: &NonlinearDrivenAssembly,
) -> Result<NonlinearDrivenSolve, String> {
    let t0 = Instant::now();
    match asm.resolved_process {
        NonlinearProcess::Auto | NonlinearProcess::Chi3 => match asm.spec.continuation {
            NonlinearContinuation::PowerRamp => solve_chi3_power_ramp(mesh, asm, t0),
            NonlinearContinuation::Arclength => solve_chi3_arclength(mesh, asm, t0),
        },
        NonlinearProcess::Chi3Xpm
        | NonlinearProcess::Shg
        | NonlinearProcess::Sfg
        | NonlinearProcess::Dfg
        | NonlinearProcess::Opa => solve_multi_harmonic_coupled(mesh, asm, t0),
    }
}

fn build_nonlinear_3d_solution(
    mesh: &Discretisation,
    h: &DrivenAssembly,
    freq_hz: f64,
    x_red: &[C],
    label: String,
    power_in_w: f64,
) -> Solution {
    let x = h.red.expand(x_red);
    let geo = geometry(mesh);
    let lamn = local_bary(3, 2);
    let nloc = lamn.len();
    let nc = mesh.n_cells();
    let nb = h.ned.nb();

    let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    let mut curl_e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];

    for ci in 0..nc {
        let gid = h.ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for kk in 0..3 {
                gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
            }
        }
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        for (a, l) in lamn.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            match &curved {
                Some(ct) => {
                    let gq = safe_curved_jacobian(ct, &lam, &gs);
                    basis3_keyed(h.order, &lam, &gq, &gid, &mut w, &mut cu);
                }
                None => basis3_keyed(h.order, &lam, &gs, &gid, &mut w, &mut cu),
            }
            let mut val = [c(0.0, 0.0); 3];
            let mut cval = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let coef = x[h.ned.cell_dofs[ci][s]];
                for kk in 0..3 {
                    val[kk] += coef * w[s][kk];
                    cval[kk] += coef * cu[s][kk];
                }
            }
            e[ci * nloc + a] = val;
            curl_e[ci * nloc + a] = cval;
        }
    }

    let big =
        e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, 2);

    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
    ];

    let mut cells = CellFields::new(mesh, 2);
    for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
        let z: Vec<C> = e.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "V/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }
    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let h_prefactor = c(0.0, 1.0 / (omega * MU0));
    for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
        let z: Vec<C> = (0..nc * nloc)
            .map(|idx| {
                let ci = idx / nloc;
                let mui = c(1.0 / h.mur_phys[ci], 0.0);
                curl_e[idx][k] * mui * h_prefactor
            })
            .collect();
        cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }

    Solution {
        label,
        fields,
        cell_fields: cell_fields_of(&cells),
        value: Some(power_in_w),
        eigenvalue: Some([freq_hz, 0.0]),
    }
}

/// Stage 3 of `Study::NonlinearDriven`: extracts `SolveResult` solutions, port powers, conversion efficiencies,
/// Manley–Rowe photon-flux conservation metrics, and bistability turning points.
pub fn nonlinear_driven_fields(
    mesh: &Discretisation,
    asm: &NonlinearDrivenAssembly,
    solved: &NonlinearDrivenSolve,
) -> Result<SolveResult, String> {
    let m_cnt = asm.harmonics.len();
    let n_steps = solved.steps.len();
    let mut solutions = Vec::new();
    let mut scalars: Vec<(String, Scalar)> = vec![
        ("process".into(), Scalar::Text(asm.resolved_process.as_str().into())),
        ("continuation".into(), Scalar::Text(asm.spec.continuation.as_str().into())),
        ("n_steps".into(), Scalar::Num(n_steps as f64)),
        ("newton_iterations".into(), Scalar::Num(solved.total_newton_iterations as f64)),
        ("factorizations".into(), Scalar::Num(solved.total_factorizations as f64)),
        ("symbolic_analyses".into(), Scalar::Num(solved.symbolic_analyses as f64)),
    ];

    for (ti, &tp_w) in solved.turning_point_powers_w.iter().enumerate() {
        scalars.push((format!("turning_point_power {} (W)", ti + 1), Scalar::Num(tp_w)));
    }
    scalars.push(("n_turning_points".into(), Scalar::Num(solved.turning_point_powers_w.len() as f64)));

    // Record step-by-step continuation curves in JSON-friendly scalars + final/representative field solutions:
    let mut curve_p_in = Vec::with_capacity(n_steps);
    let mut curve_u_j = Vec::with_capacity(n_steps);
    let mut curve_ta = Vec::with_capacity(n_steps);

    for (si, st) in solved.steps.iter().enumerate() {
        curve_p_in.push(st.power_in_w);
        curve_u_j.push(st.stored_energy_j);
        curve_ta.push(st.tangent_alpha);

        let step_tag = if n_steps == 1 { String::new() } else { format!(" step {}", si + 1) };
        scalars.push((format!("P_in{step_tag} (W)"), Scalar::Num(st.power_in_w)));
        scalars.push((format!("U_stored{step_tag} (J)"), Scalar::Num(st.stored_energy_j)));
        if st.max_carrier_density_m3 > 0.0 {
            scalars.push((format!("max_N_c{step_tag} (1/m^3)"), Scalar::Num(st.max_carrier_density_m3)));
        }
        if asm.spec.continuation == NonlinearContinuation::Arclength {
            scalars.push((format!("tangent_alpha{step_tag}"), Scalar::Num(st.tangent_alpha)));
        }

        for m in 0..m_cnt {
            let h = &asm.harmonics[m];
            let p_tot_in: f64 = st.port_p_in_w[m].iter().sum();
            let p_tot_out: f64 = st.port_p_out_w[m].iter().sum();
            scalars.push((format!("harmonic {} P_in_total{step_tag} (W)", m + 1), Scalar::Num(p_tot_in)));
            scalars.push((format!("harmonic {} P_out_total{step_tag} (W)", m + 1), Scalar::Num(p_tot_out)));
            for (pi, port) in h.ports.iter().enumerate() {
                let bo = st.port_b_out[m][pi];
                scalars.push((
                    format!("harmonic {} port '{}' P_out{step_tag} (W)", m + 1, port.name),
                    Scalar::Num(st.port_p_out_w[m][pi]),
                ));
                scalars.push((
                    format!("harmonic {} port '{}' b_out_re{step_tag} (W^0.5)", m + 1, port.name),
                    Scalar::Num(bo.re),
                ));
                scalars.push((
                    format!("harmonic {} port '{}' b_out_im{step_tag} (W^0.5)", m + 1, port.name),
                    Scalar::Num(bo.im),
                ));
            }
        }

        // Conversion efficiency & Manley-Rowe photon-flux conservation for SHG / SFG / DFG / OPA:
        if asm.resolved_process == NonlinearProcess::Shg && m_cnt == 2 {
            let w1 = 2.0 * PI * asm.harmonic_freqs_hz[0];
            let w2 = 2.0 * PI * asm.harmonic_freqs_hz[1];
            let p1_in: f64 = st.port_p_in_w[0].iter().sum();
            let p2_in: f64 = st.port_p_in_w[1].iter().sum();
            let p1_out: f64 = st.port_p_out_w[0].iter().sum();
            let p2_out: f64 = st.port_p_out_w[1].iter().sum();
            let eta = if p1_in > 0.0 { p2_out / p1_in } else { 0.0 };
            let n_in = p1_in / (HBAR * w1) + 2.0 * p2_in / (HBAR * w2);
            let n_out = p1_out / (HBAR * w1) + 2.0 * p2_out / (HBAR * w2);
            let mr_rel = if n_in > 0.0 { (n_out - n_in).abs() / n_in } else { 0.0 };
            scalars.push((format!("shg_efficiency{step_tag}"), Scalar::Num(eta)));
            scalars.push((format!("manley_rowe_rel_err{step_tag}"), Scalar::Num(mr_rel)));
        } else if matches!(asm.resolved_process, NonlinearProcess::Sfg | NonlinearProcess::Dfg | NonlinearProcess::Opa)
            && m_cnt == 3
        {
            let n_flux = |m: usize, p_vec: &[f64]| -> f64 {
                let wm = 2.0 * PI * asm.harmonic_freqs_hz[m];
                p_vec.iter().sum::<f64>() / (HBAR * wm)
            };
            let n1_in = n_flux(0, &st.port_p_in_w[0]);
            let n2_in = n_flux(1, &st.port_p_in_w[1]);
            let n3_in = n_flux(2, &st.port_p_in_w[2]);
            let n1_out = n_flux(0, &st.port_p_out_w[0]);
            let n2_out = n_flux(1, &st.port_p_out_w[1]);
            let n3_out = n_flux(2, &st.port_p_out_w[2]);
            let scale = (n1_in + n2_in + 2.0 * n3_in).max(1e-30);
            // Manley-Rowe invariants: N_1 + N_3 = const and N_2 + N_3 = const!
            let mr_13 = ((n1_out + n3_out) - (n1_in + n3_in)).abs() / scale;
            let mr_23 = ((n2_out + n3_out) - (n2_in + n3_in)).abs() / scale;
            scalars.push((format!("manley_rowe_13_rel_err{step_tag}"), Scalar::Num(mr_13)));
            scalars.push((format!("manley_rowe_23_rel_err{step_tag}"), Scalar::Num(mr_23)));
            scalars.push((format!("manley_rowe_rel_err{step_tag}"), Scalar::Num(mr_13.max(mr_23))));
        }
    }

    // Emit full 3D fields for all steps when n_steps <= 8, or for the final step (plus turning-point neighborhoods) on long sweeps:
    let emit_all = n_steps <= 8;
    for (si, st) in solved.steps.iter().enumerate() {
        if !emit_all && si + 1 != n_steps && si != 0 && si != n_steps / 2 {
            continue;
        }
        for m in 0..m_cnt {
            let f_hz = asm.harmonic_freqs_hz[m];
            let lbl = if m_cnt == 1 {
                format!("step {}: P_in = {:.4e} W (f = {:.4e} Hz)", si + 1, st.power_in_w, f_hz)
            } else {
                format!("step {} harmonic {} (f = {:.4e} Hz, P_in = {:.4e} W)", si + 1, m + 1, f_hz, st.power_in_w)
            };
            solutions.push(build_nonlinear_3d_solution(
                mesh,
                &asm.harmonics[m],
                f_hz,
                &st.x_red[m],
                lbl,
                st.power_in_w,
            ));
        }
    }
    let _ = (curve_p_in, curve_u_j, curve_ta);

    let mut diag = crate::diag::EmsDiag::new(asm.spec.nl_solver.as_str());
    let final_res = solved.steps.last().map_or(0.0, |s| s.residual_norm);
    diag.record_mode_residuals("nonlinear_harmonic_balance", vec![final_res], true, None);
    let (diagnostics, error_budget) = diag.finish();

    let d_port_dummy: Option<DrivenPortData> = None;
    let _ = d_port_dummy;

    Ok(SolveResult {
        physics: "ems".into(),
        study: "nonlinear_driven".into(),
        order: asm.harmonics[0].order as u8,
        dim: 3,
        solutions,
        scalars,
        log: vec![format!(
            "nonlinear_driven ({}, {}): {} step(s), {} Newton iterations, {} factorizations ({} symbolic)",
            asm.resolved_process.as_str(),
            asm.spec.continuation.as_str(),
            n_steps,
            solved.total_newton_iterations,
            solved.total_factorizations,
            solved.symbolic_analyses,
        )],
        seconds: solved.solve_seconds,
        dofs: asm.harmonics.iter().map(|h| h.nred).sum(),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

/// One-call 3D `Study::NonlinearDriven` solver.
pub fn solve_nonlinear_driven(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let asm = assemble_nonlinear_driven(mesh, problem)?;
    let solved = solve_nonlinear_driven_system(mesh, &asm)?;
    nonlinear_driven_fields(mesh, &asm, &solved)
}

// =============================================================================================
// Stage 1 / Stage 2 / Stage 3 of 2D Study::NonlinearModes
// =============================================================================================

#[derive(Clone, Debug)]
struct QuadPoint2d {
    wq_area: f64,
    w_basis: Vec<[f64; 2]>,
    n_shape: Vec<f64>,
}

#[derive(Clone, Debug)]
struct ActiveCell2d {
    dom: i32,
    /// Index in `asm.free` (`usize::MAX` if PEC-clamped) for each local edge basis slot `0..nloc_e`.
    free_e: Vec<usize>,
    /// Index in `asm.free` (`usize::MAX` if PEC-clamped) for each local scalar `phi` slot `0..nl`.
    free_p: Vec<usize>,
    qpts: Vec<QuadPoint2d>,
}

/// Stage 1 of 2D `Study::NonlinearModes`: linear `ModeAssemblyC`, cell material tables, active Kerr/TPA cells,
/// and static $(2 N_{\text{free}} + 2) \times (2 N_{\text{free}} + 2)$ bordered Jacobian sparsity pattern.
pub struct NonlinearModesAssembly {
    pub spec: NonlinearModesSpec,
    pub order: usize,
    pub mode_asm: ModeAssemblyC,
    pub eps_cells: Vec<EpsC>,
    pub mu_cells: Vec<MuC>,
    pub pml_cells: Vec<bool>,
    active_cells: Vec<ActiveCell2d>,
    chi3_domains: BTreeMap<i32, DomainChi3Data>,
    pub bordered_pattern: StaticRealPattern,
}

/// Converged nonlinear waveguide mode at one guided optical power $P$ (W).
#[derive(Clone, Debug)]
pub struct NonlinearModeStep {
    pub power_w: f64,
    pub neff: C,
    /// Full unscaled transverse Nedelec DOFs `et` (`net`).
    pub et: Vec<C>,
    /// Full unscaled longitudinal potential DOFs `phi = tilde_phi / k0` (`nphi`), where $E_y = i\beta\phi$.
    pub phi: Vec<C>,
    pub newton_iterations: usize,
    pub residual_norm: f64,
}

/// Stage 2 output of `Study::NonlinearModes`.
#[derive(Clone, Debug)]
pub struct NonlinearModesSolve {
    pub neff_linear: C,
    /// Full-vectorial nonlinear waveguide parameter $\gamma = \frac{\omega}{4 P_0^2}\operatorname{Re}\int \mathbf{E}_0^* \cdot \mathbf{P}^{(3,\text{Kerr})}(\mathbf{E}_0)\,\mathrm{d}A$ ($\text{W}^{-1}\text{m}^{-1}$).
    pub gamma_w_inv_m_inv: f64,
    /// Imaginary nonlinear TPA waveguide parameter $\gamma_{\text{TPA}} = -\frac{\omega}{4 P_0^2}\operatorname{Im}\int \mathbf{E}_0^* \cdot \mathbf{P}_{\text{NL}}(\mathbf{E}_0)\,\mathrm{d}A$ ($\text{W}^{-1}\text{m}^{-1}$, where $\mathrm{d}P/\mathrm{d}y = -2\gamma_{\text{TPA}} P^2$).
    pub gamma_tpa_w_inv_m_inv: f64,
    /// Effective nonlinear mode area $A_{\text{eff}} = \frac{\big(\int \|\mathbf{E}_0\|^2\,\mathrm{d}A\big)^2}{\int_{\text{NL}} \|\mathbf{E}_0\|^4\,\mathrm{d}A}$ ($\text{m}^2$).
    pub a_eff_m2: f64,
    pub steps: Vec<NonlinearModeStep>,
    pub total_newton_iterations: usize,
    pub total_factorizations: usize,
    pub symbolic_analyses: usize,
    pub solve_seconds: f64,
}

/// Assembles Stage 1 (`NonlinearModesAssembly`) for a 2D `Study::NonlinearModes` problem.
pub fn assemble_nonlinear_modes(mesh: &Discretisation, problem: &Problem) -> Result<NonlinearModesAssembly, String> {
    if mesh.dim != 2 {
        return Err("the 'nonlinear_modes' study requires a 2D cross-section mesh".into());
    }
    let Study::NonlinearModes(ref spec) = problem.study else {
        return Err("assemble_nonlinear_modes requires Study::NonlinearModes".into());
    };
    if !(spec.freq > 0.0 && spec.freq.is_finite()) {
        return Err(format!("nonlinear_modes requires positive finite freq (got {})", spec.freq));
    }
    if spec.mode_index < 1 {
        return Err("nonlinear_modes requires mode_index >= 1".into());
    }
    if spec.power_w < 0.0 || spec.powers_w.iter().any(|&p| p < 0.0) {
        return Err("nonlinear_modes optical power must be non-negative".into());
    }

    let (eps_cells, mu_cells) = cell_materials(mesh, &problem.materials)?;
    let mode_asm = assemble_modes_complex_with_bcs_and_symmetry_and_sheets(
        mesh,
        spec.freq,
        &eps_cells,
        &mu_cells,
        problem.order,
        &problem.pec,
        &problem.surface_impedance,
        &problem.sheet_2d,
        problem.imported_potential.as_deref(),
        &problem.pmc,
        problem.symmetry,
    )?;

    let pml_cells: Vec<bool> = (0..mesh.n_cells())
        .map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(femcore::props::Pml::is_active))
        .collect();

    // Extract Chi3 domains on the 2D cross-section:
    let mut chi3_domains = BTreeMap::new();
    for (&dom, an) in &problem.materials.aniso {
        let Some(ref nl) = an.nonlinear else { continue };
        let Some(ref c3) = nl.chi3 else { continue };
        c3.validate()?;
        if !c3.is_active() {
            continue;
        }
        let n0 = an
            .permittivity
            .as_ref()
            .map(|t| {
                let m = t.matrix();
                ((m[0][0] + m[1][1] + m[2][2]) / 3.0).max(1.0).sqrt()
            })
            .or_else(|| problem.materials.get(dom, "eps_r").map(|e| e.max(1.0).sqrt()))
            .unwrap_or(1.0);
        let (a_mt, b_mt) = c3.isotropic_ab(n0);
        let has_cubic = c3.has_cubic_anisotropy(n0);
        let rank4_glob = c3.rank4_global(n0, an.euler_deg);
        chi3_domains.insert(
            dom,
            DomainChi3Data { chi3: c3.clone(), n0, euler_deg: an.euler_deg, has_cubic, a_mt, b_mt, rank4_glob },
        );
    }

    let n_free = mode_asm.free.len();
    let mut free_of_full = vec![usize::MAX; mode_asm.n_tot];
    for (fi, &gd) in mode_asm.free.iter().enumerate() {
        free_of_full[gd] = fi;
    }

    let geo = geometry(mesh);
    let cd = cell_data(mesh, problem.order);
    let nl = mesh.nloc(problem.order);
    let nloc_e = if problem.order == 1 { 3 } else { 8 };
    let quad = quad_rule(2, 6);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, problem.order)).collect();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut active_cells = Vec::new();

    for ci in 0..mesh.n_cells() {
        if pml_cells[ci] {
            continue;
        }
        let dom = mesh.cell_labels[ci];
        if !chi3_domains.contains_key(&dom) {
            continue;
        }
        mesh.cell_dofs(ci, problem.order, &mut pd);
        let free_e: Vec<usize> = cd.gdofs[ci].iter().map(|&gd| free_of_full[gd]).collect();
        let free_p: Vec<usize> = pd.iter().take(nl).map(|&p| free_of_full[mode_asm.net + p]).collect();
        if free_e.iter().all(|&f| f == usize::MAX) && free_p.iter().all(|&f| f == usize::MAX) {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        let mut qpts = Vec::with_capacity(quad.w.len());
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], problem.order, &mut w, &mut curl);
            qpts.push(QuadPoint2d { wq_area: quad.w[q] * area, w_basis: w.clone(), n_shape: tab[q].0.clone() });
        }
        active_cells.push(ActiveCell2d { dom, free_e, free_p, qpts });
    }

    // Build static (2 * n_free + 2) x (2 * n_free + 2) bordered CSR pattern for (Re u, Im u, Re beta, Im beta):
    let mut c_pairs = Vec::new();
    for r in 0..n_free {
        c_pairs.push((r, r));
        for &col in mode_asm.a.row(r).0 {
            c_pairs.push((r, col as usize));
        }
        for &col in mode_asm.b.row(r).0 {
            c_pairs.push((r, col as usize));
        }
        if let Some(ref a1) = mode_asm.a1 {
            for &col in a1.row(r).0 {
                c_pairs.push((r, col as usize));
            }
        }
    }
    for cell in &active_cells {
        let all_dofs: Vec<usize> =
            cell.free_e.iter().chain(cell.free_p.iter()).copied().filter(|&f| f != usize::MAX).collect();
        for &r in &all_dofs {
            for &s in &all_dofs {
                c_pairs.push((r, s));
            }
        }
    }
    c_pairs.sort_unstable();
    c_pairs.dedup();

    let n_real = 2 * n_free;
    let mut r_pairs = Vec::with_capacity(4 * c_pairs.len() + 8 * n_real + 8);
    for &(r, s) in &c_pairs {
        for dr in 0..2 {
            for ds in 0..2 {
                r_pairs.push((2 * r + dr, 2 * s + ds));
            }
        }
    }
    for i in 0..n_real {
        for b in [n_real, n_real + 1] {
            r_pairs.push((i, b));
            r_pairs.push((b, i));
        }
    }
    for b1 in [n_real, n_real + 1] {
        for b2 in [n_real, n_real + 1] {
            r_pairs.push((b1, b2));
        }
    }
    let bordered_pattern = StaticRealPattern::from_pairs(n_real + 2, n_real + 2, r_pairs);

    Ok(NonlinearModesAssembly {
        spec: spec.clone(),
        order: problem.order,
        mode_asm,
        eps_cells,
        mu_cells,
        pml_cells,
        active_cells,
        chi3_domains,
        bordered_pattern,
    })
}

/// Computes the exact cross-section time-averaged Poynting power $P(\mathbf{u}, \beta)$ (W) on physical (non-PML) cells
/// and its gradient $(\partial P/\partial(\operatorname{Re}\mathbf{u}, \operatorname{Im}\mathbf{u}), \partial P/\partial\operatorname{Re}\beta, \partial P/\partial\operatorname{Im}\beta)$.
fn eval_2d_poynting_power_and_grad(
    mesh: &Discretisation,
    asm: &NonlinearModesAssembly,
    u_free: &[C],
    beta: C,
) -> (f64, Vec<f64>, f64, f64) {
    let ma = &asm.mode_asm;
    let n_free = ma.free.len();
    let k0 = ma.k0;
    let omega = k0 * C0;
    let mut x_tot = vec![c(0.0, 0.0); ma.n_tot];
    for (fi, &gd) in ma.free.iter().enumerate() {
        x_tot[gd] = u_free[fi];
    }

    let geo = geometry(mesh);
    let cd = cell_data(mesh, asm.order);
    let nl = mesh.nloc(asm.order);
    let nloc_e = if asm.order == 1 { 3 } else { 8 };
    let quad = quad_rule(2, 4);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, asm.order)).collect();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];

    let mut free_of_full = vec![usize::MAX; ma.n_tot];
    for (fi, &gd) in ma.free.iter().enumerate() {
        free_of_full[gd] = fi;
    }

    let mut p_val = 0.0_f64;
    let mut dp_du_real = vec![0.0_f64; 2 * n_free];
    let mut dp_db_re = 0.0_f64;
    let mut dp_db_im = 0.0_f64;

    let inv_2om_mu0 = 0.5 / (omega * MU0);
    for ci in 0..mesh.n_cells() {
        if asm.pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        mesh.cell_dofs(ci, asm.order, &mut pd);
        let m = asm.mu_cells[ci];
        // Isotropic/diagonal mu: H_x = beta/(om mu0 mu_z_wafer) * (Ez + dphi/dz), H_z = -beta/(om mu0 mu_x_wafer) * (Ex + dphi/dx)
        // Note in local_tensor (x, z, y): m.x = mu_xx, m.y = mu_zz!
        let inv_mux = (c(1.0, 0.0) / m.x).re;
        let inv_muz = (c(1.0, 0.0) / m.y).re;

        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], asm.order, &mut w, &mut curl);
            let (_, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let wq = quad.w[q] * area * inv_2om_mu0;
            let (mut ex, mut ez) = (c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let cf = x_tot[cd.gdofs[ci][s]];
                ex += cf * w[s][0];
                ez += cf * w[s][1];
            }
            let (mut dphi_dx, mut dphi_dz) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nl {
                let pf = x_tot[ma.net + pd[a]] / k0;
                dphi_dx += pf * gn[a][0];
                dphi_dz += pf * gn[a][1];
            }
            // S_y = 0.5 Re(Ez^* Hx - Ex^* Hz) = inv_2om_mu0 * Re[ beta * (inv_muz * Ez^* (Ez + dphi_dz) + inv_mux * Ex^* (Ex + dphi_dx)) ]
            let q_term = ez.conj() * (ez + dphi_dz) * inv_muz + ex.conj() * (ex + dphi_dx) * inv_mux;
            p_val += wq * (beta * q_term).re;
            dp_db_re += wq * q_term.re;
            dp_db_im += -wq * q_term.im;

            // Wirtinger derivative dP / du_f: P = wq * Re(beta * q_term); for real G(u), dG = 2 Re( (dG/du) du ) => grad_re = 2 Re(dG/du), grad_im = -2 Im(dG/du)!
            // Let F_c = wq * beta * q_term; then P = Re(F_c), so:
            // dP/d(Re u_k) = Re( dF_c/du_k + (dF_c/du_k^*)^* ), dP/d(Im u_k) = Re( i dF_c/du_k - i (dF_c/du_k^*)^* ) = -Im( dF_c/du_k - (dF_c/du_k^*)^* )!
            let wb = beta * wq;
            for s in 0..nloc_e {
                let fi = free_of_full[cd.gdofs[ci][s]];
                if fi == usize::MAX {
                    continue;
                }
                let wx = w[s][0];
                let wz = w[s][1];
                // d(F_c)/du_s = wb * (inv_muz * ez^* * wz + inv_mux * ex^* * wx)
                // d(F_c)/du_s^* = wb * (inv_muz * (ez + dphi_dz) * wz + inv_mux * (ex + dphi_dx) * wx)
                let df_du = wb * (ez.conj() * (inv_muz * wz) + ex.conj() * (inv_mux * wx));
                let df_duc = wb * ((ez + dphi_dz) * (inv_muz * wz) + (ex + dphi_dx) * (inv_mux * wx));
                let g_c = df_du + df_duc.conj();
                dp_du_real[2 * fi] += g_c.re;
                dp_du_real[2 * fi + 1] += -g_c.im;
            }
            for a in 0..nl {
                let fi = free_of_full[ma.net + pd[a]];
                if fi == usize::MAX {
                    continue;
                }
                let gnx = gn[a][0] / k0;
                let gnz = gn[a][1] / k0;
                let df_du = wb * (ez.conj() * (inv_muz * gnz) + ex.conj() * (inv_mux * gnx));
                dp_du_real[2 * fi] += df_du.re;
                dp_du_real[2 * fi + 1] += -df_du.im;
            }
        }
    }
    (p_val, dp_du_real, dp_db_re, dp_db_im)
}

/// Evaluates the 2D nonlinear waveguide mode residual and bordered Jacobian on `(Re u, Im u, Re beta, Im beta)`!
fn eval_nonlinear_mode_bordered(
    mesh: &Discretisation,
    asm: &NonlinearModesAssembly,
    u_free: &[C],
    beta: C,
    u_ref: &[C],
    p_target: f64,
    p_scale: f64,
    res_out: &mut [f64],
    mut jac_vals: Option<&mut [f64]>,
) {
    let ma = &asm.mode_asm;
    let n_free = ma.free.len();
    let n_real = 2 * n_free;
    let k0 = ma.k0;
    let k0_2 = k0 * k0;
    let omega = k0 * C0;
    let jb = c(0.0, 1.0) * beta;
    let beta2 = beta * beta;

    // 1. Linear operator residual r_lin = (A + j beta A1 + beta^2 B) u_free:
    let au = ma.a.matvec(u_free);
    let bu = ma.b.matvec(u_free);
    let a1u = ma.a1.as_ref().map(|a1| a1.matvec(u_free));

    // Also compute dr/dbeta = (j A1 + 2 beta B) u_free:
    let mut dr_dbeta = vec![c(0.0, 0.0); n_free];
    for r in 0..n_free {
        let mut val = au[r] + beta2 * bu[r];
        let mut db = bu[r] * (2.0 * beta);
        if let Some(ref a1v) = a1u {
            val += jb * a1v[r];
            db += c(0.0, 1.0) * a1v[r];
        }
        res_out[2 * r] = val.re;
        res_out[2 * r + 1] = val.im;
        dr_dbeta[r] = db;
    }

    let pat = &asm.bordered_pattern;
    if let Some(jv) = jac_vals.as_deref_mut() {
        jv.fill(0.0);
        for r in 0..n_free {
            for (&col, &v) in ma.a.row(r).0.iter().zip(ma.a.row(r).1) {
                pat.add_wirtinger_2x2(jv, r, col as usize, v, c(0.0, 0.0));
            }
            for (&col, &v) in ma.b.row(r).0.iter().zip(ma.b.row(r).1) {
                pat.add_wirtinger_2x2(jv, r, col as usize, beta2 * v, c(0.0, 0.0));
            }
            if let Some(ref a1) = ma.a1 {
                for (&col, &v) in a1.row(r).0.iter().zip(a1.row(r).1) {
                    pat.add_wirtinger_2x2(jv, r, col as usize, jb * v, c(0.0, 0.0));
                }
            }
        }
    }

    // 2. Nonlinear Kerr + TPA + FCA/FCD cross-section integral:
    // E(x) = [ Ex, i*(beta/k0)*tilde_phi, Ez ] in wafer axes (0=x, 1=y, 2=z)!
    let ib_k0 = c(0.0, 1.0) * beta / k0;
    let i_k0 = c(0.0, 1.0) / k0;
    let mut dr_dbeta_conj = vec![c(0.0, 0.0); n_free]; // non-holomorphic d(r_NL)/d(beta^*) through Ey^* = -i*(beta^*/k0)*tilde_phi^*

    for cell in &asm.active_cells {
        let Some(d3) = asm.chi3_domains.get(&cell.dom) else { continue };
        let nloc_e = cell.free_e.len();
        let nl = cell.free_p.len();
        for qp in &cell.qpts {
            let (mut ex, mut ez, mut tphi) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let fi = cell.free_e[s];
                if fi != usize::MAX {
                    ex += u_free[fi] * qp.w_basis[s][0];
                    ez += u_free[fi] * qp.w_basis[s][1];
                }
            }
            for a in 0..nl {
                let fi = cell.free_p[a];
                if fi != usize::MAX {
                    tphi += u_free[fi] * qp.n_shape[a];
                }
            }
            let ey = ib_k0 * tphi;
            let e_wafer = [ex, ey, ez];
            let (p_nl, j1_pt, j2_pt, _) = eval_chi3_point_wirtinger(d3, omega, k0, &e_wafer);

            let wq = qp.wq_area;
            let wk2 = -k0_2 * wq;
            let ib_k0_wq = c(0.0, k0 * wq) * beta; // +i * beta * k0 * wq

            // d(Ey)/d(beta) = (i/k0)*tphi, d(Ey^*)/d(beta^*) = (-i/k0)*tphi^*:
            let dey_db = i_k0 * tphi;
            let deyc_dbc = (-i_k0) * tphi.conj();
            let dp_db = [j1_pt[0][1] * dey_db, j1_pt[1][1] * dey_db, j1_pt[2][1] * dey_db];
            let dp_dbc = [j2_pt[0][1] * deyc_dbc, j2_pt[1][1] * deyc_dbc, j2_pt[2][1] * deyc_dbc];

            // Transverse rows (s in 0..nloc_e):
            for s in 0..nloc_e {
                let fs = cell.free_e[s];
                if fs == usize::MAX {
                    continue;
                }
                let wx = qp.w_basis[s][0];
                let wz = qp.w_basis[s][1];
                let cv = (p_nl[0] * wx + p_nl[2] * wz) * wk2;
                res_out[2 * fs] += cv.re;
                res_out[2 * fs + 1] += cv.im;

                dr_dbeta[fs] += (dp_db[0] * wx + dp_db[2] * wz) * wk2;
                dr_dbeta_conj[fs] += (dp_dbc[0] * wx + dp_dbc[2] * wz) * wk2;

                if let Some(jv) = jac_vals.as_deref_mut() {
                    // Cols in transverse t:
                    for t in 0..nloc_e {
                        let ft = cell.free_e[t];
                        if ft == usize::MAX {
                            continue;
                        }
                        let wtx = qp.w_basis[t][0];
                        let wtz = qp.w_basis[t][1];
                        let m1 =
                            (j1_pt[0][0] * wtx + j1_pt[0][2] * wtz) * wx + (j1_pt[2][0] * wtx + j1_pt[2][2] * wtz) * wz;
                        let m2 =
                            (j2_pt[0][0] * wtx + j2_pt[0][2] * wtz) * wx + (j2_pt[2][0] * wtx + j2_pt[2][2] * wtz) * wz;
                        pat.add_wirtinger_2x2(jv, fs, ft, m1 * wk2, m2 * wk2);
                    }
                    // Cols in longitudinal tilde_phi_b (Ey = ib_k0 * tphi, Ey^* = conj(ib_k0) * tphi^*):
                    for b in 0..nl {
                        let fb = cell.free_p[b];
                        if fb == usize::MAX {
                            continue;
                        }
                        let nb_val = qp.n_shape[b];
                        let m1 = (j1_pt[0][1] * wx + j1_pt[2][1] * wz) * (ib_k0 * nb_val * wk2);
                        let m2 = (j2_pt[0][1] * wx + j2_pt[2][1] * wz) * (ib_k0.conj() * nb_val * wk2);
                        pat.add_wirtinger_2x2(jv, fs, fb, m1, m2);
                    }
                }
            }

            // Longitudinal rows (a in 0..nl): r_NL = +i * beta * k0 * int N_a * (P_NL,y / eps0) dA
            for a in 0..nl {
                let fa = cell.free_p[a];
                if fa == usize::MAX {
                    continue;
                }
                let na = qp.n_shape[a];
                let cv = ib_k0_wq * (p_nl[1] * na);
                res_out[2 * fa] += cv.re;
                res_out[2 * fa + 1] += cv.im;

                dr_dbeta[fa] += c(0.0, k0 * wq * na) * p_nl[1] + ib_k0_wq * (dp_db[1] * na);
                dr_dbeta_conj[fa] += ib_k0_wq * (dp_dbc[1] * na);

                if let Some(jv) = jac_vals.as_deref_mut() {
                    let row_pref = ib_k0_wq * na;
                    for t in 0..nloc_e {
                        let ft = cell.free_e[t];
                        if ft == usize::MAX {
                            continue;
                        }
                        let wtx = qp.w_basis[t][0];
                        let wtz = qp.w_basis[t][1];
                        let m1 = row_pref * (j1_pt[1][0] * wtx + j1_pt[1][2] * wtz);
                        let m2 = row_pref * (j2_pt[1][0] * wtx + j2_pt[1][2] * wtz);
                        pat.add_wirtinger_2x2(jv, fa, ft, m1, m2);
                    }
                    for b in 0..nl {
                        let fb = cell.free_p[b];
                        if fb == usize::MAX {
                            continue;
                        }
                        let nb_val = qp.n_shape[b];
                        let m1 = row_pref * j1_pt[1][1] * (ib_k0 * nb_val);
                        let m2 = row_pref * j2_pt[1][1] * (ib_k0.conj() * nb_val);
                        pat.add_wirtinger_2x2(jv, fa, fb, m1, m2);
                    }
                }
            }
        }
    }

    // Add beta columns (n_real, n_real + 1) to Jacobian:
    if let Some(jv) = jac_vals.as_deref_mut() {
        for r in 0..n_free {
            let j1 = dr_dbeta[r];
            let j2 = dr_dbeta_conj[r];
            let rr = 2 * r;
            pat.add_entry(jv, rr, n_real, (j1 + j2).re);
            pat.add_entry(jv, rr, n_real + 1, -(j1 - j2).im);
            pat.add_entry(jv, rr + 1, n_real, (j1 + j2).im);
            pat.add_entry(jv, rr + 1, n_real + 1, (j1 - j2).re);
        }
    }

    // 3. Constraint row n_real: Poynting power (P(u, beta) - p_target) * p_scale = 0
    let (p_cur, dp_du, dp_db_re, dp_db_im) = eval_2d_poynting_power_and_grad(mesh, asm, u_free, beta);
    res_out[n_real] = (p_cur - p_target) * p_scale;
    if let Some(jv) = jac_vals.as_deref_mut() {
        for col in 0..n_real {
            pat.add_entry(jv, n_real, col, dp_du[col] * p_scale);
        }
        pat.add_entry(jv, n_real, n_real, dp_db_re * p_scale);
        pat.add_entry(jv, n_real, n_real + 1, dp_db_im * p_scale);
    }

    // 4. Constraint row n_real + 1: Phase gauge Im(u_ref^H u_free) = sum_k (-Im(u_ref[k]) Re(u_k) + Re(u_ref[k]) Im(u_k)) = 0
    let mut gauge = 0.0_f64;
    for k in 0..n_free {
        gauge += -u_ref[k].im * u_free[k].re + u_ref[k].re * u_free[k].im;
    }
    res_out[n_real + 1] = gauge;
    if let Some(jv) = jac_vals {
        for k in 0..n_free {
            pat.add_entry(jv, n_real + 1, 2 * k, -u_ref[k].im);
            pat.add_entry(jv, n_real + 1, 2 * k + 1, u_ref[k].re);
        }
    }
}

/// Evaluates the full-vectorial waveguide nonlinear parameter $\gamma$ ($\text{W}^{-1}\text{m}^{-1}$),
/// TPA parameter $\gamma_{\text{TPA}}$, and effective mode area $A_{\text{eff}}$ ($\text{m}^2$) on the linear eigenmode `(u_0, beta_0)`.
fn compute_waveguide_gamma_and_aeff(
    mesh: &Discretisation,
    asm: &NonlinearModesAssembly,
    u0_free: &[C],
    beta0: C,
) -> (f64, f64, f64, f64) {
    let ma = &asm.mode_asm;
    let k0 = ma.k0;
    let omega = k0 * C0;
    let (p0, _, _, _) = eval_2d_poynting_power_and_grad(mesh, asm, u0_free, beta0);
    let p0_abs = p0.abs().max(1e-30);

    // Total int |E|^2 dA over non-PML cross-section for A_eff:
    let mut x_tot = vec![c(0.0, 0.0); ma.n_tot];
    for (fi, &gd) in ma.free.iter().enumerate() {
        x_tot[gd] = u0_free[fi];
    }
    let geo = geometry(mesh);
    let cd = cell_data(mesh, asm.order);
    let nl = mesh.nloc(asm.order);
    let nloc_e = if asm.order == 1 { 3 } else { 8 };
    let quad = quad_rule(2, 6);
    let tab: Vec<(Vec<f64>, Vec<f64>)> = quad.lam.iter().map(|l| shape(l, asm.order)).collect();
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let ib_k0 = c(0.0, 1.0) * beta0 / k0;

    let mut int_e2 = 0.0_f64;
    let mut int_e4_nl = 0.0_f64;
    for ci in 0..mesh.n_cells() {
        if asm.pml_cells[ci] {
            continue;
        }
        let g = tri_g(&geo, ci);
        let area = geo.vol[ci];
        mesh.cell_dofs(ci, asm.order, &mut pd);
        let is_nl = asm.chi3_domains.contains_key(&mesh.cell_labels[ci]);
        for (q, l) in quad.lam.iter().enumerate() {
            edge_basis(l, &g, &cd.sign[ci], asm.order, &mut w, &mut curl);
            let (mut ex, mut ez, mut tphi) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..nloc_e {
                let cf = x_tot[cd.gdofs[ci][s]];
                ex += cf * w[s][0];
                ez += cf * w[s][1];
            }
            for a in 0..nl {
                tphi += x_tot[ma.net + pd[a]] * tab[q].0[a];
            }
            let ey = ib_k0 * tphi;
            let e2 = ex.norm_sqr() + ey.norm_sqr() + ez.norm_sqr();
            let wq = quad.w[q] * area;
            int_e2 += e2 * wq;
            if is_nl {
                int_e4_nl += e2 * e2 * wq;
            }
        }
    }

    let mut int_ep_kerr = c(0.0, 0.0);
    for cell in &asm.active_cells {
        let Some(d3) = asm.chi3_domains.get(&cell.dom) else { continue };
        for qp in &cell.qpts {
            let (mut ex, mut ez, mut tphi) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for s in 0..cell.free_e.len() {
                let fi = cell.free_e[s];
                if fi != usize::MAX {
                    ex += u0_free[fi] * qp.w_basis[s][0];
                    ez += u0_free[fi] * qp.w_basis[s][1];
                }
            }
            for a in 0..cell.free_p.len() {
                let fi = cell.free_p[a];
                if fi != usize::MAX {
                    tphi += u0_free[fi] * qp.n_shape[a];
                }
            }
            let ey = ib_k0 * tphi;
            let e_w = [ex, ey, ez];
            let (p_over_eps0, _, _, _) = eval_chi3_point_wirtinger(d3, omega, k0, &e_w);
            let dot_ep =
                (e_w[0].conj() * p_over_eps0[0] + e_w[1].conj() * p_over_eps0[1] + e_w[2].conj() * p_over_eps0[2])
                    * EPS0;
            int_ep_kerr += dot_ep * qp.wq_area;
        }
    }

    let pref = omega / (4.0 * p0_abs * p0_abs);
    let gamma_re = pref * int_ep_kerr.re;
    let gamma_tpa = -pref * int_ep_kerr.im;
    let a_eff = if int_e4_nl > 0.0 { (int_e2 * int_e2) / int_e4_nl } else { 0.0 };
    (p0_abs, gamma_re, gamma_tpa, a_eff)
}

/// Stage 2 of `Study::NonlinearModes`: solves linear mode `mode_index`, evaluates exact vectorial $\gamma$ and $A_{\text{eff}}$,
/// and tracks power-dependent $n_{\text{eff}}(P)$ by fixed-power bordered Newton continuation.
pub fn solve_nonlinear_modes_continuation(
    mesh: &Discretisation,
    asm: &NonlinearModesAssembly,
) -> Result<NonlinearModesSolve, String> {
    let t0 = Instant::now();
    let spec = &asm.spec;
    let ma = &asm.mode_asm;
    let k0 = ma.k0;

    let target = if spec.evanescent { c(0.0, spec.n_target) } else { c(spec.n_target, 0.0) };
    let n_ask = (spec.mode_index + 3).max(4).min(ma.free.len().saturating_sub(1).max(1));
    let lin_modes = solve_mode_pairs_complex(ma, target, n_ask)?;
    if lin_modes.neff.len() < spec.mode_index {
        return Err(format!(
            "nonlinear_modes requested mode_index = {}, but only {} linear mode(s) converged near n_target = {}",
            spec.mode_index,
            lin_modes.neff.len(),
            spec.n_target
        ));
    }
    let m_idx = spec.mode_index - 1;
    let neff0 = lin_modes.neff[m_idx];
    let beta0 = neff0 * k0;

    // Reconstruct scaled free vector u0_free from (et, k0 * phi):
    let mut x0_tot = vec![c(0.0, 0.0); ma.n_tot];
    x0_tot[..ma.net].copy_from_slice(&lin_modes.et[m_idx]);
    for (p, &pv) in lin_modes.phi[m_idx].iter().enumerate() {
        x0_tot[ma.net + p] = pv * k0;
    }
    let mut u0_free: Vec<C> = ma.free.iter().map(|&gd| x0_tot[gd]).collect();

    let (p0_raw, gamma_re, gamma_tpa, a_eff) = compute_waveguide_gamma_and_aeff(mesh, asm, &u0_free, beta0);
    if p0_raw <= 1e-30 {
        return Err("selected waveguide mode carries zero propagating Poynting power (evanescent/cutoff)".into());
    }
    // Normalize u0_free to 1 Watt Poynting power and unit phase gauge!
    let scale_1w = (1.0 / p0_raw).sqrt();
    for v in &mut u0_free {
        *v *= scale_1w;
    }
    let u0_norm = u0_free.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt().max(1e-30);
    let u_ref: Vec<C> = u0_free.iter().map(|&z| z / u0_norm).collect();

    let powers: Vec<f64> = if !spec.powers_w.is_empty() { spec.powers_w.clone() } else { vec![spec.power_w] };

    let n_free = ma.free.len();
    let n_real = 2 * n_free;
    let n_sys = n_real + 2;
    let pat = &asm.bordered_pattern;

    let mut fac: Option<ShiftedFactorization<f64>> = None;
    let mut steps = Vec::with_capacity(powers.len());
    let mut tot_iters = 0usize;

    let mut prev_p = 0.0_f64;
    let mut u_prev = u0_free.clone();
    let mut beta_prev = beta0;

    for &pw in &powers {
        femcore::progress::check()?;
        if pw == 0.0 {
            let et_z = vec![c(0.0, 0.0); ma.net];
            let phi_z = vec![c(0.0, 0.0); ma.n_tot - ma.net];
            steps.push(NonlinearModeStep {
                power_w: 0.0,
                neff: neff0,
                et: et_z,
                phi: phi_z,
                newton_iterations: 0,
                residual_norm: 0.0,
            });
            continue;
        }

        // Initial guess at power `pw`:
        let (mut u_cur, mut beta_cur) = if prev_p > 0.0 {
            let ratio = (pw / prev_p).sqrt();
            (
                u_prev.iter().map(|&z| z * ratio).collect::<Vec<C>>(),
                beta_prev + c(gamma_re * (pw - prev_p), -gamma_tpa * (pw - prev_p)),
            )
        } else {
            let amp = pw.sqrt();
            (u0_free.iter().map(|&z| z * amp).collect::<Vec<C>>(), beta0 + c(gamma_re * pw, -gamma_tpa * pw))
        };

        // Characteristic residual scale ||A u_cur||_2 for balanced constraint scaling:
        let au_scale = ma.a.matvec(&u_cur).iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt().max(1.0);
        let p_scale = au_scale / pw.max(1e-30);

        let mut res = vec![0.0; n_sys];
        let mut jv = pat.zero_values();
        let mut converged = false;
        let mut it_used = 0usize;
        let mut rel_res = f64::INFINITY;

        for it in 0..spec.max_it {
            eval_nonlinear_mode_bordered(mesh, asm, &u_cur, beta_cur, &u_ref, pw, p_scale, &mut res, Some(&mut jv));
            let r_norm = res.iter().map(|v| v * v).sum::<f64>().sqrt();
            rel_res = r_norm / au_scale;
            if rel_res <= spec.rtol || r_norm <= 1e-13 * au_scale {
                converged = true;
                it_used = it;
                break;
            }

            let j_csr = pat.into_csr(jv.clone());
            let ff = match fac.as_mut() {
                Some(ff) => {
                    ff.replace_matrices(&j_csr, None).map_err(|e| e.to_string())?;
                    ff
                }
                None => fac
                    .insert(ShiftedFactorization::new(&j_csr, None, Definiteness::General).map_err(|e| e.to_string())?),
            };
            ff.set_shift(0.0).map_err(|e| e.to_string())?;

            let mut delta: Vec<f64> = res.iter().map(|&v| -v).collect();
            ff.solve_block_in_place(&mut delta, 1).map_err(|e| e.to_string())?;
            if delta.iter().any(|v| !v.is_finite()) {
                return Err(format!("nonlinear_modes linear solve produced non-finite step at P = {pw:.4e} W"));
            }

            // Armijo backtracking on ||res||_2:
            let mut alpha = 1.0_f64;
            let mut trial_res = vec![0.0; n_sys];
            let mut accepted_ls = false;
            for _bt in 0..12 {
                let u_try: Vec<C> =
                    (0..n_free).map(|k| u_cur[k] + c(alpha * delta[2 * k], alpha * delta[2 * k + 1])).collect();
                let b_try = beta_cur + c(alpha * delta[n_real], alpha * delta[n_real + 1]);
                eval_nonlinear_mode_bordered(mesh, asm, &u_try, b_try, &u_ref, pw, p_scale, &mut trial_res, None);
                let t_norm = trial_res.iter().map(|v| v * v).sum::<f64>().sqrt();
                if t_norm.is_finite() && t_norm <= (1.0 - 1e-4 * alpha) * r_norm {
                    u_cur = u_try;
                    beta_cur = b_try;
                    accepted_ls = true;
                    break;
                }
                alpha *= 0.5;
            }
            if !accepted_ls {
                for k in 0..n_free {
                    u_cur[k] += c(alpha * delta[2 * k], alpha * delta[2 * k + 1]);
                }
                beta_cur += c(alpha * delta[n_real], alpha * delta[n_real + 1]);
            }
            it_used = it + 1;
        }

        if !converged {
            return Err(format!(
                "nonlinear_modes Newton continuation did not converge at P = {pw:.4e} W (rel_res = {rel_res:.3e})"
            ));
        }

        tot_iters += it_used;
        prev_p = pw;
        u_prev = u_cur.clone();
        beta_prev = beta_cur;

        let mut x_full = vec![c(0.0, 0.0); ma.n_tot];
        for (fi, &gd) in ma.free.iter().enumerate() {
            x_full[gd] = u_cur[fi];
        }
        let et = x_full[..ma.net].to_vec();
        let phi: Vec<C> = x_full[ma.net..].iter().map(|&v| v / k0).collect();

        steps.push(NonlinearModeStep {
            power_w: pw,
            neff: beta_cur / k0,
            et,
            phi,
            newton_iterations: it_used,
            residual_norm: rel_res,
        });
    }

    let (fz, sy) = fac.as_ref().map_or((0, 0), |ff| (ff.stats().numeric_count, ff.stats().symbolic_count));
    Ok(NonlinearModesSolve {
        neff_linear: neff0,
        gamma_w_inv_m_inv: gamma_re,
        gamma_tpa_w_inv_m_inv: gamma_tpa,
        a_eff_m2: a_eff,
        steps,
        total_newton_iterations: tot_iters,
        total_factorizations: fz,
        symbolic_analyses: sy,
        solve_seconds: t0.elapsed().as_secs_f64(),
    })
}

/// Stage 3 of `Study::NonlinearModes`: builds `SolveResult` solutions and nonlinear waveguide scalars
/// (`neff_0`, `gamma (1/(W*m))`, `gamma_tpa (1/(W*m))`, `A_eff (m^2)`, `fom_tpa`, `neff(P)`, `delta_neff(P)`).
pub fn nonlinear_modes_fields(
    mesh: &Discretisation,
    asm: &NonlinearModesAssembly,
    solved: &NonlinearModesSolve,
) -> Result<SolveResult, String> {
    let freq = asm.spec.freq;
    let k0 = asm.mode_asm.k0;
    let mut solutions = Vec::with_capacity(solved.steps.len());
    let mut scalars: Vec<(String, Scalar)> = vec![
        ("neff_linear_re".into(), Scalar::Num(solved.neff_linear.re)),
        ("neff_linear_im".into(), Scalar::Num(solved.neff_linear.im)),
        ("gamma (1/(W*m))".into(), Scalar::Num(solved.gamma_w_inv_m_inv)),
        ("gamma_tpa (1/(W*m))".into(), Scalar::Num(solved.gamma_tpa_w_inv_m_inv)),
        ("A_eff (m^2)".into(), Scalar::Num(solved.a_eff_m2)),
    ];
    if solved.gamma_tpa_w_inv_m_inv > 1e-30 {
        // Standard two-photon figure of merit FOM_TPA = n_2 / (lambda_0 * beta_TPA) = gamma / (4 pi * gamma_TPA)
        let fom = solved.gamma_w_inv_m_inv.abs() / (4.0 * PI * solved.gamma_tpa_w_inv_m_inv);
        scalars.push(("fom_tpa".into(), Scalar::Num(fom)));
    }
    scalars.push(("newton_iterations".into(), Scalar::Num(solved.total_newton_iterations as f64)));

    let n_steps = solved.steps.len();
    for (si, st) in solved.steps.iter().enumerate() {
        let mc = ModesC {
            neff: vec![st.neff],
            et: vec![st.et.clone()],
            phi: vec![st.phi.clone()],
            et_left: vec![st.et.iter().map(|z| z.conj()).collect()],
            phi_left: vec![st.phi.iter().map(|z| z.conj()).collect()],
            biorth_norms: vec![c(1.0, 0.0)],
            is_non_symmetric: asm.mode_asm.is_non_symmetric,
            residuals: vec![st.residual_norm],
            restarts: 0,
            matvecs: 0,
            converged: true,
        };
        let (fields, cell_fields) = if st.power_w > 0.0 {
            mode_fields_complex_cells(mesh, &mc, 0, freq, asm.order)
        } else {
            let n = mesh.n_dofs(asm.order);
            let zf = vec![0.0; n];
            let fields = vec![
                Field::scalar("|E|", "a.u.", zf.clone()),
                Field::scalar("|Exz|", "a.u.", zf.clone()),
                Field::scalar("|Ey|", "a.u.", zf.clone()),
                Field::scalar("Re Ex", "a.u.", zf.clone()),
                Field::scalar("Re Ez", "a.u.", zf.clone()),
                Field::scalar("Re Ey", "a.u.", zf),
            ];
            let mut cells = CellFields::new(mesh, asm.order);
            let nc_tot = mesh.n_cells() * mesh.nloc(asm.order);
            for nm in ["Ex", "Ez", "Ey"] {
                cells.push_complex(nm, "a.u.", vec![0.0; nc_tot], vec![0.0; nc_tot]);
            }
            (fields, cells)
        };
        let txt = format!("{:.8} {:+.3e}i", st.neff.re, st.neff.im);
        solutions.push(Solution {
            label: format!("P = {:.4e} W: neff = {txt}", st.power_w),
            fields,
            cell_fields: cell_fields_of(&cell_fields),
            value: Some(st.power_w),
            eigenvalue: Some([st.neff.re, st.neff.im]),
        });

        let tag = if n_steps == 1 { String::new() } else { format!(" {}", si + 1) };
        let dneff = st.neff - solved.neff_linear;
        scalars.push((format!("power{tag} (W)"), Scalar::Num(st.power_w)));
        scalars.push((format!("neff_re{tag}"), Scalar::Num(st.neff.re)));
        scalars.push((format!("neff_im{tag}"), Scalar::Num(st.neff.im)));
        scalars.push((format!("delta_neff_re{tag}"), Scalar::Num(dneff.re)));
        scalars.push((format!("delta_neff_im{tag}"), Scalar::Num(dneff.im)));
        scalars.push((format!("delta_beta_re{tag} (1/m)"), Scalar::Num(dneff.re * k0)));
    }

    let mut diag = crate::diag::EmsDiag::new("nonlinear_mode_newton");
    let final_res = solved.steps.last().map_or(0.0, |s| s.residual_norm);
    diag.record_mode_residuals("nonlinear_modes", vec![final_res], true, None);
    let (diagnostics, error_budget) = diag.finish();

    Ok(SolveResult {
        physics: "ems".into(),
        study: "nonlinear_modes".into(),
        order: asm.order as u8,
        dim: 2,
        solutions,
        scalars,
        log: vec![format!(
            "nonlinear_modes: neff(0) = {:.6}, gamma = {:.4e} 1/(W*m), A_eff = {:.4e} m^2, {} power step(s)",
            solved.neff_linear.re, solved.gamma_w_inv_m_inv, solved.a_eff_m2, n_steps
        )],
        seconds: solved.solve_seconds,
        dofs: mesh.n_dofs(asm.order),
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

/// One-call 2D `Study::NonlinearModes` solver.
pub fn solve_nonlinear_modes(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let asm = assemble_nonlinear_modes(mesh, problem)?;
    let solved = solve_nonlinear_modes_continuation(mesh, &asm)?;
    nonlinear_modes_fields(mesh, &asm, &solved)
}
