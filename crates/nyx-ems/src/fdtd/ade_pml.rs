// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-pole Drude and Lorentz Auxiliary Differential Equations (ADEs), convolutional Complex-Frequency-Shifted PML
//! (`CPML`, Roden & Gedney 2000), adiabatic conductivity `Absorber`, PEC/PMC walls, and Bloch-periodic boundaries (ADR 145).

use super::grid::{FdtdDim, Pol2d, YeeGrid, EPS0, MU0};
use super::subpixel::YeeMaterialCoeffs;
use femcore::material::{Chi3, Dispersion, Nonlinear, QpmSpec};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Outer boundary condition along one face (`-x, +x, -y, +y, -z, +z`) of the rectilinear FDTD domain.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum FdtdBcFace {
    /// Convolutional Complex-Frequency-Shifted PML (`CPML`).
    Cpml {
        #[serde(default = "default_pml_layers")]
        layers: usize,
        #[serde(default = "default_pml_order")]
        order: f64,
        #[serde(default = "default_kappa_max")]
        kappa_max: f64,
        #[serde(default = "default_alpha_max")]
        alpha_max: f64,
        #[serde(default = "default_target_r0")]
        target_reflection: f64,
    },
    /// Polynomial adiabatic matched electric + magnetic conductivity absorber (`sigma_m / mu_0 = sigma_e / eps`).
    Absorber {
        #[serde(default = "default_absorber_layers")]
        layers: usize,
        #[serde(default = "default_pml_order")]
        order: f64,
        #[serde(default = "default_target_r0")]
        target_reflection: f64,
    },
    /// Perfect electric conductor (`E_tan = 0`).
    Pec,
    /// Perfect magnetic conductor (`H_tan = 0`, mirror even E).
    Pmc,
    /// Periodic or Bloch-periodic (`F(x + L) = exp(-j k_Bloch L) F(x)`).
    Periodic {
        #[serde(default)]
        bloch_k: f64,
    },
}

impl Default for FdtdBcFace {
    fn default() -> Self {
        Self::Cpml {
            layers: default_pml_layers(),
            order: default_pml_order(),
            kappa_max: default_kappa_max(),
            alpha_max: default_alpha_max(),
            target_reflection: default_target_r0(),
        }
    }
}

fn default_pml_layers() -> usize {
    10
}
fn default_absorber_layers() -> usize {
    16
}
fn default_pml_order() -> f64 {
    3.5
}
fn default_kappa_max() -> f64 {
    5.0
}
fn default_alpha_max() -> f64 {
    0.05
}
fn default_target_r0() -> f64 {
    1e-8
}

impl FdtdBcFace {
    pub fn default_cpml(layers: usize) -> Self {
        Self::Cpml {
            layers,
            order: default_pml_order(),
            kappa_max: default_kappa_max(),
            alpha_max: default_alpha_max(),
            target_reflection: default_target_r0(),
        }
    }

    pub fn is_periodic(&self) -> bool {
        matches!(self, Self::Periodic { .. })
    }

    pub fn bloch_k(&self) -> f64 {
        match self {
            Self::Periodic { bloch_k } => *bloch_k,
            _ => 0.0,
        }
    }
}

/// Six-face boundary specification `[x_minus, x_plus, y_minus, y_plus, z_minus, z_plus]`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FdtdBoundaries {
    pub faces: [FdtdBcFace; 6],
}

impl Default for FdtdBoundaries {
    fn default() -> Self {
        Self {
            faces: [
                FdtdBcFace::default_cpml(10),
                FdtdBcFace::default_cpml(10),
                FdtdBcFace::default_cpml(10),
                FdtdBcFace::default_cpml(10),
                FdtdBcFace::default_cpml(10),
                FdtdBcFace::default_cpml(10),
            ],
        }
    }
}

impl FdtdBoundaries {
    pub fn all_pec() -> Self {
        Self { faces: std::array::from_fn(|_| FdtdBcFace::Pec) }
    }

    pub fn periodic_axes(&self) -> [bool; 3] {
        [
            self.faces[0].is_periodic() && self.faces[1].is_periodic(),
            self.faces[2].is_periodic() && self.faces[3].is_periodic(),
            self.faces[4].is_periodic() && self.faces[5].is_periodic(),
        ]
    }

    pub fn bloch_vector(&self) -> [f64; 3] {
        [self.faces[0].bloch_k(), self.faces[2].bloch_k(), self.faces[4].bloch_k()]
    }

    pub fn is_complex_bloch(&self, dim: FdtdDim) -> bool {
        let bv = self.bloch_vector();
        bv[0].abs() > 0.0 || (!matches!(dim, FdtdDim::TwoD { .. }) && bv[1].abs() > 0.0) || bv[2].abs() > 0.0
    }

    pub fn check_finite(&self) -> Result<(), String> {
        for (i, f) in self.faces.iter().enumerate() {
            match f {
                FdtdBcFace::Cpml { order, kappa_max, alpha_max, target_reflection, .. } => {
                    if !order.is_finite()
                        || *order <= 0.0
                        || !kappa_max.is_finite()
                        || *kappa_max < 1.0
                        || !alpha_max.is_finite()
                        || *alpha_max < 0.0
                        || !target_reflection.is_finite()
                        || *target_reflection <= 0.0
                        || *target_reflection >= 1.0
                    {
                        return Err(format!("FDTD CPML boundary face {i} has invalid or non-finite parameters"));
                    }
                }
                FdtdBcFace::Absorber { order, target_reflection, .. } => {
                    if !order.is_finite()
                        || *order <= 0.0
                        || !target_reflection.is_finite()
                        || *target_reflection <= 0.0
                        || *target_reflection >= 1.0
                    {
                        return Err(format!("FDTD Absorber boundary face {i} has invalid or non-finite parameters"));
                    }
                }
                FdtdBcFace::Periodic { bloch_k } => {
                    if !bloch_k.is_finite() {
                        return Err(format!("FDTD Periodic boundary face {i} has non-finite bloch_k"));
                    }
                }
                FdtdBcFace::Pec | FdtdBcFace::Pmc => {}
            }
        }
        Ok(())
    }
}

/// 1D Roden–Gedney CPML stretched-derivative profile along one coordinate axis.
#[derive(Clone, Debug)]
pub struct CpmlAxis1d {
    /// Reciprocal real stretch `1 / kappa_e[i]` at primal vertices `0..=n` (used in Ampere's curl of H).
    pub inv_kappa_e: Vec<f64>,
    pub b_e: Vec<f64>,
    pub c_e: Vec<f64>,
    /// True where `i` lies inside the `-` or `+` CPML strip for E updates.
    pub active_e: Vec<bool>,
    /// Reciprocal real stretch `1 / kappa_h[i]` at dual centers `0..n` (used in Faraday's curl of E).
    pub inv_kappa_h: Vec<f64>,
    pub b_h: Vec<f64>,
    pub c_h: Vec<f64>,
    /// True where `i` lies inside the `-` or `+` CPML strip for H updates.
    pub active_h: Vec<bool>,
    /// Optional matched adiabatic conductivity profile `sigma_e[i]` (S/m) at primal vertices `0..=n`.
    pub absorber_sigma_e: Vec<f64>,
    /// Optional matched adiabatic magnetic loss `sigma_m[i] / mu_0` (S/F = 1/(s*eps_0)) at dual centers `0..n`.
    pub absorber_sigma_m_over_mu0: Vec<f64>,
}

impl CpmlAxis1d {
    pub fn build(v: &[f64], c: &[f64], bc_minus: &FdtdBcFace, bc_plus: &FdtdBcFace, dt: f64, enabled: bool) -> Self {
        let n = c.len();
        let mut out = Self {
            inv_kappa_e: vec![1.0; n + 1],
            b_e: vec![0.0; n + 1],
            c_e: vec![0.0; n + 1],
            active_e: vec![false; n + 1],
            inv_kappa_h: vec![1.0; n],
            b_h: vec![0.0; n],
            c_h: vec![0.0; n],
            active_h: vec![false; n],
            absorber_sigma_e: vec![0.0; n + 1],
            absorber_sigma_m_over_mu0: vec![0.0; n],
        };
        if !enabled || n < 3 {
            return out;
        }
        let eta0 = (MU0 / EPS0).sqrt();

        let apply_side = |out: &mut Self, face: &FdtdBcFace, is_minus: bool| match *face {
            FdtdBcFace::Cpml { layers, order, kappa_max, alpha_max, target_reflection } => {
                let nl = layers.min(n / 2);
                if nl == 0 {
                    return;
                }
                let x_iface = if is_minus { v[nl] } else { v[n - nl] };
                let x_wall = if is_minus { v[0] } else { v[n] };
                let d_pml = (x_wall - x_iface).abs().max(1e-30);
                let sigma_max = -(order + 1.0) * target_reflection.ln() / (2.0 * eta0 * d_pml);
                let alpha_si = alpha_max * EPS0 * super::grid::C0 / d_pml;

                for i in 0..=n {
                    let depth = if is_minus { x_iface - v[i] } else { v[i] - x_iface };
                    if depth > 0.0 {
                        let xi = (depth / d_pml).clamp(0.0, 1.0);
                        let sig = sigma_max * xi.powf(order);
                        let kap = 1.0 + (kappa_max - 1.0) * xi.powf(order);
                        let alp = alpha_si * (1.0 - xi);
                        let b = (-((sig / kap) + alp) * (dt / EPS0)).exp();
                        let denom = sig * kap + kap * kap * alp;
                        let c_val = if denom > 1e-30 { sig * (b - 1.0) / denom } else { 0.0 };
                        out.inv_kappa_e[i] = 1.0 / kap;
                        out.b_e[i] = b;
                        out.c_e[i] = c_val;
                        out.active_e[i] = true;
                    }
                }
                for i in 0..n {
                    let depth = if is_minus { x_iface - c[i] } else { c[i] - x_iface };
                    if depth > 0.0 {
                        let xi = (depth / d_pml).clamp(0.0, 1.0);
                        let sig = sigma_max * xi.powf(order);
                        let kap = 1.0 + (kappa_max - 1.0) * xi.powf(order);
                        let alp = alpha_si * (1.0 - xi);
                        let b = (-((sig / kap) + alp) * (dt / EPS0)).exp();
                        let denom = sig * kap + kap * kap * alp;
                        let c_val = if denom > 1e-30 { sig * (b - 1.0) / denom } else { 0.0 };
                        out.inv_kappa_h[i] = 1.0 / kap;
                        out.b_h[i] = b;
                        out.c_h[i] = c_val;
                        out.active_h[i] = true;
                    }
                }
            }
            FdtdBcFace::Absorber { layers, order, target_reflection } => {
                let nl = layers.min(n / 2);
                if nl == 0 {
                    return;
                }
                let x_iface = if is_minus { v[nl] } else { v[n - nl] };
                let x_wall = if is_minus { v[0] } else { v[n] };
                let d_abs = (x_wall - x_iface).abs().max(1e-30);
                let sigma_max = -(order + 1.0) * target_reflection.ln() / (2.0 * eta0 * d_abs);
                for i in 0..=n {
                    let depth = if is_minus { x_iface - v[i] } else { v[i] - x_iface };
                    if depth > 0.0 {
                        let xi = (depth / d_abs).clamp(0.0, 1.0);
                        out.absorber_sigma_e[i] += sigma_max * xi.powf(order);
                    }
                }
                for i in 0..n {
                    let depth = if is_minus { x_iface - c[i] } else { c[i] - x_iface };
                    if depth > 0.0 {
                        let xi = (depth / d_abs).clamp(0.0, 1.0);
                        out.absorber_sigma_m_over_mu0[i] += (sigma_max * xi.powf(order)) / EPS0;
                    }
                }
            }
            _ => {}
        };

        apply_side(&mut out, bc_minus, true);
        apply_side(&mut out, bc_plus, false);
        out
    }

    pub fn has_cpml(&self) -> bool {
        self.active_e.iter().any(|&a| a) || self.active_h.iter().any(|&a| a)
    }

    pub fn has_absorber(&self) -> bool {
        self.absorber_sigma_e.iter().any(|&s| s > 0.0)
    }
}

/// Auxiliary CPML convolutional memory arrays (`psi_{E_a, b}` and `psi_{H_a, b}`).
///
/// Allocated cleanly so `Step` updates remain cache-aligned and zero-overhead when an axis has no CPML.
#[derive(Clone, Debug)]
pub struct CpmlState {
    pub x: CpmlAxis1d,
    pub y: CpmlAxis1d,
    pub z: CpmlAxis1d,
    pub psi_ey_x: Vec<f64>,
    pub psi_ez_x: Vec<f64>,
    pub psi_ex_y: Vec<f64>,
    pub psi_ez_y: Vec<f64>,
    pub psi_ex_z: Vec<f64>,
    pub psi_ey_z: Vec<f64>,
    pub psi_hy_x: Vec<f64>,
    pub psi_hz_x: Vec<f64>,
    pub psi_hx_y: Vec<f64>,
    pub psi_hz_y: Vec<f64>,
    pub psi_hx_z: Vec<f64>,
    pub psi_hy_z: Vec<f64>,
}

impl CpmlState {
    pub fn new(grid: &YeeGrid, bcs: &FdtdBoundaries, dt: f64) -> Self {
        let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
        let x = CpmlAxis1d::build(&grid.vx, &grid.cx, &bcs.faces[0], &bcs.faces[1], dt, true);
        let y = CpmlAxis1d::build(&grid.vy, &grid.cy, &bcs.faces[2], &bcs.faces[3], dt, !is_2d);
        let z = CpmlAxis1d::build(&grid.vz, &grid.cz, &bcs.faces[4], &bcs.faces[5], dt, true);
        let n = grid.slab_len();
        let alloc = |on: bool| if on { vec![0.0; n] } else { Vec::new() };
        Self {
            psi_ey_x: alloc(x.has_cpml()),
            psi_ez_x: alloc(x.has_cpml()),
            psi_hy_x: alloc(x.has_cpml()),
            psi_hz_x: alloc(x.has_cpml()),
            psi_ex_y: alloc(y.has_cpml()),
            psi_ez_y: alloc(y.has_cpml()),
            psi_hx_y: alloc(y.has_cpml()),
            psi_hz_y: alloc(y.has_cpml()),
            psi_ex_z: alloc(z.has_cpml()),
            psi_ey_z: alloc(z.has_cpml()),
            psi_hx_z: alloc(z.has_cpml()),
            psi_hy_z: alloc(z.has_cpml()),
            x,
            y,
            z,
        }
    }
}

// ------------------------------------------------------------------ Multi-Pole Drude / Lorentz ADEs

/// One Drude or Lorentz susceptibility pole in SI angular frequencies (rad/s) under the `exp(+j omega t)` convention:
/// - Drude (`omega_0 == 0`): `chi(omega) = -strength / (omega^2 - j gamma omega)`
/// - Lorentz (`omega_0 > 0`): `chi(omega) = +strength / (omega_0^2 - omega^2 + j gamma omega)`
///
/// Here `strength = f_k * omega_p^2` (rad^2 / s^2).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SusceptibilityPole {
    /// Resonance angular frequency `omega_0` (rad/s; `0.0` for a Drude pole).
    pub omega_0: f64,
    /// Damping rate `gamma` (rad/s, `>= 0`).
    pub gamma: f64,
    /// Oscillator strength `f_k * omega_p^2` (rad^2 / s^2, `> 0`).
    pub strength: f64,
}

/// Dispersive ADE material model attached to a 1-based domain label.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DispersiveDomainSpec {
    pub eps_inf: f64,
    pub poles: Vec<SusceptibilityPole>,
}

impl DispersiveDomainSpec {
    /// Converts a `femcore::material::Dispersion::Oscillators` model into an `FDTD` ADE pole specification.
    pub fn from_optical_dispersion(disp: &Dispersion) -> Option<Self> {
        let Dispersion::Oscillators { eps_inf, plasma, drude, lorentz, .. } = disp else {
            return None;
        };
        let wp2 = plasma * plasma;
        let mut poles = Vec::new();
        if let Some([f0, g0]) = *drude {
            if f0 * wp2 > 0.0 {
                poles.push(SusceptibilityPole { omega_0: 0.0, gamma: g0.max(0.0), strength: f0 * wp2 });
            }
        }
        for &[f, w0, g] in lorentz {
            if f * wp2 > 0.0 && w0 > 0.0 {
                poles.push(SusceptibilityPole { omega_0: w0, gamma: g.max(0.0), strength: f * wp2 });
            }
        }
        (!poles.is_empty()).then_some(Self { eps_inf: eps_inf.max(1.0), poles })
    }

    /// Evaluates exact analytical complex relative permittivity `eps_r(omega)` at frequency `freq_hz` (`exp(+j omega t)`).
    pub fn eval_eps_r(&self, freq_hz: f64) -> [f64; 2] {
        let w = 2.0 * std::f64::consts::PI * freq_hz;
        let (mut re, mut im) = (self.eps_inf, 0.0);
        for p in &self.poles {
            if p.omega_0 == 0.0 {
                // -strength / (w^2 - j gamma w) = -strength * (w^2 + j gamma w) / (w^4 + gamma^2 w^2)
                let denom = w * w * (w * w + p.gamma * p.gamma).max(1e-30);
                re -= p.strength * w * w / denom;
                im -= p.strength * p.gamma * w / denom;
            } else {
                // +strength / ((w0^2 - w^2) + j gamma w)
                let dr = p.omega_0 * p.omega_0 - w * w;
                let di = p.gamma * w;
                let denom = (dr * dr + di * di).max(1e-30);
                re += p.strength * dr / denom;
                im -= p.strength * di / denom;
            }
        }
        [re, im]
    }
}

/// Precomputed ADE pole update coefficients for one `SusceptibilityPole` at time step `dt`.
#[derive(Clone, Copy, Debug)]
pub enum CompiledPole {
    /// Trapezoidal Drude current `J^{n+1} = k_d J^n + beta_d (E^{n+1} + E^n)`.
    Drude {
        k_d: f64,
        /// `beta_d = eps_0 * strength * dt / (2 * (1 + 0.5 * gamma * dt))`
        beta_d: f64,
    },
    /// Centered second-order Lorentz polarization `P^{n+1} = c1 P^n + c2 P^{n-1} + c3 E^n`.
    Lorentz { c1: f64, c2: f64, c3: f64 },
}

impl CompiledPole {
    pub fn compile(p: &SusceptibilityPole, dt: f64) -> Self {
        if p.omega_0 <= 0.0 {
            let g = 0.5 * p.gamma * dt;
            let k_d = (1.0 - g) / (1.0 + g);
            let beta_d = (EPS0 * p.strength * 0.5 * dt) / (1.0 + g);
            Self::Drude { k_d, beta_d }
        } else {
            let g = 0.5 * p.gamma * dt;
            let w2 = p.omega_0 * p.omega_0 * dt * dt;
            let denom = 1.0 + g;
            Self::Lorentz { c1: (2.0 - w2) / denom, c2: (g - 1.0) / denom, c3: (EPS0 * p.strength * dt * dt) / denom }
        }
    }
}

/// Per-component sparse ADE bank updating all dispersive Yee edges on one component (`Ex`, `Ey`, or `Ez`).
#[derive(Clone, Debug, Default)]
pub struct ComponentAdeBank {
    /// Flat SoA indices into the electric field component slab.
    pub edge_indices: Vec<usize>,
    /// Pole range `pole_offsets[k]..pole_offsets[k+1]` for dispersive edge `k`.
    pub pole_offsets: Vec<usize>,
    pub poles: Vec<CompiledPole>,
    /// Auxiliary state slot `s1[m]` (`J^n` for Drude, `P^n` for Lorentz).
    pub s1: Vec<f64>,
    /// Auxiliary state slot `s2[m]` (unused for Drude, `P^{n-1}` for Lorentz).
    pub s2: Vec<f64>,
    /// Additional effective conductivity `sum_{Drude} beta_d` folded into the implicit `E^{n+1}` update denominator!
    pub drude_cond_add: Vec<f64>,
}

impl ComponentAdeBank {
    pub fn build(dom_arr: &[i32], disp_frac: &[f64], specs: &BTreeMap<i32, DispersiveDomainSpec>, dt: f64) -> Self {
        if specs.is_empty() {
            return Self::default();
        }
        let mut edge_indices = Vec::new();
        let mut pole_offsets = vec![0usize];
        let mut poles = Vec::new();
        let mut drude_cond_add = Vec::new();

        for (idx, &dom) in dom_arr.iter().enumerate() {
            let Some(sp) = specs.get(&dom) else { continue };
            let frac = disp_frac.get(idx).copied().unwrap_or(1.0);
            if sp.poles.is_empty() || frac <= 0.0 {
                continue;
            }
            edge_indices.push(idx);
            let mut d_add = 0.0;
            for p in &sp.poles {
                let mut sp_p = *p;
                sp_p.strength *= frac;
                let cp = CompiledPole::compile(&sp_p, dt);
                if let CompiledPole::Drude { beta_d, .. } = cp {
                    d_add += beta_d;
                }
                poles.push(cp);
            }
            drude_cond_add.push(d_add);
            pole_offsets.push(poles.len());
        }
        let n_poles = poles.len();
        Self { edge_indices, pole_offsets, poles, s1: vec![0.0; n_poles], s2: vec![0.0; n_poles], drude_cond_add }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.edge_indices.is_empty()
    }
}

/// Pre-rotated per-domain nonlinear optical constitutive parameters (`FR-41` `[nonlinear]`, `ADR 143`).
#[derive(Clone, Debug)]
pub struct NonlinearFdtdDomain {
    pub chi2_d_global: Option<[[f64; 6]; 3]>,
    pub qpm: Option<QpmSpec>,
    pub chi3: Option<Chi3>,
    pub chi3_iso_xxxx: f64,
    pub chi3_rank4_global: Option<Box<[[[[f64; 3]; 3]; 3]; 3]>>,
    pub euler_deg: [f64; 3],
    pub n0: f64,
    pub omega_0: f64,
    pub eps_min: f64,
    pub beta_tpa: f64,
    pub sigma_fca: f64,
    pub k_fcd: f64,
    pub tau_c: f64,
}

/// Active nonlinear Yee electric-field edge entry.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearEdgeEntry {
    pub idx: usize,
    pub ix: u32,
    pub iy: u32,
    pub iz: u32,
    pub dom_slot: u16,
    pub qpm_sign: f64,
}

/// Sparse nonlinear polarization & free-carrier ADE bank across `(Ex, Ey, Ez)` Yee edges (`ADR 143`, `ADR 145`).
#[derive(Clone, Debug, Default)]
pub struct NonlinearFdtdBank {
    pub domains: Vec<NonlinearFdtdDomain>,
    pub edges: [Vec<NonlinearEdgeEntry>; 3],
    /// Saved nonlinear polarization $P_{\text{NL},a}^n$ ($\text{C/m}^2$) at step $n$ for each active edge.
    pub p_nl_prev: [Vec<f64>; 3],
    /// Free-carrier density $N_c^n$ ($\text{m}^{-3}$) at step $n$ when TPA carrier dynamics are active.
    pub n_carrier: [Vec<f64>; 3],
    pub has_carriers: bool,
}

impl NonlinearFdtdBank {
    #[inline]
    pub fn is_active(&self) -> bool {
        !self.domains.is_empty()
            && (!self.edges[0].is_empty() || !self.edges[1].is_empty() || !self.edges[2].is_empty())
    }

    pub fn build(
        grid: &YeeGrid,
        mats: &YeeMaterialCoeffs,
        cpml: &CpmlState,
        nl_specs: &BTreeMap<i32, (Nonlinear, [f64; 3], f64, f64)>,
        omega_0: f64,
    ) -> Result<Self, String> {
        if nl_specs.is_empty() {
            return Ok(Self::default());
        }
        let mut domains = Vec::new();
        let mut dom_to_slot = BTreeMap::new();
        let mut has_carriers = false;

        for (&dom, (nl, euler_deg, n0, eps_min)) in nl_specs {
            if let Some(ref q) = nl.qpm {
                q.validate()?;
            }
            if let Some(ref c3) = nl.chi3 {
                c3.validate()?;
            }
            let chi2_d_global = nl.chi2.as_ref().and_then(|c2| {
                let dg = c2.d_global(*euler_deg);
                let max_d = dg.iter().flatten().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
                (max_d > 0.0).then_some(dg)
            });
            let (chi3_opt, chi3_iso_xxxx, chi3_rank4_global, beta_tpa, sigma_fca, k_fcd, tau_c) =
                if let Some(ref c3) = nl.chi3 {
                    let has_cub = c3.has_cubic_anisotropy(*n0);
                    let (a_mt, b_mt) = c3.isotropic_ab(*n0);
                    let iso_xxxx = if has_cub { 0.0 } else { (4.0 / 3.0) * (a_mt + 0.5 * b_mt) };
                    let r4 = has_cub.then(|| Box::new(c3.rank4_global(*n0, *euler_deg)));
                    let bt = c3.beta_tpa.unwrap_or(0.0).max(0.0);
                    let sf = c3.sigma_fca.unwrap_or(0.0).max(0.0);
                    let kf = c3.k_fcd.unwrap_or(0.0);
                    let tc = c3.tau_c.unwrap_or(0.0).max(0.0);
                    if bt > 0.0 && (sf > 0.0 || kf != 0.0) {
                        has_carriers = true;
                    }
                    (Some(c3.clone()), iso_xxxx, r4, bt, sf, kf, tc)
                } else {
                    (None, 0.0, None, 0.0, 0.0, 0.0, 0.0)
                };
            if chi2_d_global.is_none() && chi3_rank4_global.is_none() && chi3_iso_xxxx.abs() == 0.0 && beta_tpa == 0.0 {
                continue;
            }
            let slot = domains.len() as u16;
            dom_to_slot.insert(dom, slot);
            domains.push(NonlinearFdtdDomain {
                chi2_d_global,
                qpm: nl.qpm.clone(),
                chi3: chi3_opt,
                chi3_iso_xxxx,
                chi3_rank4_global,
                euler_deg: *euler_deg,
                n0: n0.max(1.0),
                omega_0: omega_0.max(1.0),
                eps_min: eps_min.max(1.0),
                beta_tpa,
                sigma_fca,
                k_fcd,
                tau_c,
            });
        }
        if domains.is_empty() {
            return Ok(Self::default());
        }

        let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
        let (active_x, active_y, active_z) = match grid.dim {
            FdtdDim::TwoD { polarization: Pol2d::TmEy } => (false, true, false),
            FdtdDim::TwoD { polarization: Pol2d::TeHy } => (true, false, true),
            FdtdDim::TwoD { polarization: Pol2d::Full } | FdtdDim::ThreeD => (true, true, true),
        };

        let in_abs_v = |ax: &CpmlAxis1d, i: usize| ax.active_e[i] || ax.absorber_sigma_e[i] > 0.0;
        let in_abs_c = |ax: &CpmlAxis1d, i: usize| ax.active_h[i] || ax.absorber_sigma_m_over_mu0[i] > 0.0;

        let pec_lx = mats.pec_len_x.as_slice();
        let pec_ly = mats.pec_len_y.as_slice();
        let pec_lz = mats.pec_len_z.as_slice();
        let inv_ex = mats.inv_eps_x.as_slice();
        let inv_ey = mats.inv_eps_y.as_slice();
        let inv_ez = mats.inv_eps_z.as_slice();

        let mut edges_x = Vec::new();
        let mut edges_y = Vec::new();
        let mut edges_z = Vec::new();

        if active_x {
            let iy_hi = if is_2d { 0 } else { grid.ny };
            for ix in 0..grid.nx {
                if in_abs_c(&cpml.x, ix) {
                    continue;
                }
                let x = grid.cx[ix];
                for iy in 0..=iy_hi {
                    if !is_2d && in_abs_v(&cpml.y, iy) {
                        continue;
                    }
                    let y = if is_2d { 0.0 } else { grid.vy[iy] };
                    for iz in 0..=grid.nz {
                        if in_abs_v(&cpml.z, iz) {
                            continue;
                        }
                        let g = grid.idx(ix, iy, iz);
                        if pec_lx[g] <= 0.0 || inv_ex[g] <= 0.0 {
                            continue;
                        }
                        let Some(&slot) = dom_to_slot.get(&mats.dom_ex[g]) else { continue };
                        let z = grid.vz[iz];
                        let qpm_sign = domains[slot as usize].qpm.as_ref().map_or(1.0, |q| q.poling_sign_at([x, y, z]));
                        edges_x.push(NonlinearEdgeEntry {
                            idx: g,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            dom_slot: slot,
                            qpm_sign,
                        });
                    }
                }
            }
        }

        if active_y {
            let iy_len = if is_2d { 1 } else { grid.ny };
            for ix in 0..=grid.nx {
                if in_abs_v(&cpml.x, ix) {
                    continue;
                }
                let x = grid.vx[ix];
                for iy in 0..iy_len {
                    if !is_2d && in_abs_c(&cpml.y, iy) {
                        continue;
                    }
                    let y = if is_2d { 0.0 } else { grid.cy[iy] };
                    for iz in 0..=grid.nz {
                        if in_abs_v(&cpml.z, iz) {
                            continue;
                        }
                        let g = grid.idx(ix, iy, iz);
                        if pec_ly[g] <= 0.0 || inv_ey[g] <= 0.0 {
                            continue;
                        }
                        let Some(&slot) = dom_to_slot.get(&mats.dom_ey[g]) else { continue };
                        let z = grid.vz[iz];
                        let qpm_sign = domains[slot as usize].qpm.as_ref().map_or(1.0, |q| q.poling_sign_at([x, y, z]));
                        edges_y.push(NonlinearEdgeEntry {
                            idx: g,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            dom_slot: slot,
                            qpm_sign,
                        });
                    }
                }
            }
        }

        if active_z {
            let iy_hi = if is_2d { 0 } else { grid.ny };
            for ix in 0..=grid.nx {
                if in_abs_v(&cpml.x, ix) {
                    continue;
                }
                let x = grid.vx[ix];
                for iy in 0..=iy_hi {
                    if !is_2d && in_abs_v(&cpml.y, iy) {
                        continue;
                    }
                    let y = if is_2d { 0.0 } else { grid.vy[iy] };
                    for iz in 0..grid.nz {
                        if in_abs_c(&cpml.z, iz) {
                            continue;
                        }
                        let g = grid.idx(ix, iy, iz);
                        if pec_lz[g] <= 0.0 || inv_ez[g] <= 0.0 {
                            continue;
                        }
                        let Some(&slot) = dom_to_slot.get(&mats.dom_ez[g]) else { continue };
                        let z = grid.cz[iz];
                        let qpm_sign = domains[slot as usize].qpm.as_ref().map_or(1.0, |q| q.poling_sign_at([x, y, z]));
                        edges_z.push(NonlinearEdgeEntry {
                            idx: g,
                            ix: ix as u32,
                            iy: iy as u32,
                            iz: iz as u32,
                            dom_slot: slot,
                            qpm_sign,
                        });
                    }
                }
            }
        }

        let p_nl_prev = [vec![0.0; edges_x.len()], vec![0.0; edges_y.len()], vec![0.0; edges_z.len()]];
        let n_carrier = if has_carriers {
            [vec![0.0; edges_x.len()], vec![0.0; edges_y.len()], vec![0.0; edges_z.len()]]
        } else {
            [Vec::new(), Vec::new(), Vec::new()]
        };

        Ok(Self { domains, edges: [edges_x, edges_y, edges_z], p_nl_prev, n_carrier, has_carriers })
    }
}

/// Complete 3-component ADE + nonlinear state bank (`Ex`, `Ey`, `Ez`), plus precomputed Ampere `(ca, cb)` update slabs!
#[derive(Clone, Debug)]
pub struct AdeAndAmpereCoeffs {
    pub ca_x: Vec<f64>,
    pub cb_x: Vec<f64>,
    pub ca_y: Vec<f64>,
    pub cb_y: Vec<f64>,
    pub ca_z: Vec<f64>,
    pub cb_z: Vec<f64>,
    pub da_x: Vec<f64>,
    pub db_x: Vec<f64>,
    pub da_y: Vec<f64>,
    pub db_y: Vec<f64>,
    pub da_z: Vec<f64>,
    pub db_z: Vec<f64>,
    pub ade_x: ComponentAdeBank,
    pub ade_y: ComponentAdeBank,
    pub ade_z: ComponentAdeBank,
    pub nonlinear: NonlinearFdtdBank,
}

impl AdeAndAmpereCoeffs {
    pub fn build(
        grid: &YeeGrid,
        mats: &YeeMaterialCoeffs,
        cpml: &CpmlState,
        disp_specs: &BTreeMap<i32, DispersiveDomainSpec>,
        dt: f64,
    ) -> Self {
        let n = grid.slab_len();
        let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
        let ade_x = ComponentAdeBank::build(&mats.dom_ex, &mats.disp_frac_x, disp_specs, dt);
        let ade_y = ComponentAdeBank::build(&mats.dom_ey, &mats.disp_frac_y, disp_specs, dt);
        let ade_z = ComponentAdeBank::build(&mats.dom_ez, &mats.disp_frac_z, disp_specs, dt);

        let mut drude_map_x = vec![0.0; n];
        for (&idx, &d_add) in ade_x.edge_indices.iter().zip(&ade_x.drude_cond_add) {
            drude_map_x[idx] = d_add;
        }
        let mut drude_map_y = vec![0.0; n];
        for (&idx, &d_add) in ade_y.edge_indices.iter().zip(&ade_y.drude_cond_add) {
            drude_map_y[idx] = d_add;
        }
        let mut drude_map_z = vec![0.0; n];
        for (&idx, &d_add) in ade_z.edge_indices.iter().zip(&ade_z.drude_cond_add) {
            drude_map_z[idx] = d_add;
        }

        let make_e = |inv_eps: &[f64], sig: &[f64], pec_l: &[f64], d_map: &[f64], pos_kind: usize| {
            let mut ca = vec![0.0; n];
            let mut cb = vec![0.0; n];
            for ix in 0..=grid.nx {
                let s_abs_x = if pos_kind == 0 {
                    if ix < grid.nx {
                        0.5 * (cpml.x.absorber_sigma_e[ix] + cpml.x.absorber_sigma_e[ix + 1])
                    } else {
                        0.0
                    }
                } else {
                    cpml.x.absorber_sigma_e[ix]
                };
                for iy in 0..=grid.ny {
                    let s_abs_y = if is_2d {
                        0.0
                    } else if pos_kind == 1 {
                        if iy < grid.ny {
                            0.5 * (cpml.y.absorber_sigma_e[iy] + cpml.y.absorber_sigma_e[iy + 1])
                        } else {
                            0.0
                        }
                    } else {
                        cpml.y.absorber_sigma_e[iy]
                    };
                    for iz in 0..=grid.nz {
                        let s_abs_z = if pos_kind == 2 {
                            if iz < grid.nz {
                                0.5 * (cpml.z.absorber_sigma_e[iz] + cpml.z.absorber_sigma_e[iz + 1])
                            } else {
                                0.0
                            }
                        } else {
                            cpml.z.absorber_sigma_e[iz]
                        };
                        let g = grid.idx(ix, iy, iz);
                        if pec_l[g] <= 0.0 || inv_eps[g] <= 0.0 {
                            ca[g] = 0.0;
                            cb[g] = 0.0;
                            continue;
                        }
                        let eps_abs = EPS0 / inv_eps[g];
                        // Scale adiabatic absorber conductivity with local relative permittivity for exact impedance matching!
                        let eps_r_loc = (1.0 / inv_eps[g]).max(1.0);
                        let sig_tot = sig[g] + eps_r_loc * (s_abs_x + s_abs_y + s_abs_z);
                        let half_loss = 0.5 * sig_tot * dt + 0.5 * d_map[g] * dt;
                        let denom = eps_abs + half_loss;
                        ca[g] = (eps_abs - half_loss) / denom;
                        cb[g] = dt / denom;
                    }
                }
            }
            (ca, cb)
        };

        let make_h = |inv_mu: &[f64], face_axis: usize| {
            let mut da = vec![1.0; n];
            let mut db = vec![0.0; n];
            for ix in 0..=grid.nx {
                let sm_x = if face_axis == 0 {
                    if ix == 0 {
                        cpml.x.absorber_sigma_m_over_mu0[0]
                    } else if ix >= grid.nx {
                        cpml.x.absorber_sigma_m_over_mu0[grid.nx - 1]
                    } else {
                        0.5 * (cpml.x.absorber_sigma_m_over_mu0[ix - 1] + cpml.x.absorber_sigma_m_over_mu0[ix])
                    }
                } else if ix < grid.nx {
                    cpml.x.absorber_sigma_m_over_mu0[ix]
                } else {
                    0.0
                };
                for iy in 0..=grid.ny {
                    let sm_y = if is_2d {
                        0.0
                    } else if face_axis == 1 {
                        if iy == 0 {
                            cpml.y.absorber_sigma_m_over_mu0[0]
                        } else if iy >= grid.ny {
                            cpml.y.absorber_sigma_m_over_mu0[grid.ny - 1]
                        } else {
                            0.5 * (cpml.y.absorber_sigma_m_over_mu0[iy - 1] + cpml.y.absorber_sigma_m_over_mu0[iy])
                        }
                    } else if iy < grid.ny {
                        cpml.y.absorber_sigma_m_over_mu0[iy]
                    } else {
                        0.0
                    };
                    for iz in 0..=grid.nz {
                        let sm_z = if face_axis == 2 {
                            if iz == 0 {
                                cpml.z.absorber_sigma_m_over_mu0[0]
                            } else if iz >= grid.nz {
                                cpml.z.absorber_sigma_m_over_mu0[grid.nz - 1]
                            } else {
                                0.5 * (cpml.z.absorber_sigma_m_over_mu0[iz - 1] + cpml.z.absorber_sigma_m_over_mu0[iz])
                            }
                        } else if iz < grid.nz {
                            cpml.z.absorber_sigma_m_over_mu0[iz]
                        } else {
                            0.0
                        };
                        let g = grid.idx(ix, iy, iz);
                        if inv_mu[g] <= 0.0 {
                            da[g] = 0.0;
                            db[g] = 0.0;
                            continue;
                        }
                        let mu_abs = MU0 / inv_mu[g];
                        let mu_r_loc = (1.0 / inv_mu[g]).max(1.0);
                        let sm_tot = MU0 * mu_r_loc * (sm_x + sm_y + sm_z);
                        let half_m = 0.5 * sm_tot * dt;
                        let denom = mu_abs + half_m;
                        da[g] = (mu_abs - half_m) / denom;
                        db[g] = dt / denom;
                    }
                }
            }
            (da, db)
        };

        let (ca_x, cb_x) =
            make_e(mats.inv_eps_x.as_slice(), mats.sigma_ex.as_slice(), mats.pec_len_x.as_slice(), &drude_map_x, 0);
        let (ca_y, cb_y) =
            make_e(mats.inv_eps_y.as_slice(), mats.sigma_ey.as_slice(), mats.pec_len_y.as_slice(), &drude_map_y, 1);
        let (ca_z, cb_z) =
            make_e(mats.inv_eps_z.as_slice(), mats.sigma_ez.as_slice(), mats.pec_len_z.as_slice(), &drude_map_z, 2);

        let (da_x, db_x) = make_h(mats.inv_mu_x.as_slice(), 0);
        let (da_y, db_y) = make_h(mats.inv_mu_y.as_slice(), 1);
        let (da_z, db_z) = make_h(mats.inv_mu_z.as_slice(), 2);

        Self {
            ca_x,
            cb_x,
            ca_y,
            cb_y,
            ca_z,
            cb_z,
            da_x,
            db_x,
            da_y,
            db_y,
            da_z,
            db_z,
            ade_x,
            ade_y,
            ade_z,
            nonlinear: NonlinearFdtdBank::default(),
        }
    }
}

/// Applies the ADE polarization/current source term BEFORE/AFTER updating `e_slab` from `e_n` to `e_{n+1}`.
///
/// Specifically, given `e_old = E^n` on dispersive edges before standard Ampere update (`E_tent = ca * E^n + cb * curl_H`):
/// - For Lorentz poles: computes `P^{n+1} = c1 P^n + c2 P^{n-1} + c3 E^n` and subtracts `(cb / dt) * (P^{n+1} - P^n)`.
/// - For Drude poles: subtracts `cb * 0.5 * (1 + k_d) * J^n`, then once `E^{n+1}` is finalized, advances
///   `J^{n+1} = k_d J^n + beta_d (E^{n+1} + E^n)`.
pub fn apply_ade_correction(bank: &mut ComponentAdeBank, e_slab: &mut [f64], e_old_saved: &[f64], cb: &[f64], dt: f64) {
    if bank.is_empty() {
        return;
    }
    let inv_dt = 1.0 / dt;
    for (k, &g) in bank.edge_indices.iter().enumerate() {
        let e_n = e_old_saved[k];
        let p_range = bank.pole_offsets[k]..bank.pole_offsets[k + 1];
        let mut rhs_sub = 0.0;
        for m in p_range.clone() {
            match bank.poles[m] {
                CompiledPole::Drude { k_d, .. } => {
                    rhs_sub += 0.5 * (1.0 + k_d) * bank.s1[m];
                }
                CompiledPole::Lorentz { c1, c2, c3 } => {
                    let p_n = bank.s1[m];
                    let p_nm1 = bank.s2[m];
                    let p_np1 = c1 * p_n + c2 * p_nm1 + c3 * e_n;
                    rhs_sub += (p_np1 - p_n) * inv_dt;
                    bank.s2[m] = p_n;
                    bank.s1[m] = p_np1;
                }
            }
        }
        let e_np1 = e_slab[g] - cb[g] * rhs_sub;
        e_slab[g] = e_np1;
        // Finalize Drude currents now that E^{n+1} is known
        for m in p_range {
            if let CompiledPole::Drude { k_d, beta_d } = bank.poles[m] {
                bank.s1[m] = k_d * bank.s1[m] + beta_d * (e_np1 + e_n);
            }
        }
    }
}
