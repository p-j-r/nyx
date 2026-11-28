// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! High-order curved and mixed-cell Discontinuous Galerkin Time-Domain solver (`Study::Dgtd`, `"dgtd"`, ADR 144)
//! with explicit power-of-two Local Time-Stepping (LTS).
//!
//! Structured in the three stages of ADR 9:
//! 1. [`assemble_dgtd`]: builds high-order nodal/GLL polynomial bases ($p \in \{1, 2, 3\}$) across `Tri`, `Quad`,
//!    `Tet`, `CurvedTet`, `Prism`, and `Hex` cells, pre-inverts block-local element mass matrices $M_K^{-1}$
//!    (diagonal GLL quadrature via [`femcore::linalg::matfree::Gll`] on `Quad`/`Hex` and local dense factorization
//!    on `Tri`/`Tet`/`CurvedTet`/`Prism` — zero global linear solves), forms Summation-By-Parts (SBP) skew-symmetric
//!    strong-form curl operators $D^k_K = M_K^{-1} \hat{K}^k_K$, resolves 2D/3D waveguide modal ports via
//!    [`crate::modes`], attaches penetrable zero-thickness [`Sheet2dBc`] (`FR-35`) interface jump tensors
//!    $\hat{\mathbf{n}} \times (\mathbf{H}^+ - \mathbf{H}^-) = \boldsymbol{\sigma}_{2\text{D}} \mathbf{E}_t^*$,
//!    and partitions elements into $L$ power-of-two CFL tiers ($\Delta t_\ell = 2^\ell \Delta t_0$).
//! 2. [`solve_dgtd`]: advances $(\mathbf{E}, \mathbf{H})$ and dispersive Drude/Lorentz/gyrotropic ADE states in time
//!    using conservative multi-rate hierarchical 4th-order Local Time-Stepping (`Lts`), 5-stage 4th-order low-storage
//!    Runge–Kutta (`Lserk4`), or classical 4-stage `Rk4`, with fixed-chunk deterministic parallel element updates via
//!    [`femcore::par`] (`NYX_THREADS=1` vs `8` bit-identical) and on-the-fly running DFT / time monitors.
//! 3. [`dgtd_fields`]: extracts broadband modal scattering parameters into [`DrivenPortData`], packages time-domain
//!    and LTS diagnostics into [`DgtdData`] (`DGTD_DATA_KEY`), and projects frequency-domain DFT or final time-domain
//!    electromagnetic fields onto the output [`SolveResult`].

use crate::excitation::{Excitation, ExcitationSpec};
use crate::ibc::{sheet2d_tangent_frame_3d, sorted3, EPS0, ETA0, MU0};
use crate::modes::{assemble_modes_with_bcs, solve_mode_pairs, EpsCell, C0};
use crate::problem::{symmetry_facet_flags_mixed, Problem, Study};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, CurvedTet};
use femcore::fem_mixed::{
    cell_corner_ref_coords, facet_lagrange_shape, facet_rule_mixed, lagrange_shape_into, quad_rule_mixed,
    transfinite_quad_center, IsoparametricCell, LagrangeSpace,
};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, C};
use femcore::linalg::matfree::Gll;
use femcore::linalg::prom::s_to_z_and_y_checked;
use femcore::mesh::cells::{canonical_quad_face, CellType};
use femcore::mesh::mixed::{
    cell_local_edges, cell_local_quad_faces, cell_local_tri_faces, sorted_tri_face, FacetType, MixedDiscretisation,
    QUAD_EDGES, TET_EDGES, TET_TRI_FACES, TRI_EDGES,
};
use femcore::mesh::SimplexMesh;
use femcore::props::PmlDir;
use nyx_physics::decode::Params;
use nyx_physics::{
    cell_fields_of, discretise, discretise_mixed, Ctx, DrivenPortData, Error, Field, FieldFrame, Mesh,
    MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::f64::consts::PI;
use std::time::Instant;

/// Scalar key carrying serialized [`DgtdData`] on `ems` `dgtd` results (`ADR 144`).
pub const DGTD_DATA_KEY: &str = "dgtd time-domain data";

/// Allowed choices for `flux` on `StudyDesc` `"dgtd"`.
pub const DGTD_FLUX_CHOICES: [&str; 2] = ["upwind", "central"];

/// Allowed choices for `time_stepper` on `StudyDesc` `"dgtd"`.
pub const DGTD_STEPPER_CHOICES: [&str; 3] = ["lts", "lserk4", "rk4"];

/// Allowed choices for `polarization` on `StudyDesc` `"dgtd"`.
pub const DGTD_POL_CHOICES: [&str; 3] = ["full", "tm", "te"];

/// Allowed choices for `pulse` on `StudyDesc` `"dgtd"`.
pub const DGTD_PULSE_CHOICES: [&str; 2] = ["gaussian", "ricker"];

/// Riemann numerical flux across element interfaces.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DgtdFlux {
    /// Full Rankine–Hugoniot upwind numerical flux ($\alpha = 1$).
    #[default]
    Upwind,
    /// Energy-conserving central numerical flux ($\alpha = 0$).
    Central,
}

impl DgtdFlux {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "upwind" | "" => Ok(Self::Upwind),
            "central" => Ok(Self::Central),
            other => Err(format!("dgtd study: unknown flux '{other}' (expected 'upwind' or 'central')")),
        }
    }

    #[inline]
    pub fn alpha(self) -> f64 {
        match self {
            Self::Upwind => 1.0,
            Self::Central => 0.0,
        }
    }
}

/// Explicit time integration scheme for `Study::Dgtd`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DgtdTimeStepper {
    /// Conservative power-of-two multi-rate Local Time-Stepping ($\Delta t_\ell = 2^\ell \Delta t_0$).
    #[default]
    Lts,
    /// Uniform 5-stage 4th-order Carpenter–Kennedy Low-Storage Explicit Runge–Kutta.
    Lserk4,
    /// Uniform classical 4-stage 4th-order Explicit Runge–Kutta.
    Rk4,
}

impl DgtdTimeStepper {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "lts" | "" => Ok(Self::Lts),
            "lserk4" => Ok(Self::Lserk4),
            "rk4" => Ok(Self::Rk4),
            other => Err(format!("dgtd study: unknown time_stepper '{other}' (expected 'lts', 'lserk4', or 'rk4')")),
        }
    }
}

/// Active electromagnetic field polarization in 2D or 3D.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DgtdPolarization {
    /// Full 6-component $(\mathbf{E}, \mathbf{H})$ vector Maxwell evolution.
    #[default]
    Full,
    /// 2D Transverse Magnetic ($E_z$ out-of-plane, $H_x, H_y$ in-plane).
    Tm,
    /// 2D Transverse Electric ($E_x, E_y$ in-plane, $H_z$ out-of-plane).
    Te,
}

impl DgtdPolarization {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "full" | "" => Ok(Self::Full),
            "tm" => Ok(Self::Tm),
            "te" => Ok(Self::Te),
            other => Err(format!("dgtd study: unknown polarization '{other}' (expected 'full', 'tm', or 'te')")),
        }
    }
}

/// Time-domain source pulse waveform $s(t)$.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DgtdPulse {
    /// Broadband Gaussian-modulated cosine carrier:
    /// $s(t) = A \exp(-(t - t_0)^2 / (2\sigma_t^2)) \cos(2\pi f_0 (t - t_0) + \phi)$ with $\sigma_t = 1 / (2\pi f_{\text{width}})$.
    GaussianModulated {
        f0_hz: f64,
        fwidth_hz: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        t0_s: Option<f64>,
        amplitude: f64,
        phase_rad: f64,
    },
    /// Zero-DC Ricker wavelet (second derivative of a Gaussian):
    /// $s(t) = A (1 - 2\pi^2 f_0^2 (t - t_0)^2) \exp(-\pi^2 f_0^2 (t - t_0)^2)$.
    Ricker {
        f0_hz: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        t0_s: Option<f64>,
        amplitude: f64,
    },
}

impl Default for DgtdPulse {
    fn default() -> Self {
        Self::GaussianModulated { f0_hz: 1.934e14, fwidth_hz: 4.5e13, t0_s: None, amplitude: 1.0, phase_rad: 0.0 }
    }
}

impl DgtdPulse {
    /// Center frequency $f_0$ in Hz.
    #[inline]
    pub fn f0_hz(&self) -> f64 {
        match *self {
            Self::GaussianModulated { f0_hz, .. } | Self::Ricker { f0_hz, .. } => f0_hz,
        }
    }

    /// Bandwidth parameter $f_{\text{width}}$ in Hz.
    #[inline]
    pub fn fwidth_hz(&self) -> f64 {
        match *self {
            Self::GaussianModulated { fwidth_hz, .. } => fwidth_hz,
            Self::Ricker { f0_hz, .. } => 0.5 * f0_hz,
        }
    }

    /// Pulse peak delay time $t_0$ in seconds.
    #[inline]
    pub fn t0_s(&self) -> f64 {
        match *self {
            Self::GaussianModulated { fwidth_hz, t0_s, .. } => {
                t0_s.unwrap_or_else(|| 4.0 / (2.0 * PI * fwidth_hz.max(1.0)))
            }
            Self::Ricker { f0_hz, t0_s, .. } => t0_s.unwrap_or_else(|| 2.0 / f0_hz.max(1.0)),
        }
    }

    /// Evaluates $s(t)$ at time $t$ (seconds).
    #[inline]
    pub fn eval(&self, t: f64) -> f64 {
        let dt = t - self.t0_s();
        match *self {
            Self::GaussianModulated { f0_hz, fwidth_hz, amplitude, phase_rad, .. } => {
                let sigma = 1.0 / (2.0 * PI * fwidth_hz.max(1.0));
                let arg = dt / sigma;
                if arg.abs() > 14.0 {
                    0.0
                } else {
                    amplitude * (-0.5 * arg * arg).exp() * (2.0 * PI * f0_hz * dt + phase_rad).cos()
                }
            }
            Self::Ricker { f0_hz, amplitude, .. } => {
                let u = PI * f0_hz * dt;
                let u2 = u * u;
                if u2 > 100.0 {
                    0.0
                } else {
                    amplitude * (1.0 - 2.0 * u2) * (-u2).exp()
                }
            }
        }
    }

    pub fn is_finite(&self) -> bool {
        match *self {
            Self::GaussianModulated { f0_hz, fwidth_hz, t0_s, amplitude, phase_rad } => {
                f0_hz.is_finite()
                    && fwidth_hz.is_finite()
                    && t0_s.is_none_or(f64::is_finite)
                    && amplitude.is_finite()
                    && phase_rad.is_finite()
            }
            Self::Ricker { f0_hz, t0_s, amplitude } => {
                f0_hz.is_finite() && t0_s.is_none_or(f64::is_finite) && amplitude.is_finite()
            }
        }
    }
}

/// Auxiliary Differential Equation (ADE) dispersive polarization pole for Drude, Lorentz, or gyrotropic media.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DgtdDispersionPole {
    /// Drude pole: $\partial_t \mathbf{J}_D + \gamma \mathbf{J}_D = \varepsilon_0 \omega_p^2 \mathbf{E}$.
    Drude { omega_p_rad_s: f64, gamma_rad_s: f64 },
    /// Lorentz oscillator pole:
    /// $\partial_t \mathbf{P}_L = \mathbf{J}_L$, $\partial_t \mathbf{J}_L + \gamma \mathbf{J}_L + \omega_0^2 \mathbf{P}_L = \varepsilon_0 \Delta\varepsilon\,\omega_0^2 \mathbf{E}$.
    Lorentz { omega_0_rad_s: f64, gamma_rad_s: f64, delta_eps: f64 },
    /// Magneto-optic gyrotropic cyclotron-precession pole (`FR-34`):
    /// $\partial_t \mathbf{J}_G + \gamma \mathbf{J}_G + \boldsymbol{\omega}_c \times \mathbf{J}_G = \varepsilon_0 \omega_p^2 \mathbf{E}$.
    GyrotropicDrude { omega_p_rad_s: f64, gamma_rad_s: f64, omega_c_rad_s: [f64; 3] },
}

impl DgtdDispersionPole {
    pub fn is_finite(&self) -> bool {
        match *self {
            Self::Drude { omega_p_rad_s, gamma_rad_s } => omega_p_rad_s.is_finite() && gamma_rad_s.is_finite(),
            Self::Lorentz { omega_0_rad_s, gamma_rad_s, delta_eps } => {
                omega_0_rad_s.is_finite() && gamma_rad_s.is_finite() && delta_eps.is_finite()
            }
            Self::GyrotropicDrude { omega_p_rad_s, gamma_rad_s, omega_c_rad_s } => {
                omega_p_rad_s.is_finite() && gamma_rad_s.is_finite() && omega_c_rad_s.iter().all(|v| v.is_finite())
            }
        }
    }

    #[inline]
    pub fn characteristic_rate(&self) -> f64 {
        match *self {
            Self::Drude { omega_p_rad_s, gamma_rad_s } => omega_p_rad_s.abs().max(gamma_rad_s.abs()),
            Self::Lorentz { omega_0_rad_s, gamma_rad_s, delta_eps } => {
                (omega_0_rad_s.abs() * delta_eps.abs().max(1.0).sqrt()).max(gamma_rad_s.abs())
            }
            Self::GyrotropicDrude { omega_p_rad_s, gamma_rad_s, omega_c_rad_s } => {
                let wc = (omega_c_rad_s[0].powi(2) + omega_c_rad_s[1].powi(2) + omega_c_rad_s[2].powi(2)).sqrt();
                omega_p_rad_s.abs().max(gamma_rad_s.abs()).max(wc)
            }
        }
    }
}

/// Specification of a surface Poynting flux monitor (`FluxMonitor`) across boundary or interior facets.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DgtdFluxMonitorSpec {
    pub name: String,
    pub labels: Vec<i32>,
}

/// Complete configuration for `Study::Dgtd`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DgtdConfig {
    /// Total physical simulation duration in seconds.
    pub t_final: f64,
    /// Courant–Friedrichs–Lewy stability multiplier (`0 < cfl <= 1.2`, default `0.8`).
    pub cfl: f64,
    /// Numerical Riemann flux (`Upwind` or `Central`).
    pub flux: DgtdFlux,
    /// Explicit time-stepping algorithm (`Lts`, `Lserk4`, or `Rk4`).
    pub time_stepper: DgtdTimeStepper,
    /// Maximum number of power-of-two LTS tiers $L \in [1, 12]$ (`1` forces uniform stepping).
    pub max_lts_tiers: usize,
    /// Polarization mode (`Full`, `Tm`, or `Te`).
    pub polarization: DgtdPolarization,
    /// Optional explicit polynomial basis order override $p \in \{1, 2, 3\}$ (`None` uses `problem.order`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub basis_order: Option<usize>,
    /// Pulse time waveform $s(t)$.
    pub pulse: DgtdPulse,
    /// Monitor frequencies in Hz for running DFT `FieldMonitor`, `FluxMonitor`, and `DrivenPortData` extraction.
    pub freqs_hz: Vec<f64>,
    /// 1-based index of the excited port (`1` by default).
    pub excited_port: usize,
    /// Boundary labels treated as first-order Silver–Müller absorbing boundaries.
    #[serde(default)]
    pub absorbing_labels: Vec<i32>,
    /// Probe point coordinates `[x, y, z]` (in 2D `[x, y, 0.0]`) for `FieldTimeMonitor`.
    #[serde(default)]
    pub probes: Vec<[f64; 3]>,
    /// Surface Poynting flux monitors.
    #[serde(default)]
    pub flux_monitors: Vec<DgtdFluxMonitorSpec>,
    /// Optional per-domain dispersive Drude/Lorentz/gyrotropic ADE poles.
    #[serde(default)]
    pub dispersion_poles: BTreeMap<i32, Vec<DgtdDispersionPole>>,
    /// Optional analytical cavity eigenmode indices `[m, n, p]` to initialize at $t = 0$ (for closed-cavity energy tests).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub initial_cavity_mode: Option<[usize; 3]>,
    /// Compute execution backend (`"auto"`, `"cpu"`, or `"gpu"`, `ADR 146`).
    #[serde(default)]
    pub backend: crate::problem::TdBackend,
    /// GPU floating-point precision (`"auto"`, `"f32"`, or `"f64"`, `ADR 146`).
    #[serde(default)]
    pub precision: crate::problem::TdPrecision,
    /// Opt-in activation for material `[nonlinear]` $\chi^{(2)}/\chi^{(3)}$ polarization (`false` auto-skips library blocks).
    #[serde(default)]
    pub nonlinear: bool,
    /// Minimum total nodal DOFs for `backend = "auto"` to offload uniform-step DGTD to a discrete GPU (`50_000` default).
    #[serde(default = "default_dgtd_gpu_min_dofs")]
    pub gpu_min_dofs: usize,
}

fn default_dgtd_gpu_min_dofs() -> usize {
    50_000
}

impl Default for DgtdConfig {
    fn default() -> Self {
        Self {
            t_final: 100e-15,
            cfl: 0.8,
            flux: DgtdFlux::Upwind,
            time_stepper: DgtdTimeStepper::Lts,
            max_lts_tiers: 6,
            polarization: DgtdPolarization::Full,
            basis_order: None,
            pulse: DgtdPulse::default(),
            freqs_hz: vec![1.934e14],
            excited_port: 1,
            absorbing_labels: Vec::new(),
            probes: Vec::new(),
            flux_monitors: Vec::new(),
            dispersion_poles: BTreeMap::new(),
            initial_cavity_mode: None,
            backend: crate::problem::TdBackend::Auto,
            precision: crate::problem::TdPrecision::Auto,
            nonlinear: false,
            gpu_min_dofs: default_dgtd_gpu_min_dofs(),
        }
    }
}

impl DgtdConfig {
    pub fn is_finite(&self) -> bool {
        self.t_final.is_finite()
            && self.cfl.is_finite()
            && self.pulse.is_finite()
            && self.freqs_hz.iter().all(|f| f.is_finite())
            && self.probes.iter().flatten().all(|x| x.is_finite())
            && self.dispersion_poles.values().flatten().all(DgtdDispersionPole::is_finite)
    }

    pub fn reference_freq(&self) -> f64 {
        if !self.freqs_hz.is_empty() {
            0.5 * (self.freqs_hz[0] + self.freqs_hz[self.freqs_hz.len() - 1])
        } else {
            self.pulse.f0_hz()
        }
    }
}

/// Decodes a [`DgtdConfig`] from `ProblemSpec` study parameters.
pub fn decode_dgtd(sp: &Params, ctx: &Ctx) -> Result<DgtdConfig, Error> {
    let t_final = sp.number("t_final")?;
    if !(t_final > 0.0 && t_final.is_finite()) {
        return Err(Error::invalid(format!("dgtd study requires positive finite t_final (got {t_final})")));
    }
    let cfl = sp.number("cfl")?;
    if !(cfl > 0.0 && cfl <= 2.0 && cfl.is_finite()) {
        return Err(Error::invalid(format!("dgtd study requires 0 < cfl <= 2.0 (got {cfl})")));
    }
    let flux = DgtdFlux::parse(&sp.text("flux")?).map_err(Error::invalid)?;
    let time_stepper = DgtdTimeStepper::parse(&sp.text("time_stepper")?).map_err(Error::invalid)?;
    let max_lts_raw = sp.integer("max_lts_tiers")?;
    if !(1..=12).contains(&max_lts_raw) {
        return Err(Error::invalid(format!("dgtd study requires 1 <= max_lts_tiers <= 12 (got {max_lts_raw})")));
    }
    let max_lts_tiers = max_lts_raw as usize;
    let polarization = DgtdPolarization::parse(&sp.text("polarization")?).map_err(Error::invalid)?;
    let basis_order = match sp.integer("basis_order")? {
        0 => None,
        p if (1..=3).contains(&p) => Some(p as usize),
        other => {
            return Err(Error::invalid(format!(
                "dgtd study requires basis_order in 0..=3 (0 uses problem.order; got {other})"
            )))
        }
    };
    let f0 = sp.number("f0")?;
    if !(f0 > 0.0 && f0.is_finite()) {
        return Err(Error::invalid(format!("dgtd study requires positive finite f0 (got {f0})")));
    }
    let fwidth_raw = sp.number("fwidth")?;
    if !fwidth_raw.is_finite() || fwidth_raw < 0.0 {
        return Err(Error::invalid(format!("dgtd study requires finite fwidth >= 0 (got {fwidth_raw})")));
    }
    let fwidth = if fwidth_raw > 0.0 { fwidth_raw } else { 0.25 * f0 };
    let pulse_amp = if sp.has("pulse_amplitude") { sp.number("pulse_amplitude")? } else { 1.0 };
    if !pulse_amp.is_finite() || pulse_amp < 0.0 {
        return Err(Error::invalid(format!(
            "dgtd study requires finite non-negative pulse_amplitude (got {pulse_amp})"
        )));
    }
    let pulse_kind = sp.text("pulse")?;
    let pulse = match pulse_kind.as_str() {
        "gaussian" | "" => DgtdPulse::GaussianModulated {
            f0_hz: f0,
            fwidth_hz: fwidth,
            t0_s: None,
            amplitude: pulse_amp,
            phase_rad: 0.0,
        },
        "ricker" => DgtdPulse::Ricker { f0_hz: f0, t0_s: None, amplitude: pulse_amp },
        other => {
            return Err(Error::invalid(format!(
                "dgtd study: unknown pulse '{other}' (expected 'gaussian' or 'ricker')"
            )))
        }
    };

    let explicit_freqs = sp.optional_numbers("freqs_hz")?;
    let freqs_hz = if let Some(fs) = explicit_freqs.filter(|v| !v.is_empty()) {
        for &f in &fs {
            if !(f > 0.0 && f.is_finite()) {
                return Err(Error::invalid(format!(
                    "dgtd study requires all entries of 'freqs_hz' to be positive finite numbers (got {f})"
                )));
            }
        }
        fs
    } else {
        let f_min = sp.number("f_min")?;
        let f_max = sp.number("f_max")?;
        let n_freqs = sp.integer("n_freqs")?.max(1) as usize;
        if f_min > 0.0 && f_max >= f_min {
            if n_freqs == 1 {
                vec![0.5 * (f_min + f_max)]
            } else {
                (0..n_freqs).map(|i| f_min + (f_max - f_min) * (i as f64) / ((n_freqs - 1) as f64)).collect()
            }
        } else {
            vec![f0]
        }
    };
    let excited_port = sp.integer("excited_port")?.max(1) as usize;
    let abs_sel = sp.text("absorbing_on")?;
    let mut absorbing_labels = Vec::new();
    if !abs_sel.trim().is_empty() {
        for part in abs_sel.split(',') {
            let s = part.trim();
            if !s.is_empty() {
                absorbing_labels.extend(ctx.resolver.boundaries(s)?);
            }
        }
        absorbing_labels.sort_unstable();
        absorbing_labels.dedup();
    }

    let mut flux_monitors = Vec::new();
    if sp.has("flux_monitors_on") {
        let fm_sel = sp.text("flux_monitors_on")?;
        if !fm_sel.trim().is_empty() {
            for part in fm_sel.split(',') {
                let s = part.trim();
                if !s.is_empty() {
                    let mut labels = ctx.resolver.boundaries(s)?;
                    labels.sort_unstable();
                    labels.dedup();
                    flux_monitors.push(DgtdFluxMonitorSpec { name: s.to_string(), labels });
                }
            }
        }
    }

    let backend = if sp.has("backend") {
        crate::problem::TdBackend::parse(&sp.text("backend")?).map_err(Error::invalid)?
    } else {
        crate::problem::TdBackend::Auto
    };
    let precision = if sp.has("precision") {
        crate::problem::TdPrecision::parse(&sp.text("precision")?).map_err(Error::invalid)?
    } else {
        crate::problem::TdPrecision::Auto
    };
    let nonlinear = sp.boolean("nonlinear")?;

    Ok(DgtdConfig {
        t_final,
        cfl,
        flux,
        time_stepper,
        max_lts_tiers,
        polarization,
        basis_order,
        pulse,
        freqs_hz,
        excited_port,
        absorbing_labels,
        probes: Vec::new(),
        flux_monitors,
        dispersion_poles: BTreeMap::new(),
        initial_cavity_mode: None,
        backend,
        precision,
        nonlinear,
        gpu_min_dofs: default_dgtd_gpu_min_dofs(),
    })
}

/// Frequency-domain Poynting flux spectrum recorded by a [`DgtdFluxMonitorSpec`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DgtdFluxSpectrum {
    pub name: String,
    pub freqs_hz: Vec<f64>,
    /// Unnormalized Parseval spectral energy flux $\frac{1}{2}\operatorname{Re}\int_\Gamma (\hat{\mathbf{E}}(\omega) \times \hat{\mathbf{H}}^*(\omega))\cdot\hat{\mathbf{n}}\,dS$ ($\text{W}\cdot\text{s}^2$).
    #[serde(default)]
    pub raw_flux_w_s2: Vec<f64>,
    /// Time-averaged real Poynting flux $P(\omega) = \frac{1}{2}\operatorname{Re}\int_\Gamma (\mathbf{E}(\omega) \times \mathbf{H}^*(\omega))\cdot\hat{\mathbf{n}}\,dS$
    /// normalized by the incident pulse spectrum $|\hat{s}(\omega)|^2$ (clamped to carrier peak $|\hat{s}_{\text{peak}}|^2$ for out-of-band nonlinear harmonics).
    pub flux_w: Vec<f64>,
    /// Complex tangential electric and magnetic field overlap integrals `[re, im]` normalized by $\hat{s}(\omega)$.
    pub e_avg: Vec<[f64; 2]>,
    pub h_avg: Vec<[f64; 2]>,
}

/// Time-series probe waveform recorded by `FieldTimeMonitor`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DgtdTimeSeries {
    pub position: [f64; 3],
    pub cell_index: usize,
    pub e_t: Vec<[f64; 3]>,
    pub h_t: Vec<[f64; 3]>,
}

/// Serialized time-domain diagnostics, LTS work statistics, energy history, and monitor outputs of `Study::Dgtd`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DgtdData {
    pub basis_order: usize,
    pub n_cells: usize,
    pub total_dofs: usize,
    pub dt_min_s: f64,
    pub dt_max_s: f64,
    pub h_min_m: f64,
    pub h_max_m: f64,
    pub n_lts_tiers: usize,
    pub cells_per_tier: Vec<usize>,
    pub n_macro_steps: usize,
    pub n_finest_steps: usize,
    /// Total element-stage RHS evaluations performed during `solve_dgtd`.
    pub total_element_updates: u64,
    /// Equivalent element-stage RHS evaluations that uniform global $\Delta t_0$ stepping would require.
    pub equivalent_uniform_updates: u64,
    /// Theoretical and measured LTS work speedup factor `equivalent_uniform_updates / total_element_updates`.
    pub lts_work_speedup: f64,
    pub solve_time_s: f64,
    pub times_s: Vec<f64>,
    pub energy_j: Vec<f64>,
    pub pulse_spectrum: Vec<[f64; 2]>,
    pub time_monitors: Vec<DgtdTimeSeries>,
    pub flux_monitors: Vec<DgtdFluxSpectrum>,
}

impl DgtdData {
    pub fn of(result: &SolveResult) -> Option<Self> {
        result.scalars.iter().find(|(n, _)| n == DGTD_DATA_KEY).and_then(|(_, v)| match v {
            Scalar::Text(t) => serde_json::from_str(t).ok(),
            Scalar::Num(_) => None,
        })
    }

    pub fn scalar(&self) -> (String, Scalar) {
        let text = serde_json::to_string(self).unwrap_or_default();
        (DGTD_DATA_KEY.to_string(), Scalar::Text(text))
    }
}

// ---------------------------------------------------------------------------
// High-order polynomial basis (p = 1..=3) on Tri, Quad, Tet, Prism, Hex
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn mat3_mul_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn invert_3x3(m: &[[f64; 3]; 3]) -> Result<[[f64; 3]; 3], String> {
    let c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
    let c01 = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
    let c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
    let c10 = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
    let c11 = m[0][0] * m[2][2] - m[0][2] * m[2][0];
    let c12 = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]);
    let c20 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    let c21 = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]);
    let c22 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    let det = m[0][0] * c00 + m[0][1] * c01 + m[0][2] * c02;
    if !(det.abs() > 1e-30 && det.is_finite()) {
        return Err("singular 3x3 constitutive material tensor in DGTD".into());
    }
    let inv = 1.0 / det;
    Ok([[c00 * inv, c10 * inv, c20 * inv], [c01 * inv, c11 * inv, c21 * inv], [c02 * inv, c12 * inv, c22 * inv]])
}

/// Number of local DG basis nodes on `cell_type` at polynomial degree `p` (`1..=3`).
pub fn dg_nloc(cell_type: CellType, p: usize) -> usize {
    match (cell_type, p) {
        (CellType::Tri, 1) => 3,
        (CellType::Tri, 2) => 6,
        (CellType::Tri, _) => 10,
        (CellType::Quad, _) => (p + 1) * (p + 1),
        (CellType::Tet, 1) => 4,
        (CellType::Tet, 2) => 10,
        (CellType::Tet, _) => 20,
        (CellType::Prism, 1) => 6,
        (CellType::Prism, 2) => 18,
        (CellType::Prism, _) => 40,
        (CellType::Hex, _) => (p + 1) * (p + 1) * (p + 1),
    }
}

/// Reference node coordinates of the degree-`p` DG basis on `cell_type`.
pub fn dg_ref_nodes(cell_type: CellType, p: usize, gll: &Gll) -> Vec<[f64; 3]> {
    match cell_type {
        CellType::Tri => match p {
            1 => vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            2 => vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.5, 0.0, 0.0],
                [0.0, 0.5, 0.0],
                [0.5, 0.5, 0.0],
            ],
            _ => {
                let c = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
                let mut out = c.to_vec();
                for &(i, j) in &TRI_EDGES {
                    out.push([(2.0 * c[i][0] + c[j][0]) / 3.0, (2.0 * c[i][1] + c[j][1]) / 3.0, 0.0]);
                    out.push([(c[i][0] + 2.0 * c[j][0]) / 3.0, (c[i][1] + 2.0 * c[j][1]) / 3.0, 0.0]);
                }
                out.push([1.0 / 3.0, 1.0 / 3.0, 0.0]);
                out
            }
        },
        CellType::Quad => {
            let mut out = Vec::with_capacity((p + 1) * (p + 1));
            for &eta in &gll.nodes {
                for &xi in &gll.nodes {
                    out.push([xi, eta, 0.0]);
                }
            }
            out
        }
        CellType::Tet => match p {
            1 => cell_corner_ref_coords(CellType::Tet).to_vec(),
            2 => {
                let c = cell_corner_ref_coords(CellType::Tet);
                let mut out = c.to_vec();
                for &(i, j) in &TET_EDGES {
                    out.push([0.5 * (c[i][0] + c[j][0]), 0.5 * (c[i][1] + c[j][1]), 0.5 * (c[i][2] + c[j][2])]);
                }
                out
            }
            _ => {
                let c = cell_corner_ref_coords(CellType::Tet);
                let mut out = c.to_vec();
                for &(i, j) in &TET_EDGES {
                    out.push([
                        (2.0 * c[i][0] + c[j][0]) / 3.0,
                        (2.0 * c[i][1] + c[j][1]) / 3.0,
                        (2.0 * c[i][2] + c[j][2]) / 3.0,
                    ]);
                    out.push([
                        (c[i][0] + 2.0 * c[j][0]) / 3.0,
                        (c[i][1] + 2.0 * c[j][1]) / 3.0,
                        (c[i][2] + 2.0 * c[j][2]) / 3.0,
                    ]);
                }
                for &[i, j, k] in &TET_TRI_FACES {
                    out.push([
                        (c[i][0] + c[j][0] + c[k][0]) / 3.0,
                        (c[i][1] + c[j][1] + c[k][1]) / 3.0,
                        (c[i][2] + c[j][2] + c[k][2]) / 3.0,
                    ]);
                }
                out
            }
        },
        CellType::Prism => {
            let tri_nodes = dg_ref_nodes(CellType::Tri, p, gll);
            let mut out = Vec::with_capacity(tri_nodes.len() * (p + 1));
            for &zeta in &gll.nodes {
                for tn in &tri_nodes {
                    out.push([tn[0], tn[1], zeta]);
                }
            }
            out
        }
        CellType::Hex => {
            let mut out = Vec::with_capacity((p + 1) * (p + 1) * (p + 1));
            for &zeta in &gll.nodes {
                for &eta in &gll.nodes {
                    for &xi in &gll.nodes {
                        out.push([xi, eta, zeta]);
                    }
                }
            }
            out
        }
    }
}

/// Evaluates 1D GLL Lagrange basis polynomials $\ell_a(t)$ and derivatives $\ell_a'(t)$ at $t \in [-1, 1]$.
fn gll_1d_basis(gll: &Gll, t: f64) -> (Vec<f64>, Vec<f64>) {
    let n = gll.nodes.len();
    let mut vals = vec![0.0; n];
    let mut dvals = vec![0.0; n];
    for i in 0..n {
        let mut prod = 1.0;
        let mut dsum = 0.0;
        for j in 0..n {
            if i == j {
                continue;
            }
            let denom = gll.nodes[i] - gll.nodes[j];
            prod *= (t - gll.nodes[j]) / denom;
            let mut term = 1.0 / denom;
            for k in 0..n {
                if k != i && k != j {
                    term *= (t - gll.nodes[k]) / (gll.nodes[i] - gll.nodes[k]);
                }
            }
            dsum += term;
        }
        vals[i] = prod;
        dvals[i] = dsum;
    }
    (vals, dvals)
}

/// Evaluates degree-`p` DG basis functions `vals` and reference gradients `dxi` at reference point `xi`.
pub fn dg_shape_into(
    cell_type: CellType,
    p: usize,
    gll: &Gll,
    xi: &[f64],
    vals: &mut Vec<f64>,
    dxi: &mut Vec<[f64; 3]>,
) {
    vals.clear();
    dxi.clear();
    let f3 = |l: f64| -> (f64, f64) { (0.5 * l * (3.0 * l - 1.0) * (3.0 * l - 2.0), 13.5 * l * l - 9.0 * l + 1.0) };
    let f2 = |l: f64| -> (f64, f64) { (4.5 * l * (3.0 * l - 1.0), 27.0 * l - 4.5) };

    match cell_type {
        CellType::Tri => {
            if p == 1 {
                lagrange_shape_into(CellType::Tri, LagrangeSpace::Order1, xi, vals, dxi);
            } else if p == 2 {
                lagrange_shape_into(CellType::Tri, LagrangeSpace::Order2Serendipity, xi, vals, dxi);
            } else {
                let (u, v) = (xi[0], xi[1]);
                let l = [1.0 - u - v, u, v];
                let dl: [[f64; 2]; 3] = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
                for i in 0..3 {
                    let (v3, dv3) = f3(l[i]);
                    vals.push(v3);
                    dxi.push([dv3 * dl[i][0], dv3 * dl[i][1], 0.0]);
                }
                for &(i, j) in &TRI_EDGES {
                    let (vi2, dvi2) = f2(l[i]);
                    vals.push(vi2 * l[j]);
                    dxi.push([dvi2 * l[j] * dl[i][0] + vi2 * dl[j][0], dvi2 * l[j] * dl[i][1] + vi2 * dl[j][1], 0.0]);
                    let (vj2, dvj2) = f2(l[j]);
                    vals.push(l[i] * vj2);
                    dxi.push([dl[i][0] * vj2 + l[i] * dvj2 * dl[j][0], dl[i][1] * vj2 + l[i] * dvj2 * dl[j][1], 0.0]);
                }
                vals.push(27.0 * l[0] * l[1] * l[2]);
                let g0 = 27.0 * l[1] * l[2];
                let g1 = 27.0 * l[0] * l[2];
                let g2 = 27.0 * l[0] * l[1];
                dxi.push([g1 - g0, g2 - g0, 0.0]);
            }
        }
        CellType::Quad => {
            let (lu, dlu) = gll_1d_basis(gll, xi[0]);
            let (lv, dlv) = gll_1d_basis(gll, xi[1]);
            for b in 0..lv.len() {
                for a in 0..lu.len() {
                    vals.push(lu[a] * lv[b]);
                    dxi.push([dlu[a] * lv[b], lu[a] * dlv[b], 0.0]);
                }
            }
        }
        CellType::Tet => {
            if p == 1 {
                lagrange_shape_into(CellType::Tet, LagrangeSpace::Order1, xi, vals, dxi);
            } else if p == 2 {
                lagrange_shape_into(CellType::Tet, LagrangeSpace::Order2Serendipity, xi, vals, dxi);
            } else {
                let (u, v, w) = (xi[0], xi[1], xi[2]);
                let l = [1.0 - u - v - w, u, v, w];
                let dl: [[f64; 3]; 4] = [[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
                for i in 0..4 {
                    let (v3, dv3) = f3(l[i]);
                    vals.push(v3);
                    dxi.push([dv3 * dl[i][0], dv3 * dl[i][1], dv3 * dl[i][2]]);
                }
                for &(i, j) in &TET_EDGES {
                    let (vi2, dvi2) = f2(l[i]);
                    vals.push(vi2 * l[j]);
                    dxi.push([
                        dvi2 * l[j] * dl[i][0] + vi2 * dl[j][0],
                        dvi2 * l[j] * dl[i][1] + vi2 * dl[j][1],
                        dvi2 * l[j] * dl[i][2] + vi2 * dl[j][2],
                    ]);
                    let (vj2, dvj2) = f2(l[j]);
                    vals.push(l[i] * vj2);
                    dxi.push([
                        dl[i][0] * vj2 + l[i] * dvj2 * dl[j][0],
                        dl[i][1] * vj2 + l[i] * dvj2 * dl[j][1],
                        dl[i][2] * vj2 + l[i] * dvj2 * dl[j][2],
                    ]);
                }
                for &[i, j, k] in &TET_TRI_FACES {
                    vals.push(27.0 * l[i] * l[j] * l[k]);
                    let gi = 27.0 * l[j] * l[k];
                    let gj = 27.0 * l[i] * l[k];
                    let gk = 27.0 * l[i] * l[j];
                    dxi.push([
                        gi * dl[i][0] + gj * dl[j][0] + gk * dl[k][0],
                        gi * dl[i][1] + gj * dl[j][1] + gk * dl[k][1],
                        gi * dl[i][2] + gj * dl[j][2] + gk * dl[k][2],
                    ]);
                }
            }
        }
        CellType::Prism => {
            let mut t_vals = Vec::new();
            let mut t_dxi = Vec::new();
            dg_shape_into(CellType::Tri, p, gll, &[xi[0], xi[1], 0.0], &mut t_vals, &mut t_dxi);
            let (lz, dlz) = gll_1d_basis(gll, xi[2]);
            for c_idx in 0..lz.len() {
                for a in 0..t_vals.len() {
                    vals.push(t_vals[a] * lz[c_idx]);
                    dxi.push([t_dxi[a][0] * lz[c_idx], t_dxi[a][1] * lz[c_idx], t_vals[a] * dlz[c_idx]]);
                }
            }
        }
        CellType::Hex => {
            let (lu, dlu) = gll_1d_basis(gll, xi[0]);
            let (lv, dlv) = gll_1d_basis(gll, xi[1]);
            let (lw, dlw) = gll_1d_basis(gll, xi[2]);
            for c_idx in 0..lw.len() {
                for b in 0..lv.len() {
                    for a in 0..lu.len() {
                        vals.push(lu[a] * lv[b] * lw[c_idx]);
                        dxi.push([dlu[a] * lv[b] * lw[c_idx], lu[a] * dlv[b] * lw[c_idx], lu[a] * lv[b] * dlw[c_idx]]);
                    }
                }
            }
        }
    }
}

/// Inverts an $n \times n$ symmetric positive-definite dense matrix via pivoted Gauss-Jordan elimination.
fn invert_dense_matrix(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    for k in 0..n {
        let mut piv_row = k;
        let mut piv_val = m[k * n + k].abs();
        for r in (k + 1)..n {
            let v = m[r * n + k].abs();
            if v > piv_val {
                piv_val = v;
                piv_row = r;
            }
        }
        if !(piv_val > 1e-30 && piv_val.is_finite()) {
            return Err("singular local element mass matrix in DGTD".into());
        }
        if piv_row != k {
            for col in 0..n {
                m.swap(k * n + col, piv_row * n + col);
                inv.swap(k * n + col, piv_row * n + col);
            }
        }
        let diag = m[k * n + k];
        for col in 0..n {
            m[k * n + col] /= diag;
            inv[k * n + col] /= diag;
        }
        for r in 0..n {
            if r == k {
                continue;
            }
            let factor = m[r * n + k];
            if factor.abs() > 0.0 {
                for col in 0..n {
                    m[r * n + col] -= factor * m[k * n + col];
                    inv[r * n + col] -= factor * inv[k * n + col];
                }
            }
        }
    }
    Ok(inv)
}

// ---------------------------------------------------------------------------
// Stage 1: Assembly data structures
// ---------------------------------------------------------------------------

/// Boundary or interior interface condition on a local face of an element.
#[derive(Clone, Debug)]
pub enum DgtdFaceKind {
    /// Interior face shared with neighbor cell `neighbor_cell` (or periodic counterpart with `phase`).
    Interior { neighbor_cell: usize },
    /// Perfect Electric Conductor ($\mathbf{E}_t^* = 0$).
    Pec,
    /// Perfect Magnetic Conductor ($\mathbf{H}_t^* = 0$).
    Pmc,
    /// First-order Silver–Müller absorbing boundary condition.
    SilverMuller,
    /// Modal waveguide or lumped port `port_index`.
    Port { port_index: usize },
}

/// One canonical quadrature point on a local face of element $K^-$.
#[derive(Clone, Debug)]
pub struct DgtdFaceQp {
    /// Physical position $\mathbf{x}_q$.
    pub x: [f64; 3],
    /// Surface integration weight $w_q dS_q$.
    pub w_ds: f64,
    /// Outward unit normal $\hat{\mathbf{n}}_q^-$ of element $K^-$.
    pub n_hat: [f64; 3],
    /// Non-zero basis values `(local_node_b, phi_b(xi_q^-))` on element $K^-$.
    pub trace_minus: Vec<(usize, f64)>,
    /// Non-zero basis values `(neighbor_local_node_b, phi_b(xi_q^+))` on neighbor element $K^+$.
    pub trace_plus: Vec<(usize, f64)>,
    /// Pre-multiplied mass-inverted lift weights `(local_node_a, (M_K^-1)_{ab} w_q dS_q phi_b(xi_q^-))`.
    pub lift_minus: Vec<(usize, f64)>,
    /// Optional modal electric field profile $\hat{\mathbf{e}}_{t,p}(\mathbf{x}_q)$ and magnetic profile
    /// $\hat{\mathbf{h}}_{t,p}^{\text{in}}(\mathbf{x}_q)$ at this face quadrature point when the face belongs to a port.
    pub port_mode: Option<([f64; 3], [f64; 3])>,
    /// Optional TFSF sign (`+1.0` if $K^-$ is Total-Field and $K^+$ is Scattered-Field; `-1.0` if vice versa)
    /// and background wave parameters `(k_hat, pol, origin, c_bg, y_bg)`.
    pub tfsf: Option<(f64, [f64; 3], [f64; 3], [f64; 3], f64, f64)>,
}

/// Precomputed face coupling and boundary data for one local face of element $K^-$.
#[derive(Clone, Debug)]
pub struct DgtdCellFace {
    pub kind: DgtdFaceKind,
    /// Exterior wave admittance $Y^+$ ($\text{S}$).
    pub y_plus: f64,
    /// Optional penetrable zero-thickness 2D sheet (`FR-35`):
    /// `(sigma_2d_re, u_hat, v_hat)` where `sigma_2d_re` is the $2\times 2$ real tangential sheet conductivity in Siemens.
    pub sheet_2d: Option<([[f64; 2]; 2], [f64; 3], [f64; 3])>,
    /// Optional boundary facet label (`0` for unlabelled interior faces).
    pub facet_label: i32,
    pub qps: Vec<DgtdFaceQp>,
}

/// Precomputed per-element nonlinear optical constitutive parameters (`FR-41` `[nonlinear]`, `ADR 143`).
#[derive(Clone, Debug)]
pub struct DgtdCellNonlinear {
    /// Global-frame contracted $3 \times 6$ $\chi^{(2)}$ tensor $d_{iJ}$ ($\text{m/V}$) when `[nonlinear.chi2]` is active.
    pub chi2_d_global: Option<[[f64; 6]; 3]>,
    /// Spatial QPM poling sign $s_{\text{QPM}}(\mathbf{x}_a) \in \{-1, +1\}$ at each of the $N_p$ element nodes.
    pub node_qpm_signs: Vec<f64>,
    /// Isotropic Kerr susceptibility $\chi^{(3)}_{xxxx} = \frac{4}{3}(A + \frac{1}{2}B)$ ($\text{m}^2/\text{V}^2$) when `!has_cubic`.
    pub chi3_iso_xxxx: f64,
    /// Full rotated global rank-4 $\chi^{(3)}_{ijkl}$ ($\text{m}^2/\text{V}^2$) when cubic anisotropy is active.
    pub chi3_rank4_global: Option<Box<[[[[f64; 3]; 3]; 3]; 3]>>,
    /// Two-photon absorption coefficient $\beta_{\text{TPA}}$ ($\text{m/W}$).
    pub beta_tpa: f64,
    /// Free-carrier absorption cross-section $\sigma_{\text{FCA}}$ ($\text{m}^2$).
    pub sigma_fca: f64,
    /// Free-carrier dispersion coefficient $k_{\text{FCD}}$ ($\text{m}^3$).
    pub k_fcd: f64,
    /// Free-carrier lifetime $\tau_c$ ($\text{s}$).
    pub tau_c: f64,
    /// Reference linear refractive index $n_0 = \sqrt{\varepsilon_{\text{iso}}}$.
    pub n0: f64,
    /// Reference carrier angular frequency $\omega_0$ ($\text{rad/s}$).
    pub omega_0: f64,
    /// Minimum diagonal linear permittivity $\varepsilon_{\min}$ for physical susceptibility clamping.
    pub eps_min_diag: f64,
}

/// Precomputed local operators and constitutive data for one DG element $K$.
#[derive(Clone, Debug)]
pub struct DgtdCell {
    pub cell_index: usize,
    pub cell_type: CellType,
    pub domain_label: i32,
    pub n_p: usize,
    /// Flat global DOF offset `dof_offset..dof_offset + n_p` in the global nodal arrays.
    pub dof_offset: usize,
    /// Physical coordinates of the $N_p$ local nodes.
    pub node_coords: Vec<[f64; 3]>,
    /// Quadrature weights / mass matrix $M_K$ (`len = n_p` if `is_diagonal_mass`, else `n_p * n_p`) for energy integration.
    pub is_diagonal_mass: bool,
    pub mass: Vec<f64>,
    pub inv_mass: Vec<f64>,
    /// Strong-form SBP local differentiation matrices $D^x_K, D^y_K, D^z_K$ (`n_p * n_p` row-major each).
    pub dx: Vec<f64>,
    pub dy: Vec<f64>,
    pub dz: Vec<f64>,
    /// Instantaneous relative permittivity tensor $\boldsymbol{\varepsilon}_\infty$ and its inverse $(\varepsilon_0 \boldsymbol{\varepsilon}_\infty)^{-1}$.
    pub eps_r: [[f64; 3]; 3],
    pub inv_eps: [[f64; 3]; 3],
    /// Instantaneous relative permeability tensor $\boldsymbol{\mu}_\infty$ and its inverse $(\mu_0 \boldsymbol{\mu}_\infty)^{-1}$.
    pub mu_r: [[f64; 3]; 3],
    pub inv_mu: [[f64; 3]; 3],
    /// Electric conductivity / gyrotropic coupling tensor $\boldsymbol{\sigma}$ ($\text{S/m}$).
    pub sigma_e: [[f64; 3]; 3],
    pub has_sigma_e: bool,
    /// Characteristic scalar wave admittance $Y_K = \sqrt{\varepsilon_0 \varepsilon_r / (\mu_0 \mu_r)}$.
    pub y_admittance: f64,
    /// Uniaxial PML damping rates $[\sigma_x, \sigma_y, \sigma_z]$ in $\text{s}^{-1}$ (`[0.0; 3]` outside PML).
    pub pml_sigma: [f64; 3],
    pub in_pml: bool,
    /// Dispersive ADE poles active in this cell.
    pub poles: Vec<DgtdDispersionPole>,
    /// Offset into the global auxiliary pole state array (2 vector states $[\mathbf{J}_m, \mathbf{P}_m]$ per pole per node).
    pub pole_offset: usize,
    /// Optional instantaneous $\chi^{(2)}$ / QPM / $\chi^{(3)}$ Kerr + TPA/FCA/FCD nonlinear material state (`ADR 143`).
    pub nonlinear: Option<Box<DgtdCellNonlinear>>,
    /// Element volume (3D) or area (2D), inscribed length scale $h_K = V_K / A_K$, and local stable $\Delta t_K$.
    pub measure: f64,
    pub h_scale: f64,
    pub dt_local: f64,
    /// Assigned power-of-two LTS tier $\ell \in 0..n_{\text{tiers}}$.
    pub lts_tier: usize,
    /// Local faces of this cell.
    pub faces: Vec<DgtdCellFace>,
}

/// Assembled modal or lumped port for broadband `DrivenPortData` extraction and time-domain injection.
#[derive(Clone, Debug)]
pub struct DgtdPort {
    pub name: String,
    pub labels: Vec<i32>,
    pub n_out: [f64; 3],
    pub u_hat: [f64; 3],
    pub v_hat: [f64; 3],
    /// Quadratic dispersion coefficients $[c_0, c_1, c_2]$ for $\beta_p^2(k_0^2) = c_0 + c_1 k_0^2 + c_2 k_0^4$.
    pub beta2_poly: [f64; 3],
    /// Optional explicit circuit impedance override $Z_0$ ($\Omega$).
    pub z0_override: Option<f64>,
    /// Modal wave admittance $Y_0(\omega_0) = \beta_p(\omega_0) / (\omega_0 \mu_0)$ at the carrier frequency $\omega_0$.
    pub y0_carrier: f64,
    /// Whether this port is a lumped port (`true`) or waveguide modal port (`false`).
    pub is_lumped: bool,
}

impl DgtdPort {
    /// Returns `(beta(omega), z_wave(omega), z_circuit(omega))` at angular frequency `omega` (rad/s).
    #[inline]
    pub fn impedance_at(&self, omega: f64) -> (f64, f64, f64) {
        if self.is_lumped {
            let z0 = self.z0_override.unwrap_or(50.0);
            let k0 = (omega / C0).max(1e-12);
            return (k0, z0, z0);
        }
        let k0 = omega / C0;
        let s = k0 * k0;
        let [c0, c1, c2] = self.beta2_poly;
        let beta = (c0 + c1 * s + c2 * s * s).max(1e-12).sqrt();
        let z_wave = (omega * MU0 / beta).max(1e-12);
        let z_circ = self.z0_override.unwrap_or(z_wave);
        (beta, z_wave, z_circ)
    }
}

/// Point dipole source pre-mapped to its owning cell and mass-inverted basis weights.
#[derive(Clone, Debug)]
pub struct DgtdDipoleSource {
    pub cell_index: usize,
    pub moment: [f64; 3],
    pub magnetic: bool,
    /// Pre-multiplied `(local_node_a, (M_K^-1)_{ab} phi_b(xi_dipole))`.
    pub lift_weights: Vec<(usize, f64)>,
}

/// Precomputed probe interpolation weights for a `FieldTimeMonitor`.
#[derive(Clone, Debug)]
pub struct DgtdProbe {
    pub position: [f64; 3],
    pub cell_index: usize,
    pub basis_weights: Vec<(usize, f64)>,
}

/// Precomputed surface quadrature references for a `FluxMonitor`.
#[derive(Clone, Debug)]
pub struct DgtdFluxMonitor {
    pub name: String,
    /// `(cell_index, local_face_index)` pairs belonging to this flux monitor.
    pub cell_faces: Vec<(usize, usize)>,
}

/// Stage 1 output of the DGTD solver (`ADR 9`, `ADR 144`).
#[derive(Clone, Debug)]
pub struct DgtdAssembly {
    pub dim: usize,
    pub basis_order: usize,
    pub config: DgtdConfig,
    pub mixed: MixedDiscretisation,
    pub cells: Vec<DgtdCell>,
    pub total_dofs: usize,
    pub total_pole_states: usize,
    pub dt_0: f64,
    pub dt_max: f64,
    pub h_min: f64,
    pub h_max: f64,
    pub n_lts_tiers: usize,
    pub tier_cells: Vec<Vec<usize>>,
    /// For each tier $\ell$, the 1-ring of cells in tier $\ell - 1$ that share a face with a cell in tier $\ell$.
    pub tier_finer_border: Vec<Vec<usize>>,
    pub ports: Vec<DgtdPort>,
    pub dipoles: Vec<DgtdDipoleSource>,
    pub probes: Vec<DgtdProbe>,
    pub flux_monitors: Vec<DgtdFluxMonitor>,
    pub symmetry_factor: f64,
    /// True when at least one non-PML cell has TPA-generated free-carrier FCA/FCD dynamics (`drude_carrier_step`).
    pub has_carrier_dynamics: bool,
}

// ---------------------------------------------------------------------------
// Canonical face keys, curved face evaluation, and port/sheet helpers
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum FaceKey {
    Edge2([u32; 2]),
    Tri3([u32; 3]),
    Quad4([u32; 4]),
}

fn cell_local_face_vertices(mesh: &MixedDiscretisation, ci: usize) -> Vec<(FacetType, Vec<u32>)> {
    let ct = mesh.cell_type(ci);
    let vs = mesh.cell(ci);
    let mut out = Vec::new();
    if mesh.dim == 2 {
        for &(i, j) in cell_local_edges(ct) {
            out.push((FacetType::Edge2, vec![vs[i], vs[j]]));
        }
    } else {
        for &[i, j, k] in cell_local_tri_faces(ct) {
            out.push((FacetType::Tri3, vec![vs[i], vs[j], vs[k]]));
        }
        for &[i, j, k, l] in cell_local_quad_faces(ct) {
            out.push((FacetType::Quad4, vec![vs[i], vs[j], vs[k], vs[l]]));
        }
    }
    out
}

fn canonical_face_key(ft: FacetType, fv: &[u32]) -> FaceKey {
    match ft {
        FacetType::Edge2 => FaceKey::Edge2([fv[0].min(fv[1]), fv[0].max(fv[1])]),
        FacetType::Tri3 => FaceKey::Tri3(sorted_tri_face([fv[0], fv[1], fv[2]])),
        FacetType::Quad4 => FaceKey::Quad4(canonical_quad_face([fv[0], fv[1], fv[2], fv[3]])),
    }
}

fn canonical_face_verts(key: &FaceKey) -> &[u32] {
    match key {
        FaceKey::Edge2(v) => v,
        FaceKey::Tri3(v) => v,
        FaceKey::Quad4(v) => v,
    }
}

/// Evaluates canonical face quadrature rules:
/// - `Edge2`: 1D GLL rule `Gll::new(p)` on `[-1, 1]`
/// - `Quad4`: 2D GLL rule `Gll::new(p) x Gll::new(p)` on `[-1, 1]^2`
/// - `Tri3`: symmetric triangle quadrature rule of degree `2 * p + 2`
fn canonical_face_quadrature(ft: FacetType, p: usize, gll: &Gll) -> (Vec<[f64; 2]>, Vec<f64>) {
    match ft {
        FacetType::Edge2 => {
            let pts = gll.nodes.iter().map(|&t| [t, 0.0]).collect();
            let w = gll.weights.clone();
            (pts, w)
        }
        FacetType::Quad4 => {
            let n = gll.nodes.len();
            let mut pts = Vec::with_capacity(n * n);
            let mut w = Vec::with_capacity(n * n);
            for (b, &eta) in gll.nodes.iter().enumerate() {
                for (a, &xi) in gll.nodes.iter().enumerate() {
                    pts.push([xi, eta]);
                    w.push(gll.weights[a] * gll.weights[b]);
                }
            }
            (pts, w)
        }
        FacetType::Tri3 => {
            let rule = facet_rule_mixed(FacetType::Tri3, 2 * p + 2);
            (rule.pts, rule.w)
        }
    }
}

/// Evaluates physical position `x`, surface Jacobian measure `ds`, and canonical unit normal `n_canon`
/// on a canonical face `key` at reference coordinate `xi_f`.
fn eval_canonical_face_geometry(
    mesh: &MixedDiscretisation,
    key: &FaceKey,
    xi_f: &[f64; 2],
) -> ([f64; 3], f64, [f64; 3]) {
    let cv = canonical_face_verts(key);
    let corners: Vec<[f64; 3]> = cv.iter().map(|&v| mesh.point3(v as usize)).collect();
    match *key {
        FaceKey::Edge2([v0, v1]) => {
            let eid = mesh.edge_index(v0, v1);
            let (space, nodes) = if let Some(mid) = eid.and_then(|e| mesh.edge_node(e)) {
                (LagrangeSpace::Order2Serendipity, vec![corners[0], corners[1], mid])
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(FacetType::Edge2, space, xi_f);
            let mut x = [0.0; 3];
            let mut t = [0.0; 2];
            for (a, p) in nodes.iter().enumerate() {
                x[0] += n[a] * p[0];
                x[1] += n[a] * p[1];
                t[0] += dn[a][0] * p[0];
                t[1] += dn[a][0] * p[1];
            }
            let ds = t[0].hypot(t[1]);
            let n_hat = if ds > 0.0 { [t[1] / ds, -t[0] / ds, 0.0] } else { [1.0, 0.0, 0.0] };
            (x, ds, n_hat)
        }
        FaceKey::Tri3([v0, v1, v2]) => {
            let e01 = mesh.edge_index(v0, v1);
            let e02 = mesh.edge_index(v0, v2);
            let e12 = mesh.edge_index(v1, v2);
            let any_curved = e01.and_then(|e| mesh.edge_node(e)).is_some()
                || e02.and_then(|e| mesh.edge_node(e)).is_some()
                || e12.and_then(|e| mesh.edge_node(e)).is_some();
            let (space, nodes) = if any_curved {
                let mid = |eid: Option<usize>, i: usize, j: usize| {
                    eid.and_then(|e| mesh.edge_node(e)).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ])
                };
                (
                    LagrangeSpace::Order2Serendipity,
                    vec![corners[0], corners[1], corners[2], mid(e01, 0, 1), mid(e02, 0, 2), mid(e12, 1, 2)],
                )
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(FacetType::Tri3, space, xi_f);
            let (mut x, mut t0, mut t1) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            for (a, p) in nodes.iter().enumerate() {
                for k in 0..3 {
                    x[k] += n[a] * p[k];
                    t0[k] += dn[a][0] * p[k];
                    t1[k] += dn[a][1] * p[k];
                }
            }
            let cx = cross3(t0, t1);
            let ds = norm3(cx);
            let n_hat = if ds > 0.0 { [cx[0] / ds, cx[1] / ds, cx[2] / ds] } else { [0.0, 0.0, 1.0] };
            (x, ds, n_hat)
        }
        FaceKey::Quad4(qv) => {
            let qf_id = mesh.quad_face_index(qv);
            let eids = [
                mesh.edge_index(qv[0], qv[1]),
                mesh.edge_index(qv[1], qv[2]),
                mesh.edge_index(qv[2], qv[3]),
                mesh.edge_index(qv[3], qv[0]),
            ];
            let has_qfn = qf_id.and_then(|q| mesh.quad_face_node(q)).is_some();
            let any_curved = has_qfn || eids.iter().any(|&e| e.and_then(|id| mesh.edge_node(id)).is_some());
            let (space, nodes) = if any_curved {
                let mut ns = corners.clone();
                let mut mids = [[0.0; 3]; 4];
                for (k, &(i, j)) in QUAD_EDGES.iter().enumerate() {
                    let m = eids[k].and_then(|id| mesh.edge_node(id)).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ]);
                    mids[k] = m;
                    ns.push(m);
                }
                if has_qfn {
                    ns.push(qf_id.and_then(|q| mesh.quad_face_node(q)).unwrap_or_else(|| {
                        transfinite_quad_center([corners[0], corners[1], corners[2], corners[3]], mids)
                    }));
                    (LagrangeSpace::Order2Tensor, ns)
                } else {
                    (LagrangeSpace::Order2Serendipity, ns)
                }
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(FacetType::Quad4, space, xi_f);
            let (mut x, mut t0, mut t1) = ([0.0; 3], [0.0; 3], [0.0; 3]);
            for (a, p) in nodes.iter().enumerate() {
                for k in 0..3 {
                    x[k] += n[a] * p[k];
                    t0[k] += dn[a][0] * p[k];
                    t1[k] += dn[a][1] * p[k];
                }
            }
            let cx = cross3(t0, t1);
            let ds = norm3(cx);
            let n_hat = if ds > 0.0 { [cx[0] / ds, cx[1] / ds, cx[2] / ds] } else { [0.0, 0.0, 1.0] };
            (x, ds, n_hat)
        }
    }
}

/// Maps canonical face reference coordinates `xi_f` into the reference coordinates `xi_cell` of cell `ci`,
/// along with the reference vector `xi_face_center - xi_opposite_center` pointing outward in reference space.
fn map_face_to_cell_ref(mesh: &MixedDiscretisation, ci: usize, key: &FaceKey, xi_f: &[f64; 2]) -> ([f64; 3], [f64; 3]) {
    let ct = mesh.cell_type(ci);
    let cell_vs = mesh.cell(ci);
    let corners_ref = cell_corner_ref_coords(ct);
    let f_vs = canonical_face_verts(key);
    let mut loc_idx = Vec::with_capacity(f_vs.len());
    for &fv in f_vs {
        let li = cell_vs.iter().position(|&cv| cv == fv).unwrap_or(0);
        loc_idx.push(li);
    }
    let (fn1, _) = match key {
        FaceKey::Edge2(_) => facet_lagrange_shape(FacetType::Edge2, LagrangeSpace::Order1, xi_f),
        FaceKey::Tri3(_) => facet_lagrange_shape(FacetType::Tri3, LagrangeSpace::Order1, xi_f),
        FaceKey::Quad4(_) => facet_lagrange_shape(FacetType::Quad4, LagrangeSpace::Order1, xi_f),
    };
    let mut xi_cell = [0.0; 3];
    let mut xi_f_cent = [0.0; 3];
    for (m, &li) in loc_idx.iter().enumerate() {
        let cr = corners_ref[li];
        for a in 0..3 {
            xi_cell[a] += fn1[m] * cr[a];
            xi_f_cent[a] += cr[a] / (loc_idx.len() as f64);
        }
    }
    let mut xi_opp = [0.0; 3];
    let mut n_opp = 0usize;
    for (li, cr) in corners_ref.iter().enumerate() {
        if !loc_idx.contains(&li) {
            for a in 0..3 {
                xi_opp[a] += cr[a];
            }
            n_opp += 1;
        }
    }
    if n_opp > 0 {
        for a in 0..3 {
            xi_opp[a] /= n_opp as f64;
        }
    }
    let v_out_ref = sub3(xi_f_cent, xi_opp);
    (xi_cell, v_out_ref)
}

/// Builds a canonical in-plane orthonormal frame `(u_hat, v_hat)` for a port plane with normal `n_out`.
fn canonical_port_axes(n_out: [f64; 3], dim: usize) -> ([f64; 3], [f64; 3]) {
    if dim == 2 {
        // In 2D (x, y, 0), the in-plane tangent is u_hat = (-n_y, n_x, 0) and out-of-plane is v_hat = (0, 0, 1).
        let u_raw = [-n_out[1], n_out[0], 0.0];
        let ul = norm3(u_raw).max(1e-30);
        let u_hat = [u_raw[0] / ul, u_raw[1] / ul, 0.0];
        let mut sign = 1.0;
        for k in 0..2 {
            if u_hat[k].abs() > 1e-6 {
                if u_hat[k] < 0.0 {
                    sign = -1.0;
                }
                break;
            }
        }
        return ([sign * u_hat[0], sign * u_hat[1], 0.0], [0.0, 0.0, 1.0]);
    }
    let mut sign = 1.0;
    for k in 0..3 {
        if n_out[k].abs() > 1e-6 {
            if n_out[k] < 0.0 {
                sign = -1.0;
            }
            break;
        }
    }
    let n_canon = [sign * n_out[0], sign * n_out[1], sign * n_out[2]];
    let ref_axis = if n_canon[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let rn = dot3(ref_axis, n_canon);
    let u_raw = [ref_axis[0] - rn * n_canon[0], ref_axis[1] - rn * n_canon[1], ref_axis[2] - rn * n_canon[2]];
    let ul = norm3(u_raw).max(1e-30);
    let u_hat = [u_raw[0] / ul, u_raw[1] / ul, u_raw[2] / ul];
    let v_raw = cross3(n_canon, u_hat);
    let vl = norm3(v_raw).max(1e-30);
    let v_hat = [v_raw[0] / vl, v_raw[1] / vl, v_raw[2] / vl];
    (u_hat, v_hat)
}

/// Solves a 3D `waveport` modal cross-section using `nyx_ems::modes` (or TEM fallback if all boundaries are PMC),
/// returning `(DgtdPort, modal_evaluator_triangles)` so every face quadrature point on the port gets its exact
/// modal $(\hat{\mathbf{e}}_{t,p}(\mathbf{x}_q), \hat{\mathbf{h}}_{t,p}^{\text{in}}(\mathbf{x}_q))$ profile!
struct SolvedPortProfile {
    port: DgtdPort,
    /// For 3D ports: `(tri_2d_coords, cell_2d_index)`, 2D mesh, edge coefficients `et`, `order`, and scale `inv_norm`.
    mode_2d: Option<(SimplexMesh, usize, Vec<f64>, f64)>,
    /// For 2D or uniform TEM ports: constant normalized electric polarization vector `e_const`.
    e_const: [f64; 3],
}

impl SolvedPortProfile {
    fn eval_at(&self, x: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let e_t = if let Some((ref mesh_2d, order, ref et, scale)) = self.mode_2d {
            let u = dot3(x, self.port.u_hat);
            let v = dot3(x, self.port.v_hat);
            let geo_2d = geometry(mesh_2d);
            let cd_2d = cell_data(mesh_2d, order);
            let nloc_e = if order == 1 { 3 } else { 8 };
            let mut w2 = vec![[0.0; 2]; nloc_e];
            let mut curl2 = vec![0.0; nloc_e];
            let mut best_c = 0usize;
            let mut best_l = [1.0 / 3.0; 3];
            let mut best_min_l = f64::NEG_INFINITY;
            for c_idx in 0..mesh_2d.n_cells() {
                let cv = mesh_2d.cell(c_idx);
                let p0 = mesh_2d.point(cv[0] as usize);
                let p1 = mesh_2d.point(cv[1] as usize);
                let p2 = mesh_2d.point(cv[2] as usize);
                let det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
                if det.abs() <= 1e-30 {
                    continue;
                }
                let l1 = ((u - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (v - p0[1])) / det;
                let l2 = ((p1[0] - p0[0]) * (v - p0[1]) - (u - p0[0]) * (p1[1] - p0[1])) / det;
                let l0 = 1.0 - l1 - l2;
                let min_l = l0.min(l1).min(l2);
                if min_l > best_min_l {
                    best_min_l = min_l;
                    best_c = c_idx;
                    best_l = [l0.clamp(0.0, 1.0), l1.clamp(0.0, 1.0), l2.clamp(0.0, 1.0)];
                }
                if min_l >= -1e-7 {
                    break;
                }
            }
            let g2 = tri_g(&geo_2d, best_c);
            edge_basis(&best_l, &g2, &cd_2d.sign[best_c], order, &mut w2, &mut curl2);
            let (mut eu, mut ev) = (0.0, 0.0);
            for s in 0..nloc_e {
                let coef = et[cd_2d.gdofs[best_c][s]];
                eu += coef * w2[s][0];
                ev += coef * w2[s][1];
            }
            eu *= scale;
            ev *= scale;
            [
                eu * self.port.u_hat[0] + ev * self.port.v_hat[0],
                eu * self.port.u_hat[1] + ev * self.port.v_hat[1],
                eu * self.port.u_hat[2] + ev * self.port.v_hat[2],
            ]
        } else {
            self.e_const
        };
        // Inward propagation direction is n_in = -n_out, so h_in = Y_0(omega_0) * (n_in x e_t) = -Y_0 * (n_out x e_t)
        let n_cross_e = cross3(self.port.n_out, e_t);
        let y0 = self.port.y0_carrier;
        let h_in = [-y0 * n_cross_e[0], -y0 * n_cross_e[1], -y0 * n_cross_e[2]];
        (e_t, h_in)
    }
}

#[allow(clippy::too_many_arguments)]
fn build_dgtd_port_profile(
    mesh: &MixedDiscretisation,
    problem: &Problem,
    config: &DgtdConfig,
    name: &str,
    labels: &[i32],
    mode_idx: usize,
    n_eff_guess: Option<f64>,
    z0_override: Option<f64>,
    is_lumped: bool,
    lumped_dir: Option<[f64; 3]>,
    cell_eps: &[[[f64; 3]; 3]],
    cell_mur: &[f64],
    face_owner_map: &HashMap<FaceKey, Vec<(usize, usize)>>,
) -> Result<SolvedPortProfile, String> {
    let f0 = config.pulse.f0_hz().max(1.0);
    let f_min = config.freqs_hz.first().copied().unwrap_or(f0).min(f0);
    let f_max = config.freqs_hz.last().copied().unwrap_or(f0).max(f0);
    let f_mid = 0.5 * (f_min + f_max);

    // Collect all boundary facets belonging to `labels`
    let mut port_facets = Vec::new();
    let mut seen_keys = BTreeSet::new();
    let mut n_acc = [0.0; 3];
    let mut total_area = 0.0;
    let mut avg_eps = 0.0;
    let mut avg_mur = 0.0;

    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) || problem.pec.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let ft = mesh.facet_type(f);
        let fv = mesh.facet(f);
        let key = canonical_face_key(ft, fv);
        if !seen_keys.insert(key) {
            continue;
        }
        let Some(owners) = face_owner_map.get(&key) else { continue };
        let (ci, _) = owners[0];
        let (_x_c, ds_c, mut n_c) = eval_canonical_face_geometry(mesh, &key, &[0.0, 0.0]);
        let (xi_c, v_out_ref) = map_face_to_cell_ref(mesh, ci, &key, &[0.0, 0.0]);
        let iso = IsoparametricCell::from_mesh(mesh, ci);
        let jac = iso.eval(&xi_c);
        let v_out_phys = if mesh.dim == 2 {
            [
                jac.j[0][0] * v_out_ref[0] + jac.j[0][1] * v_out_ref[1],
                jac.j[1][0] * v_out_ref[0] + jac.j[1][1] * v_out_ref[1],
                0.0,
            ]
        } else {
            mat3_mul_vec(&jac.j, v_out_ref)
        };
        if dot3(n_c, v_out_phys) < 0.0 {
            n_c = [-n_c[0], -n_c[1], -n_c[2]];
        }
        let ref_meas = match ft {
            FacetType::Edge2 => 2.0,
            FacetType::Tri3 => 0.5,
            FacetType::Quad4 => 4.0,
        };
        let area_f = ds_c * ref_meas;
        total_area += area_f;
        for k in 0..3 {
            n_acc[k] += area_f * n_c[k];
        }
        let e_iso = (cell_eps[ci][0][0] + cell_eps[ci][1][1] + cell_eps[ci][2][2]) / 3.0;
        avg_eps += area_f * e_iso.max(1.0);
        avg_mur += area_f * cell_mur[ci].max(1.0);
        port_facets.push((ft, fv.to_vec(), ci, area_f));
    }
    if port_facets.is_empty() || total_area <= 0.0 {
        return Err(format!("dgtd port '{name}' selects no active boundary facets"));
    }
    avg_eps /= total_area;
    avg_mur /= total_area;
    let nl = norm3(n_acc).max(1e-30);
    let n_out = [n_acc[0] / nl, n_acc[1] / nl, n_acc[2] / nl];
    let (u_hat, v_hat) = canonical_port_axes(n_out, mesh.dim);

    // 2D mesh or Lumped port or 3D TEM fallback helper
    let make_tem_fallback = |pol_dir: [f64; 3]| -> SolvedPortProfile {
        let mut e_dir = [
            pol_dir[0] - dot3(pol_dir, n_out) * n_out[0],
            pol_dir[1] - dot3(pol_dir, n_out) * n_out[1],
            pol_dir[2] - dot3(pol_dir, n_out) * n_out[2],
        ];
        let el = norm3(e_dir).max(1e-30);
        for k in 0..3 {
            e_dir[k] /= el;
        }
        // Normalize so int_{Gamma_p} mu_r^-1 |e_t|^2 dS = 1
        let norm_scale = (avg_mur / total_area).sqrt();
        let e_const = [e_dir[0] * norm_scale, e_dir[1] * norm_scale, e_dir[2] * norm_scale];
        let beta2_poly = [0.0, avg_eps * avg_mur, 0.0];
        let omega_0 = 2.0 * PI * f0;
        let beta_0 = omega_0 * (avg_eps * avg_mur).sqrt() / C0;
        let y0_carrier = if is_lumped { 1.0 / z0_override.unwrap_or(50.0) } else { beta_0 / (omega_0 * MU0) };
        SolvedPortProfile {
            port: DgtdPort {
                name: name.to_string(),
                labels: labels.to_vec(),
                n_out,
                u_hat,
                v_hat,
                beta2_poly,
                z0_override,
                y0_carrier,
                is_lumped,
            },
            mode_2d: None,
            e_const,
        }
    };

    if mesh.dim == 2 || is_lumped {
        let default_pol = if let Some(d) = lumped_dir {
            d
        } else if mesh.dim == 2 {
            match config.polarization {
                DgtdPolarization::Te => u_hat,
                DgtdPolarization::Tm | DgtdPolarization::Full => [0.0, 0.0, 1.0],
            }
        } else {
            v_hat
        };
        return Ok(make_tem_fallback(default_pol));
    }

    // 3D waveport: build 2D SimplexMesh by triangulating Tri3 and Quad4 facets on the port plane
    let mut verts_set = BTreeSet::new();
    let mut tris_3d: Vec<([u32; 3], usize)> = Vec::new();
    for (ft, fv, ci, _) in &port_facets {
        for &v in fv {
            verts_set.insert(v);
        }
        match ft {
            FacetType::Tri3 => tris_3d.push((sorted3(fv[0], fv[1], fv[2]), *ci)),
            FacetType::Quad4 => {
                tris_3d.push((sorted3(fv[0], fv[1], fv[2]), *ci));
                tris_3d.push((sorted3(fv[0], fv[2], fv[3]), *ci));
            }
            FacetType::Edge2 => {}
        }
    }
    let verts_3d: Vec<u32> = verts_set.into_iter().collect();
    let mut points_2d = Vec::with_capacity(verts_3d.len() * 2);
    for &v3 in &verts_3d {
        let p = mesh.point3(v3 as usize);
        points_2d.push(dot3(p, u_hat));
        points_2d.push(dot3(p, v_hat));
    }
    let v2d_of = |v3: u32| -> usize { verts_3d.binary_search(&v3).unwrap_or(0) };

    let (sym_pec, sym_pmc) = symmetry_facet_flags_mixed(mesh, problem.symmetry)?;
    let all_wp_labels: Vec<i32> = problem.waveports.iter().flat_map(|w| w.labels.iter().copied()).collect();
    let all_ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let mut pec_edges_3d = BTreeSet::new();
    let mut pmc_edges_3d = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        let is_pec = problem.pec.contains(&l)
            || sym_pec[f]
            || (mesh.exterior_facets[f]
                && !problem.pmc.contains(&l)
                && !sym_pmc[f]
                && !all_ibc_labels.contains(&l)
                && !all_wp_labels.contains(&l)
                && !config.absorbing_labels.contains(&l));
        let is_pmc = problem.pmc.contains(&l) || sym_pmc[f];
        for i in 0..fv.len() {
            let a = fv[i];
            let b = fv[(i + 1) % fv.len()];
            let e = (a.min(b), a.max(b));
            if is_pec {
                pec_edges_3d.insert(e);
            } else if is_pmc {
                pmc_edges_3d.insert(e);
            }
        }
    }

    let nc_2d = tris_3d.len();
    let mut cells_2d = Vec::with_capacity(nc_2d * 3);
    let mut cell_labels_2d = Vec::with_capacity(nc_2d);
    let mut eps_2d = Vec::with_capacity(nc_2d);
    let mut mur_2d = Vec::with_capacity(nc_2d);
    let mut edge_use_count: BTreeMap<(u32, u32), usize> = BTreeMap::new();

    for (c_idx, &(tri, ci)) in tris_3d.iter().enumerate() {
        let (i0, i1, i2) = (v2d_of(tri[0]) as u32, v2d_of(tri[1]) as u32, v2d_of(tri[2]) as u32);
        cells_2d.extend_from_slice(&[i0, i1, i2]);
        cell_labels_2d.push(c_idx as i32 + 1);
        for (a, b) in [(i0, i1), (i1, i2), (i0, i2)] {
            *edge_use_count.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
        let t3 = &cell_eps[ci];
        let proj = |a: [f64; 3], b: [f64; 3]| -> f64 {
            let mut s = 0.0;
            for r in 0..3 {
                for cc in 0..3 {
                    s += a[r] * t3[r][cc] * b[cc];
                }
            }
            s
        };
        eps_2d.push(EpsCell {
            xx: proj(u_hat, u_hat),
            xy: proj(u_hat, v_hat),
            yy: proj(v_hat, v_hat),
            zz: proj(n_out, n_out),
        });
        mur_2d.push(cell_mur[ci]);
    }

    let mut facets_2d = Vec::new();
    let mut facet_labels_2d = Vec::new();
    let mut has_pec_edge = false;
    for (&(ia, ib), &cnt) in &edge_use_count {
        let (va, vb) = (verts_3d[ia as usize], verts_3d[ib as usize]);
        let e3 = (va.min(vb), va.max(vb));
        let is_perim = cnt == 1;
        if pec_edges_3d.contains(&e3) || (is_perim && !pmc_edges_3d.contains(&e3)) {
            facets_2d.extend_from_slice(&[ia, ib]);
            facet_labels_2d.push(1);
            has_pec_edge = true;
        } else if is_perim && pmc_edges_3d.contains(&e3) {
            facets_2d.extend_from_slice(&[ia, ib]);
            facet_labels_2d.push(2);
        }
    }

    // If the cross-section is homogeneous and has no PEC walls (e.g. a pure parallel-plate / periodic TEM slab), use TEM mode
    let eps_spread =
        eps_2d.iter().map(|e| (e.xx - eps_2d[0].xx).abs().max((e.yy - eps_2d[0].yy).abs())).fold(0.0_f64, f64::max);
    if !has_pec_edge && eps_spread < 1e-6 {
        return Ok(make_tem_fallback(v_hat));
    }

    let mesh_2d = SimplexMesh::new(2, points_2d, cells_2d, cell_labels_2d, facets_2d, facet_labels_2d);
    let mode_order = problem.order.clamp(1, 2);
    let n_max =
        eps_2d.iter().zip(&mur_2d).map(|(e, &m)| ((e.xx.max(e.yy)) * m).max(1.0).sqrt()).fold(1.0_f64, f64::max);
    let (n_search, k_ask) =
        if let Some(ng) = n_eff_guess { (ng, (mode_idx + 6).max(8)) } else { (1.02 * n_max, (mode_idx + 3).max(4)) };

    let geo_2d = geometry(&mesh_2d);
    let quad2 = femcore::fem::quad_rule(2, 4);
    let tri_weights: Vec<Vec<f64>> =
        (0..nc_2d).map(|c_idx| quad2.w.iter().map(|&wq| wq * geo_2d.vol[c_idx]).collect()).collect();

    let eval_norm_mode = |et: &[f64]| -> (Vec<Vec<[f64; 3]>>, f64, f64) {
        let cd_2d = cell_data(&mesh_2d, mode_order);
        let nloc_e = if mode_order == 1 { 3 } else { 8 };
        let mut w2 = vec![[0.0; 2]; nloc_e];
        let mut curl2 = vec![0.0; nloc_e];
        let mut fields = Vec::with_capacity(nc_2d);
        let mut norm2_acc = 0.0;
        let mut u_c = 0.0;
        let mut v_c = 0.0;
        let mut area_tot = 0.0;
        for c_idx in 0..nc_2d {
            let cv = mesh_2d.cell(c_idx);
            let (p0, p1, p2) =
                (mesh_2d.point(cv[0] as usize), mesh_2d.point(cv[1] as usize), mesh_2d.point(cv[2] as usize));
            let g2 = tri_g(&geo_2d, c_idx);
            let imu = 1.0 / mur_2d[c_idx];
            let mut cell_q = Vec::with_capacity(quad2.w.len());
            for q in 0..quad2.w.len() {
                let l = &quad2.lam[q];
                let uq = l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0];
                let vq = l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1];
                let wq = tri_weights[c_idx][q];
                u_c += wq * uq;
                v_c += wq * vq;
                area_tot += wq;
                edge_basis(l, &g2, &cd_2d.sign[c_idx], mode_order, &mut w2, &mut curl2);
                let (mut eu, mut ev) = (0.0, 0.0);
                for s in 0..nloc_e {
                    let coef = et[cd_2d.gdofs[c_idx][s]];
                    eu += coef * w2[s][0];
                    ev += coef * w2[s][1];
                }
                let e3 = [eu * u_hat[0] + ev * v_hat[0], eu * u_hat[1] + ev * v_hat[1], eu * u_hat[2] + ev * v_hat[2]];
                norm2_acc += wq * imu * dot3(e3, e3);
                cell_q.push(e3);
            }
            fields.push(cell_q);
        }
        let mut inv_norm = 1.0 / norm2_acc.max(1e-300).sqrt();
        for cq in &mut fields {
            for e3 in cq {
                e3[0] *= inv_norm;
                e3[1] *= inv_norm;
                e3[2] *= inv_norm;
            }
        }
        u_c /= area_tot.max(1e-300);
        v_c /= area_tot.max(1e-300);
        let span = area_tot.sqrt().max(1e-12);
        let pref_dir = [
            u_hat[0] + std::f64::consts::FRAC_1_PI * v_hat[0],
            u_hat[1] + std::f64::consts::FRAC_1_PI * v_hat[1],
            u_hat[2] + std::f64::consts::FRAC_1_PI * v_hat[2],
        ];
        let mut m0 = 0.0;
        let mut m1 = 0.0;
        let mut peak_val = 0.0;
        let mut peak_mag2 = -1.0;
        for c_idx in 0..nc_2d {
            let cv = mesh_2d.cell(c_idx);
            let (p0, p1, p2) =
                (mesh_2d.point(cv[0] as usize), mesh_2d.point(cv[1] as usize), mesh_2d.point(cv[2] as usize));
            for q in 0..quad2.w.len() {
                let l = &quad2.lam[q];
                let uq = (l[0] * p0[0] + l[1] * p1[0] + l[2] * p2[0] - u_c) / span;
                let vq = (l[0] * p0[1] + l[1] * p1[1] + l[2] * p2[1] - v_c) / span;
                let wq = tri_weights[c_idx][q];
                let e3 = fields[c_idx][q];
                let proj = dot3(e3, pref_dir);
                m0 += wq * proj;
                m1 += wq * (uq + 0.2718281828 * vq) * proj;
                let mag2 = dot3(e3, e3);
                if mag2 > peak_mag2 + 1e-12 && proj.abs() > 1e-6 * mag2.sqrt() {
                    peak_mag2 = mag2;
                    peak_val = proj;
                }
            }
        }
        let score = if m0.abs() > 1e-4 * area_tot.sqrt() {
            m0
        } else if m1.abs() > 1e-5 * area_tot.sqrt() {
            m1
        } else {
            peak_val
        };
        if score < 0.0 {
            inv_norm = -inv_norm;
            for cq in &mut fields {
                for e3 in cq {
                    e3[0] = -e3[0];
                    e3[1] = -e3[1];
                    e3[2] = -e3[2];
                }
            }
        }
        (fields, score.abs(), inv_norm)
    };

    let solve_2d_at = |freq: f64| -> Result<Vec<(f64, Vec<Vec<[f64; 3]>>, Vec<f64>, f64)>, String> {
        let asm2 = assemble_modes_with_bcs(&mesh_2d, freq, &eps_2d, &mur_2d, mode_order, &[1], &[], &[2])?;
        if asm2.free.len() < 3 {
            return Err("too few free DOFs on port cross-section".into());
        }
        let k_actual = k_ask.min(asm2.free.len().saturating_sub(2)).max(1);
        let modes = solve_mode_pairs(&asm2, n_search, false, k_actual)?;
        if modes.neff.is_empty() {
            return Err(format!("waveport '{name}': no 2D modes converged"));
        }
        let k0 = 2.0 * PI * freq / C0;
        let mut cands = Vec::with_capacity(modes.neff.len());
        for (mi, &(n_re, n_im)) in modes.neff.iter().enumerate() {
            let beta2 = if n_re > 0.0 { (n_re * k0).powi(2) } else { -(n_im * k0).powi(2) };
            let (ef, score, scale) = eval_norm_mode(&modes.et[mi]);
            cands.push((beta2, score, ef, modes.et[mi].clone(), scale));
        }
        if let Some(ng) = n_eff_guess {
            cands.sort_by(|a, b| match (a.0 > 0.0, b.0 > 0.0) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                (true, true) => {
                    let da = (a.0.sqrt() / k0 - ng).abs();
                    let db = (b.0.sqrt() / k0 - ng).abs();
                    if (da - db).abs() > 1e-6 * ng.max(1.0) {
                        da.total_cmp(&db)
                    } else {
                        b.1.total_cmp(&a.1)
                    }
                }
                (false, false) => b.0.total_cmp(&a.0),
            });
        } else {
            cands.sort_by(|a, b| {
                let scale = a.0.abs().max(b.0.abs()).max(1.0);
                if (a.0 - b.0).abs() > 1e-5 * scale {
                    b.0.total_cmp(&a.0)
                } else {
                    b.1.total_cmp(&a.1)
                }
            });
        }
        Ok(cands.into_iter().map(|(b2, _, ef, et, sc)| (b2, ef, et, sc)).collect())
    };

    let Ok(mid_modes) = solve_2d_at(f_mid) else {
        return Ok(make_tem_fallback(v_hat));
    };
    if mode_idx > mid_modes.len() || mid_modes[mode_idx - 1].0 <= 0.0 {
        return Ok(make_tem_fallback(v_hat));
    }
    let (beta2_mid, ref_ef, ref_et, ref_scale) = mid_modes[mode_idx - 1].clone();

    let overlap_with_ref = |ef: &[Vec<[f64; 3]>]| -> f64 {
        let mut acc = 0.0;
        for (c_idx, (c1, c2)) in ef.iter().zip(&ref_ef).enumerate() {
            let imu = 1.0 / mur_2d[c_idx];
            for (q, (&v1, &v2)) in c1.iter().zip(c2).enumerate() {
                acc += tri_weights[c_idx][q] * imu * dot3(v1, v2);
            }
        }
        acc.abs()
    };

    let track_b2 = |freq: f64| -> f64 {
        let Ok(cands) = solve_2d_at(freq) else { return beta2_mid };
        let mut best_ov = -1.0;
        let mut best_b2 = beta2_mid;
        for (b2, ef, _, _) in cands {
            let ov = overlap_with_ref(&ef);
            if ov > best_ov {
                best_ov = ov;
                best_b2 = b2;
            }
        }
        best_b2
    };

    let beta2_poly = if (f_max - f_min).abs() <= 1e-9 * f_mid {
        let s_mid = (2.0 * PI * f_mid / C0).powi(2);
        [0.0, beta2_mid / s_mid.max(1e-300), 0.0]
    } else {
        let b2_min = track_b2(f_min);
        let b2_max = track_b2(f_max);
        let s0 = (2.0 * PI * f_min / C0).powi(2);
        let s1 = (2.0 * PI * f_mid / C0).powi(2);
        let s2 = (2.0 * PI * f_max / C0).powi(2);
        let w0 = b2_min / ((s0 - s1) * (s0 - s2));
        let w1 = beta2_mid / ((s1 - s0) * (s1 - s2));
        let w2 = b2_max / ((s2 - s0) * (s2 - s1));
        [w0 * s1 * s2 + w1 * s0 * s2 + w2 * s0 * s1, -(w0 * (s1 + s2) + w1 * (s0 + s2) + w2 * (s0 + s1)), w0 + w1 + w2]
    };

    let omega_0 = 2.0 * PI * f0;
    let s_0 = (omega_0 / C0).powi(2);
    let beta_0 = (beta2_poly[0] + beta2_poly[1] * s_0 + beta2_poly[2] * s_0 * s_0).max(1e-12).sqrt();
    let y0_carrier = beta_0 / (omega_0 * MU0);

    Ok(SolvedPortProfile {
        port: DgtdPort {
            name: name.to_string(),
            labels: labels.to_vec(),
            n_out,
            u_hat,
            v_hat,
            beta2_poly,
            z0_override,
            y0_carrier,
            is_lumped: false,
        },
        mode_2d: Some((mesh_2d, mode_order, ref_et, ref_scale)),
        e_const: [0.0; 3],
    })
}

/// Evaluates whether physical point `p` lies inside cell `ci` and returns its reference coordinates `xi`.
fn locate_point_in_cell(mesh: &MixedDiscretisation, ci: usize, p: [f64; 3]) -> Option<[f64; 3]> {
    let ct = mesh.cell_type(ci);
    let iso = IsoparametricCell::from_mesh(mesh, ci);
    let mut xi = match ct {
        CellType::Tri => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        CellType::Tet => [0.25, 0.25, 0.25],
        CellType::Prism => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        CellType::Quad | CellType::Hex => [0.0, 0.0, 0.0],
    };
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    for _ in 0..12 {
        let jac = iso.eval_into(&xi, &mut sn, &mut sdn);
        let r = sub3(p, jac.x);
        if (mesh.dim == 2 && r[0].hypot(r[1]) < 1e-12) || (mesh.dim == 3 && norm3(r) < 1e-12) {
            break;
        }
        // dxi = J^-1 r = (J^-T)^T r
        for k in 0..mesh.dim {
            let mut dxi_k = 0.0;
            for i in 0..mesh.dim {
                dxi_k += jac.j_inv_t[i][k] * r[i];
            }
            xi[k] += dxi_k;
        }
    }
    let tol = 1e-5;
    let inside = match ct {
        CellType::Tri => xi[0] >= -tol && xi[1] >= -tol && xi[0] + xi[1] <= 1.0 + tol,
        CellType::Quad => xi[0].abs() <= 1.0 + tol && xi[1].abs() <= 1.0 + tol,
        CellType::Tet => xi[0] >= -tol && xi[1] >= -tol && xi[2] >= -tol && xi[0] + xi[1] + xi[2] <= 1.0 + tol,
        CellType::Prism => xi[0] >= -tol && xi[1] >= -tol && xi[0] + xi[1] <= 1.0 + tol && xi[2].abs() <= 1.0 + tol,
        CellType::Hex => xi[0].abs() <= 1.0 + tol && xi[1].abs() <= 1.0 + tol && xi[2].abs() <= 1.0 + tol,
    };
    inside.then_some(xi)
}

// ---------------------------------------------------------------------------
// Stage 1: assemble_dgtd
// ---------------------------------------------------------------------------

/// Stage 1 of the DGTD solver (`ADR 9`, `ADR 144`): constructs local high-order polynomial bases ($p \in \{1, 2, 3\}$),
/// pre-inverts block-local element mass matrices $M_K^{-1}$, builds strong-form SBP curl operators, resolves
/// waveguide modal ports and zero-thickness `Sheet2D` (`FR-35`) interface jump tensors, and partitions elements
/// into power-of-two Local Time-Stepping (LTS) CFL tiers.
pub fn assemble_dgtd(mesh: &Mesh, problem: &Problem) -> Result<DgtdAssembly, String> {
    let Study::Dgtd(ref config) = problem.study else {
        return Err(format!("assemble_dgtd requires Study::Dgtd, got '{}'", problem.study.kind()));
    };
    #[cfg(not(feature = "gpu"))]
    if config.backend == crate::problem::TdBackend::Gpu {
        return Err(
            "DGTD backend 'gpu' requested, but nyx-ems was compiled without GPU support; recompile with `--features gpu`"
                .into(),
        );
    }
    let basis_order = config.basis_order.unwrap_or(problem.order).clamp(1, 3);
    let mixed = discretise_mixed(mesh).map_err(|e| e.to_string())?;
    let dim = mixed.dim;
    let nc = mixed.n_cells();
    let gll = Gll::new(basis_order);
    let f_ref = config.reference_freq().max(1.0);
    let omega_ref = 2.0 * PI * f_ref;

    // 1. Evaluate cell material tensors (permittivity, permeability, conductivity, gyrotropy, thermo-optic/EO/PE)
    let mut cell_eps_re = vec![[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]; nc];
    let mut cell_eps_im = vec![[[0.0; 3]; 3]; nc];
    let mut cell_mu_re = vec![[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]; nc];
    let mut cell_sigma = vec![[[0.0; 3]; 3]; nc];
    let mut cell_mur_scalar = vec![1.0; nc];

    let eo_pe = problem.domain_eo_pe_global();
    let cell_e_dc = problem.eval_cell_electric_field_with(nc, mixed.n_points(), dim, |ci| mixed.cell(ci));
    let cell_s_dc = problem.eval_cell_strain_with(nc, mixed.n_points(), |ci| mixed.cell(ci));

    for ci in 0..nc {
        let dom = mixed.cell_labels[ci];
        let mut eps_re = if let Some(an) = problem.materials.aniso.get(&dom) {
            if let Some(ref t) = an.permittivity {
                femcore::material::tensor_global(t, an.euler_deg)
            } else {
                let e = problem.materials.get(dom, "eps_r").unwrap_or(1.0);
                [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]]
            }
        } else {
            let e = problem.materials.get(dom, "eps_r").unwrap_or(1.0);
            [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]]
        };
        let mut eps_im = if let Some(an) = problem.materials.aniso.get(&dom) {
            if let Some(ref t) = an.permittivity_im {
                femcore::material::tensor_global(t, an.euler_deg)
            } else {
                let e_im = problem
                    .materials
                    .get(dom, "eps_r_im")
                    .unwrap_or_else(|| eps_re[0][0] * problem.materials.get(dom, "loss_tangent").unwrap_or(0.0));
                [[e_im, 0.0, 0.0], [0.0, e_im, 0.0], [0.0, 0.0, e_im]]
            }
        } else {
            let e_im = problem
                .materials
                .get(dom, "eps_r_im")
                .unwrap_or_else(|| eps_re[0][0] * problem.materials.get(dom, "loss_tangent").unwrap_or(0.0));
            [[e_im, 0.0, 0.0], [0.0, e_im, 0.0], [0.0, 0.0, e_im]]
        };

        // Thermo-optic evaluation if active
        if let Some(td) = problem.thermo_optic_domains.get(&dom) {
            let tc =
                problem.cell_temperatures.as_ref().and_then(|v| v.get(ci).copied()).or(problem.fallback_temperature);
            if let Some(tk) = tc.filter(|&t| t > 0.0 && t.is_finite()) {
                if let Some((re, im)) = td.eval_at(tk)? {
                    eps_re = re;
                    eps_im = im;
                }
            }
        }
        // Electro-optic / photoelastic perturbation if active
        if let Some((r_g, p_g)) = eo_pe.get(&dom) {
            let e_dc = cell_e_dc.as_ref().map(|v| &v[ci]);
            let s_dc = cell_s_dc.as_ref().map(|v| &v[ci]);
            let d_eps = femcore::material::coupled::delta_eps_eo_pe(&eps_re, r_g.as_ref(), p_g.as_ref(), e_dc, s_dc);
            for a in 0..3 {
                for b in 0..3 {
                    eps_re[a][b] += d_eps[a][b];
                }
            }
        }

        let mur = problem.materials.get(dom, "mu_r").unwrap_or(1.0);
        let mu_re = if let Some(an) = problem.materials.aniso.get(&dom) {
            if let Some(ref t) = an.permeability {
                femcore::material::tensor_global(t, an.euler_deg)
            } else {
                [[mur, 0.0, 0.0], [0.0, mur, 0.0], [0.0, 0.0, mur]]
            }
        } else {
            [[mur, 0.0, 0.0], [0.0, mur, 0.0], [0.0, 0.0, mur]]
        };

        let sigma_scalar = problem.materials.get(dom, "sigma").unwrap_or(0.0);
        let mut sig = [[0.0; 3]; 3];
        for a in 0..3 {
            sig[a][a] += sigma_scalar;
            for b in 0..3 {
                sig[a][b] += omega_ref * EPS0 * eps_im[a][b];
            }
        }
        // Gyrotropic skew-symmetric off-diagonal coupling g x E (FR-34)
        if let Some(an) = problem.materials.aniso.get(&dom) {
            if let Some(g) = an.gyrotropy {
                let g_rot = mat3_mul_vec(&femcore::material::rotation(an.euler_deg), g);
                // (g x E) = [g_y E_z - g_z E_y, g_z E_x - g_x E_z, g_x E_y - g_y E_x]
                let scale = omega_ref * EPS0;
                sig[0][1] -= scale * g_rot[2];
                sig[0][2] += scale * g_rot[1];
                sig[1][0] += scale * g_rot[2];
                sig[1][2] -= scale * g_rot[0];
                sig[2][0] -= scale * g_rot[1];
                sig[2][1] += scale * g_rot[0];
            }
        }

        cell_eps_re[ci] = eps_re;
        cell_eps_im[ci] = eps_im;
        cell_mu_re[ci] = mu_re;
        cell_mur_scalar[ci] = ((mu_re[0][0] + mu_re[1][1] + mu_re[2][2]) / 3.0).max(1e-6);
        cell_sigma[ci] = sig;
    }

    // 2. Collect dispersive Drude/Lorentz/gyrotropic poles and opt-in `[nonlinear]` constitutive blocks per domain
    let domain_poles: BTreeMap<i32, Vec<DgtdDispersionPole>> = config.dispersion_poles.clone();
    let mut domain_nl: BTreeMap<i32, (femcore::material::Nonlinear, [f64; 3])> = BTreeMap::new();
    if config.nonlinear {
        for (&dom, an) in &problem.materials.aniso {
            if let Some(ref nl) = an.nonlinear {
                if let Some(ref q) = nl.qpm {
                    q.validate()?;
                }
                if let Some(ref c3) = nl.chi3 {
                    c3.validate()?;
                }
                domain_nl.insert(dom, (nl.clone(), an.euler_deg));
            }
        }
        for (&dom, td) in &problem.thermo_optic_domains {
            if let std::collections::btree_map::Entry::Vacant(e) = domain_nl.entry(dom) {
                if let Some(ref nl) = td.spec.nonlinear {
                    if let Some(ref q) = nl.qpm {
                        q.validate()?;
                    }
                    if let Some(ref c3) = nl.chi3 {
                        c3.validate()?;
                    }
                    e.insert((nl.clone(), td.euler_deg));
                }
            }
        }
    }

    // 3. Bounding boxes of PML domains for uniaxial PML damping rate calculation
    let mut pml_boxes: BTreeMap<i32, ([f64; 3], [f64; 3])> = BTreeMap::new();
    for ci in 0..nc {
        let dom = mixed.cell_labels[ci];
        if let Some(pml) = problem.materials.pml.get(&dom) {
            let entry = pml_boxes.entry(dom).or_insert(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
            if let Some(b) = pml.bounds {
                *entry = b;
            } else {
                for &v in mixed.cell(ci) {
                    let p = mixed.point3(v as usize);
                    for k in 0..3 {
                        entry.0[k] = entry.0[k].min(p[k]);
                        entry.1[k] = entry.1[k].max(p[k]);
                    }
                }
            }
        }
    }

    // 4. Build face connectivity map: canonical FaceKey -> Vec<(cell_id, local_face_index)>
    let mut face_owner_map: HashMap<FaceKey, Vec<(usize, usize)>> = HashMap::with_capacity(nc * 3);
    let mut cell_local_faces: Vec<Vec<(FacetType, FaceKey)>> = Vec::with_capacity(nc);
    for ci in 0..nc {
        let lfs = cell_local_face_vertices(&mixed, ci);
        let mut cell_keys = Vec::with_capacity(lfs.len());
        for (lf_idx, (ft, fv)) in lfs.into_iter().enumerate() {
            let key = canonical_face_key(ft, &fv);
            face_owner_map.entry(key).or_default().push((ci, lf_idx));
            cell_keys.push((ft, key));
        }
        cell_local_faces.push(cell_keys);
    }

    // Map canonical FaceKey -> boundary facet label (if any)
    let mut face_label_map: HashMap<FaceKey, (i32, usize)> = HashMap::with_capacity(mixed.n_facets());
    for f in 0..mixed.n_facets() {
        let ft = mixed.facet_type(f);
        let fv = mixed.facet(f);
        let key = canonical_face_key(ft, fv);
        face_label_map.insert(key, (mixed.facet_labels[f], f));
    }

    // 5. Resolve zero-thickness 2D material sheets (FR-35, Sheet2dBc) onto canonical FaceKeys
    let mut face_sheet2d_map: HashMap<FaceKey, ([[f64; 2]; 2], [f64; 3], [f64; 3])> = HashMap::new();
    for sh in &problem.sheet_2d {
        let eval_sheet_at_verts =
            |verts: &[u32], n_raw: [f64; 3]| -> Result<([[f64; 2]; 2], [f64; 3], [f64; 3]), String> {
                let local_pot = problem.imported_potential.as_ref().and_then(|pot| {
                    if pot.is_empty() {
                        return None;
                    }
                    let mut sum = 0.0;
                    let mut cnt = 0usize;
                    for &v in verts {
                        if let Some(&pv) = pot.get(v as usize) {
                            sum += pv;
                            cnt += 1;
                        }
                    }
                    (cnt > 0).then(|| sum / (cnt as f64))
                });
                let local_temp = problem
                    .temperature_field
                    .as_ref()
                    .and_then(|tf| {
                        let mut sum = 0.0;
                        let mut cnt = 0usize;
                        for &v in verts {
                            if let Some(&tv) = tf.get(v as usize) {
                                sum += tv;
                                cnt += 1;
                            }
                        }
                        (cnt > 0).then(|| sum / (cnt as f64))
                    })
                    .or(problem.fallback_temperature);
                let s_c = sh.eval_tensor_at(f_ref, local_temp, local_pot)?;
                let s_re = [[s_c[0][0].re.max(0.0), s_c[0][1].re], [s_c[1][0].re, s_c[1][1].re.max(0.0)]];
                let (u_hat, v_hat, _) = if dim == 2 {
                    let u_raw = [-n_raw[1], n_raw[0], 0.0];
                    let ul = norm3(u_raw).max(1e-30);
                    ([u_raw[0] / ul, u_raw[1] / ul, 0.0], [0.0, 0.0, 1.0], n_raw)
                } else {
                    sheet2d_tangent_frame_3d(n_raw)
                };
                Ok((s_re, u_hat, v_hat))
            };

        for f in 0..mixed.n_facets() {
            if sh.labels.contains(&mixed.facet_labels[f]) {
                let ft = mixed.facet_type(f);
                let fv = mixed.facet(f);
                let key = canonical_face_key(ft, fv);
                let (_, _, n_c) = eval_canonical_face_geometry(&mixed, &key, &[0.0, 0.0]);
                let entry = eval_sheet_at_verts(fv, n_c)?;
                face_sheet2d_map.insert(key, entry);
            }
        }
        // Also support domain-pair interface selection `a|b` or `interface(a, b)`
        if !sh.interface_domains.is_empty() {
            for (&key, owners) in &face_owner_map {
                if owners.len() == 2 {
                    let d0 = mixed.cell_labels[owners[0].0];
                    let d1 = mixed.cell_labels[owners[1].0];
                    let matches = sh.interface_domains.iter().any(|(da, db)| {
                        (da.contains(&d0) && db.contains(&d1)) || (da.contains(&d1) && db.contains(&d0))
                    });
                    if matches {
                        let (_, _, n_c) = eval_canonical_face_geometry(&mixed, &key, &[0.0, 0.0]);
                        let entry = eval_sheet_at_verts(canonical_face_verts(&key), n_c)?;
                        face_sheet2d_map.insert(key, entry);
                    }
                }
            }
        }
    }

    // 6. Solve waveguide modal ports and lumped ports
    let mut solved_ports: Vec<SolvedPortProfile> = Vec::new();
    for wp in &problem.waveports {
        solved_ports.push(build_dgtd_port_profile(
            &mixed,
            problem,
            config,
            &wp.name,
            &wp.labels,
            wp.mode,
            wp.n_eff_guess,
            wp.z0_ohm,
            false,
            None,
            &cell_eps_re,
            &cell_mur_scalar,
            &face_owner_map,
        )?);
    }
    for lp in &problem.lumped_ports {
        solved_ports.push(build_dgtd_port_profile(
            &mixed,
            problem,
            config,
            &lp.name,
            &lp.labels,
            1,
            None,
            Some(lp.z0),
            true,
            Some(lp.direction),
            &cell_eps_re,
            &cell_mur_scalar,
            &face_owner_map,
        )?);
    }

    let (sym_pec_facets, sym_pmc_facets) = symmetry_facet_flags_mixed(&mixed, problem.symmetry)?;

    // Check for TFSF plane-wave excitation
    let tfsf_plane_wave: Option<(&ExcitationSpec, [f64; 3], [f64; 3], [f64; 3], f64, f64)> =
        problem.excitations.iter().find_map(|spec| match &spec.excitation {
            Excitation::PlaneWave { k_hat, pol, amplitude, origin, .. } => {
                let kl = norm3(*k_hat).max(1e-30);
                let kh = [k_hat[0] / kl, k_hat[1] / kl, k_hat[2] / kl];
                let p_re = [pol[0].re, pol[1].re, pol[2].re];
                let pdot = dot3(p_re, kh);
                let pt = [p_re[0] - pdot * kh[0], p_re[1] - pdot * kh[1], p_re[2] - pdot * kh[2]];
                let pl = norm3(pt).max(1e-30);
                let p_vec = [amplitude * pt[0] / pl, amplitude * pt[1] / pl, amplitude * pt[2] / pl];
                let n_bg = (spec.eps_bg.re.max(1.0) * spec.mu_bg.max(1.0)).sqrt();
                let c_bg = C0 / n_bg;
                let y_bg = (spec.eps_bg.re.max(1.0) / spec.mu_bg.max(1.0)).sqrt() / ETA0;
                Some((spec, kh, p_vec, *origin, c_bg, y_bg))
            }
            _ => None,
        });

    // 7. Assemble each DG element K: basis, M_K^-1, faces, SBP D^x, D^y, D^z, and CFL time step
    let mut cells = Vec::with_capacity(nc);
    let mut total_dofs = 0usize;
    let mut total_pole_states = 0usize;
    let mut dt_min = f64::INFINITY;
    let mut dt_max = 0.0_f64;
    let mut h_min = f64::INFINITY;
    let mut h_max = 0.0_f64;

    let mut vals_buf = Vec::new();
    let mut dxi_buf = Vec::new();
    let mut sn_buf = Vec::new();
    let mut sdn_buf = Vec::new();

    for ci in 0..nc {
        let ct = mixed.cell_type(ci);
        let dom = mixed.cell_labels[ci];
        let n_p = dg_nloc(ct, basis_order);
        let dof_offset = total_dofs;
        total_dofs += n_p;

        let poles = domain_poles.get(&dom).cloned().unwrap_or_default();
        let pole_offset = total_pole_states;
        total_pole_states += poles.len() * n_p;

        let iso = IsoparametricCell::from_mesh(&mixed, ci);
        let ref_nodes = dg_ref_nodes(ct, basis_order, &gll);
        let mut node_coords = Vec::with_capacity(n_p);
        for xi_n in &ref_nodes {
            node_coords.push(iso.eval_into(xi_n, &mut sn_buf, &mut sdn_buf).x);
        }

        // Volume mass matrix M_K and volume stiffness integrals K^{vol, k}
        let is_diagonal_mass = matches!(ct, CellType::Quad | CellType::Hex);
        let mut mass = Vec::new();
        let mut inv_mass = Vec::new();
        let mut k_vol_x = vec![0.0; n_p * n_p];
        let mut k_vol_y = vec![0.0; n_p * n_p];
        let mut k_vol_z = vec![0.0; n_p * n_p];
        let mut measure = 0.0_f64;

        if is_diagonal_mass {
            // Diagonal GLL quadrature on Quad (2D) or Hex (3D) via femcore::linalg::matfree::Gll
            mass.resize(n_p, 0.0);
            inv_mass.resize(n_p, 0.0);
            let n1 = gll.nodes.len();
            for a in 0..n_p {
                let ia = a % n1;
                let ib = (a / n1) % n1;
                let ic = if dim == 3 { a / (n1 * n1) } else { 0 };
                let w_ref = gll.weights[ia] * gll.weights[ib] * if dim == 3 { gll.weights[ic] } else { 1.0 };
                let jac = iso.eval_into(&ref_nodes[a], &mut sn_buf, &mut sdn_buf);
                let det_abs = jac.det.abs();
                if !(det_abs > 1e-30 && det_abs.is_finite()) {
                    return Err(format!("degenerate element {ci} in DGTD assembly (det = {})", jac.det));
                }
                let m_a = w_ref * det_abs;
                mass[a] = m_a;
                inv_mass[a] = 1.0 / m_a;
                measure += m_a;

                dg_shape_into(ct, basis_order, &gll, &ref_nodes[a], &mut vals_buf, &mut dxi_buf);
                for b in 0..n_p {
                    let grad_b = jac.map_grad(dxi_buf[b]);
                    k_vol_x[a * n_p + b] = m_a * grad_b[0];
                    k_vol_y[a * n_p + b] = m_a * grad_b[1];
                    k_vol_z[a * n_p + b] = m_a * grad_b[2];
                }
            }
        } else {
            // Exact quadrature and local dense mass matrix factorization on Tri, Tet, CurvedTet, Prism
            let q_deg = 2 * basis_order + if iso.curved { 2 } else { 0 };
            let rule = quad_rule_mixed(ct, q_deg);
            let mut m_dense = vec![0.0; n_p * n_p];
            let curved_tet = (ct == CellType::Tet && iso.curved).then(|| {
                let mut x10 = [[0.0; 3]; 10];
                for (idx, p3) in iso.nodes.iter().take(10).enumerate() {
                    x10[idx] = *p3;
                }
                CurvedTet::new(x10)
            });

            for (xi_q, &w_q) in rule.pts.iter().zip(&rule.w) {
                let jac = iso.eval_into(xi_q, &mut sn_buf, &mut sdn_buf);
                let det_abs = if let Some(ref ct_iso) = curved_tet {
                    let lam = [1.0 - xi_q[0] - xi_q[1] - xi_q[2], xi_q[0], xi_q[1], xi_q[2]];
                    ct_iso.jacobian(&lam).1.abs()
                } else {
                    jac.det.abs()
                };
                if !(det_abs > 1e-30 && det_abs.is_finite()) {
                    return Err(format!("degenerate element {ci} in DGTD assembly (det = {det_abs})"));
                }
                let weight = w_q * det_abs;
                measure += weight;
                dg_shape_into(ct, basis_order, &gll, xi_q, &mut vals_buf, &mut dxi_buf);
                for a in 0..n_p {
                    let wa = weight * vals_buf[a];
                    for b in 0..n_p {
                        m_dense[a * n_p + b] += wa * vals_buf[b];
                        let grad_b = jac.map_grad(dxi_buf[b]);
                        k_vol_x[a * n_p + b] += wa * grad_b[0];
                        k_vol_y[a * n_p + b] += wa * grad_b[1];
                        k_vol_z[a * n_p + b] += wa * grad_b[2];
                    }
                }
            }
            inv_mass = invert_dense_matrix(&m_dense, n_p)?;
            mass = m_dense;
        }

        // Assemble local faces and symmetric SBP surface matrices F^x, F^y, F^z
        let mut f_mat_x = vec![0.0; n_p * n_p];
        let mut f_mat_y = vec![0.0; n_p * n_p];
        let mut f_mat_z = vec![0.0; n_p * n_p];
        let mut perimeter_area = 0.0_f64;
        let mut faces = Vec::with_capacity(cell_local_faces[ci].len());

        let eps_r = cell_eps_re[ci];
        let mu_r = cell_mu_re[ci];
        let eps_iso = ((eps_r[0][0] + eps_r[1][1] + eps_r[2][2]) / 3.0).max(1e-6);
        let mu_iso = cell_mur_scalar[ci];
        let y_minus = (EPS0 * eps_iso / (MU0 * mu_iso)).sqrt();

        for &(ft, key) in &cell_local_faces[ci] {
            let Some(owners) = face_owner_map.get(&key) else { continue };
            let neighbor_cell = owners.iter().find(|&&(c_o, _)| c_o != ci).map(|&(c_o, _)| c_o);
            let (facet_label, facet_idx_opt) =
                face_label_map.get(&key).map(|&(l, fi)| (l, Some(fi))).unwrap_or((0, None));

            let port_index =
                solved_ports.iter().position(|sp| facet_label > 0 && sp.port.labels.contains(&facet_label));

            let kind = if let Some(nb) = neighbor_cell {
                if facet_label > 0 && problem.pec.contains(&facet_label) {
                    DgtdFaceKind::Pec
                } else {
                    DgtdFaceKind::Interior { neighbor_cell: nb }
                }
            } else if let Some(pi) = port_index {
                DgtdFaceKind::Port { port_index: pi }
            } else if facet_label > 0 && config.absorbing_labels.contains(&facet_label) {
                DgtdFaceKind::SilverMuller
            } else if (facet_label > 0 && problem.pmc.contains(&facet_label))
                || facet_idx_opt.is_some_and(|fi| sym_pmc_facets[fi])
            {
                DgtdFaceKind::Pmc
            } else if tfsf_plane_wave.is_some()
                && facet_label > 0
                && !problem.pec.contains(&facet_label)
                && facet_idx_opt.is_none_or(|fi| !sym_pec_facets[fi])
            {
                // When a plane_wave excitation is active without explicit ports/PML, outer non-PEC boundaries absorb
                DgtdFaceKind::SilverMuller
            } else {
                DgtdFaceKind::Pec
            };

            let y_plus = match kind {
                DgtdFaceKind::Interior { neighbor_cell: nb } => {
                    let e_nb =
                        ((cell_eps_re[nb][0][0] + cell_eps_re[nb][1][1] + cell_eps_re[nb][2][2]) / 3.0).max(1e-6);
                    let m_nb = cell_mur_scalar[nb];
                    (EPS0 * e_nb / (MU0 * m_nb)).sqrt()
                }
                DgtdFaceKind::Port { port_index: pi } => {
                    let y0_p = solved_ports[pi].port.y0_carrier;
                    if solved_ports[pi].port.is_lumped {
                        y0_p
                    } else {
                        y0_p / mu_iso
                    }
                }
                DgtdFaceKind::SilverMuller | DgtdFaceKind::Pec | DgtdFaceKind::Pmc => y_minus,
            };

            let sheet_2d = face_sheet2d_map.get(&key).copied();
            let (f_pts, f_wts) = canonical_face_quadrature(ft, basis_order, &gll);
            let mut qps = Vec::with_capacity(f_pts.len());

            // Check TFSF interface status on this face
            let tfsf_face_info = tfsf_plane_wave.and_then(|(spec, kh, p_vec, origin, c_bg, y_bg)| {
                let in_tf_minus = spec.domains.is_empty() || spec.domains.contains(&dom);
                match kind {
                    DgtdFaceKind::Interior { neighbor_cell: nb } => {
                        let dom_plus = mixed.cell_labels[nb];
                        let in_tf_plus = spec.domains.is_empty() || spec.domains.contains(&dom_plus);
                        if in_tf_minus != in_tf_plus {
                            let sign = if in_tf_minus { 1.0 } else { -1.0 };
                            Some((sign, kh, p_vec, origin, c_bg, y_bg))
                        } else {
                            None
                        }
                    }
                    DgtdFaceKind::SilverMuller if in_tf_minus => Some((1.0, kh, p_vec, origin, c_bg, y_bg)),
                    _ => None,
                }
            });

            for (xi_f, &w_f) in f_pts.iter().zip(&f_wts) {
                let (x_q, ds_q, mut n_q) = eval_canonical_face_geometry(&mixed, &key, xi_f);
                let (xi_minus, v_out_ref) = map_face_to_cell_ref(&mixed, ci, &key, xi_f);
                let jac_m = iso.eval_into(&xi_minus, &mut sn_buf, &mut sdn_buf);
                let v_out_phys = if dim == 2 {
                    [
                        jac_m.j[0][0] * v_out_ref[0] + jac_m.j[0][1] * v_out_ref[1],
                        jac_m.j[1][0] * v_out_ref[0] + jac_m.j[1][1] * v_out_ref[1],
                        0.0,
                    ]
                } else {
                    mat3_mul_vec(&jac_m.j, v_out_ref)
                };
                if dot3(n_q, v_out_phys) < 0.0 {
                    n_q = [-n_q[0], -n_q[1], -n_q[2]];
                }
                let w_ds = w_f * ds_q;
                perimeter_area += w_ds;

                dg_shape_into(ct, basis_order, &gll, &xi_minus, &mut vals_buf, &mut dxi_buf);
                let phi_minus = vals_buf.clone();
                for a in 0..n_p {
                    if phi_minus[a].abs() <= 1e-15 {
                        continue;
                    }
                    let wa = w_ds * phi_minus[a];
                    for b in 0..n_p {
                        if phi_minus[b].abs() <= 1e-15 {
                            continue;
                        }
                        let wab = wa * phi_minus[b];
                        f_mat_x[a * n_p + b] += wab * n_q[0];
                        f_mat_y[a * n_p + b] += wab * n_q[1];
                        f_mat_z[a * n_p + b] += wab * n_q[2];
                    }
                }

                let trace_minus: Vec<(usize, f64)> =
                    phi_minus.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 1e-13).collect();

                let mut lift_minus = Vec::new();
                if is_diagonal_mass {
                    for &(b, pb) in &trace_minus {
                        let lv = inv_mass[b] * w_ds * pb;
                        if lv.abs() > 1e-15 {
                            lift_minus.push((b, lv));
                        }
                    }
                } else {
                    for a in 0..n_p {
                        let mut acc = 0.0;
                        for &(b, pb) in &trace_minus {
                            acc += inv_mass[a * n_p + b] * pb;
                        }
                        let lv = acc * w_ds;
                        if lv.abs() > 1e-15 {
                            lift_minus.push((a, lv));
                        }
                    }
                }

                let trace_plus = if let DgtdFaceKind::Interior { neighbor_cell: nb } = kind {
                    let ct_nb = mixed.cell_type(nb);
                    let (xi_plus, _) = map_face_to_cell_ref(&mixed, nb, &key, xi_f);
                    dg_shape_into(ct_nb, basis_order, &gll, &xi_plus, &mut vals_buf, &mut dxi_buf);
                    vals_buf.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 1e-13).collect()
                } else {
                    Vec::new()
                };

                let port_mode = if let DgtdFaceKind::Port { port_index: pi } = kind {
                    Some(solved_ports[pi].eval_at(x_q))
                } else {
                    None
                };

                qps.push(DgtdFaceQp {
                    x: x_q,
                    w_ds,
                    n_hat: n_q,
                    trace_minus,
                    trace_plus,
                    lift_minus,
                    port_mode,
                    tfsf: tfsf_face_info,
                });
            }

            faces.push(DgtdCellFace { kind, y_plus, sheet_2d, facet_label, qps });
        }

        // Form SBP skew-symmetric strong-form differentiation matrices D^k = M_K^-1 \hat{K}^k
        let mut dx = vec![0.0; n_p * n_p];
        let mut dy = vec![0.0; n_p * n_p];
        let mut dz = vec![0.0; n_p * n_p];
        let mut k_hat_x = vec![0.0; n_p * n_p];
        let mut k_hat_y = vec![0.0; n_p * n_p];
        let mut k_hat_z = vec![0.0; n_p * n_p];
        for a in 0..n_p {
            for b in 0..n_p {
                k_hat_x[a * n_p + b] = 0.5 * (k_vol_x[a * n_p + b] - k_vol_x[b * n_p + a]) + 0.5 * f_mat_x[a * n_p + b];
                k_hat_y[a * n_p + b] = 0.5 * (k_vol_y[a * n_p + b] - k_vol_y[b * n_p + a]) + 0.5 * f_mat_y[a * n_p + b];
                k_hat_z[a * n_p + b] = 0.5 * (k_vol_z[a * n_p + b] - k_vol_z[b * n_p + a]) + 0.5 * f_mat_z[a * n_p + b];
            }
        }
        if is_diagonal_mass {
            for a in 0..n_p {
                let im = inv_mass[a];
                for b in 0..n_p {
                    dx[a * n_p + b] = im * k_hat_x[a * n_p + b];
                    dy[a * n_p + b] = im * k_hat_y[a * n_p + b];
                    dz[a * n_p + b] = im * k_hat_z[a * n_p + b];
                }
            }
        } else {
            for a in 0..n_p {
                for b in 0..n_p {
                    let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
                    for q in 0..n_p {
                        let im_aq = inv_mass[a * n_p + q];
                        sx += im_aq * k_hat_x[q * n_p + b];
                        sy += im_aq * k_hat_y[q * n_p + b];
                        sz += im_aq * k_hat_z[q * n_p + b];
                    }
                    dx[a * n_p + b] = sx;
                    dy[a * n_p + b] = sy;
                    dz[a * n_p + b] = sz;
                }
            }
        }

        // Uniaxial PML damping rate [sigma_x, sigma_y, sigma_z] at element centroid
        let mut pml_sigma = [0.0; 3];
        let mut in_pml = false;
        if let (Some(pml), Some(&(lo, hi))) = (problem.materials.pml.get(&dom), pml_boxes.get(&dom)) {
            in_pml = true;
            let mut cent = [0.0; 3];
            for p3 in &node_coords {
                for k in 0..3 {
                    cent[k] += p3[k] / (n_p as f64);
                }
            }
            let dirs = if dim == 2 { [pml.x, pml.z, PmlDir::Off] } else { [pml.x, pml.y, pml.z] };
            for axis in 0..dim {
                let d_len = (hi[axis] - lo[axis]).max(1e-12);
                let depth = match dirs[axis] {
                    PmlDir::Plus => ((cent[axis] - lo[axis]) / d_len).clamp(0.0, 1.0),
                    PmlDir::Minus => ((hi[axis] - cent[axis]) / d_len).clamp(0.0, 1.0),
                    PmlDir::Off => 0.0,
                };
                if depth > 0.0 {
                    let s_max = (pml.order + 1.0) * pml.strength.max(1.0) * C0 / (2.0 * d_len * eps_iso.sqrt());
                    pml_sigma[axis] = s_max * depth.powf(pml.order.max(1.0));
                }
            }
        }

        let mut inv_eps = invert_3x3(&eps_r)?;
        let mut inv_mu = invert_3x3(&mu_r)?;
        for a in 0..3 {
            for b in 0..3 {
                inv_eps[a][b] /= EPS0;
                inv_mu[a][b] /= MU0;
            }
        }
        let sigma_e = cell_sigma[ci];
        let has_sigma_e = sigma_e.iter().flatten().any(|&v| v.abs() > 1e-20);

        let h_scale = (measure / perimeter_area.max(1e-30)).max(1e-15);
        h_min = h_min.min(h_scale);
        h_max = h_max.max(h_scale);

        let eps_min_diag = eps_r[0][0].min(eps_r[1][1]).min(eps_r[2][2]).max(1e-3);
        let mu_min_diag = mu_r[0][0].min(mu_r[1][1]).min(mu_r[2][2]).max(1e-3);
        let c_wave = C0 / (eps_min_diag * mu_min_diag).sqrt();
        let poly_factor = ((2 * basis_order + 1) * (basis_order + 1)) as f64 / 2.0;
        let rate_wave = c_wave * poly_factor / h_scale;
        let rate_sigma =
            sigma_e.iter().flatten().map(|v| v.abs()).fold(0.0_f64, f64::max) / (2.0 * EPS0 * eps_min_diag);
        let rate_pml = pml_sigma[0].max(pml_sigma[1]).max(pml_sigma[2]);
        let rate_poles = poles.iter().map(DgtdDispersionPole::characteristic_rate).fold(0.0_f64, f64::max);
        let dt_local = config.cfl / (rate_wave + rate_sigma + rate_pml + 0.5 * rate_poles).max(1e-30);
        dt_min = dt_min.min(dt_local);
        dt_max = dt_max.max(dt_local);

        let nonlinear = if !in_pml {
            domain_nl.get(&dom).and_then(|(nl_spec, euler_deg)| {
                let n0_loc = eps_iso.max(1.0).sqrt();
                let chi2_d_global = nl_spec.chi2.as_ref().and_then(|c2| {
                    let dg = c2.d_global(*euler_deg);
                    let max_d = dg.iter().flatten().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
                    (max_d > 0.0).then_some(dg)
                });
                let node_qpm_signs = if chi2_d_global.is_some() {
                    node_coords.iter().map(|&pt| nl_spec.qpm.as_ref().map_or(1.0, |q| q.poling_sign_at(pt))).collect()
                } else {
                    Vec::new()
                };
                let (chi3_iso_xxxx, chi3_rank4_global, beta_tpa, sigma_fca, k_fcd, tau_c) =
                    if let Some(ref c3) = nl_spec.chi3 {
                        let has_cub = c3.has_cubic_anisotropy(n0_loc);
                        let (a_mt, b_mt) = c3.isotropic_ab(n0_loc);
                        let iso_xxxx = if has_cub { 0.0 } else { (4.0 / 3.0) * (a_mt + 0.5 * b_mt) };
                        let r4 = has_cub.then(|| Box::new(c3.rank4_global(n0_loc, *euler_deg)));
                        (
                            iso_xxxx,
                            r4,
                            c3.beta_tpa.unwrap_or(0.0).max(0.0),
                            c3.sigma_fca.unwrap_or(0.0).max(0.0),
                            c3.k_fcd.unwrap_or(0.0),
                            c3.tau_c.unwrap_or(0.0).max(0.0),
                        )
                    } else {
                        (0.0, None, 0.0, 0.0, 0.0, 0.0)
                    };
                let any_active = chi2_d_global.is_some()
                    || chi3_rank4_global.is_some()
                    || chi3_iso_xxxx.abs() > 0.0
                    || beta_tpa > 0.0;
                any_active.then(|| {
                    Box::new(DgtdCellNonlinear {
                        chi2_d_global,
                        node_qpm_signs,
                        chi3_iso_xxxx,
                        chi3_rank4_global,
                        beta_tpa,
                        sigma_fca,
                        k_fcd,
                        tau_c,
                        n0: n0_loc,
                        omega_0: omega_ref,
                        eps_min_diag,
                    })
                })
            })
        } else {
            None
        };

        cells.push(DgtdCell {
            cell_index: ci,
            cell_type: ct,
            domain_label: dom,
            n_p,
            dof_offset,
            node_coords,
            is_diagonal_mass,
            mass,
            inv_mass,
            dx,
            dy,
            dz,
            eps_r,
            inv_eps,
            mu_r,
            inv_mu,
            sigma_e,
            has_sigma_e,
            y_admittance: y_minus,
            pml_sigma,
            in_pml,
            poles,
            pole_offset,
            nonlinear,
            measure,
            h_scale,
            dt_local,
            lts_tier: 0,
            faces,
        });
    }

    // 8. Partition elements into L power-of-two Local Time-Stepping (LTS) tiers with 1-ring tier smoothing
    let max_tiers = if config.time_stepper == DgtdTimeStepper::Lts { config.max_lts_tiers.clamp(1, 12) } else { 1 };
    let mut raw_tiers = vec![0usize; nc];
    if max_tiers > 1 {
        for ci in 0..nc {
            let ratio = (cells[ci].dt_local / dt_min).max(1.0);
            let tier = (ratio.log2().floor() as usize).min(max_tiers - 1);
            raw_tiers[ci] = tier;
        }
        // Enforce 1-ring 2:1 tier grading: across any interior face, |tier(K^-) - tier(K^+)| <= 1
        for _pass in 0..max_tiers {
            let mut changed = false;
            for ci in 0..nc {
                let mut cap = raw_tiers[ci];
                for f in &cells[ci].faces {
                    if let DgtdFaceKind::Interior { neighbor_cell: nb } = f.kind {
                        cap = cap.min(raw_tiers[nb] + 1);
                    }
                }
                if cap < raw_tiers[ci] {
                    raw_tiers[ci] = cap;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    let n_lts_tiers = raw_tiers.iter().copied().max().unwrap_or(0) + 1;
    let mut tier_cells = vec![Vec::new(); n_lts_tiers];
    for (ci, &t) in raw_tiers.iter().enumerate() {
        cells[ci].lts_tier = t;
        tier_cells[t].push(ci);
    }
    let mut tier_finer_border = vec![Vec::new(); n_lts_tiers];
    for ell in 1..n_lts_tiers {
        let mut border = BTreeSet::new();
        for &ci in &tier_cells[ell] {
            for f in &cells[ci].faces {
                if let DgtdFaceKind::Interior { neighbor_cell: nb } = f.kind {
                    if cells[nb].lts_tier == ell - 1 {
                        border.insert(nb);
                    }
                }
            }
        }
        tier_finer_border[ell] = border.into_iter().collect();
    }

    // 9. Point dipole sources
    let mut dipoles = Vec::new();
    for spec in &problem.excitations {
        if let Excitation::Dipole { position, moment, magnetic } = &spec.excitation {
            let mut found_cell = None;
            for ci in 0..nc {
                if let Some(xi) = locate_point_in_cell(&mixed, ci, *position) {
                    found_cell = Some((ci, xi));
                    break;
                }
            }
            let (ci, xi) = found_cell.unwrap_or_else(|| {
                let mut best_c = 0usize;
                let mut best_d2 = f64::INFINITY;
                for c_idx in 0..nc {
                    for p3 in &cells[c_idx].node_coords {
                        let d2 = dot3(sub3(*position, *p3), sub3(*position, *p3));
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best_c = c_idx;
                        }
                    }
                }
                (best_c, ref_nodes_center(cells[best_c].cell_type))
            });
            dg_shape_into(cells[ci].cell_type, basis_order, &gll, &xi, &mut vals_buf, &mut dxi_buf);
            let n_p = cells[ci].n_p;
            let mut lift_weights = Vec::new();
            if cells[ci].is_diagonal_mass {
                for a in 0..n_p {
                    let w = cells[ci].inv_mass[a] * vals_buf[a];
                    if w.abs() > 1e-15 {
                        lift_weights.push((a, w));
                    }
                }
            } else {
                for a in 0..n_p {
                    let mut acc = 0.0;
                    for b in 0..n_p {
                        acc += cells[ci].inv_mass[a * n_p + b] * vals_buf[b];
                    }
                    if acc.abs() > 1e-15 {
                        lift_weights.push((a, acc));
                    }
                }
            }
            dipoles.push(DgtdDipoleSource {
                cell_index: ci,
                moment: [moment[0].re, moment[1].re, moment[2].re],
                magnetic: *magnetic,
                lift_weights,
            });
        }
    }

    // 10. FieldTimeMonitor probes
    let mut probes = Vec::with_capacity(config.probes.len());
    for &pos in &config.probes {
        let mut found = None;
        for ci in 0..nc {
            if let Some(xi) = locate_point_in_cell(&mixed, ci, pos) {
                found = Some((ci, xi));
                break;
            }
        }
        let (ci, xi) = found.unwrap_or_else(|| {
            let mut best_c = 0usize;
            let mut best_d2 = f64::INFINITY;
            for c_idx in 0..nc {
                for p3 in &cells[c_idx].node_coords {
                    let d2 = dot3(sub3(pos, *p3), sub3(pos, *p3));
                    if d2 < best_d2 {
                        best_d2 = d2;
                        best_c = c_idx;
                    }
                }
            }
            (best_c, ref_nodes_center(cells[best_c].cell_type))
        });
        dg_shape_into(cells[ci].cell_type, basis_order, &gll, &xi, &mut vals_buf, &mut dxi_buf);
        let basis_weights = vals_buf.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 1e-13).collect();
        probes.push(DgtdProbe { position: pos, cell_index: ci, basis_weights });
    }

    // 11. Surface FluxMonitors
    let mut flux_monitors = Vec::with_capacity(config.flux_monitors.len());
    for fm in &config.flux_monitors {
        let mut cell_faces = Vec::new();
        let mut seen = BTreeSet::new();
        for ci in 0..nc {
            for (lf_idx, f) in cells[ci].faces.iter().enumerate() {
                if fm.labels.contains(&f.facet_label) {
                    let key = cell_local_faces[ci][lf_idx].1;
                    if seen.insert(key) {
                        cell_faces.push((ci, lf_idx));
                    }
                }
            }
        }
        flux_monitors.push(DgtdFluxMonitor { name: fm.name.clone(), cell_faces });
    }

    let symmetry_factor = (if problem.symmetry.x_is_pec().is_some() { 2.0 } else { 1.0 })
        * (if problem.symmetry.z_is_pec().is_some() { 2.0 } else { 1.0 });
    let has_carrier_dynamics = cells
        .iter()
        .any(|c| c.nonlinear.as_ref().is_some_and(|nl| nl.beta_tpa > 0.0 && (nl.sigma_fca > 0.0 || nl.k_fcd != 0.0)));

    Ok(DgtdAssembly {
        dim,
        basis_order,
        config: config.clone(),
        mixed,
        cells,
        total_dofs,
        total_pole_states,
        dt_0: dt_min,
        dt_max,
        h_min,
        h_max,
        n_lts_tiers,
        tier_cells,
        tier_finer_border,
        ports: solved_ports.into_iter().map(|sp| sp.port).collect(),
        dipoles,
        probes,
        flux_monitors,
        symmetry_factor,
        has_carrier_dynamics,
    })
}

fn ref_nodes_center(ct: CellType) -> [f64; 3] {
    match ct {
        CellType::Tri => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        CellType::Tet => [0.25, 0.25, 0.25],
        CellType::Prism => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        CellType::Quad | CellType::Hex => [0.0, 0.0, 0.0],
    }
}

// ---------------------------------------------------------------------------
// Stage 2: solve_dgtd (Deterministic parallel RHS, Riemann fluxes, LTS/RK4)
// ---------------------------------------------------------------------------

/// Instantaneous electromagnetic, dispersive, and free-carrier state vector of the DGTD solver.
#[derive(Clone, Debug)]
pub struct DgtdState {
    pub e: Vec<[f64; 3]>,
    pub h: Vec<[f64; 3]>,
    pub j_pole: Vec<[f64; 3]>,
    pub p_pole: Vec<[f64; 3]>,
    /// Free-carrier density $N_c$ ($\text{m}^{-3}$) per nodal DOF when TPA carrier dynamics are active (`empty` otherwise).
    pub n_carrier: Vec<f64>,
}

impl DgtdState {
    pub fn zeros(n_dofs: usize, n_poles: usize) -> Self {
        Self {
            e: vec![[0.0; 3]; n_dofs],
            h: vec![[0.0; 3]; n_dofs],
            j_pole: vec![[0.0; 3]; n_poles],
            p_pole: vec![[0.0; 3]; n_poles],
            n_carrier: Vec::new(),
        }
    }

    pub fn zeros_with_carriers(n_dofs: usize, n_poles: usize, has_carriers: bool) -> Self {
        Self {
            e: vec![[0.0; 3]; n_dofs],
            h: vec![[0.0; 3]; n_dofs],
            j_pole: vec![[0.0; 3]; n_poles],
            p_pole: vec![[0.0; 3]; n_poles],
            n_carrier: if has_carriers { vec![0.0; n_dofs] } else { Vec::new() },
        }
    }
}

#[inline]
fn advance_dgtd_carriers_for_cells(asm: &DgtdAssembly, cells_slice: &[usize], state: &mut DgtdState, dt: f64) {
    if !asm.has_carrier_dynamics || state.n_carrier.is_empty() {
        return;
    }
    const TPA_RMS_SCALE: f64 = 0.816_496_580_927_726; // sqrt(2/3), matching <i_eff(t)^2> = I_opt^2 for harmonic carrier
    for &ci in cells_slice {
        let cell = &asm.cells[ci];
        let Some(ref nl) = cell.nonlinear else { continue };
        if nl.beta_tpa <= 0.0 || (nl.sigma_fca <= 0.0 && nl.k_fcd == 0.0) {
            continue;
        }
        let base = cell.dof_offset;
        let pref = TPA_RMS_SCALE * C0 * EPS0 * nl.n0;
        for a in 0..cell.n_p {
            let ea = state.e[base + a];
            let e2 = dot3(ea, ea);
            let i_eff = pref * e2;
            state.n_carrier[base + a] = femcore::material::drude_carrier_step(
                state.n_carrier[base + a],
                dt,
                nl.omega_0,
                nl.beta_tpa,
                nl.tau_c,
                i_eff,
            );
        }
    }
}

/// Stage 2 output of the DGTD solver (`ADR 9`, `ADR 144`).
#[derive(Clone, Debug)]
pub struct DgtdSolve {
    pub final_state: DgtdState,
    /// Running DFT of nodal electric field `dft_e[freq_idx][dof]` (`len = n_freqs x total_dofs`).
    pub dft_e: Vec<Vec<[C; 3]>>,
    /// Running DFT of nodal magnetic field `dft_h[freq_idx][dof]` (`len = n_freqs x total_dofs`).
    pub dft_h: Vec<Vec<[C; 3]>>,
    /// Running DFT of pulse waveform $\hat{s}(\omega_k)$ at `freqs_hz`.
    pub pulse_dft: Vec<C>,
    /// Running DFT of modal voltage $V_p(\omega_k)$ and inward modal current $I_p(\omega_k)$ per port `port_dft[freq_idx][port_idx]`.
    pub port_dft: Vec<Vec<(C, C)>>,
    /// Running DFT of surface Poynting flux monitors `flux_dft[mon_idx][freq_idx]` -> `(flux_w, e_avg, h_avg, raw_flux_w_s2)`.
    pub flux_dft: Vec<Vec<(f64, C, C, f64)>>,
    pub times_s: Vec<f64>,
    pub energy_j: Vec<f64>,
    pub probe_e: Vec<Vec<[f64; 3]>>,
    pub probe_h: Vec<Vec<[f64; 3]>>,
    pub n_macro_steps: usize,
    pub n_finest_steps: usize,
    pub total_element_updates: u64,
    pub equivalent_uniform_updates: u64,
    pub solve_time_s: f64,
    /// Execution backend actually used (`"cpu"` or `"gpu"`).
    pub backend_used: &'static str,
    /// Optional GPU adapter description when `backend_used == "gpu"`.
    pub gpu_adapter: Option<String>,
    /// Optional fallback reason emitted as a `DiagnosticWarning` (`code: "gpu_fallback_to_cpu"`).
    pub fallback_warning: Option<String>,
}

/// Evaluates the time derivative $(\partial_t \mathbf{E}_K, \partial_t \mathbf{H}_K, \partial_t \mathbf{J}_{p,K}, \partial_t \mathbf{P}_{p,K})$
/// for a single element `cell` at physical time `t_s`, reading neighbor states from `state` (with optional predicted override
/// for coarse/fine LTS neighbor cells).
#[inline]
fn eval_cell_rhs(
    asm: &DgtdAssembly,
    cell: &DgtdCell,
    state: &DgtdState,
    neighbor_override: Option<&dyn Fn(usize, usize) -> ([f64; 3], [f64; 3])>,
    t_s: f64,
    de_out: &mut [[f64; 3]],
    dh_out: &mut [[f64; 3]],
    dj_out: &mut [[f64; 3]],
    dp_out: &mut [[f64; 3]],
) {
    let n_p = cell.n_p;
    let base = cell.dof_offset;
    let e_loc = &state.e[base..base + n_p];
    let h_loc = &state.h[base..base + n_p];
    let alpha = asm.config.flux.alpha();
    let pol = asm.config.polarization;

    // 1. Volume curl: curl_h = D x H, curl_e = D x E
    let mut rhs_d = vec![[0.0; 3]; n_p];
    let mut rhs_b = vec![[0.0; 3]; n_p];
    for a in 0..n_p {
        let row = a * n_p;
        let (mut dxh_y, mut dxh_z) = (0.0, 0.0);
        let (mut dyh_x, mut dyh_z) = (0.0, 0.0);
        let (mut dzh_x, mut dzh_y) = (0.0, 0.0);
        let (mut dxe_y, mut dxe_z) = (0.0, 0.0);
        let (mut dye_x, mut dye_z) = (0.0, 0.0);
        let (mut dze_x, mut dze_y) = (0.0, 0.0);
        for b in 0..n_p {
            let dx_ab = cell.dx[row + b];
            let dy_ab = cell.dy[row + b];
            let dz_ab = cell.dz[row + b];
            let hb = h_loc[b];
            let eb = e_loc[b];
            dxh_y += dx_ab * hb[1];
            dxh_z += dx_ab * hb[2];
            dyh_x += dy_ab * hb[0];
            dyh_z += dy_ab * hb[2];
            dzh_x += dz_ab * hb[0];
            dzh_y += dz_ab * hb[1];

            dxe_y += dx_ab * eb[1];
            dxe_z += dx_ab * eb[2];
            dye_x += dy_ab * eb[0];
            dye_z += dy_ab * eb[2];
            dze_x += dz_ab * eb[0];
            dze_y += dz_ab * eb[1];
        }
        rhs_d[a] = [dyh_z - dzh_y, dzh_x - dxh_z, dxh_y - dyh_x];
        rhs_b[a] = [-(dye_z - dze_y), -(dze_x - dxe_z), -(dxe_y - dye_x)];
    }

    // 2. Surface Riemann fluxes across all local faces
    let y_m = cell.y_admittance;
    let pulse_now = asm.config.pulse.eval(t_s);
    let excited_port_idx = asm.config.excited_port.saturating_sub(1);

    for face in &cell.faces {
        let y_p = face.y_plus;
        let y_sum = (y_m + y_p).max(1e-30);
        let z_m = 1.0 / y_m.max(1e-30);
        let z_p = 1.0 / y_p.max(1e-30);
        let z_sum = (z_m + z_p).max(1e-30);

        for qp in &face.qps {
            let n = qp.n_hat;
            let mut e_m = [0.0; 3];
            let mut h_m = [0.0; 3];
            for &(b, w) in &qp.trace_minus {
                let eb = e_loc[b];
                let hb = h_loc[b];
                e_m[0] += w * eb[0];
                e_m[1] += w * eb[1];
                e_m[2] += w * eb[2];
                h_m[0] += w * hb[0];
                h_m[1] += w * hb[1];
                h_m[2] += w * hb[2];
            }

            let (mut e_p, mut h_p) = match face.kind {
                DgtdFaceKind::Interior { neighbor_cell: nb } => {
                    let nb_base = asm.cells[nb].dof_offset;
                    let mut ep = [0.0; 3];
                    let mut hp = [0.0; 3];
                    if let Some(ov) = neighbor_override {
                        for &(b, w) in &qp.trace_plus {
                            let (eb, hb) = ov(nb, b);
                            ep[0] += w * eb[0];
                            ep[1] += w * eb[1];
                            ep[2] += w * eb[2];
                            hp[0] += w * hb[0];
                            hp[1] += w * hb[1];
                            hp[2] += w * hb[2];
                        }
                    } else {
                        for &(b, w) in &qp.trace_plus {
                            let eb = state.e[nb_base + b];
                            let hb = state.h[nb_base + b];
                            ep[0] += w * eb[0];
                            ep[1] += w * eb[1];
                            ep[2] += w * eb[2];
                            hp[0] += w * hb[0];
                            hp[1] += w * hb[1];
                            hp[2] += w * hb[2];
                        }
                    }
                    (ep, hp)
                }
                DgtdFaceKind::Pec => ([-e_m[0], -e_m[1], -e_m[2]], h_m),
                DgtdFaceKind::Pmc => (e_m, [-h_m[0], -h_m[1], -h_m[2]]),
                DgtdFaceKind::SilverMuller => ([0.0; 3], [0.0; 3]),
                DgtdFaceKind::Port { port_index: pi } => {
                    if pi == excited_port_idx {
                        if let Some((e_mode, h_mode)) = qp.port_mode {
                            (
                                [pulse_now * e_mode[0], pulse_now * e_mode[1], pulse_now * e_mode[2]],
                                [pulse_now * h_mode[0], pulse_now * h_mode[1], pulse_now * h_mode[2]],
                            )
                        } else {
                            ([0.0; 3], [0.0; 3])
                        }
                    } else {
                        ([0.0; 3], [0.0; 3])
                    }
                }
            };

            // TFSF plane-wave surface injection
            if let Some((sign, kh, p_vec, origin, c_bg, y_bg)) = qp.tfsf {
                let delay = dot3(sub3(qp.x, origin), kh) / c_bg;
                let s_inc = asm.config.pulse.eval(t_s - delay);
                let e_inc = [s_inc * p_vec[0], s_inc * p_vec[1], s_inc * p_vec[2]];
                let k_cross_e = cross3(kh, e_inc);
                let h_inc = [y_bg * k_cross_e[0], y_bg * k_cross_e[1], y_bg * k_cross_e[2]];
                for k in 0..3 {
                    e_p[k] += sign * e_inc[k];
                    h_p[k] += sign * h_inc[k];
                }
            }

            let de = sub3(e_p, e_m);
            let dh = sub3(h_p, h_m);
            let n_cross_dh = cross3(n, dh);
            let n_cross_de = cross3(n, de);
            let n_cross_n_cross_de = cross3(n, n_cross_de);
            let n_cross_n_cross_dh = cross3(n, n_cross_dh);

            let (flux_e, flux_h) = if let Some((sigma_2d, u_hat, v_hat)) = face.sheet_2d {
                // Penetrable zero-thickness 2D material sheet (FR-35) Rankine-Hugoniot interface solver:
                // ((Y^- + Y^+) I_2 + sigma_2D) E_t^* = Y^- E_t^- + Y^+ E_t^+ + alpha * (n x (H^+ - H^-))
                let rhs_u = y_m * dot3(e_m, u_hat) + y_p * dot3(e_p, u_hat) + alpha * dot3(n_cross_dh, u_hat);
                let rhs_v = y_m * dot3(e_m, v_hat) + y_p * dot3(e_p, v_hat) + alpha * dot3(n_cross_dh, v_hat);
                let a00 = y_sum + sigma_2d[0][0];
                let a01 = sigma_2d[0][1];
                let a10 = sigma_2d[1][0];
                let a11 = y_sum + sigma_2d[1][1];
                let det_a = (a00 * a11 - a01 * a10).max(1e-30);
                let eu_star = (a11 * rhs_u - a01 * rhs_v) / det_a;
                let ev_star = (-a10 * rhs_u + a00 * rhs_v) / det_a;
                let e_t_star = [
                    eu_star * u_hat[0] + ev_star * v_hat[0],
                    eu_star * u_hat[1] + ev_star * v_hat[1],
                    eu_star * u_hat[2] + ev_star * v_hat[2],
                ];
                // Surface current K = sigma_2D E_t^*
                let ku = sigma_2d[0][0] * eu_star + sigma_2d[0][1] * ev_star;
                let kv = sigma_2d[1][0] * eu_star + sigma_2d[1][1] * ev_star;
                let k_vec =
                    [ku * u_hat[0] + kv * v_hat[0], ku * u_hat[1] + kv * v_hat[1], ku * u_hat[2] + kv * v_hat[2]];
                // Note: n x (n x (E^+ - E^-)) = -(E_t^+ - E_t^-), so dissipative upwind penalty has minus sign on n_cross_n_cross_de
                let c_h = y_p / y_sum;
                let c_e = alpha * (y_m * y_p) / y_sum;
                let c_k = y_m / y_sum;
                let fe = [
                    c_h * n_cross_dh[0] - c_e * n_cross_n_cross_de[0] - c_k * k_vec[0],
                    c_h * n_cross_dh[1] - c_e * n_cross_n_cross_de[1] - c_k * k_vec[1],
                    c_h * n_cross_dh[2] - c_e * n_cross_n_cross_de[2] - c_k * k_vec[2],
                ];
                let de_star = sub3(e_t_star, e_m);
                let n_cross_de_star = cross3(n, de_star);
                let fh = [-n_cross_de_star[0], -n_cross_de_star[1], -n_cross_de_star[2]];
                (fe, fh)
            } else {
                // Note: n x (n x v) = -v_t, so upwind dissipation (+alpha * (E_t^+ - E_t^-)) is (-c_de * n_cross_n_cross_de)
                let c_dh = y_p / y_sum;
                let c_de = alpha * (y_m * y_p) / y_sum;
                let fe = [
                    c_dh * n_cross_dh[0] - c_de * n_cross_n_cross_de[0],
                    c_dh * n_cross_dh[1] - c_de * n_cross_n_cross_de[1],
                    c_dh * n_cross_dh[2] - c_de * n_cross_n_cross_de[2],
                ];
                let c_e = z_p / z_sum;
                let c_h = alpha * (z_m * z_p) / z_sum;
                let fh = [
                    -c_e * n_cross_de[0] - c_h * n_cross_n_cross_dh[0],
                    -c_e * n_cross_de[1] - c_h * n_cross_n_cross_dh[1],
                    -c_e * n_cross_de[2] - c_h * n_cross_n_cross_dh[2],
                ];
                (fe, fh)
            };

            for &(a, lift_w) in &qp.lift_minus {
                rhs_d[a][0] += lift_w * flux_e[0];
                rhs_d[a][1] += lift_w * flux_e[1];
                rhs_d[a][2] += lift_w * flux_e[2];
                rhs_b[a][0] += lift_w * flux_h[0];
                rhs_b[a][1] += lift_w * flux_h[1];
                rhs_b[a][2] += lift_w * flux_h[2];
            }
        }
    }

    // 3. Point dipole sources in this cell
    for dip in &asm.dipoles {
        if dip.cell_index == cell.cell_index {
            for &(a, lw) in &dip.lift_weights {
                let src =
                    [lw * pulse_now * dip.moment[0], lw * pulse_now * dip.moment[1], lw * pulse_now * dip.moment[2]];
                if dip.magnetic {
                    rhs_b[a][0] -= src[0];
                    rhs_b[a][1] -= src[1];
                    rhs_b[a][2] -= src[2];
                } else {
                    rhs_d[a][0] -= src[0];
                    rhs_d[a][1] -= src[1];
                    rhs_d[a][2] -= src[2];
                }
            }
        }
    }

    // 4. Dispersive Drude / Lorentz / Gyrotropic ADE poles, ohmic conductivity, PML damping, and constitutive inversion
    for a in 0..n_p {
        let ea = e_loc[a];
        let ha = h_loc[a];
        if cell.has_sigma_e {
            let j_ohm = mat3_mul_vec(&cell.sigma_e, ea);
            rhs_d[a][0] -= j_ohm[0];
            rhs_d[a][1] -= j_ohm[1];
            rhs_d[a][2] -= j_ohm[2];
        }
        for (p_idx, pole) in cell.poles.iter().enumerate() {
            let slot = p_idx * n_p + a;
            let jp = state.j_pole[cell.pole_offset + slot];
            let pp = state.p_pole[cell.pole_offset + slot];
            rhs_d[a][0] -= jp[0];
            rhs_d[a][1] -= jp[1];
            rhs_d[a][2] -= jp[2];
            match *pole {
                DgtdDispersionPole::Drude { omega_p_rad_s, gamma_rad_s } => {
                    let c_e = EPS0 * omega_p_rad_s * omega_p_rad_s;
                    dj_out[slot] = [
                        c_e * ea[0] - gamma_rad_s * jp[0],
                        c_e * ea[1] - gamma_rad_s * jp[1],
                        c_e * ea[2] - gamma_rad_s * jp[2],
                    ];
                    dp_out[slot] = [0.0; 3];
                }
                DgtdDispersionPole::Lorentz { omega_0_rad_s, gamma_rad_s, delta_eps } => {
                    let w02 = omega_0_rad_s * omega_0_rad_s;
                    let c_e = EPS0 * delta_eps * w02;
                    dp_out[slot] = jp;
                    dj_out[slot] = [
                        c_e * ea[0] - gamma_rad_s * jp[0] - w02 * pp[0],
                        c_e * ea[1] - gamma_rad_s * jp[1] - w02 * pp[1],
                        c_e * ea[2] - gamma_rad_s * jp[2] - w02 * pp[2],
                    ];
                }
                DgtdDispersionPole::GyrotropicDrude { omega_p_rad_s, gamma_rad_s, omega_c_rad_s } => {
                    let c_e = EPS0 * omega_p_rad_s * omega_p_rad_s;
                    let wc_cross_j = cross3(omega_c_rad_s, jp);
                    dj_out[slot] = [
                        c_e * ea[0] - gamma_rad_s * jp[0] - wc_cross_j[0],
                        c_e * ea[1] - gamma_rad_s * jp[1] - wc_cross_j[1],
                        c_e * ea[2] - gamma_rad_s * jp[2] - wc_cross_j[2],
                    ];
                    dp_out[slot] = [0.0; 3];
                }
            }
        }

        let mut de = if let Some(ref nl) = cell.nonlinear {
            let n_c = state.n_carrier.get(base + a).copied().unwrap_or(0.0);
            let e2 = dot3(ea, ea);
            if nl.beta_tpa > 0.0 || (nl.sigma_fca > 0.0 && n_c > 0.0) {
                let sig_nl = (2.0 / 3.0) * EPS0 * EPS0 * C0 * C0 * nl.n0 * nl.n0 * nl.beta_tpa * e2
                    + EPS0 * C0 * nl.n0 * nl.sigma_fca * n_c;
                if sig_nl.is_finite() && sig_nl > 0.0 {
                    rhs_d[a][0] -= sig_nl * ea[0];
                    rhs_d[a][1] -= sig_nl * ea[1];
                    rhs_d[a][2] -= sig_nl * ea[2];
                }
            }
            let mut d_eps = [[0.0; 3]; 3];
            if let Some(ref d_glob) = nl.chi2_d_global {
                let s_qpm = nl.node_qpm_signs.get(a).copied().unwrap_or(1.0);
                let j2 = femcore::material::eval_chi2_tangent_time(d_glob, s_qpm, ea);
                for r in 0..3 {
                    for cc in 0..3 {
                        d_eps[r][cc] += j2[r][cc];
                    }
                }
            }
            if let Some(ref r4) = nl.chi3_rank4_global {
                for r in 0..3 {
                    for l_idx in 0..3 {
                        let mut acc = 0.0;
                        for j in 0..3 {
                            for k in 0..3 {
                                acc += (r4[r][j][k][l_idx] + r4[r][j][l_idx][k] + r4[r][l_idx][j][k]) * ea[j] * ea[k];
                            }
                        }
                        d_eps[r][l_idx] += acc;
                    }
                }
            } else if nl.chi3_iso_xxxx != 0.0 {
                let c3 = nl.chi3_iso_xxxx;
                for r in 0..3 {
                    d_eps[r][r] += c3 * e2;
                    for cc in 0..3 {
                        d_eps[r][cc] += 2.0 * c3 * ea[r] * ea[cc];
                    }
                }
            }
            if nl.k_fcd != 0.0 && n_c > 0.0 {
                let dfcd = 2.0 * nl.n0 * nl.k_fcd * n_c;
                d_eps[0][0] += dfcd;
                d_eps[1][1] += dfcd;
                d_eps[2][2] += dfcd;
            }
            let max_row =
                (0..3).map(|r| d_eps[r][0].abs() + d_eps[r][1].abs() + d_eps[r][2].abs()).fold(0.0_f64, f64::max);
            let cap = 0.45 * nl.eps_min_diag;
            if !max_row.is_finite() {
                d_eps = [[0.0; 3]; 3];
            } else if max_row > cap && max_row > 0.0 {
                let sc = cap / max_row;
                for r in 0..3 {
                    for cc in 0..3 {
                        d_eps[r][cc] *= sc;
                    }
                }
            }
            let mut eps_inst = cell.eps_r;
            for r in 0..3 {
                for cc in 0..3 {
                    eps_inst[r][cc] += d_eps[r][cc];
                }
            }
            match invert_3x3(&eps_inst) {
                Ok(inv_e) => {
                    let v = mat3_mul_vec(&inv_e, rhs_d[a]);
                    [v[0] / EPS0, v[1] / EPS0, v[2] / EPS0]
                }
                Err(_) => mat3_mul_vec(&cell.inv_eps, rhs_d[a]),
            }
        } else {
            mat3_mul_vec(&cell.inv_eps, rhs_d[a])
        };
        let mut dh = mat3_mul_vec(&cell.inv_mu, rhs_b[a]);

        if cell.in_pml {
            // Matched electric and magnetic uniaxial PML damping
            let s_sum = cell.pml_sigma[0] + cell.pml_sigma[1] + cell.pml_sigma[2];
            for k in 0..3 {
                let damp = s_sum - cell.pml_sigma[k];
                de[k] -= damp * ea[k];
                dh[k] -= damp * ha[k];
            }
        }

        if asm.dim == 2 {
            match pol {
                DgtdPolarization::Tm => {
                    de[0] = 0.0;
                    de[1] = 0.0;
                    dh[2] = 0.0;
                }
                DgtdPolarization::Te => {
                    de[2] = 0.0;
                    dh[0] = 0.0;
                    dh[1] = 0.0;
                }
                DgtdPolarization::Full => {}
            }
        }

        de_out[a] = de;
        dh_out[a] = dh;
    }
}

/// Evaluates the RHS for a fixed slice of active cell indices `active_cells` in deterministic fixed chunks (`femcore::par`).
struct TierRhsBuf {
    de: Vec<Vec<[f64; 3]>>,
    dh: Vec<Vec<[f64; 3]>>,
    dj: Vec<Vec<[f64; 3]>>,
    dp: Vec<Vec<[f64; 3]>>,
}

impl TierRhsBuf {
    fn new(asm: &DgtdAssembly, active_cells: &[usize]) -> Self {
        let mut de = Vec::with_capacity(active_cells.len());
        let mut dh = Vec::with_capacity(active_cells.len());
        let mut dj = Vec::with_capacity(active_cells.len());
        let mut dp = Vec::with_capacity(active_cells.len());
        for &ci in active_cells {
            let n_p = asm.cells[ci].n_p;
            let n_poles = asm.cells[ci].poles.len() * n_p;
            de.push(vec![[0.0; 3]; n_p]);
            dh.push(vec![[0.0; 3]; n_p]);
            dj.push(vec![[0.0; 3]; n_poles]);
            dp.push(vec![[0.0; 3]; n_poles]);
        }
        Self { de, dh, dj, dp }
    }
}

fn compute_cells_rhs_parallel(
    asm: &DgtdAssembly,
    active_cells: &[usize],
    state: &DgtdState,
    neighbor_override: Option<&(dyn Fn(usize, usize) -> ([f64; 3], [f64; 3]) + Sync)>,
    t_s: f64,
    out: &mut TierRhsBuf,
) {
    let n = active_cells.len();
    if n == 0 {
        return;
    }
    let chunk = 64usize;
    let nt = femcore::par::threads();
    let mut slots: Vec<(usize, &mut [[f64; 3]], &mut [[f64; 3]], &mut [[f64; 3]], &mut [[f64; 3]])> = active_cells
        .iter()
        .copied()
        .zip(out.de.iter_mut().zip(out.dh.iter_mut()).zip(out.dj.iter_mut().zip(out.dp.iter_mut())))
        .map(|(ci, ((de, dh), (dj, dp)))| {
            (ci, de.as_mut_slice(), dh.as_mut_slice(), dj.as_mut_slice(), dp.as_mut_slice())
        })
        .collect();

    femcore::par::fold_chunks_mut(&mut slots, chunk, nt, |_, slice| {
        for (ci, de, dh, dj, dp) in slice.iter_mut() {
            eval_cell_rhs(
                asm,
                &asm.cells[*ci],
                state,
                neighbor_override.map(|f| f as &dyn Fn(usize, usize) -> ([f64; 3], [f64; 3])),
                t_s,
                de,
                dh,
                dj,
                dp,
            );
        }
        0.0
    });
}

/// Computes total electromagnetic energy $W(t) = \frac{1}{2}\int_\Omega (\varepsilon_0 \mathbf{E}\cdot\boldsymbol{\varepsilon}_r\mathbf{E} + \mu_0 \mathbf{H}\cdot\boldsymbol{\mu}_r\mathbf{H})\,dV$
/// using fixed-order deterministic chunk summation (`ADR 11`).
pub fn compute_total_energy(asm: &DgtdAssembly, state: &DgtdState) -> f64 {
    let nc = asm.cells.len();
    let nt = femcore::par::threads();
    femcore::par::sum_range(nc, femcore::par::CHUNK, nt, |ci| {
        let cell = &asm.cells[ci];
        if cell.in_pml {
            return 0.0;
        }
        let n_p = cell.n_p;
        let base = cell.dof_offset;
        let e_loc = &state.e[base..base + n_p];
        let h_loc = &state.h[base..base + n_p];
        let mut w_cell = 0.0;
        if cell.is_diagonal_mass {
            for a in 0..n_p {
                let ma = cell.mass[a];
                let d_a = mat3_mul_vec(&cell.eps_r, e_loc[a]);
                let b_a = mat3_mul_vec(&cell.mu_r, h_loc[a]);
                w_cell += ma * (EPS0 * dot3(e_loc[a], d_a) + MU0 * dot3(h_loc[a], b_a));
            }
        } else {
            for a in 0..n_p {
                let d_a = mat3_mul_vec(&cell.eps_r, e_loc[a]);
                let b_a = mat3_mul_vec(&cell.mu_r, h_loc[a]);
                for b in 0..n_p {
                    let mab = cell.mass[a * n_p + b];
                    w_cell += mab * (EPS0 * dot3(e_loc[b], d_a) + MU0 * dot3(h_loc[b], b_a));
                }
            }
        }
        0.5 * w_cell
    }) * asm.symmetry_factor
}

/// Recursive power-of-two multi-rate 4th-order Local Time-Stepping (LTS) kernel (`ADR 144`).
///
/// Advances tier `ell` (and all finer tiers `0..ell`) from `t_s` to `t_s + 2^ell * dt_0`.
/// Coarser neighbors in tiers `> ell` provide their 2nd-order RK midpoint-predicted states at any sub-step time `tau`.
#[allow(clippy::too_many_arguments)]
fn advance_lts_tier(
    asm: &DgtdAssembly,
    ell: usize,
    t_s: f64,
    state: &mut DgtdState,
    u0_e: &mut [[f64; 3]],
    u0_h: &mut [[f64; 3]],
    u0_j: &mut [[f64; 3]],
    u0_p: &mut [[f64; 3]],
    k1_e: &mut [[f64; 3]],
    k1_h: &mut [[f64; 3]],
    k2_e: &mut [[f64; 3]],
    k2_h: &mut [[f64; 3]],
    tier_t0: &mut [f64],
    tier_has_k2: &mut [bool],
    rhs_bufs: &mut [TierRhsBuf],
    acc_e: &mut [Vec<Vec<[f64; 3]>>],
    acc_h: &mut [Vec<Vec<[f64; 3]>>],
    acc_j: &mut [Vec<Vec<[f64; 3]>>],
    acc_p: &mut [Vec<Vec<[f64; 3]>>],
    updates_counter: &mut u64,
) {
    let dt = asm.dt_0 * ((1usize << ell) as f64);
    let cells_ell = &asm.tier_cells[ell];
    tier_t0[ell] = t_s;
    tier_has_k2[ell] = false;

    // Save initial state u0 for all cells in tier ell
    for &ci in cells_ell {
        let cell = &asm.cells[ci];
        let r = cell.dof_offset..cell.dof_offset + cell.n_p;
        u0_e[r.clone()].copy_from_slice(&state.e[r.clone()]);
        u0_h[r.clone()].copy_from_slice(&state.h[r]);
        let npoles = cell.poles.len() * cell.n_p;
        if npoles > 0 {
            let pr = cell.pole_offset..cell.pole_offset + npoles;
            u0_j[pr.clone()].copy_from_slice(&state.j_pole[pr.clone()]);
            u0_p[pr.clone()].copy_from_slice(&state.p_pole[pr]);
        }
    }

    // Helper closure evaluating neighbor (ci_nb, loc_b) at time `tau`:
    // If `t_nb > ell`, evaluates the coarse tier's Taylor/midpoint predictor at `tau`.
    let eval_neighbor_at = |tau: f64,
                            st: &DgtdState,
                            u0e: &[[f64; 3]],
                            u0h: &[[f64; 3]],
                            k1e: &[[f64; 3]],
                            k1h: &[[f64; 3]],
                            k2e: &[[f64; 3]],
                            k2h: &[[f64; 3]],
                            t0s: &[f64],
                            has_k2: &[bool],
                            ci_nb: usize,
                            loc_b: usize|
     -> ([f64; 3], [f64; 3]) {
        let t_nb = asm.cells[ci_nb].lts_tier;
        let idx = asm.cells[ci_nb].dof_offset + loc_b;
        if t_nb > ell {
            let dtau = tau - t0s[t_nb];
            if has_k2[t_nb] {
                let dt_nb = asm.dt_0 * ((1usize << t_nb) as f64);
                let c2 = (dtau * dtau) / dt_nb.max(1e-30);
                let c1 = dtau - c2;
                (
                    [
                        u0e[idx][0] + c1 * k1e[idx][0] + c2 * k2e[idx][0],
                        u0e[idx][1] + c1 * k1e[idx][1] + c2 * k2e[idx][1],
                        u0e[idx][2] + c1 * k1e[idx][2] + c2 * k2e[idx][2],
                    ],
                    [
                        u0h[idx][0] + c1 * k1h[idx][0] + c2 * k2h[idx][0],
                        u0h[idx][1] + c1 * k1h[idx][1] + c2 * k2h[idx][1],
                        u0h[idx][2] + c1 * k1h[idx][2] + c2 * k2h[idx][2],
                    ],
                )
            } else {
                (
                    [
                        u0e[idx][0] + dtau * k1e[idx][0],
                        u0e[idx][1] + dtau * k1e[idx][1],
                        u0e[idx][2] + dtau * k1e[idx][2],
                    ],
                    [
                        u0h[idx][0] + dtau * k1h[idx][0],
                        u0h[idx][1] + dtau * k1h[idx][1],
                        u0h[idx][2] + dtau * k1h[idx][2],
                    ],
                )
            }
        } else {
            (st.e[idx], st.h[idx])
        }
    };

    if ell == 0 {
        if cells_ell.is_empty() {
            return;
        }
        // Standard 4-stage RK4 on tier 0
        let stage_weights = [dt / 6.0, dt / 3.0, dt / 3.0, dt / 6.0];
        let stage_steps = [0.5 * dt, 0.5 * dt, dt, 0.0];
        let stage_times = [t_s, t_s + 0.5 * dt, t_s + 0.5 * dt, t_s + dt];

        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let r = cell.dof_offset..cell.dof_offset + cell.n_p;
            acc_e[0][i].copy_from_slice(&u0_e[r.clone()]);
            acc_h[0][i].copy_from_slice(&u0_h[r]);
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pr = cell.pole_offset..cell.pole_offset + npoles;
                acc_j[0][i].copy_from_slice(&u0_j[pr.clone()]);
                acc_p[0][i].copy_from_slice(&u0_p[pr]);
            }
        }

        for s in 0..4 {
            let tau = stage_times[s];
            let ov = |ci_nb: usize, loc_b: usize| {
                eval_neighbor_at(tau, state, u0_e, u0_h, k1_e, k1_h, k2_e, k2_h, tier_t0, tier_has_k2, ci_nb, loc_b)
            };
            compute_cells_rhs_parallel(asm, cells_ell, state, Some(&ov), tau, &mut rhs_bufs[0]);
            *updates_counter += cells_ell.len() as u64;

            let w_acc = stage_weights[s];
            let h_next = stage_steps[s];
            for (i, &ci) in cells_ell.iter().enumerate() {
                let cell = &asm.cells[ci];
                let base = cell.dof_offset;
                let de = &rhs_bufs[0].de[i];
                let dh = &rhs_bufs[0].dh[i];
                for a in 0..cell.n_p {
                    for k in 0..3 {
                        acc_e[0][i][a][k] += w_acc * de[a][k];
                        acc_h[0][i][a][k] += w_acc * dh[a][k];
                        if s < 3 {
                            state.e[base + a][k] = u0_e[base + a][k] + h_next * de[a][k];
                            state.h[base + a][k] = u0_h[base + a][k] + h_next * dh[a][k];
                        } else {
                            state.e[base + a][k] = acc_e[0][i][a][k];
                            state.h[base + a][k] = acc_h[0][i][a][k];
                        }
                    }
                }
                let npoles = cell.poles.len() * cell.n_p;
                if npoles > 0 {
                    let pbase = cell.pole_offset;
                    let dj = &rhs_bufs[0].dj[i];
                    let dp = &rhs_bufs[0].dp[i];
                    for slot in 0..npoles {
                        for k in 0..3 {
                            acc_j[0][i][slot][k] += w_acc * dj[slot][k];
                            acc_p[0][i][slot][k] += w_acc * dp[slot][k];
                            if s < 3 {
                                state.j_pole[pbase + slot][k] = u0_j[pbase + slot][k] + h_next * dj[slot][k];
                                state.p_pole[pbase + slot][k] = u0_p[pbase + slot][k] + h_next * dp[slot][k];
                            } else {
                                state.j_pole[pbase + slot][k] = acc_j[0][i][slot][k];
                                state.p_pole[pbase + slot][k] = acc_p[0][i][slot][k];
                            }
                        }
                    }
                }
            }
        }
        advance_dgtd_carriers_for_cells(asm, cells_ell, state, dt);
        return;
    }

    // Inductive multi-rate case ell > 0:
    // Stage 1 on tier ell at t_s
    if !cells_ell.is_empty() {
        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let r = cell.dof_offset..cell.dof_offset + cell.n_p;
            acc_e[ell][i].copy_from_slice(&u0_e[r.clone()]);
            acc_h[ell][i].copy_from_slice(&u0_h[r]);
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pr = cell.pole_offset..cell.pole_offset + npoles;
                acc_j[ell][i].copy_from_slice(&u0_j[pr.clone()]);
                acc_p[ell][i].copy_from_slice(&u0_p[pr]);
            }
        }

        let ov1 = |ci_nb: usize, loc_b: usize| {
            eval_neighbor_at(t_s, state, u0_e, u0_h, k1_e, k1_h, k2_e, k2_h, tier_t0, tier_has_k2, ci_nb, loc_b)
        };
        compute_cells_rhs_parallel(asm, cells_ell, state, Some(&ov1), t_s, &mut rhs_bufs[ell]);
        *updates_counter += cells_ell.len() as u64;

        let w1 = dt / 6.0;
        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let base = cell.dof_offset;
            let de = &rhs_bufs[ell].de[i];
            let dh = &rhs_bufs[ell].dh[i];
            for a in 0..cell.n_p {
                k1_e[base + a] = de[a];
                k1_h[base + a] = dh[a];
                for k in 0..3 {
                    acc_e[ell][i][a][k] += w1 * de[a][k];
                    acc_h[ell][i][a][k] += w1 * dh[a][k];
                    state.e[base + a][k] = u0_e[base + a][k] + 0.5 * dt * de[a][k];
                    state.h[base + a][k] = u0_h[base + a][k] + 0.5 * dt * dh[a][k];
                }
            }
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pbase = cell.pole_offset;
                let dj = &rhs_bufs[ell].dj[i];
                let dp = &rhs_bufs[ell].dp[i];
                for slot in 0..npoles {
                    for k in 0..3 {
                        acc_j[ell][i][slot][k] += w1 * dj[slot][k];
                        acc_p[ell][i][slot][k] += w1 * dp[slot][k];
                        state.j_pole[pbase + slot][k] = u0_j[pbase + slot][k] + 0.5 * dt * dj[slot][k];
                        state.p_pole[pbase + slot][k] = u0_p[pbase + slot][k] + 0.5 * dt * dp[slot][k];
                    }
                }
            }
        }
    }

    // Advance finer tiers 0..ell-1 over first half-step [t_s, t_s + 0.5 * dt]
    advance_lts_tier(
        asm,
        ell - 1,
        t_s,
        state,
        u0_e,
        u0_h,
        u0_j,
        u0_p,
        k1_e,
        k1_h,
        k2_e,
        k2_h,
        tier_t0,
        tier_has_k2,
        rhs_bufs,
        acc_e,
        acc_h,
        acc_j,
        acc_p,
        updates_counter,
    );

    // Stage 2 & Stage 3 on tier ell at midpoint t_mid = t_s + 0.5 * dt (now tier ell-1 is at its exact midpoint state!)
    let t_mid = t_s + 0.5 * dt;
    if !cells_ell.is_empty() {
        let ov2 = |ci_nb: usize, loc_b: usize| {
            eval_neighbor_at(t_mid, state, u0_e, u0_h, k1_e, k1_h, k2_e, k2_h, tier_t0, tier_has_k2, ci_nb, loc_b)
        };
        compute_cells_rhs_parallel(asm, cells_ell, state, Some(&ov2), t_mid, &mut rhs_bufs[ell]);
        *updates_counter += cells_ell.len() as u64;

        let w2 = dt / 3.0;
        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let base = cell.dof_offset;
            let de = &rhs_bufs[ell].de[i];
            let dh = &rhs_bufs[ell].dh[i];
            for a in 0..cell.n_p {
                k2_e[base + a] = de[a];
                k2_h[base + a] = dh[a];
                for k in 0..3 {
                    acc_e[ell][i][a][k] += w2 * de[a][k];
                    acc_h[ell][i][a][k] += w2 * dh[a][k];
                    state.e[base + a][k] = u0_e[base + a][k] + 0.5 * dt * de[a][k];
                    state.h[base + a][k] = u0_h[base + a][k] + 0.5 * dt * dh[a][k];
                }
            }
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pbase = cell.pole_offset;
                let dj = &rhs_bufs[ell].dj[i];
                let dp = &rhs_bufs[ell].dp[i];
                for slot in 0..npoles {
                    for k in 0..3 {
                        acc_j[ell][i][slot][k] += w2 * dj[slot][k];
                        acc_p[ell][i][slot][k] += w2 * dp[slot][k];
                        state.j_pole[pbase + slot][k] = u0_j[pbase + slot][k] + 0.5 * dt * dj[slot][k];
                        state.p_pole[pbase + slot][k] = u0_p[pbase + slot][k] + 0.5 * dt * dp[slot][k];
                    }
                }
            }
        }
        tier_has_k2[ell] = true;

        // Stage 3 at t_mid
        let ov3 = |ci_nb: usize, loc_b: usize| {
            eval_neighbor_at(t_mid, state, u0_e, u0_h, k1_e, k1_h, k2_e, k2_h, tier_t0, tier_has_k2, ci_nb, loc_b)
        };
        compute_cells_rhs_parallel(asm, cells_ell, state, Some(&ov3), t_mid, &mut rhs_bufs[ell]);
        *updates_counter += cells_ell.len() as u64;

        let w3 = dt / 3.0;
        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let base = cell.dof_offset;
            let de = &rhs_bufs[ell].de[i];
            let dh = &rhs_bufs[ell].dh[i];
            for a in 0..cell.n_p {
                k2_e[base + a] = de[a];
                k2_h[base + a] = dh[a];
                for k in 0..3 {
                    acc_e[ell][i][a][k] += w3 * de[a][k];
                    acc_h[ell][i][a][k] += w3 * dh[a][k];
                    state.e[base + a][k] = u0_e[base + a][k] + dt * de[a][k];
                    state.h[base + a][k] = u0_h[base + a][k] + dt * dh[a][k];
                }
            }
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pbase = cell.pole_offset;
                let dj = &rhs_bufs[ell].dj[i];
                let dp = &rhs_bufs[ell].dp[i];
                for slot in 0..npoles {
                    for k in 0..3 {
                        acc_j[ell][i][slot][k] += w3 * dj[slot][k];
                        acc_p[ell][i][slot][k] += w3 * dp[slot][k];
                        state.j_pole[pbase + slot][k] = u0_j[pbase + slot][k] + dt * dj[slot][k];
                        state.p_pole[pbase + slot][k] = u0_p[pbase + slot][k] + dt * dp[slot][k];
                    }
                }
            }
        }
    }

    // Advance finer tiers 0..ell-1 over second half-step [t_mid, t_s + dt]
    advance_lts_tier(
        asm,
        ell - 1,
        t_mid,
        state,
        u0_e,
        u0_h,
        u0_j,
        u0_p,
        k1_e,
        k1_h,
        k2_e,
        k2_h,
        tier_t0,
        tier_has_k2,
        rhs_bufs,
        acc_e,
        acc_h,
        acc_j,
        acc_p,
        updates_counter,
    );

    // Stage 4 on tier ell at t_end = t_s + dt (now tier ell-1 is at its exact end-of-step state!)
    let t_end = t_s + dt;
    if !cells_ell.is_empty() {
        let ov4 = |ci_nb: usize, loc_b: usize| {
            eval_neighbor_at(t_end, state, u0_e, u0_h, k1_e, k1_h, k2_e, k2_h, tier_t0, tier_has_k2, ci_nb, loc_b)
        };
        compute_cells_rhs_parallel(asm, cells_ell, state, Some(&ov4), t_end, &mut rhs_bufs[ell]);
        *updates_counter += cells_ell.len() as u64;

        let w4 = dt / 6.0;
        for (i, &ci) in cells_ell.iter().enumerate() {
            let cell = &asm.cells[ci];
            let base = cell.dof_offset;
            let de = &rhs_bufs[ell].de[i];
            let dh = &rhs_bufs[ell].dh[i];
            for a in 0..cell.n_p {
                for k in 0..3 {
                    state.e[base + a][k] = acc_e[ell][i][a][k] + w4 * de[a][k];
                    state.h[base + a][k] = acc_h[ell][i][a][k] + w4 * dh[a][k];
                }
            }
            let npoles = cell.poles.len() * cell.n_p;
            if npoles > 0 {
                let pbase = cell.pole_offset;
                let dj = &rhs_bufs[ell].dj[i];
                let dp = &rhs_bufs[ell].dp[i];
                for slot in 0..npoles {
                    for k in 0..3 {
                        state.j_pole[pbase + slot][k] = acc_j[ell][i][slot][k] + w4 * dj[slot][k];
                        state.p_pole[pbase + slot][k] = acc_p[ell][i][slot][k] + w4 * dp[slot][k];
                    }
                }
            }
        }
        advance_dgtd_carriers_for_cells(asm, cells_ell, state, dt);
    }
}

/// Carpenter–Kennedy 1994 5-stage 4th-order Low-Storage Explicit Runge–Kutta (`LSERK4`) coefficients.
const LSERK4_A: [f64; 5] = [
    0.0,
    -567301805773.0 / 1357537059087.0,
    -2404267990393.0 / 2016746695238.0,
    -3550918686646.0 / 2091501179385.0,
    -1275806237668.0 / 842570457699.0,
];
const LSERK4_B: [f64; 5] = [
    1432997174477.0 / 9575080441755.0,
    5161836677717.0 / 13612068292357.0,
    1720146321549.0 / 2090206949498.0,
    3134564353537.0 / 4481467310338.0,
    2277821191437.0 / 14882151754819.0,
];
const LSERK4_C: [f64; 5] = [
    0.0,
    1432997174477.0 / 9575080441755.0,
    2526269341429.0 / 6820363962896.0,
    2006345519317.0 / 3224310063776.0,
    2802321613138.0 / 2924317926251.0,
];

/// Stage 2 of the DGTD solver (`ADR 9`, `ADR 144`): steps the electromagnetic fields over $[0, t_{\text{final}}]$
/// and accumulates running DFT and time-domain monitors.
pub fn solve_dgtd(asm: &DgtdAssembly) -> Result<DgtdSolve, String> {
    #[cfg(not(feature = "gpu"))]
    if asm.config.backend == crate::problem::TdBackend::Gpu {
        return Err(
            "DGTD backend 'gpu' requested, but nyx-ems was compiled without GPU support; recompile with `--features gpu`"
                .into(),
        );
    }

    #[allow(unused_mut)]
    let mut fallback_warning: Option<String> = None;

    #[cfg(feature = "gpu")]
    {
        let try_gpu = match asm.config.backend {
            crate::problem::TdBackend::Cpu => None,
            crate::problem::TdBackend::Gpu => Some(false),
            crate::problem::TdBackend::Auto => {
                if asm.total_dofs >= asm.config.gpu_min_dofs.max(1) {
                    Some(true)
                } else {
                    None
                }
            }
        };
        if let Some(require_discrete) = try_gpu {
            match crate::dgtd_gpu::try_solve_dgtd_gpu(asm, require_discrete)? {
                Ok(gpu_sol) => return Ok(gpu_sol),
                Err(reason) => {
                    fallback_warning = Some(reason);
                }
            }
        }
    }

    let t_start = Instant::now();
    let nc = asm.cells.len();
    let ndofs = asm.total_dofs;
    let npoles = asm.total_pole_states;
    let mut state = DgtdState::zeros_with_carriers(ndofs, npoles, asm.has_carrier_dynamics);

    // Optional initial cavity eigenmode [m, n, p] initialization for closed-cavity energy conservation verification
    if let Some([mx, my, mz]) = asm.config.initial_cavity_mode {
        let (lo, hi) = asm.mixed.bbox();
        let lx = (hi[0] - lo[0]).max(1e-12);
        let ly = (hi[1] - lo[1]).max(1e-12);
        let lz = (hi[2] - lo[2]).max(1e-12);
        let kx = (mx as f64) * PI / lx;
        let ky = (my as f64) * PI / ly;
        let kz = if asm.dim == 3 { (mz as f64) * PI / lz } else { 0.0 };
        for cell in &asm.cells {
            for (a, p) in cell.node_coords.iter().enumerate() {
                let x = p[0] - lo[0];
                let y = p[1] - lo[1];
                let z = p[2] - lo[2];
                if asm.dim == 2 {
                    match asm.config.polarization {
                        DgtdPolarization::Te => {
                            state.e[cell.dof_offset + a] =
                                [(kx * x).cos() * (ky * y).sin(), -(kx * x).sin() * (ky * y).cos(), 0.0];
                        }
                        DgtdPolarization::Tm | DgtdPolarization::Full => {
                            state.e[cell.dof_offset + a] = [0.0, 0.0, (kx * x).sin() * (ky * y).sin()];
                        }
                    }
                } else {
                    // 3D TE_{mx, my, 0} / TM_{mx, my, mz} cavity mode satisfying n x E = 0 on all 6 walls
                    state.e[cell.dof_offset + a] = [0.0, 0.0, (kx * x).sin() * (ky * y).sin() * (kz * z).cos()];
                }
            }
        }
    }

    let max_tier = asm.n_lts_tiers.saturating_sub(1);
    let substeps_per_macro = 1usize << max_tier;
    let dt_macro = asm.dt_0 * (substeps_per_macro as f64);
    let n_macro_steps = ((asm.config.t_final / dt_macro).ceil() as usize).max(1);
    let n_finest_steps = n_macro_steps * substeps_per_macro;

    let n_freqs = asm.config.freqs_hz.len();
    let n_ports = asm.ports.len();
    let n_flux = asm.flux_monitors.len();
    let mut dft_e = vec![vec![[c(0.0, 0.0); 3]; ndofs]; n_freqs];
    let mut dft_h = vec![vec![[c(0.0, 0.0); 3]; ndofs]; n_freqs];
    let mut pulse_dft = vec![c(0.0, 0.0); n_freqs];
    let mut port_dft = vec![vec![(c(0.0, 0.0), c(0.0, 0.0)); n_ports]; n_freqs];
    let mut flux_e_dft: Vec<Vec<Vec<[C; 3]>>> = asm
        .flux_monitors
        .iter()
        .map(|fm| {
            let n_qps: usize = fm.cell_faces.iter().map(|&(ci, lf)| asm.cells[ci].faces[lf].qps.len()).sum();
            vec![vec![[c(0.0, 0.0); 3]; n_qps]; n_freqs]
        })
        .collect();
    let mut flux_h_dft = flux_e_dft.clone();

    let mut times_s = Vec::with_capacity(n_macro_steps + 1);
    let mut energy_j = Vec::with_capacity(n_macro_steps + 1);
    let mut probe_e = vec![Vec::with_capacity(n_macro_steps + 1); asm.probes.len()];
    let mut probe_h = vec![Vec::with_capacity(n_macro_steps + 1); asm.probes.len()];

    let record_snapshot = |t_now: f64,
                           st: &DgtdState,
                           times: &mut Vec<f64>,
                           energies: &mut Vec<f64>,
                           pe: &mut [Vec<[f64; 3]>],
                           ph: &mut [Vec<[f64; 3]>]| {
        times.push(t_now);
        energies.push(compute_total_energy(asm, st));
        for (pi, pr) in asm.probes.iter().enumerate() {
            let base = asm.cells[pr.cell_index].dof_offset;
            let mut ep = [0.0; 3];
            let mut hp = [0.0; 3];
            for &(b, w) in &pr.basis_weights {
                for k in 0..3 {
                    ep[k] += w * st.e[base + b][k];
                    hp[k] += w * st.h[base + b][k];
                }
            }
            pe[pi].push(ep);
            ph[pi].push(hp);
        }
    };

    record_snapshot(0.0, &state, &mut times_s, &mut energy_j, &mut probe_e, &mut probe_h);

    // Pre-collect port faces for fast modal (V_p, I_p) projection
    let mut port_faces: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_ports];
    for ci in 0..nc {
        for (lf_idx, f) in asm.cells[ci].faces.iter().enumerate() {
            if let DgtdFaceKind::Port { port_index: pi } = f.kind {
                port_faces[pi].push((ci, lf_idx));
            }
        }
    }

    let accumulate_dft = |t_now: f64,
                          dt_w: f64,
                          st: &DgtdState,
                          dft_e: &mut [Vec<[C; 3]>],
                          dft_h: &mut [Vec<[C; 3]>],
                          pulse_dft: &mut [C],
                          port_dft: &mut [Vec<(C, C)>],
                          flux_e_dft: &mut [Vec<Vec<[C; 3]>>],
                          flux_h_dft: &mut [Vec<Vec<[C; 3]>>]| {
        let s_now = asm.config.pulse.eval(t_now);
        for (fi, &freq) in asm.config.freqs_hz.iter().enumerate() {
            let angle = -2.0 * PI * freq * t_now;
            let (s_ang, c_ang) = angle.sin_cos();
            let ph = c(c_ang * dt_w, s_ang * dt_w);
            pulse_dft[fi] += ph * s_now;

            let de_f = &mut dft_e[fi];
            let dh_f = &mut dft_h[fi];
            for d in 0..ndofs {
                let ed = st.e[d];
                let hd = st.h[d];
                de_f[d][0] += ph * ed[0];
                de_f[d][1] += ph * ed[1];
                de_f[d][2] += ph * ed[2];
                dh_f[d][0] += ph * hd[0];
                dh_f[d][1] += ph * hd[1];
                dh_f[d][2] += ph * hd[2];
            }

            // Port modal projections V_p(omega) and I_p(omega)
            for pi in 0..n_ports {
                let mut v_acc = 0.0;
                let mut i_acc = 0.0;
                for &(ci, lf_idx) in &port_faces[pi] {
                    let cell = &asm.cells[ci];
                    let base = cell.dof_offset;
                    let imu = 1.0 / ((cell.mu_r[0][0] + cell.mu_r[1][1] + cell.mu_r[2][2]) / 3.0).max(1e-6);
                    let face = &cell.faces[lf_idx];
                    for qp in &face.qps {
                        let Some((e_mode, _)) = qp.port_mode else { continue };
                        let mut e_m = [0.0; 3];
                        let mut h_m = [0.0; 3];
                        for &(b, w) in &qp.trace_minus {
                            let eb = st.e[base + b];
                            let hb = st.h[base + b];
                            for k in 0..3 {
                                e_m[k] += w * eb[k];
                                h_m[k] += w * hb[k];
                            }
                        }
                        // Use the upwind interface state at the port plane for clean modal separation
                        let y_m = cell.y_admittance;
                        let y_p = face.y_plus;
                        let y_sum = (y_m + y_p).max(1e-30);
                        let (e_ext, h_ext) = if pi == asm.config.excited_port.saturating_sub(1) {
                            let hm = qp.port_mode.map_or([0.0; 3], |(_, h)| h);
                            (
                                [s_now * e_mode[0], s_now * e_mode[1], s_now * e_mode[2]],
                                [s_now * hm[0], s_now * hm[1], s_now * hm[2]],
                            )
                        } else {
                            ([0.0; 3], [0.0; 3])
                        };
                        let dh = sub3(h_ext, h_m);
                        let n_cross_dh = cross3(qp.n_hat, dh);
                        let e_star = [
                            (y_m * e_m[0] + y_p * e_ext[0] + n_cross_dh[0]) / y_sum,
                            (y_m * e_m[1] + y_p * e_ext[1] + n_cross_dh[1]) / y_sum,
                            (y_m * e_m[2] + y_p * e_ext[2] + n_cross_dh[2]) / y_sum,
                        ];
                        // n_out x H^* = n_out x H^- - Y^- (E_t^* - E_t^-)
                        let n_cross_h_star = [
                            cross3(qp.n_hat, h_m)[0] - y_m * (e_star[0] - e_m[0]),
                            cross3(qp.n_hat, h_m)[1] - y_m * (e_star[1] - e_m[1]),
                            cross3(qp.n_hat, h_m)[2] - y_m * (e_star[2] - e_m[2]),
                        ];
                        // Inward current projection: for H^in = -Y_0 (n_out x E^in), n_out x H^in = +Y_0 E^in
                        v_acc += qp.w_ds * imu * dot3(e_star, e_mode);
                        i_acc += qp.w_ds * imu * dot3(n_cross_h_star, e_mode);
                    }
                }
                port_dft[fi][pi].0 += ph * v_acc;
                port_dft[fi][pi].1 += ph * i_acc;
            }

            // Surface FluxMonitors
            for (mi, fm) in asm.flux_monitors.iter().enumerate() {
                let mut q_idx = 0usize;
                for &(ci, lf_idx) in &fm.cell_faces {
                    let cell = &asm.cells[ci];
                    let base = cell.dof_offset;
                    let face = &cell.faces[lf_idx];
                    for qp in &face.qps {
                        let mut e_q = [0.0; 3];
                        let mut h_q = [0.0; 3];
                        for &(b, w) in &qp.trace_minus {
                            let eb = st.e[base + b];
                            let hb = st.h[base + b];
                            for k in 0..3 {
                                e_q[k] += w * eb[k];
                                h_q[k] += w * hb[k];
                            }
                        }
                        for k in 0..3 {
                            flux_e_dft[mi][fi][q_idx][k] += ph * e_q[k];
                            flux_h_dft[mi][fi][q_idx][k] += ph * h_q[k];
                        }
                        q_idx += 1;
                    }
                }
            }
        }
    };

    let mut total_element_updates = 0u64;
    let stages_per_step = match asm.config.time_stepper {
        DgtdTimeStepper::Lserk4 => 5u64,
        DgtdTimeStepper::Lts | DgtdTimeStepper::Rk4 => 4u64,
    };
    let equivalent_uniform_updates = (n_finest_steps as u64) * (nc as u64) * stages_per_step;

    match asm.config.time_stepper {
        DgtdTimeStepper::Lts if asm.n_lts_tiers > 1 => {
            let mut u0_e = vec![[0.0; 3]; ndofs];
            let mut u0_h = vec![[0.0; 3]; ndofs];
            let mut u0_j = vec![[0.0; 3]; npoles];
            let mut u0_p = vec![[0.0; 3]; npoles];
            let mut k1_e = vec![[0.0; 3]; ndofs];
            let mut k1_h = vec![[0.0; 3]; ndofs];
            let mut k2_e = vec![[0.0; 3]; ndofs];
            let mut k2_h = vec![[0.0; 3]; ndofs];
            let mut tier_t0 = vec![0.0; asm.n_lts_tiers];
            let mut tier_has_k2 = vec![false; asm.n_lts_tiers];

            let mut rhs_bufs: Vec<TierRhsBuf> =
                (0..asm.n_lts_tiers).map(|ell| TierRhsBuf::new(asm, &asm.tier_cells[ell])).collect();
            let mut acc_e: Vec<Vec<Vec<[f64; 3]>>> = (0..asm.n_lts_tiers)
                .map(|ell| asm.tier_cells[ell].iter().map(|&ci| vec![[0.0; 3]; asm.cells[ci].n_p]).collect())
                .collect();
            let mut acc_h = acc_e.clone();
            let mut acc_j: Vec<Vec<Vec<[f64; 3]>>> = (0..asm.n_lts_tiers)
                .map(|ell| {
                    asm.tier_cells[ell]
                        .iter()
                        .map(|&ci| vec![[0.0; 3]; asm.cells[ci].poles.len() * asm.cells[ci].n_p])
                        .collect()
                })
                .collect();
            let mut acc_p = acc_j.clone();

            for step in 0..n_macro_steps {
                let t_s = (step as f64) * dt_macro;
                if step == 0 {
                    accumulate_dft(
                        0.0,
                        0.5 * dt_macro,
                        &state,
                        &mut dft_e,
                        &mut dft_h,
                        &mut pulse_dft,
                        &mut port_dft,
                        &mut flux_e_dft,
                        &mut flux_h_dft,
                    );
                }
                advance_lts_tier(
                    asm,
                    max_tier,
                    t_s,
                    &mut state,
                    &mut u0_e,
                    &mut u0_h,
                    &mut u0_j,
                    &mut u0_p,
                    &mut k1_e,
                    &mut k1_h,
                    &mut k2_e,
                    &mut k2_h,
                    &mut tier_t0,
                    &mut tier_has_k2,
                    &mut rhs_bufs,
                    &mut acc_e,
                    &mut acc_h,
                    &mut acc_j,
                    &mut acc_p,
                    &mut total_element_updates,
                );
                let t_next = ((step + 1) as f64) * dt_macro;
                let dt_w = if step + 1 == n_macro_steps { 0.5 * dt_macro } else { dt_macro };
                accumulate_dft(
                    t_next,
                    dt_w,
                    &state,
                    &mut dft_e,
                    &mut dft_h,
                    &mut pulse_dft,
                    &mut port_dft,
                    &mut flux_e_dft,
                    &mut flux_h_dft,
                );
                record_snapshot(t_next, &state, &mut times_s, &mut energy_j, &mut probe_e, &mut probe_h);
            }
        }
        DgtdTimeStepper::Lserk4 => {
            let all_cells: Vec<usize> = (0..nc).collect();
            let mut rhs_buf = TierRhsBuf::new(asm, &all_cells);
            let mut res_e = vec![[0.0; 3]; ndofs];
            let mut res_h = vec![[0.0; 3]; ndofs];
            let mut res_j = vec![[0.0; 3]; npoles];
            let mut res_p = vec![[0.0; 3]; npoles];
            let dt = asm.dt_0;

            for step in 0..n_finest_steps {
                let t_s = (step as f64) * dt;
                if step == 0 {
                    accumulate_dft(
                        0.0,
                        0.5 * dt,
                        &state,
                        &mut dft_e,
                        &mut dft_h,
                        &mut pulse_dft,
                        &mut port_dft,
                        &mut flux_e_dft,
                        &mut flux_h_dft,
                    );
                }
                for stage in 0..5 {
                    let tau = t_s + LSERK4_C[stage] * dt;
                    compute_cells_rhs_parallel(asm, &all_cells, &state, None, tau, &mut rhs_buf);
                    total_element_updates += nc as u64;
                    let a_s = LSERK4_A[stage];
                    let b_s = LSERK4_B[stage];
                    for ci in 0..nc {
                        let cell = &asm.cells[ci];
                        let base = cell.dof_offset;
                        let de = &rhs_buf.de[ci];
                        let dh = &rhs_buf.dh[ci];
                        for a in 0..cell.n_p {
                            for k in 0..3 {
                                res_e[base + a][k] = a_s * res_e[base + a][k] + dt * de[a][k];
                                res_h[base + a][k] = a_s * res_h[base + a][k] + dt * dh[a][k];
                                state.e[base + a][k] += b_s * res_e[base + a][k];
                                state.h[base + a][k] += b_s * res_h[base + a][k];
                            }
                        }
                        let np_c = cell.poles.len() * cell.n_p;
                        if np_c > 0 {
                            let pbase = cell.pole_offset;
                            let dj = &rhs_buf.dj[ci];
                            let dp = &rhs_buf.dp[ci];
                            for slot in 0..np_c {
                                for k in 0..3 {
                                    res_j[pbase + slot][k] = a_s * res_j[pbase + slot][k] + dt * dj[slot][k];
                                    res_p[pbase + slot][k] = a_s * res_p[pbase + slot][k] + dt * dp[slot][k];
                                    state.j_pole[pbase + slot][k] += b_s * res_j[pbase + slot][k];
                                    state.p_pole[pbase + slot][k] += b_s * res_p[pbase + slot][k];
                                }
                            }
                        }
                    }
                }
                advance_dgtd_carriers_for_cells(asm, &all_cells, &mut state, dt);
                let t_next = ((step + 1) as f64) * dt;
                let dt_w = if step + 1 == n_finest_steps { 0.5 * dt } else { dt };
                accumulate_dft(
                    t_next,
                    dt_w,
                    &state,
                    &mut dft_e,
                    &mut dft_h,
                    &mut pulse_dft,
                    &mut port_dft,
                    &mut flux_e_dft,
                    &mut flux_h_dft,
                );
                if (step + 1) % substeps_per_macro == 0 || step + 1 == n_finest_steps {
                    record_snapshot(t_next, &state, &mut times_s, &mut energy_j, &mut probe_e, &mut probe_h);
                }
            }
        }
        _ => {
            // Uniform classical 4-stage RK4 (or 1-tier LTS)
            let all_cells: Vec<usize> = (0..nc).collect();
            let mut rhs_buf = TierRhsBuf::new(asm, &all_cells);
            let mut u0_e = vec![[0.0; 3]; ndofs];
            let mut u0_h = vec![[0.0; 3]; ndofs];
            let mut u0_j = vec![[0.0; 3]; npoles];
            let mut u0_p = vec![[0.0; 3]; npoles];
            let mut acc_e = vec![[0.0; 3]; ndofs];
            let mut acc_h = vec![[0.0; 3]; ndofs];
            let mut acc_j = vec![[0.0; 3]; npoles];
            let mut acc_p = vec![[0.0; 3]; npoles];
            let dt = asm.dt_0;
            let stage_weights = [dt / 6.0, dt / 3.0, dt / 3.0, dt / 6.0];
            let stage_steps = [0.5 * dt, 0.5 * dt, dt, 0.0];

            for step in 0..n_finest_steps {
                let t_s = (step as f64) * dt;
                if step == 0 {
                    accumulate_dft(
                        0.0,
                        0.5 * dt,
                        &state,
                        &mut dft_e,
                        &mut dft_h,
                        &mut pulse_dft,
                        &mut port_dft,
                        &mut flux_e_dft,
                        &mut flux_h_dft,
                    );
                }
                u0_e.copy_from_slice(&state.e);
                u0_h.copy_from_slice(&state.h);
                u0_j.copy_from_slice(&state.j_pole);
                u0_p.copy_from_slice(&state.p_pole);
                acc_e.copy_from_slice(&state.e);
                acc_h.copy_from_slice(&state.h);
                acc_j.copy_from_slice(&state.j_pole);
                acc_p.copy_from_slice(&state.p_pole);

                let stage_times = [t_s, t_s + 0.5 * dt, t_s + 0.5 * dt, t_s + dt];
                for s in 0..4 {
                    compute_cells_rhs_parallel(asm, &all_cells, &state, None, stage_times[s], &mut rhs_buf);
                    total_element_updates += nc as u64;
                    let w_acc = stage_weights[s];
                    let h_next = stage_steps[s];
                    for ci in 0..nc {
                        let cell = &asm.cells[ci];
                        let base = cell.dof_offset;
                        let de = &rhs_buf.de[ci];
                        let dh = &rhs_buf.dh[ci];
                        for a in 0..cell.n_p {
                            for k in 0..3 {
                                acc_e[base + a][k] += w_acc * de[a][k];
                                acc_h[base + a][k] += w_acc * dh[a][k];
                                if s < 3 {
                                    state.e[base + a][k] = u0_e[base + a][k] + h_next * de[a][k];
                                    state.h[base + a][k] = u0_h[base + a][k] + h_next * dh[a][k];
                                } else {
                                    state.e[base + a][k] = acc_e[base + a][k];
                                    state.h[base + a][k] = acc_h[base + a][k];
                                }
                            }
                        }
                        let np_c = cell.poles.len() * cell.n_p;
                        if np_c > 0 {
                            let pbase = cell.pole_offset;
                            let dj = &rhs_buf.dj[ci];
                            let dp = &rhs_buf.dp[ci];
                            for slot in 0..np_c {
                                for k in 0..3 {
                                    acc_j[pbase + slot][k] += w_acc * dj[slot][k];
                                    acc_p[pbase + slot][k] += w_acc * dp[slot][k];
                                    if s < 3 {
                                        state.j_pole[pbase + slot][k] = u0_j[pbase + slot][k] + h_next * dj[slot][k];
                                        state.p_pole[pbase + slot][k] = u0_p[pbase + slot][k] + h_next * dp[slot][k];
                                    } else {
                                        state.j_pole[pbase + slot][k] = acc_j[pbase + slot][k];
                                        state.p_pole[pbase + slot][k] = acc_p[pbase + slot][k];
                                    }
                                }
                            }
                        }
                    }
                }
                advance_dgtd_carriers_for_cells(asm, &all_cells, &mut state, dt);
                let t_next = ((step + 1) as f64) * dt;
                let dt_w = if step + 1 == n_finest_steps { 0.5 * dt } else { dt };
                accumulate_dft(
                    t_next,
                    dt_w,
                    &state,
                    &mut dft_e,
                    &mut dft_h,
                    &mut pulse_dft,
                    &mut port_dft,
                    &mut flux_e_dft,
                    &mut flux_h_dft,
                );
                if (step + 1) % substeps_per_macro == 0 || step + 1 == n_finest_steps {
                    record_snapshot(t_next, &state, &mut times_s, &mut energy_j, &mut probe_e, &mut probe_h);
                }
            }
        }
    }

    // Reference peak of the pulse spectrum at carrier f0 (prevents numerical noise blowup at out-of-band nonlinear harmonics like 2*f0)
    let mut s_carrier_dft = c(0.0, 0.0);
    {
        let f_c = asm.config.pulse.f0_hz();
        let dt_m = dt_macro;
        for step in 0..=n_macro_steps {
            let t_s = (step as f64) * dt_m;
            let w_t = if step == 0 || step == n_macro_steps { 0.5 * dt_m } else { dt_m };
            let ang = -2.0 * PI * f_c * t_s;
            let (sa, ca) = ang.sin_cos();
            s_carrier_dft += c(ca * w_t, sa * w_t) * asm.config.pulse.eval(t_s);
        }
    }
    let s_peak = pulse_dft.iter().map(|z| z.norm()).fold(s_carrier_dft.norm(), f64::max);
    let oob_floor = 1e-4 * s_peak;

    // Integrate FluxMonitor DFTs over surface quadrature points
    let mut flux_dft = vec![vec![(0.0, c(0.0, 0.0), c(0.0, 0.0), 0.0); n_freqs]; n_flux];
    for (mi, fm) in asm.flux_monitors.iter().enumerate() {
        for fi in 0..n_freqs {
            let s_hat = pulse_dft[fi];
            let s_norm = if s_hat.norm() >= oob_floor && s_hat.norm() > 1e-30 {
                s_hat
            } else if s_peak > 1e-30 {
                c(s_peak, 0.0)
            } else {
                c(1.0, 0.0)
            };
            let s_norm2 = s_norm.norm_sqr().max(1e-300);
            let mut raw_p_acc = 0.0;
            let mut e_acc = c(0.0, 0.0);
            let mut h_acc = c(0.0, 0.0);
            let mut area_acc = 0.0;
            let mut q_idx = 0usize;
            for &(ci, lf_idx) in &fm.cell_faces {
                let face = &asm.cells[ci].faces[lf_idx];
                for qp in &face.qps {
                    let eq = flux_e_dft[mi][fi][q_idx];
                    let hq = flux_h_dft[mi][fi][q_idx];
                    // 0.5 * Re(E x H*) . n
                    let poynting = [
                        (eq[1] * hq[2].conj() - eq[2] * hq[1].conj()).re,
                        (eq[2] * hq[0].conj() - eq[0] * hq[2].conj()).re,
                        (eq[0] * hq[1].conj() - eq[1] * hq[0].conj()).re,
                    ];
                    raw_p_acc += 0.5 * qp.w_ds * dot3(poynting, qp.n_hat);
                    // Dominant transverse field average
                    let e_dom = if asm.dim == 2 && asm.config.polarization == DgtdPolarization::Tm {
                        eq[2]
                    } else {
                        eq[0] + eq[1] + eq[2]
                    };
                    let h_dom = hq[0] + hq[1] + hq[2];
                    e_acc += (e_dom / s_norm) * qp.w_ds;
                    h_acc += (h_dom / s_norm) * qp.w_ds;
                    area_acc += qp.w_ds;
                    q_idx += 1;
                }
            }
            let inv_a = 1.0 / area_acc.max(1e-30);
            let p_acc = raw_p_acc / s_norm2;
            flux_dft[mi][fi] =
                (p_acc * asm.symmetry_factor, e_acc * inv_a, h_acc * inv_a, raw_p_acc * asm.symmetry_factor);
        }
    }

    Ok(DgtdSolve {
        final_state: state,
        dft_e,
        dft_h,
        pulse_dft,
        port_dft,
        flux_dft,
        times_s,
        energy_j,
        probe_e,
        probe_h,
        n_macro_steps,
        n_finest_steps,
        total_element_updates,
        equivalent_uniform_updates,
        solve_time_s: t_start.elapsed().as_secs_f64(),
        backend_used: "cpu",
        gpu_adapter: None,
        fallback_warning,
    })
}

// ---------------------------------------------------------------------------
// Stage 3: dgtd_fields and run_dgtd
// ---------------------------------------------------------------------------

/// Stage 3 of the DGTD solver (`ADR 9`, `ADR 144`): extracts broadband modal scattering matrices into
/// [`DrivenPortData`], packages LTS and time-monitor outputs into [`DgtdData`], and projects frequency-domain
/// DFT fields onto [`SolveResult`].
pub fn dgtd_fields(mesh: &Mesh, asm: &DgtdAssembly, sol: &DgtdSolve) -> Result<SolveResult, String> {
    let simplex = discretise(mesh).map_err(|e| e.to_string())?;
    let nc_s = simplex.n_cells();
    let nv_s = simplex.n_points();
    let out_order = problem_output_order(asm.basis_order);
    let nd_s = simplex.n_dofs(out_order);
    let nloc_s = simplex.nloc(out_order);
    let bary = local_bary(simplex.dim, out_order);

    let proj = if !asm.mixed.is_pure_simplex() {
        Some(MixedToSimplexProjector::new(&asm.mixed, &simplex).map_err(|e| e.to_string())?)
    } else {
        None
    };
    let gll = Gll::new(asm.basis_order);
    let mut vals_buf = Vec::new();
    let mut dxi_buf = Vec::new();

    let n_freqs = asm.config.freqs_hz.len().max(1);
    let mut solutions = Vec::with_capacity(n_freqs);
    let mut scalars = Vec::new();
    let s_peak_sol = sol.pulse_dft.iter().map(|z| z.norm()).fold(0.0_f64, f64::max);
    let oob_floor_sol = 1e-4 * s_peak_sol;

    for fi in 0..n_freqs {
        let freq_hz = asm.config.freqs_hz.get(fi).copied().unwrap_or_else(|| asm.config.pulse.f0_hz());
        let s_hat = sol.pulse_dft.get(fi).copied().unwrap_or(c(1.0, 0.0));
        let inv_s = if s_hat.norm() >= oob_floor_sol && s_hat.norm() > 1e-30 {
            c(1.0, 0.0) / s_hat
        } else if s_peak_sol > 1e-30 {
            c(1.0 / s_peak_sol, 0.0)
        } else {
            c(1.0, 0.0)
        };

        let mut e_cell =
            [vec![c(0.0, 0.0); nc_s * nloc_s], vec![c(0.0, 0.0); nc_s * nloc_s], vec![c(0.0, 0.0); nc_s * nloc_s]];
        let mut h_cell =
            [vec![c(0.0, 0.0); nc_s * nloc_s], vec![c(0.0, 0.0); nc_s * nloc_s], vec![c(0.0, 0.0); nc_s * nloc_s]];

        for sc in 0..nc_s {
            let orig_ci = proj.as_ref().map_or(sc, |p| p.parent_cell(sc));
            let cell = &asm.cells[orig_ci];
            let base = cell.dof_offset;
            let cv = simplex.cell(sc);
            for (i, l) in bary.iter().enumerate() {
                let mut x_pt = [0.0; 3];
                for (a, &la) in l.iter().enumerate() {
                    let p_v = simplex.point(cv[a] as usize);
                    for k in 0..simplex.dim {
                        x_pt[k] += la * p_v[k];
                    }
                }
                let xi =
                    locate_point_in_cell(&asm.mixed, orig_ci, x_pt).unwrap_or_else(|| ref_nodes_center(cell.cell_type));
                dg_shape_into(cell.cell_type, asm.basis_order, &gll, &xi, &mut vals_buf, &mut dxi_buf);
                let mut eq = [c(0.0, 0.0); 3];
                let mut hq = [c(0.0, 0.0); 3];
                if let (Some(de_f), Some(dh_f)) = (sol.dft_e.get(fi), sol.dft_h.get(fi)) {
                    for (b, &wb) in vals_buf.iter().enumerate() {
                        for k in 0..3 {
                            eq[k] += de_f[base + b][k] * wb * inv_s;
                            hq[k] += dh_f[base + b][k] * wb * inv_s;
                        }
                    }
                } else {
                    for (b, &wb) in vals_buf.iter().enumerate() {
                        for k in 0..3 {
                            eq[k] += c(sol.final_state.e[base + b][k] * wb, 0.0);
                            hq[k] += c(sol.final_state.h[base + b][k] * wb, 0.0);
                        }
                    }
                }
                let slot = sc * nloc_s + i;
                for k in 0..3 {
                    e_cell[k][slot] = eq[k];
                    h_cell[k][slot] = hq[k];
                }
            }
        }

        let split_re_im = |arr: &[C]| -> (Vec<f64>, Vec<f64>) {
            let mut re = Vec::with_capacity(arr.len());
            let mut im = Vec::with_capacity(arr.len());
            for z in arr {
                re.push(z.re);
                im.push(z.im);
            }
            (re, im)
        };
        let (ex_re, ex_im) = split_re_im(&e_cell[0]);
        let (ey_re, ey_im) = split_re_im(&e_cell[1]);
        let (ez_re, ez_im) = split_re_im(&e_cell[2]);
        let (hx_re, hx_im) = split_re_im(&h_cell[0]);
        let (hy_re, hy_im) = split_re_im(&h_cell[1]);
        let (hz_re, hz_im) = split_re_im(&h_cell[2]);

        let mut ex_n_re = node_average(&simplex, &ex_re, 1, out_order);
        let mut ex_n_im = node_average(&simplex, &ex_im, 1, out_order);
        let mut ey_n_re = node_average(&simplex, &ey_re, 1, out_order);
        let mut ey_n_im = node_average(&simplex, &ey_im, 1, out_order);
        let mut ez_n_re = node_average(&simplex, &ez_re, 1, out_order);
        let mut ez_n_im = node_average(&simplex, &ez_im, 1, out_order);
        ex_n_re.truncate(nv_s);
        ex_n_im.truncate(nv_s);
        ey_n_re.truncate(nv_s);
        ey_n_im.truncate(nv_s);
        ez_n_re.truncate(nv_s);
        ez_n_im.truncate(nv_s);
        let mag: Vec<f64> = (0..nv_s)
            .map(|i| {
                (ex_n_re[i].powi(2)
                    + ex_n_im[i].powi(2)
                    + ey_n_re[i].powi(2)
                    + ey_n_im[i].powi(2)
                    + ez_n_re[i].powi(2)
                    + ez_n_im[i].powi(2))
                .sqrt()
            })
            .collect();
        let peak = mag.iter().copied().fold(0.0_f64, f64::max).max(1e-30);
        let pad = |mut v: Vec<f64>| {
            v.resize(nd_s, 0.0);
            v
        };

        let mut cf = CellFields::new(&simplex, out_order);
        cf.push_complex("Ex", "V/m", ex_re, ex_im);
        cf.push_complex("Ey", "V/m", ey_re, ey_im);
        cf.push_complex("Ez", "V/m", ez_re, ez_im);
        cf.push_complex("Hx", "A/m", hx_re, hx_im);
        cf.push_complex("Hy", "A/m", hy_re, hy_im);
        cf.push_complex("Hz", "A/m", hz_re, hz_im);

        solutions.push(Solution {
            label: format!("f = {:.4} GHz", freq_hz / 1e9),
            value: Some(freq_hz),
            eigenvalue: Some([freq_hz, 0.0]),
            fields: vec![
                Field::scalar("|E|", "a.u.", pad(mag.iter().map(|v| v / peak).collect())),
                Field::scalar("Re Ex", "a.u.", pad(ex_n_re.iter().map(|v| v / peak).collect())),
                Field::scalar("Re Ey", "a.u.", pad(ey_n_re.iter().map(|v| v / peak).collect())),
                Field::scalar("Re Ez", "a.u.", pad(ez_n_re.iter().map(|v| v / peak).collect())),
            ],
            cell_fields: cell_fields_of(&cf),
        });
    }

    // Extract broadband DrivenPortData if ports exist
    let p_cnt = asm.ports.len();
    if p_cnt > 0 && !asm.config.freqs_hz.is_empty() {
        let q_exc = asm.config.excited_port.saturating_sub(1).min(p_cnt - 1);
        let mut z0_ohm = Vec::with_capacity(n_freqs);
        let mut s_matrix = Vec::with_capacity(n_freqs);
        let mut z_matrix = Vec::with_capacity(n_freqs);
        let mut y_matrix = Vec::with_capacity(n_freqs);

        for (fi, &freq) in asm.config.freqs_hz.iter().enumerate() {
            let omega = 2.0 * PI * freq;
            let mut z0_f = Vec::with_capacity(p_cnt);
            let mut z_waves = Vec::with_capacity(p_cnt);
            let mut a_waves = vec![c(0.0, 0.0); p_cnt];
            let mut b_waves = vec![c(0.0, 0.0); p_cnt];

            for (pi, port) in asm.ports.iter().enumerate() {
                let (_beta, z_w, z_c) = port.impedance_at(omega);
                z0_f.push([z_c, 0.0]);
                z_waves.push(z_w);
                let (v_p, i_p) = sol.port_dft[fi][pi];
                let sq_z = z_w.max(1e-30).sqrt();
                a_waves[pi] = 0.5 * (v_p / sq_z + sq_z * i_p);
                b_waves[pi] = 0.5 * (v_p / sq_z - sq_z * i_p);
            }

            let a_exc = if a_waves[q_exc].norm() > 1e-30 {
                a_waves[q_exc]
            } else {
                sol.pulse_dft[fi] / z_waves[q_exc].max(1e-30).sqrt()
            };
            let mut s_c = vec![vec![c(0.0, 0.0); p_cnt]; p_cnt];
            for r in 0..p_cnt {
                let s_rq = if a_exc.norm() > 1e-30 { b_waves[r] / a_exc } else { c(0.0, 0.0) };
                s_c[r][q_exc] = s_rq;
                if p_cnt == 2 && r != q_exc {
                    // Populate symmetric reciprocal column for 2-port passive structures when single port excited
                    s_c[q_exc][r] = s_rq;
                    s_c[r][r] = s_c[q_exc][q_exc];
                }
            }

            let z_circ_vec: Vec<f64> = z0_f.iter().map(|z| z[0]).collect();
            let (s_f, z_f, y_f, _) = s_to_z_and_y_checked(&s_c, &z_circ_vec);

            if fi == n_freqs / 2 {
                for r in 0..p_cnt {
                    for col in 0..p_cnt {
                        scalars.push((format!("|S_{}{}| (mid)", r + 1, col + 1), Scalar::Num(s_c[r][col].norm())));
                    }
                }
            }

            z0_ohm.push(z0_f);
            s_matrix.push(s_f);
            z_matrix.push(z_f);
            y_matrix.push(y_f);
        }

        let port_data = DrivenPortData {
            freqs_hz: asm.config.freqs_hz.clone(),
            ports: asm.ports.iter().map(|p| p.name.clone()).collect(),
            z0_ohm,
            s_matrix,
            z_matrix,
            y_matrix,
            hdm_solves: 1,
            rom_dim: 0,
            anchor_freqs_hz: vec![asm.config.pulse.f0_hz()],
            max_residual: 0.0,
        };
        scalars.push(port_data.scalar());
    }

    let lts_work_speedup = (sol.equivalent_uniform_updates as f64) / (sol.total_element_updates.max(1) as f64);
    scalars.push(("dt_min (s)".into(), Scalar::Num(asm.dt_0)));
    scalars.push(("dt_max (s)".into(), Scalar::Num(asm.dt_max)));
    scalars.push(("lts_tiers".into(), Scalar::Num(asm.n_lts_tiers as f64)));
    scalars.push(("lts_work_speedup".into(), Scalar::Num(lts_work_speedup)));
    if let (Some(&e0), Some(&ef)) = (sol.energy_j.first(), sol.energy_j.last()) {
        scalars.push(("energy_initial (J)".into(), Scalar::Num(e0)));
        scalars.push(("energy_final (J)".into(), Scalar::Num(ef)));
    }

    let time_monitors: Vec<DgtdTimeSeries> = asm
        .probes
        .iter()
        .enumerate()
        .map(|(pi, pr)| DgtdTimeSeries {
            position: pr.position,
            cell_index: pr.cell_index,
            e_t: sol.probe_e[pi].clone(),
            h_t: sol.probe_h[pi].clone(),
        })
        .collect();

    let flux_monitors: Vec<DgtdFluxSpectrum> = asm
        .flux_monitors
        .iter()
        .enumerate()
        .map(|(mi, fm)| {
            let flux_w = sol.flux_dft[mi].iter().map(|x| x.0).collect();
            let e_avg = sol.flux_dft[mi].iter().map(|x| [x.1.re, x.1.im]).collect();
            let h_avg = sol.flux_dft[mi].iter().map(|x| [x.2.re, x.2.im]).collect();
            let raw_flux_w_s2 = sol.flux_dft[mi].iter().map(|x| x.3).collect();
            DgtdFluxSpectrum {
                name: fm.name.clone(),
                freqs_hz: asm.config.freqs_hz.clone(),
                raw_flux_w_s2,
                flux_w,
                e_avg,
                h_avg,
            }
        })
        .collect();

    let dgtd_data = DgtdData {
        basis_order: asm.basis_order,
        n_cells: asm.cells.len(),
        total_dofs: asm.total_dofs,
        dt_min_s: asm.dt_0,
        dt_max_s: asm.dt_max,
        h_min_m: asm.h_min,
        h_max_m: asm.h_max,
        n_lts_tiers: asm.n_lts_tiers,
        cells_per_tier: asm.tier_cells.iter().map(Vec::len).collect(),
        n_macro_steps: sol.n_macro_steps,
        n_finest_steps: sol.n_finest_steps,
        total_element_updates: sol.total_element_updates,
        equivalent_uniform_updates: sol.equivalent_uniform_updates,
        lts_work_speedup,
        solve_time_s: sol.solve_time_s,
        times_s: sol.times_s.clone(),
        energy_j: sol.energy_j.clone(),
        pulse_spectrum: sol.pulse_dft.iter().map(|z| [z.re, z.im]).collect(),
        time_monitors,
        flux_monitors,
    };
    scalars.push(dgtd_data.scalar());

    let pml_doms: Vec<i32> =
        asm.cells.iter().filter(|c| c.in_pml).map(|c| c.domain_label).collect::<BTreeSet<_>>().into_iter().collect();
    scalars.push(FieldFrame { lateral: 0, normal: 2, propagation: 1, left_handed: true, outside: pml_doms }.scalar());

    let diagnostics = if sol.fallback_warning.is_some() || sol.backend_used == "gpu" {
        let mut diag = crate::diag::EmsDiag::new(if sol.backend_used == "gpu" { "dgtd_gpu" } else { "dgtd_cpu" });
        diag.iterations = Some(sol.n_finest_steps);
        if let Some(ref msg) = sol.fallback_warning {
            diag.push_warning("gpu_fallback_to_cpu", nyx_physics::DiagnosticSeverity::Warning, msg.clone());
        }
        Some(diag.finish().0)
    } else {
        None
    };

    Ok(SolveResult {
        physics: "ems".into(),
        study: "dgtd".into(),
        order: out_order as u8,
        dim: asm.dim as u8,
        dofs: asm.total_dofs,
        seconds: sol.solve_time_s,
        solutions,
        scalars,
        log: vec![],
        diagnostics,
        error_budget: None,
    })
}

fn problem_output_order(basis_order: usize) -> usize {
    basis_order.clamp(1, 2)
}

/// One-call convenience entry running all three stages ([`assemble_dgtd`], [`solve_dgtd`], [`dgtd_fields`]) of `Study::Dgtd`.
pub fn run_dgtd(mesh: &Mesh, problem: &Problem) -> Result<SolveResult, String> {
    let asm = assemble_dgtd(mesh, problem)?;
    let sol = solve_dgtd(&asm)?;
    dgtd_fields(mesh, &asm, &sol)
}
