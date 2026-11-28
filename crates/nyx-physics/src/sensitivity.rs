// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pure-data descriptors, requests, and reports for discrete adjoint sensitivities and
//! gradient-based shape, material, and topology optimization across physics crates (ADR 142).

use crate::error::Error;
use femcore::linalg::csr::CCsr;
use femcore::linalg::scalar::C64;
use femcore::topology::{HelmholtzFilter, ShapeVelocityField};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

fn default_out_port() -> usize {
    1
}

fn default_filter_radius() -> f64 {
    0.0
}

fn default_simp_penalty() -> f64 {
    3.0
}

fn default_beta() -> f64 {
    1.0
}

fn default_eta() -> f64 {
    0.5
}

fn default_val_min() -> f64 {
    1e-3
}

fn default_val_max() -> f64 {
    1.0
}

fn default_algorithm() -> String {
    "lbfgs".into()
}

fn default_max_iters() -> usize {
    50
}

fn default_gtol() -> f64 {
    1e-6
}

fn default_ftol() -> f64 {
    1e-12
}

fn default_xtol() -> f64 {
    1e-6
}

fn default_move_limit() -> f64 {
    0.2
}

/// Scalar objective differentiated by [`crate::Physics::sensitivity`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SensitivityObjective {
    /// Resonant or modal eigenfrequency `Re(f_m)` (Hz), or normalized squared mismatch
    /// `0.5 * ((Re(f_m) - target_hz) / scale)^2` when `target_hz` is given (`ems` `cavity`, `sms` `eigenfrequency`/`modes`).
    Eigenfrequency {
        #[serde(default)]
        mode: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target_hz: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        scale_hz: Option<f64>,
    },
    /// Waveguide effective index `Re(n_eff,m)` or Bloch index `Re(n_B,m)`, or squared mismatch
    /// `0.5 * (Re(n_eff,m) - target_neff)^2` when `target_neff` is given (`ems` `modes`, `modes_complex`, `bloch`).
    #[serde(alias = "neff")]
    EffectiveIndex {
        #[serde(default)]
        mode: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target_neff: Option<f64>,
    },
    /// Modal quality factor `Q = |Re(f)| / (2 |Im(f)|)`, or loss tangent `1 / Q = 2 |Im(f)| / |Re(f)|`
    /// when `inverse` is `true`.
    #[serde(alias = "q_factor")]
    QualityFactor {
        #[serde(default)]
        mode: usize,
        #[serde(default)]
        inverse: bool,
    },
    /// Scattering power transmission `|S_{out,in}|^2` (or `0.5 * (|S_{out,in}|^2 - target)^2`),
    /// evaluated at `freq_index` or averaged across all frequencies when `freq_index` is `None` (`ems` `driven`).
    #[serde(alias = "s_parameter")]
    Transmission {
        #[serde(default = "default_out_port")]
        out_port: usize,
        #[serde(default)]
        in_port: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        freq_index: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target: Option<f64>,
    },
    /// Multi-port power-splitter mismatch `sum_k 0.5 * (|S_{out_k, in}|^2 - target_k)^2` across `out_ports`
    /// (`targets` defaults to `1 / out_ports.len()` per port if empty).
    SplitterBalance {
        #[serde(default)]
        in_port: usize,
        out_ports: Vec<usize>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        targets: Vec<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        freq_index: Option<usize>,
    },
    /// Electrostatic capacitance `C_{ij}` (F) from terminal sweeps or `2 W_e / (dV)^2` on two-conductor
    /// fixed-potential studies (`dfs` `electrostatic`).
    Capacitance {
        #[serde(default)]
        term_i: usize,
        #[serde(default)]
        term_j: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target_f: Option<f64>,
    },
    /// Structural compliance `f^T u` (`sms` `static`), thermal compliance `f^T T` (`dfs` `heat_stationary`),
    /// or total electrostatic energy `0.5 * V^T K V` (`dfs` `electrostatic`).
    Compliance {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target: Option<f64>,
    },
    /// Domain-averaged temperature, smooth `L_p` peak temperature `(mean(T^p))^(1/p)`, or squared
    /// tracking error `0.5 * mean((T - target_k)^2)` (`dfs` `heat_stationary`).
    #[serde(alias = "mean_temperature")]
    TargetTemperature {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        domain: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        target_k: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        p_norm: Option<f64>,
    },
    /// Integrated field intensity `integral_{domain} |E|^2 dV` in `ems` `driven` (for focusing,
    /// mode-converter, or splitter target-region power).
    FieldPower {
        domain: String,
        #[serde(default)]
        in_port: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        freq_index: Option<usize>,
    },
}

impl SensitivityObjective {
    /// Refuses `NaN` or infinite objective parameters.
    pub fn validate_finite(&self) -> Result<(), Error> {
        let chk = |name: &str, v: Option<f64>| -> Result<(), Error> {
            if let Some(x) = v {
                if !x.is_finite() {
                    return Err(Error::NonFinite(format!("objective parameter '{name}' is not finite ({x})")));
                }
            }
            Ok(())
        };
        match self {
            Self::Eigenfrequency { target_hz, scale_hz, .. } => {
                chk("target_hz", *target_hz)?;
                chk("scale_hz", *scale_hz)?;
                if let Some(s) = scale_hz {
                    if *s <= 0.0 {
                        return Err(Error::invalid("scale_hz must be positive"));
                    }
                }
            }
            Self::EffectiveIndex { target_neff, .. } => chk("target_neff", *target_neff)?,
            Self::QualityFactor { .. } => {}
            Self::Transmission { target, .. } => chk("target", *target)?,
            Self::SplitterBalance { out_ports, targets, .. } => {
                if out_ports.is_empty() {
                    return Err(Error::invalid("splitter_balance requires at least one output port"));
                }
                if !targets.is_empty() && targets.len() != out_ports.len() {
                    return Err(Error::invalid("splitter_balance targets length must match out_ports"));
                }
                for &t in targets {
                    chk("targets[]", Some(t))?;
                }
            }
            Self::Capacitance { target_f, .. } => chk("target_f", *target_f)?,
            Self::Compliance { target } => chk("target", *target)?,
            Self::TargetTemperature { target_k, p_norm, .. } => {
                chk("target_k", *target_k)?;
                chk("p_norm", *p_norm)?;
                if let Some(p) = p_norm {
                    if *p < 1.0 {
                        return Err(Error::invalid("p_norm must be >= 1.0"));
                    }
                }
            }
            Self::FieldPower { domain, .. } => {
                if domain.trim().is_empty() {
                    return Err(Error::invalid("field_power domain must not be empty"));
                }
            }
        }
        Ok(())
    }
}

/// Design parameter differentiated by [`crate::Physics::sensitivity`] and updated by `optimize_run`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DesignParam {
    /// Named scalar CAD parameter in `geom::model::FeatureList::params`.
    ParamVar {
        name: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        min: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max: Option<f64>,
    },
    /// Named scalar material property (`eps_r`, `k`, `E`, `nu`, `rho`, `Q`, `rho_v`, ...) on `region`.
    MaterialScalar {
        region: String,
        property: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        min: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max: Option<f64>,
    },
    /// Physical propagation length `L_c` (metres, `>= 0`) of 0-based EME segment `segment` (`ems` `eme`, ADR 140 / ADR 142).
    #[serde(alias = "segment_length", alias = "cell_length")]
    EmeSegmentLength {
        #[serde(default)]
        segment: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        min: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max: Option<f64>,
    },
    /// Cellwise topology density field `rho in [0, 1]^{N_design}` on `domain`, filtered by the
    /// unstructured Helmholtz PDE filter (`filter_radius`), projected via `heaviside_project(beta, eta)`,
    /// and mapped to physical properties via SIMP / RAMP / plasmon-safe optical interpolation.
    TopologyDensity {
        domain: String,
        #[serde(default = "default_filter_radius")]
        filter_radius: f64,
        #[serde(default = "default_simp_penalty")]
        simp_penalty: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        ramp_q: Option<f64>,
        #[serde(default = "default_beta")]
        beta: f64,
        #[serde(default = "default_eta")]
        eta: f64,
        #[serde(default = "default_val_min")]
        val_min: f64,
        #[serde(default = "default_val_max")]
        val_max: f64,
        #[serde(default)]
        val_min_im: f64,
        #[serde(default)]
        val_max_im: f64,
        #[serde(default)]
        plasmon_damping: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        property: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        initial_density: Option<f64>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        densities: Vec<f64>,
    },
}

impl DesignParam {
    /// Human-readable label for logs and reports.
    pub fn label(&self) -> String {
        match self {
            Self::ParamVar { name, .. } => name.clone(),
            Self::MaterialScalar { region, property, .. } => format!("{region}.{property}"),
            Self::EmeSegmentLength { segment, .. } => format!("eme_length[{segment}]"),
            Self::TopologyDensity { domain, .. } => format!("topology:{domain}"),
        }
    }

    /// Validates parameter bounds and finiteness.
    pub fn validate_finite(&self) -> Result<(), Error> {
        let check_bounds = |label: &str, lo: Option<f64>, hi: Option<f64>| -> Result<(), Error> {
            if let Some(a) = lo {
                if !a.is_finite() {
                    return Err(Error::NonFinite(format!("min bound of '{label}' is not finite")));
                }
            }
            if let Some(b) = hi {
                if !b.is_finite() {
                    return Err(Error::NonFinite(format!("max bound of '{label}' is not finite")));
                }
            }
            if let (Some(a), Some(b)) = (lo, hi) {
                if a > b {
                    return Err(Error::invalid(format!("min bound ({a}) exceeds max bound ({b}) for '{label}'")));
                }
            }
            Ok(())
        };
        match self {
            Self::ParamVar { name, min, max } => {
                if name.trim().is_empty() {
                    return Err(Error::invalid("ParamVar name must not be empty"));
                }
                check_bounds(name, *min, *max)
            }
            Self::MaterialScalar { region, property, min, max } => {
                if region.trim().is_empty() || property.trim().is_empty() {
                    return Err(Error::invalid("MaterialScalar region and property must not be empty"));
                }
                check_bounds(&format!("{region}.{property}"), *min, *max)
            }
            Self::EmeSegmentLength { segment, min, max } => {
                let lbl = format!("eme_length[{segment}]");
                check_bounds(&lbl, *min, *max)?;
                if let Some(lo) = min {
                    if *lo < 0.0 {
                        return Err(Error::invalid(format!("min bound ({lo}) of '{lbl}' must be >= 0")));
                    }
                }
                if let Some(hi) = max {
                    if *hi < 0.0 {
                        return Err(Error::invalid(format!("max bound ({hi}) of '{lbl}' must be >= 0")));
                    }
                }
                Ok(())
            }
            Self::TopologyDensity {
                domain,
                filter_radius,
                simp_penalty,
                ramp_q,
                beta,
                eta,
                val_min,
                val_max,
                val_min_im,
                val_max_im,
                plasmon_damping,
                initial_density,
                densities,
                ..
            } => {
                if domain.trim().is_empty() {
                    return Err(Error::invalid("TopologyDensity domain must not be empty"));
                }
                for (nm, v) in [
                    ("filter_radius", *filter_radius),
                    ("simp_penalty", *simp_penalty),
                    ("beta", *beta),
                    ("eta", *eta),
                    ("val_min", *val_min),
                    ("val_max", *val_max),
                    ("val_min_im", *val_min_im),
                    ("val_max_im", *val_max_im),
                    ("plasmon_damping", *plasmon_damping),
                ] {
                    if !v.is_finite() {
                        return Err(Error::NonFinite(format!("TopologyDensity '{nm}' is not finite ({v})")));
                    }
                }
                if *filter_radius < 0.0 {
                    return Err(Error::invalid("TopologyDensity filter_radius must be >= 0"));
                }
                if *simp_penalty < 1.0 {
                    return Err(Error::invalid("TopologyDensity simp_penalty must be >= 1"));
                }
                if *beta < 0.0 {
                    return Err(Error::invalid("TopologyDensity beta must be >= 0"));
                }
                if !(*eta >= 0.0 && *eta <= 1.0) {
                    return Err(Error::invalid("TopologyDensity eta must lie in [0, 1]"));
                }
                if let Some(q) = ramp_q {
                    if !q.is_finite() || *q < 0.0 {
                        return Err(Error::invalid("TopologyDensity ramp_q must be finite and >= 0"));
                    }
                }
                if let Some(d0) = initial_density {
                    if !d0.is_finite() || !(*d0 >= 0.0 && *d0 <= 1.0) {
                        return Err(Error::invalid("TopologyDensity initial_density must lie in [0, 1]"));
                    }
                }
                if densities.iter().any(|x| !x.is_finite()) {
                    return Err(Error::NonFinite("TopologyDensity densities contain NaN or infinity".into()));
                }
                Ok(())
            }
        }
    }
}

/// Cached modal basis from a previous iterate for Hungarian mode matching (`femcore::linalg::sens::match_modes`).
#[derive(Clone, Debug)]
pub struct TrackedModeBasis {
    pub eigenvalues: Vec<C64>,
    pub eigenvectors: Vec<Vec<C64>>,
    pub mass_matrix: Option<CCsr>,
}

/// Request envelope passed to [`crate::Physics::sensitivity`].
#[derive(Clone, Serialize, Deserialize)]
pub struct SensitivityRequest {
    pub objective: SensitivityObjective,
    pub params: Vec<DesignParam>,
    /// Precomputed harmonic mesh velocities keyed by `ParamVar::name` (populated by the engine or caller).
    #[serde(skip, default)]
    pub shape_velocities: BTreeMap<String, ShapeVelocityField>,
    /// Optional pre-factorised [`HelmholtzFilter`] keyed by index in `params` (reused across optimization steps).
    #[serde(skip, default)]
    pub topology_filters: BTreeMap<usize, Arc<HelmholtzFilter>>,
    /// Optional modal basis from the preceding optimization step so mode crossings are tracked continuously.
    #[serde(skip, default)]
    pub reference_modes: Option<TrackedModeBasis>,
}

impl std::fmt::Debug for SensitivityRequest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SensitivityRequest")
            .field("objective", &self.objective)
            .field("params", &self.params)
            .field("shape_velocities", &self.shape_velocities.keys().collect::<Vec<_>>())
            .field("topology_filters_count", &self.topology_filters.len())
            .field("has_reference_modes", &self.reference_modes.is_some())
            .finish()
    }
}

impl SensitivityRequest {
    pub fn new(objective: SensitivityObjective, params: Vec<DesignParam>) -> Self {
        Self {
            objective,
            params,
            shape_velocities: BTreeMap::new(),
            topology_filters: BTreeMap::new(),
            reference_modes: None,
        }
    }

    pub fn validate_finite(&self) -> Result<(), Error> {
        self.objective.validate_finite()?;
        if self.params.is_empty() {
            return Err(Error::invalid("SensitivityRequest must specify at least one design parameter"));
        }
        for p in &self.params {
            p.validate_finite()?;
        }
        Ok(())
    }
}

/// Evaluated filter + Heaviside state of one [`DesignParam::TopologyDensity`] entry on a mesh.
pub struct EvaluatedTopologyParam {
    pub filter: Arc<HelmholtzFilter>,
    pub rho_raw: Vec<f64>,
    pub rho_tilde: Vec<f64>,
    pub rho_bar: Vec<f64>,
    pub d_rho_bar_d_tilde: Vec<f64>,
}

impl EvaluatedTopologyParam {
    /// Evaluates `rho -> rho_tilde -> rho_bar` on `simplex`, reusing `cached_filter` when supplied.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        simplex: &crate::Discretisation,
        domain_labels: &[i32],
        filter_radius: f64,
        beta: f64,
        eta: f64,
        initial_density: Option<f64>,
        densities: &[f64],
        cached_filter: Option<&Arc<HelmholtzFilter>>,
    ) -> Result<Self, Error> {
        let filter = match cached_filter {
            Some(f) => f.clone(),
            None => {
                Arc::new(HelmholtzFilter::new(simplex, filter_radius, Some(domain_labels)).map_err(Error::from_kernel)?)
            }
        };
        let nc = filter.n_design_cells();
        let rho_raw = if densities.is_empty() {
            vec![initial_density.unwrap_or(0.5); nc]
        } else if densities.len() == nc {
            densities.to_vec()
        } else {
            return Err(Error::invalid(format!(
                "TopologyDensity densities length ({}) does not match active design cell count ({nc})",
                densities.len()
            )));
        };
        let rho_tilde = filter.filter_densities(&rho_raw).map_err(Error::from_kernel)?;
        let (rho_bar, d_rho_bar_d_tilde) = femcore::topology::heaviside_project(&rho_tilde, beta, eta);
        Ok(Self { filter, rho_raw, rho_tilde, rho_bar, d_rho_bar_d_tilde })
    }

    /// Pulls physical-density sensitivities `dJ / d(rho_bar)` back to raw design sensitivities `dJ / d(rho)`.
    pub fn pullback(&self, d_j_d_rho_bar: &[f64]) -> Result<Vec<f64>, Error> {
        if d_j_d_rho_bar.len() != self.rho_bar.len() {
            return Err(Error::invalid("dJ/d(rho_bar) length mismatch in topology pullback"));
        }
        let d_j_d_tilde: Vec<f64> = d_j_d_rho_bar.iter().zip(&self.d_rho_bar_d_tilde).map(|(&g, &dr)| g * dr).collect();
        self.filter.pullback_gradient(&d_j_d_tilde).map_err(Error::from_kernel)
    }

    /// Volume-weighted physical material fraction `sum_e V_e rho_bar_e / sum_e V_e` and its exact
    /// gradient with respect to raw design variables `rho`.
    pub fn volume_fraction_and_grad(&self) -> Result<(f64, Vec<f64>), Error> {
        let v_tot = self.filter.total_volume().max(1e-30);
        let vols = self.filter.cell_volumes();
        let vf: f64 = vols.iter().zip(&self.rho_bar).map(|(&v, &r)| v * r).sum::<f64>() / v_tot;
        let dvf_dbar: Vec<f64> = vols.iter().map(|&v| v / v_tot).collect();
        let grad = self.pullback(&dvf_dbar)?;
        Ok((vf, grad))
    }
}

/// Output returned by [`crate::Physics::sensitivity`].
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SensitivityReport {
    /// Objective value `J(p)`.
    pub objective: f64,
    /// Gradient block per entry of `SensitivityRequest::params`.
    pub gradients: Vec<Vec<f64>>,
    /// Flattened gradient vector across all design variables.
    pub flat_gradient: Vec<f64>,
    /// Human-readable names of `params`.
    pub param_names: Vec<String>,
    /// Diagnostic scalars (`f_hz`, `neff_re`, `Q`, `S_mag_sq`, `capacitance_f`, `compliance`, `volume_fraction`, ...).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub aux_scalars: BTreeMap<String, f64>,
    /// Projected physical densities `rho_bar` per `TopologyDensity` parameter index.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub projected_densities: BTreeMap<usize, Vec<f64>>,
    /// Volume fraction `sum_e V_e rho_bar_e / sum_e V_e` per `TopologyDensity` parameter index.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub volume_fractions: BTreeMap<usize, f64>,
    /// Exact gradient of each `TopologyDensity` volume fraction wrt its raw `rho`, for MMA constraints.
    #[serde(skip, default)]
    pub volume_fraction_gradients: BTreeMap<usize, Vec<f64>>,
    /// Modal basis at the current iterate for continuous mode tracking across optimization steps.
    #[serde(skip, default)]
    pub mode_basis: Option<TrackedModeBasis>,
    /// Stage log lines.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub log: Vec<String>,
}

impl PartialEq for SensitivityReport {
    fn eq(&self, other: &Self) -> bool {
        self.objective == other.objective
            && self.gradients == other.gradients
            && self.flat_gradient == other.flat_gradient
            && self.param_names == other.param_names
            && self.aux_scalars == other.aux_scalars
            && self.projected_densities == other.projected_densities
            && self.volume_fractions == other.volume_fractions
            && self.log == other.log
    }
}

impl SensitivityReport {
    /// Refuses `NaN` or infinite objective or gradient entries.
    pub fn check_finite(&self) -> Result<(), Error> {
        if !self.objective.is_finite() {
            return Err(Error::NonFinite(format!("sensitivity objective is not finite ({})", self.objective)));
        }
        if let Some(pos) = self.flat_gradient.iter().position(|g| !g.is_finite()) {
            return Err(Error::NonFinite(format!("sensitivity gradient contains NaN or infinity at flat index {pos}")));
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Gradient-based optimization specification (`optimize_run`)
// ---------------------------------------------------------------------------

/// Optimization direction (`minimize` by default; `maximize` negates `(J, grad J)` internally).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OptimizeSense {
    #[default]
    Minimize,
    Maximize,
}

/// Inequality constraint for `mma` optimization (`c_i(x) <= 0`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ConstraintSpec {
    /// Bounds the volume-weighted physical density fraction `V(rho_bar)` of `params[param_index]`
    /// to `<= max_fraction` (and optionally `>= min_fraction`).
    VolumeFraction {
        #[serde(default)]
        param_index: usize,
        max_fraction: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        min_fraction: Option<f64>,
    },
}

/// Single Heaviside projection sharpness stage `(beta, max_iters)` in a topology continuation schedule.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ContinuationStage {
    pub beta: f64,
    pub max_iters: usize,
}

/// Specification for gradient-based optimization (`lbfgs`, `ncg`, or `mma`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OptimizeSpec {
    /// `"lbfgs"`, `"ncg"`, or `"mma"`.
    #[serde(default = "default_algorithm")]
    pub algorithm: String,
    #[serde(default)]
    pub sense: OptimizeSense,
    #[serde(default = "default_max_iters")]
    pub max_iters: usize,
    #[serde(default = "default_gtol")]
    pub gtol: f64,
    #[serde(default = "default_ftol")]
    pub ftol: f64,
    #[serde(default = "default_xtol")]
    pub xtol: f64,
    #[serde(default = "default_move_limit")]
    pub move_limit: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub constraints: Vec<ConstraintSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub continuation_schedule: Vec<ContinuationStage>,
}

impl Default for OptimizeSpec {
    fn default() -> Self {
        Self {
            algorithm: default_algorithm(),
            sense: OptimizeSense::Minimize,
            max_iters: default_max_iters(),
            gtol: default_gtol(),
            ftol: default_ftol(),
            xtol: default_xtol(),
            move_limit: default_move_limit(),
            constraints: Vec::new(),
            continuation_schedule: Vec::new(),
        }
    }
}

impl OptimizeSpec {
    pub fn validate_finite(&self) -> Result<(), Error> {
        let alg = self.algorithm.as_str();
        if !matches!(alg, "lbfgs" | "ncg" | "mma") {
            return Err(Error::invalid(format!(
                "unknown optimization algorithm '{alg}'; expected 'lbfgs', 'ncg', or 'mma'"
            )));
        }
        if self.max_iters == 0 {
            return Err(Error::invalid("max_iters must be >= 1"));
        }
        for (nm, v) in [("gtol", self.gtol), ("ftol", self.ftol), ("xtol", self.xtol), ("move_limit", self.move_limit)]
        {
            if !v.is_finite() || v <= 0.0 {
                return Err(Error::invalid(format!("OptimizeSpec '{nm}' must be finite and positive, got {v}")));
            }
        }
        if !self.constraints.is_empty() && alg != "mma" {
            return Err(Error::invalid("nonlinear inequality constraints require algorithm = 'mma'"));
        }
        for c in &self.constraints {
            match c {
                ConstraintSpec::VolumeFraction { max_fraction, min_fraction, .. } => {
                    if !max_fraction.is_finite() || *max_fraction <= 0.0 || *max_fraction > 1.0 {
                        return Err(Error::invalid("VolumeFraction max_fraction must lie in (0, 1]"));
                    }
                    if let Some(lo) = min_fraction {
                        if !lo.is_finite() || *lo < 0.0 || *lo > *max_fraction {
                            return Err(Error::invalid("VolumeFraction min_fraction must lie in [0, max_fraction]"));
                        }
                    }
                }
            }
        }
        for st in &self.continuation_schedule {
            if !st.beta.is_finite() || st.beta < 0.0 || st.max_iters == 0 {
                return Err(Error::invalid("ContinuationStage requires finite beta >= 0 and max_iters >= 1"));
            }
        }
        Ok(())
    }
}
