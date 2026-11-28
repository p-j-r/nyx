// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `dfs`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (names and selectors
//! resolved to mesh numbers, parameters with their defaults), or built directly by code that already has mesh numbers (tests, the
//! legacy window bridge of the engine). The envelope carries scalar material properties only; anisotropic tensors are given in
//! `materials` directly until the material files exist (phase 6).

use crate::Dfs;
use femcore::material::{bcs_surface_impedance, FilmModel, TempCurve, TempTensor};
use femcore::props::{LinearSolver, Materials, PmlDir};
use nyx_physics::decode::{
    domain_names_by_label, materials_by_label, resolve_conditions, tls_interfaces_from_spec, Params,
};
use nyx_physics::{Ctx, Error, MaterialAssign, Physics, ProblemSpec, TlsInterfaces};
use std::collections::BTreeMap;

/// Exact CODATA 2018 Stefan-Boltzmann constant `sigma_SB` in `W / (m^2 K^4)`.
pub const SIGMA_SB: f64 = 5.670_374_419e-8;

/// Coordinate topology of a domain-applied static Infinite Element (`BcKind::InfiniteElement`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IeCoordKind {
    #[default]
    Auto,
    Cartesian,
    Cylindrical,
    Spherical,
}

impl IeCoordKind {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "auto" => Ok(Self::Auto),
            "cartesian" => Ok(Self::Cartesian),
            "cylindrical" => Ok(Self::Cylindrical),
            "spherical" => Ok(Self::Spherical),
            other => Err(format!(
                "unknown infinite_element coord '{other}' (expected 'auto', 'cartesian', 'cylindrical', or 'spherical')"
            )),
        }
    }
}

/// Radial/normal coordinate-stretching profile inside a static Infinite Element domain.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IeProfile {
    /// Rational mapping `[R_in, R_out) -> [R_in, inf)` (locking-free on native `Quad`/`Prism`/`Hex`).
    #[default]
    Rational,
    /// Bounded exponential stretching `chi = exp(beta)` paired with outer Order-1 `far_field` Robin.
    Exponential,
}

impl IeProfile {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim().to_ascii_lowercase().as_str() {
            "" | "rational" => Ok(Self::Rational),
            "exponential" => Ok(Self::Exponential),
            other => Err(format!("unknown infinite_element profile '{other}' (expected 'rational' or 'exponential')")),
        }
    }
}

/// Configuration of a domain-applied static Infinite Element (`BcKind::InfiniteElement`).
#[derive(Clone, Debug, PartialEq)]
pub struct InfiniteElementSpec {
    pub coord: IeCoordKind,
    pub profile: IeProfile,
    pub x: PmlDir,
    pub y: PmlDir,
    pub z: PmlDir,
    pub center: Option<[f64; 3]>,
    pub axis: [f64; 3],
    pub r_in: Option<f64>,
    pub thickness: Option<f64>,
    pub pole_dist: Option<f64>,
    pub decay: f64,
    pub bounds_lo: Option<[f64; 3]>,
    pub bounds_hi: Option<[f64; 3]>,
}

impl Default for InfiniteElementSpec {
    fn default() -> Self {
        Self {
            coord: IeCoordKind::Auto,
            profile: IeProfile::Rational,
            x: PmlDir::Off,
            y: PmlDir::Off,
            z: PmlDir::Off,
            center: None,
            axis: [0.0, 0.0, 1.0],
            r_in: None,
            thickness: None,
            pole_dist: None,
            decay: 30.0,
            bounds_lo: None,
            bounds_hi: None,
        }
    }
}

fn parse_pml_dir(s: &str) -> Result<PmlDir, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "" | "off" => Ok(PmlDir::Off),
        "plus" => Ok(PmlDir::Plus),
        "minus" => Ok(PmlDir::Minus),
        other => Err(format!("unknown infinite_element direction '{other}' (expected 'off', 'plus', or 'minus')")),
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BcKind {
    Potential(f64),
    /// fixed electric potential for steady Ohmic conduction, V
    Voltage(f64),
    Ground,
    /// surface charge density, C/m^2
    SurfaceCharge(f64),
    /// inward normal current density `-n . J`, A/m^2 (steady conduction)
    CurrentDensity(f64),
    /// named conductor terminal for multi-RHS capacitance matrix extraction
    Terminal(String),
    /// spherical 1/r asymptotic Robin open boundary (3D electrostatics)
    Open {
        origin: Option<[f64; 3]>,
    },
    /// asymptotic open/infinite far-field Robin / Bayliss-Turkel boundary condition
    FarField {
        order: usize,
        decay_power: f64,
        center: Option<[f64; 3]>,
    },
    /// domain-applied static Infinite Element (Rational or Bounded Exponential + Order-1 Robin)
    InfiniteElement(InfiniteElementSpec),
    Temperature(f64),
    /// inward heat flux, W/m^2
    Flux(f64),
    Convection {
        h: f64,
        text: f64,
    },
    /// Stefan-Boltzmann radiative cooling/heating `-n . q = emissivity * view_factor * sigma_SB * (T^4 - t_ambient^4)`
    Radiation {
        emissivity: f64,
        t_ambient: f64,
        view_factor: f64,
    },
    /// thin superconducting sheet with sheet kinetic inductance L_k,sq (H/sq)
    LondonSheet {
        l_k_sq: f64,
    },
    /// current terminal excitation (A) and direction vector
    CurrentTerminal {
        name: String,
        current: f64,
        direction: [f64; 3],
    },
    /// superconducting fluxoid / hole loop carrying trapped flux (Wb)
    FluxLoop {
        name: String,
        flux_wb: f64,
    },
}

impl BcKind {
    pub fn name(&self) -> &'static str {
        match self {
            BcKind::Potential(_) => "Electric potential",
            BcKind::Voltage(_) => "Voltage",
            BcKind::Ground => "Ground",
            BcKind::SurfaceCharge(_) => "Surface charge",
            BcKind::CurrentDensity(_) => "Current density",
            BcKind::Terminal(_) => "Terminal",
            BcKind::Open { .. } => "Open boundary",
            BcKind::FarField { .. } => "Far-field open boundary",
            BcKind::InfiniteElement(_) => "Infinite element",
            BcKind::Temperature(_) => "Temperature",
            BcKind::Flux(_) => "Heat flux",
            BcKind::Convection { .. } => "Convective cooling",
            BcKind::Radiation { .. } => "Thermal radiation",
            BcKind::LondonSheet { .. } => "London superconducting sheet",
            BcKind::CurrentTerminal { .. } => "Current terminal",
            BcKind::FluxLoop { .. } => "Fluxoid / hole loop",
        }
    }

    pub fn is_electrostatic(&self) -> bool {
        matches!(
            self,
            BcKind::Potential(_)
                | BcKind::Ground
                | BcKind::SurfaceCharge(_)
                | BcKind::Terminal(_)
                | BcKind::Open { .. }
                | BcKind::FarField { .. }
                | BcKind::InfiniteElement(_)
        )
    }

    pub fn is_conduction(&self) -> bool {
        matches!(self, BcKind::Voltage(_) | BcKind::Potential(_) | BcKind::Ground | BcKind::CurrentDensity(_))
    }

    pub fn is_heat(&self) -> bool {
        matches!(
            self,
            BcKind::Temperature(_)
                | BcKind::Flux(_)
                | BcKind::Convection { .. }
                | BcKind::Radiation { .. }
                | BcKind::FarField { .. }
                | BcKind::InfiniteElement(_)
        )
    }

    pub fn is_magnetostatic(&self) -> bool {
        matches!(
            self,
            BcKind::LondonSheet { .. }
                | BcKind::CurrentTerminal { .. }
                | BcKind::FluxLoop { .. }
                | BcKind::Ground
                | BcKind::FarField { .. }
                | BcKind::InfiniteElement(_)
        )
    }
}

/// A condition on boundary numbers of the mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct Bc {
    pub labels: Vec<i32>,
    pub kind: BcKind,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Study {
    Electrostatic,
    /// Steady Ohmic conduction `-div(sigma(T) grad V) = 0` with Joule heating export `Q_joule = grad V . sigma grad V`
    Conduction,
    HeatStationary,
    /// `t0` is an expression in x, y, z (m) for the initial temperature
    HeatTransient {
        t_end: f64,
        n_out: usize,
        t0: String,
    },
    /// Gauge-stabilized London magnetostatics (`"auto" | "cross_section" | "planar_stream" | "mixed_3d"`)
    Magnetostatic {
        formulation: String,
    },
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Electrostatic => "electrostatic",
            Study::Conduction => "conduction",
            Study::HeatStationary => "heat_stationary",
            Study::HeatTransient { .. } => "heat_transient",
            Study::Magnetostatic { .. } => "magnetostatic",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Problem {
    /// element order, 1 or 2
    pub order: usize,
    pub study: Study,
    pub materials: Materials,
    pub domain_names: BTreeMap<i32, String>,
    pub bcs: Vec<Bc>,
    pub solver: LinearSolver,
    /// iterative solver relative residual tolerance (default `1e-10`)
    pub tol: f64,
    /// Thomas tridiagonal vertical-line AMG smoother on 3D prismatic meshes (default `true`)
    pub amg_line_smooth: bool,
    /// electrostatics formulation (`"volumetric_fem" | "planar_bem" | "auto"`, default `"volumetric_fem"`)
    pub formulation: String,
    /// 3D P2 singular conductor-sheet perimeter edge enrichment (`"none" | "quarter_point"`, default `"none"`)
    pub singular_edges: String,
    /// 1->4 panel subdivision + log-space SPD Richardson extrapolation in `planar_bem` (default `false`)
    pub richardson: bool,
    /// optional `[interfaces.{ma,ms,sa}]` TLS surface-loss parameters resolved from `ctx.library`
    pub tls_interfaces: Option<TlsInterfaces>,
    /// number of goal-oriented multi-pass h-AMR passes (`>= 1`, default `1`)
    pub adapt_passes: usize,
    /// pass-to-pass relative capacitance convergence threshold (`> 0`, default `0.001`)
    pub adapt_rel_tol: f64,
    /// consecutive passes below `adapt_rel_tol` required to stop (`>= 1`, default `1`)
    pub adapt_min_converged_passes: usize,
    /// target mesh complexity / cell growth factor per adaptive pass (`> 1`, default `1.35`)
    pub adapt_growth: f64,
    /// adaptive meshing mode (`"auto" | "sheet_2d_z" | "full_3d"`, default `"auto"`)
    pub adapt_mode: String,
    /// a posteriori cell error indicator (`"flux_jump_goal" | "prager_synge" | "hessian"`, default `"flux_jump_goal"`)
    pub adapt_estimator: String,
    /// Prager-Synge H(div) equilibrated-flux bounds (`"none" | "report" | "midpoint"`, default `"none"`)
    pub flux_bounds: String,
    /// operating frequency in Hz for AC quasi-static dielectric loss (`>= 0`, default `0.0`)
    pub freq: f64,
    /// far-field ambient temperature `T_inf` (K) for `heat_stationary` `far_field` / `infinite_element` BCs (default `293.15`)
    pub farfield_t_inf: f64,
    /// optional 1-based boundary label -> physical metal film thickness `t_m` (m) from `LabelMap`
    pub sheet_thickness_m: BTreeMap<i32, f64>,
    /// optional 1-based boundary label -> substrate over-etch trench depth `t_e` (m) from `LabelMap`
    pub sheet_etch_depth_m: BTreeMap<i32, f64>,
    /// optional per-domain temperature-dependent thermal conductivity curves `(TempTensor, euler_deg)`
    pub k_curves: BTreeMap<i32, (TempTensor, [f64; 3])>,
    /// optional per-domain temperature-dependent specific heat curves `Cp(T)`
    pub cp_curves: BTreeMap<i32, TempCurve>,
    /// optional per-domain temperature-dependent electrical conductivity curves `(TempTensor, euler_deg)`
    pub sigma_curves: BTreeMap<i32, (TempTensor, [f64; 3])>,
    /// optional imported volumetric heat source field (`W/m^3`, nodal or cellwise) from `FieldLink` `"heat_source"`
    pub heat_source: Option<Vec<f64>>,
    /// optional domain-wise conservatively projected per-cell heat source (`nc * nloc(order)`) from cross-mesh `Ctx::transferred_inputs` (`ADR 135`)
    pub heat_source_cell_nodes: Option<Vec<f64>>,
    /// optional imported temperature field (`K`, nodal or cellwise) from `FieldLink` `"temperature"`
    pub temperature_field: Option<Vec<f64>>,
    /// optional domain-filtered per-cell P2 temperature nodes (`nc * nloc(2)`) from cross-mesh `Ctx::transferred_inputs` (`ADR 135`)
    pub temperature_cell_nodes: Option<Vec<f64>>,
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem {
            order,
            study,
            materials: Materials::default(),
            domain_names: BTreeMap::new(),
            bcs: vec![],
            solver: LinearSolver::Auto,
            tol: 1e-10,
            amg_line_smooth: true,
            formulation: "volumetric_fem".into(),
            singular_edges: "none".into(),
            richardson: false,
            tls_interfaces: None,
            adapt_passes: 1,
            adapt_rel_tol: 0.001,
            adapt_min_converged_passes: 1,
            adapt_growth: 1.35,
            adapt_mode: "auto".into(),
            adapt_estimator: "flux_jump_goal".into(),
            flux_bounds: "none".into(),
            freq: 0.0,
            farfield_t_inf: 293.15,
            sheet_thickness_m: BTreeMap::new(),
            sheet_etch_depth_m: BTreeMap::new(),
            k_curves: BTreeMap::new(),
            cp_curves: BTreeMap::new(),
            sigma_curves: BTreeMap::new(),
            heat_source: None,
            heat_source_cell_nodes: None,
            temperature_field: None,
            temperature_cell_nodes: None,
        }
    }

    /// True if any condition is a non-zero Stefan-Boltzmann `Radiation` boundary condition.
    pub fn has_radiation(&self) -> bool {
        self.bcs.iter().any(
            |bc| matches!(bc.kind, BcKind::Radiation { emissivity, view_factor, .. } if emissivity * view_factor > 0.0),
        )
    }

    /// True if stationary or transient heat transfer requires nonlinear SNES Newton iteration.
    pub fn is_nonlinear_heat(&self) -> bool {
        self.has_radiation()
            || !self.k_curves.is_empty()
            || (matches!(self.study, Study::HeatTransient { .. }) && !self.cp_curves.is_empty())
    }

    pub fn with_k_curve(mut self, domain: i32, curve: TempTensor, euler_deg: [f64; 3]) -> Self {
        self.k_curves.insert(domain, (curve, euler_deg));
        self
    }

    pub fn with_cp_curve(mut self, domain: i32, curve: TempCurve) -> Self {
        self.cp_curves.insert(domain, curve);
        self
    }

    pub fn with_sigma_curve(mut self, domain: i32, curve: TempTensor, euler_deg: [f64; 3]) -> Self {
        self.sigma_curves.insert(domain, (curve, euler_deg));
        self
    }

    pub fn with_heat_source(mut self, heat_source: Vec<f64>) -> Self {
        self.heat_source = Some(heat_source);
        self
    }

    pub fn with_temperature_field(mut self, temperature_field: Vec<f64>) -> Self {
        self.temperature_field = Some(temperature_field);
        self
    }

    pub fn with_tls_interfaces(mut self, tls_interfaces: Option<TlsInterfaces>) -> Self {
        self.tls_interfaces = tls_interfaces;
        self
    }

    pub fn with_flux_bounds(mut self, flux_bounds: impl Into<String>) -> Self {
        self.flux_bounds = flux_bounds.into();
        self
    }

    pub fn with_freq(mut self, freq: f64) -> Self {
        self.freq = freq;
        self
    }

    pub fn with_farfield_t_inf(mut self, t_inf: f64) -> Self {
        self.farfield_t_inf = t_inf;
        self
    }

    /// True if any condition is an `InfiniteElement` or `materials.pml` has an active domain.
    pub fn has_infinite_element(&self) -> bool {
        self.bcs.iter().any(|bc| matches!(bc.kind, BcKind::InfiniteElement(_))) || self.materials.has_pml()
    }

    pub fn with_adapt_passes(mut self, adapt_passes: usize) -> Self {
        self.adapt_passes = adapt_passes;
        self
    }

    pub fn with_adapt_estimator(mut self, adapt_estimator: impl Into<String>) -> Self {
        self.adapt_estimator = adapt_estimator.into();
        self
    }
    pub fn with_formulation(mut self, formulation: impl Into<String>) -> Self {
        self.formulation = formulation.into();
        self
    }

    pub fn with_singular_edges(mut self, singular_edges: impl Into<String>) -> Self {
        self.singular_edges = singular_edges.into();
        self
    }

    pub fn with_richardson(mut self, richardson: bool) -> Self {
        self.richardson = richardson;
        self
    }

    pub fn with_amg_line_smooth(mut self, amg_line_smooth: bool) -> Self {
        self.amg_line_smooth = amg_line_smooth;
        self
    }

    pub fn with_tol(mut self, tol: f64) -> Self {
        self.tol = tol;
        self
    }

    /// Every number the user can type must be finite: NaN or infinity does not stop a factorisation, it flows into every field.
    pub fn check_finite(&self) -> Result<(), String> {
        let m = &self.materials;
        for (name, v) in &m.default {
            if !v.is_finite() {
                return Err(format!(
                    "material property '{name}' (default for all domains) is {v}: it must be a finite number"
                ));
            }
        }
        for (dom, props) in &m.domains {
            for (name, v) in props {
                if !v.is_finite() {
                    return Err(format!(
                        "material property '{name}' of domain {dom} is {v}: it must be a finite number"
                    ));
                }
            }
        }
        if !self.tol.is_finite() {
            return Err(format!("solver tolerance 'tol' is {}: it must be a finite number", self.tol));
        }
        if !self.freq.is_finite() {
            return Err(format!("study parameter 'freq' is {}: it must be a finite number", self.freq));
        }
        if !self.adapt_rel_tol.is_finite() {
            return Err(format!(
                "study parameter 'adapt_rel_tol' is {}: it must be a finite number",
                self.adapt_rel_tol
            ));
        }
        if !self.adapt_growth.is_finite() {
            return Err(format!("study parameter 'adapt_growth' is {}: it must be a finite number", self.adapt_growth));
        }
        if !self.farfield_t_inf.is_finite() {
            return Err(format!("far_field ambient temperature t_inf is {}: it must be finite", self.farfield_t_inf));
        }
        for (&dom, &v) in &self.sheet_thickness_m {
            if !v.is_finite() {
                return Err(format!("sheet_thickness_m of boundary {dom} is {v}: it must be a finite number"));
            }
        }
        for (&dom, &v) in &self.sheet_etch_depth_m {
            if !v.is_finite() {
                return Err(format!("sheet_etch_depth_m of boundary {dom} is {v}: it must be a finite number"));
            }
        }
        if let Some(ref cv) = self.heat_source_cell_nodes {
            if cv.iter().any(|v| !v.is_finite()) {
                return Err("projected cell heat_source field contains NaN or infinity".into());
            }
        }
        if let Some(ref tcn) = self.temperature_cell_nodes {
            if tcn.iter().any(|v| !v.is_finite()) {
                return Err("transferred cell temperature field contains NaN or infinity".into());
            }
        }
        for (i, bc) in self.bcs.iter().enumerate() {
            let vals: Vec<f64> = match &bc.kind {
                BcKind::Potential(v)
                | BcKind::Voltage(v)
                | BcKind::SurfaceCharge(v)
                | BcKind::CurrentDensity(v)
                | BcKind::Temperature(v)
                | BcKind::Flux(v) => vec![*v],
                BcKind::Convection { h, text } => vec![*h, *text],
                BcKind::Radiation { emissivity, t_ambient, view_factor } => vec![*emissivity, *t_ambient, *view_factor],
                BcKind::LondonSheet { l_k_sq } => vec![*l_k_sq],
                BcKind::CurrentTerminal { current, direction, .. } => {
                    vec![*current, direction[0], direction[1], direction[2]]
                }
                BcKind::FluxLoop { flux_wb, .. } => vec![*flux_wb],
                BcKind::Open { origin } => origin.map_or_else(Vec::new, |o| vec![o[0], o[1], o[2]]),
                BcKind::FarField { decay_power, center, .. } => {
                    let mut v = vec![*decay_power];
                    if let Some(c) = center {
                        v.extend_from_slice(c);
                    }
                    v
                }
                BcKind::InfiniteElement(ie) => {
                    let mut v = vec![ie.decay, ie.axis[0], ie.axis[1], ie.axis[2]];
                    if let Some(c) = ie.center {
                        v.extend_from_slice(&c);
                    }
                    if let Some(r) = ie.r_in {
                        v.push(r);
                    }
                    if let Some(t) = ie.thickness {
                        v.push(t);
                    }
                    if let Some(p) = ie.pole_dist {
                        v.push(p);
                    }
                    if let Some(lo) = ie.bounds_lo {
                        v.extend_from_slice(&lo);
                    }
                    if let Some(hi) = ie.bounds_hi {
                        v.extend_from_slice(&hi);
                    }
                    v
                }
                BcKind::Ground | BcKind::Terminal(_) => vec![],
            };
            if vals.iter().any(|v| !v.is_finite()) {
                return Err(format!(
                    "boundary condition {} ({} on boundaries {:?}) has a NaN or infinite value",
                    i + 1,
                    bc.kind.name(),
                    bc.labels
                ));
            }
            if let BcKind::FarField { order, decay_power, .. } = &bc.kind {
                if !(1..=2).contains(order) {
                    return Err(format!(
                        "boundary condition {} ({} on boundaries {:?}) has order {order}: it must be 1 or 2",
                        i + 1,
                        bc.kind.name(),
                        bc.labels
                    ));
                }
                if *decay_power <= 0.0 {
                    return Err(format!(
                        "boundary condition {} ({} on boundaries {:?}) has decay_power {decay_power}: it must be positive",
                        i + 1,
                        bc.kind.name(),
                        bc.labels
                    ));
                }
            }
            if let BcKind::InfiniteElement(ie) = &bc.kind {
                if ie.decay < 1.0 {
                    return Err(format!(
                        "condition {} (Infinite element on domains {:?}) has decay {}: it must be >= 1.0",
                        i + 1,
                        bc.labels,
                        ie.decay
                    ));
                }
                if let Some(r) = ie.r_in {
                    if r < 0.0 {
                        return Err(format!(
                            "condition {} (Infinite element on domains {:?}) has r_in {r}: it must be >= 0",
                            i + 1,
                            bc.labels
                        ));
                    }
                }
                if let Some(t) = ie.thickness {
                    if t < 0.0 {
                        return Err(format!(
                            "condition {} (Infinite element on domains {:?}) has thickness {t}: it must be >= 0",
                            i + 1,
                            bc.labels
                        ));
                    }
                }
                if let Some(p) = ie.pole_dist {
                    if p <= 0.0 {
                        return Err(format!(
                            "condition {} (Infinite element on domains {:?}) has pole_dist {p}: it must be > 0",
                            i + 1,
                            bc.labels
                        ));
                    }
                }
            }
            if let BcKind::Radiation { emissivity, t_ambient, view_factor } = &bc.kind {
                if !(0.0..=1.0).contains(emissivity) {
                    return Err(format!(
                        "boundary condition {} (Thermal radiation on boundaries {:?}) has emissivity {emissivity}: it must lie in [0, 1]",
                        i + 1,
                        bc.labels
                    ));
                }
                if !(0.0..=1.0).contains(view_factor) {
                    return Err(format!(
                        "boundary condition {} (Thermal radiation on boundaries {:?}) has view_factor {view_factor}: it must lie in [0, 1]",
                        i + 1,
                        bc.labels
                    ));
                }
                if *t_ambient < 0.0 {
                    return Err(format!(
                        "boundary condition {} (Thermal radiation on boundaries {:?}) has t_ambient {t_ambient}: absolute temperature must be >= 0 K",
                        i + 1,
                        bc.labels
                    ));
                }
            }
        }
        if let Some(ref hs) = self.heat_source {
            if hs.iter().any(|v| !v.is_finite()) {
                return Err("imported 'heat_source' field contains NaN or infinity".into());
            }
        }
        if let Some(ref tf) = self.temperature_field {
            if tf.iter().any(|v| !v.is_finite()) {
                return Err("imported 'temperature' field contains NaN or infinity".into());
            }
        }
        if let Study::HeatTransient { t_end, .. } = &self.study {
            if !t_end.is_finite() {
                return Err("the study's end time is NaN or infinite".into());
            }
        }
        // Tensor data of anisotropic domains, PMLs, and TLS interface definitions: scan outside quoted strings
        let text = format!("{:?} {:?} {:?}", m.aniso, m.pml, self.tls_interfaces);
        if femcore::validate::non_finite_token(&text).is_some() {
            return Err(
                "the problem description contains a NaN or infinite number (an anisotropic material, layer or interface)"
                    .into(),
            );
        }
        Ok(())
    }
}

/// Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and
/// selectors with `ctx.resolver` and refuses combinations that need the study (a heat condition in an electrostatics study).
pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error> {
    let dfs = Dfs;
    let studies = dfs.studies();
    let sd = studies
        .iter()
        .find(|s| s.kind == spec.study.kind)
        .ok_or_else(|| Error::Unsupported(format!("dfs has no study '{}'", spec.study.kind)))?;
    let sp = Params::new(&spec.study.params, &sd.params);
    let mut formulation = "volumetric_fem".to_string();
    let mut singular_edges = "none".to_string();
    let mut richardson = false;
    let mut adapt_passes = 1usize;
    let mut adapt_rel_tol = 0.001_f64;
    let mut adapt_min_converged_passes = 1usize;
    let mut adapt_growth = 1.35_f64;
    let mut adapt_mode = "auto".to_string();
    let mut adapt_estimator = "flux_jump_goal".to_string();
    let mut flux_bounds = "none".to_string();
    let mut freq = 0.0_f64;
    let mut farfield_t_inf = 293.15_f64;
    let study = match sd.kind.as_str() {
        "electrostatic" => {
            formulation = sp.text("formulation")?;
            singular_edges = sp.text("singular_edges")?;
            richardson = sp.boolean("richardson")?;
            let ap = sp.integer("adapt_passes")?;
            if ap < 1 {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_passes' must be >= 1"));
            }
            adapt_passes = ap as usize;
            adapt_rel_tol = sp.number("adapt_rel_tol")?;
            if !(adapt_rel_tol.is_finite() && adapt_rel_tol > 0.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_rel_tol' must be > 0"));
            }
            let amc = sp.integer("adapt_min_converged_passes")?;
            if amc < 1 {
                return Err(Error::invalid(
                    "study 'electrostatic': parameter 'adapt_min_converged_passes' must be >= 1",
                ));
            }
            adapt_min_converged_passes = amc as usize;
            adapt_growth = sp.number("adapt_growth")?;
            if !(adapt_growth.is_finite() && adapt_growth > 1.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_growth' must be > 1"));
            }
            adapt_mode = sp.text("adapt_mode")?;
            adapt_estimator = sp.text("adapt_estimator")?;
            flux_bounds = sp.text("flux_bounds")?;
            freq = sp.number("freq")?;
            if !(freq.is_finite() && freq >= 0.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'freq' must be >= 0"));
            }
            Study::Electrostatic
        }
        "conduction" => Study::Conduction,
        "heat_stationary" => Study::HeatStationary,
        "heat_transient" => Study::HeatTransient {
            t_end: sp.number("t_end")?,
            n_out: sp.integer("n_out")?.max(0) as usize,
            t0: sp.text("t0")?,
        },
        "magnetostatic" => Study::Magnetostatic { formulation: sp.text("formulation")? },
        other => return Err(Error::Unsupported(format!("dfs has no study '{other}'"))),
    };
    let kinds = dfs.bc_kinds();
    let mut bcs = vec![];
    for r in resolve_conditions(&spec.bcs, &kinds, ctx)? {
        let d = kinds
            .iter()
            .find(|k| k.kind == r.kind)
            .ok_or_else(|| Error::Unsupported(format!("dfs has no condition '{}'", r.kind)))?;
        let p = Params::new(&r.params, &d.params);
        let kind = match r.kind.as_str() {
            "potential" => BcKind::Potential(p.number("value")?),
            "voltage" => BcKind::Voltage(p.number("value")?),
            "ground" => BcKind::Ground,
            "surface_charge" => BcKind::SurfaceCharge(p.number("density")?),
            "current_density" => BcKind::CurrentDensity(p.number("value")?),
            "terminal" => {
                let name = p.text("name")?.trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("terminal condition requires a non-empty 'name'"));
                }
                BcKind::Terminal(name)
            }
            "open" => {
                let origin = if p.has("origin") {
                    let v = p.vector("origin")?;
                    if v.len() != 3 || v.iter().any(|x| !x.is_finite()) {
                        return Err(Error::invalid("open boundary 'origin' must be 3 finite numbers [x, y, z]"));
                    }
                    Some([v[0], v[1], v[2]])
                } else {
                    None
                };
                BcKind::Open { origin }
            }
            "far_field" => {
                let order_i = p.integer("order")?;
                if !(1..=2).contains(&order_i) {
                    return Err(Error::invalid(format!("far_field condition 'order' must be 1 or 2, got {order_i}")));
                }
                let decay_power = p.number("decay_power")?;
                if !decay_power.is_finite() || decay_power <= 0.0 {
                    return Err(Error::invalid(format!(
                        "far_field condition 'decay_power' must be positive, got {decay_power}"
                    )));
                }
                let center = if p.has("center") {
                    let cv = p.vector("center")?;
                    if cv.len() != 3 || cv.iter().any(|x| !x.is_finite()) {
                        return Err(Error::invalid(
                            "far_field condition 'center' must be 3 finite numbers [cx, cy, cz]",
                        ));
                    }
                    Some([cv[0], cv[1], cv[2]])
                } else {
                    None
                };
                if p.has("t_inf") {
                    let t_val = p.number("t_inf")?;
                    if !t_val.is_finite() {
                        return Err(Error::invalid("far_field condition 't_inf' must be a finite number"));
                    }
                    farfield_t_inf = t_val;
                }
                BcKind::FarField { order: order_i as usize, decay_power, center }
            }
            "infinite_element" => {
                let coord = IeCoordKind::parse(&p.text("coord")?).map_err(Error::invalid)?;
                let profile = IeProfile::parse(&p.text("profile")?).map_err(Error::invalid)?;
                let x = parse_pml_dir(&p.text("x")?).map_err(Error::invalid)?;
                let y = parse_pml_dir(&p.text("y")?).map_err(Error::invalid)?;
                let z = parse_pml_dir(&p.text("z")?).map_err(Error::invalid)?;
                let center = if p.has("center") {
                    let cv = p.vector("center")?;
                    if cv.len() != 3 || cv.iter().any(|v| !v.is_finite()) {
                        return Err(Error::invalid("infinite_element 'center' must be 3 finite numbers [cx, cy, cz]"));
                    }
                    Some([cv[0], cv[1], cv[2]])
                } else {
                    None
                };
                let ax_v = p.vector("axis")?;
                if ax_v.len() != 3 || ax_v.iter().any(|v| !v.is_finite()) {
                    return Err(Error::invalid("infinite_element 'axis' must be 3 finite numbers [ax, ay, az]"));
                }
                let axis = [ax_v[0], ax_v[1], ax_v[2]];
                if axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] <= 1e-30 {
                    return Err(Error::invalid("infinite_element 'axis' must be a non-zero vector"));
                }
                let r_in = p.optional_number("r_in")?.filter(|&r| r > 0.0);
                let thickness = p.optional_number("thickness")?.filter(|&t| t > 0.0);
                let pole_dist = p.optional_number("pole_dist")?.filter(|&pd| pd > 0.0);
                let decay = p.number("decay")?;
                if !decay.is_finite() || decay < 1.0 {
                    return Err(Error::invalid(format!("infinite_element 'decay' (chi) must be >= 1.0, got {decay}")));
                }
                let bounds_lo = if p.has("bounds_lo") {
                    let v = p.vector("bounds_lo")?;
                    if v.len() != 3 || v.iter().any(|x| !x.is_finite()) {
                        return Err(Error::invalid("infinite_element 'bounds_lo' must be 3 finite numbers [x, y, z]"));
                    }
                    Some([v[0], v[1], v[2]])
                } else {
                    None
                };
                let bounds_hi = if p.has("bounds_hi") {
                    let v = p.vector("bounds_hi")?;
                    if v.len() != 3 || v.iter().any(|x| !x.is_finite()) {
                        return Err(Error::invalid("infinite_element 'bounds_hi' must be 3 finite numbers [x, y, z]"));
                    }
                    Some([v[0], v[1], v[2]])
                } else {
                    None
                };
                BcKind::InfiniteElement(InfiniteElementSpec {
                    coord,
                    profile,
                    x,
                    y,
                    z,
                    center,
                    axis,
                    r_in,
                    thickness,
                    pole_dist,
                    decay,
                    bounds_lo,
                    bounds_hi,
                })
            }
            "temperature" => BcKind::Temperature(p.number("value")?),
            "flux" => BcKind::Flux(p.number("value")?),
            "convection" => BcKind::Convection { h: p.number("h")?, text: p.number("t_ext")? },
            "radiation" => {
                let emissivity = p.number("emissivity")?;
                let t_ambient = p.number("t_ambient")?;
                let view_factor = p.number("view_factor")?;
                if !emissivity.is_finite() || !(0.0..=1.0).contains(&emissivity) {
                    return Err(Error::invalid(format!(
                        "radiation condition 'emissivity' must lie in [0, 1], got {emissivity}"
                    )));
                }
                if !view_factor.is_finite() || !(0.0..=1.0).contains(&view_factor) {
                    return Err(Error::invalid(format!(
                        "radiation condition 'view_factor' must lie in [0, 1], got {view_factor}"
                    )));
                }
                if !t_ambient.is_finite() || t_ambient < 0.0 {
                    return Err(Error::invalid(format!(
                        "radiation condition 't_ambient' must be >= 0 K, got {t_ambient}"
                    )));
                }
                BcKind::Radiation { emissivity, t_ambient, view_factor }
            }
            "london_sheet" => {
                let l_k_sq_given = p.number("L_k_sq")?;
                let mat_name = p.text("material")?;
                let mat_trim = mat_name.trim();
                let l_k_sq = if !mat_trim.is_empty() {
                    let thickness = p.number("thickness")?;
                    let temperature = p.number("temperature")?;
                    let fm_str = p.text("film_model")?;
                    let film_model = FilmModel::parse(&fm_str)
                        .ok_or_else(|| Error::invalid(format!("unknown london_sheet film_model '{fm_str}'")))?;
                    let lib = ctx.library.as_ref().ok_or_else(|| {
                        Error::Unsupported(format!(
                            "london_sheet material '{mat_trim}' needs a material library, and this run has none"
                        ))
                    })?;
                    let mat = lib.require(mat_trim).map_err(|e| Error::invalid(e.to_string()))?;
                    let sc = mat.superconducting.as_ref().ok_or_else(|| {
                        Error::invalid(format!("material '{mat_trim}' has no [superconducting] block"))
                    })?;
                    let eval =
                        bcs_surface_impedance(sc, 1.0e9, temperature, thickness, film_model).map_err(Error::invalid)?;
                    eval.l_k_sq
                } else {
                    if l_k_sq_given < 0.0 {
                        return Err(Error::invalid("london_sheet L_k_sq must be non-negative"));
                    }
                    l_k_sq_given
                };
                BcKind::LondonSheet { l_k_sq }
            }
            "current_terminal" => {
                let name = p.text("name")?.trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("current_terminal condition requires a non-empty 'name'"));
                }
                let current = p.number("current")?;
                let dir_v = p.vector("direction")?;
                let direction = [
                    dir_v.first().copied().unwrap_or(1.0),
                    dir_v.get(1).copied().unwrap_or(0.0),
                    dir_v.get(2).copied().unwrap_or(0.0),
                ];
                BcKind::CurrentTerminal { name, current, direction }
            }
            "flux_loop" => {
                let name = p.text("name")?.trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("flux_loop condition requires a non-empty 'name'"));
                }
                BcKind::FluxLoop { name, flux_wb: p.number("flux_wb")? }
            }
            other => return Err(Error::Unsupported(format!("dfs has no condition '{other}'"))),
        };
        if matches!(study, Study::HeatTransient { .. }) {
            if matches!(kind, BcKind::InfiniteElement(_)) {
                return Err(Error::invalid(
                    "infinite_element cannot be used with heat_transient: frequency-independent spatial coordinate stretching distorts volumetric heat capacity (rho * Cp * det(J_IE)) across transient diffusion timescales; use unstretched (s_u = 1) geometrically graded layers of diffusion-length thickness d = sqrt(4 * D_max * t_end * ln(1/R_0)) instead",
                ));
            }
            if matches!(kind, BcKind::FarField { .. }) {
                return Err(Error::invalid(
                    "far_field asymptotic Laplace Robin boundary condition is only valid for stationary heat transfer (heat_stationary), not heat_transient; use unstretched geometrically graded diffusion-length layers d = sqrt(4 * D_max * t_end * ln(1/R_0))",
                ));
            }
        }
        let ok = match &study {
            Study::Electrostatic => kind.is_electrostatic(),
            Study::Conduction => kind.is_conduction(),
            Study::HeatStationary | Study::HeatTransient { .. } => kind.is_heat(),
            Study::Magnetostatic { .. } => kind.is_magnetostatic(),
        };
        if !ok {
            let expected = match &study {
                Study::Electrostatic => "an electrostatics",
                Study::Conduction => "a steady-conduction",
                Study::HeatStationary | Study::HeatTransient { .. } => "a heat-transfer",
                Study::Magnetostatic { .. } => "a magnetostatics",
            };
            return Err(Error::invalid(format!("'{}' is not {expected} boundary condition", kind.name())));
        }
        bcs.push(Bc { labels: r.labels, kind });
    }
    if matches!(study, Study::Electrostatic) {
        if flux_bounds != "none" && singular_edges == "quarter_point" {
            return Err(Error::invalid(
                "flux_bounds cannot be combined with singular_edges='quarter_point': contravariant Piola-mapped Raviart-Thomas fluxes have divergent L^2 energy on singular quarter-point (det J -> 0) CurvedTets",
            ));
        }
        if flux_bounds != "none" && bcs.iter().any(|bc| matches!(bc.kind, BcKind::SurfaceCharge(s) if s.abs() > 0.0)) {
            return Err(Error::invalid(
                "flux_bounds requires source-free electrostatics (rho_v == 0 and surface_charge == 0)",
            ));
        }
        if formulation == "planar_bem" {
            if bcs.iter().any(|bc| matches!(bc.kind, BcKind::FarField { .. } | BcKind::InfiniteElement(_))) {
                return Err(Error::Unsupported(
                    "planar_bem models open half-spaces via its layered Green's function and does not support volumetric 'far_field' or 'infinite_element' conditions".into(),
                ));
            }
            if flux_bounds != "none" {
                return Err(Error::Unsupported(
                    "planar_bem is a 2D boundary-element formulation and does not support 3D volumetric Prager-Synge 'flux_bounds'".into(),
                ));
            }
            if ctx.resolver.sheet_thickness_m().values().any(|&t| t > 0.0)
                || ctx.resolver.sheet_etch_depth_m().values().any(|&t| t > 0.0)
            {
                return Err(Error::Unsupported(
                    "planar_bem requires zero-thickness planar conductor sheets on a flat interface and does not support finite sheet thickness > 0 or substrate etch_depth > 0".into(),
                ));
            }
        }
    }
    let params = dfs.params();
    let pp = Params::new(&spec.params, &params);
    let solver = match pp.text("solver")?.as_str() {
        "direct" => LinearSolver::Direct,
        "iterative" => LinearSolver::Iterative,
        _ => LinearSolver::Auto,
    };
    let tol = pp.number("tol")?;
    let amg_line_smooth = pp.boolean("amg_line_smooth")?;
    let tls_ifaces = tls_interfaces_from_spec(spec, ctx);
    let tls_interfaces = (!tls_ifaces.is_empty()).then_some(tls_ifaces);
    let materials = materials_by_label(spec, ctx, &dfs.properties())?;
    let mut k_curves: BTreeMap<i32, (TempTensor, [f64; 3])> = BTreeMap::new();
    let mut cp_curves: BTreeMap<i32, TempCurve> = BTreeMap::new();
    let mut sigma_curves: BTreeMap<i32, (TempTensor, [f64; 3])> = BTreeMap::new();
    if let Some(lib) = ctx.library.as_ref() {
        let is_nonconst_curve = |c: &TempCurve| !matches!(c, TempCurve::Constant { .. });
        let is_nonconst_tensor = |t: &TempTensor| match t {
            TempTensor::Scalar(c) => is_nonconst_curve(c),
            TempTensor::Uniaxial { perp, par } => is_nonconst_curve(perp) || is_nonconst_curve(par),
        };
        for (region, assign) in &spec.materials {
            let (mat_name, euler) = match assign {
                MaterialAssign::Named(n) => (n.as_str(), [0.0; 3]),
                MaterialAssign::Oriented { material, euler_deg } => (material.as_str(), *euler_deg),
                _ => continue,
            };
            let labels = ctx.resolver.domains(region)?;
            let mat = lib.require(mat_name).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            if matches!(study, Study::HeatStationary | Study::HeatTransient { .. }) {
                if let Some(kt) = mat.thermal.conductivity_t.as_ref().filter(|t| is_nonconst_tensor(t)) {
                    for &l in &labels {
                        k_curves.insert(l, (kt.clone(), euler));
                    }
                }
                if matches!(study, Study::HeatTransient { .. }) {
                    if let Some(cpt) = mat.thermal.specific_heat_t.as_ref().filter(|c| is_nonconst_curve(c)) {
                        for &l in &labels {
                            cp_curves.insert(l, cpt.clone());
                        }
                    }
                }
            } else if matches!(study, Study::Conduction) {
                if let Some(st) = mat.electrical.conductivity_t.as_ref().filter(|t| is_nonconst_tensor(t)) {
                    for &l in &labels {
                        sigma_curves.insert(l, (st.clone(), euler));
                    }
                }
            }
        }
    }
    let heat_source = if matches!(study, Study::HeatStationary | Study::HeatTransient { .. }) {
        if let Some(hs) = ctx.inputs.get("heat_source") {
            if hs.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFinite("imported 'heat_source' field contains NaN or infinity".into()));
            }
            Some(hs.clone())
        } else {
            None
        }
    } else {
        None
    };
    let temperature_field = if matches!(study, Study::Conduction) {
        if let Some(tf) = ctx.inputs.get("temperature") {
            if tf.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFinite("imported 'temperature' field contains NaN or infinity".into()));
            }
            Some(tf.clone())
        } else {
            None
        }
    } else {
        None
    };
    if matches!(study, Study::Electrostatic)
        && flux_bounds != "none"
        && (materials.default.get("rho_v").copied().unwrap_or(0.0).abs() > 0.0
            || materials.domains.values().any(|m| m.get("rho_v").copied().unwrap_or(0.0).abs() > 0.0))
    {
        return Err(Error::invalid(
            "flux_bounds requires source-free electrostatics (rho_v == 0 and surface_charge == 0)",
        ));
    }
    Ok(Problem {
        order: spec.order as usize,
        study,
        materials,
        domain_names: domain_names_by_label(spec, ctx)?,
        bcs,
        solver,
        tol,
        amg_line_smooth,
        formulation,
        singular_edges,
        richardson,
        tls_interfaces,
        adapt_passes,
        adapt_rel_tol,
        adapt_min_converged_passes,
        adapt_growth,
        adapt_mode,
        adapt_estimator,
        flux_bounds,
        freq,
        farfield_t_inf,
        sheet_thickness_m: ctx.resolver.sheet_thickness_m(),
        sheet_etch_depth_m: ctx.resolver.sheet_etch_depth_m(),
        k_curves,
        cp_curves,
        sigma_curves,
        heat_source,
        heat_source_cell_nodes: None,
        temperature_field,
        temperature_cell_nodes: None,
    })
}
