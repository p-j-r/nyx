// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `dfs`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (names and selectors
//! resolved to mesh numbers, parameters with their defaults), or built directly by code that already has mesh numbers (tests, the
//! legacy window bridge of the engine). The envelope carries scalar material properties only; anisotropic tensors are given in
//! `materials` directly until the material files exist (phase 6).

use crate::Dfs;
use femcore::material::{bcs_surface_impedance, FilmModel};
use femcore::props::{LinearSolver, Materials};
use nyx_physics::decode::{
    domain_names_by_label, materials_by_label, resolve_conditions, tls_interfaces_from_spec, Params,
};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec, TlsInterfaces};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum BcKind {
    Potential(f64),
    Ground,
    /// surface charge density, C/m^2
    SurfaceCharge(f64),
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
    Temperature(f64),
    /// inward heat flux, W/m^2
    Flux(f64),
    Convection {
        h: f64,
        text: f64,
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
            BcKind::Ground => "Ground",
            BcKind::SurfaceCharge(_) => "Surface charge",
            BcKind::Terminal(_) => "Terminal",
            BcKind::Open { .. } => "Open boundary",
            BcKind::FarField { .. } => "Far-field open boundary",
            BcKind::Temperature(_) => "Temperature",
            BcKind::Flux(_) => "Heat flux",
            BcKind::Convection { .. } => "Convective cooling",
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
        )
    }

    pub fn is_heat(&self) -> bool {
        matches!(self, BcKind::Temperature(_) | BcKind::Flux(_) | BcKind::Convection { .. })
    }

    pub fn is_magnetostatic(&self) -> bool {
        matches!(
            self,
            BcKind::LondonSheet { .. } | BcKind::CurrentTerminal { .. } | BcKind::FluxLoop { .. } | BcKind::Ground
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
    /// optional 1-based boundary label -> physical metal film thickness `t_m` (m) from `LabelMap`
    pub sheet_thickness_m: BTreeMap<i32, f64>,
    /// optional 1-based boundary label -> substrate over-etch trench depth `t_e` (m) from `LabelMap`
    pub sheet_etch_depth_m: BTreeMap<i32, f64>,
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
            sheet_thickness_m: BTreeMap::new(),
            sheet_etch_depth_m: BTreeMap::new(),
        }
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
        for (i, bc) in self.bcs.iter().enumerate() {
            let vals: Vec<f64> = match &bc.kind {
                BcKind::Potential(v) | BcKind::SurfaceCharge(v) | BcKind::Temperature(v) | BcKind::Flux(v) => vec![*v],
                BcKind::Convection { h, text } => vec![*h, *text],
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
        }
        if let Study::HeatTransient { t_end, .. } = &self.study {
            if !t_end.is_finite() {
                return Err("the study's end time is NaN or infinite".into());
            }
        }
        if matches!(self.study, Study::Electrostatic)
            && (!self.freq.is_finite() || !self.adapt_rel_tol.is_finite() || !self.adapt_growth.is_finite())
        {
            return Err(
                "electrostatic study parameter (freq, adapt_rel_tol, or adapt_growth) is NaN or infinite".into()
            );
        }
        // tensor data of anisotropic domains
        let text = format!("{:?}", m.aniso);
        if text.contains("NaN") || text.contains("inf") {
            return Err("the problem description contains a NaN or infinite number in an anisotropic material".into());
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
    let study = match sd.kind.as_str() {
        "electrostatic" => {
            formulation = sp.text("formulation");
            singular_edges = sp.text("singular_edges");
            richardson = sp.boolean("richardson");
            let ap = sp.integer("adapt_passes");
            if ap < 1 {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_passes' must be >= 1"));
            }
            adapt_passes = ap as usize;
            adapt_rel_tol = sp.number("adapt_rel_tol");
            if !(adapt_rel_tol.is_finite() && adapt_rel_tol > 0.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_rel_tol' must be > 0"));
            }
            let amc = sp.integer("adapt_min_converged_passes");
            if amc < 1 {
                return Err(Error::invalid(
                    "study 'electrostatic': parameter 'adapt_min_converged_passes' must be >= 1",
                ));
            }
            adapt_min_converged_passes = amc as usize;
            adapt_growth = sp.number("adapt_growth");
            if !(adapt_growth.is_finite() && adapt_growth > 1.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'adapt_growth' must be > 1"));
            }
            adapt_mode = sp.text("adapt_mode");
            adapt_estimator = sp.text("adapt_estimator");
            flux_bounds = sp.text("flux_bounds");
            freq = sp.number("freq");
            if !(freq.is_finite() && freq >= 0.0) {
                return Err(Error::invalid("study 'electrostatic': parameter 'freq' must be >= 0"));
            }
            Study::Electrostatic
        }
        "heat_stationary" => Study::HeatStationary,
        "heat_transient" => Study::HeatTransient {
            t_end: sp.number("t_end"),
            n_out: sp.integer("n_out").max(0) as usize,
            t0: sp.text("t0"),
        },
        "magnetostatic" => Study::Magnetostatic { formulation: sp.text("formulation") },
        other => return Err(Error::Unsupported(format!("dfs has no study '{other}'"))),
    };
    let kinds = dfs.bc_kinds();
    let mut bcs = vec![];
    for r in resolve_conditions(&spec.bcs, &kinds, ctx)? {
        let d = kinds.iter().find(|k| k.kind == r.kind).expect("resolved against these descriptors");
        let p = Params::new(&r.params, &d.params);
        let kind = match r.kind.as_str() {
            "potential" => BcKind::Potential(p.number("value")),
            "ground" => BcKind::Ground,
            "surface_charge" => BcKind::SurfaceCharge(p.number("density")),
            "terminal" => {
                let name = p.text("name").trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("terminal condition requires a non-empty 'name'"));
                }
                BcKind::Terminal(name)
            }
            "open" => {
                let origin = if p.has("origin") {
                    let v = p.vector("origin");
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
                let order_i = p.integer("order");
                if !(1..=2).contains(&order_i) {
                    return Err(Error::invalid(format!("far_field condition 'order' must be 1 or 2, got {order_i}")));
                }
                let decay_power = p.number("decay_power");
                if !decay_power.is_finite() || decay_power <= 0.0 {
                    return Err(Error::invalid(format!(
                        "far_field condition 'decay_power' must be positive, got {decay_power}"
                    )));
                }
                let center = if p.has("center") {
                    let cv = p.vector("center");
                    if cv.len() != 3 || cv.iter().any(|x| !x.is_finite()) {
                        return Err(Error::invalid(
                            "far_field condition 'center' must be 3 finite numbers [cx, cy, cz]",
                        ));
                    }
                    Some([cv[0], cv[1], cv[2]])
                } else {
                    None
                };
                BcKind::FarField { order: order_i as usize, decay_power, center }
            }
            "temperature" => BcKind::Temperature(p.number("value")),
            "flux" => BcKind::Flux(p.number("value")),
            "convection" => BcKind::Convection { h: p.number("h"), text: p.number("t_ext") },
            "london_sheet" => {
                let l_k_sq_given = p.number("L_k_sq");
                let mat_name = p.text("material");
                let mat_trim = mat_name.trim();
                let l_k_sq = if !mat_trim.is_empty() {
                    let thickness = p.number("thickness");
                    let temperature = p.number("temperature");
                    let fm_str = p.text("film_model");
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
                let name = p.text("name").trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("current_terminal condition requires a non-empty 'name'"));
                }
                let current = p.number("current");
                let dir_v = p.vector("direction");
                let direction = [
                    dir_v.first().copied().unwrap_or(1.0),
                    dir_v.get(1).copied().unwrap_or(0.0),
                    dir_v.get(2).copied().unwrap_or(0.0),
                ];
                BcKind::CurrentTerminal { name, current, direction }
            }
            "flux_loop" => {
                let name = p.text("name").trim().to_string();
                if name.is_empty() {
                    return Err(Error::invalid("flux_loop condition requires a non-empty 'name'"));
                }
                BcKind::FluxLoop { name, flux_wb: p.number("flux_wb") }
            }
            other => return Err(Error::Unsupported(format!("dfs has no condition '{other}'"))),
        };
        let ok = match &study {
            Study::Electrostatic => kind.is_electrostatic(),
            Study::HeatStationary | Study::HeatTransient { .. } => kind.is_heat(),
            Study::Magnetostatic { .. } => kind.is_magnetostatic(),
        };
        if !ok {
            let expected = match &study {
                Study::Electrostatic => "an electrostatics",
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
            if bcs.iter().any(|bc| matches!(bc.kind, BcKind::FarField { .. })) {
                return Err(Error::Unsupported(
                    "planar_bem models open half-spaces via its layered Green's function and does not support volumetric 'far_field' boundary conditions".into(),
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
    let solver = match pp.text("solver").as_str() {
        "direct" => LinearSolver::Direct,
        "iterative" => LinearSolver::Iterative,
        _ => LinearSolver::Auto,
    };
    let tol = pp.number("tol");
    let amg_line_smooth = pp.boolean("amg_line_smooth");
    let tls_ifaces = tls_interfaces_from_spec(spec, ctx);
    let tls_interfaces = (!tls_ifaces.is_empty()).then_some(tls_ifaces);
    let materials = materials_by_label(spec, ctx, &dfs.properties())?;
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
        sheet_thickness_m: ctx.resolver.sheet_thickness_m(),
        sheet_etch_depth_m: ctx.resolver.sheet_etch_depth_m(),
    })
}
