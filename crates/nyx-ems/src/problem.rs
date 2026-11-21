// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `ems`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (region names and selectors
//! resolved to mesh numbers, parameters with their defaults), or built directly by code that already has mesh numbers (tests, the legacy
//! window bridge of the engine). The envelope carries scalar permittivity and permeability; a crystal's permittivity tensor reaches a kernel
//! through `materials.aniso` until the material files exist (phase 6).
//!
//! The study names the kernel: `Modes` (guided modes of a 2D cross-section, real or complex), `AxisymmetricModes` (a body of revolution on
//! its (r, z) half plane) or `Cavity` (3D eigenfrequencies). It is not guessed from the dimension of the mesh (ROUGH_SPOTS R10): a study
//! on a mesh of the other dimension is refused by name.

use crate::ibc::resolve_jj_l_and_ej;
use crate::Ems;
use femcore::material::{bcs_surface_impedance, FilmModel};
use femcore::props::{Materials, Periodic, Pml, PmlDir};
use nyx_physics::decode::{domain_names_by_label, materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec};
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum Study {
    /// Guided modes at `freq` (Hz) near the effective index `n_target` (imaginary if `evanescent`). Absorbing layers, loss and a crystal
    /// tilted against the guide take the complex solver.
    Modes { freq: f64, n_target: f64, evanescent: bool, n_modes: usize },
    /// Modes of a body of revolution at azimuthal order `m` (`exp(-j m theta)`, theta along wafer y) nearest `f_search` (Hz); the mesh is
    /// the (r, z) half plane with x = r. A permittivity per domain (isotropic), absorbing layers by `materials.pml`.
    AxisymmetricModes { f_search: f64, azimuthal_order: i64, n_modes: usize },
    /// Eigenfrequencies of a 3D cavity or Floquet cell nearest `f_search` (Hz; a value of zero or less takes the fundamental of a box of
    /// the size of the mesh).
    Cavity { n_modes: usize, f_search: f64 },
    /// Driven frequency-domain multi-port S/Z/Y sweep over `[f_min, f_max]` (`n_points`), solved with Hermite/Galerkin PROM
    /// (`prom_tol = Some(tol)`) or direct point-by-point HDM (`prom_tol = None`).
    Driven { f_min: f64, f_max: f64, n_points: usize, prom_tol: Option<f64> },
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Modes { .. } => "modes",
            Study::AxisymmetricModes { .. } => "axisymmetric_modes",
            Study::Cavity { .. } => "cavity",
            Study::Driven { .. } => "driven",
        }
    }

    /// Reference frequency of the study (Hz): `freq` for `Modes`, `f_search` for `AxisymmetricModes` and `Cavity`, midband for `Driven`.
    pub fn reference_freq(&self) -> f64 {
        match *self {
            Study::Modes { freq, .. } => freq,
            Study::AxisymmetricModes { f_search, .. } | Study::Cavity { f_search, .. } => f_search,
            Study::Driven { f_min, f_max, .. } => 0.5 * (f_min + f_max),
        }
    }
}

/// Kinetic-inductance / surface-impedance sheet boundary condition (`Z_s = R_s + i omega L_k_sq`).
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceImpedanceBc {
    pub labels: Vec<i32>,
    pub l_k_sq: f64,
    pub r_s: f64,
    pub material: String,
    pub thickness: f64,
    pub temperature: f64,
    pub f_ref: f64,
    pub film_model: String,
}

/// Lumped Josephson junction / nonlinear port boundary condition (`cavity` and `driven` studies).
#[derive(Clone, Debug, PartialEq)]
pub struct LumpedJjBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub l_j: f64,
    pub e_j_hz: f64,
    pub c_j: f64,
    pub e_l_hz: f64,
    pub phi_ext_rad: f64,
    pub n_g: f64,
    pub direction: [f64; 3],
    pub formulation: String,
}

/// Lumped sheet port boundary condition (`driven` study).
#[derive(Clone, Debug, PartialEq)]
pub struct LumpedPortBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub z0: f64,
    pub direction: [f64; 3],
}

/// 2D boundary-mode waveguide port boundary condition (`driven` study).
#[derive(Clone, Debug, PartialEq)]
pub struct WavePortBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub mode: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Problem {
    /// element order, 1 or 2
    pub order: usize,
    pub study: Study,
    /// `eps_r` and `mu_r` per domain (scalars, or crystal tensors in `aniso`), absorbing layers in `pml`
    pub materials: Materials,
    /// Floquet-periodic boundary pairs (cavity study)
    pub periodic: Vec<Periodic>,
    /// boundary numbers of magnetic walls (cavity study): perfect magnetic conductor, the natural condition of the curl-curl form, so
    /// the field is left free on them instead of having its tangential part set to zero. A mirror plane that keeps the modes even
    /// under the reflection (tangential E even, E normal odd) is one of these.
    pub pmc: Vec<i32>,
    /// boundary numbers of explicit PEC walls or interior 2D PEC sheets
    pub pec: Vec<i32>,
    /// superconducting / resistive surface-impedance boundary conditions
    pub surface_impedance: Vec<SurfaceImpedanceBc>,
    /// lumped Josephson junction boundary conditions (cavity / driven study)
    pub lumped_jj: Vec<LumpedJjBc>,
    /// lumped port boundary conditions (driven study)
    pub lumped_ports: Vec<LumpedPortBc>,
    /// waveguide port boundary conditions (driven study)
    pub waveports: Vec<WavePortBc>,
    /// optional 1-based domain label -> region name mapping from `ProblemSpec::materials`
    pub domain_names: BTreeMap<i32, String>,
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem {
            order,
            study,
            materials: Materials::default(),
            periodic: vec![],
            pmc: vec![],
            pec: vec![],
            surface_impedance: vec![],
            lumped_jj: vec![],
            lumped_ports: vec![],
            waveports: vec![],
            domain_names: BTreeMap::new(),
        }
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
        let study_bad = match &self.study {
            Study::Modes { freq, n_target, .. } => (!freq.is_finite())
                .then_some("frequency")
                .or((!n_target.is_finite()).then_some("target effective index")),
            Study::AxisymmetricModes { f_search, .. } | Study::Cavity { f_search, .. } => {
                (!f_search.is_finite()).then_some("search frequency")
            }
            Study::Driven { f_min, f_max, prom_tol, .. } => (!f_min.is_finite())
                .then_some("minimum frequency")
                .or((!f_max.is_finite()).then_some("maximum frequency"))
                .or(prom_tol.is_some_and(|t| !t.is_finite()).then_some("PROM tolerance")),
        };
        if let Some(what) = study_bad {
            return Err(format!("the study's {what} is NaN or infinite"));
        }
        for ibc in &self.surface_impedance {
            for (k, v) in [
                ("L_k_sq", ibc.l_k_sq),
                ("R_s", ibc.r_s),
                ("thickness", ibc.thickness),
                ("temperature", ibc.temperature),
                ("f_ref", ibc.f_ref),
            ] {
                if !v.is_finite() {
                    return Err(format!("surface_impedance parameter '{k}' is {v}: it must be a finite number"));
                }
            }
        }
        for jj in &self.lumped_jj {
            for (k, v) in [
                ("L_j", jj.l_j),
                ("E_j_Hz", jj.e_j_hz),
                ("C_j", jj.c_j),
                ("E_l_Hz", jj.e_l_hz),
                ("phi_ext_rad", jj.phi_ext_rad),
                ("n_g", jj.n_g),
                ("direction[0]", jj.direction[0]),
                ("direction[1]", jj.direction[1]),
                ("direction[2]", jj.direction[2]),
            ] {
                if !v.is_finite() {
                    return Err(format!("lumped_jj '{}' parameter '{k}' is {v}: it must be a finite number", jj.name));
                }
            }
        }
        for lp in &self.lumped_ports {
            for (k, v) in [
                ("Z_0", lp.z0),
                ("direction[0]", lp.direction[0]),
                ("direction[1]", lp.direction[1]),
                ("direction[2]", lp.direction[2]),
            ] {
                if !v.is_finite() {
                    return Err(format!(
                        "lumped_port '{}' parameter '{k}' is {v}: it must be a finite number",
                        lp.name
                    ));
                }
            }
        }
        // tensors of anisotropic domains, layers and Floquet phases: scan the printed form, outside quoted strings, for NaN and inf
        let text = format!("{:?} {:?} {:?}", m.aniso, m.pml, self.periodic);
        if text.contains("NaN") || text.contains("inf") {
            return Err(
                "the problem description contains a NaN or infinite number (an anisotropic material, a layer or a Floquet phase)"
                    .into(),
            );
        }
        Ok(())
    }
}

fn pml_dir(text: &str) -> PmlDir {
    match text {
        "plus" => PmlDir::Plus,
        "minus" => PmlDir::Minus,
        _ => PmlDir::Off,
    }
}

/// Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and
/// selectors with `ctx.resolver` and refuses combinations that need the study (a Floquet pair outside the cavity study).
pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error> {
    let ems = Ems;
    let studies = ems.studies();
    let sd = studies
        .iter()
        .find(|s| s.kind == spec.study.kind)
        .ok_or_else(|| Error::Unsupported(format!("ems has no study '{}'", spec.study.kind)))?;
    let sp = Params::new(&spec.study.params, &sd.params);
    let study = match sd.kind.as_str() {
        "modes" => Study::Modes {
            freq: sp.number("freq"),
            n_target: sp.number("n_target"),
            evanescent: sp.boolean("evanescent"),
            n_modes: sp.integer("n_modes").max(0) as usize,
        },
        "axisymmetric_modes" => Study::AxisymmetricModes {
            f_search: sp.number("f_search"),
            azimuthal_order: sp.integer("azimuthal_order"),
            n_modes: sp.integer("n_modes").max(0) as usize,
        },
        "cavity" => Study::Cavity { n_modes: sp.integer("n_modes").max(0) as usize, f_search: sp.number("f_search") },
        "driven" => {
            let f_min = sp.number("f_min");
            let f_max = sp.number("f_max");
            if !(f_min > 0.0 && f_max >= f_min) {
                return Err(Error::invalid(format!(
                    "driven study requires 0 < f_min <= f_max (got f_min = {f_min}, f_max = {f_max})"
                )));
            }
            let n_points = sp.integer("n_points").max(1) as usize;
            let prom_tol = sp.optional_number("prom_tol").filter(|&t| t > 0.0);
            Study::Driven { f_min, f_max, n_points, prom_tol }
        }
        other => return Err(Error::Unsupported(format!("ems has no study '{other}'"))),
    };
    let mut materials = materials_by_label(spec, ctx, &ems.properties())?;
    let domain_names = domain_names_by_label(spec, ctx)?;
    let mut periodic = vec![];
    let mut pmc: Vec<i32> = vec![];
    let mut pec: Vec<i32> = vec![];
    let mut surface_impedance: Vec<SurfaceImpedanceBc> = vec![];
    let mut lumped_jj: Vec<LumpedJjBc> = vec![];
    let mut lumped_ports: Vec<LumpedPortBc> = vec![];
    let mut waveports: Vec<WavePortBc> = vec![];
    let kinds = ems.bc_kinds();
    for r in resolve_conditions(&spec.bcs, &kinds, ctx)? {
        let d = kinds.iter().find(|k| k.kind == r.kind).expect("resolved against these descriptors");
        let p = Params::new(&r.params, &d.params);
        match r.kind.as_str() {
            "pml" => {
                let layer = Pml {
                    x: pml_dir(&p.text("x")),
                    y: pml_dir(&p.text("y")),
                    z: pml_dir(&p.text("z")),
                    strength: p.number("strength"),
                    order: p.number("order"),
                    re: p.number("re"),
                };
                for l in r.labels {
                    materials.pml.insert(l, layer.clone());
                }
            }
            "periodic" => {
                if !matches!(study, Study::Cavity { .. } | Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a Floquet-periodic pair belongs to the cavity study, not to '{}'",
                        study.kind()
                    )));
                }
                let dest = ctx.resolver.boundaries(&p.text("dest"))?;
                if dest.is_empty() {
                    return Err(Error::invalid(format!(
                        "the destination '{}' of a periodic pair selects nothing",
                        p.text("dest")
                    )));
                }
                periodic.push(Periodic { source: r.labels, dest, phase: p.number("phase") });
            }
            "pmc" => {
                if !matches!(study, Study::Cavity { .. } | Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a magnetic wall belongs to the cavity study, not to '{}'",
                        study.kind()
                    )));
                }
                pmc.extend(r.labels);
            }
            "pec" => {
                pec.extend(r.labels);
            }
            "surface_impedance" => {
                let mut l_k_sq = p.number("L_k_sq");
                let mut r_s = p.number("R_s");
                let material = p.text("material");
                let thickness = p.number("thickness");
                let temperature = p.number("temperature");
                let f_ref_given = p.optional_number("f_ref").unwrap_or(0.0);
                let f_ref = if f_ref_given > 0.0 { f_ref_given } else { study.reference_freq() };
                let film_model = p.text("film_model");
                let fm = FilmModel::parse(&film_model).ok_or_else(|| {
                    Error::invalid(format!(
                        "surface_impedance: unknown film_model '{film_model}' (expected 'two_sided', 'kerr' or 'thin_limit')"
                    ))
                })?;
                if !material.trim().is_empty() {
                    let lib = ctx.library.as_ref().ok_or_else(|| {
                        Error::Unsupported(format!(
                            "surface_impedance: material '{material}' needs a material library, and this run has none"
                        ))
                    })?;
                    let mat =
                        lib.require(material.trim()).map_err(|e| Error::invalid(format!("surface_impedance: {e}")))?;
                    let sc = mat.superconducting.as_ref().ok_or_else(|| {
                        Error::invalid(format!(
                            "surface_impedance: material '{}' has no [superconducting] section",
                            material.trim()
                        ))
                    })?;
                    if f_ref <= 0.0 {
                        return Err(Error::invalid(
                            "surface_impedance with a superconducting material requires a positive f_ref or study search frequency",
                        ));
                    }
                    let eval = bcs_surface_impedance(sc, f_ref, temperature, thickness, fm)
                        .map_err(|e| Error::invalid(format!("surface_impedance ('{}'): {e}", material.trim())))?;
                    l_k_sq = eval.l_k_sq;
                    r_s = eval.r_s;
                }
                if l_k_sq < 0.0 || r_s < 0.0 {
                    return Err(Error::invalid("surface_impedance requires non-negative L_k_sq and R_s"));
                }
                surface_impedance.push(SurfaceImpedanceBc {
                    labels: r.labels,
                    l_k_sq,
                    r_s,
                    material,
                    thickness,
                    temperature,
                    f_ref,
                    film_model: fm.as_str().to_string(),
                });
            }
            "lumped_jj" => {
                if !matches!(study, Study::Cavity { .. } | Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a lumped_jj condition belongs to the cavity or driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name");
                if name.trim().is_empty() {
                    return Err(Error::invalid("lumped_jj requires a non-empty 'name'"));
                }
                if lumped_jj.iter().any(|j| j.name == name) {
                    return Err(Error::invalid(format!("duplicate lumped_jj name '{name}'")));
                }
                let l_j_in = p.number("L_j");
                let e_j_in = p.number("E_j_Hz");
                let c_j = p.number("C_j");
                let e_l_hz = p.number("E_l_Hz");
                let phi_ext_rad = p.number("phi_ext_rad");
                let n_g = p.number("n_g");
                let (l_j, e_j_hz) = resolve_jj_l_and_ej(l_j_in, e_j_in);
                if l_j <= 0.0 && e_l_hz <= 0.0 && c_j <= 0.0 {
                    return Err(Error::invalid(format!(
                        "lumped_jj '{name}' requires a positive L_j, E_j_Hz, E_l_Hz, or C_j"
                    )));
                }
                let dir_vec = p.vector("direction");
                if dir_vec.len() != 3 {
                    return Err(Error::invalid(format!("lumped_jj '{name}': direction must have 3 components")));
                }
                let dnorm = (dir_vec[0] * dir_vec[0] + dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]).sqrt();
                if !(dnorm > 0.0 && dnorm.is_finite()) {
                    return Err(Error::invalid(format!("lumped_jj '{name}': direction vector must be non-zero")));
                }
                let direction = [dir_vec[0] / dnorm, dir_vec[1] / dnorm, dir_vec[2] / dnorm];
                let formulation = p.text("formulation");
                if formulation != "rank1" && formulation != "distributed_sheet" {
                    return Err(Error::invalid(format!(
                        "lumped_jj '{name}': unknown formulation '{formulation}' (expected 'rank1' or 'distributed_sheet')"
                    )));
                }
                lumped_jj.push(LumpedJjBc {
                    labels: r.labels,
                    name,
                    l_j,
                    e_j_hz,
                    c_j,
                    e_l_hz,
                    phi_ext_rad,
                    n_g,
                    direction,
                    formulation,
                });
            }
            "lumped_port" => {
                if !matches!(study, Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a lumped_port condition belongs to the driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name");
                if name.trim().is_empty() {
                    return Err(Error::invalid("lumped_port requires a non-empty 'name'"));
                }
                if lumped_ports.iter().any(|lp| lp.name == name) || waveports.iter().any(|wp| wp.name == name) {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let z0 = p.number("Z_0");
                if !(z0 > 0.0 && z0.is_finite()) {
                    return Err(Error::invalid(format!("lumped_port '{name}': Z_0 must be positive (got {z0})")));
                }
                let dir_vec = p.vector("direction");
                if dir_vec.len() != 3 {
                    return Err(Error::invalid(format!("lumped_port '{name}': direction must have 3 components")));
                }
                let dnorm = (dir_vec[0] * dir_vec[0] + dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]).sqrt();
                if !(dnorm > 0.0 && dnorm.is_finite()) {
                    return Err(Error::invalid(format!("lumped_port '{name}': direction vector must be non-zero")));
                }
                let direction = [dir_vec[0] / dnorm, dir_vec[1] / dnorm, dir_vec[2] / dnorm];
                lumped_ports.push(LumpedPortBc { labels: r.labels, name, z0, direction });
            }
            "waveport" => {
                if !matches!(study, Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a waveport condition belongs to the driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name");
                if name.trim().is_empty() {
                    return Err(Error::invalid("waveport requires a non-empty 'name'"));
                }
                if lumped_ports.iter().any(|lp| lp.name == name) || waveports.iter().any(|wp| wp.name == name) {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let mode = p.integer("mode").max(1) as usize;
                waveports.push(WavePortBc { labels: r.labels, name, mode });
            }
            other => return Err(Error::Unsupported(format!("ems has no condition '{other}'"))),
        }
    }
    pmc.sort_unstable();
    pmc.dedup();
    pec.sort_unstable();
    pec.dedup();
    Ok(Problem {
        order: spec.order as usize,
        study,
        materials,
        periodic,
        pmc,
        pec,
        surface_impedance,
        lumped_jj,
        lumped_ports,
        waveports,
        domain_names,
    })
}
