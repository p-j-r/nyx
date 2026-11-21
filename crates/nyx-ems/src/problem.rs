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
use femcore::props::{Materials, Periodic, Pml, PmlDir, TensorKind};
use nyx_physics::decode::{domain_names_by_label, materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Discretisation, Error, MaterialAssign, Physics, ProblemSpec};
use std::collections::BTreeMap;

/// Allowed values for the midplane `symmetry` parameter on `modes`, `cavity`, and `bloch`.
pub const SYMMETRY_CHOICES: [&str; 9] =
    ["none", "pec_x", "pmc_x", "pec_z", "pmc_z", "pec_x_pmc_z", "pmc_x_pec_z", "pmc_x_pmc_z", "pec_x_pec_z"];

/// Midplane PEC/PMC symmetry condition on exterior facets lying on `x = 0` and/or `z = 0`.
/// On a quarter mesh `x >= 0, z >= 0`, fundamental `TE0` is `PecXPmcZ` (`pec_x_pmc_z`) and `TM0` is `PmcXPecZ` (`pmc_x_pec_z`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Symmetry {
    #[default]
    None,
    PecX,
    PmcX,
    PecZ,
    PmcZ,
    PecXPmcZ,
    PmcXPecZ,
    PmcXPmcZ,
    PecXPecZ,
}

impl Symmetry {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "none" | "" => Ok(Self::None),
            "pec_x" => Ok(Self::PecX),
            "pmc_x" => Ok(Self::PmcX),
            "pec_z" => Ok(Self::PecZ),
            "pmc_z" => Ok(Self::PmcZ),
            "pec_x_pmc_z" => Ok(Self::PecXPmcZ),
            "pmc_x_pec_z" => Ok(Self::PmcXPecZ),
            "pmc_x_pmc_z" => Ok(Self::PmcXPmcZ),
            "pec_x_pec_z" => Ok(Self::PecXPecZ),
            other => Err(format!("unknown symmetry '{other}' (expected one of {SYMMETRY_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::PecX => "pec_x",
            Self::PmcX => "pmc_x",
            Self::PecZ => "pec_z",
            Self::PmcZ => "pmc_z",
            Self::PecXPmcZ => "pec_x_pmc_z",
            Self::PmcXPecZ => "pmc_x_pec_z",
            Self::PmcXPmcZ => "pmc_x_pmc_z",
            Self::PecXPecZ => "pec_x_pec_z",
        }
    }

    /// `Some(true)` for PEC on `x = 0`, `Some(false)` for PMC on `x = 0`, `None` otherwise.
    pub fn x_is_pec(self) -> Option<bool> {
        match self {
            Self::PecX | Self::PecXPmcZ | Self::PecXPecZ => Some(true),
            Self::PmcX | Self::PmcXPecZ | Self::PmcXPmcZ => Some(false),
            Self::None | Self::PecZ | Self::PmcZ => None,
        }
    }

    /// `Some(true)` for PEC on `z = 0`, `Some(false)` for PMC on `z = 0`, `None` otherwise.
    pub fn z_is_pec(self) -> Option<bool> {
        match self {
            Self::PecZ | Self::PmcXPecZ | Self::PecXPecZ => Some(true),
            Self::PmcZ | Self::PecXPmcZ | Self::PmcXPmcZ => Some(false),
            Self::None | Self::PecX | Self::PmcX => None,
        }
    }
}

/// Auto-detects exterior facets lying on `x = 0` (when `|xmin| <= tol` or `|xmax| <= tol`) and/or
/// `z = 0` (when `|zmin| <= tol` or `|zmax| <= tol`), returning `(sym_pec_facets, sym_pmc_facets)`.
pub fn symmetry_facet_flags(mesh: &Discretisation, symmetry: Symmetry) -> Result<(Vec<bool>, Vec<bool>), String> {
    let nf = mesh.n_facets();
    let mut sym_pec = vec![false; nf];
    let mut sym_pmc = vec![false; nf];
    if symmetry == Symmetry::None {
        return Ok((sym_pec, sym_pmc));
    }
    let (lo, hi) = mesh.bbox();
    let extent = (0..mesh.dim).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let x_idx = 0usize;
    let z_idx = if mesh.dim == 2 { 1usize } else { 2usize };
    let x_rule = symmetry.x_is_pec();
    let z_rule = symmetry.z_is_pec();
    if x_rule.is_some() && !(lo[x_idx].abs() <= tol || hi[x_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on x = 0 (mesh x spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[x_idx],
            hi[x_idx]
        ));
    }
    if z_rule.is_some() && !(lo[z_idx].abs() <= tol || hi[z_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on z = 0 (mesh z spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[z_idx],
            hi[z_idx]
        ));
    }
    let exterior = mesh.exterior_facets();
    for f in 0..nf {
        if !exterior[f] {
            continue;
        }
        let fv = mesh.facet(f);
        if let Some(is_pec) = x_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[x_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
        if let Some(is_pec) = z_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[z_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
    }
    Ok((sym_pec, sym_pmc))
}

/// Collects PEC boundary edge indices `be` and vertex indices `bn` for a 2D mesh, honoring explicit `pec`
/// facets and midplane `symmetry` (`pec_*` / `pmc_*`), while excluding `surface_impedance` and `pmc` facets.
pub(crate) fn pec_edges_and_nodes_2d_with_symmetry(
    mesh: &Discretisation,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
    symmetry: Symmetry,
) -> Result<(Vec<usize>, Vec<usize>), String> {
    let exterior = mesh.exterior_facets();
    let (sym_pec, sym_pmc) = symmetry_facet_flags(mesh, symmetry)?;
    let ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let is_pec_facet = |f: usize| -> bool {
        let l = mesh.facet_labels[f];
        pec_labels.contains(&l)
            || sym_pec[f]
            || (exterior[f] && !ibc_labels.contains(&l) && !pmc_labels.contains(&l) && !sym_pmc[f])
    };
    let mut be: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
        .collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
        .collect();
    bn.sort_unstable();
    bn.dedup();
    Ok((be, bn))
}

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
    /// Fixed-frequency 3D Floquet-Bloch modes at `freq` (Hz) near the Bloch index `n_target` in Brillouin zone `brillouin_zone`.
    Bloch { freq: f64, n_target: f64, n_modes: usize, brillouin_zone: i64 },
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
            Study::Bloch { .. } => "bloch",
            Study::Driven { .. } => "driven",
        }
    }

    /// Reference frequency of the study (Hz): `freq` for `Modes` and `Bloch`, `f_search` for `AxisymmetricModes` and `Cavity`, midband for `Driven`.
    pub fn reference_freq(&self) -> f64 {
        match *self {
            Study::Modes { freq, .. } | Study::Bloch { freq, .. } => freq,
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
    /// optional midplane PEC/PMC symmetry condition on `x = 0` and/or `z = 0`
    pub symmetry: Symmetry,
    /// global-frame real dispersive permittivity derivative `Re(d(omega * eps_r) / d(omega))` per domain label
    /// when evaluated from a material's `[optical]` dispersion model; domains absent here fall back to static `Re(eps_r)`
    pub d_omega_eps: BTreeMap<i32, [[f64; 3]; 3]>,
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
            symmetry: Symmetry::None,
            d_omega_eps: BTreeMap::new(),
        }
    }

    /// Global-frame unstretched `Re(d(omega * eps_r) / d(omega))` tensor for every cell of `mesh`,
    /// falling back to static `Re(eps_r)` on domains without an `[optical]` dispersion entry in `d_omega_eps`.
    pub fn cell_d_omega_eps(&self, mesh: &Discretisation) -> Result<Vec<[[f64; 3]; 3]>, String> {
        let static_re = self.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?;
        Ok((0..mesh.n_cells())
            .map(|ci| {
                let dom = mesh.cell_labels[ci];
                self.d_omega_eps.get(&dom).copied().unwrap_or(static_re[ci])
            })
            .collect())
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
            Study::Bloch { freq, n_target, .. } => {
                (!freq.is_finite()).then_some("frequency").or((!n_target.is_finite()).then_some("target Bloch index"))
            }
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
        for (&dom, mat) in &self.d_omega_eps {
            if mat.iter().flatten().any(|v| !v.is_finite()) {
                return Err(format!("dispersive permittivity derivative of domain {dom} contains NaN or infinity"));
            }
        }
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
        "bloch" => {
            let freq = sp.number("freq");
            if freq <= 0.0 || !freq.is_finite() {
                return Err(Error::invalid(format!("bloch study requires freq > 0 (got {freq})")));
            }
            Study::Bloch {
                freq,
                n_target: sp.number("n_target"),
                n_modes: sp.integer("n_modes").max(0) as usize,
                brillouin_zone: sp.integer("brillouin_zone"),
            }
        }
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
    let symmetry = if matches!(sd.kind.as_str(), "modes" | "cavity" | "bloch") {
        Symmetry::parse(&sp.text("symmetry")).map_err(Error::invalid)?
    } else {
        Symmetry::None
    };

    let f_ref = study.reference_freq();
    let derived_lambda = (f_ref > 0.0 && f_ref.is_finite()).then(|| crate::modes::C0 / f_ref);
    let mut ctx_eff = Ctx::new(ctx.resolver);
    ctx_eff.library = ctx.library.clone();
    ctx_eff.wavelength = ctx.wavelength.or(derived_lambda);
    ctx_eff.frequency = ctx.frequency.or_else(|| (f_ref > 0.0 && f_ref.is_finite()).then_some(f_ref));
    ctx_eff.temperature = ctx.temperature;
    ctx_eff.inputs = ctx.inputs.clone();
    ctx_eff.progress = ctx.progress.clone();

    let mut materials = materials_by_label(spec, &ctx_eff, &ems.properties())?;
    let domain_names = domain_names_by_label(spec, &ctx_eff)?;
    let temp_k = ctx_eff.temperature.or_else(|| {
        spec.params.get("temperature").and_then(serde_json::Value::as_f64).filter(|&t| t > 0.0 && t.is_finite())
    });
    let mut d_omega_eps: BTreeMap<i32, [[f64; 3]; 3]> = BTreeMap::new();
    if let (Some(lib), Some(lambda_m)) = (ctx_eff.library.as_ref(), ctx_eff.wavelength.filter(|&lam| lam <= 3e-4)) {
        for (region, assign) in &spec.materials {
            let (name, euler) = match assign {
                MaterialAssign::Named(n) => (n.as_str(), [0.0; 3]),
                MaterialAssign::Oriented { material, euler_deg } => (material.as_str(), *euler_deg),
                _ => continue,
            };
            let labels = ctx_eff.resolver.domains(region)?;
            let mat = lib.require(name).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            let Some(opt) = mat.optical.as_ref() else { continue };
            let mut diag = [0.0_f64; 3];
            for (axis, slot) in diag.iter_mut().enumerate() {
                let (n, k, dn_dlam, dk_dlam) =
                    opt.nk(axis, lambda_m, temp_k).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
                let n_tilde = femcore::linalg::csr::c(n, -k);
                let dn_tilde = femcore::linalg::csr::c(dn_dlam, -dk_dlam);
                let d_om = n_tilde * n_tilde - n_tilde * dn_tilde * (2.0 * lambda_m);
                *slot = d_om.re;
            }
            let cryst = [[diag[0], 0.0, 0.0], [0.0, diag[1], 0.0], [0.0, 0.0, diag[2]]];
            let glob = femcore::material::rotate_tensor(&cryst, &femcore::material::rotation(euler));
            for &l in &labels {
                d_omega_eps.insert(l, glob);
            }
        }
    }

    let mut periodic = vec![];
    let mut pmc: Vec<i32> = vec![];
    let mut pec: Vec<i32> = vec![];
    let mut surface_impedance: Vec<SurfaceImpedanceBc> = vec![];
    let mut lumped_jj: Vec<LumpedJjBc> = vec![];
    let mut lumped_ports: Vec<LumpedPortBc> = vec![];
    let mut waveports: Vec<WavePortBc> = vec![];
    let kinds = ems.bc_kinds();
    for r in resolve_conditions(&spec.bcs, &kinds, &ctx_eff)? {
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
                if !matches!(study, Study::Cavity { .. } | Study::Bloch { .. } | Study::Driven { .. }) {
                    return Err(Error::invalid(format!(
                        "a Floquet-periodic pair belongs to the cavity study, not to '{}'",
                        study.kind()
                    )));
                }
                let dest = ctx_eff.resolver.boundaries(&p.text("dest"))?;
                if dest.is_empty() {
                    return Err(Error::invalid(format!(
                        "the destination '{}' of a periodic pair selects nothing",
                        p.text("dest")
                    )));
                }
                periodic.push(Periodic { source: r.labels, dest, phase: p.number("phase") });
            }
            "pmc" => {
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
        symmetry,
        d_omega_eps,
    })
}
