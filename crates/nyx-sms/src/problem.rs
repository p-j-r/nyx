// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `sms`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (names and selectors
//! resolved to mesh numbers, parameters with their defaults, the imported temperature taken from `ctx.inputs`), or built directly by
//! code that already has mesh numbers (tests, the legacy window bridge of the engine).

use crate::Sms;
use femcore::props::{LinearSolver, Materials, Periodic, Pml, PmlCoord, PmlDir, PmlProfile, PmlStabilization};
use nyx_physics::decode::{materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Plane {
    Strain,
    Stress,
}

/// Piezoelectric coupling activation policy in solid mechanics (`sms`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PiezoMode {
    /// Enable coupled `(u, phi)` assembly automatically when any domain has `[piezoelectric]` data
    /// or when electrical boundary conditions (`ground`, `potential`, etc.) are present.
    #[default]
    Auto,
    /// Force coupled piezoelectric/electro-elastic assembly even with zero coupling tensor (`e = 0`).
    On,
    /// Solve pure mechanical elasticity, ignoring `[piezoelectric]` coupling tensors.
    Off,
}

/// Electrical role of a domain-level metallic `conductor` condition.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ConductorRole {
    /// Equipotential electrode whose potential floats with net free charge `charge` (0 C in eigenfrequency).
    #[default]
    Floating,
    /// Grounded electrode (`phi = 0 V`).
    Ground,
    /// Fixed-potential electrode (`phi = value`, static only when `value != 0`).
    Fixed,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BcKind {
    /// per-component constraint (x, y, z, in wafer axes) to the given displacement
    Fixed { fixed: [bool; 3], values: [f64; 3] },
    /// boundary load (force per area), in wafer axes
    Load([f64; 3]),
    /// absorbing traction for outgoing pressure and shear waves (eigenfrequency studies: complex frequencies, Q)
    LowReflecting,
    /// Grounded electrical boundary (`phi = 0 V`)
    Ground,
    /// Fixed electric potential (`phi = value` in V)
    Potential(f64),
    /// Surface free charge density (`sigma_s` in `C/m^2`, static only when non-zero)
    SurfaceCharge(f64),
    /// Equipotential floating electrode boundary with prescribed net charge `charge` (C, 0 in eigenfrequency)
    FloatingPotential { charge: f64 },
    /// Domain-level equipotential metallic body (skipped in electrical `K_{phi phi}` assembly)
    Conductor { role: ConductorRole, value: f64, charge: f64 },
}

impl BcKind {
    pub fn name(&self) -> &'static str {
        match self {
            BcKind::Fixed { .. } => "Fixed constraint",
            BcKind::Load(_) => "Boundary load",
            BcKind::LowReflecting => "Low-reflecting boundary",
            BcKind::Ground => "Ground",
            BcKind::Potential(_) => "Electric potential",
            BcKind::SurfaceCharge(_) => "Surface charge",
            BcKind::FloatingPotential { .. } => "Floating potential",
            BcKind::Conductor { .. } => "Conductor",
        }
    }

    /// True if this boundary condition applies to electrical potential `phi`.
    pub fn is_electrical(&self) -> bool {
        matches!(
            self,
            BcKind::Ground
                | BcKind::Potential(_)
                | BcKind::SurfaceCharge(_)
                | BcKind::FloatingPotential { .. }
                | BcKind::Conductor { .. }
        )
    }

    /// True if this condition selects volumetric domain labels rather than boundary facets.
    pub fn is_domain_condition(&self) -> bool {
        matches!(self, BcKind::Conductor { .. })
    }
}

/// A condition on boundary numbers of the mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct Bc {
    pub labels: Vec<i32>,
    pub kind: BcKind,
}

/// Real coordinate-stretching profile for static Infinite Elements (`sms::static`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum IeProfile {
    /// Locking-free Astley–Bettess–Zienkiewicz rational mapping `[R_in, R_out) -> [R_in, infinity)`
    /// with clamped `u = 0` outer wall (`u(infinity) = 0`).
    #[default]
    Rational,
    /// Bounded exponential stretching `u_tilde - u_0 = d * (chi^t - 1) / ln(chi)` paired with an
    /// outer 1st-order asymptotic elastic spring Robin boundary condition.
    Exponential,
}

/// Domain configuration for a static Infinite Element (`"infinite_element"` on `sms::static`).
#[derive(Clone, Debug, PartialEq)]
pub struct InfiniteElement {
    /// Underlying coordinate topology (`Cartesian`, `Cylindrical`, `Spherical`), active axes, and bounds.
    pub pml: Pml,
    /// Optional explicit pole distance `d_pole` (m) from the inner IE interface toward the source center;
    /// `None` auto-calibrates from `r_in` or distance from origin to the inner slab interface.
    pub pole_dist: Option<f64>,
    /// Real radial/normal stretching law (`Rational` or `Exponential`).
    pub profile: IeProfile,
    /// Exponential stretch ratio `chi = exp(beta) >= 1.0` (default `30.0`, used when `profile == Exponential`).
    pub decay: f64,
    /// Whether exterior outer termination facets are clamped `u = 0` (`true` by default for `Rational`,
    /// `false` by default for `Exponential` which applies a 1st-order spring Robin boundary condition).
    pub clamp_outer: bool,
}

impl Default for InfiniteElement {
    fn default() -> Self {
        Self { pml: Pml::default(), pole_dist: None, profile: IeProfile::Rational, decay: 30.0, clamp_outer: true }
    }
}

impl InfiniteElement {
    /// True if this Infinite Element domain is active.
    pub fn is_active(&self) -> bool {
        match self.pml.coord {
            PmlCoord::Cartesian => self.pml.x != PmlDir::Off || self.pml.y != PmlDir::Off || self.pml.z != PmlDir::Off,
            // Cylindrical & Spherical IEs can also infer `(r_in, thickness)` from the mesh domain when omitted (`0.0`).
            PmlCoord::Cylindrical { radial_on, axial, .. } => radial_on || axial != PmlDir::Off,
            PmlCoord::Spherical { .. } => true,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalSource {
    Uniform(f64),
    /// the temperature at the scalar dofs of the mesh and element order of the solve, exported by a heat-transfer study
    Imported(Vec<f64>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct ThermalExpansion {
    pub source: ThermalSource,
    pub t_ref: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Study {
    Static,
    /// eigenfrequencies; `f_search` > 0 finds modes near that frequency (Hz)
    Eigenfrequency {
        n_modes: usize,
        f_search: f64,
    },
    /// vibration modes of a body of revolution on its (r, z) half plane (mesh x = r); the displacement goes as `exp(-j m theta)`
    AxisymmetricEigenfrequency {
        n_modes: usize,
        f_search: f64,
        azimuthal_order: i64,
    },
    /// elastic waveguide modes on the 2D (x, z) cross-section; the displacement goes as `exp(-j q y)`
    Modes {
        n_modes: usize,
        f_search: f64,
        wavevector: f64,
    },
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Static => "static",
            Study::Eigenfrequency { .. } => "eigenfrequency",
            Study::AxisymmetricEigenfrequency { .. } => "axisymmetric_eigenfrequency",
            Study::Modes { .. } => "modes",
        }
    }

    /// Whether the study finds eigenfrequencies (`Eigenfrequency`, `AxisymmetricEigenfrequency` or `Modes`).
    pub fn is_eigen(&self) -> bool {
        matches!(self, Study::Eigenfrequency { .. } | Study::AxisymmetricEigenfrequency { .. } | Study::Modes { .. })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Problem {
    /// element order, 1 or 2
    pub order: usize,
    pub study: Study,
    pub materials: Materials,
    pub bcs: Vec<Bc>,
    pub plane: Plane,
    /// 2D out-of-plane thickness (m)
    pub thickness: f64,
    pub thermal: Option<ThermalExpansion>,
    pub solver: LinearSolver,
    /// domains that do not take part in this physics (for example air): they carry no stiffness or mass and the dofs they alone touch are
    /// held at zero
    pub inactive: Vec<i32>,
    /// Floquet-periodic boundary pairs (eigenfrequency studies)
    pub periodic: Vec<Periodic>,
    /// Static Infinite Element domains (`sms::static` open-boundary real coordinate stretching)
    pub infinite_elements: BTreeMap<i32, InfiniteElement>,
    /// Coupled piezoelectric activation mode (`Auto`, `On`, or `Off`)
    pub piezoelectric: PiezoMode,
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem {
            order,
            study,
            materials: Materials::default(),
            bcs: vec![],
            plane: Plane::Strain,
            thickness: 1.0,
            thermal: None,
            solver: LinearSolver::Auto,
            inactive: vec![],
            periodic: vec![],
            infinite_elements: BTreeMap::new(),
            piezoelectric: PiezoMode::Auto,
        }
    }

    /// True if any domain carries an active static Infinite Element.
    pub fn has_infinite_elements(&self) -> bool {
        self.infinite_elements.values().any(InfiniteElement::is_active)
    }

    /// True if any electrical boundary/domain condition (`ground`, `potential`, `surface_charge`,
    /// `floating_potential`, `conductor`) is present.
    pub fn has_electrical_bcs(&self) -> bool {
        self.bcs.iter().any(|bc| bc.kind.is_electrical())
    }

    /// Domain labels marked as equipotential `conductor` bodies (sorted, unique).
    pub fn conductor_domains(&self) -> Vec<i32> {
        let mut out: Vec<i32> = self
            .bcs
            .iter()
            .filter_map(|bc| match bc.kind {
                BcKind::Conductor { .. } => Some(bc.labels.iter().copied()),
                _ => None,
            })
            .flatten()
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// True when the coupled `(u, phi)` piezoelectric/electro-elastic solver should handle this problem.
    pub fn uses_piezoelectric(&self) -> bool {
        match self.piezoelectric {
            PiezoMode::Off => false,
            PiezoMode::On => true,
            PiezoMode::Auto => self.has_electrical_bcs(),
        }
    }

    /// Every number the user can type must be finite: NaN or infinity does not stop a factorisation or a Krylov iteration, it flows into
    /// every field. The mesh is checked by `nyx_physics::check_mesh`.
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
                BcKind::Fixed { values, .. } => values.to_vec(),
                BcKind::Load(f) => f.to_vec(),
                BcKind::LowReflecting | BcKind::Ground => vec![],
                BcKind::Potential(v) | BcKind::SurfaceCharge(v) | BcKind::FloatingPotential { charge: v } => vec![*v],
                BcKind::Conductor { value, charge, .. } => vec![*value, *charge],
            };
            if vals.iter().any(|v| !v.is_finite()) {
                return Err(format!(
                    "boundary condition {} ({} on boundaries {:?}) has a NaN or infinite value",
                    i + 1,
                    bc.kind.name(),
                    bc.labels
                ));
            }
        }
        if let Study::Eigenfrequency { f_search, .. }
        | Study::AxisymmetricEigenfrequency { f_search, .. }
        | Study::Modes { f_search, .. } = &self.study
        {
            if !f_search.is_finite() {
                return Err("the study's search frequency is NaN or infinite".into());
            }
        }
        if let Study::Modes { wavevector, .. } = &self.study {
            if !wavevector.is_finite() {
                return Err("the study's longitudinal wavevector is NaN or infinite".into());
            }
        }
        if !self.thickness.is_finite() {
            return Err("the 2D thickness is NaN or infinite".into());
        }
        // everything else (tensors, PML/IE data, thermal expansion, Floquet phases): scan the printed form for the words NaN and inf
        let text =
            format!("{:?} {:?} {:?} {:?} {:?}", m.aniso, m.pml, self.infinite_elements, self.thermal, self.periodic);
        if let Some(at) = femcore::validate::non_finite_token(&text) {
            let from = text[..at].char_indices().rev().nth(70).map_or(0, |(i, _)| i);
            return Err(format!(
                "the problem description contains a NaN or infinite number, near: ...{}",
                &text[from..(at + 5).min(text.len())]
            ));
        }
        Ok(())
    }

    /// Refuses what a wafer-frame 2D problem cannot say: there is no y axis in the plane of a cross-section.
    pub fn check_wafer(&self, dim: usize) -> Result<(), String> {
        if dim == 2 {
            if let Some((d, _)) =
                self.materials.pml.iter().find(|(_, p)| matches!(p.coord, PmlCoord::Cartesian) && p.y != PmlDir::Off)
            {
                return Err(format!(
                    "the perfectly matched layer of domain {d} absorbs along y, which is out of the plane of a 2D cross-section (x-z): use z for the direction of the surface normal"
                ));
            }
            if let Some((d, _)) = self
                .infinite_elements
                .iter()
                .find(|(_, ie)| matches!(ie.pml.coord, PmlCoord::Cartesian) && ie.pml.y != PmlDir::Off)
            {
                return Err(format!(
                    "the infinite_element of domain {d} stretches along y, which is out of the plane of a 2D cross-section (x-z): use z for the direction of the surface normal"
                ));
            }
        }
        Ok(())
    }
}

impl Problem {
    /// Rejects `infinite_element` domains on dynamic eigenfrequency studies.
    pub fn check_study_open_domains(&self) -> Result<(), Error> {
        if self.study.is_eigen() && !self.infinite_elements.is_empty() {
            return Err(Error::invalid(
                "an infinite_element domain is for the static study only (real coordinate stretching reflects propagating waves): use 'pml' for eigenfrequency studies",
            ));
        }
        Ok(())
    }

    /// Validates electrical/piezoelectric settings and boundary/domain conflicts.
    pub fn check_electrical(&self) -> Result<(), Error> {
        if self.piezoelectric == PiezoMode::Off && self.has_electrical_bcs() {
            return Err(Error::invalid(
                "electrical boundary conditions ('ground', 'potential', 'surface_charge', 'floating_potential', 'conductor') cannot be used when piezoelectric = 'off'",
            ));
        }
        let eigen = self.study.is_eigen();
        let mut seen_conductors: BTreeMap<i32, (ConductorRole, f64, f64)> = BTreeMap::new();
        for bc in &self.bcs {
            match &bc.kind {
                BcKind::Potential(v) if eigen && *v != 0.0 => {
                    return Err(Error::invalid(
                        "an eigenfrequency study takes zero electrical sources only: fixed potential must be 0 V (or use 'ground')",
                    ));
                }
                BcKind::SurfaceCharge(q) if eigen && *q != 0.0 => {
                    return Err(Error::invalid(
                        "a non-zero surface_charge boundary condition belongs to the static study",
                    ));
                }
                BcKind::FloatingPotential { charge } if eigen && *charge != 0.0 => {
                    return Err(Error::invalid(
                        "a floating_potential with non-zero charge belongs to the static study",
                    ));
                }
                BcKind::Conductor { role, value, charge } => {
                    if eigen {
                        match role {
                            ConductorRole::Fixed if *value != 0.0 => {
                                return Err(Error::invalid(
                                    "an eigenfrequency study takes zero electrical sources only: conductor fixed potential must be 0 V",
                                ));
                            }
                            ConductorRole::Floating if *charge != 0.0 => {
                                return Err(Error::invalid(
                                    "a conductor with non-zero charge belongs to the static study",
                                ));
                            }
                            _ => {}
                        }
                    }
                    for &dom in &bc.labels {
                        if self.materials.pml.contains_key(&dom) {
                            return Err(Error::invalid(format!(
                                "domain {dom} cannot be both an equipotential 'conductor' and a 'pml'"
                            )));
                        }
                        if self.infinite_elements.contains_key(&dom) {
                            return Err(Error::invalid(format!(
                                "domain {dom} cannot be both an equipotential 'conductor' and an 'infinite_element'"
                            )));
                        }
                        if let Some(&(r0, v0, c0)) = seen_conductors.get(&dom) {
                            let same = r0 == *role && (v0 - value).abs() <= 1e-12 && (c0 - charge).abs() <= 1e-12;
                            if !same {
                                return Err(Error::invalid(format!(
                                    "conflicting 'conductor' conditions assigned to domain {dom}"
                                )));
                            }
                        } else {
                            seen_conductors.insert(dom, (*role, *value, *charge));
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// The rules of the axisymmetric study that do not need the mesh (ADR 68): no boundary load, no Floquet pair, no low-reflecting boundary
    /// (its damping matrix has no `r` weight, rough spot D1), only zero constraints, and a layer cannot absorb towards the axis. A problem of
    /// another study passes. The mesh rule (the half plane lies in `x >= 0`) is checked by `cross_section::solve_axisymmetric`.
    pub fn check_axisymmetric(&self) -> Result<(), Error> {
        self.check_study_open_domains()?;
        if !matches!(self.study, Study::AxisymmetricEigenfrequency { .. }) {
            return Ok(());
        }
        if !self.periodic.is_empty() {
            return Err(Error::invalid(
                "a Floquet-periodic pair has no meaning in an axisymmetric model (no translation to repeat)",
            ));
        }
        for bc in &self.bcs {
            match &bc.kind {
                BcKind::Load(_) => return Err(Error::invalid("a boundary load belongs to the static study")),
                BcKind::LowReflecting => {
                    return Err(Error::Unsupported(
                        "a low-reflecting boundary is not available in the axisymmetric study yet (its damping matrix has no r weight, rough spot D1): use a perfectly matched layer".into(),
                    ))
                }
                BcKind::Fixed { fixed, values } => {
                    if (0..3).any(|k| fixed[k] && values[k] != 0.0) {
                        return Err(Error::invalid(
                            "an eigenfrequency study takes zero constraints only: a fixed displacement must be 0",
                        ));
                    }
                }
                BcKind::Ground
                | BcKind::Potential(_)
                | BcKind::SurfaceCharge(_)
                | BcKind::FloatingPotential { .. }
                | BcKind::Conductor { .. } => {}
            }
        }
        for (d, p) in &self.materials.pml {
            match p.coord {
                PmlCoord::Cartesian => {
                    if p.x == PmlDir::Minus {
                        return Err(Error::invalid(format!(
                            "the perfectly matched layer of domain {d} absorbs towards the axis (x = minus): a layer cannot absorb towards the axis, use plus"
                        )));
                    }
                }
                PmlCoord::Spherical { center, r_in, thickness } => {
                    let scale = r_in.abs().max(thickness.abs()).max(1.0);
                    if center[0].abs() > 1e-9 * scale {
                        return Err(Error::invalid(format!(
                            "axisymmetric spherical PML on domain {d} must have its center on the symmetry axis r = 0 (got r_center = {})",
                            center[0]
                        )));
                    }
                }
                PmlCoord::Cylindrical { center, axis, r_in, r_thickness, .. } => {
                    let scale = r_in.abs().max(r_thickness.abs()).max(1.0);
                    if center[0].abs() > 1e-9 * scale {
                        return Err(Error::invalid(format!(
                            "axisymmetric cylindrical PML on domain {d} must be centered on the symmetry axis r = 0 (got r_center = {})",
                            center[0]
                        )));
                    }
                    let ax_xy = axis[0].hypot(axis[1]);
                    let ax_n = (ax_xy * ax_xy + axis[2] * axis[2]).sqrt();
                    if ax_n > 1e-30 && ax_xy > 1e-9 * ax_n {
                        return Err(Error::invalid(format!(
                            "axisymmetric cylindrical PML on domain {d} must align its axis with z = [0, 0, +-1] (got {axis:?})"
                        )));
                    }
                }
            }
        }
        Ok(())
    }

    /// The rules of the waveguide modes study: no boundary load, and only zero constraints. Periodic lateral pairs and
    /// absorbing layers in the cross-section plane (x and z; y is refused by `check_wafer(2)`) are allowed.
    pub fn check_modes(&self) -> Result<(), Error> {
        self.check_study_open_domains()?;
        self.check_electrical()?;
        if !matches!(self.study, Study::Modes { .. }) {
            return Ok(());
        }
        for bc in &self.bcs {
            match &bc.kind {
                BcKind::Load(_) => return Err(Error::invalid("a boundary load belongs to the static study")),
                BcKind::LowReflecting => {}
                BcKind::Fixed { fixed, values } => {
                    if (0..3).any(|k| fixed[k] && values[k] != 0.0) {
                        return Err(Error::invalid(
                            "an eigenfrequency study takes zero constraints only: a fixed displacement must be 0",
                        ));
                    }
                }
                BcKind::Ground
                | BcKind::Potential(_)
                | BcKind::SurfaceCharge(_)
                | BcKind::FloatingPotential { .. }
                | BcKind::Conductor { .. } => {}
            }
        }
        Ok(())
    }
}

fn pml_dir(s: &str) -> PmlDir {
    match s {
        "plus" => PmlDir::Plus,
        "minus" => PmlDir::Minus,
        _ => PmlDir::Off,
    }
}

fn decode_coord_and_bounds(p: &Params, what: &str) -> Result<(PmlCoord, Option<([f64; 3], [f64; 3])>), Error> {
    let blo = p.optional_numbers("bounds_lo");
    let bhi = p.optional_numbers("bounds_hi");
    let bounds = match (blo, bhi) {
        (None, None) => None,
        (Some(lo), Some(hi)) => {
            if lo.len() != 3 || hi.len() != 3 {
                return Err(Error::invalid(format!(
                    "{what}: 'bounds_lo' and 'bounds_hi' must each have 3 wafer coordinates [x, y, z]"
                )));
            }
            for k in 0..3 {
                if lo[k] > hi[k] {
                    return Err(Error::invalid(format!(
                        "{what}: bounds_lo[{k}] ({}) exceeds bounds_hi[{k}] ({})",
                        lo[k], hi[k]
                    )));
                }
            }
            Some(([lo[0], lo[1], lo[2]], [hi[0], hi[1], hi[2]]))
        }
        _ => {
            return Err(Error::invalid(format!("{what}: both 'bounds_lo' and 'bounds_hi' must be given together")));
        }
    };

    let c_vec = p.vector("center");
    let center = [c_vec[0], c_vec[1], c_vec[2]];
    let a_vec = p.vector("axis");
    let axis = [a_vec[0], a_vec[1], a_vec[2]];
    let r_in = p.number("r_in");
    let thickness = p.number("thickness");

    let coord = match p.text("coord").as_str() {
        "cartesian" => PmlCoord::Cartesian,
        "cylindrical" => PmlCoord::Cylindrical {
            center,
            axis,
            r_in,
            r_thickness: thickness,
            radial_on: p.boolean("radial_on"),
            axial: pml_dir(&p.text("axial")),
            z_in: p.number("z_in"),
            z_thickness: p.number("z_thickness"),
        },
        "spherical" => PmlCoord::Spherical { center, r_in, thickness },
        other => return Err(Error::invalid(format!("{what}: unsupported coordinate system '{other}'"))),
    };
    Ok((coord, bounds))
}

/// Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and
/// selectors with `ctx.resolver`, takes the imported temperature from `ctx.inputs` and refuses combinations that need the study (a
/// perfectly matched layer in a static study, or an infinite element in a dynamic study).
pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error> {
    let sms = Sms;
    let studies = sms.studies();
    let sd = studies
        .iter()
        .find(|s| s.kind == spec.study.kind)
        .ok_or_else(|| Error::Unsupported(format!("sms has no study '{}'", spec.study.kind)))?;
    let sp = Params::new(&spec.study.params, &sd.params);
    let study = match sd.kind.as_str() {
        "static" => Study::Static,
        "eigenfrequency" => {
            Study::Eigenfrequency { n_modes: sp.integer("n_modes").max(0) as usize, f_search: sp.number("f_search") }
        }
        "axisymmetric_eigenfrequency" => Study::AxisymmetricEigenfrequency {
            n_modes: sp.integer("n_modes").max(0) as usize,
            f_search: sp.number("f_search"),
            azimuthal_order: sp.integer("azimuthal_order"),
        },
        "modes" => Study::Modes {
            n_modes: sp.integer("n_modes").max(0) as usize,
            f_search: sp.number("f_search"),
            wavevector: sp.number("wavevector"),
        },
        other => return Err(Error::Unsupported(format!("sms has no study '{other}'"))),
    };
    let eigen = study.is_eigen();
    let params = sms.params();
    let pp = Params::new(&spec.params, &params);
    let plane = if pp.text("plane") == "stress" { Plane::Stress } else { Plane::Strain };
    let thermal = match pp.text("thermal").as_str() {
        "uniform" => Some(ThermalExpansion {
            source: ThermalSource::Uniform(pp.number("temperature")),
            t_ref: pp.number("t_ref"),
        }),
        "imported" => {
            let t = ctx.inputs.get("temperature").ok_or_else(|| {
                Error::invalid("thermal expansion 'imported' needs a temperature exported by another physics (dfs.T), none was given")
            })?;
            Some(ThermalExpansion { source: ThermalSource::Imported(t.clone()), t_ref: pp.number("t_ref") })
        }
        _ => None,
    };
    let mut inactive = vec![];
    for region in pp.names("inactive") {
        inactive.extend(ctx.resolver.domains(&region)?);
    }
    inactive.sort_unstable();
    inactive.dedup();
    let mut materials = materials_by_label(spec, ctx, &sms.properties())?;
    let kinds = sms.bc_kinds();
    let mut bcs = vec![];
    let mut periodic = vec![];
    let mut infinite_elements = BTreeMap::new();
    for r in resolve_conditions(&spec.bcs, &kinds, ctx)? {
        let d = kinds.iter().find(|k| k.kind == r.kind).expect("resolved against these descriptors");
        let p = Params::new(&r.params, &d.params);
        let only_eigen = |what: &str| -> Result<(), Error> {
            if eigen {
                Ok(())
            } else {
                Err(Error::invalid(format!(
                    "{what} is used with the eigenfrequency study (it gives complex frequencies and Q)"
                )))
            }
        };
        match r.kind.as_str() {
            "fixed" => {
                let mut fixed = [false; 3];
                let mut values = [0.0; 3];
                for (k, name) in ["ux", "uy", "uz"].iter().enumerate() {
                    if p.has(name) {
                        fixed[k] = true;
                        values[k] = p.number(name);
                    }
                }
                bcs.push(Bc { labels: r.labels, kind: BcKind::Fixed { fixed, values } });
            }
            "load" => {
                if eigen {
                    return Err(Error::invalid("a boundary load belongs to the static study"));
                }
                let t = p.vector("traction");
                bcs.push(Bc { labels: r.labels, kind: BcKind::Load([t[0], t[1], t[2]]) });
            }
            "low_reflecting" => {
                only_eigen("a low-reflecting boundary")?;
                bcs.push(Bc { labels: r.labels, kind: BcKind::LowReflecting });
            }
            "pml" => {
                only_eigen("a perfectly matched layer")?;
                let (coord, bounds) = decode_coord_and_bounds(&p, "condition 'pml'")?;
                let profile = match p.text("profile").as_str() {
                    "rational" => PmlProfile::Rational { curvature: p.number("curvature") },
                    _ => PmlProfile::KappaPolynomial { kappa_max: p.number("kappa_max") },
                };
                let stabilization = match p.text("stabilization").as_str() {
                    "off" => PmlStabilization::Off,
                    "mpml" => PmlStabilization::Mpml { ratio: p.number("mpml_ratio") },
                    "adiabatic" => PmlStabilization::Adiabatic { eta_max: p.number("adiabatic_eta") },
                    _ => PmlStabilization::Auto,
                };
                let pml = Pml {
                    x: pml_dir(&p.text("x")),
                    y: pml_dir(&p.text("y")),
                    z: pml_dir(&p.text("z")),
                    strength: p.number("strength"),
                    order: p.number("order"),
                    re: p.number("re"),
                    coord,
                    bounds,
                    profile,
                    auto_scale: p.boolean("auto_scale"),
                    target_reflection: p.number("target_reflection"),
                    stabilization,
                    clamp_outer: p.boolean("clamp_outer"),
                };
                for l in r.labels {
                    materials.pml.insert(l, pml.clone());
                }
            }
            "infinite_element" => {
                if eigen {
                    return Err(Error::invalid(
                        "an infinite_element domain is for the static study only (real coordinate stretching reflects propagating waves): use 'pml' for eigenfrequency studies",
                    ));
                }
                let (coord, bounds) = decode_coord_and_bounds(&p, "condition 'infinite_element'")?;
                let profile = match p.text("profile").as_str() {
                    "exponential" => IeProfile::Exponential,
                    _ => IeProfile::Rational,
                };
                let pole_dist = p.optional_number("pole_dist").filter(|&v| v > 0.0);
                let clamp_outer =
                    if p.has("clamp_outer") { p.boolean("clamp_outer") } else { profile == IeProfile::Rational };
                let pml = Pml {
                    x: pml_dir(&p.text("x")),
                    y: pml_dir(&p.text("y")),
                    z: pml_dir(&p.text("z")),
                    coord,
                    bounds,
                    clamp_outer,
                    ..Default::default()
                };
                let ie = InfiniteElement { pml, pole_dist, profile, decay: p.number("decay"), clamp_outer };
                if !ie.is_active() {
                    return Err(Error::invalid(
                        "condition 'infinite_element' has no active stretching direction (set x/y/z for cartesian, or coord = cylindrical / spherical)",
                    ));
                }
                for l in r.labels {
                    infinite_elements.insert(l, ie.clone());
                }
            }
            "periodic" => {
                only_eigen("a Floquet-periodic pair")?;
                let dest = ctx.resolver.boundaries(&p.text("dest"))?;
                if dest.is_empty() {
                    return Err(Error::invalid(format!(
                        "periodic pair: destination '{}' selects nothing",
                        p.text("dest")
                    )));
                }
                periodic.push(Periodic { source: r.labels, dest, phase: p.number("phase") });
            }
            "ground" => bcs.push(Bc { labels: r.labels, kind: BcKind::Ground }),
            "potential" => bcs.push(Bc { labels: r.labels, kind: BcKind::Potential(p.number("value")) }),
            "surface_charge" => bcs.push(Bc { labels: r.labels, kind: BcKind::SurfaceCharge(p.number("density")) }),
            "floating_potential" => {
                bcs.push(Bc { labels: r.labels, kind: BcKind::FloatingPotential { charge: p.number("charge") } })
            }
            "conductor" => {
                let role = match p.text("role").as_str() {
                    "ground" => ConductorRole::Ground,
                    "fixed" => ConductorRole::Fixed,
                    _ => ConductorRole::Floating,
                };
                bcs.push(Bc {
                    labels: r.labels,
                    kind: BcKind::Conductor { role, value: p.number("value"), charge: p.number("charge") },
                });
            }
            other => return Err(Error::Unsupported(format!("sms has no condition '{other}'"))),
        }
    }
    let piezoelectric = match pp.text("piezoelectric").as_str() {
        "on" => PiezoMode::On,
        "off" => PiezoMode::Off,
        _ => PiezoMode::Auto,
    };
    let solver = match pp.text("solver").as_str() {
        "direct" => LinearSolver::Direct,
        "iterative" => LinearSolver::Iterative,
        _ => LinearSolver::Auto,
    };
    let problem = Problem {
        order: spec.order as usize,
        study,
        materials,
        bcs,
        plane,
        thickness: pp.number("thickness"),
        thermal,
        solver,
        inactive,
        periodic,
        infinite_elements,
        piezoelectric,
    };
    problem.check_study_open_domains()?;
    problem.check_electrical()?;
    problem.check_axisymmetric()?;
    problem.check_modes()?;
    Ok(problem)
}
