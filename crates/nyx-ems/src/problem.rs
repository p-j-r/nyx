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

use crate::Ems;
use femcore::props::{Materials, Periodic, Pml, PmlDir};
use nyx_physics::decode::{materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec};

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
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Modes { .. } => "modes",
            Study::AxisymmetricModes { .. } => "axisymmetric_modes",
            Study::Cavity { .. } => "cavity",
        }
    }
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
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem { order, study, materials: Materials::default(), periodic: vec![], pmc: vec![] }
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
        };
        if let Some(what) = study_bad {
            return Err(format!("the study's {what} is NaN or infinite"));
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
        other => return Err(Error::Unsupported(format!("ems has no study '{other}'"))),
    };
    let mut materials = materials_by_label(spec, ctx, &ems.properties())?;
    let mut periodic = vec![];
    let mut pmc: Vec<i32> = vec![];
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
                if !matches!(study, Study::Cavity { .. }) {
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
                if !matches!(study, Study::Cavity { .. }) {
                    return Err(Error::invalid(format!(
                        "a magnetic wall belongs to the cavity study, not to '{}'",
                        study.kind()
                    )));
                }
                pmc.extend(r.labels);
            }
            other => return Err(Error::Unsupported(format!("ems has no condition '{other}'"))),
        }
    }
    pmc.sort_unstable();
    pmc.dedup();
    Ok(Problem { order: spec.order as usize, study, materials, periodic, pmc })
}
