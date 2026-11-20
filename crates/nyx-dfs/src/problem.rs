// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `dfs`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (names and selectors
//! resolved to mesh numbers, parameters with their defaults), or built directly by code that already has mesh numbers (tests, the
//! legacy window bridge of the engine). The envelope carries scalar material properties only; anisotropic tensors are given in
//! `materials` directly until the material files exist (phase 6).

use crate::Dfs;
use femcore::props::{LinearSolver, Materials};
use nyx_physics::decode::{materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec};

#[derive(Clone, Debug, PartialEq)]
pub enum BcKind {
    Potential(f64),
    Ground,
    /// surface charge density, C/m^2
    SurfaceCharge(f64),
    Temperature(f64),
    /// inward heat flux, W/m^2
    Flux(f64),
    Convection {
        h: f64,
        text: f64,
    },
}

impl BcKind {
    pub fn name(&self) -> &'static str {
        match self {
            BcKind::Potential(_) => "Electric potential",
            BcKind::Ground => "Ground",
            BcKind::SurfaceCharge(_) => "Surface charge",
            BcKind::Temperature(_) => "Temperature",
            BcKind::Flux(_) => "Heat flux",
            BcKind::Convection { .. } => "Convective cooling",
        }
    }

    fn is_electrostatic(&self) -> bool {
        matches!(self, BcKind::Potential(_) | BcKind::Ground | BcKind::SurfaceCharge(_))
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
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Electrostatic => "electrostatic",
            Study::HeatStationary => "heat_stationary",
            Study::HeatTransient { .. } => "heat_transient",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Problem {
    /// element order, 1 or 2
    pub order: usize,
    pub study: Study,
    pub materials: Materials,
    pub bcs: Vec<Bc>,
    pub solver: LinearSolver,
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem { order, study, materials: Materials::default(), bcs: vec![], solver: LinearSolver::Auto }
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
                BcKind::Ground => vec![],
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
        if let Study::HeatTransient { t_end, .. } = &self.study {
            if !t_end.is_finite() {
                return Err("the study's end time is NaN or infinite".into());
            }
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
    let study = match sd.kind.as_str() {
        "electrostatic" => Study::Electrostatic,
        "heat_stationary" => Study::HeatStationary,
        "heat_transient" => Study::HeatTransient {
            t_end: sp.number("t_end"),
            n_out: sp.integer("n_out").max(0) as usize,
            t0: sp.text("t0"),
        },
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
            "temperature" => BcKind::Temperature(p.number("value")),
            "flux" => BcKind::Flux(p.number("value")),
            "convection" => BcKind::Convection { h: p.number("h"), text: p.number("t_ext") },
            other => return Err(Error::Unsupported(format!("dfs has no condition '{other}'"))),
        };
        let electro = matches!(study, Study::Electrostatic);
        if kind.is_electrostatic() != electro {
            return Err(Error::invalid(format!(
                "'{}' is not {} boundary condition",
                kind.name(),
                if electro { "an electrostatics" } else { "a heat-transfer" }
            )));
        }
        bcs.push(Bc { labels: r.labels, kind });
    }
    let params = dfs.params();
    let pp = Params::new(&spec.params, &params);
    let solver = match pp.text("solver").as_str() {
        "direct" => LinearSolver::Direct,
        "iterative" => LinearSolver::Iterative,
        _ => LinearSolver::Auto,
    };
    Ok(Problem {
        order: spec.order as usize,
        study,
        materials: materials_by_label(spec, ctx, &dfs.properties())?,
        bcs,
        solver,
    })
}
