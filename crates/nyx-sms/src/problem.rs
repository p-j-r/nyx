// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `sms`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (names and selectors
//! resolved to mesh numbers, parameters with their defaults, the imported temperature taken from `ctx.inputs`), or built directly by
//! code that already has mesh numbers (tests, the legacy window bridge of the engine). The envelope carries scalar material properties
//! only; anisotropic stiffness and expansion tensors are given in `materials.aniso` directly until the material files exist (phase 6).

use crate::Sms;
use femcore::props::{LinearSolver, Materials, Periodic, Pml, PmlDir};
use nyx_physics::decode::{materials_by_label, resolve_conditions, Params};
use nyx_physics::{Ctx, Error, Physics, ProblemSpec};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Plane {
    Strain,
    Stress,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BcKind {
    /// per-component constraint (x, y, z, in wafer axes) to the given displacement
    Fixed { fixed: [bool; 3], values: [f64; 3] },
    /// boundary load (force per area), in wafer axes
    Load([f64; 3]),
    /// absorbing traction for outgoing pressure and shear waves (eigenfrequency studies: complex frequencies, Q)
    LowReflecting,
}

impl BcKind {
    pub fn name(&self) -> &'static str {
        match self {
            BcKind::Fixed { .. } => "Fixed constraint",
            BcKind::Load(_) => "Boundary load",
            BcKind::LowReflecting => "Low-reflecting boundary",
        }
    }
}

/// A condition on boundary numbers of the mesh.
#[derive(Clone, Debug, PartialEq)]
pub struct Bc {
    pub labels: Vec<i32>,
    pub kind: BcKind,
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
                BcKind::LowReflecting => vec![],
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
        // everything else (tensors, PML data, thermal expansion, Floquet phases): scan the printed form for the words NaN and inf
        let text = format!("{:?} {:?} {:?} {:?}", m.aniso, m.pml, self.thermal, self.periodic);
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
            if let Some((d, _)) = self.materials.pml.iter().find(|(_, p)| p.y != PmlDir::Off) {
                return Err(format!(
                    "the perfectly matched layer of domain {d} absorbs along y, which is out of the plane of a 2D cross-section (x-z): use z for the direction of the surface normal"
                ));
            }
        }
        Ok(())
    }
}

impl Problem {
    /// The rules of the axisymmetric study that do not need the mesh (ADR 68): no boundary load, no Floquet pair, no low-reflecting boundary
    /// (its damping matrix has no `r` weight, rough spot D1), only zero constraints, and a layer cannot absorb towards the axis. A problem of
    /// another study passes. The mesh rule (the half plane lies in `x >= 0`) is checked by `cross_section::solve_axisymmetric`.
    pub fn check_axisymmetric(&self) -> Result<(), Error> {
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
            }
        }
        if let Some((d, _)) = self.materials.pml.iter().find(|(_, p)| p.x == PmlDir::Minus) {
            return Err(Error::invalid(format!(
                "the perfectly matched layer of domain {d} absorbs towards the axis (x = minus): a layer cannot absorb towards the axis, use plus"
            )));
        }
        Ok(())
    }

    /// The rules of the waveguide modes study: no boundary load, and only zero constraints. Periodic lateral pairs and
    /// absorbing layers in the cross-section plane (x and z; y is refused by `check_wafer(2)`) are allowed.
    pub fn check_modes(&self) -> Result<(), Error> {
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

/// Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and
/// selectors with `ctx.resolver`, takes the imported temperature from `ctx.inputs` and refuses combinations that need the study (a
/// perfectly matched layer in a static study).
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
                let pml = Pml {
                    x: pml_dir(&p.text("x")),
                    y: pml_dir(&p.text("y")),
                    z: pml_dir(&p.text("z")),
                    strength: p.number("strength"),
                    order: p.number("order"),
                    re: p.number("re"),
                };
                for l in r.labels {
                    materials.pml.insert(l, pml.clone());
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
            other => return Err(Error::Unsupported(format!("sms has no condition '{other}'"))),
        }
    }
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
    };
    problem.check_axisymmetric()?;
    problem.check_modes()?;
    Ok(problem)
}
