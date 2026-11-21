// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared helpers of the ems tests: the mesh fixtures of `femcore::testkit`, a solve helper over the typed problem, the effective index of a
//! solution, and the conversion of an old-shape `ProblemSpec` (what the use cases build) to the typed problem.
#![allow(dead_code, unused_imports)]
use femcore::progress::Progress;
pub use femcore::testkit::*;
use nyx_ems::{Problem, Study};
use nyx_physics::{Discretisation, Mesh, Solution, SolveResult};

/// An EM-wave problem with the defaults of the old `ProblemSpec::new(Physics::Emw)`: order 2, a `Modes` study at 10 GHz.
pub fn emw() -> Problem {
    Problem::new(Study::Modes { freq: 10e9, n_target: 1.0, evanescent: false, n_modes: 4 }, 2)
}

/// Solves a typed problem on a fixture mesh (what `Physics::solve` does after `decode`); the stand-in of the old `solve::solve`.
pub fn solve(mesh: &Discretisation, p: &Problem) -> Result<SolveResult, String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    nyx_ems::run(&exchange, p, &Progress::new()).map_err(|e| e.to_string())
}

/// `solve` with the polarization analysis (`nyx_analysis::polarization`, a function of the result) next to the result.
pub fn try_solve_analysed(
    mesh: &Discretisation,
    p: &Problem,
) -> Result<(SolveResult, nyx_analysis::ModeAnalysis), String> {
    let exchange = Mesh::from_solver_mesh(mesh)?;
    let result = nyx_ems::run(&exchange, p, &Progress::new()).map_err(|e| e.to_string())?;
    let analysis = nyx_analysis::polarization(&exchange, &result).map_err(|e| e.to_string())?;
    Ok((result, analysis))
}

/// The effective index (re, im) of a mode solution.
pub fn neff(s: &Solution) -> (f64, f64) {
    let e = s.eigenvalue.expect("a mode with an effective index");
    (e[0], e[1])
}

// ---------------------------------------------------------------------------------------------------------------- mesher sizing fixtures

use femcore::adapt::drive::AdaptReport;
use femcore::adapt::ops3::AdaptOptions3;
use femcore::material::M6;
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::physics_size::{adapt_tets_for, optical_wavelength, speed_of_light, PhysicsSizing};
use femcore::meshing::sizing::label_names;
use femcore::meshing::{mesh_for_sizing, MeshPlan, Sizing};
use femcore::props::Materials;
use geom::model::FeatureList;

/// What the sizing tests describe of an EM problem: element order, materials by domain number and the study. The physics crate states its
/// own `Sizing` from a problem in the envelope (`Ems::sizing`); these tests start from numbers per domain, so they state it here the same way
/// (`n^2 = mu_r * eps_r`, wavelength `c / (f n)`, headroom 1.5 for an eigenfrequency search).
pub struct Spec {
    pub order: usize,
    pub materials: Materials,
    pub study: Study,
}

impl Spec {
    pub fn emw() -> Spec {
        Spec {
            order: 2,
            materials: Materials::default(),
            study: Study::Modes { freq: 10e9, n_target: 1.0, evanescent: false, n_modes: 4 },
        }
    }

    pub fn problem(&self) -> Problem {
        Problem {
            order: self.order,
            study: self.study.clone(),
            materials: self.materials.clone(),
            periodic: vec![],
            pmc: vec![],
            surface_impedance: vec![],
        }
    }

    /// The wavelength of every domain in `domains`, and the headroom of the study.
    fn wavelengths(&self, domains: &[i32]) -> (std::collections::BTreeMap<i32, f64>, f64) {
        let (freq, headroom) = match self.study {
            Study::Modes { freq, .. } => (freq, 1.0),
            Study::Cavity { f_search, .. } => (f_search, 1.5),
            _ => panic!("a sizing test needs a mode or cavity study"),
        };
        let lam0 = speed_of_light() / freq;
        let lam = domains
            .iter()
            .map(|&d| {
                let (eps, mu) =
                    (self.materials.get(d, "eps_r").unwrap_or(1.0), self.materials.get(d, "mu_r").unwrap_or(1.0));
                (d, optical_wavelength(&[[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, eps]], mu, lam0))
            })
            .collect();
        (lam, headroom)
    }
}

/// The 2D mesh the physics asks for: the regions of the model named in a `Sizing`, as `Physics::sizing` states it.
pub fn mesh_for_problem(
    list: &FeatureList,
    plan: &MeshPlan,
    spec: &Spec,
    sz: &PhysicsSizing,
) -> Result<(MeshedGeom, Option<AdaptReport>), String> {
    let ev = femcore::meshing::geomdrive::evaluate(list)?;
    let labels = ev.label_map().map_err(|e| e.to_string())?;
    let ids: Vec<i32> = labels.domains.iter().map(|d| d.id).collect();
    let (lam, headroom) = spec.wavelengths(&ids);
    let wavelengths = labels.domains.iter().map(|d| (d.name.clone(), lam[&d.id])).collect();
    mesh_for_sizing(list, plan, &Sizing { order: spec.order as u8, headroom, wavelengths }, sz)
}

/// Adapts a tetrahedral mesh (domains by number) to the waves of `spec`.
pub fn adapt_tets(
    mesh: &SimplexMesh,
    spec: &Spec,
    sz: &PhysicsSizing,
    ceiling: f64,
    grade: f64,
    opt: &AdaptOptions3,
) -> Result<(SimplexMesh, AdaptReport), String> {
    let (lam, headroom) = spec.wavelengths(&mesh.domain_labels());
    let sizing = Sizing::by_label(spec.order as u8, headroom, &lam);
    let pml = |d: i32| spec.materials.pml.contains_key(&d);
    adapt_tets_for(mesh, &sizing, &label_names, &pml, sz, ceiling, grade, opt)
}

/// Solves a `Spec` on a fixture mesh.
pub fn solve_spec(mesh: &Discretisation, spec: &Spec) -> Result<SolveResult, String> {
    solve(mesh, &spec.problem())
}

/// The problem of a `SoiWire` on its meshed cross-section: permittivity by region name, modes near the silicon index (what the SOI use case
/// built until phase 7).
pub fn wire_problem(w: &femcore::testkit::soi::SoiWire, g: &MeshedGeom) -> Problem {
    use femcore::testkit::soi::Cladding;
    let mut p = emw();
    let n_top = match w.cladding {
        Cladding::Oxide => w.n_sio2,
        Cladding::Air => 1.0,
    };
    for (regions, n) in [(["BOX", "BOX_near"], w.n_sio2), (["Si", "Si_slab"], w.n_si), (["clad", "clad_near"], n_top)] {
        for r in regions {
            if let Some(d) = g.domain(r) {
                p.materials.set(d, "eps_r", n * n);
                p.materials.set(d, "mu_r", 1.0);
            }
        }
    }
    p.study = Study::Modes {
        freq: speed_of_light() / w.wavelength,
        n_target: w.n_target,
        evanescent: false,
        n_modes: w.n_modes,
    };
    p
}

#[allow(dead_code)]
fn _unused(_: &M6) {}
