//! nyx-dfs: diffusion physics (`dfs`, ADR 12): electrostatics and heat transfer, stationary and transient.
//!
//! The crate holds the typed problem (`problem`, decoded from the envelope by `problem::decode`), the kernels in the three stages of ADR 9
//! (`electro`: `assemble_electrostatics`, `solve_potential`, `electrostatics_fields`; `heat`: `assemble_heat`, `solve_heat_states`,
//! `heat_fields`) and the `Physics` implementation below. Moved here from `femcore::dfs` in restructure phase 4 (ADR 41).

// The kernels index several arrays with one loop variable, the form that matches the formulas (as in femcore).
#![allow(clippy::needless_range_loop)]

mod common;
pub mod electro;
pub mod heat;
pub mod problem;

pub use electro::{
    assemble_electrostatics, electrostatics_fields, solve_electrostatics, solve_potential, ElectroAssembly, EPS0,
};
pub use heat::{assemble_heat, heat_fields, solve_heat, solve_heat_states, HeatStates, HeatSystem};
pub use problem::{decode, Bc, BcKind, Problem, Study};

use femcore::progress::{self, Progress};
use nyx_physics::*;
use std::sync::Arc;

pub struct Dfs;

const WAFER: [&str; 3] = ["x", "y", "z"];

impl Physics for Dfs {
    fn id(&self) -> &'static str {
        "dfs"
    }

    fn name(&self) -> &'static str {
        "Diffusion (electrostatics, heat transfer)"
    }

    fn studies(&self) -> Vec<StudyDesc> {
        vec![
            StudyDesc {
                kind: "electrostatic".into(),
                name: "Electrostatics".into(),
                doc: "Potential and field of fixed potentials and charges: div(eps grad V) = -rho_v.".into(),
                dims: Dims::Both,
                params: vec![],
            },
            StudyDesc {
                kind: "heat_stationary".into(),
                name: "Heat transfer, stationary".into(),
                doc: "Temperature for fixed temperatures, fluxes, convection and sources: div(k grad T) + Q = 0."
                    .into(),
                dims: Dims::Both,
                params: vec![],
            },
            StudyDesc {
                kind: "heat_transient".into(),
                name: "Heat transfer, transient".into(),
                doc: "Temperature in time from an initial field: rho Cp dT/dt = div(k grad T) + Q.".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::number("t_end", "End time", "s", 1.0, "")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("n_out", "Output times", 10, 1, "number of stored time steps"),
                    ParamDesc::text(
                        "t0",
                        "Initial temperature",
                        "293.15",
                        "an expression in x, y, z (m), giving kelvin",
                    ),
                ],
            },
        ]
    }

    fn bc_kinds(&self) -> Vec<BcDesc> {
        let b = |kind: &str, name: &str, doc: &str, params: Vec<ParamDesc>| BcDesc {
            kind: kind.into(),
            name: name.into(),
            doc: doc.into(),
            applies: Applies::Boundary,
            params,
        };
        vec![
            b(
                "potential",
                "Electric potential",
                "Fixed potential (electrostatics).",
                vec![ParamDesc::number("value", "Potential", "V", 1.0, "")],
            ),
            b("ground", "Ground", "Zero potential (electrostatics).", vec![]),
            b(
                "surface_charge",
                "Surface charge",
                "Surface charge density (electrostatics).",
                vec![ParamDesc::number("density", "Surface charge density", "C/m^2", 0.0, "")],
            ),
            b(
                "temperature",
                "Temperature",
                "Fixed temperature (heat transfer).",
                vec![ParamDesc::number("value", "Temperature", "K", 293.15, "")],
            ),
            b(
                "flux",
                "Heat flux",
                "Inward heat flux (heat transfer).",
                vec![ParamDesc::number("value", "Heat flux", "W/m^2", 0.0, "positive into the body")],
            ),
            b(
                "convection",
                "Convective cooling",
                "Newton cooling to an outside temperature (heat transfer).",
                vec![
                    ParamDesc::number("h", "Heat transfer coefficient", "W/(m^2 K)", 10.0, ""),
                    ParamDesc::number("t_ext", "Outside temperature", "K", 293.15, ""),
                ],
            ),
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        vec![
            PropDesc::new("eps_r", "Relative permittivity", "", Tensor, &["electrostatic"], None),
            PropDesc::new("rho_v", "Charge density", "C/m^3", Scalar, &["electrostatic"], Some(0.0)),
            PropDesc::new("k", "Thermal conductivity", "W/(m K)", Tensor, &["heat_stationary", "heat_transient"], None),
            PropDesc::new("rho", "Density", "kg/m^3", Scalar, &["heat_transient"], None),
            PropDesc::new("Cp", "Specific heat", "J/(kg K)", Scalar, &["heat_transient"], None),
            PropDesc::new("Q", "Heat source", "W/m^3", Scalar, &["heat_stationary", "heat_transient"], Some(0.0)),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![ParamDesc::choice(
            "solver",
            "Linear solver",
            &["auto", "direct", "iterative"],
            "auto",
            "auto: sparse direct for small systems, multigrid-preconditioned CG for large 3D ones; direct: Cholesky; iterative: AMG-CG",
        )]
    }

    fn fields(&self, study: &str) -> Vec<FieldDesc> {
        match study {
            "electrostatic" => vec![
                FieldDesc::new("V", "V", &[], "electric potential"),
                FieldDesc::new("E", "V/m", &WAFER, "electric field (2D: x and z)"),
            ],
            "heat_stationary" | "heat_transient" => vec![
                FieldDesc::new("T", "K", &[], "temperature"),
                FieldDesc::new("q", "W/m^2", &WAFER, "heat flux (2D: x and z)"),
            ],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![
            FieldLink::new("dfs.T", "temperature"),
            FieldLink::new("dfs.V", "potential"),
            FieldLink::new("dfs.E", "electric_field"),
        ]
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }

    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error> {
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("dfs supports element order 1 and 2, not {}", spec.order)));
        }
        if let Some(i) = mesh.points.iter().position(|v| !v.is_finite()) {
            return Err(Error::NonFinite(format!(
                "the mesh has a non-finite coordinate (point {}, axis {}): NaN or infinity in the geometry",
                i / mesh.dim.max(1),
                i % mesh.dim.max(1)
            )));
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec) -> Sizing {
        // diffusion has no wavelength to resolve; geometry rules size the mesh
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }

    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        validate_spec(self, spec)?;
        self.validate(spec, mesh)?;
        let problem = decode(spec, ctx)?;
        run(mesh, &problem, &ctx.progress)
    }
}

/// Solves a typed problem: checks the numbers, runs the kernel of the study under `progress` (stage log, cancellation) and appends the
/// stage log to the result. `Physics::solve` is this after `decode`; code that already has mesh numbers calls it directly.
pub fn run(mesh: &Mesh, problem: &Problem, progress: &Arc<Progress>) -> Result<SolveResult, Error> {
    let _installed = progress::install(progress);
    let first_line = progress.log().len();
    if !(1..=2).contains(&problem.order) {
        return Err(Error::invalid("element order must be 1 or 2"));
    }
    if progress.is_cancelled() {
        return Err(Error::Cancelled);
    }
    if let Some(i) = mesh.points.iter().position(|v| !v.is_finite()) {
        return Err(Error::NonFinite(format!(
            "the mesh has a non-finite coordinate (point {}, axis {}): NaN or infinity in the geometry",
            i / mesh.dim.max(1),
            i % mesh.dim.max(1)
        )));
    }
    problem.check_finite().map_err(Error::NonFinite)?;
    let mut result = {
        let _s = progress::stage(format!("solving diffusion on a {}D mesh of {} cells", mesh.dim, mesh.n_cells()));
        let d = discretise(mesh)?;
        match problem.study {
            Study::Electrostatic => solve_electrostatics(&d, problem),
            _ => solve_heat(&d, problem),
        }
        .map_err(Error::from_kernel)?
    };
    check_result_finite(&result)?;
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
