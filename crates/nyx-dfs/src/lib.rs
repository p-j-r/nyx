//! nyx-dfs: diffusion physics (`dfs`, ADR 12): electrostatics and heat transfer, stationary and transient.
//!
//! SKELETON (restructure foundation). This crate declares what `dfs` is to the rest of nyx: its studies, conditions, properties, fields
//! and couplings. The kernels (`femcore::dfs`: assemble, solve, fields) move in here in phase 4, first of the three physics, after
//! phases 2 and 3 have changed the mesh type and the linear algebra underneath them. Until then `solve` answers
//! `Error::NotYetMoved`. The descriptors below are the first version, read from the current kernels; the phase 4 owner refines
//! them and bumps `docs/architecture/contracts/physics_dfs.v1.json` with a `VERSIONS.md` entry.

use nyx_physics::Mesh;
use nyx_physics::*;

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
        vec![]
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

    fn validate(&self, spec: &ProblemSpec, _mesh: &Mesh) -> Result<(), Error> {
        // the element orders the kernels have (phase 4 adds the checks that need the mesh)
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("dfs supports element order 1 and 2, not {}", spec.order)));
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec) -> Sizing {
        // diffusion has no wavelength to resolve; geometry rules size the mesh
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }

    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        Err(Error::NotYetMoved(
            "dfs: the kernels still live in femcore::dfs and move into this crate in restructure phase 4".into(),
        ))
    }
}
