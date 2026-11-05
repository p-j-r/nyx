//! nyx-sms: solid mechanics (`sms`, ADR 12): static and eigenfrequency studies of isotropic and anisotropic solids, in 2D (plane strain or
//! stress) and 3D, with absorbing layers (perfectly matched layers and low-reflecting boundaries) and Floquet-periodic pairs.
//!
//! The crate holds the typed problem (`problem`, decoded from the envelope by `problem::decode`), the kernels in the three stages of
//! ADR 9 (`solid`: static and real eigenfrequencies; `pml`, `lrb`, `floquet`: the complex-frequency variants; `axisym` and
//! `cross_section`: bodies of revolution on their half-plane) and the `Physics` implementation below. Moved here from
//! `femcore::{sms, sms_pml, sms_lrb, sms_floquet, sms_axisym, sms_kz}` in restructure phase 4 (ADR 42). `femcore::periodic` (the Floquet
//! reduction, shared with 3D EM waves) stays in femcore.

// The kernels index several arrays with one loop variable, the form that matches the formulas (as in femcore).
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::type_complexity)]

pub mod axisym;
pub mod cross_section;
pub mod floquet;
pub mod lrb;
pub mod pml;
pub mod problem;
pub mod solid;

pub use cross_section::{solve_solid_axisym, CrossModes};
pub use problem::{decode, Bc, BcKind, Plane, Problem, Study, ThermalExpansion, ThermalSource};
pub use solid::{solve_solid, Setup};

use femcore::progress::{self, Progress};
use nyx_physics::Mesh;
use nyx_physics::*;
use std::sync::Arc;

pub struct Sms;

const WAFER: [&str; 3] = ["x", "y", "z"];

fn axes_pml_params() -> Vec<ParamDesc> {
    let dir = |n: &str, l: &str| {
        ParamDesc::choice(n, l, &["off", "plus", "minus"], "off", "which side of the region the layer absorbs towards")
    };
    vec![
        dir("x", "Absorb along x"),
        dir("y", "Absorb along y (3D only)"),
        dir("z", "Absorb along z"),
        ParamDesc::number("strength", "Strength", "", 3.0, "stretch s = re - j strength t^order"),
        ParamDesc::number("order", "Order", "", 2.0, ""),
        ParamDesc::number("re", "Real part", "", 1.0, ""),
    ]
}

impl Physics for Sms {
    fn id(&self) -> &'static str {
        "sms"
    }

    fn name(&self) -> &'static str {
        "Solid mechanics"
    }

    fn studies(&self) -> Vec<StudyDesc> {
        vec![
            StudyDesc { kind: "static".into(), name: "Static".into(), doc: "Displacement and stress under fixed constraints, loads and thermal expansion.".into(), dims: Dims::Both, params: vec![] },
            StudyDesc {
                kind: "eigenfrequency".into(),
                name: "Eigenfrequencies".into(),
                doc: "Vibration modes. With absorbing layers or low-reflecting boundaries the frequencies are complex (Q).".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 0.0, "0 finds the lowest modes").range(Some(0.0), None),
                ],
            },
        ]
    }

    fn bc_kinds(&self) -> Vec<BcDesc> {
        let comp = |n: &str| {
            ParamDesc::number(
                n,
                &format!("Displacement along {}", &n[1..]),
                "m",
                0.0,
                "constrains this component to the value; absent: free",
            )
        };
        vec![
            BcDesc {
                kind: "fixed".into(),
                name: "Fixed constraint".into(),
                doc: "Constrains the listed displacement components (wafer axes) to given values.".into(),
                applies: Applies::Boundary,
                params: vec![comp("ux"), comp("uy"), comp("uz")],
            },
            BcDesc {
                kind: "load".into(),
                name: "Boundary load".into(),
                doc: "Traction (force per area) in wafer axes.".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::vector("traction", "Traction", "Pa", 3, "x, y, z")],
            },
            BcDesc {
                kind: "low_reflecting".into(),
                name: "Low-reflecting boundary".into(),
                doc: "Absorbing traction for outgoing pressure and shear waves (eigenfrequency studies).".into(),
                applies: Applies::Boundary,
                params: vec![],
            },
            BcDesc {
                kind: "pml".into(),
                name: "Perfectly matched layer".into(),
                doc: "Absorbing layer by coordinate stretching in the region.".into(),
                applies: Applies::Domain,
                params: axes_pml_params(),
            },
            BcDesc {
                kind: "periodic".into(),
                name: "Floquet-periodic pair".into(),
                doc: "The dofs on `on` are the dofs on `dest` shifted by the cell and multiplied by exp(-j phase)."
                    .into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("dest", "Destination boundary", "", "a boundary selector").required(),
                    ParamDesc::number("phase", "Phase advance", "rad", 0.0, "k . d over the period"),
                ],
            },
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        vec![
            PropDesc::new("E", "Young's modulus", "Pa", Scalar, &[], None),
            PropDesc::new("nu", "Poisson's ratio", "", Scalar, &[], None),
            PropDesc::new("rho", "Density", "kg/m^3", Scalar, &["eigenfrequency"], None),
            PropDesc::new("alpha", "Thermal expansion", "1/K", Tensor, &["static"], Some(0.0)),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![
            ParamDesc::choice("plane", "Plane assumption (2D)", &["strain", "stress"], "strain", ""),
            ParamDesc::number("thickness", "Thickness (2D plane stress)", "m", 1.0, "")
                .range(Some(f64::MIN_POSITIVE), None),
            ParamDesc::choice(
                "thermal",
                "Thermal expansion",
                &["off", "uniform", "imported"],
                "off",
                "`imported` takes the temperature exported by another physics (dfs)",
            ),
            ParamDesc::number("temperature", "Uniform temperature", "K", 293.15, "used with thermal = uniform"),
            ParamDesc::number("t_ref", "Reference temperature", "K", 293.15, "the stress-free temperature"),
            ParamDesc::choice(
                "solver",
                "Linear solver",
                &["auto", "direct", "iterative"],
                "auto",
                "auto: sparse direct for small systems, multigrid-preconditioned CG (static) or LOBPCG (lowest modes) for large ones; direct: Cholesky; iterative: multigrid",
            ),
            ParamDesc::names(
                "inactive",
                "Inactive regions",
                "regions with no stiffness or mass (for example air); the dofs only they touch are held at zero",
            ),
        ]
    }

    fn fields(&self, study: &str) -> Vec<FieldDesc> {
        match study {
            "static" => vec![
                FieldDesc::new("u", "m", &WAFER, "displacement (2D: x and z)"),
                FieldDesc::new("von Mises", "Pa", &[], "equivalent stress"),
            ],
            "eigenfrequency" => vec![FieldDesc::new("u", "a.u.", &WAFER, "displacement of a mode (normalised)")],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("sms.u", "displacement")]
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("sms.temperature", "temperature")]
    }

    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error> {
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("sms supports element order 1 and 2, not {}", spec.order)));
        }
        check_mesh(mesh)?;
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec, library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        use femcore::material::{isotropic_stiffness, stiffness_global};
        use femcore::meshing::physics_size::elastic_wavelength;
        let empty = Sizing { order: spec.order, headroom: 1.5, wavelengths: Default::default() };
        // only an eigenfrequency study with a search frequency has a wave to resolve
        let studies = self.studies();
        let Some(sd) = studies.iter().find(|s| s.kind == "eigenfrequency" && spec.study.kind == s.kind) else { return empty };
        let f = nyx_physics::decode::Params::new(&spec.study.params, &sd.params).number("f_search");
        if !(f > 0.0 && f.is_finite()) {
            return empty;
        }
        let Ok((names, m)) = nyx_physics::decode::materials_by_region(spec, library, None, &self.properties()) else {
            return empty;
        };
        let inactive: Vec<String> = nyx_physics::decode::Params::new(&spec.params, &self.params()).names("inactive");
        let mut wavelengths = std::collections::BTreeMap::new();
        for (i, region) in names.iter().enumerate() {
            let label = i as i32 + 1;
            if inactive.iter().any(|n| n == region) {
                continue;
            }
            let Some(rho) = m.get(label, "rho").filter(|r| *r > 0.0) else { continue };
            let c = match m.aniso.get(&label).and_then(|a| a.stiffness.as_ref().map(|s| stiffness_global(s, m.aniso[&label].euler_deg))) {
                Some(c) => c,
                None => match (m.get(label, "E"), m.get(label, "nu")) {
                    (Some(e), Some(nu)) if e > 0.0 && nu > -1.0 && nu < 0.5 => isotropic_stiffness(e, nu),
                    _ => continue,
                },
            };
            wavelengths.insert(region.clone(), elastic_wavelength(&c, rho, f));
        }
        Sizing { order: spec.order, headroom: 1.5, wavelengths }
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
    check_mesh(mesh)?;
    problem.check_finite().map_err(Error::NonFinite)?;
    let mut result = {
        let _s =
            progress::stage(format!("solving solid mechanics on a {}D mesh of {} cells", mesh.dim(), mesh.n_cells()));
        let d = discretise(mesh)?;
        problem.check_wafer(d.dim).map_err(Error::invalid)?;
        solve_solid(&d, problem).map_err(Error::from_kernel)?
    };
    check_result_finite(&result)?;
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
