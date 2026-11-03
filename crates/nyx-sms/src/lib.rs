//! nyx-sms: solid mechanics (`sms`, ADR 12): static and eigenfrequency studies of isotropic and anisotropic solids, in 2D (plane strain or
//! stress, with an optional out-of-plane wave number) and 3D, with absorbing layers and Floquet-periodic pairs.
//!
//! SKELETON (restructure foundation). It declares what `sms` is to the rest of nyx; the kernels (`femcore::sms`, `sms_kz`, `sms_pml`,
//! `sms_lrb`, `sms_floquet`, `sms_axisym`, `periodic`) move in here in phase 4, after phases 2 and 3. Until then `solve` answers
//! `Error::NotYetMoved`. The descriptors are the first version, read from the current kernels and `ProblemSpec`; the phase 4 owner
//! refines them and bumps `docs/architecture/contracts/physics_sms.v1.json` with a `VERSIONS.md` entry.

use nyx_physics::Mesh;
use nyx_physics::*;

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
                doc: "Vibration modes. With `kz` on a 2D cross-section the modes travel along y with that wave number; with absorbing layers or low-reflecting boundaries the frequencies are complex (Q).".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 0.0, "0 finds the lowest modes").range(Some(0.0), None),
                    ParamDesc::number("kz", "Out-of-plane wave number", "rad/m", 0.0, "2D only; positive for a wave along +y. Absent: no propagation along y"),
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
            PropDesc::new("stiffness", "Elastic matrix (Voigt)", "Pa", Tensor, &[], None),
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

    fn validate(&self, spec: &ProblemSpec, _mesh: &Mesh) -> Result<(), Error> {
        // the element orders the kernels have (phase 4 adds the checks that need the mesh)
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("sms supports element order 1 and 2, not {}", spec.order)));
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec) -> Sizing {
        // phase 4 fills `wavelengths` from the slowest plane wave of each region and the search frequency
        Sizing { order: spec.order, headroom: 1.5, wavelengths: Default::default() }
    }

    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        Err(Error::NotYetMoved(
            "sms: the kernels still live in femcore::sms* and move into this crate in restructure phase 4".into(),
        ))
    }
}
