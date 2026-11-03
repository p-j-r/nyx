//! nyx-ems: electromagnetic waves (`ems`, ADR 12): guided modes of a 2D cross-section (real or complex, with absorbing layers and
//! tilted crystals), axisymmetric modes, and the eigenfrequencies of 3D cavities and Floquet cells.
//!
//! SKELETON (restructure foundation). It declares what `ems` is to the rest of nyx; the kernels (`femcore::ems`, `ems_pml`,
//! `ems_axisym`, `ems3d`, `optical`) move in here in phase 4, after phases 2 and 3. The study is named in the problem (waveguide `modes`
//! versus `cavity`), not guessed from the mesh dimension (ROUGH_SPOTS R10). Until the kernels arrive `solve` answers
//! `Error::NotYetMoved`. The descriptors are the first version, read from the current kernels; the phase 4 owner refines them and
//! bumps `docs/architecture/contracts/physics_ems.v1.json` with a `VERSIONS.md` entry.

use nyx_physics::Mesh;
use nyx_physics::*;

pub struct Ems;

const WAFER: [&str; 3] = ["x", "y", "z"];

impl Physics for Ems {
    fn id(&self) -> &'static str {
        "ems"
    }

    fn name(&self) -> &'static str {
        "Electromagnetic waves"
    }

    fn studies(&self) -> Vec<StudyDesc> {
        vec![
            StudyDesc {
                kind: "modes".into(),
                name: "Waveguide modes".into(),
                doc: "Guided modes of a 2D cross-section at a frequency, propagating along +y as exp(-j beta y); the result is the effective index (loss is n - j k).".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::number("freq", "Frequency", "Hz", 1.934e14, "").required().range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("n_target", "Effective index near", "", 1.0, "the modes are searched around this index"),
                    ParamDesc::boolean("evanescent", "Evanescent target", false, "search an imaginary effective index"),
                    ParamDesc::integer("n_modes", "Number of modes", 4, 1, ""),
                ],
            },
            StudyDesc {
                kind: "axisymmetric_modes".into(),
                name: "Axisymmetric modes".into(),
                doc: "Whispering-gallery and other modes of a body of revolution; the mesh is the (r, z) half plane with x = r.".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::number("freq", "Frequency", "Hz", 1.934e14, "").required().range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("azimuthal_order", "Azimuthal order", 0, 0, ""),
                    ParamDesc::number("n_target", "Effective index near", "", 1.0, ""),
                    ParamDesc::integer("n_modes", "Number of modes", 4, 1, ""),
                ],
            },
            StudyDesc {
                kind: "cavity".into(),
                name: "Cavity eigenfrequencies".into(),
                doc: "Resonant frequencies of a 3D cavity or one period of a Floquet cell.".into(),
                dims: Dims::Three,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 10e9, "").required().range(Some(f64::MIN_POSITIVE), None),
                ],
            },
        ]
    }

    fn bc_kinds(&self) -> Vec<BcDesc> {
        let dir = |n: &str, l: &str| {
            ParamDesc::choice(
                n,
                l,
                &["off", "plus", "minus"],
                "off",
                "which side of the region the layer absorbs towards",
            )
        };
        vec![
            BcDesc {
                kind: "pml".into(),
                name: "Perfectly matched layer".into(),
                doc: "Absorbing layer by coordinate stretching in the region (x and z on a cross-section, x, y and z in 3D). Outer walls not otherwise described are perfect electric conductors.".into(),
                applies: Applies::Domain,
                params: vec![
                    dir("x", "Absorb along x"),
                    dir("y", "Absorb along y (3D only)"),
                    dir("z", "Absorb along z"),
                    ParamDesc::number("strength", "Strength", "", 3.0, "stretch s = re - j strength t^order"),
                    ParamDesc::number("order", "Order", "", 2.0, ""),
                    ParamDesc::number("re", "Real part", "", 1.0, ""),
                ],
            },
            BcDesc {
                kind: "periodic".into(),
                name: "Floquet-periodic pair".into(),
                doc: "The dofs on `on` are the dofs on `dest` shifted by the cell and multiplied by exp(-j phase) (cavity study).".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::text("dest", "Destination boundary", "", "a boundary selector").required(), ParamDesc::number("phase", "Phase advance", "rad", 0.0, "k . d over the period")],
            },
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        vec![
            // phase 6 turns the optical properties into dispersion models evaluated at the study wavelength
            PropDesc::new("eps_r", "Relative permittivity", "", Tensor, &[], Some(1.0)),
            PropDesc::new("mu_r", "Relative permeability", "", Scalar, &[], Some(1.0)),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }

    fn fields(&self, study: &str) -> Vec<FieldDesc> {
        match study {
            "modes" | "axisymmetric_modes" => vec![
                FieldDesc::new(
                    "E",
                    "a.u.",
                    &WAFER,
                    "electric field of a mode (x lateral, z normal, y along the propagation)",
                ),
                FieldDesc::new("|E|", "a.u.", &[], "magnitude"),
            ],
            "cavity" => vec![FieldDesc::new("E", "a.u.", &WAFER, "electric field of a mode")],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("ems.E", "optical_field")]
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("ems.temperature", "temperature")]
    }

    fn validate(&self, spec: &ProblemSpec, _mesh: &Mesh) -> Result<(), Error> {
        // the element orders the kernels have (phase 4 adds the checks that need the mesh)
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("ems supports element order 1 and 2, not {}", spec.order)));
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec) -> Sizing {
        // phase 4 fills `wavelengths` as lambda0 / n per region, n^2 = mu_r * largest eigenvalue of the permittivity
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }

    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        Err(Error::NotYetMoved(
            "ems: the kernels still live in femcore::ems* and move into this crate in restructure phase 4".into(),
        ))
    }
}
