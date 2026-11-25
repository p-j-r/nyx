// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-sms: solid mechanics (`sms`, ADR 12): static and eigenfrequency studies of isotropic and anisotropic solids, in 2D (plane strain or
//! stress) and 3D, with absorbing layers (perfectly matched layers, static infinite elements, and low-reflecting boundaries) and
//! Floquet-periodic pairs, assembled natively on mixed (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) or simplex meshes.

// The kernels index several arrays with one loop variable, the form that matches the formulas (as in femcore).
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::type_complexity)]

pub mod axisym;
pub mod cross_section;
pub mod floquet;
pub mod lrb;
pub mod piezo;
pub mod piezo_asm;
pub mod piezo_driven;
pub mod piezo_mat;
pub mod pml;
pub mod problem;
pub mod solid;
pub mod waveguide;

pub use cross_section::{
    solve_axisymmetric, solve_axisymmetric_mixed, solve_solid_axisym, solve_solid_waveguide, solve_waveguide_modes,
    solve_waveguide_modes_mixed, CrossModes,
};
pub use problem::{
    decode, Bc, BcKind, ConductorRole, IeProfile, InfiniteElement, PiezoMode, Plane, Problem, Study, ThermalExpansion,
    ThermalSource,
};
pub use solid::{solve_solid, solve_solid_mixed, Setup};

use femcore::progress::{self, Progress};
use nyx_physics::Mesh;
use nyx_physics::*;
use std::sync::Arc;

pub struct Sms;

const WAFER: [&str; 3] = ["x", "y", "z"];

fn shared_coord_params() -> Vec<ParamDesc> {
    let dir = |n: &str, l: &str| {
        ParamDesc::choice(n, l, &["off", "plus", "minus"], "off", "which side of the region stretches outward")
    };
    vec![
        ParamDesc::choice(
            "coord",
            "Coordinate topology",
            &["cartesian", "cylindrical", "spherical"],
            "cartesian",
            "coordinate system along which the layer stretches",
        ),
        dir("x", "Absorb along x"),
        dir("y", "Absorb along y (3D only)"),
        dir("z", "Absorb along z"),
        ParamDesc::vector("center", "Center", "m", 3, "center [x, y, z] for cylindrical or spherical coordinates"),
        ParamDesc::vector("axis", "Cylinder axis", "", 3, "axis vector [x, y, z] for cylindrical coordinates")
            .default_vector(vec![0.0, 0.0, 1.0]),
        ParamDesc::number("r_in", "Inner radius", "m", 0.0, "inner interface radius for cylindrical/spherical layers")
            .range(Some(0.0), None),
        ParamDesc::number(
            "thickness",
            "Radial/shell thickness",
            "m",
            0.0,
            "physical radial/spherical layer thickness (0 infers from domain mesh)",
        )
        .range(Some(0.0), None),
        ParamDesc::boolean("radial_on", "Radial stretching", true, "enable cylindrical radial stretching"),
        dir("axial", "Axial stretching (cylindrical)"),
        ParamDesc::number("z_in", "Axial inner coordinate", "m", 0.0, "inner axial interface offset along axis"),
        ParamDesc::number("z_thickness", "Axial thickness", "m", 0.0, "cylindrical endcap axial thickness")
            .range(Some(0.0), None),
        ParamDesc::numbers("bounds_lo", "Cartesian lower bounds", "m", "optional explicit [x_min, y_min, z_min]"),
        ParamDesc::numbers("bounds_hi", "Cartesian upper bounds", "m", "optional explicit [x_max, y_max, z_max]"),
    ]
}

fn axes_pml_params() -> Vec<ParamDesc> {
    let mut out = shared_coord_params();
    out.extend([
        ParamDesc::number("strength", "Strength", "", 3.0, "imaginary stretch scale psi_max")
            .range(Some(0.0), None),
        ParamDesc::number("order", "Order", "", 2.0, "polynomial stretch exponent m").range(Some(0.0), None),
        ParamDesc::number("re", "Real part", "", 1.0, "constant real stretch when kappa_max = 1"),
        ParamDesc::choice(
            "profile",
            "Stretching profile",
            &["kappa_polynomial", "rational"],
            "kappa_polynomial",
            "kappa_polynomial: bounded ABCS complex dilation (eigenfrequencies); rational: singular pole mapping",
        ),
        ParamDesc::number("kappa_max", "Max real stretch", "", 1.0, "peak real evanescent stretch kappa_max >= 1")
            .range(Some(1.0), None),
        ParamDesc::number("curvature", "Rational curvature", "", 0.9, "rational pole parameter c_rat in [0, 1)")
            .range(Some(0.0), Some(0.999_999)),
        ParamDesc::boolean("auto_scale", "Auto-calibrate scaling", false, "calibrate stretch from elastic wavelengths"),
        ParamDesc::number(
            "target_reflection",
            "Target reflection",
            "",
            1e-6,
            "target normal-incidence amplitude reflection R_0",
        )
        .range(Some(1e-15), Some(0.5)),
        ParamDesc::choice(
            "stabilization",
            "Backward-wave stabilization",
            &["auto", "off", "mpml", "adiabatic"],
            "auto",
            "auto: Christoffel pre-check enables M-PML when bulk backward waves exist; mpml: tangential damping; adiabatic: stiffness loss ramp",
        ),
        ParamDesc::number("mpml_ratio", "M-PML tangential ratio", "", 0.15, "tangential-to-normal damping ratio")
            .range(Some(0.0), Some(2.0)),
        ParamDesc::number("adiabatic_eta", "Adiabatic loss peak", "", 0.5, "peak stiffness loss factor eta_max")
            .range(Some(0.0), None),
        ParamDesc::boolean(
            "clamp_outer",
            "Clamp outer PML boundary",
            true,
            "apply u = 0 Dirichlet clamping on exterior PML termination facets",
        ),
    ]);
    out
}

fn infinite_element_params() -> Vec<ParamDesc> {
    let mut out = shared_coord_params();
    out.extend([
        ParamDesc::number(
            "pole_dist",
            "Pole distance",
            "m",
            0.0,
            "distance from inner IE interface to multipole center (0 auto-computes)",
        )
        .range(Some(0.0), None),
        ParamDesc::choice(
            "profile",
            "Infinite-element profile",
            &["rational", "exponential"],
            "rational",
            "rational: exact [R_in, R_out) -> [R_in, inf) mapping; exponential: bounded real stretch + outer 1st-order spring Robin",
        ),
        ParamDesc::number(
            "decay",
            "Exponential stretch ratio",
            "",
            30.0,
            "outer-to-inner slope ratio chi = exp(beta) for exponential profile",
        )
        .range(Some(1.0), None),
        ParamDesc::boolean(
            "clamp_outer",
            "Clamp outer boundary",
            true,
            "clamp u = 0 at infinity for rational infinite elements",
        ),
    ]);
    out
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
            StudyDesc { kind: "static".into(), name: "Static".into(), doc: "Displacement and stress under fixed constraints, loads, thermal expansion, and static infinite elements.".into(), dims: Dims::Both, params: vec![] },
            StudyDesc {
                kind: "eigenfrequency".into(),
                name: "Eigenfrequencies".into(),
                doc: "Vibration modes. With absorbing layers or low-reflecting boundaries the frequencies are complex (signed Q = Re f / (2 Im f)).".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 0.0, "0 finds the lowest modes").range(Some(0.0), None),
                ],
            },
            StudyDesc {
                kind: "axisymmetric_eigenfrequency".into(),
                name: "Axisymmetric eigenfrequencies".into(),
                doc: "Vibration modes of a body of revolution; the mesh is the (r, z) half plane with x = r, and the displacement goes as exp(-j m theta). With absorbing layers the frequencies are complex: Im f > 0 and Q = Re f / (2 Im f); the scalar 'layer fraction N' tells a resonator mode from a layer mode. Every mode also carries 'f re N (Hz)' and 'f im N (Hz)' as numbers and its energy shares outside the layers: 'radial share N', 'axial share N', 'azimuthal share N' and 'share in <region> N'. At m = 0 the torsional modes (azimuthal displacement only, azimuthal share 1) are solved together with the radial and axial ones. The modes returned are the converged eigenpairs nearest f_search: their number can be smaller or larger than n_modes.".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 0.0, "0 finds the lowest modes").range(Some(0.0), None),
                    ParamDesc::integer("azimuthal_order", "Azimuthal order", 0, -100000, "the displacement goes as exp(-j m theta); m > 0 runs along +y")
                        .range(Some(-100000.0), Some(100000.0)),
                ],
            },
            StudyDesc {
                kind: "modes".into(),
                name: "Waveguide modes".into(),
                doc: "Elastic waveguide modes on the 2D (x, z) cross-section; the displacement has three components (ux, uz, uy) and goes as exp(-j q y) with longitudinal wavevector q. With absorbing layers or low-reflecting boundaries the frequencies are complex: Im f > 0 and Q = Re f / (2 Im f); the scalar 'layer fraction N' tells a guided mode from a layer mode. Every mode outputs the fields '|u|', 'Re ux', 'Re uz', 'Re uy', '|uy|' and the scalars 'f re N (Hz)', 'f im N (Hz)', 'QN' (with loss), 'layer fraction N', 'lateral share N', 'vertical share N', 'longitudinal share N' and 'share in <region> N'.".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc::number("f_search", "Search around", "Hz", 0.0, "0 finds the lowest modes").range(Some(0.0), None),
                    ParamDesc::number(
                        "wavevector",
                        "Longitudinal wavevector",
                        "rad/m",
                        0.0,
                        "wavevector q along +y; the displacement goes as exp(-j q y)",
                    ),
                ],
            },
            StudyDesc {
                kind: "driven".into(),
                name: "Driven frequency sweep".into(),
                doc: "Driven piezoelectric frequency-domain multi-port admittance (Y), impedance (Z), and scattering (S) sweep on 2D or 3D meshes with electrical 'terminal' electrodes, accelerated by Hermite/Galerkin reduced-order modeling.".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::number("f_min", "Start frequency", "Hz", 1.0e9, "minimum sweep frequency (required unless freqs_hz is given)").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("f_max", "Stop frequency", "Hz", 2.0e9, "maximum sweep frequency (>= f_min, required unless freqs_hz is given)").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("n_points", "Frequency points", 201, 1, "number of uniformly spaced frequencies in [f_min, f_max]"),
                    ParamDesc::numbers("freqs_hz", "Explicit frequencies", "Hz", "optional explicit ascending frequency grid in Hz; overrides f_min, f_max, and n_points when non-empty"),
                    ParamDesc::number("prom_tol", "PROM residual tolerance", "", 1e-3, "relative residual tolerance for Hermite/Galerkin PROM fast sweep (0 or null disables PROM and solves point-by-point HDM)").range(Some(0.0), None),
                    ParamDesc::integer("max_hdm_solves", "Maximum HDM solves", 16, 1, "maximum number of greedy full-order HDM anchor factorizations during PROM sweep").range(Some(1.0), Some(64.0)),
                    ParamDesc::boolean("strict_prom_tol", "Strict PROM tolerance", false, "return an error if max_prom_residual still exceeds prom_tol after max_hdm_solves"),
                    ParamDesc::numbers("seed_freqs_hz", "Seed anchor frequencies", "Hz", "optional initial PROM anchor frequencies merged with [f_min, f_mid, f_max]"),
                    ParamDesc::choice(
                        "solver",
                        "Linear solver backend",
                        &["auto", "direct", "lblt", "blr"],
                        "auto",
                        "sparse linear solver backend for driven frequency-domain solves ('auto', 'direct', 'lblt', or 'blr')",
                    ),
                    ParamDesc::number(
                        "blr_tol",
                        "BLR compression tolerance",
                        "",
                        1e-6,
                        "relative low-rank compression tolerance for the 'blr' multifrontal LU backend",
                    )
                    .range(Some(0.0), None),
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
                doc: "Absorbing layer by curvilinear complex coordinate stretching in the region (eigenfrequency studies).".into(),
                applies: Applies::Domain,
                params: axes_pml_params(),
            },
            BcDesc {
                kind: "infinite_element".into(),
                name: "Infinite element".into(),
                doc: "Real coordinate-stretching far-field domain mapping [R_in, R_out) -> [R_in, infinity) for static elasticity.".into(),
                applies: Applies::Domain,
                params: infinite_element_params(),
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
            BcDesc {
                kind: "ground".into(),
                name: "Ground".into(),
                doc: "Fixes electric potential phi = 0 V on the boundary (piezoelectric / electro-elastic coupling).".into(),
                applies: Applies::Boundary,
                params: vec![],
            },
            BcDesc {
                kind: "potential".into(),
                name: "Electric potential".into(),
                doc: "Prescribes electric potential phi = value (V) on the boundary (0 V in eigenfrequency studies).".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::number("value", "Potential", "V", 1.0, "")],
            },
            BcDesc {
                kind: "surface_charge".into(),
                name: "Surface charge".into(),
                doc: "Prescribes surface free charge density sigma_s (C/m^2) on the boundary (static study).".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::number("density", "Surface charge density", "C/m^2", 0.0, "")],
            },
            BcDesc {
                kind: "floating_potential".into(),
                name: "Floating potential".into(),
                doc: "Constrains all potential dofs on the boundary to a single equipotential value with prescribed net charge (0 C in eigenfrequency studies).".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::number(
                    "charge",
                    "Total charge",
                    "C",
                    0.0,
                    "prescribed net free charge on the floating electrode",
                )],
            },
            BcDesc {
                kind: "terminal".into(),
                name: "Electrical terminal".into(),
                doc: "Multi-port driven electrical terminal boundary (applies harmonic excitation V_p = 1 V to compute multi-port admittance Y, Z, and S in the driven study).".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique port identifier").required(),
                    ParamDesc::number("Z_0", "Reference impedance", "Ohm", 50.0, "reference port impedance for S-parameters")
                        .range(Some(f64::MIN_POSITIVE), None),
                ],
            },
            BcDesc {
                kind: "conductor".into(),
                name: "Conductor".into(),
                doc: "Marks a domain as an equipotential metallic electrode (skipped in K_{phi phi} and K_{u phi}, with role floating, ground, fixed, or terminal).".into(),
                applies: Applies::Domain,
                params: vec![
                    ParamDesc::choice(
                        "role",
                        "Electrical role",
                        &["floating", "ground", "fixed", "terminal"],
                        "floating",
                        "floating: unknown equipotential with prescribed net charge; ground: phi = 0 V; fixed: phi = value; terminal: multi-port driven terminal",
                    ),
                    ParamDesc::number("value", "Potential", "V", 0.0, "electric potential when role = fixed"),
                    ParamDesc::number("charge", "Total charge", "C", 0.0, "net free charge when role = floating (static only)"),
                    ParamDesc::text("name", "Port name", "", "port identifier when role = terminal"),
                    ParamDesc::number("Z_0", "Reference impedance", "Ohm", 50.0, "reference port impedance when role = terminal")
                        .range(Some(f64::MIN_POSITIVE), None),
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
            PropDesc::new("eps_r", "Relative permittivity (clamped)", "", Tensor, &[], Some(1.0)),
            PropDesc::new("eta", "Mechanical viscosity", "Pa*s", Scalar, &["eigenfrequency"], Some(0.0)),
            PropDesc::new("piezo_e", "Piezoelectric stress-charge tensor", "C/m^2", Matrix3x6, &[], None),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![
            ParamDesc::choice("plane", "Plane assumption (2D)", &["strain", "stress"], "strain", ""),
            ParamDesc::number("thickness", "Thickness (2D plane stress)", "m", 1.0, "")
                .range(Some(f64::MIN_POSITIVE), None),
            ParamDesc::choice(
                "piezoelectric",
                "Piezoelectric coupling",
                &["auto", "on", "off"],
                "auto",
                "auto: enable coupled (u, phi) formulation when electrical boundary conditions are present; on: force coupled formulation; off: ignore piezoelectric coupling",
            ),
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
                FieldDesc::new("phi", "V", &[], "electric potential (piezoelectric / electro-elastic)"),
                FieldDesc::new("E", "V/m", &WAFER, "electric field -grad phi (piezoelectric / electro-elastic)"),
                FieldDesc::new("D", "C/m^2", &WAFER, "electric displacement D (piezoelectric / electro-elastic)"),
            ],
            "eigenfrequency" | "driven" => vec![
                FieldDesc::new("u", "a.u.", &WAFER, "displacement of a mode (normalised)"),
                FieldDesc::new("phi", "a.u.", &[], "electric potential of a mode (piezoelectric)"),
            ],
            "axisymmetric_eigenfrequency" => vec![
                FieldDesc::new("|u|", "a.u.", &[], "magnitude of the displacement of a mode, scaled to a maximum of 1"),
                FieldDesc::new(
                    "Re ux",
                    "a.u.",
                    &[],
                    "real part of the radial component, phase chosen so that the largest component is real",
                ),
                FieldDesc::new("Re uz", "a.u.", &[], "real part of the axial component"),
                FieldDesc::new("Re uy", "a.u.", &[], "real part of the azimuthal component"),
                FieldDesc::new("|uy|", "a.u.", &[], "magnitude of the azimuthal component"),
                FieldDesc::new(
                    "Re phi",
                    "a.u.",
                    &[],
                    "real part of the electric potential (coupled piezoelectric solves)",
                ),
                FieldDesc::new(
                    "|phi|",
                    "a.u.",
                    &[],
                    "magnitude of the electric potential (coupled piezoelectric solves)",
                ),
            ],
            "modes" => vec![
                FieldDesc::new("|u|", "a.u.", &[], "magnitude of the displacement of a mode, scaled to a maximum of 1"),
                FieldDesc::new(
                    "Re ux",
                    "a.u.",
                    &[],
                    "real part of the lateral (x) component, phase chosen so that the largest component is real",
                ),
                FieldDesc::new("Re uz", "a.u.", &[], "real part of the vertical (z) component"),
                FieldDesc::new("Re uy", "a.u.", &[], "real part of the longitudinal (y) component"),
                FieldDesc::new("|uy|", "a.u.", &[], "magnitude of the longitudinal (y) component"),
                FieldDesc::new("Re phi", "a.u.", &[], "real part of the electric potential of a mode (piezoelectric)"),
                FieldDesc::new("|phi|", "a.u.", &[], "magnitude of the electric potential of a mode (piezoelectric)"),
            ],
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
        if spec.study.kind == "axisymmetric_eigenfrequency" && mesh.dim() != 2 {
            return Err(Error::invalid(format!(
                "the study 'axisymmetric_eigenfrequency' needs a 2D mesh (the (r, z) half plane), this one is {}D",
                mesh.dim()
            )));
        }
        if spec.study.kind == "modes" && mesh.dim() != 2 {
            return Err(Error::invalid(format!(
                "the study 'modes' needs a 2D mesh (the (x, z) cross-section), this one is {}D",
                mesh.dim()
            )));
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec, library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        use femcore::material::{
            coupled::EPSILON_0, isotropic_stiffness, piezo_e_global, stiffness_global, tensor_global,
        };
        use femcore::meshing::physics_size::{elastic_wavelength, piezo_wave_speed_range};
        let headroom = if spec.study.kind == "driven" { 1.0 } else { 1.5 };
        let empty = Sizing { order: spec.order, headroom, wavelengths: Default::default() };
        let studies = self.studies();
        let Some(sd) = studies.iter().find(|s| {
            matches!(s.kind.as_str(), "eigenfrequency" | "axisymmetric_eigenfrequency" | "modes" | "driven")
                && spec.study.kind == s.kind
        }) else {
            return empty;
        };
        let sp = nyx_physics::decode::Params::new(&spec.study.params, &sd.params);
        let f = if spec.study.kind == "driven" {
            sp.optional_numbers("freqs_hz")
                .and_then(|fs| fs.into_iter().filter(|f| f.is_finite() && *f > 0.0).reduce(f64::max))
                .unwrap_or_else(|| sp.number("f_max"))
        } else {
            sp.number("f_search")
        };
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
            let c = match m
                .aniso
                .get(&label)
                .and_then(|a| a.stiffness.as_ref().map(|s| stiffness_global(s, m.aniso[&label].euler_deg)))
            {
                Some(c) => c,
                None => match (m.get(label, "E"), m.get(label, "nu")) {
                    (Some(e), Some(nu)) if e > 0.0 && nu > -1.0 && nu < 0.5 => isotropic_stiffness(e, nu),
                    _ => continue,
                },
            };
            let lam = if let Some(e_tb) = m.aniso.get(&label).and_then(|a| a.piezo_e.as_ref()) {
                let euler = m.aniso[&label].euler_deg;
                let e_sol = piezo_e_global(e_tb, euler);
                let eps_rel = m
                    .aniso
                    .get(&label)
                    .and_then(|a| a.permittivity.as_ref().map(|t| tensor_global(t, euler)))
                    .unwrap_or_else(|| {
                        let er = m.get(label, "eps_r").unwrap_or(1.0).max(1e-12);
                        [[er, 0.0, 0.0], [0.0, er, 0.0], [0.0, 0.0, er]]
                    });
                let mut eps_abs = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        eps_abs[a][b] = EPSILON_0 * eps_rel[a][b];
                    }
                }
                let (vs_min, _) = piezo_wave_speed_range(&c, &e_sol, &eps_abs, rho);
                vs_min / f
            } else {
                elastic_wavelength(&c, rho, f)
            };
            wavelengths.insert(region.clone(), lam);
        }
        Sizing { order: spec.order, headroom, wavelengths }
    }

    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        validate_spec(self, spec)?;
        self.validate(spec, mesh)?;
        let problem = decode(spec, ctx)?;
        run(mesh, &problem, &ctx.progress)
    }
}

/// Solves a typed problem: checks the numbers, assembles natively on `MixedDiscretisation` whenever
/// non-simplex (`Quad`/`Prism`/`Hex`), curved, or static infinite-element cells are present (while
/// preserving the `SimplexMesh` fast path on straight pure-simplex meshes), and appends the stage log.
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
    problem.check_study_open_domains()?;
    problem.check_electrical()?;
    let mut result = {
        let _s =
            progress::stage(format!("solving solid mechanics on a {}D mesh of {} cells", mesh.dim(), mesh.n_cells()));
        let mixed = discretise_mixed(mesh)?;
        problem.check_wafer(mixed.dim).map_err(Error::invalid)?;
        let use_mixed = !mixed.is_pure_simplex() || mixed.has_curved() || problem.has_infinite_elements();
        if matches!(problem.study, Study::AxisymmetricEigenfrequency { .. }) {
            if mixed.dim != 2 {
                return Err(Error::invalid("the axisymmetric study needs a 2D mesh (the (r, z) half plane)"));
            }
            problem.check_axisymmetric()?;
            if use_mixed {
                solve_axisymmetric_mixed(mesh, &mixed, problem).map_err(Error::from_kernel)?
            } else {
                let d = discretise(mesh)?;
                solve_axisymmetric(&d, problem).map_err(Error::from_kernel)?
            }
        } else if matches!(problem.study, Study::Modes { .. }) {
            if mixed.dim != 2 {
                return Err(Error::invalid("the waveguide modes study needs a 2D mesh (the (x, z) cross-section)"));
            }
            problem.check_modes()?;
            if use_mixed {
                solve_waveguide_modes_mixed(mesh, &mixed, problem).map_err(Error::from_kernel)?
            } else {
                let d = discretise(mesh)?;
                solve_waveguide_modes(&d, problem).map_err(Error::from_kernel)?
            }
        } else if matches!(problem.study, Study::Driven { .. }) {
            if use_mixed {
                piezo_driven::solve_piezo_driven_mixed(mesh, &mixed, problem).map_err(Error::from_kernel)?
            } else {
                let d = discretise(mesh)?;
                piezo_driven::solve_piezo_driven(&d, problem).map_err(Error::from_kernel)?
            }
        } else if problem.uses_piezoelectric() {
            if use_mixed {
                piezo::solve_piezo_mixed(mesh, &mixed, problem).map_err(Error::from_kernel)?
            } else {
                let d = discretise(mesh)?;
                piezo::solve_piezo(&d, problem).map_err(Error::from_kernel)?
            }
        } else if use_mixed {
            solve_solid_mixed(mesh, &mixed, problem).map_err(Error::from_kernel)?
        } else {
            let d = discretise(mesh)?;
            solve_solid(&d, problem).map_err(Error::from_kernel)?
        }
    };
    check_result_finite(&result)?;
    if !problem.uses_piezoelectric()
        && problem.piezoelectric == PiezoMode::Auto
        && problem.materials.has_piezoelectric()
    {
        result.log.push(
            "piezoelectric material constants detected with no electrical boundary conditions under piezoelectric='auto'; running pure elastic solver (set piezoelectric='on' or add electrical BCs to solve coupled (u, phi))"
                .into(),
        );
    }
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
