// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-ems: electromagnetic waves (`ems`, ADR 12): guided modes of a 2D cross-section (real or complex, with absorbing layers and
//! tilted crystals), axisymmetric modes, and the eigenfrequencies of 3D cavities and Floquet cells.
//!
//! The crate holds the typed problem (`problem`, decoded from the envelope by `problem::decode`), the kernels in the three stages of ADR 9
//! (`modes`: real waveguide modes, `modes_complex`: complex modes with PML or a tilted crystal, `axisym`: bodies of revolution,
//! `cavity`: 3D eigenfrequencies; `optical`: the modes as a field-evaluable set, for tests of the fields) and the `Physics`
//! implementation below. Moved here from `femcore::{ems, ems_pml, ems_axisym, ems3d}` in restructure phase 4 (ADR 41, 45). The study is named
//! in the problem (waveguide `modes` versus `cavity`), not guessed from the mesh dimension (ROUGH_SPOTS R10).
//!
//! The polarization of the modes is not computed here (phase 5, ADR 56): a mode study puts the complex field of every mode in its cell fields
//! and states its frame in the scalar `field frame` (`nyx_physics::FieldFrame`); `nyx_analysis::polarization` reads that.

// The kernels index several arrays with one loop variable, the form that matches the formulas (as in femcore).
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::type_complexity)]

pub mod axisym;
pub mod bloch;
pub mod cavity;
pub mod driven;
pub mod ibc;
pub mod mixed_ops;
pub mod modes;
pub mod modes_complex;
pub mod optical;
pub mod ports;
pub mod problem;
pub mod rom;

pub use axisym::{
    assemble_axisym, assemble_axisym_mixed, assemble_axisym_with_bcs, axisym_fields_cells, power_along_theta,
    solve_axisym_modes, solve_axisym_pairs, solve_axisymmetric, solve_axisymmetric_mixed, AxiAssembly, AxiMaterial,
    AxiModes,
};
pub use bloch::{
    assemble_bloch, assemble_bloch_mixed, bloch_fields_cells, solve_bloch, solve_bloch_eigen, solve_bloch_mixed,
    BlochAssembly, BlochSolve,
};
pub use cavity::{
    assemble_emw_eigen, assemble_emw_eigen_mixed, emw_eigen_fields, solve_cavity, solve_cavity_mixed,
    solve_cavity_mixed_ams, solve_emw_modes, EmwAssembly, EmwSolve, MixedCavityMasks,
};
pub use driven::{assemble_emw_driven, assemble_emw_driven_mixed, solve_driven, solve_driven_mixed, DrivenAssembly};
pub use mixed_ops::{eigs_mixed_curl_curl_ams, solve_mixed_curl_curl_iterative, EmPmlTable};
pub use modes::{
    assemble_modes, assemble_modes_mixed, assemble_modes_with_bcs, mode_fields, mode_fields_cells, n_unknowns,
    solve_mode_pairs, solve_modes, solve_modes_tensor, solve_waveguide, solve_waveguide_mixed, EpsCell, ModeAssembly,
    Modes, C0,
};
pub use modes_complex::{
    assemble_modes_complex, assemble_modes_complex_mixed, assemble_modes_complex_with_bcs, mode_fields_complex,
    mode_fields_complex_cells, solve_mode_pairs_complex, solve_modes_complex, solve_waveguide_complex,
    solve_waveguide_complex_mixed, ModeAssemblyC, ModesC,
};
pub use optical::{solve_optical, solve_optical_wafer};
pub use ports::{
    assemble_lumped_ports_3d, assemble_waveports_3d, validate_port_contacts_3d, AssembledLumpedPort, AssembledWavePort,
};
pub use problem::{
    decode, symmetry_facet_flags, LumpedJjBc, LumpedPortBc, Problem, Study, SurfaceImpedanceBc, Symmetry, WavePortBc,
    DRIVEN_SOLVER_CHOICES, IBC_MODE_CHOICES, SYMMETRY_CHOICES,
};
pub use rom::{rom_complex_eigenfrequencies, rom_eigenfrequencies};

use femcore::progress::{self, Progress};
use nyx_physics::Mesh;
use nyx_physics::*;
use std::sync::Arc;

pub struct Ems;

impl Physics for Ems {
    fn id(&self) -> &'static str {
        "ems"
    }

    fn name(&self) -> &'static str {
        "Electromagnetic waves"
    }

    fn studies(&self) -> Vec<StudyDesc> {
        let sym_param = || {
            ParamDesc::choice(
                "symmetry",
                "Midplane symmetry",
                &SYMMETRY_CHOICES,
                "none",
                "auto-detect exterior facets on x = 0 and/or z = 0 and apply PEC or PMC symmetry (on x >= 0, z >= 0: TE0 is pec_x_pmc_z, TM0 is pmc_x_pec_z)",
            )
        };
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
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "axisymmetric_modes".into(),
                name: "Axisymmetric modes".into(),
                doc: "Whispering-gallery and other modes of a body of revolution; the mesh is the (r, z) half plane with x = r.".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::number("f_search", "Search around", "Hz", 1.934e14, "the modes nearest this frequency are returned")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("azimuthal_order", "Azimuthal order", 0, -100000, "the field goes as exp(-j m theta); m > 0 runs along +y")
                        .range(Some(-100000.0), Some(100000.0)),
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
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "bloch".into(),
                name: "3D Bloch modes".into(),
                doc: "Fixed-frequency 3D Floquet-Bloch modes of a periodic waveguide or grating unit cell along +y as exp(-j k_y y), returning complex Bloch index n_B = k_y / k_0, loss alpha, and Poynting group index n_g.".into(),
                dims: Dims::Three,
                params: vec![
                    ParamDesc::number("freq", "Frequency", "Hz", 1.934e14, "").required().range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("n_target", "Bloch index near", "", 1.5, "the modes are searched around this Bloch index").required(),
                    ParamDesc::integer("n_modes", "Number of modes", 6, 1, ""),
                    ParamDesc {
                        name: "brillouin_zone".into(),
                        label: "Brillouin zone".into(),
                        unit: "".into(),
                        kind: ParamKind::Integer { default: 1, min: None, max: None },
                        required: false,
                        doc: "Brillouin zone index m shifting the target to n_target + (m - 1) * lambda_0 / Lambda".into(),
                    },
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "driven".into(),
                name: "Driven frequency sweep".into(),
                doc: "Driven frequency-domain multi-port S/Z/Y sweep on a 3D mesh with lumped and waveguide ports, accelerated by greedy Hermite/Galerkin reduced-order modeling.".into(),
                dims: Dims::Three,
                params: vec![
                    ParamDesc::number("f_min", "Start frequency", "Hz", 4.0e9, "minimum sweep frequency (required unless freqs_hz is given)").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("f_max", "Stop frequency", "Hz", 8.0e9, "maximum sweep frequency (>= f_min, required unless freqs_hz is given)").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("n_points", "Frequency points", 201, 1, "number of uniformly spaced frequencies in [f_min, f_max]"),
                    ParamDesc::numbers("freqs_hz", "Explicit frequencies", "Hz", "optional explicit ascending frequency grid in Hz; overrides f_min, f_max, and n_points when non-empty"),
                    ParamDesc::number("prom_tol", "PROM residual tolerance", "", 1e-3, "relative residual tolerance for Hermite/Galerkin PROM fast sweep (0 or null disables PROM and solves point-by-point HDM)").range(Some(0.0), None),
                    ParamDesc::integer("max_hdm_solves", "Maximum HDM solves", 16, 1, "maximum number of greedy full-order HDM anchor factorizations during PROM sweep").range(Some(1.0), Some(64.0)),
                    ParamDesc::boolean("strict_prom_tol", "Strict PROM tolerance", false, "return an error if max_prom_residual still exceeds prom_tol after max_hdm_solves"),
                    ParamDesc::numbers("seed_freqs_hz", "Seed anchor frequencies", "Hz", "optional initial PROM anchor frequencies merged with [f_min, f_mid, f_max]"),
                    ParamDesc::choice(
                        "solver",
                        "Linear solver backend",
                        &["auto", "direct", "lblt", "blr", "ams_fgmres"],
                        "auto",
                        "sparse linear solver backend for driven frequency-domain solves ('auto', 'direct', 'lblt', 'blr', or 'ams_fgmres')",
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
                doc: "Absorbing layer by coordinate stretching in the region (Cartesian, cylindrical, or spherical; kappa-polynomial or rational profile). Outer walls not otherwise described are perfect electric conductors.".into(),
                applies: Applies::Domain,
                params: vec![
                    dir("x", "Absorb along x"),
                    dir("y", "Absorb along y (3D only)"),
                    dir("z", "Absorb along z"),
                    ParamDesc::number("strength", "Strength", "", 3.0, "stretch s = re - j strength t^order"),
                    ParamDesc::number("order", "Order", "", 2.0, ""),
                    ParamDesc::number("re", "Real part", "", 1.0, ""),
                    ParamDesc::choice(
                        "coord",
                        "Coordinate system",
                        &["auto", "cartesian", "cylindrical", "spherical"],
                        "auto",
                        "coordinate topology of the absorbing layer",
                    ),
                    ParamDesc::vector(
                        "center",
                        "Center",
                        "m",
                        3,
                        "origin [x, y, z] in wafer coordinates for cylindrical or spherical layers",
                    )
                    .default_vector(vec![0.0, 0.0, 0.0]),
                    ParamDesc::vector(
                        "axis",
                        "Cylinder axis",
                        "",
                        3,
                        "unit cylinder axis [ax, ay, az] in wafer coordinates",
                    )
                    .default_vector(vec![0.0, 0.0, 1.0]),
                    ParamDesc::number(
                        "r_in",
                        "Inner radius",
                        "m",
                        0.0,
                        "inner interface radius for cylindrical or spherical shells",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "r_thickness",
                        "Radial thickness",
                        "m",
                        0.0,
                        "radial shell thickness for cylindrical layers",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::boolean(
                        "radial_on",
                        "Absorb radially",
                        true,
                        "enable radial coordinate stretching in cylindrical layers",
                    ),
                    dir("axial", "Absorb axially"),
                    ParamDesc::number(
                        "z_in",
                        "Axial start",
                        "m",
                        0.0,
                        "axial coordinate along cylinder axis where endcap stretching starts",
                    ),
                    ParamDesc::number(
                        "z_thickness",
                        "Axial thickness",
                        "m",
                        0.0,
                        "axial endcap thickness for cylindrical layers",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "thickness",
                        "Shell thickness",
                        "m",
                        0.0,
                        "shell thickness for spherical layers (or fallback radial/axial thickness)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::choice(
                        "profile",
                        "Stretching profile",
                        &["kappa_polynomial", "rational"],
                        "kappa_polynomial",
                        "coordinate-stretching law along normalized depth t in [0, 1]",
                    ),
                    ParamDesc::number(
                        "kappa_max",
                        "Maximum real stretch",
                        "",
                        1.0,
                        "peak real evanescent stretch kappa_max >= 1 for kappa_polynomial",
                    )
                    .range(Some(1.0), None),
                    ParamDesc::number(
                        "curvature",
                        "Rational curvature",
                        "",
                        0.9,
                        "pole parameter c in [0, 1) for rational stretching",
                    )
                    .range(Some(0.0), Some(0.999_999)),
                    ParamDesc::boolean(
                        "auto_scale",
                        "Auto-scale parameters",
                        false,
                        "calibrate strength and profile automatically from study wavelength and target_reflection",
                    ),
                    ParamDesc::number(
                        "target_reflection",
                        "Target reflection",
                        "",
                        1e-6,
                        "target normal-incidence amplitude reflection R_0 when auto_scale is enabled",
                    )
                    .range(Some(1e-15), Some(0.5)),
                    ParamDesc::vector(
                        "bounds_lo",
                        "Slab lower bounds",
                        "m",
                        3,
                        "optional explicit lower [x, y, z] slab bounds for Cartesian layers",
                    ),
                    ParamDesc::vector(
                        "bounds_hi",
                        "Slab upper bounds",
                        "m",
                        3,
                        "optional explicit upper [x, y, z] slab bounds for Cartesian layers",
                    ),
                ],
            },
            BcDesc {
                kind: "periodic".into(),
                name: "Floquet-periodic pair".into(),
                doc: "The dofs on `on` are the dofs on `dest` shifted by the cell and multiplied by exp(-j phase) (cavity, bloch, and driven studies).".into(),
                applies: Applies::Boundary,
                params: vec![ParamDesc::text("dest", "Destination boundary", "", "a boundary selector").required(), ParamDesc::number("phase", "Phase advance", "rad", 0.0, "k . d over the period")],
            },
            BcDesc {
                kind: "pmc".into(),
                name: "Magnetic wall (mirror plane)".into(),
                doc: "Perfect magnetic conductor on the boundary: the field is free there (tangential H zero). On a mirror plane of the structure it keeps the modes that are even under the reflection (the tangential electric field even, the normal one odd), such as the TE-like modes of a slab with the plane at its mid-plane.".into(),
                applies: Applies::Boundary,
                params: vec![],
            },
            BcDesc {
                kind: "pec".into(),
                name: "Perfect electric conductor".into(),
                doc: "Perfect electric conductor (tangential E = 0) on an interior 2D conductive sheet or outer boundary face.".into(),
                applies: Applies::Boundary,
                params: vec![],
            },
            BcDesc {
                kind: "surface_impedance".into(),
                name: "Superconducting / resistive surface impedance".into(),
                doc: "Kinetic-inductance and surface-resistance sheet impedance Z_s = R_s + i omega L_k_sq on a 2D sheet or boundary face. When `material` is given, evaluates arbitrary-purity BCS surface impedance at `(temperature, thickness, f_ref, film_model)`.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::number("L_k_sq", "Sheet kinetic inductance", "H/sq", 0.0, "sheet kinetic inductance per square").range(Some(0.0), None),
                    ParamDesc::number("R_s", "Surface resistance", "Ohm/sq", 0.0, "sheet surface resistance per square").range(Some(0.0), None),
                    ParamDesc::text("material", "Superconducting material", "", "optional material from the library with a [superconducting] block"),
                    ParamDesc::number("thickness", "Film thickness", "m", 0.0, "film thickness for BCS surface impedance").range(Some(0.0), None),
                    ParamDesc::number("temperature", "Temperature", "K", 0.01, "operating temperature for BCS surface impedance").range(Some(0.0), None),
                    ParamDesc::number("f_ref", "Reference frequency", "Hz", 0.0, "reference frequency for Z_s evaluation (0 or null uses the study frequency)").range(Some(0.0), None),
                    ParamDesc::choice("film_model", "Thin-film model", &["two_sided", "kerr", "thin_limit"], "two_sided", "thin-film current distribution model"),
                    ParamDesc::boolean("lossless", "Lossless kinetic inductance", false, "zero surface resistance R_s after BCS evaluation so 3D cavity solves use the fast real symmetric eigensolver while retaining L_k_sq"),
                    ParamDesc::choice(
                        "ibc_mode",
                        "IBC formulation mode",
                        &["direct", "perturbative", "auto"],
                        "direct",
                        "surface-impedance formulation: 'direct' adds Robin stiffness directly to the system matrix, 'perturbative' enforces strong PEC during eigensolve and computes Wheeler surface-current shifts in Stage 3, and 'auto' selects 'perturbative' for sub-grid P1 thin films",
                    ),
                ],
            },
            BcDesc {
                kind: "lumped_jj".into(),
                name: "Lumped Josephson junction".into(),
                doc: "Lumped Josephson junction / nonlinear port sheet across two conductors in a 3D cavity study, assembled as a rank-1 terminal operator (or distributed sheet) and exported in the passive Craig-Bampton CircuitRom.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique terminal identifier").required(),
                    ParamDesc::number("L_j", "Josephson inductance", "H", 0.0, "linearized Josephson inductance").range(Some(0.0), None),
                    ParamDesc::number("E_j_Hz", "Josephson energy E_J / h", "Hz", 0.0, "Josephson energy in Hz (alternative to L_j)").range(Some(0.0), None),
                    ParamDesc::number("C_j", "Junction capacitance", "F", 0.0, "intrinsic parallel junction capacitance").range(Some(0.0), None),
                    ParamDesc::number("E_l_Hz", "Parallel inductive energy E_L / h", "Hz", 0.0, "parallel superinductance energy in Hz").range(Some(0.0), None),
                    ParamDesc::number("phi_ext_rad", "External flux bias", "rad", 0.0, "external reduced flux bias phi_ext = 2 pi Phi_ext / Phi_0"),
                    ParamDesc::number("n_g", "Offset charge", "2e", 0.0, "gate offset charge in Cooper-pair units"),
                    ParamDesc::number("squid_asymmetry", "SQUID junction asymmetry", "", 0.0, "optional DC-SQUID junction asymmetry d in [0, 1)").range(Some(0.0), Some(1.0)),
                    ParamDesc::vector("direction", "Direction", "", 3, "unit vector along junction current flow").default_vector(vec![1.0, 0.0, 0.0]),
                    ParamDesc::choice("formulation", "Formulation", &["rank1", "distributed_sheet"], "rank1", "lumped terminal operator formulation"),
                ],
            },
            BcDesc {
                kind: "lumped_port".into(),
                name: "Lumped sheet port".into(),
                doc: "Lumped voltage/current sheet port with reference impedance Z_0 across two conductors in a 3D cavity or driven study.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique port identifier").required(),
                    ParamDesc::number("Z_0", "Reference impedance", "Ohm", 50.0, "reference port impedance").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::vector("direction", "Direction", "", 3, "unit vector along port current flow").default_vector(vec![1.0, 0.0, 0.0]),
                    ParamDesc::boolean("complex_damping", "Complex port damping", false, "include resistive port damping directly in 3D cavity eigensolve (when false, cavity stays real and exports port damping in CircuitRom)"),
                    ParamDesc::choice("mode", "Port mode", &["single", "cpw_even"], "single", "single-gap port ('single') or balanced two-slot CPW even mode ('cpw_even')"),
                ],
            },
            BcDesc {
                kind: "waveport".into(),
                name: "Waveguide boundary port".into(),
                doc: "2D boundary-mode waveguide port on an exterior face in a 3D driven frequency sweep.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique port identifier").required(),
                    ParamDesc::integer("mode", "Mode index", 1, 1, "1-based index of the guided/TEM boundary mode"),
                    ParamDesc::number("n_eff_guess", "Effective index guess", "", 0.0, "optional target effective index for 2D boundary-mode selection (sorts propagating modes by |n_eff - n_eff_guess|)").range(Some(0.0), None),
                    ParamDesc::number("z0_ohm", "Reference circuit impedance", "Ohm", 0.0, "optional circuit reference impedance Z_0 for DrivenPortData and S->Z/Y conversion (0 or null uses modal wave impedance)").range(Some(0.0), None),
                ],
            },
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        vec![
            // phase 6 turns the optical properties into dispersion models evaluated at the study wavelength
            PropDesc::new("eps_r", "Relative permittivity", "", Tensor, &[], Some(1.0)),
            PropDesc::new("eps_r_im", "Imaginary permittivity", "", Tensor, &[], Some(0.0)),
            PropDesc::new("loss_tangent", "Dielectric loss tangent", "", Tensor, &[], Some(0.0)),
            PropDesc::new("mu_r", "Relative permeability", "", Scalar, &[], Some(1.0)),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }

    fn fields(&self, study: &str) -> Vec<FieldDesc> {
        // The fields are scalar arrays under the flat names of the window (wafer axes: x lateral, z normal, y along the propagation). A
        // real mode (real solver) has the signed components `Ex`, `Ey`, `Ez`; a complex one (absorbing layers, tilted crystal, loss) has
        // their real parts `Re Ex`, ... after the phase is rotated so that the largest component is real. The complex components are the
        // per-cell fields `Ex`, `Ey`, `Ez` of the solution.
        match study {
            "modes" => vec![
                FieldDesc::new("|E|", "a.u.", &[], "magnitude of the electric field, scaled to a maximum of 1"),
                FieldDesc::new("|Exz|", "a.u.", &[], "magnitude of the transverse (in-plane) part"),
                FieldDesc::new("|Ey|", "a.u.", &[], "magnitude of the longitudinal part (along the propagation)"),
                FieldDesc::new("Ex", "a.u.", &[], "lateral component (real solver; `Re Ex` for a complex mode)"),
                FieldDesc::new("Ey", "a.u.", &[], "component along the propagation (real solver; `Re Ey` for a complex mode)"),
                FieldDesc::new("Ez", "a.u.", &[], "surface-normal component (real solver; `Re Ez` for a complex mode)"),
            ],
            "bloch" => vec![
                FieldDesc::new("|E|", "a.u.", &[], "magnitude of the electric field, scaled to a maximum of 1"),
                FieldDesc::new("|Exz|", "a.u.", &[], "magnitude of the transverse (in-plane) part"),
                FieldDesc::new("|Ey|", "a.u.", &[], "magnitude of the longitudinal part (along the propagation)"),
                FieldDesc::new("Re Ex", "a.u.", &[], "real part of the lateral component, phase chosen so that the largest component is real"),
                FieldDesc::new("Re Ey", "a.u.", &[], "real part of the y component (along the propagation)"),
                FieldDesc::new("Re Ez", "a.u.", &[], "real part of the surface-normal component"),
            ],
            "axisymmetric_modes" | "cavity" | "driven" => vec![
                FieldDesc::new("|E|", "a.u.", &[], "magnitude of the electric field, scaled to a maximum of 1"),
                FieldDesc::new("Re Ex", "a.u.", &[], "real part of the lateral (2D axisymmetric: radial) component, phase chosen so that the largest component is real"),
                FieldDesc::new("Re Ey", "a.u.", &[], "real part of the y component (2D axisymmetric: azimuthal)"),
                FieldDesc::new("Re Ez", "a.u.", &[], "real part of the surface-normal component"),
            ],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("ems.E", "optical_field")]
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("ems.temperature", "temperature")]
    }

    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error> {
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("ems supports element order 1 and 2, not {}", spec.order)));
        }
        check_mesh(mesh)?;
        // the study names the kernel, the mesh must be of its dimension
        let studies = self.studies();
        if let Some(sd) = studies.iter().find(|s| s.kind == spec.study.kind) {
            let ok = match sd.dims {
                Dims::Two => mesh.dim() == 2,
                Dims::Three => mesh.dim() == 3,
                Dims::Both => true,
            };
            if !ok {
                return Err(Error::invalid(format!(
                    "the study '{}' needs a {} mesh, this one is {}D",
                    sd.kind,
                    if sd.dims == Dims::Two { "2D" } else { "3D" },
                    mesh.dim()
                )));
            }
        }
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec, library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        use femcore::material::tensor_global;
        use femcore::meshing::physics_size::{optical_wavelength, speed_of_light};
        let empty = Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() };
        let studies = self.studies();
        let Some(sd) = studies.iter().find(|s| s.kind == spec.study.kind) else { return empty };
        let sp = nyx_physics::decode::Params::new(&spec.study.params, &sd.params);
        // a mode study resolves the wave at its frequency; the eigenfrequency studies find modes above the one they search around
        let (freq, headroom) = match spec.study.kind.as_str() {
            "modes" | "bloch" => (sp.number("freq").unwrap_or(0.0), 1.0),
            "cavity" | "axisymmetric_modes" => (sp.number("f_search").unwrap_or(0.0), 1.5),
            "driven" => {
                let f_hi = sp
                    .optional_numbers("freqs_hz")
                    .ok()
                    .flatten()
                    .and_then(|fs| fs.into_iter().filter(|f| f.is_finite() && *f > 0.0).reduce(f64::max))
                    .or_else(|| sp.number("f_max").ok())
                    .unwrap_or(0.0);
                (f_hi, 1.0)
            }
            _ => return empty,
        };
        if !(freq > 0.0 && freq.is_finite()) {
            return empty;
        }
        let lam0 = speed_of_light() / freq;
        let Ok((names, m)) = nyx_physics::decode::materials_by_region(spec, library, Some(lam0), &self.properties())
        else {
            return empty;
        };
        let mut wavelengths = std::collections::BTreeMap::new();
        for (i, region) in names.iter().enumerate() {
            let label = i as i32 + 1;
            let eps = match m
                .aniso
                .get(&label)
                .and_then(|a| a.permittivity.as_ref().map(|t| tensor_global(t, a.euler_deg)))
            {
                Some(t) => t,
                None => {
                    let e = m.get(label, "eps_r").unwrap_or(1.0);
                    [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]]
                }
            };
            wavelengths.insert(region.clone(), optical_wavelength(&eps, m.get(label, "mu_r").unwrap_or(1.0), lam0));
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
        let _s = progress::stage(format!(
            "solving EM waves (mode analysis) on a {}D mesh of {} cells",
            mesh.dim(),
            mesh.n_cells()
        ));
        // a layer along y has no meaning in the plane of a cross-section (x-z)
        if mesh.dim() == 2 {
            if let Some((dom, _)) = problem.materials.pml.iter().find(|(_, p)| p.y != femcore::props::PmlDir::Off) {
                return Err(Error::invalid(format!("the perfectly matched layer of domain {dom} absorbs along y, which is out of the plane of a 2D cross-section (x-z): use z for the direction of the surface normal")));
            }
        }
        let mixed = discretise_mixed(mesh)?;
        let use_mixed = !mixed.is_pure_simplex() || (mesh.dim() == 2 && !mixed.edge_nodes.is_empty());
        if use_mixed {
            if mixed.dim == 3 && matches!(problem.study, Study::Cavity { .. } | Study::Driven { .. }) {
                mixed_ops::validate_port_contacts_3d_mixed(&mixed, problem).map_err(Error::invalid)?;
            }
            match problem.study {
                Study::Modes { .. } => solve_waveguide_mixed(mesh, &mixed, problem),
                Study::AxisymmetricModes { .. } => solve_axisymmetric_mixed(mesh, &mixed, problem),
                Study::Cavity { .. } => solve_cavity_mixed(mesh, &mixed, problem),
                Study::Bloch { .. } => solve_bloch_mixed(mesh, &mixed, problem),
                Study::Driven { .. } => solve_driven_mixed(mesh, &mixed, problem),
            }
            .map_err(Error::from_kernel)?
        } else {
            let d = discretise(mesh)?;
            if d.dim == 3 && matches!(problem.study, Study::Cavity { .. } | Study::Driven { .. }) {
                validate_port_contacts_3d(&d, problem).map_err(Error::invalid)?;
            }
            match problem.study {
                Study::Modes { .. } => solve_waveguide(&d, problem),
                Study::AxisymmetricModes { .. } => solve_axisymmetric(&d, problem),
                Study::Cavity { .. } => solve_cavity(&d, problem),
                Study::Bloch { .. } => solve_bloch(&d, problem),
                Study::Driven { .. } => solve_driven(&d, problem),
            }
            .map_err(Error::from_kernel)?
        }
    };
    check_result_finite(&result)?;
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
