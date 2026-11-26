// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! nyx-dfs: diffusion physics (`dfs`, ADR 12): electrostatics and heat transfer, stationary and transient.
//!
//! The crate holds the typed problem (`problem`, decoded from the envelope by `problem::decode`), the kernels in the three stages of ADR 9
//! (`electro`: `assemble_electrostatics`, `solve_potential`, `electrostatics_fields`; `heat`: `assemble_heat`, `solve_heat_states`,
//! `heat_fields`) and the `Physics` implementation below. Moved here from `femcore::dfs` in restructure phase 4 (ADR 41).

// The kernels index several arrays with one loop variable, the form that matches the formulas (as in femcore).
#![allow(clippy::needless_range_loop)]

pub mod ac_loss;
mod common;
mod diag;
pub mod electro;
pub mod heat;
pub mod magneto;
pub mod planar_bem;
pub mod problem;

pub use electro::{
    assemble_electrostatics, electrostatics_fields, solve_electrostatics, solve_electrostatics_mixed, solve_potential,
    ElectroAssembly, EPS0,
};
pub use heat::{assemble_heat, heat_fields, solve_heat, solve_heat_mixed, solve_heat_states, HeatStates, HeatSystem};
pub use magneto::{solve_magnetostatics, solve_magnetostatics_mixed, MU0};
pub use planar_bem::{layered_greens_function, solve_planar_bem, LayeredStack};
pub use problem::{decode, Bc, BcKind, IeCoordKind, IeProfile, InfiniteElementSpec, Problem, Study};

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
        "Diffusion (electrostatics, heat transfer, magnetostatics)"
    }

    fn studies(&self) -> Vec<StudyDesc> {
        vec![
            StudyDesc {
                kind: "electrostatic".into(),
                name: "Electrostatics".into(),
                doc: "Potential, field, and multi-terminal capacitance matrix: div(eps grad V) = -rho_v.".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::choice(
                        "formulation",
                        "Formulation",
                        &["volumetric_fem", "planar_bem", "auto"],
                        "volumetric_fem",
                        "volumetric_fem: 2D/3D finite elements; planar_bem: 2D layered-dielectric boundary elements on planar conductor sheets; auto: planar_bem when applicable, otherwise volumetric_fem",
                    ),
                    ParamDesc::choice(
                        "singular_edges",
                        "Singular edge enrichment",
                        &["none", "quarter_point"],
                        "none",
                        "quarter_point: Barsoum 1/4-point mid-edge node shift on 3D P2 elements meeting open conductor sheet perimeter edges",
                    ),
                    ParamDesc::boolean(
                        "richardson",
                        "Richardson extrapolation",
                        false,
                        "in planar_bem, solve on 1->4 subdivided sheet panels and apply log-space SPD Richardson extrapolation",
                    ),
                    ParamDesc::integer(
                        "adapt_passes",
                        "Adaptive passes",
                        1,
                        1,
                        "number of goal-oriented h-adaptive mesh refinement passes (1 disables AMR)",
                    ),
                    ParamDesc::number(
                        "adapt_rel_tol",
                        "Adaptive relative tolerance",
                        "",
                        0.001,
                        "pass-to-pass relative capacitance convergence threshold",
                    )
                    .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer(
                        "adapt_min_converged_passes",
                        "Minimum converged passes",
                        1,
                        1,
                        "consecutive passes below adapt_rel_tol required to stop",
                    ),
                    ParamDesc::number(
                        "adapt_growth",
                        "Adaptive cell growth factor",
                        "",
                        1.35,
                        "target mesh complexity / cell growth factor per adaptive pass (> 1)",
                    )
                    .range(Some(1.0 + 1e-9), None),
                    ParamDesc::choice(
                        "adapt_mode",
                        "Adaptive meshing mode",
                        &["auto", "sheet_2d_z", "full_3d"],
                        "auto",
                        "auto: sheet_2d_z on 3D sheet stacks, full_3d on general 3D, 2D metric adaptation in 2D",
                    ),
                    ParamDesc::choice(
                        "adapt_estimator",
                        "Adaptive error estimator",
                        &["flux_jump_goal", "prager_synge", "hessian"],
                        "flux_jump_goal",
                        "a posteriori cell error indicator driving h-adaptation",
                    ),
                    ParamDesc::choice(
                        "flux_bounds",
                        "Complementary flux bounds",
                        &["none", "report", "midpoint"],
                        "none",
                        "Prager-Synge H(div) equilibrated-flux lower/upper capacitance bounds; midpoint sets C = 0.5 * (C^up + C^low)",
                    ),
                    ParamDesc::number(
                        "freq",
                        "Operating frequency",
                        "Hz",
                        0.0,
                        "operating frequency for AC quasi-static dielectric loss (G and Q_diel; 0 for DC)",
                    )
                    .range(Some(0.0), None),
                ],
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
            StudyDesc {
                kind: "magnetostatic".into(),
                name: "Magnetostatics (London / inductance)".into(),
                doc: "Magnetic vector potential, London supercurrent sheets and multi-terminal inductance matrix."
                    .into(),
                dims: Dims::Both,
                params: vec![ParamDesc::choice(
                    "formulation",
                    "Formulation",
                    &["auto", "cross_section", "planar_stream", "mixed_3d"],
                    "auto",
                    "auto: cross_section in 2D, mixed_3d in 3D",
                )],
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
            b("ground", "Ground", "Zero potential (electrostatics, magnetostatics).", vec![]),
            b(
                "surface_charge",
                "Surface charge",
                "Surface charge density (electrostatics).",
                vec![ParamDesc::number("density", "Surface charge density", "C/m^2", 0.0, "")],
            ),
            b(
                "terminal",
                "Terminal",
                "Named conductor terminal for capacitance matrix extraction (electrostatics).",
                vec![ParamDesc::text("name", "Terminal name", "", "unique terminal identifier").required()],
            ),
            b(
                "open",
                "Open boundary",
                "Spherical 1/r asymptotic Robin open boundary (3D electrostatics).",
                vec![ParamDesc::vector(
                    "origin",
                    "Far-field origin",
                    "m",
                    3,
                    "Origin r0 for spherical 1/r open-boundary Robin condition; defaults to mesh bounding-box centroid",
                )],
            ),
            b(
                "far_field",
                "Far-field open boundary",
                "Asymptotic open/infinite Robin and 2nd-order Bayliss-Turkel boundary condition (electrostatics, stationary heat, magnetostatics).",
                vec![
                    ParamDesc::integer(
                        "order",
                        "Asymptotic order",
                        1,
                        1,
                        "1: first-order multipole Robin condition; 2: 3D second-order symmetric Laplace-Beltrami correction",
                    )
                    .range(Some(1.0), Some(2.0)),
                    ParamDesc::number(
                        "decay_power",
                        "Decay power",
                        "",
                        1.0,
                        "asymptotic potential decay power gamma > 0 (V ~ r^-gamma; 1.0 for monopole, 2.0 for dipole)",
                    )
                    .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::vector(
                        "center",
                        "Multipole center",
                        "m",
                        3,
                        "optional multipole expansion center [cx, cy, cz]; omit to use conductor area-weighted centroid",
                    ),
                    ParamDesc::number(
                        "t_inf",
                        "Ambient far-field temperature",
                        "K",
                        293.15,
                        "far-field ambient temperature T_inf for stationary heat transfer",
                    ),
                ],
            ),
            BcDesc {
                kind: "infinite_element".into(),
                name: "Infinite element".into(),
                doc: "Real coordinate-stretched infinite-element domain mapping [R_in, R_out) -> [R_in, inf) (rational, locking-free on Quad/Prism/Hex) or bounded exponential stretching with outer far-field Robin (electrostatics, stationary heat, magnetostatics).".into(),
                applies: Applies::Domain,
                params: vec![
                    ParamDesc::choice(
                        "coord",
                        "Coordinate system",
                        &["auto", "cartesian", "cylindrical", "spherical"],
                        "auto",
                        "coordinate topology of the infinite-element domain",
                    ),
                    ParamDesc::choice(
                        "profile",
                        "Stretching profile",
                        &["rational", "exponential"],
                        "rational",
                        "rational: maps [R_in, R_out) -> [R_in, inf); exponential: bounded real stretching chi = exp(beta) paired with outer Order-1 far_field Robin",
                    ),
                    ParamDesc::choice(
                        "x",
                        "Stretch along x",
                        &["off", "plus", "minus"],
                        "off",
                        "Cartesian stretching direction along x",
                    ),
                    ParamDesc::choice(
                        "y",
                        "Stretch along y (3D only)",
                        &["off", "plus", "minus"],
                        "off",
                        "Cartesian stretching direction along y",
                    ),
                    ParamDesc::choice(
                        "z",
                        "Stretch along z",
                        &["off", "plus", "minus"],
                        "off",
                        "Cartesian stretching direction along z",
                    ),
                    ParamDesc::vector(
                        "center",
                        "Reference center",
                        "m",
                        3,
                        "center [x0, y0, z0] for cylindrical or spherical infinite elements",
                    ),
                    ParamDesc::vector(
                        "axis",
                        "Cylinder axis",
                        "",
                        3,
                        "unit axis vector [ax, ay, az] for cylindrical infinite elements",
                    )
                    .default_vector(vec![0.0, 0.0, 1.0]),
                    ParamDesc::number(
                        "r_in",
                        "Inner radius",
                        "m",
                        0.0,
                        "inner interface radius R_in (0 auto-detects from domain vertices)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "thickness",
                        "Physical layer thickness",
                        "m",
                        0.0,
                        "physical thickness d of the infinite-element shell (0 auto-detects from domain extent)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "pole_dist",
                        "Rational pole distance",
                        "m",
                        0.0,
                        "pole distance d_pole for rational stretching (0 or null uses R_in or 0.5 * L_char)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "decay",
                        "Exponential stretch ratio",
                        "",
                        30.0,
                        "bounded exponential stretch ratio chi = exp(beta) >= 1 (typically 20..50)",
                    )
                    .range(Some(1.0), None),
                    ParamDesc::vector(
                        "bounds_lo",
                        "Cartesian lower bounds",
                        "m",
                        3,
                        "optional explicit lower corner [x_min, y_min, z_min] of a Cartesian IE domain",
                    ),
                    ParamDesc::vector(
                        "bounds_hi",
                        "Cartesian upper bounds",
                        "m",
                        3,
                        "optional explicit upper corner [x_max, y_max, z_max] of a Cartesian IE domain",
                    ),
                ],
            },
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
            b(
                "london_sheet",
                "London superconducting sheet",
                "Thin superconducting film with sheet kinetic inductance L_k,sq (magnetostatics).",
                vec![
                    ParamDesc::number("L_k_sq", "Sheet kinetic inductance", "H", 0.0, "direct L_k,sq in H/sq when > 0")
                        .range(Some(0.0), None),
                    ParamDesc::text(
                        "material",
                        "Superconducting material",
                        "",
                        "material name from library when L_k_sq is 0",
                    ),
                    ParamDesc::number("thickness", "Film thickness", "m", 0.0, "film thickness for BCS evaluation")
                        .range(Some(0.0), None),
                    ParamDesc::number("temperature", "Temperature", "K", 0.01, "temperature for BCS evaluation")
                        .range(Some(0.0), None),
                    ParamDesc::choice(
                        "film_model",
                        "Thin-film model",
                        &["two_sided", "kerr", "thin_limit"],
                        "two_sided",
                        "thin-film kinetic inductance formula",
                    ),
                ],
            ),
            b(
                "current_terminal",
                "Current terminal",
                "Current excitation port on a conductor or sheet (magnetostatics).",
                vec![
                    ParamDesc::text("name", "Terminal name", "", "unique terminal identifier").required(),
                    ParamDesc::number("current", "Terminal current", "A", 1.0, "excitation current"),
                    ParamDesc::vector("direction", "Current direction", "", 3, "unit direction vector of current flow")
                        .default_vector(vec![1.0, 0.0, 0.0]),
                ],
            ),
            b(
                "flux_loop",
                "Fluxoid / hole loop",
                "Superconducting loop or hole carrying trapped flux or circulatory stream excitation (magnetostatics).",
                vec![
                    ParamDesc::text("name", "Loop name", "", "unique loop identifier").required(),
                    ParamDesc::number("flux_wb", "Trapped flux", "Wb", 0.0, "fluxoid value in Webers"),
                ],
            ),
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        vec![
            PropDesc::new("eps_r", "Relative permittivity", "", Tensor, &["electrostatic"], None),
            PropDesc::new("tan_delta", "Dielectric loss tangent", "", Tensor, &["electrostatic"], Some(0.0)),
            PropDesc::new("sigma", "Electrical conductivity", "S/m", Tensor, &["electrostatic"], Some(0.0)),
            PropDesc::new("rho_v", "Charge density", "C/m^3", Scalar, &["electrostatic"], Some(0.0)),
            PropDesc::new("k", "Thermal conductivity", "W/(m K)", Tensor, &["heat_stationary", "heat_transient"], None),
            PropDesc::new("rho", "Density", "kg/m^3", Scalar, &["heat_transient"], None),
            PropDesc::new("Cp", "Specific heat", "J/(kg K)", Scalar, &["heat_transient"], None),
            PropDesc::new("Q", "Heat source", "W/m^3", Scalar, &["heat_stationary", "heat_transient"], Some(0.0)),
            PropDesc::new("mu_r", "Relative permeability", "", Scalar, &["magnetostatic"], Some(1.0)),
        ]
    }

    fn params(&self) -> Vec<ParamDesc> {
        vec![
            ParamDesc::choice(
                "solver",
                "Linear solver",
                &["auto", "direct", "iterative"],
                "auto",
                "auto: sparse direct for small systems, multigrid-preconditioned CG for large 3D ones; direct: Cholesky; iterative: AMG-CG",
            ),
            ParamDesc::number(
                "tol",
                "Iterative solver tolerance",
                "",
                1e-10,
                "relative residual tolerance for AMG-preconditioned conjugate gradients",
            )
            .range(Some(1e-15), Some(1e-2)),
            ParamDesc::boolean(
                "amg_line_smooth",
                "Vertical-line AMG smoother",
                true,
                "auto-detect vertical prismatic columns in 3D meshes and apply tridiagonal Thomas line smoothing",
            ),
        ]
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
            "magnetostatic" => vec![
                FieldDesc::new("A", "T m", &[], "magnetic vector potential / stream magnitude"),
                FieldDesc::new("B", "T", &WAFER, "magnetic flux density (2D: x and z)"),
            ],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![
            FieldLink::new("dfs.T", "temperature"),
            FieldLink::new("dfs.V", "potential"),
            FieldLink::new("dfs.E", "electric_field"),
            FieldLink::new("dfs.B", "magnetic_field"),
        ]
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }

    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error> {
        if !(1..=2).contains(&spec.order) {
            return Err(Error::Unsupported(format!("dfs supports element order 1 and 2, not {}", spec.order)));
        }
        check_mesh(mesh)?;
        Ok(())
    }

    fn sizing(&self, spec: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        // diffusion has no wavelength to resolve; geometry rules size the mesh
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }

    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        validate_spec(self, spec)?;
        self.validate(spec, mesh)?;
        if ctx.frequency.is_none() && spec.study.kind == "electrostatic" {
            ctx.frequency =
                spec.study.params.get("freq").and_then(|v| v.as_f64()).filter(|&f| f > 0.0 && f.is_finite());
        }
        let problem = decode(spec, ctx)?;
        run(mesh, &problem, &ctx.progress)
    }
}

fn mesh_has_curved_edges(mesh: &Mesh) -> bool {
    match mesh {
        Mesh::D2(m) => !m.edge_nodes.is_empty(),
        Mesh::D3(m) => !m.edge_nodes.is_empty(),
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
    let d = discretise(mesh)?;
    if mesh_has_curved_edges(mesh) || d.has_curved() {
        match &problem.study {
            Study::HeatStationary | Study::HeatTransient { .. } if mesh.is_pure_simplex() => {
                return Err(Error::Unsupported(
                    "heat transfer does not support curved simplex edge nodes (R32)".into(),
                ));
            }
            Study::Magnetostatic { .. } => {
                return Err(Error::Unsupported("magnetostatics does not support curved edges (R32)".into()));
            }
            _ => {}
        }
    }
    if matches!(problem.study, Study::Electrostatic) {
        if !(problem.adapt_rel_tol > 0.0
            && problem.adapt_passes >= 1
            && problem.adapt_min_converged_passes >= 1
            && problem.adapt_growth > 1.0
            && problem.freq >= 0.0)
        {
            return Err(Error::invalid(
                "electrostatic study parameters require adapt_passes >= 1, adapt_min_converged_passes >= 1, adapt_rel_tol > 0, adapt_growth > 1, and freq >= 0",
            ));
        }
        if !matches!(problem.flux_bounds.as_str(), "none" | "report" | "midpoint") {
            return Err(Error::invalid(format!(
                "flux_bounds must be 'none', 'report', or 'midpoint', got '{}'",
                problem.flux_bounds
            )));
        }
        let need_rt =
            problem.flux_bounds != "none" || (problem.adapt_passes > 1 && problem.adapt_estimator == "prager_synge");
        if need_rt {
            if problem.singular_edges == "quarter_point" || (problem.order >= 2 && d.has_qp()) {
                return Err(Error::invalid(
                    "flux_bounds / prager_synge cannot be combined with singular_edges='quarter_point': contravariant Piola-mapped Raviart-Thomas fluxes have divergent L^2 energy on singular quarter-point (det J -> 0) CurvedTets",
                ));
            }
            let has_rho = problem
                .materials
                .cell_values(&d, "rho_v", Some(0.0))
                .map(|v| v.iter().any(|&x| x.abs() > 0.0))
                .unwrap_or(false);
            let has_sigma_s = problem.bcs.iter().any(|bc| matches!(bc.kind, BcKind::SurfaceCharge(s) if s.abs() > 0.0));
            if has_rho || has_sigma_s {
                return Err(Error::invalid(
                    "flux_bounds / prager_synge requires source-free electrostatics (rho_v == 0 and surface_charge == 0)",
                ));
            }
        }
        if problem.formulation == "planar_bem" {
            if problem.adapt_passes > 1 {
                return Err(Error::Unsupported(
                    "planar_bem is a single-pass 2D boundary-element formulation and does not support volumetric multi-pass AMR (`adapt_passes > 1`)".into(),
                ));
            }
            if problem.bcs.iter().any(|bc| matches!(bc.kind, BcKind::FarField { .. } | BcKind::InfiniteElement(_))) {
                return Err(Error::Unsupported(
                    "planar_bem models open half-spaces via its layered Green's function and does not support volumetric 'far_field' or 'infinite_element' conditions".into(),
                ));
            }
            if problem.flux_bounds != "none" {
                return Err(Error::Unsupported(
                    "planar_bem is a 2D boundary-element formulation and does not support 3D volumetric Prager-Synge 'flux_bounds'".into(),
                ));
            }
            if planar_bem::has_thick_sheet_or_trench_geometry(&d, problem) {
                return Err(Error::Unsupported(
                    "planar_bem requires zero-thickness planar conductor sheets on a flat interface and does not support finite sheet thickness > 0 or substrate etch_depth > 0".into(),
                ));
            }
        }
    }
    let use_mixed = !mesh.is_pure_simplex() || problem.has_infinite_element();
    let mut result = {
        let _s = progress::stage(format!("solving diffusion on a {}D mesh of {} cells", mesh.dim(), mesh.n_cells()));
        match problem.study {
            Study::Electrostatic => {
                if use_mixed {
                    solve_electrostatics_mixed(mesh, problem)
                } else {
                    solve_electrostatics(&d, problem)
                }
            }
            Study::Magnetostatic { .. } => {
                if use_mixed {
                    solve_magnetostatics_mixed(mesh, problem)
                } else {
                    solve_magnetostatics(&d, problem)
                }
            }
            Study::HeatStationary | Study::HeatTransient { .. } => {
                if use_mixed {
                    solve_heat_mixed(mesh, problem)
                } else {
                    solve_heat(&d, problem)
                }
            }
        }
        .map_err(Error::from_kernel)?
    };
    check_result_finite(&result)?;
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
