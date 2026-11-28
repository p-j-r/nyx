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
pub mod dgtd;
#[cfg(feature = "gpu")]
pub mod dgtd_gpu;
pub(crate) mod diag;
pub mod driven;
pub mod eme;
pub mod excitation;
pub mod fdtd;
pub mod floquet_port;
pub mod heat_source;
pub mod ibc;
pub mod mixed_ops;
pub mod modes;
pub mod modes_complex;
pub mod nonlinear;
pub mod optical;
pub mod ports;
pub mod problem;
pub mod rom;
pub mod sensitivity;

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
pub use dgtd::{
    assemble_dgtd, compute_total_energy, dg_nloc, dg_ref_nodes, dg_shape_into, dgtd_fields, run_dgtd, solve_dgtd,
    DgtdAssembly, DgtdCell, DgtdCellFace, DgtdConfig, DgtdData, DgtdDipoleSource, DgtdDispersionPole, DgtdFaceKind,
    DgtdFaceQp, DgtdFlux, DgtdFluxMonitor, DgtdFluxMonitorSpec, DgtdFluxSpectrum, DgtdPolarization, DgtdPort,
    DgtdProbe, DgtdPulse, DgtdSolve, DgtdState, DgtdTimeSeries, DgtdTimeStepper, DGTD_DATA_KEY, DGTD_FLUX_CHOICES,
    DGTD_POL_CHOICES, DGTD_PULSE_CHOICES, DGTD_STEPPER_CHOICES,
};
pub use driven::{
    assemble_emw_driven, assemble_emw_driven_mixed, solve_driven, solve_driven_mixed, solve_driven_mixed_with_options,
    solve_driven_with_options, DrivenAssembly, DrivenOptions,
};
pub use eme::{
    assemble_eme, conformal_bend_coord, eme_fields, eme_length_sweep, eme_mode_sweep, eme_nominal_segment_lengths,
    eme_periodicity_sweep, eme_set_spec_segment_lengths, evaluate_eme_length_sensitivity_from_data, solve_eme,
    solve_eme_cascade, EmeAssembly, EmeCascade, EmeCascadeNode, EmeData, EmeGridSpec, EmeInterfaceOverlap,
    EmeSweepPoint, EmeSweepReport, EmeUniqueSliceSummary, EME_DATA_KEY, EME_GRID_KIND_CHOICES,
};
pub use excitation::{BeamModel, Excitation, ExcitationSpec, PreparedExcitation};
pub use fdtd::{
    assemble_fdtd, fdtd_fields, find_resonances_prony, run_fdtd, solve_fdtd, DispersiveDomainSpec, FdtdAssembly,
    FdtdBcFace, FdtdBoundaries, FdtdDim, FdtdMonitorData, FdtdMonitorSpec, FdtdSolve, FdtdSourceSpec, FdtdStudyConfig,
    FdtdWaveform, GridOverrideBox, GridSpec, LayerRefinementSpec, Pol2d, PronyResonance, SubpixelMode,
    SusceptibilityPole, YeeFields, YeeGrid, FDTD_DATA_KEY,
};
pub use floquet_port::{assemble_floquet_ports_3d, AssembledFloquetPort, FloquetHarmonic};
pub use mixed_ops::{eigs_mixed_curl_curl_ams, solve_mixed_curl_curl_iterative, EmPmlTable};
pub use modes::{
    assemble_modes, assemble_modes_mixed, assemble_modes_with_bcs, mode_fields, mode_fields_cells, n_unknowns,
    solve_mode_pairs, solve_modes, solve_modes_tensor, solve_waveguide, solve_waveguide_mixed, EpsCell, ModeAssembly,
    Modes, C0,
};
pub use modes_complex::{
    assemble_modes_complex, assemble_modes_complex_mixed, assemble_modes_complex_with_bcs,
    biorthogonal_poynting_overlap, mode_fields_complex, mode_fields_complex_cells, solve_mode_pairs_complex,
    solve_modes_complex, solve_waveguide_complex, solve_waveguide_complex_mixed, ModeAssemblyC, ModesC,
};
pub use nonlinear::{
    assemble_nonlinear_driven, assemble_nonlinear_modes, eval_chi2_residual_and_jacobian,
    eval_chi3_residual_and_jacobian, nonlinear_driven_fields, nonlinear_modes_fields, solve_nonlinear_driven,
    solve_nonlinear_driven_system, solve_nonlinear_modes, solve_nonlinear_modes_continuation, Chi2CouplingMode,
    NonlinearContinuation, NonlinearDrivenAssembly, NonlinearDrivenSolve, NonlinearDrivenSpec, NonlinearDrivenStep,
    NonlinearModeStep, NonlinearModesAssembly, NonlinearModesSolve, NonlinearModesSpec, NonlinearProcess,
    NonlinearSolverKind, StaticRealPattern, CHI2_COUPLING_CHOICES, NONLINEAR_CONTINUATION_CHOICES,
    NONLINEAR_PROCESS_CHOICES, NONLINEAR_SOLVER_CHOICES,
};
pub use optical::{solve_optical, solve_optical_wafer};
pub use ports::{
    assemble_lumped_ports_3d, assemble_waveports_3d, validate_port_contacts_3d, AssembledLumpedPort, AssembledWavePort,
};
pub use problem::{
    decode, symmetry_facet_flags, DiffractionPortBc, LumpedJjBc, LumpedPortBc, Problem, Sheet2dBc, Study,
    SurfaceImpedanceBc, Symmetry, TdBackend, TdPrecision, ThermoOpticDomain, WavePortBc, DRIVEN_SOLVER_CHOICES,
    IBC_MODE_CHOICES, SYMMETRY_CHOICES, TD_BACKEND_CHOICES, TD_PRECISION_CHOICES,
};
pub use rom::{rom_complex_eigenfrequencies, rom_eigenfrequencies};

use femcore::progress::{self, Progress};
use nyx_physics::Mesh;
use nyx_physics::*;
use std::sync::Arc;

pub struct Ems;

/// Descriptor of a scattered-field excitation condition (FR-36): common parameters (name, direction, polarization, amplitude,
/// phase, background, scattered-field output) followed by the kind-specific `extra` ones.
fn exc_desc(kind: &str, name: &str, doc: &str, extra: Vec<ParamDesc>) -> BcDesc {
    let mut params = vec![ParamDesc::text("name", "Excitation name", "", "unique excitation identifier").required()];
    if kind != "dipole" {
        params.push(
            ParamDesc::vector(
                "direction",
                "Propagation direction",
                "",
                3,
                "unit propagation vector k_hat in wafer axes",
            )
            .default_vector(vec![0.0, 0.0, -1.0]),
        );
        params.push(
            ParamDesc::vector(
                "polarization",
                "Polarization (real part)",
                "",
                3,
                "electric polarization; its component along k_hat is removed",
            )
            .default_vector(vec![1.0, 0.0, 0.0]),
        );
        params.push(
            ParamDesc::vector(
                "polarization_im",
                "Polarization (imaginary part)",
                "",
                3,
                "imaginary part for elliptical polarization",
            )
            .default_vector(vec![0.0, 0.0, 0.0]),
        );
        params.push(ParamDesc::number(
            "amplitude",
            "Amplitude",
            "V/m",
            1.0,
            "field amplitude (plane wave) or focal amplitude (beams)",
        ));
        params.push(ParamDesc::number("phase", "Phase", "rad", 0.0, "phase at the origin / focus"));
    }
    params.extend(extra);
    params.push(
        ParamDesc::number(
            "eps_bg",
            "Background permittivity",
            "",
            1.0,
            "real part of the homogeneous background eps_r",
        )
        .range(Some(f64::MIN_POSITIVE), None),
    );
    params.push(ParamDesc::number(
        "eps_bg_im",
        "Background permittivity (imaginary)",
        "",
        0.0,
        "imaginary part (<= 0 for loss with exp(+i w t))",
    ));
    params
        .push(ParamDesc::number("mu_bg", "Background permeability", "", 1.0, "").range(Some(f64::MIN_POSITIVE), None));
    params.push(ParamDesc::boolean(
        "keep_scattered",
        "Keep scattered fields",
        true,
        "also emit the scattered cell fields Esx..Esz, Hsx..Hsz",
    ));
    BcDesc { kind: kind.into(), name: name.into(), doc: doc.into(), applies: Applies::Domain, params }
}

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
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "eme".into(),
                name: "Eigenmode expansion (EME)".into(),
                doc: "3D bidirectional eigenmode expansion along +y over uniform, explicit, composite, or periodic cross-section segments with BLAKE3 slice deduplication, exact cross-power interface scattering matrices, Redheffer star products, cylindrical bends, and thermo-optic coupling.".into(),
                dims: Dims::Three,
                params: vec![
                    ParamDesc::number("freq", "Frequency", "Hz", 1.934e14, "operating optical or microwave frequency")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("n_target", "Effective index near", "", 1.5, "default target effective index for 2D slice mode searches"),
                    ParamDesc::boolean("evanescent", "Evanescent target", false, "search around an imaginary effective index"),
                    ParamDesc::integer("n_modes", "Modes per slice", 6, 1, "default number of eigenmodes solved per unique 2D cross-section"),
                    ParamDesc::choice(
                        "grid_kind",
                        "Segment partition kind",
                        &EME_GRID_KIND_CHOICES,
                        "uniform",
                        "partition kind along +y ('uniform', 'explicit', 'periodic_subgrid', or 'json' for nested composite/periodic specs)",
                    ),
                    ParamDesc::integer("n_cells", "Uniform slice count", 10, 1, "number of uniform slices along +y (or per unit cell)"),
                    ParamDesc::numbers("lengths", "Explicit segment lengths", "m", "explicit positive segment lengths L_c in metres along +y"),
                    ParamDesc::integer("num_reps", "Periodic repetitions", 1, 1, "unit-cell repetition count N_rep for 'periodic_subgrid'"),
                    ParamDesc::number(
                        "bend_radius",
                        "Cylindrical bend radius",
                        "m",
                        0.0,
                        "optional in-plane bend radius R_bend in metres along x (0 disables bending)",
                    ),
                    ParamDesc::text(
                        "grid_json",
                        "Nested grid specification JSON",
                        "",
                        "optional serialized EmeGridSpec JSON ('uniform', 'explicit', 'composite', or 'periodic_subgrid') overriding flat grid parameters",
                    ),
                    ParamDesc::numbers(
                        "y_range",
                        "Propagation y bounds",
                        "m",
                        "optional [y0, y1] propagation interval in metres along +y; defaults to the mesh y-bounds",
                    ),
                    ParamDesc::integer(
                        "port_modes",
                        "Exposed terminal modes per end",
                        0,
                        0,
                        "number of end-slice modes exposed in DrivenPortData when no waveport conditions are set (0 exposes all solved modes)",
                    ),
                    ParamDesc::number(
                        "port_z0_ohm",
                        "Reference circuit impedance",
                        "Ohm",
                        50.0,
                        "reference impedance Z_0 for S-to-Z/Y conversion when no waveport z0_ohm is set",
                    )
                    .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::boolean(
                        "reconstruct_fields",
                        "Reconstruct 3D electric field",
                        true,
                        "reconstruct the 3D bidirectional modal superposition electric field on the mesh under port 1 excitation",
                    ),
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "nonlinear_driven".into(),
                name: "Nonlinear harmonic-balance & coupled-wave driven".into(),
                doc: "3D frequency-domain nonlinear optics solver: single-tone chi^(3) Kerr + TPA + FCA/FCD Harmonic Balance FEM (with power-ramp and pseudo-arclength optical bistability continuation) and multi-harmonic chi^(2) / chi^(3) XPM Coupled-Wave FEM (SHG, SFG, DFG, OPA, depleted or undepleted, with Gauss-point QPM poling).".into(),
                dims: Dims::Three,
                params: vec![
                    ParamDesc::choice(
                        "process",
                        "Nonlinear process",
                        &NONLINEAR_PROCESS_CHOICES,
                        "auto",
                        "nonlinear interaction kind ('auto', 'chi3', 'chi3_xpm', 'shg', 'sfg', 'dfg', or 'opa')",
                    ),
                    ParamDesc::number("freq", "Primary / pump frequency", "Hz", 1.934_144_89e14, "primary optical or microwave frequency f_1")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number(
                        "freq_2",
                        "Secondary / signal frequency",
                        "Hz",
                        0.0,
                        "secondary/signal frequency f_2 in Hz for non-degenerate 'sfg', 'dfg', 'opa', or 'chi3_xpm'",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number("power_w", "Input pump power", "W", 1.0, "incident pump optical power in Watts")
                        .range(Some(0.0), None),
                    ParamDesc::numbers("powers_w", "Explicit power sweep", "W", "optional sequence of incident pump powers in Watts for power_ramp continuation"),
                    ParamDesc::number("power_min_w", "Minimum sweep power", "W", 0.0, "starting power in Watts for multi-step sweeps or arclength continuation")
                        .range(Some(0.0), None),
                    ParamDesc::number("power_max_w", "Maximum sweep power", "W", 0.0, "maximum power in Watts for multi-step sweeps or arclength continuation")
                        .range(Some(0.0), None),
                    ParamDesc::integer("n_steps", "Sweep steps", 1, 1, "number of uniform power steps when powers_w is empty"),
                    ParamDesc::integer("drive_port", "Primary drive port index", 0, 0, "0-based port index driven by power_w"),
                    ParamDesc::numbers("port_powers_w", "Per-port pump powers", "W", "optional per-port incident pump powers in Watts (enables arbitrary polarization/ellipticity mixtures)"),
                    ParamDesc::numbers("port_phases_rad", "Per-port pump phases", "rad", "optional per-port incident pump phases in radians"),
                    ParamDesc::number("signal_power_w", "Signal / seed power", "W", 0.0, "incident signal/seed power in Watts for 'sfg', 'dfg', 'opa', 'chi3_xpm', or seeded 'shg'")
                        .range(Some(0.0), None),
                    ParamDesc::integer("signal_port", "Signal port index", 0, 0, "0-based port index for the secondary/signal wave"),
                    ParamDesc::number("signal_phase_rad", "Signal phase", "rad", 0.0, "incident phase of the secondary/signal wave in radians"),
                    ParamDesc::choice(
                        "continuation",
                        "Continuation method",
                        &NONLINEAR_CONTINUATION_CHOICES,
                        "power_ramp",
                        "'power_ramp' for monotonic power sweeps or 'arclength' for Keller pseudo-arclength continuation across bistable S-curve turning points",
                    ),
                    ParamDesc::choice(
                        "coupling",
                        "Harmonic coupling mode",
                        &CHI2_COUPLING_CHOICES,
                        "depleted",
                        "'depleted' solves the fully coupled multi-harmonic system with exact pump depletion and Manley-Rowe conservation; 'undepleted' omits back-action onto the pump",
                    ),
                    ParamDesc::choice(
                        "nl_solver",
                        "Nonlinear solver",
                        &NONLINEAR_SOLVER_CHOICES,
                        "newton",
                        "'newton' (sparse direct LU reusing symbolic factorization) or 'jfnk' (LU-preconditioned Jacobian-free Newton-Krylov)",
                    ),
                    ParamDesc::number("rtol", "Relative residual tolerance", "", 1e-10, "nonlinear convergence relative residual tolerance")
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("max_it", "Max Newton iterations", 40, 1, "maximum nonlinear iterations per step"),
                    ParamDesc::number("arclength_ds", "Initial arclength step", "", 0.08, "dimensionless initial step size Delta s for pseudo-arclength continuation")
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("arclength_max_steps", "Max arclength steps", 240, 1, "maximum continuation steps along the pseudo-arclength curve"),
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "nonlinear_modes".into(),
                name: "Nonlinear waveguide modes".into(),
                doc: "2D cross-section nonlinear waveguide mode analysis computing vectorial nonlinear parameter gamma (1/(W*m)), effective mode area A_eff (m^2), TPA figure of merit, and fixed-power Newton continuation n_eff(P).".into(),
                dims: Dims::Two,
                params: vec![
                    ParamDesc::number("freq", "Frequency", "Hz", 1.934_144_89e14, "operating optical frequency")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("n_target", "Effective index near", "", 1.8, "target effective index for the baseline linear mode search"),
                    ParamDesc::boolean("evanescent", "Evanescent target", false, "search around an imaginary effective index"),
                    ParamDesc::integer("mode_index", "Mode index", 1, 1, "1-based index of the guided mode to track across optical power P"),
                    ParamDesc::number("power_w", "Guided optical power", "W", 1.0, "guided Poynting power P in Watts when powers_w is empty")
                        .range(Some(0.0), None),
                    ParamDesc::numbers("powers_w", "Explicit guided powers", "W", "optional sequence of guided optical powers P >= 0 in Watts"),
                    ParamDesc::number("rtol", "Relative residual tolerance", "", 1e-10, "bordered Newton relative residual tolerance")
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::integer("max_it", "Max Newton iterations", 30, 1, "maximum Newton iterations per power point"),
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "dgtd".into(),
                name: "Discontinuous Galerkin time-domain (DGTD)".into(),
                doc: "High-order curved and mixed-cell 2D/3D Discontinuous Galerkin time-domain solver with pre-inverted element mass matrices, Rankine-Hugoniot upwind and energy-conserving central fluxes, penetrable zero-thickness 2D sheet jumps, dispersive ADE poles, instantaneous chi^(2)/chi^(3) Kerr + QPM + TPA/FCA nonlinearity, FEM-seeded waveguide modal ports, and power-of-two Local Time-Stepping (LTS).".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::number("t_final", "Simulation duration", "s", 1.0e-13, "total physical time-stepping duration")
                        .required()
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("cfl", "CFL factor", "", 0.8, "Courant-Friedrichs-Lewy local stability multiplier")
                        .range(Some(f64::MIN_POSITIVE), Some(2.0)),
                    ParamDesc::choice(
                        "flux",
                        "Numerical Riemann flux",
                        &DGTD_FLUX_CHOICES,
                        "upwind",
                        "element-interface Riemann flux ('upwind' Rankine-Hugoniot or 'central' energy-conserving)",
                    ),
                    ParamDesc::choice(
                        "time_stepper",
                        "Explicit time integrator",
                        &DGTD_STEPPER_CHOICES,
                        "lts",
                        "time integration scheme ('lts' multi-rate power-of-two local time-stepping, 'lserk4', or 'rk4')",
                    ),
                    ParamDesc::integer(
                        "max_lts_tiers",
                        "Maximum LTS tiers",
                        6,
                        1,
                        "maximum number of power-of-two local time-stepping CFL tiers L (dt_ell = 2^ell dt_0; 1 forces uniform stepping)",
                    )
                    .range(Some(1.0), Some(12.0)),
                    ParamDesc::choice(
                        "polarization",
                        "2D/3D polarization",
                        &DGTD_POL_CHOICES,
                        "full",
                        "active electromagnetic polarization ('full' 6-component, 2D 'tm', or 2D 'te')",
                    ),
                    ParamDesc::integer(
                        "basis_order",
                        "Polynomial basis order override",
                        0,
                        0,
                        "optional DG polynomial basis order p in 1..=3 (0 uses the problem element order)",
                    )
                    .range(Some(0.0), Some(3.0)),
                    ParamDesc::choice(
                        "pulse",
                        "Source pulse waveform",
                        &DGTD_PULSE_CHOICES,
                        "gaussian",
                        "broadband source pulse waveform ('gaussian' carrier or 'ricker' wavelet)",
                    ),
                    ParamDesc::number("f0", "Carrier / peak frequency", "Hz", 1.934e14, "source pulse center frequency")
                        .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number(
                        "fwidth",
                        "Pulse bandwidth",
                        "Hz",
                        0.0,
                        "Gaussian pulse bandwidth parameter f_width (0 defaults to 0.25 * f0)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "pulse_amplitude",
                        "Pulse amplitude",
                        "",
                        1.0,
                        "source pulse peak amplitude scale factor A (incident port peak power scales as A^2)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "f_min",
                        "Monitor start frequency",
                        "Hz",
                        0.0,
                        "optional start frequency for running DFT monitors (0 uses f0)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::number(
                        "f_max",
                        "Monitor stop frequency",
                        "Hz",
                        0.0,
                        "optional stop frequency for running DFT monitors (>= f_min)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::integer(
                        "n_freqs",
                        "Monitor frequency points",
                        1,
                        1,
                        "number of uniformly spaced DFT monitor frequencies in [f_min, f_max]",
                    ),
                    ParamDesc::numbers(
                        "freqs_hz",
                        "Explicit monitor frequencies",
                        "Hz",
                        "optional explicit ascending frequency grid in Hz for running DFT monitors",
                    ),
                    ParamDesc::integer(
                        "excited_port",
                        "Excited port index",
                        1,
                        1,
                        "1-based index of the driven waveport or lumped port",
                    ),
                    ParamDesc::text(
                        "absorbing_on",
                        "Silver-Muller absorbing boundaries",
                        "",
                        "optional comma-separated boundary selectors treated as first-order Silver-Muller absorbing boundaries",
                    ),
                    ParamDesc::text(
                        "flux_monitors_on",
                        "Surface Poynting flux monitors",
                        "",
                        "optional comma-separated boundary or interface selectors recorded as running DFT Poynting flux monitors",
                    ),
                    ParamDesc::choice(
                        "backend",
                        "Compute backend",
                        &TD_BACKEND_CHOICES,
                        "auto",
                        "time-domain compute backend ('auto', 'cpu', or 'gpu'; 'gpu' requires --features gpu)",
                    ),
                    ParamDesc::choice(
                        "precision",
                        "GPU floating-point precision",
                        &TD_PRECISION_CHOICES,
                        "auto",
                        "GPU floating-point precision ('auto' picks f64 for nonlinear runs when SHADER_F64 is supported, otherwise scaled f32)",
                    ),
                    ParamDesc::boolean(
                        "nonlinear",
                        "Enable nonlinear polarization",
                        false,
                        "opt-in activation for material [nonlinear] chi^(2)/chi^(3) Kerr, SHG, QPM, and TPA/FCA terms (false auto-skips library [nonlinear] blocks for linear speed)",
                    ),
                    sym_param(),
                ],
            },
            StudyDesc {
                kind: "fdtd".into(),
                name: "Finite-difference time-domain (FDTD)".into(),
                doc: "2D and 3D SIMD structured subpixel Yee FDTD solver with Farjadpour polarized dielectric averaging, Dey-Mittra conformal PEC, Drude/Lorentz ADE dispersion, instantaneous chi^(2)/chi^(3) Kerr + QPM + TPA/FCA nonlinearity, CFS-PML/Absorber/Bloch boundaries, FEM-seeded ModeSource and ModeMonitor ports, broadband DFT/flux/diffraction monitors, and Prony cavity resonance extraction.".into(),
                dims: Dims::Both,
                params: vec![
                    ParamDesc::number(
                        "wavelength_min",
                        "Minimum free-space wavelength",
                        "m",
                        1.5e-6,
                        "shortest free-space wavelength lambda_min resolved by the automatic Yee grid",
                    )
                    .range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number(
                        "min_steps_per_wave",
                        "Grid steps per wavelength",
                        "",
                        16.0,
                        "minimum Yee cells per medium wavelength lambda_min / n_max",
                    )
                    .range(Some(4.0), None),
                    ParamDesc::number(
                        "run_time",
                        "Simulation run time",
                        "s",
                        0.0,
                        "physical simulation duration in seconds (0 selects automatic pulse + transit window)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::integer(
                        "max_steps",
                        "Maximum time steps",
                        0,
                        0,
                        "optional cap on leapfrog time steps (0 = uncapped)",
                    ),
                    ParamDesc::number("courant", "Courant factor", "", 0.95, "CFL stability factor in (0, 1)")
                        .range(Some(1e-3), Some(0.999_999)),
                    ParamDesc::number(
                        "shutoff_energy_ratio",
                        "Auto-shutoff energy ratio",
                        "",
                        1e-9,
                        "early termination threshold on field energy relative to peak energy (0 disables early shutoff)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::choice(
                        "polarization_2d",
                        "2D polarization",
                        &["tm_ey", "te_hy", "full"],
                        "tm_ey",
                        "active Maxwell subsystem on 2D cross-section meshes",
                    ),
                    ParamDesc::choice(
                        "subpixel",
                        "Subpixel dielectric averaging",
                        &["polarized", "volumetric", "staircased"],
                        "polarized",
                        "anisotropic interface subpixel permittivity averaging mode",
                    ),
                    ParamDesc::boolean(
                        "conformal_pec",
                        "Dey-Mittra conformal PEC",
                        true,
                        "enable stabilized Dey-Mittra fractional edge/face updates on curved PEC boundaries",
                    ),
                    ParamDesc::number(
                        "port_z0_ohm",
                        "Reference circuit impedance",
                        "Ohm",
                        0.0,
                        "optional reference impedance Z_0 for DrivenPortData S->Z/Y conversion (0 uses modal wave impedance)",
                    )
                    .range(Some(0.0), None),
                    ParamDesc::numbers(
                        "freqs_hz",
                        "Monitor frequencies",
                        "Hz",
                        "optional monitor frequency grid in Hz applied to the default field monitor when config_json omits monitors",
                    ),
                    ParamDesc::text(
                        "config_json",
                        "FDTD configuration JSON",
                        "",
                        "optional serialized FdtdStudyConfig JSON specifying sources, monitors, layer refinements, overrides, ADE poles, and boundaries",
                    ),
                    ParamDesc::choice(
                        "backend",
                        "Compute backend",
                        &TD_BACKEND_CHOICES,
                        "auto",
                        "time-domain compute backend ('auto', 'cpu', or 'gpu'; 'gpu' requires --features gpu)",
                    ),
                    ParamDesc::choice(
                        "precision",
                        "GPU floating-point precision",
                        &TD_PRECISION_CHOICES,
                        "auto",
                        "GPU floating-point precision ('auto' picks f64 for nonlinear runs when SHADER_F64 is supported, otherwise scaled f32)",
                    ),
                    ParamDesc::boolean(
                        "nonlinear",
                        "Enable nonlinear polarization",
                        false,
                        "opt-in activation for material [nonlinear] chi^(2)/chi^(3) Kerr, SHG, QPM, and TPA/FCA terms (false auto-skips library [nonlinear] blocks for linear speed)",
                    ),
                    sym_param(),
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
                doc: "2D boundary-mode waveguide port on an exterior face in a 3D driven or EME study.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique port identifier").required(),
                    ParamDesc::integer("mode", "Mode index", 1, 1, "1-based index of the guided/TEM boundary mode"),
                    ParamDesc::number("n_eff_guess", "Effective index guess", "", 0.0, "optional target effective index for 2D boundary-mode selection (sorts propagating modes by |n_eff - n_eff_guess|)").range(Some(0.0), None),
                    ParamDesc::number("z0_ohm", "Reference circuit impedance", "Ohm", 0.0, "optional circuit reference impedance Z_0 for DrivenPortData and S->Z/Y conversion (0 or null uses modal wave impedance)").range(Some(0.0), None),
                ],
            },
            BcDesc {
                kind: "sheet_2d".into(),
                name: "Zero-thickness 2D material sheet".into(),
                doc: "Penetrable zero-thickness 2D material sheet (e.g. monolayer Kubo graphene or anisotropic 2D crystal) with 2x2 tangential surface conductivity tensor sigma_2D(omega, T, mu_c, B_0) in Siemens.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("material", "2D sheet material", "", "optional material from the library with a [sheet_2d] block"),
                    ParamDesc::choice("model", "Sheet dispersion model", &["kubo_graphene", "anisotropic", "constant"], "kubo_graphene", "surface conductivity model when material is omitted"),
                    ParamDesc::number("mu_c_eV", "Chemical potential", "eV", 0.3, "Dirac Fermi level / chemical potential mu_c in eV"),
                    ParamDesc::number("tau_s", "Relaxation time", "s", 0.5e-12, "carrier momentum relaxation time tau in seconds").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("temperature", "Temperature", "K", 300.0, "sheet temperature T in Kelvin").range(Some(0.0), None),
                    ParamDesc::number("v_f", "Fermi velocity", "m/s", 1.0e6, "Dirac Fermi velocity v_F in m/s").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("b_0", "Out-of-plane magnetic bias", "T", 0.0, "static out-of-plane magnetic flux density B_0 in Tesla (Hall effect)"),
                    ParamDesc::number("n_layers", "Number of layers", "", 1.0, "effective number of decoupled 2D layers").range(Some(f64::MIN_POSITIVE), None),
                    ParamDesc::number("theta_deg", "In-plane rotation", "deg", 0.0, "in-plane rotation angle of principal axes (u, v) in degrees"),
                    ParamDesc::number("sigma_uu_re", "Principal-u real conductivity", "S", 0.0, "inline real sheet conductivity along principal axis u").range(Some(0.0), None),
                    ParamDesc::number("sigma_uu_im", "Principal-u imag conductivity", "S", 0.0, "inline imaginary sheet conductivity along principal axis u"),
                    ParamDesc::number("sigma_vv_re", "Principal-v real conductivity", "S", 0.0, "inline real sheet conductivity along principal axis v").range(Some(0.0), None),
                    ParamDesc::number("sigma_vv_im", "Principal-v imag conductivity", "S", 0.0, "inline imaginary sheet conductivity along principal axis v"),
                    ParamDesc::number("c_ox", "Gate oxide capacitance", "F/m^2", 0.0, "gate oxide capacitance per unit area C_ox for electrostatic chemical-potential gating").range(Some(0.0), None),
                    ParamDesc::number("eps_ox", "Gate oxide relative permittivity", "", 0.0, "gate oxide relative permittivity eps_ox (used with d_ox when c_ox is 0)").range(Some(0.0), None),
                    ParamDesc::number("d_ox", "Gate oxide thickness", "m", 0.0, "gate oxide thickness d_ox in metres (used with eps_ox when c_ox is 0)").range(Some(0.0), None),
                    ParamDesc::number("v_cnp", "Charge-neutrality voltage", "V", 0.0, "gate voltage offset V_CNP at the Dirac charge-neutrality point"),
                    ParamDesc::number("v_gate", "Uniform gate voltage", "V", 0.0, "uniform gate potential bias V_gate when no spatially varying potential field is imported"),
                    ParamDesc::number("f_ref", "Reference frequency", "Hz", 0.0, "reference frequency for conductivity evaluation (0 uses the study frequency)").range(Some(0.0), None),
                ],
            },
            BcDesc {
                kind: "diffraction_port".into(),
                name: "Floquet diffraction port".into(),
                doc: "Planar unit-cell face of a transversely Floquet-periodic driven study expanded in orthonormal TE/TM Rayleigh-Floquet harmonics; every propagating order is an S-matrix channel ('name:(m,n)TE|TM', |S|^2 = diffraction efficiency), evanescent orders enter the exact DtN operator.".into(),
                applies: Applies::Boundary,
                params: vec![
                    ParamDesc::text("name", "Port name", "", "unique port identifier").required(),
                    ParamDesc::integer("m_max", "Max order along axis 1", -1, -1, "largest |m| of the Rayleigh harmonics along the first transverse axis (-1: automatic, all propagating orders at f_max plus one evanescent ring)"),
                    ParamDesc::integer("n_max", "Max order along axis 2", -1, -1, "largest |n| along the second transverse axis (-1: automatic)"),
                ],
            },
            exc_desc("plane_wave", "Plane-wave excitation", "Scattered-field plane wave E = A exp(i phase) p exp(-i k k_hat.(r - origin)) in the background (eps_bg, mu_bg); 'on' selects the scatterer (contrast) domains.", vec![
                ParamDesc::vector("origin", "Phase origin", "m", 3, "point where the incident phase equals 'phase'").default_vector(vec![0.0, 0.0, 0.0]),
            ]),
            exc_desc("gaussian_beam", "Gaussian-beam excitation", "Scattered-field Gaussian beam of waist w0 focused at 'focus' (exact complex-source/sink beam by default, or the paraxial model); 'on' selects the contrast domains.", vec![
                ParamDesc::vector("focus", "Focus", "m", 3, "beam focus").default_vector(vec![0.0, 0.0, 0.0]),
                ParamDesc::number("waist", "Waist radius w0", "m", 1e-6, "1/e field radius at the focus").range(Some(f64::MIN_POSITIVE), None),
                ParamDesc::choice("model", "Beam model", &["complex_source", "paraxial"], "complex_source", "exact complex-source/sink beam or first-order paraxial Gaussian"),
            ]),
            exc_desc("astigmatic_gaussian_beam", "Astigmatic Gaussian-beam excitation", "Scattered-field elliptical paraxial Gaussian beam with waists (waist_x, waist_y) at axial foci (focus_x, focus_y) relative to 'focus'.", vec![
                ParamDesc::vector("focus", "Reference point", "m", 3, "beam reference point on the axis").default_vector(vec![0.0, 0.0, 0.0]),
                ParamDesc::number("waist_x", "Waist along x_axis", "m", 1e-6, "").range(Some(f64::MIN_POSITIVE), None),
                ParamDesc::number("waist_y", "Waist along y", "m", 1e-6, "").range(Some(f64::MIN_POSITIVE), None),
                ParamDesc::number("focus_x", "Axial focus of x", "m", 0.0, "axial offset of the x waist from 'focus'"),
                ParamDesc::number("focus_y", "Axial focus of y", "m", 0.0, "axial offset of the y waist from 'focus'"),
                ParamDesc::vector("x_axis", "Transverse x axis", "", 3, "projected onto the plane normal to the direction").default_vector(vec![1.0, 0.0, 0.0]),
            ]),
            exc_desc("thin_lens_beam", "Focused thin-lens excitation", "Scattered-field vectorial angular-spectrum focus of an aplanatic lens of numerical aperture NA (Gauss-Legendre x trapezoid superposition of exact plane waves).", vec![
                ParamDesc::vector("focus", "Focus", "m", 3, "geometric focus").default_vector(vec![0.0, 0.0, 0.0]),
                ParamDesc::number("NA", "Numerical aperture", "", 0.5, "n_bg sin(theta_max)").range(Some(f64::MIN_POSITIVE), None),
                ParamDesc::vector("x_axis", "Transverse x axis", "", 3, "pupil x axis").default_vector(vec![1.0, 0.0, 0.0]),
                ParamDesc::number("filling", "Pupil filling factor", "", 0.0, "Gaussian pupil apodization exp(-(sin theta / sin theta_max)^2 / filling^2); 0 = uniform").range(Some(0.0), None),
                ParamDesc::integer("n_theta", "Polar quadrature points", 32, 2, ""),
                ParamDesc::integer("n_phi", "Azimuthal quadrature points", 48, 4, ""),
            ]),
            exc_desc("dipole", "Point-dipole excitation", "Hertzian electric (C m) or magnetic (A m^2) point dipole radiating into the background; it must lie outside the contrast domains selected by 'on'.", vec![
                ParamDesc::vector("position", "Position", "m", 3, "").default_vector(vec![0.0, 0.0, 0.0]),
                ParamDesc::vector("moment", "Moment (real part)", "", 3, "dipole moment, C m (electric) or A m^2 (magnetic)").default_vector(vec![1e-20, 0.0, 0.0]),
                ParamDesc::vector("moment_im", "Moment (imaginary part)", "", 3, "").default_vector(vec![0.0, 0.0, 0.0]),
                ParamDesc::boolean("magnetic", "Magnetic dipole", false, "magnetic instead of electric dipole"),
            ]),
        ]
    }

    fn properties(&self) -> Vec<PropDesc> {
        use PropKind::*;
        let mut props = vec![
            // phase 6 turns the optical properties into dispersion models evaluated at the study wavelength
            PropDesc::new("eps_r", "Relative permittivity", "", Tensor, &[], Some(1.0)),
            PropDesc::new("eps_r_im", "Imaginary permittivity", "", Tensor, &[], Some(0.0)),
            PropDesc::new("loss_tangent", "Dielectric loss tangent", "", Tensor, &[], Some(0.0)),
            PropDesc::new("mu_r", "Relative permeability", "", Scalar, &[], Some(1.0)),
            // FR-41 (ADR 143) inline chi^(3) Kerr, TPA, FCA/FCD, chi^(2) d_iJ, and QPM properties (shared by harmonic-balance and time-domain solvers):
            PropDesc::new(
                "n_2",
                "Nonlinear Kerr index",
                "m^2/W",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "chi3_xxxx",
                "Third-order susceptibility chi3_xxxx",
                "m^2/V^2",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "chi3_xxyy",
                "Third-order susceptibility chi3_xxyy",
                "m^2/V^2",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "chi3_xyxy",
                "Third-order susceptibility chi3_xyxy",
                "m^2/V^2",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "chi3_xyyx",
                "Third-order susceptibility chi3_xyyx",
                "m^2/V^2",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "b_over_a",
                "Maker-Terhune cross-polarization ratio B/A",
                "",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(1.0),
            ),
            PropDesc::new(
                "beta_tpa",
                "Two-photon absorption coefficient",
                "m/W",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "sigma_fca",
                "Free-carrier absorption cross-section",
                "m^2",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "k_fcd",
                "Free-carrier dispersion coefficient",
                "m^3",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "tau_c",
                "Free-carrier lifetime",
                "s",
                Scalar,
                &["nonlinear_driven", "nonlinear_modes", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "chi2_d",
                "Second-order nonlinear tensor d_iJ",
                "m/V",
                Matrix3x6,
                &["nonlinear_driven", "dgtd", "fdtd"],
                None,
            ),
            PropDesc::new(
                "qpm_sign",
                "Uniform QPM poling sign (+1 or -1)",
                "",
                Scalar,
                &["nonlinear_driven", "dgtd", "fdtd"],
                Some(1.0),
            ),
            PropDesc::new(
                "qpm_period",
                "Periodic QPM poling period",
                "m",
                Scalar,
                &["nonlinear_driven", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "Lambda_qpm",
                "Periodic QPM poling period alias",
                "m",
                Scalar,
                &["nonlinear_driven", "dgtd", "fdtd"],
                Some(0.0),
            ),
            PropDesc::new(
                "qpm_duty_cycle",
                "Periodic QPM positive-domain duty cycle",
                "",
                Scalar,
                &["nonlinear_driven", "dgtd", "fdtd"],
                Some(0.5),
            ),
            PropDesc::new(
                "qpm_phase_offset",
                "Periodic QPM spatial phase offset",
                "m",
                Scalar,
                &["nonlinear_driven", "dgtd", "fdtd"],
                Some(0.0),
            ),
        ];
        for i in 1..=3 {
            for j in 1..=6 {
                let key: &'static str = match (i, j) {
                    (1, 1) => "d11",
                    (1, 2) => "d12",
                    (1, 3) => "d13",
                    (1, 4) => "d14",
                    (1, 5) => "d15",
                    (1, 6) => "d16",
                    (2, 1) => "d21",
                    (2, 2) => "d22",
                    (2, 3) => "d23",
                    (2, 4) => "d24",
                    (2, 5) => "d25",
                    (2, 6) => "d26",
                    (3, 1) => "d31",
                    (3, 2) => "d32",
                    (3, 3) => "d33",
                    (3, 4) => "d34",
                    (3, 5) => "d35",
                    _ => "d36",
                };
                props.push(PropDesc::new(key, key, "m/V", Scalar, &["nonlinear_driven", "dgtd", "fdtd"], Some(0.0)));
            }
        }
        props
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
            "bloch" | "nonlinear_modes" => vec![
                FieldDesc::new("|E|", "a.u.", &[], "magnitude of the electric field, scaled to a maximum of 1"),
                FieldDesc::new("|Exz|", "a.u.", &[], "magnitude of the transverse (in-plane) part"),
                FieldDesc::new("|Ey|", "a.u.", &[], "magnitude of the longitudinal part (along the propagation)"),
                FieldDesc::new("Re Ex", "a.u.", &[], "real part of the lateral component, phase chosen so that the largest component is real"),
                FieldDesc::new("Re Ey", "a.u.", &[], "real part of the y component (along the propagation)"),
                FieldDesc::new("Re Ez", "a.u.", &[], "real part of the surface-normal component"),
            ],
            "axisymmetric_modes" | "cavity" | "driven" | "eme" | "nonlinear_driven" | "dgtd" | "fdtd" => vec![
                FieldDesc::new("|E|", "a.u.", &[], "magnitude of the electric field, scaled to a maximum of 1"),
                FieldDesc::new("Re Ex", "a.u.", &[], "real part of the lateral (2D axisymmetric: radial) component, phase chosen so that the largest component is real"),
                FieldDesc::new("Re Ey", "a.u.", &[], "real part of the y component (2D axisymmetric: azimuthal)"),
                FieldDesc::new("Re Ez", "a.u.", &[], "real part of the surface-normal component"),
            ],
            _ => vec![],
        }
    }

    fn exports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("ems.E", "optical_field"), FieldLink::new("ems.Q_opt", "heat_source")]
    }

    fn study_exports(&self, study: &str) -> Vec<FieldLink> {
        match study {
            "bloch" | "axisymmetric_modes" | "eme" | "nonlinear_driven" | "nonlinear_modes" | "dgtd" | "fdtd" => {
                vec![FieldLink::new("ems.E", "optical_field")]
            }
            _ => self.exports(),
        }
    }

    fn imports(&self) -> Vec<FieldLink> {
        vec![
            FieldLink::new("ems.temperature", "temperature"),
            FieldLink::new("ems.potential", "potential"),
            FieldLink::new("ems.electric_field", "electric_field"),
            FieldLink::new("ems.strain", "strain"),
        ]
    }

    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error> {
        let max_order = if spec.study.kind == "dgtd" { 3 } else { 2 };
        if !(1..=max_order).contains(&spec.order) {
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
            "modes" | "bloch" | "eme" | "nonlinear_modes" => (sp.number("freq").unwrap_or(0.0), 1.0),
            "nonlinear_driven" => {
                let f1 = sp.number("freq").unwrap_or(0.0);
                let f2 = sp.number("freq_2").unwrap_or(0.0);
                let proc = sp.text("process").unwrap_or_default();
                let f_hi = match proc.as_str() {
                    "shg" => 2.0 * f1,
                    "sfg" => f1 + f2.max(0.0),
                    "chi3_xpm" | "dfg" | "opa" => f1.max(f2),
                    _ => f1,
                };
                (f_hi, 1.0)
            }
            "cavity" | "axisymmetric_modes" => (sp.number("f_search").unwrap_or(0.0), 1.5),
            "driven" | "dgtd" => {
                let f_hi = sp
                    .optional_numbers("freqs_hz")
                    .ok()
                    .flatten()
                    .and_then(|fs| fs.into_iter().filter(|f| f.is_finite() && *f > 0.0).reduce(f64::max))
                    .or_else(|| sp.number("f_max").ok().filter(|&f| f > 0.0))
                    .or_else(|| sp.number("f0").ok())
                    .unwrap_or(0.0);
                (f_hi, 1.0)
            }
            "fdtd" => {
                let lam = sp.number("wavelength_min").unwrap_or(0.0);
                let f_hi = if lam > 0.0 && lam.is_finite() { speed_of_light() / lam } else { 0.0 };
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
        let mut problem = decode(spec, ctx)?;
        problem.apply_transferred_inputs(mesh, ctx)?;
        run(mesh, &problem, &ctx.progress)
    }

    fn sensitivity(
        &self,
        mesh: &Mesh,
        spec: &ProblemSpec,
        ctx: &mut Ctx,
        req: &SensitivityRequest,
    ) -> Result<SensitivityReport, Error> {
        crate::sensitivity::compute_sensitivity(self, mesh, spec, ctx, req)
    }
}

/// Solves a typed problem: checks the numbers, runs the kernel of the study under `progress` (stage log, cancellation) and appends the
/// stage log to the result. `Physics::solve` is this after `decode`; code that already has mesh numbers calls it directly.
pub fn run(mesh: &Mesh, problem: &Problem, progress: &Arc<Progress>) -> Result<SolveResult, Error> {
    let _installed = progress::install(progress);
    let first_line = progress.log().len();
    let max_order = if matches!(problem.study, Study::Dgtd(_)) { 3 } else { 2 };
    if !(1..=max_order).contains(&problem.order) {
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
                Study::Eme { .. } => solve_eme(mesh, problem),
                Study::NonlinearDriven(_) => {
                    let d = discretise(mesh)?;
                    validate_port_contacts_3d(&d, problem).map_err(Error::invalid)?;
                    solve_nonlinear_driven(&d, problem)
                }
                Study::NonlinearModes(_) => {
                    let d = discretise(mesh)?;
                    solve_nonlinear_modes(&d, problem)
                }
                Study::Dgtd(_) => run_dgtd(mesh, problem),
                Study::Fdtd(..) => run_fdtd(mesh, problem),
            }
            .map_err(Error::from_kernel)?
        } else {
            let d = discretise(mesh)?;
            if d.dim == 3
                && matches!(problem.study, Study::Cavity { .. } | Study::Driven { .. } | Study::NonlinearDriven(_))
            {
                validate_port_contacts_3d(&d, problem).map_err(Error::invalid)?;
            }
            match problem.study {
                Study::Modes { .. } => solve_waveguide(&d, problem),
                Study::AxisymmetricModes { .. } => solve_axisymmetric(&d, problem),
                Study::Cavity { .. } => solve_cavity(&d, problem),
                Study::Bloch { .. } => solve_bloch(&d, problem),
                Study::Driven { .. } => solve_driven(&d, problem),
                Study::Eme { .. } => solve_eme(mesh, problem),
                Study::NonlinearDriven(_) => solve_nonlinear_driven(&d, problem),
                Study::NonlinearModes(_) => solve_nonlinear_modes(&d, problem),
                Study::Dgtd(_) => run_dgtd(mesh, problem),
                Study::Fdtd(..) => run_fdtd(mesh, problem),
            }
            .map_err(Error::from_kernel)?
        }
    };
    check_result_finite(&result)?;
    result.log.extend(progress.log().into_iter().skip(first_line));
    Ok(result)
}
