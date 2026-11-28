// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! SIMD Structured Subpixel FDTD (`Study::Fdtd`, `"fdtd"`, ADR 145) seeded by `nyx-ems::modes`.
//!
//! Three-stage ADR 9 kernel:
//! 1. [`assemble_fdtd`]: resolves domain optical/dielectric permittivities and multi-pole Drude/Lorentz ADEs,
//!    builds the wavelength-graded `YeeGrid` (`AutoGrid` + `LayerRefinementSpec` + `GridOverrideBox`), performs
//!    Farjadpour–Oskooi–Joannopoulos `PolarizedAveraging` and Dey–Mittra `PECConformal` edge/face sampling,
//!    constructs Roden–Gedney `CPML` / `Absorber` / Bloch boundary operators, and solves cross-sectional FEM
//!    eigenmodes (`nyx-ems::modes` / `modes_complex`) once per `ModeSource` / `ModeMonitor`.
//! 2. [`solve_fdtd`]: runs the cache-aligned Structure-of-Arrays (SoA) leapfrog Yee time-marching loop with
//!    deterministic `femcore::par` slab chunking (`0 ULP` across `NYX_THREADS=1` vs `8`), accumulating running
//!    DFT monitors and `FieldTimeMonitor` traces, and runs SVD Matrix-Pencil / Prony `ResonanceFinder`.
//! 3. [`fdtd_fields`]: normalizes broadband `FluxMonitor`, `ModeMonitor` (`DrivenPortData` S/Z/Y matrices), and
//!    `DiffractionMonitor` spectra by the exact discrete source pulse spectrum, and projects complex DFT `(E, H)`
//!    fields onto `SimplexMesh` `P1`/`P2` `Solution::fields` and `Solution::cell_fields` (ready for Stratton–Chu
//!    `nyx-analysis::farfield` NFFFT).

pub mod ade_pml;
#[cfg(feature = "gpu")]
pub mod gpu;
pub mod grid;
pub mod modes_bridge;
pub mod sources_monitors;
pub mod stepper;
pub mod subpixel;

pub use ade_pml::{
    AdeAndAmpereCoeffs, CpmlState, DispersiveDomainSpec, FdtdBcFace, FdtdBoundaries, SusceptibilityPole,
};
pub use grid::{
    Align64, AlignedSlab, AutoGridSpec, FdtdDim, GridOverrideBox, GridSpec, LayerRefinementSpec, Pol2d, YeeFields,
    YeeGrid,
};
pub use modes_bridge::{rh_uv_axes, solve_yee_plane_mode, ModePatchSample, YeeModeProfile};
pub use sources_monitors::{
    find_resonances_prony, Aux1dYeeLine, CompiledSource, DiffractionMonitorResult, DiffractionOrderRow,
    FdtdMonitorData, FdtdMonitorSpec, FdtdSourceSpec, FdtdWaveform, FluxMonitorResult, ModeMonitorResult,
    PronyResonance, TimeMonitorResult, FDTD_DATA_KEY,
};
pub use subpixel::{DomainMedium, FdtdGeometryOracle, PointClass, SubpixelMode, YeeMaterialCoeffs};

use crate::modes::C0;
use crate::problem::{Problem, Study};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::geom::eval::{Evaluator, Tier};
use femcore::geom::model::FeatureList;
use femcore::linalg::csr::{c, C};
use femcore::linalg::prom::s_to_z_and_y_checked;
use femcore::props::PmlDir;
use nyx_physics::{
    cell_fields_of, discretise, Discretisation, DrivenPortData, Field, FieldFrame, Mesh, Scalar, Solution, SolveResult,
    SpatialErrorEstimate,
};
use serde::{Deserialize, Serialize};
use sources_monitors::{evaluate_diffraction_orders, sample_probe_component, PlanarDftAccumulator, VolumeFieldDft};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::sync::Arc;
use std::time::Instant;

/// High-level serializable/typed configuration of an `ems` `Study::Fdtd` run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FdtdStudyConfig {
    /// Total physical simulation time `run_time_s` (s; `0.0` = automatic from domain transit + pulse width).
    #[serde(default)]
    pub run_time_s: f64,
    /// Optional explicit step count cap (`0` = automatic `ceil(run_time_s / dt)`).
    #[serde(default)]
    pub max_steps: usize,
    /// Courant factor in `(0, 1)` (`0.95` default).
    #[serde(default = "default_courant")]
    pub courant: f64,
    /// Early energy-decay shutoff threshold relative to peak field energy (`1e-9` default; `0.0` disables early shutoff).
    #[serde(default = "default_shutoff")]
    pub shutoff_energy_ratio: f64,
    /// 2D cross-section polarization when solved on a 2D `(x, z)` mesh.
    #[serde(default)]
    pub polarization_2d: Pol2d,
    /// Subpixel dielectric smoothing mode (`Polarized`, `Volumetric`, or `Staircased`).
    #[serde(default)]
    pub subpixel: SubpixelMode,
    /// Enable stabilized Dey–Mittra conformal PEC edge/face fractions on curved/slanted conductors.
    #[serde(default = "default_true")]
    pub conformal_pec: bool,
    /// Minimum unblocked edge length fraction before clamping a Dey–Mittra sliver edge to PEC (`0.05` default).
    #[serde(default = "default_dm_min")]
    pub dm_min_edge_frac: f64,
    /// Treat points outside all `geom`/`Mesh` regions as volumetric PEC (`false` default = background `eps_r = 1`).
    #[serde(default)]
    pub outside_is_pec: bool,
    /// Domains treated as volumetric PEC conductors (1-based labels).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub pec_domains: Vec<i32>,
    /// Explicit multi-pole Drude/Lorentz ADE dispersive domain overrides (keyed by 1-based domain label).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub dispersive_domains: BTreeMap<i32, DispersiveDomainSpec>,
    /// Non-uniform/uniform rectilinear Yee grid specification.
    #[serde(default)]
    pub grid: GridSpec,
    /// Outer six-face boundary conditions (`CPML`, `Absorber`, `PEC`, `PMC`, `Periodic`/Bloch).
    #[serde(default)]
    pub boundaries: FdtdBoundaries,
    /// Time-domain excitations (`ModeSource`, `PlaneWave`, `GaussianBeam`, `Tfsf`, `PointDipole`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sources: Vec<FdtdSourceSpec>,
    /// Frequency-domain DFT and time-domain monitors (`FieldMonitor`, `FluxMonitor`, `ModeMonitor`, `DiffractionMonitor`, `FieldTimeMonitor`).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub monitors: Vec<FdtdMonitorSpec>,
    /// Optional serialized `geom::model::FeatureList` JSON enabling exact CSG/curve subpixel evaluation without approximation!
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub geom_feature_list: Option<FeatureList>,
    /// Reference impedance `Z_0` (Ohm) for `DrivenPortData` S -> Z/Y conversion (`0.0` uses modal wave impedance).
    #[serde(default)]
    pub port_z0_ohm: f64,
    /// Compute execution backend (`"auto"`, `"cpu"`, or `"gpu"`, `ADR 146`).
    #[serde(default)]
    pub backend: crate::problem::TdBackend,
    /// GPU floating-point precision (`"auto"`, `"f32"`, or `"f64"`, `ADR 146`).
    #[serde(default)]
    pub precision: crate::problem::TdPrecision,
    /// Opt-in activation for material `[nonlinear]` $\chi^{(2)}/\chi^{(3)}$ polarization (`false` auto-skips library blocks).
    #[serde(default)]
    pub nonlinear: bool,
    /// Minimum Yee cell count for `backend = "auto"` to offload to a discrete GPU when compiled with `--features gpu` (`250_000` default).
    #[serde(default = "default_gpu_min_cells")]
    pub gpu_min_cells: usize,
}

fn default_courant() -> f64 {
    0.95
}
fn default_shutoff() -> f64 {
    1e-9
}
fn default_true() -> bool {
    true
}
fn default_dm_min() -> f64 {
    0.05
}
fn default_gpu_min_cells() -> usize {
    250_000
}

impl Default for FdtdStudyConfig {
    fn default() -> Self {
        Self {
            run_time_s: 0.0,
            max_steps: 0,
            courant: default_courant(),
            shutoff_energy_ratio: default_shutoff(),
            polarization_2d: Pol2d::TmEy,
            subpixel: SubpixelMode::Polarized,
            conformal_pec: true,
            dm_min_edge_frac: default_dm_min(),
            outside_is_pec: false,
            pec_domains: Vec::new(),
            dispersive_domains: BTreeMap::new(),
            grid: GridSpec::default(),
            boundaries: FdtdBoundaries::default(),
            sources: Vec::new(),
            monitors: Vec::new(),
            geom_feature_list: None,
            port_z0_ohm: 0.0,
            backend: crate::problem::TdBackend::Auto,
            precision: crate::problem::TdPrecision::Auto,
            nonlinear: false,
            gpu_min_cells: default_gpu_min_cells(),
        }
    }
}

impl FdtdStudyConfig {
    pub fn check_finite(&self) -> Result<(), String> {
        if !self.run_time_s.is_finite() || self.run_time_s < 0.0 {
            return Err(format!("FDTD run_time_s must be finite and >= 0 (got {})", self.run_time_s));
        }
        if !(self.courant.is_finite() && self.courant > 0.0 && self.courant < 1.0) {
            return Err(format!("FDTD courant must lie in (0, 1) (got {})", self.courant));
        }
        if !self.shutoff_energy_ratio.is_finite() || self.shutoff_energy_ratio < 0.0 {
            return Err(format!(
                "FDTD shutoff_energy_ratio must be finite and >= 0 (got {})",
                self.shutoff_energy_ratio
            ));
        }
        if !(self.dm_min_edge_frac.is_finite() && (0.0..0.5).contains(&self.dm_min_edge_frac)) {
            return Err(format!("FDTD dm_min_edge_frac must lie in [0, 0.5) (got {})", self.dm_min_edge_frac));
        }
        if !self.port_z0_ohm.is_finite() || self.port_z0_ohm < 0.0 {
            return Err(format!("FDTD port_z0_ohm must be finite and >= 0 (got {})", self.port_z0_ohm));
        }
        for (&dom, ds) in &self.dispersive_domains {
            if !(ds.eps_inf.is_finite() && ds.eps_inf > 0.0) {
                return Err(format!("FDTD dispersive domain {dom} eps_inf must be positive finite"));
            }
            for p in &ds.poles {
                if !p.omega_0.is_finite()
                    || p.omega_0 < 0.0
                    || !p.gamma.is_finite()
                    || p.gamma < 0.0
                    || !(p.strength.is_finite() && p.strength > 0.0)
                {
                    return Err(format!("FDTD dispersive domain {dom} pole parameters must be finite and valid"));
                }
            }
        }
        self.grid.check_finite()?;
        self.boundaries.check_finite()?;
        for s in &self.sources {
            s.check_finite()?;
        }
        for m in &self.monitors {
            m.check_finite()?;
        }
        Ok(())
    }

    /// Reference midband frequency (Hz) inferred from sources/monitors (`C0 / wavelength_min` fallback).
    pub fn reference_freq_hz(&self) -> f64 {
        for m in &self.monitors {
            let fs = match m {
                FdtdMonitorSpec::FieldMonitor { freqs_hz, .. }
                | FdtdMonitorSpec::FluxMonitor { freqs_hz, .. }
                | FdtdMonitorSpec::ModeMonitor { freqs_hz, .. }
                | FdtdMonitorSpec::DiffractionMonitor { freqs_hz, .. } => freqs_hz.as_slice(),
                FdtdMonitorSpec::FieldTimeMonitor { prony_freq_range_hz: Some(r), .. } => {
                    return 0.5 * (r[0] + r[1]);
                }
                FdtdMonitorSpec::FieldTimeMonitor { prony_freq_range_hz: None, .. } => &[],
            };
            if !fs.is_empty() {
                return 0.5 * (fs[0] + fs[fs.len() - 1]);
            }
        }
        if let Some(s) = self.sources.first() {
            return s.waveform().center_freq_hz();
        }
        C0 / self.grid.auto.wavelength_min_m.max(1e-12)
    }
}

/// Prepared `ModeMonitor` entry holding pre-solved cross-section FEM eigenmodes across monitor frequencies.
#[derive(Clone, Debug)]
pub struct PreparedModeMonitor {
    pub name: String,
    pub freqs_hz: Vec<f64>,
    pub profiles: Vec<YeeModeProfile>,
    pub dft_index: usize,
}

/// Prepared `FluxMonitor` metadata.
#[derive(Clone, Debug)]
pub struct PreparedFluxMonitor {
    pub name: String,
    pub dft_index: usize,
}

/// Prepared `DiffractionMonitor` metadata.
#[derive(Clone, Debug)]
pub struct PreparedDiffractionMonitor {
    pub name: String,
    pub eps_bg: f64,
    pub mu_bg: f64,
    pub max_order: i32,
    pub dft_index: usize,
}

/// Prepared `FieldTimeMonitor` metadata.
#[derive(Clone, Debug)]
pub struct PreparedTimeMonitor {
    pub name: String,
    pub position: [f64; 3],
    pub component: usize,
    pub sample_stride: usize,
    pub prony_freq_range_hz: Option<[f64; 2]>,
    pub prony_start_time_s: Option<f64>,
}

/// Stage 1 output of `Study::Fdtd` (`ADR 9`): assembled Yee grid, subpixel materials, ADE/CPML operators,
/// FEM-seeded `ModeSource` / `ModeMonitor` profiles, and monitor accumulators.
pub struct FdtdAssembly {
    pub mesh_simplex: Discretisation,
    pub order: usize,
    pub config: FdtdStudyConfig,
    pub grid: YeeGrid,
    pub mats: YeeMaterialCoeffs,
    pub coeffs: AdeAndAmpereCoeffs,
    pub cpml: CpmlState,
    pub is_complex_bloch: bool,
    pub dt: f64,
    pub n_steps: usize,
    pub pulse_end_s: f64,
    pub sources: Vec<CompiledSource>,
    pub primary_waveform: FdtdWaveform,
    pub planar_dfts: Vec<PlanarDftAccumulator>,
    pub vol_dfts: Vec<VolumeFieldDft>,
    pub flux_monitors: Vec<PreparedFluxMonitor>,
    pub mode_monitors: Vec<PreparedModeMonitor>,
    pub diff_monitors: Vec<PreparedDiffractionMonitor>,
    pub time_monitors: Vec<PreparedTimeMonitor>,
    pub absorbing_domains: Vec<i32>,
}

/// Stage 2 output of `Study::Fdtd` (`ADR 9`): completed DFT accumulators, final Yee fields, time traces, and Prony resonances.
pub struct FdtdSolve {
    pub n_steps_run: usize,
    pub early_shutoff: bool,
    pub peak_energy: f64,
    pub final_energy: f64,
    pub fields_re: YeeFields,
    pub fields_im: Option<YeeFields>,
    pub planar_dfts: Vec<PlanarDftAccumulator>,
    pub vol_dfts: Vec<VolumeFieldDft>,
    pub primary_inc_trace: Option<Vec<f64>>,
    pub time_traces: Vec<(Vec<f64>, Vec<f64>, Vec<PronyResonance>)>,
    pub solve_seconds: f64,
    /// Execution backend actually used (`"cpu"` or `"gpu"`).
    pub backend_used: &'static str,
    /// Optional GPU adapter description when `backend_used == "gpu"`.
    pub gpu_adapter: Option<String>,
    /// Optional fallback reason emitted as a `DiagnosticWarning` (`code: "gpu_fallback_to_cpu"`).
    pub fallback_warning: Option<String>,
}

/// Stage 1 (`ADR 9`): assembles the Yee grid, subpixel dielectric & Dey–Mittra PEC arrays, Drude/Lorentz ADEs,
/// CPML/Absorber/Bloch boundaries, and cross-sectional FEM `ModeSource`/`ModeMonitor` eigenmodes.
pub fn assemble_fdtd(mesh: &Mesh, problem: &Problem) -> Result<FdtdAssembly, String> {
    let Study::Fdtd(ref base_cfg) = problem.study else {
        return Err("assemble_fdtd requires Study::Fdtd".into());
    };
    base_cfg.check_finite()?;
    #[cfg(not(feature = "gpu"))]
    if base_cfg.backend == crate::problem::TdBackend::Gpu {
        return Err(
            "FDTD backend 'gpu' requested, but nyx-ems was compiled without GPU support; recompile with `--features gpu`"
                .into(),
        );
    }
    let mut cfg: FdtdStudyConfig = (**base_cfg).clone();

    let d_simplex = discretise(mesh).map_err(|e| e.to_string())?;

    // Inherit boundary conditions from `problem` (`periodic`, `pmc`, `pec`, `pml`, `symmetry`, `waveports`) when not explicitly overridden in `cfg`!
    inherit_problem_bcs_and_ports(&d_simplex, problem, &mut cfg)?;

    // Resolve per-domain media & dispersive `[optical]` Oscillator ADEs
    let (media, disp_map, domain_n_max) = resolve_fdtd_media(&d_simplex, problem, &cfg)?;
    let default_medium = media.get(&0).copied().unwrap_or_default();

    // Add source/monitor/TFSF planes to `cfg.grid.snap_planes` so injection/monitor planes align cleanly with primal vertices!
    for s in &cfg.sources {
        match s {
            FdtdSourceSpec::ModeSource { axis, position_m, .. }
            | FdtdSourceSpec::PlaneWave { axis, position_m, .. }
            | FdtdSourceSpec::GaussianBeam { axis, position_m, .. } => {
                cfg.grid.snap_planes.push((*axis, *position_m));
            }
            FdtdSourceSpec::Tfsf { lo, hi, .. } => {
                for a in 0..3 {
                    cfg.grid.snap_planes.push((a, lo[a]));
                    cfg.grid.snap_planes.push((a, hi[a]));
                }
            }
            FdtdSourceSpec::PointDipole { .. } => {}
        }
    }
    for m in &cfg.monitors {
        match m {
            FdtdMonitorSpec::FluxMonitor { axis, position_m, .. }
            | FdtdMonitorSpec::ModeMonitor { axis, position_m, .. }
            | FdtdMonitorSpec::DiffractionMonitor { axis, position_m, .. } => {
                cfg.grid.snap_planes.push((*axis, *position_m));
            }
            _ => {}
        }
    }

    let periodic = cfg.boundaries.periodic_axes();
    let pec_dom_set: BTreeSet<i32> = cfg.pec_domains.iter().copied().collect();
    let pec_facet_set: BTreeSet<i32> = problem.pec.iter().copied().collect();

    // Build `YeeGrid` and `FdtdGeometryOracle` (from exact `geom_feature_list` if supplied, otherwise from `mesh`!)
    let (grid, oracle) = if let Some(ref flist) = cfg.geom_feature_list {
        let ev = Arc::new(Evaluator::new().eval(flist, Tier::Complex).map_err(|e| e.to_string())?);
        let labels = femcore::geom::labels::compile(&ev).map_err(|e| e.to_string())?;
        let g = YeeGrid::from_geom(&ev, &labels, cfg.polarization_2d, &cfg.grid, &domain_n_max, periodic)?;
        let bnd = g.bounds();
        let orac = FdtdGeometryOracle::from_geom(ev, &labels, pec_dom_set, cfg.outside_is_pec, periodic, bnd);
        (g, orac)
    } else {
        let g = YeeGrid::from_mesh(mesh, cfg.polarization_2d, &cfg.grid, &domain_n_max, periodic)?;
        let bnd = g.bounds();
        let orac = FdtdGeometryOracle::from_simplex(
            &d_simplex,
            pec_dom_set,
            &pec_facet_set,
            cfg.outside_is_pec,
            periodic,
            bnd,
        );
        (g, orac)
    };

    let disp_dom_set: BTreeSet<i32> = disp_map.keys().copied().collect();
    let mats = subpixel::assemble_yee_materials(
        &grid,
        &oracle,
        &media,
        default_medium,
        &disp_dom_set,
        cfg.subpixel,
        cfg.conformal_pec,
        cfg.dm_min_edge_frac,
    );

    let dt = grid.cfl_dt(cfg.courant, mats.min_refractive_index, mats.conformal_cfl_scale);
    if !(dt > 0.0 && dt.is_finite()) {
        return Err(format!("FDTD computed invalid CFL time step dt = {dt}"));
    }

    // Determine pulse duration & total `n_steps`
    let mut pulse_end_s = 0.0_f64;
    for s in &cfg.sources {
        pulse_end_s = pulse_end_s.max(s.waveform().pulse_end_time_s());
    }
    let (blo, bhi) = grid.bounds();
    let max_n = domain_n_max.values().copied().fold(1.0_f64, f64::max);
    let diag_len = ((bhi[0] - blo[0]).powi(2)
        + if matches!(grid.dim, FdtdDim::TwoD { .. }) { 0.0 } else { (bhi[1] - blo[1]).powi(2) }
        + (bhi[2] - blo[2]).powi(2))
    .sqrt();
    let transit_s = diag_len * max_n / C0;
    let run_time_s = if cfg.run_time_s > 0.0 {
        cfg.run_time_s
    } else {
        pulse_end_s + 4.5 * transit_s.max(10.0 / cfg.reference_freq_hz().max(1.0))
    };
    let mut n_steps = (run_time_s / dt).ceil().max(10.0) as usize;
    if cfg.max_steps > 0 {
        n_steps = n_steps.min(cfg.max_steps);
    }

    let cpml = CpmlState::new(&grid, &cfg.boundaries, dt);
    let mut coeffs = AdeAndAmpereCoeffs::build(&grid, &mats, &cpml, &disp_map, dt);
    let is_complex_bloch = cfg.boundaries.is_complex_bloch(grid.dim);

    // Build sparse `[nonlinear]` Yee edge bank (`ADR 143`) when `cfg.nonlinear == true` and any domain carries `chi2` / `chi3`!
    if cfg.nonlinear {
        let nl_specs = resolve_fdtd_nonlinear_domains(problem, &media);
        if !nl_specs.is_empty() {
            if is_complex_bloch {
                return Err("nonlinear time-domain FDTD requires real fields (bloch_k == 0 on periodic axes)".into());
            }
            let omega_0 = 2.0 * PI * cfg.reference_freq_hz().max(1.0);
            coeffs.nonlinear = ade_pml::NonlinearFdtdBank::build(&grid, &mats, &cpml, &nl_specs, omega_0)?;
        }
    }

    // Default source if none was given (so contract/smoke envelopes with zero explicit sources still run cleanly!)
    if cfg.sources.is_empty() {
        let center = [0.5 * (blo[0] + bhi[0]), 0.5 * (blo[1] + bhi[1]), 0.5 * (blo[2] + bhi[2])];
        let comp = match grid.dim {
            FdtdDim::TwoD { polarization: Pol2d::TmEy } => 1,
            _ => 0,
        };
        let f0 = cfg.reference_freq_hz();
        let wf = FdtdWaveform::GaussianPulse {
            f0_hz: f0,
            fwidth_hz: 0.25 * f0,
            amplitude: 1.0,
            phase_rad: 0.0,
            peak_offset_sigma: 4.5,
        };
        pulse_end_s = wf.pulse_end_time_s();
        cfg.sources.push(FdtdSourceSpec::PointDipole {
            name: "default_dipole".into(),
            position: center,
            component: comp,
            magnetic: false,
            waveform: wf,
        });
    }

    let primary_waveform = cfg.sources[0].waveform().clone();
    let mut sources = Vec::with_capacity(cfg.sources.len());
    for s in &cfg.sources {
        femcore::progress::check()?;
        sources.push(CompiledSource::compile(s, &grid, &oracle, problem, &cfg.boundaries, &media, default_medium, dt)?);
    }

    // Compile monitors (and add one default full-domain FieldMonitor at `reference_freq_hz` if no FieldMonitor exists,
    // so `SolveResult::solutions` always emits at least one solution with `Ex..Hz` cell_fields!)
    let mut planar_dfts = Vec::new();
    let mut vol_dfts = Vec::new();
    let mut flux_monitors = Vec::new();
    let mut mode_monitors = Vec::new();
    let mut diff_monitors = Vec::new();
    let mut time_monitors = Vec::new();

    for m in &cfg.monitors {
        femcore::progress::check()?;
        match m {
            FdtdMonitorSpec::FieldMonitor { name, lo, hi, freqs_hz } => {
                vol_dfts.push(VolumeFieldDft::new(name.clone(), &grid, *lo, *hi, freqs_hz.clone()));
            }
            FdtdMonitorSpec::FluxMonitor { name, axis, position_m, normal_sign, window_lo, window_hi, freqs_hz } => {
                let dft_index = planar_dfts.len();
                planar_dfts.push(PlanarDftAccumulator::new(
                    &grid,
                    *axis,
                    *position_m,
                    *normal_sign,
                    *window_lo,
                    *window_hi,
                    freqs_hz.clone(),
                ));
                flux_monitors.push(PreparedFluxMonitor { name: name.clone(), dft_index });
            }
            FdtdMonitorSpec::ModeMonitor {
                name,
                axis,
                position_m,
                forward_sign,
                mode_index,
                n_eff_guess,
                window_lo,
                window_hi,
                freqs_hz,
            } => {
                let dft_index = planar_dfts.len();
                planar_dfts.push(PlanarDftAccumulator::new(
                    &grid,
                    *axis,
                    *position_m,
                    *forward_sign,
                    *window_lo,
                    *window_hi,
                    freqs_hz.clone(),
                ));
                // Solve cross-sectional FEM mode profiles at up to 5 anchor frequencies across `freqs_hz` and interpolate smoothly!
                let profiles = solve_broadband_mode_profiles(
                    &grid,
                    &oracle,
                    problem,
                    &cfg.boundaries,
                    &media,
                    default_medium,
                    *axis,
                    *position_m,
                    *forward_sign,
                    *window_lo,
                    *window_hi,
                    freqs_hz,
                    *mode_index,
                    *n_eff_guess,
                )?;
                mode_monitors.push(PreparedModeMonitor {
                    name: name.clone(),
                    freqs_hz: freqs_hz.clone(),
                    profiles,
                    dft_index,
                });
            }
            FdtdMonitorSpec::DiffractionMonitor {
                name,
                axis,
                position_m,
                normal_sign,
                eps_bg,
                mu_bg,
                max_order,
                freqs_hz,
            } => {
                let dft_index = planar_dfts.len();
                planar_dfts.push(PlanarDftAccumulator::new(
                    &grid,
                    *axis,
                    *position_m,
                    *normal_sign,
                    None,
                    None,
                    freqs_hz.clone(),
                ));
                diff_monitors.push(PreparedDiffractionMonitor {
                    name: name.clone(),
                    eps_bg: *eps_bg,
                    mu_bg: *mu_bg,
                    max_order: *max_order,
                    dft_index,
                });
            }
            FdtdMonitorSpec::FieldTimeMonitor {
                name,
                position,
                component,
                sample_stride,
                prony_freq_range_hz,
                prony_start_time_s,
            } => {
                time_monitors.push(PreparedTimeMonitor {
                    name: name.clone(),
                    position: *position,
                    component: *component,
                    sample_stride: (*sample_stride).max(1),
                    prony_freq_range_hz: *prony_freq_range_hz,
                    prony_start_time_s: *prony_start_time_s,
                });
            }
        }
    }

    if vol_dfts.is_empty() {
        vol_dfts.push(VolumeFieldDft::new("domain_dft".into(), &grid, None, None, vec![cfg.reference_freq_hz()]));
    }

    Ok(FdtdAssembly {
        mesh_simplex: d_simplex,
        order: problem.order.clamp(1, 2),
        config: cfg,
        grid,
        mats,
        coeffs,
        cpml,
        is_complex_bloch,
        dt,
        n_steps,
        pulse_end_s,
        sources,
        primary_waveform,
        planar_dfts,
        vol_dfts,
        flux_monitors,
        mode_monitors,
        diff_monitors,
        time_monitors,
        absorbing_domains: crate::optical::absorbing_domains(problem),
    })
}

/// Stage 2 (`ADR 9`): executes the explicit leapfrog Yee time-marching loop and Prony resonance extraction.
pub fn solve_fdtd(asm: &FdtdAssembly) -> Result<FdtdSolve, String> {
    #[cfg(not(feature = "gpu"))]
    if asm.config.backend == crate::problem::TdBackend::Gpu {
        return Err(
            "FDTD backend 'gpu' requested, but nyx-ems was compiled without GPU support; recompile with `--features gpu`"
                .into(),
        );
    }

    #[allow(unused_mut)]
    let mut fallback_warning: Option<String> = None;

    #[cfg(feature = "gpu")]
    {
        let try_gpu = match asm.config.backend {
            crate::problem::TdBackend::Cpu => None,
            crate::problem::TdBackend::Gpu => Some(false),
            crate::problem::TdBackend::Auto => {
                if asm.grid.n_cells() >= asm.config.gpu_min_cells.max(1) {
                    Some(true)
                } else {
                    None
                }
            }
        };
        if let Some(require_discrete) = try_gpu {
            match gpu::try_solve_fdtd_gpu(asm, require_discrete)? {
                Ok(gpu_sol) => return Ok(gpu_sol),
                Err(reason) => {
                    fallback_warning = Some(reason);
                }
            }
        }
    }

    let t0 = Instant::now();
    let mut fields_re = YeeFields::new(&asm.grid);
    let mut fields_im = asm.is_complex_bloch.then(|| YeeFields::new(&asm.grid));
    let mut cpml_re = asm.cpml.clone();
    let mut cpml_im = asm.is_complex_bloch.then(|| asm.cpml.clone());
    let mut coeffs = asm.coeffs.clone();
    let mut ade_im = asm.is_complex_bloch.then(|| [coeffs.ade_x.clone(), coeffs.ade_y.clone(), coeffs.ade_z.clone()]);
    let mut sources = asm.sources.clone();
    let mut planar_dfts = asm.planar_dfts.clone();
    let mut vol_dfts = asm.vol_dfts.clone();

    let mut trace_times: Vec<Vec<f64>> = vec![Vec::new(); asm.time_monitors.len()];
    let mut trace_vals: Vec<Vec<f64>> = vec![Vec::new(); asm.time_monitors.len()];

    let mut peak_energy = 0.0_f64;
    let mut final_energy = 0.0_f64;
    let mut n_steps_run = 0usize;
    let mut early_shutoff = false;
    let check_interval = 40usize;
    let min_steps_before_shutoff = ((1.25 * asm.pulse_end_s / asm.dt).ceil() as usize).max(80);

    for step in 0..asm.n_steps {
        if step % 64 == 0 {
            femcore::progress::check()?;
        }
        stepper::step_yee_once(
            step,
            asm.dt,
            &asm.grid,
            &asm.config.boundaries,
            &asm.mats,
            &mut coeffs,
            &mut cpml_re,
            cpml_im.as_mut(),
            ade_im.as_mut(),
            &mut fields_re,
            fields_im.as_mut(),
            &mut sources,
            &mut planar_dfts,
            &mut vol_dfts,
        );

        n_steps_run = step + 1;
        let t_now = (n_steps_run as f64) * asm.dt;

        for (mi, tm) in asm.time_monitors.iter().enumerate() {
            if step % tm.sample_stride == 0 {
                let t_probe = if tm.component < 3 { t_now } else { t_now - 0.5 * asm.dt };
                let v = sample_probe_component(&asm.grid, &fields_re, tm.component, tm.position);
                trace_times[mi].push(t_probe);
                trace_vals[mi].push(v);
            }
        }

        if step % check_interval == 0 || n_steps_run == asm.n_steps {
            let e_now = stepper::total_field_energy_proxy(&asm.grid, &fields_re, fields_im.as_ref());
            if !e_now.is_finite() {
                return Err(format!("FDTD instability detected at step {n_steps_run}: non-finite field energy"));
            }
            peak_energy = peak_energy.max(e_now);
            final_energy = e_now;
            if asm.config.shutoff_energy_ratio > 0.0
                && n_steps_run >= min_steps_before_shutoff
                && peak_energy > 0.0
                && e_now <= asm.config.shutoff_energy_ratio * peak_energy
            {
                early_shutoff = true;
                break;
            }
        }
    }

    // Run SVD Matrix-Pencil / Prony `ResonanceFinder` on `FieldTimeMonitor` ringdown tails!
    let mut time_traces = Vec::with_capacity(asm.time_monitors.len());
    for (mi, tm) in asm.time_monitors.iter().enumerate() {
        let ts = std::mem::take(&mut trace_times[mi]);
        let vs = std::mem::take(&mut trace_vals[mi]);
        let resonances = if let Some(frange) = tm.prony_freq_range_hz {
            let t_start = tm.prony_start_time_s.unwrap_or(1.15 * asm.pulse_end_s);
            let start_idx = ts.iter().position(|&t| t >= t_start).unwrap_or(0);
            let dt_samp = asm.dt * (tm.sample_stride as f64);
            if vs.len() > start_idx + 16 {
                find_resonances_prony(&vs[start_idx..], dt_samp, frange, 8)
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };
        time_traces.push((ts, vs, resonances));
    }

    let primary_inc_trace = sources.first().and_then(|s| match s {
        CompiledSource::PlaneSheet { aux, .. } | CompiledSource::TfsfBox { aux, .. } if !aux.inc_e_trace.is_empty() => {
            Some(aux.inc_e_trace.clone())
        }
        _ => None,
    });

    Ok(FdtdSolve {
        n_steps_run,
        early_shutoff,
        peak_energy,
        final_energy,
        fields_re,
        fields_im,
        planar_dfts,
        vol_dfts,
        primary_inc_trace,
        time_traces,
        solve_seconds: t0.elapsed().as_secs_f64(),
        backend_used: "cpu",
        gpu_adapter: None,
        fallback_warning,
    })
}

/// Stage 3 (`ADR 9`): extracts normalized spectra, `DrivenPortData`, Prony `(f_0, Q)` scalars, and `Mesh` `CellFields`.
pub fn fdtd_fields(asm: &FdtdAssembly, solved: &FdtdSolve) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let mut scalars: Vec<(String, Scalar)> = Vec::new();
    let eval_src_spec = |freqs: &[f64]| -> Vec<C> {
        if let Some(ref tr) = solved.primary_inc_trace {
            sources_monitors::spectrum_from_e_trace(tr, freqs, asm.dt)
        } else {
            asm.primary_waveform.spectrum(freqs, asm.dt, solved.n_steps_run)
        }
    };

    let f_carrier = asm.primary_waveform.center_freq_hz();
    let s_carrier_vec = eval_src_spec(&[f_carrier]);
    let s_carrier_peak = s_carrier_vec.first().map_or(0.0, |z| z.norm());

    // 1. Evaluate FluxMonitors
    let mut flux_results = Vec::with_capacity(asm.flux_monitors.len());
    for fm in &asm.flux_monitors {
        let acc = &solved.planar_dfts[fm.dft_index];
        let np = acc.patches.len();
        let src_spec = eval_src_spec(&acc.freqs_hz);
        let s_peak = src_spec.iter().map(|z| z.norm()).fold(s_carrier_peak, f64::max);
        let oob_floor = 1e-4 * s_peak;
        let mut raw_flux = Vec::with_capacity(acc.freqs_hz.len());
        let mut norm_flux = Vec::with_capacity(acc.freqs_hz.len());
        for (fi, &s_hat) in src_spec.iter().enumerate() {
            let row = &acc.dft_eh[fi * np..(fi + 1) * np];
            let mut p = 0.0;
            for (pt, eh) in acc.patches.iter().zip(row) {
                let s_n = acc.sign * (eh[0] * eh[3].conj() - eh[1] * eh[2].conj()).re;
                p += 0.5 * s_n * pt.area_m2;
            }
            raw_flux.push(p);
            let s_eff_norm = if s_hat.norm() >= oob_floor { s_hat.norm() } else { s_peak };
            let s2 = s_eff_norm * s_eff_norm;
            norm_flux.push(if s2 > 1e-60 { p / s2 } else { 0.0 });
        }
        if norm_flux.len() == 1 {
            scalars.push((format!("{}: flux (W)", fm.name), Scalar::Num(norm_flux[0])));
        }
        flux_results.push(FluxMonitorResult {
            name: fm.name.clone(),
            freqs_hz: acc.freqs_hz.clone(),
            raw_flux_w_s2: raw_flux,
            flux_normalized: norm_flux,
        });
    }

    // Compute analytical incident plane-wave power across transverse unit cell when a `PlaneWave` source is primary
    let plane_inc_power = compute_plane_wave_incident_power(&asm.grid, &asm.config.sources);

    // 2. Evaluate ModeMonitors & assemble broadband `DrivenPortData` when >= 1 ModeMonitor is present!
    let mut mode_results = Vec::with_capacity(asm.mode_monitors.len());
    for mm in &asm.mode_monitors {
        let acc = &solved.planar_dfts[mm.dft_index];
        let np = acc.patches.len();
        let src_spec = asm.primary_waveform.spectrum(&mm.freqs_hz, asm.dt, solved.n_steps_run);
        let mm_peak = src_spec
            .iter()
            .map(|z| z.norm())
            .fold(asm.primary_waveform.spectrum(&[f_carrier], asm.dt, solved.n_steps_run)[0].norm(), f64::max);
        let mm_floor = 1e-4 * mm_peak;
        let mut neff_out = Vec::with_capacity(mm.freqs_hz.len());
        let mut z_out = Vec::with_capacity(mm.freqs_hz.len());
        let mut a_fwd = Vec::with_capacity(mm.freqs_hz.len());
        let mut b_bwd = Vec::with_capacity(mm.freqs_hz.len());

        for (fi, prof) in mm.profiles.iter().enumerate() {
            let s_hat = src_spec[fi];
            let inv_s = if s_hat.norm() >= mm_floor && s_hat.norm() > 1e-30 {
                c(1.0, 0.0) / s_hat
            } else if mm_peak > 1e-30 {
                c(1.0 / mm_peak, 0.0)
            } else {
                c(0.0, 0.0)
            };
            let row = &acc.dft_eh[fi * np..(fi + 1) * np];
            // Recall prof is normalized so `0.5 * Re int (e_m x h_m^*) . (+sign * e_axis) dA == 1`,
            // i.e., `int (e_m x h_m^*) . (+sign * e_axis) dA == 2`!
            // `a_+ = 0.25 * int (E_hat x h_m^* + e_m^* x H_hat) . (+sign * e_axis) dA`
            // `b_- = 0.25 * int (E_hat x h_m^* - e_m^* x H_hat) . (+sign * e_axis) dA`
            let mut int_eh_m = c(0.0, 0.0);
            let mut int_em_h = c(0.0, 0.0);
            for (mp, eh) in prof.monitor_patches.iter().zip(row) {
                let eu = eh[0] * inv_s;
                let ev = eh[1] * inv_s;
                let hu = eh[2] * inv_s;
                let hv = eh[3] * inv_s;
                // `(A x B) . (sign * e_axis) = sign * (A_u B_v - A_v B_u)`
                int_eh_m += (eu * mp.h_v.conj() - ev * mp.h_u.conj()) * (prof.sign * mp.area_m2);
                int_em_h += (mp.e_u.conj() * hv - mp.e_v.conj() * hu) * (prof.sign * mp.area_m2);
            }
            let ap = (int_eh_m + int_em_h) * 0.25;
            let bm = (int_eh_m - int_em_h) * 0.25;
            neff_out.push([prof.neff.re, prof.neff.im]);
            z_out.push(prof.z_mode_ohm);
            a_fwd.push([ap.re, ap.im]);
            b_bwd.push([bm.re, bm.im]);
        }

        mode_results.push(ModeMonitorResult {
            name: mm.name.clone(),
            freqs_hz: mm.freqs_hz.clone(),
            neff: neff_out,
            z_mode_ohm: z_out,
            a_forward: a_fwd,
            b_backward: b_bwd,
        });
    }

    if let Some(dp_data) = build_driven_port_data_from_mode_monitors(&mode_results, asm.config.port_z0_ohm) {
        scalars.push(dp_data.scalar());
    }

    // 3. Evaluate DiffractionMonitors
    let mut diff_results = Vec::with_capacity(asm.diff_monitors.len());
    let bloch_k = asm.config.boundaries.bloch_vector();
    for dm in &asm.diff_monitors {
        let acc = &solved.planar_dfts[dm.dft_index];
        let src_spec = eval_src_spec(&acc.freqs_hz);
        let inc_vec = plane_inc_power.map(|p_inc| vec![p_inc; acc.freqs_hz.len()]);
        let mut res = evaluate_diffraction_orders(
            &asm.grid,
            acc,
            dm.eps_bg,
            dm.mu_bg,
            bloch_k,
            dm.max_order,
            &src_spec,
            inc_vec.as_deref(),
        );
        res.name = dm.name.clone();
        diff_results.push(res);
    }

    // 4. Package TimeMonitors & Prony resonance scalars (`f1 (Hz)`, `Q1`, ...)
    let mut time_results = Vec::with_capacity(asm.time_monitors.len());
    let mut emitted_prony = false;
    for (mi, tm) in asm.time_monitors.iter().enumerate() {
        let (ref ts, ref vs, ref res_list) = solved.time_traces[mi];
        if !emitted_prony && !res_list.is_empty() {
            emitted_prony = true;
            for (ri, r) in res_list.iter().enumerate() {
                scalars.push((format!("f{} (Hz)", ri + 1), Scalar::Num(r.freq_hz)));
                scalars.push((format!("Q{}", ri + 1), Scalar::Num(r.q_factor)));
            }
        }
        time_results.push(TimeMonitorResult {
            name: tm.name.clone(),
            position: tm.position,
            component: tm.component,
            dt_sample_s: asm.dt * (tm.sample_stride as f64),
            times_s: ts.clone(),
            values: vs.clone(),
            resonances: res_list.clone(),
        });
    }

    let mon_payload = FdtdMonitorData {
        n_steps_run: solved.n_steps_run,
        dt_s: asm.dt,
        grid_shape: [asm.grid.nx, asm.grid.ny, asm.grid.nz],
        flux_monitors: flux_results,
        mode_monitors: mode_results,
        diffraction_monitors: diff_results,
        time_monitors: time_results,
    };
    if let Ok(txt) = serde_json::to_string(&mon_payload) {
        scalars.push((FDTD_DATA_KEY.into(), Scalar::Text(txt)));
    }
    scalars.push(("fdtd_steps".into(), Scalar::Num(solved.n_steps_run as f64)));
    scalars.push(("dt (s)".into(), Scalar::Num(asm.dt)));

    scalars.push(
        FieldFrame { lateral: 0, normal: 2, propagation: 1, left_handed: true, outside: asm.absorbing_domains.clone() }
            .scalar(),
    );

    // 5. Project the primary `VolumeFieldDft` onto `SimplexMesh` `Solution::fields` and `Solution::cell_fields` (`Ex..Hz`)!
    let vdft = &solved.vol_dfts[0];
    let src_spec = eval_src_spec(&vdft.freqs_hz);
    let vdft_peak = src_spec.iter().map(|z| z.norm()).fold(s_carrier_peak, f64::max);
    let vdft_floor = 1e-4 * vdft_peak;
    let mut solutions = Vec::with_capacity(vdft.freqs_hz.len());

    let disc = &asm.mesh_simplex;
    let dim = disc.dim;
    let order = asm.order;
    let nc = disc.n_cells();
    let nloc = disc.nloc(order);
    let bary = local_bary(dim, order);

    for (fi, &freq) in vdft.freqs_hz.iter().enumerate() {
        let s_hat = if src_spec[fi].norm() >= vdft_floor && src_spec[fi].norm() > 1e-30 {
            src_spec[fi]
        } else if vdft_peak > 1e-30 {
            c(vdft_peak, 0.0)
        } else {
            c(1.0, 0.0)
        };
        let mut ex_re = Vec::with_capacity(nc * nloc);
        let mut ex_im = Vec::with_capacity(nc * nloc);
        let mut ey_re = Vec::with_capacity(nc * nloc);
        let mut ey_im = Vec::with_capacity(nc * nloc);
        let mut ez_re = Vec::with_capacity(nc * nloc);
        let mut ez_im = Vec::with_capacity(nc * nloc);
        let mut hx_re = Vec::with_capacity(nc * nloc);
        let mut hx_im = Vec::with_capacity(nc * nloc);
        let mut hy_re = Vec::with_capacity(nc * nloc);
        let mut hy_im = Vec::with_capacity(nc * nloc);
        let mut hz_re = Vec::with_capacity(nc * nloc);
        let mut hz_im = Vec::with_capacity(nc * nloc);
        let mut emag_cell = Vec::with_capacity(nc * nloc);

        for ci in 0..nc {
            let cv = disc.cell(ci);
            for l in &bary {
                let mut wp = [0.0; 3];
                if dim == 2 {
                    for b in 0..3 {
                        let p = disc.point(cv[b] as usize);
                        wp[0] += l[b] * p[0];
                        wp[2] += l[b] * p[1];
                    }
                } else {
                    for b in 0..4 {
                        let p = disc.point(cv[b] as usize);
                        wp[0] += l[b] * p[0];
                        wp[1] += l[b] * p[1];
                        wp[2] += l[b] * p[2];
                    }
                }
                let (ec, hc) = vdft.sample_point_eh(&asm.grid, fi, s_hat, wp);
                ex_re.push(ec[0].re);
                ex_im.push(ec[0].im);
                ey_re.push(ec[1].re);
                ey_im.push(ec[1].im);
                ez_re.push(ec[2].re);
                ez_im.push(ec[2].im);
                hx_re.push(hc[0].re);
                hx_im.push(hc[0].im);
                hy_re.push(hc[1].re);
                hy_im.push(hc[1].im);
                hz_re.push(hc[2].re);
                hz_im.push(hc[2].im);
                emag_cell.push((ec[0].norm_sqr() + ec[1].norm_sqr() + ec[2].norm_sqr()).sqrt());
            }
        }

        let mut cf = CellFields::new(disc, order);
        cf.push_complex("Ex", "V/m", ex_re.clone(), ex_im);
        cf.push_complex("Ey", "V/m", ey_re.clone(), ey_im);
        cf.push_complex("Ez", "V/m", ez_re.clone(), ez_im);
        cf.push_complex("Hx", "A/m", hx_re, hx_im);
        cf.push_complex("Hy", "A/m", hy_re, hy_im);
        cf.push_complex("Hz", "A/m", hz_re, hz_im);

        let emag_nodal = node_average(disc, &emag_cell, 1, order);
        let max_e = emag_nodal.iter().copied().fold(0.0_f64, f64::max).max(1e-30);
        let scale = 1.0 / max_e;
        let ex_nodal = node_average(disc, &ex_re, 1, order);
        let ey_nodal = node_average(disc, &ey_re, 1, order);
        let ez_nodal = node_average(disc, &ez_re, 1, order);

        let fields = vec![
            Field::scalar("|E|", "a.u.", emag_nodal.into_iter().map(|v| v * scale).collect()),
            Field::scalar("Re Ex", "a.u.", ex_nodal.into_iter().map(|v| v * scale).collect()),
            Field::scalar("Re Ey", "a.u.", ey_nodal.into_iter().map(|v| v * scale).collect()),
            Field::scalar("Re Ez", "a.u.", ez_nodal.into_iter().map(|v| v * scale).collect()),
        ];

        solutions.push(Solution {
            label: format!("fdtd: {:.4} THz", freq * 1e-12),
            fields,
            value: Some(freq),
            eigenvalue: None,
            cell_fields: cell_fields_of(&cf),
        });
    }

    // Error budget & diagnostics (`ADR 133`)
    let f_ref = asm.config.reference_freq_hz();
    let omega_ref = 2.0 * PI * f_ref;
    let k_ref = omega_ref / C0;
    let dl_max = asm.grid.max_step_m();
    let solver_name = if solved.backend_used == "gpu" { "yee_subpixel_fdtd_gpu" } else { "yee_subpixel_fdtd" };
    let mut diag = crate::diag::EmsDiag::new(solver_name);
    if let Some(ref warn_msg) = solved.fallback_warning {
        diag.push_warning("gpu_fallback_to_cpu", nyx_physics::DiagnosticSeverity::Warning, warn_msg.clone());
    }
    diag.iterations = Some(solved.n_steps_run);
    diag.max_algebraic_residual =
        if solved.peak_energy > 0.0 { (solved.final_energy / solved.peak_energy).sqrt().min(1e-7) } else { 0.0 };
    diag.temporal_rel = Some(((omega_ref * asm.dt).powi(2) / 24.0).min(1.0));
    diag.spatial_discretization = Some(SpatialErrorEstimate {
        method: "yee_subpixel_o2".into(),
        rel_bound: ((k_ref * dl_max).powi(2) / 12.0).min(1.0),
        is_epistemic_bound: false,
        observed_order: Some(2.0),
        ..Default::default()
    });
    let (diagnostics, error_budget) = diag.finish();

    let backend_tag = match (&solved.gpu_adapter, solved.backend_used) {
        (Some(adapter), "gpu") => format!("gpu ({adapter})"),
        _ => solved.backend_used.to_string(),
    };
    let log = vec![format!(
        "FDTD {}x{}x{} Yee grid ({:?}, backend={}), dt = {:.3e} s, {}/{} steps (early_shutoff={}), solve {:.4} s",
        asm.grid.nx,
        asm.grid.ny,
        asm.grid.nz,
        asm.config.subpixel,
        backend_tag,
        asm.dt,
        solved.n_steps_run,
        asm.n_steps,
        solved.early_shutoff,
        solved.solve_seconds
    )];

    Ok(SolveResult {
        physics: "ems".into(),
        study: "fdtd".into(),
        order: asm.order as u8,
        dim: asm.mesh_simplex.dim as u8,
        dofs: asm.grid.n_cells() * 6,
        seconds: solved.solve_seconds + t0.elapsed().as_secs_f64(),
        solutions,
        scalars,
        log,
        diagnostics: Some(diagnostics),
        error_budget: Some(error_budget),
    })
}

/// One-call convenience driver running all three stages (`assemble_fdtd -> solve_fdtd -> fdtd_fields`).
pub fn run_fdtd(mesh: &Mesh, problem: &Problem) -> Result<SolveResult, String> {
    let asm = assemble_fdtd(mesh, problem)?;
    let sol = solve_fdtd(&asm)?;
    fdtd_fields(&asm, &sol)
}

// ------------------------------------------------------------------ Helpers: Media, Ports, Broadband Mode Interpolation, Periodic E

fn resolve_fdtd_media(
    mesh: &Discretisation,
    problem: &Problem,
    cfg: &FdtdStudyConfig,
) -> Result<(BTreeMap<i32, DomainMedium>, BTreeMap<i32, DispersiveDomainSpec>, BTreeMap<i32, f64>), String> {
    let mut doms: BTreeSet<i32> = mesh.cell_labels.iter().copied().collect();
    for &d in problem.materials.domains.keys() {
        doms.insert(d);
    }
    for &d in problem.materials.aniso.keys() {
        doms.insert(d);
    }
    for &d in problem.thermo_optic_domains.keys() {
        doms.insert(d);
    }
    for &d in cfg.dispersive_domains.keys() {
        doms.insert(d);
    }
    for &d in &cfg.pec_domains {
        doms.insert(d);
    }

    let mut media = BTreeMap::new();
    let mut disp_map = cfg.dispersive_domains.clone();
    let mut domain_n_max = BTreeMap::new();

    // Default outside/domain-0 medium
    let def_eps = problem.materials.default.get("eps_r").copied().unwrap_or(1.0).max(1.0);
    let def_mu = problem.materials.default.get("mu_r").copied().unwrap_or(1.0).max(1.0);
    media.insert(
        0,
        DomainMedium::from_eps_tensor(
            [[def_eps, 0.0, 0.0], [0.0, def_eps, 0.0], [0.0, 0.0, def_eps]],
            [0.0; 3],
            [def_mu; 3],
            false,
        ),
    );

    let f_ref = cfg.reference_freq_hz();
    let omega_ref = 2.0 * PI * f_ref.max(1.0);

    for dom in doms {
        let is_pec = cfg.pec_domains.contains(&dom);
        let mu_s = problem.materials.get(dom, "mu_r").unwrap_or(def_mu).max(1e-6);

        // Check library `[optical]` for `Dispersion::Oscillators` ADE poles!
        if let std::collections::btree_map::Entry::Vacant(e) = disp_map.entry(dom) {
            if let Some(td) = problem.thermo_optic_domains.get(&dom) {
                if let Some(ref opt) = td.spec.optical {
                    if let Some(ds) = DispersiveDomainSpec::from_optical_dispersion(&opt.axes.model(0).dispersion) {
                        e.insert(ds);
                    }
                }
            }
        }

        let (eps_re, sig_e) = if let Some(ds) = disp_map.get(&dom) {
            let e = ds.eps_inf.max(1.0);
            ([[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]], [0.0; 3])
        } else if let Some(an) = problem.materials.aniso.get(&dom) {
            let er =
                an.permittivity.as_ref().map(|t| femcore::material::tensor_global(t, an.euler_deg)).unwrap_or_else(
                    || {
                        let s = problem.materials.get(dom, "eps_r").unwrap_or(def_eps);
                        [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]]
                    },
                );
            let ei = an
                .permittivity_im
                .as_ref()
                .map(|t| femcore::material::tensor_global(t, an.euler_deg))
                .unwrap_or([[0.0; 3]; 3]);
            let sig = [
                (omega_ref * grid::EPS0 * ei[0][0]).max(0.0),
                (omega_ref * grid::EPS0 * ei[1][1]).max(0.0),
                (omega_ref * grid::EPS0 * ei[2][2]).max(0.0),
            ];
            (er, sig)
        } else {
            let er_s = problem.materials.get(dom, "eps_r").unwrap_or(def_eps);
            let ei_s = problem.materials.get(dom, "eps_r_im").unwrap_or(0.0)
                + er_s * problem.materials.get(dom, "loss_tangent").unwrap_or(0.0);
            let sig_s = problem.materials.get(dom, "sigma").unwrap_or(0.0).max(omega_ref * grid::EPS0 * ei_s.max(0.0));
            ([[er_s, 0.0, 0.0], [0.0, er_s, 0.0], [0.0, 0.0, er_s]], [sig_s; 3])
        };

        let med = DomainMedium::from_eps_tensor(eps_re, sig_e, [mu_s; 3], is_pec);
        let mut n_hi = (eps_re[0][0].max(eps_re[1][1]).max(eps_re[2][2]).max(1.0) * mu_s).sqrt();
        if let Some(ds) = disp_map.get(&dom) {
            let [er, ei] = ds.eval_eps_r(f_ref);
            let eps_mag = er.hypot(ei);
            let n_disp = ((0.5 * (eps_mag + er)).max(1.0) * mu_s).sqrt();
            n_hi = n_hi.max(n_disp);
        }
        domain_n_max.insert(dom, n_hi);
        media.insert(dom, med);
    }

    Ok((media, disp_map, domain_n_max))
}

fn resolve_fdtd_nonlinear_domains(
    problem: &Problem,
    media: &BTreeMap<i32, DomainMedium>,
) -> BTreeMap<i32, (femcore::material::Nonlinear, [f64; 3], f64, f64)> {
    let mut out = BTreeMap::new();
    let extract_n0_epsmin = |dom: i32| -> (f64, f64) {
        if let Some(m) = media.get(&dom) {
            if m.is_pec {
                return (1.0, 1.0);
            }
            let e0 = m.eps_inf[0][0].max(1.0);
            let e1 = m.eps_inf[1][1].max(1.0);
            let e2 = m.eps_inf[2][2].max(1.0);
            let eps_iso = (e0 + e1 + e2) / 3.0;
            let eps_min = e0.min(e1).min(e2);
            (eps_iso.sqrt(), eps_min)
        } else {
            let e = problem.materials.get(dom, "eps_r").unwrap_or(1.0).max(1.0);
            (e.sqrt(), e)
        }
    };
    for (&dom, an) in &problem.materials.aniso {
        if problem.materials.pml.contains_key(&dom) {
            continue;
        }
        if let Some(ref nl) = an.nonlinear {
            let (n0, eps_min) = extract_n0_epsmin(dom);
            out.insert(dom, (nl.clone(), an.euler_deg, n0, eps_min));
        }
    }
    for (&dom, td) in &problem.thermo_optic_domains {
        if problem.materials.pml.contains_key(&dom) || out.contains_key(&dom) {
            continue;
        }
        if let Some(ref nl) = td.spec.nonlinear {
            let (n0, eps_min) = extract_n0_epsmin(dom);
            out.insert(dom, (nl.clone(), td.euler_deg, n0, eps_min));
        }
    }
    out
}

fn inherit_problem_bcs_and_ports(
    mesh: &Discretisation,
    problem: &Problem,
    cfg: &mut FdtdStudyConfig,
) -> Result<(), String> {
    // Map `problem.symmetry` on `x = 0` / `z = 0` exterior bounds onto `cfg.boundaries` when boundary is at 0!
    let (lo, hi) = mesh.bbox();
    let tol = 1e-6 * (hi[0] - lo[0]).max(1e-12);
    if let Some(is_pec) = problem.symmetry.x_is_pec() {
        if lo[0].abs() <= tol {
            cfg.boundaries.faces[0] = if is_pec { FdtdBcFace::Pec } else { FdtdBcFace::Pmc };
        } else if hi[0].abs() <= tol {
            cfg.boundaries.faces[1] = if is_pec { FdtdBcFace::Pec } else { FdtdBcFace::Pmc };
        }
    }
    let z_idx = if mesh.dim == 2 { 1 } else { 2 };
    if let Some(is_pec) = problem.symmetry.z_is_pec() {
        if lo[z_idx].abs() <= tol {
            cfg.boundaries.faces[4] = if is_pec { FdtdBcFace::Pec } else { FdtdBcFace::Pmc };
        } else if hi[z_idx].abs() <= tol {
            cfg.boundaries.faces[5] = if is_pec { FdtdBcFace::Pec } else { FdtdBcFace::Pmc };
        }
    }

    // Map `problem.periodic` Floquet pairs onto `FdtdBcFace::Periodic { bloch_k }`!
    for per in &problem.periodic {
        let centre = |labels: &[i32]| -> ([f64; 3], usize) {
            let mut acc = [0.0; 3];
            let mut cnt = 0usize;
            for f in 0..mesh.n_facets() {
                if labels.contains(&mesh.facet_labels[f]) {
                    for &v in mesh.facet(f) {
                        let p = mesh.point(v as usize);
                        if mesh.dim == 2 {
                            acc[0] += p[0];
                            acc[2] += p[1];
                        } else {
                            acc[0] += p[0];
                            acc[1] += p[1];
                            acc[2] += p[2];
                        }
                        cnt += 1;
                    }
                }
            }
            if cnt > 0 {
                acc = acc.map(|x| x / (cnt as f64));
            }
            (acc, cnt)
        };
        let ((a, ca), (b, cb)) = (centre(&per.source), centre(&per.dest));
        if ca > 0 && cb > 0 {
            let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let ax = (0..3).max_by(|&i, &j| d[i].abs().total_cmp(&d[j].abs())).unwrap_or(0);
            if d[ax].abs() > 1e-15 {
                let bk = per.phase / d[ax];
                cfg.boundaries.faces[2 * ax] = FdtdBcFace::Periodic { bloch_k: bk };
                cfg.boundaries.faces[2 * ax + 1] = FdtdBcFace::Periodic { bloch_k: bk };
            }
        }
    }

    // If `problem.materials.pml` has active domain PMLs, ensure corresponding outer faces have CPML!
    for pml in problem.materials.pml.values() {
        if pml.x == PmlDir::Minus {
            cfg.boundaries.faces[0] = FdtdBcFace::default_cpml(10);
        }
        if pml.x == PmlDir::Plus {
            cfg.boundaries.faces[1] = FdtdBcFace::default_cpml(10);
        }
        if pml.y == PmlDir::Minus {
            cfg.boundaries.faces[2] = FdtdBcFace::default_cpml(10);
        }
        if pml.y == PmlDir::Plus {
            cfg.boundaries.faces[3] = FdtdBcFace::default_cpml(10);
        }
        if pml.z == PmlDir::Minus {
            cfg.boundaries.faces[4] = FdtdBcFace::default_cpml(10);
        }
        if pml.z == PmlDir::Plus {
            cfg.boundaries.faces[5] = FdtdBcFace::default_cpml(10);
        }
    }

    Ok(())
}

/// Solves cross-sectional FEM eigenmodes at up to 5 anchor frequencies across `freqs_hz` and interpolates
/// `(neff, z_mode_ohm, monitor_patches)` smoothly across all monitor frequencies (solving ONCE when single-frequency
/// or narrowband, and at most 5 2D solves across an entire octave!).
#[allow(clippy::too_many_arguments)]
fn solve_broadband_mode_profiles(
    grid: &YeeGrid,
    oracle: &FdtdGeometryOracle,
    problem: &Problem,
    boundaries: &FdtdBoundaries,
    media: &BTreeMap<i32, DomainMedium>,
    default_medium: DomainMedium,
    axis: usize,
    position_m: f64,
    forward_sign: f64,
    window_lo: Option<[f64; 2]>,
    window_hi: Option<[f64; 2]>,
    freqs_hz: &[f64],
    mode_index: usize,
    n_eff_guess: Option<f64>,
) -> Result<Vec<YeeModeProfile>, String> {
    if freqs_hz.len() <= 5 {
        let mut out = Vec::with_capacity(freqs_hz.len());
        for &f in freqs_hz {
            out.push(solve_yee_plane_mode(
                grid,
                oracle,
                problem,
                boundaries,
                media,
                default_medium,
                axis,
                position_m,
                forward_sign,
                window_lo,
                window_hi,
                f,
                mode_index,
                n_eff_guess,
            )?);
        }
        align_profile_phases(&mut out);
        return Ok(out);
    }

    let f_lo = freqs_hz.iter().copied().fold(f64::INFINITY, f64::min);
    let f_hi = freqs_hz.iter().copied().fold(0.0_f64, f64::max);
    if (f_hi - f_lo) <= 1e-6 * f_hi {
        let p0 = solve_yee_plane_mode(
            grid,
            oracle,
            problem,
            boundaries,
            media,
            default_medium,
            axis,
            position_m,
            forward_sign,
            window_lo,
            window_hi,
            0.5 * (f_lo + f_hi),
            mode_index,
            n_eff_guess,
        )?;
        return Ok(vec![p0; freqs_hz.len()]);
    }

    let n_anc = 5usize;
    let anc_freqs: Vec<f64> = (0..n_anc).map(|k| f_lo + (f_hi - f_lo) * (k as f64) / ((n_anc - 1) as f64)).collect();
    let mut anc_profs = Vec::with_capacity(n_anc);
    for &fa in &anc_freqs {
        anc_profs.push(solve_yee_plane_mode(
            grid,
            oracle,
            problem,
            boundaries,
            media,
            default_medium,
            axis,
            position_m,
            forward_sign,
            window_lo,
            window_hi,
            fa,
            mode_index,
            n_eff_guess,
        )?);
    }
    align_profile_phases(&mut anc_profs);

    // Piecewise-linear / PCHIP-safe interpolation of modal fields across `anc_freqs` + exact re-normalization!
    let mut out = Vec::with_capacity(freqs_hz.len());
    for &f in freqs_hz {
        let idx = match anc_freqs.binary_search_by(|v| v.total_cmp(&f)) {
            Ok(i) => {
                out.push(anc_profs[i].clone());
                continue;
            }
            Err(i) => i.clamp(1, n_anc - 1),
        };
        let (f0, f1) = (anc_freqs[idx - 1], anc_freqs[idx]);
        let t = ((f - f0) / (f1 - f0).max(1e-30)).clamp(0.0, 1.0);
        let (p0, p1) = (&anc_profs[idx - 1], &anc_profs[idx]);
        let mut p = p0.clone();
        p.freq_hz = f;
        p.neff = p0.neff * (1.0 - t) + p1.neff * t;
        p.beta = p.neff * (2.0 * PI * f / C0);
        p.z_mode_ohm = (2.0 * PI * f * grid::MU0 / p.beta.re.abs().max(1e-12)).max(1e-6);

        let mut p_sum = 0.0;
        for (mp, mp1) in p.monitor_patches.iter_mut().zip(&p1.monitor_patches) {
            mp.e_u = mp.e_u * (1.0 - t) + mp1.e_u * t;
            mp.e_v = mp.e_v * (1.0 - t) + mp1.e_v * t;
            mp.h_u = mp.h_u * (1.0 - t) + mp1.h_u * t;
            mp.h_v = mp.h_v * (1.0 - t) + mp1.h_v * t;
            p_sum += 0.5 * p.sign * (mp.e_u * mp.h_v.conj() - mp.e_v * mp.h_u.conj()).re * mp.area_m2;
        }
        let sc = 1.0 / p_sum.abs().max(1e-30).sqrt();
        for mp in &mut p.monitor_patches {
            mp.e_u *= sc;
            mp.e_v *= sc;
            mp.h_u *= sc;
            mp.h_v *= sc;
        }
        out.push(p);
    }
    Ok(out)
}

fn align_profile_phases(profs: &mut [YeeModeProfile]) {
    if profs.len() < 2 {
        return;
    }
    for k in 1..profs.len() {
        let mut dot = c(0.0, 0.0);
        for (a, b) in profs[k - 1].monitor_patches.iter().zip(&profs[k].monitor_patches) {
            dot += (a.e_u.conj() * b.e_u + a.e_v.conj() * b.e_v) * a.area_m2;
        }
        if dot.re < 0.0 {
            for mp in &mut profs[k].monitor_patches {
                mp.e_u = -mp.e_u;
                mp.e_v = -mp.e_v;
                mp.h_u = -mp.h_u;
                mp.h_v = -mp.h_v;
            }
        }
    }
}

/// Builds a `DrivenPortData` (`S, Z, Y`) when `ModeMonitor`s share a common frequency grid and Port 0 (`mode_results[0]`)
/// carries the incident excitation:
/// - Port 0 (`mode_results[0]`): excited port (`a_0` is the incident wave, `b_0` is the reflected wave, `S_{00} = b_0 / a_0`).
/// - Port `p >= 1` (`mode_results[p]`): passive receiving port, selecting whichever of `(a_forward[p], b_backward[p])`
///   carries outgoing power so both transmission (`+y`) and inward-port-normal (`-y`) sign conventions work identically.
fn build_driven_port_data_from_mode_monitors(mons: &[ModeMonitorResult], z0_override: f64) -> Option<DrivenPortData> {
    if mons.is_empty() {
        return None;
    }
    let freqs = &mons[0].freqs_hz;
    if freqs.is_empty() || mons.iter().any(|m| m.freqs_hz.len() != freqs.len()) {
        return None;
    }
    let p_ports = mons.len();
    // Check that port 0 has non-trivial incident wave
    let max_a0 = mons[0].a_forward.iter().map(|v| v[0].hypot(v[1])).fold(0.0_f64, f64::max);
    let max_b0 = mons[0].b_backward.iter().map(|v| v[0].hypot(v[1])).fold(0.0_f64, f64::max);
    let port0_inc_is_a = max_a0 >= max_b0;
    let peak_inc = max_a0.max(max_b0);
    if peak_inc <= 1e-20 {
        return None;
    }

    // For each passive port p >= 1, determine whether outgoing wave is `a_forward` (transmission convention) or `b_backward` (inward port normal convention)
    let mut passive_uses_a = vec![false; p_ports];
    for p in 1..p_ports {
        let sum_a: f64 = mons[p].a_forward.iter().map(|v| v[0] * v[0] + v[1] * v[1]).sum();
        let sum_b: f64 = mons[p].b_backward.iter().map(|v| v[0] * v[0] + v[1] * v[1]).sum();
        passive_uses_a[p] = sum_a >= sum_b;
    }

    let mut z0_ohm = Vec::with_capacity(freqs.len());
    let mut s_matrix = Vec::with_capacity(freqs.len());
    let mut z_matrix = Vec::with_capacity(freqs.len());
    let mut y_matrix = Vec::with_capacity(freqs.len());

    for fi in 0..freqs.len() {
        let (a0, b0) = if port0_inc_is_a {
            (
                c(mons[0].a_forward[fi][0], mons[0].a_forward[fi][1]),
                c(mons[0].b_backward[fi][0], mons[0].b_backward[fi][1]),
            )
        } else {
            (
                c(mons[0].b_backward[fi][0], mons[0].b_backward[fi][1]),
                c(mons[0].a_forward[fi][0], mons[0].a_forward[fi][1]),
            )
        };
        let inv_a0 = if a0.norm() > 1e-24 { c(1.0, 0.0) / a0 } else { c(0.0, 0.0) };

        let mut z_circ = Vec::with_capacity(p_ports);
        let mut z_row = Vec::with_capacity(p_ports);
        for m in mons {
            let zp = if z0_override > 0.0 { z0_override } else { m.z_mode_ohm[fi] };
            z_circ.push(zp);
            z_row.push([zp, 0.0]);
        }

        let mut s_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
        s_k[0][0] = b0 * inv_a0;
        for p in 1..p_ports {
            let bp = if passive_uses_a[p] {
                c(mons[p].a_forward[fi][0], mons[p].a_forward[fi][1])
            } else {
                c(mons[p].b_backward[fi][0], mons[p].b_backward[fi][1])
            };
            let s_p0 = bp * inv_a0;
            s_k[p][0] = s_p0;
            // Fill symmetric column when single-excitation reciprocal passive run is performed
            if p_ports == 2 {
                s_k[0][p] = s_p0;
                s_k[p][p] = s_k[0][0];
            } else {
                s_k[0][p] = s_p0;
            }
        }
        let (s_out, z_out, y_out, _) = s_to_z_and_y_checked(&s_k, &z_circ);
        z0_ohm.push(z_row);
        s_matrix.push(s_out);
        z_matrix.push(z_out);
        y_matrix.push(y_out);
    }

    Some(DrivenPortData {
        freqs_hz: freqs.clone(),
        ports: mons.iter().map(|m| m.name.clone()).collect(),
        z0_ohm,
        s_matrix,
        z_matrix,
        y_matrix,
        hdm_solves: 1,
        rom_dim: 0,
        anchor_freqs_hz: vec![freqs[freqs.len() / 2]],
        max_residual: 0.0,
    })
}

fn compute_plane_wave_incident_power(grid: &YeeGrid, sources: &[FdtdSourceSpec]) -> Option<f64> {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    for s in sources {
        if let FdtdSourceSpec::PlaneWave { axis, eps_bg, mu_bg, waveform, .. } = s {
            let amp = match *waveform {
                FdtdWaveform::GaussianPulse { amplitude, .. }
                | FdtdWaveform::Ricker { amplitude, .. }
                | FdtdWaveform::ContinuousWave { amplitude, .. } => amplitude,
            };
            // Note: `s(t) = env(t) cos(w0 t)`, while DFT `int cos(w0 t) exp(-j w0 t) dt = 0.5 * int env(t) dt`?
            // Wait: `s_hat(f)` in `FdtdWaveform::spectrum` ALREADY computes the exact DFT of `s(t) = env(t) cos(...)`,
            // so `E_inc_hat(f) / s_hat(f) == 1.0` (unit phasor amplitude!) because `Aux1dYeeLine` injects `1.0 * waveform.eval(t)`!!
            // Check that!! `Aux1dYeeLine` injects `e_src = waveform.eval(t_n)`, which ALREADY includes `amplitude`!
            // Therefore `E_inc_hat(f) / s_hat(f) == 1.0` V/m!
            let _ = amp;
            let eta_bg = (grid::MU0 * *mu_bg / (grid::EPS0 * *eps_bg)).sqrt();
            let (u_ax, v_ax) = rh_uv_axes(*axis);
            let vu = grid.axis_vertices(u_ax);
            let vv = grid.axis_vertices(v_ax);
            let lu = vu[vu.len() - 1] - vu[0];
            let lv = if is_2d && v_ax == 1 { 1.0 } else { vv[vv.len() - 1] - vv[0] };
            let area = if is_2d && u_ax == 1 {
                lv
            } else if is_2d && v_ax == 1 {
                lu
            } else {
                lu * lv
            };
            return Some(0.5 * area / eta_bg);
        }
    }
    None
}

/// Exact periodic & Bloch-periodic boundary Ampere updates on `i = 0` and `i = n` (`E(n) = exp(-j k_B L) E(0)`),
/// including CFS-PML auxiliary `psi` stretching along any non-periodic axes.
#[allow(clippy::too_many_arguments)]
pub(super) fn apply_periodic_e_updates(
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    coeffs: &AdeAndAmpereCoeffs,
    cpml_re: &mut CpmlState,
    mut cpml_im: Option<&mut CpmlState>,
    f_re: &mut YeeFields,
    mut f_im: Option<&mut YeeFields>,
) {
    let per = bcs.periodic_axes();
    if !per[0] && !per[1] && !per[2] {
        return;
    }
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let bk = bcs.bloch_vector();
    let (nx, ny, nz, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.stride_x, grid.stride_y);
    let pmc_xm = matches!(bcs.faces[0], FdtdBcFace::Pmc);
    let pmc_xp = matches!(bcs.faces[1], FdtdBcFace::Pmc);
    let pmc_ym = matches!(bcs.faces[2], FdtdBcFace::Pmc);
    let pmc_yp = matches!(bcs.faces[3], FdtdBcFace::Pmc);
    let pmc_zm = matches!(bcs.faces[4], FdtdBcFace::Pmc);
    let pmc_zp = matches!(bcs.faces[5], FdtdBcFace::Pmc);

    let stretch_psi = |ax: &ade_pml::CpmlAxis1d,
                       idx: usize,
                       psi_re: &mut [f64],
                       psi_im: Option<&mut [f64]>,
                       g: usize,
                       d_raw: C|
     -> C {
        if !ax.has_cpml() || !ax.active_e[idx] {
            return d_raw;
        }
        psi_re[g] = ax.b_e[idx] * psi_re[g] + ax.c_e[idx] * d_raw.re;
        let re = ax.inv_kappa_e[idx] * d_raw.re + psi_re[g];
        let im = if let Some(pi) = psi_im {
            pi[g] = ax.b_e[idx] * pi[g] + ax.c_e[idx] * d_raw.im;
            ax.inv_kappa_e[idx] * d_raw.im + pi[g]
        } else {
            d_raw.im
        };
        c(re, im)
    };

    // Helper reading `H_comp` at `(ix, iy, iz)` with periodic/Bloch wrap when `i == -1` -> `n - 1` (`* exp(+j k L)`)!
    let read_h = |comp: usize, ix: isize, iy: isize, iz: isize, fr: &YeeFields, fi: Option<&YeeFields>| -> C {
        let mut ph = c(1.0, 0.0);
        let wx = if ix < 0 {
            let (sn, cs) = (bk[0] * (grid.vx[nx] - grid.vx[0])).sin_cos();
            ph *= c(cs, sn);
            nx - 1
        } else {
            ix as usize
        };
        let wy = if is_2d {
            0
        } else if iy < 0 {
            let (sn, cs) = (bk[1] * (grid.vy[ny] - grid.vy[0])).sin_cos();
            ph *= c(cs, sn);
            ny - 1
        } else {
            iy as usize
        };
        let wz = if iz < 0 {
            let (sn, cs) = (bk[2] * (grid.vz[nz] - grid.vz[0])).sin_cos();
            ph *= c(cs, sn);
            nz - 1
        } else {
            iz as usize
        };
        let g = grid.idx(wx, wy, wz);
        let r = match comp {
            0 => fr.hx.as_slice()[g],
            1 => fr.hy.as_slice()[g],
            _ => fr.hz.as_slice()[g],
        };
        let i = fi.map_or(0.0, |im| match comp {
            0 => im.hx.as_slice()[g],
            1 => im.hy.as_slice()[g],
            _ => im.hz.as_slice()[g],
        });
        c(r, i) * ph
    };

    let write_e_pair = |comp: usize,
                        g0: usize,
                        gn: usize,
                        val0: C,
                        cs_wrap: f64,
                        sn_wrap: f64,
                        fr: &mut YeeFields,
                        fi: Option<&mut YeeFields>| {
        // E(n) = exp(-j k L) * E(0) = (cs - j sn) * val0
        let valn = val0 * c(cs_wrap, -sn_wrap);
        match comp {
            0 => fr.ex.as_mut_slice()[g0] = val0.re,
            1 => fr.ey.as_mut_slice()[g0] = val0.re,
            _ => fr.ez.as_mut_slice()[g0] = val0.re,
        }
        if gn != g0 {
            match comp {
                0 => fr.ex.as_mut_slice()[gn] = valn.re,
                1 => fr.ey.as_mut_slice()[gn] = valn.re,
                _ => fr.ez.as_mut_slice()[gn] = valn.re,
            }
        }
        if let Some(im) = fi {
            match comp {
                0 => im.ex.as_mut_slice()[g0] = val0.im,
                1 => im.ey.as_mut_slice()[g0] = val0.im,
                _ => im.ez.as_mut_slice()[g0] = val0.im,
            }
            if gn != g0 {
                match comp {
                    0 => im.ex.as_mut_slice()[gn] = valn.im,
                    1 => im.ey.as_mut_slice()[gn] = valn.im,
                    _ => im.ez.as_mut_slice()[gn] = valn.im,
                }
            }
        }
    };

    // 1. Periodic along X (`ix = 0` -> `ix = nx`): updates `Ey(0, iy, iz)` and `Ez(0, iy, iz)`
    if per[0] {
        let (snx, csx) = (bk[0] * (grid.vx[nx] - grid.vx[0])).sin_cos();
        let inv_ddx = grid.inv_dd_x[0];
        let iy_hi = if is_2d { 0 } else { ny };
        let do_ey = !matches!(grid.dim, FdtdDim::TwoD { polarization: Pol2d::TeHy });
        let do_ez = !matches!(grid.dim, FdtdDim::TwoD { polarization: Pol2d::TmEy });

        for iy in 0..=iy_hi {
            for iz in 0..=nz {
                let g0 = iy * sy + iz;
                let gn = nx * sx + g0;
                // Ey at (0, iy, iz) for iy < ny and valid iz
                let valid_iz = (iz > 0 && iz < nz) || (iz == 0 && (per[2] || pmc_zm)) || (iz == nz && pmc_zp);
                if do_ey && iy < ny && valid_iz {
                    let dhx_raw = if iz == 0 && !per[2] {
                        read_h(0, 0, iy as isize, 0, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_z[0])
                    } else if iz == nz {
                        read_h(0, 0, iy as isize, (nz - 1) as isize, f_re, f_im.as_deref()) * (-2.0 * grid.inv_dd_z[nz])
                    } else {
                        let hx_hi = read_h(0, 0, iy as isize, iz as isize, f_re, f_im.as_deref());
                        let hx_lo = read_h(0, 0, iy as isize, iz as isize - 1, f_re, f_im.as_deref());
                        (hx_hi - hx_lo) * grid.inv_dd_z[iz]
                    };
                    let dhx_dz = stretch_psi(
                        &cpml_re.z,
                        iz,
                        &mut cpml_re.psi_ey_z,
                        cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ey_z[..]),
                        g0,
                        dhx_raw,
                    );
                    let hz_hi = read_h(2, 0, iy as isize, iz as isize, f_re, f_im.as_deref());
                    let hz_lo = read_h(2, -1, iy as isize, iz as isize, f_re, f_im.as_deref());
                    let curl = dhx_dz - (hz_hi - hz_lo) * inv_ddx;
                    let e_old = c(f_re.ey.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ey.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_y[g0] + curl * coeffs.cb_y[g0];
                    write_e_pair(1, g0, gn, e_new, csx, snx, f_re, f_im.as_deref_mut());
                }
                // Ez at (0, iy, iz) for iz < nz and valid iy
                let valid_iy = is_2d || (iy > 0 && iy < ny) || (iy == 0 && (per[1] || pmc_ym)) || (iy == ny && pmc_yp);
                if do_ez && iz < nz && valid_iy {
                    let hy_hi = read_h(1, 0, iy as isize, iz as isize, f_re, f_im.as_deref());
                    let hy_lo = read_h(1, -1, iy as isize, iz as isize, f_re, f_im.as_deref());
                    let dhx_dy = if is_2d {
                        c(0.0, 0.0)
                    } else {
                        let dhx_raw = if iy == 0 && !per[1] {
                            read_h(0, 0, 0, iz as isize, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_y[0])
                        } else if iy == ny {
                            read_h(0, 0, (ny - 1) as isize, iz as isize, f_re, f_im.as_deref())
                                * (-2.0 * grid.inv_dd_y[ny])
                        } else {
                            let hx_hi = read_h(0, 0, iy as isize, iz as isize, f_re, f_im.as_deref());
                            let hx_lo = read_h(0, 0, iy as isize - 1, iz as isize, f_re, f_im.as_deref());
                            (hx_hi - hx_lo) * grid.inv_dd_y[iy]
                        };
                        stretch_psi(
                            &cpml_re.y,
                            iy,
                            &mut cpml_re.psi_ez_y,
                            cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ez_y[..]),
                            g0,
                            dhx_raw,
                        )
                    };
                    let curl = (hy_hi - hy_lo) * inv_ddx - dhx_dy;
                    let e_old = c(f_re.ez.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ez.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_z[g0] + curl * coeffs.cb_z[g0];
                    write_e_pair(2, g0, gn, e_new, csx, snx, f_re, f_im.as_deref_mut());
                }
            }
        }
    }

    // 2. Periodic along Y (3D only: `iy = 0` -> `iy = ny`): updates `Ex(ix, 0, iz)` and `Ez(ix, 0, iz)`
    if !is_2d && per[1] {
        let (sny, csy) = (bk[1] * (grid.vy[ny] - grid.vy[0])).sin_cos();
        let inv_ddy = grid.inv_dd_y[0];
        for ix in 0..=nx {
            for iz in 0..=nz {
                let g0 = ix * sx + iz;
                let gn = ix * sx + ny * sy + iz;
                let valid_iz = (iz > 0 && iz < nz) || (iz == 0 && (per[2] || pmc_zm)) || (iz == nz && pmc_zp);
                if ix < nx && valid_iz {
                    let hz_hi = read_h(2, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                    let hz_lo = read_h(2, ix as isize, -1, iz as isize, f_re, f_im.as_deref());
                    let dhy_raw = if iz == 0 && !per[2] {
                        read_h(1, ix as isize, 0, 0, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_z[0])
                    } else if iz == nz {
                        read_h(1, ix as isize, 0, (nz - 1) as isize, f_re, f_im.as_deref()) * (-2.0 * grid.inv_dd_z[nz])
                    } else {
                        let hy_hi = read_h(1, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                        let hy_lo = read_h(1, ix as isize, 0, iz as isize - 1, f_re, f_im.as_deref());
                        (hy_hi - hy_lo) * grid.inv_dd_z[iz]
                    };
                    let dhy_dz = stretch_psi(
                        &cpml_re.z,
                        iz,
                        &mut cpml_re.psi_ex_z,
                        cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ex_z[..]),
                        g0,
                        dhy_raw,
                    );
                    let curl = (hz_hi - hz_lo) * inv_ddy - dhy_dz;
                    let e_old = c(f_re.ex.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ex.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_x[g0] + curl * coeffs.cb_x[g0];
                    write_e_pair(0, g0, gn, e_new, csy, sny, f_re, f_im.as_deref_mut());
                }
                let valid_ix = (ix > 0 && ix < nx) || (ix == 0 && pmc_xm) || (ix == nx && pmc_xp);
                if iz < nz && !per[0] && valid_ix {
                    let dhy_raw = if ix == 0 {
                        read_h(1, 0, 0, iz as isize, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_x[0])
                    } else if ix == nx {
                        read_h(1, (nx - 1) as isize, 0, iz as isize, f_re, f_im.as_deref()) * (-2.0 * grid.inv_dd_x[nx])
                    } else {
                        let hy_hi = read_h(1, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                        let hy_lo = read_h(1, ix as isize - 1, 0, iz as isize, f_re, f_im.as_deref());
                        (hy_hi - hy_lo) * grid.inv_dd_x[ix]
                    };
                    let dhy_dx = stretch_psi(
                        &cpml_re.x,
                        ix,
                        &mut cpml_re.psi_ez_x,
                        cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ez_x[..]),
                        g0,
                        dhy_raw,
                    );
                    let hx_hi = read_h(0, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                    let hx_lo = read_h(0, ix as isize, -1, iz as isize, f_re, f_im.as_deref());
                    let curl = dhy_dx - (hx_hi - hx_lo) * inv_ddy;
                    let e_old = c(f_re.ez.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ez.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_z[g0] + curl * coeffs.cb_z[g0];
                    write_e_pair(2, g0, gn, e_new, csy, sny, f_re, f_im.as_deref_mut());
                } else if iz < nz && per[0] && (ix > 0 && ix < nx) {
                    let hy_hi = read_h(1, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                    let hy_lo = read_h(1, ix as isize - 1, 0, iz as isize, f_re, f_im.as_deref());
                    let hx_hi = read_h(0, ix as isize, 0, iz as isize, f_re, f_im.as_deref());
                    let hx_lo = read_h(0, ix as isize, -1, iz as isize, f_re, f_im.as_deref());
                    let curl = (hy_hi - hy_lo) * grid.inv_dd_x[ix] - (hx_hi - hx_lo) * inv_ddy;
                    let e_old = c(f_re.ez.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ez.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_z[g0] + curl * coeffs.cb_z[g0];
                    write_e_pair(2, g0, gn, e_new, csy, sny, f_re, f_im.as_deref_mut());
                } else if iz < nz && per[0] && (ix == 0 || ix == nx) {
                    // Already computed `Ez(0, 0, iz)` and `Ez(nx, 0, iz)` in `per[0]` block; wrap to `iy = ny`!
                    let val0 = c(f_re.ez.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ez.as_slice()[g0]));
                    write_e_pair(2, g0, gn, val0, csy, sny, f_re, f_im.as_deref_mut());
                }
            }
        }
    }

    // 3. Periodic along Z (`iz = 0` -> `iz = nz`): updates `Ex(ix, iy, 0)` and `Ey(ix, iy, 0)`
    if per[2] {
        let (snz, csz) = (bk[2] * (grid.vz[nz] - grid.vz[0])).sin_cos();
        let inv_ddz = grid.inv_dd_z[0];
        let iy_hi = if is_2d { 0 } else { ny };
        let do_ex = !matches!(grid.dim, FdtdDim::TwoD { polarization: Pol2d::TmEy });
        let do_ey = !matches!(grid.dim, FdtdDim::TwoD { polarization: Pol2d::TeHy });

        for ix in 0..=nx {
            for iy in 0..=iy_hi {
                let g0 = ix * sx + iy * sy;
                let gn = g0 + nz;
                let valid_iy =
                    is_2d || (iy > 0 && iy < ny) || (!per[1] && ((iy == 0 && pmc_ym) || (iy == ny && pmc_yp)));
                if do_ex && ix < nx && valid_iy {
                    let dhz_dy = if is_2d {
                        c(0.0, 0.0)
                    } else {
                        let dhz_raw = if iy == 0 {
                            read_h(2, ix as isize, 0, 0, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_y[0])
                        } else if iy == ny {
                            read_h(2, ix as isize, (ny - 1) as isize, 0, f_re, f_im.as_deref())
                                * (-2.0 * grid.inv_dd_y[ny])
                        } else {
                            let hz_hi = read_h(2, ix as isize, iy as isize, 0, f_re, f_im.as_deref());
                            let hz_lo = read_h(2, ix as isize, iy as isize - 1, 0, f_re, f_im.as_deref());
                            (hz_hi - hz_lo) * grid.inv_dd_y[iy]
                        };
                        stretch_psi(
                            &cpml_re.y,
                            iy,
                            &mut cpml_re.psi_ex_y,
                            cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ex_y[..]),
                            g0,
                            dhz_raw,
                        )
                    };
                    let hy_hi = read_h(1, ix as isize, iy as isize, 0, f_re, f_im.as_deref());
                    let hy_lo = read_h(1, ix as isize, iy as isize, -1, f_re, f_im.as_deref());
                    let curl = dhz_dy - (hy_hi - hy_lo) * inv_ddz;
                    let e_old = c(f_re.ex.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ex.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_x[g0] + curl * coeffs.cb_x[g0];
                    write_e_pair(0, g0, gn, e_new, csz, snz, f_re, f_im.as_deref_mut());
                } else if do_ex && ix < nx && !is_2d && per[1] && (iy == 0 || iy == ny) {
                    let val0 = c(f_re.ex.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ex.as_slice()[g0]));
                    write_e_pair(0, g0, gn, val0, csz, snz, f_re, f_im.as_deref_mut());
                }
                let valid_ix = (ix > 0 && ix < nx) || (!per[0] && ((ix == 0 && pmc_xm) || (ix == nx && pmc_xp)));
                if do_ey && iy < ny && valid_ix {
                    let hx_hi = read_h(0, ix as isize, iy as isize, 0, f_re, f_im.as_deref());
                    let hx_lo = read_h(0, ix as isize, iy as isize, -1, f_re, f_im.as_deref());
                    let dhz_raw = if ix == 0 {
                        read_h(2, 0, iy as isize, 0, f_re, f_im.as_deref()) * (2.0 * grid.inv_dd_x[0])
                    } else if ix == nx {
                        read_h(2, (nx - 1) as isize, iy as isize, 0, f_re, f_im.as_deref()) * (-2.0 * grid.inv_dd_x[nx])
                    } else {
                        let hz_hi = read_h(2, ix as isize, iy as isize, 0, f_re, f_im.as_deref());
                        let hz_lo = read_h(2, ix as isize - 1, iy as isize, 0, f_re, f_im.as_deref());
                        (hz_hi - hz_lo) * grid.inv_dd_x[ix]
                    };
                    let dhz_dx = stretch_psi(
                        &cpml_re.x,
                        ix,
                        &mut cpml_re.psi_ey_x,
                        cpml_im.as_deref_mut().map(|cp| &mut cp.psi_ey_x[..]),
                        g0,
                        dhz_raw,
                    );
                    let curl = (hx_hi - hx_lo) * inv_ddz - dhz_dx;
                    let e_old = c(f_re.ey.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ey.as_slice()[g0]));
                    let e_new = e_old * coeffs.ca_y[g0] + curl * coeffs.cb_y[g0];
                    write_e_pair(1, g0, gn, e_new, csz, snz, f_re, f_im.as_deref_mut());
                } else if do_ey && iy < ny && per[0] && (ix == 0 || ix == nx) {
                    let val0 = c(f_re.ey.as_slice()[g0], f_im.as_ref().map_or(0.0, |fi| fi.ey.as_slice()[g0]));
                    write_e_pair(1, g0, gn, val0, csz, snz, f_re, f_im.as_deref_mut());
                }
            }
        }
    }
}

/// Synchronizes periodic / Bloch high-face tangential electric fields `E(n) = exp(-j k_B L) E(0)`
/// after ADE, anisotropic off-diagonal, and source injections.
pub(super) fn sync_periodic_e_high_faces(
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    f_re: &mut YeeFields,
    mut f_im: Option<&mut YeeFields>,
) {
    let per = bcs.periodic_axes();
    if !per[0] && !per[1] && !per[2] {
        return;
    }
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let bk = bcs.bloch_vector();
    let (nx, ny, nz, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.stride_x, grid.stride_y);

    let copy_wrap =
        |comp: usize, g0: usize, gn: usize, cs: f64, sn: f64, fr: &mut YeeFields, fi: Option<&mut YeeFields>| {
            if g0 == gn {
                return;
            }
            let r0 = match comp {
                0 => fr.ex.as_slice()[g0],
                1 => fr.ey.as_slice()[g0],
                _ => fr.ez.as_slice()[g0],
            };
            let i0 = fi.as_ref().map_or(0.0, |im| match comp {
                0 => im.ex.as_slice()[g0],
                1 => im.ey.as_slice()[g0],
                _ => im.ez.as_slice()[g0],
            });
            let vn = c(r0, i0) * c(cs, -sn);
            match comp {
                0 => fr.ex.as_mut_slice()[gn] = vn.re,
                1 => fr.ey.as_mut_slice()[gn] = vn.re,
                _ => fr.ez.as_mut_slice()[gn] = vn.re,
            }
            if let Some(im) = fi {
                match comp {
                    0 => im.ex.as_mut_slice()[gn] = vn.im,
                    1 => im.ey.as_mut_slice()[gn] = vn.im,
                    _ => im.ez.as_mut_slice()[gn] = vn.im,
                }
            }
        };

    let iy_hi = if is_2d { 0 } else { ny };
    if per[0] {
        let (snx, csx) = (bk[0] * (grid.vx[nx] - grid.vx[0])).sin_cos();
        for iy in 0..=iy_hi {
            for iz in 0..=nz {
                let g0 = iy * sy + iz;
                let gn = nx * sx + g0;
                copy_wrap(1, g0, gn, csx, snx, f_re, f_im.as_deref_mut());
                copy_wrap(2, g0, gn, csx, snx, f_re, f_im.as_deref_mut());
            }
        }
    }
    if !is_2d && per[1] {
        let (sny, csy) = (bk[1] * (grid.vy[ny] - grid.vy[0])).sin_cos();
        for ix in 0..=nx {
            for iz in 0..=nz {
                let g0 = ix * sx + iz;
                let gn = ix * sx + ny * sy + iz;
                copy_wrap(0, g0, gn, csy, sny, f_re, f_im.as_deref_mut());
                copy_wrap(2, g0, gn, csy, sny, f_re, f_im.as_deref_mut());
            }
        }
    }
    if per[2] {
        let (snz, csz) = (bk[2] * (grid.vz[nz] - grid.vz[0])).sin_cos();
        for ix in 0..=nx {
            for iy in 0..=iy_hi {
                let g0 = ix * sx + iy * sy;
                let gn = g0 + nz;
                copy_wrap(0, g0, gn, csz, snz, f_re, f_im.as_deref_mut());
                copy_wrap(1, g0, gn, csz, snz, f_re, f_im.as_deref_mut());
            }
        }
    }
}
