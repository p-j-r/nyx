// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The typed problem of `ems`: what the kernels read. It is decoded from the `ProblemSpec` envelope by `decode` (region names and selectors
//! resolved to mesh numbers, parameters with their defaults), or built directly by code that already has mesh numbers (tests, the legacy
//! window bridge of the engine). The envelope carries scalar permittivity and permeability; a crystal's permittivity tensor reaches a kernel
//! through `materials.aniso` until the material files exist (phase 6).
//!
//! The study names the kernel: `Modes` (guided modes of a 2D cross-section, real or complex), `AxisymmetricModes` (a body of revolution on
//! its (r, z) half plane) or `Cavity` (3D eigenfrequencies). It is not guessed from the dimension of the mesh (ROUGH_SPOTS R10): a study
//! on a mesh of the other dimension is refused by name.

use crate::ibc::resolve_jj_l_and_ej;
use crate::Ems;
use femcore::material::{bcs_surface_impedance, FilmModel};
use femcore::props::{Materials, Periodic, Pml, PmlCoord, PmlDir, PmlProfile, TensorKind};
use nyx_physics::decode::{
    domain_names_by_label, materials_by_label, resolve_conditions, tls_interfaces_from_spec, Params,
};
use nyx_physics::{
    Ctx, Discretisation, Error, MaterialAssign, MixedDiscretisation, Physics, ProblemSpec, TlsInterfaces,
};
use std::collections::BTreeMap;

/// Allowed values for the midplane `symmetry` parameter on `modes`, `cavity`, and `bloch`.
pub const SYMMETRY_CHOICES: [&str; 9] =
    ["none", "pec_x", "pmc_x", "pec_z", "pmc_z", "pec_x_pmc_z", "pmc_x_pec_z", "pmc_x_pmc_z", "pec_x_pec_z"];

/// Midplane PEC/PMC symmetry condition on exterior facets lying on `x = 0` and/or `z = 0`.
/// On a quarter mesh `x >= 0, z >= 0`, fundamental `TE0` is `PecXPmcZ` (`pec_x_pmc_z`) and `TM0` is `PmcXPecZ` (`pmc_x_pec_z`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Symmetry {
    #[default]
    None,
    PecX,
    PmcX,
    PecZ,
    PmcZ,
    PecXPmcZ,
    PmcXPecZ,
    PmcXPmcZ,
    PecXPecZ,
}

impl Symmetry {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "none" | "" => Ok(Self::None),
            "pec_x" => Ok(Self::PecX),
            "pmc_x" => Ok(Self::PmcX),
            "pec_z" => Ok(Self::PecZ),
            "pmc_z" => Ok(Self::PmcZ),
            "pec_x_pmc_z" => Ok(Self::PecXPmcZ),
            "pmc_x_pec_z" => Ok(Self::PmcXPecZ),
            "pmc_x_pmc_z" => Ok(Self::PmcXPmcZ),
            "pec_x_pec_z" => Ok(Self::PecXPecZ),
            other => Err(format!("unknown symmetry '{other}' (expected one of {SYMMETRY_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::PecX => "pec_x",
            Self::PmcX => "pmc_x",
            Self::PecZ => "pec_z",
            Self::PmcZ => "pmc_z",
            Self::PecXPmcZ => "pec_x_pmc_z",
            Self::PmcXPecZ => "pmc_x_pec_z",
            Self::PmcXPmcZ => "pmc_x_pmc_z",
            Self::PecXPecZ => "pec_x_pec_z",
        }
    }

    /// `Some(true)` for PEC on `x = 0`, `Some(false)` for PMC on `x = 0`, `None` otherwise.
    pub fn x_is_pec(self) -> Option<bool> {
        match self {
            Self::PecX | Self::PecXPmcZ | Self::PecXPecZ => Some(true),
            Self::PmcX | Self::PmcXPecZ | Self::PmcXPmcZ => Some(false),
            Self::None | Self::PecZ | Self::PmcZ => None,
        }
    }

    /// `Some(true)` for PEC on `z = 0`, `Some(false)` for PMC on `z = 0`, `None` otherwise.
    pub fn z_is_pec(self) -> Option<bool> {
        match self {
            Self::PecZ | Self::PmcXPecZ | Self::PecXPecZ => Some(true),
            Self::PmcZ | Self::PecXPmcZ | Self::PmcXPmcZ => Some(false),
            Self::None | Self::PecX | Self::PmcX => None,
        }
    }
}

/// Auto-detects exterior facets lying on `x = 0` (when `|xmin| <= tol` or `|xmax| <= tol`) and/or
/// `z = 0` (when `|zmin| <= tol` or `|zmax| <= tol`), returning `(sym_pec_facets, sym_pmc_facets)`.
pub fn symmetry_facet_flags(mesh: &Discretisation, symmetry: Symmetry) -> Result<(Vec<bool>, Vec<bool>), String> {
    let nf = mesh.n_facets();
    let mut sym_pec = vec![false; nf];
    let mut sym_pmc = vec![false; nf];
    if symmetry == Symmetry::None {
        return Ok((sym_pec, sym_pmc));
    }
    let (lo, hi) = mesh.bbox();
    let extent = (0..mesh.dim).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let x_idx = 0usize;
    let z_idx = if mesh.dim == 2 { 1usize } else { 2usize };
    let x_rule = symmetry.x_is_pec();
    let z_rule = symmetry.z_is_pec();
    if x_rule.is_some() && !(lo[x_idx].abs() <= tol || hi[x_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on x = 0 (mesh x spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[x_idx],
            hi[x_idx]
        ));
    }
    if z_rule.is_some() && !(lo[z_idx].abs() <= tol || hi[z_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on z = 0 (mesh z spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[z_idx],
            hi[z_idx]
        ));
    }
    let exterior = mesh.exterior_facets();
    for f in 0..nf {
        if !exterior[f] {
            continue;
        }
        let fv = mesh.facet(f);
        if let Some(is_pec) = x_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[x_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
        if let Some(is_pec) = z_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[z_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
    }
    Ok((sym_pec, sym_pmc))
}

/// Auto-detects exterior facets lying on `x = 0` and/or `z = 0` on a native [`MixedDiscretisation`].
pub fn symmetry_facet_flags_mixed(
    mesh: &MixedDiscretisation,
    symmetry: Symmetry,
) -> Result<(Vec<bool>, Vec<bool>), String> {
    let nf = mesh.n_facets();
    let mut sym_pec = vec![false; nf];
    let mut sym_pmc = vec![false; nf];
    if symmetry == Symmetry::None {
        return Ok((sym_pec, sym_pmc));
    }
    let (lo, hi) = mesh.bbox();
    let extent = (0..mesh.dim).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let x_idx = 0usize;
    let z_idx = if mesh.dim == 2 { 1usize } else { 2usize };
    let x_rule = symmetry.x_is_pec();
    let z_rule = symmetry.z_is_pec();
    if x_rule.is_some() && !(lo[x_idx].abs() <= tol || hi[x_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on x = 0 (mesh x spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[x_idx],
            hi[x_idx]
        ));
    }
    if z_rule.is_some() && !(lo[z_idx].abs() <= tol || hi[z_idx].abs() <= tol) {
        return Err(format!(
            "symmetry '{}' requires an exterior boundary on z = 0 (mesh z spans [{:.3e}, {:.3e}])",
            symmetry.as_str(),
            lo[z_idx],
            hi[z_idx]
        ));
    }
    for f in 0..nf {
        if !mesh.exterior_facets[f] {
            continue;
        }
        let fv = mesh.facet(f);
        if let Some(is_pec) = x_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[x_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
        if let Some(is_pec) = z_rule {
            if fv.iter().all(|&v| mesh.point(v as usize)[z_idx].abs() <= tol) {
                if is_pec {
                    sym_pec[f] = true;
                } else {
                    sym_pmc[f] = true;
                }
            }
        }
    }
    Ok((sym_pec, sym_pmc))
}

/// Collects PEC boundary edge indices `be` and vertex indices `bn` for a 2D mesh, honoring explicit `pec`
/// facets and midplane `symmetry` (`pec_*` / `pmc_*`), while excluding `surface_impedance` and `pmc` facets.
pub(crate) fn pec_edges_and_nodes_2d_with_symmetry(
    mesh: &Discretisation,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
    symmetry: Symmetry,
) -> Result<(Vec<usize>, Vec<usize>), String> {
    let exterior = mesh.exterior_facets();
    let (sym_pec, sym_pmc) = symmetry_facet_flags(mesh, symmetry)?;
    let ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let is_pec_facet = |f: usize| -> bool {
        let l = mesh.facet_labels[f];
        pec_labels.contains(&l)
            || sym_pec[f]
            || (exterior[f] && !ibc_labels.contains(&l) && !pmc_labels.contains(&l) && !sym_pmc[f])
    };
    let mut be: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
        .collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| is_pec_facet(f))
        .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
        .collect();
    bn.sort_unstable();
    bn.dedup();
    Ok((be, bn))
}

#[derive(Clone, Debug, PartialEq)]
pub enum Study {
    /// Guided modes at `freq` (Hz) near the effective index `n_target` (imaginary if `evanescent`). Absorbing layers, loss and a crystal
    /// tilted against the guide take the complex solver.
    Modes { freq: f64, n_target: f64, evanescent: bool, n_modes: usize },
    /// Modes of a body of revolution at azimuthal order `m` (`exp(-j m theta)`, theta along wafer y) nearest `f_search` (Hz); the mesh is
    /// the (r, z) half plane with x = r. A permittivity per domain (isotropic), absorbing layers by `materials.pml`.
    AxisymmetricModes { f_search: f64, azimuthal_order: i64, n_modes: usize },
    /// Eigenfrequencies of a 3D cavity or Floquet cell nearest `f_search` (Hz; a value of zero or less takes the fundamental of a box of
    /// the size of the mesh).
    Cavity { n_modes: usize, f_search: f64 },
    /// Fixed-frequency 3D Floquet-Bloch modes at `freq` (Hz) near the Bloch index `n_target` in Brillouin zone `brillouin_zone`.
    Bloch { freq: f64, n_target: f64, n_modes: usize, brillouin_zone: i64 },
    /// Driven frequency-domain multi-port S/Z/Y sweep over `[f_min, f_max]` (`n_points`) or an explicit ascending
    /// frequency grid `freqs_hz`, solved with Hermite/Galerkin PROM (`prom_tol = Some(tol)`) or direct point-by-point
    /// HDM (`prom_tol = None`).
    Driven {
        f_min: f64,
        f_max: f64,
        n_points: usize,
        freqs_hz: Option<Vec<f64>>,
        prom_tol: Option<f64>,
        max_hdm_solves: usize,
        strict_prom_tol: bool,
        seed_freqs_hz: Vec<f64>,
        solver: String,
        blr_tol: f64,
    },
    /// 3D eigenmode expansion (`"eme"`, ADR 140) along `+y` over a segmented / periodic `EmeGridSpec` partition
    /// with BLAKE3 cross-section deduplication, cross-power interface scattering matrices, and Redheffer star products.
    Eme {
        freq: f64,
        n_target: f64,
        evanescent: bool,
        n_modes: usize,
        grid: crate::eme::EmeGridSpec,
        y_range: Option<[f64; 2]>,
        port_modes: usize,
        port_z0_ohm: f64,
        reconstruct_fields: bool,
    },
    /// 3D frequency-domain nonlinear Harmonic Balance & Coupled-Wave FEM (`"nonlinear_driven"`, ADR 143):
    /// single-tone `chi3` Kerr + TPA + FCA/FCD (with `power_ramp` and `arclength` bistability continuation)
    /// and multi-harmonic `chi2` / `chi3_xpm` (`shg`, `sfg`, `dfg`, `opa`, `depleted` or `undepleted`).
    NonlinearDriven(crate::nonlinear::NonlinearDrivenSpec),
    /// 2D cross-section fixed-power nonlinear waveguide eigenmodes (`"nonlinear_modes"`, ADR 143) computing
    /// vectorial nonlinear parameter `gamma` (`1/(W*m)`), `A_eff` (`m^2`), and power-dependent `n_eff(P)`.
    NonlinearModes(crate::nonlinear::NonlinearModesSpec),
    /// High-order curved/mixed-cell Discontinuous Galerkin Time-Domain solver (`"dgtd"`, ADR 144)
    /// with power-of-two Local Time-Stepping (LTS).
    Dgtd(crate::dgtd::DgtdConfig),
    /// 2D/3D SIMD structured subpixel FDTD (`"fdtd"`, ADR 145) on a non-uniform Yee grid with conformal
    /// dielectric/PEC averaging, Drude/Lorentz ADEs, CPML/Absorber/Bloch boundaries, FEM-seeded `ModeSource` /
    /// `ModeMonitor`, DFT/time monitors, and Prony `ResonanceFinder`.
    Fdtd(Box<crate::fdtd::FdtdStudyConfig>),
}

impl Study {
    /// The `study.kind` of the envelope.
    pub fn kind(&self) -> &'static str {
        match self {
            Study::Modes { .. } => "modes",
            Study::AxisymmetricModes { .. } => "axisymmetric_modes",
            Study::Cavity { .. } => "cavity",
            Study::Bloch { .. } => "bloch",
            Study::Driven { .. } => "driven",
            Study::Eme { .. } => "eme",
            Study::NonlinearDriven(_) => "nonlinear_driven",
            Study::NonlinearModes(_) => "nonlinear_modes",
            Study::Dgtd(_) => "dgtd",
            Study::Fdtd(..) => "fdtd",
        }
    }

    /// Reference frequency of the study (Hz): `freq` for `Modes`, `Bloch`, `Eme`, `NonlinearDriven`, and `NonlinearModes`, `f_search` for `AxisymmetricModes` and `Cavity`, midband for `Driven`, `Dgtd`, and `Fdtd`.
    pub fn reference_freq(&self) -> f64 {
        match *self {
            Study::Modes { freq, .. } | Study::Bloch { freq, .. } | Study::Eme { freq, .. } => freq,
            Study::NonlinearDriven(ref s) => s.freq,
            Study::NonlinearModes(ref s) => s.freq,
            Study::AxisymmetricModes { f_search, .. } | Study::Cavity { f_search, .. } => f_search,
            Study::Driven { f_min, f_max, ref freqs_hz, .. } => {
                if let Some(fs) = freqs_hz.as_ref().filter(|v| !v.is_empty()) {
                    0.5 * (fs[0] + fs[fs.len() - 1])
                } else {
                    0.5 * (f_min + f_max)
                }
            }
            Study::Dgtd(ref cfg) => cfg.reference_freq(),
            Study::Fdtd(ref cfg) => cfg.reference_freq_hz(),
        }
    }
}

/// Allowed values for `surface_impedance` `ibc_mode`.
pub const IBC_MODE_CHOICES: [&str; 3] = ["direct", "perturbative", "auto"];

/// Allowed values for `Study::Driven` `solver`.
pub const DRIVEN_SOLVER_CHOICES: [&str; 5] = ["auto", "direct", "lblt", "blr", "ams_fgmres"];

/// Allowed values for `Study::Fdtd` and `Study::Dgtd` `backend` (`ADR 146`).
pub const TD_BACKEND_CHOICES: [&str; 3] = ["auto", "cpu", "gpu"];

/// Allowed values for `Study::Fdtd` and `Study::Dgtd` floating-point `precision` (`ADR 146`).
pub const TD_PRECISION_CHOICES: [&str; 3] = ["auto", "f32", "f64"];

/// Execution backend selection for explicit time-domain solvers (`Study::Fdtd` and `Study::Dgtd`, `ADR 146`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TdBackend {
    /// Automatically dispatches to `wgpu` when compiled with `--features gpu`, a discrete GPU is available,
    /// the problem size meets the offload threshold, and all active physics features are GPU-supported;
    /// otherwise executes on the deterministic `f64` CPU kernel.
    #[default]
    Auto,
    /// Always executes on the deterministic `f64` CPU kernel (`0 ULP` across thread counts).
    Cpu,
    /// Explicitly requests `wgpu` GPU compute acceleration (returns an error if compiled without
    /// `--features gpu`, or falls back to CPU with a `DiagnosticWarning` if an unsupported physics
    /// feature is active).
    Gpu,
}

impl TdBackend {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "auto" | "" => Ok(Self::Auto),
            "cpu" => Ok(Self::Cpu),
            "gpu" => Ok(Self::Gpu),
            other => Err(format!("unknown time-domain backend '{other}' (expected one of {TD_BACKEND_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Cpu => "cpu",
            Self::Gpu => "gpu",
        }
    }
}

/// Floating-point arithmetic precision selection for time-domain GPU kernels (`ADR 146`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TdPrecision {
    /// Uses `f64` on GPU when nonlinear media are active and the adapter supports `SHADER_F64`;
    /// otherwise uses scaled `f32` (`E_ref` / energy normalization).
    #[default]
    Auto,
    /// Uses scaled `f32` arithmetic on GPU (`E_ref`-normalized when nonlinear media are active).
    F32,
    /// Requires native IEEE-754 `f64` GPU compute (`wgpu::Features::SHADER_F64`), falling back to
    /// CPU `f64` with a `gpu_fallback_to_cpu` warning if the adapter lacks `SHADER_F64`.
    F64,
}

impl TdPrecision {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.trim() {
            "auto" | "" => Ok(Self::Auto),
            "f32" => Ok(Self::F32),
            "f64" => Ok(Self::F64),
            other => Err(format!("unknown time-domain precision '{other}' (expected one of {TD_PRECISION_CHOICES:?})")),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::F32 => "f32",
            Self::F64 => "f64",
        }
    }
}

/// Kinetic-inductance / surface-impedance sheet boundary condition (`Z_s = R_s + i omega L_k_sq`).
#[derive(Clone, Debug, PartialEq)]
pub struct SurfaceImpedanceBc {
    pub labels: Vec<i32>,
    pub l_k_sq: f64,
    pub r_s: f64,
    pub material: String,
    pub thickness: f64,
    pub temperature: f64,
    pub f_ref: f64,
    pub film_model: String,
    pub lossless: bool,
    pub ibc_mode: String,
}

/// Lumped Josephson junction / nonlinear port boundary condition (`cavity` and `driven` studies).
#[derive(Clone, Debug, PartialEq)]
pub struct LumpedJjBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub l_j: f64,
    pub e_j_hz: f64,
    pub c_j: f64,
    pub e_l_hz: f64,
    pub phi_ext_rad: f64,
    pub n_g: f64,
    pub squid_asymmetry: Option<f64>,
    pub direction: [f64; 3],
    pub formulation: String,
}

/// Lumped sheet port boundary condition (`cavity` and `driven` studies).
#[derive(Clone, Debug, PartialEq)]
pub struct LumpedPortBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub z0: f64,
    pub direction: [f64; 3],
    pub complex_damping: bool,
    pub mode: String,
}

/// 2D boundary-mode waveguide port boundary condition (`driven` study).
#[derive(Clone, Debug, PartialEq)]
pub struct WavePortBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub mode: usize,
    pub n_eff_guess: Option<f64>,
    pub z0_ohm: Option<f64>,
}

/// Zero-thickness 2D material sheet boundary/interface condition (`"sheet_2d"`, ADR 137).
///
/// Assembles the penetrable tangential surface conductivity jump form
/// $$a_{2\text{D}}(\mathbf{w}, \mathbf{E}) = +i\omega\mu_0 \int_{\Gamma_{2\text{D}}} \mathbf{w}_t \cdot (\boldsymbol{\sigma}_{2\text{D}}(\mathbf{x})\,\mathbf{E}_t)\,\mathrm{d}S$$
/// on deduplicated interior/exterior facets with zero auxiliary DOFs and no Dirichlet clamp on $\mathbf{E}_t$.
#[derive(Clone, Debug, PartialEq)]
pub struct Sheet2dBc {
    pub labels: Vec<i32>,
    /// Optional domain-label pairs `(domains_a, domains_b)` when `bc.on` names an interior region interface
    /// (`"a|b"` or `"interface(a, b)"`) whose interior 1D edges are not emitted into `SimplexMesh::facets` by the 2D mesher.
    pub interface_domains: Vec<(Vec<i32>, Vec<i32>)>,
    /// Underlying 2D material sheet model (loaded from `[sheet_2d]` in a library material or constructed from inline BC parameters).
    pub sheet: femcore::material::Sheet2d,
    /// Optional explicit inline principal-axis complex conductivity override `[[su_re, su_im], [sv_re, sv_im]]` in Siemens.
    pub inline_sigma_uv: Option<[[f64; 2]; 2]>,
    /// Gate oxide capacitance per unit area $C_{\text{ox}}$ in $\text{F/m}^2$ ($0$ when electrostatically ungated).
    pub c_ox: f64,
    /// Charge-neutrality-point gate potential offset $V_{\text{CNP}}$ in Volts.
    pub v_cnp: f64,
    /// Uniform gate potential bias $V_{\text{gate}}$ in Volts used when no spatially varying `"potential"` field is imported.
    pub v_gate: f64,
    /// Optional reference frequency $f_{\text{ref}}$ in $\text{Hz}$ ($0$ uses the study reference frequency).
    pub f_ref: f64,
}

impl Sheet2dBc {
    /// Whether this 2D sheet has a non-symmetric conductivity tensor ($\sigma_{uv} \neq \sigma_{vu}$, e.g. Hall bias $B_0 \neq 0$).
    pub fn is_non_symmetric(&self) -> bool {
        self.inline_sigma_uv.is_none() && self.sheet.is_non_symmetric()
    }

    /// True when the Hermitian part $\frac{1}{2}(\boldsymbol{\sigma}_{2\text{D}} + \boldsymbol{\sigma}_{2\text{D}}^H)$ is non-zero
    /// (or electrostatic gating $C_{\text{ox}} > 0$ can activate conductive loss).
    pub fn has_loss(&self, freq_hz: f64) -> bool {
        if self.c_ox > 0.0 {
            return true;
        }
        self.eval_tensor_at(freq_hz, None, None).ok().is_some_and(|s| {
            let h00 = s[0][0].re;
            let h11 = s[1][1].re;
            let h01 = femcore::linalg::csr::c(0.5, 0.0) * (s[0][1] + s[1][0].conj());
            h00.abs() > 1e-15 || h11.abs() > 1e-15 || h01.norm() > 1e-15
        })
    }

    /// Evaluates the $2 \times 2$ tangential surface conductivity tensor $\boldsymbol{\sigma}_{2\text{D}}$ in Siemens (`S`)
    /// at frequency `freq_hz` ($\text{Hz}$), optional local temperature `local_temp_k` ($\text{K}$),
    /// and optional local electrostatic potential `local_potential_v` ($\text{V}$, e.g. interpolated from `dfs.V`).
    pub fn eval_tensor_at(
        &self,
        freq_hz: f64,
        local_temp_k: Option<f64>,
        local_potential_v: Option<f64>,
    ) -> Result<[[femcore::linalg::csr::C; 2]; 2], String> {
        use femcore::linalg::csr::c;
        let f_eval = if self.f_ref > 0.0 { self.f_ref } else { freq_hz };
        let omega = 2.0 * std::f64::consts::PI * f_eval.max(1.0);

        if let Some([[su_re, su_im], [sv_re, sv_im]]) = self.inline_sigma_uv {
            let su = c(su_re, su_im);
            let sv = c(sv_re, sv_im);
            let (s, co) = self.sheet.theta_rad.sin_cos();
            let c2 = co * co;
            let s2 = s * s;
            let cs = co * s;
            let s_xx = su * c2 + sv * s2;
            let s_yy = su * s2 + sv * c2;
            let s_xy = (su - sv) * cs;
            return Ok([[s_xx, s_xy], [s_xy, s_yy]]);
        }

        let mu_c_override = if self.c_ox > 0.0 {
            let v_eff = local_potential_v.unwrap_or(self.v_gate);
            let v_f = match &self.sheet.model {
                femcore::material::Sheet2dModel::KuboGraphene { v_f, .. } => *v_f,
                _ => femcore::material::DEFAULT_GRAPHENE_VF,
            };
            let (_n_s, mu_c_j) =
                femcore::material::chemical_potential_from_gate_potential(v_eff, self.v_cnp, self.c_ox, v_f)?;
            Some(mu_c_j)
        } else {
            None
        };

        let t = self.sheet.eval_tensor(omega, local_temp_k, mu_c_override, None)?;
        Ok([
            [c(t[0][0].re, t[0][0].im), c(t[0][1].re, t[0][1].im)],
            [c(t[1][0].re, t[1][0].im), c(t[1][1].re, t[1][1].im)],
        ])
    }
}

/// Floquet diffraction port on a planar unit-cell face (`driven` study, FR-36): Rayleigh harmonics `|m| <= m_max`,
/// `|n| <= n_max` (`None` = automatic: every propagating order at `f_max` plus one evanescent ring).
#[derive(Clone, Debug, PartialEq)]
pub struct DiffractionPortBc {
    pub labels: Vec<i32>,
    pub name: String,
    pub m_max: Option<i32>,
    pub n_max: Option<i32>,
}

/// Per-domain temperature-dependent optical/microwave permittivity evaluator resolved from the material library
/// for spatially varying photothermal `temperature` imports (ADR 139, FR-38).
#[derive(Clone, Debug)]
pub struct ThermoOpticDomain {
    pub spec: std::sync::Arc<femcore::material::MaterialSpec>,
    pub euler_deg: [f64; 3],
    pub wavelength_m: Option<f64>,
    pub frequency_hz: Option<f64>,
}

impl PartialEq for ThermoOpticDomain {
    fn eq(&self, other: &Self) -> bool {
        self.spec.name == other.spec.name
            && self.euler_deg == other.euler_deg
            && self.wavelength_m == other.wavelength_m
            && self.frequency_hz == other.frequency_hz
            && self.spec.to_toml() == other.spec.to_toml()
    }
}

impl ThermoOpticDomain {
    /// Evaluates global-frame `(eps_re, eps_im_loss)` at `temp_k` (K).
    pub fn eval_at(&self, temp_k: f64) -> Result<Option<(femcore::material::M3, femcore::material::M3)>, String> {
        let tk = Some(temp_k);
        if let (Some(_), Some(lam)) = (&self.spec.optical, self.wavelength_m.filter(|&l| l <= 3e-4)) {
            let eps_g = self.spec.optical_permittivity(lam, tk, self.euler_deg).map_err(|e| e.to_string())?;
            let re = [
                [eps_g[0][0].re, eps_g[0][1].re, eps_g[0][2].re],
                [eps_g[1][0].re, eps_g[1][1].re, eps_g[1][2].re],
                [eps_g[2][0].re, eps_g[2][1].re, eps_g[2][2].re],
            ];
            let im = [
                [-eps_g[0][0].im, -eps_g[0][1].im, -eps_g[0][2].im],
                [-eps_g[1][0].im, -eps_g[1][1].im, -eps_g[1][2].im],
                [-eps_g[2][0].im, -eps_g[2][1].im, -eps_g[2][2].im],
            ];
            return Ok(Some((re, im)));
        }
        if let Some(pair) = self.spec.complex_microwave_permittivity(self.frequency_hz, tk, None, false, self.euler_deg)
        {
            return Ok(Some(pair));
        }
        Ok(None)
    }

    pub fn has_loss(&self) -> bool {
        let t_ref = self.spec.reference_temperature.unwrap_or(293.15);
        self.eval_at(t_ref).ok().flatten().is_some_and(|(re, im)| is_lossy_eps_pair(&re, &im))
    }
}

/// True when `eps = re - i * im` is non-Hermitian (`eps != eps^H`), i.e. has a non-zero dissipative anti-Hermitian part
/// (distinguishing true electromagnetic absorption from lossless Hermitian gyrotropy where `re` is symmetric and `im` is skew-symmetric).
#[inline]
pub fn is_lossy_eps_pair(re: &femcore::material::M3, im: &femcore::material::M3) -> bool {
    for a in 0..3 {
        for b in 0..3 {
            if (im[a][b] + im[b][a]).abs() > 1e-15 || (re[a][b] - re[b][a]).abs() > 1e-15 {
                return true;
            }
        }
    }
    false
}

/// Equality-comparable wrapper around [`nyx_physics::TransferredField`] for cross-mesh slice sampling (`ADR 135`, `ADR 140`).
#[derive(Clone, Debug)]
pub struct ImportedTransfer(pub nyx_physics::TransferredField);

impl PartialEq for ImportedTransfer {
    fn eq(&self, other: &Self) -> bool {
        self.0.physics == other.0.physics
            && self.0.field == other.0.field
            && self.0.order == other.0.order
            && self.0.components == other.0.components
            && self.0.values == other.0.values
            && self.0.domain_names == other.0.domain_names
            && femcore::transfer::same_exchange_mesh(&self.0.mesh, &other.0.mesh)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Problem {
    /// element order, 1 or 2
    pub order: usize,
    pub study: Study,
    /// `eps_r` and `mu_r` per domain (scalars, or crystal tensors in `aniso`), absorbing layers in `pml`
    pub materials: Materials,
    /// Floquet-periodic boundary pairs (cavity study)
    pub periodic: Vec<Periodic>,
    /// boundary numbers of magnetic walls (cavity study): perfect magnetic conductor, the natural condition of the curl-curl form, so
    /// the field is left free on them instead of having its tangential part set to zero. A mirror plane that keeps the modes even
    /// under the reflection (tangential E even, E normal odd) is one of these.
    pub pmc: Vec<i32>,
    /// boundary numbers of explicit PEC walls or interior 2D PEC sheets
    pub pec: Vec<i32>,
    /// superconducting / resistive surface-impedance boundary conditions
    pub surface_impedance: Vec<SurfaceImpedanceBc>,
    /// lumped Josephson junction boundary conditions (cavity / driven study)
    pub lumped_jj: Vec<LumpedJjBc>,
    /// lumped port boundary conditions (cavity / driven study)
    pub lumped_ports: Vec<LumpedPortBc>,
    /// waveguide port boundary conditions (driven / eme study)
    pub waveports: Vec<WavePortBc>,
    /// zero-thickness 2D material sheet boundary/interface conditions
    pub sheet_2d: Vec<Sheet2dBc>,
    /// optional spatially varying electrostatic potential field V(r) (Volts) imported from `nyx-dfs` (`"potential"`)
    pub imported_potential: Option<Vec<f64>>,
    /// Floquet diffraction ports (driven study, FR-36)
    pub diffraction_ports: Vec<DiffractionPortBc>,
    /// scattered-field analytical excitations (driven study, FR-36)
    pub excitations: Vec<crate::excitation::ExcitationSpec>,
    /// optional 1-based domain label -> region name mapping from `ProblemSpec::materials`
    pub domain_names: BTreeMap<i32, String>,
    /// optional midplane PEC/PMC symmetry condition on `x = 0` and/or `z = 0`
    pub symmetry: Symmetry,
    /// global-frame real dispersive permittivity derivative `Re(d(omega * eps_r) / d(omega))` per domain label
    /// when evaluated from a material's `[optical]` dispersion model; domains absent here fall back to static `Re(eps_r)`
    pub d_omega_eps: BTreeMap<i32, [[f64; 3]; 3]>,
    /// optional `[interfaces.{ma,ms,sa}]` TLS surface-loss parameters resolved from the material library
    pub tls_interfaces: Option<TlsInterfaces>,
    /// optional spatially varying nodal or cellwise temperature field (K) imported from `dfs.T`
    pub temperature_field: Option<Vec<f64>>,
    /// optional domain-filtered per-cell temperature vector (K) populated by cross-mesh `Ctx::transferred_inputs` (`ADR 135`)
    pub cell_temperatures: Option<Vec<f64>>,
    /// optional un-resampled cross-mesh source temperature field (`ADR 135` / `ADR 140`) for exact 2D EME slice sampling
    pub transferred_temperature: Option<ImportedTransfer>,
    /// optional un-resampled cross-mesh source electrostatic potential field (`ADR 135` / `ADR 140`) for 2D EME slice `sheet_2d` gating
    pub transferred_potential: Option<ImportedTransfer>,
    /// optional uniform fallback / reference temperature (K) from `Ctx::temperature` or `ProblemSpec::params["temperature"]`
    pub fallback_temperature: Option<f64>,
    /// per-domain library thermo-optic / temperature-dependent dielectric evaluators
    pub thermo_optic_domains: BTreeMap<i32, ThermoOpticDomain>,
    /// optional spatially varying nodal DC electric field (V/m) imported from `dfs.E` or `sms.E` (`"electric_field"`)
    pub electric_field: Option<Vec<f64>>,
    /// optional domain-filtered per-cell DC electric field `[Ex, Ey, Ez]` (V/m) in wafer axes
    pub cell_electric_field: Option<Vec<[f64; 3]>>,
    /// optional un-resampled cross-mesh source DC electric field (`ADR 135` / `ADR 140` / `ADR 141`) for 2D EME slice EO sampling
    pub transferred_electric_field: Option<ImportedTransfer>,
    /// optional spatially varying nodal mechanical strain field (`m/m`, 6-component solver Voigt) imported from `sms.S` (`"strain"`)
    pub strain_field: Option<Vec<f64>>,
    /// optional domain-filtered per-cell mechanical strain `[Sxx, Syy, Szz, Sxy, Syz, Sxz]` (`m/m`, solver Voigt) in wafer axes
    pub cell_strain: Option<Vec<[f64; 6]>>,
    /// optional un-resampled cross-mesh source mechanical strain field (`ADR 135` / `ADR 140` / `ADR 141`) for 2D EME slice PE sampling
    pub transferred_strain: Option<ImportedTransfer>,
}

impl Problem {
    pub fn new(study: Study, order: usize) -> Problem {
        Problem {
            order,
            study,
            materials: Materials::default(),
            periodic: vec![],
            pmc: vec![],
            pec: vec![],
            surface_impedance: vec![],
            lumped_jj: vec![],
            lumped_ports: vec![],
            waveports: vec![],
            sheet_2d: vec![],
            imported_potential: None,
            diffraction_ports: vec![],
            excitations: vec![],
            domain_names: BTreeMap::new(),
            symmetry: Symmetry::None,
            d_omega_eps: BTreeMap::new(),
            tls_interfaces: None,
            temperature_field: None,
            cell_temperatures: None,
            transferred_temperature: None,
            transferred_potential: None,
            fallback_temperature: None,
            thermo_optic_domains: BTreeMap::new(),
            electric_field: None,
            cell_electric_field: None,
            transferred_electric_field: None,
            strain_field: None,
            cell_strain: None,
            transferred_strain: None,
        }
    }

    /// True when physical dielectric, surface-impedance, or 2D sheet electromagnetic dissipation is active.
    pub fn has_loss(&self) -> bool {
        let f_ref = self.study.reference_freq();
        self.materials.has_complex_permittivity()
            || self.surface_impedance.iter().any(|b| b.r_s > 0.0)
            || self.thermo_optic_domains.values().any(ThermoOpticDomain::has_loss)
            || self.sheet_2d.iter().any(|s| s.has_loss(f_ref))
    }

    /// True when a spatially varying imported `temperature_field` (or cross-mesh `cell_temperatures` / `transferred_temperature`) drives at least one thermo-optic domain.
    pub fn has_spatially_varying_temp(&self) -> bool {
        (self.temperature_field.is_some() || self.cell_temperatures.is_some() || self.transferred_temperature.is_some())
            && !self.thermo_optic_domains.is_empty()
    }

    /// True when an imported DC `electric_field` or mechanical `strain` perturbs at least one electro-optic (`r^S`) or photoelastic (`p`) domain.
    pub fn has_eo_or_pe_perturbation(&self) -> bool {
        let has_eo = (self.electric_field.is_some()
            || self.cell_electric_field.is_some()
            || self.transferred_electric_field.is_some())
            && self.materials.aniso.values().any(|a| a.electro_optic.is_some());
        let has_pe = (self.strain_field.is_some() || self.cell_strain.is_some() || self.transferred_strain.is_some())
            && self.materials.aniso.values().any(|a| a.photoelastic.is_some());
        has_eo || has_pe
    }

    /// True when any cell-wise permittivity perturbation (thermo-optic, electro-optic, or photoelastic) is active.
    pub fn has_cell_material_perturbation(&self) -> bool {
        self.has_spatially_varying_temp() || self.has_eo_or_pe_perturbation()
    }

    /// Pre-rotates crystal-frame electro-optic `r^S_Ik` (`6x3`) and photoelastic `p_IJ` (`6x6`) tensors
    /// into the global wafer frame (in textbook Voigt order `(xx, yy, zz, yz, xz, xy)`) per domain label.
    pub fn domain_eo_pe_global(&self) -> BTreeMap<i32, (Option<[[f64; 3]; 6]>, Option<femcore::material::M6>)> {
        let mut out = BTreeMap::new();
        for (&dom, an) in &self.materials.aniso {
            if an.electro_optic.is_some() || an.photoelastic.is_some() {
                let r_g = an
                    .electro_optic
                    .as_ref()
                    .map(|r| femcore::material::coupled::electro_optic_global(r, an.euler_deg));
                let p_g =
                    an.photoelastic.as_ref().map(|p| femcore::material::coupled::photoelastic_global(p, an.euler_deg));
                out.insert(dom, (r_g, p_g));
            }
        }
        out
    }

    /// Evaluates cell-averaged DC electric field `[Ex, Ey, Ez]` (V/m) in wafer axes across `nc` cells with vertex accessor `cell_verts`.
    pub fn eval_cell_electric_field_with<'a>(
        &self,
        nc: usize,
        nv: usize,
        dim: usize,
        mut cell_verts: impl FnMut(usize) -> &'a [u32],
    ) -> Option<Vec<[f64; 3]>> {
        if let Some(ref ce) = self.cell_electric_field {
            if ce.len() == nc {
                return Some(ce.clone());
            }
        }
        let ef = self.electric_field.as_ref()?;
        if ef.is_empty() {
            return None;
        }
        let mut out = vec![[0.0; 3]; nc];
        if ef.len() == 3 * nc {
            for ci in 0..nc {
                out[ci] = [ef[3 * ci], ef[3 * ci + 1], ef[3 * ci + 2]];
            }
            return Some(out);
        }
        if ef.len() == 2 * nc && dim == 2 {
            for ci in 0..nc {
                out[ci] = [ef[2 * ci], 0.0, ef[2 * ci + 1]];
            }
            return Some(out);
        }
        let n_nodes = if ef.len() % 3 == 0 && (dim == 3 || ef.len() / 3 >= nv && ef.len() % 2 != 0) {
            3
        } else if dim == 2 && ef.len() >= 2 * nv && ef.len() % 2 == 0 {
            2
        } else if ef.len() >= 3 * nv && ef.len() % 3 == 0 {
            3
        } else {
            return None;
        };
        for ci in 0..nc {
            let cv = cell_verts(ci);
            let inv_n = 1.0 / (cv.len().max(1) as f64);
            let mut acc = [0.0; 3];
            for &v in cv {
                let idx = v as usize;
                if n_nodes == 2 {
                    if 2 * idx + 1 < ef.len() {
                        acc[0] += ef[2 * idx];
                        acc[2] += ef[2 * idx + 1];
                    }
                } else if 3 * idx + 2 < ef.len() {
                    acc[0] += ef[3 * idx];
                    acc[1] += ef[3 * idx + 1];
                    acc[2] += ef[3 * idx + 2];
                }
            }
            out[ci] = [acc[0] * inv_n, acc[1] * inv_n, acc[2] * inv_n];
        }
        Some(out)
    }

    /// Evaluates cell-averaged mechanical strain `[Sxx, Syy, Szz, Sxy, Syz, Sxz]` (m/m, solver Voigt) in wafer axes across `nc` cells.
    pub fn eval_cell_strain_with<'a>(
        &self,
        nc: usize,
        nv: usize,
        mut cell_verts: impl FnMut(usize) -> &'a [u32],
    ) -> Option<Vec<[f64; 6]>> {
        if let Some(ref cs) = self.cell_strain {
            if cs.len() == nc {
                return Some(cs.clone());
            }
        }
        let sf = self.strain_field.as_ref()?;
        if sf.is_empty() {
            return None;
        }
        let mut out = vec![[0.0; 6]; nc];
        if sf.len() == 6 * nc {
            for ci in 0..nc {
                out[ci] = [sf[6 * ci], sf[6 * ci + 1], sf[6 * ci + 2], sf[6 * ci + 3], sf[6 * ci + 4], sf[6 * ci + 5]];
            }
            return Some(out);
        }
        if sf.len() >= 6 * nv && sf.len() % 6 == 0 {
            for ci in 0..nc {
                let cv = cell_verts(ci);
                let inv_n = 1.0 / (cv.len().max(1) as f64);
                let mut acc = [0.0; 6];
                for &v in cv {
                    let idx = v as usize;
                    if 6 * idx + 5 < sf.len() {
                        for k in 0..6 {
                            acc[k] += sf[6 * idx + k];
                        }
                    }
                }
                for k in 0..6 {
                    out[ci][k] = acc[k] * inv_n;
                }
            }
            return Some(out);
        }
        None
    }

    /// Evaluates unstretched global-frame `(eps_re, eps_im_loss)` per cell of `mesh`, applying local cell-averaged
    /// temperatures from `self.cell_temperatures` or `self.temperature_field` on any domain in `self.thermo_optic_domains`,
    /// and electro-optic / photoelastic perturbations `\Delta \boldsymbol{\varepsilon}_r = -\boldsymbol{\varepsilon}_r \Delta(\boldsymbol{\varepsilon}_r^{-1}) \boldsymbol{\varepsilon}_r`.
    pub fn cell_complex_permittivity(
        &self,
        mesh: &Discretisation,
    ) -> Result<Vec<(femcore::material::M3, femcore::material::M3)>, String> {
        let mut tens = self.materials.cell_complex_permittivity(mesh, Some(1.0))?;
        let nc = mesh.n_cells();
        let nv = mesh.n_points();
        if !self.thermo_optic_domains.is_empty() {
            if let Some(ct) = self.cell_temperatures.as_ref().filter(|v| v.len() == nc) {
                for ci in 0..nc {
                    let dom = mesh.cell_labels[ci];
                    let Some(td) = self.thermo_optic_domains.get(&dom) else { continue };
                    let tc = ct[ci];
                    if tc > 0.0 && tc.is_finite() {
                        if let Some(pair) = td.eval_at(tc)? {
                            tens[ci] = pair;
                        }
                    }
                }
            } else if let Some(tf) = self.temperature_field.as_ref().filter(|v| !v.is_empty()) {
                for ci in 0..nc {
                    let dom = mesh.cell_labels[ci];
                    let Some(td) = self.thermo_optic_domains.get(&dom) else { continue };
                    let tc = if tf.len() == nc {
                        tf[ci]
                    } else if tf.len() >= nv {
                        let cv = mesh.cell(ci);
                        cv.iter().map(|&v| tf[v as usize]).sum::<f64>() / (cv.len() as f64)
                    } else {
                        continue;
                    };
                    if tc > 0.0 && tc.is_finite() {
                        if let Some(pair) = td.eval_at(tc)? {
                            tens[ci] = pair;
                        }
                    }
                }
            }
        }
        if self.has_eo_or_pe_perturbation() {
            let eo_pe = self.domain_eo_pe_global();
            let cell_e = if let Some(ce) = self.cell_electric_field.as_ref().filter(|v| v.len() == nc) {
                Some(ce.clone())
            } else if let Some(ref ef) = self.electric_field {
                let n1 = nv;
                let n2 = mesh.n_dofs(2);
                let (ord, ncomp) = if ef.len() == 2 * n2 && mesh.dim == 2 {
                    (2, 2)
                } else if ef.len() == 3 * n2 {
                    (2, 3)
                } else if ef.len() == 2 * n1 && mesh.dim == 2 {
                    (1, 2)
                } else if ef.len() == 3 * n1 {
                    (1, 3)
                } else {
                    (0, 0)
                };
                if ord > 0 {
                    let nloc = mesh.nloc(ord);
                    let mut dofs = Vec::new();
                    let mut out = vec![[0.0; 3]; nc];
                    for ci in 0..nc {
                        mesh.cell_dofs(ci, ord, &mut dofs);
                        let mut acc = [0.0; 3];
                        for &gd in dofs.iter().take(nloc) {
                            if ncomp == 2 {
                                acc[0] += ef[2 * gd];
                                acc[2] += ef[2 * gd + 1];
                            } else {
                                acc[0] += ef[3 * gd];
                                acc[1] += ef[3 * gd + 1];
                                acc[2] += ef[3 * gd + 2];
                            }
                        }
                        let inv_n = 1.0 / (nloc as f64);
                        out[ci] = [acc[0] * inv_n, acc[1] * inv_n, acc[2] * inv_n];
                    }
                    Some(out)
                } else {
                    self.eval_cell_electric_field_with(nc, nv, mesh.dim, |ci| mesh.cell(ci))
                }
            } else {
                None
            };
            let cell_s = if let Some(cs) = self.cell_strain.as_ref().filter(|v| v.len() == nc) {
                Some(cs.clone())
            } else if let Some(ref sf) = self.strain_field {
                let n1 = nv;
                let n2 = mesh.n_dofs(2);
                let ord = if sf.len() == 6 * n2 {
                    2
                } else if sf.len() == 6 * n1 {
                    1
                } else {
                    0
                };
                if ord > 0 {
                    let nloc = mesh.nloc(ord);
                    let mut dofs = Vec::new();
                    let mut out = vec![[0.0; 6]; nc];
                    for ci in 0..nc {
                        mesh.cell_dofs(ci, ord, &mut dofs);
                        let mut acc = [0.0; 6];
                        for &gd in dofs.iter().take(nloc) {
                            for k in 0..6 {
                                acc[k] += sf[6 * gd + k];
                            }
                        }
                        let inv_n = 1.0 / (nloc as f64);
                        for k in 0..6 {
                            out[ci][k] = acc[k] * inv_n;
                        }
                    }
                    Some(out)
                } else {
                    self.eval_cell_strain_with(nc, nv, |ci| mesh.cell(ci))
                }
            } else {
                None
            };
            for ci in 0..nc {
                let dom = mesh.cell_labels[ci];
                let Some((r_g, p_g)) = eo_pe.get(&dom) else { continue };
                let e_dc = cell_e.as_ref().map(|v| &v[ci]);
                let s_dc = cell_s.as_ref().map(|v| &v[ci]);
                let d_eps =
                    femcore::material::coupled::delta_eps_eo_pe(&tens[ci].0, r_g.as_ref(), p_g.as_ref(), e_dc, s_dc);
                for a in 0..3 {
                    for b in 0..3 {
                        tens[ci].0[a][b] += d_eps[a][b];
                    }
                }
            }
        }
        Ok(tens)
    }

    /// Populates domain-aware cross-mesh `cell_temperatures` and quadratic `imported_potential` from
    /// `ctx.transferred_inputs` when the exporting physics was solved on a different mesh (`ADR 135`).
    /// When `ctx.transferred_inputs` is empty or `is_same_mesh(mesh)` is true, this is a zero-cost no-op.
    pub(crate) fn apply_transferred_inputs(&mut self, mesh: &nyx_physics::Mesh, ctx: &Ctx) -> Result<(), Error> {
        if self.fallback_temperature.is_none() {
            self.fallback_temperature = ctx.temperature;
        }
        if ctx.transferred_inputs.is_empty() {
            return Ok(());
        }
        let mixed = nyx_physics::discretise_mixed(mesh)?;
        let use_mixed = !mixed.is_pure_simplex() || (mesh.dim() == 2 && !mixed.edge_nodes.is_empty());

        if let Some(tf) = ctx.transferred_inputs.get("temperature") {
            if !tf.is_same_mesh(mesh) {
                self.transferred_temperature = Some(ImportedTransfer(tf.clone()));
                let src_disc_owned;
                let src_disc: &Discretisation = match &tf.discretisation {
                    Some(d) => d.as_ref(),
                    None => {
                        src_disc_owned = nyx_physics::discretise(&tf.mesh)?;
                        &src_disc_owned
                    }
                };
                let src_idx = femcore::transfer::TransferIndex::new(
                    femcore::transfer::SourceMesh::Simplex(src_disc),
                    &tf.domain_names,
                    femcore::transfer::TransferOptions::default(),
                )
                .map_err(|e| Error::from_kernel(e.to_string()))?;
                let src_field = femcore::transfer::SourceField::Nodal {
                    values: &tf.values,
                    ncomp: tf.components.max(1),
                    order: tf.order,
                };
                let fallback_t = ctx.temperature.unwrap_or(293.15);
                if use_mixed {
                    let nc = mixed.n_cells();
                    let dim = mixed.dim;
                    let mut ct = Vec::with_capacity(nc);
                    for ci in 0..nc {
                        let dom = mixed.cell_labels[ci];
                        let cv = mixed.cell(ci);
                        let mut xc = [0.0_f64; 3];
                        for &v in cv {
                            let p = mixed.point(v as usize);
                            for k in 0..dim {
                                xc[k] += p[k];
                            }
                        }
                        let inv_n = 1.0 / (cv.len().max(1) as f64);
                        for xk in xc.iter_mut().take(dim) {
                            *xk *= inv_n;
                        }
                        let tc = if let Some(dname) = self.domain_names.get(&dom) {
                            let ev = src_idx
                                .eval_at_points(&src_field, &xc[..dim], &[Some(dname.as_str())])
                                .map_err(|e| Error::from_kernel(e.to_string()))?;
                            if ev.found.first().copied() == Some(true) {
                                ev.values[0]
                            } else {
                                fallback_t
                            }
                        } else {
                            fallback_t
                        };
                        ct.push(tc);
                    }
                    self.cell_temperatures = Some(ct);
                } else {
                    let d = nyx_physics::discretise(mesh)?;
                    let nc = d.n_cells();
                    let nloc1 = d.nloc(1);
                    let sampled = tf.sample_cell_nodes(&d, &self.domain_names, 1)?;
                    let mut ct = Vec::with_capacity(nc);
                    for ci in 0..nc {
                        let base = ci * nloc1;
                        let mut sum = 0.0;
                        let mut cnt = 0usize;
                        for a in 0..nloc1 {
                            if sampled.found[base + a] {
                                sum += sampled.values[base + a];
                                cnt += 1;
                            }
                        }
                        ct.push(if cnt > 0 { sum / (cnt as f64) } else { fallback_t });
                    }
                    self.cell_temperatures = Some(ct);
                }
            }
        }

        if let Some(tf) = ctx.transferred_inputs.get("potential") {
            if !tf.is_same_mesh(mesh) && !self.sheet_2d.is_empty() {
                self.transferred_potential = Some(ImportedTransfer(tf.clone()));
                let src_disc_owned;
                let src_disc: &Discretisation = match &tf.discretisation {
                    Some(d) => d.as_ref(),
                    None => {
                        src_disc_owned = nyx_physics::discretise(&tf.mesh)?;
                        &src_disc_owned
                    }
                };
                let src_idx = femcore::transfer::TransferIndex::new(
                    femcore::transfer::SourceMesh::Simplex(src_disc),
                    &tf.domain_names,
                    femcore::transfer::TransferOptions::default(),
                )
                .map_err(|e| Error::from_kernel(e.to_string()))?;
                let src_field = femcore::transfer::SourceField::Nodal {
                    values: &tf.values,
                    ncomp: tf.components.max(1),
                    order: tf.order,
                };
                if use_mixed {
                    let nv = mixed.n_points();
                    let dim = mixed.dim;
                    let mut coords = Vec::with_capacity(nv * dim);
                    for v in 0..nv {
                        coords.extend_from_slice(&mixed.point(v)[..dim]);
                    }
                    let ev = src_idx
                        .eval_at_points(&src_field, &coords, &[])
                        .map_err(|e| Error::from_kernel(e.to_string()))?;
                    self.imported_potential = Some(ev.values);
                } else {
                    let d = nyx_physics::discretise(mesh)?;
                    let coords = d.dof_coordinates(2);
                    let ev = src_idx
                        .eval_at_points(&src_field, &coords, &[])
                        .map_err(|e| Error::from_kernel(e.to_string()))?;
                    let mut sheet_verts = std::collections::BTreeSet::new();
                    let mut seen_e = std::collections::BTreeSet::new();
                    for sh in &self.sheet_2d {
                        if d.dim == 2 {
                            for (va, vb) in crate::ibc::collect_sheet2d_edges_2d_simplex(&d, sh, &mut seen_e) {
                                sheet_verts.insert(va as usize);
                                sheet_verts.insert(vb as usize);
                            }
                        } else {
                            for f in 0..d.n_facets() {
                                if sh.labels.contains(&d.facet_labels[f]) {
                                    for &v in d.facet(f) {
                                        sheet_verts.insert(v as usize);
                                    }
                                }
                            }
                        }
                    }
                    if let Some(&bad_v) = sheet_verts.iter().find(|&&v| !ev.found.get(v).copied().unwrap_or(false)) {
                        return Err(Error::invalid(format!(
                            "the potential imported onto 'ems' does not cover sheet_2d vertex {bad_v} (check `active_domains` of the electrostatic mesh)"
                        )));
                    }
                    self.imported_potential = Some(ev.values);
                }
            }
        }

        if let Some(tf) = ctx.transferred_inputs.get("electric_field") {
            if !tf.is_same_mesh(mesh) {
                self.transferred_electric_field = Some(ImportedTransfer(tf.clone()));
                let src_disc_owned;
                let src_disc: &Discretisation = match &tf.discretisation {
                    Some(d) => d.as_ref(),
                    None => {
                        src_disc_owned = nyx_physics::discretise(&tf.mesh)?;
                        &src_disc_owned
                    }
                };
                let src_idx = femcore::transfer::TransferIndex::new(
                    femcore::transfer::SourceMesh::Simplex(src_disc),
                    &tf.domain_names,
                    femcore::transfer::TransferOptions::default(),
                )
                .map_err(|e| Error::from_kernel(e.to_string()))?;
                let ncomp = tf.components.max(1);
                let src_field = femcore::transfer::SourceField::Nodal { values: &tf.values, ncomp, order: tf.order };
                if use_mixed {
                    let nc = mixed.n_cells();
                    let dim = mixed.dim;
                    let mut ce = vec![[0.0; 3]; nc];
                    for ci in 0..nc {
                        let dom = mixed.cell_labels[ci];
                        if !self.materials.aniso.get(&dom).is_some_and(|a| a.electro_optic.is_some()) {
                            continue;
                        }
                        let Some(dname) = self.domain_names.get(&dom) else { continue };
                        let cv = mixed.cell(ci);
                        let mut pts = Vec::with_capacity((cv.len() + 1) * dim);
                        let mut xc = [0.0; 3];
                        for &v in cv {
                            let p = mixed.point(v as usize);
                            pts.extend_from_slice(&p[..dim]);
                            for k in 0..dim {
                                xc[k] += p[k];
                            }
                        }
                        for xk in xc.iter_mut().take(dim) {
                            *xk /= cv.len().max(1) as f64;
                        }
                        pts.extend_from_slice(&xc[..dim]);
                        let npts = cv.len() + 1;
                        let dom_filters = vec![Some(dname.as_str()); npts];
                        let ev = src_idx
                            .eval_at_points(&src_field, &pts, &dom_filters)
                            .map_err(|e| Error::from_kernel(e.to_string()))?;
                        let mut acc = [0.0; 3];
                        let mut cnt = 0usize;
                        for p_idx in 0..npts {
                            if ev.found.get(p_idx).copied() == Some(true) {
                                let base = p_idx * ncomp;
                                if ncomp == 2 && dim == 2 {
                                    acc[0] += ev.values[base];
                                    acc[2] += ev.values[base + 1];
                                } else if ncomp >= 3 {
                                    acc[0] += ev.values[base];
                                    acc[1] += ev.values[base + 1];
                                    acc[2] += ev.values[base + 2];
                                }
                                cnt += 1;
                            }
                        }
                        if cnt > 0 {
                            let inv = 1.0 / (cnt as f64);
                            ce[ci] = [acc[0] * inv, acc[1] * inv, acc[2] * inv];
                        }
                    }
                    self.cell_electric_field = Some(ce);
                } else {
                    let d = nyx_physics::discretise(mesh)?;
                    let nc = d.n_cells();
                    let nloc2 = d.nloc(2);
                    let sampled = tf.sample_cell_nodes(&d, &self.domain_names, 2)?;
                    let mut ce = vec![[0.0; 3]; nc];
                    for ci in 0..nc {
                        let base_pt = ci * nloc2;
                        let mut acc = [0.0; 3];
                        let mut cnt = 0usize;
                        for a in 0..nloc2 {
                            let pt_idx = base_pt + a;
                            if sampled.found[pt_idx] {
                                let v_base = pt_idx * ncomp;
                                if ncomp == 2 && d.dim == 2 {
                                    acc[0] += sampled.values[v_base];
                                    acc[2] += sampled.values[v_base + 1];
                                } else if ncomp >= 3 {
                                    acc[0] += sampled.values[v_base];
                                    acc[1] += sampled.values[v_base + 1];
                                    acc[2] += sampled.values[v_base + 2];
                                }
                                cnt += 1;
                            }
                        }
                        if cnt > 0 {
                            let inv = 1.0 / (cnt as f64);
                            ce[ci] = [acc[0] * inv, acc[1] * inv, acc[2] * inv];
                        }
                    }
                    self.cell_electric_field = Some(ce);
                }
            }
        }

        if let Some(tf) = ctx.transferred_inputs.get("strain") {
            if !tf.is_same_mesh(mesh) {
                self.transferred_strain = Some(ImportedTransfer(tf.clone()));
                let src_disc_owned;
                let src_disc: &Discretisation = match &tf.discretisation {
                    Some(d) => d.as_ref(),
                    None => {
                        src_disc_owned = nyx_physics::discretise(&tf.mesh)?;
                        &src_disc_owned
                    }
                };
                let src_idx = femcore::transfer::TransferIndex::new(
                    femcore::transfer::SourceMesh::Simplex(src_disc),
                    &tf.domain_names,
                    femcore::transfer::TransferOptions::default(),
                )
                .map_err(|e| Error::from_kernel(e.to_string()))?;
                let ncomp = tf.components.max(1);
                let src_field = femcore::transfer::SourceField::Nodal { values: &tf.values, ncomp, order: tf.order };
                if use_mixed {
                    let nc = mixed.n_cells();
                    let dim = mixed.dim;
                    let mut cs = vec![[0.0; 6]; nc];
                    for ci in 0..nc {
                        let dom = mixed.cell_labels[ci];
                        if !self.materials.aniso.get(&dom).is_some_and(|a| a.photoelastic.is_some()) {
                            continue;
                        }
                        let Some(dname) = self.domain_names.get(&dom) else { continue };
                        let cv = mixed.cell(ci);
                        let mut pts = Vec::with_capacity((cv.len() + 1) * dim);
                        let mut xc = [0.0; 3];
                        for &v in cv {
                            let p = mixed.point(v as usize);
                            pts.extend_from_slice(&p[..dim]);
                            for k in 0..dim {
                                xc[k] += p[k];
                            }
                        }
                        for xk in xc.iter_mut().take(dim) {
                            *xk /= cv.len().max(1) as f64;
                        }
                        pts.extend_from_slice(&xc[..dim]);
                        let npts = cv.len() + 1;
                        let dom_filters = vec![Some(dname.as_str()); npts];
                        let ev = src_idx
                            .eval_at_points(&src_field, &pts, &dom_filters)
                            .map_err(|e| Error::from_kernel(e.to_string()))?;
                        let mut acc = [0.0; 6];
                        let mut cnt = 0usize;
                        for p_idx in 0..npts {
                            if ev.found.get(p_idx).copied() == Some(true) && ncomp >= 6 {
                                let base = p_idx * ncomp;
                                for k in 0..6 {
                                    acc[k] += ev.values[base + k];
                                }
                                cnt += 1;
                            }
                        }
                        if cnt > 0 {
                            let inv = 1.0 / (cnt as f64);
                            for k in 0..6 {
                                cs[ci][k] = acc[k] * inv;
                            }
                        }
                    }
                    self.cell_strain = Some(cs);
                } else {
                    let d = nyx_physics::discretise(mesh)?;
                    let nc = d.n_cells();
                    let nloc2 = d.nloc(2);
                    let sampled = tf.sample_cell_nodes(&d, &self.domain_names, 2)?;
                    let mut cs = vec![[0.0; 6]; nc];
                    for ci in 0..nc {
                        let base_pt = ci * nloc2;
                        let mut acc = [0.0; 6];
                        let mut cnt = 0usize;
                        for a in 0..nloc2 {
                            let pt_idx = base_pt + a;
                            if sampled.found[pt_idx] && ncomp >= 6 {
                                let v_base = pt_idx * ncomp;
                                for k in 0..6 {
                                    acc[k] += sampled.values[v_base + k];
                                }
                                cnt += 1;
                            }
                        }
                        if cnt > 0 {
                            let inv = 1.0 / (cnt as f64);
                            for k in 0..6 {
                                cs[ci][k] = acc[k] * inv;
                            }
                        }
                    }
                    self.cell_strain = Some(cs);
                }
            }
        }
        Ok(())
    }

    /// Global-frame unstretched `Re(d(omega * eps_r) / d(omega))` tensor for every cell of `mesh`,
    /// falling back to static `Re(eps_r)` on domains without an `[optical]` dispersion entry in `d_omega_eps`.
    pub fn cell_d_omega_eps(&self, mesh: &Discretisation) -> Result<Vec<[[f64; 3]; 3]>, String> {
        let static_re = self.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?;
        Ok((0..mesh.n_cells())
            .map(|ci| {
                let dom = mesh.cell_labels[ci];
                self.d_omega_eps.get(&dom).copied().unwrap_or(static_re[ci])
            })
            .collect())
    }

    /// Every number the user can type must be finite: NaN or infinity does not stop a factorisation, it flows into every field.
    pub fn check_finite(&self) -> Result<(), String> {
        let m = &self.materials;
        for (name, v) in &m.default {
            if !v.is_finite() {
                return Err(format!(
                    "material property '{name}' (default for all domains) is {v}: it must be a finite number"
                ));
            }
        }
        for (dom, props) in &m.domains {
            for (name, v) in props {
                if !v.is_finite() {
                    return Err(format!(
                        "material property '{name}' of domain {dom} is {v}: it must be a finite number"
                    ));
                }
            }
        }
        let study_bad = match &self.study {
            Study::Modes { freq, n_target, .. } => (!freq.is_finite())
                .then_some("frequency")
                .or((!n_target.is_finite()).then_some("target effective index")),
            Study::Bloch { freq, n_target, .. } => {
                (!freq.is_finite()).then_some("frequency").or((!n_target.is_finite()).then_some("target Bloch index"))
            }
            Study::AxisymmetricModes { f_search, .. } | Study::Cavity { f_search, .. } => {
                (!f_search.is_finite()).then_some("search frequency")
            }
            Study::Driven { f_min, f_max, freqs_hz, prom_tol, seed_freqs_hz, blr_tol, .. } => (!f_min.is_finite())
                .then_some("minimum frequency")
                .or((!f_max.is_finite()).then_some("maximum frequency"))
                .or(prom_tol.is_some_and(|t| !t.is_finite()).then_some("PROM tolerance"))
                .or((!blr_tol.is_finite()).then_some("BLR tolerance"))
                .or(freqs_hz.as_ref().is_some_and(|v| v.iter().any(|f| !f.is_finite())).then_some("freqs_hz"))
                .or(seed_freqs_hz.iter().any(|f| !f.is_finite()).then_some("seed_freqs_hz")),
            Study::Eme { freq, n_target, grid, y_range, port_z0_ohm, .. } => (!freq.is_finite())
                .then_some("frequency")
                .or((!n_target.is_finite()).then_some("target effective index"))
                .or((!port_z0_ohm.is_finite()).then_some("port reference impedance"))
                .or(y_range.is_some_and(|r| !r[0].is_finite() || !r[1].is_finite()).then_some("y_range"))
                .or((!grid.is_finite()).then_some("EME grid specification")),
            Study::NonlinearDriven(s) => (!s.is_finite()).then_some("nonlinear_driven study parameters"),
            Study::NonlinearModes(s) => (!s.is_finite()).then_some("nonlinear_modes study parameters"),
            Study::Dgtd(ref cfg) => (!cfg.is_finite()).then_some("DGTD configuration"),
            Study::Fdtd(cfg) => {
                cfg.check_finite()?;
                None
            }
        };
        if let Some(what) = study_bad {
            return Err(format!("the study's {what} is NaN or infinite"));
        }
        for ibc in &self.surface_impedance {
            for (k, v) in [
                ("L_k_sq", ibc.l_k_sq),
                ("R_s", ibc.r_s),
                ("thickness", ibc.thickness),
                ("temperature", ibc.temperature),
                ("f_ref", ibc.f_ref),
            ] {
                if !v.is_finite() {
                    return Err(format!("surface_impedance parameter '{k}' is {v}: it must be a finite number"));
                }
            }
        }
        for jj in &self.lumped_jj {
            for (k, v) in [
                ("L_j", jj.l_j),
                ("E_j_Hz", jj.e_j_hz),
                ("C_j", jj.c_j),
                ("E_l_Hz", jj.e_l_hz),
                ("phi_ext_rad", jj.phi_ext_rad),
                ("n_g", jj.n_g),
                ("direction[0]", jj.direction[0]),
                ("direction[1]", jj.direction[1]),
                ("direction[2]", jj.direction[2]),
            ] {
                if !v.is_finite() {
                    return Err(format!("lumped_jj '{}' parameter '{k}' is {v}: it must be a finite number", jj.name));
                }
            }
            if let Some(d) = jj.squid_asymmetry {
                if !d.is_finite() {
                    return Err(format!(
                        "lumped_jj '{}' parameter 'squid_asymmetry' is {d}: it must be a finite number",
                        jj.name
                    ));
                }
            }
        }
        for lp in &self.lumped_ports {
            for (k, v) in [
                ("Z_0", lp.z0),
                ("direction[0]", lp.direction[0]),
                ("direction[1]", lp.direction[1]),
                ("direction[2]", lp.direction[2]),
            ] {
                if !v.is_finite() {
                    return Err(format!(
                        "lumped_port '{}' parameter '{k}' is {v}: it must be a finite number",
                        lp.name
                    ));
                }
            }
        }
        for wp in &self.waveports {
            if let Some(ng) = wp.n_eff_guess {
                if !ng.is_finite() {
                    return Err(format!(
                        "waveport '{}' parameter 'n_eff_guess' is {ng}: it must be a finite number",
                        wp.name
                    ));
                }
            }
            if let Some(z0) = wp.z0_ohm {
                if !z0.is_finite() {
                    return Err(format!(
                        "waveport '{}' parameter 'z0_ohm' is {z0}: it must be a finite number",
                        wp.name
                    ));
                }
            }
        }
        for sh in &self.sheet_2d {
            for (k, v) in [
                ("theta_rad", sh.sheet.theta_rad),
                ("c_ox", sh.c_ox),
                ("v_cnp", sh.v_cnp),
                ("v_gate", sh.v_gate),
                ("f_ref", sh.f_ref),
            ] {
                if !v.is_finite() {
                    return Err(format!("sheet_2d parameter '{k}' is {v}: it must be a finite number"));
                }
            }
            if let Some(uv) = &sh.inline_sigma_uv {
                if uv.iter().flatten().any(|v| !v.is_finite()) {
                    return Err("sheet_2d inline conductivity contains NaN or infinity".into());
                }
            }
            match &sh.sheet.model {
                femcore::material::Sheet2dModel::KuboGraphene { mu_c, tau, temperature, v_f, b_0, n_layers } => {
                    for (k, v) in [
                        ("mu_c", *mu_c),
                        ("tau", *tau),
                        ("temperature", *temperature),
                        ("v_f", *v_f),
                        ("b_0", *b_0),
                        ("n_layers", *n_layers),
                    ] {
                        if !v.is_finite() {
                            return Err(format!("sheet_2d parameter '{k}' is {v}: it must be a finite number"));
                        }
                    }
                }
                femcore::material::Sheet2dModel::Anisotropic { sigma_uu_0, sigma_vv_0, poles_u, poles_v, .. } => {
                    if !sigma_uu_0[0].is_finite()
                        || !sigma_uu_0[1].is_finite()
                        || !sigma_vv_0[0].is_finite()
                        || !sigma_vv_0[1].is_finite()
                    {
                        return Err("sheet_2d baseline conductivity contains NaN or infinity".into());
                    }
                    for p in poles_u.iter().chain(poles_v.iter()) {
                        if !p.omega_0.is_finite() || !p.gamma.is_finite() || !p.weight.is_finite() {
                            return Err("sheet_2d pole contains NaN or infinity".into());
                        }
                    }
                }
                femcore::material::Sheet2dModel::Tabulated {
                    omega,
                    sigma_uu_re,
                    sigma_uu_im,
                    sigma_vv_re,
                    sigma_vv_im,
                    ..
                } => {
                    if omega.iter().any(|v| !v.is_finite())
                        || sigma_uu_re.iter().any(|v| !v.is_finite())
                        || sigma_uu_im.iter().any(|v| !v.is_finite())
                        || sigma_vv_re.iter().any(|v| !v.is_finite())
                        || sigma_vv_im.iter().any(|v| !v.is_finite())
                    {
                        return Err("sheet_2d tabulated data contains NaN or infinity".into());
                    }
                }
            }
        }
        if let Some(pot) = &self.imported_potential {
            if pot.iter().any(|v| !v.is_finite()) {
                return Err("imported potential field contains NaN or infinity".into());
            }
        }
        // tensors of anisotropic domains, layers and Floquet phases: scan the printed form, outside quoted strings, for NaN and inf
        for (&dom, mat) in &self.d_omega_eps {
            if mat.iter().flatten().any(|v| !v.is_finite()) {
                return Err(format!("dispersive permittivity derivative of domain {dom} contains NaN or infinity"));
            }
        }
        if let Some(ref tf) = self.temperature_field {
            if tf.iter().any(|v| !v.is_finite()) {
                return Err("imported temperature field contains NaN or infinity".into());
            }
        }
        if let Some(ref ct) = self.cell_temperatures {
            if ct.iter().any(|v| !v.is_finite()) {
                return Err("transferred cell temperature field contains NaN or infinity".into());
            }
        }
        if let Some(ref tt) = self.transferred_temperature {
            if tt.0.values.iter().any(|v| !v.is_finite()) {
                return Err("transferred source temperature field contains NaN or infinity".into());
            }
        }
        if let Some(ref tp) = self.transferred_potential {
            if tp.0.values.iter().any(|v| !v.is_finite()) {
                return Err("transferred source potential field contains NaN or infinity".into());
            }
        }
        if let Some(ft) = self.fallback_temperature {
            if !ft.is_finite() {
                return Err("fallback temperature is NaN or infinite".into());
            }
        }
        let text = format!(
            "{:?} {:?} {:?} {:?} {:?} {:?}",
            m.aniso, m.pml, self.periodic, self.tls_interfaces, self.excitations, self.diffraction_ports
        );
        if femcore::validate::non_finite_token(&text).is_some() {
            return Err(
                "the problem description contains a NaN or infinite number (an anisotropic material, a layer, a Floquet phase or an interface)"
                    .into(),
            );
        }
        Ok(())
    }
}

fn pml_dir(text: &str) -> PmlDir {
    match text {
        "plus" => PmlDir::Plus,
        "minus" => PmlDir::Minus,
        _ => PmlDir::Off,
    }
}

/// Reads the envelope. The structural checks against the descriptors (`validate_spec`) are the engine's job; this resolves names and
/// selectors with `ctx.resolver` and refuses combinations that need the study (a Floquet pair outside the cavity study).
pub fn decode(spec: &ProblemSpec, ctx: &Ctx) -> Result<Problem, Error> {
    let ems = Ems;
    let studies = ems.studies();
    let sd = studies
        .iter()
        .find(|s| s.kind == spec.study.kind)
        .ok_or_else(|| Error::Unsupported(format!("ems has no study '{}'", spec.study.kind)))?;
    let sp = Params::new(&spec.study.params, &sd.params);
    let study = match sd.kind.as_str() {
        "modes" => Study::Modes {
            freq: sp.number("freq")?,
            n_target: sp.number("n_target")?,
            evanescent: sp.boolean("evanescent")?,
            n_modes: sp.integer("n_modes")?.max(0) as usize,
        },
        "axisymmetric_modes" => Study::AxisymmetricModes {
            f_search: sp.number("f_search")?,
            azimuthal_order: sp.integer("azimuthal_order")?,
            n_modes: sp.integer("n_modes")?.max(0) as usize,
        },
        "cavity" => Study::Cavity { n_modes: sp.integer("n_modes")?.max(0) as usize, f_search: sp.number("f_search")? },
        "bloch" => {
            let freq = sp.number("freq")?;
            if freq <= 0.0 || !freq.is_finite() {
                return Err(Error::invalid(format!("bloch study requires freq > 0 (got {freq})")));
            }
            Study::Bloch {
                freq,
                n_target: sp.number("n_target")?,
                n_modes: sp.integer("n_modes")?.max(0) as usize,
                brillouin_zone: sp.integer("brillouin_zone")?,
            }
        }
        "driven" => {
            let freqs_hz = sp.optional_numbers("freqs_hz")?;
            let (f_min, f_max, n_points) = if let Some(ref fs) = freqs_hz {
                for &f in fs {
                    if !(f > 0.0 && f.is_finite()) {
                        return Err(Error::invalid(format!(
                            "driven study requires all entries of 'freqs_hz' to be positive finite numbers (got {f})"
                        )));
                    }
                }
                if fs.windows(2).any(|w| w[1] < w[0]) {
                    return Err(Error::invalid("driven study requires 'freqs_hz' to be monotonically ascending"));
                }
                (fs[0], fs[fs.len() - 1], fs.len())
            } else {
                if !sp.has("f_min") || !sp.has("f_max") {
                    return Err(Error::invalid(
                        "driven study requires either non-empty 'freqs_hz' or both 'f_min' and 'f_max'",
                    ));
                }
                let f_min = sp.number("f_min")?;
                let f_max = sp.number("f_max")?;
                if !(f_min > 0.0 && f_max >= f_min) {
                    return Err(Error::invalid(format!(
                        "driven study requires 0 < f_min <= f_max (got f_min = {f_min}, f_max = {f_max})"
                    )));
                }
                let n_points = sp.integer("n_points")?.max(1) as usize;
                (f_min, f_max, n_points)
            };
            let prom_tol = sp.optional_number("prom_tol")?.filter(|&t| t > 0.0);
            let max_hdm_raw = sp.integer("max_hdm_solves")?;
            if !(1..=64).contains(&max_hdm_raw) {
                return Err(Error::invalid(format!(
                    "driven study requires 1 <= max_hdm_solves <= 64 (got {max_hdm_raw})"
                )));
            }
            let max_hdm_solves = max_hdm_raw as usize;
            let strict_prom_tol = sp.boolean("strict_prom_tol")?;
            let seed_freqs_hz = sp.numbers("seed_freqs_hz")?;
            for &sf in &seed_freqs_hz {
                if !(sf > 0.0 && sf.is_finite()) {
                    return Err(Error::invalid(format!(
                        "driven study requires all entries of 'seed_freqs_hz' to be positive finite numbers (got {sf})"
                    )));
                }
            }
            let solver = sp.text("solver")?;
            if !DRIVEN_SOLVER_CHOICES.contains(&solver.as_str()) {
                return Err(Error::invalid(format!(
                    "driven study: unknown solver '{solver}' (expected one of {DRIVEN_SOLVER_CHOICES:?})"
                )));
            }
            let blr_tol = sp.number("blr_tol")?;
            if !blr_tol.is_finite() || blr_tol < 0.0 {
                return Err(Error::invalid(format!("driven study requires finite blr_tol >= 0 (got {blr_tol})")));
            }
            Study::Driven {
                f_min,
                f_max,
                n_points,
                freqs_hz,
                prom_tol,
                max_hdm_solves,
                strict_prom_tol,
                seed_freqs_hz,
                solver,
                blr_tol,
            }
        }
        "eme" => {
            let freq = sp.number("freq")?;
            if !(freq > 0.0 && freq.is_finite()) {
                return Err(Error::invalid(format!("eme study requires freq > 0 (got {freq})")));
            }
            let n_target = sp.number("n_target")?;
            let evanescent = sp.boolean("evanescent")?;
            let n_modes = sp.integer("n_modes")?.max(1) as usize;
            let grid = crate::eme::decode_eme_grid(&sp)?;
            let yr = sp.numbers("y_range")?;
            let y_range = if yr.is_empty() {
                None
            } else if yr.len() == 2 {
                if !(yr[0].is_finite() && yr[1].is_finite() && yr[1] > yr[0]) {
                    return Err(Error::invalid(format!(
                        "eme study requires y_range = [y0, y1] with y0 < y1 (got [{}, {}])",
                        yr[0], yr[1]
                    )));
                }
                Some([yr[0], yr[1]])
            } else {
                return Err(Error::invalid(format!(
                    "eme study requires 'y_range' to be empty or [y0, y1] (got {} entries)",
                    yr.len()
                )));
            };
            let port_modes = sp.integer("port_modes")?.max(0) as usize;
            let port_z0_ohm = sp.number("port_z0_ohm")?;
            if !(port_z0_ohm > 0.0 && port_z0_ohm.is_finite()) {
                return Err(Error::invalid(format!(
                    "eme study requires positive finite port_z0_ohm (got {port_z0_ohm})"
                )));
            }
            let reconstruct_fields = sp.boolean("reconstruct_fields")?;
            Study::Eme {
                freq,
                n_target,
                evanescent,
                n_modes,
                grid,
                y_range,
                port_modes,
                port_z0_ohm,
                reconstruct_fields,
            }
        }
        "nonlinear_driven" => {
            let freq = sp.number("freq")?;
            if !(freq > 0.0 && freq.is_finite()) {
                return Err(Error::invalid(format!("nonlinear_driven study requires freq > 0 (got {freq})")));
            }
            let freq_2 = sp.number("freq_2")?;
            if !freq_2.is_finite() || freq_2 < 0.0 {
                return Err(Error::invalid(format!("nonlinear_driven study requires freq_2 >= 0 (got {freq_2})")));
            }
            let process = crate::nonlinear::NonlinearProcess::parse(&sp.text("process")?).map_err(Error::invalid)?;
            let continuation =
                crate::nonlinear::NonlinearContinuation::parse(&sp.text("continuation")?).map_err(Error::invalid)?;
            let coupling = crate::nonlinear::Chi2CouplingMode::parse(&sp.text("coupling")?).map_err(Error::invalid)?;
            let nl_solver =
                crate::nonlinear::NonlinearSolverKind::parse(&sp.text("nl_solver")?).map_err(Error::invalid)?;
            let power_w = sp.number("power_w")?;
            let powers_w = sp.numbers("powers_w")?;
            let power_min_w = sp.number("power_min_w")?;
            let power_max_w = sp.number("power_max_w")?;
            if power_w < 0.0 || power_min_w < 0.0 || power_max_w < 0.0 || powers_w.iter().any(|&p| p < 0.0) {
                return Err(Error::invalid("nonlinear_driven optical powers must be non-negative"));
            }
            let n_steps = sp.integer("n_steps")?.max(1) as usize;
            let drive_port = sp.integer("drive_port")?.max(0) as usize;
            let port_powers_w = sp.numbers("port_powers_w")?;
            if port_powers_w.iter().any(|&p| p < 0.0) {
                return Err(Error::invalid("nonlinear_driven 'port_powers_w' entries must be non-negative"));
            }
            let port_phases_rad = sp.numbers("port_phases_rad")?;
            let signal_power_w = sp.number("signal_power_w")?;
            if signal_power_w < 0.0 {
                return Err(Error::invalid("nonlinear_driven 'signal_power_w' must be non-negative"));
            }
            let signal_port = sp.integer("signal_port")?.max(0) as usize;
            let signal_phase_rad = sp.number("signal_phase_rad")?;
            let rtol = sp.number("rtol")?;
            if !(rtol > 0.0 && rtol.is_finite()) {
                return Err(Error::invalid(format!("nonlinear_driven requires positive finite rtol (got {rtol})")));
            }
            let max_it = sp.integer("max_it")?.max(1) as usize;
            let arclength_ds = sp.number("arclength_ds")?;
            if !(arclength_ds > 0.0 && arclength_ds.is_finite()) {
                return Err(Error::invalid("nonlinear_driven requires positive finite arclength_ds"));
            }
            let arclength_max_steps = sp.integer("arclength_max_steps")?.max(1) as usize;
            Study::NonlinearDriven(crate::nonlinear::NonlinearDrivenSpec {
                process,
                freq,
                freq_2,
                power_w,
                powers_w,
                power_min_w,
                power_max_w,
                n_steps,
                drive_port,
                port_powers_w,
                port_phases_rad,
                signal_power_w,
                signal_port,
                signal_phase_rad,
                continuation,
                coupling,
                nl_solver,
                rtol,
                max_it,
                arclength_ds,
                arclength_max_steps,
            })
        }
        "nonlinear_modes" => {
            let freq = sp.number("freq")?;
            if !(freq > 0.0 && freq.is_finite()) {
                return Err(Error::invalid(format!("nonlinear_modes study requires freq > 0 (got {freq})")));
            }
            let n_target = sp.number("n_target")?;
            let evanescent = sp.boolean("evanescent")?;
            let mode_index = sp.integer("mode_index")?.max(1) as usize;
            let power_w = sp.number("power_w")?;
            let powers_w = sp.numbers("powers_w")?;
            if power_w < 0.0 || powers_w.iter().any(|&p| p < 0.0) {
                return Err(Error::invalid("nonlinear_modes optical power must be non-negative"));
            }
            let rtol = sp.number("rtol")?;
            if !(rtol > 0.0 && rtol.is_finite()) {
                return Err(Error::invalid(format!("nonlinear_modes requires positive finite rtol (got {rtol})")));
            }
            let max_it = sp.integer("max_it")?.max(1) as usize;
            Study::NonlinearModes(crate::nonlinear::NonlinearModesSpec {
                freq,
                n_target,
                evanescent,
                mode_index,
                power_w,
                powers_w,
                rtol,
                max_it,
            })
        }
        "dgtd" => Study::Dgtd(crate::dgtd::decode_dgtd(&sp, ctx)?),
        "fdtd" => {
            let cfg_txt = sp.text("config_json")?;
            let mut cfg: crate::fdtd::FdtdStudyConfig = if cfg_txt.trim().is_empty() {
                crate::fdtd::FdtdStudyConfig::default()
            } else {
                serde_json::from_str(cfg_txt.trim())
                    .map_err(|e| Error::invalid(format!("fdtd study invalid config_json: {e}")))?
            };
            let has_cfg = !cfg_txt.trim().is_empty();
            if !has_cfg || sp.has("wavelength_min") {
                cfg.grid.auto.wavelength_min_m = sp.number("wavelength_min")?;
            }
            if !has_cfg || sp.has("min_steps_per_wave") {
                cfg.grid.auto.min_steps_per_wave = sp.number("min_steps_per_wave")?;
            }
            if !has_cfg || sp.has("run_time") {
                cfg.run_time_s = sp.number("run_time")?;
            }
            if !has_cfg || sp.has("max_steps") {
                cfg.max_steps = sp.integer("max_steps")?.max(0) as usize;
            }
            if !has_cfg || sp.has("courant") {
                cfg.courant = sp.number("courant")?;
            }
            if !has_cfg || sp.has("shutoff_energy_ratio") {
                cfg.shutoff_energy_ratio = sp.number("shutoff_energy_ratio")?;
            }
            if !has_cfg || sp.has("polarization_2d") {
                cfg.polarization_2d = match sp.text("polarization_2d")?.as_str() {
                    "tm_ey" | "" => crate::fdtd::Pol2d::TmEy,
                    "te_hy" => crate::fdtd::Pol2d::TeHy,
                    "full" => crate::fdtd::Pol2d::Full,
                    other => {
                        return Err(Error::invalid(format!(
                            "fdtd study: unknown polarization_2d '{other}' (expected 'tm_ey', 'te_hy', or 'full')"
                        )))
                    }
                };
            }
            if !has_cfg || sp.has("subpixel") {
                cfg.subpixel = match sp.text("subpixel")?.as_str() {
                    "polarized" | "" => crate::fdtd::SubpixelMode::Polarized,
                    "volumetric" => crate::fdtd::SubpixelMode::Volumetric,
                    "staircased" => crate::fdtd::SubpixelMode::Staircased,
                    other => {
                        return Err(Error::invalid(format!(
                        "fdtd study: unknown subpixel '{other}' (expected 'polarized', 'volumetric', or 'staircased')"
                    )))
                    }
                };
            }
            if !has_cfg || sp.has("conformal_pec") {
                cfg.conformal_pec = sp.boolean("conformal_pec")?;
            }
            if !has_cfg || sp.has("port_z0_ohm") {
                cfg.port_z0_ohm = sp.number("port_z0_ohm")?;
            }
            if !has_cfg || sp.has("backend") {
                let b_str = if sp.has("backend") { sp.text("backend")? } else { "auto".to_string() };
                cfg.backend = TdBackend::parse(&b_str).map_err(Error::invalid)?;
            }
            if !has_cfg || sp.has("precision") {
                let p_str = if sp.has("precision") { sp.text("precision")? } else { "auto".to_string() };
                cfg.precision = TdPrecision::parse(&p_str).map_err(Error::invalid)?;
            }
            if !has_cfg || sp.has("nonlinear") {
                cfg.nonlinear = sp.boolean("nonlinear")?;
            }
            let fs = sp.numbers("freqs_hz")?;
            if !fs.is_empty() && cfg.monitors.is_empty() {
                cfg.monitors.push(crate::fdtd::FdtdMonitorSpec::FieldMonitor {
                    name: "domain_dft".into(),
                    lo: None,
                    hi: None,
                    freqs_hz: fs,
                });
            }
            cfg.check_finite().map_err(Error::invalid)?;
            Study::Fdtd(Box::new(cfg))
        }
        other => return Err(Error::Unsupported(format!("ems has no study '{other}'"))),
    };
    let symmetry = if matches!(
        sd.kind.as_str(),
        "modes" | "cavity" | "bloch" | "driven" | "eme" | "nonlinear_driven" | "nonlinear_modes" | "dgtd" | "fdtd"
    ) {
        Symmetry::parse(&sp.text("symmetry")?).map_err(Error::invalid)?
    } else {
        Symmetry::None
    };

    let f_ref = study.reference_freq();
    let derived_lambda = (f_ref > 0.0 && f_ref.is_finite()).then(|| crate::modes::C0 / f_ref);
    let mut ctx_eff = Ctx::new(ctx.resolver);
    ctx_eff.library = ctx.library.clone();
    ctx_eff.wavelength = ctx.wavelength.or(derived_lambda);
    ctx_eff.frequency = ctx.frequency.or_else(|| (f_ref > 0.0 && f_ref.is_finite()).then_some(f_ref));
    ctx_eff.temperature = ctx.temperature;
    ctx_eff.inputs = ctx.inputs.clone();
    ctx_eff.progress = ctx.progress.clone();

    let mut materials = materials_by_label(spec, &ctx_eff, &ems.properties())?;
    let domain_names = domain_names_by_label(spec, &ctx_eff)?;
    let temp_k = ctx_eff.temperature.or_else(|| {
        spec.params.get("temperature").and_then(serde_json::Value::as_f64).filter(|&t| t > 0.0 && t.is_finite())
    });
    let mut d_omega_eps: BTreeMap<i32, [[f64; 3]; 3]> = BTreeMap::new();
    let mut thermo_optic_domains: BTreeMap<i32, ThermoOpticDomain> = BTreeMap::new();
    if let Some(lib) = ctx_eff.library.as_ref() {
        for (region, assign) in &spec.materials {
            let (name, euler) = match assign {
                MaterialAssign::Named(n) => (n.as_str(), [0.0; 3]),
                MaterialAssign::Oriented { material, euler_deg } => (material.as_str(), *euler_deg),
                _ => continue,
            };
            let labels = ctx_eff.resolver.domains(region)?;
            let mat = lib.require(name).map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
            if mat.optical.is_some() || mat.dielectric.is_some() {
                let td = ThermoOpticDomain {
                    spec: std::sync::Arc::new(mat.clone()),
                    euler_deg: euler,
                    wavelength_m: ctx_eff.wavelength,
                    frequency_hz: ctx_eff.frequency,
                };
                for &l in &labels {
                    thermo_optic_domains.insert(l, td.clone());
                }
            }
            if let (Some(opt), Some(lambda_m)) = (mat.optical.as_ref(), ctx_eff.wavelength.filter(|&lam| lam <= 3e-4)) {
                let mut diag = [0.0_f64; 3];
                for (axis, slot) in diag.iter_mut().enumerate() {
                    let (n, k, dn_dlam, dk_dlam) = opt
                        .nk(axis, lambda_m, temp_k)
                        .map_err(|e| Error::invalid(format!("region '{region}': {e}")))?;
                    let n_tilde = femcore::linalg::csr::c(n, -k);
                    let dn_tilde = femcore::linalg::csr::c(dn_dlam, -dk_dlam);
                    let d_om = n_tilde * n_tilde - n_tilde * dn_tilde * (2.0 * lambda_m);
                    *slot = d_om.re;
                }
                let cryst = [[diag[0], 0.0, 0.0], [0.0, diag[1], 0.0], [0.0, 0.0, diag[2]]];
                let glob = femcore::material::rotate_tensor(&cryst, &femcore::material::rotation(euler));
                for &l in &labels {
                    d_omega_eps.insert(l, glob);
                }
            }
        }
    }
    let temperature_field = ctx_eff.inputs.get("temperature").cloned();
    let electric_field =
        ctx_eff.inputs.get("electric_field").or_else(|| ctx_eff.inputs.get("ems.electric_field")).cloned();
    let strain_field = ctx_eff.inputs.get("strain").or_else(|| ctx_eff.inputs.get("ems.strain")).cloned();

    let mut periodic = vec![];
    let mut pmc: Vec<i32> = vec![];
    let mut pec: Vec<i32> = vec![];
    let mut surface_impedance: Vec<SurfaceImpedanceBc> = vec![];
    let mut lumped_jj: Vec<LumpedJjBc> = vec![];
    let mut lumped_ports: Vec<LumpedPortBc> = vec![];
    let mut waveports: Vec<WavePortBc> = vec![];
    let mut sheet_2d: Vec<Sheet2dBc> = vec![];
    let imported_potential = ctx.inputs.get("potential").or_else(|| ctx.inputs.get("ems.potential")).cloned();
    let mut diffraction_ports: Vec<DiffractionPortBc> = vec![];
    let mut excitations: Vec<crate::excitation::ExcitationSpec> = vec![];
    let kinds = ems.bc_kinds();
    for (bc_spec, r) in spec.bcs.iter().zip(resolve_conditions(&spec.bcs, &kinds, &ctx_eff)?) {
        let d = kinds
            .iter()
            .find(|k| k.kind == r.kind)
            .ok_or_else(|| Error::Unsupported(format!("ems has no condition '{}'", r.kind)))?;
        let p = Params::new(&r.params, &d.params);
        match r.kind.as_str() {
            "pml" => {
                let x_dir = pml_dir(&p.text("x")?);
                let y_dir = pml_dir(&p.text("y")?);
                let z_dir = pml_dir(&p.text("z")?);
                let coord_text = p.text("coord")?;
                let resolved_coord = match coord_text.as_str() {
                    "cartesian" => "cartesian",
                    "cylindrical" => "cylindrical",
                    "spherical" => "spherical",
                    "auto" => {
                        if x_dir != PmlDir::Off || y_dir != PmlDir::Off || z_dir != PmlDir::Off {
                            "cartesian"
                        } else if p.has("r_thickness")
                            || p.has("z_thickness")
                            || p.has("axial")
                            || p.has("radial_on")
                            || p.has("axis")
                        {
                            "cylindrical"
                        } else if p.has("thickness") || p.has("r_in") || p.has("center") {
                            if matches!(study, Study::Modes { .. } | Study::NonlinearModes(_)) {
                                "cylindrical"
                            } else {
                                "spherical"
                            }
                        } else {
                            "cartesian"
                        }
                    }
                    other => {
                        return Err(Error::invalid(format!(
                            "pml: unknown coord '{other}' (expected 'auto', 'cartesian', 'cylindrical', or 'spherical')"
                        )));
                    }
                };
                let center_vec = p.vector("center")?;
                if center_vec.len() != 3 {
                    return Err(Error::invalid("pml: 'center' must have 3 components"));
                }
                let center = [center_vec[0], center_vec[1], center_vec[2]];

                let coord = match resolved_coord {
                    "cartesian" => PmlCoord::Cartesian,
                    "cylindrical" => {
                        let axis_vec = if p.has("axis") {
                            p.vector("axis")?
                        } else if matches!(study, Study::Modes { .. } | Study::NonlinearModes(_)) {
                            vec![0.0, 1.0, 0.0]
                        } else {
                            vec![0.0, 0.0, 1.0]
                        };
                        if axis_vec.len() != 3 {
                            return Err(Error::invalid("pml: 'axis' must have 3 components"));
                        }
                        let anorm =
                            (axis_vec[0] * axis_vec[0] + axis_vec[1] * axis_vec[1] + axis_vec[2] * axis_vec[2]).sqrt();
                        if !(anorm > 0.0 && anorm.is_finite()) {
                            return Err(Error::invalid("pml: 'axis' must be a non-zero finite vector"));
                        }
                        let axis = [axis_vec[0] / anorm, axis_vec[1] / anorm, axis_vec[2] / anorm];
                        let r_in = p.number("r_in")?;
                        let r_thickness =
                            if p.has("r_thickness") { p.number("r_thickness")? } else { p.number("thickness")? };
                        let radial_on = p.boolean("radial_on")?;
                        let axial = pml_dir(&p.text("axial")?);
                        let z_in = p.number("z_in")?;
                        let z_thickness =
                            if p.has("z_thickness") { p.number("z_thickness")? } else { p.number("thickness")? };
                        if r_in < 0.0 || r_thickness < 0.0 || z_thickness < 0.0 {
                            return Err(Error::invalid("pml: cylindrical radii and thicknesses must be non-negative"));
                        }
                        PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness }
                    }
                    "spherical" => {
                        let r_in = p.number("r_in")?;
                        let thickness =
                            if p.has("thickness") { p.number("thickness")? } else { p.number("r_thickness")? };
                        if r_in < 0.0 || thickness < 0.0 {
                            return Err(Error::invalid("pml: spherical r_in and thickness must be non-negative"));
                        }
                        PmlCoord::Spherical { center, r_in, thickness }
                    }
                    other => return Err(Error::invalid(format!("pml: unexpected resolved coord '{other}'"))),
                };

                let profile_text = p.text("profile")?;
                let profile = match profile_text.as_str() {
                    "kappa_polynomial" => {
                        let kappa_max = p.number("kappa_max")?;
                        if !(kappa_max >= 1.0 && kappa_max.is_finite()) {
                            return Err(Error::invalid(format!(
                                "pml: kappa_max must satisfy kappa_max >= 1.0 (got {kappa_max})"
                            )));
                        }
                        PmlProfile::KappaPolynomial { kappa_max }
                    }
                    "rational" => {
                        let curvature = p.number("curvature")?;
                        if !((0.0..1.0).contains(&curvature) && curvature.is_finite()) {
                            return Err(Error::invalid(format!(
                                "pml: curvature must satisfy 0.0 <= curvature < 1.0 (got {curvature})"
                            )));
                        }
                        PmlProfile::Rational { curvature }
                    }
                    other => {
                        return Err(Error::invalid(format!(
                            "pml: unknown profile '{other}' (expected 'kappa_polynomial' or 'rational')"
                        )));
                    }
                };

                let bounds = if p.has("bounds_lo") || p.has("bounds_hi") {
                    if !(p.has("bounds_lo") && p.has("bounds_hi")) {
                        return Err(Error::invalid(
                            "pml: both 'bounds_lo' and 'bounds_hi' must be given when explicit slab bounds are specified",
                        ));
                    }
                    let lo_v = p.vector("bounds_lo")?;
                    let hi_v = p.vector("bounds_hi")?;
                    if lo_v.len() != 3 || hi_v.len() != 3 {
                        return Err(Error::invalid("pml: 'bounds_lo' and 'bounds_hi' must have 3 components"));
                    }
                    let lo = [lo_v[0], lo_v[1], lo_v[2]];
                    let hi = [hi_v[0], hi_v[1], hi_v[2]];
                    if (0..3).any(|k| lo[k] > hi[k]) {
                        return Err(Error::invalid("pml: bounds_lo must be <= bounds_hi on every axis"));
                    }
                    Some((lo, hi))
                } else {
                    None
                };

                let auto_scale = p.boolean("auto_scale")?;
                let target_reflection = p.number("target_reflection")?;
                if !(target_reflection > 0.0 && target_reflection < 1.0 && target_reflection.is_finite()) {
                    return Err(Error::invalid(format!(
                        "pml: target_reflection must lie in (0, 1) (got {target_reflection})"
                    )));
                }

                let layer = Pml {
                    x: x_dir,
                    y: y_dir,
                    z: z_dir,
                    strength: p.number("strength")?,
                    order: p.number("order")?,
                    re: p.number("re")?,
                    coord,
                    bounds,
                    profile,
                    auto_scale,
                    target_reflection,
                    ..Default::default()
                };
                for l in r.labels {
                    materials.pml.insert(l, layer.clone());
                }
            }
            "periodic" => {
                if !matches!(
                    study,
                    Study::Cavity { .. }
                        | Study::Bloch { .. }
                        | Study::Driven { .. }
                        | Study::NonlinearDriven(_)
                        | Study::Dgtd(_)
                        | Study::Fdtd(..)
                ) {
                    return Err(Error::invalid(format!(
                        "a Floquet-periodic pair belongs to the cavity study, not to '{}'",
                        study.kind()
                    )));
                }
                let dest_name = p.text("dest")?;
                let dest = ctx_eff.resolver.boundaries(&dest_name)?;
                if dest.is_empty() {
                    return Err(Error::invalid(format!(
                        "the destination '{dest_name}' of a periodic pair selects nothing"
                    )));
                }
                periodic.push(Periodic { source: r.labels, dest, phase: p.number("phase")? });
            }
            "pmc" => {
                pmc.extend(r.labels);
            }
            "pec" => {
                pec.extend(r.labels);
            }
            "surface_impedance" => {
                let mut l_k_sq = p.number("L_k_sq")?;
                let mut r_s = p.number("R_s")?;
                let material = p.text("material")?;
                let mut thickness = p.number("thickness")?;
                let temperature = p.number("temperature")?;
                let f_ref_given = p.optional_number("f_ref")?.unwrap_or(0.0);
                let f_ref = if f_ref_given > 0.0 { f_ref_given } else { study.reference_freq() };
                let lossless = p.boolean("lossless")?;
                let ibc_mode_raw = p.text("ibc_mode")?;
                let ibc_mode = if ibc_mode_raw.is_empty() { "direct".to_string() } else { ibc_mode_raw };
                if !IBC_MODE_CHOICES.contains(&ibc_mode.as_str()) {
                    return Err(Error::invalid(format!(
                        "surface_impedance: unknown ibc_mode '{ibc_mode}' (expected one of {IBC_MODE_CHOICES:?})"
                    )));
                }
                let film_model = p.text("film_model")?;
                let mut fm = FilmModel::parse(&film_model).ok_or_else(|| {
                    Error::invalid(format!(
                        "surface_impedance: unknown film_model '{film_model}' (expected 'two_sided', 'kerr' or 'thin_limit')"
                    ))
                })?;
                if !material.trim().is_empty() {
                    let lib = ctx.library.as_ref().ok_or_else(|| {
                        Error::Unsupported(format!(
                            "surface_impedance: material '{material}' needs a material library, and this run has none"
                        ))
                    })?;
                    let mat =
                        lib.require(material.trim()).map_err(|e| Error::invalid(format!("surface_impedance: {e}")))?;
                    let sc = mat.superconducting.as_ref().ok_or_else(|| {
                        Error::invalid(format!(
                            "surface_impedance: material '{}' has no [superconducting] section",
                            material.trim()
                        ))
                    })?;
                    if thickness <= 0.0 {
                        thickness = sc.nominal_thickness.unwrap_or(0.0);
                    }
                    if !p.has("film_model") {
                        fm = sc.film_model.unwrap_or(fm);
                    }
                    if f_ref <= 0.0 {
                        return Err(Error::invalid(
                            "surface_impedance with a superconducting material requires a positive f_ref or study search frequency",
                        ));
                    }
                    let eval = bcs_surface_impedance(sc, f_ref, temperature, thickness, fm)
                        .map_err(|e| Error::invalid(format!("surface_impedance ('{}'): {e}", material.trim())))?;
                    l_k_sq = eval.l_k_sq;
                    r_s = eval.r_s;
                }
                if lossless {
                    r_s = 0.0;
                }
                if l_k_sq < 0.0 || r_s < 0.0 {
                    return Err(Error::invalid("surface_impedance requires non-negative L_k_sq and R_s"));
                }
                surface_impedance.push(SurfaceImpedanceBc {
                    labels: r.labels,
                    l_k_sq,
                    r_s,
                    material,
                    thickness,
                    temperature,
                    f_ref,
                    film_model: fm.as_str().to_string(),
                    lossless,
                    ibc_mode,
                });
            }
            "lumped_jj" => {
                if !matches!(study, Study::Cavity { .. } | Study::Driven { .. } | Study::Dgtd(_)) {
                    return Err(Error::invalid(format!(
                        "a lumped_jj condition belongs to the cavity or driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name")?;
                if name.trim().is_empty() {
                    return Err(Error::invalid("lumped_jj requires a non-empty 'name'"));
                }
                if lumped_jj.iter().any(|j| j.name == name)
                    || lumped_ports.iter().any(|lp| lp.name == name)
                    || waveports.iter().any(|wp| wp.name == name)
                    || diffraction_ports.iter().any(|dp| dp.name == name)
                {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let l_j_in = p.number("L_j")?;
                let e_j_in = p.number("E_j_Hz")?;
                let c_j = p.number("C_j")?;
                let e_l_hz = p.number("E_l_Hz")?;
                let phi_ext_rad = p.number("phi_ext_rad")?;
                let n_g = p.number("n_g")?;
                let squid_asymmetry = if p.has("squid_asymmetry") { Some(p.number("squid_asymmetry")?) } else { None };
                if let Some(d) = squid_asymmetry {
                    if !d.is_finite() || !(0.0..1.0).contains(&d) {
                        return Err(Error::invalid(format!(
                            "lumped_jj '{name}': squid_asymmetry must satisfy 0.0 <= d < 1.0 (got {d})"
                        )));
                    }
                }
                let (l_j, e_j_hz) = resolve_jj_l_and_ej(l_j_in, e_j_in);
                if l_j <= 0.0 && e_l_hz <= 0.0 && c_j <= 0.0 {
                    return Err(Error::invalid(format!(
                        "lumped_jj '{name}' requires a positive L_j, E_j_Hz, E_l_Hz, or C_j"
                    )));
                }
                let dir_vec = p.vector("direction")?;
                if dir_vec.len() != 3 {
                    return Err(Error::invalid(format!("lumped_jj '{name}': direction must have 3 components")));
                }
                let dnorm = (dir_vec[0] * dir_vec[0] + dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]).sqrt();
                if !(dnorm > 0.0 && dnorm.is_finite()) {
                    return Err(Error::invalid(format!("lumped_jj '{name}': direction vector must be non-zero")));
                }
                let direction = [dir_vec[0] / dnorm, dir_vec[1] / dnorm, dir_vec[2] / dnorm];
                let formulation = p.text("formulation")?;
                if formulation != "rank1" && formulation != "distributed_sheet" {
                    return Err(Error::invalid(format!(
                        "lumped_jj '{name}': unknown formulation '{formulation}' (expected 'rank1' or 'distributed_sheet')"
                    )));
                }
                lumped_jj.push(LumpedJjBc {
                    labels: r.labels,
                    name,
                    l_j,
                    e_j_hz,
                    c_j,
                    e_l_hz,
                    phi_ext_rad,
                    n_g,
                    squid_asymmetry,
                    direction,
                    formulation,
                });
            }
            "lumped_port" => {
                if !matches!(
                    study,
                    Study::Cavity { .. } | Study::Driven { .. } | Study::NonlinearDriven(_) | Study::Dgtd(_)
                ) {
                    return Err(Error::invalid(format!(
                        "a lumped_port condition belongs to the cavity or driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name")?;
                if name.trim().is_empty() {
                    return Err(Error::invalid("lumped_port requires a non-empty 'name'"));
                }
                if lumped_ports.iter().any(|lp| lp.name == name)
                    || waveports.iter().any(|wp| wp.name == name)
                    || lumped_jj.iter().any(|j| j.name == name)
                    || diffraction_ports.iter().any(|dp| dp.name == name)
                {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let z0 = p.number("Z_0")?;
                if !(z0 > 0.0 && z0.is_finite()) {
                    return Err(Error::invalid(format!("lumped_port '{name}': Z_0 must be positive (got {z0})")));
                }
                let dir_vec = p.vector("direction")?;
                if dir_vec.len() != 3 {
                    return Err(Error::invalid(format!("lumped_port '{name}': direction must have 3 components")));
                }
                let dnorm = (dir_vec[0] * dir_vec[0] + dir_vec[1] * dir_vec[1] + dir_vec[2] * dir_vec[2]).sqrt();
                if !(dnorm > 0.0 && dnorm.is_finite()) {
                    return Err(Error::invalid(format!("lumped_port '{name}': direction vector must be non-zero")));
                }
                let direction = [dir_vec[0] / dnorm, dir_vec[1] / dnorm, dir_vec[2] / dnorm];
                let complex_damping = p.boolean("complex_damping")?;
                let mode = p.text("mode")?;
                if mode != "single" && mode != "cpw_even" {
                    return Err(Error::invalid(format!(
                        "lumped_port '{name}': unknown mode '{mode}' (expected 'single' or 'cpw_even')"
                    )));
                }
                lumped_ports.push(LumpedPortBc { labels: r.labels, name, z0, direction, complex_damping, mode });
            }
            "waveport" => {
                if !matches!(
                    study,
                    Study::Driven { .. }
                        | Study::Eme { .. }
                        | Study::NonlinearDriven(_)
                        | Study::Dgtd(_)
                        | Study::Fdtd(..)
                ) {
                    return Err(Error::invalid(format!(
                        "a waveport condition belongs to the driven or eme study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name")?;
                if name.trim().is_empty() {
                    return Err(Error::invalid("waveport requires a non-empty 'name'"));
                }
                if lumped_ports.iter().any(|lp| lp.name == name)
                    || waveports.iter().any(|wp| wp.name == name)
                    || lumped_jj.iter().any(|j| j.name == name)
                    || diffraction_ports.iter().any(|dp| dp.name == name)
                {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let mode = p.integer("mode")?.max(1) as usize;
                let n_eff_guess = if p.has("n_eff_guess") {
                    let v = p.number("n_eff_guess")?;
                    if !(v > 0.0 && v.is_finite()) {
                        return Err(Error::invalid(format!(
                            "waveport '{name}': n_eff_guess must be positive (got {v})"
                        )));
                    }
                    Some(v)
                } else {
                    None
                };
                let z0_ohm = if p.has("z0_ohm") {
                    let v = p.number("z0_ohm")?;
                    if !(v > 0.0 && v.is_finite()) {
                        return Err(Error::invalid(format!("waveport '{name}': z0_ohm must be positive (got {v})")));
                    }
                    Some(v)
                } else {
                    None
                };
                waveports.push(WavePortBc { labels: r.labels, name, mode, n_eff_guess, z0_ohm });
            }
            "sheet_2d" => {
                let mat_name = p.text("material")?;
                let mat_trim = mat_name.trim();
                let mut inline_sigma_uv = None;

                let mut sheet = if !mat_trim.is_empty() {
                    let lib = ctx.library.as_ref().ok_or_else(|| {
                        Error::Unsupported(format!(
                            "sheet_2d: material '{mat_trim}' needs a material library, and this run has none"
                        ))
                    })?;
                    let mat = lib.require(mat_trim).map_err(|e| Error::invalid(format!("sheet_2d: {e}")))?;
                    let sh = mat.sheet_2d.clone().ok_or_else(|| {
                        Error::invalid(format!("sheet_2d: material '{mat_trim}' has no [sheet_2d] section"))
                    })?;
                    let mut has_unverified = mat.source.is_empty() && sh.source.is_empty();
                    for s in mat.source.iter().chain(sh.source.iter()) {
                        if s.is_unverified() {
                            has_unverified = true;
                        }
                    }
                    if has_unverified {
                        femcore::progress::note(format!(
                            "warning [unverified_material_provenance]: 2D sheet material '{}' has unverified or missing source provenance in [sheet_2d]",
                            mat.name
                        ));
                    }
                    sh
                } else {
                    let model_raw = p.text("model")?;
                    let model_kind = if model_raw.is_empty() { "kubo_graphene" } else { model_raw.as_str() };
                    match model_kind {
                        "kubo_graphene" | "kubo" => {
                            let mu_c_ev = p.number("mu_c_eV")?;
                            let tau_s = p.number("tau_s")?;
                            let temp_bc =
                                if p.has("temperature") { p.number("temperature")? } else { temp_k.unwrap_or(300.0) };
                            let v_f = p.number("v_f")?;
                            let b_0 = p.number("b_0")?;
                            femcore::material::Sheet2d::kubo_graphene(
                                mu_c_ev * femcore::material::EV_TO_J,
                                tau_s,
                                temp_bc,
                                v_f,
                                b_0,
                            )
                            .map_err(Error::invalid)?
                        }
                        "anisotropic" | "constant" => {
                            let su_re = p.number("sigma_uu_re")?;
                            let su_im = p.number("sigma_uu_im")?;
                            let sv_re = if p.has("sigma_vv_re") { p.number("sigma_vv_re")? } else { su_re };
                            let sv_im = if p.has("sigma_vv_im") { p.number("sigma_vv_im")? } else { su_im };
                            if su_re < 0.0 || sv_re < 0.0 {
                                return Err(Error::invalid(format!(
                                    "sheet_2d requires non-negative real sheet conductivity (got sigma_uu_re={su_re}, sigma_vv_re={sv_re})"
                                )));
                            }
                            inline_sigma_uv = Some([[su_re, su_im], [sv_re, sv_im]]);
                            femcore::material::Sheet2d {
                                model: femcore::material::Sheet2dModel::Anisotropic {
                                    sigma_uu_0: [su_re, su_im],
                                    sigma_vv_0: [sv_re, sv_im],
                                    poles_u: Vec::new(),
                                    poles_v: Vec::new(),
                                    frequency_unit: "rad/s".to_string(),
                                },
                                theta_rad: 0.0,
                                t_range: None,
                                source: Vec::new(),
                            }
                        }
                        other => {
                            return Err(Error::invalid(format!(
                                "sheet_2d: unknown model '{other}' (expected 'kubo_graphene', 'anisotropic', or 'constant')"
                            )));
                        }
                    }
                };

                if p.has("theta_deg") {
                    sheet.theta_rad = p.number("theta_deg")?.to_radians();
                }
                if let femcore::material::Sheet2dModel::KuboGraphene { mu_c, tau, temperature, v_f, b_0, n_layers } =
                    &mut sheet.model
                {
                    if p.has("mu_c_eV") {
                        *mu_c = p.number("mu_c_eV")? * femcore::material::EV_TO_J;
                    }
                    if p.has("tau_s") {
                        let ts = p.number("tau_s")?;
                        if ts <= 0.0 {
                            return Err(Error::invalid(format!("sheet_2d: tau_s must be positive (got {ts})")));
                        }
                        *tau = ts;
                    }
                    if p.has("temperature") {
                        let tk = p.number("temperature")?;
                        if tk < 0.0 {
                            return Err(Error::invalid(format!(
                                "sheet_2d: temperature must be non-negative (got {tk})"
                            )));
                        }
                        *temperature = tk;
                    } else if let Some(tk) = temp_k {
                        *temperature = tk;
                    }
                    if p.has("v_f") {
                        let vf = p.number("v_f")?;
                        if vf <= 0.0 {
                            return Err(Error::invalid(format!("sheet_2d: v_f must be positive (got {vf})")));
                        }
                        *v_f = vf;
                    }
                    if p.has("b_0") {
                        *b_0 = p.number("b_0")?;
                    }
                    if p.has("n_layers") {
                        let nl = p.number("n_layers")?;
                        if nl <= 0.0 {
                            return Err(Error::invalid(format!("sheet_2d: n_layers must be positive (got {nl})")));
                        }
                        *n_layers = nl;
                    }
                    let tk_val = *temperature;
                    let _ = tk_val;
                }
                if let femcore::material::Sheet2dModel::KuboGraphene { temperature: tk_val, .. } = &sheet.model {
                    if let Some(clamp) = sheet.check_temperature(*tk_val).clamp {
                        femcore::progress::note(format!(
                            "warning [material_temperature_clamped]: sheet_2d temperature {:.3} K is outside T_range [{:.3}, {:.3}] K and was clamped to {:.3} K",
                            clamp.requested_k, clamp.bounds_k[0], clamp.bounds_k[1], clamp.clamped_k
                        ));
                    }
                }
                if p.has("sigma_uu_re") || p.has("sigma_uu_im") || p.has("sigma_vv_re") || p.has("sigma_vv_im") {
                    let su_re = p.number("sigma_uu_re")?;
                    let su_im = p.number("sigma_uu_im")?;
                    let sv_re = if p.has("sigma_vv_re") { p.number("sigma_vv_re")? } else { su_re };
                    let sv_im = if p.has("sigma_vv_im") { p.number("sigma_vv_im")? } else { su_im };
                    if su_re < 0.0 || sv_re < 0.0 {
                        return Err(Error::invalid(format!(
                            "sheet_2d requires non-negative real sheet conductivity (got sigma_uu_re={su_re}, sigma_vv_re={sv_re})"
                        )));
                    }
                    inline_sigma_uv = Some([[su_re, su_im], [sv_re, sv_im]]);
                }

                let mut c_ox = p.number("c_ox")?;
                let eps_ox = p.number("eps_ox")?;
                let d_ox = p.number("d_ox")?;
                if c_ox < 0.0 {
                    return Err(Error::invalid(format!("sheet_2d: c_ox must be >= 0 (got {c_ox})")));
                }
                if c_ox == 0.0 && (eps_ox > 0.0 || d_ox > 0.0) {
                    c_ox = femcore::material::gate_oxide_capacitance(eps_ox, d_ox).map_err(Error::invalid)?;
                }
                let v_cnp = p.number("v_cnp")?;
                let v_gate = p.number("v_gate")?;
                let f_ref_given = p.optional_number("f_ref")?.unwrap_or(0.0);
                if f_ref_given < 0.0 {
                    return Err(Error::invalid(format!("sheet_2d: f_ref must be >= 0 (got {f_ref_given})")));
                }
                let f_ref = if f_ref_given > 0.0 { f_ref_given } else { study.reference_freq() };

                let interface_domains = parse_interface_domains(&bc_spec.on, ctx_eff.resolver);
                sheet_2d.push(Sheet2dBc {
                    labels: r.labels,
                    interface_domains,
                    sheet,
                    inline_sigma_uv,
                    c_ox,
                    v_cnp,
                    v_gate,
                    f_ref,
                });
            }
            "diffraction_port" => {
                if !matches!(study, Study::Driven { .. } | Study::Dgtd(_)) {
                    return Err(Error::invalid(format!(
                        "a diffraction_port condition belongs to the driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let name = p.text("name")?;
                if name.trim().is_empty() {
                    return Err(Error::invalid("diffraction_port requires a non-empty 'name'"));
                }
                if lumped_ports.iter().any(|lp| lp.name == name)
                    || waveports.iter().any(|wp| wp.name == name)
                    || lumped_jj.iter().any(|j| j.name == name)
                    || diffraction_ports.iter().any(|dp| dp.name == name)
                {
                    return Err(Error::invalid(format!("duplicate port name '{name}'")));
                }
                let order = |key: &str| -> Result<Option<i32>, Error> {
                    let v = p.integer(key)?;
                    if v > 50 {
                        return Err(Error::invalid(format!(
                            "diffraction_port '{name}': {key} must be <= 50 (got {v})"
                        )));
                    }
                    Ok((v >= 0).then_some(v as i32))
                };
                let (m_max, n_max) = (order("m_max")?, order("n_max")?);
                diffraction_ports.push(DiffractionPortBc { labels: r.labels, name, m_max, n_max });
            }
            kind @ ("plane_wave" | "gaussian_beam" | "astigmatic_gaussian_beam" | "thin_lens_beam" | "dipole") => {
                if !matches!(study, Study::Driven { .. } | Study::Dgtd(_)) {
                    return Err(Error::invalid(format!(
                        "a '{kind}' excitation belongs to the driven study, not to '{}'",
                        study.kind()
                    )));
                }
                let spec = decode_excitation(kind, &p, r.labels)?;
                if excitations.iter().any(|e: &crate::excitation::ExcitationSpec| e.name == spec.name) {
                    return Err(Error::invalid(format!("duplicate excitation name '{}'", spec.name)));
                }
                excitations.push(spec);
            }
            other => return Err(Error::Unsupported(format!("ems has no condition '{other}'"))),
        }
    }
    let ifaces = tls_interfaces_from_spec(spec, &ctx_eff);
    let tls_interfaces = (!ifaces.is_empty()).then_some(ifaces);
    pmc.sort_unstable();
    pmc.dedup();
    pec.sort_unstable();
    pec.dedup();
    Ok(Problem {
        order: spec.order as usize,
        study,
        materials,
        periodic,
        pmc,
        pec,
        surface_impedance,
        lumped_jj,
        lumped_ports,
        waveports,
        sheet_2d,
        imported_potential,
        diffraction_ports,
        excitations,
        domain_names,
        symmetry,
        d_omega_eps,
        tls_interfaces,
        temperature_field,
        cell_temperatures: None,
        transferred_temperature: None,
        transferred_potential: None,
        fallback_temperature: temp_k,
        thermo_optic_domains,
        electric_field,
        cell_electric_field: None,
        transferred_electric_field: None,
        strain_field,
        cell_strain: None,
        transferred_strain: None,
    })
}

fn parse_interface_domains(on: &str, resolver: &dyn nyx_physics::LabelResolver) -> Vec<(Vec<i32>, Vec<i32>)> {
    let mut out = Vec::new();
    let s = on.trim();
    if let Some(inner) = s.strip_prefix("interface(").and_then(|rest| rest.strip_suffix(')')) {
        if let Some((a, b)) = inner.split_once(',') {
            if let (Ok(da), Ok(db)) = (resolver.domains(a.trim()), resolver.domains(b.trim())) {
                if !da.is_empty() && !db.is_empty() {
                    out.push((da, db));
                }
            }
        }
    } else if let Some((a, b)) = s.split_once('|') {
        if !b.contains('|') {
            if let (Ok(da), Ok(db)) = (resolver.domains(a.trim()), resolver.domains(b.trim())) {
                if !da.is_empty() && !db.is_empty() {
                    out.push((da, db));
                }
            }
        }
    }
    out
}

/// Reads one scattered-field excitation condition (`plane_wave`, `gaussian_beam`, `astigmatic_gaussian_beam`, `thin_lens_beam`,
/// `dipole`); `labels` are the contrast domains of the condition.
fn decode_excitation(kind: &str, p: &Params, labels: Vec<i32>) -> Result<crate::excitation::ExcitationSpec, Error> {
    use crate::excitation::{BeamModel, Excitation, ExcitationSpec};
    use femcore::linalg::csr::c;
    let name = p.text("name")?;
    if name.trim().is_empty() {
        return Err(Error::invalid(format!("{kind} requires a non-empty 'name'")));
    }
    let v3 = |key: &str| -> Result<[f64; 3], Error> {
        let v = p.vector(key)?;
        if v.len() != 3 || v.iter().any(|x| !x.is_finite()) {
            return Err(Error::invalid(format!("{kind} '{name}': '{key}' must have 3 finite components")));
        }
        Ok([v[0], v[1], v[2]])
    };
    let cv3 = |re: &str, im: &str| -> Result<[femcore::linalg::csr::C; 3], Error> {
        let r = v3(re)?;
        let i = v3(im)?;
        Ok([c(r[0], i[0]), c(r[1], i[1]), c(r[2], i[2])])
    };
    let num = |key: &str| -> Result<f64, Error> {
        let v = p.number(key)?;
        if !v.is_finite() {
            return Err(Error::invalid(format!("{kind} '{name}': '{key}' must be finite (got {v})")));
        }
        Ok(v)
    };
    let eps_bg = c(num("eps_bg")?, num("eps_bg_im")?);
    let mu_bg = num("mu_bg")?;
    if !(eps_bg.re > 0.0 && mu_bg > 0.0) {
        return Err(Error::invalid(format!("{kind} '{name}': eps_bg and mu_bg must be positive")));
    }
    let amplitude = if kind == "dipole" { 1.0 } else { num("amplitude")? };
    let phase = if kind == "dipole" { 0.0 } else { num("phase")? };
    let excitation = match kind {
        "plane_wave" => Excitation::PlaneWave {
            k_hat: v3("direction")?,
            pol: cv3("polarization", "polarization_im")?,
            amplitude,
            phase,
            origin: v3("origin")?,
        },
        "gaussian_beam" => {
            let model = match p.text("model")?.as_str() {
                "complex_source" => BeamModel::ComplexSource,
                "paraxial" => BeamModel::Paraxial,
                other => {
                    return Err(Error::invalid(format!(
                        "gaussian_beam '{name}': unknown model '{other}' (expected 'complex_source' or 'paraxial')"
                    )))
                }
            };
            Excitation::GaussianBeam {
                focus: v3("focus")?,
                w0: num("waist")?,
                k_hat: v3("direction")?,
                pol: cv3("polarization", "polarization_im")?,
                amplitude,
                phase,
                model,
            }
        }
        "astigmatic_gaussian_beam" => Excitation::AstigmaticGaussianBeam {
            focus: v3("focus")?,
            w0_x: num("waist_x")?,
            w0_y: num("waist_y")?,
            focus_x: num("focus_x")?,
            focus_y: num("focus_y")?,
            k_hat: v3("direction")?,
            x_axis: v3("x_axis")?,
            pol: cv3("polarization", "polarization_im")?,
            amplitude,
            phase,
        },
        "thin_lens_beam" => Excitation::ThinLensBeam {
            focus: v3("focus")?,
            na: num("NA")?,
            k_hat: v3("direction")?,
            x_axis: v3("x_axis")?,
            pol: cv3("polarization", "polarization_im")?,
            amplitude,
            phase,
            filling: num("filling")?,
            n_theta: p.integer("n_theta")?.clamp(2, 512) as usize,
            n_phi: p.integer("n_phi")?.clamp(4, 1024) as usize,
        },
        _ => Excitation::Dipole {
            position: v3("position")?,
            moment: cv3("moment", "moment_im")?,
            magnetic: p.boolean("magnetic")?,
        },
    };
    let n_bg = (eps_bg * mu_bg).sqrt().re;
    excitation.validate(n_bg).map_err(|e| Error::invalid(format!("{kind} '{name}': {e}")))?;
    Ok(ExcitationSpec {
        name,
        excitation,
        domains: labels,
        eps_bg,
        mu_bg,
        keep_scattered: p.boolean("keep_scattered")?,
    })
}
