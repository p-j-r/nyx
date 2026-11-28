// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Opt-in `wgpu` compute kernel acceleration for `Study::Dgtd` (`ADR 144`, `ADR 146`, `FR-43`).
//!
//! Dispatches uniform-step (`Lserk4`, `Rk4`, and 1-tier `Lts`) Discontinuous Galerkin Time-Domain
//! electromagnetic simulations to WGSL compute shaders operating on device-resident `f32` or
//! (when `wgpu::Features::SHADER_F64` is supported) native `f64` buffers:
//! - Pass 1 (`eval_face_qps`): evaluates interior, PEC, PMC, Silver–Müller, modal port, TFSF, and
//!   penetrable 2D sheet (`FR-35`) Riemann numerical fluxes across all element face quadrature points.
//! - Pass 2 (`eval_dof_rhs`): computes strong-form SBP volume curls ($D^x_K, D^y_K, D^z_K$), lifts
//!   face quadrature fluxes, injects point dipole sources, evaluates instantaneous nonlinear
//!   $\chi^{(2)}$/QPM/$\chi^{(3)}$/FCD susceptibility tangents and TPA/FCA conductivity losses,
//!   applies constitutive tensors ($\boldsymbol{\varepsilon}_{\text{inst}}^{-1}, \boldsymbol{\mu}^{-1}, \boldsymbol{\sigma}$),
//!   and applies uniaxial PML damping.
//! - Pass 3 (`lserk4_stage` / `rk4_stage`): advances nodal $(\mathbf{E}, \mathbf{H})$ in-place on the GPU.
//! - Pass 3b (`step_carriers`): advances TPA-generated free-carrier density $N_c$ once per macro step.
//! - Pass 4 (`accum_dof_dft`, `sample_probes`, `dgtd_energy_pass1`/`pass2`): accumulates running
//!   discrete Fourier transforms, probe time series, and normalized atomic-free two-level workgroup
//!   tree energy reductions on the device, reading back only final states and monitor outputs.
//!
//! **Future GPU & HPC Roadmap:** `wgpu`/WGSL serves as the zero-SDK portable Tier-1 compute
//! layer today and may be replaced by `CubeCL` or `rust-gpu` in the future once `f64` coverage
//! and stable toolchain compatibility mature (`femcore::linalg::device` remains backend-agnostic).
//! Dedicated HPC support and direct native bindings for NVIDIA (CUDA) and AMD (ROCm/HIP) GPUs
//! are planned for the future; if a user requires these on a priority basis, they should contact
//! the maintainers with target hardware and workload details.
//!
//! Runs requiring hierarchical multi-rate Local Time-Stepping (`n_lts_tiers > 1`), dispersive
//! Drude/Lorentz/gyrotropic ADE poles, or explicit `precision = "f64"` on adapters without
//! `SHADER_F64` return a descriptive fallback reason so [`crate::dgtd::solve_dgtd`] falls back
//! cleanly to CPU with a `DiagnosticWarning` (`code: "gpu_fallback_to_cpu"`).

use std::f64::consts::PI;
use std::time::Instant;

use bytemuck::{Pod, Zeroable};
use femcore::linalg::csr::{c, C};
use femcore::material::HBAR;
use wgpu::util::DeviceExt;

use crate::dgtd::{DgtdAssembly, DgtdFaceKind, DgtdPolarization, DgtdPulse, DgtdSolve, DgtdState, DgtdTimeStepper};
use crate::fdtd::gpu::{select_td_precision, EmsGpuContext, GpuScalar};
use crate::ibc::{EPS0, MU0};

const WG_SIZE_1D: u32 = 64;
const C0_SI: f64 = 299_792_458.0;
const TPA_RMS_SCALE: f64 = 0.816_496_580_927_726;
const DGTD_NC_REF: f64 = 1.0e20;

macro_rules! impl_pod_scalar_struct {
    ($name:ident, $sz32:expr, $sz64:expr) => {
        // SAFETY: all fields are tightly packed scalars/arrays of matching alignment with explicit padding verified by const assert.
        unsafe impl Pod for $name<f32> {}
        // SAFETY: all fields are tightly packed scalars/arrays of matching alignment with explicit padding verified by const assert.
        unsafe impl Pod for $name<f64> {}
        const _: () = assert!(std::mem::size_of::<$name<f32>>() == $sz32);
        const _: () = assert!(std::mem::size_of::<$name<f64>>() == $sz64);
    };
}

/// Per-step / per-stage uniform parameters bound at `@group(0) @binding(0)`.
/// Layout: 20 `R` scalars (`80B` in `f32`, `160B` in `f64`) + 16 `u32` (`64B`), zero implicit padding.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct DgtdStepUniforms<R> {
    dt: R,
    t_now: R,
    pulse_now: R,
    flux_alpha: R,
    rk_a: R,
    rk_b: R,
    rk_w_acc: R,
    rk_h_next: R,
    pulse_f0: R,
    pulse_fwidth: R,
    pulse_t0: R,
    pulse_amp: R,
    pulse_phase: R,
    symmetry_factor: R,
    inv_eps0: R,
    e_ref: R,
    inv_e_ref: R,
    _pad_r0: R,
    _pad_r1: R,
    _pad_r2: R,
    total_dofs: u32,
    total_qps: u32,
    n_cells: u32,
    n_freqs: u32,
    n_probes: u32,
    n_energy_wgs: u32,
    dim: u32,
    pol_mode: u32,
    rk_stage: u32,
    snap_idx: u32,
    pulse_kind: u32,
    has_nonlinear: u32,
    has_carriers: u32,
    _pad_u0: u32,
    _pad_u1: u32,
    _pad_u2: u32,
}
impl_pod_scalar_struct!(DgtdStepUniforms, 144, 224);

/// Per-cell material and energy quadrature metadata (`16 * vec4<R> + 8 * u32`).
/// Notice `eps_si_row*` holds dimensionless `eps_r` and `mu_si_row*` holds `(MU0 / EPS0) * mu_r`,
/// paired with `cell_masses` scaled by `1 / mean_measure` so nanoscale DGTD energies avoid `f32` underflow.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuCellMeta<R> {
    inv_eps_row0: [R; 4],
    inv_eps_row1: [R; 4],
    inv_eps_row2: [R; 4],
    inv_mu_row0: [R; 4],
    inv_mu_row1: [R; 4],
    inv_mu_row2: [R; 4],
    sigma_e_row0: [R; 4],
    sigma_e_row1: [R; 4],
    sigma_e_row2: [R; 4],
    /// `.xyz` = uniaxial PML damping `[s_sum - sx, s_sum - sy, s_sum - sz]`, `.w` = `1.0` if `has_sigma_e`.
    pml_damp: [R; 4],
    /// Dimensionless `eps_r` rows (used for normalized energy and nonlinear `eps_inst = eps_r + d_eps`).
    eps_si_row0: [R; 4],
    eps_si_row1: [R; 4],
    eps_si_row2: [R; 4],
    /// Impedance-scaled permeability `(MU0 / EPS0) * mu_r` rows for normalized energy evaluation.
    mu_si_row0: [R; 4],
    mu_si_row1: [R; 4],
    mu_si_row2: [R; 4],
    dof_offset: u32,
    n_p: u32,
    mass_offset: u32,
    /// `0` = in PML (skip energy), `1` = diagonal mass, `2` = dense `n_p x n_p` mass.
    energy_mode: u32,
    /// Index into `nl_cells_meta`, or `0xffff_ffff` for linear cells.
    nl_cell_slot: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}
impl_pod_scalar_struct!(GpuCellMeta, 288, 544);

/// Per-DOF volume differentiation and face-lift index metadata (`2 * vec4<R> + 8 * u32`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuDofMeta<R> {
    /// `.xyz` = E dipole moment vector; `.w` = local node QPM domain sign (`+1.0` / `-1.0`).
    dipole_e_vec: [R; 4],
    dipole_h_vec: [R; 4],
    cell_id: u32,
    cell_dof_offset: u32,
    n_p: u32,
    dmat_row_offset: u32,
    lift_start: u32,
    lift_len: u32,
    _pad0: u32,
    _pad1: u32,
}
impl_pod_scalar_struct!(GpuDofMeta, 64, 96);

/// One `(lift_w, qp_idx)` pair mapping a face quadrature point flux to a nodal DOF.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuDofLift<R> {
    lift_w: R,
    qp_idx: u32,
    _pad: u32,
}
impl_pod_scalar_struct!(GpuDofLift, 12, 16);

/// One `(weight, global_dof)` trace evaluation pair on $K^-$ or $K^+$.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuTracePair<R> {
    weight: R,
    global_dof: u32,
    _pad: u32,
}
impl_pod_scalar_struct!(GpuTracePair, 12, 16);

/// Per-face-quadrature-point metadata (`9 * vec4<R> + 8 * u32`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuFaceQpMeta<R> {
    n_hat_ym: [R; 4],
    yp_tfsf_delay: [R; 4],
    port_e_mode: [R; 4],
    port_h_mode: [R; 4],
    tfsf_e_unit: [R; 4],
    tfsf_h_unit: [R; 4],
    sheet_sigma: [R; 4],
    sheet_u_hat: [R; 4],
    sheet_v_hat: [R; 4],
    trace_minus_start: u32,
    trace_minus_len: u32,
    trace_plus_start: u32,
    trace_plus_len: u32,
    /// `0 = Interior`, `1 = Pec`, `2 = Pmc`, `3 = SilverMuller`, `4 = ExcitedPort`, `5 = PassivePort`.
    face_kind: u32,
    has_sheet2d: u32,
    has_tfsf: u32,
    _pad0: u32,
}
impl_pod_scalar_struct!(GpuFaceQpMeta, 176, 320);

/// Per-probe interpolation range in `probe_pairs`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Pod, Zeroable)]
struct GpuProbeMeta {
    pair_start: u32,
    pair_len: u32,
    _pad: [u32; 2],
}

/// Per-nonlinear-cell scaled susceptibility & carrier parameters (`26 * R + 4 * u32`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuNlCellMeta<R> {
    /// Row-major `[[f64; 6]; 3]` nonlinear `d_glob` scaled by `e_ref`.
    d_glob: [R; 18],
    /// Isotropic `chi3_xxxx` scaled by `e_ref^2`.
    chi3_iso_xxxx: R,
    eps_min_diag: R,
    sig_tpa_coeff: R,
    sig_fca_coeff: R,
    fcd_coeff: R,
    tpa_gen_coeff: R,
    decay: R,
    one_minus_decay_tau: R,
    has_chi2: u32,
    has_chi3_iso: u32,
    rank4_offset: u32,
    has_carrier: u32,
}
impl_pod_scalar_struct!(GpuNlCellMeta, 120, 224);

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn non_empty_typed_bytes<T: Pod>(slice: &[T], fallback: &[T]) -> Vec<u8> {
    if slice.is_empty() {
        bytemuck::cast_slice(fallback).to_vec()
    } else {
        bytemuck::cast_slice(slice).to_vec()
    }
}

fn build_dgtd_wgsl(use_f64: bool) -> String {
    let rt = if use_f64 { "f64" } else { "f32" };
    DGTD_WGSL_TEMPLATE.replace("__R__", rt)
}

const DGTD_WGSL_TEMPLATE: &str = r#"
struct StepUniforms {
    dt: __R__,
    t_now: __R__,
    pulse_now: __R__,
    flux_alpha: __R__,
    rk_a: __R__,
    rk_b: __R__,
    rk_w_acc: __R__,
    rk_h_next: __R__,
    pulse_f0: __R__,
    pulse_fwidth: __R__,
    pulse_t0: __R__,
    pulse_amp: __R__,
    pulse_phase: __R__,
    symmetry_factor: __R__,
    inv_eps0: __R__,
    e_ref: __R__,
    inv_e_ref: __R__,
    _pad_r0: __R__,
    _pad_r1: __R__,
    _pad_r2: __R__,
    total_dofs: u32,
    total_qps: u32,
    n_cells: u32,
    n_freqs: u32,
    n_probes: u32,
    n_energy_wgs: u32,
    dim: u32,
    pol_mode: u32,
    rk_stage: u32,
    snap_idx: u32,
    pulse_kind: u32,
    has_nonlinear: u32,
    has_carriers: u32,
    _pad_u0: u32,
    _pad_u1: u32,
    _pad_u2: u32,
};

struct CellMeta {
    inv_eps_row0: vec4<__R__>,
    inv_eps_row1: vec4<__R__>,
    inv_eps_row2: vec4<__R__>,
    inv_mu_row0: vec4<__R__>,
    inv_mu_row1: vec4<__R__>,
    inv_mu_row2: vec4<__R__>,
    sigma_e_row0: vec4<__R__>,
    sigma_e_row1: vec4<__R__>,
    sigma_e_row2: vec4<__R__>,
    pml_damp: vec4<__R__>,
    eps_si_row0: vec4<__R__>,
    eps_si_row1: vec4<__R__>,
    eps_si_row2: vec4<__R__>,
    mu_si_row0: vec4<__R__>,
    mu_si_row1: vec4<__R__>,
    mu_si_row2: vec4<__R__>,
    dof_offset: u32,
    n_p: u32,
    mass_offset: u32,
    energy_mode: u32,
    nl_cell_slot: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

struct DofMeta {
    dipole_e_vec: vec4<__R__>,
    dipole_h_vec: vec4<__R__>,
    cell_id: u32,
    cell_dof_offset: u32,
    n_p: u32,
    dmat_row_offset: u32,
    lift_start: u32,
    lift_len: u32,
    _pad0: u32,
    _pad1: u32,
};

struct DofLift {
    lift_w: __R__,
    qp_idx: u32,
    _pad: u32,
};

struct TracePair {
    weight: __R__,
    global_dof: u32,
    _pad: u32,
};

struct FaceQpMeta {
    n_hat_ym: vec4<__R__>,
    yp_tfsf_delay: vec4<__R__>,
    port_e_mode: vec4<__R__>,
    port_h_mode: vec4<__R__>,
    tfsf_e_unit: vec4<__R__>,
    tfsf_h_unit: vec4<__R__>,
    sheet_sigma: vec4<__R__>,
    sheet_u_hat: vec4<__R__>,
    sheet_v_hat: vec4<__R__>,
    trace_minus_start: u32,
    trace_minus_len: u32,
    trace_plus_start: u32,
    trace_plus_len: u32,
    face_kind: u32,
    has_sheet2d: u32,
    has_tfsf: u32,
    _pad0: u32,
};

struct ProbeMeta {
    pair_start: u32,
    pair_len: u32,
    _pad0: u32,
    _pad1: u32,
};

struct NlCellMeta {
    d_glob: array<__R__, 18>,
    chi3_iso_xxxx: __R__,
    eps_min_diag: __R__,
    sig_tpa_coeff: __R__,
    sig_fca_coeff: __R__,
    fcd_coeff: __R__,
    tpa_gen_coeff: __R__,
    decay: __R__,
    one_minus_decay_tau: __R__,
    has_chi2: u32,
    has_chi3_iso: u32,
    rank4_offset: u32,
    has_carrier: u32,
};

@group(0) @binding(0) var<uniform> u: StepUniforms;
@group(0) @binding(1) var<storage, read_write> state_e: array<vec4<__R__>>;
@group(0) @binding(2) var<storage, read_write> state_h: array<vec4<__R__>>;
@group(0) @binding(3) var<storage, read_write> de_out: array<vec4<__R__>>;
@group(0) @binding(4) var<storage, read_write> dh_out: array<vec4<__R__>>;
@group(0) @binding(5) var<storage, read_write> aux0_e: array<vec4<__R__>>;
@group(0) @binding(6) var<storage, read_write> aux0_h: array<vec4<__R__>>;
@group(0) @binding(7) var<storage, read_write> aux1_e: array<vec4<__R__>>;
@group(0) @binding(8) var<storage, read_write> aux1_h: array<vec4<__R__>>;
@group(0) @binding(9) var<storage, read_write> qp_flux_e: array<vec4<__R__>>;
@group(0) @binding(10) var<storage, read_write> qp_flux_h: array<vec4<__R__>>;

@group(1) @binding(0) var<storage, read> cells_meta: array<CellMeta>;
@group(1) @binding(1) var<storage, read> dofs_meta: array<DofMeta>;
@group(1) @binding(2) var<storage, read> dmat_xyz: array<vec4<__R__>>;
@group(1) @binding(3) var<storage, read> dof_lifts: array<DofLift>;
@group(1) @binding(4) var<storage, read> face_qps: array<FaceQpMeta>;
@group(1) @binding(5) var<storage, read> trace_pairs: array<TracePair>;
@group(1) @binding(6) var<storage, read> cell_masses: array<__R__>;

@group(2) @binding(0) var<storage, read> dft_cos_sin: array<vec2<__R__>>;
@group(2) @binding(1) var<storage, read_write> dft_e_re: array<vec4<__R__>>;
@group(2) @binding(2) var<storage, read_write> dft_e_im: array<vec4<__R__>>;
@group(2) @binding(3) var<storage, read_write> dft_h_re: array<vec4<__R__>>;
@group(2) @binding(4) var<storage, read_write> dft_h_im: array<vec4<__R__>>;
@group(2) @binding(5) var<storage, read> probes_meta: array<ProbeMeta>;
@group(2) @binding(6) var<storage, read> probe_pairs: array<TracePair>;
@group(2) @binding(7) var<storage, read_write> probe_out_e: array<vec4<__R__>>;
@group(2) @binding(8) var<storage, read_write> probe_out_h: array<vec4<__R__>>;
@group(2) @binding(9) var<storage, read_write> energy_partials: array<__R__>;
@group(2) @binding(10) var<storage, read_write> energy_history: array<__R__>;

@group(3) @binding(0) var<storage, read> nl_cells_meta: array<NlCellMeta>;
@group(3) @binding(1) var<storage, read> nl_chi3_rank4: array<__R__>;
@group(3) @binding(2) var<storage, read_write> nl_carrier: array<__R__>;

const PI_R: __R__ = __R__(3.14159265358979323846);

fn eval_pulse_at(t: __R__) -> __R__ {
    let dt = t - u.pulse_t0;
    if (u.pulse_kind == 0u) {
        let sigma = __R__(1.0) / (__R__(2.0) * PI_R * max(u.pulse_fwidth, __R__(1.0)));
        let arg = dt / sigma;
        if (abs(arg) > __R__(14.0)) {
            return __R__(0.0);
        }
        let arg_f = f32(arg);
        let ph_f = f32(__R__(2.0) * PI_R * u.pulse_f0 * dt + u.pulse_phase);
        return u.pulse_amp * __R__(exp(-0.5 * arg_f * arg_f) * cos(ph_f));
    } else {
        let arg = PI_R * u.pulse_f0 * dt;
        let u2 = arg * arg;
        if (u2 > __R__(100.0)) {
            return __R__(0.0);
        }
        let u2_f = f32(u2);
        return u.pulse_amp * (__R__(1.0) - __R__(2.0) * u2) * __R__(exp(-u2_f));
    }
}

@compute @workgroup_size(64)
fn eval_face_qps(@builtin(global_invocation_id) gid: vec3<u32>) {
    let q_idx = gid.x;
    if (q_idx >= u.total_qps) {
        return;
    }
    let qp = face_qps[q_idx];
    let n = qp.n_hat_ym.xyz;
    let y_m = qp.n_hat_ym.w;
    let y_p = qp.yp_tfsf_delay.x;
    let y_sum = max(y_m + y_p, __R__(1e-30));
    let z_m = __R__(1.0) / max(y_m, __R__(1e-30));
    let z_p = __R__(1.0) / max(y_p, __R__(1e-30));
    let z_sum = max(z_m + z_p, __R__(1e-30));
    let alpha = u.flux_alpha;

    var e_m = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    var h_m = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    for (var i: u32 = 0u; i < qp.trace_minus_len; i = i + 1u) {
        let tp = trace_pairs[qp.trace_minus_start + i];
        e_m = e_m + tp.weight * state_e[tp.global_dof].xyz;
        h_m = h_m + tp.weight * state_h[tp.global_dof].xyz;
    }

    var e_p = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    var h_p = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    if (qp.face_kind == 0u) {
        for (var i: u32 = 0u; i < qp.trace_plus_len; i = i + 1u) {
            let tp = trace_pairs[qp.trace_plus_start + i];
            e_p = e_p + tp.weight * state_e[tp.global_dof].xyz;
            h_p = h_p + tp.weight * state_h[tp.global_dof].xyz;
        }
    } else if (qp.face_kind == 1u) {
        e_p = -e_m;
        h_p = h_m;
    } else if (qp.face_kind == 2u) {
        e_p = e_m;
        h_p = -h_m;
    } else if (qp.face_kind == 4u) {
        e_p = u.pulse_now * qp.port_e_mode.xyz;
        h_p = u.pulse_now * qp.port_h_mode.xyz;
    }

    if (qp.has_tfsf != 0u) {
        let s_inc = eval_pulse_at(u.t_now - qp.yp_tfsf_delay.y);
        e_p = e_p + s_inc * qp.tfsf_e_unit.xyz;
        h_p = h_p + s_inc * qp.tfsf_h_unit.xyz;
    }

    let de = e_p - e_m;
    let dh = h_p - h_m;
    let n_cross_dh = cross(n, dh);
    let n_cross_de = cross(n, de);
    let n_cross_n_cross_de = cross(n, n_cross_de);
    let n_cross_n_cross_dh = cross(n, n_cross_dh);

    var flux_e = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    var flux_h = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    if (qp.has_sheet2d != 0u) {
        let u_hat = qp.sheet_u_hat.xyz;
        let v_hat = qp.sheet_v_hat.xyz;
        let s00 = qp.sheet_sigma.x;
        let s01 = qp.sheet_sigma.y;
        let s10 = qp.sheet_sigma.z;
        let s11 = qp.sheet_sigma.w;
        let rhs_u = y_m * dot(e_m, u_hat) + y_p * dot(e_p, u_hat) + alpha * dot(n_cross_dh, u_hat);
        let rhs_v = y_m * dot(e_m, v_hat) + y_p * dot(e_p, v_hat) + alpha * dot(n_cross_dh, v_hat);
        let a00 = y_sum + s00;
        let a11 = y_sum + s11;
        let det_a = max(a00 * a11 - s01 * s10, __R__(1e-30));
        let eu_star = (a11 * rhs_u - s01 * rhs_v) / det_a;
        let ev_star = (-s10 * rhs_u + a00 * rhs_v) / det_a;
        let e_t_star = eu_star * u_hat + ev_star * v_hat;
        let ku = s00 * eu_star + s01 * ev_star;
        let kv = s10 * eu_star + s11 * ev_star;
        let k_vec = ku * u_hat + kv * v_hat;
        let c_h = y_p / y_sum;
        let c_e = alpha * (y_m * y_p) / y_sum;
        let c_k = y_m / y_sum;
        flux_e = c_h * n_cross_dh - c_e * n_cross_n_cross_de - c_k * k_vec;
        flux_h = -cross(n, e_t_star - e_m);
    } else {
        let c_dh = y_p / y_sum;
        let c_de = alpha * (y_m * y_p) / y_sum;
        flux_e = c_dh * n_cross_dh - c_de * n_cross_n_cross_de;
        let c_e = z_p / z_sum;
        let c_h = alpha * (z_m * z_p) / z_sum;
        flux_h = -c_e * n_cross_de - c_h * n_cross_n_cross_dh;
    }

    qp_flux_e[q_idx] = vec4<__R__>(flux_e, __R__(0.0));
    qp_flux_h[q_idx] = vec4<__R__>(flux_h, __R__(0.0));
}

fn invert_eps_mul(
    r0: vec3<__R__>,
    r1: vec3<__R__>,
    r2: vec3<__R__>,
    v: vec3<__R__>,
) -> vec4<__R__> {
    let c00 = r1.y * r2.z - r1.z * r2.y;
    let c01 = -(r1.x * r2.z - r1.z * r2.x);
    let c02 = r1.x * r2.y - r1.y * r2.x;
    let det = r0.x * c00 + r0.y * c01 + r0.z * c02;
    if (abs(det) <= __R__(1e-24)) {
        return vec4<__R__>(__R__(0.0), __R__(0.0), __R__(0.0), __R__(0.0));
    }
    let c10 = -(r0.y * r2.z - r0.z * r2.y);
    let c11 = r0.x * r2.z - r0.z * r2.x;
    let c12 = -(r0.x * r2.y - r0.y * r2.x);
    let c20 = r0.y * r1.z - r0.z * r1.y;
    let c21 = -(r0.x * r1.z - r0.z * r1.x);
    let c22 = r0.x * r1.y - r0.y * r1.x;
    let inv_det = __R__(1.0) / det;
    let x = (c00 * v.x + c10 * v.y + c20 * v.z) * inv_det;
    let y = (c01 * v.x + c11 * v.y + c21 * v.z) * inv_det;
    let z = (c02 * v.x + c12 * v.y + c22 * v.z) * inv_det;
    return vec4<__R__>(x, y, z, __R__(1.0));
}

@compute @workgroup_size(64)
fn eval_dof_rhs(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs) {
        return;
    }
    let dm = dofs_meta[d];
    var dxh_y: __R__ = __R__(0.0);
    var dxh_z: __R__ = __R__(0.0);
    var dyh_x: __R__ = __R__(0.0);
    var dyh_z: __R__ = __R__(0.0);
    var dzh_x: __R__ = __R__(0.0);
    var dzh_y: __R__ = __R__(0.0);
    var dxe_y: __R__ = __R__(0.0);
    var dxe_z: __R__ = __R__(0.0);
    var dye_x: __R__ = __R__(0.0);
    var dye_z: __R__ = __R__(0.0);
    var dze_x: __R__ = __R__(0.0);
    var dze_y: __R__ = __R__(0.0);

    for (var b: u32 = 0u; b < dm.n_p; b = b + 1u) {
        let dxyz = dmat_xyz[dm.dmat_row_offset + b];
        let eb = state_e[dm.cell_dof_offset + b].xyz;
        let hb = state_h[dm.cell_dof_offset + b].xyz;
        dxh_y = dxh_y + dxyz.x * hb.y;
        dxh_z = dxh_z + dxyz.x * hb.z;
        dyh_x = dyh_x + dxyz.y * hb.x;
        dyh_z = dyh_z + dxyz.y * hb.z;
        dzh_x = dzh_x + dxyz.z * hb.x;
        dzh_y = dzh_y + dxyz.z * hb.y;

        dxe_y = dxe_y + dxyz.x * eb.y;
        dxe_z = dxe_z + dxyz.x * eb.z;
        dye_x = dye_x + dxyz.y * eb.x;
        dye_z = dye_z + dxyz.y * eb.z;
        dze_x = dze_x + dxyz.z * eb.x;
        dze_y = dze_y + dxyz.z * eb.y;
    }

    var rhs_d = vec3<__R__>(dyh_z - dzh_y, dzh_x - dxh_z, dxh_y - dyh_x);
    var rhs_b = vec3<__R__>(-(dye_z - dze_y), -(dze_x - dxe_z), -(dxe_y - dye_x));

    for (var l: u32 = 0u; l < dm.lift_len; l = l + 1u) {
        let item = dof_lifts[dm.lift_start + l];
        rhs_d = rhs_d + item.lift_w * qp_flux_e[item.qp_idx].xyz;
        rhs_b = rhs_b + item.lift_w * qp_flux_h[item.qp_idx].xyz;
    }

    rhs_d = rhs_d - u.pulse_now * dm.dipole_e_vec.xyz;
    rhs_b = rhs_b - u.pulse_now * dm.dipole_h_vec.xyz;

    let ea = state_e[d].xyz;
    let ha = state_h[d].xyz;
    let cell = cells_meta[dm.cell_id];
    if (cell.pml_damp.w > __R__(0.5)) {
        let j_ohm = vec3<__R__>(
            dot(cell.sigma_e_row0.xyz, ea),
            dot(cell.sigma_e_row1.xyz, ea),
            dot(cell.sigma_e_row2.xyz, ea)
        );
        rhs_d = rhs_d - j_ohm;
    }

    var de = vec3<__R__>(
        dot(cell.inv_eps_row0.xyz, rhs_d),
        dot(cell.inv_eps_row1.xyz, rhs_d),
        dot(cell.inv_eps_row2.xyz, rhs_d)
    );

    if (u.has_nonlinear != 0u && cell.nl_cell_slot != 0xffffffffu) {
        let nl = nl_cells_meta[cell.nl_cell_slot];
        let en = ea * u.inv_e_ref;
        let en2 = dot(en, en);
        let nc_tilde = max(nl_carrier[d], __R__(0.0));

        let sig_nl = nl.sig_tpa_coeff * en2 + nl.sig_fca_coeff * nc_tilde;
        if (sig_nl > __R__(0.0)) {
            rhs_d = rhs_d - sig_nl * ea;
        }

        var d0 = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
        var d1 = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
        var d2 = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));

        if (nl.has_chi2 != 0u) {
            let s4 = __R__(4.0) * dm.dipole_e_vec.w;
            d0 = vec3<__R__>(
                s4 * (nl.d_glob[0] * en.x + nl.d_glob[5] * en.y + nl.d_glob[4] * en.z),
                s4 * (nl.d_glob[5] * en.x + nl.d_glob[1] * en.y + nl.d_glob[3] * en.z),
                s4 * (nl.d_glob[4] * en.x + nl.d_glob[3] * en.y + nl.d_glob[2] * en.z)
            );
            d1 = vec3<__R__>(
                s4 * (nl.d_glob[6] * en.x + nl.d_glob[11] * en.y + nl.d_glob[10] * en.z),
                s4 * (nl.d_glob[11] * en.x + nl.d_glob[7] * en.y + nl.d_glob[9] * en.z),
                s4 * (nl.d_glob[10] * en.x + nl.d_glob[9] * en.y + nl.d_glob[8] * en.z)
            );
            d2 = vec3<__R__>(
                s4 * (nl.d_glob[12] * en.x + nl.d_glob[17] * en.y + nl.d_glob[16] * en.z),
                s4 * (nl.d_glob[17] * en.x + nl.d_glob[13] * en.y + nl.d_glob[15] * en.z),
                s4 * (nl.d_glob[16] * en.x + nl.d_glob[15] * en.y + nl.d_glob[14] * en.z)
            );
        }

        if (nl.rank4_offset != 0xffffffffu) {
            let r_base = nl.rank4_offset;
            var k3 = array<__R__, 9>(
                __R__(0.0), __R__(0.0), __R__(0.0),
                __R__(0.0), __R__(0.0), __R__(0.0),
                __R__(0.0), __R__(0.0), __R__(0.0)
            );
            for (var r: u32 = 0u; r < 3u; r = r + 1u) {
                for (var l: u32 = 0u; l < 3u; l = l + 1u) {
                    var acc: __R__ = __R__(0.0);
                    for (var j: u32 = 0u; j < 3u; j = j + 1u) {
                        let ej = select(select(en.x, en.y, j == 1u), en.z, j == 2u);
                        for (var k: u32 = 0u; k < 3u; k = k + 1u) {
                            let ek = select(select(en.x, en.y, k == 1u), en.z, k == 2u);
                            let t_a = nl_chi3_rank4[r_base + r * 27u + j * 9u + k * 3u + l];
                            let t_b = nl_chi3_rank4[r_base + r * 27u + j * 9u + l * 3u + k];
                            let t_c = nl_chi3_rank4[r_base + r * 27u + l * 9u + j * 3u + k];
                            acc = acc + (t_a + t_b + t_c) * ej * ek;
                        }
                    }
                    k3[r * 3u + l] = acc;
                }
            }
            d0 = d0 + vec3<__R__>(k3[0], k3[1], k3[2]);
            d1 = d1 + vec3<__R__>(k3[3], k3[4], k3[5]);
            d2 = d2 + vec3<__R__>(k3[6], k3[7], k3[8]);
        } else if (nl.has_chi3_iso != 0u) {
            let c3 = nl.chi3_iso_xxxx;
            let c3_diag = c3 * en2;
            let c3_2 = __R__(2.0) * c3;
            d0 = d0 + vec3<__R__>(c3_diag + c3_2 * en.x * en.x, c3_2 * en.x * en.y, c3_2 * en.x * en.z);
            d1 = d1 + vec3<__R__>(c3_2 * en.y * en.x, c3_diag + c3_2 * en.y * en.y, c3_2 * en.y * en.z);
            d2 = d2 + vec3<__R__>(c3_2 * en.z * en.x, c3_2 * en.z * en.y, c3_diag + c3_2 * en.z * en.z);
        }

        if (nl.fcd_coeff != __R__(0.0) && nc_tilde > __R__(0.0)) {
            let dfcd = nl.fcd_coeff * nc_tilde;
            d0.x = d0.x + dfcd;
            d1.y = d1.y + dfcd;
            d2.z = d2.z + dfcd;
        }

        let rsum0 = abs(d0.x) + abs(d0.y) + abs(d0.z);
        let rsum1 = abs(d1.x) + abs(d1.y) + abs(d1.z);
        let rsum2 = abs(d2.x) + abs(d2.y) + abs(d2.z);
        let max_row = max(max(rsum0, rsum1), rsum2);
        let cap = __R__(0.45) * nl.eps_min_diag;
        if (max_row > cap && max_row > __R__(0.0)) {
            let sc = cap / max_row;
            d0 = d0 * sc;
            d1 = d1 * sc;
            d2 = d2 * sc;
        }

        let eps_r0 = cell.eps_si_row0.xyz + d0;
        let eps_r1 = cell.eps_si_row1.xyz + d1;
        let eps_r2 = cell.eps_si_row2.xyz + d2;
        let sol = invert_eps_mul(eps_r0, eps_r1, eps_r2, rhs_d);
        if (sol.w > __R__(0.5)) {
            de = sol.xyz * u.inv_eps0;
        } else {
            de = vec3<__R__>(
                dot(cell.inv_eps_row0.xyz, rhs_d),
                dot(cell.inv_eps_row1.xyz, rhs_d),
                dot(cell.inv_eps_row2.xyz, rhs_d)
            );
        }
    }

    var dh = vec3<__R__>(
        dot(cell.inv_mu_row0.xyz, rhs_b),
        dot(cell.inv_mu_row1.xyz, rhs_b),
        dot(cell.inv_mu_row2.xyz, rhs_b)
    );

    let damp = cell.pml_damp.xyz;
    de = de - damp * ea;
    dh = dh - damp * ha;

    if (u.dim == 2u) {
        if (u.pol_mode == 1u) {
            de.x = __R__(0.0);
            de.y = __R__(0.0);
            dh.z = __R__(0.0);
        } else if (u.pol_mode == 2u) {
            de.z = __R__(0.0);
            dh.x = __R__(0.0);
            dh.y = __R__(0.0);
        }
    }

    de_out[d] = vec4<__R__>(de, __R__(0.0));
    dh_out[d] = vec4<__R__>(dh, __R__(0.0));
}

@compute @workgroup_size(64)
fn lserk4_stage(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs) {
        return;
    }
    let re = u.rk_a * aux0_e[d].xyz + u.dt * de_out[d].xyz;
    let rh = u.rk_a * aux0_h[d].xyz + u.dt * dh_out[d].xyz;
    aux0_e[d] = vec4<__R__>(re, __R__(0.0));
    aux0_h[d] = vec4<__R__>(rh, __R__(0.0));
    state_e[d] = vec4<__R__>(state_e[d].xyz + u.rk_b * re, __R__(0.0));
    state_h[d] = vec4<__R__>(state_h[d].xyz + u.rk_b * rh, __R__(0.0));
}

@compute @workgroup_size(64)
fn rk4_save_init(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs) {
        return;
    }
    let se = state_e[d];
    let sh = state_h[d];
    aux0_e[d] = se;
    aux0_h[d] = sh;
    aux1_e[d] = se;
    aux1_h[d] = sh;
}

@compute @workgroup_size(64)
fn rk4_stage(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs) {
        return;
    }
    let de = de_out[d].xyz;
    let dh = dh_out[d].xyz;
    let ae = aux1_e[d].xyz + u.rk_w_acc * de;
    let ah = aux1_h[d].xyz + u.rk_w_acc * dh;
    aux1_e[d] = vec4<__R__>(ae, __R__(0.0));
    aux1_h[d] = vec4<__R__>(ah, __R__(0.0));
    if (u.rk_stage < 3u) {
        state_e[d] = vec4<__R__>(aux0_e[d].xyz + u.rk_h_next * de, __R__(0.0));
        state_h[d] = vec4<__R__>(aux0_h[d].xyz + u.rk_h_next * dh, __R__(0.0));
    } else {
        state_e[d] = vec4<__R__>(ae, __R__(0.0));
        state_h[d] = vec4<__R__>(ah, __R__(0.0));
    }
}

@compute @workgroup_size(64)
fn step_carriers(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs || u.has_carriers == 0u) {
        return;
    }
    let cell = cells_meta[dofs_meta[d].cell_id];
    if (cell.nl_cell_slot == 0xffffffffu) {
        return;
    }
    let nl = nl_cells_meta[cell.nl_cell_slot];
    if (nl.has_carrier == 0u) {
        return;
    }
    let ea = state_e[d].xyz;
    let en = ea * u.inv_e_ref;
    let en2 = dot(en, en);
    let nc_old = nl_carrier[d];
    let nc_next = max(nc_old * nl.decay + nl.tpa_gen_coeff * en2 * en2 * nl.one_minus_decay_tau, __R__(0.0));
    nl_carrier[d] = nc_next;
}

@compute @workgroup_size(64)
fn accum_dof_dft(@builtin(global_invocation_id) gid: vec3<u32>) {
    let d = gid.x;
    if (d >= u.total_dofs) {
        return;
    }
    let ea = state_e[d].xyz;
    let ha = state_h[d].xyz;
    for (var fi: u32 = 0u; fi < u.n_freqs; fi = fi + 1u) {
        let cs = dft_cos_sin[fi];
        let out_idx = fi * u.total_dofs + d;
        dft_e_re[out_idx] = vec4<__R__>(dft_e_re[out_idx].xyz + cs.x * ea, __R__(0.0));
        dft_e_im[out_idx] = vec4<__R__>(dft_e_im[out_idx].xyz + cs.y * ea, __R__(0.0));
        dft_h_re[out_idx] = vec4<__R__>(dft_h_re[out_idx].xyz + cs.x * ha, __R__(0.0));
        dft_h_im[out_idx] = vec4<__R__>(dft_h_im[out_idx].xyz + cs.y * ha, __R__(0.0));
    }
}

@compute @workgroup_size(64)
fn sample_probes(@builtin(global_invocation_id) gid: vec3<u32>) {
    let pi = gid.x;
    if (pi >= u.n_probes) {
        return;
    }
    let pm = probes_meta[pi];
    var ep = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    var hp = vec3<__R__>(__R__(0.0), __R__(0.0), __R__(0.0));
    for (var i: u32 = 0u; i < pm.pair_len; i = i + 1u) {
        let tp = probe_pairs[pm.pair_start + i];
        ep = ep + tp.weight * state_e[tp.global_dof].xyz;
        hp = hp + tp.weight * state_h[tp.global_dof].xyz;
    }
    let out_idx = u.snap_idx * u.n_probes + pi;
    probe_out_e[out_idx] = vec4<__R__>(ep, __R__(0.0));
    probe_out_h[out_idx] = vec4<__R__>(hp, __R__(0.0));
}

var<workgroup> wg_energy: array<__R__, 64>;

@compute @workgroup_size(64)
fn dgtd_energy_pass1(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
) {
    let ci = gid.x;
    var val: __R__ = __R__(0.0);
    if (ci < u.n_cells) {
        let cell = cells_meta[ci];
        if (cell.energy_mode == 1u) {
            var w_cell: __R__ = __R__(0.0);
            for (var a: u32 = 0u; a < cell.n_p; a = a + 1u) {
                let ea = state_e[cell.dof_offset + a].xyz;
                let ha = state_h[cell.dof_offset + a].xyz;
                let d_a = vec3<__R__>(
                    dot(cell.eps_si_row0.xyz, ea),
                    dot(cell.eps_si_row1.xyz, ea),
                    dot(cell.eps_si_row2.xyz, ea)
                );
                let b_a = vec3<__R__>(
                    dot(cell.mu_si_row0.xyz, ha),
                    dot(cell.mu_si_row1.xyz, ha),
                    dot(cell.mu_si_row2.xyz, ha)
                );
                let ma = cell_masses[cell.mass_offset + a];
                w_cell = w_cell + ma * (dot(ea, d_a) + dot(ha, b_a));
            }
            val = __R__(0.5) * w_cell;
        } else if (cell.energy_mode == 2u) {
            var w_cell: __R__ = __R__(0.0);
            for (var a: u32 = 0u; a < cell.n_p; a = a + 1u) {
                let ea = state_e[cell.dof_offset + a].xyz;
                let ha = state_h[cell.dof_offset + a].xyz;
                let d_a = vec3<__R__>(
                    dot(cell.eps_si_row0.xyz, ea),
                    dot(cell.eps_si_row1.xyz, ea),
                    dot(cell.eps_si_row2.xyz, ea)
                );
                let b_a = vec3<__R__>(
                    dot(cell.mu_si_row0.xyz, ha),
                    dot(cell.mu_si_row1.xyz, ha),
                    dot(cell.mu_si_row2.xyz, ha)
                );
                let row_off = cell.mass_offset + a * cell.n_p;
                for (var b: u32 = 0u; b < cell.n_p; b = b + 1u) {
                    let eb = state_e[cell.dof_offset + b].xyz;
                    let hb = state_h[cell.dof_offset + b].xyz;
                    let mab = cell_masses[row_off + b];
                    w_cell = w_cell + mab * (dot(eb, d_a) + dot(hb, b_a));
                }
            }
            val = __R__(0.5) * w_cell;
        }
    }
    wg_energy[lid.x] = val;
    workgroupBarrier();

    var stride: u32 = 32u;
    loop {
        if (stride == 0u) {
            break;
        }
        if (lid.x < stride) {
            wg_energy[lid.x] = wg_energy[lid.x] + wg_energy[lid.x + stride];
        }
        workgroupBarrier();
        stride = stride >> 1u;
    }
    if (lid.x == 0u) {
        energy_partials[wid.x] = wg_energy[0];
    }
}

var<workgroup> wg_sum2: array<__R__, 256>;

@compute @workgroup_size(256)
fn dgtd_energy_pass2(@builtin(local_invocation_id) lid: vec3<u32>) {
    var acc: __R__ = __R__(0.0);
    var i: u32 = lid.x;
    loop {
        if (i >= u.n_energy_wgs) {
            break;
        }
        acc = acc + energy_partials[i];
        i = i + 256u;
    }
    wg_sum2[lid.x] = acc;
    workgroupBarrier();

    var stride: u32 = 128u;
    loop {
        if (stride == 0u) {
            break;
        }
        if (lid.x < stride) {
            wg_sum2[lid.x] = wg_sum2[lid.x] + wg_sum2[lid.x + stride];
        }
        workgroupBarrier();
        stride = stride >> 1u;
    }
    if (lid.x == 0u) {
        energy_history[u.snap_idx] = wg_sum2[0] * u.symmetry_factor;
    }
}
"#;

/// Carpenter–Kennedy 1994 5-stage 4th-order Low-Storage Explicit Runge–Kutta (`LSERK4`) coefficients.
const LSERK4_A: [f64; 5] = [
    0.0,
    -567301805773.0 / 1357537059087.0,
    -2404267990393.0 / 2016746695238.0,
    -3550918686646.0 / 2091501179385.0,
    -1275806237668.0 / 842570457699.0,
];
const LSERK4_B: [f64; 5] = [
    1432997174477.0 / 9575080441755.0,
    5161836677717.0 / 13612068292357.0,
    1720146321549.0 / 2090206949498.0,
    3134564353537.0 / 4481467310338.0,
    2277821191437.0 / 14882151754819.0,
];
const LSERK4_C: [f64; 5] = [
    0.0,
    1432997174477.0 / 9575080441755.0,
    2526269341429.0 / 6820363962896.0,
    2006345519317.0 / 3224310063776.0,
    2802321613138.0 / 2924317926251.0,
];

fn read_back_slice<T: Pod>(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    src: &wgpu::Buffer,
    count: usize,
) -> Result<Vec<T>, String> {
    if count == 0 {
        return Ok(Vec::new());
    }
    let byte_len = (count * std::mem::size_of::<T>()) as u64;
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("nyx-ems-dgtd-readback"),
        size: byte_len,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder =
        device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("nyx-ems-dgtd-readback-enc") });
    encoder.copy_buffer_to_buffer(src, 0, &staging, 0, byte_len);
    queue.submit(Some(encoder.finish()));

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device
        .poll(wgpu::PollType::wait_indefinitely())
        .map_err(|e| format!("wgpu device poll error during DGTD readback: {e}"))?;
    rx.recv()
        .map_err(|e| format!("wgpu map_async channel recv error: {e}"))?
        .map_err(|e| format!("wgpu buffer map_async error: {e}"))?;

    let mapped = slice.get_mapped_range().map_err(|e| format!("wgpu get_mapped_range error: {e}"))?;
    let out: Vec<T> = bytemuck::cast_slice(&mapped).to_vec();
    drop(mapped);
    staging.unmap();
    Ok(out)
}

/// Attempts to execute `Study::Dgtd` on a `wgpu` compute adapter.
///
/// Returns:
/// - `Ok(Ok(DgtdSolve))` when the simulation was executed on the GPU.
/// - `Ok(Err(fallback_reason))` when the configuration requires CPU features (such as multi-rate
///   Local Time-Stepping `n_lts_tiers > 1`, dispersive ADE poles, or explicit `precision = "f64"`
///   on adapters without `SHADER_F64`) or no suitable GPU adapter is present, allowing
///   [`crate::dgtd::solve_dgtd`] to fall back cleanly to CPU with a `DiagnosticWarning`
///   (`code: "gpu_fallback_to_cpu"`).
/// - `Err(fatal)` only on an unrecoverable internal error.
pub fn try_solve_dgtd_gpu(asm: &DgtdAssembly, require_discrete: bool) -> Result<Result<DgtdSolve, String>, String> {
    if asm.config.time_stepper == DgtdTimeStepper::Lts && asm.n_lts_tiers > 1 {
        return Ok(Err(
            "DGTD multi-rate Local Time-Stepping (n_lts_tiers > 1) requires hierarchical CPU sub-stepping; falling back to CPU DGTD"
                .into(),
        ));
    }
    if asm.total_pole_states > 0 {
        return Ok(Err("DGTD dispersive Drude/Lorentz/gyrotropic ADE poles active; falling back to CPU DGTD".into()));
    }
    if asm.total_dofs == 0 || asm.cells.is_empty() {
        return Ok(Err("DGTD assembly has zero degrees of freedom; falling back to CPU DGTD".into()));
    }

    let ctx = match EmsGpuContext::acquire_default() {
        Ok(c) => c,
        Err(e) => return Ok(Err(format!("wgpu adapter unavailable ({e}); falling back to CPU DGTD"))),
    };
    if require_discrete && !ctx.is_discrete() {
        return Ok(Err(format!(
            "auto backend requires discrete GPU, found {:?} '{}'; falling back to CPU DGTD",
            ctx.adapter_info.device_type,
            ctx.adapter_description()
        )));
    }

    let has_nonlinear = asm.cells.iter().any(|c| c.nonlinear.is_some()) || asm.has_carrier_dynamics;
    let use_f64 = match select_td_precision(&ctx, asm.config.precision, has_nonlinear) {
        Ok(v) => v,
        Err(reason) => return Ok(Err(reason)),
    };

    if use_f64 {
        solve_dgtd_on_gpu_typed::<f64>(&ctx, asm, true)
    } else {
        solve_dgtd_on_gpu_typed::<f32>(&ctx, asm, false)
    }
}

fn solve_dgtd_on_gpu_typed<R: GpuScalar>(
    ctx: &EmsGpuContext,
    asm: &DgtdAssembly,
    use_f64: bool,
) -> Result<Result<DgtdSolve, String>, String>
where
    DgtdStepUniforms<R>: Pod,
    GpuCellMeta<R>: Pod,
    GpuDofMeta<R>: Pod,
    GpuDofLift<R>: Pod,
    GpuTracePair<R>: Pod,
    GpuFaceQpMeta<R>: Pod,
    GpuNlCellMeta<R>: Pod,
{
    let t_start = Instant::now();
    let nc = asm.cells.len();
    let ndofs = asm.total_dofs;
    let n_freqs = asm.config.freqs_hz.len();
    let n_ports = asm.ports.len();
    let n_flux = asm.flux_monitors.len();
    let n_probes = asm.probes.len();
    let dt = asm.dt_0;
    let has_nonlinear = asm.cells.iter().any(|c| c.nonlinear.is_some()) || asm.has_carrier_dynamics;
    let has_carriers = asm.has_carrier_dynamics;

    // 1. Estimate peak electric field scale for f32 nonlinear normalization
    let pulse_amp_abs = match asm.config.pulse {
        DgtdPulse::GaussianModulated { amplitude, .. } | DgtdPulse::Ricker { amplitude, .. } => amplitude.abs(),
    };
    let mut est_e_peak = pulse_amp_abs.max(if asm.config.initial_cavity_mode.is_some() { 1.0 } else { 0.0 });
    for cell in &asm.cells {
        for face in &cell.faces {
            for qp in &face.qps {
                if let Some((em, _)) = qp.port_mode {
                    let em_n = (em[0] * em[0] + em[1] * em[1] + em[2] * em[2]).sqrt();
                    est_e_peak = est_e_peak.max(pulse_amp_abs * em_n);
                }
            }
        }
    }
    let (e_ref, nc_ref) = if use_f64 { (1.0_f64, 1.0_f64) } else { (est_e_peak.clamp(1.0, 1.0e10), DGTD_NC_REF) };

    // Normalize cell quadrature masses by mean non-PML cell measure and omit EPS0 inside energy shader
    // so nanoscale DGTD elements avoid f32 subnormal underflow while recovering exact Joules in host f64.
    let non_pml_count = asm.cells.iter().filter(|c| !c.in_pml).count();
    let mean_measure = if non_pml_count > 0 {
        let sum_m: f64 = asm.cells.iter().filter(|c| !c.in_pml).map(|c| c.measure).sum();
        (sum_m / (non_pml_count as f64)).max(1e-30)
    } else {
        1.0
    };
    let inv_mean_measure = 1.0 / mean_measure;
    let mu0_over_eps0 = MU0 / EPS0;

    // 2. Initial nodal state (E, H)
    let mut init_e = vec![[R::default(); 4]; ndofs];
    let init_h = vec![[R::default(); 4]; ndofs];
    if let Some([mx, my, mz]) = asm.config.initial_cavity_mode {
        let (lo, hi) = asm.mixed.bbox();
        let lx = (hi[0] - lo[0]).max(1e-12);
        let ly = (hi[1] - lo[1]).max(1e-12);
        let lz = (hi[2] - lo[2]).max(1e-12);
        let kx = (mx as f64) * PI / lx;
        let ky = (my as f64) * PI / ly;
        let kz = if asm.dim == 3 { (mz as f64) * PI / lz } else { 0.0 };
        for cell in &asm.cells {
            for (a, p) in cell.node_coords.iter().enumerate() {
                let x = p[0] - lo[0];
                let y = p[1] - lo[1];
                let z = p[2] - lo[2];
                let e_val = if asm.dim == 2 {
                    match asm.config.polarization {
                        DgtdPolarization::Te => {
                            [(kx * x).cos() * (ky * y).sin(), -(kx * x).sin() * (ky * y).cos(), 0.0]
                        }
                        DgtdPolarization::Tm | DgtdPolarization::Full => [0.0, 0.0, (kx * x).sin() * (ky * y).sin()],
                    }
                } else {
                    [0.0, 0.0, (kx * x).sin() * (ky * y).sin() * (kz * z).cos()]
                };
                init_e[cell.dof_offset + a] =
                    [R::from_f64(e_val[0]), R::from_f64(e_val[1]), R::from_f64(e_val[2]), R::default()];
            }
        }
    }

    // 3. Flatten cell metadata, nonlinear metadata, differentiation matrices, face QPs, trace pairs, and lifts
    let mut cells_meta = Vec::with_capacity(nc);
    let mut cell_masses = Vec::new();
    let mut dofs_meta = vec![GpuDofMeta::<R>::default(); ndofs];
    let mut dmat_xyz = Vec::new();
    let mut face_qps = Vec::new();
    let mut trace_pairs = Vec::new();
    let mut per_dof_lifts: Vec<Vec<GpuDofLift<R>>> = vec![Vec::new(); ndofs];
    let mut nl_cells_meta: Vec<GpuNlCellMeta<R>> = Vec::new();
    let mut nl_rank4_flat: Vec<R> = Vec::new();

    // Pre-accumulate dipole moment lift vectors per global DOF
    let mut dip_e_per_dof = vec![[0.0f64; 3]; ndofs];
    let mut dip_h_per_dof = vec![[0.0f64; 3]; ndofs];
    for dip in &asm.dipoles {
        let base = asm.cells[dip.cell_index].dof_offset;
        for &(a, lw) in &dip.lift_weights {
            let v = [lw * dip.moment[0], lw * dip.moment[1], lw * dip.moment[2]];
            if dip.magnetic {
                for k in 0..3 {
                    dip_h_per_dof[base + a][k] += v[k];
                }
            } else {
                for k in 0..3 {
                    dip_e_per_dof[base + a][k] += v[k];
                }
            }
        }
    }

    let excited_port_idx = asm.config.excited_port.saturating_sub(1);
    for (ci, cell) in asm.cells.iter().enumerate() {
        let n_p = cell.n_p;
        let base = cell.dof_offset;
        let mass_offset = cell_masses.len() as u32;
        for &m in &cell.mass {
            cell_masses.push(R::from_f64(m * inv_mean_measure));
        }
        let energy_mode = if cell.in_pml {
            0u32
        } else if cell.is_diagonal_mass {
            1u32
        } else {
            2u32
        };
        let pml_damp = if cell.in_pml {
            let s_sum = cell.pml_sigma[0] + cell.pml_sigma[1] + cell.pml_sigma[2];
            [
                R::from_f64(s_sum - cell.pml_sigma[0]),
                R::from_f64(s_sum - cell.pml_sigma[1]),
                R::from_f64(s_sum - cell.pml_sigma[2]),
                R::from_f64(if cell.has_sigma_e { 1.0 } else { 0.0 }),
            ]
        } else {
            [R::default(), R::default(), R::default(), R::from_f64(if cell.has_sigma_e { 1.0 } else { 0.0 })]
        };

        let nl_cell_slot = if let Some(ref nl) = cell.nonlinear {
            let slot = nl_cells_meta.len() as u32;
            let mut d_glob = [R::default(); 18];
            let has_chi2 = if let Some(ref dg) = nl.chi2_d_global {
                for r in 0..3 {
                    for v in 0..6 {
                        d_glob[r * 6 + v] = R::from_f64(dg[r][v] * e_ref);
                    }
                }
                1u32
            } else {
                0u32
            };
            let mut has_chi3_iso = 0u32;
            let mut rank4_offset = u32::MAX;
            if let Some(ref r4) = nl.chi3_rank4_global {
                rank4_offset = nl_rank4_flat.len() as u32;
                let e2 = e_ref * e_ref;
                for i in 0..3 {
                    for j in 0..3 {
                        for k in 0..3 {
                            for l in 0..3 {
                                nl_rank4_flat.push(R::from_f64(r4[i][j][k][l] * e2));
                            }
                        }
                    }
                }
            } else if nl.chi3_iso_xxxx != 0.0 {
                has_chi3_iso = 1u32;
            }
            let e2 = e_ref * e_ref;
            let sig_tpa_coeff =
                R::from_f64((2.0 / 3.0) * EPS0 * EPS0 * C0_SI * C0_SI * nl.n0 * nl.n0 * nl.beta_tpa * e2);
            let sig_fca_coeff = R::from_f64(EPS0 * C0_SI * nl.n0 * nl.sigma_fca * nc_ref);
            let fcd_coeff = R::from_f64(2.0 * nl.n0 * nl.k_fcd * nc_ref);
            let cell_has_carrier = u32::from(nl.beta_tpa > 0.0 && (nl.sigma_fca > 0.0 || nl.k_fcd != 0.0));
            let tpa_gen_coeff = if cell_has_carrier != 0 && nl.omega_0 > 0.0 {
                let pref = TPA_RMS_SCALE * C0_SI * EPS0 * nl.n0 * e2;
                R::from_f64(((nl.beta_tpa * pref * pref) / (2.0 * HBAR * nl.omega_0)) / nc_ref)
            } else {
                R::default()
            };
            let (decay, one_minus_decay_tau) = if nl.tau_c > 0.0 && nl.tau_c.is_finite() {
                let x = dt / nl.tau_c;
                if x > 700.0 {
                    (R::default(), R::from_f64(nl.tau_c))
                } else {
                    (R::from_f64((-x).exp()), R::from_f64(nl.tau_c * (-(-x).exp_m1())))
                }
            } else {
                (R::from_f64(1.0), R::from_f64(dt))
            };
            nl_cells_meta.push(GpuNlCellMeta {
                d_glob,
                chi3_iso_xxxx: R::from_f64(nl.chi3_iso_xxxx * e2),
                eps_min_diag: R::from_f64(nl.eps_min_diag),
                sig_tpa_coeff,
                sig_fca_coeff,
                fcd_coeff,
                tpa_gen_coeff,
                decay,
                one_minus_decay_tau,
                has_chi2,
                has_chi3_iso,
                rank4_offset,
                has_carrier: cell_has_carrier,
            });
            slot
        } else {
            u32::MAX
        };

        cells_meta.push(GpuCellMeta {
            inv_eps_row0: [
                R::from_f64(cell.inv_eps[0][0]),
                R::from_f64(cell.inv_eps[0][1]),
                R::from_f64(cell.inv_eps[0][2]),
                R::default(),
            ],
            inv_eps_row1: [
                R::from_f64(cell.inv_eps[1][0]),
                R::from_f64(cell.inv_eps[1][1]),
                R::from_f64(cell.inv_eps[1][2]),
                R::default(),
            ],
            inv_eps_row2: [
                R::from_f64(cell.inv_eps[2][0]),
                R::from_f64(cell.inv_eps[2][1]),
                R::from_f64(cell.inv_eps[2][2]),
                R::default(),
            ],
            inv_mu_row0: [
                R::from_f64(cell.inv_mu[0][0]),
                R::from_f64(cell.inv_mu[0][1]),
                R::from_f64(cell.inv_mu[0][2]),
                R::default(),
            ],
            inv_mu_row1: [
                R::from_f64(cell.inv_mu[1][0]),
                R::from_f64(cell.inv_mu[1][1]),
                R::from_f64(cell.inv_mu[1][2]),
                R::default(),
            ],
            inv_mu_row2: [
                R::from_f64(cell.inv_mu[2][0]),
                R::from_f64(cell.inv_mu[2][1]),
                R::from_f64(cell.inv_mu[2][2]),
                R::default(),
            ],
            sigma_e_row0: [
                R::from_f64(cell.sigma_e[0][0]),
                R::from_f64(cell.sigma_e[0][1]),
                R::from_f64(cell.sigma_e[0][2]),
                R::default(),
            ],
            sigma_e_row1: [
                R::from_f64(cell.sigma_e[1][0]),
                R::from_f64(cell.sigma_e[1][1]),
                R::from_f64(cell.sigma_e[1][2]),
                R::default(),
            ],
            sigma_e_row2: [
                R::from_f64(cell.sigma_e[2][0]),
                R::from_f64(cell.sigma_e[2][1]),
                R::from_f64(cell.sigma_e[2][2]),
                R::default(),
            ],
            pml_damp,
            eps_si_row0: [
                R::from_f64(cell.eps_r[0][0]),
                R::from_f64(cell.eps_r[0][1]),
                R::from_f64(cell.eps_r[0][2]),
                R::default(),
            ],
            eps_si_row1: [
                R::from_f64(cell.eps_r[1][0]),
                R::from_f64(cell.eps_r[1][1]),
                R::from_f64(cell.eps_r[1][2]),
                R::default(),
            ],
            eps_si_row2: [
                R::from_f64(cell.eps_r[2][0]),
                R::from_f64(cell.eps_r[2][1]),
                R::from_f64(cell.eps_r[2][2]),
                R::default(),
            ],
            mu_si_row0: [
                R::from_f64(mu0_over_eps0 * cell.mu_r[0][0]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[0][1]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[0][2]),
                R::default(),
            ],
            mu_si_row1: [
                R::from_f64(mu0_over_eps0 * cell.mu_r[1][0]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[1][1]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[1][2]),
                R::default(),
            ],
            mu_si_row2: [
                R::from_f64(mu0_over_eps0 * cell.mu_r[2][0]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[2][1]),
                R::from_f64(mu0_over_eps0 * cell.mu_r[2][2]),
                R::default(),
            ],
            dof_offset: base as u32,
            n_p: n_p as u32,
            mass_offset,
            energy_mode,
            nl_cell_slot,
            _pad0: 0,
            _pad1: 0,
            _pad2: 0,
        });

        for a in 0..n_p {
            let dmat_row_offset = dmat_xyz.len() as u32;
            let row = a * n_p;
            for b in 0..n_p {
                dmat_xyz.push([
                    R::from_f64(cell.dx[row + b]),
                    R::from_f64(cell.dy[row + b]),
                    R::from_f64(cell.dz[row + b]),
                    R::default(),
                ]);
            }
            let d = base + a;
            let qpm_sign = cell.nonlinear.as_ref().and_then(|nl| nl.node_qpm_signs.get(a).copied()).unwrap_or(1.0);
            dofs_meta[d] = GpuDofMeta {
                dipole_e_vec: [
                    R::from_f64(dip_e_per_dof[d][0]),
                    R::from_f64(dip_e_per_dof[d][1]),
                    R::from_f64(dip_e_per_dof[d][2]),
                    R::from_f64(qpm_sign),
                ],
                dipole_h_vec: [
                    R::from_f64(dip_h_per_dof[d][0]),
                    R::from_f64(dip_h_per_dof[d][1]),
                    R::from_f64(dip_h_per_dof[d][2]),
                    R::default(),
                ],
                cell_id: ci as u32,
                cell_dof_offset: base as u32,
                n_p: n_p as u32,
                dmat_row_offset,
                lift_start: 0,
                lift_len: 0,
                _pad0: 0,
                _pad1: 0,
            };
        }

        for face in &cell.faces {
            let (face_kind, nb_base) = match face.kind {
                DgtdFaceKind::Interior { neighbor_cell: nb } => (0u32, asm.cells[nb].dof_offset),
                DgtdFaceKind::Pec => (1u32, 0),
                DgtdFaceKind::Pmc => (2u32, 0),
                DgtdFaceKind::SilverMuller => (3u32, 0),
                DgtdFaceKind::Port { port_index: pi } => {
                    if pi == excited_port_idx {
                        (4u32, 0)
                    } else {
                        (5u32, 0)
                    }
                }
            };
            let (has_sheet2d, sheet_sigma, sheet_u_hat, sheet_v_hat) = if let Some((sig, u_hat, v_hat)) = face.sheet_2d
            {
                (
                    1u32,
                    [R::from_f64(sig[0][0]), R::from_f64(sig[0][1]), R::from_f64(sig[1][0]), R::from_f64(sig[1][1])],
                    [R::from_f64(u_hat[0]), R::from_f64(u_hat[1]), R::from_f64(u_hat[2]), R::default()],
                    [R::from_f64(v_hat[0]), R::from_f64(v_hat[1]), R::from_f64(v_hat[2]), R::default()],
                )
            } else {
                (0u32, [R::default(); 4], [R::default(); 4], [R::default(); 4])
            };

            for qp in &face.qps {
                let qp_idx = face_qps.len() as u32;
                let trace_minus_start = trace_pairs.len() as u32;
                for &(b, w) in &qp.trace_minus {
                    trace_pairs.push(GpuTracePair { global_dof: (base + b) as u32, weight: R::from_f64(w), _pad: 0 });
                }
                let trace_minus_len = qp.trace_minus.len() as u32;

                let trace_plus_start = trace_pairs.len() as u32;
                if face_kind == 0 {
                    for &(b, w) in &qp.trace_plus {
                        trace_pairs.push(GpuTracePair {
                            global_dof: (nb_base + b) as u32,
                            weight: R::from_f64(w),
                            _pad: 0,
                        });
                    }
                }
                let trace_plus_len = if face_kind == 0 { qp.trace_plus.len() as u32 } else { 0 };

                let (port_e_mode, port_h_mode) = if let Some((em, hm)) = qp.port_mode {
                    (
                        [R::from_f64(em[0]), R::from_f64(em[1]), R::from_f64(em[2]), R::default()],
                        [R::from_f64(hm[0]), R::from_f64(hm[1]), R::from_f64(hm[2]), R::default()],
                    )
                } else {
                    ([R::default(); 4], [R::default(); 4])
                };

                let (has_tfsf, tfsf_delay, tfsf_e_unit, tfsf_h_unit) =
                    if let Some((sign, kh, p_vec, origin, c_bg, y_bg)) = qp.tfsf {
                        let delay = dot3(sub3(qp.x, origin), kh) / c_bg;
                        let k_cross_p = cross3(kh, p_vec);
                        (
                            1u32,
                            R::from_f64(delay),
                            [
                                R::from_f64(sign * p_vec[0]),
                                R::from_f64(sign * p_vec[1]),
                                R::from_f64(sign * p_vec[2]),
                                R::default(),
                            ],
                            [
                                R::from_f64(sign * y_bg * k_cross_p[0]),
                                R::from_f64(sign * y_bg * k_cross_p[1]),
                                R::from_f64(sign * y_bg * k_cross_p[2]),
                                R::default(),
                            ],
                        )
                    } else {
                        (0u32, R::default(), [R::default(); 4], [R::default(); 4])
                    };

                face_qps.push(GpuFaceQpMeta {
                    n_hat_ym: [
                        R::from_f64(qp.n_hat[0]),
                        R::from_f64(qp.n_hat[1]),
                        R::from_f64(qp.n_hat[2]),
                        R::from_f64(cell.y_admittance),
                    ],
                    yp_tfsf_delay: [R::from_f64(face.y_plus), tfsf_delay, R::default(), R::default()],
                    port_e_mode,
                    port_h_mode,
                    tfsf_e_unit,
                    tfsf_h_unit,
                    sheet_sigma,
                    sheet_u_hat,
                    sheet_v_hat,
                    trace_minus_start,
                    trace_minus_len,
                    trace_plus_start,
                    trace_plus_len,
                    face_kind,
                    has_sheet2d,
                    has_tfsf,
                    _pad0: 0,
                });

                for &(a, lift_w) in &qp.lift_minus {
                    per_dof_lifts[base + a].push(GpuDofLift { qp_idx, lift_w: R::from_f64(lift_w), _pad: 0 });
                }
            }
        }
    }

    let mut dof_lifts = Vec::new();
    for (d, lifts) in per_dof_lifts.into_iter().enumerate() {
        dofs_meta[d].lift_start = dof_lifts.len() as u32;
        dofs_meta[d].lift_len = lifts.len() as u32;
        dof_lifts.extend(lifts);
    }

    let mut probes_meta = Vec::with_capacity(n_probes.max(1));
    let mut probe_pairs = Vec::new();
    for pr in &asm.probes {
        let base = asm.cells[pr.cell_index].dof_offset;
        let pair_start = probe_pairs.len() as u32;
        for &(b, w) in &pr.basis_weights {
            probe_pairs.push(GpuTracePair { global_dof: (base + b) as u32, weight: R::from_f64(w), _pad: 0 });
        }
        probes_meta.push(GpuProbeMeta { pair_start, pair_len: pr.basis_weights.len() as u32, _pad: [0; 2] });
    }

    let total_qps = face_qps.len();
    let n_steps = ((asm.config.t_final / dt).ceil() as usize).max(1);
    let n_snaps = n_steps + 1;
    let n_energy_wgs = (nc as u32).div_ceil(WG_SIZE_1D).max(1);

    // 4. Allocate device buffers
    let device = &ctx.device;
    let queue = &ctx.queue;
    let scalar_bytes = std::mem::size_of::<R>();
    let vec4_bytes = 4 * scalar_bytes;

    let make_rw_vec4_buf = |label: &str, init: &[[R; 4]]| -> wgpu::Buffer {
        device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytemuck::cast_slice(init),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
        })
    };
    let make_zero_vec4_buf = |label: &str, len: usize| -> wgpu::Buffer {
        device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: (len.max(1) * vec4_bytes) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    };
    let make_ro_buf = |label: &str, bytes: &[u8]| -> wgpu::Buffer {
        device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytes,
            usage: wgpu::BufferUsages::STORAGE,
        })
    };

    let (pulse_kind, pulse_f0, pulse_fwidth, pulse_t0, pulse_amp, pulse_phase) = match asm.config.pulse {
        DgtdPulse::GaussianModulated { f0_hz, fwidth_hz, amplitude, phase_rad, .. } => {
            (0u32, f0_hz, fwidth_hz, asm.config.pulse.t0_s(), amplitude, phase_rad)
        }
        DgtdPulse::Ricker { f0_hz, amplitude, .. } => {
            (1u32, f0_hz, 0.5 * f0_hz, asm.config.pulse.t0_s(), amplitude, 0.0)
        }
    };
    let pol_mode = match asm.config.polarization {
        DgtdPolarization::Full => 0u32,
        DgtdPolarization::Tm => 1u32,
        DgtdPolarization::Te => 2u32,
    };

    let mut uniforms = DgtdStepUniforms {
        dt: R::from_f64(dt),
        t_now: R::default(),
        pulse_now: R::from_f64(asm.config.pulse.eval(0.0)),
        flux_alpha: R::from_f64(asm.config.flux.alpha()),
        rk_a: R::default(),
        rk_b: R::default(),
        rk_w_acc: R::default(),
        rk_h_next: R::default(),
        pulse_f0: R::from_f64(pulse_f0),
        pulse_fwidth: R::from_f64(pulse_fwidth),
        pulse_t0: R::from_f64(pulse_t0),
        pulse_amp: R::from_f64(pulse_amp),
        pulse_phase: R::from_f64(pulse_phase),
        symmetry_factor: R::from_f64(asm.symmetry_factor),
        inv_eps0: R::from_f64(1.0 / EPS0),
        e_ref: R::from_f64(e_ref),
        inv_e_ref: R::from_f64(1.0 / e_ref),
        _pad_r0: R::default(),
        _pad_r1: R::default(),
        _pad_r2: R::default(),
        total_dofs: ndofs as u32,
        total_qps: total_qps as u32,
        n_cells: nc as u32,
        n_freqs: n_freqs as u32,
        n_probes: n_probes as u32,
        n_energy_wgs,
        dim: asm.dim as u32,
        pol_mode,
        rk_stage: 0,
        snap_idx: 0,
        pulse_kind,
        has_nonlinear: u32::from(has_nonlinear),
        has_carriers: u32::from(has_carriers),
        _pad_u0: 0,
        _pad_u1: 0,
        _pad_u2: 0,
    };

    let uniform_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("dgtd-uniforms"),
        contents: bytemuck::bytes_of(&uniforms),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    let state_e_buf = make_rw_vec4_buf("dgtd-state-e", &init_e);
    let state_h_buf = make_rw_vec4_buf("dgtd-state-h", &init_h);
    let de_out_buf = make_zero_vec4_buf("dgtd-de-out", ndofs);
    let dh_out_buf = make_zero_vec4_buf("dgtd-dh-out", ndofs);
    let aux0_e_buf = make_zero_vec4_buf("dgtd-aux0-e", ndofs);
    let aux0_h_buf = make_zero_vec4_buf("dgtd-aux0-h", ndofs);
    let aux1_e_buf = make_zero_vec4_buf("dgtd-aux1-e", ndofs);
    let aux1_h_buf = make_zero_vec4_buf("dgtd-aux1-h", ndofs);
    let qp_flux_e_buf = make_zero_vec4_buf("dgtd-qp-flux-e", total_qps.max(1));
    let qp_flux_h_buf = make_zero_vec4_buf("dgtd-qp-flux-h", total_qps.max(1));

    let dummy_cell = [GpuCellMeta::<R>::default()];
    let dummy_dof = [GpuDofMeta::<R>::default()];
    let dummy_dmat = [[R::default(); 4]];
    let dummy_lift = [GpuDofLift::<R>::default()];
    let dummy_qp = [GpuFaceQpMeta::<R>::default()];
    let dummy_tp = [GpuTracePair::<R>::default()];
    let dummy_scalar = [R::default(); 4];
    let dummy_probe = [GpuProbeMeta::default()];
    let dummy_nl_cell = [GpuNlCellMeta::<R>::default()];

    let cells_meta_buf = make_ro_buf("dgtd-cells-meta", &non_empty_typed_bytes(&cells_meta, &dummy_cell));
    let dofs_meta_buf = make_ro_buf("dgtd-dofs-meta", &non_empty_typed_bytes(&dofs_meta, &dummy_dof));
    let dmat_xyz_buf = make_ro_buf("dgtd-dmat-xyz", &non_empty_typed_bytes(&dmat_xyz, &dummy_dmat));
    let dof_lifts_buf = make_ro_buf("dgtd-dof-lifts", &non_empty_typed_bytes(&dof_lifts, &dummy_lift));
    let face_qps_buf = make_ro_buf("dgtd-face-qps", &non_empty_typed_bytes(&face_qps, &dummy_qp));
    let trace_pairs_buf = make_ro_buf("dgtd-trace-pairs", &non_empty_typed_bytes(&trace_pairs, &dummy_tp));
    let cell_masses_buf = make_ro_buf("dgtd-cell-masses", &non_empty_typed_bytes(&cell_masses, &dummy_scalar));

    let dft_cos_sin_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dgtd-dft-cos-sin"),
        size: (n_freqs.max(1) * 2 * scalar_bytes) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let dft_len = (n_freqs * ndofs).max(1);
    let dft_e_re_buf = make_zero_vec4_buf("dgtd-dft-e-re", dft_len);
    let dft_e_im_buf = make_zero_vec4_buf("dgtd-dft-e-im", dft_len);
    let dft_h_re_buf = make_zero_vec4_buf("dgtd-dft-h-re", dft_len);
    let dft_h_im_buf = make_zero_vec4_buf("dgtd-dft-h-im", dft_len);

    let probes_meta_buf = make_ro_buf("dgtd-probes-meta", &non_empty_typed_bytes(&probes_meta, &dummy_probe));
    let probe_pairs_buf = make_ro_buf("dgtd-probe-pairs", &non_empty_typed_bytes(&probe_pairs, &dummy_tp));
    let probe_out_len = (n_snaps * n_probes).max(1);
    let probe_out_e_buf = make_zero_vec4_buf("dgtd-probe-out-e", probe_out_len);
    let probe_out_h_buf = make_zero_vec4_buf("dgtd-probe-out-h", probe_out_len);

    let energy_partials_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dgtd-energy-partials"),
        size: (n_energy_wgs as usize * scalar_bytes).max(32) as u64,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    let energy_history_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dgtd-energy-history"),
        size: (n_snaps * scalar_bytes).max(32) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let nl_cells_meta_buf = make_ro_buf("dgtd-nl-cells", &non_empty_typed_bytes(&nl_cells_meta, &dummy_nl_cell));
    let nl_rank4_buf = make_ro_buf("dgtd-nl-rank4", &non_empty_typed_bytes(&nl_rank4_flat, &dummy_scalar));
    let nl_carrier_buf = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("dgtd-nl-carrier"),
        size: (if has_carriers { ndofs.max(1) * scalar_bytes } else { 32 }).max(32) as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // 5. Create bind groups and compute pipelines
    let bgl_entry = |binding: u32, read_only: bool, is_uniform: bool| -> wgpu::BindGroupLayoutEntry {
        wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: if is_uniform {
                    wgpu::BufferBindingType::Uniform
                } else {
                    wgpu::BufferBindingType::Storage { read_only }
                },
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }
    };

    let bgl0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("dgtd-bgl0"),
        entries: &[
            bgl_entry(0, true, true),
            bgl_entry(1, false, false),
            bgl_entry(2, false, false),
            bgl_entry(3, false, false),
            bgl_entry(4, false, false),
            bgl_entry(5, false, false),
            bgl_entry(6, false, false),
            bgl_entry(7, false, false),
            bgl_entry(8, false, false),
            bgl_entry(9, false, false),
            bgl_entry(10, false, false),
        ],
    });
    let bgl1 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("dgtd-bgl1"),
        entries: &[
            bgl_entry(0, true, false),
            bgl_entry(1, true, false),
            bgl_entry(2, true, false),
            bgl_entry(3, true, false),
            bgl_entry(4, true, false),
            bgl_entry(5, true, false),
            bgl_entry(6, true, false),
        ],
    });
    let bgl2 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("dgtd-bgl2"),
        entries: &[
            bgl_entry(0, true, false),
            bgl_entry(1, false, false),
            bgl_entry(2, false, false),
            bgl_entry(3, false, false),
            bgl_entry(4, false, false),
            bgl_entry(5, true, false),
            bgl_entry(6, true, false),
            bgl_entry(7, false, false),
            bgl_entry(8, false, false),
            bgl_entry(9, false, false),
            bgl_entry(10, false, false),
        ],
    });
    let bgl3 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("dgtd-bgl3"),
        entries: &[bgl_entry(0, true, false), bgl_entry(1, true, false), bgl_entry(2, false, false)],
    });

    let bg0 = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("dgtd-bg0"),
        layout: &bgl0,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: uniform_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: state_e_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: state_h_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: de_out_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: dh_out_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: aux0_e_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: aux0_h_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 7, resource: aux1_e_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 8, resource: aux1_h_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 9, resource: qp_flux_e_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 10, resource: qp_flux_h_buf.as_entire_binding() },
        ],
    });
    let bg1 = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("dgtd-bg1"),
        layout: &bgl1,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: cells_meta_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: dofs_meta_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: dmat_xyz_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: dof_lifts_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: face_qps_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: trace_pairs_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: cell_masses_buf.as_entire_binding() },
        ],
    });
    let bg2 = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("dgtd-bg2"),
        layout: &bgl2,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: dft_cos_sin_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: dft_e_re_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: dft_e_im_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: dft_h_re_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: dft_h_im_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: probes_meta_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: probe_pairs_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 7, resource: probe_out_e_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 8, resource: probe_out_h_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 9, resource: energy_partials_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 10, resource: energy_history_buf.as_entire_binding() },
        ],
    });
    let bg3 = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("dgtd-bg3"),
        layout: &bgl3,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: nl_cells_meta_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: nl_rank4_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: nl_carrier_buf.as_entire_binding() },
        ],
    });

    let shader_src = build_dgtd_wgsl(use_f64);
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("nyx-ems-dgtd-wgsl"),
        source: wgpu::ShaderSource::Wgsl(shader_src.into()),
    });
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("nyx-ems-dgtd-pl"),
        bind_group_layouts: &[Some(&bgl0), Some(&bgl1), Some(&bgl2), Some(&bgl3)],
        immediate_size: 0,
    });
    let make_pipeline = |entry: &str| -> wgpu::ComputePipeline {
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some(entry),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        })
    };

    let p_eval_face_qps = make_pipeline("eval_face_qps");
    let p_eval_dof_rhs = make_pipeline("eval_dof_rhs");
    let p_lserk4_stage = make_pipeline("lserk4_stage");
    let p_rk4_save_init = make_pipeline("rk4_save_init");
    let p_rk4_stage = make_pipeline("rk4_stage");
    let p_step_carriers = make_pipeline("step_carriers");
    let p_accum_dof_dft = make_pipeline("accum_dof_dft");
    let p_sample_probes = make_pipeline("sample_probes");
    let p_energy_pass1 = make_pipeline("dgtd_energy_pass1");
    let p_energy_pass2 = make_pipeline("dgtd_energy_pass2");

    let wg_dofs = (ndofs as u32).div_ceil(WG_SIZE_1D);
    let wg_qps = (total_qps as u32).div_ceil(WG_SIZE_1D);
    let wg_probes = (n_probes as u32).div_ceil(WG_SIZE_1D);

    let mut pulse_dft = vec![c(0.0, 0.0); n_freqs];
    let mut dft_cos_sin_host = vec![[R::default(); 2]; n_freqs.max(1)];
    let mut times_s = Vec::with_capacity(n_snaps);

    let enqueue_snapshot_and_dft = |t_now: f64,
                                    dt_w: f64,
                                    snap_idx: usize,
                                    u_state: &mut DgtdStepUniforms<R>,
                                    pulse_dft: &mut [C],
                                    cos_sin_host: &mut [[R; 2]],
                                    times: &mut Vec<f64>| {
        times.push(t_now);
        u_state.t_now = R::from_f64(t_now);
        u_state.pulse_now = R::from_f64(asm.config.pulse.eval(t_now));
        u_state.snap_idx = snap_idx as u32;
        queue.write_buffer(&uniform_buf, 0, bytemuck::bytes_of(u_state));

        if n_freqs > 0 {
            let s_now = asm.config.pulse.eval(t_now);
            for (fi, &freq) in asm.config.freqs_hz.iter().enumerate() {
                let angle = -2.0 * PI * freq * t_now;
                let (s_ang, c_ang) = angle.sin_cos();
                let ph = c(c_ang * dt_w, s_ang * dt_w);
                pulse_dft[fi] += ph * s_now;
                cos_sin_host[fi] = [R::from_f64(c_ang * dt_w), R::from_f64(s_ang * dt_w)];
            }
            queue.write_buffer(&dft_cos_sin_buf, 0, bytemuck::cast_slice(cos_sin_host));
        }

        let mut enc =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("dgtd-snap-dft-enc") });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("dgtd-snap-dft-pass"),
                timestamp_writes: None,
            });
            pass.set_bind_group(0, &bg0, &[]);
            pass.set_bind_group(1, &bg1, &[]);
            pass.set_bind_group(2, &bg2, &[]);
            pass.set_bind_group(3, &bg3, &[]);
            if n_freqs > 0 {
                pass.set_pipeline(&p_accum_dof_dft);
                pass.dispatch_workgroups(wg_dofs, 1, 1);
            }
            if n_probes > 0 {
                pass.set_pipeline(&p_sample_probes);
                pass.dispatch_workgroups(wg_probes, 1, 1);
            }
            pass.set_pipeline(&p_energy_pass1);
            pass.dispatch_workgroups(n_energy_wgs, 1, 1);
            pass.set_pipeline(&p_energy_pass2);
            pass.dispatch_workgroups(1, 1, 1);
        }
        queue.submit(Some(enc.finish()));
    };

    // Initial snapshot and half-weight trapezoidal DFT at t = 0
    enqueue_snapshot_and_dft(0.0, 0.5 * dt, 0, &mut uniforms, &mut pulse_dft, &mut dft_cos_sin_host, &mut times_s);

    let mut total_element_updates = 0u64;
    let use_lserk4 = asm.config.time_stepper == DgtdTimeStepper::Lserk4;
    let stage_weights = [dt / 6.0, dt / 3.0, dt / 3.0, dt / 6.0];
    let stage_steps = [0.5 * dt, 0.5 * dt, dt, 0.0];

    for step in 0..n_steps {
        let t_s = (step as f64) * dt;
        if use_lserk4 {
            for stage in 0..5 {
                let tau = t_s + LSERK4_C[stage] * dt;
                uniforms.t_now = R::from_f64(tau);
                uniforms.pulse_now = R::from_f64(asm.config.pulse.eval(tau));
                uniforms.rk_stage = stage as u32;
                uniforms.rk_a = R::from_f64(LSERK4_A[stage]);
                uniforms.rk_b = R::from_f64(LSERK4_B[stage]);
                queue.write_buffer(&uniform_buf, 0, bytemuck::bytes_of(&uniforms));

                let mut enc =
                    device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("dgtd-lserk4-stage") });
                {
                    let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: Some("dgtd-lserk4-pass"),
                        timestamp_writes: None,
                    });
                    pass.set_bind_group(0, &bg0, &[]);
                    pass.set_bind_group(1, &bg1, &[]);
                    pass.set_bind_group(2, &bg2, &[]);
                    pass.set_bind_group(3, &bg3, &[]);
                    if total_qps > 0 {
                        pass.set_pipeline(&p_eval_face_qps);
                        pass.dispatch_workgroups(wg_qps, 1, 1);
                    }
                    pass.set_pipeline(&p_eval_dof_rhs);
                    pass.dispatch_workgroups(wg_dofs, 1, 1);
                    pass.set_pipeline(&p_lserk4_stage);
                    pass.dispatch_workgroups(wg_dofs, 1, 1);
                    if has_carriers && stage == 4 {
                        pass.set_pipeline(&p_step_carriers);
                        pass.dispatch_workgroups(wg_dofs, 1, 1);
                    }
                }
                queue.submit(Some(enc.finish()));
                total_element_updates += nc as u64;
            }
        } else {
            let stage_times = [t_s, t_s + 0.5 * dt, t_s + 0.5 * dt, t_s + dt];
            for s in 0..4 {
                let tau = stage_times[s];
                uniforms.t_now = R::from_f64(tau);
                uniforms.pulse_now = R::from_f64(asm.config.pulse.eval(tau));
                uniforms.rk_stage = s as u32;
                uniforms.rk_w_acc = R::from_f64(stage_weights[s]);
                uniforms.rk_h_next = R::from_f64(stage_steps[s]);
                queue.write_buffer(&uniform_buf, 0, bytemuck::bytes_of(&uniforms));

                let mut enc =
                    device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("dgtd-rk4-stage") });
                {
                    let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: Some("dgtd-rk4-pass"),
                        timestamp_writes: None,
                    });
                    pass.set_bind_group(0, &bg0, &[]);
                    pass.set_bind_group(1, &bg1, &[]);
                    pass.set_bind_group(2, &bg2, &[]);
                    pass.set_bind_group(3, &bg3, &[]);
                    if s == 0 {
                        pass.set_pipeline(&p_rk4_save_init);
                        pass.dispatch_workgroups(wg_dofs, 1, 1);
                    }
                    if total_qps > 0 {
                        pass.set_pipeline(&p_eval_face_qps);
                        pass.dispatch_workgroups(wg_qps, 1, 1);
                    }
                    pass.set_pipeline(&p_eval_dof_rhs);
                    pass.dispatch_workgroups(wg_dofs, 1, 1);
                    pass.set_pipeline(&p_rk4_stage);
                    pass.dispatch_workgroups(wg_dofs, 1, 1);
                    if has_carriers && s == 3 {
                        pass.set_pipeline(&p_step_carriers);
                        pass.dispatch_workgroups(wg_dofs, 1, 1);
                    }
                }
                queue.submit(Some(enc.finish()));
                total_element_updates += nc as u64;
            }
        }

        let t_next = ((step + 1) as f64) * dt;
        let dt_w = if step + 1 == n_steps { 0.5 * dt } else { dt };
        enqueue_snapshot_and_dft(
            t_next,
            dt_w,
            step + 1,
            &mut uniforms,
            &mut pulse_dft,
            &mut dft_cos_sin_host,
            &mut times_s,
        );
    }

    // 6. Read back final state, DFT spectra, probe time series, and normalized energy history
    let final_e_raw: Vec<[R; 4]> = read_back_slice(device, queue, &state_e_buf, ndofs)?;
    let final_h_raw: Vec<[R; 4]> = read_back_slice(device, queue, &state_h_buf, ndofs)?;
    let mut final_state = DgtdState::zeros(ndofs, 0);
    for d in 0..ndofs {
        final_state.e[d] = [final_e_raw[d][0].to_f64(), final_e_raw[d][1].to_f64(), final_e_raw[d][2].to_f64()];
        final_state.h[d] = [final_h_raw[d][0].to_f64(), final_h_raw[d][1].to_f64(), final_h_raw[d][2].to_f64()];
    }
    if has_carriers {
        let nc_raw: Vec<R> = read_back_slice(device, queue, &nl_carrier_buf, ndofs)?;
        final_state.n_carrier = nc_raw.into_iter().map(|v| v.to_f64() * nc_ref).collect();
    }

    let mut dft_e = vec![vec![[c(0.0, 0.0); 3]; ndofs]; n_freqs];
    let mut dft_h = vec![vec![[c(0.0, 0.0); 3]; ndofs]; n_freqs];
    if n_freqs > 0 {
        let e_re_raw: Vec<[R; 4]> = read_back_slice(device, queue, &dft_e_re_buf, n_freqs * ndofs)?;
        let e_im_raw: Vec<[R; 4]> = read_back_slice(device, queue, &dft_e_im_buf, n_freqs * ndofs)?;
        let h_re_raw: Vec<[R; 4]> = read_back_slice(device, queue, &dft_h_re_buf, n_freqs * ndofs)?;
        let h_im_raw: Vec<[R; 4]> = read_back_slice(device, queue, &dft_h_im_buf, n_freqs * ndofs)?;
        for fi in 0..n_freqs {
            let off = fi * ndofs;
            for d in 0..ndofs {
                let er = e_re_raw[off + d];
                let ei = e_im_raw[off + d];
                let hr = h_re_raw[off + d];
                let hi = h_im_raw[off + d];
                dft_e[fi][d] = [
                    c(er[0].to_f64(), ei[0].to_f64()),
                    c(er[1].to_f64(), ei[1].to_f64()),
                    c(er[2].to_f64(), ei[2].to_f64()),
                ];
                dft_h[fi][d] = [
                    c(hr[0].to_f64(), hi[0].to_f64()),
                    c(hr[1].to_f64(), hi[1].to_f64()),
                    c(hr[2].to_f64(), hi[2].to_f64()),
                ];
            }
        }
    }

    let energy_raw: Vec<R> = read_back_slice(device, queue, &energy_history_buf, n_snaps)?;
    let energy_scale = EPS0 * mean_measure;
    let energy_j: Vec<f64> = energy_raw.into_iter().map(|v| v.to_f64() * energy_scale).collect();

    let mut probe_e = vec![Vec::with_capacity(n_snaps); n_probes];
    let mut probe_h = vec![Vec::with_capacity(n_snaps); n_probes];
    if n_probes > 0 {
        let pe_raw: Vec<[R; 4]> = read_back_slice(device, queue, &probe_out_e_buf, n_snaps * n_probes)?;
        let ph_raw: Vec<[R; 4]> = read_back_slice(device, queue, &probe_out_h_buf, n_snaps * n_probes)?;
        for snap in 0..n_snaps {
            for pi in 0..n_probes {
                let idx = snap * n_probes + pi;
                probe_e[pi].push([pe_raw[idx][0].to_f64(), pe_raw[idx][1].to_f64(), pe_raw[idx][2].to_f64()]);
                probe_h[pi].push([ph_raw[idx][0].to_f64(), ph_raw[idx][1].to_f64(), ph_raw[idx][2].to_f64()]);
            }
        }
    }

    // 7. Reconstruct port modal DFTs (V_p, I_p) and surface FluxMonitor DFTs from nodal (dft_e, dft_h)
    let mut port_faces: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_ports];
    for ci in 0..nc {
        for (lf_idx, f) in asm.cells[ci].faces.iter().enumerate() {
            if let DgtdFaceKind::Port { port_index: pi } = f.kind {
                port_faces[pi].push((ci, lf_idx));
            }
        }
    }

    let mut port_dft = vec![vec![(c(0.0, 0.0), c(0.0, 0.0)); n_ports]; n_freqs];
    for fi in 0..n_freqs {
        let s_hat = pulse_dft[fi];
        let de_f = &dft_e[fi];
        let dh_f = &dft_h[fi];
        for pi in 0..n_ports {
            let mut v_acc = c(0.0, 0.0);
            let mut i_acc = c(0.0, 0.0);
            for &(ci, lf_idx) in &port_faces[pi] {
                let cell = &asm.cells[ci];
                let base = cell.dof_offset;
                let imu = 1.0 / ((cell.mu_r[0][0] + cell.mu_r[1][1] + cell.mu_r[2][2]) / 3.0).max(1e-6);
                let face = &cell.faces[lf_idx];
                let y_m = cell.y_admittance;
                let y_p = face.y_plus;
                let y_sum = (y_m + y_p).max(1e-30);
                for qp in &face.qps {
                    let Some((e_mode, _)) = qp.port_mode else { continue };
                    let mut e_m = [c(0.0, 0.0); 3];
                    let mut h_m = [c(0.0, 0.0); 3];
                    for &(b, w) in &qp.trace_minus {
                        let eb = de_f[base + b];
                        let hb = dh_f[base + b];
                        for k in 0..3 {
                            e_m[k] += eb[k] * w;
                            h_m[k] += hb[k] * w;
                        }
                    }
                    let (e_ext, h_ext) = if pi == excited_port_idx {
                        let hm = qp.port_mode.map_or([0.0; 3], |(_, h)| h);
                        (
                            [s_hat * e_mode[0], s_hat * e_mode[1], s_hat * e_mode[2]],
                            [s_hat * hm[0], s_hat * hm[1], s_hat * hm[2]],
                        )
                    } else {
                        ([c(0.0, 0.0); 3], [c(0.0, 0.0); 3])
                    };
                    let dh = [h_ext[0] - h_m[0], h_ext[1] - h_m[1], h_ext[2] - h_m[2]];
                    let n = qp.n_hat;
                    let n_cross_dh =
                        [dh[2] * n[1] - dh[1] * n[2], dh[0] * n[2] - dh[2] * n[0], dh[1] * n[0] - dh[0] * n[1]];
                    let e_star = [
                        (e_m[0] * y_m + e_ext[0] * y_p + n_cross_dh[0]) / y_sum,
                        (e_m[1] * y_m + e_ext[1] * y_p + n_cross_dh[1]) / y_sum,
                        (e_m[2] * y_m + e_ext[2] * y_p + n_cross_dh[2]) / y_sum,
                    ];
                    let n_cross_hm =
                        [h_m[2] * n[1] - h_m[1] * n[2], h_m[0] * n[2] - h_m[2] * n[0], h_m[1] * n[0] - h_m[0] * n[1]];
                    let n_cross_h_star = [
                        n_cross_hm[0] - (e_star[0] - e_m[0]) * y_m,
                        n_cross_hm[1] - (e_star[1] - e_m[1]) * y_m,
                        n_cross_hm[2] - (e_star[2] - e_m[2]) * y_m,
                    ];
                    let w_imu = qp.w_ds * imu;
                    v_acc += (e_star[0] * e_mode[0] + e_star[1] * e_mode[1] + e_star[2] * e_mode[2]) * w_imu;
                    i_acc +=
                        (n_cross_h_star[0] * e_mode[0] + n_cross_h_star[1] * e_mode[1] + n_cross_h_star[2] * e_mode[2])
                            * w_imu;
                }
            }
            port_dft[fi][pi] = (v_acc, i_acc);
        }
    }

    let mut s_carrier_dft = c(0.0, 0.0);
    {
        let f_c = asm.config.pulse.f0_hz();
        for step in 0..=n_steps {
            let t_s = (step as f64) * dt;
            let w_t = if step == 0 || step == n_steps { 0.5 * dt } else { dt };
            let ang = -2.0 * PI * f_c * t_s;
            let (sa, ca) = ang.sin_cos();
            s_carrier_dft += c(ca * w_t, sa * w_t) * asm.config.pulse.eval(t_s);
        }
    }
    let s_peak = pulse_dft.iter().map(|z| z.norm()).fold(s_carrier_dft.norm(), f64::max);
    let oob_floor = 1e-4 * s_peak;

    let mut flux_dft = vec![vec![(0.0, c(0.0, 0.0), c(0.0, 0.0), 0.0); n_freqs]; n_flux];
    for (mi, fm) in asm.flux_monitors.iter().enumerate() {
        for fi in 0..n_freqs {
            let s_hat = pulse_dft[fi];
            let s_norm = if s_hat.norm() >= oob_floor && s_hat.norm() > 1e-30 {
                s_hat
            } else if s_peak > 1e-30 {
                c(s_peak, 0.0)
            } else {
                c(1.0, 0.0)
            };
            let s_norm2 = s_norm.norm_sqr().max(1e-300);
            let mut raw_p_acc = 0.0;
            let mut e_acc = c(0.0, 0.0);
            let mut h_acc = c(0.0, 0.0);
            let mut area_acc = 0.0;
            let de_f = &dft_e[fi];
            let dh_f = &dft_h[fi];
            for &(ci, lf_idx) in &fm.cell_faces {
                let cell = &asm.cells[ci];
                let base = cell.dof_offset;
                let face = &cell.faces[lf_idx];
                for qp in &face.qps {
                    let mut eq = [c(0.0, 0.0); 3];
                    let mut hq = [c(0.0, 0.0); 3];
                    for &(b, w) in &qp.trace_minus {
                        let eb = de_f[base + b];
                        let hb = dh_f[base + b];
                        for k in 0..3 {
                            eq[k] += eb[k] * w;
                            hq[k] += hb[k] * w;
                        }
                    }
                    let poynting = [
                        (eq[1] * hq[2].conj() - eq[2] * hq[1].conj()).re,
                        (eq[2] * hq[0].conj() - eq[0] * hq[2].conj()).re,
                        (eq[0] * hq[1].conj() - eq[1] * hq[0].conj()).re,
                    ];
                    raw_p_acc += 0.5 * qp.w_ds * dot3(poynting, qp.n_hat);
                    let e_dom = if asm.dim == 2 && asm.config.polarization == DgtdPolarization::Tm {
                        eq[2]
                    } else {
                        eq[0] + eq[1] + eq[2]
                    };
                    let h_dom = hq[0] + hq[1] + hq[2];
                    e_acc += (e_dom / s_norm) * qp.w_ds;
                    h_acc += (h_dom / s_norm) * qp.w_ds;
                    area_acc += qp.w_ds;
                }
            }
            let inv_a = 1.0 / area_acc.max(1e-30);
            let p_acc = raw_p_acc / s_norm2;
            flux_dft[mi][fi] =
                (p_acc * asm.symmetry_factor, e_acc * inv_a, h_acc * inv_a, raw_p_acc * asm.symmetry_factor);
        }
    }

    let stages_per_step = if use_lserk4 { 5u64 } else { 4u64 };
    let equivalent_uniform_updates = (n_steps as u64) * (nc as u64) * stages_per_step;

    Ok(Ok(DgtdSolve {
        final_state,
        dft_e,
        dft_h,
        pulse_dft,
        port_dft,
        flux_dft,
        times_s,
        energy_j,
        probe_e,
        probe_h,
        n_macro_steps: n_steps,
        n_finest_steps: n_steps,
        total_element_updates,
        equivalent_uniform_updates,
        solve_time_s: t_start.elapsed().as_secs_f64(),
        backend_used: "gpu",
        gpu_adapter: Some(format!("{}, precision={}", ctx.adapter_description(), if use_f64 { "f64" } else { "f32" })),
        fallback_warning: None,
    }))
}
