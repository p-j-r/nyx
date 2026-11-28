// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! `wgpu` WGSL GPU compute kernels for `Study::Fdtd` (`#[cfg(feature = "gpu")]`, `ADR 146`).
//!
//! Operates directly over 64-byte-aligned Structure-of-Arrays (SoA) Yee field slabs uploaded as
//! contiguous `f32` or (when `wgpu::Features::SHADER_F64` is supported) native `f64` device buffers
//! (`Ex, Ey, Ez, Hx, Hy, Hz`, fused `(ca, cb)` and `(da, db)` update coefficients, 12-slab
//! convolutional `CPML` `psi` memory arrays, non-uniform `inv_dp`/`inv_dd` tables, Dey–Mittra
//! conformal PEC fractions, off-diagonal Farjadpour `PolarizedAveraging` couplings, synchronous
//! $\chi^{(2)}$/QPM/$\chi^{(3)}$ Kerr + TPA/FCA/FCD nonlinear polarization updates with dimensionless
//! $E_{\text{ref}}$ normalization, pre-evaluated source waveforms, running DFT monitors, and
//! atomic-free two-level tree reductions). Only monitor spectra, time traces, and energy scalars
//! are read back to host memory.
//!
//! **Future GPU & HPC Roadmap:** `wgpu`/WGSL serves as the zero-SDK portable Tier-1 compute
//! layer today and may be replaced by `CubeCL` or `rust-gpu` in the future once `f64` coverage
//! and stable toolchain compatibility mature (`femcore::linalg::device` remains backend-agnostic).
//! Dedicated HPC support and direct native bindings for NVIDIA (CUDA) and AMD (ROCm/HIP) GPUs
//! are planned for the future; if a user requires these on a priority basis, they should contact
//! the maintainers with target hardware and workload details.

use super::ade_pml::FdtdBcFace;
use super::grid::{FdtdDim, Pol2d, YeeFields, C0, EPS0};
use super::modes_bridge::rh_uv_axes;
use super::sources_monitors::{find_resonances_prony, Aux1dYeeLine, CompiledSource, StaggeredStencil8};
use super::{FdtdAssembly, FdtdSolve};
use crate::problem::TdPrecision;
use bytemuck::{Pod, Zeroable};
use femcore::linalg::csr::c;
use femcore::material::HBAR;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Instant;
use wgpu::util::DeviceExt;

const TPA_RMS_SCALE: f64 = 0.816_496_580_927_726;
const FDTD_NC_REF: f64 = 1.0e20;

/// Shared `wgpu` compute device and queue for `nyx-ems` time-domain solvers (`Study::Fdtd` & `Study::Dgtd`),
/// backed by [`femcore::linalg::gpu::GpuContext::shared`] so linear-algebra and time-domain solvers share
/// one `wgpu` adapter, logical device, and command queue.
pub struct EmsGpuContext {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub adapter_info: wgpu::AdapterInfo,
}

impl EmsGpuContext {
    /// Wraps an existing `wgpu::Device` and `wgpu::Queue` (e.g. from `femcore::linalg::gpu::GpuContext`).
    pub fn from_wgpu(device: wgpu::Device, queue: wgpu::Queue) -> Self {
        let adapter_info = device.adapter_info();
        Self { device, queue, adapter_info }
    }

    /// Creates an `EmsGpuContext` sharing the underlying `wgpu::Device` and `wgpu::Queue` of `ctx`.
    pub fn from_femcore(ctx: &femcore::linalg::gpu::GpuContext) -> Self {
        Self::from_wgpu(ctx.wgpu_device().clone(), ctx.wgpu_queue().clone())
    }

    /// True when the active adapter is a dedicated discrete GPU.
    #[inline]
    pub fn is_discrete(&self) -> bool {
        self.adapter_info.device_type == wgpu::DeviceType::DiscreteGpu
    }

    /// True when the active logical GPU device exposes `wgpu::Features::SHADER_F64`.
    #[inline]
    pub fn supports_f64(&self) -> bool {
        self.device.features().contains(wgpu::Features::SHADER_F64)
    }

    /// Human-readable adapter description (`"<name> (DiscreteGpu, Vulkan)"`).
    pub fn adapter_summary(&self) -> String {
        format!("{} ({:?}, {:?})", self.adapter_info.name, self.adapter_info.device_type, self.adapter_info.backend)
    }

    /// Alias for [`Self::adapter_summary`].
    pub fn adapter_description(&self) -> String {
        self.adapter_summary()
    }

    /// Acquires the process-wide cached `EmsGpuContext`, delegating adapter selection and logical
    /// device creation to [`femcore::linalg::gpu::GpuContext::shared`] so `femcore` and `nyx-ems`
    /// share a single `wgpu::Device` and `wgpu::Queue`.
    pub fn acquire_default() -> Result<Arc<Self>, String> {
        static CTX: OnceLock<Result<Arc<EmsGpuContext>, String>> = OnceLock::new();
        CTX.get_or_init(init_default_gpu_context).clone()
    }
}

fn init_default_gpu_context() -> Result<Arc<EmsGpuContext>, String> {
    let shared = femcore::linalg::gpu::GpuContext::shared().map_err(|e| e.to_string())?;
    Ok(Arc::new(EmsGpuContext::from_femcore(&shared)))
}

/// Resolves whether a time-domain GPU run executes in `f64` (`true`) or `f32` (`false`),
/// returning an explanatory fallback reason if explicit `f64` was requested on an adapter
/// without `SHADER_F64`.
pub(crate) fn select_td_precision(
    ctx: &EmsGpuContext,
    requested: TdPrecision,
    is_nonlinear: bool,
) -> Result<bool, String> {
    match requested {
        TdPrecision::F32 => Ok(false),
        TdPrecision::F64 => {
            if ctx.supports_f64() {
                Ok(true)
            } else {
                Err(format!(
                    "time-domain precision 'f64' requested, but GPU adapter '{}' does not expose SHADER_F64; falling back to CPU f64",
                    ctx.adapter_summary()
                ))
            }
        }
        TdPrecision::Auto => Ok(is_nonlinear && ctx.supports_f64()),
    }
}

pub(crate) trait GpuScalar: Pod + Zeroable + Copy + Default {
    fn from_f64(v: f64) -> Self;
    fn to_f64(self) -> f64;
}

impl GpuScalar for f32 {
    #[inline]
    fn from_f64(v: f64) -> Self {
        v as f32
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl GpuScalar for f64 {
    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self
    }
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Pod, Zeroable)]
struct GpuFdtdParams {
    nx: u32,
    ny: u32,
    nz: u32,
    nz_pad: u32,
    stride_x: u32,
    stride_y: u32,
    slab_len: u32,
    dim_mode: u32,
    use_conf: u32,
    has_cpml_x: u32,
    has_cpml_y: u32,
    has_cpml_z: u32,
    pmc_mask: u32,
    pec_mask: u32,
    per_mask: u32,
    off_inv_dp_x: u32,
    off_inv_dp_y: u32,
    off_inv_dp_z: u32,
    off_inv_dd_x: u32,
    off_inv_dd_y: u32,
    off_inv_dd_z: u32,
    off_cpml_e_x: u32,
    off_cpml_e_y: u32,
    off_cpml_e_z: u32,
    off_cpml_h_x: u32,
    off_cpml_h_y: u32,
    off_cpml_h_z: u32,
    n_src_h_sites: u32,
    n_src_h_chans: u32,
    n_src_e_sites: u32,
    n_src_e_chans: u32,
    n_offdiag: u32,
    n_planar_patches: u32,
    n_vol_pts: u32,
    n_vol_freqs: u32,
    vol_i0: u32,
    vol_i1: u32,
    vol_j0: u32,
    vol_j1: u32,
    vol_k0: u32,
    vol_k1: u32,
    vol_freq_offset: u32,
    total_dft_freqs: u32,
    n_time_monitors: u32,
    n_energy_wgs: u32,
    n_nl_edges: u32,
    nl_has_carriers: u32,
    _pad0: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct GpuStepUniform {
    step: u32,
    _pad: [u32; 63],
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Pod, Zeroable)]
struct GpuSrcSite {
    target_offset: u32,
    item_start: u32,
    item_end: u32,
    _pad: u32,
}

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

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuSrcTerm<R> {
    w_re: R,
    w_im: R,
    chan: u32,
    _pad: u32,
}
impl_pod_scalar_struct!(GpuSrcTerm, 16, 24);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuOffDiagEntry<R> {
    b_weights: [R; 4],
    c_weights: [R; 4],
    b_offsets: [u32; 4],
    c_offsets: [u32; 4],
    target_offset: u32,
    _pad: [u32; 7],
}
impl_pod_scalar_struct!(GpuOffDiagEntry, 96, 128);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuStencil8<R> {
    idx: [u32; 8],
    w0: [R; 4],
    w1: [R; 4],
}
impl_pod_scalar_struct!(GpuStencil8, 64, 96);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuPlanarPatch<R> {
    stencils: [GpuStencil8<R>; 4],
    freq_offset: u32,
    n_freqs: u32,
    out_offset: u32,
    patch_in_mon: u32,
    n_patches_in_mon: u32,
    _pad: [u32; 3],
}
impl_pod_scalar_struct!(GpuPlanarPatch, 288, 416);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuTimeMonitor<R> {
    stencil: GpuStencil8<R>,
    is_magnetic: u32,
    sample_stride: u32,
    write_offset: u32,
    _pad: [u32; 5],
}
impl_pod_scalar_struct!(GpuTimeMonitor, 96, 128);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuNlDomain<R> {
    d_glob: [R; 18],
    chi3_iso_xxxx: R,
    eps_min: R,
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
impl_pod_scalar_struct!(GpuNlDomain, 120, 224);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuNlEdge<R> {
    qpm_sign: R,
    cb_eps_inv_dt: R,
    cb_g: R,
    grid_idx: u32,
    ix: u32,
    iy: u32,
    iz: u32,
    comp: u32,
    dom_slot: u32,
    _pad0: u32,
    _pad1: u32,
}
impl_pod_scalar_struct!(GpuNlEdge, 44, 56);

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, Zeroable)]
struct GpuNlRunParams<R> {
    e_ref: R,
    inv_e_ref: R,
    _pad0: R,
    _pad1: R,
}
impl_pod_scalar_struct!(GpuNlRunParams, 16, 32);

/// Checks whether `asm` can execute directly on the GPU or should fall back to CPU with a warning.
///
/// Returns:
/// - `Ok(Ok(FdtdSolve))` when GPU execution succeeds.
/// - `Ok(Err(reason))` when `asm` uses a feature requiring CPU fallback (or no qualifying GPU is present).
/// - `Err(err)` on a fatal runtime error (such as cancellation or non-finite blowup).
pub fn try_solve_fdtd_gpu(asm: &FdtdAssembly, require_discrete: bool) -> Result<Result<FdtdSolve, String>, String> {
    if asm.is_complex_bloch {
        return Ok(Err(
            "FDTD complex Bloch periodicity (bloch_k != 0) requires dual complex grids; falling back to CPU FDTD"
                .into(),
        ));
    }
    if !asm.coeffs.ade_x.is_empty() || !asm.coeffs.ade_y.is_empty() || !asm.coeffs.ade_z.is_empty() {
        return Ok(Err("FDTD multi-pole Drude/Lorentz ADE dispersion bank is active; falling back to CPU FDTD".into()));
    }
    if asm.vol_dfts.len() > 1 {
        return Ok(Err("FDTD multiple volumetric FieldMonitors configured; falling back to CPU FDTD".into()));
    }

    let ctx = match EmsGpuContext::acquire_default() {
        Ok(c) => c,
        Err(e) => {
            return Ok(Err(format!("wgpu compute context unavailable ({e}); falling back to CPU FDTD")));
        }
    };
    if require_discrete && !ctx.is_discrete() {
        return Ok(Err(format!(
            "wgpu adapter '{}' is not a discrete GPU; keeping backend='auto' on CPU FDTD",
            ctx.adapter_info.name
        )));
    }

    let is_nl = asm.coeffs.nonlinear.is_active();
    let use_f64 = match select_td_precision(&ctx, asm.config.precision, is_nl) {
        Ok(v) => v,
        Err(reason) => return Ok(Err(reason)),
    };

    let elem_bytes = if use_f64 { std::mem::size_of::<f64>() } else { std::mem::size_of::<f32>() };
    let max_binding = ctx.device.limits().max_storage_buffer_binding_size as u64;
    let max_buf = ctx.device.limits().max_buffer_size;
    let limit_bytes = max_binding.min(max_buf);
    let slab_len = asm.grid.slab_len();
    let coeff_bytes = (3 * slab_len * 2 * elem_bytes) as u64;
    let psi_bytes = (6 * slab_len * elem_bytes) as u64;
    let vol_pts = asm.vol_dfts.first().map_or(0, |v| v.n_pts());
    let vol_freqs = asm.vol_dfts.first().map_or(0, |v| v.freqs_hz.len());
    let vol_dft_bytes = (vol_pts * vol_freqs * 3 * 4 * elem_bytes) as u64;
    if coeff_bytes > limit_bytes || psi_bytes > limit_bytes || vol_dft_bytes > limit_bytes {
        return Ok(Err(format!(
            "FDTD device buffer requirement ({:.1} MB) exceeds wgpu storage binding limit ({:.1} MB); falling back to CPU FDTD",
            (coeff_bytes.max(psi_bytes).max(vol_dft_bytes) as f64) / (1024.0 * 1024.0),
            (limit_bytes as f64) / (1024.0 * 1024.0)
        )));
    }

    let sol = if use_f64 {
        solve_fdtd_on_gpu_typed::<f64>(&ctx, asm, true)?
    } else {
        solve_fdtd_on_gpu_typed::<f32>(&ctx, asm, false)?
    };
    Ok(Ok(sol))
}

fn non_empty_typed_bytes<T: Pod>(slice: &[T], fallback_one: &[T; 1]) -> Vec<u8> {
    if slice.is_empty() {
        bytemuck::cast_slice(fallback_one).to_vec()
    } else {
        bytemuck::cast_slice(slice).to_vec()
    }
}

#[allow(clippy::type_complexity)]
fn compile_gpu_nl<R: GpuScalar>(
    asm: &FdtdAssembly,
    use_f64: bool,
    e_ref: f64,
    nc_ref: f64,
) -> (Vec<GpuNlDomain<R>>, Vec<R>, Vec<GpuNlEdge<R>>, u32) {
    let bank = &asm.coeffs.nonlinear;
    if !bank.is_active() {
        return (Vec::new(), Vec::new(), Vec::new(), 0);
    }
    let dt = asm.dt;
    let inv_dt = 1.0 / dt.max(1e-30);
    let _ = use_f64;

    let mut doms = Vec::with_capacity(bank.domains.len());
    let mut rank4_flat: Vec<R> = Vec::new();
    for d in &bank.domains {
        let mut d_glob = [R::default(); 18];
        let mut has_chi2 = 0u32;
        if let Some(dg) = d.chi2_d_global {
            has_chi2 = 1u32;
            for row in 0..3 {
                for col in 0..6 {
                    d_glob[row * 6 + col] = R::from_f64(dg[row][col] * e_ref);
                }
            }
        }
        let mut has_chi3_iso = 0u32;
        let mut rank4_offset = u32::MAX;
        if let Some(ref r4) = d.chi3_rank4_global {
            rank4_offset = rank4_flat.len() as u32;
            let e2 = e_ref * e_ref;
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        for l in 0..3 {
                            rank4_flat.push(R::from_f64(r4[i][j][k][l] * e2));
                        }
                    }
                }
            }
        } else if d.chi3_iso_xxxx != 0.0 {
            has_chi3_iso = 1u32;
        }
        let chi3_iso_xxxx = R::from_f64(d.chi3_iso_xxxx * e_ref * e_ref);
        let sig_tpa_coeff =
            R::from_f64((2.0 / 3.0) * EPS0 * EPS0 * C0 * C0 * d.n0 * d.n0 * d.beta_tpa * (e_ref * e_ref));
        let sig_fca_coeff = R::from_f64(EPS0 * C0 * d.n0 * d.sigma_fca * nc_ref);
        let fcd_coeff = R::from_f64(2.0 * d.n0 * d.k_fcd * nc_ref);

        let has_carrier = u32::from(d.beta_tpa > 0.0 && (d.sigma_fca > 0.0 || d.k_fcd != 0.0));
        let tpa_gen_coeff = if has_carrier != 0 && d.omega_0 > 0.0 {
            let pref = TPA_RMS_SCALE * C0 * EPS0 * d.n0 * (e_ref * e_ref);
            R::from_f64(((d.beta_tpa * pref * pref) / (2.0 * HBAR * d.omega_0)) / nc_ref)
        } else {
            R::default()
        };
        let (decay, one_minus_decay_tau) = if d.tau_c > 0.0 && d.tau_c.is_finite() {
            let x = dt / d.tau_c;
            if x > 700.0 {
                (R::default(), R::from_f64(d.tau_c))
            } else {
                (R::from_f64((-x).exp()), R::from_f64(d.tau_c * (-(-x).exp_m1())))
            }
        } else {
            (R::from_f64(1.0), R::from_f64(dt))
        };

        doms.push(GpuNlDomain {
            d_glob,
            chi3_iso_xxxx,
            eps_min: R::from_f64(d.eps_min),
            sig_tpa_coeff,
            sig_fca_coeff,
            fcd_coeff,
            tpa_gen_coeff,
            decay,
            one_minus_decay_tau,
            has_chi2,
            has_chi3_iso,
            rank4_offset,
            has_carrier,
        });
    }

    let cbs = [asm.coeffs.cb_x.as_slice(), asm.coeffs.cb_y.as_slice(), asm.coeffs.cb_z.as_slice()];
    let mut edges = Vec::with_capacity(bank.edges[0].len() + bank.edges[1].len() + bank.edges[2].len());
    for comp in 0..3usize {
        for e in &bank.edges[comp] {
            let cb_g = cbs[comp][e.idx];
            edges.push(GpuNlEdge {
                qpm_sign: R::from_f64(e.qpm_sign),
                cb_eps_inv_dt: R::from_f64(cb_g * EPS0 * inv_dt),
                cb_g: R::from_f64(cb_g),
                grid_idx: e.idx as u32,
                ix: e.ix,
                iy: e.iy,
                iz: e.iz,
                comp: comp as u32,
                dom_slot: e.dom_slot as u32,
                _pad0: 0,
                _pad1: 0,
            });
        }
    }

    (doms, rank4_flat, edges, u32::from(bank.has_carriers))
}

fn solve_fdtd_on_gpu_typed<R: GpuScalar>(
    ctx: &EmsGpuContext,
    asm: &FdtdAssembly,
    use_f64: bool,
) -> Result<FdtdSolve, String>
where
    GpuSrcTerm<R>: Pod,
    GpuOffDiagEntry<R>: Pod,
    GpuStencil8<R>: Pod,
    GpuPlanarPatch<R>: Pod,
    GpuTimeMonitor<R>: Pod,
    GpuNlDomain<R>: Pod,
    GpuNlEdge<R>: Pod,
    GpuNlRunParams<R>: Pod,
{
    let t0 = Instant::now();
    let grid = &asm.grid;
    let (nx, ny, nz, nz_pad, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.nz_pad, grid.stride_x, grid.stride_y);
    let slab_len = grid.slab_len();
    let n_steps = asm.n_steps;

    let dim_mode: u32 = match grid.dim {
        FdtdDim::ThreeD => 0,
        FdtdDim::TwoD { polarization: Pol2d::TmEy } => 1,
        FdtdDim::TwoD { polarization: Pol2d::TeHy } => 2,
        FdtdDim::TwoD { polarization: Pol2d::Full } => 3,
    };

    let bcs = &asm.config.boundaries;
    let mut pmc_mask = 0u32;
    let mut pec_mask = 0u32;
    for (i, f) in bcs.faces.iter().enumerate() {
        if matches!(f, FdtdBcFace::Pmc) {
            pmc_mask |= 1u32 << i;
        }
        if matches!(f, FdtdBcFace::Pec) {
            pec_mask |= 1u32 << i;
        }
    }
    let per = bcs.periodic_axes();
    let per_mask = (u32::from(per[0])) | (u32::from(per[1]) << 1) | (u32::from(per[2]) << 2);

    // 1. Pack 1D grid and CPML tables into `Vec<[R; 4]>`
    let mut tables: Vec<[R; 4]> = Vec::new();
    let push_scalar_table = |tables: &mut Vec<[R; 4]>, arr: &[f64]| -> u32 {
        let off = tables.len() as u32;
        for &v in arr {
            tables.push([R::from_f64(v), R::default(), R::default(), R::default()]);
        }
        off
    };
    let push_cpml_e = |tables: &mut Vec<[R; 4]>, ax: &super::ade_pml::CpmlAxis1d| -> u32 {
        let off = tables.len() as u32;
        for i in 0..ax.inv_kappa_e.len() {
            tables.push([
                R::from_f64(ax.inv_kappa_e[i]),
                R::from_f64(ax.b_e[i]),
                R::from_f64(ax.c_e[i]),
                R::from_f64(if ax.active_e[i] { 1.0 } else { 0.0 }),
            ]);
        }
        off
    };
    let push_cpml_h = |tables: &mut Vec<[R; 4]>, ax: &super::ade_pml::CpmlAxis1d| -> u32 {
        let off = tables.len() as u32;
        for i in 0..ax.inv_kappa_h.len() {
            tables.push([
                R::from_f64(ax.inv_kappa_h[i]),
                R::from_f64(ax.b_h[i]),
                R::from_f64(ax.c_h[i]),
                R::from_f64(if ax.active_h[i] { 1.0 } else { 0.0 }),
            ]);
        }
        off
    };

    let off_inv_dp_x = push_scalar_table(&mut tables, &grid.inv_dp_x);
    let off_inv_dp_y = push_scalar_table(&mut tables, &grid.inv_dp_y);
    let off_inv_dp_z = push_scalar_table(&mut tables, &grid.inv_dp_z);
    let off_inv_dd_x = push_scalar_table(&mut tables, &grid.inv_dd_x);
    let off_inv_dd_y = push_scalar_table(&mut tables, &grid.inv_dd_y);
    let off_inv_dd_z = push_scalar_table(&mut tables, &grid.inv_dd_z);
    let off_cpml_e_x = push_cpml_e(&mut tables, &asm.cpml.x);
    let off_cpml_e_y = push_cpml_e(&mut tables, &asm.cpml.y);
    let off_cpml_e_z = push_cpml_e(&mut tables, &asm.cpml.z);
    let off_cpml_h_x = push_cpml_h(&mut tables, &asm.cpml.x);
    let off_cpml_h_y = push_cpml_h(&mut tables, &asm.cpml.y);
    let off_cpml_h_z = push_cpml_h(&mut tables, &asm.cpml.z);

    // 2. Pack `(ca, cb)` and `(da, db)` coefficients (`3 * slab_len` of `[R; 2]`)
    let mut coeffs_e = Vec::with_capacity(3 * slab_len);
    for (&ca, &cb) in asm.coeffs.ca_x.iter().zip(&asm.coeffs.cb_x) {
        coeffs_e.push([R::from_f64(ca), R::from_f64(cb)]);
    }
    for (&ca, &cb) in asm.coeffs.ca_y.iter().zip(&asm.coeffs.cb_y) {
        coeffs_e.push([R::from_f64(ca), R::from_f64(cb)]);
    }
    for (&ca, &cb) in asm.coeffs.ca_z.iter().zip(&asm.coeffs.cb_z) {
        coeffs_e.push([R::from_f64(ca), R::from_f64(cb)]);
    }

    let mut coeffs_h = Vec::with_capacity(3 * slab_len);
    for (&da, &db) in asm.coeffs.da_x.iter().zip(&asm.coeffs.db_x) {
        coeffs_h.push([R::from_f64(da), R::from_f64(db)]);
    }
    for (&da, &db) in asm.coeffs.da_y.iter().zip(&asm.coeffs.db_y) {
        coeffs_h.push([R::from_f64(da), R::from_f64(db)]);
    }
    for (&da, &db) in asm.coeffs.da_z.iter().zip(&asm.coeffs.db_z) {
        coeffs_h.push([R::from_f64(da), R::from_f64(db)]);
    }

    let use_conf = u32::from(asm.mats.has_conformal_pec);
    let pec_len_host: Vec<R> = if asm.mats.has_conformal_pec {
        let mut v = Vec::with_capacity(3 * slab_len);
        v.extend(asm.mats.pec_len_x.as_slice().iter().map(|&x| R::from_f64(x)));
        v.extend(asm.mats.pec_len_y.as_slice().iter().map(|&x| R::from_f64(x)));
        v.extend(asm.mats.pec_len_z.as_slice().iter().map(|&x| R::from_f64(x)));
        v
    } else {
        vec![R::from_f64(1.0); 4]
    };

    // 3. Compile all H and E source injection sites and pre-evaluate 1D temporal channels across `0..n_steps`
    let (
        src_h_sites,
        src_h_terms,
        src_h_waves,
        n_src_h_chans,
        src_e_sites,
        src_e_terms,
        src_e_waves,
        n_src_e_chans,
        primary_inc_trace,
        est_e_peak,
    ) = compile_gpu_sources::<R>(asm);

    let (e_ref, nc_ref) = if use_f64 { (1.0_f64, 1.0_f64) } else { (est_e_peak.clamp(1.0, 1.0e10), FDTD_NC_REF) };
    let (nl_domains, nl_rank4_flat, nl_edges, nl_has_carriers) = compile_gpu_nl::<R>(asm, use_f64, e_ref, nc_ref);
    let n_nl_edges = nl_edges.len();

    // 4. Compile off-diagonal Farjadpour / anisotropic entries
    let offdiag_entries = compile_gpu_offdiag::<R>(asm);

    // 5. Compile running DFT frequency phasors (`we_re, we_im, wh_re, wh_im`)
    let mut all_freqs: Vec<f64> = Vec::new();
    let mut planar_freq_offsets: Vec<u32> = Vec::with_capacity(asm.planar_dfts.len());
    for pd in &asm.planar_dfts {
        planar_freq_offsets.push(all_freqs.len() as u32);
        all_freqs.extend_from_slice(&pd.freqs_hz);
    }
    let vol_freq_offset = all_freqs.len() as u32;
    if let Some(vd) = asm.vol_dfts.first() {
        all_freqs.extend_from_slice(&vd.freqs_hz);
    }
    let total_dft_freqs = all_freqs.len();
    let mut dft_phasors: Vec<[R; 4]> = Vec::with_capacity(n_steps * total_dft_freqs.max(1));
    if total_dft_freqs > 0 {
        for step in 0..n_steps {
            let t_e = (step as f64 + 1.0) * asm.dt;
            let t_h = (step as f64 + 0.5) * asm.dt;
            for &freq in &all_freqs {
                let w = 2.0 * PI * freq;
                let (sne, cse) = (w * t_e).sin_cos();
                let (snh, csh) = (w * t_h).sin_cos();
                dft_phasors.push([
                    R::from_f64(cse * asm.dt),
                    R::from_f64(-sne * asm.dt),
                    R::from_f64(csh * asm.dt),
                    R::from_f64(-snh * asm.dt),
                ]);
            }
        }
    } else {
        dft_phasors.push([R::default(); 4]);
    }

    // Planar DFT patches
    let mut gpu_planar_patches: Vec<GpuPlanarPatch<R>> = Vec::new();
    let mut total_planar_out_vec4s = 0usize;
    for (pi, pd) in asm.planar_dfts.iter().enumerate() {
        let f_off = planar_freq_offsets[pi];
        let nf = pd.freqs_hz.len() as u32;
        let np = pd.patches.len() as u32;
        let out_off = (total_planar_out_vec4s / 2) as u32;
        for (loc_p, pt) in pd.patches.iter().enumerate() {
            let stencils = [
                to_gpu_stencil::<R>(&pt.stencils[0], pd.u_ax, slab_len),
                to_gpu_stencil::<R>(&pt.stencils[1], pd.v_ax, slab_len),
                to_gpu_stencil::<R>(&pt.stencils[2], pd.u_ax, slab_len),
                to_gpu_stencil::<R>(&pt.stencils[3], pd.v_ax, slab_len),
            ];
            gpu_planar_patches.push(GpuPlanarPatch {
                stencils,
                freq_offset: f_off,
                n_freqs: nf,
                out_offset: out_off,
                patch_in_mon: loc_p as u32,
                n_patches_in_mon: np,
                _pad: [0; 3],
            });
        }
        total_planar_out_vec4s += pd.dft_eh.len() * 2;
    }

    // Volume DFT metadata
    let (vol_i0, vol_i1, vol_j0, vol_j1, vol_k0, vol_k1, n_vol_pts, n_vol_freqs) =
        if let Some(vd) = asm.vol_dfts.first() {
            (
                vd.i_range.0 as u32,
                vd.i_range.1 as u32,
                vd.j_range.0 as u32,
                vd.j_range.1 as u32,
                vd.k_range.0 as u32,
                vd.k_range.1 as u32,
                vd.n_pts() as u32,
                vd.freqs_hz.len() as u32,
            )
        } else {
            (0, 0, 0, 0, 0, 0, 0, 0)
        };
    let total_vol_out_vec4s = (n_vol_pts as usize) * (n_vol_freqs as usize) * 3;

    // Time monitors
    let mut gpu_time_monitors: Vec<GpuTimeMonitor<R>> = Vec::with_capacity(asm.time_monitors.len());
    let mut total_time_samples = 0usize;
    let mut time_mon_counts: Vec<usize> = Vec::with_capacity(asm.time_monitors.len());
    for tm in &asm.time_monitors {
        let is_center = match tm.component {
            0 => [true, false, false],
            1 => [false, true, false],
            2 => [false, false, true],
            3 => [false, true, true],
            4 => [true, false, true],
            _ => [true, true, false],
        };
        let st8 = StaggeredStencil8::build(grid, is_center, tm.position);
        let comp_local = tm.component % 3;
        let is_mag = u32::from(tm.component >= 3);
        let stride = tm.sample_stride.max(1);
        let n_samp = n_steps.div_ceil(stride);
        gpu_time_monitors.push(GpuTimeMonitor {
            stencil: to_gpu_stencil::<R>(&st8, comp_local, slab_len),
            is_magnetic: is_mag,
            sample_stride: stride as u32,
            write_offset: total_time_samples as u32,
            _pad: [0; 5],
        });
        time_mon_counts.push(n_samp);
        total_time_samples += n_samp;
    }

    let n_energy_wgs = (slab_len as u32).div_ceil(256).max(1);
    let has_cpml_x = u32::from(asm.cpml.x.has_cpml());
    let has_cpml_y = u32::from(asm.cpml.y.has_cpml());
    let has_cpml_z = u32::from(asm.cpml.z.has_cpml());
    let any_cpml = (has_cpml_x | has_cpml_y | has_cpml_z) != 0;

    let params = GpuFdtdParams {
        nx: nx as u32,
        ny: ny as u32,
        nz: nz as u32,
        nz_pad: nz_pad as u32,
        stride_x: sx as u32,
        stride_y: sy as u32,
        slab_len: slab_len as u32,
        dim_mode,
        use_conf,
        has_cpml_x,
        has_cpml_y,
        has_cpml_z,
        pmc_mask,
        pec_mask,
        per_mask,
        off_inv_dp_x,
        off_inv_dp_y,
        off_inv_dp_z,
        off_inv_dd_x,
        off_inv_dd_y,
        off_inv_dd_z,
        off_cpml_e_x,
        off_cpml_e_y,
        off_cpml_e_z,
        off_cpml_h_x,
        off_cpml_h_y,
        off_cpml_h_z,
        n_src_h_sites: src_h_sites.len() as u32,
        n_src_h_chans,
        n_src_e_sites: src_e_sites.len() as u32,
        n_src_e_chans,
        n_offdiag: offdiag_entries.len() as u32,
        n_planar_patches: gpu_planar_patches.len() as u32,
        n_vol_pts,
        n_vol_freqs,
        vol_i0,
        vol_i1,
        vol_j0,
        vol_j1,
        vol_k0,
        vol_k1,
        vol_freq_offset,
        total_dft_freqs: total_dft_freqs as u32,
        n_time_monitors: gpu_time_monitors.len() as u32,
        n_energy_wgs,
        n_nl_edges: n_nl_edges as u32,
        nl_has_carriers,
        _pad0: 0,
    };

    // 6. Allocate wgpu buffers
    let dev = &ctx.device;
    let scalar_bytes = std::mem::size_of::<R>();
    let params_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_params"),
        contents: bytemuck::bytes_of(&params),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let step_uniforms: Vec<GpuStepUniform> =
        (0..n_steps.max(1)).map(|s| GpuStepUniform { step: s as u32, _pad: [0; 63] }).collect();
    let step_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_step_uniforms"),
        contents: bytemuck::cast_slice(&step_uniforms),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let fields_bytes = (3 * slab_len * scalar_bytes) as u64;
    let fields_e_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_fields_e"),
        size: fields_bytes.max(32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let fields_h_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_fields_h"),
        size: fields_bytes.max(32),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let coeffs_e_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_coeffs_e"),
        contents: bytemuck::cast_slice(&coeffs_e),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let coeffs_h_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_coeffs_h"),
        contents: bytemuck::cast_slice(&coeffs_h),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let psi_alloc_bytes = if any_cpml { (6 * slab_len * scalar_bytes) as u64 } else { 32 };
    let psi_e_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_psi_e"),
        size: psi_alloc_bytes,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    let psi_h_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_psi_h"),
        size: psi_alloc_bytes,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });

    let tables_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_tables"),
        contents: bytemuck::cast_slice(&tables),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let pec_len_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_pec_len"),
        contents: bytemuck::cast_slice(&pec_len_host),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let dummy_site = [GpuSrcSite::default()];
    let dummy_term = [GpuSrcTerm::<R>::default()];
    let dummy_wave = [[R::default(); 2]];
    let dummy_offdiag = [GpuOffDiagEntry::<R>::default()];
    let dummy_patch = [GpuPlanarPatch::<R>::default()];
    let dummy_tmon = [GpuTimeMonitor::<R>::default()];
    let dummy_nldom = [GpuNlDomain::<R>::default()];
    let dummy_nledge = [GpuNlEdge::<R>::default()];
    let dummy_r4 = [R::default(); 4];

    let src_h_sites_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_h_sites"),
        contents: &non_empty_typed_bytes(&src_h_sites, &dummy_site),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let src_h_terms_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_h_terms"),
        contents: &non_empty_typed_bytes(&src_h_terms, &dummy_term),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let src_h_waves_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_h_waves"),
        contents: &non_empty_typed_bytes(&src_h_waves, &dummy_wave),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let src_e_sites_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_e_sites"),
        contents: &non_empty_typed_bytes(&src_e_sites, &dummy_site),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let src_e_terms_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_e_terms"),
        contents: &non_empty_typed_bytes(&src_e_terms, &dummy_term),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let src_e_waves_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_src_e_waves"),
        contents: &non_empty_typed_bytes(&src_e_waves, &dummy_wave),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let offdiag_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_offdiag_entries"),
        contents: &non_empty_typed_bytes(&offdiag_entries, &dummy_offdiag),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let offdiag_delta_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_offdiag_delta"),
        size: ((offdiag_entries.len().max(1)) * scalar_bytes) as u64,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });

    let dft_phasors_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_dft_phasors"),
        contents: bytemuck::cast_slice(&dft_phasors),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let planar_patches_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_planar_patches"),
        contents: &non_empty_typed_bytes(&gpu_planar_patches, &dummy_patch),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let planar_dft_bytes = (total_planar_out_vec4s.max(1) * 4 * scalar_bytes) as u64;
    let planar_dft_out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_planar_dft_out"),
        size: planar_dft_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let vol_dft_bytes = (total_vol_out_vec4s.max(1) * 4 * scalar_bytes) as u64;
    let vol_dft_out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_vol_dft_out"),
        size: vol_dft_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let time_monitors_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_time_monitors"),
        contents: &non_empty_typed_bytes(&gpu_time_monitors, &dummy_tmon),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let time_out_bytes = (total_time_samples.max(1) * scalar_bytes) as u64;
    let time_out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_time_out"),
        size: time_out_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let energy_partials_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_energy_partials"),
        size: ((n_energy_wgs as usize).max(1) * scalar_bytes) as u64,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    let energy_out_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_energy_out"),
        size: 32,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let energy_readback_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_energy_readback"),
        size: 32,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    // Group 3: Nonlinear polarization & carrier buffers
    let nl_domains_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_nl_domains"),
        contents: &non_empty_typed_bytes(&nl_domains, &dummy_nldom),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let nl_rank4_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_nl_rank4"),
        contents: if nl_rank4_flat.is_empty() {
            bytemuck::cast_slice(&dummy_r4)
        } else {
            bytemuck::cast_slice(&nl_rank4_flat)
        },
        usage: wgpu::BufferUsages::STORAGE,
    });
    let nl_edges_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_nl_edges"),
        contents: &non_empty_typed_bytes(&nl_edges, &dummy_nledge),
        usage: wgpu::BufferUsages::STORAGE,
    });
    let nl_state_buf = dev.create_buffer(&wgpu::BufferDescriptor {
        label: Some("fdtd_nl_state"),
        size: ((5 * n_nl_edges.max(1)) * scalar_bytes).max(32) as u64,
        usage: wgpu::BufferUsages::STORAGE,
        mapped_at_creation: false,
    });
    let nl_run = GpuNlRunParams {
        e_ref: R::from_f64(e_ref),
        inv_e_ref: R::from_f64(1.0 / e_ref),
        _pad0: R::default(),
        _pad1: R::default(),
    };
    let nl_run_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("fdtd_nl_run"),
        contents: bytemuck::bytes_of(&nl_run),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    // 7. Compile WGSL shader and pipelines
    let wgsl_src = build_fdtd_wgsl(use_f64);
    let shader = dev.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(if use_f64 { "fdtd_wgsl_f64" } else { "fdtd_wgsl_f32" }),
        source: wgpu::ShaderSource::Wgsl(wgsl_src.into()),
    });

    let bgl0 = dev.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("fdtd_bgl0"),
        entries: &[
            bgl_uniform(0, false),
            bgl_storage(1, false),
            bgl_storage(2, false),
            bgl_storage(3, true),
            bgl_storage(4, true),
            bgl_storage(5, false),
            bgl_storage(6, true),
            bgl_storage(7, true),
        ],
    });
    let bgl1 = dev.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("fdtd_bgl1"),
        entries: &[
            bgl_uniform(0, true),
            bgl_storage(1, true),
            bgl_storage(2, true),
            bgl_storage(3, true),
            bgl_storage(4, true),
            bgl_storage(5, true),
            bgl_storage(6, true),
            bgl_storage(7, true),
            bgl_storage(8, false),
        ],
    });
    let bgl2 = dev.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("fdtd_bgl2"),
        entries: &[
            bgl_storage(0, true),
            bgl_storage(1, true),
            bgl_storage(2, false),
            bgl_storage(3, false),
            bgl_storage(4, true),
            bgl_storage(5, false),
            bgl_storage(6, false),
            bgl_storage(7, false),
        ],
    });
    let bgl3 = dev.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("fdtd_bgl3"),
        entries: &[
            bgl_storage(0, true),
            bgl_storage(1, true),
            bgl_storage(2, true),
            bgl_storage(3, false),
            bgl_uniform(4, false),
        ],
    });

    let pipeline_layout = dev.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("fdtd_pipeline_layout"),
        bind_group_layouts: &[Some(&bgl0), Some(&bgl1), Some(&bgl2), Some(&bgl3)],
        immediate_size: 0,
    });

    let make_pipe = |entry: &str| -> wgpu::ComputePipeline {
        dev.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some(entry),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        })
    };

    let pipe_update_h = make_pipe("update_h");
    let pipe_inject_h = make_pipe("inject_h");
    let pipe_offdiag_save = make_pipe("offdiag_save");
    let pipe_update_e = make_pipe("update_e");
    let pipe_offdiag_apply = make_pipe("offdiag_apply");
    let pipe_offdiag_commit = make_pipe("offdiag_commit");
    let pipe_inject_e = make_pipe("inject_e");
    let pipe_sync_per_x = make_pipe("sync_periodic_x");
    let pipe_sync_per_y = make_pipe("sync_periodic_y");
    let pipe_sync_per_z = make_pipe("sync_periodic_z");
    let pipe_clamp_pec = make_pipe("clamp_pec_e");
    let pipe_nl_save_old = make_pipe("nl_save_old");
    let pipe_nl_save_lin = make_pipe("nl_save_lin");
    let pipe_nl_newton = make_pipe("nl_newton_iter");
    let pipe_nl_commit = make_pipe("nl_commit_e");
    let pipe_nl_finish = make_pipe("nl_finish_step");
    let pipe_planar_dft = make_pipe("accum_planar_dft");
    let pipe_vol_dft = make_pipe("accum_vol_dft");
    let pipe_time_mon = make_pipe("sample_time_monitors");
    let pipe_energy_p1 = make_pipe("energy_reduce_pass1");
    let pipe_energy_p2 = make_pipe("energy_reduce_pass2");

    let make_bg0 = |label: &str, psi_buf: &wgpu::Buffer| -> wgpu::BindGroup {
        dev.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some(label),
            layout: &bgl0,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: fields_e_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: fields_h_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: coeffs_e_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: coeffs_h_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: psi_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: tables_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: pec_len_buf.as_entire_binding() },
            ],
        })
    };
    let bg0_h = make_bg0("fdtd_bg0_h", &psi_h_buf);
    let bg0_e = make_bg0("fdtd_bg0_e", &psi_e_buf);

    let step_binding_size = std::num::NonZeroU64::new(std::mem::size_of::<GpuStepUniform>() as u64);
    let bg1 = dev.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("fdtd_bg1"),
        layout: &bgl1,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                    buffer: &step_buf,
                    offset: 0,
                    size: step_binding_size,
                }),
            },
            wgpu::BindGroupEntry { binding: 1, resource: src_h_sites_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: src_h_terms_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: src_h_waves_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: src_e_sites_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: src_e_terms_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: src_e_waves_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 7, resource: offdiag_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 8, resource: offdiag_delta_buf.as_entire_binding() },
        ],
    });

    let bg2 = dev.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("fdtd_bg2"),
        layout: &bgl2,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: dft_phasors_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: planar_patches_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: planar_dft_out_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: vol_dft_out_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: time_monitors_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 5, resource: time_out_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 6, resource: energy_partials_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 7, resource: energy_out_buf.as_entire_binding() },
        ],
    });

    let bg3 = dev.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("fdtd_bg3"),
        layout: &bgl3,
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: nl_domains_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: nl_rank4_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: nl_edges_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: nl_state_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 4, resource: nl_run_buf.as_entire_binding() },
        ],
    });

    // 8. Time-marching loop batched in chunks of `check_interval = 40` steps per CommandEncoder submission!
    let wg_slab = (slab_len as u32).div_ceil(256);
    let wg_src_h = (src_h_sites.len() as u32).div_ceil(64);
    let wg_src_e = (src_e_sites.len() as u32).div_ceil(64);
    let wg_offdiag = (offdiag_entries.len() as u32).div_ceil(64);
    let wg_nl = (n_nl_edges as u32).div_ceil(64);
    let wg_planar = (gpu_planar_patches.len() as u32).div_ceil(64);
    let wg_vol = n_vol_pts.div_ceil(64);
    let wg_tmon = (gpu_time_monitors.len() as u32).div_ceil(64);
    let iy_hi = if matches!(grid.dim, FdtdDim::TwoD { .. }) { 0 } else { ny };
    let wg_per_x = (((iy_hi + 1) * (nz + 1)) as u32).div_ceil(64);
    let wg_per_y = (((nx + 1) * (nz + 1)) as u32).div_ceil(64);
    let wg_per_z = (((nx + 1) * (iy_hi + 1)) as u32).div_ceil(64);
    let need_pec_clamp = pec_mask != 0 || use_conf != 0;

    let check_interval = 40usize;
    let min_steps_before_shutoff = ((1.25 * asm.pulse_end_s / asm.dt).ceil() as usize).max(80);
    let mut peak_energy = 0.0_f64;
    let mut final_energy = 0.0_f64;
    let mut n_steps_run = 0usize;
    let mut early_shutoff = false;

    let mut batch_start = 0usize;
    while batch_start < n_steps {
        femcore::progress::check()?;
        let batch_end = (batch_start + check_interval).min(n_steps);
        let mut encoder =
            dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("fdtd_batch_encoder") });

        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("fdtd_batch_cpass"),
                timestamp_writes: None,
            });
            cpass.set_bind_group(2, &bg2, &[]);
            cpass.set_bind_group(3, &bg3, &[]);

            for step in batch_start..batch_end {
                let dyn_off = [(step * std::mem::size_of::<GpuStepUniform>()) as u32];
                cpass.set_bind_group(1, &bg1, &dyn_off);

                // 1. Update H^{n-1/2} -> H^{n+1/2}
                cpass.set_bind_group(0, &bg0_h, &[]);
                cpass.set_pipeline(&pipe_update_h);
                cpass.dispatch_workgroups(wg_slab, 1, 1);

                // 2. Inject H sources at t_n
                if wg_src_h > 0 {
                    cpass.set_pipeline(&pipe_inject_h);
                    cpass.dispatch_workgroups(wg_src_h, 1, 1);
                }

                // 3. Save nonlinear E^n and off-diagonal E^n before Ampere update
                cpass.set_bind_group(0, &bg0_e, &[]);
                if wg_nl > 0 {
                    cpass.set_pipeline(&pipe_nl_save_old);
                    cpass.dispatch_workgroups(wg_nl, 1, 1);
                }
                if wg_offdiag > 0 {
                    cpass.set_pipeline(&pipe_offdiag_save);
                    cpass.dispatch_workgroups(wg_offdiag, 1, 1);
                }

                // 4. Update E^n -> E^{n+1} (including periodic low faces i = 0)
                cpass.set_pipeline(&pipe_update_e);
                cpass.dispatch_workgroups(wg_slab, 1, 1);

                // 5. Apply off-diagonal delta E^{n+1}
                if wg_offdiag > 0 {
                    cpass.set_pipeline(&pipe_offdiag_apply);
                    cpass.dispatch_workgroups(wg_offdiag, 1, 1);
                    cpass.set_pipeline(&pipe_offdiag_commit);
                    cpass.dispatch_workgroups(wg_offdiag, 1, 1);
                }

                // 6. Inject E sources at t_{n+1/2}
                if wg_src_e > 0 {
                    cpass.set_pipeline(&pipe_inject_e);
                    cpass.dispatch_workgroups(wg_src_e, 1, 1);
                }

                // 7. Synchronize periodic high faces (i = n) and clamp PEC boundaries
                if (per_mask & 1) != 0 && wg_per_x > 0 {
                    cpass.set_pipeline(&pipe_sync_per_x);
                    cpass.dispatch_workgroups(wg_per_x, 1, 1);
                }
                if (per_mask & 2) != 0 && wg_per_y > 0 {
                    cpass.set_pipeline(&pipe_sync_per_y);
                    cpass.dispatch_workgroups(wg_per_y, 1, 1);
                }
                if (per_mask & 4) != 0 && wg_per_z > 0 {
                    cpass.set_pipeline(&pipe_sync_per_z);
                    cpass.dispatch_workgroups(wg_per_z, 1, 1);
                }
                if need_pec_clamp {
                    cpass.set_pipeline(&pipe_clamp_pec);
                    cpass.dispatch_workgroups(wg_slab, 1, 1);
                }

                // 7b. Synchronous 3-iteration nonlinear chi^(2)/chi^(3)/TPA/FCA/FCD Newton solve
                if wg_nl > 0 {
                    cpass.set_pipeline(&pipe_nl_save_lin);
                    cpass.dispatch_workgroups(wg_nl, 1, 1);
                    for _ in 0..3 {
                        cpass.set_pipeline(&pipe_nl_newton);
                        cpass.dispatch_workgroups(wg_nl, 1, 1);
                        cpass.set_pipeline(&pipe_nl_commit);
                        cpass.dispatch_workgroups(wg_nl, 1, 1);
                        if (per_mask & 1) != 0 && wg_per_x > 0 {
                            cpass.set_pipeline(&pipe_sync_per_x);
                            cpass.dispatch_workgroups(wg_per_x, 1, 1);
                        }
                        if (per_mask & 2) != 0 && wg_per_y > 0 {
                            cpass.set_pipeline(&pipe_sync_per_y);
                            cpass.dispatch_workgroups(wg_per_y, 1, 1);
                        }
                        if (per_mask & 4) != 0 && wg_per_z > 0 {
                            cpass.set_pipeline(&pipe_sync_per_z);
                            cpass.dispatch_workgroups(wg_per_z, 1, 1);
                        }
                        if need_pec_clamp {
                            cpass.set_pipeline(&pipe_clamp_pec);
                            cpass.dispatch_workgroups(wg_slab, 1, 1);
                        }
                    }
                    cpass.set_pipeline(&pipe_nl_finish);
                    cpass.dispatch_workgroups(wg_nl, 1, 1);
                }

                // 8. Accumulate planar & volumetric running DFTs and sample time probes
                if wg_planar > 0 {
                    cpass.set_pipeline(&pipe_planar_dft);
                    cpass.dispatch_workgroups(wg_planar, 1, 1);
                }
                if wg_vol > 0 && n_vol_freqs > 0 {
                    cpass.set_pipeline(&pipe_vol_dft);
                    cpass.dispatch_workgroups(wg_vol, 1, 1);
                }
                if wg_tmon > 0 {
                    cpass.set_pipeline(&pipe_time_mon);
                    cpass.dispatch_workgroups(wg_tmon, 1, 1);
                }
            }

            // 9. Fixed-order two-level tree reduction of total field energy proxy at end of batch
            cpass.set_bind_group(0, &bg0_e, &[]);
            cpass.set_pipeline(&pipe_energy_p1);
            cpass.dispatch_workgroups(n_energy_wgs, 1, 1);
            cpass.set_pipeline(&pipe_energy_p2);
            cpass.dispatch_workgroups(1, 1, 1);
        }

        encoder.copy_buffer_to_buffer(&energy_out_buf, 0, &energy_readback_buf, 0, scalar_bytes as u64);
        ctx.queue.submit(Some(encoder.finish()));

        let e_slice: Vec<R> = read_back_typed_slice(dev, &energy_readback_buf, 1)?;
        let e_now = e_slice[0].to_f64();
        n_steps_run = batch_end;
        if !e_now.is_finite() {
            return Err(format!("FDTD GPU instability detected at step {n_steps_run}: non-finite field energy"));
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
        batch_start = batch_end;
    }

    // 9. Read back ONLY monitor spectra (`planar_dft_out`, `vol_dft_out`, `time_out`)!
    let mut planar_dfts = asm.planar_dfts.clone();
    if total_planar_out_vec4s > 0 {
        let raw: Vec<R> = copy_and_read_typed(dev, &ctx.queue, &planar_dft_out_buf, total_planar_out_vec4s * 4)?;
        let mut cur = 0usize;
        for pd in &mut planar_dfts {
            for eh in &mut pd.dft_eh {
                eh[0] = c(raw[cur].to_f64(), raw[cur + 1].to_f64());
                eh[1] = c(raw[cur + 2].to_f64(), raw[cur + 3].to_f64());
                eh[2] = c(raw[cur + 4].to_f64(), raw[cur + 5].to_f64());
                eh[3] = c(raw[cur + 6].to_f64(), raw[cur + 7].to_f64());
                cur += 8;
            }
        }
    }

    let mut vol_dfts = asm.vol_dfts.clone();
    if total_vol_out_vec4s > 0 {
        let raw: Vec<R> = copy_and_read_typed(dev, &ctx.queue, &vol_dft_out_buf, total_vol_out_vec4s * 4)?;
        if let Some(vd) = vol_dfts.first_mut() {
            let mut cur = 0usize;
            for eh in &mut vd.eh_dft {
                eh[0] = c(raw[cur].to_f64(), raw[cur + 1].to_f64());
                eh[1] = c(raw[cur + 2].to_f64(), raw[cur + 3].to_f64());
                eh[2] = c(raw[cur + 4].to_f64(), raw[cur + 5].to_f64());
                eh[3] = c(raw[cur + 6].to_f64(), raw[cur + 7].to_f64());
                eh[4] = c(raw[cur + 8].to_f64(), raw[cur + 9].to_f64());
                eh[5] = c(raw[cur + 10].to_f64(), raw[cur + 11].to_f64());
                cur += 12;
            }
        }
    }

    let mut time_traces = Vec::with_capacity(asm.time_monitors.len());
    if total_time_samples > 0 {
        let raw_time: Vec<R> = copy_and_read_typed(dev, &ctx.queue, &time_out_buf, total_time_samples)?;
        for (mi, tm) in asm.time_monitors.iter().enumerate() {
            let stride = tm.sample_stride.max(1);
            let w_off = gpu_time_monitors[mi].write_offset as usize;
            let max_s = time_mon_counts[mi];
            let mut ts = Vec::with_capacity(max_s);
            let mut vs = Vec::with_capacity(max_s);
            for step in (0..n_steps_run).step_by(stride) {
                let s_idx = step / stride;
                if s_idx < max_s {
                    let t_now = ((step + 1) as f64) * asm.dt;
                    let t_probe = if tm.component < 3 { t_now } else { t_now - 0.5 * asm.dt };
                    ts.push(t_probe);
                    vs.push(raw_time[w_off + s_idx].to_f64());
                }
            }
            let resonances = if let Some(frange) = tm.prony_freq_range_hz {
                let t_start = tm.prony_start_time_s.unwrap_or(1.15 * asm.pulse_end_s);
                let start_idx = ts.iter().position(|&t| t >= t_start).unwrap_or(0);
                let dt_samp = asm.dt * (stride as f64);
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
    }

    let primary_inc = primary_inc_trace.map(|mut tr| {
        tr.truncate(n_steps_run);
        tr
    });

    Ok(FdtdSolve {
        n_steps_run,
        early_shutoff,
        peak_energy,
        final_energy,
        fields_re: YeeFields::new(grid),
        fields_im: None,
        planar_dfts,
        vol_dfts,
        primary_inc_trace: primary_inc,
        time_traces,
        solve_seconds: t0.elapsed().as_secs_f64(),
        backend_used: "gpu",
        gpu_adapter: Some(format!("{}, precision={}", ctx.adapter_summary(), if use_f64 { "f64" } else { "f32" })),
        fallback_warning: None,
    })
}

fn bgl_uniform(binding: u32, has_dynamic_offset: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset,
            min_binding_size: None,
        },
        count: None,
    }
}

fn bgl_storage(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn to_gpu_stencil<R: GpuScalar>(st: &StaggeredStencil8, comp: usize, slab_len: usize) -> GpuStencil8<R> {
    let base = (comp * slab_len) as u32;
    let mut idx = [0u32; 8];
    let mut w0 = [R::default(); 4];
    let mut w1 = [R::default(); 4];
    for k in 0..8 {
        idx[k] = base + (st.idx[k] as u32);
        if k < 4 {
            w0[k] = R::from_f64(st.w[k]);
        } else {
            w1[k - 4] = R::from_f64(st.w[k]);
        }
    }
    GpuStencil8 { idx, w0, w1 }
}

pub(crate) fn copy_and_read_typed<T: Pod>(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    src_buf: &wgpu::Buffer,
    count: usize,
) -> Result<Vec<T>, String> {
    let byte_len = (count * std::mem::size_of::<T>()) as u64;
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("ems_gpu_readback_staging"),
        size: byte_len.max(8),
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder =
        device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("ems_gpu_readback_encoder") });
    encoder.copy_buffer_to_buffer(src_buf, 0, &staging, 0, byte_len);
    queue.submit(Some(encoder.finish()));
    read_back_typed_slice(device, &staging, count)
}

pub(crate) fn read_back_typed_slice<T: Pod>(
    device: &wgpu::Device,
    staging: &wgpu::Buffer,
    count: usize,
) -> Result<Vec<T>, String> {
    let byte_len = (count * std::mem::size_of::<T>()) as u64;
    let slice = staging.slice(0..byte_len);
    let status: Arc<Mutex<Option<Result<(), String>>>> = Arc::new(Mutex::new(None));
    let status_cb = Arc::clone(&status);
    slice.map_async(wgpu::MapMode::Read, move |res| {
        if let Ok(mut guard) = status_cb.lock() {
            *guard = Some(res.map_err(|e| e.to_string()));
        }
    });
    device.poll(wgpu::PollType::wait_indefinitely()).map_err(|e| format!("wgpu device poll failed: {e}"))?;

    let map_res = status
        .lock()
        .map_err(|_| "wgpu readback mutex poisoned".to_string())?
        .take()
        .ok_or_else(|| "wgpu map_async callback did not complete".to_string())?;
    map_res?;

    let view = slice.get_mapped_range().map_err(|e| format!("wgpu get_mapped_range failed: {e}"))?;
    let elems: Vec<T> = bytemuck::cast_slice(&view).to_vec();
    drop(view);
    staging.unmap();
    Ok(elems)
}

/// Compiles all FDTD sources into deduplicated `GpuSrcSite` lists and pre-evaluates temporal waveforms
/// across `0..asm.n_steps` (including stepping `Aux1dYeeLine` for `PlaneSheet` and `TfsfBox`).
#[allow(clippy::type_complexity)]
fn compile_gpu_sources<R: GpuScalar>(
    asm: &FdtdAssembly,
) -> (
    Vec<GpuSrcSite>,
    Vec<GpuSrcTerm<R>>,
    Vec<[R; 2]>,
    u32,
    Vec<GpuSrcSite>,
    Vec<GpuSrcTerm<R>>,
    Vec<[R; 2]>,
    u32,
    Option<Vec<f64>>,
    f64,
) {
    let grid = &asm.grid;
    let coeffs = &asm.coeffs;
    let slab_len = grid.slab_len();
    let n_steps = asm.n_steps;
    let dt = asm.dt;
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });

    let stride_ax = |ax: usize| match ax {
        0 => grid.stride_x,
        1 => grid.stride_y,
        _ => 1,
    };
    let db_comp = |comp: usize| -> &[f64] {
        match comp {
            0 => &coeffs.db_x,
            1 => &coeffs.db_y,
            _ => &coeffs.db_z,
        }
    };
    let cb_comp = |comp: usize| -> &[f64] {
        match comp {
            0 => &coeffs.cb_x,
            1 => &coeffs.cb_y,
            _ => &coeffs.cb_z,
        }
    };

    let mut h_map: BTreeMap<u32, BTreeMap<u32, (f64, f64)>> = BTreeMap::new();
    let mut e_map: BTreeMap<u32, BTreeMap<u32, (f64, f64)>> = BTreeMap::new();
    let mut h_chan_waves: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut e_chan_waves: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut primary_inc_trace: Option<Vec<f64>> = None;
    let mut est_e_peak = 0.0_f64;

    let mut add_h = |target: u32, chan: u32, wre: f64, wim: f64| {
        if wre.abs() > 0.0 || wim.abs() > 0.0 {
            let entry = h_map.entry(target).or_default().entry(chan).or_insert((0.0, 0.0));
            entry.0 += wre;
            entry.1 += wim;
        }
    };
    let mut add_e = |target: u32, chan: u32, wre: f64, wim: f64| {
        if wre.abs() > 0.0 || wim.abs() > 0.0 {
            let entry = e_map.entry(target).or_default().entry(chan).or_insert((0.0, 0.0));
            entry.0 += wre;
            entry.1 += wim;
        }
    };

    for (si, src) in asm.sources.iter().enumerate() {
        match src {
            CompiledSource::Mode { profile, waveform, h_time_offset_s, inv_dp_h, inv_dd_e, .. }
            | CompiledSource::GaussianSheet { profile, waveform, h_time_offset_s, inv_dp_h, inv_dd_e, .. } => {
                let h_chan = h_chan_waves.len() as u32;
                let mut hw = Vec::with_capacity(n_steps);
                let e_chan = e_chan_waves.len() as u32;
                let mut ew = Vec::with_capacity(n_steps);
                let mut wave_max = 0.0_f64;
                for step in 0..n_steps {
                    let t_n = (step as f64) * dt;
                    let t_half = (step as f64 + 0.5) * dt;
                    let (sh_re, sh_im) = waveform.eval_analytic(t_n);
                    hw.push([sh_re, sh_im]);
                    let (se_re, se_im) = waveform.eval_analytic(t_half + *h_time_offset_s);
                    wave_max = wave_max.max(se_re.hypot(se_im));
                    ew.push([se_re, se_im]);
                }
                let e_prof_max = profile
                    .e_u_samples
                    .iter()
                    .chain(profile.e_v_samples.iter())
                    .map(|(_, z)| z.norm())
                    .fold(0.0_f64, f64::max);
                est_e_peak = est_e_peak.max(wave_max * e_prof_max);

                h_chan_waves.push(hw);
                e_chan_waves.push(ew);

                let st_a = stride_ax(profile.axis) as isize;
                let shift_h = (profile.h_cell_idx as isize - profile.plane_idx as isize) * st_a;
                let shift_e = -shift_h;
                let sg = profile.sign;
                let (u_ax, v_ax) = (profile.u_axis, profile.v_axis);

                let dbu = db_comp(u_ax);
                for &(g_ev, ev_c) in &profile.e_v_samples {
                    let g_hu = (g_ev as isize + shift_h) as usize;
                    let scale = -sg * dbu[g_hu] * *inv_dp_h;
                    add_h((u_ax * slab_len + g_hu) as u32, h_chan, scale * ev_c.re, -scale * ev_c.im);
                }
                let dbv = db_comp(v_ax);
                for &(g_eu, eu_c) in &profile.e_u_samples {
                    let g_hv = (g_eu as isize + shift_h) as usize;
                    let scale = sg * dbv[g_hv] * *inv_dp_h;
                    add_h((v_ax * slab_len + g_hv) as u32, h_chan, scale * eu_c.re, -scale * eu_c.im);
                }

                let cbu = cb_comp(u_ax);
                for &(g_hv, hv_c) in &profile.h_v_samples {
                    let g_eu = (g_hv as isize + shift_e) as usize;
                    let scale = sg * cbu[g_eu] * *inv_dd_e;
                    add_e((u_ax * slab_len + g_eu) as u32, e_chan, scale * hv_c.re, -scale * hv_c.im);
                }
                let cbv = cb_comp(v_ax);
                for &(g_hu, hu_c) in &profile.h_u_samples {
                    let g_ev = (g_hu as isize + shift_e) as usize;
                    let scale = -sg * cbv[g_ev] * *inv_dd_e;
                    add_e((v_ax * slab_len + g_ev) as u32, e_chan, scale * hu_c.re, -scale * hu_c.im);
                }
            }
            CompiledSource::PlaneSheet {
                axis, dir_sign, plane_idx, h_cell_idx, cos_pol, sin_pol, eta_bg, aux, ..
            } => {
                let mut aux_sim: Aux1dYeeLine = aux.clone();
                let h_chan = h_chan_waves.len() as u32;
                let e_chan = e_chan_waves.len() as u32;
                let mut hw = Vec::with_capacity(n_steps);
                let mut ew = Vec::with_capacity(n_steps);
                for step in 0..n_steps {
                    let t_n = (step as f64) * dt;
                    let t_half = (step as f64 + 0.5) * dt;
                    let e_1d = aux_sim.e_at_vertex(*plane_idx);
                    est_e_peak = est_e_peak.max(e_1d.abs());
                    hw.push([e_1d, 0.0]);
                    aux_sim.step_h(t_n, *eta_bg);
                    let h_1d = aux_sim.h_at_cell(*h_cell_idx);
                    ew.push([h_1d, 0.0]);
                    aux_sim.step_e(t_half, *eta_bg, dt);
                }
                if si == 0 && !aux_sim.inc_e_trace.is_empty() {
                    primary_inc_trace = Some(aux_sim.inc_e_trace);
                }
                h_chan_waves.push(hw);
                e_chan_waves.push(ew);

                let (u_ax, v_ax) = rh_uv_axes(*axis);
                let inv_dp = grid.axis_dp(*axis)[*h_cell_idx].recip();
                let inv_dd = grid.axis_dd(*axis)[*plane_idx].recip();
                let sg = *dir_sign;
                let nu = grid.axis_centers(u_ax).len();
                let nv = grid.axis_centers(v_ax).len();
                let iu_hi = if is_2d && u_ax == 1 { 0 } else { nu };
                let iv_hi = if is_2d && v_ax == 1 { 0 } else { nv };
                let to_idx = |ia: usize, iu: usize, iv: usize| {
                    let mut ijk = [0usize; 3];
                    ijk[*axis] = ia;
                    ijk[u_ax] = iu;
                    ijk[v_ax] = iv;
                    grid.idx(ijk[0], ijk[1], ijk[2])
                };

                if sin_pol.abs() > 0.0 {
                    let dbu = db_comp(u_ax);
                    for iu in 0..=iu_hi {
                        for iv in 0..nv {
                            let g = to_idx(*h_cell_idx, iu, iv);
                            add_h((u_ax * slab_len + g) as u32, h_chan, -sg * dbu[g] * *sin_pol * inv_dp, 0.0);
                        }
                    }
                    let cbv = cb_comp(v_ax);
                    for iu in 0..=iu_hi {
                        for iv in 0..nv {
                            let g = to_idx(*plane_idx, iu, iv);
                            add_e((v_ax * slab_len + g) as u32, e_chan, sg * cbv[g] * *sin_pol * inv_dd, 0.0);
                        }
                    }
                }
                if cos_pol.abs() > 0.0 {
                    let dbv = db_comp(v_ax);
                    for iu in 0..nu {
                        for iv in 0..=iv_hi {
                            let g = to_idx(*h_cell_idx, iu, iv);
                            add_h((v_ax * slab_len + g) as u32, h_chan, sg * dbv[g] * *cos_pol * inv_dp, 0.0);
                        }
                    }
                    let cbu = cb_comp(u_ax);
                    for iu in 0..nu {
                        for iv in 0..=iv_hi {
                            let g = to_idx(*plane_idx, iu, iv);
                            add_e((u_ax * slab_len + g) as u32, e_chan, sg * cbu[g] * *cos_pol * inv_dd, 0.0);
                        }
                    }
                }
            }
            CompiledSource::TfsfBox { i_lo, i_hi, axis, cos_pol, sin_pol, eta_bg, aux, .. } => {
                let prop_axis = *axis;
                let (pu, pv) = rh_uv_axes(prop_axis);
                let v_min = i_lo[prop_axis].saturating_sub(1);
                let v_max = i_hi[prop_axis] + 1;
                let n_v = v_max - v_min + 1;
                let h_chan_base = h_chan_waves.len() as u32;
                let e_chan_base = e_chan_waves.len() as u32;
                let mut h_vecs = vec![Vec::with_capacity(n_steps); n_v];
                let mut e_vecs = vec![Vec::with_capacity(n_steps); n_v];
                let mut aux_sim = aux.clone();

                for step in 0..n_steps {
                    let t_n = (step as f64) * dt;
                    let t_half = (step as f64 + 0.5) * dt;
                    for k in 0..n_v {
                        let e_1d = aux_sim.e_at_vertex(v_min + k);
                        est_e_peak = est_e_peak.max(e_1d.abs());
                        h_vecs[k].push([e_1d, 0.0]);
                    }
                    aux_sim.step_h(t_n, *eta_bg);
                    for k in 0..n_v {
                        e_vecs[k].push([aux_sim.h_at_cell(v_min + k), 0.0]);
                    }
                    aux_sim.step_e(t_half, *eta_bg, dt);
                }
                if si == 0 && !aux_sim.inc_e_trace.is_empty() {
                    primary_inc_trace = Some(aux_sim.inc_e_trace);
                }
                h_chan_waves.extend(h_vecs);
                e_chan_waves.extend(e_vecs);

                let pol_e = |comp: usize| -> f64 {
                    if comp == pu {
                        *cos_pol
                    } else if comp == pv {
                        *sin_pol
                    } else {
                        0.0
                    }
                };
                let pol_h = |comp: usize| -> f64 {
                    if comp == pu {
                        -*sin_pol
                    } else if comp == pv {
                        *cos_pol
                    } else {
                        0.0
                    }
                };

                for f_ax in 0..3 {
                    if is_2d && f_ax == 1 {
                        continue;
                    }
                    let (fu, fv) = rh_uv_axes(f_ax);
                    let fu_v_hi = if is_2d && fu == 1 { 0 } else { i_hi[fu] };
                    let fv_v_hi = if is_2d && fv == 1 { 0 } else { i_hi[fv] };
                    for &(is_hi, sg) in &[(false, 1.0_f64), (true, -1.0_f64)] {
                        let ia_e = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] };
                        let ia_h = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] - 1 };
                        let inv_dp = grid.axis_dp(f_ax)[ia_h].recip();
                        let inv_dd = grid.axis_dd(f_ax)[ia_e].recip();

                        let p_efv = pol_e(fv);
                        if p_efv.abs() > 0.0 {
                            let dbu = db_comp(fu);
                            for iu in i_lo[fu]..=fu_v_hi {
                                for iv in i_lo[fv]..i_hi[fv] {
                                    let mut ijk_e = [0usize; 3];
                                    ijk_e[f_ax] = ia_e;
                                    ijk_e[fu] = iu;
                                    ijk_e[fv] = iv;
                                    let mut ijk_h = ijk_e;
                                    ijk_h[f_ax] = ia_h;
                                    let gh = grid.idx(ijk_h[0], ijk_h[1], ijk_h[2]);
                                    let ch = h_chan_base + (ijk_e[prop_axis] - v_min) as u32;
                                    add_h((fu * slab_len + gh) as u32, ch, -sg * dbu[gh] * p_efv * inv_dp, 0.0);
                                }
                            }
                        }
                        let p_efu = pol_e(fu);
                        if p_efu.abs() > 0.0 {
                            let dbv = db_comp(fv);
                            for iu in i_lo[fu]..i_hi[fu] {
                                for iv in i_lo[fv]..=fv_v_hi {
                                    let mut ijk_e = [0usize; 3];
                                    ijk_e[f_ax] = ia_e;
                                    ijk_e[fu] = iu;
                                    ijk_e[fv] = iv;
                                    let mut ijk_h = ijk_e;
                                    ijk_h[f_ax] = ia_h;
                                    let gh = grid.idx(ijk_h[0], ijk_h[1], ijk_h[2]);
                                    let ch = h_chan_base + (ijk_e[prop_axis] - v_min) as u32;
                                    add_h((fv * slab_len + gh) as u32, ch, sg * dbv[gh] * p_efu * inv_dp, 0.0);
                                }
                            }
                        }
                        let p_hfv = pol_h(fv);
                        if p_hfv.abs() > 0.0 {
                            let cbu = cb_comp(fu);
                            for iu in i_lo[fu]..i_hi[fu] {
                                for iv in i_lo[fv]..=fv_v_hi {
                                    let mut ijk_h = [0usize; 3];
                                    ijk_h[f_ax] = ia_h;
                                    ijk_h[fu] = iu;
                                    ijk_h[fv] = iv;
                                    let mut ijk_e = ijk_h;
                                    ijk_e[f_ax] = ia_e;
                                    let ge = grid.idx(ijk_e[0], ijk_e[1], ijk_e[2]);
                                    let ch = e_chan_base + (ijk_h[prop_axis] - v_min) as u32;
                                    add_e((fu * slab_len + ge) as u32, ch, sg * cbu[ge] * p_hfv * inv_dd, 0.0);
                                }
                            }
                        }
                        let p_hfu = pol_h(fu);
                        if p_hfu.abs() > 0.0 {
                            let cbv = cb_comp(fv);
                            for iu in i_lo[fu]..=fu_v_hi {
                                for iv in i_lo[fv]..i_hi[fv] {
                                    let mut ijk_h = [0usize; 3];
                                    ijk_h[f_ax] = ia_h;
                                    ijk_h[fu] = iu;
                                    ijk_h[fv] = iv;
                                    let mut ijk_e = ijk_h;
                                    ijk_e[f_ax] = ia_e;
                                    let ge = grid.idx(ijk_e[0], ijk_e[1], ijk_e[2]);
                                    let ch = e_chan_base + (ijk_h[prop_axis] - v_min) as u32;
                                    add_e((fv * slab_len + ge) as u32, ch, -sg * cbv[ge] * p_hfu * inv_dd, 0.0);
                                }
                            }
                        }
                    }
                }
            }
            CompiledSource::Dipole { component, magnetic, weights, waveform, .. } => {
                if *magnetic {
                    let h_chan = h_chan_waves.len() as u32;
                    let mut hw = Vec::with_capacity(n_steps);
                    for step in 0..n_steps {
                        let t_n = (step as f64) * dt;
                        hw.push([waveform.eval(t_n), 0.0]);
                    }
                    h_chan_waves.push(hw);
                    let db = db_comp(*component);
                    for &(g, w) in weights {
                        add_h((*component * slab_len + g) as u32, h_chan, -db[g] * w, 0.0);
                    }
                } else {
                    let e_chan = e_chan_waves.len() as u32;
                    let mut ew = Vec::with_capacity(n_steps);
                    let mut wmax = 0.0_f64;
                    for step in 0..n_steps {
                        let t_half = (step as f64 + 0.5) * dt;
                        let s = waveform.eval(t_half);
                        wmax = wmax.max(s.abs());
                        ew.push([s, 0.0]);
                    }
                    e_chan_waves.push(ew);
                    let cb = cb_comp(*component);
                    for &(g, w) in weights {
                        est_e_peak = est_e_peak.max((cb[g] * w * wmax).abs());
                        add_e((*component * slab_len + g) as u32, e_chan, -cb[g] * w, 0.0);
                    }
                }
            }
        }
    }

    let flatten_sites = |map: BTreeMap<u32, BTreeMap<u32, (f64, f64)>>,
                         chan_waves: &[Vec<[f64; 2]>]|
     -> (Vec<GpuSrcSite>, Vec<GpuSrcTerm<R>>, Vec<[R; 2]>, u32) {
        let n_chans = chan_waves.len() as u32;
        let mut sites = Vec::with_capacity(map.len());
        let mut terms = Vec::new();
        for (target_offset, ch_map) in map {
            let item_start = terms.len() as u32;
            for (chan, (wre, wim)) in ch_map {
                terms.push(GpuSrcTerm { w_re: R::from_f64(wre), w_im: R::from_f64(wim), chan, _pad: 0 });
            }
            let item_end = terms.len() as u32;
            sites.push(GpuSrcSite { target_offset, item_start, item_end, _pad: 0 });
        }
        let mut waves = Vec::with_capacity(n_steps * (n_chans as usize).max(1));
        if n_chans > 0 {
            for step in 0..n_steps {
                for ch in 0..(n_chans as usize) {
                    waves.push([R::from_f64(chan_waves[ch][step][0]), R::from_f64(chan_waves[ch][step][1])]);
                }
            }
        } else {
            waves.push([R::default(), R::default()]);
        }
        (sites, terms, waves, n_chans)
    };

    let (h_sites, h_terms, h_waves, n_h_chans) = flatten_sites(h_map, &h_chan_waves);
    let (e_sites, e_terms, e_waves, n_e_chans) = flatten_sites(e_map, &e_chan_waves);

    (h_sites, h_terms, h_waves, n_h_chans, e_sites, e_terms, e_waves, n_e_chans, primary_inc_trace, est_e_peak)
}

fn compile_gpu_offdiag<R: GpuScalar>(asm: &FdtdAssembly) -> Vec<GpuOffDiagEntry<R>> {
    let grid = &asm.grid;
    let mats = &asm.mats;
    if !mats.has_offdiag() {
        return Vec::new();
    }
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let slab_len = grid.slab_len();
    let iex = mats.inv_eps_x.as_slice();
    let iey = mats.inv_eps_y.as_slice();
    let iez = mats.inv_eps_z.as_slice();

    let inv_iex = |g: usize| if iex[g] > 1e-12 { 1.0 / iex[g] } else { 0.0 };
    let inv_iey = |g: usize| if iey[g] > 1e-12 { 1.0 / iey[g] } else { 0.0 };
    let inv_iez = |g: usize| if iez[g] > 1e-12 { 1.0 / iez[g] } else { 0.0 };

    let mut out = Vec::with_capacity(mats.offdiag_x.len() + mats.offdiag_y.len() + mats.offdiag_z.len());

    // Ex offdiag: b = Ey (comp 1), c = Ez (comp 2)
    for e in &mats.offdiag_x {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let (b_g, b_scale) = if is_2d {
            let g0 = grid.idx(ix, 0, iz);
            let g1 = grid.idx(ix + 1, 0, iz);
            ([g0, g1, g0, g1], [0.5, 0.5, 0.0, 0.0])
        } else {
            let ym = iy.saturating_sub(1);
            let yp = iy.min(grid.ny - 1);
            (
                [grid.idx(ix, ym, iz), grid.idx(ix + 1, ym, iz), grid.idx(ix, yp, iz), grid.idx(ix + 1, yp, iz)],
                [0.25; 4],
            )
        };
        let zm = iz.saturating_sub(1);
        let zp = iz.min(grid.nz - 1);
        let c_g = [grid.idx(ix, iy, zm), grid.idx(ix + 1, iy, zm), grid.idx(ix, iy, zp), grid.idx(ix + 1, iy, zp)];
        out.push(GpuOffDiagEntry {
            b_weights: std::array::from_fn(|k| R::from_f64(e.w_ab * b_scale[k] * inv_iey(b_g[k]))),
            c_weights: std::array::from_fn(|k| R::from_f64(e.w_ac * 0.25 * inv_iez(c_g[k]))),
            b_offsets: b_g.map(|g| (slab_len + g) as u32),
            c_offsets: c_g.map(|g| (2 * slab_len + g) as u32),
            target_offset: e.idx as u32,
            _pad: [0; 7],
        });
    }

    // Ey offdiag: b = Ex (comp 0), c = Ez (comp 2)
    for e in &mats.offdiag_y {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let xm = ix.saturating_sub(1);
        let xp = ix.min(grid.nx - 1);
        let zm = iz.saturating_sub(1);
        let zp = iz.min(grid.nz - 1);
        let (b_g, b_scale) = if is_2d {
            let g0 = grid.idx(xm, 0, iz);
            let g1 = grid.idx(xp, 0, iz);
            ([g0, g1, g0, g1], [0.5, 0.5, 0.0, 0.0])
        } else {
            (
                [grid.idx(xm, iy, iz), grid.idx(xp, iy, iz), grid.idx(xm, iy + 1, iz), grid.idx(xp, iy + 1, iz)],
                [0.25; 4],
            )
        };
        let (c_g, c_scale) = if is_2d {
            let g0 = grid.idx(ix, 0, zm);
            let g1 = grid.idx(ix, 0, zp);
            ([g0, g1, g0, g1], [0.5, 0.5, 0.0, 0.0])
        } else {
            (
                [grid.idx(ix, iy, zm), grid.idx(ix, iy, zp), grid.idx(ix, iy + 1, zm), grid.idx(ix, iy + 1, zp)],
                [0.25; 4],
            )
        };
        out.push(GpuOffDiagEntry {
            b_weights: std::array::from_fn(|k| R::from_f64(e.w_ab * b_scale[k] * inv_iex(b_g[k]))),
            c_weights: std::array::from_fn(|k| R::from_f64(e.w_ac * c_scale[k] * inv_iez(c_g[k]))),
            b_offsets: b_g.map(|g| g as u32),
            c_offsets: c_g.map(|g| (2 * slab_len + g) as u32),
            target_offset: (slab_len + e.idx) as u32,
            _pad: [0; 7],
        });
    }

    // Ez offdiag: b = Ex (comp 0), c = Ey (comp 1)
    for e in &mats.offdiag_z {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let xm = ix.saturating_sub(1);
        let xp = ix.min(grid.nx - 1);
        let b_g = [grid.idx(xm, iy, iz), grid.idx(xp, iy, iz), grid.idx(xm, iy, iz + 1), grid.idx(xp, iy, iz + 1)];
        let (c_g, c_scale) = if is_2d {
            let g0 = grid.idx(ix, 0, iz);
            let g1 = grid.idx(ix, 0, iz + 1);
            ([g0, g1, g0, g1], [0.5, 0.5, 0.0, 0.0])
        } else {
            let ym = iy.saturating_sub(1);
            let yp = iy.min(grid.ny - 1);
            (
                [grid.idx(ix, ym, iz), grid.idx(ix, yp, iz), grid.idx(ix, ym, iz + 1), grid.idx(ix, yp, iz + 1)],
                [0.25; 4],
            )
        };
        out.push(GpuOffDiagEntry {
            b_weights: std::array::from_fn(|k| R::from_f64(e.w_ab * 0.25 * inv_iex(b_g[k]))),
            c_weights: std::array::from_fn(|k| R::from_f64(e.w_ac * c_scale[k] * inv_iey(c_g[k]))),
            b_offsets: b_g.map(|g| g as u32),
            c_offsets: c_g.map(|g| (slab_len + g) as u32),
            target_offset: (2 * slab_len + e.idx) as u32,
            _pad: [0; 7],
        });
    }

    out
}

fn build_fdtd_wgsl(use_f64: bool) -> String {
    let rt = if use_f64 { "f64" } else { "f32" };
    FDTD_WGSL_TEMPLATE.replace("__R__", rt)
}

const FDTD_WGSL_TEMPLATE: &str = r#"
struct GridParams {
    nx: u32,
    ny: u32,
    nz: u32,
    nz_pad: u32,
    stride_x: u32,
    stride_y: u32,
    slab_len: u32,
    dim_mode: u32,
    use_conf: u32,
    has_cpml_x: u32,
    has_cpml_y: u32,
    has_cpml_z: u32,
    pmc_mask: u32,
    pec_mask: u32,
    per_mask: u32,
    off_inv_dp_x: u32,
    off_inv_dp_y: u32,
    off_inv_dp_z: u32,
    off_inv_dd_x: u32,
    off_inv_dd_y: u32,
    off_inv_dd_z: u32,
    off_cpml_e_x: u32,
    off_cpml_e_y: u32,
    off_cpml_e_z: u32,
    off_cpml_h_x: u32,
    off_cpml_h_y: u32,
    off_cpml_h_z: u32,
    n_src_h_sites: u32,
    n_src_h_chans: u32,
    n_src_e_sites: u32,
    n_src_e_chans: u32,
    n_offdiag: u32,
    n_planar_patches: u32,
    n_vol_pts: u32,
    n_vol_freqs: u32,
    vol_i0: u32,
    vol_i1: u32,
    vol_j0: u32,
    vol_j1: u32,
    vol_k0: u32,
    vol_k1: u32,
    vol_freq_offset: u32,
    total_dft_freqs: u32,
    n_time_monitors: u32,
    n_energy_wgs: u32,
    n_nl_edges: u32,
    nl_has_carriers: u32,
    _pad0: u32,
}

struct StepUniform {
    step: u32,
}

struct SrcSite {
    target_offset: u32,
    item_start: u32,
    item_end: u32,
    _pad: u32,
}

struct SrcTerm {
    w_re: __R__,
    w_im: __R__,
    chan: u32,
    _pad: u32,
}

struct OffDiagEntry {
    b_weights: vec4<__R__>,
    c_weights: vec4<__R__>,
    b_offsets: vec4<u32>,
    c_offsets: vec4<u32>,
    target_offset: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    _pad3: u32,
    _pad4: u32,
    _pad5: u32,
    _pad6: u32,
}

struct Stencil8 {
    idx0: vec4<u32>,
    idx1: vec4<u32>,
    w0: vec4<__R__>,
    w1: vec4<__R__>,
}

struct PlanarPatch {
    st_eu: Stencil8,
    st_ev: Stencil8,
    st_hu: Stencil8,
    st_hv: Stencil8,
    freq_offset: u32,
    n_freqs: u32,
    out_offset: u32,
    patch_in_mon: u32,
    n_patches_in_mon: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

struct TimeMonitor {
    st: Stencil8,
    is_magnetic: u32,
    sample_stride: u32,
    write_offset: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    _pad3: u32,
    _pad4: u32,
}

struct NlDomain {
    d_glob: array<__R__, 18>,
    chi3_iso_xxxx: __R__,
    eps_min: __R__,
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
}

struct NlEdge {
    qpm_sign: __R__,
    cb_eps_inv_dt: __R__,
    cb_g: __R__,
    grid_idx: u32,
    ix: u32,
    iy: u32,
    iz: u32,
    comp: u32,
    dom_slot: u32,
    _pad0: u32,
    _pad1: u32,
}

struct NlRunParams {
    e_ref: __R__,
    inv_e_ref: __R__,
    _pad0: __R__,
    _pad1: __R__,
}

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<storage, read_write> fields_e: array<__R__>;
@group(0) @binding(2) var<storage, read_write> fields_h: array<__R__>;
@group(0) @binding(3) var<storage, read> coeffs_e: array<vec2<__R__>>;
@group(0) @binding(4) var<storage, read> coeffs_h: array<vec2<__R__>>;
@group(0) @binding(5) var<storage, read_write> psi_buf: array<__R__>;
@group(0) @binding(6) var<storage, read> tables: array<vec4<__R__>>;
@group(0) @binding(7) var<storage, read> pec_len: array<__R__>;

@group(1) @binding(0) var<uniform> step_u: StepUniform;
@group(1) @binding(1) var<storage, read> src_h_sites: array<SrcSite>;
@group(1) @binding(2) var<storage, read> src_h_terms: array<SrcTerm>;
@group(1) @binding(3) var<storage, read> src_h_waves: array<vec2<__R__>>;
@group(1) @binding(4) var<storage, read> src_e_sites: array<SrcSite>;
@group(1) @binding(5) var<storage, read> src_e_terms: array<SrcTerm>;
@group(1) @binding(6) var<storage, read> src_e_waves: array<vec2<__R__>>;
@group(1) @binding(7) var<storage, read> offdiag_entries: array<OffDiagEntry>;
@group(1) @binding(8) var<storage, read_write> offdiag_delta: array<__R__>;

@group(2) @binding(0) var<storage, read> dft_phasors: array<vec4<__R__>>;
@group(2) @binding(1) var<storage, read> planar_patches: array<PlanarPatch>;
@group(2) @binding(2) var<storage, read_write> planar_dft_out: array<vec4<__R__>>;
@group(2) @binding(3) var<storage, read_write> vol_dft_out: array<vec4<__R__>>;
@group(2) @binding(4) var<storage, read> time_monitors: array<TimeMonitor>;
@group(2) @binding(5) var<storage, read_write> time_out: array<__R__>;
@group(2) @binding(6) var<storage, read_write> energy_partials: array<__R__>;
@group(2) @binding(7) var<storage, read_write> energy_out: array<__R__>;

@group(3) @binding(0) var<storage, read> nl_domains: array<NlDomain>;
@group(3) @binding(1) var<storage, read> nl_chi3_rank4: array<__R__>;
@group(3) @binding(2) var<storage, read> nl_edges: array<NlEdge>;
// Flat state buffer: 0=e_old, 1=e_lin, 2=e_next, 3=p_prev, 4=n_carrier (stride = params.n_nl_edges)
@group(3) @binding(3) var<storage, read_write> nl_state: array<__R__>;
@group(3) @binding(4) var<uniform> nl_run: NlRunParams;

fn inv_dp_x(i: u32) -> __R__ { return tables[params.off_inv_dp_x + i].x; }
fn inv_dp_y(i: u32) -> __R__ { return tables[params.off_inv_dp_y + i].x; }
fn inv_dp_z(i: u32) -> __R__ { return tables[params.off_inv_dp_z + i].x; }
fn inv_dd_x(i: u32) -> __R__ { return tables[params.off_inv_dd_x + i].x; }
fn inv_dd_y(i: u32) -> __R__ { return tables[params.off_inv_dd_y + i].x; }
fn inv_dd_z(i: u32) -> __R__ { return tables[params.off_inv_dd_z + i].x; }

fn eval_stencil_e(st: Stencil8) -> __R__ {
    return st.w0.x * fields_e[st.idx0.x]
         + st.w0.y * fields_e[st.idx0.y]
         + st.w0.z * fields_e[st.idx0.z]
         + st.w0.w * fields_e[st.idx0.w]
         + st.w1.x * fields_e[st.idx1.x]
         + st.w1.y * fields_e[st.idx1.y]
         + st.w1.z * fields_e[st.idx1.z]
         + st.w1.w * fields_e[st.idx1.w];
}

fn eval_stencil_h(st: Stencil8) -> __R__ {
    return st.w0.x * fields_h[st.idx0.x]
         + st.w0.y * fields_h[st.idx0.y]
         + st.w0.z * fields_h[st.idx0.z]
         + st.w0.w * fields_h[st.idx0.w]
         + st.w1.x * fields_h[st.idx1.x]
         + st.w1.y * fields_h[st.idx1.y]
         + st.w1.z * fields_h[st.idx1.z]
         + st.w1.w * fields_h[st.idx1.w];
}

// psi_h layout (6 slabs):
// 0: psi_hy_x, 1: psi_hz_x, 2: psi_hx_y, 3: psi_hz_y, 4: psi_hx_z, 5: psi_hy_z
@compute @workgroup_size(256)
fn update_h(@builtin(global_invocation_id) gid: vec3<u32>) {
    let g = gid.x;
    let slab = params.slab_len;
    if (g >= slab) { return; }

    let sx = params.stride_x;
    let sy = params.stride_y;
    let ix = g / sx;
    let rem = g - ix * sx;
    let iy = rem / sy;
    let iz = rem - iy * sy;
    let nx = params.nx;
    let ny = params.ny;
    let nz = params.nz;
    if (ix > nx || iy > ny || iz > nz) { return; }

    let ey_off = slab;
    let ez_off = 2u * slab;
    let hy_off = slab;
    let hz_off = 2u * slab;
    let use_c = params.use_conf != 0u;

    if (params.dim_mode == 1u) {
        if (iy == 0u && ix <= nx && iz < nz) {
            var dey_dz: __R__ = (fields_e[ey_off + g + 1u] - fields_e[ey_off + g]) * inv_dp_z(iz);
            if (params.has_cpml_z != 0u) {
                let cpz = tables[params.off_cpml_h_z + iz];
                if (cpz.w > __R__(0.5)) {
                    let p_idx = 4u * slab + g;
                    let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dey_dz;
                    psi_buf[p_idx] = p_new;
                    dey_dz = cpz.x * dey_dz + p_new;
                }
            }
            let d = coeffs_h[g];
            fields_h[g] = d.x * fields_h[g] + d.y * dey_dz;
        }
        if (iy == 0u && ix < nx && iz <= nz) {
            var dey_dx: __R__ = (fields_e[ey_off + g + sx] - fields_e[ey_off + g]) * inv_dp_x(ix);
            if (params.has_cpml_x != 0u) {
                let cpx = tables[params.off_cpml_h_x + ix];
                if (cpx.w > __R__(0.5)) {
                    let p_idx = 1u * slab + g;
                    let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dey_dx;
                    psi_buf[p_idx] = p_new;
                    dey_dx = cpx.x * dey_dx + p_new;
                }
            }
            let d = coeffs_h[hz_off + g];
            fields_h[hz_off + g] = d.x * fields_h[hz_off + g] - d.y * dey_dx;
        }
        return;
    }

    if (params.dim_mode == 2u) {
        if (iy == 0u && ix < nx && iz < nz) {
            var ex0: __R__ = fields_e[g];
            var ex1: __R__ = fields_e[g + 1u];
            var ez0: __R__ = fields_e[ez_off + g];
            var ez1: __R__ = fields_e[ez_off + g + sx];
            if (use_c) {
                ex0 *= pec_len[g];
                ex1 *= pec_len[g + 1u];
                ez0 *= pec_len[ez_off + g];
                ez1 *= pec_len[ez_off + g + sx];
            }
            var dex_dz: __R__ = (ex1 - ex0) * inv_dp_z(iz);
            var dez_dx: __R__ = (ez1 - ez0) * inv_dp_x(ix);
            if (params.has_cpml_z != 0u) {
                let cpz = tables[params.off_cpml_h_z + iz];
                if (cpz.w > __R__(0.5)) {
                    let p_idx = 5u * slab + g;
                    let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dex_dz;
                    psi_buf[p_idx] = p_new;
                    dex_dz = cpz.x * dex_dz + p_new;
                }
            }
            if (params.has_cpml_x != 0u) {
                let cpx = tables[params.off_cpml_h_x + ix];
                if (cpx.w > __R__(0.5)) {
                    let p_idx = 0u * slab + g;
                    let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dez_dx;
                    psi_buf[p_idx] = p_new;
                    dez_dx = cpx.x * dez_dx + p_new;
                }
            }
            let d = coeffs_h[hy_off + g];
            fields_h[hy_off + g] = d.x * fields_h[hy_off + g] - d.y * (dex_dz - dez_dx);
        }
        return;
    }

    let is_2d = params.dim_mode == 3u;
    let iy_hi = select(ny, 0u, is_2d);

    if (ix <= nx && iy < ny && iz < nz) {
        var dez_dy: __R__ = __R__(0.0);
        if (!is_2d) {
            var ez0: __R__ = fields_e[ez_off + g];
            var ez1: __R__ = fields_e[ez_off + g + sy];
            if (use_c) {
                ez0 *= pec_len[ez_off + g];
                ez1 *= pec_len[ez_off + g + sy];
            }
            dez_dy = (ez1 - ez0) * inv_dp_y(iy);
            if (params.has_cpml_y != 0u) {
                let cpy = tables[params.off_cpml_h_y + iy];
                if (cpy.w > __R__(0.5)) {
                    let p_idx = 2u * slab + g;
                    let p_new = cpy.y * psi_buf[p_idx] + cpy.z * dez_dy;
                    psi_buf[p_idx] = p_new;
                    dez_dy = cpy.x * dez_dy + p_new;
                }
            }
        }
        var ey0: __R__ = fields_e[ey_off + g];
        var ey1: __R__ = fields_e[ey_off + g + 1u];
        if (use_c) {
            ey0 *= pec_len[ey_off + g];
            ey1 *= pec_len[ey_off + g + 1u];
        }
        var dey_dz: __R__ = (ey1 - ey0) * inv_dp_z(iz);
        if (params.has_cpml_z != 0u) {
            let cpz = tables[params.off_cpml_h_z + iz];
            if (cpz.w > __R__(0.5)) {
                let p_idx = 4u * slab + g;
                let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dey_dz;
                psi_buf[p_idx] = p_new;
                dey_dz = cpz.x * dey_dz + p_new;
            }
        }
        let d = coeffs_h[g];
        fields_h[g] = d.x * fields_h[g] - d.y * (dez_dy - dey_dz);
    }

    if (ix < nx && iy <= iy_hi && iz < nz) {
        var ex0: __R__ = fields_e[g];
        var ex1: __R__ = fields_e[g + 1u];
        var ez0: __R__ = fields_e[ez_off + g];
        var ez1: __R__ = fields_e[ez_off + g + sx];
        if (use_c) {
            ex0 *= pec_len[g];
            ex1 *= pec_len[g + 1u];
            ez0 *= pec_len[ez_off + g];
            ez1 *= pec_len[ez_off + g + sx];
        }
        var dex_dz: __R__ = (ex1 - ex0) * inv_dp_z(iz);
        var dez_dx: __R__ = (ez1 - ez0) * inv_dp_x(ix);
        if (params.has_cpml_z != 0u) {
            let cpz = tables[params.off_cpml_h_z + iz];
            if (cpz.w > __R__(0.5)) {
                let p_idx = 5u * slab + g;
                let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dex_dz;
                psi_buf[p_idx] = p_new;
                dex_dz = cpz.x * dex_dz + p_new;
            }
        }
        if (params.has_cpml_x != 0u) {
            let cpx = tables[params.off_cpml_h_x + ix];
            if (cpx.w > __R__(0.5)) {
                let p_idx = 0u * slab + g;
                let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dez_dx;
                psi_buf[p_idx] = p_new;
                dez_dx = cpx.x * dez_dx + p_new;
            }
        }
        let d = coeffs_h[hy_off + g];
        fields_h[hy_off + g] = d.x * fields_h[hy_off + g] - d.y * (dex_dz - dez_dx);
    }

    if (ix < nx && iy < ny && iz <= nz) {
        var ey0: __R__ = fields_e[ey_off + g];
        var ey1: __R__ = fields_e[ey_off + g + sx];
        if (use_c) {
            ey0 *= pec_len[ey_off + g];
            ey1 *= pec_len[ey_off + g + sx];
        }
        var dey_dx: __R__ = (ey1 - ey0) * inv_dp_x(ix);
        if (params.has_cpml_x != 0u) {
            let cpx = tables[params.off_cpml_h_x + ix];
            if (cpx.w > __R__(0.5)) {
                let p_idx = 1u * slab + g;
                let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dey_dx;
                psi_buf[p_idx] = p_new;
                dey_dx = cpx.x * dey_dx + p_new;
            }
        }
        var dex_dy: __R__ = __R__(0.0);
        if (!is_2d) {
            var ex0: __R__ = fields_e[g];
            var ex1: __R__ = fields_e[g + sy];
            if (use_c) {
                ex0 *= pec_len[g];
                ex1 *= pec_len[g + sy];
            }
            dex_dy = (ex1 - ex0) * inv_dp_y(iy);
            if (params.has_cpml_y != 0u) {
                let cpy = tables[params.off_cpml_h_y + iy];
                if (cpy.w > __R__(0.5)) {
                    let p_idx = 3u * slab + g;
                    let p_new = cpy.y * psi_buf[p_idx] + cpy.z * dex_dy;
                    psi_buf[p_idx] = p_new;
                    dex_dy = cpy.x * dex_dy + p_new;
                }
            }
        }
        let d = coeffs_h[hz_off + g];
        fields_h[hz_off + g] = d.x * fields_h[hz_off + g] - d.y * (dey_dx - dex_dy);
    }
}

@compute @workgroup_size(64)
fn inject_h(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    if (tid >= params.n_src_h_sites) { return; }
    let site = src_h_sites[tid];
    let step_base = step_u.step * params.n_src_h_chans;
    var acc: __R__ = __R__(0.0);
    for (var i = site.item_start; i < site.item_end; i = i + 1u) {
        let term = src_h_terms[i];
        let s = src_h_waves[step_base + term.chan];
        acc += term.w_re * s.x + term.w_im * s.y;
    }
    fields_h[site.target_offset] += acc;
}

fn eval_offdiag_rhs(e: OffDiagEntry) -> __R__ {
    return e.b_weights.x * fields_e[e.b_offsets.x]
         + e.b_weights.y * fields_e[e.b_offsets.y]
         + e.b_weights.z * fields_e[e.b_offsets.z]
         + e.b_weights.w * fields_e[e.b_offsets.w]
         + e.c_weights.x * fields_e[e.c_offsets.x]
         + e.c_weights.y * fields_e[e.c_offsets.y]
         + e.c_weights.z * fields_e[e.c_offsets.z]
         + e.c_weights.w * fields_e[e.c_offsets.w];
}

@compute @workgroup_size(64)
fn offdiag_save(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    if (tid >= params.n_offdiag) { return; }
    offdiag_delta[tid] = -eval_offdiag_rhs(offdiag_entries[tid]);
}

@compute @workgroup_size(64)
fn offdiag_apply(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    if (tid >= params.n_offdiag) { return; }
    offdiag_delta[tid] += eval_offdiag_rhs(offdiag_entries[tid]);
}

@compute @workgroup_size(64)
fn offdiag_commit(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    if (tid >= params.n_offdiag) { return; }
    let e = offdiag_entries[tid];
    fields_e[e.target_offset] += offdiag_delta[tid];
}

@compute @workgroup_size(256)
fn update_e(@builtin(global_invocation_id) gid: vec3<u32>) {
    let g = gid.x;
    let slab = params.slab_len;
    if (g >= slab) { return; }

    let sx = params.stride_x;
    let sy = params.stride_y;
    let ix = g / sx;
    let rem = g - ix * sx;
    let iy = rem / sy;
    let iz = rem - iy * sy;
    let nx = params.nx;
    let ny = params.ny;
    let nz = params.nz;
    if (ix > nx || iy > ny || iz > nz) { return; }

    let ey_off = slab;
    let ez_off = 2u * slab;
    let hy_off = slab;
    let hz_off = 2u * slab;

    let pmc_xm = (params.pmc_mask & 1u) != 0u;
    let pmc_xp = (params.pmc_mask & 2u) != 0u;
    let pmc_ym = (params.pmc_mask & 4u) != 0u;
    let pmc_yp = (params.pmc_mask & 8u) != 0u;
    let pmc_zm = (params.pmc_mask & 16u) != 0u;
    let pmc_zp = (params.pmc_mask & 32u) != 0u;

    let per_x = (params.per_mask & 1u) != 0u;
    let per_y = (params.per_mask & 2u) != 0u;
    let per_z = (params.per_mask & 4u) != 0u;

    let ix_lo = select(1u, 0u, pmc_xm || per_x);
    let ix_hi = select(nx - 1u, nx, pmc_xp);
    let iz_lo = select(1u, 0u, pmc_zm || per_z);
    let iz_hi = select(nz - 1u, nz, pmc_zp);

    if (params.dim_mode == 1u) {
        if (iy == 0u && ix >= ix_lo && ix <= ix_hi && iz >= iz_lo && iz <= iz_hi) {
            var dhx_dz: __R__ = __R__(0.0);
            if (iz == 0u) {
                if (per_z) {
                    dhx_dz = (fields_h[g] - fields_h[g + (nz - 1u)]) * inv_dd_z(0u);
                } else {
                    dhx_dz = __R__(2.0) * fields_h[g] * inv_dd_z(0u);
                }
            } else if (iz == nz) {
                dhx_dz = __R__(-2.0) * fields_h[g - 1u] * inv_dd_z(nz);
            } else {
                dhx_dz = (fields_h[g] - fields_h[g - 1u]) * inv_dd_z(iz);
            }
            if (params.has_cpml_z != 0u) {
                let cpz = tables[params.off_cpml_e_z + iz];
                if (cpz.w > __R__(0.5)) {
                    let p_idx = 5u * slab + g;
                    let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dhx_dz;
                    psi_buf[p_idx] = p_new;
                    dhx_dz = cpz.x * dhx_dz + p_new;
                }
            }

            var dhz_dx: __R__ = __R__(0.0);
            if (ix == 0u) {
                if (per_x) {
                    dhz_dx = (fields_h[hz_off + g] - fields_h[hz_off + g + (nx - 1u) * sx]) * inv_dd_x(0u);
                } else {
                    dhz_dx = __R__(2.0) * fields_h[hz_off + g] * inv_dd_x(0u);
                }
            } else if (ix == nx) {
                dhz_dx = __R__(-2.0) * fields_h[hz_off + g - sx] * inv_dd_x(nx);
            } else {
                dhz_dx = (fields_h[hz_off + g] - fields_h[hz_off + g - sx]) * inv_dd_x(ix);
            }
            if (params.has_cpml_x != 0u) {
                let cpx = tables[params.off_cpml_e_x + ix];
                if (cpx.w > __R__(0.5)) {
                    let p_idx = 0u * slab + g;
                    let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dhz_dx;
                    psi_buf[p_idx] = p_new;
                    dhz_dx = cpx.x * dhz_dx + p_new;
                }
            }
            let c = coeffs_e[ey_off + g];
            fields_e[ey_off + g] = c.x * fields_e[ey_off + g] + c.y * (dhx_dz - dhz_dx);
        }
        return;
    }

    if (params.dim_mode == 2u) {
        if (iy == 0u && ix < nx && iz >= iz_lo && iz <= iz_hi) {
            var dhy_dz: __R__ = __R__(0.0);
            if (iz == 0u) {
                if (per_z) {
                    dhy_dz = (fields_h[hy_off + g] - fields_h[hy_off + g + (nz - 1u)]) * inv_dd_z(0u);
                } else {
                    dhy_dz = __R__(2.0) * fields_h[hy_off + g] * inv_dd_z(0u);
                }
            } else if (iz == nz) {
                dhy_dz = __R__(-2.0) * fields_h[hy_off + g - 1u] * inv_dd_z(nz);
            } else {
                dhy_dz = (fields_h[hy_off + g] - fields_h[hy_off + g - 1u]) * inv_dd_z(iz);
            }
            if (params.has_cpml_z != 0u) {
                let cpz = tables[params.off_cpml_e_z + iz];
                if (cpz.w > __R__(0.5)) {
                    let p_idx = 4u * slab + g;
                    let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dhy_dz;
                    psi_buf[p_idx] = p_new;
                    dhy_dz = cpz.x * dhy_dz + p_new;
                }
            }
            let c = coeffs_e[g];
            fields_e[g] = c.x * fields_e[g] - c.y * dhy_dz;
        }
        if (iy == 0u && ix >= ix_lo && ix <= ix_hi && iz < nz) {
            var dhy_dx: __R__ = __R__(0.0);
            if (ix == 0u) {
                if (per_x) {
                    dhy_dx = (fields_h[hy_off + g] - fields_h[hy_off + g + (nx - 1u) * sx]) * inv_dd_x(0u);
                } else {
                    dhy_dx = __R__(2.0) * fields_h[hy_off + g] * inv_dd_x(0u);
                }
            } else if (ix == nx) {
                dhy_dx = __R__(-2.0) * fields_h[hy_off + g - sx] * inv_dd_x(nx);
            } else {
                dhy_dx = (fields_h[hy_off + g] - fields_h[hy_off + g - sx]) * inv_dd_x(ix);
            }
            if (params.has_cpml_x != 0u) {
                let cpx = tables[params.off_cpml_e_x + ix];
                if (cpx.w > __R__(0.5)) {
                    let p_idx = 1u * slab + g;
                    let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dhy_dx;
                    psi_buf[p_idx] = p_new;
                    dhy_dx = cpx.x * dhy_dx + p_new;
                }
            }
            let c = coeffs_e[ez_off + g];
            fields_e[ez_off + g] = c.x * fields_e[ez_off + g] + c.y * dhy_dx;
        }
        return;
    }

    let is_2d = params.dim_mode == 3u;
    let iy_lo = select(select(1u, 0u, pmc_ym || per_y), 0u, is_2d);
    let iy_hi = select(select(ny - 1u, ny, pmc_yp), 0u, is_2d);

    if (ix < nx && iy >= iy_lo && iy <= iy_hi && iz >= iz_lo && iz <= iz_hi) {
        var dhz_dy: __R__ = __R__(0.0);
        if (!is_2d) {
            if (iy == 0u) {
                if (per_y) {
                    dhz_dy = (fields_h[hz_off + g] - fields_h[hz_off + g + (ny - 1u) * sy]) * inv_dd_y(0u);
                } else {
                    dhz_dy = __R__(2.0) * fields_h[hz_off + g] * inv_dd_y(0u);
                }
            } else if (iy == ny) {
                dhz_dy = __R__(-2.0) * fields_h[hz_off + g - sy] * inv_dd_y(ny);
            } else {
                dhz_dy = (fields_h[hz_off + g] - fields_h[hz_off + g - sy]) * inv_dd_y(iy);
            }
            if (params.has_cpml_y != 0u) {
                let cpy = tables[params.off_cpml_e_y + iy];
                if (cpy.w > __R__(0.5)) {
                    let p_idx = 2u * slab + g;
                    let p_new = cpy.y * psi_buf[p_idx] + cpy.z * dhz_dy;
                    psi_buf[p_idx] = p_new;
                    dhz_dy = cpy.x * dhz_dy + p_new;
                }
            }
        }
        var dhy_dz: __R__ = __R__(0.0);
        if (iz == 0u) {
            if (per_z) {
                dhy_dz = (fields_h[hy_off + g] - fields_h[hy_off + g + (nz - 1u)]) * inv_dd_z(0u);
            } else {
                dhy_dz = __R__(2.0) * fields_h[hy_off + g] * inv_dd_z(0u);
            }
        } else if (iz == nz) {
            dhy_dz = __R__(-2.0) * fields_h[hy_off + g - 1u] * inv_dd_z(nz);
        } else {
            dhy_dz = (fields_h[hy_off + g] - fields_h[hy_off + g - 1u]) * inv_dd_z(iz);
        }
        if (params.has_cpml_z != 0u) {
            let cpz = tables[params.off_cpml_e_z + iz];
            if (cpz.w > __R__(0.5)) {
                let p_idx = 4u * slab + g;
                let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dhy_dz;
                psi_buf[p_idx] = p_new;
                dhy_dz = cpz.x * dhy_dz + p_new;
            }
        }
        let c = coeffs_e[g];
        fields_e[g] = c.x * fields_e[g] + c.y * (dhz_dy - dhy_dz);
    }

    if (ix >= ix_lo && ix <= ix_hi && iy < ny && iz >= iz_lo && iz <= iz_hi) {
        var dhx_dz: __R__ = __R__(0.0);
        if (iz == 0u) {
            if (per_z) {
                dhx_dz = (fields_h[g] - fields_h[g + (nz - 1u)]) * inv_dd_z(0u);
            } else {
                dhx_dz = __R__(2.0) * fields_h[g] * inv_dd_z(0u);
            }
        } else if (iz == nz) {
            dhx_dz = __R__(-2.0) * fields_h[g - 1u] * inv_dd_z(nz);
        } else {
            dhx_dz = (fields_h[g] - fields_h[g - 1u]) * inv_dd_z(iz);
        }
        if (params.has_cpml_z != 0u) {
            let cpz = tables[params.off_cpml_e_z + iz];
            if (cpz.w > __R__(0.5)) {
                let p_idx = 5u * slab + g;
                let p_new = cpz.y * psi_buf[p_idx] + cpz.z * dhx_dz;
                psi_buf[p_idx] = p_new;
                dhx_dz = cpz.x * dhx_dz + p_new;
            }
        }
        var dhz_dx: __R__ = __R__(0.0);
        if (ix == 0u) {
            if (per_x) {
                dhz_dx = (fields_h[hz_off + g] - fields_h[hz_off + g + (nx - 1u) * sx]) * inv_dd_x(0u);
            } else {
                dhz_dx = __R__(2.0) * fields_h[hz_off + g] * inv_dd_x(0u);
            }
        } else if (ix == nx) {
            dhz_dx = __R__(-2.0) * fields_h[hz_off + g - sx] * inv_dd_x(nx);
        } else {
            dhz_dx = (fields_h[hz_off + g] - fields_h[hz_off + g - sx]) * inv_dd_x(ix);
        }
        if (params.has_cpml_x != 0u) {
            let cpx = tables[params.off_cpml_e_x + ix];
            if (cpx.w > __R__(0.5)) {
                let p_idx = 0u * slab + g;
                let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dhz_dx;
                psi_buf[p_idx] = p_new;
                dhz_dx = cpx.x * dhz_dx + p_new;
            }
        }
        let c = coeffs_e[ey_off + g];
        fields_e[ey_off + g] = c.x * fields_e[ey_off + g] + c.y * (dhx_dz - dhz_dx);
    }

    if (ix >= ix_lo && ix <= ix_hi && iy >= iy_lo && iy <= iy_hi && iz < nz) {
        var dhy_dx: __R__ = __R__(0.0);
        if (ix == 0u) {
            if (per_x) {
                dhy_dx = (fields_h[hy_off + g] - fields_h[hy_off + g + (nx - 1u) * sx]) * inv_dd_x(0u);
            } else {
                dhy_dx = __R__(2.0) * fields_h[hy_off + g] * inv_dd_x(0u);
            }
        } else if (ix == nx) {
            dhy_dx = __R__(-2.0) * fields_h[hy_off + g - sx] * inv_dd_x(nx);
        } else {
            dhy_dx = (fields_h[hy_off + g] - fields_h[hy_off + g - sx]) * inv_dd_x(ix);
        }
        if (params.has_cpml_x != 0u) {
            let cpx = tables[params.off_cpml_e_x + ix];
            if (cpx.w > __R__(0.5)) {
                let p_idx = 1u * slab + g;
                let p_new = cpx.y * psi_buf[p_idx] + cpx.z * dhy_dx;
                psi_buf[p_idx] = p_new;
                dhy_dx = cpx.x * dhy_dx + p_new;
            }
        }
        var dhx_dy: __R__ = __R__(0.0);
        if (!is_2d) {
            if (iy == 0u) {
                if (per_y) {
                    dhx_dy = (fields_h[g] - fields_h[g + (ny - 1u) * sy]) * inv_dd_y(0u);
                } else {
                    dhx_dy = __R__(2.0) * fields_h[g] * inv_dd_y(0u);
                }
            } else if (iy == ny) {
                dhx_dy = __R__(-2.0) * fields_h[g - sy] * inv_dd_y(ny);
            } else {
                dhx_dy = (fields_h[g] - fields_h[g - sy]) * inv_dd_y(iy);
            }
            if (params.has_cpml_y != 0u) {
                let cpy = tables[params.off_cpml_e_y + iy];
                if (cpy.w > __R__(0.5)) {
                    let p_idx = 3u * slab + g;
                    let p_new = cpy.y * psi_buf[p_idx] + cpy.z * dhx_dy;
                    psi_buf[p_idx] = p_new;
                    dhx_dy = cpy.x * dhx_dy + p_new;
                }
            }
        }
        let c = coeffs_e[ez_off + g];
        fields_e[ez_off + g] = c.x * fields_e[ez_off + g] + c.y * (dhy_dx - dhx_dy);
    }
}

@compute @workgroup_size(64)
fn inject_e(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    if (tid >= params.n_src_e_sites) { return; }
    let site = src_e_sites[tid];
    let step_base = step_u.step * params.n_src_e_chans;
    var acc: __R__ = __R__(0.0);
    for (var i = site.item_start; i < site.item_end; i = i + 1u) {
        let term = src_e_terms[i];
        let s = src_e_waves[step_base + term.chan];
        acc += term.w_re * s.x + term.w_im * s.y;
    }
    fields_e[site.target_offset] += acc;
}

@compute @workgroup_size(64)
fn sync_periodic_x(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    let is_2d = params.dim_mode != 0u;
    let iy_hi = select(params.ny, 0u, is_2d);
    let nz1 = params.nz + 1u;
    let total = (iy_hi + 1u) * nz1;
    if (tid >= total) { return; }
    let iy = tid / nz1;
    let iz = tid - iy * nz1;
    let g0 = iy * params.stride_y + iz;
    let gn = params.nx * params.stride_x + g0;
    let slab = params.slab_len;
    fields_e[slab + gn] = fields_e[slab + g0];
    fields_e[2u * slab + gn] = fields_e[2u * slab + g0];
}

@compute @workgroup_size(64)
fn sync_periodic_y(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    let nz1 = params.nz + 1u;
    let total = (params.nx + 1u) * nz1;
    if (tid >= total) { return; }
    let ix = tid / nz1;
    let iz = tid - ix * nz1;
    let g0 = ix * params.stride_x + iz;
    let gn = g0 + params.ny * params.stride_y;
    let slab = params.slab_len;
    fields_e[gn] = fields_e[g0];
    fields_e[2u * slab + gn] = fields_e[2u * slab + g0];
}

@compute @workgroup_size(64)
fn sync_periodic_z(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tid = gid.x;
    let is_2d = params.dim_mode != 0u;
    let iy_hi = select(params.ny, 0u, is_2d);
    let ny1 = iy_hi + 1u;
    let total = (params.nx + 1u) * ny1;
    if (tid >= total) { return; }
    let ix = tid / ny1;
    let iy = tid - ix * ny1;
    let g0 = ix * params.stride_x + iy * params.stride_y;
    let gn = g0 + params.nz;
    let slab = params.slab_len;
    fields_e[gn] = fields_e[g0];
    fields_e[slab + gn] = fields_e[slab + g0];
}

@compute @workgroup_size(256)
fn clamp_pec_e(@builtin(global_invocation_id) gid: vec3<u32>) {
    let g = gid.x;
    let slab = params.slab_len;
    if (g >= slab) { return; }
    let sx = params.stride_x;
    let sy = params.stride_y;
    let ix = g / sx;
    let rem = g - ix * sx;
    let iy = rem / sy;
    let iz = rem - iy * sy;
    if (ix > params.nx || iy > params.ny || iz > params.nz) { return; }

    let ey_off = slab;
    let ez_off = 2u * slab;
    let is_2d = params.dim_mode != 0u;

    if (((params.pec_mask & 1u) != 0u && ix == 0u) || ((params.pec_mask & 2u) != 0u && ix == params.nx)) {
        fields_e[ey_off + g] = __R__(0.0);
        fields_e[ez_off + g] = __R__(0.0);
    }
    if (!is_2d && (((params.pec_mask & 4u) != 0u && iy == 0u) || ((params.pec_mask & 8u) != 0u && iy == params.ny))) {
        fields_e[g] = __R__(0.0);
        fields_e[ez_off + g] = __R__(0.0);
    }
    if (((params.pec_mask & 16u) != 0u && iz == 0u) || ((params.pec_mask & 32u) != 0u && iz == params.nz)) {
        fields_e[g] = __R__(0.0);
        fields_e[ey_off + g] = __R__(0.0);
    }
    if (params.use_conf != 0u) {
        if (pec_len[g] <= __R__(0.0)) { fields_e[g] = __R__(0.0); }
        if (pec_len[ey_off + g] <= __R__(0.0)) { fields_e[ey_off + g] = __R__(0.0); }
        if (pec_len[ez_off + g] <= __R__(0.0)) { fields_e[ez_off + g] = __R__(0.0); }
    }
}

// ---------------------------------------------------------------------------
// Nonlinear chi^(2)/QPM/chi^(3)/TPA/FCA/FCD Kernels (`Group 3`)
// ---------------------------------------------------------------------------

fn sample_edge_e_vec(ed: NlEdge) -> vec3<__R__> {
    let sx = params.stride_x;
    let sy = params.stride_y;
    let slab = params.slab_len;
    let ey_off = slab;
    let ez_off = 2u * slab;
    let is_2d = params.dim_mode != 0u;
    let ix = ed.ix;
    let iy = ed.iy;
    let iz = ed.iz;

    if (ed.comp == 0u) {
        let ex_v = fields_e[ed.grid_idx];
        var ey_v: __R__ = __R__(0.0);
        if (is_2d) {
            let g0 = ix * sx + iz;
            let g1 = (ix + 1u) * sx + iz;
            ey_v = __R__(0.5) * (fields_e[ey_off + g0] + fields_e[ey_off + g1]);
        } else {
            let ym = select(0u, iy - 1u, iy > 0u);
            let yp = min(iy, params.ny - 1u);
            let g0 = ix * sx + ym * sy + iz;
            let g1 = (ix + 1u) * sx + ym * sy + iz;
            let g2 = ix * sx + yp * sy + iz;
            let g3 = (ix + 1u) * sx + yp * sy + iz;
            ey_v = __R__(0.25) * (fields_e[ey_off + g0] + fields_e[ey_off + g1] + fields_e[ey_off + g2] + fields_e[ey_off + g3]);
        }
        let zm = select(0u, iz - 1u, iz > 0u);
        let zp = min(iz, params.nz - 1u);
        let c0 = ix * sx + iy * sy + zm;
        let c1 = (ix + 1u) * sx + iy * sy + zm;
        let c2 = ix * sx + iy * sy + zp;
        let c3 = (ix + 1u) * sx + iy * sy + zp;
        let ez_v = __R__(0.25) * (fields_e[ez_off + c0] + fields_e[ez_off + c1] + fields_e[ez_off + c2] + fields_e[ez_off + c3]);
        return vec3<__R__>(ex_v, ey_v, ez_v);
    } else if (ed.comp == 1u) {
        let ey_v = fields_e[ey_off + ed.grid_idx];
        let xm = select(0u, ix - 1u, ix > 0u);
        let xp = min(ix, params.nx - 1u);
        var ex_v: __R__ = __R__(0.0);
        if (is_2d) {
            ex_v = __R__(0.5) * (fields_e[xm * sx + iz] + fields_e[xp * sx + iz]);
        } else {
            let a0 = xm * sx + iy * sy + iz;
            let a1 = xp * sx + iy * sy + iz;
            let a2 = xm * sx + (iy + 1u) * sy + iz;
            let a3 = xp * sx + (iy + 1u) * sy + iz;
            ex_v = __R__(0.25) * (fields_e[a0] + fields_e[a1] + fields_e[a2] + fields_e[a3]);
        }
        let zm = select(0u, iz - 1u, iz > 0u);
        let zp = min(iz, params.nz - 1u);
        var ez_v: __R__ = __R__(0.0);
        if (is_2d) {
            ez_v = __R__(0.5) * (fields_e[ez_off + ix * sx + zm] + fields_e[ez_off + ix * sx + zp]);
        } else {
            let c0 = ix * sx + iy * sy + zm;
            let c1 = ix * sx + iy * sy + zp;
            let c2 = ix * sx + (iy + 1u) * sy + zm;
            let c3 = ix * sx + (iy + 1u) * sy + zp;
            ez_v = __R__(0.25) * (fields_e[ez_off + c0] + fields_e[ez_off + c1] + fields_e[ez_off + c2] + fields_e[ez_off + c3]);
        }
        return vec3<__R__>(ex_v, ey_v, ez_v);
    } else {
        let ez_v = fields_e[ez_off + ed.grid_idx];
        let xm = select(0u, ix - 1u, ix > 0u);
        let xp = min(ix, params.nx - 1u);
        let a0 = xm * sx + iy * sy + iz;
        let a1 = xp * sx + iy * sy + iz;
        let a2 = xm * sx + iy * sy + (iz + 1u);
        let a3 = xp * sx + iy * sy + (iz + 1u);
        let ex_v = __R__(0.25) * (fields_e[a0] + fields_e[a1] + fields_e[a2] + fields_e[a3]);
        var ey_v: __R__ = __R__(0.0);
        if (is_2d) {
            let b0 = ix * sx + iz;
            let b1 = ix * sx + (iz + 1u);
            ey_v = __R__(0.5) * (fields_e[ey_off + b0] + fields_e[ey_off + b1]);
        } else {
            let ym = select(0u, iy - 1u, iy > 0u);
            let yp = min(iy, params.ny - 1u);
            let b0 = ix * sx + ym * sy + iz;
            let b1 = ix * sx + yp * sy + iz;
            let b2 = ix * sx + ym * sy + (iz + 1u);
            let b3 = ix * sx + yp * sy + (iz + 1u);
            ey_v = __R__(0.25) * (fields_e[ey_off + b0] + fields_e[ey_off + b1] + fields_e[ey_off + b2] + fields_e[ey_off + b3]);
        }
        return vec3<__R__>(ex_v, ey_v, ez_v);
    }
}

fn eval_nl_scaled(
    dom: NlDomain,
    comp: u32,
    qpm_sign: __R__,
    en: vec3<__R__>,
    nc_tilde: __R__,
) -> vec3<__R__> {
    // Returns (p_tilde_comp, chi_aa, sig_nl_si) using E_ref-normalized `en` and `nc_tilde`.
    let en2 = dot(en, en);
    var p_tilde: __R__ = __R__(0.0);
    var chi_aa: __R__ = __R__(0.0);

    if (dom.has_chi2 != 0u) {
        let base = comp * 6u;
        let d0 = dom.d_glob[base + 0u];
        let d1 = dom.d_glob[base + 1u];
        let d2 = dom.d_glob[base + 2u];
        let d3 = dom.d_glob[base + 3u];
        let d4 = dom.d_glob[base + 4u];
        let d5 = dom.d_glob[base + 5u];
        let v0 = en.x * en.x;
        let v1 = en.y * en.y;
        let v2 = en.z * en.z;
        let v3 = __R__(2.0) * en.y * en.z;
        let v4 = __R__(2.0) * en.x * en.z;
        let v5 = __R__(2.0) * en.x * en.y;
        p_tilde += __R__(2.0) * qpm_sign * (d0 * v0 + d1 * v1 + d2 * v2 + d3 * v3 + d4 * v4 + d5 * v5);
        if (comp == 0u) {
            chi_aa += __R__(4.0) * qpm_sign * (d0 * en.x + d4 * en.z + d5 * en.y);
        } else if (comp == 1u) {
            chi_aa += __R__(4.0) * qpm_sign * (d1 * en.y + d3 * en.z + d5 * en.x);
        } else {
            chi_aa += __R__(4.0) * qpm_sign * (d2 * en.z + d3 * en.y + d4 * en.x);
        }
    }

    if (dom.rank4_offset != 0xffffffffu) {
        let r_base = dom.rank4_offset;
        var s: __R__ = __R__(0.0);
        var t: __R__ = __R__(0.0);
        for (var j = 0u; j < 3u; j = j + 1u) {
            let ej = select(select(en.x, en.y, j == 1u), en.z, j == 2u);
            for (var k = 0u; k < 3u; k = k + 1u) {
                let ek = select(select(en.x, en.y, k == 1u), en.z, k == 2u);
                let ejk = ej * ek;
                for (var l = 0u; l < 3u; l = l + 1u) {
                    let el = select(select(en.x, en.y, l == 1u), en.z, l == 2u);
                    let r_cjkl = nl_chi3_rank4[r_base + comp * 27u + j * 9u + k * 3u + l];
                    s += r_cjkl * ejk * el;
                }
                let r_cjkc = nl_chi3_rank4[r_base + comp * 27u + j * 9u + k * 3u + comp];
                let r_cjck = nl_chi3_rank4[r_base + comp * 27u + j * 9u + comp * 3u + k];
                let r_ccjk = nl_chi3_rank4[r_base + comp * 27u + comp * 9u + j * 3u + k];
                t += (r_cjkc + r_cjck + r_ccjk) * ejk;
            }
        }
        p_tilde += s;
        chi_aa += t;
    } else if (dom.has_chi3_iso != 0u) {
        let ec = select(select(en.x, en.y, comp == 1u), en.z, comp == 2u);
        p_tilde += dom.chi3_iso_xxxx * en2 * ec;
        chi_aa += dom.chi3_iso_xxxx * (en2 + __R__(2.0) * ec * ec);
    }

    if (dom.fcd_coeff != __R__(0.0) && nc_tilde > __R__(0.0)) {
        let dfcd = dom.fcd_coeff * nc_tilde;
        let ec = select(select(en.x, en.y, comp == 1u), en.z, comp == 2u);
        p_tilde += dfcd * ec;
        chi_aa += dfcd;
    }

    let p_cap = __R__(0.45) * dom.eps_min * sqrt(en2);
    if (abs(p_tilde) > p_cap && p_cap > __R__(0.0)) {
        p_tilde = select(-p_cap, p_cap, p_tilde >= __R__(0.0));
    }
    let tan_cap = __R__(0.45) * dom.eps_min;
    chi_aa = clamp(chi_aa, -tan_cap, tan_cap);

    var sig_nl: __R__ = dom.sig_tpa_coeff * en2 + dom.sig_fca_coeff * nc_tilde;
    if (sig_nl < __R__(0.0)) {
        sig_nl = __R__(0.0);
    }
    return vec3<__R__>(p_tilde, chi_aa, sig_nl);
}

@compute @workgroup_size(64)
fn nl_save_old(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let n = params.n_nl_edges;
    if (k >= n) { return; }
    let ed = nl_edges[k];
    nl_state[k] = fields_e[ed.comp * params.slab_len + ed.grid_idx];
}

@compute @workgroup_size(64)
fn nl_save_lin(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let n = params.n_nl_edges;
    if (k >= n) { return; }
    let ed = nl_edges[k];
    nl_state[n + k] = fields_e[ed.comp * params.slab_len + ed.grid_idx];
}

@compute @workgroup_size(64)
fn nl_newton_iter(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let n = params.n_nl_edges;
    if (k >= n) { return; }
    let ed = nl_edges[k];
    let dom = nl_domains[ed.dom_slot];
    let e_vec = sample_edge_e_vec(ed);
    let en = e_vec * nl_run.inv_e_ref;
    let en_m = select(select(en.x, en.y, ed.comp == 1u), en.z, ed.comp == 2u);
    let en_old = nl_state[k] * nl_run.inv_e_ref;
    let en_lin = nl_state[n + k] * nl_run.inv_e_ref;
    let p_prev = nl_state[3u * n + k];
    let nc_tilde = nl_state[4u * n + k];

    let eval_res = eval_nl_scaled(dom, ed.comp, ed.qpm_sign, en, nc_tilde);
    let p_m = eval_res.x;
    let chi_aa = eval_res.y;
    let sig_nl = eval_res.z;

    let res = en_m - en_lin
        + ed.cb_eps_inv_dt * (p_m - p_prev)
        + ed.cb_g * sig_nl * __R__(0.5) * (en_m + en_old);
    let jac = max(__R__(1.0) + ed.cb_eps_inv_dt * chi_aa + __R__(0.5) * ed.cb_g * sig_nl, __R__(0.25));
    let upd = en_m - res / jac;
    nl_state[2u * n + k] = upd * nl_run.e_ref;
}

@compute @workgroup_size(64)
fn nl_commit_e(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let n = params.n_nl_edges;
    if (k >= n) { return; }
    let ed = nl_edges[k];
    fields_e[ed.comp * params.slab_len + ed.grid_idx] = nl_state[2u * n + k];
}

@compute @workgroup_size(64)
fn nl_finish_step(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    let n = params.n_nl_edges;
    if (k >= n) { return; }
    let ed = nl_edges[k];
    let dom = nl_domains[ed.dom_slot];
    let en = sample_edge_e_vec(ed) * nl_run.inv_e_ref;
    let nc_tilde = nl_state[4u * n + k];
    let eval_res = eval_nl_scaled(dom, ed.comp, ed.qpm_sign, en, nc_tilde);
    nl_state[3u * n + k] = eval_res.x;

    if (params.nl_has_carriers != 0u && dom.has_carrier != 0u) {
        let en_old = nl_state[k] * nl_run.inv_e_ref;
        var en_mid = en;
        if (ed.comp == 0u) {
            en_mid.x = __R__(0.5) * (en.x + en_old);
        } else if (ed.comp == 1u) {
            en_mid.y = __R__(0.5) * (en.y + en_old);
        } else {
            en_mid.z = __R__(0.5) * (en.z + en_old);
        }
        let en2_mid = dot(en_mid, en_mid);
        let nc_next = max(nc_tilde * dom.decay + dom.tpa_gen_coeff * en2_mid * en2_mid * dom.one_minus_decay_tau, __R__(0.0));
        nl_state[4u * n + k] = nc_next;
    }
}

// ---------------------------------------------------------------------------
// DFT Monitors, Time Probes, and Energy Reduction
// ---------------------------------------------------------------------------

@compute @workgroup_size(64)
fn accum_planar_dft(@builtin(global_invocation_id) gid: vec3<u32>) {
    let pid = gid.x;
    if (pid >= params.n_planar_patches) { return; }
    let pt = planar_patches[pid];
    let eu = eval_stencil_e(pt.st_eu);
    let ev = eval_stencil_e(pt.st_ev);
    let hu = eval_stencil_h(pt.st_hu);
    let hv = eval_stencil_h(pt.st_hv);
    let step_ph_base = step_u.step * params.total_dft_freqs + pt.freq_offset;

    for (var fi = 0u; fi < pt.n_freqs; fi = fi + 1u) {
        let ph = dft_phasors[step_ph_base + fi];
        let slot = (pt.out_offset + fi * pt.n_patches_in_mon + pt.patch_in_mon) * 2u;
        planar_dft_out[slot] += vec4<__R__>(eu * ph.x, eu * ph.y, ev * ph.x, ev * ph.y);
        planar_dft_out[slot + 1u] += vec4<__R__>(hu * ph.z, hu * ph.w, hv * ph.z, hv * ph.w);
    }
}

@compute @workgroup_size(64)
fn accum_vol_dft(@builtin(global_invocation_id) gid: vec3<u32>) {
    let vid = gid.x;
    if (vid >= params.n_vol_pts) { return; }
    let ny_v = params.vol_j1 - params.vol_j0;
    let nz_v = params.vol_k1 - params.vol_k0;
    let yz_v = ny_v * nz_v;
    let loc_i = vid / yz_v;
    let rem = vid - loc_i * yz_v;
    let loc_j = rem / nz_v;
    let loc_k = rem - loc_j * nz_v;
    let ix = params.vol_i0 + loc_i;
    let iy = params.vol_j0 + loc_j;
    let iz = params.vol_k0 + loc_k;

    let sx = params.stride_x;
    let sy = params.stride_y;
    let slab = params.slab_len;
    let ey_off = slab;
    let ez_off = 2u * slab;
    let hy_off = slab;
    let hz_off = 2u * slab;

    let g000 = ix * sx + iy * sy + iz;
    let g100 = g000 + sx;
    let g001 = g000 + 1u;
    let g101 = g100 + 1u;

    let is_2d = params.dim_mode != 0u;
    var ex_c: __R__ = __R__(0.0);
    let ey_c = __R__(0.25) * (fields_e[ey_off + g000] + fields_e[ey_off + g100] + fields_e[ey_off + g001] + fields_e[ey_off + g101]);
    var ez_c: __R__ = __R__(0.0);
    let hx_c = __R__(0.5) * (fields_h[g000] + fields_h[g100]);
    var hy_c: __R__ = __R__(0.0);
    let hz_c = __R__(0.5) * (fields_h[hz_off + g000] + fields_h[hz_off + g001]);

    if (is_2d) {
        ex_c = __R__(0.5) * (fields_e[g000] + fields_e[g001]);
        ez_c = __R__(0.5) * (fields_e[ez_off + g000] + fields_e[ez_off + g100]);
        hy_c = fields_h[hy_off + g000];
    } else {
        let g010 = g000 + sy;
        let g110 = g100 + sy;
        let g011 = g010 + 1u;
        ex_c = __R__(0.25) * (fields_e[g000] + fields_e[g010] + fields_e[g001] + fields_e[g011]);
        ez_c = __R__(0.25) * (fields_e[ez_off + g000] + fields_e[ez_off + g100] + fields_e[ez_off + g010] + fields_e[ez_off + g110]);
        hy_c = __R__(0.5) * (fields_h[hy_off + g000] + fields_h[hy_off + g010]);
    }

    let step_ph_base = step_u.step * params.total_dft_freqs + params.vol_freq_offset;
    for (var fi = 0u; fi < params.n_vol_freqs; fi = fi + 1u) {
        let ph = dft_phasors[step_ph_base + fi];
        let base_slot = (fi * params.n_vol_pts + vid) * 3u;
        vol_dft_out[base_slot] += vec4<__R__>(ex_c * ph.x, ex_c * ph.y, ey_c * ph.x, ey_c * ph.y);
        vol_dft_out[base_slot + 1u] += vec4<__R__>(ez_c * ph.x, ez_c * ph.y, hx_c * ph.z, hx_c * ph.w);
        vol_dft_out[base_slot + 2u] += vec4<__R__>(hy_c * ph.z, hy_c * ph.w, hz_c * ph.z, hz_c * ph.w);
    }
}

@compute @workgroup_size(64)
fn sample_time_monitors(@builtin(global_invocation_id) gid: vec3<u32>) {
    let mid = gid.x;
    if (mid >= params.n_time_monitors) { return; }
    let tm = time_monitors[mid];
    if (step_u.step % tm.sample_stride != 0u) { return; }
    let val = select(eval_stencil_e(tm.st), eval_stencil_h(tm.st), tm.is_magnetic != 0u);
    let s_idx = step_u.step / tm.sample_stride;
    time_out[tm.write_offset + s_idx] = val;
}

var<workgroup> wg_sum: array<__R__, 256>;

@compute @workgroup_size(256)
fn energy_reduce_pass1(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
) {
    let g = gid.x;
    let l = lid.x;
    let slab = params.slab_len;
    var val: __R__ = __R__(0.0);
    if (g < slab) {
        let ex = fields_e[g];
        let ey = fields_e[slab + g];
        let ez = fields_e[2u * slab + g];
        let hx = fields_h[g];
        let hy = fields_h[slab + g];
        let hz = fields_h[2u * slab + g];
        let eta0_sq = __R__(376.730313668) * __R__(376.730313668);
        val = ex * ex + ey * ey + ez * ez + eta0_sq * (hx * hx + hy * hy + hz * hz);
    }
    wg_sum[l] = val;
    workgroupBarrier();

    for (var s = 128u; s > 0u; s = s >> 1u) {
        if (l < s) {
            wg_sum[l] += wg_sum[l + s];
        }
        workgroupBarrier();
    }
    if (l == 0u) {
        energy_partials[wid.x] = wg_sum[0];
    }
}

@compute @workgroup_size(256)
fn energy_reduce_pass2(@builtin(local_invocation_id) lid: vec3<u32>) {
    let l = lid.x;
    let n_wgs = params.n_energy_wgs;
    var acc: __R__ = __R__(0.0);
    var i = l;
    while (i < n_wgs) {
        acc += energy_partials[i];
        i += 256u;
    }
    wg_sum[l] = acc;
    workgroupBarrier();

    for (var s = 128u; s > 0u; s = s >> 1u) {
        if (l < s) {
            wg_sum[l] += wg_sum[l + s];
        }
        workgroupBarrier();
    }
    if (l == 0u) {
        energy_out[0] = wg_sum[0];
    }
}
"#;
