// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Hardware-accelerated GPU linear-algebra kernels via `wgpu` + WGSL (`FR-43a`, `ADR 146`).
//!
//! Metadata types ([`GpuDeviceInfo`], [`GpuError`], [`gpu_available`], [`probe_gpu_device`]) are
//! always compiled so callers can query GPU availability without `#[cfg]` gates. When `femcore`
//! is compiled with `--features gpu`, [`GpuContext`], [`GpuVecF32`], [`GpuSellCSigmaF32`],
//! [`GpuAmgF32`], and (when `SHADER_F64` is exposed by the adapter) [`GpuVecF64`] and
//! [`GpuSellCSigmaF64`] provide pure-Rust `wgpu` WGSL compute kernels:
//!
//! * Level-1 vector operations (`vec_axpy`, `vec_aypx`, `vec_scale`),
//! * Atomic-free fixed-order two-level binary-tree reductions (`BLOCK = 256`) for `dot` and `norm2`,
//! * Warp-coalesced `SELL-32-128` SpMV (`GpuSellCSigmaF32` / `GpuSellCSigmaF64`), and
//! * Device-resident Chebyshev `f32` AMG V-cycles (`GpuAmgF32`) driving outer `f64` PCG / FGMRES.
//!
//! **Future GPU & HPC Roadmap:** `wgpu`/WGSL serves as the zero-SDK portable Tier-1 compute
//! layer today and may be replaced by `CubeCL` or `rust-gpu` in the future once `f64` coverage
//! and stable toolchain compatibility mature (`femcore::linalg::device` keeps all Krylov and
//! eigensolver logic backend-agnostic). Dedicated HPC support and direct native bindings for
//! NVIDIA (CUDA) and AMD (ROCm/HIP) GPUs are planned for the future; if a user requires these on
//! a priority basis, they should contact the maintainers with target hardware and workload details.

use std::fmt;

#[cfg(feature = "gpu")]
use crate::linalg::amg::{AmgF32, AmgOptions, AmgStats, SmootherKind, f32_safe_scale};
#[cfg(feature = "gpu")]
use crate::linalg::csr::Csr;
#[cfg(feature = "gpu")]
use crate::linalg::device::{
    CpuVec, DEFAULT_SELL_C, DEFAULT_SELL_SIGMA, DeviceLinearOperator, DevicePreconditioner, SellCSigma, Vector,
};
#[cfg(feature = "gpu")]
use crate::linalg::linop::Preconditioner;
#[cfg(feature = "gpu")]
use bytemuck::{Pod, Zeroable};
#[cfg(feature = "gpu")]
use std::sync::{Arc, Mutex, OnceLock, mpsc};
#[cfg(feature = "gpu")]
use wgpu::util::DeviceExt;

/// Structured errors returned by GPU adapter discovery, buffer transfers, and compute dispatch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GpuError {
    /// `femcore` was built without `--features gpu`.
    FeatureDisabled,
    /// No compatible `wgpu` compute adapter was found on the host.
    NoAdapter,
    /// Requesting a logical device from the adapter failed.
    RequestDevice(String),
    /// Asynchronous staging-buffer map failed or channel closed.
    BufferMap(String),
    /// Vector or matrix dimensions did not match.
    DimensionMismatch { expected: usize, got: usize },
    /// Invalid input argument (such as non-finite matrix entries or non-32 slice height).
    InvalidArgument(&'static str),
    /// Underlying AMG hierarchy construction error.
    Solver(String),
}

impl fmt::Display for GpuError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FeatureDisabled => write!(f, "GPU backend not compiled; rebuild with `--features gpu`"),
            Self::NoAdapter => write!(f, "no usable wgpu compute adapter found"),
            Self::RequestDevice(e) => write!(f, "failed to create wgpu device: {e}"),
            Self::BufferMap(e) => write!(f, "failed to map GPU staging buffer: {e}"),
            Self::DimensionMismatch { expected, got } => {
                write!(f, "GPU dimension mismatch: expected {expected}, got {got}")
            }
            Self::InvalidArgument(m) => write!(f, "invalid GPU argument: {m}"),
            Self::Solver(e) => write!(f, "GPU solver setup error: {e}"),
        }
    }
}

impl std::error::Error for GpuError {}

/// Hardware and backend description of a selected `wgpu` compute adapter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GpuDeviceInfo {
    pub name: String,
    pub vendor: u32,
    pub device: u32,
    pub device_type: String,
    pub backend: String,
    pub is_discrete: bool,
    pub supports_f64: bool,
}

/// Returns `true` when `femcore` is compiled with `--features gpu` and a usable GPU compute
/// adapter can be initialised on this machine.
pub fn gpu_available() -> bool {
    probe_gpu_device().is_some()
}

/// Probes the highest-priority `wgpu` adapter (preferring discrete GPUs) and returns its
/// [`GpuDeviceInfo`], or `None` when `feature = "gpu"` is disabled or no adapter exists.
#[cfg(not(feature = "gpu"))]
pub fn probe_gpu_device() -> Option<GpuDeviceInfo> {
    None
}

#[cfg(feature = "gpu")]
pub fn probe_gpu_device() -> Option<GpuDeviceInfo> {
    GpuContext::shared().ok().map(|ctx| ctx.info().clone())
}

// ---------------------------------------------------------------------------------------
// wgpu Compute Context, Shaders, Vectors, SELL-C-sigma SpMV, and Chebyshev AMG
// ---------------------------------------------------------------------------------------

#[cfg(feature = "gpu")]
const MAX_WG_PER_DIM: u32 = 65_535;

#[cfg(feature = "gpu")]
const LINALG_WGSL: &str = r#"
struct VecParams {
    n: u32,
    wg_stride: u32,
    alpha: f32,
    beta: f32,
};

struct SellParams {
    nrows: u32,
    n_slices: u32,
    wg_stride: u32,
    mode: u32,
};

struct ChebParams {
    n: u32,
    wg_stride: u32,
    mode: u32,
    _pad0: u32,
    alpha: f32,
    beta: f32,
    _pad1: f32,
    _pad2: f32,
};

// ---- Unary scale: y[i] = alpha * y[i] ----
@group(0) @binding(0) var<uniform> u_params: VecParams;
@group(0) @binding(1) var<storage, read_write> u_y: array<f32>;

@compute @workgroup_size(256)
fn vec_scale(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * u_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < u_params.n) {
        u_y[idx] = u_params.alpha * u_y[idx];
    }
}

// ---- Binary axpy / aypx ----
@group(0) @binding(0) var<uniform> b_params: VecParams;
@group(0) @binding(1) var<storage, read> b_x: array<f32>;
@group(0) @binding(2) var<storage, read_write> b_y: array<f32>;

@compute @workgroup_size(256)
fn vec_axpy(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * b_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < b_params.n) {
        b_y[idx] = b_y[idx] + b_params.alpha * b_x[idx];
    }
}

@compute @workgroup_size(256)
fn vec_aypx(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * b_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < b_params.n) {
        b_y[idx] = b_x[idx] + b_params.beta * b_y[idx];
    }
}

// ---- Atomic-free fixed-order binary-tree dot product (Level 1, BLOCK = 256) ----
@group(0) @binding(0) var<uniform> d_params: VecParams;
@group(0) @binding(1) var<storage, read> d_x: array<f32>;
@group(0) @binding(2) var<storage, read> d_y: array<f32>;
@group(0) @binding(3) var<storage, read_write> d_partials: array<f32>;

var<workgroup> sdata_dot: array<f32, 256>;

@compute @workgroup_size(256)
fn tree_reduce_dot_l1(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * d_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < d_params.n) {
        sdata_dot[tid] = d_x[idx] * d_y[idx];
    } else {
        sdata_dot[tid] = 0.0;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_dot[tid] = sdata_dot[tid] + sdata_dot[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        d_partials[wg] = sdata_dot[0];
    }
}

// ---- Atomic-free single-buffer norm2_sq (Level 1) & multi-level sum reduction ----
@group(0) @binding(0) var<uniform> r_params: VecParams;
@group(0) @binding(1) var<storage, read> r_in: array<f32>;
@group(0) @binding(2) var<storage, read_write> r_out: array<f32>;

var<workgroup> sdata_red: array<f32, 256>;

@compute @workgroup_size(256)
fn tree_reduce_norm2_l1(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * r_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < r_params.n) {
        let v = r_in[idx];
        sdata_red[tid] = v * v;
    } else {
        sdata_red[tid] = 0.0;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_red[tid] = sdata_red[tid] + sdata_red[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        r_out[wg] = sdata_red[0];
    }
}

@compute @workgroup_size(256)
fn tree_reduce_sum_level(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * r_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < r_params.n) {
        sdata_red[tid] = r_in[idx];
    } else {
        sdata_red[tid] = 0.0;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_red[tid] = sdata_red[tid] + sdata_red[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        r_out[wg] = sdata_red[0];
    }
}

// ---- Warp-coalesced SELL-32-sigma SpMV ----
// mode 0: y = A * x
// mode 1: y = rhs - A * x
// mode 2: y = y + A * x
@group(0) @binding(0) var<uniform> s_params: SellParams;
@group(0) @binding(1) var<storage, read> s_slice_ptr: array<u32>;
@group(0) @binding(2) var<storage, read> s_col_idx: array<u32>;
@group(0) @binding(3) var<storage, read> s_values: array<f32>;
@group(0) @binding(4) var<storage, read> s_row_lengths: array<u32>;
@group(0) @binding(5) var<storage, read> s_perm: array<u32>;
@group(0) @binding(6) var<storage, read> s_x: array<f32>;
@group(0) @binding(7) var<storage, read> s_rhs: array<f32>;
@group(0) @binding(8) var<storage, read_write> s_y: array<f32>;

@compute @workgroup_size(256)
fn sell_spmv(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let gtid = (wid.y * s_params.wg_stride + wid.x) * 256u + lid.x;
    let slice_id = gtid >> 5u;
    let lane = gtid & 31u;
    if (slice_id >= s_params.n_slices) {
        return;
    }
    let sorted_row = slice_id * 32u + lane;
    if (sorted_row >= s_params.nrows) {
        return;
    }
    let base = s_slice_ptr[slice_id];
    let next_base = s_slice_ptr[slice_id + 1u];
    let width = (next_base - base) >> 5u;
    let rlen = s_row_lengths[sorted_row];
    var acc: f32 = 0.0;
    for (var k: u32 = 0u; k < width; k = k + 1u) {
        if (k < rlen) {
            let pos = base + k * 32u + lane;
            let col = s_col_idx[pos];
            acc = acc + s_values[pos] * s_x[col];
        }
    }
    let orig_row = s_perm[sorted_row];
    if (s_params.mode == 0u) {
        s_y[orig_row] = acc;
    } else if (s_params.mode == 1u) {
        s_y[orig_row] = s_rhs[orig_row] - acc;
    } else {
        s_y[orig_row] = s_y[orig_row] + acc;
    }
}

// ---- Chebyshev smoother vector update ----
// mode 0: step 0 with zero initial guess (d = alpha * inv_diag * res; x = d)
// mode 1: step 0 with non-zero initial guess (d = alpha * inv_diag * res; x += d)
// mode 2: step >= 1 (d = alpha * d + beta * inv_diag * res; x += d)
@group(0) @binding(0) var<uniform> c_params: ChebParams;
@group(0) @binding(1) var<storage, read> c_inv_diag: array<f32>;
@group(0) @binding(2) var<storage, read> c_res: array<f32>;
@group(0) @binding(3) var<storage, read_write> c_d: array<f32>;
@group(0) @binding(4) var<storage, read_write> c_x: array<f32>;

@compute @workgroup_size(256)
fn chebyshev_step(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * c_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx >= c_params.n) {
        return;
    }
    if (c_params.mode == 0u) {
        let di = c_params.alpha * (c_inv_diag[idx] * c_res[idx]);
        c_d[idx] = di;
        c_x[idx] = di;
    } else if (c_params.mode == 1u) {
        let di = c_params.alpha * (c_inv_diag[idx] * c_res[idx]);
        c_d[idx] = di;
        c_x[idx] = c_x[idx] + di;
    } else {
        let di = c_params.alpha * c_d[idx] + c_params.beta * (c_inv_diag[idx] * c_res[idx]);
        c_d[idx] = di;
        c_x[idx] = c_x[idx] + di;
    }
}
"#;

#[cfg(feature = "gpu")]
const LINALG_F64_WGSL: &str = r#"
struct VecParamsF64 {
    alpha: f64,
    beta: f64,
    n: u32,
    wg_stride: u32,
    _pad0: u32,
    _pad1: u32,
};

struct SellParams {
    nrows: u32,
    n_slices: u32,
    wg_stride: u32,
    mode: u32,
};

// ---- Unary scale: y[i] = alpha * y[i] ----
@group(0) @binding(0) var<uniform> u_params: VecParamsF64;
@group(0) @binding(1) var<storage, read_write> u_y: array<f64>;

@compute @workgroup_size(256)
fn vec_scale_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * u_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < u_params.n) {
        u_y[idx] = u_params.alpha * u_y[idx];
    }
}

// ---- Binary axpy / aypx ----
@group(0) @binding(0) var<uniform> b_params: VecParamsF64;
@group(0) @binding(1) var<storage, read> b_x: array<f64>;
@group(0) @binding(2) var<storage, read_write> b_y: array<f64>;

@compute @workgroup_size(256)
fn vec_axpy_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * b_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < b_params.n) {
        b_y[idx] = b_y[idx] + b_params.alpha * b_x[idx];
    }
}

@compute @workgroup_size(256)
fn vec_aypx_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * b_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < b_params.n) {
        b_y[idx] = b_x[idx] + b_params.beta * b_y[idx];
    }
}

// ---- Mixed-precision device casts for GpuAmgF32 <-> GpuVecF64 ----
@group(0) @binding(0) var<uniform> c64_params: VecParamsF64;
@group(0) @binding(1) var<storage, read> c_src_f64: array<f64>;
@group(0) @binding(2) var<storage, read_write> c_dst_f32: array<f32>;

@compute @workgroup_size(256)
fn cast_f64_to_f32(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * c64_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < c64_params.n) {
        c_dst_f32[idx] = f32(c_src_f64[idx] * c64_params.alpha);
    }
}

@group(0) @binding(0) var<uniform> c32_params: VecParamsF64;
@group(0) @binding(1) var<storage, read> c_src_f32: array<f32>;
@group(0) @binding(2) var<storage, read_write> c_dst_f64: array<f64>;

@compute @workgroup_size(256)
fn cast_f32_to_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let idx = (wid.y * c32_params.wg_stride + wid.x) * 256u + lid.x;
    if (idx < c32_params.n) {
        c_dst_f64[idx] = f64(c_src_f32[idx]) * c32_params.alpha;
    }
}

// ---- Atomic-free fixed-order binary-tree dot product (f64, BLOCK = 256) ----
@group(0) @binding(0) var<uniform> d_params: VecParamsF64;
@group(0) @binding(1) var<storage, read> d_x: array<f64>;
@group(0) @binding(2) var<storage, read> d_y: array<f64>;
@group(0) @binding(3) var<storage, read_write> d_partials: array<f64>;

var<workgroup> sdata_dot_f64: array<f64, 256>;

@compute @workgroup_size(256)
fn tree_reduce_dot_l1_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * d_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < d_params.n) {
        sdata_dot_f64[tid] = d_x[idx] * d_y[idx];
    } else {
        sdata_dot_f64[tid] = 0.0lf;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_dot_f64[tid] = sdata_dot_f64[tid] + sdata_dot_f64[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        d_partials[wg] = sdata_dot_f64[0];
    }
}

// ---- Atomic-free single-buffer norm2_sq & sum reduction (f64, BLOCK = 256) ----
@group(0) @binding(0) var<uniform> r_params: VecParamsF64;
@group(0) @binding(1) var<storage, read> r_in: array<f64>;
@group(0) @binding(2) var<storage, read_write> r_out: array<f64>;

var<workgroup> sdata_red_f64: array<f64, 256>;

@compute @workgroup_size(256)
fn tree_reduce_norm2_l1_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * r_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < r_params.n) {
        let v = r_in[idx];
        sdata_red_f64[tid] = v * v;
    } else {
        sdata_red_f64[tid] = 0.0lf;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_red_f64[tid] = sdata_red_f64[tid] + sdata_red_f64[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        r_out[wg] = sdata_red_f64[0];
    }
}

@compute @workgroup_size(256)
fn tree_reduce_sum_level_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let wg = wid.y * r_params.wg_stride + wid.x;
    let tid = lid.x;
    let idx = wg * 256u + tid;
    if (idx < r_params.n) {
        sdata_red_f64[tid] = r_in[idx];
    } else {
        sdata_red_f64[tid] = 0.0lf;
    }
    workgroupBarrier();
    var s: u32 = 128u;
    loop {
        if (s == 0u) {
            break;
        }
        if (tid < s) {
            sdata_red_f64[tid] = sdata_red_f64[tid] + sdata_red_f64[tid + s];
        }
        workgroupBarrier();
        s = s >> 1u;
    }
    if (tid == 0u) {
        r_out[wg] = sdata_red_f64[0];
    }
}

// ---- Warp-coalesced SELL-32-sigma SpMV (f64) ----
@group(0) @binding(0) var<uniform> s_params: SellParams;
@group(0) @binding(1) var<storage, read> s_slice_ptr: array<u32>;
@group(0) @binding(2) var<storage, read> s_col_idx: array<u32>;
@group(0) @binding(3) var<storage, read> s_values: array<f64>;
@group(0) @binding(4) var<storage, read> s_row_lengths: array<u32>;
@group(0) @binding(5) var<storage, read> s_perm: array<u32>;
@group(0) @binding(6) var<storage, read> s_x: array<f64>;
@group(0) @binding(7) var<storage, read> s_rhs: array<f64>;
@group(0) @binding(8) var<storage, read_write> s_y: array<f64>;

@compute @workgroup_size(256)
fn sell_spmv_f64(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
) {
    let gtid = (wid.y * s_params.wg_stride + wid.x) * 256u + lid.x;
    let slice_id = gtid >> 5u;
    let lane = gtid & 31u;
    if (slice_id >= s_params.n_slices) {
        return;
    }
    let sorted_row = slice_id * 32u + lane;
    if (sorted_row >= s_params.nrows) {
        return;
    }
    let base = s_slice_ptr[slice_id];
    let next_base = s_slice_ptr[slice_id + 1u];
    let width = (next_base - base) >> 5u;
    let rlen = s_row_lengths[sorted_row];
    var acc: f64 = 0.0lf;
    for (var k: u32 = 0u; k < width; k = k + 1u) {
        if (k < rlen) {
            let pos = base + k * 32u + lane;
            let col = s_col_idx[pos];
            acc = acc + s_values[pos] * s_x[col];
        }
    }
    let orig_row = s_perm[sorted_row];
    if (s_params.mode == 0u) {
        s_y[orig_row] = acc;
    } else if (s_params.mode == 1u) {
        s_y[orig_row] = s_rhs[orig_row] - acc;
    } else {
        s_y[orig_row] = s_y[orig_row] + acc;
    }
}
"#;

#[cfg(feature = "gpu")]
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct VecParams {
    n: u32,
    wg_stride: u32,
    alpha: f32,
    beta: f32,
}

#[cfg(feature = "gpu")]
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct VecParamsF64 {
    alpha: f64,
    beta: f64,
    n: u32,
    wg_stride: u32,
    _pad0: u32,
    _pad1: u32,
}

#[cfg(feature = "gpu")]
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct SellParams {
    nrows: u32,
    n_slices: u32,
    wg_stride: u32,
    mode: u32,
}

#[cfg(feature = "gpu")]
#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct ChebParams {
    n: u32,
    wg_stride: u32,
    mode: u32,
    _pad0: u32,
    alpha: f32,
    beta: f32,
    _pad1: f32,
    _pad2: f32,
}

#[cfg(feature = "gpu")]
#[inline]
fn workgroup_grid_256(n_items: usize) -> (u32, u32, u32) {
    let total_wg = n_items.div_ceil(256).max(1);
    if total_wg <= MAX_WG_PER_DIM as usize {
        let g = total_wg as u32;
        (g, 1, g)
    } else {
        let gx = MAX_WG_PER_DIM;
        let gy = total_wg.div_ceil(MAX_WG_PER_DIM as usize) as u32;
        (gx, gy, gx)
    }
}

#[cfg(feature = "gpu")]
struct ReductionScratch {
    cap_blocks: usize,
    buf_a: wgpu::Buffer,
    buf_b: wgpu::Buffer,
    readback: wgpu::Buffer,
}

#[cfg(feature = "gpu")]
struct GpuF64Pipelines {
    pipe_scale: wgpu::ComputePipeline,
    pipe_axpy: wgpu::ComputePipeline,
    pipe_aypx: wgpu::ComputePipeline,
    pipe_dot_l1: wgpu::ComputePipeline,
    pipe_norm2_l1: wgpu::ComputePipeline,
    pipe_sum_level: wgpu::ComputePipeline,
    pipe_sell_spmv: wgpu::ComputePipeline,
    pipe_cast_64_to_32: wgpu::ComputePipeline,
    pipe_cast_32_to_64: wgpu::ComputePipeline,
}

#[cfg(feature = "gpu")]
struct GpuContextInner {
    device: wgpu::Device,
    queue: wgpu::Queue,
    info: GpuDeviceInfo,
    dummy_f32: wgpu::Buffer,
    dummy_f64: wgpu::Buffer,
    unary_bgl: wgpu::BindGroupLayout,
    binary_bgl: wgpu::BindGroupLayout,
    reduce_dot_bgl: wgpu::BindGroupLayout,
    reduce_single_bgl: wgpu::BindGroupLayout,
    sell_bgl: wgpu::BindGroupLayout,
    cheb_bgl: wgpu::BindGroupLayout,
    pipe_scale: wgpu::ComputePipeline,
    pipe_axpy: wgpu::ComputePipeline,
    pipe_aypx: wgpu::ComputePipeline,
    pipe_dot_l1: wgpu::ComputePipeline,
    pipe_norm2_l1: wgpu::ComputePipeline,
    pipe_sum_level: wgpu::ComputePipeline,
    pipe_sell_spmv: wgpu::ComputePipeline,
    pipe_cheb: wgpu::ComputePipeline,
    f64_pipes: Option<GpuF64Pipelines>,
    reduction_scratch: Mutex<ReductionScratch>,
}

/// Shared `wgpu` compute context owning the selected device, command queue, and pre-compiled
/// WGSL linear-algebra pipelines. Cheaply cloneable (`Arc`-backed).
#[cfg(feature = "gpu")]
#[derive(Clone)]
pub struct GpuContext {
    inner: Arc<GpuContextInner>,
}

#[cfg(feature = "gpu")]
fn adapter_priority(info: &wgpu::AdapterInfo) -> u32 {
    let dev_score = match info.device_type {
        wgpu::DeviceType::DiscreteGpu => 400,
        wgpu::DeviceType::IntegratedGpu => 300,
        wgpu::DeviceType::VirtualGpu => 200,
        wgpu::DeviceType::Other => 100,
        wgpu::DeviceType::Cpu => 10,
    };
    let backend_score = match info.backend {
        wgpu::Backend::Vulkan => 40,
        wgpu::Backend::Metal => 35,
        wgpu::Backend::Dx12 => 30,
        wgpu::Backend::Gl => 10,
        wgpu::Backend::BrowserWebGpu => 5,
        wgpu::Backend::Noop => 0,
    };
    dev_score + backend_score
}

#[cfg(feature = "gpu")]
fn format_adapter_info(info: &wgpu::AdapterInfo, supports_f64: bool) -> GpuDeviceInfo {
    GpuDeviceInfo {
        name: info.name.clone(),
        vendor: info.vendor,
        device: info.device,
        device_type: format!("{:?}", info.device_type),
        backend: format!("{:?}", info.backend),
        is_discrete: info.device_type == wgpu::DeviceType::DiscreteGpu,
        supports_f64,
    }
}

#[cfg(feature = "gpu")]
fn bgl_uniform(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

#[cfg(feature = "gpu")]
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

#[cfg(feature = "gpu")]
impl GpuContext {
    /// Initialises a new `GpuContext`, preferring high-performance discrete GPUs (such as an
    /// NVIDIA discrete card over an integrated GPU) and returning [`GpuError::NoAdapter`] if no
    /// usable compute adapter is present.
    pub fn new() -> Result<Self, GpuError> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let mut adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::PRIMARY));
        if adapters.is_empty() {
            adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()));
        }
        let chosen = if let Some(best) = adapters
            .into_iter()
            .filter(|a| a.get_info().backend != wgpu::Backend::Noop)
            .max_by_key(|a| adapter_priority(&a.get_info()))
        {
            best
        } else {
            pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
                ..Default::default()
            }))
            .map_err(|_| GpuError::NoAdapter)?
        };

        let supports_f64 = chosen.features().contains(wgpu::Features::SHADER_F64);
        let info = format_adapter_info(&chosen.get_info(), supports_f64);
        let limits = chosen.limits();
        let required_features = if supports_f64 { wgpu::Features::SHADER_F64 } else { wgpu::Features::empty() };
        let (device, queue) = pollster::block_on(chosen.request_device(&wgpu::DeviceDescriptor {
            label: Some("nyx-shared-gpu"),
            required_features,
            required_limits: limits,
            ..Default::default()
        }))
        .map_err(|e| GpuError::RequestDevice(e.to_string()))?;

        Ok(Self::from_wgpu(device, queue, info))
    }

    /// Returns a process-wide cached [`GpuContext`] so repeated callers share one compiled shader
    /// module and command queue.
    pub fn shared() -> Result<Self, GpuError> {
        static GLOBAL: OnceLock<Result<GpuContext, GpuError>> = OnceLock::new();
        GLOBAL.get_or_init(Self::new).clone()
    }

    /// Wraps an existing `wgpu::Device` and `wgpu::Queue` (e.g. shared across crates) and compiles
    /// the linear-algebra WGSL pipelines on it.
    pub fn from_wgpu(device: wgpu::Device, queue: wgpu::Queue, mut info: GpuDeviceInfo) -> Self {
        let has_f64 = device.features().contains(wgpu::Features::SHADER_F64);
        info.supports_f64 = has_f64;

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("nyx-femcore-linalg-wgsl"),
            source: wgpu::ShaderSource::Wgsl(LINALG_WGSL.into()),
        });

        let unary_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-unary-bgl"),
            entries: &[bgl_uniform(0), bgl_storage(1, false)],
        });
        let binary_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-binary-bgl"),
            entries: &[bgl_uniform(0), bgl_storage(1, true), bgl_storage(2, false)],
        });
        let reduce_dot_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-reduce-dot-bgl"),
            entries: &[bgl_uniform(0), bgl_storage(1, true), bgl_storage(2, true), bgl_storage(3, false)],
        });
        let reduce_single_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-reduce-single-bgl"),
            entries: &[bgl_uniform(0), bgl_storage(1, true), bgl_storage(2, false)],
        });
        let sell_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-sell-bgl"),
            entries: &[
                bgl_uniform(0),
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
        let cheb_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("linalg-cheb-bgl"),
            entries: &[
                bgl_uniform(0),
                bgl_storage(1, true),
                bgl_storage(2, true),
                bgl_storage(3, false),
                bgl_storage(4, false),
            ],
        });

        let make_pipe_from = |mod_ref: &wgpu::ShaderModule, label: &str, bgl: &wgpu::BindGroupLayout, entry: &str| {
            let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[Some(bgl)],
                immediate_size: 0,
            });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: Some(&layout),
                module: mod_ref,
                entry_point: Some(entry),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            })
        };

        let pipe_scale = make_pipe_from(&shader, "pipe-vec-scale", &unary_bgl, "vec_scale");
        let pipe_axpy = make_pipe_from(&shader, "pipe-vec-axpy", &binary_bgl, "vec_axpy");
        let pipe_aypx = make_pipe_from(&shader, "pipe-vec-aypx", &binary_bgl, "vec_aypx");
        let pipe_dot_l1 = make_pipe_from(&shader, "pipe-dot-l1", &reduce_dot_bgl, "tree_reduce_dot_l1");
        let pipe_norm2_l1 = make_pipe_from(&shader, "pipe-norm2-l1", &reduce_single_bgl, "tree_reduce_norm2_l1");
        let pipe_sum_level = make_pipe_from(&shader, "pipe-sum-level", &reduce_single_bgl, "tree_reduce_sum_level");
        let pipe_sell_spmv = make_pipe_from(&shader, "pipe-sell-spmv", &sell_bgl, "sell_spmv");
        let pipe_cheb = make_pipe_from(&shader, "pipe-cheb-step", &cheb_bgl, "chebyshev_step");

        let f64_pipes = if has_f64 {
            let shader_f64 = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("nyx-femcore-linalg-f64-wgsl"),
                source: wgpu::ShaderSource::Wgsl(LINALG_F64_WGSL.into()),
            });
            Some(GpuF64Pipelines {
                pipe_scale: make_pipe_from(&shader_f64, "pipe-vec-scale-f64", &unary_bgl, "vec_scale_f64"),
                pipe_axpy: make_pipe_from(&shader_f64, "pipe-vec-axpy-f64", &binary_bgl, "vec_axpy_f64"),
                pipe_aypx: make_pipe_from(&shader_f64, "pipe-vec-aypx-f64", &binary_bgl, "vec_aypx_f64"),
                pipe_dot_l1: make_pipe_from(&shader_f64, "pipe-dot-l1-f64", &reduce_dot_bgl, "tree_reduce_dot_l1_f64"),
                pipe_norm2_l1: make_pipe_from(
                    &shader_f64,
                    "pipe-norm2-l1-f64",
                    &reduce_single_bgl,
                    "tree_reduce_norm2_l1_f64",
                ),
                pipe_sum_level: make_pipe_from(
                    &shader_f64,
                    "pipe-sum-level-f64",
                    &reduce_single_bgl,
                    "tree_reduce_sum_level_f64",
                ),
                pipe_sell_spmv: make_pipe_from(&shader_f64, "pipe-sell-spmv-f64", &sell_bgl, "sell_spmv_f64"),
                pipe_cast_64_to_32: make_pipe_from(&shader_f64, "pipe-cast-64-to-32", &binary_bgl, "cast_f64_to_f32"),
                pipe_cast_32_to_64: make_pipe_from(&shader_f64, "pipe-cast-32-to-64", &binary_bgl, "cast_f32_to_f64"),
            })
        } else {
            None
        };

        let dummy_f32 = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("linalg-dummy-f32"),
            contents: bytemuck::bytes_of(&[0.0f32; 4]),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        });
        let dummy_f64 = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("linalg-dummy-f64"),
            contents: bytemuck::bytes_of(&[0.0f64; 4]),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        });

        let init_blocks = 4096usize;
        let alloc_partials = |label: &str, len: usize| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: (len.max(4) * std::mem::size_of::<f64>()) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        };
        let buf_a = alloc_partials("reduce-partials-a", init_blocks);
        let buf_b = alloc_partials("reduce-partials-b", init_blocks.div_ceil(256).max(4));
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("reduce-scalar-readback"),
            size: 16,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            inner: Arc::new(GpuContextInner {
                device,
                queue,
                info,
                dummy_f32,
                dummy_f64,
                unary_bgl,
                binary_bgl,
                reduce_dot_bgl,
                reduce_single_bgl,
                sell_bgl,
                cheb_bgl,
                pipe_scale,
                pipe_axpy,
                pipe_aypx,
                pipe_dot_l1,
                pipe_norm2_l1,
                pipe_sum_level,
                pipe_sell_spmv,
                pipe_cheb,
                f64_pipes,
                reduction_scratch: Mutex::new(ReductionScratch { cap_blocks: init_blocks, buf_a, buf_b, readback }),
            }),
        }
    }

    /// Adapter metadata for the active GPU device.
    #[inline]
    pub fn info(&self) -> &GpuDeviceInfo {
        &self.inner.info
    }

    /// True when this logical GPU device supports native IEEE-754 `f64` shaders (`SHADER_F64`).
    #[inline]
    pub fn supports_f64(&self) -> bool {
        self.inner.f64_pipes.is_some()
    }

    /// Underlying `wgpu::Device` handle for sharing across physics crates (e.g. `nyx-ems`).
    #[inline]
    pub fn wgpu_device(&self) -> &wgpu::Device {
        &self.inner.device
    }

    /// Underlying `wgpu::Queue` handle for sharing across physics crates (e.g. `nyx-ems`).
    #[inline]
    pub fn wgpu_queue(&self) -> &wgpu::Queue {
        &self.inner.queue
    }

    pub(crate) fn read_buffer_f32(&self, src: &wgpu::Buffer, len: usize) -> Result<Vec<f32>, GpuError> {
        if len == 0 {
            return Ok(Vec::new());
        }
        let byte_len = (len * std::mem::size_of::<f32>()) as u64;
        let staging = self.inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gpu-vec-readback"),
            size: byte_len,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = self.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("gpu-vec-readback-enc"),
        });
        encoder.copy_buffer_to_buffer(src, 0, &staging, 0, byte_len);
        self.inner.queue.submit(Some(encoder.finish()));

        let slice = staging.slice(0..byte_len);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        let _ = self.inner.device.poll(wgpu::PollType::wait_indefinitely());
        let map_res = rx.recv().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        map_res.map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let out: Vec<f32> = bytemuck::cast_slice(&mapped).to_vec();
        drop(mapped);
        staging.unmap();
        Ok(out)
    }

    fn reduce_two_buffers(&self, x: &wgpu::Buffer, y: Option<&wgpu::Buffer>, n: usize) -> Result<f32, GpuError> {
        if n == 0 {
            return Ok(0.0);
        }
        let inner = &*self.inner;
        let nblocks = n.div_ceil(256);
        let mut guard = inner
            .reduction_scratch
            .lock()
            .map_err(|_| GpuError::BufferMap("reduction scratch mutex poisoned".to_string()))?;

        if nblocks > guard.cap_blocks {
            let new_cap = nblocks.next_power_of_two().max(256);
            guard.buf_a = inner.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("reduce-partials-a"),
                size: (new_cap * 8) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            guard.buf_b = inner.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("reduce-partials-b"),
                size: (new_cap.div_ceil(256).max(4) * 8) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            guard.cap_blocks = new_cap;
        }

        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("tree-reduce-enc"),
        });

        let (gx, gy, wg_stride) = workgroup_grid_256(n);
        let p_l1 = VecParams { n: n as u32, wg_stride, alpha: 1.0, beta: 0.0 };
        let u_l1 = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("reduce-l1-uniform"),
            contents: bytemuck::bytes_of(&p_l1),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        if let Some(yb) = y {
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-dot-l1-bg"),
                layout: &inner.reduce_dot_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_l1.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: x.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: yb.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 3, resource: guard.buf_a.as_entire_binding() },
                ],
            });
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("reduce-dot-l1-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&inner.pipe_dot_l1);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        } else {
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-norm2-l1-bg"),
                layout: &inner.reduce_single_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_l1.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: x.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: guard.buf_a.as_entire_binding() },
                ],
            });
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("reduce-norm2-l1-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&inner.pipe_norm2_l1);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }

        let mut cur_len = nblocks;
        let mut ping = true;
        while cur_len > 1 {
            let (lgx, lgy, lstride) = workgroup_grid_256(cur_len);
            let p_lvl = VecParams { n: cur_len as u32, wg_stride: lstride, alpha: 1.0, beta: 0.0 };
            let u_lvl = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("reduce-lvl-uniform"),
                contents: bytemuck::bytes_of(&p_lvl),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let (src_buf, dst_buf) = if ping { (&guard.buf_a, &guard.buf_b) } else { (&guard.buf_b, &guard.buf_a) };
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-sum-lvl-bg"),
                layout: &inner.reduce_single_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_lvl.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: src_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: dst_buf.as_entire_binding() },
                ],
            });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("reduce-sum-lvl-pass"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_sum_level);
                pass.set_bind_group(0, &bg, &[]);
                pass.dispatch_workgroups(lgx, lgy, 1);
            }
            cur_len = cur_len.div_ceil(256);
            ping = !ping;
        }

        let final_buf = if ping { &guard.buf_a } else { &guard.buf_b };
        encoder.copy_buffer_to_buffer(final_buf, 0, &guard.readback, 0, 4);
        inner.queue.submit(Some(encoder.finish()));

        let slice = guard.readback.slice(0..4);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        let _ = inner.device.poll(wgpu::PollType::wait_indefinitely());
        let map_res = rx.recv().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        map_res.map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let val = *bytemuck::from_bytes::<f32>(&mapped[0..4]);
        drop(mapped);
        guard.readback.unmap();
        Ok(val)
    }

    pub(crate) fn read_buffer_f64(&self, src: &wgpu::Buffer, len: usize) -> Result<Vec<f64>, GpuError> {
        if len == 0 {
            return Ok(Vec::new());
        }
        let byte_len = (len * std::mem::size_of::<f64>()) as u64;
        let staging = self.inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gpu-vec-f64-readback"),
            size: byte_len,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = self.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("gpu-vec-f64-readback-enc"),
        });
        encoder.copy_buffer_to_buffer(src, 0, &staging, 0, byte_len);
        self.inner.queue.submit(Some(encoder.finish()));

        let slice = staging.slice(0..byte_len);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        let _ = self.inner.device.poll(wgpu::PollType::wait_indefinitely());
        let map_res = rx.recv().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        map_res.map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let out: Vec<f64> = bytemuck::cast_slice(&mapped).to_vec();
        drop(mapped);
        staging.unmap();
        Ok(out)
    }

    fn reduce_two_buffers_f64(&self, x: &wgpu::Buffer, y: Option<&wgpu::Buffer>, n: usize) -> Result<f64, GpuError> {
        if n == 0 {
            return Ok(0.0);
        }
        let inner = &*self.inner;
        let f64p = inner
            .f64_pipes
            .as_ref()
            .ok_or(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"))?;
        let nblocks = n.div_ceil(256);
        let mut guard = inner
            .reduction_scratch
            .lock()
            .map_err(|_| GpuError::BufferMap("reduction scratch mutex poisoned".to_string()))?;

        if nblocks > guard.cap_blocks {
            let new_cap = nblocks.next_power_of_two().max(256);
            guard.buf_a = inner.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("reduce-partials-a"),
                size: (new_cap * 8) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            guard.buf_b = inner.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("reduce-partials-b"),
                size: (new_cap.div_ceil(256).max(4) * 8) as u64,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            guard.cap_blocks = new_cap;
        }

        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("tree-reduce-f64-enc"),
        });

        let (gx, gy, wg_stride) = workgroup_grid_256(n);
        let p_l1 = VecParamsF64 { alpha: 1.0, beta: 0.0, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
        let u_l1 = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("reduce-l1-f64-uniform"),
            contents: bytemuck::bytes_of(&p_l1),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        if let Some(yb) = y {
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-dot-l1-f64-bg"),
                layout: &inner.reduce_dot_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_l1.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: x.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: yb.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 3, resource: guard.buf_a.as_entire_binding() },
                ],
            });
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("reduce-dot-l1-f64-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_dot_l1);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        } else {
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-norm2-l1-f64-bg"),
                layout: &inner.reduce_single_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_l1.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: x.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: guard.buf_a.as_entire_binding() },
                ],
            });
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("reduce-norm2-l1-f64-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_norm2_l1);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }

        let mut cur_len = nblocks;
        let mut ping = true;
        while cur_len > 1 {
            let (lgx, lgy, lstride) = workgroup_grid_256(cur_len);
            let p_lvl = VecParamsF64 { alpha: 1.0, beta: 0.0, n: cur_len as u32, wg_stride: lstride, _pad0: 0, _pad1: 0 };
            let u_lvl = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("reduce-lvl-f64-uniform"),
                contents: bytemuck::bytes_of(&p_lvl),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let (src_buf, dst_buf) = if ping { (&guard.buf_a, &guard.buf_b) } else { (&guard.buf_b, &guard.buf_a) };
            let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("reduce-sum-lvl-f64-bg"),
                layout: &inner.reduce_single_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_lvl.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: src_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: dst_buf.as_entire_binding() },
                ],
            });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("reduce-sum-lvl-f64-pass"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&f64p.pipe_sum_level);
                pass.set_bind_group(0, &bg, &[]);
                pass.dispatch_workgroups(lgx, lgy, 1);
            }
            cur_len = cur_len.div_ceil(256);
            ping = !ping;
        }

        let final_buf = if ping { &guard.buf_a } else { &guard.buf_b };
        encoder.copy_buffer_to_buffer(final_buf, 0, &guard.readback, 0, 8);
        inner.queue.submit(Some(encoder.finish()));

        let slice = guard.readback.slice(0..8);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        let _ = inner.device.poll(wgpu::PollType::wait_indefinitely());
        let map_res = rx.recv().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        map_res.map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let val = *bytemuck::from_bytes::<f64>(&mapped[0..8]);
        drop(mapped);
        guard.readback.unmap();
        Ok(val)
    }
}

/// Device-resident single-precision `f32` vector stored in a `wgpu::Buffer`.
///
/// Implements [`Vector<Scalar = f64>`] so generic Krylov (`cg_vec`, `fgmres_vec`) and
/// eigensolver (`lobpcg_vec`) routines can operate directly on GPU-resident vectors.
#[cfg(feature = "gpu")]
pub struct GpuVecF32 {
    ctx: GpuContext,
    buf: wgpu::Buffer,
    len: usize,
}

#[cfg(feature = "gpu")]
impl fmt::Debug for GpuVecF32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GpuVecF32").field("len", &self.len).field("device", &self.ctx.info().name).finish()
    }
}

#[cfg(feature = "gpu")]
impl Clone for GpuVecF32 {
    fn clone(&self) -> Self {
        let byte_size = (self.len.max(4) * std::mem::size_of::<f32>()) as u64;
        let new_buf = self.ctx.inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gpu-vec-f32-clone"),
            size: byte_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        if self.len > 0 {
            let mut encoder = self.ctx.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("gpu-vec-clone-enc"),
            });
            encoder.copy_buffer_to_buffer(&self.buf, 0, &new_buf, 0, (self.len * 4) as u64);
            self.ctx.inner.queue.submit(Some(encoder.finish()));
        }
        Self { ctx: self.ctx.clone(), buf: new_buf, len: self.len }
    }
}

#[cfg(feature = "gpu")]
impl GpuVecF32 {
    /// Allocates a zero-initialised `GpuVecF32` of length `len` on `ctx`.
    pub fn zeros(ctx: &GpuContext, len: usize) -> Self {
        let alloc_elems = len.max(4);
        let zeros = vec![0.0f32; alloc_elems];
        let buf = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gpu-vec-f32"),
            contents: bytemuck::cast_slice(&zeros),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        });
        Self { ctx: ctx.clone(), buf, len }
    }

    /// Allocates and uploads a `GpuVecF32` from an `f32` host slice.
    pub fn from_f32_slice(ctx: &GpuContext, data: &[f32]) -> Self {
        if data.len() >= 4 {
            let buf = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gpu-vec-f32"),
                contents: bytemuck::cast_slice(data),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            });
            Self { ctx: ctx.clone(), buf, len: data.len() }
        } else {
            let mut pad = [0.0f32; 4];
            pad[..data.len()].copy_from_slice(data);
            let buf = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gpu-vec-f32"),
                contents: bytemuck::bytes_of(&pad),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            });
            Self { ctx: ctx.clone(), buf, len: data.len() }
        }
    }

    /// Allocates and uploads a `GpuVecF32` by narrowing an `f64` host slice to `f32`.
    pub fn from_f64_slice(ctx: &GpuContext, data: &[f64]) -> Self {
        let f32s: Vec<f32> = data.iter().map(|&v| v as f32).collect();
        Self::from_f32_slice(ctx, &f32s)
    }

    #[inline]
    pub fn context(&self) -> &GpuContext {
        &self.ctx
    }

    #[inline]
    pub fn buffer(&self) -> &wgpu::Buffer {
        &self.buf
    }

    /// Reads the GPU vector back to host memory as `Vec<f32>`.
    pub fn to_host_f32(&self) -> Result<Vec<f32>, GpuError> {
        self.ctx.read_buffer_f32(&self.buf, self.len)
    }

    /// Overwrites this GPU vector from an `f32` host slice of matching length.
    pub fn write_f32(&mut self, data: &[f32]) -> Result<(), GpuError> {
        if data.len() != self.len {
            return Err(GpuError::DimensionMismatch { expected: self.len, got: data.len() });
        }
        if !data.is_empty() {
            self.ctx.inner.queue.write_buffer(&self.buf, 0, bytemuck::cast_slice(data));
        }
        Ok(())
    }

    /// Computes the atomic-free two-level binary-tree dot product in `f32`.
    pub fn dot_f32(&self, other: &Self) -> Result<f32, GpuError> {
        if self.len != other.len {
            return Err(GpuError::DimensionMismatch { expected: self.len, got: other.len() });
        }
        if std::ptr::eq(&self.buf, &other.buf) {
            self.ctx.reduce_two_buffers(&self.buf, None, self.len)
        } else {
            self.ctx.reduce_two_buffers(&self.buf, Some(&other.buf), self.len)
        }
    }

    /// Computes the atomic-free two-level binary-tree squared 2-norm in `f32`.
    pub fn norm2_sq_f32(&self) -> Result<f32, GpuError> {
        self.ctx.reduce_two_buffers(&self.buf, None, self.len)
    }

    fn dispatch_binary(&mut self, other: &Self, alpha: f32, beta: f32, is_axpy: bool) {
        let n = self.len.min(other.len);
        if n == 0 {
            return;
        }
        let inner = &*self.ctx.inner;
        let (gx, gy, wg_stride) = workgroup_grid_256(n);
        let params = VecParams { n: n as u32, wg_stride, alpha, beta };
        let ubuf = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vec-bin-uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("vec-bin-bg"),
            layout: &inner.binary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: ubuf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: other.buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: self.buf.as_entire_binding() },
            ],
        });
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("vec-bin-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("vec-bin-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(if is_axpy { &inner.pipe_axpy } else { &inner.pipe_aypx });
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
    }

    fn dispatch_scale(&mut self, alpha: f32) {
        if self.len == 0 {
            return;
        }
        let inner = &*self.ctx.inner;
        let (gx, gy, wg_stride) = workgroup_grid_256(self.len);
        let params = VecParams { n: self.len as u32, wg_stride, alpha, beta: 0.0 };
        let ubuf = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vec-scale-uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("vec-scale-bg"),
            layout: &inner.unary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: ubuf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: self.buf.as_entire_binding() },
            ],
        });
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("vec-scale-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("vec-scale-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&inner.pipe_scale);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
    }
}

#[cfg(feature = "gpu")]
impl Vector for GpuVecF32 {
    type Scalar = f64;

    #[inline]
    fn len(&self) -> usize {
        self.len
    }

    fn dotc(&self, other: &Self) -> f64 {
        self.dot_f32(other).map_or(f64::NAN, |v| v as f64)
    }

    fn dotu(&self, other: &Self) -> f64 {
        self.dot_f32(other).map_or(f64::NAN, |v| v as f64)
    }

    fn norm2_sq(&self) -> f64 {
        self.norm2_sq_f32().map_or(f64::NAN, |v| v as f64)
    }

    fn axpy(&mut self, a: f64, x: &Self) {
        self.dispatch_binary(x, a as f32, 0.0, true);
    }

    fn scale(&mut self, a: f64) {
        self.dispatch_scale(a as f32);
    }

    fn scale_real(&mut self, r: f64) {
        self.dispatch_scale(r as f32);
    }

    fn aypx(&mut self, beta: f64, x: &Self) {
        self.dispatch_binary(x, 0.0, beta as f32, false);
    }

    fn copy_from(&mut self, other: &Self) {
        let n = self.len.min(other.len);
        if n == 0 {
            return;
        }
        let mut encoder = self.ctx.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("gpu-vec-copy-enc"),
        });
        encoder.copy_buffer_to_buffer(&other.buf, 0, &self.buf, 0, (n * 4) as u64);
        self.ctx.inner.queue.submit(Some(encoder.finish()));
    }

    fn to_host(&self) -> Vec<f64> {
        match self.to_host_f32() {
            Ok(v) => v.into_iter().map(|x| x as f64).collect(),
            Err(_) => vec![f64::NAN; self.len],
        }
    }

    fn from_host(&mut self, host: &[f64]) {
        let f32s: Vec<f32> = host.iter().map(|&v| v as f32).collect();
        if f32s.len() == self.len {
            let _ = self.write_f32(&f32s);
        } else {
            *self = Self::from_f32_slice(&self.ctx, &f32s);
        }
    }

    fn zeros_like(&self) -> Self {
        Self::zeros(&self.ctx, self.len)
    }
}

// ---------------------------------------------------------------------------------------
// GpuSellCSigmaF32: Warp-Coalesced SELL-32-128 Sparse Matrix on GPU
// ---------------------------------------------------------------------------------------

#[cfg(feature = "gpu")]
pub struct GpuSellCSigmaF32 {
    ctx: GpuContext,
    nrows: usize,
    ncols: usize,
    nnz: usize,
    n_slices: u32,
    slice_ptr_buf: wgpu::Buffer,
    col_idx_buf: wgpu::Buffer,
    values_buf: wgpu::Buffer,
    row_lengths_buf: wgpu::Buffer,
    perm_buf: wgpu::Buffer,
    uniform_spmv: wgpu::Buffer,
    gx: u32,
    gy: u32,
}

#[cfg(feature = "gpu")]
impl GpuSellCSigmaF32 {
    /// Converts a host [`Csr<f64>`] into `SELL-32-128` and uploads its structure and `f32` values
    /// to device storage buffers.
    pub fn from_csr(ctx: &GpuContext, csr: &Csr<f64>) -> Result<Self, GpuError> {
        if csr.values().iter().any(|v| !v.is_finite()) {
            return Err(GpuError::InvalidArgument("CSR matrix contains non-finite entries"));
        }
        let sell = SellCSigma::from_csr_with(csr, DEFAULT_SELL_C, DEFAULT_SELL_SIGMA)
            .map_err(|_| GpuError::InvalidArgument("failed to convert CSR to SELL-32-128"))?;
        Self::from_sell(ctx, &sell)
    }

    /// Uploads a host `SellCSigma<f64>` (`C = 32`) to device `f32` storage buffers.
    pub fn from_sell(ctx: &GpuContext, sell: &SellCSigma<f64>) -> Result<Self, GpuError> {
        if sell.c() != 32 {
            return Err(GpuError::InvalidArgument("GpuSellCSigmaF32 requires slice height C == 32"));
        }
        let mut sp_u32 = Vec::with_capacity(sell.slice_ptr().len());
        for &p in sell.slice_ptr() {
            sp_u32.push(u32::try_from(p).map_err(|_| GpuError::InvalidArgument("SELL slice_ptr overflows u32"))?);
        }
        let n_slices = (sell.slice_ptr().len().saturating_sub(1)) as u32;
        let vals_f32: Vec<f32> = sell.values().iter().map(|&v| v as f32).collect();

        let mk_u32 = |label: &str, data: &[u32]| {
            if data.is_empty() {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::bytes_of(&[0u32; 4]),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            } else {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::cast_slice(data),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            }
        };
        let mk_f32 = |label: &str, data: &[f32]| {
            if data.is_empty() {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::bytes_of(&[0.0f32; 4]),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            } else {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::cast_slice(data),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            }
        };

        let slice_ptr_buf = mk_u32("sell-slice-ptr", &sp_u32);
        let col_idx_buf = mk_u32("sell-col-idx", sell.col_idx());
        let values_buf = mk_f32("sell-values", &vals_f32);
        let row_lengths_buf = mk_u32("sell-row-lengths", sell.row_lengths());
        let perm_buf = mk_u32("sell-perm", sell.perm());

        let total_threads = (n_slices as usize) * 32;
        let (gx, gy, wg_stride) = workgroup_grid_256(total_threads);
        let sp_params = SellParams { nrows: sell.nrows() as u32, n_slices, wg_stride, mode: 0 };
        let uniform_spmv = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("sell-spmv-uniform"),
            contents: bytemuck::bytes_of(&sp_params),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        Ok(Self {
            ctx: ctx.clone(),
            nrows: sell.nrows(),
            ncols: sell.ncols(),
            nnz: sell.nnz(),
            n_slices,
            slice_ptr_buf,
            col_idx_buf,
            values_buf,
            row_lengths_buf,
            perm_buf,
            uniform_spmv,
            gx,
            gy,
        })
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    pub fn nnz(&self) -> usize {
        self.nnz
    }

    fn create_bind_group(&self, mode: u32, x_buf: &wgpu::Buffer, rhs_buf: &wgpu::Buffer, y_buf: &wgpu::Buffer) -> wgpu::BindGroup {
        let inner = &*self.ctx.inner;
        let ubuf = if mode == 0 {
            None
        } else {
            let total_threads = (self.n_slices as usize) * 32;
            let (_, _, wg_stride) = workgroup_grid_256(total_threads);
            let p = SellParams { nrows: self.nrows as u32, n_slices: self.n_slices, wg_stride, mode };
            Some(inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("sell-mode-uniform"),
                contents: bytemuck::bytes_of(&p),
                usage: wgpu::BufferUsages::UNIFORM,
            }))
        };
        let u_ref = ubuf.as_ref().unwrap_or(&self.uniform_spmv);
        inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("sell-spmv-bg"),
            layout: &inner.sell_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: u_ref.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: self.slice_ptr_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: self.col_idx_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: self.values_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: self.row_lengths_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: self.perm_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: x_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: rhs_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: y_buf.as_entire_binding() },
            ],
        })
    }

    /// Dispatches `y ← A x` on the GPU.
    pub fn spmv(&self, x: &GpuVecF32, y: &mut GpuVecF32) -> Result<(), GpuError> {
        if x.len() != self.ncols {
            return Err(GpuError::DimensionMismatch { expected: self.ncols, got: x.len() });
        }
        if y.len() != self.nrows {
            return Err(GpuError::DimensionMismatch { expected: self.nrows, got: y.len() });
        }
        if self.nrows == 0 {
            return Ok(());
        }
        let inner = &*self.ctx.inner;
        let bg = self.create_bind_group(0, &x.buf, &inner.dummy_f32, &y.buf);
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("sell-spmv-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("sell-spmv-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&inner.pipe_sell_spmv);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(self.gx, self.gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
        Ok(())
    }
}

#[cfg(feature = "gpu")]
impl DeviceLinearOperator<GpuVecF32> for GpuSellCSigmaF32 {
    #[inline]
    fn dim(&self) -> usize {
        self.nrows
    }
    fn apply(&self, x: &GpuVecF32, y: &mut GpuVecF32) {
        let _ = self.spmv(x, y);
    }
}

// ---------------------------------------------------------------------------------------
// GpuVecF64 & GpuSellCSigmaF64: Native IEEE-754 f64 Device Vectors and SpMV (SHADER_F64)
// ---------------------------------------------------------------------------------------

/// Device-resident double-precision `f64` vector stored in a `wgpu::Buffer` (`wgpu::Features::SHADER_F64`).
///
/// Implements [`Vector<Scalar = f64>`] so device-resident Krylov (`cg_vec`, `fgmres_vec`) and
/// eigensolver (`lobpcg_vec`) routines can converge to full `f64` accuracy (`< 1e-10`) on GPUs
/// that expose `SHADER_F64`, while pairing cleanly with `GpuAmgF32` as the fast inner V-cycle.
#[cfg(feature = "gpu")]
pub struct GpuVecF64 {
    ctx: GpuContext,
    buf: wgpu::Buffer,
    len: usize,
}

#[cfg(feature = "gpu")]
impl fmt::Debug for GpuVecF64 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("GpuVecF64").field("len", &self.len).field("device", &self.ctx.info().name).finish()
    }
}

#[cfg(feature = "gpu")]
impl Clone for GpuVecF64 {
    fn clone(&self) -> Self {
        let byte_size = (self.len.max(4) * std::mem::size_of::<f64>()) as u64;
        let new_buf = self.ctx.inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("gpu-vec-f64-clone"),
            size: byte_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        if self.len > 0 {
            let mut encoder = self.ctx.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("gpu-vec-f64-clone-enc"),
            });
            encoder.copy_buffer_to_buffer(&self.buf, 0, &new_buf, 0, (self.len * 8) as u64);
            self.ctx.inner.queue.submit(Some(encoder.finish()));
        }
        Self { ctx: self.ctx.clone(), buf: new_buf, len: self.len }
    }
}

#[cfg(feature = "gpu")]
impl GpuVecF64 {
    fn require_f64(ctx: &GpuContext) -> Result<(), GpuError> {
        if ctx.supports_f64() {
            Ok(())
        } else {
            Err(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"))
        }
    }

    /// Allocates a zero-initialised `GpuVecF64` of length `len` on `ctx`.
    pub fn zeros(ctx: &GpuContext, len: usize) -> Result<Self, GpuError> {
        Self::require_f64(ctx)?;
        let alloc_elems = len.max(4);
        let zeros = vec![0.0f64; alloc_elems];
        let buf = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gpu-vec-f64"),
            contents: bytemuck::cast_slice(&zeros),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
        });
        Ok(Self { ctx: ctx.clone(), buf, len })
    }

    /// Allocates and uploads a `GpuVecF64` from an `f64` host slice without narrowing.
    pub fn from_f64_slice(ctx: &GpuContext, data: &[f64]) -> Result<Self, GpuError> {
        Self::require_f64(ctx)?;
        let buf = if data.len() >= 4 {
            ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gpu-vec-f64"),
                contents: bytemuck::cast_slice(data),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            })
        } else {
            let mut pad = [0.0f64; 4];
            pad[..data.len()].copy_from_slice(data);
            ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("gpu-vec-f64"),
                contents: bytemuck::bytes_of(&pad),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            })
        };
        Ok(Self { ctx: ctx.clone(), buf, len: data.len() })
    }

    #[inline]
    pub fn context(&self) -> &GpuContext {
        &self.ctx
    }

    #[inline]
    pub fn buffer(&self) -> &wgpu::Buffer {
        &self.buf
    }

    /// Reads the GPU vector back to host memory as `Vec<f64>`.
    pub fn to_host_f64(&self) -> Result<Vec<f64>, GpuError> {
        self.ctx.read_buffer_f64(&self.buf, self.len)
    }

    /// Overwrites this GPU vector from an `f64` host slice of matching length.
    pub fn write_f64(&mut self, data: &[f64]) -> Result<(), GpuError> {
        if data.len() != self.len {
            return Err(GpuError::DimensionMismatch { expected: self.len, got: data.len() });
        }
        if !data.is_empty() {
            self.ctx.inner.queue.write_buffer(&self.buf, 0, bytemuck::cast_slice(data));
        }
        Ok(())
    }

    fn dispatch_binary(&mut self, other: &Self, alpha: f64, beta: f64, is_axpy: bool) {
        let n = self.len.min(other.len);
        if n == 0 {
            return;
        }
        let inner = &*self.ctx.inner;
        let Some(f64p) = inner.f64_pipes.as_ref() else {
            return;
        };
        let (gx, gy, wg_stride) = workgroup_grid_256(n);
        let params = VecParamsF64 { alpha, beta, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
        let ubuf = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vec-f64-bin-uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("vec-f64-bin-bg"),
            layout: &inner.binary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: ubuf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: other.buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: self.buf.as_entire_binding() },
            ],
        });
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("vec-f64-bin-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("vec-f64-bin-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(if is_axpy { &f64p.pipe_axpy } else { &f64p.pipe_aypx });
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
    }

    fn dispatch_scale(&mut self, alpha: f64) {
        if self.len == 0 {
            return;
        }
        let inner = &*self.ctx.inner;
        let Some(f64p) = inner.f64_pipes.as_ref() else {
            return;
        };
        let (gx, gy, wg_stride) = workgroup_grid_256(self.len);
        let params = VecParamsF64 { alpha, beta: 0.0, n: self.len as u32, wg_stride, _pad0: 0, _pad1: 0 };
        let ubuf = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("vec-f64-scale-uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("vec-f64-scale-bg"),
            layout: &inner.unary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: ubuf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: self.buf.as_entire_binding() },
            ],
        });
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("vec-f64-scale-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("vec-f64-scale-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_scale);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
    }
}

#[cfg(feature = "gpu")]
impl Vector for GpuVecF64 {
    type Scalar = f64;

    #[inline]
    fn len(&self) -> usize {
        self.len
    }

    fn dotc(&self, other: &Self) -> f64 {
        if self.len != other.len {
            return f64::NAN;
        }
        let res = if std::ptr::eq(&self.buf, &other.buf) {
            self.ctx.reduce_two_buffers_f64(&self.buf, None, self.len)
        } else {
            self.ctx.reduce_two_buffers_f64(&self.buf, Some(&other.buf), self.len)
        };
        res.unwrap_or(f64::NAN)
    }

    fn dotu(&self, other: &Self) -> f64 {
        self.dotc(other)
    }

    fn norm2_sq(&self) -> f64 {
        self.ctx.reduce_two_buffers_f64(&self.buf, None, self.len).unwrap_or(f64::NAN)
    }

    fn axpy(&mut self, a: f64, x: &Self) {
        self.dispatch_binary(x, a, 0.0, true);
    }

    fn scale(&mut self, a: f64) {
        self.dispatch_scale(a);
    }

    fn scale_real(&mut self, r: f64) {
        self.dispatch_scale(r);
    }

    fn aypx(&mut self, beta: f64, x: &Self) {
        self.dispatch_binary(x, 0.0, beta, false);
    }

    fn copy_from(&mut self, other: &Self) {
        let n = self.len.min(other.len);
        if n == 0 {
            return;
        }
        let mut encoder = self.ctx.inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("gpu-vec-f64-copy-enc"),
        });
        encoder.copy_buffer_to_buffer(&other.buf, 0, &self.buf, 0, (n * 8) as u64);
        self.ctx.inner.queue.submit(Some(encoder.finish()));
    }

    fn to_host(&self) -> Vec<f64> {
        self.to_host_f64().unwrap_or_else(|_| vec![f64::NAN; self.len])
    }

    fn from_host(&mut self, host: &[f64]) {
        if host.len() == self.len {
            let _ = self.write_f64(host);
        } else if let Ok(v) = Self::from_f64_slice(&self.ctx, host) {
            *self = v;
        }
    }

    fn zeros_like(&self) -> Self {
        Self::zeros(&self.ctx, self.len).unwrap_or_else(|_| self.clone())
    }
}

/// Warp-coalesced `SELL-32-128` sparse matrix in native `f64` on GPU (`wgpu::Features::SHADER_F64`).
#[cfg(feature = "gpu")]
pub struct GpuSellCSigmaF64 {
    ctx: GpuContext,
    nrows: usize,
    ncols: usize,
    nnz: usize,
    slice_ptr_buf: wgpu::Buffer,
    col_idx_buf: wgpu::Buffer,
    values_buf: wgpu::Buffer,
    row_lengths_buf: wgpu::Buffer,
    perm_buf: wgpu::Buffer,
    uniform_spmv: wgpu::Buffer,
    gx: u32,
    gy: u32,
}

#[cfg(feature = "gpu")]
impl GpuSellCSigmaF64 {
    /// Converts a host [`Csr<f64>`] into `SELL-32-128` and uploads full `f64` values to the GPU.
    pub fn from_csr(ctx: &GpuContext, csr: &Csr<f64>) -> Result<Self, GpuError> {
        if !ctx.supports_f64() {
            return Err(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"));
        }
        if csr.values().iter().any(|v| !v.is_finite()) {
            return Err(GpuError::InvalidArgument("CSR matrix contains non-finite entries"));
        }
        let sell = SellCSigma::from_csr_with(csr, DEFAULT_SELL_C, DEFAULT_SELL_SIGMA)
            .map_err(|_| GpuError::InvalidArgument("failed to convert CSR to SELL-32-128"))?;
        Self::from_sell(ctx, &sell)
    }

    /// Uploads a host `SellCSigma<f64>` (`C = 32`) to device `f64` storage buffers.
    pub fn from_sell(ctx: &GpuContext, sell: &SellCSigma<f64>) -> Result<Self, GpuError> {
        if !ctx.supports_f64() {
            return Err(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"));
        }
        if sell.c() != 32 {
            return Err(GpuError::InvalidArgument("GpuSellCSigmaF64 requires slice height C == 32"));
        }
        let mut sp_u32 = Vec::with_capacity(sell.slice_ptr().len());
        for &p in sell.slice_ptr() {
            sp_u32.push(u32::try_from(p).map_err(|_| GpuError::InvalidArgument("SELL slice_ptr overflows u32"))?);
        }
        let n_slices = (sell.slice_ptr().len().saturating_sub(1)) as u32;

        let mk_u32 = |label: &str, data: &[u32]| {
            if data.is_empty() {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::bytes_of(&[0u32; 4]),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            } else {
                ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some(label),
                    contents: bytemuck::cast_slice(data),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            }
        };
        let values_buf = if sell.values().is_empty() {
            ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("sell-f64-values"),
                contents: bytemuck::bytes_of(&[0.0f64; 4]),
                usage: wgpu::BufferUsages::STORAGE,
            })
        } else {
            ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("sell-f64-values"),
                contents: bytemuck::cast_slice(sell.values()),
                usage: wgpu::BufferUsages::STORAGE,
            })
        };

        let slice_ptr_buf = mk_u32("sell-f64-slice-ptr", &sp_u32);
        let col_idx_buf = mk_u32("sell-f64-col-idx", sell.col_idx());
        let row_lengths_buf = mk_u32("sell-f64-row-lengths", sell.row_lengths());
        let perm_buf = mk_u32("sell-f64-perm", sell.perm());

        let total_threads = (n_slices as usize) * 32;
        let (gx, gy, wg_stride) = workgroup_grid_256(total_threads);
        let sp_params = SellParams { nrows: sell.nrows() as u32, n_slices, wg_stride, mode: 0 };
        let uniform_spmv = ctx.inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("sell-f64-spmv-uniform"),
            contents: bytemuck::bytes_of(&sp_params),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        Ok(Self {
            ctx: ctx.clone(),
            nrows: sell.nrows(),
            ncols: sell.ncols(),
            nnz: sell.nnz(),
            slice_ptr_buf,
            col_idx_buf,
            values_buf,
            row_lengths_buf,
            perm_buf,
            uniform_spmv,
            gx,
            gy,
        })
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    pub fn nnz(&self) -> usize {
        self.nnz
    }

    /// Dispatches `y ← A x` in `f64` on the GPU.
    pub fn spmv(&self, x: &GpuVecF64, y: &mut GpuVecF64) -> Result<(), GpuError> {
        if x.len() != self.ncols {
            return Err(GpuError::DimensionMismatch { expected: self.ncols, got: x.len() });
        }
        if y.len() != self.nrows {
            return Err(GpuError::DimensionMismatch { expected: self.nrows, got: y.len() });
        }
        if self.nrows == 0 {
            return Ok(());
        }
        let inner = &*self.ctx.inner;
        let f64p = inner
            .f64_pipes
            .as_ref()
            .ok_or(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"))?;
        let bg = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("sell-f64-spmv-bg"),
            layout: &inner.sell_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: self.uniform_spmv.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: self.slice_ptr_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: self.col_idx_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: self.values_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: self.row_lengths_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: self.perm_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: x.buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: inner.dummy_f64.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 8, resource: y.buf.as_entire_binding() },
            ],
        });
        let mut encoder = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("sell-f64-spmv-enc"),
        });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("sell-f64-spmv-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_sell_spmv);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups(self.gx, self.gy, 1);
        }
        inner.queue.submit(Some(encoder.finish()));
        Ok(())
    }
}

#[cfg(feature = "gpu")]
impl DeviceLinearOperator<GpuVecF64> for GpuSellCSigmaF64 {
    #[inline]
    fn dim(&self) -> usize {
        self.nrows
    }
    fn apply(&self, x: &GpuVecF64, y: &mut GpuVecF64) {
        let _ = self.spmv(x, y);
    }
}

// ---------------------------------------------------------------------------------------
// GpuAmgF32: Device-Resident Chebyshev V-Cycle Preconditioner
// ---------------------------------------------------------------------------------------

#[cfg(feature = "gpu")]
struct GpuAmgFineLevel {
    nrows: usize,
    gx_vec: u32,
    gy_vec: u32,
    a_gx: u32,
    a_gy: u32,
    r_gx: u32,
    r_gy: u32,
    p_gx: u32,
    p_gy: u32,
    rhs_buf: wgpu::Buffer,
    x_buf: wgpu::Buffer,
    bg_res_from_ax: wgpu::BindGroup,
    bg_restrict: wgpu::BindGroup,
    bg_prolong_add: wgpu::BindGroup,
    bg_cheb_pre0: wgpu::BindGroup,
    bg_cheb_post0: wgpu::BindGroup,
    bg_cheb_steps: Vec<wgpu::BindGroup>,
}

#[cfg(feature = "gpu")]
struct HostTransferScratch {
    f32_in: Vec<f32>,
    coarse_rhs: Vec<f32>,
    coarse_x: Vec<f32>,
    fine_readback: wgpu::Buffer,
    coarse_readback: wgpu::Buffer,
}

/// GPU-resident Chebyshev smoothed-aggregation AMG V-cycle preconditioner (`f32` device cycle)
/// driving outer `f64` PCG / FGMRES or device-resident `GpuVecF32` Krylov solvers.
///
/// All fine-to-coarse pre-smoothing, residual, and restriction passes execute in a single
/// `wgpu::CommandEncoder` submission, the tiny coarsest level (`≤ max_coarse`, default `≤ 150`)
/// is solved via `faer` sparse LU, and all coarse-to-fine prolongation and post-smoothing passes
/// execute in a second `wgpu::CommandEncoder` submission.
#[cfg(feature = "gpu")]
pub struct GpuAmgF32 {
    ctx: GpuContext,
    host_amg: AmgF32,
    fine_levels: Vec<GpuAmgFineLevel>,
    coarse_rhs_buf: wgpu::Buffer,
    coarse_x_buf: wgpu::Buffer,
    coarse_nrows: usize,
    degree: usize,
    transfer: Mutex<HostTransferScratch>,
}

#[cfg(feature = "gpu")]
impl GpuAmgF32 {
    /// Builds the smoothed-aggregation AMG hierarchy on the host and uploads all level operators
    /// (`A_l`, `P_l`, `R_l` in `SELL-32-128`), inverse diagonals, scratch vectors, and pre-built
    /// bind groups to `ctx`.
    ///
    /// Because Gauss–Seidel is sequential, if `opts.smoother` is `GaussSeidel` the GPU cycle
    /// automatically uses a parallel degree-3 Chebyshev smoother.
    pub fn new(ctx: &GpuContext, a: &Csr<f64>, mut opts: AmgOptions) -> Result<Self, GpuError> {
        let degree = match opts.smoother {
            SmootherKind::Chebyshev { degree } => degree.max(1),
            SmootherKind::GaussSeidel => {
                opts.smoother = SmootherKind::Chebyshev { degree: 3 };
                3
            }
        };
        let host_amg = AmgF32::new(a, opts).map_err(|e| GpuError::Solver(e.to_string()))?;
        let levels = host_amg.levels_f32();
        let n_levels = levels.len();
        let inner = &*ctx.inner;

        let alloc_vec = |label: &str, len: usize| {
            let zeros = vec![0.0f32; len.max(4)];
            inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: bytemuck::cast_slice(&zeros),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            })
        };

        let coarse_nrows = levels[n_levels - 1].a.nrows;
        let coarse_rhs_buf = alloc_vec("amg-coarse-rhs", coarse_nrows);
        let coarse_x_buf = alloc_vec("amg-coarse-x", coarse_nrows);

        // Pre-allocate rhs and x buffers for every fine level so level `l` restriction/prolongation
        // can bind level `l + 1`'s `rhs_buf` and `x_buf` directly.
        let n_fine = n_levels.saturating_sub(1);
        let mut fine_rhs: Vec<wgpu::Buffer> = Vec::with_capacity(n_fine);
        let mut fine_x: Vec<wgpu::Buffer> = Vec::with_capacity(n_fine);
        for (lvl, lf) in levels.iter().take(n_fine).enumerate() {
            fine_rhs.push(alloc_vec(&format!("amg-l{lvl}-rhs"), lf.a.nrows));
            fine_x.push(alloc_vec(&format!("amg-l{lvl}-x"), lf.a.nrows));
        }

        let mut fine_levels = Vec::with_capacity(n_fine);
        // Move `fine_rhs` / `fine_x` out after creating cross-level bind groups.
        for lvl in 0..n_fine {
            let lf = &levels[lvl];
            let nrows = lf.a.nrows;
            let a_gpu = GpuSellCSigmaF32::from_csr(ctx, &lf.a_f64)?;
            let p_f64 = lf.p_f64.as_ref().ok_or(GpuError::InvalidArgument("missing P on fine AMG level"))?;
            let r_f64 = lf.r_f64.as_ref().ok_or(GpuError::InvalidArgument("missing R on fine AMG level"))?;
            let p_gpu = GpuSellCSigmaF32::from_csr(ctx, p_f64)?;
            let r_gpu = GpuSellCSigmaF32::from_csr(ctx, r_f64)?;

            let res_buf = alloc_vec(&format!("amg-l{lvl}-res"), nrows);
            let d_buf = alloc_vec(&format!("amg-l{lvl}-d"), nrows);
            let inv_diag_buf = if lf.inv_diag.is_empty() {
                alloc_vec(&format!("amg-l{lvl}-invd"), 4)
            } else {
                inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("amg-inv-diag"),
                    contents: bytemuck::cast_slice(&lf.inv_diag),
                    usage: wgpu::BufferUsages::STORAGE,
                })
            };

            let next_rhs = if lvl + 1 < n_fine { &fine_rhs[lvl + 1] } else { &coarse_rhs_buf };
            let next_x = if lvl + 1 < n_fine { &fine_x[lvl + 1] } else { &coarse_x_buf };

            let bg_res_from_ax = a_gpu.create_bind_group(1, &fine_x[lvl], &fine_rhs[lvl], &res_buf);
            let bg_restrict = r_gpu.create_bind_group(0, &res_buf, &inner.dummy_f32, next_rhs);
            let bg_prolong_add = p_gpu.create_bind_group(2, next_x, &inner.dummy_f32, &fine_x[lvl]);

            let (gx_vec, gy_vec, wg_stride) = workgroup_grid_256(nrows);
            let hi = 1.1f32 * lf.lmax.max(1e-12);
            let lo = hi / 30.0f32;
            let theta = 0.5f32 * (hi + lo);
            let delta = 0.5f32 * (hi - lo);
            let sigma = theta / delta;
            let inv_theta = 1.0f32 / theta;

            let mk_cheb_bg = |mode: u32, alpha: f32, beta: f32, src_res: &wgpu::Buffer| {
                let cp = ChebParams {
                    n: nrows as u32,
                    wg_stride,
                    mode,
                    _pad0: 0,
                    alpha,
                    beta,
                    _pad1: 0.0,
                    _pad2: 0.0,
                };
                let ubuf = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("amg-cheb-uniform"),
                    contents: bytemuck::bytes_of(&cp),
                    usage: wgpu::BufferUsages::UNIFORM,
                });
                inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("amg-cheb-bg"),
                    layout: &inner.cheb_bgl,
                    entries: &[
                        wgpu::BindGroupEntry { binding: 0, resource: ubuf.as_entire_binding() },
                        wgpu::BindGroupEntry { binding: 1, resource: inv_diag_buf.as_entire_binding() },
                        wgpu::BindGroupEntry { binding: 2, resource: src_res.as_entire_binding() },
                        wgpu::BindGroupEntry { binding: 3, resource: d_buf.as_entire_binding() },
                        wgpu::BindGroupEntry { binding: 4, resource: fine_x[lvl].as_entire_binding() },
                    ],
                })
            };

            let bg_cheb_pre0 = mk_cheb_bg(0, inv_theta, 0.0, &fine_rhs[lvl]);
            let bg_cheb_post0 = mk_cheb_bg(1, inv_theta, 0.0, &res_buf);

            let mut bg_cheb_steps = Vec::with_capacity(degree.saturating_sub(1));
            let mut rho_old = 1.0f32 / sigma;
            for _ in 1..degree {
                let rho_new = 1.0f32 / (2.0f32 * sigma - rho_old);
                let c1 = rho_new * rho_old;
                let c2 = 2.0f32 * rho_new / delta;
                bg_cheb_steps.push(mk_cheb_bg(2, c1, c2, &res_buf));
                rho_old = rho_new;
            }

            fine_levels.push((
                nrows,
                gx_vec,
                gy_vec,
                a_gpu.gx,
                a_gpu.gy,
                r_gpu.gx,
                r_gpu.gy,
                p_gpu.gx,
                p_gpu.gy,
                bg_res_from_ax,
                bg_restrict,
                bg_prolong_add,
                bg_cheb_pre0,
                bg_cheb_post0,
                bg_cheb_steps,
            ));
        }

        let assembled_fine: Vec<GpuAmgFineLevel> = fine_levels
            .into_iter()
            .zip(fine_rhs)
            .zip(fine_x)
            .map(
                |(
                    (
                        (
                            nrows,
                            gx_vec,
                            gy_vec,
                            a_gx,
                            a_gy,
                            r_gx,
                            r_gy,
                            p_gx,
                            p_gy,
                            bg_res_from_ax,
                            bg_restrict,
                            bg_prolong_add,
                            bg_cheb_pre0,
                            bg_cheb_post0,
                            bg_cheb_steps,
                        ),
                        rhs_buf,
                    ),
                    x_buf,
                )| GpuAmgFineLevel {
                    nrows,
                    gx_vec,
                    gy_vec,
                    a_gx,
                    a_gy,
                    r_gx,
                    r_gy,
                    p_gx,
                    p_gy,
                    rhs_buf,
                    x_buf,
                    bg_res_from_ax,
                    bg_restrict,
                    bg_prolong_add,
                    bg_cheb_pre0,
                    bg_cheb_post0,
                    bg_cheb_steps,
                },
            )
            .collect();

        let n0 = levels[0].a.nrows;
        let fine_readback = inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("amg-fine-readback"),
            size: (n0.max(4) * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let coarse_readback = inner.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("amg-coarse-readback"),
            size: (coarse_nrows.max(4) * 4) as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Ok(Self {
            ctx: ctx.clone(),
            host_amg,
            fine_levels: assembled_fine,
            coarse_rhs_buf,
            coarse_x_buf,
            coarse_nrows,
            degree,
            transfer: Mutex::new(HostTransferScratch {
                f32_in: vec![0.0; n0],
                coarse_rhs: vec![0.0; coarse_nrows],
                coarse_x: vec![0.0; coarse_nrows],
                fine_readback,
                coarse_readback,
            }),
        })
    }

    /// Convenience constructor building `GpuAmgF32` with a degree-`degree` Chebyshev smoother.
    pub fn new_chebyshev(ctx: &GpuContext, a: &Csr<f64>, degree: usize) -> Result<Self, GpuError> {
        Self::new(
            ctx,
            a,
            AmgOptions { smoother: SmootherKind::Chebyshev { degree: degree.max(1) }, ..Default::default() },
        )
    }

    #[inline]
    pub fn context(&self) -> &GpuContext {
        &self.ctx
    }

    #[inline]
    pub fn dim(&self) -> usize {
        self.host_amg.dim()
    }

    #[inline]
    pub fn stats(&self) -> AmgStats {
        self.host_amg.stats()
    }

    #[inline]
    pub fn chebyshev_degree(&self) -> usize {
        self.degree
    }

    fn encode_downward(&self, encoder: &mut wgpu::CommandEncoder) {
        let inner = &*self.ctx.inner;
        for fl in &self.fine_levels {
            if fl.nrows == 0 {
                continue;
            }
            // Pre-smooth step 0 (x starts at 0): d = (1/theta)*D^-1*rhs, x = d
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-pre0"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_cheb);
                pass.set_bind_group(0, &fl.bg_cheb_pre0, &[]);
                pass.dispatch_workgroups(fl.gx_vec, fl.gy_vec, 1);
            }
            for step_bg in &fl.bg_cheb_steps {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-pre-step"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_res_from_ax, &[]);
                pass.dispatch_workgroups(fl.a_gx, fl.a_gy, 1);

                pass.set_pipeline(&inner.pipe_cheb);
                pass.set_bind_group(0, step_bg, &[]);
                pass.dispatch_workgroups(fl.gx_vec, fl.gy_vec, 1);
            }
            // Restrict residual: res = rhs - A x; next_rhs = R * res
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-restrict"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_res_from_ax, &[]);
                pass.dispatch_workgroups(fl.a_gx, fl.a_gy, 1);

                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_restrict, &[]);
                pass.dispatch_workgroups(fl.r_gx, fl.r_gy, 1);
            }
        }
    }

    fn encode_upward(&self, encoder: &mut wgpu::CommandEncoder) {
        let inner = &*self.ctx.inner;
        for fl in self.fine_levels.iter().rev() {
            if fl.nrows == 0 {
                continue;
            }
            // Prolongate & add coarse correction: x += P * next_x; then post-smooth step 0
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-prolong-post0"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_prolong_add, &[]);
                pass.dispatch_workgroups(fl.p_gx, fl.p_gy, 1);

                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_res_from_ax, &[]);
                pass.dispatch_workgroups(fl.a_gx, fl.a_gy, 1);

                pass.set_pipeline(&inner.pipe_cheb);
                pass.set_bind_group(0, &fl.bg_cheb_post0, &[]);
                pass.dispatch_workgroups(fl.gx_vec, fl.gy_vec, 1);
            }
            for step_bg in &fl.bg_cheb_steps {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-post-step"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&inner.pipe_sell_spmv);
                pass.set_bind_group(0, &fl.bg_res_from_ax, &[]);
                pass.dispatch_workgroups(fl.a_gx, fl.a_gy, 1);

                pass.set_pipeline(&inner.pipe_cheb);
                pass.set_bind_group(0, step_bg, &[]);
                pass.dispatch_workgroups(fl.gx_vec, fl.gy_vec, 1);
            }
        }
    }

    fn map_staging_into(device: &wgpu::Device, staging: &wgpu::Buffer, dst: &mut [f32]) -> Result<(), GpuError> {
        if dst.is_empty() {
            return Ok(());
        }
        let byte_len = (dst.len() * 4) as u64;
        let slice = staging.slice(0..byte_len);
        let (tx, rx) = mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        let _ = device.poll(wgpu::PollType::wait_indefinitely());
        let map_res = rx.recv().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        map_res.map_err(|e| GpuError::BufferMap(e.to_string()))?;
        let mapped = slice.get_mapped_range().map_err(|e| GpuError::BufferMap(e.to_string()))?;
        dst.copy_from_slice(bytemuck::cast_slice(&mapped));
        drop(mapped);
        staging.unmap();
        Ok(())
    }

    /// Runs one GPU V-cycle directly between device-resident `GpuVecF32` vectors `r` and `z`.
    pub fn apply_gpu(&self, r: &GpuVecF32, z: &mut GpuVecF32) -> Result<(), GpuError> {
        let n = self.dim();
        if r.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: r.len() });
        }
        if z.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: z.len() });
        }
        if n == 0 {
            return Ok(());
        }

        let inner = &*self.ctx.inner;
        let mut tr = self
            .transfer
            .lock()
            .map_err(|_| GpuError::BufferMap("AMG transfer mutex poisoned".to_string()))?;

        if self.fine_levels.is_empty() {
            // Single-level (already <= max_coarse): solve directly via coarse factorisation.
            let mut enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("amg-1lvl-down"),
            });
            enc.copy_buffer_to_buffer(&r.buf, 0, &tr.coarse_readback, 0, (n * 4) as u64);
            inner.queue.submit(Some(enc.finish()));
            let HostTransferScratch { coarse_readback, coarse_rhs, coarse_x, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, coarse_readback, coarse_rhs)?;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            inner.queue.write_buffer(&z.buf, 0, bytemuck::cast_slice(coarse_x));
            return Ok(());
        }

        // 1. Single submission for fine copy + all downward levels + coarse readback copy.
        let mut down_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-down-enc"),
        });
        down_enc.copy_buffer_to_buffer(&r.buf, 0, &self.fine_levels[0].rhs_buf, 0, (n * 4) as u64);
        self.encode_downward(&mut down_enc);
        if self.coarse_nrows > 0 {
            down_enc.copy_buffer_to_buffer(
                &self.coarse_rhs_buf,
                0,
                &tr.coarse_readback,
                0,
                (self.coarse_nrows * 4) as u64,
            );
        }
        inner.queue.submit(Some(down_enc.finish()));

        // 2. Exact solve on tiny coarsest grid (<= max_coarse unknowns).
        if self.coarse_nrows > 0 {
            let HostTransferScratch { coarse_readback, coarse_rhs, coarse_x, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, coarse_readback, coarse_rhs)?;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            inner.queue.write_buffer(&self.coarse_x_buf, 0, bytemuck::cast_slice(coarse_x));
        }

        // 3. Single submission for all upward levels + copy into `z.buf`.
        let mut up_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-up-enc"),
        });
        self.encode_upward(&mut up_enc);
        up_enc.copy_buffer_to_buffer(&self.fine_levels[0].x_buf, 0, &z.buf, 0, (n * 4) as u64);
        inner.queue.submit(Some(up_enc.finish()));
        Ok(())
    }

    /// Runs one mixed-precision GPU V-cycle for an outer `f64` Krylov solver (`z ≈ A⁻¹ r`).
    pub fn apply_host_f64(&self, r: &[f64], z: &mut [f64]) -> Result<(), GpuError> {
        let n = self.dim();
        if r.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: r.len() });
        }
        if z.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: z.len() });
        }
        if n == 0 {
            return Ok(());
        }

        let (to_s, from_s) = f32_safe_scale(r);
        let inner = &*self.ctx.inner;
        let mut tr = self
            .transfer
            .lock()
            .map_err(|_| GpuError::BufferMap("AMG transfer mutex poisoned".to_string()))?;

        if self.fine_levels.is_empty() {
            for (dst, &src) in tr.coarse_rhs.iter_mut().zip(r) {
                *dst = (src * to_s) as f32;
            }
            let HostTransferScratch { coarse_rhs, coarse_x, .. } = &mut *tr;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            for (dst, &src) in z.iter_mut().zip(coarse_x.iter()) {
                *dst = (src as f64) * from_s;
            }
            return Ok(());
        }

        for (dst, &src) in tr.f32_in.iter_mut().zip(r) {
            *dst = (src * to_s) as f32;
        }
        inner.queue.write_buffer(&self.fine_levels[0].rhs_buf, 0, bytemuck::cast_slice(&tr.f32_in));

        let mut down_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-host-down-enc"),
        });
        self.encode_downward(&mut down_enc);
        if self.coarse_nrows > 0 {
            down_enc.copy_buffer_to_buffer(
                &self.coarse_rhs_buf,
                0,
                &tr.coarse_readback,
                0,
                (self.coarse_nrows * 4) as u64,
            );
        }
        inner.queue.submit(Some(down_enc.finish()));

        if self.coarse_nrows > 0 {
            let HostTransferScratch { coarse_readback, coarse_rhs, coarse_x, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, coarse_readback, coarse_rhs)?;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            inner.queue.write_buffer(&self.coarse_x_buf, 0, bytemuck::cast_slice(coarse_x));
        }

        let mut up_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-host-up-enc"),
        });
        self.encode_upward(&mut up_enc);
        up_enc.copy_buffer_to_buffer(&self.fine_levels[0].x_buf, 0, &tr.fine_readback, 0, (n * 4) as u64);
        inner.queue.submit(Some(up_enc.finish()));

        {
            let HostTransferScratch { fine_readback, f32_in, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, fine_readback, f32_in)?;
            for (dst, &src) in z.iter_mut().zip(f32_in.iter()) {
                *dst = (src as f64) * from_s;
            }
        }
        Ok(())
    }

    /// Runs one mixed-precision GPU V-cycle directly between `GpuVecF64` device vectors (`z ≈ A⁻¹ r`)
    /// without transferring fine vectors across PCIe (casting `f64 ↔ f32` on the device with
    /// dynamic norm scaling).
    pub fn apply_gpu_f64(&self, r: &GpuVecF64, z: &mut GpuVecF64) -> Result<(), GpuError> {
        let n = self.dim();
        if r.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: r.len() });
        }
        if z.len() != n {
            return Err(GpuError::DimensionMismatch { expected: n, got: z.len() });
        }
        if n == 0 {
            return Ok(());
        }
        let inner = &*self.ctx.inner;
        let f64p = inner
            .f64_pipes
            .as_ref()
            .ok_or(GpuError::InvalidArgument("GPU adapter does not support SHADER_F64"))?;

        let r_nrm = r.norm2();
        let (to_s, from_s) = if r_nrm > 0.0 && r_nrm.is_finite() && (r_nrm < 1e-10 || r_nrm > 1e10) {
            (1.0 / r_nrm, r_nrm)
        } else {
            (1.0, 1.0)
        };

        let mut tr = self
            .transfer
            .lock()
            .map_err(|_| GpuError::BufferMap("AMG transfer mutex poisoned".to_string()))?;

        if self.fine_levels.is_empty() {
            let mut enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("amg-1lvl-f64-down"),
            });
            let (gx, gy, wg_stride) = workgroup_grid_256(n);
            let p_dn = VecParamsF64 { alpha: to_s, beta: 0.0, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
            let u_dn = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("amg-cast-dn-u"),
                contents: bytemuck::bytes_of(&p_dn),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let bg_dn = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("amg-cast-dn-bg"),
                layout: &inner.binary_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_dn.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: r.buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: self.coarse_rhs_buf.as_entire_binding() },
                ],
            });
            {
                let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-cast-dn-pass"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&f64p.pipe_cast_64_to_32);
                pass.set_bind_group(0, &bg_dn, &[]);
                pass.dispatch_workgroups(gx, gy, 1);
            }
            enc.copy_buffer_to_buffer(&self.coarse_rhs_buf, 0, &tr.coarse_readback, 0, (n * 4) as u64);
            inner.queue.submit(Some(enc.finish()));
            let HostTransferScratch { coarse_readback, coarse_rhs, coarse_x, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, coarse_readback, coarse_rhs)?;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            inner.queue.write_buffer(&self.coarse_x_buf, 0, bytemuck::cast_slice(coarse_x));
            let mut up_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("amg-1lvl-f64-up"),
            });
            let p_up = VecParamsF64 { alpha: from_s, beta: 0.0, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
            let u_up = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("amg-cast-up-u"),
                contents: bytemuck::bytes_of(&p_up),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let bg_up = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("amg-cast-up-bg"),
                layout: &inner.binary_bgl,
                entries: &[
                    wgpu::BindGroupEntry { binding: 0, resource: u_up.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 1, resource: self.coarse_x_buf.as_entire_binding() },
                    wgpu::BindGroupEntry { binding: 2, resource: z.buf.as_entire_binding() },
                ],
            });
            {
                let mut pass = up_enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("amg-cast-up-pass"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&f64p.pipe_cast_32_to_64);
                pass.set_bind_group(0, &bg_up, &[]);
                pass.dispatch_workgroups(gx, gy, 1);
            }
            inner.queue.submit(Some(up_enc.finish()));
            return Ok(());
        }

        let (gx, gy, wg_stride) = workgroup_grid_256(n);
        let p_dn = VecParamsF64 { alpha: to_s, beta: 0.0, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
        let u_dn = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("amg-f64-cast-dn-u"),
            contents: bytemuck::bytes_of(&p_dn),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg_dn = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("amg-f64-cast-dn-bg"),
            layout: &inner.binary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: u_dn.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: r.buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: self.fine_levels[0].rhs_buf.as_entire_binding() },
            ],
        });

        let mut down_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-f64-down-enc"),
        });
        {
            let mut pass = down_enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("amg-f64-cast-dn-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_cast_64_to_32);
            pass.set_bind_group(0, &bg_dn, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        self.encode_downward(&mut down_enc);
        if self.coarse_nrows > 0 {
            down_enc.copy_buffer_to_buffer(
                &self.coarse_rhs_buf,
                0,
                &tr.coarse_readback,
                0,
                (self.coarse_nrows * 4) as u64,
            );
        }
        inner.queue.submit(Some(down_enc.finish()));

        if self.coarse_nrows > 0 {
            let HostTransferScratch { coarse_readback, coarse_rhs, coarse_x, .. } = &mut *tr;
            Self::map_staging_into(&inner.device, coarse_readback, coarse_rhs)?;
            self.host_amg.solve_coarse_f32(coarse_rhs, coarse_x);
            inner.queue.write_buffer(&self.coarse_x_buf, 0, bytemuck::cast_slice(coarse_x));
        }

        let p_up = VecParamsF64 { alpha: from_s, beta: 0.0, n: n as u32, wg_stride, _pad0: 0, _pad1: 0 };
        let u_up = inner.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("amg-f64-cast-up-u"),
            contents: bytemuck::bytes_of(&p_up),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let bg_up = inner.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("amg-f64-cast-up-bg"),
            layout: &inner.binary_bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: u_up.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: self.fine_levels[0].x_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: z.buf.as_entire_binding() },
            ],
        });

        let mut up_enc = inner.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("amg-f64-up-enc"),
        });
        self.encode_upward(&mut up_enc);
        {
            let mut pass = up_enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("amg-f64-cast-up-pass"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&f64p.pipe_cast_32_to_64);
            pass.set_bind_group(0, &bg_up, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        inner.queue.submit(Some(up_enc.finish()));
        Ok(())
    }
}

#[cfg(feature = "gpu")]
impl DevicePreconditioner<GpuVecF32> for GpuAmgF32 {
    fn apply(&self, r: &GpuVecF32, z: &mut GpuVecF32) {
        let _ = self.apply_gpu(r, z);
    }
}

#[cfg(feature = "gpu")]
impl DevicePreconditioner<GpuVecF64> for GpuAmgF32 {
    fn apply(&self, r: &GpuVecF64, z: &mut GpuVecF64) {
        let _ = self.apply_gpu_f64(r, z);
    }
}

#[cfg(feature = "gpu")]
impl Preconditioner<f64> for GpuAmgF32 {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        if self.apply_host_f64(r, z).is_err() {
            // Graceful host fallback if the GPU device was lost mid-solve.
            Preconditioner::apply(&self.host_amg, r, z);
        }
    }
}

#[cfg(feature = "gpu")]
impl DevicePreconditioner<CpuVec<f64>> for GpuAmgF32 {
    #[inline]
    fn apply(&self, r: &CpuVec<f64>, z: &mut CpuVec<f64>) {
        Preconditioner::apply(self, &r.0, &mut z.0);
    }
}
