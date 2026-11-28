# ADR 146: Hardware acceleration — `Vector` abstraction, opt-in `wgpu` compute kernels, and 64-bit CSR storage (`FR-43`)

* **Date:** 2026-11-27
* **Status:** Accepted
* **Applies to:** `crates/femcore/src/linalg/{device,gpu,amg,ams,mod}.rs`, `crates/femcore/Cargo.toml`, `crates/nyx-ems/src/{problem,lib,dgtd,dgtd_gpu}.rs`, `crates/nyx-ems/src/fdtd/{mod,gpu}.rs`, `crates/nyx-ems/Cargo.toml`, `crates/nyx-engine/Cargo.toml`

## Context

`nyx` targets desktop workstations and compute nodes ranging from multi-core CPUs without discrete GPUs to workstations equipped with Vulkan / DirectX 12 / Metal GPUs and large-memory nodes where sparse finite-element matrices exceed $2^{31} - 1$ nonzeros (`nnz > 4.29e9`). Prior to `FR-43`:

1. Iterative Krylov and eigensolver routines (`cg`, `fgmres`, `lobpcg`) operated directly on `&[T]` / `&mut [T]` host slices, preventing transparent execution against device-resident GPU buffers or custom vector backends without duplicating the Krylov state machines.
2. Explicit time-domain electromagnetic solvers (`Study::Fdtd` in ADR 145 and `Study::Dgtd` in ADR 144) executed exclusively on the CPU via SIMD / `femcore::par` chunked parallelism, even when a discrete GPU capable of high-bandwidth structured stencil and batched element matrix-free updates was present.
3. Classical CSR row pointer arrays used `u32` offsets (`col_ind: Vec<u32>`, `row_ptr: Vec<u32>`), bounding single-matrix storage at $4.29 \times 10^9$ nonzeros, and AMG smoothers relied on sequential lexicographical Gauss–Seidel sweeps that serialize on wide SIMD/SIMT architectures.

## Decision

We adopt a strictly opt-in, deterministic, pure-Rust hardware-acceleration architecture (`FR-43`) split into two complementary layers:

### 1. Linear algebra vector abstraction, SIMD/GPU sparse formats, polynomial smoothers, and `Csr64` (`femcore::linalg`)

1. **`Vector` trait and generic Krylov / block eigensolvers (`femcore::linalg::device`):**
   - `Vector` abstracts device-resident linear combinations (`axpy`, `aypx`, `scale`, `scale_real`, `copy_from`, `add_combination`), deterministic inner products (`dotc`, `dotu`), and norms (`norm2_sq`, `norm2`) over `Scalar` (`f64` and `C64`, plus `f32` via `GpuVecF32`).
   - `CpuVec<T>` delegates directly to `vecops` and `kernels` for 0-ULP bit-for-bit parity with the existing slice-based solvers across all `NYX_THREADS` counts, while `tree_reduce_sum`, `tree_reduce_dotc`, `tree_reduce_sum_f32`, and `tree_reduce_dot_f32` implement atomic-free fixed-order two-level binary reduction trees (`BLOCK = 256`).
   - `DeviceLinearOperator<V>` and `DevicePreconditioner<V>` decouple matrix-vector products and preconditioning from host memory, driving `cg_vec`, `fgmres_vec`, and `lobpcg_vec` while keeping existing slice-based entry points unchanged.
2. **Sliced ELLPACK (`SellCSigma`, $C = 32, \sigma = 128$) and 64-bit CSR (`Csr64`):**
   - `SellCSigma<T>` converts `Csr<T>` into column-major $C$-row slices (`DEFAULT_SELL_C = 32`) with bounded row-length sorting within $\sigma$-windows (`DEFAULT_SELL_SIGMA = 128`) and explicit `perm` / `inv_perm` tables, preserving 0-ULP SpMV parity on CPU and enabling warp-coalesced GPU SpMV.
   - `Csr64<T>` uses `indptr: Vec<u64>` and `indices: Vec<u64>` (`from_triplets64`, `spmv_checked`) for matrices with `nnz > u32::MAX` or column indices above `u32::MAX`, with checked conversion (`TryFrom<&Csr64<T>> for Csr<T>`).
3. **Parallel Chebyshev polynomial smoother and mixed-precision `AmgF32` / `AmsF32` (`femcore::linalg::{amg,ams}`):**
   - `SmootherKind::Chebyshev { degree }` computes a deterministic power-iteration spectral bound $\lambda_{\max}(D^{-1} A)$ and applies degree-$k$ Chebyshev polynomial smoothing over $[\lambda_{\max}/30, 1.1\lambda_{\max}]$ using only SpMV and diagonal scaling (zero sequential data dependencies; default remains `SmootherKind::GaussSeidel` so existing CPU goldens stay 0-ULP unchanged).
   - `AmgF32` and `AmsF32` execute the entire AMG / Hiptmair–Xu auxiliary-space Maxwell V-cycle in single precision (`f32`) with dynamic RHS norm scaling (`f32_safe_scale`), implementing `Preconditioner<f64>` to drive outer `f64` PCG / FGMRES solves to full `f64` residual tolerance at half the memory bandwidth.
4. **Opt-in shared `wgpu` backend and opportunistic `SHADER_F64` (`femcore::linalg::gpu`, `--features gpu`):**
   - Always-compiled query helpers (`GpuDeviceInfo` including `supports_f64: bool`, `GpuError`, `gpu_available`, `probe_gpu_device`) allow callers to inspect GPU and double-precision capabilities without `#[cfg]` gates.
   - Behind `--features gpu`, `GpuContext::shared()` selects the highest-priority compute adapter (preferring `DiscreteGpu`), opportunistically requests `wgpu::Features::SHADER_F64` when advertised by the adapter (`ctx.supports_f64()`), and exposes its `wgpu::Device` and `wgpu::Queue` (`wgpu_device()`, `wgpu_queue()`) so `femcore::linalg::gpu` and `nyx-ems` (`EmsGpuContext::acquire_default()`) share a single logical GPU context and command queue.
   - Alongside `GpuVecF32` and `GpuSellCSigmaF32`, `femcore::linalg::gpu` provides native double-precision `LINALG_F64_WGSL`, `GpuVecF64`, and `GpuSellCSigmaF64` when `SHADER_F64` is supported, plus a mixed-precision `DevicePreconditioner<GpuVecF64> for GpuAmgF32` bridge driving full-precision `f64` Krylov outer iterations with `f32` GPU AMG V-cycles.

### 2. Opt-in `wgpu` WGSL time-domain electromagnetic kernels (`nyx-ems`, `--features gpu`)

1. **`backend`, `precision`, and opt-in `nonlinear` study parameters:**
   - Both `Study::Fdtd` (`FdtdStudyConfig`) and `Study::Dgtd` (`DgtdConfig`) expose `ParamDesc::choice("backend", "auto", &["auto", "cpu", "gpu"])`, `ParamDesc::choice("precision", "auto", &["auto", "f32", "f64"])` (`TdPrecision`), and `ParamDesc::boolean("nonlinear", false)` (opt-in toggle, default `false`).
   - When `nonlinear == false`, assembly in both `assemble_fdtd` and `assemble_dgtd` auto-skips library `[nonlinear]` material blocks so linear simulations on materials that also catalog nonlinear constants (e.g., `silicon_nitride`, `lithium_niobate`, `silicon`) stay on fast linear/GPU paths without unintended Kerr/SHG/TPA overhead or weak-field cancellation.
   - Precision selection (`select_td_precision`):
     - `"f32"` executes on GPU in `f32` (using normalized $E_{\text{ref}}$ and $N_{c,\text{ref}} = 10^{20}\ \text{m}^{-3}$ scaling when nonlinear dynamics are active).
     - `"f64"` executes on GPU using `SHADER_F64` when supported, or triggers `"gpu_fallback_to_cpu"` (CPU `f64`) when the adapter lacks `SHADER_F64`.
     - `"auto"` selects `f64` when `nonlinear && ctx.supports_f64()`, and `f32` otherwise (keeping large linear runs at full `f32` bandwidth).
   - In default builds without `--features gpu`, `"cpu"` and `"auto"` execute the existing CPU kernels with 0-ULP bit-exact identity and zero additional dependencies (`wgpu`, `pollster`, and `bytemuck` are strictly optional behind `gpu`). Requesting `backend = "gpu"` without `--features gpu` returns a clear diagnostic error instructing the user to compile with `--features gpu`.
   - Under `--features gpu`:
     - `backend = "cpu"` always runs the CPU kernel.
     - `backend = "gpu"` dispatches to `wgpu` on any available compute adapter, falling back gracefully to CPU with a `SolverDiagnostics` warning (`DiagnosticWarning { code: "gpu_fallback_to_cpu", .. }`, ADR 133) if the problem requires multi-tier LTS (`n_lts_tiers > 1`), dispersive ADE poles, explicit `precision = "f64"` without `SHADER_F64`, or if no adapter is available.
     - `backend = "auto"` dispatches to `wgpu` when a discrete GPU (`wgpu::DeviceType::DiscreteGpu`) is available and the problem size meets the threshold (`gpu_min_cells = 250_000` for FDTD; `gpu_min_dofs = 50_000` for DGTD), and otherwise runs on the CPU.
2. **`Study::Fdtd` WGSL compute engine (`crates/nyx-ems/src/fdtd/gpu.rs`):**
   - Templated over `GpuScalar` (`f32` or native `SHADER_F64` `f64`), uploading the 64-byte-aligned SoA Yee grid arrays (`Ex, Ey, Ez, Hx, Hy, Hz`), material update tables `(ca, cb)` and `(da, db)`, non-uniform dual/primary inverse spacing tables, 12-slab CPML auxiliary $\psi$ buffers, Dey–Mittra conformal PEC masks, Farjadpour off-diagonal `PolarizedAveraging` stencils, and Group 3 nonlinear domain/edge buffers (`GpuNlDomain`, `GpuNlEdge`, rank-4 $\chi^{(3)}$, `nl_state`).
   - Executes the complete time-stepping loop on the GPU (`update_h`, `inject_h`, `offdiag_save/apply/commit`, `update_e`, `inject_e`, `sync_periodic_x/y/z`, `clamp_pec_e`, synchronous 3-iteration Newton `nl_save_old` / `nl_save_lin` / `nl_newton_iter` / `nl_commit_e` / `nl_finish_step` solving $\chi^{(2)}$ + QPM, isotropic and rank-4 $\chi^{(3)}$ Kerr, TPA/FCA conductivity losses, FCD, and `drude_carrier_step` with $E_{\text{ref}}$/$N_{c,\text{ref}}$ normalization for `f32`, `accum_planar_dft`, `accum_vol_dft`, `sample_time_monitors`, and two-level workgroup tree-reduction `energy_reduce_pass1`/`pass2`), reading back only the accumulated monitor spectra, time traces, and final field snapshots at the end of the simulation.
3. **`Study::Dgtd` WGSL compute engine (`crates/nyx-ems/src/dgtd_gpu.rs`):**
   - Implements a race-free, atomic-free nodal DG evaluation templated over `GpuScalar` (`f32` or native `SHADER_F64` `f64`):
     - Pass 1 (`eval_face_qps`, 1 thread per face quadrature point) evaluates interior, PEC, PMC, Silver–Müller, modal waveguide port, TFSF plane-wave, and penetrable 2D sheet (`FR-35`) Rankine–Hugoniot upwind/central Riemann fluxes.
     - Pass 2 (`eval_dof_rhs`, 1 thread per global nodal DOF) computes strong-form SBP volume curls ($D^x_K, D^y_K, D^z_K$), lifts surface quadrature fluxes, injects point dipoles, evaluates instantaneous nonlinear $\chi^{(2)}$/QPM, isotropic/rank-4 $\chi^{(3)}$, FCD, and TPA/FCA conductivity losses with an analytic $3\times 3$ $\boldsymbol{\varepsilon}_{\text{inst}}^{-1}$ solve ($E_{\text{ref}}$/$N_{c,\text{ref}}$-scaled in `f32`, $0.45\,\varepsilon_{\min}$ row-sum cap), and applies linear $(\boldsymbol{\varepsilon}^{-1}, \boldsymbol{\mu}^{-1}, \boldsymbol{\sigma})$ and uniaxial PML damping.
     - Pass 3 (`lserk4_stage` / `rk4_save_init` + `rk4_stage`) advances uniform-step `Lserk4`, `Rk4`, and 1-tier `Lts` states in-place on the GPU, followed by `step_carriers` once per macro step when TPA free carriers are active.
     - Pass 4 (`accum_dof_dft`, `sample_probes`, `dgtd_energy_pass1`/`pass2`) accumulates nodal DFT spectra, probe time traces, and fixed-order workgroup tree-reduced total energy on the GPU using normalized cell masses ($M_K / \bar{V}$ paired with dimensionless $\boldsymbol{\varepsilon}_r$ and $(\mu_0/\varepsilon_0)\boldsymbol{\mu}_r$, rescaled by $\varepsilon_0 \bar{V}$ in `f64` upon readback) so nanoscale meshes avoid `f32` subnormal underflow, and reconstructs modal `port_dft` (`DrivenPortData`) and surface `flux_dft` upon readback.
   - Multi-rate Local Time-Stepping (`n_lts_tiers > 1`), dispersive ADE poles, and explicit `precision = "f64"` on adapters without `SHADER_F64` fall back cleanly to CPU with a `"gpu_fallback_to_cpu"` warning.

### 3. Future GPU & HPC Roadmap

`wgpu`/WGSL serves as the zero-SDK portable Tier-1 compute layer across Vulkan, DirectX 12, and Metal today, while `femcore::linalg::device` (`Vector`, `DeviceLinearOperator`, `DevicePreconditioner`) remains strictly backend-agnostic:

* **Compute backend evolution:** `wgpu`/WGSL may be replaced by `CubeCL` or `rust-gpu` in the future once portable double-precision (`f64`) support and stable Rust toolchain integration mature.
* **HPC and vendor-native GPU backends:** Dedicated multi-node / cluster HPC support and direct native bindings for NVIDIA (CUDA) and AMD (ROCm/HIP) GPUs are planned for the future; if a user or deployment requires these capabilities on a priority basis, they should contact the maintainers with target hardware and workload details.

## Consequences

* **Zero default bloat and 0-ULP CPU preservation:** Default builds (`cargo build --release`) compile zero GPU dependencies and preserve every existing CPU golden and contract bit-for-bit.
* **Deterministic GPU reductions:** All GPU inner products and energy reductions use fixed-order workgroup shared-memory binary reduction trees (`@workgroup_size(64)` / `@workgroup_size(256)`) with zero floating-point atomics, ensuring bit-identical results across repeated GPU runs.
* **Strict crate layering:** `nyx-ems` depends only on `femcore` and `nyx-physics` (plus optional `wgpu`, `pollster`, `bytemuck` behind `gpu`), preserving the workspace dependency DAG verified by `crates/nyx-engine/tests/layering.rs`.
