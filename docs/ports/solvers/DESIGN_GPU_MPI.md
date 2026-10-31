# Design notes: GPU kernels (G1) and distributed memory (D1)

**Status: design only. Nothing in this file was implemented or measured.** The session that wrote it
had neither a GPU nor a cluster (the agreement with Paul: "list order; GPU (G1) and MPI (D1) only as documented
design notes"). Each section says what the existing `femsolve` code gives, what has to be written, and how it would be tested.
Estimates are reasoning from the measurements in `LEARNINGS.md`, not benchmarks.

---------------------------------------------------------------------------------------------------

## G1. GPU kernels for AMG-PCG and LOBPCG (port list G1; CubeCL, mixed precision)

### Which parts are worth moving

From `LEARNINGS.md` 8–9: eigensolves and AMG-PCG are **bandwidth-bound** (SpMV, Gauss–Seidel/Chebyshev sweeps, dots, axpys, block Gram matrices).
CPU streaming is 8–10 GB/s on the VM used here; a GPU offers 10–50× that, so the gain is real *only* if the data stays on the device for a whole
solve. Port in this order:

1. **SpMV** (`Csr::spmv`, `Bsr::spmv`): CSR-vector (one warp per row) for rows ≥ 16 nonzeros, SELL-C-σ (C = 32, σ = 128) for short rows (P3 of the list, item "SELL-C-sigma").
   BSR 3×3 maps to one thread per block row with the three `x` values in registers.
2. **Fused vector kernels**: `cg_update` (`vecops`), `dotc`/`norm2` (two-pass or atomic-free tree reduction), the `Basis` kernels `dots_block`, `sub_block`, `sub_block_then_dots`
   (these *are* batched GEMV/GEMM: `V^H w`, `w -= V h`); use a library GEMM when available (cuBLAS/rocBLAS through CubeCL interop) for `k ≥ 32`.
3. **AMG solve phase on the device**, setup on the host. `Amg::new` is sequential setup (aggregation, `P`, Galerkin product); it can stay on the CPU and the
   hierarchy (`A_l`, `P_l`, `R_l`, `D_l⁻¹`) is uploaded once. **Gauss–Seidel is sequential**, so the device cycle needs a parallel smoother: Chebyshev (degree 2–3 on `D⁻¹A`,
   eigenvalue bound from a few power iterations at setup) or multi-colour GS. Both are in `amg.rs`' roadmap ("a coloured or Chebyshev smoother is the route to threads"); this is the
   first piece of CPU code to write, because it also speeds up the CPU parallel path and can be validated on the CPU against the existing V-cycle (same convergence factor ±20 %).
   The coarsest level (≤150 unknowns) is a dense LU on the host or a few Jacobi sweeps on the device.
4. **Mixed precision**: V-cycle in `f32` as preconditioner for an `f64` PCG/FGMRES outer iteration (the outer residual in `f64` guarantees the answer; the
   preconditioner does not need more than 1e-3 accuracy). The module `mixed.rs` (S3) shows the pattern for direct solvers; for AMG it is a type parameter
   on the hierarchy (`Amg<f32>` built from `Csr<f32>`; needs `Scalar for f32`, currently only `f64`/`C64`).
5. **LOBPCG** (`lobpcg.rs`): per iteration it is 2–3 block SpMVs, block orthogonalisation (`B`-Gram matrices of size ≤ 3m × 3m) and a small dense eigenproblem.
   Keep the tall-skinny blocks on the device; the `3m × 3m` Rayleigh–Ritz problem goes to the host (`faer`) — a 3m of 30–300 is a few ms and 90 KB of traffic.
   Complex Hermitian needs `Complex<f32/f64>` kernels (CubeCL supports vectors of floats; complex multiplication has to be written out).

### Memory budget (reasoning)
`n = 2·10⁷` real unknowns in 3-D, 27 nonzeros per row: `A` is 6.5 GB in CSR/`f64` (4.4 GB in `f32` values); AMG adds ≈ 1.4–1.6× (operator complexity measured 1.3–1.6 on the CPU tests);
LOBPCG with block 12 needs `6·12·n·8` bytes ≈ 11.5 GB. A 24 GB card holds one such problem in `f32` preconditioner / `f64` iteration form; 16 GB needs the block reduced.
**Nobody should plan beyond that without measuring the actual card.**

### API shape that keeps the CPU path untouched
`LinearOperator<T>` and `Preconditioner<T>` are already the only interfaces the Krylov solvers use (`apply`, `apply_dotc`, `apply_block`). A device matrix implements them with
**device-resident vectors**, which the present `&[T]` signatures cannot express. So the real prerequisite is a vector abstraction:

```rust
trait Vector { type Scalar; fn dot(&self, &Self) -> Self::Scalar; fn axpy(&mut self, a, &Self); ... }
```
and generic `Cg<V: Vector>`. This is a refactor of `ksp/`, `lobpcg.rs` and `vecops.rs` (about 1 500 lines touched); do it on the CPU first (a `Vec<T>` implementation must reproduce today's bit-exact results,
which `tests/determinism.rs` checks) before any GPU code exists.

### What was *not* verified and is needed before starting
* Which GPU model, driver and whether CubeCL's backend (WGPU / CUDA / HIP) is acceptable on Paul's machines; whether `f64` throughput is usable (consumer cards run `f64` at 1/32 rate: then mixed precision is not an optimisation, it is a requirement).
* Whether reductions must be bit-reproducible on the GPU (`par.rs` promises that on the CPU). A fixed-order two-level tree reduction is possible; atomics are not.
* Break-even size: below roughly 10⁵ unknowns kernel-launch latency (≈ 10 µs × tens of launches per V-cycle) dominates; the CPU path should remain the default there.

### Test plan once a device exists
1. Kernel-level: SpMV/dots/axpy against the CPU result to 1e-14 (`f64`) / 1e-6 (`f32`), across `Csr`, `Bsr`, the 3-D Laplacian and the Nédélec matrices of `maxwell.rs`.
2. Solver-level: AMG-PCG iteration counts equal to the CPU's (Chebyshev smoother on both) on `lap3d`, then the same on the complex-symmetric PML matrices (P2).
3. Performance: GB/s of SpMV against the device's STREAM number (should be > 70 %), V-cycle time, end-to-end speed-up for 10⁶, 10⁷ unknowns **including host↔device transfers of the setup**.

---------------------------------------------------------------------------------------------------

## D1. Distributed memory (MPI) matrices and AMG (port list D1)

### What already helps
* `Schwarz` (`schwarz.rs`) is the single-process picture of a distributed solve: each subdomain is "what one rank owns", with overlap by graph layers, exact local solves and an optional coarse correction.
  Its iteration counts equal PETSc's `PCASM` (tier A), so a distributed version has an oracle to be checked against.
* `ShiftedFactorization` per node (`contour.rs`, `sweep.rs`) is embarrassingly parallel over contour nodes or frequencies; **an MPI rank per node/frequency needs no distributed matrix at all** and covers E3/W1 on a cluster of nodes with local memory big enough for one factorisation.
  This is the cheapest "cluster" step and could be done first (scatter the matrices once, gather the small dense results).
* Determinism: `par.rs` folds fixed-size chunks in index order. The same rule extends to ranks: partial sums per chunk, reduced in a fixed tree over rank order → results independent of the number of *threads* but dependent on the number of *ranks* (document it; a rank-count-independent result needs chunks that do not straddle ranks).

### Components, in order
1. **Distributed CSR**: row partition (contiguous after a graph partitioner such as METIS/ParMETIS or the existing `graph_partition` as a placeholder), local `Csr` for owned rows with columns split into owned/ghost; halo
   exchange with nonblocking send/recv, overlapped with the interior part of the SpMV. Index type: `u32` local, `u64` global (the port list notes the 32-bit limit).
2. **Distributed vector ops**: `dotc`, `norm2`, `axpy`, block `dots_block` with one `MPI_Allreduce` per block (not per vector); this is why `Basis::dots_block` and DGKS refinement (`eps/orth.rs`) are already block-oriented (one reduction per Gram–Schmidt pass).
3. **Krylov solvers**: unchanged algorithms on top of the vector abstraction described in G1 (the same refactor serves both). CG/GMRES need 1–2 allreduces per iteration; pipelined variants are possible later.
4. **Preconditioners**: block Jacobi/ILU(k) per rank (`ilu.rs` on the owned block), additive/restricted Schwarz across ranks (halo of one overlap layer per level; `Schwarz` becomes the single-rank case), field-split per field.
5. **Distributed AMG**: smoothed aggregation with *rank-local aggregation* (aggregates never cross rank boundaries; this keeps setup communication-free and costs convergence at the interface — measure), `P` block-diagonal across ranks, Galerkin product `RAP` with a
   parallel sparse matrix–matrix product (needs the halo rows of `P`), coarse levels agglomerated onto fewer ranks below a size threshold, direct (`faer`) solve of the coarsest level on one rank. Smoother: Chebyshev or hybrid GS (GS inside the rank, Jacobi across).
6. **Eigensolvers**: Krylov–Schur with a distributed basis (the dense Schur machinery is replicated on every rank — it is `m × m`), shift-and-invert needs a distributed direct solver (out of scope: MUMPS/PARDISO are explicitly "No" in the list) or an iterative inner solve with the preconditioners above (Jacobi–Davidson `jd.rs`, LOBPCG `lobpcg.rs` are the natural fits).

### Crates and tooling
`mpi` (rsmpi) needs a system MPI at build time; Debian's `libopenmpi-dev` is already required by the comparison harness, so a 2–4 rank test on one machine is possible here (`mpirun --oversubscribe -n 4`). On Windows (Paul's machine) MS-MPI is the usual choice;
rsmpi's Windows support is limited, so confirm before committing to it. An alternative that avoids system MPI altogether is a thin trait `Comm` with a thread-based implementation for tests and an rsmpi implementation for clusters.

### Test plan
1. 1-rank distributed code reproduces the serial result bit for bit; 2/4/8 ranks reproduce the serial *eigenvalues/solutions* to 1e-10 and the *iteration counts* of the serial solver within ±5 % for block-Jacobi-type preconditioners (exactly for Krylov methods without preconditioner).
2. Weak and strong scaling of SpMV/CG/AMG-PCG on the 3-D Laplacian and the Nédélec matrices; communication volume per iteration measured and compared with the surface-to-volume estimate.
3. Distributed contour solver (rank per node): eigenvalue sets equal to the serial `contour_eigs` result to 1e-12.

### What was not decided
Whether the target is a few large-memory nodes (threads + the contour/frequency-sweep parallelism above may already be enough) or a real cluster. **This decides whether D1 should be done at all**
(port list question 9, still open): the answer determines if the distributed-vector refactor is worth its cost.
