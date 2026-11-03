# 8. The ported solver library is merged but not wired into the physics

Date: 2026-10-31

## Status

Accepted as an interim state (commit `7a9a85e`, "Merge the solver and mesher ports into femcore"). The swap-in is open work.

## Context

`femsolve` (written in a separate repository, edition 2024, no access to nyx) provides one implementation of Krylov solvers, eigensolvers (Krylov-Schur, LOBPCG, Jacobi-Davidson, contour, Chebyshev filtering, NLEIGS, TOAR), preconditioners (ILU/ILUT/IC, Schwarz, AMG, AMS/ADS, field-split), sparse direct solvers with symbolic reuse and a thread policy, time steppers, optimisation, reduced-order sweeps and matrix-free operators, for real and complex scalars. nyx's physics runs on older modules (`eigen`, `lobpcg`, `amg`, `linalg`, `csparse`, `sparse`) that duplicate parts of it.

## Decision

* The library is moved in as `femcore::solvers` with its tests (`tests/solvers/*`), studies and CLI tools as examples (`fs_eps`, `fs_gen`, `fs_ksp`, `fs_contour`, `fs_snes`). The PETSc/SLEPc comparison harness stays in the original repository because nyx ships no Python.
* **No kernel calls it yet.** The integration order is in `docs/ports/solvers/NYX_INTEGRATION.md` (adapters, reproduce, compare, thread policy, swap the eigensolver behind `ShiftInvertEigenSolver`, keep the old one selectable). Its statements about nyx are unverified because its author could not see nyx.
* The stage boundary of ADR 9 is where a swap happens: "the backend is chosen in the solve stage, so a `femcore::solvers` adapter replaces one function per kernel".
* The library keeps determinism (fixed chunks, results bit-identical for every thread count; `NYX_THREADS` selects the pool through `crate::par`) and the `portable-kernels` feature, which trades the GEMM-based restart multiply for CPU-independent last bits.
* Licence: BSD-2-Clause with PETSc/SLEPc notices; compatibility with nyx's licence is unchecked.

## Consequences

* Duplicate AMG, LOBPCG, Arnoldi and four sparse-matrix types live side by side (R02).
* Any claim about the ported solvers' speed or robustness is about the library on its own test problems and the Linux VM it was written on, not about nyx runs.
* faer 0.24 is used by both stacks; the ported library touches it in two files only (`faer_glue.rs`, `direct.rs`).

## Code

`crates/femcore/src/solvers/**`, `examples/fs_*.rs`, `tests/solvers/**`, `docs/ports/solvers/*`; old stack `amg.rs`, `eigen.rs`, `lobpcg.rs`, `linalg.rs`, `csparse.rs`, `sparse.rs`.

## Contract

`contracts/kernel_stages`: each kernel's `solve_*` stage is the swap point and stays composition-equivalent to the one-call solver (`tests/stages.rs`).
