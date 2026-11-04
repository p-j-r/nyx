# 11. Determinism: fixed chunks, ordered reductions, one thread pool

Date: 2026-10-31

## Status

Accepted. The rule is stated in `CLAUDE.md` ("the same input gives the same numbers") and in `solvers/par.rs`.

## Context

Parallel floating-point sums depend on the order of addition. Results that change with the thread count make goldens, caches and comparisons between runs worthless.

## Decision

* Parallel kernels cut work into **fixed-size chunks independent of the thread count**, compute partials per chunk and fold them **in chunk order**. A result is bit-identical for every thread count and equal to the serial run.
* `femcore::par` owns the one shared rayon pool; `NYX_THREADS=N` sets its size (`1` makes everything serial). `solvers::par` runs on that pool; small problems stay serial because waking a worker costs 50 to 150 microseconds.
* `NYX_PROFILE=1` echoes the solver log.
* `portable-kernels` removes the one known source of cross-machine drift (the faer GEMM restart multiply picks FMA/blocking by CPU); it never mattered between thread counts.

## Consequences

* `meshers::par` uses `std::thread::scope` and `FEMMESH_THREADS` instead of the pool (R03). The chunk rule is the same, so results agree, but the thread count is controlled by a different variable.
* Determinism is across threads, not across machines: see ADR 10.
* Tests: `tests/determinism.rs`, `tests/solvers/determinism`, and bit-identical reruns in the mesher tests.

## Code

`crates/femcore/src/par.rs`, `solvers/par.rs`.

## Contract

Covered by the existing determinism tests; `contracts/kernel_stages` does not repeat them.

## Amendment 2026-11-04 (mesh phase)

The mesher's own `par` module is gone: its fixed-chunk reductions (`CHUNK`, `fold_chunks_mut`, `sum_range`) moved into `femcore::par`, and the metric adaptation now takes its thread count from `femcore::par::threads()` (`NYX_THREADS`, at most 32); `FEMMESH_THREADS` no longer exists. Results are bit-identical for any thread count by the chunk rule, so nothing but the variable name changes (R03 is resolved for the mesher; the solver library keeps `solvers/par.rs`).
