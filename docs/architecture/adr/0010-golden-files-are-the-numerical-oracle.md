# 10. Golden files are the numerical oracle, and are not regenerated to hide a difference

Date: 2026-11-02

## Status

Accepted (rule in `CLAUDE.md`, a047506). Open: the reference platform and tolerance (R22).

## Context

Porting the kernels to wafer axes (ADR 1) needed an oracle that did not depend on the code being changed. The kernels were run through `solve::solve` before each phase and their outputs frozen.

## Decision

* `crates/femcore/tests/golden_frames/*.json` (read by `tests/golden_frames.rs`, 24 tests) and `crates/geom/golden/*.json` with `*.expected.json` (read by `tests/goldens.rs`; closed-form areas) are the numerical oracle.
* `NYX_WRITE_GOLDEN=1` (femcore) and `GEOM_WRITE_GOLDEN=1` (geom) rewrite them; do so only for a kernel known to be right, and say so in the commit. Never regenerate to make a failing test pass before the cause is understood.
* The golden comparison accepts an eigenvector with all its signed fields negated, because the overall sign of a mode is round-off when two entries tie in modulus.
* Tolerances are per case (1e-9 for `dfs`, 1e-8 to 1e-7 for 2D EM and solid, 1e-6 for Brillouin and rings).

## Consequences

* Seven golden tests fail on main on this machine (Windows, Rust 1.99.1) by ~1e-6 relative, identical with `NYX_THREADS=1`, and already at `35cceba`, the commit that wrote them. The goldens carry the platform of their freezing (CPU, FMA, GEMM blocking). Until the reference platform or a tolerance is decided, new contract tests in this work do not depend on these goldens.
* The ring goldens were written with the kernel of `1f2cb58` and the Brillouin goldens with a kernel that still propagated along -y on the mirrored crystal; they check +y indirectly (displacement along y is mirrored in the files).
* `solid_pml::low_reflecting_boundary_gives_the_analytic_damped_bar_modes` fails before and after the port and has no golden.

## Code

`crates/femcore/tests/golden_frames.rs`, `golden_frames/`, `crates/geom/tests/goldens.rs`, `crates/geom/golden/`.

## Contract

`contracts/frames` and `contracts/kernel_stages` complement the goldens with checks that do not need frozen floating-point values.
