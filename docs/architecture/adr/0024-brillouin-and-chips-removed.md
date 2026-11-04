# 24. Stimulated Brillouin scattering and the ring chips are removed

Date: 2026-11-03

## Status

Accepted by paul on 2026-11-03 (restructure plan, goal 8). Supersedes in part ADR 12 (the code `sbs`), ADR 21 (the Brillouin use of the analysis) and ADR 22 (the four ring chips). The code is kept in the tag `sbs-last` of the `nyx2` repository.

Amendment 2026-11-04: on paul's instruction the out-of-plane wave-number solid model (`ProblemSpec::out_of_plane`, `sms_kz::solve_solid_kz*`, its 4 goldens, `tests/solid_kz.rs`, the app checkbox) is removed as well. The cross-section solid machinery in `sms_kz.rs` stays because the axisymmetric model (`sms_axisym`, kept by G1) is built on it; it is renamed and moved in phase 4. The Floquet-vs-cross-section agreement test went with `solid_kz.rs`; `solid_eigen_floquet` remains as a golden.

## Context

nyx started as a Brillouin application. The restructure narrows it to a general finite-element application with a few physics per project. Brillouin scattering needed the optical modes, the elastic modes with a phase-matched out-of-plane wave number, photoelastic and moving-boundary overlaps and gain; it spread over `sbs.rs`, `sbs_ring.rs`, `chips.rs`, a `Physics::Brillouin` variant, a `BrillouinSpec` inside `ProblemSpec`, an application tab, five golden files and about 5000 lines of tests and studies. Five of the seven failing golden tests were Brillouin goldens.

## Decision

* `Physics::Brillouin`, `BrillouinSpec`, `PhotoSpec`/`PhotoKind`, `sbs.rs`, `sbs_ring.rs`, `chips.rs`, `tests/brillouin.rs`, the five `brillouin_*` golden files, the chip and rib examples and the Brillouin tab and example of the application are deleted.
* What other code needed from `sbs.rs` stays: the optical mode set with its field sampler, the vacuum constants and the 24-point quadrature rule move to `femcore::optical` (used by `analysis` and `polarization`, ADR 21).
* The kz solid kernel (`sms_kz`) and the axisymmetric kernels stay: they are ordinary studies (gate G1, default keep).
* Brillouin may return later as an optional physics crate through the `Physics` trait (ADR 25); the tag `sbs-last` holds the last working state.

## Consequences

* `ProblemSpec` shrinks to four physics; `problem_spec.v1.json` was rewritten on purpose (`VERSIONS.md`).
* The golden failures drop from seven to two (`golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`, both optical, R22).
* `geom` keeps its golden feature lists of the chip shapes (`crates/geom/golden`): they test the geometry kernel, not Brillouin.
* The `wg800` study in the notes of the earlier work is only reproducible from `sbs-last`.

## Code

Removed: `crates/femcore/src/{sbs,sbs_ring,chips}.rs`, `tests/brillouin.rs`, `examples/{chip_report,chip_digits,rib_shapes,rib_sweep,wg800}.rs`. Added: `crates/femcore/src/optical.rs`.

## Contract

`problem_spec` (rewritten, version 1.1).
