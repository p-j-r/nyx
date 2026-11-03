# 21. Polarization of an optical mode is computed from its complex fields, in the solvers' local frame

Date: 2026-10-31

## Status

Accepted (commit `4107cfa`, "Analysis tab, ring and 3D polarization"; `e60ad11` for the first version). `docs/photonic_use_cases.md` is the description.

Superseded in part by ADR 24: the Brillouin study no longer produces an analysis; the optical mode set moved to `femcore::optical`.

## Context

Users of waveguide mode solvers want to know whether a mode is quasi-TE, quasi-TM, hybrid or circular. The plotted fields are real and lose the phase between components, so they cannot give the Stokes parameters of hybrid or leaky modes. Degenerate modes (the square 300 nm wire in oxide has an exactly degenerate TE/TM pair) come out of the eigensolver as arbitrary mixtures, so their individual polarizations mean nothing.

## Decision

* `analysis` and `polarization` work from the **complex fields the solver found at quadrature points** (with area weights), not from plotted fields. Every optical mode solver returns `SolveResult::analysis` next to its fields: the 2D mode solver (real, complex with PML or a tilted crystal), the 3D eigenmode solver, the Brillouin study of a straight guide (pump and Stokes) and of a ring (whispering-gallery modes).
* Per mode: integrated Stokes vector, a point on the Poincare sphere, degree of polarization, ellipticity and orientation, the TE and longitudinal energy shares, and a class (quasi-TE, quasi-TM, linear +-45, elliptical, circular, hybrid).
* **Degenerate modes** (effective indices within 5e-4, or frequencies within 1e-6 relative) form a group; the analysis reports the group Stokes vector (independent of the mixture) and the resolved pair of definite x/y polarization.
* The **sense of rotation** is defined by what can be seen (`S3 > 0`: the E vector turns counter-clockwise as time advances), not by a time convention. Left-handed axis triples (rings, wafer-frame 3D) flip the sign convention.
* The axes of the analysis are the local frame of the mode solvers: x lateral (quasi-TE), **y the surface normal (quasi-TM), z the propagation axis**. In the wafer frame these are x, z and y.

## Consequences

* A second naming of axes lives in the analysis, its documentation and the Analysis tab (x, y normal, z propagation), beside the public wafer frame of ADR 1 (R28).
* `SolveResult` carries an `Option<Analysis>`, so the result type depends on the analysis module (R27).
* `nyx-store::runs::metrics_of` reads `te_fraction` from the analysis rows.
* The class thresholds are heuristics; no closed-form benchmark other than the slab and the square wire exists (tests in `tests/photonics.rs`, `polarization` unit tests).

## Code

`crates/femcore/src/analysis.rs`, `polarization.rs`, `ems.rs`, `ems_pml.rs`, `ems3d.rs`, `sbs.rs`, `sbs_ring.rs`; `crates/app/src/tabs.rs` (Analysis tab); `docs/photonic_use_cases.md`. History: `e60ad11`, `4107cfa`.

## Contract

None of its own. `contracts/solve_result` notes that `analysis` is part of the result.
