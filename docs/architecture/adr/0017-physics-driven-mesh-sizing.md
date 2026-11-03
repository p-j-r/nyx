# 17. Element sizes come from the physics of the problem

Date: 2026-11-02

## Status

Accepted as opt-in (`e5fad41`, `12042b4`; app option `e950f05`, checkbox `aa513ea`). Defaults are starting values calibrated on smooth modes only.

## Context

A mesh that resolves a wave needs a fixed number of points per wavelength in each material; a mesh for a high-index core must be finer than one for the cladding, and elastic modes depend on the slowest anisotropic plane wave. Users were choosing element sizes by hand.

## Decision

`physics_size` derives an element size per domain from the study: wavelength of EM modes and eigenfrequencies (n times shorter in index-n material), slowest anisotropic plane wave for solid eigenfrequencies, Brillouin (optical and acoustic), with **12 points per wavelength for order 1 and 6 for order 2**, PML relaxation and eigenfrequency headroom. `meshing::mesh_for_problem` intersects this metric with the geometry metric (2D); `physics_size::adapt_tets` does the same graded in 3D. Estimators (`estimate::Estimator::{combine, e_field}`) treat E components as one vector field with one metric.

## Consequences

* Only order 1 and 2 are covered; points per wavelength for higher orders are unmeasured (`MISSING_FEATURES.md` section 1).
* No test runs a Brillouin or ring study on a sized mesh; sharp corners, thin layers and high contrast are uncalibrated.
* The app shows the option as one checkbox on the model editor (`GeometrySpec::physics_sizing`), with no report of the sizes it used and no per-domain override.
* `physics_size` reads `ProblemSpec`, so the mesher depends on the problem types, and a change in `Study`/`BrillouinSpec` field names reaches the mesher.

## Code

`crates/femcore/src/physics_size.rs`, `meshing.rs`, `estimate.rs`, `crates/app/src/state.rs`. History: `e5fad41`, `12042b4`, `e950f05`, `aa513ea`.

## Contract

`contracts/problem_spec` (the field names `physics_size` reads: `Study::{Modes, Eigenfrequency}`, `BrillouinSpec::wavelength`, `Materials`).
