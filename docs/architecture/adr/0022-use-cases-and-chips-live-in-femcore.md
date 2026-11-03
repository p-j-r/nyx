# 22. Ready-made problems (SOI, Bragg, the four ring chips) live in femcore, built from geom models with named regions

Date: 2026-10-31

## Status

Accepted by practice (`e60ad11`, `4107cfa`; the chips predate them, `6de7458`). Re-expressed on geom by `08dd007` and `ffa1bf6`.

## Context

The physics modules need realistic, reproducible inputs to be tested and demonstrated: a 300 nm SOI wire at 1550 nm, symmetric and asymmetric Bragg gratings, and the silicon-on-sapphire rib waveguides of four ring chips (300TE, 300TM, 500TE, 500TM) as Brillouin problems with crystal orientations, absorbing layers and photoelastic data.

## Decision

* `femcore::usecases::{soi, bragg}` and `femcore::chips` return a model (a geom `FeatureList` with parameters), the `ProblemSpec` and the mesh plan; they call `meshing::mesh_model` / `mesh_evaluated` with `MeshPlan::simple`.
* Regions have **names** (`sapphire`, `slab`, `rib`, `air`, `sapphire.pml_xmin`, ...); the integers come from `LabelMap` (ADR 2), replacing the hard-coded domain numbers (chips 1-4, PML 11-15, 21-22, 41-45).
* Crystal orientations are written in the wafer frame (ADR 1): silicon axes `[p, q, z]`, r-plane sapphire `[a, c x a, c]` with the projection of the c axis at angle psi to the waveguide.
* The examples of the app (numbers 16 to 19) and the `examples/` studies (`chip_report`, `soi_modes`, `bragg_gratings`, `wg800`) call these functions; `docs/chip_report.md` and `docs/photonic_use_cases.md` record the numbers they produced on the machine that ran them.

## Consequences

* Application content lives in the solver library (R15): `femcore` cannot be built without the chips; `ai/nyx-session` calls `usecases::soi` (`SoiWire`, `slab_neff`) for its only template, so the MCP interface also depends on the application content.
* Their tests (`tests/usecases.rs`, `photonics.rs`, `bragg_stack.rs`) are the only end-to-end checks of whole studies besides the goldens.
* Reference values in docs are the authors' measurements on one machine (R22).

## Code

`crates/femcore/src/usecases/{mod,soi,bragg}.rs`, `chips.rs`, `crates/app/src/examples.rs`, `crates/femcore/examples/{chip_report,chip_digits,soi_modes,bragg_gratings,rib_sweep,rib_shapes,wg800}.rs`, `docs/chip_report.md`, `docs/photonic_use_cases.md`.

## Contract

None of its own; a changed region name breaks `usecases` tests.
