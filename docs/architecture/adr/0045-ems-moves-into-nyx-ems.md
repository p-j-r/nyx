# 45. The electromagnetic kernels move into `nyx-ems`; what they share with the analysis and the solid kernels stays in `femcore`

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-ems`). Carries out ADR 41 for `ems` (the second physics after `dfs`). Closes ROUGH_SPOTS R10 (the EM study is named, not guessed from the dimension of the mesh). Builds on ADR 25 (the `Physics` trait), ADR 9 (three stages) and ADR 12 (codes).

## Context

`femcore` had four EM modules: `ems` (real full-vector modes of a 2D cross-section), `ems_pml` (the same in complex arithmetic: absorbing layers, loss, a crystal tilted against the guide), `ems_axisym` (bodies of revolution) and `ems3d` (3D cavities and Floquet cells). They read the old `ProblemSpec` and dispatched on the dimension of the mesh. Three things outside them used their internals: `optical` (the field sampler of the polarization analysis builds on the edge-element basis and on the complex material tensors, and its `solve_optical` called the mode solvers), the solid kernels (`sms_pml`, `sms_kz` used the PML stretch of `ems_pml`), and the use cases (`usecases/{soi,bragg}` called `solve::solve`). `femcore` cannot depend on `nyx-ems` (the crate depends on `femcore`), and `analysis`, `polarization` and `optical` stay in `femcore` until phase 5 because the old `SolveResult` carries an `Analysis`.

## Decision

1. **Crate layout.** `nyx-ems` holds `problem` (typed `Problem`, `Study`, `decode`), `modes` (real modes), `modes_complex`, `axisym`, `cavity` (the `ems3d` kernel), `optical` (modes as the field-evaluable set, `solve_optical[_wafer]`) and, in `lib.rs`, the `Ems` physics, `run` and `run_analysed`. Kernels were copied and edited only for types (ADR 41, item 2): same formulas, same order of operations, which is why the goldens keep their numbers.
2. **The study names the kernel.** `modes` (2D), `axisymmetric_modes` (2D, the (r, z) half plane) and `cavity` (3D). `Physics::validate` and the kernels refuse a study on a mesh of the other dimension by name ("the study 'cavity' needs a 3D mesh"). The dispatch on the mesh dimension survives only in the bridge (ADR 47), because the old problem has no kind for it.
3. **What stays in `femcore`, and why.** `edge` (Whitney and degree-2 Nedelec basis on a triangle, per-cell numbering and signs), `pml` (`stretch`, `cell_stretches`) and, in `optical`, the vacuum constants including `C0`, `EpsC`/`MuC` and `cell_materials[_signed]`. They are what the field sampler and the solid kernels share with the EM kernels, so they sit below both. `periodic` (Floquet pairs), `section`, `material`, `fields`, `fem`, `analysis`, `polarization` also stay. `section::axes_of_3d_periodic` is the same function as `axes_of_3d_problem` taking the Floquet pairs alone. `Analysis::summary` is what `attach` adds to scalars and log, for results that are not the old `SolveResult`.
4. **Polarization stays in the result for now.** The 2D and cavity studies compute the polarization with `femcore::analysis` (as before the move) and put the scalars `polarization N` and the table into the new result; `run_analysed` returns the `Analysis` next to it for the window. That keeps the goldens, which record those scalars, identical. Phase 5 moves the computation into `nyx-analysis`, removes the scalars from the result, and rewrites the goldens deliberately (ADR 10).
5. **Descriptors (contract `physics_ems` 1.1).** The fields are scalar arrays under the flat names of the window and of the goldens (`|E|`, `|Exz|`, `|Ey|`, `Ex`, `Ey`, `Ez`; `Re Ex`, ... for complex modes and cavities), because the old results were not a vector field with components and the goldens compare these names. The axisymmetric study takes `f_search` (Hz), a signed `azimuthal_order` and `n_modes`; `n_target` is gone (the kernel searches by frequency). The shape of a result is pinned by `physics_ems_result` 1.0.
6. **Goldens and tests.** Five goldens (`emw_pec_straight`, `emw_pec_tilted`, `emw_pml_x`, `emw_soi_wire`, `emw_soi_hybrid`) moved unchanged to `crates/nyx-ems/tests/golden/`; the first three keep their tolerances (1e-7). The two silicon-wire files never reproduced on this machine (R22); gate G2 (option B) sets their tolerances to 5e-5 relative on n_eff and 1e-3 of the field maximum, with the class of the polarization compared as text. The ems tests of `emw2d`, `emw3d`, `pml`, `aniso_optics`, `bragg_stack`, `photonics`, `physics_meshing`, the SOI tests of `meshing` and the ems test of `anisotropic`, `frame_roundtrip`, `stages`, `adapt3d`, `wafer_names` and `ring_frame` moved with the kernels, with the same numbers and tolerances. `femcore::solve` answers an error for `Physics::Emw`.

## Consequences

* Positive: `ems` is a crate with its own problem type; the window keeps working through the bridge; the stage API is contract-tested in the crate; every moved golden passes with unchanged numbers (the two SOI files at G2 tolerances).
* Negative: `nyx-ems` depends on `femcore::analysis` until phase 5; four small shared modules sit in `femcore` only because of the crate direction; the envelope carries scalar permittivity only (crystal tensors reach the kernel through the typed problem until the material files, phase 6).
* The imports `crate::ems_pml::cell_stretches` in `sms_pml.rs` and `sms_kz.rs` are now `crate::pml::cell_stretches`; the sms thread's copies of those files need the same one-line change at the merge.
* Open: `Physics::sizing` still returns the placeholder (the wavelength per region is read by `physics_size::spec_sizing` from the old spec); the adapters go when the window is rewired (phase 7).

## Code

`crates/nyx-ems/src/*`, `crates/femcore/src/{edge,pml,optical}.rs`, `crates/nyx-engine/src/legacy/ems.rs`.

## Contract

`physics_ems.v1.json` 1.1, `physics_ems_result.v1.json` 1.0, the stage-API test in `crates/nyx-ems/tests/contracts.rs`.
