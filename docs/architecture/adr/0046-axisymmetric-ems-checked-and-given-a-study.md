# 46. The axisymmetric EM kernel is checked against closed forms and gets a study (gate G1)

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-ems`). Answers the G1 action for the `ems_axisym` half (paul: keep the axisymmetric kernels, but check their implementation before they move). The solid half (`sms_axisym`) is the sms thread's.

## Context

`ems_axisym` solves the modes of a body of revolution at azimuthal order `m` (`exp(-j m theta)`, theta along wafer +y) on the (r, z) half plane: `E_r, E_z` in the degree-2 Nedelec space of the waveguide solver and `E_theta` in P2, a Hermitian pencil `curl-curl x = (omega/c)^2 eps x`. Before this phase it had no study: nothing in `ProblemSpec` could ask for it, and it was reached only by tests. Paul was not sure of the implementation.

## What was verified

All numbers are in `crates/nyx-ems/tests/axisym.rs` (run with `--nocapture` for the printed errors), second-order elements:

* **Coaxial cavity, PEC, `a < r < b`, `0 < z < h`, m = 0**: the TEM modes have `f = p c / (2h)` for any `a`, `b`. Relative errors 1.4e-6 (p = 1) and 8e-6 (p = 2). The same cavity at m = 1 has nothing at the TEM frequency.
* **PEC cylinder cavity, TM modes at m = 1, 2, 3**: `f = c x_{m,1} / (2 pi R)` with `x_{m,1}` the first zero of the Bessel function `J_m`. This is the check the coaxial case cannot make: it exercises the `m / r` terms of the curl and the `E_theta` unknown. A thin rod of radius `a` on the axis keeps the mesh off the axis: relative errors 1.1e-5 (m = 2) and 2.7e-5 (m = 3); at m = 1 the error is 1.3e-3 with a 2 mm rod and 3.5e-4 with a 1 mm rod, so the rest at m = 1 is the rod (about `(a/R)^2` times a factor), not the discretisation.
* **Direction**: for a silicon wire in silica on a 20 mm ring, `m > 0` and `m < 0` give the same frequency and the power along theta (wafer y) is positive for `m > 0` and negative for `m < 0`, equal in size to 1e-3.

## What was not verified

* A dielectric whispering-gallery mode against an analytic or independent reference (only the symmetry in `m` and the sign of the power are checked); the Bessel-zero test is the closest independent check of the azimuthal terms.
* Absorbing layers in the axisymmetric kernel (`AxiMaterial::with_pml`, stretch of `r` with `r~`/`r` for the azimuthal component) are exercised by the code path of the `pml` BC but not compared with a closed form or with a larger domain. Open boundaries are not modelled otherwise: the walls are perfect conductors, adequate for well-confined modes in a large enough domain.
* Order 1 elements, large `m` (the mesh must resolve `m / r`), meshes that touch the axis (not supported; the mesh must stay off `r = 0`).

## Decision

1. The kernel moves unchanged (same formulas) into `nyx_ems::axisym`.
2. It gets a study, `axisymmetric_modes` (parameters `f_search`, signed `azimuthal_order`, `n_modes`), a fields stage (`axisym_fields_cells`: `|E|` and `Re Ex` (r), `Re Ey` (theta), `Re Ez` (z) with the phase rotated so the largest component is real, plus complex per-cell components), and a solve (`solve_axisymmetric`): permittivity and permeability per domain from `eps_r` and `mu_r`, absorbing layers from the `pml` condition (`x` is the r direction, `z` the z direction). `AxiMaterial::from_materials` builds the material data. Results: one solution per mode with its frequency (`value`), scalars `f N (Hz)` and `Q N` when the layers make it complex.
3. The G1 note in the status file is satisfied by the list above; the gaps are open items for whoever uses the kernel for a real resonator.

## Consequences

* The kernel is reachable from the envelope and from the engine; it has no goldens (there is nothing frozen from the old code to compare with) and its tests are closed forms.
* The window cannot ask for it (the old problem has no axisymmetric EM study); it appears with the engine commands in phase 7.

## Code

`crates/nyx-ems/src/axisym.rs`, `crates/nyx-ems/tests/axisym.rs`.

## Contract

`physics_ems.v1.json` 1.1 (the `axisymmetric_modes` study), `physics_ems_result.v1.json`.
