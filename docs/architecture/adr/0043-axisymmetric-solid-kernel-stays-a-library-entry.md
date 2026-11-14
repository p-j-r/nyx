# 43. Bodies of revolution stay a library entry of `nyx-sms`

Date: 2026-11-04

## Status

Amended by ADR 68 and ADR 69: the library entry stays, the study `axisymmetric_eigenfrequency` is added (ADR 68), a mesh may touch the axis and the layers are cylindrical (ADR 69). Accepted (restructure phase 4, bookmark `restructure-4-sms`). Follows the gate answer G1 (keep the axisymmetric kernels, with a verification note) and ADR 24 (the out-of-plane wave-number model is removed).

## Context

`sms_axisym` (stiffness and mass of a body of revolution on its (r, z) half-plane at azimuthal order `m`) and `sms_kz` (the cross-section assembly, damping, absorbing layers and eigensolve it is built on) were reachable only through `solve_solid_axisym`, called by two test files; no study of the envelope and no window control selects them. paul decided that `sms_kz` went (the axis is y now) and that the axisymmetric kernel stays, with its implementation checked.

## Decision

* The two files move into `nyx-sms` as `axisym` and `cross_section` (the latter was `sms_kz`; its type `KzModes` is `CrossModes` and its field `kz` is `m_az`, since there is no wave number any more). They read the typed `Problem` (an `Eigenfrequency` study; the azimuthal order is an argument) and give a `nyx_physics::SolveResult` with the scalar fields `|u|`, `ux`, `uz`, `uy`, `|uy|` and complex components per cell, plus the `CrossModes` (complex frequencies and mode shapes) the checks read.
* They stay **a library entry** (`nyx_sms::solve_solid_axisym`), not a study of the envelope. A study `axisymmetric` with the parameter `m` would change the descriptors and the sizing (the mesh is a half-plane with `x = r`) and is a decision for the owner of the window; nothing needs it before phase 7.
* **Verification** (G1, phase 3, numbers from phase 3; the tests below pass on the moved code with their thresholds): the breathing mode of a 1 m ring (0.02 m square section) against `sqrt(E/rho)/(2 pi R)` to 3.4e-5, in-plane flexure m = 2 and m = 3 against Euler-Bernoulli ring theory to 1.9e-4 and 5.0e-4 (theory itself is off by that much); the tests are `nyx-sms/tests/axisym_closed_form.rs`. The frame checks of `ring_frame.rs` (crystal ring independent of the azimuth of the c axis, the sign of theta along wafer +y, torsion as `uy`) and the frozen frequencies `ring_model_axes` pass unchanged.

## Consequences

* The ring tests keep running on the moved code; the electromagnetic ring (`ems_axisym`) is the ems thread's.
* A mesh that touches the axis `r = 0` is refused as before. *(Superseded by ADR 69: a mesh may touch the axis, with regularity conditions on the axis dofs.)*

## Code

`crates/nyx-sms/src/{axisym,cross_section}.rs`, `crates/nyx-sms/tests/{axisym_closed_form,ring_frame}.rs`.
