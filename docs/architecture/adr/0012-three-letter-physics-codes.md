# 12. Physics modules carry three-letter codes

Date: 2026-11-03

## Status

Accepted (Paul's choice; commit `d2e0f74`).

## Context

Module names mixed physics and implementation history (`scalar`, `solid`, `emw`, `brillouin`, `ring_brillouin`) and the same words meant different things in different crates (`solvers::scalar` is the solver library's scalar trait).

## Decision

| Code | Physics | Modules |
|---|---|---|
| `dfs` | scalar diffusion (electrostatics, heat) | `dfs` |
| `sms` | solid mechanics | `sms`, `sms_pml`, `sms_floquet`, `sms_lrb`, `sms_kz`, `sms_axisym` |
| `ems` | EM waves | `ems`, `ems_pml`, `ems3d`, `ems_axisym` |
| `sbs` | stimulated Brillouin scattering | `sbs`, `sbs_ring` |

New physics modules are named the same way: three letters, then a suffix for the variant. The public `Physics` enum is unchanged (`Electrostatics`, `Heat`, `Solid`, `Emw`, `Brillouin`); the codes name modules, not the type.

## Consequences

* The rename left comments, timer labels and test file names with the old prefixes (`emw3d`, `emw_pml`, `solid_pml`) (R13).
* `periodic.rs` and `analysis.rs`/`polarization.rs` are not under a code. `periodic` is used by `sms_floquet` and `ems3d`, so it belongs to no single physics.
* The codes are unrelated to the three-letter names of the kernel stages.

## Code

`crates/femcore/src/{dfs,sms*,ems*,sbs*}.rs`, `lib.rs` module list, `docs/KERNEL_PORT_STATUS.md` table.

## Contract

None needed beyond compilation; the module list is part of `contracts/kernel_stages`.
