# 67. Axis conditions and the cylindrical stretch of the axisymmetric solid kernel

Date: 2026-11-14

## Status

Draft (capability R6, workstream A, bookmark `capability-r6-sms-axis`). Part 1 (K1, the axis) is written here. **Part 2 (K2, the cylindrical stretch of the absorbing layer) is added by workstream B**; its "What was verified" lines go below part 1's. ADR 66 (the study) is workstream C/D's. Amends ADR 43 ("a mesh that touches the axis is refused as before" no longer holds). Plan: `SMS_AXISYM_PLAN.md` section 4.3.

## Context

Benchmark A1 (a disk on a post on a substrate) has bodies that touch the axis `r = 0`. ADR 43 refused such a mesh. The refusal was a check `r <= 0` at a quadrature point, which the interior quadrature rules (`quad_rule`: three points of degree 2, six of degree 4, none on an edge) can only fail on a degenerate cell, so assembly already ran on a mesh that touches the axis. What was missing is the regularity of `u(r, z) exp(-i m theta)` at `r = 0`: without essential conditions on the axis dofs the `1/r` terms of the energy act only as a penalty and a mode may carry a spurious displacement on the axis.

## Decision (K1)

* The axis dofs are the dofs (vertices and P2 edge midpoints) with `r <= 1e-9` of the extent of the mesh (`axisym::axis_tolerance`). `axisym::axis_constraints` returns the components to hold at zero and `assemble_cross_section` adds them to the fixed set before the reduction, so a boundary condition and an axis condition on one dof are the same dof held once.
* Derivation: finiteness at `r -> 0` of `e_tt = (u_r - i m u_t) / r`, `g_zt = du_t/dz - i m u_z / r` and `g_rt = du_t/dr - u_t / r - i m u_r / r` gives `u_r = i m u_t`, `u_z = 0` for `m != 0`, and then `(1 - m^2) u_t = 0` from `g_rt`:

  | `m` | held at zero on the axis |
  |---|---|
  | 0 | `u_r`, `u_theta` (`u_z` free) |
  | +-1 | `u_z`; `u_r` and `u_theta` stay free (the tie `u_r = i m u_theta` is not imposed, `ROUGH_SPOTS.md` R6a) |
  | \|m\| >= 2 | all three |

* An azimuthal order that is not an integer is refused when the mesh touches the axis (a single-valued displacement has an integer order). A mesh with a vertex at `r < -1e-9` of its extent is refused by cell; so is a cell whose three vertices lie on the axis (no area). Meshes off the axis are untouched: no existing number changes (`ring_model_axes` golden, `axisym_closed_form.rs`, `ring_frame.rs` pass unchanged).
* The log of the result says which components were held on how many axis dofs.
* **Shared function signature touched:** `CrossAssembly` gains the public field `axis_dofs: usize`. Nothing else in the public stage API changes.

## Consequences

* Solid disks, posts and substrates can be solved; the radial contour modes of a thin disk now have a closed-form check.
* For `m = +-1` the axis is regular only in the weak sense (the tie is a penalty); rigid translation and tilt, which satisfy it exactly, come out at zero frequency.

## What was verified

Tests `crates/nyx-sms/tests/axisym_axis.rs` (second-order elements; the printed errors are the measured ones, see the status section of `docs/RESTRUCTURE_STATUS.md`): T1 torsion of a free solid cylinder (`u_theta = J1(kr)`, zeros of `J2`), T2 radial plane-strain modes with `u_z = 0` on the ends (`u_r = J1(kr)`, `k R J0 / J1 = (1 - 2 nu) / (1 - nu)`), T3 the radial contour modes of a thin disk (`lambda J0 / J1 = 1 - nu`) and `m = 2` stability under refinement, plus `m = 1` rigid modes and the refusal of a mesh that crosses the axis.

## What was not verified

* The exact tie for `m = +-1`: converge-rate of an `m = +-1` mode that touches the axis was not measured beyond the rigid modes.
* Anisotropic materials on the axis: the kernel averages the stiffness about the axis (`average_about_axis`), as before.

## Code

`crates/nyx-sms/src/{axisym,cross_section}.rs`, `crates/nyx-sms/tests/axisym_axis.rs`.
