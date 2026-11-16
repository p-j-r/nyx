# 69. Axis conditions and the cylindrical stretch of the axisymmetric solid kernel

Date: 2026-11-14

## Status

Proposed (capability R6, the merge of the three workstreams: K1 axis conditions on `capability-r6-sms-axis`, K2 cylindrical stretch on `capability-r6-sms-pml`, joined on `capability-r6-sms-axisym`). Amends ADR 43 ("a mesh that touches the axis is refused as before" no longer holds, and the layers of the axisymmetric solid are now the cylindrical stretch). ADR 68 is the study that reaches this kernel. Numbered 69 because the A3 branch (`restructure-10-a3-snowflake`) already holds ADR 66 and 67. Plan: `SMS_AXISYM_PLAN.md` sections 4.3 and 8.

## Context

**Axis (K1).** Benchmark A1 (a disk on a post on a substrate) has bodies that touch the axis `r = 0`. ADR 43 refused such a mesh. The refusal was a check `r <= 0` at a quadrature point, which the interior quadrature rules (`quad_rule`: three points of degree 2, six of degree 4, none on an edge) can only fail on a degenerate cell, so assembly already ran on a mesh that touches the axis. What was missing is the regularity of `u(r, z) exp(-i m theta)` at `r = 0`: without essential conditions on the axis dofs the `1/r` terms of the energy act only as a penalty and a mode may carry a spurious displacement on the axis.

**Absorbing layers (K2).** `nyx_sms::axisym` assembled the absorbing layers of the axisymmetric solid with the stretches `s_r` and `s_z` only: the gradients were divided by them and the cell integral multiplied by `s_r s_z`, but the `1/r` terms of the strains (`e_tt`, `g_zt`, `g_rt`, the `m / r` terms) and the weight `r` kept the real radius. In a layer that stretches `r` this is not the analytic continuation of the equations, so the layer is not perfectly matched in the cylindrical case (finding G2 of the plan). The electromagnetic kernel does it right (`AxiMaterial::with_pml` in `nyx-ems/src/axisym.rs`, ADR 46: `r~` and `s_theta = r~ / r`). The solid kernel had no layer test at all before this one.

## Decision (K1): the axis

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

## Decision (K2): the cylindrical stretch

1. **Same stretch as the EM kernel.** Inside a layer whose `x` direction (the radial one) is stretched, the radius becomes the complex `r~ = int s_r dr` counted from the interface of the layer. With `s_theta = r~ / r` every `1/r` of the strains becomes `1/(s_theta r)` and the weight `r dr dz` becomes `s_theta r s_r s_z dr dz`; in the mass the weight `r` becomes `s_theta r` and the determinant is `s_r s_z s_theta`. `s_z` and `s_r` on the gradients stay.
2. **Where.** `nyx_sms::axisym::cylindrical_stretches(mesh, pml)` returns `[s_r, s_z, s_theta]` per cell: the first two from `femcore::pml::cell_stretches`, the third from `femcore::pml::radial_coordinate` at the cell centre (the integral of the profile, exactly the formula that was inline in `AxiMaterial::with_pml`; moved to femcore as a pure move, `nyx-ems` now calls it, its numbers are unchanged). `assemble_stiffness_axisym` and `assemble_mass_axisym` keep their signatures and read the third component of the stretch they are given (`None` is the unstretched ring); `assemble_cross_section` builds the stretch with `cylindrical_stretches`. A cell outside a layer with a radial stretch has `s_theta = 1` exactly, so every model without such a layer assembles to the same matrices as before.
3. **The test function is not conjugated in the stretch.** The weak form is the analytic continuation of `int conj(eps(v)) : D : eps(u) r dA` in which the conjugation reverses the phase `exp(-i m theta)` of the test function and nothing else. So inside a layer the `m / r` terms of the test side are `+i m / (s_theta r)`, not the conjugate of the trial side's `-i m / (s_theta r)` (that would conjugate `s_theta`). The code uses `beta_h = i m / (s_theta r)`; for `s_theta = 1` it is the old `conj(beta)`. The matrices are complex symmetric inside a layer, as the 2D PML kernel's already are, and the shift-and-invert eigensolve is the one the PML kernel uses.
4. **`s_theta` is the value at the cell centre, as in the EM kernel.** In a layer `r~(r)` is therefore `(r~_c / r_c) r` inside a cell, not exactly the radial integral: the error is of order `h / r` per cell and falls with the mesh (measured below). Evaluating `s_theta` at the quadrature points would remove it; not done so that the two kernels read alike (plan step 4).

## The deliberate difference from Bindel and Govindjee

As the plan reads the paper (Bindel and Govindjee, their Sec. 3, after the weak form), the Jacobian factor `r` of the axisymmetric case "should not be transformed": the layer stretches the derivatives and the `1/r` of the strains are left to the real radius, which is what nyx did before this change. K2 uses the analytic continuation `r -> r~` (Collino and Monk; the EM kernel's choice) because it is the perfectly matched form for the cylindrical wave equation: the plan expects the two to differ by a reflection at the layer interface of order `lambda / r`. I did not have the paper at hand while writing this; the quotation is the plan's. If paul prefers the paper's form it is one change (`s_theta = 1`), which the T6 test of this file already does through the public stage pieces; no flag was added.

T6 (`tests/axisym_pml.rs`) measures what it makes, on a thin annular wheel on a ring-shaped post on a substrate, absorbing layers left, right and below, order-2 profile, `m = 0`, wheel mode near `f = 0.086` (units `E = rho = 1`, `nu = 0.3`; the wheel sits at `r` 12 to 16, the layers at `r` 2 to 8 and 20 to 26):

| wheel, layers | `Q`, `s_theta = r~ / r` | `Q`, `s_theta = 1` | relative difference |
|---|---|---|---|
| `r` 12 to 16, thickness 3, strength 3 | 121.84 | 122.38 | 4.5e-3 |
| `r` 12 to 16, thickness 6, strength 3 | 122.18 | 121.95 | 1.9e-3 |
| `r` 9 to 13, thickness 3, strength 3 | 76.90 | 77.55 | 8.4e-3 |
| `r` 9 to 13, thickness 3, strength 6 | 77.15 | 78.13 | 1.3e-2 |

So on this model the two forms agree to about 1 %, more where the layers are nearer the axis, and the difference is of the order of the insensitivity to the layer itself (next section), not larger. The standard candle A1 has the disk and its layers nearer the axis and a different wave-length to radius ratio; its report must state which form was used and may repeat T6 on it.

## Consequences

* Solid disks, posts and substrates can be solved; the radial contour modes of a thin disk now have a closed-form check.
* For `m = +-1` the axis is regular only in the weak sense (the tie is a penalty); rigid translation and tilt, which satisfy it exactly, come out at zero frequency.

## What was verified

**K1, the axis.** Tests `crates/nyx-sms/tests/axisym_axis.rs` (second-order elements; the printed errors are the measured ones, see the status section of `docs/RESTRUCTURE_STATUS.md`): T1 torsion of a free solid cylinder (`u_theta = J1(kr)`, zeros of `J2`), T2 radial plane-strain modes with `u_z = 0` on the ends (`u_r = J1(kr)`, `k R J0 / J1 = (1 - 2 nu) / (1 - nu)`), T3 the radial contour modes of a thin disk (`lambda J0 / J1 = 1 - nu`) and `m = 2` stability under refinement, plus `m = 1` rigid modes and the refusal of a mesh that crosses the axis.

**K2, the stretch.**

All numbers are printed by `crates/nyx-sms/tests/axisym_pml.rs` (`-- --nocapture`), second-order elements; T4 to T6 use off-axis meshes, T5b and T6b the topology of A1 (all touching the axis). Units `E = rho = 1`, `nu = 0.3` unless said.

* **Unit stretch changes nothing** (T4a): a layer with strength 0 and `re = 1` on a ring gives the unstretched frequencies to 1e-9 at `m = 0` and 2. All pre-existing sms and ems tests pass unchanged (the old numbers are bit-for-bit as far as the tests can see: no golden was touched, `ring_model_axes` included).
* **A real radial stretch is the longer ring** (T4b): `re = 1.5` on a ring of width 0.3 at `r = 1` gives the frequencies of the ring of width 0.45, `m = 0` and 3, worst relative error 1.1e-2, 7.2e-4, 3.0e-4 for 6 x 2, 12 x 4 and 24 x 8 cells. This is the per-cell `s_theta` of decision 4: it falls with the mesh and is an upper bound for the effect at a real layer.
* **Large radius** (T4): a clamped strip `r` in `[R0, R0 + L]` with a `+x` layer (strength 0.2, order 1), `R0 = 1000 L`, `m = 0`, against the plane-strain PML kernel on the same strip at the origin: the complex frequencies of the first four in-plane modes agree to 2.1e-4, 6.7e-5, 2.3e-5, 8.1e-5 relative (`L / R0 = 1e-3`; the torsional modes are left out by their `u_theta` share).
* **Layer insensitivity** (T5): the wheel mode near 0.086, `Q = 121.84` with layers of thickness 3 and strength 3 (647 cells), 122.18 with thickness 6 (1161 cells) and 122.19 with strength 6: changes of 2.9e-3 in `Q` and 2.9e-5 in `Re f`, against the 5 % the plan asks for. The wheel holds 95 to 96 % of the mass of the mode.
* **K1 and K2 together, a disk on a post on a substrate, all touching the axis** (T5b, T6b, merged tree; the disk is `r` 0 to 8, 0.4 thick, the post `r` 0 to 1, the substrate `r` 0 to 20 and 3 deep, layers to the right and below, `E = rho = 1`, `nu = 0.3`, `m = 0`; the contour mode is within 10 % of the free disk's closed form): `Q = 2125.4` with layers of thickness 3 and strength 3 (830 cells), 2028.0 with thickness 6 (1310 cells), 2081.9 with strength 6; `Re f` changes by 3.2e-5 at most, `Q` by 4.6 % for the doubled thickness and 2.0 % for the doubled strength (the first is under the 5 % of the plan, with little margin on this coarse mesh). The untransformed form (`s_theta = 1`) gives `Q = 2126.1`, 2043.7, 2103.6: relative differences of 3.1e-4, 7.8e-3 and 1.0e-2, of the order of the insensitivity to the layer itself. So the two forms differ by about 1 % in `Q` here as they do on the annular wheel (T6 above), and the form of Bindel and Govindjee is not needed to reproduce their anchor loss to that accuracy.
* **The study on a solid disk that touches the axis** (`axisym_study.rs`): a thin free disk (`R = 1 mm`, `t = 20 um`) through the envelope, the typed problem and the library entry gives the same frequencies, and the contour mode is `1.774722e6 Hz` against the closed form `1.774765e6 Hz` (2.4e-5 relative); the log names the axis condition.
* **The whole path with the axis and the layers** (engine example 17 and the MCP test): a polysilicon-like disk on a post on a substrate, built by commands, gives the contour mode at 4.505e8 Hz (free disk 4.619e8) with `Q = 1733` and a layer fraction of 0.001, while the modes of the layers have layer fractions of 0.5 to 0.6 and `Q` of 2; on a mesh refined by 0.7 the `Q` of that mode changes by 9 % (the example is deliberately coarse).


## What was not verified

* The exact tie for `m = +-1`: converge-rate of an `m = +-1` mode that touches the axis was not measured beyond the rigid modes.
* Anisotropic materials on the axis: the kernel averages the stiffness about the axis (`average_about_axis`), as before.

* **A layer next to the axis.** T4 to T6 are off-axis and T5b and T6b have their layers far from it; whether `s_theta` inside a layer behaves near `r -> 0` does not arise, because a layer with a radial stretch lies away from the axis and a layer towards the axis (`x = minus`) is refused by the study (ADR 68).
* A complex stretch at `m != 0` is checked only through the unit and real stretch above and by T5/T6 at `m = 0`; there is no closed form for a complex stretch at `m != 0`.
* The `x = minus` direction of the kernel (the left layers of T5) is exercised only through the insensitivity of `Q`, with no separate check.
* The size of the reflection at the interface, which is the quantity the choice between the two forms is about, is not measured directly: T6 compares quality factors of one mode, not the reflection coefficient of a plane wave.
* The 3D-equivalent and the low-reflecting boundary (finding D1 of the plan: `lrb::damping_matrix` has no `r` weight) are untouched.

## Code

`crates/nyx-sms/src/{axisym,cross_section}.rs`, `crates/nyx-sms/tests/axisym_axis.rs`. Part 2: `crates/nyx-sms/src/axisym.rs` (`cylindrical_stretches`, the stiffness and mass), `crates/nyx-sms/src/cross_section.rs` (one line), `crates/femcore/src/pml.rs` (`radial_coordinate`, moved from `crates/nyx-ems/src/axisym.rs`), `crates/nyx-sms/tests/axisym_pml.rs`.

## Contract

None of its own: `CrossAssembly` gains the public fields `axis_dofs` and `pml_cells` (the latter from ADR 68); `kernel_stages` 1.4 (ADR 68) covers the stage API.
