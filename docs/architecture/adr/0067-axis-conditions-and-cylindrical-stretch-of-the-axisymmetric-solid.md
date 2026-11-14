# 67 (part 2: K2). The cylindrical stretch of the axisymmetric solid kernel

Date: 2026-11-14

This file is the K2 half of ADR 67 (bookmark `capability-r6-sms-pml`), written under the same file name as the K1 half (bookmark `capability-r6-sms-axis`, axis conditions) so that the two merge as an add/add: at the merge, keep part 1 as it is and append everything below under the heading "K2, the cylindrical stretch", with one Status, Date, Code and Contract block for the whole ADR. The K1 half is not covered here.

## Context

`nyx_sms::axisym` assembled the absorbing layers of the axisymmetric solid with the stretches `s_r` and `s_z` only: the gradients were divided by them and the cell integral multiplied by `s_r s_z`, but the `1/r` terms of the strains (`e_tt`, `g_zt`, `g_rt`, the `m / r` terms) and the weight `r` kept the real radius. In a layer that stretches `r` this is not the analytic continuation of the equations, so the layer is not perfectly matched in the cylindrical case (finding G2 of the plan). The electromagnetic kernel does it right (`AxiMaterial::with_pml` in `nyx-ems/src/axisym.rs`, ADR 46: `r~` and `s_theta = r~ / r`). The solid kernel had no layer test at all before this one.

## Decision (K2)

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

## What was verified

All numbers are printed by `crates/nyx-sms/tests/axisym_pml.rs` (`-- --nocapture`), second-order elements, off-axis meshes only. Units `E = rho = 1`, `nu = 0.3` unless said.

* **Unit stretch changes nothing** (T4a): a layer with strength 0 and `re = 1` on a ring gives the unstretched frequencies to 1e-9 at `m = 0` and 2. All pre-existing sms and ems tests pass unchanged (the old numbers are bit-for-bit as far as the tests can see: no golden was touched, `ring_model_axes` included).
* **A real radial stretch is the longer ring** (T4b): `re = 1.5` on a ring of width 0.3 at `r = 1` gives the frequencies of the ring of width 0.45, `m = 0` and 3, worst relative error 1.1e-2, 7.2e-4, 3.0e-4 for 6 x 2, 12 x 4 and 24 x 8 cells. This is the per-cell `s_theta` of decision 4: it falls with the mesh and is an upper bound for the effect at a real layer.
* **Large radius** (T4): a clamped strip `r` in `[R0, R0 + L]` with a `+x` layer (strength 0.2, order 1), `R0 = 1000 L`, `m = 0`, against the plane-strain PML kernel on the same strip at the origin: the complex frequencies of the first four in-plane modes agree to 2.1e-4, 6.7e-5, 2.3e-5, 8.1e-5 relative (`L / R0 = 1e-3`; the torsional modes are left out by their `u_theta` share).
* **Layer insensitivity** (T5): the wheel mode near 0.086, `Q = 121.84` with layers of thickness 3 and strength 3 (647 cells), 122.18 with thickness 6 (1161 cells) and 122.19 with strength 6: changes of 2.9e-3 in `Q` and 2.9e-5 in `Re f`, against the 5 % the plan asks for. The wheel holds 95 to 96 % of the mass of the mode.

## What was not verified

* **A mesh that touches the axis.** Everything above is off-axis (K1 is not on `main`; the wheel is annular instead of a disk). Whether `s_theta` inside a layer behaves near `r -> 0` does not arise: a layer with a radial stretch lies away from the axis, and a layer towards the axis (`x = minus`) is refused by the study of the plan.
* A complex stretch at `m != 0` is checked only through the unit and real stretch above and by T5/T6 at `m = 0`; there is no closed form for a complex stretch at `m != 0`.
* The `x = minus` direction of the kernel (the left layers of T5) is exercised only through the insensitivity of `Q`, with no separate check.
* The size of the reflection at the interface, which is the quantity the choice between the two forms is about, is not measured directly: T6 compares quality factors of one mode, not the reflection coefficient of a plane wave.
* The 3D-equivalent and the low-reflecting boundary (finding D1 of the plan: `lrb::damping_matrix` has no `r` weight) are untouched.

## Code

`crates/nyx-sms/src/axisym.rs` (`cylindrical_stretches`, the stiffness and mass), `crates/nyx-sms/src/cross_section.rs` (one line), `crates/femcore/src/pml.rs` (`radial_coordinate`, moved from `crates/nyx-ems/src/axisym.rs`), `crates/nyx-sms/tests/axisym_pml.rs`.
