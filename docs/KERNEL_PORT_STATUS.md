# Kernel port to wafer axes: status

Plan: `nyx_frames/KERNEL_PORT_PLAN.md` in the project files. Oracle: `crates/femcore/tests/golden_frames/*.json`, frozen from the kernels before they were ported (`tests/golden_frames.rs`; `NYX_WRITE_GOLDEN=1` rewrites them, only do that on a kernel known to be right).

| Phase | Kernels | State |
|---|---|---|
| 0 | `section` module, golden files | done |
| 1 | `scalar` (electrostatics, heat) | done, adapter bypassed, golden files match to 1e-9 |
| 1 pilot / 3 | `emw`, `emw_pml` (real, tilted-crystal and PML paths), `Physics::Emw` on a cross-section | done, golden files match to 1e-8 / 1e-7 |
| 2 | `solid`, `solid_pml`, `solid_floquet`, `solid_lrb` on a cross-section without `out_of_plane` | done, 7 golden cases match to 1e-8 / 1e-7 (the low-reflecting kernel has no golden: its test fails on this machine before and after) |
| 4 | `solid_kz` | open |
| 5 | `brillouin` | open; its optical part goes through `frame_adapter::to_wafer_frame` (transitional) |
| 6 | axisymmetric family | open |
| 7 | delete the adapter | open |

The solid kernels work in the Voigt order of the axes (x, z, y), `section::plane_stiffness` (so their in-plane block is (xx, zz, xz) and plane stress condenses yy); boundary components, loads, expansion tensors, PML directions and field names go through `section` too.

`section` is the one place with axis numbers: `plane_tensor` (x, z, y; for diffusion and in-plane elasticity), `local_tensor` and `LOCAL_AXES` (lateral, normal, propagation; for the mode solvers), `component_names`, `wafer_point`. Everything propagates along `+y` (decided 2026-11-02): the mode solvers carry `exp(-j beta y)` (`PROP_SIGN = +1`; their local triple (x, z, y) is left-handed, which only reverses the reading of the polarization sense, `TimeConvention::MinusI`). Neff did not change, only the sign of `Ey` (golden files regenerated for that). The Brillouin kernels still use `LEGACY_PROP_SIGN = -1` through `brillouin::solve_optical` until phase 5, when the pump moves to `+y` as well.

Known: `solid_pml::low_reflecting_boundary_gives_the_analytic_damped_bar_modes` fails on this machine's cloud container before and after the port (Re f1 = 4e-7, a rigid-body mode); not caused by the port.
