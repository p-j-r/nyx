# Kernel port to wafer axes: status

Plan: `nyx_frames/KERNEL_PORT_PLAN.md` in the project files. Oracle: `crates/femcore/tests/golden_frames/*.json`, frozen from the kernels before they were ported (`tests/golden_frames.rs`; `NYX_WRITE_GOLDEN=1` rewrites them, only do that on a kernel known to be right).

| Phase | Kernels | State |
|---|---|---|
| 0 | `section` module, golden files | done |
| 1 | `scalar` (electrostatics, heat) | done, adapter bypassed, golden files match to 1e-9 |
| 1 pilot / 3 | `emw`, `emw_pml` (real, tilted-crystal and PML paths), `Physics::Emw` on a cross-section | done, golden files match to 1e-8 / 1e-7 |
| 2 | `solid`, `solid_pml`, `solid_floquet`, `solid_lrb` on a cross-section without `out_of_plane` | done, 7 golden cases match to 1e-8 / 1e-7 (the low-reflecting kernel has no golden: its test fails on this machine before and after) |
| 4 | `solid_kz` (`solve_solid_kz`: a cross-section with `out_of_plane`) | done, 4 golden cases (tilted sapphire, component constraints, PML along x and z) match to 1e-8 / 1e-7; the model-frame `solve_solid_kz_full` stays until phase 5 for Brillouin |
| 5 | `brillouin` (straight waveguide: optical modes, photoelastic overlap, moving boundary, gains, viscosity) | done, 3 golden cases match to 1e-6 (silicon wire, cubic and tilted crystal, backward and forward); the ring (`ring_radius`) still goes through the adapter until phase 6 |
| 6 | axisymmetric family (`solid_axisym` through `solve_solid_axisym`, `emw_axisym`, `ring_brillouin`, the ring in `brillouin::solve_ring`) | done, 2 ring Brillouin golden cases match to 1e-6 (the Q of a lossless ring mode is round-off and is not compared); the ring frame tests pass; `m > 0` is propagation along +y (signed test) |
| 7 | delete the adapter | open |

The solid kernels work in the Voigt order of the axes (x, z, y), `section::plane_stiffness` (so their in-plane block is (xx, zz, xz) and plane stress condenses yy); boundary components, loads, expansion tensors, PML directions and field names go through `section` too.

`section` is the one place with axis numbers: `plane_tensor` (x, z, y; for diffusion and in-plane elasticity), `local_tensor` and `LOCAL_AXES` (lateral, normal, propagation; for the mode solvers), `component_names`, `wafer_point`. Everything propagates along `+y` (decided 2026-11-02): the mode solvers carry `exp(-j beta y)` (`PROP_SIGN = +1`; their local triple (x, z, y) is left-handed, which only reverses the reading of the polarization sense, `TimeConvention::MinusI`). Neff did not change, only the sign of `Ey` (golden files regenerated for that). The Brillouin kernels still use `LEGACY_PROP_SIGN = -1` through `brillouin::solve_optical` until phase 5, when the pump moves to `+y` as well.

The Brillouin kernels run in the axes (x, z, y) of the cross-section with the wave along `+y` (optical fields, elastic modes, photoelastic tensor `brillouin::photo_from_spec`, viscosity in `plane_stiffness` order); `mirror_axis` is a mesh axis (0 = x, 1 = z), so it needed no translation, only its documentation. The golden files were written with the kernel that still propagated along -y, on the mirror image of the crystal (half turn about y, which even-rank tensors cannot tell from a mirror), so they check the +y direction as well: gains, Q, frequencies and group indices are the same, and only the sign of the displacement component along y differs (the mirror image of the displacement, written into the files).

The ring kernels run in the ring's local wafer frame (r, theta, z) = (x, y, z), which is right-handed with no sign flip: the displacement blocks are (ur, uz, utheta) = (ux, uz, uy), the crystal data are averaged about the ring axis z, the layers are x (radial) and z. `exp(-i m theta)` with `m > 0` is a wave along +theta = +y: the spatial dependence follows from a rigid translation (test in `tests/ring_frame.rs`), and the time convention e^{j omega t} (the same as the straight kernels) from the sign of the optical power along theta, `emw_axisym::power_along_theta`, which is positive for `m > 0` and negative for `m < 0` at the same frequency. The ring goldens were written with the kernels of main (`1f2cb58`).

The kz kernel carries the wave `exp(-j k y)` with displacement blocks (x, z, y); the equations did not change (a reflected third axis with a reflected wave number is the same eigenproblem), only the stiffness order, the constraint components and the names. The golden comparison accepts an eigenvector with all its signed fields negated: when two entries of a mode tie in modulus its overall sign is round-off.

Known: `solid_pml::low_reflecting_boundary_gives_the_analytic_damped_bar_modes` fails on this machine's cloud container before and after the port (Re f1 = 4e-7, a rigid-body mode); not caused by the port.
