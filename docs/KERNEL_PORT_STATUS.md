# Kernel port to wafer axes: status

Plan: `nyx_frames/KERNEL_PORT_PLAN.md` in the project files. Oracle: `crates/femcore/tests/golden_frames/*.json`, frozen from the kernels before they were ported (`tests/golden_frames.rs`; `NYX_WRITE_GOLDEN=1` rewrites them, only do that on a kernel known to be right).

| Phase | Kernels | State |
|---|---|---|
| 0 | `section` module, golden files | done |
| 1 | `dfs` (electrostatics, heat) | done, adapter bypassed, golden files match to 1e-9 |
| 1 pilot / 3 | `ems`, `ems_pml` (real, tilted-crystal and PML paths), `Physics::Emw` on a cross-section | done, golden files match to 1e-8 / 1e-7 |
| 2 | `sms`, `sms_pml`, `sms_floquet`, `sms_lrb` on a cross-section without `out_of_plane` | done, 7 golden cases match to 1e-8 / 1e-7 (the low-reflecting kernel has no golden: its test fails on this machine before and after) |
| 4 | `sms_kz` (`solve_solid_kz`: a cross-section with `out_of_plane`) | done, 4 golden cases (tilted sapphire, component constraints, PML along x and z) match to 1e-8 / 1e-7; the model-frame `solve_solid_kz_full` stays until phase 5 for Brillouin |
| 5 | `sbs` (straight waveguide: optical modes, photoelastic overlap, moving boundary, gains, viscosity) | done, 3 golden cases match to 1e-6 (silicon wire, cubic and tilted crystal, backward and forward); the ring (`ring_radius`) still goes through the adapter until phase 6 |
| 6 | axisymmetric family (`sms_axisym` through `solve_solid_axisym`, `ems_axisym`, `sbs_ring`, the ring in `brillouin::solve_ring`) | done, 2 ring Brillouin golden cases match to 1e-6 (the Q of a lossless ring mode is round-off and is not compared); the ring frame tests pass; `m > 0` is propagation along +y (signed test) |
| 7 | delete the adapter (`frame_adapter` removed: `to_model_frame`, `rename_wafer_fields`, `WAFER_TO_MODEL`, `Kernel`, the model-frame tests; `check_wafer_input` and `axes_of_3d_problem` moved into `section`) | done |

The solid kernels work in the Voigt order of the axes (x, z, y), `section::plane_stiffness` (so their in-plane block is (xx, zz, xz) and plane stress condenses yy); boundary components, loads, expansion tensors, PML directions and field names go through `section` too.

`section` is the one place with axis numbers: `plane_tensor` (x, z, y; for diffusion and in-plane elasticity), `local_tensor` and `LOCAL_AXES` (lateral, normal, propagation; for the mode solvers), `component_names`, `wafer_point`. Everything propagates along `+y` (decided 2026-11-02): the mode solvers carry `exp(-j beta y)` (`PROP_SIGN = +1`; their local triple (x, z, y) is left-handed, which only reverses the reading of the polarization sense, `TimeConvention::MinusI`). Neff did not change, only the sign of `Ey` (golden files regenerated for that). The Brillouin kernels still use `LEGACY_PROP_SIGN = -1` through `brillouin::solve_optical` until phase 5, when the pump moves to `+y` as well.

The Brillouin kernels run in the axes (x, z, y) of the cross-section with the wave along `+y` (optical fields, elastic modes, photoelastic tensor `brillouin::photo_from_spec`, viscosity in `plane_stiffness` order); `mirror_axis` is a mesh axis (0 = x, 1 = z), so it needed no translation, only its documentation. The golden files were written with the kernel that still propagated along -y, on the mirror image of the crystal (half turn about y, which even-rank tensors cannot tell from a mirror), so they check the +y direction as well: gains, Q, frequencies and group indices are the same, and only the sign of the displacement component along y differs (the mirror image of the displacement, written into the files).

The ring kernels run in the ring's local wafer frame (r, theta, z) = (x, y, z), which is right-handed with no sign flip: the displacement blocks are (ur, uz, utheta) = (ux, uz, uy), the crystal data are averaged about the ring axis z, the layers are x (radial) and z. `exp(-i m theta)` with `m > 0` is a wave along +theta = +y: the spatial dependence follows from a rigid translation (test in `tests/ring_frame.rs`), and the time convention e^{j omega t} (the same as the straight kernels) from the sign of the optical power along theta, `emw_axisym::power_along_theta`, which is positive for `m > 0` and negative for `m < 0` at the same frequency. The ring goldens were written with the kernels of main (`1f2cb58`).

The kz kernel carries the wave `exp(-j k y)` with displacement blocks (x, z, y); the equations did not change (a reflected third axis with a reflected wave number is the same eigenproblem), only the stiffness order, the constraint components and the names. The golden comparison accepts an eigenvector with all its signed fields negated: when two entries of a mode tie in modulus its overall sign is round-off.

Known: `solid_pml::low_reflecting_boundary_gives_the_analytic_damped_bar_modes` fails on this machine's cloud container before and after the port (Re f1 = 4e-7, a rigid-body mode); not caused by the port.

## Physics codes (2026-11-03)

Every physics module carries a three-letter code, chosen by Paul:

| Code | Physics | Modules |
|---|---|---|
| `dfs` | scalar diffusion (electrostatics, heat) | `dfs` (was `scalar`) |
| `sms` | solid mechanics | `sms`, `sms_pml`, `sms_floquet`, `sms_lrb`, `sms_kz`, `sms_axisym` (were `solid*`) |
| `ems` | EM waves | `ems`, `ems_pml`, `ems3d`, `ems_axisym` (were `emw*`) |
| `sbs` | stimulated Brillouin scattering | `sbs`, `sbs_ring` (were `brillouin`, `ring_brillouin`) |

`solvers::scalar` (the solver library's scalar trait) is unrelated and keeps its name.

## Kernel stages: assemble, solve, fields (2026-11-03)

Every kernel is three functions; `solve_*` is their composition, and the golden files are unchanged, so the split changed no number.

1. **assemble**: matrices on the free dofs, constraints, loads, the material data later stages need.
2. **solve**: the raw result (complex eigenpairs, a solution vector, the states of a transient run). The backend is chosen here (`spec.solver`, or the eigensolver the stage calls), so a `femcore::solvers` adapter replaces one function per kernel.
3. **fields**: the plot fields (node-averaged, as before) and `Solution::cell_fields`: the same components per cell and local node, not averaged, with wafer names and complex values where the quantity is complex (`fields::CellFields`, indexed `cell * nloc + a`). The node average of a cell field is the plot field of the same name (checked on every golden case).

| Kernel | assemble | solve | fields |
|---|---|---|---|
| dfs electrostatics | `assemble_electrostatics` | `solve_potential` | `electrostatics_fields` (`V`, `Ex`, `Ez`) |
| dfs heat | `assemble_heat` | `solve_heat_states` | `heat_fields` (`T`, `qx`, `qz`) |
| ems 2D real | `assemble_modes` | `solve_mode_pairs` | `mode_fields_cells` (`Ex`, `Ez`, `Ey`) |
| ems 2D complex / PML / tilted | `assemble_modes_complex` | `solve_mode_pairs_complex` | `mode_fields_complex_cells` (complex `Ex`, `Ez`, `Ey`) |
| sms kz and axisymmetric | `assemble_cross_section` | `solve_cross_section_modes` | `cross_section_fields` (complex `ux`, `uz`, `uy`) |
| ems axisymmetric | `assemble_axisym` | `solve_axisym_pairs` | `mode_samples` (unchanged) |
| sms plain static / eigen | inline in `solve_solid` (`setup`, `assemble_solid`, loads) | static: `solve_spd`; eigen: `solve_real_modes` (AMG + LOBPCG kept) | `fields_of` (`ux`, `uz`, stresses, `von Mises`) |
| sms PML, Floquet, low-reflecting | inline | inline | `complex_mode_fields` (complex `ux`, `uz`) |
| ems 3D | inline in `solve_emw_eigen` | inline | cell `Ex`, `Ey`, `Ez` (complex) |

The last two rows have the fields stage and the cell fields, but their assemble and solve are still one function each (they share large local state with the Floquet reduction, the low-reflecting damping and the polarization analysis); they are the next to cut.
