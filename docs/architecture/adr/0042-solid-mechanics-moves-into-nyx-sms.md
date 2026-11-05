# 42. Solid mechanics moves into `nyx-sms`

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-sms`). Applies the pattern of ADR 41 (checklist of six steps) to `sms`; `dfs` was the reference.

## Context

Solid mechanics lived in `femcore` as `sms` (static and real eigenfrequencies), `sms_pml` (absorbing layers), `sms_lrb` (low-reflecting boundaries), `sms_floquet` (periodic cells), `sms_axisym` and `sms_kz` (bodies of revolution). They read `femcore::problem::ProblemSpec`, returned the old `SolveResult` (flat fields, a `disp` array) and were called by `femcore::solve::solve`, by about fifteen test files and by three benchmark examples. The crate `nyx-sms` was a skeleton whose `solve` answered `NotYetMoved`.

## Decision

1. **Typed problem.** `nyx_sms::Problem` has the field names of the old spec that the solid kernels read (`order`, `study`, `materials`, `bcs`, `plane`, `thickness`, `thermal`, `solver`, `inactive`, `periodic`), so the kernels changed only in types. `Study` is `Static` or `Eigenfrequency { n_modes, f_search }`; `BcKind` is `Fixed`, `Load`, `LowReflecting`; a perfectly matched layer is `materials.pml` and a Floquet pair is `periodic`, as before. `ThermalSource` is `Uniform(T)` or `Imported(Vec<f64>)`: the temperature per dof that another physics exported. The old `FromHeat` plus `SolveState::heat_t` is the bridge's business (below); the kernel has no state.
2. **`decode`** reads the envelope: the physics options `plane`, `thickness`, `thermal` (`off`, `uniform`, `imported`), `temperature`, `t_ref`, `inactive` (region names to domain numbers) and the new `solver`; `fixed` constrains the components that are given (`ux`, `uy`, `uz`, absent means free), `load` is a traction vector, `pml` is a domain condition, `periodic` takes the source from `on` and the destination from its own selector `dest`. `imported` takes the role `temperature` from `ctx.inputs` and refuses a run without it. `pml`, `low_reflecting` and `periodic` belong to the eigenfrequency study and `load` to the static one; the other combination is refused in `decode` with a message that names the study (the old kernel refused only the PML one, and ignored the others silently).
3. **Module layout.** `solid` (was `sms`: setup, stiffness, static, real modes, fields), `pml`, `lrb`, `floquet`, `axisym`, `cross_section` (was `sms_kz`, see ADR 43). Every kernel keeps its three public stages (ADR 9); `solve_solid` composes them and chooses the variant exactly as `femcore::sms::solve_solid` did. `run(mesh, &Problem, &Arc<Progress>)` checks the numbers (`Problem::check_finite`, the same messages as the old `validate::check_inputs`, including the refusal of a layer that absorbs along y on a cross-section), solves under the progress handle and appends the stage log. `Physics::solve` is `validate_spec`, `validate`, `decode`, `run`.
4. **Results.** The displacement of a real study (static, real modes) is one vector field `u` with the components of the mesh (`x`, `z` on a cross-section; `x`, `y`, `z` in 3D); its values are the displacement per dof, interleaved, so `Solution::flat` still gives `|u|`, `ux`, `uz`, bit for bit what the old fields were. `von Mises` and the stresses `sxx`... are scalar fields. A complex-frequency mode (layer, Floquet pair, low-reflecting boundary) has the scalar fields `|u|` (modulus, largest component real and of unit modulus) and `Re ux`, `Re uz` (`Re uy` in 3D), and the complex components per cell. They are scalars because a vector `Re u` would add a flat name `|Re u|` that the pre-restructure results and the golden files do not have. The window bridge reads `disp` from `u` or from `Re u{x,z,y}`. The descriptor's `fields` therefore lists the real-mode field `u`; the shape of every variant is pinned by `physics_sms_result.v1.json`.
5. **Descriptors** (`physics_sms` 1.2): the physics option `solver` is added (it was `ProblemSpec::solver`); the property `stiffness` is removed because an elastic matrix has no scalar form and the envelope carries scalars only (ADR 41 item 2): anisotropic stiffness and expansion reach the kernel through `materials.aniso` of the typed problem until the material files of phase 6.
6. **Bridge.** `nyx_engine::legacy::sms` converts the old problem (same mesh numbers), maps `ThermalSource::FromHeat` to `Imported(state.heat_t)` and keeps the old message when there is none (the length check against the mesh is in the kernel), runs `nyx_sms::run`, converts the result back with `result_back` and fills `disp` (ADR 41 item 6).
7. **Tests.** Goldens `solid_*` (six) and `ring_model_axes` moved to `crates/nyx-sms/tests/golden/` and pass without being rewritten, at the same tolerances (1e-8 static, 1e-7 complex). `solid_static_thermal_from_heat` needs dfs and sms and lives with the coupled test in `nyx-engine/tests`. The solid tests of `femcore/tests` moved with their numbers and assertions: `solid`, `solid_pml`, `anisotropic` (solid half), `stages` (the four solid stage tests), `ring_frame` (all but the optical ring test), `axisym_closed_form` (the ring half), `amg_solve`, `cancel_and_log`, `determinism`, `finite_guard` (input checks), `frame`, `frame_roundtrip`, `geom_meshes` (plate with a hole), `golden_frames` (solid cases) and `split` (whole: its oracle is a solid bar). `femcore` keeps what is not solid. New: envelope against typed problem bit for bit, the Floquet pair and the layer through the envelope, an imported temperature against the same uniform one, refusals, the stage API, the result-shape exemplar.
8. **Benchmarks** `eig_beam`, `mem_profile`, `profile3d` drive the solid kernel through the bridge and moved to `nyx-engine/examples` with `bench` (ADR 41 item 9); `tools/baseline.sh` runs them from there.

## Consequences

* Positive: `nyx-sms` is a physics crate on the same footing as `nyx-dfs`; the engine can solve `sms` problems and chain `dfs` into `sms` through `Ctx::inputs`; goldens unchanged.
* Negative: the descriptors promise `u` while a lossy mode gives `Re u{x,z}` (point 4); `sizing` for `sms` is still the skeleton's (no wavelengths), so `physics_size::spec_sizing` stays until the ems move is done (ADR 41, open item 3).
* Open: whether a complex mode should offer a vector `Re u` once the flat naming is replaced in phase 7; native non-simplex elements (ADR 41 item 5).

## Code

`crates/nyx-sms/src/*`, `crates/nyx-engine/src/legacy/sms.rs`, tests as listed.

## Contract

`physics_sms.v1.json` 1.2, `physics_sms_result.v1.json` 1.0 (new), `kernel_stages` 1.3.
