# 68. The axisymmetric solid gets a study

Date: 2026-11-14

## Status

Proposed (capability R6, the merge of the three workstreams on `capability-r6-sms-axisym`). Amends ADR 43: the library entry stays, the study is added. The axis conditions and the cylindrical stretch of the kernel are ADR 69. Numbered 68 because the A3 branch (`restructure-10-a3-snowflake`) already holds ADR 66 and 67.

## Context

Benchmark A1 (Bindel and Govindjee 2005, a disk on a post anchored in a substrate) was blocked on "no axisymmetric `sms`". The kernel exists and is tested (ADR 43: `solve_solid_axisym`, ring closed forms, crystal-ring frame checks), but it was a library function: `Sms::studies()` listed `static` and `eigenfrequency` only, and `run` always called `solve_solid`. ADR 43 left the study to the owner of the window ("nothing needs it before phase 7"); `CAPABILITY_PLAN.md` R6 and the A1 report ask for it, and paul stated that the support exists and only needs wiring. `nyx-ems` is the precedent (ADR 46: study `axisymmetric_modes`).

Everything downstream is generic: `nyx-engine` never names a physics, `set_study` takes any kind the descriptors list, `add_pml` plus `add_pml_conditions` write `x`/`z` layers of a 2D model, the window draws the Study tab from the descriptors, and the MCP tools are the engine commands (ADR 57). So the change is inside `nyx-sms` and the contracts of `sms`.

## Decision

A study of `sms`, parallel to the `axisymmetric_modes` study of `ems`.

* **The study.** `kind` `axisymmetric_eigenfrequency`, name "Axisymmetric eigenfrequencies", `Dims::Two`, parameters `n_modes` (integer, default 6, min 1), `f_search` (Hz, default 0 = lowest modes, range >= 0) and `azimuthal_order` (integer `m`, range -100000 to 100000, same wording and range as `ems`; the displacement goes as `exp(-j m theta)`, `m > 0` runs along +y). The mesh is the (r, z) half plane with `x = r`; the displacement components are `ux` radial, `uz` axial, `uy` azimuthal (the order of `fixed`).
* **The typed problem.** `Study::AxisymmetricEigenfrequency { n_modes, f_search, azimuthal_order: i64 }`; `Study::kind()` returns the string, `Study::is_eigen()` is true for it, `decode` reads it, `check_finite` covers `f_search`.
* **What the study refuses** (`Problem::check_axisymmetric`, called by `decode` and by `run`, so the typed problem is held to the same rules as the envelope): a mesh that is not 2D (`Sms::validate` and `run`); a mesh with a point at `x < -1e-9` of its width ("the (r, z) half plane must lie in x >= 0", `solve_axisymmetric`); a boundary load; a Floquet-periodic pair (no translation to repeat); a low-reflecting boundary, with code `unsupported`, because `lrb::damping_matrix` has no `r` weight (rough spot D1: the damping of the axisymmetric kernel is wrong today and is not made reachable through a new door); a layer along `y` (`check_wafer`) or with `x = minus` ("a layer cannot absorb towards the axis"); a fixed value other than zero. Absorbing layers (`pml` with `x = plus` and `z`) are allowed.
* **What it ignores.** `plane`, `thickness`, `thermal` and `solver` are not read; when a non-default `plane` or `thickness` is present the log says so. `inactive` is honoured. `Sms::sizing` treats the new kind like `eigenfrequency` (a wavelength per region from `f_search`), and `Sms::fields` lists `|u|`, `Re ux`, `Re uz`, `Re uy`, `|uy|`.
* **The result.** One `Solution` per mode sorted by `Re f`; `value = Re f` (Hz) and **`eigenvalue = [Re f, Im f]`** on every mode (also without layers, where `Im f` is rounding noise of the complex solver, below 1e-8 of `Re f` in the test). Node fields `|u|`, `Re ux`, `Re uz`, `Re uy`, `|uy|` (phase rotated so the largest component is real, unit maximum). Cell fields are the complex `ux`, `uz`, `uy` as before. Scalars per mode `N`: `fN (Hz)` (a number when lossless, the text `re +imi` when lossy), `QN = Re f / (2 |Im f|)` when lossy, and **`layer fraction N`**, the share of `int rho |u|^2 r` that lies in absorbing layers (0 without layers; from the real mass weights `rho r` of the layer cells, not from the stretched ones). It lets a host or a test tell a resonator mode from a layer mode, the complaint of findings A1 and A4. The kernel gives `Im f > 0` for the layers of `sms` (measured in `tests/axisym_study.rs`); the study doc says so.
* **Names (rough spot D2).** The node fields of the shared stage `cross_section_fields` were `ux`, `uz`, `uy` although they are real parts of complex modes; they are now `Re ux`, `Re uz`, `Re uy`, as the 2D PML kernel and `nyx-ems` name them (contract `kernel_stages` 1.4 says so). The per-cell complex fields keep their names.
* **Entries.** `cross_section::solve_axisymmetric(mesh, problem)` is the study's entry (the name `nyx-ems` uses); `solve_solid_axisym(mesh, problem, m)` stays for the ring tests. `assemble_cross_section` accepts both `Study::Eigenfrequency` (with `m` from the argument) and the new study. `run` dispatches on the study; `solve_solid` forwards the new study to `solve_axisymmetric`.
* The `Physics` trait, `ProblemSpec`, `SolveResult`, the registry and the engine commands do not change shape (`physics_api` stays 1.1).

## Consequences

* The study is reachable from the envelope, the engine, the window (Study tab) and MCP with no change outside `nyx-sms`; example project and MCP proof are the engine step of the plan.
* A body that touches the axis is solved with the regularity conditions of ADR 69 and its layers are the cylindrical stretch of ADR 69, so a disk on a post on a substrate is a model of this study (engine example 17, the MCP proof). The study's own mesh rule stays `x >= 0`.
* Not covered, and recorded: the low-reflecting boundary with `r` weight (D1), the azimuthal tie at `m = +-1` on the axis, a complex shift (`f_search` stays real, enough for Q of 1e3 to 1e5), the `periodic` Floquet pair.
* Contracts: `physics_sms` 1.3 (additive: the study and its fields), `physics_sms_result` 1.1 (additive: the shapes `axisymmetric_eigenfrequency` and `axisymmetric_eigenfrequency_pml`), `kernel_stages` 1.4 (additive, with the rename of D2).

## What was verified

`crates/nyx-sms/tests/axisym_study.rs`: on a solid disk that touches the axis the study gives the contour mode of the closed form to 2.4e-5; the engine example `17_disk_anchor_loss` solves through engine commands (test T10), `convergence_check` compares the complex eigenvalues (T11) and the MCP server over stdio describes the study and returns `Q1` for the same model (T12); the envelope (`Sms::solve` of a `ProblemSpec`), the typed problem (`nyx_sms::run`) and the library entry (`solve_solid_axisym` with the same `m`) give identical frequencies at `m = 0` and `m = 2`; the refusals above each fail with their own message and a good problem passes; an absorbing layer on the outer part of a ring section gives complex frequencies with `Q = Re f / (2 |Im f|)` and layer fractions between 0 and 1 (0.79 to 0.96 for the four lowest modes of the test ring); descriptors, fields and sizing. `stages.rs` checks that the stages compose to `solve_axisymmetric`, `determinism.rs` that two runs give the same bits (a layer, `m = 1`), `cancel_and_log.rs` the stage log and `cancelled`, `contracts.rs` the descriptors exemplar, the result-shape exemplar (two new entries) and an example spec of the new kind.

## What was not verified

The quality factor of a model that matches a published anchor-loss number (benchmark A1 is a later step: `crates/nyx-engine/tests/benchmark_a1.rs` is an ignored skeleton until paul has agreed the material and the criteria); the sign convention of `Im f` is observed, not derived; the window's Study tab was not looked at (it follows from `forms.rs` being descriptor-driven, ADR 57).

## Code

`crates/nyx-sms/src/{problem,lib,cross_section,solid}.rs`, `crates/nyx-sms/tests/{axisym_study,stages,contracts,determinism,cancel_and_log}.rs`, `crates/nyx-engine/tests/{examples,make_examples,benchmark_a1}.rs`, `crates/nyx-mcp/tests/stdio.rs`, `examples/17_disk_anchor_loss.nyx.json`, `docs/architecture/contracts/{physics_sms,physics_sms_result}.v1.json`, `docs/architecture/contracts/VERSIONS.md`.

## Contract

`physics_sms` 1.3, `physics_sms_result` 1.1, `kernel_stages` 1.4.
