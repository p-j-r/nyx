# 80. Modes carry what identifies them, and a convergence check follows a mode

Date: 2026-11-15

## Status

Proposed (capability R11, bookmark `capability-r11-sms-mode-scalars`, on `capability-r6-sms-axisym`). Amends ADR 68 (the result of the axisymmetric study gains scalars) and ADR 61 (`convergence_check` gains a way to pair solutions). Numbers 70 to 79 are reserved by other capability plans.

## Context

Benchmark A1 (Bindel and Govindjee 2005) was run through the MCP server as a product, and four things in the way were not physics:

1. A lossy mode's frequency was only the text `fN (Hz)` (`"6.93e8 +2.56e7i"`), which a host had to parse, while `QN` was a number.
2. Nothing said which family a mode belongs to. At `m = 0` the torsional modes (azimuthal displacement only) are solved together with the radial and axial ones, a flexural mode sits near the contour mode, and `layer fraction N` only separates the modes of the layers. The benchmark summed squares of node values, which is not the energy and needs every field read back.
3. `convergence_check` compares solution `i` with solution `i`. The order of the complex modes of an open resonator changes between meshes (a layer mode moves past a resonator mode), so the check reported a relative difference near 1 for a mode that had moved by 1e-4.
4. The study did not say that `m = 0` holds torsional modes, or that the number of modes returned can differ from `n_modes`.

## Decision

**Scalars per mode of `axisymmetric_eigenfrequency`** (stage `cross_section_fields`), additive:

* `f re N (Hz)` and `f im N (Hz)`: the two parts of the eigenvalue as numbers, for every mode, with or without loss. The scalar `fN (Hz)` keeps its form (a number without loss, text with loss; `physics_sms_result` 1.1), because removing it breaks hosts for no gain.
* `radial share N`, `axial share N`, `azimuthal share N`: the split of `int rho |u|^2 r dr dz` over the three displacement components, taken over the cells outside the absorbing layers. They sum to 1.
* `share in domain <label> N` for each domain outside the layers that has mass: the same integral per domain. They sum to 1. The kernel knows domain numbers only; the engine replaces `domain <label>` by the region name after the solve (`solve::name_domains`), so a host sees `share in disk N`. A caller of the kernel without an engine sees the numbers.

The integrals are computed by `axisym::mode_energy_axisym` in one pass over the cells (the integrand and quadrature of `assemble_mass_axisym`, no stretch), not by one mass matrix per domain. They are mass-weighted, so they are shares of the kinetic energy of the physical displacement. `layer fraction N` is unchanged.

**`convergence_check` pairs solutions** (`convergence::Matching`), additive:

* `match`: `"index"` (default: nothing changes for existing callers) or `"nearest"`. With `nearest` each solution of the present mesh is paired with the nearest eigenvalue of the finer mesh in the complex plane (distance over the larger modulus), and a pair is kept only when the two are each other's nearest (mutual nearest neighbours). No threshold is asked of the caller.
* `track: { f }`: only the solution whose eigenvalue has its real part nearest `f` is compared, against the nearest eigenvalue of the finer mesh. It implies `nearest`; with `match: "index"` it is refused.
* With pairing, the rows are the eigenvalue, the value and the per-mode scalars of each pair, named with the index of the mode on the present mesh; `pairs` lists `{ coarse, refined, distance }`, `unpaired` the solutions of each mesh without a partner; `worst_relative` and `converged` cover the paired modes only, and the note says so. Scalars that do not belong to a mode are not compared in this form.
* A per-mode scalar is recognised by its name: the one-based mode number at the end, before a unit in brackets if there is one (`Q3`, `layer fraction 3`, `f re 3 (Hz)`). This is a naming convention of results, not a field of `SolveResult`.

**Descriptions.** The study doc states the new scalars, that `m = 0` includes torsional modes, and that the modes returned are the converged eigenpairs nearest `f_search` (the block Krylov-Schur run returns the pairs whose residual is below 1e-6: fewer than `n_modes` when some did not converge, more when a block converged past it). This is how the solver behaves, not a defect, so it is described and not changed. The guide card says that Q converges slower than f, to refine at re-entrant corners, and to compare meshes with `match` or `track`.

**`benchmark_a1.rs` is removed** from `crates/nyx-engine/tests`. Benchmarks are not nyx software (paul, 2026-11-14): no benchmark driver and no test asserting a paper's number lives in the repository. The closed-form tests and example 17 remain the regression guards of the study.

## Consequences

* A host picks a mode by numbers: small `layer fraction`, `radial share` and `share in disk` large, nearest a frequency. No field needs to be read.
* Pairing by eigenvalue alone can pair two different modes that cross between meshes; the mutual rule drops most of these and `distance` shows how far a pair moved. Pairing by field overlap would be stricter and needs both fields on one mesh; not done.
* The convention for per-mode scalar names is implicit. A scalar of another physics that ends in a number within the range of the solutions is taken as a per-mode scalar when pairing is asked for. Recorded in `ROUGH_SPOTS.md` (R11a).
* `convergence_check` still refines only the global element size and its rules, not the sizes of `add_size` features (ADR 61 as it was). On the model that asked for the pairing, the A1 disk, the resolution is set by such sizes, so the check follows the right mode but tests the bulk mesh only. The tool description and the guide card say so; scaling size features is recorded as R11b and not done here.
* The same scalars are the natural next step for `ems` (`axisymmetric_modes`, `cavity`), and for the 2D and 3D complex eigenfrequency kernels of `sms` (`pml.rs`, `lrb.rs`), which still give only the text frequency. Not done here.
* Contracts: `physics_sms_result` 1.2, `physics_sms` 1.4, `engine_commands` and `mcp_tools` next minor, `kernel_stages` 1.5 (all additive).

## What was verified

`crates/nyx-sms/tests/axisym_mode_scalars.rs`: `f re` and `f im` equal the eigenvalue for lossless and lossy modes; the shares sum to 1; the torsional mode of a solid cylinder has azimuthal share 1 and the radial plane-strain mode radial share 1 (thresholds in the test); a two-domain body's domain shares sum to 1; a layer has no share. `crates/nyx-engine/tests/examples.rs`: on example 17 the shares appear under the region names and the contour mode is radial and in the disk; `convergence_check` with `nearest` pairs or lists every mode, `track` gives one mode, `index` is unchanged, bad arguments are refused. `crates/nyx-mcp/tests/stdio.rs`: the scalars and the two new arguments over stdio. The A1 reference case rerun on a build of this bookmark gives the numbers of the build before it (recorded with the benchmark, outside the repository).

## What was not verified

Pairing on a case where two resonator modes cross between meshes; the shares for `|m| >= 1` against a closed form (the sums are checked, the split is not); the window (it shows the scalars of a result generically).

## Code

`crates/nyx-sms/src/{axisym,cross_section,lib}.rs`, `crates/nyx-sms/tests/axisym_mode_scalars.rs`, `crates/nyx-engine/src/{convergence,commands,solve}.rs`, `crates/nyx-engine/cards/nyx-instructions.md`, `crates/nyx-engine/tests/examples.rs`, `crates/nyx-mcp/tests/stdio.rs`, `docs/architecture/contracts/`.

## Contract

`physics_sms_result` 1.2, `physics_sms` 1.4, `engine_commands` 1.6, `mcp_tools` 2.3, `kernel_stages` 1.5.
