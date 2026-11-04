# Restructure: status

Handoff file between working sessions (see `CLAUDE.md`). The plan is `docs/architecture/RESTRUCTURE_PLAN.md`; read it first. **Update your own section before you finish**; do not edit another session's section.

Repository: `nyx2` (jj, colocated git, full history of `nyx`, **no remote, never push**; ADR 27). The old `nyx` repository is untouched.
Last updated: 2026-11-04, by the foundation session.

## Phases

| Phase | Name | State | Bookmark | Notes |
|---|---|---|---|---|
| 0 | Baselines | **done** | `restructure-0-baselines` | tags `pre-restructure`, `sbs-last`; numbers below; `tools/baseline.sh` |
| 1 | Delete sbs, chips | **done** | `restructure-1-delete-sbs-chips` | G1 answered (keep axisymmetric, check them; `sms_kz` kept); ADR 24 |
| F | Foundation (workspace layout, crate skeletons, `Physics` trait, conventions) | **done** | `restructure-foundation` (tip) | ADR 25, 26, 27; start every parallel workspace from this bookmark |
| 2 | femcore: one mesh model, all meshing | ready | `restructure-2-mesh` | see "Parallel work" |
| 3 | femcore: one linear algebra | ready (gate G2 before it ends) | `restructure-3-linalg` | |
| 3x | G2 experiment (read-only) | ready | `restructure-3x-g2-experiment` | |
| 4 | nyx-physics content, kernels into the physics crates | waits for 2 and 3 | | the interface exists (F); gate G3 answered |
| 5 | nyx-analysis | waits for 4 | | |
| 6 | Materials | **implemented** (+ extends, more models, provenance, checked library), awaiting paul's answers (defaults in its section) | `restructure-6-materials` | ADR 29, 30, 34, 35; `docs/MATERIALS.md` |
| 7 | nyx-engine, MCP parity, thin app | waits for 4, 5, 6 | | gate G5; 7a (connect merge) is ready now |
| 8 | File format | waits for 7 | | gate G6 |
| 9 | Close-out | waits for 8 | | |

Bookmarks stack: `main` (= `arch-model-adrs` tip f74c12ea) <- `restructure-0-baselines` <- `restructure-1-delete-sbs-chips` <- `restructure-foundation`.

## How to start a working session

```
cd C:\Users\robin-adm\Documents\opomech\nyx2
jj workspace add ../nyx2-<name> -r restructure-foundation     # <name> from the table below
cd ../nyx2-<name>
jj bookmark create restructure-<n>-<name> -r @                # or: jj new first, then create the bookmark on your own commit
```

* All workspaces share the build directory `..\nyx2-target` (`.cargo/config.toml`). The first build waits for the lock if another session is compiling; registry dependencies are not recompiled. **Disk**: C: had 17 GB free; the shared target is about 3 GB. Do not create extra target directories.
* Run long builds and tests in the background; the machine is shared by all sessions, so timings taken while another session compiles are worthless (see "Benchmark noise").
* Never push, never touch `main`, never use the `default` workspace of another session.

## Interfaces fixed by the foundation (parallel sessions may rely on them)

* **`nyx-physics`** (ADR 25): `Physics` trait (`id name studies bc_kinds properties params fields exports imports validate sizing solve`), `Ctx`, `LabelResolver`, `ProblemSpec` envelope, descriptors, `SolveResult`, `Error` (with `Cancelled`, `NotYetMoved`), `Registry` (`solve_order` of one-way chains), `validate_spec`, `schema::physics_schema`, `contract::{descriptor_snapshot, check_exemplar}`, `canonical_json`.
* **`nyx_physics::Mesh`** is `femcore::mesh::Mesh`, a type alias of `SimplexMesh`. Phase 2 changes what it is; code outside femcore uses only the alias.
* **Layering test** `crates/nyx-engine/tests/layering.rs` holds the allowed dependencies and the temporary edges; a session that adds a crate or a dependency changes the table in the same commit.
* **Contracts** (ADR 19): `problem_spec.v2.json`, `solve_result.v2.json`, `physics_descriptors.v1.json`, `physics_{dfs,sms,ems}.v1.json` are new; `problem_spec.v1.json` (old `ProblemSpec`) stays until phase 4. Exemplars are compared as values.
* **Test rule** (ADR 26): never compare JSON text or key order.
* **Descriptors of dfs/sms/ems are a first version** read from the current kernels. The phase 4 owner refines them.

## Parallel work ready now

Each row is one working session in its own jj workspace of `nyx2`, started from `restructure-foundation`.

| Session | Plan | Crate / files | jj workspace | Bookmark | Gate questions for paul |
|---|---|---|---|---|---|
| mesh | phase 2 | `femcore`: `mesh view geomdrive mesher mesher3d refine3d meshing meshing_pml pml3 physics_size estimate estimate3 cellmesh meshers/**` and their tests; contract `mesh` | `nyx2-mesh` | `restructure-2-mesh` | R05: keep `mesher3d`/`refine3d` as test fixtures only? (plan default yes). Confirm "domains and boundaries numbered from 1 everywhere". Where does the neutral `Sizing` live (proposal: `femcore::meshing::Sizing`, and `nyx_physics::Sizing` becomes a `pub use` of it)? |
| linalg | phase 3 | `femcore`: kernels `dfs sms sms_* ems ems_*` solve stages, old `sparse csparse linalg amg eigen lobpcg`, `solvers/**` | `nyx2-linalg` | `restructure-3-linalg` | **G2** (tolerance and reference platform; the experiment below informs it). Is "benchmarks within 10%" best-of-3 on an idle machine? Check the axisymmetric kernels against closed forms before moving them (G1 note) |
| g2 | phase 3 (read-only) | no source change: runs the two failing `ems` goldens with `portable-kernels` and with a tighter eigensolver tolerance, writes the finding in this file | `nyx2-g2` | `restructure-3x-g2-experiment` | none to start; at the end: which rule (plan section 7 options A, B, C) |
| materials | phase 6 | `femcore::material` (+ `toml` dependency, a `materials/` library outside the code), `physics_size` untouched | `nyx2-materials` | `restructure-6-materials` | G4 answered (TOML). Shipped library location (`materials/*.toml` at the repository root?). Import of the refractiveindex.info format: yes or no? Which dispersion models first (Sellmeier, Cauchy, tabulated are the plan's minimum)? |
| connect | phase 7a | merge `crates/app/src/connect.rs` and `crates/nyx-connect` into one `nyx-connect`; `app` depends on it | `nyx2-connect` | `restructure-7a-connect` | none |
| dsl | phase 9 (early, optional) | `docs/architecture/workspace.dsl`, `diagrams.md`: add the new crates, validate with Structurizr | `nyx2-dsl` | `restructure-9a-structurizr` | May the session run Docker or download the Structurizr CLI? (blocked by permissions last time) |

Waits for the above: phase 4 (needs 2 and 3), then 5 (needs 4), 7 (needs 4, 5, 6), 8, 9. **Phase 4 can fan out after its first step**: one session lands `dfs` (and any change to `nyx-physics` it needs), then `sms` and `ems` are separate files and can be moved by two sessions at once.

Not in the plan, needs paul's go: **engine core early.** `Project`, the command set and the job runner of `nyx-engine` do not depend on the kernels (solve is stubbed by `NotYetMoved`). Starting them before phase 4 would overlap with phase 2 only in `app/src/state.rs` and the meshing calls. The plan puts phase 7 after 4; say if you want it earlier.

## Baselines

**Machine**: Windows 11, 20 logical processors, 31.7 GB, Rust 1.99.1 pinned. Produced with `tools/baseline.sh`.

### Test counts

| State | passed | failed | ignored | Notes |
|---|---|---|---|---|
| `pre-restructure` (f74c12ea, code of main a047506), `cargo test --release --workspace --no-fail-fast`, run 2026-11-03 | 694 | **24** | 17 | 7 are the known golden failures; **17 more are an artefact**: with `--workspace` the app switches on `serde_json/preserve_order`, and the golden comparison compared key order. `-p femcore` alone gave the 7 (CLAUDE.md of that date said 632 / 7 for a different selection). The `ai/` workspace was separate (its tests passed) |
| phase 1 + foundation (`restructure-foundation`), same command | 764 | 2 | 9 | the only failures are `golden_soi_wire_modes` and `golden_soi_wire_with_sloped_wall_hybrid_modes` (R22/G2). The root workspace now includes the MCP crates, so these counts cover them. Fixed on the way: goldens compare sorted key lists, the MCP tool-list contract compares values |

Failing tests before phase 1 (the 7): `golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`, `golden_brillouin_wire_in_a_cubic_crystal`, `golden_brillouin_wire_in_a_tilted_crystal`, `golden_brillouin_forward_scattering_in_a_tilted_crystal`, `golden_ring_brillouin_cubic_crystal`, `golden_ring_brillouin_tilted_crystal`. Five were deleted with Brillouin.

### Benchmark timings (wall seconds, `NYX_PROFILE=1`, quiet machine)

Arguments are fixed in `tools/baseline.sh`. `bench` lists seven solves; the 3D eigenfrequency one dominates.

| Case | pre (run 2) | phase 1 (run 2) | Note |
|---|---|---|---|
| `bench` total | 46.4 | 46.5 | |
| electrostatics 2D, 38,697 dofs | 0.09 | 0.09 | |
| solid static 2D, 77,394 dofs | 0.27 | 0.27 | |
| solid 6 eigenfrequencies 2D | 1.78 | 1.72 | |
| solid static 3D, 205,479 dofs (AMG-CG) | 2.17 | 2.17 | |
| solid 6 eigenfrequencies 3D | 33.7 | 33.8 | |
| solid static 3D direct | 4.84 | 4.88 | |
| solid static 3D AMG-CG | 2.16 | 2.16 | |
| EM modes 2D (Nedelec 2), 10,059 dofs | 0.89 | 0.89 | |
| `profile3d 16 16 16 p2` | 2.18 | 2.22 | includes process start |
| `mem_profile 12 12 12 p2 static auto` | 1.28 | 1.30 | |
| `cavity_bench 0.012 2` | 4.34 | 4.43 | |
| `solve_scaling` | 14.7 | 13.1 | block width 8: 14.9 ms per vector (0.32x of one) in both |

**Benchmark noise.** The first runs, taken while other work compiled or ran on the machine, gave 63 and 83 s for `bench` and 49.3 / 59.7 s for the 3D eigensolve (against 33.7 once the machine was idle). The 2D cases are stable to a few percent, the 3D eigensolve moves by up to 40 % with machine load. A phase compares with the **best of three runs on an idle machine** (`BENCH_ONLY=1 bash tools/baseline.sh <dir>`); the plan's "within 10 %" is only meaningful that way. Raw logs of the five runs are in `C:\Users\robin-adm\Documents\opomech\nyx2-baseline-out` (outside the repository).

### Contract suites

All pass on `restructure-foundation`: `cargo test --release --workspace --test contracts` and `-p nyx-engine --test layering`.

## Foundation: what was done and decided (this session)

1. **Repository.** `nyx2` is cloned from `nyx` with the **full history** of the `arch-model-adrs` tip, not a clean import: the goldens, the provenance of the ports, the commit hashes cited by the ADRs and `blame` stay valid. Remote removed. Tags `pre-restructure`, `sbs-last`.
2. **Phase 0**: `tools/baseline.sh`; numbers above.
3. **Phase 1**: deleted `sbs.rs`, `sbs_ring.rs`, `chips.rs`, `Physics::Brillouin`, `BrillouinSpec`, `PhotoSpec`, `tests/brillouin.rs`, the five Brillouin goldens, the chip and rib examples, the app's Brillouin tab, example and template. Kept what other code needed from `sbs.rs` as the new `femcore::optical` (optical mode set, field sampler, `EPS0`/`MU0`, `fine_rule`); its field evaluator now computes only E (D and H were used by the Brillouin coupling only). Kept `usecases/`, `sms_kz`, `ems_axisym`, `sms_axisym` (G1 answered: axisymmetric ones to be checked). The app's examples 15 to 18 became 14 to 17 and the "Brillouin gain chip" template became "Strip on sapphire". `problem_spec.v1.json` rewritten on purpose (1.1). `geom` keeps its chip-shape goldens (geometry tests). ADR 24.
4. **Layout**: `ai/` crates moved to `crates/` and into the root workspace (the plan had this in phase 7); `nyx-connect` is still duplicated in `app/src/connect.rs` (merge = session `connect`). `.cargo/config.toml` shares the target directory. ADR 26.
5. **New crates**: `nyx-physics` (the interface), `nyx-dfs`, `nyx-sms`, `nyx-ems` (descriptors, `solve` = `NotYetMoved`), `nyx-analysis` (the `Analysis` trait only), `nyx-engine` (registry from cargo features, layering test). rustfmt-formatted (settings as `geom`) and clippy-clean. ADR 25.
6. **Conventions** in `CLAUDE.md`: test rule on JSON order, contract rules, workspace per session, shared target directory.

### Deviations from the plan (and why)

* The `Physics` trait has one more method, `params()` (descriptors of the envelope's physics-level options); the plan listed `params` in the envelope but had no descriptor for it.
* `ProblemSpec::order` is validated as 1 to 9 by the envelope; the physics decides what it supports (the mesher feature chain wants order above 2 later).
* `ai/` joined the root workspace now, not in phase 7; reason: one lockfile and one `--workspace` run make the baselines meaningful and let the layering test see everything.
* `Sizing` is defined in `nyx-physics` for now; the mesh session decides where the neutral type lives (see its gate question).

### Findings the next sessions should know

* **Hidden coupling by cargo features.** Any test that compares JSON key order or text differs between `-p crate` and `--workspace` builds (`preserve_order`). Fixed in the goldens and the MCP contract; do not reintroduce.
* **Benchmark noise** as above; the machine is shared.
* **Disk**: 17 GB free on C:. `nyx-jj/target` is 29 GB and is not used by `nyx2`; paul may delete it to make room.
* The CLAUDE.md of the old repository counted 632 passed / 7 failed; with the whole workspace the counts are the ones in the table.

## Mesher features (`FEATURE_DEPENDENCIES.md`) against the plan

That file lives uncommitted in the old checkout (`nyx/docs/ports/meshers/FEATURE_DEPENDENCIES.md`, R24); copy it into `nyx2` when the mesher work starts. They come after the restructure and only where they fit. Tweaks the plan needs, for paul's information:

* **Wait for phase 2's file moves.** H2, H5, T1 to T4, C0 to C3, E1, E3, P3 edit `ops3.rs`, `geom3.rs`, `estimate3.rs`, `physics_size.rs`, `pml3.rs`, `cellmesh.rs`, `mesh.rs`, all of which phase 2 moves or rewrites. jj does not follow renames: a concurrent edit to the old path becomes a modify/delete conflict. **The mesh session should land a pure "move the files" commit first** (no logic changes, a bookmark of its own), so mesher streams can start from it. Items that touch only tests and docs (H1, H3, V1, T5, H4's fixtures) can start earlier.
* **X1 (native mixed-cell interface) is phase 2's main product** and P4/P5 (adapter from `PmlBox::cells` to the solvers) are then phase 3/4 work: the kernels read the native mesh through the split until they take cells directly. The mesh contract must state the node order per cell type and order.
* **X2/C4 (stored curved nodes)**: ADR 18 defers them. Phase 2's mesh model should reserve `nodes_per_cell(type, order)` and keep `mesh` contract v1 "no curved-node storage" until C4; no restructure change needed.
* **X3 (Nedelec above 2)**: the envelope no longer hard-codes order 1 or 2; `ems` descriptors and `Physics::validate` carry the limit.
* **X4/X5/X6** (driven study, time domain, PML along general normals) are new study kinds or `pml` condition parameters of the physics crates: they fit the trait; no plan change.
* **E7 (Brillouin and ring studies with sized meshes)**: Brillouin is gone; the ring (axisymmetric) part remains with `ems_axisym`. **H4 (four-chip baselines)**: `chips.rs` is gone (tag `sbs-last`), the geometries survive as `crates/geom/golden` chip shapes; rebuild the fixtures from those with generic ems/sms studies.
* **A1 to A3 (app routes)** become engine commands in phase 7, not `app/src/state.rs` changes. **C0/E5/E6 (`physics_size.rs`)** target the neutral `Sizing` and the mesher-side `PhysicsSizing` (points per wavelength), not the physics crates: calibration stays on the mesher side, the physics only supplies wavelengths.
* **H6 (Mmg/Omega_h comparison)** and D1/D2 stay paul's decisions.

## Gate answers from paul (2026-11-03 evening, in the thread)

* **G1**: keep the axisymmetric kernels (`ems_axisym`, `sms_axisym`) **with a note that their implementation must be checked**, paul is not sure of it. Action for the phase 3 and 4 owners: before moving them, review them against closed forms (the earlier checks quoted in the pre-restructure README: whispering-gallery modes, breathing and flexural modes of a ring) and say in this file what was verified. About `sms_kz`, paul does not know why it exists after the axis changes: it is the solid-mechanics kernel for a 2D cross-section with an out-of-plane wave number k (three displacement components, `u(x, z) exp(-j k y)`, wafer axes since the kernel port), i.e. the phonon waveguide dispersion, not Brillouin-specific. It stays (nothing was deleted); paul may still ask for its removal.
* **G3**: "sure": crate names as in the plan (`nyx-physics`, `nyx-dfs`, `nyx-sms`, `nyx-ems`, `nyx-analysis`, `nyx-engine`).
* **G4**: "sure": TOML for materials.

## Open questions for paul

Gates of the plan not yet answered: **G2** (tolerance and reference platform; needed before phase 3 ends), **G5**, **G6**. New: the session list above (engine core early? Docker for Structurizr?), and whether to free `nyx-jj/target`.

## Log

* 2026-11-03: architecture model, 23 ADRs, 28 rough spots, contract tests (bookmark `arch-model-adrs`), restructure plan written.
* 2026-11-04: foundation session: `nyx2` created, phases 0 and 1, workspace layout, `nyx-physics` and the skeletons, ADRs 24 to 27, parallel list.

## Materials session (phase 6), written by the materials session, 2026-11-04

State: **implemented, awaiting paul's answers to three questions** (defaults below). Bookmark `restructure-6-materials` (workspace `nyx2-materials`), cut from `restructure-foundation`. Nothing pushed.

Done:
* `femcore::material` is now a folder: `mod.rs` (the old API, unchanged except the new variant `Stiffness::General(M6)`), `units.rs`, `reader.rs` (TOML sections, unknown keys with suggestions, plain errors), `elastic.rs`, `tensor.rs`, `optics.rs`, `spec.rs` (`MaterialSpec`, `Library`), `schema.rs`. `toml` (feature `preserve_order`) is a `femcore` dependency.
* Elastic: class and constants (isotropic, cubic, hexagonal, trigonal, tetragonal, orthotropic, monoclinic y or z, triclinic), reduced matrix (upper triangle or full), compliance (constants or matrix; the shear factors are derived), E/nu, lambda/mu, K/G; symmetry, stability (Mandel eigenvalues, Born criteria for cubic), class consistency and derived-entry checks. Tensors: 1, 3, 6 entries, 3x3, `{perp, par}`, with positivity checks.
* Optics: constant, Sellmeier, Cauchy, tabulated n,k (PCHIP or linear), per axis (ordinary/extraordinary or x/y/z), declared wavelength unit, required validity range for formula models, errors outside the range, `dn_dT`, complex index `n - jk`, `d(n)/d(lambda)`, group index, permittivity tensor rotated into the wafer frame.
* Library `materials/{silicon,silica,sapphire,air}.toml` (loaded at run time; `Library::load_dir`). Silica and sapphire use the Malitson coefficients quoted from the refractiveindex.info database; silicon the Li (293 K) table. Tests pin the library to the in-code `sapphire()` and `silicon()`.
* ADR 29 and 30 (not 28: the connect session took it), `docs/MATERIALS.md` (user guide; its TOML examples are tested), schema `material_schema.v1.json`, `material_kitchen_sink.v1.{toml,json}`, `VERSIONS.md` and `contracts/README.md` rows, `CODE_TO_ADR.md` rows (`check_code_to_adr.sh` passes).
* Tests: 20 unit tests in the module, `crates/femcore/tests/materials.rs` (28: published values, all stiffness forms, round trips, every error message), `tests/material_contract.rs` (3). Layering test unchanged and passing.

Defaults I used for the questions to paul (each is a one-line change):
1. **Library location**: `materials/*.toml` at the repository root, loaded at run time (`Library::load_dir(path)`); the code holds no copy. Change: move the folder and the path in `tests/materials.rs`.
2. **refractiveindex.info import**: not implemented; `source` records provenance. The sellmeier model accepts `lambda0` (the resonance wavelengths, as that database lists them), so coefficients paste without conversion.
3. **Dispersion models**: constant, Sellmeier, Cauchy, tabulated first. Lorentz, Drude, Schott, Herzberger are additions to the `Dispersion` enum.

Left for others / later:
* `sapphire()` and `silicon()` are still in `material/mod.rs` and used by the app and tests; the plan says remove them. Removing needs the app to load the library (engine, phase 7). The app can already use `MaterialSpec::anisotropy(euler)` and `legacy_scalars()`.
* Wiring evaluation into a study (`ems` property descriptor `optical.eps`, a wavelength sweep in the engine) is phase 4 and 7.
* Numbers marked "approximate" in `materials/silicon.toml` (dn_dT) and the mechanical, thermal values of silica and sapphire are typical handbook values; check before quantitative use.
* `femcore --test contracts` and the other sessions' files were not touched.

Test run of this session (`cargo test --release --workspace --no-fail-fast`, Windows, other sessions compiling at the same time): 809 passed, 9 ignored, 3 failed: the two known `golden_frames` failures (R22) and `nyx-store workspace::tests::dashboard_rules_hide_after_three_days_cap_pin_and_never_delete`, which passes when run alone (an unrelated, probably order- or id-dependent flake in a crate this phase did not touch; noted for whoever owns `nyx-store`).

### Materials follow-up (same day): `extends`, more models, provenance, checked library

* **`extends`** (ADR 34): a file may start from a library material (`extends = "silicon"`); top-level keys and whole tables (`[elastic]`, `[thermal]`, `[electrical]`, `[optical]`) of the file replace the base's, never mixed key by key; name not inherited; checks run on the merged result; missing base and loops give plain errors (chain shown). `Library::load_dir` resolves chains, `Library::read(text, name)` resolves one text. Example `materials/examples/silicon_measured.toml`, tests `tests/material_extends.rs`, guide section tested. Note: the coordinator first asked for field-level overrides (one axis, one dn/dT), then for the whole-table semantics proposed to paul; I implemented whole-table replacement (restate the table to change one axis).
* **Models** (ADR 35): `schott`, `herzberger`, `drude`, `lorentz`, `drude_lorentz` (frequency units eV, meV, THz, GHz, Hz, rad/s, 1/cm; e^{j omega t}); same rules (declared units, required range, error outside). Drude-Lorentz reproduces the Rakic tables of refractiveindex.info for Au, Ag, Cu, Al to 2e-4 (rows compared; 5-digit tables), which also checks the sign convention. Unit invariance (Herzberger um vs nm, eV vs THz vs rad/s) and finite-difference derivatives are tested.
* **Provenance** (ADR 35): `source` as text, `{cite, url, doi, verified}` or a list, at the top and in every table and optical axis; kept by the loader and written back; a test requires a source on every file and table of the library and the words approximate or unverified on any entry without a `verified` date. Schema is now 1.1 (additive; `extends`, sources, new model keys); VERSIONS rows added; the first kitchen-sink exemplar was rewritten on purpose (sources) and a second one (`material_kitchen_sink_models`) added.
* **Library checked against the literature (2026-11-04)**. Corrected: silicon stiffness 165.7/63.9/79.6 to 165.6/63.9/79.5 GPa (Hopcroft, Nix, Kenny 2010); silicon thermal conductivity 130 to 148 W/(m K) (approximate: sources 130 to 156); silicon dn/dT 1.8e-4 to 1.87e-4 1/K (Frey et al. 2006 via a secondary compilation; Komma 2012 at 1550 nm); silica E 73.1 to 72.7 GPa, nu 0.17 to 0.16, density 2203 to 2201, specific heat 703 to 770, expansion 0.55e-6 to 0.52e-6 (Corning 7980 data), plus dn/dT 9.7e-6 (Corning, 1530 nm); sapphire parallel conductivity 25.2 to 25.8, parallel expansion 5.3e-6 to 5.6e-6, parallel permittivity 11.6 to 11.5, specific heat 760 to 756 (supplier data sheets; they disagree, so flagged approximate). Confirmed: sapphire stiffness (Wachtman, Tefft, Lam 1960: c11 4.968, c12 1.636, c13 1.109, c14 -0.235, c33 4.981, c44 1.474 x 1e12 dyn/cm^2), the silica, sapphire (o and e) Sellmeier coefficients and the silicon table (read from the database files with their references). Added: gold, silver, copper, aluminium (Rakic LD), silicon_nitride (Luke 2015), lithium_niobate (Zelmon 1997); optical data only, as their source tables say.
* **Flagged as approximate or unverified in the files**: silicon thermal conductivity, specific heat and permittivity; silica permittivity (3.6 to 3.9 range, midpoint used); sapphire thermal, expansion and permittivity; air (all values; air is treated as vacuum, n = 1). Silicon nitride is a film (the index depends on deposition). dn/dT of sapphire and lithium niobate is not in the library.
* ioffe.ru (the usual source for semiconductor data) has an invalid TLS certificate; I did not bypass the check (the tool refused it) and used other sources for silicon. If paul wants Ioffe values, they can be added by hand.
* The legacy `silicon()` and `sapphire()` in `material/mod.rs` keep the earlier rounded numbers (tests of other sessions and goldens use them); the pin test in `tests/materials.rs` allows the small differences and documents them. The ADR numbers 34 and 35 were chosen after 31 to 33 were taken by the mesh and linalg sessions; renumber if another session took them in the meantime.
