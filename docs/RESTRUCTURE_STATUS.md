# Restructure: status

Handoff file between working sessions (see `CLAUDE.md`). The plan is `docs/architecture/RESTRUCTURE_PLAN.md`; read it first. **Update your own section before you finish**; do not edit another session's section.

Repository: `nyx2` (jj, colocated git, full history of `nyx`, **no remote, never push**; ADR 27). The old `nyx` repository is untouched.
Last updated: 2026-11-04, by the foundation session.

## Phases

| Phase | Name | State | Bookmark | Notes |
|---|---|---|---|---|
| 0 | Baselines | **done** | `restructure-0-baselines` | tags `pre-restructure`, `sbs-last`; numbers below; `tools/baseline.sh` |
| 1 | Delete sbs, chips | **done** | `restructure-1-delete-sbs-chips` | G1 answered (keep axisymmetric, check them; `sms_kz` kept); ADR 24 |
| 1b | Remove `sms_kz` (out-of-plane model), on paul's instruction | **done** | `restructure-1b-remove-sms-kz` (on top of `restructure-foundation`) | 755 pass, 9 ignored, 2 expected golden failures; merge it before phase 3 touches `sms_kz.rs` |
| F | Foundation (workspace layout, crate skeletons, `Physics` trait, conventions) | **done** | `restructure-foundation` (tip) | ADR 25, 26, 27; start every parallel workspace from this bookmark |
| 2 | femcore: one mesh model, all meshing | ready | `restructure-2-mesh` | see "Parallel work" |
| 3 | femcore: one linear algebra | ready (gate G2 before it ends) | `restructure-3-linalg` | |
| 3x | G2 experiment (read-only) | ready | `restructure-3x-g2-experiment` | |
| 4 | nyx-physics content, kernels into the physics crates | waits for 2 and 3 | | the interface exists (F); gate G3 answered |
| 5 | nyx-analysis | waits for 4 | | |
| 6 | Materials | ready (G4 answered: TOML) | `restructure-6-materials` | |
| 7 | nyx-engine, MCP parity, thin app | waits for 4, 5, 6 | | gate G5; 7a (connect merge) is ready now |
| 8 | File format | waits for 7 | | gate G6 |
| 9 | Close-out | waits for 8 | | |

Bookmarks stack: `main` (= `arch-model-adrs` tip f74c12ea) <- `restructure-0-baselines` <- `restructure-1-delete-sbs-chips` <- `restructure-foundation`.

Bookmarks that exist after the machine restart (2026-11-04, none is merged into another yet; whoever integrates decides the order): `restructure-1b-remove-sms-kz` and `restructure-1c-gate-answers` (children of the foundation; 1c holds every answered gate and should be the base of new sessions), `restructure-2-mesh`, `restructure-2b-mesh3d-adapter`, `restructure-3-linalg`, `restructure-3x-g2-experiment`, `restructure-4e-engine-core` (engine core early, approved), `restructure-6-materials`, `restructure-7a-connect`. Not created yet: `restructure-9a-structurizr`. Each session records its own state in its own section below.

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

**Engine core early: approved (2026-11-04), thread started from `restructure-1c-gate-answers`.** `Project`, the command set and the job runner of `nyx-engine` do not depend on the kernels (solve is stubbed by `NotYetMoved`). Starting them before phase 4 would overlap with phase 2 only in `app/src/state.rs` and the meshing calls. The plan puts phase 7 after 4; say if you want it earlier.

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

* **G1**: keep the axisymmetric kernels (`ems_axisym`, `sms_axisym`) **with a note that their implementation must be checked**, paul is not sure of it. Action for the phase 3 and 4 owners: before moving them, review them against closed forms (the earlier checks quoted in the pre-restructure README: whispering-gallery modes, breathing and flexural modes of a ring) and say in this file what was verified. About `sms_kz`, paul does not know why it exists after the axis changes: it is the solid-mechanics kernel for a 2D cross-section with an out-of-plane wave number k (three displacement components, `u(x, z) exp(-j k y)`, wafer axes since the kernel port), i.e. the phonon waveguide dispersion, not Brillouin-specific. **Update 2026-11-04: paul said to remove it** (the axis is along y now, so the code is obsolete): done on bookmark `restructure-1b-remove-sms-kz` (ADR 24 amendment). What went: `ProblemSpec::out_of_plane`, `solve_solid_kz*`, `assemble_stiffness_kz`, the 4 `solid_kz_*` goldens, `tests/solid_kz.rs` (its Floquet-vs-cross-section test too), the app checkbox, the example's `out_of_plane`. What stays: the shared cross-section machinery in `sms_kz.rs`, because `sms_axisym` is built on it; phase 4 renames and moves it into `nyx-sms`.
* **G2** (2026-11-04): loose tolerances everywhere, no reference platform (option B); recorded by the g2 session.
* **G5** (2026-11-04): the MCP may manipulate **every tool geom offers** (the whole feature-list editing surface), not a minimal set. The plan's "first set (strip, rib, ring)" is only a build order and the acceptance test: those three must be rebuildable through MCP; the tool set is generated from geom's feature descriptors so nothing is left out. Phase 7 owner: derive the tools from geom, do not hand-pick.
* **G6** (2026-11-04): the JSON layout of the project file is **not fixed**; it is derived from the `ProblemSpec` envelope, `Project` and `FeatureList` in phase 8 (the plan's default stands). No layout decision is needed before then.
* **Structurizr** (2026-11-04): paul allows downloading the Structurizr CLI (session `dsl`, bookmark `restructure-9a-structurizr`, optional and low priority).
* **G3**: "sure": crate names as in the plan (`nyx-physics`, `nyx-dfs`, `nyx-sms`, `nyx-ems`, `nyx-analysis`, `nyx-engine`).
* **G4**: "sure": TOML for materials.

## Engine core early (session `nyx2-engine`, bookmark `restructure-4e-engine-core`, 2026-11-04)

Started from `restructure-1c-gate-answers`. ADRs 37, 38, 39. **State: done for the part that does not depend on phases 2 to 4; the rest waits (below).**

**Dependency check (asked by paul: "start the engine core if it has no dependencies").** From the plan and the code: `Project`, the commands, the job runner and the store wiring do not depend on the kernels or on the linear algebra (phase 3). The old `app::Project` could not simply move: it is built on the old `femcore::problem::ProblemSpec` (replaced by the envelope in phase 4) and calls `femcore::{meshing, mesher, geomdrive, physics_size}` (phase 2 replaces and moves those). The engine therefore uses the new envelope and `geom` directly, and has exactly two couplings: **the mesher** (one function, `FemcoreMesher::mesh` in `meshing.rs`; phase 2 rewrites it) and **`Physics::solve`** (answers `NotYetMoved` for the skeleton physics until phase 4). No file owned by the phase 2, 3 or 6 sessions was touched.

**Built.** `Project` (feature list, mesh settings, one envelope per physics, active), JSON with sorted keys; `Engine` with caches keyed by content (stale results are reported, nothing is invalidated by hand); 64 commands (`Engine::catalog()` / `Engine::call(name, json)`), each with a JSON Schema: project, geometry edits, physics/study/material/condition edits, mesh, solve, jobs, results, store; a job runner (thread per job, `Progress`, cancel, `info`, `try_take`, `wait`, panic becomes a failed job); the coupled solve chain (order from `Registry::solve_order`, imports filled from exports); store wiring (create or open a workflow folder, commit versions, every solve, failed or cancelled run is a row of the run table with metrics, log and result files); the geometry descriptors (G5): one `add_*` command per geom operation (23), generated from `nyx_engine::geomdesc`, with a test that the descriptors equal the serde fields and that the strip, rib and ring can be built, evaluated and meshed through commands alone. `nyx-store`: specs may be JSON, and the caller supplies the card function (ADR 39); additive, the app is unchanged.

**Tests.** `cargo test --release -p nyx-engine`: contracts 2, engine 10 (strip, rib and ring through commands, stand-in physics for solve/chain/cancel/store/reopen, refusals leave the project unchanged), geomdesc 7, layering 2, registry 4. Whole workspace (`cargo test --release --workspace --no-fail-fast`): 775 pass, 9 ignored, 2 fail, the two expected golden failures (R22); the base `restructure-1c-gate-answers` has 755 pass, so the 20 new tests are the difference. clippy `-D warnings` and rustfmt clean for `nyx-engine`. `nyx-store`: one new test. New contract exemplars: `engine_commands.v1.json`, `geom_operations.v1.json`, `engine_project.v1.json` (provisional, version 0.1; see `VERSIONS.md`).

**What waits (and on what).**
* Solving with the real physics: phase 4 (`NotYetMoved` is reported and recorded as a failed run). Nothing in the engine changes then.
* `physics_sizing` of the mesh: phase 2 (neutral `Sizing`); the engine answers `Unsupported`. 3D models (feature list `dim` 3) do not mesh yet (also the case in the old app); a 2D model can be swept along y with `sweep` layers.
* Rewiring `app` and `nyx-mcp` to `Engine::call` and dissolving `nyx-session`: phase 7 proper. The MCP server can become a loop over `catalog()` and `call()`; it has to add the handshake, the `nyx.project_hash` header and the assumptions ledger of the old session if wanted.
* File names in the store (`experiment.nyx.xml` holds JSON for an engine project), removal of the XML branch and of `femcore` from `nyx-store`: phase 7/8. The final project JSON layout and its `format` key: phase 8 (G6).
* Not done: `set_material` takes a material name or scalar properties (the TOML files and dispersion models are phase 6, `femcore::material`); the engine does not read material files yet. Analyses (`nyx-analysis`) are not wired: phase 5. The batch runner and HDF5 export of the plan's engine list are not built (HDF5 stays an optional export; batch is `Engine::call` in a loop).

**For the other sessions.** Phase 2 owner: change `FemcoreMesher::mesh` only; names and labels come from `geom::labels::LabelMap`, which already numbers domains from 1. Phase 4 owners: your physics are exercised by `solve` through `Ctx` (resolver, progress, imported fields by role); `describe_physics` and the schema of each command come from your descriptors, so refining a descriptor changes `engine_commands.v1.json` (rewrite it with a `VERSIONS.md` line). Phase 7 owner: start from `Engine`; do not write a second project type.

**Findings.** `geom` depends on serde only and a test pins that (ADR 3), so the descriptors live in the engine, not in `geom`. A running test executable in another session locked the shared target directory once (LNK1104); the coordinator suggested a private `CARGO_TARGET_DIR`, which was not needed here.

## Open questions for paul

Nothing is open for paul from the foundation any more: the plan's gates G1 to G6 are answered, and the coordinator reports the rest settled on 2026-11-04: mesh defaults accepted (`mesher3d`/`refine3d` kept as a fixture with a README, numbering from 1, `Sizing` in `femcore::meshing`), materials settled in its own thread (further models added, sources verified, unverified values flagged), mesher decision D1 = skip Mmg, D2 answered (the geom 3D finding went to paul), and **engine core early = yes** (a thread "Build engine core early" started from `restructure-1c-gate-answers`).

## Log

* 2026-11-03: architecture model, 23 ADRs, 28 rough spots, contract tests (bookmark `arch-model-adrs`), restructure plan written.
* 2026-11-04: foundation session: `nyx2` created, phases 0 and 1, workspace layout, `nyx-physics` and the skeletons, ADRs 24 to 27, parallel list.

## Phase 2 (mesh session, workspace `nyx2-mesh`, bookmark `restructure-2-mesh`)

State 2026-11-04: **step 1 to 3 done, the rest not started.**

Done (jj commits on the bookmark, in order):
1. `a4f8aa69` pure file move into `femcore::{mesh, meshing, adapt}` (ADR 31 lists the new paths); only `mod` lines, one import and a lint allow changed; old paths kept as a temporary re-export block.
2. `1fe9e8aa` `CODE_TO_ADR.md` rows for the moved files (`check_code_to_adr.sh` passes).
3. `5f36745b` all callers (femcore, tests, examples, app, nyx-session) use the new paths; the re-export block is gone; `lib.rs` module map updated.
4. `cdaa9b43` ADR 31.

Checks at `cdaa9b43`: `cargo build --release --workspace --all-targets` ok; `cargo test --release --workspace --no-fail-fast`: **764 passed, 2 failed, 9 ignored**, identical to the foundation baseline (the failures are the two known `golden_frames` `ems` tests, R22/G2). No logic changed; no golden, contract exemplar or benchmark touched (benchmarks not re-run: a pure move cannot change them, they are measured at the end of the phase).

Not done yet (in the plan's order): one native mesh (`mesh::native`/`cells` become the only mesh, `SimplexMesh` a view), delete `mesh::cellmesh` and `mesh::interop`, numbering from 1 everywhere, one PML generator (`meshing::pml` and `meshing::pml_box` still both exist), neutral `Sizing` in `femcore::meshing` with `nyx_physics::Sizing` re-exporting it, `adapt::par` onto `femcore::par`, `mesher3d`/`refine3d` as test fixtures only, rewritten `mesh_*` contracts and schemas (ADR 19), ADRs for those decisions (number 32 and up; coordinate with the linalg thread), best-of-3 benchmarks on an idle machine.

Gate questions for paul (open; defaults will be used and kept easy to reverse): R05 keep `mesher3d`/`refine3d` as test fixtures only (default yes); domains and boundaries numbered from 1 everywhere (default yes); neutral `Sizing` in `femcore::meshing::Sizing` re-exported by `nyx_physics` (default yes).

Notes for other sessions: any file that names the old paths (`femcore::geomdrive`, `femcore::meshers::*`, `crate::view`, `crate::estimate*`, and so on) must use `femcore::meshing::geomdrive`, `femcore::mesh::native` etc.; rebase onto `restructure-2-mesh` or apply the same path rewrite.

Update (mesh session, later the same day): commit `1de99a33` made `mesher3d` and `refine3d` a **test fixture** (`crates/femcore/tests/common/solid3d`, included by `tests/emw3d.rs` and `examples/cavity_bench.rs`); their own tests and `examples/mesh3d_bench.rs` were deleted. The fixture stays because the PEC sphere cavity (closed form, `emw3d`) and the `cavity_bench` baseline need a sphere mesh, which `extrude_y` cannot make. Rest of the plan (single native mesh, 1-based labels, one PML generator, `Sizing`, contracts, benchmarks) is written as ADR 32 (Proposed) with five commit-sized steps; **none of it is implemented yet**. Tests after the fixture move: `emw3d` 6/6; full suite not re-run since `cdaa9b43` (764/2/9).

### Mesher features against the plan (mesh session, 2026-11-04)

* **D1 (paul): the Mmg / Omega_h comparison (H6) is dropped.** Marked in `docs/ports/meshers/MISSING_FEATURES.md`.
* **D2: what 3D geometry geom has and what meshes it.**
  * `geom` describes 3D as a *prismatic layer stack* (`crates/geom/src/stack.rs`, `scene.rs`): `Box`, `Cylinder`, `Extrude`, `Cone`, `Revolve` along one shared sweep axis (`scene.rs:246-255`), `Transform` only about that axis (`scene.rs:345`), `Pattern`, `Union`/`Difference`/`Intersect`/`Overlay`; `Sphere` and `Torus` exist as analytic bodies but only inside such a host (`stack.rs:203`). It evaluates volumes, named faces and edges, and has analytic surfaces (`FaceGeom::Quadric`, `Complex::facet_project`, `complex.rs:409`; `Evaluated::classify(p)`, `eval.rs:149`). A rotation about another axis or two different sweep axes is refused.
  * **Nothing meshes a general 3D geom model.** `meshing::geomdrive::mesh_geom_2d` returns an error for `dim != 2` (`crates/femcore/src/meshing/geomdrive.rs:160`); `geom::tess` is display-only (`tess.rs:1`); `app/src/state.rs:16` says so: 3D problems are made by sweeping the 2D section along y (`MeshCache::mesh_extruded` -> `meshing::mesher::extrude_y`). So an extruded/prismatic geom solid is covered; a sphere, torus or revolved body is not.
  * **`mesher3d`/`refine3d` (fixture `crates/femcore/tests/common/solid3d`) are the closest existing piece, but not wired and not enough alone.** They stuff a BCC lattice from their own `Solid` enum (sphere, box, cylinder, extruded outline), not from a geom model. The missing adapter is small in concept: geom already offers `Evaluated::classify` (point to region) and `Complex::facet_project` (snap to an analytic face), which is what the stuffing needs in place of `Solid::sdf`/`face_at`. The adaptation side (`adapt::geom3::FacetGeometry`, `drive3`) can then improve the start mesh. Estimated work: a `geom -> Solid-like oracle` adapter, face/domain numbers from geom's `LabelMap`, tests on a sphere-in-box and a torus against closed forms. Not part of the restructure; decision for paul whether to promote the fixture back into `femcore::meshing` when that work starts (the fixture README says how).

### Phase 2 state at the end of the mesh session (2026-11-04)

ADR 32 steps 1 to 7 are done on `restructure-2-mesh` (tip = the commit that carries this text); one jj commit per step.

| Step | Commit | What |
|---|---|---|
| move | `a4f8aa69` (before the sms_kz rebase), callers `5f36745b`, ADR 31 | `femcore::{mesh, meshing, adapt}` |
| 1 | `01597c35` | boundary numbers 1-based everywhere; mesh contract 1.1 |
| 2 | `ae77dffe` | `SimplexMesh::from_cells*` take `mesh::cells::CellMesh<D>`; nyx `CellMesh`/`ElementType`, `cells_to_nyx`, `PmlBox.cells` deleted |
| 3, 4 | `f9f3af74` | `mesh::native` -> `adapt::work` (private working mesh of the adaptation), `mesh::interop` dissolved; ADR 32 amended; removal of `SimplexMesh` is a phase 4 deliverable (plan updated) |
| 5 | `db231977` | one PML generator (`pml_box::build_pml`), physical vertex and cell numbering kept and tested (answers the FEATURE_DEPENDENCIES P3 question: yes) |
| 6 | `9ff426a9` | neutral `femcore::meshing::Sizing`, re-exported as `nyx_physics::Sizing`; `physics_size` takes it, `ProblemSpec` only via `spec_sizing` |
| 7 | `c62349da` | `mesh_cells.v1.json` exemplar and contract test; VERSIONS and README rows |

* **Checks**: `cargo build --release --workspace --all-targets` OK; `cargo test --release --workspace --no-fail-fast`: **744 passed, 2 failed (the known `ems` goldens), 9 ignored** (the cell-model exemplar test, once it compares values, brings 743 to 744). The count moved from 764 by: -15 (tests of the deleted `tests/mesher3d.rs`), -9 (the `sms_kz` removal, `restructure-1b-remove-sms-kz`, rebased under this stack), +2 new tests (PML numbering, `Sizing`) and +1 (the cell-model exemplar). `check_code_to_adr.sh` passes; clippy `-D warnings` on the new crates is clean.
* **Benchmarks** (`BENCH_ONLY=1 bash tools/baseline.sh`, nothing else compiling, three runs of mine and two interleaved control runs of the parent commit `ef9106ef` on the same machine state). Mesh phase against its parent, best of the runs: `bench` total 82.9 s vs 85.4 s, 3D eigensolve 59.5 s vs 61.3 s, 3D direct solve 11.55 s vs 11.75 s, 2D eigensolve 1.71 s vs 1.70 s, `cavity_bench` 4.83 s vs 4.84 s, `solve_scaling` 18.3 s vs 18.2 s: **no difference beyond noise (within 3 %)**. **Caution**: on today's machine state both are about 1.8 times slower than the phase 0 numbers in the table above (3D eigensolve 59 to 70 s against 33.7 s), with no compile running; whatever changed (power plan, thermal state, background services after the restart) the old absolute numbers cannot be reproduced now, so phase 3 and later must compare against a control run of their parent commit, not against the phase 0 table.
* **Not done / for others**: `mesh::SimplexMesh` removal (phase 4); the `ProblemSpec` adapters in `physics_size` and `mesh_for_problem` (phase 4); a `geom -> mesher3d` adapter for general 3D shapes (separate thread, starts from the fixture `tests/common/solid3d`); `meshing::pml` was renamed `pml_geom` after the phase.
* **Practical findings**: `nyx2-mesh-target` (private build directory, about 3 GB) can be deleted when the phase is accepted; a shared target directory is locked by running test binaries (LNK1104) and by parallel compiles, so a private `CARGO_TARGET_DIR` per workspace is worth the disk. Line endings: several docs are CRLF in the repository (`RESTRUCTURE_PLAN.md`, `CODE_TO_ADR.md` before my edits); keep them when editing.

### Phase 2 after acceptance (paul: "do the renaming, then delete what needs to be")

* `meshing::pml` is now `meshing::pml_geom` (commit `f3deedcf`).
* The mesher's own `par` module is gone: `CHUNK`, `fold_chunks_mut`, `sum_range` live in `femcore::par`; `FEMMESH_THREADS` no longer exists, the adaptation uses `NYX_THREADS` (commit `20736971`, ADR 11 amendment).
* The `mesh_cells` contract test compared JSON text and failed under `--workspace` (`preserve_order`); it now compares values (`3b753030`). Suite after these commits: **744 passed, 2 failed (the known `ems` goldens), 9 ignored**; build, `check_code_to_adr.sh`, clippy on the new crates clean.
* Rough spots R01, R03 (mesher part), R04, R05, R09 marked resolved in `ROUGH_SPOTS.md`; `README.md` and `docs/HOW_IT_WORKS.md` describe the new layout.
* **Deleted**: the private build folder `nyx2-mesh-target` and the benchmark logs `nyx2-mesh-bench` (both outside the repository, mine); the temporary control workspace `nyx2-mesh-ctl` (already forgotten). **Audit of dead code**: no compiler warning in `femcore`; six public functions have no caller anywhere (`CellMesh::as_simplex`, `mesh::all_cell_dofs`, `fmsh::load_dataset`, `fmsh::save_steps`, `pml_box::layers_for_size`, `adapt::work::signed_measure`); they are file-format or conversion API and were kept (`as_simplex` is the bridge a later adaptation on `CellMesh` would use); `tests/common/solid3d` kept by decision.

## Phase 3 (linalg), session `linalg` (updated 2026-11-04, bookmark `restructure-3-linalg`, rebased onto `restructure-1b-remove-sms-kz`)

**State: the solve stages of every kernel run on the new linear algebra; the old modules are gone; open items below.** G2 was answered by paul (option B: loose tolerances everywhere, no reference platform; ADR 33). `EigenOptions::tol` stays configurable (`NYX_EIG_TOL` overrides it for experiments).

**Done**
* One matrix type `solvers::csr::Csr<T>` (+ `Triplets`/`Coo`/`CCoo`, aliases `RCsr`/`CCsr`/`C`/`c`), one door for the kernels `solvers::driver` (`DirectSolver<T>`, `eigs_near` and the wrappers `eigs_sym`, `eigs_sym_indefinite`, `eigs_complex(_projected)`, `eigs_general`, `lobpcg_guarded`, dense helpers, `EigenOptions`), progress and cancel in the driver. Elasticity multigrid moved to `solvers::amg_elastic` (hierarchy still on a private working matrix type, the old `sparse::Csr`/`Csr32`). The real non-symmetric 2D mode pencil is solved by `solvers::arnoldi_real` (the old real block Arnoldi, ported onto `Csr`).
* **Deleted**: `sparse.rs`, `csparse.rs`, `linalg.rs`, `amg.rs`, `eigen.rs`, `lobpcg.rs`, the bridge and the `NYX_OLD_LINALG` switch, examples `lu_scaling`, `pool_bench`. (Their unit tests went with them: the workspace count drops by about 9 tests.)
* ADR 33 (supersedes ADR 8), contract `linalg_options.v1.json` (+ `VERSIONS.md`), `tests/linalg_contract.rs` (8 tests incl. a converged complex pair), `tests/axisym_closed_form.rs`, `CODE_TO_ADR.md` rows, ADR 8 marked superseded.
* **Axisymmetric kernels verified against closed forms (G1)**: `ems_axisym` on a PEC coaxial cavity, TEM modes `f = p c / 2h`: relative error 1.4e-6 (p=1), 8e-6 (p=2); `sms_axisym` on a 1 m ring, 0.02 m square section: breathing mode 3.4e-5 from `sqrt(E/rho)/(2 pi R)`, in-plane flexure m=2: 1.9e-4, m=3: 5.0e-4 from Euler-Bernoulli ring theory (which has an error of that size itself). Both look right.
* **Complex modes of the real mode pencil are not dropped** (paul: some modes are leaky). The real solver now reports a *converged complex Ritz pair among the wanted eigenvalues* (`EigPairs::complex_found`, `Modes::complex_wanted`); `ems::solve_emw` then repeats the study in `ems_pml::solve_emw_complex` (complex vectors, complex n_eff). How real and complex are told apart is in ADR 33 and `arnoldi_real.rs`: the Ritz values come from the real Schur/Hessenberg form of the projected matrix, a real value has imaginary part exactly 0 and a complex one comes in a conjugate pair; the only cutoff, `|Im| > 10 tol |theta|`, is below the accuracy asked of the value itself and only guards two nearly equal real values that the projection merged; and a complex pair counts only if it has converged (residual below `tol`) and ranks among the nearest `nev + 2`. Contract test `a_converged_complex_pair_of_a_real_pencil_is_reported_and_never_dropped` (eigenvalues `1 +- 0.5i` and real ones). Truly complex operators (PML, lossy materials) use the complex path as before.
* **Golden change (paul: adopt the extra converged mode for now)**: `emw_pml_x.json` has 4 modes instead of 3. The new 4th mode: `neff = 1.730216 -5.076e-6i` (Q 1.70e5, quasi-TM z, P 1.000). Modes 1 to 3 are unchanged (neff agree to 1e-12 in the real part). Why: the old solver delivered only 3 pairs the kernel accepts (the kernel drops pairs it did not converge); I did not trace which limit of the old solver left the 4th out; the library converges it at the kernel's tolerance 1e-10. It is **tolerance-sensitive**: with `NYX_EIG_TOL=1e-13` the study returns 3 modes again, so it is marginal, not proven genuine or spurious. Details and the way back: `crates/femcore/tests/golden_frames/GOLDEN_CHANGES.md`.
* `solid_kz_pml_z` and the other kz goldens: **moot**, `sms_kz` was removed on `restructure-1b-remove-sms-kz` (rebased onto it).

**Tests on this tip** (femcore targets run after the rebase: golden_frames, stages, ring_frame, axisym_closed_form, solid_pml, linalg_contract, solid, contracts, emw2d): all pass except the 2 known SOI goldens (R22). Workspace run before the rebase and the last changes: 762 passed, 4 failed (the 2 SOI goldens, `emw_pml_x`, `solid_kz_pml_z`; the last two are now resolved), 9 ignored; doctests failed to start only because another session replaced `libfemcore` in the shared target. **Not re-run after the rebase**: the whole workspace and `photonics`.

**Kernel files touched (for the merge with phase 2 and with the kz removal)**: `crates/femcore/src/`: `dfs.rs`, `ems.rs`, `ems_pml.rs`, `ems_axisym.rs`, `ems3d.rs`, `sms.rs`, `sms_axisym.rs`, `sms_floquet.rs`, `sms_kz.rs` (imports and the eigensolve call of the axisymmetric path only), `sms_lrb.rs`, `sms_pml.rs`, `fem.rs`, `periodic.rs`, `lib.rs`, `optical.rs`/`analysis.rs`/`polarization.rs`/`problem.rs`/`fields.rs` (imports only); tests `stages.rs`, `pml.rs`, `ring_frame.rs`, `aniso_optics.rs`. The edits are about matrix types (`Csr`, `CCsr`, `.nrows()`, `indptr()`/`values()`), imports and solver calls; none touch mesh types.

**Benchmarks (NOT the required best-of-3 on an idle machine; the machine was shared)**, `bench`, old stack against new driver at the same load, before the `DirectSolver` pattern fix: 2D direct solves were 35 to 50 % slower (electrostatics 0.08 to 0.11 s, static 2D 0.26 to 0.39 s), 2D eigenfrequencies faster (1.73 to 1.06 s), 3D eigenfrequencies faster (58 to 41 s under load; idle old baseline 33.7 s), 3D direct +13 %, AMG-CG unchanged. The pattern fix (lower triangle built without sorting, full pattern lazy) is committed but **not re-measured**. Required before accepting: best of 3, idle machine, `BENCH_ONLY=1 bash tools/baseline.sh <dir>`.

**Open / left**
1. **Benchmarks on an idle machine** are still outstanding: every number in this section was taken while other sessions were running (the same 3D direct solve took 4.7 s in one run and 11.7 s in another); a quiet window is needed for best of 3 of every stage (`BENCH_ONLY=1 bash tools/baseline.sh <dir>`, `complex_bench`, `cavity_bench`). Not timed at all: `sms_pml`, axisym and `ems3d` complex at realistic size.
2. **R26 done**: the stage functions of `sms` (static, real modes), `sms_lrb` and `ems3d` are public and in contract `kernel_stages` 1.2 (`tests/contracts.rs`).
3. **Rename `solvers` to `linalg`** (plan section 2) not done: it rewrites every `crate::solvers::` import and would collide with the mesh thread's edits in the same kernel files; do it at the merge point, in one mechanical commit.
4. **Done**: `amg_elastic` builds its hierarchy on `Csr<f64>`; only `Csr32` (single-precision copy for the cycles) is private to the module.
5. The two PML-study mode counts depend on the eigen tolerance (above). `ems::mode_fields_cells`: `Ez` is not flipped with `Et` by the sign convention (latent gauge bug, not from this phase).
6. History: the commits of this bookmark below the tip are marked conflicted by jj after the rebase (the conflicts were resolved at the tip); `jj squash` of the chain into one commit would tidy it, left to the coordinator.

### Solver path per stage: old against new (full table with the reasons in ADR 33)

Evidence from the code at `restructure-foundation` and at this bookmark plus `bench`, `cavity_bench`, `profile3d`, `mem_profile`. **Specialisations unchanged**: Cholesky/LU direct solves, elasticity AMG with P2 coarse level and the AMG-CG/LOBPCG selection thresholds and fallbacks, dense transient eigensolve, real block Arnoldi for the 2D mode pencil. **Changed to a library algorithm**: shift-invert (real and complex) is now restarted block Krylov-Schur, LOBPCG is the library's. **Possible regressions**: small direct solves (+35 to 50 % before the pattern fix, not re-measured); large complex problems (PML, Floquet, lrb, axisym, ems_pml, ems3d complex), where `max_dim` is ignored and multi-RHS solves are sequential, not timed; LOBPCG 3D (faster than old under load, idle baseline unmatched); nested dissection unused for large 3D direct. **Faster**: shift-invert 2D (1.73 to 1.06 s), ems3d real (4.4 to 3.3 s, best of 3, load unknown).

### Parallel block solves in the complex stages (added after paul's question)

The library's default of sequential triangular solves comes from a faer 0.24.4 finding (its own parallel multi-RHS solve is pathologically slow); the old complex code instead ran the sequential solve on column groups in separate threads. Reinstated as `EigenOptions::par_solve_min_n` (contract `linalg_options` 1.1; default 5000 for the complex wrappers, 0 elsewhere, `NYX_PAR_SOLVE_MIN_N` overrides). `complex_bench` (`ems_pml` guide), old / new sequential / new parallel: 6.8k unknowns 0.175 / 0.178 / 0.175 s; 19.8k 0.539 / 0.588 / 0.515 s; 35.5k 0.963 / 1.033 / 0.914 s. Parallel is kept as default. Not covered: `sms_pml`, axisym, `ems3d` complex, sizes above 35k. Details: ADR 33, last section.

### Ordering rule (measured, not adopted)

Nested dissection for large 3D direct solves was benchmarked (`NYX_ORDERING=auto|amd` on `bench`, two runs each, loaded machine) and lost in every pair: 3D direct 205k unknowns 13.5 / 17.8 s with `Auto` against 4.7 / 10.9 s with AMD; small 2D solves 60 to 80 % slower with `Auto`. The rule stays AMD for all sizes (as the old stack); the best AMD run (4.68 s) matches the phase 0 idle baseline (4.84 s). Details: ADR 33.

### Final state of phase 3 on this tip (2026-11-04)

* **Whole workspace on the tip, isolated target directory** (`CARGO_TARGET_DIR` of its own, so no other session's `femcore` leaks into the doctests): 755 passed, 2 failed, 9 ignored, `photonics` and the doctests included. The 2 failures are the known SOI wire goldens (`golden_soi_wire_modes`, `golden_soi_wire_with_sloped_wall_hybrid_modes`, R22). (Phase 1 + foundation baseline: 764 passed, 2 failed, 9 ignored; the difference is the unit tests of the deleted old modules and of `sms_kz` and the kz goldens removed on `restructure-1b`, plus this phase's new tests.)
* **Measured on a loaded machine** (other sessions running; to be repeated idle): small direct solves are back at the phase 0 baseline (electrostatics 2D 0.07 to 0.08 s against 0.09; solid static 2D 0.22 to 0.25 s against 0.27); 2D eigenfrequencies 0.83 s against 1.72; EM modes 2D 0.94 s in the one quiet run against 0.89 baseline, about 2.0 s under load; 3D direct 4.68 s in the best run against 4.84 baseline (11.7 s under load); 3D eigenfrequencies 29.5 s in the best run against 33.7 baseline (43 to 49 s under load); `cavity_bench` 3.3 s against 4.4; complex PML solves (`complex_bench`) 5 % faster than the old stack with column-parallel solves.
* **Decisions recorded in ADR 33**: AMD ordering for every size (nested dissection lost in the measurement), column-parallel block solves for complex stages (`par_solve_min_n`), real Arnoldi for the 2D mode pencil with complex pairs re-solved in the complex path, tolerance per G2 option B.
* **Still open**: idle-machine best-of-3 (above); `solvers` to `linalg` rename (left for the merge point, mechanical); a real Schur form in the library would let the 2D mode pencil use the library solver (not needed now).

## Phase 2b: 3D mesh adapter (workspace `nyx2-mesh3d`, bookmark `restructure-2b-mesh3d-adapter`)

State 2026-11-04: **done on the bookmark, not merged.** Based on `restructure-2-mesh` at `f9f3af74` (ADR 32 steps 1 to 3); rebase onto its later steps. ADR 40.

* What: `meshing::{mesher3d, refine3d}` are library code again, driven by geom through the trait `Domain3` (`meshing::geomdrive3::GeomDomain3`: `Evaluated::classify` for the region at a point, `Complex::facet_project` for the exact closest point of a face, `LabelMap` for the integers). `meshing::mesh_model` now meshes a 3D geom model (Simple mesher, inline layers; adaptive and structured PML are refused with a message). Cells carry domain numbers and boundary facets (outer faces and faces between regions) the boundary numbers of the model's label map, both 1-based. The mesher returns the native `CellMesh<3>`; `MeshedGeom.mesh` is `SimplexMesh::from_cells` of it. `Solid`, `sdf` and `face_at` are gone; the old creases fitting is dropped (see ADR 40).
* Files touched outside the new ones: `meshing/mod.rs` (three `pub mod` lines and a 3D dispatch at the top of `mesh_evaluated_report`), `tests/emw3d.rs` and `examples/cavity_bench.rs` (sphere from geom: `geomdrive3::sphere_cavity`), `tests/contracts.rs` (one test appended), ADR index, `CODE_TO_ADR.md`, `contracts/VERSIONS.md` (mesh 1.2). Deleted: `tests/common/solid3d/`, `meshing/README.md`.
* **No change in geom.** The restrictions it states stay: a curved body must sit inside a prismatic host and not touch anything, a lone sphere or torus has no host (`stack.rs:203`), mixed sweep axes are refused (`scene.rs:255`). A PEC cavity is a sphere in a *hole* host box, which geom accepts.
* Tests (all pass, own target dir): `tests/mesh3d.rs` 6/6 (sphere in a box: two domains, 7 boundary labels, box faces exact, sphere vertices on the sphere to 1e-9, conformity of the facets; torus: volume 0.9 % and area 0.2 % off the closed form at h 0.1; cavity in a hole host; a `Size` rule refines a region; 2D-only mesher options refused), `emw3d` 6/6 (the PEC sphere within 1 % of the closed form on a geom-made mesh), the contract test `mesh_geom3d_numbers_cells_and_facets_with_the_model_label_map`. The rest of the suite was **not** re-run.
* Not done: thread-count determinism test and benchmarks for the new mesher; crease fitting; 3D adaptation and wavelength sizing; mixed-axis and free-standing curved bodies (geom work).
