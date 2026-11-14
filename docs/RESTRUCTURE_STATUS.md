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
| 4 | nyx-physics content, kernels into the physics crates | **done** (awaiting acceptance) | `restructure-4-physics` | the interface exists (F); gate G3 answered |
| 5 | nyx-analysis | waits for 4 | | |
| 6 | Materials | **implemented** (+ extends, more models, provenance, checked library), awaiting paul's answers (defaults in its section) | `restructure-6-materials` | ADR 29, 30, 34, 35; `docs/MATERIALS.md` |
| 7 | nyx-engine, MCP parity, thin app | **implemented**, awaiting acceptance (see its section at the end) | `restructure-7-engine-frontends` (on the merge of 5, 6, 7a) | ADR 57 to 60; 7a (connect merge) is merged into it |
| 8 | File format | **implemented**, awaiting acceptance (see its section at the end) | `restructure-8-format` (on `main` = phase 7) | ADR 63; tag `last-old-format` |
| 9 | Close-out | **implemented**, awaiting acceptance (see its section at the end) | `restructure-9-closeout` (on `main` = phase 8) | ADR 64, 65 |

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

## Phase 7a: connect merge (session `connect`, bookmark `restructure-7a-connect`) - done

`app/src/connect.rs` is deleted; `crates/nyx-connect` is the only copy (R07 resolved, ADR 28). `app` depends on it; the layering table needed no change (`nyx -> nyx-connect` was already allowed).

* **Crate**: modules `env`, `hosts`, `config`, `server`, `report`; `--json` on the command; `examples/fake_mcp.rs` is the stand-in server for the handshake tests (built by `cargo test`, not shipped).
* **Decisions taken in the merge** (details in ADR 28): backup is `<name>.nyx-backup`, made once, never overwritten (the old time-stamped copy could be overwritten within a second); Claude Code's entry has `type`, `args`, `env: {}`; keys the user added to the `nyx` entry survive a reconnect; `find_server` = `NYX_MCP`, then the folder of the running program (that is all a shipped build does); only in a cargo build tree (`is_cargo_build_tree`, detected by the `deps` and `.fingerprint` folders, not by `debug_assertions`, because developers run `--release`) it also tries the program's folder and four parents (five folders) and their `target/release`; the Connect menu now lists every installed host (five), not only Claude Code and Desktop.
* **Contracts** (new, `VERSIONS.md` and the README table updated): `connect_hosts` v1 and `connect_report` v1 (`connect_report.v1.json` examples plus `connect_report.v1.schema.json`, generated by `report::schema()`; the examples and the real program's output are validated against it).
* **Tests**: `cargo test --release -p nyx-connect`: 23 behaviour + 6 contract tests pass. Whole workspace `cargo test --release --workspace --no-fail-fast`: 782 passed, 2 failed (the two known `ems` goldens, R22/G2), 9 ignored. Layering test and `check_code_to_adr.sh` (564 files, 0 without a row) pass. The app builds (`cargo build --release -p nyx`). The Connect menu was **not** opened by hand in the window.
* **Docs**: ADR 28, `ADR_INDEX.md`, `CODE_TO_ADR.md` (the `app/src/connect.rs` row is gone), `ROUGH_SPOTS.md` R07 marked resolved, ADR 15 notes the partial supersession.
* **Left for others**: `workspace.dsl` line 218 ("Connect menu duplicates this logic") and the `nyx.app -> nyx.connect` relationship are stale; that is the `dsl` session's file. Phase 7 proper still has `nyx-mcp`/`nyx-session` calling nothing in `nyx-connect` (the table allows the edge when the engine arrives).
* **Watch out when merging bookmarks**: `Cargo.lock` (new `nyx -> nyx-connect` edge) and `crates/app/Cargo.toml`/`app.rs` (Connect menu hunk only) are the only shared files I touched. A stale-fingerprint quirk made cargo reuse an old `nyx-connect` rlib once in the shared target dir; touching the sources fixed it.

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

## Phase 4 (physics session, workspace `nyx2-phase4`, bookmarks `restructure-4-physics`, `restructure-4a-base`)

Updated 2026-11-04. **State: the merge, the rename, the shared groundwork, `dfs` and the removal of `SimplexMesh` from the physics interface are done (bookmark `restructure-4-physics`; `restructure-4a-base` is the commit before the mesh change); `sms` and `ems` are not moved yet (two threads).** Read ADR 41 first: it is the pattern and the checklist for the next two physics.

**Done**
* **Merge** of `restructure-2-mesh` and `restructure-3-linalg` (conflicts only in docs and `femcore/src/lib.rs`), bookmark `restructure-4-physics`. `femcore::solvers` renamed `femcore::linalg` (module, `tests/linalg`, imports, docs). `restructure-2b-mesh3d-adapter` rebased onto the merge (ADR 40, one test-file conflict resolved); `restructure-4e-engine-core` merged in (engine core with the stale mesher paths fixed). The line with all of it is `restructure-4a-base`.
* **Shared groundwork** (ADR 41): `femcore::props` (materials, PML, periodic pairs, linear solver; `problem` re-exports), `fem::{check_boundary_labels, cell_values_at_nodes}`, `femcore::testkit` (mesh fixtures; `tests/common` re-exports), `nyx_physics::{decode, testkit}` (typed parameters, name and selector resolution, the golden harness), `Field::scalar/vector`, `Solution::flat_fields/flat`, `Discretisation` and `discretise`, `check_result_finite`, `femcore::adapt::PlotFields` (the adaptation loops take any result), `nyx_engine::legacy` (bridge for the window and the MCP session).
* **dfs moved** into `crates/nyx-dfs` (electrostatics, stationary and transient heat): typed `Problem` and `decode`, the three stages, `run`, `Physics::solve`, descriptors with the new `solver` option (contract `physics_dfs` 1.1). `femcore::dfs` is deleted; `femcore::solve` answers an error for the dfs arms. `app` and `nyx-session` call `nyx_engine::legacy::solve[_with]`. The engine solves `dfs` for real through `Physics::solve` (`tests/engine.rs`, `dfs_solves_the_strip_through_the_engine`).
* **Tests moved with the kernel, numbers unchanged**: the four dfs goldens (files moved unchanged, tolerance 1e-9), `scalar.rs` (9 closed-form tests), 12 more from `amg_solve`, `anisotropic`, `finite_guard`, `frame`, `frame_roundtrip`, `geom_meshes`, `wafer_names`, 3 stage tests, 2 meshing tests; contracts: envelope against typed problem bit for bit, refusals, cancel, stage API. Tests that need two physics are in `nyx-engine/tests/legacy_thermal.rs` (heat feeding thermal expansion, with the golden `solid_static_thermal_from_heat`, whose file is still in `femcore/tests/golden_frames`). The old-shape `solve_result` contract tests are in `nyx-engine/tests/legacy_contract.rs`. `bench` and `stage_log` moved to `nyx-engine/examples` (`tools/baseline.sh` runs `bench` from there).
* **Checks**: `cargo test --release --workspace --no-fail-fast`: **774 passed, 2 failed, 9 ignored** (the two failures are the known SOI-wire goldens, R22; baseline of the merged line before the dfs work: 771 / 2 / 9). `cargo clippy -D warnings` and `cargo fmt` clean for the new crates. `check_code_to_adr.sh` passes. Layering test passes (new temporary edge `nyx-session -> nyx-engine`).

**Rules for the sms and ems threads**
* Do not edit what the other physics owns. `femcore/src/lib.rs` (your `mod` lines), `femcore/src/solve.rs` (your arms), `nyx-engine/src/legacy/mod.rs` (your arm; your file is `legacy/<id>.rs`), `docs/architecture/{ADR_INDEX,CODE_TO_ADR}.md`, `contracts/{README,VERSIONS}.md` and this file are the shared hot spots: edit your own lines only, rebase often.
* `femcore::periodic` (Floquet pairs) is used by `sms_floquet` and `ems3d`: it stays in femcore; neither thread moves it. The same for `section`, `material`, `fields`, `fem`, `optical` (until phase 5) and `polarization`/`analysis` (phase 5).
* Name the mesh only as `nyx_physics::Mesh` and `nyx_physics::Discretisation`; take `Discretisation` through `nyx_physics::discretise(mesh)` at the top of `solve`.
* The golden harness compares through flat names; a vector (`E`, `u`) is one `Field` with components, `Solution::flat_fields` expands it. A solid result's `disp` array is the vector field `u`; the bridge (`legacy::result_back`) leaves `disp: None`, so the sms thread must fill it in its own `legacy/sms.rs` (and `neff` comes from `Solution::eigenvalue`).

**Decisions and deviations**
* **Legacy bridge** (coordinator agreed): the window and the MCP session still build the old `ProblemSpec`; they reach moved kernels through `nyx_engine::legacy`. It goes with `femcore::problem` and `femcore::solve` in phase 7.
* **`SimplexMesh` removal** (coordinator agreed to this reading): **done** (commit f326636a): `Physics` takes `femcore::mesh::Mesh`, the native cell mesh in 2D or 3D (`mesh/exchange.rs`), and each solve builds the numbering `Discretisation` with `nyx_physics::discretise` (a unit test pins that the round trip reproduces the solver mesh exactly); the struct `SimplexMesh` itself stays in femcore for the mesher and the window until phase 7; non-simplex cells are still cut into simplices so goldens keep their numbers; native quadrilaterals, prisms and hexahedra are a later, golden-changing step. ADR 41 item 5.
* The adaptation loops became generic (`PlotFields`); a signature change of `solve_adaptively[3]`, no change in behaviour.
* Contract `physics_dfs` 1.1: the physics option `solver` replaces `ProblemSpec::solver`.

**Open / next**
1. `sms` and `ems` (two threads from `restructure-4a-base`); each deletes its `femcore` module and arm, moves its goldens and tests, bumps its `physics_<id>` exemplar.
2. (done) `SimplexMesh` removal from the interface. Left for phase 7: rename the struct to `Discretisation` and let the mesher produce `Mesh` directly when the window is rewired.
3. `Physics::sizing` for sms and ems, then delete `physics_size::spec_sizing` and `mesh_for_problem`.
4. Rewrite the contracts for the envelope (`kernel_stages`, `problem_spec`) once all three physics are in; delete `problem_spec.v1`.
5. Idle-machine benchmarks against a control run of the parent commit, best of three (not taken yet; the machine is shared).

## Phase 4: ems (workspace `nyx2-ems`, bookmark `restructure-4-ems`, started from the tip of `restructure-4-physics` = a5ba6825)

Updated 2026-11-05. **State: `ems`, `ems_pml`, `ems_axisym` and `ems3d` are in `crates/nyx-ems`; the tests and goldens moved with them. ADR 45, 46, 47.** Read ADR 45 first.

**Done**
* **Kernels moved** (formulas and order of operations unchanged): `nyx_ems::{modes, modes_complex, axisym, cavity, optical}`, typed `Problem`/`Study` (`modes`, `axisymmetric_modes`, `cavity`: the study names the kernel, R10) with `decode`, `run`/`run_analysed`, `Physics::solve`, descriptors refined (`physics_ems` 1.1), result-shape exemplar `physics_ems_result` 1.0. `femcore::{ems, ems_pml, ems_axisym, ems3d}` are deleted; `femcore::solve` answers an error for `Emw`.
* **Stays in femcore, new modules**: `edge` (2D edge-element basis), `pml` (`stretch`, `cell_stretches`), and in `optical` the complex material tensors (`EpsC`, `MuC`, `cell_materials`) and `C0`; `section::axes_of_3d_periodic`, `Analysis::summary`. Reason: `optical`/`analysis`/`polarization` stay until phase 5 and `femcore` cannot depend on `nyx-ems`; the solid PML kernels share `pml`. **For the sms merge**: `sms_pml.rs` and `sms_kz.rs` import `crate::pml::cell_stretches` instead of `crate::ems_pml::cell_stretches`.
* **Polarization**: still computed in `nyx-ems` through `femcore::analysis`; the scalars `polarization N` and the table are in the result (the goldens record them, so they keep their numbers); `run_analysed` returns the `Analysis` for the window. Phase 5 removes the scalars and rewrites those goldens on purpose.
* **Bridge**: `nyx_engine::legacy::ems` (2D + `Modes` = waveguide modes, 3D + `Eigenfrequency` = cavity, `analysis` set on the old result). **Use cases** take their solver as an argument (`usecases::Solver`), ADR 47.
* **Goldens moved unchanged** to `crates/nyx-ems/tests/golden/`: `emw_pec_straight`, `emw_pec_tilted`, `emw_pml_x` (1e-7, as before), and the two SOI wire files `emw_soi_wire`, `emw_soi_hybrid`, which never reproduced on this machine (R22) and now use the G2 tolerances (n_eff 5e-5 relative, fields 1e-3 of their maximum, polarization class compared as text). The entry of `GOLDEN_CHANGES.md` about `emw_pml_x` moved with them.
* **Tests moved, same numbers and tolerances**: to `nyx-ems/tests/`: `emw2d`, `emw3d`, `pml`, `aniso_optics`, `bragg_stack`, `photonics`, `physics_meshing`, `meshing` (the four SOI tests), `moved` (anisotropic guide, frame guide, stage composition, wafer names, Floquet cell on an adapted mesh), `axisym` (coax cavity and the whispering-gallery power direction), `golden`, `contracts` (envelope against typed bit for bit incl. a Floquet pair, refusals by name, cancel, stage API, result shape). To `nyx-engine/tests/usecases.rs`: the use-case tests (they need the old result with its analysis). Examples `soi_modes`, `bragg_gratings`, `cavity_bench`, `complex_bench` are in `nyx-engine/examples`; `tools/baseline.sh` runs `cavity_bench` from there.
* **G1 (axisymmetric kernel checked)**: coax TEM modes 1.4e-6 and 8e-6; PEC cylinder TM modes at m = 1, 2, 3 against the Bessel zeros 3.5e-4 (1 mm axis rod; 1.3e-3 with a 2 mm rod, so the rest is the rod), 1.1e-5, 2.7e-5; `m` and `-m` give the same frequency and opposite power along y. **Not verified**: a dielectric whispering-gallery mode against a reference, the absorbing layers of the axisymmetric kernel, order 1, large `m`. Details in ADR 46. The kernel now has a study (`axisymmetric_modes`) and a fields stage.
* **Checks** (isolated target directory, because the shared one had a locked test executable from another session): `cargo test --release --workspace --no-fail-fast`: **783 passed, 1 failed, 9 ignored**. The failure is `geom --test props union_area_is_associative_and_commutative` (a randomised property test of `geom` on a 1e-9 sliver polygon, 7.565216931 against 7.565216924; `geom` is untouched by this phase). The two SOI golden failures of the baseline are gone (G2 tolerances). `cargo clippy -D warnings` and `cargo fmt` clean for `nyx-ems`. Layering test and `check_code_to_adr.sh` pass.

**Merge notes for the phase 4 thread**
* Shared hot spots touched (my lines only): `femcore/src/lib.rs` (the `ems*` mods replaced by `edge`, `pml`, and the doc line), `femcore/src/solve.rs` (the `Emw` arm), `femcore/tests/contracts.rs` (three lines of `kernel_stage_api_is_public`), `femcore/tests/golden_frames.rs` (the four ems tests, `crystal`, `C0`), `femcore/tests/{stages,anisotropic,frame_roundtrip,ring_frame,adapt3d,axisym_closed_form,meshing}.rs` (one ems test each), `nyx-engine/src/legacy/mod.rs` (one arm), ADR index, `CODE_TO_ADR.md`, `contracts/{README,VERSIONS}.md`. `femcore/tests/wafer_names.rs` held only the ems test and is deleted. Expect conflicts in the sms thread's edits of the same test files.
* `CLAUDE.md` still says `nyx-ems` is a skeleton and lists the femcore examples; update it after both physics are in.

**Open / next**
1. `Physics::sizing` for ems (wavelength per region), then delete the ems part of `physics_size::spec_sizing`.
2. Phase 5: move the polarization out of the result into `nyx-analysis`; the goldens lose the `polarization N` scalars then.
3. Phase 6: the envelope carries scalar permittivity only; crystals through material files.
4. Idle-machine benchmarks (`cavity_bench`, `complex_bench`) against a control run: not taken; the machine is shared.

## Phase 4: sms (workspace `nyx2-sms`, bookmark `restructure-4-sms`, based on `restructure-4-physics` at `f326636a`)

Updated 2026-11-04. **State: `sms` is moved into `nyx-sms`; goldens unchanged.** ADR 42 (the move, the typed problem, the result shape), ADR 43 (bodies of revolution stay a library entry), ADR 44 (what stays in femcore and the few shared lines touched).

**Done**
* `crates/nyx-sms`: `problem` (typed `Problem`, `decode`), `solid` (was `sms`), `pml`, `lrb`, `floquet`, `axisym`, `cross_section` (was `sms_kz`), `run`, `Physics::solve`. Kernels edited for types only. `femcore::{sms, sms_pml, sms_lrb, sms_floquet, sms_axisym, sms_kz}` are deleted; `femcore::solve` answers an error for `Physics::Solid`.
* Bridge `nyx_engine::legacy::sms`: old problem to typed, `FromHeat` to `Imported(state.heat_t)`, result back, **`disp` filled** from the vector field `u` (real studies) or `Re u{x,z,y}` (complex modes). `neff` stays `None` for sms. The arm is in `legacy::solve_with`.
* Descriptors `physics_sms` 1.2 (option `solver`; property `stiffness` removed), new exemplar `physics_sms_result.v1.json` 1.0, `kernel_stages` 1.3.
* **Goldens moved unchanged and passing**: `solid_static_sapphire_strain`, `solid_static_sapphire_stress_thermal`, `solid_eigen_sapphire`, `solid_eigen_floquet`, `solid_eigen_pml_x`, `solid_eigen_pml_z`, `ring_model_axes` to `crates/nyx-sms/tests/golden/`; `solid_static_thermal_from_heat` to `crates/nyx-engine/tests/golden/`. Nothing was regenerated.
* **Tests moved** to `crates/nyx-sms/tests`: `solid`, `solid_pml`, `anisotropic` (solid half), `stages` (4 tests), `ring_frame` (5 of 6), `axisym_closed_form` (ring), `amg_solve`, `cancel_and_log`, `determinism`, `finite_guard` (+ the finite result check of the new result), `frame` (2), `frame_roundtrip` (1), `geom_meshes` (plate with a hole), `golden`, `split` (whole file: a solid bar is its oracle). New in `contracts.rs`: envelope against typed problem (bits), layer and Floquet pair through the envelope, imported temperature, refusals, stage API, result-shape exemplar. femcore keeps the ems halves of those files. Benchmarks `eig_beam`, `mem_profile`, `profile3d` moved to `nyx-engine/examples` (`tools/baseline.sh` updated).
* **Checks** (`cargo test --release --workspace --no-fail-fast`, shared machine): **787 passed, 2 failed, 7 ignored**; the two failures are the known SOI-wire goldens (R22). Baseline of `f326636a`: 775 / 2 / 9 (the ignored count differs because the runs were taken on different targets; not investigated). `cargo clippy -D warnings` and `cargo fmt` clean for `nyx-sms` and `nyx-engine`; `check_code_to_adr.sh` passes.

**Axisymmetric kernel check (G1)**: the ring tests (`axisym_closed_form.rs`, `ring_frame.rs`) pass on the moved code with their thresholds (breathing mode 2e-3, flexure 2e-2); the figures quoted in phase 3 were not re-measured.

**For the merge (see ADR 44)**
* Shared lines I edited: `femcore/src/{lib,solve,progress,validate,fem,ems3d,h5io}.rs` and `linalg/amg_elastic/mod.rs` (`sms::lap` became `progress::lap_since`; if the ems move deletes or edits `ems3d.rs`, take its version and call `progress::lap_since`), `nyx-engine/src/legacy/mod.rs` (one `pub mod`, one arm), ADR_INDEX, CODE_TO_ADR, contracts README and VERSIONS, `tools/baseline.sh`.
* `nyx-sms` calls `femcore::ems_pml::cell_stretches`: it must stay public in femcore (or move to `props`).
* Test files edited on both sides: `femcore/tests/{golden_frames,stages,anisotropic,ring_frame,axisym_closed_form,frame_roundtrip,contracts}.rs` (my hunks removed the solid parts only; ems parts untouched).
* `CLAUDE.md` still calls `nyx-sms` a skeleton; left for the merge since it says the same of `nyx-ems`.

**Open / not done**
1. `Physics::sizing` for sms is the skeleton's (no wavelengths); `physics_size::spec_sizing` stays.
2. Axisymmetric solid study is not in the envelope (ADR 43).
3. Anisotropic stiffness cannot be given through the envelope until the material files (phase 6).
4. Descriptor `fields` for the eigenfrequency study lists `u`; a lossy mode has `Re ux`, `Re uz` (ADR 42 point 4).
5. Benchmarks of the solid examples not re-timed (shared machine).

### Phase 4 close (physics session, 2026-11-05): merge of sms and ems, checks, benchmarks

* **Merge** of `restructure-4-sms` (e13f322f) and `restructure-4-ems` (10baedf0) into `restructure-4-physics`. Conflicts were the ones ADR 44 and the ems notes predicted: `cell_stretches` is `femcore::pml::cell_stretches` (the ems thread made it; `nyx-sms` imports it from there); the `lap` change touched `ems3d.rs`, which moved into `nyx-ems`; each thread had left the other physics' half in the shared femcore test files (`anisotropic`, `axisym_closed_form`, `contracts`, `frame_roundtrip`, `golden_frames`, `stages`): both halves had moved into their crates, so they were dropped and the emptied files deleted. `femcore::solve` had no live arm left and is deleted; `nyx_engine::legacy` runs all three physics. `femcore/src` no longer holds a physics kernel (`ems3d`, `sms_kz`, `sms_pml` deleted). `tools/baseline.sh` runs `bench`, `profile3d`, `mem_profile` and `cavity_bench` from `nyx-engine` (the examples that drive every physics).
* **Checks on the tip**: `cargo test --release --workspace --no-fail-fast`: **790 passed, 0 failed, 9 ignored**. The two SOI-wire goldens pass at the G2 tolerances (the ems thread's work: 5e-5 relative on n_eff, 1e-3 on field maxima), every other golden is the original file at its original tolerance. `cargo clippy -D warnings` clean for the six new crates, layering test passes, `check_code_to_adr.sh` passes.
* **Benchmarks** (plan section 5 item 4: best of three, control run of the parent commit `718cd218` = the merge of phases 2 and 3 with the rename, interleaved control/new, `BENCH_ONLY=1 tools/baseline.sh`; the machine is about 1.8x slower than at the phase 0 baseline and was not otherwise idle, so single values vary, e.g. 3D eigensolve 41 to 66 s). Best of three, control against phase 4 tip, wall seconds:

| Case | control | phase 4 | change |
|---|---|---|---|
| `bench` total | 63.0 | 62.3 | -1 % |
| electrostatics 2D (38,697 dofs) | 0.08 | 0.09 | +0.01 s (see below) |
| solid static 2D | 0.23 | 0.22 | -4 % |
| solid 6 eigenfrequencies 2D | 0.83 | 0.81 | -2 % |
| solid static 3D | 2.29 | 2.29 | 0 |
| solid 6 eigenfrequencies 3D | 41.91 | 40.99 | -2 % |
| solid static 3D direct | 11.26 | 11.23 | 0 |
| solid static 3D AMG-CG | 3.56 | 3.62 | +2 % |
| EM modes 2D (Nedelec 2) | 1.94 | 1.98 | +2 % |
| `profile3d 16 16 16 p2` | 2.15 | 2.15 | 0 |
| `mem_profile 12 12 12 p2 static auto` | 1.26 | 1.33 | +5 % (a 1.3 s process, start-up dominated) |
| `cavity_bench 0.012 2` | 4.03 | 4.12 | +2 % |
| `solve_scaling` | 18.1 | 18.2 | +1 % |

  Everything is inside the 10 % rule. The one visible difference is electrostatics 2D, 0.08 to 0.09 s: the window bridge now converts the solver mesh to the exchange mesh and each solve builds its numbering again (`Mesh::from_solver_mesh`, `discretise`), a fixed cost of about 10 ms on a 77,000-cell mesh that does not grow with the solve. Raw logs: `C:\Users\robin-adm\Documents\opomech\nyx2-bench` (outside the repository; delete when accepted). The control workspace `nyx2-phase4-ctl` was forgotten and its folder removed; its build products are in the shared `nyx2-target`.
* **Plan definition of done** (section 3, phase 4): the kernels of `dfs`, `sms`, `ems` are in their crates with typed problems decoded from the envelope (R10: the EM study is named in the study; R12: property names are validated against the descriptors by `validate_spec`; R19: the explicit `Ctx` replaces the thread-local); `registry.get(id).solve(...)` solves (engine test `dfs_solves_the_strip_through_the_engine`, the coupled chain, the use-case tests); nothing outside the physics crates and `nyx-engine` names a physics (layering test). Not done, on purpose: the old window problem (`femcore::problem`) and the bridge stay until phase 7; `problem_spec.v1.json` stays until then, `problem_spec.v2` and the `physics_*` exemplars are the envelope contracts; `physics_size::spec_sizing` and `mesh_for_problem` (adapters of the old problem) stay until the window is rewired; native quadrilateral/prism/hexahedron elements (golden-changing) are a later step.

## Phase 4m: adaptive meshing, wavelength sizing and structured PML in 3D (workspace `nyx2-mesh3d-adapt`, bookmark `restructure-4m-mesh3d-features`)

State 2026-11-05: **done on the bookmark, not merged.** Started from `restructure-4-physics` at `a5ba6825`, rebased onto the closed phase 4 (`45d172f7`; the only conflicts were the shared doc tables, resolved by keeping both sides). Paul: all three are essential in 3D. Stayed inside `femcore::{meshing, adapt}`; no physics crate was edited (one test file is new in `nyx-ems`). ADRs 50 to 55 (50 is the first free number above the sms/ems ones).

**What works now** (3D geom models, any dimension through `meshing::request::mesh_request`; `MeshPlan` as before):
* **Wavelength sizing** (ADR 50): `meshing::mesh_for_sizing(list, plan, &Sizing, &PhysicsSizing)`. The simple mesher gets one element size per region (a ball with short waves is refined, the air around it is not); the adaptive mesher intersects the wave metric with the geometry metric. `mesh_for_problem` stays 2D.
* **Adaptive meshing** (ADR 51): `Mesher::Adaptive` in 3D = lattice start mesh + `adapt3` with the exact geom faces as the oracle (`meshing::adapt3d`); every vertex of a refined sphere is on the sphere, the box faces stay exact, labels and interface facets are kept.
* **Structured PML** (ADR 52): `PmlMode::Structured` in 3D (`meshing::pml_geom3`): the physical window is meshed alone, `build_pml` appends prisms and hexahedra, layer cells carry geom's piece numbers, outer faces geom's face numbers; physical vertex and cell ids are kept (contract test); native `CellMesh<3>` and the tetrahedral cut are both available.
* **Maxwell residual estimator** (ADR 53): `adapt::maxwell` (`CurlField` trait, `maxwell_indicators`, `indicator_metrics`) and `estimate3::adapt_loop3` (the loop with the estimator and the adaptation as arguments; `solve_adaptively3` is a call of it). PEC box cavity, lowest mode: order 1 reaches 5.2e-4 on 333 tetrahedra (uniform blocks 5.3e-3 on 696; the Hessian loop was worse than uniform), order 2 reaches 8.0e-6 on 1402 (uniform 6.7e-5 on 1470). Test in `crates/nyx-ems/tests/maxwell_estimator.rs` because it needs the cavity kernel.
* **Data** (ADR 54): `MeshRequest` / `MeshReport` with JSON Schemas, exemplars in `docs/architecture/contracts/` (`mesh_request*.json`, `mesh_report*.json`), versions in `VERSIONS.md` (mesh 1.3, mesh_request 1.0, mesh_report 1.0, adapt_maxwell 1.0).

**Faults found in code that was already there, and fixed** (they were in the way, details in ADR 51 and 55):
* the lattice mesher made zero-volume slivers wherever a model plane lies on a lattice plane (3300 cells of a layer stack: 140 slivers; they turned the PML prisms inside out): vertices a few ulp apart are welded;
* its meshes depended on the thread count (cell order by block; unstable dedup of equal cut ids): now identical for 1, 3 and 8 threads (test with child processes);
* `adapt3` could fold the mesh over itself: `edge_removal_pass` turned ring triangles inside out and duplicated edges, `flip_pass` did not know the faces and edges made earlier in its own pass. `Work3::check` (`FEMMESH_CHECK=1`) now finds overlaps; before, a box of volume 64 came out with 65.95;
* `adapt3` never removed dead tetrahedra from its arrays and repeated sweeps that changed a handful of cells (`Work3::compact_tets`, `AdaptOptions3::sweep_frac`; default behaviour unchanged);
* `pml_box::build_pml`'s `orient` copied all the points of the mesh for every layer cell (quadratic): 34 s to 0.9 s for 600 000 cells.

**Checks.** `cargo test --release --workspace --no-fail-fast` on the rebased tree: **811 passed, 0 failed, 11 ignored** (phase 4 closed at 790 / 0 / 9; the difference is the new tests, 2 more ignored studies). New tests: `femcore` `mesh3d_features` (12), `contracts` (+4), `adapt::maxwell` unit tests (2), `nyx-ems` `maxwell_estimator` (3 + 2 studies, ignored). `check_code_to_adr.sh` passes.

**Benchmarks** (`crates/femcore/examples/mesh3d_bench.rs`, run best of three, interleaved with a control run of the parent commit `a5ba6825` in a second workspace with its own target directory, nothing else compiling; wall seconds; the solver benchmarks of `tools/baseline.sh` do not touch changed code except `cavity_bench`, measured here):

| Case | control (parent) | this bookmark | Note |
|---|---|---|---|
| lattice, sphere in a box, h 0.1 (794 016 cells) | 0.219 | 0.205 | noise |
| lattice, layer stack, h 0.05 | 1.164 | 1.222 | +5 %, inside the spread of the runs (1.16 to 1.38 / 1.22 to 1.27); 3 816 slivers fewer (6 250 216 to 6 246 400 cells) |
| sphere cavity 0.05 / 0.004 (105 120 cells) | 0.045 | 0.041 | same mesh cell for cell |
| `add_pml_3d`, block 40 layers, 605 820 cells | 34.3 | 0.86 | the quadratic copy in `orient` |
| `adapt3` on a block, graded (17 425 / 17 419 cells, qmin 0.588 both) | 9.40 | 8.48 | -10 % |
| `cavity_bench 0.012 2` (whole example) | 3.204 | 3.170 | |

Features that the parent refuses (no control), best of three: wavelength sizing with the simple mesher 0.050 s for 10 608 cells; adaptive sphere 0.6 x 0.5 5.1 s for 40 500 cells; adaptive plus waves 0.68 s for 3 399 cells; structured PML on a 479 232-cell mesh 4.2 s (inline layers of the same model: 6.1 s for 831 408 cells); Maxwell indicators on 105 120 cells 0.18 s. A strongly graded adaptation (waves 0.15 in a ball, start 0.8, 17 000 tetrahedra) takes about 30 s: `adapt3` is serial and rebuilds its incidence per sweep.

**Not done / for others**
* The engine and the window do not call `mesh_request`/`mesh_for_sizing` yet (`FemcoreMesher::mesh` still answers Unsupported for physics sizing); that is a one-function change in `nyx-engine`.
* Physics crates should implement `CurlField` for their edge-element solutions (`nyx-ems`'s `cavity::{Ned3, basis3}` is what the test adapter uses) and report `Sizing`; the estimator skips absorbing layers and works for real isotropic media only.
* `adapt3` speed (threads, incidence per sweep), surface flips on curved patches, curved creases from the exact edges (polylines of the start mesh are used), goal-oriented adaptivity: `docs/ports/meshers/FEATURE_DEPENDENCIES.md` (copied into this repository, annotated) lists them.
* Structured layers need a box-shaped physical window; round windows and non-box bodies stay open (P6).
* Practical: cargo derives the same artifact hash for the same crate in every jj workspace, so workspaces that share `../nyx2-target` overwrite each other's `libfemcore` (my first test run linked another session's library and failed on missing symbols). I used private target directories `nyx2-target-adapt` (about 4 GB) and `nyx2-target-adapt-ctl`; both can be deleted, and so can the workspace `nyx2-mesh3d-adapt-ctl` (`jj workspace forget nyx2-mesh3d-adapt-ctl`).

## Phase 5: nyx-analysis (workspace `nyx2-phase5`, bookmark `restructure-5-analysis`, started from `restructure-4m-mesh3d-features`)

Updated 2026-11-05. **State: done on the bookmark, not merged.** ADR 56. Read it first.

**Merge step.** `restructure-4m-mesh3d-features` is a linear descendant of the closed phase 4 (`restructure-4-physics` = `45d172f7`), so the line is a fast-forward and phase 5 is cut from the 4m tip (`3ecdd978`); no bookmark was moved. `restructure-4-sms`, `restructure-4-ems`, `restructure-4a-base`, `restructure-4e-engine-core`, `restructure-2b-mesh3d-adapter`, `restructure-3-linalg` and `restructure-2-mesh` are all ancestors of it. **Not in the line** (other phases, each cut from older commits, left alone): `restructure-6-materials`, `restructure-7a-connect`, `restructure-3x-g2-experiment`, `restructure-3-linalg-history`. Control run on the clean 4m tip (own workspace, own target directory): `cargo test --release --workspace --no-fail-fast` 811 passed, 0 failed, 11 ignored.

**Done**
* `nyx-analysis` holds `polarization` (was `femcore::polarization`) and `modes` (was `femcore::analysis`, its struct `Analysis` is `ModeAnalysis` now); both files are deleted from `femcore` with `optical::mode_samples`. Their unit tests moved with them. New `samples` reads the modes of a `SolveResult` as weighted samples; `polarization(mesh, result)` and the `PolarizationAnalysis` implementation of the `Analysis` trait are the entry points. Layering: `nyx-analysis` may depend on `femcore` (mesh geometry, quadrature).
* **A result carries what the analysis needs, not the analysis.** `nyx_physics::FieldFrame` (new `frame.rs`) is the frame of the vector fields: wafer axes of the lateral, normal and propagation direction, handedness, absorbing domains left out; the ems mode and cavity studies write it as the scalar `field frame`. The complex field is in the cell fields (already so); the cavity solutions got their complex frequency as `eigenvalue`.
* **The `polarization N` scalars are gone** from every ems result, with `Solved`, `run_analysed`, `with_analysis`, the table in the log, `axisym::mode_samples` and the `optical_set_*` helpers. `nyx_ems::run` returns the result.
* **Window / MCP**: `nyx_engine::legacy::solve_with_analysis` returns the old-shape result and the `ModeAnalysis`; `ResultView::analysis` replaces `SolveResult::analysis` of the old window result; the window's log keeps the table; `Solution::te_fraction` (set by the bridge) feeds the use cases (`fundamental_neff`) and the run-table metric. `nyx_engine::analysis` re-exports the crate.
* **R28** resolved: samples are `lat`/`nor`/`prop`, documentation maps them to the wafer axes. The labels of the Analysis tab are unchanged.
* **Goldens (deliberate, listed in `crates/nyx-ems/tests/golden/GOLDEN_CHANGES.md`)**: the five mode files lose `polarization N` and gain `field frame`; the removed lines live on unchanged in `polarization.json` and the analysis is checked against them (class as text, numbers to 1.5e-3, 5e-3 for the silicon-wire files). **Also found and fixed**: the real mode solver did not flip `E_y` with the transverse field when it makes the largest component positive, so the relative sign of `E_y` between modes was arbitrary; the analysis of a degenerate group reads it. Fixed in the plot and cell field; the `Ey` probe sign changes in 3 + 1 + 1 modes of three files, nothing else changes (`emw_pec_tilted`, `emw_pml_x` unchanged).
* Contract `physics_ems_result` 1.1 (exemplar rewritten on purpose, `VERSIONS.md`). ADR 56, ADR index, `CODE_TO_ADR.md`, `ROUGH_SPOTS.md` (R28), `workspace.dsl` note, `CLAUDE.md`, `docs/photonic_use_cases.md` updated.

**Checks.** `cargo test --release --workspace --no-fail-fast`: **815 passed, 0 failed, 11 ignored** (the closed line has 811; 4 new tests are in `nyx-analysis/tests/result_samples.rs` and 1 for the frame text in `nyx-physics`, so 5 added; the other side of the sum was not itemised). `cargo clippy -D warnings` and `cargo fmt` clean for the six new crates; layering test and `check_code_to_adr.sh` pass. The Analysis-tab behaviour is tested through the same assertions as before: `nyx-ems/tests/{photonics,emw3d,golden}.rs`, `nyx-engine/tests/usecases.rs` (degenerate pair resolved, hybrid wire tilted); the slow Bragg test (`--ignored`) was not run. The egui window was compiled and its own 24 tests pass; the Analysis tab was **not looked at on screen**.

**Benchmarks** (best of three, control = the clean 4m tip in its own workspace and target directory, runs interleaved control/new, whole machine otherwise idle apart from nothing compiling; the first cavity_bench pair is dropped because the control compiled it). Wall seconds: EM modes 2D (Nedelec 2, `bench`) 2.18 control against 1.90 (-13 %, the kernel no longer builds the analysis); solid 6 eigenfrequencies 2D 0.81 against 0.82; 3D eigensolve 41.32 against 41.58; `bench` total (process wall) 60.0 against 61.0 (+1.6 %, inside the spread of 60.0 to 63.3 of the control); `cavity_bench 0.012 2` (whole example) 3.89 against 3.82 (-2 %; two pairs after the dropped one). Nothing beyond noise; the solvers are untouched. Logs: `nyx2-phase5-bench` next to the repository.

**Open / for others**
1. The engine has no `analyse` command yet (the plan lists it in phase 7); `PolarizationAnalysis` and `FieldFrame` are what it needs. `nyx-store`'s `metrics_of` still reads `Solution::te_fraction` of the old result (phase 7).
2. `legacy::solve_with_analysis`, `Solution::te_fraction` and `femcore::problem` go together in phase 7.
3. A physics with other vector fields (sms displacement) gets an analysis by stating a frame; none is written.
4. The overall phase table above still says "waits for 4" for phase 5 (the table is not mine to edit).
5. Private build directories to delete when accepted: `nyx2-target-phase5` and `nyx2-target-phase5-ctl` (about 4 GB each), workspaces `nyx2-phase5-ctl` (`jj workspace forget nyx2-phase5-ctl`), logs `nyx2-phase5-*.log` and `nyx2-phase5-bench`.

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

## Phase 7: engine, thin front-ends, MCP parity (workspace `nyx2-phase7`, bookmark `restructure-7-engine-frontends`, 2026-11-05)

Started from the merge of `restructure-5-analysis`, `restructure-6-materials` and `restructure-7a-connect` (conflicts were in the docs tables, `Cargo.lock` and the Connect-menu hunk of `app.rs`; line endings of the old app files were CRLF, new files are LF). State: **implemented, awaiting acceptance.** ADR 57 to 60.

**Done.**
* **Engine.** The mesher takes a `MeshRequest` (2D and 3D, adaptive, structured PML, element sizes from the waves; `ems` and `sms` state their own `Sizing`, `Physics::sizing` takes the library). Named, oriented and project materials resolve in the engine (`Ctx.library`, `Ctx.wavelength`, `MaterialAssign::Oriented`, `Project.materials`, `Library::with_texts`). New commands: `get_guide`, `list_analyses`, `analyse`, `list_materials`, `describe_material`, `add_material`, `remove_material`, `add_pml_conditions`, `set_geometry`, `set_problem`. The default solve includes the physics the active one imports from. `facade` re-exports what a front-end draws with. `auto_store`: the MCP server binds the project to a workflow at the first solve.
* **MCP.** `nyx-mcp` is a loop over `Engine::catalog()`/`call()`; all 23 geometry operations are tools from the descriptors; `nyx-session` is deleted. Contract `mcp_tools.v2.json` (breaking), `engine_commands.v1.json` 1.1. An end-to-end stdio test builds a wire from tools, solves, analyses and finds the run in the stored workflow.
* **Window.** On `Engine`; Study, Physics and Materials tabs drawn from descriptors (`forms.rs`); examples are files (`examples/*.nyx.json`, 16, written through commands by `tests/make_examples.rs`); the app depends on `nyx-engine` and `nyx-connect` only. App tests include a headless frame of every tab on five examples.
* **Deleted.** `nyx_engine::legacy`, `femcore::{problem, usecases, xmlio, h5io}` (and the `hdf5`/`phdf5`/`quick-xml` dependencies), the XML project format, `worker.rs`, built-in material presets, `nyx-store`'s dependency on `femcore` (ADR 59), `complex_bench`, `soi_modes`, `bragg_gratings`, `usecases.rs`, `legacy_*` tests. The layering table has no temporary edge. The benchmark programs now call `Physics::solve` on the generated meshes.

**Checks.** Merge point before any change: 901 pass, 0 fail, 11 ignored. After the phase: see the line at the end of this section. Benchmarks (best of 3, control built fresh from `restructure-5-analysis` in `nyx2-phase7-ctl`, same machine, alternating): `bench` 33.2 s against 34.0 s (after a first rewrite of `bench` that skipped the two eigenfrequency runs was corrected to do the same work as the control); `profile3d 16 16 16 p2` 1.60 s against 1.60 s; `mem_profile 12 12 12 p2 static auto` 0.72 s against 0.74 s; `cavity_bench 0.012 2` 3.19 s against 3.17 s. No difference beyond noise.

**Decisions to confirm (paul).**
1. The one-call `recipe_run` (SOI strip) and `check_run` are gone; the guide shows the command sequence. A generic `mesh_convergence` command is not written. Say if you want one.
2. The Bragg grating examples are dropped (swept layers with `domain_offset` have unnamed domains, so materials cannot be assigned by name); naming them is a mesh follow-up.
3. HDF5 output is not ported (optional engine export, later).
4. Examples use the library materials, so their numbers differ slightly from the old window's rounded ones.
5. Physics-sized meshes in the SOS example are slow (the `13_silicon_on_sapphire` test takes about 60 s); not compared with the old window.

**For phase 8.** Spec files are still named `experiment.nyx.xml` / `specs/*.nyx.xml` and hold JSON; the project JSON has no `format` key; `docs/architecture/workspace.dsl`, `diagrams.md` and `README.md` of the repository still describe the old structure (phase 9). `xmlio` is gone, so there is no XML left to delete.

**Not touched.** The window itself was not seen on screen in this session (the scripted run exited cleanly but produced no screenshot); please look at it once.

**Phase 7 result:** `cargo test --release --workspace --no-fail-fast`: 887 pass, 0 fail, 10 ignored (merge point 901 / 0 / 11; the difference is the deleted use-case, legacy and old-contract tests, the new engine, example, MCP, app and store tests). `check_code_to_adr.sh` and the contract and layering suites pass.

**Follow-up after review (2026-11-05, paul).** Decision 1 answered: convergence is essential, written (`convergence_check`, ADR 61: engine command, MCP tool, Study-tab button). Decision 3 answered: HDF5 stays central; the writer is ported to the new types (`nyx-engine::h5io`, feature `hdf5`, ADR 62): project, mesh and all fields of every physics, exact round trip, 2.5 to 2.9 times smaller than JSON, `results.h5` in every recorded run, commands `export_results` and `results_file_info`. XDMF is written too (recovered from the old writer, `nyx-engine::xdmf`; every reference verified by a test). Open: (a) paul's "fast xml format": is XDMF what was meant, or XML for the project file (the plan says JSON, phase 8)? (b) make `hdf5` a default feature (needs the HDF5 library everywhere; here `HDF5_DIR=C:/Users/robin-adm/anaconda3/Library` and its `bin` on PATH build and run it); (c) the window cannot open a results file yet; (d) the MPI-IO writer is not ported. Builds with the feature use a separate target directory in these notes (`nyx2-target-phase7-h5`).

**Result after the follow-up:** default build 889 pass, 0 fail, 10 ignored (workspace); build with `hdf5`: engine, app and MCP tests 75 pass, 0 fail (incl. `tests/h5.rs`: exact round trip, XDMF references, `results.h5` in every run).

## Phase 8: file format (workspace `nyx2-phase8`, bookmark `restructure-8-format`, 2026-11-12)

Started from `main` = `restructure-7-engine-frontends` (74851193, phase 7 merged by paul on 2026-11-05 to 2026-11-12). The commit is tagged `last-old-format`. State: **implemented, awaiting acceptance.** ADR 63.

**Done.**
* **`format` key (value 1) at three levels:** `Project`, `geom::model::FeatureList`, `nyx_physics::ProblemSpec`. `Project::from_json` is the file boundary and is strict (key required in the project, the geometry and every problem; a different value is refused); a command given a problem or feature list without the key reads it as the current format. No migration code; a pre-phase-8 file fails with a message that says what is missing. The layout itself is the one phase 7 derived from the envelope (G6): no other key was renamed or moved.
* **No XML left that is not output.** `xmlio` and the XML project format were already gone in phase 7. Now: spec files are `experiment.nyx.json` and `specs/<hash>.nyx.json` (were `.nyx.xml`), `strip_paths_xml` is `strip_paths`, the window's open dialog offers JSON, the XML-compatibility check in the geom tests is gone. The only XML written is the XDMF description next to HDF5 results (ADR 62).
* **No bulk data as JSON.** The JSON copy of the fields in a build without `hdf5` (`result_<physics>.json`) is removed. A run then records metrics and log only and `log.txt` says why. With `hdf5` every run keeps `results.h5` and `results.xdmf` as before.
* **Contracts.** `engine_project.v2.json` and `problem_spec.v3.json` (the old files removed; `VERSIONS.md`), the geom golden feature lists and the 16 `examples/*.nyx.json` carry the key. `*.expected.json` and every numeric golden are unchanged. `store_workflow.v1` and `store_run.v1` are unchanged (the store does not know the layout, ADR 59).
* New tests: `crates/nyx-engine/tests/format.rs` (4), one in `geom/tests/selectors.rs`; the spec contract test and the store test of the run folder check the new rules.

**Checks.** `cargo test --release --workspace --no-fail-fast` (default build): 894 pass, 0 fail, 10 ignored (phase 7: 889 / 0 / 10; +5 new tests). With `hdf5` (`HDF5_DIR=C:/Users/robin-adm/anaconda3/Library`, its `bin` on PATH): `nyx-engine`, `nyx-mcp`, `nyx-store` all pass, including `tests/h5.rs`. Contract and layering suites and `check_code_to_adr.sh` pass. No solver code changed, so no benchmark was re-run. `clippy -D warnings` on the new crates stops at an unused import in `nyx-ems/tests/physics_meshing.rs` and `cargo fmt --check` at old formatting in `nyx-engine/examples`; both exist on `main` and were left alone (phase 9, R23).

**Open for paul (carried from phase 7, ADR 62/63).**
1. ~~Make `hdf5` a default feature?~~ Answered yes (2026-11-12) and done; see "Decisions by paul".
2. ~~The window cannot open a results file yet.~~ Done 2026-11-12 (`open_results`).
3. ~~The MPI-IO writer is not ported.~~ Ported 2026-11-12 behind feature `phdf5`; needs an MPI-enabled HDF5 library and an MPI SDK to build, which this machine lacks: please build `-p nyx-engine --features phdf5` once on a machine that has them.
4. "Fast xml format": is XDMF what was meant? The project file is JSON by plan; XML for it would be a reversal.
5. Projects written before this phase (the old examples, a saved workspace of the phase 7 build) do not open; they are regenerated or re-saved by hand.

**Decisions by paul (2026-11-12).** `hdf5` is a default feature of `nyx-engine` (done: `Cargo.toml`, `.cargo/config.toml` sets `HDF5_DIR`; the DLL must be on PATH). The window opens a results file (**done**: `open_results`, Open dialog and drop; 901 pass, 0 fail, 10 ignored). The MPI-IO writer is ported (**done, but the MPI part is not compiled or run here**: no MPI SDK or parallel HDF5 on this machine; the partition logic is tested with three emulated processes; 904 pass, 0 fail, 10 ignored). XDMF is as far as XML goes; the project file stays JSON.

## Phase 9: close-out (workspace `nyx2-phase9`, bookmark `restructure-9-closeout`, 2026-11-13)

Started from `main` = `restructure-8-format` (0e2a5e12, phase 8 accepted and merged by paul on 2026-11-12). State: **implemented, awaiting acceptance.** ADR 64, 65.

**Done.**
* **Lints and formatting (R23, ADR 65).** `cargo clippy --workspace --all-targets -- -D warnings` is clean (it was about 75 warnings: renamed egui calls, dead code in the app, auto-fixes, `type_complexity` and `large_enum_variant` allowed with a reason, three new lints of this toolchain). `rustfmt` is applied to `geom`, `nyx-physics`, `nyx-dfs`, `nyx-sms`, `nyx-ems`, `nyx-analysis`, `nyx-engine` (`cargo fmt --check` passes on them) and to no other crate. `rust-version` is 1.99.
* **Guard tests.** `femcore/tests/expr_grammar.rs` (R06: the two expression grammars accept and refuse the same corpus and agree on every value) and `nyx-engine/tests/layering.rs::geom_depends_on_serde_and_optional_rayon_only` (R18).
* **ADR 64** (R08): `hash_bytes` (BLAKE3) is the only persisted hash, all other hashes are cache keys. **ADR 65**: ADR 4, 22, 23 superseded; 5, 14, 15 marked as implemented or amended; the formatting and lint rule; the records of the earlier design phase are left as history (R14); the parked items.
* **Model and docs.** `workspace.dsl` and `diagrams.md` rewritten for the final structure (crate layering exactly as the layering test; engine and femcore components; two dynamic views); `CODE_TO_ADR.md` lost the 16 rows of deleted files; `ROUGH_SPOTS.md` has the resolutions; `README.md` of `docs/architecture`, `ADR_INDEX.md` and `CLAUDE.md` are up to date.

**Not done (needs paul).**
1. ~~`workspace.dsl` is not validated with Structurizr.~~ Done 2026-11-13 (paul agreed to the download): the Structurizr CLI (about 99 MB zip, in `opomech/structurizr-cli`, outside the repository) validates it with exit 0 and exports all six views to Mermaid; a deliberately broken copy is refused with the line number. No slip was found. `diagrams.md` stays hand-written.
2. The licence check of the BSD-2-Clause ports against nyx's own licence (never done; parked in the plan "until the end"): nothing may be published before it.
3. Run-time plug-ins, two-way coupling and Brillouin as an optional physics are parked as in the plan.
4. `-p nyx-engine --features phdf5` (the MPI-IO writer, phase 8 follow-up) has still not been built on a machine with MPI.

**Checks.** `cargo test --release --workspace --no-fail-fast`: 907 pass, 0 fail, 10 ignored (phase 8: 904; the three new tests are the two expression-grammar tests and the geom dependency test). `cargo clippy --workspace --all-targets -- -D warnings` clean; `cargo fmt --check` clean for the seven formatted crates; `check_code_to_adr.sh` passes. No solver numerics changed (lint edits only), so no benchmark was re-run.

## Phase 9c: the full C4 set (workspace `nyx2-c4`, bookmark `restructure-9b-c4`, 2026-11-13)

Started from `main` = `restructure-9a-structurizr` (a6219af1). State: **implemented, awaiting acceptance.** No ADR (documentation and tooling only).

**Done.** `docs/architecture/c4/`: context, containers (one per crate), components (91, every container, a view each; femcore also in three parts) and code pages (one per component: files, uses, used by, class diagram, table of public items), two dynamic views, exports (Mermaid, PlantUML, JSON), a browser viewer, a memo and a README. The master `docs/architecture/workspace.dsl` keeps the hand-written model and includes the generated components, relationships and views. `c4/tools/gen.pl` reads the Rust source and `c4/spec/components.spec`; `build.sh` validates with the Structurizr CLI and exports; `make_standalone.sh` builds the standalone copy `Documents/opomech/nyx-c4/` (outside the repository). Checks: the DSL validates (exit 0); all 110 Mermaid diagrams render in headless Chrome; every `.rs` file of `src/` belongs to one component; the crate-to-crate references found in the source are all in the layering table. Rebuild with `bash docs/architecture/c4/tools/build.sh` when files or crates change.

**Not done.** The worked 1550 nm use case with a dynamic view (next step); a C4 deployment view; tests, examples and benches are outside the model.

## Capability R6, workstream A: axis conditions of the axisymmetric solid (workspace `nyx2-sms-axis`, bookmark `capability-r6-sms-axis`, 2026-11-14)

Started from `main` = `restructure-9b-c4` (152eb4ed). Plan: `SMS_AXISYM_PLAN.md` steps 1 and 2 (K1, gate G3). State: **implemented, awaiting merge with workstreams B (K2 cylindrical PML) and C (K3 fields and the study).** The other plan steps are not in this branch.

**Done.** A mesh may touch the axis `r = 0`. `axisym::axis_constraints` holds the components regularity forces to zero on the axis dofs (`m = 0`: `u_r`, `u_theta`; `|m| = 1`: `u_z`; `|m| >= 2`: all three) and `assemble_cross_section` adds them to the fixed set; a cell with a vertex at `r < -1e-9` of the mesh extent is refused (the old per-point `r <= 0` check stays and only fires on a zero-area cell, because the quadrature rules are interior: confirmed in `femcore::fem::quad_rule`). The log says how many axis dofs were held and which components. Derivation in the header of `axisym.rs`, decision in ADR 67 part 1, gap recorded as `ROUGH_SPOTS.md` R6a.

**Shared signature touched.** `CrossAssembly` gained one public field, `axis_dofs: usize` (built in one place, `assemble_cross_section`). No function signature changed. `cross_section_fields` gained only a log line.

**Checks.** Baseline `cargo test --release --workspace --no-fail-fast` before the tests below existed: 907 pass, 10 ignored, 0 fail (it ran while the K1 edit was being written, so treat it as the documented baseline reproduced). `cargo test --release -p nyx-sms`: all green with the new `tests/axisym_axis.rs` (6 tests); `ring_model_axes` golden, `axisym_closed_form.rs`, `ring_frame.rs`, `stages.rs`, `contracts.rs` pass untouched. `clippy -D warnings` and `cargo fmt -p nyx-sms --check` clean; `check_code_to_adr.sh` passes. Measured errors (P2): T1 torsion of a solid cylinder 1.3e-5 and 9.9e-5; T2 radial modes with `u_z = 0` on the ends 6.0e-8 and 1.2e-5; T3 thin disk (t/R = 0.02) radial contour modes 1.5e-5 and 7.7e-5; `m = 2` lowest frequency changes 1.1e-4 under one refinement; `m = 1` free cylinder has two rigid modes at 1e-13 and 1.3e-4 Hz against an elastic 2282 Hz. Thresholds are 5e-4, 1e-4, 5e-4, 1e-2.

**Open points.** R6a (the `m = +-1` tie is a penalty, not imposed). The ADR numbers were re-assigned to 68 and 69 in workstream D, which also did the docs step.

## Capability R6, workstream D: the merge, the engine and MCP proof, the docs (workspace `nyx2-sms-axisym`, bookmark `capability-r6-sms-axisym`, 2026-11-14)

Started from the three workstream bookmarks on `main` = `restructure-9b-c4` (152eb4ed): A `capability-r6-sms-axis` (862042c7, K1 axis), B `capability-r6-sms-pml` (66761af8, K2 cylindrical PML), C `capability-r6-sms-study` (aeccff2c, K3 and the study). Plan: `SMS_AXISYM_PLAN.md` steps 6, 7 and 9 (step 8, the A1 run through MCP, is a later thread). State: **implemented, awaiting acceptance.** ADR 68, 69.

**Done.**
* **The merge** (an octopus merge of the three bookmarks). Conflicts: `axisym.rs` (header, A and B), `cross_section.rs` (`CrossAssembly` and its constructor, A, B and C: the struct keeps both `axis_dofs` and `pml_cells`), and the ADR 67 file that A and B had both written (add/add; joined into one ADR). Nothing else conflicted. No golden changed; `ring_model_axes` passes untouched.
* **The one contract consequence:** the `physics_sms_result` exemplar entry `axisymmetric_eigenfrequency_pml` changed only in its mode label (Q 2.751 to 2.475 of the one-ring example), because C wrote it with the old PML. Rewritten with `NYX_WRITE_CONTRACT=1` inside the existing 1.1 line, which says so. `mcp_tools` and `engine_commands` did not move.
* **ADR numbers.** `restructure-10-a3-snowflake` already has ADR 66 (sweep layers) and 67 (pmc wall) and nothing uses 68 or higher, so ours are **68 (the study, was 66)** and **69 (axis conditions and cylindrical stretch, was 67)**. Renumbered in the sources, ADRs, contracts, rough spots and C4 spec. The A3 line will conflict with this one on `ADR_INDEX.md`, `CODE_TO_ADR.md`, `ROUGH_SPOTS.md` (it adds R29 to R31; ours are R6a to R6e) and the status file, not on ADR numbers.
* **Tests held to off-axis meshes or the old PML, strengthened:** `axisym_pml.rs` T5b and T6b repeat T5 and T6 on a solid disk on a post on a substrate that all touch the axis (A1's topology); `axisym_study.rs` solves a thin solid disk on the axis through the envelope, the typed problem and the library entry (three-way equal) and against the contour closed form (2.4e-5). Numbers are in ADR 69.
* **Step 6.** `examples/17_disk_anchor_loss.nyx.json` (built by commands, `make_examples.rs`; the other 16 files are byte-identical on regeneration), T10 (`examples.rs`: solves, every mode has an eigenvalue and a layer fraction, the contour mode at 4.505e8 Hz has Q = 1733 and a layer fraction of 0.001, layer modes have Q near 2), T11 (`convergence_check` rows for the real and imaginary parts), T12 (`nyx-mcp/tests/stdio.rs`: `describe_physics sms` lists the study with its three parameters, the tool sequence solves and `result_info` carries `Q1`, `layer fraction 1`, `f1 (Hz)`).
* **T13** `crates/nyx-engine/tests/benchmark_a1.rs`: `#[ignore]`, builds A1's model by commands with the paper's dimensions; `CRITERIA_AGREED = false`, so the criteria are printed and not asserted until paul agrees them. Indicative run with the assumed material (E 160 GPa, nu 0.22, rho 2330) and a substrate 30 um wide and deep, layers 7.4 um: f = 712.87 MHz (0.4 % from 715.6), Q = 5883 (6 % from 6250), layer fraction 0.000; Q moves by more than 5 % between h = 0.5 and 0.25 um at the disk, and by less than 5 % when the layers double. **Not a result**: the material is assumed and the mesh criterion is not met.
* **Step 7 docs.** ADR 68 and 69 finalised, ADR 43 status, `ADR_INDEX.md`, `CODE_TO_ADR.md` (check passes), C4 spec and `workspace.dsl` (container text and the AICommand description; the generated dynamic README is rewritten by the generator, so the anchor-loss run is named in the view description instead), `build.sh` and `make_standalone.sh` run (`nyx-c4` is updated), `CLAUDE.md`, `HOW_IT_WORKS.md` (new paragraph), `nyx-instructions` guide card (one step, 2647 characters), `ROUGH_SPOTS.md` (R6b to R6e), `VERSIONS.md`. `CAPABILITY_PLAN.md` rows (uncommitted file of the default workspace) are marked done.

**Checks.** `cargo test --release --workspace --no-fail-fast`: **932 pass, 0 fail, 11 ignored** (baseline 907 / 10 / 0; the extra ignored one is `benchmark_a1`). `cargo clippy --release --workspace --all-targets --no-deps -- -D warnings`: clean. `cargo fmt --check` for the seven formatted crates: clean. Contract and layering suites pass; `bash docs/architecture/check_code_to_adr.sh`: 0 files without a row. Build folder `../nyx2-sms-axisym-target`.

**Open points for paul.**
1. A1 material constants and pass criteria (plan section 10 item 3): `CRITERIA_AGREED` stays false until answered; step 8 (the MCP run, `benchmarks/` folder, `benchmarks/README.md`) is a later thread. Nothing in `benchmarks/` was committed.
2. The cylindrical stretch r to r~ is kept (the default). Measured difference to the paper's untransformed form: about 1 % in Q (B: 0.2 to 1.3 % on an annular wheel; D: 0.03 to 1.0 % on a disk on the axis).
3. Q of an A1-type model is mesh-sensitive at coarse meshes (R6e): 4.6 % for doubled layers on the T5b model, 9 % for a 0.7 refinement of example 17.
4. Follow-ups, not done: low-reflecting boundary with `r` weight (R6b), exact `m = +-1` tie (R6a), complex shift (R6c), `s_theta` at the quadrature points (R6d).
5. The ADR numbers (68, 69) and the merge with the A3 line (the files named above).
