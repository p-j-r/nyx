# 131. VVUQ Phase 0 P0 soundness, security, and Layer 1 formal verification gates

Date: 2026-11-25

## Status

Accepted (VVUQ Phases 0–1; renumbered from provisional `0125` in `docs/vvuq/VVUQ_VERIFICATION_PLAN_EGM_REFINED.md` because ADRs `0125–0130` are already allocated on `feat/nyx-roadmap`). Extends ADR 10, ADR 19, ADR 26, ADR 61, ADR 91, ADR 112, and ADR 120.

## Context

`docs/vvuq/VVUQ_VERIFICATION_PLAN.md` and `docs/vvuq/VVUQ_VERIFICATION_PLAN_EGM_REFINED.md` define a 5-layer verification, validation, and uncertainty quantification (VVUQ) roadmap for `nyx`. Before numerical code verification (Layer 2) or solution error budgets (Layer 3) can be trusted, Phase 0 (`P0-01`..`P0-15`) and Phase 1 (Layer 1 Rust formal methods & 3-tier CI) close verified `unsafe` UB risks, silent `NaN`/degenerate-mesh masking, unvalidated untrusted inputs (GDSII, HDF5, `.nyx` archives, MCP framing, Monty script timeouts, HTTPS material fetch SSRF), and non-deterministic `Problem` keys or subprocess thread propagation.

## Decision

* **Phase 0 (`P0-01`..`P0-15`):**
  * Replace unsafe `RawCols::at` and unvalidated `LineSmoother` indexing with lifetime-bound tile slicing (`RawCols<'a, T>::chunk_mut`) and construction-time index/slice validation.
  * Enforce release-mode bounds and `u32::try_from` checks across `Triplets` and `Csr` constructors (`femcore::linalg::csr`).
  * Extend `nyx_physics::validate::check_result_finite` across `fields`, `cell_fields` (`re`/`im`), `eigenvalue`, and `value` while permitting valid lossless `Q = f64::INFINITY` on result scalars.
  * Fix `curve3` zero/NaN determinant masking, implement `SimplexMesh::check()` + `from_cells_split_with` tet orientation + `mesher3d` post-weld `det3 > 0`, remove `has_sheets` gating on `adapt/ops3.rs` overlap checks, and default `strict_mesh = true` across `GeomMeshOpts`, `MeshParams`, `MeshRequest`, and `MeshSettings` (with `default_true`/`is_true` serde helpers). `SimplexMesh::check()` and `extrude_with` tet orientation also uncovered and closed two latent bugs: `pml_geom.rs` failing to clear pre-existing facets before repopulating 2D `r = 0` boundaries, and `nyx-ems` `Ned3` tie-breaking flipping edge signs on 1-layer periodic extrusions (`Ned3::key_pairs` / `basis3_keyed`).
  * Enforce GDSII hierarchy/record/polygon resource caps (`geom::gds`), HDF5 `max_solutions` & byte-budget guards (`nyx-storage`, `nyx-engine::h5io`), deterministic HDF5 root attributes (no `created_utc`), `.nys`/`.xdmf` archive exclusions (`nyx-store::archive`), `nyx-script` timeout/`printer!` output caps, `nyx-mcp` bounded framing/`sync_channel`/`-32002`, subprocess `NYX_THREADS` + TOML-hashed `Engine::problem_key` + `set_project` check reset, `driven_prom` assertion, contract README / `CODE_TO_ADR` / portable `gates.rs`, and HTTPS allowlist/SSRF/Windows-reserved-name hardening in `material_fetch.rs` + script-sandbox refusal.
* **Phase 1 (Layer 1 Formal Methods & 3-Tier CI):**
  * Stand up `.github/workflows/vvuq-ci.yml`, `tools/vvuq_ci.sh`, and `tools/vvuq_ci.ps1` implementing Tier 1/2/3 verification gates.
  * Lower `#[cfg(miri)]` parallel dispatch thresholds (`femcore::linalg::{par, amg_elastic}`), add `#[cfg(kani)]` harnesses `K1–K3` (`femcore`, `nyx-storage`, `geom::gds`), and add `crates/femcore/tests/portable_kernels.rs` verifying `< 1e-13` differential parity between default SIMD and `--features portable-kernels`.
  * Deploy 5 unified `bolero` fuzz/regression targets (`fuzz_gds`, `fuzz_nyx_storage`, `fuzz_script_trace` with explicit `rounds > 500` error in `trace::guard`, `fuzz_path_and_archive`, `fuzz_mcp_framing` over `nyx-mcp::lib`), expand `proptest` across `geom` (`props`, `bvh`, fixing a thick-sheet lateral-seam bug in `topo3::build`), `femcore` (`vvuq_props`), and `nyx-physics` (`vvuq_props`), keeping `bolero`/`proptest` strictly in `[dev-dependencies]`.
  * Achieve `finite_guard` parity across `dfs`/`sms`/`ems` (`femcore::validate::non_finite_token`), add `SolverError::NotFactored` and typed `Result` returns in `femcore::linalg::direct` + `nyx-physics::decode::Params` + `nyx-ems` + `nyx-analysis`, and complete `// SAFETY:` + scope-guard audits on `app::console` and `nyx-install::platform::windows`.

## Consequences

* `mesh_request` (`1.4 -> 1.5`), `engine_commands` (`1.31 -> 1.32`), and `script_api` (`1.2 -> 1.3`) evolve additively while maintaining backward compatibility for saved projects and HDF5 files.
* All 17 `examples/*.nyx.json` pass under `strict_mesh = true`.

## Code

`docs/vvuq/*`, `crates/{geom,femcore,nyx-physics,nyx-storage,nyx-store,nyx-engine,nyx-script,nyx-mcp,nyx-ems,nyx-install}/*`, `.github/workflows/*`, `tools/vvuq_ci.*`.
