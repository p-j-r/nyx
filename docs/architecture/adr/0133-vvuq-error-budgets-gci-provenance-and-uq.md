# 133. VVUQ Phases 3–5: Error budgets, 3-grid GCI, validation, determinism, provenance, and `vvuq_report.v1`

Date: 2026-11-26

## Status

Accepted (VVUQ Phases 3–5 / Layers 3–5; renumbered from provisional `0127` in `docs/vvuq/VVUQ_VERIFICATION_PLAN_EGM_REFINED.md` because ADRs `0125–0130` were already allocated, Phase 0–1 took `0131`, and Phase 2 took `0132`). Extends ADR 19, ADR 53, ADR 61, ADR 62, ADR 91, ADR 99, ADR 107, ADR 120, ADR 131, and ADR 132.

## Context

Layers 3–5 close the end-to-end verification and reproducibility loop on top of Layers 1–2:
* **Layer 3 (Solution Verification & Error Budgets):** Every `SolveResult` needs algebraic residual floor, conditional 1-norm `rcond`, solver convergence warnings, spatial discretization error (`Prager-Synge`, 2-grid delta, or 3-grid Roache GCI), existing AMR wiring (`femcore::adapt::maxwell`, `femcore::adapt::estimate3`), and unconditional surfacing through HDF5 worker round-trips, `result_info`, MCP, and GUI.
* **Layer 4 (Physical Validation & Lightweight UQ):** Published experimental and cross-code validation suites in `benchmarks/validation/` enforced via the explicit ASME V&V 20 `AND` rule keeping Prager–Synge epistemic bounds strictly separate from aleatoric $1\sigma$ uncertainties, plus active-study material provenance (`Source::flags_uncertainty()`, `T_range`, `TempCurve::eval_checked`) and zero-solve cQED / LHS-Sobol QMC `sweep_run` UQ.
* **Layer 5 (Determinism, Provenance & CI):** Two-tier (`<= 60k` 0-ULP vs `> 60k` bounded spectral/residual) determinism, 0-dep `nyx-store` `BuildProvenance`/`RuntimeProvenance`, deterministic HDF5 provenance without `created_utc`, `/provenance/materials_json` drift checks, unsigned BLAKE3-linked `vvuq_report.v1.json`, and 3-tier CI matrices.

## Decision

* **Unified diagnostics & error budgets (`solve_result` `2.0 -> 2.1`):** Extend `nyx_physics::SolveResult` additively with `diagnostics: Option<SolverDiagnostics>` and `error_budget: Option<NumericalErrorBudget>` (`#[serde(default, skip_serializing_if = "Option::is_none")]`), plus `DiagnosticSeverity`, `DiagnosticWarning`, and `SpatialErrorEstimate`, persisted under `/results/<physics>/{diagnostics_json,error_budget_json}` in HDF5 layout v4.
* **Roache GCI & `BoundCertified` (`femcore::vvuq`, `nyx-engine::{convergence,validated}`):** Add `roache_gci` supporting worker-isolated opt-in `grids: 3` and inverted coarsening ladders, `R11a` structured `mode_index` metadata, `R11b` low-$r < 1.3$ fallback warnings, backward-compatible `"state"` with additive `"tier"` elevating Prager–Synge solves to `bound_certified`, and pragmatic `adapt::maxwell` (`ems::cavity`) + `adapt::estimate3` (`sms` `R6f`) wiring.
* **Validation, material validity & lightweight UQ:** Enforce ASME V&V 20 `AND` checks in `benchmarks/validation/` and 17/17 `examples.rs`, surface active-study material/temperature warnings without `material_schema` v2.1 bloat, and provide zero-solve cQED UQ + streaming LHS/Sobol `SweepCheckpoint`.
* **Two-tier determinism, provenance & `vvuq_export`:** Enforce subprocess thread determinism (`<= 60k` bit-exact / `> 60k` `< 1e-10`), add `nyx-store` provenance (`store_run` `1.0 -> 1.1`), `nyx-engine/build.rs`, deterministic HDF5 root provenance attrs, and `vvuq_export` emitting `vvuq_report.v1.json`.

## Consequences

* Workspace layering (`crates/nyx-engine/tests/layering.rs`) stays unchanged (`nyx-store` and `nyx-storage` keep zero workspace dependencies).
* Legacy HDF5 results files and existing JSON contracts remain backward compatible via `#[serde(default)]` optional fields.

## Code

`crates/femcore/src/vvuq.rs`, `crates/nyx-physics/src/{result,lib}.rs`, `crates/femcore/src/{linalg,adapt,material}/*`, `crates/nyx-{dfs,sms,ems}/src/*`, `crates/nyx-analysis/src/*`, `crates/nyx-store/src/runs.rs`, `crates/nyx-engine/{build.rs,src/*}`, `crates/app/src/{app,tabs}.rs`, `benchmarks/validation/*`.
