# Contract versions

One line per contract. Bump the version here in the same commit that changes the exemplar (ADR 19). Major = breaking, minor = additive.

| Contract | Exemplar / snapshot | Version | Changed | Why |
|---|---|---|---|---|
| problem_spec | `problem_spec.v1.json` | 1.1 | 2026-11-03 | `Physics::Brillouin`, `BrillouinSpec` and the `brillouin` key removed (phase 1, ADR 24). Deliberate rewrite; no reader of the old file is kept. Superseded by `problem_spec.v2.json` in phase 4 |
| store_workflow | `store_workflow.v1.json` | 1.0 | 2026-11-03 | first snapshot; the format number inside the file is `nyx = 1` |
| store_run | `store_run.v1.json` | 1.0 | 2026-11-03 | first snapshot |
| mcp_tools | `mcp_tools.v1.json` | 1.0 | 2026-11-03 | first snapshot, 9 tools |
| geom_model | `crates/geom/golden/*.json` | 1.0 | 2026-11-03 | no exemplar of its own; the goldens are the exemplars; the file carries no format key yet |
| mesh | (tests only) `.fmsh` version byte 2 | 1.0 | 2026-11-03 | label bases, node order, no curved nodes |
| frames | (tests only) `section` constants | 1.0 | 2026-11-03 | wafer frame, +y |
| kernel_stages | (tests only) public stage functions | 1.1 | 2026-11-03 | the list of one-call solvers loses `sbs::solve_brillouin` (ADR 24); stage functions unchanged |
| solve_result | (tests only) structure of `SolveResult` | 1.0 | 2026-11-03 | not serialisable (R27) |
| problem_spec | `problem_spec.v2.json` | 2.0 | 2026-11-04 | the envelope `{ physics, order, study{kind,params}, materials, bcs[{on,kind,params}], params }` of `nyx-physics` (ADR 25). Coexists with v1 until phase 4 deletes the old `ProblemSpec` |
| solve_result | `solve_result.v2.json` | 2.0 | 2026-11-04 | `SolveResult`, `Solution`, `Field`, `CellField`, `Scalar` of `nyx-physics` are serialisable (R27); no `analysis` member |
| physics_descriptors | `physics_descriptors.v1.json` | 1.0 | 2026-11-04 | serde shape of the descriptor types (`StudyDesc`, `BcDesc`, `PropDesc`, `ParamDesc`, `FieldDesc`, `FieldLink`, `Sizing`) |
| physics_dfs | `physics_dfs.v1.json` | 1.0 | 2026-11-04 | first descriptors of `dfs`, read from the current kernels; the phase 4 owner refines them (bump here) |
| physics_sms | `physics_sms.v1.json` | 1.0 | 2026-11-04 | same, `sms` |
| physics_ems | `physics_ems.v1.json` | 1.0 | 2026-11-04 | same, `ems` |
| physics_api | (tests only) the `Physics` trait surface, registry rules, error codes | 1.0 | 2026-11-04 | `crates/nyx-physics/tests/contracts.rs`; changing a method is a bump |
| material_schema | `material_schema.v1.json` | 1.0 | 2026-11-04 | JSON Schema of material file format 1 (`femcore::material::schema`, ADR 29, 30); test `crates/femcore/tests/material_contract.rs` checks it lists exactly the keys the loader accepts |
| material_example | `material_kitchen_sink.v1.toml` + `.json` | 1.0 | 2026-11-04 | one file using every section of format 1 and the resolved `MaterialSpec` it must read to (SI values, Voigt matrix in Pa) |
| material_schema | `material_schema.v1.json` | 1.1 | 2026-11-04 | additive: `extends` (ADR 34), `source` as text, table or list at the top and in every table (ADR 35), models `schott`, `herzberger`, `drude`, `lorentz`, `drude_lorentz` with `frequency_unit`, `plasma_frequency`, `damping`, `drude_weight`, `eps_inf`, `oscillators`. Files of schema 1.0 stay valid. Exemplar rewritten on purpose |
| material_example | `material_kitchen_sink.v1.json` | 1.1 | 2026-11-04 | the resolved material now carries the sources; `source` is a list of tables instead of strings (`MaterialSpec.source: Vec<Source>`). Rewritten on purpose |
| material_example | `material_kitchen_sink_models.v1.toml` + `.json` | 1.0 | 2026-11-04 | second file: the Schott, Herzberger and Drude-Lorentz models on three axes, table-form source |

The file name carries the major version (`.v1.`). A breaking change adds `.v2.` next to `.v1.`, keeps the old file for the reader of the previous major, and updates the test.
