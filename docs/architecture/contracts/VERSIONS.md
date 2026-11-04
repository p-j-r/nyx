# Contract versions

One line per contract. Bump the version here in the same commit that changes the exemplar (ADR 19). Major = breaking, minor = additive.

| Contract | Exemplar / snapshot | Version | Changed | Why |
|---|---|---|---|---|
| problem_spec | `problem_spec.v1.json` | 1.2 | 2026-11-04 | `Physics::Brillouin`, `BrillouinSpec` and the `brillouin` key removed (phase 1, ADR 24); 1.2: the `out_of_plane` key removed (sms_kz model, ADR 24 amendment). Deliberate rewrite; no reader of the old file is kept. Superseded by `problem_spec.v2.json` in phase 4 |
| store_workflow | `store_workflow.v1.json` | 1.0 | 2026-11-03 | first snapshot; the format number inside the file is `nyx = 1` |
| store_run | `store_run.v1.json` | 1.0 | 2026-11-03 | first snapshot |
| mcp_tools | `mcp_tools.v1.json` | 1.0 | 2026-11-03 | first snapshot, 9 tools |
| geom_model | `crates/geom/golden/*.json` | 1.0 | 2026-11-03 | no exemplar of its own; the goldens are the exemplars; the file carries no format key yet |
| mesh | (tests only) `.fmsh` version byte 2 | 1.2 | 2026-11-04 | boundary-element numbers of the nyx `CellMesh` are 1-based like the domain numbers (ADR 32 step 1); `SimplexMesh` labels unchanged, `from_cells` no longer shifts. 1.2 (ADR 40): additive, a 3D geom model meshes through `mesh_model` with the label-map numbering (test `mesh_geom3d_numbers_cells_and_facets_with_the_model_label_map`) |
| mesh_cells | `mesh_cells.v1.json` | 1.0 | 2026-11-04 | the one mesh model (ADR 32): node order per cell type, faces, label base 1 for domains and boundaries, `.fmsh` version 2, how the solver view is cut from it; the nyx `CellMesh` it replaces is gone |
| frames | (tests only) `section` constants | 1.0 | 2026-11-03 | wafer frame, +y |
| kernel_stages | (tests only) public stage functions | 1.1 | 2026-11-03 | the list of one-call solvers loses `sbs::solve_brillouin` (ADR 24); stage functions unchanged |
| solve_result | (tests only) structure of `SolveResult` | 1.0 | 2026-11-03 | not serialisable (R27) |
| problem_spec | `problem_spec.v2.json` | 2.0 | 2026-11-04 | the envelope `{ physics, order, study{kind,params}, materials, bcs[{on,kind,params}], params }` of `nyx-physics` (ADR 25). Coexists with v1 until phase 4 deletes the old `ProblemSpec` |
| solve_result | `solve_result.v2.json` | 2.0 | 2026-11-04 | `SolveResult`, `Solution`, `Field`, `CellField`, `Scalar` of `nyx-physics` are serialisable (R27); no `analysis` member |
| physics_descriptors | `physics_descriptors.v1.json` | 1.0 | 2026-11-04 | serde shape of the descriptor types (`StudyDesc`, `BcDesc`, `PropDesc`, `ParamDesc`, `FieldDesc`, `FieldLink`, `Sizing`) |
| physics_dfs | `physics_dfs.v1.json` | 1.0 | 2026-11-04 | first descriptors of `dfs`, read from the current kernels; the phase 4 owner refines them (bump here) |
| physics_sms | `physics_sms.v1.json` | 1.1 | 2026-11-04 | first descriptors of `sms`; 1.1: the `kz` study parameter removed with the out-of-plane model |
| physics_ems | `physics_ems.v1.json` | 1.0 | 2026-11-04 | same, `ems` |
| physics_api | (tests only) the `Physics` trait surface, registry rules, error codes | 1.0 | 2026-11-04 | `crates/nyx-physics/tests/contracts.rs`; changing a method is a bump |
| engine_commands | `engine_commands.v1.json` | 1.0 | 2026-11-04 | the catalog of `nyx-engine` commands: name, title, description, access, JSON Schema of the arguments (what the MCP server turns into tools); includes one `add_*` command per geometry operation |
| geom_operations | `geom_operations.v1.json` | 1.0 | 2026-11-04 | descriptors of the 23 geom operations (`nyx_engine::geomdesc`): fields, types, defaults, generated schema (ADR 38) |
| engine_project | `engine_project.v1.json` | 0.1 | 2026-11-04 | an example project written by commands alone. **Provisional**: the final layout is derived in phase 8 (gate G6); no `format` key yet |

The file name carries the major version (`.v1.`). A breaking change adds `.v2.` next to `.v1.`, keeps the old file for the reader of the previous major, and updates the test.
| linalg_options | `linalg_options.v1.json` | 1.1 | 2026-11-04 | 1.1 adds `par_solve_min_n` (default 0, serde default so 1.0 files still read); 1.0 was the first snapshot: `EigenOptions { tol, block, ncv }` and defaults of `femcore::solvers::driver` (ADR 33); API and closed-form checks in `crates/femcore/tests/linalg_contract.rs` |
| kernel_stages | (tests only) public stage functions | 1.2 | 2026-11-04 | the stage functions of plain `sms`, `sms_lrb` and `ems3d` are public (R26, phase 3); additive |
