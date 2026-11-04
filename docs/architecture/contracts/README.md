# Interface contracts

The rule is ADR 19. Each contract is an exemplar file here plus a test that generates the same bytes from the code.

| Contract | Pins | Test | Run |
|---|---|---|---|
| `problem_spec` | serde shape of `ProblemSpec` (one spec per physics), XML encoding, leniency, `Physics` order | `crates/femcore/tests/contracts.rs` `problem_spec_*`, `physics_enum_*` | `cargo test --release -p femcore --test contracts` |
| `solve_result` | structure of a `SolveResult`, order refusal, `"cancelled"` | same file, `solve_*`, `a_cancelled_solve_*` | same |
| `frames` | wafer-frame constants, name tables, tensor permutations, PML axes | same file, `frames_*` | same |
| `mesh` | label bases, node order of quads, no curved-node storage, `.fmsh` version and round trip | same file, `mesh_*` | same |
| `kernel_stages` | public stage functions of every kernel family (R26 closed in phase 3) | same file, `kernel_stage_api_is_public`; composition in `tests/stages.rs` | same |
| `geom_model` | geom depends on serde only; golden feature lists round-trip key for key | `crates/geom/tests/contracts.rs`; `tests/goldens.rs` | `cargo test --release -p geom --test contracts` |
| `store` | `workflow.json`, `runs.jsonl` line, format gate, metric names, concept card | `crates/nyx-store/tests/contracts.rs` | `cargo test --release -p nyx-store --test contracts` (`-- --include-ignored` shows R25) |
| `mcp_tools` | tool names, annotations, arguments, handshake, handle format | `crates/nyx-mcp/tests/contracts.rs` | `cargo test --release -p nyx-mcp --test contracts` |
| `problem_spec` v2 | serde shape of the problem envelope (`nyx-physics`), one spec per physics | `crates/nyx-physics/tests/contracts.rs` `problem_spec_v2_*` | `cargo test --release -p nyx-physics --test contracts` |
| `solve_result` v2 | serde shape of a result | same, `solve_result_v2_*` | same |
| `physics_descriptors`, `physics_api` | descriptor types; trait surface, registry rules, error codes | same | same |
| `physics_dfs`, `physics_sms`, `physics_ems` | the descriptors each physics declares (studies, conditions, properties, fields, couplings, generated JSON Schema) | `crates/nyx-{dfs,sms,ems}/tests/contracts.rs` | `cargo test --release -p nyx-dfs -p nyx-sms -p nyx-ems --test contracts` |
| `layering` (tests only) | which workspace crate may depend on which | `crates/nyx-engine/tests/layering.rs` | `cargo test --release -p nyx-engine --test layering` |
| `linalg_options`, driver | `EigenOptions` shape and defaults; the driver surface the kernels call; closed-form direct solve, Laplacian and non-symmetric eigenvalues, dense helpers, cancel | `crates/femcore/tests/linalg_contract.rs` | `cargo test --release -p femcore --test linalg_contract` |

Rewriting an exemplar: `NYX_WRITE_CONTRACT=1 cargo test ... --test contracts <test name>`. Do it only with a `VERSIONS.md` entry. None of these tests depends on frozen floating-point values, so they pass on a machine where the golden files fail (R22).

Files: `problem_spec.v1.json` (old `ProblemSpec`, until phase 4), `problem_spec.v2.json`, `solve_result.v2.json`, `physics_descriptors.v1.json`, `physics_{dfs,sms,ems}.v1.json`, `store_workflow.v1.json`, `store_run.v1.json`, `mcp_tools.v1.json`; versions in `VERSIONS.md`. Exemplars are compared as values, never as text: the order of the keys of a JSON object is not part of a contract (ADR 26).
