# Interface contracts

The rule is ADR 19. Each contract is an exemplar file here plus a test that generates the same bytes from the code.

| Contract | Pins | Test | Run |
|---|---|---|---|
| `problem_spec` | serde shape of `ProblemSpec` (one spec per physics), XML encoding, leniency, `Physics` order | `crates/femcore/tests/contracts.rs` `problem_spec_*`, `physics_enum_*` | `cargo test --release -p femcore --test contracts` |
| `solve_result` | structure of a `SolveResult`, order refusal, `"cancelled"` | same file, `solve_*`, `a_cancelled_solve_*` | same |
| `frames` | wafer-frame constants, name tables, tensor permutations, PML axes | same file, `frames_*` | same |
| `mesh` | label bases, node order of quads, no curved-node storage, `.fmsh` version and round trip | same file, `mesh_*` | same |
| `kernel_stages` | public stage functions of 7 kernel families | same file, `kernel_stage_api_is_public`; composition in `tests/stages.rs` | same |
| `geom_model` | geom depends on serde only; golden feature lists round-trip key for key | `crates/geom/tests/contracts.rs`; `tests/goldens.rs` | `cargo test --release -p geom --test contracts` |
| `store` | `workflow.json`, `runs.jsonl` line, format gate, metric names, concept card | `crates/nyx-store/tests/contracts.rs` | `cargo test --release -p nyx-store --test contracts` (`-- --include-ignored` shows R25) |
| `mcp_tools` | tool names, annotations, arguments, handshake, handle format | `ai/nyx-mcp/tests/contracts.rs` | `cd ai && cargo test --release -p nyx-mcp --test contracts` |

Rewriting an exemplar: `NYX_WRITE_CONTRACT=1 cargo test ... --test contracts <test name>`. Do it only with a `VERSIONS.md` entry. None of these tests depends on frozen floating-point values, so they pass on a machine where the golden files fail (R22).

Files: `problem_spec.v1.json`, `store_workflow.v1.json`, `store_run.v1.json`, `mcp_tools.v1.json`; versions in `VERSIONS.md`.
