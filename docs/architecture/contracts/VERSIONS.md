# Contract versions

One line per contract. Bump the version here in the same commit that changes the exemplar (ADR 19). Major = breaking, minor = additive.

| Contract | Exemplar / snapshot | Version | Changed | Why |
|---|---|---|---|---|
| problem_spec | `problem_spec.v1.json` | 1.0 | 2026-11-03 | first snapshot, taken from main `a047506` |
| store_workflow | `store_workflow.v1.json` | 1.0 | 2026-11-03 | first snapshot; the format number inside the file is `nyx = 1` |
| store_run | `store_run.v1.json` | 1.0 | 2026-11-03 | first snapshot |
| mcp_tools | `mcp_tools.v1.json` | 1.0 | 2026-11-03 | first snapshot, 9 tools |
| geom_model | `crates/geom/golden/*.json` | 1.0 | 2026-11-03 | no exemplar of its own; the goldens are the exemplars; the file carries no format key yet |
| mesh | (tests only) `.fmsh` version byte 2 | 1.0 | 2026-11-03 | label bases, node order, no curved nodes |
| frames | (tests only) `section` constants | 1.0 | 2026-11-03 | wafer frame, +y |
| kernel_stages | (tests only) public stage functions | 1.0 | 2026-11-03 | 7 of 10 kernel families public (R26) |
| solve_result | (tests only) structure of `SolveResult` | 1.0 | 2026-11-03 | not serialisable (R27) |

The file name carries the major version (`.v1.`). A breaking change adds `.v2.` next to `.v1.`, keeps the old file for the reader of the previous major, and updates the test.
