# 19. Interfaces that outlive a run are pinned by exemplar files and contract tests, and carry a version

Date: 2026-11-03

## Status

Accepted by this change (the tests and exemplars are in the same commit). Adding a `format` key to `Project` and `FeatureList` is **proposed, not done** (it alters every saved file).

## Context

Saved projects, `specs/` files, `workflow.json`, `runs.jsonl`, geom feature lists, `.fmsh` files and the MCP tool list are read by code that was not built together with the code that wrote them: later versions of nyx, nyx-store reading specs as untyped values, AI hosts that registered `nyx` weeks ago. serde reads leniently (unknown keys ignored, `#[serde(default)]` keys filled in), so a renamed or dropped field does not fail; it silently loses a setting (ROUGH_SPOTS R20). Only `workflow.json`, `.fmsh` and the HDF5 file had a version gate.

## Decision

1. Each interface has a **contract**: a checked-in exemplar under `docs/architecture/contracts/<name>.v<N>.json` (or a snapshot built from the running program for the MCP tool list) and a test that produces the same bytes from the code (`crates/femcore/tests/contracts.rs`, `crates/geom/tests/contracts.rs`, `crates/nyx-store/tests/contracts.rs`, `ai/nyx-mcp/tests/contracts.rs`). The exemplar is the schema; a separate JSON Schema would need a validator dependency and a second thing to keep in step.
2. A contract test fails on **any** change of the interface, additive or not. The failure message points here. Changing the exemplar (`NYX_WRITE_CONTRACT=1`, same convention as the goldens in ADR 10) is a deliberate act that goes together with an entry in `docs/architecture/contracts/VERSIONS.md`:
   * **additive** (new optional key with a default, new tool, new optional argument): bump the minor, keep reading old files;
   * **breaking** (rename, removal, changed meaning or unit, new required key, changed node order or label base): bump the major, add a reader for the previous major or state that old files are refused, and say so in the commit.
3. Behaviour that the interfaces *rely on* but that is not a data shape is pinned too: the label bases of meshes, the node order of quadrilaterals, the `.fmsh` version byte, the wafer-frame constants of `section`, the string `"cancelled"`, the set of public kernel-stage functions.
4. Contract tests do **not** use frozen floating-point results; those are the goldens of ADR 10 and currently fail on this machine (R22). Contract tests must pass on any platform.
5. A contract that exposes a known defect is written anyway and `#[ignore]`d with the rough-spot number in the reason (`store_card_reads_current_project`, R25). `cargo test -- --include-ignored` shows it failing.
6. Stable identity rule, recorded once: spec and file identity inside nyx-store is blake3 over the bytes; AI-session handles are FNV-1a over the inputs; geom caching uses geom's own content signature. They are not interchangeable.

## Consequences

* A reviewer sees an interface change as a diff of one JSON file plus a line in `VERSIONS.md`.
* The exemplars are not a JSON Schema; they do not say which keys are optional. The tests that exercise leniency (`problem_spec_reading_is_lenient_in_both_directions`, `store_run_record_exemplar_is_stable_and_old_lines_still_read`) document what is optional today.
* Proposed next steps, each its own change: a `format` key on `Project` and `FeatureList` read by a migration function; a contract crate holding `ProblemSpec`, `SolveResult`, `Project` and `xmlio` (R15, R16, R27); a contract for the archive manifest and for the `.nyx.xml` encoding of a `Project`.

## Code

`crates/femcore/tests/contracts.rs`, `crates/geom/tests/contracts.rs`, `crates/nyx-store/tests/contracts.rs`, `ai/nyx-mcp/tests/contracts.rs`, `docs/architecture/contracts/*`.

## Contract

This ADR is the rule; `docs/architecture/contracts/README.md` lists the contracts and how to run them.
