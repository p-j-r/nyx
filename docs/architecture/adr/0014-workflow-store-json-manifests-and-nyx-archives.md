# 14. Work is stored as workflows: JSON manifests, an append-only run log and deterministic `.nyx` archives

Date: 2026-11-01

## Status

Accepted by implementation (`93321f5` nyx-store, `3dd760f` design and mockups, `37e5871` wired into the app). The open questions of `docs/WORKFLOWS_DESIGN.md` section 11 (should runs of different specs mix in one experiment?) are open.

## Context

Researchers iterate: change a width, re-run, compare. Autosaved project files in `nyx_runs` did not say which work belonged together, could not be searched, and could not be shared as one unit.

## Decision

* The unit of work is a **workflow**: a folder holding `workflow.json` (manifest: experiments, tags, status), `history.jsonl` (one line per version of the set-up), `notes.md`, `specs/` and one folder per experiment. Workflows live under a workspace folder (`workspace.json`, default `Documents/nyx`).
* **JSON, not XML,** for manifests and runs: faster, appendable, trivial for MCP. The experiment's **spec stays the project XML** (`.nyx.xml`), so existing projects and their round-trip tests keep working.
* Runs are **append-only** lines in `runs.jsonl` (`RunRecord`: id, time, seconds, status, spec hash, version of the set-up, parameters, metrics, tags, star, note, artifacts, error); a digest is derived from them. Metrics come from `SolveResult` (`runs::metrics_of`): n_eff, first scalar, TE fraction, dof count, every numeric scalar in snake_case.
* A `.nyx` file is **one deterministic small ZIP** with an `ArchiveManifest` (format, kind, level, nyx version, files); import is a safe extraction (no path escapes). Three share levels: recipe (default), results, meshes.
* The **concept card** (family, dimensions, physics, study, wavelength, materials inferred from refractive index) is generated from the spec, never typed; search uses plain words with a synonym table and numbers with units.
* Everything is plain serde with ordered maps so files are byte-for-byte reproducible. blake3 hashes identify specs.
* Nothing is deleted by rules; "Done" leaves the dashboard at once; a 12-card cap and drafts that never touch disk.

## Consequences

* `nyx-store` knows `ProblemSpec` and `Project` only through untyped values read tolerantly, so a newer spec still opens, and a changed layout is silently unseen (R25 is exactly this).
* `RunRecord.femcore` is filled from the wrong crate's version (R21).
* The store depends on all of `femcore` for two items (R15).
* `workflow.json` carries its format under the key `nyx` (`FORMAT = 1`) and `WorkflowDir::open` refuses a newer one with "made with a newer nyx"; the archive manifest has a `format` too. Its `version` key is the number of the current set-up version, not a format number. No migration code exists.

## Code

`crates/nyx-store/src/*`, `crates/app/src/{shell,wfview,sheets}.rs`, `docs/WORKFLOWS_DESIGN.md`. History: `93321f5`, `3dd760f`, `959d1d5`, `37e5871`, `2b5b66a`, `b5e5bfb`.

## Contract

`contracts/store` (`crates/nyx-store/tests/contracts.rs`): exemplar `workflow.json` and `runs.jsonl` line must round-trip, an old run line without the later keys must still read, a newer manifest must be refused, the metric names of a `SolveResult` are pinned, and the concept-card contract is written and `#[ignore]`d pending R25. The archive manifest is not covered yet.
