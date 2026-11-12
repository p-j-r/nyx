# 63. One `format` key, JSON for small files, HDF5 for bulk data

Date: 2026-11-12

## Status

Accepted by the plan (phase 8) and by paul's storage decision of 2026-11-05 (HDF5 is central; project and spec files stay small JSON; no bulk data as JSON). Replaces the provisional layout of the project file and the `engine_project` exemplar (ADR 37 named it provisional); completes ADR 62.

## Context

Until phase 8 the saved layout was provisional: no file said which layout it was. Spec files were still called `experiment.nyx.xml` and `specs/<hash>.nyx.xml` although they have held JSON since phase 7 (the XML project format is gone with `xmlio`), `nyx-store` had a function `strip_paths_xml` that parsed JSON, and a build without HDF5 wrote the fields of a result as `result_<physics>.json`.

## Decision

* **`format` key, first value 1, written once at three levels:** `Project` (`nyx-engine`, `PROJECT_FORMAT`), `FeatureList` (`geom::model::FORMAT`) and `ProblemSpec` (`nyx_physics::PROBLEM_FORMAT`). The three are separate because each is also read alone: a feature list by `set_geometry` and the geom goldens, a problem by `set_problem` and the physics contracts. A value other than 1 is refused with a plain message by `FeatureList::validate`, `validate_spec` and `Project::check_formats`.
* **No migration code.** `Project::from_json` (the file boundary) requires the key in the project, in its geometry and in every problem; a file of the earlier layout is refused with a message that names the missing key and says older files are not converted. A command that is given a problem or feature list *without* the key (an MCP call, a test) reads it as the current format; only files are strict. The last commit that wrote the old layout is tagged `last-old-format` (74851193).
* **Spec files are JSON by name.** `experiment.nyx.json`, `specs/<hash>.nyx.json`, `strip_paths` (was `strip_paths_xml`). The `nyx` archive accepts exactly these names. There is no XML code left in the workspace except the XDMF *description* file written next to HDF5 results (ADR 62), which is XML by the XDMF specification and is output only. The project and spec files are not XML.
* **No bulk data as JSON.** `result_<physics>.json` is gone. A recorded run keeps `results.h5` (+ `results.xdmf`) in a build with feature `hdf5`; a build without it records the metrics and the log only, and says so in `log.txt`. Run records, `workflow.json`, `history.jsonl`, specs, material files and descriptors stay small text. Golden files of tests are inputs and oracles, not stored results.
* **Contract exemplars.** `engine_project.v2.json` (was `.v1.`, provisional) and `problem_spec.v3.json` (was `.v2.`) carry the key; the old files are removed (`VERSIONS.md`). `store_workflow` and `store_run` do not change: the store holds spec text and a hash and does not know the layout (ADR 59). The golden feature lists of `geom` and the 16 example projects were rewritten with the key and nothing else; the numbers in `*.expected.json` are untouched.

## Consequences

* A later change of layout bumps the number, adds a reader or a refusal in one place and a `.v<n+1>.` exemplar.
* Project files written before phase 8 (the examples of phase 7, a saved workspace) do not open. They are regenerated (`make_examples`) or re-saved from a window built before; there is no converter.
* **`hdf5` is a default feature of `nyx-engine`** (paul, 2026-11-12), so the app, the MCP server and every test build have it; `.cargo/config.toml` sets `HDF5_DIR` for the development machine (an environment value wins) and the HDF5 DLL must be on PATH at run time. `--no-default-features --features dfs,sms,ems` builds without it, and a run then stores no fields.
* Decided by paul on 2026-11-12, to do: the window opens a results file; the MPI-IO writer is ported (for GPU and multicore use); XDMF stays the XML side, the project file stays JSON.

## Code

`crates/geom/src/model.rs`, `crates/nyx-physics/src/{spec,validate}.rs`, `crates/nyx-engine/src/{project,store_link}.rs`, `crates/nyx-engine/tests/format.rs`, `crates/nyx-store/src/{flow,archive}.rs`, `docs/architecture/contracts/{engine_project.v2,problem_spec.v3}.json`.
