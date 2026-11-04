# 39. The store takes spec text (JSON or XML) and a card function from its caller

Date: 2026-11-04

## Status

Accepted. Partly implements the phase 7 item "`nyx-store` loses its `femcore` dependency" (R21, R25): the card is no longer made by the store for new projects. The dependency itself stays until the app stops writing XML.

## Context

`WorkflowDir` parsed every spec as `nyx_project` XML through `femcore::xmlio` and made the search card from the old `Project` shape (`active: "Emw"`, `specs[].study`). The engine's project is JSON with a different layout, and a store that knows what a project looks like is the coupling R21 and R25 describe.

## Decision

* A spec text that starts with `{` is JSON, anything else is the old XML. Version diffs, freezing, hashing and the share-time path stripping (`strip_paths_xml`) accept both.
* `WorkflowDir` has a `card_of: fn(&Value) -> Card`; `create_with` and `open_with` set it, `create`, `create_as` and `open` keep the old default. The engine passes `store_link::card_of`, so the card of a new project is made by the code that knows the project.
* Run metrics are made by the engine too (`store_link::metrics_of`) and passed in the `RunRecord`; `runs::metrics_of` (old `SolveResult`) stays for the app until phase 7 rewires it.
* The file names are not changed yet: a JSON spec is stored as `experiment.nyx.xml` and `specs/<hash>.nyx.xml`. Renaming is a format decision (phase 8, gate G6) and touches the `.nyx` archive whitelist and the contract `store_workflow`.

## Consequences

* The app and the old MCP session keep working unchanged; the dashboard can list an engine workflow with its generated card.
* Phase 7 proper deletes the XML branch, `card::from_value`, `runs::metrics_of` and the `femcore` dependency; phase 8 renames the files.

## Code

`crates/nyx-store/src/{flow,archive}.rs`, `crates/nyx-engine/src/store_link.rs`.
