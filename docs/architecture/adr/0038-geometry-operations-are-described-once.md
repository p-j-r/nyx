# 38. Geometry operations are described once, and every tool is generated from the description

Date: 2026-11-04

## Status

Accepted (gate G5, paul 2026-11-04: the MCP exposes every tool geom offers; strip, rib and ring are the acceptance test, not the limit).

## Context

The MCP had one hand-written template (`soi_strip`). Gate G5 says nothing may be left out and nothing hand-picked. `geom` has 23 operations (`Op`) with serde forms, but no description a tool or a form can be generated from.

## Decision

* `nyx_engine::geomdesc` holds, for each `Op` kind, a descriptor: title, what it does, group (shape, build, labels, boundary), dimensions, and its fields with type, required flag and default. From it: a JSON Schema per operation, the tool name (`add_rect`, `add_rounded_rect`, ...), `build_op` (tool arguments to `Op`: features by **name**, selectors as **text**, defaults filled, unknown or missing fields reported with the valid names) and `args_of` (an `Op` back to the same form).
* The descriptors live in the engine, **not in `geom`**, because `geom` depends on serde only (ADR 3, pinned by `crates/geom/tests/contracts.rs`) and the descriptors need `serde_json`.
* Completeness is enforced: `kind_of` is an exhaustive match on `Op`, so adding a variant does not compile until it is described; a test builds every kind from its descriptor and checks that the descriptor's field names equal the serde field names; a test checks that the engine has an `add_*` command per descriptor; `build_op` goes through the same serde reader that loads a saved feature list.
* Beyond the `add_*` tools the engine offers the edits of a feature list as commands (`update_feature`, `remove_feature`, `rename_feature`, `suppress_feature`, `tag_feature`, `set_param`, `remove_param`) and `geometry_names` (regions and boundary names before any mesh exists). Every edit is a transaction: the list is validated after the edit and left untouched if the edit is refused.
* The steps of a `Transform` and the `kind` of a `Pattern` are written in the serde form of `TStep` and `PatternKind` (documented in the schema); they are the only nested values not flattened to text.

## Consequences

* A new `geom` operation costs one descriptor and shows up as a tool, in the catalog contract and as a source for the window's forms.
* The acceptance test (rebuild a strip, a rib and a ring through commands, evaluate and mesh them) is `crates/nyx-engine/tests/engine.rs`.

## Code

`crates/nyx-engine/src/{geomdesc,geometry,commands}.rs`, `crates/nyx-engine/tests/geomdesc.rs`.

## Contract

`geom_operations.v1.json`.
