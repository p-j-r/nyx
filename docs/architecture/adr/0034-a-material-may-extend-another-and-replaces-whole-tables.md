# 34. A material may extend another one and replaces whole tables

Date: 2026-11-04

## Status

Accepted (paul's question of 2026-11-04: does adding measured values need a whole new material or only replace some values of an existing one? Only replace). Additive to material file format 1 (schema 1.1): files without `extends` read as before. Builds on ADR 29.

## Context

A researcher with an own ellipsometry table for silicon should not have to copy `silicon.toml` and edit it, which would also copy the stiffness and thermal data and lose the link to the base.

## Decision

* A material file may contain `extends = "<name>"`. The base is a material of the same `Library`, found by file name (ignoring case); `Library::load_dir` resolves the chain between the files of a folder, `Library::read(text, name)` resolves one text (for a file outside the folder, for example one embedded in a project) against the loaded library.
* **Merge rule** (`material::merge_over`): every top-level key of the file replaces the base's; a table (`[elastic]`, `[thermal]`, `[electrical]`, `[optical]`) is replaced as a whole and never mixed key by key. So a half-replaced stiffness (the constants of one class with the matrix of another), a unit key that follows a value it was not written for (`density_unit`), or a Sellmeier term mixed into a table cannot occur. `name` and `extends` are not inherited; the name is the file's own or its file name.
* **Validation runs on the merged result** with the same code as for any file: symmetry, stability, range, no extrapolation, units, source. An error in a derived file says that the combined result was checked and names the base.
* **Errors**: a missing base names the key and suggests the closest library names; a loop is reported as its chain (`a -> b -> c -> a`); a base that itself fails is reported with its message; reading a file with `extends` without a library is refused with a pointer to `Library`.
* Overriding one axis or one coefficient inside a table is deliberately not possible: restate the table. This keeps the rule explainable in one sentence; field-level merging can be added later if it proves necessary.

## Consequences

* Schema 1.1 adds `extends` (additive). The resolved `MaterialSpec` carries no trace of `extends`; `to_toml` writes the combined material.
* Example `materials/examples/silicon_measured.toml`; tests in `crates/femcore/tests/material_extends.rs`.

## Code

`crates/femcore/src/material/spec.rs` (`Library`, `resolve`, `merge_over`), `materials/examples/`, `crates/femcore/tests/material_extends.rs`, `docs/MATERIALS.md`.

## Contract

`material_schema.v1.json` (version 1.1).
