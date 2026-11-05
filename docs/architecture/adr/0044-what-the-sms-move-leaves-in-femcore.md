# 44. What the sms move leaves in femcore

Date: 2026-11-04

## Status

Accepted (restructure phase 4, bookmark `restructure-4-sms`). Records the seams between the sms and ems moves, which ran in parallel.

## Context

The solid kernels shared a few things with the electromagnetic ones, and the two moves were done by separate sessions that must merge without either rewriting the other's files.

## Decision

Stays in `femcore`, used by `nyx-sms` as it is:

* `periodic` (Floquet reduction of dofs: `match_vertices`, `nodal_relations`, `Reduction`, `congruence`): shared with `ems3d`; neither move touches it.
* `ems_pml::cell_stretches` (the complex stretch factors of a perfectly matched layer per cell): used by `sms::pml` and `sms::cross_section`. It is public in `femcore` and **must stay public there or move to a shared place** (`props` is the candidate) when `ems_pml` moves into `nyx-ems`; `nyx-sms` must then change one import. `section::mesh_pml_dirs` it relies on stays.
* `section`, `material`, `fields`, `fem`, `optical`, `polarization`, `analysis`, `props`: as ADR 41 says.

Small changes in shared files, each a few lines, listed so the merge can be checked:

* `progress::lap_since(&mut Instant, label)` replaces `sms::lap`, which `fem.rs`, `ems3d.rs` and `linalg/amg_elastic` called for stage timings (six lines in `ems3d.rs`, three in `fem.rs`, seven in `amg_elastic/mod.rs`). If the ems move deletes `ems3d.rs` or edits those lines, take its version and call `progress::lap_since`.
* `validate::non_finite_token` is `pub` (the typed problem's finite check scans tensors and layers for the words NaN and inf).
* `solve.rs` answers an error for `Physics::Solid` as it does for the dfs arms; `lib.rs` loses the six `sms*` modules. `femcore::problem` and the old `ProblemSpec` stay for the window until phase 7.
* `h5io`'s round-trip test builds a solid-shaped result by hand (the kernel is no longer in femcore), like its heat one.

## Consequences

* The three crates `nyx-dfs`, `nyx-sms`, `nyx-ems` depend on `femcore` only; the layering test is unchanged by this move.
* Until the ems move lands, `femcore::solve` still runs `ems`.

## Code

`crates/femcore/src/{progress,validate,solve,lib,h5io}.rs`, `crates/nyx-sms/src/pml.rs`.
