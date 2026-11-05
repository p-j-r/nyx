# 59. The store knows nothing of a project: no femcore, no XML, the engine makes cards and metrics

Date: 2026-11-05

## Status

Accepted. Completes ADR 39 and resolves ROUGH_SPOTS R21 and R25 as the plan asked ("do not fix them in the old structure").

## Context

ADR 39 let the engine pass its own card function to `WorkflowDir`, but the store still parsed `nyx_project` XML through `femcore::xmlio`, made cards from the old `Project` shape and made run metrics from the old `SolveResult`.

## Decision

* `nyx-store` has no dependency on `femcore` or on any workspace crate. A spec is JSON text. `WorkflowDir::create`/`open` use `card::plain`, which reads only a `title`; programs that know projects pass `create_with`/`open_with` with `nyx_engine::store_link::card_of`. Run metrics are made by `nyx_engine::store_link::metrics_of`. `card::from_value`, `card::from_xml`, `runs::metrics_of` and the XML branch of path stripping (`strip_paths_xml`, which keeps its name, handles JSON) are deleted.
* The files keep their names (`experiment.nyx.xml`, `specs/<hash>.nyx.xml`) although they hold JSON; renaming them changes the `.nyx` archive whitelist and is the file-format decision of phase 8 (gate G6).
* No migration: workflows stored with XML specs do not open. This follows constraint 9 of the plan (no migration of old projects).
* `femcore` loses `xmlio` and the `quick-xml` dependency; the HDF5 writer (`h5io`, features `hdf5-metno` and `phdf5`) is deleted with the old result type it was written for.

## Consequences

* The dashboard cards of an engine workflow come from the engine's card (physics, study, wavelength, materials, parameters of the geometry).
* HDF5 results are not written by this build. The plan keeps HDF5 as an optional engine export; porting it to `nyx_physics::SolveResult` and the native mesh is a follow-up.

## Code

`crates/nyx-store/src/{card,flow,archive,runs}.rs`, `crates/nyx-engine/src/store_link.rs`.
