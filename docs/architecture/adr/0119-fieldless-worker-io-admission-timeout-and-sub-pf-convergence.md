# 119. Field-less worker HDF5 I/O, per-mode read budgets, memory admission timeouts, and sub-picofarad linear convergence checks

Date: 2026-11-22

## Status

Accepted. Amends ADR 61, 62, 86, 107, and 112.

## Context

Large 3D multi-mode solves (`100k+` tetrahedra, 10+ modes) suffered from three operational bottlenecks: (1) writing and round-tripping multi-gigabyte field arrays from isolated workers when callers only needed derived cQED analyses, while whole-dataset read caps rejected valid multi-mode files whose individual modes fit comfortably in memory; (2) unbounded memory-admission polling in headless batch/script runs when RAM remained low; and (3) `convergence_check` clamping `< 1e-12` capacitance values to a zero relative-error floor and squaring `adapt_scale` when refining meshes.

## Decision

* **Per-mode HDF5 budget & field-less `save_fields: false` + `post_analyses` (`nyx-storage`, `nyx-engine`):**
  * Enforced per-mode row-slice read budgets (`shape[1..] * elem_bytes <= max_dataset_bytes`, default `512 MiB`) alongside an overall `max_total_result_bytes` cap (`8 GiB` default, `NYX_H5_MAX_BYTES` env override) in `h5io::check_packed_budget`.
  * Added `save_fields: bool` (default `true`) and `post_analyses: Vec<Value>` across `QueueOpts`, `solve_run`, `queue_submit`, `sweep_run`, `worker::serve_job`, `shared_queue`, `h5io::write_results_with_options`, and `Stored::fields_saved` (`@fields_saved = 0` HDF5 attribute). When `save_fields == false`, `Engine::run_post_analyses` executes requested analyses on full in-memory fields before stripping mesh/field datasets, skips `xdmf::write_xdmf`, and keeps `100k+`-tet `results.h5` files under `256 KiB`.
* **Bounded memory admission timeout (`nyx-engine::estimate`, `app::runfile`):**
  * Added `memory_wait_timeout_s: Option<f64>` (`>= 0.0`) on `solve_run`, `queue_submit`, `sweep_run`, and `Engine::set_memory_wait_timeout`, obeying priority `explicit memory_wait_timeout_s > NYX_MEMORY_WAIT_TIMEOUT_S > caller default` (`30 s` for batch CLI `nyx run` / `run_script`; `None` on bare `Engine::new()`). Exceeding the wait aborts promptly with an error naming `need`, `free`, `total`, and `timeout_s`.
* **Sub-picofarad relative convergence, `quantities` filtering, and linear `adapt_scale` (`nyx-engine::convergence`):**
  * Replaced the `scale < 1e-12` relative-error floor with unit-fraction absolute comparison and a `scale <= 1e-30` zero guard so sub-pF (`fF`–`aF`) capacitances compare true relative errors.
  * Added validated `quantities: Option<Vec<String>>` filtering to `convergence_check`.
  * Stopped multiplying `adapt_scale` by `factor` in `convergence::refined` (scaling `points_per_wavelength /= factor` and `dz_min *= factor` instead) so `h`, `Op::Size`, `dz_min`, and `wave_size` refine by `factor^1` instead of `factor^2`.

## Consequences

* Isolated worker sweeps and multi-mode cQED jobs avoid multi-gigabyte HDF5 round-trips, fail fast under memory starvation in unattended batches, and perform accurate linear mesh-convergence checks on femtofarad-scale superconducting capacitances.

## Code

`crates/nyx-storage/src/lib.rs`, `crates/nyx-engine/src/{h5io,xdmf,store_link,estimate,jobs,worker,queue,shared_queue,engine,sweep,commands,convergence}.rs`, `crates/app/src/runfile.rs`, `crates/nyx-engine/tests/{h5,worker,admission,convergence_sizes,squadds_fr_e2e}.rs`, `docs/architecture/contracts/engine_commands.v1.json`.
