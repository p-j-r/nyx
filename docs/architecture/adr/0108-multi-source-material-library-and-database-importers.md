# 108. Multi-source material library and on-demand database importers

Date: 2026-11-21

## Status

Accepted.

## Context

`nyx` previously shipped ten room-temperature material TOML files (`air`, `aluminium`, `copper`, `gold`, `lithium_niobate`, `sapphire`, `silica`, `silicon`, `silicon_nitride`, `silver`). Real photonics, phononics/MEMS, cryogenic, and superconducting device simulations routinely depend on author-, fabrication-process-, and temperature-specific literature curves (for example Li 1980 vs. Salzberg 1957 vs. Palik for silicon optics; Luke 2015 LPCVD vs. Philipp 1973 stoichiometric vs. PECVD silicon nitride; McSkimin 1953 elastic stiffness; 4 K cryogenic thermal conductivity and expansion; and normal-state / kinetic-inductance thin films). Meanwhile, bundling tens of thousands of uncalibrated 0 K DFT or raw dispersion records (such as the full `refractiveindex.info`, `AFLOW`, `Materials Project`, or `NumBAT` catalogs) inside the repository would bloat distribution size and invite silent misuse without verification.

## Decision

* **Multi-source curated TOML library (`materials/*.toml`, 38 materials):**
  * Preserve all existing base material names and pinned numerical properties so existing models and regression tests remain unchanged, while upgrading base files with cited `[elastic]`, `[thermal]`, `[electrical]` sections and a one-line comparative `description`.
  * Present conflicting authors, deposition processes, temperatures, and thin-film superconducting variants as distinct TOML files (e.g. `silicon_salzberg`, `silicon_palik`, `silicon_mcskimin`, `silicon_4k`, `silicon_nitride_philipp`, `silicon_nitride_pecvd`, `gallium_arsenide_papatryfonos`, `gold_johnson`, `silver_johnson`, `lithium_niobate_mgo`, `silica_thermal_box`, `quartz_alpha`, `polysilicon_polymumps`, `aluminum_nitride`, `aluminum_nitride_piezomumps`, `scaln_30`, `diamond`, `silica_4k`, `sapphire_4k`, `lithium_niobate_4k`, `copper_4k`, and `_sc`-suffixed normal-state superconducting thin films `aluminium_sc_emft`, `tantalum_alpha_sc`, `niobium_sc`, `titanium_nitride_sc`, `niobium_titanium_nitride_sc`, `granular_aluminium_sc`), each carrying a one-line `description` explaining how it differs from sibling variants.
* **Schema & uncertainty hygiene:**
  * Keep material file `format = 1` unchanged (`MaterialContract`). Approximate thin-film/handbook parameters carry `"APPROXIMATE"`, and automatic database imports carry `"UNVERIFIED"` in `[source].cite` (checked via `Source::flags_uncertainty()`).
* **Pure-Rust on-demand importers (`femcore::material::import`, `nyx-engine::material_fetch`):**
  * Zero external crate dependencies added (respecting offline builds). Pure-Rust parsers convert `refractiveindex.info` YAML (Formulas 1, 2, 3, 4, 5, 6 and tabulated `(n, k)`), `AFLOW` AFLUX JSON (`ael_stiffness_tensor`, `agl_thermal_conductivity_300K`, `agl_thermal_expansion_300K`, `density`), `Materials Project` summary/elasticity JSON (`ieee_format` + `dielectric.e_total`), and `NumBAT` JSONC (preserving Brillouin photoelastic `p_ij` / acoustic viscosity `eta_ij` in `notes` while stripping placeholder DOIs).
  * Expose 24 curated `REMOTE_CATALOG` presets (`si_green_2008`, `si_franta_10k`, `aflow_3c_sic`, `numbat_si_poulton`, `nb_weaver`, etc.) without bundling them in `materials/`.
* **Engine, CLI, and GUI surfaces:**
  * Engine commands `list_remote_materials` (`Read`), `import_material` (`Edit`), and `fetch_material` (`Expensive`), with `folder` guarded by `path_guard::PATH_ARGS` and skipped during static `nyx check` planning.
  * `list_materials` and `describe_material` expose each material's comparative `description`.
  * CLI `nyx material {list,show,import,fetch}` in `crates/app/src/matcmd.rs` using `curl`/`curl.exe` over HTTPS (5 MiB cap, 20 s timeout) without adding HTTP crate dependencies.
  * Window Materials tab displays the material `description` and provides a collapsed "Download or import material" panel without introducing `Color32::*` literals.

## Not in this decision

* Complex superconductive Mattis–Bardeen / kinetic-inductance constitutive surface impedance models remain outside `format = 1` here (`_sc` files provide normal-state RF/optical/mechanical/thermal baselines).
* Materials Project authenticated REST calls (`X-API-KEY`) are not performed directly by `nyx`: users download MP JSON externally and convert it locally with `nyx material import mp` or `import_material`.

## Consequences

* Users can compare multiple published models directly inside the bundled 38-material library and pull additional curves from `refractiveindex.info`, `AFLOW`, `Materials Project`, or `NumBAT` on demand into a project or `materials/` folder.
* All converted files validate against `MaterialSpec` (and optional `extends` inheritance) at import time and clearly flag `"UNVERIFIED"` database provenance.

## Code

`materials/*`, `crates/femcore/src/material/import.rs`, `crates/femcore/tests/material_import.rs`, `crates/nyx-engine/src/material_fetch.rs`, `crates/nyx-engine/src/{commands,commandlist,materials}.rs`, `crates/app/src/{matcmd,main,state,app,tabs}.rs`, `docs/MATERIALS.md`.
