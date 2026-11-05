# 35. Provenance is first-class; Drude, Lorentz, Schott and Herzberger models; the library is checked against the literature

Date: 2026-11-04

## Status

Accepted (requests of paul relayed on 2026-11-04). Additive to material file format 1 (schema 1.1). Extends ADR 29 and 30.

## Decision

* **Provenance.** The top of a file and each table (`[elastic]`, `[thermal]`, `[electrical]`, `[optical]`, each optical axis) may carry `source`: a citation text, a table `{ cite, url, doi, verified }`, or a list of those (`material::source`). `verified` is the date `YYYY-MM-DD` on which the values were compared with the publication; `url` must be http(s), `doi` must start with `10.`. The loader keeps the sources in `MaterialSpec`, `Elastic`, `Thermal`, `Electrical`, `Optical` and `AxisModel`, and `to_toml` writes them back.
* **Library rule, enforced by a test** (`material_library.rs`): every shipped file has a top-level source and every table that holds data has a source (an optical axis may share the one of `[optical]`); a source without a verification date must say "approximate" or "unverified" in its text.
* **More dispersion models**: `schott` (`n^2 = a0 + a1 x^2 + a2 x^-2 + ...`), `herzberger` (its 0.028 um^2 constant is converted to the declared wavelength unit), and the oscillator family `drude`, `lorentz`, `drude_lorentz` in the `e^{j omega t}` convention: `eps = eps_inf - f0 Wp^2/(w^2 - j g0 w) + sum f_k Wp^2/(w0_k^2 - w^2 + j g_k w)`, frequencies in a declared `frequency_unit` (`eV`, `meV`, `THz`, `GHz`, `Hz`, `rad/s`, `1/cm`; `Hz` units are converted with `2 pi`). This is the form of Rakic et al. 1998. Same rules as the other formula models: declared wavelength unit, required `range`, error outside the range, per axis, derivative, thermo-optic shift.
* **Library checked against sources** (2026-11-04). Corrected: silicon stiffness to Hopcroft 2010 (c11 165.6, c12 63.9, c44 79.5 GPa; the code's `silicon()` still carries 165.7/63.9/79.6), silicon thermal conductivity 130 to 148 W/(m K) (flagged approximate: sources range 130 to 156), silicon dn/dT 1.8e-4 to 1.87e-4 1/K (Frey 2006), silica to Corning 7980 data (E 72.7 GPa, nu 0.16, density 2201, specific heat 770, expansion 0.52e-6; dn/dT 9.7e-6 from the Corning table), sapphire parallel conductivity 25.2 to 25.8, parallel expansion 5.3e-6 to 5.6e-6, parallel permittivity 11.6 to 11.5, specific heat 760 to 756 (supplier data, flagged approximate; suppliers disagree). Confirmed unchanged: sapphire stiffness (Wachtman 1960), the Sellmeier data of silica, sapphire (both rays) and the silicon table. Added: gold, silver, copper, aluminium (Rakic-LD, reproduced against the published tables to 2e-4), silicon nitride (Luke 2015), lithium niobate (Zelmon 1997). Not verifiable in this session and therefore flagged: air properties, silica permittivity (a range), silicon permittivity (11.7 or 11.9).

## Consequences

* The in-code `sapphire()` and `silicon()` differ slightly from the library; the pin test in `materials.rs` allows that and says why. They go when the engine loads the library (phase 7).
* Metals other than the four, temperature-dependent Sellmeier models and dn/dT for sapphire are not in the library.
* The site ioffe.ru (a common source for semiconductor data) has an invalid certificate; its values were not used and the certificate check was not bypassed.

## Code

`crates/femcore/src/material/{source,optics,units}.rs`, `materials/*.toml`, `crates/femcore/tests/{material_library,materials}.rs`.

## Contract

`material_schema.v1.json` (1.1), `material_kitchen_sink_models.v1.toml` and `.json`.
