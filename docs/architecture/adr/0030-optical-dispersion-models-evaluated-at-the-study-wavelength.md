# 30. Optical properties are dispersion models evaluated at the study wavelength

Date: 2026-11-04

## Status

Accepted (plan section 2, "Materials"; design in `DECISIONS_ROUND_3.md` section 4). Defaults used while paul has not answered the open questions: models constant, Sellmeier, Cauchy and tabulated first; no refractiveindex.info importer.

## Context

Optical materials were one number per domain. Photonic work needs `n(lambda)`, group index, loss and, for heated devices, `dn/dT`, per crystal axis, and must never produce a number outside the range where the data hold.

## Decision

* `[optical]` of a material holds models, not numbers: `constant` (n or eps, optional k), `sellmeier` (`n^2 = A + sum B x^2/(x^2 - C)`, `C` as squared resonance wavelengths or `lambda0` as wavelengths), `cauchy` (`n = sum c_i x^(-2 i)`), `tabulated` (`[wavelength, n]` or `[wavelength, n, k]` rows; monotone cubic (PCHIP) or linear interpolation of n and k against wavelength). A uniaxial crystal has `[optical.ordinary]` and `[optical.extraordinary]` (unique axis along crystal z), a biaxial one `[optical.x]`, `.y`, `.z`.
* **The wavelength unit is declared** (`wavelength_unit = "um"`) and is required whenever a coefficient, a table or a range uses it, so that published coefficients are pasted without conversion.
* **No silent extrapolation.** Formula models need a validity `range`; a table's range is its extent. A wavelength outside the range is an `OpticsError::OutOfRange` whose message names the range. A table may extend beyond its data only if the file says `extrapolate = "linear"`. At load, a Sellmeier resonance inside the range, a non-positive n over the range, an unsorted table, a negative k are refused.
* **Convention**: time dependence `e^{j omega t}`, complex index `n - j k`, relative permittivity `(n - j k)^2`. `Optical::index`, `d_index_d_lambda` (analytic for formulas, spline derivative for tables), `group_index`, `permittivity_diag` and `permittivity_global` (rotated into the wafer frame with the Bunge angles) take the wavelength in metres and an optional temperature.
* **Thermo-optic**: `dn_dT` (1/K) per axis shifts n by `dn_dT (T - T_ref)`, `reference_temperature` default 293.15 K; without a temperature argument nothing shifts.
* The study does not see a formula: the engine evaluates the material at the wavelength (a sweep loops in the engine) and hands the tensor to the physics, which declares the property it needs in its descriptors. Dispersive cavities (a nonlinear eigenproblem) are later work.

## Consequences

* Lorentz, Drude, Schott and Herzberger models are not implemented; `Dispersion` is an enum, so they are additions.
* Closed-form tests: fused silica (Malitson 1965: n(587.56 nm) = 1.458464, n(1.55 um) = 1.4440), sapphire (ordinary 1.7682 and extraordinary 1.7600 at 589.3 nm; 1.7462 and 1.7383 at 1.55 um), silicon (Li, 3.4757 at 1.55 um).

## Code

`crates/femcore/src/material/optics.rs`, `materials/{silicon,silica,sapphire,air}.toml`, `crates/femcore/tests/materials.rs`.

## Contract

As ADR 29 (the schema covers `[optical]`).
