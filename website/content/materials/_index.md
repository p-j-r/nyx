+++
title = "Materials Library"
description = "Unit-explicit format = 2 TOML constitutive models, crystallographic Z-X-Z Euler tensors, dispersion formulas, Mattis–Bardeen superconductors, and remote fetching."
template = "materials.html"
+++

Every material in `nyx` is a plain-text TOML file (`format = 2`, validated against `docs/architecture/contracts/material_schema.v2.json`) in which **every physical quantity carries an explicit unit** (`"496.8 GPa"`, `"2500 kg/m^3"`, `"1395 ueV"`). Missing units, positive-definiteness violations, or wavelengths evaluated outside a formula's declared `range = [min, max]` are rejected immediately with the exact table and key name — silent extrapolation is forbidden.

## Constitutive Sections & Crystallography

- **`[elastic]`:** Isotropic (`E, nu`, `lambda, mu`, or `K, G`), symmetry-class (`isotropic`, `cubic`, `hexagonal`, `trigonal`, `tetragonal`, `orthotropic`, `monoclinic`, `triclinic`), full 6×6 stiffness matrix, or compliance matrix in textbook Voigt order `1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy` with engineering shear strain. Crystal regions are rotated into the wafer frame by Bunge `Z-X-Z` intrinsic Euler angles.
- **`[optical]`:** Supports `constant`, `sellmeier`, `cauchy`, `tabulated` (monotone cubic `pchip` or `linear`), `schott`, `herzberger`, `rii_formula_4`, `rii_formula_5`, `rii_formula_6`, `drude`, `lorentz`, and `drude_lorentz`, plus linear `dn_dT`, cryogenic Bose–Einstein `dn_dT_model`, uniaxial (`[optical.ordinary]` / `[optical.extraordinary]`), and biaxial (`x`/`y`/`z`) dispersion.
- **`[thermal]` & `[electrical]`:** Scalar, diagonal, 6-component, or `{ perp, par }` tensors with continuous temperature curves (`tabulated_T`, NIST `log10_poly`, or `einstein_shift`).
- **Microwave & Quantum (`[dielectric]`, `[superconducting]`, `[interfaces.{ma,ms,sa}]`):** Complex microwave permittivity $\tilde{\varepsilon}_r(T) = \varepsilon_r(1 - i \tan\delta(T))$ with optional TLS saturation (`tls_tan_delta_0`, `tls_E_c`, `tls_T_0`, `tls_beta`), oxide interface participation thickness and loss tangent, and Mattis–Bardeen superconducting thin-film parameters (`tc`, `gap_0`, `lambda_L`, `xi_0`, `sigma_n`, `dynes_gamma`, `r_res`, `nominal_thickness`).
- **Coupled Tensors:** Full crystal-frame `[piezoelectric]` (`e` in `C/m^2` or `d` in `pC/N` with thermodynamic check $C^D - C^E = e^T (\varepsilon^S)^{-1} e$), `[photoelastic]` ($p_{IJ}$), `[electro_optic]` ($r_{IJ}$), and positive semi-definite `[viscoelastic]` dynamic viscosity $\eta_{IJ}$ for Akhiezer/Kelvin–Voigt $Q_m$ damping.

## Inheritance (`extends`) & Provenance

Every table carries a `source` citation and optional `verified` date; unverified estimates must carry `"APPROXIMATE"` or `"UNVERIFIED"`. Custom or measured materials inherit cleanly via `extends = "silicon"` without duplicating tensors: each overridden table replaces its parent table atomically so units and parameters never mix across sources.

## Bundled Domains & CLI Evaluator

The shipped library spans five domains — **base optical/RF/thermal** (`air`, `aluminium`, `copper`, `gold`, `lithium_niobate`, `silica`, `silicon`, `silicon_nitride`, `silver`), **photonics author/process variants** (`silicon_salzberg`, `silicon_palik`, `silicon_nitride_philipp`, `silicon_nitride_pecvd`, `lithium_niobate_mgo`, `gallium_arsenide`, `gallium_arsenide_papatryfonos`, `gold_johnson`, `silver_johnson`), **phononics, crystals & MEMS** (`silicon_mcskimin`, `silica_thermal_box`, `quartz_alpha`, `polysilicon_polymumps`, `aluminum_nitride`, `aluminum_nitride_piezomumps`, `scaln_30`, `diamond`), **cryogenic 4 K baselines** (`silicon_4k`, `silica_4k`, `lithium_niobate_4k`, `copper_4k`), and **superconducting films** (`aluminium_sc_emft`, `niobium_sc`, `titanium_nitride_sc`, `niobium_titanium_nitride_sc`, `granular_aluminium_sc`, `nbn_sc`, `beta_tin_sc`).

```sh
nyx material list
nyx material show niobium_sc
nyx material eval niobium_sc --temp 0.02 --freq 5e9 --thickness 100e-9
nyx material bcs niobium_sc --freq 5e9 --temp 0.02 --thickness 100e-9
nyx material fetch rii main/Si/nk/Green-2008 --extends silicon --name silicon_green --save
```
