# Writing a material

A material is one small text file in the TOML format, for example `materials/silicon.toml`. Nothing is assumed: write the unit of every number. If something is wrong, nyx says which material, which table and which key, and what to do. The format is version 2 (`format = 2` on the first line; legacy `format = 1` files continue to load through a read-only compatibility shim). The formal description is the JSON Schema `docs/architecture/contracts/material_schema.v2.json`; `docs/architecture/contracts/material_kitchen_sink.v2.toml` is one file that uses every section. Decisions: ADR 29, 30, 101 and 102.

## The shortest useful file

```toml
format = 2
name = "glass"
kind = "dielectric"
density = "2500 kg/m^3"

[elastic]
unit = "GPa"
E = 70
nu = 0.22

[optical]
model = "constant"
n = 1.5
```

## Quantities and units

Write `"496.8 GPa"`, or `{ value = 496.8, unit = "GPa" }`, or a bare number when the table declares the unit (`unit = "GPa"` in `[elastic]`, `conductivity_unit = "W/(m*K)"` next to a thermal value). Dimensionless numbers (Poisson's ratio, relative permittivity, loss tangent, photoelastic $p_{IJ}$) are bare.

| Kind | Examples |
|---|---|
| stiffness / stress | `Pa kPa MPa GPa TPa` (residual stress may be signed) |
| density | `kg/m^3`, `g/cm^3` |
| length (wavelength, thickness) | `m mm um nm pm angstrom` (`µm` also) |
| thermal conductivity | `W/(m*K)` |
| specific heat | `J/(kg*K)` |
| per kelvin (expansion, dn/dT) | `1/K`, `ppm/K`, `1e-6/K` |
| compliance | `1/GPa`, `1/TPa`, `1/Pa`, `m^2/N` |
| electrical conductivity / resistivity | `S/m`, `ohm*m`, `uohm*cm` |
| temperature | `K`, `mK`, `uK`, `degC` |
| frequency | `Hz`, `kHz`, `MHz`, `GHz`, `THz` |
| energy (BCS gap) | `eV`, `meV`, `ueV`, `J` |
| surface resistance / sheet inductance | `ohm`, `mohm`, `uohm`, `nohm` / `H`, `nH`, `pH` |
| piezoelectric / electro-optic | `C/m^2` ($e_{iJ}$), `pC/N`, `pm/V` ($d_{iJ}$, $r_{IJ}$) |
| dynamic viscosity | `Pa*s`, `mPa*s` |

## Top level

`format` (required, `2` or legacy `1`), `name`, `description`, `kind` (`dielectric`, `metal`, `superconductor`, `semiconductor`, `piezoelectric`, `fluid`), `point_group` (one of the 32 crystallographic Hermann–Mauguin classes, e.g. `"m-3m"`, `"3m"`, `"6mm"`, `"4/mmm"`), `crystal_setting` (IEEE Std 176-1987 axis convention string), `reference_temperature` (default `"293.15 K"`), `residual_stress` (signed thin-film stress), `source` (where the numbers come from, see "Provenance"), `notes`, `density`, `extends` (see "Starting from another material"). The orientation of a crystal is not in the material: the project places it with Bunge Euler angles.

## `[elastic]`: the stiffness

The stiffness is the 6x6 matrix in the textbook Voigt order **1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy** with engineering shear strain, in the crystal frame. Give it in exactly one of these ways:

1. **Class and constants.** `class` is `isotropic`, `cubic`, `hexagonal`, `trigonal`, `tetragonal`, `orthotropic`, `monoclinic` (with `unique_axis = "y"` (default) or `"z"`) or `triclinic`. The constants are `cIJ` with I at most J; the class says which are independent:

   | Class | Independent constants |
   |---|---|
   | isotropic | c11 c12 |
   | cubic | c11 c12 c44 |
   | hexagonal | c11 c12 c13 c33 c44 |
   | trigonal (32, 3m, -3m) | c11 c12 c13 c14 c33 c44 |
   | tetragonal (4mm, 422, -42m, 4/mmm) | c11 c12 c13 c33 c44 c66 |
   | orthotropic | c11 c12 c13 c22 c23 c33 c44 c55 c66 |
   | monoclinic, y | the nine above and c15 c25 c35 c46 |
   | monoclinic, z | the nine above and c16 c26 c36 c45 |
   | triclinic | all 21 |

   You may also write a dependent entry (for example `c66` of a hexagonal crystal); it must agree with the others.
2. **The matrix**: `matrix = [...]`, the full 6x6 or the upper triangle as rows of 6, 5, 4, 3, 2, 1 numbers. A full matrix must be symmetric. Add `class = "cubic"` to have nyx check the matrix against that class.
3. **The compliance**: `compliance = [...]` (same shapes; it is the inverse matrix of the stiffness in Voigt form, `compliance_unit = "1/TPa"`), or `class` with `s11 s12 ...` and `compliance_unit`.
4. **Isotropic pairs**: `E` and `nu`, or `lambda` and `mu`, or `K` and `G`.

```toml
[elastic]                  # sapphire, by class
class = "trigonal"
unit = "GPa"
c11 = 496.8
c12 = 163.6
c13 = 110.9
c14 = -23.5
c33 = 498.1
c44 = 147.4
```

```toml
[elastic]                  # the same silicon as a matrix (upper triangle)
unit = "GPa"
matrix = [
  [165.7, 63.9, 63.9,  0.0,  0.0,  0.0],
  [       165.7, 63.9, 0.0,  0.0,  0.0],
  [              165.7, 0.0,  0.0,  0.0],
  [                     79.6, 0.0,  0.0],
  [                           79.6, 0.0],
  [                                 79.6],
]
```

What is checked: symmetry (it names the two entries that differ), stability (positive definite; for a cubic crystal the message gives c11 > |c12|, c11 + 2 c12 > 0, c44 > 0), consistency with the declared class (it names the entry and the relation), the range of Poisson's ratio, unknown keys.

## `[thermal]` and `[electrical]`: tensors and temperature curves

`thermal`: `conductivity`, `specific_heat`, `expansion`, plus optional continuous temperature curves `conductivity_T`, `specific_heat_T`, `expansion_T` (each using `model = "tabulated_T"` with `data = [[T, val], ...]`, `model = "log10_poly"` with NIST logarithmic polynomial `coefficients` and `range = [T_min, T_max]`, or `model = "einstein_shift"`). `[elastic]` also accepts `c_T = { c11 = ..., c12 = ..., c44 = ... }`. `electrical`: `permittivity` (relative, low frequency), `permeability` (relative), `conductivity`. A tensor is given as one number (isotropic), three (x y z), six (xx yy zz yz xz xy), a 3x3 array, or `{ perp = .., par = .. }` for a crystal with its unique axis along crystal z.

```toml
[thermal]
conductivity = { perp = 23.0, par = 25.2 }
conductivity_unit = "W/(m*K)"
```

Permittivity and thermal conductivity must be positive in every direction, electrical conductivity not negative; expansion may be negative.

## `[optical]`: dispersion

Optical properties are models evaluated at the study wavelength. Write the wavelength unit once, `wavelength_unit = "um"`; the coefficients, the `range` and a table are in that unit. The complex index is `n - j k` (loss is a positive `k`).

| `model` | Keys | Form |
|---|---|---|
| `constant` | `n` or `eps`, optional `k` | no dispersion |
| `sellmeier` | `B`, and `C` (squared resonance wavelengths) or `lambda0` (resonance wavelengths), optional `A` (default 1) | `n^2 = A + sum B x^2 / (x^2 - C)` |
| `cauchy` | `coefficients = [A, B, C, ...]` | `n = A + B/x^2 + C/x^4 + ...` |
| `tabulated` | `data = [[wavelength, n], ...]` (or 3-column `[wavelength, n, k]`), optional separate `k_data = [[wavelength, k], ...]`, optional `interpolation = "linear"` (default monotone cubic `"pchip"`) | table |
| `schott` | `coefficients = [a0, a1, a2, a3, a4, a5]` | `n^2 = a0 + a1 x^2 + a2 x^-2 + a3 x^-4 + a4 x^-6 + a5 x^-8` (Schott glass catalogue) |
| `herzberger` | `coefficients = [A, B, C, D, E, F]` | `n = A + B L + C L^2 + D x^2 + E x^4 + F x^6`, `L = 1/(x^2 - 0.028 um^2)` (the 0.028 is converted to your wavelength unit) |
| `rii_formula_4` | `coefficients = [c0..c8, ...]` | RefractiveIndex.INFO Formula 4 (exact power-law poles and polynomial terms) |
| `rii_formula_5` | `coefficients = [c0, c1, c2, ...]` | RefractiveIndex.INFO Formula 5 (`n = c0 + sum c_i x^c_{i+1}`) |
| `rii_formula_6` | `coefficients = [c0, c1, c2, ...]` | RefractiveIndex.INFO Formula 6 (`n - 1 = c0 + sum B_i / (C_i - x^-2)`) |
| `drude` | `frequency_unit`, `plasma_frequency`, `damping`, optional `eps_inf` (default 1), `drude_weight` (default 1) | `eps = eps_inf - f0 Wp^2 / (w^2 - j g w)` (free carriers, metals) |
| `lorentz` | `frequency_unit`, `plasma_frequency`, `oscillators = [[f, resonance, damping], ...]`, optional `eps_inf` | `eps = eps_inf + sum f Wp^2 / (w0^2 - w^2 + j g w)` |
| `drude_lorentz` | both of the above | the sum of both (Rakic et al. 1998 write metals this way) |

Rules: `range = [min, max]` is required for every formula model, i.e. everything but `constant` and `tabulated` (published coefficients hold only there); a table's range is its extent. **A wavelength outside the range is an error**, never an extrapolation; `extrapolate = "linear"` in a table allows it on purpose. `dn_dT = "1.87e-4 1/K"` adds a linear thermo-optic shift `dn_dT (T - reference_temperature)` (default reference 293.15 K) when a temperature is given, or `dn_dT_model = { model = "einstein", ... }` evaluates a cryogenic Bose–Einstein index shift down to 4 K. A uniaxial crystal uses the tables `[optical.ordinary]` and `[optical.extraordinary]` (unique axis along crystal z), a biaxial one `[optical.x]`, `[optical.y]`, `[optical.z]`; each holds its own model.

```toml
[optical]
wavelength_unit = "um"
range = [0.21, 6.7]
model = "sellmeier"          # fused silica, Malitson 1965
B = [0.6961663, 0.4079426, 0.8974794]
lambda0 = [0.0684043, 0.1162414, 9.896161]
```

## `[dielectric]`, `[interfaces]` and `[superconducting]`: microwave & cQAD models

`[dielectric]` separates microwave/RF permittivity and dielectric loss tangent from `[electrical]` DC conductivity: `permittivity` (scalar or tensor), `loss_tangent` (scalar or tensor $\tan\delta \ge 0$), optional temperature curve `loss_tangent_T`, and optional two-level-system (TLS) saturation parameters (`tls_tan_delta_0`, `tls_E_c`, `tls_T_0`, `tls_beta`), yielding $\tilde{\varepsilon}_r(T) = \varepsilon_r (1 - i\tan\delta(T))$ in `ems` microwave studies. `[interfaces.{ma,ms,sa}]` records thin-interface thickness `t` and `loss_tangent` for metal–air, metal–substrate, and substrate–air participation analyses.

`[superconducting]` parameterizes the Mattis–Bardeen complex conductivity $\sigma(\omega, T) = \sigma_1(\omega, T) - i\sigma_2(\omega, T)$ and surface impedance $Z_s(\omega, T, d) = R_s + i\omega L_{s,\square}$ (used by the `ems` `surface_impedance` boundary condition and `nyx material bcs`):

```toml
format = 2
name = "nb_film_example"
kind = "superconductor"
point_group = "m-3m"

[superconducting]
tc = "9.2 K"
gap_0 = "1395 ueV"
lambda_L = "33.3 nm"
xi_0 = "39.0 nm"
sigma_n = "5.5e7 S/m"
dynes_gamma = "0.5 ueV"
r_res = "8.0 nOhm/sq"
x_qp_floor = 1.0e-7
nominal_thickness = "100 nm"
film_model = "kerr"

[interfaces.ma]
oxide = "Nb2O5"
thickness = "4.5 nm"
eps_r = 33.0
tan_delta_tls_0 = 8.0e-4
e_crit = "12.0 V/m"
```

## Coupled-field & damping tensors: `[piezoelectric]`, `[photoelastic]`, `[electro_optic]`, `[viscoelastic]`

`format = 2` validates and stores four coupled-field and acoustic-damping tensors in the crystal frame (rotated by the region's Bunge Euler angles):

* `[piezoelectric]`: `form = "e"` ($e_{iJ}$ in `C/m^2`) or `form = "d"` ($d_{iJ}$ in `pC/N` or `pm/V`), given by crystallographic symmetry keys (`e15`, `e22`, `e31`, `e33`, `d11`, `d14`, ...) or `matrix = [[...], [...], [...]]` ($3\times 6$), plus `clamped_permittivity` ($\varepsilon^S / \varepsilon_0$) or `free_permittivity` ($\varepsilon^T / \varepsilon_0$). `MaterialSpec` converts between $e$ and $d$ with the Voigt shear weight and enforces the thermodynamic link $C^D - C^E = e^T (\varepsilon^S)^{-1} e$ and $\varepsilon^T - \varepsilon^S = e (C^E)^{-1} e^T$.
* `[photoelastic]`: elasto-optic tensor $p_{IJ}$ ($6\times 6$, dimensionless), given by `class` (`isotropic`, `cubic`, `trigonal`, `hexagonal`, `orthotropic`, `triclinic`) with keys `p11`, `p12`, `p13`, `p14`, `p31`, `p33`, `p41`, `p44`, or `matrix`.
* `[electro_optic]`: Pockels tensor $r_{IJ}$ ($6\times 3$, default `unit = "pm/V"`), with `condition = "clamped" | "free"` and keys `r13`, `r22`, `r33`, `r41`, `r42`, `r51`, or `matrix`.
* `[viscoelastic]`: Kelvin–Voigt / Akhiezer dynamic viscosity tensor $\eta_{IJ}$ ($6\times 6$, default `unit = "mPa*s"`, checked positive semi-definite), assembled into the viscous damping matrix $\mathbf{C}_{\text{visc}}$ by `sms` eigenfrequency solves to return mechanical quality factors $Q_m$.

## The library

`materials/*.toml` is loaded at run time with `Library::load_dir`; a material is found by its file name, ignoring case. Shipped: 40 `format = 2` materials across five domains (`materials/examples/silicon_measured.toml` remains at `format = 1` as the legacy compatibility fixture). Every table of every file names its source; values that could not be checked against a primary publication carry `"APPROXIMATE"` or `"UNVERIFIED"` in that source (the library test enforces this for the shipped files).

## Provenance: where a number comes from

The top of a file and each table (`[elastic]`, `[thermal]`, `[electrical]`, `[dielectric]`, `[superconducting]`, `[interfaces.{ma,ms,sa}]`, `[piezoelectric]`, `[photoelastic]`, `[electro_optic]`, `[viscoelastic]`, `[optical]`, and each axis of `[optical]` if it differs) may carry `source`: a citation text, a table `{ cite = "...", url = "https://...", doi = "10.xxxx/...", verified = "2026-11-04" }`, or a list of those. `verified` is the date on which someone compared the values with the publication; an entry without it must say "approximate" or "unverified" in its text (the library test enforces this for the shipped files).

## Starting from another material

To change a few values of a library material (your own measured index, a different density) write a short file that starts from it:

```toml extends
format = 2
extends = "silicon"
source = "Own ellipsometry, 2026-11 (example)"

[optical]
wavelength_unit = "um"
range = [1.5, 1.6]
model = "tabulated"
data = [[1.50, 3.4799, 0.0], [1.55, 3.4761, 0.0], [1.60, 3.4722, 0.0]]
```

Rules: the base is looked up by name in the same library (a missing base or a loop is reported with the chain); every top-level key you write replaces the base's; each table (`[elastic]`, `[thermal]`, `[electrical]`, `[dielectric]`, `[superconducting]`, `[interfaces]`, `[piezoelectric]`, `[photoelastic]`, `[electro_optic]`, `[viscoelastic]`, `[optical]`) is replaced as a whole, never mixed key by key, so a half-replaced stiffness cannot occur and a unit never follows a value it was not written for; the `name` is not inherited; all checks run on the combined result. Read such a file with `library.read(text, "name")`; `Library::load_dir` resolves `extends` between the files of the folder. `materials/examples/silicon_measured.toml` is an example.

## From Rust

```rust
let lib = femcore::material::Library::load_dir(std::path::Path::new("materials"))?;
let sapphire = lib.require("sapphire")?;
let aniso = sapphire.anisotropy([0.0, 90.0, 0.0]);                 // for the existing kernels
let eps = sapphire.optical_permittivity(1.55e-6, Some(300.0), [0.0, 90.0, 0.0])?;   // 3x3 complex, wafer frame
```

## Bundled multi-source library (40 materials)

Because optical, mechanical, thermal, coupled-field and superconducting thin-film parameters depend strongly on author dataset, deposition process and temperature, conflicting sources are shipped as separate named TOML files rather than overwritten silently. Every shipped material carries a one-line `description` stating what distinguishes it from sibling variants.

| Domain | Material names | Key differences |
|---|---|---|
| Base optical / RF / thermal (10) | `air`, `aluminium`, `copper`, `gold`, `lithium_niobate`, `sapphire`, `silica`, `silicon`, `silicon_nitride`, `silver` | Room-temperature `format = 2` baselines with `[dielectric]`, `[photoelastic]`, `[viscoelastic]`, `[piezoelectric]` and `[electro_optic]` where applicable (`silicon`: Li 1980 + Biegelsen $p_{IJ}$ + Helme-King $\eta_{IJ}$; `lithium_niobate`: Zelmon 1997 + Weis-Gaylord $e_{iJ}, p_{IJ}, r_{IJ}$; `sapphire`: Malitson-Dodge + Gladden + Krupka microwave loss) |
| Photonics author / process variants (9) | `silicon_salzberg`, `silicon_palik`, `silicon_nitride_philipp`, `silicon_nitride_pecvd`, `lithium_niobate_mgo`, `gallium_arsenide`, `gallium_arsenide_papatryfonos`, `gold_johnson`, `silver_johnson` | `silicon_salzberg` (Salzberg-Villa 1957 IR) & `silicon_palik` (Palik tabular 0.25-20 um) extend `silicon`; `silicon_nitride_philipp` & `silicon_nitride_pecvd`; `lithium_niobate_mgo` (`rii_formula_4`, 5 mol% MgO Gayer 2008); `gallium_arsenide` (`-43m` with $e_{14}, r_{41}, p_{IJ}$) vs. `gallium_arsenide_papatryfonos`; `gold_johnson` / `silver_johnson` |
| Phononics, crystals & MEMS (8) | `silicon_mcskimin`, `silica_thermal_box`, `quartz_alpha`, `polysilicon_polymumps`, `aluminum_nitride`, `aluminum_nitride_piezomumps`, `scaln_30`, `diamond` | `silicon_mcskimin` (McSkimin 1953 ultrasonics 298 K); `silica_thermal_box`; `quartz_alpha` (point group `32`, Bechmann $d_{11}, d_{14}$, Narasimhamurty $p_{IJ}$, Lamb $\eta_{IJ}$); `polysilicon_polymumps`; `aluminum_nitride` (`6mm`, Tsubouchi $e_{iJ}$, Davydov $p_{IJ}$) vs. `aluminum_nitride_piezomumps` ($d_{33}, d_{31}$) vs. `scaln_30` (30% ScAlN, ~2x piezoelectric); `diamond` |
| Cryogenic 4 K (5) | `silicon_4k`, `silica_4k`, `sapphire_4k`, `lithium_niobate_4k`, `copper_4k` | Cryogenic 4 K baselines with continuous `conductivity_T` (`tabulated_T` / NIST `log10_poly`), cryogenic `[dielectric]` $\tan\delta$, and `[interfaces.sa]` |
| Superconducting films & metals (8) | `aluminium_sc_emft`, `tantalum_alpha_sc`, `niobium_sc`, `titanium_nitride_sc`, `niobium_titanium_nitride_sc`, `granular_aluminium_sc`, `nbn_sc`, `beta_tin_sc` | Full `[superconducting]` Mattis–Bardeen / kinetic-inductance tables ($T_c$, $\Delta_0$, $\lambda_{L,0}$, $\xi_0$, $\rho_n$, $\ell$, $d$, `regime`, `subgap_loss_floor`) plus `[interfaces.{ma,ms}]` (`nbn_sc`: 16 K NbN film; `beta_tin_sc`: 3.722 K tetragonal $\beta$-Sn) |

## Evaluating, searching and downloading materials

`nyx` evaluates temperature-, frequency- and thickness-dependent properties directly from the CLI, engine scripts, MCP and the GUI (**Materials** tab), and converts remote records from `refractiveindex.info`, `AFLOW`, `Materials Project` and `NumBAT` on demand into validated `format = 2` TOML files marked `"UNVERIFIED"` in `[source].cite`:

```sh
# List bundled materials (with format=2 capability tags) or curated remote presets
nyx material list
nyx material list --remote --domain superconducting

# Inspect any bundled material or evaluate its constitutive state at (T, f, d)
nyx material show niobium_sc
nyx material eval niobium_sc --temp 0.02 --freq 5e9 --thickness 100e-9
nyx material bcs niobium_sc --freq 5e9 --temp 0.02 --thickness 100e-9

# Search curated presets and live AFLOW / Materials Project endpoints
nyx material search silicon
nyx material search Si --live --source aflow

# Download a preset, a refractiveindex.info / AFLOW / NumBAT entry, or a Materials Project ID (uses --api-key or MP_API_KEY via curl stdin)
nyx material fetch si_green_2008 --out si_green.toml
nyx material fetch rii main/Si/nk/Green-2008 --extends silicon --name silicon_green --save
nyx material fetch mp mp-149 --api-key "$MP_API_KEY" --name silicon_mp --out silicon_mp.toml

# Convert a local download offline (supports rii YAML, aflow JSON, mp JSON, numbat JSONC)
nyx material import mp mp_149_si.json --name silicon_mp --extends silicon --out silicon_mp.toml
```

Within engine command scripts, MCP sessions and the GUI (**Materials > Download or import material** and **Live (T, f, d) evaluator**), the equivalent engine commands are `evaluate_material`, `bcs_surface_impedance`, `list_remote_materials`, `search_remote_materials`, `import_material` and `fetch_material`.
