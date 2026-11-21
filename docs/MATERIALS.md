# Writing a material

A material is one small text file in the TOML format, for example `materials/silicon.toml`. Nothing is assumed: write the unit of every number. If something is wrong, nyx says which material, which table and which key, and what to do. The format is version 1 (`format = 1` on the first line). The formal description is the JSON Schema `docs/architecture/contracts/material_schema.v1.json`; `docs/architecture/contracts/material_kitchen_sink.v1.toml` is one file that uses every section. Decisions: ADR 29 and 30.

## The shortest useful file

```toml
format = 1
name = "glass"
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

Write `"496.8 GPa"`, or `{ value = 496.8, unit = "GPa" }`, or a bare number when the table declares the unit (`unit = "GPa"` in `[elastic]`, `conductivity_unit = "W/(m*K)"` next to a thermal value). Dimensionless numbers (Poisson's ratio, relative permittivity) are bare.

| Kind | Examples |
|---|---|
| stiffness | `Pa kPa MPa GPa TPa` |
| density | `kg/m^3`, `g/cm^3` |
| length (wavelength) | `m mm um nm pm angstrom` (`µm` also) |
| thermal conductivity | `W/(m*K)` |
| specific heat | `J/(kg*K)` |
| per kelvin (expansion, dn/dT) | `1/K`, `ppm/K`, `1e-6/K` |
| compliance | `1/GPa`, `1/TPa`, `1/Pa`, `m^2/N` |
| electrical conductivity | `S/m` |
| temperature | `K`, `degC` |

## Top level

`format` (required), `name`, `description`, `source` (where the numbers come from, see "Provenance"), `notes`, `density`, `extends` (see "Starting from another material"). The orientation of a crystal is not in the material: the project places it with Bunge Euler angles.

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

## `[thermal]` and `[electrical]`: tensors

`thermal`: `conductivity`, `specific_heat`, `expansion`. `electrical`: `permittivity` (relative, low frequency), `permeability` (relative), `conductivity`. A tensor is given as one number (isotropic), three (x y z), six (xx yy zz yz xz xy), a 3x3 array, or `{ perp = .., par = .. }` for a crystal with its unique axis along crystal z.

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
| `tabulated` | `data = [[wavelength, n], ...]` or with a third column `k`; optional `interpolation = "linear"` (default monotone cubic `"pchip"`) | table |
| `schott` | `coefficients = [a0, a1, a2, a3, a4, a5]` | `n^2 = a0 + a1 x^2 + a2 x^-2 + a3 x^-4 + a4 x^-6 + a5 x^-8` (Schott glass catalogue) |
| `herzberger` | `coefficients = [A, B, C, D, E, F]` | `n = A + B L + C L^2 + D x^2 + E x^4 + F x^6`, `L = 1/(x^2 - 0.028 um^2)` (the 0.028 is converted to your wavelength unit) |
| `drude` | `frequency_unit`, `plasma_frequency`, `damping`, optional `eps_inf` (default 1), `drude_weight` (default 1) | `eps = eps_inf - f0 Wp^2 / (w^2 - j g w)` (free carriers, metals) |
| `lorentz` | `frequency_unit`, `plasma_frequency`, `oscillators = [[f, resonance, damping], ...]`, optional `eps_inf` | `eps = eps_inf + sum f Wp^2 / (w0^2 - w^2 + j g w)` |
| `drude_lorentz` | both of the above | the sum of both (Rakic et al. 1998 write metals this way) |

Rules: `range = [min, max]` is required for every formula model, i.e. everything but `constant` and `tabulated` (published coefficients hold only there); a table's range is its extent. **A wavelength outside the range is an error**, never an extrapolation; `extrapolate = "linear"` in a table allows it on purpose. `dn_dT = "1.87e-4 1/K"` adds a thermo-optic shift `dn_dT (T - reference_temperature)` (default reference 293.15 K) when a temperature is given. A uniaxial crystal uses the tables `[optical.ordinary]` and `[optical.extraordinary]` (unique axis along crystal z), a biaxial one `[optical.x]`, `[optical.y]`, `[optical.z]`; each holds its own model.

```toml
[optical]
wavelength_unit = "um"
range = [0.21, 6.7]
model = "sellmeier"          # fused silica, Malitson 1965
B = [0.6961663, 0.4079426, 0.8974794]
lambda0 = [0.0684043, 0.1162414, 9.896161]
```

## The library

`materials/*.toml` is loaded at run time with `Library::load_dir`; a material is found by its file name, ignoring case. Shipped: `silicon`, `silica`, `sapphire`, `air` (mechanical, thermal, electrical and optical data), `silicon_nitride` (LPCVD film: density, isotropic stiffness and Luke 2015 optical dispersion), `lithium_niobate` (Weis and Gaylord 1985 density and trigonal constant-E stiffness, Zelmon 1997 ordinary and extraordinary optical dispersion) and, with optical data only, `gold`, `silver`, `copper`, `aluminium` (Lorentz-Drude fits of Rakic et al. 1998, 0.25 to 6 um for gold, to 12 um for silver and copper, 0.06 to 248 um for aluminium). `materials/examples/` holds files that are not part of the library, for example `silicon_measured.toml`. Every table of every file names its source; values that could not be checked against a publication say "approximate" or "unverified" in that source (silicon thermal conductivity and permittivity, silica permittivity, the sapphire thermal, expansion and permittivity values, silicon nitride thin-film stiffness, air).

## Provenance: where a number comes from

The top of a file and each table (`[elastic]`, `[thermal]`, `[electrical]`, `[optical]`, and each axis of `[optical]` if it differs) may carry `source`: a citation text, a table `{ cite = "...", url = "https://...", doi = "10.xxxx/...", verified = "2026-11-04" }`, or a list of those. `verified` is the date on which someone compared the values with the publication; an entry without it must say "approximate" or "unverified" in its text (the library test enforces this for the shipped files). The loader keeps the sources (`MaterialSpec.source`, `Elastic.source`, `Thermal.source`, `Electrical.source`, `Optical.source`, `AxisModel.source`) and writes them back.

## Starting from another material

To change a few values of a library material (your own measured index, a different density) write a short file that starts from it:

```toml extends
format = 1
extends = "silicon"
source = "Own ellipsometry, 2026-11 (example)"

[optical]
wavelength_unit = "um"
range = [1.5, 1.6]
model = "tabulated"
data = [[1.50, 3.4799, 0.0], [1.55, 3.4761, 0.0], [1.60, 3.4722, 0.0]]
```

Rules: the base is looked up by name in the same library (a missing base or a loop is reported with the chain); every top-level key you write replaces the base's; a table (`[elastic]`, `[thermal]`, `[electrical]`, `[optical]`) is replaced as a whole, never mixed key by key, so a half-replaced stiffness cannot occur and a unit never follows a value it was not written for; the `name` is not inherited; all checks run on the combined result. Read such a file with `library.read(text, "name")`; `Library::load_dir` resolves `extends` between the files of the folder. `materials/examples/silicon_measured.toml` is an example.

## From Rust

```rust
let lib = femcore::material::Library::load_dir(std::path::Path::new("materials"))?;
let sapphire = lib.require("sapphire")?;
let aniso = sapphire.anisotropy([0.0, 90.0, 0.0]);                 // for the existing kernels
let eps = sapphire.optical_permittivity(1.55e-6, Some(300.0), [0.0, 90.0, 0.0])?;   // 3x3 complex, wafer frame
```
