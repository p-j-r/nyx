# Material specification: code

C4 level 4 (code) of the component **Material specification** (`femcore.material`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> TOML material files with explicit units, tensors, dispersion models evaluated at the wavelength, Zimmermann BCS superconducting conductivity, provenance and the library loader.

**Decisions:** [ADR 29](../../../adr/0029-material-files-are-toml-with-explicit-units.md), [ADR 30](../../../adr/0030-optical-dispersion-models-evaluated-at-the-study-wavelength.md), [ADR 34](../../../adr/0034-a-material-may-extend-another-and-replaces-whole-tables.md), [ADR 35](../../../adr/0035-material-provenance-oscillator-models-and-the-checked-library.md), [ADR 101](../../../adr/0101-sms-waveguide-modes-and-floquet-pml.md), [ADR 104](../../../adr/0104-bcs-surface-impedance-and-lumped-josephson-ports.md)  
**Size:** 10 files, 4406 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/material/bcs.rs` | 721 |  |
| `crates/femcore/src/material/elastic.rs` | 784 |  |
| `crates/femcore/src/material/mod.rs` | 402 |  |
| `crates/femcore/src/material/optics.rs` | 916 |  |
| `crates/femcore/src/material/reader.rs` | 269 |  |
| `crates/femcore/src/material/schema.rs` | 192 |  |
| `crates/femcore/src/material/source.rs` | 136 |  |
| `crates/femcore/src/material/spec.rs` | 594 |  |
| `crates/femcore/src/material/tensor.rs` | 171 |  |
| `crates/femcore/src/material/units.rs` | 221 |  |

## What it uses

Nothing outside the component (the standard library and external crates are not listed).

## Used by

| Component | Modules | Items |
|---|---|---|
| [dfs: Physics adapter](../dfs/adapter.md) | material::spec | Library |
| [dfs: Shared helpers](../dfs/common.md) | material | M3 |
| [dfs: Electrostatics kernel](../dfs/electro.md) | material | M3 |
| [dfs: Heat kernel](../dfs/heat.md) | material | M3 |
| [dfs: Typed problem](../dfs/problem.md) | material, material::bcs | FilmModel, bcs_surface_impedance |
| [ems: Physics adapter](../ems/adapter.md) | material, material::spec | Library, tensor_global |
| [ems: Thin-film IBC and Josephson ROM](../ems/ibc.md) | material | tensor_global |
| [ems: Typed problem](../ems/problem.md) | material, material::bcs | FilmModel, bcs_surface_impedance |
| [engine: Commands](../engine/commands.md) | material, material::spec | Anisotropy, Library |
| [engine: Job runner and queue](../engine/jobs.md) | material::spec | Library |
| [engine: Materials](../engine/materials.md) | material, material::bcs, material::spec | EV_TO_J, FilmModel, Library, MaterialSpec, bcs_surface_impedance, surface_impedance_from_params |
| [engine: Results files](../engine/resultsIO.md) | material::spec | Library |
| [Finite-element basics](./fem.md) | material | Anisotropy, M3, M6 |
| [Meshing facade](./meshing.md) | material | M6 |
| [Frame and sections](./section.md) | material | M3, M6 |
| [physics: Problem envelope](../physics/envelope.md) | material, material::elastic, material::spec | Class, Library, Tensor |
| [physics: Physics trait and registry](../physics/trait.md) | material::spec | Library |
| [sms: Physics adapter](../sms/adapter.md) | material, material::spec | Library, isotropic_stiffness, stiffness_global |
| [sms: Axisymmetric kernel](../sms/axisym.md) | material | M6 |
| [sms: Cross-section study kernel](../sms/crosssection.md) | material | M6 |
| [sms: Low-reflecting boundary](../sms/lrb.md) | material | M6 |
| [sms: Solid kernel](../sms/solid.md) | material | M3, M6 |
| [sms: Waveguide kernel](../sms/waveguide.md) | material | M6 |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace material_bcs_rs {
    class FilmModel {
        <<enumeration>>
        TwoSided
        Kerr
        ThinLimit
        +parse(s) Option~FilmModel~
        +as_str() str
    }
    class BcsEval {
        <<struct>>
        +sigma Complex64
        +sigma_n f64
        +delta_j f64
        +lambda_eff f64
        +z_s Complex64
        +r_s f64
        +l_k_sq f64
        +l_eff f64
        +more_1_fields
    }
    class fns_material_bcs_rs {
        <<module material/bcs.rs>>
        +bcs_gap(tc, gap_0_j, temperature) Result
        +surface_impedance_from_params(tc, gap_0_j, lambda_l, xi_0, mean_free_path, sigma_n, freq_hz, temperature, thickness, film_model) Result
        +bcs_surface_impedance(sc, freq_hz, temperature, thickness, film_model) Result
        +K_B const
        +HBAR const
        +EV_TO_J const
        +MU0 const
        +more_1_items
    }
  }
  namespace material_elastic_rs {
    class Unique {
        <<enumeration>>
        Y
        Z
    }
    class Class {
        <<enumeration>>
        Isotropic
        Cubic
        Hexagonal
        Trigonal
        Tetragonal
        Orthotropic
        Monoclinic
        Triclinic
        +name() str
        +parse(s, unique) Option~Class~
        +constants() Vec
    }
    class Elastic {
        <<struct>>
        +c M6
        +class Option~Class~
        +exact_constants bool
        +source Vec~Source~
        +from_matrix(c) Result
        +from_constants(class, constants) Result
        +isotropic(e, nu) Result
        +compliance() M6
        +detect_class() Class
        +global(euler_deg) M6
        +eigenvalue_range() tuple
        +read(s) Result
        +more_1_methods()
    }
    class fns_material_elastic_rs {
        <<module material/elastic.rs>>
        ~sym_eigenvalues(a) slice
        +invert6(m) Option~M6~
        ~NAMES const
    }
  }
  namespace material_mod_rs {
    class Stiffness {
        <<enumeration>>
        Cubic
        Hexagonal
        Trigonal
        Orthotropic
        General
        +name() str
        +constants_mut() Vec
    }
    class Tensor {
        <<enumeration>>
        Uniaxial
        Diagonal
        Symmetric
        +name() str
        +matrix() M3
        +constants_mut() Vec
    }
    class Anisotropy {
        <<struct>>
        +stiffness Option~Stiffness~
        +conductivity Option~Tensor~
        +permittivity Option~Tensor~
        +expansion Option~Tensor~
        +euler_deg slice
    }
    class fns_material_mod_rs {
        <<module material/mod.rs>>
        +rotation(euler_deg) M3
        +euler_from_rotation(m) slice
        +rotate_tensor(t, r) M3
        +stiffness_global(s, euler_deg) M6
        +isotropic_stiffness(e, nu) M6
        +tensor_global(t, euler_deg) M3
        +iso_tensor(v) M3
        +is_positive_definite(c) bool
        +sapphire() tuple
        +silicon() Stiffness
        +M3 type
        +M6 type
    }
  }
  namespace material_optics_rs {
    class Interp {
        <<enumeration>>
        Pchip
        Linear
    }
    class Extrapolate {
        <<enumeration>>
        Error
        Linear
    }
    class Dispersion {
        <<enumeration>>
        Constant
        Sellmeier
        Cauchy
        Tabulated
        Schott
        Herzberger
        Oscillators
        +name() str
    }
    class AxisModel {
        <<struct>>
        +dispersion Dispersion
        +source Vec~Source~
        +dn_dt f64
    }
    class Axes {
        <<enumeration>>
        Isotropic
        Uniaxial
        Biaxial
        +model(axis) AxisModel
        +models() Vec
    }
    class Optical {
        <<struct>>
        +source Vec~Source~
        +wavelength_unit String
        +unit_m f64
        +range Option
        +reference_temperature f64
        +axes Axes
        +is_anisotropic() bool
        +range_m() Option
        +nk(axis, lambda_m, temp_k) Result
        +index(axis, lambda_m, temp_k) Result
        +d_index_d_lambda(axis, lambda_m, temp_k) Result
        +group_index(axis, lambda_m, temp_k) Result
        +permittivity_diag(lambda_m, temp_k) ResultComplex64
        +permittivity_global(lambda_m, temp_k, euler_deg) ResultComplex64
        +more_3_methods()
    }
    class OpticsError {
        <<enumeration>>
        OutOfRange
        NotFinite
        BadWavelength
    }
    class fns_material_optics_rs {
        <<module material/optics.rs>>
        +MODEL_KEYS const
        +OPTICAL_KEYS const
    }
  }
  namespace material_reader_rs {
    class MaterialError {
        <<struct>>
        +material String
        +path String
        +message String
    }
    class Tree {
        <<enumeration>>
        Num
        List
        Uniaxial
    }
    class Section {
        <<struct>>
        +material String
        +path String
        +new(material, path, table) Self
        +err(key, message) MaterialError
        +has(key) bool
        +keys() Vec~String~
        +value(key) Option
        +string(key) Result
        +bool(key) Result
        +number(key) Result
        +more_6_methods()
    }
    class fns_material_reader_rs {
        <<module material/reader.rs>>
        ~edit_distance(a, b) usize
    }
  }
  namespace material_schema_rs {
    class fns_material_schema_rs {
        <<module material/schema.rs>>
        +material_schema() Value
        +SCHEMA_VERSION const
    }
  }
  namespace material_source_rs {
    class Source {
        <<struct>>
        +cite String
        +url Option~String~
        +doi Option~String~
        +verified Option~String~
        +text(cite) Source
        +is_verified() bool
        +flags_uncertainty() bool
    }
    class fns_material_source_rs {
        <<module material/source.rs>>
        +read_sources(s, key) Result
        +write_sources(v) Option~Value~
    }
  }
  namespace material_spec_rs {
    class Thermal {
        <<struct>>
        +source Vec~Source~
        +conductivity Option~TensorProp~
        +specific_heat Option~f64~
        +expansion Option~TensorProp~
    }
    class Electrical {
        <<struct>>
        +source Vec~Source~
        +permittivity Option~TensorProp~
        +permeability Option~f64~
        +conductivity Option~TensorProp~
    }
    class Superconducting {
        <<struct>>
        +source Vec~Source~
        +tc f64
        +gap_0 Option~f64~
        +lambda_l Option~f64~
        +xi_0 Option~f64~
        +mean_free_path Option~f64~
        +sigma_n Option~f64~
        +effective_gap_0() f64
    }
    class MaterialSpec {
        <<struct>>
        +name String
        +description Option~String~
        +source Vec~Source~
        +notes Option~String~
        +density Option~f64~
        +elastic Option~Elastic~
        +thermal Thermal
        +electrical Electrical
        +more_2_fields
        +from_toml_named(text, fallback_name) Result
        +from_table(table, fallback_name) Result
        +from_toml(text) Result
        +load(path) Result
        +to_toml() String
        +anisotropy(euler_deg) Anisotropy
        +legacy_scalars() BTreeMap
        +optical_permittivity(lambda_m, temp_k, euler_deg) ResultComplex64
    }
    class Library {
        <<struct>>
        +load_dir(dir) Result
        +read(text, fallback_name) Result
        +with_texts(texts) Result
        +read_table(text, fallback_name) Result
        +get(name) Option
        +require(name) Result
        +names() Vec~String~
        +len() usize
        +more_1_methods()
    }
    class Keys {
        <<struct>>
        +top slice
        +elastic slice
        +thermal slice
        +electrical slice
        +optical slice
        +model slice
        +superconducting slice
    }
    class fns_material_spec_rs {
        <<module material/spec.rs>>
        +merge_over(base, over) Table
        +keys() Keys
        +FORMAT const
        +SUPERCONDUCTING_KEYS const
    }
  }
  namespace material_tensor_rs {
    class TensorProp {
        <<enumeration>>
        Scalar
        Diagonal
        Uniaxial
        Symmetric
        +matrix() M3
        +principal_values() slice
        +as_scalar() Option~f64~
        +to_legacy() Tensor
        +check(req, label) Result
        +from_tree(s, key, tree) Result
        +read(s, key, dim, req, label) Result
        +write(unit) Value
    }
    class Requirement {
        <<enumeration>>
        PositiveDefinite
        NonNegative
        SymmetricOnly
    }
  }
  namespace material_units_rs {
    class Dim {
        <<enumeration>>
        Dimensionless
        Pressure
        Density
        Length
        Temperature
        ThermalConductivity
        SpecificHeat
        PerKelvin
        Compliance
        ElectricalConductivity
        Frequency
        Energy
        +describe() str
    }
    class fns_material_units_rs {
        <<module material/units.rs>>
        +lookup(unit) Option
        +to_si(value, unit, expect) Result
        +split_quantity(s) Result
    }
  }
  Class --> Unique
  Elastic --> Class
  Elastic --> Source
  Anisotropy --> Stiffness
  Anisotropy --> Tensor
  Dispersion --> Interp
  Dispersion --> Extrapolate
  AxisModel --> Dispersion
  AxisModel --> Source
  Axes --> AxisModel
  Optical --> Source
  Optical --> Axes
  Thermal --> Source
  Thermal --> TensorProp
  Electrical --> Source
  Electrical --> TensorProp
  Superconducting --> Source
  MaterialSpec --> Source
  MaterialSpec --> Elastic
  MaterialSpec --> Thermal
  MaterialSpec --> Electrical
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| const | `K_B` | `material/bcs.rs:18` | `f64` —  Boltzmann constant `k_B` in J/K (SI exact). |
| const | `HBAR` | `material/bcs.rs:20` | `f64` —  Reduced Planck constant `hbar` in J*s. |
| const | `EV_TO_J` | `material/bcs.rs:22` | `f64` —  Elementary charge / eV-to-Joule conversion factor in J/eV (SI exact). |
| const | `MU0` | `material/bcs.rs:24` | `f64` —  Vacuum magnetic permeability `mu_0` in H/m. |
| const | `BCS_GAP_RATIO` | `material/bcs.rs:26` | `f64` —  Default weak-coupling BCS zero-temperature gap ratio `Delta(0) / (k_B T_c)`. |
| enum | `FilmModel` | `material/bcs.rs:75` | variants: `TwoSided`, `Kerr`, `ThinLimit` —  Thin-film boundary model mapping bulk complex conductivity `sigma` and film thickness `d` |
| struct | `BcsEval` | `material/bcs.rs:106` | fields: `sigma`, `sigma_n`, `delta_j`, `lambda_eff`, `z_s`, `r_s`, `l_k_sq`, `l_eff`, `g_eff` —  Result of a BCS complex conductivity and thin-film surface impedance evaluation. |
| fn | `bcs_gap` | `material/bcs.rs:173` | `pub fn bcs_gap(tc: f64, gap_0_j: Option<f64>, temperature: f64) -> Result<f64, String>` —  Solves the exact self-consistent BCS gap equation for `Delta(T)` in Joules. |
| fn | `surface_impedance_from_params` | `material/bcs.rs:488` | `pub fn surface_impedance_from_params( tc: f64, gap_0_j: Option<f64>, lambda_l: Option<f64>, xi_0: Option<f64>, mean_free_path: Option<f64>, sigma_n: Option<f64>, freq_hz: f64, temperature: f64, thickness: f64, film_model: FilmModel, ) -> Result<BcsEval, String>` —  Evaluates arbitrary-purity BCS complex conductivity and thin-film surface impedance from explicit parameters. |
| fn | `bcs_surface_impedance` | `material/bcs.rs:616` | `pub fn bcs_surface_impedance( sc: &Superconducting, freq_hz: f64, temperature: f64, thickness: f64, film_model: FilmModel, ) -> Result<BcsEval, String>` —  Evaluates arbitrary-purity BCS complex conductivity and thin-film surface impedance for a `[superconducting]` block. |
| fn | `FilmModel::parse` | `material/bcs.rs:86` | `pub fn parse(s: &str) -> Option<FilmModel>` |
| fn | `FilmModel::as_str` | `material/bcs.rs:95` | `pub fn as_str(self) -> &'static str` |
| enum | `Unique` | `material/elastic.rs:14` | variants: `Y`, `Z` |
| enum | `Class` | `material/elastic.rs:21` | variants: `Isotropic`, `Cubic`, `Hexagonal`, `Trigonal`, `Tetragonal`, `Orthotropic`, `Monoclinic`, `Triclinic` —  Crystal class of the stiffness (the elastic symmetry group, not the full point group). |
| const | `NAMES` | `material/elastic.rs:54` | `[[&str` |
| fn | `sym_eigenvalues` | `material/elastic.rs:230` | `pub(super) fn sym_eigenvalues<const N: usize>(mut a: [[f64; N]; N]) -> [f64` —  Eigenvalues of a symmetric matrix by cyclic Jacobi rotations, ascending. |
| fn | `invert6` | `material/elastic.rs:276` | `pub fn invert6(m: &M6) -> Option<M6>` —  Inverse of a 6x6 matrix by Gauss-Jordan elimination with partial pivoting; `None` if singular. |
| struct | `Elastic` | `material/elastic.rs:329` | fields: `c`, `class`, `exact_constants`, `source` |
| fn | `Class::name` | `material/elastic.rs:85` | `pub fn name(self) -> &'static str` |
| fn | `Class::parse` | `material/elastic.rs:98` | `pub fn parse(s: &str, unique: Unique) -> Option<Class>` |
| fn | `Class::constants` | `material/elastic.rs:180` | `pub fn constants(self) -> Vec<&'static str>` —  Names of the independent constants, in the order of the class definition. |
| fn | `Elastic::from_matrix` | `material/elastic.rs:342` | `pub fn from_matrix(c: M6) -> Result<Elastic, String>` |
| fn | `Elastic::from_constants` | `material/elastic.rs:348` | `pub fn from_constants(class: Class, constants: &[(&str, f64)]) -> Result<Elastic, String>` —  Class and independent constants in Pa, for example `(Class::Cubic, [("c11", 165.7e9), ...])`. |
| fn | `Elastic::isotropic` | `material/elastic.rs:361` | `pub fn isotropic(e: f64, nu: f64) -> Result<Elastic, String>` —  Isotropic from Young's modulus and Poisson's ratio. |
| fn | `Elastic::compliance` | `material/elastic.rs:373` | `pub fn compliance(&self) -> M6` |
| fn | `Elastic::detect_class` | `material/elastic.rs:378` | `pub fn detect_class(&self) -> Class` —  The highest-symmetry class whose constants reproduce the matrix (to 1e-6 of its largest entry). |
| fn | `Elastic::global` | `material/elastic.rs:392` | `pub fn global(&self, euler_deg: [f64; 3]) -> M6` —  Stiffness in the solver's Voigt order (xx, yy, zz, xy, yz, xz) after rotating the crystal with Bunge angles (degrees). |
| fn | `Elastic::eigenvalue_range` | `material/elastic.rs:397` | `pub fn eigenvalue_range(&self) -> (f64, f64)` —  Smallest and largest eigenvalue of the stiffness in the normalised (Mandel) form, Pa. |
| fn | `Elastic::read` | `material/elastic.rs:403` | `pub fn read(s: &mut Section) -> Result<Elastic, MaterialError>` —  Reads the `[elastic]` table. |
| fn | `Elastic::write` | `material/elastic.rs:566` | `pub fn write(&self) -> toml::Table` —  The `[elastic]` table written from this value: constants if the file gave them, otherwise the matrix (upper triangle), in GPa. |
| type | `M3` | `material/mod.rs:35` |  |
| type | `M6` | `material/mod.rs:36` |  |
| enum | `Stiffness` | `material/mod.rs:40` | variants: `Cubic`, `Hexagonal`, `Trigonal`, `Orthotropic`, `General` |
| enum | `Tensor` | `material/mod.rs:151` | variants: `Uniaxial`, `Diagonal`, `Symmetric` —  A symmetric second-rank tensor property in the crystal frame. |
| struct | `Anisotropy` | `material/mod.rs:183` | fields: `stiffness`, `conductivity`, `permittivity`, `expansion`, `euler_deg` |
| fn | `rotation` | `material/mod.rs:192` | `pub fn rotation(euler_deg: [f64; 3]) -> M3` |
| fn | `euler_from_rotation` | `material/mod.rs:207` | `pub fn euler_from_rotation(m: &M3) -> [f64` —  Bunge angles (degrees) of a rotation matrix whose columns are the crystal axes in the global |
| fn | `rotate_tensor` | `material/mod.rs:231` | `pub fn rotate_tensor(t: &M3, r: &M3) -> M3` |
| fn | `stiffness_global` | `material/mod.rs:284` | `pub fn stiffness_global(s: &Stiffness, euler_deg: [f64; 3]) -> M6` —  Stiffness in the global frame (solver Voigt order). |
| fn | `isotropic_stiffness` | `material/mod.rs:289` | `pub fn isotropic_stiffness(e: f64, nu: f64) -> M6` —  Isotropic stiffness (solver Voigt order). |
| fn | `tensor_global` | `material/mod.rs:303` | `pub fn tensor_global(t: &Tensor, euler_deg: [f64; 3]) -> M3` |
| fn | `iso_tensor` | `material/mod.rs:307` | `pub fn iso_tensor(v: f64) -> M3` |
| fn | `is_positive_definite` | `material/mod.rs:312` | `pub fn is_positive_definite(c: &M6) -> bool` —  Positive definiteness check by Cholesky on the 6x6 matrix. |
| fn | `sapphire` | `material/mod.rs:334` | `pub fn sapphire() -> (Stiffness, Tensor, Tensor, Tensor, f64)` —  Typical room-temperature crystal data for the material library (GPa converted to Pa). |
| fn | `silicon` | `material/mod.rs:344` | `pub fn silicon() -> Stiffness` |
| fn | `Stiffness::name` | `material/mod.rs:52` | `pub fn name(&self) -> &'static str` |
| fn | `Stiffness::constants_mut` | `material/mod.rs:126` | `pub fn constants_mut(&mut self) -> Vec<(&'static str, &mut f64)>` —  The constants as `(label, value)` for editing. |
| fn | `Tensor::name` | `material/mod.rs:159` | `pub fn name(&self) -> &'static str` |
| fn | `Tensor::matrix` | `material/mod.rs:166` | `pub fn matrix(&self) -> M3` |
| fn | `Tensor::constants_mut` | `material/mod.rs:173` | `pub fn constants_mut(&mut self) -> Vec<(&'static str, &mut f64)>` |
| enum | `Interp` | `material/optics.rs:17` | variants: `Pchip`, `Linear` |
| enum | `Extrapolate` | `material/optics.rs:23` | variants: `Error`, `Linear` |
| enum | `Dispersion` | `material/optics.rs:32` | variants: `Constant`, `Sellmeier`, `Cauchy`, `Tabulated`, `Schott`, `Herzberger`, `Oscillators` —  One dispersion model. Parameters are in the file's wavelength unit (`x` below is the wavelength in that unit). |
| struct | `AxisModel` | `material/optics.rs:67` | fields: `dispersion`, `source`, `dn_dt` |
| enum | `Axes` | `material/optics.rs:77` | variants: `Isotropic`, `Uniaxial`, `Biaxial` |
| struct | `Optical` | `material/optics.rs:85` | fields: `source`, `wavelength_unit`, `unit_m`, `range`, `reference_temperature`, `axes` |
| enum | `OpticsError` | `material/optics.rs:101` | variants: `OutOfRange`, `NotFinite`, `BadWavelength` |
| const | `MODEL_KEYS` | `material/optics.rs:590` | `&[&str]` |
| const | `OPTICAL_KEYS` | `material/optics.rs:847` | `&[&str]` —  Keys of the `[optical]` table itself (the model keys of an isotropic material live there too). |
| fn | `Dispersion::name` | `material/optics.rs:49` | `pub fn name(&self) -> &'static str` |
| fn | `Axes::model` | `material/optics.rs:297` | `pub fn model(&self, axis: usize) -> &AxisModel` |
| fn | `Axes::models` | `material/optics.rs:310` | `pub fn models(&self) -> Vec<&AxisModel>` |
| fn | `Optical::is_anisotropic` | `material/optics.rs:320` | `pub fn is_anisotropic(&self) -> bool` |
| fn | `Optical::range_m` | `material/optics.rs:325` | `pub fn range_m(&self) -> Option<(f64, f64)>` —  The validity range in metres. |
| fn | `Optical::nk` | `material/optics.rs:330` | `pub fn nk(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<(f64, f64, f64, f64), OpticsError>` —  `(n, k, dn/dlambda, dk/dlambda)` for crystal axis 0, 1 or 2 (derivatives per metre). |
| fn | `Optical::index` | `material/optics.rs:359` | `pub fn index(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<Complex64, OpticsError>` —  Complex index `n - j k`. |
| fn | `Optical::d_index_d_lambda` | `material/optics.rs:365` | `pub fn d_index_d_lambda(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<Complex64, OpticsError>` —  `d(n - j k)/d lambda` per metre. |
| fn | `Optical::group_index` | `material/optics.rs:371` | `pub fn group_index(&self, axis: usize, lambda_m: f64, temp_k: Option<f64>) -> Result<f64, OpticsError>` —  Group index `n - lambda dn/dlambda` (real part). |
| fn | `Optical::permittivity_diag` | `material/optics.rs:377` | `pub fn permittivity_diag(&self, lambda_m: f64, temp_k: Option<f64>) -> Result<[Complex64` —  Relative permittivity along the three crystal axes, `(n - j k)^2`. |
| fn | `Optical::permittivity_global` | `material/optics.rs:387` | `pub fn permittivity_global(&self, lambda_m: f64, temp_k: Option<f64>, euler_deg: [f64; 3]) -> Result<[[Complex64` —  Relative permittivity tensor in the global frame for a crystal rotated by the Bunge angles (degrees). |
| fn | `Optical::read` | `material/optics.rs:403` | `pub fn read(s: &mut Section) -> Result<Optical, MaterialError>` |
| fn | `Optical::validate` | `material/optics.rs:477` | `pub fn validate(&self) -> Result<(), (String, String)>` —  Checks that every model is usable over the declared range. |
| fn | `Optical::write` | `material/optics.rs:544` | `pub fn write(&self) -> toml::Table` |
| struct | `MaterialError` | `material/reader.rs:12` | fields: `material`, `path`, `message` —  A problem found in a material file or description, in plain words. |
| enum | `Tree` | `material/reader.rs:33` | variants: `Num`, `List`, `Uniaxial` —  A quantity tree after unit conversion: numbers in SI, lists, or the uniaxial pair. |
| struct | `Section` | `material/reader.rs:39` | fields: `material`, `path` |
| fn | `edit_distance` | `material/reader.rs:46` | `pub(crate) fn edit_distance(a: &str, b: &str) -> usize` |
| fn | `Section::new` | `material/reader.rs:61` | `pub fn new(material: &str, path: &str, table: &'a toml::Table) -> Self` |
| fn | `Section::err` | `material/reader.rs:65` | `pub fn err(&self, key: &str, message: impl Into<String>) -> MaterialError` |
| fn | `Section::has` | `material/reader.rs:76` | `pub fn has(&self, key: &str) -> bool` |
| fn | `Section::keys` | `material/reader.rs:80` | `pub fn keys(&self) -> Vec<String>` |
| fn | `Section::value` | `material/reader.rs:84` | `pub fn value(&mut self, key: &str) -> Option<&'a toml::Value>` |
| fn | `Section::string` | `material/reader.rs:92` | `pub fn string(&mut self, key: &str) -> Result<Option<String>, MaterialError>` |
| fn | `Section::bool` | `material/reader.rs:100` | `pub fn bool(&mut self, key: &str) -> Result<Option<bool>, MaterialError>` |
| fn | `Section::number` | `material/reader.rs:117` | `pub fn number(&mut self, key: &str) -> Result<Option<f64>, MaterialError>` —  A plain number (no unit); a string such as "0.17" is accepted too. |
| fn | `Section::number_of` | `material/reader.rs:124` | `pub fn number_of(&self, key: &str, v: &toml::Value) -> Result<f64, MaterialError>` |
| fn | `Section::numbers` | `material/reader.rs:135` | `pub fn numbers(&mut self, key: &str) -> Result<Option<Vec<f64>>, MaterialError>` |
| fn | `Section::quantities` | `material/reader.rs:209` | `pub fn quantities(&mut self, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<Tree>, MaterialError>` —  A quantity tree (number, list, or uniaxial pair) converted to SI. |
| fn | `Section::quantity` | `material/reader.rs:216` | `pub fn quantity(&mut self, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<f64>, MaterialError>` —  One quantity in SI. |
| fn | `Section::sub` | `material/reader.rs:224` | `pub fn sub(&mut self, key: &str) -> Result<Option<Section<'a>>, MaterialError>` |
| fn | `Section::finish` | `material/reader.rs:236` | `pub fn finish(&self, known: &[&str]) -> Result<(), MaterialError>` —  Reports the first key nobody read, with the closest known key as a suggestion. |
| const | `SCHEMA_VERSION` | `material/schema.rs:9` | `&str` |
| fn | `material_schema` | `material/schema.rs:68` | `pub fn material_schema() -> Value` —  The JSON Schema (draft 2020-12) of format 1. |
| struct | `Source` | `material/source.rs:11` | fields: `cite`, `url`, `doi`, `verified` |
| fn | `read_sources` | `material/source.rs:93` | `pub fn read_sources(s: &mut Section, key: &str) -> Result<Vec<Source>, MaterialError>` —  Reads `key` of the section as a list of sources (empty if the key is absent). |
| fn | `write_sources` | `material/source.rs:102` | `pub fn write_sources(v: &[Source]) -> Option<toml::Value>` —  The TOML value for a list of sources, or `None` if it is empty. |
| fn | `Source::text` | `material/source.rs:20` | `pub fn text(cite: &str) -> Source` |
| fn | `Source::is_verified` | `material/source.rs:25` | `pub fn is_verified(&self) -> bool` —  True if the entry has a verification date. |
| fn | `Source::flags_uncertainty` | `material/source.rs:30` | `pub fn flags_uncertainty(&self) -> bool` —  True if the citation text itself flags the data as unchecked. |
| const | `FORMAT` | `material/spec.rs:20` | `i64` —  The version of the material file format this code reads and writes. |
| struct | `Thermal` | `material/spec.rs:23` | fields: `source`, `conductivity`, `specific_heat`, `expansion` |
| struct | `Electrical` | `material/spec.rs:35` | fields: `source`, `permittivity`, `permeability`, `conductivity` |
| struct | `Superconducting` | `material/spec.rs:47` | fields: `source`, `tc`, `gap_0`, `lambda_l`, `xi_0`, `mean_free_path`, `sigma_n` |
| struct | `MaterialSpec` | `material/spec.rs:73` | fields: `name`, `description`, `source`, `notes`, `density`, `elastic`, `thermal`, `electrical`, `optical`, `superconducting` |
| const | `SUPERCONDUCTING_KEYS` | `material/spec.rs:93` | `&[&str]` |
| struct | `Library` | `material/spec.rs:367` |  A folder of material files, loaded at run time. Materials are found by file name (without `.toml`), ignoring case. |
| fn | `merge_over` | `material/spec.rs:556` | `pub fn merge_over(base: &toml::Table, over: &toml::Table) -> toml::Table` —  `over` laid over `base`: every top-level key of `over` replaces the one of `base`, and a table (`elastic`, `thermal`, `electrical`, |
| struct | `Keys` | `material/spec.rs:573` | fields: `top`, `elastic`, `thermal`, `electrical`, `optical`, `model`, `superconducting` —  The keys each table of a material file accepts (the contract test compares them with the JSON Schema). |
| fn | `keys` | `material/spec.rs:583` | `pub fn keys() -> Keys` |
| fn | `Superconducting::effective_gap_0` | `material/spec.rs:67` | `pub fn effective_gap_0(&self) -> f64` —  Zero-temperature BCS energy gap Delta(0) in J: `gap_0` if given, otherwise the weak-coupling |
| fn | `MaterialSpec::from_toml_named` | `material/spec.rs:118` | `pub fn from_toml_named(text: &str, fallback_name: &str) -> Result<MaterialSpec, MaterialError>` —  Reads a material from TOML text. `fallback_name` is used when the file has no `name`. |
| fn | `MaterialSpec::from_table` | `material/spec.rs:128` | `pub fn from_table(table: &toml::Table, fallback_name: &str) -> Result<MaterialSpec, MaterialError>` —  Reads a material from a parsed table (already merged with its base if it had `extends`). |
| fn | `MaterialSpec::from_toml` | `material/spec.rs:218` | `pub fn from_toml(text: &str) -> Result<MaterialSpec, MaterialError>` |
| fn | `MaterialSpec::load` | `material/spec.rs:223` | `pub fn load(path: &Path) -> Result<MaterialSpec, MaterialError>` —  Reads a file; the file name (without `.toml`) is the name if the file has none. |
| fn | `MaterialSpec::to_toml` | `material/spec.rs:235` | `pub fn to_toml(&self) -> String` —  The material as TOML text, in the same format (comments of an original file are not kept). |
| fn | `MaterialSpec::anisotropy` | `material/spec.rs:310` | `pub fn anisotropy(&self, euler_deg: [f64; 3]) -> Anisotropy` —  The material in the older `Anisotropy` form (stiffness as a general matrix, tensors), placed with Bunge angles in degrees. |
| fn | `MaterialSpec::legacy_scalars` | `material/spec.rs:321` | `pub fn legacy_scalars(&self) -> BTreeMap<String, f64>` —  Scalar properties under the older names (`rho`, `Cp`, `mu_r`, and `eps_r`, `k`, `alpha`, `E`, `nu` where the data are isotropic). |
| fn | `MaterialSpec::optical_permittivity` | `material/spec.rs:345` | `pub fn optical_permittivity(&self, lambda_m: f64, temp_k: Option<f64>, euler_deg: [f64; 3]) -> Result<[[Complex64` —  Relative permittivity tensor at a wavelength (m) and temperature (K) in the frame of a crystal rotated by Bunge angles. |
| fn | `Library::load_dir` | `material/spec.rs:385` | `pub fn load_dir(dir: &Path) -> Result<Library, Vec<MaterialError>>` —  Loads every `*.toml` file of the folder (not its subfolders). All files with problems are reported together. |
| fn | `Library::read` | `material/spec.rs:440` | `pub fn read(&self, text: &str, fallback_name: &str) -> Result<MaterialSpec, MaterialError>` —  Reads material text that may use `extends`, looking the base up in this library (for a file outside the library folder, |
| fn | `Library::with_texts` | `material/spec.rs:446` | `pub fn with_texts(&self, texts: &BTreeMap<String, String>) -> Result<Library, Vec<MaterialError>>` —  A copy of the library with more materials given as TOML text by name (the materials a project carries). A text may extend a material |
| fn | `Library::read_table` | `material/spec.rs:467` | `pub fn read_table(&self, text: &str, fallback_name: &str) -> Result<(MaterialSpec, toml::Table), MaterialError>` —  `read`, with the table after `extends` was resolved. |
| fn | `Library::get` | `material/spec.rs:485` | `pub fn get(&self, name: &str) -> Option<&MaterialSpec>` |
| fn | `Library::require` | `material/spec.rs:490` | `pub fn require(&self, name: &str) -> Result<&MaterialSpec, MaterialError>` —  Like `get`, with an error that lists the available names. |
| fn | `Library::names` | `material/spec.rs:498` | `pub fn names(&self) -> Vec<String>` |
| fn | `Library::len` | `material/spec.rs:502` | `pub fn len(&self) -> usize` |
| fn | `Library::is_empty` | `material/spec.rs:506` | `pub fn is_empty(&self) -> bool` |
| enum | `TensorProp` | `material/tensor.rs:15` | variants: `Scalar`, `Diagonal`, `Uniaxial`, `Symmetric` |
| enum | `Requirement` | `material/tensor.rs:26` | variants: `PositiveDefinite`, `NonNegative`, `SymmetricOnly` —  What a physical tensor must satisfy. |
| fn | `TensorProp::matrix` | `material/tensor.rs:36` | `pub fn matrix(&self) -> M3` |
| fn | `TensorProp::principal_values` | `material/tensor.rs:45` | `pub fn principal_values(&self) -> [f64` |
| fn | `TensorProp::as_scalar` | `material/tensor.rs:50` | `pub fn as_scalar(&self) -> Option<f64>` —  The scalar if the tensor is isotropic. |
| fn | `TensorProp::to_legacy` | `material/tensor.rs:58` | `pub fn to_legacy(&self) -> super::Tensor` —  The tensor in the older `Anisotropy` form. |
| fn | `TensorProp::check` | `material/tensor.rs:68` | `pub fn check(&self, req: Requirement, label: &str) -> Result<(), String>` —  Checks the requirement; `label` names the property in the message. |
| fn | `TensorProp::from_tree` | `material/tensor.rs:85` | `pub fn from_tree(s: &Section, key: &str, tree: &Tree) -> Result<TensorProp, MaterialError>` |
| fn | `TensorProp::read` | `material/tensor.rs:126` | `pub fn read(s: &mut Section, key: &str, dim: Dim, req: Requirement, label: &str) -> Result<Option<TensorProp>, MaterialError>` —  Reads `key` of the section, converts to SI, checks the requirement. |
| fn | `TensorProp::write` | `material/tensor.rs:134` | `pub fn write(&self, unit: f64) -> toml::Value` —  The value as TOML, numbers divided by `unit` (SI factor of the unit written next to it). |
| enum | `Dim` | `material/units.rs:9` | variants: `Dimensionless`, `Pressure`, `Density`, `Length`, `Temperature`, `ThermalConductivity`, `SpecificHeat`, `PerKelvin`, `Compliance`, `ElectricalConductivity`, `Frequency`, `Energy` —  The kind of quantity an entry holds; a unit of another kind is refused. |
| fn | `lookup` | `material/units.rs:85` | `pub fn lookup(unit: &str) -> Option<(Dim, f64)>` —  Dimension and factor to SI of a unit name. `None` if the unit is unknown. |
| fn | `to_si` | `material/units.rs:142` | `pub fn to_si(value: f64, unit: &str, expect: Dim) -> Result<f64, String>` —  `value` in `unit` converted to SI, if `unit` is a unit of dimension `expect`. |
| fn | `split_quantity` | `material/units.rs:160` | `pub fn split_quantity(s: &str) -> Result<(f64, String), String>` —  Splits `"496.8 GPa"` (or `"496.8GPa"`) into the number and the unit text. |
| fn | `Dim::describe` | `material/units.rs:30` | `pub fn describe(self) -> &'static str` |

Trait implementations in this component: `OpticsError: Display`, `OpticsError: Error`, `MaterialError: Display`, `MaterialError: Error`
