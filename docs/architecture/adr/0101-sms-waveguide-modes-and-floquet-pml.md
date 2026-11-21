# 101. Solid waveguide modes, Floquet-periodic absorbing layers, and lossless convergence floor

Date: 2026-11-21

## Status

Accepted. Amends ADR 42 (`nyx-sms`), ADR 61 (`convergence_check`), ADR 68 and 80 (cross-section solid studies and mode identification scalars), ADR 35 (material library), and ADR 97 (`IM_FLOOR` in `convergence_check`).

## Context

Four coupled issues surfaced when modelling phononic and optomechanical waveguides (such as a 300 nm silicon strip on R-plane sapphire at longitudinal wavevector $q \approx 2.35 \times 10^7\ \text{rad/m}$):

1. **Floquet-periodic phase dropped by the `sms` PML kernel.** In `nyx-sms::solid::solve_variant`, any problem with absorbing layers (`spec.materials.has_pml()`) dispatches immediately to `pml::solve_solid_pml`, bypassing `floquet::solve_floquet`. However, `pml::assemble_solid_pml` reduced `K` and `M` only with `free_dofs(n, &st.fixed)` and never consulted `spec.periodic`. On a 3D 1-layer swept cell with lateral PMLs and a `ymin -> ymax` Floquet pair (`phase = q * L_y`), the periodic boundaries were left open and free, making the eigenfrequencies independent of `phase`. In addition, `assemble_solid_pml` called `spec.materials.cell_values(mesh, "rho", None)`, which failed when an `inactive` region (such as `air`) had no mass density `"rho"`.
2. **No native 2.5D elastic waveguide study on 2D cross-sections.** While `nyx-ems` provides a 2D `modes` study (`exp(-j beta y)`) and `nyx-sms` provides `axisymmetric_eigenfrequency` (`exp(-j m theta)` on the `(r, z)` half-plane via `cross_section.rs`), finding the elastic modes of a translationally invariant waveguide required extruding the 2D cross-section into a 1-layer 3D swept mesh with a Floquet pair. The 3D swept mesh doubles the nodal DOFs (both `y = 0` and `y = L_y` planes before reduction), introduces a finite-difference approximation $O((q L_y)^2)$ along $y$, and requires 3D structured PMLs.
3. **Coupled-beta fixed-point rounding on lossless tilted crystals exceeded `IM_FLOOR`.** On a lossless waveguide with a tilted anisotropic crystal (such as R-plane sapphire, `euler_deg = [0.0, 57.6, 0.0]`), `nyx-ems` iterates the coupled longitudinal propagation constant until $|\Delta \beta / \beta| < \text{COUPLED\_BETA\_TOL} = 10^{-8}$. That fixed-point tolerance leaves a residual numerical $|\text{Im}(n_{\text{eff}}) / \text{Re}(n_{\text{eff}})|$ between $10^{-9}$ and $10^{-8}$ ($Q \sim 10^8$), just above the previous `IM_FLOOR = 1e-9` ($Q \ge 5 \times 10^8$) in `nyx-engine::convergence`, causing `convergence_check` to fail on the unflagged `scalar Q 1` noise ratio. Moreover, `noise_modes` only inspected `eigenvalue[i].im`, missing solvers (`nyx-sms::pml`) that report `scalar Q<n>` without `eigenvalue`.
4. **Missing elastic tables on `lithium_niobate` and `silicon_nitride`.** Both shipped material files had `[optical]` dispersion tables only, preventing optomechanical or acoustic simulations with `set_material`.

## Decision

1. **Floquet-periodic congruence reduction and inactive mass skip in `nyx-sms::pml`.**
   * `assemble_solid_pml` skips domains listed in `spec.inactive` when looking up `"rho"` (assigning `0.0`, matching `solid::cell_rho`), and requires `"rho"` on active domains only.
   * When `spec.periodic.is_empty()`, `assemble_solid_pml`, `solve_pml_modes`, and `pml_fields` keep the exact `free_dofs` + `submatrix(&free, &free)` path so all non-periodic PML results and goldens remain bit-identical.
   * When `!spec.periodic.is_empty()`, `assemble_solid_pml` matches periodic boundary vertices with `femcore::periodic::{match_vertices, nodal_relations}`, builds `Reduction::build(n, &relations, &st.fixed)`, and reduces the complex stiffness and mass matrices via `congruence(&pm, &kfull)` and `congruence(&pm, &mfull)` ($P^H K P$ and $P^H M P$). Because `congruence` conjugates the prolongation entries on the row (test) side (`exp(+i phase)`) without conjugating the interior matrix entries, the PML coordinate stretches remain complex-symmetric while the Bloch phase acts hermitianly across the periodic pair. `solve_pml_modes` filters unconverged Krylov vectors (`residual < 1e-6`), and `pml_fields` expands reduced eigenvectors via `red.expand(...)`.

2. **Native 2.5D elastic waveguide `modes` study in `nyx-sms`.**
   * `crates/nyx-sms/src/waveguide.rs` provides `assemble_stiffness_waveguide(mesh, c6, order, q, stretch)` and `mode_energy_waveguide(mesh, rho, order, class, n_classes, modes)` on a 2D `(x, z)` cross-section mesh with three displacement component blocks `(0: u_x, 1: u_z, 2: u_y)` matching the cross-section Voigt permutation `(xx, zz, yy, xz, zy, xy)` of `femcore::section::plane_stiffness(c, 2)` (ADR 1).
   * With displacement $\mathbf{u}(x, z)\,e^{-i q y}$ at longitudinal wavevector $q$ ($\text{rad/m}$), the strain operator splits into $B = B_0 + (-i q) B_1$ on the trial side and $B_0 + (+i q) B_1$ on the test side; inside a 2D PML (`femcore::pml::cell_stretches`), in-plane derivatives carry $1/s_x$ and $1/s_z$ with Jacobian factor $s_x s_z$ unconjugated on both sides.
   * `CrossSection::Waveguide(f64)` extends `cross_section.rs` alongside `CrossSection::Axisymmetric(f64)`. Unlike the axisymmetric model, `Waveguide` preserves full 3D anisotropic crystal orientation without rotational averaging (`average_about_axis`), imposes no $r = 0$ axis constraints, uses planar Cartesian mass matrices and PML stretches, and supports lateral 1D Floquet-periodic boundary pairs via `Reduction::build` and `congruence`.
   * `Sms::studies()` exposes `kind: "modes"` ("Waveguide modes", `Dims::Two`, parameters `n_modes`, `f_search`, `wavevector` in $\text{rad/m}$), with fields `|u|`, `Re ux`, `Re uz`, `Re uy`, `|uy|` and per-mode identification scalars `f re N (Hz)`, `f im N (Hz)`, `fN (Hz)`, `QN` (when lossy), `layer fraction N`, `lateral share N`, `vertical share N`, `longitudinal share N`, and `share in domain <label> N` (renamed to `share in <region> N` by the engine).

3. **Lossless convergence imaginary floor and `scalar Q` noise detection in `nyx-engine::convergence`.**
   * `IM_FLOOR` is raised from `1e-9` to `1e-8` ($Q_{\text{floor}} = 1 / (2 \cdot \text{IM\_FLOOR}) = 5 \times 10^7$), covering the `COUPLED_BETA_TOL = 1e-8` fixed-point termination tolerance of tilted-crystal optical waveguide solves while keeping any physical loss ($Q \le 10^7$) strictly counted.
   * `noise_modes` also flags mode `i` when `scalar Q<n>` (or `scalar Q <n> ...`) is present without `eigenvalue[i].im` and exceeds `q_floor` on both the coarse and refined meshes.

4. **Elastic properties in `materials/lithium_niobate.toml` and `materials/silicon_nitride.toml`.**
   * `lithium_niobate.toml` gains `density = "4640 kg/m^3"` and trigonal (class `3m`) constant-$E$ `[elastic]` constants (`c11 = 203.0`, `c12 = 57.3`, `c13 = 75.2`, `c14 = 8.5`, `c33 = 242.4`, `c44 = 59.5` GPa; Weis and Gaylord, *Appl. Phys. A* 37, 191–203, 1985).
   * `silicon_nitride.toml` gains `density = "3100 kg/m^3"` and approximate stoichiometric LPCVD thin-film isotropic `[elastic]` constants (`E = 290.0` GPa, `nu = 0.27`). `femcore::material::Elastic::read` canonically rounds isotropic derived $(c_{11}, c_{12})$ to their GPa representation `(x / 1e9) * 1e9` so `MaterialSpec::to_toml()` round-trips bit-identically when sent to a worker process.

## Consequences

* 3D swept Floquet cells with absorbing layers honour the Bloch phase advance `phase` and allow inactive domains without `"rho"`.
* Phononic and acousto-optic waveguide cross-sections can be solved directly in 2D with `study = "modes"` and `wavevector = q` (with or without 2D `add_pml`), running roughly $3\times$ to $5\times$ faster than a 3D swept cell and avoiding longitudinal finite-difference error.
* Lossless tilted-crystal optical waveguides pass `convergence_check` with coupled-beta rounding automatically marked as `noise: true`.
* Contracts: `physics_sms` 1.5, `physics_sms_result` 1.3, `kernel_stages` 1.6 (all additive).

## Code

`crates/nyx-sms/src/{pml,waveguide,cross_section,problem,lib,solid}.rs`, `crates/nyx-sms/tests/{pml_floquet,waveguide_modes,contracts,stages}.rs`, `crates/nyx-engine/src/convergence.rs`, `crates/nyx-engine/tests/validated.rs`, `crates/femcore/src/material/elastic.rs`, `materials/{lithium_niobate,silicon_nitride}.toml`, `crates/femcore/tests/material_library.rs`, `docs/MATERIALS.md`, `docs/architecture/contracts/{physics_sms,physics_sms_result}.v1.json`, `docs/architecture/contracts/VERSIONS.md`.

## Contract

`physics_sms` 1.5, `physics_sms_result` 1.3, `kernel_stages` 1.6.
