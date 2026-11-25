# 129. Coupled piezoelectric FEM in `nyx-sms`

Date: 2026-11-25

## Status

Accepted. Extends ADR 42 (`nyx-sms` architecture), ADR 43 & ADR 68 & ADR 69 (axisymmetric solid mechanics), ADR 101 (2.5D waveguide modes), ADR 106 (driven frequency-domain PROM & multi-port scattering), ADR 109 (material `format = 2` coupled tensors), ADR 113 (inline tensor schema), ADR 118 (`NetworkAnalysis`), and ADR 128 (`MixedDiscretisation` & curvilinear PML/IE).

## Context

Prior to this decision, `nyx-sms` solved only uncoupled linear elasticity `(u_x, u_y, u_z)`, while `femcore::material` already parsed `[piezoelectric]` TOML sections (`e` / `d`, `format = 2`). Coupled electromechanical resonator, SAW/BAW filter, actuator, and disk-resonator workflows required:

1. Exact IEEE/textbook-to-solver Voigt tensor conventions across isotropic, hexagonal `6mm`, and trigonal `3m` crystals under arbitrary Bunge Z-X-Z Euler rotations (`C^E`, `C^D`, `e`, `eps^S`, `eta`), plus thermodynamic self-consistency validation (`|C^D - C^E - e^T (eps^S)^{-1} e|_F / |C^E|_F`) on material loading.
2. Coupled `(u, phi)` 2D/3D static, eigenfrequency, 2.5D waveguide (`modes`), axisymmetric (`6mm`/`inf m`), and multi-port frequency-domain `driven` admittance solvers on both `SimplexMesh` and native `MixedDiscretisation` (`Quad`/`Prism`/`Hex`), supporting conductors, floating electrodes (`Reduction` hub DOFs), PMLs, `infinite_element` domains, and thermal expansion.
3. Extreme block scale disparity (`K_uu ~ 1e11 Pa` vs `K_pp ~ 1e-11 F/m`, `22` orders of magnitude) and indefinite saddle-point structure `[[K_uu, K_{u phi}], [K_{phi u}, -K_{phi phi}]]`, which break standard direct factorizations without equilibration and break naive algebraic multigrid without block Schur preconditioning.
4. Equivalent-circuit resonator parameter extraction (`piezo_resonator` mBVD model, `f_s`, `f_p`, `k_t^2`, `k_eff^2`, unwrapped phase-slope `Q_s`/`Q_p`).

## Decision

* **Shared constitutive, equilibration, and PROM infrastructure (`femcore`, `nyx-physics`, `nyx-engine`):**
  * Fixed solver-to-textbook Voigt conversion on material specification (`stiffness_constant_e`, `stiffness_constant_d`, non-double-permuted `viscosity`), added `piezo_e_global`, `plane_piezo_e`, `plane_piezo_e_textbook`, and `rotate_m6_textbook`, enforced `<= 5%` Frobenius thermodynamic consistency on load, cleared inherited `[dielectric]` fields on `merge_subtable` overrides, and enabled rotated `[piezoelectric]`, `[photoelastic]`, and `[viscoelastic]` evaluation in `nyx-engine::materials`.
  * Extracted 8-pass symmetric Ruiz equilibration (`femcore::linalg::equilibrate`, supporting `RCsr`/`CCsr` affine terms) and multi-port Hermite + isometric-QR PROM (`femcore::linalg::prom`) from `nyx-ems::driven` for shared reuse across `nyx-ems` and `nyx-sms`.
  * Upgraded `femcore::periodic` with `congruence_fast`, `restrict_rhs`, inhomogeneous Dirichlet lifting, and weighted DSU cycle/collision validation, and added zero-shift guard (`theta.norm() < 1e-15`) in `driver::eigs_near`.
  * Extended `nyx-physics` (`PropKind::Matrix3x6`, `PropValue::Matrix3x6([[f64; 6]; 3])`, `MaterialAssign::Oriented` at `oneOf[2]` in JSON Schema, `"piezo_e"`, `"eps_r"`, `"eta"`, and shared `"piezo"` decode guard) while skipping non-scalar `Matrix3x6` in `app::tabs` scalar UI grids.
* **2D/3D coupled static, eigenfrequency & Schur solvers (`nyx-sms::{piezo_mat,piezo_asm,piezo}`):**
  * `Sms` exposes `piezoelectric` (`"auto" | "on" | "off"`), `"eps_r"` / `"eta"` / `"piezo_e"`, and electrical BCs (`ground`, `potential`, `surface_charge`, `floating_potential`, `conductor`, `terminal`), activating coupled `(u, phi)` assembly when `[piezoelectric]` data pairs with electrical BCs/permittivity under `"auto"`, or unconditionally under `"on"`.
  * Dual-backend assembly on `SimplexMesh` and `MixedDiscretisation` handles exact 3x3 `Plane::Stress` electromechanical static condensation `(D_uu, e_pu, eps_ps)`, `PmlJacobian` / `ie_rational_jacobian` / `ie_exponential_jacobian`, piezo-stiffened Christoffel M-PML backward-wave checks, connected-component gauge BFS for ungrounded electrostatic components, and electrode hub `Reduction` DOFs.
  * Power-of-2 block equilibration + 8-pass Ruiz scaling balances `[K_uu, K_{u phi}; K_{phi u}, -K_{phi phi}]` before sparse solve. Static 3D simplex problems with `solver = "iterative"` (or `"auto"` above `120,000` DOFs) use Right-FGMRES with upper block-triangular `PiezoSchurPc` (`AmgSystem` displacement V-cycle + `M_phi = K_pp + diag(sum_k K_{pu,ik}^2 / K_{uu,kk})` solved by SPD `DirectSolver` or line-smoothed `AmgSystem` above `50,000` `phi` DOFs when hub-free). Eigenfrequency studies eliminate `phi` implicitly inside `SchurShiftInvertOp<T>` (real/complex Krylov-Schur) and `SchurQepShiftInvertOp` (Kelvin-Voigt + LRB QEP).
* **Driven multi-port admittance PROM & mBVD resonator analysis (`nyx-sms::piezo_driven`, `nyx-analysis::piezo_resonator`):**
  * Registers `Study::Driven` and `"terminal"` BCs on `Sms`, computes frequency-domain port admittance `Y(omega)` via static condensation against non-conjugate `T_PI = T_IP^T`, extracts clamped capacitance `C_0`, scattering `S`, and impedance `Z`, supports `"auto" | "direct" | "lblt" | "blr"` backends, exports `DrivenPortData`, and accelerates sweeps using shared multi-port `femcore::linalg::prom`.
  * Adds `piezo_resonator` analysis in `nyx-analysis` (`schur_reduced_admittance` for multi-port terminations, `f_s`/`f_p`, unwrapped phase-slope `Q_s`/`Q_p`, IEEE `k_t^2` / `k_eff^2`, and bounded reflection-space Levenberg–Marquardt 6-parameter mBVD fit `(C_0, R_0, R_s, C_m, L_m, R_m)`).
* **2.5D Waveguide & Axisymmetric `6mm`/`inf m` extensions (`nyx-sms::{waveguide,axisym,cross_section}`):**
  * Extends `modes` and `axisymmetric_eigenfrequency` to 4-block `(u_x, u_z, u_y, phi)` / `(u_r, u_z, u_theta, phi)` complex assemblies on both `SimplexMesh` and `MixedDiscretisation` with explicit non-conjugated `K_{phi u}`, reusing `SchurShiftInvertOp<C>` and QEP solvers. Axisymmetry enforces `<= 5%` rotational invariance around `z` (diagnostic warning `> 1e-3`, error `> 5%`), `phi = 0` regularity on `r = 0` for `|m| >= 1` (`phi` free at `m = 0`), floating electrode restriction to `m = 0`, and cylindrical `PmlJacobian` (`s_theta`).

## Consequences

* `nyx-sms` solves static, modal, waveguide, axisymmetric `6mm`, and driven multi-port piezoelectric BAW, SAW, Lamb, MEMS, and disk resonators across P1/P2 `SimplexMesh` and native `MixedDiscretisation` (`Quad`/`Prism`/`Hex`), preserving bit-for-bit regression stability on all pure-elastic golden suites.
* Large 3D static piezoelectric systems converge in $O(1)$ iterations under `PiezoSchurPc` + Right-FGMRES, while `nyx-ems::driven` and `nyx-sms::piezo_driven` share a single verified equilibration and multi-port Hermite PROM stack.

## Code

`crates/femcore/src/linalg/{equilibrate,prom,driver,mod}.rs`, `crates/femcore/src/material/{mod,spec,coupled}.rs`, `crates/femcore/src/{periodic,props,section}.rs`, `crates/femcore/src/meshing/physics_size.rs`, `crates/femcore/tests/material_library.rs`, `crates/nyx-physics/src/{desc,spec,schema,validate,decode}.rs`, `crates/nyx-physics/tests/inline_tensor.rs`, `crates/nyx-ems/src/driven.rs`, `crates/nyx-sms/src/{piezo_mat,piezo_asm,piezo,piezo_driven,axisym,waveguide,cross_section,pml,problem,solid,lib}.rs`, `crates/nyx-sms/tests/{piezo,piezo_driven,piezo_waveguide,piezo_axisym,contracts}.rs`, `crates/nyx-analysis/src/{piezo_resonator,network,lib}.rs`, `crates/nyx-engine/src/{analyse,commands,estimate,materials}.rs`, `crates/nyx-engine/tests/{materials,piezo_driven_network}.rs`, `crates/app/src/tabs.rs`, `materials/aluminum_nitride_piezomumps.toml`.

## Contract

`physics_api` 1.7, `physics_sms` 1.8 (`physics_sms.v1.json`), `physics_sms_result` 1.5 (`physics_sms_result.v1.json`), `kernel_stages` 1.8, `physics_dfs` 1.7 (`physics_dfs.v1.json`), `physics_ems` 1.8 (`physics_ems.v1.json`), `engine_commands` 1.31 (`engine_commands.v1.json`), and `docs/api/*.md`.
