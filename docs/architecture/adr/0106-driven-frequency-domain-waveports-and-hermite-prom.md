# 106. Driven frequency-domain solver, waveports, and Hermite Tangential Galerkin PROM

Date: 2026-11-21

## Status

Accepted (superconducting qubit milestone M5). Extends the 3D `ems` physics of ADR 45 and ADR 104.

## Context

Characterizing microwave transmission $S_{21}(\omega)$ and black-box quantization (BBQ) driving-point admittance $Y_{\text{JJ}}(\omega)$ across hundreds of frequency points around high-$Q$ qubit and readout resonances with a direct 3D FEM solve per frequency is prohibitively slow and fails to exploit the smooth rational/Meromorphic structure of the Maxwell operator. Moreover, exterior `waveport` boundaries carry frequency-dependent modal propagation constants $\beta_p(\omega)$ and impedances $Z_{0,p}(\omega)$ that cannot be lumped into a frequency-independent damping matrix $\mathbf{C}_r$.

## Decision

* **`Study::Driven` (`f_min`, `f_max`, `n_points`, `prom_tol`, `crates/nyx-ems/src/{lib,problem,driven}.rs`).** Adds `study = "driven"` (`Dims::Three`). Setting `prom_tol` to a positive tolerance (default `1e-3`) enables greedy PROM sweep; setting `prom_tol` to `null` or `<= 0.0` runs exact point-by-point direct HDM solves on every frequency sample. `pec`, `pmc`, `surface_impedance`, and `lumped_jj` remain fully supported inside `Study::Driven`.
* **Port conditions (`crates/nyx-ems/src/ports.rs`).** Adds `BcKind::LumpedPort { name, z0, direction }` (rank-1 voltage/current terminal) and `BcKind::WavePort { name, mode }`. Each `waveport` builds a local orthonormal 2D plane triangulation on the port face inheriting adjacent 3D cell $\varepsilon_r/\mu_r$ and conductor edges, solves 2D waveguide modes via `nyx_ems::modes`, enforces smooth eigenvector sign/phase alignment across anchor frequencies ($\operatorname{Re}\langle \mathbf{e}_{t,p}(\omega_k), \mathbf{e}_{t,p}(\omega_{\text{ref}})\rangle \ge 0$), fits $\beta_p^2(k_0^2)$ quadratically for dispersive sweeps, and assembles separate sparse 3D waveport matrices $\mathbf{W}_p$ and modal vectors $\mathbf{b}_p$ so dispersive $\sum_p i\beta_p(\omega)\mathbf{W}_{p,r}$ stays explicit outside constant $\mathbf{C}_r$.
* **Hermite Tangential Galerkin PROM (`crates/nyx-ems/src/driven.rs`).** At each greedily chosen anchor frequency $\omega_a$, `solve_driven` applies symmetric Ruiz diagonal equilibration, factorizes $\tilde{\mathbf{A}}(\omega_a)$ **once** via `DirectSolver::new(&a_anchor, false)`, solves value, 1st, and 2nd frequency derivatives ($\mathbf{X}, \mathbf{X}', \mathbf{X}''$) against that single LU factorization, augments the real orthonormal projection basis $\mathbf{Q}$ and isometric QR image-space basis $\mathbf{U}$ for microsecond full-space residual certification over all `n_points`, extracts multi-port $\mathbf{S}(\omega), \mathbf{Z}(\omega), \mathbf{Y}(\omega)$ into `DrivenPortData`, and stores representative 3D field solutions.

## Consequences

* A 500-point driven sweep across coupled 3D resonator–feedline–junction circuits requires `<= 3` HDM factorizations (`<= 12` limit), projects all 500 points in `~17 ms` (`< 1 s`), agrees with point-by-point direct HDM to `< 2e-13` relative $S_{21}$ error, and matches `cavity` BBQ eigenfrequencies to `< 0.002%`.
* `physics_ems.v1.json` and `physics_ems_result.v1.json` gain `driven`, `lumped_port`, and `waveport` descriptors and result shapes (additive).

## Code

`crates/nyx-physics/src/cqed.rs` (`DrivenPortData`), `crates/nyx-ems/src/{driven,ports,lib,problem}.rs`; test `crates/nyx-ems/tests/driven_prom.rs`.
