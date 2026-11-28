# ADR 137: Gyrotropic and non-symmetric 3x3 electromagnetic tensors

Status: Accepted

## Context

Magneto-optic and gyromagnetic photonic and microwave devices—including Faraday rotators, non-reciprocal Ce:YIG-on-silicon waveguide phase shifters, and ferrite circulators/isolators—break Lorentz reciprocity through an applied static magnetic bias $\mathbf{B}_0$ or remnant magnetization $\mathbf{M}_0$. Under the engine's harmonic time convention $e^{+j\omega t}$, their relative permittivity and permeability tensors take the non-symmetric (Hermitian when lossless) form

$$\boldsymbol{\varepsilon}_r(\omega) = \boldsymbol{\varepsilon}_{\mathrm{sym}}(\omega) + i\,\mathbf{G}_\times(\omega), \qquad \boldsymbol{\mu}_r(\omega) = \boldsymbol{\mu}_{\mathrm{sym}}(\omega) + i\,\mathbf{M}_\times(\omega),$$

where $\mathbf{g} = (g_x, g_y, g_z)$ and $\mathbf{g}_\mu = (g_{\mu,x}, g_{\mu,y}, g_{\mu,z})$ are axial gyrotropy vectors (Voigt vectors) acting via the cross product $\mathbf{G}_\times \mathbf{E} = \mathbf{g} \times \mathbf{E}$:

$$\mathbf{G}_\times = \begin{bmatrix} 0 & -g_z & +g_y \\ +g_z & 0 & -g_x \\ -g_y & +g_x & 0 \end{bmatrix}, \qquad \varepsilon_{xy} = -i g_z,\quad \varepsilon_{yx} = +i g_z,\quad \varepsilon_{xz} = +i g_y,\quad \varepsilon_{zx} = -i g_y,\quad \varepsilon_{yz} = -i g_x,\quad \varepsilon_{zy} = +i g_x.$$

Because $\mathbf{g}$ is an axial vector (pseudovector), under an orthogonal crystal-to-wafer Bunge Euler transformation $\mathbf{R}$ it transforms as $\mathbf{g}_{\mathrm{global}} = (\det \mathbf{R})\,\mathbf{R}\,\mathbf{g}_{\mathrm{local}}$, preserving the skew-symmetric tensor equivalence $\mathbf{R}\,\mathbf{G}_{\times,\mathrm{local}}\,\mathbf{R}^T = (\det \mathbf{R})[\mathbf{R}\mathbf{g}_{\mathrm{local}}]_\times$.

Previously, `EpsC` stored only the six upper-triangle entries of a complex-symmetric permittivity tensor, `MuC` stored only three diagonal permeability entries, `modes_complex` assumed $\varepsilon_{ab} = \varepsilon_{ba}$ in the longitudinal-transverse waveguide coupling blocks, `ShiftedFactorization` dispatched complex matrices with matching sparsity pattern and diagonal to complex-symmetric $L B L^T$ (`SymLblt`), and `modes_complex` / `SparseDrivenSolver` used the complex-symmetric left-eigenvector shortcut $\mathbf{y} = \bar{\mathbf{x}}$. When $\boldsymbol{\varepsilon}$ or $\boldsymbol{\mu}$ has a skew-symmetric gyrotropic part ($A^T \neq A$), complex-symmetric $L B L^T$ factorization and $\mathbf{y} = \bar{\mathbf{x}}$ are invalid.

At the same time, reciprocal dielectric, superconducting, and PML simulations rely on the speed and memory efficiency of `SymLblt` and the bit-exactness of existing golden oracles.

## Decision

1. **Material representation and axial-vector rotation (`femcore::material`, `femcore::props`, `femcore::optical`)**:
   - Extend `EpsC` and `MuC` to full 9-entry complex $3 \times 3$ tensors (`xx, yy, zz, xy, yx, xz, zx, yz, zy: C`) with `is_non_symmetric() -> bool` (`|T_{ab} - T_{ba}| > 1e-15`), `from_matrix`, `to_matrix`, and exact complex $3 \times 3$ matrix inversion `MuC::inv_matrix()`.
   - When permeability is scalar isotropic (`mu_r`), retain the exact scalar division path `c(1.0 / mur, 0.0)` so all reciprocal golden tests remain 0-ULP bit-identical.
   - Add `rotate_axial_vector(r, g)` ($\mathbf{g}_{\mathrm{global}} = (\det \mathbf{R})\,\mathbf{R}\,\mathbf{g}_{\mathrm{local}}$), `skew_from_axial`, and `axial_from_skew` in `femcore::material::tensor`.
   - Support `[optical.gyrotropy]` in TOML material files with `model = "constant"`, `"lorentz"` (multi-oscillator magneto-optic Lorentz poles along unit axis $\hat{\mathbf{u}}$), and `"tabulated"` ($\lambda \mapsto (g_x, g_y, g_z)$), as well as `[dielectric]` / `[electrical]` and inline problem material properties `gyrotropy` / `gyrotropy_im`, `mu_gyrotropy` / `mu_gyrotropy_im`, and full 9-entry complex tensor matrices in `ems`.

2. **Three-stage non-symmetric assembly and solver dispatch (`nyx-ems`, `femcore::linalg`)**:
   - Track `is_non_symmetric: bool` in Stage 1 across 2D waveguide mode assembly (`ModeAssemblyC`), 3D cavity/eigenmode assembly (`EmwAssembly`), 3D Floquet-Bloch assembly (`BlochAssembly`), 3D driven frequency-domain assembly (`DrivenAssembly`), and native mixed-cell materials (`MixedCellMaterials`).
   - In 2D waveguide `modes_complex`, route problems with `is_non_symmetric()` or off-diagonal $\boldsymbol{\mu}$ to the full 9-entry element kernel `add_waveguide_complex_point_full`, which assembles distinct $(\mathbf{e}_t, \phi)$ (`-mg` using $\varepsilon_{02}, \varepsilon_{12}$) and $(\phi, \mathbf{e}_t)$ (`+mgt` using $\varepsilon_{20}, \varepsilon_{21}$) longitudinal-transverse coupling blocks as well as transverse-longitudinal magnetic $\boldsymbol{\mu}^{-1}$ curl coupling (`a1_ee`), while keeping the reciprocal fast path untouched when `!is_non_symmetric && !has_off_diag_mu`. Support negative target effective index (`n_target < 0`) to select backward-propagating ($-\beta$) branches of the quadratic eigenvalue problem directly.
   - When `is_non_symmetric == false`, preserve `SymLblt` complex-symmetric factorization and the $\mathbf{y} = \bar{\mathbf{x}}$ shortcut without alteration.
   - When `is_non_symmetric == true`:
     - Automatically dispatch Stage 2 shifted factorizations (`DirectSolver::new_general_lu`, `SparseDrivenSolver::new_with_symmetry`, `ShiftedBlrLu`, `AmsFgmresSolver`) to general complex LU (`LuColamd`) instead of `SymLblt`.
     - Reject `solver = "lblt"` in `ems::driven` with a clear diagnostic error directing the caller to `"auto"`, `"direct"`, `"blr"`, or `"ams_fgmres"`.
     - Compute true left eigenvectors $\mathbf{A}^H \mathbf{y} = \bar{\lambda}\,\mathbf{B}^H \mathbf{y}$ via `femcore::linalg::sens::left_eigenvector` in `modes_complex`, normalize them biorthogonally ($\mathbf{y}_i^H \mathbf{B} \mathbf{x}_i = 1$), and expose `ModesC::{et_left, phi_left, biorth_norms}` and `biorthogonal_poynting_overlap`.
     - Report `PromOperator::is_symmetric(&self) -> bool { !self.is_non_symmetric }` on `DrivenAssembly` so multi-port driven sweeps and PROM projections solve every active port column without symmetrizing $S_{ij} \neq S_{ji}$.

## Consequences

- **Positive**: Faraday rotation, non-reciprocal waveguide phase shift ($\Delta\beta_{\mathrm{NR}} = \beta_+ - \beta_-$), and multi-port magneto-optic / ferrite isolators and circulators ($|S_{21}| \neq |S_{12}|$, $\mathbf{S}^H \mathbf{S} = \mathbf{I}$) are solved natively on both simplex and mixed-cell meshes with $P_1/P_2$ Nédélec elements and curvilinear PMLs.
- **Positive**: Zero runtime or memory overhead and 0-ULP bit-identity on all reciprocal problems (`is_non_symmetric == false`).
- **Negative**: Non-symmetric factorizations require general complex LU (`LuColamd`) and an explicit adjoint triangular solve for true left eigenvectors, approximately doubling factorization memory compared to `SymLblt` when gyrotropy is active.

## Code

- `crates/femcore/src/material/{mod.rs,tensor.rs,optics.rs,coupled.rs,spec.rs}`
- `crates/femcore/src/props.rs`
- `crates/femcore/src/optical.rs`
- `crates/femcore/src/linalg/{direct.rs,driver.rs,prom.rs}`
- `crates/nyx-physics/src/{validate.rs,decode.rs}`
- `crates/nyx-ems/src/{modes.rs,modes_complex.rs,cavity.rs,bloch.rs,driven.rs,mixed_ops.rs}`
- `crates/nyx-ems/tests/gyrotropic.rs`
