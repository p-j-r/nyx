# 110. Periodic 3D lattice conformance, per-layer region `sweep`, and fixed-frequency 3D Bloch modes with dispersive Poynting group index

Date: 2026-11-21

## Status

Accepted. Amends ADR 40, 45, 51, 66, 67, and 80.

## Context

Subwavelength grating (SWG) waveguides, photonic-crystal nanobeams, and corrugated Bragg couplers are periodic along the wafer propagation axis $+y$ with pitch $\Lambda$, while engineered dispersion requires group velocity $v_g$ and group index $n_g = c_0 / v_g$ from a single fixed-frequency solve rather than finite-differencing two eigenfrequency solves. Previously:

1. `nyx-ems` solved fixed-wavevector 3D `cavity` eigenfrequencies $\omega(\mathbf{k})$, requiring root-finding sweeps to locate $k_y(\omega)$ at fixed laser wavelength $\lambda_0$ and failing to capture radiative complex-$k_y$ leaky grating modes above the cladding light line.
2. Second-order Nédélec (`Ned3`) periodic congruence (`cavity::ned_relations`) requires periodic source and destination faces to match vertex-for-vertex, edge-for-edge, and triangle-for-triangle (`face_id[&sorted3(vt)]`), which independent BCC lattice warping or Rivara refinement could break when local `add_size` rules touched a periodic boundary.
3. `SweepLayer` supported only uniform `domain_offset` or whole-layer `fill`, preventing multi-region per-layer remapping (such as turning grating teeth into cladding in the etched half-period) and 3D `geom` model extrusion.

## Decision

* **Complex quadratic shift-invert eigensolver (`femcore::linalg::driver::eigs_quadratic_complex`):**
  * Linearizes $(\mathbf{A}_0 + \nu \mathbf{A}_1 + \nu^2 \mathbf{A}_2)\mathbf{u} = \mathbf{0}$ around complex shift $\sigma \in \mathbb{C}$ via `QepShiftInvert` of dimension $2n$, factoring $\mathbf{P}(\sigma) = \mathbf{A}_0 + \sigma \mathbf{A}_1 + \sigma^2 \mathbf{A}_2$ once with `ShiftedFactorization`, honoring cancellation and warm starts, and reporting true quadratic residuals.
* **3-stage 3D Floquet–Bloch study (`nyx-ems::bloch`, `Study::Bloch`):**
  * Propagates along wafer $+y$ as $\mathbf{E}(x,y,z) = \mathbf{e}(x,y,z) e^{-j k_y y}$ with periodic cell amplitude $\mathbf{e}(x, y + \Lambda, z) = \mathbf{e}(x, y, z)$ (`phase = 0` on the longitudinal $y$-periodic pair, optional Floquet phase on transverse periodic pairs) and unconjugated transverse $(x, z)$ PML stretches paired with Hermitian congruence $\mathbf{P}^H (\cdot) \mathbf{P}$.
  * Solves the quadratic pencil in normalized Bloch index $\nu = n_B = k_y / k_0$ shifted around $\sigma = n_{\text{target}} + (m - 1)\lambda_0 / \Lambda$ (`brillouin_zone = m`), reconstructs physical $\mathbf{E}$ and $\mathbf{H}$ on all 10 P2 tetrahedral nodes, and populates `FieldFrame::outside` with 3D PML domains.
  * Enables `pmc` across all `ems` studies and adds midplane `symmetry` (`none`, `pec_x`, `pmc_x`, `pec_z`, `pmc_z`, `pec_x_pmc_z`, `pmc_x_pec_z`, `pmc_x_pmc_z`, `pec_x_pec_z`) on `modes`, `cavity`, and `bloch`.
* **Parent-to-child 3D periodic meshing & per-layer `SweepLayer::regions` (`femcore::meshing`, `nyx-engine`):**
  * `GeomDomain3` extracts axis-aligned `PeriodicAxis3 { axis, lo, hi }` pairs; `mesher3d` quantizes periodic axes to integer cell counts ($\ge 2$), restricts parent-plane corner warping to exact `lo`, mirrors `warped[child] = warped[parent] + shift`, synchronizes periodic equivalence classes across flat-tet undo rounds and sliver welds, inherits parent cut bisections, and uses canonical tie-break IDs `(canon_id, id)` in `quad_tris` and `split_prism` for identical quad diagonals; `refine3d::refine_to_size_periodic` and `curve3` synchronize edge bisections and curved mid-edge projections parent-to-child; `adapt3d` freezes 3D periodic face patches.
  * `SweepLayer` adds optional `regions: Option<BTreeMap<String, String>>` for per-layer region remapping in 2D (`extrude_y_with_remap`) and 3D (`sweep_model_3d`).
* **Single-solve dispersive Poynting $v_g$, $n_g$, and mode scalars (`nyx-ems`):**
  * Evaluates time-averaged longitudinal power flow $P_y = \frac{1}{2}\operatorname{Re}\int (\mathbf{E} \times \mathbf{H}^*)_y\,d\Omega$ and Brillouin dispersive electromagnetic energy $U_{\text{EM}} = \frac{1}{4}\int \left(\varepsilon_0 \mathbf{E}^* \frac{\partial(\omega\boldsymbol{\varepsilon}_r)}{\partial\omega}\mathbf{E} + \mu_0 \mathbf{H}^* \boldsymbol{\mu}_r \mathbf{H}\right) d\Omega$ over non-PML cells using $\frac{\partial(\omega\varepsilon_r)}{\partial\omega} = \operatorname{Re}(\tilde{n}^2 - 2\lambda \tilde{n}\,d\tilde{n}/d\lambda)$ from `Optical::nk`, emitting `v_g N (m/s)`, `n_g N`, `Ex/Ey/Ez share N`, `share in domain <label> N`, and `layer fraction N` across 2D `modes`/`modes_complex` and 3D `bloch`/`cavity`.

## Consequences

* Fixed-wavelength 3D periodic waveguide and leaky grating dispersion, radiative loss $\alpha$ (dB/cm), and dispersive group index $n_g$ are obtained in one eigensolve without finite-difference frequency pairs.
* Arbitrary local `add_size` rules on cores and grating teeth preserve exact `Ned3` face-and-edge periodic congruence.
* Legacy `nyx-ems` golden files remain untouched (`nyx_physics::testkit` ignores additive FR-3 scalars when absent from legacy golden expectations).

## Code

`crates/femcore/src/linalg/driver.rs`, `crates/femcore/src/optical.rs`, `crates/femcore/src/meshing/{mod,geomdrive3,mesher3d,refine3d,curve3,mesher}.rs`, `crates/femcore/tests/{linalg_contract,periodic_mesh3d}.rs`, `crates/nyx-ems/src/{bloch,modes,modes_complex,cavity,rom,axisym,optical,problem,lib}.rs`, `crates/nyx-ems/tests/{bloch3d,group_index_and_scalars,contracts}.rs`, `crates/nyx-physics/src/testkit.rs`, `crates/nyx-analysis/src/samples.rs`, `crates/nyx-engine/src/{project,preview,meshing,commands}.rs`, `crates/nyx-engine/tests/{swept_cell,periodic_photonics_e2e}.rs`.
