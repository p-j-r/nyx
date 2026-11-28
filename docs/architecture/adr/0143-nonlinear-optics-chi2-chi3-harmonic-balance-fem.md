# ADR 144: Nonlinear optics ($\chi^{(2)}$ and $\chi^{(3)}$ Kerr, TPA/FCA, Harmonic Balance FEM, and nonlinear waveguide modes)

## Status

Accepted (2026-11-27). Implements FR-41.

## Context

Integrated photonic and quantum-optical devices rely on bulk second- and third-order optical nonlinearities: quasi-phase-matched second-harmonic, sum-frequency, difference-frequency, and optical parametric generation ($\chi^{(2)}$ SHG/SFG/DFG/OPA) in non-centrosymmetric crystals such as $\text{LiNbO}_3$, and Kerr self-/cross-phase modulation ($\chi^{(3)}$ SPM/XPM), optical bistability in high-$Q$ microresonators ($\text{Si}_3\text{N}_4$), and two-photon plus free-carrier absorption/dispersion (TPA + FCA/FCD) in semiconductor waveguides ($\text{Si}$).

To model these phenomena accurately without external scripting glue or scalar polarization approximations:

1. Material files (`format = 2`, ADR 109) need a unit-checked, provenance-tracked `[nonlinear]` section that enforces crystallographic `PointGroup` selection rules (forbidding bulk $\chi^{(2)}$ in centrosymmetric crystals unless explicitly overridden), rotates rank-3 $d_{iJ}$ and cubic rank-4 $\chi^{(3)}_{ijkl}$ tensors via Bunge Z-X-Z Euler angles, represents spatial quasi-phase-matching (`QpmSpec`) poling domains, and exposes shared pointwise time-domain polarization evaluators for transient solvers.
2. Frequency-domain 3D simulations must capture full cross-polarization vector Kerr coupling ($\mathbf{E}\cdot\mathbf{E}^*$ and conjugate $\mathbf{E}\cdot\mathbf{E}$ Maker–Terhune terms), steady-state TPA-generated free-carrier density $N_c(\mathbf{E})$, multi-harmonic chromatic dispersion $n(\omega_m)$, depleted pump back-action satisfying Manley–Rowe photon-flux conservation, and saddle-node turning points along S-shaped Kerr cavity bistability curves where power-controlled continuation fails ($\det\mathbf{J} = 0$).
3. Waveguide cross-section design requires the full-vectorial nonlinear parameter $\gamma$ ($\text{W}^{-1}\text{m}^{-1}$), nonlinear effective area $A_{\text{eff}}$, TPA figure of merit, and fixed-power continuation of $n_{\text{eff}}(P)$.

## Decision

### 1. Material `format = 2` `[nonlinear]` schema, symmetry validation, and time-domain evaluators (`femcore::material::nonlinear`)

Add `Nonlinear` (`chi2: Option<Chi2>`, `chi3: Option<Chi3>`, `qpm: Option<QpmSpec>`, `t_range`, `source`) to `MaterialSpec` and `Anisotropy`:

- **`[nonlinear.chi2]`**: Contracted $3 \times 6$ optical $d_{iJ}$ tensor (`(xx, yy, zz, yz, xz, xy)`, SI $\text{m/V}$, accepting `"pm/V"` units or named coefficients `d33`, `d31`, `d15`, `d22`, `d14`, `d36` expanded by crystal class and optional `kleinman = true`). `validate_chi2_point_group` rejects non-zero $d_{iJ}$ on all 11 centrosymmetric point groups and chiral cubic `432` by Neumann's principle (unless `broken_centrosymmetry = true`) and checks non-centrosymmetric symmetry relations (`3m`, `6mm`, `4mm`, `-43m`, `23`, `-42m`, `-6m2`, `mm2`). `Chi2::d_global(euler_deg)` expands $d_{iJ} \to d_{ijk}$, applies the Bunge Z-X-Z orthogonal rotation, and contracts back to $3 \times 6$.
- **`[nonlinear.qpm]`**: Spatial poling sign (`sign = +1 | -1`) or periodic poling grating (`period` / `Lambda_qpm`, `axis`, `duty_cycle`, `phase_offset`) evaluated pointwise via `QpmSpec::poling_sign_at([x, y, z])`.
- **`[nonlinear.chi3]`**: Supports isotropic Kerr index `n_2` ($\text{m}^2/\text{W}$) or explicit tensor entries `chi3_xxxx`, `chi3_xxyy`, `chi3_xyxy`, `chi3_xyyx` ($\text{m}^2/\text{V}^2$) and Maker–Terhune ratio `b_over_a` (default `1.0`), yielding the full vector Kerr polarization
  $$\mathbf{P}^{(3,\text{Kerr})}(\omega) = \varepsilon_0\Big[A(\mathbf{E}\cdot\mathbf{E}^*)\mathbf{E} + \tfrac{1}{2}B(\mathbf{E}\cdot\mathbf{E})\mathbf{E}^*\Big]$$
  plus rotated cubic (`m-3m` / `-43m`) anisotropy $\delta\chi^{(3)} = \chi_{xxxx} - (\chi_{xxyy} + \chi_{xyxy} + \chi_{xyyx})$, two-photon absorption `beta_tpa` ($\text{m/W}$), free-carrier absorption cross-section `sigma_fca` ($\text{m}^2$), signed free-carrier dispersion `k_fcd` ($\text{m}^3$), and recombination lifetime `tau_c` ($\text{s}$).
- **Shared time-domain evaluators**: Export `eval_chi2_polarization_time`, `eval_chi2_tangent_time`, `eval_chi3_kerr_time`, `eval_chi3_envelope_time`, and unconditionally positive semi-analytic `drude_carrier_step` so transient `dgtd` and `fdtd` kernels consume the identical `[nonlinear]` TOML definitions.
- Ship literature-sourced `[nonlinear]` sections in `materials/lithium_niobate.toml` (`3m`, Shoji 1997), `materials/silicon_nitride.toml` (Ikeda 2008), and `materials/silicon.toml` (`m-3m`, Bristow 2007 & Soref–Bennett 1987).

### 2. Three-stage `Study::NonlinearDriven` (`"nonlinear_driven"`) in `nyx-ems::nonlinear`

Implement `assemble_nonlinear_driven`, `solve_nonlinear_driven_system`, and `nonlinear_driven_fields` on 3D Nédélec (`P1`/`P2`, straight and curved) meshes with `waveport` and `lumped_port` excitations (mapping RMS port power $P_{\text{in}}$ to peak phasor drive $\sqrt{2 P_{\text{in}}}$):

- **Single-tone $\chi^{(3)}$ Kerr + TPA + FCA/FCD Harmonic Balance FEM (`process = "chi3"`)**:
  - Assembles the real $2 N_{\text{red}} \times 2 N_{\text{red}}$ split $(\operatorname{Re}\hat{\mathbf{y}}, \operatorname{Im}\hat{\mathbf{y}})$ residual $\hat{\mathbf{F}}(\hat{\mathbf{y}}, \alpha) = \hat{\mathbf{A}}(\omega)\hat{\mathbf{y}} - k_0^2 \hat{\mathbf{g}}_{\text{NL}}(\hat{\mathbf{y}}) - \alpha\hat{\mathbf{b}}_{\text{ref}} = \mathbf{0}$ with steady-state carrier density $N_c(\mathbf{E}) = \frac{\tau_c \beta_{\text{TPA}}}{2\hbar\omega} I(\mathbf{E})^2$ and exact Wirtinger Fréchet Jacobian $\mathbf{J} = (\partial\hat{\mathbf{F}}/\partial\hat{\mathbf{y}}, \partial\hat{\mathbf{F}}/\partial\hat{\mathbf{y}}^*)$ (honoring Floquet phase ratios and Ruiz equilibration factors).
  - Reuses a static CSR sparsity pattern (`StaticRealPattern`) across `femcore::linalg::snes::newton` (and LU-preconditioned `snes::jfnk`).
  - Supports `continuation = "power_ramp"` and bordered $(2 N_{\text{red}} + 1) \times (2 N_{\text{red}} + 1)$ Keller **pseudo-arclength continuation** (`continuation = "arclength"`), stepping through fold singularities ($\det\mathbf{J} = 0$) with a single symbolic factorization and extracting saddle-node turning-point powers via cubic Hermite interpolation of $\dot{\alpha}(s) = 0$.
- **Multi-harmonic $\chi^{(2)}$ (`shg`, `sfg`, `dfg`, `opa`) & two-tone `chi3_xpm` Coupled-Wave FEM**:
  - Evaluates library `[optical]` dispersion independently at each harmonic wavelength $\lambda_m = c_0 / f_m$ (`clone_problem_at_freq`) and samples `QpmSpec::poling_sign_at` at every 3D Gauss quadrature point.
  - Supports fast `undepleted` sequential solves and exact `depleted` coupled block-Newton (`snes::newton` / `snes::jfnk`) with full parametric back-action ($\mathbf{P}_1 = 2\varepsilon_0 d : \mathbf{E}_2\mathbf{E}_1^*$ for SHG; $\mathbf{P}_1 = 2\varepsilon_0 d : \mathbf{E}_3\mathbf{E}_2^*, \mathbf{P}_2 = 2\varepsilon_0 d : \mathbf{E}_3\mathbf{E}_1^*$ for 3-wave mixing), reporting `shg_efficiency` and discrete Manley–Rowe photon-flux conservation errors.

### 3. Three-stage 2D `Study::NonlinearModes` (`"nonlinear_modes"`) in `nyx-ems::nonlinear`

Implement `assemble_nonlinear_modes`, `solve_nonlinear_modes_continuation`, and `nonlinear_modes_fields` on 2D waveguide cross-sections:

- Normalizes the linear eigenmode $(\mathbf{e}_t, \tilde{\boldsymbol{\phi}}, \beta_0)$ to unit cross-section Poynting power and evaluates the full-vectorial nonlinear waveguide parameter $\gamma = \frac{\omega}{4 P_0^2}\operatorname{Re}\int \mathbf{E}_0^* \cdot \mathbf{P}^{(3,\text{Kerr})}(\mathbf{E}_0)\,\mathrm{d}A$, TPA parameter $\gamma_{\text{TPA}}$, effective mode area $A_{\text{eff}}$, and `fom_tpa`.
- Tracks power-dependent complex effective index $n_{\text{eff}}(P)$ via bordered $(2 N_{\text{free}} + 2) \times (2 N_{\text{free}} + 2)$ fixed-power + phase-gauge Newton continuation with exact $\partial/\partial\beta$ and conjugate $\partial/\partial\beta^*$ derivatives on $E_y = i(\beta / k_0)\tilde{\phi}$.

## Consequences

- **Positive**: Depleted SHG/SFG/DFG/OPA, QPM periodic poling, vector Kerr SPM/XPM, silicon TPA+FCA saturation, high-$Q$ Kerr microresonator bistability S-curves, and 2D waveguide $\gamma$ / $n_{\text{eff}}(P)$ run natively in pure Rust with quadratic Newton convergence, single-symbolic-analysis reuse, and bit-exact `NYX_THREADS` determinism.
- **Trade-offs**: Multi-harmonic coupled Newton scales with the number of active harmonic blocks ($2\times$ DOFs for SHG/XPM, $3\times$ DOFs for 3-wave mixing); weak-conversion sweeps can select `coupling = "undepleted"` or `nl_solver = "jfnk"` to avoid coupled block factorizations.

## Code

- `crates/femcore/src/material/{nonlinear.rs,mod.rs,spec.rs,schema.rs,reader.rs,units.rs}`
- `materials/{lithium_niobate.toml,silicon_nitride.toml,silicon.toml}`
- `crates/nyx-physics/src/decode.rs`
- `crates/nyx-engine/src/materials.rs`
- `crates/nyx-ems/src/{nonlinear.rs,dgtd.rs,fdtd,problem.rs,lib.rs,ports.rs,sensitivity.rs}`
- `crates/femcore/tests/{nonlinear_material.rs,material_contract.rs,material_library.rs}`
- `crates/nyx-ems/tests/{fr41_nonlinear.rs,contracts.rs,finite_guard.rs}`
- `crates/nyx-engine/tests/phase3_integration.rs`

## Contract

- `material_schema` `2.2` (`docs/architecture/contracts/material_schema.v2.json`)
- `material_example` `2.1` (`docs/architecture/contracts/material_kitchen_sink.v2.{toml,json}`)
- `physics_ems` `1.11` (`docs/architecture/contracts/physics_ems.v1.json`)
- `physics_ems_result` `1.5` (`docs/architecture/contracts/physics_ems_result.v1.json`)
- `kernel_stages` `1.11` (`crates/nyx-ems/tests/contracts.rs`)
