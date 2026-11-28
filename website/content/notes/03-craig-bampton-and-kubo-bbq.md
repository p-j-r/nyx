+++
title = "Craig–Bampton Hierarchical Quantisation & Two-Pole Kubo Black-Box Equivalence"
date = 2026-11-24
description = "Port-interior Craig–Bampton decomposition, companion QEP and Schur-reduced Black-Box Quantisation, 4-way SQuADDS equivalence to 1e-6%, and the rank-1 port-exclusion trade-off."

[extra]
note_number = "03"
crates = "nyx-analysis::hamiltonian · nyx-analysis::network · nyx-ems::rom"
reading_time = "10 min"
+++

Normal-mode Fock quantisation (`basis = "normal_fock"`) works well for weakly anharmonic transmons when truncated at $6\text{--}8$ levels per mode, because the cosine nonlinearity $\cos\hat{\varphi}$ only mixes nearby Fock states. On a Fluxonium, Cooper-pair box, or flux-biased DC-SQUID dissipator, however, multi-well wavefunctions and steep inductive barriers require dozens of local port states per junction; expanding the full tensor product of cavity and junction modes explodes exponentially.

In `crates/nyx-analysis/src/hamiltonian.rs`, `nyx` resolves this via **Mode B (`basis = "craig_bampton"`)** and cross-checks it against **Schur-complement Black-Box Quantisation (`bbq_modes`)** and **two-pole Kubo linear-response synthesis** in `crates/nyx-analysis/src/network.rs`.

## Nonlinear-Port + Linear-Cavity Craig–Bampton Basis

Starting from the 3D electromagnetic `CircuitRom` extracted by `nyx-ems::cavity`, let $\mathbf{B}\_{r,\mathrm{nl}}$ select the $J$ nonlinear Josephson ports from the $r$-dimensional reduced-order space $\mathbb{R}^r$. `build_linear_cb_decomposition` decomposes $\mathbb{R}^r$ into:

1. **Nonlinear port lift $T\_p$** satisfying $\mathbf{B}\_{r,\mathrm{nl}}^T T\_p = I\_J$, which isolates the local port flux and charge operators across each junction.
2. **Short-circuit interior cavity modes $T\_c$** satisfying $\mathbf{B}\_{r,\mathrm{nl}}^T T\_c = 0$, obtained by diagonalising the linear cavity blocks with the junction ports shorted.

Each 1D nonlinear port $j$ is pre-diagonalised independently before coupling to the cavity modes:

- When $E\_{L,\mathrm{eff}} > 1\,\text{kHz}$ (Fluxonium / shunted loops), the 1D port Hamiltonian is diagonalised in a `40`-state harmonic oscillator basis (`n_flux_basis = 40`), with $\exp(i\hat{\varphi}\_j)$ evaluated once via Padé-13 `expm`.
- When $E\_{L,\mathrm{eff}} \le 1\,\text{kHz}$ (charge-sensitive transmon / Cooper-pair box / DC-SQUID), it uses a `31`-state charge-periodic basis (`n_charge_cut = 15`).

Projecting $(n\_j, \varphi\_j)$ onto the lowest `port_trunc` (default `8`) eigenstates of each port reduces the full coupled Hilbert space (capped at `4096`) so tightly that a full DC-SQUID flux sweep executes in under `100 us` per flux point.{% <sidenote id="sn-labeling" label="1"> %}Dressed eigenvalues are mapped back to bare product states `{|0...0>, |1_a>, |2_a>, |1_a 1_b>}` using globally bijective greedy descending bipartite matching on `|<n | psi_k>|^2` to extract signed `f01_Hz`, `alpha_Hz`, and `chi_matrix_Hz`.{% </sidenote> %}

## Companion QEP, Schur BBQ & Taubin `S_21` Circle Fits

In parallel, `crates/nyx-analysis/src/network.rs` provides linearised frequency-domain and complex-eigenmode reductions on the same `CircuitRom` or `ems::driven` port data:

- **Companion QEP:** Linearises $(K\_{\mathrm{tot}} + i k\_0 C\_{\mathrm{tot}} - k\_0^2 M\_{\mathrm{tot}}) \mathbf{z} = \mathbf{0}$ into a balanced $2r \times 2r$ complex matrix via `femcore::linalg::dense::schur`, returning complex eigenfrequencies $f\_{\mathrm{re}} + i f\_{\mathrm{im}}$, external quality factors $Q\_c$, and linewidths $\kappa = 2 |f\_{\mathrm{im}}|$.
- **Schur-Complement BBQ (`bbq_modes`):** Partitions ports into active junction ports $A$ and $50\,\Omega$-terminated feedline ports $T$, eliminates $T$ exactly via

$$
Y\_{\mathrm{red}}(\omega) = Y\_{AA}(\omega) - Y\_{AT}(\omega)\left(Y\_{TT}(\omega) + Y\_{\mathrm{load}}\right)^{-1} Y\_{TA}(\omega)
$$

locates Foster parallel-resonance zeros $\mathrm{Im}[Y\_{JJ}(\omega\_m)] = 0$, and extracts $C\_{\Sigma,m} = \frac{1}{2}\,\mathrm{Im}[Y\_{JJ}'(\omega\_m)]$, $L\_{\mathrm{eff},m}$, and Purcell lifetimes $T\_{1,\mathrm{Purcell}}$.
- **Taubin `S_21` Circle Fit (`s21_fit`):** Fits driven scattering curves on the complex plane via a $3 \times 3$ centered Taubin algebraic circle fit plus Gauss–Newton arctan phase extraction. On `C0.3` ($\lambda/4$ CPW hanger) and `12`-slab `S2.1` (`r0105`, `r1192`), `CircuitRom` QEP and 2-port `ems::driven` Taubin agree to **`Delta f_r = 0.0000%–1.67e-6%`** and **`Delta kappa = 0.0042%–0.257%`**.

## 4-Way `QubitCavity` Equivalence & Cryogenic Hardware

On `12` coupled transmon+resonator chips (`6` $\lambda/4$ `CLT` + `6` $\lambda/2$ `NCap`, `149.8k–221.8k` tets on identical 3D meshes), `docs/SQUADDS_BENCHMARK_REPORT.md` compared closed `craig_bampton` (`M1`), passive-port QEP + `craig_bampton` (`M2`), and 3-port `ems::driven` Schur BBQ + exact two-pole Kubo susceptibility (`M3`):

- **`M1` vs. `M2` (`12/12` PASS):** `med |Delta f_01| = 1.38e-6%` (max `0.0120%`), `med |Delta alpha| = 2.65e-6%` (max `1.918%`), and `CLT` `|Delta f_r| <= 9.89e-4%`.
- **`M1` vs. `M3` Kubo BBQ:** `med |Delta f_01| = 5.0e-6%` (`<= 0.0001%` on `11/12` designs), `med |Delta alpha| = 0.0056%` (`<= 0.0941%` on `12/12`), and **Thomas–Reiche–Kuhn (TRK) oscillator-strength sum rules** (`16 E_C |<k|n|m>|^2 / f_mk`) satisfied to **`1.000000000000 +/- 2.7e-10`** across all `12` designs.

{% <solver_log cmd="cargo test --release -p nyx-analysis --test squadds_smoke -- --test-threads=1"> %}
transmon_cross_3net_mode_g2_far_field_and_i_edge_thickness_smoke ... ok
dissipator_dc_squid_craig_bampton_vs_rom_port_flux_sweep_smoke ... ok
  4-mode S4.3 Dissipator ROM (qubit, dissipator_jj, readout, filter):
  25-point native squid_flux_sweep (10.32 -> 2.93 GHz in < 0.25 s)
  matches point-by-point ej_hz_override craig_bampton eigenvalues to < 1e-12
muninn_3d_tls_participation_nine_param_oxide_and_purcell_smoke ... ok
{% </solver_log> %}

On full cryogenic hardware suites (`Phase 1`, `30/30 ok`):

- **`S4.1` `WM1` (6-qubit MIT-LL chip):** With mask-verified `+0.9 um` etch bias (`W1` PEC), `ems::cavity` hits **`1.32%` RMSRE** on `f_r` (`-1.52%` to `+1.43%` across `Q1..Q6`), **`2.3x` closer to mK hardware** than the published `3.0%` SQuADDS Ansys HFSS reference.
- **`S4.2` `MUNINN`:** Reproduces exact experimental coherence ranking (**`Ancilla 12.82 ms > Qubit 11.22 ms > Env 8.19 ms`** matching `32 > 31 > 28 us`) and solves the full `1,085,076`-tet (`1.402M` Nédélec DOF) 3-qubit chip in `82.55 s` (`17.51 GiB`), extracting the 21-pair `ZZ` coupling matrix.
- **`S4.3` `Dissipator`:** Solves a `464,520`-tet 4-mode `craig_bampton` model in `46.56 s` and sweeps `100` DC-SQUID flux points (`675`-dim Hilbert space) in **`1.10 s`** (`11.0 ms/pt`), resolving the `2g = 166.8 MHz` avoided crossing cleanly.

## Honest Disclosure: The `rank1` Port-Exclusion Trade-Off

In `crates/nyx-ems/src/cavity.rs` (`L1280–1310`, `Core Bug 2`), refining a `rank1` `lumped_port` or `lumped_jj` sheet inserts interior/lateral-vacuum port vertices where `(G^T b_j)_v == 0` into `cond_vertices`, dropping their zero-frequency gradient columns from `g_red` and stalling `krylov_schur`. Enforcing an unrefined `30 um` port-exclusion box prevents every solver hang (`163/163` EM `ok`), but forces P1 elements across narrow `NCap` port gaps to overestimate port capacitance (`S2.2` `kappa` HT median `435.75%` in raw QEP vs. `41.59%` when `C_c` is extracted via `dfs::electro` in `hybrid` mode; `WM1` `alpha` RMSRE `27.32%`). Similarly, on `1/12` `M3` runs (`s2_3_p1_clt_03`), capping `complex_solve` to `8` vertical slabs (`<= 222k` tets to stay under `19.71 GiB`) under-resolved sub-kHz `CLT` hanger coupling.
