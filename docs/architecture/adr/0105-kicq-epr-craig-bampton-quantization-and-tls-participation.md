# 105. KICQ-EPR, hierarchical Craig–Bampton circuit quantization, and edge-regularized TLS participation

Date: 2026-11-21

## Status

Accepted (superconducting qubit milestones M4.2–M4.4). Extends the post-solve `Analysis` plugin architecture of ADR 56 inside `crates/nyx-analysis`.

## Context

Quantizing 3D superconducting circuits from full-wave EM results requires three dedicated post-solve analyses:
1. KICQ Energy Participation Ratio (`epr`) that tracks kinetic inductance separately from linear superinductance and allows instant sweeps over junction inductance without re-running 3D FEM.
2. Full quantum Hamiltonian diagonalization (`hamiltonian`) handling both weakly anharmonic transmon–cavity systems (multi-mode Fock cosine expansion) and strongly nonlinear/multi-well circuits like Fluxonium at $\Phi_{\text{ext}} = 0.5\Phi_0$ where linearized normal-mode $\omega_p^2$ goes negative, in pure Rust with microsecond/millisecond flux/charge sweeps.
3. Surface dielectric TLS loss participation (`tls_participation`) that agrees with explicit nanometre-scale oxide integrals in 2D (`MA`, `MS`, `SA`) and avoids logarithmic $\ln(w/h_{\text{edge}})$ divergence on 3D zero-thickness sheets under `convergence_check`.

## Decision

* **`EprAnalysis` (`id = "epr"`, `crates/nyx-analysis/src/epr.rs`).** Consumes `CircuitRom`, supports instantaneous `< 20 us` $r\times r$ generalized eigenvalue re-projection under `lj_override` or `ej_hz_override`, splits parallel linear superinductance $E_L$ from nonlinear Josephson energy $E_J$, reports signed modal participations $p_{mj}$, sign matrix $s_{mj}$, phase non-collinearity, $\varphi_{\text{zpf},mj}$, perturbative $\alpha_m$ and $\chi_{mn}$, kinetic inductance fractions $\alpha_{k,m}$, and `equipartition_error` (`< 1e-12`).
* **`HamiltonianAnalysis` (`id = "hamiltonian"`, `crates/nyx-analysis/src/hamiltonian.rs`).** Implements two native Rust quantization modes on `CircuitRom` using `femcore::linalg::{expm, dense::hermitian_eig}`:
  - `basis = "normal_fock"`: Multi-mode normal Fock product basis (`fock_trunc`), evaluating exact displacement exponentials $\exp(i\hat{\varphi}_{mj})$ minus quadratic parts;
  - `basis = "craig_bampton"`: Decomposes the passive ROM into port lift $\mathbf{T}_p$ ($\mathbf{B}_r^T \mathbf{T}_p = \mathbf{I}_J$) and short-circuit interior cavity modes $\mathbf{T}_c$ ($\mathbf{B}_r^T \mathbf{T}_c = 0$), solves each nonlinear 1D port independently (`E_L > 1e3 Hz` in a 40-state positive-definite bare $E_L\text{–}E_C$ harmonic-oscillator basis with cached `expm(i phi)` so $\omega_p$ never goes imaginary at half-flux; `E_L <= 1e3 Hz` in a 31-state charge-periodic $|n\rangle$ rotor basis with offset charge `n_g`), projects $(\hat{n}_j, \hat{\varphi}_j)$ onto the lowest `port_trunc` states, diagonalizes the coupled port–cavity tensor product, and runs full `phi_ext_sweep` / `ng_sweep` curves over `[0, 1]` in `< 10 ms`.
  - Both modes assign dressed eigenstate labels $|0\dots 0\rangle, |1_m\rangle, |2_m\rangle, |1_m 1_n\rangle$ via globally bijective greedy descending bipartite matching on $|\langle \mathbf{n} | \psi_k\rangle|^2$ to prevent state collisions near avoided crossings.
* **`TlsParticipationAnalysis` (`id = "tls_participation"`, `crates/nyx-analysis/src/tls.rs`).** Evaluates `p_ma`, `p_ms`, `p_sa`, `loss_tan_total`, `q_tls`, and `t1_tls_s` from either 2D `dfs`/`ems` fields (reading `vacuum` on `MA`/`SA` and `substrate` on `MS`, matching explicit $2\text{ nm}$ oxide integrals `< 0.03%`) or 3D zero-thickness sheet solves by splitting off a perimeter edge collar (`edge_cutoff` $R_0$, extracting convergent $I_{\text{edge}}$ from the 3D volume energy inside $r < R_0$ combined with Woods et al. 2019 finite-thickness corner weights so results converge `< 1.3%` as $h\to 0$ and vary `< 1.1%` with $R_0$).

## Consequences

* `nyx-analysis` stays lightweight in the dependency graph (production dependencies remain strictly `femcore` and `nyx-physics`; `geom`/`nyx-dfs`/`nyx-ems` appear only in `[dev-dependencies]` for benchmarks, preserving `crates/nyx-engine/tests/layering.rs`).
* Quantitative benchmarks Test 2 (2D & 3D TLS), Test 4 (Park et al. 2026 35 nm disordered Nb $Q_1+R_1$), and Test 5 (native 3D Fluxonium at $\Phi_{\text{ext}} = 0.5\Phi_0$ + $< 6\text{ ms}$ 100-point flux sweep) execute and pass directly in Rust.

## Code

`crates/nyx-analysis/{Cargo.toml,src/{epr,hamiltonian,tls,lib}.rs}`, `crates/nyx-engine/src/analyse.rs`; tests `crates/nyx-analysis/tests/{tls_benchmark,park2026_transmon_benchmark,fluxonium_benchmark}.rs`, `crates/nyx-engine/tests/cqed_validated.rs`.
