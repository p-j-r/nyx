+++
title = "ASME V&V 20 Benchmarks"
description = "Rigorous ASME V&V 20 validation, SQuADDS superconducting circuit benchmarks, Roache GCI convergence, and transparent disclosure of known engine limitations."
template = "benchmarks.html"
+++

`nyx` treats verification and validation as reproducible software tests inside `crates/nyx-engine/tests/validation_*.rs` and `benchmarks/validation/*.json`, gated continuously by `tools/vvuq_ci.ps1` / `tools/vvuq_ci.sh`. Benchmarks are reported with both what passes cleanly and where current architectural limits sit.

## ASME V&V 20 Explicit `AND` Rule

Every validation target defines an experimental or cross-code reference value $D$, input parameter uncertainty $u_{\text{input}}$, reference uncertainty $u_D$, tolerance ceiling $\mathrm{tol} = \mathrm{tol}_{\mathrm{rel}} |D| + \mathrm{tol}_{\mathrm{abs}}$, and coverage factor $k$ (typically $k = 2.0$). Numerical discretisation uncertainty $u_{\text{num}}$ is extracted directly from 2-grid or 3-grid Roache Grid Convergence Index (GCI) refinement. A benchmark passes if and only if both conditions hold simultaneously:

$$|S - D| \le \mathrm{tol}_{\mathrm{rel}} |D| + \mathrm{tol}_{\mathrm{abs}} \qquad \text{AND} \qquad |S - D| \le k \sqrt{u_{\text{num}}^2 + u_{\text{input}}^2 + u_D^2}$$

Where `nyx` supplies deterministic Prager–Synge complementary capacitance brackets $[S_{\text{low}}, S_{\text{up}}]$, the epistemic half-width $\delta_{\text{ep}} = \frac{1}{2}(S_{\text{up}} - S_{\text{low}})$ is **never folded into Gaussian variance**; instead the test enforces $|S_{\text{mid}} - D| \le \delta_{\text{ep}} + k \sqrt{u_{\text{input}}^2 + u_D^2}$.

## SQuADDS cQED Comprehensive Benchmark (`418` Runs)

An unmodified baseline (`9480093a`, zero core patches) was evaluated across `418` runs against Ansys Q3D, Ansys HFSS, `FasterCap`, and three fabricated cryogenic multi-qubit chips (`WM1`, `MUNINN`, and `Dissipator`), completing `410/418` (`98.1%`) overall and `100.0%` across Phase 0 Smoke (`32/32`), Phase 1 Cryogenic Hardware (`30/30`), and Phase 2 Core (`164/164`) within a `19.71 GiB` peak commit ceiling:

- **Planar cQED Electrostatics (`S1.1–S1.3`, `60/60 ok`):** `dfs::electro` (`order = 2`, quarter-point edge grading `AF-2` + spherical-harmonic `far_field(l=2)` `AF-3`) reaches `1.11%` unweighted (`1.04%` Horvitz–Thompson) median self-capacitance error against Ansys Q3D (`0.67%` on `NCap`, `0.59%` on `CapNInterdigitalTee`, `1.53%` on `TransmonCross`), `MAE = 1.369 fF`, `-0.06%` signed bias, and **bitwise-exact (`0.000e0 < 1e-13`) multi-threaded AMG determinism** (`3.08x–3.67x` speedup at `NYX_THREADS=8`).
- **Certified Prager–Synge Bounds & 3-Layer Planar BEM:** Equilibrated $H(\mathrm{div})$ `flux_bounds="midpoint"` strictly brackets **`68/69` (`98.6%`)** of Q3D matrix entries across `S1.1–S1.4`. On narrow-gap `S1.4` geometries, `planar_bem` solves stratified wafers **`76x` faster** (`1.03 s`) at `4.04%` self-capacitance median error.
- **Cryogenic Hardware (`WM1`, `MUNINN`, `Dissipator`):** With mask-verified `+0.9 um` lithographic etch bias, `ems::cavity` PEC (`WM1 W1`) predicts all six readout resonators to **`1.32%` RMSRE** (`2.3x` closer to mK measurements than the `3.0%` published HFSS reference). `MUNINN` reproduces exact experimental coherence ordering (`Ancilla 12.82 ms > Qubit 11.22 ms > Env 8.19 ms`) and solves the full `1,085,076`-tet (`1.402M`-DOF) 3-qubit chip in `82.55 s`. `Dissipator` executes a 100-point DC-SQUID flux sweep (`675`-dim Craig–Bampton Hilbert space) in `1.10 s`.
- **4-Way `QubitCavity` Equivalence (`12/12 ok`):** Closed Craig–Bampton (`M1`), passive-port QEP (`M2`), and 3-port `ems::driven` Schur BBQ + two-pole Kubo (`M3`) agree to `med |Delta f_01| = 1.38e-6%` (`M2`) and `5.0e-6%` (`M3`), `med |Delta alpha| = 0.0000%` (`M2`) and `0.0056%` (`M3`), with Thomas–Reiche–Kuhn oscillator-strength sum rules satisfied to `1.000000000000 +/- 2.7e-10`.

## Honest Disclosures: Known Limitations on Baseline `HEAD`

1. **Finite-Die `FR-11` Thickness Guard (Bug 1):** `crates/nyx-analysis/src/tls.rs` falsely detects vertical `Si|Vac` outer die sidewalls as etched trenches when `sub_box` sits inside a larger vacuum box, blocking `metal_thickness_m` edge corrections on finite-die stacks.
2. **Refined `rank1` Port Exclusion Trade-Off (Bug 2):** Refining `rank1` `lumped_port`/`lumped_jj` sheets places interior port vertices into `cond_vertices`, dropping zero-gradient columns and stalling `krylov_schur`. An unrefined `30 um` port-exclusion box prevents hangs (`163/163` EM `ok`) at the cost of coarse P1 port capacitance (`kappa` HT median `435.75%` vs `41.59%` in `hybrid`).
3. **Multi-Layer `z_coarsen` & Narrow-Gap Steiner Locking (Bugs 3 & 4):** `z_coarsen: true` in `sheet3d.rs` can collapse exterior boundary labels; dense `200–300 nm` GDS fillets across `2–3 um` gaps exhaust `max_steiner_multiplier = 20` in Ruppert 2D refinement (`365.6%` volumetric error on `S1.4`, avoided cleanly by `planar_bem` at `4.04%`).
4. **3D `ems` Parallel `faer` Thread Drift (Bug 5):** Whereas `dfs` AMG is bitwise deterministic (`0.000e0`), 3D `ems` parallel `faer` `LDL^T`/`LU` exhibits `3.59e-10` to `1.15e-7` relative drift between `NYX_THREADS=1` and `8`.
