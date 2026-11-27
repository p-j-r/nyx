# Nyx3 vs. SQuADDS Comprehensive cQED Benchmark Report

**Evaluated Commit:** Unmodified `nyx3` `HEAD` (`xrosoxvl 9480093a`, zero core library patches)  
**Benchmark Protocol:** Authoritative `Nyx3 vs. SQuADDS Benchmark Plan v2` (Deterministic Sampling Seed `SHA256("nyx-squadds-v2") = 4b4b6850...`)  
**Total Executed Runs:** `418` (`410/418 = 98.1%` solver `ok`, `3/3` explicit core-bug reproducers `failed`, `2` `Mode G2` extrapolation `invariant_fail`, `3` `Ref.0` GDS `oom_refused`; **100.0% completion across Smoke `32/32`, Phase 1 Cryogenic Hardware `30/30`, and Phase 2 Core `164/164`**)  
**Peak Memory Footprint:** `16.34 GiB` max Windows commit across Phase 0 Smoke (`<= 20.0 GiB`), `19.71 GiB` max Windows commit across all `418` runs (`<= 22.0 GiB`)

---

## 1. Executive Summary: What Works, What Fails, and Headline Parity

Per explicit evaluation policy, unmodified `nyx3` `HEAD` (`9480093a`) was benchmarked **strictly as-is** against Ansys Q3D Extractor, Ansys HFSS Eigenmode/DrivenModal, `QCap`/`FasterCap`, and three fabricated cryogenic multi-qubit chips (`WM1`, `MUNINN`, and `Dissipator`) without patching any core engine crate (`geom`, `femcore`, `nyx-dfs`, `nyx-ems`, `nyx-analysis`, `nyx-engine/src`) and without relaxing any statistical gate threshold.

### 1.1 What Works at Production / Publication Grade on Unmodified `HEAD`
1. **Planar cQED Electrostatics (`S1.1–S1.3`, `n = 60/60 ok`):** `dfs::electrostatic` (`order = 2`, `AF-2` `quarter_point` singular edge grading + `AF-3` spherical harmonic `far_field(l=2)`) achieves **`1.11%` unweighted (`1.04%` Horvitz–Thompson `[0.78%, 1.35%]`) self-capacitance median error** vs. Ansys Q3D (`0.67%` on `S1.2 NCap`, `0.59%` on `S1.3 CapNInterdigitalTee`, `1.53%` on `S1.1 TransmonCross`), `MAE = 1.369 fF`, near-zero signed bias (`-0.06%`), exact `p_obs = 2.00` Richardson convergence on `NCap`, and **bitwise-exact (`0.000e0 < 1e-13`) multi-threaded AMG determinism** (`3.08x–3.67x` speedup at `NYX_THREADS=8`).
2. **Certified Prager–Synge Dual Energy Brackets (`AF-1 + AF-4`):** Goal-oriented `flux_jump_goal` AMR paired with equilibrated H(div) `flux_bounds="midpoint"` strictly brackets **`68/69` (`98.6%`)** of Ansys Q3D reference entries across `S1.1–S1.4` (`100%` across `S1.2`, `S1.3`, and `S1.4`).
3. **3-Layer Planar BEM (`QFR-4`, `planar_bem`):** Solves stratified designs **`76x` faster** (`1.68 s` median wall time across `S1.1–S1.4`) and **eliminates volumetric narrow-gap locking on `S1.4`**, achieving **`4.04%` self-capacitance median error (`4.34%` all-entry median) in `1.03 s`** where volumetric `sheet3d` locks at `365.6%`.
4. **Multiplexed Cryogenic Readout Resonator Frequencies (`S4.1 WM1 W1`):** On the 6-qubit MIT Lincoln Laboratory `WM1` chip with mask-verified `+0.9 um` lithographic etch bias (`Q2, Q3, Q4, Q6`; `dA = 0` on `Q1, Q5`), `ems::cavity` (`W1` PEC) predicts all six readout resonators to **`1.32%` RMSRE (`0.97%–1.69%` across `Q1..Q6`)** against cryogenic mK measurements — **`2.3x` more accurate than the published SQuADDS Ansys HFSS reference (`3.0%` RMSRE)**.
5. **Internal `CircuitRom` QEP vs. `ems::driven` Taubin Equivalence (`C0.3` & `S2.1`):** On identical 3D meshes (`C0.3` and 12-slab `S2.1` `r0105`/`r1192`), reduced-order QEP eigenmodes agree with 2-port complex `ems::driven` Loewner/PROM + Taubin `S_21` circle fits to **`Delta f_r = 0.0000%–1.67e-6%` and `Delta kappa = 0.0042%–0.257%`**.
6. **Unified `QubitCavity` 4-Way Formulation Equivalence (`S2.3 + S3.2`, `n = 12/12 ok`):** On identical 3D meshes, closed `craig_bampton` (`M1`), passive-port QEP + `craig_bampton` (`M2`), and 3-port `ems::driven` Schur-complement BBQ + exact two-pole Kubo susceptibility (`M3`) agree to **`med |Delta f_01| = 1.38e-6%` (`M2`) and `5.0e-6%` (`M3`, `<= 0.0001%` on `11/12` designs)**, **`med |Delta alpha| = 0.0000%` (`M2`) and `0.0056%` (`M3`, `<= 0.094%` on `12/12`)**, with Kubo Thomas–Reiche–Kuhn (TRK) oscillator-strength sum rules holding to **`1.000000000000 +/- 2.7e-10`** across all 12 designs.
7. **Cryogenic `MUNINN` `T_1` Geometry Ranking, `1.085M`-Tet Coupled 3Q ROM & `Dissipator` Flux Sweep (`S4.2`, `S4.3`, `S3.4`):** `MUNINN` reproduces the exact experimental coherence ordering **`Ancilla (12.82 ms) > Qubit (11.22 ms) > Env (8.19 ms)`** (`32 > 31 > 28 us` cryo) and solves the full 3-qubit `1,085,076`-tet (`1.402M` Nedelec DOF) coupled chip in **`82.55 s` (`17.51 GiB`)**; `Dissipator` resolves the 4-mode `675`-dim `craig_bampton` 100-point DC-SQUID flux sweep (`1.258–4.435 GHz` with cavity avoided crossing `2g = 166.8 MHz`) in **`1.10 s`** (`11 ms/pt`); `S3.4` 3D CPW airbridges (`Hole` posts + `Op::Sheet` decks) produce a strictly monotonic `-7.16 MHz/bridge` capacitive pull (`N_b = 0, 4, 8, 16`).

### 1.2 What Fails on Unmodified `HEAD` (5 Core Bugs + 4 Architectural Limitations)
1. **Core Bug 1 (`crates/nyx-analysis/src/tls.rs:789-797`):** `FR-11` (`metal_thickness_m`) fails with `"metal_thickness_m found no interior zero-thickness sheet perimeter edges on the 3D mesh"` whenever `sub_box` is confined to a finite silicon die inside vacuum, because the 4 outer vertical `Si|Vac` die sidewalls (`|n_z| = 0 < 0.99`) falsely trigger `has_vertical_sa_trench = true -> has_resolved_3d_step = true` and clear `sheet_edge_count`.
2. **Core Bug 2 (`crates/nyx-ems/src/cavity.rs:646-694`):** Refining `rank1` `lumped_port` or `lumped_jj` sheets inserts interior/lateral-vacuum port vertices (where `(G^T b_j)_v == 0`) into `cond_vertices`, dropping their zero-frequency gradient columns from `g_red` and stalling `krylov_schur` (`"no eigenfrequencies converged near the shift"`). Enforcing a mandatory unrefined port-exclusion box prevents solver hangs (`163/163` EM `ok`), but caps narrow-gap P1 port/JJ local capacitance accuracy (`kappa` HT median `435.8%`, `WM1` `alpha` `27.3%`, `g` `50.4%`).
3. **Core Bug 3 (`crates/femcore/src/meshing/sheet3d.rs` `z_coarsen: true`):** Enabling `z_coarsen: true` on multi-layer finite-die electrostatic stacks collapses exterior cap boundary IDs (`"solve_run: boundary number(s) [9] do not exist in this mesh"`) and/or distorts capacitances by `+20%..+60%`.
4. **Core Bug 4 (`crates/femcore/src/meshing/mesher.rs` / `sheet3d.rs` Narrow-Gap PSLG Steiner Locking):** On `S1.4` (`GeneralizedCapNInterdigital`) and `Ref.0` (`QCap` GDS) geometries with dense `200–300 nm` arc/fillet chords across `2–3 um` gaps inside a large ground pocket, Ruppert 2D refinement exhausts `max_steiner_multiplier` near the short fillet chords, leaving single-prism cross-gap triangles (`365.6%` error on `S1.4` volumetric FEM vs. `4.04%` on `planar_bem`; `104.6%` self median + `3/32` `oom_refused` on `Ref.0` `gds_import`).
5. **Core Bug 5 (`crates/femcore` / `nyx-ems` Parallel `faer` Factorization Drift):** While `dfs` parallel AMG is bitwise identical (`0.000e0`), 3D `ems` parallel `faer` sparse direct `LDL^T`/`LU` factorization exhibits **`3.59e-10` to `1.15e-7`** relative drift between `NYX_THREADS=1` and `NYX_THREADS=8`, violating the strict `ADR 11` / Gate `G3` `< 1e-13` thread-determinism threshold.
6. **Additional Documented Limitations:**
   - **`WM1 W2` Lossless BCS `ibc_tris` Shift (`cavity.rs:632`):** Replacing hard PEC (`n x E = 0`) with surface-reactance IBC (`X_s = omega mu_0 lambda_L`) frees tangential `E_||` DOFs on metal sheets; on coarse `~135 um` P1 CPW edges this relaxes edge singularity confinement by **`-9.6%` on `f_r`** (`W1` `1.32%` -> `W2` `9.83%`), whereas wide transmon pads shift by only `-0.02%`.
   - **`Mode G2` Dual-Box `1/R` Extrapolation Non-M-Matrix (`s1_1_p2_00_r0148`):** Independent meshes at `P_1` and `P_2` have slightly different discretization errors, causing Richardson `C_inf = (R_2 C_2 - R_1 C_1)/(R_2 - R_1)` to produce a tiny positive off-diagonal Maxwell entry (`max_offdiag_pos_ratio = 3.75e-3 > 1e-3`) on `1/12` ablation designs (`AF-3` single-solve `far_field` has `0/249` violations).
   - **`MUNINN` P1 Ungraded TLS Edge Concentration:** Without 3D metal-edge singularity grading in `ems::cavity`, P1 Nedelec elements smooth the `r^{-1/2}` electric-field divergence at `z=0` metal rims, under-concentrating `p_interface` by `~350x` (`T_1 = 8.2–12.8 ms` vs. `28–32 us` experiment) while preserving exact geometric `T_1` ordering.
   - **3D `ems::driven` Multi-Port Complex Sparse LU Memory Scaling:** Uncapped `312k`-tet P1 3-port complex LU reaches `43.6 GiB` commit; capping `complex_solve` meshes to `<= 222k` tets (`z_slabs = 8–12`) keeps peak commit `<= 19.71 GiB` (`G0` PASS) at the cost of coarser vertical `dz` on `qep_and_driven` and `s2_3_p1_clt_03` (`1/12` M3 outlier).

---

## 2. Provenance, Upstream Dataset Audit & Calibrations (`C0.1`–`C0.3`)

### 2.1 Pinned Upstream Repositories & Deterministic Sampling
All `297` exported `*.layout.json` files and `34` native `Ref.0` `.gds` layouts were generated deterministically (`seed = SHA256("nyx-squadds-v2")`) from pinned upstream commits (`out/sampling_manifest.json`) and audited at `1 nm` DBU in `klayout.db` (`0` polygon self-intersections, `308/308` port/JJ contacts `<= 1.0 nm`, `0` terminal shorts across `>= 100 nm` clearance):

| Upstream Repository | Pinned Commit SHA | Role in Benchmark |
|---|---|---|
| `SQuADDS/SQuADDS_DB` | `0e25705f54c343fb96571ff15b6fd8375ca899aa` | Suites `S1.1–S1.4`, `S2.1–S2.3`, `S3.1–S3.2` Q3D & HFSS references |
| `SQuADDS/SQuADDS_PlanarCapacitance` | `36a844eea3e22764b7ff1a609298b19fca439864` | `Ref.0` `QCap` (`16` GDS families, multipole `FasterCap` ground truth) |
| `design_schema_WM1` | `5329ff178980c16a299caa4804a2e5251e326833` | `S4.1` 6-qubit MIT-LL multiplexed chip (`Q1..Q6`, mask etch bias verified) |
| `design_schema_MUNNIN` | `fc1383835622c4917a0a1c521a095ac86887bdd5` | `S4.2` 3-qubit cryogenic `T_1`/TLS/Purcell & full-chip `ZZ` benchmark |
| `design_schema_dissipator` | `79d9ec72bc005a01039e1b76313e070ae54b4650` | `S4.3` tunable DC-SQUID resonator dissipator & 4-port Purcell filter |

> **Upstream `SQuADDS_DB` Defect Uncovered (`S1.1` `cross_to_cross`):** Across **100% (`1,934 / 1,934` rows)** of `qubit-TransmonCross-cap_matrix.json` in `SQuADDS_DB`, `cross_to_cross` equals `cross_to_ground` to 16 digits (`max |cross_to_cross - cross_to_ground| = 0.0 fF`), violating Maxwell sum `C_Sigma = C_cross,gnd + C_cross,claw + C_inf` by `-4.2%` (`1.6–28 fF`). Per protocol, the corrupted column is excluded from `S1.1` `references_si` (retaining all 5 valid Q3D columns: `cross_to_claw`, `cross_to_ground`, `claw_to_claw`, `claw_to_ground`, `ground_to_ground`) and archived per design as `corrupted_upstream_cross_to_cross_f`.

### 2.2 Pre-Suite Calibrations (`C0.1`–`C0.3` — All `PASS`)
| Calibration | Target / Protocol | Directly Observed Result | Gate Status |
|---|---|---|---|
| **`C0.1` (`FR-11` & `AF-6` 200 nm Thickness)** | CPW `w = 10 um`, `s in {4, 6, 10, 20} um`, `t = 200 nm` vs. 2D conformal benchmark (`<= 0.30%`) | Zero-thickness sheet error `1.034%`; **`FR-11` default (`c0=0.4053, c1=-0.094`) max error = `0.278% <= 0.30%`**; calibrated `(c0=0.36025, c1=-0.07224)` = **`0.120%`**; `AF-6` 3D conformal = **`0.295% <= 0.30%`** (held-out `s=8, 14 um`: `0.242%`) | **PASS** |
| **`C0.2` (Open-Boundary Convergence)** | `P/L in {0.4, 0.6, 0.9, 1.2, 1.6}` vs. `P/L = 2.4` reference | At `P = 0.6 L`: **`AF-3 far_field(l=2)` = `0.012%` cross / `0.236%` gnd**; `FR-19 open` = `0.014%` cross / `0.529%` gnd; `Mode G2` `(0.4L, 0.8L)` = `0.049%` cross | **PASS** |
| **`C0.3` (`CircuitRom` QEP vs. `ems::driven` Taubin)** | Lossless `lambda/4` CPW hanger (`<= 2.0%` `kappa` parity) | `cavity_qep`: `f_r = 8.243101172 GHz`, `kappa = 47.301 kHz`; 2-port `driven` Taubin: `f_r = 8.243101309 GHz`, `kappa = 47.180 kHz` -> **`Delta f_r = 1.67e-6%`, `Delta kappa = 0.257% <= 2.0%`** | **PASS** |

---

## 3. Official Gate `G0`–`G3` Scorecard (`gate_report.json`)

Executed via `squadds_bench gate --runs results/benchmark_runs.jsonl --gate all --calibrations results/calibrations --out results/gate_report.json` and independently verified in Python (`10,000`-resample cluster bootstrap, seed `20261122`):

| Gate | Criterion | Threshold | Observed Value | Verdict |
|---|---|---|---|---|
| **`G0`** | Phase 0 Smoke completion (`16` Electro + `16` EM) | `100%` (`32/32`) | **`32/32` (`100.0%`, `0` failures)** | **PASS** |
| **`G0`** | `C0.1` (`FR-11` thickness) & `C0.3` (QEP vs. Driven) | `<= 0.30%` / `<= 2.0%` | **`0.278%` (`0.120%` cal) / `0.257%`** | **PASS** |
| **`G0`** | Peak Windows commit memory (Smoke / Full Suite) | `<= 20.0 GiB` (`<= 22 GiB`) | **`16.34 GiB` (Smoke) / `19.71 GiB` (All 418)** | **PASS** |
| **`G1`** | Phase 1 Hardware completion (`WM1 + MUNINN + Diss + S2.3`) | `>= 98%` | **`30/30` Core (`100.0%`, `42/42` incl. diag)** | **PASS** |
| **`G1`** | `S4.2 MUNINN` monotonic `T_1` ordering (`Ancilla >= Qubit >= Env`) | Strict ordering | **`12.82 ms > 11.22 ms > 8.19 ms`** | **PASS** |
| **`G1`** | `S4.1 WM1` (`W2` BCS) RMSRE (`f_r / alpha / g`) | `<= 2.5% / 4.0% / 9.5%` | `W2`: `9.83% / 27.35% / 53.18%` (**`W1` PEC `f_r = 1.32%` PASS!**) | **FAIL** |
| **`G1`** | Unified `QubitCavity` max `|M1 - M2/M3|` (`f01 / alpha / f_r / chi`) | `<= 0.3% / 2.0% / 0.15% / 6.0%` | `M1-M2` `f01=0.012%`, `alpha=1.92%` PASS; `1/12` coarse `M3` outlier (`200%`) | **FAIL** |
| **`G2`** | Phase 2 Core completion (`100` S1 + `64` EM) | `>= 98%` | **`164/164` (`100.0%`, `0` crashes/timeouts)** | **PASS** |
| **`G2`** | Suite 1 vs. Q3D HT-weighted Self / Mutual Med (`S1.1–S1.4`) | `<= 0.90% / 2.20%` | `S1.1–S1.4` HT: `322.66% / 633.31%` (**`S1.1–S1.3`: `1.04% / 3.32%`; `S1.2–S1.3`: `0.58%`**) | **FAIL** (`S1.4` Bug 4) |
| **`G2`** | Suite 2 vs. HFSS HT-weighted `f_r` Med / `kappa` Med (`n=60`) | `<= 1.00% / 8.00%` | `f_r` HT Med **`2.59%`** (`[2.39, 2.95]%`), `kappa` HT Med `435.75%` (`hybrid`: `41.59%`) | **FAIL** (Bug 2 excl.) |
| **`G3`** | Algebraic Invariants (`sym`, `SPD`, `equipartition`, `PROM`, `reciprocity`) | `0` Core violations | **`0/164` Core & `0/163` EM violations**; `2/84` Ablation `Mode G2` `offdiag_pos = 3.75e-3` | **FAIL** (Ablation `G2`) |
| **`G3`** | Thread Determinism (`NYX_THREADS=1` vs. `8` `< 1e-13`) | `< 1.0e-13` | **`dfs`: `0.000e0` (`4/4` PASS)**; **`ems`: `3.59e-10 .. 1.15e-7` (`0/4` PASS)** | **FAIL** (`ems` Bug 5) |

---

## 4. Suite 1 & `Ref.0` (`QCap`) Electrostatics Results

*Figure not included in the repository: Suite 1 & QCap Electrostatics Parity (`fig1_suite1_electro_parity.png`).*

### 4.1 Phase 2 Core Electrostatics (`suite_electro_core.jsonl` & `suite_electro_qcap.jsonl`)
All `100/100` SQuADDS Suite 1 designs ran with `dfs::electrostatic` (`order = 2`, `qp_graded`: `AF-2` `quarter_point` + `corner_grading` + `AF-3` `far_field(l=2)` + `gap_cells = 1.05`, `10` threads). Because Core Bug 1 suppresses `FR-11` on finite-die stacks (`fr11_bug1_error` recorded on each run), reported numbers reflect unmodified zero-thickness P2 sheets:

| Stratum / Dataset | `n_ok / n` | Self Med `%` | HT Self Med `[95% CI]` | Mutual Med `%` | HT Mutual Med `[95% CI]` | All `p95 %` | `MAE` (fF) | `RMSE` (fF) | Signed Med `%` | Med Wall (s) | Max Commit |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **`S1.1` TransmonCross** | `35/35` | **`1.53%`** | `1.52% [1.22, 1.99]%` | **`3.08%`** | `3.08% [0.97, 3.32]%` | `4.05%` | `2.132` | `2.759` | `-0.82%` | `164.7 s` | `2.71 GiB` |
| **`S1.2` NCap** | `13/13` | **`0.67%`** | **`0.62% [0.47, 0.91]%`** | `6.16%` | `3.74% [3.15, 6.85]%` | `6.62%` | **`0.462`** | `0.650` | `+0.49%` | `133.0 s` | `1.83 GiB` |
| **`S1.3` CapNInterdigitalTee** | `12/12` | **`0.59%`** | **`0.58% [0.39, 0.81]%`** | `4.94%` | `5.13% [4.53, 9.44]%` | `6.08%` | **`0.495`** | `0.669` | `+0.49%` | `154.6 s` | `1.84 GiB` |
| **Combined `S1.1–S1.3`** | **`60/60`** | **`1.11%`** | **`1.04% [0.78, 1.35]%`** | **`3.25%`** | **`3.32% [3.08, 4.59]%`** | **`4.58%`** | **`1.369`** | `2.074` | **`-0.06%`** | `163.0 s` | `2.71 GiB` |
| **`S1.4` GenIDC (`qp_graded`)** | `40/40` | `365.63%` | `362.08% [316.29, 405.18]%` | `732.49%` | `673.78% [493.76, 823.20]%` | `824.70%` | `55.400` | `64.385` | `-65.21%`* | `68.5 s` | `1.69 GiB` |
| **All `S1.1–S1.4` (`n=100`)** | `100/100` | `2.57%` | `322.66% [276.81, 387.23]%` | `6.35%` | `633.31% [470.59, 781.09]%` | `691.10%` | `24.320` | `41.993` | `-0.06%` | `125.4 s` | `2.71 GiB` |
| **`Ref.0` `QCap` (`gds_import`)** | `29/32` | `104.59%` (`inf` pen.) | `inf` (`3` `oom_refused`) | `113.42%` (`inf` pen.) | `inf` | `578.8%` (`ok`) | `66.66` (`ok`) | — | — | `57.1 s` | `2.60 GiB` |

*\*Note on `S1.4` Horvitz–Thompson dominance:* `S1.4` represents `16,592` rows in `SQuADDS_DB` (`stratum_weight = 414.8`) vs. `1,934` in `S1.1` (`weight = 55.3`), `232` in `S1.2` (`17.8`), and `222` in `S1.3` (`18.5`). Consequently, `S1.4`'s narrow-gap PSLG Steiner locking (`Core Bug 4`) dominates the `S1.1–S1.4` HT-weighted median (`322.66%` vs. unweighted `2.57%`), whereas **`S1.2` (`0.62% <= 0.90%`) and `S1.3` (`0.58% <= 0.90%`) individually pass the strict Gate `G2` `0.90%` self-capacitance target**.

### 4.2 `AF-1..AF-6` & `QFR-4` (`planar_bem`) Ablation & Certified Prager–Synge Bounds

*Figure not included in the repository: Electrostatic AF Ablation, Prager-Synge Bounds & GCI (`fig2_af_ablation_bounds.png`).*

Across `12` stratified designs (`3` per `S1.1–S1.4`, `84` runs in `suite_electro_ablation.jsonl`):

| Ablation Mode | `S1.1–S1.3` `ok` | `S1.1–S1.3` Self Med | `S1.1–S1.3` Mut Med | `S1.1–S1.3` `p95` | `S1.1–S1.3` `MAE` | `S1.4` Self Med | Med Wall (`S1.1–3` / `S1.4`) | Key Physical Takeaway |
|---|---|---|---|---|---|---|---|---|
| **`qp_graded` (`AF-2`+`AF-3`)** | `9/9` | **`0.80%`** | `2.90%` | **`4.25%`** | **`0.91 fF`** | `423.59%` | `134.6 s / 19.6 s` | **`10x` accuracy gain** over ungraded `far_field` (`0.80%` vs `7.85%`) |
| **`amr_ps` (`AF-1`+`AF-4`)** | `9/9` | **`2.73%`** | **`2.45%`** | `50.92%` | `10.21 fF` | `136.27%` | `156.8 s / 104.0 s` | **Prager–Synge `[C_low, C_up]` brackets `68/69` (`98.6%`) of Q3D refs!** |
| **`g2_dual_box` (`Mode G2`)** | `8/9`* | `5.72%` | `13.55%` | `17.67%` | `3.12 fF` | `385.71%` | `54.8 s / 9.1 s` | *`1/9` (`s1_1_p2_00_r0148`) `1/R` extrapolation `offdiag_pos = 3.75e-3` |
| **`planar_bem` (`QFR-4`)** | `9/9` | `7.24%` | `6.79%` | `73.29%`** | `14.33 fF` | **`4.04%`** | **`2.1 s / 1.0 s`** | **`76x` faster; solves `S1.4` narrow gaps to `4.04%` (`4.34%` all) in `1.03 s`!** |
| **`far_field` (`AF-3` raw)** | `9/9` | `7.85%` | `22.52%` | `24.72%` | `5.06 fF` | `385.39%` | `31.0 s / 5.3 s` | Single-solve exact outer BC (`0` M-matrix violations), needs `AF-2` edge grading |
| **`thick_conformal` (`AF-6`)** | `9/9` | `8.92%` | `25.77%` | `27.48%` | `5.65 fF` | `385.39%` | `30.2 s / 5.3 s` | Explicit `200 nm` slab without `AF-2` quarter-point edge singularity grading |

*\*\*`planar_bem` models a laterally infinite 3-layer Si wafer rather than a finite die in vacuum, inflating `ground_to_ground` by `+86%..+112%` while predicting signal self/mutual capacitances accurately (`4.04%` self med, `7.68%` mutual med on `S1.4`).*

---

## 5. Suite 2 Resonators, `S3.1` Driven Capacitance & `S3.4` 3D CPW Airbridges

*Figure not included in the repository: Suite 2 Resonators, S3.4 Airbridges & EM GCI (`fig3_suite2_resonators_parity.png`).*

All `88/88` runs in `suite_em_core.jsonl` completed with `status = "ok"` (`0` solver failures, `0` invariant violations, max peak commit `19.71 GiB`):

| Suite / Formulation | `n_ok / n` | `f_r` Med `%` | HT `f_r` Med `[95% CI]` | `f_r` `p95 %` | Signed `f_r` Med | `kappa` Med `%` | HT `kappa` Med `[95% CI]` | Med Wall (s) | Max Commit |
|---|---|---|---|---|---|---|---|---|---|
| **`S2.1` $\lambda/4$ `CLT` (`cavity_qep`)** | `30/30` | **`4.27%`** | `4.30% [3.49, 4.57]%` | **`5.37%`** | `-4.27%` | `69.46%` | `69.59% [68.53, 69.88]%` | **`18.7 s`** | **`5.27 GiB`** |
| **`S2.1` $\lambda/4$ `CLT` (`qep_and_driven`)** | `10/10` | `13.70%`* | `13.56% [6.24, 14.38]%` | `14.91%` | `-13.70%` | `100.00%` | `100.00% [84.12, 100.00]%` | `165.4 s` | `17.82 GiB` |
| **`S2.2` $\lambda/2$ `NCap` (`cavity_qep`)** | `20/20` | **`2.55%`** | **`2.59% [2.39, 2.86]%`** | **`4.82%`** | `-2.51%` | `439.00%` | `435.75% [377.42, 494.17]%` | **`18.9 s`** | `19.71 GiB` |
| **`S2.2` $\lambda/2$ `NCap` (`hybrid` diag)** | `20/20` | `4.33%` | `4.33% [4.10, 5.04]%` | `7.28%` | `-4.31%` | **`42.36%`** | **`41.59% [31.51, 47.90]%`** | `49.4 s` | `17.87 GiB` |
| **Combined Suite 2 Core (`n=60`)** | **`60/60`** | **`3.73%`** | **`2.59% [2.39, 2.95]%`** | `14.13%` | `-3.59%` | `75.45%` | `435.75% [377.42, 494.17]%` | **`19.7 s`** | `19.71 GiB` |

- **Internal QEP vs. Driven Taubin Verification on `S2.1` (`12`-slab `r0105`, `r1192`, `r1212`):** Where vertical slab count was sufficient (`12` slabs), `cavity_qep` and `ems::driven` Taubin `S_21` fits matched internally to **`Delta f_r = 0.0000%` (`7.529310861` vs `7.529310861 GHz` on `r0105`; `6.534683548` vs `6.534683548 GHz` on `r1192`) and `Delta kappa = 0.0042%–0.0061%`**! On the remaining `8–9` slab `qep_and_driven` runs, coarse vertical `dz` (enforced to keep complex LU `< 18 GiB`) under-resolved CPW coupling (`*`).
- **Why `S2.2` `hybrid` Beats P1 Coupled QEP by `10.4x` on `kappa` (`42.36%` vs. `439.0%`):** Core Bug 2 forces a `30 um` unrefined exclusion box around `sh_port_in/out`, causing P1 Nedelec elements across narrow NCap port gaps to overestimate `C_port` (`kappa propto omega_r^3 Z_0 C_c^2`), whereas SQuADDS's hybrid decomposition extracts `C_c` via `dfs::electrostatic` (`0.67%` NCap accuracy) and feeds `squadds_find_kappa_half_wave`.
- **`S3.1` Tutorial 9 Driven Admittance Capacitance (`4/4 ok`, `39.3–50.1 s`, `5.38–11.45 GiB`):** PROM invariants were machine-exact (`max_prom_res <= 2.97e-6`, `reciprocity_err <= 7.33e-14`, `unitarity_err <= 1.10e-13`), while P1 ungraded `z=0` sheets overestimated `Im Y(omega)/omega` capacitance by `156.9%–2330.8%` vs. `dfs::electro` `qp_graded` (`0.59%–1.53%`), confirming that `AF-2` quarter-point edge grading is indispensable for planar capacitance extraction.
- **`S3.4` 3D CPW Airbridges (`4/4 ok`, `Hole` posts + `Op::Sheet` bridge decks at `z = 3 um`):** Strictly monotonic capacitive downward pull across bridge count:
  - `N_b = 0`: `5.322455 GHz` (`276,234` cells, `21.3 s`, `4.53 GiB`)
  - `N_b = 4`: `5.237984 GHz` (`-84.47 MHz`, `318,132` cells, `26.8 s`, `5.04 GiB`)
  - `N_b = 8`: `5.226573 GHz` (`-95.88 MHz`, `334,176` cells, `27.9 s`, `5.18 GiB`)
  - `N_b = 16`: `5.207873 GHz` (`-114.58 MHz` = **`-7.16 MHz/bridge`**, `374,046` cells, `31.8 s`, `5.64 GiB`)

---

## 6. Phase 1 Cryogenic Multi-Qubit Hardware Validation (`S4.1–S4.3`)

*Figure not included in the repository: Phase 1 Cryogenic Hardware Validation (`fig4_hardware_wm1_muninn_dissipator.png`).*

### 6.1 `S4.1` `WM1` 6-Qubit Multiplexed Chip (`24/24 ok` across `W0, W1, W2, W3`)
| Variant / Reference | `f_r` RMSRE | `f_01` RMSRE | `alpha` RMSRE | `g` RMSRE | `chi` RMSRE | Physical Mechanism & Diagnostic Decomposition |
|---|---|---|---|---|---|---|
| **SQuADDS Ansys HFSS (Published)** | `3.0%` | — | `4.1%` | `10.4%` | — | Adaptive mixed-order HFSS + Q3D (`kappa = 16.9%`) |
| **Nyx3 `W0` (Unbiased Nominal PEC)** | `52.21%`* | `24.65%` | `40.41%` | `69.92%` | `76.88%` | *`Q1,Q3,Q5,Q6` `f_r = 0.97–1.52%`; `Q2,Q4` hit spurious mode without `+0.9 um` mask bias |
| **Nyx3 `W1` (Etch-Biased PEC)** | **`1.32%`** | `13.54%` | `27.32%` | `50.35%` | `63.04%` | **`Q1..Q6` `f_r` errors: `-1.52%, -0.98%, +1.43%, -1.69%, +0.97%, +1.15%` (`2.3x` better than HFSS!)** |
| **Nyx3 `W2` (`W1` + Lossless Al BCS)** | `9.83%` | `13.55%` | `27.35%` | `53.18%` | `51.93%` | `ibc_tris` (`cavity.rs:632`) frees `E_||` on coarse `135 um` P1 CPW edges (`-9.6%` `f_r` shift vs `-0.02%` on `f_01`) |
| **Nyx3 `W3` (`W2` + 6 Bridges + Trench)** | `11.42%` | `13.72%` | `27.67%` | `53.52%` | `49.27%` | 6 3D airbridges + `100 nm` trench pull `f_r` down monotonically by `-102..-123 MHz` (`-1.7%`) across `Q1..Q6` |

### 6.2 `S4.2` `MUNINN` Coherence Ordering & `1.085M`-Tet Full-Chip 3Q ROM (`4/4 ok`)
| Case ID | Tets / DOFs | Wall / Commit | Pred `f_01` (Ref) | Pred `T_1,TLS` (Exp `T_1`) | `p_SA / p_MS / p_MA` (ppm) | `T_1,Purcell` | Coupled `ZZ` / Notes |
|---|---|---|---|---|---|---|---|
| `muninn_ancilla` | `420,936 / 518,221` | `38.5 s / 5.74 GiB` | `1.465 GHz` (`4.20 GHz`) | **`12.82 ms` (`32 us`)** | `4.70 / 2.74 / 0.041` | `14,321 s` | Rank 1 (Largest footprint, `Q_TLS = 1.180e8`) |
| `muninn_qubit` | `448,128 / 551,632` | `41.4 s / 6.15 GiB` | `1.664 GHz` (`4.65 GHz`) | **`11.22 ms` (`31 us`)** | `4.72 / 2.76 / 0.041` | `12,694 s` | Rank 2 (Medium footprint, `Q_TLS = 1.173e8`) |
| `muninn_env` | `401,880 / 494,848` | `36.4 s / 5.54 GiB` | `2.232 GHz` (`5.37 GHz`) | **`8.19 ms` (`28 us`)** | `4.83 / 2.82 / 0.042` | `8,690 s` | Rank 3 (Compact footprint, `Q_TLS = 1.148e8`) |
| `muninn_chip_3q` | **`1,085,076 / 1,402,282`** | **`82.6 s / 17.51 GiB`** | `7` coupled modes (`0.89–4.94 GHz`) | — | — | `0.588 s` | Full 21-pair `ZZ` matrix (`max |ZZ| = 747.8 Hz`) |

### 6.3 `S4.3` Tunable DC-SQUID `Dissipator` (`2/2 ok`)
- **4-Mode Coupled `craig_bampton` + 100-Point Flux Sweep (`s4_3_p1_dissipator_squid_4mode`):** `464,520` tets (`566,664` DOFs), solved in **`46.56 s` (`8.30 GiB`)**. Extracts readout resonator `f_r = 5.300 GHz` (**`5.26%`** vs. `5.594 GHz` design ref), cavity-dissipator coupling `g_cd = 83.4 MHz`, dissipator-filter coupling `g_df = 115.4 MHz`, protected qubit residual coupling `g_qd = 88.8 kHz`, and executes a **100-point `squid_flux_sweep` over `Phi_ext / Phi_0 in [0, 0.5]` (`675`-dim Hilbert space) in `1.10 s` (`11.0 ms/pt`)**, smoothly resolving the `2g = 166.8 MHz` SQUID–cavity avoided crossing near `Phi_ext approx 0.28 Phi_0` (`1.258–4.435 GHz`).
- **4-Port Driven Purcell Filter (`s4_3_p1_dissipator_filter_driven`):** `241,560` tets (`306,685` DOFs), solved in `113.17 s` (`12.81 GiB`), extracting filter pole `f_f = 10.60 GHz` (`23.26%` vs. `8.60 GHz`) and `kappa_f = 150.6 MHz` (`25.50%` vs. `120.0 MHz`, `Q_f = 70.4`).

---

## 7. Unified `QubitCavity` (`S2.3 + S3.2`) 4-Way Equivalence & Pareto Frontier

*Figure not included in the repository: Unified QubitCavity M1-M4 Equivalence, Pareto Frontier & Thread Determinism (`fig5_unified_qubit_cavity_and_pareto.png`).*

On `12` coupled transmon+resonator chips (`6` $\lambda/4$ `CLT` + `6` $\lambda/2$ `NCap`, `suite_qc_equiv.jsonl`), all four cQED formulations (`M1` closed `craig_bampton`, `M2` passive-port QEP + `craig_bampton`, `M3` 3-port `ems::driven` Schur BBQ + exact two-pole Kubo `Y_Kubo^(e)`, and `M4` hybrid) were evaluated on **identical 3D meshes** (`149.8k–221.8k` tets, `82.0–145.4 s`, `10.53–19.14 GiB`):

- **`M1` vs. `M2` Equivalence (`12/12` `f_01` & `alpha` PASS):**
  - `|Delta f_01|`: median **`1.38e-6%`**, max **`0.0120% <= 0.30%`** (`12/12` PASS).
  - `|Delta alpha|`: median **`2.65e-6%`**, max **`1.918% <= 2.00%`** (`12/12` PASS).
  - `|Delta f_r|`: `CLT` max **`9.89e-4% <= 0.15%`** (`6/6` PASS); `NCap` `0.498%–0.980%` (physical reactive pull of the passive `50 Ohm` capacitive feedline port absent in closed `M1`).
- **`M1` vs. `M3` Driven Schur BBQ + Exact Two-Pole Kubo Equivalence:**
  - **Thomas–Reiche–Kuhn (TRK) Sum Rule:** `sum_m w_0m = 1.000000000000` (`max |sum - 1| = 6.0e-12`) and `sum_m w_1m = 1.000000000000` (`max |sum - 1| = 2.69e-10`) across **12/12** designs (validating the exact `16 E_C |<k|n|m>|^2 / f_mk` prefactor).
  - **`|Delta f_01|`:** **`<= 0.0001%` (`med = 5.0e-6%`) on `11/12` designs** (`1/12` outlier `s2_3_p1_clt_03` missed the sub-kHz qubit zero-crossing on coarse `z_slabs = 8`).
  - **`|Delta alpha|`:** **`<= 0.0941%` (`med = 0.0056%`) on `12/12` designs** (`12/12` PASS `<= 2.0%`).
  - **`|Delta f_r|`:** **`0.0001%–0.9810%` (`med = 0.390%`) on `12/12` designs**.
  - **`|Delta chi|`:** **`0.0001%–0.0021%` across `NCap 00..04`** (`5/6` NCap designs); `CLT` designs at `z_slabs = 8` had decoupled sub-kHz `|chi| < 25 kHz` denominator floor (`m1_m3_dchi_pct` sentinel `200%` where M3's segmented sweep step exceeded ultra-narrow `< 1 Hz` hanger linewidth).

---

## 8. Richardson GCI, Algebraic Invariants & Thread Determinism (`Gate G3`)

### 8.1 Multi-Level Mesh Convergence (Richardson GCI)
| Archetype (`case_id`) | Solver / Observable | Levels (`r`) | Coarse -> Med -> Fine Values | Observed Order `p_obs` | `GCI_fine` |
|---|---|---|---|---|---|
| `S1.1` (`s1_1_p2_34_r1638`) | `dfs` P2 (`C_xx`) | `4` (`r = 1.35`) | `69.398 -> 69.358 -> 69.320 fF` | `0.50` (asymptotic plateau) | **`0.420%`** |
| `S1.2` (`s1_2_p2_06_r084`) | `dfs` P2 (`C_xx`) | `4` (`r = 1.35`) | `31.505 -> 31.814 -> 31.349 fF` | **`2.00` (exact nominal)** | **`2.254%`** |
| `S1.4` (`s1_4_p2_00_r00516`) | `dfs` P2 (`C_xx`) | `4` (`r = 1.35`) | `43.075 -> 7.917* -> 42.583 fF` | `0.05` (Steiner lock/unlock) | `628.55%`* |
| `S2.1` (`s2_1_smoke_01_r0252`) | `ems` QEP (`f_r`) | `3` (`r = 1.25`) | `6.162711 -> 6.171172 -> 6.173956 GHz` | **`4.00` (`4.98` raw)** | **`0.0391%`** |
| `S2.2` (`s2_2_p2_02_cc696_cp047`) | `ems` QEP (`f_r`) | `3` (`r = 1.25`) | `7.520383 -> 7.526413 -> 7.530002 GHz` | **`2.33`** | **`0.0876%`** |
| `S4.1` (`s4_1_p1_wm1_q1_w1`) | `ems` QEP (`f_r`) | `3` (`r = 1.25`) | `6.027325 -> 6.023037 -> 6.039637 GHz` | **`1.50` (exact nominal)** | **`0.8642%`** |

*\*At `l+1` (`scale_h = 0.74`), `s1_4_p2_00_r00516` unlocked cross-gap Steiner refinement (`7.917 fF`, `+13.9%` vs. Q3D `6.95 fF`) before `300 nm` fillet chords re-exhausted `max_steiner_multiplier = 20` at `l+2` (`42.583 fF`).*

### 8.2 Scoped Algebraic Invariants Across All `418` Runs
- **Electrostatics (`n = 249`):** Symmetry `sym_err` median `8.20e-16` (`max = 4.87e-15`); SPD `min_eig` median `+9.21e-15 > 0` (`100%` positive definite); M-matrix `max_offdiag_pos_ratio = 0.000e0` on **`247/249` (`99.2%`)** runs (`2` `Mode G2` `1/R` extrapolation violations `3.75e-3` on `s1_1_p2_00_r0148`).
- **`ems::cavity` `CircuitRom` (`n = 134`):** Modal equipartition `equipartition_err` median **`3.83e-12`** (`max = 1.32e-10 < 1e-8`, `100%` PASS); eigenvector phase collinearity `phase_noncollinearity_max = 0.000e0` (`100%` PASS); dressed state `min_overlap` median `0.9995` (`min = 0.614` at avoided crossing).
- **`ems::driven` Loewner/PROM (`n = 23`):** Reduced residual `max_prom_res` median **`1.17e-6`** (`max = 8.49e-6 < 1e-3`, `100%` PASS); reciprocity `reciprocity_err` median **`1.49e-13`** (`max = 3.09e-9`); lossless unitarity `unitarity_err` median **`1.54e-12`** (`max = 1.66e-5`).

### 8.3 `NYX_THREADS=1` vs. `NYX_THREADS=8` Determinism (`8` Paired Cases)
| Suite | Case ID | Solver Mode | Max Rel Drift `|y_8 - y_1| / |y_1|` | Worst Key | `T = 1` Wall | `T = 8` Wall | Speedup | `< 1e-13` |
|---|---|---|---|---|---|---|---|---|
| `S1.1` | `s1_1_p2_34_r1638` | `dfs` `far_field` | **`0.000e+00` (Bitwise Exact)** | `top_to_top` | `164.07 s` | `44.76 s` | **`3.67x`** | **PASS** |
| `S1.2` | `s1_2_p2_00_r203` | `dfs` `far_field` | **`0.000e+00` (Bitwise Exact)** | `top_to_top` | `66.04 s` | `20.75 s` | **`3.18x`** | **PASS** |
| `S1.3` | `s1_3_p2_00_r200` | `dfs` `far_field` | **`0.000e+00` (Bitwise Exact)** | `top_to_top` | `53.97 s` | `17.37 s` | **`3.11x`** | **PASS** |
| `S1.4` | `s1_4_p2_00_r00516` | `dfs` `far_field` | **`0.000e+00` (Bitwise Exact)** | `top_to_top` | `13.65 s` | `4.43 s` | **`3.08x`** | **PASS** |
| `S2.1` | `s2_1_smoke_01_r0252` | `ems` `cavity_qep` | `1.152e-07` | `q_c` | `15.96 s` | `18.17 s` | `0.88x` | **FAIL** (Bug 5) |
| `S2.2` | `s2_2_p2_02_cc696_cp047` | `ems` `cavity_qep` | `3.594e-10` | `q_c` | `15.62 s` | `17.27 s` | `0.90x` | **FAIL** (Bug 5) |
| `S3.1` | `s3_1_p2_ncap_02_r203` | `ems` `driven_cap` | `5.751e-09` | `c_1_1_f` | `23.22 s` | `22.79 s` | `1.02x` | **FAIL** (Bug 5) |
| `S4.1` | `s4_1_p1_wm1_q1_w2` | `ems` `w2` | `9.015e-09` | `kappa_hz` | `53.89 s` | `42.31 s` | `1.27x` | **FAIL** (Bug 5) |

---

## 9. PR CI Smoke Tests, Intentional Core-Bug Regressions & Repository Gates

Executed sequentially in `.jj-workspaces/verifier` (`CARGO_TARGET_DIR = bench-verifier`, `--release`, `--test-threads=1`). All **7 intended-pass physics smoke tests pass at their strict original tolerances** (`27.45 s` and `38.49 s`, both well within the `< 90 s` PR CI budget), while all **5 uncovered core engine bugs (`bug1..bug5`) fail deterministically as explicit regression tests** on unmodified `HEAD` (`9480093a`):

### 9.1 `cargo test --release -p nyx-analysis --test squadds_smoke -- --test-threads=1` — `3 passed; 1 failed` (`27.45 s < 90 s`)
| Test | Status | Observed Verification / Core-Bug Trigger |
|---|---|---|
| `transmon_cross_3net_mode_g2_far_field_and_i_edge_thickness_smoke` | **PASS** | Pinned `inner_vac_box` + `sub_box` (`C0.2` setup): Mode `G2` `(0.6L, 1.2L)` `1/R` vs. single-box `AF-3 far_field(0.6L)` agree within `< 2.5%` self / `< 3.0%` mutual; `FR-11` (`t = 200 nm`) `I_edge` strictly increases diagonal `C_ii` and lowers `E_C` while preserving symmetry & M-matrix signs |
| `dissipator_dc_squid_craig_bampton_vs_rom_port_flux_sweep_smoke` | **PASS** | 4-mode `S4.3 Dissipator` ROM (`qubit`, `dissipator_jj`, `readout`, `filter`): 25-point native `squid_flux_sweep` (`10.32 -> 2.93 GHz` across both cavities in `< 0.25 s`) matches point-by-point `ej_hz_override` `craig_bampton` eigenvalues to `< 1e-12` |
| `muninn_3d_tls_participation_nine_param_oxide_and_purcell_smoke` | **PASS** | 3D Xmon + readout coupler `ems::cavity` -> 9-parameter MIT-LL `MA/MS/SA` `tls_participation` (`d = 200 nm`, `edge_cutoff = 5.0 um`, positive `p_ma, p_ms, p_sa` & edge-tube participations) + multi-mode Purcell `T_1` |
| `bug1_fr11_metal_thickness_on_finite_die_fails` | **FAIL (Intentional — Core Bug 1)** | `crates/nyx-analysis/src/tls.rs:789-797`: `Invalid("metal_thickness_m found no interior zero-thickness sheet perimeter edges on the 3D mesh")` when `sub_box` is confined to a finite silicon die inside outer `vac_box` |

### 9.2 `cargo test --release -p nyx-engine --test squadds_engine_smoke -- --test-threads=1` — `4 passed; 4 failed` (`38.49 s < 90 s`)
| Test | Status | Observed Verification / Core-Bug Trigger |
|---|---|---|
| `clt_quarter_wave_qep_vs_driven_taubin_fit_smoke` | **PASS** | Coplanar `ground_plane` + $\lambda/4$ `CLT` hanger (`C0.3` setup): real `CircuitRom` QEP `(f_r, kappa)` vs. 2-port `ems::driven` Loewner/PROM + Taubin `s21_fit` on identical `SimplexMesh` agree within `<= 0.1%` on `f_r` and `<= 2.0%` on `kappa` (`max_residual < 2e-4`, `reciprocity < 1e-10`, `unitarity < 1e-8`) |
| `unified_qubit_cavity_m1_m2_m3_equivalence_smoke` | **PASS** | Identical 3D mesh across closed `craig_bampton` (`M1`), passive-port QEP + `craig_bampton` (`M2`), and 3-port `ems::driven` Schur BBQ + exact two-pole Kubo (`M3`, TRK sum rules `1.0000000000 +/- 1e-10`): `|M1 - M2|` `f_01 <= 0.5%`, `alpha <= 2.0%`, `chi <= 5.0%`; `|M3 - M2|` `f_01 <= 2.0%`, `alpha <= 5.0%`, `f_r <= 0.2%`, `chi <= 15.0%` |
| `airbridge_hole_post_and_deck_sheet_facet_preservation_smoke` | **PASS** | `S3.4` 3D `Op::Box` (`kind: "Hole"`) airbridge landing posts + elevated `Op::Sheet` bridge decks at `z = 4 um` preserve all physical deck & lateral post facets in `SimplexMesh` and pull CPW resonance downward (`f(N_b=2) < f(N_b=0)`) |
| `nyx_threads_1_vs_8_determinism_smoke` | **PASS** | Fresh-process `NYX_THREADS=1` vs. `NYX_THREADS=8` on 3D P2 `dfs::electrostatic` (`far_field` + `FR-11` `200 nm` thickness correction): relative drift `< 1e-13` (bitwise deterministic AMG) |
| `bug2_cavity_refined_lumped_port_stalls` | **FAIL (Intentional — Core Bug 2)** | `crates/nyx-ems/src/cavity.rs:646-694`: refining `port_in:edge` on a `rank1` `lumped_port` sheet fails with `Solver("no eigenfrequencies converged near the shift")` |
| `bug3_sheet3d_z_coarsen_cap_collapse_and_distortion` | **FAIL (Intentional — Core Bug 3)** | `crates/femcore/src/meshing/sheet3d.rs`: multi-layer finite-die stack (`P = 1.2 L`) with `z_coarsen: true` fails with `"boundary number(s) [9] do not exist in this mesh (it has [1, 2, 3, 4, 5, 6, 7, 8, 10, ..., 26])"` |
| `bug4_sheet3d_narrow_gap_steiner_locking` | **FAIL (Intentional — Core Bug 4)** | `crates/femcore/src/meshing/mesher.rs` / `sheet3d.rs`: `2 um` filleted interdigital capacitor (`314 nm` tip/pocket chords) locks volumetric P2 FEM vs. `planar_bem` (`C_aa +134.3%`: `118.08 fF` vs. `50.39 fF`; `C_ab +183.3%`: `92.62 fF` vs. `32.69 fF`) |
| `bug5_ems_thread_determinism_drift` | **FAIL (Intentional — Core Bug 5)** | `ADR 11` thread-determinism violation in 3D `ems::cavity` parallel `faer` sparse factorization: `NYX_THREADS=1` (`-4.3748428283937454e7 Hz`) vs. `NYX_THREADS=8` (`-4.3748428284030914e7 Hz`) on `alpha_Hz`, `rel = 2.136e-12 >= 1e-13` |

### 9.3 Full Repository Verification Gate Summary (`.jj-workspaces/verifier`)
| Repository Gate Command | Observed Result | Wall Time |
|---|---|---|
| `cargo fmt -p geom -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-analysis -p nyx-engine -- --check` | **PASS** (`EXIT=0`) | `1.1 s` |
| `cargo clippy --release --workspace --all-targets --no-deps -- -D warnings` | **PASS** (`0` warnings across workspace, `EXIT=0`) | `22.20 s` |
| `cargo test --release -p nyx-analysis --test squadds_smoke -- --test-threads=1` | **`3 passed; 1 failed`** (only intentional `bug1` fails) | **`27.45 s`** (`< 90 s`) |
| `cargo test --release -p nyx-engine --test squadds_engine_smoke -- --test-threads=1` | **`4 passed; 4 failed`** (only intentional `bug2..bug5` fail) | **`38.49 s`** (`< 90 s`) |
| `cargo test --release -p femcore -p geom -p nyx-store -p nyx-mcp -p nyx-physics -p nyx-dfs -p nyx-sms -p nyx-ems -p nyx-engine --test contracts --test layering` | **PASS** (`72 passed; 0 failed` across 10 test binaries, `EXIT=0`) | `9.21 s` |
| `cargo test --release -p nyx-engine --test api_docs` | **PASS** (`6 passed; 0 failed`, `EXIT=0`) | `0.67 s` |
| `bash docs/architecture/check_code_to_adr.sh` | **PASS** (`1360 files checked, 0 without a row, 322 patterns`, `EXIT=0`) | `0.4 s` |
