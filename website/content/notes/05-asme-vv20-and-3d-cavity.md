+++
title = "ASME V&V 20 Epistemic Separation, 3D Curved Cavities & Thread Determinism"
date = 2026-11-26
description = "3D Nédélec curl-curl cavities with gradient projection and un-averaged CellFields slicing, ASME V&V 20 validation gates, and honest disclosure of parallel faer thread drift."

[extra]
note_number = "05"
crates = "nyx-ems::cavity · femcore::vvuq · nyx-analysis::view3d"
reading_time = "9 min"
+++

Three-dimensional electromagnetic cavity and resonator solvers must solve the curl-curl eigenproblem

$$
\nabla \times \left(\mu^{-1} \nabla \times \mathbf{E}\right) = \left(\frac{2\pi f}{c\_0}\right)^2 \varepsilon \mathbf{E}
$$

without letting its infinite-dimensional nullspace $\nabla H^1(\Omega)$ pollute physical modes near $f > 0$. In `crates/nyx-ems/src/cavity.rs`, `nyx` discretises $\mathbf{E}$ using tangentially continuous Nédélec elements of the first kind on tetrahedra (`order = 1` Whitney edge functions; `order = 2` with `20` DOFs per tetrahedron: `12` edge/gradient functions and `8` face functions, integrating curved `CurvedTet` cells with a degree-6 quadrature rule) and projects out the discrete gradient nullspace inside `projected_shift_invert_complex`.

## 3D Open Spherical-PML Cavity & Un-Averaged `CellFields`

In `examples/20_cavity_3d_spherical_pml.nyx.json`, a `25 mm` radius air sphere is enclosed by a `10 mm` 5-layer spherical PML shell (`r_in = 0.025 m`, `kappa_max = 2.0`, `strength = 16.480`, `order = 1`). Solving `3,622` Nédélec DOFs around `f_search = 6.0 GHz` completes in `0.68 s` (`135` restarts, max algebraic residual `1.054e-8`), recovering the triply-degenerate fundamental leaky resonance at `1.642024e9 + 1.700781e9i Hz` (`Q = 0.483`) and a fourth mode at `2.824384e9 + 2.557435e9i Hz` (`Q = 0.552`).

{{ <plate src="figures/cavity_3d_cutaway.png" number="PLATE 04" alt="3D orthographic cutaway rendering of four modes in an open spherical PML cavity" caption="cutaway3d rendering of 20_cavity_3d_spherical_pml (3,622 DOFs, 0.68 s): triply-degenerate fundamental leaky resonances (modes 1–3, f = 1.642 + 1.701i GHz) and mode 4 (2.824 + 2.557i GHz)." meta="20_cavity_3d_spherical_pml · cutaway3d · 3,622 DOFs" /> }}

{{ <plate src="figures/cavity_3d_slice.png" number="PLATE 05" alt="Planar z-slice rendering of un-averaged CellFields inside the 3D spherical PML cavity" caption="slice3d plane=\"z\" cross-section of the same 4 modes evaluated directly on un-averaged per-cell CellFields, preserving sharp discontinuities across region boundaries." meta="20_cavity_3d_spherical_pml · slice3d plane=\"z\"" /> }}

Notice why `figure` (`crates/nyx-analysis/src/view3d.rs`) stores and evaluates `CellFields` per cell node (`cell_fields[c * nloc + a]`) rather than averaging fields across vertices before plotting: across dielectric or PML interfaces, normal components $E\_\perp$ are genuinely discontinuous, and smoothing across shared vertices smears boundary fields. Both `cutaway3d` and `slice3d` evaluate scalar expressions (`|E|`, `re(Ex)`, `abs2(E)`) directly inside each tetrahedron (subdivided into 8 sub-tetrahedra for P2 cells).

{% <solver_log cmd="nyx run render_site_figures.nys  # 20_cavity_3d_spherical_pml"> %}
ems: cavity, order=1, dim=3, dofs=3622, seconds=0.6795
  solver_kind: projected_shift_invert_complex, converged: true, restarts: 135
  max_algebraic_residual: 1.0544204088520082e-08
  mode 1: 1.642024e9 + 1.701e9i Hz  (Q = 4.827e-1, layer_frac = 0.9067, res = 5.610e-09)
  mode 2: 1.642024e9 + 1.701e9i Hz  (Q = 4.827e-1, layer_frac = 0.9067, res = 1.054e-08)
  mode 3: 1.642024e9 + 1.701e9i Hz  (Q = 4.827e-1, layer_frac = 0.9067, res = 8.801e-09)
  mode 4: 2.824384e9 + 2.557e9i Hz  (Q = 5.522e-1, layer_frac = 0.7033, res = 2.225e-11)
{% </solver_log> %}

## ASME V&V 20 Validation & Epistemic Separation

`crates/femcore/src/vvuq.rs` and `crates/nyx-engine/tests/validation_*.rs` enforce ASME V&V 20 validation across electrostatics, microwave cavities, photonics, MEMS anchor loss, and NAFEMS thermal benchmarks (`T1–T4` plus Timoshenko's coupled Si/Al bimorph).{% <sidenote id="sn-mms" label="1"> %}To verify PDE assembly at the code level without Cargo feature thrashing, `femcore::vvuq` uses a zero-sized monomorphised `SourceField<const D: usize>` trait (`UniformSource` in production; `FnSource` in MMS integration tests) sharing 100% of production quadrature and Piola mappings.{% </sidenote> %}

For every target $D$, the comparison error $E = S - D$ must satisfy both the domain tolerance ceiling and the expanded validation uncertainty simultaneously (`AND` rule):

$$
|E| \le \mathrm{tol}\_{\mathrm{rel}}|D| + \mathrm{tol}\_{\mathrm{abs}} \qquad \text{AND} \qquad |E| \le k \sqrt{u\_{\mathrm{num}}^2 + u\_{\mathrm{input}}^2 + u\_D^2}
$$

where $u\_{\mathrm{num}}$ is computed from 2-grid or 3-grid Roache GCI refinement, while deterministic Prager–Synge intervals $[S\_{\mathrm{low}}, S\_{\mathrm{up}}]$ keep their epistemic half-width $\delta\_{\mathrm{ep}} = \frac{1}{2}(S\_{\mathrm{up}} - S\_{\mathrm{low}})$ strictly outside the Gaussian sum: $|S\_{\mathrm{mid}} - D| \le \delta\_{\mathrm{ep}} + k \sqrt{u\_{\mathrm{input}}^2 + u\_D^2}$.

Representative validated baselines in `benchmarks/validation/`:

- **3D Lossy Copper Cavity (`40 mm × 20 mm × 30 mm`, $\sigma = 5.8\times 10^7\,\text{S/m}$):** `TE101` eigenfrequency $f\_0 = 6.245676\,\text{GHz}$, Pozar unloaded surface-loss quality factor $Q\_0 = 10,806.04$, and 2-port `CircuitRom` QEP vs. `ems::driven` Taubin $S\_{21}$ parity.
- **NAFEMS Thermal & Thermomechanical (`T1–T4` + Bimorph):** `T1` hot-spot membrane stress (`50.50 MPa` at Point D), `T2` radiating bar (`927.008 K`), `T3` BDF2 transient conduction (`36.603 C` at `t = 32 s`), `T4` steady 2D convection (`18.254 C`), and Timoshenko Si/Al bimorph deflection (`0.5855 mm`).

## Honest Disclosures: Parallel `faer` Thread Drift & Open-PML Mode Selection

Two solver subtleties discovered during `Gate G3` and VVUQ verification merit explicit engineering disclosure:

1. **3D `ems` Parallel `faer` Factorisation Drift (Core Bug 5):** While `dfs` parallel smoothed-aggregation AMG is **bitwise identical (`0.000e+00`)** between `NYX_THREADS=1` and `NYX_THREADS=8` across all `4/4` Suite 1 benchmarks (achieving `3.08x–3.67x` speedups), 3D `ems` parallel `faer` `LDL^T`/`LU` sparse direct factorisation exhibits **`3.59e-10` to `1.15e-7`** relative drift between `1` and `8` threads (`docs/SQUADDS_BENCHMARK_REPORT.md` §8.3), failing the strict `< 1e-13` `ADR 11` determinism gate.
2. **Krylov Mode Filtering in Open PML Studies:** In PML-bounded open geometries, continuous Berenger layer modes can enter the shift-and-invert window alongside physical cavity or waveguide modes (`VVUQ_VERIFICATION_PLAN_EGM_REFINED.md` §P1.2). Asking the Krylov eigensolver for extra modes (`+4` over-ask) and filtering on `layer_fraction` ensures physical resonances are never displaced from the returned mode list.
