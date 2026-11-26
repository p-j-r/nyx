# Nyx Phase 4 Validation Benchmark Datasets

Versioned, cited external experimental and cross-code validation datasets used by `crates/nyx-engine/tests/validation_*.rs`.

## ASME V&V 20 Acceptance Rule

Every benchmark entry in this directory defines:
- `reference_value` ($D$): Experimental measurement, published cross-code oracle (Ansys Q3D/HFSS, COMSOL Multiphysics, AWS Palace, MPB/MEEP/Lumerical MODE, `scqubits`), or standardized analytical benchmark target (NAFEMS T1–T4).
- `u_input`: Documented $1\sigma$ aleatoric uncertainty from material constants, fabrication tolerances, or geometric etch bias.
- `u_d`: Documented $1\sigma$ experimental or reference-solver benchmark uncertainty.
- `tol_rel`, `tol_abs`: Explicit domain acceptance ceiling ($\mathrm{tol} = \mathrm{tol}_{\mathrm{rel}} |D| + \mathrm{tol}_{\mathrm{abs}}$).
- `coverage_k`: ASME V&V 20 coverage factor ($k = 2.0$ unless tighter $k = 1.0$ is set).

Validation suites compute $u_{\text{num}}$ directly from mesh refinement (2-grid/3-grid Roache GCI) or dual bounds, form:
$$E = S - D, \qquad u_{\text{aleatoric}} = \sqrt{u_{\text{input}}^2 + u_D^2}, \qquad u_{\text{val}} = \sqrt{u_{\text{num}}^2 + u_{\text{input}}^2 + u_D^2}$$
and enforce the **explicit `AND` rule** with strict epistemic/aleatoric separation:
1. $\lvert E \rvert \le \mathrm{tol}_{\mathrm{rel}} \lvert D \rvert + \mathrm{tol}_{\mathrm{abs}}$ **AND**
2. $\lvert E \rvert \le k \, u_{\text{val}}$ (or for deterministic Prager–Synge intervals $[S_{\text{low}}, S_{\text{up}}]$, $\lvert S_{\text{mid}} - D \rvert \le \delta_{\text{ep}} + k \, u_{\text{aleatoric}}$ without folding $\delta_{\text{ep}}$ into Gaussian variance).

## Dataset Inventory

| File | Domain | Primary Citations |
|---|---|---|
| `cqed_park2026.json` | cQED Experimental | Park et al., *npj Quantum Inf.* (2026) 100 nm Nb CPW & 35 nm disordered Nb `Q1+R1` |
| `cqed_palace_q3d_squadds.json` | cQED Cross-Code | AWS Palace (`spheres`/`cylinder`), Ansys Q3D / SQuADDS planar transmon & resonator circuits |
| `cqed_scqubits_tls.json` | cQED Spectra & TLS | Frozen Python/spectral `scqubits` Fluxonium+resonator & Transmon spectra; Calusine (2018) / Woods (2019) / Place (2021) cryogenic $Q_{\text{TLS}}$ |
| `photonics_benchmarks.json` | Photonics & Bloch | 450×230 nm Si-on-sapphire, Chrostowski (2015) / Lumerical/MPB SOI wires, Bragg stack, 3D SWG Bloch & leaky grating |
| `solid_mems_piezo_benchmarks.json` | MEMS, Anchor Loss, Piezo | Bindel & Govindjee (2005) / COMSOL PML disk anchor loss, Akhiezer viscoelastic $Q_m$, Tsubouchi (1985) AlN BAW FBAR |
| `thermal_nafems_benchmarks.json` | Thermal & Thermomechanical | NAFEMS T1 (Hot-spot), T2 (Radiation), T3 (1D Transient), T4 (2D Convection), Timoshenko (1925) Si–sapphire bimorph |
