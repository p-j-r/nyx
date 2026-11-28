+++
title = "nyx — Finite-Element Simulation for Photonics, Phononics & Quantum Circuits"
description = "Pure-Rust finite-element workbench for waveguide modes, elastic anchor loss, 3D cavities, electrostatics, heat transfer, and superconducting quantum circuit quantisation."

[extra]
hero_kicker = "RUST 1.99 WORKBENCH · v0.1.0 · AGPL-3.0-OR-LATER"
hero_title = "Finite-element physics from first principles, in pure Rust."
hero_lead = "Waveguide modes, elastic anchor loss, 3D curved-element cavities, Prager–Synge dual electrostatic bounds, and hierarchical superconducting circuit quantisation — driven identically from a desktop window, Monty scripts, headless batch CLI, or an MCP server."
bibtex = """@software{robin_nyx_2026,
  author  = {Robin, Paul Joseph},
  title   = {{nyx}: Finite-Element Simulation for Photonics, Phononics, Piezoelectrics, and Superconducting Quantum Circuits},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/p-j-r/nyx}
}"""

[[extra.slides]]
src = "figures/soi_wire_modes.png"
title = "SOI Rib Waveguide · Vector P2 EM (`15_soi_wire`)"
caption = "986 cells, 515 points, 2,015 DOFs, solved in 0.08 s; neff 1 = 2.671811 (ng = 4.1969, 94.5% in Si); neff 2 = 2.333849 (ng = 4.6294, 85.1% in Si); sparse factorisation (6,819 unknowns, 150,676 non-zeros): 0.04 s."
alt = "Vector P2 optical mode gallery of a 500 by 300 nm silicon-on-insulator waveguide showing TE0 and TM0 electric field distributions at 1550 nm."
study = "15_soi_wire · ems::modes"

[[extra.slides]]
src = "figures/wgm_cyl_pml.png"
title = "Cylindrical WGM Whispering Gallery · Axisymmetric PML (`18_wgm_cylindrical_pml`)"
caption = "5,225 DOFs, solved in 0.38 s; residual 8.716e-14; mode 1 = 2.40057, mode 2 = 2.41462 + 6.65e-9i, mode 3 = 2.38500 + 3.85e-8i."
alt = "Four-mode gallery of a high-index dielectric cylinder surrounded by an outward cylindrical PML annulus."
study = "18_wgm_cylindrical_pml · ems::modes"

[[extra.slides]]
src = "figures/anchor_loss_modes.png"
title = "MEMS Anchor Loss · 6-Mode Elastic PML (`11_anchor_loss`)"
caption = "6 coupled elastic-PML eigenmodes from 183.703 Hz (Q = 14,444.6) to 3.096 kHz (Q = 119.3)."
alt = "Six-panel deformed elastic eigenfrequency gallery of a beam clamped to a Cartesian PML substrate."
study = "11_anchor_loss · sms::eigenfrequency"

[[extra.slides]]
src = "figures/cavity_3d_cutaway.png"
title = "3D Spherical-PML Cavity · Curved P2 Cutaway (`20_cavity_3d_spherical_pml`)"
caption = "4 damped 3D cavity modes via curved Nedelec elements, 1.642 GHz (Q = 0.5) to 2.824 GHz (Q = 0.6)."
alt = "3D orthographic cutaway rendering of four open spherical-PML electromagnetic cavity modes."
study = "20_cavity_3d_spherical_pml · ems::cavity"

[[extra.slides]]
src = "figures/cavity_3d_slice.png"
title = "3D Cavity Un-Averaged P2 Field Slice (`20_cavity_3d_spherical_pml`)"
caption = "Interactive slice3d plane=\"z\" cutaway preserving inter-element discontinuities."
alt = "3D planar z-slice rendering of un-averaged P2 electric fields inside the open spherical-PML resonator."
study = "20_cavity_3d_spherical_pml · ems::cavity"
+++

`nyx` is a desktop and headless finite-element workbench written in Rust `1.99.1` (`AGPL-3.0-only OR LicenseRef-nyx-Commercial`). All numerics — geometry evaluation, constrained Delaunay and 3D isosurface-stuffing meshing, finite-element assembly, sparse direct factorisations (`faer`), smoothed-aggregation algebraic multigrid, complex Krylov eigensolvers, and pure-Rust HDF5 storage (`hdf5-pure` + `zlib-rs`) — execute without any external C solver or Python runtime dependency.

## One engine, four front-ends

Every action in `nyx` maps onto a single versioned catalog of **114 JSON-Schema commands** (`nyx api list`), ensuring identical numerical results across:

1. **The desktop window** (`nyx`, built on `eframe`/`egui` with `wgpu` rendering and live field probing).
2. **Typed Monty scripts** (`*.nys`, Python syntax executed in an isolated `nyx-script-worker` sandbox).
3. **Headless CLI & shared queue** (`nyx --batch`, `nyx check`, `nyx run`, `nyx queue submit`).
4. **AI host integration** (`nyx-mcp` and `nyx-connect` for Antigravity, Claude Code, Claude Desktop, Cursor, VS Code, and Windsurf).

## Core physics & verification pillars

- **Electromagnetic waves (`nyx-ems`):** Full-vector Nédélec(2)+P2 cross-section waveguide modes (`exp(-j beta y)` in the wafer frame), axisymmetric rings, Floquet–Bloch periodic lattices, driven multi-port Loewner/PROM frequency sweeps, surface-impedance sheets, and 3D cavity eigenmodes with gradient projection and curved tetrahedral geometry.
- **Solid mechanics (`nyx-sms`):** Linear anisotropic elasticity, 2.5D elastic waveguide modes (`exp(-j q y)`), axisymmetric eigenfrequencies with `r = 0` regularity enforcement, Kelvin–Voigt/Akhiezer viscoelastic damping, low-reflecting boundaries, and Cartesian/cylindrical/spherical perfectly matched layers (`Q = Re f / (2 |Im f|)`).
- **Diffusion & static fields (`nyx-dfs`):** P1/P2 electrostatics (`2D`/`3D` plus 3-layer planar BEM `planar_bem`), stationary and transient heat transfer (exact modal/matrix exponential for small systems, BDF2 with AMG above 2,500 free unknowns), magnetostatics, and AC surface loss.
- **Verification & cQED analysis (`nyx-analysis` & `femcore::vvuq`):** Equilibrated `H(div)` Raviart–Thomas Prager–Synge dual capacitance bounds (`[C_low, C_up]`), Stokes/Poincaré optical polarisation classification, Craig–Bampton hierarchical circuit quantisation, Schur-complement Black-Box Quantisation, Taubin `S_21` circle fits, and ASME V&V 20 validation suites.
