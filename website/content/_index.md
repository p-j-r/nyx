+++
title = "nyx — Finite-Element Simulation for Photonics, Phononics & Quantum Circuits"
description = "Pure-Rust finite-element workbench for waveguide modes, elastic anchor loss, 3D cavities, electrostatics, heat transfer, and superconducting quantum circuit quantisation."

[extra]
kicker = "RUST 1.99.1 WORKBENCH · v0.1.0 · AGPL-3.0-OR-LATER"
hero_kicker = "RUST 1.99.1 WORKBENCH · v0.1.0 · AGPL-3.0-OR-LATER"
hero_title = "Finite-element physics from first principles, in pure Rust."
hero_subtitle = "Waveguide modes, elastic anchor loss, 3D curved-element cavities, Prager–Synge dual electrostatic bounds, and hierarchical superconducting circuit quantisation — driven identically from a desktop window, Monty scripts, headless batch CLI, or an MCP server."
hero_lead = "Waveguide modes, elastic anchor loss, 3D curved-element cavities, Prager–Synge dual electrostatic bounds, and hierarchical superconducting circuit quantisation — driven identically from a desktop window, Monty scripts, headless batch CLI, or an MCP server."
bibtex = """@software{robin_nyx_2026,
  author  = {Robin, Paul Joseph},
  title   = {{nyx}: Finite-Element Simulation for Photonics, Phononics, Piezoelectrics, and Superconducting Quantum Circuits},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/p-j-r/nyx}
}"""

[[extra.slides]]
plate_id = "PLATE 01.1"
image = "figures/soi_wire_modes.png"
badge = "PHOTONICS · 15_soi_wire"
title = "SOI Wire Waveguide · Vector Nédélec(2)+P2 (`15_soi_wire`)"
caption = "500 × 300 nm silicon wire at 1550 nm in the wafer frame: TE0 (neff = 2.671811, ng = 4.1969, 94.5% power in Si) and TM0 (neff = 2.333849, ng = 4.6294, 85.1% in Si)."
alt = "Vector P2 optical mode gallery of a 500 by 300 nm silicon-on-insulator waveguide showing TE0 and TM0 electric field distributions at 1550 nm."
command = "nyx --batch examples/15_soi_wire.nyx.json"
metrics = "986 cells · 2,015 DOFs · 0.08 s"
solver_log = [
  "986 cells, 515 points",
  "ems: solved in 0.08 s, 2015 dofs, 2 solution(s)",
  "  neff 1 = 2.671811   (n_g = 4.1969, share in Si = 94.5%)",
  "  neff 2 = 2.333849   (n_g = 4.6294, share in Si = 85.1%)",
  "  factorising sparse matrix (6819 unknowns, 150676 non-zeros): 0.04 s"
]

[[extra.slides]]
plate_id = "PLATE 01.2"
image = "figures/wgm_cyl_pml.png"
badge = "OPEN RADIATION · CYLINDRICAL PML"
title = "Whispering-Gallery Cylinder · Axisymmetric PML (`18_wgm_cylindrical_pml`)"
caption = "5-layer structured quadrilateral cylindrical PML annulus around a 1.2 µm core: max algebraic residual 8.716e-14 across 4 eigenmodes in 0.38 s."
alt = "Four-mode gallery of a high-index dielectric cylinder surrounded by an outward cylindrical PML annulus."
command = "nyx run render_site_figures.nys --workspace ./run  # 18_wgm_cylindrical_pml"
metrics = "5,225 DOFs · 0.38 s · res 8.716e-14"
solver_log = [
  "ems: modes, order=2, dim=2, dofs=5225, seconds=0.3793",
  "  solver_kind: shift_invert_complex, converged: true, restarts: 64",
  "  mode 1: neff = 2.400572 - 1.559e-10i  (core = 0.9900, res = 1.375e-15)",
  "  mode 2: neff = 2.414624 + 6.653e-09i  (core = 0.9899, res = 4.588e-14)",
  "  mode 3: neff = 2.385004 + 3.849e-08i  (core = 0.9932, res = 6.400e-14)"
]

[[extra.slides]]
plate_id = "PLATE 01.3"
image = "figures/anchor_loss_modes.png"
badge = "SOLID MECHANICS · ANCHOR LOSS"
title = "MEMS Cantilever on Elastic PML Substrate (`11_anchor_loss`)"
caption = "Cartesian 5-layer elastic PML separating high-Q flexural modes (183.703 Hz, Q = 14,444.6) from leaky continuous substrate-PML modes (Q = 1.46–1.61)."
alt = "Six-panel deformed elastic eigenfrequency gallery of a beam clamped to a Cartesian PML substrate."
command = "nyx run render_site_figures.nys --workspace ./run  # 11_anchor_loss"
metrics = "4,832 DOFs · 0.086 s · unstable_modes = 0"
solver_log = [
  "sms: eigenfrequency, order=2, dim=2, dofs=4832, seconds=0.0858",
  "  mode 1:  183.703 Hz + 6.359e-3i Hz   (Q = 14444.57, layer_frac = 2.71e-7)",
  "  mode 2: 1132.207 Hz + 0.4566i Hz    (Q =  1239.95, layer_frac = 9.09e-5)",
  "  mode 3: 2286.157 Hz + 710.22i Hz    (Q =     1.61, layer_frac = 0.2715)",
  "  mode 6: 3096.439 Hz + 12.979i Hz    (Q =   119.28, layer_frac = 5.75e-3)"
]

[[extra.slides]]
plate_id = "PLATE 01.4"
image = "figures/cavity_3d_cutaway.png"
badge = "3D CAVITY · SPHERICAL PML"
title = "3D Open Spherical-PML Cavity · Gradient-Projected Cutaway"
caption = "Orthographic cutaway3d view of 4 damped 3D cavity modes via projected_shift_invert_complex around 6.0 GHz."
alt = "3D orthographic cutaway rendering of four open spherical-PML electromagnetic cavity modes."
command = "nyx run render_site_figures.nys --workspace ./run  # 20_cavity_3d_spherical_pml"
metrics = "3,622 DOFs · 0.68 s · 135 restarts"
solver_log = [
  "ems: cavity, order=1, dim=3, dofs=3622, seconds=0.6795",
  "  solver_kind: projected_shift_invert_complex, converged: true",
  "  mode 1..3: 1.642024e9 + 1.701e9i Hz  (Q = 0.483, layer_frac = 0.9067)",
  "  mode 4:    2.824384e9 + 2.557e9i Hz  (Q = 0.552, layer_frac = 0.7033)"
]

[[extra.slides]]
plate_id = "PLATE 01.5"
image = "figures/cavity_3d_slice.png"
badge = "FIELD RENDERING · UN-AVERAGED P2"
title = "3D Cavity Planar Slice · Un-Averaged Per-Cell Fields"
caption = "slice3d plane=\"z\" evaluates |E| directly inside each tetrahedron without smoothing normal discontinuities across interfaces."
alt = "3D planar z-slice rendering of un-averaged P2 electric fields inside the open spherical-PML resonator."
command = "figure(study=\"cavity\", view=\"slice3d\", plane=\"z\", expr=\"|E|\")"
metrics = "view3d · CellFields · sharp interface jumps"
solver_log = [
  "figure: 4 modes rendered to cavity_3d_slice.png",
  "  view = slice3d, plane = z, expr = |E|",
  "  rendering un-averaged per-cell CellFields (interface jumps preserved)"
]
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
