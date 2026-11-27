# nyx

A desktop application for finite-element photonics and phononics (waveguide modes, effective index, polarization, vibration modes, anchor loss,
cavities, heat and electrostatics), written in Rust. The numerics are pure Rust (`femcore`, three physics crates on top of it); no Python is needed to
build or run anything (tests and CI use Python 3.12 with `h5py` and `numpy` to check the results files against a second HDF5 reader). There are three ways to drive one engine: the **window**, **scripts** (`*.nys`, Python syntax) and an **MCP server** for AI hosts.
Every action is an engine command with a JSON Schema (114 of them, one page each in [`docs/api`](docs/api)), so the three front-ends do the same things.

```
crates/geom           geometry kernel: parametric wafer-frame models, named regions, selectors
crates/femcore        the FEM core: meshes, meshing and adaptation, linear algebra (sparse, eigen, multigrid), finite-element basics, materials
crates/nyx-physics    the Physics trait, the problem envelope, descriptors, SolveResult, JSON Schemas
crates/nyx-dfs        diffusion: electrostatics, heat transfer (stationary, transient)
crates/nyx-sms        solid mechanics: static, eigenfrequencies, axisymmetric eigenfrequencies, waveguide modes; fixed, load, low-reflecting, PML, Floquet conditions
crates/nyx-ems        EM waves: waveguide modes, axisymmetric modes, cavity eigenfrequencies; PML, Floquet, magnetic wall
crates/nyx-analysis   analyses of a result: polarization (Stokes, Poincare sphere), figures of modes (SVG, PNG)
crates/nyx-store      workspaces, workflows, experiments and run records (JSON)
crates/nyx-storage    the pure-Rust HDF5 results file engine (read and write, packed layout)
crates/nyx-engine     the Engine: project, commands, jobs and queue, sweeps, results files (HDF5), guide cards
crates/nyx-script     scripts (`*.nys`) run by Monty in a worker process
crates/nyx-mcp        the MCP server (four script tools) ; crates/nyx-connect connects it to AI apps
crates/nyx-install    installer and packager: nyx-setup, nyx-uninstall, nyx-pack (.nyxpack, zip, dmg)
crates/app            the egui application (binary `nyx`)
docs/                 HOW_IT_WORKS.md, docs/architecture (ADRs, C4 model, contracts), docs/api (command reference), docs/ports (provenance of ported code)
```

## Licence

**nyx is dual licensed: you may use it under the GNU Affero General Public License, version 3 (see [LICENSE](LICENSE)), or under a commercial licence, whichever you choose.** (`SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial`)

Copyright (C) 2026 p-j-r. A commercial licence is available from the copyright holder for use that the AGPL does not fit. Parts that carry their own licence keep it and are not relicensed: `crates/geom` and the ported solver and mesher modules of `femcore` (BSD-2-Clause, see `docs/ports`), and the third-party crates in `Cargo.lock`. [NOTICE](NOTICE) lists them.

## Citing nyx

If you use `nyx` in a paper, thesis, or technical report, please cite the repository (also available via GitHub's **Cite this repository** button from [`CITATION.cff`](CITATION.cff)):

> Simulations were performed using `nyx` (v0.1.0, https://github.com/p-j-r/nyx).

```bibtex
@software{robin_nyx_2026,
  author  = {Robin, Paul Joseph},
  title   = {{nyx}: Finite-Element Simulation for Photonics, Phononics, Piezoelectrics, and Superconducting Quantum Circuits},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/p-j-r/nyx}
}
```

## Build and run

```
cargo run --release -p nyx                          # the application
cargo test --release --workspace                    # all tests (a few minutes); `-- --ignored` adds the slow ones
cargo clippy --release --workspace --all-targets --no-deps -- -D warnings
bash tools/vvuq_ci.sh --tier 1                      # local VVUQ fast gate (PowerShell mirror: .\tools\vvuq_ci.ps1 -Tier 1)
bash tools/package_dist.sh                          # release programs in dist/<name>/ with LICENSE, NOTICE, nyx-setup, and zip
```

* Rust 1.99.1 (pinned in `rust-toolchain.toml`). Always build with `--release`: the solvers are unusably slow otherwise. The first build takes several
  minutes (egui, wgpu, faer).
* Results files (`results.h5`, cargo feature `hdf5`, on by default) use 100% pure-Rust `hdf5-pure` + `zlib-rs` (ADR 107): **no external `HDF5_DIR`,
  C library, or `hdf5.dll` is required** to build or run. 
* A nyx install is `nyx`, `nyx-worker` (solves in a child process) and `nyx-script-worker` (runs scripts) side by side, plus `nyx-mcp`, `nyx-connect`
  (AI hosts), and `nyx-uninstall`. `cargo build --release -p nyx-script` before `cargo test -p nyx-mcp`.
* Cross-OS release bundles (`.zip`, self-extracting `nyx-setup`, `.nyxpack`, and macOS `.dmg`) are built by `nyx-pack` locally or via
  [`.github/workflows/release.yml`](.github/workflows/release.yml) on Windows, Linux, and macOS (`arm64` + `x86_64`).
* Several checkouts or workspaces on one machine can share a build directory by setting `CARGO_TARGET_DIR`; `CLAUDE.md` has the conventions and the checks.

[`docs/HOW_IT_WORKS.md`](docs/HOW_IT_WORKS.md) explains in plain language what the program does from mesh to picture and what every line of the solver log
means; every source file of `femcore` starts with a short description of what it does.

## Using nyx

**The window** (`nyx`, or `nyx FILE` to open one). File > New workflow (Ctrl+N) and Open (Ctrl+O), Edit > Search or run anything (Ctrl+K), a Dashboard,
Settings, Run > Compute (F5) and Stop (Esc), Help > About nyx and licences. The Study, Physics and Materials tabs are drawn from the descriptors of the
physics. A project holds the geometry (one `geom` model, regions and boundaries named, selectors such as `outer & plane(x = 0)`), materials per domain
(a library plus the project's own), one problem per physics, and mesh settings. A 3D project is the 2D section swept along y. `nyx --size 1400x860 --script
"..."` drives the window from the command line (for screenshots).

**The command line.**

```
nyx --batch model.nyx.json [--physics ems,sms]   # solve without a window; prints the summary and the stage log
nyx check FILE.nys ; nyx run FILE.nys [--workspace DIR]   # scripts (a script may only use the one folder it is given)
nyx api list ; nyx api show COMMAND [--json]     # the command reference, same text as docs/api
nyx queue submit WORKSPACE PROJECT.nyx.json [PHYSICS...] [--priority 1|2|3] [--hold]   # the queue shared by several programs
```

**Scripts and AI hosts.** A `*.nys` script calls the engine commands (`nyx.project_open(path=...)`, `nyx.solve_run(...)`, ...); samples are in
[`docs/scripting-samples`](docs/scripting-samples). The MCP server has four tools (`script_run`, `script_wait`, `script_check`, `api`) over one script
session; `nyx-connect list | connect <app|all> | disconnect | test | instructions | snippet` (and the Connect menu) wire it into Claude Code,
Claude Desktop, Cursor, VS Code, Windsurf and Antigravity. What costs time with an AI host is in [`docs/HOST_SETUP.md`](docs/HOST_SETUP.md).

**Jobs, sweeps and the queue.** Meshing and solving run as jobs in worker processes (cancellable, with a memory estimate); `sweep_run` records one run per
point of a parameter sweep, and a warm start can carry the modes of one point to the next; the queue holds work for several programs on one machine.

**Examples** (`examples/`, `01_capacitor` to `20_cavity_3d_spherical_pml`): one or more per physics, among them a WR-90 waveguide, a PML slab, a box cavity,
a silicon-on-insulator wire (`15_soi_wire`) and its hybrid version, silicon on sapphire, anchor loss of a disk, and PML examples in cylindrical, axisymmetric and 3D spherical form (`18` to `20`). They are found next to the program or in a
parent folder.

## What is solved

| Physics | Studies and features |
|---|---|
| `dfs` | **Electrostatics** 2D/3D (per-domain eps_r, space charge, potential, ground, surface charge, energy, capacitance). **Heat transfer**, stationary and transient, 2D/3D (k, rho, Cp, Q, temperature, flux, convection; exact semi-discrete solution, BDF2 for large models) |
| `sms` | **Static** (plane stress or strain, 3D), **eigenfrequencies**, **axisymmetric eigenfrequencies** (bodies of revolution on the (r, z) half plane, axis conditions, cylindrical absorbing layers), **waveguide modes** (`modes`, 2.5D elastic waves `exp(-j q y)` on the 2D (x, z) cross-section at wavevector `q`); fixed constraints, boundary loads, thermal expansion; low-reflecting boundaries, PML, Floquet-periodic pairs; per mode frequency parts and energy shares per component and region |
| `ems` | **Waveguide modes** (full-vector Nedelec(2)+P2 on the cross-section, complex neff and Q with PML, tilted anisotropic crystals), **axisymmetric modes**, **cavity eigenfrequencies** (Nedelec edge elements of degree 1 and 2 on tetrahedra, optional curved cells), PML, Floquet-periodic pairs (Bloch dispersion), magnetic walls |

Elements are P1/P2 Lagrange on triangles and tetrahedra (Nedelec for the waves). Stiffness, thermal expansion, conductivity, permittivity and the
photoelastic tensor are full tensors with Euler-angle orientation (`crates/femcore/tests/anisotropic.rs`). The solvers are checked against
closed forms (waveguides, cavities, beams, plates, bars) and against golden files (`crates/nyx-{dfs,sms,ems}/tests/golden`, `crates/geom/golden`); the physics
crates' tests say what each is checked against. Brillouin scattering and the chips built on it were removed
(`docs/architecture/adr/0024-brillouin-and-chips-removed.md`; the last state with them is the tag `sbs-last`).

**Analyses.** Every optical mode study returns the Stokes parameters and the Poincare-sphere point of each mode, computed from the complex field the solver
found (degree of polarization, ellipticity and orientation, TE and longitudinal shares, a class, and the resolved pair for degenerate modes): the `analyse`
command and the Analysis tab (`nyx-analysis`, ADR 56). The older write-up of the 300 nm SOI use cases is kept in [`docs/history/`](docs/history/README.md); its code references predate the restructure.

**Figures of modes.** The `figure` command writes a gallery or a components figure of the modes of a result as SVG or PNG (ADR 99): each mode is put in a
fixed phase (the in-plane motion largest in Re, sign fixed, peak scaled to 1; the reply reports how far the propagation component is in quadrature), drawn on the
cross-section (the first y-plane of a 3D result) with the outline of the regions and arrows of the in-plane motion; a displacement deforms the section, other
fields do not. The window, scripts and the MCP server share it. An example is [`docs/figures/example_modes_gallery.png`](docs/figures/example_modes_gallery.png).

**Also in the engine** (each has a page in [`docs/api`](docs/api)): a `magnetostatic` study and infinite elements in `dfs`; `driven` studies and piezoelectric coupling (ADR 129) in `sms`; Bloch and `driven` studies, lumped ports, wave ports, surface impedance and lumped Josephson junctions in `ems`; cQED analyses in `nyx-analysis` (capacitance matrix, EPR, Hamiltonian, TLS participation, network, piezo resonator, `uq_cqed`); `convergence_check` (grid-convergence index) and `vvuq_export` for verification and validation (ADR 131 to 133, `benchmarks/validation`, `bash tools/vvuq_ci.sh --tier N`); `import_gds` for GDS layouts; `fetch_material` for remote materials. The installer and packager (`nyx-setup`, `nyx-uninstall`, `nyx-pack`) are described under Build and run.

## Axes: the wafer frame

Projects are described in the axes of the wafer: x and y span the wafer, **z is the surface normal and y is the direction of propagation** of a waveguide
(waves go as `exp(-j beta y)`). A 2D cross-section is the x-z plane, mesh coordinates `(x, z)`. Crystal orientations (Euler Z-X-Z angles), tensors, conditions
and the names of result fields (`ux, uy, uz`, `Ez`, ...) all refer to x, y, z. `femcore::section` is the one place that knows axis numbers; the solid kernels use
the Voigt order (x, z, y). Units are SI (metres) everywhere in the library and the tools.

## Meshing

The application meshes `geom` models: 2D cross-sections (constrained Delaunay with Ruppert refinement, the `spade` crate) and 3D models with a
lattice mesher (`femcore::meshing::mesher3d`: isosurface stuffing after Labelle and Shewchuk, interfaces cut exactly, local size per body or face with
`add_size`, `min_layers` for thin films). The mesh setting `curved` gives the 3D cells on a curved surface quadratic geometry (ADR 82, 83); with it off every
number is the one of the straight mesh. Structured PML cells stay native quadrilaterals, prisms and hexahedra (ADR 128). Other cells that are not tetrahedra (hexahedra, prisms, quadrilaterals given as a `CellMesh`) are cut into simplices with no extra
vertices (`SplitMode::Minimal`: 6 tetrahedra per hexahedron, 3 per prism; the older `Star` split has a vertex at every face and cell centre). Metric-based
adaptation lives in `femcore::adapt`.

## Linear solvers

Small systems use a sparse Cholesky or LU factorisation (`faer`). Large ones use conjugate gradients preconditioned by smoothed-aggregation algebraic
multigrid (own Rust code: rigid-body modes as near-nullspace, Chebyshev smoothing, single-precision matrices inside the cycle, and for quadratic elements the
linear subspace as first coarse space). `Automatic` switches multigrid on at 25,000 free unknowns in 3D and 400,000 in 2D (`femcore::fem::AMG_THRESHOLD*`) and falls back
to the factorisation if it does not converge. The lowest eigenfrequencies use LOBPCG with the multigrid cycle as preconditioner above 80,000 free unknowns when a
probe shows multigrid converging fast. Eigenproblems of the waves use a block shift-and-invert Krylov solver written for this project (it finds repeated
eigenvalues such as degenerate modes). The code of the ported solver library (Krylov methods, block Krylov-Schur, Jacobi-Davidson, LOBPCG and contour
eigensolvers, AMG, AMS, ADS, ILU, Schwarz, time steppers, for real and complex systems) is `femcore::linalg`; its provenance, licences and the plan for swapping
it in are in [`docs/ports/README.md`](docs/ports/README.md).

Measured earlier on a 20-core machine (3D, not re-run for this page): a drilled plate (190,000 unknowns, static, P2) takes 27.6 s with the direct solver and 1.7 s
with multigrid; a block of 681,000 unknowns takes 341 s against 7.1 s. Very thin parts (two element layers) make multigrid slow with quadratic elements, and the
eigenfrequency study then chooses the factorisation.

## Speed, threads and Cancel

* All parallel loops run on one shared thread pool (`femcore::par`, rayon): as many threads as logical processors (at most 32); `NYX_THREADS=N` changes that
  (1 = serial). Sparse factorisations can use a smaller pool (`NYX_FACTOR_THREADS`). The same input gives the same numbers for any thread count.
* A running solve is stopped with Stop (Esc) or `job_cancel`; the solvers check at every iteration and at every multigrid level. Earlier results stay untouched.
* Every solve records a log of what it did (`What the solver did` in the Results tab, `log.txt` in a run folder, printed by `--batch`). `NYX_PROFILE=1` echoes it to
  the console while it runs; `NYX_FRAME_LOG=frames.csv` writes the time of every window frame to that CSV file (the variable only has an effect with a `.csv` name).

## Projects, runs and results files

Project and spec files are JSON (`*.nyx.json`, with `format` 1: a file without it or of another format is refused, nothing is converted). Workspaces
(default `Documents\nyx`) hold workflows, experiments and recorded runs. A recorded run keeps its results in its folder as `results.h5` (HDF5; the project, the
mesh and every field of every solution, compressed, read back bit for bit) and `results.xdmf` (see below), with `log.txt`; no JSON copy of the fields is written.
`export_results` writes the same file anywhere and `open_results` reads it back as a project with its results.

**`results.h5`**, layout 4 (ADR 62, 99): `/mesh` (points, cells, facets, labels), `/project`, and per physics `/results/<physics>/{points,cells}` (the
nodes and the node of every cell node, in nyx's cell node order named in the attribute `cell_order`), `solution_k/fields/<name>` (units as attributes) and
`cell_fields` (the complex fields per cell node: the value at node `cells[c, a]` is `cell_fields[c * nloc + a]`). There is no second copy of the nodes or cells
for a viewer.

**`results.xdmf`** (XML, next to the `.h5`): the light index of that heavy data, for **nyx and its tools only, not for ParaView or VisIt** (the file says so: their
order of the quadratic cell nodes differs from nyx's, stated in `CellNodeOrder`). It names `points`, `cells` and the stored fields as they are, with their
shapes; several solutions become a temporal collection. A test resolves every reference against the HDF5 file.

**Parallel output** (cargo feature `phdf5`) is not available: the feature is a stub with no MPI code, and `export_results` with `parallel=true` returns an error ("parallel HDF5 (MPI-IO) is not available in this build"). Results files are always written serially.

## Limits

* Linear materials only (no contact or nonlinear solids). Absorbing layers (complex coordinate stretching) are available for the wave modes and for solid
  eigenfrequencies; low-reflecting boundaries need a search frequency (they make the problem quadratic) and cannot be combined with periodic conditions yet.
* The 2D figure kinds (`gallery`, `components`) are drawn for mode studies with a vector field (displacement or electric field); sweep curves and the
  axisymmetric half plane are not drawn. The 3D kinds (`slice3d`, `cutaway3d`, `isosurface3d`, `mesh3d`, see `docs/api/figure.md`) draw an expression of a 3D result.
* Decisions, and why each was taken, are in the ADRs (`docs/architecture/adr/`); the plans, status logs and rough-spot lists of earlier work are kept for reference in [`docs/history/`](docs/history/README.md) and are not current.

## Where to read next

[`docs/HOW_IT_WORKS.md`](docs/HOW_IT_WORKS.md) (what the program does), [`docs/architecture/README.md`](docs/architecture/README.md) (ADRs, the C4 model, the
interface contracts), [`docs/api`](docs/api) (every command), [`CLAUDE.md`](CLAUDE.md) (layout, build, conventions and the version-control rules for work in this repository).
