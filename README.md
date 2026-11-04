# nyx

A desktop application for finite-element simulations, written in Rust.
The numerics are pure Rust (`femcore`); no Python is needed to build or run anything. The
solvers are checked against closed forms and self-consistency tests (see "What is solved").

```
nyx/
  crates/femcore   the solver library (meshes, sparse + eigen solvers, physics)
  crates/app       the egui application (binary `nyx`)
  docs/            HOW_IT_WORKS.md, docs/architecture (ADRs, contracts, restructure plan) and docs/ports (provenance of the ported libraries)
```

## Build and run

```
cargo run --release -p nyx                 # start the application
cargo run --release -p nyx -- project.nyx.xml   # ... and open a project or HDF5 results file
cargo test --release --workspace              # about 500 tests (about 200 of them for the ported solvers) against closed forms and self-consistency checks, ~2.5 min
cargo test --release --workspace -- --ignored # the slow tests (the solver library, the Bragg grating comparison ~1.5 min)
cargo clippy --release --workspace --all-targets -- -D warnings   # clean
```

Rust 1.84 or newer (checked with 1.99). The first build takes several minutes (egui, wgpu and
faer are compiled with optimisations).

**docs/HOW_IT_WORKS.md** explains in plain language what the program does from mesh to picture and what every
line of the solver log means; every source file of `femcore` starts with a short description of what it does.

## Speed, threads and Cancel

* All parallel loops run on one shared thread pool (`femcore::par`, rayon). It has as many threads as logical
  processors (at most 32); `NYX_THREADS=N` changes that (1 = serial). Sparse factorisations use a smaller pool
  of 8 threads (`NYX_FACTOR_THREADS`) because they stop scaling there. The same input gives the same numbers
  every time.
* A running solve can be stopped with **Cancel** next to Compute (or Esc). The solvers check at every iteration and
  at every multigrid level; a started matrix factorisation has to finish first. Earlier results stay untouched.
* Every solve records a log of what it did (`What the solver did` in the Results tab, `.log.txt` next to autosaved
  runs, printed by `--batch`). `NYX_PROFILE=1` echoes it to the console while it runs; `NYX_FRAME_LOG=1` prints
  the time of every window frame.
* NaN or infinite input (a material, a coordinate, a load) is refused with a message that names it.
* Measured on a 20-thread machine: the Brillouin example of the earlier README (23,240 cells, 141,000 unknowns) took
  about 43 s (it is gone, see `sbs-last`); a
  3.1-million-unknown 3D solid solves in about 40 s in 13 GB; a million-triangle result rotates at ~45 ms per frame.

## Meshes with hexahedra and prisms

Cells that are not tetrahedra are cut into tetrahedra with no extra vertices (`SplitMode::Minimal`:
6 tetrahedra per hexahedron, 3 per prism, cut along the diagonal from the lowest-numbered vertex of
every face, so neighbours always agree). The older star split (a vertex at every face and cell centre,
24 and 14 tetrahedra) is kept as `SplitMode::Star`. The minimal split has far fewer cells and unknowns
(a quarter of the cells for hexahedra), and the AMG-CG solve needs several times fewer iterations on it.
Meshes of mixed cells are described by `femcore::mesh::cells::CellMesh` (the one mesh model, labels from 1); `SimplexMesh::from_cells_split` cuts them into the triangles or tetrahedra the solvers use.

## Nanophotonics cross-sections

The example "Silicon strip on sapphire at 1550 nm" is a 450 nm x 230 nm silicon strip on sapphire
(n = 3.45 / 1.746 / 1, PML on all sides).
* **Optical mode:** neff = 2.31747 (mesh-converged; 2.31753 at 40 nm, 2.31747 at 15 nm element size in the
  strip); a coarse mesh in the automatic test agrees with a mesh twice as fine to 1e-3.
* **Acoustic modes with an out-of-plane wave number** (Study > Eigenfrequencies > "Out-of-plane wave number"):
  three displacement components on the cross-section, `u(x, y) exp(-i kz z)`, full rotated stiffness of
  silicon and sapphire (Euler 60, 57.6, 0 degrees), PML in the sapphire for leaky modes (complex frequency,
  Q), domains that carry no solid (air) switched off for that physics (`inactive`). Checked against Young's
  rod velocity (3 %), Euler-Bernoulli bending (6 %), a clamped plate shear mode (0.2 %) and the analytic
  complex frequency of a bar with a complex stretch (Q = 5, 0.3 %). The example searches around 9.47 GHz at
  kz = 460 rad/m, the forward-scattering setting.
* **Not yet:** large 3D photonic-crystal and suspended-membrane models need 3D vector
  electromagnetic eigenmodes (Nedelec elements), Floquet periodic conditions and low-reflecting boundaries
  at that scale.

Things found by running this geometry: the shape mesher left slivers where shapes share edges (coordinates a
rounding error apart), which made the electromagnetic problem singular; shapes are now made conforming and
triangulated in normalised coordinates. `serde_json` read floats with 1 ulp error, which showed up as
drifting geometry in saved projects; it now round-trips exactly. The complex mode solver is 4x faster
(parallel orthogonalisation, exact residuals without assembling vectors) and omits the clustered radiation
modes inside absorbing layers that would not converge, instead of running to its limit for them.

## Brillouin scattering (removed)

The Brillouin study (optical and phase-matched elastic modes, photoelastic and moving-boundary gain, rings) and the four ring chips were
removed in the restructure (`docs/architecture/adr/0024-brillouin-and-chips-removed.md`). The last state that has them is the tag `sbs-last`.
The kz solid kernel and the axisymmetric kernels that it was built on stay as ordinary studies.

## Ported solver and mesher libraries

`femcore::solvers` (Krylov methods, block Krylov-Schur / Jacobi-Davidson / LOBPCG / contour eigensolvers, AMG / AMS / ADS / ILU /
Schwarz preconditioners, time steppers, nonlinear solvers, optimisers, for real and complex systems) and `femcore::meshers`
(metric-based 2D mesh adaptation, structured PML meshes with quadrilaterals, prisms and hexahedra, the lossless `.fmsh` file format) were
merged from the `inspired-rust-port` repository (BSD-2-Clause). They have their own tests (`tests/linalg`, `tests/meshers.rs`, the unit
tests of both modules) and are **not yet wired into the physics modules**, which still use the older `eigen`, `lobpcg`, `amg` and `linalg`.
`meshing::pml_box::pml_domains` turns the PML regions into `Pml` domains; `adapt::work` converts between the solver mesh and the adaptation's working mesh.
Provenance, licences, what the authors measured and the plan for swapping the solvers in: `docs/ports/README.md`.

## Polarization analysis and the photonic use cases

Every optical mode solver (2D real and complex, 3D eigenmodes and Floquet cells) returns the **Stokes
parameters and the Poincare-sphere point of each mode**, computed from the complex field the solver found: degree of polarization, ellipticity and
orientation, TE and longitudinal shares, a class (quasi-TE, quasi-TM, linear, elliptical, circular, hybrid), and for degenerate modes the group Stokes vector and the
resolved pair (the eigensolver returns an arbitrary mixture of degenerate modes). The application shows it in the **Analysis** tab (sphere, mode list, table, CSV and SVG
export); the table is also in the solver log and the batch output. Modules `femcore::polarization` and `femcore::analysis`.

Use cases (Examples 16 to 19, modules `femcore::usecases::{soi, bragg}`, studies `soi_modes` and `bragg_gratings` in `crates/nyx-engine/examples`):
the modes of a wire on **300 nm SOI** at 1550 nm (TE0 2.6718, TM0 2.3337, against the infinite slab's 3.049 and 2.616; a sloped wall mixes the nearly degenerate pair of a 300 nm
wire into two states tilted by 45 degrees), and **symmetric and asymmetric Bragg gratings** in that wire as 3D Floquet cells at the band edge (stop band 5.8 nm for
80 nm teeth on both sides, 3.2 nm for one side, coarse mesh). Numbers, checks and limits: `docs/photonic_use_cases.md`.

## Axes: the wafer frame

Projects are described in the axes of the wafer (there is no other frame): x and y span the
wafer, **z is the surface normal and y is the direction of propagation** of a waveguide, so a 2D cross-section is the **x-z plane**
(mesh coordinates `(x, z)`, `y` out of the page; the frame is right-handed, so a wave travelling into the page travels along `+y`, and the out-of-plane wave number `k` of an elastic study is positive for it). Crystal
orientations (Euler Z-X-Z angles of the crystal axes in x, y, z), photoelastic tensors, permittivity, constraint and load components and
the names of the result fields (`ux, uy, uz`, `szz`, `Ez`, ...) all refer to x, y, z; the kernels read these axes as they are (`femcore::section` is the one place with axis numbers; golden files in `crates/femcore/tests/golden_frames` froze the results of the earlier kernels, see `docs/KERNEL_PORT_STATUS.md`). There is no other public frame (a `frame` key in an old project file is ignored). 3D problems are unaffected.
Crystal orientations are written in these axes (the ring chips that once showed it are gone with ADR 24; silicon axes `[p, q, z]`, sapphire axes `[a, c x a, c]`, c tilted 57.6 degrees from z towards the
propagation direction at azimuth psi = 0; the wafer-frame stiffness equals a direct fourth-rank rotation of the crystal tensor in
every entry, test `shear_components_of_a_rotated_stiffness_are_in_the_solver_order`).

**A bug found on the way:** the shear components of every rotated stiffness were stored in the wrong order of the solver's Voigt vector
(xz, xy, yz instead of xy, yz, xz). Isotropic and cubic crystals in the axes of the model hide it; sapphire, orthotropic materials and
crystals rotated about other axes were affected. Fixed, with a regression test against a direct fourth-rank rotation.

## Anisotropy everywhere, including the optics

Stiffness (6x6, rotated), thermal expansion, conductivity, permittivity and the photoelastic tensor (rank 4, rotated) are full tensors in every solver.
New: **optical modes of crystals tilted against the propagation axis.** A permittivity with `xz`/`yz` terms (sapphire with its c axis
in the plane of the wafer's normal and the guide) couples the longitudinal to the transverse field in the electric term, and the 2D problem
becomes quadratic in the propagation constant, `(A0 + j b A1 + b^2 B) x = 0`. It is solved exactly (fixed-point iteration on `b` in the
linear term with Aitken acceleration, residual below 1e-10), also with PML (`eps_ij -> det(S) eps_ij / (s_i s_j)`, which is also used in 3D now). Checked against the 3D
eigenmode solver on a periodic cell of a PEC guide filled with a tilted uniaxial crystal (eps_o = 4, eps_e = 6.5, tilt 40 degrees): the four lowest modes
agree to 0.04 to 0.4 % (the 3D mesh is coarse), while leaving out the coupling makes them differ by 0.8, 5.9 and 7.3 % (test `aniso_optics.rs`). The moving-boundary
integral uses the anisotropic form built from `M = [[eps_TT - eps_Tn eps_nn^-1 eps_nT, eps_Tn eps_nn^-1], [eps_nn^-1 eps_nT, -eps_nn^-1]]` acting on the continuous
`(E_tangential, D_normal)`, and the photoelastic term uses `-(1/eps0) conj(D_s) (p:S) D_p`, exact for any symmetric permittivity (both reduce to the isotropic forms).
Ring optics use diagonal (r, z, theta) tensors: a crystal is averaged over rotations about the axis, as for the elastic and photoelastic tensors.

## Saved runs

Every Compute writes its results to a folder (Study tab > Run folder; default `Documents\nyx_runs`): `temp_untitled_<date>_<time>.hdf5` with the
XDMF description (builds with HDF5 only), the project as `temp_untitled_<date>_<time>.nyx.xml` and the log of what the solver did as
`temp_untitled_<date>_<time>.log.txt`, so nothing is lost before a project is saved. The file names use the local time.

## Linear solvers

Small systems use a sparse Cholesky factorisation. Large ones use conjugate gradients preconditioned by
smoothed-aggregation algebraic multigrid (own Rust code): rigid-body modes as near-nullspace, Chebyshev
smoothing with 3x3 node blocks, single-precision matrices inside the cycle, and for quadratic elements the
linear subspace as first coarse space. `Automatic` switches on at about 30,000 free unknowns in 3D (400,000
in 2D) and falls back to the factorisation if multigrid does not converge. It serves static solid mechanics,
steady and transient heat (one hierarchy for all time steps), electrostatics and the lowest eigenfrequencies
(LOBPCG with the multigrid cycle as preconditioner, no factorisation; used when a probe shows multigrid
converging fast, otherwise the factorisation). Measured on this machine (20 cores), 3D:

| model | unknowns | direct | multigrid |
|---|---|---|---|
| drilled plate, static, P2 | 189,804 | 27.6 s | 1.7 s (15 iterations) |
| block, static, P2 | 680,943 | 340.8 s | 7.1 s (16 iterations) |
| sapphire cantilever, static, P2 | 280,755 | 67.1 s | 14 s (189 iterations) |
| drilled plate, 6 eigenfrequencies, P2 | 189,804 | 43.6 s | 16.8 s |

Very thin parts (two element layers, such as the sapphire plate) make multigrid slow with quadratic
elements; it still wins for static solves, and the eigenfrequency study then chooses the factorisation.

## 3D mesher (Solids)

The Solids geometry kind has been removed from the application (every geometry is now a `geom` model); the 3D solid mesher is kept only as the test fixture `crates/femcore/tests/common/solid3d` (it was `femcore::mesher3d`) and builds tetrahedral meshes of boxes, spheres, cylinders and extruded outlines, painted
in order (a later body overrides earlier ones, domain 0 is a hole). The method is isosurface stuffing
(Labelle and Shewchuk, 2007, as implemented in Quartet): fill the solid with the tetrahedra of a
body-centred cubic lattice, move lattice vertices that lie within alpha of an edge length of the surface onto
it, and cut the few tetrahedra that remain across it. Every step is local, so it runs on all cores: a sphere
with a hole of 25.7 million tetrahedra took 1.3 s (20 threads), 1.7 million took 0.11 s. The design follows
what makes the fast open-source meshers fast (lattice or Delaunay kernels, parallel blocks, no global
search structures; HXT, Quartet), not their code.

**Fitted interfaces and edges.** A material interface between bodies of different domains is cut
exactly: the faces of the bodies are extra cutting surfaces, so a planar interface at x = 0.37 on a
lattice of 0.1 gives domain volumes exact to rounding, no cell straddles it, and a curved interface (a ball in
a block) converges faster than the element size. A sharp edge or corner of a body that is not on the lattice
(a re-entrant edge of an L-shaped solid, the rim of a plate) gets vertices moved onto it by alternating
projections onto the faces that meet there, so the volume of the L-shape is exact (`snap_features`, on by default).

**Local size and thin bodies.** A body can carry its own element size (`Body3d.h`). The mesh is then refined
inside it by longest-edge (Rivara) bisection: a cell whose longest edge is larger than 1.3 times the wanted size
is split together with every cell around its terminal edge, so the mesh stays conforming and grades
smoothly into the coarse part. New vertices on a curved boundary or interface are moved onto the true surface
(Newton steps on the signed distance), and boundary facets and labels are split with their cells. A body
thinner than `min_layers` (default 2) elements along an axis makes the lattice finer along that axis (at most 4
times finer than h), so a film of 0.05 on a lattice of 0.2 still holds two layers.

Checked: boxes exact (volume, face areas), sphere volume converging with h and a closed surface,
all cells positively oriented and conforming, smallest dihedral angle above 8 degrees on the sphere,
boundary numbers, a drilled plate whose peak von Mises stress matches the textbook stress
concentration factor within about 1 percent; exact volumes for off-lattice planar interfaces and edges;
refined bodies with exact volumes, conforming meshes, interface vertices on the sphere to 1e-8, and the
longest edge in a refined body within the requested size. Limits: refinement is by bisection of the
lattice mesh (it cannot find a feature the lattice missed, which is what `min_layers` is for), thin gaps
between bodies are not detected, at most 63 body faces take part in interface fitting, and curved
sharp edges (the rim of a cylinder) are fitted by the same projections only where the faces are planes,
cylinders or spheres.

## Results files and projects: HDF5, XDMF and XML

Each simulation can write its results after Compute (Study tab > Results file) and `File > Open`
reads such files back. The writer is `femcore::h5io`.

* **`results.h5`** (HDF5): `/mesh` (points, cells, facets and labels), `/spec` (the problem as JSON
  bytes), `/results/solution_k/fields/<name>` (and `disp`) with units as attributes, and `/viz`
  (the same data in the shape ParaView needs: quadratic cells in VTK node order, 3-component vectors).
* **`results.xdmf`** (XML written with `quick-xml`, next to the `.h5`): the light-data description
  that ParaView and VisIt open; the heavy data stay in the HDF5 file. Several solutions (time steps)
  become a temporal collection. Checked by reading the files back with meshio (triangle6/tetra10,
  fields equal to the HDF5 data, time values correct); the midside node order is also tested in Rust.
* **`project.nyx.xml`**: the whole model (geometry, materials, boundary conditions, studies,
  output settings) as readable XML, written and read with `quick-xml`; `.nyx.json` still works.
  Round-trip tested on all examples.

**Dual-mode I/O:** one writer, two drivers, chosen per simulation (and per `--output` flag).

* *Serial HDF5* (feature `hdf5`): needs only the HDF5 library. `set HDF5_DIR=C:/path/to/hdf5`, then
  `cargo build --release -p nyx --features hdf5`. `dist/` has a ready build with its DLLs.
* *Parallel HDF5, MPI-IO* (feature `phdf5`, implies `hdf5`): every MPI process makes the same
  collective HDF5 calls and writes its own block of rows of every dataset, so the file equals a
  serial one. Builds without the feature contain no MPI code and report that parallel output is not
  available. **The `phdf5` code could not be compiled or run on this machine** (no MPI SDK, and the
  installed HDF5 has no MPI support). What is tested: the row partitioning, and a replay in which three
  emulated ranks each write only their rows into one file and the result equals the serial file.
* Headless runs for a workstation or cluster, no window: `nyx --batch model.nyx.xml
  --output serial out.h5`, or `mpiexec -n 16 nyx --batch model.nyx.xml --output parallel out.h5`.
  Every process solves the whole problem (the solvers themselves are not distributed over MPI) and
  writes its own rows of the output.

## Using the application

* **File > New (Ctrl+N)** starts a simulation from scratch: choose the physics and 2D or 3D (a 3D project is the 2D section swept along y); you get a plate as a `geom` model, default materials and no boundary conditions, then fill in the tabs below and press Compute (F5).
* **Geometry** - one `geom` model (regions and boundaries are named; a selector such as `outer & plane(x = 0)` picks boundaries), an element size, and optionally a sweep along y into 3D layers. The rectangle, block, shapes and solids kinds no longer exist; the removed ones in an old project are named in the error when it is opened.
  **File > Open** (or drag a file in) reads saved projects (`.nyx.xml`, `.nyx.json`) and HDF5 results.
* **Materials** - per domain, with a small library (Si, SiO2, Si3N4, sapphire (isotropic), Al,
  Cu, steel, water, FR-4, air). Values are shared between the physics that use them.
* **Physics** - boundary conditions per physics. *Pick on mesh* lets you click boundaries in the
  viewport. Solid mechanics also has plane stress/strain, thickness and thermal expansion (from a
  uniform temperature or from the last heat-transfer result).
* **Study** and **Compute** (F5) - runs in the background; results appear with a field selector,
  colour maps, deformation, a slider over modes or time steps, a summary plot and the numbers.
* **File > Export** writes every dof with all fields as CSV. **Python (optional)** runs a script
  of yours (any `.py` file you choose; none is shipped with the program) with the project JSON and the result CSV as arguments.

Examples (menu *Examples*) cover each physics. For automated checks the app can be driven from
the command line and takes screenshots:
`nyx --size 1400x860 --script "example 1; compute; shot out.png; quit"`.

## What is solved (and how it was checked)

| Physics | Features | Checked against |
|---|---|---|
| Electrostatics 2D/3D | per-domain eps_r, space charge, potential/ground/surface charge, energy, capacitance | exact linear potential between plates (2D and 3D), series-capacitor closed form |
| Heat transfer 2D/3D | per-domain k, rho, Cp, Q; temperature, flux, convection; stationary and transient (exact semi-discrete solution, BDF2 for large models); initial T as an expression | exact 1D profiles (convection, volumetric source), decaying sine, lumped cooling |
| Solid mechanics 2D (plane stress/strain) and 3D | per-domain E, nu, rho, alpha; fixed (per component, with prescribed values), boundary loads, thermal expansion; static, eigenfrequencies; stress components and von Mises | Euler-Bernoulli plus shear tip deflection and eigenfrequencies of a cantilever (2D and 3D), uniaxial stress at the cell centres, free thermal expansion, layers in series |
| EM waves, mode analysis (2D cross-section) | full-vector Nedelec(2)+P2 formulation, inhomogeneous media, PEC walls, real or evanescent targets | closed forms of filled and empty WR-90 (propagating and evanescent modes, Ez shape), a dielectric-slab transcendental equation |
| PML (modes) | complex stretching along x/y, complex neff and Q | no-PML limit equals the real solver; real stretch = wider guide; constant complex stretch = analytic complex neff; absorbing slab keeps the guided mode |
| Anisotropy | cubic/trigonal/... stiffness with rotation, tensor k, eps, alpha | isotropic limit of the cubic path equals the isotropic solution; cubic plate modulus vs the compliance formula at 0 and 45 degrees; sapphire invariance under 120 degrees; stress-free anisotropic thermal expansion; anisotropic capacitance; TE10/TE01 modes in a biaxial filling |
| EM waves, 3D eigenmodes | Nedelec edge elements of degree 1 and 2 on tetrahedra, inhomogeneous and complex media, PEC walls, PML stretching, gradient (null space) modes removed by projection | box cavity frequencies, sphere cavity modes, first-order convergence rate, real stretch = longer cavity |
| Floquet-periodic boundaries | phase-shifted constraints between matching boundaries for EM modes (2D and 3D), 2D/3D solids and solids with an out-of-plane wavenumber; Bloch dispersion | waveguide dispersion, a 3D rod against the 2D out-of-plane wavenumber solution |
| Low-reflecting boundaries (solids) | impedance damping of the boundary (quadratic eigenproblem in the frequency, linearised) | analytic damped bar: frequency and Q |
| Solids with out-of-plane wavenumber kz | 2D cross-section, 3 displacement components, complex Hermitian stiffness, PML, LRB, Floquet | Euler-Bernoulli beam, Young's rod, a clamped plate |

The solvers are P1/P2 Lagrange on triangles and tetrahedra (Nedelec degree 2 for the waves),
with sparse Cholesky/LU from `faer`, and a block shift-and-invert Krylov eigensolver written for
this project (it finds repeated eigenvalues such as the degenerate modes of a square membrane).

## Limits

* Linear materials only (no contact or nonlinear solids). Perfectly matched layers are available for the waveguide modes and for solid eigenfrequencies (complex coordinate stretching, complex neff or frequency); low-reflecting boundaries need a search frequency (they make the problem quadratic) and cannot be combined with periodic conditions yet. Anisotropic data is supported: crystal stiffness (cubic, hexagonal, trigonal, orthotropic) with Euler-angle orientation, and tensor conductivity, permittivity and thermal expansion. For the waveguide modes the crystal must be oriented with the propagation axis z as a principal axis.
* Meshes: the app meshes `geom` models (2D, optionally swept along y). Quadrilaterals, hexahedra and prisms given as a `CellMesh` are split into triangles/tetrahedra (exact for straight-sided cells; a hexahedron becomes 6 tetrahedra (24 with the old star split, see above)). Pyramids are not supported (view only). The 2D mesher (`femcore::meshing::mesher`) meshes 2D rectangles, circles and polygons (overlapping shapes, holes, per-shape element size, boundaries numbered per shape edge) with constrained Delaunay triangulation and Ruppert refinement (the `spade` crate), and can extrude the 2D mesh through layers into 3D tetrahedra. the 3D solid mesher (boxes, spheres, cylinders, extruded outlines; see the 3D mesher section) is a test fixture under `crates/femcore/tests/common/solid3d`; the application no longer offers it.
* Dense eigendecomposition is used for transient heat up to about 2500 free dofs.
