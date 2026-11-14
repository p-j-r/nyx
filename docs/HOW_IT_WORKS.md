# How nyx works, from a click on Compute to a picture

This is a tour of what the program does, in the order it does it, with the name of the file in
`crates/femcore/src` (the solver library) or `crates/app/src` (the window) that does each step. No
background in finite elements is assumed; every term is explained where it first appears. The numbers
quoted were measured on the development machine (20 logical processors, 31 GB of memory); they are
examples, not promises.

## 1. The big picture

```
  what you enter                        what the program does                             what you get
  ─────────────                         ───────────────────────                           ────────────
  geometry  ──► mesh  ──► small matrices per triangle/tetrahedron ──► one big sparse matrix ──► solve ──► fields
  materials                (element matrices, §4)                      (assembly, §5)           (§6, §7)    stresses
  boundary conditions ───────────────────────────────────────────────►  constraints (§5)                   frequencies
  study (static, eigenfrequencies, modes, ...)                                                             gains (§9)
```

Everything physical (a bending beam, the capacitance of two plates, the optical modes of a waveguide) is
described by a **partial differential equation**: a rule that links a field (displacement, temperature,
electric field) to its own rate of change in space. A computer cannot handle a field that has a value at
every point, so the **finite element method** replaces it by a limited number of unknowns: the field is
cut into small pieces (triangles in 2D, tetrahedra in 3D, together the **mesh**), and inside each piece the
field is assumed to follow a simple shape fixed by a few numbers (its values at the corners, and for
"P2" elements also at the middle of each edge). The unknowns are those numbers. The equation, applied
piece by piece and added up, becomes a large system of ordinary equations

```
K u = f          (statics: stiffness matrix K, unknowns u, forces f)
K u = ω² M u     (vibrations: mass matrix M, frequencies ω)
```

and the rest of the program is about building `K` and `M` quickly and solving these systems quickly.

## 2. What you give it: the problem (`problem.rs`)

A **problem** (`ProblemSpec`) says: which physics (electrostatics, heat transfer, solid mechanics, EM
waves, Brillouin scattering); the element order (1 = values at the corners, 2 = also at edge middles,
more accurate); the **materials** per domain (numbers such as Young's modulus E, or full 6×6 stiffness
tables of a crystal with its orientation angles); the **boundary conditions** (which boundary is fixed,
loaded, held at a potential, ...); and the **study** (static solve, eigenfrequencies, waveguide modes, a
transient run over time). A project file (`.nyx.xml`) is exactly this plus the geometry.

There is one public coordinate convention, the *wafer frame*: a waveguide cross-section is the x-z plane, y is the
direction the light travels, z is the surface normal. The 2D solver kernels read the problem in these axes and name their
result fields in them (`uz`, `szz`, ...); `femcore::section` is the one place with axis numbers (a 2D mesh is the x-z plane,
the out-of-plane axis is y, and everything propagates along +y).

## 3. The mesh (`mesh/`, `meshing/`, `adapt/`)

A `SimplexMesh` is a list of points, a list of cells (3 points per triangle, 4 per tetrahedron) with a
**domain number** each, and a list of **boundary facets** (edges in 2D, triangles in 3D) with a
**boundary number** each. Domain numbers say where a material applies; boundary numbers say where a
condition applies (for a rectangle 1 = left, 2 = bottom, 3 = top,
4 = right).

Meshes come from four places: generated rectangles and blocks (`mesh/mod.rs`); 2D shapes (rectangles, circles,
polygons, with holes and per-shape element size) triangulated with refinement by `meshing/mesher.rs`, optionally
extruded into 3D layers; 3D solids (boxes, spheres, cylinders, extrusions) filled with well-shaped tetrahedra by
`tests/common/solid3d` (isosurface stuffing, with local refinement; a test fixture since phase 2, not part of the library). Meshes given as mixed cells (`mesh/cells.rs`) are accepted too: cells that are not
simplices (quadrilaterals, hexahedra, prisms) are cut into triangles/tetrahedra without adding points.

**Unknown numbering.** For P1 the unknowns are the mesh points. For P2 they are the points followed by one
per edge (its midpoint). A vector problem such as elasticity has one block of these per component: all x
displacements, then all y, then all z.

## 4. Element matrices (`fem.rs`, `material.rs`)

On one triangle or tetrahedron the field is a sum of **shape functions** times the unknowns. Shape
functions are written with barycentric coordinates `l_i` (the weights with which a point is a mixture of
the cell's corners): P1 uses `l_i`; P2 uses `l_i(2 l_i − 1)` at corners and `4 l_i l_j` at edge middles.
Cells are straight-sided, so the gradients of the `l_i` are constant over a cell (`fem::geometry`).

An **element matrix** is an integral over the cell, evaluated by **quadrature** (a few well-chosen points
with weights; exact for polynomials up to degree 2 or 4, `fem::quad_rule`). For elasticity
(`sms.rs`) it is `∫ Bᵀ D B dV`: `B` turns the unknowns into strains (how much the material is stretched
and sheared), `D` is the material's stiffness table (6×6, in the order xx, yy, zz, xy, yz, xz, with
engineering shear strains), so the integral is the elastic energy stored per unit of unknowns. For
crystals `D` is rotated into the wafer axes by the orientation angles (`material.rs`, Z-X-Z Euler angles).
The mass matrix is `∫ ρ Nᵀ N dV`. Other physics have their own versions (`dfs.rs`: diffusion matrices
for electric potential and temperature; `ems*.rs`: curl-curl matrices with edge elements, see §8).

Variants of the elasticity element matrix: out-of-plane wave number (`sms_kz.rs`: all three displacement
components live on a 2D cross-section with a travelling-wave factor `exp(−i kz z)`), bodies of revolution
(`sms_axisym.rs`), complex coordinate stretching for absorbing layers (`sms_pml.rs`), and boundary
damping that absorbs outgoing waves (`sms_lrb.rs`).

**Axisymmetric eigenfrequency study (`nyx-sms`, ADR 68 and 69).** A body of revolution is a 2D model whose x is the
radius r; the displacement goes as `exp(−i m θ)` with the azimuthal order `m` and three components (radial, axial,
azimuthal). The study `axisymmetric_eigenfrequency` takes `n_modes`, `f_search` and `azimuthal_order`. A mesh may touch
the axis `r = 0`: there the kernel holds the components that regularity forces to zero (`m = 0`: radial and azimuthal;
`|m| = 1`: axial; `|m| ≥ 2`: all three), and the log line `axis r = 0: N dofs, K held at zero` says how many. Absorbing
layers (on +x, +z, −z, never towards the axis) are the *cylindrical* stretch: the radius becomes the complex `r̃` and the
`1/r` terms and the weight `r` carry `r̃/r`, as in the electromagnetic ring. The frequencies are then complex; the result
carries `fN (Hz)`, `QN = Re f / (2|Im f|)` and `layer fraction N`, the share of the mass that lies in the layers (modes
with a large fraction belong to the layers, not to the resonator, and have a Q near 1).

## 5. Assembly and boundary conditions (`fem::assemble_parallel`, `sparse.rs`, `periodic.rs`)

**Assembly** loops over all cells, computes each element matrix and adds its entries into the global matrix
at the rows/columns of the cell's unknowns (a point shared by several cells collects contributions from all of
them). The program collects entries as (row, column, value) triplets (`Coo`) and sums duplicates when it
converts to compressed rows (`Csr`). The cell loop runs on all processor threads, each handling a block of cells and
folding its finished sub-blocks into a running sum so that memory stays modest (a 3-million-unknown solid problem
peaks at about 13 GB).

**Boundary conditions.** A *load* (force per area) adds a boundary integral to the right-hand side `f`. A
*prescribed value* (a fixed face, a potential, a temperature) is applied by elimination: the constrained
unknowns are removed from the system and their known values move to the right-hand side
(`fem::solve_dirichlet`; the multigrid path instead replaces those rows by scaled identity rows). *Periodic*
(Floquet) conditions tie the unknowns of one boundary to those of the opposite boundary times a phase factor
(`periodic.rs`): the system is reduced to the free unknowns by a sparse change of variables. *Absorbing layers*
(`PML`) make coordinates complex in a layer so that waves entering it decay instead of reflecting; this makes
the matrices complex, and the frequencies or effective indices complex too (the imaginary part gives the loss and
the quality factor Q).

## 6. Solving `K u = f`: direct or multigrid (`linalg.rs`, `amg.rs`, `fem::solve_spd`)

*Direct* means a **sparse factorisation** (`K = L Lᵀ`, Cholesky, or LU for general matrices; `faer` does it).
It is exact (to rounding), reusable for many right-hand sides, but its memory and time grow quickly with the
size of 3D problems because the factors fill in. It cannot be interrupted once started (a cancel waits for it).

*Multigrid* (`amg.rs`) is an iterative method for large 3D problems. The idea: a simple smoothing iteration
(Gauss–Seidel / Chebyshev) removes the *rough* part of the error quickly but leaves the *smooth* part; the
smooth part can be represented on a much coarser mesh, solved there (recursively), and the correction carried
back. Here the coarse levels are built from the matrix itself (**smoothed aggregation**: neighbouring
unknowns are grouped, the coarse space is made to reproduce rigid-body motions exactly, the interpolation is
smoothed) rather than from a second mesh. One such "V-cycle" is then used as a **preconditioner** for the
conjugate-gradient iteration (`amg::pcg`), which typically converges in 15–30 iterations to a residual of
1e-10 whatever the size. A 3.1-million-unknown block takes about 40 s in total, of which the multigrid set-up and the
iterations are roughly 15 s each.

*Auto* (the default) picks multigrid for static solves above about 30,000 unknowns in 3D (400,000 in 2D) and the
factorisation otherwise; if multigrid does not converge it falls back to the factorisation and says so in the log.

## 7. Eigenproblems `K u = λ M u` (`eigen.rs`, `lobpcg.rs`, `csparse.rs`)

Free vibrations, cavity resonances and waveguide modes are eigenproblems: find the special pairs (λ, u).

* **Shift-and-invert Krylov** (`eigen.rs`, real symmetric; `csparse.rs`, complex). To find the eigenvalues
  nearest a chosen *shift* σ, the matrix `(K − σM)` is factorised once and the method repeatedly applies
  `(K − σM)⁻¹ M` to a block of vectors; the eigenvectors of interest dominate that operator, so a small
  set of vectors (the "Krylov space") is enough to describe them. Each step prints the worst residual (how
  far the current estimates are from solving the equation), and the loop ends when it is below the goal.
  Starting from a block of several vectors finds repeated eigenvalues (the two bending directions of a square
  beam, for instance). Search near zero uses a shift just below the spectrum (`f_search = 0`).
* **LOBPCG** (`lobpcg.rs`). For large 3D problems a factorisation of `K − σM` is too expensive, so the lowest modes are
  found iteratively with multigrid as preconditioner. It needs no factorisation. *Auto* tries it only above
  80,000 unknowns in 3D, 400,000 in 2D (below that the factorisation is as fast), probes whether multigrid is a good preconditioner,
  and **gives up early** when the residual stops improving (a slender part can defeat multigrid): the stage
  log then says "LOBPCG stalled ...; used a sparse factorisation instead".
* **Complex problems with absorbing layers or tilted crystals** (`csparse.rs`, `ems_pml.rs`). The matrices are
  complex symmetric. When a crystal's axes are tilted against the propagation direction, the unknown
  propagation constant β appears in the matrix itself; `ems_pml.rs` solves this by repeating the eigensolve
  with the latest β until β stops changing (relative change below `COUPLED_BETA_TOL` = 1e-8, which sits just above
  the eigensolver's own noise of ~1e-9; a stricter test used to run all 10 repeats of 13 s each for nothing).
  Low-reflecting boundaries make the eigenproblem quadratic in ω; `sms_lrb.rs` rewrites it as a linear problem of
  twice the size.

## 8. Electromagnetic waves (`ems.rs`, `ems_pml.rs`, `ems_axisym.rs`, `ems3d.rs`)

For light in a waveguide cross-section the unknown is the electric field, which has a component along
the guide and two across it. The program uses the standard **mixed formulation**: the transverse field is
expanded in Nédélec *edge* elements (unknowns on edges, tangential component continuous across cells — the
right continuity for E, which avoids false "spurious" modes) and the longitudinal field in ordinary
Lagrange elements. The eigenvalue is `−β²`; the effective index is `n_eff = β / k0`. Perfect conductor walls, complex
materials and PML are handled; `ems_axisym.rs` treats rings (azimuthal order `m`) and `ems3d.rs` closed 3D cavities
(`curl(1/μ curl E) = (2πf/c)² ε E`; the set of gradient fields, which has eigenvalue zero and is not physical, is
projected out).

## 9. Brillouin scattering (`sbs.rs`, `sbs_ring.rs`, `chips.rs`)

Stimulated Brillouin scattering couples two light waves (pump and Stokes) and a sound wave. The chain for a
waveguide cross-section:

1. Solve the **optical modes** of the cross-section (§8) at the chosen wavelength: effective indices, fields E, D, H.
2. For a chosen pump/Stokes pair, the sound wave must carry the momentum difference: its out-of-plane wave number is
   `q = β_pump + β_Stokes` (backward scattering) or `β_pump − β_Stokes` (forward).
3. Solve the **elastic modes** at that `q` (`sms_kz.rs`; frequencies are complex with PML in the substrate, the
   imaginary part gives the radiation loss Q).
4. Combine optical and elastic fields in two overlap integrals: the **photoelastic** one (the strain of the
   sound wave changes the permittivity, which scatters light; `p:S` with the material's photoelastic tensor) and
   the **moving-boundary** one (the displaced interface scatters light), each normalised by the optical powers.
   With the conventions chosen for this project the strain is conjugated (`StrainConvention::Conjugate`, the
   default; `AsIs` uses the strain as it is).
5. The **gain** per mode follows from the overlaps, the mode's Q and group velocities. The report lists, per
   elastic mode: frequency, Q, share of the energy in the body and in the PML, MB, PE and total gain, and g0.

`chips.rs` builds the four ring-chip cross-sections (silicon rib on sapphire) and `sbs_ring.rs`
does the same for whole rings (axisymmetric). Example 15 of the Examples menu is the 300TE chip.

## 10. From the solution to the picture (`fields.rs`, `app/src/viewport.rs`, `scene.rs`)

A solution is a vector of unknowns. To plot stresses, the program takes the displacement gradient at the
nodes of every cell, forms the stress there (`σ = D ε`), then **averages over the cells that share a node**
and computes von Mises stress. `Display` keeps the surface of
the mesh (the cells in 2D, the boundary facets in 3D; for P2 each triangle is cut in four so the midpoint values
show). The window colours the vertices from the field and colour map, projects them to the screen, sorts triangles
back to front, and hands one mesh to the graphics library. `scene.rs` keeps the parts that did not change between
frames, so panning a million-triangle result costs ~15 ms and rotating it ~45 ms.

## 11. What the log says

Every solve writes a plain-language log of what it did (also in the window under "What the solver did", in the
`.log.txt` saved next to each autosaved run, and printed by `nyx --batch`). Lines are indented below the
stage they belong to; a stage line gets its time when it ends. A real example, 3D cube, 6 modes, multigrid
LOBPCG:

```
solving Solid mechanics on a 3D mesh of 16464 cells: 7.40 s       <- the whole study
  setup (materials, constraints): 0.001 s                          <- materials per cell, fixed boundaries
  assembling element matrices (16464 cells, 16 threads): 0.19 s    <- §5, K
  stiffness assembly: 0.187 s
  assembling element matrices (16464 cells, 16 threads): 0.03 s    <- the mass matrix M
  building the multigrid hierarchy (73167 unknowns): 0.27 s        <- §6, coarse levels
    multigrid levels (unknowns): [73167, 10125] + coarsest 750
  LOBPCG eigensolver with multigrid preconditioner (6 modes, 70644 unknowns): 6.38 s
    iteration 0: worst residual of the 6 wanted modes 8.4e-1 (goal 1e-8)
    iteration 10: worst residual of the 6 wanted modes 6.0e-5 (goal 1e-8)
```

And from the Brillouin example (excerpt; the full log is 60 lines): the dominating costs are the factorisations
(`factorising the complex matrix (161625 unknowns); this cannot be interrupted: 4-10 s`), one per optical
solve and one for the elastic solve; the tilted-crystal refinement now stops after one step:

```
optical modes of the cross-section: 16.05 s
  complex eigenmodes (1 wanted, 161625 unknowns): 5.11 s
  ...
  tilted-crystal coupling, optical mode 1: refinement step 1, relative change of the propagation constant 2.7e-9 (stop below 1e-8)
pump 1 -> Stokes 1: phase matching, elastic modes and gains: 27.76 s
  complex eigenmodes (24 wanted, 101934 unknowns): 24.85 s      <- the elastic modes at the phase-matched q
    Krylov space 216 vectors, worst residual of the 24 wanted modes 1.2e-9 (goal 1e-9)
    Krylov space of 224 vectors: orthogonalising 5.7 s, solves 5.8 s, small eigenproblems 1.7 s, total 23.8 s
```

Reading hints: a residual that falls by a factor 10–1000 every few steps is healthy; a residual that stays
flat means the method is struggling (the program will stop it or fall back). `this cannot be interrupted`
marks a factorisation: Cancel takes effect when it ends. Lines starting with "multigrid did not converge"
or "LOBPCG stalled" mean a fallback to the direct method happened.

## 12. Where results go

* In the window: the field selector, the colour map, deformation, the slider over modes/time steps, the summary
  table and the log. Hover shows the value under the cursor.
* **Autosave** (Study tab, on by default): every Compute writes, into the runs folder
  (`Documents\nyx_runs`), `temp_untitled_YYYYMMDD_HHMMSS.nyx.xml` (the project), `.log.txt` (the log) and, in builds
  with HDF5, `.hdf5` + `.xdmf` (all fields; open the `.xdmf` in ParaView). A results file can also be named in the Study tab or
  with `--output serial FILE.h5`.
* **Headless** (`nyx --batch project.nyx.xml`): solves without a window, prints the scalars and the log.
* Export: *File > Export* writes every unknown with all fields as CSV; *Python (optional)* runs a script of your own
  with the project JSON and the CSV as arguments (no `.py` file ships with the program).

## 13. Threads, Cancel and diagnostics

* Parallel loops all run on **one shared thread pool** (`par.rs`; rayon). Its size is the number of logical processors
  (at most 32); set `NYX_THREADS=4` to leave cores to other programs or `NYX_THREADS=1` for a serial run.
  Sparse factorisations use their own, smaller pool of 8 (`NYX_FACTOR_THREADS`) because they get slower beyond that.
* Results are reproducible: the work is cut into pieces that depend only on the thread count and the pieces are added
  in a fixed order, so the same run gives the same numbers every time (different thread counts can differ in the last
  digits).
* **Cancel** (toolbar button or Esc): the solver checks a flag at the top of every iteration, every multigrid level, every
  refinement step and every Brillouin pair and then returns "cancelled"; earlier results stay. Assembly of a matrix and a
  started factorisation finish first. Only one job runs at a time; mesh builds and file opens requested meanwhile wait.
* Diagnostics: `NYX_PROFILE=1` echoes the log to the console as it is written; `NYX_FRAME_LOG=1` prints the
  time of every window frame and which drawing passes ran; `NYX_LOBPCG_TRACE=1` and `NYX_AMG_*` (see `amg.rs`) are for solver work.
* **Input checks** (`validate.rs`): NaN or infinite numbers in materials, coordinates, boundary values or study
  parameters are refused with a message naming the item; a solution containing NaN/infinity is refused too.

## 14. Module guide (crates/femcore/src)

| file | what it does |
|---|---|
| `solve.rs` | the one entry point; checks, dispatches to the physics, renames fields, adds the log |
| `problem.rs` | problem and result types (see §2) |
| `validate.rs` | refuses NaN/infinity before and after a solve |
| `progress.rs` | cancel flag, current stage, stage log (thread-local, so solver functions keep their signatures) |
| `par.rs` | the shared thread pool and the factorisation pool |
| `mesh/mod.rs` | simplex meshes (the solver view), unknown numbering, structured generators, relabelling by boxes |
| `meshing/mesher.rs` | 2D shapes → triangles (constrained Delaunay + Ruppert refinement), extrusion to 3D |
| `tests/common/solid3d` | (test fixture) 3D solids → tetrahedra (isosurface stuffing; fitted interfaces and edges; bisection refinement) |
| `mesh/view.rs` | surface of any mesh for drawing, including cells the solvers cannot use |
| `mesh/cells.rs` | the one mesh model: mixed cells (triangles, quadrilaterals, tetrahedra, prisms, hexahedra), cut into simplices for the solvers; `adapt/` holds the metric adaptation |
| `fem.rs` | quadrature, shape functions, scalar matrices, constraints, direct-vs-multigrid choice |
| `fields.rs` | gradients at nodes, node averaging, plotting surface |
| `material.rs` | anisotropic stiffness/permittivity/conductivity/expansion tensors, crystal rotation |
| `expr.rs` | the small expression language for initial temperatures |
| `sparse.rs`, `csparse.rs` | real and complex sparse matrices, complex Arnoldi eigensolver |
| `linalg.rs` | sparse Cholesky/LU, dense eigensolvers |
| `amg.rs` | smoothed-aggregation multigrid and preconditioned conjugate gradients |
| `eigen.rs`, `lobpcg.rs` | shift-and-invert Krylov eigensolvers; multigrid LOBPCG with stall guard |
| `dfs.rs` | electrostatics and heat transfer (stationary, transient: exact for small systems, BDF2 for large) |
| `sms.rs` | linear elasticity: static, eigenfrequencies, stresses, thermal expansion |
| `sms_kz.rs`, `sms_axisym.rs` | elastic modes with out-of-plane wave number; bodies of revolution |
| `sms_pml.rs`, `sms_lrb.rs` | anchor loss through PML; low-reflecting boundaries (quadratic eigenproblem) |
| `sms_floquet.rs`, `periodic.rs` | periodic cells (phononic crystals), Floquet constraints |
| `ems.rs`, `ems_pml.rs` | 2D waveguide modes; complex modes, PML, tilted crystals |
| `ems_axisym.rs`, `ems3d.rs` | ring resonator modes; 3D cavities |
| `sbs.rs`, `sbs_ring.rs`, `chips.rs` | Brillouin overlaps and gain; rings; the four chips |
| `xmlio.rs`, `h5io.rs` | project XML; HDF5/XDMF results |

Window crate (crates/app/src): `main.rs` (command line), `app.rs` (state, jobs, layout), `tabs.rs` (inspector),
`state.rs` (project), `worker.rs` (background jobs), `viewport.rs` + `scene.rs` (drawing and picking),
`examples.rs`, `batch.rs`, `script.rs` (scripted runs), `colormap.rs`, `theme.rs`.

## 15. Running the tests and the benchmarks

```
cargo test --release --workspace                 # about 160 tests, ~25 s
cargo test --release --workspace -- --ignored    # 7 slow Brillouin tests, ~70 s
cargo clippy --release --workspace --all-targets -- -D warnings
```

Examples (`cargo run --release -p femcore --example NAME`; `stage_log` and `bench` are in `nyx-engine`, since they drive every physics): `stage_log` (the logs quoted above), `bench` (timings on larger
models), `profile3d` (stage times of a 3D static solve), `mem_profile` (peak memory by stage), `pool_bench` (matrix-vector
product and multigrid iteration cost), `lu_scaling` (factorisation speed against threads), `eig_beam` (eigenfrequencies of a
beam with each solver), `chip_report` / `chip_digits` (the Brillouin chips; `chip_digits` prints 12 digits for before/after
comparisons). `nyx --size 1400x860 --script "example 1;
compute; shot out.png; quit"` drives the window (commands in `script.rs`).
