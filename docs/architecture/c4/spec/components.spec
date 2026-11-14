# Component spec of the nyx C4 model (level 3).
#
# One section per container (= crate): [key | crate directory under crates/ | Rust crate name].
# One line per component:  id | Title | files | ADRs | description
#   files: space separated, relative to the crate's src/. `dir/` takes every file under it, `file.rs#a,b,Type::m` takes only
#   those items of a file shared by two components (free items by name, methods as Type::method).
# `edge A -> B | label` adds a relationship the source scan cannot see (components by id, `key.id` across containers).
# tools/gen.pl reads this, scans the Rust source, and writes the Structurizr fragments and every code page.
# Every source file of a crate must belong to exactly one component; gen.pl fails otherwise.

[app | app | nyx]
window    | Window | main.rs app.rs state.rs | 23, 57 | The application object around one Engine: window layout, the run-time state (mesh on show, result on show, plot settings) and the frame update loop.
shell     | Shell and dashboard | shell.rs dashboard.rs wfview.rs sheets.rs examples.rs | 14, 23 | The tab bar, the dashboard of workflows with live status, the workflow view (concept, build, analysis), the sheets (new, search, versions, share, settings) and the example starting points.
inspector | Inspector tabs | tabs.rs forms.rs geomui.rs | 23, 57 | The tabs beside the viewport: geometry editor, materials, physics, study, results, analysis. Forms are drawn from the physics descriptors; the geometry editor generates one form per feature kind.
viewport  | Viewport | viewport.rs scene.rs colormap.rs | 23 | The 2D/3D viewport: camera, drawing of meshes and result fields, picking, colour maps, and what is kept between frames.
style     | Look and widgets | theme.rs kit.rs | 23 | The nyx look (palette, fonts, egui visuals, light and dark) and the small brand widgets (cards, chips, buttons, sparklines, sheets).
headless  | Scripts and batch | batch.rs script.rs pyscript.rs | 23 | Runs without a person at the window: scripted runs with screenshots (--script), batch runs for clusters, and the optional user Python script.

[mcp | nyx-mcp | nyx_mcp]
protocol  | Protocol loop | main.rs#main,usage,serve,at_least,session_id,default_dir,NAME,VERSION,PROTOCOLS,Server,Server::new,Server::header,Server::initialize,Server::handle | 15, 57 | Newline-delimited JSON-RPC 2.0 over stdin/stdout: initialize with protocol negotiation, method dispatch, notifications, logs on stderr. No async runtime and no SDK.
tools     | Tool bridge | main.rs#tool_list,destructive,render,Server::run,Server::tools_call | 15, 57 | The tool list is the engine's command catalog (with annotations: read-only, destructive); a tool call is Engine::call, its reply rendered as text.
edge mcp.tools -> engine.commands | One tool per command: catalog and call
edge mcp.protocol -> engine.registry | Server instructions text

[connect | nyx-connect | nyx_connect]
hosts     | Host table | hosts.rs env.rs | 28 | Where each AI host (Claude Code, Claude Desktop, ...) keeps its MCP config and what an entry looks like; the home folder, application-data folders and the operating system.
config    | Host config | config.rs | 28 | Reads and rewrites a host's JSON config: status of the nyx entry, connect, disconnect, list.
server    | Server probe | server.rs | 28 | Finds the nyx-mcp program, builds the snippet for a host that is not in the table, and does a real MCP handshake against it.
report    | Report | report.rs lib.rs | 28, 19 | The machine-readable form of what nyx-connect finds and does, with its JSON Schema.
cli       | Command line | main.rs | 28 | The command-line twin of the Connect menu: list, connect, disconnect, test, optionally as JSON.

[engine | nyx-engine | nyx_engine]
commands  | Commands | commands.rs engine.rs guide.rs | 57, 38, 15 | The Engine and its catalog: every action by name with JSON arguments, JSON Schema and JSON reply; the dispatcher; the guide cards a host reads on demand. The window and the MCP server use the same entries.
project   | Project | project.rs | 63 | The saved state: feature list, mesh settings, one ProblemSpec per physics, project materials. The format key is checked at the file boundary.
geometry  | Geometry operations | geometry.rs geomdesc.rs | 38, 3 | Every edit of a geom feature list as a function that does it or says why not, and the descriptors of the geom operations (they live here because geom depends on serde only).
jobs      | Job runner | jobs.rs solve.rs convergence.rs | 61 | Mesh and solve jobs on their own threads with progress and cancel, solving physics in the order of the one-way couplings, mesh convergence on a scratch engine.
meshing   | Meshing seam | meshing.rs pml.rs | 54, 50 | MeshRequest to MeshedModel (mesh, names, label resolver); the only caller of the femcore mesher; absorbing layers across the two halves of a project.
materials | Materials | materials.rs | 58 | The material library plus the materials a project carries, resolved at the study wavelength.
analyse   | Analyses | analyse.rs | 56 | Runs the analyses of nyx-analysis on the latest result of a physics.
resultsIO | Results files | h5io.rs xdmf.rs results_file.rs | 62, 63 | HDF5 results (serial or MPI-IO), the XDMF description, and opening a results file as a project with its results.
storeLink | Store link | store_link.rs | 39, 59 | Binds the project to a workflow of the workspace, commits versions, records runs (card and metrics made here).
registry  | Registry and facade | lib.rs | 25 | The physics compiled in (cargo features dfs, sms, ems) and the facade that re-exports what a front-end draws with.
edge engine.project -> @projectFile | Saves and loads
edge engine.resultsIO -> @resultsFile | Writes and opens
edge engine.materials -> @materialLibrary | Loads

[physics | nyx-physics | nyx_physics]
trait        | Physics trait and registry | physics.rs registry.rs error.rs lib.rs | 25 | The Physics trait (everything the rest of nyx may ask of a physics), the Registry of the physics a build contains, and the plain-data error every method returns.
envelope     | Problem envelope | spec.rs decode.rs validate.rs | 25, 41, 63 | The one shape every physics is given (ProblemSpec), its decoding into typed parameters with descriptor defaults and name resolution, and its structural validation.
descriptors  | Descriptors and schemas | desc.rs schema.rs | 25 | What a physics tells nyx about itself (studies, conditions, parameters) and the JSON Schemas generated from it for the window, the MCP server and the documentation.
result       | Result and frame | result.rs frame.rs | 25, 56 | SolveResult as plain serialisable data, and the frame a result states for the analysis of its vector fields.
testsupport  | Test support | testkit.rs contract.rs | 10, 19 | The golden-file harness and the contract-test helpers shared by every physics crate.

[dfs | nyx-dfs | nyx_dfs]
adapter  | Physics adapter | lib.rs | 41, 25 | The Physics implementation of dfs: descriptors, validation, solve entry.
problem  | Typed problem | problem.rs | 41 | The typed dfs problem the kernels read, decoded from the ProblemSpec envelope.
heat     | Heat kernel | heat.rs | 41, 9 | Heat transfer on P1/P2 simplices, stationary and transient.
electro  | Electrostatics kernel | electro.rs | 41, 9 | Electrostatics on P1/P2 simplices with fixed-potential and surface-charge conditions.
common   | Shared helpers | common.rs | 41 | Per-cell gradients and values shared by the heat and electrostatics kernels.

[sms | nyx-sms | nyx_sms]
adapter       | Physics adapter | lib.rs | 42, 25, 68 | The Physics implementation of sms with three studies (static, eigenfrequency, axisymmetric eigenfrequency): descriptors, validation, solve entry.
problem       | Typed problem | problem.rs | 42, 68 | The typed sms problem the kernels read, decoded from the envelope.
solid         | Solid kernel | solid.rs | 42, 9 | Linear elastic solid mechanics on P1/P2 triangles and tetrahedra: static solve and eigenfrequencies.
crosssection  | Axisymmetric study kernel | cross_section.rs | 43, 68, 69 | Axisymmetric study kernel: constraints and axis conditions, absorbing layers, complex eigensolve and fields on the (r, z) half plane (the shared cross-section stages).
axisym        | Axisymmetric kernel | axisym.rs | 43, 69 | Elastic modes of a body of revolution on its half-plane cross-section, with the axis regularity conditions and the cylindrical PML stretch.
floquet       | Floquet kernel | floquet.rs | 42 | Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals).
pml           | Absorbing-layer kernel | pml.rs | 42 | Eigenfrequencies of elastic bodies with perfectly matched layers.
lrb           | Low-reflecting boundary | lrb.rs | 42 | Low-reflecting boundaries for solid mechanics.

[ems | nyx-ems | nyx_ems]
adapter   | Physics adapter | lib.rs | 45, 25 | The Physics implementation of ems: descriptors, validation, solve entry.
problem   | Typed problem | problem.rs | 45 | The typed ems problem the kernels read, decoded from the envelope.
modes     | Mode kernel | modes.rs | 45, 9 | Full-vector mode analysis of a 2D waveguide cross-section (inhomogeneous media allowed).
modescx   | Complex-mode kernel | modes_complex.rs | 45 | Complex-valued mode analysis: lossy media and perfectly matched layers.
axisym    | Axisymmetric kernel | axisym.rs | 46 | Electromagnetic eigenmodes of a body of revolution (rings, disks, coaxial cavities).
cavity    | 3D cavity kernel | cavity.rs | 45 | Electromagnetic eigenmodes of 3D cavities and resonators.
optical   | Optical mode sets | optical.rs | 45, 21 | The mode solvers' results as an OpticalSet, with the material data needed to evaluate their fields.

[analysis | nyx-analysis | nyx_analysis]
api           | Analysis API | lib.rs | 56 | The entry points an engine calls: an analysis is a function of a SolveResult and its mesh.
samples       | Field samples | samples.rs | 56 | The modes of a SolveResult as weighted samples of the field, which is what the polarization analysis integrates.
polarization  | Polarization | polarization.rs | 56, 21 | Stokes parameters, Poincare sphere and the classification of an optical mode.
modes         | Mode report | modes.rs | 56 | One report per mode, a table, the Poincare picture and the degenerate groups.

[store | nyx-store | nyx_store]
model     | Data model and card | model.rs card.rs | 14, 39 | The workflow manifest, experiments, run digest and the concept card generated from an experiment's spec.
workspace | Workspace | workspace.rs lib.rs | 14 | The folder behind the dashboard: settings, the list of workflows (manifests only) and the sections of the dashboard.
flow      | Workflow folder | flow.rs | 14, 59 | One workflow folder: workflow.json, history.jsonl, notes, specs and one folder per experiment.
versions  | Versions | versions.rs | 14, 59 | Every change to an experiment's set-up kept in history.jsonl with the hash of its spec.
runs      | Runs | runs.rs | 64 | Append-only run records in runs.jsonl, the digest derived from them and the filter language of the runs table.
search    | Search | search.rs | 14 | Search over the manifests: plain words with a synonym table and numbers with units matched against card dimensions.
archive   | Archives | archive.rs | 14 | .nyx files: a workflow as one deterministic ZIP and the safe way back in.
util      | Utilities | util.rs | 14 | UTC time without a time crate, time-ordered ids, atomic file writes, hashes and slugs.
edge store.flow -> @workspaceStore | Reads and writes workflow folders
edge store.runs -> @workspaceStore | Appends runs.jsonl

[geom | geom | geom]
featurelist  | Feature list | model.rs expr.rs print.rs error.rs lib.rs | 2, 3, 63 | The serialisable feature list (single source of truth for code and GUI), its operations and validation, parameter expressions, printing as Rust builder calls, and the error type.
evaluation   | Evaluation | eval.rs diff.rs scene.rs | 2, 6 | The Evaluator turns a feature list into an evaluated model at a chosen tier with a cache; scenes instantiate profiles with transforms and patterns; change detection between evaluations for sweeps.
kernel2d     | 2D kernel | prim.rs curve.rs arrange.rs complex.rs | 3, 18 | Primitives, planar curves, the planar arrangement of all boundary pieces and the evaluated region complex.
kernel3d     | 3D kernel | stack.rs planar.rs topo3.rs | 3, 1 | Prismatic models as a layer stack, general planar-faced solids, and the edges and vertices of a 3D model.
queries      | Queries and selectors | query.rs select.rs labels.rs sim.rs | 2, 3, 7, 18 | Exact queries on an evaluated model, serialisable selectors, the LabelMap between names and the integer labels solvers consume, and simulation operations as data.
tessellation | Display tessellation | tess.rs | 3 | Tessellation for the viewport and picking only, never for meshing.
math         | Math helpers | math.rs | 3 | Small vector, box and similarity-transform helpers.

[femcore | femcore | femcore]
mesh        | Mesh | mesh/mod.rs mesh/cells.rs mesh/exchange.rs mesh/view.rs mesh/quality.rs | 4, 5, 32 | The native mixed-cell mesh, the exchange mesh the physics take, simplex meshes with their numbering of unknowns, draw-only meshes and quality statistics.
meshio      | Mesh files | mesh/io/ | 16 | The compact lossless .fmsh format for meshes and results and its entropy coder.
fixtures    | Fixtures and test kit | mesh/fixtures.rs testkit/ | 41, 10, 60 | Test meshes and the SOI wire cross-section used by the tests of every physics crate.
meshing     | Meshing facade | meshing/mod.rs meshing/request.rs meshing/sizing.rs meshing/physics_size.rs meshing/geomdrive.rs meshing/mesher.rs | 6, 17, 54 | MeshPlan and mesh_model (the one way to mesh a geom model), the serialisable request and report, sizing from the physics, and the 2D triangulation core.
meshing3d   | 3D meshing | meshing/mesher3d.rs meshing/refine3d.rs meshing/geomdrive3.rs meshing/adapt3d.rs | 40, 50, 51, 55 | The lattice tetrahedral mesher on Domain3, Rivara refinement, the adapter from a 3D geom model, and adaptive meshing driven by exact faces.
structpml   | Structured absorbing layers | meshing/pml_box.rs meshing/pml_geom.rs meshing/pml_geom3.rs meshing/pml3.rs | 7, 52, 5 | Structured PML around a box-shaped or geom-model mesh: prisms, hexahedra, quads, labels.
meshcheck   | Meshing validation harness | meshing/validate/ | 7 | A small independent scalar Helmholtz solver with PML used to validate meshes.
metric      | Metric fields | adapt/mod.rs adapt/metric.rs adapt/field.rs adapt/field3.rs adapt/work.rs adapt/timer.rs | 7, 5 | Riemannian metrics in 2D and 3D, metric fields, the simplex working mesh and per-phase timers.
adaptops    | Adaptation operators | adapt/ops.rs adapt/ops3.rs adapt/topo.rs adapt/topo3.rs | 7, 51 | Edge split, collapse, flip and smoothing on the working meshes in 2D and 3D, with the vertex classification.
adaptdrive  | Adaptation drivers | adapt/drive.rs adapt/drive3.rs adapt/geom.rs adapt/geom3.rs | 7, 51 | Passes of operators until the mesh conforms to the metric, with geometry oracles saying where a boundary vertex may go.
estimators  | Error estimators | adapt/estimate.rs adapt/estimate3.rs adapt/maxwell.rs | 7, 17, 53 | Error estimation and the solve-estimate-adapt loop; the Maxwell residual estimator for edge elements.
lacore      | Matrices and vectors | linalg/mod.rs linalg/csr.rs linalg/csr/ linalg/bsr.rs linalg/basis.rs linalg/vecops.rs linalg/scalar.rs linalg/kernels.rs linalg/par.rs linalg/rng.rs linalg/linop.rs linalg/mm.rs linalg/bundle.rs linalg/dense/ linalg/faer_glue.rs linalg/maxwell.rs | 33, 8, 11 | The one matrix type Csr over a generic scalar, block storage, multivector kernels, operator interfaces, dense helpers, Matrix Market, problem bundles and test-problem generators.
ladriver    | Linear algebra driver | linalg/driver.rs linalg/problem.rs linalg/pc.rs | 33 | The door the physics kernels use: a direct solve, eigenpairs near a shift, and the solver interface for models.
ladirect    | Direct solvers | linalg/direct.rs linalg/sparse_lu.rs linalg/blr.rs linalg/lowrank.rs linalg/mixed.rs linalg/ordering.rs | 33, 8 | Shifted sparse direct solves with symbolic reuse, LU with a user ordering, block low-rank factorisation, mixed precision, nested dissection.
laksp       | Krylov solvers | linalg/ksp/ | 33, 8 | CG, MINRES, GMRES, BiCGStab, COCR, QMR, IDR(s), GCRO-DR for real and complex systems.
laeigen     | Eigensolvers | linalg/eps/ linalg/lobpcg.rs linalg/jd.rs linalg/toar.rs linalg/nep.rs linalg/nleigs.rs linalg/chebfsi.rs linalg/contour.rs linalg/arnoldi_real.rs linalg/sens.rs | 33, 8, 11 | Krylov-Schur, LOBPCG, Jacobi-Davidson, quadratic and nonlinear eigenproblems, contour integral, Chebyshev filtering, real Arnoldi for the 2D modes, sensitivities.
laprecond   | Preconditioners | linalg/amg.rs linalg/amg_elastic/ linalg/ams.rs linalg/ads.rs linalg/schwarz.rs linalg/fieldsplit.rs linalg/ilu.rs | 33, 8 | Smoothed-aggregation multigrid (scalar and elastic), auxiliary-space Maxwell and divergence solvers, Schwarz, field-split and incomplete factorisations.
laother     | Time, nonlinear, sweeps | linalg/ts/ linalg/snes.rs linalg/opt.rs linalg/expm.rs linalg/sweep.rs linalg/matfree.rs | 33, 8 | Time integration, Newton and optimisation, matrix exponential action, reduced-order frequency sweeps and matrix-free high-order operators.
fem         | Finite-element basics | fem.rs edge.rs fields.rs pml.rs periodic.rs props.rs | 9, 4, 45, 69 | Lagrange P1/P2 elements and quadrature, the 2D edge-element basis, derived fields, PML coordinate stretching, Floquet pairs and material properties per domain.
material    | Material specification | material/ | 29, 30, 34, 35 | TOML material files with explicit units, tensors, dispersion models evaluated at the wavelength, provenance and the library loader.
optical     | Optical mode sets | optical.rs | 21, 9, 45 | Optical modes of a 2D cross-section as a field-evaluable set with the material data to evaluate them.
section     | Frame and sections | section.rs | 1 | The one place that knows the wafer axes (x lateral, z normal, waves along +y).
expr        | Expressions | expr.rs | 65, 2 | A small safe expression evaluator for user input; a copy of geom's kept equal by a corpus test.
runtime     | Runtime services | par.rs progress.rs validate.rs lib.rs | 9, 11 | The shared worker pool, stage log and cancel, and refusal of NaN and infinity.
view femcore | meshing | Meshing and adaptation | mesh meshio fixtures meshing meshing3d structpml meshcheck metric adaptops adaptdrive estimators
view femcore | linalg | Linear algebra | lacore ladriver ladirect laksp laeigen laprecond laother
view femcore | fem | Finite elements, materials and services | fem material optical section expr runtime mesh
