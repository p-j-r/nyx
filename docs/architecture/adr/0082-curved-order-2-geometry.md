# 82. Curved (isoparametric, quadratic) tetrahedra for the `ems` cavity study; amends ADR 18

Date: 2026-11-15

## Status

Implemented on `capability-r12-curved-elements` (default off; paul decides later whether it becomes the default). Amends ADR 18, which deferred curved nodes.

## Context

The 3D lattice mesher (ADR 40) fits boundary vertices to the exact surface but the cells are straight, so a sphere is a polyhedron and the geometric error is O(h^2) whatever the element order (ADR 18). The A7 simple-cubic cross-check against MPB (`benchmarks/johnson-joannopoulos-2001`, part 2) failed its pre-declared limits by exactly that error: it falls as h^2 and its sign flips between a dielectric and an air sphere.

## Decision

* **Storage.** `CellMesh<D>::edge_nodes: Vec<([u32; 2], [f64; D])>`: for an edge (sorted vertex pair) the position of its mid-edge geometry node, sorted by pair, empty by default. An edge not listed is straight (its node is the midpoint). The solver mesh (`mesh::SimplexMesh`) carries the same list by edge index (`edge_node`, `cell_is_curved`, `cell_geometry_nodes`, `edge_node_list`). Round trips (`Mesh::from_solver_mesh`, `discretise`, `from_cells`) keep it; `from_cells` drops it when cells had to be split; refinement and adaptation start from lists that are empty. HDF5 mesh files carry two optional datasets (`/mesh/edge_node_pairs`, `/mesh/edge_node_points`); a file without them is a straight mesh. XDMF output keeps straight cells.
* **Placing the nodes** (`femcore::meshing::curve3`). For every edge of a boundary or interface facet `Domain3::project_edge` gives the surface position of the midpoint; it is kept when it moved by more than `1e-9` and less than `0.25` of the edge length (a planar face adds nothing). Every curved cell must keep `det J > 0.1 det J_straight` at its ten nodes and at the points of a degree-6 rule; if not, the offsets of its curved edges are halved (every cell around an edge sees the same offset), and after four halvings the edges are straightened. `CurveReport` counts curved edges and cells, relaxed and straightened edges and the smallest scaled Jacobian. An edge whose ends lie on two faces that place its midpoint differently (a crease) stays straight (`Domain3::edge_on_crease`). Independent of the thread count.
* **Switch.** `MeshSettings.curved`, `MeshRequest.curved`, `GeomMeshOpts.curved`, `set_mesh` argument `curved`; default false and not serialised when false, so every existing project, golden and example is byte-identical. After an adaptation or structured layers the edges are curved last (`curve_solver_mesh`).
* **Kernel.** `femcore::fem::CurvedTet` (ten nodes, quadratic Lagrange map, `jacobian(lam) -> (J^-T grad_ref lambda_a, det J)`), `quad_rule_collapsed(degree)` (Gauss on a collapsed cube, any degree). `nyx-ems::cavity::assemble_maxwell` integrates cells with a curved edge with a degree-6 rule and, at each point, the gradients `g_a(q)` of the barycentric functions in place of the constant ones: `basis3` is unchanged, because `lam_a g_b - lam_b g_a` is then exactly the covariant Piola image of the Whitney function and `2 g_a x g_b` equals `J curl_ref / det J`; the order-2 gradient and face functions follow by the product rule. Cells with no curved edge keep the affine path and the same arithmetic (the matrices of a straight mesh are unchanged). The discrete gradient, numbering, periodic relations and reduction are topological and unchanged. `emw_eigen_fields` evaluates `E` with the same per-point gradients.
* **Other physics.** `sms` and `dfs` ignore `edge_nodes` (they see straight cells); only the `ems` cavity study uses them.

## Consequences

* Closed forms (tests on the bookmark): the volume of a meshed unit sphere is 2.1e-5 off at h = r/4 (straight 1.3e-2) and falls at about h^4.5 against h^2 for straight; the lowest triplet of a PEC sphere cavity converges at rate 2.96 (curved) against 2.10 (straight) over h = r/5 to r/7, and at h = r/3.5 the error is 3.5e-5 against 6.4e-3. A constant field is represented exactly on a quadratic image of a box, and the curl of the basis equals the curl of its field by finite differences.
* Rough spots (`ROUGH_SPOTS.md` R32 to R35): plotted quadratic nodes stay at the straight positions; sms and dfs see straight cells; creases and a curved face meeting a periodic face stay straight (`capability-r7-periodic-cells` must give matching curved edges on a paired face); curved cells inside a PML are allowed but untested; the volume integrals of `nyx-analysis` (mode volume, layer fraction) use straight measures.
* Contracts: `mesh_cells` 1.1, `mesh_request` 1.1, `mesh_report` 1.1, `engine_commands` 1.8 (1.6 is A3, 1.7 is R11), `mcp_tools` 2.4 (2.3 is R11); `results_file` additive (optional datasets).

## Code

`crates/femcore/src/fem.rs` (`CurvedTet`, `quad_rule_collapsed`), `mesh/cells.rs`, `mesh/mod.rs`, `mesh/exchange.rs`, `meshing/curve3.rs`, `meshing/geomdrive3.rs`, `crates/nyx-ems/src/cavity.rs`, `crates/nyx-engine/src/h5io.rs`.
