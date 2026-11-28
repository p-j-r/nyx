# 135. Per-physics meshing, active-domain subsets and domain-aware cross-mesh field transfer

Date: 2026-11-27

## Status

Accepted. Extends ADR 9 (three-stage kernels), ADR 11 (deterministic parallel reductions), ADR 25 (the `Physics` trait and `Ctx`), ADR 54 (the
mesh request), ADR 61 (the solve chain), ADR 82 (curved cells), ADR 90 (whole-array reductions), ADR 128 (native mixed cells) and ADR 131
(no library panics, finite guards).

## Context

Until now one project had one mesh: every physics of a one-way chain was solved on it, and a field exported by one physics was handed to the
next one as the raw per-dof array of that shared mesh (`Ctx::inputs`). Coupled problems do not want one mesh. An optical mode wants a fine
mesh in the core and the cladding; a thermal solve is smooth and wants a coarse one; a solid-mechanics solve wants the solid regions and its
elastic absorbing layers only, never the air. The reductions (`integrate`, `field_op`, `overlap`) refused outright to combine fields solved on
different meshes.

Transferring a field between two meshes of one geometry has one trap that matters more than accuracy: a field that is discontinuous across a
material interface (a stress, a flux, a permittivity-weighted quantity, a temperature-dependent material value) is smeared when the target
point near the interface is read in a source cell of the neighbouring material. Two non-matching discretisations of one interface disagree
by up to the chord deviation, so a naive point location reads the wrong side for a band of target nodes.

## Decision

* **Per-physics mesh settings.** `Project::mesh_by_physics: BTreeMap<String, MeshSettings>` (serde default, skipped when empty, so every
  existing project file is unchanged) gives a physics settings of its own; `Project::mesh_for(physics)` falls back to `Project::mesh`.
  `MeshSettings::active_domains: Option<Vec<String>>` names the regions to keep (names or `name*` patterns of the geometry). `FemcoreMesher`
  meshes the whole model and `MeshedModel::restrict_to_active` drops the other cells (`femcore::transfer::restrict_to_domains`): vertices are
  renumbered in their old order, labelled facets that are faces of kept cells keep their labels, curved edges and quad-face nodes are carried
  over, and **domain labels and the `LabelMap` are kept**, so every mesh of the project shares the canonical `geom` region and boundary names.
* **Solve chain.** `solve::run_chain_models` meshes a physics with settings of its own (cached by the hash of canonical geometry, settings
  and sizing, inside the chain and process-wide through weak references) and solves it on that mesh with a view of the project whose `mesh`
  is its own (the checks and the AMR loops read it). `JobOutput::Solved::models` carries those meshes and the engine stores each result on the
  mesh it was solved on, so a reduction reads it on exactly its own cells. A project with per-physics meshes is solved in-process (a worker
  process returns one mesh). The problem key of a physics with its own mesh includes its settings.
* **Imports across meshes.** `nyx_physics::TransferredField` (the source `Arc<Mesh>`, its discretisation, canonical domain and boundary names,
  order, components, per-dof values and cell fields; `is_same_mesh` as the fast-path test) and the additive `Ctx::transferred_inputs`. When the
  exporting physics was solved on the same mesh (`Arc` identity or equal vertices, cells and labels) `Ctx::inputs` is filled exactly as
  before, at zero cost. Otherwise the field is sampled at the dofs of the importing mesh, domain-filtered, and put into `Ctx::inputs` (so a
  kernel that only reads `inputs` works unchanged) and the source itself into `Ctx::transferred_inputs`. A target node whose region the
  source mesh does not have is an error naming the count, never a silent default.
* **`femcore::transfer`.** One cell-AABB BVH per canonical domain name and one over all cells, over `SimplexMesh` (`Tri`, `Tet`, curved
  `CurvedTet`) and `MixedDiscretisation` (`Tri`, `Quad`, `Tet`, `Prism`, `Hex`, straight or curved). Three rules:
  1. *Strict domain-label filtering:* a query carries the canonical name of its target cell's region and is evaluated only in source cells of
     the same name. A field that jumps 12:1 across a silicon / oxide interface keeps both values at the interface.
  2. *Direct sampling and conservative projection:* `eval_at_points`, `eval_scalar`, `eval_vector`, `sample_cell_nodes`, `sample_dofs` and
     `sample_quadrature` evaluate `P1`/`P2` nodal fields, per-cell `CellFields`, per-cell constants and mixed Lagrange vectors directly.
     `project_scalar_conservative` solves the domain-wise `M_tgt Q_tgt = int N_tgt Q_src dV` (consistent mass by Jacobi-PCG, or HRZ lumped)
     on an independent continuous space per domain, with the right-hand side integrated on the target cells from the directly sampled
     source, then a constant correction per domain makes `int Q_tgt dV` equal the source integral (taken by the source quadrature) to
     rounding.
  3. *Boundary-facet transfer:* `BoundaryIndex` keeps one facet BVH per canonical boundary name and evaluates at the closest point of the
     matching facets, in the adjacent cell of a chosen side for an interface.
  Point location tries the BVH candidates in ascending cell number (the first containing cell wins); curved and mixed cells are inverted by
  Newton iteration; a point outside its domain is snapped to the nearest cell of the same domain within `snap_rel` times the source
  diagonal (branch-and-bound, ties to the lower cell number). Queries run in fixed chunks on the shared pool and are concatenated in chunk
  order; the projection scatters sequentially, so results do not depend on the thread count.
* **Reductions.** `reduce::prepare` no longer refuses leaves of results on different meshes: the first leaf fixes the primary mesh, and a
  leaf of another mesh is sampled directly at the degree-4 quadrature points of every primary cell (for `integrate`, `overlap` and
  `confinement`, preserving the exact polynomial order of the source even on a coarser or lower-order primary mesh) and at the local nodes
  of every primary cell (for `extremum`, `probe`, `line_cut` and boundary integrals), domain-filtered; per-cell values are sampled at
  primary cell centres, region indicators are evaluated on the primary mesh. The answers of `field_op`, `integrate`, `overlap` and the other
  reductions carry `cross_mesh_transfer: true` and the snapped and missing point counts. A primary point whose region the other mesh lacks
  (an active-domain subset) reads 0.

## Consequences

* Existing projects, results files, contract exemplars and kernels are unchanged: the new project and mesh-settings fields are skipped when
  empty, `Ctx::new` initialises `transferred_inputs` empty, and the same-mesh import path is the old code path.
* The crate layering holds: `femcore::transfer` depends on nothing outside femcore, `nyx-physics` uses it through femcore, and only
  `nyx-engine` drives meshing and transfer.
* Volume reductions (`integrate`, `overlap`, `confinement`) evaluate transferred fields directly at the degree-4 quadrature points of the
  primary mesh, so a `P2` source field is integrated exactly on every primary cell contained in a source cell even when the primary mesh is
  `P1`. Point/extremum/boundary queries read the primary cell-node interpolant of the transferred field.
* A worker-isolated or shared-queue solve of a project with per-physics meshes is refused when isolation is asked for explicitly, and runs
  in-process otherwise.

## Code

`crates/femcore/src/{transfer,lib}.rs`, `crates/nyx-physics/src/{physics,lib}.rs`, `crates/nyx-engine/src/{project,meshing,solve,reduce,engine,jobs,sweep,amr}.rs`,
`crates/nyx-engine/tests/cross_mesh.rs`.
