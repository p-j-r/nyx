# Meshing facade: code

C4 level 4 (code) of the component **Meshing facade** (`femcore.meshing`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> MeshPlan and mesh_model (the one way to mesh a geom model), the serialisable request and report, sizing from the physics, and the 2D triangulation core.

**Decisions:** [ADR 6](../../../adr/0006-one-meshing-interface.md), [ADR 17](../../../adr/0017-physics-driven-mesh-sizing.md), [ADR 54](../../../adr/0054-mesh-request-and-report-as-data.md)  
**Size:** 9 files, 4475 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/meshing/geomdrive.rs` | 363 | From a `geom` model to a mesh and to the integers the solvers take, in the wafer frame throughout: x lateral, y the direction of propagation, z the surface normal. A 2D model is the x-z cross-section and its mesh has the columns `(x, z)`... |
| `crates/femcore/src/meshing/mesher.rs` | 1352 | The triangulation core shared by the geom path (`geomdrive`): chains of boundary points become a constrained Delaunay triangulation with Ruppert refinement (via the `spade` crate), cells labelled by a domain function; and the sweep of a ... |
| `crates/femcore/src/meshing/mod.rs` | 595 | The one way to mesh a `geom` model: a `MeshPlan` says how fine and with which mesher, `mesh_model` returns a `MeshedGeom` with the names, labels and exact geometry of the model, whichever mesher made the cells. |
| `crates/femcore/src/meshing/physics_size.rs` | 221 | Element sizes from the physics: how fine a mesh must be for the waves the study looks for, before any solution exists. The core takes the neutral `meshing::Sizing` a physics reports (`cell_sizes_for`, `vertex_sizes_for`, `metrics_for`, `... |
| `crates/femcore/src/meshing/request.rs` | 417 | A mesh asked for and reported in data (ADR 54): the serialisable request that the engine, the MCP tools and the window can build, run with [`mesh_request`], and the report that comes back, with JSON Schemas for both (draft 2020-12). Noth... |
| `crates/femcore/src/meshing/sheet3d_amr.rs` | 361 |  |
| `crates/femcore/src/meshing/sheet3d_curve.rs` | 433 |  |
| `crates/femcore/src/meshing/singularity.rs` | 640 |  |
| `crates/femcore/src/meshing/sizing.rs` | 93 | The neutral request a physics makes of the mesher: what to resolve, with no knowledge of any physics or of `ProblemSpec`. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `geom` (geom)](../../3-components/geom.md) | complex, curve, diff, eval, labels, math, model, query +2 | Aabb, ChangeKind, Complex, Curve, DiscretizeOpts, EdgeGeom, EntityKind, Evaluated, Evaluator, FaceGeom +10 |
| [Adaptation drivers](./adaptdrive.md) | adapt::drive, adapt::drive3, adapt::geom | AdaptOptions, AdaptReport, Geometry2, PolylineGeometry, adapt2_with, adapt3 |
| [Adaptation operators](./adaptops.md) | adapt::ops3, adapt::topo, adapt::topo3 | AdaptOptions3, Sig, TopoOptions3, Work2 |
| [Finite-element basics](./fem.md) | props | Pml, PmlDir |
| [Material specification](./material.md) | material | M6 |
| [Mesh](./mesh.md) | mesh | SimplexMesh |
| [3D meshing](./meshing3d.md) | meshing::adapt3d, meshing::curve3, meshing::geomdrive3, meshing::mesher3d, meshing::sheet3d | CurveReport, Domain3, EdgeNodes, GeomDomain3, Sizes3, adapt_to_metrics3, apply_sheet_post_enrichment, build_sized_overlay_mesh_2d, extrude_overlay_2d, extrude_overlay_2d_with_plan +3 |
| [Metric fields](./metric.md) | adapt, adapt::field, adapt::field3, adapt::metric, adapt::work | BackgroundMetric2, BackgroundMetric3, MetricField, SimplexMesh, Sym, bbox_of, intersect, isotropic, simplex_from_adapt, simplex_to_adapt +1 |
| [Runtime services](./runtime.md) | progress | note |
| [Structured absorbing layers](./structpml.md) | meshing::pml_geom | PmlMode, PmlReport, structured_pml |

## Used by

| Component | Modules | Items |
|---|---|---|
| [dfs: Electrostatics kernel](../dfs/electro.md) | meshing::singularity | DEFAULT_QUARTER_POINT_XI, apply_quarter_point_nodes, sheet_rim_vertices |
| [ems: Physics adapter](../ems/adapter.md) | meshing::physics_size | optical_wavelength, speed_of_light |
| [engine: Meshing seam](../engine/meshing.md) | meshing, meshing::geomdrive, meshing::mesher, meshing::physics_size, meshing::request, meshing::sheet3d_amr | Adaptive, Layer, MeshReport, MeshRequest, MeshedGeom, MesherRequest, PhysicsSizing, RemappedLayer +12 |
| [engine: Registry and facade](../engine/registry.md) | meshing::geomdrive |  |
| [Error estimators](./estimators.md) | meshing, meshing::geomdrive | Adaptive, GeomMeshOpts, MeshedGeom, adapt_to_field, adapt_to_metrics, evaluate, geometry_floor, mesh_geom_2d |
| [Fixtures and test kit](./fixtures.md) | meshing, meshing::geomdrive, meshing::mesher | GeomMeshOpts, Layer, MeshPlan, MeshedGeom, block, evaluate, extrude_y, mesh_evaluated +3 |
| [3D meshing](./meshing3d.md) | meshing, meshing::geomdrive, meshing::mesher, meshing::physics_size, meshing::sheet3d_curve, meshing::singularity +1 | Adaptive, Chain, DEFAULT_QUARTER_POINT_XI, Fill, GeomMeshOpts, MeshParams, MeshPlan, MeshedGeom +16 |
| [Structured absorbing layers](./structpml.md) | meshing::geomdrive | MeshedGeom |
| [physics: Descriptors and schemas](../physics/descriptors.md) | meshing::sizing | Sizing |
| [sms: Physics adapter](../sms/adapter.md) | meshing::physics_size | elastic_wavelength |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace meshing_geomdrive_rs {
    class GeomMeshOpts {
        <<struct>>
        +h f64
        +max_dev Option~f64~
        +min_angle_deg f64
        +curved bool
        +gap_cells Option~f64~
        +z_coarsen bool
        +max_cells Option~usize~
        +dz_min Option~f64~
        +more_7_fields
        +new(h) Self
    }
    class MeshedGeom {
        <<struct>>
        +mesh SimplexMesh
        +labels LabelMap
        +eval Arc~Evaluated~
        +warnings Vec~String~
        +domain(region) Option~i32~
        +domains(pattern) Vec~i32~
        +boundary(name) Vec~i32~
        +select(selector) Result
        +pml_materials(strength) BTreeMap
    }
    class SizeMap {
        <<struct>>
        +region_h Vec~f64~
        +edge_rule Vec
    }
    class Plate {
        <<struct>>
        +mesh SimplexMesh
        +domain i32
        +xmin Vec~i32~
        +xmax Vec~i32~
        +zmin Vec~i32~
        +zmax Vec~i32~
    }
    class Block {
        <<struct>>
        +mesh SimplexMesh
        +domain i32
        +xmin Vec~i32~
        +xmax Vec~i32~
        +zmin Vec~i32~
        +zmax Vec~i32~
        +ymin Vec~i32~
        +ymax Vec~i32~
    }
    class fns_meshing_geomdrive_rs {
        <<module meshing/geomdrive.rs>>
        +evaluate(list) Result
        ~size_map(ev, h) SizeMap
        +mesh_geom_2d(ev, opts) Result
        +select_labels(ev, labels, selector) Result
        +plate(lx, lz, h) Result
        +block(lx, ly, lz, h, ny) Result
        +sort_nodes(m) SimplexMesh
    }
  }
  namespace meshing_mesher_rs {
    class Layer {
        <<struct>>
        +thickness f64
        +n usize
        +domain_offset i32
        +fill Option~i32~
    }
    class Chain {
        <<struct>>
        +pts Vec
        +label i32
        +part usize
    }
    class Fill {
        <<struct>>
        +lo slice
        +hi slice
        +h f64
        +owner BoxdynFnf642bool
    }
    class MeshParams {
        <<struct>>
        +h f64
        +h_fine_min f64
        +min_angle_deg f64
        +lo slice
        +hi slice
        +max_steiner_multiplier Option~usize~
        +strict_mesh bool
        +max_cells Option~usize~
    }
    class MeshChainsReport {
        <<struct>>
        +n_init usize
        +added_vertices usize
        +refinement_complete bool
        +healed_caps usize
        +dropped_slivers usize
        +warnings Vec~String~
    }
    class fns_meshing_mesher_rs {
        <<module meshing/mesher.rs>>
        ~dist_point_segment(p, a, b) f64
        ~mesh_chains(chains, fills, mp, domain_at) Result
        +extrude_y(m2, layers, y0) Result
        +extrude_y_with_remap(m2, layers, y0) Result
        +sweep_model_3d(ev, labels, h_global, max_dev, min_angle_deg, sizes, layers) Result
        +RemappedLayer type
    }
  }
  namespace meshing_mod_rs {
    class Adaptive {
        <<struct>>
        +scale f64
        +grade f64
        +max_cells Option~usize~
    }
    class Mesher {
        <<enumeration>>
        Simple
        Adaptive
    }
    class MeshPlan {
        <<struct>>
        +opts GeomMeshOpts
        +mesher Mesher
        +pml PmlMode
        +simple(h) Self
        +adaptive(h, scale) Self
    }
    class GeomOracle2 {
        <<struct>>
        +from_work(ev, labels, w) Self
    }
    class Reuse {
        <<enumeration>>
        Hit
        Meshed
    }
    class CacheStats {
        <<struct>>
        +requests usize
        +hits usize
        +eval_seconds f64
        +mesh_seconds f64
    }
    class MeshCache {
        <<struct>>
        +stats CacheStats
        +new(capacity) Self
        +mesh(list, plan) Result
        +mesh_extruded(list, plan, layers) Result
    }
    class fns_meshing_mod_rs {
        <<module meshing/mod.rs>>
        +mesh_model(list, plan) Result
        +mesh_evaluated(ev, plan) Result
        +mesh_evaluated_report(ev, plan) Result
        +mesh_for_sizing(list, plan, sizing, sz) Result
        +wave_sizes_by_region(ev, labels, sizing, sz) Result
        +a_priori_metrics(g, opts, a) Result
        +geometry_floor(g, opts, a) Result
        +adapt_to_metrics(g, metrics, a) Result
        +adapt_to_field(g, field, a) Result
        ~periodic_labels(g) Vec~i32~
    }
  }
  namespace meshing_physics_size_rs {
    class PhysicsSizing {
        <<struct>>
        +points_per_wavelength Option~f64~
        +pml_relax f64
        +eigen_headroom f64
        +ppw(order) f64
    }
    class fns_meshing_physics_size_rs {
        <<module meshing/physics_size.rs>>
        +wave_speed_range(c, rho) tuple
        +cell_sizes_for(sizing, resolve, cell_labels, pml, sz) Result
        +vertex_sizes_for(sizing, resolve, mesh, pml, sz) Result
        +metrics_for(sizing, resolve, mesh, pml, sz, ceiling) Result
        +adapt_tets_for(mesh, sizing, resolve, pml, sz, ceiling, grade, opt) Result
        +optical_wavelength(eps, mu, lam0) f64
        +speed_of_light() f64
        +elastic_wavelength(c, rho, f) f64
    }
  }
  namespace meshing_request_rs {
    class MesherRequest {
        <<enumeration>>
        Simple
        Adaptive
    }
    class WavesRequest {
        <<struct>>
        +sizing Sizing
        +points_per_wavelength Option~f64~
        +pml_relax f64
        +eigen_headroom f64
    }
    class MeshRequest {
        <<struct>>
        +h f64
        +max_dev Option~f64~
        +mesher MesherRequest
        +structured_pml bool
        +waves Option~WavesRequest~
        +curved bool
        +gap_cells Option~f64~
        +z_coarsen bool
        +more_9_fields
        +plan() MeshPlan
    }
    class DomainReport {
        <<struct>>
        +id i32
        +name String
        +cells usize
        +wave_size Option~f64~
    }
    class AdaptSummary {
        <<struct>>
        +iterations usize
        +converged bool
        +capped bool
        +q_min f64
        +cells_start usize
        +cells_end usize
    }
    class MeshReport {
        <<struct>>
        +dim usize
        +vertices usize
        +cells usize
        +facets usize
        +domains Vec~DomainReport~
        +adaptation Option~AdaptSummary~
        +structured_layers Option
        +curved Option~CurvedSummary~
        +more_1_fields
    }
    class CurvedSummary {
        <<struct>>
        +edges usize
        +cells usize
        +min_scaled_jacobian f64
        +from_mesh(m) Option~Self~
    }
    class fns_meshing_request_rs {
        <<module meshing/request.rs>>
        +mesh_request(list, req) Result
        +request_schema() Value
        +report_schema() Value
    }
  }
  namespace meshing_sheet3d_amr_rs {
    class OverlayOracle2 {
        <<struct>>
        +ov Complex
        +from_work(ov, w) Self
        +pinned_vertices(ov, m2) Vec~u32~
    }
    class OverlaySheetFeatures {
        <<struct>>
        +rim_segments Vec
        +corners Vec
        +plane_axes slice
    }
    class fns_meshing_sheet3d_amr_rs {
        <<module meshing/sheet3d_amr.rs>>
        +extract_overlay_sheet_features(cx) Option~OverlaySheetFeatures~
        +adapt_overlay_2d(ov, m2, metrics, a) Result
        +adapt_overlay_2d_field(ov, m2, field, a) Result
        +build_overlay_mesh_2d(ev, plan) Result
        +sheet_rim_vertices_3d(mesh, sheet_labels) Vec~bool~
        +adapt_sheet_full_3d(g, metrics, a) Result
        +adapt_sheet_full_3d_with_plan(g, metrics, a, plan, waves) Result
    }
  }
  namespace meshing_sheet3d_curve_rs {
    class Sheet3dTransitionRef {
        <<struct>>
        +map_v slice
        +rank slice
    }
    class fns_meshing_sheet3d_curve_rs {
        <<module meshing/sheet3d_curve.rs>>
        +relax_edge_nodes(pos, cells, cand) tuple
        +project_sheet3d_curved_edges(cx, m2, bnd_edges2, z_planes, used_v, slab_layer, level_of_zp, transitions, mesh) CurveReport
        +apply_sheet3d_curved_edges(g) CurveReport
    }
  }
  namespace meshing_singularity_rs {
    class SheetCorner {
        <<struct>>
        +vertex usize
        +pos slice
        +omega_rad f64
        +alpha f64
        +h_edge f64
        +is_sharp() bool
        +h_at(p) f64
        +h_min() f64
    }
    class QuarterPointReport {
        <<struct>>
        +rim_vertices usize
        +shifted_edges usize
        +qp_cells usize
        +clamped_curved_edges usize
        +unshifted_edges usize
        +min_interior_jacobian_ratio f64
    }
    class fns_meshing_singularity_rs {
        <<module meshing/singularity.rs>>
        +morrison_lewis_exponent(omega_rad, eps_ratio) f64
        +corner_singularity_exponent(omega_rad, eps_ratio) f64
        +corner_grading_h(d, h_edge, alpha) f64
        +babuska_corner_h(d, h_edge, alpha) f64
        +corner_graded_h_at(p, corners, base_h) f64
        +coupled_corner_dz(base_dz, corners) f64
        +measure_sheet_corner_angles_m2(m2, sheet_covered_cells) Vec
        +measure_overlay_sheet_corners(cx, z_s) Vec
        +collect_overlay_sheet_corners(cx, z_s, ov_edge_h, eps_ratio) Vec~SheetCorner~
        +analytical_qp_jacobian_ratio(n_rim_in_tet, lam_out, xi_qp) f64
        +more_2_functions()
        +DEFAULT_QUARTER_POINT_XI const
        +CORNER_GRADING_P const
        +CORNER_RADIUS_FACTOR const
        +CORNER_H_MIN_FACTOR const
        +more_1_items
    }
  }
  namespace meshing_sizing_rs {
    class Sizing {
        <<struct>>
        +order u8
        +headroom f64
        +wavelengths BTreeMap
        +by_label(order, headroom, wavelengths) Sizing
        +expand_wildcards(names) Sizing
        +resolve(domain) Result
    }
    class fns_meshing_sizing_rs {
        <<module meshing/sizing.rs>>
        +label_names(name) Option~i32~
    }
  }
  Mesher --> Adaptive
  MeshPlan --> GeomMeshOpts
  MeshPlan --> Mesher
  MeshCache --> CacheStats
  WavesRequest --> Sizing
  MeshRequest --> MesherRequest
  MeshRequest --> WavesRequest
  MeshReport --> DomainReport
  MeshReport --> AdaptSummary
  MeshReport --> CurvedSummary
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `GeomMeshOpts` | `meshing/geomdrive.rs:29` | fields: `h`, `max_dev`, `min_angle_deg`, `curved`, `gap_cells`, `z_coarsen`, `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, `max_steiner_multiplier`, `strict_mesh`, ... —  Sizes of a mesh made from a model. |
| struct | `MeshedGeom` | `meshing/geomdrive.rs:85` | fields: `mesh`, `labels`, `eval`, `warnings` —  A meshed model: the mesh, the compiled names, and the evaluation they came from. |
| fn | `evaluate` | `meshing/geomdrive.rs:157` | `pub fn evaluate(list: &FeatureList) -> Result<Evaluated, String>` —  Evaluates a feature list to the tier the mesher needs. |
| struct | `SizeMap` | `meshing/geomdrive.rs:163` | fields: `region_h`, `edge_rule` —  The element sizes a model asks for: per region (the global size, the `Size` rules on regions, and thickness / layers in an absorbing |
| fn | `size_map` | `meshing/geomdrive.rs:168` | `pub(crate) fn size_map(ev: &Evaluated, h: f64) -> SizeMap` |
| fn | `mesh_geom_2d` | `meshing/geomdrive.rs:194` | `pub fn mesh_geom_2d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String>` —  Meshes a 2D model (the x-z cross-section) with triangles. |
| fn | `select_labels` | `meshing/geomdrive.rs:281` | `pub fn select_labels(ev: &Evaluated, labels: &LabelMap, selector: &str) -> Result<Vec<i32>, String>` —  Labels of whatever a selector text picks in an evaluated model: boundary labels for edges, domain labels for regions. |
| struct | `Plate` | `meshing/geomdrive.rs:298` | fields: `mesh`, `domain`, `xmin`, `xmax`, `zmin`, `zmax` —  A rectangle `[0, lx] x [0, lz]` of the x-z plane as one domain, meshed with element size `h`, with the boundary labels of its four sides |
| fn | `plate` | `meshing/geomdrive.rs:307` | `pub fn plate(lx: f64, lz: f64, h: f64) -> Result<Plate, String>` |
| struct | `Block` | `meshing/geomdrive.rs:320` | fields: `mesh`, `domain`, `xmin`, `xmax`, `zmin`, `zmax`, `ymin`, `ymax` —  A box `[0, lx] x [0, ly] x [0, lz]`: the plate swept along y through `ny` layers (`extrude_y`). The vertices of the section are put in coordinate order first (`sort_nodes`), so that opposite side walls can be paired by a Floquet condition. The side walls keep the labels of the |
| fn | `block` | `meshing/geomdrive.rs:331` | `pub fn block(lx: f64, ly: f64, lz: f64, h: f64, ny: usize) -> Result<Block, String>` |
| fn | `sort_nodes` | `meshing/geomdrive.rs:342` | `pub fn sort_nodes(m: &SimplexMesh) -> SimplexMesh` —  The mesh with its vertices renumbered in the order of their coordinates (x first, then the second axis, equal to a tolerance of |
| fn | `GeomMeshOpts::new` | `meshing/geomdrive.rs:62` | `pub fn new(h: f64) -> Self` |
| fn | `MeshedGeom::domain` | `meshing/geomdrive.rs:96` | `pub fn domain(&self, region: &str) -> Option<i32>` —  The domain label of a region. |
| fn | `MeshedGeom::domains` | `meshing/geomdrive.rs:102` | `pub fn domains(&self, pattern: &str) -> Vec<i32>` —  Domain labels of regions by name. A list is comma separated; `"sapphire*"` is `sapphire` and every region whose name starts with |
| fn | `MeshedGeom::boundary` | `meshing/geomdrive.rs:122` | `pub fn boundary(&self, name: &str) -> Vec<i32>` —  Boundary labels by name (`"rib.top"`, `"air\|rib"`; a face split into pieces has several) or tag/group. |
| fn | `MeshedGeom::select` | `meshing/geomdrive.rs:127` | `pub fn select(&self, selector: &str) -> Result<Vec<i32>, String>` —  Boundary labels of whatever a selector text picks (`"outer & plane(z = -3e-6)"`). |
| fn | `MeshedGeom::pml_materials` | `meshing/geomdrive.rs:133` | `pub fn pml_materials(&self, strength: f64) -> BTreeMap<i32, Pml>` —  The absorbing layers of the model as solver stretching, keyed by domain: `strength` is the imaginary stretch at the outer edge, the |
| struct | `Layer` | `meshing/mesher.rs:11` | fields: `thickness`, `n`, `domain_offset`, `fill` |
| struct | `Chain` | `meshing/mesher.rs:23` | fields: `pts`, `label`, `part` —  One closed polyline edge chain: sample points and the boundary label of every segment. |
| fn | `dist_point_segment` | `meshing/mesher.rs:30` | `pub(crate) fn dist_point_segment(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64` |
| struct | `Fill` | `meshing/mesher.rs:39` | fields: `lo`, `hi`, `h`, `owner` —  One lattice fill: a region that wants elements of size `h` (finer than the global size) inside the box `lo..hi`; `owner` says whether a point is |
| struct | `MeshParams` | `meshing/mesher.rs:46` | fields: `h`, `h_fine_min`, `min_angle_deg`, `lo`, `hi`, `max_steiner_multiplier`, `strict_mesh`, `max_cells` |
| struct | `MeshChainsReport` | `meshing/mesher.rs:64` | fields: `n_init`, `added_vertices`, `refinement_complete`, `healed_caps`, `dropped_slivers`, `warnings` |
| fn | `mesh_chains` | `meshing/mesher.rs:361` | `pub(crate) fn mesh_chains( chains: &[Chain], fills: &[Fill], mp: &MeshParams, domain_at: &dyn Fn([f64; 2]) -> Option<i32>, ) -> Result<(SimplexMesh, MeshChainsReport), String>` —  The shared tail of the meshers: constrained Delaunay triangulation of the chains (conforming, labelled), lattice points, Ruppert refinement, |
| type | `RemappedLayer` | `meshing/mesher.rs:690` |  A sweep [`Layer`] paired with an optional domain-label remap `(source_domain -> target_domain)` for that layer. |
| fn | `extrude_y` | `meshing/mesher.rs:696` | `pub fn extrude_y(m2: &SimplexMesh, layers: &[Layer], y0: f64) -> Result<SimplexMesh, String>` —  Sweeps the x-z cross-section of the wafer frame (a 2D mesh with the columns `(x, z)`) along the propagation axis y through the given |
| fn | `extrude_y_with_remap` | `meshing/mesher.rs:703` | `pub fn extrude_y_with_remap( m2: &SimplexMesh, layers: &[RemappedLayer], y0: f64, ) -> Result<SimplexMesh, String>` —  Like [`extrude_y`], with an optional per-layer domain map `(source_domain -> target_domain)` applied when `l.fill` is `None`. |
| fn | `sweep_model_3d` | `meshing/mesher.rs:782` | `pub fn sweep_model_3d( ev: &geom::eval::Evaluated, labels: &geom::labels::LabelMap, h_global: f64, max_dev: Option<f64>, min_angle_deg: f64, sizes: &crate::meshing::geomdrive3::Sizes3, layers: &[RemappedLayer], ) -> Result<SimplexMesh, String>` —  Sweeps a 3D `geom` model along the propagation axis `y` through `layers`: extracts the combined 2D `(x, z)` cross-section |
| struct | `Adaptive` | `meshing/mod.rs:55` | fields: `scale`, `grade`, `max_cells` —  Settings of the adaptive mesher. |
| enum | `Mesher` | `meshing/mod.rs:73` | variants: `Simple`, `Adaptive` —  Which mesher makes the cells. |
| struct | `MeshPlan` | `meshing/mod.rs:80` | fields: `opts`, `mesher`, `pml` —  How a model is meshed: the sizes (as for `mesh_geom_2d`) and the mesher. |
| fn | `mesh_model` | `meshing/mod.rs:106` | `pub fn mesh_model(list: &FeatureList, plan: &MeshPlan) -> Result<MeshedGeom, String>` —  Meshes a model. |
| fn | `mesh_evaluated` | `meshing/mod.rs:111` | `pub fn mesh_evaluated(ev: Evaluated, plan: &MeshPlan) -> Result<MeshedGeom, String>` —  Meshes an evaluated model. |
| fn | `mesh_evaluated_report` | `meshing/mod.rs:116` | `pub fn mesh_evaluated_report(ev: Evaluated, plan: &MeshPlan) -> Result<(MeshedGeom, Option<AdaptReport>), String>` —  Like `mesh_evaluated`, with the report of the adaptation when there was one. |
| fn | `mesh_for_sizing` | `meshing/mod.rs:140` | `pub fn mesh_for_sizing(list: &FeatureList, plan: &MeshPlan, sizing: &Sizing, sz: &PhysicsSizing) -> Result<(MeshedGeom, Option<AdaptReport>), String>` —  Meshes a model of either dimension for the neutral `Sizing` a physics reports (ADR 50): the sizes the model asks for, and in |
| fn | `wave_sizes_by_region` | `meshing/mod.rs:205` | `pub fn wave_sizes_by_region(ev: &Evaluated, labels: &LabelMap, sizing: &Sizing, sz: &PhysicsSizing) -> Result<Vec<f64>, String>` —  The element size the waves of `sizing` need in every region of a 3D model, by region index; `f64::INFINITY` where a region carries no |
| fn | `a_priori_metrics` | `meshing/mod.rs:285` | `pub fn a_priori_metrics(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String>` —  The isotropic metric at every vertex of a meshed model that the model itself asks for: the size of the region(s) the vertex belongs |
| fn | `geometry_floor` | `meshing/mod.rs:291` | `pub fn geometry_floor(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<2>>, String>` —  The sizes the geometry insists on whatever the solution does: `Size` rules on regions and edges, the layers of an absorbing region, and |
| fn | `adapt_to_metrics` | `meshing/mod.rs:379` | `pub fn adapt_to_metrics(g: MeshedGeom, metrics: &[Sym<2>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String>` —  Adapts a meshed model to a metric given at every vertex of its mesh. The metrics are graded and interpolated on the old mesh; the |
| fn | `adapt_to_field` | `meshing/mod.rs:393` | `pub fn adapt_to_field(g: MeshedGeom, field: &dyn MetricField<2>, a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String>` —  Adapts a meshed model to a metric field defined anywhere in the plane, for example the final metric of an earlier, slightly different |
| fn | `periodic_labels` | `meshing/mod.rs:412` | `pub(crate) fn periodic_labels(g: &MeshedGeom) -> Vec<i32>` —  Boundary labels of the edges (2D) or faces (3D) that take part in a periodic pair of the model. Their chains or patches stay |
| struct | `GeomOracle2` | `meshing/mod.rs:449` |  The exact geometry as the oracle of the adaptation: where a vertex of a boundary or interface chain may go is the closest point on |
| enum | `Reuse` | `meshing/mod.rs:511` | variants: `Hit`, `Meshed` —  How a cached request was answered. |
| struct | `CacheStats` | `meshing/mod.rs:521` | fields: `requests`, `hits`, `eval_seconds`, `mesh_seconds` —  Counters of a `MeshCache`. |
| struct | `MeshCache` | `meshing/mod.rs:536` | fields: `stats` —  Meshes a series of related models (a parameter sweep, a refinement study, the same model for several physics) without repeating work. |
| fn | `MeshPlan::simple` | `meshing/mod.rs:89` | `pub fn simple(h: f64) -> Self` —  The simple mesher at element size `h`. |
| fn | `MeshPlan::adaptive` | `meshing/mod.rs:94` | `pub fn adaptive(h: f64, scale: f64) -> Self` —  The adaptive mesher: start at element size `h`, adapt to `scale * h` (and to the curvature of the edges at that scale). |
| fn | `GeomOracle2::from_work` | `meshing/mod.rs:456` | `pub fn from_work(ev: &'a Evaluated, labels: &LabelMap, w: &Work2) -> Self` |
| fn | `MeshCache::new` | `meshing/mod.rs:552` | `pub fn new(capacity: usize) -> Self` —  A cache that keeps the `capacity` most recent distinct requests. |
| fn | `MeshCache::mesh` | `meshing/mod.rs:562` | `pub fn mesh(&mut self, list: &FeatureList, plan: &MeshPlan) -> Result<(MeshedGeom, Reuse), String>` —  The meshed model for `list` and `plan`, from the cache when it is there. |
| fn | `MeshCache::mesh_extruded` | `meshing/mod.rs:589` | `pub fn mesh_extruded(&mut self, list: &FeatureList, plan: &MeshPlan, layers: &[crate::meshing::mesher::Layer]) -> Result<(MeshedGeom, Reuse), String>` —  The section of `list` swept along y through `layers`. The section is the cached part: changing only the layers re-extrudes it. |
| struct | `PhysicsSizing` | `meshing/physics_size.rs:30` | fields: `points_per_wavelength`, `pml_relax`, `eigen_headroom` —  How finely the waves are resolved. |
| fn | `wave_speed_range` | `meshing/physics_size.rs:56` | `pub fn wave_speed_range(c: &M6, rho: f64) -> (f64, f64)` —  Slowest and fastest plane-wave speed of a medium with stiffness `c` (solver Voigt order xx, yy, zz, xy, yz, xz; engineering shear) and |
| fn | `cell_sizes_for` | `meshing/physics_size.rs:93` | `pub fn cell_sizes_for( sizing: &Sizing, resolve: &dyn Fn(&str) -> Option<i32>, cell_labels: &[i32], pml: &dyn Fn(i32) -> bool, sz: &PhysicsSizing, ) -> Result<Vec<f64>, String>` —  Target element size of every cell: the shortest wavelength of its domain, divided by the headroom and by the points per wavelength, and |
| fn | `vertex_sizes_for` | `meshing/physics_size.rs:116` | `pub fn vertex_sizes_for( sizing: &Sizing, resolve: &dyn Fn(&str) -> Option<i32>, mesh: &SimplexMesh, pml: &dyn Fn(i32) -> bool, sz: &PhysicsSizing, ) -> Result<Vec<f64>, String>` —  Target element size at every vertex: the smallest of the cells that meet there. A vertex of cells without a wave gets `f64::INFINITY`. |
| fn | `metrics_for` | `meshing/physics_size.rs:135` | `pub fn metrics_for<const D: usize>( sizing: &Sizing, resolve: &dyn Fn(&str) -> Option<i32>, mesh: &SimplexMesh, pml: &dyn Fn(i32) -> bool, sz: &PhysicsSizing, ceiling: f64, ) -> Result<Vec<Sym<D>>, String>` —  Isotropic metric of the physics size at every vertex of a mesh of any dimension (`D` = 2 or 3). Vertices without a wave get `ceiling` (a |
| fn | `adapt_tets_for` | `meshing/physics_size.rs:153` | `pub fn adapt_tets_for( mesh: &SimplexMesh, sizing: &Sizing, resolve: &dyn Fn(&str) -> Option<i32>, pml: &dyn Fn(i32) -> bool, sz: &PhysicsSizing, ceiling: f64, grade: f64, opt: &crate::adapt::ops3::AdaptOptions3, ) -> Result<(SimplexMesh, crate::adapt::drive::AdaptReport), String>` —  Adapts a tetrahedral mesh to a neutral `Sizing`: the sizes of `vertex_sizes_for` (capped at `ceiling`), graded with slope `grade`, |
| fn | `optical_wavelength` | `meshing/physics_size.rs:176` | `pub fn optical_wavelength(eps: &[[f64; 3]; 3], mu: f64, lam0: f64) -> f64` —  Free-space wavelength `lam0` shortened by the largest refractive index of a medium with the relative permittivity tensor `eps` and the |
| fn | `speed_of_light` | `meshing/physics_size.rs:182` | `pub fn speed_of_light() -> f64` —  Vacuum speed of light (m/s). |
| fn | `elastic_wavelength` | `meshing/physics_size.rs:188` | `pub fn elastic_wavelength(c: &M6, rho: f64, f: f64) -> f64` —  Shortest wavelength of the elastic waves of a medium with stiffness `c` (solver Voigt order) and density `rho` at frequency `f`: the |
| fn | `PhysicsSizing::ppw` | `meshing/physics_size.rs:48` | `pub fn ppw(&self, order: usize) -> f64` —  Points per wavelength for an element order. |
| enum | `MesherRequest` | `meshing/request.rs:18` | variants: `Simple`, `Adaptive` —  Which mesher makes the cells. |
| struct | `WavesRequest` | `meshing/request.rs:49` | fields: `sizing`, `points_per_wavelength`, `pml_relax`, `eigen_headroom` —  The waves of a study, as the physics reports them (`Sizing`) and how finely they are resolved. |
| struct | `MeshRequest` | `meshing/request.rs:64` | fields: `h`, `max_dev`, `mesher`, `structured_pml`, `waves`, `curved`, `gap_cells`, `z_coarsen`, `max_cells`, `dz_min`, `max_z_slabs`, `z_growth_ratio`, ... —  A mesh asked for. |
| struct | `DomainReport` | `meshing/request.rs:146` | fields: `id`, `name`, `cells`, `wave_size` —  One domain of the meshed model. |
| struct | `AdaptSummary` | `meshing/request.rs:157` | fields: `iterations`, `converged`, `capped`, `q_min`, `cells_start`, `cells_end` —  What the adaptation did. |
| struct | `MeshReport` | `meshing/request.rs:169` | fields: `dim`, `vertices`, `cells`, `facets`, `domains`, `adaptation`, `structured_layers`, `curved`, `warnings` —  A mesh reported. |
| struct | `CurvedSummary` | `meshing/request.rs:191` | fields: `edges`, `cells`, `min_scaled_jacobian` —  The curved geometry of a 3D mesh (ADR 82). |
| fn | `mesh_request` | `meshing/request.rs:293` | `pub fn mesh_request(list: &FeatureList, req: &MeshRequest) -> Result<(MeshedGeom, MeshReport), String>` —  Meshes `list` as `req` asks (any dimension) and reports what came out. |
| fn | `request_schema` | `meshing/request.rs:317` | `pub fn request_schema() -> Value` —  JSON Schema of [`MeshRequest`]. |
| fn | `report_schema` | `meshing/request.rs:379` | `pub fn report_schema() -> Value` —  JSON Schema of [`MeshReport`]. |
| fn | `MeshRequest::plan` | `meshing/request.rs:118` | `pub fn plan(&self) -> MeshPlan` —  The plan this request names. |
| fn | `CurvedSummary::from_mesh` | `meshing/request.rs:204` | `pub fn from_mesh(m: &crate::mesh::SimplexMesh) -> Option<Self>` —  Summarizes the curved and/or quarter-point geometry of a 3D mesh, evaluating quarter-point cells |
| struct | `OverlayOracle2` | `meshing/sheet3d_amr.rs:40` | fields: `ov` —  Exact 2D geometry oracle over `LayerStack::overlay` in the `(x, y)` sweep plane. |
| struct | `OverlaySheetFeatures` | `meshing/sheet3d_amr.rs:180` | fields: `rim_segments`, `corners`, `plane_axes` —  Extracted 2D sheet perimeter segments, sharp corners, and sweep-plane world axes from a 3D stack `Complex`. |
| fn | `extract_overlay_sheet_features` | `meshing/sheet3d_amr.rs:190` | `pub fn extract_overlay_sheet_features(cx: &Complex) -> Option<OverlaySheetFeatures>` —  Extracts 2D sheet rim segments and sharp corners from a 3D prismatic stack `Complex`. |
| fn | `adapt_overlay_2d` | `meshing/sheet3d_amr.rs:262` | `pub fn adapt_overlay_2d( ov: &Complex, m2: &SimplexMesh, metrics: &[Sym<2>], a: &Adaptive, ) -> Result<(SimplexMesh, AdaptReport), String>` —  Adapts a 2D overlay triangulation `m2` over `ov` (`st.overlay`) to vertex metrics `metrics`, |
| fn | `adapt_overlay_2d_field` | `meshing/sheet3d_amr.rs:280` | `pub fn adapt_overlay_2d_field( ov: &Complex, m2: &SimplexMesh, field: &dyn MetricField<2>, a: &Adaptive, ) -> Result<(SimplexMesh, AdaptReport), String>` —  Adapts a 2D overlay triangulation `m2` over `ov` (`st.overlay`) to a general [`MetricField<2>`]. |
| fn | `build_overlay_mesh_2d` | `meshing/sheet3d_amr.rs:309` | `pub fn build_overlay_mesh_2d(ev: &Evaluated, plan: &MeshPlan) -> Result<SimplexMesh, String>` —  Builds the initial 2D conformal triangulation `m2` over `ev.complex().stack.overlay` (`cell_labels` |
| fn | `sheet_rim_vertices_3d` | `meshing/sheet3d_amr.rs:315` | `pub fn sheet_rim_vertices_3d(mesh: &SimplexMesh, sheet_labels: &[i32]) -> Vec<bool>` —  Identifies vertices of a 3D `SimplexMesh` lying on open 1D perimeter rims of `sheet_labels` |
| fn | `adapt_sheet_full_3d` | `meshing/sheet3d_amr.rs:344` | `pub fn adapt_sheet_full_3d(g: MeshedGeom, metrics: &[Sym<3>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String>` —  Adapts a 3D [`MeshedGeom`] containing interior 2D `Op::Sheet` boundaries (`l0 == l1` or `l0 != l1`) |
| fn | `adapt_sheet_full_3d_with_plan` | `meshing/sheet3d_amr.rs:350` | `pub fn adapt_sheet_full_3d_with_plan( g: MeshedGeom, metrics: &[Sym<3>], a: &Adaptive, plan: &MeshPlan, waves: Option<(&Sizing, &PhysicsSizing)>, ) -> Result<(MeshedGeom, AdaptReport), String>` —  Same as [`adapt_sheet_full_3d`], and re-applies `AF-5` (`curved`) + `AF-2` (`quarter_point`) |
| fn | `OverlayOracle2::from_work` | `meshing/sheet3d_amr.rs:48` | `pub fn from_work(ov: &'a Complex, w: &Work2) -> Self` —  Builds the overlay oracle from a classified 2D working mesh `w` over `ov`. |
| fn | `OverlayOracle2::pinned_vertices` | `meshing/sheet3d_amr.rs:129` | `pub fn pinned_vertices(ov: &Complex, m2: &SimplexMesh) -> Vec<u32>` —  Mesh vertex indices of `m2` that coincide with `ov.vertices` in `(x, y)` (pinned as `Corner`). |
| struct | `Sheet3dTransitionRef` | `meshing/sheet3d_curve.rs:42` | fields: `map_v`, `rank` —  Borrowed transition mapping for one 2D vertical coarsening step `M_m -> M_{m+1}` in `sheet3d`. |
| fn | `relax_edge_nodes` | `meshing/sheet3d_curve.rs:50` | `pub fn relax_edge_nodes( pos: &(dyn Fn(u32) -> [f64; 3] + Sync), cells: &[[u32; 4]], mut cand: EdgeNodes, ) -> (EdgeNodes, CurveReport)` —  Validates and relaxes candidate mid-edge geometry nodes `cand` (sorted by `[u32; 2]` vertex pair) |
| fn | `project_sheet3d_curved_edges` | `meshing/sheet3d_curve.rs:65` | `pub fn project_sheet3d_curved_edges( cx: &Complex, m2: &SimplexMesh, bnd_edges2: &[(u32, u32, u32)], z_planes: &[f64], used_v: &[u32], slab_layer: &[Option<usize>], level_of_zp: &[usize], transitions: &[Sheet3dTransitionRef<'_>], mesh: &mut SimplexMesh, ) -> CurveReport` —  Projects midpoints of 2D boundary/interface edges `bnd_edges2` lying on `Curve::Circle` or |
| fn | `apply_sheet3d_curved_edges` | `meshing/sheet3d_curve.rs:239` | `pub fn apply_sheet3d_curved_edges(g: &mut MeshedGeom) -> CurveReport` —  Reconstructs the prismatic 2D column/slab topology of a `MeshedGeom` built by |
| const | `DEFAULT_QUARTER_POINT_XI` | `meshing/singularity.rs:29` | `f64` —  Default Barsoum quarter-point parameter (`0.25` from the singular rim vertex). |
| const | `CORNER_GRADING_P` | `meshing/singularity.rs:32` | `f64` —  Polynomial degree `p = 2` for Babuska corner grading (`p + 0.5 = 2.5`). |
| const | `CORNER_RADIUS_FACTOR` | `meshing/singularity.rs:35` | `f64` —  Corner grading influence radius `R_c / h_edge = 2.0`. |
| const | `CORNER_H_MIN_FACTOR` | `meshing/singularity.rs:38` | `f64` —  Mandatory corner grading floor `h_min / h_edge = 0.05` (prevents runaway subdivision). |
| const | `MORRISON_LEWIS_TABLE` | `meshing/singularity.rs:48` | `[(f64, f64)` —  Tabulated Morrison-Lewis (1976) 3D Laplace corner singularity exponents `(omega_rad, alpha)` for a flat |
| fn | `morrison_lewis_exponent` | `meshing/singularity.rs:76` | `pub fn morrison_lewis_exponent(omega_rad: f64, eps_ratio: f64) -> f64` —  Piecewise-linearly interpolated Morrison-Lewis 3D corner singularity exponent `alpha(omega, eps_1 / eps_2)` |
| fn | `corner_singularity_exponent` | `meshing/singularity.rs:98` | `pub fn corner_singularity_exponent(omega_rad: f64, eps_ratio: f64) -> f64` —  Alias for [`morrison_lewis_exponent`]. |
| fn | `corner_grading_h` | `meshing/singularity.rs:106` | `pub fn corner_grading_h(d: f64, h_edge: f64, alpha: f64) -> f64` —  Babuska's optimal radical corner grading law: |
| fn | `babuska_corner_h` | `meshing/singularity.rs:126` | `pub fn babuska_corner_h(d: f64, h_edge: f64, alpha: f64) -> f64` —  Alias for [`corner_grading_h`]. |
| struct | `SheetCorner` | `meshing/singularity.rs:132` | fields: `vertex`, `pos`, `omega_rad`, `alpha`, `h_edge` —  A detected 2D conductive sheet corner on `st.overlay` or `m2`. |
| fn | `corner_graded_h_at` | `meshing/singularity.rs:163` | `pub fn corner_graded_h_at(p: [f64; 2], corners: &[SheetCorner], base_h: f64) -> f64` —  Evaluates the corner-graded 2D element size at point `p` across all sharp corners in `corners`, |
| fn | `coupled_corner_dz` | `meshing/singularity.rs:185` | `pub fn coupled_corner_dz(base_dz: f64, corners: &[SheetCorner]) -> f64` —  Coupled vertical slab step `dz` at a sheet plane when `corner_grading` is active: |
| fn | `measure_sheet_corner_angles_m2` | `meshing/singularity.rs:201` | `pub fn measure_sheet_corner_angles_m2(m2: &SimplexMesh, sheet_covered_cells: &[bool]) -> Vec<Option<f64>>` —  Measures the interior sheet opening angle `omega in (0, 2*pi)` at every vertex of a 2D overlay triangulation `m2` |
| fn | `measure_overlay_sheet_corners` | `meshing/singularity.rs:267` | `pub fn measure_overlay_sheet_corners(cx: &Complex, z_s: f64) -> Vec<Option<f64>>` —  Measures the interior sheet opening angle `omega in (0, 2*pi)` at every `st.overlay` vertex `0..ov.vertices.len()` |
| fn | `collect_overlay_sheet_corners` | `meshing/singularity.rs:337` | `pub fn collect_overlay_sheet_corners( cx: &Complex, z_s: f64, ov_edge_h: &[f64], eps_ratio: f64, ) -> Vec<SheetCorner>` —  Collects all active sharp [`SheetCorner`]s on sheet plane `z_s` of `cx`, pairing each corner with its |
| fn | `analytical_qp_jacobian_ratio` | `meshing/singularity.rs:384` | `pub fn analytical_qp_jacobian_ratio(n_rim_in_tet: usize, lam_out: f64, xi_qp: f64) -> f64` —  Exact analytical 3D Jacobian determinant ratio `det(J(lam)) / det(J_straight)` for a tetrahedron whose |
| fn | `sheet_rim_vertices` | `meshing/singularity.rs:397` | `pub fn sheet_rim_vertices(mesh: &SimplexMesh, sheet_labels: &[i32]) -> Vec<bool>` —  Identifies all 1D open perimeter vertices (`rim_verts[v] == true`) of interior conductive sheet facets |
| struct | `QuarterPointReport` | `meshing/singularity.rs:426` | fields: `rim_vertices`, `shifted_edges`, `qp_cells`, `clamped_curved_edges`, `unshifted_edges`, `min_interior_jacobian_ratio` —  Summary report returned by [`apply_quarter_point_nodes`]. |
| fn | `apply_quarter_point_nodes` | `meshing/singularity.rs:468` | `pub fn apply_quarter_point_nodes( mesh: &mut SimplexMesh, rim_verts: &[bool], xi_qp: f64, curved_rim_edges: &BTreeSet<[u32; 2]>, ) -> QuarterPointReport` —  Applies 3D Barsoum quarter-point mid-edge geometry nodes (`EdgeNodes`) to `mesh`: |
| fn | `SheetCorner::is_sharp` | `meshing/singularity.rs:145` | `pub fn is_sharp(&self) -> bool` —  True when the interior angle deviates from a smooth `180 deg` straight edge (`\|omega - pi\| > 0.05 rad`). |
| fn | `SheetCorner::h_at` | `meshing/singularity.rs:150` | `pub fn h_at(&self, p: [f64; 2]) -> f64` —  Graded element size `h(d)` at point `p` induced by this corner. |
| fn | `SheetCorner::h_min` | `meshing/singularity.rs:156` | `pub fn h_min(&self) -> f64` —  Minimum floor size `h_min = 0.05 * h_edge` at the corner apex (`d = 0`). |
| struct | `Sizing` | `meshing/sizing.rs:14` | fields: `order`, `headroom`, `wavelengths` —  What the mesher needs from a physics to size the elements before any solution exists. |
| fn | `label_names` | `meshing/sizing.rs:55` | `pub fn label_names(name: &str) -> Option<i32>` —  Resolver for meshes without names: the key of a region is its decimal domain label. |
| fn | `Sizing::by_label` | `meshing/sizing.rs:25` | `pub fn by_label(order: u8, headroom: f64, wavelengths: &BTreeMap<i32, f64>) -> Sizing` —  A request for a mesh that only has domain numbers: the keys are the decimal labels. |
| fn | `Sizing::expand_wildcards` | `meshing/sizing.rs:32` | `pub fn expand_wildcards(&self, names: &[&str]) -> Sizing` —  The same request with a name that ends in `*` replaced by the regions it stands for: `Air*` is `Air` and its absorbing pieces |
| fn | `Sizing::resolve` | `meshing/sizing.rs:44` | `pub fn resolve(&self, domain: &dyn Fn(&str) -> Option<i32>) -> Result<BTreeMap<i32, f64>, String>` —  The wavelengths by domain label, the names resolved by `domain`. A region the model does not have is an error. |

Trait implementations in this component: `Adaptive: Default`, `MeshPlan: From`, `GeomOracle2: Geometry2`, `MeshCache: Default`, `PhysicsSizing: Default`, `OverlayOracle2: Geometry2`, `QuarterPointReport: Default`
