# 3D meshing: code

C4 level 4 (code) of the component **3D meshing** (`femcore.meshing3d`) in the container **femcore** (crate `femcore`). [Component view of the container](../../3-components/femcore.md) · [All code pages](../README.md)

> The lattice tetrahedral mesher on Domain3, Rivara refinement, the adapter from a 3D geom model, adaptive meshing driven by exact faces, face sizes and max_dev, and the curved (quadratic) edges of boundary cells.

**Decisions:** [ADR 40](../../../adr/0040-3d-geom-models-meshed-by-the-lattice-mesher.md), [ADR 50](../../../adr/0050-wavelength-sizing-in-3d.md), [ADR 51](../../../adr/0051-adaptive-meshing-of-3d-geom-models.md), [ADR 55](../../../adr/0055-lattice-mesher-welds-slivers-and-is-independent-of-the-thread-count.md), [ADR 82](../../../adr/0082-curved-order-2-geometry.md), [ADR 83](../../../adr/0083-face-size-rules-and-max-dev-in-the-lattice-mesher.md)  
**Size:** 5 files, 2090 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/femcore/src/meshing/adapt3d.rs` | 320 | Adaptive meshing of a 3D `geom` model (ADR 51): the metric adaptation of `adapt::drive3` driven by the exact faces of the model. |
| `crates/femcore/src/meshing/curve3.rs` | 276 | Curved edges for 3D meshes (ADR 82): the mid-edge geometry nodes that make the tetrahedra around a curved surface quadratic (isoparametric) instead of straight. |
| `crates/femcore/src/meshing/geomdrive3.rs` | 338 | From a 3D `geom` model to a tetrahedral mesh: the adapter between `geom::Evaluated` and the lattice mesher (`mesher3d`), which only knows the trait `Domain3`. |
| `crates/femcore/src/meshing/mesher3d.rs` | 912 | Tetrahedral meshes of 3D models: a lattice fill that is fitted to the exact geometry of any model that can say which region a point is in. |
| `crates/femcore/src/meshing/refine3d.rs` | 244 | Local refinement of a tetrahedral mesh by longest-edge (Rivara) bisection. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `geom` (geom)](../../3-components/geom.md) | eval, labels, math, model, select | Aabb, EntityKind, Evaluated, LabelMap, Model, RegionKind |
| [Adaptation drivers](./adaptdrive.md) | adapt::drive, adapt::drive3, adapt::geom3 | AdaptReport, FacetGeometry, Geometry3, adapt3_with |
| [Adaptation operators](./adaptops.md) | adapt::ops3, adapt::topo, adapt::topo3 | AdaptOptions3, Sig, TopoOptions3, Work3 |
| [Finite-element basics](./fem.md) | fem | CurvedTet, quad_rule_collapsed |
| [Mesh](./mesh.md) | mesh, mesh::cells | CellMesh, CellType, SimplexMesh |
| [Meshing facade](./meshing.md) | meshing, meshing::geomdrive | Adaptive, GeomMeshOpts, MeshedGeom, evaluate, periodic_labels |
| [Metric fields](./metric.md) | adapt::field3, adapt::metric, adapt::work | BackgroundMetric3, SimplexMesh, Sym, isotropic, simplex_to_adapt |
| [Runtime services](./runtime.md) | par, progress | check, lap, threads |

## Used by

| Component | Modules | Items |
|---|---|---|
| [Meshing facade](./meshing.md) | meshing::curve3, meshing::geomdrive3 | Sizes3, scaled_jacobian_min |
| [Structured absorbing layers](./structpml.md) | meshing::adapt3d | FaceTable, interface_facets |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace meshing_adapt3d_rs {
    class FaceTable {
        <<struct>>
        +new(ev, labels) Result
        +nearest_of(p, cands) Option~i32~
        +interface_label(p, a, b) i32
    }
    class GeomOracle3 {
        <<struct>>
        +from_work(table, w) Self
    }
    class fns_meshing_adapt3d_rs {
        <<module meshing/adapt3d.rs>>
        +interface_facets(point, cells, cell_label, table) tuple
        +adapt_to_metrics3(g, metrics, a) Result
        +vertex_sizes3(g, opts, a, geometry_only) Result
        +a_priori_metrics3(g, opts, a) Result
        +geometry_floor3(g, opts, a) Result
    }
  }
  namespace meshing_curve3_rs {
    class CurveReport {
        <<struct>>
        +curved_edges usize
        +curved_cells usize
        +relaxed usize
        +straightened usize
        +min_scaled_jacobian f64
    }
    class fns_meshing_curve3_rs {
        <<module meshing/curve3.rs>>
        +curve_edges(pos, cells, facets, dom) tuple
        +scaled_jacobian_min(mesh, c) f64
        +curve_boundary_edges(mesh, dom) CurveReport
        +curve_solver_mesh(mesh, dom) CurveReport
        +EdgeNodes type
        +MIN_SCALED_JACOBIAN const
    }
  }
  namespace meshing_geomdrive3_rs {
    class GeomDomain3 {
        <<struct>>
        +new(ev, labels, h) Result
        +new_with(ev, labels, h, sizes) Result
        +with_max_dev(max_dev) Self
        +check_periodic_faces() Result
        +region_sizes() slice
        +is_meshed(r) bool
        +coarsest_meshed_size() f64
    }
    class Sizes3 {
        <<struct>>
        +region_h Vec~f64~
        +skip Vec~bool~
    }
    class fns_meshing_geomdrive3_rs {
        <<module meshing/geomdrive3.rs>>
        +mesh_geom_3d(ev, opts) Result
        +mesh_geom_3d_report(ev, opts) Result
        +mesh_geom_3d_with(ev, opts, sizes) Result
        +sphere_cavity(r, h) Result
        +sphere_cavity_with(r, opts) Result
    }
  }
  namespace meshing_mesher3d_rs {
    class Domain3 {
        <<interface>>
        +region_at(p) Option~u32~
        +n_regions() usize
        +region_box(r) tuple
        +region_size(r) Option~f64~
        +region_label(r) i32
        +boundary_label(p, a, b) i32
        +project_edge(mid, pa, pb) Optionf64
        +edge_on_crease(_pa, _pb) bool
        +more_2_methods()
    }
    class Spec3d {
        <<struct>>
        +h f64
        +alpha f64
        +min_layers f64
        +new(h) Self
    }
    class Report3d {
        <<struct>>
        +lattice_vertices usize
        +warped usize
        +cut_points usize
        +cells usize
        +refined usize
        +curve Option~CurveReport~
    }
    class fns_meshing_mesher3d_rs {
        <<module meshing/mesher3d.rs>>
        +mesh_3d(dom, spec) Result
        +MAX_SURFACES const
    }
  }
  namespace meshing_refine3d_rs {
    class EdgeKind {
        <<enumeration>>
        Interior
        Boundary
        Interface
    }
    class RefineReport {
        <<struct>>
        +rounds usize
        +bisections usize
        +cells_before usize
        +cells_after usize
    }
    class fns_meshing_refine3d_rs {
        <<module meshing/refine3d.rs>>
        +refine_to_size(mesh, size, project, ratio, max_cells) Result
    }
  }
  GeomDomain3 ..|> Domain3
  Report3d --> CurveReport
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `FaceTable` | `meshing/adapt3d.rs:35` |  The faces of a model as the adaptation needs them: the domain numbers on their two sides, their boundary number and their box. |
| struct | `GeomOracle3` | `meshing/adapt3d.rs:124` |  The exact geometry as the oracle of the 3D adaptation. |
| fn | `interface_facets` | `meshing/adapt3d.rs:159` | `pub fn interface_facets(point: &dyn Fn(u32) -> [f64; 3], cells: &[u32], cell_label: &[i32], table: &FaceTable) -> (Vec<[u32` —  Faces between cells of different domains, labelled, oriented out of the lower-numbered domain. |
| fn | `adapt_to_metrics3` | `meshing/adapt3d.rs:199` | `pub fn adapt_to_metrics3(g: MeshedGeom, metrics: &[Sym<3>], a: &Adaptive) -> Result<(MeshedGeom, AdaptReport), String>` —  Adapts a meshed 3D model to a metric given at every vertex of its mesh: graded, the faces of the model as the oracle, the labels of |
| fn | `vertex_sizes3` | `meshing/adapt3d.rs:236` | `pub fn vertex_sizes3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive, geometry_only: bool) -> Result<Vec<f64>, String>` —  The element size every vertex of a meshed 3D model is asked for by the model itself: the size of the regions it belongs to (the |
| fn | `a_priori_metrics3` | `meshing/adapt3d.rs:287` | `pub fn a_priori_metrics3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<3>>, String>` —  The isotropic metric the model itself asks for at every vertex (see [`vertex_sizes3`]). |
| fn | `geometry_floor3` | `meshing/adapt3d.rs:292` | `pub fn geometry_floor3(g: &MeshedGeom, opts: &GeomMeshOpts, a: &Adaptive) -> Result<Vec<Sym<3>>, String>` —  The sizes the geometry insists on whatever the solution does, as metrics: the floor of an estimator metric. |
| fn | `FaceTable::new` | `meshing/adapt3d.rs:44` | `pub fn new(ev: &'a Evaluated, labels: &LabelMap) -> Result<Self, String>` |
| fn | `FaceTable::nearest_of` | `meshing/adapt3d.rs:79` | `pub fn nearest_of(&self, p: [f64; 3], cands: impl Iterator<Item = usize>) -> Option<i32>` —  The boundary number of the nearest of the faces `cands` to `p`. |
| fn | `FaceTable::interface_label` | `meshing/adapt3d.rs:102` | `pub fn interface_label(&self, p: [f64; 3], a: i32, b: i32) -> i32` —  The boundary number of the face a facet with centroid `p` between the domains `a` and `b` lies on: the nearest face between exactly |
| fn | `GeomOracle3::from_work` | `meshing/adapt3d.rs:132` | `pub fn from_work(table: &'a FaceTable<'a>, w: &Work3) -> Self` |
| type | `EdgeNodes` | `meshing/curve3.rs:24` |  The curved edges: sorted vertex pair and the position of the mid-edge node. |
| struct | `CurveReport` | `meshing/curve3.rs:28` | fields: `curved_edges`, `curved_cells`, `relaxed`, `straightened`, `min_scaled_jacobian` —  What `curve_boundary_edges` did. |
| const | `MIN_SCALED_JACOBIAN` | `meshing/curve3.rs:48` | `f64` —  A cell must keep this fraction of its straight Jacobian determinant everywhere. |
| fn | `curve_edges` | `meshing/curve3.rs:69` | `pub fn curve_edges( pos: &(dyn Fn(u32) -> [f64; 3] + Sync), cells: &[[u32; 4]], facets: &[[u32; 3]], dom: &dyn Domain3, ) -> (EdgeNodes, CurveReport)` —  The core: curved edges for the tetrahedra `cells` with vertex positions `pos`, from the edges of the facets `facets` (boundary |
| fn | `scaled_jacobian_min` | `meshing/curve3.rs:223` | `pub fn scaled_jacobian_min(mesh: &SimplexMesh, c: usize) -> f64` —  The smallest ratio of the Jacobian determinant to the straight cell's over cell `c` of a solver mesh (ten nodes and a degree-6 rule). |
| fn | `curve_boundary_edges` | `meshing/curve3.rs:252` | `pub fn curve_boundary_edges(mesh: &mut CellMesh<3>, dom: &dyn Domain3) -> CurveReport` —  Curves the edges of the boundary and interface facets of a tetrahedral `CellMesh` (replacing its list of curved edges). Meshes with |
| fn | `curve_solver_mesh` | `meshing/curve3.rs:265` | `pub fn curve_solver_mesh(mesh: &mut SimplexMesh, dom: &dyn Domain3) -> CurveReport` —  The same for a solver mesh (after an adaptation, say): replaces its curved edges. |
| struct | `GeomDomain3` | `meshing/geomdrive3.rs:25` |  A `geom` model as the geometry of the lattice mesher. |
| struct | `Sizes3` | `meshing/geomdrive3.rs:279` | fields: `region_h`, `skip` —  What a caller asks of the 3D mesh beyond the model's own sizes. |
| fn | `mesh_geom_3d` | `meshing/geomdrive3.rs:290` | `pub fn mesh_geom_3d(ev: Evaluated, opts: &GeomMeshOpts) -> Result<MeshedGeom, String>` —  Meshes a 3D model with tetrahedra at element size `opts.h`: cells labelled by region (domain numbers), boundary facets by geom face. |
| fn | `mesh_geom_3d_report` | `meshing/geomdrive3.rs:295` | `pub fn mesh_geom_3d_report(ev: Evaluated, opts: &GeomMeshOpts) -> Result<(MeshedGeom, Report3d, CellMesh<3>), String>` —  Like `mesh_geom_3d`, with the report of the mesher and the native mesh (all tetrahedra) the solver mesh was made from. |
| fn | `mesh_geom_3d_with` | `meshing/geomdrive3.rs:300` | `pub fn mesh_geom_3d_with(ev: Evaluated, opts: &GeomMeshOpts, sizes: &Sizes3) -> Result<(MeshedGeom, Report3d, CellMesh<3>), String>` —  Like `mesh_geom_3d_report`, with the extra requests of `sizes` (ADR 50, 52). |
| fn | `sphere_cavity` | `meshing/geomdrive3.rs:326` | `pub fn sphere_cavity(r: f64, h: f64) -> Result<MeshedGeom, String>` —  A spherical cavity of radius `r` centred at the origin, meshed at element size `h`: one domain (`cavity`) and one boundary, the surface of |
| fn | `sphere_cavity_with` | `meshing/geomdrive3.rs:331` | `pub fn sphere_cavity_with(r: f64, opts: &GeomMeshOpts) -> Result<MeshedGeom, String>` —  `sphere_cavity` with the mesh options of the caller (`curved`, `max_dev`). |
| fn | `GeomDomain3::new` | `meshing/geomdrive3.rs:40` | `pub fn new(ev: &'a Evaluated, labels: &'a LabelMap, h: f64) -> Result<Self, String>` |
| fn | `GeomDomain3::new_with` | `meshing/geomdrive3.rs:46` | `pub fn new_with(ev: &'a Evaluated, labels: &'a LabelMap, h: f64, sizes: &Sizes3) -> Result<Self, String>` —  Like `new`, with the extra requests of `sizes`: a finer element size per region (the wavelength of a study, ADR 50) and regions that |
| fn | `GeomDomain3::with_max_dev` | `meshing/geomdrive3.rs:112` | `pub fn with_max_dev(mut self, max_dev: Option<f64>) -> Self` —  Adds the chord limit `max_dev` (metres) of the request to every curved face: the element size at which a chord stays within it |
| fn | `GeomDomain3::check_periodic_faces` | `meshing/geomdrive3.rs:125` | `pub fn check_periodic_faces(&self) -> Result<(), String>` —  A face with a size request that a periodic face is within reach of: refused (the two faces of a pair would be refined differently |
| fn | `GeomDomain3::region_sizes` | `meshing/geomdrive3.rs:149` | `pub fn region_sizes(&self) -> &[f64]` —  The element size of each geom region (global size lowered by `Size` rules, absorbing layers and requested sizes). |
| fn | `GeomDomain3::is_meshed` | `meshing/geomdrive3.rs:154` | `pub fn is_meshed(&self, r: u32) -> bool` —  Whether a region contributes cells. |
| fn | `GeomDomain3::coarsest_meshed_size` | `meshing/geomdrive3.rs:159` | `pub fn coarsest_meshed_size(&self) -> f64` —  The largest element size any meshed region asks for. |
| trait | `Domain3` | `meshing/mesher3d.rs:31` | methods: `region_at`, `n_regions`, `region_box`, `region_size`, `region_label`, `boundary_label`, `project_edge`, `edge_on_crease`, `size_at`, `asks_size_at` —  The geometry the mesher fits to. Regions are numbered `0..n_regions()`; a point outside every region is in the void. |
| struct | `Spec3d` | `meshing/mesher3d.rs:64` | fields: `h`, `alpha`, `min_layers` —  Settings of the 3D mesher. |
| struct | `Report3d` | `meshing/mesher3d.rs:83` | fields: `lattice_vertices`, `warped`, `cut_points`, `cells`, `refined`, `curve` —  A report of how the mesh was made. |
| const | `MAX_SURFACES` | `meshing/mesher3d.rs:97` | `usize` —  Largest number of cutting surfaces: the solid and the regions. |
| fn | `mesh_3d` | `meshing/mesher3d.rs:303` | `pub fn mesh_3d(dom: &dyn Domain3, spec: &Spec3d) -> Result<(CellMesh<3>, Report3d), String>` —  Meshes the model with tetrahedra. Cells carry `region_label`, boundary facets (outer faces and faces between regions) carry |
| fn | `Spec3d::new` | `meshing/mesher3d.rs:76` | `pub fn new(h: f64) -> Self` |
| enum | `EdgeKind` | `meshing/refine3d.rs:18` | variants: `Interior`, `Boundary`, `Interface` —  Where an edge that is being split lies. |
| struct | `RefineReport` | `meshing/refine3d.rs:27` | fields: `rounds`, `bisections`, `cells_before`, `cells_after` |
| fn | `refine_to_size` | `meshing/refine3d.rs:72` | `pub fn refine_to_size( mesh: &CellMesh<3>, size: &(dyn Fn(P3) -> f64 + Sync), project: &(dyn Fn(P3, P3, P3, EdgeKind) -> P3 + Sync), ratio: f64, max_cells: usize, ) -> Result<(CellMesh<3>, RefineReport), String>` —  Refines `mesh` until no cell has an edge longer than `ratio` times `size(centroid)`. |

Trait implementations in this component: `GeomOracle3: Geometry3`, `CurveReport: Default`, `GeomDomain3: Domain3`
