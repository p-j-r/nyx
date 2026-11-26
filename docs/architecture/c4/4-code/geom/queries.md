# Queries and selectors: code

C4 level 4 (code) of the component **Queries and selectors** (`geom.queries`) in the container **geom** (crate `geom`). [Component view of the container](../../3-components/geom.md) · [All code pages](../README.md)

> Exact queries on an evaluated model, serialisable selectors, the LabelMap between names and the integer labels solvers consume, flat 2D BVH with depth spill stack, and simulation operations as data.

**Decisions:** [ADR 2](../../../adr/0002-geometry-as-named-feature-list.md), [ADR 3](../../../adr/0003-geom-is-a-lean-separate-kernel.md), [ADR 7](../../../adr/0007-metric-adaptation-and-structured-pml.md), [ADR 18](../../../adr/0018-curved-elements-deferred.md), [ADR 132](../../../adr/0132-vvuq-asme-six-pillar-code-verification.md)  
**Size:** 5 files, 2976 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/geom/src/bvh.rs` | 246 |  |
| `crates/geom/src/labels.rs` | 254 | `LabelMap`: the compiled, deterministic mapping between names and the integer labels the solvers consume. Domain label = region (1-based in overlay order, pins honoured and checked); boundary label = one per facet (ordered by provenance,... |
| `crates/geom/src/query.rs` | 405 | Exact queries on an evaluated 2D model: measures, curve evaluation with first and second derivatives, curvature radii, chord / facet deviation, edge discretisation by deviation and size, midpoint projection (where a quadratic element put... |
| `crates/geom/src/select.rs` | 1145 | Selectors: small serialisable rules that resolve against an evaluated complex to typed entity sets (regions, facets, vertices). The text form (`name(core.top) \| interface(Si, BOX)`) parses to the serde form and prints back deterministica... |
| `crates/geom/src/sim.rs` | 926 | Simulation operations as *data*: absorbing-layer regions, mesh size rules, periodic face pairs (verified congruent), and the projection oracle. geom describes which entities and how; femcore decides what a condition means and applies it. |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Evaluation](./evaluation.md) | eval, scene | Evaluated, Node, PmlInfo, RegionX, axis_perm |
| [Feature list](./featurelist.md) | error, expr, model | Axis, Expr, GeomError, GeomResult, Op, PmlCoordKind, PmlGrowMode, PmlRole, RegionKind |
| [2D kernel](./kernel2d.md) | complex, curve, prim | Body, Complex, Curve, Edge, Face, FaceGeom, Region, curve_kind |
| [Math helpers](./math.md) | math | Aabb, V2 |

## Used by

| Component | Modules | Items |
|---|---|---|
| [engine: Geometry operations](../engine/geometry.md) | select | Selector |
| [engine: Meshing seam](../engine/meshing.md) | labels, select, sim | EntityKind, EntitySet, LabelMap, PmlRegion, Selector |
| [femcore: Fixtures and test kit](../femcore/fixtures.md) | select | Selector |
| [femcore: Meshing facade](../femcore/meshing.md) | labels, query, select | DiscretizeOpts, EntityKind, LabelMap, Selector |
| [femcore: 3D meshing](../femcore/meshing3d.md) | labels, query, select, sim | DiscretizeOpts, EntityKind, LabelMap, PeriodicPair |
| [Evaluation](./evaluation.md) | sim | SimData |
| [Feature list](./featurelist.md) | select | Selector |
| [2D kernel](./kernel2d.md) | bvh | AabbBvh2, Bbox2, ComplexBvh |
| [3D kernel](./kernel3d.md) | bvh | ComplexBvh |
| [Display tessellation](./tessellation.md) | query | DiscretizeOpts |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace bvh_rs {
    class Bbox2 {
        <<struct>>
        +lo V2
        +hi V2
        +new(lo, hi) Self
        +from_aabb_xy(b) Self
        +is_empty() bool
        +grow(p)
        +union(o) Bbox2
        +pad(d) Bbox2
        +contains_point(p, pad) bool
        +overlaps(o, pad) bool
        +more_2_methods()
    }
    class BvhNode2 {
        <<struct>>
        +bbox Bbox2
        +left_or_start u32
        +right u32
        +count u32
    }
    class AabbBvh2 {
        <<struct>>
        +nodes Vec~BvhNode2~
        +items Vec~u32~
        +item_boxes Vec~Bbox2~
        +build(item_boxes) Self
        +is_empty() bool
        +query_point(p, pad, visit)
        +query_box(target, pad, visit)
        +query_ray_x(p, pad, visit)
    }
    class ComplexBvh {
        <<struct>>
        +regions AabbBvh2
        +loops AabbBvh2
        +edges AabbBvh2
        +loop_edge_bvhs Vec
    }
  }
  namespace labels_rs {
    class DomainLabel {
        <<struct>>
        +id i32
        +name String
        +region u32
        +kind RegionKind
        +pinned bool
    }
    class BoundaryLabel {
        <<struct>>
        +id i32
        +name String
        +entities Vec~u32~
        +pinned bool
    }
    class Group {
        <<struct>>
        +name String
        +rule String
        +labels Vec~i32~
    }
    class LabelMap {
        <<struct>>
        +domains Vec~DomainLabel~
        +boundaries Vec~BoundaryLabel~
        +groups Vec~Group~
        +sheet_thickness_m BTreeMap
        +sheet_etch_depth_m BTreeMap
        +boundary(name) Vec~i32~
        +name_of_boundary(id) Option
        +domain(name) Option~i32~
        +name_of_domain(id) Option
        +edges_of_label(id) Vec~u32~
        +label_of_edge(e) Option~i32~
    }
    class Evaluated {
        <<impl block>>
        +label_map() GeomResult~LabelMap~
        +select_labels(sel) GeomResult
        +length(edge) f64
        +volume(region) f64
        +interface_area(a, b) f64
        +edge_eval(e, t) V3
        +edge_tangent(e, t) V3
        +edge_curvature(e, t) f64
        +more_20_methods()
    }
    class fns_labels_rs {
        <<module labels.rs>>
        +compile(ev) GeomResult~LabelMap~
    }
  }
  namespace query_rs {
    class DiscretizeOpts {
        <<struct>>
        +max_dev f64
        +max_len f64
    }
    class Fractions {
        <<struct>>
        +u f64
        +v f64
    }
    class Complex {
        <<impl block>>
        +discretize_edge(e, o) Vec
    }
  }
  namespace select_rs {
    class Cmp {
        <<enumeration>>
        Lt
        Le
        Gt
        Ge
    }
    class BoxMode {
        <<enumeration>>
        Within
        Touching
    }
    class SurfKind {
        <<enumeration>>
        Plane
        Cylinder
        Cone
        Sphere
        Torus
        Line
        Circle
        Ellipse
    }
    class NormalSpec {
        <<enumeration>>
        Axis
        Vector
    }
    class Selector {
        <<enumeration>>
        Name
        Tag
        Region
        Boundary
        Outer
        Interface
        From
        Kind
        Radius
        Area
        Plane
        Normal
        more_8_variants
        +parse(text) GeomResult~Selector~
        +atoms() usize
        +or(o) Selector
        +and(o) Selector
        +minus(o) Selector
    }
    class EntityKind {
        <<enumeration>>
        Regions
        Faces
        Edges
        Vertices
    }
    class EntitySet {
        <<struct>>
        +kind EntityKind
        +ids Vec~u32~
        +new(kind, ids) EntitySet
    }
    class Resolved {
        <<struct>>
        +set EntitySet
        +warnings Vec~String~
    }
  }
  namespace sim_rs {
    class PmlRegion {
        <<struct>>
        +region u32
        +parent Option~u32~
        +window Aabb
        +sides slice
        +thickness f64
        +layers f64
        +grading f64
        +coord PmlCoordKind
        +more_7_fields
    }
    class SizeRule {
        <<struct>>
        +feature FeatureId
        +entities EntitySet
        +h Option~f64~
        +max_dev Option~f64~
        +min_layers Option~f64~
    }
    class PeriodicPair {
        <<struct>>
        +feature FeatureId
        +a EntitySet
        +b EntitySet
        +translation V3
        +checked_congruent bool
    }
    class EntityRef {
        <<enumeration>>
        Region
        Face
        Edge
        Vertex
    }
    class SimData {
        <<struct>>
        +pml Vec~PmlRegion~
        +size Vec~SizeRule~
        +periodic Vec~PeriodicPair~
    }
    class Oracle {
        <<struct>>
        +project(p, e) Option~V3~
    }
  }
  BvhNode2 --> Bbox2
  AabbBvh2 --> BvhNode2
  AabbBvh2 --> Bbox2
  ComplexBvh --> AabbBvh2
  LabelMap --> DomainLabel
  LabelMap --> BoundaryLabel
  LabelMap --> Group
  Selector --> SurfKind
  Selector --> Cmp
  Selector --> NormalSpec
  Selector --> BoxMode
  EntitySet --> EntityKind
  Resolved --> EntitySet
  SizeRule --> EntitySet
  PeriodicPair --> EntitySet
  SimData --> PmlRegion
  SimData --> SizeRule
  SimData --> PeriodicPair
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `Bbox2` | `bvh.rs:12` | fields: `lo`, `hi` —  2D axis-aligned bounding box in the planar curve `(u, v)` frame. |
| struct | `BvhNode2` | `bvh.rs:97` | fields: `bbox`, `left_or_start`, `right`, `count` —  One node in a flat `Vec`-backed 2D AABB BVH. |
| struct | `AabbBvh2` | `bvh.rs:110` | fields: `nodes`, `items`, `item_boxes` —  Flat, deterministic 2D AABB bounding volume hierarchy. |
| struct | `ComplexBvh` | `bvh.rs:239` | fields: `regions`, `loops`, `edges`, `loop_edge_bvhs` —  Spatial acceleration hierarchy attached to a 2D [`crate::complex::Complex`] (in the 2D curve `(u, v)` frame). |
| fn | `Bbox2::new` | `bvh.rs:21` | `pub fn new(lo: V2, hi: V2) -> Self` |
| fn | `Bbox2::from_aabb_xy` | `bvh.rs:26` | `pub fn from_aabb_xy(b: &Aabb) -> Self` |
| fn | `Bbox2::is_empty` | `bvh.rs:35` | `pub fn is_empty(&self) -> bool` |
| fn | `Bbox2::grow` | `bvh.rs:40` | `pub fn grow(&mut self, p: V2)` |
| fn | `Bbox2::union` | `bvh.rs:48` | `pub fn union(&self, o: &Bbox2) -> Bbox2` |
| fn | `Bbox2::pad` | `bvh.rs:62` | `pub fn pad(&self, d: f64) -> Bbox2` |
| fn | `Bbox2::contains_point` | `bvh.rs:71` | `pub fn contains_point(&self, p: V2, pad: f64) -> bool` |
| fn | `Bbox2::overlaps` | `bvh.rs:76` | `pub fn overlaps(&self, o: &Bbox2, pad: f64) -> bool` |
| fn | `Bbox2::overlaps_ray_x` | `bvh.rs:85` | `pub fn overlaps_ray_x(&self, p: V2, pad: f64) -> bool` —  True if this box overlaps the horizontal ray `[p[0], +inf) x {p[1]}` within `pad`. |
| fn | `Bbox2::center` | `bvh.rs:90` | `pub fn center(&self) -> V2` |
| fn | `AabbBvh2::build` | `bvh.rs:118` | `pub fn build(item_boxes: Vec<Bbox2>) -> Self` —  Builds a deterministic top-down median-split 2D BVH over `item_boxes`. |
| fn | `AabbBvh2::is_empty` | `bvh.rs:165` | `pub fn is_empty(&self) -> bool` |
| fn | `AabbBvh2::query_point` | `bvh.rs:220` | `pub fn query_point(&self, p: V2, pad: f64, visit: impl FnMut(u32))` —  Visits all item indices whose 2D AABB (padded by `pad`) contains `p`. |
| fn | `AabbBvh2::query_box` | `bvh.rs:226` | `pub fn query_box(&self, target: &Bbox2, pad: f64, visit: impl FnMut(u32))` —  Visits all item indices whose 2D AABB overlaps `target` within `pad`. |
| fn | `AabbBvh2::query_ray_x` | `bvh.rs:232` | `pub fn query_ray_x(&self, p: V2, pad: f64, visit: impl FnMut(u32))` —  Visits all item indices whose 2D AABB overlaps the horizontal ray `[p[0], +inf) x {p[1]}` within `pad`. |
| struct | `DomainLabel` | `labels.rs:12` | fields: `id`, `name`, `region`, `kind`, `pinned` |
| struct | `BoundaryLabel` | `labels.rs:21` | fields: `id`, `name`, `entities`, `pinned` |
| struct | `Group` | `labels.rs:31` | fields: `name`, `rule`, `labels` —  A named group of boundary labels and the rule that produced it. |
| struct | `LabelMap` | `labels.rs:38` | fields: `domains`, `boundaries`, `groups`, `sheet_thickness_m`, `sheet_etch_depth_m` |
| fn | `compile` | `labels.rs:88` | `pub fn compile(ev: &Evaluated) -> GeomResult<LabelMap>` |
| fn | `LabelMap::boundary` | `labels.rs:52` | `pub fn boundary(&self, name: &str) -> Vec<i32>` —  All boundary labels whose name is `name` (a face split into fragments has several). |
| fn | `LabelMap::name_of_boundary` | `labels.rs:61` | `pub fn name_of_boundary(&self, id: i32) -> Option<&str>` |
| fn | `LabelMap::domain` | `labels.rs:64` | `pub fn domain(&self, name: &str) -> Option<i32>` |
| fn | `LabelMap::name_of_domain` | `labels.rs:67` | `pub fn name_of_domain(&self, id: i32) -> Option<&str>` |
| fn | `LabelMap::edges_of_label` | `labels.rs:70` | `pub fn edges_of_label(&self, id: i32) -> Vec<u32>` |
| fn | `LabelMap::label_of_edge` | `labels.rs:73` | `pub fn label_of_edge(&self, e: u32) -> Option<i32>` |
| fn | `Evaluated::label_map` | `labels.rs:234` | `pub fn label_map(&self) -> GeomResult<LabelMap>` |
| fn | `Evaluated::select_labels` | `labels.rs:239` | `pub fn select_labels(&self, sel: &Selector) -> GeomResult<Vec<i32>>` —  Integer labels (boundary labels for facets, domain labels for regions) a selector resolves to. |
| struct | `DiscretizeOpts` | `query.rs:22` | fields: `max_dev`, `max_len` |
| struct | `Fractions` | `query.rs:29` | fields: `u`, `v` —  Parametric position on an entity: one fraction along an edge (`v` is 0), two for a face. |
| fn | `Evaluated::length` | `query.rs:38` | `pub fn length(&self, edge: u32) -> f64` —  Length of an edge. |
| fn | `Evaluated::volume` | `query.rs:42` | `pub fn volume(&self, region: u32) -> f64` —  Area of a region (volume in 3D). |
| fn | `Evaluated::interface_area` | `query.rs:46` | `pub fn interface_area(&self, a: u32, b: u32) -> f64` —  Measure of the interface between two regions: length of the shared edges in 2D, area in 3D. |
| fn | `Evaluated::edge_eval` | `query.rs:54` | `pub fn edge_eval(&self, e: u32, t: f64) -> V3` |
| fn | `Evaluated::edge_tangent` | `query.rs:58` | `pub fn edge_tangent(&self, e: u32, t: f64) -> V3` |
| fn | `Evaluated::edge_curvature` | `query.rs:63` | `pub fn edge_curvature(&self, e: u32, t: f64) -> f64` —  Signed curvature along the canonical direction (positive: turns left). |
| fn | `Evaluated::edge_closest` | `query.rs:67` | `pub fn edge_closest(&self, e: u32, p: V3) -> (V3, f64)` —  Closest point of the edge (between its end points) and its parameter. |
| fn | `Evaluated::min_radius_of_curvature` | `query.rs:75` | `pub fn min_radius_of_curvature(&self, e: u32) -> f64` —  Smallest radius of curvature on the edge (infinite for a line); a mesher derives `n = pi / acos(1 - d / r)` chords from it. |
| fn | `Evaluated::chord_deviation` | `query.rs:94` | `pub fn chord_deviation(&self, e: u32, ta: f64, tb: f64) -> f64` —  Largest distance between the curve on `[ta, tb]` and its chord (exact for lines and circles, sampled for ellipses). |
| fn | `Evaluated::facet_deviation` | `query.rs:99` | `pub fn facet_deviation(&self, e: u32, pts: &[V3]) -> f64` —  Largest distance from the polyline through `pts` (straight facets) to the exact edge, sampled along each segment. |
| fn | `Evaluated::point_at_mid` | `query.rs:113` | `pub fn point_at_mid(&self, e: u32, ta: f64, tb: f64) -> (f64, V3)` —  The point of the exact edge nearest the midpoint of the chord between parameters `ta` and `tb` (the mid-edge node of a |
| fn | `Evaluated::project_midpoint` | `query.rs:128` | `pub fn project_midpoint(&self, e: u32, a: V3, b: V3) -> V3` —  The point of the exact edge nearest the chord midpoint of two positions on (or near) the edge; valid for chords up to a half turn. |
| fn | `Evaluated::discretize_edge` | `query.rs:134` | `pub fn discretize_edge(&self, e: u32, o: DiscretizeOpts) -> Vec<(f64, V3)>` —  Points `(t, p)` along the edge, end points included, such that the chord deviation is at most `max_dev` and each chord is at |
| fn | `Evaluated::fractions_of` | `query.rs:143` | `pub fn fractions_of(&self, ent: EntityRef, p: V3) -> Option<Fractions>` —  Fractions of the point of the entity closest to `p`. Edges: `u` in `[0, 1]` along `[t0, t1]`. Returns `None` for entities |
| fn | `Evaluated::position_at` | `query.rs:203` | `pub fn position_at(&self, ent: EntityRef, fr: Fractions) -> Option<V3>` —  The exact point of the entity at the given fractions. |
| fn | `Evaluated::topology_signature` | `query.rs:255` | `pub fn topology_signature(&self) -> u64` —  Hash of the combinatorial structure only: regions (kinds), edges (curve kinds, end vertices, sides, origins and names), |
| fn | `Evaluated::overlay_geometry_hash` | `query.rs:336` | `pub fn overlay_geometry_hash(&self) -> Option<u64>` —  Hash of the xy geometry of the overlay of a layer stack (curves and edge ranges): unchanged when only z extents move. |
| fn | `Complex::discretize_edge` | `query.rs:352` | `pub fn discretize_edge(&self, e: u32, o: DiscretizeOpts) -> Vec<(f64, V3)>` |
| enum | `Cmp` | `select.rs:19` | variants: `Lt`, `Le`, `Gt`, `Ge` |
| enum | `BoxMode` | `select.rs:27` | variants: `Within`, `Touching` |
| enum | `SurfKind` | `select.rs:33` | variants: `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`, `Line`, `Circle`, `Ellipse` |
| enum | `NormalSpec` | `select.rs:45` | variants: `Axis`, `Vector` |
| enum | `Selector` | `select.rs:55` | variants: `Name`, `Tag`, `Region`, `Boundary`, `Outer`, `Interface`, `From`, `Kind`, `Radius`, `Area`, `Plane`, `Normal`, ... |
| enum | `EntityKind` | `select.rs:80` | variants: `Regions`, `Faces`, `Edges`, `Vertices` |
| struct | `EntitySet` | `select.rs:88` | fields: `kind`, `ids` |
| struct | `Resolved` | `select.rs:103` | fields: `set`, `warnings` |
| fn | `EntitySet::new` | `select.rs:95` | `pub fn new(kind: EntityKind, mut ids: Vec<u32>) -> EntitySet` |
| fn | `Selector::parse` | `select.rs:495` | `pub fn parse(text: &str) -> GeomResult<Selector>` |
| fn | `Selector::atoms` | `select.rs:498` | `pub fn atoms(&self) -> usize` |
| fn | `Selector::or` | `select.rs:505` | `pub fn or(self, o: Selector) -> Selector` |
| fn | `Selector::and` | `select.rs:508` | `pub fn and(self, o: Selector) -> Selector` |
| fn | `Selector::minus` | `select.rs:511` | `pub fn minus(self, o: Selector) -> Selector` |
| fn | `Evaluated::resolve` | `select.rs:993` | `pub fn resolve(&self, sel: &Selector) -> GeomResult<Resolved>` —  Resolves a selector to an entity set plus warnings (an empty result is a warning, not an error). |
| fn | `Evaluated::resolve_inner` | `select.rs:997` | `pub(crate) fn resolve_inner(&self, sel: &Selector, no_id: bool) -> GeomResult<Resolved>` |
| fn | `Evaluated::suggest` | `select.rs:1008` | `pub fn suggest(&self, picked: &EntitySet) -> Vec<Selector>` —  Selectors, simplest first, that each resolve to exactly `picked`; empty if none does. |
| struct | `PmlRegion` | `sim.rs:13` | fields: `region`, `parent`, `window`, `sides`, `thickness`, `layers`, `grading`, `coord`, `mode`, `center`, `axis`, `r_in`, ... |
| struct | `SizeRule` | `sim.rs:34` | fields: `feature`, `entities`, `h`, `max_dev`, `min_layers` |
| struct | `PeriodicPair` | `sim.rs:43` | fields: `feature`, `a`, `b`, `translation`, `checked_congruent` |
| enum | `EntityRef` | `sim.rs:52` | variants: `Region`, `Face`, `Edge`, `Vertex` |
| struct | `SimData` | `sim.rs:60` | fields: `pml`, `size`, `periodic` |
| struct | `Oracle` | `sim.rs:867` |  |
| fn | `Evaluated::pml_regions` | `sim.rs:67` | `pub fn pml_regions(&self) -> &[PmlRegion]` |
| fn | `Evaluated::size_rules` | `sim.rs:70` | `pub fn size_rules(&self) -> &[SizeRule]` |
| fn | `Evaluated::periodic_pairs` | `sim.rs:73` | `pub fn periodic_pairs(&self) -> &[PeriodicPair]` |
| fn | `Evaluated::compute_sim` | `sim.rs:77` | `pub(crate) fn compute_sim(&mut self) -> GeomResult<()>` |
| fn | `Evaluated::oracle` | `sim.rs:210` | `pub fn oracle(&self) -> Oracle<'_>` —  Projection oracle: closest point on an exact entity (shaped to be wrapped by femcore's `Geometry2::project`). |
| fn | `Evaluated::validate_pml_ray_invariance` | `sim.rs:632` | `pub fn validate_pml_ray_invariance(&self) -> GeomResult<()>` —  Validates Teixeira-Chew / Collino-Monk analytic-continuation ray invariance (`d_u eps == 0`) |
| fn | `Oracle::project` | `sim.rs:873` | `pub fn project(&self, p: V3, e: EntityRef) -> Option<V3>` —  Closest point of the entity to `p`; `None` for an unknown id or an entity kind the model does not have. |

Trait implementations in this component: `Selector: Display`
