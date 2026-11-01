# geom: handoff brief for the session that builds the geometry crate

Audience: the session that will write `crates/geom` (separate from the nyx integration session). Read this file, then
`docs/GEOMETRY_ARCHITECTURE.md` (the design and the reasons), then the code pointers in section 12.
Written 2026-10-31 against nyx main `54b8829`. The architecture doc says *what and why*; this file says *what to deliver, to what contract, and how it will be checked*.
Names in Rust sketches are proposals; the semantics, ids, ordering and error behaviour are the contract.

---

## 0. Decisions already made by Paul (do not reopen)

1. `geom` is a **separate crate**, written by this session, integrated into nyx later by the integration session.
2. **Import is by mesh only** (femcore already accepts cell lists). No STEP, no NURBS, no `Import` leaf, no `geom-step`. This supersedes section 7.6
   and Phase 5 of the architecture doc.
3. **Priorities**: efficiency, speed, accuracy; **lean** (roughly 6k lines of Rust, `serde` the only required dependency).
4. **No manual numbering**: names and rule-based selectors are the interface; ints are a compiled output. Everything stays reachable from code and GUI.
5. **Code-first and GUI are equivalent**: the serialisable feature list is the single source of truth.
6. The **curved-element refactor of the solvers is not this session's job**. The integration session does it after receiving the code and benchmarks it
   against closed forms. This session's job is to make that refactor possible (section 6).

---

## 1. Scope

In scope: feature list + parameters + expressions, evaluation tiers and caching, analytic B-rep (planes, cylinders, cones, spheres, tori; lines,
circles, ellipses, procedural surface-surface curves), the n-ary arrangement with CSG labelling into a region complex, provenance and names,
selectors, `LabelMap`, projection/closest-point oracle, edge discretisation by deviation and size, layer-stack detection, display tessellation,
builder API, serde round trip, tests, `geom_bench`.

Out of scope (do not build): meshing of any kind (no tetrahedra, no triangulation for FEM; display tessellation only), STEP/IGES/NURBS/import,
fillets/chamfers/shells/offsets in 3D, sketch constraint solver, assemblies, healing, any dependency on `femcore`, `faer`, the GUI or physics, boundary-condition
*types* (geom resolves *which* entities, femcore decides what a condition means), solver changes.

---

## 2. Crate conventions (match the repository)

- Location `crates/geom`, add to `[workspace] members`. `edition = "2021"`, `rust-version = "1.84"` (workspace values).
- Dependencies: `serde` (derive). Optional feature `parallel` = `rayon` (off by default). Dev-dependencies only: `serde_json`, `proptest`.
  No `femcore` dependency, ever (dependency direction is `app -> femcore -> geom`).
- Units: SI metres, `f64`. The model may be nanometre-sized at micrometre offsets; all tolerances are relative to the model size (section 7.1).
- Every source file starts with a short `//!` explaining what it does (femcore convention).
- `cargo clippy --release --workspace --all-targets -- -D warnings` must be clean (the repo keeps it clean).
- Errors: `Result<_, GeomError>` where `GeomError: Display` and `Into<String>`; femcore's APIs use `Result<_, String>`. Messages are plain language and name
  the feature(s) involved ("`core.side` touches `clad.side` tangentially ..."). NaN or infinite input is refused with a message that names the input
  (repo rule, see README "NaN or infinite input ... is refused").
- **Determinism**: same input gives bit-identical output, independent of thread count (repo guarantee; `tests/determinism.rs` exists in femcore). No
  `HashMap` iteration order may leak into ids or output; sort or use `BTreeMap`/indexed vectors.
- **No panics on user input**: a random feature list returns `Err`, never panics (property test, section 9).
- Serialisation must round-trip through `serde_json` with `float_roundtrip` **and** through `femcore::xmlio` (it writes any serde value as XML:
  objects with a `t` attribute, keys that are not valid XML names become `<entry key=..>`, arrays as `<item>`). Practical consequences: prefer enums
  and structs over maps with non-string keys; a number-or-expression field is an `untagged` enum of `f64` and `String`; no `u128`, no NaN literals.
- Windows is the development machine; no platform-specific code or paths.

---

## 3. Data model (summary; reasons in the architecture doc)

```text
FeatureList { params: Vec<Param>, features: Vec<Feature> }
Feature     { id: FeatureId(u32, stable), name: String, op: Op, suppressed: bool, origin: Gui | Code, tags: Vec<String> }
Param       { name, value: Expr }                         Expr = f64 | String      (reuse the grammar of femcore/src/expr.rs, minus x y z)
Op          Rect|Circle|Ellipse|Polygon|RoundedRect (2D profiles) | Box|Cylinder|Cone|Sphere|Torus
            | Extrude{profile, z0, z1 | layers, taper?} | Revolve{profile, axis, angles}
            | Transform{input, translate/rotate/mirror/uniform scale} | Pattern{input, linear|polar, n}
            | Union[ids] | Difference{a, [b]} | Intersect[ids]
            | Overlay{regions: Vec<RegionDef>}            (root; list order is painter priority, later wins)
            | Pad/Pml{target, sides, thickness, layers, grading} | Cut{plane, keep} | Periodic{sel_a, sel_b, translation}
            | Size{sel, h?, max_deviation?, min_layers?} | Tag{sel, name}
RegionDef   { name, input: FeatureId, kind: Material | Hole | Pml | Void, pinned_id: Option<i32> }
```

References point to earlier features only (flat list = DAG in topological order; no cycles possible). Normalisation (flatten n-ary Boolean nodes,
drop empties by AABB, deduplicate identical operands via hash-consing, `A cap A = A`) runs before evaluation.

Evaluation tiers (each cheaper tier never builds a costlier one): **Bounds** (AABB per node, params only), **Classify** (inside/sdf-bound from the tree, no
B-rep), **Complex** (arrangement, provenance, selectors), **Display** (tessellation with face ids). Per-node results are cached under a content hash of the
op parameters (after expression evaluation) and the input hashes. A parameter change must recompute only the dependent nodes (tested).

---

## 4. Behaviour spec details

### 4.1 Roles and auto-names

| Primitive | Roles (fixed, documented, never renamed) |
|---|---|
| Rect (2D profile) | `left bottom top right` (the existing mesher's order; keep for legacy lifting) |
| Circle / Ellipse | `edge` |
| Polygon | `e0 e1 ...` in vertex order, or user-supplied names |
| RoundedRect | `left bottom top right` plus `c_bl c_br c_tl c_tr` for the corner arcs |
| Box | `xmin xmax ymin ymax zmin zmax` |
| Cylinder / Cone | `side bottom top` |
| Sphere / Torus | `surface` |
| Extrude(profile) | `side:<profile role>`, `bottom`, `top` (`bottom` at the lower z) |
| Revolve(profile) | `swept:<profile role>`, `start`, `end` (the last two absent for a full turn) |

Transforms and patterns keep roles; pattern instances add `[i]` to the feature name (`hole[3].side`).
Auto-name of a boundary face: `<feature>.<role>`; of an interface between regions: `A|B` with the two **region** names in sorted order
(`BOX|Si`, so the name does not depend on face orientation). Fragments of one original face share the name. Edges and vertices follow the same
scheme (edge name = `<face name>/<face name>` for the two faces meeting at it, sorted).

### 4.2 Overlay semantics

Regions are ordered; where several contain a point, the **last** wins (this is the current painter rule). `Hole` regions remove cells (the current
"domain 0 or negative"). `Void` is the implicit outside (id 0 / no region). A `Pml` kind region is a normal region that femcore treats specially.
Interfaces exist between any two different regions (including Void); an interface face has two sides and belongs to the complex once.

### 4.3 Selector grammar (text form must parse to the serde form and back, deterministic printing)

```text
sel   := term (('|' | '&' | '-') term)*   | '!' term | '(' sel ')'
term  := name(<auto-name>)  | tag(<tag>)  | region(<name>)  | boundary(sel_region) | outer | interface(<A>,<B>) | interface(<A>,*)
       | from(<feature>,<role>) | kind(plane|cylinder|cone|sphere|torus|line|circle|ellipse) | radius(<expr>[,<tol>]) | area(<cmp> <expr>)
       | plane(x|y|z = <expr>[,<tol>]) | normal(<±x|±y|±z|vector>[,<angle tol>]) | extreme(<±axis>) | in_box(lo, hi [, within|touching])
       | adjacent(sel) | id(<int>)          # id(...) is the legacy atom: addresses a compiled int label
```

Results are typed (`Faces | Edges | Regions | Vertices`); mixing types in `|`, `&`, `-` is an error with a clear message. `|`, `&`, `-` bind left to right,
`!` complements within the entity type of its operand. Names are matched case-sensitively. A selector that resolves to nothing returns an **empty set plus a
warning object**, not an error (the GUI shows it red; the integration session decides what a physics layer does).

Also required: `suggest(&Evaluated, picked: &EntitySet) -> Vec<Selector>`, returning candidates ordered by simplicity that each resolve to *exactly* the
picked set (auto-name, tag, name+predicate, geometric predicate, box), or an empty list if none does.

### 4.4 LabelMap (compile step)

`compile(&Evaluated) -> LabelMap` mints ints **deterministically**: domain label = region id (1-based, Overlay order; pinned ids honoured and checked for
collisions), boundary label = face id (1-based, ordered by provenance key then centroid; 2D: edges). Output holds `name <-> int` both ways, the
entity count per label, and the rule that produced a named group. A selector resolves to `Vec<i32>` of boundary or domain labels via the map.
Legacy projects: `lift_legacy(Geometry)` (done by the integration session, but geom must expose `pinned_id` to make it possible) maps old numbering to pins:
2D rect edges `1..4 = left,bottom,top,right` continued across shapes in order; polygon edges each their own number; circle one; 3D
`Solid::Box` 1..6 as its doc says (xmin, xmax, ymin, ymax, zmin, zmax; note the docs of `SimplexMesh::block` disagree, see section 8.3), cylinder `1 side 2 base 3 top`, sphere 1, extrusion `edges..., bottom, top`.

### 4.5 Simulation ops: geom describes, femcore applies

`Periodic`, `Pml/Pad`, `Cut`, `Size` produce **data** in the evaluated model, not meshes:
`PeriodicPair { face_a, face_b, translation, checked_congruent: bool }`, `PmlRegion { region, box, sides, thickness, layers, grading }`, `SizeRule { entities, h, max_dev, min_layers }`,
and `Cut` actually changes the geometry (keeps one side of a plane, adds the cut face named `<feature>.cut`, tagged `symmetry`).
`Periodic` must verify the two faces are congruent under the translation (same area, same boundary curves up to `tol`) and fail with a named error if not.

### 4.6 Layer-stack detection

`layer_stack(&Evaluated) -> Option<LayerStack>` returns `Some` when the 3D model is a stack along one axis of prisms over 2D partitions:
`LayerStack { axis, layers: Vec<{z0, z1, partition: Complex2D, region_of: Vec<RegionId>}>, overlay: Complex2D }` where `overlay` is the **common refinement** of
all layers' 2D partitions (every layer boundary is an edge of it) and each overlay face maps to one region per layer. Must also report the
face/edge provenance mapping so side faces of the extruded mesh can be named from overlay edges and caps from layer tops/bottoms.
Typical cases that must be detected: SOI/rib strips, rib on sapphire with PML boxes, the Bragg cell (layers with and without teeth), holes through a slab.
Not detected (returns `None`, general 3D path): tapered sweeps, spheres, tori, rotated primitives.

---

## 5. API contract needed by femcore (sketch)

```rust
// build / edit
FeatureList::new(); builder: Model::new_2d()/new_3d(); .param(..); .rect(..) ... returns FeatureId; list.to_rust() -> String;
list.with(&[("w", 450e-9)]) -> FeatureList                      // parameter override, no mutation
// evaluate
Evaluated = model.eval(Tier::Complex)?;                           // cached; Tier::{Bounds, Classify, Complex, Display}
eval.bounds() -> Aabb;     eval.classify(p) -> Option<RegionId>;  eval.sdf_bound(p, margin) -> f64
// topology and geometry (ids are dense u32, deterministic)
eval.regions() / faces() / edges() / vertices();  face.surface(), face.loops(), face.sides() -> (Option<RegionId>, Option<RegionId>), face.name(), face.provenance()
surface.eval(u,v), .normal(p), .closest(p) -> (point, uv), .curvatures(p) -> (k1, k2);   curve.eval(t), .tangent(t), .curvature(t), .closest(p) -> (point, t)
eval.area(face), eval.length(edge), eval.volume(region), eval.interface_area(a, b)           // analytic
// meshing support (no meshing inside geom)
eval.discretize_edge(edge, DiscretizeOpts{max_dev, max_len}) -> Vec<(f64 /*t*/, [f64;3])>   // deterministic; both neighbouring faces use the same cached points
eval.surfaces_for_stuffing() -> Vec<{surface, face_aabb, faces}>                              // plane/cyl/cone/sphere/torus + local extent
eval.layer_stack() -> Option<LayerStack>
eval.size_rules() / periodic_pairs() / pml_regions()
oracle = eval.oracle();  oracle.project(p, EntityRef::{Face(id)|Edge(id)|Vertex(id)}) -> [f64;3]
// naming
Selector (serde) + Selector::parse(&str)/to_string(); eval.resolve(&sel) -> EntitySet (+ warnings); eval.suggest(&picked) -> Vec<Selector>
eval.label_map() -> LabelMap (serde);   label_map.boundary("Si|BOX") -> Vec<i32>;  .name_of_boundary(i32); .domain("Si") -> i32
// display
eval.tessellate(max_dev) -> DisplayMesh { points, tris, tri_face: Vec<u32>, segs, seg_edge: Vec<u32> }   // for viewport and picking
```

`oracle.project` is shaped to be wrapped by femcore's existing `meshers::adapt::geom::Geometry2::project(p, chain)` (analytic projector per boundary label).

---

## 6. Preparing for the curved-element refactor (what geom must expose; the refactor itself is not yours)

Why it is needed: femcore cells are straight-sided and constant-Jacobian (`fem::geometry` returns one `vol` and one barycentric-gradient set per
cell, called from about 20 sites in about 15 files), so curved boundaries carry an O(h^2) geometric error whatever the element order. The integration session will add a
per-cell geometry provider (`Affine | Curved`) and gradients per quadrature point, with the covariant Piola map for the Nedelec edge elements, **using
only the following from geom**. Provide all of them, tested:

1. **Exact closest point and normal** on every face and edge (already in the API) with **first and second derivatives**: `curve.eval/tangent/curvature`,
   `surface.normal/curvatures`. A quadratic geometry element needs the on-curve midpoint of an edge and a consistent normal.
2. **`curve.point_at_mid(t_a, t_b)`** and **`surface.project_midpoint(a, b)`** helpers: given the parameters (or positions) of two nodes on an edge/face, the point of the
   exact entity nearest their chord midpoint (this is where the P2 mid-edge node goes). Must be monotone and stable for periodic parameters (circle seam, cylinder theta seam).
3. **Association-friendly ids**: `FaceId/EdgeId/VertexId` dense, deterministic, stable under unrelated edits of other features. femcore builds
   `facet_face: Vec<u32>` and `node_on: Vec<(u8,u32)>` from them; geom does not know meshes.
4. **Deviation metrics for verification**: `chord_deviation(edge, t_a, t_b)` and `facet_deviation(face, p0, p1, p2)` (max distance of a straight edge/triangle
   from the exact entity). Used to set the sizing rule `max_dev` and to assert in tests that deviation-driven discretisation hits its target.
5. **Orientation**: every face has an outward normal defined by its region side convention (document it: normal points from the first listed region to the second).
6. **Curvature radius queries for sizing**: `min_radius_of_curvature(face|edge)` so a mesher can derive `n = pi / acos(1 - d/r)` chords.
7. **Closed-form geometric oracles** (used as benchmark ground truth later): exact area/perimeter/volume of every primitive and composed region, `analytic` (planar
   polygons, circle arcs, cylinders, spheres, tori), not tessellated. Property to test: sum of region volumes equals union volume to `1e-12` relative.

What the integration session will do afterwards, so you can see how your outputs are used: (a) a driver per mesher route that discretises edges once, meshes, then
projects boundary nodes onto exact entities and rejects inverted cells; (b) deviation-driven sizing; (c) benchmark of straight-cell vs deviation-sized vs
curved-P2 elements against the closed forms in section 8.4; (d) only if (c) shows a need, the solver-wide curved-cell refactor.

---

## 7. Robustness and numerical rules

### 7.1 Tolerance model
`L` = bounding-box diagonal of the model, `tol = 1e-9 * L`, one value per model, passed in a context, not a global. Angles use `1e-10` rad for "same orientation".
Work in coordinates centred and scaled by `L` internally if it helps; results are returned in model units.

### 7.2 Interning (core of coincidence handling)
Planes (unit normal sign-canonical, offset), cylinders (axis line, radius), cones, spheres, tori are hash-consed to a `SurfId` within `tol`; edges are interned by (surface pair, parameter range).
Two coincident faces have equal `SurfId`; "coplanar" is an integer compare. For same-surface overlap decide sides from normal orientation, not by sampling
an offset point that would straddle the other primitive. This is what makes `conform_chains` (femcore `mesher.rs`, ~140 lines of tolerance repair) unnecessary.

### 7.3 Predicates
2D orientation/in-circle: exact predicates (`robust` crate is already in `Cargo.lock` via `spade`; it may be a dependency of geom). 3D: tolerance-based on analytic parameters.

### 7.4 Refuse, don't guess
Grazing contact (distance in `(0, k*tol]` with `k` a small documented constant, or surfaces tangent within angle tolerance but not coincident) is an error
naming both features and what to change. The same for zero-thickness slivers (a region thinner than `tol`), self-overlapping profiles, and too many
fragments (a documented cap with a message, replacing silent blow-up).

### 7.5 Always-on self-checks (debug and tests)
Each edge used by exactly two face sides (one on the outer boundary); loops closed; Euler relation per shell; positive region volumes; sum of region volumes =
union volume; interface area equal from both sides.

---

## 8. Acceptance material

### 8.1 Golden models (deliver as JSON `FeatureList` files in `crates/geom/golden/`, each with an `expected.json` of geometric invariants)

Each golden file must evaluate to the listed regions and face-name set, with areas/volumes matching closed forms to `1e-10` relative. The integration session
builds these into meshes and compares **solver outputs** against the current code (so geom does not need femcore for that).

| Model | Source in nyx | Facts to reproduce |
|---|---|---|
| SOI strip, oxide cladding | `usecases/soi.rs` `SoiWire::default()` | three regions BOX/Si/clad (current domains 1, 2, 3), 500x300 nm core on 2 um BOX, 1 um margin |
| SOI with rib slab, sloped walls, left wall 75 deg | `SoiWire` fields `slab`, `sidewall_deg`, `sidewall_left_deg` | polygon core; breaks mirror symmetry |
| Four ring-chip cross-sections incl. PML regions | `chips.rs` (`chip_shapes`, domains 1-4, 11-15, 21-22, 41-45) | the PML combination domains come from `Pml` features, no hand numbers |
| Bragg cell: symmetric, asymmetric, uniform | `usecases/bragg.rs` | layer stack (teeth present/absent), identical end faces, `Periodic` pair |
| Plate with a circular hole | `tests/mesher.rs` `plate_with_a_hole_...` | Kirsch stress concentration (factor 3) |
| Meshed disk | `tests/mesher.rs` `poisson_on_a_meshed_disk_...` | exact area `pi r^2` |
| Annulus, concentric cylinders, concentric spheres | new | capacitance/Lame closed forms (8.4) |
| Sphere in a box (2 regions), cylinder cavity, torus | new | volumes `4/3 pi r^3`, `pi r^2 h`, `2 pi^2 R r^2` |
| Box with a slanted plane cut; two skew cylinders | new | exercises `SurfSurf` and refusal path |

### 8.2 Existing reference numbers the integration session will hold the new path to
SOI 300 nm strip at default sizes: 1312 cells, 9077 unknowns, TE neff 2.67179 (TE share 0.989), TM neff 2.33374 (docs/photonic_use_cases.md);
Bragg coarse mesh: 6696 cells, 44,838 unknowns, symmetric band edges 194.680 and 195.415 THz, asymmetric 196.236 and 196.644 THz; chip numbers in `docs/chip_report.md`.
Mesh counts will differ with a new mesher driver, so the gate is mesh-converged physical results within the stated mesh sensitivity, not identical meshes.

### 8.3 Known inconsistency to leave as a test, not to "fix" silently
`mesher3d.rs` documents `Solid::Box` boundary numbering as "the numbering of `SimplexMesh::block`" (1 x=min, 2 x=max, 3 y=min, 4 y=max, 5 z=min, 6 z=max), while
`mesh.rs` documents `block` as 1 x=0, 2 y=0, 3 z=0, 4 z=lz, 5 y=ly, 6 x=lx. Roles replace both; legacy pins follow the **code's actual behaviour**, which the integration
session will check by running the existing tests, not geom.

### 8.4 Closed-form problems for the later benchmark (geometry supplied by geom; solves done by the integration session). Verify each constant before use.

| Physics | Geometry | Closed form |
|---|---|---|
| 2D EM, PEC circular guide | disc radius R | cutoffs `k_c R`: TE11 1.841184, TM01 2.404826, TE21 3.054237, TM11 3.831706 |
| 2D EM, step-index fibre | core in cladding, circular | exact characteristic equation (Bessel/Hankel) |
| 3D EM, PEC cylinder cavity | cylinder R, h | `f_010 = 2.404826 c / (2 pi R)`, independent of h |
| 3D EM, PEC spherical cavity | sphere R | lowest TM: `kR = 2.7437`, lowest TE: `kR = 4.4934` |
| Electrostatics / heat | concentric cylinders / spheres | `C = 2 pi eps L / ln(b/a)`, `C = 4 pi eps a b / (b - a)` |
| Solid, plane strain | thick-walled cylinder under pressure | Lame solution |
| Solid | plate with hole | Kirsch (exists as a test) |

Protocol (integration session): refine `h` by 2 repeatedly; report error and observed order for (i) straight cells as today, (ii) deviation-driven sizing, (iii) curved cells.
Hypothesis to test, not a result: straight cells saturate near O(h^2) for P2 on curved boundaries; curved P2 restores a higher order.

---

## 9. Test plan for geom

- Unit and closed-form: areas/volumes/perimeters of all primitives and compositions; intersection curves of every supported surface pair against closed forms.
- Invariants (7.5) after every operation in randomised trees.
- Property tests (`proptest`): union volume associativity and commutativity; `classify` agrees with the B-rep (point inside region iff inside the shell) on random points;
  `A cap A = A`; name stability (adding an unrelated feature or reordering independent features does not change any existing auto-name); random feature lists never panic.
- Selector tests: parse/print round trip; every example in section 4.3; `suggest` returns only selectors that resolve exactly to the picked set; empty-result warnings.
- Serde: JSON round trip equality (bit-exact floats); XML-compatibility check by a test that writes and reads the value through the same quick-xml approach (a small local copy of the encoding rules, or a dev-dependency on a shared helper; the integration session will add the real `xmlio` test).
- Builder equivalence: `build(list.to_rust())` equals `list` (a generated test that compiles the printed code is optional; a pure data check is required).
- Determinism and thread independence (with `parallel` on).
- Tolerance edge cases: nm models at 1e-6 offsets, coincident faces, T-junctions, grazing contact refused with a named error.
- Cache: changing one parameter recomputes only dependent nodes (count evaluations).

## 10. Performance (targets to verify, none measured yet)
`examples/geom_bench.rs`: evaluate a ~50-feature model (target under about 5 ms), resolve selectors over ~500 faces (under about 1 ms), parameter sweep of 100 values showing cache hits,
a 200-hole pattern for pair-finding cost. If pair-finding dominates, add a face-AABB sort-and-sweep or grid (not a BVH first). Record numbers in the crate README.

## 11. Suggested order of work (early integration value first)

| Step | Deliver | Done when |
|---|---|---|
| A | Surfaces/curves + 2D profiles + 2D arrangement + overlay + names + selectors + `LabelMap` + serde + builder | 2D goldens (SOI, chips, plate, disc) pass invariants and closed forms; femcore can already start the 2D driver |
| B | Params/expressions, tiers, hash cache, normalisation, `to_rust`, tessellation, oracle, deviation metrics, curvature queries (section 6) | sweep and cache tests pass; section 6 items tested |
| C | 3D primitives, sweeps, transforms, patterns, 3D arrangement (all five surfaces, `SurfSurf`), `Cut`, `Periodic`, `Pml`, `Size` data | 3D goldens pass; refusal tests pass |
| D | Layer-stack detection, `surfaces_for_stuffing` with face extents (replaces the 64-surface cap on the mesher side), `geom_bench` | Bragg cell detected as a stack; bench recorded |

Deliver after each step with the golden files and a CHANGELOG line, so the integration session can start on A while B-D are written.

## 12. Code pointers (read these in nyx)

| What | Where |
|---|---|
| Current 2D mesher, conforming hack, label rules | `crates/femcore/src/mesher.rs` (`Shape`, `Part`, `conform_chains`, end of `mesh_2d`) |
| Current 3D mesher, surfaces, creases, cutting surfaces | `crates/femcore/src/mesher3d.rs` (`Solid`, `Surf`, `FaceSurf`, `mesh_3d`) |
| The mesh contract with solvers | `crates/femcore/src/mesh.rs` (`SimplexMesh`: `cell_labels`, `facet_labels`, `relabel`) |
| Where ints are consumed | `crates/femcore/src/problem.rs` (`Bc.labels`, `Materials` per-domain ints, `Periodic`, `inactive`, `optical_pml`) |
| Geometry oracle seam | `crates/femcore/src/meshers/adapt/geom.rs` (`Geometry2::project`) |
| Structured PML builder | `crates/femcore/src/meshers/adapt/pml.rs`, `meshers/interop.rs` |
| Use cases to reproduce | `crates/femcore/src/usecases/{soi,bragg}.rs`, `chips.rs`, `docs/photonic_use_cases.md`, `docs/chip_report.md` |
| Expression grammar to mirror | `crates/femcore/src/expr.rs` |
| Project save format | `crates/femcore/src/xmlio.rs`, `crates/app/src/state.rs` (`Geometry`, `GeometrySpec`, `Region`) |
| Element geometry (curved-cell target) | `crates/femcore/src/fem.rs` (`geometry`, `quad_rule`) |

## 13. Integration plan after delivery (for the other session; recorded here so both agree)

1. Add `geom` as a `femcore` dependency; `lift_legacy()` for `Shapes`/`Solids`/`Rectangle`/`Block`.
2. Drivers: 2D (edge discretisation once, `mesh_2d` core reused), layer stack (overlay + `extrude`), 3D stuffing fed by `surfaces_for_stuffing`; association arrays; `LabelMap` in results.
3. Name-keyed `Bc`/`Materials`/`Periodic` resolved by one `resolve(&mut ProblemSpec, &LabelMap)` before `solve`. Check that `Bc.labels` matching is not a linear scan per facet.
4. Golden comparisons (8.1, 8.2). 5. Closed-form benchmark and the curved-cell decision (section 6, 8.4). 6. GUI feature-list panel and pick-to-selector.

## 14. Open items and risks to report back on
- Arrangement robustness for coincident and grazing faces is the main risk; report any configuration the interning cannot decide, with a minimal feature list.
- Fragment-name convention (a split face keeps one name; a single fragment needs a predicate): confirm it feels right once the first goldens exist.
- If total size passes about 8k lines before step C is done, stop and report which module is growing and why (the lean requirement is part of the acceptance).
- Anything in this file that conflicts with the code you read: trust the code, note the conflict in your CHANGELOG.
