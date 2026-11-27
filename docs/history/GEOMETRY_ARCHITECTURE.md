# nyx geometry: proposal for a lean CSG-over-B-rep kernel

Update 2026-10-31: Paul decided that mesh import is enough, so section 7.6 (STEP) and Phase 5 are **dropped**; and the `geom` crate will be written by a
separate session from `docs/GEOM_HANDOFF.md`, with the curved-element refactor done later by the integration session. Where the two files differ, the handoff wins.

Status: design proposal, no code written. Based on a read of `femcore` at main `54b8829` (`mesher.rs`, `mesher3d.rs`, `mesh.rs`,
`problem.rs`, `usecases/*`, `chips.rs`, `meshers/adapt/*`, `app/src/state.rs`, `app/src/tabs.rs`) and of the six reference systems.
Nothing here has been benchmarked; every number marked "target" is a goal to be checked in Phase 1, not a measurement.

---

## 0. Summary of recommendations

1. **One geometry language for 2D and 3D**: a flat, serialisable *feature list* whose entries reference earlier entries (a DAG by
   reference, section 5). Leaves are primitives and profile sweeps, inner nodes are transforms, patterns and Booleans, and the root is an
   ordered, named **Overlay** of regions. This replaces `Rectangle`, `Block`, `Shapes`, `Solids`, `extrude` and the bbox `relabel` with one thing.
2. **Evaluate it into a *region complex*, not a single solid.** FEM needs a partition of space into named regions with *shared* interface faces
   (one face, two sides), not a pile of overlapping solids. The central algorithm is one n-ary **arrangement** (common refinement of all
   primitive boundaries) plus a CSG expression evaluated on the per-primitive inside bits of every fragment. Union, difference and intersection
   are the same algorithm, so there is only one Boolean to get right.
3. **B-rep with analytic geometry only**: planes, cylinders, cones, spheres, tori; lines, circles, ellipses and a procedural surface-surface
   curve. No NURBS, no fillets, no STEP in v1. That covers every geometry in the repo today and the ones we are likely to add (strips, ribs,
   sloped walls, discs, rings, holes, photonic-crystal slabs, spheres) and keeps the kernel at roughly 6k lines.
4. **"Exact" means**: curved boundaries stay analytic until the mesher asks; topology decisions are integer comparisons on *interned* geometry
   (coincident faces share one surface id), not float fuzz; every mesh boundary node is projected onto its exact face/edge. It does *not* mean
   arbitrary-precision arithmetic. Near-tangent input is refused with a message naming the two features, not meshed into garbage.
5. **Names, not numbers.** Every face, edge and region descends from exactly one primitive role (`core.top`, `Si|BOX`), so persistent naming
   falls out of provenance. Boundary conditions, materials, mesh sizes, periodic pairs and PML are attached through **selectors** (a small rule
   language: names, tags, interface/outer, plane/normal/extreme, boxes, set algebra). Ints still exist, but as a *compiled output*
   (`LabelMap`), addressable from code and GUI, never typed by the user.
6. **Code and GUI edit the same data.** The feature list *is* the model; the Rust builder only constructs it, the GUI form edits it, and "Copy
   as Rust" prints it. Parameters are named and shared, so sweeps are one call and GUI sliders drive the same thing.
7. **Mesher handoff keeps our best meshers and adds an association.** The 2D CDT mesher, the layered extrusion and the isosurface-stuffing
   mesher stay. They consume exact edge/face data from the kernel and return the existing `SimplexMesh` plus a compact face/edge association
   per facet and boundary node. Solvers are unchanged.
8. **A 2.5D layer-stack fast path.** Almost everything photonic here (SOI, rib on sapphire, Bragg cell, photonic-crystal slab) is a stack of
   prisms over 2D partitions. Detecting that and meshing "2D overlay of all layers, then extrude" is the fastest, best-conditioned and
   periodic-friendly route, and the stuffing mesher is kept for the truly 3D rest.
9. **Efficiency from your CSG notes, adapted**: per-node AABBs, lazy tiers, flatten/simplify, content-hash caching and n-ary operations are in;
   a general graph library, per-subtree BVHs and operation-reordering heuristics are not, until a profile asks (section 6).
10. **Do the cheap, high-value part first (Phase 0)**: names, selectors and a `LabelMap` on top of the *existing* meshers, before any kernel
    exists. That already removes the "look up the boundary number in the GUI" problem.

---

## 1. What the current geometry code does, and where it hurts

What is good and should be kept: isosurface stuffing in `mesher3d.rs` (parallel, local, good dihedral angles, fits creases and interfaces),
the constrained-Delaunay 2D mesher working in normalised coordinates, the layered `extrude` that gives identical end meshes (so Floquet pairs
match), `SimplexMesh` as the one contract with the solvers, XML serialisation of any serde value (`xmlio.rs`), and the small safe `expr.rs`.

What hurts, with where it lives:

| # | Problem | Evidence |
|---|---------|----------|
| 1 | Three unrelated geometry languages plus a patch layer. `Geometry::{Rectangle, Block, Shapes, Solids}`, `extrude`, and `Region` bbox relabels. 2D `Shape` (3 kinds) and 3D `Solid` (4 kinds) share only the `Extrusion` bridge. Two separate GUI editors (`tabs.rs` `shapes_editor` ~l.588, `solids_editor` ~l.427). | `state.rs:15`, `mesher.rs:15`, `mesher3d.rs:35` |
| 2 | Boundary numbers are creation order, so inserting a shape renumbers everything after it. Even "a box" has several orderings in the docs: rectangle L,B,T,R; `SimplexMesh::block` 1 x=0, 2 y=0, 3 z=0, 4 z=lz, 5 y=ly, 6 x=lx; `Solid::Box` documented as x-min, x-max, y-min, y-max, z-min, z-max "the numbering of `SimplexMesh::block`", which contradicts the block doc (I read the comments, I did not run it). | `mesher.rs` header, `mesh.rs` header, `mesher3d.rs` header |
| 3 | Boundary labels are assigned *after* meshing by nearest chain segment (2D) or nearest signed distance (3D), with "last shape wins" tie rules. At coincident faces the label depends on tolerance and order; the 2D code has an error path "belongs to no shape edge". | end of `mesh_2d`; `Field::boundary_label` |
| 4 | Conformity is repaired, not guaranteed: `conform_chains` (~140 lines, tolerance `1e-9*scale`) merges near-equal vertices and inserts T-junctions. The README records that slivers from shared edges made the EM problem singular. A debug block with hard-coded coordinates (`3.5e-6, 9e-6`, `NYX_MESH_DEBUG`) is still in the mesher. | `mesher.rs:301` |
| 5 | Curves are polygons. A circle is `max(12, ceil(2*pi*r/h))` chords, controlled by `h` only, not by deviation. `refine()` may split a constrained boundary segment at its chord midpoint, which is inside the circle (inferred from `keep_constraint_edges` + `refine`; not tested). | `mesher.rs` `mesh_2d` |
| 6 | Cells are straight-sided (HOW_IT_WORKS section 4: constant gradients), so any curved surface carries an O(h^2) geometric error whatever the element order. For P2 and Nedelec-2 this is the accuracy ceiling for discs, rings and spheres. | `fem::geometry` |
| 7 | Primitives cannot be rotated or placed freely. `Cylinder.axis` is 0/1/2, boxes are axis-aligned, there is no transform node, no pattern, no revolve. | `Solid` |
| 8 | A hard cap on interface surfaces (`MAX_SURFACES = 64`, `cl: [i8; 64]` per lattice vertex) and a "too many body faces" error. | `mesher3d.rs` |
| 9 | Domain numbers are hard-coded in use cases: SOI 1/2/3; chips 1-4, PML 11-15, 21-22, 41-45. Boundary conditions elsewhere are raw `Bc.labels: Vec<i32>`. | `soi.rs`, `chips.rs`, `problem.rs` |
| 10 | The existing bbox relabel tests *centroids* of cells/facets against a box. It has no notion of "which face" so a box that clips two faces relabels both. It is the partial fix already in place and the right instinct, but at the wrong level. | `mesh.rs` `relabel` |

Seam that already exists and should become the contract: `meshers::adapt::geom::Geometry2::project(p, chain)` is a geometry oracle
("closest admissible point on this chain"), with analytic overrides per boundary label. A B-rep edge/face *is* that oracle.

---

## 2. Requirements (restated from the request)

- CSG-tree variant for building geometry, an internal B-rep for exactness.
- Efficiency, speed, accuracy first. Not bloated.
- No manual domain or boundary numbering, but still addressable from code and GUI.
- Code-first and GUI are two complete routes to the same feature.
- Pipeline: tree -> B-rep -> our meshers -> our solvers.

---

## 3. What to take from each reference (and what not to)

I fetched the AnvilCAD, next.brep.io, DeclaraCAD meshing and Fornjot pages. I did not fetch the Manifold and opencascade-rs repositories; what I say about them (and about Gmsh and OpenCASCADE) is from general knowledge and is worth a quick check before anyone relies on it.

| System | Take | Leave |
|---|---|---|
| **Manifold** (mesh Booleans) | Provenance through Booleans (every output triangle knows its input), always-valid output by construction, n-ary/batch operations, parallelism. | Mesh-only representation: it cannot give exact curved boundaries or analytic projection, which is the whole point here. |
| **opencascade-rs** (OCCT bindings) | Proof that exact B-rep Booleans on quadrics are a solved but very large problem. | A C++ toolchain, long builds and packaging pain on Windows, an LGPL dependency, a meshing module aimed at display. Its naming is by topology index, the classic topological-naming problem. Could be an *optional importer* later, never the core. |
| **DeclaraCAD meshing** | The separation geometry -> hypotheses -> mesh. | Its own docs select edges as `e[0], e[2], e[5], e[17]`: that is the COMSOL numbering pain in code form. Our sizing/BC rules attach to selectors instead. |
| **next.brep.io** | A checklist of what a full kernel contains (sketch solver, fillets, healing, STEP, sheet metal). | Almost all of it. It is a project self-assessment ("no kernel or application tests or benchmarks were run" by its own page), so I treat it as a feature list, not as a dependency candidate. |
| **Fornjot** | Principles: reliability over features, a friendly code-first API, clear errors. | The project is archived (README, 19 Jun 2026: "shut down, goals never reached"). |
| **AnvilCAD** (Rust, own kernels) | A feature list with suppress/rollback, global variables with equations, configurations (= named parameter sets), an agent/CLI surface. | Sketch constraint solving, assemblies, mates, fillets, DFM checks, PDM: scope creep for a simulation tool. Parameters plus expressions replace constrained sketches. |
| **Gmsh/OCC `fragment`, physical groups** (prior art, not fetched) | `fragment` is exactly our arrangement: it makes interfaces conforming and shared. Named physical groups are the right output concept. | Selection by entity tag number and bounding box queries, i.e. the pain being removed. |

---

## 4. Architecture overview

```
 code (Rust builder)  --\                                           /--> display tessellation (viewport, picking)
                         >--> FeatureList (serde, named params) --> eval (lazy, cached DAG)
 GUI (form + viewport) --/                                          \--> RegionComplex (B-rep, provenance, tags)
                                                                              |
                         Selectors (rules) ---------- resolve -------------->|--> LabelMap (names <-> ints), sizing rules,
                                                                              |    periodic pairs, PML, BC/material targets
                                                                              v
                       mesh drivers:  2D CDT | layer-stack (2D overlay + extrude) | 3D stuffing     (femcore, existing code)
                                                                              |
                                           SimplexMesh (unchanged) + Association (facet->face, node->edge/face/vertex)
                                                                              |
                                                          solvers (unchanged), results carry LabelMap
```

Crates and layering (dependency arrow points at what is depended on):

```
app  ->  femcore  ->  geom          (geom depends only on serde; optional rayon behind a feature)
```

`geom` knows nothing about meshes, `SimplexMesh`, physics or the GUI. It exposes: the feature list, evaluation, the region complex, selectors,
a projection oracle, exact queries (area, volume, bounds, classification) and a display tessellation with face ids. `femcore` owns the drivers that
turn a region complex into a `SimplexMesh`, so the meshers keep their code and tests.

---

## 5. The tree: a CSG variant, stored as a flat feature list

### 5.1 Node set (the "variant")

```text
Param            name = value | expression                 (global, named, shared by code and GUI)
Profile (2D)     Rect, Circle, Ellipse, Polygon, RoundedRect, Arc-polyline wire; every segment has a role name
Primitive (3D)   Box, Cylinder, Cone, Sphere, Torus
Sweep            Extrude(profile, range | layers, taper?), Revolve(profile, axis, angles)
Place            Transform(node, translate/rotate/mirror/uniform scale), Pattern(linear | polar, n)
Combine          Union[n-ary], Difference(a, [b..]), Intersect[n-ary]        (shape algebra inside one body)
Simulation ops   Overlay (root), Pad/Pml(box or region, thickness, layers), Cut(plane, keep side) = symmetry,
                 Periodic(sel_a, sel_b, translation), Size(sel, h | max_deviation | layers)
```

The variant relative to textbook CSG: (a) sweeps are leaves, so extruded profiles keep their profile edges as named side faces; (b) the root is
not a single solid but an **ordered Overlay of named regions** (the painter's order you already use, now explicit, with "hole" as a region
kind); (c) simulation ops (PML, periodic, symmetry, sizing) are first-class nodes, not post-hoc patches; (d) the same node set works for a 2D
profile partition and a 3D solid partition.

### 5.2 Storage: flat list, DAG by reference

```text
Feature { id: FeatureId (stable, u32), name: String, op: Op, suppressed: bool }     // Op refers to earlier FeatureIds
FeatureList { params: Vec<Param>, features: Vec<Feature> }
```

A flat list whose references point backwards is already a DAG, in topological order for free. Reasons to prefer it over nested trees or a graph
library: a GUI form edits one row; XML/JSON diffs cleanly in git; undo is a list edit; suppress/rollback (AnvilCAD style) is trivial; a profile
used by two extrudes, or a pattern instance, is one node referenced twice. Verdict on graphs is in section 6.2.

### 5.3 Evaluation

Evaluation is demand-driven in tiers; asking for a cheaper tier never builds a more expensive one:

| Tier | Needs | Gives |
|---|---|---|
| Bounds | params only | AABB per node (GUI bbox, lattice/box size for the mesher, cheap culling) |
| Classify | tree only (no B-rep) | `inside(p)` per primitive and for the whole expression, signed-distance bound near surfaces; this is what stuffing already needs |
| Complex | arrangement | `RegionComplex` with provenance, names, tags, selectors |
| Display | complex | triangles/polylines carrying face ids, for the viewport and picking |

Each node result is cached under a content hash (op parameters after expression evaluation, plus the hashes of its inputs). Changing `w` only
re-evaluates the nodes that depend on it. Meshing dominates cost so this matters most for parameter sweeps and live GUI edits.

---

## 6. Efficiency: your CSG notes mapped onto this design, and the DAG question

### 6.1 Each optimisation, and what I recommend

| Idea | Verdict | How it applies here |
|---|---|---|
| Bounding volumes per subtree | **Yes, AABB per node.** | Computed in the Bounds tier. Disjoint operands make a Boolean trivial; `inside(p)` skips a subtree whose padded AABB does not contain `p`; the arrangement only tests face pairs with overlapping AABBs. |
| BVH over subtrees | **Not for the tree**, maybe for faces. | Models have tens of nodes and tens to a few hundred faces; a sort-and-sweep or uniform grid over face AABBs is enough. Add a BVH only if a profile of a large pattern (hundreds of holes) shows pair-finding dominating. |
| Lazy evaluation | **Yes, as the tier system above.** | The tree stays symbolic until a consumer needs bounds, classification, a complex or a display mesh. |
| Spatial subdivision (octree/grid) | **Already there, keep it.** | `mesher3d` evaluates only lattice cells that can contain boundary ("local and independent"). The new part is feeding it per-cell *local* surface lists from face AABBs, which also removes the `MAX_SURFACES = 64` cap. |
| Flatten and simplify | **Yes.** | Flatten nested unions/intersections into n-ary nodes; remove empties by AABB (`A cap B = empty` when disjoint, `A \ B = A`); `A cap A = A` and duplicate operands removed through hash-consing. Cheap and done at feature-list normalisation. |
| Operation reordering | **Not initially.** | With tens of operands, ordering is noise. Do the safe part (difference operands are culled by AABB against the minuend before anything else) and revisit with measurements. |
| Caching | **Yes, mandatory.** | Content-hash per node, described in 5.3. |
| Batch/n-ary Booleans | **Native, not an optimisation.** | The arrangement takes all primitives at once and evaluates the CSG expression on the inside bits of each fragment, so `union([a,b,c,d])` is one pass by construction. Binary operators are sugar. |
| Robust predicates and coplanar handling | **Yes, and this is the hard part.** | See section 7.3: interned surfaces make coincidence an id comparison; one tolerance scale per model; exact orientation predicates (the `robust` crate is already in `Cargo.lock` through `spade`) for 2D decisions. Your note is right that most failures are floating point and topology failures, not maths. |

### 6.2 Is a graph worth it? Recommendation: a DAG yes, a graph library no

- **What is worth having**: sharing (a profile or sub-assembly referenced more than once), hash-consing (identical sub-expressions are one node,
  so `A cap A = A` and duplicate detection are free), and per-node caching. All three come from the flat list with backward references.
- **What is not worth it**: `petgraph` or a general graph structure. There are no cycles to detect (references point backwards), topological
  order is the list order, and dependency invalidation is "recompute nodes whose input hashes changed". A graph crate adds a dependency and an
  API to learn for no capability the list lacks. A node-graph GUI (Grasshopper style) is a *view* that could be drawn from the same list later; it
  does not change storage.
- **When to revisit**: if sub-assembly instancing of whole region complexes with different parameters becomes common (many instances of one
  parametric cell), promote "FeatureList reference" to a node. The list format allows this without migration.

---

## 7. The B-rep kernel

### 7.1 Representation (struct-of-arrays, u32 ids, no pointers)

```text
Surface   Plane | Cylinder | Cone | Sphere | Torus          (frame + parameters; closest_point, normal, uv <-> xyz)
Curve     Line | Circle | Ellipse | SurfSurf{a, b, range}    (eval, closest_point; SurfSurf solves the 2-surface system with Newton, marched adaptively)
Vertex    point
Edge      curve, [t0, t1], v0, v1
Face      surface, orientation, loops: Vec<Vec<EdgeUse{edge, forward}>>, provenance, name
Region    name, kind (material | hole | pml | ...), shell(s) of face-sides, priority
```

2D models use the same structure one dimension down (profile partition: regions are planar faces, interfaces are edges), so selectors and
names behave identically. `SurfSurf` is the escape hatch for intersections that are not conics (for example skew cylinder-cylinder); it is
exactly how `mesher3d` already treats cutting surfaces (`bisect_surface`, `newton`).

Not in v1, by design: NURBS, offsets, fillets, chamfers, shells, healing, STEP. Rounded *profile* corners are in (they are planar arcs).

### 7.2 The one algorithm: arrangement + CSG labelling

1. Collect the faces of all primitives (after transforms). Intern surfaces (7.3).
2. For each pair of faces with overlapping AABBs, compute intersection curves (closed form for plane/plane, plane/quadric, coaxial quadrics,
   sphere/sphere; `SurfSurf` otherwise). Split edges at intersection vertices; imprint curves on both faces in each face's parameter space
   (planar: local xy; cylinder and cone: (theta, z); sphere and torus: two angles), build loops, find face fragments.
3. For each fragment compute its inside bits for every primitive (exact `sdf` sign at an interior sample point, on both sides of the fragment)
   and evaluate the CSG expression and the Overlay priority on each side. A fragment where the two sides give different regions is a **boundary
   or interface face**; with the same region on both sides it is dropped.
4. Where a fragment lies *on* another primitive's face (shared interned surface), decide sides from normal orientation, not from sampling at an
   offset that would straddle the other primitive.
5. Build shells per region, validate (7.4), attach provenance (8.1).

Every output face descends from exactly one input face, even at Boolean cuts (the cut face of B in A minus B is B's face, flipped). That is why
naming survives, and is the same idea Manifold calls tracking the original ID.

### 7.3 Tolerance and coincidence

- One scale per model, `L` = bounding-box diagonal, one `tol = 1e-9 * L` (the 2D mesher already works this way, and nanometre structures at
  micrometre offsets are the reason). All kernel decisions use it; the mesher keeps its normalisation.
- **Interning**: planes (unit normal sign-canonical, offset), cylinders, spheres and so on are hash-consed to a `SurfId` within `tol`. Two
  coincident faces get the same id, so "are these coplanar" is an integer compare. Edges are interned by (surface pair, range). This removes the
  need for anything like `conform_chains`.
- 2D orientation/in-circle tests use exact predicates; 3D decisions are tolerance-based on analytic parameters.
- **Refuse, don't guess**: grazing contact (distance within a small multiple of `tol` but not zero) is reported as an error naming both
  features ("core.side touches clad.side tangentially: move by more than tol or make them coincident"). Even the next.brep.io self-assessment lists near-tangent
  and nearly degenerate intersections as its remaining sensitive case, so refusing with a clear message is the honest behaviour.

### 7.4 Self-checks (cheap, always on in debug and in tests)

Every edge used by exactly two face-sides (one on the outer boundary); loops closed; Euler relation per shell; region volumes positive; sum of
region volumes equals the volume of the union of all Overlay regions; each interface area appears on both neighbours with equal value. Volumes
and areas are analytic (planar polygons, cylinders, spheres, tori), so these double as correctness oracles against closed forms.

### 7.5 What is exact and what is not

Exact (to rounding in f64): the curve and surface definitions, intersection curves between the supported pairs, areas and volumes, the
projection of any point onto any face or edge, and all topology that is decided by interned ids. Not exact: the arrangement itself in f64, the
`SurfSurf` curves (Newton to tolerance), and anything outside the supported surface set.

### 7.6 Import (three levels, very different costs)

| Level | What | Cost | Recommendation |
|---|---|---|---|
| Mesh import | Existing cell lists (`CellMesh`, `.fmsh`) go straight to the solvers, no kernel involved. | Already there | Keep. Imported meshes are not parametric and have no names beyond their labels. |
| **STEP, analytic subset** | Read `MANIFOLD_SOLID_BREP` / `ADVANCED_FACE` whose surfaces are plane, cylinder, cone, sphere, torus and whose edges are line, circle, ellipse. These map one-to-one onto the kernel's own types, so there is no conversion theory, only a text parser (ISO 10303-21 syntax is simple) and entity mapping. | Medium: parser ~600 lines, mapping ~500, plus a point-in-solid test (~300) because an imported solid has no closed-form signed distance. My estimate, unverified. | Worth doing **after** v1, as a separate optional crate (`geom-step`). Faces with B-spline/NURBS surfaces are refused with a clear message naming the face. Imported faces get roles `face[i]` plus any STEP names, and enter the Overlay as opaque leaf primitives. |
| STEP with NURBS | Free-form surfaces, trimmed NURBS, healing. | Large: it is the part that makes real kernels big. | Do not build. If ever needed, take it as a mesh (tessellate elsewhere, import at the first level) or an optional external converter. |

What to do now so import stays cheap: keep the node set open to an `Import(brep)` leaf, keep topology independent of how it was built (no
"this face came from a Box" assumptions outside provenance), and keep classification a trait (`inside(p)`), not a hard-wired signed-distance formula.

---

## 8. Names and selection (no numbers to look up)

### 8.1 Provenance gives stable names

Every primitive and profile declares **roles**:

| Primitive | Roles |
|---|---|
| Box | `xmin xmax ymin ymax zmin zmax` |
| Cylinder | `side bottom top` |
| Cone / Sphere / Torus | `side bottom top` / `surface` / `surface` |
| Extrude of a profile | `side:<profile edge role>`, `bottom`, `top` (profile edges are named: Rect gives `left bottom top right`, Polygon gives `e0 e1 ...` or user names) |
| Revolve | `swept:<profile edge role>`, `start`, `end` |

A face's **auto-name** is `<feature>.<role>` for outer boundary faces (`core.top`), and `A|B` for interfaces between regions (`Si|BOX`).
Fragments of one original face keep the same name; the group is addressed by name, a single fragment by name plus a geometric predicate.
Insertion, reordering or re-parametrising other features does not change a name, which is exactly what the current creation-order numbers cannot
promise.

Users add **tags** to features or to selector results (`tag("input")`) for physics roles that are not geometric.

### 8.2 Selector language

Selectors are small serialisable expressions; they resolve against an evaluated complex to a set of faces, edges or regions.

```text
atoms     name(core.top)   tag(input)   region(Si)   interface(Si, BOX)   outer   interface(*)
          kind(cylinder)   radius(1.2e-6, tol)   area(>1e-12)
          plane(z = 0)     normal(+z, tol=5deg)     extreme(+z)        in_box(lo, hi, within | touching)
          adjacent(sel)    boundary(region)          from(feature, role)
algebra   a | b     a & b     a - b     !a
```

Examples for the repo's own cases:

| Intent | Selector |
|---|---|
| Fixed support on the anchor | `name(anchor.zmin)` |
| Silicon/oxide interface (moving-boundary term in Brillouin) | `interface(Si, BOX) \| interface(Si, clad)` |
| All outer walls of the simulation window (PEC in the SOI case) | `outer - tag(pml_outer)` |
| Floquet pair of a Bragg cell | `periodic(name(cell.zmin), name(cell.zmax), translate(0,0,period))` |
| Top of the substrate wherever it is | `region(sapphire).boundary & extreme(+z)` |
| Mesh size on the silicon and its interface | `size(region(Si) \| interface(Si,*), h_core)` |

`in_box` has explicit semantics (entirely inside vs touching), fixing the centroid ambiguity noted in problem 10.

### 8.3 Ints still exist, as compiled output

To leave solvers unchanged, `Bc.labels: Vec<i32>`, `Materials` per-domain ints and `Periodic.source/dest` stay. The compile step:

- **domain label = region id**, **boundary label = face id** (one int per B-rep face, minted deterministically by sorting provenance keys then
  centroids). A selector resolves to a set of face ids, which is exactly a `Vec<i32>` of boundary labels. No solver code changes.
- `LabelMap { domains: name <-> int, boundaries: name <-> int + selector hits }` is returned next to the mesh, stored in results files, and shown
  in the GUI Labels panel (name, int, mesh entity count, the rule that produced it). From code: `lm.boundary("Si|BOX") -> Vec<i32>`,
  `lm.name_of_boundary(7)`.
- A user can **pin** an id (`region("Si").id(2)`) for legacy compatibility or external tools. Legacy projects keep working because `Selector::Id(i32)`
  is an atom, and legacy `Geometry` variants load through a `lift()` that builds the equivalent feature list (Shapes and Solids map almost 1:1).
- To verify before relying on it: how `Bc.labels` is matched in assembly. If it is a linear `contains` per facet, hundreds of face ids from a
  photonic-crystal slab would be slow; a bitset or hash lookup is a small change in the solver glue.

### 8.4 GUI picking produces rules, never ids

Click a face: the GUI computes the simplest selector that resolves to exactly the picked set, trying in order auto-name, tag, name plus
predicate, then geometric predicate (plane/normal/box), shows it as text with "resolves to N faces (M mesh facets)", and stores *the rule*. A rule
that later resolves to zero is flagged red instead of silently applying to nothing. Multi-select offers "same normal", "same region interface",
"all faces of this feature" as one-click growths.

---

## 9. Code-first and GUI as two routes to the same feature

- The `FeatureList` is the single source of truth and the saved form (inside `.nyx.xml`, via the existing serde -> XML path).
- **Builder**: thin wrappers that push `Op` variants and return a `FeatureId` handle. No logic lives in the builder, so anything the GUI can
  do, code can do, and the reverse. Loops and conditionals in code just produce more features.
- **GUI**: one list panel (reorder, suppress, rename, edit parameters through expression fields) plus viewport picking. One form per `Op`
  variant, driven by a small descriptor table (about 20 lines per op) instead of the two hand-written editors that exist today.
- **Round trip**: `FeatureList -> Rust code` is a printer (data to text, easy); `code -> FeatureList` is just running the builder. The invariant to
  test: `build(print(list)) == list` and serde round-trip equality.
- **Parameters** are named, global, expression-valued (reusing `expr.rs`). Sweeps are `model.with(&[("w", 450e-9)])`; GUI sliders set the same
  names. This replaces the per-use-case structs (`SoiWire`, `ChipSpec`, `BraggGrating`) as the place where geometry parameters live.
- **Generated blocks**: features created by a code loop are marked `origin: Code` and shown grouped in the GUI, editable through their
  parameters, not individually (otherwise the next code run overwrites hand edits).
- Python users: the saved tree is JSON/XML, so a 100-line pure-Python emitter can write models with no extra dependency in the library. Not part
  of the core.

Illustrative sketch (not compiled; names are proposals). The 300 nm SOI strip from `usecases/soi.rs`:

```rust
let mut m = Model::new_2d();                       // x lateral, y surface normal
m.param("w", 500e-9).param("t", 300e-9).param("h_core", 40e-9);
let clad = m.rect("clad", ("-w/2-1e-6", "-2e-6"), ("w/2+1e-6", "t+1e-6"));
let box_ = m.rect("box",  ("-w/2-1e-6", "-2e-6"), ("w/2+1e-6", "0"));
let core = m.rect("core", ("-w/2", "0"), ("w/2", "t"));
m.overlay([("SiO2", clad), ("BOX", box_), ("Si", core)]);        // painter order = list order
m.size(sel!(region(Si) | interface(Si, *)), "h_core");
let mesh = m.mesh()?;                              // SimplexMesh + LabelMap
let walls = mesh.labels.select(sel!(outer))?;      // Vec<i32>, for Bc{labels: walls, ..}
```

---

## 10. From B-rep to mesh to solver

### 10.1 What the kernel hands to each mesher

| Route | Chosen when | Input from the kernel | Existing code reused |
|---|---|---|---|
| **2D** | 2D model | each edge discretised *once* by deviation and size (shared by both neighbouring faces, so conformity is by construction); face loops; region ids | `mesh_2d` CDT, Ruppert refinement, normalisation |
| **Layer stack (2.5D)** | 3D model that is a stack of prisms over 2D partitions (SOI, ribs, Bragg, photonic-crystal slabs, suspended membranes) | the 2D *overlay of all layers* (common refinement, so every layer boundary is a mesh edge), z ranges, per-layer region ids | `extrude` (3 tets per prism, identical end meshes so periodic pairs match) |
| **3D general** | spheres, tori, rotated or tapered bodies | analytic surfaces (the planes, cylinders and spheres `mesher3d` already knows plus cone and torus), CSG classification, interface faces from the complex, size rules | isosurface stuffing, crease snapping, `refine3d` |
| **Structured PML** | `Pad`/`Pml` feature | the physical mesh's boundary on a box | `meshers::adapt::pml`: it already builds prisms/hexes and labels `(region, material)`; the kernel just supplies the box and thicknesses, replacing the hand-numbered 11-15 / 41-45 in `chips.rs` |

Drivers live in `femcore`; the kernel never produces tetrahedra.

### 10.2 The association (small, additive)

Alongside `SimplexMesh` (unchanged fields, so the ~500 tests are not disturbed) each mesher returns:

```text
facet_face: Vec<u32>          // B-rep face (2D: edge) of every boundary facet
node_on:    Vec<(u8, u32)>    // Interior | Face(id) | Edge(id) | Vertex(id) for every node
```

This enables, with no further geometry knowledge in the solver or adaptation code:
(1) post-refinement projection of boundary nodes onto the *exact* face/edge (the `Geometry2` oracle, now implemented by the kernel, also for
`refine3d` and the 2D adaptation loop);
(2) curved elements (section 11);
(3) selecting results by rule after the solve ("integrate the flux through `interface(Si,BOX)`");
(4) the Labels panel counts.

### 10.3 Other handoff rules

- **Periodic**: `Periodic(sel_a, sel_b, translation)` verifies the faces are congruent under the translation; the mesher meshes face A and copies it to B (layer stack gets this for free; the 3D path meshes one and mirrors the other).
- **Sizing** is a field built from `Size` rules (per selector: max `h`, max deviation, minimum layers across thin bodies, replacing `min_layers`), evaluated per entity, not per global `h` plus `Part.h`.
- **Curvature**: a circle of radius `r` with allowed deviation `d` gets `n = pi / acos(1 - d/r)` points, so geometric error is stated, not an accident of `h`.

### 10.4 Solver side

`LabelMap` resolves names in `ProblemSpec` (add name-keyed `materials.set_named("Si", ..)` and `Bc { on: Selector }` next to the int forms); a single
`resolve(&mut ProblemSpec, &LabelMap)` before `solve::solve` fills the existing int fields. Result files store the map so a saved field is
readable without the project.

---

## 11. Accuracy plan

| Source of error | Today | Proposal |
|---|---|---|
| Boundary vertex position | Chord midpoints can end up off the curve after refinement (2D); 3D warps vertices onto the surface (good). | Project every boundary node onto its exact face/edge using the association; check no cell inverts. |
| Chord count | `h`-driven only. | Deviation-driven (10.3). |
| Straight-sided cells | O(h^2) geometric error on curved boundaries even for P2/Nedelec-2. | Step 1 (cheap, in the mesher): deviation-driven sizing on curved faces so the faceting error is a stated tolerance. Step 2 (only if a benchmark says it is needed): curved P2 cells, which is a **solver-wide refactor, not a local patch** (see 11.1). Benchmark first: disc or ring resonance against the Bessel-function closed form. |

### 11.1 What curved cells would really change (corrected after reading the element code)

`fem::geometry(mesh)` returns one volume and one set of barycentric gradients per cell, and it is called from about 20 sites in about 15 files
(`fem`, `solid*`, `emw*`, `brillouin`, `ring_brillouin`, `fields`). Every element routine relies on those being constant, which is also why the
quadrature rules are exact. A curved (quadratic-map) cell has a position-dependent Jacobian, so each routine would need gradients and `det J` per
quadrature point, and integrands stop being polynomial. Lagrange elements (scalar, solid) are the easy half. The Nedelec edge elements of the
optical solvers additionally need the covariant Piola map, so the EM code changes non-trivially. Only the thin shell of cells touching curved faces
needs it (a small fraction), via a per-cell `Affine | Curved` geometry provider, but the refactor still touches every physics module. Axisymmetric
solvers already treat the revolve direction exactly; this concerns curved cross-sections and spheres. I therefore recommend deferring it behind
step 1 and a measurement. As a rough scale (not measured): an inscribed polygon with chord `h = r/10` loses about `(h/r)^2/6 ~ 1.7e-3` of the
disc area, so a resonance frequency shifts by order 1e-3 relative.
| Interface fidelity | Fitted in 3D, repaired in 2D. | Conforming by construction in both. |
| Booleans | Tolerance-implicit. | Interned coincidence, refusal on tangency, volume/area self-checks. |

I would not claim a convergence order for the curved-element step before it is measured.

---

## 12. Keeping it lean

Budget (lines of Rust, excluding tests), order of magnitude only:

| Module | Lines |
|---|---|
| surfaces and curves (+ projection, intersection pairs) | ~900 |
| topology arena, validity checks | ~500 |
| profiles and sketch wires (named roles, rounded corners) | ~500 |
| primitives, sweeps, transforms, patterns | ~700 |
| arrangement + CSG labelling | ~1800 |
| selectors, provenance, LabelMap | ~700 |
| feature list, params, eval tiers, cache, serde | ~700 |
| display tessellation with face ids | ~400 |
| **total** | **~6200** |

Rules to stay lean: dependencies are `serde` (and `rayon` behind a `parallel` feature, off by default); no NURBS, STEP or fillets until a
concrete model needs them; the kernel does no meshing; one Boolean algorithm; descriptors instead of per-op bespoke GUI code. The new code
replaces `conform_chains`, the per-shape label logic, `Solid::face_at`/`faces_with_extent`, the `Region` relabel path and one of the two GUI
editors, which offsets a good part of the added size.

---

## 13. Migration plan

Sizes are relative (S/M/L), not calendar estimates.

| Phase | Content | Exit criterion |
|---|---|---|
| **0** (S) | Names and selectors on the *existing* meshers: give `Part`/`Body3d` names, emit provenance (facet -> part, edge) from `mesh_2d`/`mesh_3d` (the `Chain{part,label}` data is already there), add `Selector`, `LabelMap`, name-keyed `Bc`/materials, persistent resolution. | Examples 1-19 run with name-based BCs and materials; no number typed in the GUI for them. |
| **1** (L) | `geom` crate: surfaces/curves, profiles, 2D arrangement, feature list + params + eval tiers + serde, display tessellation, 2D mesher driver (replaces `conform_chains`). `lift()` for legacy `Shapes`. `geom_bench` example. | SOI wire, the four chips and the 2D Bragg cross-section reproduce the current mesh quality; n_eff of the SOI TE mode stays within the mesh-convergence band already documented (2.67179); self-checks green. First real timings recorded. |
| **2** (L) | 3D: plane/cylinder/sphere/cone/torus arrangement, layer-stack driver, stuffing driver fed by kernel surfaces with local surface lists (no 64 cap), `Pml`, `Periodic`, `Cut`, `Size`. | Bragg gratings (examples 18-19) through the layer-stack path with identical band-edge frequencies to the documented values within mesh tolerance; one genuinely 3D case (sphere in box) through stuffing. |
| **3** (M) | GUI: feature-list panel, pick -> selector, Labels panel, "Copy as Rust", parameter sliders; retire the two old editors. | Same model built in code and in GUI compares equal as a `FeatureList`. |
| **4** (M) | Association-driven projection after refinement, deviation-driven sizing, 3D adaptive refinement through the oracle. Curved P2 cells only if the disc/ring benchmark shows step 1 is not enough (solver-wide, section 11.1). | Disc/ring closed-form benchmark shows the measured effect; no regression elsewhere. |
| **5** (optional) | STEP import of the analytic subset as a separate `geom-step` crate (section 7.6). NURBS: not planned. | A real user model needs it. |

Guard rails throughout: golden tests that rebuild each existing use case through the new path and compare solver outputs; property tests
(`proptest` as a dev-dependency only) for associativity of union volume, `classify` versus B-rep agreement on random points, and
name stability under unrelated edits.

---

## 14. Risks and open questions

- **Arrangement robustness** is the real risk (coincident and grazing faces). Mitigation: interning, refusal with a named error, restricting v1
  to the listed surface types, and the self-checks. If it slips, the fallback is the existing `Spec3d` painter path behind the same selectors
  (Phase 0 already delivers the naming value).
- **Fragment naming**: when a boolean splits `Si.top` into two faces, the name covers both; users who need one fragment must add a predicate.
  I think that is the right default, but it is a convention to confirm.
- **Curved cells need solver changes** (section 11). Decision needed before Phase 4.
- **Sweeps with taper** (sloped side walls in 3D) are not prismatic, so they leave the layer-stack path. The 2D cross-section already handles the
  common sloped-wall case (`sidewall_deg`); 3D tapered extrusions go through stuffing.
- Not measured: all performance targets (build of a ~50-feature model under about 5 ms, selector resolution of ~500 faces under about 1 ms) are
  guesses to be checked in Phase 1.

Decisions I would like from you:
1. A separate `geom` crate (my recommendation) or a module inside `femcore`?
2. Is STEP/NURBS import ever needed? Default: no.
3. Whether curved P2 cells (a solver-wide refactor, section 11.1) are worth considering at all, after deviation-driven sizing has been measured.
4. Selectors as a text syntax as well as the Rust builder (it makes GUI display and saved files readable)? Default: yes.
5. Start with Phase 0 now, since it delivers the "no looking up numbers" win independently of the kernel?
