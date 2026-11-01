# Changelog

## Revolve

* `Revolve` (builders `revolve`, `revolve_along`): a full turn of a `Rect` or `Polygon` profile about the sweep axis through `center`; the
  profile's `(x, y)` are radius (>= 0) and height along the axis. Segments off the axis become faces `swept:<profile role>`: cone sides,
  cylinder walls and flat rings (new face types `CylSide`, `Annulus`, and a `flip` flag on cone sides so the normal always points out of the
  solid), with a circle edge at every corner off the axis. Exact volume by Pappus, signed distance, projection and parametrisation. Same
  clear-of-everything rule as the other curved bodies. A profile that crosses the axis or itself is refused. Golden: stepped bell; tests: tube,
  bell, revolve about x, refusals, `to_rust` round trip.
* The remaining work is in [`../TODO.md`](../TODO.md).

## Any sweep axis, cone bodies, Pml in the polyhedral path

* **Cylinders, cones and extrusions along x or y.** `Cylinder.axis` is honoured, `Extrude` gained `axis` (serde default z; the profile's
  `(x, y)` are the two other axes in cyclic order), builders `cylinder_along`, `extrude_along`, `cone_along`. All prismatic solids of a model must
  share one sweep axis (mixed axes are refused naming both). Internally the model is built in the frame `(axis+1, axis+2, axis)` so the layer stack
  is unchanged; every value the API returns is converted back (`Complex::wv/iv`): bounds, positions, normals, boxes, projections, fractions,
  tessellation, stuffing surfaces, `Pml` sides and window. Box roles keep world names (`xmin ... zmax`).
* **`Cone`** (frustum, either radius may be 0): a body with faces `side`, `bottom`, `top` and circle edges where they meet (seam vertex at angle 0),
  exact volume, area, distance, projection (in the half plane `(rho, z)`), parametrisation. Clear-of-everything rule as for spheres and tori.
  Goldens/tests: diaphragm in a box, pointed horn along x, touching floor refused.
* **`Pml` in the polyhedral path** (rotated / tapered / slanted models): pad pieces are polyhedral boxes, regions split as in 2D and the prismatic 3D path;
  pieces with no face are dropped.
* Fixed on the way: the stuffing plane normals of lateral faces were not converted to world coordinates.

### Still not implemented

`Revolve` of general profiles, a sphere / torus / cone / cylinder cut by a plane or touching another solid (the faces would have to be split
along intersection curves), skew or tilted cylinders, curved solids together with rotated / tapered / slanted polyhedra, bodies that touch the
face of their host (a cone standing on a floor).

## Polyhedral path (rotations, tapers, slanted cuts)

* **General planar-faced solids** (`planar.rs`): boxes in any orientation, `Rect`/`Polygon` extrusions (concave profiles are triangulated into
  convex prisms whose internal faces vanish in the arrangement), tapered extrusions of convex profiles (`Extrude.taper_deg`, new field with a
  serde default), rotations about any axis (`Transform`, `Pattern` polar), mirrors, scale, patterns, and `Cut` with any normal. Selected
  automatically when a model uses a rotation about x or y, a taper or a slanted cut.
* **Arrangement**: plane by plane. Every plane that carries a face is arranged in 2D (reusing the 2D arrangement) from the cross-section of each
  solid, with the composite class (region above, region below); cells whose two sides differ are the faces. Edges are the cell boundaries merged
  by 3D line (so T-junctions between planes are resolved) and vertices their end points; region volumes follow from the divergence theorem.
  Faces on one plane with the same pair of regions are one face (possibly several components). The closure, volume and edge-incidence checks
  hold on 2000 random rotated / cut models (`tests/props.rs`).
* Goldens added: box with a slanted cut (5 faces, 9 edges, 6 vertices), tapered frustum. `tests/planar.rs`: volumes against closed forms, Euler
  numbers, rotated cubes inside and across other solids, concave beam, pattern of tilted bars, refusals.
* Not in this path: curved solids (cylinder, sphere, torus) together with rotated / tapered / slanted geometry, `Pml` (refused with a message).

### Still missing (needs curved-surface intersection curves)

Cone and frustum solids, `Revolve`, a sphere / torus / cylinder cut by a plane or by another curved solid, skew or tilted cylinders, cylinders
along x or y, curved solids mixed with rotated polyhedra. These need the face-splitting of curved surfaces along `SurfSurf` curves in the surface
parametrisation; they are refused with a message that names the features, not approximated.

## Edges, vertices, curved bodies, refusal catalogue

* **3D edges and vertices** (`topo3.rs`): a fragment-level B-rep derived from the layer stack. Vertical edges run along overlay vertices
  through runs of layers where at least two distinct faces meet; horizontal edges are overlay edges at a level where at least two faces
  meet (lateral faces below and above, caps on both sides); vertices are overlay vertices at levels. Edges are shared by two faces, three or
  more at a junction of regions. Names: sorted face names joined by `/`. A new check verifies that every edge is shared by at least two
  faces and that the edges around every face add up to its boundary (box: V - E + F = 2). Selectors (`name(a/b)`), oracle, fractions, signature use them.
* **Curved bodies**: `Sphere` and `Torus` (axis z), each one closed face `surface` without edges, exact volume, area, distance, projection and
  parametrisation. They sit inside a prismatic host and must be clear of every other face (exact distance for spheres, bounding box for tori)
  and of each other (disjoint, or one sphere strictly inside another). Anything that touches or crosses is refused with both names.
  Goldens: sphere in a box, concentric spheres, torus in a box.
* `tests/refusals.rs`: one test per refusal class (undecidable 2D geometry, invalid input, unsupported 3D content, touching curved bodies,
  labels / selectors / periodic pairs); every message must name the features involved.
* Still not implemented: intersection curves between curved surfaces and planes or each other (`SurfSurf`), cones, revolves, tapered sweeps,
  rotations about x or y, cylinders along x or y, slanted cuts. These are refused, not approximated.

## Naming schemes

* `NameScheme::{Feature, Region, Numbered}` on the `Overlay` (serde default `Feature`, so existing files are unchanged), `Model::set_naming`, printed by `to_rust`. Works for 2D edges and 3D faces.

## Steps C and D (prismatic 3D, layer stacks) - delivered in part, see "Not delivered"

* 3D feature lists: `Box`, `Cylinder` (axis z), `Extrude` of any 2D profile, Transform (translate, rotate about z, mirror, uniform
  scale), Pattern, Booleans, `Overlay`, `Pml` (26 side combinations in 3D), `Cut` (axis-aligned planes; faces `<feature>.cut`, tag
  `symmetry`), `Periodic` (verified congruent by sampling both sides), `Size`, `Tag`. Roles: Box `xmin xmax ymin ymax zmin zmax`, Cylinder
  `side bottom top`, Extrude `side:<role> bottom top`.
* Evaluation as a **layer stack** (`stack.rs`): z levels, one 2D arrangement per layer, and the *overlay*, the common refinement of all layers
  (the 2D arrangement run with composite per-layer classes). `Evaluated::layer_stack()` returns the layers (z range, partition, region of every
  overlay cell) and the overlay. Detected for every 3D model this build can represent; non-prismatic content is refused (below).
* 3D faces derived from the stack: lateral faces = overlay edges swept over runs of layers (merged while the sides and name agree), caps =
  overlay cells at levels where the region below differs from the region above. Names follow the 2D rules (`A|B`, `<feature>.<role>`).
  Volumes, face areas, projection (`oracle().project(.., Face(f))`), fractions (two per face), `classify`, selectors and `LabelMap` all work on faces.
* Self-check for 3D: every face separates two sides, volumes are positive, and every region is closed (the outward area vectors sum to
  zero in x, y and z, to a few tolerances).
* `surfaces_for_stuffing()` (planes and cylinders, coincident faces share one surface, with the extent of its faces), 3D display tessellation,
  `diff` kind `ZOnly` (same topology, same xy overlay geometry, different z).
* Golden models: three Bragg periods (symmetric, asymmetric, uniform) with a Floquet pair, cylinder cavity, concentric cylinders.
* Property test: random 3D box / cylinder models never panic or fail the closure check.

### Not delivered (needs the general 3D arrangement)

* Sphere, Cone, Torus, Revolve, tapered sweeps, rotations about x or y, cylinders along x or y, slanted cut planes, skew cylinders and the
  `SurfSurf` intersection curves. A model that needs them is refused with a message (`Cut ... slanted`, `Cylinder ... axis`, `must be the input
  of an Extrude`). The goldens "sphere in a box", "torus", "concentric spheres", "box with a slanted cut" and "two skew cylinders" are therefore
  not provided; the slanted cut and skew cases exist only as refusals.
* 3D edges and vertices of the B-rep (the 2D edges and vertices of the overlay are available through `layer_stack().overlay`), 3D curve
  queries and 3D `discretize_edge`; `adjacent(...)` between 3D faces; the Euler relation (replaced by the closure check).
* Because everything prismatic is also a layer stack, `layer_stack()` never returns `None` for a model that evaluates.

## Step B (tiers, cache, queries, sweep support) - delivered

* Normalisation of region expressions (flatten, drop empties, `A cap A = A`, bounding-box-disjoint intersections and difference
  operands) and pruning of primitives nothing refers to; the primitive and pair-intersection caches are keyed by content hashes
  (`Evaluator::stats()` counts builds and hits; a one-parameter change rebuilds one primitive and the pairs that involve it).
* `FeatureList::to_rust` (builder code; round trip checked as pure data by a small interpreter in `tests/common/interp.rs`).
* Queries: `length`, `volume`, `interface_area`, curve evaluation / tangent / curvature, `min_radius_of_curvature`, `chord_deviation`,
  `facet_deviation`, `point_at_mid`, `project_midpoint` (seam safe), `discretize_edge` (deviation and size, shared by both neighbours).
* Display tessellation with region / edge ids (`Evaluated::tessellate`).
* Addendum 1: `topology_signature`, stable ids within a signature (edge order is by first origin and the fraction along the origin segment, not by
  position), `fractions_of` / `position_at` (round trip to 1e-12 L), `diff -> ChangeReport`, documented continuous parametrisation (see `query.rs`).
* Property tests: random feature lists never panic or fail a self-check (6000 cases run), union area associativity / commutativity,
  `A cap A = A`, unrelated features keep existing names.
* Bug found by the property tests and fixed: the collinear re-merge removed the wrong edge when the second fragment had the lower index.
* New refusal: a primitive thinner than `64 tol` (about `6e-8 L`) is refused by name; before it produced an inconsistent loop.
* Conflict with the brief: the tolerance is relative to the model size, so adding a distant feature changes `tol`; input within a few
  `1e-9 L` of degenerate can therefore change meaning (the name-stability property is tested on well-separated input).
* Measured (release, 4-core VM): 49-feature model 0.77 ms cold, 4-term selector over 250 edges 6 us, 200-hole plate 0.94 ms,
  100-point width sweep 70 ms. The pair cache saves little at this size because the arrangement dominates; see the README.

## Step A (2D kernel) - delivered

* Added: curves (line, circle, ellipse; canonical parametrisation, exact Green-theorem areas, closest point, curve-curve intersection
  with tangency and near-coincidence reporting), 2D profiles (Rect, RoundedRect, Circle, Ellipse, Polygon) with fixed role names,
  Transform / Pattern / Union / Difference / Intersect, `Overlay` (painter order, Hole/Void/Pml kinds, pinned ids), `Pml` (window +
  sides -> physical region and one `<region>.pml_<sides>` piece per side combination, empty pieces dropped), `Cut` (`<feature>.cut`,
  tag `symmetry`), `Periodic` (verified congruent), `Size`, `Tag`.
* n-ary planar arrangement with interned curves and vertices, owner-oriented side decisions, tolerance `1e-9 L`, refusal band
  `64 tol` (tangent contact, near-coincident curves, corners that nearly touch, self-crossing profiles, thin slivers) with errors that name the features.
* Region complex: regions with loops (holes assigned to outer loops), edges (left to right normal convention), vertices, auto-names
  (`feature.role`, `A|B`, `a/b`), provenance, always-on self checks (`Complex::check`).
* Selectors: text form with deterministic printing and serde form, `resolve` with typed sets and warnings, `suggest`.
* `LabelMap`: deterministic domain and boundary ints, pins, name <-> int in both directions, groups with their rule.
* Oracle (`project` on edges / vertices / regions), serde round trip (JSON bit-exact; XML value encoding checked), builder.
* Golden models in `golden/`: SOI strip, SOI rib with unequal sloped walls, the four ring chips with PML, Bragg cross-section,
  plate with hole, disc, annulus, concentric discs; closed-form areas to 1e-10.
* Conflict with the brief / design choices (see README "Decisions"): the brief's `Pml{target,...}` is implemented as a window plus
  region splitting; edge ids are ordered by (first origin, fraction along the origin segment) rather than by centroid so they do not
  depend on parameter values (addendum item 2).
