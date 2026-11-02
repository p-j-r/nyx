# geom

A lean CSG-over-B-rep geometry kernel for FEM: a serialisable feature list (parameters, profiles, Booleans, an ordered `Overlay` of
named regions) is evaluated into a **region complex** with shared interface edges, provenance-based names, rule-based selectors and a
compiled `LabelMap`. Code and GUI edit the same data. No meshing happens here.

Status: **steps A and B (2D) are delivered; step C/D are delivered for prismatic, polyhedral and (clear-standing) curved 3D models; what is still refused is listed in `tests/refusals.rs` and in [`../docs/GEOMETRY_ARCHITECTURE.md`](../docs/GEOMETRY_ARCHITECTURE.md)**; formerly: prismatic 3D models (boxes, z-cylinders, extruded profiles, layer stacks); spheres, tori, cones, revolves and slanted or skew geometry are not (see `CHANGELOG.md`); see [`CHANGELOG.md`](CHANGELOG.md) and the brief in [`../docs/GEOM_HANDOFF.md`](../docs/GEOM_HANDOFF.md).

```bash
cargo test --release -p geom
cargo clippy -p geom --all-targets -- -D warnings
GEOM_WRITE_GOLDEN=1 cargo test --release -p geom --test goldens   # regenerate golden/*.json from tests/common
```

## The wafer frame

Everything is in the frame of the wafer, with no conversions anywhere: **x lateral, y the direction of propagation, z the surface normal**
(right-handed). A 2D model is the **x-z cross-section**: its points are `(x, z)` and its world coordinates `(x, 0, z)`. A 3D extrusion
defaults to the propagation axis, y, with the profile `(x, z)`.

* `Extrude { profile, axis (default Y), from, to, taper_deg }` (older files with `z0`/`z1` still load); builders `extrude(name, profile, y0, y1)`,
  `extrude_tapered`, `extrude_along(name, axis, ..)`. Caps are named by the world axis: `ymin`/`ymax` for the default.
* Rotations are right-handed about the named world axis (`Rotate`, `Polar`), 2D ones about y turn the x-z plane. A 2D `Translate`/`Mirror`
  out of the x-z plane and a 2D `Cut` normal with a y component are errors.
* `Pml` sides are `xmin xmax zmin zmax` in 2D (`ymin ymax` too in 3D); pieces are `<region>.pml_zmin`, `.pml_xmax_zmin`, ...
* Selectors: `plane(z = 0)`, `extreme(+z)`, `normal(+z)`, `in_box([x0, 0, z0], [x1, 0, z1])` (2D ignores y).
* Internal design (the only one): a prismatic model is built in a *sweep frame* (profile plane first, sweep axis last), a permutation of the world axes
  fixed by `scene::axis_perm`; `Complex` keeps its data in that frame and its fields are `pub(crate)`. Every public query (`bbox`, vertices, `classify`,
  `project`, `facet_normals`, `Oracle`, the `Evaluated` accessors) answers in world coordinates, and the test `a_y_sweep_and_a_z_sweep_answer_alike`
  pins that the choice of sweep axis does not change any answer.

```rust
use geom::{eval::Tier, model::*, select::Selector};
let mut m = Model::new_2d();
m.param("w", 500e-9).param("t", 300e-9);
// rect(name, (x_min, z_min), (x_max, z_max)): the x-z cross-section
let clad = m.rect("clad", ("-w/2-1e-6", "-2e-6"), ("w/2+1e-6", "t+1e-6"));
let bx   = m.rect("box",  ("-w/2-1e-6", "-2e-6"), ("w/2+1e-6", 0.0));
let core = m.rect("core", ("-w/2", 0.0), ("w/2", "t"));
m.overlay(&[("clad", clad), ("BOX", bx), ("Si", core)]);       // painter order: later wins
let ev = m.eval(Tier::Complex)?;
let lm = ev.label_map()?;                                       // names <-> ints
let si_box = lm.boundary("BOX|Si");                             // Vec<i32>
let walls = ev.select_labels(&Selector::parse("outer")?)?;
```

## 3D: three paths

Prismatic models (boxes, cylinders and extrusions along one axis) use the layer stack below. Models with a rotation about x or y, a tapered extrusion or a
slanted `Cut` use the polyhedral path (`planar.rs`): convex planar solids combined by Booleans, arranged plane by plane. Spheres and tori are
edge-less bodies added to the prismatic stack. Curved surfaces that would have to be split along intersection curves are refused.

## 3D: layer stacks, edges and curved bodies

A 3D model is evaluated as a stack of layers (cut at every z level where a solid starts or ends); each layer is a 2D arrangement, and the
*overlay* is the common refinement of all of them, so every layer boundary is an edge of it. `Evaluated::layer_stack()` gives the layers
(range along the sweep axis, 2D partition, the region of every overlay cell) and the overlay: exactly what a "2D overlay mesh + extrude" driver needs, with the
names of side faces from the overlay edges and of caps from the levels. Faces are fragments (a face split by another layer's features keeps one name).
Boundary faces are `<feature>.<role>` with roles `xmin ... zmax` (Box), `side bottom top` (Cylinder), `side:<role> bottom top` (Extrude; `bottom`/`top` are the caps of the sweep axis, `ymin`/`ymax` along y);
interfaces are `A|B`; the faces of a `Cut` are `<feature>.cut`. The normal of a face points from `first` to `second` region
(lateral: left to right of the swept overlay edge; caps: below to above). 3D edges and vertices are derived from the face fragments
(`edges3()`, `vertices()`): an edge is shared by two or more faces and named by them (`box.xmin/box.ymin`). `Sphere` and `Torus` bodies
are added as edge-less closed faces when they are clear of everything else (see `CHANGELOG.md` for the exact rule).

## Decisions and conventions

* **Tolerance**: `tol = 1e-9 * L`, `L` the diagonal of the model; undecidable band `64 tol` is refused, not guessed (see `arrange.rs`).
* **Edges** run along the canonical direction of their curve (lines: angle in (-90 deg, 90 deg]; circles and ellipses counter-clockwise from
  the major axis, i.e. from angle 0 for circles). `left` / `right` are the regions on either side; the edge normal points from `left` to `right`.
  A `None` side is the void (outside, or a removed hole).
* **Naming scheme** (`Model::set_naming`, stored in the `Overlay`): `Feature` (default) names a boundary after the shape that drew it,
  `box.bottom`; `Region` after the material it bounds, `BOX.bottom`; `Numbered` after the region and a running number, `BOX.3`, which needs no
  roles at all. Interfaces are always `A|B`. Whatever the scheme, `from(feature, role)` still finds a face by the shape that produced it.
* **Names** (default scheme): boundary edges (one side void) `<feature instance>.<role>` taken from the origin that belongs to the inside region;
  interfaces `A|B` with the region names sorted. Roles are fixed per primitive (see `prim.rs`). Pattern instances append `[i]`.
* **Ids** are dense and deterministic; edges are ordered by (first origin feature id, instance, role, fraction along the origin segment),
  which does not depend on parameter values.
* `A & extreme(+y)` evaluates `extreme` inside the set `A`.
* `Pml { window, sides, thickness, layers, grading }` splits every material region into the physical part and one
  `<region>.pml_<sides>` region per side combination (kind `Pml`); `Evaluated::pml_regions()` carries the data. Pieces that end up empty
  (painted over) are dropped.

## Limits worth knowing

* Prismatic 3D only (see above); everything else is refused with a message.
* Tolerances are relative to the model size, so features closer than `64 * 1e-9 L` (about `6e-8 L`) to something else are refused; in a layer stack
  the z levels follow the same rule.

* A line whose direction passes through the vertical flips its canonical direction (so `left` and `right` swap); a model swept through
  that angle changes its topology signature.
* Ellipse intersections are found numerically (720 samples per curve); two roots closer than a sample interval can be missed.
* `robust` predicates are not used: all decisions are tolerance based on interned geometry.

## Performance (release, one core; `cargo run --release --example geom_bench`)

| Case | Time |
|---|---|
| evaluate a 49-feature model (48 regions, 250 edges), cold cache | 0.77 ms |
| Classify tier only | 0.15 ms |
| resolve a 4-term selector over 250 edges | 6 us |
| label map | 0.1 ms |
| 100 evaluations of a width sweep (shared cache) | 70 ms (primitives: 99 built, 4701 reused; pairs: 893 computed, 19907 reused) |
| 200-hole plate, cold | 0.94 ms |
| one Bragg period (3D, two layers), cold | 0.12 ms |
| 100-hole photonic-crystal slab (3D: 106 faces, 104 overlay edges), cold | 1.1 ms |

The pair cache removes the intersection work of unchanged pairs, but at these sizes the rest of the arrangement (classification, loops) is
the same cost as the intersections, so a warm cache is only slightly faster than a cold one. Caching the arrangement per connected group
is the next step if a profile of a large model asks for it.

## Parameter sweeps

`Evaluated::topology_signature()` hashes the combinatorial structure (names and provenance included, parameter values excluded). Equal
signatures guarantee identical region / edge / vertex id sets in identical order. `fractions_of` / `position_at` give parametric
positions of boundary nodes, `diff(old, new)` classifies a change as `None`, `Morph` (same signature) or `Remesh`; `ZOnly` needs layer
stacks (step C/D). Edge direction is canonical per curve kind, so a line rotated through the vertical swaps `left`/`right` (a signature change).
