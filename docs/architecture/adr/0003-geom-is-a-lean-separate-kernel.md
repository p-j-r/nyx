# 3. geom is a lean, separate, analytic kernel that depends on nothing of nyx

Date: 2026-10-31

## Status

Accepted. Paul's decisions 1, 2, 3 and 6 in `docs/GEOM_HANDOFF.md` section 0; the STEP/Phase 5 plan in `GEOMETRY_ARCHITECTURE.md` was dropped on 2026-10-31.

## Context

A full CAD kernel was out of scope. The geometries nyx needs (strips, ribs, sloped walls, discs, rings, holes, photonic-crystal slabs, spheres) are analytic, and almost every photonic structure is a stack of prisms over 2D partitions.

## Decision

* `geom` is a separate crate, roughly 6k lines in the plan, with `serde` as the only required dependency (`rayon` optional, dev-dependencies `serde_json`, `proptest`). **No dependency on `femcore`, `faer`, the GUI or physics, ever.** The dependency direction is `app -> femcore -> geom`.
* Analytic geometry only: planes, cylinders, cones, spheres, tori; lines, circles, ellipses and procedural surface-surface curves. No NURBS, fillets, shells, sketch solver or assemblies.
* **Import is by mesh only.** No STEP, IGES or `Import` leaf.
* A 2.5D **layer-stack fast path** handles prismatic 3D models (mesh the 2D overlay of all layers, then extrude), which is also the periodic-friendly route.
* geom resolves *which* entities a rule picks; femcore decides what a condition means. Meshing is not in geom (display tessellation only).
* The curved-element refactor of the solvers is not geom's job. geom only makes it possible (exact queries, closest-point oracle).

## Consequences

* geom duplicates the small expression grammar of `femcore::expr` on purpose (rough spot R06).
* There is no 3D geom meshing pipeline: 3D adaptation starts from a tet mesh (extrusion or import) (`MISSING_FEATURES.md` section 3).
* The meshers use geom through `femcore::geomdrive` and `meshing::GeomOracle2`, never the reverse.
* The source is about 10.4k lines (`crates/geom/src`), well above the 6k of the plan; in practice the leanness that held is the dependency list, not the line count.

## Code

`crates/geom` (all of it), `crates/femcore/src/geomdrive.rs`, `meshing.rs`. History: `08dd007`.

## Contract

`contracts/geom_model` (see ADR 2). The dependency rule is checked by `crates/geom/tests/contracts.rs`, test `geom_depends_on_serde_only`.
