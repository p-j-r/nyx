# 2. Geometry is a serialisable feature list addressed by names and selectors

Date: 2026-10-31

## Status

Accepted. Decisions 1 to 5 of `docs/GEOM_HANDOFF.md` section 0 are Paul's ("do not reopen").

## Context

The earlier code had three unrelated geometry languages (`Rectangle`, `Block`, `Shapes`, `Solids` plus `extrude` and bbox relabels), boundary numbers that depended on creation order, labels assigned after meshing by nearest chain segment, and domain numbers hard-coded in the use cases (SOI 1/2/3, chips 1-4, PML 11-15, 21-22, 41-45). See `docs/GEOMETRY_ARCHITECTURE.md` section 1 for the ten specific problems.

## Decision

* The model of a geometry is one flat, serialisable **feature list** (`geom::model::FeatureList`) with named parameters and expressions. The Rust builder, the GUI form and "Copy as Rust" all edit or print this same data: code-first and GUI are equivalent.
* A model evaluates to a **region complex**: named regions with shared interface faces, produced by one n-ary arrangement plus a CSG expression, so union, difference and intersection are one algorithm.
* **Names, not numbers.** Every region, edge and face descends from one primitive role. Boundary conditions, materials, mesh sizes, periodic pairs and PML attach through **selectors** (a small serialisable rule language). Integers exist only as a compiled output, `LabelMap`, which is what the solvers read.
* Parameters are shared by name, so sweeps are one call and sliders drive the same data.

## Consequences

* All geometry in the app is a geom model (`ffa1bf6` removed rectangle, block, shapes, solids). Use cases build geom models (`usecases/soi.rs`, `bragg.rs`, `chips.rs`).
* Boundary conditions can stay valid when a shape is inserted; the GUI can show "click a boundary, get the rules that pick it" (`9b49399`).
* `Bc.labels: Vec<i32>` remains in `ProblemSpec`; it is the compiled form and is stable only for one evaluation of one model (a limitation of the interface, see contract `problem_spec`).
* Persistent naming makes the saved model, not the mesh, the artefact that is versioned and hashed (`nyx-store`, `nyx-session` handles).

## Code

`crates/geom/src/{model,eval,labels,select,sim,print,diff}.rs`, `crates/femcore/src/geomdrive.rs`, `crates/app/src/{geomui,state}.rs`. History: `08dd007`, `ffa1bf6`, `4a3cb55`, `5ebb4dd`, `bfa09c1`, `9b49399`.

## Contract

`contracts/geom_model`: exemplar feature lists in `crates/geom/golden/*.json` already evaluate to closed-form expectations (`tests/goldens.rs`); `crates/geom/tests/contracts.rs` adds the byte-exact JSON round trip of the exemplars and the dependency rule of ADR 3.
