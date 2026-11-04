# 31. femcore holds meshes, meshing and adaptation in three module trees

Date: 2026-11-04

## Status

Accepted (phase 2 of the restructure, step 1). Supersedes ADR 7 in part (where the ported mesher lives); the one-native-mesh decision (ADR 5) is carried out in later steps of phase 2.

## Context

Mesh code was spread over fifteen top-level modules of `femcore` (`mesh`, `view`, `cellmesh`, `geomdrive`, `mesher`, `mesher3d`, `refine3d`, `meshing`, `meshing_pml`, `pml3`, `physics_size`, `estimate`, `estimate3`) and a ported `meshers/` tree with its own mesh types, quality, fixtures, I/O, adaptation, validation and PML generator. Two feature streams (the mesher features and phase 2 itself) were going to edit the same files, and jj does not follow renames, so a move made later would turn every concurrent edit into a modify/delete conflict.

## Decision

* The files move, without any logic change, into three module trees (the first commit of the phase, `restructure-2-mesh` a4f8aa69):
  * `femcore::mesh`: `mod.rs` (`SimplexMesh`, was `mesh.rs`), `view`, `cellmesh` (the nyx `CellMesh`, legacy), `native` (the dimension-generic mesh, was `meshers/mesh.rs`), `cells` (was `meshers/cellmesh.rs`), `interop` (legacy), `io` (`.fmsh`), `quality`, `fixtures`.
  * `femcore::meshing`: `mod.rs` (`MeshPlan`, `mesh_model`, was `meshing.rs`), `geomdrive`, `mesher`, `mesher3d`, `refine3d`, `pml` (was `meshing_pml.rs`), `pml3`, `pml_box` (the structured PML generator, was `meshers/adapt/pml.rs`), `physics_size`, `validate` (the independent Helmholtz checker).
  * `femcore::adapt`: the metric adaptation (`metric`, `field`, `topo`, `ops`, `drive`, 2D and 3D), `estimate`, `estimate3`, `par`, `timer`.
* The old module paths do not survive: callers were rewritten in the same phase (second commit) and no re-export remains.
* `mesh::native` and `mesh::cells` stay beside the legacy types until the later steps of the phase delete `mesh::cellmesh` and `mesh::interop` and make the native mesh the only one (ADR 5). The move changes no behaviour; the tests, goldens and contract exemplars are unchanged.

## Consequences

* A mesher feature stream (the `FEATURE_DEPENDENCIES.md` items) can start from the move commit; its files have their final paths.
* `docs/ports/meshers/*` still describe the old `meshers::` paths; the provenance and licence notes remain valid.
* `meshers::par` is still its own module (`adapt::par`) until it is replaced by `femcore::par` (ADR 11).

## Code

`crates/femcore/src/{mesh,meshing,adapt}/`, `crates/femcore/src/lib.rs`.

## Contract

Contract `mesh` (`crates/femcore/tests/contracts.rs`, `mesh_*`) is unchanged by this step.
