# 67. `pmc`: a magnetic wall for `ems` cavities and Floquet cells

Date: 2026-11-14

## Status

Accepted (benchmark A3, bookmark `restructure-10-a3-snowflake`). Extends the cavity kernel of ADR 45.

## Context

The `ems` cavity study puts a perfect electric conductor on every exterior boundary that is not a Floquet face (ADR 45: "Outer walls not otherwise described are perfect electric conductors"). A slab with a mirror plane at its mid-plane has two parity classes. The TE-like modes of a photonic crystal slab, the ones the literature calls the even modes, have the tangential electric field even under the reflection, so the mirror plane is a perfect magnetic conductor for them. With only PEC walls the even modes cannot be isolated, and the cell has twice the volume it needs.

A magnetic wall is the natural condition of the curl-curl form: nothing is constrained on it.

## Decision

* `ems` has a boundary condition `pmc` ("Magnetic wall (mirror plane)", no parameters, cavity study only). The facets with its boundary numbers are left out of the PEC set of the cavity kernel; the edge dofs and the nodal gradient space on them are free, as inside the domain.
* `Problem` gains `pmc: Vec<i32>`; `decode` fills it and refuses the condition outside the cavity study.
* `physics_ems` descriptors 1.2 (additive).

## Consequences

* A mirror plane halves the cell and selects the even modes: checked against the closed form of a box whose upper wall is magnetic (`a_magnetic_wall_is_the_mirror_plane_of_the_doubled_box`, the modes of the doubled PEC box with `sin(p pi z / 2d)`, p odd).
* The odd modes need the mirror plane electric: that is the default, so they are the same project with `pmc` removed on the half cell.
* `sms` has no such condition: a mirror plane for the solid is a fixed normal displacement (`fixed` with the component that is normal).

## Code

`crates/nyx-ems/src/problem.rs`, `crates/nyx-ems/src/cavity.rs` (`assemble_emw_eigen`), `crates/nyx-ems/src/lib.rs` (`bc_kinds`), test in `crates/nyx-ems/tests/emw3d.rs`.
