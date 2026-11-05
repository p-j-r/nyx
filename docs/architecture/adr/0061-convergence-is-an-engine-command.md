# 61. Mesh convergence is an engine command

Date: 2026-11-05

## Status

Accepted. Replaces the `mesh_converge` check of the old session tools (ADR 15), which existed only for the silicon wire.

## Context

A number from a finite-element solve means little without a second mesh. The old `check_run` repeated a solve at 0.6 times the element sizes, but only for the SOI wire use case. After ADR 57 no tool can do it for a project of any physics.

## Decision

`convergence_check {physics, factor = 0.6, tolerance?}` solves the physics of the current project on a scratch engine, then on a project whose mesh is finer by `factor` (element size, adaptive size factor and chord deviation times `factor`, points per wavelength divided by it, the layers of a sweep with `1 / factor` times as many cells), and compares every eigenvalue (real and imaginary part), every value that orders the solutions and every numeric scalar. The reply has a row per quantity (coarse, refined, difference, relative difference, pass or fail when a tolerance was given), both cell counts and unknown counts, the worst relative difference and, with a tolerance, `converged`. The project, its mesh and its results are not changed. It is an expensive command (two solves).

## Consequences

* A host or the window can say whether a quoted number is converged for any physics; the guide cards tell a host to run it before quoting.
* The difference of two meshes is an indication, not an error bound; a smaller `factor` shows whether it is still falling.

## Code

`crates/nyx-engine/src/convergence.rs`, test `a_finer_mesh_confirms_or_refutes_a_number_through_convergence_check` in `crates/nyx-engine/tests/examples.rs`.
