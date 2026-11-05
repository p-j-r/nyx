# 60. The examples are project files; the use cases are gone

Date: 2026-11-05

## Status

Accepted. Implements "the library contains no built-in devices or use cases; devices are `.nyx` files outside the program" (plan constraint 3) and the deletion of `usecases/` (plan phase 7).

## Context

The window opened with 18 projects built by code (`app/examples.rs`), two of which (the silicon-on-insulator wire and the Bragg grating) came from `femcore::usecases`, which the MCP session also called. Their geometry, materials and conditions lived in Rust, so changing an example meant recompiling, and the tests of other crates used the use cases as fixtures.

## Decision

* The examples are files, `examples/NN_name.nyx.json` at the repository root, found next to the program, in a parent folder of it or at `NYX_EXAMPLES` (`app/src/examples.rs`). They are written by `crates/nyx-engine/tests/make_examples.rs` (run it with `-- --ignored`), which builds each one **through engine commands** (`add_rect`, `add_overlay`, `add_physics`, `set_material`, `add_condition`, `add_pml_conditions`, ...), so an example is also a worked use of the interface. 16 examples: capacitor, cantilever, heated plate, thermal bar (heat feeding thermal expansion), block cantilever, WR-90 guide, loaded guide, sapphire cantilever (rotated crystal from the library), PML slab, plate with a hole, anchor loss (PML substrate), drilled plate, silicon strip on sapphire (optical and acoustic modes with PML, named materials), box cavity, SOI wire and its hybrid-mode variant.
* `crates/nyx-engine/tests/examples.rs` opens, checks, meshes and solves every example and compares the ones with a closed form (cantilever tip against Euler-Bernoulli, WR-90 cutoff, TE101 frequency of the box, the documented TE0 index of the SOI wire). These tests replace `usecases.rs`, `legacy_thermal.rs` and `legacy_contract.rs`.
* The two Bragg grating examples are not carried over: a swept layer with a `domain_offset` has domain numbers that carry no name, so a material cannot be assigned to them by name. Naming the domains of swept layers is open (a follow-up for the mesh owner).
* `femcore::testkit::soi` keeps the cross-section of the SOI wire as a fixture for the meshing and golden tests of other crates; `nyx-ems/tests/common` builds the typed problem from it.
* The benchmark programs (`bench`, `profile3d`, `mem_profile`, `cavity_bench`, `stage_log`, `eig_beam`) drive the physics through the `Physics` trait on generated meshes (`examples/common`); `complex_bench`, `soi_modes` and `bragg_gratings` are deleted with the use cases they studied.

## Consequences

* A change to an example is a regenerated file. The program ships with its `examples` and `materials` folders, not with devices in code.
* The window's New dialog and the start-from cards of the dashboard point at example files by position in the sorted list.

## Code

`examples/*.nyx.json`, `crates/nyx-engine/tests/{make_examples,examples}.rs`, `crates/nyx-engine/examples/*`, `crates/app/src/examples.rs`.
