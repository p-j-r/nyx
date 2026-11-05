# 58. The engine resolves named materials from the library and the project

Date: 2026-11-05

## Status

Accepted. Implements the wiring of the material files (ADR 29, 30, 34, 35) into the studies that phases 4 and 6 left to phase 7.

## Context

The envelope accepted `MaterialAssign::Named` since phase 4, but every physics answered "material files arrive with phase 6". Phase 6 built the library and the models (`femcore::material`); nothing used them. The window had its own presets and crystal editors, which the new window could not carry.

## Decision

* `MaterialAssign` has a third form, `Oriented { material, euler_deg }`: a named material placed in the wafer frame with Bunge angles. `Named`, `Oriented` and `Props` are told apart by their JSON shape (text, object with `material`, object of numbers).
* `Ctx` carries the material `library` (`femcore::material::Library`) and the free-space `wavelength` of the study. `nyx_physics::decode::materials_by_label` resolves a named material for the properties the physics lists: scalars under the descriptor names (`legacy_scalars`), the crystal data (`Anisotropy`) where the material is anisotropic or rotated, and for a physics that reads `eps_r` the permittivity from the optical model at the wavelength (isotropic: a scalar; birefringent: a tensor in the wafer frame). A material that is lossy at that wavelength is refused by the real-mode solver with a message that points to absorbing layers.
* The **engine** supplies both: `Engine::library()` is the folder library (`materials/` next to the program or in a parent folder, or `NYX_MATERIALS`) with the materials the **project** carries (`Project.materials`: name to TOML text, which may `extends` a library material; added with `add_material`, read at once so a mistake is reported there). The wavelength is `c / f` of the first study parameter of the physics with the unit Hz.
* `Physics::sizing` gets the library (`sizing(&self, spec, library)`) and `ems` and `sms` now state their wavelengths (`optical_wavelength`, `elastic_wavelength` in `femcore::meshing::physics_size`), so `MeshSettings::physics_sizing` works for named materials and anisotropic crystals. The mesher takes a `MeshRequest` (ADR 54): `FemcoreMesher` meshes 2D and 3D models, simple or adaptive, with structured absorbing layers and element sizes from the waves; a 2D model can still be swept along y.
* Commands: `list_materials`, `describe_material`, `add_material`, `remove_material`; `set_material` takes the three forms. `add_pml_conditions` writes one `pml` condition per absorbing piece of the geometry (the pieces and their sides come from the `Pml` feature), so the geometry and the physics cannot disagree.

## Consequences

* A person or a host picks silicon, silica or sapphire by name; an example of a rotated sapphire block (r-plane, Bunge angles) is a project file.
* The numbers of the examples differ slightly from the old window's: the library values were checked against the literature (ADR 35) and are not the earlier rounded ones. The goldens of the physics crates are unaffected (they use inline numbers).
* The `silicon()` and `sapphire()` helpers of `femcore::material` stay for tests of other crates; no program uses them.

## Code

`crates/nyx-physics/src/{spec,physics,decode}.rs`, `crates/nyx-engine/src/{materials,meshing,pml}.rs`, `crates/femcore/src/material/spec.rs` (`Library::with_texts`), `crates/femcore/src/meshing/physics_size.rs`, `crates/nyx-{ems,sms}/src/lib.rs` (`sizing`), `crates/nyx-engine/tests/{materials,examples}.rs`.
