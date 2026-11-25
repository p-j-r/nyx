# Project: code

C4 level 4 (code) of the component **Project** (`engine.project`) in the container **nyx-engine** (crate `nyx_engine`). [Component view of the container](../../3-components/engine.md) · [All code pages](../README.md)

> The saved state: feature list, mesh settings, one ProblemSpec per physics, project materials. The format key is checked at the file boundary.

**Decisions:** [ADR 63](../../../adr/0063-the-format-key-and-no-json-for-bulk-data.md)  
**Size:** 1 file, 318 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-engine/src/project.rs` | 318 |  |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [container `femcore` (femcore)](../../3-components/femcore.md) | meshing::request | PmlModeRequest |
| [container `geom` (geom)](../../3-components/geom.md) | model | FORMAT, FeatureList |
| [container `physics` (nyx-physics)](../../3-components/physics.md) | error, registry, spec, validate | Error, FORMAT, ProblemSpec, Registry, canonical_json, validate_spec |
| container `projectFile` (Project file) |  (declared: Saves and loads) |  |
| [Registry and facade](./registry.md) | crate root | physics_for |

## Used by

| Component | Modules | Items |
|---|---|---|
| [app: Inspector tabs](../app/inspector.md) | crate root, project | Auto, Inline, Structured, SweepLayer |
| [app: Window](../app/window.md) | project | MeshSettings |
| [Commands](./commands.md) | crate root, project | Auto, MeshSettings, Project |
| [Job runner and queue](./jobs.md) | project | Project |
| [Meshing seam](./meshing.md) | crate root, project | MeshSettings, PmlModeRequest, Project, SweepLayer |
| [Registry and facade](./registry.md) | crate root, project | MeshSettings, PmlModeRequest, Project, SweepLayer |
| [Results files](./resultsIO.md) | project | Project |
| [Store link](./storeLink.md) | project | Project |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class SweepLayer {
      <<struct>>
      +thickness f64
      +n usize
      +domain_offset i32
      +fill Option~String~
      +regions Option
  }
  class MeshSettings {
      <<struct>>
      +h f64
      +sweep Option
      +adapt_scale Option~f64~
      +physics_sizing bool
      +points_per_wavelength Option~f64~
      +pml_mode PmlModeRequest
      +structured_pml Option~bool~
      +max_dev Option~f64~
      +more_14_fields
      +resolved_pml_mode() PmlModeRequest
  }
  class Project {
      <<struct>>
      +format u32
      +title String
      +geometry FeatureList
      +mesh MeshSettings
      +physics Vec~ProblemSpec~
      +active Option~String~
      +materials BTreeMap
      +new(title, dim) Project
      +spec(physics) Option
      +spec_mut(physics) Option
      +need_spec(physics) Result
      +default_targets() Vec~String~
      +to_json() String
      +from_json(text) Result
      +check_formats() Result
      +more_1_methods()
  }
  class fns_project_rs {
      <<module project.rs>>
      +PROJECT_FORMAT const
  }
  MeshSettings --> SweepLayer
  Project --> MeshSettings
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `SweepLayer` | `project.rs:25` | fields: `thickness`, `n`, `domain_offset`, `fill`, `regions` —  One layer of a model swept along y into a 3D mesh (the 2D cross-section is extruded with `n` cells across `thickness`). |
| struct | `MeshSettings` | `project.rs:42` | fields: `h`, `sweep`, `adapt_scale`, `physics_sizing`, `points_per_wavelength`, `pml_mode`, `structured_pml`, `max_dev`, `curved`, `gap_cells`, `z_coarsen`, `adapt_grade`, ... —  How the geometry is meshed. |
| const | `PROJECT_FORMAT` | `project.rs:156` | `u32` —  The file format of a project (ADR 63). |
| struct | `Project` | `project.rs:163` | fields: `format`, `title`, `geometry`, `mesh`, `physics`, `active`, `materials` |
| fn | `MeshSettings::resolved_pml_mode` | `project.rs:146` | `pub fn resolved_pml_mode(&self) -> PmlModeRequest` —  Resolves effective [`PmlModeRequest`], honoring legacy `structured_pml: Some(true \| false)` |
| fn | `Project::new` | `project.rs:183` | `pub fn new(title: &str, dim: u8) -> Project` —  An empty project: no features, no physics. |
| fn | `Project::spec` | `project.rs:195` | `pub fn spec(&self, physics: &str) -> Option<&ProblemSpec>` |
| fn | `Project::spec_mut` | `project.rs:199` | `pub fn spec_mut(&mut self, physics: &str) -> Option<&mut ProblemSpec>` |
| fn | `Project::need_spec` | `project.rs:204` | `pub fn need_spec(&mut self, physics: &str) -> Result<&mut ProblemSpec, Error>` —  The spec of a physics the project must have, with the list of those it has in the message when it has not. |
| fn | `Project::default_targets` | `project.rs:215` | `pub fn default_targets(&self) -> Vec<String>` —  The codes to solve when none is named: the active physics, else all of them. |
| fn | `Project::to_json` | `project.rs:223` | `pub fn to_json(&self) -> String` —  The project as JSON text with sorted keys. |
| fn | `Project::from_json` | `project.rs:230` | `pub fn from_json(text: &str) -> Result<Project, Error>` —  Reads a project file. The `format` key must be present (in the project, the feature list and every problem) and supported: a file of |
| fn | `Project::check_formats` | `project.rs:253` | `pub fn check_formats(&self) -> Result<(), Error>` —  Refuses a format this program does not read, in the project, its feature list and its problems. |
| fn | `Project::check` | `project.rs:273` | `pub fn check(&self, registry: &Registry) -> Result<(), Error>` —  Everything that can be checked without a mesh: the feature list, the element size, the physics against the registry and each problem |

Trait implementations in this component: `MeshSettings: Default`
