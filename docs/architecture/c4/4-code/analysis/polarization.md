# Polarization: code

C4 level 4 (code) of the component **Polarization** (`analysis.polarization`) in the container **nyx-analysis** (crate `nyx_analysis`). [Component view of the container](../../3-components/analysis.md) · [All code pages](../README.md)

> Stokes parameters, Poincare sphere and the classification of an optical mode.

**Decisions:** [ADR 56](../../../adr/0056-polarization-analysis-is-a-function-of-a-result.md), [ADR 21](../../../adr/0021-polarization-analysis-from-complex-fields.md)  
**Size:** 1 file, 936 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/nyx-analysis/src/polarization.rs` | 936 | Polarization analysis of optical modes: Stokes parameters, Poincaré sphere and a classification of the mode. |

## What it uses

Nothing outside the component (the standard library and external crates are not listed).

## Used by

| Component | Modules | Items |
|---|---|---|
| [Mode report](./modes.md) | polarization |  |
| [Field samples](./samples.md) | polarization | Cx, Sample, TimeConvention |
| [engine: Analyses](../engine/analyse.md) | polarization |  |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  class TimeConvention {
      <<enumeration>>
      PlusJ
      MinusI
  }
  class Stokes {
      <<struct>>
      +s0 f64
      +s1 f64
      +s2 f64
      +s3 f64
      +from_fields(lat, nor, conv) Stokes
      +normalized() slice
      +dop() f64
      +psi_deg() f64
      +chi_deg() f64
  }
  class Sample {
      <<struct>>
      +lat Cx
      +nor Cx
      +prop Cx
      +w f64
  }
  class PolClass {
      <<enumeration>>
      QuasiTe
      QuasiTm
      LinearPlus45
      LinearMinus45
      Elliptical
      Circular
      Hybrid
      +label() str
  }
  class ClassOptions {
      <<struct>>
      +p_min f64
      +circ_min f64
      +ell_min f64
      +time TimeConvention
  }
  class ModeReport {
      <<struct>>
      +integrated Stokes
      +sphere slice
      +dop f64
      +chi_deg f64
      +psi_deg f64
      +te_fraction f64
      +long_fraction f64
      +class PolClass
      +summary_line() String
  }
  class SphereState {
      <<struct>>
      +label String
      +s slice
      +dashed bool
  }
  class WireLine {
      <<struct>>
      +pts Vec
      +front bool
      +width f64
      +opacity f64
  }
  class WireLabel {
      <<struct>>
      +text String
      +x f64
      +y f64
      +front bool
  }
  class SphereDrawing {
      <<struct>>
      +size f64
      +centre tuple
      +radius f64
      +lines Vec~WireLine~
      +labels Vec~WireLabel~
  }
  class ModeMarker {
      <<struct>>
      +label String
      +centre tuple
      +tip tuple
      +front bool
      +dashed bool
      +colour usize
  }
  class fns_polarization_rs {
      <<module polarization.rs>>
      +sphere_distance(a, b) f64
      +local_stokes(samples, conv) Vec~Stokes~
      +classify(integrated, opts) PolClass
      +analyze(samples, opts) Option~ModeReport~
      +degenerate_groups(neff, tol) Vec
      +p1_vertex_weights(points, cells, cell_labels, skip) Vec~f64~
      +poincare_svg(points, size) String
      +poincare_modes_svg(states, size) String
      +sphere_drawing(size) SphereDrawing
      +mode_markers(states, size) Vec~ModeMarker~
      +Cx type
      +MODE_COLOURS const
  }
  ClassOptions --> TimeConvention
  ModeReport --> Stokes
  ModeReport --> PolClass
  SphereDrawing --> WireLine
  SphereDrawing --> WireLabel
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| type | `Cx` | `polarization.rs:36` |  Complex number as [re, im]. |
| enum | `TimeConvention` | `polarization.rs:48` | variants: `PlusJ`, `MinusI` —  Time dependence of the complex field. |
| struct | `Stokes` | `polarization.rs:68` | fields: `s0`, `s1`, `s2`, `s3` |
| fn | `sphere_distance` | `polarization.rs:126` | `pub fn sphere_distance(a: [f64; 3], b: [f64; 3]) -> f64` —  Angle in radians between two points of the Poincaré sphere (pi for orthogonal polarizations). |
| struct | `Sample` | `polarization.rs:137` | fields: `lat`, `nor`, `prop`, `w` —  The field of a mode at one point of the cross-section and the area (or volume share) that point stands for. |
| enum | `PolClass` | `polarization.rs:146` | variants: `QuasiTe`, `QuasiTm`, `LinearPlus45`, `LinearMinus45`, `Elliptical`, `Circular`, `Hybrid` |
| struct | `ClassOptions` | `polarization.rs:183` | fields: `p_min`, `circ_min`, `ell_min`, `time` —  Limits of the classification, all on the integrated Stokes vector of the mode. |
| struct | `ModeReport` | `polarization.rs:200` | fields: `integrated`, `sphere`, `dop`, `chi_deg`, `psi_deg`, `te_fraction`, `long_fraction`, `class` |
| fn | `local_stokes` | `polarization.rs:233` | `pub fn local_stokes(samples: &[Sample], conv: TimeConvention) -> Vec<Stokes>` —  Stokes parameters of every sample, in the order of `samples` (for plots on the sphere and maps of the cross-section). |
| fn | `classify` | `polarization.rs:237` | `pub fn classify(integrated: &Stokes, opts: &ClassOptions) -> PolClass` |
| fn | `analyze` | `polarization.rs:264` | `pub fn analyze(samples: &[Sample], opts: &ClassOptions) -> Option<ModeReport>` —  The report of one mode; `None` if the weighted intensity is zero or not finite. |
| fn | `degenerate_groups` | `polarization.rs:292` | `pub fn degenerate_groups(neff: &[(f64, f64)], tol: f64) -> Vec<Vec<usize>>` —  Groups of modes whose effective indices differ by at most `tol` (consecutive in sorted order of the real part). The |
| fn | `p1_vertex_weights` | `polarization.rs:322` | `pub fn p1_vertex_weights(points: &[f64], cells: &[u32], cell_labels: &[i32], skip: &[i32]) -> Vec<f64>` —  Area weight of every vertex of a 2D triangle mesh (a third of each adjacent triangle), skipping triangles whose |
| fn | `poincare_svg` | `polarization.rs:483` | `pub fn poincare_svg(points: &[(f64, f64, f64, f64)], size: f64) -> String` —  A Poincar\u{e9} sphere picture in the style of a Bloch sphere (translucent ball, light grid, bold equator, axes with |
| struct | `SphereState` | `polarization.rs:528` | fields: `label`, `s`, `dashed` —  One mode on the sphere: a label, its (s1, s2, s3) and whether it belongs to a degenerate group (drawn dashed |
| fn | `poincare_modes_svg` | `polarization.rs:536` | `pub fn poincare_modes_svg(states: &[SphereState], size: f64) -> String` —  The same sphere with one numbered marker per mode on a line from the centre to its integrated state, like a Bloch vector: a fully |
| struct | `WireLine` | `polarization.rs:578` | fields: `pts`, `front`, `width`, `opacity` —  A polyline of the sphere picture in picture coordinates (x to the right, y down, the picture is a square of side `size`). |
| struct | `WireLabel` | `polarization.rs:588` | fields: `text`, `x`, `y`, `front` —  A label at a point of the picture. |
| struct | `SphereDrawing` | `polarization.rs:597` | fields: `size`, `centre`, `radius`, `lines`, `labels` —  The empty Poincare sphere as lines and labels, for programs that draw it themselves (the SVG functions draw the same picture). |
| struct | `ModeMarker` | `polarization.rs:607` | fields: `label`, `centre`, `tip`, `front`, `dashed`, `colour` —  One mode on the picture: its tip, the colour (index into `MODE_COLOURS`) and whether it ends on the surface (`len` = 1) or inside. |
| const | `MODE_COLOURS` | `polarization.rs:617` | `[&str` —  The colours of the mode markers, in the order of the numbers (as `#rrggbb`). |
| fn | `sphere_drawing` | `polarization.rs:645` | `pub fn sphere_drawing(size: f64) -> SphereDrawing` —  The picture of the empty sphere in a square of side `size`. |
| fn | `mode_markers` | `polarization.rs:683` | `pub fn mode_markers(states: &[SphereState], size: f64) -> Vec<ModeMarker>` —  The markers of `states` on the picture of `sphere_drawing(size)`; states with a non-finite coordinate are skipped. |
| fn | `Stokes::from_fields` | `polarization.rs:76` | `pub fn from_fields(lat: Cx, nor: Cx, conv: TimeConvention) -> Stokes` |
| fn | `Stokes::normalized` | `polarization.rs:94` | `pub fn normalized(&self) -> [f64` —  (s1, s2, s3) divided by S0; zeros when there is no intensity. A point of a single sample has unit length. |
| fn | `Stokes::dop` | `polarization.rs:103` | `pub fn dop(&self) -> f64` —  Degree of polarization in [0, 1]. |
| fn | `Stokes::psi_deg` | `polarization.rs:109` | `pub fn psi_deg(&self) -> f64` —  Orientation of the major axis in degrees, in (-90, 90], measured from x towards y. |
| fn | `Stokes::chi_deg` | `polarization.rs:115` | `pub fn chi_deg(&self) -> f64` —  Ellipticity angle in degrees, in [-45, 45]: 0 linear, +45 CCW circular, -45 CW circular. Uses the direction of |
| fn | `PolClass::label` | `polarization.rs:166` | `pub fn label(&self) -> &'static str` |
| fn | `ModeReport::summary_line` | `polarization.rs:216` | `pub fn summary_line(&self) -> String` |

Trait implementations in this component: `ClassOptions: Default`
