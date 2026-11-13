# Feature list: code

C4 level 4 (code) of the component **Feature list** (`geom.featurelist`) in the container **geom** (crate `geom`). [Component view of the container](../../3-components/geom.md) · [All code pages](../README.md)

> The serialisable feature list (single source of truth for code and GUI), its operations and validation, parameter expressions, printing as Rust builder calls, and the error type.

**Decisions:** [ADR 2](../../../adr/0002-geometry-as-named-feature-list.md), [ADR 3](../../../adr/0003-geom-is-a-lean-separate-kernel.md), [ADR 63](../../../adr/0063-the-format-key-and-no-json-for-bulk-data.md)  
**Size:** 5 files, 1410 lines.

## Files

| File | Lines | What it says about itself |
|---|---:|---|
| `crates/geom/src/error.rs` | 42 | The error type of the crate: plain-language messages that name the features involved. Converts into `String` because the femcore APIs use `Result<_, String>`. |
| `crates/geom/src/expr.rs` | 337 | Parameter expressions: a number or a small safe formula over named parameters. The grammar mirrors nyx `femcore::expr` (numbers, `+ - * / ^ **`, parentheses, the functions `sin cos tan exp log ln log10 sqrt abs sinh cosh tanh atan min ma... |
| `crates/geom/src/lib.rs` | 23 | `geom`: a lean CSG-over-B-rep geometry kernel for the nyx FEM stack. |
| `crates/geom/src/model.rs` | 848 | The serialisable feature list (the single source of truth for code and GUI), its operations, validation, and the `Model` builder that only pushes `Op`s and returns `FeatureId` handles. Nothing here evaluates geometry. |
| `crates/geom/src/print.rs` | 160 | `FeatureList::to_rust`: prints the builder calls that recreate a list (the "Copy as Rust" of the GUI). One statement per feature, `let fN = m.method(...)` where `N` is the feature id; operations without a builder shorthand fall back to `... |

## What it uses

| Target | Modules | Items |
|---|---|---|
| [Queries and selectors](./queries.md) | select | Selector |

## Used by

| Component | Modules | Items |
|---|---|---|
| [engine: Commands](../engine/commands.md) | model | FeatureList, Op |
| [engine: Geometry operations](../engine/geometry.md) | error, expr, model | Expr, Feature, FeatureList, GeomError, Op, Origin, Param, check_name |
| [engine: Meshing seam](../engine/meshing.md) | model | FeatureList |
| [engine: Project](../engine/project.md) | model | FORMAT, FeatureList |
| [engine: Registry and facade](../engine/registry.md) | crate root |  |
| [engine: Store link](../engine/storeLink.md) | expr, model | Expr, Op |
| [femcore: Error estimators](../femcore/estimators.md) | model | FeatureList |
| [femcore: Fixtures and test kit](../femcore/fixtures.md) | model | FeatureId, Model |
| [femcore: Meshing facade](../femcore/meshing.md) | model | FeatureList, Model, RegionKind |
| [femcore: 3D meshing](../femcore/meshing3d.md) | model | Model, RegionKind |
| [femcore: Structured absorbing layers](../femcore/structpml.md) | model | RegionKind |
| [Evaluation](./evaluation.md) | error, expr, model | Env, Expr, GeomError, GeomResult |
| [2D kernel](./kernel2d.md) | error, model | FeatureId, GeomError, GeomResult, NameScheme, RegionKind |
| [3D kernel](./kernel3d.md) | error, expr, model | Axis, Env, Expr, GeomError, GeomResult, NameScheme, RegionKind |
| [Queries and selectors](./queries.md) | error, expr, model | Axis, Expr, GeomError, GeomResult, Op, RegionKind |

## Class diagram

Public types with their public fields (or variants) and methods, free functions per file as a `module` box; arrows are field types and `impl Trait`. Capped at 36 types and a few members each; the table below is complete.

```mermaid
classDiagram
  namespace error_rs {
    class GeomError {
        <<struct>>
        +_0 String
    }
    class String {
        <<impl block>>
    }
    class fns_error_rs {
        <<module error.rs>>
        +GeomResult type
    }
  }
  namespace expr_rs {
    class Expr {
        <<enumeration>>
        Num
        Text
        +source() String
    }
    class Env {
        <<struct>>
        +vars Vec
        +get(name) Option~f64~
        +eval(e, what) GeomResult~f64~
    }
    class fns_expr_rs {
        <<module expr.rs>>
        +eval_text(text, env) Result
    }
  }
  namespace model_rs {
    class Axis {
        <<enumeration>>
        X
        Y
        Z
        +index() usize
        +name() str
    }
    class Origin {
        <<enumeration>>
        Gui
        Code
    }
    class NameScheme {
        <<enumeration>>
        Feature
        Region
        Numbered
    }
    class RegionKind {
        <<enumeration>>
        Material
        Hole
        Pml
        Void
    }
    class Param {
        <<struct>>
        +name String
        +value Expr
    }
    class RegionDef {
        <<struct>>
        +name String
        +input FeatureId
        +kind RegionKind
        +pinned_id Option~i32~
    }
    class LabelPin {
        <<struct>>
        +sel Selector
        +id i32
    }
    class Side {
        <<enumeration>>
        XMin
        XMax
        YMin
        YMax
        ZMin
        ZMax
        +axis() usize
        +sign() i8
        +name() str
    }
    class TStep {
        <<enumeration>>
        Translate
        Rotate
        Mirror
        Scale
    }
    class PatternKind {
        <<enumeration>>
        Linear
        Polar
    }
    class Keep {
        <<enumeration>>
        Above
        Below
    }
    class Op {
        <<enumeration>>
        Rect
        Circle
        Ellipse
        Polygon
        RoundedRect
        Box
        Cylinder
        Sphere
        Cone
        Revolve
        Torus
        Extrude
        more_11_variants
        +kind() str
        +exprs() Vec
        +inputs() Vec~FeatureId~
    }
    class Feature {
        <<struct>>
        +id FeatureId
        +name String
        +op Op
        +suppressed bool
        +origin Origin
        +tags Vec~String~
    }
    class FeatureList {
        <<struct>>
        +format u32
        +dim u8
        +params Vec~Param~
        +features Vec~Feature~
        +new(dim) FeatureList
        +feature(id) Option
        +feature_by_name(name) Option
        +with(set) FeatureList
        +validate() GeomResult
        +to_rust() String
    }
    class Model {
        <<struct>>
        +list FeatureList
        +new_2d() Model
        +new_3d() Model
        +from_list(list) Model
        +set_origin(o) Self
        +param(name, value) Self
        +push(name, op) FeatureId
        +set_naming(scheme) Self
        +tag_feature(id, tag) Self
        +more_33_methods()
    }
    class fns_model_rs {
        <<module model.rs>>
        +check_name(kind, name) GeomResult
        +FeatureId type
        +FORMAT const
    }
  }
  namespace print_rs {
    class fns_print_rs {
        <<module print.rs>>
        +to_rust(l) String
    }
  }
  GeomError --> String
  Expr --> String
  Env --> String
  Param --> String
  Param --> Expr
  RegionDef --> String
  RegionDef --> RegionKind
  TStep --> Expr
  TStep --> Axis
  PatternKind --> Expr
  PatternKind --> Axis
  Op --> Expr
  Op --> String
  Op --> Axis
  Op --> TStep
  Op --> PatternKind
  Op --> RegionDef
  Op --> LabelPin
  Op --> NameScheme
  Op --> Side
  Op --> Keep
  Feature --> String
  Feature --> Op
  Feature --> Origin
  FeatureList --> Param
  FeatureList --> Feature
  Model --> FeatureList
```

## Public items

| Kind | Item | File:line | Signature or summary |
|---|---|---|---|
| struct | `GeomError` | `error.rs:7` | fields: `_0` |
| type | `GeomResult` | `error.rs:9` | `<T>` |
| enum | `Expr` | `expr.rs:11` | variants: `Num`, `Text` —  A number or an expression string; serialises as a bare JSON number or string (`untagged`). |
| fn | `eval_text` | `expr.rs:288` | `pub fn eval_text(text: &str, env: &dyn Fn(&str) -> Option<f64>) -> Result<f64, String>` —  Evaluates `text` with names resolved by `env`; the result must be finite. |
| struct | `Env` | `expr.rs:300` | fields: `vars` —  Named parameters evaluated in order; later parameters may use earlier ones. |
| fn | `Expr::source` | `expr.rs:44` | `pub fn source(&self) -> String` —  Source text for printing (numbers use Rust's shortest round-trip form). |
| fn | `Env::get` | `expr.rs:305` | `pub fn get(&self, name: &str) -> Option<f64>` |
| fn | `Env::eval` | `expr.rs:309` | `pub fn eval(&self, e: &Expr, what: &str) -> GeomResult<f64>` |
| type | `FeatureId` | `model.rs:9` |  |
| enum | `Axis` | `model.rs:21` | variants: `X`, `Y`, `Z` |
| enum | `Origin` | `model.rs:37` | variants: `Gui`, `Code` |
| enum | `NameScheme` | `model.rs:44` | variants: `Feature`, `Region`, `Numbered` —  How boundary faces (a region against the void) are named. Interfaces between two regions are always `A\|B`. |
| enum | `RegionKind` | `model.rs:55` | variants: `Material`, `Hole`, `Pml`, `Void` |
| struct | `Param` | `model.rs:66` | fields: `name`, `value` |
| struct | `RegionDef` | `model.rs:72` | fields: `name`, `input`, `kind`, `pinned_id` |
| struct | `LabelPin` | `model.rs:82` | fields: `sel`, `id` —  Pins the integer boundary label of everything a selector resolves to (all fragments share the int). |
| enum | `Side` | `model.rs:88` | variants: `XMin`, `XMax`, `YMin`, `YMax`, `ZMin`, `ZMax` |
| enum | `TStep` | `model.rs:114` | variants: `Translate`, `Rotate`, `Mirror`, `Scale` |
| enum | `PatternKind` | `model.rs:135` | variants: `Linear`, `Polar` |
| enum | `Keep` | `model.rs:148` | variants: `Above`, `Below` |
| enum | `Op` | `model.rs:155` | variants: `Rect`, `Circle`, `Ellipse`, `Polygon`, `RoundedRect`, `Box`, `Cylinder`, `Sphere`, `Cone`, `Revolve`, `Torus`, `Extrude`, ... |
| struct | `Feature` | `model.rs:381` | fields: `id`, `name`, `op`, `suppressed`, `origin`, `tags` |
| const | `FORMAT` | `model.rs:393` | `u32` —  The file format of a feature list (ADR 63). Written in every list; a list of another format is refused, there is no migration. |
| struct | `FeatureList` | `model.rs:400` | fields: `format`, `dim`, `params`, `features` |
| fn | `check_name` | `model.rs:413` | `pub fn check_name(kind: &str, name: &str) -> GeomResult<()>` |
| struct | `Model` | `model.rs:519` | fields: `list` —  Builder: pushes operations and returns handles. Holds no logic beyond name checks. |
| fn | `Axis::index` | `model.rs:28` | `pub fn index(self) -> usize` |
| fn | `Axis::name` | `model.rs:31` | `pub fn name(self) -> &'static str` |
| fn | `Side::axis` | `model.rs:98` | `pub fn axis(self) -> usize` |
| fn | `Side::sign` | `model.rs:101` | `pub fn sign(self) -> i8` |
| fn | `Side::name` | `model.rs:108` | `pub fn name(self) -> &'static str` |
| fn | `Op::kind` | `model.rs:296` | `pub fn kind(&self) -> &'static str` |
| fn | `Op::exprs` | `model.rs:325` | `pub fn exprs(&self) -> Vec<&Expr>` —  Every number-or-expression field of the operation, in a fixed order. |
| fn | `Op::inputs` | `model.rs:367` | `pub fn inputs(&self) -> Vec<FeatureId>` —  Features this operation refers to (all must come earlier in the list). |
| fn | `FeatureList::new` | `model.rs:426` | `pub fn new(dim: u8) -> FeatureList` |
| fn | `FeatureList::feature` | `model.rs:430` | `pub fn feature(&self, id: FeatureId) -> Option<&Feature>` |
| fn | `FeatureList::feature_by_name` | `model.rs:434` | `pub fn feature_by_name(&self, name: &str) -> Option<&Feature>` |
| fn | `FeatureList::with` | `model.rs:439` | `pub fn with(&self, set: &[(&str, f64)]) -> FeatureList` —  Parameter override without mutation: a copy with the listed parameters set to the given values. |
| fn | `FeatureList::validate` | `model.rs:451` | `pub fn validate(&self) -> GeomResult<()>` —  Structural checks: unique ids and names, references to earlier features only, one overlay, valid names. |
| fn | `FeatureList::to_rust` | `model.rs:512` | `pub fn to_rust(&self) -> String` —  Rust builder code that recreates this list (see `Model`). |
| fn | `Model::new_2d` | `model.rs:533` | `pub fn new_2d() -> Model` |
| fn | `Model::new_3d` | `model.rs:536` | `pub fn new_3d() -> Model` |
| fn | `Model::from_list` | `model.rs:539` | `pub fn from_list(list: FeatureList) -> Model` |
| fn | `Model::set_origin` | `model.rs:544` | `pub fn set_origin(&mut self, o: Origin) -> &mut Self` —  Features added after this call are marked as coming from the GUI. |
| fn | `Model::param` | `model.rs:549` | `pub fn param(&mut self, name: &str, value: impl Into<Expr>) -> &mut Self` |
| fn | `Model::push` | `model.rs:554` | `pub fn push(&mut self, name: &str, op: Op) -> FeatureId` |
| fn | `Model::set_naming` | `model.rs:569` | `pub fn set_naming(&mut self, scheme: NameScheme) -> &mut Self` —  Chooses how boundary faces are named (see `NameScheme`); applies to the model's overlay. |
| fn | `Model::tag_feature` | `model.rs:578` | `pub fn tag_feature(&mut self, id: FeatureId, tag: &str) -> &mut Self` |
| fn | `Model::suppress` | `model.rs:584` | `pub fn suppress(&mut self, id: FeatureId, on: bool) -> &mut Self` |
| fn | `Model::rect` | `model.rs:591` | `pub fn rect( &mut self, name: &str, min: (impl Into<Expr>, impl Into<Expr>), max: (impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` |
| fn | `Model::rounded_rect` | `model.rs:599` | `pub fn rounded_rect( &mut self, name: &str, min: (impl Into<Expr>, impl Into<Expr>), max: (impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::box_` | `model.rs:608` | `pub fn box_( &mut self, name: &str, min: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), max: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` |
| fn | `Model::cylinder` | `model.rs:616` | `pub fn cylinder( &mut self, name: &str, base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>, h: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::extrude_tapered` | `model.rs:625` | `pub fn extrude_tapered( &mut self, name: &str, profile: FeatureId, y0: impl Into<Expr>, y1: impl Into<Expr>, taper_deg: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::cylinder_along` | `model.rs:639` | `pub fn cylinder_along( &mut self, name: &str, axis: Axis, base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>, h: impl Into<Expr>, ) -> FeatureId` —  A cylinder along any axis (`axis` x, y or z); all prismatic solids of a model must share one axis. |
| fn | `Model::extrude_along` | `model.rs:650` | `pub fn extrude_along( &mut self, name: &str, axis: Axis, profile: FeatureId, a0: impl Into<Expr>, a1: impl Into<Expr>, ) -> FeatureId` —  A profile swept along any axis, from `a0` to `a1` on that axis (see `Op::Extrude`). |
| fn | `Model::cone` | `model.rs:660` | `pub fn cone( &mut self, name: &str, base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), r0: impl Into<Expr>, r1: impl Into<Expr>, h: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::cone_along` | `model.rs:670` | `pub fn cone_along( &mut self, name: &str, axis: Axis, base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), r0: impl Into<Expr>, r1: impl Into<Expr>, h: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::revolve` | `model.rs:682` | `pub fn revolve( &mut self, name: &str, profile: FeatureId, center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` —  Turns a `Rect` / `Polygon` profile once about the z axis through `center` (see `Op::Revolve`). |
| fn | `Model::revolve_along` | `model.rs:690` | `pub fn revolve_along( &mut self, name: &str, axis: Axis, profile: FeatureId, center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` |
| fn | `Model::sphere` | `model.rs:699` | `pub fn sphere( &mut self, name: &str, center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::torus` | `model.rs:707` | `pub fn torus( &mut self, name: &str, center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), major: impl Into<Expr>, minor: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::extrude` | `model.rs:717` | `pub fn extrude(&mut self, name: &str, profile: FeatureId, y0: impl Into<Expr>, y1: impl Into<Expr>) -> FeatureId` —  Sweeps an `(x, z)` profile along y from `y0` to `y1` (the wafer frame: cross-section in x-z, propagation along y). |
| fn | `Model::circle` | `model.rs:720` | `pub fn circle(&mut self, name: &str, center: (impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>) -> FeatureId` |
| fn | `Model::ellipse` | `model.rs:723` | `pub fn ellipse( &mut self, name: &str, center: (impl Into<Expr>, impl Into<Expr>), radii: (impl Into<Expr>, impl Into<Expr>), angle_deg: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::polygon` | `model.rs:732` | `pub fn polygon<A: Into<Expr> + Clone, B: Into<Expr> + Clone>(&mut self, name: &str, pts: &[(A, B)]) -> FeatureId` |
| fn | `Model::polygon_named` | `model.rs:736` | `pub fn polygon_named<A: Into<Expr> + Clone, B: Into<Expr> + Clone>( &mut self, name: &str, pts: &[(A, B)], edge_names: &[&str], ) -> FeatureId` |
| fn | `Model::transform` | `model.rs:745` | `pub fn transform(&mut self, name: &str, input: FeatureId, steps: Vec<TStep>) -> FeatureId` |
| fn | `Model::translate` | `model.rs:748` | `pub fn translate( &mut self, name: &str, input: FeatureId, d: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` |
| fn | `Model::linear_pattern` | `model.rs:756` | `pub fn linear_pattern( &mut self, name: &str, input: FeatureId, step: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), n: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::union` | `model.rs:765` | `pub fn union(&mut self, name: &str, ids: &[FeatureId]) -> FeatureId` |
| fn | `Model::difference` | `model.rs:768` | `pub fn difference(&mut self, name: &str, a: FeatureId, b: &[FeatureId]) -> FeatureId` |
| fn | `Model::intersect` | `model.rs:771` | `pub fn intersect(&mut self, name: &str, ids: &[FeatureId]) -> FeatureId` |
| fn | `Model::overlay` | `model.rs:775` | `pub fn overlay(&mut self, regions: &[(&str, FeatureId)]) -> FeatureId` —  Material regions in painter order (later wins). |
| fn | `Model::overlay_defs` | `model.rs:782` | `pub fn overlay_defs(&mut self, regions: Vec<RegionDef>, pins: Vec<LabelPin>) -> FeatureId` |
| fn | `Model::overlay_full` | `model.rs:786` | `pub fn overlay_full( &mut self, regions: &[(&str, FeatureId, RegionKind, Option<i32>)], pins: &[(&str, i32)], ) -> FeatureId` —  Regions with kinds and pinned ids, and boundary pins given as `(selector text, label)`. |
| fn | `Model::pml` | `model.rs:801` | `pub fn pml( &mut self, name: &str, window: FeatureId, sides: &[Side], thickness: impl Into<Expr>, layers: impl Into<Expr>, grading: impl Into<Expr>, ) -> FeatureId` |
| fn | `Model::cut` | `model.rs:822` | `pub fn cut(&mut self, name: &str, normal: (f64, f64, f64), offset: impl Into<Expr>, keep: Keep) -> FeatureId` |
| fn | `Model::periodic` | `model.rs:825` | `pub fn periodic( &mut self, name: &str, a: Selector, b: Selector, translation: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>), ) -> FeatureId` |
| fn | `Model::size` | `model.rs:834` | `pub fn size( &mut self, name: &str, sel: Selector, h: Option<Expr>, max_deviation: Option<Expr>, min_layers: Option<Expr>, ) -> FeatureId` |
| fn | `Model::tag` | `model.rs:844` | `pub fn tag(&mut self, name: &str, sel: Selector, tag: &str) -> FeatureId` |
| fn | `to_rust` | `print.rs:36` | `pub fn to_rust(l: &FeatureList) -> String` |

Trait implementations in this component: `GeomError: Display`, `GeomError: Error`, `String: From`, `GeomError: From`, `GeomError: From`, `Expr: From`, `Expr: From`, `Expr: From`, `Expr: From`, `Expr: From`
