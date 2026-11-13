//! The serialisable feature list (the single source of truth for code and GUI), its operations, validation, and the
//! `Model` builder that only pushes `Op`s and returns `FeatureId` handles. Nothing here evaluates geometry.

use crate::error::{GeomError, GeomResult};
use crate::expr::Expr;
use crate::select::Selector;
use serde::{Deserialize, Serialize};

pub type FeatureId = u32;

fn z_axis() -> Axis {
    Axis::Z
}

/// The default sweep of an `Extrude` is y, the propagation direction of the wafer frame.
fn y_axis() -> Axis {
    Axis::Y
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    pub fn index(self) -> usize {
        self as usize
    }
    pub fn name(self) -> &'static str {
        ["x", "y", "z"][self as usize]
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Origin {
    Gui,
    Code,
}

/// How boundary faces (a region against the void) are named. Interfaces between two regions are always `A|B`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum NameScheme {
    /// `<feature>.<role>`: after the shape that drew the face (`box.bottom`)
    #[default]
    Feature,
    /// `<region>.<role>`: after the material / region it bounds (`BOX.bottom`)
    Region,
    /// `<region>.<n>`: the region name and a running number in a deterministic order (`BOX.3`); needs no knowledge of roles
    Numbered,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum RegionKind {
    Material,
    /// removes its cells (the old "domain 0 or negative")
    Hole,
    /// a normal region that the physics layer treats as an absorbing layer
    Pml,
    /// the implicit outside, listed explicitly (behaves like a hole)
    Void,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Param {
    pub name: String,
    pub value: Expr,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RegionDef {
    pub name: String,
    pub input: FeatureId,
    pub kind: RegionKind,
    /// keep an integer domain label from a legacy project
    pub pinned_id: Option<i32>,
}

/// Pins the integer boundary label of everything a selector resolves to (all fragments share the int).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LabelPin {
    pub sel: Selector,
    pub id: i32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Side {
    XMin,
    XMax,
    YMin,
    YMax,
    ZMin,
    ZMax,
}

impl Side {
    pub fn axis(self) -> usize {
        self as usize / 2
    }
    pub fn sign(self) -> i8 {
        if (self as usize).is_multiple_of(2) {
            -1
        } else {
            1
        }
    }
    pub fn name(self) -> &'static str {
        ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"][self as usize]
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TStep {
    Translate([Expr; 3]),
    /// rotation by `deg` degrees about the axis through the point `about` (2D: axis Z)
    Rotate {
        axis: Axis,
        deg: Expr,
        about: [Expr; 3],
    },
    /// reflection across the plane `axis = at`
    Mirror {
        axis: Axis,
        at: Expr,
    },
    /// uniform scale about the point `about`
    Scale {
        factor: Expr,
        about: [Expr; 3],
    },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum PatternKind {
    Linear {
        step: [Expr; 3],
    },
    /// `total_deg` split evenly over `n` instances (360 for a full ring)
    Polar {
        axis: Axis,
        center: [Expr; 3],
        total_deg: Expr,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Keep {
    /// the side the normal points to
    Above,
    Below,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Op {
    Rect {
        min: [Expr; 2],
        max: [Expr; 2],
    },
    Circle {
        center: [Expr; 2],
        r: Expr,
    },
    Ellipse {
        center: [Expr; 2],
        radii: [Expr; 2],
        angle_deg: Expr,
    },
    /// counter-clockwise or clockwise; `names` gives the roles of the edges `i -> i+1` (default `e0 e1 ...`)
    Polygon {
        pts: Vec<[Expr; 2]>,
        names: Vec<String>,
    },
    RoundedRect {
        min: [Expr; 2],
        max: [Expr; 2],
        r: Expr,
    },
    /// 3D: axis-aligned box (roles `xmin xmax ymin ymax zmin zmax`)
    Box {
        min: [Expr; 3],
        max: [Expr; 3],
    },
    /// 3D: cylinder from `base` along +axis (this build: axis z only); roles `side bottom top`
    Cylinder {
        base: [Expr; 3],
        axis: Axis,
        r: Expr,
        h: Expr,
    },
    /// 3D: sphere (one face `surface`; must not touch any other face)
    Sphere {
        center: [Expr; 3],
        r: Expr,
    },
    /// 3D: cone frustum along the sweep axis from `base`, radius `r0` at the bottom and `r1` at the top (either may be 0); faces `side`,
    /// `bottom`, `top`, circle edges where they meet. Must lie clear of every other face
    Cone {
        base: [Expr; 3],
        #[serde(default = "z_axis")]
        axis: Axis,
        r0: Expr,
        r1: Expr,
        h: Expr,
    },
    /// 3D: a full turn of a `Rect` or `Polygon` profile about the sweep axis through `center`. The profile's `(x, y)` are the radius
    /// (must be >= 0) and the height along the axis. Faces `swept:<profile role>`: cone sides, cylinder walls and flat rings, with a
    /// circle edge at every profile corner off the axis. Must lie clear of every other face.
    Revolve {
        profile: FeatureId,
        #[serde(default = "z_axis")]
        axis: Axis,
        center: [Expr; 3],
    },
    /// 3D: torus with axis z (one face `surface`; must not touch any other face)
    Torus {
        center: [Expr; 3],
        major: Expr,
        minor: Expr,
    },
    /// 3D: a 2D profile feature swept along `axis` (default y, the propagation direction) from `from` to `to`; roles
    /// `side:<profile role>` and the caps `ymin`, `ymax` (`bottom`, `top` for a sweep along z). The profile's two coordinates are
    /// the other axes: `(x, z)` for a sweep along y, `(x, y)` for z and `(y, z)` for x.
    Extrude {
        profile: FeatureId,
        #[serde(default = "y_axis")]
        axis: Axis,
        #[serde(alias = "z0")]
        from: Expr,
        #[serde(alias = "z1")]
        to: Expr,
        /// wall angle from the sweep direction in degrees (positive narrows towards `to`); makes the solid polyhedral
        #[serde(default)]
        taper_deg: Option<Expr>,
    },
    Transform {
        input: FeatureId,
        steps: Vec<TStep>,
    },
    Pattern {
        input: FeatureId,
        kind: PatternKind,
        n: Expr,
    },
    Union {
        ids: Vec<FeatureId>,
    },
    Difference {
        a: FeatureId,
        b: Vec<FeatureId>,
    },
    Intersect {
        ids: Vec<FeatureId>,
    },
    /// the root: region order is painter priority, later wins
    Overlay {
        regions: Vec<RegionDef>,
        pins: Vec<LabelPin>,
        #[serde(default)]
        naming: NameScheme,
    },
    /// absorbing layers of `thickness` around the window `window` on `sides` (empty: all); splits every matching region
    /// (`regions` empty: all) into the physical part and one `<region>.pml_<sides>` piece per side combination
    Pml {
        window: FeatureId,
        sides: Vec<Side>,
        thickness: Expr,
        layers: Expr,
        grading: Expr,
        regions: Vec<String>,
    },
    /// keeps one side of the plane `normal . p = offset`; the cut face is named `<feature>.cut` and tagged `symmetry`
    Cut {
        normal: [Expr; 3],
        offset: Expr,
        keep: Keep,
    },
    Periodic {
        sel_a: Selector,
        sel_b: Selector,
        translation: [Expr; 3],
    },
    Size {
        sel: Selector,
        h: Option<Expr>,
        max_deviation: Option<Expr>,
        min_layers: Option<Expr>,
    },
    Tag {
        sel: Selector,
        name: String,
    },
}

impl Op {
    pub fn kind(&self) -> &'static str {
        match self {
            Op::Rect { .. } => "Rect",
            Op::Circle { .. } => "Circle",
            Op::Ellipse { .. } => "Ellipse",
            Op::Polygon { .. } => "Polygon",
            Op::RoundedRect { .. } => "RoundedRect",
            Op::Box { .. } => "Box",
            Op::Sphere { .. } => "Sphere",
            Op::Cone { .. } => "Cone",
            Op::Revolve { .. } => "Revolve",
            Op::Torus { .. } => "Torus",
            Op::Cylinder { .. } => "Cylinder",
            Op::Extrude { .. } => "Extrude",
            Op::Transform { .. } => "Transform",
            Op::Pattern { .. } => "Pattern",
            Op::Union { .. } => "Union",
            Op::Difference { .. } => "Difference",
            Op::Intersect { .. } => "Intersect",
            Op::Overlay { .. } => "Overlay",
            Op::Pml { .. } => "Pml",
            Op::Cut { .. } => "Cut",
            Op::Periodic { .. } => "Periodic",
            Op::Size { .. } => "Size",
            Op::Tag { .. } => "Tag",
        }
    }

    /// Every number-or-expression field of the operation, in a fixed order.
    pub fn exprs(&self) -> Vec<&Expr> {
        match self {
            Op::Rect { min, max } => min.iter().chain(max).collect(),
            Op::RoundedRect { min, max, r } => min.iter().chain(max).chain([r]).collect(),
            Op::Circle { center, r } => center.iter().chain([r]).collect(),
            Op::Box { min, max } => min.iter().chain(max).collect(),
            Op::Sphere { center, r } => center.iter().chain([r]).collect(),
            Op::Cone { base, r0, r1, h, .. } => base.iter().chain([r0, r1, h]).collect(),
            Op::Revolve { center, .. } => center.iter().collect(),
            Op::Torus { center, major, minor } => center.iter().chain([major, minor]).collect(),
            Op::Cylinder { base, r, h, .. } => base.iter().chain([r, h]).collect(),
            Op::Extrude { from: z0, to: z1, taper_deg, .. } => {
                [Some(z0), Some(z1), taper_deg.as_ref()].into_iter().flatten().collect()
            }
            Op::Ellipse { center, radii, angle_deg } => center.iter().chain(radii).chain([angle_deg]).collect(),
            Op::Polygon { pts, .. } => pts.iter().flatten().collect(),
            Op::Transform { steps, .. } => steps
                .iter()
                .flat_map(|s| match s {
                    TStep::Translate(v) => v.iter().collect::<Vec<_>>(),
                    TStep::Rotate { deg, about, .. } => about.iter().chain([deg]).collect(),
                    TStep::Mirror { at, .. } => vec![at],
                    TStep::Scale { factor, about } => about.iter().chain([factor]).collect(),
                })
                .collect(),
            Op::Pattern { kind, n, .. } => match kind {
                PatternKind::Linear { step } => step.iter().chain([n]).collect(),
                PatternKind::Polar { center, total_deg, .. } => center.iter().chain([total_deg, n]).collect(),
            },
            Op::Pml { thickness, layers, grading, .. } => vec![thickness, layers, grading],
            Op::Cut { normal, offset, .. } => normal.iter().chain([offset]).collect(),
            Op::Periodic { translation, .. } => translation.iter().collect(),
            Op::Size { h, max_deviation, min_layers, .. } => {
                [h, max_deviation, min_layers].into_iter().flatten().collect()
            }
            Op::Union { .. } | Op::Difference { .. } | Op::Intersect { .. } | Op::Overlay { .. } | Op::Tag { .. } => {
                vec![]
            }
        }
    }

    /// Features this operation refers to (all must come earlier in the list).
    pub fn inputs(&self) -> Vec<FeatureId> {
        match self {
            Op::Transform { input, .. } | Op::Pattern { input, .. } => vec![*input],
            Op::Extrude { profile, .. } | Op::Revolve { profile, .. } => vec![*profile],
            Op::Union { ids } | Op::Intersect { ids } => ids.clone(),
            Op::Difference { a, b } => std::iter::once(*a).chain(b.iter().copied()).collect(),
            Op::Overlay { regions, .. } => regions.iter().map(|r| r.input).collect(),
            Op::Pml { window, .. } => vec![*window],
            _ => vec![],
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Feature {
    pub id: FeatureId,
    pub name: String,
    pub op: Op,
    #[serde(default)]
    pub suppressed: bool,
    pub origin: Origin,
    #[serde(default)]
    pub tags: Vec<String>,
}

/// The file format of a feature list (ADR 63). Written in every list; a list of another format is refused, there is no migration.
pub const FORMAT: u32 = 1;

fn current_format() -> u32 {
    FORMAT
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FeatureList {
    /// the file format, `1` (ADR 63)
    #[serde(default = "current_format")]
    pub format: u32,
    /// 2 or 3
    pub dim: u8,
    pub params: Vec<Param>,
    pub features: Vec<Feature>,
}

/// Characters that would be ambiguous in auto-names (`feature.role`, `A|B`, `a/b`) or in the selector text form.
const FORBIDDEN: &str = ".|/,()[]&!:=\"'<>*";

pub fn check_name(kind: &str, name: &str) -> GeomResult<()> {
    if name.trim().is_empty() {
        return Err(GeomError(format!("a {kind} has an empty name")));
    }
    if name != name.trim() || name.chars().any(|c| c.is_control() || FORBIDDEN.contains(c)) {
        return Err(GeomError(format!(
            "{kind} name '{name}' contains a character from {FORBIDDEN:?}, a control character or edge spaces; those are reserved for auto-names and selectors"
        )));
    }
    Ok(())
}

impl FeatureList {
    pub fn new(dim: u8) -> FeatureList {
        FeatureList { format: FORMAT, dim, params: vec![], features: vec![] }
    }

    pub fn feature(&self, id: FeatureId) -> Option<&Feature> {
        self.features.iter().find(|f| f.id == id)
    }

    pub fn feature_by_name(&self, name: &str) -> Option<&Feature> {
        self.features.iter().find(|f| f.name == name)
    }

    /// Parameter override without mutation: a copy with the listed parameters set to the given values.
    pub fn with(&self, set: &[(&str, f64)]) -> FeatureList {
        let mut l = self.clone();
        for (n, v) in set {
            match l.params.iter_mut().find(|p| p.name == *n) {
                Some(p) => p.value = Expr::Num(*v),
                None => l.params.push(Param { name: n.to_string(), value: Expr::Num(*v) }),
            }
        }
        l
    }

    /// Structural checks: unique ids and names, references to earlier features only, one overlay, valid names.
    pub fn validate(&self) -> GeomResult<()> {
        if self.format != FORMAT {
            return Err(GeomError(format!(
                "feature list format {} is not supported (this program reads format {FORMAT})",
                self.format
            )));
        }
        if self.dim != 2 && self.dim != 3 {
            return Err(GeomError(format!("dimension must be 2 or 3, not {}", self.dim)));
        }
        let mut seen: Vec<(FeatureId, &str)> = vec![];
        for p in &self.params {
            check_name("parameter", &p.name)?;
        }
        let mut overlays = 0;
        for f in &self.features {
            check_name("feature", &f.name)?;
            if seen.iter().any(|(i, _)| *i == f.id) {
                return Err(GeomError(format!("feature id {} is used twice", f.id)));
            }
            if let Some((_, n)) = seen.iter().find(|(_, n)| *n == f.name) {
                return Err(GeomError(format!("two features are named '{n}'; names must be unique")));
            }
            for i in f.op.inputs() {
                if !seen.iter().any(|(s, _)| *s == i) {
                    return Err(GeomError(format!(
                        "feature '{}' refers to feature id {i}, which does not come before it in the list",
                        f.name
                    )));
                }
            }
            for t in &f.tags {
                check_name("tag", t)?;
            }
            match &f.op {
                Op::Overlay { regions, .. } => {
                    overlays += 1;
                    for (k, r) in regions.iter().enumerate() {
                        check_name("region", &r.name)?;
                        if regions[..k].iter().any(|o| o.name == r.name) {
                            return Err(GeomError(format!("two regions are named '{}'", r.name)));
                        }
                    }
                }
                Op::Tag { name, .. } => check_name("tag", name)?,
                Op::Polygon { names, .. } => {
                    for n in names {
                        check_name("polygon edge", n)?;
                    }
                }
                _ => {}
            }
            seen.push((f.id, &f.name));
        }
        if overlays > 1 {
            return Err(GeomError("a feature list may have only one Overlay".into()));
        }
        Ok(())
    }

    /// Rust builder code that recreates this list (see `Model`).
    pub fn to_rust(&self) -> String {
        crate::print::to_rust(self)
    }
}

/// Builder: pushes operations and returns handles. Holds no logic beyond name checks.
#[derive(Clone, Debug)]
pub struct Model {
    pub list: FeatureList,
    next_id: FeatureId,
    origin: Origin,
}

fn e2<A: Into<Expr>, B: Into<Expr>>(p: (A, B)) -> [Expr; 2] {
    [p.0.into(), p.1.into()]
}
fn e3<A: Into<Expr>, B: Into<Expr>, C: Into<Expr>>(p: (A, B, C)) -> [Expr; 3] {
    [p.0.into(), p.1.into(), p.2.into()]
}

impl Model {
    pub fn new_2d() -> Model {
        Model { list: FeatureList::new(2), next_id: 1, origin: Origin::Code }
    }
    pub fn new_3d() -> Model {
        Model { list: FeatureList::new(3), next_id: 1, origin: Origin::Code }
    }
    pub fn from_list(list: FeatureList) -> Model {
        let next_id = list.features.iter().map(|f| f.id).max().unwrap_or(0) + 1;
        Model { list, next_id, origin: Origin::Code }
    }
    /// Features added after this call are marked as coming from the GUI.
    pub fn set_origin(&mut self, o: Origin) -> &mut Self {
        self.origin = o;
        self
    }

    pub fn param(&mut self, name: &str, value: impl Into<Expr>) -> &mut Self {
        self.list.params.push(Param { name: name.to_string(), value: value.into() });
        self
    }

    pub fn push(&mut self, name: &str, op: Op) -> FeatureId {
        let id = self.next_id;
        self.next_id += 1;
        self.list.features.push(Feature {
            id,
            name: name.to_string(),
            op,
            suppressed: false,
            origin: self.origin,
            tags: vec![],
        });
        id
    }

    /// Chooses how boundary faces are named (see `NameScheme`); applies to the model's overlay.
    pub fn set_naming(&mut self, scheme: NameScheme) -> &mut Self {
        for f in &mut self.list.features {
            if let Op::Overlay { naming, .. } = &mut f.op {
                *naming = scheme;
            }
        }
        self
    }

    pub fn tag_feature(&mut self, id: FeatureId, tag: &str) -> &mut Self {
        if let Some(f) = self.list.features.iter_mut().find(|f| f.id == id) {
            f.tags.push(tag.to_string());
        }
        self
    }
    pub fn suppress(&mut self, id: FeatureId, on: bool) -> &mut Self {
        if let Some(f) = self.list.features.iter_mut().find(|f| f.id == id) {
            f.suppressed = on;
        }
        self
    }

    pub fn rect(
        &mut self,
        name: &str,
        min: (impl Into<Expr>, impl Into<Expr>),
        max: (impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Rect { min: e2(min), max: e2(max) })
    }
    pub fn rounded_rect(
        &mut self,
        name: &str,
        min: (impl Into<Expr>, impl Into<Expr>),
        max: (impl Into<Expr>, impl Into<Expr>),
        r: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::RoundedRect { min: e2(min), max: e2(max), r: r.into() })
    }
    pub fn box_(
        &mut self,
        name: &str,
        min: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        max: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Box { min: e3(min), max: e3(max) })
    }
    pub fn cylinder(
        &mut self,
        name: &str,
        base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        r: impl Into<Expr>,
        h: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Cylinder { base: e3(base), axis: Axis::Z, r: r.into(), h: h.into() })
    }
    pub fn extrude_tapered(
        &mut self,
        name: &str,
        profile: FeatureId,
        y0: impl Into<Expr>,
        y1: impl Into<Expr>,
        taper_deg: impl Into<Expr>,
    ) -> FeatureId {
        self.push(
            name,
            Op::Extrude { profile, axis: Axis::Y, from: y0.into(), to: y1.into(), taper_deg: Some(taper_deg.into()) },
        )
    }
    /// A cylinder along any axis (`axis` x, y or z); all prismatic solids of a model must share one axis.
    pub fn cylinder_along(
        &mut self,
        name: &str,
        axis: Axis,
        base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        r: impl Into<Expr>,
        h: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Cylinder { base: e3(base), axis, r: r.into(), h: h.into() })
    }
    /// A profile swept along any axis, from `a0` to `a1` on that axis (see `Op::Extrude`).
    pub fn extrude_along(
        &mut self,
        name: &str,
        axis: Axis,
        profile: FeatureId,
        a0: impl Into<Expr>,
        a1: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Extrude { profile, axis, from: a0.into(), to: a1.into(), taper_deg: None })
    }
    pub fn cone(
        &mut self,
        name: &str,
        base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        r0: impl Into<Expr>,
        r1: impl Into<Expr>,
        h: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Cone { base: e3(base), axis: Axis::Z, r0: r0.into(), r1: r1.into(), h: h.into() })
    }
    pub fn cone_along(
        &mut self,
        name: &str,
        axis: Axis,
        base: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        r0: impl Into<Expr>,
        r1: impl Into<Expr>,
        h: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Cone { base: e3(base), axis, r0: r0.into(), r1: r1.into(), h: h.into() })
    }
    /// Turns a `Rect` / `Polygon` profile once about the z axis through `center` (see `Op::Revolve`).
    pub fn revolve(
        &mut self,
        name: &str,
        profile: FeatureId,
        center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Revolve { profile, axis: Axis::Z, center: e3(center) })
    }
    pub fn revolve_along(
        &mut self,
        name: &str,
        axis: Axis,
        profile: FeatureId,
        center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Revolve { profile, axis, center: e3(center) })
    }
    pub fn sphere(
        &mut self,
        name: &str,
        center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        r: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Sphere { center: e3(center), r: r.into() })
    }
    pub fn torus(
        &mut self,
        name: &str,
        center: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        major: impl Into<Expr>,
        minor: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Torus { center: e3(center), major: major.into(), minor: minor.into() })
    }
    /// Sweeps an `(x, z)` profile along y from `y0` to `y1` (the wafer frame: cross-section in x-z, propagation along y).
    pub fn extrude(&mut self, name: &str, profile: FeatureId, y0: impl Into<Expr>, y1: impl Into<Expr>) -> FeatureId {
        self.push(name, Op::Extrude { profile, axis: Axis::Y, from: y0.into(), to: y1.into(), taper_deg: None })
    }
    pub fn circle(&mut self, name: &str, center: (impl Into<Expr>, impl Into<Expr>), r: impl Into<Expr>) -> FeatureId {
        self.push(name, Op::Circle { center: e2(center), r: r.into() })
    }
    pub fn ellipse(
        &mut self,
        name: &str,
        center: (impl Into<Expr>, impl Into<Expr>),
        radii: (impl Into<Expr>, impl Into<Expr>),
        angle_deg: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Ellipse { center: e2(center), radii: e2(radii), angle_deg: angle_deg.into() })
    }
    pub fn polygon<A: Into<Expr> + Clone, B: Into<Expr> + Clone>(&mut self, name: &str, pts: &[(A, B)]) -> FeatureId {
        let pts = pts.iter().cloned().map(e2).collect();
        self.push(name, Op::Polygon { pts, names: vec![] })
    }
    pub fn polygon_named<A: Into<Expr> + Clone, B: Into<Expr> + Clone>(
        &mut self,
        name: &str,
        pts: &[(A, B)],
        edge_names: &[&str],
    ) -> FeatureId {
        let pts = pts.iter().cloned().map(e2).collect();
        self.push(name, Op::Polygon { pts, names: edge_names.iter().map(|s| s.to_string()).collect() })
    }
    pub fn transform(&mut self, name: &str, input: FeatureId, steps: Vec<TStep>) -> FeatureId {
        self.push(name, Op::Transform { input, steps })
    }
    pub fn translate(
        &mut self,
        name: &str,
        input: FeatureId,
        d: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Transform { input, steps: vec![TStep::Translate(e3(d))] })
    }
    pub fn linear_pattern(
        &mut self,
        name: &str,
        input: FeatureId,
        step: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
        n: impl Into<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Pattern { input, kind: PatternKind::Linear { step: e3(step) }, n: n.into() })
    }
    pub fn union(&mut self, name: &str, ids: &[FeatureId]) -> FeatureId {
        self.push(name, Op::Union { ids: ids.to_vec() })
    }
    pub fn difference(&mut self, name: &str, a: FeatureId, b: &[FeatureId]) -> FeatureId {
        self.push(name, Op::Difference { a, b: b.to_vec() })
    }
    pub fn intersect(&mut self, name: &str, ids: &[FeatureId]) -> FeatureId {
        self.push(name, Op::Intersect { ids: ids.to_vec() })
    }
    /// Material regions in painter order (later wins).
    pub fn overlay(&mut self, regions: &[(&str, FeatureId)]) -> FeatureId {
        let regions = regions
            .iter()
            .map(|(n, i)| RegionDef { name: n.to_string(), input: *i, kind: RegionKind::Material, pinned_id: None })
            .collect();
        self.push("overlay", Op::Overlay { regions, pins: vec![], naming: NameScheme::Feature })
    }
    pub fn overlay_defs(&mut self, regions: Vec<RegionDef>, pins: Vec<LabelPin>) -> FeatureId {
        self.push("overlay", Op::Overlay { regions, pins, naming: NameScheme::Feature })
    }
    /// Regions with kinds and pinned ids, and boundary pins given as `(selector text, label)`.
    pub fn overlay_full(
        &mut self,
        regions: &[(&str, FeatureId, RegionKind, Option<i32>)],
        pins: &[(&str, i32)],
    ) -> FeatureId {
        let regions = regions
            .iter()
            .map(|(n, i, k, p)| RegionDef { name: n.to_string(), input: *i, kind: *k, pinned_id: *p })
            .collect();
        let pins = pins
            .iter()
            .map(|(s, id)| LabelPin { sel: Selector::parse(s).expect("a valid selector in overlay_full"), id: *id })
            .collect();
        self.push("overlay", Op::Overlay { regions, pins, naming: NameScheme::Feature })
    }
    pub fn pml(
        &mut self,
        name: &str,
        window: FeatureId,
        sides: &[Side],
        thickness: impl Into<Expr>,
        layers: impl Into<Expr>,
        grading: impl Into<Expr>,
    ) -> FeatureId {
        self.push(
            name,
            Op::Pml {
                window,
                sides: sides.to_vec(),
                thickness: thickness.into(),
                layers: layers.into(),
                grading: grading.into(),
                regions: vec![],
            },
        )
    }
    pub fn cut(&mut self, name: &str, normal: (f64, f64, f64), offset: impl Into<Expr>, keep: Keep) -> FeatureId {
        self.push(name, Op::Cut { normal: e3(normal), offset: offset.into(), keep })
    }
    pub fn periodic(
        &mut self,
        name: &str,
        a: Selector,
        b: Selector,
        translation: (impl Into<Expr>, impl Into<Expr>, impl Into<Expr>),
    ) -> FeatureId {
        self.push(name, Op::Periodic { sel_a: a, sel_b: b, translation: e3(translation) })
    }
    pub fn size(
        &mut self,
        name: &str,
        sel: Selector,
        h: Option<Expr>,
        max_deviation: Option<Expr>,
        min_layers: Option<Expr>,
    ) -> FeatureId {
        self.push(name, Op::Size { sel, h, max_deviation, min_layers })
    }
    pub fn tag(&mut self, name: &str, sel: Selector, tag: &str) -> FeatureId {
        self.push(name, Op::Tag { sel, name: tag.to_string() })
    }
}
