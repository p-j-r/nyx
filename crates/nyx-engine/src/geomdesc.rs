// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Descriptors of the geom feature operations (they live here, not in `geom`, because `geom` depends on serde only: ADR 3): for every `Op` its kind, what it does, its fields with their types and defaults, and the JSON
//! Schema built from them. The MCP server generates one tool per descriptor (restructure gate G5: every tool geom offers, none hand-picked),
//! the window can build its forms from the same table, and `build_op` turns the friendly arguments of a tool call (selectors as text,
//! features by name) into the `Op` the feature list stores.
//!
//! Completeness is enforced twice: `kind_of` matches every `Op` variant, so a new variant does not compile until it is listed here, and
//! the tests check that the field names of every descriptor equal the serde field names of a sample operation of that kind.

use geom::model::{FeatureList, Op};
use geom::select::Selector;
use nyx_physics::Error;
use serde_json::{json, Map, Value};

type GeomResult<T> = Result<T, Error>;

fn inv(text: String) -> Error {
    Error::Invalid(text)
}

/// How the value of one field is written in a tool call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ty {
    /// a number or an expression text over the model parameters (`"w/2"`)
    Expr,
    /// a fixed number of expressions (a point, a size)
    Exprs(usize),
    /// a list of points with this many expressions each
    Points(usize),
    /// a list of closed rings, each a list of points with this many expressions
    Rings(usize),
    /// `"X"`, `"Y"` or `"Z"`
    Axis,
    /// another feature: its name (or its id)
    Feature,
    /// a list of features
    Features,
    Text,
    TextList,
    /// a selector in its text form (`outer & plane(x = 0)`)
    Selector,
    /// a list of `xmin xmax ymin ymax zmin zmax`
    Sides,
    /// `"Above"` or `"Below"`
    Keep,
    /// `"Feature"`, `"Region"` or `"Numbered"`
    Naming,
    /// regions in painter order: `{ name, input (feature), kind, pinned_id }`
    Regions,
    /// label pins: `{ sel (selector text), id }`
    Pins,
    /// transformation steps in the serde form of `TStep` (`{"Translate": [dx, dy, dz]}`, `{"Rotate": {...}}`)
    Steps,
    /// the serde form of `PatternKind` (`{"Linear": {"step": [..]}}` or `{"Polar": {...}}`)
    Pattern,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FieldDesc {
    pub name: &'static str,
    pub ty: Ty,
    pub doc: &'static str,
    /// must be given; a field that is not required has `default`
    pub required: bool,
    pub default: Value,
}

/// What an operation makes, for grouping tools.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Group {
    /// a 2D or 3D primitive
    Shape,
    /// moves, copies, sweeps or combines other features
    Build,
    /// names the regions and the boundaries and sets the mesh sizes
    Labels,
    /// absorbing layers, cuts, periodicity
    Boundary,
}

impl Group {
    pub fn name(self) -> &'static str {
        match self {
            Group::Shape => "shape",
            Group::Build => "build",
            Group::Labels => "naming",
            Group::Boundary => "boundary",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct OpDesc {
    /// the `Op` variant name (`"Rect"`)
    pub kind: &'static str,
    pub title: &'static str,
    pub doc: &'static str,
    pub group: Group,
    /// the dimensions of feature list it can be added to
    pub dims: &'static [u8],
    pub fields: Vec<FieldDesc>,
}

fn req(name: &'static str, ty: Ty, doc: &'static str) -> FieldDesc {
    FieldDesc { name, ty, doc, required: true, default: Value::Null }
}

fn opt(name: &'static str, ty: Ty, doc: &'static str, default: Value) -> FieldDesc {
    FieldDesc { name, ty, doc, required: false, default }
}

const SWEEP_AXIS_DOC: &str = "the axis of the sweep";

/// Every operation, in a fixed order (shapes, builders, naming, boundary).
pub fn all() -> Vec<OpDesc> {
    use Group::{Boundary, Build, Labels, Shape};
    use Ty::{
        Axis, Expr, Exprs, Feature, Features, Keep, Naming, Pattern, Pins, Points, Regions, Rings, Selector, Sides,
        Steps, Text, TextList,
    };
    let d = |kind, title, doc, group, dims, fields| OpDesc { kind, title, doc, group, dims, fields };
    vec![
        d("Rect", "Rectangle", "A rectangle in the x-z cross-section, given by its two corners.", Shape, &[2], vec![
            req("min", Exprs(2), "lower-left corner (x, z) in metres"),
            req("max", Exprs(2), "upper-right corner (x, z) in metres"),
        ]),
        d("Circle", "Circle", "A disc in the x-z cross-section.", Shape, &[2], vec![
            req("center", Exprs(2), "centre (x, z)"),
            req("r", Expr, "radius"),
        ]),
        d("Ellipse", "Ellipse", "An ellipse in the x-z cross-section, rotated by an angle.", Shape, &[2], vec![
            req("center", Exprs(2), "centre (x, z)"),
            req("radii", Exprs(2), "semi-axes along the rotated x and z"),
            opt("angle_deg", Expr, "rotation in degrees", json!(0.0)),
        ]),
        d("Polygon", "Polygon", "A closed polygon in the x-z cross-section; counter-clockwise or clockwise. Edge i goes from point i to point i+1.", Shape, &[2], vec![
            req("pts", Points(2), "the corners (x, z)"),
            opt("names", TextList, "role names of the edges (default e0, e1, ...), used in boundary names", json!([])),
        ]),
        d("Polygons", "Polygon batch", "A batch of outer polygon rings and optional hole rings in 2D or 3D, with Ramer-Douglas-Peucker closed-ring simplification (shared corners pinned) and minimum hole area filtering.", Shape, &[2, 3], vec![
            req("polygons", Rings(2), "outer polygon rings, each a list of 2D corners in metres"),
            opt("holes", Rings(2), "hole polygon rings to subtract", json!([])),
            opt("simplify_tol", Expr, "Ramer-Douglas-Peucker simplification tolerance in metres", json!(50e-9)),
            opt("min_hole_area", Expr, "drop hole rings with unsigned shoelace area below this threshold in m^2", json!(0.0)),
        ]),
        d("RoundedRect", "Rounded rectangle", "A rectangle with quarter-circle corners.", Shape, &[2], vec![
            req("min", Exprs(2), "lower-left corner (x, z)"),
            req("max", Exprs(2), "upper-right corner (x, z)"),
            req("r", Expr, "corner radius"),
        ]),
        d("Box", "Box", "An axis-aligned box (3D). Face roles: xmin xmax ymin ymax zmin zmax.", Shape, &[3], vec![
            req("min", Exprs(3), "lower corner (x, y, z)"),
            req("max", Exprs(3), "upper corner (x, y, z)"),
        ]),
        d("Cylinder", "Cylinder", "A cylinder (3D) from `base` along +z. Face roles: side bottom top.", Shape, &[3], vec![
            req("base", Exprs(3), "centre of the bottom face"),
            opt("axis", Axis, "this build: z only", json!("Z")),
            req("r", Expr, "radius"),
            req("h", Expr, "height"),
        ]),
        d("Sphere", "Sphere", "A sphere (3D). One face, `surface`; it must not touch any other face.", Shape, &[3], vec![
            req("center", Exprs(3), "centre"),
            req("r", Expr, "radius"),
        ]),
        d("Cone", "Cone frustum", "A cone frustum (3D) along the axis from `base`; either radius may be 0. Faces side, bottom, top.", Shape, &[3], vec![
            req("base", Exprs(3), "centre of the bottom face"),
            opt("axis", Axis, "the axis of the cone", json!("Z")),
            req("r0", Expr, "radius at the bottom"),
            req("r1", Expr, "radius at the top"),
            req("h", Expr, "height"),
        ]),
        d("Revolve", "Revolve", "A full turn of a Rect or Polygon profile about the axis through `center` (3D). The profile's (x, y) are the radius (>= 0) and the height.", Build, &[3], vec![
            req("profile", Feature, "the 2D profile feature"),
            opt("axis", Axis, "the axis of revolution", json!("Z")),
            req("center", Exprs(3), "a point on the axis"),
        ]),
        d("Torus", "Torus", "A torus with axis z (3D). One face, `surface`.", Shape, &[3], vec![
            req("center", Exprs(3), "centre"),
            req("major", Expr, "major radius"),
            req("minor", Expr, "minor radius"),
        ]),
        d("Extrude", "Extrude", "Sweeps a 2D profile feature along an axis (default y, the propagation direction) into a 3D solid. Roles side:<profile role> and the caps ymin, ymax.", Build, &[3], vec![
            req("profile", Feature, "the 2D profile feature"),
            opt("axis", Axis, SWEEP_AXIS_DOC, json!("Y")),
            req("from", Expr, "start coordinate along the axis"),
            req("to", Expr, "end coordinate along the axis"),
            opt("taper_deg", Expr, "wall angle in degrees (positive narrows towards `to`); makes the solid polyhedral", Value::Null),
        ]),
        d("Sheet", "2D sheet in 3D", "Embeds a 2D profile or boolean feature as a zero-thickness conductive sheet at plane z (3D).", Build, &[3], vec![
            req("input", Feature, "the 2D profile or boolean feature"),
            opt("z", Expr, "z coordinate of the sheet plane in metres", json!(0.0)),
        ]),
        d("Transform", "Transform", "Moves, rotates, mirrors or scales a feature; the steps apply in order.", Build, &[2, 3], vec![
            req("input", Feature, "the feature to transform"),
            req("steps", Steps, "list of steps in the serde form of TStep"),
        ]),
        d("Pattern", "Pattern", "n copies of a feature, in a line or around an axis.", Build, &[2, 3], vec![
            req("input", Feature, "the feature to repeat"),
            req("kind", Pattern, "Linear { step } or Polar { axis, center, total_deg }"),
            req("n", Expr, "number of copies"),
        ]),
        d("Union", "Union", "The union of features.", Build, &[2, 3], vec![req("ids", Features, "the features to join")]),
        d("Difference", "Difference", "The feature `a` with the features `b` cut away.", Build, &[2, 3], vec![
            req("a", Feature, "the feature kept"),
            req("b", Features, "the features removed"),
        ]),
        d("Intersect", "Intersection", "The common part of features.", Build, &[2, 3], vec![req("ids", Features, "the features to intersect")]),
        d("Overlay", "Regions (overlay)", "The root of the model: names the regions, in painter order (a later region wins where they overlap). A list has one overlay. Pins fix the integer label of selected boundaries.", Labels, &[2, 3], vec![
            req("regions", Regions, "regions in painter order: { name, input (feature name), kind (Material | Hole | Pml | Void), pinned_id }"),
            opt("pins", Pins, "label pins: { sel (selector text), id }", json!([])),
            opt("naming", Naming, "how boundary faces against the void are named", json!("Feature")),
        ]),
        d("Pml", "Absorbing layers", "Absorbing layers of `thickness` around a window on some sides; splits every matching region into the physical part and one piece per side combination.", Boundary, &[2, 3], vec![
            req("window", Feature, "the feature that bounds the physical window"),
            opt("sides", Sides, "sides that get a layer (empty: all)", json!([])),
            req("thickness", Expr, "layer thickness"),
            req("layers", Expr, "number of element layers across it"),
            req("grading", Expr, "order of the polynomial profile"),
            opt("regions", TextList, "regions to split (empty: all)", json!([])),
        ]),
        d("Cut", "Cut by a plane", "Keeps one side of the plane normal . p = offset. The cut face is named <feature>.cut and tagged `symmetry`.", Boundary, &[2, 3], vec![
            req("normal", Exprs(3), "plane normal"),
            req("offset", Expr, "plane offset along the normal"),
            req("keep", Keep, "which side to keep"),
        ]),
        d("Periodic", "Periodic pair", "Declares two boundary sets as periodic images of each other.", Boundary, &[2, 3], vec![
            req("sel_a", Selector, "first boundary set"),
            req("sel_b", Selector, "second boundary set"),
            req("translation", Exprs(3), "the translation from A to B"),
        ]),
        d("Size", "Mesh size", "Sets the element size near the selected boundaries or regions.", Labels, &[2, 3], vec![
            req("sel", Selector, "where the size applies"),
            opt("h", Expr, "element size", Value::Null),
            opt("max_deviation", Expr, "largest distance of a chord from its curve", Value::Null),
            opt("min_layers", Expr, "fewest element layers across a thin region", Value::Null),
        ]),
        d("Tag", "Tag", "Gives a selection a name (a boundary group) to refer to later.", Labels, &[2, 3], vec![
            req("sel", Selector, "what to tag"),
            req("name", Text, "the tag"),
        ]),
    ]
}

/// The kind name of an operation. Exhaustive on purpose: a new `Op` variant does not compile until it is added here and in `all`.
pub fn kind_of(op: &Op) -> &'static str {
    match op {
        Op::Rect { .. } => "Rect",
        Op::Circle { .. } => "Circle",
        Op::Ellipse { .. } => "Ellipse",
        Op::Polygon { .. } => "Polygon",
        Op::Polygons { .. } => "Polygons",
        Op::RoundedRect { .. } => "RoundedRect",
        Op::Box { .. } => "Box",
        Op::Cylinder { .. } => "Cylinder",
        Op::Sphere { .. } => "Sphere",
        Op::Cone { .. } => "Cone",
        Op::Revolve { .. } => "Revolve",
        Op::Torus { .. } => "Torus",
        Op::Extrude { .. } => "Extrude",
        Op::Sheet { .. } => "Sheet",
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

pub fn find(kind: &str) -> Option<OpDesc> {
    all().into_iter().find(|d| d.kind == kind)
}

/// The tool name for an operation: `add_` and the snake_case kind (`add_rounded_rect`).
pub fn tool_name(kind: &str) -> String {
    let mut s = String::from("add_");
    for (i, c) in kind.chars().enumerate() {
        if c.is_ascii_uppercase() {
            if i > 0 {
                s.push('_');
            }
            s.push(c.to_ascii_lowercase());
        } else {
            s.push(c);
        }
    }
    s
}

// ------------------------------------------------------------------------------------------------------------- JSON Schema

fn expr_schema() -> Value {
    json!({ "oneOf": [{ "type": "number" }, { "type": "string" }], "description": "a number in metres (or the unit of the field) or an expression over the model parameters" })
}

fn ty_schema(t: Ty) -> Value {
    match t {
        Ty::Expr => expr_schema(),
        Ty::Exprs(n) => json!({ "type": "array", "items": expr_schema(), "minItems": n, "maxItems": n }),
        Ty::Points(n) => {
            json!({ "type": "array", "items": { "type": "array", "items": expr_schema(), "minItems": n, "maxItems": n }, "minItems": 3 })
        }
        Ty::Rings(n) => {
            json!({ "type": "array", "items": { "type": "array", "items": { "type": "array", "items": expr_schema(), "minItems": n, "maxItems": n }, "minItems": 3 } })
        }
        Ty::Axis => json!({ "enum": ["X", "Y", "Z"] }),
        Ty::Feature => {
            json!({ "oneOf": [{ "type": "string" }, { "type": "integer" }], "description": "the name of an earlier feature (or its id)" })
        }
        Ty::Features => json!({ "type": "array", "items": { "oneOf": [{ "type": "string" }, { "type": "integer" }] } }),
        Ty::Text => json!({ "type": "string" }),
        Ty::TextList => json!({ "type": "array", "items": { "type": "string" } }),
        Ty::Selector => json!({ "type": "string", "description": "selector text, e.g. `outer & plane(x = 0)`" }),
        Ty::Sides => json!({ "type": "array", "items": { "enum": ["XMin", "XMax", "YMin", "YMax", "ZMin", "ZMax"] } }),
        Ty::Keep => json!({ "enum": ["Above", "Below"] }),
        Ty::Naming => json!({ "enum": ["Feature", "Region", "Numbered"] }),
        Ty::Regions => json!({
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": { "type": "string" },
                    "input": { "oneOf": [{ "type": "string" }, { "type": "integer" }] },
                    "kind": { "enum": ["Material", "Hole", "Pml", "Void"], "default": "Material" },
                    "pinned_id": { "type": ["integer", "null"] }
                },
                "required": ["name", "input"],
                "additionalProperties": false
            }
        }),
        Ty::Pins => json!({
            "type": "array",
            "items": {
                "type": "object",
                "properties": { "sel": { "type": "string" }, "id": { "type": "integer" } },
                "required": ["sel", "id"],
                "additionalProperties": false
            }
        }),
        Ty::Steps => json!({
            "type": "array",
            "items": { "type": "object", "minProperties": 1, "maxProperties": 1,
                "description": "{\"Translate\": [dx, dy, dz]} | {\"Rotate\": {axis, deg, about}} | {\"Mirror\": {axis, at}} | {\"Scale\": {factor, about}}" }
        }),
        Ty::Pattern => json!({ "type": "object", "minProperties": 1, "maxProperties": 1,
            "description": "{\"Linear\": {\"step\": [dx, dy, dz]}} | {\"Polar\": {axis, center, total_deg}}" }),
    }
}

/// The JSON Schema of the arguments of one operation (the `inputSchema` of its tool, without the fields every add tool shares).
pub fn schema(d: &OpDesc) -> Value {
    let mut props = Map::new();
    let mut required = vec![];
    for f in &d.fields {
        let mut s = ty_schema(f.ty);
        if let Value::Object(o) = &mut s {
            o.insert("title".into(), json!(f.doc));
            if f.required && matches!(f.ty, Ty::Rings(_)) {
                o.insert("minItems".into(), json!(1));
            }
            if !f.required {
                o.insert("default".into(), f.default.clone());
            }
        }
        props.insert(f.name.into(), s);
        if f.required {
            required.push(json!(f.name));
        }
    }
    json!({ "type": "object", "properties": Value::Object(props), "required": required, "additionalProperties": false })
}

// ------------------------------------------------------------------------------------------------------------- building an Op

fn feature_ref(list: &FeatureList, v: &Value, what: &str) -> GeomResult<Value> {
    match v {
        Value::Number(n) if n.is_u64() => {
            let id = n.as_u64().unwrap() as u32;
            if list.feature(id).is_none() {
                return Err(inv(format!("{what}: there is no feature with id {id}")));
            }
            Ok(json!(id))
        }
        Value::String(s) => match list.feature_by_name(s) {
            Some(f) => Ok(json!(f.id)),
            None => Err(inv(format!(
                "{what}: there is no feature named '{s}' (the features are: {})",
                list.features.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(", ")
            ))),
        },
        other => Err(inv(format!("{what}: a feature is given by its name or id, not {other}"))),
    }
}

fn selector_value(v: &Value, what: &str) -> GeomResult<Value> {
    match v {
        Value::String(s) => {
            let sel = Selector::parse(s).map_err(|e| inv(format!("{what}: {e}")))?;
            serde_json::to_value(sel).map_err(|e| inv(e.to_string()))
        }
        other if other.is_object() => Ok(other.clone()),
        other => Err(inv(format!("{what}: a selector is written as text, e.g. \"outer & plane(x = 0)\", not {other}"))),
    }
}

fn convert(list: &FeatureList, f: &FieldDesc, v: &Value, kind: &str) -> GeomResult<Value> {
    let what = format!("{kind}.{}", f.name);
    let items = |v: &Value| -> GeomResult<Vec<Value>> {
        v.as_array().cloned().ok_or_else(|| inv(format!("{what} must be a list")))
    };
    Ok(match f.ty {
        Ty::Feature => feature_ref(list, v, &what)?,
        Ty::Features => Value::Array(items(v)?.iter().map(|x| feature_ref(list, x, &what)).collect::<GeomResult<_>>()?),
        Ty::Selector => selector_value(v, &what)?,
        Ty::Regions => {
            let mut out = vec![];
            for r in items(v)? {
                let o = r
                    .as_object()
                    .ok_or_else(|| inv(format!("{what}: each region is an object {{name, input, kind}}")))?;
                let mut o = o.clone();
                let input = o
                    .get("input")
                    .ok_or_else(|| inv(format!("{what}: a region needs `input` (the feature it is painted from)")))?;
                let id = feature_ref(list, input, &what)?;
                o.insert("input".into(), id);
                o.entry("kind").or_insert(json!("Material"));
                o.entry("pinned_id").or_insert(Value::Null);
                out.push(Value::Object(o));
            }
            Value::Array(out)
        }
        Ty::Pins => {
            let mut out = vec![];
            for p in items(v)? {
                let mut o =
                    p.as_object().cloned().ok_or_else(|| inv(format!("{what}: each pin is an object {{sel, id}}")))?;
                let sel = o.get("sel").ok_or_else(|| inv(format!("{what}: a pin needs `sel`")))?;
                let sel = selector_value(sel, &what)?;
                o.insert("sel".into(), sel);
                out.push(Value::Object(o));
            }
            Value::Array(out)
        }
        Ty::Rings(_) => {
            let arr = items(v)?;
            // Accept both a list of rings `[[[x, y], ...], ...]` and a single ring `[[x, y], ...]`.
            let is_single_ring = arr
                .first()
                .and_then(Value::as_array)
                .and_then(|pt| pt.first())
                .is_some_and(|coord| coord.is_number() || coord.is_string());
            if is_single_ring {
                Value::Array(vec![Value::Array(arr)])
            } else {
                Value::Array(arr)
            }
        }
        _ => v.clone(),
    })
}

/// Builds the operation `kind` from the arguments of a tool call. Features are named by name (or id), selectors are text; optional
/// fields take their defaults; unknown or missing fields are reported with the names the operation has. The result is checked by serde,
/// which is the same reader that loads a saved feature list.
pub fn build_op(list: &FeatureList, kind: &str, args: &Value) -> GeomResult<Op> {
    let d = find(kind).ok_or_else(|| {
        inv(format!(
            "there is no operation '{kind}' (they are: {})",
            all().iter().map(|d| d.kind).collect::<Vec<_>>().join(", ")
        ))
    })?;
    let empty = Map::new();
    let given = match args {
        Value::Null => &empty,
        Value::Object(m) => m,
        other => return Err(inv(format!("the arguments of {kind} must be an object, not {other}"))),
    };
    for k in given.keys() {
        if !d.fields.iter().any(|f| f.name == k) {
            return Err(inv(format!(
                "{kind} has no field '{k}' (it has: {})",
                d.fields.iter().map(|f| f.name).collect::<Vec<_>>().join(", ")
            )));
        }
    }
    let mut out = Map::new();
    for f in &d.fields {
        match given.get(f.name) {
            Some(v) if !(v.is_null() && !f.required) => {
                out.insert(f.name.into(), convert(list, f, v, kind)?);
            }
            _ if f.required => {
                return Err(inv(format!("{kind} needs `{}`: {}", f.name, f.doc)));
            }
            _ => {
                out.insert(f.name.into(), f.default.clone());
            }
        }
    }
    let mut tagged = Map::new();
    tagged.insert(kind.into(), Value::Object(out));
    serde_json::from_value::<Op>(Value::Object(tagged)).map_err(|e| inv(format!("{kind}: {e}")))
}

/// The arguments of an existing operation in the form `build_op` takes (features by name, selectors as text): what an editing tool shows
/// and what `update_feature` merges into.
pub fn args_of(list: &FeatureList, op: &Op) -> Value {
    let kind = kind_of(op);
    let d = find(kind).expect("every kind has a descriptor");
    let v = serde_json::to_value(op).expect("an Op serialises");
    let mut fields = v.get(kind).and_then(Value::as_object).cloned().unwrap_or_default();
    let name_of = |id: &Value| -> Value {
        id.as_u64().and_then(|i| list.feature(i as u32)).map(|f| json!(f.name)).unwrap_or_else(|| id.clone())
    };
    for f in &d.fields {
        let Some(val) = fields.get(f.name).cloned() else { continue };
        let nv = match f.ty {
            Ty::Feature => name_of(&val),
            Ty::Features => Value::Array(val.as_array().into_iter().flatten().map(&name_of).collect()),
            Ty::Selector => {
                serde_json::from_value::<Selector>(val.clone()).map(|s| json!(s.to_string())).unwrap_or(val)
            }
            Ty::Regions => Value::Array(
                val.as_array()
                    .into_iter()
                    .flatten()
                    .map(|r| {
                        let mut r = r.clone();
                        if let Some(i) = r.get("input").cloned() {
                            r["input"] = name_of(&i);
                        }
                        r
                    })
                    .collect(),
            ),
            Ty::Pins => Value::Array(
                val.as_array()
                    .into_iter()
                    .flatten()
                    .map(|p| {
                        let mut p = p.clone();
                        if let Some(s) = p.get("sel").cloned().and_then(|s| serde_json::from_value::<Selector>(s).ok())
                        {
                            p["sel"] = json!(s.to_string());
                        }
                        p
                    })
                    .collect(),
            ),
            _ => val,
        };
        fields.insert(f.name.into(), nv);
    }
    Value::Object(fields)
}
