// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! JSON Schemas generated from descriptors, so the MCP server, the window and the documentation read the same truth. The schemas are
//! draft 2020-12 documents; `physics_schema` describes the whole `ProblemSpec` of one physics.

use crate::desc::{ParamDesc, ParamKind, PropKind};
use crate::physics::Physics;
use serde_json::{json, Map, Value};

fn param_schema(d: &ParamDesc) -> Value {
    let mut s = Map::new();
    match &d.kind {
        ParamKind::Number { default, min, max } => {
            s.insert("type".into(), json!("number"));
            s.insert("default".into(), json!(default));
            if let Some(m) = min {
                s.insert("minimum".into(), json!(m));
            }
            if let Some(m) = max {
                s.insert("maximum".into(), json!(m));
            }
        }
        ParamKind::Integer { default, min, max } => {
            s.insert("type".into(), json!("integer"));
            s.insert("default".into(), json!(default));
            if let Some(m) = min {
                s.insert("minimum".into(), json!(m));
            }
            if let Some(m) = max {
                s.insert("maximum".into(), json!(m));
            }
        }
        ParamKind::Bool { default } => {
            s.insert("type".into(), json!("boolean"));
            s.insert("default".into(), json!(default));
        }
        ParamKind::Choice { options, default } => {
            s.insert("type".into(), json!("string"));
            s.insert("enum".into(), json!(options));
            s.insert("default".into(), json!(default));
        }
        ParamKind::Vector { len, default } => {
            s.insert("type".into(), json!("array"));
            s.insert("items".into(), json!({ "type": "number" }));
            s.insert("minItems".into(), json!(len));
            s.insert("maxItems".into(), json!(len));
            s.insert("default".into(), json!(default));
        }
        ParamKind::Text { default } => {
            s.insert("type".into(), json!("string"));
            s.insert("default".into(), json!(default));
        }
        ParamKind::Names { default } => {
            s.insert("type".into(), json!("array"));
            s.insert("items".into(), json!({ "type": "string" }));
            s.insert("default".into(), json!(default));
        }
    }
    let mut desc = d.label.clone();
    if !d.unit.is_empty() {
        desc.push_str(&format!(" [{}]", d.unit));
    }
    if !d.doc.is_empty() {
        desc.push_str(&format!(". {}", d.doc));
    }
    s.insert("description".into(), json!(desc));
    Value::Object(s)
}

/// The schema of an object made of `descs`: unknown keys are refused, required ones listed.
pub fn params_schema(descs: &[ParamDesc]) -> Value {
    let mut props = Map::new();
    for d in descs {
        props.insert(d.name.clone(), param_schema(d));
    }
    let required: Vec<&str> = descs.iter().filter(|d| d.required).map(|d| d.name.as_str()).collect();
    json!({ "type": "object", "properties": props, "required": required, "additionalProperties": false })
}

/// The schema of a `ProblemSpec` for `p`: its study kinds, condition kinds, property names and physics options.
pub fn physics_schema(p: &dyn Physics) -> Value {
    let studies: Vec<Value> = p.studies().iter().map(|s| json!({ "type": "object", "properties": { "kind": { "const": s.kind }, "params": params_schema(&s.params) }, "required": ["kind"], "additionalProperties": false, "description": s.doc })).collect();
    let bcs: Vec<Value> = p.bc_kinds().iter().map(|k| json!({ "type": "object", "properties": { "on": { "type": "string", "description": "a region name or a boundary selector" }, "kind": { "const": k.kind }, "params": params_schema(&k.params) }, "required": ["on", "kind"], "additionalProperties": false, "description": k.doc })).collect();
    let mut props = Map::new();
    for d in p.properties() {
        let desc = format!("{} [{}]", d.label, d.unit);
        let prop_schema = match d.kind {
            PropKind::Tensor => json!({
                "description": desc,
                "oneOf": [
                    { "type": "number" },
                    { "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3 },
                    { "type": "array", "items": { "type": "number" }, "minItems": 6, "maxItems": 6 },
                    {
                        "type": "array",
                        "items": { "type": "array", "items": { "type": "number" }, "minItems": 3, "maxItems": 3 },
                        "minItems": 3,
                        "maxItems": 3
                    }
                ]
            }),
            _ => json!({ "type": "number", "description": desc }),
        };
        props.insert(d.name.clone(), prop_schema);
    }
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "title": format!("nyx problem for {} ({})", p.name(), p.id()),
        "type": "object",
        "properties": {
            "physics": { "const": p.id() },
            "order": { "type": "integer", "minimum": 1, "maximum": 9 },
            "study": { "oneOf": studies },
            "materials": { "type": "object", "additionalProperties": { "oneOf": [ { "type": "string", "description": "a material by name" }, { "type": "object", "properties": props, "additionalProperties": false } ] }, "description": "material by region name" },
            "bcs": { "type": "array", "items": { "oneOf": bcs } },
            "params": params_schema(&p.params()),
        },
        "required": ["physics", "order", "study"],
        "additionalProperties": false
    })
}
