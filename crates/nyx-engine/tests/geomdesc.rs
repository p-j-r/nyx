// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Contract of `geomdesc` (ADR 38): every operation has a descriptor, the descriptor's fields are exactly the serde fields of the
//! operation, `build_op` accepts the friendly form (names, selector text) and round-trips through `args_of`, and refusals name what is wrong.
use geom::model::{FeatureList, Model, Op};
use nyx_engine::geomdesc::*;
use serde_json::{json, Value};
use std::collections::BTreeSet;

fn sample(t: Ty) -> Value {
    match t {
        Ty::Expr => json!(1.0),
        Ty::Exprs(n) => json!(vec![1.0; n]),
        Ty::Points(_) => json!([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]),
        Ty::Rings(_) => json!([[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]]),
        Ty::Axis => json!("Z"),
        Ty::Feature => json!("p"),
        Ty::Features => json!(["p"]),
        Ty::Text => json!("t"),
        Ty::TextList => json!([]),
        Ty::Selector => json!("outer"),
        Ty::Sides => json!([]),
        Ty::Keep => json!("Above"),
        Ty::Naming => json!("Feature"),
        Ty::Regions => json!([{ "name": "a", "input": "p" }]),
        Ty::Pins => json!([{ "sel": "outer", "id": 1 }]),
        Ty::Steps => json!([{ "Translate": [1.0, 0.0, 0.0] }]),
        Ty::Pattern => json!({ "Linear": { "step": [1.0, 0.0, 0.0] } }),
        Ty::PmlMode => json!("Outward"),
        Ty::PmlCoord => json!("Cartesian"),
    }
}

fn base() -> FeatureList {
    let mut m = Model::new_2d();
    m.rect("p", (0.0, 0.0), (1.0, 1.0));
    m.list
}

fn required_args(d: &OpDesc) -> Value {
    Value::Object(d.fields.iter().filter(|f| f.required).map(|f| (f.name.to_string(), sample(f.ty))).collect())
}

fn all_args(d: &OpDesc) -> Value {
    Value::Object(d.fields.iter().map(|f| (f.name.to_string(), sample(f.ty))).collect())
}

#[test]
fn every_kind_is_described_once_with_its_serde_fields() {
    let list = base();
    let all = all();
    let kinds: BTreeSet<_> = all.iter().map(|d| d.kind).collect();
    assert_eq!(kinds.len(), all.len(), "a kind is described twice");
    assert_eq!(all.len(), 25, "the number of Op variants");
    for d in &all {
        let op = build_op(&list, d.kind, &required_args(d)).unwrap_or_else(|e| panic!("{}: {e}", d.kind));
        assert_eq!(kind_of(&op), d.kind);
        let op_all = build_op(&list, d.kind, &all_args(d)).unwrap_or_else(|e| panic!("{}: {e}", d.kind));
        let v = serde_json::to_value(&op_all).unwrap();
        let got: BTreeSet<String> = v[d.kind].as_object().unwrap().keys().cloned().collect();
        let want: BTreeSet<String> = d.fields.iter().map(|f| f.name.to_string()).collect();
        assert_eq!(got, want, "{}: the descriptor fields differ from the serde fields", d.kind);
        // the schema lists the required fields and no other
        let s = schema(d);
        let req: BTreeSet<String> =
            s["required"].as_array().unwrap().iter().map(|x| x.as_str().unwrap().to_string()).collect();
        let want_req: BTreeSet<String> = d.fields.iter().filter(|f| f.required).map(|f| f.name.to_string()).collect();
        assert_eq!(req, want_req);
    }
}

#[test]
fn args_of_gives_back_what_build_op_took() {
    let list = base();
    for d in all() {
        let args = required_args(&d);
        let op = build_op(&list, d.kind, &args).unwrap();
        let back = args_of(&list, &op);
        let op2 = build_op(&list, d.kind, &back).unwrap_or_else(|e| panic!("{}: {e}", d.kind));
        assert_eq!(op, op2, "{}", d.kind);
    }
}

#[test]
fn defaults_fill_optional_fields_and_a_null_means_default() {
    let list = base();
    let op = build_op(&list, "Extrude", &json!({ "profile": "p", "from": 0.0, "to": 1.0 })).unwrap();
    match op {
        Op::Extrude { axis, taper_deg, .. } => {
            assert_eq!(axis, geom::model::Axis::Y);
            assert!(taper_deg.is_none());
        }
        _ => panic!(),
    }
    build_op(&list, "Ellipse", &json!({ "center": [0, 0], "radii": [1, 2], "angle_deg": null })).unwrap();
}

#[test]
fn refusals_say_what_is_wrong_in_plain_words() {
    let list = base();
    let e = build_op(&list, "Rect", &json!({ "min": [0, 0] })).unwrap_err().to_string();
    assert!(e.contains("needs `max`"), "{e}");
    let e = build_op(&list, "Rect", &json!({ "min": [0, 0], "max": [1, 1], "radius": 2 })).unwrap_err().to_string();
    assert!(e.contains("no field 'radius'") && e.contains("min, max"), "{e}");
    let e = build_op(&list, "Union", &json!({ "ids": ["nope"] })).unwrap_err().to_string();
    assert!(e.contains("no feature named 'nope'") && e.contains("p"), "{e}");
    let e = build_op(&list, "Tag", &json!({ "sel": "plane(", "name": "x" })).unwrap_err().to_string();
    assert!(e.starts_with("Tag.sel"), "{e}");
    let e = build_op(&list, "Hexagon", &json!({})).unwrap_err().to_string();
    assert!(e.contains("Rect"), "{e}");
    let e = build_op(&list, "Rect", &json!({ "min": [0, 0], "max": [1] })).unwrap_err().to_string();
    assert!(e.starts_with("Rect:"), "{e}");
}

#[test]
fn tool_names_are_snake_case() {
    assert_eq!(tool_name("Rect"), "add_rect");
    assert_eq!(tool_name("RoundedRect"), "add_rounded_rect");
    let names: BTreeSet<String> = all().iter().map(|d| tool_name(d.kind)).collect();
    assert_eq!(names.len(), all().len());
}

#[test]
fn a_strip_a_rib_and_a_ring_can_be_rebuilt_from_descriptors_alone() {
    // the acceptance test of gate G5: the strip, rib and ring are built only through build_op and evaluate
    let mut list = FeatureList::new(2);
    let add = |list: &mut FeatureList, name: &str, kind: &str, args: Value| {
        let op = build_op(list, kind, &args).unwrap();
        let id = list.features.iter().map(|f| f.id).max().map_or(1, |m| m + 1);
        list.features.push(geom::model::Feature {
            id,
            name: name.into(),
            op,
            suppressed: false,
            origin: geom::model::Origin::Code,
            tags: vec![],
        });
    };
    add(&mut list, "box", "Rect", json!({ "min": [-3e-6, -3e-6], "max": [3e-6, 0.0] }));
    add(&mut list, "slab", "Rect", json!({ "min": [-3e-6, 0.0], "max": [3e-6, 1e-7] }));
    add(&mut list, "rib", "Rect", json!({ "min": [-2.5e-7, 1e-7], "max": [2.5e-7, 3e-7] }));
    add(&mut list, "ring_out", "Circle", json!({ "center": [0.0, 0.0], "r": 5e-6 }));
    add(&mut list, "ring_in", "Circle", json!({ "center": [0.0, 0.0], "r": 4.5e-6 }));
    add(&mut list, "ring", "Difference", json!({ "a": "ring_out", "b": ["ring_in"] }));
    add(
        &mut list,
        "regions",
        "Overlay",
        json!({ "regions": [
        { "name": "BOX", "input": "box" }, { "name": "slab", "input": "slab" }, { "name": "rib", "input": "rib" } ] }),
    );
    list.validate().unwrap();
    geom::eval::Evaluator::new().eval(&list, geom::eval::Tier::Complex).unwrap();
}

/// The exemplar `docs/architecture/contracts/geom_operations.v1.json` pins the descriptors (ADR 38, ADR 19). Rewrite it with
/// `NYX_WRITE_CONTRACT=1` together with a line in `contracts/VERSIONS.md`.
#[test]
fn the_descriptors_match_their_contract_exemplar() {
    fn sorted(v: &Value) -> Value {
        match v {
            Value::Object(m) => {
                let mut keys: Vec<_> = m.keys().collect();
                keys.sort();
                Value::Object(keys.into_iter().map(|k| (k.clone(), sorted(&m[k]))).collect())
            }
            Value::Array(a) => Value::Array(a.iter().map(sorted).collect()),
            o => o.clone(),
        }
    }
    let snapshot = Value::Array(
        all()
            .iter()
            .map(|d| {
                json!({
                    "kind": d.kind, "tool": tool_name(d.kind), "title": d.title, "doc": d.doc, "group": d.group.name(), "dims": d.dims,
                    "fields": d.fields.iter().map(|f| json!({ "name": f.name, "type": format!("{:?}", f.ty), "required": f.required, "default": f.default, "doc": f.doc })).collect::<Vec<_>>(),
                    "schema": schema(d),
                })
            })
            .collect(),
    );
    let got = sorted(&snapshot);
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/architecture/contracts/geom_operations.v1.json");
    if std::env::var("NYX_WRITE_CONTRACT").is_ok_and(|v| v == "1" || v == "geom_operations") {
        std::fs::write(&path, serde_json::to_string_pretty(&got).unwrap() + "\n").unwrap();
        return;
    }
    let want: Value = serde_json::from_str(
        &std::fs::read_to_string(&path).expect("exemplar missing; write it with NYX_WRITE_CONTRACT=1"),
    )
    .unwrap();
    assert_eq!(
        sorted(&want),
        got,
        "the descriptors changed; if on purpose, rewrite the exemplar and add a VERSIONS.md line (ADR 19)"
    );
}
