// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 VVUQ verification tests for `femcore::material`:
//! - `T_range` alignment across `TempCurve`, `spec.rs`, `schema.rs`, and library TOMLs (`additionalProperties: false`);
//! - `TempCurve::eval_checked`, `TempCurve::eval_strict`, and `TempTensor::eval_*_checked` non-silent clamping;
//! - `AxisModel` / `Optical::nk_checked` / `Optical::permittivity_global_checked` `dn_dT_range` clamping metadata;
//! - `Source::flags_uncertainty()`, `is_approximate()`, and `is_unverified()`.

use femcore::material::schema::material_schema;
use femcore::material::state_t::PointGroup;
use femcore::material::{Library, MaterialSpec, Source, TempClampInfo, TempCurve, TempTensor, TEMP_CURVE_KEYS};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::path::Path;

fn repo_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().parent().unwrap()
}

fn toml_to_json(v: &toml::Value) -> Value {
    match v {
        toml::Value::String(s) => Value::String(s.clone()),
        toml::Value::Integer(i) => json!(*i),
        toml::Value::Float(f) => json!(*f),
        toml::Value::Boolean(b) => Value::Bool(*b),
        toml::Value::Datetime(d) => Value::String(d.to_string()),
        toml::Value::Array(a) => Value::Array(a.iter().map(toml_to_json).collect()),
        toml::Value::Table(t) => Value::Object(t.iter().map(|(k, x)| (k.clone(), toml_to_json(x))).collect()),
    }
}

fn resolve_ref<'a>(root: &'a Value, schema: &'a Value) -> &'a Value {
    if let Some(r) = schema.get("$ref").and_then(Value::as_str) {
        if let Some(def) = r.strip_prefix("#/$defs/") {
            return &root["$defs"][def];
        }
    }
    schema
}

fn json_kind(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "boolean",
        Value::Number(n) if n.is_i64() || n.is_u64() => "integer",
        Value::Number(_) => "number",
        Value::String(_) => "string",
        Value::Array(_) => "array",
        Value::Object(_) => "object",
    }
}

fn validate_json_schema(root: &Value, schema: &Value, v: &Value, at: &str) -> Result<(), String> {
    let schema = resolve_ref(root, schema);
    if let Some(c) = schema.get("const") {
        if c != v {
            return Err(format!("{at}: expected const {c}, got {v}"));
        }
    }
    if let Some(e) = schema.get("enum") {
        if !e.as_array().unwrap().contains(v) {
            return Err(format!("{at}: {v} not in enum {e}"));
        }
    }
    if let Some(t) = schema.get("type") {
        let allowed: Vec<&str> = match t {
            Value::String(s) => vec![s.as_str()],
            Value::Array(a) => a.iter().map(|x| x.as_str().unwrap()).collect(),
            _ => vec![],
        };
        let k = json_kind(v);
        if !allowed.iter().any(|a| *a == k || (*a == "number" && k == "integer")) {
            return Err(format!("{at}: type {k}, expected {allowed:?}"));
        }
    }
    if let (Some(x), Some(m)) = (v.as_f64(), schema.get("minimum").and_then(Value::as_f64)) {
        if x < m {
            return Err(format!("{at}: {x} < minimum {m}"));
        }
    }
    if let (Some(x), Some(m)) = (v.as_f64(), schema.get("exclusiveMinimum").and_then(Value::as_f64)) {
        if x <= m {
            return Err(format!("{at}: {x} <= exclusiveMinimum {m}"));
        }
    }
    if let Some(one) = schema.get("oneOf") {
        let mut errs = Vec::new();
        let mut ok = 0usize;
        for alt in one.as_array().unwrap() {
            match validate_json_schema(root, alt, v, at) {
                Ok(()) => ok += 1,
                Err(e) => errs.push(e),
            }
        }
        if ok != 1 {
            return Err(format!("{at}: matched {ok} of oneOf alternatives ({})", errs.join(" | ")));
        }
    }
    if let Some(any) = schema.get("anyOf") {
        let mut errs = Vec::new();
        let mut ok = 0usize;
        for alt in any.as_array().unwrap() {
            match validate_json_schema(root, alt, v, at) {
                Ok(()) => ok += 1,
                Err(e) => errs.push(e),
            }
        }
        if ok == 0 {
            return Err(format!("{at}: matched 0 of anyOf alternatives ({})", errs.join(" | ")));
        }
    }
    if let Value::Object(obj) = v {
        for r in schema.get("required").and_then(Value::as_array).into_iter().flatten() {
            let k = r.as_str().unwrap();
            if !obj.contains_key(k) {
                return Err(format!("{at}: missing required key '{k}'"));
            }
        }
        let props = schema.get("properties").and_then(Value::as_object);
        let forbid_extra = schema.get("additionalProperties") == Some(&Value::Bool(false));
        for (k, val) in obj {
            if let Some(ps) = props.and_then(|p| p.get(k)) {
                validate_json_schema(root, ps, val, &format!("{at}.{k}"))?;
            } else if forbid_extra {
                return Err(format!("{at}: unexpected property '{k}' under additionalProperties=false"));
            }
        }
    }
    if let Value::Array(arr) = v {
        if let Some(min_i) = schema.get("minItems").and_then(Value::as_u64) {
            if (arr.len() as u64) < min_i {
                return Err(format!("{at}: len {} < minItems {min_i}", arr.len()));
            }
        }
        if let Some(max_i) = schema.get("maxItems").and_then(Value::as_u64) {
            if (arr.len() as u64) > max_i {
                return Err(format!("{at}: len {} > maxItems {max_i}", arr.len()));
            }
        }
        if let Some(items) = schema.get("items") {
            if items.is_object() {
                for (i, elem) in arr.iter().enumerate() {
                    validate_json_schema(root, items, elem, &format!("{at}[{i}]"))?;
                }
            }
        }
    }
    Ok(())
}

#[test]
fn t_range_schema_closes_temp_curve_and_model_and_validates_toml_materials() {
    let schema = material_schema();
    assert_eq!(schema["$defs"]["temp_curve"]["additionalProperties"], json!(false));
    assert_eq!(schema["$defs"]["model"]["additionalProperties"], json!(false));

    let tc_props: BTreeSet<String> = schema["$defs"]["temp_curve"]["properties"]
        .as_object()
        .unwrap()
        .keys()
        .cloned()
        .collect();
    let expected_tc: BTreeSet<String> = TEMP_CURVE_KEYS.iter().map(|s| (*s).to_string()).collect();
    assert_eq!(tc_props, expected_tc, "TEMP_CURVE_KEYS must match $defs.temp_curve");
    assert!(tc_props.contains("T_range"));
    assert!(!tc_props.contains("t_range"));

    let mut samples = vec![repo_root().join("docs/architecture/contracts/material_kitchen_sink.v2.toml")];
    let mut lib_files: Vec<_> = std::fs::read_dir(repo_root().join("materials"))
        .unwrap()
        .filter_map(|e| e.ok().map(|x| x.path()))
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("toml"))
        .collect();
    lib_files.sort();
    samples.extend(lib_files);
    for path in samples {
        let text = std::fs::read_to_string(&path).unwrap();
        let raw: toml::Table = text.parse().unwrap();
        // Format-1 TOML library files use the v1 schema; validate all format = 2 TOMLs against material_schema().
        if raw.get("format").and_then(|v| v.as_integer()) != Some(2) {
            continue;
        }
        let json_val = toml_to_json(&toml::Value::Table(raw));
        validate_json_schema(&schema, &schema, &json_val, path.file_name().unwrap().to_str().unwrap())
            .unwrap_or_else(|e| panic!("{} failed schema validation: {e}", path.display()));
    }

    // Lowercase `t_range` inside log10_poly must fail BOTH JSON Schema (additionalProperties: false) AND MaterialSpec loader.
    let bad_toml = r#"
format = 2
name = "bad_lowercase_t_range"
[thermal]
specific_heat_T = { form = "log10_poly", t_range = [4.0, 300.0], coeffs = [1.0, 0.2], unit = "J/(kg*K)" }
"#;
    let bad_json = toml_to_json(&toml::Value::Table(bad_toml.parse::<toml::Table>().unwrap()));
    let schema_err = validate_json_schema(&schema, &schema, &bad_json, "bad").unwrap_err();
    assert!(schema_err.contains("additionalProperties") || schema_err.contains("oneOf"), "{schema_err}");
    let load_err = MaterialSpec::from_toml(bad_toml).unwrap_err().to_string();
    assert!(load_err.contains("T_range"), "{load_err}");
}

#[test]
fn temp_curve_eval_checked_and_eval_strict_report_out_of_range_clamps() {
    let dir = repo_root().join("materials");
    let lib = Library::load_dir(&dir).unwrap();
    let cu4k = lib.require("copper_4k").unwrap();

    let cp_curve = cu4k.thermal.specific_heat_t.as_ref().expect("copper_4k specific_heat_T");
    assert_eq!(cp_curve.t_range(), Some([4.0, 300.0]));

    // Endpoints and interior are NOT clamped.
    for t_ok in [4.0, 77.0, 300.0] {
        let checked = cp_curve.eval_checked(t_ok);
        assert_eq!(checked.clamp, None, "T={t_ok} K should not clamp");
        assert_eq!(checked.value, cp_curve.eval(t_ok));
        assert_eq!(cp_curve.eval_strict(t_ok).unwrap(), checked.value);
    }

    // Below T_range clamps to 4.0 K and triggers eval_strict Err.
    let below = cp_curve.eval_checked(1.5);
    assert_eq!(
        below.clamp,
        Some(TempClampInfo { requested_k: 1.5, clamped_k: 4.0, bounds_k: [4.0, 300.0] })
    );
    assert_eq!(below.value, cp_curve.eval(4.0));
    assert_eq!(cp_curve.eval_strict(1.5).unwrap_err(), below.clamp.unwrap());

    // Above T_range clamps to 300.0 K.
    let above = cp_curve.eval_checked(450.0);
    assert_eq!(
        above.clamp,
        Some(TempClampInfo { requested_k: 450.0, clamped_k: 300.0, bounds_k: [4.0, 300.0] })
    );
    assert_eq!(above.value, cp_curve.eval(300.0));

    // Tabulated TempTensor (`conductivity_T` [4.0, 293.15]) also reports clamping accurately.
    let (k_ok, k_clamps_ok) = cu4k.thermal_conductivity_at_checked(Some(4.0)).unwrap();
    assert!(k_clamps_ok.is_empty());
    assert!((k_ok.as_scalar().unwrap() - 642.3).abs() < 1e-9);

    let (k_hi, k_clamps_hi) = cu4k.thermal_conductivity_at_checked(Some(350.0)).unwrap();
    assert_eq!(
        k_clamps_hi,
        vec![(None, TempClampInfo { requested_k: 350.0, clamped_k: 293.15, bounds_k: [4.0, 293.15] })]
    );
    assert!((k_hi.as_scalar().unwrap() - 401.0).abs() < 1e-9);

    // Uniaxial TempTensor reports distinct perp/par clamping.
    let uniaxial = TempTensor::Uniaxial {
        perp: TempCurve::Tabulated {
            t: vec![10.0, 100.0],
            values: vec![1.0, 2.0],
            t_ref: None,
            source: Vec::new(),
        },
        par: TempCurve::Log10Poly {
            t_range: [20.0, 80.0],
            coeffs: vec![0.5, 0.1],
            t_ref: None,
            source: Vec::new(),
        },
    };
    let (_, u_clamps) = uniaxial.eval_tensor_checked(5.0);
    assert_eq!(u_clamps.len(), 2);
    assert_eq!(u_clamps[0].0, Some("perp"));
    assert_eq!(u_clamps[0].1.clamped_k, 10.0);
    assert_eq!(u_clamps[1].0, Some("par"));
    assert_eq!(u_clamps[1].1.clamped_k, 20.0);
}

#[test]
fn optical_axis_model_dn_dt_range_clamping_and_source_uncertainty_flags() {
    let mat = MaterialSpec::from_toml(
        r#"
format = 2
name = "bounded_thermo_optic"
source = { cite = "Verified unit test", verified = "2026-11-26" }

[optical]
wavelength_unit = "um"
reference_temperature = "293.15 K"
model = "constant"
n = 2.0
dn_dT = "1.0e-4 1/K"
dn_dT_range = [100.0, 400.0]
source = { cite = "APPROXIMATE thermo-optic range exemplar" }
"#,
    )
    .unwrap();
    let opt = mat.optical.as_ref().unwrap();
    let in_range = opt.nk_checked(0, 1.55e-6, Some(200.0)).unwrap();
    assert_eq!(in_range.clamp, None);

    let cold = opt.nk_checked(0, 1.55e-6, Some(40.0)).unwrap();
    assert_eq!(
        cold.clamp,
        Some(TempClampInfo { requested_k: 40.0, clamped_k: 100.0, bounds_k: [100.0, 400.0] })
    );
    assert_eq!(cold.value.0, opt.nk(0, 1.55e-6, Some(100.0)).unwrap().0);

    let (_, eps_clamps) = mat.optical_permittivity_checked(1.55e-6, Some(520.0), [0.0; 3]).unwrap();
    assert_eq!(eps_clamps.len(), 3);
    assert_eq!(eps_clamps[0].1.clamped_k, 400.0);

    // Check Source helpers.
    let s_approx = Source::text("APPROXIMATE room-T estimate");
    let s_unver = Source::text("UNVERIFIED handbook table");
    let s_both = Source::text("APPROXIMATE and UNVERIFIED placeholder");
    let s_verified = Source { cite: "Peer reviewed paper".into(), url: None, doi: None, verified: Some("2026-11-26".into()) };

    assert!(s_approx.flags_uncertainty() && s_approx.is_approximate() && !s_approx.is_unverified());
    assert!(s_unver.flags_uncertainty() && s_unver.is_unverified() && !s_unver.is_approximate());
    assert!(s_both.flags_uncertainty() && s_both.is_approximate() && s_both.is_unverified());
    assert!(!s_verified.flags_uncertainty() && s_verified.is_verified());
    let _ = PointGroup::Isotropic;
}
