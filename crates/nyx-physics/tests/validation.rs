// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Structural validation of a problem against the descriptors, and the JSON Schemas generated from them.
mod common;
use common::*;
use nyx_physics::schema::{params_schema, physics_schema};
use nyx_physics::*;
use serde_json::json;

fn good() -> ProblemSpec {
    let mut s = ProblemSpec::new("tst", 2, "modes");
    s.study.params = json!({ "freq": 2e14, "n_modes": 3, "kind": "b", "v": [1.0, 2.0, 3.0], "skip": ["air"] });
    s.materials.insert("core".into(), MaterialAssign::Props([("eps_r".to_string(), 12.0)].into()));
    s.materials.insert("clad".into(), MaterialAssign::Named("silica".into()));
    s.bcs.push(BcSpec { on: "outer".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) });
    s.params = json!({ "plane": "stress" });
    s
}

fn err(spec: &ProblemSpec) -> String {
    validate_spec(&Mock::new("tst"), spec).unwrap_err().to_string()
}

#[test]
fn a_well_formed_problem_passes() {
    validate_spec(&Mock::new("tst"), &good()).unwrap();
}

#[test]
fn the_physics_order_study_and_conditions_are_checked() {
    let mut s = good();
    s.physics = "other".into();
    assert!(err(&s).contains("was given to"));
    let mut s = good();
    s.order = 0;
    assert!(err(&s).contains("element order must be"));
    let mut s = good();
    s.study.kind = "cavity".into();
    let e = err(&s);
    assert!(e.contains("no study 'cavity'") && e.contains("modes"), "{e}");
    let mut s = good();
    s.bcs[0].kind = "glue".into();
    assert!(err(&s).contains("no condition 'glue'"));
    let mut s = good();
    s.bcs[0].on = " ".into();
    assert!(err(&s).contains("does not say where"));
}

#[test]
fn parameters_are_checked_for_presence_type_and_range() {
    let mut s = good();
    s.study.params = json!({ "n_modes": 3 });
    assert!(err(&s).contains("'freq'") && err(&s).contains("required"));
    let mut s = good();
    s.study.params = json!({ "freq": 0.5 });
    assert!(err(&s).contains("outside the range"), "{}", err(&s));
    let mut s = good();
    s.study.params = json!({ "freq": "fast" });
    assert!(err(&s).contains("finite number"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "n_modes": 0 });
    assert!(err(&s).contains("outside the range"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "n_modes": 2.5 });
    assert!(err(&s).contains("integer"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "kind": "c" });
    assert!(err(&s).contains("not one of"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "v": [1.0, 2.0] });
    assert!(err(&s).contains("3 finite numbers"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "evanescent": 1 });
    assert!(err(&s).contains("true or false"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "skip": [1] });
    assert!(err(&s).contains("region names"));
    let mut s = good();
    s.study.params = json!({ "freq": 1e9, "frequncy": 2 });
    let e = err(&s);
    assert!(e.contains("unknown parameter 'frequncy'") && e.contains("freq"), "{e}");
    let mut s = good();
    s.study.params = json!([1, 2]);
    assert!(err(&s).contains("must be an object"));
    let mut s = good();
    s.params = json!({ "plane": "strain", "colour": "red" });
    assert!(err(&s).contains("physics options") && err(&s).contains("colour"));
}

#[test]
fn material_property_names_must_be_ones_the_physics_reads() {
    let mut s = good();
    s.materials.insert("core".into(), MaterialAssign::Props([("eps".to_string(), 12.0)].into()));
    let e = err(&s);
    assert!(e.contains("region 'core'") && e.contains("'eps'") && e.contains("eps_r"), "{e}");
    let mut s = good();
    s.materials.insert("core".into(), MaterialAssign::Props([("E".to_string(), f64::NAN)].into()));
    assert!(err(&s).contains("finite"));
}

#[test]
fn the_json_schema_describes_the_problem_the_physics_accepts() {
    let m = Mock::new("tst");
    let schema = physics_schema(&m);
    assert_eq!(schema["properties"]["physics"]["const"], json!("tst"));
    assert_eq!(schema["properties"]["order"]["minimum"], json!(1));
    let study = &schema["properties"]["study"]["oneOf"][0];
    assert_eq!(study["properties"]["kind"]["const"], json!("modes"));
    let freq = &study["properties"]["params"]["properties"]["freq"];
    assert_eq!(freq["type"], json!("number"));
    assert_eq!(freq["minimum"], json!(1.0));
    assert_eq!(study["properties"]["params"]["required"], json!(["freq"]));
    assert_eq!(study["properties"]["params"]["additionalProperties"], json!(false));
    assert_eq!(study["properties"]["params"]["properties"]["kind"]["enum"], json!(["a", "b"]));
    assert_eq!(study["properties"]["params"]["properties"]["v"]["minItems"], json!(3));
    assert_eq!(study["properties"]["params"]["properties"]["skip"]["items"]["type"], json!("string"));
    assert_eq!(schema["properties"]["bcs"]["items"]["oneOf"][0]["properties"]["kind"]["const"], json!("fixed"));
    assert!(schema["properties"]["materials"]["additionalProperties"]["oneOf"][1]["properties"].get("E").is_some());
    assert!(freq["description"].as_str().unwrap().contains("[Hz]"));
    // an empty parameter list is an object that takes no keys
    assert_eq!(params_schema(&[])["additionalProperties"], json!(false));
}

#[test]
fn the_schema_and_the_validator_agree_on_the_example() {
    // every parameter the schema lists is one the validator accepts, and the other way round
    let m = Mock::new("tst");
    let names: Vec<String> = m.studies()[0].params.iter().map(|p| p.name.clone()).collect();
    let schema = params_schema(&m.studies()[0].params);
    let listed: Vec<String> = schema["properties"].as_object().unwrap().keys().cloned().collect();
    let mut a = names.clone();
    let mut b = listed;
    a.sort();
    b.sort();
    assert_eq!(a, b);
}

#[test]
fn optional_null_parameters_use_default_while_required_null_is_refused() {
    let mut s = good();
    // optional parameter `n_modes` is passed explicitly as JSON `null`
    s.study.params = json!({ "freq": 2e14, "n_modes": null });
    validate_spec(&Mock::new("tst"), &s).unwrap();
    // required parameter `freq` passed as `null` must still be rejected
    s.study.params = json!({ "freq": null });
    assert!(err(&s).contains("'freq'") && err(&s).contains("required"));
}
