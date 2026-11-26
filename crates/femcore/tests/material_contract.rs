// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Contract of the material file formats:
//! - Format 1 compatibility exemplars (`material_schema.v1.json`, `material_kitchen_sink.v1.json`,
//!   `material_kitchen_sink_models.v1.json`) MUST remain unchanged.
//! - Format 2 exemplars (`material_schema.v2.json`, `material_kitchen_sink.v2.json`) and agreement of
//!   the JSON Schemas with the loader keys in `spec.rs`. Rewrite v2 exemplars with `NYX_WRITE_CONTRACT=1`
//!   and a `VERSIONS.md` entry.

use femcore::material::schema::{material_schema, material_schema_v1, SCHEMA_VERSION, SCHEMA_VERSION_V1};
use femcore::material::spec::{keys as loader_keys, FORMAT, FORMAT_V1};
use femcore::material::MaterialSpec;
use serde_json::Value;
use std::collections::BTreeSet;
use std::path::PathBuf;

fn contracts() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/architecture/contracts")
}

fn check(name: &str, got: &Value) {
    let path = contracts().join(name);
    if std::env::var("NYX_WRITE_CONTRACT").is_ok() {
        std::fs::write(&path, serde_json::to_string_pretty(got).unwrap() + "\n").unwrap();
    }
    let want: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{name}: {e}"))).unwrap();
    assert_eq!(&want, got, "{name} differs from the code; rewrite it on purpose with NYX_WRITE_CONTRACT=1 and a VERSIONS.md line");
}

fn check_immutable_v1(name: &str, got: &Value) {
    let path = contracts().join(name);
    let want: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{name}: {e}"))).unwrap();
    assert_eq!(&want, got, "v1 contract {name} must stay backwards-compatible!");
}

#[test]
fn schema_exemplars() {
    check_immutable_v1("material_schema.v1.json", &material_schema_v1());
    assert_eq!(material_schema_v1()["x-schema-version"], SCHEMA_VERSION_V1);
    assert_eq!(material_schema_v1()["x-file-format"], FORMAT_V1);

    check("material_schema.v2.json", &material_schema());
    assert_eq!(material_schema()["x-schema-version"], SCHEMA_VERSION);
    assert_eq!(material_schema()["x-file-format"], FORMAT);
}

#[test]
fn kitchen_sink_v1_resolves_to_the_exemplar() {
    let text = std::fs::read_to_string(contracts().join("material_kitchen_sink.v1.toml")).unwrap();
    let m = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}"));
    check_immutable_v1("material_kitchen_sink.v1.json", &serde_json::to_value(&m).unwrap());
}

#[test]
fn model_kitchen_sink_v1_resolves_to_the_exemplar() {
    let text = std::fs::read_to_string(contracts().join("material_kitchen_sink_models.v1.toml")).unwrap();
    let m = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}"));
    check_immutable_v1("material_kitchen_sink_models.v1.json", &serde_json::to_value(&m).unwrap());
}

#[test]
fn kitchen_sink_v2_resolves_to_the_exemplar() {
    let text = std::fs::read_to_string(contracts().join("material_kitchen_sink.v2.toml")).unwrap();
    let m = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}"));
    check("material_kitchen_sink.v2.json", &serde_json::to_value(&m).unwrap());
}

fn keys(v: &Value) -> BTreeSet<String> {
    v["properties"].as_object().unwrap().keys().filter(|k| !is_const_key(k)).cloned().collect()
}

fn is_const_key(k: &str) -> bool {
    let b = k.as_bytes();
    if b.len() == 3
        && matches!(b[0], b'c' | b's' | b'e' | b'd' | b'p' | b'r')
        && (b'1'..=b'6').contains(&b[1])
        && (b'1'..=b'6').contains(&b[2])
    {
        return true;
    }
    b.len() == 5
        && k.starts_with("eta")
        && (b'1'..=b'6').contains(&b[3])
        && (b'1'..=b'6').contains(&b[4])
}

fn set(list: &[&str]) -> BTreeSet<String> {
    list.iter().filter(|k| !is_const_key(k)).map(|k| k.to_string()).collect()
}

#[test]
fn the_schema_lists_exactly_the_keys_the_loader_knows() {
    let k = loader_keys();

    // Verify v1 schema parity
    let s1 = material_schema_v1();
    assert_eq!(keys(&s1), set(k.top_v1), "v1 top level");
    assert_eq!(keys(&s1["properties"]["elastic"]), set(k.elastic), "v1 [elastic]");
    assert_eq!(keys(&s1["properties"]["electrical"]), set(k.electrical), "v1 [electrical]");

    // Verify v2 schema parity across all sections
    let s2 = material_schema();
    assert_eq!(keys(&s2), set(k.top), "v2 top level");
    assert_eq!(keys(&s2["properties"]["elastic"]), set(k.elastic), "v2 [elastic]");
    assert_eq!(keys(&s2["properties"]["thermal"]), set(k.thermal), "v2 [thermal]");
    assert_eq!(keys(&s2["properties"]["electrical"]), set(k.electrical), "v2 [electrical]");
    assert_eq!(keys(&s2["properties"]["dielectric"]), set(k.dielectric), "v2 [dielectric]");
    assert_eq!(keys(&s2["properties"]["superconducting"]), set(k.superconducting), "v2 [superconducting]");
    assert_eq!(keys(&s2["properties"]["piezoelectric"]), set(k.piezoelectric), "v2 [piezoelectric]");
    assert_eq!(keys(&s2["properties"]["photoelastic"]), set(k.photoelastic), "v2 [photoelastic]");
    assert_eq!(keys(&s2["properties"]["electro_optic"]), set(k.electro_optic), "v2 [electro_optic]");
    assert_eq!(keys(&s2["properties"]["viscoelastic"]), set(k.viscoelastic), "v2 [viscoelastic]");
    assert_eq!(keys(&s2["properties"]["interfaces"]), set(k.interfaces), "v2 [interfaces]");
    assert_eq!(keys(&s2["$defs"]["interface_layer"]), set(k.interface_layer), "v2 interface_layer");
    let mut optical = set(k.optical);
    optical.extend(set(k.model));
    assert_eq!(keys(&s2["properties"]["optical"]), optical, "v2 [optical]");
    assert_eq!(keys(&s2["$defs"]["model"]), set(k.model), "v2 model table");
    assert_eq!(keys(&s2["$defs"]["temp_curve"]), set(k.temp_curve), "v2 temp_curve table");
}
