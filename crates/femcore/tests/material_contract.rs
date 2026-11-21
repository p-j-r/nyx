// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Contract of the material file format (ADR 30): the JSON Schema exemplar, the kitchen-sink file and its resolved form, and the
//! agreement of the schema with the keys the loader accepts. Rewrite an exemplar with `NYX_WRITE_CONTRACT=1` and a line in
//! `docs/architecture/contracts/VERSIONS.md`.

use femcore::material::schema::{material_schema, SCHEMA_VERSION};
use femcore::material::spec::FORMAT;
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

#[test]
fn schema_exemplar() {
    check("material_schema.v1.json", &material_schema());
    assert_eq!(material_schema()["x-schema-version"], SCHEMA_VERSION);
    assert_eq!(material_schema()["x-file-format"], FORMAT);
}

#[test]
fn kitchen_sink_resolves_to_the_exemplar() {
    let text = std::fs::read_to_string(contracts().join("material_kitchen_sink.v1.toml")).unwrap();
    let m = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}"));
    check("material_kitchen_sink.v1.json", &serde_json::to_value(&m).unwrap());
}

fn keys(v: &Value) -> BTreeSet<String> {
    v["properties"].as_object().unwrap().keys().filter(|k| !is_const_key(k)).cloned().collect()
}

fn is_const_key(k: &str) -> bool {
    let b = k.as_bytes();
    b.len() == 3 && (b[0] == b'c' || b[0] == b's') && (b'1'..=b'6').contains(&b[1]) && (b'1'..=b'6').contains(&b[2])
}

fn set(list: &[&str]) -> BTreeSet<String> {
    list.iter().filter(|k| !is_const_key(k)).map(|k| k.to_string()).collect()
}

#[test]
fn the_schema_lists_exactly_the_keys_the_loader_knows() {
    use femcore::material::spec::{keys as loader_keys};
    let s = material_schema();
    let k = loader_keys();
    assert_eq!(keys(&s), set(k.top), "top level");
    assert_eq!(keys(&s["properties"]["elastic"]), set(k.elastic), "[elastic]");
    assert_eq!(keys(&s["properties"]["thermal"]), set(k.thermal), "[thermal]");
    assert_eq!(keys(&s["properties"]["electrical"]), set(k.electrical), "[electrical]");
    let mut optical = set(k.optical);
    optical.extend(set(k.model));
    assert_eq!(keys(&s["properties"]["optical"]), optical, "[optical]");
    assert_eq!(keys(&s["properties"]["superconducting"]), set(k.superconducting), "[superconducting]");
    assert_eq!(keys(&s["$defs"]["model"]), set(k.model), "model table");
}

#[test]
fn model_kitchen_sink_resolves_to_the_exemplar() {
    let text = std::fs::read_to_string(contracts().join("material_kitchen_sink_models.v1.toml")).unwrap();
    let m = MaterialSpec::from_toml(&text).unwrap_or_else(|e| panic!("{e}"));
    check("material_kitchen_sink_models.v1.json", &serde_json::to_value(&m).unwrap());
}
