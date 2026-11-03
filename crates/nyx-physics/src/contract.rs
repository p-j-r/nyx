//! Helpers for the contract tests of every physics crate (ADR 19): the snapshot of a physics' descriptors, and the exemplar file check.
//!
//! A contract exemplar is a JSON file in `docs/architecture/contracts/`; the test builds the same JSON from the code and compares the two as
//! values (the order of the keys of an object is not part of a contract). `NYX_WRITE_CONTRACT=1` rewrites the file; do it only together with
//! a line in `contracts/VERSIONS.md`.

use crate::physics::Physics;
use crate::schema::physics_schema;
use crate::spec::canonical_json;
use serde_json::{json, Value};
use std::path::PathBuf;

/// Everything a physics tells the rest of nyx, as one JSON value: the exemplar `physics_<id>.v1.json` is this value.
pub fn descriptor_snapshot(p: &dyn Physics) -> Value {
    let studies = p.studies();
    let fields: serde_json::Map<String, Value> =
        studies.iter().map(|s| (s.kind.clone(), serde_json::to_value(p.fields(&s.kind)).unwrap())).collect();
    json!({
        "id": p.id(),
        "name": p.name(),
        "studies": serde_json::to_value(&studies).unwrap(),
        "bc_kinds": serde_json::to_value(p.bc_kinds()).unwrap(),
        "properties": serde_json::to_value(p.properties()).unwrap(),
        "params": serde_json::to_value(p.params()).unwrap(),
        "fields": Value::Object(fields),
        "exports": serde_json::to_value(p.exports()).unwrap(),
        "imports": serde_json::to_value(p.imports()).unwrap(),
        "schema": physics_schema(p),
    })
}

pub fn contract_dir(manifest_dir: &str) -> PathBuf {
    PathBuf::from(manifest_dir).join("../../docs/architecture/contracts")
}

/// Compares `got` with the exemplar file `name` (call with `env!("CARGO_MANIFEST_DIR")`). Panics with the instruction to version the change.
pub fn check_exemplar(manifest_dir: &str, name: &str, got: &Value) {
    let got = canonical_json(got);
    let path = contract_dir(manifest_dir).join(name);
    if std::env::var("NYX_WRITE_CONTRACT").is_ok_and(|v| v == "1" || name.starts_with(&v)) {
        std::fs::write(&path, serde_json::to_string_pretty(&got).unwrap() + "\n").unwrap();
        return;
    }
    let text = std::fs::read_to_string(&path).unwrap_or_else(|_| {
        panic!("missing exemplar {path:?}; write it with NYX_WRITE_CONTRACT=1 and add a VERSIONS.md line")
    });
    let want: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(canonical_json(&want), got, "the contract exemplar {name} changed; if on purpose, rewrite it with NYX_WRITE_CONTRACT=1 and add a VERSIONS.md line (ADR 19)");
}
