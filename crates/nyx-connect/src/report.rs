//! The machine-readable form of what `nyx-connect` finds and does (`nyx-connect list|connect|disconnect|test --json`), and its JSON Schema.
//!
//! Three reports, all with `{"schema": "nyx-connect.<kind>", "version": 1, ...}`:
//! * `list`: the server path and every host with its status;
//! * `change`: the outcome of a `connect` or `disconnect`, one row per host (success or the plain-language error);
//! * `test`: the answer of the handshake with the server, or why there was none.
//!
//! The version is the major of the contract `connect_report` (`docs/architecture/contracts/connect_report.v1.*`): an added optional key keeps
//! it, a renamed, removed or re-typed key raises it. [`schema`] is generated here and compared to the checked-in schema file by the contract test.
use crate::config::{Done, HostInfo, Status};
use crate::server::ServerInfo;
use serde_json::{json, Value};
use std::path::Path;

/// The major version written into every report.
pub const VERSION: u64 = 1;

fn head(kind: &str) -> serde_json::Map<String, Value> {
    let mut m = serde_json::Map::new();
    m.insert("schema".into(), json!(format!("nyx-connect.{kind}")));
    m.insert("version".into(), json!(VERSION));
    m
}

/// `{"state": ..., "found"?: ..., "why"?: ...}`; `state` is one of `not_installed`, `not_connected`, `connected`, `stale`, `unreadable`.
pub fn status_json(s: &Status) -> Value {
    match s {
        Status::NotInstalled => json!({ "state": "not_installed" }),
        Status::NotConnected => json!({ "state": "not_connected" }),
        Status::Connected => json!({ "state": "connected" }),
        Status::Stale { found } => json!({ "state": "stale", "found": found }),
        Status::Unreadable { why } => json!({ "state": "unreadable", "why": why }),
    }
}

pub fn host_json(h: &HostInfo) -> Value {
    json!({ "id": h.id, "name": h.name, "config": h.config.to_string_lossy(), "after_connect": h.after_connect, "status": status_json(&h.status) })
}

/// The `list` report. `server_found` says whether `server` is a file.
pub fn list_report(server: &Path, hosts: &[HostInfo]) -> Value {
    let mut m = head("list");
    m.insert("server".into(), json!(server.to_string_lossy()));
    m.insert("server_found".into(), json!(server.is_file()));
    m.insert("hosts".into(), Value::Array(hosts.iter().map(host_json).collect()));
    Value::Object(m)
}

/// One row of a `change` report.
pub fn change_row(host: &str, result: &Result<Done, String>) -> Value {
    match result {
        Ok(d) => json!({ "host": host, "ok": true, "changed": d.changed, "config": d.config.to_string_lossy(), "backup": d.backup.as_ref().map(|b| b.to_string_lossy().into_owned()) }),
        Err(e) => json!({ "host": host, "ok": false, "error": e }),
    }
}

/// The `change` report; `action` is `connect` or `disconnect`.
pub fn change_report(action: &str, rows: Vec<Value>) -> Value {
    let mut m = head("change");
    m.insert("action".into(), json!(action));
    m.insert("results".into(), Value::Array(rows));
    Value::Object(m)
}

/// The `test` report.
pub fn test_report(result: &Result<ServerInfo, String>) -> Value {
    let mut m = head("test");
    match result {
        Ok(i) => {
            m.insert("ok".into(), json!(true));
            m.insert("server".into(), json!({ "name": i.name, "version": i.version, "protocol": i.protocol, "tools": i.tools }));
        }
        Err(e) => {
            m.insert("ok".into(), json!(false));
            m.insert("error".into(), json!(e));
        }
    }
    Value::Object(m)
}

/// The JSON Schema (2020-12) of the three reports: a `oneOf` over `#/$defs/list`, `#/$defs/change` and `#/$defs/test`.
/// Only the keywords `type`, `const`, `enum`, `required`, `properties`, `additionalProperties`, `items`, `oneOf` and `$ref` are used, so a small
/// validator is enough (the contract test has one).
pub fn schema() -> Value {
    let string = || json!({ "type": "string" });
    let nullable_string = || json!({ "type": ["string", "null"] });
    let kind = |k: &str| json!({ "const": format!("nyx-connect.{k}") });
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "nyx-connect.report.v1",
        "title": "nyx-connect report, version 1",
        "oneOf": [{ "$ref": "#/$defs/list" }, { "$ref": "#/$defs/change" }, { "$ref": "#/$defs/test" }],
        "$defs": {
            "status": {
                "oneOf": [
                    { "type": "object", "required": ["state"], "additionalProperties": false,
                      "properties": { "state": { "enum": ["not_installed", "not_connected", "connected"] } } },
                    { "type": "object", "required": ["state", "found"], "additionalProperties": false,
                      "properties": { "state": { "const": "stale" }, "found": string() } },
                    { "type": "object", "required": ["state", "why"], "additionalProperties": false,
                      "properties": { "state": { "const": "unreadable" }, "why": string() } }
                ]
            },
            "host": {
                "type": "object",
                "required": ["id", "name", "config", "after_connect", "status"],
                "additionalProperties": false,
                "properties": {
                    "id": { "enum": ["claude-code", "claude-desktop", "cursor", "vscode", "windsurf", "antigravity"] },
                    "name": string(), "config": string(), "after_connect": string(),
                    "status": { "$ref": "#/$defs/status" }
                }
            },
            "list": {
                "type": "object",
                "required": ["schema", "version", "server", "server_found", "hosts"],
                "additionalProperties": false,
                "properties": {
                    "schema": kind("list"), "version": { "const": VERSION },
                    "server": string(), "server_found": { "type": "boolean" },
                    "hosts": { "type": "array", "items": { "$ref": "#/$defs/host" } }
                }
            },
            "row": {
                "oneOf": [
                    { "type": "object", "required": ["host", "ok", "changed", "config", "backup"], "additionalProperties": false,
                      "properties": { "host": string(), "ok": { "const": true }, "changed": { "type": "boolean" }, "config": string(), "backup": nullable_string() } },
                    { "type": "object", "required": ["host", "ok", "error"], "additionalProperties": false,
                      "properties": { "host": string(), "ok": { "const": false }, "error": string() } }
                ]
            },
            "change": {
                "type": "object",
                "required": ["schema", "version", "action", "results"],
                "additionalProperties": false,
                "properties": {
                    "schema": kind("change"), "version": { "const": VERSION },
                    "action": { "enum": ["connect", "disconnect"] },
                    "results": { "type": "array", "items": { "$ref": "#/$defs/row" } }
                }
            },
            "test": {
                "oneOf": [
                    { "type": "object", "required": ["schema", "version", "ok", "server"], "additionalProperties": false,
                      "properties": {
                          "schema": kind("test"), "version": { "const": VERSION }, "ok": { "const": true },
                          "server": { "type": "object", "required": ["name", "version", "protocol", "tools"], "additionalProperties": false,
                                      "properties": { "name": string(), "version": string(), "protocol": string(), "tools": { "type": "array", "items": string() } } }
                      } },
                    { "type": "object", "required": ["schema", "version", "ok", "error"], "additionalProperties": false,
                      "properties": { "schema": kind("test"), "version": { "const": VERSION }, "ok": { "const": false }, "error": string() } }
                ]
            }
        }
    })
}
