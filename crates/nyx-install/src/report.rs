// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Machine-readable JSON report (`--json`, contract `install_report.v1`) emitted by `nyx-setup` and `nyx-uninstall`.

use crate::manifest::InstallScope;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HostOutcome {
    pub id: String,
    pub action: String,
    pub status: String,
    pub config_path: Option<String>,
    pub backup_path: Option<String>,
    pub detail: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DoctorCheck {
    pub name: String,
    pub tier: String,
    pub status: String,
    pub detail: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InstallReport {
    pub schema: String,
    pub version: u32,
    pub ok: bool,
    pub action: String,
    pub dry_run: bool,
    pub scope: InstallScope,
    pub install_root: String,
    pub manifest_path: Option<String>,
    pub nyx_version: Option<String>,
    pub previous_version: Option<String>,
    pub files_written: usize,
    pub files_removed: usize,
    pub path_updated: bool,
    pub desktop_updated: bool,
    pub file_assoc_updated: bool,
    pub ai_hosts: Vec<HostOutcome>,
    pub doctor_checks: Vec<DoctorCheck>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub next_steps: Vec<String>,
}

impl InstallReport {
    pub fn new(
        action: &str,
        scope: InstallScope,
        install_root: impl Into<String>,
        dry_run: bool,
    ) -> Self {
        Self {
            schema: format!("nyx-install.{action}"),
            version: 1,
            ok: true,
            action: action.to_string(),
            dry_run,
            scope,
            install_root: install_root.into(),
            manifest_path: None,
            nyx_version: None,
            previous_version: None,
            files_written: 0,
            files_removed: 0,
            path_updated: false,
            desktop_updated: false,
            file_assoc_updated: false,
            ai_hosts: vec![],
            doctor_checks: vec![],
            warnings: vec![],
            errors: vec![],
            next_steps: vec![],
        }
    }
}

/// JSON Schema (`install_report.v1.schema.json`) describing CLI `--json` reports.
pub fn report_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.dev/contracts/install_report.v1.schema.json",
        "title": "install_report.v1",
        "description": "Machine-readable report printed by `nyx-setup <action> --json`.",
        "type": "object",
        "additionalProperties": false,
        "required": [
            "schema", "version", "ok", "action", "dry_run", "scope", "install_root",
            "manifest_path", "nyx_version", "previous_version", "files_written", "files_removed",
            "path_updated", "desktop_updated", "file_assoc_updated",
            "ai_hosts", "doctor_checks", "warnings", "errors", "next_steps"
        ],
        "properties": {
            "schema": { "enum": ["nyx-install.install", "nyx-install.upgrade", "nyx-install.repair", "nyx-install.uninstall", "nyx-install.status", "nyx-install.doctor"] },
            "version": { "const": 1 },
            "ok": { "type": "boolean" },
            "action": { "enum": ["install", "upgrade", "repair", "uninstall", "status", "doctor"] },
            "dry_run": { "type": "boolean" },
            "scope": { "enum": ["user", "system", "custom"] },
            "install_root": { "type": "string" },
            "manifest_path": { "type": ["string", "null"] },
            "nyx_version": { "type": ["string", "null"] },
            "previous_version": { "type": ["string", "null"] },
            "files_written": { "type": "integer" },
            "files_removed": { "type": "integer" },
            "path_updated": { "type": "boolean" },
            "desktop_updated": { "type": "boolean" },
            "file_assoc_updated": { "type": "boolean" },
            "ai_hosts": { "type": "array", "items": { "$ref": "#/$defs/host_outcome" } },
            "doctor_checks": { "type": "array", "items": { "$ref": "#/$defs/doctor_check" } },
            "warnings": { "type": "array", "items": { "type": "string" } },
            "errors": { "type": "array", "items": { "type": "string" } },
            "next_steps": { "type": "array", "items": { "type": "string" } }
        },
        "$defs": {
            "host_outcome": {
                "type": "object",
                "additionalProperties": false,
                "required": ["id", "action", "status", "config_path", "backup_path", "detail"],
                "properties": {
                    "id": { "type": "string" },
                    "action": { "type": "string" },
                    "status": { "type": "string" },
                    "config_path": { "type": ["string", "null"] },
                    "backup_path": { "type": ["string", "null"] },
                    "detail": { "type": ["string", "null"] }
                }
            },
            "doctor_check": {
                "type": "object",
                "additionalProperties": false,
                "required": ["name", "tier", "status", "detail"],
                "properties": {
                    "name": { "type": "string" },
                    "tier": { "enum": ["required", "optional"] },
                    "status": { "enum": ["ok", "warn", "fail", "info"] },
                    "detail": { "type": "string" }
                }
            }
        }
    })
}
