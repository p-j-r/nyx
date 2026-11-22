// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pinned interface contracts for `nyx-install` (ADR 19, ADR 122):
//! * `install_manifest.v1`: `install-manifest.json` receipt + JSON Schema.
//! * `install_report.v1`: machine-readable `nyx-setup <action> --json` report + JSON Schema.
//!
//! Rewrite intentionally with:
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-install --test contracts`

use nyx_install::manifest::{
    manifest_schema, AiHostRecord, DesktopEntry, FileAssociation, FileBase, FileRole,
    InstallLayout, InstallManifest, InstallScope, ManifestFile, ManifestSymlink, PathEdit,
    TargetOs, MANIFEST_FORMAT, MANIFEST_SCHEMA,
};
use nyx_install::report::{report_schema, DoctorCheck, HostOutcome, InstallReport};
use serde_json::{json, Value};
use std::path::PathBuf;

fn contract_file(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/architecture/contracts")
        .join(name)
}

fn check_exemplar(name: &str, got: &Value) {
    let path = contract_file(name);
    if std::env::var_os("NYX_WRITE_CONTRACT").is_some() {
        std::fs::write(&path, serde_json::to_string_pretty(got).unwrap() + "\n").unwrap();
        return;
    }
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let want: Value = serde_json::from_str(&text).unwrap();
    assert_eq!(
        &want, got,
        "{name} drifted from pinned contract exemplar (update VERSIONS.md and run with NYX_WRITE_CONTRACT=1)"
    );
}

fn validate(root: &Value, schema: &Value, inst: &Value, at: &str, errs: &mut Vec<String>) {
    const KNOWN: [&str; 14] = [
        "$schema",
        "$id",
        "title",
        "$defs",
        "$ref",
        "oneOf",
        "const",
        "enum",
        "type",
        "required",
        "properties",
        "additionalProperties",
        "items",
        "description",
    ];
    let obj = schema.as_object().expect("schema must be object");
    for k in obj.keys() {
        assert!(
            KNOWN.contains(&k.as_str()),
            "unsupported JSON Schema keyword {k}"
        );
    }
    if let Some(r) = obj.get("$ref") {
        let name = r
            .as_str()
            .unwrap()
            .strip_prefix("#/$defs/")
            .expect("only #/$defs refs");
        return validate(root, &root["$defs"][name], inst, at, errs);
    }
    if let Some(alts) = obj.get("oneOf") {
        let n = alts
            .as_array()
            .unwrap()
            .iter()
            .filter(|a| {
                let mut e = vec![];
                validate(root, a, inst, at, &mut e);
                e.is_empty()
            })
            .count();
        if n != 1 {
            errs.push(format!(
                "{at}: matched {n} of oneOf alternatives (expected 1)"
            ));
        }
    }
    if let Some(c) = obj.get("const") {
        if c != inst {
            errs.push(format!("{at}: expected const {c}, got {inst}"));
        }
    }
    if let Some(e) = obj.get("enum") {
        if !e.as_array().unwrap().contains(inst) {
            errs.push(format!("{at}: {inst} not in enum {e}"));
        }
    }
    if let Some(t) = obj.get("type") {
        let types: Vec<&str> = match t {
            Value::String(s) => vec![s.as_str()],
            Value::Array(a) => a.iter().map(|v| v.as_str().unwrap()).collect(),
            _ => panic!("bad type"),
        };
        let ok = types.iter().any(|t| match *t {
            "object" => inst.is_object(),
            "array" => inst.is_array(),
            "string" => inst.is_string(),
            "boolean" => inst.is_boolean(),
            "null" => inst.is_null(),
            "integer" => inst.is_i64() || inst.is_u64(),
            "number" => inst.is_number(),
            other => panic!("unknown schema type {other}"),
        });
        if !ok {
            errs.push(format!("{at}: expected {types:?}, got {inst}"));
        }
    }
    if let Some(m) = inst.as_object() {
        for k in obj
            .get("required")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            if !m.contains_key(k.as_str().unwrap()) {
                errs.push(format!("{at}: missing key {k}"));
            }
        }
        let props = obj.get("properties").and_then(Value::as_object);
        for (k, v) in m {
            match props.and_then(|p| p.get(k)) {
                Some(s) => validate(root, s, v, &format!("{at}/{k}"), errs),
                None if obj.get("additionalProperties") == Some(&json!(false)) => {
                    errs.push(format!("{at}: unexpected key {k}"));
                }
                None => {}
            }
        }
    }
    if let (Some(items), Some(arr)) = (obj.get("items"), inst.as_array()) {
        for (i, v) in arr.iter().enumerate() {
            validate(root, items, v, &format!("{at}/{i}"), errs);
        }
    }
}

fn manifest_examples() -> Vec<(&'static str, Value)> {
    let win = InstallManifest {
        format: MANIFEST_FORMAT,
        schema: MANIFEST_SCHEMA.into(),
        nyx_version: "0.1.0".into(),
        installer_version: "0.1.0".into(),
        install_id: "019a0000001000100".into(),
        target_triple: "x86_64-pc-windows-msvc".into(),
        os: TargetOs::Windows,
        scope: InstallScope::User,
        layout: InstallLayout::Split,
        install_root: r"C:\Users\alice\AppData\Local\Programs\nyx".into(),
        bin_dir: r"C:\Users\alice\AppData\Local\Programs\nyx\bin".into(),
        resource_dir: r"C:\Users\alice\AppData\Local\Programs\nyx".into(),
        manifest_path: r"C:\Users\alice\AppData\Local\Programs\nyx\install-manifest.json".into(),
        installed_at: "2026-11-22T18:00:00Z".into(),
        updated_at: "2026-11-22T18:00:00Z".into(),
        previous_version: None,
        files: vec![
            ManifestFile {
                rel_path: "nyx.exe".into(),
                base: FileBase::Bin,
                role: FileRole::Binary,
                size: 72_967_680,
                blake3: "b3:0000000000000000000000000000000000000000000000000000000000000001"
                    .into(),
                mode: 0o755,
            },
            ManifestFile {
                rel_path: "THIRD_PARTY_NOTICES.txt".into(),
                base: FileBase::Resource,
                role: FileRole::License,
                size: 4096,
                blake3: "b3:0000000000000000000000000000000000000000000000000000000000000002"
                    .into(),
                mode: 0o644,
            },
        ],
        symlinks: vec![],
        path_edits: vec![PathEdit::WindowsRegistry {
            hive: "HKCU".into(),
            subkey: "Environment".into(),
            value_name: "Path".into(),
            entry_added: r"C:\Users\alice\AppData\Local\Programs\nyx\bin".into(),
            created_value: false,
            preserved_expand_sz: true,
        }],
        desktop_entries: vec![
            DesktopEntry::WindowsShortcut {
                path:
                    r"C:\Users\alice\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\nyx.lnk"
                        .into(),
            },
            DesktopEntry::WindowsUninstallRegistry {
                hive: "HKCU".into(),
                subkey: r"Software\Microsoft\Windows\CurrentVersion\Uninstall\nyx".into(),
            },
        ],
        file_associations: vec![FileAssociation::WindowsRegistryKey {
            hive: "HKCU".into(),
            subkey: r"Software\Classes\nyx.Archive.1".into(),
        }],
        ai_hosts: vec![AiHostRecord {
            host_id: "antigravity".into(),
            config_path: r"C:\Users\alice\.gemini\config\mcp_config.json".into(),
            mcp_exe_path: r"C:\Users\alice\AppData\Local\Programs\nyx\bin\nyx-mcp.exe".into(),
            connected_by_installer: true,
            previous_state: "not_connected".into(),
            backup_path: Some(r"C:\Users\alice\.gemini\config\mcp_config.json.nyx-backup".into()),
        }],
    };
    let linux = InstallManifest {
        format: MANIFEST_FORMAT,
        schema: MANIFEST_SCHEMA.into(),
        nyx_version: "0.1.0".into(),
        installer_version: "0.1.0".into(),
        install_id: "019a0000002000100".into(),
        target_triple: "x86_64-unknown-linux-gnu".into(),
        os: TargetOs::Linux,
        scope: InstallScope::User,
        layout: InstallLayout::Split,
        install_root: "/home/alice/.local/share/nyx".into(),
        bin_dir: "/home/alice/.local/share/nyx/bin".into(),
        resource_dir: "/home/alice/.local/share/nyx".into(),
        manifest_path: "/home/alice/.local/share/nyx/install-manifest.json".into(),
        installed_at: "2026-11-22T18:00:00Z".into(),
        updated_at: "2026-11-22T18:05:00Z".into(),
        previous_version: Some("0.0.9".into()),
        files: vec![],
        symlinks: vec![ManifestSymlink {
            link_path: "/home/alice/.local/bin/nyx".into(),
            target_path: "/home/alice/.local/share/nyx/bin/nyx".into(),
        }],
        path_edits: vec![PathEdit::UnixShellProfile {
            profile_path: "/home/alice/.bashrc".into(),
            marker_block_id: "nyx installer".into(),
            created_file: false,
        }],
        desktop_entries: vec![DesktopEntry::FreedesktopDesktop {
            path: "/home/alice/.local/share/applications/nyx.desktop".into(),
        }],
        file_associations: vec![FileAssociation::FreedesktopMimeFile {
            path: "/home/alice/.local/share/mime/packages/nyx.xml".into(),
        }],
        ai_hosts: vec![],
    };
    vec![
        ("windows_user", serde_json::to_value(win).unwrap()),
        ("linux_user", serde_json::to_value(linux).unwrap()),
    ]
}

fn report_examples() -> Vec<(&'static str, Value)> {
    let mut inst = InstallReport::new(
        "install",
        InstallScope::User,
        r"C:\Users\alice\AppData\Local\Programs\nyx",
        false,
    );
    inst.manifest_path =
        Some(r"C:\Users\alice\AppData\Local\Programs\nyx\install-manifest.json".into());
    inst.nyx_version = Some("0.1.0".into());
    inst.files_written = 50;
    inst.path_updated = true;
    inst.desktop_updated = true;
    inst.file_assoc_updated = true;
    inst.ai_hosts.push(HostOutcome {
        id: "cursor".into(),
        action: "connected".into(),
        status: "connected".into(),
        config_path: Some(r"C:\Users\alice\.cursor\mcp.json".into()),
        backup_path: Some(r"C:\Users\alice\.cursor\mcp.json.nyx-backup".into()),
        detail: Some("Restart Cursor to use nyx.".into()),
    });
    inst.doctor_checks.push(DoctorCheck {
        name: "mcp_handshake".into(),
        tier: "required".into(),
        status: "ok".into(),
        detail: "nyx-mcp 0.1.0 (4 tools)".into(),
    });

    let mut uninst = InstallReport::new(
        "uninstall",
        InstallScope::User,
        "/home/alice/.local/share/nyx",
        false,
    );
    uninst.nyx_version = Some("0.1.0".into());
    uninst.files_removed = 50;
    uninst.path_updated = true;
    uninst.next_steps.push(
        "Your workflows and simulation results in /home/alice/Documents/nyx were preserved.".into(),
    );

    vec![
        ("install_ok", serde_json::to_value(inst).unwrap()),
        ("uninstall_ok", serde_json::to_value(uninst).unwrap()),
    ]
}

#[test]
fn install_manifest_v1_schema_and_exemplars_match_contract() {
    let schema = manifest_schema();
    check_exemplar("install_manifest.v1.schema.json", &schema);
    let ex = manifest_examples();
    let doc: serde_json::Map<String, Value> = ex
        .iter()
        .map(|(k, v)| ((*k).to_string(), v.clone()))
        .collect();
    check_exemplar("install_manifest.v1.json", &Value::Object(doc));
    for (name, inst) in ex {
        let mut errs = vec![];
        validate(&schema, &schema, &inst, "", &mut errs);
        assert!(errs.is_empty(), "{name}: {}", errs.join("\n"));
    }
    // Validator must reject corrupted schema input.
    let mut bad = manifest_examples()[0].1.clone();
    bad["format"] = json!(99);
    let mut errs = vec![];
    validate(&schema, &schema, &bad, "", &mut errs);
    assert!(
        !errs.is_empty(),
        "validator must reject wrong format version"
    );
}

#[test]
fn install_report_v1_schema_and_exemplars_match_contract() {
    let schema = report_schema();
    check_exemplar("install_report.v1.schema.json", &schema);
    let ex = report_examples();
    let doc: serde_json::Map<String, Value> = ex
        .iter()
        .map(|(k, v)| ((*k).to_string(), v.clone()))
        .collect();
    check_exemplar("install_report.v1.json", &Value::Object(doc));
    for (name, inst) in ex {
        let mut errs = vec![];
        validate(&schema, &schema, &inst, "", &mut errs);
        assert!(errs.is_empty(), "{name}: {}", errs.join("\n"));
    }
}
