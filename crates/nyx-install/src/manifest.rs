// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Versioned install receipt (`install-manifest.json`, `format = 1`, contract `install_manifest.v1`).

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

pub const MANIFEST_FORMAT: u32 = 1;
pub const MANIFEST_SCHEMA: &str = "nyx-install.manifest";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TargetOs {
    Windows,
    Linux,
    Macos,
}

impl TargetOs {
    pub fn current() -> Self {
        if cfg!(target_os = "windows") {
            TargetOs::Windows
        } else if cfg!(target_os = "macos") {
            TargetOs::Macos
        } else {
            TargetOs::Linux
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            TargetOs::Windows => "windows",
            TargetOs::Linux => "linux",
            TargetOs::Macos => "macos",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InstallScope {
    User,
    System,
    Custom,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InstallLayout {
    Split,
    Fhs,
    MacosApp,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FileBase {
    Bin,
    Resource,
    Share,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FileRole {
    Binary,
    Uninstaller,
    Asset,
    License,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ManifestFile {
    pub rel_path: String,
    pub base: FileBase,
    pub role: FileRole,
    pub size: u64,
    pub blake3: String,
    pub mode: u32,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ManifestSymlink {
    pub link_path: String,
    pub target_path: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PathEdit {
    WindowsRegistry {
        hive: String,
        subkey: String,
        value_name: String,
        entry_added: String,
        created_value: bool,
        preserved_expand_sz: bool,
    },
    UnixShellProfile {
        profile_path: String,
        marker_block_id: String,
        created_file: bool,
    },
    MacosPathsD {
        file_path: String,
    },
    LinuxProfileD {
        file_path: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum DesktopEntry {
    WindowsShortcut { path: String },
    WindowsUninstallRegistry { hive: String, subkey: String },
    FreedesktopDesktop { path: String },
    FreedesktopIcon { path: String },
    FreedesktopMime { path: String },
    MacosAppBundle { path: String },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum FileAssociation {
    WindowsRegistryKey {
        hive: String,
        subkey: String,
    },
    WindowsRegistryValue {
        hive: String,
        subkey: String,
        value_name: String,
        set_value: String,
        previous_value: Option<String>,
        created_key: bool,
    },
    FreedesktopMimeFile {
        path: String,
    },
    MacosLaunchServices {
        app_path: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AiHostRecord {
    pub host_id: String,
    pub config_path: String,
    pub mcp_exe_path: String,
    pub connected_by_installer: bool,
    pub previous_state: String,
    pub backup_path: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InstallManifest {
    pub format: u32,
    pub schema: String,
    pub nyx_version: String,
    pub installer_version: String,
    pub install_id: String,
    pub target_triple: String,
    pub os: TargetOs,
    pub scope: InstallScope,
    pub layout: InstallLayout,
    pub install_root: String,
    pub bin_dir: String,
    pub resource_dir: String,
    pub manifest_path: String,
    pub installed_at: String,
    pub updated_at: String,
    pub previous_version: Option<String>,
    pub files: Vec<ManifestFile>,
    pub symlinks: Vec<ManifestSymlink>,
    pub path_edits: Vec<PathEdit>,
    pub desktop_entries: Vec<DesktopEntry>,
    pub file_associations: Vec<FileAssociation>,
    pub ai_hosts: Vec<AiHostRecord>,
}

/// Formats `secs` since Unix epoch as `YYYY-MM-DDTHH:MM:SSZ` (matching `nyx_store::util::iso`).
pub fn iso_utc(secs: i64) -> String {
    let z = secs.div_euclid(86_400) + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let year = if m <= 2 { y + 1 } else { y };
    let r = secs.rem_euclid(86_400);
    format!(
        "{year:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        r / 3600,
        r % 3600 / 60,
        r % 60
    )
}

pub fn utc_now_iso() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    iso_utc(secs)
}

/// Generates a short deterministic-friendly install ID (hex timestamp + pid + sequence).
pub fn new_install_id() -> String {
    static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis() as u64)
        .unwrap_or(0);
    let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    format!(
        "{:011x}{:04x}{:02x}",
        ms & 0xfff_ffff_ffff,
        std::process::id() & 0xffff,
        n & 0xff
    )
}

/// Resolves the receipt path for an installation.
/// On macOS, receipts MUST reside outside `Nyx.app` so `_CodeSignature/CodeResources` is never invalidated.
pub fn resolve_manifest_path(
    os: TargetOs,
    scope: InstallScope,
    install_root: &Path,
    home: Option<&Path>,
) -> PathBuf {
    match os {
        TargetOs::Windows | TargetOs::Linux => install_root.join("install-manifest.json"),
        TargetOs::Macos => match scope {
            InstallScope::System if install_root == Path::new("/Applications/Nyx.app") => {
                PathBuf::from("/Library/Application Support/nyx/install/install-manifest.json")
            }
            InstallScope::User => {
                let expected = home.map(|h| h.join("Applications").join("Nyx.app"));
                if expected.as_deref() == Some(install_root) {
                    if let Some(h) = home {
                        return h
                            .join("Library/Application Support/nyx/install/install-manifest.json");
                    }
                }
                install_root
                    .parent()
                    .unwrap_or(install_root)
                    .join(".nyx-install")
                    .join("install-manifest.json")
            }
            _ => install_root
                .parent()
                .unwrap_or(install_root)
                .join(".nyx-install")
                .join("install-manifest.json"),
        },
    }
}

/// Discovers an existing `install-manifest.json` from a binary inside an installed tree or OS receipt path.
pub fn discover_manifest_from_root(
    os: TargetOs,
    scope: InstallScope,
    install_root: &Path,
    home: Option<&Path>,
) -> Option<PathBuf> {
    if os == TargetOs::Macos {
        let mut cands = vec![];
        if let Some(p) = install_root
            .parent()
            .map(|d| d.join(".nyx-install").join("install-manifest.json"))
        {
            cands.push(p);
        }
        cands.push(PathBuf::from(
            "/Library/Application Support/nyx/install/install-manifest.json",
        ));
        if let Some(h) = home {
            cands.push(h.join("Library/Application Support/nyx/install/install-manifest.json"));
        }
        for c in cands {
            if c.is_file() {
                if let Ok(text) = std::fs::read_to_string(&c) {
                    if let Ok(m) = serde_json::from_str::<InstallManifest>(&text) {
                        if Path::new(&m.install_root) == install_root {
                            return Some(c);
                        }
                    }
                }
            }
        }
    }
    let direct = resolve_manifest_path(os, scope, install_root, home);
    direct.is_file().then_some(direct)
}

/// Full path where `ManifestFile` lives on disk.
/// `FileBase::Share` resolves to `<install_root>/Contents` for `MacosApp` bundles and `<install_root>/share/nyx` otherwise.
pub fn file_abs_path(manifest: &InstallManifest, file: &ManifestFile) -> PathBuf {
    match file.base {
        FileBase::Bin => Path::new(&manifest.bin_dir).join(&file.rel_path),
        FileBase::Resource => Path::new(&manifest.resource_dir).join(&file.rel_path),
        FileBase::Share => {
            if manifest.layout == InstallLayout::MacosApp {
                Path::new(&manifest.install_root)
                    .join("Contents")
                    .join(&file.rel_path)
            } else {
                Path::new(&manifest.install_root)
                    .join("share/nyx")
                    .join(&file.rel_path)
            }
        }
    }
}

/// JSON Schema (`install_manifest.v1.schema.json`) describing `install-manifest.json`.
pub fn manifest_schema() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.dev/contracts/install_manifest.v1.schema.json",
        "title": "install_manifest.v1",
        "description": "Transactional install receipt written by nyx-setup and verified during upgrade, repair, doctor, and uninstall.",
        "type": "object",
        "additionalProperties": false,
        "required": [
            "format", "schema", "nyx_version", "installer_version", "install_id", "target_triple",
            "os", "scope", "layout", "install_root", "bin_dir", "resource_dir", "manifest_path",
            "installed_at", "updated_at", "previous_version", "files", "symlinks",
            "path_edits", "desktop_entries", "file_associations", "ai_hosts"
        ],
        "properties": {
            "format": { "const": 1 },
            "schema": { "const": MANIFEST_SCHEMA },
            "nyx_version": { "type": "string" },
            "installer_version": { "type": "string" },
            "install_id": { "type": "string" },
            "target_triple": { "type": "string" },
            "os": { "enum": ["windows", "linux", "macos"] },
            "scope": { "enum": ["user", "system", "custom"] },
            "layout": { "enum": ["split", "fhs", "macos_app"] },
            "install_root": { "type": "string" },
            "bin_dir": { "type": "string" },
            "resource_dir": { "type": "string" },
            "manifest_path": { "type": "string" },
            "installed_at": { "type": "string" },
            "updated_at": { "type": "string" },
            "previous_version": { "type": ["string", "null"] },
            "files": { "type": "array", "items": { "$ref": "#/$defs/manifest_file" } },
            "symlinks": { "type": "array", "items": { "$ref": "#/$defs/manifest_symlink" } },
            "path_edits": { "type": "array", "items": { "$ref": "#/$defs/path_edit" } },
            "desktop_entries": { "type": "array", "items": { "$ref": "#/$defs/desktop_entry" } },
            "file_associations": { "type": "array", "items": { "$ref": "#/$defs/file_association" } },
            "ai_hosts": { "type": "array", "items": { "$ref": "#/$defs/ai_host" } }
        },
        "$defs": {
            "manifest_file": {
                "type": "object",
                "additionalProperties": false,
                "required": ["rel_path", "base", "role", "size", "blake3", "mode"],
                "properties": {
                    "rel_path": { "type": "string" },
                    "base": { "enum": ["bin", "resource", "share"] },
                    "role": { "enum": ["binary", "uninstaller", "asset", "license"] },
                    "size": { "type": "integer" },
                    "blake3": { "type": "string" },
                    "mode": { "type": "integer" }
                }
            },
            "manifest_symlink": {
                "type": "object",
                "additionalProperties": false,
                "required": ["link_path", "target_path"],
                "properties": {
                    "link_path": { "type": "string" },
                    "target_path": { "type": "string" }
                }
            },
            "path_edit": {
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "hive", "subkey", "value_name", "entry_added", "created_value", "preserved_expand_sz"],
                        "properties": {
                            "kind": { "const": "windows_registry" },
                            "hive": { "type": "string" },
                            "subkey": { "type": "string" },
                            "value_name": { "type": "string" },
                            "entry_added": { "type": "string" },
                            "created_value": { "type": "boolean" },
                            "preserved_expand_sz": { "type": "boolean" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "profile_path", "marker_block_id", "created_file"],
                        "properties": {
                            "kind": { "const": "unix_shell_profile" },
                            "profile_path": { "type": "string" },
                            "marker_block_id": { "type": "string" },
                            "created_file": { "type": "boolean" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "file_path"],
                        "properties": {
                            "kind": { "const": "macos_paths_d" },
                            "file_path": { "type": "string" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "file_path"],
                        "properties": {
                            "kind": { "const": "linux_profile_d" },
                            "file_path": { "type": "string" }
                        }
                    }
                ]
            },
            "desktop_entry": {
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "windows_shortcut" }, "path": { "type": "string" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "hive", "subkey"],
                        "properties": {
                            "kind": { "const": "windows_uninstall_registry" },
                            "hive": { "type": "string" },
                            "subkey": { "type": "string" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "freedesktop_desktop" }, "path": { "type": "string" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "freedesktop_icon" }, "path": { "type": "string" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "freedesktop_mime" }, "path": { "type": "string" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "macos_app_bundle" }, "path": { "type": "string" } }
                    }
                ]
            },
            "file_association": {
                "oneOf": [
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "hive", "subkey"],
                        "properties": {
                            "kind": { "const": "windows_registry_key" },
                            "hive": { "type": "string" },
                            "subkey": { "type": "string" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "hive", "subkey", "value_name", "set_value", "previous_value", "created_key"],
                        "properties": {
                            "kind": { "const": "windows_registry_value" },
                            "hive": { "type": "string" },
                            "subkey": { "type": "string" },
                            "value_name": { "type": "string" },
                            "set_value": { "type": "string" },
                            "previous_value": { "type": ["string", "null"] },
                            "created_key": { "type": "boolean" }
                        }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "path"],
                        "properties": { "kind": { "const": "freedesktop_mime_file" }, "path": { "type": "string" } }
                    },
                    {
                        "type": "object",
                        "additionalProperties": false,
                        "required": ["kind", "app_path"],
                        "properties": { "kind": { "const": "macos_launch_services" }, "app_path": { "type": "string" } }
                    }
                ]
            },
            "ai_host": {
                "type": "object",
                "additionalProperties": false,
                "required": ["host_id", "config_path", "mcp_exe_path", "connected_by_installer", "previous_state", "backup_path"],
                "properties": {
                    "host_id": { "type": "string" },
                    "config_path": { "type": "string" },
                    "mcp_exe_path": { "type": "string" },
                    "connected_by_installer": { "type": "boolean" },
                    "previous_state": { "type": "string" },
                    "backup_path": { "type": ["string", "null"] }
                }
            }
        }
    })
}
