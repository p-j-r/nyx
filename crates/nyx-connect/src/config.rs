// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Reading and rewriting a host's JSON config: the status of the `nyx` entry, adding it, removing it.
//!
//! Rules that hold for every host (pinned by `tests/behaviour.rs`):
//! * only the `nyx` entry is touched; every other key and server stays, in its order (`serde_json` `preserve_order`);
//! * a file that is not plain JSON (comments, a syntax error) is never written;
//! * the first time nyx changes a file, the original is copied next to it as `<name>.nyx-backup`; that copy is never overwritten, so it is
//!   always the file as it was before nyx touched it;
//! * the new file is written to a temporary file in the same folder and renamed over the old one, so a failure leaves no half config.
use crate::env::Env;
use crate::hosts::{hosts, Host};
use crate::ENTRY;
use serde_json::{Map, Value};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Status {
    /// The host's config folder does not exist.
    NotInstalled,
    /// The host is there, nyx is not in its config.
    NotConnected,
    /// nyx is in the config and points at this exe.
    Connected,
    /// nyx is in the config but points somewhere else (an older build, another folder).
    Stale { found: String },
    /// The config exists but is not plain JSON (comments, a syntax error), so nyx will not touch it.
    Unreadable { why: String },
}

#[derive(Clone, Debug)]
pub struct HostInfo {
    pub id: &'static str,
    pub name: &'static str,
    pub config: PathBuf,
    pub after_connect: &'static str,
    pub status: Status,
}

/// What a [`connect`] or [`disconnect`] did.
#[derive(Clone, Debug)]
pub struct Done {
    pub host: &'static str,
    pub config: PathBuf,
    /// The copy of the file as it was before nyx first changed it, when the file existed then.
    pub backup: Option<PathBuf>,
    pub changed: bool,
}

pub fn status(host: &Host, exe: &Path) -> Status {
    let cfg = match read_config(&host.config) {
        Ok(c) => c,
        Err(why) => return Status::Unreadable { why },
    };
    let entry = cfg.as_ref().and_then(|c| c.get(host.shape.key())).and_then(|s| s.get(ENTRY));
    match entry {
        Some(e) => {
            let found = e.get("command").and_then(Value::as_str).unwrap_or("").to_string();
            if same_path(&found, exe) {
                Status::Connected
            } else {
                Status::Stale { found }
            }
        }
        None if cfg.is_some() || host.marker.is_dir() => Status::NotConnected,
        None => Status::NotInstalled,
    }
}

/// Every host with its status, for a settings panel.
pub fn list(env: &Env, exe: &Path) -> Vec<HostInfo> {
    hosts(env).into_iter().map(|h| HostInfo { id: h.id, name: h.name, after_connect: h.after_connect, status: status(&h, exe), config: h.config.clone() }).collect()
}

fn same_path(a: &str, b: &Path) -> bool {
    let b = b.to_string_lossy();
    if cfg!(windows) {
        a.eq_ignore_ascii_case(&b)
    } else {
        a == b
    }
}

/// `Ok(None)`: no file or an empty file. `Err`: present but not plain JSON, or not an object.
fn read_config(path: &Path) -> Result<Option<Value>, String> {
    let text = match fs::read_to_string(path) {
        Ok(t) => t,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(format!("cannot read {}: {e}", path.display())),
    };
    let text = text.trim_start_matches('\u{feff}');
    if text.trim().is_empty() {
        return Ok(None);
    }
    match serde_json::from_str::<Value>(text) {
        Ok(v) if v.is_object() => Ok(Some(v)),
        Ok(_) => Err(format!("{} is JSON but not an object; not changed", path.display())),
        Err(e) => Err(format!("{} is not plain JSON ({e}); it may contain comments. Not changed: add the entry by hand (see the snippet)", path.display())),
    }
}

fn find<'a>(hs: &'a [Host], id: &str) -> Result<&'a Host, String> {
    hs.iter().find(|h| h.id == id).ok_or_else(|| format!("no host '{id}'; valid: {}", hs.iter().map(|h| h.id).collect::<Vec<_>>().join(", ")))
}

/// Adds or updates the `nyx` entry. Other servers and settings in the file are kept, in their order, and so are keys the user added to the
/// nyx entry itself.
pub fn connect(env: &Env, id: &str, exe: &Path) -> Result<Done, String> {
    let hs = hosts(env);
    let host = find(&hs, id)?;
    if !exe.is_file() {
        return Err(format!("the nyx-mcp program is not at {}; build it first (cargo build --release -p nyx-mcp)", exe.display()));
    }
    let mut cfg = read_config(&host.config)?.unwrap_or_else(|| Value::Object(Map::new()));
    let key = host.shape.key();
    let servers = cfg.as_object_mut().expect("checked to be an object").entry(key).or_insert_with(|| Value::Object(Map::new()));
    let servers = servers.as_object_mut().ok_or_else(|| format!("'{key}' in {} is not an object; not touching it", host.config.display()))?;
    let want = host.shape.entry_over(servers.get(ENTRY), exe);
    if servers.get(ENTRY) == Some(&want) {
        return Ok(Done { host: host.id, config: host.config.clone(), backup: existing_backup(&host.config), changed: false });
    }
    servers.insert(ENTRY.to_string(), want);
    let backup = write_with_backup(&host.config, &cfg)?;
    Ok(Done { host: host.id, config: host.config.clone(), backup, changed: true })
}

/// Removes the `nyx` entry and leaves everything else.
pub fn disconnect(env: &Env, id: &str) -> Result<Done, String> {
    let hs = hosts(env);
    let host = find(&hs, id)?;
    let unchanged = |host: &Host| Done { host: host.id, config: host.config.clone(), backup: existing_backup(&host.config), changed: false };
    let Some(mut cfg) = read_config(&host.config)? else {
        return Ok(unchanged(host));
    };
    let removed = cfg.get_mut(host.shape.key()).and_then(Value::as_object_mut).map(|s| s.shift_remove(ENTRY).is_some()).unwrap_or(false);
    if !removed {
        return Ok(unchanged(host));
    }
    let backup = write_with_backup(&host.config, &cfg)?;
    Ok(Done { host: host.id, config: host.config.clone(), backup, changed: true })
}

fn backup_path(path: &Path) -> PathBuf {
    let mut name = path.file_name().unwrap_or_default().to_os_string();
    name.push(".nyx-backup");
    path.with_file_name(name)
}

fn existing_backup(path: &Path) -> Option<PathBuf> {
    Some(backup_path(path)).filter(|b| b.is_file())
}

/// Copies the old file to `<name>.nyx-backup` unless that copy exists, then replaces it by writing a temp file next to it and renaming.
/// Returns the backup, if there is one.
fn write_with_backup(path: &Path, cfg: &Value) -> Result<Option<PathBuf>, String> {
    let dir = path.parent().ok_or("config path has no folder")?;
    fs::create_dir_all(dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    let backup = backup_path(path);
    if path.is_file() && !backup.exists() {
        fs::copy(path, &backup).map_err(|e| format!("cannot back up {}: {e}", path.display()))?;
    }
    let mut tmp_name = path.file_name().unwrap_or_default().to_os_string();
    tmp_name.push(".nyx-tmp");
    let tmp = path.with_file_name(tmp_name);
    let text = serde_json::to_string_pretty(cfg).map_err(|e| e.to_string())? + "\n";
    fs::write(&tmp, text).map_err(|e| format!("cannot write {}: {e}", tmp.display()))?;
    fs::rename(&tmp, path).map_err(|e| {
        let _ = fs::remove_file(&tmp);
        format!("cannot replace {}: {e}", path.display())
    })?;
    Ok(existing_backup(path))
}
