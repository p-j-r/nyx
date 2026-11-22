// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Lifecycle integration with `nyx-connect` for install, upgrade, repair, and surgical uninstall.

use crate::manifest::{AiHostRecord, InstallScope};
use crate::report::HostOutcome;
use nyx_connect::{
    connect, disconnect_matching, hosts, is_cargo_build_tree, list, plain_path, snippet,
    test_server_within, Env, HostInfo, Status,
};
use std::path::Path;
use std::time::Duration;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConnectPolicy {
    Detected,
    All,
    None,
    Explicit(Vec<String>),
}

impl ConnectPolicy {
    pub fn parse(spec: &str) -> Self {
        match spec.trim().to_ascii_lowercase().as_str() {
            "" | "detected" | "auto" => Self::Detected,
            "all" => Self::All,
            "none" | "off" => Self::None,
            other => Self::Explicit(
                other
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect(),
            ),
        }
    }
}

pub fn status_label(st: &Status) -> &'static str {
    match st {
        Status::NotInstalled => "not_installed",
        Status::NotConnected => "not_connected",
        Status::Connected => "connected",
        Status::Stale { .. } => "stale",
        Status::Unreadable { .. } => "unreadable",
    }
}

fn same_exe(a: &str, b: &Path) -> bool {
    let b_s = plain_path(b.to_path_buf()).to_string_lossy().into_owned();
    if cfg!(windows) {
        a.eq_ignore_ascii_case(&b_s)
    } else {
        a == b_s
    }
}

/// Applies AI-host MCP connections for install or upgrade while honoring unreadable JSONC files,
/// local cargo developer configs, and user-disconnected hosts on upgrades.
pub fn apply_connections(
    env: Option<&Env>,
    scope: InstallScope,
    new_mcp_exe: &Path,
    old_mcp_exe: Option<&Path>,
    old_records: &[AiHostRecord],
    policy: &ConnectPolicy,
    dry_run: bool,
) -> (
    Vec<AiHostRecord>,
    Vec<HostOutcome>,
    Vec<String>,
    Vec<String>,
) {
    let mut records = Vec::new();
    let mut outcomes = Vec::new();
    let mut warnings = Vec::new();
    let mut next_steps = Vec::new();

    if matches!(policy, ConnectPolicy::None) {
        return (records, outcomes, warnings, next_steps);
    }

    if scope == InstallScope::System {
        next_steps.push(format!(
            "System install skipped per-user AI host configs; run `{} connect all` as your normal user account.",
            new_mcp_exe
                .with_file_name(if cfg!(windows) { "nyx-connect.exe" } else { "nyx-connect" })
                .display()
        ));
        return (records, outcomes, warnings, next_steps);
    }

    let detected_env;
    let env_ref = match env {
        Some(e) => e,
        None => {
            let Ok(e) = Env::detect() else {
                warnings.push("Could not detect user home directory for AI host MCP setup.".into());
                return (records, outcomes, warnings, next_steps);
            };
            detected_env = e;
            &detected_env
        }
    };

    let clean_mcp = plain_path(new_mcp_exe.to_path_buf());
    let host_infos: Vec<HostInfo> = list(env_ref, &clean_mcp);

    for info in host_infos {
        let prev_label = status_label(&info.status).to_string();
        let old_rec = old_records.iter().find(|r| r.host_id == info.id);
        let should_connect = match policy {
            ConnectPolicy::None => false,
            ConnectPolicy::Explicit(ids) => ids.iter().any(|id| id == info.id),
            ConnectPolicy::All => !matches!(
                info.status,
                Status::NotInstalled | Status::Unreadable { .. }
            ),
            ConnectPolicy::Detected => match &info.status {
                Status::NotConnected => {
                    // Upgrade safety: if the installer previously connected it and the user explicitly removed it, leave it off.
                    !old_rec.is_some_and(|r| r.connected_by_installer)
                }
                Status::Connected => true,
                Status::Stale { found } => {
                    let points_to_old_install = old_mcp_exe.is_some_and(|old| same_exe(found, old));
                    let is_dev_tree = is_cargo_build_tree(Path::new(found));
                    points_to_old_install || !is_dev_tree
                }
                Status::NotInstalled | Status::Unreadable { .. } => false,
            },
        };

        if let Status::Unreadable { why } = &info.status {
            if !matches!(policy, ConnectPolicy::None) {
                warnings.push(format!(
                    "{}: config not modified ({why}). Paste snippet:\n{}",
                    info.name,
                    snippet(&clean_mcp)
                ));
            }
            outcomes.push(HostOutcome {
                id: info.id.to_string(),
                action: "skipped".to_string(),
                status: "unreadable".to_string(),
                config_path: Some(info.config.display().to_string()),
                backup_path: None,
                detail: Some(why.clone()),
            });
            continue;
        }

        if !should_connect {
            outcomes.push(HostOutcome {
                id: info.id.to_string(),
                action: "skipped".to_string(),
                status: prev_label,
                config_path: Some(info.config.display().to_string()),
                backup_path: None,
                detail: None,
            });
            continue;
        }

        if dry_run {
            outcomes.push(HostOutcome {
                id: info.id.to_string(),
                action: "would_connect".to_string(),
                status: "connected".to_string(),
                config_path: Some(info.config.display().to_string()),
                backup_path: None,
                detail: None,
            });
            continue;
        }

        match connect(env_ref, info.id, &clean_mcp) {
            Ok(done) => {
                let backup = done.backup.as_ref().map(|p| p.display().to_string());
                records.push(AiHostRecord {
                    host_id: info.id.to_string(),
                    config_path: done.config.display().to_string(),
                    mcp_exe_path: clean_mcp.display().to_string(),
                    connected_by_installer: true,
                    previous_state: prev_label,
                    backup_path: backup.clone(),
                });
                outcomes.push(HostOutcome {
                    id: info.id.to_string(),
                    action: if done.changed {
                        "connected"
                    } else {
                        "unchanged"
                    }
                    .to_string(),
                    status: "connected".to_string(),
                    config_path: Some(done.config.display().to_string()),
                    backup_path: backup,
                    detail: Some(info.after_connect.to_string()),
                });
            }
            Err(e) => {
                warnings.push(format!("Could not connect {}: {e}", info.name));
                outcomes.push(HostOutcome {
                    id: info.id.to_string(),
                    action: "failed".to_string(),
                    status: prev_label,
                    config_path: Some(info.config.display().to_string()),
                    backup_path: None,
                    detail: Some(e),
                });
            }
        }
    }

    (records, outcomes, warnings, next_steps)
}

/// Surgically removes `"nyx"` entries from AI host configs ONLY where configured `"command"` matches `expected_mcp_exe`.
pub fn remove_matching_connections(
    env: Option<&Env>,
    expected_mcp_exe: &Path,
    dry_run: bool,
) -> (Vec<HostOutcome>, Vec<String>) {
    let mut outcomes = Vec::new();
    let mut warnings = Vec::new();

    let detected_env;
    let env_ref = match env {
        Some(e) => e,
        None => {
            let Ok(e) = Env::detect() else {
                return (outcomes, warnings);
            };
            detected_env = e;
            &detected_env
        }
    };

    let clean_mcp = plain_path(expected_mcp_exe.to_path_buf());
    for h in hosts(env_ref) {
        let st = nyx_connect::status(&h, &clean_mcp);
        match st {
            Status::Connected => {
                if dry_run {
                    outcomes.push(HostOutcome {
                        id: h.id.to_string(),
                        action: "would_disconnect".to_string(),
                        status: "not_connected".to_string(),
                        config_path: Some(h.config.display().to_string()),
                        backup_path: None,
                        detail: None,
                    });
                    continue;
                }
                match disconnect_matching(env_ref, h.id, &clean_mcp) {
                    Ok(done) => outcomes.push(HostOutcome {
                        id: h.id.to_string(),
                        action: if done.changed {
                            "disconnected"
                        } else {
                            "unchanged"
                        }
                        .to_string(),
                        status: "not_connected".to_string(),
                        config_path: Some(done.config.display().to_string()),
                        backup_path: done.backup.map(|p| p.display().to_string()),
                        detail: None,
                    }),
                    Err(e) => warnings.push(format!("Could not disconnect {}: {e}", h.name)),
                }
            }
            Status::Stale { found } => {
                outcomes.push(HostOutcome {
                    id: h.id.to_string(),
                    action: "preserved_other".to_string(),
                    status: "stale".to_string(),
                    config_path: Some(h.config.display().to_string()),
                    backup_path: None,
                    detail: Some(format!(
                        "Left entry pointing to other installation: {found}"
                    )),
                });
            }
            Status::Unreadable { why } => {
                warnings.push(format!(
                    "{} config not modified on uninstall: {why}",
                    h.name
                ));
            }
            Status::NotInstalled | Status::NotConnected => {}
        }
    }
    (outcomes, warnings)
}

/// Probes `nyx-mcp` (and transitively `nyx-script-worker`) via a real JSON-RPC `initialize` + `tools/list` handshake.
pub fn smoke_test_mcp(mcp_exe: &Path) -> Result<nyx_connect::ServerInfo, String> {
    test_server_within(mcp_exe, Duration::from_secs(10))
}
