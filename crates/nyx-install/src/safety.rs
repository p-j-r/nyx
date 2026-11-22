// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Path-safety guardrails, poison-prefix denylist, and BLAKE3-verified non-recursive uninstall removal.
//!
//! `std::fs::remove_dir_all` is deliberately banned on install roots so a corrupted manifest or hostile
//! `--prefix` can never wipe user home folders or `Documents/nyx` workspaces.

use crate::manifest::{file_abs_path, FileRole, InstallManifest, InstallScope, TargetOs};
use anyhow::{bail, Result};
use std::collections::BTreeSet;
use std::path::{Component, Path, PathBuf};

/// Normalizes a path lexically (collapsing `.` and `..` without touching disk/symlinks) and strips Windows `\\?\`.
pub fn normalize_lexical(path: &Path) -> PathBuf {
    let stripped = nyx_connect::plain_path(path.to_path_buf());
    let mut out = PathBuf::new();
    for comp in stripped.components() {
        match comp {
            Component::Prefix(p) => out.push(p.as_os_str()),
            Component::RootDir => out.push(comp.as_os_str()),
            Component::CurDir => {}
            Component::ParentDir => {
                if !matches!(
                    out.components().next_back(),
                    Some(Component::RootDir | Component::Prefix(_)) | None
                ) {
                    out.pop();
                }
            }
            Component::Normal(s) => out.push(s),
        }
    }
    if out.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        out
    }
}

fn path_key(os: TargetOs, p: &Path) -> String {
    let norm = normalize_lexical(p);
    let s = norm.to_string_lossy().replace('\\', "/");
    let trimmed = if s.len() > 1 && !(os == TargetOs::Windows && s.len() == 3 && s.ends_with(":/"))
    {
        s.trim_end_matches('/')
    } else {
        s.as_str()
    };
    if os == TargetOs::Windows {
        trimmed.to_ascii_lowercase()
    } else {
        trimmed.to_string()
    }
}

fn push_env_path(out: &mut Vec<PathBuf>, val: Option<PathBuf>) {
    if let Some(p) = val.filter(|p| !p.as_os_str().is_empty()) {
        out.push(p);
    }
}

/// Validates that `install_root` is safe to install into or uninstall from.
/// Rejects filesystem roots, home directories, `Documents`, `Documents/nyx`, system directories,
/// and bare OS program roots without a `nyx` / `Nyx.app` subdirectory.
pub fn validate_install_root(
    root: &Path,
    os: TargetOs,
    _scope: InstallScope,
    home: Option<&Path>,
    local_appdata: Option<&Path>,
    program_files: Option<&Path>,
    allow_existing_non_nyx: bool,
) -> Result<()> {
    if root.as_os_str().is_empty() || root.as_os_str().to_string_lossy().contains('\0') {
        bail!("install path is empty or contains NUL");
    }
    if !root.is_absolute() {
        bail!("install path must be absolute: {}", root.display());
    }
    let norm = normalize_lexical(root);
    let key = path_key(os, &norm);

    if norm.parent().is_none() {
        bail!("refusing to target filesystem root {}", norm.display());
    }
    if os == TargetOs::Windows {
        let b = key.as_bytes();
        if b.len() <= 3 && b.get(1) == Some(&b':') {
            bail!("refusing to target drive root {}", norm.display());
        }
    }

    let mut forbidden: Vec<PathBuf> = vec![];
    push_env_path(&mut forbidden, home.map(Path::to_path_buf));
    if let Some(h) = home {
        forbidden.push(h.join("Documents"));
        forbidden.push(h.join("Documents").join("nyx"));
        forbidden.push(h.join("Documents").join("nyx_runs"));
        forbidden.push(h.join("Desktop"));
        forbidden.push(h.join("Downloads"));
        forbidden.push(h.join(".local"));
        forbidden.push(h.join(".local").join("share"));
        forbidden.push(h.join("Applications"));
        forbidden.push(h.join("Library"));
        forbidden.push(h.join("Library").join("Application Support"));
    }
    push_env_path(&mut forbidden, local_appdata.map(Path::to_path_buf));
    if let Some(la) = local_appdata {
        forbidden.push(la.join("Programs"));
    }
    push_env_path(&mut forbidden, program_files.map(Path::to_path_buf));
    for p in [
        "/usr",
        "/usr/local",
        "/usr/local/lib",
        "/bin",
        "/sbin",
        "/lib",
        "/etc",
        "/var",
        "/opt",
        "/Applications",
        "/Library",
        "/System",
        "/tmp",
    ] {
        forbidden.push(PathBuf::from(p));
    }
    if let Ok(win_root) = std::env::var("SystemRoot") {
        forbidden.push(PathBuf::from(&win_root));
        forbidden.push(PathBuf::from(win_root).join("System32"));
    }

    for f in forbidden {
        if path_key(os, &f) == key {
            bail!(
                "refusing protected or ambiguous install root '{}'; use a dedicated subdirectory such as '{}'",
                norm.display(),
                f.join("nyx").display()
            );
        }
    }

    // Never install into or delete a live nyx workspace directory containing `workspace.json`.
    if norm.join("workspace.json").is_file() {
        bail!(
            "refusing '{}' because it is a nyx user workspace (contains workspace.json)",
            norm.display()
        );
    }

    if !allow_existing_non_nyx && norm.exists() {
        if !norm.is_dir() {
            bail!(
                "target path exists and is not a directory: {}",
                norm.display()
            );
        }
        let has_manifest = norm.join("install-manifest.json").is_file()
            || home.is_some_and(|h| {
                os == TargetOs::Macos
                    && (norm == h.join("Applications/Nyx.app")
                        || norm == Path::new("/Applications/Nyx.app"))
            })
            || norm
                .parent()
                .is_some_and(|p| p.join(".nyx-install/install-manifest.json").is_file());
        let has_nyx_bin = norm.join("bin/nyx").is_file()
            || norm.join("bin/nyx.exe").is_file()
            || norm.join("Contents/MacOS/nyx").is_file();
        let non_empty = std::fs::read_dir(&norm)
            .map(|mut rd| rd.next().is_some())
            .unwrap_or(false);
        if non_empty && !has_manifest && !has_nyx_bin {
            bail!(
                "target directory '{}' already exists and is not an existing nyx installation (pass --force or choose an empty folder)",
                norm.display()
            );
        }
    }

    Ok(())
}

/// Validates a relative path stored in a bundle or manifest against traversal, drive prefixes,
/// alternate data streams, NUL bytes, and reserved DOS device names.
pub fn validate_rel_path(rel: &str) -> Result<PathBuf> {
    if rel.is_empty() || rel.contains('\0') {
        bail!("empty or NUL-containing entry path");
    }
    if rel.starts_with('/') || rel.starts_with('\\') || rel.contains(':') {
        bail!("unsafe entry path '{rel}'");
    }
    let normalized = rel.replace('\\', "/");
    let mut out = PathBuf::new();
    for seg in normalized.split('/') {
        if seg.is_empty() || seg == "." || seg == ".." {
            bail!("invalid path segment in '{rel}'");
        }
        let upper = seg.split('.').next().unwrap_or(seg).to_ascii_uppercase();
        if matches!(
            upper.as_str(),
            "CON"
                | "PRN"
                | "AUX"
                | "NUL"
                | "COM1"
                | "COM2"
                | "COM3"
                | "COM4"
                | "COM5"
                | "COM6"
                | "COM7"
                | "COM8"
                | "COM9"
                | "LPT1"
                | "LPT2"
                | "LPT3"
                | "LPT4"
                | "LPT5"
                | "LPT6"
                | "LPT7"
                | "LPT8"
                | "LPT9"
        ) {
            bail!("reserved device name in entry path '{rel}'");
        }
        out.push(seg);
    }
    Ok(out)
}

pub fn file_blake3_hex(path: &Path) -> std::io::Result<String> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut hasher = blake3::Hasher::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("b3:{}", hasher.finalize().to_hex()))
}

#[derive(Clone, Debug, Default)]
pub struct RemoveSummary {
    pub removed_files: usize,
    pub skipped_modified: Vec<PathBuf>,
    pub deferred_uninstaller: Option<PathBuf>,
}

/// Removes ONLY manifest-recorded files (verifying BLAKE3 digests unless `force`), then removes
/// empty manifest directories bottom-up via `std::fs::remove_dir` (never `remove_dir_all`).
pub fn safe_remove_manifest_files(
    manifest: &InstallManifest,
    force: bool,
    defer_uninstaller: bool,
) -> Result<RemoveSummary> {
    let install_root = PathBuf::from(&manifest.install_root);
    let mut summary = RemoveSummary::default();
    let mut dirs_to_prune: BTreeSet<PathBuf> = BTreeSet::new();
    dirs_to_prune.insert(install_root.clone());
    dirs_to_prune.insert(PathBuf::from(&manifest.bin_dir));
    dirs_to_prune.insert(PathBuf::from(&manifest.resource_dir));

    for entry in &manifest.files {
        validate_rel_path(&entry.rel_path)?;
        let abs = file_abs_path(manifest, entry);
        let mut cur = abs.parent();
        while let Some(d) = cur {
            if d.starts_with(&install_root) {
                dirs_to_prune.insert(d.to_path_buf());
                if d == install_root {
                    break;
                }
            } else {
                break;
            }
            cur = d.parent();
        }
        if defer_uninstaller && entry.role == FileRole::Uninstaller {
            summary.deferred_uninstaller = Some(abs);
            continue;
        }
        if !abs.exists() && std::fs::symlink_metadata(&abs).is_err() {
            continue;
        }
        if !force {
            if let Ok(actual) = file_blake3_hex(&abs) {
                if actual != entry.blake3 {
                    summary.skipped_modified.push(abs);
                    continue;
                }
            }
        }
        std::fs::remove_file(&abs)?;
        summary.removed_files += 1;
    }

    // Remove the manifest itself if inside install_root, or parent receipt dir if outside (macOS).
    let manifest_p = PathBuf::from(&manifest.manifest_path);
    if manifest_p.is_file() {
        let _ = std::fs::remove_file(&manifest_p);
        if let Some(mp) = manifest_p.parent() {
            if mp != install_root {
                let _ = std::fs::remove_dir(mp);
            }
        }
    }

    // Bottom-up directory cleanup: deepest paths first, ignoring DirectoryNotEmpty.
    let mut sorted_dirs: Vec<PathBuf> = dirs_to_prune.into_iter().collect();
    sorted_dirs.sort_by(|a, b| {
        b.components()
            .count()
            .cmp(&a.components().count())
            .then(b.cmp(a))
    });
    for d in sorted_dirs {
        let _ = std::fs::remove_dir(&d);
    }

    Ok(summary)
}
