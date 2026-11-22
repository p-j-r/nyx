// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Trait-abstracted OS platform backends for Windows, Linux, and macOS plus a pure-Rust hermetic fake backend.

pub mod linux;
pub mod macos;
pub mod windows;

use crate::manifest::{
    resolve_manifest_path, DesktopEntry, FileAssociation, InstallLayout, InstallManifest,
    InstallScope, ManifestFile, ManifestSymlink, PathEdit, TargetOs,
};
use crate::report::DoctorCheck;
use anyhow::Result;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

#[derive(Clone, Debug)]
pub struct HostPaths {
    pub home: Option<PathBuf>,
    pub appdata: Option<PathBuf>,
    pub local_appdata: Option<PathBuf>,
    pub program_files: Option<PathBuf>,
    pub program_data: Option<PathBuf>,
    pub xdg_data_home: Option<PathBuf>,
    pub shell: Option<String>,
    pub display: Option<String>,
    pub wayland_display: Option<String>,
    pub desktop_dir: Option<PathBuf>,
    pub documents_dir: Option<PathBuf>,
    pub temp_dir: PathBuf,
}

impl HostPaths {
    pub fn detect() -> Self {
        let env_path = |k: &str| {
            std::env::var_os(k)
                .filter(|v| !v.is_empty())
                .map(PathBuf::from)
        };
        let home = env_path("HOME").or_else(|| env_path("USERPROFILE"));
        Self {
            desktop_dir: home.as_ref().map(|h| h.join("Desktop")),
            documents_dir: env_path("NYX_WORKSPACE")
                .or_else(|| home.as_ref().map(|h| h.join("Documents").join("nyx"))),
            appdata: env_path("APPDATA"),
            local_appdata: env_path("LOCALAPPDATA"),
            program_files: env_path("ProgramFiles"),
            program_data: env_path("ProgramData"),
            xdg_data_home: env_path("XDG_DATA_HOME"),
            home,
            shell: std::env::var("SHELL").ok().filter(|s| !s.is_empty()),
            display: std::env::var("DISPLAY").ok().filter(|s| !s.is_empty()),
            wayland_display: std::env::var("WAYLAND_DISPLAY")
                .ok()
                .filter(|s| !s.is_empty()),
            temp_dir: std::env::temp_dir(),
        }
    }

    pub fn sandboxed(root: &Path, os: TargetOs) -> Self {
        let home = root.join("home");
        Self {
            home: Some(home.clone()),
            appdata: Some(home.join("AppData/Roaming")),
            local_appdata: Some(home.join("AppData/Local")),
            program_files: Some(root.join("ProgramFiles")),
            program_data: Some(root.join("ProgramData")),
            xdg_data_home: Some(home.join(".local/share")),
            shell: Some(
                if os == TargetOs::Macos {
                    "/bin/zsh"
                } else {
                    "/bin/bash"
                }
                .to_string(),
            ),
            display: Some(":0".to_string()),
            wayland_display: None,
            desktop_dir: Some(home.join("Desktop")),
            documents_dir: Some(home.join("Documents").join("nyx")),
            temp_dir: root.join("tmp"),
        }
    }

    pub fn connect_env(&self, os: TargetOs) -> Option<nyx_connect::Env> {
        let home = self.home.clone()?;
        Some(nyx_connect::Env {
            home,
            appdata: self.appdata.clone(),
            local_appdata: self.local_appdata.clone(),
            os: match os {
                TargetOs::Windows => nyx_connect::Os::Windows,
                TargetOs::Macos => nyx_connect::Os::Mac,
                TargetOs::Linux => nyx_connect::Os::Linux,
            },
        })
    }
}

#[derive(Clone, Debug)]
pub struct PlatformOpts {
    pub os: TargetOs,
    pub scope: InstallScope,
    pub layout: InstallLayout,
    pub add_path: bool,
    pub desktop: bool,
    pub desktop_explicit: bool,
    pub desktop_shortcut: bool,
    pub file_assoc: bool,
    pub force: bool,
    pub allow_root_user: bool,
    pub dry_run: bool,
    pub testing_reg_root: Option<String>,
}

impl PlatformOpts {
    pub fn default_for(os: TargetOs, scope: InstallScope) -> Self {
        Self {
            os,
            scope,
            layout: match os {
                TargetOs::Macos => InstallLayout::MacosApp,
                TargetOs::Windows | TargetOs::Linux => InstallLayout::Split,
            },
            add_path: true,
            desktop: true,
            desktop_explicit: false,
            desktop_shortcut: false,
            file_assoc: true,
            force: false,
            allow_root_user: false,
            dry_run: false,
            testing_reg_root: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct LayoutPaths {
    pub layout: InstallLayout,
    pub install_root: PathBuf,
    pub bin_dir: PathBuf,
    pub resource_dir: PathBuf,
    pub share_dir: PathBuf,
    pub link_dir: Option<PathBuf>,
    pub manifest_path: PathBuf,
}

pub trait PlatformBackend {
    fn os(&self) -> TargetOs;
    fn default_install_root(&self, scope: InstallScope, paths: &HostPaths) -> PathBuf;
    fn resolve_layout(&self, root: &Path, scope: InstallScope, paths: &HostPaths) -> LayoutPaths;
    fn check_elevation(&self, scope: InstallScope, allow_root_user: bool) -> Result<()>;
    fn harden_process(&self) {}

    fn install_symlinks(
        &self,
        bin_dir: &Path,
        link_dir: Option<&Path>,
        bin_names: &[String],
        force: bool,
        dry_run: bool,
    ) -> Result<Vec<ManifestSymlink>>;
    fn remove_symlinks(&self, symlinks: &[ManifestSymlink], dry_run: bool) -> Result<usize>;

    fn apply_path(
        &self,
        bin_dir: &Path,
        old_bin_dir: Option<&Path>,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<PathEdit>, Vec<String>)>;
    fn remove_path(&self, edits: &[PathEdit], opts: &PlatformOpts) -> Result<bool>;

    fn install_desktop(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<DesktopEntry>, Vec<String>)>;
    fn remove_desktop(&self, entries: &[DesktopEntry], opts: &PlatformOpts) -> Result<bool>;

    fn install_file_assoc(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<FileAssociation>, Vec<String>)>;
    fn remove_file_assoc(&self, assocs: &[FileAssociation], opts: &PlatformOpts) -> Result<bool>;

    fn finalize_bundle(
        &self,
        _manifest: &InstallManifest,
        _opts: &PlatformOpts,
        _paths: &HostPaths,
    ) -> Result<(Vec<ManifestFile>, Vec<String>)> {
        Ok((vec![], vec![]))
    }

    fn prepare_uninstall_hook(
        &self,
        _manifest: &InstallManifest,
        _opts: &PlatformOpts,
    ) -> Result<()> {
        Ok(())
    }

    fn doctor_checks(
        &self,
        manifest: Option<&InstallManifest>,
        paths: &HostPaths,
    ) -> Vec<DoctorCheck>;
}

pub fn native_backend() -> Box<dyn PlatformBackend> {
    if cfg!(target_os = "windows") {
        Box::new(windows::WindowsBackend)
    } else if cfg!(target_os = "macos") {
        Box::new(macos::MacosBackend)
    } else {
        Box::new(linux::LinuxBackend)
    }
}

pub fn compute_standard_layout(
    os: TargetOs,
    root: &Path,
    scope: InstallScope,
    paths: &HostPaths,
) -> LayoutPaths {
    match os {
        TargetOs::Windows | TargetOs::Linux => {
            let link_dir = if os == TargetOs::Linux && !matches!(scope, InstallScope::Custom) {
                match scope {
                    InstallScope::User => paths.home.as_ref().map(|h| h.join(".local/bin")),
                    InstallScope::System => Some(PathBuf::from("/usr/local/bin")),
                    InstallScope::Custom => None,
                }
            } else {
                None
            };
            LayoutPaths {
                layout: InstallLayout::Split,
                install_root: root.to_path_buf(),
                bin_dir: root.join("bin"),
                resource_dir: root.to_path_buf(),
                share_dir: root.join("share/nyx"),
                link_dir,
                manifest_path: resolve_manifest_path(os, scope, root, paths.home.as_deref()),
            }
        }
        TargetOs::Macos => {
            let link_dir = match scope {
                InstallScope::User => paths.home.as_ref().map(|h| h.join(".local/bin")),
                InstallScope::System => Some(PathBuf::from("/usr/local/bin")),
                InstallScope::Custom => None,
            };
            let contents = root.join("Contents");
            LayoutPaths {
                layout: InstallLayout::MacosApp,
                install_root: root.to_path_buf(),
                bin_dir: contents.join("MacOS"),
                resource_dir: contents.join("Resources"),
                share_dir: contents,
                link_dir,
                manifest_path: resolve_manifest_path(os, scope, root, paths.home.as_deref()),
            }
        }
    }
}

/// Shared Unix atomic symlink replacement honoring foreign files unless `force`.
pub fn create_symlink_atomic(
    target: &Path,
    link_path: &Path,
    force: bool,
    dry_run: bool,
) -> Result<bool> {
    if let Ok(meta) = std::fs::symlink_metadata(link_path) {
        if meta.file_type().is_symlink() {
            if let Ok(existing_target) = std::fs::read_link(link_path) {
                if existing_target == target {
                    return Ok(false);
                }
                let looks_like_nyx = existing_target.to_string_lossy().contains("nyx")
                    || existing_target.to_string_lossy().contains("Nyx.app");
                if !looks_like_nyx && !force {
                    anyhow::bail!(
                        "refusing to overwrite existing foreign symlink '{}' -> '{}' without --force",
                        link_path.display(),
                        existing_target.display()
                    );
                }
            }
        } else if !force {
            anyhow::bail!(
                "refusing to overwrite existing non-nyx file '{}' without --force",
                link_path.display()
            );
        }
    }
    if dry_run {
        return Ok(true);
    }
    if let Some(p) = link_path.parent() {
        std::fs::create_dir_all(p)?;
    }
    let tmp = link_path.with_extension(format!("nyx-lnk-{}", std::process::id()));
    let _ = std::fs::remove_file(&tmp);
    #[cfg(unix)]
    std::os::unix::fs::symlink(target, &tmp)?;
    #[cfg(windows)]
    {
        // On Windows in cross-OS tests we avoid requiring SeCreateSymbolicLinkPrivilege when not needed.
        std::fs::write(&tmp, format!("symlink -> {}", target.display()))?;
    }
    if let Err(e) = std::fs::rename(&tmp, link_path) {
        let _ = std::fs::remove_file(&tmp);
        return Err(e.into());
    }
    Ok(true)
}

pub fn prune_empty_xdg_parents(file_path: &Path) {
    let mut cur = file_path.parent();
    while let Some(dir) = cur {
        let name = dir.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let is_xdg_subdir = matches!(
            name,
            "apps" | "hicolor" | "icons" | "applications" | "packages" | "mime" | "share"
        ) || (name.contains('x')
            && name.chars().all(|c| c.is_ascii_digit() || c == 'x'));
        if !is_xdg_subdir {
            break;
        }
        if std::fs::remove_dir(dir).is_err() {
            break;
        }
        cur = dir.parent();
    }
}

pub fn remove_manifest_symlinks(symlinks: &[ManifestSymlink], dry_run: bool) -> Result<usize> {
    let mut removed = 0usize;
    for s in symlinks {
        let lp = Path::new(&s.link_path);
        let Ok(meta) = std::fs::symlink_metadata(lp) else {
            continue;
        };
        if meta.file_type().is_symlink() {
            if let Ok(t) = std::fs::read_link(lp) {
                if t != Path::new(&s.target_path) {
                    continue;
                }
            }
        }
        if !dry_run && std::fs::remove_file(lp).is_ok() {
            removed += 1;
        }
    }
    Ok(removed)
}

/// Shared marked-block editor for Unix shell profiles (`~/.zshrc`, `~/.bashrc`, `~/.profile`, `config.fish`).
/// Preserves symlinks (writing in place into symlink target rather than replacing the symlink with a regular file).
pub fn upsert_marked_block(
    path: &Path,
    marker_id: &str,
    block_body: &str,
    create_if_missing: bool,
    dry_run: bool,
) -> Result<(bool, bool)> {
    let start_marker = format!("# >>> {marker_id} >>>");
    let end_marker = format!("# <<< {marker_id} <<<");
    let full_block = format!("{start_marker}\n{block_body}\n{end_marker}\n");

    let existed = path.exists() || std::fs::symlink_metadata(path).is_ok();
    if !existed && !create_if_missing {
        return Ok((false, false));
    }

    let existing = if existed {
        std::fs::read_to_string(path).unwrap_or_default()
    } else {
        String::new()
    };
    let new_text = replace_or_append_block(&existing, &start_marker, &end_marker, &full_block);
    if existing == new_text && existed {
        return Ok((false, false));
    }
    if dry_run {
        return Ok((true, !existed));
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    // Preserve symlink target if `path` is a symlink (write directly to the resolved target path).
    let target_write = if std::fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink())
    {
        path.canonicalize().unwrap_or_else(|_| path.to_path_buf())
    } else {
        path.to_path_buf()
    };
    std::fs::write(target_write, new_text)?;
    Ok((true, !existed))
}

pub fn remove_marked_block(
    path: &Path,
    marker_id: &str,
    delete_if_created_empty: bool,
    dry_run: bool,
) -> Result<bool> {
    if !path.exists() {
        return Ok(false);
    }
    let start_marker = format!("# >>> {marker_id} >>>");
    let end_marker = format!("# <<< {marker_id} <<<");
    let existing = std::fs::read_to_string(path)?;
    let cleaned = strip_marked_block(&existing, &start_marker, &end_marker);
    if cleaned == existing {
        return Ok(false);
    }
    if dry_run {
        return Ok(true);
    }
    if delete_if_created_empty
        && cleaned.trim().is_empty()
        && !std::fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink())
    {
        let _ = std::fs::remove_file(path);
    } else {
        let target_write =
            if std::fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink()) {
                path.canonicalize().unwrap_or_else(|_| path.to_path_buf())
            } else {
                path.to_path_buf()
            };
        std::fs::write(target_write, cleaned)?;
    }
    Ok(true)
}

fn replace_or_append_block(
    existing: &str,
    start_marker: &str,
    end_marker: &str,
    full_block: &str,
) -> String {
    let stripped = strip_marked_block(existing, start_marker, end_marker);
    let mut out = stripped;
    if !out.is_empty() && !out.ends_with('\n') {
        out.push('\n');
    }
    out.push_str(full_block);
    out
}

fn strip_marked_block(existing: &str, start_marker: &str, end_marker: &str) -> String {
    let mut out = String::with_capacity(existing.len());
    let mut skipping = false;
    for line in existing.split_inclusive('\n') {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed == start_marker {
            skipping = true;
            continue;
        }
        if skipping {
            if trimmed == end_marker {
                skipping = false;
            }
            continue;
        }
        out.push_str(line);
    }
    out
}

/// In-memory state recorded by [`FakePlatformBackend`] for hermetic multi-OS unit and fault-injection tests.
#[derive(Clone, Debug, Default)]
pub struct FakeState {
    pub win_reg: BTreeMap<String, (u32, String)>,
    pub elevated: bool,
    pub path_applied: usize,
    pub path_removed: usize,
    pub desktop_applied: usize,
    pub desktop_removed: usize,
    pub assoc_applied: usize,
    pub assoc_removed: usize,
    pub uninstall_hook_called: bool,
    pub inject_fail_phase: Option<&'static str>,
}

/// Hermetic cross-OS backend that simulates Windows, Linux, or macOS transactions in a sandbox directory
/// without touching host registry hives, `/Applications`, or `$HOME`.
pub struct FakePlatformBackend {
    pub target_os: TargetOs,
    pub state: Arc<Mutex<FakeState>>,
}

impl FakePlatformBackend {
    pub fn new(target_os: TargetOs) -> Self {
        Self {
            target_os,
            state: Arc::new(Mutex::new(FakeState::default())),
        }
    }
}

impl PlatformBackend for FakePlatformBackend {
    fn os(&self) -> TargetOs {
        self.target_os
    }

    fn default_install_root(&self, scope: InstallScope, paths: &HostPaths) -> PathBuf {
        match self.target_os {
            TargetOs::Windows => windows::default_root(scope, paths),
            TargetOs::Linux => linux::default_root(scope, paths),
            TargetOs::Macos => macos::default_root(scope, paths),
        }
    }

    fn resolve_layout(&self, root: &Path, scope: InstallScope, paths: &HostPaths) -> LayoutPaths {
        compute_standard_layout(self.target_os, root, scope, paths)
    }

    fn check_elevation(&self, scope: InstallScope, _allow_root_user: bool) -> Result<()> {
        if scope == InstallScope::System && !self.state.lock().unwrap().elevated {
            anyhow::bail!("administrator / root privileges are required for --system installation");
        }
        Ok(())
    }

    fn install_symlinks(
        &self,
        bin_dir: &Path,
        link_dir: Option<&Path>,
        bin_names: &[String],
        force: bool,
        dry_run: bool,
    ) -> Result<Vec<ManifestSymlink>> {
        let Some(ld) = link_dir else {
            return Ok(vec![]);
        };
        let mut syms = Vec::new();
        for name in bin_names {
            let target = bin_dir.join(name);
            let link = ld.join(name);
            create_symlink_atomic(&target, &link, force, dry_run)?;
            syms.push(ManifestSymlink {
                link_path: link.display().to_string(),
                target_path: target.display().to_string(),
            });
        }
        Ok(syms)
    }

    fn remove_symlinks(&self, symlinks: &[ManifestSymlink], dry_run: bool) -> Result<usize> {
        remove_manifest_symlinks(symlinks, dry_run)
    }

    fn apply_path(
        &self,
        bin_dir: &Path,
        old_bin_dir: Option<&Path>,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<PathEdit>, Vec<String>)> {
        if let Some(fail) = self.state.lock().unwrap().inject_fail_phase {
            if fail == "path" {
                anyhow::bail!("injected failure in apply_path");
            }
        }
        if !opts.add_path || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        let mut st = self.state.lock().unwrap();
        st.path_applied += 1;
        match self.target_os {
            TargetOs::Windows => {
                let key = "HKCU\\Environment\\Path".to_string();
                let (ty, existing) = st.win_reg.get(&key).cloned().unwrap_or((
                    2, /* REG_EXPAND_SZ */
                    "%SystemRoot%\\system32".to_string(),
                ));
                let old_s = old_bin_dir.map(|p| p.display().to_string());
                let new_s = bin_dir.display().to_string();
                let (updated, changed, _) =
                    windows::edit_path_string(&existing, &new_s, old_s.as_deref());
                if changed {
                    st.win_reg.insert(key, (ty, updated));
                }
                Ok((
                    vec![PathEdit::WindowsRegistry {
                        hive: "HKCU".into(),
                        subkey: "Environment".into(),
                        value_name: "Path".into(),
                        entry_added: new_s,
                        created_value: false,
                        preserved_expand_sz: ty == 2,
                    }],
                    vec![],
                ))
            }
            TargetOs::Linux | TargetOs::Macos => {
                drop(st);
                let rc = if self.target_os == TargetOs::Macos {
                    paths.home.as_ref().map(|h| h.join(".zshrc"))
                } else {
                    paths.home.as_ref().map(|h| h.join(".bashrc"))
                };
                if let Some(p) = rc {
                    let snippet = format!(
                        "case \":${{PATH}}:\" in\n    *\":{}:\"*) ;;\n    *) export PATH=\"{}:$PATH\" ;;\nesac",
                        bin_dir.display(),
                        bin_dir.display()
                    );
                    let (_, created) =
                        upsert_marked_block(&p, "nyx installer", &snippet, true, false)?;
                    return Ok((
                        vec![PathEdit::UnixShellProfile {
                            profile_path: p.display().to_string(),
                            marker_block_id: "nyx installer".into(),
                            created_file: created,
                        }],
                        vec![],
                    ));
                }
                Ok((vec![], vec![]))
            }
        }
    }

    fn remove_path(&self, edits: &[PathEdit], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        let mut st = self.state.lock().unwrap();
        st.path_removed += 1;
        let mut changed = false;
        for e in edits {
            match e {
                PathEdit::WindowsRegistry { entry_added, .. } => {
                    let key = "HKCU\\Environment\\Path".to_string();
                    if let Some((ty, cur)) = st.win_reg.get(&key).cloned() {
                        let (updated, c) = windows::remove_from_path_string(&cur, entry_added);
                        if c {
                            st.win_reg.insert(key, (ty, updated));
                            changed = true;
                        }
                    }
                }
                PathEdit::UnixShellProfile {
                    profile_path,
                    marker_block_id,
                    created_file,
                } => {
                    if remove_marked_block(
                        Path::new(profile_path),
                        marker_block_id,
                        *created_file,
                        false,
                    )? {
                        changed = true;
                    }
                }
                PathEdit::MacosPathsD { .. } | PathEdit::LinuxProfileD { .. } => {}
            }
        }
        Ok(changed)
    }

    fn install_desktop(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<DesktopEntry>, Vec<String>)> {
        if let Some(fail) = self.state.lock().unwrap().inject_fail_phase {
            if fail == "desktop" {
                anyhow::bail!("injected failure in install_desktop");
            }
        }
        self.state.lock().unwrap().desktop_applied += 1;
        if !opts.desktop || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        match self.target_os {
            TargetOs::Linux => linux::write_xdg_desktop(manifest, opts, paths),
            TargetOs::Windows => Ok((
                vec![DesktopEntry::WindowsUninstallRegistry {
                    hive: "HKCU".into(),
                    subkey: "Software\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\nyx".into(),
                }],
                vec![],
            )),
            TargetOs::Macos => Ok((
                vec![DesktopEntry::MacosAppBundle {
                    path: manifest.install_root.clone(),
                }],
                vec![],
            )),
        }
    }

    fn remove_desktop(&self, entries: &[DesktopEntry], opts: &PlatformOpts) -> Result<bool> {
        self.state.lock().unwrap().desktop_removed += 1;
        if opts.dry_run {
            return Ok(false);
        }
        let mut removed = false;
        for e in entries {
            match e {
                DesktopEntry::FreedesktopDesktop { path }
                | DesktopEntry::FreedesktopIcon { path }
                | DesktopEntry::FreedesktopMime { path }
                | DesktopEntry::WindowsShortcut { path } => {
                    let p = Path::new(path);
                    if p.is_file() && std::fs::remove_file(p).is_ok() {
                        prune_empty_xdg_parents(p);
                        removed = true;
                    }
                }
                DesktopEntry::WindowsUninstallRegistry { .. }
                | DesktopEntry::MacosAppBundle { .. } => {
                    removed = true;
                }
            }
        }
        Ok(removed)
    }

    fn install_file_assoc(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<FileAssociation>, Vec<String>)> {
        if let Some(fail) = self.state.lock().unwrap().inject_fail_phase {
            if fail == "file_assoc" {
                anyhow::bail!("injected failure in install_file_assoc");
            }
        }
        self.state.lock().unwrap().assoc_applied += 1;
        if !opts.file_assoc || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        match self.target_os {
            TargetOs::Linux => linux::write_xdg_mime(manifest, opts, paths),
            TargetOs::Macos => Ok((
                vec![FileAssociation::MacosLaunchServices {
                    app_path: manifest.install_root.clone(),
                }],
                vec![],
            )),
            TargetOs::Windows => Ok((
                vec![FileAssociation::WindowsRegistryKey {
                    hive: "HKCU".into(),
                    subkey: "Software\\Classes\\nyx.Archive.1".into(),
                }],
                vec![],
            )),
        }
    }

    fn remove_file_assoc(&self, assocs: &[FileAssociation], opts: &PlatformOpts) -> Result<bool> {
        self.state.lock().unwrap().assoc_removed += 1;
        if opts.dry_run {
            return Ok(false);
        }
        let mut removed = false;
        for a in assocs {
            if let FileAssociation::FreedesktopMimeFile { path } = a {
                let p = Path::new(path);
                if p.is_file() && std::fs::remove_file(p).is_ok() {
                    prune_empty_xdg_parents(p);
                    removed = true;
                }
            } else {
                removed = true;
            }
        }
        Ok(removed)
    }

    fn finalize_bundle(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        _paths: &HostPaths,
    ) -> Result<(Vec<ManifestFile>, Vec<String>)> {
        if self.target_os == TargetOs::Macos && !opts.dry_run {
            let added = macos::ensure_bundle_metadata(
                Path::new(&manifest.install_root),
                &manifest.nyx_version,
            )?;
            return Ok((added, vec![]));
        }
        Ok((vec![], vec![]))
    }

    fn prepare_uninstall_hook(
        &self,
        _manifest: &InstallManifest,
        _opts: &PlatformOpts,
    ) -> Result<()> {
        self.state.lock().unwrap().uninstall_hook_called = true;
        Ok(())
    }

    fn doctor_checks(
        &self,
        _manifest: Option<&InstallManifest>,
        _paths: &HostPaths,
    ) -> Vec<DoctorCheck> {
        vec![DoctorCheck {
            name: "platform_backend".into(),
            tier: "required".into(),
            status: "ok".into(),
            detail: format!("simulated {} backend ready", self.target_os.as_str()),
        }]
    }
}
