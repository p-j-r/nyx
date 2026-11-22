// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Linux platform driver: XDG / FHS layouts, atomic `~/.local/bin` / `/usr/local/bin` symlinks,
//! marked shell profile blocks, FreeDesktop `.desktop` / `hicolor` PNG icons / `shared-mime-info` XML,
//! and two-tier `dlopen` GUI readiness diagnostics.

use super::{
    compute_standard_layout, create_symlink_atomic, prune_empty_xdg_parents,
    remove_manifest_symlinks, remove_marked_block, upsert_marked_block, HostPaths, LayoutPaths,
    PlatformBackend, PlatformOpts,
};
use crate::manifest::{
    DesktopEntry, FileAssociation, InstallManifest, InstallScope, ManifestSymlink, PathEdit,
    TargetOs,
};
use crate::report::DoctorCheck;
#[cfg(unix)]
use anyhow::bail;
use anyhow::Result;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub const MARK_PNG: &[u8] = include_bytes!("../../../../docs/brand/nyx-mark.png");
const MARKER_ID: &str = "nyx installer";

pub fn default_root(scope: InstallScope, paths: &HostPaths) -> PathBuf {
    match scope {
        InstallScope::System => PathBuf::from("/usr/local/lib/nyx"),
        InstallScope::User | InstallScope::Custom => paths
            .xdg_data_home
            .clone()
            .or_else(|| paths.home.as_ref().map(|h| h.join(".local/share")))
            .unwrap_or_else(|| PathBuf::from("/tmp/nyx"))
            .join("nyx"),
    }
}

pub fn render_desktop_entry(bin_dir: &Path) -> String {
    let nyx_bin = bin_dir.join("nyx");
    format!(
        "[Desktop Entry]\n\
         Type=Application\n\
         Version=1.5\n\
         Name=nyx\n\
         GenericName=Finite-Element Simulation\n\
         Comment=Photonics, phononics, electrostatics and electromagnetics finite-element solver\n\
         TryExec={exe}\n\
         Exec=\"{exe}\" %f\n\
         Icon=nyx\n\
         Terminal=false\n\
         Categories=Science;Physics;Engineering;NumericalAnalysis;\n\
         MimeType=application/x-nyx-project;application/x-nyx-archive;text/x-nyx-script;\n\
         StartupWMClass=nyx\n\
         Keywords=FEM;photonics;phononics;Maxwell;waveguide;cavity;cQED;\n",
        exe = nyx_bin.display()
    )
}

pub fn render_mime_xml() -> &'static str {
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
     <mime-info xmlns=\"http://www.freedesktop.org/standards/shared-mime-info\">\n\
       <mime-type type=\"application/x-nyx-project\">\n\
         <comment>nyx Simulation Project</comment>\n\
         <sub-class-of type=\"application/json\"/>\n\
         <glob pattern=\"*.nyx.json\" weight=\"60\"/>\n\
       </mime-type>\n\
       <mime-type type=\"application/x-nyx-archive\">\n\
         <comment>nyx Workflow Archive</comment>\n\
         <sub-class-of type=\"application/zip\"/>\n\
         <glob pattern=\"*.nyx\"/>\n\
       </mime-type>\n\
       <mime-type type=\"text/x-nyx-script\">\n\
         <comment>nyx Simulation Script</comment>\n\
         <sub-class-of type=\"text/x-python\"/>\n\
         <glob pattern=\"*.nys\"/>\n\
       </mime-type>\n\
     </mime-info>\n"
}

fn xdg_data_root(opts: &PlatformOpts, paths: &HostPaths, install_root: &Path) -> PathBuf {
    match opts.scope {
        InstallScope::System => PathBuf::from("/usr/local/share"),
        InstallScope::User => paths
            .xdg_data_home
            .clone()
            .or_else(|| paths.home.as_ref().map(|h| h.join(".local/share")))
            .unwrap_or_else(|| PathBuf::from("/tmp")),
        InstallScope::Custom => install_root.join("share"),
    }
}

pub fn write_xdg_desktop(
    manifest: &InstallManifest,
    opts: &PlatformOpts,
    paths: &HostPaths,
) -> Result<(Vec<DesktopEntry>, Vec<String>)> {
    let mut entries = Vec::new();
    let warnings = Vec::new();
    let xdg_root = xdg_data_root(opts, paths, Path::new(&manifest.install_root));
    // Auto-skip on headless Linux unless `--desktop` was passed or a desktop session / applications folder exists.
    if opts.scope == InstallScope::User
        && !opts.desktop_explicit
        && paths.display.is_none()
        && paths.wayland_display.is_none()
        && !xdg_root.join("applications").is_dir()
    {
        return Ok((entries, warnings));
    }

    let desktop_path = xdg_root.join("applications/nyx.desktop");
    if let Some(p) = desktop_path.parent() {
        std::fs::create_dir_all(p)?;
    }
    let exec_dir = if opts.scope == InstallScope::System {
        PathBuf::from("/usr/local/bin")
    } else if let Some(h) = &paths.home {
        if opts.scope == InstallScope::User {
            h.join(".local/bin")
        } else {
            PathBuf::from(&manifest.bin_dir)
        }
    } else {
        PathBuf::from(&manifest.bin_dir)
    };
    std::fs::write(&desktop_path, render_desktop_entry(&exec_dir))?;
    entries.push(DesktopEntry::FreedesktopDesktop {
        path: desktop_path.display().to_string(),
    });

    if let Ok(img) = image::load_from_memory(MARK_PNG) {
        let rgba = img.to_rgba8();
        for size in [32u32, 64, 128, 256, 512] {
            let icon_path = xdg_root.join(format!("icons/hicolor/{size}x{size}/apps/nyx.png"));
            if let Some(p) = icon_path.parent() {
                std::fs::create_dir_all(p)?;
            }
            let resized =
                image::imageops::resize(&rgba, size, size, image::imageops::FilterType::Lanczos3);
            if resized.save(&icon_path).is_ok() {
                entries.push(DesktopEntry::FreedesktopIcon {
                    path: icon_path.display().to_string(),
                });
            }
        }
    }

    best_effort_refresh(&xdg_root);
    Ok((entries, warnings))
}

pub fn write_xdg_mime(
    manifest: &InstallManifest,
    opts: &PlatformOpts,
    paths: &HostPaths,
) -> Result<(Vec<FileAssociation>, Vec<String>)> {
    let xdg_root = xdg_data_root(opts, paths, Path::new(&manifest.install_root));
    if opts.scope == InstallScope::User
        && !opts.desktop_explicit
        && paths.display.is_none()
        && paths.wayland_display.is_none()
        && !xdg_root.join("applications").is_dir()
    {
        return Ok((vec![], vec![]));
    }
    let xml_path = xdg_root.join("mime/packages/nyx.xml");
    if let Some(p) = xml_path.parent() {
        std::fs::create_dir_all(p)?;
    }
    std::fs::write(&xml_path, render_mime_xml())?;
    best_effort_refresh(&xdg_root);
    Ok((
        vec![FileAssociation::FreedesktopMimeFile {
            path: xml_path.display().to_string(),
        }],
        vec![],
    ))
}

fn best_effort_refresh(xdg_root: &Path) {
    for (tool, arg) in [
        ("update-desktop-database", xdg_root.join("applications")),
        ("update-mime-database", xdg_root.join("mime")),
        ("gtk-update-icon-cache", xdg_root.join("icons/hicolor")),
    ] {
        let mut cmd = Command::new(tool);
        if tool == "gtk-update-icon-cache" {
            cmd.args(["-q", "-f", "-t"]);
        }
        cmd.arg(&arg).stdout(Stdio::null()).stderr(Stdio::null());
        if let Ok(mut child) = cmd.spawn() {
            let start = std::time::Instant::now();
            while start.elapsed() < std::time::Duration::from_secs(2) {
                match child.try_wait() {
                    Ok(Some(_)) => break,
                    Ok(None) => std::thread::sleep(std::time::Duration::from_millis(20)),
                    Err(_) => break,
                }
            }
            let _ = child.kill();
        }
    }
}

pub fn posix_path_block(bin_or_link_dir: &str) -> String {
    format!(
        "case \":${{PATH}}:\" in\n\
         \x20   *\":{dir}:\"*) ;;\n\
         \x20   *) export PATH=\"{dir}:$PATH\" ;;\n\
         esac",
        dir = bin_or_link_dir
    )
}

pub fn fish_path_block(bin_or_link_dir: &str) -> String {
    format!("fish_add_path -g -p \"{bin_or_link_dir}\"")
}

/// Targets existing shell rc files without ever creating `~/.bash_profile` from scratch.
pub fn select_shell_profiles(home: &Path, shell: Option<&str>) -> Vec<(PathBuf, bool, bool)> {
    // (path, is_fish, create_if_missing)
    let mut targets = Vec::new();
    let shell_name = shell
        .and_then(|s| Path::new(s).file_name())
        .and_then(|n| n.to_str())
        .unwrap_or("");
    let bashrc = home.join(".bashrc");
    let profile = home.join(".profile");
    let bash_profile = home.join(".bash_profile");
    let zshrc = home.join(".zshrc");
    let fish_cfg = home.join(".config/fish/config.fish");

    if bashrc.exists() || std::fs::symlink_metadata(&bashrc).is_ok() {
        targets.push((bashrc, false, false));
    }
    if bash_profile.exists() || std::fs::symlink_metadata(&bash_profile).is_ok() {
        targets.push((bash_profile, false, false));
    }
    if profile.exists()
        || std::fs::symlink_metadata(&profile).is_ok()
        || shell_name == "bash"
        || shell_name == "sh"
    {
        targets.push((profile, false, shell_name == "bash" || shell_name == "sh"));
    }
    if zshrc.exists() || std::fs::symlink_metadata(&zshrc).is_ok() || shell_name == "zsh" {
        targets.push((zshrc, false, shell_name == "zsh"));
    }
    if fish_cfg.exists() || std::fs::symlink_metadata(&fish_cfg).is_ok() || shell_name == "fish" {
        targets.push((fish_cfg, true, shell_name == "fish"));
    }
    targets
}

pub struct LinuxBackend;

impl PlatformBackend for LinuxBackend {
    fn os(&self) -> TargetOs {
        TargetOs::Linux
    }

    fn default_install_root(&self, scope: InstallScope, paths: &HostPaths) -> PathBuf {
        default_root(scope, paths)
    }

    fn resolve_layout(&self, root: &Path, scope: InstallScope, paths: &HostPaths) -> LayoutPaths {
        compute_standard_layout(TargetOs::Linux, root, scope, paths)
    }

    fn check_elevation(&self, scope: InstallScope, allow_root_user: bool) -> Result<()> {
        #[cfg(unix)]
        {
            let euid = unsafe { libc::geteuid() };
            match scope {
                InstallScope::System if euid != 0 => {
                    bail!("root privileges are required for --system installation (run with sudo or omit --system for ~/.local install)");
                }
                InstallScope::User if euid == 0 && !allow_root_user => {
                    bail!("refusing --user install while running as root (would create root-owned files in user home; pass --system or --allow-root-user)");
                }
                _ => {}
            }
        }
        let _ = (scope, allow_root_user);
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
        _old_bin_dir: Option<&Path>,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<PathEdit>, Vec<String>)> {
        if !opts.add_path {
            return Ok((vec![], vec![]));
        }
        let mut edits = Vec::new();
        let warnings = Vec::new();
        match opts.scope {
            InstallScope::System => {
                // Standard /usr/local/bin is already on default PATH.
                Ok((edits, warnings))
            }
            InstallScope::User => {
                let Some(home) = &paths.home else {
                    return Ok((edits, warnings));
                };
                let path_expr = "$HOME/.local/bin";
                for (profile_path, is_fish, create_if_missing) in
                    select_shell_profiles(home, paths.shell.as_deref())
                {
                    let snippet = if is_fish {
                        fish_path_block(path_expr)
                    } else {
                        posix_path_block(path_expr)
                    };
                    let (changed, created) = upsert_marked_block(
                        &profile_path,
                        MARKER_ID,
                        &snippet,
                        create_if_missing,
                        opts.dry_run,
                    )?;
                    if changed {
                        edits.push(PathEdit::UnixShellProfile {
                            profile_path: profile_path.display().to_string(),
                            marker_block_id: MARKER_ID.to_string(),
                            created_file: created,
                        });
                    }
                }
                Ok((edits, warnings))
            }
            InstallScope::Custom => {
                let _ = bin_dir;
                Ok((edits, warnings))
            }
        }
    }

    fn remove_path(&self, edits: &[PathEdit], opts: &PlatformOpts) -> Result<bool> {
        let mut changed = false;
        for edit in edits {
            match edit {
                PathEdit::UnixShellProfile {
                    profile_path,
                    marker_block_id,
                    created_file,
                } => {
                    if remove_marked_block(
                        Path::new(profile_path),
                        marker_block_id,
                        *created_file,
                        opts.dry_run,
                    )? {
                        changed = true;
                    }
                }
                PathEdit::LinuxProfileD { file_path }
                    if Path::new(file_path).is_file()
                        && !opts.dry_run
                        && std::fs::remove_file(file_path).is_ok() =>
                {
                    changed = true;
                }
                _ => {}
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
        if !opts.desktop || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        write_xdg_desktop(manifest, opts, paths)
    }

    fn remove_desktop(&self, entries: &[DesktopEntry], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        let mut removed = false;
        for entry in entries {
            if let DesktopEntry::FreedesktopDesktop { path }
            | DesktopEntry::FreedesktopIcon { path }
            | DesktopEntry::FreedesktopMime { path } = entry
            {
                let p = Path::new(path);
                if p.is_file() && std::fs::remove_file(p).is_ok() {
                    prune_empty_xdg_parents(p);
                    removed = true;
                }
            }
        }
        if let Some(DesktopEntry::FreedesktopDesktop { path }) = entries.first() {
            if let Some(xdg) = Path::new(path).parent().and_then(Path::parent) {
                best_effort_refresh(xdg);
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
        if !opts.file_assoc || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        write_xdg_mime(manifest, opts, paths)
    }

    fn remove_file_assoc(&self, assocs: &[FileAssociation], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        let mut changed = false;
        let mut refresh_dir: Option<PathBuf> = None;
        for a in assocs {
            if let FileAssociation::FreedesktopMimeFile { path } = a {
                let p = Path::new(path);
                refresh_dir = p
                    .parent()
                    .and_then(Path::parent)
                    .and_then(Path::parent)
                    .map(Path::to_path_buf);
                if p.is_file() && std::fs::remove_file(p).is_ok() {
                    prune_empty_xdg_parents(p);
                    changed = true;
                }
            }
        }
        if let Some(xdg) = refresh_dir {
            best_effort_refresh(&xdg);
        }
        Ok(changed)
    }

    fn doctor_checks(
        &self,
        manifest: Option<&InstallManifest>,
        paths: &HostPaths,
    ) -> Vec<DoctorCheck> {
        let mut checks = Vec::new();
        if let Some(m) = manifest {
            let bin_dir = Path::new(&m.bin_dir);
            let ok = ["nyx", "nyx-mcp", "nyx-script-worker", "nyx-connect"]
                .iter()
                .all(|b| bin_dir.join(b).is_file());
            checks.push(DoctorCheck {
                name: "binaries".into(),
                tier: "required".into(),
                status: if ok { "ok" } else { "fail" }.into(),
                detail: if ok {
                    format!("core binaries present in {}", bin_dir.display())
                } else {
                    format!("missing one or more binaries in {}", bin_dir.display())
                },
            });
        }
        let has_display = paths.display.is_some() || paths.wayland_display.is_some();
        #[cfg(target_os = "linux")]
        {
            let gl_ok = dlopen_probe(&["libGL.so.1", "libEGL.so.1"]);
            let x_or_wl_ok = dlopen_probe(&["libwayland-client.so.0", "libX11.so.6"]);
            checks.push(DoctorCheck {
                name: "gui_runtime_libs".into(),
                tier: "optional".into(),
                status: if gl_ok && x_or_wl_ok {
                    "ok"
                } else if has_display {
                    "warn"
                } else {
                    "info"
                }
                .into(),
                detail: format!(
                    "GL/EGL present: {gl_ok}; Wayland/X11 present: {x_or_wl_ok} (headless CLI and MCP work without GUI libs)"
                ),
            });
        }
        #[cfg(not(target_os = "linux"))]
        {
            checks.push(DoctorCheck {
                name: "gui_runtime_libs".into(),
                tier: "optional".into(),
                status: if has_display { "ok" } else { "info" }.into(),
                detail: "Linux runtime check skipped on non-Linux host".into(),
            });
        }
        checks
    }
}

#[cfg(target_os = "linux")]
fn dlopen_probe(names: &[&str]) -> bool {
    use std::ffi::CString;
    for name in names {
        let Ok(c) = CString::new(*name) else { continue };
        let h = unsafe { libc::dlopen(c.as_ptr(), libc::RTLD_LAZY | libc::RTLD_LOCAL) };
        if !h.is_null() {
            unsafe {
                libc::dlclose(h);
            }
            return true;
        }
    }
    false
}
