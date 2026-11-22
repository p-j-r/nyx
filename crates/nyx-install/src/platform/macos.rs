// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! macOS platform driver: signed-bundle-preserving `Nyx.app` layout, external `Application Support`
//! manifest receipt, pure-Rust `.icns` + `Info.plist` generators, `~/.local/bin` / `/usr/local/bin`
//! symlinks, `/etc/paths.d/nyx`, and `lsregister` LaunchServices integration.

use super::{
    compute_standard_layout, create_symlink_atomic,
    linux::{fish_path_block, posix_path_block, MARK_PNG},
    remove_manifest_symlinks, remove_marked_block, upsert_marked_block, HostPaths, LayoutPaths,
    PlatformBackend, PlatformOpts,
};
use crate::manifest::{
    DesktopEntry, FileAssociation, FileBase, FileRole, InstallManifest, InstallScope, ManifestFile,
    ManifestSymlink, PathEdit, TargetOs,
};
use crate::report::DoctorCheck;
use crate::safety::file_blake3_hex;
#[cfg(unix)]
use anyhow::bail;
use anyhow::Result;
use std::path::{Path, PathBuf};
#[cfg(target_os = "macos")]
use std::process::Command;

const MARKER_ID: &str = "nyx installer";

pub fn default_root(scope: InstallScope, paths: &HostPaths) -> PathBuf {
    match scope {
        InstallScope::System => PathBuf::from("/Applications/Nyx.app"),
        InstallScope::User | InstallScope::Custom => paths
            .home
            .clone()
            .unwrap_or_else(|| PathBuf::from("/tmp"))
            .join("Applications")
            .join("Nyx.app"),
    }
}

/// Pure-Rust Apple `.icns` builder producing embedded PNG icon chunks (`ic07` 128, `ic08` 256, `ic09` 512, `ic10` 1024).
pub fn build_icns(mark_png: &[u8]) -> Result<Vec<u8>> {
    let img = image::load_from_memory(mark_png)?.to_rgba8();
    let specs: [([u8; 4], u32); 4] = [
        (*b"ic07", 128),
        (*b"ic08", 256),
        (*b"ic09", 512),
        (*b"ic10", 1024),
    ];
    let mut chunks = Vec::new();
    for (tag, dim) in specs {
        let resized =
            image::imageops::resize(&img, dim, dim, image::imageops::FilterType::Lanczos3);
        let mut png_bytes = Vec::new();
        let dyn_img = image::DynamicImage::ImageRgba8(resized);
        dyn_img.write_to(
            &mut std::io::Cursor::new(&mut png_bytes),
            image::ImageFormat::Png,
        )?;
        let chunk_len = (8 + png_bytes.len()) as u32;
        chunks.extend_from_slice(&tag);
        chunks.extend_from_slice(&chunk_len.to_be_bytes());
        chunks.extend_from_slice(&png_bytes);
    }
    let total_len = (8 + chunks.len()) as u32;
    let mut out = Vec::with_capacity(total_len as usize);
    out.extend_from_slice(b"icns");
    out.extend_from_slice(&total_len.to_be_bytes());
    out.extend_from_slice(&chunks);
    Ok(out)
}

pub fn render_info_plist(
    version: &str,
    exe_name: &str,
    bundle_name: &str,
    bundle_id: &str,
    include_docs: bool,
) -> String {
    let doc_section = if include_docs {
        "  <key>CFBundleDocumentTypes</key>\n\
         \x20 <array>\n\
         \x20   <dict>\n\
         \x20     <key>CFBundleTypeName</key><string>nyx Simulation Project</string>\n\
         \x20     <key>CFBundleTypeRole</key><string>Editor</string>\n\
         \x20     <key>LSHandlerRank</key><string>Owner</string>\n\
         \x20     <key>CFBundleTypeExtensions</key><array><string>nyx.json</string></array>\n\
         \x20   </dict>\n\
         \x20   <dict>\n\
         \x20     <key>CFBundleTypeName</key><string>nyx Workflow Archive</string>\n\
         \x20     <key>CFBundleTypeRole</key><string>Editor</string>\n\
         \x20     <key>LSHandlerRank</key><string>Owner</string>\n\
         \x20     <key>CFBundleTypeExtensions</key><array><string>nyx</string></array>\n\
         \x20   </dict>\n\
         \x20   <dict>\n\
         \x20     <key>CFBundleTypeName</key><string>nyx Simulation Script</string>\n\
         \x20     <key>CFBundleTypeRole</key><string>Editor</string>\n\
         \x20     <key>LSHandlerRank</key><string>Owner</string>\n\
         \x20     <key>CFBundleTypeExtensions</key><array><string>nys</string></array>\n\
         \x20   </dict>\n\
         \x20   <dict>\n\
         \x20     <key>CFBundleTypeName</key><string>nyx Results File (HDF5)</string>\n\
         \x20     <key>CFBundleTypeRole</key><string>Viewer</string>\n\
         \x20     <key>LSHandlerRank</key><string>Alternate</string>\n\
         \x20     <key>CFBundleTypeExtensions</key><array><string>h5</string></array>\n\
         \x20   </dict>\n\
         \x20 </array>\n"
    } else {
        ""
    };
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <!DOCTYPE plist PUBLIC \"-//Apple//DTD PLIST 1.0//EN\" \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">\n\
         <plist version=\"1.0\">\n\
         <dict>\n\
           <key>CFBundlePackageType</key><string>APPL</string>\n\
           <key>CFBundleName</key><string>{bundle_name}</string>\n\
           <key>CFBundleDisplayName</key><string>{bundle_name}</string>\n\
           <key>CFBundleIdentifier</key><string>{bundle_id}</string>\n\
           <key>CFBundleVersion</key><string>{version}</string>\n\
           <key>CFBundleShortVersionString</key><string>{version}</string>\n\
           <key>CFBundleExecutable</key><string>{exe_name}</string>\n\
           <key>CFBundleIconFile</key><string>Nyx</string>\n\
           <key>NSHighResolutionCapable</key><true/>\n\
           <key>LSMinimumSystemVersion</key><string>12.0</string>\n\
         {doc_section}\
         </dict>\n\
         </plist>\n"
    )
}

/// Ensures unsigned local/CI `.app` bundles carry `Contents/Info.plist` and `Contents/Resources/Nyx.icns`,
/// returning `ManifestFile` records for any newly synthesized metadata files so uninstall cleans the entire bundle.
/// Skips mutating any bundle that already has `Contents/_CodeSignature`.
pub fn ensure_bundle_metadata(app_root: &Path, version: &str) -> Result<Vec<ManifestFile>> {
    let contents = app_root.join("Contents");
    if contents.join("_CodeSignature").exists() {
        return Ok(vec![]);
    }
    let res_dir = contents.join("Resources");
    std::fs::create_dir_all(&res_dir)?;
    let mut created = Vec::new();
    let plist_path = contents.join("Info.plist");
    if !plist_path.is_file() {
        let bytes = render_info_plist(version, "nyx", "nyx", "dev.pjr.nyx", true).into_bytes();
        std::fs::write(&plist_path, &bytes)?;
        created.push(ManifestFile {
            rel_path: "Info.plist".to_string(),
            base: FileBase::Share,
            role: FileRole::Asset,
            size: bytes.len() as u64,
            blake3: file_blake3_hex(&plist_path)?,
            mode: 0o644,
        });
    }
    let icns_path = res_dir.join("Nyx.icns");
    if !icns_path.is_file() {
        let icns = build_icns(MARK_PNG)?;
        std::fs::write(&icns_path, &icns)?;
        created.push(ManifestFile {
            rel_path: "Nyx.icns".to_string(),
            base: FileBase::Resource,
            role: FileRole::Asset,
            size: icns.len() as u64,
            blake3: file_blake3_hex(&icns_path)?,
            mode: 0o644,
        });
    }
    Ok(created)
}

#[cfg(target_os = "macos")]
fn lsregister_bin() -> &'static str {
    "/System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister"
}

pub struct MacosBackend;

impl PlatformBackend for MacosBackend {
    fn os(&self) -> TargetOs {
        TargetOs::Macos
    }

    fn default_install_root(&self, scope: InstallScope, paths: &HostPaths) -> PathBuf {
        default_root(scope, paths)
    }

    fn resolve_layout(&self, root: &Path, scope: InstallScope, paths: &HostPaths) -> LayoutPaths {
        compute_standard_layout(TargetOs::Macos, root, scope, paths)
    }

    fn check_elevation(&self, scope: InstallScope, allow_root_user: bool) -> Result<()> {
        #[cfg(unix)]
        {
            let euid = unsafe { libc::geteuid() };
            match scope {
                InstallScope::System if euid != 0 => {
                    bail!("administrator/root privileges are required for /Applications --system install (run with sudo or omit --system for ~/Applications install)");
                }
                InstallScope::User if euid == 0 && !allow_root_user => {
                    bail!("refusing --user install while running as root (pass --system or --allow-root-user)");
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
        match opts.scope {
            InstallScope::System => {
                // `/usr/local/bin` is on `/etc/paths`; if symlinked outside `/usr/local/bin` or custom `--prefix`, register `/etc/paths.d/nyx`.
                let _ = bin_dir;
                Ok((edits, vec![]))
            }
            InstallScope::User => {
                let Some(home) = &paths.home else {
                    return Ok((edits, vec![]));
                };
                let snippet = posix_path_block("$HOME/.local/bin");
                for (profile_path, is_fish, create) in [
                    (home.join(".zprofile"), false, false),
                    (home.join(".zshrc"), false, true),
                    (home.join(".bash_profile"), false, false),
                    (home.join(".bashrc"), false, false),
                    (home.join(".config/fish/config.fish"), true, false),
                ] {
                    let body = if is_fish {
                        fish_path_block("$HOME/.local/bin")
                    } else {
                        snippet.clone()
                    };
                    let (changed, created) =
                        upsert_marked_block(&profile_path, MARKER_ID, &body, create, opts.dry_run)?;
                    if changed {
                        edits.push(PathEdit::UnixShellProfile {
                            profile_path: profile_path.display().to_string(),
                            marker_block_id: MARKER_ID.to_string(),
                            created_file: created,
                        });
                    }
                }
                Ok((edits, vec![]))
            }
            InstallScope::Custom => Ok((edits, vec![])),
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
                PathEdit::MacosPathsD { file_path }
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
        _paths: &HostPaths,
    ) -> Result<(Vec<DesktopEntry>, Vec<String>)> {
        if opts.dry_run || !opts.desktop {
            return Ok((vec![], vec![]));
        }
        Ok((
            vec![DesktopEntry::MacosAppBundle {
                path: manifest.install_root.clone(),
            }],
            vec![],
        ))
    }

    fn remove_desktop(&self, _entries: &[DesktopEntry], _opts: &PlatformOpts) -> Result<bool> {
        Ok(true)
    }

    fn install_file_assoc(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        _paths: &HostPaths,
    ) -> Result<(Vec<FileAssociation>, Vec<String>)> {
        if !opts.file_assoc || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        let app = Path::new(&manifest.install_root);
        #[cfg(target_os = "macos")]
        {
            if Path::new(lsregister_bin()).is_file() {
                let _ = Command::new(lsregister_bin())
                    .args(["-f", &app.to_string_lossy()])
                    .output();
            }
        }
        Ok((
            vec![FileAssociation::MacosLaunchServices {
                app_path: app.display().to_string(),
            }],
            vec![],
        ))
    }

    fn remove_file_assoc(&self, _assocs: &[FileAssociation], _opts: &PlatformOpts) -> Result<bool> {
        Ok(true)
    }

    fn finalize_bundle(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        _paths: &HostPaths,
    ) -> Result<(Vec<ManifestFile>, Vec<String>)> {
        if opts.dry_run {
            return Ok((vec![], vec![]));
        }
        #[allow(unused_mut)]
        let mut warnings = Vec::new();
        let app = Path::new(&manifest.install_root);
        let added = ensure_bundle_metadata(app, &manifest.nyx_version)?;
        #[cfg(target_os = "macos")]
        {
            if !app.join("Contents/_CodeSignature").exists() {
                // Apple Silicon kills unsigned copied Mach-O binaries; apply inside-out ad-hoc codesign cleanly.
                let macos_dir = app.join("Contents/MacOS");
                if let Ok(rd) = std::fs::read_dir(&macos_dir) {
                    for e in rd.filter_map(Result::ok) {
                        let _ = Command::new("codesign")
                            .args(["--force", "--sign", "-"])
                            .arg(e.path())
                            .output();
                    }
                }
                let status = Command::new("codesign")
                    .args(["--force", "--sign", "-"])
                    .arg(app)
                    .output();
                if status.is_err() {
                    warnings.push(format!("ad-hoc codesign failed for {}", app.display()));
                }
            }
        }
        Ok((added, warnings))
    }

    fn prepare_uninstall_hook(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
    ) -> Result<()> {
        if opts.dry_run {
            return Ok(());
        }
        #[cfg(target_os = "macos")]
        {
            let app = Path::new(&manifest.install_root);
            if Path::new(lsregister_bin()).is_file() {
                let _ = Command::new(lsregister_bin())
                    .args(["-u", &app.to_string_lossy()])
                    .output();
            }
        }
        let _ = manifest;
        Ok(())
    }

    fn doctor_checks(
        &self,
        manifest: Option<&InstallManifest>,
        _paths: &HostPaths,
    ) -> Vec<DoctorCheck> {
        let mut checks = Vec::new();
        if let Some(m) = manifest {
            let app = Path::new(&m.install_root);
            let signed = app.join("Contents/_CodeSignature").exists();
            checks.push(DoctorCheck {
                name: "macos_bundle_seal".into(),
                tier: "optional".into(),
                status: if signed { "ok" } else { "info" }.into(),
                detail: if signed {
                    format!("code signature present at {}", app.display())
                } else {
                    format!(
                        "unsigned or ad-hoc bundle {}; if quarantined run `xattr -dr com.apple.quarantine \"{}\"`",
                        app.display(),
                        app.display()
                    )
                },
            });
            #[cfg(target_os = "macos")]
            {
                let ver = Command::new("codesign")
                    .args(["--verify", "--deep", "--strict"])
                    .arg(app)
                    .output();
                if let Ok(o) = ver {
                    if !o.status.success() {
                        checks.push(DoctorCheck {
                            name: "codesign_verify".into(),
                            tier: "required".into(),
                            status: "fail".into(),
                            detail: String::from_utf8_lossy(&o.stderr).trim().to_string(),
                        });
                    }
                }
            }
        }
        checks
    }
}
