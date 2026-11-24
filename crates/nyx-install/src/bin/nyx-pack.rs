// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Release packager CLI (`nyx-pack`): generates `THIRD_PARTY_NOTICES.txt`, packs streaming `NYXPK001` archives,
//! appends 8-byte-aligned 96B `NYXINST1` trailers onto `nyx-setup`, and emits `.zip`, dir trees, `.icns`, and `.dmg` stagings.

use anyhow::{bail, Context, Result};
use nyx_install::bundle::{
    binary_names, collect_pack_spec, exe_filename, target_triple_guess, write_dir_tree,
    write_nyxpack, write_self_extracting, write_zip_bundle, PackCollectOpts,
};
use nyx_install::manifest::TargetOs;
use nyx_install::platform::macos::{build_icns, render_info_plist};
use std::path::{Path, PathBuf};

const USAGE_TEXT: &str = "usage: nyx-pack --target-dir DIR --repo-root DIR [OPTIONS]\n\
     Options:\n\
       --out-dir DIR          Output directory (optional if explicit output paths given)\n\
       --stub PATH            Bare nyx-setup executable (bundled as nyx-uninstall and used for --self-extract)\n\
       --self-extract PATH    Write self-extracting nyx-setup executable with 96B NYXINST1 trailer\n\
       --bundle PATH          Write standalone .nyxpack payload\n\
       --zip PATH             Write standalone .zip distribution archive\n\
       --dir PATH             Stage unpacked directory layout\n\
       --dmg-staging DIR      Stage dual-path macOS .dmg root (Nyx.app + Set Up CLI & AI Hosts.app + Applications)\n\
       --notices PATH         Use pre-generated THIRD_PARTY_NOTICES.txt\n\
       --target TRIPLE        Cargo target triple override\n\
       --version VER          nyx version override\n\
       --omit-worker          Omit standalone nyx-worker from packed binaries\n\
       -V, --version          Print version and exit\n\
       -h, --help             Print this help message and exit";

fn usage(exit_code: i32) -> ! {
    if exit_code == 0 {
        println!("{USAGE_TEXT}");
    } else {
        eprintln!("{USAGE_TEXT}");
    }
    std::process::exit(exit_code);
}

fn copy_dir_recursive(src: &Path, dst: &Path) -> Result<()> {
    std::fs::create_dir_all(dst)?;
    for entry in std::fs::read_dir(src)? {
        let entry = entry?;
        let sp = entry.path();
        let dp = dst.join(entry.file_name());
        if sp.is_dir() {
            copy_dir_recursive(&sp, &dp)?;
        } else {
            std::fs::copy(&sp, &dp)?;
        }
    }
    Ok(())
}

fn write_blake3_sidecar(path: &Path) -> Result<()> {
    let hex = nyx_install::safety::file_blake3_hex(path)?;
    let bare = hex.strip_prefix("b3:").unwrap_or(&hex);
    let fname = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("artifact");
    let mut sidecar = path.as_os_str().to_os_string();
    sidecar.push(".blake3");
    std::fs::write(PathBuf::from(sidecar), format!("{bare}  {fname}\n"))?;
    Ok(())
}

fn stage_macos_dmg_root(
    staging: &Path,
    target_dir: &Path,
    repo_root: &Path,
    stub: &Path,
    notices_txt: &str,
    version: &str,
    include_worker: bool,
) -> Result<()> {
    std::fs::create_dir_all(staging)?;
    let nyx_app = staging.join("Nyx.app");
    let macos_dir = nyx_app.join("Contents/MacOS");
    let res_dir = nyx_app.join("Contents/Resources");
    std::fs::create_dir_all(&macos_dir)?;
    std::fs::create_dir_all(&res_dir)?;

    for stem in binary_names(include_worker) {
        let src = target_dir.join(stem);
        let dst = macos_dir.join(stem);
        std::fs::copy(&src, &dst).with_context(|| format!("copying {}", src.display()))?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let _ = std::fs::set_permissions(&dst, std::fs::Permissions::from_mode(0o755));
        }
    }
    std::fs::copy(stub, macos_dir.join("nyx-uninstall"))?;

    copy_dir_recursive(&repo_root.join("materials"), &res_dir.join("materials"))?;
    copy_dir_recursive(&repo_root.join("examples"), &res_dir.join("examples"))?;
    std::fs::copy(repo_root.join("LICENSE"), res_dir.join("LICENSE"))?;
    std::fs::copy(repo_root.join("NOTICE"), res_dir.join("NOTICE"))?;
    std::fs::write(res_dir.join("THIRD_PARTY_NOTICES.txt"), notices_txt)?;
    let icns = build_icns(&std::fs::read(repo_root.join("docs/brand/nyx-mark.png"))?)?;
    std::fs::write(res_dir.join("Nyx.icns"), &icns)?;
    std::fs::write(
        nyx_app.join("Contents/Info.plist"),
        render_info_plist(version, "nyx", "nyx", "dev.pjr.nyx", true),
    )?;

    // Companion helper bundle in .dmg root that auto-detects sibling `../Nyx.app`.
    let helper_app = staging.join("Set Up CLI & AI Hosts.app");
    std::fs::create_dir_all(helper_app.join("Contents/MacOS"))?;
    std::fs::create_dir_all(helper_app.join("Contents/Resources"))?;
    std::fs::copy(stub, helper_app.join("Contents/MacOS/nyx-setup"))?;
    std::fs::write(helper_app.join("Contents/Resources/Nyx.icns"), &icns)?;
    std::fs::write(
        helper_app.join("Contents/Info.plist"),
        render_info_plist(
            version,
            "nyx-setup",
            "Set Up CLI & AI Hosts",
            "dev.pjr.nyx.setup",
            false,
        ),
    )?;

    #[cfg(unix)]
    {
        let apps_link = staging.join("Applications");
        let _ = std::fs::remove_file(&apps_link);
        std::os::unix::fs::symlink("/Applications", apps_link)?;
    }
    Ok(())
}

fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        usage(2);
    }
    if matches!(args[0].as_str(), "-h" | "--help") {
        usage(0);
    }
    if matches!(args[0].as_str(), "-V" | "--version") {
        println!("nyx-pack {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    let mut target_dir: Option<PathBuf> = None;
    let mut repo_root: Option<PathBuf> = None;
    let mut out_dir: Option<PathBuf> = None;
    let mut stub: Option<PathBuf> = None;
    let mut self_extract: Option<PathBuf> = None;
    let mut bundle: Option<PathBuf> = None;
    let mut zip_out: Option<PathBuf> = None;
    let mut dir_out: Option<PathBuf> = None;
    let mut dmg_staging: Option<PathBuf> = None;
    let mut notices_path: Option<PathBuf> = None;
    let mut target_triple = target_triple_guess().to_string();
    let mut version = env!("CARGO_PKG_VERSION").to_string();
    let mut include_worker = true;

    let mut i = 0;
    while i < args.len() {
        let next = |i: &mut usize, flag: &str| -> Result<String> {
            *i += 1;
            args.get(*i)
                .cloned()
                .with_context(|| format!("{flag} requires a value"))
        };
        match args[i].as_str() {
            "--target-dir" => target_dir = Some(PathBuf::from(next(&mut i, "--target-dir")?)),
            "--repo-root" => repo_root = Some(PathBuf::from(next(&mut i, "--repo-root")?)),
            "--out-dir" => out_dir = Some(PathBuf::from(next(&mut i, "--out-dir")?)),
            "--stub" => stub = Some(PathBuf::from(next(&mut i, "--stub")?)),
            "--self-extract" => self_extract = Some(PathBuf::from(next(&mut i, "--self-extract")?)),
            "--bundle" => bundle = Some(PathBuf::from(next(&mut i, "--bundle")?)),
            "--zip" => zip_out = Some(PathBuf::from(next(&mut i, "--zip")?)),
            "--dir" => dir_out = Some(PathBuf::from(next(&mut i, "--dir")?)),
            "--dmg-staging" => dmg_staging = Some(PathBuf::from(next(&mut i, "--dmg-staging")?)),
            "--notices" => notices_path = Some(PathBuf::from(next(&mut i, "--notices")?)),
            "--target" => target_triple = next(&mut i, "--target")?,
            "--version" => version = next(&mut i, "--version")?,
            "--omit-worker" => include_worker = false,
            other => bail!("unknown argument '{other}'"),
        }
        i += 1;
    }

    let target_dir = target_dir.context("--target-dir is required")?;
    let repo_root = repo_root.context("--repo-root is required")?;
    let os = if target_triple.contains("windows") {
        TargetOs::Windows
    } else if target_triple.contains("apple") || target_triple.contains("darwin") {
        TargetOs::Macos
    } else {
        TargetOs::Linux
    };

    let default_stub = target_dir.join(exe_filename("nyx-setup", os));
    let stub_path = stub.unwrap_or(default_stub);
    let notices_txt = match notices_path {
        Some(p) => {
            std::fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?
        }
        None => {
            nyx_install::notices::generate_third_party_notices(&repo_root, Some(&target_triple))?
        }
    };

    let (spec, tmp_notices) = collect_pack_spec(&PackCollectOpts {
        bin_source_dir: &target_dir,
        repo_root: &repo_root,
        uninstaller_stub: &stub_path,
        notices_text: Some(&notices_txt),
        os,
        target_triple: &target_triple,
        nyx_version: &version,
        include_worker,
    })?;

    if let Some(d) = dir_out {
        write_dir_tree(&spec, &d)?;
        println!("dir: {}", d.display());
    }
    if let Some(z) = zip_out {
        write_zip_bundle(&spec, &z)?;
        write_blake3_sidecar(&z)?;
        println!("zip: {}", z.display());
    }
    if let Some(st) = dmg_staging {
        stage_macos_dmg_root(
            &st,
            &target_dir,
            &repo_root,
            &stub_path,
            &notices_txt,
            &version,
            include_worker,
        )?;
        println!("dmg_staging: {}", st.display());
    }

    if bundle.is_some() || self_extract.is_some() || out_dir.is_some() {
        let pack_path = bundle.clone().unwrap_or_else(|| {
            out_dir
                .as_ref()
                .map(|o| o.join(format!("nyx-{version}-{target_triple}.nyxpack")))
                .unwrap_or_else(|| {
                    std::env::temp_dir().join(format!("nyx-pack-{}.nyxpack", std::process::id()))
                })
        });
        let written = write_nyxpack(&spec, &pack_path)?;
        if bundle.is_some() || out_dir.is_some() {
            write_blake3_sidecar(&pack_path)?;
            println!("bundle: {}", pack_path.display());
        }
        if let Some(sfx) = self_extract.or_else(|| {
            out_dir.as_ref().map(|o| {
                let stem = format!("nyx-setup-{version}-{target_triple}");
                o.join(exe_filename(&stem, os))
            })
        }) {
            write_self_extracting(
                &stub_path,
                &pack_path,
                &sfx,
                &written.toc_bytes,
                written.toc_blake3,
                written.payload_blake3,
            )?;
            write_blake3_sidecar(&sfx)?;
            println!("self_extract: {}", sfx.display());
            if bundle.is_none() && out_dir.is_none() {
                let _ = std::fs::remove_file(&pack_path);
            }
        }
    }

    if let Some(tmp) = tmp_notices {
        let _ = std::fs::remove_file(tmp);
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("nyx-pack: {e}");
        std::process::exit(1);
    }
}
