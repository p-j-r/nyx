// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Windows native Registry (`REG_EXPAND_SZ`), ARP, file-association, and COM `.lnk` integration tests
//! scoped inside an isolated scratch key `HKCU\Software\nyx-install-test-<pid>` (cleaned via RAII `Drop`).
#![cfg(windows)]

use nyx_install::manifest::{
    DesktopEntry, FileAssociation, InstallLayout, InstallManifest, InstallScope, PathEdit, TargetOs,
};
use nyx_install::platform::windows::win32::{
    create_shortcut_lnk, delete_registry_tree, edit_path_in_subkey, install_file_associations,
    install_windows_desktop_and_arp, query_registry_value, remove_file_associations,
    remove_path_in_subkey, remove_windows_desktop_and_arp,
};
use nyx_install::platform::{HostPaths, PlatformOpts};
use std::fs;
use std::path::PathBuf;

struct RegSandbox {
    base: String,
    temp_dir: PathBuf,
}

impl RegSandbox {
    fn new(name: &str) -> Self {
        let pid = std::process::id();
        let base = format!(r"Software\nyx-install-test-{name}-{pid}");
        delete_registry_tree(false, &base);
        let temp_dir = std::env::temp_dir().join(format!("nyx-win-int-{name}-{pid}"));
        let _ = fs::remove_dir_all(&temp_dir);
        fs::create_dir_all(&temp_dir).unwrap();
        Self { base, temp_dir }
    }
}

impl Drop for RegSandbox {
    fn drop(&mut self) {
        delete_registry_tree(false, &self.base);
        let _ = fs::remove_dir_all(&self.temp_dir);
    }
}

#[test]
fn real_win32_registry_path_preserves_reg_expand_sz_and_arp_and_assoc_cleanup() {
    let sb = RegSandbox::new("reg");
    let env_sub = format!(r"{}\Environment", sb.base);

    // Pre-seed REG_EXPAND_SZ (`%SystemRoot%\System32`).
    let first = edit_path_in_subkey(
        false,
        &env_sub,
        "Path",
        r"%SystemRoot%\System32",
        None,
        false,
    )
    .unwrap()
    .unwrap();
    let _ = first;
    let added = edit_path_in_subkey(
        false,
        &env_sub,
        "Path",
        r"C:\Users\Alice\AppData\Local\Programs\nyx\bin",
        None,
        false,
    )
    .unwrap()
    .unwrap();
    let (ty, raw) = query_registry_value(false, &env_sub, "Path")
        .unwrap()
        .unwrap();
    assert_eq!(ty, 2, "must preserve REG_EXPAND_SZ (type 2)");
    assert!(
        raw.contains(r"%SystemRoot%\System32"),
        "must not expand %%SystemRoot%%: {raw}"
    );
    assert!(raw.contains(r"Programs\nyx\bin"));

    if let PathEdit::WindowsRegistry {
        entry_added,
        created_value,
        ..
    } = added
    {
        assert!(
            remove_path_in_subkey(false, &env_sub, "Path", &entry_added, created_value, false)
                .unwrap()
        );
    }
    let (ty2, raw2) = query_registry_value(false, &env_sub, "Path")
        .unwrap()
        .unwrap();
    assert_eq!(ty2, 2);
    assert_eq!(raw2, r"%SystemRoot%\System32");

    // Verify ARP and OpenWith file-association writes + clean removal inside scratch hive.
    let manifest = InstallManifest {
        format: 1,
        schema: "nyx-install.manifest".into(),
        nyx_version: "0.1.0".into(),
        installer_version: "0.1.0".into(),
        install_id: "test".into(),
        target_triple: "x86_64-pc-windows-msvc".into(),
        os: TargetOs::Windows,
        scope: InstallScope::User,
        layout: InstallLayout::Split,
        install_root: r"C:\Users\Alice\AppData\Local\Programs\nyx".into(),
        bin_dir: r"C:\Users\Alice\AppData\Local\Programs\nyx\bin".into(),
        resource_dir: r"C:\Users\Alice\AppData\Local\Programs\nyx".into(),
        manifest_path: r"C:\Users\Alice\AppData\Local\Programs\nyx\install-manifest.json".into(),
        installed_at: "2026-11-22T00:00:00Z".into(),
        updated_at: "2026-11-22T00:00:00Z".into(),
        previous_version: None,
        files: vec![],
        symlinks: vec![],
        path_edits: vec![],
        desktop_entries: vec![],
        file_associations: vec![],
        ai_hosts: vec![],
    };
    let mut opts = PlatformOpts::default_for(TargetOs::Windows, InstallScope::User);
    opts.testing_reg_root = Some(sb.base.clone());
    opts.desktop = false; // ARP tested via registry; COM .lnk tested separately below.

    let paths = HostPaths::sandboxed(&sb.temp_dir, TargetOs::Windows);
    let (desk, _) = install_windows_desktop_and_arp(&manifest, &opts, &paths).unwrap();
    assert!(matches!(
        &desk[0],
        DesktopEntry::WindowsUninstallRegistry { .. }
    ));
    let (_, icon_val) =
        query_registry_value(false, &format!(r"{}\Uninstall\nyx", sb.base), "DisplayIcon")
            .unwrap()
            .unwrap();
    assert!(
        !icon_val.starts_with('"') && icon_val.ends_with(",0"),
        "DisplayIcon must be unquoted: {icon_val}"
    );
    let (_, q_uninst) = query_registry_value(
        false,
        &format!(r"{}\Uninstall\nyx", sb.base),
        "QuietUninstallString",
    )
    .unwrap()
    .unwrap();
    assert!(
        q_uninst.starts_with('"'),
        "QuietUninstallString must quote exe: {q_uninst}"
    );
    assert!(remove_windows_desktop_and_arp(&desk).unwrap());
    assert!(
        query_registry_value(false, &format!(r"{}\Uninstall\nyx", sb.base), "DisplayName")
            .unwrap()
            .is_none()
    );

    let (assocs, _) = install_file_associations(&manifest, &opts).unwrap();
    assert!(assocs
        .iter()
        .any(|a| matches!(a, FileAssociation::WindowsRegistryKey { .. })));
    assert!(
        query_registry_value(false, &format!(r"{}\Classes\.nyx", sb.base), "")
            .unwrap()
            .is_some_and(|(_, v)| v == "nyx.Archive.1")
    );
    assert!(
        query_registry_value(false, &format!(r"{}\Classes\.json", sb.base), "")
            .unwrap()
            .is_none(),
        ".json Default must NEVER be touched"
    );
    assert!(remove_file_associations(&assocs).unwrap());
    assert!(
        query_registry_value(false, &format!(r"{}\Classes\nyx.Archive.1", sb.base), "")
            .unwrap()
            .is_none()
    );
}

#[test]
fn real_win32_com_creates_valid_lnk_shortcut_file() {
    let sb = RegSandbox::new("lnk");
    let fake_exe = sb.temp_dir.join("nyx.exe");
    fs::write(&fake_exe, b"MZ").unwrap();
    let lnk = sb.temp_dir.join("nyx.lnk");
    create_shortcut_lnk(&lnk, &fake_exe, &sb.temp_dir, "nyx test shortcut").unwrap();
    let bytes = fs::read(&lnk).unwrap();
    assert!(bytes.len() > 76, "lnk too small");
    assert_eq!(
        &bytes[0..4],
        &[0x4C, 0x00, 0x00, 0x00],
        "missing ShellLink LNK header magic"
    );
}
