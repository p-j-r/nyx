// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Windows platform driver: `HKCU`/`HKLM` `Path` (`REG_EXPAND_SZ`), `WM_SETTINGCHANGE`,
//! Add/Remove Programs (`Uninstall\nyx`), COM `.lnk` shortcuts, safe `OpenWithProgids`,
//! `SetDefaultDllDirectories`, `AttachConsole`, and pure-Rust Mode A / Mode B self-deletion.

use super::{compute_standard_layout, HostPaths, LayoutPaths, PlatformBackend, PlatformOpts};
use crate::manifest::{
    DesktopEntry, FileAssociation, InstallManifest, InstallScope, ManifestSymlink, PathEdit,
    TargetOs,
};
use crate::report::DoctorCheck;
#[cfg(windows)]
use anyhow::bail;
use anyhow::Result;
use std::path::{Path, PathBuf};

pub const ARP_SUBKEY: &str = "Software\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\nyx";
pub const AUMID: &str = "dev.pjr.nyx";

pub fn default_root(scope: InstallScope, paths: &HostPaths) -> PathBuf {
    match scope {
        InstallScope::System => paths
            .program_files
            .clone()
            .unwrap_or_else(|| PathBuf::from(r"C:\Program Files"))
            .join("nyx"),
        InstallScope::User | InstallScope::Custom => paths
            .local_appdata
            .clone()
            .or_else(|| paths.home.as_ref().map(|h| h.join("AppData/Local")))
            .unwrap_or_else(|| PathBuf::from(r"C:\Users\Default\AppData\Local"))
            .join("Programs")
            .join("nyx"),
    }
}

/// Normalizes a single Windows `Path` segment solely for case/slash/quote-insensitive comparison.
pub fn norm_win_path(seg: &str) -> String {
    let s = seg.trim().trim_matches('"').trim();
    let s = s.strip_prefix(r"\\?\").unwrap_or(s);
    let mut s = s.replace('/', "\\");
    while s.ends_with('\\') && !(s.len() == 3 && s.as_bytes().get(1) == Some(&b':')) && s.len() > 1
    {
        s.pop();
    }
    s.to_lowercase()
}

/// Appends `add_entry` (preserving all existing segments verbatim) and optionally strips `remove_old`.
/// Returns `(new_raw_value, changed, has_env_tokens)`.
pub fn edit_path_string(
    raw: &str,
    add_entry: &str,
    remove_old: Option<&str>,
) -> (String, bool, bool) {
    let target_norm = norm_win_path(add_entry);
    let old_norm = remove_old.map(norm_win_path);
    let mut parts: Vec<String> = Vec::new();
    let mut already = false;
    let mut removed_any = false;

    for raw_seg in raw.split(';') {
        let trimmed = raw_seg.trim();
        if trimmed.is_empty() {
            continue;
        }
        let seg_norm = norm_win_path(trimmed);
        if let Some(ref on) = old_norm {
            if &seg_norm == on && on != &target_norm {
                removed_any = true;
                continue;
            }
        }
        if seg_norm == target_norm {
            if already {
                removed_any = true;
                continue;
            }
            already = true;
        }
        parts.push(trimmed.to_string());
    }

    let mut changed = removed_any;
    if !already {
        parts.push(add_entry.trim().to_string());
        changed = true;
    }
    let joined = parts.join(";");
    let has_tokens = joined.contains('%');
    (joined, changed, has_tokens)
}

/// Removes every segment matching `entry` while preserving all remaining segments verbatim.
pub fn remove_from_path_string(raw: &str, entry: &str) -> (String, bool) {
    let target_norm = norm_win_path(entry);
    let mut parts: Vec<String> = Vec::new();
    let mut changed = false;
    for raw_seg in raw.split(';') {
        let trimmed = raw_seg.trim();
        if trimmed.is_empty() {
            continue;
        }
        if norm_win_path(trimmed) == target_norm {
            changed = true;
            continue;
        }
        parts.push(trimmed.to_string());
    }
    (parts.join(";"), changed)
}

/// Unquoted `<root>\bin\nyx.exe,0` (`DisplayIcon` is parsed by `PathParseIconLocationW`; quotes break Windows Settings icons).
pub fn display_icon_value(nyx_exe: &Path) -> String {
    format!(
        "{},0",
        nyx_connect::plain_path(nyx_exe.to_path_buf()).display()
    )
}

/// Quoted command string passed to `CreateProcessW` (mitigates CWE-428 unquoted service/uninstall paths).
pub fn quote_command_exe(exe: &Path, args: &str) -> String {
    let clean = nyx_connect::plain_path(exe.to_path_buf());
    if args.is_empty() {
        format!("\"{}\"", clean.display())
    } else {
        format!("\"{}\" {args}", clean.display())
    }
}

pub struct WindowsBackend;

impl PlatformBackend for WindowsBackend {
    fn os(&self) -> TargetOs {
        TargetOs::Windows
    }

    fn default_install_root(&self, scope: InstallScope, paths: &HostPaths) -> PathBuf {
        default_root(scope, paths)
    }

    fn resolve_layout(&self, root: &Path, scope: InstallScope, paths: &HostPaths) -> LayoutPaths {
        compute_standard_layout(TargetOs::Windows, root, scope, paths)
    }

    fn harden_process(&self) {
        #[cfg(windows)]
        win32::harden_dll_search_order();
    }

    fn check_elevation(&self, scope: InstallScope, _allow_root_user: bool) -> Result<()> {
        if scope == InstallScope::System {
            #[cfg(windows)]
            if !win32::is_process_elevated() {
                bail!("administrator privileges are required for --system installation (re-run elevated or omit --system for per-user install)");
            }
        }
        Ok(())
    }

    fn install_symlinks(
        &self,
        _bin_dir: &Path,
        _link_dir: Option<&Path>,
        _bin_names: &[String],
        _force: bool,
        _dry_run: bool,
    ) -> Result<Vec<ManifestSymlink>> {
        // Windows uses `<InstallRoot>\bin` directly on PATH, no per-binary symlinks needed.
        Ok(vec![])
    }

    fn remove_symlinks(&self, _symlinks: &[ManifestSymlink], _dry_run: bool) -> Result<usize> {
        Ok(0)
    }

    fn apply_path(
        &self,
        bin_dir: &Path,
        old_bin_dir: Option<&Path>,
        opts: &PlatformOpts,
        _paths: &HostPaths,
    ) -> Result<(Vec<PathEdit>, Vec<String>)> {
        if !opts.add_path || opts.dry_run {
            return Ok((vec![], vec![]));
        }
        #[cfg(windows)]
        {
            let clean_bin = nyx_connect::plain_path(bin_dir.to_path_buf());
            let clean_old = old_bin_dir.map(|p| nyx_connect::plain_path(p.to_path_buf()));
            let subkey = opts.testing_reg_root.as_deref().map_or_else(
                || {
                    if opts.scope == InstallScope::System {
                        "SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment"
                            .to_string()
                    } else {
                        "Environment".to_string()
                    }
                },
                |prefix| format!("{prefix}\\Environment"),
            );
            let use_hklm = opts.scope == InstallScope::System && opts.testing_reg_root.is_none();
            let edit = win32::edit_path_in_subkey(
                use_hklm,
                &subkey,
                "Path",
                &clean_bin.to_string_lossy(),
                clean_old.as_ref().map(|p| p.to_string_lossy()).as_deref(),
                opts.testing_reg_root.is_none(),
            )?;
            Ok((edit.into_iter().collect(), vec![]))
        }
        #[cfg(not(windows))]
        {
            let _ = (bin_dir, old_bin_dir);
            Ok((vec![], vec![]))
        }
    }

    fn remove_path(&self, edits: &[PathEdit], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        #[cfg(windows)]
        {
            let mut changed = false;
            for edit in edits {
                if let PathEdit::WindowsRegistry {
                    hive,
                    subkey,
                    value_name,
                    entry_added,
                    created_value,
                    ..
                } = edit
                {
                    let use_hklm = hive.eq_ignore_ascii_case("HKLM");
                    if win32::remove_path_in_subkey(
                        use_hklm,
                        subkey,
                        value_name,
                        entry_added,
                        *created_value,
                        opts.testing_reg_root.is_none(),
                    )? {
                        changed = true;
                    }
                }
            }
            Ok(changed)
        }
        #[cfg(not(windows))]
        {
            let _ = edits;
            Ok(false)
        }
    }

    fn install_desktop(
        &self,
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<DesktopEntry>, Vec<String>)> {
        if opts.dry_run {
            return Ok((vec![], vec![]));
        }
        #[cfg(windows)]
        {
            win32::install_windows_desktop_and_arp(manifest, opts, paths)
        }
        #[cfg(not(windows))]
        {
            let _ = (manifest, paths);
            Ok((vec![], vec![]))
        }
    }

    fn remove_desktop(&self, entries: &[DesktopEntry], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        #[cfg(windows)]
        {
            win32::remove_windows_desktop_and_arp(entries)
        }
        #[cfg(not(windows))]
        {
            let _ = entries;
            Ok(false)
        }
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
        #[cfg(windows)]
        {
            win32::install_file_associations(manifest, opts)
        }
        #[cfg(not(windows))]
        {
            let _ = manifest;
            Ok((vec![], vec![]))
        }
    }

    fn remove_file_assoc(&self, assocs: &[FileAssociation], opts: &PlatformOpts) -> Result<bool> {
        if opts.dry_run {
            return Ok(false);
        }
        #[cfg(windows)]
        {
            win32::remove_file_associations(assocs)
        }
        #[cfg(not(windows))]
        {
            let _ = assocs;
            Ok(false)
        }
    }

    fn doctor_checks(
        &self,
        manifest: Option<&InstallManifest>,
        _paths: &HostPaths,
    ) -> Vec<DoctorCheck> {
        let mut checks = Vec::new();
        let sys_root = std::env::var("SystemRoot").unwrap_or_else(|_| r"C:\Windows".to_string());
        let vcruntime = Path::new(&sys_root)
            .join("System32")
            .join("vcruntime140.dll");
        checks.push(DoctorCheck {
            name: "vcruntime140".into(),
            tier: "required".into(),
            status: if vcruntime.is_file() { "ok" } else { "warn" }.into(),
            detail: if vcruntime.is_file() {
                format!("found {}", vcruntime.display())
            } else {
                format!(
                    "{} not found (install Microsoft Visual C++ 2015-2022 Redistributable if nyx fails to start)",
                    vcruntime.display()
                )
            },
        });
        if let Some(m) = manifest {
            let bin_dir = Path::new(&m.bin_dir);
            let missing_bins: Vec<&str> = [
                "nyx.exe",
                "nyx-mcp.exe",
                "nyx-script-worker.exe",
                "nyx-connect.exe",
            ]
            .into_iter()
            .filter(|b| !bin_dir.join(b).is_file())
            .collect();
            checks.push(DoctorCheck {
                name: "binaries".into(),
                tier: "required".into(),
                status: if missing_bins.is_empty() {
                    "ok"
                } else {
                    "fail"
                }
                .into(),
                detail: if missing_bins.is_empty() {
                    format!("all runtime binaries present in {}", bin_dir.display())
                } else {
                    format!("missing binaries: {}", missing_bins.join(", "))
                },
            });
        }
        checks
    }
}

pub fn attach_parent_console_if_needed() {
    #[cfg(windows)]
    win32::attach_console_if_unredirected();
}

pub fn prepare_self_delete_sync_reaper(
    uninstaller_exe: &Path,
    install_root: &Path,
    wait_pid: u32,
) -> Result<()> {
    #[cfg(windows)]
    {
        win32::spawn_reap_final(uninstaller_exe, install_root, wait_pid)
    }
    #[cfg(not(windows))]
    {
        let _ = (uninstaller_exe, install_root, wait_pid);
        Ok(())
    }
}

pub fn run_reap_final(wait_pid: u32, install_root: &Path) -> Result<()> {
    #[cfg(windows)]
    {
        win32::execute_reap_final(wait_pid, install_root)
    }
    #[cfg(not(windows))]
    {
        let _ = (wait_pid, install_root);
        Ok(())
    }
}

pub fn is_inside_dir(exe: &Path, dir: &Path) -> bool {
    let e = norm_win_path(&exe.to_string_lossy());
    let d = norm_win_path(&dir.to_string_lossy());
    e == d || e.starts_with(&format!("{d}\\"))
}

/// Returns `true` only when `subkey` names an application-owned `nyx` registry subtree safe fo
/// recursive `RegDeleteTreeW` deletion (never empty, top-level hive roots, `Software`, or shared
/// Windows keys such as `Software\Microsoft` or `Software\Classes`).
pub fn is_safe_nyx_registry_subtree(subkey: &str) -> bool {
    if subkey.is_empty() || subkey.contains('\0') {
        return false;
    }
    let norm = subkey
        .replace('/', "\\")
        .trim_matches('\\')
        .to_ascii_lowercase();
    if norm.is_empty()
        || norm
            .split('\\')
            .any(|p| p.is_empty() || p == "." || p == "..")
    {
        return false;
    }
    let parts: Vec<&str> = norm.split('\\').collect();
    if parts.len() < 2 {
        return false;
    }
    if norm == ARP_SUBKEY.to_ascii_lowercase() {
        return true;
    }
    // Isolated test sandbox roots (`Software\nyx-install-test-...`) and subkeys beneath them.
    if parts[0] == "software" && parts[1].starts_with("nyx-install-test-") {
        return true;
    }
    // Testing or per-app ARP uninstall subtree (`<prefix>\Uninstall\nyx`).
    if norm.ends_with("\\uninstall\\nyx") {
        return true;
    }
    // Exclusively owned ProgIDs (`...\Classes\nyx.*`) or `...\Classes\Applications\nyx.exe`.
    if let Some(idx) = parts.iter().position(|&p| p == "classes") {
        if let Some(&child) = parts.get(idx + 1) {
            if child.starts_with("nyx.")
                || (child == "applications" && parts.get(idx + 2) == Some(&"nyx.exe"))
            {
                return true;
            }
        }
    }
    false
}

/// Guards `execute_reap_final` against operating on bare drive roots or top-level Windows OS
/// directories when scheduling `MoveFileExW(..., MOVEFILE_DELAY_UNTIL_REBOOT)` or directory removal.
pub fn is_safe_reap_root(install_root: &Path) -> bool {
    let s = norm_win_path(&install_root.to_string_lossy());
    if s.is_empty() || s.contains('\0') {
        return false;
    }
    let trimmed = s.trim_end_matches('\\');
    if trimmed.len() <= 2 || trimmed.ends_with(':') {
        return false;
    }
    let Some(leaf) = install_root.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    if leaf.is_empty() {
        return false;
    }
    !matches!(
        trimmed,
        r"c:\windows"
            | r"c:\program files"
            | r"c:\program files (x86)"
            | r"c:\users"
            | r"c:\programdata"
    ) && install_root.parent().is_some()
}

/// Evacuates a running installed `nyx-uninstall.exe` into `%TEMP%` and relaunches GUI uninstall (Mode A).
pub fn evacuate_to_temp_and_relaunch(current_exe: &Path, args: &[String]) -> Result<bool> {
    #[cfg(windows)]
    {
        win32::evacuate_to_temp_and_relaunch(current_exe, args)
    }
    #[cfg(not(windows))]
    {
        let _ = (current_exe, args);
        Ok(false)
    }
}

#[cfg(windows)]
pub mod win32 {
    use super::*;
    use std::ffi::c_void;
    use std::os::windows::ffi::OsStrExt;
    use std::ptr::{null, null_mut};

    const HKEY_CURRENT_USER: isize = -0x7FFF_FFFFisize; // 0x80000001
    const HKEY_LOCAL_MACHINE: isize = -0x7FFF_FFFEisize; // 0x80000002
    const KEY_READ: u32 = 0x20019;
    const KEY_WRITE: u32 = 0x20006;
    const KEY_WOW64_64KEY: u32 = 0x0100;
    const REG_SZ: u32 = 1;
    const REG_EXPAND_SZ: u32 = 2;
    const REG_DWORD: u32 = 4;
    const ERROR_SUCCESS: i32 = 0;
    const ERROR_FILE_NOT_FOUND: i32 = 2;
    const HWND_BROADCAST: isize = 0xFFFF;
    const WM_SETTINGCHANGE: u32 = 0x001A;
    const SMTO_ABORTIFHUNG: u32 = 0x0002;
    const SHCNE_ASSOCCHANGED: i32 = 0x0800_0000;
    const SHCNF_IDLIST: u32 = 0x0000;
    const LOAD_LIBRARY_SEARCH_SYSTEM32: u32 = 0x0000_0800;
    const STD_OUTPUT_HANDLE: u32 = (-11i32) as u32;
    const STD_ERROR_HANDLE: u32 = (-12i32) as u32;
    const ATTACH_PARENT_PROCESS: u32 = 0xFFFF_FFFF;
    const TOKEN_QUERY: u32 = 0x0008;
    const TOKEN_ELEVATION_CLASS: u32 = 20;
    const FILE_TYPE_UNKNOWN: u32 = 0;
    const MOVEFILE_DELAY_UNTIL_REBOOT: u32 = 0x0000_0004;

    #[repr(C)]
    struct Guid {
        data1: u32,
        data2: u16,
        data3: u16,
        data4: [u8; 8],
    }

    const CLSID_SHELL_LINK: Guid = Guid {
        data1: 0x00021401,
        data2: 0x0000,
        data3: 0x0000,
        data4: [0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46],
    };
    const IID_ISHELL_LINK_W: Guid = Guid {
        data1: 0x000214F9,
        data2: 0x0000,
        data3: 0x0000,
        data4: [0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46],
    };
    const IID_IPERSIST_FILE: Guid = Guid {
        data1: 0x0000010B,
        data2: 0x0000,
        data3: 0x0000,
        data4: [0xC0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46],
    };
    const IID_IPROPERTY_STORE: Guid = Guid {
        data1: 0x886D8EEB,
        data2: 0x8CF2,
        data3: 0x4446,
        data4: [0x8D, 0x02, 0xCD, 0xBA, 0x1D, 0xBD, 0xCF, 0x99],
    };

    #[repr(C)]
    struct PropertyKey {
        fmtid: Guid,
        pid: u32,
    }

    const PKEY_APP_USER_MODEL_ID: PropertyKey = PropertyKey {
        fmtid: Guid {
            data1: 0x9F4C2855,
            data2: 0x9F79,
            data3: 0x4B39,
            data4: [0xA8, 0xD0, 0xE1, 0xD4, 0x2D, 0xE1, 0xD5, 0xF3],
        },
        pid: 5,
    };

    #[repr(C)]
    struct PropVariant {
        vt: u16,
        w_reserved1: u16,
        w_reserved2: u16,
        w_reserved3: u16,
        ptr_val: *const u16,
        extra: usize,
    }

    #[repr(C)]
    struct StartupInfoW {
        cb: u32,
        reserved: *mut u16,
        desktop: *mut u16,
        title: *mut u16,
        x: u32,
        y: u32,
        x_size: u32,
        y_size: u32,
        x_chars: u32,
        y_chars: u32,
        fill: u32,
        flags: u32,
        show: u16,
        reserved2: u16,
        reserved2_ptr: *mut u8,
        std_in: *mut c_void,
        std_out: *mut c_void,
        std_err: *mut c_void,
    }

    #[repr(C)]
    struct ProcessInformation {
        process: *mut c_void,
        thread: *mut c_void,
        pid: u32,
        tid: u32,
    }

    // Compile-time ABI layout invariants matching the Windows SDK `#[repr(C)]` definitions.
    const _: () = {
        assert!(std::mem::size_of::<Guid>() == 16);
        assert!(std::mem::align_of::<Guid>() == 4);
        assert!(std::mem::size_of::<PropertyKey>() == 20);
        assert!(std::mem::align_of::<PropertyKey>() == 4);
        assert!(
            std::mem::size_of::<PropVariant>()
                == if cfg!(target_pointer_width = "64") {
                    24
                } else {
                    16
                }
        );
        assert!(std::mem::size_of::<StartupInfoW>() > 0);
        assert!(std::mem::size_of::<ProcessInformation>() > 0);
    };

    #[link(name = "kernel32")]
    extern "system" {
        fn SetDefaultDllDirectories(flags: u32) -> i32;
        fn GetStdHandle(n_std_handle: u32) -> *mut c_void;
        fn GetFileType(handle: *mut c_void) -> u32;
        fn AttachConsole(pid: u32) -> i32;
        fn GetCurrentProcess() -> *mut c_void;
        fn OpenProcess(access: u32, inherit: i32, pid: u32) -> *mut c_void;
        fn WaitForSingleObject(h: *mut c_void, ms: u32) -> u32;
        fn CloseHandle(h: *mut c_void) -> i32;
        fn MoveFileExW(src: *const u16, dst: *const u16, flags: u32) -> i32;
        fn CreateProcessW(
            application: *const u16,
            command_line: *mut u16,
            process_attributes: *const c_void,
            thread_attributes: *const c_void,
            inherit_handles: i32,
            flags: u32,
            environment: *const c_void,
            directory: *const u16,
            startup: *const StartupInfoW,
            info: *mut ProcessInformation,
        ) -> i32;
    }

    #[link(name = "advapi32")]
    extern "system" {
        fn OpenProcessToken(proc_h: *mut c_void, access: u32, token_out: *mut *mut c_void) -> i32;
        fn GetTokenInformation(
            token_h: *mut c_void,
            info_class: u32,
            info: *mut c_void,
            len: u32,
            ret_len: *mut u32,
        ) -> i32;
        fn RegOpenKeyExW(
            hkey: isize,
            subkey: *const u16,
            opts: u32,
            sam: u32,
            result: *mut isize,
        ) -> i32;
        fn RegCreateKeyExW(
            hkey: isize,
            subkey: *const u16,
            reserved: u32,
            class: *const u16,
            opts: u32,
            sam: u32,
            sec: *const c_void,
            result: *mut isize,
            disp: *mut u32,
        ) -> i32;
        fn RegQueryValueExW(
            hkey: isize,
            value_name: *const u16,
            reserved: *mut u32,
            ty: *mut u32,
            data: *mut u8,
            data_len: *mut u32,
        ) -> i32;
        fn RegSetValueExW(
            hkey: isize,
            value_name: *const u16,
            reserved: u32,
            ty: u32,
            data: *const u8,
            data_len: u32,
        ) -> i32;
        fn RegDeleteValueW(hkey: isize, value_name: *const u16) -> i32;
        fn RegDeleteTreeW(hkey: isize, subkey: *const u16) -> i32;
        fn RegDeleteKeyW(hkey: isize, subkey: *const u16) -> i32;
        fn RegQueryInfoKeyW(
            hkey: isize,
            class: *mut u16,
            class_len: *mut u32,
            reserved: *mut u32,
            subkeys: *mut u32,
            max_subkey: *mut u32,
            max_class: *mut u32,
            values: *mut u32,
            max_val_name: *mut u32,
            max_val_len: *mut u32,
            sec_desc: *mut u32,
            last_write: *mut u64,
        ) -> i32;
        fn RegCloseKey(hkey: isize) -> i32;
    }

    #[link(name = "user32")]
    extern "system" {
        fn SendMessageTimeoutW(
            hwnd: isize,
            msg: u32,
            wparam: usize,
            lparam: isize,
            flags: u32,
            timeout: u32,
            result: *mut usize,
        ) -> usize;
    }

    #[link(name = "ole32")]
    extern "system" {
        fn CoInitializeEx(reserved: *const c_void, coinit: u32) -> i32;
        fn CoUninitialize();
        fn CoCreateInstance(
            rclsid: *const Guid,
            outer: *mut c_void,
            clsctx: u32,
            riid: *const Guid,
            ppv: *mut *mut c_void,
        ) -> i32;
    }

    #[link(name = "shell32")]
    extern "system" {
        fn SHChangeNotify(event_id: i32, flags: u32, item1: *const c_void, item2: *const c_void);
    }

    fn wide(s: &str) -> Vec<u16> {
        s.encode_utf16().chain(Some(0)).collect()
    }

    fn wide_path(p: &Path) -> Vec<u16> {
        p.as_os_str().encode_wide().chain(Some(0)).collect()
    }

    pub fn harden_dll_search_order() {
        // SAFETY: `LOAD_LIBRARY_SEARCH_SYSTEM32` is a valid `DWORD` flag fo
        // `kernel32!SetDefaultDllDirectories`; takes no pointers and returns `BOOL`.
        unsafe {
            let _ = SetDefaultDllDirectories(LOAD_LIBRARY_SEARCH_SYSTEM32);
        }
    }

    pub fn attach_console_if_unredirected() {
        // SAFETY: `GetStdHandle` uses valid `STD_*_HANDLE` constants. `GetFileType` is only called
        // after checking `!h.is_null()` and `h != INVALID_HANDLE_VALUE` (`-1isize`).
        // `AttachConsole(ATTACH_PARENT_PROCESS)` uses the valid Win32 sentinel `0xFFFF_FFFF`.
        unsafe {
            let h_out = GetStdHandle(STD_OUTPUT_HANDLE);
            let h_err = GetStdHandle(STD_ERROR_HANDLE);
            let has_out = !h_out.is_null()
                && h_out != (-1isize as *mut c_void)
                && GetFileType(h_out) != FILE_TYPE_UNKNOWN;
            let has_err = !h_err.is_null()
                && h_err != (-1isize as *mut c_void)
                && GetFileType(h_err) != FILE_TYPE_UNKNOWN;
            if !has_out && !has_err {
                let _ = AttachConsole(ATTACH_PARENT_PROCESS);
            }
        }
    }

    pub fn is_process_elevated() -> bool {
        // SAFETY: `GetCurrentProcess()` returns a valid pseudo-handle. `token` is checked non-null
        // (and != `-1`) before passing to `GetTokenInformation` (with an aligned `u32` buffer and
        // exact byte size) and is unconditionally closed with `CloseHandle(token)` before return.
        unsafe {
            let mut token: *mut c_void = null_mut();
            if OpenProcessToken(GetCurrentProcess(), TOKEN_QUERY, &mut token) == 0
                || token.is_null()
                || token == (-1isize as *mut c_void)
            {
                return false;
            }
            let mut elevated: u32 = 0;
            let mut ret_len: u32 = 0;
            let ok = GetTokenInformation(
                token,
                TOKEN_ELEVATION_CLASS,
                (&mut elevated as *mut u32).cast(),
                std::mem::size_of::<u32>() as u32,
                &mut ret_len,
            );
            CloseHandle(token);
            ok != 0 && elevated != 0
        }
    }

    fn root_hkey(use_hklm: bool) -> isize {
        if use_hklm {
            HKEY_LOCAL_MACHINE
        } else {
            HKEY_CURRENT_USER
        }
    }

    pub fn broadcast_environment_change() {
        let w_env = wide("Environment");
        let mut res: usize = 0;
        // SAFETY: `w_env` is a NUL-terminated UTF-16 buffer alive for the entire synchronous
        // `SendMessageTimeoutW` call (`WM_SETTINGCHANGE` expects an `LPCWSTR` in `LPARAM`), and
        // `&mut res` is a valid aligned `PDWORD_PTR` output slot.
        unsafe {
            let _ = SendMessageTimeoutW(
                HWND_BROADCAST,
                WM_SETTINGCHANGE,
                0,
                w_env.as_ptr() as isize,
                SMTO_ABORTIFHUNG,
                5000,
                &mut res,
            );
        }
    }

    fn reg_query_string(hkey: isize, value_name: &str) -> Result<Option<(u32, String)>> {
        if hkey == 0 {
            bail!("RegQueryValueExW called with null HKEY");
        }
        let w_val = wide(value_name);
        let val_ptr = if value_name.is_empty() {
            null()
        } else {
            w_val.as_ptr()
        };
        let mut ty: u32 = 0;
        let mut bytes: u32 = 0;
        // SAFETY: `hkey` is an open registry key handle, `val_ptr` is null or NUL-terminated UTF-16
        // alive across the call, `&mut ty` and `&mut bytes` are valid aligned `u32` out-pointers,
        // and `lpData` is null for a size query.
        let rc = unsafe {
            RegQueryValueExW(hkey, val_ptr, null_mut(), &mut ty, null_mut(), &mut bytes)
        };
        if rc == ERROR_FILE_NOT_FOUND {
            return Ok(None);
        }
        if rc != ERROR_SUCCESS {
            bail!("RegQueryValueExW failed ({rc}) for '{value_name}'");
        }
        let n_u16 = (bytes as usize).div_ceil(2) + 1;
        let byte_cap = n_u16
            .checked_mul(2)
            .and_then(|b| u32::try_from(b).ok())
            .ok_or_else(|| anyhow::anyhow!("registry value size overflows u32 for '{value_name}'"))?;
        let mut buf = vec![0u16; n_u16];
        let mut actual_bytes = byte_cap;
        // SAFETY: `buf` is uniquely owned and has `byte_cap` writable bytes matching `actual_bytes`;
        // `val_ptr`, `&mut ty`, and `&mut actual_bytes` remain valid for the duration of the call.
        let rc2 = unsafe {
            RegQueryValueExW(
                hkey,
                val_ptr,
                null_mut(),
                &mut ty,
                buf.as_mut_ptr().cast(),
                &mut actual_bytes,
            )
        };
        if rc2 != ERROR_SUCCESS {
            bail!("RegQueryValueExW read failed ({rc2}) for '{value_name}'");
        }
        let u16_len = (actual_bytes as usize) / 2;
        let slice = &buf[..u16_len.min(buf.len())];
        let end = slice.iter().position(|&c| c == 0).unwrap_or(slice.len());
        Ok(Some((ty, String::from_utf16_lossy(&slice[..end]))))
    }

    fn reg_set_string(hkey: isize, value_name: &str, ty: u32, text: &str) -> Result<()> {
        if hkey == 0 {
            bail!("RegSetValueExW called with null HKEY");
        }
        let w_val = wide(value_name);
        let w_data = wide(text);
        let bytes = w_data
            .len()
            .checked_mul(2)
            .and_then(|b| u32::try_from(b).ok())
            .ok_or_else(|| anyhow::anyhow!("registry string exceeds u32 byte length"))?;
        // SAFETY: `hkey` is a valid open registry key, value name pointer is null or NUL-terminated
        // UTF-16, and `w_data` is NUL-terminated UTF-16 valid for `bytes` readable bytes.
        let rc = unsafe {
            RegSetValueExW(
                hkey,
                if value_name.is_empty() {
                    null()
                } else {
                    w_val.as_ptr()
                },
                0,
                ty,
                w_data.as_ptr().cast(),
                bytes,
            )
        };
        if rc != ERROR_SUCCESS {
            bail!("RegSetValueExW failed ({rc}) for '{value_name}'");
        }
        Ok(())
    }

    fn reg_set_dword(hkey: isize, value_name: &str, v: u32) -> Result<()> {
        if hkey == 0 {
            bail!("RegSetValueExW DWORD called with null HKEY");
        }
        let w_val = wide(value_name);
        let bytes = v.to_le_bytes();
        // SAFETY: `w_val` is NUL-terminated UTF-16 and `bytes` is a 4-byte stack array matching
        // `cbData = 4` for `REG_DWORD`.
        let rc = unsafe { RegSetValueExW(hkey, w_val.as_ptr(), 0, REG_DWORD, bytes.as_ptr(), 4) };
        if rc != ERROR_SUCCESS {
            bail!("RegSetValueExW DWORD failed ({rc}) for '{value_name}'");
        }
        Ok(())
    }

    pub fn edit_path_in_subkey(
        use_hklm: bool,
        subkey: &str,
        value_name: &str,
        add_entry: &str,
        old_entry: Option<&str>,
        broadcast: bool,
    ) -> Result<Option<PathEdit>> {
        if subkey.is_empty() {
            bail!("edit_path_in_subkey refused empty subkey");
        }
        let w_sub = wide(subkey);
        let mut hkey: isize = 0;
        let mut disp: u32 = 0;
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined `HKEY_*` constant, `w_sub` is a
        // non-empty NUL-terminated UTF-16 string, and `&mut hkey` / `&mut disp` are valid pointers.
        let rc = unsafe {
            RegCreateKeyExW(
                root_hkey(use_hklm),
                w_sub.as_ptr(),
                0,
                null(),
                0,
                KEY_READ | KEY_WRITE | KEY_WOW64_64KEY,
                null(),
                &mut hkey,
                &mut disp,
            )
        };
        if rc != ERROR_SUCCESS {
            bail!("RegCreateKeyExW failed ({rc}) for {subkey}");
        }

        let res = (|| -> Result<Option<PathEdit>> {
            let existing = reg_query_string(hkey, value_name)?;
            let created_value = existing.is_none();
            let (orig_type, orig_raw) = existing.unwrap_or((REG_EXPAND_SZ, String::new()));
            let (new_raw, changed, has_tokens) = edit_path_string(&orig_raw, add_entry, old_entry);
            let write_type = if orig_type == REG_EXPAND_SZ || has_tokens || created_value {
                REG_EXPAND_SZ
            } else {
                orig_type
            };
            if changed {
                reg_set_string(hkey, value_name, write_type, &new_raw)?;
            }
            Ok(Some(PathEdit::WindowsRegistry {
                hive: if use_hklm { "HKLM" } else { "HKCU" }.into(),
                subkey: subkey.to_string(),
                value_name: value_name.to_string(),
                entry_added: add_entry.to_string(),
                created_value,
                preserved_expand_sz: write_type == REG_EXPAND_SZ,
            }))
        })();
        // SAFETY: `hkey` was successfully opened by `RegCreateKeyExW` above and is closed once.
        unsafe {
            RegCloseKey(hkey);
        }
        if broadcast {
            broadcast_environment_change();
        }
        res
    }

    pub fn remove_path_in_subkey(
        use_hklm: bool,
        subkey: &str,
        value_name: &str,
        entry: &str,
        created_value: bool,
        broadcast: bool,
    ) -> Result<bool> {
        if subkey.is_empty() {
            return Ok(false);
        }
        let w_sub = wide(subkey);
        let mut hkey: isize = 0;
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is NUL-terminated
        // UTF-16, and `&mut hkey` is a valid out-pointer.
        let rc = unsafe {
            RegOpenKeyExW(
                root_hkey(use_hklm),
                w_sub.as_ptr(),
                0,
                KEY_READ | KEY_WRITE | KEY_WOW64_64KEY,
                &mut hkey,
            )
        };
        if rc == ERROR_FILE_NOT_FOUND {
            return Ok(false);
        }
        if rc != ERROR_SUCCESS {
            bail!("RegOpenKeyExW failed ({rc}) for {subkey}");
        }
        let res = (|| -> Result<bool> {
            let Some((orig_type, orig_raw)) = reg_query_string(hkey, value_name)? else {
                return Ok(false);
            };
            let (updated, changed) = remove_from_path_string(&orig_raw, entry);
            if !changed {
                return Ok(false);
            }
            if updated.trim().is_empty() && created_value {
                let w_val = wide(value_name);
                // SAFETY: `hkey` is open with `KEY_WRITE` and `w_val` is NUL-terminated UTF-16.
                unsafe {
                    RegDeleteValueW(
                        hkey,
                        if value_name.is_empty() {
                            null()
                        } else {
                            w_val.as_ptr()
                        },
                    );
                }
            } else {
                let write_type = if orig_type == REG_EXPAND_SZ || updated.contains('%') {
                    REG_EXPAND_SZ
                } else {
                    orig_type
                };
                reg_set_string(hkey, value_name, write_type, &updated)?;
            }
            Ok(true)
        })();
        // SAFETY: `hkey` was opened by `RegOpenKeyExW` above and is closed once.
        unsafe {
            RegCloseKey(hkey);
        }
        if broadcast && res.as_ref().is_ok_and(|&c| c) {
            broadcast_environment_change();
        }
        res
    }

    pub fn query_registry_value(
        use_hklm: bool,
        subkey: &str,
        value_name: &str,
    ) -> Result<Option<(u32, String)>> {
        if subkey.is_empty() {
            return Ok(None);
        }
        let w_sub = wide(subkey);
        let mut hkey: isize = 0;
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is NUL-terminated
        // UTF-16, and `&mut hkey` is a valid out-pointer.
        let rc = unsafe {
            RegOpenKeyExW(
                root_hkey(use_hklm),
                w_sub.as_ptr(),
                0,
                KEY_READ | KEY_WOW64_64KEY,
                &mut hkey,
            )
        };
        if rc == ERROR_FILE_NOT_FOUND {
            return Ok(None);
        }
        if rc != ERROR_SUCCESS {
            bail!("RegOpenKeyExW failed ({rc})");
        }
        let r = reg_query_string(hkey, value_name);
        // SAFETY: `hkey` was opened by `RegOpenKeyExW` above and is closed once.
        unsafe {
            RegCloseKey(hkey);
        }
        r
    }

    pub fn delete_registry_tree(use_hklm: bool, subkey: &str) {
        // Scope guard: refuse to pass empty or non-nyx registry subtrees to `RegDeleteTreeW`.
        if !is_safe_nyx_registry_subtree(subkey) {
            return;
        }
        let w_sub = wide(subkey);
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is non-null
        // NUL-terminated UTF-16, and `is_safe_nyx_registry_subtree` verified `subkey` is an
        // exclusive nyx ARP/ProgID/test-sandbox subtree.
        unsafe {
            let _ = RegDeleteTreeW(root_hkey(use_hklm), w_sub.as_ptr());
        }
    }

    // SAFETY precondition: `obj` must be a non-null, pointer-aligned live COM interface pointe
    // whose first word is its vtable pointer.
    unsafe fn vtable_ptr(obj: *mut c_void) -> *const *const c_void {
        assert!(!obj.is_null(), "COM interface pointer must be non-null");
        assert!(
            (obj as usize).is_multiple_of(std::mem::align_of::<*const *const c_void>()),
            "COM interface pointer must be pointer-aligned"
        );
        // SAFETY: Caller guarantees `obj` points to a valid COM object with a leading vptr.
        let vt = unsafe { *(obj as *const *const *const c_void) };
        assert!(!vt.is_null(), "COM vtable pointer must be non-null");
        assert!(
            (vt as usize).is_multiple_of(std::mem::align_of::<*const c_void>()),
            "COM vtable pointer must be pointer-aligned"
        );
        vt
    }

    // SAFETY precondition: `vt` is a valid COM vtable for the target interface and `slot` is within
    // that interface's documented vtable length, matching function pointer signature `F`.
    #[inline]
    unsafe fn com_method<F: Copy>(vt: *const *const c_void, slot: usize) -> F {
        const {
            assert!(
                std::mem::size_of::<F>() == std::mem::size_of::<*const c_void>(),
                "COM vtable function pointer size mismatch"
            );
        }
        // SAFETY: `vt` is non-null/aligned and caller guarantees `slot` is in-bounds for the COM
        // interface; we also assert the loaded slot pointer is non-null before transmuting.
        let entry = unsafe { *vt.add(slot) };
        assert!(!entry.is_null(), "COM vtable slot {slot} is null");
        // SAFETY: `F` has pointer size (`size_of::<F>() == size_of::<*const c_void>()`), `entry`
        // is non-null, and caller guarantees `extern "system"` COM signature parity.
        unsafe { std::mem::transmute_copy::<*const c_void, F>(&entry) }
    }

    /// Creates a Windows `.lnk` Shell Link using pure-Rust COM `IShellLinkW`, `IPropertyStore` (AUMID), and `IPersistFile`.
    pub fn create_shortcut_lnk(
        lnk_path: &Path,
        target_exe: &Path,
        working_dir: &Path,
        description: &str,
    ) -> Result<()> {
        if let Some(p) = lnk_path.parent() {
            std::fs::create_dir_all(p)?;
        }
        let clean_target = nyx_connect::plain_path(target_exe.to_path_buf());
        let clean_wd = nyx_connect::plain_path(working_dir.to_path_buf());
        let clean_lnk = nyx_connect::plain_path(lnk_path.to_path_buf());

        // SAFETY: All COM calls run on the current thread inside a balanced `CoInitializeEx` /
        // `CoUninitialize` scope. Every returned COM pointer (`shell_link`, `prop_store`, `persist`)
        // is checked non-null and released via `IUnknown::Release` before `CoUninitialize`; all
        // `wide*` UTF-16 buffers outlive the COM method calls that borrow them, and vtable slots
        // match the documented Windows SDK `IUnknown` / `IShellLinkW` / `IPropertyStore` /
        // `IPersistFile` layouts (`shobjidl_core.h`, `propsys.h`, `objidl.h`).
        unsafe {
            let hr_init = CoInitializeEx(null(), 2 /* COINIT_APARTMENTTHREADED */);
            let need_uninit = hr_init == 0 || hr_init == 1;
            let res = (|| -> Result<()> {
                let mut shell_link: *mut c_void = null_mut();
                let hr = CoCreateInstance(
                    &CLSID_SHELL_LINK,
                    null_mut(),
                    1, /* CLSCTX_INPROC_SERVER */
                    &IID_ISHELL_LINK_W,
                    &mut shell_link,
                );
                if hr < 0 || shell_link.is_null() {
                    bail!("CoCreateInstance(ShellLink) failed: 0x{:08X}", hr as u32);
                }

                let vt = vtable_ptr(shell_link);
                // `IUnknown`: 0=QueryInterface, 2=Release; `IShellLinkW`: 7=SetDescription,
                // 9=SetWorkingDirectory, 17=SetIconLocation, 20=SetPath.
                let query_iface: unsafe extern "system" fn(
                    *mut c_void,
                    *const Guid,
                    *mut *mut c_void,
                ) -> i32 = com_method(vt, 0);
                let release_sl: unsafe extern "system" fn(*mut c_void) -> u32 = com_method(vt, 2);
                let set_desc: unsafe extern "system" fn(*mut c_void, *const u16) -> i32 =
                    com_method(vt, 7);
                let set_wd: unsafe extern "system" fn(*mut c_void, *const u16) -> i32 =
                    com_method(vt, 9);
                let set_icon: unsafe extern "system" fn(*mut c_void, *const u16, i32) -> i32 =
                    com_method(vt, 17);
                let set_path: unsafe extern "system" fn(*mut c_void, *const u16) -> i32 =
                    com_method(vt, 20);

                let w_target = wide_path(&clean_target);
                let w_wd = wide_path(&clean_wd);
                let w_desc = wide(description);
                let w_lnk = wide_path(&clean_lnk);

                let work = (|| -> Result<()> {
                    if set_path(shell_link, w_target.as_ptr()) < 0 {
                        bail!("IShellLinkW::SetPath failed");
                    }
                    if set_wd(shell_link, w_wd.as_ptr()) < 0 {
                        bail!("IShellLinkW::SetWorkingDirectory failed");
                    }
                    if set_icon(shell_link, w_target.as_ptr(), 0) < 0 {
                        bail!("IShellLinkW::SetIconLocation failed");
                    }
                    if set_desc(shell_link, w_desc.as_ptr()) < 0 {
                        bail!("IShellLinkW::SetDescription failed");
                    }

                    // Set PKEY_AppUserModel_ID so taskbar grouping matches the running window.
                    let mut prop_store: *mut c_void = null_mut();
                    if query_iface(shell_link, &IID_IPROPERTY_STORE, &mut prop_store) >= 0
                        && !prop_store.is_null()
                    {
                        let pvt = vtable_ptr(prop_store);
                        // `IPropertyStore`: 2=Release, 6=SetValue, 7=Commit.
                        let rel_ps: unsafe extern "system" fn(*mut c_void) -> u32 =
                            com_method(pvt, 2);
                        let set_val: unsafe extern "system" fn(
                            *mut c_void,
                            *const PropertyKey,
                            *const PropVariant,
                        ) -> i32 = com_method(pvt, 6);
                        let commit: unsafe extern "system" fn(*mut c_void) -> i32 =
                            com_method(pvt, 7);
                        let w_aumid = wide(AUMID);
                        let pv = PropVariant {
                            vt: 31, // VT_LPWSTR
                            w_reserved1: 0,
                            w_reserved2: 0,
                            w_reserved3: 0,
                            ptr_val: w_aumid.as_ptr(),
                            extra: 0,
                        };
                        let _ = set_val(prop_store, &PKEY_APP_USER_MODEL_ID, &pv);
                        let _ = commit(prop_store);
                        rel_ps(prop_store);
                    }

                    let mut persist: *mut c_void = null_mut();
                    let hr_pf = query_iface(shell_link, &IID_IPERSIST_FILE, &mut persist);
                    if hr_pf < 0 || persist.is_null() {
                        bail!(
                            "QueryInterface(IPersistFile) failed: 0x{:08X}",
                            hr_pf as u32
                        );
                    }
                    let fvt = vtable_ptr(persist);
                    // `IPersistFile`: 2=Release, 6=Save.
                    let rel_pf: unsafe extern "system" fn(*mut c_void) -> u32 = com_method(fvt, 2);
                    let save_pf: unsafe extern "system" fn(*mut c_void, *const u16, i32) -> i32 =
                        com_method(fvt, 6);
                    let hr_save = save_pf(persist, w_lnk.as_ptr(), 1);
                    rel_pf(persist);
                    if hr_save < 0 {
                        bail!("IPersistFile::Save failed: 0x{:08X}", hr_save as u32);
                    }
                    Ok(())
                })();
                release_sl(shell_link);
                work
            })();
            if need_uninit {
                CoUninitialize();
            }
            res
        }
    }

    pub fn install_windows_desktop_and_arp(
        manifest: &InstallManifest,
        opts: &PlatformOpts,
        paths: &HostPaths,
    ) -> Result<(Vec<DesktopEntry>, Vec<String>)> {
        let mut entries = Vec::new();
        let warnings = Vec::new();
        let use_hklm = opts.scope == InstallScope::System && opts.testing_reg_root.is_none();
        let hive_name = if use_hklm { "HKLM" } else { "HKCU" }.to_string();
        let subkey = opts.testing_reg_root.as_ref().map_or_else(
            || ARP_SUBKEY.to_string(),
            |p| format!("{p}\\Uninstall\\nyx"),
        );

        let bin_dir = Path::new(&manifest.bin_dir);
        let nyx_exe = bin_dir.join("nyx.exe");
        let uninst_exe = bin_dir.join("nyx-uninstall.exe");
        let total_bytes: u64 = manifest.files.iter().map(|f| f.size).sum();
        let est_kib = total_bytes.div_ceil(1024).min(u32::MAX as u64) as u32;
        let date_str = manifest
            .updated_at
            .get(..10)
            .unwrap_or("2026-01-01")
            .replace('-', "");

        let w_sub = wide(&subkey);
        let mut hkey: isize = 0;
        let mut disp: u32 = 0;
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is NUL-terminated
        // UTF-16, and `&mut hkey` / `&mut disp` are valid out-pointers.
        let rc = unsafe {
            RegCreateKeyExW(
                root_hkey(use_hklm),
                w_sub.as_ptr(),
                0,
                null(),
                0,
                KEY_READ | KEY_WRITE | KEY_WOW64_64KEY,
                null(),
                &mut hkey,
                &mut disp,
            )
        };
        if rc != ERROR_SUCCESS {
            bail!("failed to write Uninstall registry key {subkey} ({rc})");
        }
        let r = (|| -> Result<()> {
            reg_set_string(hkey, "DisplayName", REG_SZ, "nyx")?;
            reg_set_string(hkey, "DisplayVersion", REG_SZ, &manifest.nyx_version)?;
            reg_set_string(hkey, "Publisher", REG_SZ, "p-j-r")?;
            reg_set_string(
                hkey,
                "InstallLocation",
                REG_SZ,
                &nyx_connect::plain_path(PathBuf::from(&manifest.install_root)).to_string_lossy(),
            )?;
            reg_set_string(hkey, "InstallDate", REG_SZ, &date_str)?;
            reg_set_string(hkey, "DisplayIcon", REG_SZ, &display_icon_value(&nyx_exe))?;
            reg_set_string(
                hkey,
                "UninstallString",
                REG_SZ,
                &quote_command_exe(&uninst_exe, ""),
            )?;
            reg_set_string(
                hkey,
                "QuietUninstallString",
                REG_SZ,
                &quote_command_exe(&uninst_exe, "uninstall --yes --quiet"),
            )?;
            reg_set_string(
                hkey,
                "ModifyPath",
                REG_SZ,
                &quote_command_exe(&uninst_exe, "repair"),
            )?;
            reg_set_dword(hkey, "EstimatedSize", est_kib)?;
            reg_set_dword(hkey, "NoModify", 1)?;
            reg_set_dword(hkey, "NoRepair", 0)?;
            Ok(())
        })();
        // SAFETY: `hkey` was opened by `RegCreateKeyExW` above and is closed once.
        unsafe {
            RegCloseKey(hkey);
        }
        r?;
        entries.push(DesktopEntry::WindowsUninstallRegistry {
            hive: hive_name,
            subkey,
        });

        if opts.desktop && opts.testing_reg_root.is_none() {
            let working_dir = paths
                .home
                .as_ref()
                .map(|h| h.join("Documents"))
                .unwrap_or_else(|| PathBuf::from(&manifest.install_root));
            let start_programs = if use_hklm {
                paths
                    .program_data
                    .as_ref()
                    .map(|p| p.join("Microsoft\\Windows\\Start Menu\\Programs"))
            } else {
                paths
                    .appdata
                    .as_ref()
                    .map(|p| p.join("Microsoft\\Windows\\Start Menu\\Programs"))
            };
            if let Some(dir) = start_programs {
                let lnk = dir.join("nyx.lnk");
                create_shortcut_lnk(
                    &lnk,
                    &nyx_exe,
                    &working_dir,
                    "nyx — Finite-element photonics, phononics & electromagnetics",
                )?;
                entries.push(DesktopEntry::WindowsShortcut {
                    path: lnk.display().to_string(),
                });
            }
            if opts.desktop_shortcut {
                if let Some(desk) = &paths.desktop_dir {
                    let lnk = desk.join("nyx.lnk");
                    create_shortcut_lnk(
                        &lnk,
                        &nyx_exe,
                        &working_dir,
                        "nyx — Finite-element photonics, phononics & electromagnetics",
                    )?;
                    entries.push(DesktopEntry::WindowsShortcut {
                        path: lnk.display().to_string(),
                    });
                }
            }
        }
        Ok((entries, warnings))
    }

    pub fn remove_windows_desktop_and_arp(entries: &[DesktopEntry]) -> Result<bool> {
        let mut changed = false;
        for entry in entries {
            match entry {
                DesktopEntry::WindowsShortcut { path } => {
                    if Path::new(path).is_file() && std::fs::remove_file(path).is_ok() {
                        changed = true;
                    }
                }
                DesktopEntry::WindowsUninstallRegistry { hive, subkey } => {
                    let use_hklm = hive.eq_ignore_ascii_case("HKLM");
                    delete_registry_tree(use_hklm, subkey);
                    changed = true;
                }
                _ => {}
            }
        }
        Ok(changed)
    }

    fn write_reg_key_value(
        use_hklm: bool,
        subkey: &str,
        value_name: &str,
        value: &str,
        record_prev: bool,
    ) -> Result<(bool, Option<String>)> {
        if subkey.is_empty() {
            bail!("write_reg_key_value refused empty subkey");
        }
        let w_sub = wide(subkey);
        let mut hkey: isize = 0;
        let mut disp: u32 = 0;
        // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is NUL-terminated
        // UTF-16, and `&mut hkey` / `&mut disp` are valid out-pointers.
        let rc = unsafe {
            RegCreateKeyExW(
                root_hkey(use_hklm),
                w_sub.as_ptr(),
                0,
                null(),
                0,
                KEY_READ | KEY_WRITE | KEY_WOW64_64KEY,
                null(),
                &mut hkey,
                &mut disp,
            )
        };
        if rc != ERROR_SUCCESS {
            bail!("RegCreateKeyExW failed ({rc}) for {subkey}");
        }
        let created_key = disp == 1; // REG_CREATED_NEW_KEY
        let prev = if record_prev {
            reg_query_string(hkey, value_name)
                .ok()
                .flatten()
                .map(|(_, s)| s)
        } else {
            None
        };
        let res = reg_set_string(hkey, value_name, REG_SZ, value);
        // SAFETY: `hkey` was opened by `RegCreateKeyExW` above and is closed once.
        unsafe {
            RegCloseKey(hkey);
        }
        res?;
        Ok((created_key, prev))
    }

    pub fn install_file_associations(
        manifest: &InstallManifest,
        opts: &PlatformOpts,
    ) -> Result<(Vec<FileAssociation>, Vec<String>)> {
        let use_hklm = opts.scope == InstallScope::System && opts.testing_reg_root.is_none();
        let hive = if use_hklm { "HKLM" } else { "HKCU" }.to_string();
        let classes_base = opts.testing_reg_root.as_ref().map_or_else(
            || "Software\\Classes".to_string(),
            |p| format!("{p}\\Classes"),
        );
        let nyx_exe = Path::new(&manifest.bin_dir).join("nyx.exe");
        let open_cmd = quote_command_exe(&nyx_exe, "\"%1\"");
        let icon_val = display_icon_value(&nyx_exe);

        let mut records = Vec::new();
        for (prog_id, friendly) in [
            ("nyx.Archive.1", "nyx Workflow Archive"),
            ("nyx.Script.1", "nyx Simulation Script"),
            ("nyx.Project.1", "nyx Simulation Project"),
            ("nyx.Results.1", "nyx Results File"),
        ] {
            let root = format!("{classes_base}\\{prog_id}");
            write_reg_key_value(use_hklm, &root, "", friendly, false)?;
            write_reg_key_value(
                use_hklm,
                &format!("{root}\\DefaultIcon"),
                "",
                &icon_val,
                false,
            )?;
            write_reg_key_value(
                use_hklm,
                &format!("{root}\\shell\\open\\command"),
                "",
                &open_cmd,
                false,
            )?;
            records.push(FileAssociation::WindowsRegistryKey {
                hive: hive.clone(),
                subkey: root,
            });
        }

        // Own extensions `.nyx` and `.nys` cleanly set (Default) + OpenWithProgids.
        for (ext, prog_id) in [(".nyx", "nyx.Archive.1"), (".nys", "nyx.Script.1")] {
            let ext_key = format!("{classes_base}\\{ext}");
            let (created_key, prev) = write_reg_key_value(use_hklm, &ext_key, "", prog_id, true)?;
            records.push(FileAssociation::WindowsRegistryValue {
                hive: hive.clone(),
                subkey: ext_key.clone(),
                value_name: String::new(),
                set_value: prog_id.to_string(),
                previous_value: prev,
                created_key,
            });
            let ow_key = format!("{ext_key}\\OpenWithProgids");
            let (ow_created, _) = write_reg_key_value(use_hklm, &ow_key, prog_id, "", false)?;
            records.push(FileAssociation::WindowsRegistryValue {
                hive: hive.clone(),
                subkey: ow_key,
                value_name: prog_id.to_string(),
                set_value: String::new(),
                previous_value: None,
                created_key: ow_created,
            });
        }

        // Shared extensions `.json` and `.h5`: register ONLY OpenWithProgids (never touch `(Default)` or `UserChoice`).
        for (ext, prog_id) in [(".json", "nyx.Project.1"), (".h5", "nyx.Results.1")] {
            let ow_key = format!("{classes_base}\\{ext}\\OpenWithProgids");
            let (ow_created, _) = write_reg_key_value(use_hklm, &ow_key, prog_id, "", false)?;
            records.push(FileAssociation::WindowsRegistryValue {
                hive: hive.clone(),
                subkey: ow_key,
                value_name: prog_id.to_string(),
                set_value: String::new(),
                previous_value: None,
                created_key: ow_created,
            });
        }

        let app_key = format!("{classes_base}\\Applications\\nyx.exe");
        write_reg_key_value(
            use_hklm,
            &format!("{app_key}\\shell\\open\\command"),
            "",
            &open_cmd,
            false,
        )?;
        for ext in [".nyx", ".nys", ".json", ".h5"] {
            write_reg_key_value(
                use_hklm,
                &format!("{app_key}\\SupportedTypes"),
                ext,
                "",
                false,
            )?;
        }
        records.push(FileAssociation::WindowsRegistryKey {
            hive,
            subkey: app_key,
        });

        if opts.testing_reg_root.is_none() {
            // SAFETY: `SHCNE_ASSOCCHANGED | SHCNF_IDLIST` with null item pointers is the
            // documented Win32 shell notification indicating file-association updates.
            unsafe {
                SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, null(), null());
            }
        }
        Ok((records, vec![]))
    }

    pub fn remove_file_associations(assocs: &[FileAssociation]) -> Result<bool> {
        let mut changed = false;
        // First pass: restore / remove values carefully, leaving shared parent keys intact unless empty.
        for a in assocs {
            if let FileAssociation::WindowsRegistryValue {
                hive,
                subkey,
                value_name,
                set_value,
                previous_value,
                created_key,
            } = a
            {
                if subkey.is_empty() || !subkey.contains('\\') {
                    continue;
                }
                let use_hklm = hive.eq_ignore_ascii_case("HKLM");
                let w_sub = wide(subkey);
                let mut hkey: isize = 0;
                // SAFETY: `root_hkey(use_hklm)` is a valid predefined HKEY, `w_sub` is a non-empty
                // NUL-terminated UTF-16 subkey path, and `&mut hkey` is a valid out-pointer.
                if unsafe {
                    RegOpenKeyExW(
                        root_hkey(use_hklm),
                        w_sub.as_ptr(),
                        0,
                        KEY_READ | KEY_WRITE | KEY_WOW64_64KEY,
                        &mut hkey,
                    )
                } == ERROR_SUCCESS
                {
                    let cur = reg_query_string(hkey, value_name)
                        .ok()
                        .flatten()
                        .map(|(_, s)| s);
                    if cur.as_deref() == Some(set_value.as_str()) || set_value.is_empty() {
                        if let Some(prev) = previous_value {
                            let _ = reg_set_string(hkey, value_name, REG_SZ, prev);
                        } else {
                            let w_val = wide(value_name);
                            // SAFETY: `hkey` is open with `KEY_WRITE`; value name pointer is null
                            // for `(Default)` or NUL-terminated UTF-16 `w_val`.
                            unsafe {
                                RegDeleteValueW(
                                    hkey,
                                    if value_name.is_empty() {
                                        null()
                                    } else {
                                        w_val.as_ptr()
                                    },
                                );
                            }
                        }
                        changed = true;
                    }
                    let mut n_sub: u32 = 0;
                    let mut n_vals: u32 = 0;
                    // SAFETY: `hkey` is open with `KEY_READ`; only `&mut n_sub` and `&mut n_vals`
                    // are non-null, matching Win32 optional out-parameters for `RegQueryInfoKeyW`.
                    let info_ok = unsafe {
                        RegQueryInfoKeyW(
                            hkey,
                            null_mut(),
                            null_mut(),
                            null_mut(),
                            &mut n_sub,
                            null_mut(),
                            null_mut(),
                            &mut n_vals,
                            null_mut(),
                            null_mut(),
                            null_mut(),
                            null_mut(),
                        )
                    } == ERROR_SUCCESS;
                    // SAFETY: `hkey` was opened by `RegOpenKeyExW` above and is closed once.
                    unsafe {
                        RegCloseKey(hkey);
                    }
                    if *created_key && info_ok && n_sub == 0 && n_vals == 0 {
                        // SAFETY: `RegDeleteKeyW` only deletes an empty leaf key (`n_sub == 0 &&
                        // n_vals == 0`) that `install_file_associations` created (`subkey` verified
                        // non-empty and nested); `w_sub` is NUL-terminated UTF-16.
                        unsafe {
                            let _ = RegDeleteKeyW(root_hkey(use_hklm), w_sub.as_ptr());
                        }
                    }
                }
            }
        }
        // Second pass: remove exclusively-owned ProgID / Applications subtrees.
        for a in assocs {
            if let FileAssociation::WindowsRegistryKey { hive, subkey } = a {
                let use_hklm = hive.eq_ignore_ascii_case("HKLM");
                delete_registry_tree(use_hklm, subkey);
                changed = true;
            }
        }
        // SAFETY: Documented `SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, NULL, NULL)` call.
        unsafe {
            SHChangeNotify(SHCNE_ASSOCCHANGED, SHCNF_IDLIST, null(), null());
        }
        Ok(changed)
    }

    fn quote_win32_arg(word: &std::ffi::OsStr, out: &mut Vec<u16>) {
        const BACKSLASH: u16 = b'\\' as u16;
        const QUOTE: u16 = b'"' as u16;
        out.push(QUOTE);
        let mut backslashes = 0usize;
        for c in word.encode_wide() {
            if c == BACKSLASH {
                backslashes += 1;
                continue;
            }
            if c == QUOTE {
                out.extend(std::iter::repeat_n(BACKSLASH, backslashes * 2 + 1));
            } else {
                out.extend(std::iter::repeat_n(BACKSLASH, backslashes));
            }
            backslashes = 0;
            out.push(c);
        }
        out.extend(std::iter::repeat_n(BACKSLASH, backslashes * 2));
        out.push(QUOTE);
    }

    fn spawn_win32_detached(program: &Path, args: &[std::ffi::OsString]) -> Result<()> {
        let mut line: Vec<u16> = vec![];
        quote_win32_arg(program.as_os_str(), &mut line);
        for a in args {
            line.push(' ' as u16);
            quote_win32_arg(a, &mut line);
        }
        line.push(0);
        let app: Vec<u16> = program.as_os_str().encode_wide().chain(Some(0)).collect();
        let temp_dir = std::env::temp_dir();
        let w_cwd = wide_path(&temp_dir);
        // SAFETY: `StartupInfoW` is a `#[repr(C)]` POD struct of integers and pointers where zero
        // is valid; `cb` is set to the exact struct byte size immediately below.
        let mut si: StartupInfoW = unsafe { std::mem::zeroed() };
        si.cb = std::mem::size_of::<StartupInfoW>() as u32;
        // SAFETY: `ProcessInformation` is a `#[repr(C)]` POD struct of handles and IDs where zero
        // is valid before `CreateProcessW` initializes it.
        let mut pi: ProcessInformation = unsafe { std::mem::zeroed() };
        // DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW
        let flags = 0x0000_0008 | 0x0000_0200 | 0x0800_0000;
        // SAFETY: `app`, `line` (mutable buffer required by `CreateProcessW`), and `w_cwd` are
        // NUL-terminated UTF-16 buffers alive for the duration of the call; `&si` and `&mut pi`
        // are aligned valid pointers matching the Windows `STARTUPINFOW` / `PROCESS_INFORMATION` ABI.
        let ok = unsafe {
            CreateProcessW(
                app.as_ptr(),
                line.as_mut_ptr(),
                null(),
                null(),
                0,
                flags,
                null(),
                w_cwd.as_ptr(),
                &si,
                &mut pi,
            )
        };
        if ok == 0 {
            return Err(std::io::Error::last_os_error().into());
        }
        // SAFETY: `CreateProcessW` succeeded (`ok != 0`), so any non-null `pi.process` and
        // `pi.thread` handles are valid OS handles owned by us and closed exactly once.
        unsafe {
            if !pi.process.is_null() {
                CloseHandle(pi.process);
            }
            if !pi.thread.is_null() {
                CloseHandle(pi.thread);
            }
        }
        Ok(())
    }

    fn temp_reaper_copy(src_exe: &Path) -> Result<PathBuf> {
        let temp_dir = std::env::temp_dir();
        let copy_path = temp_dir.join(format!("nyx-uninstall-reaper-{}.exe", std::process::id()));
        std::fs::copy(src_exe, &copy_path)?;
        Ok(copy_path)
    }

    pub fn spawn_reap_final(
        uninstaller_exe: &Path,
        install_root: &Path,
        wait_pid: u32,
    ) -> Result<()> {
        let _ = std::env::set_current_dir(std::env::temp_dir());
        let src = if uninstaller_exe.is_file() {
            uninstaller_exe.to_path_buf()
        } else {
            std::env::current_exe()?
        };
        let copy_path = temp_reaper_copy(&src)?;
        let args = vec![
            std::ffi::OsString::from("--_reap-final"),
            std::ffi::OsString::from("--wait-pid"),
            std::ffi::OsString::from(wait_pid.to_string()),
            std::ffi::OsString::from("--root"),
            install_root.as_os_str().to_os_string(),
        ];
        spawn_win32_detached(&copy_path, &args)
    }

    pub fn execute_reap_final(wait_pid: u32, install_root: &Path) -> Result<()> {
        if !is_safe_reap_root(install_root) {
            bail!(
                "refusing reap_final on unsafe install root '{}'",
                install_root.display()
            );
        }
        let _ = std::env::set_current_dir(std::env::temp_dir());
        // SAFETY: `OpenProcess(SYNCHRONIZE, FALSE, wait_pid)` takes scalar arguments; `h` is
        // checked non-null (and != `INVALID_HANDLE_VALUE`) before waiting and closing via `CloseHandle`.
        unsafe {
            let h = OpenProcess(0x0010_0000 /* SYNCHRONIZE */, 0, wait_pid);
            if !h.is_null() && h != (-1isize as *mut c_void) {
                WaitForSingleObject(h, 10_000);
                CloseHandle(h);
            }
        }
        let bin_uninst = install_root.join("bin").join("nyx-uninstall.exe");
        for attempt in 0..40 {
            if !bin_uninst.exists() {
                break;
            }
            if std::fs::remove_file(&bin_uninst).is_ok() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(25 + attempt * 5));
        }
        let _ = std::fs::remove_dir(install_root.join("bin"));
        let _ = std::fs::remove_dir(install_root);

        if bin_uninst.is_file() && is_process_elevated() {
            let w_uninst = wide_path(&bin_uninst);
            let w_bin = wide_path(&install_root.join("bin"));
            let w_root = wide_path(install_root);
            // SAFETY: `install_root` passed `is_safe_reap_root`; `w_uninst`, `w_bin`, and `w_root`
            // are NUL-terminated UTF-16 buffers alive for the synchronous `MoveFileExW` calls, with
            // `lpNewFileName = NULL` and `MOVEFILE_DELAY_UNTIL_REBOOT`.
            unsafe {
                let _ = MoveFileExW(w_uninst.as_ptr(), null(), MOVEFILE_DELAY_UNTIL_REBOOT);
                let _ = MoveFileExW(w_bin.as_ptr(), null(), MOVEFILE_DELAY_UNTIL_REBOOT);
                let _ = MoveFileExW(w_root.as_ptr(), null(), MOVEFILE_DELAY_UNTIL_REBOOT);
            }
        }

        // Sweep older temporary reapers left in %TEMP%.
        let me = std::env::current_exe().ok();
        if let Ok(rd) = std::fs::read_dir(std::env::temp_dir()) {
            for e in rd.filter_map(Result::ok) {
                let p = e.path();
                let name = e.file_name().to_string_lossy().into_owned();
                if name.starts_with("nyx-uninstall-reaper-")
                    && name.ends_with(".exe")
                    && me.as_ref() != Some(&p)
                {
                    let _ = std::fs::remove_file(p);
                }
            }
        }
        Ok(())
    }

    pub fn evacuate_to_temp_and_relaunch(current_exe: &Path, args: &[String]) -> Result<bool> {
        let _ = std::env::set_current_dir(std::env::temp_dir());
        let copy_path = temp_reaper_copy(current_exe)?;
        let mut child_args: Vec<std::ffi::OsString> = Vec::new();
        for a in args {
            child_args.push(std::ffi::OsString::from(a));
        }
        if !args.iter().any(|a| a == "--prefix") {
            if let Some(root) = current_exe.parent().and_then(Path::parent) {
                child_args.push(std::ffi::OsString::from("uninstall"));
                child_args.push(std::ffi::OsString::from("--prefix"));
                child_args.push(root.as_os_str().to_os_string());
            }
        }
        spawn_win32_detached(&copy_path, &child_args)?;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::{is_safe_nyx_registry_subtree, is_safe_reap_root, ARP_SUBKEY};
    use std::path::Path;

    #[test]
    fn registry_tree_scope_guard_allows_only_nyx_subtrees() {
        assert!(is_safe_nyx_registry_subtree(ARP_SUBKEY));
        assert!(is_safe_nyx_registry_subtree(
            r"Software\nyx-install-test-reg-1234"
        ));
        assert!(is_safe_nyx_registry_subtree(
            r"Software\nyx-install-test-reg-1234\Uninstall\nyx"
        ));
        assert!(is_safe_nyx_registry_subtree(
            r"Software\Classes\nyx.Archive.1"
        ));
        assert!(is_safe_nyx_registry_subtree(
            r"Software\Classes\Applications\nyx.exe"
        ));

        for unsafe_key in [
            "",
            "Software",
            r"Software\Microsoft",
            r"Software\Microsoft\Windows\CurrentVersion\Uninstall",
            r"Software\Classes",
            r"Software\Classes\.json",
            r"Software\Classes\Applications",
            "Environment",
            r"Software\nyx\..\Microsoft",
            "Software\0bad",
        ] {
            assert!(
                !is_safe_nyx_registry_subtree(unsafe_key),
                "guard must reject {unsafe_key:?}"
            );
        }
    }

    #[test]
    fn reap_root_scope_guard_rejects_drive_and_system_roots() {
        assert!(is_safe_reap_root(Path::new(
            r"C:\Users\Alice\AppData\Local\Programs\nyx"
        )));
        assert!(is_safe_reap_root(Path::new(r"C:\Program Files\nyx")));

        for bad in [
            "",
            ".",
            r"C:\",
            "C:",
            r"C:\Windows",
            r"C:\Program Files",
            r"C:\Program Files (x86)",
            r"C:\Users",
        ] {
            assert!(
                !is_safe_reap_root(Path::new(bad)),
                "guard must reject reap root {bad:?}"
            );
        }
    }
}
