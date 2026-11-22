// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tier-1 hermetic multi-OS fault-injection, safety, trailer, `nyx-connect`, and rollback tests.

use anyhow::Result;
use nyx_connect::{connect, status as connect_status, Status};
use nyx_install::bundle::{
    locate_embedded_trailer, write_nyxpack, write_self_extracting, IntegrityVerifier, PackSpec,
    PayloadSource, PayloadVerifier, StagedItem, CODEC_ZLIB_RS, TRAILER_LEN, TRAILER_MAGIC,
};
use nyx_install::connect_hook::ConnectPolicy;
use nyx_install::manifest::{
    resolve_manifest_path, FileBase, FileRole, InstallManifest, InstallScope, TargetOs,
};
use nyx_install::platform::linux::MARK_PNG;
use nyx_install::platform::macos::build_icns;
use nyx_install::platform::windows::{
    display_icon_value, edit_path_string, quote_command_exe, remove_from_path_string,
};
use nyx_install::platform::{upsert_marked_block, FakePlatformBackend, HostPaths};
use nyx_install::process::{ProcessInspector, QueueStatus, RunningProc};
use nyx_install::safety::{file_blake3_hex, validate_install_root, validate_rel_path};
use nyx_install::tx::{run_doctor, run_install_or_upgrade, run_repair, run_uninstall, TxOptions};
use std::fs;
use std::io::Cursor;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

struct TempSandbox {
    root: PathBuf,
}

impl TempSandbox {
    fn new(name: &str) -> Self {
        static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let id = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "nyx-inst-hermetic-{name}-{}-{id}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(root.join("home")).unwrap();
        fs::create_dir_all(root.join("tmp")).unwrap();
        Self { root }
    }
}

impl Drop for TempSandbox {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

#[derive(Default)]
struct FakeProcessInspector {
    procs: Mutex<Vec<RunningProc>>,
    terminated: Mutex<Vec<u32>>,
    queue: Mutex<QueueStatus>,
    paused_calls: Mutex<usize>,
    resumed_calls: Mutex<usize>,
}

impl ProcessInspector for FakeProcessInspector {
    fn running_in_root(&self, root: &Path) -> Vec<RunningProc> {
        self.procs
            .lock()
            .unwrap()
            .iter()
            .filter(|p| p.exe_path.starts_with(root))
            .cloned()
            .collect()
    }

    fn terminate_pid(&self, pid: u32) -> Result<()> {
        self.terminated.lock().unwrap().push(pid);
        self.procs.lock().unwrap().retain(|p| p.pid != pid);
        Ok(())
    }

    fn shared_queue_status(&self, _workspace: &Path) -> QueueStatus {
        self.queue.lock().unwrap().clone()
    }

    fn pause_shared_queue(&self, _workspace: &Path) -> Result<bool> {
        *self.paused_calls.lock().unwrap() += 1;
        Ok(true)
    }

    fn resume_shared_queue(&self, _workspace: &Path) -> Result<()> {
        *self.resumed_calls.lock().unwrap() += 1;
        Ok(())
    }
}

fn build_test_nyxpack(
    dir: &Path,
    os: TargetOs,
    version: &str,
    files: &[(&str, FileBase, FileRole, &[u8])],
) -> PathBuf {
    let stage = dir.join(format!("src_{version}"));
    fs::create_dir_all(&stage).unwrap();
    let mut items = Vec::new();
    for (idx, (rel, base, role, content)) in files.iter().enumerate() {
        let f = stage.join(format!("item_{idx}.bin"));
        fs::write(&f, content).unwrap();
        items.push(StagedItem {
            source_path: f,
            rel_path: (*rel).to_string(),
            base: *base,
            role: *role,
            mode: if *base == FileBase::Bin { 0o755 } else { 0o644 },
        });
    }
    let spec = PackSpec {
        nyx_version: version.to_string(),
        target_triple: "test-target".into(),
        os,
        items,
    };
    let pack = dir.join(format!("bundle_{version}.nyxpack"));
    write_nyxpack(&spec, &pack).unwrap();
    pack
}

fn base_tx_opts(prefix: PathBuf) -> TxOptions {
    let mut o = TxOptions::default_user();
    o.prefix = Some(prefix);
    o.scope = InstallScope::User;
    o.run_mcp_smoke = false;
    o.connect_policy = ConnectPolicy::None;
    o
}

#[test]
fn hermetic_install_upgrade_repair_uninstall_across_all_three_oses() {
    for os in [TargetOs::Windows, TargetOs::Linux, TargetOs::Macos] {
        let sb = TempSandbox::new(os.as_str());
        let paths = HostPaths::sandboxed(&sb.root, os);
        let docs_nyx = paths.documents_dir.as_ref().unwrap();
        fs::create_dir_all(docs_nyx).unwrap();
        fs::write(docs_nyx.join("workflow.txt"), b"user project data").unwrap();

        let prefix = if os == TargetOs::Macos {
            sb.root.join("Apps/Nyx.app")
        } else {
            sb.root.join("install_root")
        };

        let bin_name = if os == TargetOs::Windows {
            "nyx.exe"
        } else {
            "nyx"
        };
        let mcp_name = if os == TargetOs::Windows {
            "nyx-mcp.exe"
        } else {
            "nyx-mcp"
        };
        let uninst_name = if os == TargetOs::Windows {
            "nyx-uninstall.exe"
        } else {
            "nyx-uninstall"
        };

        let p1 = build_test_nyxpack(
            &sb.root,
            os,
            "0.1.0",
            &[
                (bin_name, FileBase::Bin, FileRole::Binary, b"nyx-v1"),
                (mcp_name, FileBase::Bin, FileRole::Binary, b"mcp-v1"),
                (
                    uninst_name,
                    FileBase::Bin,
                    FileRole::Uninstaller,
                    b"uninst-v1",
                ),
                (
                    "materials/si.toml",
                    FileBase::Resource,
                    FileRole::Asset,
                    b"name = \"silicon\"\n",
                ),
                (
                    "obsolete.txt",
                    FileBase::Resource,
                    FileRole::Asset,
                    b"old asset",
                ),
            ],
        );
        let backend = FakePlatformBackend::new(os);
        let inspector = FakeProcessInspector::default();
        let opts = base_tx_opts(prefix.clone());

        let r1 = run_install_or_upgrade(
            "install",
            &PayloadSource::NyxPack { path: p1 },
            &opts,
            &backend,
            &inspector,
            &paths,
            &mut |_| {},
        )
        .unwrap();
        assert_eq!(r1.action, "install");
        assert!(r1.ok);

        let manifest_path =
            resolve_manifest_path(os, InstallScope::Custom, &prefix, paths.home.as_deref());
        assert!(
            manifest_path.is_file(),
            "manifest missing at {}",
            manifest_path.display()
        );
        if os == TargetOs::Macos {
            // Manifest MUST reside outside Nyx.app!
            assert!(
                !manifest_path.starts_with(&prefix),
                "macOS manifest must sit outside Nyx.app"
            );
            assert!(prefix.join("Contents/Info.plist").is_file());
            assert!(prefix.join("Contents/Resources/Nyx.icns").is_file());
        }

        // Upgrade to 0.2.0: removes obsolete.txt and updates nyx binary.
        let p2 = build_test_nyxpack(
            &sb.root,
            os,
            "0.2.0",
            &[
                (bin_name, FileBase::Bin, FileRole::Binary, b"nyx-v2"),
                (mcp_name, FileBase::Bin, FileRole::Binary, b"mcp-v2"),
                (
                    uninst_name,
                    FileBase::Bin,
                    FileRole::Uninstaller,
                    b"uninst-v2",
                ),
                (
                    "materials/si.toml",
                    FileBase::Resource,
                    FileRole::Asset,
                    b"name = \"silicon\"\n",
                ),
            ],
        );
        let r2 = run_install_or_upgrade(
            "install",
            &PayloadSource::NyxPack { path: p2.clone() },
            &opts,
            &backend,
            &inspector,
            &paths,
            &mut |_| {},
        )
        .unwrap();
        assert_eq!(r2.action, "upgrade");
        assert_eq!(r2.previous_version.as_deref(), Some("0.1.0"));
        assert!(
            !prefix.join("obsolete.txt").exists(),
            "obsolete file should be removed on upgrade"
        );

        // Corrupt a binary and verify `doctor` detects it, then `repair` restores it.
        let m: InstallManifest =
            serde_json::from_str(&fs::read_to_string(&manifest_path).unwrap()).unwrap();
        let nyx_abs = PathBuf::from(&m.bin_dir).join(bin_name);
        fs::write(&nyx_abs, b"CORRUPTED").unwrap();
        let doc_bad = run_doctor(&opts, &backend, &paths).unwrap();
        assert!(!doc_bad.ok);
        let _rep = run_repair(
            &PayloadSource::NyxPack { path: p2 },
            &opts,
            &backend,
            &inspector,
            &paths,
            &mut |_| {},
        )
        .unwrap();
        assert_eq!(fs::read(&nyx_abs).unwrap(), b"nyx-v2");
        assert!(run_doctor(&opts, &backend, &paths).unwrap().ok);

        // Create temporary dir + solve_history and uninstall with --purge-cache.
        fs::create_dir_all(paths.temp_dir.join("nyx-scratch")).unwrap();
        fs::write(docs_nyx.join("solve_history.jsonl"), b"{}\n").unwrap();
        let mut un_opts = opts.clone();
        un_opts.purge_cache = true;
        let run_u = run_uninstall(&un_opts, &backend, &inspector, &paths, &mut |_| {}).unwrap();
        assert!(run_u.ok);
        assert!(
            !prefix.exists(),
            "{}: install_root {} still exists after uninstall",
            os.as_str(),
            prefix.display()
        );
        assert!(!paths.temp_dir.join("nyx-scratch").exists());
        assert!(!docs_nyx.join("solve_history.jsonl").exists());
        assert_eq!(
            fs::read(docs_nyx.join("workflow.txt")).unwrap(),
            b"user project data"
        );
    }
}

#[test]
fn rollback_journal_restores_previous_state_on_mid_phase_failure() {
    for fail_phase in ["path", "desktop", "file_assoc"] {
        let sb = TempSandbox::new(fail_phase);
        let paths = HostPaths::sandboxed(&sb.root, TargetOs::Linux);
        let prefix = sb.root.join("inst");
        let p1 = build_test_nyxpack(
            &sb.root,
            TargetOs::Linux,
            "0.1.0",
            &[("nyx", FileBase::Bin, FileRole::Binary, b"GOOD_V1")],
        );
        let p2 = build_test_nyxpack(
            &sb.root,
            TargetOs::Linux,
            "0.2.0",
            &[("nyx", FileBase::Bin, FileRole::Binary, b"NEW_V2")],
        );

        let backend = FakePlatformBackend::new(TargetOs::Linux);
        let inspector = FakeProcessInspector::default();
        let opts = base_tx_opts(prefix.clone());
        run_install_or_upgrade(
            "install",
            &PayloadSource::NyxPack { path: p1 },
            &opts,
            &backend,
            &inspector,
            &paths,
            &mut |_| {},
        )
        .unwrap();

        // Inject failure on upgrade in Phase 5.
        backend.state.lock().unwrap().inject_fail_phase = Some(fail_phase);
        let err = run_install_or_upgrade(
            "upgrade",
            &PayloadSource::NyxPack { path: p2 },
            &opts,
            &backend,
            &inspector,
            &paths,
            &mut |_| {},
        );
        assert!(err.is_err(), "expected failure at {fail_phase}");
        assert_eq!(
            fs::read(prefix.join("bin/nyx")).unwrap(),
            b"GOOD_V1",
            "rollback must restore v1 binary"
        );
        assert!(
            fs::read_dir(&sb.root)
                .unwrap()
                .filter_map(Result::ok)
                .all(|e| !e.file_name().to_string_lossy().contains(".journal.jsonl")),
            "journal must be cleaned up after rollback"
        );
    }
}

#[test]
fn uninstall_preserves_user_modified_files_unless_forced() {
    let sb = TempSandbox::new("user_mod");
    let paths = HostPaths::sandboxed(&sb.root, TargetOs::Linux);
    let prefix = sb.root.join("inst");
    let p1 = build_test_nyxpack(
        &sb.root,
        TargetOs::Linux,
        "0.1.0",
        &[
            ("nyx", FileBase::Bin, FileRole::Binary, b"BIN"),
            (
                "materials/custom.toml",
                FileBase::Resource,
                FileRole::Asset,
                b"orig = 1\n",
            ),
        ],
    );
    let backend = FakePlatformBackend::new(TargetOs::Linux);
    let inspector = FakeProcessInspector::default();
    let opts = base_tx_opts(prefix.clone());
    run_install_or_upgrade(
        "install",
        &PayloadSource::NyxPack { path: p1 },
        &opts,
        &backend,
        &inspector,
        &paths,
        &mut |_| {},
    )
    .unwrap();

    // User modifies materials/custom.toml.
    let mod_file = prefix.join("materials/custom.toml");
    fs::write(&mod_file, b"user_edited = true\n").unwrap();

    let rep = run_uninstall(&opts, &backend, &inspector, &paths, &mut |_| {}).unwrap();
    assert!(
        mod_file.is_file(),
        "modified file must be preserved when --force is false"
    );
    assert!(
        !prefix.join("bin/nyx").exists(),
        "unmodified binary must still be removed"
    );
    assert!(rep.warnings.iter().any(|w| w.contains("custom.toml")));
}

#[test]
fn poison_prefixes_and_unsafe_archive_entries_are_rejected() {
    let sb = TempSandbox::new("poison");
    let home = sb.root.join("home");
    let docs_nyx = home.join("Documents/nyx");
    fs::create_dir_all(&docs_nyx).unwrap();
    fs::write(docs_nyx.join("workspace.json"), b"{}").unwrap();

    for bad in [
        Path::new("/"),
        &home,
        &docs_nyx,
        Path::new("/usr/local"),
        Path::new("/Applications"),
    ] {
        assert!(
            validate_install_root(
                bad,
                TargetOs::Linux,
                InstallScope::User,
                Some(&home),
                None,
                None,
                false
            )
            .is_err(),
            "should reject {}",
            bad.display()
        );
    }
    #[cfg(windows)]
    {
        assert!(validate_install_root(
            Path::new(r"C:\"),
            TargetOs::Windows,
            InstallScope::System,
            Some(&home),
            None,
            Some(Path::new(r"C:\Program Files")),
            false
        )
        .is_err());
    }

    // Non-empty unrelated directory without --force refused.
    let occupied = sb.root.join("occupied");
    fs::create_dir_all(&occupied).unwrap();
    fs::write(occupied.join("precious.dat"), b"123").unwrap();
    assert!(validate_install_root(
        &occupied,
        TargetOs::Linux,
        InstallScope::Custom,
        Some(&home),
        None,
        None,
        false
    )
    .is_err());

    for bad_rel in [
        "",
        "../escape",
        "a/../b",
        "/abs/file",
        r"\win\abs",
        "foo:ads",
        "CON.txt",
        "nul",
        "bad\0name",
    ] {
        assert!(
            validate_rel_path(bad_rel).is_err(),
            "must reject unsafe rel path '{bad_rel}'"
        );
    }
    assert!(validate_rel_path("materials/si.toml").is_ok());
}

#[test]
fn pe_authenticode_aware_trailer_locates_payload_even_with_cert_table_appended() {
    let sb = TempSandbox::new("pe_trailer");
    let p1 = build_test_nyxpack(
        &sb.root,
        TargetOs::Windows,
        "0.1.0",
        &[(
            "nyx.exe",
            FileBase::Bin,
            FileRole::Binary,
            b"MZ-PAYLOAD-BYTES",
        )],
    );
    let payload_bytes = fs::read(&p1).unwrap();
    let payload_b3 = *blake3::hash(&payload_bytes).as_bytes();
    let toc_len = u32::from_le_bytes(payload_bytes[8..12].try_into().unwrap()) as usize;
    let toc_b3 = *blake3::hash(&payload_bytes[12..12 + toc_len]).as_bytes();

    // Build a minimal valid PE32+ stub with DataDirectory[SECURITY] slot at offset e_lfanew + 24 + 112 + 32.
    let stub = sb.root.join("stub.exe");
    let mut pe = vec![0u8; 512];
    pe[0..2].copy_from_slice(b"MZ");
    pe[60..64].copy_from_slice(&0x80u32.to_le_bytes()); // e_lfanew = 0x80
    pe[0x80..0x84].copy_from_slice(b"PE\0\0");
    let opt = 0x80 + 24;
    pe[opt..opt + 2].copy_from_slice(&0x20Bu16.to_le_bytes()); // PE32+
    pe[opt + 108..opt + 112].copy_from_slice(&16u32.to_le_bytes()); // NumberOfRvaAndSizes = 16
    fs::write(&stub, &pe).unwrap();

    let sfx = sb.root.join("signed_setup.exe");
    let trailer = write_self_extracting(
        &stub,
        &p1,
        &sfx,
        &payload_bytes[12..12 + toc_len],
        toc_b3,
        payload_b3,
    )
    .unwrap();
    assert_eq!(trailer.codec, CODEC_ZLIB_RS);

    // Unsigned EOF - 96 works and is 8-byte aligned.
    let mut signed_bytes = fs::read(&sfx).unwrap();
    assert_eq!(
        signed_bytes.len() % 8,
        0,
        "trailer end MUST be 8-byte aligned for Authenticode"
    );
    let cert_offset = signed_bytes.len() as u32;
    let fake_cert = vec![0xAAu8; 256];
    let dd4 = opt + 112 + 4 * 8;
    signed_bytes[dd4..dd4 + 4].copy_from_slice(&cert_offset.to_le_bytes());
    signed_bytes[dd4 + 4..dd4 + 8].copy_from_slice(&(fake_cert.len() as u32).to_le_bytes());
    signed_bytes.extend_from_slice(&fake_cert);

    let (found, trailer_start) = locate_embedded_trailer(&mut Cursor::new(&signed_bytes))
        .unwrap()
        .expect("should find NYXINST1 before Authenticode cert table");
    assert_eq!(trailer_start as usize, cert_offset as usize - TRAILER_LEN);
    assert_eq!(found, trailer);
    assert_eq!(&TRAILER_MAGIC[..], b"NYXINST1");

    // Tamper 1 payload byte -> IntegrityVerifier must catch it!
    signed_bytes[512 + 12 + toc_len] ^= 0xFF;
    let actual_b3 = *blake3::hash(&signed_bytes[512..512 + payload_bytes.len()]).as_bytes();
    assert!(IntegrityVerifier
        .verify_payload_stream(&found, &actual_b3)
        .is_err());
}

#[test]
fn windows_path_quoting_and_linux_symlinks_and_macos_icns_invariants() {
    // 1. Windows Path edits preserve %SystemRoot% verbatim and handle case/quotes/slashes.
    let (added, changed, has_pct) = edit_path_string(
        r"%SystemRoot%\system32;C:\Tools\",
        r"C:\Users\Alice Smith\AppData\Local\Programs\nyx\bin",
        None,
    );
    assert!(changed && has_pct);
    assert!(added.starts_with(r"%SystemRoot%\system32;"));
    let (_, changed2, _) = edit_path_string(
        &added,
        r"c:/users/alice smith/appdata/local/programs/nyx/bin\",
        None,
    );
    assert!(
        !changed2,
        "case/slash equivalent PATH entry must be idempotent"
    );
    let (removed, rem_changed) = remove_from_path_string(
        &added,
        r"C:\Users\Alice Smith\AppData\Local\Programs\nyx\bin",
    );
    assert!(rem_changed);
    assert_eq!(removed, r"%SystemRoot%\system32;C:\Tools\");

    // DisplayIcon is unquoted; UninstallString is strictly quoted.
    let nyx_p = Path::new(r"C:\Program Files\nyx\bin\nyx.exe");
    assert_eq!(
        display_icon_value(nyx_p),
        r"C:\Program Files\nyx\bin\nyx.exe,0"
    );
    assert_eq!(
        quote_command_exe(nyx_p, "--yes"),
        r#""C:\Program Files\nyx\bin\nyx.exe" --yes"#
    );

    // 2. Linux profile block preserves symlink target without overwriting symlink node.
    let sb = TempSandbox::new("dotfile");
    let real_rc = sb.root.join("dotfiles/zshrc");
    fs::create_dir_all(real_rc.parent().unwrap()).unwrap();
    fs::write(&real_rc, "# existing zsh config\n").unwrap();
    let sym_rc = sb.root.join("home/.zshrc");
    #[cfg(unix)]
    std::os::unix::fs::symlink(&real_rc, &sym_rc).unwrap();
    #[cfg(not(unix))]
    fs::write(&sym_rc, "# existing zsh config\n").unwrap();
    upsert_marked_block(
        &sym_rc,
        "nyx installer",
        "export PATH=\"$HOME/.local/bin:$PATH\"",
        false,
        false,
    )
    .unwrap();
    assert!(fs::read_to_string(&sym_rc)
        .unwrap()
        .contains("# >>> nyx installer >>>"));
    assert!(!sb.root.join("home/.bash_profile").exists());

    // 3. macOS pure-Rust .icns generator starts with `icns` + valid PNG chunks `ic07`..`ic10`.
    let icns = build_icns(MARK_PNG).unwrap();
    assert_eq!(&icns[0..4], b"icns");
    let total = u32::from_be_bytes(icns[4..8].try_into().unwrap()) as usize;
    assert_eq!(total, icns.len());
    assert!(icns.windows(4).any(|w| w == b"ic07") && icns.windows(4).any(|w| w == b"ic10"));
}

#[test]
fn connect_lifecycle_respects_user_choice_and_surgical_uninstall() {
    let sb = TempSandbox::new("connect_life");
    let os = TargetOs::Linux;
    let paths = HostPaths::sandboxed(&sb.root, os);
    let env = paths.connect_env(os).unwrap();
    fs::create_dir_all(sb.root.join("home/.cursor")).unwrap();
    fs::create_dir_all(sb.root.join("home/.claude")).unwrap();

    // Point claude-code at an external/foreign nyx-mcp installation.
    let other_mcp = sb.root.join("other").join("nyx-mcp");
    fs::create_dir_all(other_mcp.parent().unwrap()).unwrap();
    fs::write(&other_mcp, b"other").unwrap();
    connect(&env, "claude-code", &other_mcp).unwrap();

    // Install with ConnectPolicy::Explicit(["cursor", "claude-code"]) -> wait: test Detected vs surgical remove:
    // First connect cursor via Detected (cursor has .cursor dir, claude-code points to other_mcp).
    // Note `other_mcp` does NOT have deps/.fingerprint, so `Detected` treats non-cargo Stale unless we also connect explicitly,
    // wait: if cursor is connected to installed mcp and claude-code is re-pointed to `other_mcp` after install:
    let prefix = sb.root.join("inst");
    let pack = build_test_nyxpack(
        &sb.root,
        os,
        "0.1.0",
        &[
            ("nyx", FileBase::Bin, FileRole::Binary, b"NYX"),
            ("nyx-mcp", FileBase::Bin, FileRole::Binary, b"MCP"),
            ("nyx-uninstall", FileBase::Bin, FileRole::Uninstaller, b"UN"),
        ],
    );
    let backend = FakePlatformBackend::new(os);
    let inspector = FakeProcessInspector::default();
    let mut opts = base_tx_opts(prefix.clone());
    opts.connect_policy = ConnectPolicy::Explicit(vec!["cursor".into()]);
    let r = run_install_or_upgrade(
        "install",
        &PayloadSource::NyxPack { path: pack },
        &opts,
        &backend,
        &inspector,
        &paths,
        &mut |_| {},
    )
    .unwrap();
    assert!(r.ok);
    let inst_mcp = prefix.join("bin").join("nyx-mcp");
    let cursor_h = nyx_connect::hosts(&env)
        .into_iter()
        .find(|h| h.id == "cursor")
        .unwrap();
    let claude_h = nyx_connect::hosts(&env)
        .into_iter()
        .find(|h| h.id == "claude-code")
        .unwrap();
    assert_eq!(connect_status(&cursor_h, &inst_mcp), Status::Connected);
    assert_eq!(connect_status(&claude_h, &other_mcp), Status::Connected);

    // Uninstall removes cursor (pointing to `inst_mcp`) while leaving claude-code (`other_mcp`) untouched!
    let ur = run_uninstall(&opts, &backend, &inspector, &paths, &mut |_| {}).unwrap();
    assert!(ur.ok);
    assert_eq!(connect_status(&cursor_h, &inst_mcp), Status::NotConnected);
    assert_eq!(connect_status(&claude_h, &other_mcp), Status::Connected);
}

#[test]
fn busy_process_and_open_queue_require_force_before_mutating() {
    let sb = TempSandbox::new("busy_proc");
    let paths = HostPaths::sandboxed(&sb.root, TargetOs::Linux);
    let prefix = sb.root.join("inst");
    let pack = build_test_nyxpack(
        &sb.root,
        TargetOs::Linux,
        "0.1.0",
        &[("nyx", FileBase::Bin, FileRole::Binary, b"A")],
    );

    let backend = FakePlatformBackend::new(TargetOs::Linux);
    let inspector = FakeProcessInspector::default();
    inspector.procs.lock().unwrap().push(RunningProc {
        pid: 4242,
        exe_path: prefix.join("bin/nyx"),
        name: "nyx".into(),
    });
    let opts = base_tx_opts(prefix.clone());
    let err = run_install_or_upgrade(
        "install",
        &PayloadSource::NyxPack { path: pack.clone() },
        &opts,
        &backend,
        &inspector,
        &paths,
        &mut |_| {},
    )
    .unwrap_err();
    assert_eq!(nyx_install::tx::exit_code_of(&err), 4);
    assert!(!prefix.exists());

    // With force = true, inspector terminates PID and pauses/resumes queue.
    *inspector.queue.lock().unwrap() = QueueStatus {
        active_runner: true,
        open_items: 1,
        paused: false,
    };
    let mut forced = opts;
    forced.force = true;
    let ok = run_install_or_upgrade(
        "install",
        &PayloadSource::NyxPack { path: pack },
        &forced,
        &backend,
        &inspector,
        &paths,
        &mut |_| {},
    )
    .unwrap();
    assert!(ok.ok);
    assert_eq!(*inspector.terminated.lock().unwrap(), vec![4242]);
    assert_eq!(*inspector.paused_calls.lock().unwrap(), 1);
    assert_eq!(*inspector.resumed_calls.lock().unwrap(), 1);
    let _ = file_blake3_hex(&prefix.join("bin/nyx")).unwrap();
}
