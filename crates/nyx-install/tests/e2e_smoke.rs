// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Tier-2 real end-to-end installer test:
//! 1. Invokes `nyx-pack` over real release binaries into a self-extracting `nyx-setup` with 96B `NYXINST1` trailer.
//! 2. Runs `install --prefix <tmp>/inst --no-path --no-desktop --no-file-assoc --connect none --yes --json`.
//! 3. Executes installed `nyx --version`, `nyx --batch .../01_capacitor.nyx.json` from an empty temp CWD, `nyx run hello.nys`,
//!    and `nyx_connect::test_server_within` against installed `nyx-mcp`.
//! 4. Executes installed `nyx-uninstall uninstall --yes --json` and verifies full cleanup + `Documents/nyx` preservation.

use nyx_install::bundle::exe_filename;
use nyx_install::manifest::TargetOs;
use nyx_install::report::InstallReport;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{Duration, Instant};

fn exe_dir_from_pack() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-pack"))
        .parent()
        .expect("nyx-pack parent dir")
        .to_path_buf()
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap()
}

#[test]
fn pack_install_batch_script_mcp_and_uninstall_preserve_user_data() {
    let target_dir = exe_dir_from_pack();
    let os = TargetOs::current();
    for stem in [
        "nyx",
        "nyx-mcp",
        "nyx-worker",
        "nyx-script-worker",
        "nyx-connect",
        "nyx-setup",
    ] {
        let p = target_dir.join(exe_filename(stem, os));
        assert!(
            p.is_file(),
            "missing compiled binary {} (build workspace binaries first: cargo build --release -p nyx -p nyx-engine -p nyx-mcp -p nyx-script -p nyx-connect -p nyx-install)",
            p.display()
        );
    }

    let tmp_root = std::env::temp_dir().join(format!("nyx-e2e-smoke-{}", std::process::id()));
    let _ = fs::remove_dir_all(&tmp_root);
    let inst_root = if os == TargetOs::Macos {
        tmp_root.join("inst/Nyx.app")
    } else {
        tmp_root.join("inst")
    };
    let empty_cwd = tmp_root.join("empty_cwd");
    let fake_home = tmp_root.join("fake_home");
    let docs_nyx = fake_home.join("Documents/nyx");
    fs::create_dir_all(&empty_cwd).unwrap();
    fs::create_dir_all(&docs_nyx).unwrap();
    fs::write(docs_nyx.join("keep.txt"), b"must survive uninstall").unwrap();

    // 1. Pack self-extracting installer.
    let sfx_path = tmp_root.join(exe_filename("nyx-setup-e2e", os));
    let pack_out = Command::new(env!("CARGO_BIN_EXE_nyx-pack"))
        .args([
            "--target-dir",
            &target_dir.to_string_lossy(),
            "--repo-root",
            &repo_root().to_string_lossy(),
            "--stub",
            env!("CARGO_BIN_EXE_nyx-setup"),
            "--self-extract",
            &sfx_path.to_string_lossy(),
        ])
        .output()
        .expect("run nyx-pack");
    assert!(
        pack_out.status.success(),
        "nyx-pack failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&pack_out.stdout),
        String::from_utf8_lossy(&pack_out.stderr)
    );
    assert!(sfx_path.is_file());

    // 2. Run self-extracting installer into `<tmp>/inst`.
    let install_prefix_arg = if os == TargetOs::Macos {
        tmp_root.join("inst")
    } else {
        inst_root.clone()
    };
    let inst_out = Command::new(&sfx_path)
        .current_dir(&empty_cwd)
        .env("HOME", &fake_home)
        .env("USERPROFILE", &fake_home)
        .env("NYX_WORKSPACE", &docs_nyx)
        .args([
            "install",
            "--prefix",
            &install_prefix_arg.to_string_lossy(),
            "--no-path",
            "--no-desktop",
            "--no-file-assoc",
            "--connect",
            "none",
            "--yes",
            "--json",
        ])
        .output()
        .expect("run sfx install");
    assert!(
        inst_out.status.success(),
        "install failed (code {:?}):\nstdout: {}\nstderr: {}",
        inst_out.status.code(),
        String::from_utf8_lossy(&inst_out.stdout),
        String::from_utf8_lossy(&inst_out.stderr)
    );
    let rep: InstallReport = serde_json::from_slice(&inst_out.stdout).unwrap_or_else(|e| {
        panic!(
            "invalid install JSON ({e}): {}",
            String::from_utf8_lossy(&inst_out.stdout)
        )
    });
    assert!(rep.ok && rep.action == "install");

    let (bin_dir, res_dir) = if os == TargetOs::Macos {
        (
            inst_root.join("Contents/MacOS"),
            inst_root.join("Contents/Resources"),
        )
    } else {
        (inst_root.join("bin"), inst_root.clone())
    };
    let installed_nyx = bin_dir.join(exe_filename("nyx", os));
    let installed_mcp = bin_dir.join(exe_filename("nyx-mcp", os));
    let installed_uninst = bin_dir.join(exe_filename("nyx-uninstall", os));
    assert!(installed_nyx.is_file());
    assert!(res_dir.join("THIRD_PARTY_NOTICES.txt").is_file());

    // 3a. Headless `nyx --version`.
    let ver_out = Command::new(&installed_nyx)
        .current_dir(&empty_cwd)
        .arg("--version")
        .output()
        .unwrap();
    assert!(ver_out.status.success());
    assert!(String::from_utf8_lossy(&ver_out.stdout).starts_with("nyx "));

    // 3b. Headless `nyx --batch` on bundled example from empty temp CWD (no repo checkout!).
    let cap_example = res_dir.join("examples/01_capacitor.nyx.json");
    let batch_out = Command::new(&installed_nyx)
        .current_dir(&empty_cwd)
        .env("HOME", &fake_home)
        .env("USERPROFILE", &fake_home)
        .env("NYX_WORKSPACE", &docs_nyx)
        .args(["--batch", &cap_example.to_string_lossy()])
        .output()
        .unwrap();
    assert!(
        batch_out.status.success(),
        "installed `nyx --batch` failed from empty CWD:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&batch_out.stdout),
        String::from_utf8_lossy(&batch_out.stderr)
    );

    // 3c. Isolated `.nys` script execution (proves nyx-script-worker & set_self_as_worker).
    let nys_path = empty_cwd.join("smoke.nys");
    fs::write(
        &nys_path,
        "nyx.project_new(title=\"e2e_nys\", dim=2)\ninfo = nyx.project_info()\nprint(\"TITLE:\" + str(info[\"title\"]))\n",
    )
    .unwrap();
    let run_out = Command::new(&installed_nyx)
        .current_dir(&empty_cwd)
        .env("HOME", &fake_home)
        .env("USERPROFILE", &fake_home)
        .args([
            "run",
            &nys_path.to_string_lossy(),
            "--workspace",
            &empty_cwd.to_string_lossy(),
        ])
        .output()
        .unwrap();
    assert!(
        run_out.status.success(),
        "installed `nyx run smoke.nys` failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&run_out.stdout),
        String::from_utf8_lossy(&run_out.stderr)
    );
    assert!(String::from_utf8_lossy(&run_out.stdout).contains("TITLE:e2e_nys"));

    // 3d. MCP server handshake probe on installed `nyx-mcp`.
    let mcp_info = nyx_connect::test_server_within(&installed_mcp, Duration::from_secs(10))
        .expect("installed nyx-mcp handshake must succeed");
    assert!(mcp_info.tools.contains(&"script_run".to_string()));

    // 4. Uninstall using installed `nyx-uninstall` itself (exercises Windows Mode B --_reap-final!).
    let uninst_out = Command::new(&installed_uninst)
        .current_dir(&empty_cwd)
        .env("HOME", &fake_home)
        .env("USERPROFILE", &fake_home)
        .env("NYX_WORKSPACE", &docs_nyx)
        .args(["uninstall", "--yes", "--json"])
        .output()
        .expect("run installed nyx-uninstall");
    assert!(
        uninst_out.status.success(),
        "uninstall failed:\nstdout: {}\nstderr: {}",
        String::from_utf8_lossy(&uninst_out.stdout),
        String::from_utf8_lossy(&uninst_out.stderr)
    );
    let urep: InstallReport = serde_json::from_slice(&uninst_out.stdout).unwrap();
    assert!(urep.ok && urep.action == "uninstall");

    let t0 = Instant::now();
    while inst_root.exists() && t0.elapsed() < Duration::from_secs(6) {
        std::thread::sleep(Duration::from_millis(50));
    }
    assert!(
        !Path::new(&inst_root).exists(),
        "install root {} still exists after uninstall + reaper wait",
        inst_root.display()
    );
    assert_eq!(
        fs::read(docs_nyx.join("keep.txt")).unwrap(),
        b"must survive uninstall"
    );

    let _ = fs::remove_dir_all(&tmp_root);
}
