// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Supply-chain & dependency hygiene tests for `nyx-install` (CLI default feature set + offline SPDX notices).

use nyx_install::bundle::target_triple_guess;
use nyx_install::notices;
use std::path::PathBuf;
use std::process::Command;

#[test]
fn default_installer_crate_has_zero_c_sys_or_wgpu_dependencies() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());

    // Normal runtime dependency tree must contain ZERO native C *-sys wrappers, C builders, or GUI/wgpu crates.
    // (`windows-sys` is pure-Rust Win32 FFI bindings pulled unconditionally on Windows by `zip >= 9`.)
    let normal_out = Command::new(&cargo)
        .current_dir(&manifest_dir)
        .args(["tree", "--offline", "-e", "normal", "-p", "nyx-install"])
        .output()
        .expect("run cargo tree -e normal");
    assert!(normal_out.status.success());
    let normal_tree = String::from_utf8_lossy(&normal_out.stdout);
    let has_c_sys = normal_tree
        .lines()
        .any(|line| line.contains("-sys v") && !line.contains("windows-sys v"));
    assert!(
        !has_c_sys,
        "forbidden runtime dep '-sys v':\n{normal_tree}"
    );
    for forbidden in [
        " cc ",
        " cmake ",
        " bindgen ",
        " pkg-config ",
        " vcpkg ",
        "wgpu v",
        "eframe v",
    ] {
        assert!(
            !normal_tree.contains(forbidden),
            "forbidden runtime dep '{forbidden}':\n{normal_tree}"
        );
    }

    // Build + normal tree must also forbid native C library wrappers (allowing embed-resource's vswhom-sys/windows-sys & blake3 cc).
    let out = Command::new(&cargo)
        .current_dir(&manifest_dir)
        .args([
            "tree",
            "--offline",
            "-e",
            "normal,build",
            "-p",
            "nyx-install",
        ])
        .output()
        .expect("run cargo tree");
    assert!(
        out.status.success(),
        "cargo tree failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let tree = String::from_utf8_lossy(&out.stdout);
    for forbidden in [
        " cmake ",
        " bindgen ",
        " pkg-config ",
        " vcpkg ",
        "hdf5-sys",
        "libhdf5-sys",
        "libz-sys",
        "zlib-ng-sys",
        "openssl-sys",
        "wgpu v",
        "eframe v",
    ] {
        assert!(
            !tree.contains(forbidden),
            "forbidden dep '{forbidden}' in default nyx-install:\n{tree}"
        );
    }
    assert!(
        tree.contains("zlib-rs "),
        "expected pure-Rust zlib-rs in default nyx-install:\n{tree}"
    );
    assert!(
        tree.contains("blake3 "),
        "expected blake3 in nyx-install:\n{tree}"
    );
}

#[test]
fn offline_spdx_notice_generator_covers_lockfile_and_expression_operators() {
    assert!(notices::validate_spdx_expression("Apache-2.0 OR GPL-2.0-only").is_ok());
    assert!(notices::validate_spdx_expression(
        "(MIT OR Apache-2.0) AND OFL-1.1 AND Ubuntu-font-1.0"
    )
    .is_ok());
    assert!(notices::validate_spdx_expression("Apache-2.0 WITH LLVM-exception").is_ok());

    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap();
    let triple = target_triple_guess();
    let txt = notices::generate_third_party_notices(&repo_root, Some(triple))
        .expect("generate THIRD_PARTY_NOTICES.txt offline");
    assert!(txt.contains("AGPL-3.0-only OR LicenseRef-nyx-Commercial"));
    assert!(txt.contains(notices::CORRESPONDING_SOURCE_URL));
    assert!(txt.contains("blake3 "));
    assert!(txt.contains("serde "));
}
