// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Compile-time build metadata for `nyx-engine` provenance (`RunRecord`, `results.h5`, and `vvuq_report`).
//!
//! Emits `NYX_BUILD_GIT_COMMIT`, `NYX_BUILD_GIT_DIRTY`, `NYX_BUILD_RUSTC_VERSION`, `NYX_BUILD_TARGET`,
//! `NYX_BUILD_PROFILE`, and `NYX_BUILD_CARGO_LOCK_BLAKE3` with existence-guarded `cargo:rerun-if-changed`
//! directives so missing VCS files never cause spurious rebuilds.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn sanitize(s: &str) -> String {
    s.replace(['\r', '\n'], " ").trim().to_string()
}

fn watch_if_exists(p: &Path) {
    if p.exists() {
        println!("cargo:rerun-if-changed={}", p.display());
    }
}

fn resolve_git_dir(ws_root: &Path) -> Option<PathBuf> {
    let dot_git = ws_root.join(".git");
    if dot_git.is_dir() {
        return Some(dot_git);
    }
    if dot_git.is_file() {
        let text = fs::read_to_string(&dot_git).ok()?;
        let rel = text.lines().find_map(|l| l.trim().strip_prefix("gitdir:"))?.trim();
        let path = PathBuf::from(rel);
        return Some(if path.is_absolute() { path } else { ws_root.join(path) });
    }
    None
}

fn resolve_jj_repo_dir(ws_root: &Path) -> Option<PathBuf> {
    let dot_jj_repo = ws_root.join(".jj").join("repo");
    if dot_jj_repo.is_dir() {
        return Some(dot_jj_repo);
    }
    if dot_jj_repo.is_file() {
        let text = fs::read_to_string(&dot_jj_repo).ok()?;
        let rel = text.trim();
        if rel.is_empty() {
            return None;
        }
        let path = PathBuf::from(rel);
        let base = ws_root.join(".jj");
        return Some(if path.is_absolute() { path } else { base.join(path) });
    }
    None
}

fn run_trimmed(program: &str, args: &[&str], cwd: &Path) -> Option<String> {
    let out = Command::new(program).args(args).current_dir(cwd).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = sanitize(&String::from_utf8_lossy(&out.stdout));
    (!s.is_empty()).then_some(s)
}

fn cargo_lock_blake3(ws_root: &Path) -> String {
    let lock = ws_root.join("Cargo.lock");
    if lock.is_file() {
        println!("cargo:rerun-if-changed={}", lock.display());
        if let Ok(bytes) = fs::read(&lock) {
            return format!("b3:{}", blake3::hash(&bytes).to_hex());
        }
    }
    "unknown".to_string()
}

fn watch_vcs(ws_root: &Path) {
    if let Some(git_dir) = resolve_git_dir(ws_root) {
        let head = git_dir.join("HEAD");
        watch_if_exists(&head);
        watch_if_exists(&git_dir.join("index"));
        if let Ok(text) = fs::read_to_string(&head) {
            if let Some(ref_rel) = text.trim().strip_prefix("ref:") {
                watch_if_exists(&git_dir.join(ref_rel.trim()));
            }
        }
    }
    let jj_checkout = ws_root.join(".jj").join("working_copy").join("checkout");
    watch_if_exists(&jj_checkout);
    if let Some(jj_repo) = resolve_jj_repo_dir(ws_root) {
        watch_if_exists(&jj_repo.join("op_heads").join("heads"));
    }
}

fn git_or_jj_provenance(ws_root: &Path) -> (String, String) {
    let git_commit = run_trimmed("git", &["rev-parse", "--short=12", "HEAD"], ws_root);
    let git_dirty = if git_commit.is_some() {
        Command::new("git")
            .args(["status", "--porcelain"])
            .current_dir(ws_root)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| if String::from_utf8_lossy(&o.stdout).trim().is_empty() { "false" } else { "true" })
    } else {
        None
    };

    if let Some(commit) = git_commit {
        let dirty = git_dirty.unwrap_or("false").to_string();
        return (commit, dirty);
    }

    let jj_commit =
        run_trimmed("jj", &["log", "--no-graph", "-r", "@", "-n", "1", "-T", "commit_id.short(12)"], ws_root);
    let jj_dirty = run_trimmed(
        "jj",
        &["log", "--no-graph", "-r", "@", "-n", "1", "-T", "if(empty, \"false\", \"true\")"],
        ws_root,
    );
    (jj_commit.unwrap_or_else(|| "unknown".to_string()), jj_dirty.unwrap_or_else(|| "false".to_string()))
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap_or_else(|_| ".".to_string());
    let manifest_path = Path::new(&manifest_dir);
    let ws_root = manifest_path.join("../..").canonicalize().unwrap_or_else(|_| manifest_path.to_path_buf());

    watch_vcs(&ws_root);
    let cargo_lock_hash = cargo_lock_blake3(&ws_root);
    let (git_commit, git_dirty) = git_or_jj_provenance(&ws_root);

    let rustc = env::var("RUSTC").unwrap_or_else(|_| "rustc".to_string());
    let rustc_version = run_trimmed(&rustc, &["--version"], &ws_root).unwrap_or_else(|| "unknown".to_string());
    let target = sanitize(&env::var("TARGET").unwrap_or_else(|_| "unknown".to_string()));
    let profile = sanitize(&env::var("PROFILE").unwrap_or_else(|_| "unknown".to_string()));

    println!("cargo:rustc-env=NYX_BUILD_GIT_COMMIT={git_commit}");
    println!("cargo:rustc-env=NYX_BUILD_GIT_DIRTY={git_dirty}");
    println!("cargo:rustc-env=NYX_BUILD_RUSTC_VERSION={rustc_version}");
    println!("cargo:rustc-env=NYX_BUILD_TARGET={target}");
    println!("cargo:rustc-env=NYX_BUILD_PROFILE={profile}");
    println!("cargo:rustc-env=NYX_BUILD_CARGO_LOCK_BLAKE3={cargo_lock_hash}");
}
