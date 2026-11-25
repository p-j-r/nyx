// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unified `bolero` fuzz and regression tests for `path_guard::confine`, `material_fetch` validation, and `nyx_store::archive` ZIP checks.

use bolero::TypeGenerator;
use nyx_engine::facade::store::archive::{check_name, import, inspect, Limits};
use nyx_engine::material_fetch::{validate_material_name, validate_remote_url};
use nyx_engine::path_guard::confine;
use std::io::{Cursor, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, ZipWriter};

static SEQ: AtomicUsize = AtomicUsize::new(0);

const WINDOWS_RESERVED: &[&str] = &[
    "con", "prn", "aux", "nul", "com1", "com2", "com3", "com4", "com5", "com6", "com7", "com8", "com9", "lpt1", "lpt2",
    "lpt3", "lpt4", "lpt5", "lpt6", "lpt7", "lpt8", "lpt9",
];

fn shared_root() -> &'static (PathBuf, PathBuf) {
    static ROOT: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    ROOT.get_or_init(|| {
        let dir = std::env::temp_dir().join(format!("nyx_fuzz_path_guard_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create root");
        let canon = dir.canonicalize().expect("canonicalize root");
        (canon.clone(), canon)
    })
}

fn fresh_temp_dir(tag: &str) -> PathBuf {
    let id = SEQ.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("nyx_fuzz_{tag}_{}_{id}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn build_zip(entries: &[(&str, &[u8])]) -> Vec<u8> {
    let mut zw = ZipWriter::new(Cursor::new(Vec::new()));
    let opts = SimpleFileOptions::default().compression_method(CompressionMethod::Deflated);
    for (name, data) in entries {
        if zw.start_file(*name, opts).is_ok() {
            let _ = zw.write_all(data);
        }
    }
    zw.finish().map(|c| c.into_inner()).unwrap_or_default()
}

fn exercise_archive_bytes(bytes: &[u8], lim: &Limits) {
    let tmp = fresh_temp_dir("arch");
    let nyx_path = tmp.join("sample.nyx");
    let inbox = tmp.join("inbox");
    std::fs::create_dir_all(&inbox).expect("create inbox");
    std::fs::write(&nyx_path, bytes).expect("write archive");
    let _ = inspect(&nyx_path, lim);
    if let Ok(wf) = import(&nyx_path, &inbox, lim) {
        let canon_inbox = inbox.canonicalize().expect("canon inbox");
        let canon_wf = wf.root.canonicalize().expect("canon wf");
        assert!(canon_wf.starts_with(&canon_inbox), "imported workflow escaped inbox: {canon_wf:?}");
    }
    let _ = std::fs::remove_dir_all(&tmp);
}

fn is_windows_reserved_stem(comp: &str) -> bool {
    let stem = comp.split('.').next().unwrap_or("");
    let lower = stem.to_ascii_lowercase();
    WINDOWS_RESERVED.contains(&lower.as_str())
}

#[derive(Debug, TypeGenerator)]
struct CraftedArchive {
    entry_name: String,
    payload: Vec<u8>,
    include_fake_manifest: bool,
    extra_entries: u8,
}

#[test]
fn fuzz_path_guard_and_material_validation_no_panic() {
    let (root, canon_root) = shared_root();
    bolero::check!().with_iterations(512).for_each(|bytes: &[u8]| {
        let s = String::from_utf8_lossy(bytes);
        let s: &str = match s.char_indices().nth(512) {
            Some((idx, _)) => &s[..idx],
            None => &s,
        };

        if let Ok(resolved) = confine(root, s) {
            assert!(
                resolved.starts_with(canon_root),
                "confined path {resolved:?} escaped root {canon_root:?} for input {s:?}"
            );
        }

        if let Ok(name) = validate_material_name(s) {
            assert!(!name.is_empty() && name.len() <= 64);
            assert!(name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-'));
            assert!(!is_windows_reserved_stem(name));
        }

        let prefixed = format!("https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/{s}");
        for cand_url in [s, prefixed.as_str()] {
            if let Ok(v) = validate_remote_url(cand_url) {
                assert!(v.url.starts_with("https://"));
                let rest = &v.url["https://".len()..];
                let authority = rest.split(['/', '?', '#']).next().unwrap_or("");
                assert!(!authority.is_empty() && !authority.contains('@') && !authority.contains(':'));
                assert!(!authority.contains('\\'));
                if v.host == "raw.githubusercontent.com" {
                    assert!(!v.url.contains('\\'));
                }
                assert!(!v.host.is_empty());
            }
        }

        if check_name(s).is_ok() {
            assert!(!s.is_empty() && s.len() <= 240);
            assert!(!s.contains('\\') && !s.contains(':') && !s.contains('\0'));
            assert!(!s.starts_with('/'));
            for part in s.split('/') {
                assert!(part != "." && part != ".." && !part.is_empty());
                assert!(!is_windows_reserved_stem(part));
            }
        }
    });
}

#[test]
fn fuzz_archive_raw_and_crafted_zips_no_panic() {
    let tight = Limits { max_entries: 16, max_file: 16 * 1024, max_total: 32 * 1024, max_ratio: 100 };

    bolero::check!().with_iterations(256).for_each(|bytes: &[u8]| {
        exercise_archive_bytes(bytes, &tight);
    });

    bolero::check!().with_type::<CraftedArchive>().with_iterations(256).for_each(|case: &CraftedArchive| {
        let safe_name: String = case.entry_name.chars().take(120).collect();
        let mut entries: Vec<(String, Vec<u8>)> = vec![(safe_name, case.payload.iter().copied().take(4096).collect())];
        if case.include_fake_manifest {
            entries.push((
                "nyx-archive.json".into(),
                br#"{"format":1,"kind":"workflow","level":"Recipe","nyx_version":"0.0.0","files":[]}"#.to_vec(),
            ));
        }
        for i in 0..usize::from(case.extra_entries % 20) {
            entries.push((format!("extra_{i}.txt"), b"payload".to_vec()));
        }
        let refs: Vec<(&str, &[u8])> = entries.iter().map(|(n, b)| (n.as_str(), b.as_slice())).collect();
        let zip_bytes = build_zip(&refs);
        exercise_archive_bytes(&zip_bytes, &tight);
    });
}

#[test]
fn regression_hostile_paths_urls_reserved_names_and_zip_slip() {
    let (root, canon_root) = shared_root();
    let outside = fresh_temp_dir("outside");

    // Traversal, absolute paths, outside dirs, NUL, and empty paths rejected; anything accepted stays under root.
    for hostile in [
        "../outside.txt",
        "sub/../../outside.txt",
        "/etc/passwd",
        outside.join("secret.txt").to_str().unwrap(),
        "bad\0name",
        "",
    ] {
        assert!(confine(root, hostile).is_err(), "confine must reject {hostile:?}");
    }
    #[cfg(unix)]
    {
        let link_path = root.join("escape_link");
        let _ = std::fs::remove_file(&link_path);
        std::os::unix::fs::symlink(&outside, &link_path).expect("create escape symlink");
        assert!(confine(root, "escape_link/secret.txt").is_err(), "confine must reject symlink pointing outside root");
        let _ = std::fs::remove_file(&link_path);
    }
    for win_or_weird in ["C:relative.nyx", "D:\\abs\\file", r"\\server\share\file", "stream:hidden.txt", "CON"] {
        if let Ok(p) = confine(root, win_or_weird) {
            assert!(p.starts_with(canon_root), "accepted path {p:?} escaped {canon_root:?}");
        }
    }
    let ok_path = confine(root, "nested/good_file.txt").expect("valid confined path");
    assert!(ok_path.starts_with(canon_root));

    // Material name + URL validation (cross-platform Windows reserved names & colon/traversal bans).
    for bad_name in ["", "con", "PRN", "AUX", "nul", "Com1", "lpt9", "a/b", "a\\b", "a:b", "..", &"x".repeat(65)] {
        assert!(validate_material_name(bad_name).is_err(), "must reject material name {bad_name:?}");
    }
    assert_eq!(validate_material_name("SiO2_thin-01").unwrap(), "SiO2_thin-01");

    for bad_url in [
        "http://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/Si.yml",
        "https://evil.example.com/Si.toml",
        "https://raw.githubusercontent.com/p-j-r/nyx-materials/main/Si.toml",
        "https://user@raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/Si.yml",
        "https://raw.githubusercontent.com:443/polyanskiy/refractiveindex.info-database/main/Si.yml",
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/../../secret.toml",
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/%2e%2e/secret.toml",
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main\\Si.yml",
        "https://raw.githubusercontent.com/org/repo/main/Si.toml?token=1",
    ] {
        assert!(validate_remote_url(bad_url).is_err(), "must reject URL {bad_url:?}");
    }
    let ok_url = validate_remote_url(
        "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/master/database/data/main/Si/nk/Aspnes.yml",
    )
    .expect("valid allowlisted curated HTTPS URL");
    assert_eq!(ok_url.host, "raw.githubusercontent.com");
    let ok_ri = validate_remote_url("https://refractiveindex.info/database/data/main/Si/nk/Aspnes.yml")
        .expect("valid allowlisted host URL");
    assert_eq!(ok_ri.host, "refractiveindex.info");

    // Archive entry name + ZIP slip / Windows reserved names / invalid archive rejection.
    for bad_entry in
        ["../escape", "e1/../../escape", "C:rel", "C:\\x", "a:b", "e1/CON.txt", "e1/aux", "e1/LPT1.yaml", "/abs", ""]
    {
        assert!(check_name(bad_entry).is_err(), "check_name must reject {bad_entry:?}");
    }

    let tmp = fresh_temp_dir("hostile_zip");
    let inbox = tmp.join("inbox");
    std::fs::create_dir_all(&inbox).unwrap();
    let lim = Limits::default();

    let zip_slip = build_zip(&[
        ("nyx-archive.json", br#"{"format":1,"kind":"workflow","level":"Recipe","nyx_version":"0.1.0","files":[]}"#),
        ("../../evil.txt", b"owned"),
    ]);
    let p_slip = tmp.join("slip.nyx");
    std::fs::write(&p_slip, &zip_slip).unwrap();
    assert!(import(&p_slip, &inbox, &lim).is_err());

    let p_bad = tmp.join("corrupt.nyx");
    std::fs::write(&p_bad, b"PK\x03\x04not_a_real_zip_central_dir").unwrap();
    assert!(import(&p_bad, &inbox, &lim).is_err());
    assert!(!Path::new(&tmp.join("evil.txt")).exists());

    let _ = std::fs::remove_dir_all(&outside);
    let _ = std::fs::remove_dir_all(&tmp);
}
