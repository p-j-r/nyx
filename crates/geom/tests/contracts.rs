//! Interface contracts of geom (docs/architecture/adr/0002 and 0003).
//!
//! * `geom_depends_on_serde_only`: the dependency direction `app -> femcore -> geom` is kept by this test, not by convention.
//! * `golden_feature_lists_round_trip_key_for_key`: a saved feature list read and written again is the same JSON value. `goldens.rs` already
//!   checks that the list equals the builder's; this also catches a field that serde would silently drop or default (ROUGH_SPOTS R20).
//! * `feature_list_top_level_keys`: the keys saved files carry today. There is no format-version key yet (ADR 19 proposes `format`).
use geom::model::FeatureList;
use serde_json::Value;
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn geom_depends_on_serde_only() {
    let text = std::fs::read_to_string(root().join("Cargo.toml")).unwrap();
    let mut in_deps = false;
    let mut deps = vec![];
    for line in text.lines() {
        let l = line.trim();
        if l.starts_with('[') {
            in_deps = l == "[dependencies]";
            continue;
        }
        if in_deps && !l.is_empty() && !l.starts_with('#') {
            if let Some((name, _)) = l.split_once('=') {
                deps.push(name.trim().to_string());
            }
        }
    }
    deps.sort();
    assert_eq!(deps, vec!["rayon".to_string(), "serde".to_string()], "geom may depend on serde and (optionally) rayon only; no femcore, faer or GUI (ADR 3)");
    assert!(text.contains("rayon = { version = \"1\", optional = true }"), "rayon must stay optional");
}

fn golden_lists() -> Vec<(String, String)> {
    let mut out = vec![];
    for e in std::fs::read_dir(root().join("golden")).unwrap() {
        let p = e.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().to_string();
        if name.ends_with(".json") && !name.ends_with(".expected.json") {
            out.push((name, std::fs::read_to_string(&p).unwrap()));
        }
    }
    out.sort();
    out
}

#[test]
fn golden_feature_lists_round_trip_key_for_key() {
    let lists = golden_lists();
    assert!(lists.len() >= 20, "found {} golden feature lists", lists.len());
    for (name, text) in lists {
        let original: Value = serde_json::from_str(&text).unwrap();
        let list: FeatureList = serde_json::from_str(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
        let again: Value = serde_json::to_value(&list).unwrap();
        assert_eq!(original, again, "{name}: a key was dropped, added or defaulted by the round trip; that is a change of the saved format");
    }
}

#[test]
fn feature_list_top_level_keys() {
    for (name, text) in golden_lists() {
        let v: Value = serde_json::from_str(&text).unwrap();
        let mut keys: Vec<&str> = v.as_object().unwrap().keys().map(|s| s.as_str()).collect();
        keys.sort();
        for k in ["dim", "features", "params"] {
            assert!(keys.contains(&k), "{name}: missing top-level key {k} in {keys:?}");
        }
    }
}

#[test]
fn every_feature_carries_the_fields_that_nyx_store_and_the_gui_read() {
    for (name, text) in golden_lists() {
        let v: Value = serde_json::from_str(&text).unwrap();
        for f in v["features"].as_array().unwrap() {
            for k in ["id", "name", "op", "suppressed", "origin", "tags"] {
                assert!(f.get(k).is_some(), "{name}: feature without `{k}`");
            }
            assert!(f["op"].as_object().map(|o| o.len() == 1).unwrap_or(f["op"].is_string()), "{name}: `op` is an externally tagged enum");
        }
    }
}
