//! The path guard (productivity plan S4, ADR [89]): with a root set, no command reads or writes a file outside it, and no command of the catalog
//! takes a path-like argument the guard does not know about.

use nyx_engine::path_guard::PATH_ARGS;
use nyx_engine::Engine;
use serde_json::json;
use std::path::PathBuf;

fn temp_root(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-guard-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A command added to the catalog with an argument that names a file or a folder must be listed in `PATH_ARGS`, or a script could use it to leave the root.
#[test]
fn every_path_like_argument_of_the_catalog_is_guarded() {
    let e = Engine::new();
    let mut unguarded = Vec::new();
    for c in e.catalog() {
        let Some(props) = c.schema.get("properties").and_then(|p| p.as_object()) else { continue };
        for (key, def) in props {
            let doc = def.get("description").and_then(|d| d.as_str()).unwrap_or("").to_lowercase();
            // whole words of the name (`profile` is not a file), and words in the description that say a file or a folder is meant
            let words: Vec<String> = key.to_lowercase().split('_').map(str::to_string).collect();
            let pathlike = ["path", "file", "folder", "dir", "directory", "filename", "workspace"]
                .iter()
                .any(|w| words.iter().any(|k| k == w))
                || doc.contains("file path")
                || doc.contains("folder")
                || doc.contains("directory");
            if pathlike && !PATH_ARGS.contains(&key.as_str()) {
                unguarded.push(format!("{}.{key} ({doc})", c.name));
            }
        }
    }
    assert!(
        unguarded.is_empty(),
        "path-like arguments the guard does not know: {unguarded:?}; add the name to PATH_ARGS"
    );
    // and the listed names are really used
    let named: Vec<String> = e
        .catalog()
        .iter()
        .flat_map(|c| {
            c.schema
                .get("properties")
                .and_then(|p| p.as_object())
                .map(|p| p.keys().cloned().collect::<Vec<_>>())
                .unwrap_or_default()
        })
        .collect();
    for k in PATH_ARGS {
        assert!(named.iter().any(|n| n == k), "PATH_ARGS names `{k}` but no command takes it");
    }
}

#[test]
fn a_command_cannot_write_outside_the_root() {
    let root = temp_root("write");
    let outside = temp_root("write-outside");
    let mut e = Engine::new();
    e.call("project_new", &json!({ "dim": 2, "title": "guard" })).unwrap();
    e.set_path_root(Some(&root)).unwrap();

    // a relative path lands in the root
    e.call("project_save", &json!({ "path": "model.nyx.json" })).unwrap();
    assert!(root.join("model.nyx.json").exists());

    // out by .., by an absolute path
    let up = e.call("project_save", &json!({ "path": "../escaped.nyx.json" })).unwrap_err().to_string();
    assert!(up.contains("project_save") && up.contains("path"), "{up}");
    let target = outside.join("abs.nyx.json");
    assert!(e.call("project_save", &json!({ "path": target.to_str().unwrap() })).is_err());
    assert!(!target.exists(), "the command ran although the path was refused");
    assert!(!root.parent().unwrap().join("escaped.nyx.json").exists());

    // a folder argument is guarded the same way
    assert!(e.call("store_open", &json!({ "folder": outside.to_str().unwrap() })).is_err());

    // reading: a file outside the root is not opened
    std::fs::write(outside.join("secret.nyx.json"), "{}").unwrap();
    assert!(e.call("project_open", &json!({ "path": outside.join("secret.nyx.json").to_str().unwrap() })).is_err());

    // and with the guard off the same command works
    e.set_path_root(None).unwrap();
    e.call("project_save", &json!({ "path": target.to_str().unwrap() })).unwrap();
    assert!(target.exists());

    let _ = std::fs::remove_dir_all(&root);
    let _ = std::fs::remove_dir_all(&outside);
}
