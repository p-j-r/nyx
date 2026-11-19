//! The API reference (productivity plan S2, ADR 57): `docs/api/` is generated from `Engine::catalog()`, one page per command and an index, and
//! must be exactly what the generator writes now, so no page can disagree with a command and no command can lack a page. Every example on a
//! page is a line that `nyx check` accepts.
//!
//! Rewrite after a deliberate change of a command (the interface version in `docs/architecture/contracts/VERSIONS.md` changes with it):
//! `NYX_WRITE_API=1 cargo test --release -p nyx-engine --test api_docs`

use nyx_engine::apiref;
use nyx_engine::commandlist;
use nyx_engine::Engine;
use std::collections::BTreeSet;
use std::path::PathBuf;

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/api")
}

fn lf(s: &str) -> String {
    s.replace("\r\n", "\n")
}

#[test]
fn the_reference_in_docs_api_is_what_the_catalog_generates() {
    let pages = apiref::pages(&Engine::new());
    let dir = dir();
    if std::env::var_os("NYX_WRITE_API").is_some() {
        std::fs::create_dir_all(&dir).unwrap();
        for (name, text) in &pages {
            std::fs::write(dir.join(name), text).unwrap();
        }
        // a page of a command that no longer exists goes
        for f in std::fs::read_dir(&dir).unwrap().map(|f| f.unwrap().path()) {
            let name = f.file_name().unwrap().to_string_lossy().into_owned();
            if name.ends_with(".md") && !pages.contains_key(&name) {
                std::fs::remove_file(f).unwrap();
            }
        }
        return;
    }
    let mut wrong = vec![];
    for (name, text) in &pages {
        match std::fs::read_to_string(dir.join(name)) {
            Ok(have) if lf(&have) == *text => {}
            Ok(_) => wrong.push(format!("{name}: differs from what the catalog generates")),
            Err(_) => wrong.push(format!("{name}: missing")),
        }
    }
    let have: BTreeSet<String> = std::fs::read_dir(&dir)
        .expect("docs/api")
        .map(|f| f.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".md"))
        .collect();
    for n in have {
        if !pages.contains_key(&n) {
            wrong.push(format!("{n}: no such command in the catalog"));
        }
    }
    assert!(
        wrong.is_empty(),
        "the API reference is out of date; regenerate with NYX_WRITE_API=1 cargo test --release -p nyx-engine --test api_docs\n{}",
        wrong.join("\n")
    );
}

#[test]
fn every_command_has_a_page_the_index_lists_and_an_example_that_passes_the_check() {
    let e = Engine::new();
    let cmds = e.catalog();
    let pages = apiref::pages(&e);
    assert_eq!(pages.len(), cmds.len() + 1, "one page per command and the index");
    let index = &pages["README.md"];
    for c in &cmds {
        assert!(pages.contains_key(&format!("{}.md", c.name)), "{}", c.name);
        assert!(index.contains(&format!("[`{0}`]({0}.md)", c.name)), "the index lacks {}", c.name);
        let page = &pages[&format!("{}.md", c.name)];
        assert!(page.contains("## Arguments") && page.contains("## Example"), "{}", c.name);
        // the example is a list line that the static check accepts
        let (lines, bad) = commandlist::parse(&apiref::example_line(c));
        assert!(bad.is_empty() && lines.len() == 1, "{}: {bad:?}", c.name);
        let problems = e.check_static(&lines);
        assert!(problems.is_empty(), "{}: its example does not pass `nyx check`: {:?}", c.name, problems);
    }
}

#[test]
fn a_page_shows_the_arguments_their_types_and_where_the_unit_is_stated() {
    let e = Engine::new();
    let pages = apiref::pages(&e);
    let rect = &pages["add_rect.md"];
    assert!(rect.contains("| `min` | list of 2 (number or text) | yes |"), "{rect}");
    assert!(rect.contains("lower-left corner (x, z) in metres"), "the unit the schema states is on the page");
    let phys = &pages["describe_physics.md"];
    assert!(phys.contains("`dfs`, `ems`, `sms`"), "an enumeration lists its choices: {phys}");
    let sweep = &pages["set_mesh.md"];
    assert!(sweep.contains("`sweep[].thickness`"), "the fields of a list of objects are rows too");
    assert!(pages["project_info.md"].contains("None."), "a command with no arguments says so");
    assert!(pages["solve_run.md"].contains("access `expensive`"));
}

#[test]
fn the_nearest_names_are_the_hint_of_a_mistyped_command() {
    let names: Vec<String> = Engine::new().catalog().into_iter().map(|c| c.name).collect();
    assert_eq!(apiref::nearest(&names, "solve_runn")[0], "solve_run");
    assert!(apiref::nearest(&names, "zzzzzzzzzzzzzzzz").is_empty());
}

#[test]
fn every_argument_of_every_command_says_what_it_is() {
    // the reference is only as good as the schemas: an argument with no title or description is a gap in the catalog, found here
    // (the fields of a list of objects are covered by the description of the list when it has one)
    let mut bare = vec![];
    for c in Engine::new().catalog() {
        for r in apiref::rows(&c) {
            if r.doc.trim().is_empty() && !r.name.contains('.') && !r.name.contains('[') {
                bare.push(format!("{}.{}", c.name, r.name));
            }
        }
    }
    assert!(bare.is_empty(), "arguments with no description: {bare:?}");
}

#[test]
fn the_pages_do_not_depend_on_the_features_of_the_build() {
    // a schema keeps the order of insertion in a build where a crate enables `preserve_order`; the pages sort, so the same text comes out
    let e = Engine::new();
    let a = apiref::pages(&e);
    let mut cmds = e.catalog();
    for c in &mut cmds {
        // reverse every object of the schema: the text of the page must not change
        fn rev(v: &mut serde_json::Value) {
            match v {
                serde_json::Value::Object(m) => {
                    let mut pairs: Vec<(String, serde_json::Value)> = std::mem::take(m).into_iter().collect();
                    pairs.reverse();
                    for (k, mut x) in pairs {
                        rev(&mut x);
                        m.insert(k, x);
                    }
                }
                serde_json::Value::Array(a) => a.iter_mut().for_each(rev),
                _ => {}
            }
        }
        rev(&mut c.schema);
        assert_eq!(apiref::command_page(c), a[&format!("{}.md", c.name)], "{}", c.name);
    }
}
