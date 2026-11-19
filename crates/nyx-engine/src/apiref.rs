//! The API reference (productivity plan S2): one Markdown page per engine command and an index, generated from `Engine::catalog()` and from
//! nothing else, so the reference cannot disagree with the commands (ADR 57: the catalog is the one list). A page has the command's title,
//! description and access, a table of its arguments (flattened from the JSON Schema: type, whether it is required, the unit or meaning the schema
//! states, the default) and an example line in the format of a command list (`nyx check`, ADR 89). The example is made from the schema alone
//! (required arguments with placeholder values) and is checked to pass `nyx check`. `tests/api_docs.rs` regenerates the pages and compares them
//! with `docs/api/`; `nyx api list` and `nyx api show` print the same text.

use crate::commands::{Access, CommandDesc};
use crate::engine::Engine;
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;

/// What an access level means, in words.
pub fn access_words(a: Access) -> &'static str {
    match a {
        Access::Read => "reads the project or the results",
        Access::Edit => "changes the project",
        Access::Expensive => "starts a mesh build or a solve, or writes files",
    }
}

fn access_name(a: Access) -> &'static str {
    match a {
        Access::Read => "read",
        Access::Edit => "edit",
        Access::Expensive => "expensive",
    }
}

/// The commands whose names are closest to `name` (at most three, within an edit distance of three): the hint of a mistyped command.
pub fn nearest(names: &[String], name: &str) -> Vec<String> {
    let mut near: Vec<(usize, &String)> =
        names.iter().map(|n| (edit_distance(n, name), n)).filter(|(d, _)| *d <= 3).collect();
    near.sort();
    near.into_iter().take(3).map(|(_, n)| n.clone()).collect()
}

fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for i in 1..=a.len() {
        let mut prev = row[0];
        row[0] = i;
        for j in 1..=b.len() {
            let keep = row[j];
            row[j] = (row[j] + 1).min(row[j - 1] + 1).min(prev + usize::from(a[i - 1] != b[j - 1]));
            prev = keep;
        }
    }
    row[b.len()]
}

// ------------------------------------------------------------------------------------------------------------------ arguments

/// One argument of a command, as the page's table shows it.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    /// the name; a field of an object argument is `outer.inner`, an item of a list `outer[].inner`
    pub name: String,
    pub ty: String,
    pub required: bool,
    /// the title and description the schema gives: where the unit or the meaning is stated
    pub doc: String,
    pub default: Option<String>,
}

fn types_of(s: &Value) -> Vec<String> {
    if let Some(alts) = s.get("oneOf").and_then(Value::as_array) {
        let mut out: Vec<String> = vec![];
        for a in alts {
            for t in types_of(a) {
                if !out.contains(&t) {
                    out.push(t);
                }
            }
        }
        return out;
    }
    match s.get("type") {
        Some(Value::String(t)) => vec![type_name(t, s)],
        Some(Value::Array(a)) => a.iter().filter_map(Value::as_str).map(|t| type_name(t, s)).collect(),
        _ if s.get("enum").is_some() => vec!["text".into()],
        _ => vec!["any".into()],
    }
}

fn type_name(t: &str, s: &Value) -> String {
    match t {
        "string" => "text".into(),
        "integer" => "integer".into(),
        "number" => "number".into(),
        "boolean" => "true or false".into(),
        "object" => "object".into(),
        "null" => "null".into(),
        "array" => {
            let item = s.get("items").map(|i| types_of(i).join(" or ")).unwrap_or_else(|| "any".into());
            match (s.get("minItems").and_then(Value::as_u64), s.get("maxItems").and_then(Value::as_u64)) {
                (Some(a), Some(b)) if a == b => format!("list of {a} ({item})"),
                _ => format!("list of {item}"),
            }
        }
        other => other.to_string(),
    }
}

fn plain(v: &Value) -> String {
    v.as_str().map(String::from).unwrap_or_else(|| v.to_string())
}

fn doc_of(s: &Value) -> String {
    let mut parts: Vec<String> = vec![];
    for k in ["title", "description"] {
        if let Some(t) = s.get(k).and_then(Value::as_str).filter(|t| !t.is_empty()) {
            parts.push(t.to_string());
        }
    }
    // what each item of a list is, when the list itself says nothing more than its shape
    if let Some(item) = s.get("items").and_then(|i| i.get("description")).and_then(Value::as_str) {
        parts.push(format!("each item: {item}"));
    }
    parts.join(". ")
}

fn flatten(schema: &Value, prefix: &str, depth: usize, out: &mut Vec<Row>) {
    let Some(props) = schema.get("properties").and_then(Value::as_object) else { return };
    let required: Vec<&str> =
        schema.get("required").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str).collect();
    // sorted by name: a schema map keeps insertion order in a build where any crate enables `preserve_order` and sorts otherwise (ADR 26)
    let mut props: Vec<(&String, &Value)> = props.iter().collect();
    props.sort_by_key(|(n, _)| n.as_str());
    for (name, s) in props {
        let full = format!("{prefix}{name}");
        let mut ty = types_of(s).join(" or ");
        if let Some(allowed) = s.get("enum").and_then(Value::as_array) {
            let shown: Vec<String> = allowed.iter().map(|a| format!("`{}`", plain(a))).collect();
            ty = format!("{ty}: {}", shown.join(", "));
        }
        out.push(Row {
            name: full.clone(),
            ty,
            required: required.contains(&name.as_str()),
            doc: doc_of(s),
            default: s.get("default").map(plain),
        });
        if depth < 3 {
            if s.get("properties").is_some() {
                flatten(s, &format!("{full}."), depth + 1, out);
            }
            if let Some(items) = s.get("items").filter(|i| i.get("properties").is_some()) {
                flatten(items, &format!("{full}[]."), depth + 1, out);
            }
        }
    }
}

/// The arguments of a command, flattened, in the order of the schema (sorted by name: the schema is a map).
pub fn rows(c: &CommandDesc) -> Vec<Row> {
    let mut out = vec![];
    flatten(&c.schema, "", 0, &mut out);
    out
}

// ------------------------------------------------------------------------------------------------------------------ example

fn placeholder(name: &str, s: &Value, depth: usize) -> Value {
    if let Some(d) = s.get("default") {
        return d.clone();
    }
    if let Some(alts) = s.get("oneOf").and_then(Value::as_array) {
        return alts.first().map(|a| placeholder(name, a, depth)).unwrap_or(Value::Null);
    }
    if let Some(e) = s.get("enum").and_then(Value::as_array).and_then(|e| e.first()) {
        return e.clone();
    }
    let t = match s.get("type") {
        Some(Value::String(t)) => t.as_str(),
        Some(Value::Array(a)) => a.iter().filter_map(Value::as_str).find(|t| *t != "null").unwrap_or("null"),
        _ => "",
    };
    match t {
        "string" => json!(format!("<{name}>")),
        "integer" => json!(s.get("minimum").and_then(Value::as_i64).unwrap_or(0).max(0)),
        "number" => json!(s.get("minimum").and_then(Value::as_f64).unwrap_or(0.0).max(0.0)),
        "boolean" => json!(false),
        "array" => {
            let n = s.get("minItems").and_then(Value::as_u64).unwrap_or(0) as usize;
            let item = s.get("items").cloned().unwrap_or(Value::Null);
            Value::Array((0..n).map(|_| placeholder(name, &item, depth + 1)).collect())
        }
        "object" => object_example(s, depth + 1),
        _ => Value::Null,
    }
}

fn object_example(schema: &Value, depth: usize) -> Value {
    let mut o = Map::new();
    if depth > 3 {
        return Value::Object(o);
    }
    let props = schema.get("properties").and_then(Value::as_object);
    for r in schema.get("required").and_then(Value::as_array).into_iter().flatten().filter_map(Value::as_str) {
        let s = props.and_then(|p| p.get(r)).cloned().unwrap_or(Value::Null);
        o.insert(r.to_string(), placeholder(r, &s, depth));
    }
    Value::Object(o)
}

/// The `args` of an example call: every required argument with a placeholder of its type (`<name>` for text, the first choice of an enumeration,
/// zero for a number, a list of the length the schema asks for).
pub fn example_args(c: &CommandDesc) -> Value {
    object_example(&c.schema, 0)
}

/// The example as one line of a command list.
pub fn example_line(c: &CommandDesc) -> String {
    // keys in a fixed order, whatever features the build has
    let args = nyx_physics::canonical_json(&example_args(c));
    format!("{{\"cmd\":{},\"args\":{args}}}", json!(c.name))
}

// ------------------------------------------------------------------------------------------------------------------ pages

fn cell(s: &str) -> String {
    s.replace('|', "\\|").replace(['\n', '\r'], " ")
}

/// The Markdown page of one command.
pub fn command_page(c: &CommandDesc) -> String {
    let mut s = String::new();
    s.push_str(&format!("# `{}`\n\n", c.name));
    s.push_str(&format!("**{}** · access `{}`: {}\n\n", c.title, access_name(c.access), access_words(c.access)));
    if !c.description.trim().is_empty() {
        s.push_str(c.description.trim());
        s.push_str("\n\n");
    }
    s.push_str("## Arguments\n\n");
    let rows = rows(c);
    if rows.is_empty() {
        s.push_str("None.\n\n");
    } else {
        s.push_str("| Argument | Type | Required | Default | What it is |\n|---|---|---|---|---|\n");
        for r in &rows {
            s.push_str(&format!(
                "| `{}` | {} | {} | {} | {} |\n",
                r.name,
                cell(&r.ty),
                if r.required { "yes" } else { "no" },
                r.default.as_deref().map(|d| format!("`{}`", cell(d))).unwrap_or_default(),
                cell(&r.doc)
            ));
        }
        s.push('\n');
    }
    s.push_str("Lengths are in metres and every unit is SI unless the argument says otherwise.\n\n");
    s.push_str("## Example\n\nA line of a command list (`nyx check`, `nyx run`), made from the schema with placeholder values (`<...>` is text to replace):\n\n```json\n");
    s.push_str(&example_line(c));
    s.push_str("\n```\n");
    s
}

/// The index page: every command with its access and title, and the hash of the catalog the pages were made from.
pub fn index_page(cmds: &[CommandDesc]) -> String {
    let hash = nyx_store::util::hash_bytes(
        nyx_physics::canonical_json(&serde_json::to_value(cmds).unwrap_or(Value::Null)).to_string().as_bytes(),
    );
    let mut s = String::new();
    s.push_str("# nyx engine commands\n\n");
    s.push_str("Generated from `Engine::catalog()` by `nyx_engine::apiref` and checked by `crates/nyx-engine/tests/api_docs.rs`: **do not edit**. After a deliberate change of a command, rewrite with `NYX_WRITE_API=1 cargo test --release -p nyx-engine --test api_docs` (and the interface version, `engine_commands` in `docs/architecture/contracts/VERSIONS.md`, changes with it).\n\n");
    s.push_str(&format!(
        "Catalog hash: `{hash}` · {} commands. `nyx api list` and `nyx api show <command>` print the same text.\n\n",
        cmds.len()
    ));
    s.push_str("| Command | Access | What it does |\n|---|---|---|\n");
    for c in cmds {
        s.push_str(&format!("| [`{0}`]({0}.md) | {1} | {2} |\n", c.name, access_name(c.access), cell(&c.title)));
    }
    s
}

/// Every page of the reference by file name (`<command>.md`, `README.md`), sorted.
pub fn pages(e: &Engine) -> BTreeMap<String, String> {
    let mut cmds = e.catalog();
    cmds.sort_by(|a, b| a.name.cmp(&b.name));
    let mut out = BTreeMap::new();
    for c in &cmds {
        out.insert(format!("{}.md", c.name), command_page(c));
    }
    out.insert("README.md".into(), index_page(&cmds));
    out
}
