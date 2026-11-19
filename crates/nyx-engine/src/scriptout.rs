//! Export as script (productivity plan S7b, ADR [94]): the project as a script that rebuilds it. A person who built a model in the window
//! gives it to an agent (or to themselves, as text) and the agent runs it with `params` to change a value, so a hand-built model is the
//! start of a script and not a JSON blob.
//!
//! The engine builds a list of `(command, arguments)`, replays it on a scratch engine and offers it only when the result is the project; then it
//! renders the list as Python. Three forms, the first that verifies: `commands` (`project_new`, `set_param`, one `add_<kind>` per feature,
//! `set_mesh`, `add_material`, `set_problem`, `set_active`), `geometry_json` (the same with the geometry as one `set_geometry`: exact ids) and
//! `project_json` (one `project_load` of the whole project: always exact, not readable).
//!
//! "Is the project" for `commands` means: the same description of the geometry (`geometry::describe`: names instead of ids, so a model with
//! removed features and id gaps replays to the same model, and the origin of a feature, `Gui` or `Code`, is not compared: it changes no result)
//! and the same project in everything else. The other two forms must be exactly equal.

use nyx_physics::{canonical_json, Error};
use serde_json::{json, Map, Value};

use crate::geomdesc;
use crate::geometry;
use crate::project::{MeshSettings, Project};
use crate::Engine;

type Cmd = (String, Value);

/// The script and how it was made.
pub struct Exported {
    pub script: String,
    /// `commands`, `geometry_json` or `project_json`
    pub form: &'static str,
    pub commands: usize,
}

/// The project as a script that rebuilds it (it does not solve).
pub fn export(p: &Project) -> Result<Exported, Error> {
    let tiers: [(&'static str, Result<Vec<Cmd>, Error>); 3] = [
        ("commands", build(p, false)),
        ("geometry_json", build(p, true)),
        ("project_json", Ok(vec![("project_load".to_string(), json!({ "json": p.to_json() }))])),
    ];
    for (form, cmds) in tiers {
        let Ok(cmds) = cmds else { continue };
        if replays_to(p, &cmds, form == "commands") {
            return Ok(Exported { script: render(p, form, &cmds), form, commands: cmds.len() });
        }
    }
    Err(Error::invalid("this project cannot be written as a script: no form of it replays to the same project"))
}

fn mesh_args(m: &MeshSettings) -> Value {
    let d = MeshSettings::default();
    let mut a = Map::new();
    let mut put = |k: &str, cur: Value, def: Value| {
        if cur != def {
            a.insert(k.to_string(), cur);
        }
    };
    put("h", json!(m.h), json!(d.h));
    put("adapt_scale", json!(m.adapt_scale), json!(d.adapt_scale));
    put("physics_sizing", json!(m.physics_sizing), json!(d.physics_sizing));
    put("points_per_wavelength", json!(m.points_per_wavelength), json!(d.points_per_wavelength));
    put("structured_pml", json!(m.structured_pml), json!(d.structured_pml));
    put("max_dev", json!(m.max_dev), json!(d.max_dev));
    put("sweep", json!(m.sweep), json!(d.sweep));
    Value::Object(a)
}

/// The commands that rebuild `p` (the geometry as `set_geometry` when `geometry_json`).
fn build(p: &Project, geometry_json: bool) -> Result<Vec<Cmd>, Error> {
    let mut c: Vec<Cmd> = vec![("project_new".into(), json!({ "title": p.title, "dim": p.geometry.dim }))];
    let desc = geometry::describe(&p.geometry);
    if geometry_json {
        c.push((
            "set_geometry".into(),
            json!({ "geometry": serde_json::to_value(&p.geometry).map_err(|e| Error::invalid(e.to_string()))? }),
        ));
    } else {
        for q in desc["params"].as_array().into_iter().flatten() {
            c.push(("set_param".into(), json!({ "name": q["name"], "value": q["value"] })));
        }
        for f in desc["features"].as_array().into_iter().flatten() {
            let kind = f["kind"].as_str().ok_or_else(|| Error::invalid("a feature without a kind"))?;
            let mut args = f["args"].as_object().cloned().unwrap_or_default();
            args.insert("name".into(), f["name"].clone());
            c.push((geomdesc::tool_name(kind), Value::Object(args)));
            if f["suppressed"].as_bool() == Some(true) {
                c.push(("suppress_feature".into(), json!({ "name": f["name"] })));
            }
            for t in f["tags"].as_array().into_iter().flatten() {
                c.push(("tag_feature".into(), json!({ "name": f["name"], "tag": t })));
            }
        }
    }
    let mesh = mesh_args(&p.mesh);
    if mesh.as_object().is_some_and(|m| !m.is_empty()) {
        c.push(("set_mesh".into(), mesh));
    }
    for (name, toml) in &p.materials {
        c.push(("add_material".into(), json!({ "name": name, "toml": toml })));
    }
    for s in &p.physics {
        c.push(("set_problem".into(), json!({ "physics": s.physics, "problem": serde_json::to_value(s).map_err(|e| Error::invalid(e.to_string()))? })));
    }
    if let Some(a) = &p.active {
        c.push(("set_active".into(), json!({ "physics": a })));
    }
    Ok(c)
}

/// Replays `cmds` on a scratch engine; true when what comes out is `p` (see the module comment for what that means).
fn replays_to(p: &Project, cmds: &[Cmd], described: bool) -> bool {
    let mut e = Engine::new();
    for (cmd, args) in cmds {
        if e.call(cmd, args).is_err() {
            return false;
        }
    }
    let q = e.project();
    if described {
        let rest = |x: &Project| {
            let mut v = serde_json::to_value(x).expect("a Project serialises");
            v.as_object_mut().map(|o| o.remove("geometry"));
            canonical_json(&v)
        };
        geometry::describe(&q.geometry) == geometry::describe(&p.geometry) && rest(q) == rest(p)
    } else {
        canonical_json(&serde_json::to_value(q).expect("a Project serialises"))
            == canonical_json(&serde_json::to_value(p).expect("a Project serialises"))
    }
}

// ---- rendering ----

/// Python keywords: an argument of such a name cannot be written `name=value`.
const KEYWORDS: &[&str] = &[
    "False", "None", "True", "and", "as", "assert", "async", "await", "break", "class", "continue", "def", "del",
    "elif", "else", "except", "finally", "for", "from", "global", "if", "import", "in", "is", "lambda", "nonlocal",
    "not", "or", "pass", "raise", "return", "try", "while", "with", "yield",
];

fn is_ident(s: &str) -> bool {
    !s.is_empty()
        && !s.starts_with(|c: char| c.is_ascii_digit())
        && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
        && !KEYWORDS.contains(&s)
}

fn text_literal(s: &str) -> String {
    // a long text of lines (a material file) reads better between triple quotes; anything that could end or escape them stays a JSON string
    if s.contains('\n')
        && !s.contains("\"\"\"")
        && !s.contains('\\')
        && !s.ends_with('"')
        && !s.chars().any(|c| c.is_control() && c != '\n')
    {
        return format!("\"\"\"{s}\"\"\"");
    }
    serde_json::to_string(s).expect("a string prints")
}

/// A JSON value as Python source. A dict or list that does not fit a line is laid out over lines, `indent` spaces in.
fn literal(v: &Value, indent: usize) -> String {
    let inline = match v {
        Value::Null => return "None".into(),
        Value::Bool(b) => return if *b { "True".into() } else { "False".into() },
        Value::Number(n) => {
            // an integer stays an integer; a float is written with the shortest text that reads back to the same value
            return if let Some(i) = n.as_i64() {
                i.to_string()
            } else if let Some(u) = n.as_u64() {
                u.to_string()
            } else {
                format!("{:?}", n.as_f64().unwrap_or(f64::NAN))
            };
        }
        Value::String(s) => return text_literal(s),
        Value::Array(a) => format!("[{}]", a.iter().map(|x| literal(x, 0)).collect::<Vec<_>>().join(", ")),
        Value::Object(o) => format!(
            "{{{}}}",
            sorted(o)
                .iter()
                .map(|(k, x)| format!("{}: {}", serde_json::to_string(k).expect("a key prints"), literal(x, 0)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    };
    if inline.len() + indent <= 110 && !inline.contains('\n') {
        return inline;
    }
    let pad = " ".repeat(indent + 4);
    match v {
        Value::Array(a) => format!(
            "[\n{}\n{}]",
            a.iter().map(|x| format!("{pad}{},", literal(x, indent + 4))).collect::<Vec<_>>().join("\n"),
            " ".repeat(indent)
        ),
        Value::Object(o) => format!(
            "{{\n{}\n{}}}",
            sorted(o)
                .iter()
                .map(|(k, x)| format!(
                    "{pad}{}: {},",
                    serde_json::to_string(k).expect("a key prints"),
                    literal(x, indent + 4)
                ))
                .collect::<Vec<_>>()
                .join("\n"),
            " ".repeat(indent)
        ),
        _ => inline,
    }
}

fn sorted(o: &Map<String, Value>) -> Vec<(&String, &Value)> {
    let mut v: Vec<_> = o.iter().collect();
    v.sort_by(|a, b| a.0.cmp(b.0));
    v
}

/// The names of the arguments in the order a reader expects: `name` first, then the descriptor's order for a geometry operation, else alphabetical.
fn arg_order(cmd: &str, args: &Map<String, Value>) -> Vec<String> {
    let mut keys: Vec<String> = args.keys().cloned().collect();
    keys.sort();
    let fields: Vec<&'static str> = geomdesc::all()
        .into_iter()
        .find(|d| geomdesc::tool_name(d.kind) == cmd)
        .map(|d| d.fields.iter().map(|f| f.name).collect())
        .unwrap_or_default();
    keys.sort_by_key(|k| (k != "name", fields.iter().position(|f| f == k).unwrap_or(usize::MAX)));
    keys
}

fn call_line(cmd: &str, args: &Value) -> String {
    let empty = Map::new();
    let a = args.as_object().unwrap_or(&empty);
    if cmd == "set_param" {
        // every parameter reads `params`, so `script_run params={"w": 3e-7}` is the variant with nothing edited
        let name = a.get("name").and_then(Value::as_str).unwrap_or_default();
        return format!(
            "nyx.set_param(name={}, value=params.get({}, {}))",
            serde_json::to_string(name).expect("a name prints"),
            serde_json::to_string(name).expect("a name prints"),
            literal(&a["value"], 0)
        );
    }
    let order = arg_order(cmd, a);
    if order.iter().all(|k| is_ident(k)) {
        let body = order.iter().map(|k| format!("{k}={}", literal(&a[k], 4))).collect::<Vec<_>>().join(", ");
        if body.len() + cmd.len() + 6 <= 110 && !body.contains('\n') {
            return format!("nyx.{cmd}({body})");
        }
        let inner = order.iter().map(|k| format!("    {k}={},", literal(&a[k], 4))).collect::<Vec<_>>().join("\n");
        return format!("nyx.{cmd}(\n{inner}\n)");
    }
    format!("nyx.call({}, {})", serde_json::to_string(cmd).expect("a name prints"), literal(args, 0))
}

fn render(p: &Project, form: &str, cmds: &[Cmd]) -> String {
    let title: String = p.title.chars().map(|c| if c.is_control() { ' ' } else { c }).collect();
    let mut s = format!(
        "# {title} -- exported from nyx (project_script, ADR 94; form: {form}).\n# It rebuilds the project command by command and does not solve: add nyx.solve_run() where you want one.\n# Over MCP run it with script_run allow_destructive true (it starts with project_new, which discards the open project).\n# Each model parameter reads `params`: script_run params={{\"w\": 3e-7}} gives that variant with nothing edited (under `nyx run` params is empty and the values are the ones below).\n"
    );
    if form == "project_json" {
        s.push_str("# This model could not be written command by command, so it is one project_load of the whole project (exact, not readable).\n");
    }
    if form == "geometry_json" {
        s.push_str("# The geometry is one set_geometry (its overlay pins or ids must be kept exactly), the rest is written command by command.\n");
    }
    s.push('\n');
    for (cmd, args) in cmds {
        s.push_str(&call_line(cmd, args));
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small model with a parameter, an id gap, a suppressed and a tagged feature, a mesh option.
    fn model() -> Engine {
        let mut e = Engine::new();
        let mut c = |cmd: &str, args: Value| {
            e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}: {er}"));
        };
        c("project_new", json!({ "title": "fallbacks", "dim": 2 }));
        c("set_param", json!({ "name": "w", "value": 5e-7 }));
        c("add_rect", json!({ "name": "box", "min": [-2e-6, -2e-6], "max": [2e-6, 0.0] }));
        c("add_rect", json!({ "name": "junk", "min": [0.0, 0.0], "max": [1e-7, 1e-7] }));
        c("add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", 2.2e-7] }));
        c("remove_feature", json!({ "name": "junk" }));
        c("tag_feature", json!({ "name": "core", "tag": "wire" }));
        c("set_mesh", json!({ "h": 2e-7 }));
        e
    }

    /// The two fallback forms are not reached by any real project, so each is driven here directly: it replays to the project, exactly, and
    /// renders as a script that says which form it is.
    #[test]
    fn the_fallback_forms_replay_exactly_and_say_which_they_are() {
        let e = model();
        let p = e.project();
        let geometry = build(p, true).unwrap();
        assert!(geometry.iter().any(|(c, _)| c == "set_geometry") && !geometry.iter().any(|(c, _)| c == "add_rect"));
        assert!(replays_to(p, &geometry, false), "geometry_json replays to the exact project (ids kept)");
        let s = render(p, "geometry_json", &geometry);
        assert!(s.contains("form: geometry_json") && s.contains("nyx.set_geometry("), "{s}");

        let whole = vec![("project_load".to_string(), json!({ "json": p.to_json() }))];
        assert!(replays_to(p, &whole, false), "project_json replays to the exact project");
        let s = render(p, "project_json", &whole);
        assert!(
            s.contains("form: project_json") && s.contains("nyx.project_load(") && s.contains("not readable"),
            "{s}"
        );
    }

    /// The exactness check is real: a list that rebuilds something else is refused, and the readable form is the one chosen when it holds.
    #[test]
    fn a_wrong_command_list_is_refused_and_the_readable_form_is_chosen() {
        let e = model();
        let p = e.project();
        let mut wrong = build(p, false).unwrap();
        wrong.retain(|(c, _)| c != "tag_feature");
        assert!(!replays_to(p, &wrong, true), "a list without the tag is not the project");
        let r = export(p).unwrap();
        assert_eq!(r.form, "commands");
    }
}
