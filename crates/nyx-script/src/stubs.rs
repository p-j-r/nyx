//! The type stub of the `nyx` module (S4): `Engine::catalog()` written as Python declarations. Monty type-checks a script against it before the
//! script runs, so a misspelled command, a misspelled or missing argument fails in milliseconds with the line named, and an editor that reads the
//! same text (`nyx script stubs` writes it as `nyx.pyi`) offers the commands and their documentation while the script is written.
//!
//! Only names and required-ness are checked, not value types: the catalog's schemas say `number` for a field that also takes an expression
//! text (`"w/2"`), and a static check that refused a script the engine accepts would be worse than none. Values are `Any`; the engine checks them.

use nyx_engine::commands::CommandDesc;
use nyx_engine::Engine;
use serde_json::Value;

const KEYWORDS: &[&str] = &[
    "False", "None", "True", "and", "as", "assert", "async", "await", "break", "class", "continue", "def", "del", "elif", "else", "except", "finally", "for",
    "from", "global", "if", "import", "in", "is", "lambda", "nonlocal", "not", "or", "pass", "raise", "return", "try", "while", "with", "yield",
];

/// One docstring line must not break the triple quotes or the indentation.
fn doc_line(s: &str) -> String {
    s.replace("\"\"\"", "'''").replace('\\', "\\\\").replace(['\n', '\r'], " ")
}

fn method(c: &CommandDesc) -> String {
    let props = c.schema.get("properties").and_then(Value::as_object);
    let required: Vec<&str> = c.schema.get("required").and_then(Value::as_array).map(|r| r.iter().filter_map(Value::as_str).collect()).unwrap_or_default();
    let mut params = Vec::new();
    // a name that is a Python keyword cannot be written as `name=value`: the script uses `nyx.call("command", **{"from": 1})` and the stub accepts it
    let mut open = c.schema.get("additionalProperties").and_then(Value::as_bool).unwrap_or(false);
    let mut docs = Vec::new();
    for (key, def) in props.into_iter().flatten() {
        let d = def.get("description").and_then(Value::as_str).unwrap_or("");
        docs.push(format!("{key}: {}", doc_line(d)));
        if KEYWORDS.contains(&key.as_str()) || !key.chars().all(|ch| ch.is_alphanumeric() || ch == '_') {
            open = true;
        } else if required.contains(&key.as_str()) {
            params.push(format!("{key}: Any"));
        }
    }
    // required first, then the optional ones (a keyword-only parameter may follow another without a default, but keeping the order is clearer)
    for (key, _) in props.into_iter().flatten() {
        if !required.contains(&key.as_str()) && !KEYWORDS.contains(&key.as_str()) && key.chars().all(|ch| ch.is_alphanumeric() || ch == '_') {
            params.push(format!("{key}: Any = ..."));
        }
    }
    // `*` must be followed by a named parameter, so a command with none takes `**more` alone
    let sig = match (params.is_empty(), open) {
        (true, false) => String::new(),
        (true, true) => ", **more: Any".to_string(),
        (false, false) => format!(", *, {}", params.join(", ")),
        (false, true) => format!(", *, {}, **more: Any", params.join(", ")),
    };
    let mut s = format!("    def {}(self{sig}) -> Any:\n", c.name);
    s.push_str(&format!("        \"\"\"{}", doc_line(&c.title)));
    for d in docs {
        s.push_str(&format!("\n\n        {d}"));
    }
    s.push_str("\"\"\"\n        ...\n");
    s
}

/// The stub of `nyx` for the commands of `engine`: a class with one method per command and the global `nyx` of that class.
pub fn nyx_pyi(engine: &Engine) -> String {
    let mut s = String::from(
        "# nyx.pyi -- generated from Engine::catalog(); do not edit. `nyx` is the object every nyx script has: one method per engine command.\nfrom typing import Any\n\nclass _Nyx:\n",
    );
    s.push_str("    def call(self, name: str, args: dict[str, Any] = ..., /, **kwargs: Any) -> Any:\n        \"\"\"Any command by name, with its arguments as keywords or as one dict.\"\"\"\n        ...\n");
    for c in engine.catalog() {
        s.push_str(&method(&c));
    }
    s.push_str("\nnyx: _Nyx\n");
    s.push_str(&crate::trace::stub());
    s
}
