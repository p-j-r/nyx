//! Typed parameters of the tools: one table per tool that produces the JSON Schema an AI host sees, validates the arguments, fills defaults
//! and records every assumed value in a ledger. The same table could drive a GUI form or a script API; nothing here knows about MCP.

use serde_json::{json, Map, Value};

/// How much a parameter moves the answer. High ones are asked for when missing, the others are assumed and listed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Matters {
    Low,
    Medium,
    High,
}

impl Matters {
    pub fn as_str(self) -> &'static str {
        match self {
            Matters::Low => "low",
            Matters::Medium => "medium",
            Matters::High => "high",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Ty {
    Num { min: f64, max: f64 },
    Int { min: i64, max: i64 },
    Bool,
    Enum(&'static [&'static str]),
    /// a handle of the given kind letter (`g`, `m`, `p`, `r`)
    Handle(char),
    /// a list of strings, each from the set
    List(&'static [&'static str]),
}

#[derive(Clone, Copy, Debug)]
pub enum Policy {
    Required,
    Optional,
    /// assumed when missing, and said so in the result
    Default { json: &'static str, why: &'static str },
    /// asked for when missing (the result is `needs_input`); `use_defaults: true` accepts the fallback instead
    Ask { json: &'static str, why: &'static str, suggest: &'static str },
}

#[derive(Clone, Copy, Debug)]
pub struct Param {
    pub name: &'static str,
    pub ty: Ty,
    pub unit: &'static str,
    pub doc: &'static str,
    pub policy: Policy,
    pub matters: Matters,
}

/// What a tool may do to the world, as reported to the host (MCP tool annotations).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Access {
    /// reads session state only
    Read,
    /// creates artifacts in the session folder, never touches user files
    SessionWrite,
    /// as `SessionWrite`, and may run for seconds to minutes
    Expensive,
}

#[derive(Clone, Copy, Debug)]
pub struct ToolDef {
    pub name: &'static str,
    pub title: &'static str,
    pub description: &'static str,
    pub params: &'static [Param],
    pub access: Access,
}

/// A value the user did not give, recorded so the final answer can say what was assumed.
#[derive(Clone, Debug, PartialEq)]
pub struct Assumed {
    pub name: String,
    pub value: Value,
    pub why: String,
    pub matters: &'static str,
}

impl Assumed {
    pub fn to_json(&self) -> Value {
        json!({ "name": self.name, "value": self.value, "why": self.why, "matters": self.matters })
    }
}

/// An error the AI can act on: a short code, a message, and the data needed to fix the call (valid names, ranges).
#[derive(Clone, Debug)]
pub struct ToolError {
    pub code: &'static str,
    pub message: String,
    pub details: Map<String, Value>,
}

impl ToolError {
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        ToolError { code, message: message.into(), details: Map::new() }
    }
    pub fn with(mut self, key: &str, value: impl Into<Value>) -> Self {
        self.details.insert(key.to_string(), value.into());
        self
    }
    pub fn to_json(&self) -> Value {
        let mut m = Map::new();
        m.insert("error".into(), json!(self.code));
        m.insert("message".into(), json!(self.message));
        for (k, v) in &self.details {
            m.insert(k.clone(), v.clone());
        }
        Value::Object(m)
    }
}

/// Why a tool did not return a result: a real error, or a question for the user.
#[derive(Clone, Debug)]
pub enum Stop {
    Err(ToolError),
    NeedsInput(Value),
}

impl From<ToolError> for Stop {
    fn from(e: ToolError) -> Self {
        Stop::Err(e)
    }
}

pub type ToolResult = Result<Value, Stop>;

/// Arguments after validation: every parameter has its value or is absent only if optional.
#[derive(Clone, Debug)]
pub struct Bound {
    vals: Map<String, Value>,
    pub assumed: Vec<Assumed>,
}

impl Bound {
    pub fn num(&self, name: &str) -> f64 {
        self.vals.get(name).and_then(Value::as_f64).unwrap_or(f64::NAN)
    }
    pub fn opt_num(&self, name: &str) -> Option<f64> {
        self.vals.get(name).and_then(Value::as_f64)
    }
    pub fn int(&self, name: &str) -> i64 {
        self.vals.get(name).and_then(Value::as_i64).unwrap_or(0)
    }
    pub fn opt_int(&self, name: &str) -> Option<i64> {
        self.vals.get(name).and_then(Value::as_i64)
    }
    pub fn boolean(&self, name: &str) -> bool {
        self.vals.get(name).and_then(Value::as_bool).unwrap_or(false)
    }
    pub fn str(&self, name: &str) -> &str {
        self.vals.get(name).and_then(Value::as_str).unwrap_or("")
    }
    pub fn opt_str(&self, name: &str) -> Option<&str> {
        self.vals.get(name).and_then(Value::as_str)
    }
    pub fn list(&self, name: &str) -> Vec<String> {
        self.vals.get(name).and_then(Value::as_array).map(|a| a.iter().filter_map(|v| v.as_str().map(String::from)).collect()).unwrap_or_default()
    }
    pub fn get(&self, name: &str) -> Option<&Value> {
        self.vals.get(name)
    }
}

fn unit_note(p: &Param) -> String {
    if p.unit.is_empty() {
        String::new()
    } else {
        format!(" [{}]", p.unit)
    }
}

impl Param {
    fn property(&self) -> Value {
        let mut doc = format!("{}{}", self.doc, unit_note(self));
        match self.policy {
            Policy::Default { json, why } => doc.push_str(&format!(" Default {json} ({why}); matters: {}.", self.matters.as_str())),
            Policy::Ask { json, why, .. } => doc.push_str(&format!(" Asked if missing (matters: {}); fallback {json} ({why}) with use_defaults.", self.matters.as_str())),
            _ => {}
        }
        let mut p = match self.ty {
            Ty::Num { min, max } => json!({"type": "number", "minimum": min, "maximum": max}),
            Ty::Int { min, max } => json!({"type": "integer", "minimum": min, "maximum": max}),
            Ty::Bool => json!({"type": "boolean"}),
            Ty::Enum(v) => json!({"type": "string", "enum": v}),
            Ty::Handle(k) => json!({"type": "string", "pattern": format!("^{k}_[0-9a-f]+$")}),
            Ty::List(v) => json!({"type": "array", "items": {"type": "string", "enum": v}}),
        };
        p["description"] = json!(doc);
        p
    }
}

impl ToolDef {
    pub fn input_schema(&self) -> Value {
        let mut props = Map::new();
        let mut required = vec![];
        let mut asks = false;
        for p in self.params {
            props.insert(p.name.to_string(), p.property());
            match p.policy {
                Policy::Required => required.push(p.name),
                Policy::Ask { .. } => asks = true,
                _ => {}
            }
        }
        if asks {
            props.insert("use_defaults".into(), json!({"type": "boolean", "description": "Accept the fallback for every parameter that would otherwise be asked from the user."}));
        }
        json!({"type": "object", "properties": props, "required": required, "additionalProperties": false})
    }

    pub fn annotations(&self) -> Value {
        json!({
            "title": self.title,
            "readOnlyHint": self.access == Access::Read,
            "destructiveHint": false,
            "idempotentHint": true,
            "openWorldHint": false
        })
    }

    /// Checks the arguments against the table, fills defaults (recording them) and returns the bound values.
    pub fn bind(&self, args: &Value) -> Result<Bound, Stop> {
        let empty = Map::new();
        let given = match args {
            Value::Null => &empty,
            Value::Object(m) => m,
            _ => return Err(ToolError::new("bad_arguments", "arguments must be a JSON object").into()),
        };
        let valid: Vec<&str> = self.params.iter().map(|p| p.name).collect();
        let asks = self.params.iter().any(|p| matches!(p.policy, Policy::Ask { .. }));
        for k in given.keys() {
            if !valid.contains(&k.as_str()) && !(asks && k == "use_defaults") {
                return Err(ToolError::new("unknown_param", format!("'{k}' is not a parameter of {}", self.name)).with("got", k.as_str()).with("valid", valid).into());
            }
        }
        let use_defaults = given.get("use_defaults").and_then(Value::as_bool).unwrap_or(false);
        let mut vals = Map::new();
        let mut assumed = vec![];
        for p in self.params {
            match given.get(p.name) {
                Some(Value::Null) | None => match p.policy {
                    Policy::Required => {
                        return Err(ToolError::new("missing_param", format!("'{}' is required by {}", p.name, self.name)).with("param", p.name).with("doc", p.doc).into());
                    }
                    Policy::Optional => {}
                    Policy::Default { json, why } => {
                        let v: Value = serde_json::from_str(json).expect("default is valid JSON");
                        assumed.push(Assumed { name: p.name.into(), value: v.clone(), why: why.into(), matters: p.matters.as_str() });
                        vals.insert(p.name.into(), v);
                    }
                    Policy::Ask { json, why, suggest } => {
                        let v: Value = serde_json::from_str(json).expect("default is valid JSON");
                        if use_defaults {
                            assumed.push(Assumed { name: p.name.into(), value: v.clone(), why: why.into(), matters: p.matters.as_str() });
                            vals.insert(p.name.into(), v);
                        } else {
                            return Err(Stop::NeedsInput(json!({
                                "status": "needs_input",
                                "tool": self.name,
                                "slot": p.name,
                                "doc": p.doc,
                                "unit": p.unit,
                                "why": why,
                                "matters": p.matters.as_str(),
                                "suggest": serde_json::from_str::<Value>(suggest).unwrap_or(Value::Null),
                                "next": format!("Ask the user for '{}', then call {} again with it; or pass use_defaults: true to accept {json}.", p.name, self.name)
                            })));
                        }
                    }
                },
                Some(v) => {
                    vals.insert(p.name.into(), check_value(p, v)?);
                }
            }
        }
        Ok(Bound { vals, assumed })
    }
}

fn check_value(p: &Param, v: &Value) -> Result<Value, Stop> {
    let wrong = |want: &str| -> Stop { ToolError::new("wrong_type", format!("'{}' must be {want}", p.name)).with("param", p.name).with("got", v.clone()).into() };
    match p.ty {
        Ty::Num { min, max } => {
            let x = v.as_f64().ok_or_else(|| wrong("a number"))?;
            if !x.is_finite() || x < min || x > max {
                return Err(ToolError::new("out_of_range", format!("'{}' = {x} is outside [{min}, {max}]{}", p.name, unit_note(p)))
                    .with("param", p.name)
                    .with("got", x)
                    .with("min", min)
                    .with("max", max)
                    .with("unit", p.unit)
                    .into());
            }
            Ok(json!(x))
        }
        Ty::Int { min, max } => {
            let x = v.as_i64().ok_or_else(|| wrong("an integer"))?;
            if x < min || x > max {
                return Err(ToolError::new("out_of_range", format!("'{}' = {x} is outside [{min}, {max}]", p.name)).with("param", p.name).with("min", min).with("max", max).into());
            }
            Ok(json!(x))
        }
        Ty::Bool => Ok(json!(v.as_bool().ok_or_else(|| wrong("true or false"))?)),
        Ty::Enum(set) => {
            let s = v.as_str().ok_or_else(|| wrong("a string"))?;
            if !set.contains(&s) {
                return Err(ToolError::new("bad_value", format!("'{s}' is not a valid {}", p.name)).with("param", p.name).with("valid", set).into());
            }
            Ok(json!(s))
        }
        Ty::Handle(k) => {
            let s = v.as_str().ok_or_else(|| wrong("a handle string"))?;
            let ok = s.len() > 2 && s.starts_with(k) && s.as_bytes()[1] == b'_' && s[2..].chars().all(|c| c.is_ascii_hexdigit());
            if !ok {
                return Err(ToolError::new("bad_handle", format!("'{s}' is not a handle of kind {k}_")).with("param", p.name).into());
            }
            Ok(json!(s))
        }
        Ty::List(set) => {
            let a = v.as_array().ok_or_else(|| wrong("a list of strings"))?;
            let mut out = vec![];
            for it in a {
                let s = it.as_str().ok_or_else(|| wrong("a list of strings"))?;
                if !set.contains(&s) {
                    return Err(ToolError::new("bad_value", format!("'{s}' is not a valid entry of {}", p.name)).with("param", p.name).with("valid", set).into());
                }
                out.push(json!(s));
            }
            Ok(Value::Array(out))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: &[Param] = &[
        Param { name: "t", ty: Ty::Num { min: 1e-8, max: 1e-5 }, unit: "m", doc: "thickness", policy: Policy::Required, matters: Matters::High },
        Param { name: "w", ty: Ty::Num { min: 1e-8, max: 1e-4 }, unit: "m", doc: "width", policy: Policy::Ask { json: "5e-7", why: "typical", suggest: "[4.5e-7,5e-7]" }, matters: Matters::High },
        Param { name: "clad", ty: Ty::Enum(&["air", "oxide"]), unit: "", doc: "cladding", policy: Policy::Default { json: "\"air\"", why: "assumed" }, matters: Matters::Medium },
    ];
    const T: ToolDef = ToolDef { name: "t", title: "T", description: "d", params: P, access: Access::Read };

    #[test]
    fn missing_ask_param_is_a_question_and_use_defaults_accepts_the_fallback() {
        assert!(matches!(T.bind(&json!({"t": 3e-7})), Err(Stop::NeedsInput(_))));
        let b = T.bind(&json!({"t": 3e-7, "use_defaults": true})).unwrap();
        assert_eq!(b.num("w"), 5e-7);
        assert_eq!(b.str("clad"), "air");
        assert_eq!(b.assumed.iter().map(|a| a.name.as_str()).collect::<Vec<_>>(), ["w", "clad"]);
    }

    #[test]
    fn errors_name_the_valid_choices() {
        let Err(Stop::Err(e)) = T.bind(&json!({"t": 3e-7, "w": 5e-7, "clad": "glass"})) else { panic!() };
        assert_eq!(e.code, "bad_value");
        assert_eq!(e.details["valid"], json!(["air", "oxide"]));
        let Err(Stop::Err(e)) = T.bind(&json!({"t": 3e-7, "w": 5e-7, "colour": 1})) else { panic!() };
        assert_eq!(e.code, "unknown_param");
        let Err(Stop::Err(e)) = T.bind(&json!({"t": 3.0, "w": 5e-7})) else { panic!() };
        assert_eq!(e.code, "out_of_range");
        let Err(Stop::Err(e)) = T.bind(&json!({"w": 5e-7})) else { panic!() };
        assert_eq!(e.code, "missing_param");
    }

    #[test]
    fn the_schema_lists_required_parameters_and_units() {
        let s = T.input_schema();
        assert_eq!(s["required"], json!(["t"]));
        assert!(s["properties"]["t"]["description"].as_str().unwrap().contains("[m]"));
        assert!(s["properties"]["use_defaults"].is_object());
        assert_eq!(s["additionalProperties"], json!(false));
    }
}
