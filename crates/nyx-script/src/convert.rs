// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! JSON to and from the values of a script (S4). Written by hand: a script value is a `MontyObject`, an engine value is a `serde_json::Value`,
//! and a number keeps its kind (an integer stays an integer, everything else is a float).

use monty_types::{CallArgs, MontyObject, ObjectRef};
use serde_json::Value;

/// JSON to a script value.
pub fn to_object(v: &Value) -> MontyObject {
    match v {
        Value::Null => MontyObject::none(),
        Value::Bool(b) => MontyObject::bool(*b),
        Value::Number(n) => match n.as_i64() {
            Some(i) => MontyObject::int(i),
            None => MontyObject::float(n.as_f64().unwrap_or(f64::NAN)),
        },
        Value::String(s) => MontyObject::string(s.as_str()),
        Value::Array(a) => MontyObject::list(a.iter().map(to_object)),
        Value::Object(o) => MontyObject::dict(o.iter().map(|(k, x)| (MontyObject::string(k.as_str()), to_object(x)))),
    }
}

/// A script value to JSON; a value JSON cannot hold (a function, a set, a dict key that is not text, an integer past 64 bits) is an error.
pub fn from_object(o: ObjectRef) -> Result<Value, String> {
    match o.type_name() {
        "NoneType" => Ok(Value::Null),
        "bool" => Ok(Value::Bool(o.as_bool().unwrap_or(false))),
        "int" => o.as_int().map(Value::from).ok_or_else(|| "an integer larger than 64 bits cannot be sent to the engine (use a float)".to_string()),
        "float" => {
            let f = o.as_float().unwrap_or(f64::NAN);
            serde_json::Number::from_f64(f).map(Value::Number).ok_or_else(|| format!("{f} is not a JSON number"))
        }
        "str" => Ok(Value::String(o.as_str().unwrap_or_default().to_string())),
        "list" | "tuple" => o.items().unwrap_or_default().into_iter().map(from_object).collect::<Result<Vec<_>, _>>().map(Value::Array),
        "dict" => {
            let mut m = serde_json::Map::new();
            for (k, x) in o.pairs().unwrap_or_default() {
                let key = k.as_str().ok_or_else(|| format!("a dict key of type {} cannot be sent to the engine (keys are text)", k.type_name()))?;
                m.insert(key.to_string(), from_object(x)?);
            }
            Ok(Value::Object(m))
        }
        t => Err(format!("a value of type {t} cannot be sent to the engine")),
    }
}

/// The arguments of one command call as the engine's JSON object: the keyword arguments, or one positional dict. `skip` leading positional
/// arguments (the command name of `nyx.call`) are not part of it.
pub fn args_json(args: &CallArgs, skip: usize) -> Result<Value, String> {
    let kw: Vec<_> = args.kwargs().collect();
    let pos: Vec<_> = args.args().skip(skip).collect();
    match (pos.len(), kw.is_empty()) {
        (0, _) => {
            let mut m = serde_json::Map::new();
            for (k, x) in kw {
                m.insert(k.as_str().ok_or("a keyword name is not text")?.to_string(), from_object(x)?);
            }
            Ok(Value::Object(m))
        }
        (1, true) if pos[0].type_name() == "dict" => from_object(pos[0]),
        _ => Err("call a command with keyword arguments (nyx.set_param(name=\"w\", value=5e-7)) or with one dict".to_string()),
    }
}
