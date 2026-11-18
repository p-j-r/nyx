//! The Rhai binding (S4 spike): `nyx::<command>(#{...})` for every command of the catalog, and `nyx::call(name, map)`.
//!
//! An engine error is thrown as a Rhai error, so `try`/`catch` recovers. The host sets operation, depth and size limits and registers no
//! file, network, clock or random function (the `no_time` feature removes the clock).

use std::cell::{Cell, RefCell};
use std::rc::Rc;

use nyx_engine::Engine;
use rhai::{Dynamic, EvalAltResult, Map, Module, Position};
use serde_json::Value;

use crate::{Limits, Outcome};

type Shared<'a> = Rc<RefCell<&'a mut Engine>>;

fn fail(msg: String) -> Box<EvalAltResult> {
    EvalAltResult::ErrorRuntime(Dynamic::from(msg), Position::NONE).into()
}

fn send(engine: &Shared, calls: &Rc<Cell<u64>>, name: &str, args: Dynamic) -> Result<Dynamic, Box<EvalAltResult>> {
    let args: Value = from_script(&args).map_err(|e| fail(format!("{name}: arguments are not JSON: {e}")))?;
    calls.set(calls.get() + 1);
    let reply = engine.borrow_mut().call(name, &args).map_err(|e| fail(e.to_string()))?;
    Ok(to_script(&reply))
}

fn build<'a>(engine: &'a mut Engine, limits: &Limits, printed: &Rc<RefCell<Vec<String>>>) -> (rhai::Engine, Rc<Cell<u64>>) {
    let names: Vec<String> = engine.catalog().into_iter().map(|c| c.name).collect();
    let shared: Shared<'a> = Rc::new(RefCell::new(engine));
    let calls = Rc::new(Cell::new(0u64));
    // the closures hold `shared` for as long as the Rhai engine lives; the caller drops the Rhai engine before the borrow of `engine` ends
    let shared: Rc<RefCell<&'static mut Engine>> = unsafe { std::mem::transmute(shared) };
    let mut module = Module::new();
    for n in names {
        let (e, c, nm) = (shared.clone(), calls.clone(), n.clone());
        module.set_native_fn(n, move |args: Map| send(&e, &c, &nm, Dynamic::from_map(args)));
    }
    let (e, c) = (shared.clone(), calls.clone());
    module.set_native_fn("call", move |name: &str, args: Map| send(&e, &c, name, Dynamic::from_map(args)));
    let mut rh = rhai::Engine::new();
    rh.register_static_module("nyx", module.into());
    rh.register_fn("fmt", |x: f64, d: i64| crate::fmt_fixed(x, d));
    rh.register_fn("fmt", |x: i64, d: i64| crate::fmt_fixed(x as f64, d));
    rh.set_max_operations(limits.steps);
    rh.set_max_call_levels(limits.depth);
    rh.set_max_expr_depths(limits.depth * 2, limits.depth);
    // Rhai has no heap limit: strings, arrays and maps are bounded instead
    rh.set_max_string_size(limits.heap_bytes / 16);
    rh.set_max_array_size(limits.heap_bytes / 64);
    rh.set_max_map_size(limits.heap_bytes / 128);
    let p = printed.clone();
    rh.on_print(move |s| p.borrow_mut().push(s.to_string()));
    let p = printed.clone();
    rh.on_debug(move |s, _, _| p.borrow_mut().push(s.to_string()));
    (rh, calls)
}

/// Parses the script and nothing else.
pub fn parse(src: &str) -> Result<(), String> {
    rhai::Engine::new().compile(src).map(|_| ()).map_err(|e| e.to_string())
}

/// Runs a script. Its value is the script's last expression.
pub fn run(engine: &mut Engine, _name: &str, src: &str, limits: &Limits) -> Result<Outcome, String> {
    let printed = Rc::new(RefCell::new(Vec::new()));
    let (rh, calls) = build(engine, limits, &printed);
    let result = rh.eval::<Dynamic>(src);
    drop(rh);
    let v = result.map_err(|e| e.to_string())?;
    let value: Value = from_script(&v)?;
    let printed = printed.borrow().clone();
    Ok(Outcome { value, printed, calls: calls.get() })
}

/// JSON to a script value without serde: serde_json's `arbitrary_precision` (turned on for the whole build by `starlark`) makes serde write numbers as maps.
pub fn to_script(v: &Value) -> Dynamic {
    match v {
        Value::Null => Dynamic::UNIT,
        Value::Bool(b) => Dynamic::from(*b),
        Value::Number(n) => match n.as_i64() {
            Some(i) => Dynamic::from(i),
            None => Dynamic::from(n.as_f64().unwrap_or(f64::NAN)),
        },
        Value::String(s) => Dynamic::from(s.clone()),
        Value::Array(a) => Dynamic::from_array(a.iter().map(to_script).collect()),
        Value::Object(o) => {
            let mut m = Map::new();
            for (k, x) in o {
                m.insert(k.as_str().into(), to_script(x));
            }
            Dynamic::from_map(m)
        }
    }
}

/// A script value to JSON; a value JSON cannot hold (a function, a custom type) is an error.
pub fn from_script(d: &Dynamic) -> Result<Value, String> {
    if d.is_unit() {
        return Ok(Value::Null);
    }
    if let Ok(b) = d.as_bool() {
        return Ok(Value::Bool(b));
    }
    if let Ok(i) = d.as_int() {
        return Ok(Value::from(i));
    }
    if let Ok(f) = d.as_float() {
        return serde_json::Number::from_f64(f).map(Value::Number).ok_or_else(|| format!("{f} is not a JSON number"));
    }
    if d.is_string() {
        return Ok(Value::String(d.clone().into_string().unwrap_or_default()));
    }
    if d.is_array() {
        let a = d.clone().into_array().unwrap_or_default();
        return a.iter().map(from_script).collect::<Result<Vec<_>, _>>().map(Value::Array);
    }
    if d.is_map() {
        let m = d.clone().cast::<Map>();
        let mut o = serde_json::Map::new();
        for (k, x) in &m {
            o.insert(k.to_string(), from_script(x)?);
        }
        return Ok(Value::Object(o));
    }
    Err(format!("a value of type {} cannot be sent to the engine", d.type_name()))
}
