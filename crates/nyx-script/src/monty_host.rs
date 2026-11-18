//! The Monty binding (S4 spike): a script written in Python syntax calls `nyx.<command>(name=..., ...)` for every command of the catalog, and
//! `nyx.call(name, **args)`. Monty runs the script on its own VM and pauses it at each call; this host answers with `Engine::call`.
//!
//! An engine error is raised in the script as a `RuntimeError`, so `try`/`except` recovers. The sandbox has no file, network, clock or
//! random access except through functions the host registers, and the host registers none of those.

use std::time::Duration;

use monty::{MontyRun, RunProgress};
use monty_types::{
    CallArgs, CompileOptions, DateTimeSource, ExcType, ExtFunctionResult, MontyException, MontyObject, MontyUuid, NameLookupResult,
    ObjectRef, OsFunctionCall, OsPolicy, PrintWriter, ProcessTime, RandomSeed, RandomStart, ResourceLimits, ResourceTracker, SleepMode,
};
use nyx_engine::Engine;
use serde_json::Value;

use crate::{Limits, Outcome};

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

/// A script value to JSON; a value JSON cannot hold (a function, a set, a non-text dict key) is an error.
pub fn from_object(o: ObjectRef) -> Result<Value, String> {
    match o.type_name() {
        "NoneType" => Ok(Value::Null),
        "bool" => Ok(Value::Bool(o.as_bool().unwrap_or(false))),
        "int" => match o.as_int() {
            Some(i) => Ok(Value::from(i)),
            None => Err("an integer larger than 64 bits cannot be sent to the engine (use a float)".to_string()),
        },
        "float" => {
            let f = o.as_float().unwrap_or(f64::NAN);
            serde_json::Number::from_f64(f).map(Value::Number).ok_or_else(|| format!("{f} is not a JSON number"))
        }
        "str" => Ok(Value::String(o.as_str().unwrap_or_default().to_string())),
        "list" | "tuple" => o.items().unwrap_or_default().into_iter().map(from_object).collect::<Result<Vec<_>, _>>().map(Value::Array),
        "dict" => {
            let mut m = serde_json::Map::new();
            for (k, x) in o.pairs().unwrap_or_default() {
                let key = k.as_str().ok_or_else(|| format!("a dict key of type {} cannot be sent to the engine", k.type_name()))?;
                m.insert(key.to_string(), from_object(x)?);
            }
            Ok(Value::Object(m))
        }
        t => Err(format!("a value of type {t} cannot be sent to the engine")),
    }
}

/// The arguments of one call as the engine's JSON object: keyword arguments, or one positional dict.
fn args_json(args: &CallArgs, skip: usize) -> Result<Value, String> {
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

fn dispatch(engine: &mut Engine, name: &str, args: &CallArgs, calls: &mut u64) -> Result<Value, String> {
    *calls += 1;
    if name == "call" {
        let cmd = args.arg(0).and_then(|a| a.as_str().map(str::to_string)).ok_or("nyx.call needs the command name first")?;
        let a = args_json(args, 1)?;
        return engine.call(&cmd, &a).map_err(|e| e.to_string());
    }
    let a = args_json(args, 0)?;
    engine.call(name, &a).map_err(|e| e.to_string())
}

fn limits_of(l: &Limits) -> ResourceLimits {
    ResourceLimits {
        max_feed_duration: Some(Duration::from_secs_f64(l.seconds)),
        max_recursion_depth: l.depth,
        max_suspensions: l.calls,
        ..ResourceLimits::default()
    }
}

/// Parses and compiles the script and nothing else.
pub fn parse(name: &str, src: &str) -> Result<(), String> {
    MontyRun::new(src.to_owned(), name, vec!["nyx".to_string()], CompileOptions::default()).map(|_| ()).map_err(|e| e.to_string())
}

/// Runs a script. Its value is the script's last expression.
pub fn run(engine: &mut Engine, name: &str, src: &str, limits: &Limits) -> Result<Outcome, String> {
    let mut names: Vec<String> = engine.catalog().into_iter().map(|c| c.name).collect();
    names.push("call".to_string());
    let nyx = MontyObject::module(
        "nyx",
        MontyUuid::from_u128(0x6e79_785f_7363_7269_7074),
        names.iter().map(|n| (MontyObject::string(n.as_str()), MontyObject::function(n.as_str(), None))),
    );
    let known: std::collections::HashSet<&str> = names.iter().map(String::as_str).collect();
    // a script sees no real clock and no OS entropy: a run must give the same numbers every time. `time.time()` and the like ask the host,
    // which refuses them below; `random` starts from a fixed seed (a script may call random.seed(n) itself); sleeping returns at once.
    let policy = OsPolicy {
        datetime: DateTimeSource::CallHost,
        sleep: SleepMode::Zero,
        process_time: ProcessTime::Zero,
        random_start: RandomStart::Seed(RandomSeed::Int(0.into())),
        ..OsPolicy::default()
    };
    let runner = MontyRun::new(src.to_owned(), name, vec!["nyx".to_string()], CompileOptions::default())
        .map_err(|e| e.to_string())?
        .with_os_policy(policy);
    let mut out = String::new();
    let cap = Some(limits.heap_bytes);
    let tracker = ResourceTracker::new(limits_of(limits));
    let mut calls = 0u64;
    let mut progress = runner.start(vec![nyx], tracker, PrintWriter::CollectString(&mut out, cap)).map_err(|e| e.to_string())?;
    let value = loop {
        match progress {
            RunProgress::Complete(v) => break from_object(v.as_ref())?,
            RunProgress::FunctionCall(call) => {
                if calls >= limits.calls as u64 {
                    let e = MontyException::new(ExcType::RuntimeError, Some(format!("the script sent more than {} engine commands", limits.calls)));
                    return Err(call.abort(e, PrintWriter::CollectString(&mut out, cap)).err().map(|e| e.to_string()).unwrap_or_default());
                }
                let fname = call.function_name.clone();
                let result = if !known.contains(fname.as_str()) {
                    ExtFunctionResult::Error(MontyException::new(ExcType::NameError, Some(format!("`{fname}` is not available in nyx scripts"))))
                } else {
                    match dispatch(engine, &fname, &call.args, &mut calls) {
                        Ok(v) => ExtFunctionResult::Return(to_object(&v)),
                        Err(e) => ExtFunctionResult::Error(MontyException::new(ExcType::RuntimeError, Some(format!("{fname}: {e}")))),
                    }
                };
                progress = call.resume(result, PrintWriter::CollectString(&mut out, cap)).map_err(|e| e.to_string())?;
            }
            // a name the script uses but nothing defines: Python's own NameError, with its line
            RunProgress::NameLookup(l) => {
                progress = l.resume(NameLookupResult::Undefined, PrintWriter::CollectString(&mut out, cap)).map_err(|e| e.to_string())?;
            }
            // files, imports of modules that are not built in, the clock: refused, with Python's own exception (PermissionError with the path,
            // ModuleNotFoundError) so the script can catch it and the person sees what was refused
            RunProgress::OsCall(o) => {
                let result = match &o.function_call {
                    OsFunctionCall::Import(_) => o.function_call.no_handler_result(),
                    f if f.fs_primary_path().is_some() => o.function_call.no_handler_result(),
                    _ => ExtFunctionResult::Error(MontyException::new(ExcType::RuntimeError, Some("not available in nyx scripts: no clock, no environment variables (a run must give the same numbers every time)".to_string()))),
                };
                progress = o.resume(result, PrintWriter::CollectString(&mut out, cap)).map_err(|e| e.to_string())?;
            }
            _ => return Err("the script awaited something (async code is not supported in nyx scripts)".to_string()),
        }
    };
    let printed: Vec<String> = out.lines().map(str::to_string).collect();
    Ok(Outcome { value, printed, calls })
}
