//! Running and checking a script (S4). `run` starts one worker process through monty-pool, feeds it the script and answers each pause (an engine
//! command, an unknown name, a file or clock request) until the script completes or a limit ends it. The engine never leaves this process.
//!
//! What a script can do, and what it gets when it tries more:
//!
//! * engine commands: `nyx.<command>(...)`, answered by `Engine::call`; an engine error is raised in the script as `RuntimeError`, which it can catch
//! * files: none. `open`, `pathlib` and the `os` calls raise `PermissionError`; the commands that read or write files (`project_save`, ...) work
//!   only inside `Options::root`, enforced by the engine's path guard
//! * clock and random numbers: the clock is refused (`RuntimeError`), `random` starts from the fixed seed 0, `sleep` returns at once, so the same
//!   script and inputs give the same numbers (ADR 11)
//! * time, memory, depth, commands, output: `Limits`; the worker process is stopped past them and the engine goes on

use std::sync::{Arc, Mutex};
use std::time::Duration;

use monty_pool::{on_print_sync, Pool, PoolConfig, PoolError, ReplConfig, ResumeValue, TurnEvent};
use monty_types::{
    DateTimeSource, ExcType, MontyException, MontyObject, MontyUuid, NameLookupResult, NamedValues, OsPolicy, ProcessTime, RandomSeed, RandomStart,
    ResourceLimits, SleepMode,
};
use nyx_engine::Engine;
use serde_json::Value;

use crate::convert::{args_json, from_object, to_object};
use crate::{find_worker, stubs, Limits, Options, Outcome};

/// The clock is refused, `random` is seeded, sleeping returns at once: a run gives the same numbers every time.
fn os_policy() -> OsPolicy {
    OsPolicy {
        datetime: DateTimeSource::CallHost,
        sleep: SleepMode::Zero,
        process_time: ProcessTime::Zero,
        random_start: RandomStart::Seed(RandomSeed::Int(0.into())),
        ..OsPolicy::default()
    }
}

fn repl_config(name: &str, engine: &Engine, limits: &Limits, type_check: bool) -> ReplConfig {
    let rl = ResourceLimits::default()
        .max_feed_duration(Duration::from_secs_f64(limits.seconds))
        .max_memory(limits.memory_bytes)
        .max_recursion_depth(limits.depth)
        .max_suspensions(limits.calls);
    ReplConfig {
        script_name: name.to_string(),
        limits: Some(rl),
        type_check,
        type_check_stubs: Some(stubs::nyx_pyi(engine)),
        print_flush_interval: Some(Duration::ZERO),
        os_policy: os_policy(),
        ..ReplConfig::default()
    }
}

/// The `nyx` module a script sees: one function per command of the catalog and `call`.
fn nyx_module(engine: &Engine) -> NamedValues {
    let mut names: Vec<String> = engine.catalog().into_iter().map(|c| c.name).collect();
    names.push("call".to_string());
    let module = MontyObject::module(
        "nyx",
        MontyUuid::from_u128(0x6e79_785f_7363_7269_7074),
        names.iter().map(|n| (MontyObject::string(n.as_str()), MontyObject::function(n.as_str(), None))),
    );
    let mut inputs = NamedValues::new();
    inputs.push("nyx", module);
    inputs
}

fn dispatch(engine: &mut Engine, name: &str, args: &monty_types::CallArgs) -> Result<Value, String> {
    if name == "call" {
        let cmd = args.arg(0).and_then(|a| a.as_str().map(str::to_string)).ok_or("nyx.call needs the command name first")?;
        let a = args_json(args, 1)?;
        return engine.call(&cmd, &a).map_err(|e| e.to_string());
    }
    let a = args_json(args, 0)?;
    engine.call(name, &a).map_err(|e| e.to_string())
}

/// What the pool's error means to the person who wrote the script.
fn explain(e: PoolError, limits: &Limits) -> String {
    match e {
        PoolError::Runtime(exc) => exc.to_string(),
        PoolError::Typing(t) => format!("the script failed the static check, nothing ran:\n{t}"),
        PoolError::Timeout { timeout } => format!("the script computed for more than {timeout:?} and was stopped (limit `seconds` = {} s of script time; engine time is not counted)", limits.seconds),
        PoolError::Crashed { cause, .. } => format!(
            "the script worker stopped ({cause:?}). Most likely the script used more than {} MiB of memory (limit `memory_bytes`) or overflowed the stack; the engine was not affected",
            limits.memory_bytes >> 20
        ),
        other => other.to_string(),
    }
}

/// An exception the script can catch for a name the host refuses.
fn refuse(kind: ExcType, text: impl Into<String>) -> ResumeValue {
    ResumeValue::Error(MontyException::new(kind, Some(text.into())))
}

async fn drive(engine: &mut Engine, name: &str, src: &str, limits: &Limits, worker: std::path::PathBuf, execute: bool) -> Result<Outcome, String> {
    let mut cfg = PoolConfig::subprocess(worker);
    cfg.min_processes = 0;
    cfg.max_processes = 1;
    cfg.max_checkouts_per_worker = Some(1);
    // a hang the sandbox's own limits cannot see (the time between two commands, counted without the engine's answers)
    cfg.request_timeout = Some(Duration::from_secs_f64(limits.seconds + 10.0));
    let pool = Pool::new(cfg).await.map_err(|e| explain(e, limits))?;
    let repl = repl_config(name, engine, limits, true);
    let mut session = pool.checkout(&repl).await.map_err(|e| explain(e, limits))?;

    let out = Arc::new(Mutex::new(String::new()));
    let cap = limits.output_bytes;
    let sink = out.clone();
    let mut on_print = on_print_sync(move |_stream, text| {
        let mut o = sink.lock().unwrap();
        if o.len() < cap {
            o.push_str(text);
        }
    });

    // `execute` false: the type check is the point and the script must not run. The check happens before the first instruction; a script that
    // passes it is stopped at its first engine command, which `check` never answers.
    let mut event = session.feed(src, nyx_module(engine), vec![], false, &mut on_print).await.map_err(|e| explain(e, limits))?;
    let mut calls = 0u64;
    let result = loop {
        match event {
            TurnEvent::Complete(v) => break from_object(v.as_ref()),
            TurnEvent::FunctionCall { function_name, args, .. } => {
                if !execute {
                    break Ok(Value::Null);
                }
                calls += 1;
                let reply = match dispatch(engine, &function_name, &args) {
                    Ok(v) => ResumeValue::Return(to_object(&v)),
                    Err(e) => refuse(ExcType::RuntimeError, format!("{function_name}: {e}")),
                };
                event = session.resume(reply, &mut on_print).await.map_err(|e| explain(e, limits))?;
            }
            // a name nothing defines: Python's own NameError with the line
            TurnEvent::NameLookup { .. } => {
                event = session.resume_name_lookup(NameLookupResult::Undefined, &mut on_print).await.map_err(|e| explain(e, limits))?;
            }
            // files, imports of modules that are not built in, the clock: refused. The sandbox raises PermissionError (with the path) or
            // ModuleNotFoundError itself; the clock gets a sentence that says why.
            TurnEvent::OsCall { function_name, .. } => {
                let clock = function_name.contains("time") || function_name.contains("date") || function_name.contains("now") || function_name.contains("today");
                let reply = if clock {
                    refuse(ExcType::RuntimeError, "not available in nyx scripts: no clock (a run must give the same numbers every time)")
                } else {
                    ResumeValue::NotHandled
                };
                event = session.resume(reply, &mut on_print).await.map_err(|e| explain(e, limits))?;
            }
            TurnEvent::ResolveFutures { .. } => {
                return Err("the script awaited something (async code is not supported in nyx scripts)".to_string());
            }
        }
    };
    let _ = session.finish().await;
    pool.close().await;
    let printed: Vec<String> = out.lock().unwrap().lines().map(str::to_string).collect();
    Ok(Outcome { value: result?, printed, calls })
}

fn with_runtime<T>(f: impl std::future::Future<Output = T>) -> Result<T, String> {
    let rt = tokio::runtime::Builder::new_current_thread().enable_all().build().map_err(|e| format!("no async runtime: {e}"))?;
    Ok(rt.block_on(f))
}

fn worker_of(opts: &Options) -> Result<std::path::PathBuf, String> {
    opts.worker.clone().or_else(find_worker).ok_or_else(|| "the script worker program `nyx-script-worker` was not found next to this program (build it with `cargo build --release -p nyx-script`, or name it in NYX_SCRIPT_WORKER)".to_string())
}

/// Runs a script. Its value is its last expression. The engine's path guard is on for the run and restored after it.
pub fn run(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<Outcome, String> {
    let worker = worker_of(opts)?;
    let before = engine.path_root().map(std::path::Path::to_path_buf);
    engine.set_path_root(Some(&opts.root))?;
    let r = with_runtime(drive(engine, name, src, limits, worker, true));
    let _ = engine.set_path_root(before.as_deref());
    r?
}

/// Parses and type-checks a script against the commands of `engine` and runs nothing: a misspelled command or argument, a missing argument, a
/// syntax error. `Ok` when the script is fine.
pub fn check(engine: &mut Engine, name: &str, src: &str, limits: &Limits, opts: &Options) -> Result<(), String> {
    let worker = worker_of(opts)?;
    with_runtime(drive(engine, name, src, limits, worker, false))?.map(|_| ())
}
