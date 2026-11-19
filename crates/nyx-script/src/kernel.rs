//! The kernel (productivity plan S7, ADR [94]): one Monty session in a worker process that outlives a script. `run` feeds a script to it and
//! answers each pause (an engine command, a name nothing defines, a clock or file request) until the script completes, fails, or the time the
//! caller gave is spent; `wait` continues a script that was left waiting; `check` type-checks a script against the commands and the variables of
//! the session without running it. `host::run` and `host::check` (a person's `nyx run` and `nyx check`) are one script on a kernel of their own.
//!
//! * variables, functions and imports of one script are there in the next (Monty feeds are a REPL); a worker that dies (memory, stack) is
//!   replaced from the dump taken after the last script that finished, so they are back
//! * a command that waits for a job (`solve_run`, `convergence_check`, `job_wait`, `queue_wait`) is served as: start the job without waiting,
//!   then poll it; when the caller's time is spent the script stays paused at the call and the answer is `Step::Running`. The polling is the
//!   host's, not the script's: it is not a command the script sent, so it does not count against `Limits::calls`
//! * `nyx.aio.<command>` is the same four commands as futures (`await asyncio.gather(...)`); a future left pending when a script ends has its
//!   job cancelled
//! * the project-discarding commands (`is_destructive`) are refused unless the run says they are allowed

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use monty_pool::{on_print_sync, Checkout, Pool, PoolConfig, PoolError, ReplConfig, ResumeValue, TurnEvent};
use monty_types::{
    CallArgs, DateTimeSource, ExcType, MontyException, MontyObject, MontyUuid, NameLookupResult, NamedValues, OsPolicy, ProcessTime, RandomSeed, RandomStart,
    ResourceLimits, SleepMode,
};
use nyx_engine::Engine;
use serde_json::{json, Value};

use crate::convert::{args_json, from_object, to_object};
use crate::{find_worker, stubs, trace, Limits, Options, Outcome};

/// The commands that wait for a job. Each is served by starting the job and polling it, and each has an `nyx.aio` twin.
pub const POLLABLE: &[&str] = &["solve_run", "convergence_check", "job_wait", "queue_wait"];

/// How often the host asks the engine whether a job has finished.
const POLL: Duration = Duration::from_millis(20);

/// Commands that discard the project or a part of it.
pub fn is_destructive(name: &str) -> bool {
    matches!(name, "project_new" | "project_load" | "project_open" | "store_open") || name.starts_with("remove_")
}

/// What a script that is not finished is waiting for.
#[derive(Debug)]
pub struct Progress {
    /// what the script printed since the last answer
    pub printed: Vec<String>,
    pub calls: u64,
    /// the command and job it waits for
    pub waiting: String,
}

/// The answer to `run` and `wait`.
#[derive(Debug)]
pub enum Step {
    Done(Outcome),
    Running(Progress),
}

/// A job the host polls for a script.
struct Poll {
    cmd: String,
    job: u64,
    started: Instant,
    /// the call's own timeout (`job_wait`, `queue_wait`): at it the call returns the job's entry with `waiting: true`
    timeout: Option<Duration>,
}

enum Fut {
    Ready(Result<Value, String>),
    Poll(Poll),
}

enum Paused {
    Call(Poll),
    Futures(Vec<u32>),
}

enum Next {
    Event(TurnEvent),
    PollCall(Poll),
    Settle(Vec<u32>),
    Resume,
}

enum Fail {
    /// the script failed; the session is fine
    Soft(String),
    /// the worker is gone or out of step: the session is replaced
    Fatal(String),
}

enum Act {
    Value(Result<Value, String>),
    Poll(Poll),
    Future,
}

enum Polled {
    Pending,
    Done(Result<Value, String>),
}

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

fn repl_config(engine: &Engine, limits: &Limits) -> ReplConfig {
    let rl = ResourceLimits::default()
        .max_feed_duration(Duration::from_secs_f64(limits.seconds))
        .max_memory(limits.memory_bytes)
        .max_recursion_depth(limits.depth)
        .max_suspensions(limits.calls);
    ReplConfig {
        script_name: "script".to_string(),
        limits: Some(rl),
        type_check: true,
        type_check_stubs: Some(stubs::nyx_pyi(engine)),
        print_flush_interval: Some(Duration::ZERO),
        os_policy: os_policy(),
        ..ReplConfig::default()
    }
}

/// The `nyx` module a script sees: one function per command of the catalog, `call`, and `aio` with the pollable commands as async functions.
fn nyx_module(engine: &Engine) -> NamedValues {
    let mut names: Vec<String> = engine.catalog().into_iter().map(|c| c.name).collect();
    names.push("call".to_string());
    let mut attrs: Vec<(MontyObject, MontyObject)> = names.iter().map(|n| (MontyObject::string(n.as_str()), MontyObject::function(n.as_str(), None))).collect();
    let aio = MontyObject::module(
        "nyx.aio",
        MontyUuid::from_u128(0x006e_7978_5f61_696f),
        POLLABLE.iter().map(|n| (MontyObject::string(*n), MontyObject::async_function(format!("aio_{n}"), None))),
    );
    attrs.push((MontyObject::string("aio"), aio));
    let module = MontyObject::module("nyx", MontyUuid::from_u128(0x6e79_785f_7363_7269_7074), attrs);
    let mut inputs = NamedValues::new();
    inputs.push("nyx", module);
    inputs
}

fn params_input(params: &Value) -> NamedValues {
    let mut i = NamedValues::new();
    i.push("params", to_object(params));
    i
}

/// What the pool's error means to the person who wrote the script.
pub(crate) fn explain(e: PoolError, limits: &Limits) -> String {
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

fn reply(name: &str, r: Result<Value, String>) -> ResumeValue {
    match r {
        Ok(v) => ResumeValue::Return(to_object(&v)),
        Err(e) => refuse(ExcType::RuntimeError, format!("{name}: {e}")),
    }
}

/// The command a call names, and its arguments as the engine's JSON (`nyx.call("cmd", {...})` names it first).
fn command_of(name: &str, args: &CallArgs) -> Result<(String, Value), String> {
    if name == "call" {
        let cmd = args.arg(0).and_then(|a| a.as_str().map(str::to_string)).ok_or("nyx.call needs the command name first")?;
        Ok((cmd, args_json(args, 1)?))
    } else {
        Ok((name.to_string(), args_json(args, 0)?))
    }
}

/// Starts the job of a command that waits for one, or says it is a plain command (`None`).
fn start_poll(engine: &mut Engine, cmd: &str, a: &Value) -> Result<Option<Poll>, String> {
    let mk = |job: u64, timeout: Option<Duration>| Poll { cmd: cmd.to_string(), job, started: Instant::now(), timeout };
    match cmd {
        "solve_run" | "convergence_check" => {
            if a.get("wait").and_then(Value::as_bool) == Some(false) {
                return Ok(None);
            }
            let mut b = a.clone();
            b["wait"] = json!(false);
            let r = engine.call(cmd, &b).map_err(|e| e.to_string())?;
            let job = r.get("job").and_then(Value::as_u64).ok_or_else(|| format!("{cmd} did not start a job: {r}"))?;
            Ok(Some(mk(job, None)))
        }
        "job_wait" | "queue_wait" => {
            let Some(job) = a.get("id").and_then(Value::as_u64) else { return Ok(None) }; // the engine says what is wrong with the arguments
            let timeout = match a.get("timeout_s").and_then(Value::as_f64) {
                Some(t) if t >= 0.0 && t.is_finite() => Some(Duration::from_secs_f64(t)),
                Some(_) => return Ok(None),
                None if cmd == "queue_wait" => return Ok(None),
                None => None,
            };
            Ok(Some(mk(job, timeout)))
        }
        _ => Ok(None),
    }
}

/// Asks the engine once whether the job has finished (applying its outcome when it has).
fn poll_once(engine: &mut Engine, p: &Poll) -> Polled {
    match engine.call("job_wait", &json!({ "id": p.job, "timeout_s": 0.0 })) {
        Ok(v) if v.get("waiting").and_then(Value::as_bool) == Some(true) => {
            if let Some(t) = p.timeout {
                if p.started.elapsed() >= t {
                    return Polled::Done(Ok(v));
                }
            } else if let Some(why) = engine.job_blocked(p.job) {
                return Polled::Done(Err(format!("job {} cannot start: {why}", p.job)));
            }
            Polled::Pending
        }
        Ok(v) => Polled::Done(Ok(v)),
        Err(e) => Polled::Done(Err(e.to_string())),
    }
}

/// A Monty session that outlives a script.
pub struct Kernel {
    rt: Arc<tokio::runtime::Runtime>,
    pool: Pool,
    repl: ReplConfig,
    limits: Limits,
    opts: Options,
    session: Option<Checkout>,
    /// the session after the prelude and after each script that finished: what a replacement worker starts from
    dump: Option<Vec<u8>>,
    out: Arc<Mutex<String>>,
    paused: Option<Paused>,
    futs: HashMap<u32, Fut>,
    calls: u64,
    allow_destructive: bool,
    /// the engine's path root before the run that is still going, to put back
    root_before: Option<Option<PathBuf>>,
}

macro_rules! printer {
    ($k:expr) => {{
        let sink = $k.out.clone();
        let cap = $k.limits.output_bytes;
        on_print_sync(move |_stream, text| {
            let mut o = sink.lock().unwrap();
            if o.len() < cap {
                o.push_str(text);
            }
        })
    }};
}

impl Kernel {
    /// A kernel for `engine`'s catalog. No process is started until the first script.
    pub fn new(engine: &Engine, limits: Limits, opts: Options) -> Result<Kernel, String> {
        let worker = opts.worker.clone().or_else(find_worker).ok_or_else(|| {
            "the script worker program `nyx-script-worker` was not found next to this program (build it with `cargo build --release -p nyx-script`, or name it in NYX_SCRIPT_WORKER)".to_string()
        })?;
        let rt = tokio::runtime::Builder::new_current_thread().enable_all().build().map_err(|e| format!("no async runtime: {e}"))?;
        let mut cfg = PoolConfig::subprocess(worker);
        cfg.min_processes = 0;
        // the kernel's session and the throw-away one of `check`
        cfg.max_processes = 2;
        cfg.max_checkouts_per_worker = Some(1);
        // a hang the sandbox's own limits cannot see (the time between two commands, counted without the engine's answers)
        cfg.request_timeout = Some(Duration::from_secs_f64(limits.seconds + 10.0));
        let pool = rt.block_on(Pool::new(cfg)).map_err(|e| explain(e, &limits))?;
        let repl = repl_config(engine, &limits);
        Ok(Kernel {
            rt: Arc::new(rt),
            pool,
            repl,
            limits,
            opts,
            session: None,
            dump: None,
            out: Arc::new(Mutex::new(String::new())),
            paused: None,
            futs: HashMap::new(),
            calls: 0,
            allow_destructive: true,
            root_before: None,
        })
    }

    /// The name the static check gives the script in its messages (`file.nys:2:1`). It is the name of the session, so it is set before the first
    /// script: `nyx run` and `nyx check` name their file; the MCP server's kernel, which serves many scripts, keeps a neutral one.
    pub fn set_script_name(&mut self, name: &str) {
        self.repl.script_name = name.to_string();
    }

    /// What the script that has not finished is waiting for.
    pub fn waiting(&self) -> Option<String> {
        self.paused.as_ref().map(|p| match p {
            Paused::Call(c) => format!("{} (job {})", c.cmd, c.job),
            Paused::Futures(ids) => format!("{} future(s) started with nyx.aio", ids.len()),
        })
    }

    fn take_printed(&self) -> Vec<String> {
        let mut o = self.out.lock().unwrap();
        let lines = o.lines().map(str::to_string).collect();
        o.clear();
        lines
    }

    /// Feeds `src` to the session. It returns when the script is finished, or when `wait` is spent and it is waiting for a job.
    pub fn run(&mut self, engine: &mut Engine, name: &str, src: &str, params: Option<&Value>, wait: Option<Duration>, allow_destructive: bool) -> Result<Step, String> {
        // a branch on a value of the point in a function given as a formula (S6): refused before anything runs, with the line. It reads this
        // script only: a function defined by an earlier script and given as a formula here is not seen
        trace::guard(name, src)?;
        if let Some(w) = self.waiting() {
            return Err(format!("a script is still running, waiting for {w}: continue it with script_wait, or abandon it with script_wait cancel"));
        }
        let rt = self.rt.clone();
        self.calls = 0;
        self.allow_destructive = allow_destructive;
        self.out.lock().unwrap().clear();
        self.root_before = Some(engine.path_root().map(std::path::Path::to_path_buf));
        engine.set_path_root(Some(&self.opts.root))?;
        // the runs the script makes keep its text and hash (S4)
        engine.set_script(Some((name.to_string(), src.to_string())));
        let deadline = wait.map(|w| Instant::now() + w);
        let params = params.cloned().unwrap_or_else(|| json!({}));
        let r = rt.block_on(self.start(engine, name, src, &params, deadline));
        self.after(engine, &r);
        r
    }

    /// Continues the script that is waiting; at the end of `wait` it may still be.
    pub fn wait(&mut self, engine: &mut Engine, wait: Option<Duration>) -> Result<Step, String> {
        if self.paused.is_none() {
            return Err("no script is running".to_string());
        }
        let rt = self.rt.clone();
        let deadline = wait.map(|w| Instant::now() + w);
        let Some(mut s) = self.session.take() else { return Err("the script's session is gone".to_string()) };
        let r = rt.block_on(async {
            let r = self.advance(&mut s, engine, Next::Resume, deadline).await;
            self.settle_session(s, r).await
        });
        self.after(engine, &r);
        r
    }

    /// Abandons the script that is waiting: its job is cancelled and the session starts again from the last script that finished.
    pub fn cancel(&mut self, engine: &mut Engine) -> Result<(), String> {
        let Some(p) = self.paused.take() else { return Err("no script is running".to_string()) };
        if let Paused::Call(c) = &p {
            let _ = engine.call("job_cancel", &json!({ "id": c.job }));
        }
        self.cancel_futures(engine);
        if let Some(s) = self.session.take() {
            let rt = self.rt.clone();
            let _ = rt.block_on(s.terminate());
        }
        self.restore_engine(engine);
        self.out.lock().unwrap().clear();
        Ok(())
    }

    /// Parses and type-checks `src` against the commands and the variables of the session and runs nothing that matters: the session it is
    /// checked in is a copy of the kernel's, thrown away after.
    pub fn check(&mut self, engine: &mut Engine, name: &str, src: &str) -> Result<(), String> {
        trace::guard(name, src)?;
        let rt = self.rt.clone();
        rt.block_on(async {
            let mut s = self.open_session(engine).await?;
            let mut on_print = printer!(self);
            // the check happens before the first instruction; a script that passes it is stopped at its first engine command, which `check`
            // never answers
            let r = s.feed(src, params_input(&json!({})), vec![], false, &mut on_print).await;
            let _ = s.finish().await;
            self.out.lock().unwrap().clear();
            match r {
                Ok(_) => Ok(()),
                // the part of the script before its first command ran, with no `params` (it may read one): what that raises is not a fault of
                // the script, only what the static check and the parser find is
                Err(PoolError::Runtime(exc)) if exc.exc_type() != ExcType::SyntaxError => Ok(()),
                Err(e) => Err(explain(e, &self.limits)),
            }
        })
    }

    // ---- the session ----

    /// A new session: restored from the last dump, or fresh with the prelude (the symbolic values and the `nyx` object that traces a function
    /// given as a formula, S6). The prelude makes no engine command, is not type-checked (the stub declares its names) and is a feed of its own,
    /// so the script's tracebacks keep the script's line numbers.
    async fn open_session(&mut self, engine: &Engine) -> Result<Checkout, String> {
        let mut c = self.pool.checkout(&self.repl).await.map_err(|e| explain(e, &self.limits))?;
        let mut on_print = printer!(self);
        if let Some(d) = &self.dump {
            c.restore(d.clone(), vec![], &mut on_print).await.map_err(|e| explain(e, &self.limits))?;
            return Ok(c);
        }
        match c.feed(trace::prelude(engine), nyx_module(engine), vec![], true, &mut on_print).await.map_err(|e| explain(e, &self.limits))? {
            TurnEvent::Complete(_) => {}
            _ => return Err("internal error: the script prelude paused".to_string()),
        }
        self.dump = c.dump().await.ok();
        Ok(c)
    }

    async fn start(&mut self, engine: &mut Engine, _name: &str, src: &str, params: &Value, deadline: Option<Instant>) -> Result<Step, String> {
        let mut s = match self.session.take() {
            Some(s) => s,
            None => self.open_session(engine).await?,
        };
        let mut on_print = printer!(self);
        let first = s.feed(src, params_input(params), vec![], false, &mut on_print).await;
        let r = match first {
            Ok(ev) => self.advance(&mut s, engine, Next::Event(ev), deadline).await,
            Err(e) => Err(self.pool_fail(e)),
        };
        self.settle_session(s, r).await
    }

    /// Puts the session back (a script that finished also updates the dump), or drops it when the worker is gone.
    async fn settle_session(&mut self, mut s: Checkout, r: Result<Step, Fail>) -> Result<Step, String> {
        match r {
            Ok(step) => {
                if matches!(step, Step::Done(_)) {
                    if let Ok(d) = s.dump().await {
                        self.dump = Some(d);
                    }
                }
                self.session = Some(s);
                Ok(step)
            }
            Err(Fail::Soft(m)) => {
                self.session = Some(s);
                Err(m)
            }
            Err(Fail::Fatal(m)) => {
                let _ = s.terminate().await;
                self.paused = None;
                let back = if self.dump.is_some() { "; the variables of the scripts that finished are back in a new worker" } else { "" };
                Err(format!("{m}{back}"))
            }
        }
    }

    /// Once a script has finished, failed or been abandoned: what it started is cancelled and the engine's globals are put back.
    fn after(&mut self, engine: &mut Engine, r: &Result<Step, String>) {
        if !matches!(r, Ok(Step::Running(_))) {
            self.paused = None;
            self.cancel_futures(engine);
            self.restore_engine(engine);
        }
    }

    fn restore_engine(&mut self, engine: &mut Engine) {
        engine.set_script(None);
        if let Some(before) = self.root_before.take() {
            let _ = engine.set_path_root(before.as_deref());
        }
    }

    fn cancel_futures(&mut self, engine: &mut Engine) {
        for (_, f) in self.futs.drain() {
            if let Fut::Poll(p) = f {
                let _ = engine.call("job_cancel", &json!({ "id": p.job }));
            }
        }
    }

    fn pool_fail(&self, e: PoolError) -> Fail {
        let soft = matches!(e, PoolError::Runtime(_) | PoolError::Typing(_));
        let m = explain(e, &self.limits);
        if soft {
            Fail::Soft(m)
        } else {
            Fail::Fatal(m)
        }
    }

    // ---- the driver ----

    /// Answers pauses until the script completes, fails, or has to wait past `deadline`.
    async fn advance(&mut self, s: &mut Checkout, engine: &mut Engine, start: Next, deadline: Option<Instant>) -> Result<Step, Fail> {
        let mut on_print = printer!(self);
        let mut next = start;
        loop {
            next = match next {
                Next::Event(TurnEvent::Complete(v)) => {
                    let value = from_object(v.as_ref()).map_err(Fail::Soft)?;
                    return Ok(Step::Done(Outcome { value, printed: self.take_printed(), calls: self.calls }));
                }
                Next::Event(TurnEvent::FunctionCall { function_name, args, call_id, .. }) => {
                    self.calls += 1;
                    match self.act(engine, &function_name, call_id, &args) {
                        Act::Value(r) => Next::Event(s.resume(reply(&function_name, r), &mut on_print).await.map_err(|e| self.pool_fail(e))?),
                        Act::Future => Next::Event(s.resume(ResumeValue::Future, &mut on_print).await.map_err(|e| self.pool_fail(e))?),
                        Act::Poll(p) => Next::PollCall(p),
                    }
                }
                Next::Event(TurnEvent::ResolveFutures { pending_call_ids, .. }) => Next::Settle(pending_call_ids),
                // a name nothing defines: Python's own NameError with the line
                Next::Event(TurnEvent::NameLookup { .. }) => {
                    Next::Event(s.resume_name_lookup(NameLookupResult::Undefined, &mut on_print).await.map_err(|e| self.pool_fail(e))?)
                }
                // files, imports of modules that are not built in, the clock: refused. The sandbox raises PermissionError (with the path) or
                // ModuleNotFoundError itself; the clock gets a sentence that says why.
                Next::Event(TurnEvent::OsCall { function_name, .. }) => {
                    let clock = function_name.contains("time") || function_name.contains("date") || function_name.contains("now") || function_name.contains("today");
                    let r = if clock {
                        refuse(ExcType::RuntimeError, "not available in nyx scripts: no clock (a run must give the same numbers every time)")
                    } else {
                        ResumeValue::NotHandled
                    };
                    Next::Event(s.resume(r, &mut on_print).await.map_err(|e| self.pool_fail(e))?)
                }
                Next::PollCall(p) => match self.poll_until(engine, &p, deadline).await {
                    Some(r) => Next::Event(s.resume(reply(&p.cmd, r), &mut on_print).await.map_err(|e| self.pool_fail(e))?),
                    None => {
                        self.paused = Some(Paused::Call(p));
                        return Ok(Step::Running(self.progress()));
                    }
                },
                Next::Settle(ids) => match self.settle(engine, &ids, deadline).await {
                    Some(results) => Next::Event(s.resume_futures(results, &mut on_print).await.map_err(|e| self.pool_fail(e))?),
                    None => {
                        self.paused = Some(Paused::Futures(ids));
                        return Ok(Step::Running(self.progress()));
                    }
                },
                Next::Resume => match self.paused.take() {
                    Some(Paused::Call(p)) => Next::PollCall(p),
                    Some(Paused::Futures(ids)) => Next::Settle(ids),
                    None => return Err(Fail::Soft("no script is running".to_string())),
                },
            };
        }
    }

    fn progress(&self) -> Progress {
        Progress { printed: self.take_printed(), calls: self.calls, waiting: self.waiting().unwrap_or_default() }
    }

    /// One command of a script.
    fn act(&mut self, engine: &mut Engine, name: &str, call_id: u32, args: &CallArgs) -> Act {
        let (is_aio, base) = match name.strip_prefix("aio_") {
            Some(b) => (true, b),
            None => (false, name),
        };
        let (cmd, a) = match command_of(base, args) {
            Ok(x) => x,
            Err(e) => return self.value_or_future(is_aio, call_id, Err(e)),
        };
        if !self.allow_destructive && is_destructive(&cmd) {
            let e = format!("{cmd} discards the project or a part of it: it is refused unless script_run is called with allow_destructive true");
            return self.value_or_future(is_aio, call_id, Err(e));
        }
        match start_poll(engine, &cmd, &a) {
            Ok(Some(p)) if is_aio => {
                self.futs.insert(call_id, Fut::Poll(p));
                Act::Future
            }
            Ok(Some(p)) => Act::Poll(p),
            Ok(None) => {
                let r = engine.call(&cmd, &a).map_err(|e| e.to_string());
                self.value_or_future(is_aio, call_id, r)
            }
            Err(e) => self.value_or_future(is_aio, call_id, Err(e)),
        }
    }

    fn value_or_future(&mut self, is_aio: bool, call_id: u32, r: Result<Value, String>) -> Act {
        if is_aio {
            self.futs.insert(call_id, Fut::Ready(r));
            Act::Future
        } else {
            Act::Value(r)
        }
    }

    /// Polls one job until it has finished, its own timeout passed, or `deadline` (`None` when it was not finished by then).
    async fn poll_until(&mut self, engine: &mut Engine, p: &Poll, deadline: Option<Instant>) -> Option<Result<Value, String>> {
        loop {
            if let Polled::Done(r) = poll_once(engine, p) {
                return Some(r);
            }
            if deadline.is_some_and(|d| Instant::now() >= d) {
                return None;
            }
            tokio::time::sleep(POLL).await;
        }
    }

    /// Answers the futures that can be answered; waits until at least one can, or `deadline` (`None`).
    async fn settle(&mut self, engine: &mut Engine, ids: &[u32], deadline: Option<Instant>) -> Option<Vec<(u32, ResumeValue)>> {
        loop {
            let mut out: Vec<(u32, ResumeValue)> = vec![];
            for id in ids {
                let done = match self.futs.get(id) {
                    Some(Fut::Ready(_)) => match self.futs.remove(id) {
                        Some(Fut::Ready(r)) => Some(("nyx.aio".to_string(), r)),
                        _ => None,
                    },
                    Some(Fut::Poll(p)) => match poll_once(engine, p) {
                        Polled::Done(r) => {
                            let cmd = p.cmd.clone();
                            self.futs.remove(id);
                            Some((cmd, r))
                        }
                        Polled::Pending => None,
                    },
                    None => Some(("nyx.aio".to_string(), Err("this future is not known to the host".to_string()))),
                };
                if let Some((cmd, r)) = done {
                    out.push((*id, reply(&cmd, r)));
                }
            }
            if !out.is_empty() {
                return Some(out);
            }
            if deadline.is_some_and(|d| Instant::now() >= d) {
                return None;
            }
            tokio::time::sleep(POLL).await;
        }
    }
}

impl Drop for Kernel {
    fn drop(&mut self) {
        let rt = self.rt.clone();
        if let Some(s) = self.session.take() {
            let _ = rt.block_on(s.finish());
        }
        rt.block_on(self.pool.close());
    }
}
