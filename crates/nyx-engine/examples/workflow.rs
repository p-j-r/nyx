//! Replays a list of engine commands against an `Engine` and times each one (productivity plan P0). The list is a JSON-lines file, one command per
//! line, the format that `nyx check` and `nyx run` take over later (plan step S1); the workflows of `tools/workflows/` are written in it.
//!
//!   cargo run --release -p nyx-engine --example workflow -- FILE.jsonl [--csv OUT.csv]
//!   (run it from the repository root: the lists name `examples/*.nyx.json` and the material library is looked up from there)
//!
//! A line is `{"cmd": "<engine command>", "args": {...}}` and nothing else is needed: `cmd` is a command of `Engine::catalog()` and `args` its
//! arguments, exactly as an MCP host would send them. Blank lines and lines starting with `#` are skipped. These optional keys control the replay
//! only (they are not engine arguments, and a front-end that runs the list for real ignores them):
//!
//! * `note`: a label shown in the table;
//! * `job_arg`: the name of an argument that is set to the id of the last job a reply named (`solve_run` with `wait: false` returns `{job}`);
//! * `after_stage`: wait until the stage of that job contains this text before issuing the command (the way to cancel "during assembly");
//!   `delay_s` waits that many seconds more inside the stage. The run stops with an error if the job ends first, if a `job_cancel` finds no running
//!   job, or if an `expect_error` command succeeds, so a cancel that did not happen cannot pass for one;
//! * `expect_error`: an error reply is an expected outcome (a cancelled job), not the end of the run.
//!
//! `NYX_NO_ECOQOS=1` (Windows) opts the process out of power throttling, a diagnostic for the slowdown after the first solve.
//!
//! Output: one row per command with its wall time; the time from a `job_cancel` to the end of the `job_wait` that follows it is the time from
//! the cancel request to the job's `Cancelled` state. Wall times of a single run drift by up to 1.8x on this machine between sessions
//! (`NYX2_PENDING_CHANGES.md`, restructure status): compare best of three against a control build.

use nyx_engine::Engine;
use serde_json::{json, Value};
use std::io::Write;
use std::time::{Duration, Instant};

/// Diagnostic, off unless `NYX_NO_ECOQOS` is set: asks Windows not to throttle this process (`PROCESS_POWER_THROTTLING_EXECUTION_SPEED` switched
/// off). Without it a long run in a background process is moved to the efficiency cores after its first seconds and every later solve is 2 to 3
/// times slower (`docs/PRODUCTIVITY_STATUS.md`, finding 1). The engine itself does not do this; the program that hosts it would.
#[cfg(windows)]
fn opt_out_of_power_throttling() {
    #[repr(C)]
    struct State {
        version: u32,
        control_mask: u32,
        state_mask: u32,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetCurrentProcess() -> *mut std::ffi::c_void;
        fn SetProcessInformation(process: *mut std::ffi::c_void, class: i32, info: *const State, size: u32) -> i32;
    }
    // class 4 = ProcessPowerThrottling, version 1, control mask 1 = execution speed, state mask 0 = throttling off
    let st = State { version: 1, control_mask: 1, state_mask: 0 };
    // SAFETY: the pseudo handle of the current process and a correctly sized, initialised struct.
    let ok = unsafe { SetProcessInformation(GetCurrentProcess(), 4, &st, std::mem::size_of::<State>() as u32) };
    eprintln!("power throttling opt-out: {}", if ok != 0 { "set" } else { "refused" });
}

fn main() {
    #[cfg(windows)]
    if std::env::var_os("NYX_NO_ECOQOS").is_some() {
        opt_out_of_power_throttling();
    }
    let a: Vec<String> = std::env::args().skip(1).collect();
    let Some(file) = a.first() else {
        eprintln!("usage: workflow FILE.jsonl [--csv OUT.csv]");
        std::process::exit(2);
    };
    let csv = a.iter().position(|x| x == "--csv").and_then(|i| a.get(i + 1));
    let text = std::fs::read_to_string(file).unwrap_or_else(|e| {
        eprintln!("{file}: {e}");
        std::process::exit(2);
    });
    let mut out = csv.map(|p| std::fs::File::create(p).expect("the csv file"));
    if let Some(f) = out.as_mut() {
        let _ = writeln!(f, "line,cmd,note,seconds,since_cancel_s,stage,status,reply");
    }
    let mut engine = Engine::new();
    for n in &engine.library_notes {
        eprintln!("library: {n}");
    }
    let (mut last_job, mut cancelled_at): (Option<u64>, Option<Instant>) = (None, None);
    let t_all = Instant::now();
    println!("{:>4}  {:<18} {:<34} {:>9} {:>9}  status", "line", "command", "note", "seconds", "to cancel");
    for (n, raw) in text.lines().enumerate() {
        let raw = raw.trim();
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let line: Value = serde_json::from_str(raw).unwrap_or_else(|e| fail(n + 1, &format!("not JSON: {e}")));
        let cmd = line["cmd"].as_str().unwrap_or_else(|| fail(n + 1, "no `cmd`")).to_string();
        let mut args = line.get("args").cloned().unwrap_or(json!({}));
        let note = line["note"].as_str().unwrap_or("").to_string();
        if let Some(name) = line["job_arg"].as_str() {
            args[name] =
                json!(last_job.unwrap_or_else(|| fail(n + 1, "`job_arg`: no job was started before this line")));
        }
        let mut stage = String::new();
        if let Some(want) = line["after_stage"].as_str() {
            let job = last_job.unwrap_or_else(|| fail(n + 1, "`after_stage`: no job was started before this line"));
            let t = Instant::now();
            loop {
                let info = engine.jobs().info(job).unwrap_or_else(|| fail(n + 1, "`after_stage`: the job is gone"));
                stage = info.stage.clone();
                if info.state != nyx_engine::JobState::Running {
                    fail(
                        n + 1,
                        &format!("`after_stage`: the job ended before stage '{want}' (last stage seen: '{stage}')"),
                    );
                }
                if stage.contains(want) {
                    break;
                }
                if t.elapsed() > Duration::from_secs(900) {
                    fail(n + 1, &format!("`after_stage`: stage '{want}' not reached in 15 minutes"));
                }
                std::thread::sleep(Duration::from_millis(2));
            }
            if let Some(d) = line["delay_s"].as_f64() {
                // wait `d` more seconds inside the matched stage; the stage the job is in then is the one reported
                std::thread::sleep(Duration::from_secs_f64(d));
                let info = engine.jobs().info(job).unwrap_or_else(|| fail(n + 1, "`after_stage`: the job is gone"));
                if info.state != nyx_engine::JobState::Running {
                    fail(n + 1, "`delay_s`: the job ended during the delay");
                }
                stage = info.stage;
            }
        }
        let t0 = Instant::now();
        if cmd == "job_cancel" {
            cancelled_at = Some(t0);
        }
        let r = engine.call(&cmd, &args);
        let seconds = t0.elapsed().as_secs_f64();
        let mut since_cancel = String::new();
        if cmd == "job_wait" {
            if let Some(t) = cancelled_at.take() {
                since_cancel = format!("{:.3}", t.elapsed().as_secs_f64());
            }
        }
        if cmd == "job_cancel" && r.as_ref().is_ok_and(|v| v["cancelled"] != true) {
            fail(n + 1, "`job_cancel`: the job was no longer running, nothing was cancelled");
        }
        if line["expect_error"].as_bool() == Some(true) && r.is_ok() {
            fail(n + 1, "`expect_error`: the command succeeded (a cancelled job completed instead)");
        }
        let (status, reply) = match &r {
            Ok(v) => {
                if let Some(j) = v.get("job").and_then(Value::as_u64) {
                    last_job = Some(j);
                }
                ("ok".to_string(), v.to_string())
            }
            Err(e) => {
                if line["expect_error"].as_bool() != Some(true) {
                    println!("{:>4}  {:<18} {:<34} {:>9.3}  FAILED: {e}", n + 1, cmd, note, seconds);
                    std::process::exit(1);
                }
                (format!("error {}", e.code()), e.to_string())
            }
        };
        println!(
            "{:>4}  {:<18} {:<34} {:>9.3} {:>9}  {}{}",
            n + 1,
            cmd,
            note,
            seconds,
            since_cancel,
            status,
            if stage.is_empty() { String::new() } else { format!("  (stage then: {stage})") }
        );
        if let Some(f) = out.as_mut() {
            let q = |s: &str| format!("\"{}\"", s.replace('"', "\"\""));
            let short: String = reply.chars().take(300).collect();
            let _ = writeln!(
                f,
                "{},{},{},{:.4},{},{},{},{}",
                n + 1,
                cmd,
                q(&note),
                seconds,
                since_cancel,
                q(&stage),
                q(&status),
                q(&short)
            );
        }
    }
    println!("total {:.2} s", t_all.elapsed().as_secs_f64());
}

fn fail(line: usize, what: &str) -> ! {
    eprintln!("line {line}: {what}");
    std::process::exit(2);
}
