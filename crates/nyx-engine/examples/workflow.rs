//! Replays a list of engine commands against an `Engine` and times each one (productivity plan P0). The list is a JSON-lines file, one command per
//! line, the format `nyx check` and `nyx run` take (plan step S1, `nyx_engine::commandlist`, which this example uses to replay); the workflows of `tools/workflows/` are written in it.
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

use nyx_engine::commandlist;
use nyx_engine::Engine;
use std::io::Write;
use std::time::Instant;

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
    // the list is read, checked and replayed by the engine's command list module (the same code `nyx run` uses)
    let (lines, problems) = commandlist::parse(&text);
    if !problems.is_empty() {
        for p in &problems {
            eprintln!("{p}");
        }
        std::process::exit(2);
    }
    let mut out = csv.map(|p| std::fs::File::create(p).expect("the csv file"));
    if let Some(f) = out.as_mut() {
        let _ = writeln!(f, "line,cmd,note,seconds,since_cancel_s,stage,status,reply");
    }
    let mut engine = Engine::new();
    for n in &engine.library_notes {
        eprintln!("library: {n}");
    }
    let t_all = Instant::now();
    println!("{:>4}  {:<18} {:<34} {:>9} {:>9}  status", "line", "command", "note", "seconds", "to cancel");
    let result = commandlist::run(&mut engine, &lines, &mut |d| {
        let note = d.line.note.clone().unwrap_or_default();
        let since_cancel = d.since_cancel.map(|s| format!("{s:.3}")).unwrap_or_default();
        let (status, reply) = match d.reply {
            Ok(v) => ("ok".to_string(), v.to_string()),
            Err(e) if d.line.expect_error => (format!("error {}", e.code()), e.to_string()),
            Err(e) => (format!("FAILED: {e}"), e.to_string()),
        };
        println!(
            "{:>4}  {:<18} {:<34} {:>9.3} {:>9}  {}{}",
            d.line.number,
            d.line.cmd,
            note,
            d.seconds,
            since_cancel,
            status,
            if d.stage.is_empty() { String::new() } else { format!("  (stage then: {})", d.stage) }
        );
        if let Some(f) = out.as_mut() {
            let q = |s: &str| format!("\"{}\"", s.replace('"', "\"\""));
            let short: String = reply.chars().take(300).collect();
            let _ = writeln!(
                f,
                "{},{},{},{:.4},{},{},{},{}",
                d.line.number,
                d.line.cmd,
                q(&note),
                d.seconds,
                since_cancel,
                q(d.stage),
                q(&status),
                q(&short)
            );
        }
    });
    if let Err(e) = result {
        // a command that failed was printed as a row already; a replay key that could not be honoured is reported here
        eprintln!("{e}");
        std::process::exit(if e.message.starts_with('`') { 2 } else { 1 });
    }
    println!("total {:.2} s", t_all.elapsed().as_secs_f64());
}
