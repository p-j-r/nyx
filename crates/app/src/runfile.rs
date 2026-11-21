// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Command lists from the command line (productivity plan S1), without a window:
//!
//! ```text
//! nyx check FILE [--json]
//! nyx run FILE [--workspace DIR] [--no-check] [--full] [--json] [--seconds N] [--memory MB] [--calls N]
//! ```
//!
//! FILE is a JSON-lines list of engine commands (`{"cmd": ..., "args": {...}}`, see `nyx_engine::commandlist`). `check` reads it, checks every
//! command against the catalog and a dry run of the whole list against a copy of the engine, and prints what the list would mesh and solve; it
//! starts no job and writes nothing. `run` checks first (unless `--no-check`) and then sends the commands to the engine one after another, the same
//! engine commands the window and the MCP server use. With `--workspace` the solves are recorded as runs of a workflow in that folder.
//! Exit codes: 0 done, 1 a command failed, 2 the list is wrong or the usage is.
//!
//! A file ending in `.nys` is a **script** (productivity plan S4): Python syntax, run by Monty in a worker process (`nyx-script-worker`, which has to
//! sit next to this program). `check` parses it and checks every command and argument name against the catalog without running anything; `run` does
//! that first and then runs it. A script's files are confined to `--workspace` (else to the script's folder); `--seconds`, `--memory` (megabytes) and
//! `--calls` raise the strict default limits for one run.

use nyx_engine::commandlist::{self, Problem};
use nyx_engine::Engine;
use serde_json::json;

pub type Exit = (i32, String);

fn usage(what: &str) -> Exit {
    (2, format!("{what}\nusage: nyx check FILE [--json]\n       nyx run FILE [--workspace DIR] [--no-check] [--full] [--json] [--seconds N] [--memory MB] [--calls N]\n       (FILE.nys is a script, any other file a JSON-lines command list; the limits apply to scripts)"))
}

struct Opts {
    file: String,
    workspace: Option<String>,
    no_check: bool,
    full: bool,
    json: bool,
}

fn options(args: &[String]) -> Result<Opts, Exit> {
    let mut o = Opts { file: String::new(), workspace: None, no_check: false, full: false, json: false };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--workspace" => {
                i += 1;
                o.workspace = Some(args.get(i).ok_or_else(|| usage("--workspace needs a folder"))?.clone());
            }
            "--no-check" => o.no_check = true,
            // limits of a script run; their values are read by `script_limits`
            "--seconds" | "--memory" | "--calls" => {
                i += 1;
                if i >= args.len() {
                    return Err(usage(&format!("{} needs a number", args[i - 1])));
                }
            }
            "--full" => o.full = true,
            "--json" => o.json = true,
            a if a.starts_with("--") => return Err(usage(&format!("unknown option {a}"))),
            a if o.file.is_empty() => o.file = a.to_string(),
            a => return Err(usage(&format!("one file only, not also {a}"))),
        }
        i += 1;
    }
    if o.file.is_empty() {
        return Err(usage("no file"));
    }
    Ok(o)
}

fn read(file: &str) -> Result<String, Exit> {
    std::fs::read_to_string(file)
        .map(|s| s.strip_prefix('\u{feff}').map(str::to_string).unwrap_or(s))
        .map_err(|e| (2, format!("{file}: {e}")))
}

fn print_problems(problems: &[Problem]) {
    for p in problems {
        eprintln!("{p}");
    }
}

/// A script (a `*.nys` file) is Python syntax run by Monty in a worker process; any other file is a JSON-lines command list.
fn is_script(file: &str) -> bool {
    file.ends_with(".nys")
}

/// The limits of a script run: the strict defaults, raised per run by `--seconds N`, `--memory MB` and `--calls N`.
fn script_limits(args: &[String]) -> Result<nyx_script::Limits, Exit> {
    let mut l = nyx_script::Limits::default();
    let mut i = 0;
    while i < args.len() {
        let want = |i: usize| args.get(i + 1).ok_or_else(|| usage(&format!("{} needs a number", args[i])));
        match args[i].as_str() {
            "--seconds" => l.seconds = want(i)?.parse().map_err(|_| usage("--seconds needs a number"))?,
            "--memory" => l.memory_bytes = want(i)?.parse::<usize>().map_err(|_| usage("--memory needs a number of megabytes"))? << 20,
            "--calls" => l.calls = want(i)?.parse().map_err(|_| usage("--calls needs a whole number"))?,
            _ => {
                i += 1;
                continue;
            }
        }
        i += 2;
    }
    Ok(l)
}

/// Where a script may read and write files: the workspace when one is given, else the folder of the script.
fn script_root(o: &Opts) -> std::path::PathBuf {
    match &o.workspace {
        Some(w) => std::path::PathBuf::from(w),
        None => std::path::Path::new(&o.file).parent().filter(|p| !p.as_os_str().is_empty()).map(std::path::Path::to_path_buf).unwrap_or_else(|| std::path::PathBuf::from(".")),
    }
}

fn check_script(o: &Opts, text: &str, args: &[String]) -> Result<(), Exit> {
    let mut engine = Engine::new();
    let opts = nyx_script::Options::new(script_root(o));
    nyx_script::host::check(&mut engine, &o.file, text, &script_limits(args)?, &opts).map_err(|e| (2, format!("{}\n{}: the script failed the check; nothing was run", e, o.file)))?;
    println!("{}: parses and matches the nyx commands; nothing was run (a script's plan, which commands it sends and what they cost, depends on its results and is known only by running it)", o.file);
    Ok(())
}

fn run_script(o: &Opts, text: &str, args: &[String]) -> Result<(), Exit> {
    let mut engine = Engine::new();
    for n in &engine.library_notes {
        eprintln!("library: {n}");
    }
    if let Some(w) = &o.workspace {
        engine.set_workspace(std::path::PathBuf::from(w), true);
    }
    let opts = nyx_script::Options::new(script_root(o));
    let started = std::time::Instant::now();
    // a script that fails the static check is a wrong script (2, as for a bad command list); one that fails while running is a failed run (1)
    let out = nyx_script::host::run(&mut engine, &o.file, text, &script_limits(args)?, &opts).map_err(|e| (if e.starts_with("the script failed the static check") { 2 } else { 1 }, e))?;
    if o.json {
        println!("{}", json!({ "value": out.value, "printed": out.printed, "calls": out.calls, "seconds": started.elapsed().as_secs_f64() }));
    } else {
        for l in &out.printed {
            println!("{l}");
        }
        if !out.value.is_null() {
            println!("{}", serde_json::to_string_pretty(&out.value).unwrap_or_default());
        }
        println!("{} engine commands in {:.2} s", out.calls, started.elapsed().as_secs_f64());
    }
    Ok(())
}

/// `nyx check FILE`.
pub fn check(args: &[String]) -> Result<(), Exit> {
    let o = options(args)?;
    let text = read(&o.file)?;
    if is_script(&o.file) {
        return check_script(&o, &text, args);
    }
    let engine = Engine::new();
    let r = engine.check_list(&text);
    if o.json {
        println!("{}", serde_json::to_string_pretty(&r).unwrap_or_default());
    } else if r.ok() {
        println!("{}: {} commands, hash {}", o.file, r.commands, r.hash);
        for s in &r.steps {
            println!("{:>4}  {:<18} {:<11} {}", s.line, s.cmd, s.kind, s.summary);
        }
        println!(
            "plan: {} mesh build(s) and {} solve job(s); nothing was meshed, solved or written (cost estimates are plan step Q2)",
            r.meshes, r.solves
        );
    }
    if r.ok() {
        Ok(())
    } else {
        if !o.json {
            print_problems(&r.problems);
        }
        Err((2, format!("{}: {} problem(s); nothing was run", o.file, r.problems.len())))
    }
}

/// `nyx run FILE`.
pub fn run(args: &[String]) -> Result<(), Exit> {
    let o = options(args)?;
    let text = read(&o.file)?;
    if is_script(&o.file) {
        return run_script(&o, &text, args);
    }
    let mut engine = Engine::new();
    if !o.no_check {
        let r = engine.check_list(&text);
        if !r.ok() {
            print_problems(&r.problems);
            return Err((2, format!("{}: {} problem(s); nothing was run (--no-check skips the check)", o.file, r.problems.len())));
        }
    }
    let (lines, problems) = commandlist::parse(&text);
    if !problems.is_empty() {
        print_problems(&problems);
        return Err((2, format!("{}: {} problem(s); nothing was run", o.file, problems.len())));
    }
    for n in &engine.library_notes {
        eprintln!("library: {n}");
    }
    if let Some(w) = &o.workspace {
        engine.set_workspace(std::path::PathBuf::from(w), true);
    }
    let started = std::time::Instant::now();
    let result = commandlist::run(&mut engine, &lines, &mut |d| {
        if o.json {
            let mut v = json!({ "line": d.line.number, "cmd": d.line.cmd, "seconds": d.seconds });
            match d.reply {
                Ok(r) => v["reply"] = r.clone(),
                Err(e) => {
                    v["error"] = json!({ "code": e.code(), "message": e.to_string() });
                }
            }
            println!("{v}");
            return;
        }
        let (status, reply) = match d.reply {
            Ok(v) => ("ok".to_string(), v.to_string()),
            Err(e) => (format!("error {}", e.code()), e.to_string()),
        };
        let note = d.line.note.as_deref().unwrap_or("");
        println!("{:>4}  {:<18} {:>8.3} s  {}{}", d.line.number, d.line.cmd, d.seconds, status, if note.is_empty() { String::new() } else { format!("  ({note})") });
        let shown = if o.full || reply.chars().count() <= 240 { reply } else { format!("{}...", reply.chars().take(240).collect::<String>()) };
        println!("        {shown}");
    });
    if !o.json {
        println!("{} commands in {:.2} s", lines.len(), started.elapsed().as_secs_f64());
    }
    result.map_err(|e| (1, e.to_string()))
}
