//! Command lists from the command line (productivity plan S1), without a window:
//!
//! ```text
//! nyx check FILE [--json]
//! nyx run FILE [--workspace DIR] [--no-check] [--full] [--json]
//! ```
//!
//! FILE is a JSON-lines list of engine commands (`{"cmd": ..., "args": {...}}`, see `nyx_engine::commandlist`). `check` reads it, checks every
//! command against the catalog and a dry run of the whole list against a copy of the engine, and prints what the list would mesh and solve; it
//! starts no job and writes nothing. `run` checks first (unless `--no-check`) and then sends the commands to the engine one after another, the same
//! engine commands the window and the MCP server use. With `--workspace` the solves are recorded as runs of a workflow in that folder.
//! Exit codes: 0 done, 1 a command failed, 2 the list is wrong or the usage is.

use nyx_engine::commandlist::{self, Problem};
use nyx_engine::Engine;
use serde_json::json;

pub type Exit = (i32, String);

fn usage(what: &str) -> Exit {
    (2, format!("{what}\nusage: nyx check FILE [--json]\n       nyx run FILE [--workspace DIR] [--no-check] [--full] [--json]"))
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
    std::fs::read_to_string(file).map_err(|e| (2, format!("{file}: {e}")))
}

fn print_problems(problems: &[Problem]) {
    for p in problems {
        eprintln!("{p}");
    }
}

/// `nyx check FILE`.
pub fn check(args: &[String]) -> Result<(), Exit> {
    let o = options(args)?;
    let text = read(&o.file)?;
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
