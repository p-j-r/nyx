//! `nyx-connect`: the command-line twin of the app's Connect menu. It writes nothing unless you say connect/disconnect.
//! `--json` prints the versioned report of `nyx_connect::report` instead of text (contract `connect_report`).
use nyx_connect::report::{change_report, change_row, list_report, test_report};
use nyx_connect::{claude_code_command, connect, disconnect, list, server_exe, snippet, test_server, Env, Status};

fn usage() {
    println!(
        "nyx-connect: connect the nyx MCP server to your AI apps\n\n  nyx-connect list                 show every supported app and whether nyx is connected\n  nyx-connect connect <app|all>    add nyx to one app, or to every installed app ('all')\n  nyx-connect disconnect <app>     remove nyx from one app\n  nyx-connect test                 start the server and list its tools\n  nyx-connect snippet              print the JSON to paste into any other app\n\nAdd --json to list, connect, disconnect or test for a machine-readable report.\nApps: claude-code, claude-desktop, cursor, vscode, windsurf. Set NYX_MCP to use a nyx-mcp program in another place."
    );
}

fn describe(s: &Status) -> String {
    match s {
        Status::NotInstalled => "not installed".into(),
        Status::NotConnected => "installed, nyx not connected".into(),
        Status::Connected => "connected".into(),
        Status::Stale { found } => format!("points at another program ({found}); run connect to update"),
        Status::Unreadable { why } => format!("cannot edit its config: {why}"),
    }
}

fn print_json(v: &serde_json::Value) {
    println!("{}", serde_json::to_string_pretty(v).unwrap_or_default());
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let json = args.iter().any(|a| a == "--json");
    args.retain(|a| a != "--json");
    let exe = server_exe();
    let env = match Env::detect() {
        Ok(e) => e,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };
    let code = match args.first().map(String::as_str) {
        Some("list") => {
            let hosts = list(&env, &exe);
            if json {
                print_json(&list_report(&exe, &hosts));
            } else {
                println!("server: {}{}", exe.display(), if exe.is_file() { "" } else { "  (NOT FOUND, build it first)" });
                for h in hosts {
                    println!("{:15} {:16} {}", h.id, h.name, describe(&h.status));
                }
            }
            0
        }
        Some("connect") => match args.get(1).map(String::as_str) {
            None => {
                eprintln!("connect which app? (see 'nyx-connect list', or 'all')");
                2
            }
            Some(which) => {
                let ids: Vec<&str> = if which == "all" {
                    let l = list(&env, &exe);
                    l.iter().filter(|h| !matches!(h.status, Status::NotInstalled | Status::Unreadable { .. })).map(|h| h.id).collect()
                } else {
                    vec![which]
                };
                let results: Vec<_> = ids.iter().map(|id| (*id, connect(&env, id, &exe))).collect();
                if json {
                    print_json(&change_report("connect", results.iter().map(|(id, r)| change_row(id, r)).collect()));
                } else {
                    for (id, r) in &results {
                        match r {
                            Ok(d) if d.changed => println!("{id}: connected ({}){}", d.config.display(), d.backup.as_ref().map(|b| format!(", original kept as {}", b.display())).unwrap_or_default()),
                            Ok(_) => println!("{id}: already connected"),
                            Err(e) => eprintln!("{id}: {e}"),
                        }
                    }
                    if ids.is_empty() {
                        println!("no installed app found; use 'nyx-connect snippet' for the JSON");
                    }
                }
                i32::from(results.iter().any(|(_, r)| r.is_err()))
            }
        },
        Some("disconnect") => match args.get(1) {
            None => {
                eprintln!("disconnect which app?");
                2
            }
            Some(id) => {
                let r = disconnect(&env, id);
                if json {
                    print_json(&change_report("disconnect", vec![change_row(id, &r)]));
                } else {
                    match &r {
                        Ok(d) if d.changed => println!("{id}: removed ({})", d.config.display()),
                        Ok(_) => println!("{id}: nyx was not connected"),
                        Err(e) => eprintln!("{id}: {e}"),
                    }
                }
                i32::from(r.is_err())
            }
        },
        Some("test") => {
            let r = test_server(&exe);
            if json {
                print_json(&test_report(&r));
            } else {
                match &r {
                    Ok(i) => println!("{} {} (protocol {}), {} tools: {}", i.name, i.version, i.protocol, i.tools.len(), i.tools.join(", ")),
                    Err(e) => eprintln!("{e}"),
                }
            }
            i32::from(r.is_err())
        }
        Some("snippet") => {
            println!("{}\n\nClaude Code command:\n{}", snippet(&exe), claude_code_command(&exe));
            0
        }
        Some("--help") | Some("-h") | None => {
            usage();
            0
        }
        Some(o) => {
            eprintln!("unknown command '{o}'");
            usage();
            2
        }
    };
    std::process::exit(code);
}
