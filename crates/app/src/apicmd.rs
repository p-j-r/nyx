//! The API reference from the command line (productivity plan S2): `nyx api list [--json]` lists the engine commands and `nyx api show COMMAND
//! [--json]` prints one command's page, the same text as `docs/api/COMMAND.md`. Both are generated from `Engine::catalog()` by
//! `nyx_engine::apiref` (ADR 57); nothing here knows a command. Exit code 2 for a command that does not exist (with the nearest names) or a
//! usage error.

use nyx_engine::apiref;
use nyx_engine::Engine;
use serde_json::json;

pub type Exit = (i32, String);

fn usage(what: &str) -> Exit {
    (2, format!("{what}\nusage: nyx api list [--json]\n       nyx api show COMMAND [--json]"))
}

pub fn run(args: &[String]) -> Result<(), Exit> {
    let json = args.iter().any(|a| a == "--json");
    let rest: Vec<&String> = args.iter().filter(|a| *a != "--json").collect();
    if let Some(o) = rest.iter().find(|a| a.starts_with("--")) {
        return Err(usage(&format!("unknown option {o}")));
    }
    let mut cmds = Engine::new().catalog();
    cmds.sort_by(|a, b| a.name.cmp(&b.name));
    match rest.first().map(|s| s.as_str()) {
        Some("list") if rest.len() == 1 => {
            if json {
                let v: Vec<_> = cmds
                    .iter()
                    .map(|c| json!({ "name": c.name, "title": c.title, "access": c.access }))
                    .collect();
                println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
            } else {
                let width = cmds.iter().map(|c| c.name.len()).max().unwrap_or(0);
                for c in &cmds {
                    println!("{:<width$}  {:<9}  {}", c.name, serde_json::to_value(c.access).unwrap_or_default().as_str().unwrap_or(""), c.title);
                }
            }
            Ok(())
        }
        Some("show") if rest.len() == 2 => {
            let name = rest[1].as_str();
            let Some(c) = cmds.iter().find(|c| c.name == name) else {
                let names: Vec<String> = cmds.iter().map(|c| c.name.clone()).collect();
                let near = apiref::nearest(&names, name);
                let hint = if near.is_empty() { "`nyx api list` lists them".to_string() } else { format!("did you mean {}?", near.join(", ")) };
                return Err((2, format!("there is no command '{name}': {hint}")));
            };
            if json {
                let v = json!({
                    "name": c.name, "title": c.title, "description": c.description, "access": c.access, "schema": c.schema,
                    "example": { "cmd": c.name, "args": apiref::example_args(c) },
                });
                println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
            } else {
                print!("{}", apiref::command_page(c));
            }
            Ok(())
        }
        Some("stubs") if rest.len() == 1 => {
            // the type stub of the `nyx` object of scripts, for an editor (save it as nyx.pyi next to the script) and for `nyx check`
            print!("{}", nyx_script::stubs::nyx_pyi(&Engine::new()));
            Ok(())
        }
        Some("list") | Some("show") => Err(usage("wrong number of arguments")),
        Some(other) => Err(usage(&format!("unknown subcommand '{other}'"))),
        None => Err(usage("no subcommand")),
    }
}
