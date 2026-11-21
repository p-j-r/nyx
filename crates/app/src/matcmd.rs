// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Command-line material library & remote database manager (`nyx material ...`, ADR 101):
//! `nyx material list [--remote] [--domain DOMAIN] [--json]`
//! `nyx material show NAME [--json]`
//! `nyx material import FORMAT FILE [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]`
//! `nyx material fetch SOURCE [TARGET] [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]`

use nyx_engine::apiref;
use nyx_engine::material_fetch::{default_materials_folder, default_name_from_target};
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pub type Exit = (i32, String);

const USAGE: &str = "\
usage: nyx material list [--remote] [--domain DOMAIN] [--json]
       nyx material show NAME [--json]
       nyx material import FORMAT FILE [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]
       nyx material fetch SOURCE [TARGET] [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]";

fn usage(what: &str) -> Exit {
    (2, format!("{what}\n{USAGE}"))
}

pub fn run(args: &[String]) -> Result<(), Exit> {
    let Some((sub, rest)) = args.split_first() else {
        return Err(usage("missing subcommand"));
    };
    if matches!(sub.as_str(), "--help" | "-h" | "help") {
        println!("{USAGE}");
        return Ok(());
    }
    match sub.as_str() {
        "list" => cmd_list(rest),
        "show" => cmd_show(rest),
        "import" => cmd_import(rest),
        "fetch" => cmd_fetch(rest),
        other => Err(usage(&format!("unknown material subcommand '{other}'"))),
    }
}

fn take_flag_value(args: &[String], idx: &mut usize, flag: &str) -> Result<String, Exit> {
    if *idx + 1 >= args.len() {
        return Err(usage(&format!("{flag} requires a value")));
    }
    *idx += 1;
    Ok(args[*idx].clone())
}

fn cmd_list(args: &[String]) -> Result<(), Exit> {
    let mut remote = false;
    let mut json_out = false;
    let mut domain: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--remote" => remote = true,
            "--json" => json_out = true,
            "--domain" => domain = Some(take_flag_value(args, &mut i, "--domain")?),
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos => return Err(usage(&format!("unexpected argument '{pos}'"))),
        }
        i += 1;
    }

    let engine = Engine::new();
    if remote {
        let rows = engine.list_remote_materials(domain.as_deref()).map_err(|e| (2, e.to_string()))?;
        if json_out {
            println!("{}", serde_json::to_string_pretty(&rows).unwrap_or_default());
        } else if let Some(arr) = rows.as_array() {
            let w_id = arr.iter().filter_map(|r| r["id"].as_str().map(str::len)).max().unwrap_or(0);
            let w_dom = arr.iter().filter_map(|r| r["domain"].as_str().map(str::len)).max().unwrap_or(0);
            let w_fmt = arr.iter().filter_map(|r| r["format"].as_str().map(str::len)).max().unwrap_or(0);
            for r in arr {
                println!(
                    "{:<w_id$}  {:<w_dom$}  {:<w_fmt$}  {}  {}",
                    r["id"].as_str().unwrap_or(""),
                    r["domain"].as_str().unwrap_or(""),
                    r["format"].as_str().unwrap_or(""),
                    r["target"].as_str().unwrap_or(""),
                    r["summary"].as_str().unwrap_or("")
                );
            }
        }
        return Ok(());
    }

    if domain.is_some() {
        return Err(usage("--domain is only supported with --remote"));
    }
    let rows = engine.list_materials().map_err(|e| (1, e.to_string()))?;
    if json_out {
        println!("{}", serde_json::to_string_pretty(&rows).unwrap_or_default());
    } else if let Some(arr) = rows.as_array() {
        let width = arr.iter().filter_map(|r| r["name"].as_str().map(str::len)).max().unwrap_or(0);
        for r in arr {
            let mut tags = vec![];
            if r["optical"] == true {
                tags.push("optical");
            }
            if r["elastic"] == true {
                tags.push("elastic");
            }
            let desc = r["description"].as_str().unwrap_or("");
            println!("{:<width$}  {:<16}  {desc}", r["name"].as_str().unwrap_or(""), tags.join(","));
        }
    }
    Ok(())
}

fn cmd_show(args: &[String]) -> Result<(), Exit> {
    let mut json_out = false;
    let mut name: Option<String> = None;
    for a in args {
        match a.as_str() {
            "--json" => json_out = true,
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos if name.is_none() => name = Some(pos.to_string()),
            pos => return Err(usage(&format!("unexpected argument '{pos}'"))),
        }
    }
    let Some(name) = name else {
        return Err(usage("missing material NAME"));
    };
    let engine = Engine::new();
    match engine.describe_material(&name) {
        Ok(d) => {
            if json_out {
                println!("{}", serde_json::to_string_pretty(&d).unwrap_or_default());
            } else if let Some(toml) = d["toml"].as_str() {
                print!("{toml}");
            }
            Ok(())
        }
        Err(_) => {
            let names: Vec<String> = engine
                .list_materials()
                .ok()
                .and_then(|v| v.as_array().map(|a| a.iter().filter_map(|x| x["name"].as_str().map(str::to_string)).collect()))
                .unwrap_or_default();
            let near = apiref::nearest(&names, &name);
            let hint = if near.is_empty() {
                "`nyx material list` lists installed materials".to_string()
            } else {
                format!("did you mean {}?", near.join(", "))
            };
            Err((2, format!("no material '{name}': {hint}")))
        }
    }
}

struct WriteOpts {
    pos: Vec<String>,
    name: Option<String>,
    extends: Option<String>,
    description: Option<String>,
    save: bool,
    out: Option<PathBuf>,
    json: bool,
}

fn parse_write_opts(args: &[String], max_pos: usize) -> Result<WriteOpts, Exit> {
    let mut o = WriteOpts { pos: vec![], name: None, extends: None, description: None, save: false, out: None, json: false };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--name" => o.name = Some(take_flag_value(args, &mut i, "--name")?),
            "--extends" => o.extends = Some(take_flag_value(args, &mut i, "--extends")?),
            "--description" => o.description = Some(take_flag_value(args, &mut i, "--description")?),
            "--save" => o.save = true,
            "--out" => o.out = Some(PathBuf::from(take_flag_value(args, &mut i, "--out")?)),
            "--json" => o.json = true,
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos => {
                if o.pos.len() >= max_pos {
                    return Err(usage(&format!("unexpected argument '{pos}'")));
                }
                o.pos.push(pos.to_string());
            }
        }
        i += 1;
    }
    if o.save && o.out.is_some() {
        return Err(usage("pass either --save or --out PATH, not both"));
    }
    Ok(o)
}

fn is_dir_target(path: &Path) -> bool {
    if path.is_dir() {
        return true;
    }
    let s = path.to_string_lossy();
    s.ends_with('/') || s.ends_with('\\')
}

fn resolve_folder_arg(save: bool, out: Option<&Path>) -> Result<(Option<PathBuf>, Option<PathBuf>), Exit> {
    if save {
        let dir = default_materials_folder().ok_or_else(|| {
            (1, "no materials/ folder was found (set NYX_MATERIALS or pass --out PATH)".to_string())
        })?;
        return Ok((Some(dir), None));
    }
    if let Some(p) = out {
        if is_dir_target(p) {
            return Ok((Some(p.to_path_buf()), None));
        }
        return Ok((None, Some(p.to_path_buf())));
    }
    Ok((None, None))
}

fn emit_result(mut v: Value, file_out: Option<&Path>, json_out: bool) -> Result<(), Exit> {
    let emitted = v["imported_toml"].as_str().or_else(|| v["toml"].as_str()).unwrap_or("").to_string();
    if let Some(dest) = file_out {
        if let Some(parent) = dest.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent).map_err(|e| (1, format!("cannot create '{}': {e}", parent.display())))?;
        }
        std::fs::write(dest, &emitted).map_err(|e| (1, format!("cannot write '{}': {e}", dest.display())))?;
        if let Some(obj) = v.as_object_mut() {
            obj.insert("saved_path".into(), json!(dest.to_string_lossy().into_owned()));
        }
    }
    if json_out {
        println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
    } else if let Some(saved) = v["saved_path"].as_str() {
        println!("wrote {saved}");
    } else if !emitted.is_empty() {
        print!("{emitted}");
    }
    Ok(())
}

fn cmd_import(args: &[String]) -> Result<(), Exit> {
    let o = parse_write_opts(args, 2)?;
    if o.pos.len() != 2 {
        return Err(usage("import requires FORMAT and FILE"));
    }
    let format = &o.pos[0];
    let file = Path::new(&o.pos[1]);
    let content = std::fs::read_to_string(file)
        .map_err(|e| (1, format!("cannot read '{}': {e}", file.display())))?;
    let name = o.name.unwrap_or_else(|| default_name_from_target(&o.pos[1]));
    let (folder_arg, file_out) = resolve_folder_arg(o.save, o.out.as_deref())?;
    let mut engine = Engine::new();
    let v = engine
        .import_material(
            &name,
            format,
            &content,
            o.description.as_deref(),
            o.extends.as_deref(),
            false,
            folder_arg.as_deref(),
        )
        .map_err(|e| (1, e.to_string()))?;
    emit_result(v, file_out.as_deref(), o.json)
}

fn cmd_fetch(args: &[String]) -> Result<(), Exit> {
    let o = parse_write_opts(args, 2)?;
    if o.pos.is_empty() {
        return Err(usage("fetch requires SOURCE (curated preset id, `rii`, `aflow`, or `numbat`)"));
    }
    let source = &o.pos[0];
    let target = o.pos.get(1).map(String::as_str);
    let (folder_arg, file_out) = resolve_folder_arg(o.save, o.out.as_deref())?;
    let mut engine = Engine::new();
    let v = engine
        .fetch_material(
            source,
            target,
            o.name.as_deref(),
            o.description.as_deref(),
            o.extends.as_deref(),
            false,
            folder_arg.as_deref(),
        )
        .map_err(|e| (1, e.to_string()))?;
    emit_result(v, file_out.as_deref(), o.json)
}
