// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Command-line material library & remote database manager (`nyx material ...`, ADR 101, ADR 102):
//! `nyx material list [--remote] [--domain DOMAIN] [--json]`
//! `nyx material show NAME [--json]`
//! `nyx material eval NAME [--temp K] [--freq HZ] [--wavelength M] [--e-field V_PER_M] [--thickness M] [--film-model MODEL] [--json]`
//! `nyx material bcs [MATERIAL] --freq HZ --temp K [--thickness M] [--film-model MODEL] [--tc K] [--gap-0 J] [--lambda-l M] [--xi-0 M] [--ell M] [--sigma-n S_PER_M] [--dynes-gamma J] [--r-res OHM_SQ] [--x-qp-floor X] [--json]`
//! `nyx material search QUERY [--domain DOMAIN] [--source SOURCE] [--live] [--api-key KEY] [--json]`
//! `nyx material import FORMAT FILE [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]`
//! `nyx material fetch SOURCE [TARGET] [--name NAME] [--extends BASE] [--description TEXT] [--api-key KEY] [--save | --out PATH] [--json]`

use nyx_engine::apiref;
use nyx_engine::material_fetch::{default_materials_folder, default_name_from_target};
use nyx_engine::Engine;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

pub type Exit = (i32, String);

const USAGE: &str = "\
usage: nyx material list [--remote] [--domain DOMAIN] [--json]
       nyx material show NAME [--json]
       nyx material eval NAME [--temp K] [--freq HZ] [--wavelength M] [--e-field V_PER_M] [--thickness M] [--film-model MODEL] [--json]
       nyx material bcs [MATERIAL] --freq HZ --temp K [--thickness M] [--film-model MODEL] [--tc K] [--gap-0 J] [--lambda-l M] [--xi-0 M] [--ell M] [--sigma-n S_PER_M] [--dynes-gamma J] [--r-res OHM_SQ] [--x-qp-floor X] [--json]
       nyx material search QUERY [--domain DOMAIN] [--source SOURCE] [--live] [--api-key KEY] [--json]
       nyx material import FORMAT FILE [--name NAME] [--extends BASE] [--description TEXT] [--save | --out PATH] [--json]
       nyx material fetch SOURCE [TARGET] [--name NAME] [--extends BASE] [--description TEXT] [--api-key KEY] [--save | --out PATH] [--json]";

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
        "eval" => cmd_eval(rest),
        "bcs" => cmd_bcs(rest),
        "search" => cmd_search(rest),
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
            if r["dielectric"] == true {
                tags.push("dielectric");
            }
            if r["superconducting"] == true {
                tags.push("sc");
            }
            if r["piezoelectric"] == true {
                tags.push("piezo");
            }
            if r["photoelastic"] == true {
                tags.push("pe");
            }
            if r["electro_optic"] == true {
                tags.push("eo");
            }
            if r["viscoelastic"] == true {
                tags.push("visc");
            }
            let desc = r["description"].as_str().unwrap_or("");
            println!("{:<width$}  {:<24}  {desc}", r["name"].as_str().unwrap_or(""), tags.join(","));
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
    api_key: Option<String>,
    save: bool,
    out: Option<PathBuf>,
    json: bool,
}

fn parse_write_opts(args: &[String], max_pos: usize) -> Result<WriteOpts, Exit> {
    let mut o = WriteOpts {
        pos: vec![],
        name: None,
        extends: None,
        description: None,
        api_key: None,
        save: false,
        out: None,
        json: false,
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--name" => o.name = Some(take_flag_value(args, &mut i, "--name")?),
            "--extends" => o.extends = Some(take_flag_value(args, &mut i, "--extends")?),
            "--description" => o.description = Some(take_flag_value(args, &mut i, "--description")?),
            "--api-key" => o.api_key = Some(take_flag_value(args, &mut i, "--api-key")?),
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

fn parse_f64_flag(args: &[String], idx: &mut usize, flag: &str) -> Result<f64, Exit> {
    let s = take_flag_value(args, idx, flag)?;
    s.parse::<f64>().map_err(|_| usage(&format!("{flag} expects a number, got '{s}'")))
}

fn cmd_eval(args: &[String]) -> Result<(), Exit> {
    let mut name: Option<String> = None;
    let mut temp: Option<f64> = None;
    let mut freq: Option<f64> = None;
    let mut wavelength: Option<f64> = None;
    let mut e_field: Option<f64> = None;
    let mut thickness: Option<f64> = None;
    let mut film_model: Option<String> = None;
    let mut json_out = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--temp" | "-T" => temp = Some(parse_f64_flag(args, &mut i, "--temp")?),
            "--freq" | "-f" => freq = Some(parse_f64_flag(args, &mut i, "--freq")?),
            "--wavelength" | "-w" => wavelength = Some(parse_f64_flag(args, &mut i, "--wavelength")?),
            "--e-field" => e_field = Some(parse_f64_flag(args, &mut i, "--e-field")?),
            "--thickness" | "-d" => thickness = Some(parse_f64_flag(args, &mut i, "--thickness")?),
            "--film-model" => film_model = Some(take_flag_value(args, &mut i, "--film-model")?),
            "--json" => json_out = true,
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos if name.is_none() => name = Some(pos.to_string()),
            pos => return Err(usage(&format!("unexpected argument '{pos}'"))),
        }
        i += 1;
    }
    let Some(name) = name else {
        return Err(usage("eval requires material NAME"));
    };
    let engine = Engine::new();
    let v = engine
        .evaluate_material(
            &name,
            temp,
            freq,
            wavelength,
            e_field,
            thickness,
            film_model.as_deref(),
            None,
        )
        .map_err(|e| (1, e.to_string()))?;
    if json_out {
        println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
    } else {
        println!("{} (format {}, T = {} K)", v["name"].as_str().unwrap_or(&name), v["format"], v["temperature_k"]);
        if let Some(th) = v.get("thermal") {
            println!("  thermal: Cp = {} J/(kg*K), k = {}", th["specific_heat_j_kgk"], th["conductivity_w_mk"]);
        }
        if let Some(dl) = v.get("dielectric") {
            println!("  dielectric: eps_re = {}, tan_delta_xx = {:e}", dl["eps_re"], dl["tan_delta_xx"].as_f64().unwrap_or(0.0));
        }
        if let Some(opt) = v.get("optical") {
            println!("  optical: n = {}, k = {}", opt["n"], opt["k"]);
        }
        if let Some(sc) = v.get("superconducting") {
            println!(
                "  superconducting ({}): R_s = {:.4e} Ohm/sq, L_k = {:.4e} H/sq, lambda_eff = {:.2e} m, Delta(T) = {:.4} meV",
                sc["film_model"].as_str().unwrap_or("kerr"),
                sc["r_s_ohm_sq"].as_f64().unwrap_or(0.0),
                sc["l_k_sq_h_sq"].as_f64().unwrap_or(0.0),
                sc["lambda_eff_m"].as_f64().unwrap_or(0.0),
                sc["delta_t_mev"].as_f64().unwrap_or(0.0),
            );
        }
    }
    Ok(())
}

fn cmd_bcs(args: &[String]) -> Result<(), Exit> {
    let mut material: Option<String> = None;
    let mut freq: Option<f64> = None;
    let mut temp: Option<f64> = None;
    let mut thickness: Option<f64> = None;
    let mut film_model: Option<String> = None;
    let mut tc: Option<f64> = None;
    let mut gap_0: Option<f64> = None;
    let mut lambda_l: Option<f64> = None;
    let mut xi_0: Option<f64> = None;
    let mut ell: Option<f64> = None;
    let mut sigma_n: Option<f64> = None;
    let mut dynes_gamma: Option<f64> = None;
    let mut r_res: Option<f64> = None;
    let mut x_qp_floor: Option<f64> = None;
    let mut json_out = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--freq" | "-f" => freq = Some(parse_f64_flag(args, &mut i, "--freq")?),
            "--temp" | "-T" => temp = Some(parse_f64_flag(args, &mut i, "--temp")?),
            "--thickness" | "-d" => thickness = Some(parse_f64_flag(args, &mut i, "--thickness")?),
            "--film-model" => film_model = Some(take_flag_value(args, &mut i, "--film-model")?),
            "--tc" => tc = Some(parse_f64_flag(args, &mut i, "--tc")?),
            "--gap-0" => gap_0 = Some(parse_f64_flag(args, &mut i, "--gap-0")?),
            "--lambda-l" => lambda_l = Some(parse_f64_flag(args, &mut i, "--lambda-l")?),
            "--xi-0" => xi_0 = Some(parse_f64_flag(args, &mut i, "--xi-0")?),
            "--ell" | "--mean-free-path" => ell = Some(parse_f64_flag(args, &mut i, "--ell")?),
            "--sigma-n" => sigma_n = Some(parse_f64_flag(args, &mut i, "--sigma-n")?),
            "--dynes-gamma" => dynes_gamma = Some(parse_f64_flag(args, &mut i, "--dynes-gamma")?),
            "--r-res" => r_res = Some(parse_f64_flag(args, &mut i, "--r-res")?),
            "--x-qp-floor" => x_qp_floor = Some(parse_f64_flag(args, &mut i, "--x-qp-floor")?),
            "--json" => json_out = true,
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos if material.is_none() => material = Some(pos.to_string()),
            pos => return Err(usage(&format!("unexpected argument '{pos}'"))),
        }
        i += 1;
    }
    let Some(f_hz) = freq else {
        return Err(usage("bcs requires --freq HZ"));
    };
    let Some(t_k) = temp else {
        return Err(usage("bcs requires --temp K"));
    };
    let mut obj = serde_json::Map::new();
    obj.insert("frequency".into(), json!(f_hz));
    obj.insert("temperature".into(), json!(t_k));
    if let Some(m) = material {
        obj.insert("material".into(), json!(m));
    }
    if let Some(v) = thickness {
        obj.insert("thickness".into(), json!(v));
    }
    if let Some(v) = film_model {
        obj.insert("film_model".into(), json!(v));
    }
    if let Some(v) = tc {
        obj.insert("tc".into(), json!(v));
    }
    if let Some(v) = gap_0 {
        obj.insert("gap_0".into(), json!(v));
    }
    if let Some(v) = lambda_l {
        obj.insert("lambda_L".into(), json!(v));
    }
    if let Some(v) = xi_0 {
        obj.insert("xi_0".into(), json!(v));
    }
    if let Some(v) = ell {
        obj.insert("ell".into(), json!(v));
    }
    if let Some(v) = sigma_n {
        obj.insert("sigma_n".into(), json!(v));
    }
    if let Some(v) = dynes_gamma {
        obj.insert("dynes_gamma".into(), json!(v));
    }
    if let Some(v) = r_res {
        obj.insert("r_res".into(), json!(v));
    }
    if let Some(v) = x_qp_floor {
        obj.insert("x_qp_floor".into(), json!(v));
    }
    let engine = Engine::new();
    let v = engine.bcs_surface_impedance_cmd(&Value::Object(obj)).map_err(|e| (1, e.to_string()))?;
    if json_out {
        println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
    } else {
        println!(
            "BCS ({}, f = {:.3e} Hz, T = {:.4} K, d = {:.2e} m): R_s = {:.6e} Ohm/sq, L_k = {:.6e} H/sq, lambda_eff = {:.3e} m, Delta(T) = {:.4} meV",
            v["film_model"].as_str().unwrap_or("kerr"),
            v["frequency_hz"].as_f64().unwrap_or(f_hz),
            v["temperature_k"].as_f64().unwrap_or(t_k),
            v["thickness_m"].as_f64().unwrap_or(0.0),
            v["r_s_ohm_sq"].as_f64().unwrap_or(0.0),
            v["l_k_sq_h_sq"].as_f64().unwrap_or(0.0),
            v["lambda_eff_m"].as_f64().unwrap_or(0.0),
            v["delta_t_mev"].as_f64().unwrap_or(0.0),
        );
    }
    Ok(())
}

fn cmd_search(args: &[String]) -> Result<(), Exit> {
    let mut query: Option<String> = None;
    let mut domain: Option<String> = None;
    let mut source: Option<String> = None;
    let mut live = false;
    let mut api_key: Option<String> = None;
    let mut json_out = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--domain" => domain = Some(take_flag_value(args, &mut i, "--domain")?),
            "--source" => source = Some(take_flag_value(args, &mut i, "--source")?),
            "--live" => live = true,
            "--api-key" => api_key = Some(take_flag_value(args, &mut i, "--api-key")?),
            "--json" => json_out = true,
            flag if flag.starts_with('-') => return Err(usage(&format!("unknown option {flag}"))),
            pos if query.is_none() => query = Some(pos.to_string()),
            pos => return Err(usage(&format!("unexpected argument '{pos}'"))),
        }
        i += 1;
    }
    let Some(q) = query else {
        return Err(usage("search requires QUERY"));
    };
    let engine = Engine::new();
    let v = engine
        .search_remote_materials(&q, domain.as_deref(), source.as_deref(), live, api_key.as_deref())
        .map_err(|e| (1, e.to_string()))?;
    if json_out {
        println!("{}", serde_json::to_string_pretty(&v).unwrap_or_default());
    } else {
        if let Some(arr) = v["catalog_matches"].as_array() {
            for r in arr {
                println!(
                    "{:<28}  {:<15}  {:<6}  {}  {}",
                    r["id"].as_str().unwrap_or(""),
                    r["domain"].as_str().unwrap_or(""),
                    r["format"].as_str().unwrap_or(""),
                    r["target"].as_str().unwrap_or(""),
                    r["summary"].as_str().unwrap_or("")
                );
            }
        }
        if let Some(arr) = v["live_matches"].as_array() {
            for r in arr {
                println!(
                    "[live:{}] {:<24}  {}  {}",
                    r["format"].as_str().unwrap_or(""),
                    r["id"].as_str().unwrap_or(""),
                    r["target"].as_str().unwrap_or(""),
                    r["summary"].as_str().unwrap_or("")
                );
            }
        }
    }
    Ok(())
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
        .fetch_material_with_key(
            source,
            target,
            o.name.as_deref(),
            o.description.as_deref(),
            o.extends.as_deref(),
            false,
            folder_arg.as_deref(),
            o.api_key.as_deref(),
        )
        .map_err(|e| (1, e.to_string()))?;
    emit_result(v, file_out.as_deref(), o.json)
}
