//! Headless runs without a window, for workstations and clusters:
//!
//! ```text
//! nyx --batch model.nyx.json [--physics ems,sms]
//! ```
//!
//! The project is read, meshed and solved by the engine, the same code as the window and the MCP server, and the summary and the stage log
//! of every solved physics are printed. Without `--physics` the active physics is solved (with those it takes fields from).

use nyx_engine::facade::physics::Scalar;
use nyx_engine::Engine;
use serde_json::json;

pub fn run(args: &[String]) -> Result<(), String> {
    let mut path = None;
    let mut physics: Option<Vec<String>> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--physics" => {
                i += 1;
                let v = args.get(i).ok_or("--physics needs the code of a physics (dfs, sms, ems), or several separated by commas")?;
                physics = Some(v.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect());
            }
            a if a.starts_with("--") => return Err(format!("unknown option {a}")),
            a => path = Some(a.to_string()),
        }
        i += 1;
    }
    let path = path.ok_or("usage: nyx --batch PROJECT.nyx.json [--physics dfs|sms|ems[,...]]")?;
    let mut engine = Engine::new();
    engine.call("project_open", &json!({ "path": path })).map_err(|e| e.to_string())?;
    let mesh = engine.mesh_now().map_err(|e| e.to_string())?;
    println!("{} cells, {} points", mesh["mesh"]["cells"], engine.model().map(|m| m.mesh().n_points()).unwrap_or(0));
    engine.solve_now(physics).map_err(|e| e.to_string())?;
    let mut ids: Vec<String> = engine.project().physics.iter().map(|s| s.physics.clone()).collect();
    ids.retain(|id| engine.result(id).is_some());
    for id in ids {
        let (result, _) = engine.result(&id).expect("listed");
        println!("{id}: solved in {:.2} s, {} dofs, {} solution(s)", result.seconds, result.dofs, result.solutions.len());
        for (n, v) in &result.scalars {
            match v {
                Scalar::Num(x) => println!("  {n} = {x}"),
                Scalar::Text(t) => println!("  {n} = {t}"),
            }
        }
        if !result.log.is_empty() {
            println!("what the solver did:");
            for l in &result.log {
                println!("  {l}");
            }
        }
    }
    Ok(())
}
