// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Runs an optional user Python script: `python script.py project.json [result.csv]`. Python is never needed by the application itself.
//! The script runs on its own thread; its output comes back to the window as a message.

use crate::app::Msg;
use eframe::egui;
use std::path::{Path, PathBuf};
use std::sync::mpsc::Sender;

pub fn spawn(script: PathBuf, project_json: String, csv: Option<String>, tx: Sender<(u64, Msg)>, ctx: egui::Context) {
    std::thread::spawn(move || {
        let run = || -> Result<String, String> {
            let dir = std::env::temp_dir().join(format!("nyx-script-{}", std::process::id()));
            std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            let pj = dir.join("project.json");
            std::fs::write(&pj, project_json).map_err(|e| e.to_string())?;
            let mut args: Vec<std::ffi::OsString> = vec![script.clone().into(), pj.clone().into()];
            if let Some(c) = csv {
                let cp = dir.join("result.csv");
                std::fs::write(&cp, c).map_err(|e| e.to_string())?;
                args.push(cp.into());
            }
            let mut candidates: Vec<(String, Vec<String>)> = vec![];
            if let Ok(p) = std::env::var("NYX_PYTHON") {
                candidates.push((p, vec![]));
            }
            candidates.push(("python".into(), vec![]));
            candidates.push(("py".into(), vec!["-3".into()]));
            candidates.push(("python3".into(), vec![]));
            for (exe, pre) in candidates {
                let out = std::process::Command::new(&exe).args(&pre).args(&args).current_dir(script.parent().unwrap_or(Path::new("."))).output();
                match out {
                    Ok(o) => {
                        let mut text = String::from_utf8_lossy(&o.stdout).into_owned();
                        let err = String::from_utf8_lossy(&o.stderr);
                        if !err.trim().is_empty() {
                            text += "\n--- stderr ---\n";
                            text += &err;
                        }
                        if !o.status.success() {
                            text += &format!("\n(exit status {})", o.status);
                        }
                        return Ok(text);
                    }
                    Err(_) => continue,
                }
            }
            Err("Python was not found. It is optional: it is only used to run your own scripts (set NYX_PYTHON to the interpreter to use).".into())
        };
        let _ = tx.send((0, Msg::Script(run())));
        ctx.request_repaint();
    });
}
