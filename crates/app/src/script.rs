// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Scripted runs for testing and screenshots: `nyx --script "example 1; compute; shot out.png; quit"`.
//!
//! Commands (separated by `;`): `example N`, `new es|heat|solid|emw [3d]`, `open PATH`, `physics es|heat|solid|emw|dfs|sms|ems`,
//! `tab geometry|materials|physics|study|results|analysis`, `compute`, `field NAME`,
//! `sol N`, `cancel_after SECONDS` (cancels the running solve then), `rotate DYAW DPITCH`, `pan DX DY`, `zoom F`, `view mesh|result`, `pick N`, `study eigen|static|transient|stationary`,
//! `save FILE.nyx.json`, `shot FILE`, `wait`, `quit`.

use crate::app::App;
use crate::state::*;
use eframe::egui;
use serde_json::json;
use std::sync::atomic::Ordering;
use std::time::Instant;

impl App {
    pub fn script_step(&mut self, ctx: &egui::Context) {
        if self.script.is_empty() {
            return;
        }
        ctx.request_repaint();
        if self.cancel_at.is_some_and(|t| Instant::now() >= t) {
            self.cancel_at = None;
            self.cancel_job();
        }
        if self.feedback.is_busy() && !self.script_free {
            self.feedback_tick(ctx);
            ctx.request_repaint();
            return;
        }
        if ((self.busy.is_some() || !self.deferred.is_empty() || self.sweep.as_ref().is_some_and(|a| !a.complete()) || self.feedback.is_busy()) && !self.script_free) || self.shot_pending.load(Ordering::SeqCst) {
            return;
        }
        if self.settle > 0 {
            self.settle -= 1;
            return;
        }
        let Some(cmd) = self.script.pop_front() else { return };
        let cmd = cmd.trim().to_string();
        let (name, arg) = cmd.split_once(' ').map(|(a, b)| (a.to_string(), b.trim().to_string())).unwrap_or((cmd.clone(), String::new()));
        self.settle = 4;
        match name.as_str() {
            "example" => {
                if let Ok(i) = arg.parse::<usize>() {
                    if self.shell.active() {
                        if let Some(e) = crate::examples::list().into_iter().nth(i) {
                            let t = e.title.split('(').next().unwrap_or("Example").trim().to_string();
                            self.new_draft(&t, crate::shell::Start::Example(e.path));
                        }
                    } else {
                        self.load_example(i);
                    }
                }
            }
            "new" => {
                let words: Vec<&str> = arg.split_whitespace().collect();
                let physics = short_code(words.first().copied().unwrap_or(""));
                let three_d = words.contains(&"3d");
                let what = self.physics_name(physics);
                if self.shell.active() {
                    let title = format!("{} ({})", what, if three_d { "3D" } else { "2D" });
                    self.new_draft(&title, crate::shell::Start::Blank { physics: physics.into(), three_d });
                } else {
                    self.new_project(physics, three_d, &what);
                }
            }
            "open" => {
                self.shell.sheet = crate::shell::Sheet::None;
                self.open_path(arg.into());
            }
            "physics" => {
                let code = short_code(&arg);
                if self.shell.active() {
                    let _ = self.ensure_active_workflow();
                }
                if self.engine.project().spec(code).is_none() {
                    let _ = self.cmd("add_physics", json!({ "physics": code }));
                }
                self.cmd("set_active", json!({ "physics": code }));
                self.rules_dirty = true;
                self.picking = None;
            }
            "tab" => {
                let target_tab = match arg.as_str() {
                    "geometry" => Tab::Geometry,
                    "materials" => Tab::Materials,
                    "physics" => Tab::Physics,
                    "study" => Tab::Study,
                    "results" => Tab::Results,
                    "analysis" => Tab::Analysis,
                    _ => Tab::Geometry,
                };
                self.tab = target_tab;
                if self.shell.active() {
                    self.shell.sheet = crate::shell::Sheet::None;
                    let _ = self.ensure_active_workflow();
                    let shown = if target_tab == Tab::Analysis { Tab::Results } else { target_tab };
                    if let Some(w) = self.shell.current_mut() {
                        w.tab = shown;
                        w.sub = crate::shell::Sub::Main;
                    }
                }
            }
            "compute" => {
                if self.shell.active() {
                    self.shell.sheet = crate::shell::Sheet::None;
                    let _ = self.ensure_active_workflow();
                    self.shell_compute();
                } else {
                    self.compute();
                }
            }
            // `edit h xFACTOR | VALUE | reset`: the element size of the Geometry form, as typed (not applied); the script then waits
            // deterministically until the check and trial mesh finish (`feedback.is_busy()`), plus two settle frames for repaint.
            "edit" => {
                let words: Vec<&str> = arg.split_whitespace().collect();
                if let ["h", v] = words.as_slice() {
                    let now = self.engine.project().mesh.h;
                    self.mesh_edit.h = if *v == "reset" { now } else if let Some(f) = v.strip_prefix('x').and_then(|f| f.parse::<f64>().ok()) { self.mesh_edit.h * f } else { v.parse().unwrap_or(now) };
                    self.geometry_dirty = self.mesh_edit != self.engine.project().mesh || self.geom_edit != self.engine.project().geometry;
                    self.edits_changed(ctx);
                    self.settle = 2;
                }
            }
            // `check TOLERANCE [FACTOR]`: a mesh-convergence check of the active physics (the script waits for it); the state it leaves is drawn under Results
            "check" => {
                let mut words = arg.split_whitespace().filter_map(|w| w.parse::<f64>().ok());
                if let Some(t) = words.next() {
                    self.conv_tol = t;
                }
                if let Some(f) = words.next() {
                    self.conv_factor = f;
                }
                self.check_convergence();
            }
            // shell layout: `wf example N`, `page dashboard|settings|N`, `section concept|build|analysis`, `sheet new|palette|versions|share|about|none`; looks for a shot (not saved): `theme light|dark|system`, `contrast standard|high|auto`, `motion on|off`
            "wf" => {
                if let Some(e) = arg.strip_prefix("example ").and_then(|n| n.trim().parse::<usize>().ok()).and_then(|i| crate::examples::list().into_iter().nth(i)) {
                    let t = e.title.split('(').next().unwrap_or("Example").trim().to_string();
                    self.new_draft(&t, crate::shell::Start::Example(e.path));
                } else if let Some(p) = arg.strip_prefix("open ") {
                    // `wf open PATH`: a project file as a workflow (for a shot of a model the examples do not have)
                    self.new_draft("Opened model", crate::shell::Start::Example(p.trim().into()));
                }
            }
            "page" => match arg.as_str() {
                "dashboard" => {
                    self.flush();
                    self.refresh_entries();
                    self.shell.sheet = crate::shell::Sheet::None;
                    self.shell.page = crate::shell::Page::Dashboard;
                }
                "settings" => {
                    self.flush();
                    self.shell.sheet = crate::shell::Sheet::None;
                    self.shell.page = crate::shell::Page::Settings;
                }
                n => {
                    if let Some(id) = n.parse::<usize>().ok().and_then(|i| self.shell.tabs.get(i).cloned()) {
                        self.shell.sheet = crate::shell::Sheet::None;
                        self.show_workflow(&id);
                    }
                }
            },
            // `section NAME`: a tab of the left panel of the open workflow (ADR 124). `concept`, `build` and `analysis` are the names of the older window and
            // keep working for the scripts of `tools/ux-shots.sh`: geometry, physics and results. `runs` is the page of the table and the plots.
            "section" => {
                use crate::shell::Sub;
                let (tab, sub) = match arg.as_str() {
                    "materials" => (Tab::Materials, Sub::Main),
                    "physics" | "build" => (Tab::Physics, Sub::Main),
                    "study" | "sweep" => (Tab::Study, Sub::Main),
                    "results" | "analysis" => (Tab::Results, Sub::Main),
                    "runs" => (Tab::Results, Sub::Runs),
                    _ => (Tab::Geometry, Sub::Main),
                };
                if let Some(w) = self.shell.current_mut() {
                    w.tab = tab;
                    w.sub = sub;
                    w.sweep.shown = arg == "sweep";
                }
                self.tab = tab;
            }
            // `theme light|dark`: the look of the window for the shots of `tools/ux-shots.sh` (not saved to the workspace settings)
            "theme" => {
                if let Some(ws) = &mut self.shell.ws {
                    ws.settings.theme = arg.clone();
                }
                self.sync_look(ctx);
            }
            // `contrast standard|high|auto`, `motion on|off` (reduce motion), settings for a shot, not saved
            "contrast" => {
                if let Some(ws) = &mut self.shell.ws {
                    ws.settings.contrast = arg.clone();
                }
                self.sync_look(ctx);
            }
            // `canvas dark|theme`, `colour material|region`, `isolate N|none` (the legend row N, 1-based): settings for a shot, not saved
            "canvas" => {
                if let Some(ws) = &mut self.shell.ws {
                    ws.settings.canvas = arg.clone();
                }
            }
            "colour" => {
                if let Some(ws) = &mut self.shell.ws {
                    ws.settings.region_colouring = arg.clone();
                }
            }
            // `wait off|on`: with `off` the script carries on while a job runs, so that a running job and the queue can be shot (`on`, the default,
            // waits for the job to end before the next command); `queue pause|resume|open|close`: the engine's queue and its popover
            "wait" => self.script_free = arg == "off",
            // `cvd protan|deutan|tritan|none`: the shots that follow show the window as a person with that colour-vision difference sees it (cvd.rs);
            // `scale F`: the interface scale (1.5 and 2.0 are the 150 and 200 percent of Windows), not saved
            "cvd" => self.shot_cvd = crate::cvd::Kind::parse(&arg),
            "scale" => {
                if let Ok(f) = arg.parse::<f32>() {
                    let clamped = f.clamp(0.5, 3.0);
                    self.shell.script_scale = Some(clamped);
                    self.apply_zoom_factor(ctx, clamped);
                }
            }
            // `tick N M ...`: ticks those runs (1-based, in the order of the table's data) in the Analysis table, for a shot of the ringed runs
            "tick" => {
                if let Some(w) = self.shell.current_mut() {
                    for n in arg.split_whitespace().filter_map(|x| x.parse::<usize>().ok()) {
                        if let Some(r) = n.checked_sub(1).and_then(|i| w.runs.get(i)) {
                            w.selected.insert(r.id.clone());
                        }
                    }
                }
            }
            // `sweep show` opens the Sweep section of the Build page; `sweep run PARAM FROM TO N` fills it in and starts the sweep (for shots)
            "sweep" => {
                let words: Vec<&str> = arg.split_whitespace().collect();
                let Some(w) = self.shell.current_mut() else { return };
                w.sweep.shown = true;
                if let ["run", param, from, to, n] = words.as_slice() {
                    w.sweep.param = param.to_string();
                    w.sweep.from = from.to_string();
                    w.sweep.to = to.to_string();
                    w.sweep.n = n.parse().unwrap_or(5);
                    let (id, form) = (w.id.clone(), w.sweep.clone());
                    self.start_sweep(&id, &form);
                }
            }
            "queue" => match arg.as_str() {
                "pause" => self.engine.jobs().pause(),
                "resume" => self.engine.jobs().resume(),
                "open" => self.show_queue = true,
                _ => self.show_queue = false,
            },
            // `legend all`: opens the list of every legend row (what a click on "N more" does)
            "legend" => self.vp.legend_all = arg == "all",
            "isolate" => self.vp.isolate = arg.parse::<usize>().ok().and_then(|n| n.checked_sub(1)),
            "motion" => {
                if let Some(ws) = &mut self.shell.ws {
                    ws.settings.reduce_motion = arg == "on";
                }
                self.sync_look(ctx);
            }
            // `done N`: marks the workflow of tab N done, so that it moves to Earlier on the Dashboard (for the shot of that table)
            "done" => {
                if let Some(id) = arg.parse::<usize>().ok().and_then(|i| self.shell.tabs.get(i).cloned()) {
                    self.edit_workflow(&id, |d| d.set_status(nyx_engine::facade::store::Status::Done));
                    self.refresh_entries();
                }
            }
            "sheet" => {
                use crate::shell::Sheet;
                self.shell.sheet = match arg.as_str() {
                    "new" => Sheet::New,
                    "palette" => {
                        self.shell.palette_q.clear();
                        Sheet::Palette
                    }
                    "versions" => Sheet::Versions,
                    "share" => Sheet::Share,
                    "about" => Sheet::About,
                    // `sheet script` opens the script dock (ADR 124), the way the key does
                    "script" => {
                        self.nys.open = false;
                        self.open_nys_script();
                        Sheet::None
                    }
                    _ => Sheet::None,
                };
            }
            "cancel_after" => self.cancel_at = arg.parse::<f64>().ok().map(|s| Instant::now() + std::time::Duration::from_secs_f64(s)),
            "field" => {
                if let Some(i) = self.current {
                    let rv = &mut self.results[i];
                    if let Some(p) = rv.solution().fields.iter().position(|f| f.name == arg) {
                        rv.field = p;
                        rv.values_cache = None;
                    }
                }
            }
            "sol" => {
                if let (Some(i), Ok(n)) = (self.current, arg.parse::<usize>()) {
                    let rv = &mut self.results[i];
                    rv.sol = n.min(rv.sols.len() - 1);
                }
            }
            "save" => {
                if self.push_edits() && self.cmd("project_save", json!({ "path": arg })).is_some() {
                    self.status = format!("Saved {arg}");
                }
            }
            "view" => {
                self.show_results = arg == "result" || arg == "field";
                // `view field`: the Results tab (where the validated state of the result is laid out) with the field on show
                if arg == "field" {
                    if let Some(w) = self.shell.current_mut() {
                        w.tab = Tab::Results;
                        w.sub = crate::shell::Sub::Main;
                    }
                    self.tab = Tab::Results;
                }
            }
            // camera moves, for timing and screenshots: `rotate DYAW DPITCH` (radians), `pan DX DY` (pixels), `zoom FACTOR`
            "rotate" => {
                let v: Vec<f32> = arg.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if v.len() == 2 {
                    self.vp.cam.yaw += v[0];
                    self.vp.cam.pitch = (self.vp.cam.pitch + v[1]).clamp(-1.55, 1.55);
                }
            }
            "pan" => {
                let v: Vec<f32> = arg.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                if v.len() == 2 {
                    self.vp.cam.pan += egui::Vec2::new(v[0], v[1]);
                }
            }
            "zoom" => {
                if let Ok(f) = arg.parse::<f32>() {
                    self.vp.cam.zoom = (self.vp.cam.zoom * f).clamp(0.02, 1e5);
                }
            }
            "pick" => self.picking = arg.parse().ok(),
            "study" => {
                // eigen | static | transient | stationary: the study of that name in the active physics, with its defaults
                if let Some(active) = self.active() {
                    let kinds: Vec<String> = self.engine.registry().get(&active).map(|p| p.studies().iter().map(|s| s.kind.clone()).collect()).unwrap_or_default();
                    let want = match arg.as_str() {
                        "eigen" => "eigenfrequency",
                        "transient" => "heat_transient",
                        "stationary" => "heat_stationary",
                        other => other,
                    };
                    if let Some(k) = kinds.iter().find(|k| k.as_str() == want || k.starts_with(want)).cloned() {
                        self.cmd("set_study", json!({ "physics": active, "kind": k, "params": if k == "eigenfrequency" { json!({ "n_modes": 6, "f_search": 0.0 }) } else if k == "heat_transient" { json!({ "t_end": 20000.0, "n_out": 20, "t0": "293.15" }) } else { json!(null) } }));
                    }
                }
            }
            "shot" => {
                let path = arg.clone();
                let flag = self.shot_pending.clone();
                let cvd = self.shot_cvd;
                flag.store(true, Ordering::SeqCst);
                let t_shot = Instant::now();
                ctx.request_screenshot(move |img| {
                    if std::env::var_os("NYX_FRAME_LOG").is_some() {
                        eprintln!("SHOT delivered after {:.1} ms", t_shot.elapsed().as_secs_f64() * 1e3);
                    }
                    let [w, h] = img.size;
                    let mut raw = Vec::with_capacity(w * h * 4);
                    for p in &img.pixels {
                        raw.extend_from_slice(&[p.r(), p.g(), p.b(), p.a()]);
                    }
                    if let Some(k) = cvd {
                        crate::cvd::simulate_rgba(&mut raw, k);
                    }
                    if let Some(buf) = image::RgbaImage::from_raw(w as u32, h as u32, raw) {
                        let _ = buf.save(&path);
                    }
                    flag.store(false, Ordering::SeqCst);
                });
            }
            "quit" => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
            _ => {}
        }
    }
}

/// The code of a physics from the short names the script has always used (`es`, `heat`, `solid`, `emw`) or from its own code.
fn short_code(s: &str) -> &'static str {
    match s {
        "es" | "heat" | "dfs" => "dfs",
        "solid" | "sms" => "sms",
        _ => "ems",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feedback::SizeView;
    use crate::shell::{Page, Sheet};
    use std::time::{Duration, Instant};

    #[test]
    fn script_step_drives_adr124_shell_commands() {
        let root = std::env::temp_dir().join(format!("nyx-script-shell-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let mut a = App::with_context(egui::Context::default());
        a.feedback.slow_limit = Some(f64::INFINITY);
        a.attach_workspace(root.clone());
        assert!(a.shell.active());
        assert_eq!(a.shell.page, Page::Dashboard);

        let step = |app: &mut App, cmds: &[&str]| {
            app.script = cmds.iter().map(|s| (*s).to_string()).collect();
            let t0 = Instant::now();
            while !app.script.is_empty() || app.busy.is_some() || !app.deferred.is_empty() {
                app.settle = 0;
                app.poll();
                app.run_queued();
                let ctx = app.ctx.clone();
                app.script_step(&ctx);
                assert!(t0.elapsed() < Duration::from_secs(30), "script timed out: {:?}", app.script);
                std::thread::sleep(Duration::from_millis(5));
            }
        };

        step(&mut a, &["new emw", "tab physics"]);
        let wid_new = match &a.shell.page {
            Page::Workflow(id) => id.clone(),
            other => panic!("expected workflow page after `new emw`, got {other:?}"),
        };
        assert_eq!(a.active().as_deref(), Some("ems"));
        assert_eq!(a.shell.open[&wid_new].tab, Tab::Physics);

        step(&mut a, &["physics dfs"]);
        assert_eq!(a.active().as_deref(), Some("dfs"));

        step(&mut a, &["page settings"]);
        assert_eq!(a.shell.page, Page::Settings);
        assert_eq!(a.shell.sheet, Sheet::None);

        step(&mut a, &["sheet palette"]);
        assert_eq!(a.shell.sheet, Sheet::Palette);
        step(&mut a, &["sheet about"]);
        assert_eq!(a.shell.sheet, Sheet::About);
        step(&mut a, &["sheet new"]);
        assert_eq!(a.shell.sheet, Sheet::New);
        step(&mut a, &["sheet none"]);
        assert_eq!(a.shell.sheet, Sheet::None);

        step(&mut a, &["example 0", "tab geometry"]);
        let wid_ex = match &a.shell.page {
            Page::Workflow(id) => id.clone(),
            other => panic!("expected example workflow page, got {other:?}"),
        };
        assert_ne!(wid_ex, wid_new);
        assert_eq!(a.shell.tabs.len(), 2);
        assert_eq!(a.shell.open[&wid_ex].tab, Tab::Geometry);

        // `compute` through the shell starts/records a workflow solve and returns idle cleanly
        step(&mut a, &["compute"]);
        assert!(!a.results.is_empty(), "expected solve result in shell");

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn script_edit_waits_deterministically_for_trial_mesh_feedback() {
        let root = std::env::temp_dir().join(format!("nyx-script-edit-wait-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let mut a = App::with_context(egui::Context::default());
        a.feedback.slow_limit = Some(f64::INFINITY);
        a.attach_workspace(root.clone());

        // Open example 0 in the shell and finish initial rebuild_mesh
        let ctx = a.ctx.clone();
        a.script = ["example 0"].into_iter().map(str::to_string).collect();
        let t0 = Instant::now();
        while !a.script.is_empty() || a.busy.is_some() || !a.deferred.is_empty() {
            a.settle = 0;
            a.poll();
            a.run_queued();
            a.script_step(&ctx);
            assert!(t0.elapsed() < Duration::from_secs(20));
            std::thread::sleep(Duration::from_millis(5));
        }

        // Queue `edit h x0.6` followed by `tab results` with zero extra settle delay
        a.script = ["edit h x0.6", "tab results"].into_iter().map(str::to_string).collect();
        a.settle = 0;
        a.script_step(&ctx);
        assert!(a.feedback.is_busy(), "feedback timer/trial mesh must be busy right after edit");

        // While inside the REST window (first 50 ms), `tab results` MUST NOT advance even after many fast ticks
        for _ in 0..20 {
            a.poll();
            a.script_step(&ctx);
            assert_eq!(a.script.len(), 1, "`tab results` must wait while feedback is busy");
            assert_eq!(a.shell.current_mut().unwrap().tab, Tab::Geometry);
            std::thread::sleep(Duration::from_millis(2));
        }

        // Once the REST delay and background trial mesh complete, `tab results` executes with Ready size line
        let wait_t = Instant::now();
        while !a.script.is_empty() {
            a.poll();
            a.script_step(&ctx);
            assert!(wait_t.elapsed() < Duration::from_secs(10), "edit wait timed out");
            std::thread::sleep(Duration::from_millis(5));
        }
        assert!(
            matches!(a.feedback.size, SizeView::Ready(_)),
            "expected Ready trial-mesh feedback before advancing, got {:?}",
            a.feedback.size
        );
        assert_eq!(a.shell.current_mut().unwrap().tab, Tab::Results);

        let _ = std::fs::remove_dir_all(&root);
    }
}
