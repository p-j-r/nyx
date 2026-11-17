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
        if (self.busy.is_some() && !self.script_free) || self.shot_pending.load(Ordering::SeqCst) {
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
                    self.load_example(i);
                }
            }
            "new" => {
                let words: Vec<&str> = arg.split_whitespace().collect();
                let physics = short_code(words.first().copied().unwrap_or(""));
                let three_d = words.contains(&"3d");
                let what = self.physics_name(physics);
                self.new_project(physics, three_d, &what);
            }
            "open" => self.open_path(arg.into()),
            "physics" => {
                self.cmd("set_active", json!({ "physics": short_code(&arg) }));
                self.rules_dirty = true;
            }
            "tab" => {
                self.tab = match arg.as_str() {
                    "geometry" => Tab::Geometry,
                    "materials" => Tab::Materials,
                    "physics" => Tab::Physics,
                    "study" => Tab::Study,
                    "results" => Tab::Results,
                    "analysis" => Tab::Analysis,
                    _ => Tab::Geometry,
                }
            }
            "compute" => self.shell_compute(),
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
                    self.shell.page = crate::shell::Page::Dashboard;
                }
                "settings" => self.shell.page = crate::shell::Page::Settings,
                n => {
                    if let Some(id) = n.parse::<usize>().ok().and_then(|i| self.shell.tabs.get(i).cloned()) {
                        self.show_workflow(&id);
                    }
                }
            },
            "section" => {
                if let Some(w) = self.shell.current_mut() {
                    w.sub = crate::shell::Sub::Main;
                    w.part = match arg.as_str() {
                        "build" => crate::shell::Part::Build,
                        "analysis" => crate::shell::Part::Analysis,
                        _ => crate::shell::Part::Concept,
                    };
                }
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
            "queue" => match arg.as_str() {
                "pause" => self.engine.jobs().pause(),
                "resume" => self.engine.jobs().resume(),
                "open" => self.show_queue = true,
                _ => self.show_queue = false,
            },
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
                    "palette" => Sheet::Palette,
                    "versions" => Sheet::Versions,
                    "share" => Sheet::Share,
                    "about" => Sheet::About,
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
            "view" => self.show_results = arg == "result",
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
