//! Scripted runs for testing and screenshots: `nyx --script "example 1; compute; shot out.png; quit"`.
//!
//! Commands (separated by `;`): `example N`, `new es|heat|solid|emw [3d]`, `open PATH`, `physics es|heat|solid|emw`,
//! `tab geometry|materials|physics|study|results|analysis`, `compute`, `field NAME`,
//! `sol N`, `cancel_after SECONDS` (cancels the running solve then), `rotate DYAW DPITCH`, `pan DX DY`, `zoom F`, `view mesh|result`, `pick N`, `study eigen|static|transient|stationary`,
//! `output serial|parallel FILE`, `save FILE.xml|json`, `shot FILE`, `wait`, `quit`.

use crate::app::App;
use crate::state::*;
use eframe::egui;
use femcore::problem::*;
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
        if self.busy.is_some() || self.shot_pending.load(Ordering::SeqCst) {
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
                let physics = match words.first().copied().unwrap_or("") {
                    "es" => Physics::Electrostatics,
                    "heat" => Physics::Heat,
                    "solid" => Physics::Solid,
                    _ => Physics::Emw,
                };
                let three_d = words.contains(&"3d");
                self.new_project(crate::examples::blank(physics, three_d), physics.name());
            }
            "open" => self.open_path(arg.into()),
            "physics" => {
                self.proj.active = match arg.as_str() {
                    "es" => Physics::Electrostatics,
                    "heat" => Physics::Heat,
                    "solid" => Physics::Solid,
                    _ => Physics::Emw,
                }
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
            // shell layout: `wf example N`, `page dashboard|settings|N`, `section concept|build|analysis`, `sheet new|palette|versions|share`
            "wf" => {
                if let Some(i) = arg.strip_prefix("example ").and_then(|n| n.trim().parse::<usize>().ok()).filter(|i| *i < crate::examples::NAMES.len()) {
                    let t = crate::examples::NAMES[i].split('(').next().unwrap_or("Example").trim().to_string();
                    self.new_draft(&t, crate::examples::by_index(i));
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
            "sheet" => {
                use crate::shell::Sheet;
                self.shell.sheet = match arg.as_str() {
                    "new" => Sheet::New,
                    "palette" => Sheet::Palette,
                    "versions" => Sheet::Versions,
                    "share" => Sheet::Share,
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
                    rv.sol = n.min(rv.result.solutions.len() - 1);
                }
            }
            "output" => {
                let phys = self.proj.active;
                let sp = self.proj.spec_mut(phys);
                sp.output.mode = if arg.starts_with("parallel") { femcore::problem::OutputMode::Hdf5Parallel } else { femcore::problem::OutputMode::Hdf5Serial };
                sp.output.path = arg.trim_start_matches("parallel ").trim_start_matches("serial ").to_string();
            }
            "save" => {
                self.proj.geometry = self.geom_edit.clone();
                let text = if arg.ends_with(".xml") { self.proj.to_xml() } else { serde_json::to_string_pretty(&self.proj).map_err(|e| e.to_string()) };
                match text.and_then(|t| std::fs::write(&arg, t).map_err(|e| e.to_string())) {
                    Ok(()) => self.status = format!("Saved {arg}"),
                    Err(e) => self.set_error(e),
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
                let phys = self.proj.active;
                let s = self.proj.spec_mut(phys);
                s.study = match arg.as_str() {
                    "eigen" => Study::Eigenfrequency { n_modes: 6, f_search: 0.0 },
                    "static" => Study::Static,
                    "stationary" => Study::Stationary,
                    "transient" => Study::Transient { t_end: 20000.0, n_out: 20, t0: "293.15".into() },
                    _ => s.study.clone(),
                };
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
