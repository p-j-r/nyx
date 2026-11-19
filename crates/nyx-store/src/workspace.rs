//! The workspace: the folder behind the dashboard, its settings (`workspace.json`), the list of workflows (manifests only) and the
//! dashboard retention rules.
//!
//! Nothing here deletes a workflow. "Disappears from the dashboard" means moved to *Earlier*, where search still finds it. The only files
//! housekeeping removes are the heavy artifacts of old failed runs, and of unstarred runs beyond the newest N when the person chose to
//! keep only N (the default is to keep everything, as in the Settings screen).

use crate::model::*;
use crate::runs::RunStatus;
use crate::util::*;
use crate::flow::WorkflowDir;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    /// hide a workflow from the dashboard after this many days without an edit, run or note
    pub hide_after_days: u32,
    /// cards on the dashboard at most (pinned and recent)
    pub dashboard_cap: usize,
    /// keep the heavy files of the newest N unstarred runs per experiment; `None` keeps all of them
    pub keep_heavy: Option<usize>,
    /// remove heavy files of failed and cancelled runs after this many days; `None` never
    pub failed_heavy_days: Option<u32>,
    pub share_default: ShareLevel,
    pub strip_paths: bool,
    /// "ask first" or "never": shared scripts are never run silently
    pub run_shared_scripts: String,
    pub run_at_once: usize,
    pub length_unit: String,
    pub frequency_unit: String,
    /// "system", "light" or "dark"
    pub theme: String,
    pub assistant: bool,
    /// what the assistant may do when it is switched on
    pub assistant_level: String,
    /// the viewport and plot background: "dark" (the default, also in a light window) or "theme" (follows the look)
    pub canvas: String,
    /// "auto" (the system High Contrast flag where the platform has one, else standard), "standard" or "high"
    pub contrast: String,
    /// no animation: no sliding progress block, no easing
    pub reduce_motion: bool,
    /// what a region's colour follows: "material" (the default) or "region" (its number)
    pub region_colouring: String,
    /// the window zoom set with Ctrl and +/-, 1.0 = 100 %
    pub zoom: f32,
    /// workflow ids open as tabs, in order
    pub open_tabs: Vec<String>,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            hide_after_days: 3,
            dashboard_cap: 12,
            keep_heavy: None,
            failed_heavy_days: None,
            share_default: ShareLevel::Recipe,
            strip_paths: true,
            run_shared_scripts: "ask first".into(),
            run_at_once: 1,
            length_unit: "nm".into(),
            frequency_unit: "THz".into(),
            theme: "system".into(),
            assistant: false,
            assistant_level: "Ask before changing a setup".into(),
            open_tabs: vec![],
            canvas: "dark".into(),
            contrast: "auto".into(),
            reduce_motion: false,
            region_colouring: "material".into(),
            zoom: 1.0,
        }
    }
}

/// One workflow on disk as the dashboard sees it.
#[derive(Clone, Debug)]
pub struct Entry {
    pub dir: PathBuf,
    pub wf: Workflow,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Sections {
    /// indices into the entries given to `sections`
    pub running: Vec<usize>,
    pub recent: Vec<usize>,
    pub earlier: Vec<usize>,
}

pub struct Workspace {
    pub root: PathBuf,
    pub settings: Settings,
}

/// `Documents\nyx` of the user (the current folder when there is no home).
pub fn default_root() -> PathBuf {
    std::env::var("USERPROFILE").or_else(|_| std::env::var("HOME")).map(|h| PathBuf::from(h).join("Documents").join("nyx")).unwrap_or_else(|_| PathBuf::from("nyx"))
}

impl Workspace {
    /// Opens (and creates) the workspace folder; a missing or unreadable `workspace.json` gives the default settings.
    pub fn open(root: &Path) -> Result<Workspace, String> {
        std::fs::create_dir_all(root).map_err(|e| io_err(root, e))?;
        let settings = std::fs::read_to_string(root.join("workspace.json")).ok().and_then(|t| serde_json::from_str(&t).ok()).unwrap_or_default();
        Ok(Workspace { root: root.to_path_buf(), settings })
    }

    pub fn save_settings(&self) -> Result<(), String> {
        let p = self.root.join("workspace.json");
        write_atomic(&p, serde_json::to_string_pretty(&self.settings).map_err(|e| e.to_string())?.as_bytes()).map_err(|e| io_err(&p, e))
    }

    /// Every workflow folder with a readable manifest, newest id first. Only `workflow.json` is read.
    pub fn list(&self) -> Vec<Entry> {
        let Ok(rd) = std::fs::read_dir(&self.root) else { return vec![] };
        let mut v: Vec<Entry> = rd
            .filter_map(|e| e.ok())
            .filter(|e| e.path().join("workflow.json").is_file())
            .filter_map(|e| {
                let text = std::fs::read_to_string(e.path().join("workflow.json")).ok()?;
                let wf: Workflow = serde_json::from_str(&text).ok()?;
                Some(Entry { dir: e.path(), wf })
            })
            .collect();
        v.sort_by(|a, b| b.wf.id.cmp(&a.wf.id));
        v
    }

    pub fn find(&self, id: &str) -> Option<Entry> {
        self.list().into_iter().find(|e| e.wf.id == id)
    }

    /// Dashboard sections for `entries` at time `now`:
    /// * *running*: the ids in `running`, always shown first;
    /// * *recent*: pinned (always, up to the cap) then not-done workflows touched within `hide_after_days`, newest first, up to
    ///   `dashboard_cap` in all; the oldest unpinned fall to *earlier*;
    /// * *earlier*: everything else, including workflows marked done. Viewing never counts as touching.
    pub fn sections(&self, entries: &[Entry], running: &[String], now: i64) -> Sections {
        let s = &self.settings;
        let window = s.hide_after_days as i64 * 86_400;
        let mut sec = Sections::default();
        let mut cand: Vec<usize> = vec![];
        for (i, e) in entries.iter().enumerate() {
            if running.contains(&e.wf.id) {
                sec.running.push(i);
            } else if e.wf.pinned || (e.wf.status == Status::Active && now - e.wf.touched_secs() <= window) {
                cand.push(i);
            } else {
                sec.earlier.push(i);
            }
        }
        // pinned first, then by recency
        cand.sort_by_key(|&i| (!entries[i].wf.pinned, std::cmp::Reverse(entries[i].wf.touched_secs())));
        let room = s.dashboard_cap.saturating_sub(sec.running.len());
        let (keep, fall) = cand.split_at(room.min(cand.len()));
        sec.recent = keep.to_vec();
        sec.earlier.extend_from_slice(fall);
        sec.earlier.sort_by_key(|&i| std::cmp::Reverse(entries[i].wf.touched_secs()));
        sec
    }

    /// Removes heavy artifacts (everything under `runs/<id>/` except the log) according to the settings. Returns the bytes freed.
    /// Starred runs and every run record are kept.
    pub fn housekeeping(&self, w: &mut WorkflowDir, now: i64) -> Result<u64, String> {
        let mut freed = 0;
        let ids: Vec<String> = w.wf.experiments.iter().map(|e| e.id.clone()).collect();
        for eid in ids {
            let (rs, _) = w.runs(&eid);
            let ok_unstarred: Vec<&crate::runs::RunRecord> = rs.iter().filter(|r| r.status == RunStatus::Ok && !r.star).collect();
            let mut doomed: Vec<&str> = vec![];
            if let Some(n) = self.settings.keep_heavy {
                let cut = ok_unstarred.len().saturating_sub(n);
                doomed.extend(ok_unstarred[..cut].iter().map(|r| r.id.as_str()));
            }
            if let Some(days) = self.settings.failed_heavy_days {
                doomed.extend(rs.iter().filter(|r| r.status != RunStatus::Ok && !r.star && parse_iso(&r.t).is_some_and(|t| now - t > days as i64 * 86_400)).map(|r| r.id.as_str()));
            }
            for rid in doomed {
                let dir = w.run_dir(&eid, rid);
                let Ok(rd) = std::fs::read_dir(&dir) else { continue };
                for f in rd.filter_map(|f| f.ok()) {
                    let p = f.path();
                    if p.file_name().is_some_and(|n| n == "log.txt") || !p.is_file() {
                        continue;
                    }
                    freed += f.metadata().map(|m| m.len()).unwrap_or(0);
                    std::fs::remove_file(&p).map_err(|e| io_err(&p, e))?;
                }
            }
        }
        Ok(freed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow::tests::{spec_json, tmp};
    use crate::runs::RunRecord;

    fn entry(title: &str, age_days: i64, now: i64) -> Entry {
        let mut wf = Workflow::new(title, "gui");
        wf.touched = iso(now - age_days * 86_400);
        Entry { dir: PathBuf::from(title), wf }
    }

    #[test]
    fn dashboard_rules_hide_after_three_days_cap_pin_and_never_delete() {
        let now = 1_800_000_000;
        let ws = Workspace { root: PathBuf::new(), settings: Settings { dashboard_cap: 3, ..Default::default() } };
        let mut es = vec![entry("fresh", 0, now), entry("two days", 2, now), entry("old", 10, now), entry("done", 0, now), entry("pinned old", 40, now), entry("third", 1, now), entry("running old", 30, now)];
        es[3].wf.status = Status::Done;
        es[4].wf.pinned = true;
        let s = ws.sections(&es, &[es[6].wf.id.clone()], now);
        let names = |v: &[usize]| v.iter().map(|&i| es[i].wf.title.as_str()).collect::<Vec<_>>();
        assert_eq!(names(&s.running), ["running old"]);
        // cap 3 with one running leaves room for 2: the pinned one first, then the freshest
        assert_eq!(names(&s.recent), ["pinned old", "fresh"]);
        // everything else is still listed under Earlier
        let mut earlier = names(&s.earlier);
        earlier.sort();
        assert_eq!(earlier, ["done", "old", "third", "two days"]);
        assert_eq!(s.running.len() + s.recent.len() + s.earlier.len(), es.len(), "nothing is dropped");
        // edge: exactly three days is still recent, a second more is not
        let e3 = entry("edge", 3, now);
        assert_eq!(ws.sections(std::slice::from_ref(&e3), &[], now).recent.len(), 1);
        assert_eq!(ws.sections(&[e3], &[], now + 1).earlier.len(), 1);
    }

    #[test]
    fn the_workspace_lists_manifests_and_keeps_settings() {
        let root = tmp("ws");
        let mut ws = Workspace::open(&root).unwrap();
        assert_eq!(ws.settings.hide_after_days, 3);
        let a = WorkflowDir::create(&root, "First", "gui", &spec_json(5e-7, 4)).unwrap();
        std::thread::sleep(std::time::Duration::from_millis(1100));
        let b = WorkflowDir::create(&root, "Second", "gui", &spec_json(6e-7, 4)).unwrap();
        std::fs::create_dir_all(root.join("not-a-workflow")).unwrap();
        let l = ws.list();
        assert_eq!(l.len(), 2);
        assert_eq!(l[0].wf.id, b.wf.id, "newest first");
        assert!(ws.find(&a.wf.id).is_some());
        ws.settings.hide_after_days = 7;
        ws.settings.open_tabs = vec![a.wf.id.clone()];
        ws.save_settings().unwrap();
        assert_eq!(Workspace::open(&root).unwrap().settings, ws.settings);
        std::fs::write(root.join("workspace.json"), "{ broken").unwrap();
        assert_eq!(Workspace::open(&root).unwrap().settings, Settings::default(), "a damaged settings file does not stop nyx");
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn housekeeping_frees_old_heavy_files_but_keeps_logs_stars_and_records() {
        let root = tmp("hk");
        let mut ws = Workspace::open(&root).unwrap();
        ws.settings.keep_heavy = Some(1);
        let mut d = WorkflowDir::create(&root, "hk", "gui", &spec_json(5e-7, 4)).unwrap();
        let spec = d.read_spec("e1").unwrap();
        for i in 0..3 {
            let r = d.record_run("e1", RunRecord { sec: 0.1, ..Default::default() }, &spec).unwrap();
            let dir = d.run_dir("e1", &r.id);
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("results.h5"), vec![0u8; 100]).unwrap();
            std::fs::write(dir.join("log.txt"), "log").unwrap();
            if i == 0 {
                d.patch_run("e1", &r.id, crate::runs::RunPatch { star: Some(true), ..Default::default() }).unwrap();
            }
        }
        let freed = ws.housekeeping(&mut d, now_secs()).unwrap();
        // 3 runs: run 1 starred (kept), runs 2 and 3 unstarred, newest 1 kept -> run 2's heavy file goes
        assert_eq!(freed, 100);
        let (rs, _) = d.runs("e1");
        assert_eq!(rs.len(), 3);
        assert!(!d.run_dir("e1", &rs[1].id).join("results.h5").exists());
        assert!(d.run_dir("e1", &rs[1].id).join("log.txt").exists());
        assert!(d.run_dir("e1", &rs[0].id).join("results.h5").exists());
        assert!(d.run_dir("e1", &rs[2].id).join("results.h5").exists());
        let _ = std::fs::remove_dir_all(&root);
    }
}
