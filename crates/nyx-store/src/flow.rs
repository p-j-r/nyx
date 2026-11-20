// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! One workflow folder: `workflow.json`, `history.jsonl`, `notes.md`, `specs/` and one folder per experiment.
//!
//! ```text
//! <slug>_<id>/                  the folder never moves when the title changes
//!   workflow.json               the manifest (KB)
//!   history.jsonl               versions of the set-up
//!   notes.md                    free text
//!   specs/<hash>.nyx.json        frozen specs, content addressed (versions and runs point here)
//!   analysis/                   scripts (never run until the person approves them)
//!   <experiment id>/
//!     experiment.nyx.json        the current spec (a `nyx_project`, as the app saves it)
//!     runs.jsonl                append-only run records
//!     runs/<run id>/            heavy and deletable: log.txt, results.h5, ...
//! ```
//!
//! Every write of a whole file is write-then-rename. The manifest's run digest is derived from `runs.jsonl` and rebuilt when it is missing
//! or stale, so a crash never loses a run.

use crate::model::*;
use crate::runs::{self, RunPatch, RunRecord};
use crate::util::*;
use crate::versions::{self, Version, VersionKind};
use std::path::{Path, PathBuf};

/// Automatic versions closer than this on the same experiment (and with no run on the earlier one) are merged.
pub const MERGE_SECS: i64 = 300;

#[derive(Debug)]
pub struct WorkflowDir {
    pub root: PathBuf,
    pub wf: Workflow,
    /// Makes the search card of an experiment from its spec value. The default reads only the title; the engine sets its own
    /// (`create_with`, `open_with`), because the store does not know what a project looks like (R21, R25).
    pub card_of: CardFn,
}

/// A function from the value of a spec to its card.
pub type CardFn = fn(&serde_json::Value) -> Card;

/// The value of a spec text: the project as JSON.
fn parse_spec(text: &str) -> Result<serde_json::Value, String> {
    serde_json::from_str(text).map_err(|e| format!("not a nyx project: {e}"))
}

impl WorkflowDir {
    // ------------------------------------------------------------------------------------------------------ paths

    pub fn manifest_path(&self) -> PathBuf {
        self.root.join("workflow.json")
    }
    pub fn history_path(&self) -> PathBuf {
        self.root.join("history.jsonl")
    }
    pub fn notes_path(&self) -> PathBuf {
        self.root.join("notes.md")
    }
    pub fn spec_path(&self, eid: &str) -> PathBuf {
        self.root.join(eid).join("experiment.nyx.json")
    }
    pub fn runs_path(&self, eid: &str) -> PathBuf {
        self.root.join(eid).join("runs.jsonl")
    }
    pub fn run_dir(&self, eid: &str, rid: &str) -> PathBuf {
        self.root.join(eid).join("runs").join(rid)
    }
    pub fn frozen_path(&self, hash: &str) -> PathBuf {
        self.root.join("specs").join(format!("{}.nyx.json", hash.replace(':', "-")))
    }

    // ------------------------------------------------------------------------------------------------------ create, open, save

    /// Creates the folder `<workspace>/<slug>_<id>` with one experiment holding `spec_json` (the first version).
    pub fn create(workspace: &Path, title: &str, origin: &str, spec_json: &str) -> Result<WorkflowDir, String> {
        Self::create_as(workspace, None, title, origin, spec_json)
    }

    /// Like `create` with a chosen id (a draft that has lived in a tab keeps the id it was shown under).
    pub fn create_as(workspace: &Path, id: Option<&str>, title: &str, origin: &str, spec_json: &str) -> Result<WorkflowDir, String> {
        Self::create_with(workspace, id, title, origin, spec_json, crate::card::plain)
    }

    /// Like `create_as` with the function that makes cards from specs.
    pub fn create_with(workspace: &Path, id: Option<&str>, title: &str, origin: &str, spec_json: &str, card_of: CardFn) -> Result<WorkflowDir, String> {
        let mut wf = Workflow::new(title, origin);
        if let Some(id) = id {
            wf.id = id.to_string();
        }
        let root = workspace.join(format!("{}_{}", slug(title), wf.id));
        std::fs::create_dir_all(&root).map_err(|e| io_err(&root, e))?;
        wf.experiments.clear();
        let mut d = WorkflowDir { root, wf, card_of };
        d.add_experiment("Experiment 1", spec_json, None)?;
        d.write_gitignore()?;
        Ok(d)
    }

    pub fn open(root: &Path) -> Result<WorkflowDir, String> {
        Self::open_with(root, crate::card::plain)
    }

    /// Like `open` with the function that makes cards from specs.
    pub fn open_with(root: &Path, card_of: CardFn) -> Result<WorkflowDir, String> {
        let p = root.join("workflow.json");
        let text = std::fs::read_to_string(&p).map_err(|e| io_err(&p, e))?;
        let wf: Workflow = serde_json::from_str(&text).map_err(|e| io_err(&p, format!("not a nyx workflow: {e}")))?;
        if wf.nyx > FORMAT {
            return Err(format!("{}: made with a newer nyx (format {}, this one reads {FORMAT})", p.display(), wf.nyx));
        }
        Ok(WorkflowDir { root: root.to_path_buf(), wf, card_of })
    }

    pub fn save(&self) -> Result<(), String> {
        let text = serde_json::to_string_pretty(&self.wf).map_err(|e| e.to_string())?;
        write_atomic(&self.manifest_path(), text.as_bytes()).map_err(|e| io_err(&self.manifest_path(), e))
    }

    /// Edited, ran or noted: the dashboard's three-day clock restarts. Looking at a workflow never calls this.
    pub fn touch(&mut self) -> Result<(), String> {
        self.wf.touched = iso(now_secs());
        self.save()
    }

    fn write_gitignore(&self) -> Result<(), String> {
        let p = self.root.join(".gitignore");
        if p.exists() {
            return Ok(());
        }
        write_atomic(&p, b"# meshes and heavy results stay out of git\n*/runs/*/*.h5\n*/runs/*/*.hdf5\n*/runs/*/*.mesh\n*.tmp*\n").map_err(|e| io_err(&p, e))
    }

    // ------------------------------------------------------------------------------------------------------ experiments

    fn freeze(&self, text: &str) -> Result<String, String> {
        let h = hash_bytes(text.as_bytes());
        let p = self.frozen_path(&h);
        if !p.exists() {
            write_atomic(&p, text.as_bytes()).map_err(|e| io_err(&p, e))?;
        }
        Ok(h)
    }

    pub fn read_spec(&self, eid: &str) -> Result<String, String> {
        let p = self.spec_path(eid);
        std::fs::read_to_string(&p).map_err(|e| io_err(&p, e))
    }

    pub fn read_frozen(&self, hash: &str) -> Result<String, String> {
        let p = self.frozen_path(hash);
        std::fs::read_to_string(&p).map_err(|e| io_err(&p, e))
    }

    /// Adds an experiment (a fork when `parent` is given) and records the first version of it.
    pub fn add_experiment(&mut self, title: &str, spec_json: &str, parent: Option<&str>) -> Result<String, String> {
        let value = parse_spec(spec_json)?;
        let id = format!("e{}", self.wf.experiments.len() + 1);
        let hash = self.freeze(spec_json)?;
        write_atomic(&self.spec_path(&id), spec_json.as_bytes()).map_err(|e| io_err(&self.spec_path(&id), e))?;
        self.wf.experiments.push(Experiment {
            id: id.clone(),
            title: title.into(),
            parent: parent.map(String::from),
            links: vec![],
            spec_hash: hash.clone(),
            card: (self.card_of)(&value),
            runs: Digest::default(),
            model: None,
        });
        let summary = match parent {
            Some(p) => format!("Forked from {p}"),
            None => "Created".into(),
        };
        self.push_version(&id, hash, None, summary, vec![], false)?;
        self.touch()?;
        Ok(id)
    }

    /// Copies an experiment's current spec and runs nothing: the fork remembers its parent.
    pub fn fork(&mut self, eid: &str, title: &str) -> Result<String, String> {
        let text = self.read_spec(eid)?;
        self.add_experiment(title, &text, Some(eid))
    }

    pub fn link(&mut self, from: &str, to: &str, label: &str) -> Result<(), String> {
        if self.wf.experiment(to).is_none() || from == to {
            return Err("link two different experiments of this workflow".into());
        }
        let e = self.wf.experiment_mut(from).ok_or("no such experiment")?;
        e.links.retain(|l| l.to != to);
        e.links.push(Link { to: to.into(), label: label.into() });
        self.touch()
    }

    // ------------------------------------------------------------------------------------------------------ versions

    pub fn history(&self) -> Vec<Version> {
        std::fs::read_to_string(self.history_path()).map(|t| versions::parse_history(&t)).unwrap_or_default()
    }

    fn write_history(&self, v: &[Version]) -> Result<(), String> {
        write_atomic(&self.history_path(), versions::render_history(v).as_bytes()).map_err(|e| io_err(&self.history_path(), e))
    }

    fn push_version(&mut self, eid: &str, hash: String, prev: Option<String>, summary: String, changes: Vec<versions::Change>, merge: bool) -> Result<u32, String> {
        let mut h = self.history();
        let now = iso(now_secs());
        if merge {
            if let Some(last) = h.last().cloned() {
                let recent = parse_iso(&last.t).is_some_and(|t| now_secs() - t < MERGE_SECS);
                let has_runs = runs::read(&self.runs_path(eid)).0.iter().any(|r| r.version == last.n);
                if last.kind == VersionKind::Auto && last.experiment == eid && last.prev_hash.is_some() && recent && !has_runs {
                    // re-diff against the spec the merged version started from
                    let (prev, changes, summary) = match last.prev_hash.as_deref().and_then(|p| self.read_frozen(p).ok()) {
                        Some(old) => {
                            let ch = versions::diff(&parse_spec(&old)?, &parse_spec(&self.read_frozen(&hash)?)?);
                            let s = versions::summarize(&ch);
                            (last.prev_hash.clone(), ch, s)
                        }
                        None => (prev, changes, summary),
                    };
                    let n = last.n;
                    *h.last_mut().unwrap() = Version { n, t: now, kind: VersionKind::Auto, name: None, experiment: eid.into(), spec_hash: hash, prev_hash: prev, summary, changes };
                    self.write_history(&h)?;
                    return Ok(n);
                }
            }
        }
        let n = h.last().map(|v| v.n).unwrap_or(0) + 1;
        h.push(Version { n, t: now, kind: VersionKind::Auto, name: None, experiment: eid.into(), spec_hash: hash, prev_hash: prev, summary, changes });
        self.write_history(&h)?;
        self.wf.version = n;
        Ok(n)
    }

    /// Stores the current spec of `eid` when it differs from the saved one: freezes it, regenerates the card, adds (or merges into) an
    /// automatic version and touches the workflow. Returns the version number, `None` when nothing changed.
    pub fn commit_spec(&mut self, eid: &str, spec_json: &str) -> Result<Option<u32>, String> {
        let new_hash = hash_bytes(spec_json.as_bytes());
        let old_hash = self.wf.experiment(eid).ok_or_else(|| format!("no experiment {eid}"))?.spec_hash.clone();
        if new_hash == old_hash {
            return Ok(None);
        }
        let new_value = parse_spec(spec_json)?;
        let old_value = self.read_frozen(&old_hash).ok().and_then(|x| parse_spec(&x).ok());
        let changes = old_value.as_ref().map(|o| versions::diff(o, &new_value)).unwrap_or_default();
        let summary = versions::summarize(&changes);
        self.freeze(spec_json)?;
        write_atomic(&self.spec_path(eid), spec_json.as_bytes()).map_err(|e| io_err(&self.spec_path(eid), e))?;
        let e = self.wf.experiment_mut(eid).unwrap();
        e.spec_hash = new_hash.clone();
        e.card = (self.card_of)(&new_value);
        let n = self.push_version(eid, new_hash, Some(old_hash), summary, changes, true)?;
        self.wf.version = self.history().last().map(|v| v.n).unwrap_or(n);
        self.touch()?;
        Ok(Some(n))
    }

    /// Gives version `n` a name; named versions are kept as they are.
    pub fn name_version(&mut self, n: u32, name: &str) -> Result<(), String> {
        let mut h = self.history();
        let v = h.iter_mut().find(|v| v.n == n).ok_or("no such version")?;
        v.kind = VersionKind::Named;
        v.name = Some(name.trim().to_string()).filter(|s| !s.is_empty());
        self.write_history(&h)?;
        self.touch()
    }

    /// The spec text of version `n`.
    pub fn version_spec(&self, n: u32) -> Result<String, String> {
        let h = self.history();
        let v = h.iter().find(|v| v.n == n).ok_or("no such version")?;
        self.read_frozen(&v.spec_hash)
    }

    /// Makes version `n`'s spec the current one of its experiment, as a new version ("Restore this version (adds v13)").
    pub fn restore(&mut self, n: u32) -> Result<(String, String), String> {
        let h = self.history();
        let v = h.iter().find(|v| v.n == n).ok_or("no such version")?.clone();
        let text = self.read_frozen(&v.spec_hash)?;
        let eid = v.experiment.clone();
        let cur = self.wf.experiment(&eid).ok_or("that experiment is gone")?.spec_hash.clone();
        if cur != v.spec_hash {
            let (old, new) = (parse_spec(&self.read_frozen(&cur)?)?, parse_spec(&text)?);
            let changes = versions::diff(&old, &new);
            write_atomic(&self.spec_path(&eid), text.as_bytes()).map_err(|e| io_err(&self.spec_path(&eid), e))?;
            let e = self.wf.experiment_mut(&eid).unwrap();
            e.spec_hash = v.spec_hash.clone();
            e.card = (self.card_of)(&new);
            self.push_version(&eid, v.spec_hash.clone(), Some(cur), format!("Restored v{n}"), changes, false)?;
            self.touch()?;
        }
        Ok((eid, text))
    }

    // ------------------------------------------------------------------------------------------------------ runs

    pub fn runs(&self, eid: &str) -> (Vec<RunRecord>, usize) {
        runs::read(&self.runs_path(eid))
    }

    /// Adds a run to `eid`: freezes the spec it solved, appends the line and refreshes the digest.
    pub fn record_run(&mut self, eid: &str, mut r: RunRecord, spec_json: &str) -> Result<RunRecord, String> {
        r.spec = self.freeze(spec_json)?;
        let (existing, _) = self.runs(eid);
        if r.id.is_empty() {
            r.id = runs::next_id(&existing);
        }
        if r.t.is_empty() {
            r.t = iso(now_secs());
        }
        if r.femcore.is_empty() {
            r.femcore = env!("CARGO_PKG_VERSION").into();
        }
        r.version = self.wf.version;
        runs::append(&self.runs_path(eid), &r)?;
        // one manifest write for the digest and the touch (each write is a temporary file and a rename: 2 to 5 ms on a machine with a scanner)
        self.wf.touched = iso(now_secs());
        self.rebuild_digest(eid)?;
        Ok(r)
    }

    pub fn patch_run(&mut self, eid: &str, rid: &str, p: RunPatch) -> Result<(), String> {
        runs::append_patch(&self.runs_path(eid), rid, p)?;
        self.rebuild_digest(eid)?;
        self.touch()
    }

    pub fn rebuild_digest(&mut self, eid: &str) -> Result<(), String> {
        let (rs, _) = self.runs(eid);
        let d = runs::digest(&rs);
        self.wf.experiment_mut(eid).ok_or("no such experiment")?.runs = d;
        self.save()
    }

    /// True when the digest of any experiment no longer matches `runs.jsonl` (a crash between the append and the manifest write).
    pub fn digest_is_stale(&self) -> bool {
        self.wf.experiments.iter().any(|e| runs::digest(&self.runs(&e.id).0) != e.runs)
    }

    /// The main metric of an experiment as a series over its ok runs, for the dashboard sparkline: `neff_re` or the first metric.
    pub fn sparkline(&self, eid: &str) -> Vec<f64> {
        let (rs, _) = self.runs(eid);
        let ok: Vec<&RunRecord> = rs.iter().filter(|r| r.status == runs::RunStatus::Ok).collect();
        let key = if ok.iter().any(|r| r.metrics.contains_key("neff_re")) { Some("neff_re".to_string()) } else { ok.iter().find_map(|r| r.metrics.keys().find(|k| !matches!(k.as_str(), "dofs" | "n_solutions")).cloned()) };
        let Some(key) = key else { return vec![] };
        let mut v: Vec<f64> = ok.iter().filter_map(|r| r.metrics.get(&key).copied()).collect();
        let skip = v.len().saturating_sub(40);
        v.drain(..skip);
        v
    }

    // ------------------------------------------------------------------------------------------------------ workflow properties

    pub fn set_title(&mut self, t: &str) -> Result<(), String> {
        self.wf.title = t.trim().to_string();
        self.touch()
    }
    pub fn set_tags(&mut self, tags: Vec<String>) -> Result<(), String> {
        self.wf.tags = tags;
        self.touch()
    }
    pub fn set_experiment_title(&mut self, eid: &str, t: &str) -> Result<(), String> {
        self.wf.experiment_mut(eid).ok_or("no such experiment")?.title = t.trim().to_string();
        self.touch()
    }
    /// Pinning and marking done are bookkeeping: they do not restart the three-day clock.
    pub fn set_pinned(&mut self, p: bool) -> Result<(), String> {
        self.wf.pinned = p;
        self.save()
    }
    pub fn set_status(&mut self, s: Status) -> Result<(), String> {
        self.wf.status = s;
        self.save()
    }
    pub fn notes(&self) -> String {
        std::fs::read_to_string(self.notes_path()).unwrap_or_default()
    }
    pub fn set_notes(&mut self, text: &str) -> Result<(), String> {
        write_atomic(&self.notes_path(), text.as_bytes()).map_err(|e| io_err(&self.notes_path(), e))?;
        self.touch()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use serde_json::json;

    /// A minimal project (JSON) with a title, for tests that need a spec: `w` and `n_modes` make two specs differ.
    pub fn spec_json(w: f64, n_modes: usize) -> String {
        json!({ "title": format!("rectangle {w:e} m, {n_modes} modes"), "width": w, "n_modes": n_modes }).to_string()
    }

    pub fn tmp(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("nyx-store-{tag}-{}-{}", std::process::id(), short_id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn a_workflow_is_created_edited_run_and_reopened() {
        let ws = tmp("flow");
        let mut d = WorkflowDir::create(&ws, "SOI TE modes sweep", "gui", &spec_json(5e-7, 4)).unwrap();
        assert!(d.root.file_name().unwrap().to_str().unwrap().starts_with("soi-te-modes-sweep_"));
        assert_eq!(d.wf.version, 1);
        assert_eq!(d.history().len(), 1);
        let e1 = d.wf.experiments[0].id.clone();
        assert!(d.wf.experiments[0].card.summary.contains("rectangle"), "{}", d.wf.experiments[0].card.summary);

        // an unchanged spec is not a new version; a changed one is, and a quick second change merges into it
        assert_eq!(d.commit_spec(&e1, &spec_json(5e-7, 4)).unwrap(), None);
        assert_eq!(d.commit_spec(&e1, &spec_json(5.5e-7, 4)).unwrap(), Some(2));
        assert_eq!(d.commit_spec(&e1, &spec_json(5.5e-7, 6)).unwrap(), Some(2), "merged");
        let h = d.history();
        assert_eq!(h.len(), 2);
        assert!(h[1].changes.iter().any(|c| c.path.ends_with("width")) && h[1].changes.iter().any(|c| c.path.ends_with("n_modes")), "{:?}", h[1].changes);

        // a run on v2 stops later edits merging into it
        let rec = RunRecord { sec: 0.11, params: [("w".to_string(), 5.5e-7)].into(), metrics: [("neff_re".to_string(), 2.44)].into(), ..Default::default() };
        let spec_now = d.read_spec(&e1).unwrap();
        let rec = d.record_run(&e1, rec, &spec_now).unwrap();
        assert!(rec.id.starts_with("r0001_") && rec.version == 2 && rec.spec.starts_with("b3:"));
        assert_eq!(d.commit_spec(&e1, &spec_json(6e-7, 6)).unwrap(), Some(3));
        assert_eq!(d.history().len(), 3);
        assert_eq!(d.wf.experiments[0].runs.count, 1);
        assert!((d.wf.experiments[0].runs.metrics["neff_re"].max - 2.44).abs() < 1e-12);

        // naming, restoring (adds a version), fork and link
        d.name_version(2, "Baseline").unwrap();
        let (eid, text) = d.restore(2).unwrap();
        assert_eq!(eid, e1);
        assert_eq!(hash_bytes(text.as_bytes()), d.wf.experiments[0].spec_hash);
        assert_eq!(d.history().last().unwrap().summary, "Restored v2");
        let e2 = d.fork(&e1, "TE modes vs thickness").unwrap();
        assert_eq!(d.wf.experiment(&e2).unwrap().parent.as_deref(), Some("e1"));
        d.link(&e2, &e1, "validates").unwrap();
        assert!(d.link(&e1, &e1, "x").is_err());

        // reopen from disk: same manifest; a manifest that missed a run is detected and repaired
        let reopened = WorkflowDir::open(&d.root).unwrap();
        assert_eq!(reopened.wf, d.wf);
        let mut stale = WorkflowDir::open(&d.root).unwrap();
        stale.wf.experiments[0].runs = Digest::default();
        assert!(stale.digest_is_stale());
        stale.rebuild_digest(&e1).unwrap();
        assert!(!stale.digest_is_stale());
        assert!(d.root.join(".gitignore").exists());
        assert_eq!(d.sparkline(&e1), [2.44]);
        let _ = std::fs::remove_dir_all(&ws);
    }

    fn json_card(v: &serde_json::Value) -> Card {
        Card { summary: format!("engine card: {}", v["title"].as_str().unwrap_or("")), ..Default::default() }
    }

    #[test]
    fn json_specs_and_a_card_function_of_the_caller_work_without_the_old_project_shape() {
        let ws = tmp("json");
        let spec = |t: &str| serde_json::json!({ "title": t, "out": r"C:\Users\me\r.h5" }).to_string();
        let mut d = WorkflowDir::create_with(&ws, None, "json", "ai", &spec("a"), json_card).unwrap();
        let e = d.wf.experiments[0].id.clone();
        assert_eq!(d.wf.experiments[0].card.summary, "engine card: a");
        assert_eq!(d.commit_spec(&e, &spec("b")).unwrap(), Some(2));
        assert_eq!(d.wf.experiments[0].card.summary, "engine card: b");
        assert!(d.history()[1].changes.iter().any(|c| c.path.ends_with("title")), "{:?}", d.history()[1].changes);
        let reopened = WorkflowDir::open_with(&d.root, json_card).unwrap();
        assert_eq!(reopened.wf, d.wf);
        let (out, n) = crate::archive::strip_paths(&spec("b")).unwrap();
        assert_eq!(n, 1);
        assert!(out.contains("r.h5") && !out.contains("Users"));
        assert!(d.commit_spec(&e, "{not json").is_err());
        let _ = std::fs::remove_dir_all(&ws);
    }

    #[test]
    fn a_newer_format_is_refused_not_misread() {
        let ws = tmp("newer");
        let d = WorkflowDir::create(&ws, "x", "gui", &spec_json(5e-7, 4)).unwrap();
        let text = std::fs::read_to_string(d.manifest_path()).unwrap().replace("\"nyx\": 1", "\"nyx\": 99");
        std::fs::write(d.manifest_path(), text).unwrap();
        assert!(WorkflowDir::open(&d.root).unwrap_err().contains("newer nyx"));
        let _ = std::fs::remove_dir_all(&ws);
    }
}
