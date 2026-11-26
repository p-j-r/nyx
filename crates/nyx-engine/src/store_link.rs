// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The engine's persistence, through `nyx-store`: a project is the spec of one experiment of a workflow folder, every solve is a run in its
//! run table, and the search card and the run metrics are made here, from the project and the results, because the store does not know what
//! a project or a result looks like (R21, R25). `StoreLink` is the project's binding to a workflow folder.

use crate::project::Project;
use nyx_physics::{Error, Scalar, SolveResult};
use nyx_store::runs::{self, Artifact, RunRecord, RunStatus};
use nyx_store::{Card, Dim, GeoCard, WorkflowDir};
use std::collections::BTreeMap;
use std::path::Path;

fn store_err(e: String) -> Error {
    Error::invalid(e)
}

/// The search card of a project (`nyx_store::CardFn`): what the dashboard shows and search finds, generated, never typed.
pub fn card_of(v: &serde_json::Value) -> Card {
    let Ok(p) = serde_json::from_value::<Project>(v.clone()) else {
        return Card { summary: "unreadable project".into(), ..Card::default() };
    };
    let regions: Vec<String> = p
        .geometry
        .features
        .iter()
        .filter_map(|f| match &f.op {
            geom::model::Op::Overlay { regions, .. } => {
                Some(regions.iter().map(|r| r.name.clone()).collect::<Vec<_>>())
            }
            _ => None,
        })
        .flatten()
        .collect();
    let active = p.active.clone().or_else(|| p.physics.first().map(|s| s.physics.clone())).unwrap_or_default();
    let spec = p.spec(&active);
    let study = spec.map(|s| s.study.kind.clone()).unwrap_or_else(|| "none".into());
    let wavelength_m = spec.and_then(|s| s.study.params.get("wavelength").and_then(|w| w.as_f64()));
    let mut materials = vec![];
    for s in &p.physics {
        for (region, m) in &s.materials {
            materials.push(match m {
                nyx_physics::MaterialAssign::Named(n) => format!("{region}: {n}"),
                nyx_physics::MaterialAssign::Oriented { material, euler_deg } => {
                    format!("{region}: {material} at {euler_deg:?} deg")
                }
                nyx_physics::MaterialAssign::Props(_) | nyx_physics::MaterialAssign::PropsTensor(_) => {
                    format!("{region}: inline properties")
                }
            });
        }
    }
    let dims: Vec<Dim> = p
        .geometry
        .params
        .iter()
        .filter_map(|q| match q.value {
            geom::expr::Expr::Num(m) => Some(Dim { name: q.name.clone(), of: "parameter".into(), m }),
            _ => None,
        })
        .collect();
    let swept = p.geometry.dim == 2 && p.mesh.sweep.as_ref().is_some_and(|s| !s.is_empty());
    let effective_dim = if swept { 3 } else { p.geometry.dim };
    let mut keywords: Vec<String> = p.physics.iter().map(|s| s.physics.clone()).collect();
    keywords.extend(p.physics.iter().map(|s| s.study.kind.clone()));
    keywords.extend(regions.iter().map(|r| r.to_lowercase()));
    keywords.extend(p.geometry.features.iter().map(|f| crate::geomdesc::kind_of(&f.op).to_lowercase()));
    if swept {
        keywords.push("3d".into());
        keywords.push("swept".into());
    }
    keywords.sort();
    keywords.dedup();
    let wl = wavelength_m.map(|w| format!(", {:.0} nm", w * 1e9)).unwrap_or_default();
    let dim_text = if swept { "3D (swept)".to_string() } else { format!("{effective_dim}D") };
    let summary = format!(
        "{}, {} feature{}, {} region{}, {} {}{}",
        dim_text,
        p.geometry.features.len(),
        if p.geometry.features.len() == 1 { "" } else { "s" },
        regions.len(),
        if regions.len() == 1 { "" } else { "s" },
        if active.is_empty() { "no physics" } else { &active },
        study,
        wl
    );
    Card {
        physics: active,
        study,
        dim: effective_dim,
        wavelength_m,
        geometry: GeoCard {
            family: "feature list".into(),
            summary: format!("{} features, regions {}", p.geometry.features.len(), regions.join(", ")),
            dims,
        },
        materials,
        keywords,
        summary,
    }
}

/// The numbers of one result that go in the run table: the first solution's value and eigenvalue, the counts, and every numeric scalar the
/// physics reported (name made into `snake_case`). Formatted text scalars are never parsed.
pub fn metrics_of(res: &SolveResult) -> BTreeMap<String, f64> {
    let mut m = BTreeMap::new();
    if let Some(s) = res.solutions.first() {
        if let Some(v) = s.value {
            m.insert("value_1".into(), v);
        }
        if let Some([re, im]) = s.eigenvalue {
            m.insert("eigenvalue_re".into(), re);
            m.insert("eigenvalue_im".into(), im);
        }
    }
    m.insert("n_solutions".into(), res.solutions.len() as f64);
    m.insert("dofs".into(), res.dofs as f64);
    m.insert("seconds".into(), res.seconds);
    for (name, s) in &res.scalars {
        if let Scalar::Num(x) = s {
            if x.is_finite() {
                let key: String = name
                    .chars()
                    .map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' })
                    .collect();
                let key = key.split('_').filter(|s| !s.is_empty()).collect::<Vec<_>>().join("_");
                if !key.is_empty() {
                    m.entry(key).or_insert(*x);
                }
            }
        }
    }
    m
}

/// How a run ended.
pub enum RunOutcome<'a> {
    Solved(&'a [SolveResult]),
    Failed(&'a Error),
}

/// The binding of a project to one experiment of a workflow folder.
pub struct StoreLink {
    pub dir: WorkflowDir,
    pub experiment: String,
    /// the script (name, text) the runs recorded now are made by, when a script drives the engine: its text is kept with each run (S4)
    pub script: Option<(String, String)>,
}

impl StoreLink {
    /// Makes a new workflow folder in `workspace` holding `project` as its first experiment.
    pub fn create(workspace: &Path, project: &Project, origin: &str) -> Result<StoreLink, Error> {
        let dir = WorkflowDir::create_with(workspace, None, &project.title, origin, &project.to_json(), card_of)
            .map_err(store_err)?;
        let experiment = dir.wf.experiments[0].id.clone();
        Ok(StoreLink { dir, experiment, script: None })
    }

    /// Opens a workflow folder and reads the project of its first experiment (or of `experiment`).
    pub fn open(root: &Path, experiment: Option<&str>) -> Result<(StoreLink, Project), Error> {
        let dir = WorkflowDir::open_with(root, card_of).map_err(store_err)?;
        let experiment = match experiment {
            Some(e) => dir
                .wf
                .experiment(e)
                .map(|e| e.id.clone())
                .ok_or_else(|| Error::invalid(format!("no experiment {e}")))?,
            None => dir
                .wf
                .experiments
                .first()
                .map(|e| e.id.clone())
                .ok_or_else(|| Error::invalid("the workflow has no experiment"))?,
        };
        let text = dir.read_spec(&experiment).map_err(store_err)?;
        let project = Project::from_json(&text)?;
        Ok((StoreLink { dir, experiment, script: None }, project))
    }

    /// Stores the project when it differs from the saved one; returns the new version number, `None` when nothing changed.
    pub fn commit(&mut self, project: &Project) -> Result<Option<u32>, Error> {
        self.dir.commit_spec(&self.experiment, &project.to_json()).map_err(store_err)
    }

    /// Adds a run to the table: the project as it was solved (frozen), the metrics of the results and the stage log. The heavy files
    /// (`results.h5` with the project, the mesh and all fields, and `results.xdmf`, in a build with HDF5; `log.txt`) go in the run's folder, which can be deleted without losing the record. The fields are never written as JSON (ADR 62, 63).
    pub fn record(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<&femcore::mesh::SimplexMesh>,
    ) -> Result<RunRecord, Error> {
        self.record_with_options(project, outcome, seconds, log, mesh, true, &[])
    }

    /// `record` with explicit `fields_saved` and optional `post_analyses` results.
    #[allow(clippy::too_many_arguments)]
    pub fn record_with_options(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<&femcore::mesh::SimplexMesh>,
        fields_saved: bool,
        post_analyses: &[serde_json::Value],
    ) -> Result<RunRecord, Error> {
        self.commit(project)?;
        self.record_frozen(project, outcome, seconds, log, mesh, None, fields_saved, post_analyses)
    }

    /// Adds a run for a project that is not the experiment's own (a point of a sweep): the run freezes the project it solved, and the
    /// experiment's project, its versions and its history are left as they are.
    pub fn record_point(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<&femcore::mesh::SimplexMesh>,
        reuse: Option<&serde_json::Value>,
    ) -> Result<RunRecord, Error> {
        self.record_frozen(project, outcome, seconds, log, mesh, reuse, true, &[])
    }

    /// `record_point` with explicit `fields_saved` and optional `post_analyses` results.
    #[allow(clippy::too_many_arguments)]
    pub fn record_point_with_options(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<&femcore::mesh::SimplexMesh>,
        reuse: Option<&serde_json::Value>,
        fields_saved: bool,
        post_analyses: &[serde_json::Value],
    ) -> Result<RunRecord, Error> {
        self.record_frozen(project, outcome, seconds, log, mesh, reuse, fields_saved, post_analyses)
    }

    #[allow(clippy::too_many_arguments)]
    fn record_frozen(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
        mesh: Option<&femcore::mesh::SimplexMesh>,
        reuse: Option<&serde_json::Value>,
        fields_saved: bool,
        post_analyses: &[serde_json::Value],
    ) -> Result<RunRecord, Error> {
        let (existing, _) = self.dir.runs(&self.experiment);
        let id = runs::next_id(&existing);
        let mut rec = RunRecord {
            id: id.clone(),
            sec: seconds,
            femcore: format!("nyx-engine {}", env!("CARGO_PKG_VERSION")),
            reuse: reuse.cloned(),
            ..Default::default()
        };
        for p in project.geometry.params.iter().filter_map(|p| match p.value {
            geom::expr::Expr::Num(v) => Some((p.name.clone(), v)),
            _ => None,
        }) {
            rec.params.insert(p.0, p.1);
        }
        let mut log: Vec<String> = log.to_vec();
        let run_dir = self.dir.run_dir(&self.experiment, &id);
        let mut artifacts = vec![];
        #[allow(unused_mut)]
        let mut extra: Vec<Artifact> = vec![];
        let mut write = |name: &str, kind: &str, bytes: Vec<u8>| -> Result<(), Error> {
            std::fs::create_dir_all(&run_dir).map_err(|e| Error::invalid(format!("{}: {e}", run_dir.display())))?;
            let path = run_dir.join(name);
            nyx_store::util::write_atomic(&path, &bytes)
                .map_err(|e| Error::invalid(format!("{}: {e}", path.display())))?;
            artifacts.push(Artifact {
                kind: kind.into(),
                path: format!("runs/{id}/{name}"),
                bytes: bytes.len() as u64,
            });
            Ok(())
        };
        match outcome {
            RunOutcome::Solved(results) => {
                rec.status = RunStatus::Ok;
                for r in results {
                    let prefix = if results.len() > 1 { format!("{}_", r.physics) } else { String::new() };
                    for (k, v) in metrics_of(r) {
                        rec.metrics.insert(format!("{prefix}{k}"), v);
                    }
                }
                // the fields are never written as JSON text (ADR 62, 63): without the feature a run keeps its metrics and its log only
                #[cfg(not(feature = "hdf5"))]
                let log_note = Some(
                    "the fields are not stored: this build has no HDF5 support (cargo feature `hdf5`)".to_string(),
                );
                #[cfg(feature = "hdf5")]
                let log_note: Option<String> = None;
                if let Some(n) = log_note {
                    log.push(n);
                }
                #[cfg(feature = "hdf5")]
                if let Some(mesh) = mesh {
                    std::fs::create_dir_all(&run_dir)
                        .map_err(|e| Error::invalid(format!("{}: {e}", run_dir.display())))?;
                    let path = run_dir.join("results.h5");
                    let refs: Vec<(String, &SolveResult)> = results.iter().map(|r| (r.physics.clone(), r)).collect();
                    crate::h5io::write_results_with_options(
                        &path,
                        project,
                        mesh,
                        &refs,
                        crate::h5io::FieldStorageLayout::L2ShuffleLzfPacked,
                        fields_saved,
                    )
                    .map_err(Error::invalid)?;
                    if fields_saved {
                        let x = crate::xdmf::write_xdmf(&path, mesh, &refs).map_err(Error::invalid)?;
                        extra.push(Artifact {
                            kind: "result".into(),
                            path: format!("runs/{id}/results.xdmf"),
                            bytes: std::fs::metadata(&x).map(|m| m.len()).unwrap_or(0),
                        });
                    }
                    extra.push(Artifact {
                        kind: "result".into(),
                        path: format!("runs/{id}/results.h5"),
                        bytes: crate::h5io::file_bytes(&path),
                    });
                }
                #[cfg(not(feature = "hdf5"))]
                let _ = (mesh, fields_saved);
                if !post_analyses.is_empty() {
                    let bytes = serde_json::to_vec(post_analyses).map_err(|e| Error::invalid(e.to_string()))?;
                    write("post_analyses.json", "analysis", bytes)?;
                }
            }
            RunOutcome::Failed(e) => {
                rec.status = if matches!(e, Error::Cancelled) { RunStatus::Cancelled } else { RunStatus::Failed };
                rec.error = Some(e.to_string());
            }
        }
        if !log.is_empty() {
            write("log.txt", "log", format!("{}\n", log.join("\n")).into_bytes())?;
        }
        // the script that made this run, with its hash in the file name (`script-b3:...` would not be a Windows file name: the hash is hex)
        if let Some((_, text)) = &self.script {
            let hash = nyx_store::util::hash_bytes(text.as_bytes());
            write(&format!("script-{}.nys", hash.trim_start_matches("b3:")), "script", text.clone().into_bytes())?;
        }
        artifacts.extend(extra);
        rec.artifacts = artifacts;
        let spec = project.to_json();
        self.dir.record_run(&self.experiment, rec, &spec).map_err(store_err)
    }

    /// The run table of the bound experiment, oldest first.
    pub fn runs(&self) -> Vec<RunRecord> {
        self.dir.runs(&self.experiment).0
    }
}

#[cfg(test)]
mod tests {
    use super::card_of;
    use crate::project::{Project, SweepLayer};
    use serde_json::json;

    #[test]
    fn swept_2d_projects_report_effective_3d_card_and_search_keywords() {
        let base = Project::new("Waveguide", 2);
        let plain_val: serde_json::Value = serde_json::from_str(&base.to_json()).unwrap();
        let plain = card_of(&plain_val);
        assert_eq!(plain.dim, 2);
        assert!(plain.summary.starts_with("2D, "), "{}", plain.summary);
        assert!(!plain.keywords.contains(&"swept".to_string()));

        let mut empty_sweep = base.clone();
        empty_sweep.mesh.sweep = Some(vec![]);
        let empty_val: serde_json::Value = serde_json::from_str(&empty_sweep.to_json()).unwrap();
        let empty_card = card_of(&empty_val);
        assert_eq!(empty_card.dim, 2);
        assert!(empty_card.summary.starts_with("2D, "), "{}", empty_card.summary);

        let mut swept = base;
        swept.mesh.sweep =
            Some(vec![SweepLayer { thickness: 2.2e-7, n: 4, domain_offset: 0, fill: None, regions: None }]);
        let swept_val: serde_json::Value = serde_json::from_str(&swept.to_json()).unwrap();
        let card = card_of(&swept_val);
        assert_eq!(card.dim, 3);
        assert!(card.summary.starts_with("3D (swept), "), "{}", card.summary);
        assert!(card.keywords.contains(&"3d".to_string()));
        assert!(card.keywords.contains(&"swept".to_string()));

        let native_3d = Project::new("Block", 3);
        let native_val = serde_json::to_value(&native_3d).unwrap_or_else(|_| json!({}));
        let native_card = card_of(&native_val);
        assert_eq!(native_card.dim, 3);
        assert!(native_card.summary.starts_with("3D, "), "{}", native_card.summary);
        assert!(!native_card.summary.contains("swept"));
    }

    #[test]
    fn record_with_options_records_script_nys_and_xdmf_consistently() {
        let ws = std::env::temp_dir().join(format!("nyx-store-link-art-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&ws);
        std::fs::create_dir_all(&ws).unwrap();

        let mut p = Project::new("ArtTest", 2);
        p.mesh.h = 0.25;
        let (mut link, _) = (super::StoreLink::create(&ws, &p, "ai").unwrap(), ());
        link.script = Some(("demo.nys".into(), "nyx.project_info()".into()));
        let mesh = femcore::mesh::SimplexMesh::new(
            2,
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            vec![0, 1, 2],
            vec![1],
            vec![0, 1, 1, 2, 2, 0],
            vec![1, 1, 1],
        );
        let res = nyx_physics::SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            order: 1,
            dim: 2,
            dofs: 3,
            seconds: 0.01,
            solutions: vec![nyx_physics::Solution { label: "s0".into(), value: Some(1.0), ..Default::default() }],
            scalars: vec![],
            log: vec!["ok".into()],
            ..Default::default()
        };

        let r_full = link
            .record_with_options(
                &p,
                super::RunOutcome::Solved(std::slice::from_ref(&res)),
                0.01,
                &["log".into()],
                Some(&mesh),
                true,
                &[],
            )
            .unwrap();
        assert!(r_full.artifacts.iter().any(|a| a.kind == "script" && a.path.ends_with(".nys")));
        if cfg!(feature = "hdf5") {
            assert!(r_full.artifacts.iter().any(|a| a.path.ends_with("/results.xdmf")));
            assert!(r_full.artifacts.iter().any(|a| a.path.ends_with("/results.h5")));
        }

        let r_light = link
            .record_with_options(
                &p,
                super::RunOutcome::Solved(std::slice::from_ref(&res)),
                0.01,
                &["log".into()],
                Some(&mesh),
                false,
                &[],
            )
            .unwrap();
        assert!(r_light.artifacts.iter().any(|a| a.kind == "script" && a.path.ends_with(".nys")));
        assert!(!r_light.artifacts.iter().any(|a| a.path.ends_with("/results.xdmf")));

        let _ = std::fs::remove_dir_all(&ws);
    }
}
