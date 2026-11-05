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
                nyx_physics::MaterialAssign::Oriented { material, euler_deg } => format!("{region}: {material} at {euler_deg:?} deg"),
                nyx_physics::MaterialAssign::Props(_) => format!("{region}: inline properties"),
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
    let mut keywords: Vec<String> = p.physics.iter().map(|s| s.physics.clone()).collect();
    keywords.extend(p.physics.iter().map(|s| s.study.kind.clone()));
    keywords.extend(regions.iter().map(|r| r.to_lowercase()));
    keywords.extend(p.geometry.features.iter().map(|f| crate::geomdesc::kind_of(&f.op).to_lowercase()));
    keywords.sort();
    keywords.dedup();
    let wl = wavelength_m.map(|w| format!(", {:.0} nm", w * 1e9)).unwrap_or_default();
    let summary = format!(
        "{}D, {} features, {} regions, {} {}{}",
        p.geometry.dim,
        p.geometry.features.len(),
        regions.len(),
        if active.is_empty() { "no physics" } else { &active },
        study,
        wl
    );
    Card {
        physics: active,
        study,
        dim: p.geometry.dim,
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
}

impl StoreLink {
    /// Makes a new workflow folder in `workspace` holding `project` as its first experiment.
    pub fn create(workspace: &Path, project: &Project, origin: &str) -> Result<StoreLink, Error> {
        let dir = WorkflowDir::create_with(workspace, None, &project.title, origin, &project.to_json(), card_of)
            .map_err(store_err)?;
        let experiment = dir.wf.experiments[0].id.clone();
        Ok(StoreLink { dir, experiment })
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
        Ok((StoreLink { dir, experiment }, project))
    }

    /// Stores the project when it differs from the saved one; returns the new version number, `None` when nothing changed.
    pub fn commit(&mut self, project: &Project) -> Result<Option<u32>, Error> {
        self.dir.commit_spec(&self.experiment, &project.to_json()).map_err(store_err)
    }

    /// Adds a run to the table: the project as it was solved (frozen), the metrics of the results and the stage log. The heavy files
    /// (`result_<physics>.json`, `log.txt`) go in the run's folder, which can be deleted without losing the record.
    pub fn record(
        &mut self,
        project: &Project,
        outcome: RunOutcome,
        seconds: f64,
        log: &[String],
    ) -> Result<RunRecord, Error> {
        self.commit(project)?;
        let (existing, _) = self.dir.runs(&self.experiment);
        let id = runs::next_id(&existing);
        let mut rec = RunRecord {
            id: id.clone(),
            sec: seconds,
            femcore: format!("nyx-engine {}", env!("CARGO_PKG_VERSION")),
            ..Default::default()
        };
        for p in project.geometry.params.iter().filter_map(|p| match p.value {
            geom::expr::Expr::Num(v) => Some((p.name.clone(), v)),
            _ => None,
        }) {
            rec.params.insert(p.0, p.1);
        }
        let run_dir = self.dir.run_dir(&self.experiment, &id);
        let mut artifacts = vec![];
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
                    let text = serde_json::to_vec(r).map_err(|e| Error::invalid(e.to_string()))?;
                    write(&format!("result_{}.json", r.physics), "result", text)?;
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
        rec.artifacts = artifacts;
        let spec = project.to_json();
        self.dir.record_run(&self.experiment, rec, &spec).map_err(store_err)
    }

    /// The run table of the bound experiment, oldest first.
    pub fn runs(&self) -> Vec<RunRecord> {
        self.dir.runs(&self.experiment).0
    }
}
