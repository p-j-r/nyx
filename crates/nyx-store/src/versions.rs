//! Versions: every change to an experiment's set-up is kept in `history.jsonl` (one line per version) with the hash of its spec; the
//! specs themselves sit content-addressed in `specs/`, so any version can be diffed, restored or forked from.
//!
//! Automatic versions made within five minutes of each other on the same experiment, with no run recorded on the earlier one, are merged
//! into one (the intermediate states are not kept); a named version is never merged into and nothing else is ever removed.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum VersionKind {
    #[default]
    Auto,
    Named,
}

/// One changed value, by dotted path into the saved project.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Change {
    pub path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub from: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub to: Option<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Section {
    Concept,
    Build,
    Analysis,
}

impl Section {
    pub fn name(self) -> &'static str {
        match self {
            Section::Concept => "CONCEPT",
            Section::Build => "BUILD",
            Section::Analysis => "ANALYSIS",
        }
    }
}

impl Change {
    /// Geometry and materials are Concept; display settings are Analysis; everything else of the spec (study, solver, boundaries,
    /// element order, output) is Build.
    pub fn section(&self) -> Section {
        let p = self.path.as_str();
        let mut parts = p.split('.');
        match (parts.next(), parts.next(), parts.next()) {
            (Some("geometry"), _, _) | (Some("active"), _, _) | (Some("specs"), _, Some("materials")) => Section::Concept,
            (Some("colormap"), _, _) => Section::Analysis,
            _ => Section::Build,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Version {
    pub n: u32,
    pub t: String,
    pub kind: VersionKind,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    pub experiment: String,
    pub spec_hash: String,
    /// the spec this one was made from (`None` for the first of an experiment)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prev_hash: Option<String>,
    pub summary: String,
    pub changes: Vec<Change>,
}

fn fmt_leaf(v: &Value) -> String {
    match v {
        Value::Number(n) => n.as_f64().map(|f| if f != 0.0 && (f.abs() < 1e-3 || f.abs() >= 1e6) { format!("{f:e}") } else { format!("{f}") }).unwrap_or_else(|| n.to_string()),
        Value::String(s) => s.clone(),
        o => o.to_string(),
    }
}

fn flatten(v: &Value, path: &str, out: &mut BTreeMap<String, String>) {
    match v {
        Value::Object(m) => {
            if m.is_empty() {
                out.insert(path.into(), "{}".into());
            }
            for (k, x) in m {
                flatten(x, &if path.is_empty() { k.clone() } else { format!("{path}.{k}") }, out);
            }
        }
        Value::Array(a) => {
            if a.is_empty() {
                out.insert(path.into(), "[]".into());
            }
            for (i, x) in a.iter().enumerate() {
                flatten(x, &if path.is_empty() { i.to_string() } else { format!("{path}.{i}") }, out);
            }
        }
        leaf => {
            out.insert(path.into(), fmt_leaf(leaf));
        }
    }
}

/// Every leaf that differs between two project values, in path order. The paths of spec arrays use the physics position, so the
/// same path in two versions is the same setting.
pub fn diff(a: &Value, b: &Value) -> Vec<Change> {
    let (mut fa, mut fb) = (BTreeMap::new(), BTreeMap::new());
    flatten(a, "", &mut fa);
    flatten(b, "", &mut fb);
    let mut out = vec![];
    for (k, va) in &fa {
        match fb.get(k) {
            Some(vb) if vb == va => {}
            Some(vb) => out.push(Change { path: k.clone(), from: Some(va.clone()), to: Some(vb.clone()) }),
            None => out.push(Change { path: k.clone(), from: Some(va.clone()), to: None }),
        }
    }
    for (k, vb) in &fb {
        if !fa.contains_key(k) {
            out.push(Change { path: k.clone(), from: None, to: Some(vb.clone()) });
        }
    }
    out.sort_by(|x, y| x.path.cmp(&y.path));
    out
}

/// One line for a list of versions: the first two changes, then how many more.
pub fn summarize(changes: &[Change]) -> String {
    if changes.is_empty() {
        return "no change".into();
    }
    let one = |c: &Change| match (&c.from, &c.to) {
        (Some(a), Some(b)) => format!("{} {} to {}", c.path, a, b),
        (None, Some(b)) => format!("{} added: {}", c.path, b),
        (Some(a), None) => format!("{} removed: {}", c.path, a),
        _ => c.path.clone(),
    };
    let mut s = changes.iter().take(2).map(one).collect::<Vec<_>>().join("; ");
    if changes.len() > 2 {
        s += &format!("; and {} more", changes.len() - 2);
    }
    s
}

/// The text of the "what changed" panel: changes grouped under CONCEPT, BUILD, ANALYSIS.
pub fn render(changes: &[Change]) -> String {
    let mut s = String::new();
    for sec in [Section::Concept, Section::Build, Section::Analysis] {
        s += &format!("{}\n", sec.name());
        let mine: Vec<&Change> = changes.iter().filter(|c| c.section() == sec).collect();
        if mine.is_empty() {
            s += "  no changes\n";
        }
        for c in mine {
            match (&c.from, &c.to) {
                (Some(a), Some(b)) => s += &format!("  {:<34}{} -> {}\n", c.path, a, b),
                (None, Some(b)) => s += &format!("  {:<34}added: {}\n", c.path, b),
                (Some(a), None) => s += &format!("  {:<34}removed (was {})\n", c.path, a),
                _ => {}
            }
        }
    }
    s
}

pub fn parse_history(text: &str) -> Vec<Version> {
    text.lines().filter_map(|l| serde_json::from_str::<Version>(l).ok()).collect()
}

pub fn render_history(v: &[Version]) -> String {
    v.iter().map(|x| serde_json::to_string(x).unwrap_or_default() + "\n").collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn diff_lists_changed_added_and_removed_leaves_by_section() {
        let a = json!({"geometry": {"w": 4.5e-7, "h": 2e-7}, "specs": [{"study": {"Modes": {"n_modes": 4}}, "materials": {"default": {"eps_r": 1.0}}}], "colormap": "Viridis"});
        let b = json!({"geometry": {"w": 5e-7, "h": 2e-7, "extra": 1}, "specs": [{"study": {"Modes": {"n_modes": 6}}, "materials": {"default": {"eps_r": 2.0}}}], "colormap": "Viridis"});
        let d = diff(&a, &b);
        let paths: Vec<_> = d.iter().map(|c| (c.path.as_str(), c.section())).collect();
        assert_eq!(
            paths,
            [("geometry.extra", Section::Concept), ("geometry.w", Section::Concept), ("specs.0.materials.default.eps_r", Section::Concept), ("specs.0.study.Modes.n_modes", Section::Build)]
        );
        assert_eq!(d[1].from.as_deref(), Some("4.5e-7"));
        assert_eq!(d[1].to.as_deref(), Some("5e-7"));
        assert!(summarize(&d).ends_with("and 2 more"));
        assert!(render(&d).contains("ANALYSIS\n  no changes"));
        assert!(diff(&a, &a).is_empty());
    }
}
