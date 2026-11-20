// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Run records: one append-only line per run in `runs.jsonl`, the digest derived from them, and the small filter language of the runs table.
//!
//! A run is one Compute. Appending is one `write` of one line on a file opened for append, so a window and a detached worker can both
//! add runs without a lock, and a crash can only lose a half-written last line (which readers skip and count). Star, tag and note edits are
//! further lines `{"id":..,"patch":{..}}`, folded into the record when read; nothing is ever rewritten.

use crate::model::{Digest, Range};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum RunStatus {
    #[default]
    Ok,
    Failed,
    Cancelled,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Artifact {
    pub kind: String,
    /// relative to the experiment folder
    pub path: String,
    #[serde(default)]
    pub bytes: u64,
}

/// What a mesh-convergence check said about the problem a run solved (ADR 91, productivity plan P6). Plain data: the store knows nothing of a project.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Convergence {
    /// the ratio of the finer element size to the run's own
    pub factor: f64,
    /// the largest relative difference accepted
    pub tolerance: f64,
    /// what was compared and how the solutions were matched, in words ("every eigenvalue, value and numeric scalar (4), solutions matched by index")
    pub observable: String,
    /// every compared quantity was within the tolerance
    pub passed: bool,
    /// the largest relative difference found
    #[serde(default)]
    pub worst_relative: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct RunRecord {
    pub id: String,
    pub t: String,
    pub sec: f64,
    pub status: RunStatus,
    /// hash of the frozen spec this run solved
    pub spec: String,
    pub femcore: String,
    /// version of the workflow's set-up this run was made on
    #[serde(default)]
    pub version: u32,
    pub params: BTreeMap<String, f64>,
    pub metrics: BTreeMap<String, f64>,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub star: bool,
    #[serde(default)]
    pub note: String,
    #[serde(default)]
    pub artifacts: Vec<Artifact>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// the convergence check of the problem this run solved, when one was made (optional: runs written before it existed load without it)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub convergence: Option<Convergence>,
    /// how the run's eigensolves were started when a sweep asked for warm starts (productivity plan P5b): which point's modes, how many
    /// restarts, and that the result is not confirmed until it is solved cold (optional: other runs carry none)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reuse: Option<serde_json::Value>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RunPatch {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub star: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tags: Option<Vec<String>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub convergence: Option<Convergence>,
}

#[derive(Serialize, Deserialize)]
struct PatchLine {
    id: String,
    patch: RunPatch,
}

fn append_line(path: &Path, line: &str) -> Result<(), String> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(|e| crate::util::io_err(dir, e))?;
    }
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path).map_err(|e| crate::util::io_err(path, e))?;
    // one write call, newline included: concurrent appenders cannot interleave inside a line
    f.write_all(format!("{line}\n").as_bytes()).map_err(|e| crate::util::io_err(path, e))
}

pub fn append(path: &Path, r: &RunRecord) -> Result<(), String> {
    append_line(path, &serde_json::to_string(r).map_err(|e| e.to_string())?)
}

pub fn append_patch(path: &Path, id: &str, patch: RunPatch) -> Result<(), String> {
    append_line(path, &serde_json::to_string(&PatchLine { id: id.into(), patch }).map_err(|e| e.to_string())?)
}

/// The runs in file order with patches applied, and the number of lines that could not be read (a torn last line, say).
pub fn read(path: &Path) -> (Vec<RunRecord>, usize) {
    let Ok(text) = std::fs::read_to_string(path) else { return (vec![], 0) };
    parse(&text)
}

pub fn parse(text: &str) -> (Vec<RunRecord>, usize) {
    let mut runs: Vec<RunRecord> = vec![];
    let mut bad = 0;
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        if let Ok(p) = serde_json::from_str::<PatchLine>(line) {
            if let Some(r) = runs.iter_mut().find(|r| r.id == p.id) {
                if let Some(s) = p.patch.star {
                    r.star = s;
                }
                if let Some(t) = p.patch.tags {
                    r.tags = t;
                }
                if let Some(n) = p.patch.note {
                    r.note = n;
                }
                if let Some(c) = p.patch.convergence {
                    r.convergence = Some(c);
                }
            }
        } else if let Ok(r) = serde_json::from_str::<RunRecord>(line) {
            runs.push(r);
        } else {
            bad += 1;
        }
    }
    (runs, bad)
}

/// Counts and ranges of `runs` (failed and cancelled runs count but add no ranges).
pub fn digest(runs: &[RunRecord]) -> Digest {
    let mut d = Digest::default();
    let mut ps: BTreeMap<String, Option<Range>> = BTreeMap::new();
    let mut ms: BTreeMap<String, Option<Range>> = BTreeMap::new();
    for r in runs {
        d.count += 1;
        match r.status {
            RunStatus::Ok => d.ok += 1,
            _ => d.failed += 1,
        }
        if d.last.as_deref().is_none_or(|l| r.t.as_str() > l) {
            d.last = Some(r.t.clone());
        }
        if r.star {
            d.starred.push(r.id.clone());
        }
        if r.status != RunStatus::Ok {
            continue;
        }
        for (k, v) in &r.params {
            Range::add(ps.entry(k.clone()).or_default(), *v);
        }
        for (k, v) in &r.metrics {
            Range::add(ms.entry(k.clone()).or_default(), *v);
        }
    }
    d.params = ps.into_iter().filter_map(|(k, v)| v.map(|v| (k, v))).collect();
    d.metrics = ms.into_iter().filter_map(|(k, v)| v.map(|v| (k, v))).collect();
    d
}

/// The next run id for `runs`: `r0007_7hb` (number, then three characters of entropy so two writers cannot collide).
pub fn next_id(runs: &[RunRecord]) -> String {
    let n = runs.iter().filter_map(|r| r.id.strip_prefix('r').and_then(|s| s.split('_').next()).and_then(|s| s.parse::<u32>().ok())).max().unwrap_or(0) + 1;
    format!("r{n:04}_{}", crate::util::run_suffix())
}

// ------------------------------------------------------------------------------------------------------ the filter language

/// A number with an optional unit suffix, in SI: `450 nm`, `1.55um`, `193 THz`, `2e-7`.
pub fn parse_quantity(s: &str) -> Option<f64> {
    let s = s.trim();
    let split = s.char_indices().find(|(i, c)| c.is_alphabetic() && !(matches!(c, 'e' | 'E') && s[i + 1..].starts_with(|n: char| n.is_ascii_digit() || n == '-' || n == '+'))).map(|(i, _)| i).unwrap_or(s.len());
    let (n, unit) = s.split_at(split);
    let v: f64 = n.trim().parse().ok()?;
    let (mul, div) = match unit.trim() {
        "" | "m" | "Hz" | "s" => (1.0, 1.0),
        "mm" => (1.0, 1e3),
        "um" | "µm" | "μm" => (1.0, 1e6),
        "nm" => (1.0, 1e9),
        "pm" => (1.0, 1e12),
        "kHz" => (1e3, 1.0),
        "MHz" => (1e6, 1.0),
        "GHz" => (1e9, 1.0),
        "THz" => (1e12, 1.0),
        "ms" => (1.0, 1e3),
        "us" | "µs" => (1.0, 1e6),
        "ns" => (1.0, 1e9),
        _ => return None,
    };
    Some(v * mul / div)
}

#[derive(Clone, Debug, PartialEq)]
enum Cmp {
    Gt(f64),
    Ge(f64),
    Lt(f64),
    Le(f64),
    Eq(f64),
    Ne(f64),
    Between(f64, f64),
    Is(String),
}

#[derive(Clone, Debug, PartialEq)]
struct Clause {
    field: String,
    cmp: Cmp,
}

/// A parsed filter such as `params.w between 450 nm and 550 nm and metrics.neff_re > 2.4 and star`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Filter {
    clauses: Vec<Clause>,
}

impl Filter {
    /// Fields: `params.NAME`, `metrics.NAME`, `version`, `sec`; `status = ok|failed|cancelled`, `tag = NAME`, `star` and `note` alone.
    pub fn parse(text: &str) -> Result<Filter, String> {
        let toks: Vec<&str> = text.split_whitespace().collect();
        let mut clauses = vec![];
        let mut i = 0;
        // a value is one token, or two when the second is a unit
        let value = |i: &mut usize| -> Result<f64, String> {
            let t = toks.get(*i).ok_or("a value is missing at the end")?;
            *i += 1;
            if t.parse::<f64>().is_err() {
                return parse_quantity(t).ok_or_else(|| format!("'{t}' is not a number"));
            }
            // a plain number: a unit may follow as the next word
            match toks.get(*i).filter(|u| parse_quantity(&format!("1{u}")).is_some() && **u != "and" && **u != "between") {
                Some(u) => {
                    *i += 1;
                    parse_quantity(&format!("{t} {u}")).ok_or_else(|| format!("'{t} {u}' is not a quantity"))
                }
                None => t.parse::<f64>().map_err(|_| format!("'{t}' is not a number")),
            }
        };
        while i < toks.len() {
            let field = toks[i].to_string();
            i += 1;
            if matches!(field.as_str(), "star" | "starred") {
                clauses.push(Clause { field: "star".into(), cmp: Cmp::Is("true".into()) });
            } else {
                let op = *toks.get(i).ok_or_else(|| format!("'{field}' needs a comparison such as > 2.4"))?;
                i += 1;
                let cmp = match op {
                    ">" => Cmp::Gt(value(&mut i)?),
                    ">=" => Cmp::Ge(value(&mut i)?),
                    "<" => Cmp::Lt(value(&mut i)?),
                    "<=" => Cmp::Le(value(&mut i)?),
                    "!=" => Cmp::Ne(value(&mut i)?),
                    "=" | "==" => {
                        if matches!(field.as_str(), "status" | "tag" | "note") {
                            let w = toks.get(i).ok_or("a word is missing at the end")?;
                            i += 1;
                            Cmp::Is(w.to_string())
                        } else {
                            Cmp::Eq(value(&mut i)?)
                        }
                    }
                    "between" => {
                        let a = value(&mut i)?;
                        if toks.get(i) != Some(&"and") {
                            return Err("'between' needs 'A and B'".into());
                        }
                        i += 1;
                        Cmp::Between(a, value(&mut i)?)
                    }
                    o => return Err(format!("'{o}' is not a comparison (use > >= < <= = != between)")),
                };
                clauses.push(Clause { field, cmp });
            }
            match toks.get(i) {
                None => {}
                Some(&"and") => i += 1,
                Some(t) => return Err(format!("expected 'and' but found '{t}'")),
            }
        }
        Ok(Filter { clauses })
    }

    pub fn is_empty(&self) -> bool {
        self.clauses.is_empty()
    }

    pub fn matches(&self, r: &RunRecord) -> bool {
        self.clauses.iter().all(|c| {
            let num = if let Some(k) = c.field.strip_prefix("params.") {
                r.params.get(k).copied()
            } else if let Some(k) = c.field.strip_prefix("metrics.") {
                r.metrics.get(k).copied()
            } else if c.field == "version" {
                Some(r.version as f64)
            } else if c.field == "sec" {
                Some(r.sec)
            } else {
                None
            };
            match (&c.cmp, num) {
                (Cmp::Is(w), _) => match c.field.as_str() {
                    "status" => format!("{:?}", r.status).eq_ignore_ascii_case(w),
                    "tag" => r.tags.iter().any(|t| t.eq_ignore_ascii_case(w)),
                    "note" => r.note.to_lowercase().contains(&w.to_lowercase()),
                    "star" => r.star,
                    _ => false,
                },
                (_, None) => false,
                (Cmp::Gt(x), Some(v)) => v > *x,
                (Cmp::Ge(x), Some(v)) => v >= *x,
                (Cmp::Lt(x), Some(v)) => v < *x,
                (Cmp::Le(x), Some(v)) => v <= *x,
                (Cmp::Eq(x), Some(v)) => (v - x).abs() <= 1e-9 * x.abs().max(1e-300) + 1e-300,
                (Cmp::Ne(x), Some(v)) => (v - x).abs() > 1e-9 * x.abs().max(1e-300) + 1e-300,
                (Cmp::Between(a, b), Some(v)) => v >= a.min(*b) && v <= a.max(*b),
            }
        })
    }
}

/// CSV of runs: id, time, status, version, then every parameter and metric column (the union, sorted).
pub fn to_csv(runs: &[&RunRecord]) -> String {
    let mut pk: Vec<&String> = runs.iter().flat_map(|r| r.params.keys()).collect();
    pk.sort();
    pk.dedup();
    let mut mk: Vec<&String> = runs.iter().flat_map(|r| r.metrics.keys()).collect();
    mk.sort();
    mk.dedup();
    let mut s = String::from("run,time,status,version,star,seconds");
    for k in &pk {
        s += &format!(",params.{k}");
    }
    for k in &mk {
        s += &format!(",metrics.{k}");
    }
    s += "\n";
    for r in runs {
        s += &format!("{},{},{:?},{},{},{}", r.id, r.t, r.status, r.version, r.star, r.sec);
        for k in &pk {
            s += &format!(",{}", r.params.get(*k).map(|v| format!("{v:e}")).unwrap_or_default());
        }
        for k in &mk {
            s += &format!(",{}", r.metrics.get(*k).map(|v| format!("{v:e}")).unwrap_or_default());
        }
        s += "\n";
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(id: &str, w: f64, neff: f64) -> RunRecord {
        RunRecord { id: id.into(), t: "2026-11-01T00:00:00Z".into(), sec: 0.1, params: [("w".to_string(), w)].into(), metrics: [("neff_re".to_string(), neff)].into(), ..Default::default() }
    }

    #[test]
    fn quantities_with_units() {
        assert_eq!(parse_quantity("450 nm"), Some(450e-9));
        assert_eq!(parse_quantity("1.55um"), Some(1.55e-6));
        assert_eq!(parse_quantity("2e-7"), Some(2e-7));
        assert_eq!(parse_quantity("1.5e3Hz"), Some(1500.0));
        assert_eq!(parse_quantity("193 THz"), Some(193e12));
        assert_eq!(parse_quantity("3 bananas"), None);
    }

    #[test]
    fn the_filter_language_selects_runs() {
        let rs = [run("r1", 4.5e-7, 2.38), run("r2", 5.0e-7, 2.41), run("r3", 5.5e-7, 2.45), run("r4", 6.0e-7, 2.5)];
        let pick = |q: &str| -> Vec<String> { let f = Filter::parse(q).unwrap(); rs.iter().filter(|r| f.matches(r)).map(|r| r.id.clone()).collect() };
        assert_eq!(pick("params.w between 450 nm and 550 nm and metrics.neff_re > 2.4"), ["r2", "r3"]);
        assert_eq!(pick("params.w between 4.5e-7 and 5.5e-7"), ["r1", "r2", "r3"]);
        assert_eq!(pick("metrics.neff_re >= 2.5"), ["r4"]);
        assert_eq!(pick("params.w < 500 nm"), ["r1"]);
        assert_eq!(pick("").len(), 4);
        assert!(pick("params.nothing > 1").is_empty());
        assert!(Filter::parse("params.w ?? 3").is_err());
        assert!(Filter::parse("params.w between 1").is_err());
        assert!(Filter::parse("params.w >").is_err());
    }

    #[test]
    fn append_only_lines_patches_and_a_torn_last_line() {
        let dir = std::env::temp_dir().join(format!("nyx-store-runs-{}", crate::util::short_id()));
        let path = dir.join("runs.jsonl");
        let mut a = run("r0001_abc", 5e-7, 2.4);
        a.version = 3;
        append(&path, &a).unwrap();
        append(&path, &run("r0002_def", 5.1e-7, 2.41)).unwrap();
        append_patch(&path, "r0001_abc", RunPatch { star: Some(true), tags: Some(vec!["good".into()]), note: Some("keep".into()), ..Default::default() }).unwrap();
        // a crash in the middle of a line
        std::fs::OpenOptions::new().append(true).open(&path).unwrap().write_all(b"{\"id\":\"r0003_x\",\"t\":\"20").unwrap();
        let (runs, bad) = read(&path);
        assert_eq!(runs.len(), 2);
        assert_eq!(bad, 1);
        assert!(runs[0].star && runs[0].tags == ["good"] && runs[0].note == "keep" && !runs[1].star);
        assert_eq!(next_id(&runs).split('_').next(), Some("r0003"));
        let d = digest(&runs);
        assert_eq!((d.count, d.ok, d.failed), (2, 2, 0));
        assert_eq!(d.starred, ["r0001_abc"]);
        assert_eq!(d.params["w"].n, 2);
        assert!((d.metrics["neff_re"].max - 2.41).abs() < 1e-12);
        assert!(to_csv(&runs.iter().collect::<Vec<_>>()).starts_with("run,time,status,version,star,seconds,params.w,metrics.neff_re\n"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_convergence_verdict_is_a_patch_and_old_lines_load_without_it() {
        let dir = std::env::temp_dir().join(format!("nyx-store-conv-{}", crate::util::short_id()));
        let path = dir.join("runs.jsonl");
        append(&path, &run("r0001_a", 5e-7, 2.4)).unwrap();
        append(&path, &run("r0002_b", 5e-7, 2.4)).unwrap();
        let c = Convergence { factor: 0.6, tolerance: 1e-3, observable: "eigenvalues, matched by index".into(), passed: true, worst_relative: 2.5e-4 };
        append_patch(&path, "r0001_a", RunPatch { convergence: Some(c.clone()), ..Default::default() }).unwrap();
        // a later check of the same problem replaces the earlier verdict; a star patch afterwards leaves it alone
        let failed = Convergence { passed: false, worst_relative: 0.02, ..c.clone() };
        append_patch(&path, "r0001_a", RunPatch { convergence: Some(failed.clone()), ..Default::default() }).unwrap();
        append_patch(&path, "r0001_a", RunPatch { star: Some(true), ..Default::default() }).unwrap();
        let (runs, bad) = read(&path);
        assert_eq!(bad, 0);
        assert_eq!(runs[0].convergence, Some(failed));
        assert!(runs[0].star);
        assert_eq!(runs[1].convergence, None);
        // a record without a verdict does not write the key, so every earlier line and exemplar is unchanged
        assert!(!serde_json::to_string(&runs[1]).unwrap().contains("convergence"));
        let old = r#"{"id":"r1","t":"2026-10-31T00:00:00Z","sec":1.0,"status":"ok","spec":"b3:1","femcore":"0.1.0","params":{},"metrics":{}}"#;
        assert!(parse(old).0[0].convergence.is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
