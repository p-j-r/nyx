//! Search over the manifests: plain words with a synonym table, and numbers with units matched against card dimensions, wavelengths
//! and the parameter ranges of runs. No index and no engine: a few thousand manifests scan in milliseconds.
//!
//! "300 nm sapphire plasmon" becomes the numeric term 3e-7 m and the word terms "sapphire" (also "silicon on sapphire", "Al2O3") and "plasmon".
//! A workflow is shown when at least half of the terms match; the score is the share matched. The query box says how it was understood.

use crate::model::{Experiment, Workflow};
use crate::runs::{parse_quantity, RunRecord};
use crate::workspace::Entry;

/// Tolerance for a numeric term: 2 %.
pub const TOL: f64 = 0.02;

const STOP: &[&str] = &["a", "an", "the", "of", "on", "in", "and", "for", "with", "to", "from", "at", "by"];

/// A word and what else it should find. `sos` is read as silicon on sapphire (SOI, silicon on insulator, is separate).
const SYNONYMS: &[(&str, &[&str])] = &[
    ("sos", &["silicon on sapphire", "sapphire", "al2o3"]),
    ("sapphire", &["al2o3", "alumina", "silicon on sapphire"]),
    ("al2o3", &["sapphire", "alumina"]),
    ("soi", &["silicon on insulator", "sio2", "oxide"]),
    ("silica", &["sio2", "oxide", "glass"]),
    ("sio2", &["silica", "oxide"]),
    ("silicon", &["si"]),
    ("nitride", &["si3n4", "silicon nitride"]),
    ("sin", &["si3n4", "silicon nitride"]),
    ("waveguide", &["rib", "strip", "wire", "guide"]),
    ("rib", &["waveguide"]),
    ("neff", &["effective index", "modes"]),
    ("phonon", &["solid mechanics", "elastic"]),
    ("thermal", &["heat", "temperature"]),
    ("optical", &["em waves", "electromagnetic"]),
];

#[derive(Clone, Debug, PartialEq)]
pub enum Term {
    Word(String),
    /// a quantity in SI with the text the person typed
    Number(f64, String),
}

#[derive(Clone, Debug, Default)]
pub struct Query {
    pub terms: Vec<Term>,
}

impl Query {
    pub fn parse(text: &str) -> Query {
        let toks: Vec<String> = text.split_whitespace().map(|t| t.trim_matches(|c: char| c == ',' || c == '"' || c == '\u{201c}' || c == '\u{201d}').to_lowercase()).filter(|t| !t.is_empty()).collect();
        let mut terms = vec![];
        let mut i = 0;
        while i < toks.len() {
            let t = &toks[i];
            // "300 nm" as two tokens, or "300nm" as one
            if let Some(v) = parse_quantity(t).filter(|_| t.chars().any(|c| c.is_ascii_digit())) {
                if t.chars().all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | '-' | '+')) {
                    if let Some(u) = toks.get(i + 1).filter(|u| parse_quantity(&format!("1{u}")).is_some()) {
                        terms.push(Term::Number(parse_quantity(&format!("{t} {u}")).unwrap_or(v), format!("{t} {u}")));
                        i += 2;
                        continue;
                    }
                    terms.push(Term::Number(v, t.clone()));
                } else {
                    terms.push(Term::Number(v, t.clone()));
                }
            } else if !STOP.contains(&t.as_str()) {
                terms.push(Term::Word(t.clone()));
            }
            i += 1;
        }
        Query { terms }
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    /// How the query was understood, as short phrases for the "Understood as" panel.
    pub fn understood(&self) -> Vec<String> {
        self.terms
            .iter()
            .map(|t| match t {
                Term::Number(v, s) => format!("{} {s}", if *v < 1e-3 { "size" } else { "value" }),
                Term::Word(w) => match SYNONYMS.iter().find(|(k, _)| k == w) {
                    Some((_, syn)) => format!("{w} (also {})", syn.join(", ")),
                    None => format!("word {w}"),
                },
            })
            .collect()
    }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() <= TOL * b.abs().max(1e-300)
}

fn haystack(wf: &Workflow) -> String {
    let mut h = format!("{} {} {}", wf.title, wf.tags.join(" "), format!("{:?}", wf.status));
    for e in &wf.experiments {
        h += &format!(" {} {} {} {} {} {}", e.title, e.card.summary, e.card.keywords.join(" "), e.card.materials.join(" "), e.card.physics, e.card.geometry.family);
    }
    h.to_lowercase()
}

fn number_matches(e: &Experiment, v: f64) -> Option<String> {
    for d in &e.card.geometry.dims {
        if near(d.m, v) {
            return Some(format!("{} {} {:.0} nm", d.of, d.name, d.m * 1e9));
        }
    }
    if let Some(w) = e.card.wavelength_m {
        if near(w, v) {
            return Some(format!("wavelength {:.0} nm", w * 1e9));
        }
    }
    for (k, r) in &e.runs.params {
        if (r.min <= v * (1.0 + TOL) && r.max >= v * (1.0 - TOL)) || near(r.min, v) || near(r.max, v) {
            return Some(format!("runs sweep {k} over {v:.3e}"));
        }
    }
    None
}

#[derive(Clone, Debug)]
pub struct Hit {
    /// index into the entries searched
    pub index: usize,
    pub score: f64,
    /// what matched, for the line under the title
    pub reasons: Vec<String>,
}

/// Workflows matching `q`, best first (ties: newest touched first). An empty query matches nothing.
pub fn search(q: &Query, entries: &[Entry]) -> Vec<Hit> {
    if q.is_empty() {
        return vec![];
    }
    let mut hits = vec![];
    for (index, e) in entries.iter().enumerate() {
        let h = haystack(&e.wf);
        let mut matched = 0;
        let mut reasons = vec![];
        for t in &q.terms {
            match t {
                Term::Word(w) => {
                    let syn = SYNONYMS.iter().find(|(k, _)| k == w).map(|(_, s)| *s).unwrap_or(&[]);
                    if h.contains(w.as_str()) {
                        matched += 1;
                        reasons.push(w.clone());
                    } else if let Some(s) = syn.iter().find(|s| h.contains(**s)) {
                        matched += 1;
                        reasons.push(format!("{w} ({s})"));
                    }
                }
                Term::Number(v, s) => {
                    if let Some(r) = e.wf.experiments.iter().find_map(|x| number_matches(x, *v)) {
                        matched += 1;
                        reasons.push(format!("{s}: {r}"));
                    }
                }
            }
        }
        let score = matched as f64 / q.terms.len() as f64;
        if matched > 0 && score >= 0.5 {
            hits.push(Hit { index, score, reasons });
        }
    }
    hits.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap().then(entries[b.index].wf.touched.cmp(&entries[a.index].wf.touched)));
    hits
}

#[derive(Clone, Debug)]
pub struct RunHit {
    pub index: usize,
    pub experiment: String,
    pub run: String,
    pub what: String,
}

/// Runs whose parameters match a numeric term of `q`, in the first `limit` workflows of `hits` (reads their `runs.jsonl`).
pub fn run_hits(q: &Query, entries: &[Entry], hits: &[Hit], limit: usize) -> Vec<RunHit> {
    let nums: Vec<(f64, &String)> = q.terms.iter().filter_map(|t| if let Term::Number(v, s) = t { Some((*v, s)) } else { None }).collect();
    let mut out = vec![];
    if nums.is_empty() {
        return out;
    }
    for h in hits.iter().take(limit) {
        let e = &entries[h.index];
        for x in &e.wf.experiments {
            let (rs, _) = crate::runs::read(&e.dir.join(&x.id).join("runs.jsonl"));
            for r in rs.iter().take(2000) {
                if let Some((k, v)) = matching_param(r, &nums) {
                    out.push(RunHit { index: h.index, experiment: x.id.clone(), run: r.id.clone(), what: format!("{k} {:.0} nm", v * 1e9) });
                    if out.len() >= 12 {
                        return out;
                    }
                }
            }
        }
    }
    out
}

fn matching_param<'a>(r: &'a RunRecord, nums: &[(f64, &String)]) -> Option<(&'a String, f64)> {
    r.params.iter().find(|(_, v)| nums.iter().any(|(n, _)| near(**v, *n))).map(|(k, v)| (k, *v))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::*;
    use std::path::PathBuf;

    fn entry(title: &str, keywords: &[&str], dims: &[(&str, f64)], tags: &[&str]) -> Entry {
        let mut wf = Workflow::new(title, "gui");
        wf.tags = tags.iter().map(|s| s.to_string()).collect();
        let card = Card {
            physics: "Emw".into(),
            keywords: keywords.iter().map(|s| s.to_string()).collect(),
            geometry: GeoCard { family: "2D shapes".into(), summary: "frame".into(), dims: dims.iter().map(|(n, m)| Dim { name: n.to_string(), of: "domain 1".into(), m: *m }).collect() },
            summary: "frame on sapphire".into(),
            ..Default::default()
        };
        wf.experiments.push(Experiment { id: "e1".into(), title: "frame".into(), parent: None, links: vec![], spec_hash: String::new(), card, runs: Digest::default(), model: None });
        Entry { dir: PathBuf::new(), wf }
    }

    #[test]
    fn numbers_with_units_and_words_with_synonyms() {
        let q = Query::parse("plasmon scattering on a 300 nm SOS frame");
        assert_eq!(q.terms.len(), 5, "{:?}", q.terms);
        assert!(q.terms.contains(&Term::Number(3e-7, "300 nm".into())));
        assert!(Query::parse("300nm").terms.contains(&Term::Number(3e-7, "300nm".into())));
        assert_eq!(Query::parse("1550").terms, [Term::Number(1550.0, "1550".into())]);
        let es = vec![
            entry("Plasmon scattering, Ag disk on SOS frame", &["silicon", "sapphire"], &[("width", 3.0e-7)], &["plasmon"]),
            entry("Nanoantenna gap sweep", &["sapphire"], &[("gap", 3.05e-7)], &["plasmon"]),
            entry("Slot waveguide", &["silicon", "silica"], &[("gap", 1e-7)], &[]),
        ];
        let hits = search(&q, &es);
        assert_eq!(hits.iter().map(|h| h.index).collect::<Vec<_>>(), [0, 1], "the slot waveguide matches nothing");
        assert!(hits[0].score > hits[1].score);
        assert!(hits[0].reasons.iter().any(|r| r.contains("300 nm")));
        // 305 nm is within 2 % of 300 nm? 1.7 %: yes
        assert!(hits[1].reasons.iter().any(|r| r.contains("gap")));
        assert!(search(&Query::parse("sos"), &es).iter().any(|h| h.index == 0), "SOS finds silicon on sapphire through the table");
        assert!(search(&Query::parse("soi"), &es).is_empty() || !search(&Query::parse("soi"), &es).iter().any(|h| h.index == 0 || h.index == 1), "SOI is not SOS");
        assert!(search(&Query::parse(""), &es).is_empty());
        let u = q.understood().join(" | ");
        assert!(u.contains("size 300 nm") && u.contains("sos (also silicon on sapphire"), "{u}");
    }
}
