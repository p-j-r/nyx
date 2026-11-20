// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The data model: the workflow manifest (`workflow.json`), its experiments, the auto-generated concept card, the run digest and
//! the share levels. Everything is plain serde and ordered (`BTreeMap`) so files are byte-for-byte reproducible.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const FORMAT: u32 = 1;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Status {
    #[default]
    Active,
    Done,
}

/// What a `.nyx` file carries (Share screen). `Recipe` is the default and the recommended one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ShareLevel {
    /// geometry, materials, mesh and solver settings, run table, versions, notes, tags
    #[default]
    Recipe,
    /// the recipe plus the heavy result files (fields)
    Results,
    /// the recipe, results and mesh files
    Meshes,
}

impl ShareLevel {
    pub const ALL: [ShareLevel; 3] = [ShareLevel::Recipe, ShareLevel::Results, ShareLevel::Meshes];
    pub fn name(self) -> &'static str {
        match self {
            ShareLevel::Recipe => "Recipe only",
            ShareLevel::Results => "Recipe and results",
            ShareLevel::Meshes => "With meshes",
        }
    }
}

/// A named dimension read from the experiment's spec, in metres.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Dim {
    pub name: String,
    pub of: String,
    pub m: f64,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct GeoCard {
    pub family: String,
    pub summary: String,
    pub dims: Vec<Dim>,
}

/// What search sees of an experiment. Generated from its spec (`card::from_spec`), never typed.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Card {
    pub physics: String,
    pub study: String,
    pub dim: u8,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wavelength_m: Option<f64>,
    pub geometry: GeoCard,
    pub materials: Vec<String>,
    pub keywords: Vec<String>,
    /// one line, e.g. "2D shapes, 3 parts, EM waves modes, 1550 nm"
    pub summary: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Range {
    pub min: f64,
    pub max: f64,
    pub n: usize,
}

impl Range {
    pub fn add(r: &mut Option<Range>, v: f64) {
        match r {
            Some(r) => {
                r.min = r.min.min(v);
                r.max = r.max.max(v);
                r.n += 1;
            }
            None => *r = Some(Range { min: v, max: v, n: 1 }),
        }
    }
}

/// Counts and ranges of an experiment's runs; derived from `runs.jsonl` and rebuilt when stale.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Digest {
    pub count: usize,
    pub ok: usize,
    pub failed: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last: Option<String>,
    pub params: BTreeMap<String, Range>,
    pub metrics: BTreeMap<String, Range>,
    pub starred: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Link {
    pub to: String,
    pub label: String,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Experiment {
    pub id: String,
    pub title: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub links: Vec<Link>,
    /// hash of `experiment.nyx.json`
    pub spec_hash: String,
    pub card: Card,
    pub runs: Digest,
    /// a recorded surrogate model for this experiment's runs; no model code exists yet, so this stays `None`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub model: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Workflow {
    pub nyx: u32,
    pub kind: String,
    pub id: String,
    pub title: String,
    pub status: Status,
    pub pinned: bool,
    /// "gui", "batch", "ai", "imported"
    pub origin: String,
    pub created: String,
    /// last edit, run or note; merely viewing never changes it
    pub touched: String,
    pub tags: Vec<String>,
    #[serde(default)]
    pub share: ShareLevel,
    /// the current version number of the set-up (`history.jsonl`)
    #[serde(default)]
    pub version: u32,
    /// shared scripts stay off until the person approves them
    #[serde(default)]
    pub scripts_approved: bool,
    pub experiments: Vec<Experiment>,
}

impl Workflow {
    pub fn new(title: &str, origin: &str) -> Workflow {
        let now = crate::util::iso(crate::util::now_secs());
        Workflow {
            nyx: FORMAT,
            kind: "workflow".into(),
            id: crate::util::short_id(),
            title: title.to_string(),
            status: Status::Active,
            pinned: false,
            origin: origin.into(),
            created: now.clone(),
            touched: now,
            tags: vec![],
            share: ShareLevel::Recipe,
            version: 0,
            scripts_approved: false,
            experiments: vec![],
        }
    }

    pub fn experiment(&self, id: &str) -> Option<&Experiment> {
        self.experiments.iter().find(|e| e.id == id)
    }

    pub fn experiment_mut(&mut self, id: &str) -> Option<&mut Experiment> {
        self.experiments.iter_mut().find(|e| e.id == id)
    }

    pub fn touched_secs(&self) -> i64 {
        crate::util::parse_iso(&self.touched).unwrap_or(0)
    }

    pub fn run_count(&self) -> usize {
        self.experiments.iter().map(|e| e.runs.count).sum()
    }

    /// One line shown on dashboard cards: the card summary of the last experiment.
    pub fn summary(&self) -> String {
        self.experiments.last().map(|e| e.card.summary.clone()).unwrap_or_default()
    }
}
