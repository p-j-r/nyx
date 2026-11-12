//! The problem envelope: the one shape every physics is given. Each physics owns the typed form of `study.params`, `bcs[].params` and
//! `params` and validates them against its descriptors (`validate`); nothing physics-specific lives here.
//!
//! Regions and boundaries are named by `geom` names and selectors (ADR 2); the engine compiles them to labels through a
//! `LabelResolver` when it solves. The `format` key (first value 1) is written in every problem and checked by `validate_spec` (ADR 63).

use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;

/// What the material of a region is.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum MaterialAssign {
    /// the name of a material of the project or of a library (a TOML file, ADR 29)
    Named(String),
    /// a material of the project or of the library placed in the wafer frame: the Bunge angles (degrees) of its crystal axes
    Oriented { material: String, euler_deg: [f64; 3] },
    /// scalar properties given inline by descriptor name (`E`, `nu`, `rho`, ...)
    Props(BTreeMap<String, f64>),
}

/// The study: a `kind` the physics lists in `studies()` and its parameters.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StudySpec {
    pub kind: String,
    #[serde(default, skip_serializing_if = "Value::is_null")]
    pub params: Value,
}

/// A boundary or domain condition: `on` is a selector over the geometry (`outer & plane(x = 0)`) or a region name.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BcSpec {
    pub on: String,
    pub kind: String,
    #[serde(default, skip_serializing_if = "Value::is_null")]
    pub params: Value,
}

/// The file format of a problem (ADR 63). Written in every problem; a problem of another format is refused, there is no migration.
pub const FORMAT: u32 = 1;

fn current_format() -> u32 {
    FORMAT
}

/// A problem for one physics.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProblemSpec {
    /// the file format, `1` (ADR 63)
    #[serde(default = "current_format")]
    pub format: u32,
    /// the three-letter code of the physics (`dfs`, `sms`, `ems`)
    pub physics: String,
    /// element order (1 and 2 today; the physics says which it supports)
    pub order: u8,
    pub study: StudySpec,
    /// material by region name
    #[serde(default)]
    pub materials: BTreeMap<String, MaterialAssign>,
    #[serde(default)]
    pub bcs: Vec<BcSpec>,
    /// physics-level options (`ParamDesc`s of `Physics::params`)
    #[serde(default, skip_serializing_if = "Value::is_null")]
    pub params: Value,
}

impl ProblemSpec {
    pub fn new(physics: &str, order: u8, study_kind: &str) -> ProblemSpec {
        ProblemSpec {
            format: FORMAT,
            physics: physics.into(),
            order,
            study: StudySpec { kind: study_kind.into(), params: Value::Null },
            materials: BTreeMap::new(),
            bcs: vec![],
            params: Value::Null,
        }
    }
}

/// The same JSON with every object's keys in sorted order. Contract exemplars are written and compared in this form because `serde_json`
/// keeps insertion order when any crate of the build enables `preserve_order` (the app does) and sorts otherwise.
pub fn canonical_json(v: &Value) -> Value {
    match v {
        Value::Object(m) => {
            let mut keys: Vec<&String> = m.keys().collect();
            keys.sort();
            let mut out = serde_json::Map::new();
            for k in keys {
                out.insert(k.clone(), canonical_json(&m[k]));
            }
            Value::Object(out)
        }
        Value::Array(a) => Value::Array(a.iter().map(canonical_json).collect()),
        other => other.clone(),
    }
}
