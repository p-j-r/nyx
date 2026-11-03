//! What a solve returns: plain, serialisable data (flat arrays), so a result can be stored, sent to another process or compared with a
//! golden file. It does not carry an analysis: analyses are functions of a result (`nyx-analysis`).

use serde::{Deserialize, Serialize};

/// A field with one or more components, one value per dof (P1: vertices; P2: vertices then edge midpoints).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Field {
    pub name: String,
    pub unit: String,
    /// wafer-axis names of the components (`x`, `y`, `z`); empty for a scalar
    pub components: Vec<String>,
    /// `n_dofs` values for a scalar, `n_dofs * components.len()` for a vector, interleaved by dof
    pub values: Vec<f64>,
}

/// A field before it is averaged over the cells at a node, one value per cell node (`nloc` per cell).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CellField {
    pub name: String,
    pub unit: String,
    pub nloc: usize,
    pub re: Vec<f64>,
    /// empty for a real field
    pub im: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum Scalar {
    Num(f64),
    Text(String),
}

/// One solution of a study: a time step, an eigenmode, a frequency.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Solution {
    pub label: String,
    pub fields: Vec<Field>,
    /// the real value that orders the solutions: time (s), frequency (Hz) or mode number
    pub value: Option<f64>,
    /// a complex eigenvalue in the study's own unit (effective index, complex frequency), as (re, im)
    pub eigenvalue: Option<[f64; 2]>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cell_fields: Vec<CellField>,
}

impl Solution {
    pub fn field(&self, name: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.name == name)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct SolveResult {
    /// the code of the physics that produced it and the study kind
    pub physics: String,
    pub study: String,
    pub order: u8,
    /// 2 or 3
    pub dim: u8,
    pub dofs: usize,
    pub seconds: f64,
    pub solutions: Vec<Solution>,
    pub scalars: Vec<(String, Scalar)>,
    /// the plain-language stage log of the run
    pub log: Vec<String>,
}

impl SolveResult {
    /// The first numeric scalar whose name starts with `name`.
    pub fn scalar(&self, name: &str) -> Option<f64> {
        self.scalars.iter().find(|(n, _)| n.starts_with(name)).and_then(|(_, v)| match v {
            Scalar::Num(x) => Some(*x),
            Scalar::Text(_) => None,
        })
    }
}
