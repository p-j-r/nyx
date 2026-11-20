// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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

    /// The values of a field by its flat name (`V`, `|E|`, `Ex`; see `flat_fields`).
    pub fn flat(&self, name: &str) -> Option<Vec<f64>> {
        self.flat_fields().into_iter().find(|(n, _, _)| n == name).map(|(_, _, v)| v)
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

impl Field {
    /// A scalar field, one value per dof.
    pub fn scalar(name: impl Into<String>, unit: impl Into<String>, values: Vec<f64>) -> Field {
        Field { name: name.into(), unit: unit.into(), components: vec![], values }
    }

    /// A vector field from one array per component (`comps[j][dof]`), stored interleaved by dof.
    pub fn vector(name: impl Into<String>, unit: impl Into<String>, names: &[&str], comps: &[Vec<f64>]) -> Field {
        assert_eq!(names.len(), comps.len(), "one array per component name");
        let n = comps.first().map_or(0, Vec::len);
        let mut values = Vec::with_capacity(n * comps.len());
        for i in 0..n {
            for c in comps {
                values.push(c[i]);
            }
        }
        Field {
            name: name.into(),
            unit: unit.into(),
            components: names.iter().map(|s| s.to_string()).collect(),
            values,
        }
    }

    /// The array of component `j` of a vector field (the whole array for a scalar).
    pub fn component(&self, j: usize) -> Vec<f64> {
        let m = self.components.len().max(1);
        self.values.chunks(m).map(|p| p[j]).collect()
    }
}

impl Solution {
    /// The fields in the flat form the window and the golden files use: a scalar keeps its name; a vector `E` with components `x`, `z`
    /// becomes `|E|` (the norm), `Ex`, `Ez`, in that order. This is the naming of the pre-restructure results (ADR 41).
    pub fn flat_fields(&self) -> Vec<(String, String, Vec<f64>)> {
        let mut out = vec![];
        for f in &self.fields {
            if f.components.is_empty() {
                out.push((f.name.clone(), f.unit.clone(), f.values.clone()));
            } else {
                let comps: Vec<Vec<f64>> = (0..f.components.len()).map(|j| f.component(j)).collect();
                let n = comps[0].len();
                let norm: Vec<f64> = (0..n).map(|i| comps.iter().map(|c| c[i] * c[i]).sum::<f64>().sqrt()).collect();
                out.push((format!("|{}|", f.name), f.unit.clone(), norm));
                for (name, c) in f.components.iter().zip(comps) {
                    out.push((format!("{}{}", f.name, name), f.unit.clone(), c));
                }
            }
        }
        out
    }
}

/// The per-cell fields of a kernel (`femcore::fields::CellFields`) as plain serialisable data.
pub fn cell_fields_of(cf: &femcore::fields::CellFields) -> Vec<CellField> {
    cf.fields
        .iter()
        .map(|f| CellField {
            name: f.name.clone(),
            unit: f.unit.clone(),
            nloc: cf.nloc,
            re: f.re.clone(),
            im: f.im.clone(),
        })
        .collect()
}

/// The inverse of `cell_fields_of`, for code that still reads `CellFields` (the window until phase 7, the analysis until phase 5).
pub fn cell_fields_back(order: usize, fields: &[CellField]) -> femcore::fields::CellFields {
    let nloc = fields.first().map_or(0, |f| f.nloc);
    femcore::fields::CellFields {
        order,
        nloc,
        fields: fields
            .iter()
            .map(|f| femcore::fields::CellField {
                name: f.name.clone(),
                unit: f.unit.clone(),
                re: f.re.clone(),
                im: f.im.clone(),
            })
            .collect(),
    }
}

impl femcore::adapt::PlotFields for SolveResult {
    fn plot_fields(&self) -> Vec<Vec<(String, Vec<f64>)>> {
        self.solutions.iter().map(|s| s.flat_fields().into_iter().map(|(n, _, v)| (n, v)).collect()).collect()
    }
}
