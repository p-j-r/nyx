// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Descriptors: what a physics tells the rest of nyx about itself. The window builds forms from them, the MCP server builds JSON
//! Schemas from them, the mesher reads the sizing hint, and the engine validates a problem against them. All plain serde data.

use serde::{Deserialize, Serialize};

/// Which meshes a study runs on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Dims {
    /// a 2D cross-section (x-z plane of the wafer frame)
    Two,
    /// a 3D mesh
    Three,
    Both,
}

/// The type and range of one parameter.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ParamKind {
    Number {
        default: f64,
        min: Option<f64>,
        max: Option<f64>,
    },
    Integer {
        default: i64,
        min: Option<i64>,
        max: Option<i64>,
    },
    Bool {
        default: bool,
    },
    Choice {
        options: Vec<String>,
        default: String,
    },
    /// a fixed-length list of numbers (a force vector)
    Vector {
        len: usize,
        default: Vec<f64>,
    },
    /// free text (an expression in x, y, z)
    Text {
        default: String,
    },
    /// a list of region names (the regions that take no part in this physics)
    Names {
        default: Vec<String>,
    },
    /// a variable-length list of numbers (e.g. an explicit frequency list)
    Numbers {
        default: Vec<f64>,
    },
}

/// One named parameter of a study, a boundary condition or the physics.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ParamDesc {
    pub name: String,
    pub label: String,
    /// SI unit, empty for dimensionless
    pub unit: String,
    pub kind: ParamKind,
    /// a required parameter has no usable default: it must be given
    pub required: bool,
    pub doc: String,
}

impl ParamDesc {
    fn make(name: &str, label: &str, unit: &str, kind: ParamKind, doc: &str) -> ParamDesc {
        ParamDesc { name: name.into(), label: label.into(), unit: unit.into(), kind, required: false, doc: doc.into() }
    }
    pub fn number(name: &str, label: &str, unit: &str, default: f64, doc: &str) -> ParamDesc {
        Self::make(name, label, unit, ParamKind::Number { default, min: None, max: None }, doc)
    }
    pub fn integer(name: &str, label: &str, default: i64, min: i64, doc: &str) -> ParamDesc {
        Self::make(name, label, "", ParamKind::Integer { default, min: Some(min), max: None }, doc)
    }
    pub fn boolean(name: &str, label: &str, default: bool, doc: &str) -> ParamDesc {
        Self::make(name, label, "", ParamKind::Bool { default }, doc)
    }
    pub fn choice(name: &str, label: &str, options: &[&str], default: &str, doc: &str) -> ParamDesc {
        Self::make(
            name,
            label,
            "",
            ParamKind::Choice { options: options.iter().map(|s| s.to_string()).collect(), default: default.into() },
            doc,
        )
    }
    pub fn vector(name: &str, label: &str, unit: &str, len: usize, doc: &str) -> ParamDesc {
        Self::make(name, label, unit, ParamKind::Vector { len, default: vec![0.0; len] }, doc)
    }
    /// Sets the default vector for a `ParamKind::Vector` or `ParamKind::Numbers` descriptor.
    pub fn default_vector(mut self, vals: Vec<f64>) -> ParamDesc {
        match &mut self.kind {
            ParamKind::Vector { len, default } => {
                assert_eq!(vals.len(), *len, "default vector length mismatch");
                *default = vals;
            }
            ParamKind::Numbers { default } => {
                *default = vals;
            }
            _ => {}
        }
        self
    }
    pub fn text(name: &str, label: &str, default: &str, doc: &str) -> ParamDesc {
        Self::make(name, label, "", ParamKind::Text { default: default.into() }, doc)
    }
    pub fn names(name: &str, label: &str, doc: &str) -> ParamDesc {
        Self::make(name, label, "", ParamKind::Names { default: vec![] }, doc)
    }
    pub fn numbers(name: &str, label: &str, unit: &str, doc: &str) -> ParamDesc {
        Self::make(name, label, unit, ParamKind::Numbers { default: vec![] }, doc)
    }
    /// Marks the parameter as one that must be given.
    pub fn required(mut self) -> ParamDesc {
        self.required = true;
        self
    }
    /// Sets the allowed range of a number or an integer.
    pub fn range(mut self, lo: Option<f64>, hi: Option<f64>) -> ParamDesc {
        match &mut self.kind {
            ParamKind::Number { min, max, .. } => {
                *min = lo;
                *max = hi;
            }
            ParamKind::Integer { min, max, .. } => {
                *min = lo.map(|v| v as i64);
                *max = hi.map(|v| v as i64);
            }
            _ => {}
        }
        self
    }
}

/// A kind of study (stationary, eigenfrequency, modes, ...).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StudyDesc {
    /// the key written in a problem: `study.kind`
    pub kind: String,
    pub name: String,
    pub doc: String,
    pub dims: Dims,
    pub params: Vec<ParamDesc>,
}

/// Whether a boundary condition sits on a boundary (facet) selection or fills a domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Applies {
    Boundary,
    Domain,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BcDesc {
    pub kind: String,
    pub name: String,
    pub doc: String,
    pub applies: Applies,
    pub params: Vec<ParamDesc>,
}

/// How a material property is given.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PropKind {
    Scalar,
    /// a symmetric tensor: one value, three principal values or six entries (phase 6)
    Tensor,
    /// a 3x6 coupling matrix (e.g. piezoelectric stress-charge tensor e_iJ)
    Matrix3x6,
    /// a model evaluated at the study wavelength (Sellmeier, tabulated, ...; phase 6)
    Dispersive,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PropDesc {
    /// the name used in a material (`E`, `nu`, `rho`, `eps_r`, ...)
    pub name: String,
    pub label: String,
    pub unit: String,
    pub kind: PropKind,
    /// the study kinds that need it; empty means every study of the physics
    pub needed_by: Vec<String>,
    /// the value used when a material does not give one; `None` means it must be given
    pub default: Option<f64>,
}

impl PropDesc {
    pub fn new(
        name: &str,
        label: &str,
        unit: &str,
        kind: PropKind,
        needed_by: &[&str],
        default: Option<f64>,
    ) -> PropDesc {
        PropDesc {
            name: name.into(),
            label: label.into(),
            unit: unit.into(),
            kind,
            needed_by: needed_by.iter().map(|s| s.to_string()).collect(),
            default,
        }
    }
}

/// A field a study produces. Components are named in wafer axes (`x` lateral, `y` propagation, `z` normal); empty is a scalar.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FieldDesc {
    pub name: String,
    pub unit: String,
    pub components: Vec<String>,
    pub doc: String,
}

impl FieldDesc {
    pub fn new(name: &str, unit: &str, components: &[&str], doc: &str) -> FieldDesc {
        FieldDesc {
            name: name.into(),
            unit: unit.into(),
            components: components.iter().map(|s| s.to_string()).collect(),
            doc: doc.into(),
        }
    }
}

/// One end of a one-way coupling between two physics (restructure plan section 2). A physics *exports* fields under a `role` and
/// *imports* a role it needs; an import is satisfied by any export with the same role (`temperature`, `heat_source`, ...).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FieldLink {
    /// the field, qualified with the physics code (`dfs.T`)
    pub field: String,
    /// the exchanged quantity
    pub role: String,
}

impl FieldLink {
    pub fn new(field: &str, role: &str) -> FieldLink {
        FieldLink { field: field.into(), role: role.into() }
    }
}

/// What the mesher needs from a physics to size the elements before any solution exists. Defined by the mesher (phase 2), so that the
/// meshing code knows no physics; the physics crates fill it in (`Physics::sizing`).
pub use femcore::meshing::Sizing;
