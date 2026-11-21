// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Reading a TOML material file: a `Section` wraps one table, remembers which keys were read, converts quantities with their units
//! to SI and reports every key nobody read ("unknown key", with a suggestion). Errors name the material and the key path.

use super::units::{self, Dim};
use std::collections::BTreeSet;
use std::fmt;

/// A problem found in a material file or description, in plain words.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaterialError {
    pub material: String,
    /// where in the file, for example `elastic.c44`; empty for the whole file
    pub path: String,
    pub message: String,
}

impl fmt::Display for MaterialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "material '{}'", self.material)?;
        if !self.path.is_empty() {
            write!(f, ", [{}]", self.path)?;
        }
        write!(f, ": {}", self.message)
    }
}

impl std::error::Error for MaterialError {}

/// A quantity tree after unit conversion: numbers in SI, lists, or the uniaxial pair.
#[derive(Clone, Debug, PartialEq)]
pub enum Tree {
    Num(f64),
    List(Vec<Tree>),
    Uniaxial(f64, f64),
}

pub struct Section<'a> {
    pub material: String,
    pub path: String,
    table: &'a toml::Table,
    used: BTreeSet<String>,
}

pub(crate) fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for i in 1..=a.len() {
        let mut cur = vec![i];
        for j in 1..=b.len() {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            cur.push((prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + cost));
        }
        prev = cur;
    }
    prev[b.len()]
}

impl<'a> Section<'a> {
    pub fn new(material: &str, path: &str, table: &'a toml::Table) -> Self {
        Section { material: material.to_string(), path: path.to_string(), table, used: BTreeSet::new() }
    }

    pub fn err(&self, key: &str, message: impl Into<String>) -> MaterialError {
        let path = if key.is_empty() {
            self.path.clone()
        } else if self.path.is_empty() {
            key.to_string()
        } else {
            format!("{}.{}", self.path, key)
        };
        MaterialError { material: self.material.clone(), path, message: message.into() }
    }

    pub fn has(&self, key: &str) -> bool {
        self.table.contains_key(key)
    }

    pub fn keys(&self) -> Vec<String> {
        self.table.keys().cloned().collect()
    }

    pub fn value(&mut self, key: &str) -> Option<&'a toml::Value> {
        let v = self.table.get(key);
        if v.is_some() {
            self.used.insert(key.to_string());
        }
        v
    }

    pub fn string(&mut self, key: &str) -> Result<Option<String>, MaterialError> {
        match self.value(key) {
            None => Ok(None),
            Some(toml::Value::String(s)) => Ok(Some(s.clone())),
            Some(_) => Err(self.err(key, "must be text in quotes")),
        }
    }

    pub fn bool(&mut self, key: &str) -> Result<Option<bool>, MaterialError> {
        match self.value(key) {
            None => Ok(None),
            Some(toml::Value::Boolean(b)) => Ok(Some(*b)),
            Some(_) => Err(self.err(key, "must be true or false")),
        }
    }

    fn num_of(v: &toml::Value) -> Option<f64> {
        match v {
            toml::Value::Float(f) => Some(*f),
            toml::Value::Integer(i) => Some(*i as f64),
            _ => None,
        }
    }

    /// A plain number (no unit); a string such as "0.17" is accepted too.
    pub fn number(&mut self, key: &str) -> Result<Option<f64>, MaterialError> {
        match self.value(key) {
            None => Ok(None),
            Some(v) => self.number_of(key, v).map(Some),
        }
    }

    pub fn number_of(&self, key: &str, v: &toml::Value) -> Result<f64, MaterialError> {
        let x = match v {
            toml::Value::String(s) => s.trim().parse::<f64>().ok(),
            other => Self::num_of(other),
        };
        match x {
            Some(x) if x.is_finite() => Ok(x),
            _ => Err(self.err(key, "must be a finite number")),
        }
    }

    pub fn numbers(&mut self, key: &str) -> Result<Option<Vec<f64>>, MaterialError> {
        match self.value(key) {
            None => Ok(None),
            Some(toml::Value::Array(a)) => a.iter().map(|v| self.number_of(key, v)).collect::<Result<Vec<_>, _>>().map(Some),
            Some(_) => Err(self.err(key, "must be a list of numbers, for example [1.0, 2.0]")),
        }
    }

    /// The unit declared next to `key` as `<key>_unit`, or the default of the table.
    fn unit_for(&mut self, key: &str, default_unit: Option<&str>) -> Result<Option<String>, MaterialError> {
        let k = format!("{key}_unit");
        if let Some(u) = self.string(&k)? {
            return Ok(Some(u));
        }
        Ok(default_unit.map(str::to_string))
    }

    fn convert(&self, key: &str, x: f64, unit: &str, dim: Dim) -> Result<f64, MaterialError> {
        units::to_si(x, unit, dim).map_err(|m| self.err(key, m))
    }

    fn tree(&self, key: &str, v: &toml::Value, dim: Dim, unit: Option<&str>) -> Result<Tree, MaterialError> {
        match v {
            toml::Value::String(s) => {
                let (x, u) = units::split_quantity(s).map_err(|m| self.err(key, m))?;
                if u.is_empty() {
                    return self.bare(key, x, dim, unit);
                }
                Ok(Tree::Num(self.convert(key, x, &u, dim)?))
            }
            toml::Value::Integer(_) | toml::Value::Float(_) => self.bare(key, Self::num_of(v).unwrap(), dim, unit),
            toml::Value::Array(a) => Ok(Tree::List(a.iter().map(|e| self.tree(key, e, dim, unit)).collect::<Result<_, _>>()?)),
            toml::Value::Table(t) => {
                if let Some(val) = t.get("value") {
                    let u = t.get("unit").and_then(|u| u.as_str()).or(unit);
                    for k in t.keys() {
                        if k != "value" && k != "unit" {
                            return Err(self.err(key, format!("the inline table has an unknown entry '{k}'; use {{ value = ..., unit = \"...\" }}")));
                        }
                    }
                    return self.tree(key, val, dim, u);
                }
                if t.contains_key("perp") && t.contains_key("par") && t.len() == 2 {
                    let one = |k: &str| -> Result<f64, MaterialError> {
                        match self.tree(key, &t[k], dim, unit)? {
                            Tree::Num(x) => Ok(x),
                            _ => Err(self.err(key, format!("'{k}' must be a single value"))),
                        }
                    };
                    return Ok(Tree::Uniaxial(one("perp")?, one("par")?));
                }
                Err(self.err(key, "a table here must be { value = ..., unit = \"...\" } or { perp = ..., par = ... }"))
            }
            _ => Err(self.err(key, "must be a number, a quantity such as \"496.8 GPa\", or a list of them")),
        }
    }

    fn bare(&self, key: &str, x: f64, dim: Dim, unit: Option<&str>) -> Result<Tree, MaterialError> {
        if !x.is_finite() {
            return Err(self.err(key, "must be a finite number"));
        }
        if dim == Dim::Dimensionless {
            return Ok(Tree::Num(x));
        }
        match unit {
            Some(u) => Ok(Tree::Num(self.convert(key, x, u, dim)?)),
            None => Err(self.err(
                key,
                format!("has no unit; write it with the unit, for example \"{x} {}\", or declare {key}_unit = \"...\" (units are never assumed)", example_unit(dim)),
            )),
        }
    }

    /// A quantity tree (number, list, or uniaxial pair) converted to SI.
    pub fn quantities(&mut self, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<Tree>, MaterialError> {
        let Some(v) = self.value(key) else { return Ok(None) };
        let unit = self.unit_for(key, default_unit)?;
        self.tree(key, v, dim, unit.as_deref()).map(Some)
    }

    /// One quantity in SI.
    pub fn quantity(&mut self, key: &str, dim: Dim, default_unit: Option<&str>) -> Result<Option<f64>, MaterialError> {
        match self.quantities(key, dim, default_unit)? {
            None => Ok(None),
            Some(Tree::Num(x)) => Ok(Some(x)),
            Some(_) => Err(self.err(key, "must be a single value here")),
        }
    }

    pub fn sub(&mut self, key: &str) -> Result<Option<Section<'a>>, MaterialError> {
        match self.value(key) {
            None => Ok(None),
            Some(toml::Value::Table(t)) => {
                let path = if self.path.is_empty() { key.to_string() } else { format!("{}.{key}", self.path) };
                Ok(Some(Section::new(&self.material, &path, t)))
            }
            Some(_) => Err(self.err(key, format!("must be a table, written as [{key}]"))),
        }
    }

    /// Reports the first key nobody read, with the closest known key as a suggestion.
    pub fn finish(&self, known: &[&str]) -> Result<(), MaterialError> {
        let unknown: Vec<&String> = self.table.keys().filter(|k| !self.used.contains(*k)).collect();
        if unknown.is_empty() {
            return Ok(());
        }
        let k = unknown[0];
        let best = known.iter().map(|n| (edit_distance(k, n), *n)).min();
        let hint = match best {
            Some((d, n)) if d <= 2 || (d <= 3 && k.len() > 6) => format!("; did you mean '{n}'?"),
            _ if known.is_empty() => String::new(),
            _ => format!("; known keys here: {}", known.join(", ")),
        };
        let many = if unknown.len() > 1 { format!(" ({} unknown keys in all)", unknown.len()) } else { String::new() };
        Err(self.err("", format!("unknown key '{k}'{hint}{many}")))
    }
}

fn example_unit(dim: Dim) -> &'static str {
    match dim {
        Dim::Pressure => "GPa",
        Dim::Density => "kg/m^3",
        Dim::Length => "um",
        Dim::Temperature => "K",
        Dim::ThermalConductivity => "W/(m*K)",
        Dim::SpecificHeat => "J/(kg*K)",
        Dim::PerKelvin => "1/K",
        Dim::Compliance => "1/GPa",
        Dim::ElectricalConductivity => "S/m",
        Dim::Dimensionless => "",
        Dim::Frequency => "eV",
        Dim::FrequencyHz => "GHz",
        Dim::Energy => "ueV",
        Dim::SheetResistance => "nOhm",
        Dim::SheetInductance => "pH/sq",
        Dim::PiezoelectricStress => "C/m^2",
        Dim::PiezoelectricStrain => "pC/N",
        Dim::ElectroOptic => "pm/V",
        Dim::DynamicViscosity => "mPa*s",
        Dim::ElectricField => "V/m",
    }
}
