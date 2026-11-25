// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! A material as a human-written TOML file: `MaterialSpec` supporting `format = 2` (multi-physics:
//! temperature-dependent thermal/dielectric/optical curves, BCS superconductors with sub-gap mK loss,
//! surface oxide TLS interfaces, piezoelectric Legendre duals, photoelastic, electro-optic, and
//! viscoelastic tensors) with a full read-only compatibility shim for `format = 1` files and
//! representation-aware deep merge for `extends`.

use super::bcs::{bcs_surface_impedance, BcsEval, FilmModel, Superconducting};
use super::coupled::{
    resolve_piezo_legendre_duals, Dielectric, ElectroOptic, Interfaces, Photoelastic, Piezoelectric, Viscoelastic,
};
use super::elastic::Elastic;
use super::optics::{self, Optical, OpticsError};
use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
use super::state_t::{CrystalSetting, MaterialKind, PointGroup, TempCurve, TempTensor};
use super::tensor::{Requirement, TensorProp};
use super::units::Dim;
use super::{rotate_tensor, rotation, Anisotropy, Stiffness, Tensor, M3, M6};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

/// Current primary version of the material file format (`format = 2`).
pub const FORMAT: i64 = 2;
/// Legacy v1 file format version still supported read-only via dual reader.
pub const FORMAT_V1: i64 = 1;

fn default_format() -> i64 {
    FORMAT
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Thermal {
    #[serde(default)]
    pub source: Vec<Source>,
    /// W/(m K) at reference temperature (or evaluated from `conductivity_t`)
    pub conductivity: Option<TensorProp>,
    /// J/(kg K) at reference temperature (or evaluated from `specific_heat_t`)
    pub specific_heat: Option<f64>,
    /// 1/K at reference temperature (or evaluated from `expansion_t`)
    pub expansion: Option<TensorProp>,
    /// Temperature-dependent thermal conductivity `k(T)`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conductivity_t: Option<TempTensor>,
    /// Temperature-dependent specific heat `c_p(T)`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specific_heat_t: Option<TempCurve>,
    /// Temperature-dependent instantaneous thermal expansion coefficient `alpha(T)`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expansion_t: Option<TempTensor>,
    /// Integrated thermal contraction `(L(T) - L_ref) / L_ref`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contraction: Option<TempTensor>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Electrical {
    #[serde(default)]
    pub source: Vec<Source>,
    /// relative, at low frequency
    pub permittivity: Option<TensorProp>,
    /// relative
    pub permeability: Option<f64>,
    /// S/m
    pub conductivity: Option<TensorProp>,
    /// dimensionless dielectric loss tangent `tan(delta)`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loss_tangent: Option<TensorProp>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MaterialSpec {
    #[serde(skip, default = "default_format")]
    pub format: i64,
    pub name: String,
    pub description: Option<String>,
    /// references the numbers come from
    pub source: Vec<Source>,
    pub notes: Option<String>,
    /// kg/m^3
    pub density: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<MaterialKind>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub point_group: Option<PointGroup>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub crystal_setting: Option<CrystalSetting>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_temperature: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub residual_stress: Option<TensorProp>,
    pub elastic: Option<Elastic>,
    pub thermal: Thermal,
    pub electrical: Electrical,
    pub optical: Option<Optical>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dielectric: Option<Dielectric>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub superconducting: Option<Superconducting>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub piezoelectric: Option<Piezoelectric>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub photoelastic: Option<Photoelastic>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub electro_optic: Option<ElectroOptic>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub viscoelastic: Option<Viscoelastic>,
    #[serde(default, skip_serializing_if = "Interfaces::is_empty")]
    pub interfaces: Interfaces,
}

impl Default for MaterialSpec {
    fn default() -> Self {
        Self {
            format: FORMAT,
            name: String::new(),
            description: None,
            source: Vec::new(),
            notes: None,
            density: None,
            kind: None,
            point_group: None,
            crystal_setting: None,
            reference_temperature: None,
            residual_stress: None,
            elastic: None,
            thermal: Thermal::default(),
            electrical: Electrical::default(),
            optical: None,
            dielectric: None,
            superconducting: None,
            piezoelectric: None,
            photoelastic: None,
            electro_optic: None,
            viscoelastic: None,
            interfaces: Interfaces::default(),
        }
    }
}

pub const TOP_KEYS: &[&str] = &[
    "format",
    "name",
    "description",
    "source",
    "notes",
    "extends",
    "kind",
    "point_group",
    "crystal_setting",
    "reference_temperature",
    "reference_temperature_unit",
    "density",
    "density_unit",
    "residual_stress",
    "residual_stress_unit",
    "elastic",
    "thermal",
    "electrical",
    "dielectric",
    "superconducting",
    "piezoelectric",
    "photoelastic",
    "electro_optic",
    "viscoelastic",
    "interfaces",
    "optical",
];

pub const TOP_KEYS_V1: &[&str] = &[
    "format",
    "name",
    "description",
    "source",
    "notes",
    "extends",
    "density",
    "density_unit",
    "elastic",
    "thermal",
    "electrical",
    "optical",
];

pub const ELASTIC_KEYS: &[&str] = &["source", "unit", "compliance_unit", "class", "unique_axis", "matrix", "compliance", "E", "nu", "lambda", "mu", "K", "G", "c11", "c12", "c44"];
pub const THERMAL_KEYS: &[&str] = &[
    "source",
    "conductivity",
    "conductivity_unit",
    "conductivity_T",
    "conductivity_T_unit",
    "specific_heat",
    "specific_heat_unit",
    "specific_heat_T",
    "specific_heat_T_unit",
    "expansion",
    "expansion_unit",
    "expansion_T",
    "expansion_T_unit",
    "contraction",
    "contraction_unit",
];
pub const ELECTRICAL_KEYS: &[&str] = &["source", "permittivity", "permeability", "conductivity", "conductivity_unit"];

fn positive(s: &Section, key: &str, v: Option<f64>, what: &str) -> Result<Option<f64>, MaterialError> {
    match v {
        Some(x) if !x.is_finite() || x <= 0.0 => Err(s.err(key, format!("{what} must be positive (found {x})"))),
        other => Ok(other),
    }
}

fn is_temp_curve_value(v: &toml::Value) -> bool {
    match v {
        toml::Value::Table(t) => {
            if t.contains_key("form") || t.contains_key("T") || t.contains_key("T_range") || t.contains_key("theta_e") {
                return true;
            }
            if t.contains_key("perp") && t.contains_key("par") && t.len() == 2 {
                return matches!(&t["perp"], toml::Value::Table(pt) if pt.contains_key("form") || pt.contains_key("T") || pt.contains_key("T_range") || pt.contains_key("theta_e"))
                    || matches!(&t["par"], toml::Value::Table(pt) if pt.contains_key("form") || pt.contains_key("T") || pt.contains_key("T_range") || pt.contains_key("theta_e"));
            }
            false
        }
        _ => false,
    }
}

fn m3_to_tensor(m: &M3) -> Tensor {
    let diag = m[0][1].abs() < 1e-12 && m[0][2].abs() < 1e-12 && m[1][2].abs() < 1e-12;
    if diag {
        if (m[0][0] - m[1][1]).abs() < 1e-12 * m[0][0].abs().max(1.0) {
            Tensor::Uniaxial { perp: m[0][0], par: m[2][2] }
        } else {
            Tensor::Diagonal { x: m[0][0], y: m[1][1], z: m[2][2] }
        }
    } else {
        Tensor::Symmetric {
            xx: m[0][0],
            yy: m[1][1],
            zz: m[2][2],
            xy: m[0][1],
            yz: m[1][2],
            xz: m[0][2],
        }
    }
}

impl MaterialSpec {
    /// Reads a material from TOML text. `fallback_name` is used when the file has no `name`.
    pub fn from_toml_named(text: &str, fallback_name: &str) -> Result<MaterialSpec, MaterialError> {
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        let table: toml::Table = text.parse().map_err(|e: toml::de::Error| MaterialError {
            material: fallback_name.to_string(),
            path: String::new(),
            message: format!("the file is not valid TOML: {e}"),
        })?;
        Self::from_table(&table, fallback_name)
    }

    /// Reads a material from a parsed table (already merged with its base if it had `extends`).
    pub fn from_table(table: &toml::Table, fallback_name: &str) -> Result<MaterialSpec, MaterialError> {
        let name = table.get("name").and_then(|v| v.as_str()).unwrap_or(fallback_name).to_string();
        let mut s = Section::new(&name, "", table);
        if let Some(base) = table.get("extends") {
            let base = base.as_str().unwrap_or("?");
            return Err(s.err(
                "extends",
                format!("extends = \"{base}\" can only be read together with its base material: load the file through Library::load_dir or Library::read, which find \"{base}\" in the library"),
            ));
        }

        let fmt = match s.value("format") {
            Some(toml::Value::Integer(1)) => FORMAT_V1,
            Some(toml::Value::Integer(2)) => FORMAT,
            Some(toml::Value::Integer(n)) => return Err(s.err("format", format!("this file is format {n}, but this version of nyx reads format 1 or 2"))),
            Some(_) => return Err(s.err("format", "must be a whole number; write format = 2 (or format = 1)")),
            None => return Err(s.err("format", "is missing; the first line of a material file is format = 2 (or format = 1)")),
        };
        s.string("name")?;
        let description = s.string("description")?;
        let notes = s.string("notes")?;
        let source = read_sources(&mut s, "source")?;

        let kind = match s.string("kind")? {
            Some(ref text) => Some(MaterialKind::parse(text).ok_or_else(|| {
                s.err("kind", format!("unknown material kind '{text}'; expected crystalline, polycrystalline, amorphous, superconductor, semiconductor, dielectric, metal, polymer or fluid"))
            })?),
            None => None,
        };
        let point_group = match s.string("point_group")? {
            Some(ref text) => Some(PointGroup::parse(text).ok_or_else(|| {
                s.err("point_group", format!("unknown Hermann-Mauguin point_group '{text}' (use e.g. \"m-3m\", \"-3m\", \"3m\", \"6mm\", \"isotropic\")"))
            })?),
            None => None,
        };
        let crystal_setting = match s.string("crystal_setting")? {
            Some(ref text) => Some(CrystalSetting::parse(text).ok_or_else(|| {
                s.err("crystal_setting", format!("unknown crystal_setting '{text}'; expected \"IEEE-1987\" or \"IRE-1949\""))
            })?),
            None => None,
        };
        let raw_rt = s.quantity("reference_temperature", Dim::Temperature, None)?;
        let reference_temperature = positive(&s, "reference_temperature", raw_rt, "reference_temperature")?;
        let t_ref = reference_temperature.unwrap_or(293.15);

        let density = s.quantity("density", Dim::Density, None)?;
        let density = positive(&s, "density", density, "the density")?;
        let residual_stress = TensorProp::read(&mut s, "residual_stress", Dim::Pressure, Requirement::SymmetricOnly, "residual_stress")?;

        let elastic = match s.sub("elastic")? {
            None => None,
            Some(mut e) => {
                let el = Elastic::read(&mut e)?;
                e.finish(ELASTIC_KEYS)?;
                Some(el)
            }
        };

        let mut thermal = Thermal::default();
        if let Some(mut t) = s.sub("thermal")? {
            thermal.source = read_sources(&mut t, "source")?;
            if t.value("conductivity").is_some_and(is_temp_curve_value) {
                let curve = TempTensor::read(&mut t, "conductivity", Dim::ThermalConductivity, None)?.unwrap();
                let val = curve.eval_tensor(t_ref);
                val.check(Requirement::PositiveDefinite, "thermal conductivity").map_err(|m| t.err("conductivity", m))?;
                thermal.conductivity = Some(val);
                thermal.conductivity_t = Some(curve);
            } else {
                thermal.conductivity = TensorProp::read(&mut t, "conductivity", Dim::ThermalConductivity, Requirement::PositiveDefinite, "thermal conductivity")?;
            }
            if let Some(curve) = TempTensor::read(&mut t, "conductivity_T", Dim::ThermalConductivity, None)? {
                let val = curve.eval_tensor(t_ref);
                val.check(Requirement::PositiveDefinite, "thermal.conductivity_T").map_err(|m| t.err("conductivity_T", m))?;
                if thermal.conductivity.is_none() {
                    thermal.conductivity = Some(val);
                }
                thermal.conductivity_t = Some(curve);
            }

            if t.value("specific_heat").is_some_and(is_temp_curve_value) {
                let curve = TempCurve::read(&mut t, "specific_heat", Dim::SpecificHeat, None)?.unwrap();
                let v = curve.eval(t_ref);
                if !v.is_finite() || v <= 0.0 {
                    return Err(t.err("specific_heat", format!("specific_heat({t_ref} K) must be positive (found {v})")));
                }
                thermal.specific_heat = Some(v);
                thermal.specific_heat_t = Some(curve);
            } else {
                let cp = t.quantity("specific_heat", Dim::SpecificHeat, None)?;
                thermal.specific_heat = positive(&t, "specific_heat", cp, "the specific heat")?;
            }
            if let Some(curve) = TempCurve::read(&mut t, "specific_heat_T", Dim::SpecificHeat, None)? {
                let v = curve.eval(t_ref);
                if !v.is_finite() || v <= 0.0 {
                    return Err(t.err("specific_heat_T", format!("specific_heat_T({t_ref} K) must be positive (found {v})")));
                }
                if thermal.specific_heat.is_none() {
                    thermal.specific_heat = Some(v);
                }
                thermal.specific_heat_t = Some(curve);
            }

            if t.value("expansion").is_some_and(is_temp_curve_value) {
                let curve = TempTensor::read(&mut t, "expansion", Dim::PerKelvin, None)?.unwrap();
                let val = curve.eval_tensor(t_ref);
                val.check(Requirement::SymmetricOnly, "thermal expansion").map_err(|m| t.err("expansion", m))?;
                thermal.expansion = Some(val);
                thermal.expansion_t = Some(curve);
            } else {
                thermal.expansion = TensorProp::read(&mut t, "expansion", Dim::PerKelvin, Requirement::SymmetricOnly, "thermal expansion")?;
            }
            if let Some(curve) = TempTensor::read(&mut t, "expansion_T", Dim::PerKelvin, None)? {
                let val = curve.eval_tensor(t_ref);
                val.check(Requirement::SymmetricOnly, "thermal.expansion_T").map_err(|m| t.err("expansion_T", m))?;
                if thermal.expansion.is_none() {
                    thermal.expansion = Some(val);
                }
                thermal.expansion_t = Some(curve);
            }

            thermal.contraction = TempTensor::read(&mut t, "contraction", Dim::Dimensionless, None)?;
            t.finish(THERMAL_KEYS)?;
        }

        let mut electrical = Electrical::default();
        if let Some(mut t) = s.sub("electrical")? {
            electrical.source = read_sources(&mut t, "source")?;
            electrical.permittivity = TensorProp::read(&mut t, "permittivity", Dim::Dimensionless, Requirement::PositiveDefinite, "relative permittivity")?;
            let mu = t.number("permeability")?;
            electrical.permeability = positive(&t, "permeability", mu, "the relative permeability")?;
            electrical.conductivity = TensorProp::read(&mut t, "conductivity", Dim::ElectricalConductivity, Requirement::NonNegative, "electrical conductivity")?;
            if electrical.conductivity.is_none() {
                electrical.conductivity = TensorProp::read(&mut t, "sigma", Dim::ElectricalConductivity, Requirement::NonNegative, "electrical conductivity")?;
            }
            electrical.loss_tangent = TensorProp::read(&mut t, "loss_tangent", Dim::Dimensionless, Requirement::NonNegative, "dielectric loss tangent")?;
            if electrical.loss_tangent.is_none() {
                electrical.loss_tangent = TensorProp::read(&mut t, "tan_delta", Dim::Dimensionless, Requirement::NonNegative, "dielectric loss tangent")?;
            }
            const ELECTRICAL_PARSE_KEYS: &[&str] = &[
                "source",
                "permittivity",
                "permeability",
                "conductivity",
                "conductivity_unit",
                "sigma",
                "sigma_unit",
                "loss_tangent",
                "tan_delta",
            ];
            t.finish(ELECTRICAL_PARSE_KEYS)?;
        }

        let dielectric = s.sub("dielectric")?.map(|mut d| Dielectric::read(&mut d)).transpose()?;
        let superconducting = s.sub("superconducting")?.map(|mut sc| Superconducting::read(&mut sc)).transpose()?;
        let piezoelectric = s.sub("piezoelectric")?.map(|mut pz| Piezoelectric::read(&mut pz, point_group)).transpose()?;
        let photoelastic = s.sub("photoelastic")?.map(|mut pe| Photoelastic::read(&mut pe, point_group)).transpose()?;
        let electro_optic = s.sub("electro_optic")?.map(|mut eo| ElectroOptic::read(&mut eo, point_group)).transpose()?;
        let viscoelastic = s.sub("viscoelastic")?.map(|mut ve| Viscoelastic::read(&mut ve, point_group)).transpose()?;
        let interfaces = s.sub("interfaces")?.map(|mut iface| Interfaces::read(&mut iface)).transpose()?.unwrap_or_default();

        if let Some(pg) = point_group.filter(|pg| pg.forbids_piezo_and_pockels()) {
            if piezoelectric.is_some() {
                return Err(s.err(
                    "piezoelectric",
                    format!("point_group \"{}\" is centrosymmetric (or chiral cubic 432) and forbids linear piezoelectricity by Neumann's principle", pg.as_str()),
                ));
            }
            if electro_optic.is_some() {
                return Err(s.err(
                    "electro_optic",
                    format!("point_group \"{}\" forbids the linear electro-optic (Pockels) effect by Neumann's principle", pg.as_str()),
                ));
            }
        }

        let optical = match s.sub("optical")? {
            None => None,
            Some(mut o) => {
                let opt = Optical::read(&mut o)?;
                let mut known: Vec<&str> = optics::OPTICAL_KEYS.to_vec();
                known.extend(optics::MODEL_KEYS);
                o.finish(&known)?;
                Some(opt)
            }
        };

        if fmt == FORMAT_V1 {
            // Retain v1 surface compatibility while still permitting [superconducting] from sibling branch
            let v1_plus: Vec<&str> = TOP_KEYS_V1.iter().copied().chain(std::iter::once("superconducting")).collect();
            s.finish(&v1_plus)?;
        } else {
            s.finish(TOP_KEYS)?;
        }

        Ok(MaterialSpec {
            format: fmt,
            name,
            description,
            source,
            notes,
            density,
            kind,
            point_group,
            crystal_setting,
            reference_temperature,
            residual_stress,
            elastic,
            thermal,
            electrical,
            optical,
            dielectric,
            superconducting,
            piezoelectric,
            photoelastic,
            electro_optic,
            viscoelastic,
            interfaces,
        })
    }

    pub fn from_toml(text: &str) -> Result<MaterialSpec, MaterialError> {
        Self::from_toml_named(text, "unnamed")
    }

    /// Reads a file; the file name (without `.toml`) is the name if the file has none.
    pub fn load(path: &Path) -> Result<MaterialSpec, MaterialError> {
        let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("unnamed").to_string();
        let text = std::fs::read_to_string(path).map_err(|e| MaterialError { material: stem.clone(), path: String::new(), message: format!("cannot read {}: {e}", path.display()) })?;
        Self::from_toml_named(&text, &stem).map_err(|mut e| {
            if e.material != stem {
                e.material = format!("{} (file {stem}.toml)", e.material);
            }
            e
        })
    }

    /// The material as TOML text in `format = 2` (comments of an original file are not kept).
    pub fn to_toml(&self) -> String {
        let mut t = toml::Table::new();
        t.insert("format".into(), toml::Value::Integer(FORMAT));
        t.insert("name".into(), self.name.clone().into());
        if let Some(d) = &self.description {
            t.insert("description".into(), d.clone().into());
        }
        if let Some(k) = self.kind {
            t.insert("kind".into(), k.as_str().into());
        }
        if let Some(pg) = self.point_group {
            t.insert("point_group".into(), pg.as_str().into());
        }
        if let Some(cs) = self.crystal_setting {
            t.insert("crystal_setting".into(), cs.as_str().into());
        }
        if let Some(rt) = self.reference_temperature {
            t.insert("reference_temperature".into(), format!("{rt} K").into());
        }
        if let Some(v) = write_sources(&self.source) {
            t.insert("source".into(), v);
        }
        if let Some(n) = &self.notes {
            t.insert("notes".into(), n.clone().into());
        }
        if let Some(d) = self.density {
            t.insert("density".into(), format!("{d} kg/m^3").into());
        }
        put_tensor(&mut t, "residual_stress", self.residual_stress.as_ref(), Some("Pa"));
        if let Some(e) = &self.elastic {
            t.insert("elastic".into(), toml::Value::Table(e.write()));
        }
        let mut th = toml::Table::new();
        if let Some(v) = write_sources(&self.thermal.source) {
            th.insert("source".into(), v);
        }
        put_tensor(&mut th, "conductivity", self.thermal.conductivity.as_ref(), Some("W/(m*K)"));
        if let Some(ct) = &self.thermal.conductivity_t {
            th.insert("conductivity_T".into(), ct.to_toml(Some("W/(m*K)")));
        }
        if let Some(c) = self.thermal.specific_heat {
            th.insert("specific_heat".into(), format!("{c} J/(kg*K)").into());
        }
        if let Some(cpt) = &self.thermal.specific_heat_t {
            th.insert("specific_heat_T".into(), cpt.to_toml(Some("J/(kg*K)")));
        }
        put_tensor(&mut th, "expansion", self.thermal.expansion.as_ref(), Some("1/K"));
        if let Some(at) = &self.thermal.expansion_t {
            th.insert("expansion_T".into(), at.to_toml(Some("1/K")));
        }
        if let Some(con) = &self.thermal.contraction {
            th.insert("contraction".into(), con.to_toml(None));
        }
        if !th.is_empty() {
            t.insert("thermal".into(), toml::Value::Table(th));
        }
        let mut el = toml::Table::new();
        if let Some(v) = write_sources(&self.electrical.source) {
            el.insert("source".into(), v);
        }
        put_tensor(&mut el, "permittivity", self.electrical.permittivity.as_ref(), None);
        if let Some(m) = self.electrical.permeability {
            el.insert("permeability".into(), toml::Value::Float(m));
        }
        put_tensor(&mut el, "conductivity", self.electrical.conductivity.as_ref(), Some("S/m"));
        put_tensor(&mut el, "loss_tangent", self.electrical.loss_tangent.as_ref(), None);
        if !el.is_empty() {
            t.insert("electrical".into(), toml::Value::Table(el));
        }
        if let Some(dl) = &self.dielectric {
            let mut tbl = toml::Table::new();
            dl.write(&mut tbl);
            if !tbl.is_empty() {
                t.insert("dielectric".into(), toml::Value::Table(tbl));
            }
        }
        if let Some(sc) = &self.superconducting {
            let mut tbl = toml::Table::new();
            sc.write(&mut tbl);
            t.insert("superconducting".into(), toml::Value::Table(tbl));
        }
        if let Some(pz) = &self.piezoelectric {
            let mut tbl = toml::Table::new();
            pz.write(&mut tbl);
            t.insert("piezoelectric".into(), toml::Value::Table(tbl));
        }
        if let Some(pe) = &self.photoelastic {
            let mut tbl = toml::Table::new();
            pe.write(&mut tbl);
            t.insert("photoelastic".into(), toml::Value::Table(tbl));
        }
        if let Some(eo) = &self.electro_optic {
            let mut tbl = toml::Table::new();
            eo.write(&mut tbl);
            t.insert("electro_optic".into(), toml::Value::Table(tbl));
        }
        if let Some(ve) = &self.viscoelastic {
            let mut tbl = toml::Table::new();
            ve.write(&mut tbl);
            t.insert("viscoelastic".into(), toml::Value::Table(tbl));
        }
        if !self.interfaces.is_empty() {
            let mut tbl = toml::Table::new();
            self.interfaces.write(&mut tbl);
            t.insert("interfaces".into(), toml::Value::Table(tbl));
        }
        if let Some(o) = &self.optical {
            t.insert("optical".into(), toml::Value::Table(o.write()));
        }
        toml::to_string(&t).expect("a table of plain values always serialises")
    }

    // ----- use by the solvers -----

    /// Clamped (constant-strain) static/microwave relative permittivity `eps^S` in the crystal frame.
    pub fn permittivity_clamped(&self) -> Option<M3> {
        let dl_c = self
            .dielectric
            .as_ref()
            .and_then(|d| d.eps_clamped.as_ref().or(d.eps_r.as_ref()))
            .map(TensorProp::matrix);
        if dl_c.is_some() {
            return dl_c;
        }
        let dl_f = self.dielectric.as_ref().and_then(|d| d.eps_free.as_ref()).map(TensorProp::matrix);
        if let (Some(ef), Some(el), Some(pz)) = (dl_f, &self.elastic, &self.piezoelectric) {
            if let Ok((ec, _, _, _)) = resolve_piezo_legendre_duals(&el.c, pz, None, Some(&ef)) {
                return Some(ec);
            }
        }
        dl_f.or_else(|| self.electrical.permittivity.as_ref().map(TensorProp::matrix))
    }

    /// Free (constant-stress) static/microwave relative permittivity `eps^T` in the crystal frame.
    pub fn permittivity_free(&self) -> Option<M3> {
        let dl_f = self.dielectric.as_ref().and_then(|d| d.eps_free.as_ref()).map(TensorProp::matrix);
        if dl_f.is_some() {
            return dl_f;
        }
        let ec = self
            .dielectric
            .as_ref()
            .and_then(|d| d.eps_clamped.as_ref().or(d.eps_r.as_ref()))
            .or(self.electrical.permittivity.as_ref())
            .map(TensorProp::matrix);
        if let (Some(ec_m), Some(el), Some(pz)) = (ec, &self.elastic, &self.piezoelectric) {
            if let Ok((_, ef, _, _)) = resolve_piezo_legendre_duals(&el.c, pz, Some(&ec_m), None) {
                return Some(ef);
            }
        }
        ec
    }

    /// Short-circuit (constant-electric-field) stiffness `C^E` in textbook Voigt order `(xx, yy, zz, yz, xz, xy)` (`Pa`).
    pub fn stiffness_constant_e(&self) -> Option<M6> {
        self.elastic.as_ref().map(|e| e.c)
    }

    /// Open-circuit (constant-displacement) stiffness `C^D = C^E + e^T (eps^S)^-1 e` in textbook Voigt order `(xx, yy, zz, yz, xz, xy)` (`Pa`).
    pub fn stiffness_constant_d(&self) -> Option<M6> {
        let el = self.elastic.as_ref()?;
        if let Some(pz) = &self.piezoelectric {
            let ec = self
                .dielectric
                .as_ref()
                .and_then(|d| d.eps_clamped.as_ref().or(d.eps_r.as_ref()))
                .or(self.electrical.permittivity.as_ref())
                .map(TensorProp::matrix);
            let ef = self.dielectric.as_ref().and_then(|d| d.eps_free.as_ref()).map(TensorProp::matrix);
            if let Ok((_, _, _, cd_voigt)) = resolve_piezo_legendre_duals(&el.c, pz, ec.as_ref(), ef.as_ref()) {
                return Some(cd_voigt);
            }
        }
        Some(el.c)
    }

    /// Evaluates global-frame microwave complex permittivity `(eps_re, eps_im_loss)` where
    /// `eps_tilde = eps_re - i eps_im_loss` (`exp(+i omega t)`) and `eps_im_loss >= 0`.
    pub fn complex_microwave_permittivity(
        &self,
        freq_hz: Option<f64>,
        temp_k: Option<f64>,
        e_field_v_per_m: Option<f64>,
        prefer_free: bool,
        euler_deg: [f64; 3],
    ) -> Option<(M3, M3)> {
        let (mut re_cryst, mut im_cryst) = if let Some(dl) = &self.dielectric {
            if let Some((re, im)) = dl.complex_permittivity(freq_hz, temp_k, e_field_v_per_m, prefer_free) {
                (re, im)
            } else {
                let re = if prefer_free { self.permittivity_free()? } else { self.permittivity_clamped()? };
                (re, [[0.0; 3]; 3])
            }
        } else {
            let re = if prefer_free { self.permittivity_free()? } else { self.permittivity_clamped()? };
            (re, [[0.0; 3]; 3])
        };
        if prefer_free && self.dielectric.as_ref().and_then(|d| d.eps_free.as_ref()).is_none() {
            if let Some(free_re) = self.permittivity_free() {
                re_cryst = free_re;
            }
        }
        if !im_cryst.iter().flatten().any(|&v| v > 0.0) {
            if let Some(td) = &self.electrical.loss_tangent {
                let td_m = td.matrix();
                for i in 0..3 {
                    for j in 0..3 {
                        let mut s1 = 0.0;
                        let mut s2 = 0.0;
                        for k in 0..3 {
                            s1 += re_cryst[i][k] * td_m[k][j];
                            s2 += td_m[i][k] * re_cryst[k][j];
                        }
                        im_cryst[i][j] = 0.5 * (s1 + s2);
                    }
                }
            }
        }
        let r = rotation(euler_deg);
        Some((rotate_tensor(&re_cryst, &r), rotate_tensor(&im_cryst, &r)))
    }

    /// Evaluates thermal conductivity `k(T)` (W/(m*K)) at `temp_k`, defaulting to reference temperature.
    pub fn thermal_conductivity_at(&self, temp_k: Option<f64>) -> Option<TensorProp> {
        match (temp_k, &self.thermal.conductivity_t) {
            (Some(t), Some(curve)) => Some(curve.eval_tensor(t)),
            _ => self.thermal.conductivity.clone(),
        }
    }

    /// Evaluates specific heat `c_p(T)` (J/(kg*K)) at `temp_k`, defaulting to reference temperature.
    pub fn specific_heat_at(&self, temp_k: Option<f64>) -> Option<f64> {
        match (temp_k, &self.thermal.specific_heat_t) {
            (Some(t), Some(curve)) => Some(curve.eval(t)),
            _ => self.thermal.specific_heat,
        }
    }

    /// Evaluates thermal expansion `alpha(T)` (1/K) at `temp_k`, defaulting to reference temperature.
    pub fn thermal_expansion_at(&self, temp_k: Option<f64>) -> Option<TensorProp> {
        match (temp_k, &self.thermal.expansion_t) {
            (Some(t), Some(curve)) => Some(curve.eval_tensor(t)),
            _ => self.thermal.expansion.clone(),
        }
    }

    /// Convenience BCS surface impedance evaluator using `[superconducting]` defaults when omitted.
    pub fn surface_impedance_at(
        &self,
        freq_hz: f64,
        temp_k: f64,
        thickness: Option<f64>,
        film_model: Option<FilmModel>,
    ) -> Result<BcsEval, MaterialError> {
        let sc = self.superconducting.as_ref().ok_or_else(|| MaterialError {
            material: self.name.clone(),
            path: "superconducting".into(),
            message: "this material has no [superconducting] section".into(),
        })?;
        let d = thickness.or(sc.nominal_thickness).unwrap_or(100e-9);
        let fm = film_model.or(sc.film_model).unwrap_or_default();
        bcs_surface_impedance(sc, freq_hz, temp_k, d, fm).map_err(|message| MaterialError {
            material: self.name.clone(),
            path: "superconducting".into(),
            message,
        })
    }

    /// The material in the `Anisotropy` form placed with Bunge angles in degrees.
    pub fn anisotropy(&self, euler_deg: [f64; 3]) -> Anisotropy {
        let perm_cryst = self.permittivity_clamped();
        let perm_im = self
            .dielectric
            .as_ref()
            .and_then(|d| d.complex_permittivity(None, self.reference_temperature, None, false))
            .map(|(_, im)| im)
            .filter(|im| im.iter().flatten().any(|&v| v > 0.0))
            .or_else(|| {
                let re_cryst = perm_cryst?;
                let td_m = self.electrical.loss_tangent.as_ref()?.matrix();
                let mut im = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        let mut s1 = 0.0;
                        let mut s2 = 0.0;
                        for k in 0..3 {
                            s1 += re_cryst[i][k] * td_m[k][j];
                            s2 += td_m[i][k] * re_cryst[k][j];
                        }
                        im[i][j] = 0.5 * (s1 + s2);
                    }
                }
                im.iter().flatten().any(|&v| v > 0.0).then_some(im)
            });
        let ce_voigt = self.elastic.as_ref().map(|e| e.c);
        let piezo_e = self
            .piezoelectric
            .as_ref()
            .and_then(|pz| pz.stress_charge(ce_voigt.as_ref()).ok());
        Anisotropy {
            stiffness: self.elastic.as_ref().map(|e| Stiffness::General(e.c)),
            conductivity: self
                .thermal
                .conductivity
                .as_ref()
                .or(self.electrical.conductivity.as_ref())
                .map(TensorProp::to_legacy),
            permittivity: perm_cryst.map(|m| m3_to_tensor(&m)),
            expansion: self.thermal.expansion.as_ref().map(TensorProp::to_legacy),
            euler_deg,
            permittivity_im: perm_im.map(|m| m3_to_tensor(&m)),
            viscosity: self.viscoelastic.as_ref().map(|v| v.matrix),
            piezo_e,
            photoelastic: self.photoelastic.as_ref().map(|pe| pe.matrix),
        }
    }

    /// Scalar properties under the older names (`rho`, `Cp`, `mu_r`, and `eps_r`, `k`, `alpha`, `E`, `nu` where the data are isotropic).
    pub fn legacy_scalars(&self) -> BTreeMap<String, f64> {
        let mut m = BTreeMap::new();
        let mut put = |k: &str, v: Option<f64>| {
            if let Some(v) = v {
                m.insert(k.to_string(), v);
            }
        };
        put("rho", self.density);
        put("Cp", self.thermal.specific_heat);
        put("mu_r", self.electrical.permeability);
        put(
            "eps_r",
            self.electrical.permittivity.as_ref().and_then(TensorProp::as_scalar).or_else(|| {
                self.dielectric
                    .as_ref()
                    .and_then(|d| d.eps_clamped.as_ref().or(d.eps_free.as_ref()).or(d.eps_r.as_ref()))
                    .and_then(TensorProp::as_scalar)
            }),
        );
        put("sigma", self.electrical.conductivity.as_ref().and_then(TensorProp::as_scalar));
        put("k", self.thermal.conductivity.as_ref().and_then(TensorProp::as_scalar));
        put("alpha", self.thermal.expansion.as_ref().and_then(TensorProp::as_scalar));
        if let Some(e) = &self.elastic {
            if e.detect_class() == super::Class::Isotropic {
                let (lam, mu) = (e.c[0][1], 0.5 * (e.c[0][0] - e.c[0][1]));
                put("E", Some(mu * (3.0 * lam + 2.0 * mu) / (lam + mu)));
                put("nu", Some(lam / (2.0 * (lam + mu))));
            }
        }
        put(
            "tan_delta",
            self.dielectric
                .as_ref()
                .and_then(|dl| dl.loss_tangent.as_ref())
                .or(self.electrical.loss_tangent.as_ref())
                .and_then(TensorProp::as_scalar),
        );
        if let Some(dl) = &self.dielectric {
            put("tls_tan_delta_0", dl.tls_tan_delta_0);
        }
        if let Some(sc) = &self.superconducting {
            put("Tc", Some(sc.tc));
        }
        m
    }

    /// Relative permittivity tensor at a wavelength (m) and temperature (K) in the frame of a crystal rotated by Bunge angles.
    pub fn optical_permittivity(&self, lambda_m: f64, temp_k: Option<f64>, euler_deg: [f64; 3]) -> Result<[[Complex64; 3]; 3], MaterialError> {
        let o = self.optical.as_ref().ok_or_else(|| MaterialError {
            material: self.name.clone(),
            path: "optical".into(),
            message: "this material has no [optical] section, so it has no refractive index; add one (for example model = \"constant\", n = 1.45)".into(),
        })?;
        o.permittivity_global(lambda_m, temp_k, euler_deg).map_err(|e: OpticsError| MaterialError { material: self.name.clone(), path: "optical".into(), message: e.to_string() })
    }
}

fn put_tensor(t: &mut toml::Table, key: &str, v: Option<&TensorProp>, unit: Option<&str>) {
    if let Some(v) = v {
        t.insert(key.into(), v.write(1.0));
        if let Some(u) = unit {
            t.insert(format!("{key}_unit"), u.into());
        }
    }
}

/// A folder of material files, loaded at run time. Materials are found by file name (without `.toml`), ignoring case.
/// A file may start from another material of the library with `extends = "name"` (see `merge_over`).
#[derive(Clone, Debug, Default)]
pub struct Library {
    materials: BTreeMap<String, MaterialSpec>,
    /// the tables after `extends` was resolved, which later files and `read` extend in turn
    tables: BTreeMap<String, toml::Table>,
}

/// The library names closest to `name`: all that tie for the smallest edit distance, which must be at most 3.
fn nearest<'a>(name: &str, names: impl Iterator<Item = &'a String>) -> Vec<&'a String> {
    let scored: Vec<(usize, &String)> = names.map(|n| (super::reader::edit_distance(name, n), n)).collect();
    let best = scored.iter().map(|(d, _)| *d).min().unwrap_or(usize::MAX);
    if best > 3 {
        return Vec::new();
    }
    scored.into_iter().filter(|(d, _)| *d == best).map(|(_, n)| n).collect()
}

impl Library {
    /// Loads every `*.toml` file of the folder (not its subfolders). All files with problems are reported together.
    pub fn load_dir(dir: &Path) -> Result<Library, Vec<MaterialError>> {
        let entries = std::fs::read_dir(dir).map_err(|e| {
            vec![MaterialError { material: dir.display().to_string(), path: String::new(), message: format!("cannot read the material folder: {e}") }]
        })?;
        let mut files: Vec<_> = entries.filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.extension().and_then(|x| x.to_str()) == Some("toml")).collect();
        files.sort();
        let (mut raw, mut errors) = (BTreeMap::new(), Vec::new());
        for f in files {
            let stem = f.file_stem().and_then(|s| s.to_str()).unwrap_or_default().to_string();
            let text = match std::fs::read_to_string(&f) {
                Ok(t) => t,
                Err(e) => {
                    errors.push(MaterialError { material: stem, path: String::new(), message: format!("cannot read {}: {e}", f.display()) });
                    continue;
                }
            };
            let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
            match text.parse::<toml::Table>() {
                Ok(t) => {
                    raw.insert(stem.to_ascii_lowercase(), (stem, t));
                }
                Err(e) => errors.push(MaterialError { material: stem, path: String::new(), message: format!("the file is not valid TOML: {e}") }),
            }
        }
        let (mut materials, mut tables) = (BTreeMap::new(), BTreeMap::new());
        for key in raw.keys().cloned().collect::<Vec<_>>() {
            let stem = raw[&key].0.clone();
            let merged = resolve(&key, &raw, &mut Vec::new(), &mut tables);
            let result = merged.and_then(|t| {
                let extends = raw[&key].1.get("extends").and_then(|v| v.as_str()).map(str::to_string);
                MaterialSpec::from_table(&t, &stem).map_err(|mut e| {
                    if let Some(b) = extends {
                        e.message.push_str(&format!(" (this material extends \"{b}\"; the combined result is checked)"));
                    }
                    if e.material != stem {
                        e.material = format!("{} (file {stem}.toml)", e.material);
                    }
                    e
                })
            });
            match result {
                Ok(m) => {
                    materials.insert(key, m);
                }
                Err(e) => errors.push(e),
            }
        }
        if errors.is_empty() {
            Ok(Library { materials, tables })
        } else {
            Err(errors)
        }
    }

    /// Reads material text that may use `extends`, looking the base up in this library (for a file outside the library folder,
    /// for example one embedded in a project).
    pub fn read(&self, text: &str, fallback_name: &str) -> Result<MaterialSpec, MaterialError> {
        self.read_table(text, fallback_name).map(|(m, _)| m)
    }

    /// A copy of the library with more materials given as TOML text by name (the materials a project carries). A text may extend a material
    /// of the library or one given before it (in name order). All texts with problems are reported together.
    pub fn with_texts(&self, texts: &BTreeMap<String, String>) -> Result<Library, Vec<MaterialError>> {
        let mut lib = self.clone();
        let mut errors = Vec::new();
        for (name, text) in texts {
            match lib.read_table(text, name) {
                Ok((m, t)) => {
                    let key = name.to_ascii_lowercase();
                    lib.materials.insert(key.clone(), m);
                    lib.tables.insert(key, t);
                }
                Err(e) => errors.push(e),
            }
        }
        if errors.is_empty() {
            Ok(lib)
        } else {
            Err(errors)
        }
    }

    /// `read`, with the table after `extends` was resolved.
    pub fn read_table(&self, text: &str, fallback_name: &str) -> Result<(MaterialSpec, toml::Table), MaterialError> {
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        let table: toml::Table = text.parse().map_err(|e: toml::de::Error| MaterialError {
            material: fallback_name.to_string(),
            path: String::new(),
            message: format!("the file is not valid TOML: {e}"),
        })?;
        let merged = match table.get("extends") {
            None => table,
            Some(b) => {
                let base = b.as_str().ok_or_else(|| MaterialError { material: fallback_name.to_string(), path: "extends".into(), message: "must be the name of a library material in quotes".into() })?;
                let bt = self.tables.get(&base.to_ascii_lowercase()).ok_or_else(|| missing_base(fallback_name, base, self.tables.keys()))?;
                merge_over(bt, &table)
            }
        };
        let spec = MaterialSpec::from_table(&merged, fallback_name)?;
        Ok((spec, merged))
    }

    pub fn get(&self, name: &str) -> Option<&MaterialSpec> {
        self.materials.get(&name.to_ascii_lowercase())
    }

    /// Like `get`, with an error that lists the available names.
    pub fn require(&self, name: &str) -> Result<&MaterialSpec, MaterialError> {
        self.get(name).ok_or_else(|| MaterialError {
            material: name.to_string(),
            path: String::new(),
            message: format!("not in the library; available: {}", self.names().join(", ")),
        })
    }

    pub fn names(&self) -> Vec<String> {
        self.materials.keys().cloned().collect()
    }

    pub fn len(&self) -> usize {
        self.materials.len()
    }

    pub fn is_empty(&self) -> bool {
        self.materials.is_empty()
    }
}

fn missing_base<'a>(who: &str, base: &str, names: impl Iterator<Item = &'a String> + Clone) -> MaterialError {
    let close = nearest(&base.to_ascii_lowercase(), names.clone());
    let hint = if close.is_empty() { String::new() } else { format!("; did you mean {}?", close.iter().map(|n| format!("\"{n}\"")).collect::<Vec<_>>().join(" or ")) };
    MaterialError {
        material: who.to_string(),
        path: "extends".into(),
        message: format!("extends = \"{base}\", but the library has no material of that name{hint} (available: {})", names.cloned().collect::<Vec<_>>().join(", ")),
    }
}

/// The table of `key` with its `extends` chain resolved. `stack` is the chain being resolved, to detect loops.
fn resolve(key: &str, raw: &BTreeMap<String, (String, toml::Table)>, stack: &mut Vec<String>, done: &mut BTreeMap<String, toml::Table>) -> Result<toml::Table, MaterialError> {
    if let Some(t) = done.get(key) {
        return Ok(t.clone());
    }
    let (stem, table) = &raw[key];
    let Some(base) = table.get("extends") else {
        done.insert(key.to_string(), table.clone());
        return Ok(table.clone());
    };
    let err = |m: String| MaterialError { material: stem.clone(), path: "extends".into(), message: m };
    let base = base.as_str().ok_or_else(|| err("must be the name of a library material in quotes".into()))?;
    let bkey = base.to_ascii_lowercase();
    stack.push(key.to_string());
    if stack.contains(&bkey) {
        let mut chain: Vec<String> = stack.iter().map(|k| raw[k].0.clone()).collect();
        chain.push(raw[&bkey].0.clone());
        stack.pop();
        return Err(err(format!("the extends chain forms a loop: {}. A material cannot start from itself, directly or through others", chain.join(" -> "))));
    }
    if !raw.contains_key(&bkey) {
        stack.pop();
        return Err(missing_base(stem, base, raw.keys()));
    }
    let base_table = resolve(&bkey, raw, stack, done);
    stack.pop();
    let base_table = base_table.map_err(|e| err(format!("its base \"{base}\" has a problem: {}", e.message)))?;
    let merged = merge_over(&base_table, table);
    done.insert(key.to_string(), merged.clone());
    Ok(merged)
}

fn is_elastic_rep_key(k: &str) -> bool {
    if matches!(k, "matrix" | "compliance" | "class" | "unique_axis" | "E" | "nu" | "lambda" | "mu" | "K" | "G") {
        return true;
    }
    let b = k.as_bytes();
    b.len() == 3 && (b[0] == b'c' || b[0] == b's') && (b'1'..=b'6').contains(&b[1]) && (b'1'..=b'6').contains(&b[2])
}

fn strip_unit_if_unqualified(out: &mut toml::Table, over: &toml::Table, k: &str) {
    let ukey = format!("{k}_unit");
    if !over.contains_key(&ukey) {
        out.remove(&ukey);
    }
}

fn merge_subtable(section: &str, base: &toml::Table, over: &toml::Table) -> toml::Table {
    let mut r = base.clone();

    if section == "elastic" && over.keys().any(|k| is_elastic_rep_key(k)) {
        r.retain(|k, _| !is_elastic_rep_key(k) && k != "unit" && k != "compliance_unit");
    }

    if section == "optical" || section.starts_with("optical.") {
        // Reset conflicting optical model keys when changing model or axis structure
        let model_changed = over.get("model").is_some_and(|m| r.get("model") != Some(m));
        if model_changed {
            for k in ["n", "k", "eps", "A", "B", "C", "lambda0", "coefficients", "data", "interpolation", "extrapolate", "frequency_unit", "eps_inf", "plasma_frequency", "damping", "drude_weight", "oscillators"] {
                r.remove(k);
            }
        }
        let iso_over = over.contains_key("model");
        let axes_over = over.contains_key("ordinary") || over.contains_key("extraordinary") || over.contains_key("x") || over.contains_key("y") || over.contains_key("z");
        if iso_over && section == "optical" {
            for k in ["ordinary", "extraordinary", "x", "y", "z"] {
                r.remove(k);
            }
        } else if axes_over && section == "optical" {
            for k in optics::MODEL_KEYS {
                r.remove(*k);
            }
        }
    }

    if section == "piezoelectric" {
        let form_changed = over.get("form").is_some_and(|f| r.get("form") != Some(f));
        let rep_changed = over.contains_key("matrix") || over.keys().any(|k| k.starts_with('e') || k.starts_with('d'));
        if form_changed || rep_changed {
            r.remove("matrix");
            r.retain(|k, _| !((k.starts_with('e') || k.starts_with('d')) && k.len() == 3 && k.as_bytes()[1].is_ascii_digit()));
        }
    }

    if matches!(section, "photoelastic" | "viscoelastic" | "electro_optic" | "electro_optic.clamped" | "electro_optic.free") {
        let rep_changed = over.contains_key("matrix") || over.keys().any(|k| k.len() == 3 && (k.starts_with('p') || k.starts_with('r')) || (k.len() == 5 && k.starts_with("eta")));
        if rep_changed {
            r.remove("matrix");
            r.retain(|k, _| !((k.starts_with('p') || k.starts_with('r')) && k.len() == 3 && k.as_bytes()[1].is_ascii_digit()) && !k.starts_with("eta"));
        }
    }

    for (k, v) in over {
        let sub_path = format!("{section}.{k}");
        if let (Some(toml::Value::Table(bt)), toml::Value::Table(ot)) = (r.get(k), v) {
            // Only recurse into structured sub-tables or compatible temp-curve form; if `form` changed, replace whole curve
            if ot.get("form").is_some() && bt.get("form") != ot.get("form") {
                r.insert(k.clone(), v.clone());
            } else {
                r.insert(k.clone(), toml::Value::Table(merge_subtable(&sub_path, bt, ot)));
            }
        } else {
            r.insert(k.clone(), v.clone());
        }
        strip_unit_if_unqualified(&mut r, over, k);
    }
    r
}

/// `over` laid over `base`:
/// - When `over` (or `base` if `over` omits `format`) is `format = 1`: legacy whole-table replacement (ADR 30).
/// - When `over` is `format = 2`: representation-aware deep merge across nested tables.
///
/// The merged table carries `max(base.format, over.format)` so a legacy `format = 1` child extending a
/// `format = 2` library material inherits all library sections without schema-key rejection.
pub fn merge_over(base: &toml::Table, over: &toml::Table) -> toml::Table {
    let base_fmt = base.get("format").and_then(|v| v.as_integer()).unwrap_or(FORMAT_V1);
    let over_fmt = over.get("format").and_then(|v| v.as_integer()).unwrap_or(base_fmt);

    if over_fmt == FORMAT_V1 {
        let mut r = base.clone();
        r.remove("extends");
        r.remove("name");
        for (k, v) in over {
            if k == "extends" {
                continue;
            }
            r.insert(k.clone(), v.clone());
            if k == "density" && !over.contains_key("density_unit") {
                r.remove("density_unit");
            }
        }
        r.insert("format".into(), toml::Value::Integer(base_fmt.max(over_fmt)));
        return r;
    }

    let mut r = base.clone();
    r.remove("extends");
    r.remove("name");
    for (k, v) in over {
        if k == "extends" {
            continue;
        }
        if let (Some(toml::Value::Table(bt)), toml::Value::Table(ot)) = (r.get(k), v) {
            r.insert(k.clone(), toml::Value::Table(merge_subtable(k, bt, ot)));
        } else {
            r.insert(k.clone(), v.clone());
        }
        strip_unit_if_unqualified(&mut r, over, k);
    }
    r.insert("format".into(), toml::Value::Integer(base_fmt.max(over_fmt)));
    r
}

/// The keys each table of a material file accepts (the contract test compares them with the JSON Schema).
pub struct Keys {
    pub top: &'static [&'static str],
    pub top_v1: &'static [&'static str],
    pub elastic: &'static [&'static str],
    pub thermal: &'static [&'static str],
    pub electrical: &'static [&'static str],
    pub dielectric: &'static [&'static str],
    pub superconducting: &'static [&'static str],
    pub piezoelectric: &'static [&'static str],
    pub photoelastic: &'static [&'static str],
    pub electro_optic: &'static [&'static str],
    pub viscoelastic: &'static [&'static str],
    pub interfaces: &'static [&'static str],
    pub interface_layer: &'static [&'static str],
    pub optical: &'static [&'static str],
    pub model: &'static [&'static str],
}

pub fn keys() -> Keys {
    Keys {
        top: TOP_KEYS,
        top_v1: TOP_KEYS_V1,
        elastic: ELASTIC_KEYS,
        thermal: THERMAL_KEYS,
        electrical: ELECTRICAL_KEYS,
        dielectric: super::coupled::DIELECTRIC_KEYS,
        superconducting: super::bcs::SC_KEYS,
        piezoelectric: super::coupled::PIEZOELECTRIC_KEYS,
        photoelastic: super::coupled::PHOTOELASTIC_KEYS,
        electro_optic: super::coupled::ELECTRO_OPTIC_KEYS,
        viscoelastic: super::coupled::VISCOELASTIC_KEYS,
        interfaces: super::coupled::INTERFACES_KEYS,
        interface_layer: super::coupled::IFACE_KEYS,
        optical: optics::OPTICAL_KEYS,
        model: optics::MODEL_KEYS,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linbo3_anisotropy_and_stiffness_duals_use_textbook_voigt() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
        let lib = Library::load_dir(&dir).unwrap();
        let ln = lib.get("lithium_niobate").unwrap();
        let an = ln.anisotropy([0.0, 0.0, 0.0]);
        let e = an.piezo_e.expect("LiNbO3 piezo_e");
        // Weis & Gaylord (1985) 3m symmetry in textbook Voigt (xx, yy, zz, yz, xz, xy):
        // e15 == e24 == 3.7638, e16 == -e22, e14 == 0
        assert!((e[0][4] - 3.7638).abs() < 1e-3, "e15 = {}", e[0][4]);
        assert!((e[0][4] - e[1][3]).abs() < 1e-12, "e15 vs e24: {} vs {}", e[0][4], e[1][3]);
        assert!((e[0][5] + e[1][1]).abs() < 1e-12, "e16 vs -e22: {} vs {}", e[0][5], -e[1][1]);
        assert!(e[0][3].abs() < 1e-12, "spurious e14 = {}", e[0][3]);

        // Both C^E and C^D are in textbook Voigt (C44 = C[3][3] = C[4][4])
        let ce = ln.stiffness_constant_e().unwrap();
        let cd = ln.stiffness_constant_d().unwrap();
        assert!((ce[3][3] - ce[4][4]).abs() < 1e-3);
        assert!((cd[3][3] - cd[4][4]).abs() < 1e-3);
        assert!(cd[3][3] > 80e9 && cd[3][3] < 100e9, "C44^D = {}", cd[3][3]);

        // Viscosity stored in textbook Voigt rotates once via rotate_m6_global -> solver Voigt (xx,yy,zz,xy,yz,xz)
        let aln = lib.get("aluminum_nitride").unwrap();
        let eta_tb = aln.anisotropy([0.0, 0.0, 0.0]).viscosity.unwrap();
        let eta_sol = super::super::rotate_m6_global(&eta_tb, [0.0, 0.0, 0.0]);
        // AlN: eta11=0.80 mPa*s, eta12=0.25 mPa*s -> eta66 (xy) = 0.275 mPa*s, eta44=eta55 (yz,xz) = 0.25 mPa*s
        assert!((eta_sol[3][3] - 0.275e-3).abs() < 1e-12, "solver xy eta = {}", eta_sol[3][3]);
        assert!((eta_sol[4][4] - 0.25e-3).abs() < 1e-12, "solver yz eta = {}", eta_sol[4][4]);
        assert!((eta_sol[5][5] - 0.25e-3).abs() < 1e-12, "solver xz eta = {}", eta_sol[5][5]);
    }
}
