//! A material as a human-written TOML file: `MaterialSpec` (density, elastic stiffness, thermal and electrical tensors, optical
//! dispersion), reading with plain-language errors, writing back, and `Library`, a folder of such files loaded at run time.
//! The format is version `FORMAT` (key `format = 1`); the schema is in `schema.rs` and documented in `docs/MATERIALS.md`.

use super::elastic::Elastic;
use super::optics::{self, Optical, OpticsError};
use super::reader::{MaterialError, Section};
use super::tensor::{Requirement, TensorProp};
use super::units::Dim;
use super::{Anisotropy, Stiffness};
use num_complex::Complex64;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

/// The version of the material file format this code reads and writes.
pub const FORMAT: i64 = 1;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Thermal {
    /// W/(m K)
    pub conductivity: Option<TensorProp>,
    /// J/(kg K)
    pub specific_heat: Option<f64>,
    /// 1/K
    pub expansion: Option<TensorProp>,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Electrical {
    /// relative, at low frequency
    pub permittivity: Option<TensorProp>,
    /// relative
    pub permeability: Option<f64>,
    /// S/m
    pub conductivity: Option<TensorProp>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MaterialSpec {
    pub name: String,
    pub description: Option<String>,
    /// references the numbers come from
    pub source: Vec<String>,
    pub notes: Option<String>,
    /// kg/m^3
    pub density: Option<f64>,
    pub elastic: Option<Elastic>,
    pub thermal: Thermal,
    pub electrical: Electrical,
    pub optical: Option<Optical>,
}

const TOP_KEYS: &[&str] = &["format", "name", "description", "source", "notes", "density", "density_unit", "elastic", "thermal", "electrical", "optical"];
const ELASTIC_KEYS: &[&str] = &["unit", "compliance_unit", "class", "unique_axis", "matrix", "compliance", "E", "nu", "lambda", "mu", "K", "G", "c11", "c12", "c44"];
const THERMAL_KEYS: &[&str] = &["conductivity", "conductivity_unit", "specific_heat", "specific_heat_unit", "expansion", "expansion_unit"];
const ELECTRICAL_KEYS: &[&str] = &["permittivity", "permeability", "conductivity", "conductivity_unit"];

fn positive(s: &Section, key: &str, v: Option<f64>, what: &str) -> Result<Option<f64>, MaterialError> {
    match v {
        Some(x) if x <= 0.0 => Err(s.err(key, format!("{what} must be positive (found {x})"))),
        other => Ok(other),
    }
}

impl MaterialSpec {
    /// Reads a material from TOML text. `fallback_name` is used when the file has no `name`.
    pub fn from_toml_named(text: &str, fallback_name: &str) -> Result<MaterialSpec, MaterialError> {
        let table: toml::Table = text.parse().map_err(|e: toml::de::Error| MaterialError {
            material: fallback_name.to_string(),
            path: String::new(),
            message: format!("the file is not valid TOML: {e}"),
        })?;
        let name = table.get("name").and_then(|v| v.as_str()).unwrap_or(fallback_name).to_string();
        let mut s = Section::new(&name, "", &table);

        match s.value("format") {
            Some(toml::Value::Integer(FORMAT)) => {}
            Some(toml::Value::Integer(n)) => return Err(s.err("format", format!("this file is format {n}, but this version of nyx reads format {FORMAT}"))),
            Some(_) => return Err(s.err("format", "must be a whole number; write format = 1")),
            None => return Err(s.err("format", "is missing; the first line of a material file is format = 1 (the version of this file format)")),
        }
        s.string("name")?;
        let description = s.string("description")?;
        let notes = s.string("notes")?;
        let source = match s.value("source") {
            None => Vec::new(),
            Some(toml::Value::String(x)) => vec![x.clone()],
            Some(toml::Value::Array(a)) => a
                .iter()
                .map(|v| v.as_str().map(str::to_string).ok_or_else(|| s.err("source", "must be text, or a list of texts")))
                .collect::<Result<_, _>>()?,
            Some(_) => return Err(s.err("source", "must be text, or a list of texts")),
        };
        let density = s.quantity("density", Dim::Density, None)?;
        let density = positive(&s, "density", density, "the density")?;

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
            thermal.conductivity = TensorProp::read(&mut t, "conductivity", Dim::ThermalConductivity, Requirement::PositiveDefinite, "thermal conductivity")?;
            let cp = t.quantity("specific_heat", Dim::SpecificHeat, None)?;
            thermal.specific_heat = positive(&t, "specific_heat", cp, "the specific heat")?;
            thermal.expansion = TensorProp::read(&mut t, "expansion", Dim::PerKelvin, Requirement::SymmetricOnly, "thermal expansion")?;
            t.finish(THERMAL_KEYS)?;
        }

        let mut electrical = Electrical::default();
        if let Some(mut t) = s.sub("electrical")? {
            electrical.permittivity = TensorProp::read(&mut t, "permittivity", Dim::Dimensionless, Requirement::PositiveDefinite, "relative permittivity")?;
            let mu = t.number("permeability")?;
            electrical.permeability = positive(&t, "permeability", mu, "the relative permeability")?;
            electrical.conductivity = TensorProp::read(&mut t, "conductivity", Dim::ElectricalConductivity, Requirement::NonNegative, "electrical conductivity")?;
            t.finish(ELECTRICAL_KEYS)?;
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
        s.finish(TOP_KEYS)?;
        Ok(MaterialSpec { name, description, source, notes, density, elastic, thermal, electrical, optical })
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

    /// The material as TOML text, in the same format (comments of an original file are not kept).
    pub fn to_toml(&self) -> String {
        let mut t = toml::Table::new();
        t.insert("format".into(), toml::Value::Integer(FORMAT));
        t.insert("name".into(), self.name.clone().into());
        if let Some(d) = &self.description {
            t.insert("description".into(), d.clone().into());
        }
        match self.source.len() {
            0 => {}
            1 => {
                t.insert("source".into(), self.source[0].clone().into());
            }
            _ => {
                t.insert("source".into(), toml::Value::Array(self.source.iter().map(|x| x.clone().into()).collect()));
            }
        }
        if let Some(n) = &self.notes {
            t.insert("notes".into(), n.clone().into());
        }
        if let Some(d) = self.density {
            t.insert("density".into(), format!("{d} kg/m^3").into());
        }
        if let Some(e) = &self.elastic {
            t.insert("elastic".into(), toml::Value::Table(e.write()));
        }
        let mut th = toml::Table::new();
        put_tensor(&mut th, "conductivity", self.thermal.conductivity.as_ref(), Some("W/(m*K)"));
        if let Some(c) = self.thermal.specific_heat {
            th.insert("specific_heat".into(), format!("{c} J/(kg*K)").into());
        }
        put_tensor(&mut th, "expansion", self.thermal.expansion.as_ref(), Some("1/K"));
        if !th.is_empty() {
            t.insert("thermal".into(), toml::Value::Table(th));
        }
        let mut el = toml::Table::new();
        put_tensor(&mut el, "permittivity", self.electrical.permittivity.as_ref(), None);
        if let Some(m) = self.electrical.permeability {
            el.insert("permeability".into(), toml::Value::Float(m));
        }
        put_tensor(&mut el, "conductivity", self.electrical.conductivity.as_ref(), Some("S/m"));
        if !el.is_empty() {
            t.insert("electrical".into(), toml::Value::Table(el));
        }
        if let Some(o) = &self.optical {
            t.insert("optical".into(), toml::Value::Table(o.write()));
        }
        toml::to_string(&t).expect("a table of plain values always serialises")
    }

    // ----- use by the solvers -----

    /// The material in the older `Anisotropy` form (stiffness as a general matrix, tensors), placed with Bunge angles in degrees.
    pub fn anisotropy(&self, euler_deg: [f64; 3]) -> Anisotropy {
        Anisotropy {
            stiffness: self.elastic.as_ref().map(|e| Stiffness::General(e.c)),
            conductivity: self.thermal.conductivity.as_ref().map(TensorProp::to_legacy),
            permittivity: self.electrical.permittivity.as_ref().map(TensorProp::to_legacy),
            expansion: self.thermal.expansion.as_ref().map(TensorProp::to_legacy),
            euler_deg,
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
        put("eps_r", self.electrical.permittivity.as_ref().and_then(TensorProp::as_scalar));
        put("k", self.thermal.conductivity.as_ref().and_then(TensorProp::as_scalar));
        put("alpha", self.thermal.expansion.as_ref().and_then(TensorProp::as_scalar));
        if let Some(e) = &self.elastic {
            if e.detect_class() == super::Class::Isotropic {
                let (lam, mu) = (e.c[0][1], 0.5 * (e.c[0][0] - e.c[0][1]));
                put("E", Some(mu * (3.0 * lam + 2.0 * mu) / (lam + mu)));
                put("nu", Some(lam / (2.0 * (lam + mu))));
            }
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
#[derive(Clone, Debug, Default)]
pub struct Library {
    materials: BTreeMap<String, MaterialSpec>,
}

impl Library {
    /// Loads every `*.toml` file of the folder (not its subfolders). All files with problems are reported together.
    pub fn load_dir(dir: &Path) -> Result<Library, Vec<MaterialError>> {
        let entries = std::fs::read_dir(dir).map_err(|e| {
            vec![MaterialError { material: dir.display().to_string(), path: String::new(), message: format!("cannot read the material folder: {e}") }]
        })?;
        let mut files: Vec<_> = entries.filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.extension().and_then(|x| x.to_str()) == Some("toml")).collect();
        files.sort();
        let (mut materials, mut errors) = (BTreeMap::new(), Vec::new());
        for f in files {
            let key = f.file_stem().and_then(|s| s.to_str()).unwrap_or_default().to_ascii_lowercase();
            match MaterialSpec::load(&f) {
                Ok(m) => {
                    materials.insert(key, m);
                }
                Err(e) => errors.push(e),
            }
        }
        if errors.is_empty() {
            Ok(Library { materials })
        } else {
            Err(errors)
        }
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

/// The keys each table of a material file accepts (the contract test compares them with the JSON Schema).
pub struct Keys {
    pub top: &'static [&'static str],
    pub elastic: &'static [&'static str],
    pub thermal: &'static [&'static str],
    pub electrical: &'static [&'static str],
    pub optical: &'static [&'static str],
    pub model: &'static [&'static str],
}

pub fn keys() -> Keys {
    Keys { top: TOP_KEYS, elastic: ELASTIC_KEYS, thermal: THERMAL_KEYS, electrical: ELECTRICAL_KEYS, optical: optics::OPTICAL_KEYS, model: optics::MODEL_KEYS }
}
