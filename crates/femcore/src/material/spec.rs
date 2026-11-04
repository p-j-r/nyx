//! A material as a human-written TOML file: `MaterialSpec` (density, elastic stiffness, thermal and electrical tensors, optical
//! dispersion), reading with plain-language errors, writing back, and `Library`, a folder of such files loaded at run time.
//! The format is version `FORMAT` (key `format = 1`); the schema is in `schema.rs` and documented in `docs/MATERIALS.md`.

use super::elastic::Elastic;
use super::optics::{self, Optical, OpticsError};
use super::reader::{MaterialError, Section};
use super::source::{read_sources, write_sources, Source};
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
    #[serde(default)]
    pub source: Vec<Source>,
    /// W/(m K)
    pub conductivity: Option<TensorProp>,
    /// J/(kg K)
    pub specific_heat: Option<f64>,
    /// 1/K
    pub expansion: Option<TensorProp>,
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
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MaterialSpec {
    pub name: String,
    pub description: Option<String>,
    /// references the numbers come from
    pub source: Vec<Source>,
    pub notes: Option<String>,
    /// kg/m^3
    pub density: Option<f64>,
    pub elastic: Option<Elastic>,
    pub thermal: Thermal,
    pub electrical: Electrical,
    pub optical: Option<Optical>,
}

const TOP_KEYS: &[&str] = &["format", "name", "description", "source", "notes", "extends", "density", "density_unit", "elastic", "thermal", "electrical", "optical"];
const ELASTIC_KEYS: &[&str] = &["source", "unit", "compliance_unit", "class", "unique_axis", "matrix", "compliance", "E", "nu", "lambda", "mu", "K", "G", "c11", "c12", "c44"];
const THERMAL_KEYS: &[&str] = &["source", "conductivity", "conductivity_unit", "specific_heat", "specific_heat_unit", "expansion", "expansion_unit"];
const ELECTRICAL_KEYS: &[&str] = &["source", "permittivity", "permeability", "conductivity", "conductivity_unit"];

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

        match s.value("format") {
            Some(toml::Value::Integer(FORMAT)) => {}
            Some(toml::Value::Integer(n)) => return Err(s.err("format", format!("this file is format {n}, but this version of nyx reads format {FORMAT}"))),
            Some(_) => return Err(s.err("format", "must be a whole number; write format = 1")),
            None => return Err(s.err("format", "is missing; the first line of a material file is format = 1 (the version of this file format)")),
        }
        s.string("name")?;
        let description = s.string("description")?;
        let notes = s.string("notes")?;
        let source = read_sources(&mut s, "source")?;
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
            thermal.source = read_sources(&mut t, "source")?;
            thermal.conductivity = TensorProp::read(&mut t, "conductivity", Dim::ThermalConductivity, Requirement::PositiveDefinite, "thermal conductivity")?;
            let cp = t.quantity("specific_heat", Dim::SpecificHeat, None)?;
            thermal.specific_heat = positive(&t, "specific_heat", cp, "the specific heat")?;
            thermal.expansion = TensorProp::read(&mut t, "expansion", Dim::PerKelvin, Requirement::SymmetricOnly, "thermal expansion")?;
            t.finish(THERMAL_KEYS)?;
        }

        let mut electrical = Electrical::default();
        if let Some(mut t) = s.sub("electrical")? {
            electrical.source = read_sources(&mut t, "source")?;
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
        if let Some(v) = write_sources(&self.source) {
            t.insert("source".into(), v);
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
        if let Some(v) = write_sources(&self.thermal.source) {
            th.insert("source".into(), v);
        }
        put_tensor(&mut th, "conductivity", self.thermal.conductivity.as_ref(), Some("W/(m*K)"));
        if let Some(c) = self.thermal.specific_heat {
            th.insert("specific_heat".into(), format!("{c} J/(kg*K)").into());
        }
        put_tensor(&mut th, "expansion", self.thermal.expansion.as_ref(), Some("1/K"));
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
        MaterialSpec::from_table(&merged, fallback_name)
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

/// `over` laid over `base`: every top-level key of `over` replaces the one of `base`, and a table (`elastic`, `thermal`, `electrical`,
/// `optical`) is replaced as a whole, never mixed key by key, so a half-replaced stiffness cannot occur. `name` is not inherited; a
/// unit key (`density_unit`) of the base does not follow a value that replaces the base's (units are never assumed).
pub fn merge_over(base: &toml::Table, over: &toml::Table) -> toml::Table {
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
    r
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
