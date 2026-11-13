//! The materials of a session: the library of material files (`materials/*.toml`, ADR 29 to 35) and the materials a project carries as TOML
//! text, which may extend a library material. A problem names a material of either (`MaterialAssign::Named` or `Oriented`); the engine
//! hands the merged library to the physics through `Ctx` and works out the wavelength an optical model is evaluated at from the study.

use crate::engine::Engine;
use femcore::material::{Library, MaterialSpec};
use nyx_physics::{Error, Physics, ProblemSpec};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// The folders searched for material files, first match wins: `NYX_MATERIALS` (a path list), `materials` next to the program, and a
/// `materials` folder of any parent of the program's folder or of the current folder that holds `.toml` files (a build tree).
pub fn default_dirs() -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = vec![];
    if let Some(p) = std::env::var_os("NYX_MATERIALS") {
        v.extend(std::env::split_paths(&p));
    }
    let has_files = |d: &Path| {
        std::fs::read_dir(d)
            .map(|r| r.filter_map(|e| e.ok()).any(|e| e.path().extension().and_then(|x| x.to_str()) == Some("toml")))
            .unwrap_or(false)
    };
    let mut starts: Vec<PathBuf> = vec![];
    if let Ok(exe) = std::env::current_exe() {
        if let Some(d) = exe.parent() {
            starts.push(d.to_path_buf());
        }
    }
    if let Ok(cwd) = std::env::current_dir() {
        starts.push(cwd);
    }
    for s in starts {
        let mut d: Option<&Path> = Some(&s);
        while let Some(dir) = d {
            let m = dir.join("materials");
            if has_files(&m) && !v.contains(&m) {
                v.push(m);
                break;
            }
            d = dir.parent();
        }
    }
    v
}

/// Loads the first folder that has material files; an empty library (and the reason) when none loads.
pub fn load_library(dirs: &[PathBuf]) -> (Library, Vec<String>) {
    let mut notes = vec![];
    for d in dirs {
        match Library::load_dir(d) {
            Ok(l) => return (l, notes),
            Err(errs) => notes.extend(errs.iter().map(|e| format!("{}: {e}", d.display()))),
        }
    }
    if dirs.is_empty() {
        notes.push(
            "no material folder found (set NYX_MATERIALS, or keep a `materials` folder next to the program)".into(),
        );
    }
    (Library::default(), notes)
}

/// The frequency (Hz) of a study, from the first parameter of the study descriptor with the unit Hz (given, else its default); the free-space
/// wavelength follows. `None` for a study without a frequency.
pub fn study_wavelength(p: &dyn Physics, spec: &ProblemSpec) -> Option<f64> {
    use nyx_physics::ParamKind;
    let studies = p.studies();
    let sd = studies.iter().find(|s| s.kind == spec.study.kind)?;
    let d = sd.params.iter().find(|d| d.unit == "Hz")?;
    let f = spec.study.params.get(&d.name).and_then(Value::as_f64).or(match d.kind {
        ParamKind::Number { default, .. } => Some(default),
        _ => None,
    })?;
    (f > 0.0 && f.is_finite()).then(|| 299_792_458.0 / f)
}

fn describe(m: &MaterialSpec) -> Value {
    json!({
        "name": m.name,
        "density": m.density,
        "elastic": m.elastic.as_ref().map(|e| json!({ "class": format!("{:?}", e.detect_class()) })),
        "optical": m.optical.is_some(),
        "thermal": m.thermal.conductivity.is_some() || m.thermal.specific_heat.is_some() || m.thermal.expansion.is_some(),
        "electrical": m.electrical.permittivity.is_some() || m.electrical.permeability.is_some() || m.electrical.conductivity.is_some(),
        "toml": m.to_toml(),
    })
}

impl Engine {
    /// The library a solve sees: the folder's materials and the project's own, which may extend them.
    pub fn library(&self) -> Result<Arc<Library>, Error> {
        if self.project.materials.is_empty() {
            return Ok(self.base_library.clone());
        }
        let lib = self.base_library.with_texts(&self.project.materials).map_err(|errs| {
            Error::invalid(format!(
                "the materials of the project do not read: {}",
                errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; ")
            ))
        })?;
        Ok(Arc::new(lib))
    }

    /// Replaces the folder library (tests, or a front-end that lets the person pick another folder).
    pub fn set_base_library(&mut self, library: Library) {
        self.base_library = Arc::new(library);
    }

    pub fn list_materials(&self) -> Result<Value, Error> {
        let lib = self.library()?;
        Ok(json!(lib
            .names()
            .iter()
            .filter_map(|n| lib.get(n))
            .map(|m| json!({
                "name": m.name,
                "in_project": self.project.materials.keys().any(|k| k.eq_ignore_ascii_case(&m.name)),
                "optical": m.optical.is_some(),
                "elastic": m.elastic.is_some(),
            }))
            .collect::<Vec<_>>()))
    }

    pub fn describe_material(&self, name: &str) -> Result<Value, Error> {
        let lib = self.library()?;
        let m = lib.require(name).map_err(|e| Error::invalid(e.to_string()))?;
        Ok(describe(m))
    }

    /// Adds (or replaces) a material the project carries, given as TOML text. It is read now, so a mistake is reported here.
    pub fn add_material(&mut self, name: &str, toml_text: &str) -> Result<Value, Error> {
        if name.trim().is_empty() || name.contains(['/', '\\', '*', ',']) {
            return Err(Error::invalid("a material name is plain text without slashes, commas or stars"));
        }
        let mut texts = self.project.materials.clone();
        texts.insert(name.to_string(), toml_text.to_string());
        let lib = self
            .base_library
            .with_texts(&texts)
            .map_err(|errs| Error::invalid(errs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join("; ")))?;
        let m = lib.get(name).expect("just added");
        let v = describe(m);
        self.project.materials = texts;
        Ok(v)
    }

    pub fn remove_material(&mut self, name: &str) -> Result<(), Error> {
        if self.project.materials.remove(name).is_none() {
            return Err(Error::invalid(format!("the project carries no material '{name}'")));
        }
        Ok(())
    }
}
