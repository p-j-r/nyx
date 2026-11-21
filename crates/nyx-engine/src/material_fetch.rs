// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! On-demand material database catalog, pure-Rust conversion and HTTPS download (`refractiveindex.info`, `AFLOW`, `Materials Project`, `NumBAT`, ADR 101).

use crate::engine::Engine;
use crate::materials::describe;
use femcore::material::import::{
    catalog_entry, import_to_toml, resolve_url, validate_against_library, ImportFormat, ImportOptions, REMOTE_CATALOG,
};
use nyx_physics::Error;
use serde_json::{json, Value};
use std::path::{Path, PathBuf};
use std::process::Command;

/// One curated remote preset exposed to front-ends (`nyx` GUI and CLI) without leaking `femcore` types.
#[derive(Clone, Debug)]
pub struct RemoteEntryInfo {
    pub id: &'static str,
    pub domain: &'static str,
    pub format: &'static str,
    pub format_label: &'static str,
    pub target: &'static str,
    pub extends: Option<&'static str>,
    pub summary: &'static str,
    pub url: Option<String>,
}

/// Returns every curated remote preset together with its resolved HTTPS source URL (when directly downloadable without an API key).
pub fn remote_catalog() -> Vec<RemoteEntryInfo> {
    REMOTE_CATALOG
        .iter()
        .map(|c| RemoteEntryInfo {
            id: c.id,
            domain: c.domain,
            format: c.format.id(),
            format_label: c.format.label(),
            target: c.target,
            extends: c.extends,
            summary: c.summary,
            url: resolve_url(c.format, c.target).ok(),
        })
        .collect()
}

/// First writable/discovered `materials/` directory for CLI `--save` or GUI save-to-library toggles.
pub fn default_materials_folder() -> Option<PathBuf> {
    crate::materials::default_dirs().into_iter().next()
}

/// Converts raw external database text into validated `format = 1` material TOML text (`rii`, `aflow`, `mp`, `numbat`).
pub fn convert_to_toml(
    format: &str,
    content: &str,
    name: &str,
    description: Option<&str>,
    extends: Option<&str>,
    source_url: Option<&str>,
) -> Result<String, Error> {
    let fmt = ImportFormat::parse(format).map_err(|e| Error::invalid(e.to_string()))?;
    let opts = ImportOptions { name, description, extends, source_url };
    import_to_toml(fmt, content, &opts).map_err(|e| Error::invalid(e.to_string()))
}

/// Resolves a format + database coordinate into its HTTPS URL (for GUI hints or diagnostics).
pub fn resolve_remote_url(format: &str, target: &str) -> Result<String, Error> {
    let fmt = ImportFormat::parse(format).map_err(|e| Error::invalid(e.to_string()))?;
    resolve_url(fmt, target).map_err(|e| Error::invalid(e.to_string()))
}

/// Downloads a UTF-8 text payload over HTTPS using system `curl` (`--max-time 20`, 5 MiB cap).
pub fn fetch_url_text(url: &str) -> Result<String, Error> {
    let u = url.trim();
    if !u.starts_with("https://") {
        return Err(Error::invalid(format!("only https:// URLs may be downloaded, got '{u}'")));
    }
    let bin = if cfg!(windows) { "curl.exe" } else { "curl" };
    let out = Command::new(bin)
        .args(["--fail", "--silent", "--show-error", "--location", "--max-time", "20", "--max-filesize", "5242880", u])
        .output()
        .map_err(|e| {
            Error::invalid(format!(
                "could not run {bin} to download '{u}': {e} (download the file manually and use `nyx material import`)"
            ))
        })?;
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        let msg = stderr.trim();
        return Err(Error::invalid(if msg.is_empty() {
            format!("download of '{u}' failed with status {}", out.status)
        } else {
            format!("download of '{u}' failed: {msg}")
        }));
    }
    let text = String::from_utf8(out.stdout)
        .map_err(|e| Error::invalid(format!("downloaded payload from '{u}' is not valid UTF-8: {e}")))?;
    if text.trim().is_empty() {
        return Err(Error::invalid(format!("downloaded payload from '{u}' is empty")));
    }
    Ok(text)
}

fn validate_material_name(name: &str) -> Result<&str, Error> {
    let n = name.trim();
    if n.is_empty() || n.contains(['/', '\\', '*', ',']) || n.contains("..") {
        return Err(Error::invalid("a material name must be plain text without '/', '\\', '*', ',' or '..'"));
    }
    Ok(n)
}

/// Derives a safe material filename/identifier stem from a database entry path or query.
pub fn default_name_from_target(target: &str) -> String {
    let tail = target.trim().trim_end_matches('/').rsplit('/').next().unwrap_or(target.trim());
    let stem = tail
        .strip_suffix(".yml")
        .or_else(|| tail.strip_suffix(".yaml"))
        .or_else(|| tail.strip_suffix(".jsonc"))
        .or_else(|| tail.strip_suffix(".json"))
        .unwrap_or(tail);
    let mut slug = String::with_capacity(stem.len());
    let mut prev_us = false;
    for ch in stem.chars() {
        if ch.is_ascii_alphanumeric() {
            slug.push(ch.to_ascii_lowercase());
            prev_us = false;
        } else if !prev_us && !slug.is_empty() {
            slug.push('_');
            prev_us = true;
        }
    }
    while slug.ends_with('_') {
        slug.pop();
    }
    if slug.is_empty() {
        "imported_material".into()
    } else {
        slug
    }
}

fn save_toml_to_folder(folder: &Path, name: &str, toml_text: &str) -> Result<PathBuf, Error> {
    std::fs::create_dir_all(folder)
        .map_err(|e| Error::invalid(format!("cannot create folder '{}': {e}", folder.display())))?;
    let path = folder.join(format!("{name}.toml"));
    std::fs::write(&path, toml_text).map_err(|e| Error::invalid(format!("cannot write '{}': {e}", path.display())))?;
    Ok(path)
}

fn finalize_import(
    engine: &mut Engine,
    name: &str,
    toml_text: &str,
    add_to_project: bool,
    folder: Option<&Path>,
    extra: &[(&str, Value)],
) -> Result<Value, Error> {
    let mut v = if add_to_project {
        engine.add_material(name, toml_text)?
    } else {
        let lib = engine.library()?;
        let spec = validate_against_library(&lib, name, toml_text).map_err(|e| Error::invalid(e.to_string()))?;
        describe(&spec)
    };
    if let Some(obj) = v.as_object_mut() {
        obj.insert("added_to_project".into(), Value::Bool(add_to_project));
        obj.insert("imported_toml".into(), Value::String(toml_text.to_string()));
        if let Some(dir) = folder {
            let saved = save_toml_to_folder(dir, name, toml_text)?;
            obj.insert("saved_path".into(), Value::String(saved.to_string_lossy().into_owned()));
        }
        for (k, val) in extra {
            obj.insert((*k).to_string(), val.clone());
        }
    }
    Ok(v)
}

impl Engine {
    /// Lists curated remote materials in `REMOTE_CATALOG`, optionally filtered by `domain` (`photonics`, `phononics`, `superconducting`, `all`).
    pub fn list_remote_materials(&self, domain: Option<&str>) -> Result<Value, Error> {
        let filter = domain.map(str::trim).filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("all"));
        if let Some(d) = filter {
            if !["photonics", "phononics", "superconducting"].iter().any(|k| k.eq_ignore_ascii_case(d)) {
                return Err(Error::invalid(format!(
                    "unknown domain '{d}'; use photonics, phononics, superconducting or all"
                )));
            }
        }
        let rows: Vec<Value> = remote_catalog()
            .into_iter()
            .filter(|c| filter.is_none_or(|d| c.domain.eq_ignore_ascii_case(d)))
            .map(|c| {
                json!({
                    "id": c.id,
                    "domain": c.domain,
                    "format": c.format,
                    "format_label": c.format_label,
                    "target": c.target,
                    "extends": c.extends,
                    "summary": c.summary,
                    "url": c.url,
                })
            })
            .collect();
        Ok(Value::Array(rows))
    }

    /// Converts raw database text (`rii`, `aflow`, `mp`, `numbat`), validates it against the library and adds/saves as requested.
    #[allow(clippy::too_many_arguments)]
    pub fn import_material(
        &mut self,
        name: &str,
        format: &str,
        content: &str,
        description: Option<&str>,
        extends: Option<&str>,
        add_to_project: bool,
        folder: Option<&Path>,
    ) -> Result<Value, Error> {
        let name = validate_material_name(name)?;
        let fmt = ImportFormat::parse(format).map_err(|e| Error::invalid(e.to_string()))?;
        let opts = ImportOptions { name, description, extends, source_url: None };
        let toml_text = import_to_toml(fmt, content, &opts).map_err(|e| Error::invalid(e.to_string()))?;
        finalize_import(self, name, &toml_text, add_to_project, folder, &[("format", json!(fmt.id()))])
    }

    /// Downloads a remote database entry or catalog preset over HTTPS, converts it into validated TOML, and adds/saves as requested.
    #[allow(clippy::too_many_arguments)]
    pub fn fetch_material(
        &mut self,
        source: &str,
        target: Option<&str>,
        name: Option<&str>,
        description: Option<&str>,
        extends: Option<&str>,
        add_to_project: bool,
        folder: Option<&Path>,
    ) -> Result<Value, Error> {
        let src = source.trim();
        let tgt_opt = target.map(str::trim).filter(|s| !s.is_empty());

        let preset = if src.eq_ignore_ascii_case("catalog") {
            let id =
                tgt_opt.ok_or_else(|| Error::invalid("`target` preset id is required when `source` is 'catalog'"))?;
            Some(catalog_entry(id).ok_or_else(|| Error::invalid(format!("unknown remote preset '{id}'")))?)
        } else if tgt_opt.is_none() {
            catalog_entry(src)
        } else {
            None
        };

        let (fmt, resolved_target, resolved_extends, resolved_desc, default_name, preset_id) = if let Some(entry) =
            preset
        {
            (
                entry.format,
                entry.target.to_string(),
                extends.or(entry.extends),
                description.or(Some(entry.summary)),
                entry.id.to_string(),
                Some(entry.id),
            )
        } else {
            let fmt = ImportFormat::parse(src).map_err(|e| Error::invalid(e.to_string()))?;
            let tgt = tgt_opt.ok_or_else(|| {
                Error::invalid(format!(
                    "`target` is required for source '{src}' (or pass a curated preset id from `list_remote_materials` as `source`)"
                ))
            })?;
            (fmt, tgt.to_string(), extends, description, default_name_from_target(tgt), None)
        };

        let final_name =
            validate_material_name(name.map(str::trim).filter(|s| !s.is_empty()).unwrap_or(&default_name))?.to_string();
        let url = resolve_url(fmt, &resolved_target).map_err(|e| Error::invalid(e.to_string()))?;
        let raw = fetch_url_text(&url)?;
        let opts = ImportOptions {
            name: &final_name,
            description: resolved_desc,
            extends: resolved_extends,
            source_url: Some(&url),
        };
        let toml_text = import_to_toml(fmt, &raw, &opts).map_err(|e| Error::invalid(e.to_string()))?;
        let mut extra = vec![("format", json!(fmt.id())), ("fetched_url", json!(url))];
        if let Some(id) = preset_id {
            extra.push(("preset", json!(id)));
        }
        finalize_import(self, &final_name, &toml_text, add_to_project, folder, &extra)
    }
}
