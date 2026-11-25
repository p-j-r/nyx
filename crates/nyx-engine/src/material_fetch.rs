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

fn is_installed_bundle_dir(dir: &Path) -> bool {
    let in_app_bundle = dir.ancestors().any(|a| {
        a.file_name().and_then(|n| n.to_str()) == Some("Contents")
            && a.parent().and_then(|p| p.extension()).and_then(|e| e.to_str()) == Some("app")
    });
    if in_app_bundle {
        return true;
    }
    dir.join("install-manifest.json").is_file()
        || dir.parent().is_some_and(|p| p.join("install-manifest.json").is_file())
        || dir.parent().and_then(Path::parent).is_some_and(|p| p.join("install-manifest.json").is_file())
}

fn dir_is_writable(dir: &Path) -> bool {
    let probe = dir.join(format!(".nyx-write-probe-{}", std::process::id()));
    match std::fs::OpenOptions::new().write(true).create_new(true).open(&probe) {
        Ok(f) => {
            drop(f);
            let _ = std::fs::remove_file(&probe);
            true
        }
        Err(_) => false,
    }
}

/// First writable `materials/` directory for CLI `--save` or GUI save-to-library toggles.
/// Falls back to [`crate::materials::user_materials_dir`] (`<workspace>/materials`) when the bundled
/// library sits inside an installed package/bundle or is read-only.
pub fn default_materials_folder() -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("NYX_MATERIALS") {
        if let Some(first) = std::env::split_paths(&p).next() {
            if first.is_dir() && dir_is_writable(&first) {
                return Some(first);
            }
        }
    }
    if let Some(first) = crate::materials::default_dirs().into_iter().next() {
        if !is_installed_bundle_dir(&first) && dir_is_writable(&first) {
            return Some(first);
        }
    }
    crate::materials::user_materials_dir()
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

const WINDOWS_RESERVED_NAMES: &[&str] = &[
    "con", "prn", "aux", "nul", "com1", "com2", "com3", "com4", "com5", "com6", "com7", "com8", "com9", "lpt1", "lpt2",
    "lpt3", "lpt4", "lpt5", "lpt6", "lpt7", "lpt8", "lpt9",
];

const RAW_GITHUB_PREFIXES: &[&str] = &["/polyanskiy/refractiveindex.info-database/", "/bjornsturmberg/NumBAT/"];

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedUrl {
    pub url: String,
    pub host: String,
}

/// Validates that `url` is a clean `https://` URL targeting an allowlisted remote material host (no userinfo, custom port, backslash, or path traversal).
pub fn validate_remote_url(url: &str) -> Result<ValidatedUrl, Error> {
    let u = url.trim();
    if u.is_empty() || u.chars().any(|c| c.is_ascii_control() || c.is_ascii_whitespace()) {
        return Err(Error::invalid(format!(
            "invalid remote URL '{u}': empty or contains whitespace/control characters"
        )));
    }
    let Some(rest) = u.strip_prefix("https://") else {
        return Err(Error::invalid(format!("only https:// URLs may be downloaded, got '{u}'")));
    };
    let auth_end = rest.find(['/', '?', '#']).unwrap_or(rest.len());
    let authority = &rest[..auth_end];
    let suffix = &rest[auth_end..];
    if authority.is_empty()
        || authority.contains('@')
        || authority.contains('\\')
        || authority.contains(':')
        || authority.ends_with('.')
        || suffix.starts_with('#')
    {
        return Err(Error::invalid(format!("invalid or unsafe URL authority in '{u}'")));
    }
    let host = authority.to_ascii_lowercase();
    let allowed_host = matches!(
        host.as_str(),
        "api.materialsproject.org"
            | "aflowlib.org"
            | "aflow.org"
            | "aflowlib.duke.edu"
            | "refractiveindex.info"
            | "raw.githubusercontent.com"
    );
    if !allowed_host {
        return Err(Error::invalid(format!("remote URL host '{host}' is not on the HTTPS allowlist")));
    }
    if host == "raw.githubusercontent.com" {
        let path = suffix.split(['?', '#']).next().unwrap_or("");
        let lower_path = path.to_ascii_lowercase();
        if path.contains('\\')
            || path.split('/').any(|seg| seg == ".." || seg == ".")
            || lower_path.contains("%2e")
            || lower_path.contains("%2f")
            || lower_path.contains("%5c")
            || !RAW_GITHUB_PREFIXES.iter().any(|p| path.starts_with(p))
        {
            return Err(Error::invalid(format!(
                "raw.githubusercontent.com URL path '{path}' is not pinned to a curated material repository"
            )));
        }
    }
    Ok(ValidatedUrl { url: u.to_string(), host })
}

fn build_curl_args(url: &str, with_stdin_config: bool) -> Vec<String> {
    let mut args = vec![
        "--fail".to_string(),
        "--silent".to_string(),
        "--show-error".to_string(),
        "--proto".to_string(),
        "=https".to_string(),
        "--max-redirs".to_string(),
        "0".to_string(),
        "--max-time".to_string(),
        "20".to_string(),
        "--max-filesize".to_string(),
        "5242880".to_string(),
    ];
    if with_stdin_config {
        args.push("--config".to_string());
        args.push("-".to_string());
    }
    args.push("--".to_string());
    args.push(url.to_string());
    args
}

/// Downloads a UTF-8 text payload over HTTPS using system `curl` (`--proto =https --max-redirs 0 --max-time 20`, 5 MiB cap).
pub fn fetch_url_text(url: &str) -> Result<String, Error> {
    fetch_url_text_with_api_key(url, None)
}

/// Downloads a UTF-8 text payload over HTTPS using system `curl` (`--proto =https --max-redirs 0 --max-time 20`, 5 MiB cap),
/// passing an optional `X-API-KEY` header (only for `api.materialsproject.org`) via `curl --config -` on stdin so the key never appears in `argv` or logs.
pub fn fetch_url_text_with_api_key(url: &str, api_key: Option<&str>) -> Result<String, Error> {
    use std::io::Write;
    use std::process::Stdio;
    let parsed = validate_remote_url(url)?;
    let u = parsed.url.as_str();
    let bin = if cfg!(windows) { "curl.exe" } else { "curl" };
    let key_clean = api_key.map(str::trim).filter(|s| !s.is_empty());
    if let Some(k) = key_clean {
        if parsed.host != "api.materialsproject.org" {
            return Err(Error::invalid(format!(
                "API key header X-API-KEY may only be sent to api.materialsproject.org, not '{}'",
                parsed.host
            )));
        }
        if k.contains(['\r', '\n']) {
            return Err(Error::invalid("API key must not contain newline characters"));
        }
    }
    let args = build_curl_args(u, key_clean.is_some());
    let mut cmd = Command::new(bin);
    cmd.args(&args).stdout(Stdio::piped()).stderr(Stdio::piped());
    if key_clean.is_some() {
        cmd.stdin(Stdio::piped());
    }
    let mut child = cmd.spawn().map_err(|e| {
        Error::invalid(format!(
            "could not run {bin} to download '{u}': {e} (download the file manually and use `nyx material import`)"
        ))
    })?;
    if let Some(k) = key_clean {
        if let Some(mut stdin) = child.stdin.take() {
            let escaped = k.replace('\\', "\\\\").replace('"', "\\\"");
            let _ = writeln!(stdin, "header = \"X-API-KEY: {escaped}\"");
        }
    }
    let out = child
        .wait_with_output()
        .map_err(|e| Error::invalid(format!("could not wait on {bin} while downloading '{u}': {e}")))?;
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

pub fn validate_material_name(name: &str) -> Result<&str, Error> {
    let n = name.trim();
    if n.is_empty() || n.len() > 64 || !n.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-') {
        return Err(Error::invalid(
            "a material name must be 1..=64 ASCII characters [A-Za-z0-9_-]+ (no ':', path separators, or '..')",
        ));
    }
    if WINDOWS_RESERVED_NAMES.contains(&n.to_ascii_lowercase().as_str()) {
        return Err(Error::invalid(format!(
            "'{n}' is a reserved Windows device name and cannot be used as a material name"
        )));
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
        self.fetch_material_with_key(source, target, name, description, extends, add_to_project, folder, None)
    }

    /// Downloads a remote database entry or catalog preset over HTTPS (using `api_key` or `MP_API_KEY` via `curl --config -` stdin for Materials Project),
    /// converts it into validated TOML, and adds/saves as requested.
    #[allow(clippy::too_many_arguments)]
    pub fn fetch_material_with_key(
        &mut self,
        source: &str,
        target: Option<&str>,
        name: Option<&str>,
        description: Option<&str>,
        extends: Option<&str>,
        add_to_project: bool,
        folder: Option<&Path>,
        api_key: Option<&str>,
    ) -> Result<Value, Error> {
        if self.path_root().is_some() {
            return Err(Error::invalid(
                "fetch_material is refused in sandboxed script mode (path_root is set); use `import_material` instead",
            ));
        }
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
        let env_mp_key = std::env::var("MP_API_KEY").ok();
        let api_key_clean = api_key.map(str::trim).filter(|s| !s.is_empty());
        let env_key_clean = env_mp_key.as_deref().map(str::trim).filter(|s| !s.is_empty());
        let any_key_present = api_key_clean.is_some() || env_key_clean.is_some();
        let target_is_url = resolved_target.contains("://");
        if target_is_url {
            if fmt == ImportFormat::MaterialsProject {
                return Err(Error::invalid(
                    "arbitrary URL targets are not allowed for source 'mp'; provide a Materials Project id such as 'mp-149'",
                ));
            }
            if any_key_present {
                return Err(Error::invalid(
                    "arbitrary URL targets are refused when an API key (`api_key` or `MP_API_KEY`) is present",
                ));
            }
        }
        let effective_key = if fmt == ImportFormat::MaterialsProject { api_key_clean.or(env_key_clean) } else { None };
        let url = if fmt == ImportFormat::MaterialsProject {
            if effective_key.is_none() {
                return Err(Error::invalid(
                    "Materials Project live download requires an API key (`api_key` argument, `--api-key` CLI flag, or `MP_API_KEY` environment variable)",
                ));
            }
            let mp_id = resolved_target.trim();
            if mp_id.is_empty() || !mp_id.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_') {
                return Err(Error::invalid(format!("invalid Materials Project id '{mp_id}'")));
            }
            format!(
                "https://api.materialsproject.org/materials/summary/?material_ids={mp_id}&_fields=material_id,formula_pretty,symmetry,density,elasticity,dielectric,piezoelectric"
            )
        } else {
            resolve_url(fmt, &resolved_target).map_err(|e| Error::invalid(e.to_string()))?
        };
        let validated = validate_remote_url(&url)?;
        let raw = fetch_url_text_with_api_key(&validated.url, effective_key)?;
        let opts = ImportOptions {
            name: &final_name,
            description: resolved_desc,
            extends: resolved_extends,
            source_url: Some(&validated.url),
        };
        let toml_text = import_to_toml(fmt, &raw, &opts).map_err(|e| Error::invalid(e.to_string()))?;
        let mut extra = vec![("format", json!(fmt.id())), ("fetched_url", json!(validated.url))];
        if let Some(id) = preset_id {
            extra.push(("preset", json!(id)));
        }
        finalize_import(self, &final_name, &toml_text, add_to_project, folder, &extra)
    }

    /// Searches the curated `REMOTE_CATALOG` and (when `live` is true) live AFLOW and Materials Project HTTPS endpoints.
    pub fn search_remote_materials(
        &self,
        query: &str,
        domain: Option<&str>,
        source: Option<&str>,
        live: bool,
        api_key: Option<&str>,
    ) -> Result<Value, Error> {
        if live && self.path_root().is_some() {
            return Err(Error::invalid(
                "search_remote_materials with live=true is refused in sandboxed script mode (path_root is set)",
            ));
        }
        let q = query.trim();
        let q_low = q.to_ascii_lowercase();
        let dom_filter = domain.map(str::trim).filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("all"));
        let src_filter = source.map(str::trim).filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("all"));

        let catalog_matches: Vec<Value> = remote_catalog()
            .into_iter()
            .filter(|c| dom_filter.is_none_or(|d| c.domain.eq_ignore_ascii_case(d)))
            .filter(|c| src_filter.is_none_or(|s| c.format.eq_ignore_ascii_case(s)))
            .filter(|c| {
                q_low.is_empty()
                    || c.id.to_ascii_lowercase().contains(&q_low)
                    || c.target.to_ascii_lowercase().contains(&q_low)
                    || c.summary.to_ascii_lowercase().contains(&q_low)
                    || c.extends.is_some_and(|e| e.to_ascii_lowercase().contains(&q_low))
                    || c.domain.to_ascii_lowercase().contains(&q_low)
                    || c.format.to_ascii_lowercase().contains(&q_low)
            })
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

        let mut live_matches: Vec<Value> = Vec::new();
        let mut live_warnings: Vec<String> = Vec::new();

        if live && !q.is_empty() && q.chars().all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_') {
            let want_aflow = src_filter.is_none_or(|s| s.eq_ignore_ascii_case("aflow"));
            if want_aflow {
                let aflow_url =
                    format!("https://aflowlib.duke.edu/search/API/?species({q}),ael_voigt_bulk_modulus_B,$paging(1,8)");
                match fetch_url_text(&aflow_url) {
                    Ok(body) => {
                        if let Ok(Value::Array(arr)) = serde_json::from_str::<Value>(&body) {
                            for item in arr.into_iter().take(8) {
                                let compound = item.get("compound").and_then(Value::as_str).unwrap_or(q);
                                let auid = item.get("auid").and_then(Value::as_str).unwrap_or("");
                                let aurl = item.get("aurl").and_then(Value::as_str).unwrap_or("");
                                let target_url = if !aurl.is_empty() {
                                    format!("https://{}/?format=json", aurl.trim_start_matches("aflowlib.duke.edu:"))
                                } else {
                                    aflow_url.clone()
                                };
                                live_matches.push(json!({
                                    "format": "aflow",
                                    "id": if auid.is_empty() { compound.to_string() } else { format!("{compound}_{auid}") },
                                    "target": target_url,
                                    "summary": format!("AFLOW live DFT entry for {compound} ({auid})"),
                                }));
                            }
                        }
                    }
                    Err(e) => live_warnings.push(format!("AFLOW live search unavailable: {e}")),
                }
            }

            let want_mp = src_filter.is_none_or(|s| s.eq_ignore_ascii_case("mp"));
            let env_mp = std::env::var("MP_API_KEY").ok();
            let mp_key = api_key.map(str::trim).filter(|s| !s.is_empty()).or(env_mp.as_deref());
            if want_mp {
                if let Some(key) = mp_key {
                    let mp_url = if q.starts_with("mp-") {
                        format!(
                            "https://api.materialsproject.org/materials/summary/?material_ids={q}&_fields=material_id,formula_pretty,symmetry,density&_limit=8"
                        )
                    } else {
                        format!(
                            "https://api.materialsproject.org/materials/summary/?formula={q}&_fields=material_id,formula_pretty,symmetry,density&_limit=8"
                        )
                    };
                    match fetch_url_text_with_api_key(&mp_url, Some(key)) {
                        Ok(body) => {
                            if let Ok(v) = serde_json::from_str::<Value>(&body) {
                                if let Some(arr) = v.get("data").and_then(Value::as_array) {
                                    for item in arr.iter().take(8) {
                                        let mid = item.get("material_id").and_then(Value::as_str).unwrap_or(q);
                                        let formula = item.get("formula_pretty").and_then(Value::as_str).unwrap_or(q);
                                        let pg = item
                                            .get("symmetry")
                                            .and_then(|s| s.get("point_group"))
                                            .and_then(Value::as_str)
                                            .unwrap_or("");
                                        live_matches.push(json!({
                                            "format": "mp",
                                            "id": format!("mp_{}_{}", formula.to_ascii_lowercase(), mid.replace('-', "_")),
                                            "target": mid,
                                            "summary": format!("Materials Project {formula} ({mid}, point group {pg})"),
                                        }));
                                    }
                                }
                            }
                        }
                        Err(e) => live_warnings.push(format!("Materials Project live search failed: {e}")),
                    }
                } else if src_filter.is_some_and(|s| s.eq_ignore_ascii_case("mp")) {
                    live_warnings
                        .push("Materials Project live search skipped: provide `api_key` or set `MP_API_KEY`".into());
                }
            }
        }

        Ok(json!({
            "query": q,
            "catalog_matches": catalog_matches,
            "live_matches": live_matches,
            "live_warnings": live_warnings,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allowlist_rejects_untrusted_hosts_and_unsafe_urls() {
        assert!(validate_remote_url("https://api.materialsproject.org/materials/summary/?material_ids=mp-149").is_ok());
        assert!(validate_remote_url("https://aflowlib.org/API/aflux/?species(Si)").is_ok());
        assert!(validate_remote_url("https://aflow.org/API/aflux/?species(Si)").is_ok());
        assert!(validate_remote_url("https://aflowlib.duke.edu/search/API/?species(Si)").is_ok());
        assert!(validate_remote_url("https://refractiveindex.info/data/main/Si/nk/Green-2008.yml").is_ok());
        assert!(validate_remote_url(
            "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/main/database/data/main/Si/nk/Green-2008.yml"
        )
        .is_ok());
        assert!(validate_remote_url(
            "https://raw.githubusercontent.com/bjornsturmberg/NumBAT/master/backend/material_data/Si_2021_Poulton.json"
        )
        .is_ok());

        for bad in [
            "http://api.materialsproject.org/materials/summary/",
            "https://evil.example/steal",
            "https://127.0.0.1/admin",
            "https://169.254.169.254/latest/meta-data/",
            "https://api.materialsproject.org@evil.example/",
            "https://api.materialsproject.org:8443/materials/summary/",
            "https://api.materialsproject.org\\@evil.example/",
            "https://raw.githubusercontent.com/evil/repo/main/payload.yml",
            "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/../../evil/repo/secret.yml",
            "https://raw.githubusercontent.com/polyanskiy/refractiveindex.info-database/%2e%2e/secret.yml",
        ] {
            assert!(validate_remote_url(bad).is_err(), "expected rejection for {bad}");
        }
    }

    #[test]
    fn curl_args_disable_redirects_and_force_https() {
        let args = build_curl_args("https://aflow.org/API/aflux/?species(Si)", false);
        assert!(!args.iter().any(|a| a == "--location" || a == "-L"), "{args:?}");
        assert!(args.windows(2).any(|w| w == ["--proto", "=https"]), "{args:?}");
        assert!(args.windows(2).any(|w| w == ["--max-redirs", "0"]), "{args:?}");

        let mp_args = build_curl_args("https://api.materialsproject.org/materials/summary/", true);
        assert!(mp_args.windows(2).any(|w| w == ["--max-redirs", "0"]), "{mp_args:?}");
        assert!(mp_args.windows(2).any(|w| w == ["--config", "-"]), "{mp_args:?}");
    }

    #[test]
    fn api_key_is_never_sent_to_non_mp_hosts() {
        let err = fetch_url_text_with_api_key("https://aflow.org/API/aflux/?species(Si)", Some("secret-token"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("api.materialsproject.org"), "{err}");
    }

    #[test]
    fn validate_material_name_rejects_reserved_windows_names_and_punctuation() {
        assert_eq!(validate_material_name("silicon_4k-custom").unwrap(), "silicon_4k-custom");
        for bad in [
            "",
            "foo:bar",
            "foo/bar",
            "foo\\bar",
            "..",
            "a.b",
            "con",
            "CON",
            "Prn",
            "aux",
            "NUL",
            "com1",
            "COM9",
            "lpt1",
            "LPT9",
            &"a".repeat(65),
        ] {
            assert!(validate_material_name(bad).is_err(), "expected rejection for {bad:?}");
        }
    }
}
