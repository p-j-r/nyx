use crate::sanitize::{contains_denied, is_excluded_material_file, sanitize_material_text};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;
use toml::Value;

const CAPABILITY_KEYS: &[&str] = &[
    "elastic",
    "thermal",
    "dielectric",
    "electrical",
    "optical",
    "piezoelectric",
    "photoelastic",
    "electro_optic",
    "viscoelastic",
    "superconducting",
    "sheet_2d",
    "nonlinear",
    "interfaces",
];

#[derive(Debug, Serialize)]
pub struct MaterialSource {
    pub cite: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub doi: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verified: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct MaterialEntry {
    pub id: String,
    pub name: String,
    pub kind: String,
    pub point_group: Option<String>,
    pub reference_temperature: Option<String>,
    pub density: Option<String>,
    pub description: String,
    pub capabilities: Vec<String>,
    pub sources: Vec<MaterialSource>,
}

#[derive(Debug, Serialize)]
pub struct MaterialsCatalog {
    pub count: usize,
    pub kinds: Vec<String>,
    pub materials: Vec<MaterialEntry>,
}

pub fn generate_materials_json(repo: &Path, out_site: &Path) -> Result<usize, String> {
    let mat_dir = repo.join("materials");
    let entries = fs::read_dir(&mat_dir)
        .map_err(|e| format!("failed to read {}: {e}", mat_dir.display()))?;

    let mut raw_tables: BTreeMap<String, toml::map::Map<String, Value>> = BTreeMap::new();
    for entry in entries.filter_map(Result::ok) {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let Some(fname) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if !fname.ends_with(".toml") || is_excluded_material_file(fname) {
            continue;
        }
        let stem = fname.trim_end_matches(".toml").to_string();
        let text = fs::read_to_string(&path)
            .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
        let val: Value = text
            .parse::<Value>()
            .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;
        let table = val
            .as_table()
            .cloned()
            .ok_or_else(|| format!("{} root is not a TOML table", path.display()))?;
        raw_tables.insert(stem, table);
    }

    let mut materials = Vec::with_capacity(raw_tables.len());
    let mut kinds_set = BTreeSet::new();

    for id in raw_tables.keys() {
        let resolved = resolve_extends(id, &raw_tables, 0)?;
        let name = resolved
            .get("name")
            .and_then(Value::as_str)
            .unwrap_or(id.as_str())
            .to_string();
        let kind = resolved
            .get("kind")
            .and_then(Value::as_str)
            .unwrap_or("unspecified")
            .to_string();
        kinds_set.insert(kind.clone());

        let point_group = resolved
            .get("point_group")
            .and_then(Value::as_str)
            .map(str::to_string);
        let reference_temperature = resolved
            .get("reference_temperature")
            .and_then(Value::as_str)
            .map(str::to_string);
        let density = resolved
            .get("density")
            .and_then(Value::as_str)
            .map(str::to_string);

        let raw_desc = resolved
            .get("description")
            .and_then(Value::as_str)
            .unwrap_or("");
        let description = sanitize_material_text(raw_desc);

        let mut capabilities = Vec::new();
        for &cap in CAPABILITY_KEYS {
            if resolved.contains_key(cap) {
                capabilities.push(cap.to_string());
            }
        }

        let sources = collect_sources(&resolved);

        // Fail closed if any denied term survived into the material entry
        let entry = MaterialEntry {
            id: id.clone(),
            name,
            kind,
            point_group,
            reference_temperature,
            density,
            description,
            capabilities,
            sources,
        };
        let check_json = serde_json::to_string(&entry)
            .map_err(|e| format!("failed to serialize material '{id}': {e}"))?;
        if let Some(pat) = contains_denied(&check_json) {
            return Err(format!(
                "material '{id}' retained denied term '{pat}' after sanitisation"
            ));
        }
        materials.push(entry);
    }

    let catalog = MaterialsCatalog {
        count: materials.len(),
        kinds: kinds_set.into_iter().collect(),
        materials,
    };

    let data_dir = out_site.join("data");
    fs::create_dir_all(&data_dir)
        .map_err(|e| format!("failed to create {}: {e}", data_dir.display()))?;
    let out_path = data_dir.join("materials.json");
    let json = serde_json::to_string_pretty(&catalog)
        .map_err(|e| format!("failed to serialize materials.json: {e}"))?;
    fs::write(&out_path, format!("{json}\n"))
        .map_err(|e| format!("failed to write {}: {e}", out_path.display()))?;

    Ok(catalog.count)
}

fn resolve_extends(
    id: &str,
    tables: &BTreeMap<String, toml::map::Map<String, Value>>,
    depth: usize,
) -> Result<toml::map::Map<String, Value>, String> {
    if depth > 8 {
        return Err(format!("material inheritance cycle at '{id}'"));
    }
    let child = tables
        .get(id)
        .ok_or_else(|| format!("missing material '{id}'"))?;
    if let Some(parent_id) = child.get("extends").and_then(Value::as_str) {
        if !tables.contains_key(parent_id) {
            return Err(format!(
                "material '{id}' extends excluded or missing parent '{parent_id}'"
            ));
        }
        let mut merged = resolve_extends(parent_id, tables, depth + 1)?;
        // Per ADR 34: child replaces whole tables / top-level keys of parent.
        for (k, v) in child {
            if k == "extends" {
                continue;
            }
            merged.insert(k.clone(), v.clone());
        }
        Ok(merged)
    } else {
        Ok(child.clone())
    }
}

fn collect_sources(table: &toml::map::Map<String, Value>) -> Vec<MaterialSource> {
    let mut out = Vec::new();
    let mut seen_cites = BTreeSet::new();

    if let Some(src) = table.get("source") {
        push_source_value(src, &mut out, &mut seen_cites);
    }
    for (k, v) in table {
        if k == "source" {
            continue;
        }
        if let Some(sub) = v.as_table() {
            if let Some(src) = sub.get("source") {
                push_source_value(src, &mut out, &mut seen_cites);
            }
            for sub_v in sub.values() {
                if let Some(nested) = sub_v.as_table() {
                    if let Some(src) = nested.get("source") {
                        push_source_value(src, &mut out, &mut seen_cites);
                    }
                }
            }
        }
    }
    out
}

fn push_source_value(
    val: &Value,
    out: &mut Vec<MaterialSource>,
    seen: &mut BTreeSet<String>,
) {
    match val {
        Value::String(s) => {
            let cite = sanitize_material_text(s);
            if !cite.is_empty() && seen.insert(cite.clone()) {
                out.push(MaterialSource {
                    cite,
                    doi: None,
                    url: None,
                    verified: None,
                });
            }
        }
        Value::Table(t) => {
            if let Some(raw_cite) = t.get("cite").and_then(Value::as_str) {
                let cite = sanitize_material_text(raw_cite);
                if !cite.is_empty() && seen.insert(cite.clone()) {
                    out.push(MaterialSource {
                        cite,
                        doi: t.get("doi").and_then(Value::as_str).map(str::to_string),
                        url: t.get("url").and_then(Value::as_str).map(str::to_string),
                        verified: t.get("verified").and_then(Value::as_str).map(str::to_string),
                    });
                }
            }
        }
        Value::Array(arr) => {
            for item in arr {
                push_source_value(item, out, seen);
            }
        }
        _ => {}
    }
}
