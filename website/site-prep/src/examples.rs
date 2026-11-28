use crate::sanitize::{contains_denied, is_excluded_example_file};
use serde::Serialize;
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

#[derive(Debug, Serialize)]
pub struct ExampleEntry {
    pub id: String,
    pub slug: String,
    pub title: String,
    pub category: String,
    pub physics: String,
    pub study: String,
    pub description: String,
    pub file: String,
    pub command: String,
    pub domains: Vec<String>,
    pub plate: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct ExamplesCatalog {
    pub count: usize,
    pub categories: Vec<String>,
    pub examples: Vec<ExampleEntry>,
}

pub fn generate_examples_json(repo: &Path, out_site: &Path) -> Result<usize, String> {
    let ex_dir = repo.join("examples");
    let entries = fs::read_dir(&ex_dir)
        .map_err(|e| format!("failed to read {}: {e}", ex_dir.display()))?;

    let mut files: Vec<_> = entries
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.ends_with(".nyx.json") && !is_excluded_example_file(n))
        })
        .collect();
    files.sort();

    let mut examples = Vec::with_capacity(files.len());
    let mut categories_set = BTreeSet::new();

    for path in files {
        let fname = path
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| format!("invalid filename {}", path.display()))?;
        let id = fname.trim_end_matches(".nyx.json").to_string();
        let slug = id.replace('_', "-");

        let raw = fs::read_to_string(&path)
            .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
        let doc: Value = serde_json::from_str(&raw)
            .map_err(|e| format!("failed to parse {}: {e}", path.display()))?;

        let title = doc
            .get("title")
            .and_then(Value::as_str)
            .unwrap_or(&id)
            .to_string();
        let dim = doc
            .get("geometry")
            .and_then(|g| g.get("dim"))
            .and_then(Value::as_u64)
            .unwrap_or(2);

        let physics = doc
            .get("active")
            .and_then(Value::as_str)
            .unwrap_or("ems")
            .to_string();

        let phys_list = doc.get("physics").and_then(Value::as_array);
        let active_phys_obj = phys_list.and_then(|arr| {
            arr.iter()
                .find(|p| p.get("physics").and_then(Value::as_str) == Some(physics.as_str()))
                .or_else(|| arr.first())
        });

        let study = active_phys_obj
            .and_then(|p| p.get("study"))
            .and_then(|s| s.get("kind"))
            .and_then(Value::as_str)
            .unwrap_or("unknown")
            .to_string();

        let category = match physics.as_str() {
            "dfs" => "Diffusion & Thermal",
            "sms" => "Solid Mechanics",
            "ems" => "Electromagnetics",
            _ => "Multiphysics",
        }
        .to_string();
        categories_set.insert(category.clone());

        let mut domains: Vec<String> = active_phys_obj
            .and_then(|p| p.get("materials"))
            .and_then(Value::as_object)
            .map(|m| m.keys().cloned().collect())
            .unwrap_or_default();
        domains.sort();

        let all_physics_codes: Vec<&str> = phys_list
            .map(|arr| {
                arr.iter()
                    .filter_map(|p| p.get("physics").and_then(Value::as_str))
                    .collect()
            })
            .unwrap_or_default();
        let bc_count = active_phys_obj
            .and_then(|p| p.get("bcs"))
            .and_then(Value::as_array)
            .map(Vec::len)
            .unwrap_or(0);
        let order = active_phys_obj
            .and_then(|p| p.get("order"))
            .and_then(Value::as_u64);

        let description = build_description(
            dim,
            &physics,
            &study,
            &all_physics_codes,
            &domains,
            bc_count,
            order,
        );

        let plate = match id.as_str() {
            "15_soi_wire" => Some("figures/soi_wire_modes.png".to_string()),
            "18_wgm_cylindrical_pml" => Some("figures/wgm_cyl_pml.png".to_string()),
            "11_anchor_loss" => Some("figures/anchor_loss_modes.png".to_string()),
            "20_cavity_3d_spherical_pml" => Some("figures/cavity_3d_cutaway.png".to_string()),
            _ => None,
        };

        let entry = ExampleEntry {
            id: id.clone(),
            slug,
            title,
            category,
            physics,
            study,
            description,
            file: format!("examples/{fname}"),
            command: format!("nyx --batch examples/{fname}"),
            domains,
            plate,
        };
        let serialized = serde_json::to_string(&entry)
            .map_err(|e| format!("failed to serialize example '{id}': {e}"))?;
        if let Some(pat) = contains_denied(&serialized) {
            return Err(format!("example '{id}' contains denied term '{pat}'"));
        }
        examples.push(entry);
    }

    let catalog = ExamplesCatalog {
        count: examples.len(),
        categories: categories_set.into_iter().collect(),
        examples,
    };

    let data_dir = out_site.join("data");
    fs::create_dir_all(&data_dir)
        .map_err(|e| format!("failed to create {}: {e}", data_dir.display()))?;
    let out_path = data_dir.join("examples.json");
    let json = serde_json::to_string_pretty(&catalog)
        .map_err(|e| format!("failed to serialize examples.json: {e}"))?;
    fs::write(&out_path, format!("{json}\n"))
        .map_err(|e| format!("failed to write {}: {e}", out_path.display()))?;

    Ok(catalog.count)
}

fn build_description(
    dim: u64,
    physics: &str,
    study: &str,
    all_physics: &[&str],
    domains: &[String],
    bc_count: usize,
    order: Option<u64>,
) -> String {
    let coupled = if all_physics.len() > 1 {
        format!(" (coupled: {})", all_physics.join(" + "))
    } else {
        String::new()
    };
    let order_note = order.map(|o| format!(", P{o} elements")).unwrap_or_default();
    let domain_note = if domains.is_empty() {
        String::new()
    } else {
        format!(" across regions [{}]", domains.join(", "))
    };
    let bc_note = if bc_count == 1 {
        " with 1 boundary condition".to_string()
    } else if bc_count > 1 {
        format!(" with {bc_count} boundary conditions")
    } else {
        String::new()
    };
    format!("{dim}D {physics} `{study}` study{coupled}{order_note}{domain_note}{bc_note}.")
}
