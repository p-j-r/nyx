use crate::sanitize::{contains_denied, sanitize_api_text};
use serde::Serialize;
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

const EXPECTED_COMMAND_COUNT: usize = 114;
const EXPECTED_CATALOG_HASH: &str = "b3:60697136d596a5b047f70971b85b0f54";

#[derive(Debug, Serialize)]
pub struct ApiArgument {
    pub name: String,
    #[serde(rename = "type")]
    pub arg_type: String,
    pub required: bool,
    pub default: Option<String>,
    pub description: String,
}

#[derive(Debug, Serialize)]
pub struct ApiCommand {
    pub name: String,
    pub category: String,
    pub access: String,
    pub summary: String,
    pub description: String,
    pub signature: String,
    pub arguments: Vec<ApiArgument>,
    pub example: Option<String>,
    pub doc_file: String,
}

#[derive(Debug, Serialize)]
pub struct ApiCatalog {
    pub command_count: usize,
    pub catalog_hash: String,
    pub categories: Vec<String>,
    pub commands: Vec<ApiCommand>,
}

pub fn generate_api_json(repo: &Path, out_site: &Path) -> Result<(usize, String), String> {
    let api_dir = repo.join("docs/api");
    let readme_path = api_dir.join("README.md");
    let readme = fs::read_to_string(&readme_path)
        .map_err(|e| format!("failed to read {}: {e}", readme_path.display()))?;

    let catalog_hash = parse_catalog_hash(&readme)?;
    if catalog_hash != EXPECTED_CATALOG_HASH {
        return Err(format!(
            "docs/api/README.md catalog hash mismatch: got '{catalog_hash}', expected '{EXPECTED_CATALOG_HASH}'"
        ));
    }

    let entries = fs::read_dir(&api_dir)
        .map_err(|e| format!("failed to read {}: {e}", api_dir.display()))?;
    let mut md_files: Vec<_> = entries
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.extension().and_then(|e| e.to_str()) == Some("md")
                && p.file_name().and_then(|n| n.to_str()) != Some("README.md")
        })
        .collect();
    md_files.sort();

    let mut commands = Vec::with_capacity(md_files.len());
    let mut categories_set = BTreeSet::new();

    for path in md_files {
        let fname = path
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| format!("invalid filename {}", path.display()))?;
        let raw = fs::read_to_string(&path)
            .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
        let cmd = parse_command_page(fname, &raw)?;
        categories_set.insert(cmd.category.clone());
        commands.push(cmd);
    }

    if commands.len() != EXPECTED_COMMAND_COUNT {
        return Err(format!(
            "expected {EXPECTED_COMMAND_COUNT} commands in docs/api/*.md, found {}",
            commands.len()
        ));
    }

    let catalog = ApiCatalog {
        command_count: commands.len(),
        catalog_hash: catalog_hash.clone(),
        categories: categories_set.into_iter().collect(),
        commands,
    };

    let json = serde_json::to_string_pretty(&catalog)
        .map_err(|e| format!("failed to serialize api.json: {e}"))?;
    if let Some(pat) = contains_denied(&json) {
        return Err(format!(
            "api.json retained denied term '{pat}' after sanitisation"
        ));
    }

    let data_dir = out_site.join("data");
    fs::create_dir_all(&data_dir)
        .map_err(|e| format!("failed to create {}: {e}", data_dir.display()))?;
    let out_path = data_dir.join("api.json");
    fs::write(&out_path, format!("{json}\n"))
        .map_err(|e| format!("failed to write {}: {e}", out_path.display()))?;

    Ok((catalog.command_count, catalog_hash))
}

fn parse_catalog_hash(readme: &str) -> Result<String, String> {
    for line in readme.lines() {
        if let Some(rest) = line.strip_prefix("Catalog hash: `") {
            if let Some((hash, _)) = rest.split_once('`') {
                return Ok(hash.to_string());
            }
        }
    }
    Err("could not locate 'Catalog hash:' in docs/api/README.md".to_string())
}

fn parse_command_page(fname: &str, raw: &str) -> Result<ApiCommand, String> {
    let lines: Vec<&str> = raw.lines().collect();
    let name = lines
        .first()
        .and_then(|l| l.strip_prefix("# `"))
        .and_then(|l| l.strip_suffix('`'))
        .ok_or_else(|| format!("{fname}: missing '# `<name>`' header"))?
        .to_string();

    let access_line = lines
        .iter()
        .find(|l| l.starts_with("**") && l.contains("access `"))
        .ok_or_else(|| format!("{fname}: missing title/access line"))?;

    let (title_part, access_part) = access_line
        .split_once(" · access `")
        .ok_or_else(|| format!("{fname}: malformed access line"))?;
    let summary = sanitize_api_text(title_part.trim_matches('*').trim());
    let access = access_part
        .split('`')
        .next()
        .unwrap_or("read")
        .to_string();

    let args_idx = lines
        .iter()
        .position(|l| l.trim() == "## Arguments")
        .ok_or_else(|| format!("{fname}: missing '## Arguments'"))?;
    let access_idx = lines
        .iter()
        .position(|l| *l == *access_line)
        .unwrap_or(1);

    let desc_lines: Vec<String> = lines[access_idx + 1..args_idx]
        .iter()
        .map(|l| sanitize_api_text(l.trim()))
        .filter(|l| !l.is_empty())
        .collect();
    let description = desc_lines.join(" ");

    let example_idx = lines
        .iter()
        .position(|l| l.trim() == "## Example")
        .unwrap_or(lines.len());

    let mut arguments = Vec::new();
    for line in &lines[args_idx + 1..example_idx] {
        let trimmed = line.trim();
        if !trimmed.starts_with("| `") {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('|').collect();
        // | `name` | type | required | default | doc | -> 7 parts with leading/trailing empty
        if cols.len() < 6 {
            continue;
        }
        let arg_name = cols[1].trim().trim_matches('`').to_string();
        let arg_type = cols[2].trim().replace("\\|", "|");
        let required = cols[3].trim() == "yes";
        let raw_default = cols[4].trim().trim_matches('`');
        let default = if raw_default.is_empty() {
            None
        } else {
            Some(sanitize_api_text(raw_default))
        };
        let doc = sanitize_api_text(cols[5].trim());
        arguments.push(ApiArgument {
            name: arg_name,
            arg_type,
            required,
            default,
            description: doc,
        });
    }

    let signature = format_signature(&name, &arguments);

    let mut example = None;
    let mut in_json = false;
    for line in &lines[example_idx..] {
        let trimmed = line.trim();
        if trimmed == "```json" {
            in_json = true;
            continue;
        }
        if in_json {
            if trimmed == "```" {
                break;
            }
            if !trimmed.is_empty() {
                example = Some(sanitize_api_text(trimmed));
            }
        }
    }

    let category = categorize_command(&name).to_string();

    Ok(ApiCommand {
        name: name.clone(),
        category,
        access,
        summary,
        description,
        signature,
        arguments,
        example,
        doc_file: format!("docs/api/{fname}"),
    })
}

fn format_signature(name: &str, args: &[ApiArgument]) -> String {
    // Only include top-level arguments (exclude nested `foo.bar` / `foo[].bar` rows)
    let top: Vec<&ApiArgument> = args
        .iter()
        .filter(|a| !a.name.contains('.') && !a.name.contains('['))
        .collect();
    if top.is_empty() {
        return format!("{name}()");
    }
    let req: Vec<String> = top
        .iter()
        .filter(|a| a.required)
        .map(|a| format!("{}: {}", a.name, a.arg_type))
        .collect();
    let opt: Vec<String> = top
        .iter()
        .filter(|a| !a.required)
        .map(|a| format!("[{}: {}]", a.name, a.arg_type))
        .collect();
    let mut parts = req;
    parts.extend(opt);
    format!("{name}({})", parts.join(", "))
}

fn categorize_command(name: &str) -> &'static str {
    match name {
        "add_absorbing_layer"
        | "add_box"
        | "add_chip_stack"
        | "add_circle"
        | "add_cone"
        | "add_cut"
        | "add_cylinder"
        | "add_difference"
        | "add_ellipse"
        | "add_extrude"
        | "add_intersect"
        | "add_overlay"
        | "add_pattern"
        | "add_periodic"
        | "add_pml"
        | "add_polygon"
        | "add_polygons"
        | "add_rect"
        | "add_revolve"
        | "add_rounded_rect"
        | "add_sheet"
        | "add_size"
        | "add_sphere"
        | "add_tag"
        | "add_torus"
        | "add_transform"
        | "add_union"
        | "geometry_names"
        | "import_gds"
        | "list_operations"
        | "mesh_build"
        | "pml_diagnostics"
        | "remove_feature"
        | "remove_param"
        | "rename_feature"
        | "set_geometry"
        | "set_mesh"
        | "set_param"
        | "suppress_feature"
        | "tag_feature"
        | "update_feature" => "Geometry & Meshing",

        "add_condition"
        | "add_material"
        | "add_physics"
        | "add_pml_conditions"
        | "bcs_surface_impedance"
        | "describe_material"
        | "describe_physics"
        | "evaluate_material"
        | "fetch_material"
        | "import_material"
        | "list_materials"
        | "list_physics"
        | "list_remote_materials"
        | "remove_condition"
        | "remove_material"
        | "remove_physics"
        | "search_remote_materials"
        | "set_active"
        | "set_material"
        | "set_physics_params"
        | "set_problem"
        | "set_study" => "Physics & Materials",

        "convergence_check" | "solve_run" | "sweep_cancel" | "sweep_run" | "sweep_table" => {
            "Solving & Sweeps"
        }

        "analyse"
        | "confinement"
        | "extremum"
        | "field"
        | "field_op"
        | "figure"
        | "integrate"
        | "line_cut"
        | "list_analyses"
        | "overlap"
        | "probe"
        | "result_field"
        | "result_info"
        | "vvuq_export" => "Analyses & Fields",

        "check_project"
        | "export_results"
        | "open_results"
        | "project_info"
        | "project_json"
        | "project_load"
        | "project_new"
        | "project_open"
        | "project_save"
        | "project_script"
        | "results_file_info"
        | "store_open"
        | "store_runs"
        | "store_save" => "Projects & Storage",

        "job_cancel" | "job_list" | "job_status" | "job_wait" | "poll_jobs" | "queue_cancel"
        | "queue_get" | "queue_hold" | "queue_list" | "queue_pause" | "queue_release"
        | "queue_resume" | "queue_set_priority" | "queue_submit" | "queue_wait" => "Jobs & Queue",

        "get_guide" | "list_commands" | "list_features" => "Discovery & Reference",

        _ => "Engine",
    }
}
