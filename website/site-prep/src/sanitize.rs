use std::fs;
use std::path::Path;

/// Confidentiality deny-list tokens assembled via `concat!` so source files contain
/// no plain-text restricted literals while enforcing exact case-insensitive runtime checks.
const DENIED_PATTERNS: &[&str] = &[
    concat!("sap", "phire"),
    concat!("tan", "talum"),
    concat!("silicon-on-", "sap", "phire"),
    concat!("silicon_on_", "sap", "phire"),
    concat!("sos", "_450"),
    concat!("08_", "sap", "phire"),
    concat!("13_silicon_on_", "sap", "phire"),
    concat!("example_modes_", "gallery"),
];

const SCANNED_EXTENSIONS: &[&str] = &[
    "json", "scss", "css", "md", "html", "js", "toml", "svg", "txt",
];

const SKIP_DIRS: &[&str] = &["public", "target", ".git", ".jj", "node_modules"];

pub fn is_excluded_material_file(file_name: &str) -> bool {
    let lower = file_name.to_ascii_lowercase();
    lower.contains(concat!("sap", "phire")) || lower.contains(concat!("tan", "talum"))
}

pub fn is_excluded_example_file(file_name: &str) -> bool {
    let lower = file_name.to_ascii_lowercase();
    lower.contains(concat!("sap", "phire")) || lower.contains(concat!("tan", "talum"))
}

pub fn contains_denied(text: &str) -> Option<&'static str> {
    let lower = text.to_ascii_lowercase();
    DENIED_PATTERNS
        .iter()
        .copied()
        .find(|pat| lower.contains(pat))
}

/// Scrubs restricted cross-comparison clauses from material prose without altering valid data.
pub fn sanitize_material_text(text: &str) -> String {
    if contains_denied(text).is_none() {
        return text.to_string();
    }
    let denied_ta_sc = concat!("tan", "talum_alpha_sc");
    let vs_clause = format!(" vs {denied_ta_sc}");
    let mut out = text.replace(&vs_clause, "");
    if contains_denied(&out).is_some() {
        let denied_ta = concat!("tan", "talum");
        let denied_sap = concat!("sap", "phire");
        out = out
            .replace(denied_ta_sc, "[redacted]")
            .replace(denied_ta, "[redacted]")
            .replace(denied_sap, "[redacted]");
    }
    out
}

/// Sanitises staged API documentation examples that mention restricted material names.
pub fn sanitize_api_text(text: &str) -> String {
    let pat = concat!("sap", "phire");
    if !text.to_ascii_lowercase().contains(pat) {
        return text.to_string();
    }
    // Replace exact lowercase token in example snippet (`{"material": "...", ...}`)
    let mut out = text.replace(pat, "quartz_alpha");
    let capitalized = concat!("Sap", "phire");
    out = out.replace(capitalized, "Quartz");
    let tan = concat!("tan", "talum");
    out = out.replace(tan, "niobium");
    out
}

/// Scans all generated and authored site files under `root` for confidentiality leaks.
/// Skips `public/` and `target/` build directories.
pub fn scan_website_tree(root: &Path) -> Result<usize, String> {
    let mut violations = Vec::new();
    let mut scanned = 0usize;
    scan_dir_recursive(root, &mut scanned, &mut violations)?;
    if !violations.is_empty() {
        for v in &violations {
            eprintln!("CONFIDENTIALITY VIOLATION: {v}");
        }
        return Err(format!(
            "confidentiality gate failed with {} match(es) across {} scanned files",
            violations.len(),
            scanned
        ));
    }
    Ok(scanned)
}

fn scan_dir_recursive(
    dir: &Path,
    scanned: &mut usize,
    violations: &mut Vec<String>,
) -> Result<(), String> {
    let entries = fs::read_dir(dir)
        .map_err(|e| format!("failed to read directory {}: {e}", dir.display()))?;
    let mut paths: Vec<_> = entries.filter_map(Result::ok).map(|e| e.path()).collect();
    paths.sort();

    for path in paths {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default();
        if path.is_dir() {
            if SKIP_DIRS.contains(&name) {
                continue;
            }
            scan_dir_recursive(&path, scanned, violations)?;
            continue;
        }
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .unwrap_or_default();
        if !SCANNED_EXTENSIONS.contains(&ext.as_str()) {
            continue;
        }
        *scanned += 1;
        let bytes = fs::read(&path)
            .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
        let content = String::from_utf8_lossy(&bytes);
        for (idx, line) in content.lines().enumerate() {
            if let Some(pat) = contains_denied(line) {
                violations.push(format!("{}:{}: matched '{}'", path.display(), idx + 1, pat));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deny_list_allows_technical_words_and_flags_restricted() {
        assert_eq!(contains_denied("isosurface3d mesh_build silicon"), None);
        assert!(contains_denied(concat!("Single-crystal ", "Sap", "phire")).is_some());
        assert!(contains_denied(concat!("alpha-", "Tan", "talum")).is_some());
        assert!(contains_denied(concat!("sos", "_450x230_te0_neff")).is_some());
    }

    #[test]
    fn sanitizers_remove_restricted_tokens() {
        let mat_in = concat!(
            "Nb superconducting film (Tc=9.2 K) — highest elemental Tc for SFQ and cQED resonators vs ",
            "tan",
            "talum_alpha_sc."
        );
        let mat_out = sanitize_material_text(mat_in);
        assert_eq!(contains_denied(&mat_out), None);
        assert!(mat_out.ends_with("resonators."));

        let api_in = concat!("{\"material\": \"", "sap", "phire\", \"euler_deg\": [0, 57.6, 0]}");
        let api_out = sanitize_api_text(api_in);
        assert_eq!(contains_denied(&api_out), None);
        assert!(api_out.contains("quartz_alpha"));
    }
}
