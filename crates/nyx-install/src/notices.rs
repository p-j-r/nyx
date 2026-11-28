// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Generates `THIRD_PARTY_NOTICES.txt` offline from `cargo metadata` and validates SPDX licence expressions.

use anyhow::{bail, Context, Result};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

pub const CORRESPONDING_SOURCE_URL: &str = "https://github.com/p-j-r/nyx";

const KNOWN_SPDX_IDS: &[&str] = &[
    "0BSD",
    "AGPL-3.0",
    "AGPL-3.0-only",
    "AGPL-3.0-or-later",
    "Apache-2.0",
    "BSD-2-Clause",
    "BSD-3-Clause",
    "BSL-1.0",
    "BlueOak-1.0.0",
    "CC0-1.0",
    "CDLA-Permissive-2.0",
    "GPL-2.0",
    "GPL-2.0-only",
    "GPL-2.0-or-later",
    "GPL-3.0-only",
    "GPL-3.0-or-later",
    "ISC",
    "LGPL-2.1-only",
    "LGPL-3.0-only",
    "LicenseRef-nyx-Commercial",
    "MIT",
    "MIT-0",
    "MPL-2.0",
    "NCSA",
    "OFL-1.1",
    "OpenSSL",
    "Python-2.0.1",
    "Ubuntu-font-1.0",
    "Unicode-3.0",
    "Unicode-DFS-2016",
    "Unlicense",
    "Zlib",
    "bzip2-1.0.6",
];

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CrateLicense {
    pub name: String,
    pub version: String,
    pub license: String,
    pub repository: String,
}

/// Parses an SPDX expression with `OR` / `AND` / `/` / `WITH` / parentheses and checks every identifier against `KNOWN_SPDX_IDS`.
pub fn validate_spdx_expression(expr: &str) -> Result<Vec<String>> {
    let normalized = expr.replace(['(', ')'], " ").replace('/', " OR ");
    let tokens: Vec<&str> = normalized.split_whitespace().collect();
    if tokens.is_empty() {
        bail!("empty SPDX license expression");
    }
    let mut ids = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let tok = tokens[i];
        match tok {
            "OR" | "AND" => {
                i += 1;
                continue;
            }
            "WITH" => {
                // Skip the following exception token (e.g. LLVM-exception).
                i += 2;
                continue;
            }
            id => {
                let clean = id.trim_end_matches('+');
                if !KNOWN_SPDX_IDS.contains(&clean) && !clean.starts_with("LicenseRef-") {
                    bail!("unrecognised SPDX license token '{id}' in expression '{expr}'");
                }
                ids.push(id.to_string());
                i += 1;
            }
        }
    }
    if ids.is_empty() {
        bail!("no SPDX identifiers in expression '{expr}'");
    }
    Ok(ids)
}

/// Queries `cargo metadata --format-version 1 --offline` (filtered by `target_triple` if given) and collects all third-party packages.
pub fn collect_third_party_licenses(
    repo_root: &Path,
    target_triple: Option<&str>,
) -> Result<Vec<CrateLicense>> {
    let cargo = std::env::var("CARGO").unwrap_or_else(|_| "cargo".to_string());
    let manifest_path = repo_root.join("Cargo.toml");
    let mut cmd = Command::new(&cargo);
    cmd.arg("metadata")
        .arg("--format-version")
        .arg("1")
        .arg("--offline")
        .arg("--manifest-path")
        .arg(&manifest_path);
    if let Some(triple) = target_triple.filter(|s| !s.trim().is_empty()) {
        cmd.arg("--filter-platform").arg(triple);
    }
    let out = cmd
        .output()
        .with_context(|| format!("running cargo metadata in {}", repo_root.display()))?;
    if !out.status.success() {
        bail!(
            "cargo metadata failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let v: Value = serde_json::from_slice(&out.stdout).context("parsing cargo metadata JSON")?;
    let pkgs = v["packages"].as_array().context("missing packages array")?;

    let mut unique: BTreeMap<(String, String), CrateLicense> = BTreeMap::new();
    for p in pkgs {
        // Workspace crates have `source: null`.
        if p.get("source").is_none_or(Value::is_null) {
            continue;
        }
        let name = p["name"].as_str().unwrap_or_default().to_string();
        let version = p["version"].as_str().unwrap_or_default().to_string();
        let license_str = p
            .get("license")
            .and_then(Value::as_str)
            .filter(|s| !s.trim().is_empty())
            .map(str::to_string)
            .or_else(|| {
                p.get("license_file")
                    .and_then(Value::as_str)
                    .filter(|s| !s.trim().is_empty())
                    .map(|f| format!("See bundled license_file: {f}"))
            });
        let Some(license) = license_str else {
            bail!("crate {name} {version} has neither license nor license_file");
        };
        if !license.starts_with("See bundled license_file:") {
            validate_spdx_expression(&license)?;
        }
        let repository = p
            .get("repository")
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_string();
        unique.insert(
            (name.clone(), version.clone()),
            CrateLicense {
                name,
                version,
                license,
                repository,
            },
        );
    }
    Ok(unique.into_values().collect())
}

/// Renders deterministic `THIRD_PARTY_NOTICES.txt` text from resolved crate metadata.
pub fn render_third_party_notices(crates: &[CrateLicense], target_triple: Option<&str>) -> String {
    let mut out = String::new();
    out.push_str("nyx — Third-Party Software Notices and License Summary\n");
    out.push_str("======================================================\n\n");
    out.push_str("Copyright (C) 2026 p-j-r\n");
    out.push_str("SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial\n\n");
    out.push_str("This distribution of nyx is dual licensed under the GNU Affero General Public\n");
    out.push_str(
        "License, version 3 (see LICENSE), or a commercial license available from p-j-r.\n",
    );
    out.push_str(&format!(
        "Corresponding Source (AGPL-3.0 Section 6): {CORRESPONDING_SOURCE_URL}\n"
    ));
    out.push_str(
        "Ports and fonts keeping their own licences (BSD-2-Clause, SIL OFL 1.1) are detailed\n",
    );
    out.push_str("in NOTICE in this directory.\n\n");
    if let Some(triple) = target_triple {
        out.push_str(&format!(
            "Resolved Cargo dependency graph target: {triple}\n"
        ));
    }
    out.push_str(&format!("Third-party crates ({}):\n\n", crates.len()));
    for c in crates {
        if c.repository.is_empty() {
            out.push_str(&format!("- {} {} — {}\n", c.name, c.version, c.license));
        } else {
            out.push_str(&format!(
                "- {} {} — {} ({})\n",
                c.name, c.version, c.license, c.repository
            ));
        }
    }
    out
}

pub fn generate_third_party_notices(
    repo_root: &Path,
    target_triple: Option<&str>,
) -> Result<String> {
    let crates = collect_third_party_licenses(repo_root, target_triple)?;
    Ok(render_third_party_notices(&crates, target_triple))
}
