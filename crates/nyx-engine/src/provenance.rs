// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Deterministic build, runtime, solver-diagnostic, and material provenance helpers for `nyx-engine` (`store_link`, `h5io`, `results_file`).

use crate::project::Project;
use femcore::material::Library;
use nyx_physics::{canonical_json, DiagnosticSeverity, MaterialAssign, SolveResult};
use nyx_store::runs::{BuildProvenance, RunDiagnosticsSummary, RuntimeProvenance};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};

/// Resolved material provenance entry stored in `/provenance/materials_json` inside `results.h5`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MaterialProvenance {
    /// BLAKE3 digest (`b3:...`) of the resolved canonical material TOML (or canonical JSON for inline `Props`).
    pub blake3: String,
    /// Where the material came from (`"library"`, `"project"`, `"inline"`, or `"missing"`).
    pub source: String,
    /// Bibliographic source citations attached to the material specification.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub citations: Vec<String>,
    /// Canonical JSON snapshot of the resolved material specification or inline properties.
    pub params_snapshot: Value,
}

/// Returns compile-time build provenance captured by `crates/nyx-engine/build.rs`.
pub fn build_provenance() -> BuildProvenance {
    BuildProvenance {
        nyx_version: env!("CARGO_PKG_VERSION").to_string(),
        git_commit: env!("NYX_BUILD_GIT_COMMIT").to_string(),
        git_dirty: env!("NYX_BUILD_GIT_DIRTY") == "true",
        rustc_version: env!("NYX_BUILD_RUSTC_VERSION").to_string(),
        target_triple: env!("NYX_BUILD_TARGET").to_string(),
        profile: env!("NYX_BUILD_PROFILE").to_string(),
        cargo_lock_blake3: env!("NYX_BUILD_CARGO_LOCK_BLAKE3").to_string(),
    }
}

/// Resolved `NYX_THREADS` worker pool size used by deterministic solver kernels.
pub fn resolved_nyx_threads() -> usize {
    femcore::par::threads()
}

/// Resolved `NYX_FACTOR_THREADS` pool size used by sparse factorisations.
pub fn resolved_nyx_factor_threads() -> usize {
    femcore::par::factor_threads()
}

/// Returns execution runtime provenance (`created_utc` belongs strictly in `RunRecord`, never in deterministic HDF5).
pub fn runtime_provenance(created_utc: impl Into<String>) -> RuntimeProvenance {
    RuntimeProvenance {
        os: std::env::consts::OS.to_string(),
        cpu_arch: std::env::consts::ARCH.to_string(),
        nyx_threads: resolved_nyx_threads(),
        nyx_factor_threads: resolved_nyx_factor_threads(),
        created_utc: created_utc.into(),
    }
}

/// Condenses per-physics [`nyx_physics::SolverDiagnostics`] into a compact [`RunDiagnosticsSummary`].
pub fn summarize_diagnostics(results: &[SolveResult]) -> Option<RunDiagnosticsSummary> {
    let mut results_with_diagnostics = 0usize;
    let mut all_converged = true;
    let mut info_count = 0usize;
    let mut warning_count = 0usize;
    let mut error_count = 0usize;
    let mut worst_rank = 0u8;
    let mut max_algebraic_residual = 0.0_f64;
    let mut min_rcond_estimate: Option<f64> = None;
    let mut cholesky_fallbacks = 0usize;
    let mut code_set = BTreeSet::new();
    let mut warning_seen = BTreeSet::new();
    let mut warnings = Vec::new();

    for r in results {
        let Some(d) = &r.diagnostics else { continue };
        results_with_diagnostics += 1;
        all_converged &= d.converged;
        if d.max_algebraic_residual.is_finite() {
            max_algebraic_residual = max_algebraic_residual.max(d.max_algebraic_residual);
        }
        if let Some(rc) = d.rcond_estimate.filter(|x| x.is_finite()) {
            min_rcond_estimate = Some(min_rcond_estimate.map_or(rc, |prev| prev.min(rc)));
        }
        cholesky_fallbacks += d.cholesky_fallbacks;
        for w in &d.warnings {
            match w.severity {
                DiagnosticSeverity::Info => {
                    info_count += 1;
                    worst_rank = worst_rank.max(1);
                }
                DiagnosticSeverity::Warning => {
                    warning_count += 1;
                    worst_rank = worst_rank.max(2);
                }
                DiagnosticSeverity::Error => {
                    error_count += 1;
                    worst_rank = worst_rank.max(3);
                }
            }
            code_set.insert(w.code.clone());
            let line = format!("{}: [{}] {}", r.physics, w.code, w.message);
            if warning_seen.insert(line.clone()) {
                warnings.push(line);
            }
        }
    }

    if results_with_diagnostics == 0 {
        return None;
    }

    let worst_severity = match worst_rank {
        3 => "error",
        2 => "warning",
        1 => "info",
        _ => "none",
    }
    .to_string();

    Some(RunDiagnosticsSummary {
        results_with_diagnostics,
        all_converged,
        info_count,
        warning_count,
        error_count,
        worst_severity,
        max_algebraic_residual,
        min_rcond_estimate,
        cholesky_fallbacks,
        warning_codes: code_set.into_iter().collect(),
        warnings,
    })
}

fn inline_entry(assign: &MaterialAssign) -> MaterialProvenance {
    let snap = serde_json::to_value(assign).map(|v| canonical_json(&v)).unwrap_or(Value::Null);
    MaterialProvenance {
        blake3: nyx_store::util::hash_bytes(snap.to_string().as_bytes()),
        source: "inline".to_string(),
        citations: Vec::new(),
        params_snapshot: snap,
    }
}

/// Resolves deterministic material provenance for all materials referenced by `project` against a pre-merged `library`.
pub fn resolve_materials_with_merged_library(
    project: &Project,
    library: &Library,
) -> BTreeMap<String, MaterialProvenance> {
    let mut out = BTreeMap::new();
    let mut named_refs: BTreeSet<String> = project.materials.keys().cloned().collect();

    for spec in &project.physics {
        for (region, assign) in &spec.materials {
            match assign {
                MaterialAssign::Named(name) => {
                    named_refs.insert(name.clone());
                }
                MaterialAssign::Oriented { material, .. } => {
                    named_refs.insert(material.clone());
                }
                MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_) => {
                    out.insert(format!("inline:{}:{region}", spec.physics), inline_entry(assign));
                }
            }
        }
    }

    for name in named_refs {
        let in_project = project.materials.keys().any(|k| k.eq_ignore_ascii_case(&name));
        if let Some(m) = library.get(&name) {
            let key = if m.name.is_empty() { name.clone() } else { m.name.clone() };
            let in_proj = in_project || project.materials.keys().any(|k| k.eq_ignore_ascii_case(&key));
            let blake3 = nyx_store::util::hash_bytes(m.to_toml().as_bytes());
            let citations = m.source.iter().map(|s| s.cite.clone()).collect();
            let params_snapshot = serde_json::to_value(m).map(|v| canonical_json(&v)).unwrap_or(Value::Null);
            out.insert(
                key,
                MaterialProvenance {
                    blake3,
                    source: if in_proj { "project" } else { "library" }.to_string(),
                    citations,
                    params_snapshot,
                },
            );
        } else {
            let blake3 = project
                .materials
                .iter()
                .find(|(k, _)| k.eq_ignore_ascii_case(&name))
                .map(|(_, txt)| nyx_store::util::hash_bytes(txt.as_bytes()))
                .unwrap_or_else(|| "missing".to_string());
            out.insert(
                name,
                MaterialProvenance {
                    blake3,
                    source: if in_project { "project" } else { "missing" }.to_string(),
                    citations: Vec::new(),
                    params_snapshot: Value::Null,
                },
            );
        }
    }

    out
}

/// Resolves deterministic material provenance for `project`, loading [`crate::materials::default_dirs`] if `base_library` is `None`.
pub fn resolve_materials_for_project(
    project: &Project,
    base_library: Option<&Library>,
) -> BTreeMap<String, MaterialProvenance> {
    let loaded;
    let base = match base_library {
        Some(l) => l,
        None => {
            loaded = crate::materials::load_library(&crate::materials::default_dirs()).0;
            &loaded
        }
    };
    if project.materials.is_empty() {
        return resolve_materials_with_merged_library(project, base);
    }
    match base.with_texts(&project.materials) {
        Ok(merged) => resolve_materials_with_merged_library(project, &merged),
        Err(_) => resolve_materials_with_merged_library(project, base),
    }
}

/// Compares stored `/provenance/materials_json` entries against the local merged material library and returns explicit drift warnings.
pub fn check_material_drift(
    stored: &BTreeMap<String, MaterialProvenance>,
    project: &Project,
    local_merged_library: &Library,
) -> Vec<String> {
    let mut warnings = Vec::new();
    let mut inline_now: BTreeMap<String, MaterialProvenance> = BTreeMap::new();
    for spec in &project.physics {
        for (region, assign) in &spec.materials {
            if matches!(assign, MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_)) {
                inline_now.insert(format!("inline:{}:{region}", spec.physics), inline_entry(assign));
            }
        }
    }

    for (name, prev) in stored {
        if prev.source == "inline" {
            match inline_now.get(name) {
                Some(cur) if cur.blake3 == prev.blake3 => {}
                Some(cur) => warnings.push(format!(
                    "material drift: inline specification '{name}' changed ({} -> {})",
                    prev.blake3, cur.blake3
                )),
                None => warnings.push(format!(
                    "material drift: inline specification '{name}' ({}) is missing from project",
                    prev.blake3
                )),
            }
            continue;
        }
        match local_merged_library.get(name) {
            Some(cur) => {
                let cur_hash = nyx_store::util::hash_bytes(cur.to_toml().as_bytes());
                if cur_hash != prev.blake3 {
                    warnings.push(format!(
                        "material drift: '{name}' differs from local material library (stored {}, local {cur_hash})",
                        prev.blake3
                    ));
                }
            }
            None => warnings.push(format!(
                "material drift: stored material '{name}' ({}) is missing from local material library",
                prev.blake3
            )),
        }
    }

    warnings
}
