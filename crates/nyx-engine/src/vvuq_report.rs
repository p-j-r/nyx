// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Unsigned, BLAKE3-linked canonical VVUQ provenance and error-budget report (`vvuq_report.v1`).
//!
//! Links cryptographic BLAKE3 digests of canonical project JSON, canonical solver mesh,
//! resolved active materials (TOML or inline canonical JSON), and optional saved/linked `results.h5`
//! with build provenance, `SimplexMesh::check()` integrity summary, per-physics
//! `SolverDiagnostics` & `NumericalErrorBudget`, `Validated` state/tier verdicts, and active
//! material provenance/clamping warnings.

use crate::engine::Engine;
use crate::project::Project;
use crate::validated::{Check, Validated, ValidationTier};
use femcore::mesh::SimplexMesh;
use nyx_physics::{
    audit_active_materials, canonical_json, Ctx, DiagnosticSeverity, DiagnosticWarning, Error, LabelResolver,
    MaterialAssign, NumericalErrorBudget, SolverDiagnostics,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// Contract format identifier for [`VvuqReportV1`].
pub const VVUQ_REPORT_FORMAT: &str = "vvuq_report.v1";

fn sanitize_opt_f64(v: Option<f64>) -> Option<f64> {
    v.filter(|x| x.is_finite())
}

fn opt_env_nonempty(v: Option<&'static str>) -> Option<String> {
    v.map(str::trim).filter(|s| !s.is_empty()).map(str::to_string)
}

fn parse_opt_bool(v: Option<&'static str>) -> Option<bool> {
    match v.map(str::trim) {
        Some("true") | Some("1") | Some("yes") => Some(true),
        Some("false") | Some("0") | Some("no") => Some(false),
        _ => None,
    }
}

fn fallback_target_triple() -> String {
    format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS)
}

fn is_material_warning_code(code: &str) -> bool {
    matches!(code, "material_approximate" | "material_unverified" | "temperature_clamped")
        || code.starts_with("material_")
}

struct NoopResolver;

impl LabelResolver for NoopResolver {
    fn domains(&self, _region: &str) -> Result<Vec<i32>, Error> {
        Ok(Vec::new())
    }

    fn boundaries(&self, _selector: &str) -> Result<Vec<i32>, Error> {
        Ok(Vec::new())
    }
}

/// Computes the canonical BLAKE3 digest (`b3:<32 hex>`) of a [`Project`]'s sorted canonical JSON.
pub fn canonical_project_blake3(project: &Project) -> String {
    nyx_store::util::hash_bytes(project.to_json().as_bytes())
}

/// Computes the canonical BLAKE3 digest (`b3:<32 hex>`) of a solver [`SimplexMesh`] across
/// dimension, coordinates, cell connectivity/labels, facet connectivity/labels, curved mid-edge
/// nodes, and Barsoum quarter-point edge indices.
pub fn canonical_mesh_blake3(mesh: &SimplexMesh) -> String {
    let mut buf = Vec::with_capacity(
        32 + mesh.points.len() * 8 + mesh.cells.len() * 4 + mesh.cell_labels.len() * 4 + mesh.facets.len() * 4,
    );
    buf.extend_from_slice(b"nyx-simplex-mesh-v1\0");
    buf.extend_from_slice(&(mesh.dim as u64).to_le_bytes());

    buf.extend_from_slice(&(mesh.points.len() as u64).to_le_bytes());
    for &x in &mesh.points {
        buf.extend_from_slice(&x.to_le_bytes());
    }

    buf.extend_from_slice(&(mesh.cells.len() as u64).to_le_bytes());
    for &c in &mesh.cells {
        buf.extend_from_slice(&c.to_le_bytes());
    }

    buf.extend_from_slice(&(mesh.cell_labels.len() as u64).to_le_bytes());
    for &l in &mesh.cell_labels {
        buf.extend_from_slice(&l.to_le_bytes());
    }

    buf.extend_from_slice(&(mesh.facets.len() as u64).to_le_bytes());
    for &f in &mesh.facets {
        buf.extend_from_slice(&f.to_le_bytes());
    }

    buf.extend_from_slice(&(mesh.facet_labels.len() as u64).to_le_bytes());
    for &l in &mesh.facet_labels {
        buf.extend_from_slice(&l.to_le_bytes());
    }

    let curved = mesh.edge_node_list();
    buf.extend_from_slice(&(curved.len() as u64).to_le_bytes());
    for ([a, b], pos) in curved {
        buf.extend_from_slice(&a.to_le_bytes());
        buf.extend_from_slice(&b.to_le_bytes());
        for c in pos {
            buf.extend_from_slice(&c.to_le_bytes());
        }
    }

    let qp = mesh.qp_edge_list();
    buf.extend_from_slice(&(qp.len() as u64).to_le_bytes());
    for [a, b] in qp {
        buf.extend_from_slice(&a.to_le_bytes());
        buf.extend_from_slice(&b.to_le_bytes());
    }

    nyx_store::util::hash_bytes(&buf)
}

/// Build and runtime provenance captured via `build.rs` env vars (with clean fallbacks).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BuildMetadata {
    pub nyx_version: String,
    pub git_commit: Option<String>,
    pub git_dirty: Option<bool>,
    pub rustc_version: Option<String>,
    pub target_triple: String,
    pub profile: String,
    pub cargo_lock_blake3: Option<String>,
    pub os: String,
    pub arch: String,
    pub nyx_threads: usize,
    pub nyx_factor_threads: usize,
}

impl BuildMetadata {
    /// Captures build metadata from compile-time `NYX_BUILD_*` environment variables (with safe
    /// fallbacks before `build.rs` is merged) plus active runtime thread configuration.
    pub fn current() -> Self {
        Self {
            nyx_version: env!("CARGO_PKG_VERSION").to_string(),
            git_commit: opt_env_nonempty(option_env!("NYX_BUILD_GIT_COMMIT")),
            git_dirty: parse_opt_bool(option_env!("NYX_BUILD_GIT_DIRTY")),
            rustc_version: opt_env_nonempty(option_env!("NYX_BUILD_RUSTC_VERSION")),
            target_triple: opt_env_nonempty(option_env!("NYX_BUILD_TARGET")).unwrap_or_else(fallback_target_triple),
            profile: opt_env_nonempty(option_env!("NYX_BUILD_PROFILE")).unwrap_or_else(|| {
                if cfg!(debug_assertions) {
                    "debug".to_string()
                } else {
                    "release".to_string()
                }
            }),
            cargo_lock_blake3: opt_env_nonempty(option_env!("NYX_BUILD_CARGO_LOCK_BLAKE3")),
            os: std::env::consts::OS.to_string(),
            arch: std::env::consts::ARCH.to_string(),
            nyx_threads: femcore::par::threads(),
            nyx_factor_threads: femcore::par::factor_threads(),
        }
    }
}

/// Summary of [`SimplexMesh::check()`] structural, orientation, manifold, and curved/QP integrity checks.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MeshCheckSummary {
    pub ok: bool,
    pub error: Option<String>,
    pub dim: usize,
    pub n_points: usize,
    pub n_cells: usize,
    pub n_facets: usize,
    pub n_edges: usize,
    pub n_curved_edges: usize,
    pub n_qp_edges: usize,
}

impl MeshCheckSummary {
    /// Evaluates [`SimplexMesh::check()`] and summarizes mesh counts and integrity status.
    pub fn from_mesh(mesh: &SimplexMesh) -> Self {
        let check = mesh.check();
        Self {
            ok: check.is_ok(),
            error: check.err().map(|e| e.to_string()),
            dim: mesh.dim,
            n_points: mesh.n_points(),
            n_cells: mesh.n_cells(),
            n_facets: mesh.n_facets(),
            n_edges: mesh.n_edges(),
            n_curved_edges: mesh.edge_node_list().len(),
            n_qp_edges: mesh.qp_edge_list().len(),
        }
    }
}

/// Active material provenance or temperature-clamping warning scoped to a physics problem.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct MaterialWarningEntry {
    pub physics: String,
    pub code: String,
    pub severity: DiagnosticSeverity,
    pub message: String,
}

/// Cell counts recorded by a mesh-convergence check.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValidatedCheckCells {
    pub coarse: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub medium: Option<u64>,
    pub refined: Option<u64>,
}

/// Verdict summary of a mesh-convergence check (`validated::Check`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValidatedCheckSummary {
    pub factor: f64,
    pub tolerance: f64,
    pub observable: String,
    pub matching: String,
    pub passed: bool,
    pub grids: usize,
    pub worst_relative: Option<f64>,
    pub worst_quantity: String,
    pub cells: ValidatedCheckCells,
    pub seconds: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r_eff: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub observed_order: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gci_21: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gci_regime: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_algebraic_residual: Option<f64>,
}

impl ValidatedCheckSummary {
    fn from_check(c: &Check) -> Self {
        Self {
            factor: c.factor,
            tolerance: c.tolerance,
            observable: c.observable(),
            matching: c.matching.clone(),
            passed: c.passed,
            grids: c.grids,
            worst_relative: sanitize_opt_f64(Some(c.worst_relative)),
            worst_quantity: c.worst_quantity.clone(),
            cells: ValidatedCheckCells { coarse: c.cells_coarse, medium: c.cells_medium, refined: c.cells_refined },
            seconds: if c.seconds.is_finite() { c.seconds } else { 0.0 },
            r_eff: sanitize_opt_f64(c.r_eff),
            observed_order: sanitize_opt_f64(c.observed_order),
            gci_21: sanitize_opt_f64(c.gci_21),
            gci_regime: c.gci_regime.clone(),
            max_algebraic_residual: sanitize_opt_f64(c.max_algebraic_residual),
        }
    }
}

/// Per-physics validation verdict (`state`, `tier`, `reason`, optional `check`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ValidatedReportEntry {
    pub state: String,
    pub tier: String,
    pub reason: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub check: Option<ValidatedCheckSummary>,
}

impl ValidatedReportEntry {
    fn from_validated(v: &Validated) -> Self {
        Self {
            state: v.state.as_str().to_string(),
            tier: v.tier.as_str().to_string(),
            reason: v.reason.clone(),
            check: v.check.as_ref().map(ValidatedCheckSummary::from_check),
        }
    }
}

/// Compact `(state, tier)` pair per physics inside top-level [`ValidatedSummary`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValidatedTierPair {
    pub state: String,
    pub tier: String,
}

/// Project-wide validation summary (weakest tier across all configured physics).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValidatedSummary {
    pub state: String,
    pub tier: String,
    pub reason: String,
    pub by_physics: BTreeMap<String, ValidatedTierPair>,
}

/// Per-physics solve, solver diagnostics, error budget, and validation status.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PhysicsVvuqEntry {
    pub study: String,
    pub order: u8,
    pub solved: bool,
    pub stale: bool,
    pub dofs: Option<usize>,
    pub seconds: Option<f64>,
    pub diagnostics: Option<SolverDiagnostics>,
    pub error_budget: Option<NumericalErrorBudget>,
    pub validated: ValidatedReportEntry,
}

/// Canonical unsigned `vvuq_report.v1.json` document linking digests, provenance, mesh check,
/// solver diagnostics, numerical error budgets, validation tiers, and material warnings.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct VvuqReportV1 {
    pub format: String,
    pub generated_utc: String,
    pub project_title: String,
    pub project_blake3: String,
    pub mesh_blake3: Option<String>,
    pub materials_blake3: BTreeMap<String, String>,
    pub results_h5_blake3: Option<String>,
    pub build: BuildMetadata,
    pub mesh_check: Option<MeshCheckSummary>,
    pub material_warnings: Vec<MaterialWarningEntry>,
    pub validated: ValidatedSummary,
    pub physics: BTreeMap<String, PhysicsVvuqEntry>,
}

impl VvuqReportV1 {
    /// Serializes the report as canonical sorted-key pretty JSON text.
    pub fn to_canonical_json(&self) -> Result<String, Error> {
        let val = serde_json::to_value(self).map_err(|e| Error::invalid(format!("vvuq_report serialize: {e}")))?;
        serde_json::to_string_pretty(&canonical_json(&val))
            .map_err(|e| Error::invalid(format!("vvuq_report print: {e}")))
    }
}

/// Returns the JSON Schema (`vvuq_report.v1.schema.json`) describing [`VvuqReportV1`].
pub fn vvuq_report_schema() -> Value {
    let opt_num = json!({ "type": ["number", "null"] });
    let opt_int = json!({ "type": ["integer", "null"], "minimum": 0 });
    let opt_str = json!({ "type": ["string", "null"] });
    let opt_bool = json!({ "type": ["boolean", "null"] });
    let state_enum = json!({ "enum": ["exploratory", "checked", "check_not_affordable"] });
    let tier_enum = json!({
        "enum": [
            "bound_certified",
            "gci_verified",
            "two_grid_checked",
            "residual_checked",
            "exploratory",
            "not_affordable",
            "non_convergent"
        ]
    });

    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://nyx.dev/contracts/vvuq_report.v1.schema.json",
        "title": VVUQ_REPORT_FORMAT,
        "description": "Unsigned canonical VVUQ provenance and numerical error-budget report emitted by `vvuq_export`.",
        "type": "object",
        "additionalProperties": false,
        "required": [
            "format",
            "generated_utc",
            "project_title",
            "project_blake3",
            "mesh_blake3",
            "materials_blake3",
            "results_h5_blake3",
            "build",
            "mesh_check",
            "material_warnings",
            "validated",
            "physics"
        ],
        "properties": {
            "format": { "const": VVUQ_REPORT_FORMAT },
            "generated_utc": { "type": "string" },
            "project_title": { "type": "string" },
            "project_blake3": { "type": "string" },
            "mesh_blake3": opt_str,
            "materials_blake3": {
                "type": "object",
                "additionalProperties": { "type": "string" }
            },
            "results_h5_blake3": opt_str,
            "build": { "$ref": "#/$defs/build_metadata" },
            "mesh_check": {
                "oneOf": [
                    { "type": "null" },
                    { "$ref": "#/$defs/mesh_check" }
                ]
            },
            "material_warnings": {
                "type": "array",
                "items": { "$ref": "#/$defs/material_warning" }
            },
            "validated": { "$ref": "#/$defs/validated_summary" },
            "physics": {
                "type": "object",
                "additionalProperties": { "$ref": "#/$defs/physics_entry" }
            }
        },
        "$defs": {
            "build_metadata": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "nyx_version",
                    "git_commit",
                    "git_dirty",
                    "rustc_version",
                    "target_triple",
                    "profile",
                    "cargo_lock_blake3",
                    "os",
                    "arch",
                    "nyx_threads",
                    "nyx_factor_threads"
                ],
                "properties": {
                    "nyx_version": { "type": "string" },
                    "git_commit": opt_str,
                    "git_dirty": opt_bool,
                    "rustc_version": opt_str,
                    "target_triple": { "type": "string" },
                    "profile": { "type": "string" },
                    "cargo_lock_blake3": opt_str,
                    "os": { "type": "string" },
                    "arch": { "type": "string" },
                    "nyx_threads": { "type": "integer", "minimum": 1 },
                    "nyx_factor_threads": { "type": "integer", "minimum": 1 }
                }
            },
            "mesh_check": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "ok",
                    "error",
                    "dim",
                    "n_points",
                    "n_cells",
                    "n_facets",
                    "n_edges",
                    "n_curved_edges",
                    "n_qp_edges"
                ],
                "properties": {
                    "ok": { "type": "boolean" },
                    "error": opt_str,
                    "dim": { "type": "integer", "enum": [2, 3] },
                    "n_points": { "type": "integer", "minimum": 0 },
                    "n_cells": { "type": "integer", "minimum": 0 },
                    "n_facets": { "type": "integer", "minimum": 0 },
                    "n_edges": { "type": "integer", "minimum": 0 },
                    "n_curved_edges": { "type": "integer", "minimum": 0 },
                    "n_qp_edges": { "type": "integer", "minimum": 0 }
                }
            },
            "material_warning": {
                "type": "object",
                "additionalProperties": false,
                "required": ["physics", "code", "severity", "message"],
                "properties": {
                    "physics": { "type": "string" },
                    "code": { "type": "string" },
                    "severity": { "enum": ["info", "warning", "error"] },
                    "message": { "type": "string" }
                }
            },
            "diagnostic_warning": {
                "type": "object",
                "additionalProperties": false,
                "required": ["code", "severity", "message"],
                "properties": {
                    "code": { "type": "string" },
                    "severity": { "enum": ["info", "warning", "error"] },
                    "message": { "type": "string" }
                }
            },
            "solver_diagnostics": {
                "type": "object",
                "additionalProperties": false,
                "required": ["solver_kind", "converged", "max_algebraic_residual", "cholesky_fallbacks"],
                "properties": {
                    "solver_kind": { "type": "string" },
                    "converged": { "type": "boolean" },
                    "iterations": { "type": "integer", "minimum": 0 },
                    "restarts": { "type": "integer", "minimum": 0 },
                    "max_algebraic_residual": { "type": "number" },
                    "per_mode_residuals": {
                        "type": "array",
                        "items": { "type": "number" }
                    },
                    "rcond_estimate": { "type": "number" },
                    "cholesky_fallbacks": { "type": "integer", "minimum": 0 },
                    "warnings": {
                        "type": "array",
                        "items": { "$ref": "#/$defs/diagnostic_warning" }
                    }
                }
            },
            "spatial_error_estimate": {
                "type": "object",
                "additionalProperties": false,
                "required": ["method", "rel_bound", "is_epistemic_bound"],
                "properties": {
                    "method": { "type": "string" },
                    "rel_bound": { "type": "number" },
                    "is_epistemic_bound": { "type": "boolean" },
                    "observed_order": { "type": "number" },
                    "refinement_ratio": { "type": "number" },
                    "safety_factor": { "type": "number" },
                    "extrapolated_scalars_only": { "type": "boolean" }
                }
            },
            "numerical_error_budget": {
                "type": "object",
                "additionalProperties": false,
                "required": ["algebraic_rel"],
                "properties": {
                    "algebraic_rel": { "type": "number" },
                    "spatial_discretization": { "$ref": "#/$defs/spatial_error_estimate" },
                    "temporal_rel": { "type": "number" },
                    "model_truncation_rel": { "type": "number" },
                    "combined_num_uncertainty_rel": { "type": "number" },
                    "alg_below_disc_ratio": { "type": "number" }
                }
            },
            "validated_check_cells": {
                "type": "object",
                "additionalProperties": false,
                "required": ["coarse", "refined"],
                "properties": {
                    "coarse": opt_int,
                    "medium": opt_int,
                    "refined": opt_int
                }
            },
            "validated_check": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "factor",
                    "tolerance",
                    "observable",
                    "matching",
                    "passed",
                    "grids",
                    "worst_relative",
                    "worst_quantity",
                    "cells",
                    "seconds"
                ],
                "properties": {
                    "factor": { "type": "number" },
                    "tolerance": { "type": "number" },
                    "observable": { "type": "string" },
                    "matching": { "type": "string" },
                    "passed": { "type": "boolean" },
                    "grids": { "type": "integer", "minimum": 2 },
                    "worst_relative": opt_num,
                    "worst_quantity": { "type": "string" },
                    "cells": { "$ref": "#/$defs/validated_check_cells" },
                    "seconds": { "type": "number" },
                    "r_eff": { "type": "number" },
                    "observed_order": { "type": "number" },
                    "gci_21": { "type": "number" },
                    "gci_regime": { "type": "string" },
                    "max_algebraic_residual": { "type": "number" }
                }
            },
            "validated_entry": {
                "type": "object",
                "additionalProperties": false,
                "required": ["state", "tier", "reason"],
                "properties": {
                    "state": state_enum,
                    "tier": tier_enum,
                    "reason": { "type": "string" },
                    "check": { "$ref": "#/$defs/validated_check" }
                }
            },
            "validated_tier_pair": {
                "type": "object",
                "additionalProperties": false,
                "required": ["state", "tier"],
                "properties": {
                    "state": state_enum,
                    "tier": tier_enum
                }
            },
            "validated_summary": {
                "type": "object",
                "additionalProperties": false,
                "required": ["state", "tier", "reason", "by_physics"],
                "properties": {
                    "state": state_enum,
                    "tier": tier_enum,
                    "reason": { "type": "string" },
                    "by_physics": {
                        "type": "object",
                        "additionalProperties": { "$ref": "#/$defs/validated_tier_pair" }
                    }
                }
            },
            "physics_entry": {
                "type": "object",
                "additionalProperties": false,
                "required": [
                    "study",
                    "order",
                    "solved",
                    "stale",
                    "dofs",
                    "seconds",
                    "diagnostics",
                    "error_budget",
                    "validated"
                ],
                "properties": {
                    "study": { "type": "string" },
                    "order": { "type": "integer", "minimum": 1 },
                    "solved": { "type": "boolean" },
                    "stale": { "type": "boolean" },
                    "dofs": opt_int,
                    "seconds": opt_num,
                    "diagnostics": {
                        "oneOf": [
                            { "type": "null" },
                            { "$ref": "#/$defs/solver_diagnostics" }
                        ]
                    },
                    "error_budget": {
                        "oneOf": [
                            { "type": "null" },
                            { "$ref": "#/$defs/numerical_error_budget" }
                        ]
                    },
                    "validated": { "$ref": "#/$defs/validated_entry" }
                }
            }
        }
    })
}

/// Deterministic canonical exemplar matching `docs/architecture/contracts/vvuq_report.v1.json`.
pub fn vvuq_report_exemplar() -> VvuqReportV1 {
    let mut materials_blake3 = BTreeMap::new();
    materials_blake3.insert("copper_4k".into(), "b3:11111111111111111111111111111111".into());
    let mut by_physics = BTreeMap::new();
    by_physics.insert("dfs".into(), ValidatedTierPair { state: "checked".into(), tier: "gci_verified".into() });
    let mut physics = BTreeMap::new();
    physics.insert(
        "dfs".into(),
        PhysicsVvuqEntry {
            study: "heat_stationary".into(),
            order: 2,
            solved: true,
            stale: false,
            dofs: Some(425),
            seconds: Some(0.015),
            diagnostics: Some(SolverDiagnostics {
                solver_kind: "direct_cholesky".into(),
                converged: true,
                iterations: None,
                restarts: None,
                max_algebraic_residual: 1.25e-15,
                per_mode_residuals: Vec::new(),
                rcond_estimate: Some(1.4e-4),
                cholesky_fallbacks: 0,
                warnings: vec![DiagnosticWarning {
                    code: "material_approximate".into(),
                    severity: DiagnosticSeverity::Warning,
                    message: "region 'bar' (material 'copper_4k' [thermal]): approximate RRR=100 copper fit".into(),
                }],
            }),
            error_budget: Some(NumericalErrorBudget {
                algebraic_rel: 1.25e-15,
                spatial_discretization: Some(nyx_physics::SpatialErrorEstimate {
                    method: "roache_gci_3grid".into(),
                    rel_bound: 2.4e-5,
                    is_epistemic_bound: false,
                    observed_order: Some(3.92),
                    refinement_ratio: Some(1.667),
                    safety_factor: Some(1.25),
                    extrapolated_scalars_only: None,
                }),
                temporal_rel: None,
                model_truncation_rel: None,
                combined_num_uncertainty_rel: Some(2.4e-5),
                alg_below_disc_ratio: Some(5.2e-11),
            }),
            validated: ValidatedReportEntry {
                state: "checked".into(),
                tier: "gci_verified".into(),
                reason: "3-grid Roache GCI check passed: element size times 0.6, GCI_21 = 2.40e-05 (p_obs = 3.92)"
                    .into(),
                check: Some(ValidatedCheckSummary {
                    factor: 0.6,
                    tolerance: 1e-3,
                    observable:
                        "1 quantity (eigenvalues, values, numeric scalars), solutions matched by index, 3-grid GCI"
                            .into(),
                    matching: "index".into(),
                    passed: true,
                    grids: 3,
                    worst_relative: Some(1.9e-5),
                    worst_quantity: "value[0]".into(),
                    cells: ValidatedCheckCells { coarse: Some(64), medium: Some(178), refined: Some(494) },
                    seconds: 0.08,
                    r_eff: Some(1.667),
                    observed_order: Some(3.92),
                    gci_21: Some(2.4e-5),
                    gci_regime: Some("asymptotic".into()),
                    max_algebraic_residual: Some(1.25e-15),
                }),
            },
        },
    );

    VvuqReportV1 {
        format: VVUQ_REPORT_FORMAT.into(),
        generated_utc: "2026-11-26T00:00:00Z".into(),
        project_title: "Cryogenic Thermal Bar".into(),
        project_blake3: "b3:00000000000000000000000000000000".into(),
        mesh_blake3: Some("b3:22222222222222222222222222222222".into()),
        materials_blake3,
        results_h5_blake3: Some("b3:33333333333333333333333333333333".into()),
        build: BuildMetadata {
            nyx_version: "0.1.0".into(),
            git_commit: Some("0123456789abcdef".into()),
            git_dirty: Some(false),
            rustc_version: Some("rustc 1.99.1".into()),
            target_triple: "x86_64-unknown-linux-gnu".into(),
            profile: "release".into(),
            cargo_lock_blake3: Some("b3:44444444444444444444444444444444".into()),
            os: "linux".into(),
            arch: "x86_64".into(),
            nyx_threads: 4,
            nyx_factor_threads: 4,
        },
        mesh_check: Some(MeshCheckSummary {
            ok: true,
            error: None,
            dim: 2,
            n_points: 45,
            n_cells: 64,
            n_facets: 32,
            n_edges: 108,
            n_curved_edges: 0,
            n_qp_edges: 0,
        }),
        material_warnings: vec![MaterialWarningEntry {
            physics: "dfs".into(),
            code: "material_approximate".into(),
            severity: DiagnosticSeverity::Warning,
            message: "region 'bar' (material 'copper_4k' [thermal]): approximate RRR=100 copper fit".into(),
        }],
        validated: ValidatedSummary {
            state: "checked".into(),
            tier: "gci_verified".into(),
            reason: "weakest validated tier across configured physics (dfs)".into(),
            by_physics,
        },
        physics,
    }
}

impl Engine {
    fn collect_materials_blake3(&self) -> Result<BTreeMap<String, String>, Error> {
        let lib = self.library()?;
        let mut out = BTreeMap::new();
        for spec in &self.project.physics {
            for (region, assign) in &spec.materials {
                match assign {
                    MaterialAssign::Named(name) | MaterialAssign::Oriented { material: name, .. } => {
                        if !out.contains_key(name) {
                            if let Some(mat) = lib.get(name) {
                                out.insert(name.clone(), nyx_store::util::hash_bytes(mat.to_toml().as_bytes()));
                            } else if let Some(raw) = self.project.materials.get(name) {
                                out.insert(name.clone(), nyx_store::util::hash_bytes(raw.as_bytes()));
                            }
                        }
                    }
                    MaterialAssign::Props(_) | MaterialAssign::PropsTensor(_) => {
                        let key = format!("inline:{}:{region}", spec.physics);
                        let v = serde_json::to_value(assign)
                            .map_err(|e| Error::invalid(format!("material spec '{key}': {e}")))?;
                        out.insert(key, nyx_store::util::hash_bytes(canonical_json(&v).to_string().as_bytes()));
                    }
                }
            }
            for bc in &spec.bcs {
                if let Some(name) =
                    bc.params.get("material").and_then(Value::as_str).map(str::trim).filter(|s| !s.is_empty())
                {
                    if !out.contains_key(name) {
                        if let Some(mat) = lib.get(name) {
                            out.insert(name.to_string(), nyx_store::util::hash_bytes(mat.to_toml().as_bytes()));
                        } else if let Some(raw) = self.project.materials.get(name) {
                            out.insert(name.to_string(), nyx_store::util::hash_bytes(raw.as_bytes()));
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    fn collect_material_warnings(&self) -> Result<Vec<MaterialWarningEntry>, Error> {
        let lib = self.library()?;
        let mut seen: BTreeSet<(String, String, String)> = BTreeSet::new();
        let mut out = Vec::new();
        let mut push = |physics: &str, w: &DiagnosticWarning| {
            if is_material_warning_code(&w.code)
                && seen.insert((physics.to_string(), w.code.clone(), w.message.clone()))
            {
                out.push(MaterialWarningEntry {
                    physics: physics.to_string(),
                    code: w.code.clone(),
                    severity: w.severity,
                    message: w.message.clone(),
                });
            }
        };

        let noop = NoopResolver;
        for spec in &self.project.physics {
            let stored_fresh = self.results.get(&spec.physics).filter(|_| !self.is_stale(&spec.physics));
            if let Some(stored) = stored_fresh {
                if let Some(diag) = &stored.result.diagnostics {
                    for w in &diag.warnings {
                        push(&spec.physics, w);
                    }
                    continue;
                }
            }
            if let Some(p) = self.registry.get(&spec.physics) {
                let mut ctx = Ctx::new(&noop);
                ctx.library = Some(lib.clone());
                ctx.wavelength = crate::materials::study_wavelength(p, spec);
                ctx.frequency = crate::materials::study_frequency(p, spec);
                let props = p.properties();
                let temps = crate::material_warnings::collect_pre_solve_temperatures(spec, &ctx);
                if let Ok(audit) = audit_active_materials(spec, &ctx, &props, &temps) {
                    for w in audit.all_warnings() {
                        push(&spec.physics, &w);
                    }
                }
            }
        }
        Ok(out)
    }

    /// Returns the linked `results.h5` path when still valid for the current project and results,
    /// and allowed by the active [`path_guard`](crate::path_guard) root.
    pub(crate) fn linked_results_h5(&self) -> Option<PathBuf> {
        let (stamp, path) = self.last_results_h5.as_ref()?;
        if *stamp != self.stamp() || !path.is_file() {
            return None;
        }
        if self.project.physics.iter().any(|s| self.results.contains_key(&s.physics) && self.is_stale(&s.physics)) {
            return None;
        }
        if let Some(root) = &self.path_root {
            let s = path.to_str()?;
            crate::path_guard::confine(root, s).ok()?;
        }
        Some(path.clone())
    }

    /// Builds an unsigned [`VvuqReportV1`] for the current project, mesh, materials, results, and checks.
    pub fn vvuq_report(&self, results_h5: Option<&Path>) -> Result<VvuqReportV1, Error> {
        let project_blake3 = canonical_project_blake3(&self.project);
        let (mesh_blake3, mesh_check) = if let Some(model) = self.model() {
            let mesh = model.solver_mesh();
            (Some(canonical_mesh_blake3(&mesh)), Some(MeshCheckSummary::from_mesh(&mesh)))
        } else {
            (None, None)
        };

        let materials_blake3 = self.collect_materials_blake3()?;
        let h5_candidate = results_h5.map(Path::to_path_buf).or_else(|| self.linked_results_h5());
        let results_h5_blake3 = match h5_candidate {
            Some(p) if p.is_file() => {
                if let Some(root) = &self.path_root {
                    let s = p.to_str().ok_or_else(|| Error::invalid("results.h5 path is not valid UTF-8"))?;
                    crate::path_guard::confine(root, s)
                        .map_err(|why| Error::invalid(format!("vvuq_export: `results_h5`: {why}")))?;
                }
                Some(nyx_store::util::hash_file(&p).map_err(|e| Error::invalid(format!("{}: {e}", p.display())))?)
            }
            _ => None,
        };

        let material_warnings = self.collect_material_warnings()?;

        let mut physics = BTreeMap::new();
        let mut by_physics = BTreeMap::new();
        let mut weakest: Option<(String, Validated)> = None;
        for spec in &self.project.physics {
            let id = &spec.physics;
            let stored = self.results.get(id);
            let stale = stored.is_none_or(|_| self.is_stale(id));
            let val = self.validated(id);
            by_physics.insert(
                id.clone(),
                ValidatedTierPair { state: val.state.as_str().to_string(), tier: val.tier.as_str().to_string() },
            );
            let take = weakest.as_ref().is_none_or(|(_, cur)| val.tier.weakness() > cur.tier.weakness());
            if take {
                weakest = Some((id.clone(), val.clone()));
            }
            physics.insert(
                id.clone(),
                PhysicsVvuqEntry {
                    study: spec.study.kind.clone(),
                    order: spec.order,
                    solved: stored.is_some(),
                    stale,
                    dofs: stored.map(|s| s.result.dofs),
                    seconds: stored.map(|s| s.result.seconds),
                    diagnostics: stored.and_then(|s| s.result.diagnostics.clone()),
                    error_budget: stored.and_then(|s| s.result.error_budget.clone()),
                    validated: ValidatedReportEntry::from_validated(&val),
                },
            );
        }

        let validated = match weakest {
            Some((id, w)) => ValidatedSummary {
                state: w.state.as_str().to_string(),
                tier: w.tier.as_str().to_string(),
                reason: if physics.len() == 1 {
                    w.reason
                } else {
                    format!("weakest validated tier across configured physics ({id}): {}", w.reason)
                },
                by_physics,
            },
            None => ValidatedSummary {
                state: "exploratory".into(),
                tier: ValidationTier::Exploratory.as_str().to_string(),
                reason: "the project has no physics configured".into(),
                by_physics,
            },
        };

        Ok(VvuqReportV1 {
            format: VVUQ_REPORT_FORMAT.into(),
            generated_utc: nyx_store::util::iso(nyx_store::util::now_secs()),
            project_title: self.project.title.clone(),
            project_blake3,
            mesh_blake3,
            materials_blake3,
            results_h5_blake3,
            build: BuildMetadata::current(),
            mesh_check,
            material_warnings,
            validated,
            physics,
        })
    }

    /// Exports an unsigned canonical `vvuq_report.v1.json` file to `path` (confined to [`Engine::path_root`]
    /// when set), optionally exporting a companion `results.h5` file beside the report first.
    pub fn vvuq_export(&mut self, path: &str, export_h5: bool) -> Result<Value, Error> {
        let out_path = crate::path_guard::check_path(self.path_root.as_deref(), "vvuq_export", "path", path)?;
        if let Some(parent) = out_path.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent).map_err(|e| Error::invalid(format!("{}: {e}", parent.display())))?;
        }
        let mut h5_path = self.linked_results_h5();
        if h5_path.is_none() && export_h5 {
            let companion_raw = out_path.with_file_name("results.h5");
            let companion_str = companion_raw
                .to_str()
                .ok_or_else(|| Error::invalid("vvuq_export: companion results.h5 path is not valid UTF-8"))?;
            let companion =
                crate::path_guard::check_path(self.path_root.as_deref(), "vvuq_export", "path", companion_str)?;
            let snapshot = self.results_snapshot()?;
            crate::results_file::write_snapshot(&snapshot, &companion.to_string_lossy(), false, true)?;
            self.last_results_h5 = Some((self.stamp(), companion.clone()));
            h5_path = Some(companion);
        }

        let report = self.vvuq_report(h5_path.as_deref())?;
        let text = report.to_canonical_json()?;
        let bytes = format!("{text}\n");
        nyx_store::util::write_atomic(&out_path, bytes.as_bytes())
            .map_err(|e| Error::invalid(format!("{}: {e}", out_path.display())))?;

        let report_val = serde_json::to_value(&report).map_err(|e| Error::invalid(format!("vvuq_report json: {e}")))?;
        let mut reply = json!({
            "path": out_path.to_string_lossy(),
            "format": VVUQ_REPORT_FORMAT,
            "bytes": bytes.len() as u64,
            "report": report_val,
        });
        if let Some(hp) = h5_path {
            reply["results_h5_path"] = json!(hp.to_string_lossy());
        }
        Ok(reply)
    }
}
