// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Artifact schemas (`*.layout.json`, `sampling_manifest.json`, and `benchmark_runs.jsonl`) and failure-taxonomy
//! classification for the SQuADDS cQED benchmark runner (`Plan §4.5` & `§5.1`).

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Dielectric and sample-holder vertical layer stack (`Plan §2.2.A` & `§4.5`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StackSpec {
    #[serde(default = "default_h_sub", alias = "t_sub_m", alias = "substrate_thickness_m")]
    pub h_sub_m: f64,
    #[serde(default = "default_h_air_top", alias = "h_vac_top_m", alias = "sample_holder_top_m")]
    pub h_air_top_m: f64,
    #[serde(default = "default_h_air_bot", alias = "h_vac_bot_m", alias = "sample_holder_bottom_m")]
    pub h_air_bot_m: f64,
    #[serde(default)]
    pub z_min_m: Option<f64>,
    #[serde(default)]
    pub z_max_m: Option<f64>,
    #[serde(default = "default_eps_sub", alias = "eps_r_sub")]
    pub eps_sub: f64,
    #[serde(default = "default_t_metal", alias = "metal_thickness_m")]
    pub t_metal_m: f64,
    #[serde(default)]
    pub trench_depth_m: f64,
}

fn default_h_sub() -> f64 {
    750.0e-6
}
fn default_h_air_top() -> f64 {
    890.0e-6
}
fn default_h_air_bot() -> f64 {
    900.0e-6
}
fn default_eps_sub() -> f64 {
    11.45
}
fn default_t_metal() -> f64 {
    200.0e-9
}

impl Default for StackSpec {
    fn default() -> Self {
        Self {
            h_sub_m: default_h_sub(),
            h_air_top_m: default_h_air_top(),
            h_air_bot_m: default_h_air_bot(),
            z_min_m: None,
            z_max_m: None,
            eps_sub: default_eps_sub(),
            t_metal_m: default_t_metal(),
            trench_depth_m: 0.0,
        }
    }
}

impl StackSpec {
    /// Lower z-coordinate of the bottom vacuum cavity below `z = -h_sub_m`.
    /// Handles both conventions cleanly: if `z_min_m` is explicitly supplied, uses it; otherwise if
    /// `h_air_bot_m` was exported as `sample_holder_bottom = 1650 um` (total depth from `z = 0`) vs
    /// `900 um` (thickness below `-h_sub_m`), checks `z_min_m` first and treats `h_air_bot_m` as the
    /// vacuum thickness below the wafer backplane unless `h_air_bot_m` equals `1650e-6` with `750e-6` Si.
    pub fn vac_bot_zmin(&self) -> f64 {
        if let Some(z_min) = self.z_min_m {
            return z_min;
        }
        if (self.h_sub_m - 750.0e-6).abs() < 1.0e-9 && (self.h_air_bot_m - 1650.0e-6).abs() < 1.0e-9 {
            -1650.0e-6
        } else {
            -(self.h_sub_m + self.h_air_bot_m.max(0.0))
        }
    }

    /// Thickness of the lower vacuum layer `>= 0`.
    pub fn vac_bot_thickness(&self) -> f64 {
        (-self.h_sub_m - self.vac_bot_zmin()).max(0.0)
    }

    /// Upper z-coordinate of the top air/vacuum cavity above `z = 0`.
    pub fn vac_top_zmax(&self) -> f64 {
        self.z_max_m.unwrap_or(self.h_air_top_m)
    }
}

/// Planar or elevated 2D conductive/port/junction polygon sheet (`Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SheetSpec {
    pub name: String,
    #[serde(default)]
    pub role: String,
    #[serde(default, alias = "net")]
    pub terminal: Option<String>,
    #[serde(default, alias = "z")]
    pub z_m: f64,
    #[serde(default = "default_simplify_tol", alias = "simplify_tol")]
    pub simplify_tol_m: f64,
    #[serde(default)]
    pub arc_tol_m: f64,
    #[serde(default)]
    pub edge_h_m: Option<f64>,
    #[serde(default)]
    pub face_h_m: Option<f64>,
    #[serde(default)]
    pub thickness_m: Option<f64>,
    #[serde(default)]
    pub etch_depth_m: Option<f64>,
    #[serde(default)]
    pub polygons: Vec<Vec<[f64; 2]>>,
    #[serde(default)]
    pub holes: Vec<Vec<[f64; 2]>>,
}

fn default_simplify_tol() -> f64 {
    50.0e-9
}

/// 3D rectangular hole landing post for CPW airbridge crossovers (`Plan §3.3.D` & `§4.5`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HolePostSpec {
    pub name: String,
    #[serde(alias = "min")]
    pub min_m: [f64; 3],
    #[serde(alias = "max")]
    pub max_m: [f64; 3],
}

/// Lumped port or Josephson junction specification attached to a sheet (`Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PortSpec {
    #[serde(default)]
    pub name: Option<String>,
    pub sheet: String,
    /// `"lumped_port"` | `"lumped_jj"`
    pub kind: String,
    #[serde(default)]
    pub z0_ohm: Option<f64>,
    #[serde(default = "default_port_dir")]
    pub direction: [f64; 3],
    #[serde(default)]
    pub l_j_h: Option<f64>,
    #[serde(default)]
    pub e_j_hz: Option<f64>,
    #[serde(default)]
    pub c_j_f: Option<f64>,
    #[serde(default)]
    pub e_l_hz: Option<f64>,
    #[serde(default, alias = "squid_asymmetry")]
    pub squid_d: Option<f64>,
    #[serde(default)]
    pub mode: Option<String>,
}

fn default_port_dir() -> [f64; 3] {
    [0.0, 1.0, 0.0]
}

impl PortSpec {
    pub fn port_name(&self) -> String {
        self.name.clone().unwrap_or_else(|| self.sheet.clone())
    }
}

/// Physics-agnostic layout snapshot (`*.layout.json`, `Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct LayoutSpec {
    #[serde(default = "default_schema_version")]
    pub schema_version: u32,
    pub case_id: String,
    pub suite_id: String,
    #[serde(default)]
    pub tier: Option<String>,
    #[serde(default, alias = "stratum", alias = "family", alias = "topology")]
    pub cluster_id: Option<String>,
    #[serde(default)]
    pub stratum_weight: Option<f64>,
    #[serde(default)]
    pub n_stratum_pop: Option<usize>,
    #[serde(default)]
    pub n_stratum_sample: Option<usize>,
    #[serde(default = "default_die_bbox", alias = "bbox_m")]
    pub die_bbox_m: [f64; 4],
    #[serde(default)]
    pub stack: StackSpec,
    #[serde(default)]
    pub sheets: Vec<SheetSpec>,
    #[serde(default)]
    pub hole_posts: Vec<HolePostSpec>,
    #[serde(default)]
    pub ports: Vec<PortSpec>,
    #[serde(default)]
    pub wire_bonds_logged: Value,
    #[serde(default)]
    pub references_si: BTreeMap<String, f64>,
    #[serde(default)]
    pub metal_thickness_m: Option<f64>,
    #[serde(default)]
    pub etch_bias_m: Option<f64>,
    #[serde(default)]
    pub squid_asymmetry: Option<f64>,
    #[serde(default)]
    pub e_j_hz: Option<f64>,
    #[serde(default)]
    pub c_j_f: Option<f64>,
    #[serde(default)]
    pub l_j_h: Option<f64>,
    #[serde(default)]
    pub f_search_hz: Option<f64>,
    #[serde(default)]
    pub n_modes: Option<usize>,
    #[serde(default)]
    pub h_global_m: Option<f64>,
    #[serde(default)]
    pub min_gap_m: Option<f64>,
    #[serde(default)]
    pub resonator_type: Option<String>,
    #[serde(default)]
    pub wm1_level: Option<String>,
    #[serde(default)]
    pub mode: Option<String>,
    #[serde(default)]
    pub gds_path: Option<String>,
    #[serde(default)]
    pub gds_layers: Option<BTreeMap<String, [i16; 2]>>,
    #[serde(default)]
    pub gds_pin_layer: Option<[i16; 2]>,
    #[serde(default)]
    pub gds_cell: Option<String>,
    #[serde(default)]
    pub coupler_layout_path: Option<String>,
    #[serde(default)]
    pub qubit_layout_path: Option<String>,
    #[serde(default)]
    pub cavity_layout_path: Option<String>,
    #[serde(flatten)]
    pub extra: Map<String, Value>,
}

fn default_schema_version() -> u32 {
    1
}

fn default_die_bbox() -> [f64; 4] {
    [-500.0e-6, -500.0e-6, 500.0e-6, 500.0e-6]
}

impl LayoutSpec {
    /// Effective physical metal film thickness in metres (`t_metal_m`).
    pub fn effective_metal_thickness_m(&self) -> f64 {
        self.metal_thickness_m.unwrap_or(self.stack.t_metal_m)
    }

    /// Characteristic lateral chip dimension `L_chip = max(W_chip, H_chip)`.
    pub fn chip_length_m(&self) -> f64 {
        let w = (self.die_bbox_m[2] - self.die_bbox_m[0]).abs();
        let h = (self.die_bbox_m[3] - self.die_bbox_m[1]).abs();
        w.max(h).max(100.0e-6)
    }

    /// Horvitz–Thompson stratum weight (`N_stratum / n_stratum` when population/sample counts are given, else `stratum_weight` or `1.0`).
    pub fn ht_weight(&self) -> f64 {
        if let Some(w) = self.stratum_weight.filter(|&w| w.is_finite() && w > 0.0) {
            return w;
        }
        if let (Some(n_pop), Some(n_samp)) = (self.n_stratum_pop, self.n_stratum_sample) {
            if n_pop > 0 && n_samp > 0 {
                return (n_pop as f64) / (n_samp as f64);
            }
        }
        1.0
    }

    /// Cluster identifier for cluster bootstrap resampling.
    pub fn effective_cluster_id(&self) -> String {
        if let Some(ref c) = self.cluster_id {
            if !c.is_empty() {
                return c.clone();
            }
        }
        format!("{}:{}", self.suite_id, self.case_id)
    }
}

/// Stage wall-clock breakdown in seconds (`Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct StageTimings {
    pub arrange: f64,
    pub mesh: f64,
    pub assemble: f64,
    pub solve: f64,
    pub post: f64,
    pub total: f64,
}

/// Process peak memory counters in bytes (`Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct MemoryBytes {
    pub peak_ws: u64,
    pub peak_commit: u64,
}

/// Scoped algebraic and physical invariants (`Plan §4.5` & `§5.1.3`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct InvariantsRecord {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sym_err: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_eig: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_offdiag_pos_ratio: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub equipartition_err: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_noncollinearity_max: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_prom_res: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reciprocity_err: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unitarity_err: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_overlap: Option<f64>,
}

/// Canonical per-run JSONL record (`benchmark_runs.jsonl`, `Plan §4.5`).
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RunRecord {
    pub case_id: String,
    pub suite_id: String,
    pub tier: String,
    #[serde(default)]
    pub mode: String,
    #[serde(default)]
    pub cluster_id: String,
    #[serde(default = "default_weight")]
    pub stratum_weight: f64,
    /// `"ok" | "mesh_fail" | "oom_refused" | "solver_nonconv" | "prom_tol_exceeded" | "invariant_fail" | "timeout"`
    pub status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub failure_detail: Option<String>,
    pub cells: u64,
    pub dofs: u64,
    pub order: u8,
    pub threads: usize,
    pub timings_s: StageTimings,
    pub memory_bytes: MemoryBytes,
    pub invariants: InvariantsRecord,
    pub predictions_si: BTreeMap<String, f64>,
    pub references_si: BTreeMap<String, f64>,
    pub rel_errors_pct: BTreeMap<String, f64>,
    pub abs_errors_unit: BTreeMap<String, f64>,
    #[serde(default, skip_serializing_if = "Value::is_null")]
    pub extra: Value,
}

fn default_weight() -> f64 {
    1.0
}

/// Unified small-coupling relative-error denominator floors (`Plan §5.1.2`).
pub fn denom_floor_si(metric_key: &str) -> f64 {
    let k = metric_key.to_ascii_lowercase();
    if k.starts_with("c_")
        || k.ends_with("_f")
        || k.contains("to_cross")
        || k.contains("to_claw")
        || k.contains("to_ground")
        || k.contains("to_top")
        || k.contains("to_bottom")
        || k.contains("to_south")
        || k.contains("to_north")
        || k.contains("cap")
    {
        0.5e-15 // |C_ref| >= 0.5 fF
    } else if k.contains("chi") {
        25.0e3 // |chi_ref| / 2pi >= 25 kHz
    } else if k.contains("kappa") {
        20.0e3 // kappa_ref / 2pi >= 20 kHz
    } else {
        1.0e-30
    }
}

/// Physical unit scale (`fF`, `MHz`, `kHz`, `us`) for `abs_errors_unit` (`Plan §5.1.1`).
pub fn unit_scale_and_suffix(metric_key: &str) -> (f64, &'static str) {
    let k = metric_key.to_ascii_lowercase();
    if k.starts_with("c_")
        || k.ends_with("_f")
        || k.contains("to_cross")
        || k.contains("to_claw")
        || k.contains("to_ground")
        || k.contains("to_top")
        || k.contains("to_bottom")
        || k.contains("to_south")
        || k.contains("to_north")
        || k.contains("cap")
    {
        (1.0e15, "fF")
    } else if k.contains("chi") || (k.contains("kappa") && !k.contains("kappa_f")) {
        (1.0e-3, "kHz")
    } else if k.contains("t1") {
        (1.0e6, "us")
    } else if k.contains("hz")
        || k.starts_with("f_")
        || k.starts_with("f01")
        || k.starts_with("fr")
        || k.contains("alpha")
        || k.starts_with("g_")
        || k == "g"
        || k.contains("kappa_f")
    {
        (1.0e-6, "MHz")
    } else {
        (1.0, "SI")
    }
}

/// Populates `rel_errors_pct` (respecting unified denominator floors) and `abs_errors_unit` (retaining 100% of entries),
/// plus normalized Frobenius capacitance error `c_frobenius_rel_pct` when multiple capacitance entries exist.
pub fn populate_errors(rec: &mut RunRecord) {
    rec.rel_errors_pct.clear();
    rec.abs_errors_unit.clear();

    let mut frob_diff_sq = 0.0_f64;
    let mut frob_ref_sq = 0.0_f64;
    let mut n_cap_terms = 0usize;

    for (key, &y_ref) in &rec.references_si {
        if !y_ref.is_finite() {
            continue;
        }
        let Some(&y_pred) = rec.predictions_si.get(key) else {
            continue;
        };
        if !y_pred.is_finite() {
            continue;
        }
        let diff = y_pred - y_ref;
        let (scale, suffix) = unit_scale_and_suffix(key);
        rec.abs_errors_unit.insert(format!("{key}_{suffix}"), diff.abs() * scale);

        if suffix == "fF" && !key.contains("raw") {
            frob_diff_sq += diff * diff;
            frob_ref_sq += y_ref * y_ref;
            n_cap_terms += 1;
        }

        let floor = denom_floor_si(key);
        if y_ref.abs() >= floor && y_ref.abs() > 0.0 {
            let rel_pct = 100.0 * diff.abs() / y_ref.abs();
            rec.rel_errors_pct.insert(key.clone(), rel_pct);
        }
    }

    // Also allow explicit raw capacitance keys `raw_<key>` to compare against `<key>` reference
    let raw_keys: Vec<(String, f64)> = rec
        .predictions_si
        .iter()
        .filter_map(|(k, &v)| k.strip_prefix("raw_").map(|base| (base.to_string(), v)))
        .collect();
    for (base, y_raw) in raw_keys {
        if let Some(&y_ref) = rec.references_si.get(&base) {
            if y_ref.is_finite() && y_raw.is_finite() {
                let (scale, suffix) = unit_scale_and_suffix(&base);
                rec.abs_errors_unit.insert(format!("raw_{base}_{suffix}"), (y_raw - y_ref).abs() * scale);
                let floor = denom_floor_si(&base);
                if y_ref.abs() >= floor && y_ref.abs() > 0.0 {
                    rec.rel_errors_pct.insert(format!("raw_{base}"), 100.0 * (y_raw - y_ref).abs() / y_ref.abs());
                }
            }
        }
    }

    if n_cap_terms >= 2 && frob_ref_sq > 0.0 {
        let frob_pct = 100.0 * (frob_diff_sq.sqrt() / frob_ref_sq.sqrt());
        rec.rel_errors_pct.insert("c_frobenius_rel".into(), frob_pct);
    }
}

/// Classifies an engine/mesher/solver error string onto the canonical `Plan §4.5` taxonomy.
pub fn classify_failure(msg: &str) -> (&'static str, String) {
    let m = msg.to_ascii_lowercase();
    let status = if m.contains("memory_wait_timeout")
        || (m.contains("memory")
            && (m.contains("need") || m.contains("refuse") || m.contains("budget") || m.contains("oom")))
        || m.contains("max_cells")
    {
        "oom_refused"
    } else if m.contains("timeout") || m.contains("timed out") {
        "timeout"
    } else if m.contains("prom") && (m.contains("residual") || m.contains("tol")) {
        "prom_tol_exceeded"
    } else if m.contains("invariant") {
        "invariant_fail"
    } else if m.contains("mesh")
        || m.contains("ruppert")
        || m.contains("steiner")
        || m.contains("pslg")
        || m.contains("arrange")
        || m.contains("polygon")
        || m.contains("overlay")
        || m.contains("sliver")
        || m.contains("jacobian")
    {
        "mesh_fail"
    } else {
        "solver_nonconv"
    };
    (status, msg.to_string())
}

/// Sanitizes Python-emitted bare `NaN`, `Infinity`, and `-Infinity` tokens outside JSON string literals into `null`
/// so `serde_json::from_str` parses `Option<f64>` fields cleanly without failing on non-RFC-8259 float tokens.
pub fn parse_json_lenient_nan<T: serde::de::DeserializeOwned>(txt: &str) -> Result<T, serde_json::Error> {
    if !txt.contains("NaN") && !txt.contains("Infinity") {
        return serde_json::from_str(txt);
    }
    let bytes = txt.as_bytes();
    let mut out = String::with_capacity(txt.len());
    let mut i = 0usize;
    let mut in_str = false;
    let mut esc = false;
    let is_ident_char = |b: u8| b.is_ascii_alphanumeric() || b == b'_';

    while i < bytes.len() {
        let b = bytes[i];
        if in_str {
            out.push(b as char);
            if esc {
                esc = false;
            } else if b == b'\\' {
                esc = true;
            } else if b == b'"' {
                in_str = false;
            }
            i += 1;
            continue;
        }
        if b == b'"' {
            in_str = true;
            out.push('"');
            i += 1;
            continue;
        }
        let prev_ident = i > 0 && is_ident_char(bytes[i - 1]);
        if !prev_ident {
            if bytes[i..].starts_with(b"NaN") && (i + 3 == bytes.len() || !is_ident_char(bytes[i + 3])) {
                out.push_str("null");
                i += 3;
                continue;
            }
            if bytes[i..].starts_with(b"-Infinity") && (i + 9 == bytes.len() || !is_ident_char(bytes[i + 9])) {
                out.push_str("null");
                i += 9;
                continue;
            }
            if bytes[i..].starts_with(b"Infinity") && (i + 8 == bytes.len() || !is_ident_char(bytes[i + 8])) {
                out.push_str("null");
                i += 8;
                continue;
            }
        }
        out.push(b as char);
        i += 1;
    }
    serde_json::from_str(&out)
}

/// One resolved manifest item ready to execute.
#[derive(Clone, Debug)]
pub struct ManifestItem {
    pub case_id: String,
    pub suite_id: String,
    pub tier: String,
    pub layout_path: PathBuf,
    pub mode_override: Option<String>,
    pub cluster_id: Option<String>,
    pub stratum_weight: Option<f64>,
}

fn tier_matches(entry_tier: &str, requested_tier: &str) -> bool {
    let e = entry_tier.trim().to_ascii_lowercase();
    let r = requested_tier.trim().to_ascii_lowercase();
    if r == "all" || e == r {
        return true;
    }
    match r.as_str() {
        "p1" | "phase1" | "phase_1" => matches!(e.as_str(), "p1" | "phase1" | "phase_1" | "core_p1" | "track_a_p1"),
        "p1_diag" => e == "p1_diag",
        "p2" | "phase2" | "phase_2" => matches!(e.as_str(), "p2" | "phase2" | "phase_2" | "core_p2" | "track_a_p2"),
        "core" => matches!(
            e.as_str(),
            "core" | "p1" | "p2" | "phase1" | "phase2" | "phase_1" | "phase_2" | "core_p1" | "core_p2"
        ),
        "qcap" | "ref0" | "ref.0" => matches!(e.as_str(), "qcap" | "ref0" | "ref.0" | "p3_qcap" | "phase3_qcap"),
        "airbridge" | "airbridges" | "s3.4" => {
            matches!(e.as_str(), "airbridge" | "airbridges" | "s3.4" | "p3_airbridge")
        }
        "ablation" => matches!(e.as_str(), "ablation" | "p2" | "phase2" | "core"),
        _ => false,
    }
}

fn suite_matches(entry_suite: &str, filter_suite: Option<&str>) -> bool {
    let Some(f) = filter_suite else {
        return true;
    };
    if f.trim().is_empty() || f.eq_ignore_ascii_case("all") {
        return true;
    }
    let es = entry_suite.to_ascii_lowercase();
    f.split(',').any(|tok| {
        let t = tok.trim().to_ascii_lowercase();
        !t.is_empty() && (es == t || es.starts_with(&t) || t.starts_with(&es))
    })
}

fn resolve_layout_file(inputs_dir: &Path, manifest_dir: &Path, rel_or_abs: &str) -> Option<PathBuf> {
    let p = PathBuf::from(rel_or_abs);
    if p.is_absolute() && p.is_file() {
        return Some(p);
    }
    for base in [inputs_dir, manifest_dir, inputs_dir.parent().unwrap_or(inputs_dir)] {
        let cand = base.join(&p);
        if cand.is_file() {
            return Some(cand);
        }
    }
    None
}

fn collect_layout_files_recursive(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|x| x.path())).collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            collect_layout_files_recursive(&path, out);
        } else if path.file_name().and_then(|s| s.to_str()).is_some_and(|n| n.ends_with(".layout.json")) {
            out.push(path);
        }
    }
}

fn parse_manifest_object(
    obj: &Map<String, Value>,
    inherited_tier: Option<&str>,
    inherited_suite: Option<&str>,
    inputs_dir: &Path,
    manifest_dir: &Path,
    file_index: &BTreeMap<String, PathBuf>,
) -> Option<ManifestItem> {
    let rel_path = obj
        .get("layout_path")
        .or_else(|| obj.get("layout"))
        .or_else(|| obj.get("path"))
        .or_else(|| obj.get("file"))
        .and_then(Value::as_str);

    let case_id_opt = obj.get("case_id").or_else(|| obj.get("id")).and_then(Value::as_str);
    if rel_path.is_none() && case_id_opt.is_none() {
        return None;
    }

    let resolved = rel_path.and_then(|rp| resolve_layout_file(inputs_dir, manifest_dir, rp)).or_else(|| {
        case_id_opt
            .and_then(|cid| file_index.get(&format!("{cid}.layout.json")).or_else(|| file_index.get(cid)).cloned())
    })?;

    let mut case_id = case_id_opt.unwrap_or("").to_string();
    let mut suite_id = obj
        .get("suite_id")
        .or_else(|| obj.get("dataset_id"))
        .or_else(|| obj.get("suite"))
        .and_then(Value::as_str)
        .or(inherited_suite)
        .unwrap_or("")
        .to_string();
    let mut tier = obj
        .get("tier")
        .or_else(|| obj.get("phase"))
        .and_then(Value::as_str)
        .or(inherited_tier)
        .unwrap_or("smoke")
        .to_string();

    if case_id.is_empty() || suite_id.is_empty() {
        if let Ok(txt) = std::fs::read_to_string(&resolved) {
            if let Ok(lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                if case_id.is_empty() {
                    case_id = lay.case_id;
                }
                if suite_id.is_empty() {
                    suite_id = lay.suite_id;
                }
                if let Some(t) = lay.tier {
                    if inherited_tier.is_none() && obj.get("tier").is_none() {
                        tier = t;
                    }
                }
            }
        }
    }

    let mode_override = obj.get("mode").or_else(|| obj.get("wm1_level")).and_then(Value::as_str).map(String::from);
    let cluster_id = obj
        .get("cluster_id")
        .or_else(|| obj.get("stratum"))
        .or_else(|| obj.get("family"))
        .and_then(Value::as_str)
        .map(String::from);
    let stratum_weight = obj.get("stratum_weight").or_else(|| obj.get("ht_weight")).and_then(Value::as_f64);

    Some(ManifestItem { case_id, suite_id, tier, layout_path: resolved, mode_override, cluster_id, stratum_weight })
}

fn walk_manifest_value(
    val: &Value,
    inherited_tier: Option<&str>,
    inherited_suite: Option<&str>,
    inputs_dir: &Path,
    manifest_dir: &Path,
    file_index: &BTreeMap<String, PathBuf>,
    out: &mut Vec<ManifestItem>,
) {
    match val {
        Value::Array(arr) => {
            for item in arr {
                walk_manifest_value(item, inherited_tier, inherited_suite, inputs_dir, manifest_dir, file_index, out);
            }
        }
        Value::Object(obj) => {
            if let Some(item) =
                parse_manifest_object(obj, inherited_tier, inherited_suite, inputs_dir, manifest_dir, file_index)
            {
                out.push(item);
                return;
            }
            let next_suite =
                obj.get("suite_id").or_else(|| obj.get("dataset_id")).and_then(Value::as_str).or(inherited_suite);
            let next_tier = obj.get("tier").and_then(Value::as_str).or(inherited_tier);

            for (k, child) in obj {
                let k_low = k.to_ascii_lowercase();
                if k_low == "ablation_recommended_12_ids" || k_low == "provenance" || k_low == "tier_counts" {
                    continue;
                }
                let child_tier = if matches!(
                    k_low.as_str(),
                    "smoke"
                        | "p1"
                        | "p1_diag"
                        | "p2"
                        | "phase1"
                        | "phase2"
                        | "phase_1"
                        | "phase_2"
                        | "core"
                        | "qcap"
                        | "airbridge"
                        | "airbridges"
                        | "ablation"
                ) {
                    Some(k.as_str())
                } else {
                    next_tier
                };
                let child_suite = if k_low.starts_with("s1.")
                    || k_low.starts_with("s2.")
                    || k_low.starts_with("s3.")
                    || k_low.starts_with("s4.")
                    || k_low.starts_with("ref")
                {
                    Some(k.as_str())
                } else {
                    next_suite
                };
                walk_manifest_value(child, child_tier, child_suite, inputs_dir, manifest_dir, file_index, out);
            }
        }
        Value::String(s) if s.ends_with(".layout.json") => {
            if let Some(resolved) = resolve_layout_file(inputs_dir, manifest_dir, s) {
                if let Ok(txt) = std::fs::read_to_string(&resolved) {
                    if let Ok(lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                        out.push(ManifestItem {
                            case_id: lay.case_id,
                            suite_id: lay.suite_id,
                            tier: inherited_tier.map(String::from).or(lay.tier).unwrap_or_else(|| "smoke".into()),
                            layout_path: resolved,
                            mode_override: lay.mode.or(lay.wm1_level),
                            cluster_id: lay.cluster_id,
                            stratum_weight: lay.stratum_weight,
                        });
                    }
                }
            }
        }
        _ => {}
    }
}

/// Loads and filters manifest items for `tier` and optional `suite` filter.
pub fn load_manifest_items(
    manifest_path: &Path,
    inputs_dir: &Path,
    tier: &str,
    suite_filter: Option<&str>,
) -> Result<Vec<ManifestItem>, String> {
    let mut all_files = Vec::new();
    collect_layout_files_recursive(inputs_dir, &mut all_files);
    let mut file_index = BTreeMap::new();
    for f in &all_files {
        if let Some(name) = f.file_name().and_then(|s| s.to_str()) {
            file_index.insert(name.to_string(), f.clone());
        }
    }

    let mut raw_items = Vec::new();
    let mut ablation_rec_ids: Vec<String> = Vec::new();
    if manifest_path.is_file() {
        let txt = std::fs::read_to_string(manifest_path)
            .map_err(|e| format!("failed to read manifest '{}': {e}", manifest_path.display()))?;
        let val: Value = parse_json_lenient_nan(&txt)
            .map_err(|e| format!("failed to parse JSON manifest '{}': {e}", manifest_path.display()))?;
        if let Some(arr) = val.get("ablation_recommended_12_ids").and_then(Value::as_array) {
            ablation_rec_ids = arr.iter().filter_map(Value::as_str).map(String::from).collect();
        }
        let manifest_dir = manifest_path.parent().unwrap_or(inputs_dir);
        walk_manifest_value(&val, None, None, inputs_dir, manifest_dir, &file_index, &mut raw_items);
    } else {
        for f in &all_files {
            if let Ok(txt) = std::fs::read_to_string(f) {
                if let Ok(lay) = parse_json_lenient_nan::<LayoutSpec>(&txt) {
                    raw_items.push(ManifestItem {
                        case_id: lay.case_id,
                        suite_id: lay.suite_id,
                        tier: lay.tier.unwrap_or_else(|| "smoke".into()),
                        layout_path: f.clone(),
                        mode_override: lay.mode.or(lay.wm1_level),
                        cluster_id: lay.cluster_id,
                        stratum_weight: lay.stratum_weight,
                    });
                }
            }
        }
    }

    let has_explicit_qcap = raw_items.iter().any(|it| tier_matches(&it.tier, "qcap"));
    let has_explicit_airbridge = raw_items.iter().any(|it| tier_matches(&it.tier, "airbridge"));

    // Deduplicate by (case_id, tier, mode_override, layout_path)
    let mut seen = std::collections::BTreeSet::new();
    let mut filtered: Vec<ManifestItem> = raw_items
        .into_iter()
        .filter(|it| {
            let s_low = it.suite_id.to_ascii_lowercase();
            let tier_ok = if tier.eq_ignore_ascii_case("qcap") {
                if has_explicit_qcap {
                    tier_matches(&it.tier, tier)
                } else {
                    s_low.contains("qcap") || s_low.starts_with("ref")
                }
            } else if tier.eq_ignore_ascii_case("airbridge") || tier.eq_ignore_ascii_case("airbridges") {
                if has_explicit_airbridge {
                    tier_matches(&it.tier, tier)
                } else {
                    s_low.starts_with("s3.4") || s_low.contains("airbridge")
                }
            } else if tier.eq_ignore_ascii_case("p1") {
                if tier_matches(&it.tier, tier) {
                    true
                } else if it.tier.eq_ignore_ascii_case("core") {
                    s_low.starts_with("s4.") || s_low.starts_with("s2.3") || s_low.starts_with("s3.2")
                } else {
                    false
                }
            } else if tier.eq_ignore_ascii_case("p2") {
                if tier_matches(&it.tier, tier) {
                    true
                } else if it.tier.eq_ignore_ascii_case("core") {
                    s_low.starts_with("s1.")
                        || s_low.starts_with("s2.1")
                        || s_low.starts_with("s2.2")
                        || s_low.starts_with("s3.1")
                } else {
                    false
                }
            } else {
                tier_matches(&it.tier, tier)
            };
            tier_ok && suite_matches(&it.suite_id, suite_filter)
        })
        .filter(|it| {
            let key = (
                it.case_id.clone(),
                it.tier.to_ascii_lowercase(),
                it.mode_override.clone().unwrap_or_default(),
                it.layout_path.clone(),
            );
            seen.insert(key)
        })
        .collect();

    // If `--tier ablation` was requested and no explicit `"ablation"` tier existed in the manifest,
    // use `ablation_recommended_12_ids` when present, or deterministically pick 3 stratified designs per S1.1..S1.4.
    if tier.eq_ignore_ascii_case("ablation") && filtered.iter().all(|it| !it.tier.eq_ignore_ascii_case("ablation")) {
        if !ablation_rec_ids.is_empty() {
            let mut by_id: BTreeMap<String, ManifestItem> = BTreeMap::new();
            for mut it in filtered {
                it.tier = "ablation".into();
                it.mode_override = None;
                by_id.insert(it.case_id.clone(), it);
            }
            let mut ab_sel = Vec::new();
            for cid in &ablation_rec_ids {
                if let Some(it) = by_id.remove(cid) {
                    ab_sel.push(it);
                }
            }
            filtered = ab_sel;
        } else {
            let mut s1_only: Vec<ManifestItem> =
                filtered.into_iter().filter(|it| it.suite_id.to_ascii_lowercase().starts_with("s1.")).collect();
            s1_only.sort_by(|a, b| a.suite_id.cmp(&b.suite_id).then_with(|| a.case_id.cmp(&b.case_id)));
            let mut by_suite: BTreeMap<String, Vec<ManifestItem>> = BTreeMap::new();
            for mut it in s1_only {
                it.tier = "ablation".into();
                let key = it.suite_id.split(':').next().unwrap_or(&it.suite_id).trim().to_uppercase();
                by_suite.entry(key).or_default().push(it);
            }
            let mut ab_sel = Vec::new();
            for (_s, group) in by_suite {
                if group.len() <= 3 {
                    ab_sel.extend(group);
                } else {
                    let n = group.len();
                    for idx in [0, n / 2, n - 1] {
                        ab_sel.push(group[idx].clone());
                    }
                }
            }
            filtered = ab_sel;
        }
    }

    Ok(filtered)
}
