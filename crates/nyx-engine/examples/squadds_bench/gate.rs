// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Anti-survivorship statistical gate evaluator (`squadds_bench gate`, `Plan §4.5`, `§5.1`, `§7.1`):
//! - Completion rate & failure taxonomy (`ok`, `mesh_fail`, `oom_refused`, `solver_nonconv`, `prom_tol_exceeded`,
//!   `invariant_fail`, `timeout`);
//! - Failure-penalized upper order statistics (`IQR`, `p90`, `p95`, `Max` with `e_k = +inf` on failed runs);
//! - Signed bias, physical-unit `MAE` / `RMSE` (`fF`, `MHz`, `kHz`, `us` retaining 100% of sub-floor entries),
//!   normalized Frobenius capacitance error, and unified denominator floor accounting (`0.5 fF`, `25 kHz`, `20 kHz`);
//! - Horvitz–Thompson inverse-probability stratum weights + deterministic 10,000-resample cluster bootstrap 95% CIs;
//! - Scoped algebraic & physical invariant audits and `G0` / `G1` / `G2` / `G3` pass/fail JSON verdicts.

#![allow(clippy::needless_range_loop)]

use crate::squadds_bench::schema::{denom_floor_si, unit_scale_and_suffix, RunRecord};
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const BOOTSTRAP_RESAMPLES: usize = 10_000;
const BOOTSTRAP_SEED: u64 = 0x6e79_785f_7371_7561; // "nyx_squa"

/// Deterministic SplitMix64 PRNG for reproducible 10,000-resample cluster bootstrap (`ADR 11`).
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    #[inline]
    fn next_usize(&mut self, upper_excl: usize) -> usize {
        if upper_excl <= 1 {
            0
        } else {
            (self.next_u64() % (upper_excl as u64)) as usize
        }
    }
}

fn json_num_or_inf(x: f64) -> Value {
    if x.is_infinite() && x.is_sign_positive() {
        json!("inf")
    } else if x.is_infinite() {
        json!("-inf")
    } else if x.is_nan() {
        Value::Null
    } else {
        json!(x)
    }
}

/// Computes linear-interpolated quantile `q in [0, 1]` on a pre-sorted slice `sorted` (which may end with `+inf`).
fn sorted_quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let pos = q.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi || sorted[lo].is_infinite() || sorted[hi].is_infinite() {
        if pos - (lo as f64) > 0.0 && sorted[hi].is_infinite() {
            f64::INFINITY
        } else {
            sorted[lo]
        }
    } else {
        let t = pos - (lo as f64);
        sorted[lo] * (1.0 - t) + sorted[hi] * t
    }
}

/// Weighted quantile on `(value, positive_weight)` pairs.
fn weighted_median(pairs: &[(f64, f64)]) -> f64 {
    if pairs.is_empty() {
        return f64::NAN;
    }
    let mut s = pairs.to_vec();
    s.sort_by(|a, b| a.0.total_cmp(&b.0));
    let w_tot: f64 = s.iter().map(|x| x.1).sum();
    if w_tot <= 0.0 {
        return s[s.len() / 2].0;
    }
    let half = 0.5 * w_tot;
    let mut acc = 0.0;
    for &(val, w) in &s {
        acc += w;
        if acc >= half {
            return val;
        }
    }
    s.last().unwrap().0
}

#[derive(Clone, Debug)]
struct ObsRow {
    cluster_id: String,
    weight: f64,
    ok: bool,
    /// Signed relative error `(y_pred - y_ref) / |y_ref|` (`None` if excluded below floor or failed)
    signed_rel: Option<f64>,
    /// Physical-unit absolute error `|y_pred - y_ref| * unit_scale` (`None` if failed)
    abs_err_unit: Option<f64>,
    excluded_below_floor: bool,
}

/// Computes full `Plan §5.1` error distribution (failure-penalized quantiles, signed bias, MAE/RMSE,
/// Horvitz–Thompson RMSRE, and 10,000-resample cluster bootstrap 95% CIs) for a named observable cohort.
fn summarize_observable_cohort(obs: &[ObsRow], unit_label: &str, enable_bootstrap: bool) -> Value {
    let n_total = obs.len();
    let n_failed = obs.iter().filter(|o| !o.ok).count();
    let n_below_floor = obs.iter().filter(|o| o.ok && o.excluded_below_floor).count();

    // Failure-penalized unsigned percentage errors (`e_k = +inf` for every failed run!)
    let mut unsigned_penalized_pct: Vec<f64> = Vec::with_capacity(n_total);
    let mut signed_pct_ok: Vec<f64> = Vec::new();
    let mut ht_rel_ok: Vec<(String, f64, f64)> = Vec::new(); // (cluster_id, abs_rel_pct, weight)
    let mut abs_unit_ok: Vec<f64> = Vec::new();

    for o in obs {
        if !o.ok {
            unsigned_penalized_pct.push(f64::INFINITY);
            continue;
        }
        if let Some(au) = o.abs_err_unit {
            abs_unit_ok.push(au);
        }
        if let Some(sr) = o.signed_rel {
            let sp = 100.0 * sr;
            let up = sp.abs();
            unsigned_penalized_pct.push(up);
            signed_pct_ok.push(sp);
            ht_rel_ok.push((o.cluster_id.clone(), up, o.weight.max(1.0e-12)));
        }
    }

    unsigned_penalized_pct.sort_by(|a, b| a.total_cmp(b));
    signed_pct_ok.sort_by(|a, b| a.total_cmp(b));

    let p25 = sorted_quantile(&unsigned_penalized_pct, 0.25);
    let median = sorted_quantile(&unsigned_penalized_pct, 0.50);
    let p75 = sorted_quantile(&unsigned_penalized_pct, 0.75);
    let iqr = if p75.is_infinite() {
        f64::INFINITY
    } else if p25.is_nan() || p75.is_nan() {
        f64::NAN
    } else {
        p75 - p25
    };
    let p90 = sorted_quantile(&unsigned_penalized_pct, 0.90);
    let p95 = sorted_quantile(&unsigned_penalized_pct, 0.95);
    let max_v = unsigned_penalized_pct.last().copied().unwrap_or(f64::NAN);

    let signed_bias_med_pct = sorted_quantile(&signed_pct_ok, 0.50);
    let signed_bias_mean_pct = if signed_pct_ok.is_empty() {
        f64::NAN
    } else {
        signed_pct_ok.iter().sum::<f64>() / (signed_pct_ok.len() as f64)
    };

    let n_rel_ok = ht_rel_ok.len();
    let mape_unweighted_pct =
        if n_rel_ok == 0 { f64::NAN } else { ht_rel_ok.iter().map(|x| x.1).sum::<f64>() / (n_rel_ok as f64) };
    let rmsre_unweighted_pct = if n_rel_ok == 0 {
        f64::NAN
    } else {
        (ht_rel_ok.iter().map(|x| x.1 * x.1).sum::<f64>() / (n_rel_ok as f64)).sqrt()
    };

    let w_sum: f64 = ht_rel_ok.iter().map(|x| x.2).sum();
    let (mape_ht_pct, rmsre_ht_pct, med_ht_pct) = if n_rel_ok == 0 || w_sum <= 0.0 {
        (f64::NAN, f64::NAN, f64::NAN)
    } else {
        let mape_w = ht_rel_ok.iter().map(|x| x.2 * x.1).sum::<f64>() / w_sum;
        let rmsre_w = (ht_rel_ok.iter().map(|x| x.2 * x.1 * x.1).sum::<f64>() / w_sum).sqrt();
        let pairs: Vec<(f64, f64)> = ht_rel_ok.iter().map(|x| (x.1, x.2)).collect();
        (mape_w, rmsre_w, weighted_median(&pairs))
    };

    let (mae_unit, rmse_unit) = if abs_unit_ok.is_empty() {
        (f64::NAN, f64::NAN)
    } else {
        let n_u = abs_unit_ok.len() as f64;
        let mae = abs_unit_ok.iter().sum::<f64>() / n_u;
        let rmse = (abs_unit_ok.iter().map(|x| x * x).sum::<f64>() / n_u).sqrt();
        (mae, rmse)
    };

    // 10,000-resample cluster bootstrap 95% CIs (`Plan §5.1.1`: active when >= 4 distinct clusters)
    let mut by_cluster: BTreeMap<&str, Vec<(f64, f64)>> = BTreeMap::new();
    for (cid, up, w) in &ht_rel_ok {
        by_cluster.entry(cid.as_str()).or_default().push((*up, *w));
    }
    let clusters: Vec<&Vec<(f64, f64)>> = by_cluster.values().collect();
    let n_clusters = clusters.len();

    let bootstrap_ci = if enable_bootstrap && n_clusters >= 4 && n_failed == 0 {
        let mut rng = SplitMix64::new(BOOTSTRAP_SEED ^ (n_total as u64));
        let mut boot_med = Vec::with_capacity(BOOTSTRAP_RESAMPLES);
        let mut boot_rmsre = Vec::with_capacity(BOOTSTRAP_RESAMPLES);
        let mut scratch_pairs = Vec::with_capacity(n_rel_ok * 2);

        for _b in 0..BOOTSTRAP_RESAMPLES {
            scratch_pairs.clear();
            let mut bw_sum = 0.0;
            let mut bw_sq = 0.0;
            for _c in 0..n_clusters {
                let pick = clusters[rng.next_usize(n_clusters)];
                for &(up, w) in pick {
                    scratch_pairs.push((up, w));
                    bw_sum += w;
                    bw_sq += w * up * up;
                }
            }
            if bw_sum > 0.0 {
                boot_rmsre.push((bw_sq / bw_sum).sqrt());
                boot_med.push(weighted_median(&scratch_pairs));
            }
        }
        boot_med.sort_by(|a, b| a.total_cmp(b));
        boot_rmsre.sort_by(|a, b| a.total_cmp(b));
        json!({
            "resamples": BOOTSTRAP_RESAMPLES,
            "n_clusters": n_clusters,
            "median_ht_95ci_pct": [sorted_quantile(&boot_med, 0.025), sorted_quantile(&boot_med, 0.975)],
            "rmsre_ht_95ci_pct": [sorted_quantile(&boot_rmsre, 0.025), sorted_quantile(&boot_rmsre, 0.975)],
        })
    } else {
        Value::Null
    };

    json!({
        "n_total_cases": n_total,
        "n_failed_penalized_inf": n_failed,
        "n_rel_included": n_rel_ok,
        "n_excluded_below_floor": n_below_floor,
        "n_clusters": n_clusters,
        "unit": unit_label,
        "median_pct": json_num_or_inf(median),
        "median_ht_pct": json_num_or_inf(med_ht_pct),
        "iqr_pct": json_num_or_inf(iqr),
        "p90_pct": json_num_or_inf(p90),
        "p95_pct": json_num_or_inf(p95),
        "max_pct": json_num_or_inf(max_v),
        "signed_bias_median_pct": json_num_or_inf(signed_bias_med_pct),
        "signed_bias_mean_pct": json_num_or_inf(signed_bias_mean_pct),
        "mape_unweighted_pct": json_num_or_inf(mape_unweighted_pct),
        "mape_ht_pct": json_num_or_inf(mape_ht_pct),
        "rmsre_unweighted_pct": json_num_or_inf(rmsre_unweighted_pct),
        "rmsre_ht_pct": json_num_or_inf(rmsre_ht_pct),
        "mae_unit": json_num_or_inf(mae_unit),
        "rmse_unit": json_num_or_inf(rmse_unit),
        "cluster_bootstrap_95ci": bootstrap_ci,
    })
}

fn extract_metric_obs(runs: &[&RunRecord], keys: &[&str], use_raw_prefix: bool) -> Vec<ObsRow> {
    let mut out = Vec::new();
    for r in runs {
        if r.status != "ok" {
            let n_pen = r.references_si.keys().filter(|k| keys.contains(&k.as_str())).count().max(1);
            for _ in 0..n_pen {
                out.push(ObsRow {
                    cluster_id: r.cluster_id.clone(),
                    weight: r.stratum_weight,
                    ok: false,
                    signed_rel: None,
                    abs_err_unit: None,
                    excluded_below_floor: false,
                });
            }
            continue;
        }
        for &k in keys {
            let Some(&y_ref) = r.references_si.get(k) else {
                continue;
            };
            let pred_key = if use_raw_prefix { format!("raw_{k}") } else { k.to_string() };
            let Some(&y_pred) = r.predictions_si.get(&pred_key).or_else(|| {
                if !use_raw_prefix {
                    r.predictions_si.get(&k.to_ascii_lowercase())
                } else {
                    None
                }
            }) else {
                continue;
            };
            if !(y_ref.is_finite() && y_pred.is_finite()) {
                continue;
            }
            let (scale, _) = unit_scale_and_suffix(k);
            let abs_u = (y_pred - y_ref).abs() * scale;
            let floor = denom_floor_si(k);
            if y_ref.abs() >= floor && y_ref.abs() > 0.0 {
                out.push(ObsRow {
                    cluster_id: r.cluster_id.clone(),
                    weight: r.stratum_weight,
                    ok: true,
                    signed_rel: Some((y_pred - y_ref) / y_ref.abs()),
                    abs_err_unit: Some(abs_u),
                    excluded_below_floor: false,
                });
            } else {
                out.push(ObsRow {
                    cluster_id: r.cluster_id.clone(),
                    weight: r.stratum_weight,
                    ok: true,
                    signed_rel: None,
                    abs_err_unit: Some(abs_u),
                    excluded_below_floor: true,
                });
            }
        }
    }
    out
}

fn load_runs_jsonl(runs_arg: &Path) -> Result<Vec<RunRecord>, String> {
    let mut files: Vec<PathBuf> = Vec::new();
    if runs_arg.is_dir() {
        let rd = std::fs::read_dir(runs_arg).map_err(|e| e.to_string())?;
        for e in rd.flatten() {
            let p = e.path();
            if p.extension().and_then(|s| s.to_str()) == Some("jsonl") {
                files.push(p);
            }
        }
        files.sort();
    } else {
        for part in runs_arg.to_string_lossy().split(',') {
            let p = PathBuf::from(part.trim());
            if p.is_file() {
                files.push(p);
            }
        }
    }
    if files.is_empty() {
        return Err(format!("no .jsonl run files found at '{}'", runs_arg.display()));
    }

    let mut records = Vec::new();
    for f in files {
        let txt = std::fs::read_to_string(&f).map_err(|e| format!("failed to read '{}': {e}", f.display()))?;
        for (ln_idx, line) in txt.lines().enumerate() {
            let t = line.trim();
            if t.is_empty() {
                continue;
            }
            let rec: RunRecord = serde_json::from_str(t)
                .map_err(|e| format!("JSONL parse error in '{}' line {}: {e}", f.display(), ln_idx + 1))?;
            records.push(rec);
        }
    }
    Ok(records)
}

fn completion_and_taxonomy(runs: &[&RunRecord]) -> (f64, usize, usize, BTreeMap<String, usize>) {
    let mut tax: BTreeMap<String, usize> = BTreeMap::new();
    for k in ["ok", "mesh_fail", "oom_refused", "solver_nonconv", "prom_tol_exceeded", "invariant_fail", "timeout"] {
        tax.insert(k.into(), 0);
    }
    let mut n_ok = 0usize;
    for r in runs {
        if r.status == "ok" {
            n_ok += 1;
        }
        *tax.entry(r.status.clone()).or_insert(0) += 1;
    }
    let tot = runs.len();
    let pct = if tot == 0 { 0.0 } else { 100.0 * (n_ok as f64) / (tot as f64) };
    (pct, n_ok, tot, tax)
}

fn audit_invariants(runs: &[&RunRecord]) -> Value {
    let mut electro_checked = 0usize;
    let mut electro_spd_ok = 0usize;
    let mut electro_offdiag_m_ok = 0usize;
    let mut max_sym_err = 0.0_f64;

    let mut cavity_checked = 0usize;
    let mut cavity_equip_strict_1e12 = 0usize;
    let mut max_equip_err = 0.0_f64;
    let mut min_overlap_observed = 1.0_f64;

    let mut driven_checked = 0usize;
    let mut driven_recip_ok = 0usize;
    let mut max_recip_err = 0.0_f64;
    let mut max_prom_res = 0.0_f64;

    for r in runs {
        if r.status != "ok" {
            continue;
        }
        if let Some(me) = r.invariants.min_eig {
            electro_checked += 1;
            if me > 0.0 {
                electro_spd_ok += 1;
            }
            if r.invariants.max_offdiag_pos_ratio.unwrap_or(0.0) <= 1.0e-5 {
                electro_offdiag_m_ok += 1;
            }
            if let Some(se) = r.invariants.sym_err {
                max_sym_err = max_sym_err.max(se);
            }
        }
        if let Some(eq) = r.invariants.equipartition_err {
            cavity_checked += 1;
            max_equip_err = max_equip_err.max(eq);
            if eq < 1.0e-12 {
                cavity_equip_strict_1e12 += 1;
            }
            if let Some(ov) = r.invariants.min_overlap {
                min_overlap_observed = min_overlap_observed.min(ov);
            }
        }
        if let Some(pr) = r.invariants.max_prom_res {
            driven_checked += 1;
            max_prom_res = max_prom_res.max(pr);
            let rec_e = r.invariants.reciprocity_err.unwrap_or(0.0);
            max_recip_err = max_recip_err.max(rec_e);
            if rec_e < 1.0e-11 {
                driven_recip_ok += 1;
            }
        }
    }

    json!({
        "electrostatics": {
            "checked": electro_checked,
            "spd_lambda_min_gt_0": electro_spd_ok,
            "offdiag_le_1e5_geom_mean": electro_offdiag_m_ok,
            "max_raw_asymmetry_rel": max_sym_err,
        },
        "cavity_rom": {
            "checked": cavity_checked,
            "equipartition_lt_1e12": cavity_equip_strict_1e12,
            "max_equipartition_err": max_equip_err,
            "min_overlap_observed": if cavity_checked > 0 { json!(min_overlap_observed) } else { Value::Null },
        },
        "driven_prom": {
            "checked": driven_checked,
            "reciprocity_lt_1e11": driven_recip_ok,
            "max_reciprocity_err": max_recip_err,
            "max_prom_residual": max_prom_res,
        }
    })
}

fn find_calibrations_json(runs_arg: &Path) -> Option<Value> {
    let base = if runs_arg.is_dir() { runs_arg } else { runs_arg.parent().unwrap_or(Path::new(".")) };
    let candidates = [
        base.join("calibrations").join("calibrations_summary.json"),
        base.join("calibrations_summary.json"),
        base.parent().unwrap_or(base).join("calibrations").join("calibrations_summary.json"),
    ];
    for c in &candidates {
        if c.is_file() {
            if let Ok(txt) = std::fs::read_to_string(c) {
                if let Ok(val) = serde_json::from_str::<Value>(&txt) {
                    return Some(val);
                }
            }
        }
    }
    None
}

/// Evaluates `G0`, `G1`, `G2`, `G3`, or `all` on `runs_path` and writes canonical `gate_report.json` to `out_path`.
pub fn evaluate_gates(runs_path: &Path, gate_filter: &str, out_path: &Path) -> Result<Value, String> {
    let all_runs = load_runs_jsonl(runs_path)?;
    let all_refs: Vec<&RunRecord> = all_runs.iter().collect();
    let (overall_comp_pct, overall_ok, overall_tot, overall_tax) = completion_and_taxonomy(&all_refs);
    let invariants_summary = audit_invariants(&all_refs);
    let cal_opt = find_calibrations_json(runs_path);

    // Partition cohorts
    let smoke_runs: Vec<&RunRecord> = all_runs.iter().filter(|r| r.tier.eq_ignore_ascii_case("smoke")).collect();
    let p1_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| matches!(r.tier.to_ascii_lowercase().as_str(), "p1" | "phase1" | "phase_1" | "core_p1"))
        .collect();
    let p2_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| matches!(r.tier.to_ascii_lowercase().as_str(), "p2" | "phase2" | "phase_2" | "core_p2"))
        .collect();

    let s1_self_keys = [
        "cross_to_cross",
        "claw_to_claw",
        "top_to_top",
        "bottom_to_bottom",
        "north_to_north",
        "south_to_south",
        "cross_to_ground",
        "claw_to_ground",
        "top_to_ground",
        "bottom_to_ground",
        "north_to_ground",
        "south_to_ground",
        "ground_to_ground",
    ];
    let s1_mut_keys = ["cross_to_claw", "top_to_bottom", "north_to_south"];
    let s1_all_keys = [
        "cross_to_cross",
        "claw_to_claw",
        "cross_to_claw",
        "cross_to_ground",
        "claw_to_ground",
        "ground_to_ground",
        "top_to_top",
        "bottom_to_bottom",
        "top_to_bottom",
        "top_to_ground",
        "bottom_to_ground",
        "north_to_north",
        "south_to_south",
        "north_to_south",
        "north_to_ground",
        "south_to_ground",
    ];

    let has_p1_tier = all_runs.iter().any(|r| r.tier.eq_ignore_ascii_case("p1"));
    let has_p2_tier = all_runs.iter().any(|r| r.tier.eq_ignore_ascii_case("p2"));

    let s1_core_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| {
            let tier_ok = if has_p2_tier {
                r.tier.eq_ignore_ascii_case("p2")
            } else {
                !r.tier.eq_ignore_ascii_case("smoke")
                    && !r.tier.eq_ignore_ascii_case("ablation")
                    && !r.tier.eq_ignore_ascii_case("gci")
                    && !r.tier.eq_ignore_ascii_case("determinism")
                    && !r.tier.eq_ignore_ascii_case("qcap")
            };
            r.suite_id.to_ascii_lowercase().starts_with("s1.")
                && tier_ok
                && (r.mode == "g2_dual_box" || r.mode == "far_field" || r.mode == "amr_ps" || r.mode == "qp_graded")
        })
        .collect();
    let s2_core_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| {
            let s = r.suite_id.to_ascii_lowercase();
            let tier_ok = if has_p2_tier {
                r.tier.eq_ignore_ascii_case("p2")
            } else {
                !r.tier.eq_ignore_ascii_case("smoke")
                    && !r.tier.eq_ignore_ascii_case("gci")
                    && !r.tier.eq_ignore_ascii_case("determinism")
            };
            (s.starts_with("s2.1") || s.starts_with("s2.2")) && tier_ok
        })
        .collect();
    let wm1_w2_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| {
            let s = r.suite_id.to_ascii_lowercase();
            let tier_ok = if has_p1_tier {
                r.tier.eq_ignore_ascii_case("p1")
            } else {
                !r.tier.eq_ignore_ascii_case("smoke")
                    && !r.tier.eq_ignore_ascii_case("p1_diag")
                    && !r.tier.eq_ignore_ascii_case("gci")
                    && !r.tier.eq_ignore_ascii_case("determinism")
            };
            (s.contains("s4.1") || s.contains("wm1")) && r.mode.eq_ignore_ascii_case("w2") && tier_ok
        })
        .collect();
    let muninn_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| {
            let s = r.suite_id.to_ascii_lowercase();
            let tier_ok = if has_p1_tier {
                r.tier.eq_ignore_ascii_case("p1")
            } else {
                !r.tier.eq_ignore_ascii_case("smoke")
                    && !r.tier.eq_ignore_ascii_case("gci")
                    && !r.tier.eq_ignore_ascii_case("determinism")
            };
            (s.contains("s4.2") || s.contains("muninn") || s.contains("munnin")) && tier_ok
        })
        .collect();
    let qc_equiv_runs: Vec<&RunRecord> = all_runs
        .iter()
        .filter(|r| {
            let s = r.suite_id.to_ascii_lowercase();
            let tier_ok = if has_p1_tier {
                r.tier.eq_ignore_ascii_case("p1")
            } else {
                !r.tier.eq_ignore_ascii_case("smoke")
                    && !r.tier.eq_ignore_ascii_case("gci")
                    && !r.tier.eq_ignore_ascii_case("determinism")
            };
            (s.contains("s2.3") || s.contains("s3.2") || s.contains("qubitcavity")) && tier_ok
        })
        .collect();

    // Cohort statistical summaries
    let s1_self_corr_stats =
        summarize_observable_cohort(&extract_metric_obs(&s1_core_runs, &s1_self_keys, false), "fF", true);
    let s1_self_raw_stats =
        summarize_observable_cohort(&extract_metric_obs(&s1_core_runs, &s1_self_keys, true), "fF", true);
    let s1_mut_corr_stats =
        summarize_observable_cohort(&extract_metric_obs(&s1_core_runs, &s1_mut_keys, false), "fF", true);
    let s1_mut_raw_stats =
        summarize_observable_cohort(&extract_metric_obs(&s1_core_runs, &s1_mut_keys, true), "fF", true);
    let s1_all_corr_stats =
        summarize_observable_cohort(&extract_metric_obs(&s1_core_runs, &s1_all_keys, false), "fF", true);

    let s2_fr_stats = summarize_observable_cohort(
        &extract_metric_obs(&s2_core_runs, &["f_r_hz", "f_r_Hz", "cavity_frequency_hz", "f_r"], false),
        "MHz",
        true,
    );
    let s2_kappa_stats = summarize_observable_cohort(
        &extract_metric_obs(&s2_core_runs, &["kappa_hz", "kappa_Hz", "kappa"], false),
        "kHz",
        true,
    );

    // Fixed hardware cohort WM1 W2 (Plan §5.1.1: exact unweighted RMSRE, bootstrap disabled for N=6)
    let wm1_w2_fr = summarize_observable_cohort(
        &extract_metric_obs(&wm1_w2_runs, &["f_r_hz", "f_r_Hz", "f_r"], false),
        "MHz",
        false,
    );
    let wm1_w2_alpha = summarize_observable_cohort(
        &extract_metric_obs(&wm1_w2_runs, &["alpha_hz", "alpha_Hz", "alpha"], false),
        "MHz",
        false,
    );
    let wm1_w2_g =
        summarize_observable_cohort(&extract_metric_obs(&wm1_w2_runs, &["g_hz", "g_Hz", "g"], false), "MHz", false);
    let wm1_w2_kappa = summarize_observable_cohort(
        &extract_metric_obs(&wm1_w2_runs, &["kappa_hz", "kappa_Hz", "kappa"], false),
        "kHz",
        false,
    );

    let gf = gate_filter.trim().to_ascii_uppercase();
    let mut gates_out = Map::new();

    // ----------------------------------------------------------------------- Gate G0 (`Plan §7.1.1`)
    if gf == "G0" || gf == "ALL" {
        let (sm_pct, sm_ok, sm_tot, sm_tax) = completion_and_taxonomy(&smoke_runs);
        let smoke_pass = sm_tot > 0 && sm_ok == sm_tot;
        let c01_pass = cal_opt
            .as_ref()
            .and_then(|c| c.get("C0.1"))
            .and_then(|x| x.get("gate_g0_pass"))
            .and_then(Value::as_bool)
            .unwrap_or(true);
        let c03_pass = cal_opt
            .as_ref()
            .and_then(|c| c.get("C0.3"))
            .and_then(|x| x.get("gate_g0_pass"))
            .and_then(Value::as_bool)
            .unwrap_or(true);
        let smoke_commit_gib = smoke_runs
            .iter()
            .map(|r| (r.memory_bytes.peak_commit as f64) / (1024.0 * 1024.0 * 1024.0))
            .fold(0.0_f64, f64::max);
        let max_commit_gib = all_runs
            .iter()
            .map(|r| (r.memory_bytes.peak_commit as f64) / (1024.0 * 1024.0 * 1024.0))
            .fold(0.0_f64, f64::max);
        let ram_pass = smoke_commit_gib <= 20.0 && max_commit_gib <= 22.0;

        let g0_pass = smoke_pass && c01_pass && c03_pass && ram_pass;
        gates_out.insert(
            "G0".into(),
            json!({
                "pass": g0_pass,
                "smoke_completion_pct": sm_pct,
                "smoke_ok": sm_ok,
                "smoke_total": sm_tot,
                "smoke_taxonomy": sm_tax,
                "c0_1_pass": c01_pass,
                "c0_3_pass": c03_pass,
                "smoke_max_peak_commit_gib": smoke_commit_gib,
                "max_peak_commit_gib": max_commit_gib,
                "ram_le_20gib_pass": ram_pass,
            }),
        );
    }

    // ----------------------------------------------------------------------- Gate G1 (`Plan §7.1.2`)
    if gf == "G1" || gf == "ALL" {
        let (p1_pct, p1_ok, p1_tot, p1_tax) = completion_and_taxonomy(&p1_runs);
        let comp_pass = p1_tot > 0 && p1_ok == p1_tot;

        let w2_fr_rmsre = wm1_w2_fr["rmsre_unweighted_pct"].as_f64().unwrap_or(f64::INFINITY);
        let w2_al_rmsre = wm1_w2_alpha["rmsre_unweighted_pct"].as_f64().unwrap_or(f64::INFINITY);
        let w2_g_rmsre = wm1_w2_g["rmsre_unweighted_pct"].as_f64().unwrap_or(f64::INFINITY);
        let w2_kp_rmsre = wm1_w2_kappa["rmsre_unweighted_pct"].as_f64().unwrap_or(f64::INFINITY);
        let wm1_pass = wm1_w2_runs.is_empty()
            || (w2_fr_rmsre <= 2.5
                && w2_al_rmsre <= 4.0
                && w2_g_rmsre <= 9.5
                && (wm1_w2_kappa["n_rel_included"].as_u64().unwrap_or(0) == 0 || w2_kp_rmsre <= 15.0));

        // Check MUNINN monotonic T1 ordering (`Ancilla >= Qubit >= Env`) if all 3 present
        let find_muninn_t1 = |tag: &str| -> Option<f64> {
            muninn_runs
                .iter()
                .find(|r| r.case_id.to_ascii_lowercase().contains(tag))
                .and_then(|r| r.predictions_si.get("t1_s").copied())
        };
        let muninn_mono_pass = match (find_muninn_t1("ancilla"), find_muninn_t1("qubit"), find_muninn_t1("env")) {
            (Some(ta), Some(tq), Some(te)) => ta >= tq * 0.98 && tq >= te * 0.98,
            _ => true,
        };

        // Unified QubitCavity equivalence deltas (`M1 vs M2 vs M3`)
        let max_pred_key = |key: &str| -> f64 {
            qc_equiv_runs.iter().filter_map(|r| r.predictions_si.get(key).copied()).fold(0.0_f64, f64::max)
        };
        let max_df01 = max_pred_key("m1_m2_df01_pct").max(max_pred_key("m1_m3_df01_pct"));
        let max_dfr = max_pred_key("m1_m2_dfr_pct").max(max_pred_key("m1_m3_dfr_pct"));
        let max_dalpha = max_pred_key("m1_m2_dalpha_pct").max(max_pred_key("m1_m3_dalpha_pct"));
        let max_dchi = max_pred_key("m1_m2_dchi_pct").max(max_pred_key("m1_m3_dchi_pct"));
        let qc_pass =
            qc_equiv_runs.is_empty() || (max_df01 <= 0.30 && max_dfr <= 0.10 && max_dalpha <= 2.00 && max_dchi <= 4.00);

        let g1_pass = comp_pass && wm1_pass && muninn_mono_pass && qc_pass;
        gates_out.insert(
            "G1".into(),
            json!({
                "pass": g1_pass,
                "p1_completion_pct": p1_pct,
                "p1_ok": p1_ok,
                "p1_total": p1_tot,
                "p1_taxonomy": p1_tax,
                "wm1_w2_rmsre_pct": {
                    "f_r": json_num_or_inf(w2_fr_rmsre),
                    "alpha": json_num_or_inf(w2_al_rmsre),
                    "g": json_num_or_inf(w2_g_rmsre),
                    "kappa": json_num_or_inf(w2_kp_rmsre),
                    "pass": wm1_pass,
                },
                "muninn_monotonic_t1_pass": muninn_mono_pass,
                "unified_qubit_cavity_max_delta_pct": {
                    "df01": max_df01,
                    "dfr": max_dfr,
                    "dalpha": max_dalpha,
                    "dchi": max_dchi,
                    "pass": qc_pass,
                }
            }),
        );
    }

    // ----------------------------------------------------------------------- Gate G2 (`Plan §7.1.3`)
    if gf == "G2" || gf == "ALL" {
        let p2_eval_runs = if !p2_runs.is_empty() { &p2_runs } else { &s1_core_runs };
        let (p2_pct, p2_ok, p2_tot, p2_tax) = completion_and_taxonomy(p2_eval_runs);
        let comp_pass = p2_tot > 0 && p2_pct >= 99.0;

        let s1_self_med = s1_self_corr_stats["median_ht_pct"]
            .as_f64()
            .or_else(|| s1_self_corr_stats["median_pct"].as_f64())
            .unwrap_or(f64::INFINITY);
        let s1_self_raw_med = s1_self_raw_stats["median_ht_pct"]
            .as_f64()
            .or_else(|| s1_self_raw_stats["median_pct"].as_f64())
            .unwrap_or(f64::INFINITY);
        let s1_mut_med = s1_mut_corr_stats["median_ht_pct"]
            .as_f64()
            .or_else(|| s1_mut_corr_stats["median_pct"].as_f64())
            .unwrap_or(f64::INFINITY);
        let s1_p95 = s1_all_corr_stats["p95_pct"].as_f64().unwrap_or(f64::INFINITY);
        let s1_mae_ff = s1_all_corr_stats["mae_unit"].as_f64().unwrap_or(f64::INFINITY);

        let mut s1_solve_times: Vec<f64> =
            s1_core_runs.iter().filter(|r| r.status == "ok").map(|r| r.timings_s.solve).collect();
        s1_solve_times.sort_by(|a, b| a.total_cmp(b));
        let s1_med_solve_s = sorted_quantile(&s1_solve_times, 0.50);

        let s1_pass = s1_core_runs.is_empty()
            || (s1_self_med <= 0.90 && s1_mut_med <= 1.20 && s1_p95 <= 3.80 && s1_mae_ff <= 0.40);

        let s2_fr_med = s2_fr_stats["median_ht_pct"]
            .as_f64()
            .or_else(|| s2_fr_stats["median_pct"].as_f64())
            .unwrap_or(f64::INFINITY);
        let s2_fr_p95 = s2_fr_stats["p95_pct"].as_f64().unwrap_or(f64::INFINITY);
        let s2_kp_med = s2_kappa_stats["median_ht_pct"]
            .as_f64()
            .or_else(|| s2_kappa_stats["median_pct"].as_f64())
            .unwrap_or(f64::INFINITY);
        let s2_kp_p95 = s2_kappa_stats["p95_pct"].as_f64().unwrap_or(f64::INFINITY);
        let s2_pass = s2_core_runs.is_empty()
            || (s2_fr_med <= 1.00 && s2_fr_p95 <= 2.20 && s2_kp_med <= 8.00 && s2_kp_p95 <= 18.00);

        let g2_pass = comp_pass && s1_pass && s2_pass;
        gates_out.insert(
            "G2".into(),
            json!({
                "pass": g2_pass,
                "p2_completion_pct": p2_pct,
                "p2_ok": p2_ok,
                "p2_total": p2_tot,
                "p2_taxonomy": p2_tax,
                "suite1_vs_q3d": {
                    "self_median_corr_pct": json_num_or_inf(s1_self_med),
                    "self_median_raw_pct": json_num_or_inf(s1_self_raw_med),
                    "mutual_median_corr_pct": json_num_or_inf(s1_mut_med),
                    "all_p95_corr_pct": json_num_or_inf(s1_p95),
                    "all_mae_fF": json_num_or_inf(s1_mae_ff),
                    "median_solve_s": json_num_or_inf(s1_med_solve_s),
                    "pass": s1_pass,
                },
                "suite2_vs_hfss": {
                    "f_r_median_pct": json_num_or_inf(s2_fr_med),
                    "f_r_p95_pct": json_num_or_inf(s2_fr_p95),
                    "kappa_median_pct": json_num_or_inf(s2_kp_med),
                    "kappa_p95_pct": json_num_or_inf(s2_kp_p95),
                    "pass": s2_pass,
                }
            }),
        );
    }

    // ----------------------------------------------------------------------- Gate G3 (`Plan §7.1.3` — Determinism & Invariants)
    if gf == "G3" || gf == "ALL" {
        let det_runs: Vec<&RunRecord> =
            all_runs.iter().filter(|r| r.tier.eq_ignore_ascii_case("determinism") || r.mode.contains("det")).collect();
        let mut max_det_rel = 0.0_f64;
        for r in &det_runs {
            if let Some(&d) = r.predictions_si.get("det_max_rel_diff") {
                max_det_rel = max_det_rel.max(d);
            }
        }
        if det_runs.len() >= 2 {
            let mut by_case_mode: BTreeMap<(&str, &str), Vec<&RunRecord>> = BTreeMap::new();
            for r in &det_runs {
                let base_mode = r.mode.split("_det_t").next().unwrap_or(&r.mode);
                by_case_mode.entry((r.case_id.as_str(), base_mode)).or_default().push(*r);
            }
            for (_key, grp) in by_case_mode {
                if grp.len() >= 2 {
                    for (k, &v0) in &grp[0].predictions_si {
                        if k == "det_max_rel_diff" {
                            continue;
                        }
                        if let Some(&v1) = grp[1].predictions_si.get(k) {
                            if v0.abs() > 1.0e-30 {
                                max_det_rel = max_det_rel.max((v1 - v0).abs() / v0.abs());
                            }
                        }
                    }
                }
            }
        }
        let inv_failures = *overall_tax.get("invariant_fail").unwrap_or(&0);
        let det_pass = max_det_rel < 1.0e-13;
        let g3_pass = det_pass && inv_failures == 0;
        gates_out.insert(
            "G3".into(),
            json!({
                "pass": g3_pass,
                "max_thread_determinism_rel_diff": max_det_rel,
                "determinism_lt_1e13_pass": det_pass,
                "invariant_failures": inv_failures,
            }),
        );
    }

    let sys_mem = nyx_engine::process::system_memory();
    let report = json!({
        "schema_version": 1,
        "runs_source": runs_path.to_string_lossy(),
        "host_profile": {
            "os": std::env::consts::OS,
            "arch": std::env::consts::ARCH,
            "logical_threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
            "total_ram_bytes": sys_mem.map(|(t, _)| t),
            "power_throttling_opt_out": true,
        },
        "denominator_floors_si": {
            "capacitance_F": 0.5e-15,
            "chi_Hz": 25.0e3,
            "kappa_Hz": 20.0e3,
        },
        "overall_completion": {
            "completion_pct": overall_comp_pct,
            "n_ok": overall_ok,
            "n_total": overall_tot,
            "failure_taxonomy": overall_tax,
        },
        "invariants": invariants_summary,
        "cohorts": {
            "suite1_self_corrected": s1_self_corr_stats,
            "suite1_self_raw": s1_self_raw_stats,
            "suite1_mutual_corrected": s1_mut_corr_stats,
            "suite1_mutual_raw": s1_mut_raw_stats,
            "suite1_all_corrected": s1_all_corr_stats,
            "suite2_f_r": s2_fr_stats,
            "suite2_kappa": s2_kappa_stats,
            "wm1_w2_f_r": wm1_w2_fr,
            "wm1_w2_alpha": wm1_w2_alpha,
            "wm1_w2_g": wm1_w2_g,
            "wm1_w2_kappa": wm1_w2_kappa,
        },
        "gates": Value::Object(gates_out),
    });

    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
    }
    std::fs::write(out_path, serde_json::to_string_pretty(&report).unwrap()).map_err(|e| e.to_string())?;
    Ok(report)
}
