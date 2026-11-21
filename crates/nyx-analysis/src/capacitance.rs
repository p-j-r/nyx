// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Multi-terminal capacitance (`capacitance_matrix`) and inductance (`inductance_matrix`) analyses of `dfs` results.
//!
//! `CapacitanceMatrixAnalysis` reads `CapacitanceMatrixData` from `SolveResult::scalars` (produced by an electrostatic
//! study of `dfs` with one or more `terminal` boundary conditions) and returns the Maxwell capacitance matrix (F), the
//! mutual/SPICE capacitance matrix (F), and the charging energy matrix `E_C / h = (e^2 / 2h) C_M^(-1)` (Hz).
//!
//! When optional `params` are supplied to [`capacitance_matrix_with_params`] (or `Analysis::run_with_params`), it also
//! performs:
//! - Ground-terminal projection (`"eliminate": "drop"`) or floating-reference Schur elimination
//!   (`"eliminate": "floating_schur"`, `C_ss - C_sg C_gg^-1 C_gs`), reporting `reduced_terminals`,
//!   `maxwell_reduced_F`, `mutual_reduced_F`, `ec_reduced_Hz`, `condition_number`, `null_modes`, and `has_null_mode`;
//! - Single-pad grounded (`Xmon`) and differential two-pad floating (`TransmonPocket`, Eqs. 17–18 of
//!   `arXiv:2611.22057v1`) lumped qubit extraction (`C_Sigma_F`, `E_C_Hz`, `C_eff_F`, `beta`, and transmon `ej_over_ec`,
//!   `f01_Hz`, `anharmonicity_Hz`, `n_zpf` when `E_J_Hz` is given);
//! - Quarter-wave, half-wave, and lumped readout/bus resonator ZPF voltage (`C_r_F`, `V_rms_V`) and pairwise
//!   capacitive coupling `g_Hz` and dispersive shift `chi_Hz`.
//!
//! `InductanceMatrixAnalysis` reads `InductanceMatrixData` from `SolveResult::scalars` (produced by a magnetostatic
//! study of `dfs`) and returns the multi-terminal inductance matrix `L` (H, or H/m for 2D `cross_section` models) and
//! the active magnetostatic formulation (`"cross_section"`, `"planar_stream"`, or `"mixed_3d"`).

use crate::epr::{E_CHARGE, PLANCK_H};
use crate::Analysis;
use femcore::linalg::driver::dense_eigh;
use nyx_physics::{CapacitanceMatrixData, Error, InductanceMatrixData, Mesh, SolveResult};
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

/// Relative eigenvalue floor below which a symmetric capacitance mode is classified as a gauge/floating null mode.
const NULL_MODE_REL_TOL: f64 = 1.0e-10;
/// Finite ceiling used for `condition_number` when an exact zero/floating null mode is present (keeping JSON finite).
const COND_CAP: f64 = 1.0e16;

/// Extracts the typed [`CapacitanceMatrixData`] payload from `result`.
pub fn capacitance_matrix(result: &SolveResult) -> Result<CapacitanceMatrixData, Error> {
    CapacitanceMatrixData::of(result).ok_or_else(|| {
        Error::invalid(
            "the result carries no capacitance matrix (run an electrostatic study of dfs with terminal conditions)",
        )
    })
}

/// True when `result` carries a multi-terminal capacitance matrix payload.
pub fn capacitance_matrix_applies(result: &SolveResult) -> bool {
    CapacitanceMatrixData::of(result).is_some()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EliminateMode {
    Drop,
    FloatingSchur,
}

impl EliminateMode {
    fn as_str(self) -> &'static str {
        match self {
            Self::Drop => "drop",
            Self::FloatingSchur => "floating_schur",
        }
    }
}

struct ComputedQubit {
    name: String,
    pads: Vec<String>,
    c_sigma_f: f64,
    e_c_hz: f64,
    c_eff_f: BTreeMap<String, f64>,
    beta: BTreeMap<String, f64>,
    e_j_hz: Option<f64>,
    f01_hz: Option<f64>,
    n_zpf: Option<f64>,
    json: Value,
}

struct ComputedResonator {
    name: String,
    terminal: String,
    f_r_hz: f64,
    v_rms_v: f64,
    json: Value,
}

fn mutual_from_maxwell(maxwell: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = maxwell.len();
    let mut mutual = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut row_sum = 0.0;
        for (k, &cik) in maxwell[i].iter().enumerate().take(n) {
            row_sum += cik;
            if i != k {
                mutual[i][k] = -cik;
            }
        }
        mutual[i][i] = row_sum;
    }
    mutual
}

/// Eigendecomposes symmetric PSD `maxwell` (`m x m`), returning `(ec_matrix_hz, condition_number, null_modes)`.
#[allow(clippy::needless_range_loop)]
fn spectral_ec_and_cond(maxwell: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, f64, usize), Error> {
    let m = maxwell.len();
    if m == 0 {
        return Ok((Vec::new(), 1.0, 0));
    }
    let mut flat = Vec::with_capacity(m * m);
    for row in maxwell {
        flat.extend_from_slice(row);
    }
    let (eigs, vecs) = dense_eigh(&flat, m).map_err(Error::from_kernel)?;
    let lam_max = eigs.iter().copied().fold(0.0_f64, f64::max);
    let lam_min = eigs.iter().copied().fold(f64::INFINITY, f64::min);
    if lam_max <= 0.0 || lam_min < -1.0e-8 * lam_max {
        return Err(Error::invalid("reduced Maxwell capacitance matrix is not positive semi-definite"));
    }

    let null_thresh = NULL_MODE_REL_TOL * lam_max;
    let null_modes = eigs.iter().filter(|&&lam| lam <= null_thresh).count();
    let denom = lam_min.max(lam_max / COND_CAP);
    let condition_number = (lam_max / denom).min(COND_CAP);

    let scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);
    let mut ec = vec![vec![0.0; m]; m];
    for (k, &lam) in eigs.iter().enumerate().take(m) {
        if lam > null_thresh {
            let inv_l = scale / lam;
            let vk = &vecs[k];
            for i in 0..m {
                for j in 0..m {
                    ec[i][j] += inv_l * vk[i] * vk[j];
                }
            }
        }
    }
    for i in 0..m {
        for j in (i + 1)..m {
            let sym = 0.5 * (ec[i][j] + ec[j][i]);
            ec[i][j] = sym;
            ec[j][i] = sym;
        }
    }
    Ok((ec, condition_number, null_modes))
}

fn parse_ground_names(params_map: &Map<String, Value>) -> Result<Vec<String>, Error> {
    match params_map.get("ground") {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::String(s)) => Ok(vec![s.clone()]),
        Some(Value::Array(arr)) => {
            let mut out = Vec::with_capacity(arr.len());
            for (i, item) in arr.iter().enumerate() {
                let s = item
                    .as_str()
                    .ok_or_else(|| Error::invalid(format!("ground[{i}] must be a terminal name string")))?;
                out.push(s.to_string());
            }
            Ok(out)
        }
        _ => Err(Error::invalid("parameter 'ground' must be a terminal name string or array of strings")),
    }
}

fn parse_eliminate_mode(params_map: &Map<String, Value>) -> Result<EliminateMode, Error> {
    match params_map.get("eliminate") {
        None | Some(Value::Null) => Ok(EliminateMode::Drop),
        Some(Value::String(s)) => match s.as_str() {
            "drop" => Ok(EliminateMode::Drop),
            "floating_schur" => Ok(EliminateMode::FloatingSchur),
            other => {
                Err(Error::invalid(format!("parameter 'eliminate' must be 'drop' or 'floating_schur', not '{other}'")))
            }
        },
        _ => Err(Error::invalid("parameter 'eliminate' must be 'drop' or 'floating_schur'")),
    }
}

fn parse_qubits(
    q_val: &Value,
    all_terminals: &[String],
    ground_set: &BTreeSet<String>,
    reduced_terminals: &[String],
    red_idx_of: &BTreeMap<&str, usize>,
    c_red: &[Vec<f64>],
) -> Result<Vec<ComputedQubit>, Error> {
    let arr = q_val.as_array().ok_or_else(|| Error::invalid("parameter 'qubits' must be an array of objects"))?;
    let mut seen_names = BTreeSet::new();
    let mut out = Vec::with_capacity(arr.len());
    let ec_scale = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H);

    for (qi, item) in arr.iter().enumerate() {
        let obj = item.as_object().ok_or_else(|| Error::invalid(format!("qubits[{qi}] must be an object")))?;
        for k in obj.keys() {
            if !matches!(k.as_str(), "name" | "pads" | "E_J_Hz") {
                return Err(Error::invalid(format!(
                    "unknown field '{k}' in qubits[{qi}] (expected: name, pads, E_J_Hz)"
                )));
            }
        }
        let name = obj
            .get("name")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| Error::invalid(format!("qubits[{qi}].name must be a non-empty string")))?
            .to_string();
        if !seen_names.insert(name.clone()) {
            return Err(Error::invalid(format!("duplicate qubit name '{name}'")));
        }

        let pads_arr = obj
            .get("pads")
            .and_then(Value::as_array)
            .ok_or_else(|| Error::invalid(format!("qubit '{name}' requires 'pads' as [str] or [str, str]")))?;
        if !(1..=2).contains(&pads_arr.len()) {
            return Err(Error::invalid(format!(
                "qubit '{name}' 'pads' must have 1 (grounded) or 2 (differential) terminal names, got {}",
                pads_arr.len()
            )));
        }
        let mut pads = Vec::with_capacity(pads_arr.len());
        let mut pad_indices = Vec::with_capacity(pads_arr.len());
        for (pi, pv) in pads_arr.iter().enumerate() {
            let p_name =
                pv.as_str().ok_or_else(|| Error::invalid(format!("qubit '{name}' pads[{pi}] must be a string")))?;
            if !all_terminals.iter().any(|t| t == p_name) {
                return Err(Error::invalid(format!(
                    "qubit '{name}' pad '{p_name}' does not match any terminal in {all_terminals:?}"
                )));
            }
            if ground_set.contains(p_name) {
                return Err(Error::invalid(format!(
                    "qubit '{name}' pad '{p_name}' cannot also be an eliminated ground terminal"
                )));
            }
            let &idx = red_idx_of.get(p_name).ok_or_else(|| {
                Error::invalid(format!("qubit '{name}' pad '{p_name}' is not present in reduced terminals"))
            })?;
            pads.push(p_name.to_string());
            pad_indices.push(idx);
        }
        if pad_indices.len() == 2 && pad_indices[0] == pad_indices[1] {
            return Err(Error::invalid(format!("qubit '{name}' two-pad specification must name two distinct pads")));
        }

        let e_j_hz = match obj.get("E_J_Hz") {
            None | Some(Value::Null) => None,
            Some(v) => Some(v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("qubit '{name}' E_J_Hz must be a finite positive frequency in Hz"))
            })?),
        };

        let mut c_eff_map = BTreeMap::new();
        let mut beta_map = BTreeMap::new();
        let (topology, c_sigma_f, extra_pad_caps) = if pad_indices.len() == 1 {
            let s = pad_indices[0];
            let c_ss = c_red[s][s];
            if !(c_ss > 0.0 && c_ss.is_finite()) {
                return Err(Error::invalid(format!("qubit '{name}' self-capacitance C_ss must be positive")));
            }
            for (k_idx, term_k) in reduced_terminals.iter().enumerate() {
                if k_idx == s {
                    continue;
                }
                let c_eff_k = c_red[s][k_idx].abs();
                c_eff_map.insert(term_k.clone(), c_eff_k);
                beta_map.insert(term_k.clone(), c_eff_k / c_ss);
            }
            ("grounded_1pad", c_ss, None)
        } else {
            let (s, b) = (pad_indices[0], pad_indices[1]);
            let c_ss = c_red[s][s];
            let c_bb = c_red[b][b];
            let c_sb = c_red[s][b];
            if !(c_ss > 0.0 && c_bb > 0.0) {
                return Err(Error::invalid(format!("qubit '{name}' pad self-capacitances must be positive")));
            }
            let sum_diag = c_ss + c_bb;
            let denom = sum_diag + 2.0 * c_sb;
            let det = c_ss * c_bb - c_sb * c_sb;
            let c_sigma = if denom > 1.0e-14 * sum_diag {
                if det <= 0.0 {
                    return Err(Error::invalid(format!(
                        "qubit '{name}' 2-pad Maxwell submatrix is not positive definite"
                    )));
                }
                det / denom
            } else if denom >= -1.0e-12 * sum_diag && c_sb < 0.0 {
                // Pure isolated two-pad floating dipole limit (C_ss = C_bb = -C_sb, zero shunt to ground)
                -c_sb
            } else {
                return Err(Error::invalid(format!(
                    "qubit '{name}' 2-pad Maxwell entries (C_ss={c_ss:e}, C_bb={c_bb:e}, C_sb={c_sb:e}) violate passivity"
                )));
            };

            for (k_idx, term_k) in reduced_terminals.iter().enumerate() {
                if k_idx == s || k_idx == b {
                    continue;
                }
                let c_sk = c_red[s][k_idx].abs();
                let c_bk = c_red[b][k_idx].abs();
                let c_eff_k = if denom > 1.0e-14 * sum_diag {
                    ((c_bb + c_sb) * c_sk - (c_ss + c_sb) * c_bk).abs() / denom
                } else {
                    0.5 * (c_sk - c_bk).abs()
                };
                c_eff_map.insert(term_k.clone(), c_eff_k);
                beta_map.insert(term_k.clone(), c_eff_k / c_sigma);
            }
            ("differential_2pad", c_sigma, Some((c_ss, c_bb, c_sb)))
        };

        let e_c_hz = ec_scale / c_sigma_f;
        let mut q_json = Map::new();
        q_json.insert("name".into(), json!(name));
        q_json.insert("topology".into(), json!(topology));
        q_json.insert("pads".into(), json!(pads));
        if let Some((c_ss, c_bb, c_sb)) = extra_pad_caps {
            q_json.insert("C_ss_F".into(), json!(c_ss));
            q_json.insert("C_bb_F".into(), json!(c_bb));
            q_json.insert("C_sb_F".into(), json!(c_sb));
        }
        q_json.insert("C_Sigma_F".into(), json!(c_sigma_f));
        q_json.insert("E_C_Hz".into(), json!(e_c_hz));
        q_json.insert("C_eff_F".into(), json!(c_eff_map));
        q_json.insert("beta".into(), json!(beta_map));

        let (mut f01_opt, mut nzpf_opt) = (None, None);
        if let Some(ej) = e_j_hz {
            let ej_over_ec = ej / e_c_hz;
            let f01 = (8.0 * ej * e_c_hz).sqrt() - e_c_hz;
            let alpha = -e_c_hz;
            let nzpf = (ej / (32.0 * e_c_hz)).powf(0.25);
            q_json.insert("E_J_Hz".into(), json!(ej));
            q_json.insert("ej_over_ec".into(), json!(ej_over_ec));
            q_json.insert("f01_Hz".into(), json!(f01));
            q_json.insert("anharmonicity_Hz".into(), json!(alpha));
            q_json.insert("n_zpf".into(), json!(nzpf));
            f01_opt = Some(f01);
            nzpf_opt = Some(nzpf);
        }

        out.push(ComputedQubit {
            name,
            pads,
            c_sigma_f,
            e_c_hz,
            c_eff_f: c_eff_map,
            beta: beta_map,
            e_j_hz,
            f01_hz: f01_opt,
            n_zpf: nzpf_opt,
            json: Value::Object(q_json),
        });
    }
    Ok(out)
}

fn parse_resonators(
    r_val: &Value,
    all_terminals: &[String],
    ground_set: &BTreeSet<String>,
    red_idx_of: &BTreeMap<&str, usize>,
) -> Result<Vec<ComputedResonator>, Error> {
    let arr = r_val.as_array().ok_or_else(|| Error::invalid("parameter 'resonators' must be an array of objects"))?;
    let mut seen_names = BTreeSet::new();
    let mut out = Vec::with_capacity(arr.len());
    let hbar = PLANCK_H / (2.0 * PI);

    for (ri, item) in arr.iter().enumerate() {
        let obj = item.as_object().ok_or_else(|| Error::invalid(format!("resonators[{ri}] must be an object")))?;
        for k in obj.keys() {
            if !matches!(k.as_str(), "name" | "terminal" | "f_r_Hz" | "Z0_ohm" | "resonator_type" | "C_r_F") {
                return Err(Error::invalid(format!(
                    "unknown field '{k}' in resonators[{ri}] (expected: name, terminal, f_r_Hz, Z0_ohm, resonator_type, C_r_F)"
                )));
            }
        }
        let name = obj
            .get("name")
            .and_then(Value::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .ok_or_else(|| Error::invalid(format!("resonators[{ri}].name must be a non-empty string")))?
            .to_string();
        if !seen_names.insert(name.clone()) {
            return Err(Error::invalid(format!("duplicate resonator name '{name}'")));
        }

        let terminal = obj
            .get("terminal")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid(format!("resonator '{name}' requires 'terminal' string")))?
            .to_string();
        if !all_terminals.iter().any(|t| t == &terminal) {
            return Err(Error::invalid(format!(
                "resonator '{name}' terminal '{terminal}' does not match any terminal in {all_terminals:?}"
            )));
        }
        if ground_set.contains(&terminal) || !red_idx_of.contains_key(terminal.as_str()) {
            return Err(Error::invalid(format!(
                "resonator '{name}' terminal '{terminal}' cannot be an eliminated ground terminal"
            )));
        }

        let f_r_hz = obj
            .get("f_r_Hz")
            .and_then(Value::as_f64)
            .filter(|&x| x.is_finite() && x > 0.0)
            .ok_or_else(|| Error::invalid(format!("resonator '{name}' requires finite positive 'f_r_Hz'")))?;

        let z0_ohm = match obj.get("Z0_ohm") {
            None | Some(Value::Null) => 50.0,
            Some(v) => v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' 'Z0_ohm' must be a finite positive impedance in Ohm"))
            })?,
        };

        let resonator_type = match obj.get("resonator_type") {
            None | Some(Value::Null) => "quarter_wave",
            Some(Value::String(s)) => match s.as_str() {
                "quarter_wave" | "half_wave" | "lumped" => s.as_str(),
                other => {
                    return Err(Error::invalid(format!(
                        "resonator '{name}' 'resonator_type' must be 'quarter_wave', 'half_wave', or 'lumped', not '{other}'"
                    )))
                }
            },
            _ => return Err(Error::invalid(format!("resonator '{name}' 'resonator_type' must be a string"))),
        };

        let c_r_override = match obj.get("C_r_F") {
            None | Some(Value::Null) => None,
            Some(v) => Some(v.as_f64().filter(|&x| x.is_finite() && x > 0.0).ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' 'C_r_F' must be a finite positive capacitance in F"))
            })?),
        };

        let omega_r = 2.0 * PI * f_r_hz;
        let c_r_f = match resonator_type {
            "quarter_wave" => PI / (4.0 * omega_r * z0_ohm),
            "half_wave" => PI / (2.0 * omega_r * z0_ohm),
            "lumped" => c_r_override.ok_or_else(|| {
                Error::invalid(format!("resonator '{name}' with resonator_type='lumped' requires positive 'C_r_F'"))
            })?,
            _ => unreachable!(),
        };
        let v_rms_v = (hbar * omega_r / (2.0 * c_r_f)).sqrt();

        let r_json = json!({
            "name": name,
            "terminal": terminal,
            "resonator_type": resonator_type,
            "f_r_Hz": f_r_hz,
            "Z0_ohm": z0_ohm,
            "C_r_F": c_r_f,
            "V_rms_V": v_rms_v,
        });
        out.push(ComputedResonator { name, terminal, f_r_hz, v_rms_v, json: r_json });
    }
    Ok(out)
}

/// Runs the multi-terminal capacitance analysis with optional ground elimination, qubit reduction, and resonator
/// dispersive coupling parameters (`params`).
///
/// When `params` is `Value::Null` or an empty object `{}`, returns the exact 4-field base JSON
/// (`terminals`, `maxwell_matrix_F`, `mutual_matrix_F`, `ec_matrix_Hz`) unchanged.
pub fn capacitance_matrix_with_params(result: &SolveResult, params: &Value) -> Result<Value, Error> {
    let cap = capacitance_matrix(result)?;
    let params_map = match params {
        Value::Null => {
            return Ok(json!({
                "terminals": cap.terminals,
                "maxwell_matrix_F": cap.maxwell_matrix_f,
                "mutual_matrix_F": cap.mutual_matrix_f,
                "ec_matrix_Hz": cap.ec_matrix_hz,
            }));
        }
        Value::Object(m) if m.is_empty() => {
            return Ok(json!({
                "terminals": cap.terminals,
                "maxwell_matrix_F": cap.maxwell_matrix_f,
                "mutual_matrix_F": cap.mutual_matrix_f,
                "ec_matrix_Hz": cap.ec_matrix_hz,
            }));
        }
        Value::Object(m) => m,
        _ => return Err(Error::invalid("capacitance_matrix params must be a JSON object or null")),
    };

    for k in params_map.keys() {
        if !matches!(k.as_str(), "ground" | "eliminate" | "qubits" | "resonators") {
            return Err(Error::invalid(format!(
                "unknown capacitance_matrix parameter '{k}' (expected: ground, eliminate, qubits, resonators)"
            )));
        }
    }

    let n = cap.terminals.len();
    if n == 0 || cap.maxwell_matrix_f.len() != n || cap.maxwell_matrix_f.iter().any(|r| r.len() != n) {
        return Err(Error::invalid("CapacitanceMatrixData dimensions are empty or inconsistent"));
    }
    if cap.maxwell_matrix_f.iter().flatten().any(|x| !x.is_finite()) {
        return Err(Error::NonFinite("Maxwell capacitance matrix contains non-finite entries".into()));
    }

    let ground_raw = parse_ground_names(params_map)?;
    let eliminate = parse_eliminate_mode(params_map)?;
    let mut ground_idx_set = BTreeSet::new();
    for g in &ground_raw {
        let idx = cap.terminals.iter().position(|t| t == g).ok_or_else(|| {
            Error::invalid(format!("ground terminal '{g}' does not match any terminal in {:?}", cap.terminals))
        })?;
        ground_idx_set.insert(idx);
    }

    let ground_indices: Vec<usize> = ground_idx_set.iter().copied().collect();
    let ground_terminals: Vec<String> = ground_indices.iter().map(|&i| cap.terminals[i].clone()).collect();
    let ground_name_set: BTreeSet<String> = ground_terminals.iter().cloned().collect();
    let sig_indices: Vec<usize> = (0..n).filter(|i| !ground_idx_set.contains(i)).collect();
    if sig_indices.is_empty() {
        return Err(Error::invalid("ground cannot eliminate all terminals; at least one signal terminal must remain"));
    }

    let m = sig_indices.len();
    let reduced_terminals: Vec<String> = sig_indices.iter().map(|&i| cap.terminals[i].clone()).collect();
    let mut maxwell_reduced_f = vec![vec![0.0; m]; m];
    for a in 0..m {
        for b in a..m {
            let (ia, ib) = (sig_indices[a], sig_indices[b]);
            let sym = 0.5 * (cap.maxwell_matrix_f[ia][ib] + cap.maxwell_matrix_f[ib][ia]);
            maxwell_reduced_f[a][b] = sym;
            maxwell_reduced_f[b][a] = sym;
        }
    }

    if eliminate == EliminateMode::FloatingSchur && !ground_indices.is_empty() {
        let kg = ground_indices.len();
        let mut flat_gg = Vec::with_capacity(kg * kg);
        for &ig in &ground_indices {
            for &jg in &ground_indices {
                flat_gg.push(0.5 * (cap.maxwell_matrix_f[ig][jg] + cap.maxwell_matrix_f[jg][ig]));
            }
        }
        let (eigs_gg, vecs_gg) = dense_eigh(&flat_gg, kg).map_err(Error::from_kernel)?;
        let lam_max_gg = eigs_gg.iter().copied().fold(0.0_f64, f64::max);
        let lam_min_gg = eigs_gg.iter().copied().fold(f64::INFINITY, f64::min);
        if lam_max_gg <= 0.0 || lam_min_gg < -1.0e-8 * lam_max_gg {
            return Err(Error::invalid("ground submatrix C_gg is not positive semi-definite in floating_schur"));
        }
        let thresh_gg = NULL_MODE_REL_TOL * lam_max_gg;
        let mut c_sg = vec![vec![0.0; kg]; m];
        for (a, &ia) in sig_indices.iter().enumerate() {
            for (p, &ig) in ground_indices.iter().enumerate() {
                c_sg[a][p] = 0.5 * (cap.maxwell_matrix_f[ia][ig] + cap.maxwell_matrix_f[ig][ia]);
            }
        }
        for (k, &lam) in eigs_gg.iter().enumerate().take(kg) {
            if lam > thresh_gg {
                let vk = &vecs_gg[k];
                let w: Vec<f64> = (0..m).map(|a| (0..kg).map(|p| c_sg[a][p] * vk[p]).sum()).collect();
                for a in 0..m {
                    for b in a..m {
                        let sub = (w[a] * w[b]) / lam;
                        maxwell_reduced_f[a][b] -= sub;
                        if a != b {
                            maxwell_reduced_f[b][a] -= sub;
                        }
                    }
                }
            }
        }
    }

    let mutual_reduced_f = mutual_from_maxwell(&maxwell_reduced_f);
    let (ec_reduced_hz, condition_number, null_modes) = spectral_ec_and_cond(&maxwell_reduced_f)?;
    let has_null_mode = null_modes > 0;

    let red_idx_of: BTreeMap<&str, usize> =
        reduced_terminals.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();

    let qubits_opt = match params_map.get("qubits") {
        None | Some(Value::Null) => None,
        Some(v) => Some(parse_qubits(
            v,
            &cap.terminals,
            &ground_name_set,
            &reduced_terminals,
            &red_idx_of,
            &maxwell_reduced_f,
        )?),
    };

    let resonators_opt = match params_map.get("resonators") {
        None | Some(Value::Null) => None,
        Some(v) => Some(parse_resonators(v, &cap.terminals, &ground_name_set, &red_idx_of)?),
    };

    let mut out = Map::new();
    out.insert("terminals".into(), json!(cap.terminals));
    out.insert("maxwell_matrix_F".into(), json!(cap.maxwell_matrix_f));
    out.insert("mutual_matrix_F".into(), json!(cap.mutual_matrix_f));
    out.insert("ec_matrix_Hz".into(), json!(cap.ec_matrix_hz));
    out.insert("eliminate".into(), json!(eliminate.as_str()));
    out.insert("ground_terminals".into(), json!(ground_terminals));
    out.insert("reduced_terminals".into(), json!(reduced_terminals));
    out.insert("maxwell_reduced_F".into(), json!(maxwell_reduced_f));
    out.insert("mutual_reduced_F".into(), json!(mutual_reduced_f));
    out.insert("ec_reduced_Hz".into(), json!(ec_reduced_hz));
    out.insert("condition_number".into(), json!(condition_number));
    out.insert("null_modes".into(), json!(null_modes));
    out.insert("has_null_mode".into(), json!(has_null_mode));

    if let (Some(qubits), Some(resonators)) = (&qubits_opt, &resonators_opt) {
        let mut couplings = Vec::with_capacity(qubits.len() * resonators.len());
        for q in qubits {
            for r in resonators {
                if q.pads.iter().any(|p| p == &r.terminal) {
                    return Err(Error::invalid(format!(
                        "resonator '{}' terminal '{}' cannot be one of qubit '{}' pads {:?}",
                        r.name, r.terminal, q.name, q.pads
                    )));
                }
                let c_eff = q.c_eff_f.get(&r.terminal).copied().unwrap_or(0.0);
                let beta = q.beta.get(&r.terminal).copied().unwrap_or_else(|| {
                    if q.c_sigma_f > 0.0 {
                        c_eff / q.c_sigma_f
                    } else {
                        0.0
                    }
                });
                let mut c_obj = Map::new();
                c_obj.insert("qubit".into(), json!(q.name));
                c_obj.insert("resonator".into(), json!(r.name));
                c_obj.insert("terminal".into(), json!(r.terminal));
                c_obj.insert("C_eff_F".into(), json!(c_eff));
                c_obj.insert("beta".into(), json!(beta));
                c_obj.insert("V_rms_V".into(), json!(r.v_rms_v));
                if let (Some(_ej), Some(f01), Some(nzpf)) = (q.e_j_hz, q.f01_hz, q.n_zpf) {
                    let g_hz = 2.0 * E_CHARGE * beta * nzpf * r.v_rms_v / PLANCK_H;
                    let delta_hz = f01 - r.f_r_hz;
                    let chi_hz = if delta_hz.abs() > 1.0e-12 * r.f_r_hz.max(1.0)
                        && (delta_hz - q.e_c_hz).abs() > 1.0e-12 * q.e_c_hz.max(1.0)
                    {
                        Some(-(g_hz * g_hz) / delta_hz * (q.e_c_hz / (delta_hz - q.e_c_hz)))
                    } else {
                        None
                    };
                    c_obj.insert("n_zpf".into(), json!(nzpf));
                    c_obj.insert("g_Hz".into(), json!(g_hz));
                    c_obj.insert("Delta_Hz".into(), json!(delta_hz));
                    c_obj.insert("chi_Hz".into(), json!(chi_hz));
                }
                couplings.push(Value::Object(c_obj));
            }
        }
        out.insert("couplings".into(), Value::Array(couplings));
    }

    if let Some(qubits) = qubits_opt {
        out.insert("qubits".into(), Value::Array(qubits.into_iter().map(|q| q.json).collect()));
    }
    if let Some(resonators) = resonators_opt {
        out.insert("resonators".into(), Value::Array(resonators.into_iter().map(|r| r.json).collect()));
    }

    Ok(Value::Object(out))
}

/// Multi-terminal Maxwell, mutual, and charging-energy capacitance matrix analysis behind the [`Analysis`] trait.
pub struct CapacitanceMatrixAnalysis;

impl Analysis for CapacitanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "capacitance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal capacitance and charging energy matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        capacitance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        capacitance_matrix_with_params(result, &Value::Null)
    }

    fn run_with_params(&self, _mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        capacitance_matrix_with_params(result, params)
    }
}

/// Extracts the typed [`InductanceMatrixData`] payload from `result`.
pub fn inductance_matrix(result: &SolveResult) -> Result<InductanceMatrixData, Error> {
    InductanceMatrixData::of(result)
        .ok_or_else(|| Error::invalid("the result carries no inductance matrix (run a magnetostatic study of dfs)"))
}

/// True when `result` carries a magnetostatic inductance matrix payload.
pub fn inductance_matrix_applies(result: &SolveResult) -> bool {
    InductanceMatrixData::of(result).is_some()
}

/// Multi-terminal magnetostatic inductance matrix analysis behind the [`Analysis`] trait.
pub struct InductanceMatrixAnalysis;

impl Analysis for InductanceMatrixAnalysis {
    fn id(&self) -> &'static str {
        "inductance_matrix"
    }

    fn name(&self) -> &'static str {
        "Multi-terminal inductance matrix"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        inductance_matrix_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        let ind = inductance_matrix(result)?;
        Ok(json!({
            "terminals": ind.terminals,
            "inductance_matrix_H": ind.inductance_matrix_h,
            "formulation": ind.formulation,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacitance_and_inductance_analyses_round_trip_to_json() {
        let sm = femcore::testkit::rect(1.0, 0.5, 2, 2);
        let mesh = Mesh::from_solver_mesh(&sm).unwrap();

        let cap = CapacitanceMatrixData {
            terminals: vec!["pad1".into(), "pad2".into()],
            maxwell_matrix_f: vec![vec![80e-15, -30e-15], vec![-30e-15, 85e-15]],
            mutual_matrix_f: vec![vec![50e-15, 30e-15], vec![30e-15, 55e-15]],
            ec_matrix_hz: vec![vec![2.8e8, 1.0e8], vec![1.0e8, 2.6e8]],
        };
        let ind = InductanceMatrixData {
            terminals: vec!["loop1".into()],
            inductance_matrix_h: vec![vec![150e-12]],
            formulation: "mixed_3d".into(),
        };
        let res = SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            scalars: vec![cap.scalar(), ind.scalar()],
            ..Default::default()
        };

        assert!(CapacitanceMatrixAnalysis.applies_to(&res));
        assert!(InductanceMatrixAnalysis.applies_to(&res));

        let cap_json = CapacitanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(cap_json["terminals"], json!(["pad1", "pad2"]));
        assert_eq!(cap_json["maxwell_matrix_F"][0][0].as_f64().unwrap(), 80e-15);
        assert_eq!(cap_json["mutual_matrix_F"][0][1].as_f64().unwrap(), 30e-15);
        assert_eq!(cap_json["ec_matrix_Hz"][1][1].as_f64().unwrap(), 2.6e8);

        let ind_json = InductanceMatrixAnalysis.run(&mesh, &res).unwrap();
        assert_eq!(ind_json["terminals"], json!(["loop1"]));
        assert_eq!(ind_json["inductance_matrix_H"][0][0].as_f64().unwrap(), 150e-12);
        assert_eq!(ind_json["formulation"], json!("mixed_3d"));

        let empty = SolveResult::default();
        assert!(!CapacitanceMatrixAnalysis.applies_to(&empty));
        assert!(!InductanceMatrixAnalysis.applies_to(&empty));
        assert!(CapacitanceMatrixAnalysis.run(&mesh, &empty).is_err());
        assert!(InductanceMatrixAnalysis.run(&mesh, &empty).is_err());
    }
}
