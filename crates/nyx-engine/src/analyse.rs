// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analyses of results as engine commands: a result is solved once, an analysis is a function of it and of the mesh it was solved on
//! (`nyx-analysis`, ADR 56). The window draws the typed `ModeAnalysis`; the MCP server gets the same through `analyse`, as plain JSON.

use crate::engine::Engine;
use nyx_analysis::{
    tls_participation_quasi_electrostatic, uq_cqed_with_companion_and_params, Analysis, CapacitanceMatrixAnalysis,
    EprAnalysis, HamiltonianAnalysis, InductanceMatrixAnalysis, ModeAnalysis, NetworkAnalysis, PiezoResonatorAnalysis,
    PolarizationAnalysis, TlsParticipationAnalysis, UqCqedAnalysis,
};
use nyx_physics::{CapacitanceMatrixData, CircuitRom, Error};
use serde_json::{json, Value};

/// The analyses of this build.
pub fn analyses() -> Vec<Box<dyn Analysis>> {
    vec![
        Box::new(PolarizationAnalysis),
        Box::new(CapacitanceMatrixAnalysis),
        Box::new(InductanceMatrixAnalysis),
        Box::new(EprAnalysis),
        Box::new(HamiltonianAnalysis),
        Box::new(TlsParticipationAnalysis),
        Box::new(NetworkAnalysis),
        Box::new(PiezoResonatorAnalysis),
        Box::new(UqCqedAnalysis),
    ]
}

fn parse_f64_vec(v: &Value, name: &str) -> Result<Vec<f64>, Error> {
    let arr = v.as_array().ok_or_else(|| Error::invalid(format!("parameter '{name}' must be a list of numbers")))?;
    let mut out = Vec::with_capacity(arr.len());
    for (i, x) in arr.iter().enumerate() {
        let f = x
            .as_f64()
            .filter(|n| n.is_finite())
            .ok_or_else(|| Error::invalid(format!("parameter '{name}[{i}]' must be a finite number")))?;
        out.push(f);
    }
    Ok(out)
}

fn floating_pair_weights(cap_opt: Option<&CapacitanceMatrixData>, i_minus: usize, i_plus: usize) -> (f64, f64) {
    if let Some(cap) = cap_opt {
        let m = &cap.maxwell_matrix_f;
        if i_minus < m.len() && i_plus < m.len() && i_plus < m[i_minus].len() && i_minus < m[i_plus].len() {
            let s_m = (m[i_minus][i_minus] + m[i_minus][i_plus]).max(0.0);
            let s_p = (m[i_plus][i_plus] + m[i_plus][i_minus]).max(0.0);
            if s_m + s_p > 1e-30 {
                return (-s_p / (s_m + s_p), s_m / (s_m + s_p));
            }
        }
    }
    (-0.5, 0.5)
}

impl Engine {
    /// The mesh a result was solved on, and the result itself (`Stored.model`, ADR 90).
    pub(crate) fn analysis_inputs(
        &self,
        physics: &str,
    ) -> Result<(std::sync::Arc<crate::meshing::MeshedModel>, std::sync::Arc<nyx_physics::SolveResult>), Error> {
        let stored = self
            .results
            .get(physics)
            .ok_or_else(|| Error::invalid(format!("there is no result for '{physics}' yet; solve first")))?;
        Ok((stored.model.clone(), stored.result.clone()))
    }

    fn run_hybrid_quasi_electrostatic_tls(
        &self,
        physics: &str,
        model: &std::sync::Arc<crate::meshing::MeshedModel>,
        result: &std::sync::Arc<nyx_physics::SolveResult>,
        params: &Value,
    ) -> Result<Value, Error> {
        let ep_opt =
            params.get("electrostatic_physics").or_else(|| params.get("electro_physics")).filter(|v| !v.is_null());
        let (electro_phys, electro_model, electro_result) = if let Some(ep_val) = ep_opt {
            let ep =
                ep_val.as_str().map(str::trim).filter(|s| !s.is_empty()).ok_or_else(|| {
                    Error::invalid("parameter 'electrostatic_physics' must be a non-empty physics code")
                })?;
            let (em, er) = self.analysis_inputs(ep)?;
            if er.study != "electrostatic" {
                return Err(Error::invalid(format!(
                    "electrostatic_physics '{ep}' holds a '{}' result, expected an 'electrostatic' study of dfs",
                    er.study
                )));
            }
            (ep.to_string(), em, er)
        } else {
            if result.study != "electrostatic" {
                return Err(Error::invalid(format!(
                    "modal_voltages without `electrostatic_physics` requires an 'electrostatic' result (got '{}')",
                    result.study
                )));
            }
            (physics.to_string(), model.clone(), result.clone())
        };

        let n_sol = electro_result.solutions.len();
        if n_sol == 0 {
            return Err(Error::invalid(format!(
                "electrostatic result '{electro_phys}' has no saved electric-field solutions (solve with save_fields: true)"
            )));
        }

        let explicit_mv = match params.get("modal_voltages").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => {
                let rows = v
                    .as_array()
                    .ok_or_else(|| Error::invalid("parameter 'modal_voltages' must be a list of voltage vectors"))?;
                if rows.is_empty() {
                    return Err(Error::invalid("parameter 'modal_voltages' must not be empty"));
                }
                let mut mv = Vec::with_capacity(rows.len());
                for (m, r) in rows.iter().enumerate() {
                    mv.push(parse_f64_vec(r, &format!("modal_voltages[{m}]"))?);
                }
                Some(mv)
            }
        };
        let explicit_mf = match params.get("modal_freqs_hz").filter(|v| !v.is_null()) {
            None => None,
            Some(v) => Some(parse_f64_vec(v, "modal_freqs_hz")?),
        };

        let (modal_voltages, modal_freqs_hz) = if let Some(mv) = explicit_mv {
            let mf = if let Some(mf) = explicit_mf {
                mf
            } else if let Some(rom) = CircuitRom::of(result) {
                if rom.modes.len() < mv.len() {
                    return Err(Error::invalid(format!(
                        "modal_voltages has {} mode(s), but CircuitRom only has {} mode(s); pass `modal_freqs_hz` explicitly",
                        mv.len(),
                        rom.modes.len()
                    )));
                }
                rom.modes.iter().take(mv.len()).map(|m| m.freq_hz).collect()
            } else if result.solutions.len() >= mv.len()
                && result.solutions.iter().take(mv.len()).all(|s| s.eigenvalue.is_some())
            {
                result.solutions.iter().take(mv.len()).map(|s| s.eigenvalue.unwrap()[0].abs()).collect()
            } else {
                return Err(Error::invalid(
                    "pass `modal_freqs_hz` alongside `modal_voltages` when the target result has no eigenmode frequencies",
                ));
            };
            (mv, mf)
        } else {
            let rom = CircuitRom::of(result).ok_or_else(|| {
                Error::invalid(
                    "hybrid quasi-electrostatic tls_participation requires CircuitRom on the target physics (add lumped_jj or lumped_port on ems cavity) or explicit `modal_voltages`",
                )
            })?;
            if rom.modes.is_empty() {
                return Err(Error::invalid("CircuitRom has no converged modes"));
            }
            let cap_opt = CapacitanceMatrixData::of(&electro_result);
            let find_term = |t_val: &Value| -> Result<usize, Error> {
                if let Some(idx) = t_val.as_u64() {
                    let i = idx as usize;
                    if i < n_sol {
                        return Ok(i);
                    }
                    return Err(Error::invalid(format!(
                        "terminal index {i} is out of bounds for {n_sol} electrostatic solution(s)"
                    )));
                }
                if let Some(name) = t_val.as_str() {
                    if let Some(ref cap) = cap_opt {
                        if let Some(pos) = cap.terminals.iter().position(|t| t == name) {
                            if pos < n_sol {
                                return Ok(pos);
                            }
                        }
                        return Err(Error::invalid(format!(
                            "terminal '{name}' not found in electrostatic terminals {:?}",
                            cap.terminals
                        )));
                    }
                }
                Err(Error::invalid("port_terminals entry must be a terminal name string or 0-based solution index"))
            };

            let n_ports = rom.ports.len();
            let mut port_map: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_ports];
            if let Some(pt_val) = params.get("port_terminals").filter(|v| !v.is_null()) {
                let pt_obj = pt_val.as_object().ok_or_else(|| {
                    Error::invalid(
                        "parameter 'port_terminals' must be an object mapping port names to a terminal name or [minus_pad, plus_pad]",
                    )
                })?;
                for (p_idx, port) in rom.ports.iter().enumerate() {
                    if let Some(spec) = pt_obj.get(&port.name).filter(|v| !v.is_null()) {
                        if let Some(arr) = spec.as_array() {
                            if arr.len() != 2 {
                                return Err(Error::invalid(format!(
                                    "port_terminals['{}'] array must have 2 entries [minus_pad, plus_pad]",
                                    port.name
                                )));
                            }
                            let im = find_term(&arr[0])?;
                            let ip = find_term(&arr[1])?;
                            let (wm, wp) = floating_pair_weights(cap_opt.as_ref(), im, ip);
                            port_map[p_idx] = vec![(im, wm), (ip, wp)];
                        } else {
                            port_map[p_idx] = vec![(find_term(spec)?, 1.0)];
                        }
                    }
                }
            } else {
                let mut named_matches = 0usize;
                if let Some(ref cap) = cap_opt {
                    for (p_idx, port) in rom.ports.iter().enumerate() {
                        if let Some(pos) = cap.terminals.iter().position(|t| t == &port.name) {
                            if pos < n_sol {
                                port_map[p_idx] = vec![(pos, 1.0)];
                                named_matches += 1;
                            }
                        }
                    }
                }
                if named_matches == 0 {
                    let jj_indices: Vec<usize> = rom
                        .ports
                        .iter()
                        .enumerate()
                        .filter_map(|(i, p)| (p.kind == "lumped_jj").then_some(i))
                        .collect();
                    if n_ports == n_sol {
                        for (p_idx, pm) in port_map.iter_mut().enumerate().take(n_ports) {
                            *pm = vec![(p_idx, 1.0)];
                        }
                    } else if n_ports == 1 && n_sol == 2 {
                        let (w0, w1) = floating_pair_weights(cap_opt.as_ref(), 0, 1);
                        port_map[0] = vec![(0, w0), (1, w1)];
                    } else if !jj_indices.is_empty() && jj_indices.len() == n_sol {
                        for (sol_i, &p_idx) in jj_indices.iter().enumerate() {
                            port_map[p_idx] = vec![(sol_i, 1.0)];
                        }
                    } else if jj_indices.len() == 1 && n_sol == 2 {
                        let (w0, w1) = floating_pair_weights(cap_opt.as_ref(), 0, 1);
                        port_map[jj_indices[0]] = vec![(0, w0), (1, w1)];
                    } else {
                        return Err(Error::invalid(format!(
                            "cannot auto-map {} CircuitRom port(s) to {} electrostatic terminal solution(s); specify `port_terminals` or `modal_voltages` in params",
                            n_ports, n_sol
                        )));
                    }
                }
            }

            let mut mv = Vec::new();
            let mut mf = Vec::new();
            for (m_idx, mode_row) in rom.modes.iter().enumerate() {
                let mut best_norm2 = 0.0_f64;
                let mut best_phasor = [1.0_f64, 0.0_f64];
                for (p_idx, &[vr, vi]) in mode_row.v_mj.iter().enumerate() {
                    if p_idx < port_map.len() && !port_map[p_idx].is_empty() {
                        let n2 = vr * vr + vi * vi;
                        if n2 > best_norm2 {
                            best_norm2 = n2;
                            best_phasor = [vr, vi];
                        }
                    }
                }
                let b_norm = best_phasor[0].hypot(best_phasor[1]);
                let mut w_vec = vec![0.0_f64; n_sol];
                if b_norm > 1e-30 {
                    let (br, bi) = (best_phasor[0] / b_norm, best_phasor[1] / b_norm);
                    for (p_idx, &[vr, vi]) in mode_row.v_mj.iter().enumerate() {
                        if p_idx < port_map.len() {
                            let v_signed = vr * br + vi * bi;
                            for &(sol_i, coeff) in &port_map[p_idx] {
                                w_vec[sol_i] += coeff * v_signed;
                            }
                        }
                    }
                }
                if w_vec.iter().any(|&x| x.abs() > 1e-30) {
                    let freq = explicit_mf.as_ref().and_then(|fs| fs.get(m_idx).copied()).unwrap_or(mode_row.freq_hz);
                    mv.push(w_vec);
                    mf.push(freq);
                }
            }
            if mv.is_empty() {
                return Err(Error::invalid(
                    "all CircuitRom modes have zero terminal voltage on the mapped electrostatic terminals",
                ));
            }
            (mv, mf)
        };

        let modes = tls_participation_quasi_electrostatic(
            electro_model.mesh(),
            &electro_result,
            &modal_voltages,
            &modal_freqs_hz,
            params,
        )?;
        let sel = match params.get("solution").filter(|v| !v.is_null()) {
            None => 0usize,
            Some(v) => {
                let s =
                    v.as_u64().ok_or_else(|| Error::invalid("parameter 'solution' must be a non-negative integer"))?
                        as usize;
                if s >= modes.len() {
                    return Err(Error::invalid(format!(
                        "solution index {s} is out of bounds (hybrid quasi-electrostatic TLS evaluated {} mode(s))",
                        modes.len()
                    )));
                }
                s
            }
        };

        let mut v = modes[sel].to_json();
        let combined_stale = self.is_stale(physics) || (electro_phys != physics && self.is_stale(&electro_phys));
        let validated = if electro_phys != physics {
            self.validated_of(&[physics.to_string(), electro_phys.clone()])
        } else {
            self.validated_json(physics)
        };
        if let Some(o) = v.as_object_mut() {
            o.insert("quasi_electrostatic".into(), json!(true));
            o.insert("electrostatic_physics".into(), json!(electro_phys));
            o.insert("modes".into(), json!(modes.iter().map(|m| m.to_json()).collect::<Vec<_>>()));
            o.insert("analysis".into(), json!("tls_participation"));
            o.insert("physics".into(), json!(physics));
            o.insert("stale".into(), json!(combined_stale));
            o.insert("age_s".into(), json!(self.result_age(physics)));
            o.insert("validated".into(), validated);
        }
        Ok(v)
    }

    /// The polarization analysis of the modes of `physics`, typed (what the Analysis tab draws).
    pub fn polarization(&self, physics: &str) -> Result<ModeAnalysis, Error> {
        let (model, result) = self.analysis_inputs(physics)?;
        nyx_analysis::polarization(model.mesh(), &result)
    }

    /// Whether the polarization applies to the latest result of `physics`.
    pub fn polarization_applies(&self, physics: &str) -> bool {
        self.results.get(physics).is_some_and(|s| nyx_analysis::polarization_applies(&s.result))
    }

    /// Lists the analyses and says which apply to which result.
    pub fn list_analyses(&self) -> Value {
        json!(analyses()
            .iter()
            .map(|a| json!({
                "analysis": a.id(),
                "name": a.name(),
                "applies_to": self.results.iter().filter(|(_, s)| a.applies_to(&s.result)).map(|(k, _)| k.clone()).collect::<Vec<_>>(),
                "stale": self.results.iter().filter(|(_, s)| a.applies_to(&s.result)).map(|(k, _)| (k.clone(), self.is_stale(k))).collect::<std::collections::BTreeMap<_, _>>(),
            }))
            .collect::<Vec<_>>())
    }

    /// Runs one analysis on the latest result of `physics` and returns its JSON.
    pub fn analyse(&self, analysis: &str, physics: &str) -> Result<Value, Error> {
        self.analyse_with(analysis, physics, &Value::Null)
    }

    /// Runs one analysis on the latest result of `physics` with optional `params` and returns its JSON.
    pub fn analyse_with(&self, analysis: &str, physics: &str, params: &Value) -> Result<Value, Error> {
        let all = analyses();
        let a = all.iter().find(|a| a.id() == analysis).ok_or_else(|| {
            Error::invalid(format!(
                "there is no analysis '{analysis}' (they are: {})",
                all.iter().map(|a| a.id()).collect::<Vec<_>>().join(", ")
            ))
        })?;
        let (model, result) = self.analysis_inputs(physics)?;
        if analysis == "tls_participation"
            && (params
                .get("electrostatic_physics")
                .or_else(|| params.get("electro_physics"))
                .is_some_and(|v| !v.is_null())
                || params.get("modal_voltages").is_some_and(|v| !v.is_null()))
        {
            return self.run_hybrid_quasi_electrostatic_tls(physics, &model, &result, params);
        }
        if analysis == "uq_cqed"
            && params
                .get("electrostatic_physics")
                .or_else(|| params.get("electro_physics"))
                .is_some_and(|v| !v.is_null())
        {
            let ep_val = params.get("electrostatic_physics").or_else(|| params.get("electro_physics")).unwrap();
            let ep =
                ep_val.as_str().map(str::trim).filter(|s| !s.is_empty()).ok_or_else(|| {
                    Error::invalid("parameter 'electrostatic_physics' must be a non-empty physics code")
                })?;
            let (em, er) = self.analysis_inputs(ep)?;
            let mut v = uq_cqed_with_companion_and_params(Some(model.mesh()), &result, Some((em.mesh(), &er)), params)?
                .to_json();
            let combined_stale = self.is_stale(physics) || (ep != physics && self.is_stale(ep));
            let validated = if ep != physics {
                self.validated_of(&[physics.to_string(), ep.to_string()])
            } else {
                self.validated_json(physics)
            };
            if let Some(o) = v.as_object_mut() {
                o.insert("analysis".into(), json!(analysis));
                o.insert("physics".into(), json!(physics));
                o.insert("electrostatic_physics".into(), json!(ep));
                o.insert("stale".into(), json!(combined_stale));
                o.insert("age_s".into(), json!(self.result_age(physics)));
                o.insert("validated".into(), validated);
            }
            return Ok(v);
        }
        if !a.applies_to(&result) {
            return Err(Error::invalid(format!(
                "the {analysis} analysis does not apply to the {} result of '{physics}'",
                result.study
            )));
        }
        let mut v = a.run_with_params(model.mesh(), &result, params)?;
        if let Some(o) = v.as_object_mut() {
            o.insert("analysis".into(), json!(analysis));
            o.insert("physics".into(), json!(physics));
            o.insert("stale".into(), json!(self.is_stale(physics)));
            o.insert("age_s".into(), json!(self.result_age(physics)));
            o.insert("validated".into(), self.validated_json(physics));
        }
        Ok(v)
    }
}
