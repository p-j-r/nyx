// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Analyses of results as engine commands: a result is solved once, an analysis is a function of it and of the mesh it was solved on
//! (`nyx-analysis`, ADR 56). The window draws the typed `ModeAnalysis`; the MCP server gets the same through `analyse`, as plain JSON.

use crate::engine::Engine;
use nyx_analysis::{
    optomechanics_with_meshes, tls_participation_quasi_electrostatic, uq_cqed_with_companion_and_params, Analysis,
    BeamOverlapAnalysis, CapacitanceMatrixAnalysis, EprAnalysis, FarFieldAnalysis, HamiltonianAnalysis,
    InductanceMatrixAnalysis, ModeAnalysis, NetworkAnalysis, OptomechanicsAnalysis, OptomechanicsDomainMat,
    OptomechanicsOptions, PiezoResonatorAnalysis, PolarizationAnalysis, TlsParticipationAnalysis, UqCqedAnalysis,
};
use nyx_physics::{CapacitanceMatrixData, CircuitRom, Ctx, Error};
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
        Box::new(FarFieldAnalysis),
        Box::new(BeamOverlapAnalysis),
        Box::new(OptomechanicsAnalysis),
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
        if analysis == "optomechanics"
            && params.get("mechanical_result").is_none_or(Value::is_null)
            && params.get("optical_result").is_none_or(Value::is_null)
        {
            return self.run_optomechanics_analysis(physics, &model, &result, params);
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

    /// Evaluates a zero-re-solve analytical post-solve EME sweep (`length`, `periodicity`, or `modes`)
    /// on the cached [`nyx_ems::EmeData`] of `physics` (defaulting to `"ems"`).
    pub fn eme_sweep(&self, args: &Value) -> Result<Value, Error> {
        let physics = args.get("physics").and_then(Value::as_str).unwrap_or("ems");
        let mode = args
            .get("mode")
            .and_then(Value::as_str)
            .ok_or_else(|| Error::invalid("`mode` is required ('length', 'periodicity', or 'modes')"))?;
        let (_model, result) = self.analysis_inputs(physics)?;
        if result.study != "eme" {
            return Err(Error::invalid(format!(
                "eme_sweep requires a solved 'eme' result on '{physics}', found study '{}'",
                result.study
            )));
        }
        let data = nyx_ems::EmeData::of(&result)
            .ok_or_else(|| Error::invalid(format!("result of '{physics}' does not carry cached EME cascade data")))?;

        let report = match mode {
            "length" => {
                let seg = match args.get("segment").filter(|v| !v.is_null()) {
                    None => 0usize,
                    Some(v) => {
                        v.as_u64().ok_or_else(|| Error::invalid("`segment` must be a non-negative integer"))? as usize
                    }
                };
                let lengths_val =
                    args.get("lengths").ok_or_else(|| Error::invalid("`lengths` is required when mode = 'length'"))?;
                let lengths = parse_f64_vec(lengths_val, "lengths")?;
                nyx_ems::eme_length_sweep(&data, seg, &lengths).map_err(Error::invalid)?
            }
            "periodicity" => {
                let subgrid = match args.get("subgrid").filter(|v| !v.is_null()) {
                    None => 0usize,
                    Some(v) => {
                        v.as_u64().ok_or_else(|| Error::invalid("`subgrid` must be a non-negative integer"))? as usize
                    }
                };
                let reps_arr = args.get("reps").and_then(Value::as_array).ok_or_else(|| {
                    Error::invalid("`reps` must be a non-empty list of positive integers when mode = 'periodicity'")
                })?;
                let mut reps = Vec::with_capacity(reps_arr.len());
                for (i, rv) in reps_arr.iter().enumerate() {
                    let r = rv
                        .as_u64()
                        .filter(|&n| n >= 1)
                        .ok_or_else(|| Error::invalid(format!("parameter 'reps[{i}]' must be an integer >= 1")))?
                        as usize;
                    reps.push(r);
                }
                nyx_ems::eme_periodicity_sweep(&data, subgrid, &reps).map_err(Error::invalid)?
            }
            "modes" => {
                let mc_arr = args.get("mode_counts").and_then(Value::as_array).ok_or_else(|| {
                    Error::invalid("`mode_counts` must be a non-empty list of positive integers when mode = 'modes'")
                })?;
                let mut mode_counts = Vec::with_capacity(mc_arr.len());
                for (i, mv) in mc_arr.iter().enumerate() {
                    let m = mv.as_u64().filter(|&n| n >= 1).ok_or_else(|| {
                        Error::invalid(format!("parameter 'mode_counts[{i}]' must be an integer >= 1"))
                    })? as usize;
                    mode_counts.push(m);
                }
                nyx_ems::eme_mode_sweep(&data, &mode_counts).map_err(Error::invalid)?
            }
            other => {
                return Err(Error::invalid(format!(
                    "unknown eme_sweep mode '{other}' (expected 'length', 'periodicity', or 'modes')"
                )));
            }
        };

        let mut v = serde_json::to_value(&report).map_err(|e| Error::invalid(e.to_string()))?;
        if let Some(o) = v.as_object_mut() {
            o.insert("physics".into(), json!(physics));
            o.insert("stale".into(), json!(self.is_stale(physics)));
            o.insert("age_s".into(), json!(self.result_age(physics)));
            o.insert("validated".into(), self.validated_json(physics));
        }
        Ok(v)
    }

    fn run_optomechanics_analysis(
        &self,
        physics: &str,
        model: &std::sync::Arc<crate::meshing::MeshedModel>,
        result: &std::sync::Arc<nyx_physics::SolveResult>,
        params: &Value,
    ) -> Result<Value, Error> {
        let (opt_phys, opt_model, opt_result, mech_phys, mech_model, mech_result) = if result.physics == "sms" {
            let op = params
                .get("optical_physics")
                .and_then(Value::as_str)
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .unwrap_or("ems");
            let (om, or) = self.analysis_inputs(op)?;
            (op.to_string(), om, or, physics.to_string(), model.clone(), result.clone())
        } else {
            let mp = params
                .get("mechanical_physics")
                .and_then(Value::as_str)
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .unwrap_or("sms");
            let (mm, mr) = self.analysis_inputs(mp)?;
            (physics.to_string(), model.clone(), result.clone(), mp.to_string(), mm, mr)
        };

        let mut opts = OptomechanicsOptions::from_params(params)?;
        for d in &opt_model.labels().domains {
            opts.opt_domain_names.entry(d.id).or_insert_with(|| d.name.clone());
        }
        for b in &opt_model.labels().boundaries {
            opts.opt_boundary_names.entry(b.id).or_insert_with(|| b.name.clone());
        }
        for d in &mech_model.labels().domains {
            opts.mech_domain_names.entry(d.id).or_insert_with(|| d.name.clone());
        }
        for b in &mech_model.labels().boundaries {
            opts.mech_boundary_names.entry(b.id).or_insert_with(|| b.name.clone());
        }

        let explicit_mat_keys: std::collections::BTreeSet<String> = params
            .get("domain_materials")
            .and_then(Value::as_object)
            .map(|o| o.keys().cloned().collect())
            .unwrap_or_default();

        let lib = self.library().ok();
        if let Some(opt_spec) = self.project.spec(&opt_phys) {
            if opts.optical_freq_hz.is_none() {
                if let Some(f) = opt_spec
                    .study
                    .params
                    .get("freq")
                    .or_else(|| opt_spec.study.params.get("f_search"))
                    .or_else(|| opt_spec.params.get("freq"))
                    .or_else(|| opt_spec.params.get("f_search"))
                    .and_then(Value::as_f64)
                    .filter(|f| *f > 0.0 && f.is_finite())
                {
                    opts.optical_freq_hz = Some(f);
                }
            }
            let wl = opts.optical_freq_hz.filter(|f| *f > 0.0).map(|f| femcore::optical::C0 / f);
            let mut ctx = Ctx::new(opt_model.as_ref());
            ctx.library = lib.clone();
            ctx.wavelength = wl;
            let props = self.registry.get(&opt_phys).map(|p| p.properties()).unwrap_or_default();
            if let Ok(mats) = nyx_physics::decode::materials_by_label(opt_spec, &ctx, &props) {
                for d in &opt_model.labels().domains {
                    if explicit_mat_keys.contains(&d.name) {
                        continue;
                    }
                    let eps_r = if let Some(an) = mats.aniso.get(&d.id) {
                        if let Some(ref p) = an.permittivity {
                            femcore::material::tensor_global(p, an.euler_deg)
                        } else {
                            femcore::material::iso_tensor(mats.get(d.id, "eps_r").unwrap_or(1.0))
                        }
                    } else {
                        femcore::material::iso_tensor(mats.get(d.id, "eps_r").unwrap_or(1.0))
                    };
                    let pe = mats.aniso.get(&d.id).and_then(|an| {
                        an.photoelastic
                            .as_ref()
                            .map(|p| femcore::material::coupled::photoelastic_global(p, an.euler_deg))
                    });
                    let entry = opts.domain_materials.entry(d.name.clone()).or_insert_with(|| OptomechanicsDomainMat {
                        eps_r,
                        photoelastic: pe,
                        rho: 0.0,
                    });
                    entry.eps_r = eps_r;
                    if entry.photoelastic.is_none() {
                        entry.photoelastic = pe;
                    }
                }
            }
        }

        if let Some(mech_spec) = self.project.spec(&mech_phys) {
            let mut ctx = Ctx::new(mech_model.as_ref());
            ctx.library = lib.clone();
            let props = self.registry.get(&mech_phys).map(|p| p.properties()).unwrap_or_default();
            if let Ok(mats) = nyx_physics::decode::materials_by_label(mech_spec, &ctx, &props) {
                for d in &mech_model.labels().domains {
                    if explicit_mat_keys.contains(&d.name) {
                        continue;
                    }
                    let rho = mats.get(d.id, "rho").unwrap_or(0.0);
                    let pe = mats.aniso.get(&d.id).and_then(|an| {
                        an.photoelastic
                            .as_ref()
                            .map(|p| femcore::material::coupled::photoelastic_global(p, an.euler_deg))
                    });
                    let entry = opts.domain_materials.entry(d.name.clone()).or_insert_with(|| OptomechanicsDomainMat {
                        eps_r: femcore::material::iso_tensor(1.0),
                        photoelastic: pe,
                        rho,
                    });
                    if rho > 0.0 {
                        entry.rho = rho;
                    }
                    if entry.photoelastic.is_none() {
                        entry.photoelastic = pe;
                    }
                }
            }
        }

        let rep = optomechanics_with_meshes(opt_model.mesh(), &opt_result, mech_model.mesh(), &mech_result, &opts)?;
        let mut v = rep.to_json();
        let combined_stale = self.is_stale(&opt_phys) || self.is_stale(&mech_phys);
        let validated = if opt_phys != mech_phys {
            self.validated_of(&[opt_phys.clone(), mech_phys.clone()])
        } else {
            self.validated_json(physics)
        };
        if let Some(o) = v.as_object_mut() {
            o.insert("analysis".into(), json!("optomechanics"));
            o.insert("physics".into(), json!(physics));
            o.insert("optical_physics".into(), json!(opt_phys));
            o.insert("mechanical_physics".into(), json!(mech_phys));
            o.insert("stale".into(), json!(combined_stale));
            o.insert("age_s".into(), json!(self.result_age(physics)));
            o.insert("validated".into(), validated);
        }
        Ok(v)
    }
}
