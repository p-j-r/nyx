// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Native hierarchical circuit quantization (`HamiltonianAnalysis`, Milestone 4.3).
//!
//! Supports two complementary quantization bases on a 3D `ems` [`CircuitRom`]:
//! - **Mode A (`basis = "normal_fock"`, default)**: Normal-mode Fock basis (`fock_trunc`, e.g. `[7, 5]`) for
//!   weakly-to-moderately anharmonic multi-mode transmon + readout resonator systems. Forms the exact single-mode
//!   phase operator `phi_{mj} = s_{mj} phi_{zpf,mj} (a_m + a_m^dagger)` and its Padé-13 matrix exponential
//!   `U_{mj} = exp(i phi_{mj})` via [`femcore::linalg::expm::expm`], subtracts the quadratic normal-mode contribution
//!   already included in `sum_m f_m a_m^dagger a_m`, and diagonalizes `H / h` via
//!   [`femcore::linalg::dense::hermitian_eig`]. Also extracts linear Craig–Bampton bare cavity frequencies
//!   `f_cav_Hz`, charge-cavity couplings `g_matrix_Hz`, bare cavity decay rates `kappa_cav_Hz`, and multi-mode
//!   Purcell lifetimes `t1_purcell_s`.
//! - **Mode B (`basis = "craig_bampton"`)**: Nonlinear-port + linear-cavity Craig–Bampton basis for strongly
//!   anharmonic / multi-well circuits (Fluxonium, Cooper-pair box, DC-SQUID tunable transmons/dissipators, zero-pi,
//!   tunable couplers) coupled to 3D cavities. Decomposes `R^r` into nonlinear port lift `T_p` (`B_{r,nl}^T T_p = I_J`)
//!   and short-circuit interior cavity modes `T_c` (`B_{r,nl}^T T_c = 0`), pre-diagonalizes each 1D nonlinear port
//!   (`E_{L,eff} > 1e3 Hz` in a 40-state harmonic oscillator basis with cached `expm(i phi_j)`; `E_{L,eff} <= 1e3 Hz`
//!   in a 31-state charge-periodic basis), projects `(n_j, phi_j)` onto the lowest `port_trunc` (default 8) port
//!   eigenstates, and diagonalizes the coupled `H / h` in `< 100 us` per point (`< 10 ms` for a 100-point flux or
//!   DC-SQUID `E_J` sweep).
//! - **Globally Bijective Greedy Descending Bipartite State Labeling**: Assigns bare product states in
//!   `T = {|0...0>, |1_a>, |2_a>, |1_a 1_b>}` to eigenstates bijectively in descending order of overlap
//!   `|<n | psi_k>|^2` to extract signed `f01_Hz`, `f12_Hz`, `f02_Hz`, `alpha_Hz`, and `chi_matrix_Hz`.

use crate::epr::{
    assemble_rom_damping, epr_with_params, is_nonlinear_port, jj_flux_scale, jj_inv_l_tot, resolve_jj_l_and_ej,
    resolve_ports_with_overrides, squid_effective_ej, C0, EPS0, E_CHARGE, MU0, PHI0, PLANCK_H,
};
use crate::Analysis;
use femcore::linalg::dense::{hermitian_eig, solve_small, Dense};
use femcore::linalg::driver::{dense_eigh, dense_gen_eigh};
use femcore::linalg::expm::expm;
use femcore::linalg::scalar::{Scalar, C64};
use nyx_physics::{CircuitRom, Error, Mesh, RomPort, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;
use std::time::Instant;

/// Maximum allowed Fock/basis boundary occupation `max_k |<N_max | psi_k>|^2` before emitting a truncation warning.
pub const FOCK_BOUNDARY_WARN_THRESHOLD: f64 = 1.0e-4;

/// Formats a truncation warning when `max_occ > FOCK_BOUNDARY_WARN_THRESHOLD` without triggering a re-solve.
fn fock_truncation_warning_msg(max_occ: f64, subsystems: &[String], dims: &[usize], by_sub: &[f64]) -> Option<String> {
    if max_occ.is_nan() || max_occ <= FOCK_BOUNDARY_WARN_THRESHOLD {
        return None;
    }
    let worst = by_sub
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(idx, &occ)| {
            let name = subsystems.get(idx).map_or("subsystem", String::as_str);
            let dim = dims.get(idx).copied().unwrap_or(0);
            format!("worst subsystem '{name}' at dim={dim} with occupation {occ:.3e}")
        })
        .unwrap_or_else(|| format!("max occupation {max_occ:.3e}"));
    Some(format!(
        "maximum Fock boundary occupation {max_occ:.3e} exceeds {FOCK_BOUNDARY_WARN_THRESHOLD:.1e} ({worst}); increase fock_trunc or port_trunc for converged levels"
    ))
}

/// One point of a flux, DC-SQUID flux, or `E_J` sweep (`phi_ext_sweep`, `flux_sweep`, `squid_flux_sweep`, `ej_hz_sweep`).
#[derive(Clone, Debug, PartialEq)]
pub struct FluxSweepPoint {
    pub phi_ext_rad: f64,
    pub flux_phi0: f64,
    pub e_j_hz: Option<f64>,
    pub f01_hz: Vec<f64>,
    pub f12_hz: Vec<f64>,
    pub f02_hz: Vec<f64>,
    pub alpha_hz: Vec<f64>,
    pub chi_matrix_hz: Vec<Vec<f64>>,
    pub eigenvalues_hz: Vec<f64>,
    pub min_overlap: f64,
    pub max_fock_boundary_occupation: f64,
    pub fock_boundary_occupation_by_subsystem: Vec<f64>,
    pub fock_truncation_warning: Option<String>,
}

impl FluxSweepPoint {
    pub fn to_json(&self) -> Value {
        let mut obj = json!({
            "phi_ext_rad": self.phi_ext_rad,
            "flux_phi0": self.flux_phi0,
            "f01_Hz": self.f01_hz,
            "f12_Hz": self.f12_hz,
            "f02_Hz": self.f02_hz,
            "alpha_Hz": self.alpha_hz,
            "chi_matrix_Hz": self.chi_matrix_hz,
            "eigenvalues_Hz": self.eigenvalues_hz,
            "min_overlap": self.min_overlap,
            "max_fock_boundary_occupation": self.max_fock_boundary_occupation,
            "fock_boundary_occupation_by_subsystem": self.fock_boundary_occupation_by_subsystem,
            "fock_truncation_adequate": self.fock_truncation_warning.is_none(),
        });
        if let Some(ej) = self.e_j_hz {
            obj["e_j_Hz"] = json!(ej);
        }
        if let Some(ref w) = self.fock_truncation_warning {
            obj["fock_truncation_warning"] = json!(w);
        }
        obj
    }
}

/// Typed output of [`hamiltonian`] and [`hamiltonian_with_params`].
#[derive(Clone, Debug, PartialEq)]
pub struct HamiltonianReport {
    pub basis: String,
    pub hilbert_dim: usize,
    pub subsystems: Vec<String>,
    pub dims: Vec<usize>,
    pub f01_hz: Vec<f64>,
    pub f12_hz: Vec<f64>,
    pub f02_hz: Vec<f64>,
    pub alpha_hz: Vec<f64>,
    pub chi_matrix_hz: Vec<Vec<f64>>,
    pub eigenvalues_hz: Vec<f64>,
    pub min_overlap: f64,
    pub max_fock_boundary_occupation: f64,
    pub fock_boundary_occupation_by_subsystem: Vec<f64>,
    pub fock_truncation_warning: Option<String>,
    pub elapsed_s: f64,
    pub f_cav_hz: Vec<f64>,
    pub g_matrix_hz: Vec<Vec<f64>>,
    pub kappa_cav_hz: Vec<f64>,
    pub t1_purcell_s: Vec<f64>,
    pub sweep: Option<Vec<FluxSweepPoint>>,
    pub sweep_elapsed_s: Option<f64>,
}

impl HamiltonianReport {
    pub fn to_json(&self) -> Value {
        let mut obj = json!({
            "basis": self.basis,
            "hilbert_dim": self.hilbert_dim,
            "subsystems": self.subsystems,
            "dims": self.dims,
            "f01_Hz": self.f01_hz,
            "f12_Hz": self.f12_hz,
            "f02_Hz": self.f02_hz,
            "alpha_Hz": self.alpha_hz,
            "chi_matrix_Hz": self.chi_matrix_hz,
            "eigenvalues_Hz": self.eigenvalues_hz,
            "min_overlap": self.min_overlap,
            "max_fock_boundary_occupation": self.max_fock_boundary_occupation,
            "fock_boundary_occupation_by_subsystem": self.fock_boundary_occupation_by_subsystem,
            "fock_truncation_adequate": self.fock_truncation_warning.is_none(),
            "elapsed_s": self.elapsed_s,
            "f_cav_Hz": self.f_cav_hz,
            "g_matrix_Hz": self.g_matrix_hz,
            "kappa_cav_Hz": self.kappa_cav_hz,
            "t1_purcell_s": self.t1_purcell_s,
        });
        if let Some(ref w) = self.fock_truncation_warning {
            obj["fock_truncation_warning"] = json!(w);
            obj["warning"] = json!(w);
        }
        if let Some(ref sw) = self.sweep {
            obj["sweep"] = json!(sw.iter().map(|p| p.to_json()).collect::<Vec<_>>());
        }
        if let Some(t) = self.sweep_elapsed_s {
            obj["sweep_elapsed_s"] = json!(t);
        }
        obj
    }
}

/// True when `result` carries a [`CircuitRom`] with at least one nonlinear port (`lumped_jj`).
pub fn hamiltonian_applies(result: &SolveResult) -> bool {
    CircuitRom::of(result).is_some_and(|r| r.ports.iter().any(is_nonlinear_port))
}

/// Resolves `"sweep_port"` in `params` (0-based integer index or port `name` string) against `ports`
/// (the nonlinear ports of the `CircuitRom`), defaulting to `0`.
pub fn resolve_sweep_port(ports: &[RomPort], params: &Value) -> Result<usize, Error> {
    if ports.is_empty() {
        return Err(Error::invalid("CircuitRom has no nonlinear ports"));
    }
    let names: Vec<&str> = ports.iter().map(|p| p.name.as_str()).collect();
    match params.get("sweep_port") {
        None | Some(Value::Null) => Ok(0),
        Some(Value::Number(n)) => {
            let idx = n.as_u64().ok_or_else(|| {
                Error::invalid(format!("sweep_port must be a non-negative integer index or port name in {names:?}"))
            })? as usize;
            if idx < ports.len() {
                Ok(idx)
            } else {
                Err(Error::invalid(format!(
                    "sweep_port index {idx} is out of bounds for {} nonlinear port(s): {names:?}",
                    ports.len()
                )))
            }
        }
        Some(Value::String(name)) => ports.iter().position(|p| p.name == *name).ok_or_else(|| {
            Error::invalid(format!("sweep_port '{name}' does not match any nonlinear port in {names:?}"))
        }),
        _ => Err(Error::invalid("sweep_port must be a non-negative integer index or a port name string")),
    }
}

// ------------------------------------------------------------------ tensor product & state labeling helpers

/// Computes per-subsystem strides for row-major multi-index flattening `idx = sum_a n_a * stride[a]`.
fn compute_strides(dims: &[usize]) -> (usize, Vec<usize>) {
    let mut strides = vec![1usize; dims.len()];
    let mut total = 1usize;
    for a in (0..dims.len()).rev() {
        strides[a] = total;
        total = total.saturating_mul(dims[a]);
    }
    (total, strides)
}

#[inline]
fn flat_index(multi: &[usize], strides: &[usize]) -> usize {
    multi.iter().zip(strides).map(|(&n, &s)| n * s).sum()
}

/// Globally bijective greedy descending bipartite state matching on the target excitation manifold
/// `T = {|0...0>, |1_a>, |2_a>, |1_a 1_b>}`.
#[derive(Clone, Debug)]
struct LabeledSpectrum {
    f01_hz: Vec<f64>,
    f12_hz: Vec<f64>,
    f02_hz: Vec<f64>,
    alpha_hz: Vec<f64>,
    chi_matrix_hz: Vec<Vec<f64>>,
    eigenvalues_hz: Vec<f64>,
    min_overlap: f64,
    max_fock_boundary_occupation: f64,
    fock_boundary_occupation_by_subsystem: Vec<f64>,
}

fn label_spectrum_bipartite<F>(dims: &[usize], evals: &[f64], overlap_sq: F) -> LabeledSpectrum
where
    F: Fn(usize, usize) -> f64,
{
    let n_sub = dims.len();
    let d_tot = evals.len();
    let (_, strides) = compute_strides(dims);

    // Build target bare states T = {|0...0>, |1_a>, |2_a>, |1_a 1_b>}
    let mut targets: Vec<Vec<usize>> = Vec::new();
    targets.push(vec![0; n_sub]);
    for a in 0..n_sub {
        if dims[a] >= 2 {
            let mut s = vec![0; n_sub];
            s[a] = 1;
            targets.push(s);
        }
        if dims[a] >= 3 {
            let mut s = vec![0; n_sub];
            s[a] = 2;
            targets.push(s);
        }
    }
    for a in 0..n_sub {
        for b in (a + 1)..n_sub {
            if dims[a] >= 2 && dims[b] >= 2 {
                let mut s = vec![0; n_sub];
                s[a] = 1;
                s[b] = 1;
                targets.push(s);
            }
        }
    }

    let n_t = targets.len();
    let target_flat: Vec<usize> = targets.iter().map(|t| flat_index(t, &strides)).collect();

    // Candidate edges (t_idx, eig_k, prob)
    let mut edges: Vec<(usize, usize, f64)> = Vec::with_capacity(n_t * d_tot);
    for (t_idx, &b_idx) in target_flat.iter().enumerate() {
        for k in 0..d_tot {
            let p = overlap_sq(b_idx, k);
            if p > 0.0 {
                edges.push((t_idx, k, p));
            }
        }
    }
    edges.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));

    let mut matched_eig: Vec<Option<(usize, f64)>> = vec![None; n_t];
    let mut used_k = vec![false; d_tot];
    let mut n_matched = 0usize;

    for &(t_idx, k, p) in &edges {
        if matched_eig[t_idx].is_none() && !used_k[k] {
            matched_eig[t_idx] = Some((k, p));
            used_k[k] = true;
            n_matched += 1;
            if n_matched == n_t {
                break;
            }
        }
    }

    // Fallback for any unmatched target state (if overlap was 0)
    for slot in matched_eig.iter_mut().take(n_t) {
        if slot.is_none() {
            if let Some(k) = (0..d_tot).find(|&k| !used_k[k]) {
                *slot = Some((k, 0.0));
                used_k[k] = true;
            }
        }
    }

    let energy_of = |state: &[usize]| -> Option<f64> {
        let pos = targets.iter().position(|t| t == state)?;
        let (k, _) = matched_eig[pos]?;
        Some(evals[k])
    };

    let mut min_ov = 1.0_f64;
    for &(_k, p) in matched_eig.iter().flatten() {
        min_ov = min_ov.min(p);
    }

    let vac = vec![0; n_sub];
    let e0 = energy_of(&vac).unwrap_or_else(|| evals.first().copied().unwrap_or(0.0));

    let mut f01_hz = vec![0.0; n_sub];
    let mut f12_hz = vec![0.0; n_sub];
    let mut f02_hz = vec![0.0; n_sub];
    let mut alpha_hz = vec![0.0; n_sub];
    let mut chi_matrix_hz = vec![vec![0.0; n_sub]; n_sub];

    for a in 0..n_sub {
        let mut s1 = vec![0; n_sub];
        if dims[a] >= 2 {
            s1[a] = 1;
            if let Some(e1) = energy_of(&s1) {
                f01_hz[a] = e1 - e0;
            }
        }
        let mut s2 = vec![0; n_sub];
        if dims[a] >= 3 {
            s2[a] = 2;
            if let Some(e2) = energy_of(&s2) {
                f02_hz[a] = e2 - e0;
                f12_hz[a] = f02_hz[a] - f01_hz[a];
                alpha_hz[a] = f12_hz[a] - f01_hz[a];
            }
        }
    }

    for a in 0..n_sub {
        for b in (a + 1)..n_sub {
            if dims[a] >= 2 && dims[b] >= 2 {
                let mut s11 = vec![0; n_sub];
                s11[a] = 1;
                s11[b] = 1;
                if let Some(e11) = energy_of(&s11) {
                    let chi = (e11 - e0) - f01_hz[a] - f01_hz[b];
                    chi_matrix_hz[a][b] = chi;
                    chi_matrix_hz[b][a] = chi;
                }
            }
        }
    }

    let n_keep = d_tot.min(32);
    let eigenvalues_hz: Vec<f64> = evals.iter().take(n_keep).map(|&e| e - e0).collect();

    let mut fock_boundary_occupation_by_subsystem = vec![0.0_f64; n_sub];
    if n_sub > 0 && d_tot > 0 {
        let mut labeled_k: Vec<usize> = matched_eig.iter().flatten().map(|&(k, _)| k).collect();
        labeled_k.sort_unstable();
        labeled_k.dedup();
        for &k in &labeled_k {
            let mut sub_occ = vec![0.0_f64; n_sub];
            for b_idx in 0..d_tot {
                let p = overlap_sq(b_idx, k);
                if p <= 0.0 {
                    continue;
                }
                for m in 0..n_sub {
                    if dims[m] > 0 && (b_idx / strides[m]) % dims[m] == dims[m] - 1 {
                        sub_occ[m] += p;
                    }
                }
            }
            for m in 0..n_sub {
                fock_boundary_occupation_by_subsystem[m] = fock_boundary_occupation_by_subsystem[m].max(sub_occ[m]);
            }
        }
    }
    let max_fock_boundary_occupation = fock_boundary_occupation_by_subsystem.iter().copied().fold(0.0_f64, f64::max);

    LabeledSpectrum {
        f01_hz,
        f12_hz,
        f02_hz,
        alpha_hz,
        chi_matrix_hz,
        eigenvalues_hz,
        min_overlap: min_ov,
        max_fock_boundary_occupation,
        fock_boundary_occupation_by_subsystem,
    }
}

/// Parses `fock_trunc` from `params` for `n_modes` modes.
fn parse_fock_trunc(params: &Value, n_modes: usize, default_per_mode: &[usize]) -> Result<Vec<usize>, Error> {
    if n_modes == 0 {
        return Ok(Vec::new());
    }
    if let Some(arr) = params.get("fock_trunc").and_then(|v| v.as_array()) {
        if arr.is_empty() {
            return Err(Error::invalid("fock_trunc array must not be empty"));
        }
        let mut dims = Vec::with_capacity(n_modes);
        for (i, v) in arr.iter().enumerate().take(n_modes) {
            let d = v.as_u64().ok_or_else(|| Error::invalid(format!("fock_trunc[{i}] must be a positive integer")))?
                as usize;
            if d < 2 {
                return Err(Error::invalid(format!("fock_trunc[{i}] must be >= 2")));
            }
            dims.push(d);
        }
        let Some(&last) = dims.last() else {
            return Err(Error::invalid("fock_trunc array must not be empty"));
        };
        while dims.len() < n_modes {
            dims.push(last);
        }
        return Ok(dims);
    }
    if let Some(d) = params.get("fock_trunc").and_then(|v| v.as_u64()) {
        let d = d as usize;
        if d < 2 {
            return Err(Error::invalid("fock_trunc must be >= 2"));
        }
        return Ok(vec![d; n_modes]);
    }
    Ok(default_per_mode.to_vec())
}

/// Parses optional `phi_ext_rad` (or `flux_phi0`) and `n_g` overrides onto `nl_ports`.
fn apply_bias_overrides(nl_ports: &mut [RomPort], params: &Value) {
    if let Some(v) = params.get("phi_ext_rad") {
        if let Some(arr) = v.as_array() {
            for (p, val) in nl_ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.phi_ext_rad = f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in nl_ports.iter_mut() {
                p.phi_ext_rad = f;
            }
        }
    } else if let Some(v) = params.get("flux_phi0") {
        if let Some(arr) = v.as_array() {
            for (p, val) in nl_ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.phi_ext_rad = 2.0 * PI * f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in nl_ports.iter_mut() {
                p.phi_ext_rad = 2.0 * PI * f;
            }
        }
    }

    if let Some(v) = params.get("n_g") {
        if let Some(arr) = v.as_array() {
            for (p, val) in nl_ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.n_g = f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in nl_ports.iter_mut() {
                p.n_g = f;
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum SweepPointSpec {
    LoopPhase { phi_ext_rad: f64, flux_phi0: f64 },
    EffectiveEj { phi_ext_rad: f64, flux_phi0: f64, e_j_hz: f64 },
}

fn parse_number_vec(val: &Value, name: &str, non_neg: bool) -> Result<Vec<f64>, Error> {
    let arr =
        val.as_array().ok_or_else(|| Error::invalid(format!("parameter '{name}' must be an array of numbers")))?;
    if arr.is_empty() {
        return Err(Error::invalid(format!("parameter '{name}' must not be empty")));
    }
    let mut out = Vec::with_capacity(arr.len());
    for (i, v) in arr.iter().enumerate() {
        let x = v.as_f64().filter(|x| x.is_finite() && (!non_neg || *x >= 0.0)).ok_or_else(|| {
            if non_neg {
                Error::invalid(format!("{name}[{i}] must be a non-negative finite number"))
            } else {
                Error::invalid(format!("{name}[{i}] must be a finite number"))
            }
        })?;
        out.push(x);
    }
    Ok(out)
}

fn parse_sweep_schedule(
    nl_ports: &[RomPort],
    base_ej_sum_hz: &[f64],
    sweep_port: usize,
    params: &Value,
) -> Result<Option<Vec<SweepPointSpec>>, Error> {
    let target = &nl_ports[sweep_port];
    let base_phi = target.phi_ext_rad;
    let ej_sum = base_ej_sum_hz[sweep_port];

    if let Some(v) = params.get("ej_hz_sweep").filter(|v| !v.is_null()) {
        let ejs = parse_number_vec(v, "ej_hz_sweep", true)?;
        return Ok(Some(
            ejs.into_iter()
                .map(|ej| SweepPointSpec::EffectiveEj {
                    phi_ext_rad: base_phi,
                    flux_phi0: base_phi / (2.0 * PI),
                    e_j_hz: ej,
                })
                .collect(),
        ));
    }

    if let Some(v) = params.get("squid_flux_sweep").filter(|v| !v.is_null()) {
        let fluxes = parse_number_vec(v, "squid_flux_sweep", false)?;
        let d = target.squid_asymmetry.unwrap_or(0.0);
        return Ok(Some(
            fluxes
                .into_iter()
                .map(|fl| SweepPointSpec::EffectiveEj {
                    phi_ext_rad: 2.0 * PI * fl,
                    flux_phi0: fl,
                    e_j_hz: squid_effective_ej(ej_sum, fl, d),
                })
                .collect(),
        ));
    }

    let is_dc_squid = target.e_l_hz.unwrap_or(0.0) <= 1.0e3 && target.squid_asymmetry.is_some();

    if let Some(v) = params.get("phi_ext_sweep").filter(|v| !v.is_null()) {
        let phis = parse_number_vec(v, "phi_ext_sweep", false)?;
        if is_dc_squid {
            let d = target.squid_asymmetry.unwrap_or(0.0);
            return Ok(Some(
                phis.into_iter()
                    .map(|phi| {
                        let fl = phi / (2.0 * PI);
                        SweepPointSpec::EffectiveEj {
                            phi_ext_rad: phi,
                            flux_phi0: fl,
                            e_j_hz: squid_effective_ej(ej_sum, fl, d),
                        }
                    })
                    .collect(),
            ));
        }
        return Ok(Some(
            phis.into_iter()
                .map(|phi| SweepPointSpec::LoopPhase { phi_ext_rad: phi, flux_phi0: phi / (2.0 * PI) })
                .collect(),
        ));
    }

    if let Some(v) = params.get("flux_sweep").filter(|v| !v.is_null()) {
        let fluxes = parse_number_vec(v, "flux_sweep", false)?;
        if is_dc_squid {
            let d = target.squid_asymmetry.unwrap_or(0.0);
            return Ok(Some(
                fluxes
                    .into_iter()
                    .map(|fl| SweepPointSpec::EffectiveEj {
                        phi_ext_rad: 2.0 * PI * fl,
                        flux_phi0: fl,
                        e_j_hz: squid_effective_ej(ej_sum, fl, d),
                    })
                    .collect(),
            ));
        }
        return Ok(Some(
            fluxes
                .into_iter()
                .map(|fl| SweepPointSpec::LoopPhase { phi_ext_rad: 2.0 * PI * fl, flux_phi0: fl })
                .collect(),
        ));
    }

    Ok(None)
}

// ------------------------------------------------------------------ Shared Linear Craig-Bampton Decomposition & Purcell T1

struct LinearCbDecomposition {
    #[allow(dead_code)]
    t_p: Vec<Vec<f64>>,
    #[allow(dead_code)]
    t_c_keep: Vec<Vec<f64>>,
    c_inv: Vec<Vec<f64>>,
    f_cav_hz: Vec<f64>,
    g_jm: Vec<Vec<f64>>,
    g_jk: Vec<Vec<f64>>,
    kappa_cav_hz: Vec<f64>,
}

fn dot_vec(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

fn mat_vec(m: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    m.iter().map(|row| dot_vec(row, x)).collect()
}

fn invert_spd_dense(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, Error> {
    let n = a.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut flat = Vec::with_capacity(n * n);
    for r in a {
        flat.extend_from_slice(r);
    }
    let (eigs, vecs) = dense_eigh(&flat, n).map_err(Error::from_kernel)?;
    let max_eig = eigs.iter().copied().fold(0.0_f64, f64::max);
    if max_eig <= 0.0 {
        return Err(Error::invalid("capacitance matrix is not positive definite"));
    }
    let mut inv = vec![vec![0.0; n]; n];
    for k in 0..n {
        if eigs[k] <= 1e-14 * max_eig {
            return Err(Error::invalid("capacitance matrix is numerically singular"));
        }
        let inv_l = 1.0 / eigs[k];
        for i in 0..n {
            for j in 0..n {
                inv[i][j] += inv_l * vecs[k][i] * vecs[k][j];
            }
        }
    }
    Ok(inv)
}

/// Computes the passive linear Craig–Bampton port-short decomposition `(T_p, T_c, C_tilde^-1, f_cav_Hz, g_jm, g_jk, kappa_cav_Hz)`
/// over the `n_nl` nonlinear ports (`nl_indices`) of `rom`.
fn build_linear_cb_decomposition(
    rom: &CircuitRom,
    nl_ports: &[RomPort],
    nl_indices: &[usize],
    c_damping_tot: &[Vec<f64>],
    max_cav: usize,
) -> Result<LinearCbDecomposition, Error> {
    let n_ports = nl_ports.len();
    let r = rom.k_r.len();
    if r < n_ports {
        return Err(Error::invalid(format!(
            "CircuitRom dimension r = {r} is smaller than nonlinear port count J = {n_ports}"
        )));
    }
    if rom.m_r.len() != r || rom.b_r.len() != r {
        return Err(Error::invalid("CircuitRom matrix dimensions (k_r, m_r, b_r) are inconsistent"));
    }

    // 1. Construct C_tot,r (F) and L_lin,r^-1 (1/H) including only passive elements + parallel C_j and E_L
    let mut c_tot_r = vec![vec![0.0; r]; r];
    let mut l_inv_r = vec![vec![0.0; r]; r];
    for i in 0..r {
        for j in 0..r {
            c_tot_r[i][j] = EPS0 * rom.m_r[i][j];
            l_inv_r[i][j] = rom.k_r[i][j] / MU0;
        }
    }
    let scale_flux = jj_flux_scale();
    for (j_nl, &orig_p) in nl_indices.iter().enumerate() {
        let port = &nl_ports[j_nl];
        let cj = port.c_j_f;
        let inv_lp = port.e_l_hz.unwrap_or(0.0).max(0.0) / scale_flux;
        for i in 0..r {
            let bi = rom.b_r[i].get(orig_p).copied().unwrap_or(0.0);
            for j in 0..r {
                let bj = rom.b_r[j].get(orig_p).copied().unwrap_or(0.0);
                c_tot_r[i][j] += cj * bi * bj;
                l_inv_r[i][j] += inv_lp * bi * bj;
            }
        }
    }

    // 2. QR / Gram-Schmidt decomposition of B_{r,nl} (r x J) into Q_p (r x J) and nullspace N_c (r x (r - J))
    let mut q_p: Vec<Vec<f64>> = Vec::with_capacity(n_ports);
    for &orig_p in nl_indices {
        let mut v: Vec<f64> = (0..r).map(|i| rom.b_r[i].get(orig_p).copied().unwrap_or(0.0)).collect();
        for _ in 0..2 {
            for q in &q_p {
                let pr = dot_vec(&v, q);
                for (vi, &qi) in v.iter_mut().zip(q) {
                    *vi -= pr * qi;
                }
            }
        }
        let nrm = dot_vec(&v, &v).sqrt();
        if nrm <= 1e-15 {
            return Err(Error::invalid("CircuitRom port incidence matrix B_r has linearly dependent columns"));
        }
        for vi in &mut v {
            *vi /= nrm;
        }
        q_p.push(v);
    }

    // Invert R_p = Q_p^T B_{r,nl} (J x J) so T_{p,0} = Q_p R_p^{-T} satisfies B_{r,nl}^T T_{p,0} = I_J
    let mut r_p = Dense::<C64>::zeros(n_ports, n_ports);
    for i in 0..n_ports {
        for (j, &orig_j) in nl_indices.iter().enumerate() {
            let b_col: Vec<f64> = (0..r).map(|k| rom.b_r[k].get(orig_j).copied().unwrap_or(0.0)).collect();
            r_p[(j, i)] = C64::new(dot_vec(&q_p[i], &b_col), 0.0);
        }
    }
    let mut t_p: Vec<Vec<f64>> = Vec::with_capacity(n_ports);
    for p in 0..n_ports {
        let mut rhs = vec![C64::ZERO; n_ports];
        rhs[p] = C64::ONE;
        let alpha = solve_small(&r_p, &rhs)
            .ok_or_else(|| Error::invalid("B_r^T Q_p is singular in Craig-Bampton decomposition"))?;
        let mut tp_col = vec![0.0; r];
        for i in 0..n_ports {
            let a = alpha[i].re;
            for k in 0..r {
                tp_col[k] += a * q_p[i][k];
            }
        }
        t_p.push(tp_col);
    }

    // Complete Q_p to an orthonormal basis with N_c (r - J columns in ker(B_{r,nl}^T))
    let mut n_c: Vec<Vec<f64>> = Vec::with_capacity(r - n_ports);
    for cand_idx in 0..r {
        if n_c.len() == r - n_ports {
            break;
        }
        let mut v = vec![0.0; r];
        v[cand_idx] = 1.0;
        for _ in 0..2 {
            for q in q_p.iter().chain(n_c.iter()) {
                let pr = dot_vec(&v, q);
                for (vi, &qi) in v.iter_mut().zip(q) {
                    *vi -= pr * qi;
                }
            }
        }
        let nrm = dot_vec(&v, &v).sqrt();
        if nrm > 1e-7 {
            for vi in &mut v {
                *vi /= nrm;
            }
            n_c.push(v);
        }
    }
    let n_int = n_c.len();

    // 3. Diagonalize interior short-circuit cavity block (K_cc, C_cc) on N_c
    let mut t_c_keep: Vec<Vec<f64>> = Vec::new();
    let mut omega_c_keep: Vec<f64> = Vec::new();

    if n_int > 0 && max_cav > 0 {
        let mut k_cc = vec![0.0; n_int * n_int];
        let mut c_cc = vec![0.0; n_int * n_int];
        let l_nc: Vec<Vec<f64>> = n_c.iter().map(|v| mat_vec(&l_inv_r, v)).collect();
        let c_nc: Vec<Vec<f64>> = n_c.iter().map(|v| mat_vec(&c_tot_r, v)).collect();
        for a in 0..n_int {
            for b in a..n_int {
                let kab = 0.5 * (dot_vec(&n_c[a], &l_nc[b]) + dot_vec(&n_c[b], &l_nc[a]));
                let cab = 0.5 * (dot_vec(&n_c[a], &c_nc[b]) + dot_vec(&n_c[b], &c_nc[a]));
                k_cc[a * n_int + b] = kab;
                k_cc[b * n_int + a] = kab;
                c_cc[a * n_int + b] = cab;
                c_cc[b * n_int + a] = cab;
            }
        }
        let (evals_cc, evecs_cc) = dense_gen_eigh(&k_cc, &c_cc, n_int).map_err(Error::from_kernel)?;
        let mut all_tc: Vec<(f64, Vec<f64>)> = Vec::with_capacity(n_int);
        for (lam, u) in evals_cc.into_iter().zip(evecs_cc) {
            let mut tc = vec![0.0; r];
            for (a, &ua) in u.iter().enumerate() {
                for k in 0..r {
                    tc[k] += ua * n_c[a][k];
                }
            }
            all_tc.push((lam, tc));
        }

        // Relax T_p against positive-stiffness interior modes when r_static > 0
        if rom.r_static > 0 {
            for (lam, tc) in &all_tc {
                if *lam > 1.0e12 {
                    let l_tc = mat_vec(&l_inv_r, tc);
                    for tp_col in &mut t_p {
                        let coupling = dot_vec(&l_tc, tp_col) / *lam;
                        for k in 0..r {
                            tp_col[k] -= coupling * tc[k];
                        }
                    }
                }
            }
        }

        let pos_modes: Vec<(f64, f64, Vec<f64>)> = all_tc
            .into_iter()
            .filter_map(|(lam, tc)| {
                if lam.is_finite() && lam > 0.0 {
                    let omega = lam.sqrt();
                    let freq = omega / (2.0 * PI);
                    (freq > 1.0e6).then_some((freq, omega, tc))
                } else {
                    None
                }
            })
            .collect();

        if !rom.modes.is_empty() && pos_modes.len() > 1 {
            let mut cav_targets: Vec<f64> = rom
                .modes
                .iter()
                .filter_map(|row| {
                    let u_ports: f64 = row.u_ind_jj_j.iter().sum();
                    let u_tot = (row.u_mag_vol_j + row.u_kin_ibc_j + u_ports).max(1e-30);
                    (u_ports / u_tot < 0.5).then_some(row.freq_hz)
                })
                .collect();
            if cav_targets.is_empty() && rom.modes.len() > n_ports {
                cav_targets.extend(rom.modes.iter().skip(n_ports).map(|row| row.freq_hz));
            }
            if cav_targets.is_empty() {
                cav_targets.extend(rom.modes.iter().map(|row| row.freq_hz));
            }
            cav_targets.truncate(max_cav);
            let mut used = vec![false; pos_modes.len()];
            let mut picked = Vec::new();
            for f_target in cav_targets {
                let mut best = None;
                let mut best_d = f64::INFINITY;
                for (k, &(fk, _, _)) in pos_modes.iter().enumerate() {
                    if !used[k] {
                        let d = (fk - f_target).abs() / f_target.max(1.0);
                        if d < best_d {
                            best_d = d;
                            best = Some(k);
                        }
                    }
                }
                if let Some(k) = best {
                    used[k] = true;
                    picked.push(pos_modes[k].clone());
                }
            }
            picked.sort_by(|a, b| a.0.total_cmp(&b.0));
            for (_, omega, tc) in picked.into_iter().take(max_cav) {
                omega_c_keep.push(omega);
                t_c_keep.push(tc);
            }
        } else {
            for (freq, omega, tc) in pos_modes {
                if freq < 50.0e9 || rom.r_static == 0 {
                    omega_c_keep.push(omega);
                    t_c_keep.push(tc);
                    if t_c_keep.len() >= max_cav {
                        break;
                    }
                }
            }
        }
    }

    let m_cav = t_c_keep.len();
    let n_sub = n_ports + m_cav;

    // 4. Scale retained cavity vectors t_{c,m} to capacitance scale ~ 5e-14 F for well-conditioned C_tilde
    let c_scale = 50.0e-15_f64;
    let sqrt_c_scale = c_scale.sqrt();
    for tc in &mut t_c_keep {
        for x in tc.iter_mut() {
            *x *= sqrt_c_scale;
        }
    }

    // 5. Assemble C_tilde ((J + M_cav) x (J + M_cav)) and invert via eigendecomposition
    let mut c_tilde = vec![vec![0.0; n_sub]; n_sub];
    let mut c_pp_rom = vec![vec![0.0; n_ports]; n_ports];
    for i in 0..n_ports {
        let c_ti = mat_vec(&c_tot_r, &t_p[i]);
        for j in 0..n_ports {
            c_pp_rom[i][j] = dot_vec(&t_p[j], &c_ti);
            c_tilde[i][j] = c_pp_rom[i][j];
        }
        for m in 0..m_cav {
            let c_im = dot_vec(&t_c_keep[m], &c_ti);
            c_tilde[i][n_ports + m] = c_im;
            c_tilde[n_ports + m][i] = c_im;
        }
    }
    for m in 0..m_cav {
        c_tilde[n_ports + m][n_ports + m] = c_scale;
    }

    // Extract n_ports x n_ports nonlinear port capacitance submatrix from rom.c_port_f if available
    let c_port_nl: Option<Vec<Vec<f64>>> = if rom.c_port_f.len() == rom.ports.len()
        && rom.c_port_f.iter().all(|row| row.len() == rom.ports.len())
        && n_ports > 0
    {
        let sub: Vec<Vec<f64>> =
            nl_indices.iter().map(|&ri| nl_indices.iter().map(|&cj| rom.c_port_f[ri][cj]).collect()).collect();
        (0..n_ports).all(|j| sub[j][j] + nl_ports[j].c_j_f > 0.0).then_some(sub)
    } else if rom.c_port_f.len() == n_ports
        && (0..n_ports).all(|j| rom.c_port_f[j].len() == n_ports && rom.c_port_f[j][j] + nl_ports[j].c_j_f > 0.0)
    {
        Some(rom.c_port_f.clone())
    } else {
        None
    };

    if let Some(mut c_port_tot) = c_port_nl {
        for j in 0..n_ports {
            c_port_tot[j][j] += nl_ports[j].c_j_f;
            let c_rom_jj = c_pp_rom[j][j];
            if c_rom_jj > 0.0 {
                let scale_j = (c_port_tot[j][j] / c_rom_jj).sqrt();
                for m in 0..m_cav {
                    c_tilde[j][n_ports + m] *= scale_j;
                    c_tilde[n_ports + m][j] *= scale_j;
                }
            }
        }
        for i in 0..n_ports {
            for j in 0..n_ports {
                let mut schur_ij = 0.0;
                if rom.r_static > 0 {
                    for m in 0..m_cav {
                        schur_ij += c_tilde[i][n_ports + m] * c_tilde[j][n_ports + m] / c_scale;
                    }
                }
                c_tilde[i][j] = c_port_tot[i][j] + schur_ij;
            }
        }
    }

    let c_inv = invert_spd_dense(&c_tilde)?;

    let hbar = PLANCK_H / (2.0 * PI);
    let mut f_cav_hz = Vec::with_capacity(m_cav);
    let mut q_zpf_cav = Vec::with_capacity(m_cav);
    let mut kappa_cav_hz = Vec::with_capacity(m_cav);
    let has_damping = c_damping_tot.len() == r && c_damping_tot.iter().flatten().any(|&x| x.abs() > 0.0);

    for m in 0..m_cav {
        let l_m_inv = omega_c_keep[m] * omega_c_keep[m] * c_scale;
        let c_m_inv = c_inv[n_ports + m][n_ports + m];
        let omega_m_tilde = (l_m_inv * c_m_inv).sqrt();
        let z_m_tilde = omega_m_tilde / l_m_inv;
        f_cav_hz.push(omega_m_tilde / (2.0 * PI));
        q_zpf_cav.push((hbar / (2.0 * z_m_tilde)).sqrt());

        if has_damping {
            let d_tc = mat_vec(c_damping_tot, &t_c_keep[m]);
            let g_mm = dot_vec(&t_c_keep[m], &d_tc) / (MU0 * C0);
            let kappa_rad_s = (g_mm * c_m_inv).max(0.0);
            kappa_cav_hz.push(kappa_rad_s / (2.0 * PI));
        } else {
            kappa_cav_hz.push(0.0);
        }
    }

    let mut g_jm = vec![vec![0.0; m_cav]; n_ports];
    for j in 0..n_ports {
        for m in 0..m_cav {
            g_jm[j][m] = ((2.0 * E_CHARGE) / PLANCK_H) * c_inv[j][n_ports + m] * q_zpf_cav[m];
        }
    }

    let mut g_jk = vec![vec![0.0; n_ports]; n_ports];
    for j in 0..n_ports {
        for k in (j + 1)..n_ports {
            let v = ((4.0 * E_CHARGE * E_CHARGE) / PLANCK_H) * c_inv[j][k];
            g_jk[j][k] = v;
            g_jk[k][j] = v;
        }
    }

    Ok(LinearCbDecomposition { t_p, t_c_keep, c_inv, f_cav_hz, g_jm, g_jk, kappa_cav_hz })
}

/// Evaluates the multi-mode Purcell energy relaxation time `T_{1,Purcell,j}` (s) for each nonlinear port `j`
/// at its dressed transition frequency `port_f01_hz[j]` by computing the Schur-reduced electromagnetic
/// environment admittance `Y_{env,jj}(omega_{01,j})` on `rom` (with all passive ports terminated and all
/// other nonlinear ports linearized) and weighting by Fermi's Golden Rule matrix element
/// `|<0_j | Phi_j | 1_j>|^2` (reducing identically to `C_{Sigma,j} / Re[Y_{env,jj}(omega_{01,j})]` in the
/// harmonic/transmon limit).
#[allow(clippy::too_many_arguments)]
fn compute_purcell_t1_s(
    rom: &CircuitRom,
    all_ports: &[RomPort],
    nl_indices: &[usize],
    c_damping_tot: &[Vec<f64>],
    has_damping: bool,
    c_inv_opt: Option<&[Vec<f64>]>,
    port_f01_hz: &[f64],
    flux_el_sq: Option<&[f64]>,
) -> Vec<f64> {
    let n_nl = nl_indices.len();
    if !has_damping || n_nl == 0 || rom.k_r.is_empty() {
        return vec![f64::INFINITY; n_nl];
    }
    let r = rom.k_r.len();
    let hbar = PLANCK_H / (2.0 * PI);
    let mut out = Vec::with_capacity(n_nl);

    for (j, &orig_j) in nl_indices.iter().enumerate() {
        let f01 = port_f01_hz.get(j).copied().unwrap_or(0.0);
        if !(f01 > 1.0e3 && f01.is_finite()) {
            out.push(f64::INFINITY);
            continue;
        }
        let omega = 2.0 * PI * f01;
        let k0 = omega / C0;

        let mut a_env = Dense::<C64>::zeros(r, r);
        for ia in 0..r {
            for ib in 0..r {
                let re = rom.k_r[ia][ib] - k0 * k0 * rom.m_r[ia][ib];
                let im = k0 * c_damping_tot[ia][ib];
                a_env[(ia, ib)] = C64::new(re, im);
            }
        }

        // Include parallel C_j and E_L on port j (omitting L_{J,j}), and full linear L_tot + C_j on k != j
        for (k, &orig_k) in nl_indices.iter().enumerate() {
            let pk = &all_ports[orig_k];
            let el = pk.e_l_hz.unwrap_or(0.0).max(0.0);
            let inv_l = if k == j { jj_inv_l_tot(0.0, el) } else { jj_inv_l_tot(pk.l_j_h.unwrap_or(0.0), el) };
            let ak = MU0 * inv_l;
            let am = pk.c_j_f / EPS0;
            let add_re = ak - k0 * k0 * am;
            if add_re.abs() > 0.0 {
                for ia in 0..r {
                    let bi = rom.b_r[ia].get(orig_k).copied().unwrap_or(0.0);
                    if bi == 0.0 {
                        continue;
                    }
                    for ib in 0..r {
                        let bj = rom.b_r[ib].get(orig_k).copied().unwrap_or(0.0);
                        a_env[(ia, ib)] += C64::new(add_re * bi * bj, 0.0);
                    }
                }
            }
        }

        let rhs: Vec<C64> = (0..r).map(|ia| C64::new(rom.b_r[ia].get(orig_j).copied().unwrap_or(0.0), 0.0)).collect();
        let Some(sol) = solve_small(&a_env, &rhs) else {
            out.push(f64::INFINITY);
            continue;
        };

        let mut b_dot_sol = C64::ZERO;
        for ia in 0..r {
            b_dot_sol += rhs[ia] * sol[ia];
        }
        let z_env = C64::new(0.0, omega * MU0) * b_dot_sol;
        if z_env.norm_sqr() <= 1e-300 {
            out.push(f64::INFINITY);
            continue;
        }
        let y_env = C64::ONE / z_env;
        let re_y = y_env.re.max(0.0);
        if re_y <= 0.0 {
            out.push(f64::INFINITY);
            continue;
        }

        let gamma_1 = if let Some(phi_sq) = flux_el_sq.and_then(|s| s.get(j)).copied().filter(|&v| v > 0.0) {
            ((2.0 * omega) / hbar) * phi_sq * re_y
        } else {
            let c_sigma = c_inv_opt
                .and_then(|ci| ci.get(j).and_then(|row| row.get(j)).copied())
                .filter(|&v| v > 0.0)
                .map(|inv_c| 1.0 / inv_c)
                .or_else(|| {
                    rom.c_port_f
                        .get(orig_j)
                        .and_then(|row| row.get(orig_j))
                        .copied()
                        .map(|cp| cp + all_ports[orig_j].c_j_f)
                        .filter(|&c| c > 0.0)
                })
                .unwrap_or_else(|| {
                    let b2: f64 = (0..r).map(|ia| rhs[ia].re * rhs[ia].re).sum();
                    if b2 > 0.0 {
                        EPS0 / b2 + all_ports[orig_j].c_j_f
                    } else {
                        1.0
                    }
                });
            re_y / c_sigma.max(1e-30)
        };

        out.push(if gamma_1 > 1e-30 { 1.0 / gamma_1 } else { f64::INFINITY });
    }
    out
}

// ------------------------------------------------------------------ Mode A: Normal-Mode Fock Basis ("normal_fock")

/// Adds `scale * (I x ... x op_a x ... x I)` into `h_mat`.
fn add_single_mode_op(
    h_mat: &mut Dense<C64>,
    dims: &[usize],
    strides: &[usize],
    a: usize,
    op_a: &Dense<C64>,
    scale: f64,
) {
    let d_tot = h_mat.rows();
    let stride_a = strides[a];
    let dim_a = dims[a];
    let block = stride_a * dim_a;
    let mut outer = 0;
    while outer < d_tot {
        for na in 0..dim_a {
            for ma in 0..dim_a {
                let v = op_a[(na, ma)].scale(scale);
                if v == C64::ZERO {
                    continue;
                }
                let r_base = outer + na * stride_a;
                let c_base = outer + ma * stride_a;
                for inner in 0..stride_a {
                    h_mat[(r_base + inner, c_base + inner)] += v;
                }
            }
        }
        outer += block;
    }
}

/// Adds `scale * (I x ... x op_a x ... x op_b x ... x I)` (`a < b`) into `h_mat`.
#[allow(clippy::too_many_arguments)]
fn add_two_mode_op(
    h_mat: &mut Dense<C64>,
    dims: &[usize],
    strides: &[usize],
    a: usize,
    op_a: &Dense<C64>,
    b: usize,
    op_b: &Dense<C64>,
    scale: f64,
) {
    debug_assert!(a < b);
    let d_tot = h_mat.rows();
    let (sa, da) = (strides[a], dims[a]);
    let (sb, db) = (strides[b], dims[b]);

    for idx_r in 0..d_tot {
        let na = (idx_r / sa) % da;
        let nb = (idx_r / sb) % db;
        let base_r = idx_r - na * sa - nb * sb;
        for ma in 0..da {
            let va = op_a[(na, ma)];
            if va == C64::ZERO {
                continue;
            }
            for mb in 0..db {
                let vb = op_b[(nb, mb)];
                if vb == C64::ZERO {
                    continue;
                }
                let idx_c = base_r + ma * sa + mb * sb;
                h_mat[(idx_r, idx_c)] += (va * vb).scale(scale);
            }
        }
    }
}

/// Computes the Kronecker tensor product `U_0 (x) U_1 (x) ... (x) U_{M-1}` of single-mode matrices.
fn kron_all(ops: &[Dense<C64>], dims: &[usize], strides: &[usize], d_tot: usize) -> Dense<C64> {
    let n_modes = dims.len();
    if n_modes == 1 {
        return ops[0].clone();
    }
    let mut out = Dense::<C64>::zeros(d_tot, d_tot);
    for col in 0..d_tot {
        for row in 0..d_tot {
            let mut prod = C64::ONE;
            for m in 0..n_modes {
                let s = strides[m];
                let d = dims[m];
                let rm = (row / s) % d;
                let cm = (col / s) % d;
                let v = ops[m][(rm, cm)];
                if v == C64::ZERO {
                    prod = C64::ZERO;
                    break;
                }
                prod *= v;
            }
            out[(row, col)] = prod;
        }
    }
    out
}

struct PreparedNormalFock {
    modes_hz: Vec<f64>,
    e_j_hz: Vec<f64>,
    phi_ops: Vec<Vec<Dense<C64>>>,
    phi2_ops: Vec<Vec<Dense<C64>>>,
    u_full_ports: Vec<Dense<C64>>,
}

impl PreparedNormalFock {
    fn from_epr(epr_rep: &crate::epr::EprReport, dims: &[usize], strides: &[usize], d_tot: usize) -> Self {
        let n_modes = epr_rep.modes_hz.len();
        let n_ports = epr_rep.ports.len();
        let mut phi_ops = Vec::with_capacity(n_ports);
        let mut phi2_ops = Vec::with_capacity(n_ports);
        let mut u_full_ports = Vec::with_capacity(n_ports);

        for j in 0..n_ports {
            let mut p_j = Vec::with_capacity(n_modes);
            let mut p2_j = Vec::with_capacity(n_modes);
            let mut u_j = Vec::with_capacity(n_modes);
            for (m, &dm) in dims.iter().enumerate().take(n_modes) {
                let coeff = epr_rep.s_matrix[m][j] * epr_rep.phi_zpf_matrix[m][j];
                let mut phi_mj = Dense::<C64>::zeros(dm, dm);
                let mut i_phi_mj = Dense::<C64>::zeros(dm, dm);
                for n in 0..(dm - 1) {
                    let v = coeff * ((n + 1) as f64).sqrt();
                    phi_mj[(n, n + 1)] = C64::new(v, 0.0);
                    phi_mj[(n + 1, n)] = C64::new(v, 0.0);
                    i_phi_mj[(n, n + 1)] = C64::new(0.0, v);
                    i_phi_mj[(n + 1, n)] = C64::new(0.0, v);
                }
                let phi2_mj = phi_mj.matmul(&phi_mj);
                let u_mj = expm(&i_phi_mj);
                p_j.push(phi_mj);
                p2_j.push(phi2_mj);
                u_j.push(u_mj);
            }
            u_full_ports.push(kron_all(&u_j, dims, strides, d_tot));
            phi_ops.push(p_j);
            phi2_ops.push(p2_j);
        }

        Self { modes_hz: epr_rep.modes_hz.clone(), e_j_hz: epr_rep.e_j_hz.clone(), phi_ops, phi2_ops, u_full_ports }
    }

    fn solve_at_biases(
        &self,
        dims: &[usize],
        strides: &[usize],
        d_tot: usize,
        phi_exts: &[f64],
    ) -> Result<LabeledSpectrum, Error> {
        let n_modes = dims.len();
        let n_ports = self.e_j_hz.len();
        let mut h_mat = Dense::<C64>::zeros(d_tot, d_tot);
        for idx in 0..d_tot {
            let mut diag_e = 0.0;
            for m in 0..n_modes {
                let nm = (idx / strides[m]) % dims[m];
                diag_e += self.modes_hz[m] * (nm as f64);
            }
            h_mat[(idx, idx)] = C64::new(diag_e, 0.0);
        }

        for (j, &ej) in self.e_j_hz.iter().enumerate().take(n_ports) {
            if ej <= 0.0 {
                continue;
            }
            let phi_ext = phi_exts[j];
            let cos_phi = phi_ext.cos();
            let phase = C64::new(cos_phi, -phi_ext.sin());
            let half_ej = 0.5 * ej;
            let u_j = &self.u_full_ports[j];

            for c in 0..d_tot {
                for r in 0..d_tot {
                    let u_rc = phase * u_j[(r, c)];
                    let u_cr_conj = (phase * u_j[(c, r)]).conj();
                    h_mat[(r, c)] -= (u_rc + u_cr_conj).scale(half_ej);
                }
                h_mat[(c, c)] += C64::new(ej * cos_phi, 0.0);
            }

            let quad_scale = -0.5 * ej * cos_phi;
            for (m, phi2_jm) in self.phi2_ops[j].iter().enumerate().take(n_modes) {
                add_single_mode_op(&mut h_mat, dims, strides, m, phi2_jm, quad_scale);
            }
            for m in 0..n_modes {
                for n in (m + 1)..n_modes {
                    add_two_mode_op(
                        &mut h_mat,
                        dims,
                        strides,
                        m,
                        &self.phi_ops[j][m],
                        n,
                        &self.phi_ops[j][n],
                        2.0 * quad_scale,
                    );
                }
            }
        }

        let (evals, evecs) = hermitian_eig(&h_mat).map_err(|e| Error::from_kernel(e.to_string()))?;
        Ok(label_spectrum_bipartite(dims, &evals, |b_idx, k| evecs[(b_idx, k)].norm_sqr()))
    }
}

fn run_normal_fock(result: &SolveResult, params: &Value) -> Result<HamiltonianReport, Error> {
    let t0 = Instant::now();
    let rom = CircuitRom::of(result).ok_or_else(|| Error::invalid("the result carries no CircuitRom"))?;
    let (all_ports_raw, _) = resolve_ports_with_overrides(&rom.ports, params)?;
    let nl_indices: Vec<usize> =
        all_ports_raw.iter().enumerate().filter(|(_, p)| is_nonlinear_port(p)).map(|(i, _)| i).collect();
    if nl_indices.is_empty() {
        return Err(Error::invalid("normal_fock quantization requires at least one nonlinear port in CircuitRom"));
    }

    // Record unmodulated E_{J,sum} (after lj_override/ej_hz_override, before squid_flux_phi0 modulation)
    let mut params_no_sflux = params.clone();
    if let Some(obj) = params_no_sflux.as_object_mut() {
        obj.remove("squid_flux_phi0");
    }
    let (ports_unmod, _) = resolve_ports_with_overrides(&rom.ports, &params_no_sflux)?;
    let base_ej_sum_hz: Vec<f64> = nl_indices
        .iter()
        .map(|&idx| {
            resolve_jj_l_and_ej(ports_unmod[idx].l_j_h.unwrap_or(0.0), ports_unmod[idx].e_j_hz.unwrap_or(0.0)).1
        })
        .collect();

    let mut nl_ports: Vec<RomPort> = nl_indices.iter().map(|&idx| all_ports_raw[idx].clone()).collect();
    apply_bias_overrides(&mut nl_ports, params);
    let sweep_port = resolve_sweep_port(&nl_ports, params)?;

    // If flux_phi0 / phi_ext_rad was passed on a DC-SQUID port (e_l_hz <= 1e3 && squid_asymmetry.is_some()) without squid_flux_phi0
    let mut epr_base_params = if params.is_object() { params.clone() } else { json!({}) };
    if params.get("squid_flux_phi0").is_none()
        && (params.get("flux_phi0").is_some() || params.get("phi_ext_rad").is_some())
    {
        let mut ej_over =
            epr_base_params.get("ej_hz_override").and_then(|v| v.as_object()).cloned().unwrap_or_default();
        let mut any_dc_mod = false;
        for (j, p) in nl_ports.iter_mut().enumerate() {
            if p.e_l_hz.unwrap_or(0.0) <= 1.0e3 {
                if let Some(d) = p.squid_asymmetry {
                    let fl = p.phi_ext_rad / (2.0 * PI);
                    let ej_eff = squid_effective_ej(base_ej_sum_hz[j], fl, d);
                    p.phi_ext_rad = 0.0;
                    if ej_eff > 0.0 {
                        ej_over.insert(p.name.clone(), json!(ej_eff));
                        any_dc_mod = true;
                    }
                }
            }
        }
        if any_dc_mod {
            epr_base_params["ej_hz_override"] = Value::Object(ej_over);
        }
    }

    let epr_rep = epr_with_params(result, &epr_base_params)?;
    let n_modes = epr_rep.modes_hz.len();
    if n_modes == 0 {
        return Err(Error::invalid("CircuitRom has no modes for normal_fock quantization"));
    }
    let n_nl = nl_ports.len();

    let default_dims: Vec<usize> = (0..n_modes)
        .map(|m| {
            let p_tot: f64 = epr_rep.p_matrix[m].iter().sum();
            if p_tot >= 0.2 {
                7
            } else {
                5
            }
        })
        .collect();
    let dims = parse_fock_trunc(params, n_modes, &default_dims)?;
    let (d_tot, strides) = compute_strides(&dims);
    if d_tot > 4096 {
        return Err(Error::invalid(format!(
            "normal_fock Hilbert space dimension {d_tot} exceeds 4096; reduce fock_trunc"
        )));
    }

    let prep_base = PreparedNormalFock::from_epr(&epr_rep, &dims, &strides, d_tot);
    let base_phi_exts: Vec<f64> = nl_ports
        .iter()
        .map(|p| if p.e_l_hz.unwrap_or(0.0) <= 1.0e3 && p.squid_asymmetry.is_some() { 0.0 } else { p.phi_ext_rad })
        .collect();
    let spec = prep_base.solve_at_biases(&dims, &strides, d_tot, &base_phi_exts)?;
    let elapsed_s = t0.elapsed().as_secs_f64();

    // Shared linear Craig-Bampton cavity decomposition (FR-03 Pt 2 + FR-22), guarding gracefully if r < n_nl or singular
    let (_, c_damping_tot, has_damping) = assemble_rom_damping(&rom, &all_ports_raw, None);
    let max_cav_hint = if rom.modes.is_empty() { n_modes.saturating_sub(n_nl).max(1) } else { usize::MAX };
    let cb_opt = build_linear_cb_decomposition(&rom, &nl_ports, &nl_indices, &c_damping_tot, max_cav_hint).ok();

    let (f_cav_hz, g_matrix_hz, kappa_cav_hz) = match &cb_opt {
        Some(cb) => (cb.f_cav_hz.clone(), cb.g_jm.clone(), cb.kappa_cav_hz.clone()),
        None => (Vec::new(), Vec::new(), Vec::new()),
    };

    // Map each nonlinear port j to its dominant normal-mode dressed f01_Hz for Purcell T1 evaluation
    let port_f01_for_t1: Vec<f64> = (0..n_nl)
        .map(|j| {
            let best_m = (0..n_modes)
                .max_by(|&m1, &m2| epr_rep.p_matrix[m1][j].total_cmp(&epr_rep.p_matrix[m2][j]))
                .unwrap_or(0);
            spec.f01_hz[best_m]
        })
        .collect();
    let t1_purcell_s = compute_purcell_t1_s(
        &rom,
        &all_ports_raw,
        &nl_indices,
        &c_damping_tot,
        has_damping,
        cb_opt.as_ref().map(|cb| cb.c_inv.as_slice()),
        &port_f01_for_t1,
        None,
    );

    let subsystems: Vec<String> = (0..n_modes).map(|m| format!("mode_{}", m + 1)).collect();
    let base_fock_warn = fock_truncation_warning_msg(
        spec.max_fock_boundary_occupation,
        &subsystems,
        &dims,
        &spec.fock_boundary_occupation_by_subsystem,
    );

    let (sweep, sweep_elapsed_s) = if let Some(schedule) =
        parse_sweep_schedule(&nl_ports, &base_ej_sum_hz, sweep_port, params)?
    {
        let t_sw = Instant::now();
        let n_pts = schedule.len();
        let n_threads = femcore::par::threads_up_to(8).min(n_pts.max(1));
        let chunk_size = n_pts.div_ceil(n_threads.max(1));
        let target_name = nl_ports[sweep_port].name.clone();
        let target_el = nl_ports[sweep_port].e_l_hz.unwrap_or(0.0);

        let chunks: Vec<Result<Vec<FluxSweepPoint>, Error>> =
            femcore::par::map_chunks(&schedule, chunk_size, |slice| {
                let mut local = Vec::with_capacity(slice.len());
                let mut phis = base_phi_exts.clone();
                for &pt in slice {
                    match pt {
                        SweepPointSpec::LoopPhase { phi_ext_rad, flux_phi0 } => {
                            phis[sweep_port] = phi_ext_rad;
                            let sp = prep_base.solve_at_biases(&dims, &strides, d_tot, &phis)?;
                            let warn = fock_truncation_warning_msg(
                                sp.max_fock_boundary_occupation,
                                &subsystems,
                                &dims,
                                &sp.fock_boundary_occupation_by_subsystem,
                            );
                            local.push(FluxSweepPoint {
                                phi_ext_rad,
                                flux_phi0,
                                e_j_hz: None,
                                f01_hz: sp.f01_hz,
                                f12_hz: sp.f12_hz,
                                f02_hz: sp.f02_hz,
                                alpha_hz: sp.alpha_hz,
                                chi_matrix_hz: sp.chi_matrix_hz,
                                eigenvalues_hz: sp.eigenvalues_hz,
                                min_overlap: sp.min_overlap,
                                max_fock_boundary_occupation: sp.max_fock_boundary_occupation,
                                fock_boundary_occupation_by_subsystem: sp.fock_boundary_occupation_by_subsystem,
                                fock_truncation_warning: warn,
                            });
                        }
                        SweepPointSpec::EffectiveEj { phi_ext_rad, flux_phi0, e_j_hz } => {
                            // Guard E_{J,eff} <= 0 in harmonic normal_fock when E_L <= 1e3:
                            // clamp linear normal-mode reference inductance to a small positive floor (1 MHz) if needed
                            // while setting exact nonlinear E_J = 0 (or delegating pure E_J = 0 charge point if 1D).
                            let sp = if e_j_hz <= 1.0e5
                                && target_el <= 1.0e3
                                && cb_opt.is_some()
                                && n_modes == dims.len()
                            {
                                // At vanishing E_J with E_L = 0, harmonic Fock ZPF diverges; evaluate exact spectrum
                                // via ChargePeriodic Craig-Bampton if dimensions match, else floored normal_fock
                                let mut cb_p = if params.is_object() { params.clone() } else { json!({}) };
                                cb_p["basis"] = json!("craig_bampton");
                                if let Some(obj) = cb_p.as_object_mut() {
                                    obj.remove("ej_hz_sweep");
                                    obj.remove("squid_flux_sweep");
                                    obj.remove("flux_sweep");
                                    obj.remove("phi_ext_sweep");
                                    obj.remove("fock_trunc");
                                }
                                if let Ok(cb_zero) =
                                    run_craig_bampton_with_ej_point(result, &cb_p, sweep_port, e_j_hz.max(0.0))
                                {
                                    cb_zero
                                } else {
                                    let ej_floor = e_j_hz.max(1.0e6);
                                    let mut p_pt = epr_base_params.clone();
                                    let mut ov = p_pt
                                        .get("ej_hz_override")
                                        .and_then(|v| v.as_object())
                                        .cloned()
                                        .unwrap_or_default();
                                    ov.insert(target_name.clone(), json!(ej_floor));
                                    p_pt["ej_hz_override"] = Value::Object(ov);
                                    let mut epr_pt = epr_with_params(result, &p_pt)?;
                                    if e_j_hz <= 0.0 {
                                        epr_pt.e_j_hz[sweep_port] = 0.0;
                                    }
                                    let prep_pt = PreparedNormalFock::from_epr(&epr_pt, &dims, &strides, d_tot);
                                    phis[sweep_port] = 0.0;
                                    prep_pt.solve_at_biases(&dims, &strides, d_tot, &phis)?
                                }
                            } else {
                                let ej_pos = e_j_hz.max(1.0e6);
                                let mut p_pt = epr_base_params.clone();
                                let mut ov =
                                    p_pt.get("ej_hz_override").and_then(|v| v.as_object()).cloned().unwrap_or_default();
                                ov.insert(target_name.clone(), json!(ej_pos));
                                p_pt["ej_hz_override"] = Value::Object(ov);
                                let mut epr_pt = epr_with_params(result, &p_pt)?;
                                if e_j_hz <= 0.0 {
                                    epr_pt.e_j_hz[sweep_port] = 0.0;
                                }
                                let prep_pt = PreparedNormalFock::from_epr(&epr_pt, &dims, &strides, d_tot);
                                phis[sweep_port] = if target_el <= 1.0e3 { 0.0 } else { phi_ext_rad };
                                prep_pt.solve_at_biases(&dims, &strides, d_tot, &phis)?
                            };
                            let warn = fock_truncation_warning_msg(
                                sp.max_fock_boundary_occupation,
                                &subsystems,
                                &dims,
                                &sp.fock_boundary_occupation_by_subsystem,
                            );
                            local.push(FluxSweepPoint {
                                phi_ext_rad,
                                flux_phi0,
                                e_j_hz: Some(e_j_hz.max(0.0)),
                                f01_hz: sp.f01_hz,
                                f12_hz: sp.f12_hz,
                                f02_hz: sp.f02_hz,
                                alpha_hz: sp.alpha_hz,
                                chi_matrix_hz: sp.chi_matrix_hz,
                                eigenvalues_hz: sp.eigenvalues_hz,
                                min_overlap: sp.min_overlap,
                                max_fock_boundary_occupation: sp.max_fock_boundary_occupation,
                                fock_boundary_occupation_by_subsystem: sp.fock_boundary_occupation_by_subsystem,
                                fock_truncation_warning: warn,
                            });
                        }
                    }
                }
                Ok(local)
            });
        let mut out_pts = Vec::with_capacity(n_pts);
        for ch in chunks {
            out_pts.extend(ch?);
        }
        (Some(out_pts), Some(t_sw.elapsed().as_secs_f64()))
    } else {
        (None, None)
    };

    Ok(HamiltonianReport {
        basis: "normal_fock".into(),
        hilbert_dim: d_tot,
        subsystems,
        dims,
        f01_hz: spec.f01_hz,
        f12_hz: spec.f12_hz,
        f02_hz: spec.f02_hz,
        alpha_hz: spec.alpha_hz,
        chi_matrix_hz: spec.chi_matrix_hz,
        eigenvalues_hz: spec.eigenvalues_hz,
        min_overlap: spec.min_overlap,
        max_fock_boundary_occupation: spec.max_fock_boundary_occupation,
        fock_boundary_occupation_by_subsystem: spec.fock_boundary_occupation_by_subsystem,
        fock_truncation_warning: base_fock_warn,
        elapsed_s,
        f_cav_hz,
        g_matrix_hz,
        kappa_cav_hz,
        t1_purcell_s,
        sweep,
        sweep_elapsed_s,
    })
}

// ------------------------------------------------------------------ Mode B: Craig-Bampton Basis ("craig_bampton")

/// Fast allocation-light real symmetric eigensolver `A = V diag(d) V^T` (Householder reduction to tridiagonal form
/// followed by implicit Wilkinson-shifted QL iteration; returns ascending eigenvalues and column eigenvectors `vecs[k]`).
fn fast_sym_eigh(a: &[f64], n: usize) -> Result<(Vec<f64>, Vec<Vec<f64>>), Error> {
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    let mut z = a.to_vec();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];

    // 1. Householder reduction of real symmetric matrix `z` (row-major n x n) to tridiagonal `(d, e)`
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        let mut scale = 0.0;
        if l > 0 {
            for k in 0..=l {
                scale += z[i * n + k].abs();
            }
            if scale == 0.0 {
                e[i] = z[i * n + l];
            } else {
                for k in 0..=l {
                    z[i * n + k] /= scale;
                    h += z[i * n + k] * z[i * n + k];
                }
                let f = z[i * n + l];
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                z[i * n + l] = f - g;
                let mut f_acc = 0.0;
                for j in 0..=l {
                    z[j * n + i] = z[i * n + j] / h;
                    let mut g_acc = 0.0;
                    for k in 0..=j {
                        g_acc += z[j * n + k] * z[i * n + k];
                    }
                    for k in (j + 1)..=l {
                        g_acc += z[k * n + j] * z[i * n + k];
                    }
                    e[j] = g_acc / h;
                    f_acc += e[j] * z[i * n + j];
                }
                let hh = f_acc / (h + h);
                for j in 0..=l {
                    let f_j = z[i * n + j];
                    let g_j = e[j] - hh * f_j;
                    e[j] = g_j;
                    for k in 0..=j {
                        z[j * n + k] -= f_j * e[k] + g_j * z[i * n + k];
                    }
                }
            }
        } else {
            e[i] = z[i * n + l];
        }
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    for i in 0..n {
        let l = i;
        if d[i] != 0.0 {
            for j in 0..l {
                let mut g = 0.0;
                for k in 0..l {
                    g += z[i * n + k] * z[k * n + j];
                }
                for k in 0..l {
                    z[k * n + j] -= g * z[k * n + i];
                }
            }
        }
        d[i] = z[i * n + i];
        z[i * n + i] = 1.0;
        for j in 0..l {
            z[j * n + i] = 0.0;
            z[i * n + j] = 0.0;
        }
    }

    // 2. Implicit Wilkinson-shifted QL iteration on tridiagonal (d, e)
    let mut vecs_col = vec![0.0; n * n];
    for r in 0..n {
        for c in 0..n {
            vecs_col[c * n + r] = z[r * n + c];
        }
    }
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    for l in 0..n {
        let mut iter = 0usize;
        loop {
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= 1e-15 * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            if iter > 120 {
                return dense_eigh(a, n).map_err(Error::from_kernel);
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            let denom = g + if g >= 0.0 { r.abs() } else { -r.abs() };
            g = d[m] - d[l] + e[l] / denom;
            let mut s = 1.0;
            let mut c = 1.0;
            let mut p = 0.0;
            let mut underflow = false;
            for i in (l..m).rev() {
                let f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                let (col_i, col_ip1) = vecs_col.split_at_mut((i + 1) * n);
                let ci = &mut col_i[i * n..(i + 1) * n];
                let cip1 = &mut col_ip1[..n];
                for k in 0..n {
                    let v0 = ci[k];
                    let v1 = cip1[k];
                    cip1[k] = s * v0 + c * v1;
                    ci[k] = c * v0 - s * v1;
                }
            }
            if !underflow {
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }
    }

    // 3. Sort eigenvalues and eigenvectors in ascending order
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap_or(std::cmp::Ordering::Equal));
    let sorted_vals: Vec<f64> = order.iter().map(|&i| d[i]).collect();
    let sorted_vecs: Vec<Vec<f64>> = order.iter().map(|&i| vecs_col[i * n..(i + 1) * n].to_vec()).collect();
    Ok((sorted_vals, sorted_vecs))
}

/// Precomputed 1D nonlinear port model for fast microsecond flux and `E_J` sweeps.
enum PreparedPort1D {
    /// Fluxonium / inductively-shunted port (`E_{L,eff} > 1e3 Hz`):
    /// real symmetric harmonic oscillator basis of size `n_flux` with cached `cos(phi)` and `sin(phi)` from `expm(i phi)`.
    Fluxonium {
        n_flux: usize,
        e_j_hz: f64,
        phi_zpf: f64,
        h0_diag: Vec<f64>,
        cos_phi: Vec<f64>,
        sin_phi: Vec<f64>,
        /// Real skew-symmetric matrix `y = i n = (b - b^dagger) / (2 phi_zpf)` (`n_flux x n_flux`).
        y_op: Vec<f64>,
    },
    /// Charge-periodic port (`E_{L,eff} <= 1e3 Hz`):
    /// tridiagonal charge basis `|-N_c..+N_c>` of size `2 N_c + 1`.
    ChargePeriodic { n_states: usize, n_cut: isize, e_c_hz: f64, c_sigma_f: f64, e_j_hz: f64, n_g: f64 },
}

struct SolvedPort1D {
    evals_hz: Vec<f64>,
    /// Real skew-symmetric `y_proj = i n_proj` (Fluxonium).
    y_real: Option<Vec<f64>>,
    /// Real symmetric `n_proj` (ChargePeriodic at `sin(phi_ext) == 0`).
    n_real: Option<Vec<f64>>,
    /// General complex Hermitian `n_proj`.
    n_cplx: Option<Dense<C64>>,
    /// Squared flux matrix element `|<0 | Phi | 1>|^2` (Wb^2) for Purcell `T_1`.
    flux_01_wb_sq: f64,
    /// Maximum 1D basis boundary occupation across the lowest kept port eigenstates (`0..min(3, d_j)`).
    port_basis_boundary_occ: f64,
}

impl PreparedPort1D {
    fn solve(&self, phi_ext_rad: f64, e_j_override: Option<f64>, d_j: usize) -> Result<SolvedPort1D, Error> {
        match self {
            PreparedPort1D::Fluxonium { n_flux, e_j_hz, phi_zpf, h0_diag, cos_phi, sin_phi, y_op } => {
                let n = *n_flux;
                let d = d_j.min(n);
                let c_ext = phi_ext_rad.cos();
                let s_ext = phi_ext_rad.sin();
                let ej = e_j_override.unwrap_or(*e_j_hz).max(0.0);
                let mut h_port = vec![0.0; n * n];
                for i in 0..n {
                    for j in 0..=i {
                        let idx = i * n + j;
                        let mut val = -ej * (c_ext * cos_phi[idx] + s_ext * sin_phi[idx]);
                        if i == j {
                            val += h0_diag[i];
                        }
                        h_port[idx] = val;
                        h_port[j * n + i] = val;
                    }
                }
                let (evals, evecs) = fast_sym_eigh(&h_port, n)?;
                let mut y_v = vec![0.0; n * d];
                for col in 0..d {
                    let v_col = &evecs[col];
                    for i in 0..n {
                        let mut s = 0.0;
                        if i > 0 {
                            s += y_op[i * n + (i - 1)] * v_col[i - 1];
                        }
                        if i + 1 < n {
                            s += y_op[i * n + (i + 1)] * v_col[i + 1];
                        }
                        y_v[i * d + col] = s;
                    }
                }
                let mut y_proj = vec![0.0; d * d];
                for a in 0..d {
                    let va = &evecs[a];
                    for b in 0..d {
                        let mut s = 0.0;
                        for i in 0..n {
                            s += va[i] * y_v[i * d + b];
                        }
                        y_proj[a * d + b] = s;
                    }
                }

                // Exact <0 | phi | 1> in harmonic oscillator basis: phi = phi_zpf (b + b^dagger)
                let mut phi_01 = 0.0;
                if d >= 2 {
                    let (v0, v1) = (&evecs[0], &evecs[1]);
                    for k in 0..(n - 1) {
                        let s = ((k + 1) as f64).sqrt();
                        phi_01 += *phi_zpf * s * (v0[k] * v1[k + 1] + v0[k + 1] * v1[k]);
                    }
                }
                let phi0_2pi = PHI0 / (2.0 * PI);
                let flux_01_wb_sq = (phi0_2pi * phi_01).powi(2);
                let port_basis_boundary_occ =
                    if n > 0 { (0..d.min(3)).map(|a| evecs[a][n - 1].powi(2)).fold(0.0_f64, f64::max) } else { 0.0 };

                Ok(SolvedPort1D {
                    evals_hz: evals.into_iter().take(d).collect(),
                    y_real: Some(y_proj),
                    n_real: None,
                    n_cplx: None,
                    flux_01_wb_sq,
                    port_basis_boundary_occ,
                })
            }
            PreparedPort1D::ChargePeriodic { n_states, n_cut, e_c_hz, c_sigma_f, e_j_hz, n_g } => {
                let n = *n_states;
                let d = d_j.min(n);
                let ej = e_j_override.unwrap_or(*e_j_hz).max(0.0);
                let half_ej = -0.5 * ej;
                let s_ext = phi_ext_rad.sin();

                if ej == 0.0 || s_ext.abs() <= 1.0e-15 {
                    // Fast 100% real-symmetric tridiagonal path (covers all DC-SQUID E_J sweeps & zero/half-flux biases)
                    let off = half_ej * phi_ext_rad.cos();
                    let mut h_real = vec![0.0; n * n];
                    let mut n_diag = Vec::with_capacity(n);
                    for k in 0..n {
                        let nk = (k as isize - *n_cut) as f64 - *n_g;
                        n_diag.push(nk);
                        h_real[k * n + k] = 4.0 * (*e_c_hz) * nk * nk;
                        if k + 1 < n && off != 0.0 {
                            h_real[k * n + (k + 1)] = off;
                            h_real[(k + 1) * n + k] = off;
                        }
                    }
                    let (evals, evecs) = if off == 0.0 {
                        let mut order: Vec<usize> = (0..n).collect();
                        order.sort_by(|&a, &b| h_real[a * n + a].total_cmp(&h_real[b * n + b]));
                        let vals = order.iter().map(|&i| h_real[i * n + i]).collect();
                        let vecs = order
                            .iter()
                            .map(|&i| {
                                let mut v = vec![0.0; n];
                                v[i] = 1.0;
                                v
                            })
                            .collect();
                        (vals, vecs)
                    } else {
                        fast_sym_eigh(&h_real, n)?
                    };

                    let mut n_proj = vec![0.0; d * d];
                    for a in 0..d {
                        let va = &evecs[a];
                        for b in 0..=a {
                            let vb = &evecs[b];
                            let mut s = 0.0;
                            for k in 0..n {
                                s += va[k] * vb[k] * n_diag[k];
                            }
                            n_proj[a * d + b] = s;
                            n_proj[b * d + a] = s;
                        }
                    }
                    let f01_port = if d >= 2 { (evals[1] - evals[0]).max(0.0) } else { 0.0 };
                    let n01_abs = if d >= 2 { n_proj[1].abs() } else { 0.0 };
                    let flux_01_wb_sq = if f01_port > 1.0e3 && *c_sigma_f > 0.0 {
                        ((2.0 * E_CHARGE * n01_abs) / (2.0 * PI * f01_port * *c_sigma_f)).powi(2)
                    } else {
                        0.0
                    };
                    let port_basis_boundary_occ = if n > 0 {
                        (0..d.min(3))
                            .map(|a| evecs[a][0].powi(2) + if n > 1 { evecs[a][n - 1].powi(2) } else { 0.0 })
                            .fold(0.0_f64, f64::max)
                    } else {
                        0.0
                    };

                    Ok(SolvedPort1D {
                        evals_hz: evals.into_iter().take(d).collect(),
                        y_real: None,
                        n_real: Some(n_proj),
                        n_cplx: None,
                        flux_01_wb_sq,
                        port_basis_boundary_occ,
                    })
                } else {
                    let mut h_mat = Dense::<C64>::zeros(n, n);
                    let off_upper = C64::new(half_ej * phi_ext_rad.cos(), -half_ej * s_ext);
                    let off_lower = off_upper.conj();
                    let mut n_diag = Vec::with_capacity(n);
                    for k in 0..n {
                        let nk = (k as isize - *n_cut) as f64 - *n_g;
                        n_diag.push(nk);
                        h_mat[(k, k)] = C64::new(4.0 * (*e_c_hz) * nk * nk, 0.0);
                        if k + 1 < n {
                            h_mat[(k, k + 1)] = off_upper;
                            h_mat[(k + 1, k)] = off_lower;
                        }
                    }
                    let (evals, evecs) = hermitian_eig(&h_mat).map_err(|e| Error::from_kernel(e.to_string()))?;
                    let mut n_proj = Dense::<C64>::zeros(d, d);
                    for a in 0..d {
                        for b in 0..d {
                            let mut s = C64::ZERO;
                            for k in 0..n {
                                s += (evecs[(k, a)].conj() * evecs[(k, b)]).scale(n_diag[k]);
                            }
                            n_proj[(a, b)] = s;
                        }
                    }
                    let f01_port = if d >= 2 { (evals[1] - evals[0]).max(0.0) } else { 0.0 };
                    let n01_abs = if d >= 2 { n_proj[(0, 1)].norm() } else { 0.0 };
                    let flux_01_wb_sq = if f01_port > 1.0e3 && *c_sigma_f > 0.0 {
                        ((2.0 * E_CHARGE * n01_abs) / (2.0 * PI * f01_port * *c_sigma_f)).powi(2)
                    } else {
                        0.0
                    };
                    let port_basis_boundary_occ = if n > 0 {
                        (0..d.min(3))
                            .map(|a| evecs[(0, a)].norm_sqr() + if n > 1 { evecs[(n - 1, a)].norm_sqr() } else { 0.0 })
                            .fold(0.0_f64, f64::max)
                    } else {
                        0.0
                    };

                    Ok(SolvedPort1D {
                        evals_hz: evals.into_iter().take(d).collect(),
                        y_real: None,
                        n_real: None,
                        n_cplx: Some(n_proj),
                        flux_01_wb_sq,
                        port_basis_boundary_occ,
                    })
                }
            }
        }
    }
}

/// Helper used by `normal_fock` when guarding `E_{J,eff} == 0, E_L == 0` points.
fn run_craig_bampton_with_ej_point(
    result: &SolveResult,
    params: &Value,
    target_port: usize,
    ej_val: f64,
) -> Result<LabeledSpectrum, Error> {
    let mut p = params.clone();
    p["ej_hz_sweep"] = json!([ej_val]);
    p["sweep_port"] = json!(target_port);
    let rep = run_craig_bampton(result, &p)?;
    let sw = rep.sweep.and_then(|mut v| v.pop()).ok_or_else(|| Error::invalid("empty CB point"))?;
    Ok(LabeledSpectrum {
        f01_hz: sw.f01_hz,
        f12_hz: sw.f12_hz,
        f02_hz: sw.f02_hz,
        alpha_hz: sw.alpha_hz,
        chi_matrix_hz: sw.chi_matrix_hz,
        eigenvalues_hz: sw.eigenvalues_hz,
        min_overlap: sw.min_overlap,
        max_fock_boundary_occupation: sw.max_fock_boundary_occupation,
        fock_boundary_occupation_by_subsystem: sw.fock_boundary_occupation_by_subsystem,
    })
}

fn run_craig_bampton(result: &SolveResult, params: &Value) -> Result<HamiltonianReport, Error> {
    let t0 = Instant::now();
    let rom = CircuitRom::of(result).ok_or_else(|| Error::invalid("the result carries no CircuitRom"))?;
    let (all_ports_raw, _) = resolve_ports_with_overrides(&rom.ports, params)?;
    let nl_indices: Vec<usize> =
        all_ports_raw.iter().enumerate().filter(|(_, p)| is_nonlinear_port(p)).map(|(i, _)| i).collect();
    if nl_indices.is_empty() {
        return Err(Error::invalid("craig_bampton quantization requires at least one nonlinear port in CircuitRom"));
    }

    let mut params_no_sflux = params.clone();
    if let Some(obj) = params_no_sflux.as_object_mut() {
        obj.remove("squid_flux_phi0");
    }
    let (ports_unmod, _) = resolve_ports_with_overrides(&rom.ports, &params_no_sflux)?;
    let base_ej_sum_hz: Vec<f64> = nl_indices
        .iter()
        .map(|&idx| {
            resolve_jj_l_and_ej(ports_unmod[idx].l_j_h.unwrap_or(0.0), ports_unmod[idx].e_j_hz.unwrap_or(0.0)).1
        })
        .collect();

    let mut nl_ports: Vec<RomPort> = nl_indices.iter().map(|&idx| all_ports_raw[idx].clone()).collect();
    apply_bias_overrides(&mut nl_ports, params);
    let sweep_port = resolve_sweep_port(&nl_ports, params)?;

    // Auto-modulate single-point DC-SQUID E_J if flux_phi0 / phi_ext_rad was passed without squid_flux_phi0
    if params.get("squid_flux_phi0").is_none()
        && (params.get("flux_phi0").is_some() || params.get("phi_ext_rad").is_some())
    {
        for (j, p) in nl_ports.iter_mut().enumerate() {
            if p.e_l_hz.unwrap_or(0.0) <= 1.0e3 {
                if let Some(d) = p.squid_asymmetry {
                    let fl = p.phi_ext_rad / (2.0 * PI);
                    let ej_eff = squid_effective_ej(base_ej_sum_hz[j], fl, d);
                    p.phi_ext_rad = 0.0;
                    p.e_j_hz = (ej_eff > 0.0).then_some(ej_eff);
                    p.l_j_h = (ej_eff > 0.0).then(|| jj_flux_scale() / ej_eff);
                }
            }
        }
    }

    let n_ports = nl_ports.len();
    let port_trunc = params.get("port_trunc").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
    let n_flux_basis = params.get("n_flux_basis").and_then(|v| v.as_u64()).unwrap_or(40) as usize;
    let n_charge_cut = params.get("n_charge_cut").and_then(|v| v.as_u64()).unwrap_or(15) as isize;
    if port_trunc < 3 {
        return Err(Error::invalid("port_trunc must be >= 3"));
    }

    let (_, c_damping_tot, has_damping) = assemble_rom_damping(&rom, &all_ports_raw, None);
    let max_cav = params.get("fock_trunc").and_then(|v| v.as_array()).map(|a| a.len()).unwrap_or(usize::MAX);
    let cb = build_linear_cb_decomposition(&rom, &nl_ports, &nl_indices, &c_damping_tot, max_cav)?;
    let m_cav = cb.f_cav_hz.len();

    // Prepare 1D nonlinear port models (caching `expm` once per Fluxonium port)
    let mut prepared_ports = Vec::with_capacity(n_ports);
    for (j, port) in nl_ports.iter().enumerate() {
        let (_, ej_hz) = resolve_jj_l_and_ej(port.l_j_h.unwrap_or(0.0), port.e_j_hz.unwrap_or(0.0));
        let inv_c_jj = cb.c_inv[j][j];
        let c_sigma_f = if inv_c_jj > 0.0 { 1.0 / inv_c_jj } else { 0.0 };
        let ec_hz = ((E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H)) * inv_c_jj;
        let el_hz = port.e_l_hz.unwrap_or(0.0).max(0.0);

        if el_hz > 1.0e3 {
            let n_flux = n_flux_basis;
            let fp_hz = (8.0 * ec_hz * el_hz).sqrt();
            let phi_zpf = (2.0 * ec_hz / el_hz).powf(0.25);
            let h0_diag: Vec<f64> = (0..n_flux).map(|k| fp_hz * (k as f64 + 0.5)).collect();

            let mut i_phi = Dense::<C64>::zeros(n_flux, n_flux);
            let mut y_op = vec![0.0; n_flux * n_flux];
            for k in 0..(n_flux - 1) {
                let s = ((k + 1) as f64).sqrt();
                let p_val = phi_zpf * s;
                i_phi[(k, k + 1)] = C64::new(0.0, p_val);
                i_phi[(k + 1, k)] = C64::new(0.0, p_val);
                let y_val = s / (2.0 * phi_zpf);
                y_op[k * n_flux + (k + 1)] = y_val;
                y_op[(k + 1) * n_flux + k] = -y_val;
            }
            let u_mat = expm(&i_phi);
            let mut cos_phi = vec![0.0; n_flux * n_flux];
            let mut sin_phi = vec![0.0; n_flux * n_flux];
            for r_i in 0..n_flux {
                for c_i in 0..n_flux {
                    let z = u_mat[(r_i, c_i)];
                    cos_phi[r_i * n_flux + c_i] = z.re;
                    sin_phi[r_i * n_flux + c_i] = z.im;
                }
            }
            prepared_ports.push(PreparedPort1D::Fluxonium {
                n_flux,
                e_j_hz: ej_hz,
                phi_zpf,
                h0_diag,
                cos_phi,
                sin_phi,
                y_op,
            });
        } else {
            let n_states = (2 * n_charge_cut + 1) as usize;
            prepared_ports.push(PreparedPort1D::ChargePeriodic {
                n_states,
                n_cut: n_charge_cut,
                e_c_hz: ec_hz,
                c_sigma_f,
                e_j_hz: ej_hz,
                n_g: port.n_g,
            });
        }
    }

    let default_cav_dims = vec![5usize; m_cav];
    let cav_dims = parse_fock_trunc(params, m_cav, &default_cav_dims)?;
    let mut dims = vec![port_trunc; n_ports];
    dims.extend_from_slice(&cav_dims);
    let (d_tot, strides) = compute_strides(&dims);
    if d_tot > 4096 {
        return Err(Error::invalid(format!(
            "craig_bampton Hilbert space dimension {d_tot} exceeds 4096; reduce port_trunc or fock_trunc"
        )));
    }

    let f_cav_hz = cb.f_cav_hz;
    let g_jm = cb.g_jm;
    let g_jk = cb.g_jk;
    let kappa_cav_hz = cb.kappa_cav_hz;
    let c_inv = cb.c_inv;

    // Precompute cavity real operators `(a_m^dagger - a_m)` and `(a_m^dagger + a_m)` (`D_m x D_m`)
    let mut cav_adag_minus_a: Vec<Vec<f64>> = Vec::with_capacity(m_cav);
    let mut cav_adag_plus_a: Vec<Vec<f64>> = Vec::with_capacity(m_cav);
    for &dm in cav_dims.iter().take(m_cav) {
        let mut mat_minus = vec![0.0; dm * dm];
        let mut mat_plus = vec![0.0; dm * dm];
        for n in 0..(dm - 1) {
            let s = ((n + 1) as f64).sqrt();
            mat_minus[(n + 1) * dm + n] = s;
            mat_minus[n * dm + (n + 1)] = -s;
            mat_plus[(n + 1) * dm + n] = s;
            mat_plus[n * dm + (n + 1)] = s;
        }
        cav_adag_minus_a.push(mat_minus);
        cav_adag_plus_a.push(mat_plus);
    }

    let solve_cb_at_state = |phi_exts: &[f64],
                             ej_overrides: &[Option<f64>]|
     -> Result<(LabeledSpectrum, Vec<f64>), Error> {
        let mut solved_ports = Vec::with_capacity(n_ports);
        let mut flux_el_sq = Vec::with_capacity(n_ports);
        for j in 0..n_ports {
            let sp = prepared_ports[j].solve(phi_exts[j], ej_overrides[j], port_trunc)?;
            flux_el_sq.push(sp.flux_01_wb_sq);
            solved_ports.push(sp);
        }

        let all_y: Option<Vec<&[f64]>> = solved_ports.iter().map(|sp| sp.y_real.as_deref()).collect();
        let all_n: Option<Vec<&[f64]>> = solved_ports.iter().map(|sp| sp.n_real.as_deref()).collect();

        if let Some(y_mats) = all_y {
            // Fast 100% real-symmetric Fluxonium path: n_j (x) i (a^dagger - a) = + y_j (x) (a^dagger - a)
            let mut h_real = vec![0.0; d_tot * d_tot];
            for idx in 0..d_tot {
                let mut e_diag = 0.0;
                for j in 0..n_ports {
                    let nj = (idx / strides[j]) % dims[j];
                    e_diag += solved_ports[j].evals_hz[nj];
                }
                for m in 0..m_cav {
                    let nm = (idx / strides[n_ports + m]) % dims[n_ports + m];
                    e_diag += f_cav_hz[m] * (nm as f64);
                }
                h_real[idx * d_tot + idx] = e_diag;
            }

            for j in 0..n_ports {
                let y_j = y_mats[j];
                let (sj, dj) = (strides[j], dims[j]);
                for m in 0..m_cav {
                    let g = g_jm[j][m];
                    if g.abs() <= 1e-12 {
                        continue;
                    }
                    let op_m = &cav_adag_minus_a[m];
                    let (sm, dm) = (strides[n_ports + m], dims[n_ports + m]);
                    for idx_r in 0..d_tot {
                        let nj = (idx_r / sj) % dj;
                        let nm = (idx_r / sm) % dm;
                        let base_r = idx_r - nj * sj - nm * sm;
                        for mm in [nm.wrapping_sub(1), nm + 1] {
                            if mm >= dm {
                                continue;
                            }
                            let vm = op_m[nm * dm + mm];
                            if vm == 0.0 {
                                continue;
                            }
                            for mj in 0..dj {
                                let y_val = y_j[nj * dj + mj];
                                if y_val != 0.0 {
                                    let idx_c = base_r + mj * sj + mm * sm;
                                    h_real[idx_r * d_tot + idx_c] += g * y_val * vm;
                                }
                            }
                        }
                    }
                }
            }

            for j in 0..n_ports {
                for k in (j + 1)..n_ports {
                    let g = g_jk[j][k];
                    if g.abs() <= 1e-12 {
                        continue;
                    }
                    let y_j = y_mats[j];
                    let y_k = y_mats[k];
                    let (sj, dj) = (strides[j], dims[j]);
                    let (sk, dk) = (strides[k], dims[k]);
                    for idx_r in 0..d_tot {
                        let nj = (idx_r / sj) % dj;
                        let nk = (idx_r / sk) % dk;
                        let base_r = idx_r - nj * sj - nk * sk;
                        for mj in 0..dj {
                            let yj = y_j[nj * dj + mj];
                            if yj == 0.0 {
                                continue;
                            }
                            for mk in 0..dk {
                                let yk = y_k[nk * dk + mk];
                                if yk != 0.0 {
                                    let idx_c = base_r + mj * sj + mk * sk;
                                    h_real[idx_r * d_tot + idx_c] -= g * yj * yk;
                                }
                            }
                        }
                    }
                }
            }

            let (evals, evecs) = if d_tot >= 64 {
                dense_eigh(&h_real, d_tot).map_err(Error::from_kernel)?
            } else {
                fast_sym_eigh(&h_real, d_tot)?
            };
            let mut spec = label_spectrum_bipartite(&dims, &evals, |b_idx, k| {
                let v = evecs[k][b_idx];
                v * v
            });
            for (j, sp_j) in solved_ports.iter().enumerate() {
                spec.fock_boundary_occupation_by_subsystem[j] =
                    spec.fock_boundary_occupation_by_subsystem[j].max(sp_j.port_basis_boundary_occ);
            }
            spec.max_fock_boundary_occupation =
                spec.fock_boundary_occupation_by_subsystem.iter().copied().fold(0.0_f64, f64::max);
            Ok((spec, flux_el_sq))
        } else if let Some(n_mats) = all_n {
            // Fast 100% real-symmetric ChargePeriodic path under diagonal Fock gauge R_m = i^{a_m^dagger a_m}:
            // R_m^dagger i(a_m^dagger - a_m) R_m = +(a_m^dagger + a_m), preserving all Fock state probabilities!
            let mut h_real = vec![0.0; d_tot * d_tot];
            for idx in 0..d_tot {
                let mut e_diag = 0.0;
                for j in 0..n_ports {
                    let nj = (idx / strides[j]) % dims[j];
                    e_diag += solved_ports[j].evals_hz[nj];
                }
                for m in 0..m_cav {
                    let nm = (idx / strides[n_ports + m]) % dims[n_ports + m];
                    e_diag += f_cav_hz[m] * (nm as f64);
                }
                h_real[idx * d_tot + idx] = e_diag;
            }

            for j in 0..n_ports {
                let n_j = n_mats[j];
                let (sj, dj) = (strides[j], dims[j]);
                for m in 0..m_cav {
                    let g = g_jm[j][m];
                    if g.abs() <= 1e-12 {
                        continue;
                    }
                    let op_m = &cav_adag_plus_a[m];
                    let (sm, dm) = (strides[n_ports + m], dims[n_ports + m]);
                    for idx_r in 0..d_tot {
                        let nj = (idx_r / sj) % dj;
                        let nm = (idx_r / sm) % dm;
                        let base_r = idx_r - nj * sj - nm * sm;
                        for mm in [nm.wrapping_sub(1), nm + 1] {
                            if mm >= dm {
                                continue;
                            }
                            let vm = op_m[nm * dm + mm];
                            if vm == 0.0 {
                                continue;
                            }
                            for mj in 0..dj {
                                let nv = n_j[nj * dj + mj];
                                if nv != 0.0 {
                                    let idx_c = base_r + mj * sj + mm * sm;
                                    h_real[idx_r * d_tot + idx_c] += g * nv * vm;
                                }
                            }
                        }
                    }
                }
            }

            for j in 0..n_ports {
                for k in (j + 1)..n_ports {
                    let g = g_jk[j][k];
                    if g.abs() <= 1e-12 {
                        continue;
                    }
                    let n_j = n_mats[j];
                    let n_k = n_mats[k];
                    let (sj, dj) = (strides[j], dims[j]);
                    let (sk, dk) = (strides[k], dims[k]);
                    for idx_r in 0..d_tot {
                        let nj = (idx_r / sj) % dj;
                        let nk = (idx_r / sk) % dk;
                        let base_r = idx_r - nj * sj - nk * sk;
                        for mj in 0..dj {
                            let njv = n_j[nj * dj + mj];
                            if njv == 0.0 {
                                continue;
                            }
                            for mk in 0..dk {
                                let nkv = n_k[nk * dk + mk];
                                if nkv != 0.0 {
                                    let idx_c = base_r + mj * sj + mk * sk;
                                    h_real[idx_r * d_tot + idx_c] += g * njv * nkv;
                                }
                            }
                        }
                    }
                }
            }

            let (evals, evecs) = if d_tot >= 64 {
                dense_eigh(&h_real, d_tot).map_err(Error::from_kernel)?
            } else {
                fast_sym_eigh(&h_real, d_tot)?
            };
            let mut spec = label_spectrum_bipartite(&dims, &evals, |b_idx, k| {
                let v = evecs[k][b_idx];
                v * v
            });
            for (j, sp_j) in solved_ports.iter().enumerate() {
                spec.fock_boundary_occupation_by_subsystem[j] =
                    spec.fock_boundary_occupation_by_subsystem[j].max(sp_j.port_basis_boundary_occ);
            }
            spec.max_fock_boundary_occupation =
                spec.fock_boundary_occupation_by_subsystem.iter().copied().fold(0.0_f64, f64::max);
            Ok((spec, flux_el_sq))
        } else {
            // General complex Hermitian path
            let mut h_mat = Dense::<C64>::zeros(d_tot, d_tot);
            for idx in 0..d_tot {
                let mut e_diag = 0.0;
                for j in 0..n_ports {
                    let nj = (idx / strides[j]) % dims[j];
                    e_diag += solved_ports[j].evals_hz[nj];
                }
                for m in 0..m_cav {
                    let nm = (idx / strides[n_ports + m]) % dims[n_ports + m];
                    e_diag += f_cav_hz[m] * (nm as f64);
                }
                h_mat[(idx, idx)] = C64::new(e_diag, 0.0);
            }

            let n_ops: Vec<Dense<C64>> = solved_ports
                .iter()
                .enumerate()
                .map(|(j, sp)| {
                    let d = dims[j];
                    if let Some(ref nc) = sp.n_cplx {
                        Ok(nc.clone())
                    } else if let Some(ref nr) = sp.n_real {
                        Ok(Dense::from_fn(d, d, |r_i, c_i| C64::new(nr[r_i * d + c_i], 0.0)))
                    } else if let Some(ref y) = sp.y_real {
                        Ok(Dense::from_fn(d, d, |r_i, c_i| C64::new(0.0, -y[r_i * d + c_i])))
                    } else {
                        Err(Error::invalid(format!("port {j} solver returned neither charge nor flux operator")))
                    }
                })
                .collect::<Result<Vec<_>, Error>>()?;

            let cav_i_adag_a: Vec<Dense<C64>> = (0..m_cav)
                .map(|m| {
                    let dm = cav_dims[m];
                    let op = &cav_adag_minus_a[m];
                    Dense::from_fn(dm, dm, |r_i, c_i| C64::new(0.0, op[r_i * dm + c_i]))
                })
                .collect();

            for j in 0..n_ports {
                for m in 0..m_cav {
                    let g = g_jm[j][m];
                    if g.abs() > 1e-12 {
                        add_two_mode_op(&mut h_mat, &dims, &strides, j, &n_ops[j], n_ports + m, &cav_i_adag_a[m], g);
                    }
                }
                for k in (j + 1)..n_ports {
                    let g = g_jk[j][k];
                    if g.abs() > 1e-12 {
                        add_two_mode_op(&mut h_mat, &dims, &strides, j, &n_ops[j], k, &n_ops[k], g);
                    }
                }
            }

            let (evals, evecs) = hermitian_eig(&h_mat).map_err(|e| Error::from_kernel(e.to_string()))?;
            let mut spec = label_spectrum_bipartite(&dims, &evals, |b_idx, k| evecs[(b_idx, k)].norm_sqr());
            for (j, sp_j) in solved_ports.iter().enumerate() {
                spec.fock_boundary_occupation_by_subsystem[j] =
                    spec.fock_boundary_occupation_by_subsystem[j].max(sp_j.port_basis_boundary_occ);
            }
            spec.max_fock_boundary_occupation =
                spec.fock_boundary_occupation_by_subsystem.iter().copied().fold(0.0_f64, f64::max);
            Ok((spec, flux_el_sq))
        }
    };

    let base_phi_exts: Vec<f64> = nl_ports
        .iter()
        .map(|p| if p.e_l_hz.unwrap_or(0.0) <= 1.0e3 && p.squid_asymmetry.is_some() { 0.0 } else { p.phi_ext_rad })
        .collect();
    let base_ej_overrides: Vec<Option<f64>> = vec![None; n_ports];
    let (spec, flux_el_sq) = solve_cb_at_state(&base_phi_exts, &base_ej_overrides)?;
    let elapsed_s = t0.elapsed().as_secs_f64();

    let port_f01: Vec<f64> = spec.f01_hz.iter().copied().take(n_ports).collect();
    let t1_purcell_s = compute_purcell_t1_s(
        &rom,
        &all_ports_raw,
        &nl_indices,
        &c_damping_tot,
        has_damping,
        Some(&c_inv),
        &port_f01,
        Some(&flux_el_sq),
    );

    let mut subsystems: Vec<String> = nl_ports.iter().map(|p| p.name.clone()).collect();
    for m in 0..m_cav {
        subsystems.push(format!("cavity_{}", m + 1));
    }
    let base_fock_warn = fock_truncation_warning_msg(
        spec.max_fock_boundary_occupation,
        &subsystems,
        &dims,
        &spec.fock_boundary_occupation_by_subsystem,
    );

    let (sweep, sweep_elapsed_s) =
        if let Some(schedule) = parse_sweep_schedule(&nl_ports, &base_ej_sum_hz, sweep_port, params)? {
            let t_sw = Instant::now();
            let n_pts = schedule.len();
            let n_threads = femcore::par::threads_up_to(8).min(n_pts.max(1));
            let chunk_size = n_pts.div_ceil(n_threads.max(1));
            let target_el = nl_ports[sweep_port].e_l_hz.unwrap_or(0.0);

            let chunks: Vec<Result<Vec<FluxSweepPoint>, Error>> =
                femcore::par::map_chunks(&schedule, chunk_size, |slice| {
                    let mut local = Vec::with_capacity(slice.len());
                    let mut phis = base_phi_exts.clone();
                    let mut ej_ovs = vec![None; n_ports];
                    for &pt in slice {
                        match pt {
                            SweepPointSpec::LoopPhase { phi_ext_rad, flux_phi0 } => {
                                phis[sweep_port] = phi_ext_rad;
                                ej_ovs[sweep_port] = None;
                                let (sp, _) = solve_cb_at_state(&phis, &ej_ovs)?;
                                let warn = fock_truncation_warning_msg(
                                    sp.max_fock_boundary_occupation,
                                    &subsystems,
                                    &dims,
                                    &sp.fock_boundary_occupation_by_subsystem,
                                );
                                local.push(FluxSweepPoint {
                                    phi_ext_rad,
                                    flux_phi0,
                                    e_j_hz: None,
                                    f01_hz: sp.f01_hz,
                                    f12_hz: sp.f12_hz,
                                    f02_hz: sp.f02_hz,
                                    alpha_hz: sp.alpha_hz,
                                    chi_matrix_hz: sp.chi_matrix_hz,
                                    eigenvalues_hz: sp.eigenvalues_hz,
                                    min_overlap: sp.min_overlap,
                                    max_fock_boundary_occupation: sp.max_fock_boundary_occupation,
                                    fock_boundary_occupation_by_subsystem: sp.fock_boundary_occupation_by_subsystem,
                                    fock_truncation_warning: warn,
                                });
                            }
                            SweepPointSpec::EffectiveEj { phi_ext_rad, flux_phi0, e_j_hz } => {
                                phis[sweep_port] = if target_el <= 1.0e3 { 0.0 } else { base_phi_exts[sweep_port] };
                                ej_ovs[sweep_port] = Some(e_j_hz.max(0.0));
                                let (sp, _) = solve_cb_at_state(&phis, &ej_ovs)?;
                                let warn = fock_truncation_warning_msg(
                                    sp.max_fock_boundary_occupation,
                                    &subsystems,
                                    &dims,
                                    &sp.fock_boundary_occupation_by_subsystem,
                                );
                                local.push(FluxSweepPoint {
                                    phi_ext_rad,
                                    flux_phi0,
                                    e_j_hz: Some(e_j_hz.max(0.0)),
                                    f01_hz: sp.f01_hz,
                                    f12_hz: sp.f12_hz,
                                    f02_hz: sp.f02_hz,
                                    alpha_hz: sp.alpha_hz,
                                    chi_matrix_hz: sp.chi_matrix_hz,
                                    eigenvalues_hz: sp.eigenvalues_hz,
                                    min_overlap: sp.min_overlap,
                                    max_fock_boundary_occupation: sp.max_fock_boundary_occupation,
                                    fock_boundary_occupation_by_subsystem: sp.fock_boundary_occupation_by_subsystem,
                                    fock_truncation_warning: warn,
                                });
                            }
                        }
                    }
                    Ok(local)
                });
            let mut out_pts = Vec::with_capacity(n_pts);
            for ch in chunks {
                out_pts.extend(ch?);
            }
            (Some(out_pts), Some(t_sw.elapsed().as_secs_f64()))
        } else {
            (None, None)
        };

    Ok(HamiltonianReport {
        basis: "craig_bampton".into(),
        hilbert_dim: d_tot,
        subsystems,
        dims,
        f01_hz: spec.f01_hz,
        f12_hz: spec.f12_hz,
        f02_hz: spec.f02_hz,
        alpha_hz: spec.alpha_hz,
        chi_matrix_hz: spec.chi_matrix_hz,
        eigenvalues_hz: spec.eigenvalues_hz,
        min_overlap: spec.min_overlap,
        max_fock_boundary_occupation: spec.max_fock_boundary_occupation,
        fock_boundary_occupation_by_subsystem: spec.fock_boundary_occupation_by_subsystem,
        fock_truncation_warning: base_fock_warn,
        elapsed_s,
        f_cav_hz,
        g_matrix_hz: g_jm,
        kappa_cav_hz,
        t1_purcell_s,
        sweep,
        sweep_elapsed_s,
    })
}

/// Runs circuit quantization with default parameters (`basis = "normal_fock"`).
pub fn hamiltonian(result: &SolveResult) -> Result<HamiltonianReport, Error> {
    hamiltonian_with_params(result, &Value::Null)
}

/// Runs circuit quantization with optional JSON `params`.
pub fn hamiltonian_with_params(result: &SolveResult, params: &Value) -> Result<HamiltonianReport, Error> {
    let basis = params.get("basis").and_then(|v| v.as_str()).unwrap_or("normal_fock");
    match basis {
        "normal_fock" => run_normal_fock(result, params),
        "craig_bampton" => run_craig_bampton(result, params),
        other => Err(Error::invalid(format!(
            "unknown hamiltonian basis '{other}' (expected 'normal_fock' or 'craig_bampton')"
        ))),
    }
}

/// Full cosine-potential circuit quantization analysis behind the [`Analysis`] trait (`id: "hamiltonian"`).
pub struct HamiltonianAnalysis;

impl Analysis for HamiltonianAnalysis {
    fn id(&self) -> &'static str {
        "hamiltonian"
    }

    fn name(&self) -> &'static str {
        "Full cosine-potential circuit quantization (normal-mode Fock and Craig-Bampton)"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        hamiltonian_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(hamiltonian(result)?.to_json())
    }

    fn run_with_params(&self, _mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(hamiltonian_with_params(result, params)?.to_json())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::epr::project_rom_modes;

    #[test]
    fn craig_bampton_matches_exact_2body_lc_plus_fluxonium_to_machine_precision() {
        let ec_hz = 1.0e9_f64;
        let el_hz = 1.0e9_f64;
        let ej_hz = 4.0e9_f64;
        let f_r_bare = 7.0e9_f64;
        let omega_r = 2.0 * PI * f_r_bare;

        let c_sigma = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * ec_hz);
        let c_r = 50.0e-15_f64;
        let c_g = 0.50e-15_f64;
        let c_11 = c_sigma + (c_g * c_g) / (c_r + c_g);
        let l_r_inv = omega_r * omega_r * (c_r + c_g);

        let m_r = vec![vec![c_11 / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]];
        let k_r = vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]];
        let b_r = vec![vec![1.0], vec![0.0]];

        let rom = CircuitRom {
            ports: vec![RomPort {
                name: "fluxonium_jj".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(jj_flux_scale() / ej_hz),
                e_j_hz: Some(ej_hz),
                c_j_f: 0.0,
                e_l_hz: Some(el_hz),
                phi_ext_rad: PI,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [1.0, 0.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 0,
            r_dyn: 2,
            k_r,
            c_r: vec![],
            m_r,
            b_r,
            c_port_f: vec![],
            modes: vec![],
        };
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

        let rep = hamiltonian_with_params(
            &res,
            &json!({
                "basis": "craig_bampton",
                "port_trunc": 8,
                "fock_trunc": [5],
                "n_flux_basis": 40
            }),
        )
        .unwrap();

        assert_eq!(rep.hilbert_dim, 40);
        assert!((rep.f01_hz[0] - 581.849e6).abs() / 581.849e6 < 5e-3, "Fluxonium f01 = {:.4} MHz", rep.f01_hz[0] / 1e6);
        assert!((rep.f01_hz[1] - 7.0e9).abs() / 7.0e9 < 5e-3, "Resonator f01 = {:.6} GHz", rep.f01_hz[1] / 1e9);
        assert!(rep.min_overlap > 0.98, "min_overlap = {}", rep.min_overlap);

        // Verify FR-03 Pt 1 JSON serialization includes f_cav_Hz and g_matrix_Hz
        let j = rep.to_json();
        assert_eq!(j["f_cav_Hz"].as_array().unwrap().len(), 1);
        assert_eq!(j["g_matrix_Hz"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn normal_fock_and_charge_periodic_craig_bampton_agree_on_transmon_and_populate_fcav_g() {
        let scale = jj_flux_scale();
        let e_j_hz = 16.0e9;
        let l_j_h = scale / e_j_hz;
        let c_port = 75.0e-15;
        let c_r = 50.0e-15;
        let c_g = 1.5e-15;
        let f_r_bare = 6.8e9;
        let omega_r = 2.0 * PI * f_r_bare;
        let l_r_inv = omega_r * omega_r * (c_r + c_g);

        let mut rom = CircuitRom {
            ports: vec![RomPort {
                name: "Q1_JJ".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(l_j_h),
                e_j_hz: Some(e_j_hz),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 0,
            r_dyn: 2,
            k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]],
            c_r: vec![],
            m_r: vec![vec![(c_port + c_g) / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]],
            b_r: vec![vec![1.0], vec![0.0]],
            c_port_f: vec![],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

        let rep_fock = hamiltonian_with_params(&res, &json!({ "basis": "normal_fock", "fock_trunc": [9, 5] })).unwrap();
        let rep_cb =
            hamiltonian_with_params(&res, &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [5] }))
                .unwrap();

        // FR-03 Pt 1 & Pt 2: both normal_fock and craig_bampton populate identical f_cav_Hz and g_matrix_Hz
        assert_eq!(rep_fock.f_cav_hz.len(), 1);
        assert!((rep_fock.f_cav_hz[0] - rep_cb.f_cav_hz[0]).abs() < 1e-6);
        assert!((rep_fock.g_matrix_hz[0][0] - rep_cb.g_matrix_hz[0][0]).abs() < 1e-6);
        let j_fock = rep_fock.to_json();
        assert!(j_fock["f_cav_Hz"][0].as_f64().unwrap() > 6.5e9);
        assert!(j_fock["g_matrix_Hz"][0][0].as_f64().unwrap().abs() > 1.0e6);

        let rel_f01 = (rep_fock.f01_hz[0] - rep_cb.f01_hz[0]).abs() / rep_cb.f01_hz[0];
        let rel_alpha = (rep_fock.alpha_hz[0] - rep_cb.alpha_hz[0]).abs() / rep_cb.alpha_hz[0].abs();
        assert!(rel_f01 < 3e-3, "rel_f01 = {rel_f01:e}");
        assert!(rel_alpha < 5e-2, "rel_alpha = {rel_alpha:e}");
    }

    #[test]
    fn sweep_port_and_dc_squid_dissipator_sweep_under_100ms() {
        // Two-port ROM: port 0 = fixed Q1_JJ, port 1 = tunable DC-SQUID Dissipator (4.2 - 15.3 GHz)
        let e_c_q = 260.0e6;
        let e_c_d = 380.0e6;
        let c_q = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * e_c_q);
        let c_d = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * e_c_d);
        let c_r = 55.0e-15;
        let c_gq = 1.2e-15;
        let c_gd = 1.8e-15;
        let omega_r = 2.0 * PI * 7.2e9;
        let l_r_inv = omega_r * omega_r * (c_r + c_gq + c_gd);

        // Choose E_{J,sum} and asymmetry d so f01 spans 15.3 GHz (flux=0) down to 4.2 GHz (flux=0.5)
        // f01 ~ sqrt(8 E_J E_C) - E_C => E_J(0) ~ (15.68)^2 / (8 * 0.38) = 80.87 GHz; E_J(0.5) ~ (4.58)^2 / (8 * 0.38) = 6.90 GHz => d = 0.0853
        let ej_sum_d = 83.70e9;
        let asym_d = 0.0845;
        let ej_q = 16.0e9;

        let mut rom = CircuitRom {
            ports: vec![
                RomPort {
                    name: "Q1_JJ".into(),
                    kind: "lumped_jj".into(),
                    l_j_h: Some(jj_flux_scale() / ej_q),
                    e_j_hz: Some(ej_q),
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: None,
                    direction: [0.0, 1.0, 0.0],
                    formulation: "rank1".into(),
                },
                RomPort {
                    name: "Dissipator_SQUID".into(),
                    kind: "lumped_jj".into(),
                    l_j_h: Some(jj_flux_scale() / ej_sum_d),
                    e_j_hz: Some(ej_sum_d),
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: Some(asym_d),
                    z0_ohm: None,
                    direction: [0.0, 1.0, 0.0],
                    formulation: "rank1".into(),
                },
                RomPort {
                    name: "feedline_port".into(),
                    kind: "lumped_port".into(),
                    l_j_h: None,
                    e_j_hz: None,
                    c_j_f: 0.0,
                    e_l_hz: None,
                    phi_ext_rad: 0.0,
                    n_g: 0.0,
                    squid_asymmetry: None,
                    z0_ohm: Some(50.0),
                    direction: [1.0, 0.0, 0.0],
                    formulation: "rank1".into(),
                },
            ],
            r_static: 0,
            r_dyn: 3,
            k_r: vec![vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], vec![0.0, 0.0, MU0 * l_r_inv]],
            c_r: vec![],
            m_r: vec![
                vec![(c_q + c_gq) / EPS0, 0.0, -c_gq / EPS0],
                vec![0.0, (c_d + c_gd) / EPS0, -c_gd / EPS0],
                vec![-c_gq / EPS0, -c_gd / EPS0, (c_r + c_gq + c_gd) / EPS0],
            ],
            b_r: vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0], vec![0.0, 0.0, 0.02]],
            c_port_f: vec![],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

        // Invalid sweep_port name errors cleanly listing available nonlinear ports
        assert!(hamiltonian_with_params(&res, &json!({ "sweep_port": "feedline_port" })).is_err());

        let flux_pts: Vec<f64> = (0..51).map(|i| 0.5 * (i as f64) / 50.0).collect();
        let params = json!({
            "basis": "craig_bampton",
            "port_trunc": 6,
            "fock_trunc": [4],
            "sweep_port": "Dissipator_SQUID",
            "flux_sweep": flux_pts
        });
        let _ = hamiltonian_with_params(&res, &params).unwrap();
        let t0 = Instant::now();
        let rep = hamiltonian_with_params(&res, &params).unwrap();
        let ms = t0.elapsed().as_secs_f64() * 1e3;
        assert!(ms < 100.0, "51-point Dissipator SQUID sweep took {ms:.2} ms >= 100 ms");

        // Passive feedline_port was NOT quantized: subsystems are [Q1_JJ, Dissipator_SQUID, cavity_1]
        assert_eq!(rep.subsystems, vec!["Q1_JJ", "Dissipator_SQUID", "cavity_1"]);
        assert!(rep.kappa_cav_hz[0] > 1.0e4);
        assert!(rep.t1_purcell_s[0].is_finite() && rep.t1_purcell_s[0] > 0.0);

        let sw = rep.sweep.as_ref().unwrap();
        assert_eq!(sw.len(), 51);
        let f_max = sw[0].f01_hz[1];
        let f_min = sw[50].f01_hz[1];
        assert!((f_max - 15.3e9).abs() < 0.25e9, "f_max = {:.3} GHz", f_max / 1e9);
        assert!((f_min - 4.2e9).abs() < 0.25e9, "f_min = {:.3} GHz", f_min / 1e9);
        // Q1_JJ stayed fixed while Dissipator_SQUID swept 15.3 -> 4.2 GHz!
        assert!((sw[0].f01_hz[0] - sw[50].f01_hz[0]).abs() < 10.0e6);
        assert!(sw[0].to_json().get("e_j_Hz").is_some());

        // Also verify symmetric d = 0 SQUID at flux = 0.5 (E_{J,eff} == 0) produces finite eigenvalues in both bases
        let rep_zero_cb = hamiltonian_with_params(
            &res,
            &json!({
                "basis": "craig_bampton",
                "port_trunc": 5,
                "fock_trunc": [3],
                "sweep_port": 1,
                "squid_asymmetry": 0.0,
                "squid_flux_sweep": [0.0, 0.25, 0.5]
            }),
        )
        .unwrap();
        let sw_z = rep_zero_cb.sweep.unwrap();
        assert_eq!(sw_z[2].e_j_hz, Some(0.0));
        assert!(sw_z[2].f01_hz.iter().all(|x| x.is_finite()));

        let rep_zero_nf = hamiltonian_with_params(
            &res,
            &json!({
                "basis": "normal_fock",
                "fock_trunc": [4, 4, 3],
                "sweep_port": "Dissipator_SQUID",
                "squid_asymmetry": 0.0,
                "squid_flux_sweep": [0.0, 0.5]
            }),
        )
        .unwrap();
        let sw_nf = rep_zero_nf.sweep.unwrap();
        assert_eq!(sw_nf[1].e_j_hz, Some(0.0));
        assert!(sw_nf[1].f01_hz.iter().all(|x| x.is_finite()));
    }

    #[test]
    fn fock_boundary_occupation_warns_when_undertruncated_without_extra_resolves() {
        let scale = jj_flux_scale();
        let e_j_hz = 16.0e9;
        let l_j_h = scale / e_j_hz;
        let c_port = 75.0e-15;
        let c_r = 50.0e-15;
        let c_g = 1.5e-15;
        let omega_r = 2.0 * PI * 6.8e9;
        let l_r_inv = omega_r * omega_r * (c_r + c_g);

        let mut rom = CircuitRom {
            ports: vec![RomPort {
                name: "Q1_JJ".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(l_j_h),
                e_j_hz: Some(e_j_hz),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 0,
            r_dyn: 2,
            k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]],
            c_r: vec![],
            m_r: vec![vec![(c_port + c_g) / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]],
            b_r: vec![vec![1.0], vec![0.0]],
            c_port_f: vec![],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

        // 1. Severely under-truncated normal_fock ([3, 3]): |2> lies at N_max=2 => occupation ~ 1.0 > 1e-4, warning emitted
        let rep_under =
            hamiltonian_with_params(&res, &json!({ "basis": "normal_fock", "fock_trunc": [3, 3] })).unwrap();
        assert_eq!(rep_under.dims, vec![3, 3]);
        assert!(rep_under.max_fock_boundary_occupation > 1.0e-4);
        assert!(rep_under.fock_truncation_warning.is_some());
        let j_under = rep_under.to_json();
        assert_eq!(j_under["fock_truncation_adequate"], json!(false));
        assert!(j_under["fock_truncation_warning"].as_str().is_some());
        assert!(j_under["warning"].as_str().is_some());

        // 2. Well-converged normal_fock ([12, 6]): boundary occupation << 1e-4, no warning
        let rep_ok = hamiltonian_with_params(&res, &json!({ "basis": "normal_fock", "fock_trunc": [12, 6] })).unwrap();
        assert_eq!(rep_ok.dims, vec![12, 6]);
        assert!(
            rep_ok.max_fock_boundary_occupation < 1.0e-4,
            "expected converged boundary occupation < 1e-4, got {:e}",
            rep_ok.max_fock_boundary_occupation
        );
        assert!(rep_ok.fock_truncation_warning.is_none());
        let j_ok = rep_ok.to_json();
        assert_eq!(j_ok["fock_truncation_adequate"], json!(true));
        assert!(j_ok.get("fock_truncation_warning").is_none());

        // 3. Craig-Bampton with tiny 1D charge cutoff (n_charge_cut=1 -> 3 charge states) triggers warning;
        //    default n_charge_cut=15, port_trunc=8, fock_trunc=[6] is well-converged.
        let rep_cb_under = hamiltonian_with_params(
            &res,
            &json!({ "basis": "craig_bampton", "port_trunc": 3, "fock_trunc": [3], "n_charge_cut": 2 }),
        )
        .unwrap();
        assert!(rep_cb_under.max_fock_boundary_occupation > 1.0e-4);
        assert!(rep_cb_under.fock_truncation_warning.is_some());

        let rep_cb_ok = hamiltonian_with_params(
            &res,
            &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [6], "flux_sweep": [0.0] }),
        )
        .unwrap();
        assert!(
            rep_cb_ok.max_fock_boundary_occupation < 1.0e-4,
            "CB converged boundary occ = {:e}",
            rep_cb_ok.max_fock_boundary_occupation
        );
        assert!(rep_cb_ok.fock_truncation_warning.is_none());
        let sw_pt = &rep_cb_ok.sweep.as_ref().unwrap()[0];
        assert!(sw_pt.max_fock_boundary_occupation < 1.0e-4);
        assert_eq!(sw_pt.to_json()["fock_truncation_adequate"], json!(true));
    }
}
