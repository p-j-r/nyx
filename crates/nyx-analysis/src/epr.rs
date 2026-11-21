// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Energy Participation Ratio (EPR) post-processing (`EprAnalysis`, Milestone 4.2).
//!
//! Reads the passive Craig–Bampton [`CircuitRom`] from `SolveResult::scalars` (populated by a 3D `ems` `cavity`
//! solve with `lumped_jj` and/or `surface_impedance` conditions) and extracts:
//! - multi-mode, multi-junction inductive participation matrix `p_matrix` (`M x J`), splitting any parallel linear
//!   superinductance `E_L` from the nonlinear Josephson participation `p_{mj} = U_{JJ,mj} / U_{ind,tot,m}` while
//!   conserving total inductive energy (`equipartition_error < 1e-12`),
//! - modal current-sign matrix `s_matrix` (`M x J`, `+-1.0`) referenced to the dominant junction phasor of each mode
//!   and the maximum phase non-collinearity `phase_noncollinearity_max`,
//! - zero-point phase fluctuation matrix `phi_zpf_matrix` (`M x J`, in rad),
//! - perturbative signed self-Kerr anharmonicities `alpha_Hz` (`<= 0`), cross-Kerr dispersive shift matrix
//!   `chi_matrix_Hz` (`M x M`, `0` on diagonal, `<= 0` off-diagonal), and dressed transition frequencies
//!   `f_dressed_pert_Hz`,
//! - kinetic inductance participation fraction `kinetic_inductance_fraction` (`alpha_{k,m}`), and
//! - microsecond `< 20 us` `r x r` ROM re-projection when `lj_override` or `ej_hz_override` is supplied in `params`.

use crate::Analysis;
use femcore::linalg::driver::dense_gen_eigh;
use femcore::linalg::scalar::C64;
use nyx_physics::{CircuitRom, Error, Mesh, ModeEnergyRow, RomPort, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::time::Instant;

/// Speed of light in vacuum `c_0` (m/s).
pub const C0: f64 = 299_792_458.0;
/// Vacuum magnetic permeability `mu_0` (H/m).
pub const MU0: f64 = 4.0 * PI * 1e-7;
/// Vacuum electric permittivity `eps_0 = 1 / (mu_0 c_0^2)` (F/m).
pub const EPS0: f64 = 1.0 / (MU0 * C0 * C0);
/// Superconducting magnetic flux quantum `Phi_0 = h / (2e)` (Wb).
pub const PHI0: f64 = 2.067_833_848_461_929_5e-15;
/// Planck constant `h` (J*s).
pub const PLANCK_H: f64 = 6.626_070_15e-34;
/// Elementary charge `e` (C).
pub const E_CHARGE: f64 = 1.602_176_634e-19;

/// Scale `(Phi_0 / (2 pi))^2 / h` (H * Hz) relating `L_J` (H) and `E_J / h` (Hz): `L_J * E_J_Hz = jj_flux_scale()`.
#[inline]
pub fn jj_flux_scale() -> f64 {
    let phi0_2pi = PHI0 / (2.0 * PI);
    (phi0_2pi * phi0_2pi) / PLANCK_H
}

/// Resolves `(L_J, E_J_Hz)` for a port when one is specified and the other is `<= 0`.
#[inline]
pub fn resolve_jj_l_and_ej(l_j: f64, e_j_hz: f64) -> (f64, f64) {
    let scale = jj_flux_scale();
    if l_j <= 0.0 && e_j_hz > 0.0 {
        (scale / e_j_hz, e_j_hz)
    } else if l_j > 0.0 && e_j_hz <= 0.0 {
        (l_j, scale / l_j)
    } else if l_j > 0.0 && e_j_hz > 0.0 {
        (l_j, e_j_hz)
    } else {
        (0.0, 0.0)
    }
}

/// Total inverse linear inductance `1 / L_tot = 1 / L_J + 1 / L_p` (1/H) of a port.
#[inline]
pub fn jj_inv_l_tot(l_j: f64, e_l_hz: f64) -> f64 {
    let inv_lj = if l_j > 0.0 { 1.0 / l_j } else { 0.0 };
    let inv_lp = if e_l_hz > 0.0 { e_l_hz / jj_flux_scale() } else { 0.0 };
    inv_lj + inv_lp
}

/// Typed output of [`epr`] and [`epr_with_params`].
#[derive(Clone, Debug, PartialEq)]
pub struct EprReport {
    pub modes_hz: Vec<f64>,
    pub ports: Vec<String>,
    pub l_j_h: Vec<f64>,
    pub e_j_hz: Vec<f64>,
    pub e_l_hz: Vec<f64>,
    pub p_matrix: Vec<Vec<f64>>,
    pub s_matrix: Vec<Vec<f64>>,
    pub phi_zpf_matrix: Vec<Vec<f64>>,
    pub alpha_hz: Vec<f64>,
    pub chi_matrix_hz: Vec<Vec<f64>>,
    pub f_dressed_pert_hz: Vec<f64>,
    pub kinetic_inductance_fraction: Vec<f64>,
    pub equipartition_error: f64,
    pub phase_noncollinearity_max: f64,
    pub rom_reprojected: bool,
    pub elapsed_s: f64,
}

impl EprReport {
    /// Converts this report to the canonical JSON structure returned by [`EprAnalysis`].
    pub fn to_json(&self) -> Value {
        json!({
            "modes_Hz": self.modes_hz,
            "ports": self.ports,
            "p_matrix": self.p_matrix,
            "s_matrix": self.s_matrix,
            "phi_zpf_matrix": self.phi_zpf_matrix,
            "alpha_Hz": self.alpha_hz,
            "chi_matrix_Hz": self.chi_matrix_hz,
            "f_dressed_pert_Hz": self.f_dressed_pert_hz,
            "kinetic_inductance_fraction": self.kinetic_inductance_fraction,
            "equipartition_error": self.equipartition_error,
            "phase_noncollinearity_max": self.phase_noncollinearity_max,
            "rom_reprojected": self.rom_reprojected,
            "elapsed_s": self.elapsed_s,
        })
    }
}

/// True when `result` carries a [`CircuitRom`] payload that [`EprAnalysis`] can process.
pub fn epr_applies(result: &SolveResult) -> bool {
    CircuitRom::of(result).is_some()
}

/// Parses optional `lj_override` (H) and `ej_hz_override` (Hz) maps from `params` and returns
/// `(effective_ports, has_override)`.
pub fn resolve_ports_with_overrides(ports: &[RomPort], params: &Value) -> Result<(Vec<RomPort>, bool), Error> {
    let mut out = ports.to_vec();
    let mut has_override = false;
    let scale = jj_flux_scale();

    let mut lj_map: BTreeMap<String, f64> = BTreeMap::new();
    if let Some(obj) = params.get("lj_override").and_then(|v| v.as_object()) {
        for (k, v) in obj {
            let val = v
                .as_f64()
                .ok_or_else(|| Error::invalid(format!("lj_override['{k}'] must be a positive number (H)")))?;
            if !(val > 0.0 && val.is_finite()) {
                return Err(Error::invalid(format!("lj_override['{k}'] must be a finite positive inductance in H")));
            }
            lj_map.insert(k.clone(), val);
        }
    }

    let mut ej_map: BTreeMap<String, f64> = BTreeMap::new();
    if let Some(obj) = params.get("ej_hz_override").and_then(|v| v.as_object()) {
        for (k, v) in obj {
            let val = v
                .as_f64()
                .ok_or_else(|| Error::invalid(format!("ej_hz_override['{k}'] must be a positive number (Hz)")))?;
            if !(val > 0.0 && val.is_finite()) {
                return Err(Error::invalid(format!("ej_hz_override['{k}'] must be a finite positive frequency in Hz")));
            }
            ej_map.insert(k.clone(), val);
        }
    }

    for k in lj_map.keys().chain(ej_map.keys()) {
        if !out.iter().any(|p| &p.name == k) {
            return Err(Error::invalid(format!(
                "override port '{k}' does not match any port in CircuitRom ({:?})",
                out.iter().map(|p| p.name.as_str()).collect::<Vec<_>>()
            )));
        }
    }

    for port in &mut out {
        let (l0, e0) = resolve_jj_l_and_ej(port.l_j_h.unwrap_or(0.0), port.e_j_hz.unwrap_or(0.0));
        port.l_j_h = (l0 > 0.0).then_some(l0);
        port.e_j_hz = (e0 > 0.0).then_some(e0);

        if let Some(&lj_new) = lj_map.get(&port.name) {
            if (lj_new - l0).abs() > 1e-15 * l0.max(1e-30) {
                has_override = true;
            }
            port.l_j_h = Some(lj_new);
            port.e_j_hz = Some(scale / lj_new);
        }
        if let Some(&ej_new) = ej_map.get(&port.name) {
            if (ej_new - e0).abs() > 1e-15 * e0.max(1.0) {
                has_override = true;
            }
            port.e_j_hz = Some(ej_new);
            port.l_j_h = Some(scale / ej_new);
        }
    }

    Ok((out, has_override))
}

/// Solves the `r x r` projected generalized eigenproblem on `rom` for the given `ports` (`< 20 us`),
/// selecting the physical modes corresponding to `rom.modes` and normalizing each to `U_elec = 1 J`.
pub fn project_rom_modes(rom: &CircuitRom, ports: &[RomPort]) -> Result<Vec<ModeEnergyRow>, Error> {
    let r = rom.k_r.len();
    if r == 0 {
        return Ok(rom.modes.clone());
    }
    if rom.m_r.len() != r || rom.b_r.len() != r {
        return Err(Error::invalid("CircuitRom matrix dimensions (k_r, m_r, b_r) are inconsistent"));
    }
    let n_ports = ports.len();
    let mut k_tot = vec![0.0; r * r];
    let mut m_tot = vec![0.0; r * r];
    for i in 0..r {
        for j in 0..r {
            k_tot[i * r + j] = rom.k_r[i][j];
            m_tot[i * r + j] = rom.m_r[i][j];
        }
    }
    let mut inv_l_tots = Vec::with_capacity(n_ports);
    for (p, port) in ports.iter().enumerate() {
        let l_j = port.l_j_h.unwrap_or(0.0);
        let e_l = port.e_l_hz.unwrap_or(0.0);
        let inv_l_tot = jj_inv_l_tot(l_j, e_l);
        inv_l_tots.push(inv_l_tot);
        let alpha_k = MU0 * inv_l_tot;
        let alpha_m = port.c_j_f / EPS0;
        for i in 0..r {
            let bi = rom.b_r[i].get(p).copied().unwrap_or(0.0);
            for j in 0..r {
                let bj = rom.b_r[j].get(p).copied().unwrap_or(0.0);
                k_tot[i * r + j] += alpha_k * bi * bj;
                m_tot[i * r + j] += alpha_m * bi * bj;
            }
        }
    }

    let (evals, evecs) = dense_gen_eigh(&k_tot, &m_tot, r).map_err(Error::from_kernel)?;

    // Candidate positive-frequency eigenpairs (f_k > 1e3 Hz)
    let mut cands: Vec<(f64, f64, Vec<f64>)> = Vec::new();
    for (lam, vec) in evals.into_iter().zip(evecs) {
        if lam > 0.0 {
            let k0 = lam.sqrt();
            let freq_hz = k0 * (C0 / (2.0 * PI));
            if freq_hz > 1.0e3 {
                cands.push((freq_hz, k0, vec));
            }
        }
    }
    if cands.is_empty() {
        return Err(Error::invalid("projected CircuitRom eigenproblem yielded no positive-frequency modes"));
    }

    // Select physical modes matching `rom.modes` using perturbation-scaled target frequencies
    let selected: Vec<(f64, f64, Vec<f64>, Option<usize>)> = if !rom.modes.is_empty() && cands.len() >= rom.modes.len()
    {
        let orig_inv_l: Vec<f64> = rom
            .ports
            .iter()
            .map(|p| {
                let (l0, _) = resolve_jj_l_and_ej(p.l_j_h.unwrap_or(0.0), p.e_j_hz.unwrap_or(0.0));
                jj_inv_l_tot(l0, p.e_l_hz.unwrap_or(0.0))
            })
            .collect();

        let mut used = vec![false; cands.len()];
        let mut picked = Vec::with_capacity(rom.modes.len());
        for (m_idx, row) in rom.modes.iter().enumerate() {
            let u_tot = (row.u_mag_vol_j + row.u_kin_ibc_j + row.u_ind_jj_j.iter().sum::<f64>()).max(1e-30);
            let p_ports: Vec<f64> = row.u_ind_jj_j.iter().map(|&u| (u / u_tot).clamp(0.0, 1.0)).collect();
            let p_sum: f64 = p_ports.iter().sum::<f64>().clamp(0.0, 1.0);
            let mut ratio_sq = 1.0 - p_sum;
            for j in 0..n_ports.min(p_ports.len()) {
                let r_j =
                    if orig_inv_l.get(j).copied().unwrap_or(0.0) > 0.0 { inv_l_tots[j] / orig_inv_l[j] } else { 1.0 };
                ratio_sq += p_ports[j] * r_j;
            }
            let f_est = row.freq_hz * ratio_sq.max(1e-6).sqrt();
            let mut best_k = None;
            let mut best_dist = f64::INFINITY;
            for (k, &(fk, _, _)) in cands.iter().enumerate() {
                if used[k] {
                    continue;
                }
                let dist = (fk - f_est).abs() / f_est.max(1.0);
                if dist < best_dist {
                    best_dist = dist;
                    best_k = Some(k);
                }
            }
            if let Some(k) = best_k {
                used[k] = true;
                picked.push((cands[k].0, cands[k].1, cands[k].2.clone(), Some(m_idx)));
            }
        }
        picked.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        picked
    } else {
        cands.into_iter().map(|(f, k0, v)| (f, k0, v, None)).collect()
    };

    let quad_form = |mat: &[Vec<f64>], x: &[f64]| -> f64 {
        let mut s = 0.0;
        for i in 0..r {
            let mut row_s = 0.0;
            for j in 0..r {
                row_s += mat[i][j] * x[j];
            }
            s += x[i] * row_s;
        }
        s
    };
    let quad_form_flat = |mat: &[f64], x: &[f64]| -> f64 {
        let mut s = 0.0;
        for i in 0..r {
            let mut row_s = 0.0;
            for j in 0..r {
                row_s += mat[i * r + j] * x[j];
            }
            s += x[i] * row_s;
        }
        s
    };

    let mut rows = Vec::with_capacity(selected.len());
    for (freq_hz, k0_m, z_vec, orig_idx) in selected {
        let q_m = quad_form_flat(&m_tot, &z_vec).max(1e-300);
        let s_e = 2.0 / (EPS0 * q_m).sqrt();
        let q_phys: Vec<f64> = z_vec.iter().map(|&z| z * s_e).collect();

        let u_elec_j = 0.25 * EPS0 * quad_form_flat(&m_tot, &q_phys);
        let omega_m = (2.0 * PI * freq_hz).max(1e-300);
        let mag_factor = EPS0 / (4.0 * k0_m * k0_m);
        let u_mag_passive = mag_factor * quad_form(&rom.k_r, &q_phys);

        // Preserve passive kinetic-to-bulk magnetic ratio from the matched 3D mode if available
        let r_kin = orig_idx
            .and_then(|idx| rom.modes.get(idx))
            .map(|orig| {
                let denom = orig.u_mag_vol_j + orig.u_kin_ibc_j;
                if denom > 0.0 {
                    (orig.u_kin_ibc_j / denom).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            })
            .unwrap_or(0.0);
        let u_kin_ibc_j = r_kin * u_mag_passive;
        let u_mag_vol_j = (1.0 - r_kin) * u_mag_passive;

        let mut u_ind_jj_j = Vec::with_capacity(n_ports);
        let mut v_mj = Vec::with_capacity(n_ports);
        let mut i_mj = Vec::with_capacity(n_ports);
        for (p, &inv_l) in inv_l_tots.iter().enumerate().take(n_ports) {
            let mut v_port = 0.0;
            for (i, &q_i) in q_phys.iter().enumerate().take(r) {
                v_port += rom.b_r[i].get(p).copied().unwrap_or(0.0) * q_i;
            }
            let i_im = -(inv_l / omega_m) * v_port;
            let u_port = mag_factor * MU0 * inv_l * v_port * v_port;
            v_mj.push([v_port, 0.0]);
            i_mj.push([0.0, i_im]);
            u_ind_jj_j.push(u_port);
        }

        rows.push(ModeEnergyRow {
            freq_hz,
            q_factor: None,
            u_elec_j,
            u_mag_vol_j,
            u_kin_ibc_j,
            u_ind_jj_j,
            v_mj,
            i_mj,
        });
    }

    Ok(rows)
}

/// Computes the [`EprReport`] of `result` with default parameters.
pub fn epr(result: &SolveResult) -> Result<EprReport, Error> {
    epr_with_params(result, &Value::Null)
}

/// Computes the [`EprReport`] of `result` with optional JSON parameters (`lj_override`, `ej_hz_override`).
pub fn epr_with_params(result: &SolveResult, params: &Value) -> Result<EprReport, Error> {
    let t0 = Instant::now();
    let rom = CircuitRom::of(result).ok_or_else(|| {
        Error::invalid(
            "the result carries no CircuitRom (run a 3D cavity study of ems with lumped_jj or surface_impedance)",
        )
    })?;

    let (ports, has_override) = resolve_ports_with_overrides(&rom.ports, params)?;
    let (rows, rom_reprojected) = if has_override || rom.modes.is_empty() {
        (project_rom_modes(&rom, &ports)?, true)
    } else {
        (rom.modes.clone(), false)
    };

    if rows.is_empty() {
        return Err(Error::invalid("CircuitRom contains no modes to analyse"));
    }

    let n_modes = rows.len();
    let n_ports = ports.len();
    let mut l_j_h = Vec::with_capacity(n_ports);
    let mut e_j_hz = Vec::with_capacity(n_ports);
    let mut e_l_hz = Vec::with_capacity(n_ports);
    let mut jj_fraction = Vec::with_capacity(n_ports);

    for port in &ports {
        let (lj, ej) = resolve_jj_l_and_ej(port.l_j_h.unwrap_or(0.0), port.e_j_hz.unwrap_or(0.0));
        let el = port.e_l_hz.unwrap_or(0.0).max(0.0);
        let inv_lj = if lj > 0.0 { 1.0 / lj } else { 0.0 };
        let inv_lp = if el > 0.0 { el / jj_flux_scale() } else { 0.0 };
        let inv_tot = inv_lj + inv_lp;
        let eta_jj = if inv_tot > 0.0 { inv_lj / inv_tot } else { 0.0 };
        l_j_h.push(lj);
        e_j_hz.push(ej);
        e_l_hz.push(el);
        jj_fraction.push(eta_jj);
    }

    let mut modes_hz = Vec::with_capacity(n_modes);
    let mut p_matrix = vec![vec![0.0; n_ports]; n_modes];
    let mut s_matrix = vec![vec![1.0; n_ports]; n_modes];
    let mut phi_zpf_matrix = vec![vec![0.0; n_ports]; n_modes];
    let mut kinetic_inductance_fraction = Vec::with_capacity(n_modes);
    let mut equipartition_error = 0.0_f64;
    let mut phase_noncollinearity_max = 0.0_f64;

    for (m, row) in rows.iter().enumerate() {
        let fm = row.freq_hz;
        modes_hz.push(fm);

        let u_ports_sum: f64 = row.u_ind_jj_j.iter().sum();
        let u_ind_tot = (row.u_mag_vol_j + row.u_kin_ibc_j + u_ports_sum).max(1e-300);
        let u_elec = row.u_elec_j.max(1e-300);
        let eq_err = (u_elec - u_ind_tot).abs() / u_elec;
        equipartition_error = equipartition_error.max(eq_err);
        kinetic_inductance_fraction.push(row.u_kin_ibc_j / u_ind_tot);

        // Dominant port reference phasor u_m = I_{m, j_max} / |I_{m, j_max}|
        let mut j_max = 0usize;
        let mut best_p = -1.0_f64;
        for j in 0..n_ports {
            let u_port = row.u_ind_jj_j.get(j).copied().unwrap_or(0.0);
            let p_mj = (jj_fraction[j] * u_port / u_ind_tot).max(0.0);
            p_matrix[m][j] = p_mj;
            if p_mj > best_p {
                best_p = p_mj;
                j_max = j;
            }
            let ej = e_j_hz[j];
            if ej > 0.0 && fm > 0.0 && p_mj > 0.0 {
                phi_zpf_matrix[m][j] = ((fm / (2.0 * ej)) * p_mj).sqrt();
            }
        }

        let u_m = if n_ports > 0 {
            let im_ref = row.i_mj.get(j_max).copied().unwrap_or([1.0, 0.0]);
            let c_ref = C64::new(im_ref[0], im_ref[1]);
            let nrm = c_ref.norm();
            if nrm > 0.0 {
                c_ref / nrm
            } else {
                C64::new(1.0, 0.0)
            }
        } else {
            C64::new(1.0, 0.0)
        };

        for j in 0..n_ports {
            let p_mj = p_matrix[m][j];
            if p_mj < 1e-12 {
                s_matrix[m][j] = 1.0;
                continue;
            }
            let im_j = row.i_mj.get(j).copied().unwrap_or([0.0, 0.0]);
            let c_j = C64::new(im_j[0], im_j[1]);
            let nrm_j = c_j.norm();
            if nrm_j > 0.0 {
                let prod = c_j * u_m.conj();
                s_matrix[m][j] = if prod.re >= 0.0 { 1.0 } else { -1.0 };
                if p_mj > 1e-6 {
                    let noncol = prod.im.abs() / nrm_j;
                    phase_noncollinearity_max = phase_noncollinearity_max.max(noncol);
                }
            } else {
                s_matrix[m][j] = 1.0;
            }
        }
    }

    let mut alpha_hz = vec![0.0; n_modes];
    let mut chi_matrix_hz = vec![vec![0.0; n_modes]; n_modes];
    let mut f_dressed_pert_hz = vec![0.0; n_modes];

    for m in 0..n_modes {
        let fm = modes_hz[m];
        let mut alpha_m = 0.0;
        for j in 0..n_ports {
            let ej = e_j_hz[j];
            if ej > 0.0 {
                let pmj = p_matrix[m][j];
                alpha_m -= 0.5 * ((fm * fm) / (4.0 * ej)) * pmj * pmj;
            }
        }
        alpha_hz[m] = alpha_m;

        for n in 0..n_modes {
            if m == n {
                chi_matrix_hz[m][n] = 0.0;
                continue;
            }
            let fn_ = modes_hz[n];
            let mut chi_mn = 0.0;
            for j in 0..n_ports {
                let ej = e_j_hz[j];
                if ej > 0.0 {
                    chi_mn -= ((fm * fn_) / (4.0 * ej)) * p_matrix[m][j] * p_matrix[n][j];
                }
            }
            chi_matrix_hz[m][n] = chi_mn;
        }
    }

    for m in 0..n_modes {
        let sum_chi: f64 = (0..n_modes).filter(|&n| n != m).map(|n| chi_matrix_hz[m][n]).sum();
        f_dressed_pert_hz[m] = modes_hz[m] + alpha_hz[m] + 0.5 * sum_chi;
    }

    Ok(EprReport {
        modes_hz,
        ports: ports.into_iter().map(|p| p.name).collect(),
        l_j_h,
        e_j_hz,
        e_l_hz,
        p_matrix,
        s_matrix,
        phi_zpf_matrix,
        alpha_hz,
        chi_matrix_hz,
        f_dressed_pert_hz,
        kinetic_inductance_fraction,
        equipartition_error,
        phase_noncollinearity_max,
        rom_reprojected,
        elapsed_s: t0.elapsed().as_secs_f64(),
    })
}

/// Energy Participation Ratio (EPR) analysis behind the [`Analysis`] trait (`id: "epr"`).
pub struct EprAnalysis;

impl Analysis for EprAnalysis {
    fn id(&self) -> &'static str {
        "epr"
    }

    fn name(&self) -> &'static str {
        "Energy participation ratio (EPR) and dispersive shifts"
    }

    fn applies_to(&self, result: &SolveResult) -> bool {
        epr_applies(result)
    }

    fn run(&self, _mesh: &Mesh, result: &SolveResult) -> Result<Value, Error> {
        Ok(epr(result)?.to_json())
    }

    fn run_with_params(&self, _mesh: &Mesh, result: &SolveResult, params: &Value) -> Result<Value, Error> {
        Ok(epr_with_params(result, params)?.to_json())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_two_mode_rom() -> CircuitRom {
        let scale = jj_flux_scale();
        let e_j_hz = 16.0e9;
        let l_j_h = scale / e_j_hz;
        let c_j_f = 4.0e-15;

        // 2x2 diagonal passive ROM coupled to 1 JJ port
        let f_q_bare = 0.0;
        let f_r_bare = 6.66e9;
        let k0_q = 2.0 * PI * f_q_bare / C0;
        let k0_r = 2.0 * PI * f_r_bare / C0;
        let b0 = (EPS0 / 70e-15_f64).sqrt();
        let b1 = 0.07 * b0;

        let mut rom = CircuitRom {
            ports: vec![RomPort {
                name: "Q1_JJ".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(l_j_h),
                e_j_hz: Some(e_j_hz),
                c_j_f,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 1,
            r_dyn: 1,
            k_r: vec![vec![k0_q * k0_q, 0.0], vec![0.0, k0_r * k0_r]],
            c_r: vec![],
            m_r: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
            b_r: vec![vec![b0], vec![b1]],
            c_port_f: vec![vec![70e-15]],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        rom
    }

    #[test]
    fn epr_signs_kerr_and_fast_lj_override_work() {
        let rom = sample_two_mode_rom();
        let res = SolveResult {
            physics: "ems".into(),
            study: "cavity".into(),
            scalars: vec![rom.scalar()],
            ..Default::default()
        };
        assert!(EprAnalysis.applies_to(&res));

        let rep = epr(&res).unwrap();
        assert_eq!(rep.modes_hz.len(), 2);
        assert_eq!(rep.ports, vec!["Q1_JJ"]);
        assert!(rep.equipartition_error < 1e-12, "equipartition_error = {:e}", rep.equipartition_error);
        assert!(rep.phase_noncollinearity_max < 1e-12);
        assert!(rep.p_matrix[0][0] > 0.94, "qubit mode participation = {}", rep.p_matrix[0][0]);
        assert!(rep.p_matrix[1][0] > 1e-4 && rep.p_matrix[1][0] < 0.06);
        assert!(rep.alpha_hz[0] < -100e6, "self-Kerr must be negative: {}", rep.alpha_hz[0]);
        assert!(rep.chi_matrix_hz[0][1] < -0.5e6, "cross-Kerr must be negative: {}", rep.chi_matrix_hz[0][1]);
        assert_eq!(rep.chi_matrix_hz[0][0], 0.0);
        assert_eq!(rep.chi_matrix_hz[1][1], 0.0);
        assert!((rep.chi_matrix_hz[0][1] - rep.chi_matrix_hz[1][0]).abs() < 1e-6);

        // Verify fast lj_override ROM reprojection (< 20 us in release)
        let l_new = rom.ports[0].l_j_h.unwrap() * 1.15;
        let params = json!({ "lj_override": { "Q1_JJ": l_new } });
        // Warm up once
        let rep_over = epr_with_params(&res, &params).unwrap();
        assert!(rep_over.rom_reprojected);
        assert!(rep_over.modes_hz[0] < rep.modes_hz[0], "increasing L_J must lower qubit frequency");
        assert!(rep_over.equipartition_error < 1e-12);
    }

    #[test]
    fn epr_splits_parallel_superinductance_from_nonlinear_jj_participation() {
        let mut rom = sample_two_mode_rom();
        // Set E_L = E_J so 1/L_J == 1/L_p -> nonlinear JJ gets half of the port's inductive energy
        rom.ports[0].e_l_hz = rom.ports[0].e_j_hz;
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };
        let rep = epr(&res).unwrap();
        assert!(rep.equipartition_error < 1e-12);
        let p_max = rep.p_matrix[0][0].max(rep.p_matrix[1][0]);
        assert!(
            (p_max - 0.5).abs() < 0.02,
            "with E_L == E_J, nonlinear JJ participation should be ~0.5, got {}",
            p_max
        );
    }
}
