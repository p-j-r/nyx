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
//!   [`femcore::linalg::dense::hermitian_eig`].
//! - **Mode B (`basis = "craig_bampton"`)**: Nonlinear-port + linear-cavity Craig–Bampton basis for strongly
//!   anharmonic / multi-well circuits (Fluxonium, Cooper-pair box, zero-pi, tunable couplers) coupled to 3D cavities.
//!   Decomposes `R^r` into port lift `T_p` (`B_r^T T_p = I_J`) and short-circuit interior cavity modes `T_c`
//!   (`B_r^T T_c = 0`), pre-diagonalizes each 1D nonlinear port (`E_{L,eff} > 1e3 Hz` in a 40-state harmonic
//!   oscillator basis with cached `expm(i phi_j)`; `E_{L,eff} <= 1e3 Hz` in a 31-state charge-periodic basis),
//!   projects `(n_j, phi_j)` onto the lowest `port_trunc` (default 8) port eigenstates, and diagonalizes the coupled
//!   `H / h` in `< 100 us` per flux point (`< 10 ms` for a 100-point flux sweep).
//! - **Globally Bijective Greedy Descending Bipartite State Labeling**: Assigns bare product states in
//!   `T = {|0...0>, |1_a>, |2_a>, |1_a 1_b>}` to eigenstates bijectively in descending order of overlap
//!   `|<n | psi_k>|^2` to extract signed `f01_Hz`, `f12_Hz`, `f02_Hz`, `alpha_Hz`, and `chi_matrix_Hz`.

use crate::epr::{
    epr_with_params, jj_flux_scale, resolve_jj_l_and_ej, resolve_ports_with_overrides, EPS0, E_CHARGE, MU0, PLANCK_H,
};
use crate::Analysis;
use femcore::linalg::dense::{hermitian_eig, Dense};
use femcore::linalg::driver::{dense_eigh, dense_gen_eigh};
use femcore::linalg::expm::expm;
use femcore::linalg::scalar::{Scalar, C64};
use nyx_physics::{CircuitRom, Error, Mesh, RomPort, SolveResult};
use serde_json::{json, Value};
use std::f64::consts::PI;
use std::time::Instant;

/// One point of a flux sweep (`phi_ext_sweep` or `flux_sweep`).
#[derive(Clone, Debug, PartialEq)]
pub struct FluxSweepPoint {
    pub phi_ext_rad: f64,
    pub flux_phi0: f64,
    pub f01_hz: Vec<f64>,
    pub f12_hz: Vec<f64>,
    pub f02_hz: Vec<f64>,
    pub alpha_hz: Vec<f64>,
    pub chi_matrix_hz: Vec<Vec<f64>>,
    pub eigenvalues_hz: Vec<f64>,
    pub min_overlap: f64,
}

impl FluxSweepPoint {
    pub fn to_json(&self) -> Value {
        json!({
            "phi_ext_rad": self.phi_ext_rad,
            "flux_phi0": self.flux_phi0,
            "f01_Hz": self.f01_hz,
            "f12_Hz": self.f12_hz,
            "f02_Hz": self.f02_hz,
            "alpha_Hz": self.alpha_hz,
            "chi_matrix_Hz": self.chi_matrix_hz,
            "eigenvalues_Hz": self.eigenvalues_hz,
            "min_overlap": self.min_overlap,
        })
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
    pub elapsed_s: f64,
    pub f_cav_hz: Vec<f64>,
    pub g_matrix_hz: Vec<Vec<f64>>,
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
            "elapsed_s": self.elapsed_s,
        });
        if let Some(ref sw) = self.sweep {
            obj["sweep"] = json!(sw.iter().map(|p| p.to_json()).collect::<Vec<_>>());
        }
        if let Some(t) = self.sweep_elapsed_s {
            obj["sweep_elapsed_s"] = json!(t);
        }
        obj
    }
}

/// True when `result` carries a [`CircuitRom`] with at least one nonlinear port.
pub fn hamiltonian_applies(result: &SolveResult) -> bool {
    CircuitRom::of(result).is_some_and(|r| !r.ports.is_empty())
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
    for &(k, p) in matched_eig.iter().flatten() {
        let _ = k;
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

    LabeledSpectrum { f01_hz, f12_hz, f02_hz, alpha_hz, chi_matrix_hz, eigenvalues_hz, min_overlap: min_ov }
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
        let last = *dims.last().unwrap();
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

/// Parses optional `phi_ext_rad` (or `flux_phi0`) and `n_g` overrides onto `ports`.
fn apply_bias_overrides(ports: &mut [RomPort], params: &Value) {
    if let Some(v) = params.get("phi_ext_rad") {
        if let Some(arr) = v.as_array() {
            for (p, val) in ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.phi_ext_rad = f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in ports.iter_mut() {
                p.phi_ext_rad = f;
            }
        }
    } else if let Some(v) = params.get("flux_phi0") {
        if let Some(arr) = v.as_array() {
            for (p, val) in ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.phi_ext_rad = 2.0 * PI * f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in ports.iter_mut() {
                p.phi_ext_rad = 2.0 * PI * f;
            }
        }
    }

    if let Some(v) = params.get("n_g") {
        if let Some(arr) = v.as_array() {
            for (p, val) in ports.iter_mut().zip(arr) {
                if let Some(f) = val.as_f64() {
                    p.n_g = f;
                }
            }
        } else if let Some(f) = v.as_f64() {
            for p in ports.iter_mut() {
                p.n_g = f;
            }
        }
    }
}

/// Parses an optional flux sweep array from `params` (`phi_ext_sweep` in rad, or `flux_sweep` in `Phi_0`).
fn parse_flux_sweep(params: &Value) -> Option<Vec<f64>> {
    if let Some(arr) = params.get("phi_ext_sweep").and_then(|v| v.as_array()) {
        let pts: Vec<f64> = arr.iter().filter_map(|v| v.as_f64()).collect();
        if !pts.is_empty() {
            return Some(pts);
        }
    }
    if let Some(arr) = params.get("flux_sweep").and_then(|v| v.as_array()) {
        let pts: Vec<f64> = arr.iter().filter_map(|v| v.as_f64()).map(|f| 2.0 * PI * f).collect();
        if !pts.is_empty() {
            return Some(pts);
        }
    }
    None
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

fn run_normal_fock(result: &SolveResult, params: &Value) -> Result<HamiltonianReport, Error> {
    let t0 = Instant::now();
    let rom = CircuitRom::of(result).ok_or_else(|| Error::invalid("the result carries no CircuitRom"))?;
    let epr_rep = epr_with_params(result, params)?;
    let n_modes = epr_rep.modes_hz.len();
    if n_modes == 0 {
        return Err(Error::invalid("CircuitRom has no modes for normal_fock quantization"));
    }

    let (mut ports, _) = resolve_ports_with_overrides(&rom.ports, params)?;
    apply_bias_overrides(&mut ports, params);
    let n_ports = ports.len();

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

    // Precompute single-mode operators phi_{mj}, phi_{mj}^2, and U_{mj} = expm(i phi_{mj})
    let mut phi_ops: Vec<Vec<Dense<C64>>> = Vec::with_capacity(n_ports);
    let mut phi2_ops: Vec<Vec<Dense<C64>>> = Vec::with_capacity(n_ports);
    let mut u_full_ports: Vec<Dense<C64>> = Vec::with_capacity(n_ports);

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
        u_full_ports.push(kron_all(&u_j, &dims, &strides, d_tot));
        phi_ops.push(p_j);
        phi2_ops.push(p2_j);
    }

    let solve_at_biases = |phi_exts: &[f64]| -> Result<LabeledSpectrum, Error> {
        let mut h_mat = Dense::<C64>::zeros(d_tot, d_tot);
        // Bare linear harmonic mode energies sum_m f_m a_m^dagger a_m
        for idx in 0..d_tot {
            let mut diag_e = 0.0;
            for m in 0..n_modes {
                let nm = (idx / strides[m]) % dims[m];
                diag_e += epr_rep.modes_hz[m] * (nm as f64);
            }
            h_mat[(idx, idx)] = C64::new(diag_e, 0.0);
        }

        // Add nonlinear cosine correction for each port j:
        // - [ (E_J / 2)(e^{-i phi_ext} U_j + e^{+i phi_ext} U_j^dagger) - E_J cos(phi_ext) I + (1/2) E_J cos(phi_ext) phi_j^2 ]
        for j in 0..n_ports {
            let ej = epr_rep.e_j_hz[j];
            if ej <= 0.0 {
                continue;
            }
            let phi_ext = phi_exts[j];
            let cos_phi = phi_ext.cos();
            let phase = C64::new(cos_phi, -phi_ext.sin());
            let half_ej = 0.5 * ej;
            let u_j = &u_full_ports[j];

            for c in 0..d_tot {
                for r in 0..d_tot {
                    let u_rc = phase * u_j[(r, c)];
                    let u_cr_conj = (phase * u_j[(c, r)]).conj();
                    h_mat[(r, c)] -= (u_rc + u_cr_conj).scale(half_ej);
                }
                h_mat[(c, c)] += C64::new(ej * cos_phi, 0.0);
            }

            let quad_scale = -0.5 * ej * cos_phi;
            for (m, phi2_jm) in phi2_ops[j].iter().enumerate().take(n_modes) {
                add_single_mode_op(&mut h_mat, &dims, &strides, m, phi2_jm, quad_scale);
            }
            for m in 0..n_modes {
                for n in (m + 1)..n_modes {
                    add_two_mode_op(
                        &mut h_mat,
                        &dims,
                        &strides,
                        m,
                        &phi_ops[j][m],
                        n,
                        &phi_ops[j][n],
                        2.0 * quad_scale,
                    );
                }
            }
        }

        let (evals, evecs) = hermitian_eig(&h_mat).map_err(|e| Error::from_kernel(e.to_string()))?;
        Ok(label_spectrum_bipartite(&dims, &evals, |b_idx, k| evecs[(b_idx, k)].norm_sqr()))
    };

    let base_phi_exts: Vec<f64> = ports.iter().map(|p| p.phi_ext_rad).collect();
    let spec = solve_at_biases(&base_phi_exts)?;
    let elapsed_s = t0.elapsed().as_secs_f64();

    let (sweep, sweep_elapsed_s) = if let Some(sweep_pts) = parse_flux_sweep(params) {
        let t_sw = Instant::now();
        let mut out_pts = Vec::with_capacity(sweep_pts.len());
        let mut phis = base_phi_exts.clone();
        for phi in sweep_pts {
            phis[0] = phi;
            let sp = solve_at_biases(&phis)?;
            out_pts.push(FluxSweepPoint {
                phi_ext_rad: phi,
                flux_phi0: phi / (2.0 * PI),
                f01_hz: sp.f01_hz,
                f12_hz: sp.f12_hz,
                f02_hz: sp.f02_hz,
                alpha_hz: sp.alpha_hz,
                chi_matrix_hz: sp.chi_matrix_hz,
                eigenvalues_hz: sp.eigenvalues_hz,
                min_overlap: sp.min_overlap,
            });
        }
        (Some(out_pts), Some(t_sw.elapsed().as_secs_f64()))
    } else {
        (None, None)
    };

    let subsystems = (0..n_modes).map(|m| format!("mode_{}", m + 1)).collect();
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
        elapsed_s,
        f_cav_hz: Vec::new(),
        g_matrix_hz: Vec::new(),
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

    // 2. Implicit Wilkinson-shifted QL iteration on tridiagonal (d, e) accumulating eigenvectors in rows of `z_t`
    // Store eigenvectors in column-major layout `vecs_col[col * n + row]` for contiguous inner loop updates
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

fn dot_vec(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

fn mat_vec(m: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    m.iter().map(|row| dot_vec(row, x)).collect()
}

/// Precomputed 1D nonlinear port model for fast microsecond flux sweeps.
enum PreparedPort1D {
    /// Fluxonium / inductively-shunted port (`E_{L,eff} > 1e3 Hz`):
    /// real symmetric harmonic oscillator basis of size `n_flux` with cached `cos(phi)` and `sin(phi)` from `expm(i phi)`.
    Fluxonium {
        n_flux: usize,
        e_j_hz: f64,
        h0_diag: Vec<f64>,
        cos_phi: Vec<f64>,
        sin_phi: Vec<f64>,
        /// Real skew-symmetric matrix `y = i n = (b - b^dagger) / (2 phi_zpf)` (`n_flux x n_flux`).
        y_op: Vec<f64>,
    },
    /// Charge-periodic port (`E_{L,eff} <= 1e3 Hz`):
    /// tridiagonal charge basis `|-N_c..+N_c>` of size `2 N_c + 1`.
    ChargePeriodic { n_states: usize, n_cut: isize, e_c_hz: f64, e_j_hz: f64, n_g: f64 },
}

/// Solves a 1D port at `(phi_ext_rad, n_g_override)` and returns its lowest `d_j` eigenvalues (Hz)
/// and the projected `d_j x d_j` operator:
/// - for `Fluxonium`: real symmetric projected `y_proj = i n_proj` (`is_y_real = true`),
/// - for `ChargePeriodic`: complex Hermitian projected `n_proj` (`is_y_real = false`).
struct SolvedPort1D {
    evals_hz: Vec<f64>,
    y_real: Option<Vec<f64>>,
    n_cplx: Option<Dense<C64>>,
}

impl PreparedPort1D {
    fn solve(&self, phi_ext_rad: f64, d_j: usize) -> Result<SolvedPort1D, Error> {
        match self {
            PreparedPort1D::Fluxonium { n_flux, e_j_hz, h0_diag, cos_phi, sin_phi, y_op } => {
                let n = *n_flux;
                let d = d_j.min(n);
                let c_ext = phi_ext_rad.cos();
                let s_ext = phi_ext_rad.sin();
                let ej = *e_j_hz;
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
                // Project real skew-symmetric tridiagonal y_op onto lowest d eigenvectors: y_proj = V^T y_op V
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
                Ok(SolvedPort1D { evals_hz: evals.into_iter().take(d).collect(), y_real: Some(y_proj), n_cplx: None })
            }
            PreparedPort1D::ChargePeriodic { n_states, n_cut, e_c_hz, e_j_hz, n_g } => {
                let n = *n_states;
                let d = d_j.min(n);
                let mut h_mat = Dense::<C64>::zeros(n, n);
                let half_ej = -0.5 * (*e_j_hz);
                let off_upper = C64::new(half_ej * phi_ext_rad.cos(), -half_ej * phi_ext_rad.sin());
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
                Ok(SolvedPort1D { evals_hz: evals.into_iter().take(d).collect(), y_real: None, n_cplx: Some(n_proj) })
            }
        }
    }
}

/// Builds the Craig-Bampton decomposition `(T_p, T_c)` and inverse capacitance matrix `C_tilde^-1`,
/// returning `(subsystems, dims, port_models, f_cav_dressed_hz, g_jm_hz, ec_jk_hz, g_mn_hz, all_fluxonium)`.
fn run_craig_bampton(result: &SolveResult, params: &Value) -> Result<HamiltonianReport, Error> {
    let t0 = Instant::now();
    let rom = CircuitRom::of(result).ok_or_else(|| Error::invalid("the result carries no CircuitRom"))?;
    if rom.ports.is_empty() {
        return Err(Error::invalid("craig_bampton quantization requires at least one nonlinear port in CircuitRom"));
    }
    let (mut ports, _) = resolve_ports_with_overrides(&rom.ports, params)?;
    apply_bias_overrides(&mut ports, params);

    let n_ports = ports.len();
    let r = rom.k_r.len();
    if r < n_ports {
        return Err(Error::invalid(format!("CircuitRom dimension r = {r} is smaller than port count J = {n_ports}")));
    }

    let port_trunc = params.get("port_trunc").and_then(|v| v.as_u64()).unwrap_or(8) as usize;
    let n_flux_basis = params.get("n_flux_basis").and_then(|v| v.as_u64()).unwrap_or(40) as usize;
    let n_charge_cut = params.get("n_charge_cut").and_then(|v| v.as_u64()).unwrap_or(15) as isize;
    if port_trunc < 3 {
        return Err(Error::invalid("port_trunc must be >= 3"));
    }

    // 1. Construct C_tot,r (F) and L_lin,r^-1 (1/H)
    let mut c_tot_r = vec![vec![0.0; r]; r];
    let mut l_inv_r = vec![vec![0.0; r]; r];
    for i in 0..r {
        for j in 0..r {
            c_tot_r[i][j] = EPS0 * rom.m_r[i][j];
            l_inv_r[i][j] = rom.k_r[i][j] / MU0;
        }
    }
    let scale_flux = jj_flux_scale();
    for (p, port) in ports.iter().enumerate() {
        let cj = port.c_j_f;
        let inv_lp = port.e_l_hz.unwrap_or(0.0).max(0.0) / scale_flux;
        for i in 0..r {
            let bi = rom.b_r[i][p];
            for j in 0..r {
                let bj = rom.b_r[j][p];
                c_tot_r[i][j] += cj * bi * bj;
                l_inv_r[i][j] += inv_lp * bi * bj;
            }
        }
    }

    // 2. QR / Gram-Schmidt decomposition of B_r (r x J) into Q_p (r x J) and nullspace N_c (r x (r - J))
    let mut q_p: Vec<Vec<f64>> = Vec::with_capacity(n_ports);
    for p in 0..n_ports {
        let mut v: Vec<f64> = (0..r).map(|i| rom.b_r[i][p]).collect();
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

    // Invert R_p = Q_p^T B_r (J x J) so T_{p,0} = Q_p R_p^{-T} satisfies B_r^T T_{p,0} = I_J
    let mut r_p = Dense::<C64>::zeros(n_ports, n_ports);
    for i in 0..n_ports {
        for j in 0..n_ports {
            let b_col: Vec<f64> = (0..r).map(|k| rom.b_r[k][j]).collect();
            r_p[(j, i)] = C64::new(dot_vec(&q_p[i], &b_col), 0.0);
        }
    }
    let mut t_p: Vec<Vec<f64>> = Vec::with_capacity(n_ports);
    for p in 0..n_ports {
        let mut rhs = vec![C64::ZERO; n_ports];
        rhs[p] = C64::ONE;
        let alpha = femcore::linalg::dense::solve_small(&r_p, &rhs)
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

    // Complete Q_p to an orthonormal basis with N_c (r - J columns in ker(B_r^T))
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

    if n_int > 0 {
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

        // Relax T_p against positive-stiffness interior modes when r_static > 0 (3D FEM static-inductive lift relaxation)
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

        // Select physical short-circuit cavity modes (above 1 MHz)
        let pos_modes: Vec<(f64, f64, Vec<f64>)> = all_tc
            .into_iter()
            .filter_map(|(lam, tc)| {
                if lam > 0.0 {
                    let omega = lam.sqrt();
                    let freq = omega / (2.0 * PI);
                    (freq > 1.0e6).then_some((freq, omega, tc))
                } else {
                    None
                }
            })
            .collect();

        let max_cav = params.get("fock_trunc").and_then(|v| v.as_array()).map(|a| a.len()).unwrap_or(usize::MAX);

        if !rom.modes.is_empty() && pos_modes.len() > 1 {
            // Identify cavity resonator modes in `rom.modes` (modes with total port inductive participation < 0.5)
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
            picked.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
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

    // Reconcile 3D H^1 `c_port_f` when provided
    let has_exact_c_port = rom.c_port_f.len() == n_ports
        && (0..n_ports).all(|j| rom.c_port_f[j].len() == n_ports && rom.c_port_f[j][j] + ports[j].c_j_f > 0.0);
    if has_exact_c_port {
        let mut c_port_tot = rom.c_port_f.clone();
        for j in 0..n_ports {
            c_port_tot[j][j] += ports[j].c_j_f;
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

    // 6. Prepare 1D nonlinear port models (caching `expm` once per Fluxonium port)
    let mut prepared_ports = Vec::with_capacity(n_ports);
    let mut all_fluxonium = true;
    for (j, port) in ports.iter().enumerate() {
        let (_, ej_hz) = resolve_jj_l_and_ej(port.l_j_h.unwrap_or(0.0), port.e_j_hz.unwrap_or(0.0));
        let ec_hz = ((E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H)) * c_inv[j][j];
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
            prepared_ports.push(PreparedPort1D::Fluxonium { n_flux, e_j_hz: ej_hz, h0_diag, cos_phi, sin_phi, y_op });
        } else {
            all_fluxonium = false;
            let n_states = (2 * n_charge_cut + 1) as usize;
            prepared_ports.push(PreparedPort1D::ChargePeriodic {
                n_states,
                n_cut: n_charge_cut,
                e_c_hz: ec_hz,
                e_j_hz: ej_hz,
                n_g: port.n_g,
            });
        }
    }

    // 7. Cavity parameters and bilinear coupling rates
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

    let hbar = PLANCK_H / (2.0 * PI);
    let mut f_cav_hz = Vec::with_capacity(m_cav);
    let mut q_zpf_cav = Vec::with_capacity(m_cav);
    for m in 0..m_cav {
        let l_m_inv = omega_c_keep[m] * omega_c_keep[m] * c_scale;
        let c_m_inv = c_inv[n_ports + m][n_ports + m];
        let omega_m_tilde = (l_m_inv * c_m_inv).sqrt();
        let z_m_tilde = omega_m_tilde / l_m_inv;
        f_cav_hz.push(omega_m_tilde / (2.0 * PI));
        q_zpf_cav.push((hbar / (2.0 * z_m_tilde)).sqrt());
    }

    // Bilinear port-cavity coupling rates g_{jm} (Hz)
    let mut g_jm = vec![vec![0.0; m_cav]; n_ports];
    for j in 0..n_ports {
        for m in 0..m_cav {
            g_jm[j][m] = ((2.0 * E_CHARGE) / PLANCK_H) * c_inv[j][n_ports + m] * q_zpf_cav[m];
        }
    }

    // Bilinear port-port capacitive coupling coefficients (4 e^2 / h) (C^-1)_{jk} (Hz)
    let mut g_jk = vec![vec![0.0; n_ports]; n_ports];
    for j in 0..n_ports {
        for k in (j + 1)..n_ports {
            let v = ((4.0 * E_CHARGE * E_CHARGE) / PLANCK_H) * c_inv[j][k];
            g_jk[j][k] = v;
            g_jk[k][j] = v;
        }
    }

    // Precompute cavity real skew-symmetric operators `(a_m^dagger - a_m)` (`D_m x D_m`)
    let mut cav_adag_minus_a: Vec<Vec<f64>> = Vec::with_capacity(m_cav);
    for &dm in cav_dims.iter().take(m_cav) {
        let mut mat = vec![0.0; dm * dm];
        for n in 0..(dm - 1) {
            let s = ((n + 1) as f64).sqrt();
            // (a^dagger - a)_{n+1, n} = +sqrt(n+1), (a^dagger - a)_{n, n+1} = -sqrt(n+1)
            mat[(n + 1) * dm + n] = s;
            mat[n * dm + (n + 1)] = -s;
        }
        cav_adag_minus_a.push(mat);
    }

    let solve_cb_at_biases = |phi_exts: &[f64]| -> Result<LabeledSpectrum, Error> {
        let mut solved_ports = Vec::with_capacity(n_ports);
        for j in 0..n_ports {
            solved_ports.push(prepared_ports[j].solve(phi_exts[j], port_trunc)?);
        }

        if all_fluxonium {
            // Fast 100% real-symmetric path:
            // n_j (x) i (a_m^dagger - a_m) = (-i y_j) (x) i (a_m^dagger - a_m) = + y_j (x) (a_m^dagger - a_m)
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

            // Add port-cavity coupling: g_{jm} y_j (x) (a_m^dagger - a_m)
            for j in 0..n_ports {
                let y_j = solved_ports[j].y_real.as_ref().unwrap();
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
                        // (a_m^dagger - a_m) is only non-zero for mm = nm +- 1
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

            // Add port-port coupling if J > 1: g_{jk} n_j (x) n_k = - g_{jk} y_j (x) y_k
            for j in 0..n_ports {
                for k in (j + 1)..n_ports {
                    let g = g_jk[j][k];
                    if g.abs() <= 1e-12 {
                        continue;
                    }
                    let y_j = solved_ports[j].y_real.as_ref().unwrap();
                    let y_k = solved_ports[k].y_real.as_ref().unwrap();
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

            let (evals, evecs) = fast_sym_eigh(&h_real, d_tot)?;
            Ok(label_spectrum_bipartite(&dims, &evals, |b_idx, k| {
                let v = evecs[k][b_idx];
                v * v
            }))
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
                    if let Some(ref nc) = sp.n_cplx {
                        nc.clone()
                    } else {
                        let y = sp.y_real.as_ref().unwrap();
                        let d = dims[j];
                        Dense::from_fn(d, d, |r_i, c_i| C64::new(0.0, -y[r_i * d + c_i]))
                    }
                })
                .collect();

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
            Ok(label_spectrum_bipartite(&dims, &evals, |b_idx, k| evecs[(b_idx, k)].norm_sqr()))
        }
    };

    let base_phi_exts: Vec<f64> = ports.iter().map(|p| p.phi_ext_rad).collect();
    let spec = solve_cb_at_biases(&base_phi_exts)?;
    let elapsed_s = t0.elapsed().as_secs_f64();

    let (sweep, sweep_elapsed_s) = if let Some(sweep_pts) = parse_flux_sweep(params) {
        let t_sw = Instant::now();
        let n_pts = sweep_pts.len();
        let n_threads = femcore::par::threads_up_to(8).min(n_pts.max(1));
        let chunk_size = n_pts.div_ceil(n_threads.max(1));
        let chunks: Vec<Result<Vec<FluxSweepPoint>, Error>> =
            femcore::par::map_chunks(&sweep_pts, chunk_size, |slice| {
                let mut local = Vec::with_capacity(slice.len());
                let mut phis = base_phi_exts.clone();
                for &phi in slice {
                    phis[0] = phi;
                    let sp = solve_cb_at_biases(&phis)?;
                    local.push(FluxSweepPoint {
                        phi_ext_rad: phi,
                        flux_phi0: phi / (2.0 * PI),
                        f01_hz: sp.f01_hz,
                        f12_hz: sp.f12_hz,
                        f02_hz: sp.f02_hz,
                        alpha_hz: sp.alpha_hz,
                        chi_matrix_hz: sp.chi_matrix_hz,
                        eigenvalues_hz: sp.eigenvalues_hz,
                        min_overlap: sp.min_overlap,
                    });
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

    let mut subsystems: Vec<String> = ports.into_iter().map(|p| p.name).collect();
    for m in 0..m_cav {
        subsystems.push(format!("cavity_{}", m + 1));
    }

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
        elapsed_s,
        f_cav_hz,
        g_matrix_hz: g_jm,
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
        // Exact 2-body LC resonator (f_r = 7.0 GHz) + Fluxonium (E_C/h = 1.0 GHz, E_L/h = 1.0 GHz, E_J/h = 4.0 GHz)
        // at Phi_ext = 0.5 Phi_0 (phi_ext = pi)
        let ec_hz = 1.0e9_f64;
        let el_hz = 1.0e9_f64;
        let ej_hz = 4.0e9_f64;
        let f_r_bare = 7.0e9_f64;
        let omega_r = 2.0 * PI * f_r_bare;

        let c_sigma = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * ec_hz);
        let c_r = 50.0e-15_f64;
        let c_g = 0.50e-15_f64;
        // Choose c_q so that Schur complement (C^-1)_{11} = 1 / c_sigma => exact E_C / h = 1.0 GHz
        let c_11 = c_sigma + (c_g * c_g) / (c_r + c_g);
        let l_r_inv = omega_r * omega_r * (c_r + c_g);

        // Build synthetic 2x2 CircuitRom with exact 2-body capacitance and reluctance matrices:
        // Coordinate 0: Fluxonium port flux Phi_1; Coordinate 1: Resonator flux Phi_2
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
        // Reference bare Fluxonium (E_C=1.0, E_L=1.0, E_J=4.0, flux=0.5, cutoff=40) has f_01 = 581.849 MHz
        assert!((rep.f01_hz[0] - 581.849e6).abs() / 581.849e6 < 5e-3, "Fluxonium f01 = {:.4} MHz", rep.f01_hz[0] / 1e6);
        assert!((rep.f01_hz[1] - 7.0e9).abs() / 7.0e9 < 5e-3, "Resonator f01 = {:.6} GHz", rep.f01_hz[1] / 1e9);
        assert!(rep.min_overlap > 0.98, "min_overlap = {}", rep.min_overlap);
    }

    #[test]
    fn normal_fock_and_charge_periodic_craig_bampton_agree_on_transmon() {
        let scale = jj_flux_scale();
        let e_j_hz = 16.0e9;
        let l_j_h = scale / e_j_hz;
        let c_port = 75.0e-15;
        let b0 = (EPS0 / c_port).sqrt();

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
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 1,
            r_dyn: 0,
            k_r: vec![vec![0.0]],
            c_r: vec![],
            m_r: vec![vec![1.0]],
            b_r: vec![vec![b0]],
            c_port_f: vec![vec![c_port]],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

        let rep_fock = hamiltonian_with_params(&res, &json!({ "basis": "normal_fock", "fock_trunc": [10] })).unwrap();
        let rep_cb = hamiltonian_with_params(&res, &json!({ "basis": "craig_bampton", "port_trunc": 8 })).unwrap();

        let rel_f01 = (rep_fock.f01_hz[0] - rep_cb.f01_hz[0]).abs() / rep_cb.f01_hz[0];
        let rel_alpha = (rep_fock.alpha_hz[0] - rep_cb.alpha_hz[0]).abs() / rep_cb.alpha_hz[0].abs();
        assert!(
            rel_f01 < 3e-3,
            "normal_fock f01 ({:.4} GHz) vs charge-basis CB ({:.4} GHz): rel diff {rel_f01:e}",
            rep_fock.f01_hz[0] / 1e9,
            rep_cb.f01_hz[0] / 1e9
        );
        assert!(
            rel_alpha < 5e-2,
            "normal_fock alpha ({:.2} MHz) vs charge-basis CB ({:.2} MHz): rel diff {rel_alpha:e}",
            rep_fock.alpha_hz[0] / 1e6,
            rep_cb.alpha_hz[0] / 1e6
        );
    }
}
