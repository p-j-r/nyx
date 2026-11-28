// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Cross-sectional FEM eigenmode bridge (`nyx-ems::modes` & `nyx-ems::modes_complex`) for FDTD `ModeSource`
//! unidirectional injection and `ModeMonitor` broadband modal S-parameter decomposition (ADR 145).

use super::ade_pml::{FdtdBcFace, FdtdBoundaries};
use super::grid::{FdtdDim, Pol2d, YeeGrid, C0, MU0};
use super::subpixel::{DomainMedium, FdtdGeometryOracle, PointClass};
use crate::modes::{assemble_modes_with_bcs_and_symmetry, solve_mode_pairs, EpsCell};
use crate::modes_complex::{assemble_modes_complex_with_bcs, solve_mode_pairs_complex};
use crate::problem::{Problem, Symmetry};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::fem::{geometry, shape};
use femcore::linalg::csr::{c, C};
use femcore::mesh::SimplexMesh;
use femcore::optical::{EpsC, MuC};
use std::collections::BTreeMap;
use std::f64::consts::PI;

/// Solved and Yee-interpolated cross-sectional waveguide mode at one reference frequency `freq_hz`.
#[derive(Clone, Debug)]
pub struct YeeModeProfile {
    pub freq_hz: f64,
    /// Propagation axis (`0 = x`, `1 = y`, `2 = z`).
    pub axis: usize,
    /// Nominal propagation sign (`+1.0` along `+axis`, `-1.0` along `-axis`).
    pub sign: f64,
    /// Primal plane index `i_plane` on `axis` where tangential `E` lives (`v_axis[i_plane]`).
    pub plane_idx: usize,
    /// Exact plane coordinate `v_axis[plane_idx]` (m).
    pub plane_pos_m: f64,
    /// Dual cell index `i_h` adjacent to `plane_idx` for TF/SF magnetic injection (`plane_idx - 1` for `+1`, `plane_idx` for `-1`).
    pub h_cell_idx: usize,
    /// Complex effective index `n_eff` (`Re(n_eff) > 0`, `Im(n_eff) <= 0` passive) from `nyx-ems::modes` / `modes_complex`.
    pub neff: C,
    /// Complex propagation constant `beta = k0 * neff` (rad/m).
    pub beta: C,
    /// Modal wave impedance `Z_mode = omega * mu_0 / Re(beta)` (Ohm).
    pub z_mode_ohm: f64,
    /// Tangential axis indices `(u_axis, v_axis)` chosen so `(e_u, e_v, e_axis)` is right-handed (`e_u x e_v = +e_axis`)!
    pub u_axis: usize,
    pub v_axis: usize,
    /// Sampled tangential `E_u`, `E_v` (V/m, complex `[re, im]`) at every staggered Yee tangential E-edge on the port window,
    /// normalized so forward discrete cross-power `0.5 * Re sum (E x H^*) . (+sign * e_axis) dA == 1 W`.
    pub e_u_samples: Vec<(usize, C)>,
    pub e_v_samples: Vec<(usize, C)>,
    /// Sampled tangential `H_u`, `H_v` (A/m, complex `[re, im]`) at adjacent staggered Yee H-faces (`h_cell_idx`),
    /// corresponding to a `+sign` forward-propagating wave of unit power (`1 W`) at `plane_pos_m`.
    pub h_u_samples: Vec<(usize, C)>,
    pub h_v_samples: Vec<(usize, C)>,
    /// Co-located modal `(e_u, e_v, h_u, h_v, area_m2)` on the dual face centers `(u_c, v_c)` of the port window
    /// for accurate `ModeMonitor` overlap integrals `(a_fwd, b_bwd)`.
    pub monitor_patches: Vec<ModePatchSample>,
}

/// One quadrature cell patch `(iu, iv)` on the port plane `(u_axis, v_axis)` carrying co-located modal `(E_t, H_t)` for `ModeMonitor`.
#[derive(Clone, Copy, Debug)]
pub struct ModePatchSample {
    pub iu: usize,
    pub iv: usize,
    pub pos_3d: [f64; 3],
    pub area_m2: f64,
    /// Forward (`+sign * e_axis`) modal tangential electric and magnetic fields (`e_u`, `e_v`, `h_u`, `h_v`).
    pub e_u: C,
    pub e_v: C,
    pub h_u: C,
    pub h_v: C,
}

/// Right-handed transverse axes `(u, v)` for propagation axis `axis` (`e_u x e_v == +e_axis`).
#[inline]
pub fn rh_uv_axes(axis: usize) -> (usize, usize) {
    match axis {
        0 => (1, 2), // y x z = +x
        1 => (2, 0), // z x x = +y
        _ => (0, 1), // x x y = +z
    }
}

/// Solves a 2D cross-sectional eigenproblem via `nyx-ems::modes` / `modes_complex` on the specified Yee port plane
/// and interpolates the `mode_index`-th (`1`-based) guided mode onto the staggered Yee grid!
#[allow(clippy::too_many_arguments)]
pub fn solve_yee_plane_mode(
    grid: &YeeGrid,
    oracle: &FdtdGeometryOracle,
    problem: &Problem,
    boundaries: &FdtdBoundaries,
    media: &BTreeMap<i32, DomainMedium>,
    default_medium: DomainMedium,
    axis: usize,
    position_m: f64,
    sign: f64,
    window_lo: Option<[f64; 2]>,
    window_hi: Option<[f64; 2]>,
    freq_hz: f64,
    mode_index: usize,
    n_eff_guess: Option<f64>,
) -> Result<YeeModeProfile, String> {
    if !(freq_hz.is_finite() && freq_hz > 0.0) {
        return Err(format!("FDTD mode solver requires positive finite freq_hz (got {freq_hz})"));
    }
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    if is_2d && axis == 1 {
        return Err(
            "2D FDTD in wafer (x, z) plane cannot place a ModeSource/ModeMonitor normal to out-of-plane axis y".into(),
        );
    }
    let n_ax = grid.axis_centers(axis).len();
    if n_ax < 2 {
        return Err(format!("FDTD grid along axis {axis} needs at least 2 cells for ModeSource/ModeMonitor"));
    }
    // Choose interior primal plane `1 ..= n_ax - 1` nearest `position_m` so both `plane_idx - 1` and `plane_idx` exist!
    let raw_p = grid.nearest_vertex(axis, position_m);
    let plane_idx = raw_p.clamp(1, n_ax - 1);
    let plane_pos_m = grid.axis_vertices(axis)[plane_idx];
    let dir_sign = if sign >= 0.0 { 1.0 } else { -1.0 };
    let h_cell_idx = if dir_sign > 0.0 { plane_idx - 1 } else { plane_idx };

    let (u_ax, v_ax) = rh_uv_axes(axis);
    let vu_all = grid.axis_vertices(u_ax);
    let vv_all = grid.axis_vertices(v_ax);
    let cu_all = grid.axis_centers(u_ax);
    let cv_all = grid.axis_centers(v_ax);
    let dpu_all = grid.axis_dp(u_ax);
    let dpv_all = grid.axis_dp(v_ax);

    let (iu0, iu1) = window_cell_range(vu_all, window_lo.map(|w| w[0]), window_hi.map(|w| w[0]));
    let (iv0, iv1) = if is_2d {
        // In 2D `(x, z)`, one of `(u_ax, v_ax)` is `1` (`y`, where `ny == 1`, `vy = [-0.5, 0.5]`).
        (0, 1)
    } else {
        window_cell_range(vv_all, window_lo.map(|w| w[1]), window_hi.map(|w| w[1]))
    };
    let (iu0, iu1, iv0, iv1) = if is_2d && u_ax == 1 {
        let (jv0, jv1) = window_cell_range(
            vv_all,
            window_lo.map(|w| w[0]).or(window_lo.map(|w| w[1])),
            window_hi.map(|w| w[0]).or(window_hi.map(|w| w[1])),
        );
        (0, 1, jv0, jv1)
    } else {
        (iu0, iu1, iv0, iv1)
    };

    // Build 2D cross-section `SimplexMesh` in `(s0, s1)` where solver local axes `(s0, s1, s2)` satisfy:
    // Recall `nyx-ems::modes` uses LEFT-HANDED `(s0, s1, s2)`!
    // Wait: if we map `(s0, s1) = (v_ax, u_ax)` and `s2 = +axis`, let's check handedness!
    // Since `(e_u, e_v, e_axis)` is RIGHT-HANDED (`e_u x e_v = +e_axis`), swapping the two transverse axes
    // `(e_v, e_u, e_axis)` makes `e_v x e_u = -e_axis`, which is LEFT-HANDED — IDENTICAL to `(x, z, +y)` in `modes.rs`!
    // Let's double check! For `axis == 1` (`+y`), `rh_uv_axes(1) = (2, 0) = (z, x)`, so `(s0, s1) = (v_ax, u_ax) = (0, 2) = (x, z)`!!
    // BINGO!! Setting `(s0_ax, s1_ax) = (v_ax, u_ax)` reproduces `(x, z)` identically when `axis == 1` (`+y`),
    // and generalizes `modes.rs`'s exact curl formulas to ANY axis `0, 1, 2` with ZERO sign changes!
    let (s0_ax, s1_ax) = (v_ax, u_ax);
    let (is0_lo, is0_hi) = (iv0, iv1);
    let (is1_lo, is1_hi) = (iu0, iu1);
    let vs0_all = grid.axis_vertices(s0_ax);
    let vs1_all = grid.axis_vertices(s1_ax);
    let vs0 = &vs0_all[is0_lo..=is0_hi];
    let vs1 = &vs1_all[is1_lo..=is1_hi];
    let outer_flags = [is0_lo == 0, is0_hi + 1 == vs0_all.len(), is1_lo == 0, is1_hi + 1 == vs1_all.len()];

    let (mesh_2d, cell_eps_3x3, cell_mu_r, max_n_slice, pmc_labels_2d, pec_labels_2d) = build_cross_section_mesh_2d(
        grid,
        oracle,
        boundaries,
        media,
        default_medium,
        axis,
        plane_pos_m,
        dir_sign,
        s0_ax,
        s1_ax,
        vs0,
        vs1,
        outer_flags,
    )?;

    let n_target = n_eff_guess.unwrap_or(max_n_slice * 0.98).max(1.0);
    let order = problem.order.clamp(1, 2);
    let k_req = (mode_index + 3).max(4);
    let k0 = 2.0 * PI * freq_hz / C0;
    let omega = 2.0 * PI * freq_hz;

    // Solve cross-section eigenproblem via `nyx-ems::modes` (or `modes_complex` if anisotropic/lossy)
    let mut use_complex = false;
    let mut eps_real_cells = Vec::with_capacity(mesh_2d.n_cells());
    for ep in &cell_eps_3x3 {
        if ep[0][2].abs() > 1e-12 || ep[1][2].abs() > 1e-12 || (ep[0][1] - ep[1][0]).abs() > 1e-12 {
            use_complex = true;
        }
        eps_real_cells.push(EpsCell { xx: ep[0][0], xy: ep[0][1], yy: ep[1][1], zz: ep[2][2] });
    }

    let (neff, et_vec, phi_vec) = if !use_complex {
        let asm = assemble_modes_with_bcs_and_symmetry(
            &mesh_2d,
            freq_hz,
            &eps_real_cells,
            &cell_mu_r,
            order,
            &pec_labels_2d,
            &[],
            &pmc_labels_2d,
            Symmetry::None,
        )?;
        let k_solve = k_req.min(asm.a.nrows().saturating_sub(1)).max(1);
        let modes = solve_mode_pairs(&asm, n_target, false, k_solve)?;
        if modes.neff.is_empty() {
            return Err("FDTD ModeSource/ModeMonitor: 2D FEM eigensolver found zero modes on port cross-section".into());
        }
        // Sort converged propagating modes by descending Re(neff) (matching `WavePortBc` in `nyx-ems::driven`)
        let mut idxs: Vec<usize> = (0..modes.neff.len()).collect();
        idxs.sort_by(|&a, &b| modes.neff[b].0.abs().total_cmp(&modes.neff[a].0.abs()));
        let pick = idxs.get(mode_index.saturating_sub(1)).copied().unwrap_or(idxs[0]);
        let nf = c(modes.neff[pick].0.abs(), -modes.neff[pick].1.abs());
        let et_c: Vec<C> = modes.et[pick].iter().map(|&v| c(v, 0.0)).collect();
        let phi_c: Vec<C> = modes.phi[pick].iter().map(|&v| c(v, 0.0)).collect();
        (nf, et_c, phi_c)
    } else {
        let eps_c: Vec<EpsC> = cell_eps_3x3
            .iter()
            .map(|m| {
                EpsC::symmetric(
                    c(m[0][0], 0.0),
                    c(m[0][1], 0.0),
                    c(m[1][1], 0.0),
                    c(m[2][2], 0.0),
                    c(m[0][2], 0.0),
                    c(m[1][2], 0.0),
                )
            })
            .collect();
        let mu_c: Vec<MuC> = cell_mu_r.iter().map(|&mr| MuC::isotropic(c(mr, 0.0))).collect();
        let asm = assemble_modes_complex_with_bcs(
            &mesh_2d,
            freq_hz,
            &eps_c,
            &mu_c,
            order,
            &pec_labels_2d,
            &[],
            &pmc_labels_2d,
        )?;
        let k_solve = k_req.min(asm.a.nrows().saturating_sub(1)).max(1);
        let mc = solve_mode_pairs_complex(&asm, c(n_target, 0.0), k_solve)?;
        if mc.neff.is_empty() {
            return Err("FDTD ModeSource/ModeMonitor: complex 2D FEM eigensolver found zero modes".into());
        }
        let mut idxs: Vec<usize> = (0..mc.neff.len()).collect();
        idxs.sort_by(|&a, &b| mc.neff[b].re.abs().total_cmp(&mc.neff[a].re.abs()));
        let pick = idxs.get(mode_index.saturating_sub(1)).copied().unwrap_or(idxs[0]);
        let mut nf = mc.neff[pick];
        if nf.re < 0.0 {
            nf = -nf;
        }
        if nf.im > 0.0 {
            nf.im = -nf.im;
        }
        (nf, mc.et[pick].clone(), mc.phi[pick].clone())
    };

    let beta_plus = neff * k0;
    // Build fast evaluator over `mesh_2d` returning `(E_{s0}, E_{s1}, H_{s0}, H_{s1})` for `+axis` wave!
    // Remember `(s0, s1, +axis)` is identical to `(x, z, +y)` in `eme.rs`:
    //   `e_s0 = er0`, `e_s1 = er1`,
    //   `h_s0 = (-j / (omega * mu)) * (+j * beta * (dp1 + er1))`
    //   `h_s1 = (-j / (omega * mu)) * (-j * beta * (dp0 + er0))`
    // Notice for a `-axis` forward wave (`dir_sign == -1`), reversing propagation direction (`beta -> -beta`, `phi -> -phi`)
    // keeps `E_t = (E_u, E_v)` unchanged and flips `H_t -> -H_t` so `(E_t x H_t^*) . (-e_axis) > 0`!
    let sampler = SliceFemSampler::new(&mesh_2d, order, omega, beta_plus, &cell_mu_r, &et_vec, &phi_vec);

    // Fix deterministic canonical phase (largest transverse electric component real and positive)
    let mut best_amp = c(0.0, 0.0);
    let mut patches_raw = Vec::with_capacity((iu1 - iu0) * (iv1 - iv0));
    for iu in iu0..iu1 {
        for iv in iv0..iv1 {
            let u_c = cu_all[iu];
            let v_c = cv_all[iv];
            let (s0, s1) = (v_c, u_c);
            let (es0, es1, hs0, hs1) = sampler.eval_eh([s0, s1]);
            // Map `(s0, s1) = (v_ax, u_ax)` -> `(u, v)`: `e_u = es1`, `e_v = es0`, and for `dir_sign` wave: `h_u = dir_sign * hs1`, `h_v = dir_sign * hs0`!
            // Check right-handed Poynting along `dir_sign * e_axis`:
            // `(E x H^*) . (dir_sign * e_axis) = dir_sign * (e_u * h_v^* - e_v * h_u^*)`
            // `= dir_sign * (es1 * dir_sign * hs0^* - es0 * dir_sign * hs1^*) = es1 * hs0^* - es0 * hs1^*`!
            // Compare with `eme.rs` line 3126: `el1.conj() * h0 - el0.conj() * h1` (`Re` is identical!).
            let eu = es1;
            let ev = es0;
            let hu = hs1 * dir_sign;
            let hv = hs0 * dir_sign;
            if eu.norm_sqr() > best_amp.norm_sqr() {
                best_amp = eu;
            }
            if ev.norm_sqr() > best_amp.norm_sqr() {
                best_amp = ev;
            }
            let area = dpu_all[iu] * dpv_all[iv];
            let mut p3 = [0.0; 3];
            p3[axis] = plane_pos_m;
            p3[u_ax] = u_c;
            p3[v_ax] = v_c;
            patches_raw.push(ModePatchSample { iu, iv, pos_3d: p3, area_m2: area, e_u: eu, e_v: ev, h_u: hu, h_v: hv });
        }
    }

    let phase_rot = if best_amp.norm() > 1e-30 { best_amp.conj() / best_amp.norm() } else { c(1.0, 0.0) };
    let mut p_sum = 0.0;
    for mp in &patches_raw {
        // P = 0.5 * Re int (E x H^*) . (dir_sign * e_axis) dA = 0.5 * dir_sign * Re(e_u * h_v^* - e_v * h_u^*) dA
        let s_n = dir_sign * (mp.e_u * mp.h_v.conj() - mp.e_v * mp.h_u.conj()).re;
        p_sum += 0.5 * s_n * mp.area_m2;
    }
    let flip_h = if p_sum < 0.0 { -1.0 } else { 1.0 };
    let p_abs = p_sum.abs().max(1e-30);
    let scale = phase_rot * (1.0 / p_abs.sqrt());

    for mp in &mut patches_raw {
        mp.e_u *= scale;
        mp.e_v *= scale;
        mp.h_u *= scale * flip_h;
        mp.h_v *= scale * flip_h;
    }

    // Sample staggered Yee `E_u`, `E_v` on `plane_idx` and `H_u`, `H_v` on `h_cell_idx`!
    let mut e_u_samples = Vec::new();
    let mut e_v_samples = Vec::new();
    let mut h_u_samples = Vec::new();
    let mut h_v_samples = Vec::new();

    // Helper mapping `(i_axis, i_u, i_v)` -> `(ix, iy, iz)` -> `grid.idx`
    let to_grid_idx = |ia: usize, iu: usize, iv: usize| -> usize {
        let mut ijk = [0usize; 3];
        ijk[axis] = ia;
        ijk[u_ax] = iu;
        ijk[v_ax] = iv;
        grid.idx(ijk[0], ijk[1], ijk[2])
    };

    // E_u lives at `(v_axis[plane_idx], c_u[iu], v_v[iv])` for `iu in iu0..iu1`, `iv in iv0..=iv1`
    let iv_e_hi = if is_2d && v_ax == 1 { 0 } else { iv1 };
    for iu in iu0..iu1 {
        let u = cu_all[iu];
        for iv in iv0..=iv_e_hi {
            let v = if is_2d && v_ax == 1 { 0.0 } else { vv_all[iv] };
            let (_, es1, _, _) = sampler.eval_eh([v, u]);
            let val = es1 * scale;
            if val.norm() > 1e-18 {
                e_u_samples.push((to_grid_idx(plane_idx, iu, iv), val));
            }
        }
    }

    // E_v lives at `(v_axis[plane_idx], v_u[iu], c_v[iv])` for `iu in iu0..=iu1`, `iv in iv0..iv1`
    let iu_e_hi = if is_2d && u_ax == 1 { 0 } else { iu1 };
    for iu in iu0..=iu_e_hi {
        let u = if is_2d && u_ax == 1 { 0.0 } else { vu_all[iu] };
        for iv in iv0..iv1 {
            let v = cv_all[iv];
            let (es0, _, _, _) = sampler.eval_eh([v, u]);
            let val = es0 * scale;
            if val.norm() > 1e-18 {
                e_v_samples.push((to_grid_idx(plane_idx, iu, iv), val));
            }
        }
    }

    // H_u lives at `(c_axis[h_cell_idx], v_u[iu], c_v[iv])` for `iu in iu0..=iu1`, `iv in iv0..iv1`
    for iu in iu0..=iu_e_hi {
        let u = if is_2d && u_ax == 1 { 0.0 } else { vu_all[iu] };
        for iv in iv0..iv1 {
            let v = cv_all[iv];
            let (_, _, _, hs1) = sampler.eval_eh([v, u]);
            let val = hs1 * (dir_sign * flip_h) * scale;
            if val.norm() > 1e-18 {
                h_u_samples.push((to_grid_idx(h_cell_idx, iu, iv), val));
            }
        }
    }

    // H_v lives at `(c_axis[h_cell_idx], c_u[iu], v_v[iv])` for `iu in iu0..iu1`, `iv in iv0..=iv1`
    for iu in iu0..iu1 {
        let u = cu_all[iu];
        for iv in iv0..=iv_e_hi {
            let v = if is_2d && v_ax == 1 { 0.0 } else { vv_all[iv] };
            let (_, _, hs0, _) = sampler.eval_eh([v, u]);
            let val = hs0 * (dir_sign * flip_h) * scale;
            if val.norm() > 1e-18 {
                h_v_samples.push((to_grid_idx(h_cell_idx, iu, iv), val));
            }
        }
    }

    let z_mode_ohm = (omega * MU0 / beta_plus.re.abs().max(1e-12)).max(1e-6);

    Ok(YeeModeProfile {
        freq_hz,
        axis,
        sign: dir_sign,
        plane_idx,
        plane_pos_m,
        h_cell_idx,
        neff,
        beta: beta_plus,
        z_mode_ohm,
        u_axis: u_ax,
        v_axis: v_ax,
        e_u_samples,
        e_v_samples,
        h_u_samples,
        h_v_samples,
        monitor_patches: patches_raw,
    })
}

fn window_cell_range(v: &[f64], lo_opt: Option<f64>, hi_opt: Option<f64>) -> (usize, usize) {
    let n = v.len() - 1;
    let lo = lo_opt.unwrap_or(v[0]);
    let hi = hi_opt.unwrap_or(v[n]);
    let mut i0 = 0usize;
    while i0 + 1 < n && v[i0 + 1] <= lo + 1e-15 {
        i0 += 1;
    }
    let mut i1 = n;
    while i1 > i0 + 1 && v[i1 - 1] >= hi - 1e-15 {
        i1 -= 1;
    }
    (i0, i1.max(i0 + 1))
}

/// Builds a 2D `SimplexMesh` on the cross-section `(vs0, vs1)` at `axis = plane_pos_m`, sampling the `FdtdGeometryOracle`
/// (with sub-grid refinement in 2D so curved/thin waveguide cores are captured accurately by `nyx-ems::modes`).
#[allow(clippy::too_many_arguments)]
fn build_cross_section_mesh_2d(
    grid: &YeeGrid,
    oracle: &FdtdGeometryOracle,
    boundaries: &FdtdBoundaries,
    media: &BTreeMap<i32, DomainMedium>,
    default_medium: DomainMedium,
    axis: usize,
    plane_pos_m: f64,
    dir_sign: f64,
    s0_ax: usize,
    s1_ax: usize,
    vs0_raw: &[f64],
    vs1_raw: &[f64],
    outer_flags: [bool; 4],
) -> Result<(SimplexMesh, Vec<[[f64; 3]; 3]>, Vec<f64>, f64, Vec<i32>, Vec<i32>), String> {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    // Refine 2x along coarse axes (`<= 24` vertices) so P2 FEM mode shapes are well-resolved without inflating dense grids
    let refine_1d = |v: &[f64], dummy: bool| -> Vec<f64> {
        if dummy || v.len() > 24 {
            return v.to_vec();
        }
        let mut out = Vec::with_capacity(2 * v.len() - 1);
        for w in v.windows(2) {
            out.push(w[0]);
            out.push(0.5 * (w[0] + w[1]));
        }
        if let Some(&last) = v.last() {
            out.push(last);
        }
        out
    };
    let s0_dummy = is_2d && s0_ax == 1;
    let s1_dummy = is_2d && s1_ax == 1;
    let vs0 = refine_1d(vs0_raw, s0_dummy);
    let vs1 = refine_1d(vs1_raw, s1_dummy);

    let v_axis = grid.axis_vertices(axis);
    let dp_nudge = grid.axis_dp(axis).first().copied().unwrap_or(1e-9);
    let sample_axis_pos =
        (plane_pos_m + dir_sign * 1e-4 * dp_nudge).clamp(v_axis[0] + 1e-12, v_axis[v_axis.len() - 1] - 1e-12);

    let n0 = vs0.len() - 1;
    let n1 = vs1.len() - 1;
    let nv0 = vs0.len();
    let mut pts = Vec::with_capacity(2 * nv0 * vs1.len());
    for j in 0..=n1 {
        for i in 0..=n0 {
            pts.push(vs0[i]);
            pts.push(vs1[j]);
        }
    }
    let v_id = |i: usize, j: usize| -> u32 { (j * nv0 + i) as u32 };

    let mut cells = Vec::with_capacity(6 * n0 * n1);
    let mut cell_labels = Vec::with_capacity(2 * n0 * n1);
    let mut cell_eps = Vec::with_capacity(2 * n0 * n1);
    let mut cell_mu = Vec::with_capacity(2 * n0 * n1);
    let mut max_n = 1.0_f64;

    let perm_axes = [s0_ax, s1_ax, axis];
    for j in 0..n1 {
        for i in 0..n0 {
            let (v00, v10, v11, v01) = (v_id(i, j), v_id(i + 1, j), v_id(i + 1, j + 1), v_id(i, j + 1));
            let tris = if 0.5 * (vs0[i] + vs0[i + 1]) < 0.0 {
                [[v10, v11, v01], [v10, v01, v00]]
            } else {
                [[v00, v10, v11], [v00, v11, v01]]
            };
            for tri in tris {
                let c0 = (pts[2 * tri[0] as usize] + pts[2 * tri[1] as usize] + pts[2 * tri[2] as usize]) / 3.0;
                let c1 =
                    (pts[2 * tri[0] as usize + 1] + pts[2 * tri[1] as usize + 1] + pts[2 * tri[2] as usize + 1]) / 3.0;
                let mut p3 = [0.0; 3];
                p3[axis] = sample_axis_pos;
                p3[s0_ax] = if s0_dummy { 0.0 } else { c0 };
                p3[s1_ax] = if s1_dummy { 0.0 } else { c1 };
                let (dom, med) = match oracle.classify(p3) {
                    PointClass::Domain(d) => (d, media.get(&d).copied().unwrap_or(default_medium)),
                    PointClass::Pec => (
                        1,
                        DomainMedium { eps_inf: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], ..default_medium },
                    ),
                    PointClass::Outside => (1, default_medium),
                };
                let mut ep_loc = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        ep_loc[a][b] = med.eps_inf[perm_axes[a]][perm_axes[b]];
                    }
                }
                let mr = med.mu_r[axis];
                let n_here = (ep_loc[0][0].max(ep_loc[1][1]).max(ep_loc[2][2]) * mr).max(1.0).sqrt();
                max_n = max_n.max(n_here);

                cells.extend_from_slice(&tri);
                cell_labels.push(dom.max(1));
                cell_eps.push(ep_loc);
                cell_mu.push(mr);
            }
        }
    }

    // Boundary facets: label 1 = PEC walls, label 2 = PMC walls:
    // - On dummy 2D `y` walls: `TmEy` -> PEC (`1`), `TeHy` -> PMC (`2`).
    // - On outer YeeGrid boundaries (`outer_flags == true`): `FdtdBcFace::Pmc` -> PMC (`2`), otherwise PEC (`1`).
    // - On interior sub-window cuts (`outer_flags == false`, e.g. adjacent multi-port splitter windows): natural PMC (`2`),
    //   matching `nyx-ems::ports::assemble_waveports_3d` on non-PEC interior waveport perimeter edges!
    let dummy_is_pmc = matches!(grid.dim, FdtdDim::TwoD { polarization: Pol2d::TeHy });
    let side_label = |ax: usize, is_dummy: bool, is_outer: bool, is_hi: bool| -> i32 {
        if is_dummy {
            if dummy_is_pmc {
                2
            } else {
                1
            }
        } else if is_outer {
            let face_idx = 2 * ax + usize::from(is_hi);
            if matches!(boundaries.faces[face_idx], FdtdBcFace::Pmc) {
                2
            } else {
                1
            }
        } else {
            2
        }
    };
    let lab_s1_lo = side_label(s1_ax, s1_dummy, outer_flags[2], false);
    let lab_s1_hi = side_label(s1_ax, s1_dummy, outer_flags[3], true);
    let lab_s0_lo = side_label(s0_ax, s0_dummy, outer_flags[0], false);
    let lab_s0_hi = side_label(s0_ax, s0_dummy, outer_flags[1], true);

    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for i in 0..n0 {
        facets.extend_from_slice(&[v_id(i, 0), v_id(i + 1, 0)]);
        facet_labels.push(lab_s1_lo);
        facets.extend_from_slice(&[v_id(i + 1, n1), v_id(i, n1)]);
        facet_labels.push(lab_s1_hi);
    }
    for j in 0..n1 {
        facets.extend_from_slice(&[v_id(0, j + 1), v_id(0, j)]);
        facet_labels.push(lab_s0_lo);
        facets.extend_from_slice(&[v_id(n0, j), v_id(n0, j + 1)]);
        facet_labels.push(lab_s0_hi);
    }

    let any_pmc = [lab_s1_lo, lab_s1_hi, lab_s0_lo, lab_s0_hi].contains(&2);
    let any_pec = [lab_s1_lo, lab_s1_hi, lab_s0_lo, lab_s0_hi].contains(&1);
    let mesh = SimplexMesh::new(2, pts, cells, cell_labels, facets, facet_labels);
    let pmc_labels = if any_pmc { vec![2] } else { Vec::new() };
    let pec_labels = if any_pec { vec![1] } else { Vec::new() };
    Ok((mesh, cell_eps, cell_mu, max_n, pmc_labels, pec_labels))
}

/// Evaluates transverse `(E_{s0}, E_{s1}, H_{s0}, H_{s1})` at arbitrary 2D points `[s0, s1]` on a solved `SimplexMesh`.
struct SliceFemSampler<'a> {
    mesh: &'a SimplexMesh,
    geo: femcore::fem::Geometry,
    cd: femcore::edge::CellData,
    order: usize,
    omega: f64,
    beta: C,
    mu_r: &'a [f64],
    et: &'a [C],
    phi: &'a [C],
    tris: Vec<([f64; 2], [f64; 2], [f64; 2])>,
}

impl<'a> SliceFemSampler<'a> {
    fn new(
        mesh: &'a SimplexMesh,
        order: usize,
        omega: f64,
        beta: C,
        mu_r: &'a [f64],
        et: &'a [C],
        phi: &'a [C],
    ) -> Self {
        let tris = (0..mesh.n_cells())
            .map(|ci| {
                let cv = mesh.cell(ci);
                let p = |k: usize| {
                    let q = mesh.point(cv[k] as usize);
                    [q[0], q[1]]
                };
                (p(0), p(1), p(2))
            })
            .collect();
        Self { mesh, geo: geometry(mesh), cd: cell_data(mesh, order), order, omega, beta, mu_r, et, phi, tris }
    }

    fn eval_eh(&self, q: [f64; 2]) -> (C, C, C, C) {
        let mut found: Option<(usize, [f64; 3])> = None;
        for (ci, &(p0, p1, p2)) in self.tris.iter().enumerate() {
            let det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
            if det.abs() <= 1e-30 {
                continue;
            }
            let l1 = ((q[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (q[1] - p0[1])) / det;
            let l2 = ((p1[0] - p0[0]) * (q[1] - p0[1]) - (q[0] - p0[0]) * (p1[1] - p0[1])) / det;
            let l0 = 1.0 - l1 - l2;
            if l0 >= -1e-7 && l1 >= -1e-7 && l2 >= -1e-7 {
                found = Some((ci, [l0.clamp(0.0, 1.0), l1.clamp(0.0, 1.0), l2.clamp(0.0, 1.0)]));
                break;
            }
        }
        let Some((ci, lam)) = found else {
            return (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        };

        let nloc_e = if self.order == 1 { 3 } else { 8 };
        let nl = self.mesh.nloc(self.order);
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        let g = tri_g(&self.geo, ci);
        edge_basis(&lam, &g, &self.cd.sign[ci], self.order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (mut es0, mut es1) = (c(0.0, 0.0), c(0.0, 0.0));
        for s in 0..nloc_e {
            let cf = self.et[self.cd.gdofs[ci][s]];
            es0 += cf * w[s][0];
            es1 += cf * w[s][1];
        }
        let (_, dn) = shape(&lam, self.order);
        let mut pd = Vec::with_capacity(6);
        self.mesh.cell_dofs(ci, self.order, &mut pd);
        let (mut dp0, mut dp1) = (c(0.0, 0.0), c(0.0, 0.0));
        for a in 0..nl {
            let gn0: f64 = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
            let gn1: f64 = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            let pv = self.phi[pd[a]];
            dp0 += pv * gn0;
            dp1 += pv * gn1;
        }
        let j_beta = c(0.0, 1.0) * self.beta;
        let inv_mu = 1.0 / self.mu_r[ci].max(1e-12);
        let h_pref = c(0.0, -inv_mu / (self.omega * MU0));
        let hs0 = h_pref * (j_beta * (dp1 + es1));
        let hs1 = h_pref * (-j_beta * (dp0 + es0));
        (es0, es1, hs0, hs1)
    }
}
