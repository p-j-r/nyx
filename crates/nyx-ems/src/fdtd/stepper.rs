// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Deterministic SIMD-friendly Structure-of-Arrays (SoA) Yee Maxwell time-marching engine (ADR 145).
//!
//! All slab updates partition outer `x`-slabs into fixed-size chunks via `femcore::par`, guaranteeing
//! 0-ULP bitwise-identical field histories and DFT monitors across `NYX_THREADS=1` vs `8`.

use super::ade_pml::{
    apply_ade_correction, AdeAndAmpereCoeffs, ComponentAdeBank, CpmlState, FdtdBcFace, FdtdBoundaries,
    NonlinearEdgeEntry, NonlinearFdtdBank, NonlinearFdtdDomain,
};
use super::grid::{FdtdDim, Pol2d, YeeFields, YeeGrid, C0, EPS0};
use super::modes_bridge::rh_uv_axes;
use super::sources_monitors::{Aux1dYeeLine, CompiledSource, PlanarDftAccumulator, VolumeFieldDft};
use super::subpixel::{OffDiagEntry, YeeMaterialCoeffs};

/// Fixed slab chunk size for deterministic `femcore::par` multi-threading (`ADR 11`).
const SLAB_CHUNK: usize = 4;

/// Advances one complete leapfrog Yee time step `n -> n + 1`:
/// 1. `H^{n - 1/2} -> H^{n + 1/2}` (Faraday curl + CPML `psi_H` + TF/SF magnetic corrections + Bloch/Periodic wraps)
/// 2. `E^n -> E^{n + 1}` (Ampere curl + CPML `psi_E` + ADE Drude/Lorentz poles + anisotropic off-diagonal `W` + TF/SF & dipole sources + PEC/PMC/Bloch BCs)
/// 3. Accumulates all running DFT and time-series monitors at `t_e = (n + 1) * dt`, `t_h = (n + 0.5) * dt`.
#[allow(clippy::too_many_arguments)]
pub fn step_yee_once(
    step_idx: usize,
    dt: f64,
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    mats: &YeeMaterialCoeffs,
    coeffs: &mut AdeAndAmpereCoeffs,
    cpml_re: &mut CpmlState,
    mut cpml_im: Option<&mut CpmlState>,
    ade_im: Option<&mut [ComponentAdeBank; 3]>,
    fields_re: &mut YeeFields,
    mut fields_im: Option<&mut YeeFields>,
    sources: &mut [CompiledSource],
    planar_dfts: &mut [PlanarDftAccumulator],
    vol_dfts: &mut [VolumeFieldDft],
) {
    let t_n = (step_idx as f64) * dt;
    let t_half = (step_idx as f64 + 0.5) * dt;
    let t_np1 = (step_idx as f64 + 1.0) * dt;

    // 1. Update H^{n - 1/2} -> H^{n + 1/2} on real (and optional Bloch imaginary) grids
    update_h_slab(grid, mats, coeffs, cpml_re, fields_re);
    if let (Some(cp_i), Some(f_i)) = (cpml_im.as_deref_mut(), fields_im.as_deref_mut()) {
        update_h_slab(grid, mats, coeffs, cp_i, f_i);
    }

    // Step 1D auxiliary lines for PlaneWave / TFSF sources and apply TF/SF magnetic corrections at `t_n`
    for src in sources.iter_mut() {
        apply_source_h(src, grid, coeffs, fields_re, t_n, dt);
    }

    // 2. Save dispersive `E^n` (and nonlinear `E^n` if active) before Ampere update, advance `E^n -> E^{n + 1}`, and apply ADE + off-diagonal + sources + nonlinear polarization
    let saved_re = save_dispersive_e(&coeffs.ade_x, &coeffs.ade_y, &coeffs.ade_z, fields_re);
    let saved_nl_re = coeffs.nonlinear.is_active().then(|| save_nonlinear_e(&coeffs.nonlinear, fields_re));
    let saved_im = match (&ade_im, fields_im.as_deref()) {
        (Some(banks), Some(f_i)) => Some(save_dispersive_e(&banks[0], &banks[1], &banks[2], f_i)),
        _ => None,
    };

    // Save previous E for off-diagonal polarized/anisotropic coupling (`D_perp ~ eps_0 / inv_eps * E^n`)
    let offdiag_rhs_re = mats.has_offdiag().then(|| compute_offdiag_e_delta(grid, mats, fields_re));
    let offdiag_rhs_im = match (mats.has_offdiag(), fields_im.as_deref()) {
        (true, Some(f_i)) => Some(compute_offdiag_e_delta(grid, mats, f_i)),
        _ => None,
    };

    update_e_slab(grid, bcs, coeffs, cpml_re, fields_re);
    if let (Some(cp_i), Some(f_i)) = (cpml_im.as_deref_mut(), fields_im.as_deref_mut()) {
        update_e_slab(grid, bcs, coeffs, cp_i, f_i);
    }
    super::apply_periodic_e_updates(grid, bcs, coeffs, cpml_re, cpml_im, fields_re, fields_im.as_deref_mut());

    if !coeffs.ade_x.is_empty() {
        apply_ade_correction(&mut coeffs.ade_x, fields_re.ex.as_mut_slice(), &saved_re[0], &coeffs.cb_x, dt);
    }
    if !coeffs.ade_y.is_empty() {
        apply_ade_correction(&mut coeffs.ade_y, fields_re.ey.as_mut_slice(), &saved_re[1], &coeffs.cb_y, dt);
    }
    if !coeffs.ade_z.is_empty() {
        apply_ade_correction(&mut coeffs.ade_z, fields_re.ez.as_mut_slice(), &saved_re[2], &coeffs.cb_z, dt);
    }
    if let (Some(banks), Some(f_i), Some(sv_i)) = (ade_im, fields_im.as_deref_mut(), saved_im) {
        if !banks[0].is_empty() {
            apply_ade_correction(&mut banks[0], f_i.ex.as_mut_slice(), &sv_i[0], &coeffs.cb_x, dt);
        }
        if !banks[1].is_empty() {
            apply_ade_correction(&mut banks[1], f_i.ey.as_mut_slice(), &sv_i[1], &coeffs.cb_y, dt);
        }
        if !banks[2].is_empty() {
            apply_ade_correction(&mut banks[2], f_i.ez.as_mut_slice(), &sv_i[2], &coeffs.cb_z, dt);
        }
    }

    // Apply off-diagonal subpixel / anisotropic delta correction (`Delta E_a = W_{ab} (D_b^{n+1} - D_b^n) + W_{ac} (D_c^{n+1} - D_c^n)`)
    if let Some(prev_d) = offdiag_rhs_re {
        apply_offdiag_e_delta(grid, mats, fields_re, &prev_d);
    }
    if let (Some(prev_di), Some(f_i)) = (offdiag_rhs_im, fields_im.as_deref_mut()) {
        apply_offdiag_e_delta(grid, mats, f_i, &prev_di);
    }

    // Inject TF/SF electric corrections and soft dipole sources at `t_{n + 1/2}`
    for src in sources.iter_mut() {
        apply_source_e(src, grid, coeffs, fields_re, t_half, dt);
    }

    // Apply instantaneous chi^(2)/QPM/chi^(3) Kerr + TPA/FCA/FCD nonlinear polarization update (`ADR 143`)
    if let Some(ref e_old_nl) = saved_nl_re {
        super::sync_periodic_e_high_faces(grid, bcs, fields_re, fields_im.as_deref_mut());
        enforce_e_pec_clamp(grid, bcs, mats, fields_re, fields_im.as_deref_mut());
        apply_nonlinear_e_update(grid, bcs, mats, coeffs, fields_re, e_old_nl, dt);
    }

    // Synchronize periodic / Bloch high faces (`i = n`) and clamp PEC boundaries before DFT sampling
    super::sync_periodic_e_high_faces(grid, bcs, fields_re, fields_im.as_deref_mut());
    enforce_e_pec_clamp(grid, bcs, mats, fields_re, fields_im.as_deref_mut());

    // 3. Accumulate running DFT monitors at `t_e = t_{n+1}`, `t_h = t_{n+1/2}`
    let f_im_ref = fields_im.as_deref();
    for mon in planar_dfts.iter_mut() {
        mon.accumulate(grid, fields_re, f_im_ref, t_np1, t_half, dt);
    }
    for vmon in vol_dfts.iter_mut() {
        vmon.accumulate(grid, fields_re, f_im_ref, t_np1, t_half, dt);
    }
}

/// Computes discrete electromagnetic field energy proxy `sum (E^2 + eta_0^2 H^2)` for auto-shutoff detection,
/// reduced in deterministic fixed-chunk order via `femcore::par::sum_range`.
pub fn total_field_energy_proxy(grid: &YeeGrid, fields_re: &YeeFields, fields_im: Option<&YeeFields>) -> f64 {
    let n = grid.slab_len();
    let ex = fields_re.ex.as_slice();
    let ey = fields_re.ey.as_slice();
    let ez = fields_re.ez.as_slice();
    let hx = fields_re.hx.as_slice();
    let hy = fields_re.hy.as_slice();
    let hz = fields_re.hz.as_slice();
    let eta0_sq = 376.730_313_668_f64 * 376.730_313_668_f64;
    let nt = femcore::par::threads();
    femcore::par::sum_range(n, 2048, nt, |g| {
        let mut s =
            ex[g] * ex[g] + ey[g] * ey[g] + ez[g] * ez[g] + eta0_sq * (hx[g] * hx[g] + hy[g] * hy[g] + hz[g] * hz[g]);
        if let Some(fi) = fields_im {
            let exi = fi.ex.as_slice()[g];
            let eyi = fi.ey.as_slice()[g];
            let ezi = fi.ez.as_slice()[g];
            let hxi = fi.hx.as_slice()[g];
            let hyi = fi.hy.as_slice()[g];
            let hzi = fi.hz.as_slice()[g];
            s += exi * exi + eyi * eyi + ezi * ezi + eta0_sq * (hxi * hxi + hyi * hyi + hzi * hzi);
        }
        s
    })
}

// ------------------------------------------------------------------ Inner SIMD Slab Kernels (`H` & `E`)

fn update_h_slab(
    grid: &YeeGrid,
    mats: &YeeMaterialCoeffs,
    coeffs: &AdeAndAmpereCoeffs,
    cpml: &mut CpmlState,
    fields: &mut YeeFields,
) {
    let (nx, ny, nz, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.stride_x, grid.stride_y);
    let slab_elems = sx;
    let ex = fields.ex.as_slice();
    let ey = fields.ey.as_slice();
    let ez = fields.ez.as_slice();
    let pec_lx = mats.pec_len_x.as_slice();
    let pec_ly = mats.pec_len_y.as_slice();
    let pec_lz = mats.pec_len_z.as_slice();
    let use_conf = mats.has_conformal_pec;

    let has_cpml_x = cpml.x.has_cpml();
    let has_cpml_y = cpml.y.has_cpml();
    let has_cpml_z = cpml.z.has_cpml();

    let run_par = femcore::par::threads() > 1 && grid.slab_len() >= 16_384 && nx >= 8;

    match grid.dim {
        FdtdDim::TwoD { polarization: Pol2d::TmEy } => {
            // Only Hx (ix in 0..=nx, iz in 0..nz) and Hz (ix in 0..nx, iz in 0..=nz) depend on Ey!
            // Hx: dEy_dz = (ey[ix, 0, iz + 1] - ey[ix, 0, iz]) * inv_dp_z[iz]; curl_ex = -dEy_dz => Hx = da_x Hx + db_x dEy_dz_s
            {
                let hx_all = fields.hx.as_mut_slice();
                let chunk_sz = SLAB_CHUNK * slab_elems;
                if has_cpml_z {
                    let psi_z = &mut cpml.psi_hx_z;
                    for ix in 0..=nx {
                        let base = ix * sx;
                        for iz in 0..nz {
                            let g = base + iz;
                            let dey_dz = (ey[g + 1] - ey[g]) * grid.inv_dp_z[iz];
                            let dey_s = if cpml.z.active_h[iz] {
                                psi_z[g] = cpml.z.b_h[iz] * psi_z[g] + cpml.z.c_h[iz] * dey_dz;
                                cpml.z.inv_kappa_h[iz] * dey_dz + psi_z[g]
                            } else {
                                dey_dz
                            };
                            hx_all[g] = coeffs.da_x[g] * hx_all[g] + coeffs.db_x[g] * dey_s;
                        }
                    }
                } else {
                    femcore::par::for_chunks_mut(hx_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min((nx + 1).saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            let g_base = ix * sx;
                            let c_base = loc_x * slab_elems;
                            for iz in 0..nz {
                                let g = g_base + iz;
                                let dey_dz = (ey[g + 1] - ey[g]) * grid.inv_dp_z[iz];
                                chunk[c_base + iz] = coeffs.da_x[g] * chunk[c_base + iz] + coeffs.db_x[g] * dey_dz;
                            }
                        }
                    });
                }
            }
            // Hz: dEy_dx = (ey[ix + 1, 0, iz] - ey[ix, 0, iz]) * inv_dp_x[ix]; curl_ez = +dEy_dx => Hz = da_z Hz - db_z dEy_dx_s
            {
                let hz_all = fields.hz.as_mut_slice();
                let chunk_sz = SLAB_CHUNK * slab_elems;
                if has_cpml_x {
                    let psi_x = &mut cpml.psi_hz_x;
                    for ix in 0..nx {
                        let base = ix * sx;
                        let inv_dx = grid.inv_dp_x[ix];
                        let act_x = cpml.x.active_h[ix];
                        let (bx, cx_p, ikx) = (cpml.x.b_h[ix], cpml.x.c_h[ix], cpml.x.inv_kappa_h[ix]);
                        for iz in 0..=nz {
                            let g = base + iz;
                            let dey_dx = (ey[g + sx] - ey[g]) * inv_dx;
                            let dey_s = if act_x {
                                psi_x[g] = bx * psi_x[g] + cx_p * dey_dx;
                                ikx * dey_dx + psi_x[g]
                            } else {
                                dey_dx
                            };
                            hz_all[g] = coeffs.da_z[g] * hz_all[g] - coeffs.db_z[g] * dey_s;
                        }
                    }
                } else {
                    femcore::par::for_chunks_mut(hz_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            let g_base = ix * sx;
                            let c_base = loc_x * slab_elems;
                            let inv_dx = grid.inv_dp_x[ix];
                            for iz in 0..=nz {
                                let g = g_base + iz;
                                let dey_dx = (ey[g + sx] - ey[g]) * inv_dx;
                                chunk[c_base + iz] = coeffs.da_z[g] * chunk[c_base + iz] - coeffs.db_z[g] * dey_dx;
                            }
                        }
                    });
                }
            }
        }
        FdtdDim::TwoD { polarization: Pol2d::TeHy } => {
            // Only Hy (ix in 0..nx, iz in 0..nz) depends on (Ex, Ez):
            // curl_ey = dEx_dz - dEz_dx => Hy = da_y Hy - db_y (dEx_dz_s - dEz_dx_s)
            let hy_all = fields.hy.as_mut_slice();
            if has_cpml_x || has_cpml_z {
                let psi_x = &mut cpml.psi_hy_x;
                let psi_z = &mut cpml.psi_hy_z;
                for ix in 0..nx {
                    let base = ix * sx;
                    let inv_dx = grid.inv_dp_x[ix];
                    let act_x = has_cpml_x && cpml.x.active_h[ix];
                    for iz in 0..nz {
                        let g = base + iz;
                        let (ex0, ex1, ez0, ez1) = if use_conf {
                            (
                                ex[g] * pec_lx[g],
                                ex[g + 1] * pec_lx[g + 1],
                                ez[g] * pec_lz[g],
                                ez[g + sx] * pec_lz[g + sx],
                            )
                        } else {
                            (ex[g], ex[g + 1], ez[g], ez[g + sx])
                        };
                        let dex_dz = (ex1 - ex0) * grid.inv_dp_z[iz];
                        let dez_dx = (ez1 - ez0) * inv_dx;
                        let dex_s = if has_cpml_z && cpml.z.active_h[iz] {
                            psi_z[g] = cpml.z.b_h[iz] * psi_z[g] + cpml.z.c_h[iz] * dex_dz;
                            cpml.z.inv_kappa_h[iz] * dex_dz + psi_z[g]
                        } else {
                            dex_dz
                        };
                        let dez_s = if act_x {
                            psi_x[g] = cpml.x.b_h[ix] * psi_x[g] + cpml.x.c_h[ix] * dez_dx;
                            cpml.x.inv_kappa_h[ix] * dez_dx + psi_x[g]
                        } else {
                            dez_dx
                        };
                        hy_all[g] = coeffs.da_y[g] * hy_all[g] - coeffs.db_y[g] * (dex_s - dez_s);
                    }
                }
            } else {
                let chunk_sz = SLAB_CHUNK * slab_elems;
                femcore::par::for_chunks_mut(hy_all, chunk_sz, run_par, |ci, chunk| {
                    let ix0 = ci * SLAB_CHUNK;
                    let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                    for loc_x in 0..n_ix {
                        let ix = ix0 + loc_x;
                        let g_base = ix * sx;
                        let c_base = loc_x * slab_elems;
                        let inv_dx = grid.inv_dp_x[ix];
                        for iz in 0..nz {
                            let g = g_base + iz;
                            let dex_dz = (ex[g + 1] - ex[g]) * grid.inv_dp_z[iz];
                            let dez_dx = (ez[g + sx] - ez[g]) * inv_dx;
                            chunk[c_base + iz] =
                                coeffs.da_y[g] * chunk[c_base + iz] - coeffs.db_y[g] * (dex_dz - dez_dx);
                        }
                    }
                });
            }
        }
        FdtdDim::TwoD { polarization: Pol2d::Full } | FdtdDim::ThreeD => {
            let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
            // Full 3D (or 2.5D) Hx, Hy, Hz update!
            // 1. Hx (`ix in 0..=nx, iy in 0..ny, iz in 0..nz`)
            {
                let hx_all = fields.hx.as_mut_slice();
                let chunk_sz = SLAB_CHUNK * slab_elems;
                if has_cpml_y || has_cpml_z {
                    let psi_y = &mut cpml.psi_hx_y;
                    let psi_z = &mut cpml.psi_hx_z;
                    for ix in 0..=nx {
                        for iy in 0..ny {
                            let base = ix * sx + iy * sy;
                            let inv_dy = grid.inv_dp_y[iy];
                            let act_y = !is_2d && has_cpml_y && cpml.y.active_h[iy];
                            for iz in 0..nz {
                                let g = base + iz;
                                let dez_dy = if is_2d {
                                    0.0
                                } else {
                                    let d = if use_conf {
                                        ez[g + sy] * pec_lz[g + sy] - ez[g] * pec_lz[g]
                                    } else {
                                        ez[g + sy] - ez[g]
                                    } * inv_dy;
                                    if act_y {
                                        psi_y[g] = cpml.y.b_h[iy] * psi_y[g] + cpml.y.c_h[iy] * d;
                                        cpml.y.inv_kappa_h[iy] * d + psi_y[g]
                                    } else {
                                        d
                                    }
                                };
                                let d_raw = if use_conf {
                                    ey[g + 1] * pec_ly[g + 1] - ey[g] * pec_ly[g]
                                } else {
                                    ey[g + 1] - ey[g]
                                } * grid.inv_dp_z[iz];
                                let dey_dz = if has_cpml_z && cpml.z.active_h[iz] {
                                    psi_z[g] = cpml.z.b_h[iz] * psi_z[g] + cpml.z.c_h[iz] * d_raw;
                                    cpml.z.inv_kappa_h[iz] * d_raw + psi_z[g]
                                } else {
                                    d_raw
                                };
                                hx_all[g] = coeffs.da_x[g] * hx_all[g] - coeffs.db_x[g] * (dez_dy - dey_dz);
                            }
                        }
                    }
                } else {
                    femcore::par::for_chunks_mut(hx_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min((nx + 1).saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            for iy in 0..ny {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                let inv_dy = grid.inv_dp_y[iy];
                                for iz in 0..nz {
                                    let g = g_base + iz;
                                    let dez_dy = if is_2d { 0.0 } else { (ez[g + sy] - ez[g]) * inv_dy };
                                    let dey_dz = (ey[g + 1] - ey[g]) * grid.inv_dp_z[iz];
                                    chunk[c_base + iz] =
                                        coeffs.da_x[g] * chunk[c_base + iz] - coeffs.db_x[g] * (dez_dy - dey_dz);
                                }
                            }
                        }
                    });
                }
            }

            // 2. Hy (`ix in 0..nx, iy in 0..=ny, iz in 0..nz`)
            {
                let hy_all = fields.hy.as_mut_slice();
                let iy_hi = if is_2d { 0 } else { ny };
                let chunk_sz = SLAB_CHUNK * slab_elems;
                if has_cpml_x || has_cpml_z {
                    let psi_x = &mut cpml.psi_hy_x;
                    let psi_z = &mut cpml.psi_hy_z;
                    for ix in 0..nx {
                        let inv_dx = grid.inv_dp_x[ix];
                        let act_x = has_cpml_x && cpml.x.active_h[ix];
                        for iy in 0..=iy_hi {
                            let base = ix * sx + iy * sy;
                            for iz in 0..nz {
                                let g = base + iz;
                                let dex_raw = if use_conf {
                                    ex[g + 1] * pec_lx[g + 1] - ex[g] * pec_lx[g]
                                } else {
                                    ex[g + 1] - ex[g]
                                } * grid.inv_dp_z[iz];
                                let dez_raw = if use_conf {
                                    ez[g + sx] * pec_lz[g + sx] - ez[g] * pec_lz[g]
                                } else {
                                    ez[g + sx] - ez[g]
                                } * inv_dx;
                                let dex_dz = if has_cpml_z && cpml.z.active_h[iz] {
                                    psi_z[g] = cpml.z.b_h[iz] * psi_z[g] + cpml.z.c_h[iz] * dex_raw;
                                    cpml.z.inv_kappa_h[iz] * dex_raw + psi_z[g]
                                } else {
                                    dex_raw
                                };
                                let dez_dx = if act_x {
                                    psi_x[g] = cpml.x.b_h[ix] * psi_x[g] + cpml.x.c_h[ix] * dez_raw;
                                    cpml.x.inv_kappa_h[ix] * dez_raw + psi_x[g]
                                } else {
                                    dez_raw
                                };
                                hy_all[g] = coeffs.da_y[g] * hy_all[g] - coeffs.db_y[g] * (dex_dz - dez_dx);
                            }
                        }
                    }
                } else {
                    femcore::par::for_chunks_mut(hy_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            let inv_dx = grid.inv_dp_x[ix];
                            for iy in 0..=iy_hi {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                for iz in 0..nz {
                                    let g = g_base + iz;
                                    let dex_dz = (ex[g + 1] - ex[g]) * grid.inv_dp_z[iz];
                                    let dez_dx = (ez[g + sx] - ez[g]) * inv_dx;
                                    chunk[c_base + iz] =
                                        coeffs.da_y[g] * chunk[c_base + iz] - coeffs.db_y[g] * (dex_dz - dez_dx);
                                }
                            }
                        }
                    });
                }
            }

            // 3. Hz (`ix in 0..nx, iy in 0..ny, iz in 0..=nz`)
            {
                let hz_all = fields.hz.as_mut_slice();
                let chunk_sz = SLAB_CHUNK * slab_elems;
                if has_cpml_x || has_cpml_y {
                    let psi_x = &mut cpml.psi_hz_x;
                    let psi_y = &mut cpml.psi_hz_y;
                    for ix in 0..nx {
                        let inv_dx = grid.inv_dp_x[ix];
                        let act_x = has_cpml_x && cpml.x.active_h[ix];
                        for iy in 0..ny {
                            let base = ix * sx + iy * sy;
                            let inv_dy = grid.inv_dp_y[iy];
                            let act_y = !is_2d && has_cpml_y && cpml.y.active_h[iy];
                            for iz in 0..=nz {
                                let g = base + iz;
                                let dey_raw = if use_conf {
                                    ey[g + sx] * pec_ly[g + sx] - ey[g] * pec_ly[g]
                                } else {
                                    ey[g + sx] - ey[g]
                                } * inv_dx;
                                let dey_dx = if act_x {
                                    psi_x[g] = cpml.x.b_h[ix] * psi_x[g] + cpml.x.c_h[ix] * dey_raw;
                                    cpml.x.inv_kappa_h[ix] * dey_raw + psi_x[g]
                                } else {
                                    dey_raw
                                };
                                let dex_dy = if is_2d {
                                    0.0
                                } else {
                                    let d = if use_conf {
                                        ex[g + sy] * pec_lx[g + sy] - ex[g] * pec_lx[g]
                                    } else {
                                        ex[g + sy] - ex[g]
                                    } * inv_dy;
                                    if act_y {
                                        psi_y[g] = cpml.y.b_h[iy] * psi_y[g] + cpml.y.c_h[iy] * d;
                                        cpml.y.inv_kappa_h[iy] * d + psi_y[g]
                                    } else {
                                        d
                                    }
                                };
                                hz_all[g] = coeffs.da_z[g] * hz_all[g] - coeffs.db_z[g] * (dey_dx - dex_dy);
                            }
                        }
                    }
                } else {
                    femcore::par::for_chunks_mut(hz_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            let inv_dx = grid.inv_dp_x[ix];
                            for iy in 0..ny {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                let inv_dy = grid.inv_dp_y[iy];
                                for iz in 0..=nz {
                                    let g = g_base + iz;
                                    let dey_dx = (ey[g + sx] - ey[g]) * inv_dx;
                                    let dex_dy = if is_2d { 0.0 } else { (ex[g + sy] - ex[g]) * inv_dy };
                                    chunk[c_base + iz] =
                                        coeffs.da_z[g] * chunk[c_base + iz] - coeffs.db_z[g] * (dey_dx - dex_dy);
                                }
                            }
                        }
                    });
                }
            }
        }
    }
}

fn update_e_slab(
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    coeffs: &AdeAndAmpereCoeffs,
    cpml: &mut CpmlState,
    fields: &mut YeeFields,
) {
    let (nx, ny, nz, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.stride_x, grid.stride_y);
    let slab_elems = sx;
    let hx = fields.hx.as_slice();
    let hy = fields.hy.as_slice();
    let hz = fields.hz.as_slice();

    let has_cpml_x = cpml.x.has_cpml();
    let has_cpml_y = cpml.y.has_cpml();
    let has_cpml_z = cpml.z.has_cpml();
    let pmc_xm = matches!(bcs.faces[0], FdtdBcFace::Pmc);
    let pmc_xp = matches!(bcs.faces[1], FdtdBcFace::Pmc);
    let pmc_ym = matches!(bcs.faces[2], FdtdBcFace::Pmc);
    let pmc_yp = matches!(bcs.faces[3], FdtdBcFace::Pmc);
    let pmc_zm = matches!(bcs.faces[4], FdtdBcFace::Pmc);
    let pmc_zp = matches!(bcs.faces[5], FdtdBcFace::Pmc);

    let run_par = femcore::par::threads() > 1 && grid.slab_len() >= 16_384 && nx >= 8;

    match grid.dim {
        FdtdDim::TwoD { polarization: Pol2d::TmEy } => {
            // Ey(ix, 0, iz): curl_hy = dHx_dz - dHz_dx
            let ey_all = fields.ey.as_mut_slice();
            let ix_lo = if pmc_xm { 0 } else { 1 };
            let ix_hi = if pmc_xp { nx } else { nx - 1 };
            let iz_lo = if pmc_zm { 0 } else { 1 };
            let iz_hi = if pmc_zp { nz } else { nz - 1 };

            if has_cpml_x || has_cpml_z || pmc_xm || pmc_xp || pmc_zm || pmc_zp {
                let psi_x = &mut cpml.psi_ey_x;
                let psi_z = &mut cpml.psi_ey_z;
                for ix in ix_lo..=ix_hi {
                    let base = ix * sx;
                    let inv_ddx = grid.inv_dd_x[ix];
                    let act_x = has_cpml_x && cpml.x.active_e[ix];
                    for iz in iz_lo..=iz_hi {
                        let g = base + iz;
                        let dhx_raw = if iz == 0 {
                            2.0 * hx[g] * grid.inv_dd_z[0]
                        } else if iz == nz {
                            -2.0 * hx[g - 1] * grid.inv_dd_z[nz]
                        } else {
                            (hx[g] - hx[g - 1]) * grid.inv_dd_z[iz]
                        };
                        let dhz_raw = if ix == 0 {
                            2.0 * hz[g] * inv_ddx
                        } else if ix == nx {
                            -2.0 * hz[g - sx] * inv_ddx
                        } else {
                            (hz[g] - hz[g - sx]) * inv_ddx
                        };
                        let dhx_dz = if has_cpml_z && cpml.z.active_e[iz] {
                            psi_z[g] = cpml.z.b_e[iz] * psi_z[g] + cpml.z.c_e[iz] * dhx_raw;
                            cpml.z.inv_kappa_e[iz] * dhx_raw + psi_z[g]
                        } else {
                            dhx_raw
                        };
                        let dhz_dx = if act_x {
                            psi_x[g] = cpml.x.b_e[ix] * psi_x[g] + cpml.x.c_e[ix] * dhz_raw;
                            cpml.x.inv_kappa_e[ix] * dhz_raw + psi_x[g]
                        } else {
                            dhz_raw
                        };
                        ey_all[g] = coeffs.ca_y[g] * ey_all[g] + coeffs.cb_y[g] * (dhx_dz - dhz_dx);
                    }
                }
            } else {
                let chunk_sz = SLAB_CHUNK * slab_elems;
                femcore::par::for_chunks_mut(ey_all, chunk_sz, run_par, |ci, chunk| {
                    let ix0 = ci * SLAB_CHUNK;
                    let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                    for loc_x in 0..n_ix {
                        let ix = ix0 + loc_x;
                        if ix == 0 || ix >= nx {
                            continue;
                        }
                        let g_base = ix * sx;
                        let c_base = loc_x * slab_elems;
                        let inv_ddx = grid.inv_dd_x[ix];
                        for iz in 1..nz {
                            let g = g_base + iz;
                            let dhx_dz = (hx[g] - hx[g - 1]) * grid.inv_dd_z[iz];
                            let dhz_dx = (hz[g] - hz[g - sx]) * inv_ddx;
                            chunk[c_base + iz] =
                                coeffs.ca_y[g] * chunk[c_base + iz] + coeffs.cb_y[g] * (dhx_dz - dhz_dx);
                        }
                    }
                });
            }
        }
        FdtdDim::TwoD { polarization: Pol2d::TeHy } => {
            // Ex: curl_hx = -dHy_dz; Ez: curl_hz = +dHy_dx
            {
                let ex_all = fields.ex.as_mut_slice();
                let iz_lo = if pmc_zm { 0 } else { 1 };
                let iz_hi = if pmc_zp { nz } else { nz - 1 };
                let psi_z = &mut cpml.psi_ex_z;
                for ix in 0..nx {
                    let base = ix * sx;
                    for iz in iz_lo..=iz_hi {
                        let g = base + iz;
                        let dhy_raw = if iz == 0 {
                            2.0 * hy[g] * grid.inv_dd_z[0]
                        } else if iz == nz {
                            -2.0 * hy[g - 1] * grid.inv_dd_z[nz]
                        } else {
                            (hy[g] - hy[g - 1]) * grid.inv_dd_z[iz]
                        };
                        let dhy_dz = if has_cpml_z && cpml.z.active_e[iz] {
                            psi_z[g] = cpml.z.b_e[iz] * psi_z[g] + cpml.z.c_e[iz] * dhy_raw;
                            cpml.z.inv_kappa_e[iz] * dhy_raw + psi_z[g]
                        } else {
                            dhy_raw
                        };
                        ex_all[g] = coeffs.ca_x[g] * ex_all[g] - coeffs.cb_x[g] * dhy_dz;
                    }
                }
            }
            {
                let ez_all = fields.ez.as_mut_slice();
                let ix_lo = if pmc_xm { 0 } else { 1 };
                let ix_hi = if pmc_xp { nx } else { nx - 1 };
                let psi_x = &mut cpml.psi_ez_x;
                for ix in ix_lo..=ix_hi {
                    let base = ix * sx;
                    let inv_ddx = grid.inv_dd_x[ix];
                    let act_x = has_cpml_x && cpml.x.active_e[ix];
                    for iz in 0..nz {
                        let g = base + iz;
                        let dhy_raw = if ix == 0 {
                            2.0 * hy[g] * inv_ddx
                        } else if ix == nx {
                            -2.0 * hy[g - sx] * inv_ddx
                        } else {
                            (hy[g] - hy[g - sx]) * inv_ddx
                        };
                        let dhy_dx = if act_x {
                            psi_x[g] = cpml.x.b_e[ix] * psi_x[g] + cpml.x.c_e[ix] * dhy_raw;
                            cpml.x.inv_kappa_e[ix] * dhy_raw + psi_x[g]
                        } else {
                            dhy_raw
                        };
                        ez_all[g] = coeffs.ca_z[g] * ez_all[g] + coeffs.cb_z[g] * dhy_dx;
                    }
                }
            }
        }
        FdtdDim::TwoD { polarization: Pol2d::Full } | FdtdDim::ThreeD => {
            let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
            let has_pmc = pmc_xm || pmc_xp || pmc_ym || pmc_yp || pmc_zm || pmc_zp;

            // 1. Ex (`ix in 0..nx`, `iy in iy_lo..=iy_hi`, `iz in iz_lo..=iz_hi`): curl_hx = dHz_dy - dHy_dz
            {
                let ex_all = fields.ex.as_mut_slice();
                if has_cpml_y || has_cpml_z || has_pmc || is_2d {
                    let iy_lo = if is_2d || pmc_ym { 0 } else { 1 };
                    let iy_hi = if is_2d {
                        0
                    } else if pmc_yp {
                        ny
                    } else {
                        ny - 1
                    };
                    let iz_lo = if pmc_zm { 0 } else { 1 };
                    let iz_hi = if pmc_zp { nz } else { nz - 1 };
                    let psi_y = &mut cpml.psi_ex_y;
                    let psi_z = &mut cpml.psi_ex_z;
                    for ix in 0..nx {
                        for iy in iy_lo..=iy_hi {
                            let base = ix * sx + iy * sy;
                            let inv_ddy = grid.inv_dd_y[iy];
                            let act_y = !is_2d && has_cpml_y && cpml.y.active_e[iy];
                            for iz in iz_lo..=iz_hi {
                                let g = base + iz;
                                let dhz_dy = if is_2d {
                                    0.0
                                } else {
                                    let d = if iy == 0 {
                                        2.0 * hz[g] * inv_ddy
                                    } else if iy == ny {
                                        -2.0 * hz[g - sy] * inv_ddy
                                    } else {
                                        (hz[g] - hz[g - sy]) * inv_ddy
                                    };
                                    if act_y {
                                        psi_y[g] = cpml.y.b_e[iy] * psi_y[g] + cpml.y.c_e[iy] * d;
                                        cpml.y.inv_kappa_e[iy] * d + psi_y[g]
                                    } else {
                                        d
                                    }
                                };
                                let dhy_raw = if iz == 0 {
                                    2.0 * hy[g] * grid.inv_dd_z[0]
                                } else if iz == nz {
                                    -2.0 * hy[g - 1] * grid.inv_dd_z[nz]
                                } else {
                                    (hy[g] - hy[g - 1]) * grid.inv_dd_z[iz]
                                };
                                let dhy_dz = if has_cpml_z && cpml.z.active_e[iz] {
                                    psi_z[g] = cpml.z.b_e[iz] * psi_z[g] + cpml.z.c_e[iz] * dhy_raw;
                                    cpml.z.inv_kappa_e[iz] * dhy_raw + psi_z[g]
                                } else {
                                    dhy_raw
                                };
                                ex_all[g] = coeffs.ca_x[g] * ex_all[g] + coeffs.cb_x[g] * (dhz_dy - dhy_dz);
                            }
                        }
                    }
                } else {
                    let chunk_sz = SLAB_CHUNK * slab_elems;
                    femcore::par::for_chunks_mut(ex_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            for iy in 1..ny {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                let inv_ddy = grid.inv_dd_y[iy];
                                for iz in 1..nz {
                                    let g = g_base + iz;
                                    let dhz_dy = (hz[g] - hz[g - sy]) * inv_ddy;
                                    let dhy_dz = (hy[g] - hy[g - 1]) * grid.inv_dd_z[iz];
                                    chunk[c_base + iz] =
                                        coeffs.ca_x[g] * chunk[c_base + iz] + coeffs.cb_x[g] * (dhz_dy - dhy_dz);
                                }
                            }
                        }
                    });
                }
            }

            // 2. Ey (`ix in ix_lo..=ix_hi`, `iy in 0..ny`, `iz in iz_lo..=iz_hi`): curl_hy = dHx_dz - dHz_dx
            {
                let ey_all = fields.ey.as_mut_slice();
                if has_cpml_x || has_cpml_z || has_pmc {
                    let ix_lo = if pmc_xm { 0 } else { 1 };
                    let ix_hi = if pmc_xp { nx } else { nx - 1 };
                    let iz_lo = if pmc_zm { 0 } else { 1 };
                    let iz_hi = if pmc_zp { nz } else { nz - 1 };
                    let psi_x = &mut cpml.psi_ey_x;
                    let psi_z = &mut cpml.psi_ey_z;
                    for ix in ix_lo..=ix_hi {
                        let inv_ddx = grid.inv_dd_x[ix];
                        let act_x = has_cpml_x && cpml.x.active_e[ix];
                        for iy in 0..ny {
                            let base = ix * sx + iy * sy;
                            for iz in iz_lo..=iz_hi {
                                let g = base + iz;
                                let dhx_raw = if iz == 0 {
                                    2.0 * hx[g] * grid.inv_dd_z[0]
                                } else if iz == nz {
                                    -2.0 * hx[g - 1] * grid.inv_dd_z[nz]
                                } else {
                                    (hx[g] - hx[g - 1]) * grid.inv_dd_z[iz]
                                };
                                let dhz_raw = if ix == 0 {
                                    2.0 * hz[g] * inv_ddx
                                } else if ix == nx {
                                    -2.0 * hz[g - sx] * inv_ddx
                                } else {
                                    (hz[g] - hz[g - sx]) * inv_ddx
                                };
                                let dhx_dz = if has_cpml_z && cpml.z.active_e[iz] {
                                    psi_z[g] = cpml.z.b_e[iz] * psi_z[g] + cpml.z.c_e[iz] * dhx_raw;
                                    cpml.z.inv_kappa_e[iz] * dhx_raw + psi_z[g]
                                } else {
                                    dhx_raw
                                };
                                let dhz_dx = if act_x {
                                    psi_x[g] = cpml.x.b_e[ix] * psi_x[g] + cpml.x.c_e[ix] * dhz_raw;
                                    cpml.x.inv_kappa_e[ix] * dhz_raw + psi_x[g]
                                } else {
                                    dhz_raw
                                };
                                ey_all[g] = coeffs.ca_y[g] * ey_all[g] + coeffs.cb_y[g] * (dhx_dz - dhz_dx);
                            }
                        }
                    }
                } else {
                    let chunk_sz = SLAB_CHUNK * slab_elems;
                    femcore::par::for_chunks_mut(ey_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            if ix == 0 {
                                continue;
                            }
                            let inv_ddx = grid.inv_dd_x[ix];
                            for iy in 0..ny {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                for iz in 1..nz {
                                    let g = g_base + iz;
                                    let dhx_dz = (hx[g] - hx[g - 1]) * grid.inv_dd_z[iz];
                                    let dhz_dx = (hz[g] - hz[g - sx]) * inv_ddx;
                                    chunk[c_base + iz] =
                                        coeffs.ca_y[g] * chunk[c_base + iz] + coeffs.cb_y[g] * (dhx_dz - dhz_dx);
                                }
                            }
                        }
                    });
                }
            }

            // 3. Ez (`ix in ix_lo..=ix_hi`, `iy in iy_lo..=iy_hi`, `iz in 0..nz`): curl_hz = dHy_dx - dHx_dy
            {
                let ez_all = fields.ez.as_mut_slice();
                if has_cpml_x || has_cpml_y || has_pmc || is_2d {
                    let ix_lo = if pmc_xm { 0 } else { 1 };
                    let ix_hi = if pmc_xp { nx } else { nx - 1 };
                    let iy_lo = if is_2d || pmc_ym { 0 } else { 1 };
                    let iy_hi = if is_2d {
                        0
                    } else if pmc_yp {
                        ny
                    } else {
                        ny - 1
                    };
                    let psi_x = &mut cpml.psi_ez_x;
                    let psi_y = &mut cpml.psi_ez_y;
                    for ix in ix_lo..=ix_hi {
                        let inv_ddx = grid.inv_dd_x[ix];
                        let act_x = has_cpml_x && cpml.x.active_e[ix];
                        for iy in iy_lo..=iy_hi {
                            let base = ix * sx + iy * sy;
                            let inv_ddy = grid.inv_dd_y[iy];
                            let act_y = !is_2d && has_cpml_y && cpml.y.active_e[iy];
                            for iz in 0..nz {
                                let g = base + iz;
                                let dhy_raw = if ix == 0 {
                                    2.0 * hy[g] * inv_ddx
                                } else if ix == nx {
                                    -2.0 * hy[g - sx] * inv_ddx
                                } else {
                                    (hy[g] - hy[g - sx]) * inv_ddx
                                };
                                let dhy_dx = if act_x {
                                    psi_x[g] = cpml.x.b_e[ix] * psi_x[g] + cpml.x.c_e[ix] * dhy_raw;
                                    cpml.x.inv_kappa_e[ix] * dhy_raw + psi_x[g]
                                } else {
                                    dhy_raw
                                };
                                let dhx_dy = if is_2d {
                                    0.0
                                } else {
                                    let d = if iy == 0 {
                                        2.0 * hx[g] * inv_ddy
                                    } else if iy == ny {
                                        -2.0 * hx[g - sy] * inv_ddy
                                    } else {
                                        (hx[g] - hx[g - sy]) * inv_ddy
                                    };
                                    if act_y {
                                        psi_y[g] = cpml.y.b_e[iy] * psi_y[g] + cpml.y.c_e[iy] * d;
                                        cpml.y.inv_kappa_e[iy] * d + psi_y[g]
                                    } else {
                                        d
                                    }
                                };
                                ez_all[g] = coeffs.ca_z[g] * ez_all[g] + coeffs.cb_z[g] * (dhy_dx - dhx_dy);
                            }
                        }
                    }
                } else {
                    let chunk_sz = SLAB_CHUNK * slab_elems;
                    femcore::par::for_chunks_mut(ez_all, chunk_sz, run_par, |ci, chunk| {
                        let ix0 = ci * SLAB_CHUNK;
                        let n_ix = (chunk.len() / slab_elems).min(nx.saturating_sub(ix0));
                        for loc_x in 0..n_ix {
                            let ix = ix0 + loc_x;
                            if ix == 0 {
                                continue;
                            }
                            let inv_ddx = grid.inv_dd_x[ix];
                            for iy in 1..ny {
                                let g_base = ix * sx + iy * sy;
                                let c_base = loc_x * slab_elems + iy * sy;
                                let inv_ddy = grid.inv_dd_y[iy];
                                for iz in 0..nz {
                                    let g = g_base + iz;
                                    let dhy_dx = (hy[g] - hy[g - sx]) * inv_ddx;
                                    let dhx_dy = (hx[g] - hx[g - sy]) * inv_ddy;
                                    chunk[c_base + iz] =
                                        coeffs.ca_z[g] * chunk[c_base + iz] + coeffs.cb_z[g] * (dhy_dx - dhx_dy);
                                }
                            }
                        }
                    });
                }
            }
        }
    }
}

fn save_dispersive_e(
    bx: &ComponentAdeBank,
    by: &ComponentAdeBank,
    bz: &ComponentAdeBank,
    fields: &YeeFields,
) -> [Vec<f64>; 3] {
    let pick = |bank: &ComponentAdeBank, slab: &[f64]| -> Vec<f64> {
        if bank.is_empty() {
            Vec::new()
        } else {
            bank.edge_indices.iter().map(|&g| slab[g]).collect()
        }
    };
    [pick(bx, fields.ex.as_slice()), pick(by, fields.ey.as_slice()), pick(bz, fields.ez.as_slice())]
}

// ------------------------------------------------------------------ Off-Diagonal Anisotropic / Subpixel Coupling

/// Evaluates `W_{ab} D_{r,b} + W_{ac} D_{r,c}` (`D_r ~ E / inv_eps`) at every off-diagonal interface/crystal edge.
fn compute_offdiag_e_delta(grid: &YeeGrid, mats: &YeeMaterialCoeffs, fields: &YeeFields) -> [Vec<f64>; 3] {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let ex = fields.ex.as_slice();
    let ey = fields.ey.as_slice();
    let ez = fields.ez.as_slice();
    let iex = mats.inv_eps_x.as_slice();
    let iey = mats.inv_eps_y.as_slice();
    let iez = mats.inv_eps_z.as_slice();

    let dr_x = |ix: usize, iy: usize, iz: usize| -> f64 {
        let g = grid.idx(ix, iy, iz);
        if iex[g] > 1e-12 {
            ex[g] / iex[g]
        } else {
            0.0
        }
    };
    let dr_y = |ix: usize, iy: usize, iz: usize| -> f64 {
        let g = grid.idx(ix, iy, iz);
        if iey[g] > 1e-12 {
            ey[g] / iey[g]
        } else {
            0.0
        }
    };
    let dr_z = |ix: usize, iy: usize, iz: usize| -> f64 {
        let g = grid.idx(ix, iy, iz);
        if iez[g] > 1e-12 {
            ez[g] / iez[g]
        } else {
            0.0
        }
    };

    let eval_x = |e: &OffDiagEntry| -> f64 {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let dy_avg = if is_2d {
            0.5 * (dr_y(ix, 0, iz) + dr_y(ix + 1, 0, iz))
        } else {
            let ym = iy.saturating_sub(1);
            let yp = iy.min(grid.ny - 1);
            0.25 * (dr_y(ix, ym, iz) + dr_y(ix + 1, ym, iz) + dr_y(ix, yp, iz) + dr_y(ix + 1, yp, iz))
        };
        let zm = iz.saturating_sub(1);
        let zp = iz.min(grid.nz - 1);
        let dz_avg = 0.25 * (dr_z(ix, iy, zm) + dr_z(ix + 1, iy, zm) + dr_z(ix, iy, zp) + dr_z(ix + 1, iy, zp));
        e.w_ab * dy_avg + e.w_ac * dz_avg
    };

    let eval_y = |e: &OffDiagEntry| -> f64 {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let xm = ix.saturating_sub(1);
        let xp = ix.min(grid.nx - 1);
        let dx_avg = if is_2d {
            0.5 * (dr_x(xm, 0, iz) + dr_x(xp, 0, iz))
        } else {
            0.25 * (dr_x(xm, iy, iz) + dr_x(xp, iy, iz) + dr_x(xm, iy + 1, iz) + dr_x(xp, iy + 1, iz))
        };
        let zm = iz.saturating_sub(1);
        let zp = iz.min(grid.nz - 1);
        let dz_avg = if is_2d {
            0.5 * (dr_z(ix, 0, zm) + dr_z(ix, 0, zp))
        } else {
            0.25 * (dr_z(ix, iy, zm) + dr_z(ix, iy, zp) + dr_z(ix, iy + 1, zm) + dr_z(ix, iy + 1, zp))
        };
        e.w_ab * dx_avg + e.w_ac * dz_avg
    };

    let eval_z = |e: &OffDiagEntry| -> f64 {
        let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
        let xm = ix.saturating_sub(1);
        let xp = ix.min(grid.nx - 1);
        let dx_avg = 0.25 * (dr_x(xm, iy, iz) + dr_x(xp, iy, iz) + dr_x(xm, iy, iz + 1) + dr_x(xp, iy, iz + 1));
        let dy_avg = if is_2d {
            0.5 * (dr_y(ix, 0, iz) + dr_y(ix, 0, iz + 1))
        } else {
            let ym = iy.saturating_sub(1);
            let yp = iy.min(grid.ny - 1);
            0.25 * (dr_y(ix, ym, iz) + dr_y(ix, yp, iz) + dr_y(ix, ym, iz + 1) + dr_y(ix, yp, iz + 1))
        };
        e.w_ab * dx_avg + e.w_ac * dy_avg
    };

    [
        mats.offdiag_x.iter().map(eval_x).collect(),
        mats.offdiag_y.iter().map(eval_y).collect(),
        mats.offdiag_z.iter().map(eval_z).collect(),
    ]
}

fn apply_offdiag_e_delta(grid: &YeeGrid, mats: &YeeMaterialCoeffs, fields: &mut YeeFields, prev: &[Vec<f64>; 3]) {
    let curr = compute_offdiag_e_delta(grid, mats, fields);
    // Simultaneous symmetric Farjadpour/anisotropic update `E^{n+1} - E^n = W (D_r^{n+1} - D_r^n)`
    let ex_mut = fields.ex.as_mut_slice();
    for ((entry, &c_val), &p_val) in mats.offdiag_x.iter().zip(&curr[0]).zip(&prev[0]) {
        ex_mut[entry.idx] += c_val - p_val;
    }
    let ey_mut = fields.ey.as_mut_slice();
    for ((entry, &c_val), &p_val) in mats.offdiag_y.iter().zip(&curr[1]).zip(&prev[1]) {
        ey_mut[entry.idx] += c_val - p_val;
    }
    let ez_mut = fields.ez.as_mut_slice();
    for ((entry, &c_val), &p_val) in mats.offdiag_z.iter().zip(&curr[2]).zip(&prev[2]) {
        ez_mut[entry.idx] += c_val - p_val;
    }
}

// ------------------------------------------------------------------ Nonlinear Constitutive Update (`chi2`, `QPM`, `chi3`, `TPA/FCA/FCD`)

fn save_nonlinear_e(nl: &NonlinearFdtdBank, fields: &YeeFields) -> [Vec<f64>; 3] {
    let pick = |edges: &[NonlinearEdgeEntry], slab: &[f64]| -> Vec<f64> { edges.iter().map(|e| slab[e.idx]).collect() };
    [
        pick(&nl.edges[0], fields.ex.as_slice()),
        pick(&nl.edges[1], fields.ey.as_slice()),
        pick(&nl.edges[2], fields.ez.as_slice()),
    ]
}

#[inline]
fn sample_edge_e_vec(grid: &YeeGrid, fields: &YeeFields, comp: usize, e: &NonlinearEdgeEntry) -> [f64; 3] {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let ex = fields.ex.as_slice();
    let ey = fields.ey.as_slice();
    let ez = fields.ez.as_slice();
    let (ix, iy, iz) = (e.ix as usize, e.iy as usize, e.iz as usize);
    match comp {
        0 => {
            let ex_v = ex[e.idx];
            let ey_v = if is_2d {
                0.5 * (ey[grid.idx(ix, 0, iz)] + ey[grid.idx(ix + 1, 0, iz)])
            } else {
                let ym = iy.saturating_sub(1);
                let yp = iy.min(grid.ny - 1);
                0.25 * (ey[grid.idx(ix, ym, iz)]
                    + ey[grid.idx(ix + 1, ym, iz)]
                    + ey[grid.idx(ix, yp, iz)]
                    + ey[grid.idx(ix + 1, yp, iz)])
            };
            let zm = iz.saturating_sub(1);
            let zp = iz.min(grid.nz - 1);
            let ez_v = 0.25
                * (ez[grid.idx(ix, iy, zm)]
                    + ez[grid.idx(ix + 1, iy, zm)]
                    + ez[grid.idx(ix, iy, zp)]
                    + ez[grid.idx(ix + 1, iy, zp)]);
            [ex_v, ey_v, ez_v]
        }
        1 => {
            let ey_v = ey[e.idx];
            let xm = ix.saturating_sub(1);
            let xp = ix.min(grid.nx - 1);
            let ex_v = if is_2d {
                0.5 * (ex[grid.idx(xm, 0, iz)] + ex[grid.idx(xp, 0, iz)])
            } else {
                0.25 * (ex[grid.idx(xm, iy, iz)]
                    + ex[grid.idx(xp, iy, iz)]
                    + ex[grid.idx(xm, iy + 1, iz)]
                    + ex[grid.idx(xp, iy + 1, iz)])
            };
            let zm = iz.saturating_sub(1);
            let zp = iz.min(grid.nz - 1);
            let ez_v = if is_2d {
                0.5 * (ez[grid.idx(ix, 0, zm)] + ez[grid.idx(ix, 0, zp)])
            } else {
                0.25 * (ez[grid.idx(ix, iy, zm)]
                    + ez[grid.idx(ix, iy, zp)]
                    + ez[grid.idx(ix, iy + 1, zm)]
                    + ez[grid.idx(ix, iy + 1, zp)])
            };
            [ex_v, ey_v, ez_v]
        }
        _ => {
            let ez_v = ez[e.idx];
            let xm = ix.saturating_sub(1);
            let xp = ix.min(grid.nx - 1);
            let ex_v = 0.25
                * (ex[grid.idx(xm, iy, iz)]
                    + ex[grid.idx(xp, iy, iz)]
                    + ex[grid.idx(xm, iy, iz + 1)]
                    + ex[grid.idx(xp, iy, iz + 1)]);
            let ey_v = if is_2d {
                0.5 * (ey[grid.idx(ix, 0, iz)] + ey[grid.idx(ix, 0, iz + 1)])
            } else {
                let ym = iy.saturating_sub(1);
                let yp = iy.min(grid.ny - 1);
                0.25 * (ey[grid.idx(ix, ym, iz)]
                    + ey[grid.idx(ix, yp, iz)]
                    + ey[grid.idx(ix, ym, iz + 1)]
                    + ey[grid.idx(ix, yp, iz + 1)])
            };
            [ex_v, ey_v, ez_v]
        }
    }
}

#[inline]
fn eval_nl_pol_and_diag_tan(
    dom: &NonlinearFdtdDomain,
    comp: usize,
    qpm_sign: f64,
    e_vec: [f64; 3],
    n_c: f64,
) -> (f64, f64, f64) {
    let e2 = e_vec[0] * e_vec[0] + e_vec[1] * e_vec[1] + e_vec[2] * e_vec[2];
    if !e2.is_finite() {
        return (0.0, 0.0, 0.0);
    }
    let mut p_comp = 0.0;
    let mut chi_aa = 0.0;

    if let Some(ref d_glob) = dom.chi2_d_global {
        let p2 = femcore::material::eval_chi2_polarization_time(d_glob, qpm_sign, e_vec);
        let j2 = femcore::material::eval_chi2_tangent_time(d_glob, qpm_sign, e_vec);
        p_comp += p2[comp];
        chi_aa += j2[comp][comp];
    }

    if let Some(ref r4) = dom.chi3_rank4_global {
        let mut s = 0.0;
        let mut t = 0.0;
        for j in 0..3 {
            for k in 0..3 {
                let ejk = e_vec[j] * e_vec[k];
                for l in 0..3 {
                    s += r4[comp][j][k][l] * ejk * e_vec[l];
                }
                t += (r4[comp][j][k][comp] + r4[comp][j][comp][k] + r4[comp][comp][j][k]) * ejk;
            }
        }
        p_comp += EPS0 * s;
        chi_aa += t;
    } else if let Some(ref c3) = dom.chi3 {
        if dom.chi3_iso_xxxx != 0.0 {
            let p3 = femcore::material::eval_chi3_kerr_time(c3, dom.n0, dom.euler_deg, e_vec);
            p_comp += p3[comp];
            chi_aa += dom.chi3_iso_xxxx * (e2 + 2.0 * e_vec[comp] * e_vec[comp]);
        }
    }

    if dom.k_fcd != 0.0 && n_c > 0.0 {
        let dfcd = 2.0 * dom.n0 * dom.k_fcd * n_c;
        p_comp += EPS0 * dfcd * e_vec[comp];
        chi_aa += dfcd;
    }

    let p_cap = 0.45 * EPS0 * dom.eps_min * e2.sqrt();
    if !p_comp.is_finite() {
        p_comp = 0.0;
    } else if p_comp.abs() > p_cap && p_cap > 0.0 {
        p_comp = p_comp.signum() * p_cap;
    }
    let tan_cap = 0.45 * dom.eps_min;
    if !chi_aa.is_finite() {
        chi_aa = 0.0;
    } else {
        chi_aa = chi_aa.clamp(-tan_cap, tan_cap);
    }

    let sig_nl = if dom.beta_tpa > 0.0 || (dom.sigma_fca > 0.0 && n_c > 0.0) {
        (2.0 / 3.0) * EPS0 * EPS0 * C0 * C0 * dom.n0 * dom.n0 * dom.beta_tpa * e2
            + EPS0 * C0 * dom.n0 * dom.sigma_fca * n_c
    } else {
        0.0
    };
    (p_comp, chi_aa, if sig_nl.is_finite() && sig_nl > 0.0 { sig_nl } else { 0.0 })
}

/// Solves the centered nonlinear Yee Ampere update via 3 synchronous Jacobi-Newton iterations (`0 ULP` across threads),
/// then updates `p_nl_prev` and advances TPA free-carrier densities via `femcore::material::drude_carrier_step`.
fn apply_nonlinear_e_update(
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    mats: &YeeMaterialCoeffs,
    coeffs: &mut AdeAndAmpereCoeffs,
    fields_re: &mut YeeFields,
    e_old_nl: &[Vec<f64>; 3],
    dt: f64,
) {
    let e_lin = save_nonlinear_e(&coeffs.nonlinear, fields_re);
    let mut e_next: [Vec<f64>; 3] = [
        vec![0.0; coeffs.nonlinear.edges[0].len()],
        vec![0.0; coeffs.nonlinear.edges[1].len()],
        vec![0.0; coeffs.nonlinear.edges[2].len()],
    ];
    let inv_dt = 1.0 / dt.max(1e-30);
    let run_par = femcore::par::threads() > 1;
    let cbs = [coeffs.cb_x.as_slice(), coeffs.cb_y.as_slice(), coeffs.cb_z.as_slice()];

    for _iter in 0..3 {
        for comp in 0..3 {
            let edges = &coeffs.nonlinear.edges[comp];
            if edges.is_empty() {
                continue;
            }
            let cb_slab = cbs[comp];
            let p_prev = &coeffs.nonlinear.p_nl_prev[comp];
            let nc_arr = &coeffs.nonlinear.n_carrier[comp];
            let elin_c = &e_lin[comp];
            let eold_c = &e_old_nl[comp];
            let doms = &coeffs.nonlinear.domains;

            femcore::par::for_chunks_mut(&mut e_next[comp], 256, run_par && edges.len() >= 512, |ci, chunk| {
                let k0 = ci * 256;
                for (loc, out_v) in chunk.iter_mut().enumerate() {
                    let k = k0 + loc;
                    let entry = &edges[k];
                    let dom = &doms[entry.dom_slot as usize];
                    let n_c = nc_arr.get(k).copied().unwrap_or(0.0);
                    let e_vec = sample_edge_e_vec(grid, fields_re, comp, entry);
                    let ea_m = e_vec[comp];
                    let (p_m, chi_aa, sig_nl) = eval_nl_pol_and_diag_tan(dom, comp, entry.qpm_sign, e_vec, n_c);
                    let cb_g = cb_slab[entry.idx];
                    let res = ea_m - elin_c[k]
                        + (cb_g * inv_dt) * (p_m - p_prev[k])
                        + cb_g * sig_nl * 0.5 * (ea_m + eold_c[k]);
                    let jac = (1.0 + (cb_g * inv_dt) * EPS0 * chi_aa + 0.5 * cb_g * sig_nl).max(0.25);
                    let upd = ea_m - res / jac;
                    *out_v = if upd.is_finite() { upd } else { elin_c[k] };
                }
            });
        }

        // Commit synchronous iterate to `fields_re` and resynchronize periodic/PEC faces
        for (entry, &v) in coeffs.nonlinear.edges[0].iter().zip(&e_next[0]) {
            fields_re.ex.as_mut_slice()[entry.idx] = v;
        }
        for (entry, &v) in coeffs.nonlinear.edges[1].iter().zip(&e_next[1]) {
            fields_re.ey.as_mut_slice()[entry.idx] = v;
        }
        for (entry, &v) in coeffs.nonlinear.edges[2].iter().zip(&e_next[2]) {
            fields_re.ez.as_mut_slice()[entry.idx] = v;
        }
        super::sync_periodic_e_high_faces(grid, bcs, fields_re, None);
        enforce_e_pec_clamp(grid, bcs, mats, fields_re, None);
    }

    // Store converged `P_NL^{n+1}` and advance TPA free-carrier density `N_c` via `drude_carrier_step`
    const TPA_RMS_SCALE: f64 = 0.816_496_580_927_726;
    let has_carriers = coeffs.nonlinear.has_carriers;
    let doms = &coeffs.nonlinear.domains;
    for comp in 0..3 {
        let edges = &coeffs.nonlinear.edges[comp];
        if edges.is_empty() {
            continue;
        }
        let p_prev = &mut coeffs.nonlinear.p_nl_prev[comp];
        let nc_arr = &mut coeffs.nonlinear.n_carrier[comp];
        let eold_c = &e_old_nl[comp];
        for (k, entry) in edges.iter().enumerate() {
            let dom = &doms[entry.dom_slot as usize];
            let n_c = nc_arr.get(k).copied().unwrap_or(0.0);
            let e_vec = sample_edge_e_vec(grid, fields_re, comp, entry);
            let (p_np1, _, _) = eval_nl_pol_and_diag_tan(dom, comp, entry.qpm_sign, e_vec, n_c);
            p_prev[k] = p_np1;
            if has_carriers && dom.beta_tpa > 0.0 && (dom.sigma_fca > 0.0 || dom.k_fcd != 0.0) {
                let mut e_mid = e_vec;
                e_mid[comp] = 0.5 * (e_vec[comp] + eold_c[k]);
                let e2_mid = e_mid[0] * e_mid[0] + e_mid[1] * e_mid[1] + e_mid[2] * e_mid[2];
                let i_eff = TPA_RMS_SCALE * C0 * EPS0 * dom.n0 * e2_mid;
                nc_arr[k] = femcore::material::drude_carrier_step(n_c, dt, dom.omega_0, dom.beta_tpa, dom.tau_c, i_eff);
            }
        }
    }
}

// ------------------------------------------------------------------ Source Injection (`H` & `E` TF/SF + Dipoles)

fn h_comp_mut(fields: &mut YeeFields, comp: usize) -> &mut [f64] {
    match comp {
        0 => fields.hx.as_mut_slice(),
        1 => fields.hy.as_mut_slice(),
        _ => fields.hz.as_mut_slice(),
    }
}

fn e_comp_mut(fields: &mut YeeFields, comp: usize) -> &mut [f64] {
    match comp {
        0 => fields.ex.as_mut_slice(),
        1 => fields.ey.as_mut_slice(),
        _ => fields.ez.as_mut_slice(),
    }
}

fn db_comp(coeffs: &AdeAndAmpereCoeffs, comp: usize) -> &[f64] {
    match comp {
        0 => &coeffs.db_x,
        1 => &coeffs.db_y,
        _ => &coeffs.db_z,
    }
}

fn cb_comp(coeffs: &AdeAndAmpereCoeffs, comp: usize) -> &[f64] {
    match comp {
        0 => &coeffs.cb_x,
        1 => &coeffs.cb_y,
        _ => &coeffs.cb_z,
    }
}

fn apply_source_h(
    src: &mut CompiledSource,
    grid: &YeeGrid,
    coeffs: &AdeAndAmpereCoeffs,
    fields: &mut YeeFields,
    t_n: f64,
    _dt: f64,
) {
    let stride_ax = |ax: usize| match ax {
        0 => grid.stride_x,
        1 => grid.stride_y,
        _ => 1,
    };
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });

    match src {
        CompiledSource::Mode { profile, waveform, inv_dp_h, .. }
        | CompiledSource::GaussianSheet { profile, waveform, inv_dp_h, .. } => {
            let (s_re, s_im) = waveform.eval_analytic(t_n);
            if s_re == 0.0 && s_im == 0.0 {
                return;
            }
            let st_a = stride_ax(profile.axis) as isize;
            let shift = (profile.h_cell_idx as isize - profile.plane_idx as isize) * st_a;
            let sg = profile.sign;
            let (u_ax, v_ax) = (profile.u_axis, profile.v_axis);

            // H_u(h_cell_idx) -= sg * db_u * E_v_inc(plane_idx) * inv_dp_h
            {
                let dbu = db_comp(coeffs, u_ax);
                let hu_arr = h_comp_mut(fields, u_ax);
                for &(g_ev, ev_c) in &profile.e_v_samples {
                    let g_hu = (g_ev as isize + shift) as usize;
                    let e_inc = ev_c.re * s_re - ev_c.im * s_im;
                    hu_arr[g_hu] -= sg * dbu[g_hu] * e_inc * *inv_dp_h;
                }
            }
            // H_v(h_cell_idx) += sg * db_v * E_u_inc(plane_idx) * inv_dp_h
            {
                let dbv = db_comp(coeffs, v_ax);
                let hv_arr = h_comp_mut(fields, v_ax);
                for &(g_eu, eu_c) in &profile.e_u_samples {
                    let g_hv = (g_eu as isize + shift) as usize;
                    let e_inc = eu_c.re * s_re - eu_c.im * s_im;
                    hv_arr[g_hv] += sg * dbv[g_hv] * e_inc * *inv_dp_h;
                }
            }
        }
        CompiledSource::PlaneSheet { axis, dir_sign, plane_idx, h_cell_idx, cos_pol, sin_pol, eta_bg, aux, .. } => {
            // Note: `aux.e` currently holds `E_1d^n`! Inject `E_1d^n` into 3D `H` BEFORE stepping `aux.h` to `n + 1/2`!
            let e_1d = aux.e_at_vertex(*plane_idx);
            if e_1d != 0.0 {
                let (u_ax, v_ax) = rh_uv_axes(*axis);
                let eu_inc = *cos_pol * e_1d;
                let ev_inc = *sin_pol * e_1d;
                let inv_dp = grid.axis_dp(*axis)[*h_cell_idx].recip();
                let sg = *dir_sign;
                let nu = grid.axis_centers(u_ax).len();
                let nv = grid.axis_centers(v_ax).len();
                let iu_hi = if is_2d && u_ax == 1 { 0 } else { nu };
                let iv_hi = if is_2d && v_ax == 1 { 0 } else { nv };

                let to_idx = |ia: usize, iu: usize, iv: usize| {
                    let mut ijk = [0usize; 3];
                    ijk[*axis] = ia;
                    ijk[u_ax] = iu;
                    ijk[v_ax] = iv;
                    grid.idx(ijk[0], ijk[1], ijk[2])
                };

                if ev_inc != 0.0 {
                    let dbu = db_comp(coeffs, u_ax);
                    let hu_arr = h_comp_mut(fields, u_ax);
                    for iu in 0..=iu_hi {
                        for iv in 0..nv {
                            let g = to_idx(*h_cell_idx, iu, iv);
                            hu_arr[g] -= sg * dbu[g] * ev_inc * inv_dp;
                        }
                    }
                }
                if eu_inc != 0.0 {
                    let dbv = db_comp(coeffs, v_ax);
                    let hv_arr = h_comp_mut(fields, v_ax);
                    for iu in 0..nu {
                        for iv in 0..=iv_hi {
                            let g = to_idx(*h_cell_idx, iu, iv);
                            hv_arr[g] += sg * dbv[g] * eu_inc * inv_dp;
                        }
                    }
                }
            }
            aux.step_h(t_n, *eta_bg);
        }
        CompiledSource::TfsfBox { i_lo, i_hi, axis, cos_pol, sin_pol, eta_bg, aux, .. } => {
            // `aux.e` holds `E_1d^n`! Apply 6-face (4-edge in 2D) TF/SF H-corrections, THEN advance `aux.step_h`!
            apply_tfsf_box_h(grid, coeffs, fields, *i_lo, *i_hi, *axis, *cos_pol, *sin_pol, aux);
            aux.step_h(t_n, *eta_bg);
        }
        CompiledSource::Dipole { component, magnetic: true, weights, waveform, .. } => {
            let s = waveform.eval(t_n);
            if s != 0.0 {
                let db = db_comp(coeffs, *component);
                let harr = h_comp_mut(fields, *component);
                for &(g, w) in weights.iter() {
                    harr[g] -= db[g] * s * w;
                }
            }
        }
        CompiledSource::Dipole { magnetic: false, .. } => {}
    }
}

fn apply_source_e(
    src: &mut CompiledSource,
    grid: &YeeGrid,
    coeffs: &AdeAndAmpereCoeffs,
    fields: &mut YeeFields,
    t_half: f64,
    dt: f64,
) {
    let stride_ax = |ax: usize| match ax {
        0 => grid.stride_x,
        1 => grid.stride_y,
        _ => 1,
    };
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });

    match src {
        CompiledSource::Mode { profile, waveform, h_time_offset_s, inv_dd_e, .. }
        | CompiledSource::GaussianSheet { profile, waveform, h_time_offset_s, inv_dd_e, .. } => {
            let (s_re, s_im) = waveform.eval_analytic(t_half + *h_time_offset_s);
            if s_re == 0.0 && s_im == 0.0 {
                return;
            }
            let st_a = stride_ax(profile.axis) as isize;
            let shift_to_e = (profile.plane_idx as isize - profile.h_cell_idx as isize) * st_a;
            let sg = profile.sign;
            let (u_ax, v_ax) = (profile.u_axis, profile.v_axis);

            // E_u(plane_idx) += sg * cb_u * H_v_inc(h_cell_idx) * inv_dd_e
            {
                let cbu = cb_comp(coeffs, u_ax);
                let eu_arr = e_comp_mut(fields, u_ax);
                for &(g_hv, hv_c) in &profile.h_v_samples {
                    let g_eu = (g_hv as isize + shift_to_e) as usize;
                    let h_inc = hv_c.re * s_re - hv_c.im * s_im;
                    eu_arr[g_eu] += sg * cbu[g_eu] * h_inc * *inv_dd_e;
                }
            }
            // E_v(plane_idx) -= sg * cb_v * H_u_inc(h_cell_idx) * inv_dd_e
            {
                let cbv = cb_comp(coeffs, v_ax);
                let ev_arr = e_comp_mut(fields, v_ax);
                for &(g_hu, hu_c) in &profile.h_u_samples {
                    let g_ev = (g_hu as isize + shift_to_e) as usize;
                    let h_inc = hu_c.re * s_re - hu_c.im * s_im;
                    ev_arr[g_ev] -= sg * cbv[g_ev] * h_inc * *inv_dd_e;
                }
            }
        }
        CompiledSource::PlaneSheet { axis, dir_sign, plane_idx, h_cell_idx, cos_pol, sin_pol, eta_bg, aux, .. } => {
            // `aux.h` was already stepped to `t_{n + 1/2}` at the end of `apply_source_h`!
            let h_1d = aux.h_at_cell(*h_cell_idx);
            if h_1d != 0.0 {
                let (u_ax, v_ax) = rh_uv_axes(*axis);
                let hu_inc = -*sin_pol * h_1d;
                let hv_inc = *cos_pol * h_1d;
                let inv_dd = grid.axis_dd(*axis)[*plane_idx].recip();
                let sg = *dir_sign;
                let nu = grid.axis_centers(u_ax).len();
                let nv = grid.axis_centers(v_ax).len();
                let iu_hi = if is_2d && u_ax == 1 { 0 } else { nu };
                let iv_hi = if is_2d && v_ax == 1 { 0 } else { nv };

                let to_idx = |ia: usize, iu: usize, iv: usize| {
                    let mut ijk = [0usize; 3];
                    ijk[*axis] = ia;
                    ijk[u_ax] = iu;
                    ijk[v_ax] = iv;
                    grid.idx(ijk[0], ijk[1], ijk[2])
                };

                if hv_inc != 0.0 {
                    let cbu = cb_comp(coeffs, u_ax);
                    let eu_arr = e_comp_mut(fields, u_ax);
                    for iu in 0..nu {
                        for iv in 0..=iv_hi {
                            let g = to_idx(*plane_idx, iu, iv);
                            eu_arr[g] += sg * cbu[g] * hv_inc * inv_dd;
                        }
                    }
                }
                if hu_inc != 0.0 {
                    let cbv = cb_comp(coeffs, v_ax);
                    let ev_arr = e_comp_mut(fields, v_ax);
                    for iu in 0..=iu_hi {
                        for iv in 0..nv {
                            let g = to_idx(*plane_idx, iu, iv);
                            ev_arr[g] -= sg * cbv[g] * hu_inc * inv_dd;
                        }
                    }
                }
            }
            aux.step_e(t_half, *eta_bg, dt);
        }
        CompiledSource::TfsfBox { i_lo, i_hi, axis, cos_pol, sin_pol, eta_bg, aux, .. } => {
            apply_tfsf_box_e(grid, coeffs, fields, *i_lo, *i_hi, *axis, *cos_pol, *sin_pol, aux);
            aux.step_e(t_half, *eta_bg, dt);
        }
        CompiledSource::Dipole { component, magnetic: false, weights, waveform, .. } => {
            let s = waveform.eval(t_half);
            if s != 0.0 {
                let cb = cb_comp(coeffs, *component);
                let earr = e_comp_mut(fields, *component);
                for &(g, w) in weights.iter() {
                    earr[g] -= cb[g] * s * w;
                }
            }
        }
        CompiledSource::Dipole { magnetic: true, .. } => {}
    }
}

/// Applies exact 6-face (4-face in 2D) TF/SF magnetic corrections using cyclic right-handed `(f_ax, fu, fv)`!
#[allow(clippy::too_many_arguments)]
fn apply_tfsf_box_h(
    grid: &YeeGrid,
    coeffs: &AdeAndAmpereCoeffs,
    fields: &mut YeeFields,
    i_lo: [usize; 3],
    i_hi: [usize; 3],
    prop_axis: usize,
    cos_pol: f64,
    sin_pol: f64,
    aux: &Aux1dYeeLine,
) {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let (pu, pv) = rh_uv_axes(prop_axis);

    // Evaluates incident `E_comp` at `(ijk_e, comp)` where `ijk_e` has vertex indices along axes `!= comp` and cell index along `comp`!
    // Notice `E_inc` only has nonzero components along `pu` (`cos_pol`) and `pv` (`sin_pol`), which BOTH sit at PRIMAL VERTICES along `prop_axis`!!
    // Wait: check that!! Because `pu != prop_axis` and `pv != prop_axis`, BOTH `E_{pu}` and `E_{pv}` sit at PRIMAL VERTEX `ijk[prop_axis]`!
    // There is ZERO spatial interpolation along `prop_axis` — `aux.e_at_vertex(ijk[prop_axis])` is EXACT!!
    let eval_e_inc = |comp: usize, ijk: [usize; 3]| -> f64 {
        let e1d = aux.e_at_vertex(ijk[prop_axis]);
        if comp == pu {
            cos_pol * e1d
        } else if comp == pv {
            sin_pol * e1d
        } else {
            0.0
        }
    };

    for f_ax in 0..3 {
        if is_2d && f_ax == 1 {
            continue;
        }
        let (fu, fv) = rh_uv_axes(f_ax);
        let fu_v_hi = if is_2d && fu == 1 { 0 } else { i_hi[fu] };
        let fv_v_hi = if is_2d && fv == 1 { 0 } else { i_hi[fv] };

        for &(is_hi, sg) in &[(false, 1.0_f64), (true, -1.0_f64)] {
            let ia_e = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] };
            let ia_h = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] - 1 };
            let inv_dp = grid.axis_dp(f_ax)[ia_h].recip();

            // 1. H_{fu}(ia_h, iu in i_lo[fu]..=fu_v_hi, iv in i_lo[fv]..i_hi[fv]) -= sg * db_{fu} * E_{fv}(ia_e, iu, iv) * inv_dp
            {
                let dbu = db_comp(coeffs, fu);
                let hu_arr = h_comp_mut(fields, fu);
                for iu in i_lo[fu]..=fu_v_hi {
                    for iv in i_lo[fv]..i_hi[fv] {
                        let mut ijk_e = [0usize; 3];
                        ijk_e[f_ax] = ia_e;
                        ijk_e[fu] = iu;
                        ijk_e[fv] = iv;
                        let ev_val = eval_e_inc(fv, ijk_e);
                        if ev_val != 0.0 {
                            let mut ijk_h = ijk_e;
                            ijk_h[f_ax] = ia_h;
                            let gh = grid.idx(ijk_h[0], ijk_h[1], ijk_h[2]);
                            hu_arr[gh] -= sg * dbu[gh] * ev_val * inv_dp;
                        }
                    }
                }
            }

            // 2. H_{fv}(ia_h, iu in i_lo[fu]..i_hi[fu], iv in i_lo[fv]..=fv_v_hi) += sg * db_{fv} * E_{fu}(ia_e, iu, iv) * inv_dp
            {
                let dbv = db_comp(coeffs, fv);
                let hv_arr = h_comp_mut(fields, fv);
                for iu in i_lo[fu]..i_hi[fu] {
                    for iv in i_lo[fv]..=fv_v_hi {
                        let mut ijk_e = [0usize; 3];
                        ijk_e[f_ax] = ia_e;
                        ijk_e[fu] = iu;
                        ijk_e[fv] = iv;
                        let eu_val = eval_e_inc(fu, ijk_e);
                        if eu_val != 0.0 {
                            let mut ijk_h = ijk_e;
                            ijk_h[f_ax] = ia_h;
                            let gh = grid.idx(ijk_h[0], ijk_h[1], ijk_h[2]);
                            hv_arr[gh] += sg * dbv[gh] * eu_val * inv_dp;
                        }
                    }
                }
            }
        }
    }
}

/// Applies exact 6-face (4-face in 2D) TF/SF electric corrections using cyclic right-handed `(f_ax, fu, fv)`!
#[allow(clippy::too_many_arguments)]
fn apply_tfsf_box_e(
    grid: &YeeGrid,
    coeffs: &AdeAndAmpereCoeffs,
    fields: &mut YeeFields,
    i_lo: [usize; 3],
    i_hi: [usize; 3],
    prop_axis: usize,
    cos_pol: f64,
    sin_pol: f64,
    aux: &Aux1dYeeLine,
) {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let (pu, pv) = rh_uv_axes(prop_axis);

    // Notice `H_inc` only has nonzero components along `pu` (`-sin_pol`) and `pv` (`+cos_pol`),
    // which BOTH sit at DUAL CELL CENTER `ijk[prop_axis]` (`c_{prop_axis}[ic]`)!!
    // Again, ZERO spatial interpolation along `prop_axis` — `aux.h_at_cell(ijk[prop_axis])` is ALGEBRAICALLY EXACT!!
    let eval_h_inc = |comp: usize, ijk: [usize; 3]| -> f64 {
        let h1d = aux.h_at_cell(ijk[prop_axis]);
        if comp == pu {
            -sin_pol * h1d
        } else if comp == pv {
            cos_pol * h1d
        } else {
            0.0
        }
    };

    for f_ax in 0..3 {
        if is_2d && f_ax == 1 {
            continue;
        }
        let (fu, fv) = rh_uv_axes(f_ax);
        let fu_v_hi = if is_2d && fu == 1 { 0 } else { i_hi[fu] };
        let fv_v_hi = if is_2d && fv == 1 { 0 } else { i_hi[fv] };

        for &(is_hi, sg) in &[(false, 1.0_f64), (true, -1.0_f64)] {
            let ia_e = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] };
            let ia_h = if is_hi { i_hi[f_ax] } else { i_lo[f_ax] - 1 };
            let inv_dd = grid.axis_dd(f_ax)[ia_e].recip();

            // 1. E_{fu}(ia_e, iu in i_lo[fu]..i_hi[fu], iv in i_lo[fv]..=fv_v_hi) += sg * cb_{fu} * H_{fv}(ia_h, iu, iv) * inv_dd
            {
                let cbu = cb_comp(coeffs, fu);
                let eu_arr = e_comp_mut(fields, fu);
                for iu in i_lo[fu]..i_hi[fu] {
                    for iv in i_lo[fv]..=fv_v_hi {
                        let mut ijk_h = [0usize; 3];
                        ijk_h[f_ax] = ia_h;
                        ijk_h[fu] = iu;
                        ijk_h[fv] = iv;
                        let hv_val = eval_h_inc(fv, ijk_h);
                        if hv_val != 0.0 {
                            let mut ijk_e = ijk_h;
                            ijk_e[f_ax] = ia_e;
                            let ge = grid.idx(ijk_e[0], ijk_e[1], ijk_e[2]);
                            eu_arr[ge] += sg * cbu[ge] * hv_val * inv_dd;
                        }
                    }
                }
            }

            // 2. E_{fv}(ia_e, iu in i_lo[fu]..=fu_v_hi, iv in i_lo[fv]..i_hi[fv]) -= sg * cb_{fv} * H_{fu}(ia_h, iu, iv) * inv_dd
            {
                let cbv = cb_comp(coeffs, fv);
                let ev_arr = e_comp_mut(fields, fv);
                for iu in i_lo[fu]..=fu_v_hi {
                    for iv in i_lo[fv]..i_hi[fv] {
                        let mut ijk_h = [0usize; 3];
                        ijk_h[f_ax] = ia_h;
                        ijk_h[fu] = iu;
                        ijk_h[fv] = iv;
                        let hu_val = eval_h_inc(fu, ijk_h);
                        if hu_val != 0.0 {
                            let mut ijk_e = ijk_h;
                            ijk_e[f_ax] = ia_e;
                            let ge = grid.idx(ijk_e[0], ijk_e[1], ijk_e[2]);
                            ev_arr[ge] -= sg * cbv[ge] * hu_val * inv_dd;
                        }
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------ Outer & Conformal PEC Clamping

fn enforce_e_pec_clamp(
    grid: &YeeGrid,
    bcs: &FdtdBoundaries,
    mats: &YeeMaterialCoeffs,
    fields_re: &mut YeeFields,
    fields_im: Option<&mut YeeFields>,
) {
    let is_2d = matches!(grid.dim, FdtdDim::TwoD { .. });
    let (nx, ny, nz, sx, sy) = (grid.nx, grid.ny, grid.nz, grid.stride_x, grid.stride_y);
    let iy_hi = if is_2d { 0 } else { ny };

    let clamp_one = |f: &mut YeeFields| {
        if matches!(bcs.faces[0], FdtdBcFace::Pec) {
            for iy in 0..=iy_hi {
                for iz in 0..=nz {
                    let g = iy * sy + iz;
                    f.ey.as_mut_slice()[g] = 0.0;
                    f.ez.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if matches!(bcs.faces[1], FdtdBcFace::Pec) {
            for iy in 0..=iy_hi {
                for iz in 0..=nz {
                    let g = nx * sx + iy * sy + iz;
                    f.ey.as_mut_slice()[g] = 0.0;
                    f.ez.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if !is_2d && matches!(bcs.faces[2], FdtdBcFace::Pec) {
            for ix in 0..=nx {
                for iz in 0..=nz {
                    let g = ix * sx + iz;
                    f.ex.as_mut_slice()[g] = 0.0;
                    f.ez.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if !is_2d && matches!(bcs.faces[3], FdtdBcFace::Pec) {
            for ix in 0..=nx {
                for iz in 0..=nz {
                    let g = ix * sx + ny * sy + iz;
                    f.ex.as_mut_slice()[g] = 0.0;
                    f.ez.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if matches!(bcs.faces[4], FdtdBcFace::Pec) {
            for ix in 0..=nx {
                for iy in 0..=iy_hi {
                    let g = ix * sx + iy * sy;
                    f.ex.as_mut_slice()[g] = 0.0;
                    f.ey.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if matches!(bcs.faces[5], FdtdBcFace::Pec) {
            for ix in 0..=nx {
                for iy in 0..=iy_hi {
                    let g = ix * sx + iy * sy + nz;
                    f.ex.as_mut_slice()[g] = 0.0;
                    f.ey.as_mut_slice()[g] = 0.0;
                }
            }
        }
        if mats.has_conformal_pec {
            for (v, &l) in f.ex.as_mut_slice().iter_mut().zip(mats.pec_len_x.as_slice()) {
                if l <= 0.0 {
                    *v = 0.0;
                }
            }
            for (v, &l) in f.ey.as_mut_slice().iter_mut().zip(mats.pec_len_y.as_slice()) {
                if l <= 0.0 {
                    *v = 0.0;
                }
            }
            for (v, &l) in f.ez.as_mut_slice().iter_mut().zip(mats.pec_len_z.as_slice()) {
                if l <= 0.0 {
                    *v = 0.0;
                }
            }
        }
    };
    clamp_one(fields_re);
    if let Some(fi) = fields_im {
        clamp_one(fi);
    }
}
