// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-42 Track B: SIMD Structured Subpixel FDTD (`Study::Fdtd`, ADR 145) verification suite:
//! 1. Broadband (`1.20 - 1.70 um`) dielectric slab Fresnel transmission/reflection vs exact Airy formula (`< 1e-3`),
//!    Lorentz ADE dispersive slab verification, and periodic grating `DiffractionMonitor` Rayleigh orders.
//! 2. FEM-seeded `ModeSource` / `ModeMonitor` (`nyx-ems::modes`) 3D waveguide width-step and `1 x 2` Y-splitter
//!    broadband S-parameters (`DrivenPortData`) vs full-wave 3D `ems::driven` (`< 1e-3` in `|S_ij|^2`).
//! 3. Empty 2D/3D `TFSF` box scattered-field leakage (`< 1e-10`) and 3D `PointDipole` DFT `CellField`s feeding
//!    `nyx-analysis::farfield::FarFieldAnalysis` Stratton–Chu NFFFT (`D_max = 1.50`).
//! 4. Conformal Farjadpour `Polarized` subpixel dielectric smoothing and Dey–Mittra `conformal_pec` curved-interface
//!    accuracy improvement vs `Staircased` FDTD on `geom` `FeatureList` curved cavities.
//! 5. SVD Matrix-Pencil / Prony `ResonanceFinder` on 2D and 3D lossy dielectric-loaded cavities vs `ems::cavity`
//!    (`< 0.5%` in both `f_0` and `Q`).
//! 6. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` 0-ULP bitwise-identical fields and monitor spectra (`slab_len > 16384`).

use femcore::geom::model::Model;
use femcore::linalg::csr::c;
use femcore::mesh::SimplexMesh;
use femcore::meshing::mesher::{extrude_y_with_remap, Layer};
use nyx_analysis::{Analysis, FarFieldAnalysis};
use nyx_ems::{
    run, run_fdtd, DispersiveDomainSpec, FdtdBcFace, FdtdBoundaries, FdtdMonitorData, FdtdMonitorSpec, FdtdSourceSpec,
    FdtdStudyConfig, FdtdWaveform, GridOverrideBox, Pol2d, Problem, Study, SubpixelMode, SusceptibilityPole,
    WavePortBc, C0,
};
use nyx_physics::{DrivenPortData, Mesh, Scalar};
use serde_json::json;
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;

const EPS0: f64 = 8.854_187_812_8e-12;

/// Structured 2D `(x, z)` `SimplexMesh` over explicit `x_grid` and `z_grid` coordinates.
fn structured_2d_mesh(x_grid: &[f64], z_grid: &[f64], mut domain_fn: impl FnMut(f64, f64) -> i32) -> SimplexMesh {
    let nx = x_grid.len() - 1;
    let nz = z_grid.len() - 1;
    let nv_x = nx + 1;
    let mut pts = Vec::with_capacity(2 * (nx + 1) * (nz + 1));
    for &z in z_grid {
        for &x in x_grid {
            pts.push(x);
            pts.push(z);
        }
    }
    let vid = |ix: usize, iz: usize| -> u32 { (iz * nv_x + ix) as u32 };
    let mut cells = Vec::with_capacity(6 * nx * nz);
    let mut cell_labels = Vec::with_capacity(2 * nx * nz);
    for iz in 0..nz {
        let zc = 0.5 * (z_grid[iz] + z_grid[iz + 1]);
        for ix in 0..nx {
            let xc = 0.5 * (x_grid[ix] + x_grid[ix + 1]);
            let dom = domain_fn(xc, zc);
            let (v00, v10, v11, v01) = (vid(ix, iz), vid(ix + 1, iz), vid(ix + 1, iz + 1), vid(ix, iz + 1));
            if xc < 0.0 {
                cells.extend_from_slice(&[v10, v11, v01, v10, v01, v00]);
            } else {
                cells.extend_from_slice(&[v00, v10, v11, v00, v11, v01]);
            }
            cell_labels.push(dom);
            cell_labels.push(dom);
        }
    }
    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    for iz in 0..nz {
        facets.extend_from_slice(&[vid(0, iz), vid(0, iz + 1)]);
        facet_labels.push(1);
        facets.extend_from_slice(&[vid(nx, iz), vid(nx, iz + 1)]);
        facet_labels.push(2);
    }
    for ix in 0..nx {
        facets.extend_from_slice(&[vid(ix, 0), vid(ix + 1, 0)]);
        facet_labels.push(3);
        facets.extend_from_slice(&[vid(ix, nz), vid(ix + 1, nz)]);
        facet_labels.push(4);
    }
    SimplexMesh::new(2, pts, cells, cell_labels, facets, facet_labels)
}

/// Conforming 6-tet-per-hex Kuhn 3D `SimplexMesh` over `[lo, hi]` with per-hex domain classifier.
fn kuhn_3d_mesh(lo: [f64; 3], hi: [f64; 3], n: [usize; 3], mut domain_fn: impl FnMut([f64; 3]) -> i32) -> SimplexMesh {
    let idx = |i: usize, j: usize, k: usize| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let mut pts = Vec::with_capacity(3 * (n[0] + 1) * (n[1] + 1) * (n[2] + 1));
    for k in 0..=n[2] {
        let z = lo[2] + (hi[2] - lo[2]) * (k as f64) / (n[2] as f64);
        for j in 0..=n[1] {
            let y = lo[1] + (hi[1] - lo[1]) * (j as f64) / (n[1] as f64);
            for i in 0..=n[0] {
                let x = lo[0] + (hi[0] - lo[0]) * (i as f64) / (n[0] as f64);
                pts.extend_from_slice(&[x, y, z]);
            }
        }
    }
    let perms = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
    let mut cells = Vec::with_capacity(24 * n[0] * n[1] * n[2]);
    let mut cell_labels = Vec::with_capacity(6 * n[0] * n[1] * n[2]);
    for k in 0..n[2] {
        let zc = lo[2] + (hi[2] - lo[2]) * (k as f64 + 0.5) / (n[2] as f64);
        for j in 0..n[1] {
            let yc = lo[1] + (hi[1] - lo[1]) * (j as f64 + 0.5) / (n[1] as f64);
            for i in 0..n[0] {
                let xc = lo[0] + (hi[0] - lo[0]) * (i as f64 + 0.5) / (n[0] as f64);
                let dom = domain_fn([xc, yc, zc]);
                for perm in perms {
                    let mut v = [0usize; 3];
                    let mut tet = [0u32; 4];
                    tet[0] = idx(i, j, k) as u32;
                    for (s, ax) in perm.into_iter().enumerate() {
                        v[ax] = 1;
                        tet[s + 1] = idx(i + v[0], j + v[1], k + v[2]) as u32;
                    }
                    let pt = |u: u32| [pts[3 * u as usize], pts[3 * u as usize + 1], pts[3 * u as usize + 2]];
                    let (p0, p1, p2, p3) = (pt(tet[0]), pt(tet[1]), pt(tet[2]), pt(tet[3]));
                    let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                    let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                    let e3 = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
                    let det = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
                        + e1[1] * (e2[2] * e3[0] - e2[0] * e3[2])
                        + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
                    if det < 0.0 {
                        tet.swap(2, 3);
                    }
                    cells.extend_from_slice(&tet);
                    cell_labels.push(dom);
                }
            }
        }
    }
    let mut face_counts: HashMap<[u32; 3], usize> = HashMap::new();
    for tet in cells.as_chunks::<4>().0 {
        for skip in 0..4 {
            let mut f = [0u32; 3];
            let mut m = 0;
            for (q, &vq) in tet.iter().enumerate() {
                if q != skip {
                    f[m] = vq;
                    m += 1;
                }
            }
            f.sort_unstable();
            *face_counts.entry(f).or_insert(0) += 1;
        }
    }
    let mut ext: Vec<[u32; 3]> = face_counts.into_iter().filter(|&(_, c)| c == 1).map(|(f, _)| f).collect();
    ext.sort_unstable();
    let mut facets = Vec::with_capacity(3 * ext.len());
    let mut flabels = Vec::with_capacity(ext.len());
    let tol = 1e-8 * (hi[0] - lo[0]).abs().max(1e-12);
    for f in ext {
        let c3 = [0, 1, 2]
            .map(|ax| (pts[3 * f[0] as usize + ax] + pts[3 * f[1] as usize + ax] + pts[3 * f[2] as usize + ax]) / 3.0);
        let lbl = if (c3[0] - lo[0]).abs() <= tol {
            1
        } else if (c3[0] - hi[0]).abs() <= tol {
            2
        } else if (c3[1] - lo[1]).abs() <= tol {
            3
        } else if (c3[1] - hi[1]).abs() <= tol {
            4
        } else if (c3[2] - lo[2]).abs() <= tol {
            5
        } else {
            6
        };
        facets.extend_from_slice(&f);
        flabels.push(lbl);
    }
    SimplexMesh::new(3, pts, cells, cell_labels, facets, flabels)
}

/// Exact normal-incidence Fabry–Perot Airy `(T, R)` of a slab of thickness `d_m` and complex refractive index `n_tilde = n - j k` in air.
fn airy_slab_tr(freq_hz: f64, d_m: f64, n_tilde: femcore::linalg::csr::C) -> (f64, f64) {
    let one = c(1.0, 0.0);
    let r01 = (one - n_tilde) / (one + n_tilde);
    let r10 = -r01;
    let t01_t10 = one - r01 * r01;
    let k0 = 2.0 * PI * freq_hz / C0;
    // exp(-j k0 n_tilde d): with n_tilde = n - j kappa (kappa >= 0), -j k0 (n - j kappa) d = -k0 kappa d - j k0 n d
    let arg = c(0.0, -k0 * d_m) * n_tilde;
    let exp_m1 = arg.exp();
    let exp_m2 = (arg * 2.0).exp();
    let denom = one + r01 * r10 * exp_m2;
    let t_amp = t01_t10 * exp_m1 / denom;
    let r_amp = (r01 + r10 * exp_m2) / denom;
    (t_amp.norm_sqr(), r_amp.norm_sqr())
}

#[test]
fn broadband_slab_fresnel_vs_airy_and_diffraction_orders() {
    // 11 wavelengths across 1.20 um .. 1.70 um
    let lambdas: Vec<f64> = (0..=10).map(|k| 1.20e-6 + 0.05e-6 * (k as f64)).collect();
    let freqs: Vec<f64> = lambdas.iter().map(|&lam| C0 / lam).collect();
    let f_lo = *freqs.iter().min_by(|a, b| a.total_cmp(b)).unwrap();
    let f_hi = *freqs.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    let wf = FdtdWaveform::broadband_for_range(f_lo, f_hi);

    // --- Part A: Lossless dielectric slab (`d = 0.50 um`, `eps_r = 4.0`, `n = 2.0`) vs exact Airy (`< 1e-3`) ---
    let d_slab = 0.50e-6;
    let z_lo_slab = -0.25e-6;
    let z_hi_slab = 0.25e-6;
    let x_per = [0.0, 0.025e-6, 0.05e-6];
    let z_dom = [-2.20e-6, -1.00e-6, z_lo_slab, z_hi_slab, 1.00e-6, 2.20e-6];
    let sm_slab = structured_2d_mesh(&x_per, &z_dom, |_x, z| if z > z_lo_slab && z < z_hi_slab { 2 } else { 1 });
    let mesh_slab = Mesh::from_solver_mesh(&sm_slab).unwrap();

    let mut cfg = FdtdStudyConfig {
        polarization_2d: Pol2d::TmEy,
        subpixel: SubpixelMode::Polarized,
        shutoff_energy_ratio: 1e-11,
        boundaries: FdtdBoundaries {
            faces: [
                FdtdBcFace::Periodic { bloch_k: 0.0 },
                FdtdBcFace::Periodic { bloch_k: 0.0 },
                FdtdBcFace::Pmc,
                FdtdBcFace::Pmc,
                FdtdBcFace::Cpml { layers: 22, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
                FdtdBcFace::Cpml { layers: 22, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
            ],
        },
        ..Default::default()
    };
    cfg.grid.auto.wavelength_min_m = 1.20e-6;
    cfg.grid.auto.min_steps_per_wave = 96.0;
    cfg.grid.snap_planes.push((2, z_lo_slab));
    cfg.grid.snap_planes.push((2, z_hi_slab));
    cfg.sources.push(FdtdSourceSpec::PlaneWave {
        name: "pw_inc".into(),
        axis: 2,
        position_m: -1.00e-6,
        direction: 1.0,
        pol_angle_rad: PI * 0.5, // E_v = E_y (`Pol2d::TmEy`: axis=2 -> u=0 (x), v=1 (y))
        eps_bg: 1.0,
        mu_bg: 1.0,
        waveform: wf.clone(),
    });
    cfg.monitors.push(FdtdMonitorSpec::DiffractionMonitor {
        name: "trans".into(),
        axis: 2,
        position_m: 1.00e-6,
        normal_sign: 1.0,
        eps_bg: 1.0,
        mu_bg: 1.0,
        max_order: 1,
        freqs_hz: freqs.clone(),
    });
    cfg.monitors.push(FdtdMonitorSpec::DiffractionMonitor {
        name: "refl".into(),
        axis: 2,
        position_m: -1.25e-6,
        normal_sign: -1.0,
        eps_bg: 1.0,
        mu_bg: 1.0,
        max_order: 1,
        freqs_hz: freqs.clone(),
    });

    let mut prob_slab = Problem::new(Study::Fdtd(Box::new(cfg.clone())), 1);
    prob_slab.materials.set(1, "eps_r", 1.0);
    prob_slab.materials.set(2, "eps_r", 4.0);
    let res_slab = run_fdtd(&mesh_slab, &prob_slab).expect("FDTD Airy slab solve");
    let mdata = FdtdMonitorData::of(&res_slab).expect("FdtdMonitorData present");
    let d_tx = &mdata.diffraction_monitors[0];
    let d_rx = &mdata.diffraction_monitors[1];

    let mut max_t_err = 0.0_f64;
    let mut max_r_err = 0.0_f64;
    for (fi, &f) in freqs.iter().enumerate() {
        let (t_exact, r_exact) = airy_slab_tr(f, d_slab, c(2.0, 0.0));
        let t_got = d_tx.total_efficiency[fi];
        let r_got = d_rx.total_efficiency[fi];
        max_t_err = max_t_err.max((t_got - t_exact).abs());
        max_r_err = max_r_err.max((r_got - r_exact).abs());
        assert!(
            (t_got + r_got - 1.0).abs() < 1e-3,
            "slab power conservation violated at f={f:.4e}: T={t_got:.6}, R={r_got:.6}"
        );
    }
    assert!(
        max_t_err < 1e-3 && max_r_err < 1e-3,
        "broadband Airy slab error exceeded 1e-3: max |dT| = {max_t_err:.3e}, max |dR| = {max_r_err:.3e}"
    );

    // --- Part B: Lorentz ADE dispersive slab vs exact complex-index Airy formula ---
    let mut cfg_ade = cfg.clone();
    let omega_0 = 2.0 * PI * (C0 / 0.85e-6);
    let gamma = 0.08 * omega_0;
    let eps_inf = 2.25;
    let delta_eps = 1.10;
    cfg_ade.dispersive_domains.insert(
        2,
        DispersiveDomainSpec {
            eps_inf,
            poles: vec![SusceptibilityPole { omega_0, gamma, strength: delta_eps * omega_0 * omega_0 }],
        },
    );
    let mut prob_ade = Problem::new(Study::Fdtd(Box::new(cfg_ade)), 1);
    prob_ade.materials.set(1, "eps_r", 1.0);
    prob_ade.materials.set(2, "eps_r", eps_inf);
    let res_ade = run_fdtd(&mesh_slab, &prob_ade).expect("FDTD Lorentz ADE slab solve");
    let mdata_ade = FdtdMonitorData::of(&res_ade).unwrap();
    for (fi, &f) in freqs.iter().enumerate() {
        let w = 2.0 * PI * f;
        // exp(+j w t): eps(w) = eps_inf + strength / (omega_0^2 - w^2 + j gamma w)
        let eps_w = c(eps_inf, 0.0) + c(delta_eps * omega_0 * omega_0, 0.0) / c(omega_0 * omega_0 - w * w, gamma * w);
        let mut n_c = eps_w.sqrt();
        if n_c.re < 0.0 {
            n_c = -n_c;
        }
        let (t_ex, r_ex) = airy_slab_tr(f, d_slab, n_c);
        let t_fdtd = mdata_ade.diffraction_monitors[0].total_efficiency[fi];
        let r_fdtd = mdata_ade.diffraction_monitors[1].total_efficiency[fi];
        assert!(
            (t_fdtd - t_ex).abs() < 2.5e-3 && (r_fdtd - r_ex).abs() < 2.5e-3,
            "Lorentz ADE slab mismatch at f={f:.4e}: FDTD (T={t_fdtd:.5}, R={r_fdtd:.5}) vs exact (T={t_ex:.5}, R={r_ex:.5})"
        );
    }

    // --- Part C: 2D Periodic Binary Dielectric Grating (`Lambda_x = 1.80 um > lambda`, Rayleigh `m = -1, 0, +1`) ---
    let pitch_x = 1.80e-6;
    let nx_g = 36usize;
    let x_grat: Vec<f64> = (0..=nx_g).map(|i| -0.5 * pitch_x + pitch_x * (i as f64) / (nx_g as f64)).collect();
    let z_grat = [-2.20e-6, -1.00e-6, -0.20e-6, 0.20e-6, 1.00e-6, 2.20e-6];
    let sm_grat =
        structured_2d_mesh(
            &x_grat,
            &z_grat,
            |x, z| {
                if z > -0.20e-6 && z < 0.20e-6 && x.abs() < 0.45e-6 {
                    2
                } else {
                    1
                }
            },
        );
    let mesh_grat = Mesh::from_solver_mesh(&sm_grat).unwrap();

    let mut cfg_grat = cfg;
    cfg_grat.grid.auto.min_steps_per_wave = 26.0;
    cfg_grat.grid.snap_planes = vec![(2, -0.20e-6), (2, 0.20e-6), (0, -0.45e-6), (0, 0.45e-6)];
    if let FdtdMonitorSpec::DiffractionMonitor { ref mut max_order, .. } = cfg_grat.monitors[0] {
        *max_order = 2;
    }
    if let FdtdMonitorSpec::DiffractionMonitor { ref mut max_order, .. } = cfg_grat.monitors[1] {
        *max_order = 2;
    }
    let mut prob_grat = Problem::new(Study::Fdtd(Box::new(cfg_grat)), 1);
    prob_grat.materials.set(1, "eps_r", 1.0);
    prob_grat.materials.set(2, "eps_r", 4.0);
    let res_grat = run_fdtd(&mesh_grat, &prob_grat).expect("FDTD periodic grating solve");
    let m_grat = FdtdMonitorData::of(&res_grat).unwrap();
    let tx_g = &m_grat.diffraction_monitors[0];
    let rx_g = &m_grat.diffraction_monitors[1];

    // Check midband `lambda = 1.45 um` (`fi = 5`): `m = +-1` orders must propagate symmetrically and total T + R == 1
    let fi_mid = 5;
    let eta_p1 =
        tx_g.orders_by_freq[fi_mid].iter().find(|o| o.m == 1 && o.n == 0).expect("m=+1 order").efficiency_total;
    let eta_m1 =
        tx_g.orders_by_freq[fi_mid].iter().find(|o| o.m == -1 && o.n == 0).expect("m=-1 order").efficiency_total;
    assert!(eta_p1 > 0.04, "expected non-trivial 1st-order grating diffraction efficiency, got {eta_p1:.5}");
    assert!((eta_p1 - eta_m1).abs() < 1e-5, "symmetric grating must give eta_+1 == eta_-1: {eta_p1} vs {eta_m1}");
    let tot_power = tx_g.total_efficiency[fi_mid] + rx_g.total_efficiency[fi_mid];
    assert!((tot_power - 1.0).abs() < 5e-3, "grating Rayleigh orders must conserve power: T + R = {tot_power:.5}");
}

fn tag_end_waveports(mesh: &mut SimplexMesh, y0: f64, y1: f64) {
    let tol_y = 1e-8 * (y1 - y0).abs().max(1e-12);
    let exterior = mesh.exterior_facets();
    for (f, &is_ext) in exterior.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let fv = mesh.facet(f);
        let ys = [mesh.point(fv[0] as usize)[1], mesh.point(fv[1] as usize)[1], mesh.point(fv[2] as usize)[1]];
        if ys.iter().all(|&y| (y - y0).abs() <= tol_y) {
            mesh.facet_labels[f] = 201;
        } else if ys.iter().all(|&y| (y - y1).abs() <= tol_y) {
            mesh.facet_labels[f] = 202;
        }
    }
}

fn tag_splitter_waveports(mesh: &mut SimplexMesh, y0: f64, y1: f64) {
    let tol_y = 1e-8 * (y1 - y0).abs().max(1e-12);
    let tol_x = 1e-12;
    let exterior = mesh.exterior_facets();
    for (f, &is_ext) in exterior.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let fv = mesh.facet(f);
        let p0 = mesh.point(fv[0] as usize);
        let p1 = mesh.point(fv[1] as usize);
        let p2 = mesh.point(fv[2] as usize);
        let xc = (p0[0] + p1[0] + p2[0]) / 3.0;
        let ys = [p0[1], p1[1], p2[1]];
        if ys.iter().all(|&y| (y - y0).abs() <= tol_y) {
            mesh.facet_labels[f] = 201;
        } else if ys.iter().all(|&y| (y - y1).abs() <= tol_y) {
            mesh.facet_labels[f] = if xc < 0.0 { 202 } else { 203 };
        }
    }
    // Tag interior `x = 0` symmetry-plane triangles with label `99` (`problem.pmc = vec![99]`) so `ports.rs`
    // places natural PMC (`n x H = 0`) on the `x = 0` boundary edge shared between `p2` and `p3`!
    let mut seen_mid = std::collections::BTreeSet::new();
    for ci in 0..mesh.n_cells() {
        let cv = mesh.cell(ci);
        for [a, b, d] in [[cv[0], cv[1], cv[2]], [cv[0], cv[1], cv[3]], [cv[0], cv[2], cv[3]], [cv[1], cv[2], cv[3]]] {
            if mesh.point(a as usize)[0].abs() <= tol_x
                && mesh.point(b as usize)[0].abs() <= tol_x
                && mesh.point(d as usize)[0].abs() <= tol_x
            {
                let mut key = [a, b, d];
                key.sort_unstable();
                if seen_mid.insert(key) {
                    mesh.facets.extend_from_slice(&[a, b, d]);
                    mesh.facet_labels.push(99);
                }
            }
        }
    }
}

#[test]
fn fem_seeded_mode_source_and_monitor_waveguide_step_and_splitter_vs_driven() {
    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let progress = femcore::progress::Progress::new();

    // --- Part A: 3D SOI-contrast strictly single-mode waveguide width-step (`w1 = 0.22 um -> w2 = 0.36 um`, `eps_core = 12.11`, `eps_clad = 2.085`) ---
    // Note: in lateral PEC box `W_box = 0.84 um` (`x in [-0.42, +0.42] um`), optical width `L_opt = 0.36 * 3.48 + 0.48 * 1.444 = 1.946 um < 1.5 * 1.55 um = 2.325 um`,
    // keeping ALL higher-order symmetric core & cladding modes (`TM_2`, `n_z >= 1`) strictly evanescent in both solvers!
    let x_grid = [
        -0.42e-6, -0.34e-6, -0.26e-6, -0.21e-6, -0.18e-6, -0.145e-6, -0.11e-6, -0.07e-6, -0.035e-6, 0.0, 0.035e-6,
        0.07e-6, 0.11e-6, 0.145e-6, 0.18e-6, 0.21e-6, 0.26e-6, 0.34e-6, 0.42e-6,
    ];
    let z_grid = [0.0, 0.10e-6, 0.20e-6];
    // Domain 1: SiO2 cladding (`eps = 2.085`), Domain 2: Si inner core `|x| < 0.11 um` (`eps = 12.11`),
    // Domain 3: step shoulder `0.11 < |x| < 0.18 um` (SiO2 in input half, remapped to Si in output half).
    let sec_step = structured_2d_mesh(&x_grid, &z_grid, |xc, _zc| {
        if xc.abs() < 0.11e-6 {
            2
        } else if xc.abs() < 0.18e-6 {
            3
        } else {
            1
        }
    });

    // 1. Solve reference 3D `ems::driven` on `y in [0, 1.60 um]` (step at `y = 0.80 um`)
    let l_half = 0.80e-6;
    let mut remap_wide = BTreeMap::new();
    remap_wide.insert(3, 2);
    let drv_layers = vec![
        (Layer { thickness: l_half, n: 8, domain_offset: 0, fill: None }, BTreeMap::new()),
        (Layer { thickness: l_half, n: 8, domain_offset: 0, fill: None }, remap_wide.clone()),
    ];
    let mut sm_drv = extrude_y_with_remap(&sec_step, &drv_layers, 0.0).unwrap();
    tag_end_waveports(&mut sm_drv, 0.0, 2.0 * l_half);
    let ex_drv = Mesh::from_solver_mesh(&sm_drv).unwrap();

    let mut p_drv = Problem::new(
        Study::Driven {
            f_min: freq,
            f_max: freq,
            n_points: 1,
            freqs_hz: Some(vec![freq]),
            prom_tol: None,
            max_hdm_solves: 1,
            strict_prom_tol: false,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p_drv.materials.set(1, "eps_r", 2.085);
    p_drv.materials.set(2, "eps_r", 12.11);
    p_drv.materials.set(3, "eps_r", 2.085);
    p_drv.waveports = vec![
        WavePortBc { labels: vec![201], name: "in".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![202], name: "out".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];
    let res_drv = run(&ex_drv, &p_drv, &progress).expect("3D driven SOI step solve");
    let dp_drv = DrivenPortData::of(&res_drv).unwrap();
    let s11_drv_sq = c(dp_drv.s_matrix[0][0][0][0], dp_drv.s_matrix[0][0][0][1]).norm_sqr();
    let s21_drv_sq = c(dp_drv.s_matrix[0][1][0][0], dp_drv.s_matrix[0][1][0][1]).norm_sqr();

    // 2. Solve 3D `Study::Fdtd` with 20-layer CPML lead extensions `y in [-1.35 um, +2.55 um]`, `ModeSource` at `y = -0.65 um`,
    //    and `ModeMonitor`s at `y = 0.0` (`in`) and `y = 1.60 um` (`out`)!
    let fdtd_layers = vec![
        (Layer { thickness: 2.15e-6, n: 4, domain_offset: 0, fill: None }, BTreeMap::new()),
        (Layer { thickness: 1.75e-6, n: 4, domain_offset: 0, fill: None }, remap_wide),
    ];
    let sm_fdtd = extrude_y_with_remap(&sec_step, &fdtd_layers, -1.35e-6).unwrap();
    let ex_fdtd = Mesh::from_solver_mesh(&sm_fdtd).unwrap();

    let wf_wg = FdtdWaveform::GaussianPulse {
        f0_hz: freq,
        fwidth_hz: 0.20 * freq,
        amplitude: 1.0,
        phase_rad: 0.0,
        peak_offset_sigma: 4.5,
    };
    let mut cfg_wg = FdtdStudyConfig {
        subpixel: SubpixelMode::Polarized,
        shutoff_energy_ratio: 1e-8,
        port_z0_ohm: 50.0,
        boundaries: FdtdBoundaries {
            faces: [
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
                FdtdBcFace::Cpml { layers: 20, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
                FdtdBcFace::Cpml { layers: 20, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
            ],
        },
        ..Default::default()
    };
    cfg_wg.grid.auto.wavelength_min_m = lam0;
    cfg_wg.grid.auto.min_steps_per_wave = 26.0;
    cfg_wg.grid.uniform_dl = Some([0.0, 0.0, 0.10e-6]);
    cfg_wg.grid.layer_refinements.push(nyx_ems::LayerRefinementSpec {
        axis: 0,
        bounds_m: [-0.20e-6, 0.20e-6],
        max_step_m: 0.012e-6,
        edge_step_m: None,
        edge_band_m: None,
    });
    cfg_wg.grid.snap_planes = vec![(0, -0.18e-6), (0, -0.11e-6), (0, 0.11e-6), (0, 0.18e-6), (1, 0.80e-6)];
    cfg_wg.sources.push(FdtdSourceSpec::ModeSource {
        name: "src_in".into(),
        axis: 1,
        position_m: -0.65e-6,
        direction: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: None,
        window_hi: None,
        waveform: wf_wg.clone(),
    });
    cfg_wg.monitors.push(FdtdMonitorSpec::ModeMonitor {
        name: "in".into(),
        axis: 1,
        position_m: 0.0,
        forward_sign: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: None,
        window_hi: None,
        freqs_hz: vec![freq],
    });
    cfg_wg.monitors.push(FdtdMonitorSpec::ModeMonitor {
        name: "out".into(),
        axis: 1,
        position_m: 2.0 * l_half,
        forward_sign: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: None,
        window_hi: None,
        freqs_hz: vec![freq],
    });

    let mut p_fdtd_wg = Problem::new(Study::Fdtd(Box::new(cfg_wg)), 2);
    p_fdtd_wg.materials.set(1, "eps_r", 2.085);
    p_fdtd_wg.materials.set(2, "eps_r", 12.11);
    p_fdtd_wg.materials.set(3, "eps_r", 2.085);
    let res_fdtd_wg = run_fdtd(&ex_fdtd, &p_fdtd_wg).expect("3D FDTD SOI step solve");
    let dp_fdtd = DrivenPortData::of(&res_fdtd_wg).expect("FDTD DrivenPortData emitted");
    let s11_fdtd_sq = c(dp_fdtd.s_matrix[0][0][0][0], dp_fdtd.s_matrix[0][0][0][1]).norm_sqr();
    let s21_fdtd_sq = c(dp_fdtd.s_matrix[0][1][0][0], dp_fdtd.s_matrix[0][1][0][1]).norm_sqr();

    assert!(
        (s21_fdtd_sq - s21_drv_sq).abs() < 1e-3 && (s11_fdtd_sq - s11_drv_sq).abs() < 1e-3,
        "FDTD vs Driven SOI waveguide step mismatch >= 1e-3: FDTD (|S21|^2={s21_fdtd_sq:.5}, |S11|^2={s11_fdtd_sq:.5}) vs Driven (|S21|^2={s21_drv_sq:.5}, |S11|^2={s11_drv_sq:.5})"
    );

    // --- Part B: 3D Symmetric `1 x 2` Y-Splitter (`p1` -> `p2`, `p3`) vs 3D `ems::driven` (`< 1e-3` in `|S_ij|^2`) ---
    // Lateral PEC box `x in [-0.42, +0.42] um` (`eps_clad = 2.25`, `eps_core = 6.25`: `(0.75 * 1.55 / 0.42)^2 = 7.66 > 6.25`,
    // so `TM_2` is deeply evanescent (`|gamma_y| > 6.8 um^-1`, `< 1e-6` power at `1.00 um` ports)):
    //   `|x| < 0.10 um` -> dom 2 (stem center: core in stem & junction, clad in output arms)
    //   `0.10 <= |x| < 0.25 um` -> dom 4 (clad in stem, core in junction AND output arms)
    //   `|x| >= 0.25 um` -> dom 1 (cladding everywhere)
    let x_split = [
        -0.42e-6, -0.335e-6, -0.25e-6, -0.20e-6, -0.15e-6, -0.10e-6, -0.05e-6, 0.0, 0.05e-6, 0.10e-6, 0.15e-6, 0.20e-6,
        0.25e-6, 0.335e-6, 0.42e-6,
    ];
    let z_split = [0.0, 0.10e-6, 0.20e-6];
    let sec_split = structured_2d_mesh(&x_split, &z_split, |xc, _zc| {
        let ax = xc.abs();
        if ax < 0.10e-6 {
            2
        } else if ax < 0.25e-6 {
            4
        } else {
            1
        }
    });
    let mut rm_stem = BTreeMap::new();
    rm_stem.insert(4, 1);
    let mut rm_junc = BTreeMap::new();
    rm_junc.insert(4, 2);
    let mut rm_arms = BTreeMap::new();
    rm_arms.insert(2, 1);
    rm_arms.insert(4, 2);

    let y_tot = 2.20e-6;
    let split_layers_drv = vec![
        (Layer { thickness: 0.80e-6, n: 8, domain_offset: 0, fill: None }, rm_stem.clone()),
        (Layer { thickness: 0.40e-6, n: 5, domain_offset: 0, fill: None }, rm_junc.clone()),
        (Layer { thickness: 1.00e-6, n: 10, domain_offset: 0, fill: None }, rm_arms.clone()),
    ];
    let mut sm_sp_drv = extrude_y_with_remap(&sec_split, &split_layers_drv, 0.0).unwrap();
    tag_splitter_waveports(&mut sm_sp_drv, 0.0, y_tot);
    let ex_sp_drv = Mesh::from_solver_mesh(&sm_sp_drv).unwrap();

    let mut p_sp_drv = Problem::new(
        Study::Driven {
            f_min: freq,
            f_max: freq,
            n_points: 1,
            freqs_hz: Some(vec![freq]),
            prom_tol: None,
            max_hdm_solves: 1,
            strict_prom_tol: false,
            seed_freqs_hz: vec![],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p_sp_drv.materials.set(1, "eps_r", 2.25);
    p_sp_drv.materials.set(2, "eps_r", 6.25);
    p_sp_drv.pmc = vec![99];
    p_sp_drv.waveports = vec![
        WavePortBc { labels: vec![201], name: "p1".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![202], name: "p2".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![203], name: "p3".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];
    let res_sp_drv = run(&ex_sp_drv, &p_sp_drv, &progress).expect("3D driven Y-splitter solve");
    let dp_sp_drv = DrivenPortData::of(&res_sp_drv).unwrap();
    let sp11_drv = c(dp_sp_drv.s_matrix[0][0][0][0], dp_sp_drv.s_matrix[0][0][0][1]).norm_sqr();
    let sp21_drv = c(dp_sp_drv.s_matrix[0][1][0][0], dp_sp_drv.s_matrix[0][1][0][1]).norm_sqr();
    let sp31_drv = c(dp_sp_drv.s_matrix[0][2][0][0], dp_sp_drv.s_matrix[0][2][0][1]).norm_sqr();

    let split_layers_fdtd = vec![
        (Layer { thickness: 2.15e-6, n: 4, domain_offset: 0, fill: None }, rm_stem),
        (Layer { thickness: 0.40e-6, n: 4, domain_offset: 0, fill: None }, rm_junc),
        (Layer { thickness: 1.95e-6, n: 4, domain_offset: 0, fill: None }, rm_arms),
    ];
    let sm_sp_fdtd = extrude_y_with_remap(&sec_split, &split_layers_fdtd, -1.35e-6).unwrap();
    let ex_sp_fdtd = Mesh::from_solver_mesh(&sm_sp_fdtd).unwrap();

    let mut cfg_sp = FdtdStudyConfig {
        subpixel: SubpixelMode::Polarized,
        shutoff_energy_ratio: 1e-8,
        port_z0_ohm: 50.0,
        boundaries: FdtdBoundaries {
            faces: [
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
                FdtdBcFace::Cpml { layers: 20, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
                FdtdBcFace::Cpml { layers: 20, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-8 },
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
            ],
        },
        ..Default::default()
    };
    cfg_sp.grid.auto.wavelength_min_m = lam0;
    cfg_sp.grid.auto.min_steps_per_wave = 36.0;
    cfg_sp.grid.uniform_dl = Some([0.0, 0.0, 0.10e-6]);
    cfg_sp.grid.layer_refinements.push(nyx_ems::LayerRefinementSpec {
        axis: 0,
        bounds_m: [-0.27e-6, 0.27e-6],
        max_step_m: 0.010e-6,
        edge_step_m: None,
        edge_band_m: None,
    });
    cfg_sp.grid.snap_planes =
        vec![(0, -0.25e-6), (0, -0.10e-6), (0, 0.0), (0, 0.10e-6), (0, 0.25e-6), (1, 0.80e-6), (1, 1.20e-6)];
    cfg_sp.sources.push(FdtdSourceSpec::ModeSource {
        name: "src_p1".into(),
        axis: 1,
        position_m: -0.65e-6,
        direction: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: None,
        window_hi: None,
        waveform: wf_wg,
    });
    cfg_sp.monitors.push(FdtdMonitorSpec::ModeMonitor {
        name: "p1".into(),
        axis: 1,
        position_m: 0.0,
        forward_sign: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: None,
        window_hi: None,
        freqs_hz: vec![freq],
    });
    // Note on `axis = 1` (`+y`): `rh_uv_axes(1) = (2, 0)` -> `u = z`, `v = x`!
    // So `window_lo = [z_min, x_min]`, `window_hi = [z_max, x_max]`!
    cfg_sp.monitors.push(FdtdMonitorSpec::ModeMonitor {
        name: "p2".into(),
        axis: 1,
        position_m: y_tot,
        forward_sign: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: Some([0.0, -0.42e-6]),
        window_hi: Some([0.20e-6, 0.0]),
        freqs_hz: vec![freq],
    });
    cfg_sp.monitors.push(FdtdMonitorSpec::ModeMonitor {
        name: "p3".into(),
        axis: 1,
        position_m: y_tot,
        forward_sign: 1.0,
        mode_index: 1,
        n_eff_guess: None,
        window_lo: Some([0.0, 0.0]),
        window_hi: Some([0.20e-6, 0.42e-6]),
        freqs_hz: vec![freq],
    });

    let mut p_sp_fdtd = Problem::new(Study::Fdtd(Box::new(cfg_sp)), 2);
    p_sp_fdtd.materials.set(1, "eps_r", 2.25);
    p_sp_fdtd.materials.set(2, "eps_r", 6.25);
    let res_sp_fdtd = run_fdtd(&ex_sp_fdtd, &p_sp_fdtd).expect("3D FDTD Y-splitter solve");
    let dp_sp_fdtd = DrivenPortData::of(&res_sp_fdtd).unwrap();
    let sp11_fdtd = c(dp_sp_fdtd.s_matrix[0][0][0][0], dp_sp_fdtd.s_matrix[0][0][0][1]).norm_sqr();
    let sp21_fdtd = c(dp_sp_fdtd.s_matrix[0][1][0][0], dp_sp_fdtd.s_matrix[0][1][0][1]).norm_sqr();
    let sp31_fdtd = c(dp_sp_fdtd.s_matrix[0][2][0][0], dp_sp_fdtd.s_matrix[0][2][0][1]).norm_sqr();

    assert!(
        (sp21_fdtd - sp31_fdtd).abs() < 1e-5,
        "symmetric Y-splitter arms must have equal transmission: {sp21_fdtd:.8} vs {sp31_fdtd:.8}"
    );
    assert!(
        (sp11_fdtd - sp11_drv).abs() < 1e-3
            && (sp21_fdtd - sp21_drv).abs() < 1e-3
            && (sp31_fdtd - sp31_drv).abs() < 1e-3,
        "FDTD vs Driven Y-splitter mismatch >= 1e-3: FDTD ({sp11_fdtd:.5}, {sp21_fdtd:.5}, {sp31_fdtd:.5}) vs Driven ({sp11_drv:.5}, {sp21_drv:.5}, {sp31_drv:.5})"
    );
}

#[test]
fn empty_tfsf_box_leakage_and_dipole_farfield_nffft() {
    // --- Part A: 3D Empty TFSF Box Scattered-Field Leakage (`< 1e-10`) ---
    // Domain `[-1.50, +1.50]^3 um` with `6`-layer CPML (`|x|,|y|,|z| > 0.90 um`), keeping TFSF box `[-0.30, +0.30]^3 um`
    // and all SF probes (`+-0.50 um`) strictly inside the unstretched vacuum interior!
    let sm3 = kuhn_3d_mesh([-1.50e-6; 3], [1.50e-6; 3], [6, 6, 6], |_| 1);
    let m3 = Mesh::from_solver_mesh(&sm3).unwrap();
    let f0 = C0 / 1.50e-6;
    let wf = FdtdWaveform::GaussianPulse {
        f0_hz: f0,
        fwidth_hz: 0.25 * f0,
        amplitude: 1.0,
        phase_rad: 0.0,
        peak_offset_sigma: 4.5,
    };

    let cpml6 = FdtdBcFace::Cpml { layers: 6, order: 3.5, kappa_max: 4.0, alpha_max: 0.05, target_reflection: 1e-7 };
    let mut cfg_tfsf = FdtdStudyConfig {
        max_steps: 260,
        shutoff_energy_ratio: 0.0,
        boundaries: FdtdBoundaries { faces: [cpml6; 6] },
        ..Default::default()
    };
    cfg_tfsf.grid.auto.wavelength_min_m = 1.50e-6;
    cfg_tfsf.grid.auto.min_steps_per_wave = 15.0;
    cfg_tfsf.sources.push(FdtdSourceSpec::Tfsf {
        name: "tfsf_vac".into(),
        lo: [-0.30e-6, -0.30e-6, -0.30e-6],
        hi: [0.30e-6, 0.30e-6, 0.30e-6],
        axis: 2,
        direction: 1.0,
        pol_angle_rad: 0.37, // oblique polarization exciting both Ex and Ey / Hx and Hy!
        eps_bg: 1.0,
        mu_bg: 1.0,
        waveform: wf,
    });
    // Monitor inside TF region (`[0, 0, 0]`) and 4 monitors in SF region outside all faces
    cfg_tfsf.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
        name: "tf_center".into(),
        position: [0.0, 0.0, 0.0],
        component: 0,
        sample_stride: 1,
        prony_freq_range_hz: None,
        prony_start_time_s: None,
    });
    for (nm, pos) in [
        ("sf_xp", [0.50e-6, 0.0, 0.0]),
        ("sf_yp", [0.0, 0.50e-6, 0.0]),
        ("sf_zp", [0.0, 0.0, 0.50e-6]),
        ("sf_zm", [0.0, 0.0, -0.50e-6]),
    ] {
        cfg_tfsf.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
            name: nm.into(),
            position: pos,
            component: 0,
            sample_stride: 1,
            prony_freq_range_hz: None,
            prony_start_time_s: None,
        });
    }

    let mut prob_tfsf = Problem::new(Study::Fdtd(Box::new(cfg_tfsf)), 1);
    prob_tfsf.materials.set(1, "eps_r", 1.0);
    let res_tfsf = run_fdtd(&m3, &prob_tfsf).expect("3D empty TFSF solve");
    let md_tfsf = FdtdMonitorData::of(&res_tfsf).unwrap();
    let peak_tf = md_tfsf.time_monitors[0].values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    assert!(peak_tf > 0.5, "TF interior must carry incident wave (got peak_tf = {peak_tf:.4e})");
    for smon in &md_tfsf.time_monitors[1..] {
        let peak_sf = smon.values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let leak_rel = peak_sf / peak_tf;
        assert!(
            leak_rel < 1e-10,
            "empty TFSF box scattered-field leakage at '{}' exceeded 1e-10: {leak_rel:.3e}",
            smon.name
        );
    }

    // --- Part B: 3D Hertzian `PointDipole` -> DFT `FieldMonitor` `CellField`s -> `FarFieldAnalysis` NFFFT (`D_max ~ 1.5`) ---
    // Vacuum box `[-0.90, +0.90]^3` m with `12 x 12 x 12` Kuhn tets so `b = 0.45` m (`3` hexes from center) aligns with interior tet faces!
    let sm_dip = kuhn_3d_mesh([-0.90; 3], [0.90; 3], [12, 12, 12], |_| 1);
    let m_dip = Mesh::from_solver_mesh(&sm_dip).unwrap();
    let f_dip = C0; // lambda = 1.0 m
    let mut cfg_dip = FdtdStudyConfig { shutoff_energy_ratio: 1e-8, ..Default::default() };
    cfg_dip.grid.auto.wavelength_min_m = 1.0;
    cfg_dip.grid.auto.min_steps_per_wave = 20.0;
    cfg_dip.sources.push(FdtdSourceSpec::PointDipole {
        name: "ez_dipole".into(),
        position: [0.0, 0.0, 0.0],
        component: 2, // z-oriented electric dipole
        magnetic: false,
        waveform: FdtdWaveform::GaussianPulse {
            f0_hz: f_dip,
            fwidth_hz: 0.25 * f_dip,
            amplitude: 1.0,
            phase_rad: 0.0,
            peak_offset_sigma: 4.5,
        },
    });
    cfg_dip.monitors.push(FdtdMonitorSpec::FieldMonitor {
        name: "nffft_vol".into(),
        lo: None,
        hi: None,
        freqs_hz: vec![f_dip],
    });
    let mut prob_dip = Problem::new(Study::Fdtd(Box::new(cfg_dip)), 2);
    prob_dip.materials.set(1, "eps_r", 1.0);
    let res_dip = run_fdtd(&m_dip, &prob_dip).expect("3D dipole FDTD solve");

    let nffft_params = json!({
        "surface": { "kind": "box", "lo": [-0.45, -0.45, -0.45], "hi": [0.45, 0.45, 0.45], "sampling": "faces" },
        "n_theta": 37,
        "n_phi": 24
    });
    let ff_out =
        FarFieldAnalysis.run_with_params(&m_dip, &res_dip, &nffft_params).expect("FarFieldAnalysis on FDTD CellFields");
    let d_max = ff_out["D_max"].as_f64().expect("D_max number");
    let theta_deg = ff_out["D_max_theta_deg"].as_f64().expect("D_max_theta_deg");
    assert!(
        (d_max - 1.5).abs() < 0.04,
        "FDTD Hertzian dipole Stratton-Chu NFFFT expected D_max ~ 1.50, got {d_max:.4}"
    );
    assert!(
        (theta_deg - 90.0).abs() < 1e-6,
        "z-dipole peak radiation must lie in the equatorial plane theta = 90 deg, got {theta_deg}"
    );
}

#[test]
fn conformal_subpixel_dielectric_and_pec_convergence() {
    // --- Part A: Dey-Mittra `conformal_pec` vs staircased PEC on a 2D `Pol2d::TeHy` circular PEC cylinder (`R = 0.40 m`) ---
    // Exact analytical TE_{11} resonance of circular PEC cylinder of radius R: `k R = j'_{1,1} = 1.841_183_781_340_659`.
    let r_pec = 0.40;
    let f_exact_pec = C0 * 1.841_183_781_340_659 / (2.0 * PI * r_pec);

    let mut mdl_pec = Model::new_2d();
    let c_pec = mdl_pec.circle("cav_circle", (0.0, 0.0), r_pec);
    mdl_pec.overlay(&[("cav", c_pec)]);
    let fl_pec = mdl_pec.list;
    let box_2d = structured_2d_mesh(&[-0.45, 0.0, 0.45], &[-0.45, 0.0, 0.45], |_, _| 1);
    let m_box_2d = Mesh::from_solver_mesh(&box_2d).unwrap();

    let run_circ_pec = |conformal_pec: bool| -> f64 {
        let mut cfg = FdtdStudyConfig {
            polarization_2d: Pol2d::TeHy,
            subpixel: SubpixelMode::Staircased,
            conformal_pec,
            outside_is_pec: true,
            max_steps: 2400,
            shutoff_energy_ratio: 0.0,
            boundaries: FdtdBoundaries { faces: [FdtdBcFace::Pec; 6] },
            geom_feature_list: Some(fl_pec.clone()),
            ..Default::default()
        };
        let dx = 0.042; // ~ 9.5 cells per radius (non-commensurate with R = 0.40!)
        let nx = 22usize;
        let vx: Vec<f64> = (0..=nx).map(|i| -0.5 * (nx as f64) * dx + (i as f64) * dx).collect();
        cfg.grid.explicit_vertices = Some([vx.clone(), vec![-0.5, 0.5], vx]);
        cfg.sources.push(FdtdSourceSpec::PointDipole {
            name: "hy_exc".into(),
            position: [0.14, 0.0, 0.07],
            component: 1,
            magnetic: true,
            waveform: FdtdWaveform::GaussianPulse {
                f0_hz: f_exact_pec,
                fwidth_hz: 0.20 * f_exact_pec,
                amplitude: 1.0,
                phase_rad: 0.0,
                peak_offset_sigma: 4.5,
            },
        });
        cfg.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
            name: "hy_probe".into(),
            position: [-0.13, 0.0, 0.09],
            component: 4, // Hy
            sample_stride: 1,
            prony_freq_range_hz: Some([0.85 * f_exact_pec, 1.15 * f_exact_pec]),
            prony_start_time_s: None,
        });
        let mut p = Problem::new(Study::Fdtd(Box::new(cfg)), 1);
        p.materials.set(1, "eps_r", 1.0);
        let res = run_fdtd(&m_box_2d, &p).expect("circular PEC cavity FDTD");
        let md = FdtdMonitorData::of(&res).unwrap();
        md.time_monitors[0].resonances[0].freq_hz
    };

    let f_conf = run_circ_pec(true);
    let f_stair_pec = run_circ_pec(false);
    let err_conf = (f_conf - f_exact_pec).abs() / f_exact_pec;
    let err_stair_pec = (f_stair_pec - f_exact_pec).abs() / f_exact_pec;
    assert!(
        err_conf < 5e-3 && err_conf * 2.0 < err_stair_pec,
        "Dey-Mittra conformal PEC must beat staircased PEC by >2x: err_conf={err_conf:.4e} (f={f_conf:.5e}) vs err_stair={err_stair_pec:.4e} (f={f_stair_pec:.5e}), exact={f_exact_pec:.5e}"
    );

    // --- Part B: Farjadpour `SubpixelMode::Polarized` vs `Staircased` across sub-cell-shifted dielectric rod radii (`eps_r = 6.25`) ---
    // Sweeping `r` across half a coarse cell (`N = 20`, `dx = 0.05 m`) eliminates accidental single-radius staircased zero-crossings!
    let run_rod = |radius: f64, n_cells: usize, mode: SubpixelMode| -> f64 {
        let mut mdl_rod = Model::new_2d();
        let r_air = mdl_rod.rect("air_box", (-0.5, -0.5), (0.5, 0.5));
        let r_rod = mdl_rod.circle("rod_cyl", (0.0, 0.0), radius);
        mdl_rod.overlay(&[("air", r_air), ("rod", r_rod)]);
        let fl_rod = mdl_rod.list;

        let vx: Vec<f64> = (0..=n_cells).map(|i| -0.5 + (i as f64) / (n_cells as f64)).collect();
        let f_guess = 1.30e8;
        let mut cfg = FdtdStudyConfig {
            polarization_2d: Pol2d::TeHy,
            subpixel: mode,
            conformal_pec: false,
            max_steps: (95 * n_cells).max(2000),
            shutoff_energy_ratio: 0.0,
            boundaries: FdtdBoundaries { faces: [FdtdBcFace::Pec; 6] },
            geom_feature_list: Some(fl_rod),
            ..Default::default()
        };
        cfg.grid.explicit_vertices = Some([vx.clone(), vec![-0.5, 0.5], vx]);
        cfg.sources.push(FdtdSourceSpec::PointDipole {
            name: "hy_src".into(),
            position: [0.09, 0.0, 0.04],
            component: 1,
            magnetic: true,
            waveform: FdtdWaveform::GaussianPulse {
                f0_hz: f_guess,
                fwidth_hz: 0.22 * f_guess,
                amplitude: 1.0,
                phase_rad: 0.0,
                peak_offset_sigma: 4.5,
            },
        });
        cfg.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
            name: "hy_prb".into(),
            position: [-0.08, 0.0, 0.06],
            component: 4,
            sample_stride: 1,
            prony_freq_range_hz: Some([0.90e8, 1.75e8]),
            prony_start_time_s: None,
        });
        let mut p = Problem::new(Study::Fdtd(Box::new(cfg)), 1);
        p.materials.set(1, "eps_r", 1.0);
        p.materials.set(2, "eps_r", 6.25);
        let res = run_fdtd(&m_box_2d, &p).expect("dielectric rod FDTD");
        FdtdMonitorData::of(&res).unwrap().time_monitors[0].resonances[0].freq_hz
    };

    let radii = [0.195, 0.205, 0.215];
    let mut sq_pol = 0.0;
    let mut sq_str = 0.0;
    for &r in &radii {
        let f_ref = run_rod(r, 84, SubpixelMode::Polarized);
        let f_pol = run_rod(r, 20, SubpixelMode::Polarized);
        let f_str = run_rod(r, 20, SubpixelMode::Staircased);
        sq_pol += ((f_pol - f_ref) / f_ref).powi(2);
        sq_str += ((f_str - f_ref) / f_ref).powi(2);
    }
    let rms_pol = (sq_pol / (radii.len() as f64)).sqrt();
    let rms_str = (sq_str / (radii.len() as f64)).sqrt();
    assert!(
        rms_pol < rms_str * 0.65 && rms_pol < 8e-3,
        "Polarized subpixel averaging must improve coarse-grid TeHy dielectric rod RMS error over Staircased: rms_pol={rms_pol:.4e} vs rms_str={rms_str:.4e}"
    );
}

#[test]
fn prony_resonance_finder_2d_and_3d_cavity_vs_ems_cavity() {
    let progress = femcore::progress::Progress::new();

    // --- Part A: 2D Lossy Dielectric-Loaded PEC Cavity (`Pol2d::TmEy`) vs Extruded 3D `ems::cavity` (`< 0.5%` in `f_0` and `Q`) ---
    // Cavity `[0, 0.04] x [0, 0.03]` m with off-center dielectric post `[0.012, 0.028] x [0.008, 0.022]` (`eps_r = 4.0`, `tan_delta = 0.008`)
    let xg: Vec<f64> = (0..=10).map(|i| 0.004 * (i as f64)).collect();
    let zg = [0.0, 0.004, 0.008, 0.0115, 0.015, 0.0185, 0.022, 0.026, 0.030];
    let sm2 = structured_2d_mesh(&xg, &zg, |xc, zc| {
        if (0.012..=0.028).contains(&xc) && (0.008..=0.022).contains(&zc) {
            2
        } else {
            1
        }
    });
    let m2 = Mesh::from_solver_mesh(&sm2).unwrap();

    // Extrude 1 thin cell (`dy = 0.004 m`) along `y` with PEC on all 6 walls (`E_x = E_z = 0` on `y = 0, dy`)
    // so 3D `ems::cavity` solves the exact 2D `TmEy` (`E = y_hat Ey(x, z)`, `k_y = 0`) mode!
    let layers_1 = vec![(Layer { thickness: 0.004, n: 1, domain_offset: 0, fill: None }, BTreeMap::new())];
    let sm2_ext = extrude_y_with_remap(&sm2, &layers_1, 0.0).unwrap();
    let m2_ext = Mesh::from_solver_mesh(&sm2_ext).unwrap();

    let mut p_cav2 = Problem::new(Study::Cavity { n_modes: 1, f_search: 3.8e9 }, 2);
    p_cav2.materials.set(1, "eps_r", 1.0);
    p_cav2.materials.set(2, "eps_r", 4.0);
    p_cav2.materials.set(2, "eps_r_im", 4.0 * 0.008);
    let r_cav2 = run(&m2_ext, &p_cav2, &progress).expect("2D-extruded cavity solve");
    let f_fem_2d = r_cav2.solutions[0].value.unwrap();
    let q_fem_2d = r_cav2
        .scalars
        .iter()
        .find_map(|(k, v)| match (k.as_str(), v) {
            ("Q1", Scalar::Num(q)) => Some(*q),
            _ => None,
        })
        .expect("Q1 scalar");

    let mut cfg_2d = FdtdStudyConfig {
        polarization_2d: Pol2d::TmEy,
        subpixel: SubpixelMode::Polarized,
        max_steps: 4500,
        shutoff_energy_ratio: 0.0,
        boundaries: FdtdBoundaries { faces: [FdtdBcFace::Pec; 6] },
        ..Default::default()
    };
    cfg_2d.grid.auto.wavelength_min_m = C0 / f_fem_2d;
    cfg_2d.grid.auto.min_steps_per_wave = 28.0;
    cfg_2d.grid.snap_planes = vec![(0, 0.012), (0, 0.028), (2, 0.008), (2, 0.022)];
    cfg_2d.sources.push(FdtdSourceSpec::PointDipole {
        name: "ey_src".into(),
        position: [0.018, 0.0, 0.014],
        component: 1,
        magnetic: false,
        waveform: FdtdWaveform::GaussianPulse {
            f0_hz: f_fem_2d,
            fwidth_hz: 0.12 * f_fem_2d,
            amplitude: 1.0,
            phase_rad: 0.0,
            peak_offset_sigma: 4.5,
        },
    });
    // Center `prony_freq_range_hz` symmetrically around `f_fem_2d` so `reference_freq_hz() == f_fem_2d` (`sigma = omega_0 eps_0 eps_im`)!
    cfg_2d.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
        name: "ey_ringdown".into(),
        position: [0.021, 0.0, 0.016],
        component: 1,
        sample_stride: 1,
        prony_freq_range_hz: Some([0.82 * f_fem_2d, 1.18 * f_fem_2d]),
        prony_start_time_s: None,
    });
    let mut p_fdtd_2d = Problem::new(Study::Fdtd(Box::new(cfg_2d)), 1);
    p_fdtd_2d.materials.set(1, "eps_r", 1.0);
    p_fdtd_2d.materials.set(2, "eps_r", 4.0);
    p_fdtd_2d.materials.set(2, "eps_r_im", 4.0 * 0.008);
    let r_fdtd_2d = run_fdtd(&m2, &p_fdtd_2d).expect("2D cavity FDTD Prony solve");
    let pr2 = &FdtdMonitorData::of(&r_fdtd_2d).unwrap().time_monitors[0].resonances[0];
    let df_2d = (pr2.freq_hz - f_fem_2d).abs() / f_fem_2d;
    let dq_2d = (pr2.q_factor - q_fem_2d).abs() / q_fem_2d;
    assert!(
        df_2d < 5e-3 && dq_2d < 5e-3,
        "2D Prony vs ems::cavity exceeded 0.5%: Prony (f0={:.6e}, Q={:.3}) vs FEM (f0={f_fem_2d:.6e}, Q={q_fem_2d:.3}) [df={df_2d:.4e}, dQ={dq_2d:.4e}]",
        pr2.freq_hz,
        pr2.q_factor
    );

    // --- Part B: 3D Inhomogeneously Dielectric-Loaded Lossy PEC Cavity vs 3D `ems::cavity` (`< 0.5%` in `f_0` and `Q`) ---
    // Box `[0, 0.04] x [0, 0.03] x [0, 0.02]` m with central dielectric block `[0.01, 0.03] x [0.0075, 0.0225] x [0, 0.02]` (`8 x 8 x 4` Kuhn tets)
    let sm3_cav = kuhn_3d_mesh([0.0, 0.0, 0.0], [0.04, 0.03, 0.02], [8, 8, 4], |c3| {
        if (0.01..=0.03).contains(&c3[0]) && (0.0075..=0.0225).contains(&c3[1]) {
            2
        } else {
            1
        }
    });
    let m3_cav = Mesh::from_solver_mesh(&sm3_cav).unwrap();

    let mut p_cav3 = Problem::new(Study::Cavity { n_modes: 1, f_search: 3.0e9 }, 2);
    p_cav3.materials.set(1, "eps_r", 1.0);
    p_cav3.materials.set(2, "eps_r", 4.0);
    p_cav3.materials.set(2, "eps_r_im", 4.0 * 0.010);
    let r_cav3 = run(&m3_cav, &p_cav3, &progress).expect("3D inhomogeneous cavity FEM solve");
    let f_fem_3d = r_cav3.solutions[0].value.unwrap();
    let q_fem_3d = r_cav3
        .scalars
        .iter()
        .find_map(|(k, v)| match (k.as_str(), v) {
            ("Q1", Scalar::Num(q)) => Some(*q),
            _ => None,
        })
        .expect("3D cavity Q1 scalar");

    let mut cfg_3d = FdtdStudyConfig {
        subpixel: SubpixelMode::Polarized,
        max_steps: 5200,
        shutoff_energy_ratio: 0.0,
        boundaries: FdtdBoundaries { faces: [FdtdBcFace::Pec; 6] },
        ..Default::default()
    };
    cfg_3d.grid.auto.wavelength_min_m = C0 / f_fem_3d;
    cfg_3d.grid.auto.min_steps_per_wave = 44.0;
    cfg_3d.grid.snap_planes = vec![(0, 0.01), (0, 0.03), (1, 0.0075), (1, 0.0225)];
    cfg_3d.sources.push(FdtdSourceSpec::PointDipole {
        name: "ez_cav_src".into(),
        position: [0.018, 0.014, 0.010],
        component: 2,
        magnetic: false,
        waveform: FdtdWaveform::GaussianPulse {
            f0_hz: f_fem_3d,
            fwidth_hz: 0.12 * f_fem_3d,
            amplitude: 1.0,
            phase_rad: 0.0,
            peak_offset_sigma: 4.5,
        },
    });
    cfg_3d.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
        name: "ez_ringdown".into(),
        position: [0.022, 0.016, 0.010],
        component: 2,
        sample_stride: 1,
        prony_freq_range_hz: Some([0.82 * f_fem_3d, 1.18 * f_fem_3d]),
        prony_start_time_s: None,
    });
    let mut p_fdtd_3d = Problem::new(Study::Fdtd(Box::new(cfg_3d)), 1);
    p_fdtd_3d.materials.set(1, "eps_r", 1.0);
    p_fdtd_3d.materials.set(2, "eps_r", 4.0);
    p_fdtd_3d.materials.set(2, "eps_r_im", 4.0 * 0.010);
    let r_fdtd_3d = run_fdtd(&m3_cav, &p_fdtd_3d).expect("3D cavity FDTD Prony solve");
    let pr3 = &FdtdMonitorData::of(&r_fdtd_3d).unwrap().time_monitors[0].resonances[0];
    let df_3d = (pr3.freq_hz - f_fem_3d).abs() / f_fem_3d;
    let dq_3d = (pr3.q_factor - q_fem_3d).abs() / q_fem_3d;
    assert!(
        df_3d < 5e-3 && dq_3d < 5e-3,
        "3D Prony vs ems::cavity exceeded 0.5%: Prony (f0={:.6e}, Q={:.3}) vs FEM (f0={f_fem_3d:.6e}, Q={q_fem_3d:.3}) [df={df_3d:.4e}, dQ={dq_3d:.4e}]",
        pr3.freq_hz,
        pr3.q_factor
    );
}

#[test]
fn fdtd_solve_is_bitwise_identical_across_thread_counts() {
    if let Ok(out_path) = std::env::var("NYX_FR42B_DET_OUT") {
        let vals = run_determinism_3d_workload();
        let mut bytes = Vec::with_capacity(8 * vals.len());
        for v in vals {
            bytes.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        std::fs::write(out_path, bytes).expect("write determinism payload");
        return;
    }

    let exe = std::env::current_exe().expect("current test binary");
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-fr42b-det-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args(["--exact", "fdtd_solve_is_bitwise_identical_across_thread_counts", "--nocapture"])
            .env("NYX_THREADS", nt)
            .env("NYX_FR42B_DET_OUT", &out_path)
            .output()
            .expect("spawn FR-42b determinism subprocess");
        assert!(out.status.success(), "child NYX_THREADS={nt} failed: {}", String::from_utf8_lossy(&out.stderr));
        let bytes = std::fs::read(&out_path).expect("read child bytes");
        let _ = std::fs::remove_file(&out_path);
        bytes
    };

    let b1 = run_child("1");
    let b8 = run_child("8");
    assert!(!b1.is_empty());
    assert_eq!(b1, b8, "3D FDTD fields/monitors differed between NYX_THREADS=1 and NYX_THREADS=8");
}

fn run_determinism_3d_workload() -> Vec<f64> {
    // `28 x 28 x 28` grid (`slab_len = 29^3 = 24_389 > 16_384` and `nx = 28 >= 8`) with non-CPML PEC/Periodic faces
    // so `femcore::par::for_chunks_mut` and `femcore::par::sum_range` multi-threaded paths execute actively!
    let sm = kuhn_3d_mesh([-0.5e-6; 3], [0.5e-6; 3], [3, 3, 3], |c3| {
        if c3[0].abs() < 0.2e-6 && c3[2].abs() < 0.2e-6 {
            2
        } else {
            1
        }
    });
    let mesh = Mesh::from_solver_mesh(&sm).unwrap();
    let v28: Vec<f64> = (0..=28).map(|i| -0.5e-6 + 1.0e-6 * (i as f64) / 28.0).collect();
    let f0 = C0 / 1.50e-6;

    let mut cfg = FdtdStudyConfig {
        max_steps: 90,
        shutoff_energy_ratio: 1e-12,
        subpixel: SubpixelMode::Polarized,
        boundaries: FdtdBoundaries {
            faces: [
                FdtdBcFace::Periodic { bloch_k: 1.2e6 },
                FdtdBcFace::Periodic { bloch_k: 1.2e6 },
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
                FdtdBcFace::Pec,
            ],
        },
        ..Default::default()
    };
    cfg.grid.explicit_vertices = Some([v28.clone(), v28.clone(), v28]);
    cfg.grid.override_boxes.push(GridOverrideBox { lo: [-0.2e-6; 3], hi: [0.2e-6; 3], dl: [None, None, None] });
    cfg.sources.push(FdtdSourceSpec::PointDipole {
        name: "det_dip".into(),
        position: [0.05e-6, -0.04e-6, 0.03e-6],
        component: 2,
        magnetic: false,
        waveform: FdtdWaveform::GaussianPulse {
            f0_hz: f0,
            fwidth_hz: 0.30 * f0,
            amplitude: 1.0,
            phase_rad: 0.2,
            peak_offset_sigma: 4.0,
        },
    });
    cfg.monitors.push(FdtdMonitorSpec::FluxMonitor {
        name: "det_flux".into(),
        axis: 1,
        position_m: 0.15e-6,
        normal_sign: 1.0,
        window_lo: None,
        window_hi: None,
        freqs_hz: vec![0.9 * f0, f0, 1.1 * f0],
    });
    cfg.monitors.push(FdtdMonitorSpec::FieldTimeMonitor {
        name: "det_time".into(),
        position: [-0.08e-6, 0.11e-6, -0.05e-6],
        component: 2,
        sample_stride: 1,
        prony_freq_range_hz: Some([0.5 * f0, 1.8 * f0]),
        prony_start_time_s: Some(0.0),
    });
    cfg.monitors.push(FdtdMonitorSpec::FieldMonitor { name: "det_vol".into(), lo: None, hi: None, freqs_hz: vec![f0] });

    let mut prob = Problem::new(Study::Fdtd(Box::new(cfg)), 1);
    prob.materials.set(1, "eps_r", 1.0);
    prob.materials.set(2, "eps_r", 4.5);
    prob.materials.set(2, "sigma", 2.0 * PI * f0 * EPS0 * 0.02);
    let res = run_fdtd(&mesh, &prob).expect("3D determinism FDTD solve");
    let md = FdtdMonitorData::of(&res).unwrap();

    let mut out: Vec<f64> = Vec::new();
    for s in &res.solutions {
        for f in &s.fields {
            out.extend_from_slice(&f.values);
        }
        for cf in &s.cell_fields {
            out.extend_from_slice(&cf.re);
            out.extend_from_slice(&cf.im);
        }
    }
    for fm in &md.flux_monitors {
        out.extend_from_slice(&fm.raw_flux_w_s2);
        out.extend_from_slice(&fm.flux_normalized);
    }
    for tm in &md.time_monitors {
        out.extend_from_slice(&tm.values);
        for r in &tm.resonances {
            out.push(r.freq_hz);
            out.push(r.q_factor);
        }
    }
    out
}
