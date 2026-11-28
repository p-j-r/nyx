// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! FR-37 FEM-native Eigenmode Expansion (`Study::Eme`, ADR 140) verification suite:
//! 1. Directional coupler crossover length `L_pi` & sinusoidal power-transfer curve vs 2D `ems::modes`
//!    supermode beat length (`< 0.2%`) and full-wave 3D `ems::driven` (`< 1%`), plus BLAKE3 slice deduplication.
//! 2. Dielectric waveguide step transition & continuous 3D adiabatic linear taper vs 3D `ems::driven` (`< 0.5%`),
//!    plus zero-re-solve `eme_mode_sweep` modal truncation convergence.
//! 3. `1,000`-period DBR Bragg grating via `O(log N_rep)` binary-exponentiation `eme_periodicity_sweep`
//!    vs exact analytical `2 x 2` Abeles transfer-matrix stopband & passband (`< 1e-5` in `< 1 s`).
//! 4. Cylindrical bent-waveguide metric (`conformal_bend_coord` & `bend_radius`), gyrotropic biorthogonal
//!    unitarity, and `nyx_analysis::network` Touchstone `.sNp` compatibility.

use femcore::linalg::csr::{c, C};
use femcore::mesh::SimplexMesh;
use femcore::meshing::mesher::{extrude_y_with_remap, Layer};
use nyx_ems::{
    conformal_bend_coord, eme_length_sweep, eme_mode_sweep, eme_periodicity_sweep, run, solve_waveguide, EmeData,
    EmeGridSpec, Problem, Study, WavePortBc, C0,
};
use nyx_physics::{DrivenPortData, Mesh, Scalar};
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::time::Instant;

/// Builds a structured 2D `(x, z)` cross-section mesh over explicit `x_grid` and `z_grid` breakpoints,
/// labeling each rectangular cell `(ix, iz)` via `domain_fn(xc, zc)` and marking exterior boundary edges:
/// `1 = xmin`, `2 = xmax`, `3 = zmin`, `4 = zmax`.
fn build_structured_2d_section(
    x_grid: &[f64],
    z_grid: &[f64],
    mut domain_fn: impl FnMut(f64, f64) -> i32,
) -> SimplexMesh {
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
            cells.extend_from_slice(&[v00, v10, v11, v00, v11, v01]);
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

/// Splits end-face facets (`y = y0` and `y = y1`) of a 3D `SimplexMesh` at `x = x_mid` into 4 waveport labels:
/// `101` (`y = y0, x <= x_mid`), `102` (`y = y0, x >= x_mid`), `103` (`y = y1, x <= x_mid`), `104` (`y = y1, x >= x_mid`).
fn tag_split_coupler_end_ports(mesh: &mut SimplexMesh, y0: f64, y1: f64, x_mid: f64) {
    let tol_y = 1e-8 * (y1 - y0).abs().max(1e-12);
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
        let on_y0 = [p0[1], p1[1], p2[1]].iter().all(|&y| (y - y0).abs() <= tol_y);
        let on_y1 = [p0[1], p1[1], p2[1]].iter().all(|&y| (y - y1).abs() <= tol_y);
        if on_y0 {
            mesh.facet_labels[f] = if xc < x_mid { 101 } else { 102 };
        } else if on_y1 {
            mesh.facet_labels[f] = if xc < x_mid { 103 } else { 104 };
        }
    }
}

/// Tags full end-face facets (`y = y0` -> `201`, `y = y1` -> `202`) on a 3D `SimplexMesh`.
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

#[test]
fn directional_coupler_beat_length_and_3d_driven_agreement() {
    // Dual-core dielectric-loaded PEC waveguide directional coupler at lambda0 = 1.55 um (f = C0 / 1.55e-6 Hz):
    // Height d = 0.24 um (< lambda0 / (2 n_core) = 0.31 um) with z-invariant cores selects pure TE_{m,0}
    // (`E = z_hat E_z(x)`, `phi == 0`) modes on standard PEC walls (`pmc` empty!), while 8 x-cells per half-face
    // (`W_half = 1.85 um`) give 15 P2 TE_{m,0} DOFs per split port so 4-mode Arnoldi converges in 2 cycles!
    // Cladding eps_clad = 2.89 (n = 1.70), two symmetric cores eps_core = 6.25 (n = 2.50):
    // Core 1 in x in [-0.66, -0.21] um, Core 2 in x in [+0.21, +0.66] um (gap = 0.42 um).
    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let x_grid = [
        -1.85e-6, -1.50e-6, -1.18e-6, -0.90e-6, -0.66e-6, -0.51e-6, -0.36e-6, -0.21e-6, 0.0, 0.21e-6, 0.36e-6, 0.51e-6,
        0.66e-6, 0.90e-6, 1.18e-6, 1.50e-6, 1.85e-6,
    ];
    let z_grid = [0.0, 0.08e-6, 0.16e-6, 0.24e-6];

    let sec2d = build_structured_2d_section(&x_grid, &z_grid, |xc, _zc| {
        let in_core1 = xc > -0.66e-6 && xc < -0.21e-6;
        let in_core2 = xc > 0.21e-6 && xc < 0.66e-6;
        if in_core1 || in_core2 {
            2
        } else {
            1
        }
    });

    // 1. Reference 2D supermode solve via `ems::modes`
    let mut p_modes = Problem::new(Study::Modes { freq, n_target: 2.35, evanescent: false, n_modes: 2 }, 2);
    p_modes.materials.set(1, "eps_r", 2.89);
    p_modes.materials.set(2, "eps_r", 6.25);
    let r_modes = solve_waveguide(&sec2d, &p_modes).expect("2D coupler supermode solve");
    // Fundamental even and odd supermodes are the top two effective indices
    let mut neffs: Vec<f64> = r_modes.solutions.iter().filter_map(|s| s.eigenvalue.map(|ev| ev[0])).collect();
    neffs.sort_by(|a, b| b.total_cmp(a));
    let (n_even, n_odd) = (neffs[0], neffs[1]);
    assert!(n_even > n_odd + 0.01, "expected split supermodes, got {neffs:?}");
    let l_pi_modes = lam0 / (2.0 * (n_even - n_odd));

    // 2. Build 3D coupler mesh of nominal length L_3db = 0.5 * L_pi with 4 split waveports (p1, p2 at y=0; p3, p4 at y=L)
    let l_3db = 0.5 * l_pi_modes;
    let ny_3db = ((l_3db / 0.09e-6).ceil() as usize).max(40);
    let layers = vec![(Layer { thickness: l_3db, n: ny_3db, domain_offset: 0, fill: None }, BTreeMap::new())];
    let mut mesh3_3db = extrude_y_with_remap(&sec2d, &layers, 0.0).expect("extrude 3D coupler");
    tag_split_coupler_end_ports(&mut mesh3_3db, 0.0, l_3db, 0.0);
    let ex3_3db = Mesh::from_solver_mesh(&mesh3_3db).expect("exchange mesh");

    let waveports = vec![
        WavePortBc { labels: vec![101], name: "p1".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![102], name: "p2".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![103], name: "p3".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![104], name: "p4".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];

    // Solve EME with 6 segments & 5 modes (`TE_{1..5, 0}` guided + box modes) to verify BLAKE3 deduplication (`unique_slices == 1`)
    // and compare directly against full-wave 3D `ems::driven`!
    let mut p_eme_multi = Problem::new(
        Study::Eme {
            freq,
            n_target: 1.95,
            evanescent: false,
            n_modes: 5,
            grid: EmeGridSpec::Uniform { n_cells: 6, n_modes: None, bend_radius: None },
            y_range: None,
            port_modes: 0,
            port_z0_ohm: 50.0,
            reconstruct_fields: true,
        },
        2,
    );
    p_eme_multi.materials.set(1, "eps_r", 2.89);
    p_eme_multi.materials.set(2, "eps_r", 6.25);
    p_eme_multi.waveports = waveports.clone();

    let progress = femcore::progress::Progress::new();
    let res_eme_multi = run(&ex3_3db, &p_eme_multi, &progress).expect("EME multi-segment coupler solve");
    let data_multi = EmeData::of(&res_eme_multi).expect("EmeData cached");
    assert_eq!(
        data_multi.unique_slices.len(),
        1,
        "BLAKE3 deduplication must collapse all 6 identical coupler slices to 1 unique slice"
    );

    // Verify 6-segment and 1-segment 5-mode EME give identical S-matrices at L = L_3db!
    let mut p_eme_5single = p_eme_multi.clone();
    p_eme_5single.study = Study::Eme {
        freq,
        n_target: 1.95,
        evanescent: false,
        n_modes: 5,
        grid: EmeGridSpec::Uniform { n_cells: 1, n_modes: None, bend_radius: None },
        y_range: None,
        port_modes: 0,
        port_z0_ohm: 50.0,
        reconstruct_fields: false,
    };
    let res_eme_5single = run(&ex3_3db, &p_eme_5single, &progress).expect("EME 1-segment 5-mode coupler solve");
    let dp_multi = DrivenPortData::of(&res_eme_multi).unwrap();
    let dp_5single = DrivenPortData::of(&res_eme_5single).unwrap();
    for r in 0..4 {
        for col in 0..4 {
            let diff = (dp_multi.s_matrix[0][r][col][0] - dp_5single.s_matrix[0][r][col][0])
                .hypot(dp_multi.s_matrix[0][r][col][1] - dp_5single.s_matrix[0][r][col][1]);
            assert!(diff < 1e-11, "6-segment vs 1-segment EME mismatch at ({r},{col}): {diff:.3e}");
        }
    }

    // 2-supermode 1-segment EME length sweep across `L in [0, 2 * L_pi]` with zero re-solves to extract pure supermode beat length `L_pi`
    let mut p_eme = p_eme_multi.clone();
    p_eme.study = Study::Eme {
        freq,
        n_target: 2.35,
        evanescent: false,
        n_modes: 2,
        grid: EmeGridSpec::Uniform { n_cells: 1, n_modes: None, bend_radius: None },
        y_range: None,
        port_modes: 0,
        port_z0_ohm: 50.0,
        reconstruct_fields: false,
    };
    let res_eme = run(&ex3_3db, &p_eme, &progress).expect("EME 1-segment 2-supermode coupler solve");
    let eme_data = EmeData::of(&res_eme).expect("EmeData");

    // Fine analytical length sweep around L_pi to extract EME crossover length L_pi_eme
    let n_pts = 401usize;
    let lengths: Vec<f64> = (0..n_pts).map(|i| (i as f64) * (2.0 * l_pi_modes) / ((n_pts - 1) as f64)).collect();
    let sweep = eme_length_sweep(&eme_data, 0, &lengths).expect("eme_length_sweep");
    assert!(sweep.elapsed_s < 0.5, "zero-re-solve length sweep must be sub-second, took {} s", sweep.elapsed_s);

    // Parabolic peak refinement on cross-port power |S_{3,0}(L)|^2 (`p4` is index 3, `p3` bar port is index 2)
    let mut best_k = 0usize;
    let mut best_p_cross = -1.0_f64;
    for (k, pt) in sweep.points.iter().enumerate() {
        assert!(pt.passivity_max_sigma <= 1.0 + 1e-8, "non-passive coupler S at k={k}: {}", pt.passivity_max_sigma);
        let p_cross = pt.power_from_port0[3];
        if k <= n_pts * 3 / 4 && p_cross > best_p_cross {
            best_p_cross = p_cross;
            best_k = k;
        }
    }
    let l_pi_eme = if best_k > 0 && best_k + 1 < n_pts {
        let y_m = sweep.points[best_k - 1].power_from_port0[3];
        let y_0 = sweep.points[best_k].power_from_port0[3];
        let y_p = sweep.points[best_k + 1].power_from_port0[3];
        let dl = lengths[1] - lengths[0];
        let denom = (y_m - 2.0 * y_0 + y_p).abs().max(1e-30);
        lengths[best_k] + 0.5 * dl * (y_m - y_p) / denom
    } else {
        lengths[best_k]
    };
    let rel_l_pi_err = (l_pi_eme - l_pi_modes).abs() / l_pi_modes;
    assert!(
        rel_l_pi_err < 2e-3,
        "EME directional coupler L_pi = {l_pi_eme:.6e} m vs 2D modes L_pi = {l_pi_modes:.6e} m (rel err = {rel_l_pi_err:.4e} >= 0.2%)"
    );
    assert!(best_p_cross > 0.92, "expected peak crossover power > 0.92 at L_pi, got {best_p_cross:.5}");

    // 3. Full-wave 3D `ems::driven` reference on the exact same 3D coupler mesh at L = L_3db!
    let mut p_driven = Problem::new(
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
    p_driven.materials.set(1, "eps_r", 2.89);
    p_driven.materials.set(2, "eps_r", 6.25);
    p_driven.waveports = waveports;

    let res_driven = run(&ex3_3db, &p_driven, &progress).expect("3D driven coupler solve");
    let dp_drv = DrivenPortData::of(&res_driven).expect("DrivenPortData");

    let p_bar_eme = c(dp_5single.s_matrix[0][2][0][0], dp_5single.s_matrix[0][2][0][1]).norm_sqr();
    let p_crs_eme = c(dp_5single.s_matrix[0][3][0][0], dp_5single.s_matrix[0][3][0][1]).norm_sqr();
    let p_bar_drv = c(dp_drv.s_matrix[0][2][0][0], dp_drv.s_matrix[0][2][0][1]).norm_sqr();
    let p_crs_drv = c(dp_drv.s_matrix[0][3][0][0], dp_drv.s_matrix[0][3][0][1]).norm_sqr();

    assert!(
        (p_bar_eme - p_bar_drv).abs() < 0.01 && (p_crs_eme - p_crs_drv).abs() < 0.01,
        "EME vs 3D driven coupler power mismatch >= 1%: EME (bar={p_bar_eme:.4}, cross={p_crs_eme:.4}) vs Driven (bar={p_bar_drv:.4}, cross={p_crs_drv:.4})"
    );
}

#[test]
fn waveguide_step_and_adiabatic_taper_match_3d_driven_within_half_percent() {
    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let progress = femcore::progress::Progress::new();

    // --- Part A: Symmetric waveguide width step transition (w1 = 0.54 um -> w2 = 0.68 um, TE_{m,0} PEC slab) ---
    // 16 x-cells (`W = 3.60 um`) give 31 P2 TE_{m,0} DOFs so warm-started 4-mode Arnoldi in `ports.rs` converges in 2 cycles!
    let x_grid = [
        -1.80e-6, -1.50e-6, -1.20e-6, -0.92e-6, -0.68e-6, -0.34e-6, -0.305e-6, -0.27e-6, 0.0, 0.27e-6, 0.305e-6,
        0.34e-6, 0.68e-6, 0.92e-6, 1.20e-6, 1.50e-6, 1.80e-6,
    ];
    let z_grid = [0.0, 0.08e-6, 0.16e-6, 0.24e-6];
    // Domain 1: cladding (eps = 2.56), Domain 2: inner core |x| < 0.27 um (eps = 6.25 in both halves),
    // Domain 3: shoulder 0.27 < |x| < 0.34 um (cladding in half 1, remapped to core 2 in half 2!).
    let sec_step = build_structured_2d_section(&x_grid, &z_grid, |xc, _zc| {
        if xc.abs() < 0.27e-6 {
            2
        } else if xc.abs() < 0.34e-6 {
            3
        } else {
            1
        }
    });
    let l_half = 0.60e-6;
    let mut remap_wide = BTreeMap::new();
    remap_wide.insert(3, 2);
    let step_layers = vec![
        (Layer { thickness: l_half, n: 5, domain_offset: 0, fill: None }, BTreeMap::new()),
        (Layer { thickness: l_half, n: 5, domain_offset: 0, fill: None }, remap_wide),
    ];
    let mut mesh3_step = extrude_y_with_remap(&sec_step, &step_layers, 0.0).expect("extrude step transition");
    tag_end_waveports(&mut mesh3_step, 0.0, 2.0 * l_half);
    let ex_step = Mesh::from_solver_mesh(&mesh3_step).unwrap();

    let step_wps = vec![
        WavePortBc { labels: vec![201], name: "in".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![202], name: "out".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];

    // Shift n_target = 2.12 keeps all 7 propagating TE_{1..7, 0} modes (4 symmetric + 3 antisymmetric) strictly closer to sigma than lambda = 0!
    let mut p_step_eme = Problem::new(
        Study::Eme {
            freq,
            n_target: 2.12,
            evanescent: false,
            n_modes: 7,
            grid: EmeGridSpec::Explicit {
                lengths: vec![l_half, l_half],
                n_modes: None,
                bend_radius: None,
                mesh_lengths: None,
            },
            y_range: None,
            port_modes: 1,
            port_z0_ohm: 50.0,
            reconstruct_fields: true,
        },
        2,
    );
    p_step_eme.materials.set(1, "eps_r", 2.56);
    p_step_eme.materials.set(2, "eps_r", 6.25);
    p_step_eme.materials.set(3, "eps_r", 2.56);
    p_step_eme.waveports = step_wps.clone();

    let res_step_eme = run(&ex_step, &p_step_eme, &progress).expect("EME step solve");
    let step_data = EmeData::of(&res_step_eme).unwrap();
    assert_eq!(step_data.unique_slices.len(), 2, "step transition has 2 unique cross-sections");
    assert_eq!(step_data.unique_slices[0].n_modes, 7, "all 7 TE_{{1..7,0}} modes must converge on slice 0");
    assert_eq!(step_data.unique_slices[1].n_modes, 7, "all 7 TE_{{1..7,0}} modes must converge on slice 1");
    // Exercise modal truncation sweep (`eme_mode_sweep`) across guided + box radiation modes
    let mode_sweep = eme_mode_sweep(&step_data, &[1, 3, 5, 7]).expect("eme_mode_sweep");
    assert_eq!(mode_sweep.points.len(), 4);
    let t_step_eme = mode_sweep.points.last().unwrap().power_from_port0[1];

    let mut p_step_drv = Problem::new(
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
    p_step_drv.materials = p_step_eme.materials.clone();
    p_step_drv.waveports = step_wps;
    let res_step_drv = run(&ex_step, &p_step_drv, &progress).expect("3D driven step solve");
    let dp_step_drv = DrivenPortData::of(&res_step_drv).unwrap();
    let t_step_drv = c(dp_step_drv.s_matrix[0][1][0][0], dp_step_drv.s_matrix[0][1][0][1]).norm_sqr();

    assert!(
        (t_step_eme - t_step_drv).abs() < 5e-3,
        "step waveguide transmission mismatch >= 0.5%: T_EME = {t_step_eme:.5} vs T_driven = {t_step_drv:.5}"
    );

    // --- Part B: Continuous 3D sloped adiabatic linear taper (exercising general marching-tets + non-matching TransferIndex overlap!) ---
    // Extrude a structured 2D mesh along y in [0, L_taper], then smoothly deform core vertices x(y) linearly from w1 = 0.46 um to w2 = 0.86 um!
    let x_base = [
        -1.80e-6, -1.45e-6, -1.12e-6, -0.82e-6, -0.55e-6, -0.23e-6, -0.115e-6, 0.0, 0.115e-6, 0.23e-6, 0.55e-6,
        0.82e-6, 1.12e-6, 1.45e-6, 1.80e-6,
    ];
    let sec_taper_base =
        build_structured_2d_section(&x_base, &z_grid, |xc, _zc| if xc.abs() < 0.23e-6 { 2 } else { 1 });
    let l_taper = 2.40e-6;
    let taper_layers = vec![(Layer { thickness: l_taper, n: 14, domain_offset: 0, fill: None }, BTreeMap::new())];
    let mut mesh3_taper = extrude_y_with_remap(&sec_taper_base, &taper_layers, 0.0).expect("extrude taper base");
    // Smoothly warp interior x coordinates (`|x| < 1.80e-6`) so the core widens linearly from half-width 0.23 um at y=0 to 0.43 um at y=L_taper!
    let w_outer = 1.80e-6;
    let hw0 = 0.23e-6;
    let hw1 = 0.43e-6;
    for v in 0..mesh3_taper.n_points() {
        let x0 = mesh3_taper.points[3 * v];
        let y = mesh3_taper.points[3 * v + 1];
        let s = (y / l_taper).clamp(0.0, 1.0);
        let hw_y = hw0 + (hw1 - hw0) * s;
        let ax = x0.abs();
        let x_new = if ax <= hw0 + 1e-12 {
            x0 * (hw_y / hw0)
        } else if ax < w_outer - 1e-12 {
            let t_clad = (ax - hw0) / (w_outer - hw0);
            x0.signum() * (hw_y + t_clad * (w_outer - hw_y))
        } else {
            x0
        };
        mesh3_taper.points[3 * v] = x_new;
    }
    mesh3_taper.check().expect("sloped 3D adiabatic taper mesh is valid");
    tag_end_waveports(&mut mesh3_taper, 0.0, l_taper);
    let ex_taper = Mesh::from_solver_mesh(&mesh3_taper).unwrap();

    let taper_wps = vec![
        WavePortBc { labels: vec![201], name: "in".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { labels: vec![202], name: "out".into(), mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];

    let mut p_taper_eme = Problem::new(
        Study::Eme {
            freq,
            n_target: 2.05,
            evanescent: false,
            n_modes: 6,
            grid: EmeGridSpec::Uniform { n_cells: 12, n_modes: None, bend_radius: None },
            y_range: None,
            port_modes: 1,
            port_z0_ohm: 50.0,
            reconstruct_fields: false,
        },
        2,
    );
    p_taper_eme.materials.set(1, "eps_r", 2.56);
    p_taper_eme.materials.set(2, "eps_r", 6.25);
    p_taper_eme.waveports = taper_wps.clone();

    let res_taper_eme = run(&ex_taper, &p_taper_eme, &progress).expect("EME continuous sloped taper solve");
    let taper_eme_data = EmeData::of(&res_taper_eme).unwrap();
    // Because the sidewalls slope continuously, marching-tets produces 12 distinct cross-sections!
    assert_eq!(taper_eme_data.unique_slices.len(), 12);
    let dp_taper_eme = DrivenPortData::of(&res_taper_eme).unwrap();
    let t_taper_eme = c(dp_taper_eme.s_matrix[0][1][0][0], dp_taper_eme.s_matrix[0][1][0][1]).norm_sqr();

    let mut p_taper_drv = Problem::new(
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
    p_taper_drv.materials = p_taper_eme.materials.clone();
    p_taper_drv.waveports = taper_wps;
    let res_taper_drv = run(&ex_taper, &p_taper_drv, &progress).expect("3D driven taper solve");
    let dp_taper_drv = DrivenPortData::of(&res_taper_drv).unwrap();
    let t_taper_drv = c(dp_taper_drv.s_matrix[0][1][0][0], dp_taper_drv.s_matrix[0][1][0][1]).norm_sqr();

    assert!(
        t_taper_drv > 0.95 && (t_taper_eme - t_taper_drv).abs() < 5e-3,
        "adiabatic taper EME vs 3D driven transmission mismatch >= 0.5%: T_EME = {t_taper_eme:.5} vs T_driven = {t_taper_drv:.5}"
    );
}

/// Exact analytical 1D Abeles / Airy reflectance `R = |S11|^2` and transmittance `T = |S21|^2`
/// of `(n_H, L_H; n_L, L_L)^{n_rep}` (with left port in `n_H` and right port in `n_L`), computed
/// with periodic renormalization so high-`n_rep` Bragg stopbands (`(n_H/n_L)^1000 ~ 10^166`) never overflow `f64`!
#[allow(clippy::too_many_arguments)]
fn abeles_dbr_rt(freq: f64, n_in: f64, n_out: f64, n_h: f64, l_h: f64, n_l: f64, l_l: f64, n_rep: usize) -> (f64, f64) {
    let k0 = 2.0 * PI * freq / C0;
    let layer_m = |n: f64, d: f64| -> [[C; 2]; 2] {
        let phi = k0 * n * d;
        let (sp, cp) = phi.sin_cos();
        [[c(cp, 0.0), c(0.0, sp / n)], [c(0.0, n * sp), c(cp, 0.0)]]
    };
    // Rescale M / ||M||_inf at each squaring step: R = |num_r / denom|^2 is projective (invariant to scalar multiples of M!),
    // and T = 1 - R for lossless real refractive indices (checked directly against unscaled det(M) == 1 when scale_log is moderate).
    let mul2_norm = |a: ([[C; 2]; 2], f64), b: ([[C; 2]; 2], f64)| -> ([[C; 2]; 2], f64) {
        let (ma, la) = a;
        let (mb, lb) = b;
        let mut prod = [
            [ma[0][0] * mb[0][0] + ma[0][1] * mb[1][0], ma[0][0] * mb[0][1] + ma[0][1] * mb[1][1]],
            [ma[1][0] * mb[0][0] + ma[1][1] * mb[1][0], ma[1][0] * mb[0][1] + ma[1][1] * mb[1][1]],
        ];
        let s = prod.iter().flatten().map(|z| z.norm()).fold(0.0_f64, f64::max).max(1e-300);
        for r in &mut prod {
            for v in r {
                *v /= s;
            }
        }
        (prod, la + lb + s.ln())
    };
    let mut base = mul2_norm((layer_m(n_h, l_h), 0.0), (layer_m(n_l, l_l), 0.0));
    let mut tot = ([[c(1.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), c(1.0, 0.0)]], 0.0);
    let mut e = n_rep;
    while e > 0 {
        if (e & 1) == 1 {
            tot = mul2_norm(tot, base);
        }
        e >>= 1;
        if e > 0 {
            base = mul2_norm(base, base);
        }
    }
    let (m, log_scale) = tot;
    let denom = m[0][0] * n_in + m[0][1] * (n_in * n_out) + m[1][0] + m[1][1] * n_out;
    let num_r = m[0][0] * n_in + m[0][1] * (n_in * n_out) - m[1][0] - m[1][1] * n_out;
    let r = (num_r / denom).norm_sqr();
    let t = if log_scale > 300.0 {
        0.0
    } else {
        let inv_scale = (-log_scale).exp();
        (n_out / n_in) * (c(2.0 * n_in * inv_scale, 0.0) / denom).norm_sqr()
    };
    (r, t)
}

#[test]
fn dbr_bragg_grating_1000_periods_binary_exponentiation_matches_abeles_tmm() {
    // Parallel-plate TEM waveguide (PEC at z = 0, H; PMC at x = 0, W) so the fundamental TEM mode is dispersionless
    // with exact n_eff == sqrt(eps_r)! One single unit cell [0, L_H + L_L] is meshed in 3D, and EME evaluates
    // N_rep = 1 .. 1,000 periods via O(log N_rep) Redheffer binary exponentiation!
    let lam_bragg = 1.55e-6;
    let (n_h, n_l) = (2.20_f64, 1.50_f64);
    let l_h = lam_bragg / (4.0 * n_h);
    let l_l = lam_bragg / (4.0 * n_l);

    let sec = build_structured_2d_section(&[0.0, 0.20e-6, 0.40e-6], &[0.0, 0.15e-6, 0.30e-6], |_, _| 1);
    let layers = vec![
        (Layer { thickness: l_h, n: 1, domain_offset: 0, fill: Some(1) }, BTreeMap::new()),
        (Layer { thickness: l_l, n: 1, domain_offset: 0, fill: Some(2) }, BTreeMap::new()),
    ];
    let mesh3_uc = extrude_y_with_remap(&sec, &layers, 0.0).expect("extrude 1-period DBR");
    let ex_uc = Mesh::from_solver_mesh(&mesh3_uc).unwrap();

    let t_start = Instant::now();
    let progress = femcore::progress::Progress::new();

    // Test across Bragg center, stopband edge, and passband frequencies!
    for &detune in &[1.00, 1.08, 1.18, 1.28, 1.38] {
        let freq = (C0 / lam_bragg) * detune;
        let mut p = Problem::new(
            Study::Eme {
                freq,
                n_target: 2.0,
                evanescent: false,
                n_modes: 1,
                grid: EmeGridSpec::PeriodicSubgrid {
                    unit_cell: Box::new(EmeGridSpec::Explicit {
                        lengths: vec![l_h, l_l],
                        n_modes: None,
                        bend_radius: None,
                        mesh_lengths: None,
                    }),
                    num_reps: 1000,
                    unit_length: None,
                },
                y_range: None,
                port_modes: 1,
                port_z0_ohm: 50.0,
                reconstruct_fields: false,
            },
            1,
        );
        p.materials.set(1, "eps_r", n_h * n_h);
        p.materials.set(2, "eps_r", n_l * n_l);
        // x-walls (labels 1 and 2) are PMC so the waveguide supports an exact plane-wave TEM mode!
        p.pmc = vec![1, 2];

        let res = run(&ex_uc, &p, &progress).expect("1000-period DBR EME solve");
        let data = EmeData::of(&res).expect("EmeData");
        assert_eq!(data.unique_slices.len(), 2, "1000-period DBR must solve only 2 unique slices");

        // Note: unit cell starts in n_H (left terminal) and ends in n_L (right terminal)!
        let reps_list = [1usize, 5, 25, 100, 500, 1000];
        let sweep = eme_periodicity_sweep(&data, 0, &reps_list).expect("eme_periodicity_sweep");
        for (&n_rep, pt) in reps_list.iter().zip(&sweep.points) {
            let r_eme = pt.power_from_port0[0];
            let t_eme = pt.power_from_port0[1];
            let (r_exact, t_exact) = abeles_dbr_rt(freq, n_h, n_l, n_h, l_h, n_l, l_l, n_rep);
            assert!(
                (r_eme - r_exact).abs() < 1e-5 && (t_eme - t_exact).abs() < 1e-5,
                "DBR mismatch at detune={detune}, N_rep={n_rep}: EME (R={r_eme:.8}, T={t_eme:.8}) vs Abeles (R={r_exact:.8}, T={t_exact:.8})"
            );
            assert!((r_eme + t_eme - 1.0).abs() < 1e-10, "unitarity drift at N_rep={n_rep}: R + T = {}", r_eme + t_eme);
        }
    }

    let elapsed = t_start.elapsed().as_secs_f64();
    assert!(
        elapsed < 1.0,
        "5-frequency x 1,000-period DBR binary-exponentiation EME sweep must finish in < 1 s, took {elapsed:.3} s"
    );
}

#[test]
fn bent_waveguide_and_gyrotropic_biorthogonal_eme_and_network_touchstone() {
    // 1. Validate conformal_bend_coord formula and singularity guard
    let u = conformal_bend_coord(0.5e-6, 10.0e-6).unwrap();
    assert!((u - 10.0e-6 * (1.05_f64).ln()).abs() < 1e-18);
    assert!(conformal_bend_coord(-10.0e-6, 10.0e-6).is_err(), "must reject x = -R_bend singularity");

    let lam0 = 1.55e-6;
    let freq = C0 / lam0;
    let x_grid = [-0.90e-6, -0.50e-6, -0.22e-6, 0.0, 0.22e-6, 0.50e-6, 0.90e-6];
    let z_grid = [-0.55e-6, -0.20e-6, 0.20e-6, 0.55e-6];
    let sec =
        build_structured_2d_section(
            &x_grid,
            &z_grid,
            |xc, zc| {
                if xc.abs() < 0.22e-6 && zc.abs() < 0.20e-6 {
                    2
                } else {
                    1
                }
            },
        );
    let layers = vec![(Layer { thickness: 2.0e-6, n: 2, domain_offset: 0, fill: None }, BTreeMap::new())];
    let mesh3 = extrude_y_with_remap(&sec, &layers, 0.0).unwrap();
    let ex3 = Mesh::from_solver_mesh(&mesh3).unwrap();
    let progress = femcore::progress::Progress::new();

    // Straight vs bent waveguide (`R_bend = 6 um`): outward centrifugal shift increases fundamental effective index!
    let solve_with_bend = |r_bend: Option<f64>| -> SolveResult {
        let mut p = Problem::new(
            Study::Eme {
                freq,
                n_target: 2.2,
                evanescent: false,
                n_modes: 2,
                grid: EmeGridSpec::Composite {
                    subgrids: vec![
                        EmeGridSpec::Explicit {
                            lengths: vec![0.5e-6],
                            n_modes: None,
                            bend_radius: Some(vec![None]),
                            mesh_lengths: None,
                        },
                        EmeGridSpec::Explicit {
                            lengths: vec![1.0e-6],
                            n_modes: None,
                            bend_radius: Some(vec![r_bend]),
                            mesh_lengths: None,
                        },
                        EmeGridSpec::Explicit {
                            lengths: vec![0.5e-6],
                            n_modes: None,
                            bend_radius: Some(vec![None]),
                            mesh_lengths: None,
                        },
                    ],
                },
                y_range: None,
                port_modes: 1,
                port_z0_ohm: 50.0,
                reconstruct_fields: false,
            },
            2,
        );
        p.materials.set(1, "eps_r", 2.10);
        p.materials.set(2, "eps_r", 6.25);
        run(&ex3, &p, &progress).expect("bent EME solve")
    };
    use nyx_physics::SolveResult;
    let res_straight = solve_with_bend(None);
    let res_bent = solve_with_bend(Some(6.0e-6));
    let d_straight = EmeData::of(&res_straight).unwrap();
    let d_bent = EmeData::of(&res_bent).unwrap();
    assert_eq!(d_straight.unique_slices.len(), 1);
    assert_eq!(d_bent.unique_slices.len(), 2, "straight + bent + straight reuses the straight slice via BLAKE3");
    assert!(
        d_bent.unique_slices[1].neff_fwd[0][0] > d_straight.unique_slices[0].neff_fwd[0][0] + 1e-4,
        "cylindrical bend metric must shift guided effective index upward: bent={} vs straight={}",
        d_bent.unique_slices[1].neff_fwd[0][0],
        d_straight.unique_slices[0].neff_fwd[0][0]
    );

    // 2. Lossless Faraday gyrotropic (`is_non_symmetric`) EME cascade preserves unitarity & runs through NetworkAnalysis Touchstone!
    let mut p_gyro = Problem::new(
        Study::Eme {
            freq,
            n_target: 2.2,
            evanescent: false,
            n_modes: 2,
            grid: EmeGridSpec::Uniform { n_cells: 2, n_modes: None, bend_radius: None },
            y_range: None,
            port_modes: 2,
            port_z0_ohm: 50.0,
            reconstruct_fields: true,
        },
        2,
    );
    p_gyro.materials.set(1, "eps_r", 2.10);
    p_gyro.materials.set(2, "eps_r", 6.25);
    // Add lossless Faraday gyrotropy g = [0, 0.15, 0] along propagation axis +y -> skew-Hermitian imaginary off-diagonal eps_xz = -i g_y, eps_zx = +i g_y
    p_gyro.materials.set_gyrotropy(2, [0.0, 0.15, 0.0]);
    let res_gyro = run(&ex3, &p_gyro, &progress).expect("gyrotropic EME solve");
    let sigma_max = res_gyro
        .scalars
        .iter()
        .find_map(|(k, v)| match (k.as_str(), v) {
            ("passivity_max_sigma", Scalar::Num(s)) => Some(*s),
            _ => None,
        })
        .unwrap();
    assert!(
        (sigma_max - 1.0).abs() < 1e-8,
        "lossless gyrotropic EME cascade must be unitary, got sigma_max = {sigma_max}"
    );

    // Verify `nyx_analysis::network` produces valid Touchstone v1.1 `.s4p` directly from `Study::Eme` result!
    let net_rep = nyx_analysis::network_with_params(&res_gyro, &serde_json::json!({ "touchstone_format": "RI" }))
        .expect("NetworkAnalysis on EME result");
    let snp = net_rep.touchstone_snp.expect("Touchstone text");
    assert!(snp.contains("# Hz S RI R 50"), "unexpected Touchstone header: {snp}");
}
