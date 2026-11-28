// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration and verification tests for FR-42 Track A: High-Order Curved/Mixed-Cell Local Time-Stepping
//! Discontinuous Galerkin Time-Domain solver (`Study::Dgtd`, `"dgtd"`, ADR 144).

mod common;

use femcore::linalg::csr::c;
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use nyx_ems::{
    assemble_dgtd, dgtd_fields, run_dgtd, solve_dgtd, DgtdConfig, DgtdDispersionPole, DgtdFlux, DgtdFluxMonitorSpec,
    DgtdPolarization, DgtdPulse, DgtdTimeStepper, Excitation, ExcitationSpec, Problem, Sheet2dBc, Study, WavePortBc,
    C0, DGTD_DATA_KEY,
};
use nyx_physics::{DrivenPortData, Mesh};

const ETA0: f64 = 376.730_313_668;

/// Builds a 2D structured triangle mesh on `[0, lx] x [0, ly]` where the first `nx_fine` columns have width `dx_fine`
/// (domain 1) and the remaining `nx_coarse` columns have width `dx_coarse` (domain 2).
/// Boundary labels: `1` (`x = 0`), `2` (`x = L_x`), `3` (`y = 0`), `4` (`y = L_y`), `5` (interior interface `x = nx_fine * dx_fine`).
fn build_multiscale_2d_mesh(
    nx_fine: usize,
    dx_fine: f64,
    nx_coarse: usize,
    dx_coarse: f64,
    ny: usize,
    ly: f64,
) -> Mesh {
    let nx = nx_fine + nx_coarse;
    let mut xs = Vec::with_capacity(nx + 1);
    xs.push(0.0);
    for i in 0..nx_fine {
        xs.push((i + 1) as f64 * dx_fine);
    }
    let x_split = nx_fine as f64 * dx_fine;
    for i in 0..nx_coarse {
        xs.push(x_split + (i + 1) as f64 * dx_coarse);
    }
    let dy = ly / (ny as f64);

    let mut points = Vec::with_capacity((nx + 1) * (ny + 1) * 2);
    for j in 0..=ny {
        let y = (j as f64) * dy;
        for &x in &xs {
            points.push(x);
            points.push(y);
        }
    }
    let vid = |i: usize, j: usize| -> u32 { (j * (nx + 1) + i) as u32 };

    let mut cells = Vec::with_capacity(nx * ny * 6);
    let mut cell_labels = Vec::with_capacity(nx * ny * 2);
    for j in 0..ny {
        for i in 0..nx {
            let dom = if i < nx_fine { 1 } else { 2 };
            let (v00, v10, v01, v11) = (vid(i, j), vid(i + 1, j), vid(i, j + 1), vid(i + 1, j + 1));
            cells.extend_from_slice(&[v00, v10, v11]);
            cell_labels.push(dom);
            cells.extend_from_slice(&[v00, v11, v01]);
            cell_labels.push(dom);
        }
    }

    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    // x = 0 (label 1) and x = Lx (label 2)
    for j in 0..ny {
        facets.extend_from_slice(&[vid(0, j), vid(0, j + 1)]);
        facet_labels.push(1);
        facets.extend_from_slice(&[vid(nx, j), vid(nx, j + 1)]);
        facet_labels.push(2);
        if nx_fine > 0 && nx_fine < nx {
            facets.extend_from_slice(&[vid(nx_fine, j), vid(nx_fine, j + 1)]);
            facet_labels.push(5);
        }
    }
    // y = 0 (label 3) and y = Ly (label 4)
    for i in 0..nx {
        facets.extend_from_slice(&[vid(i, 0), vid(i + 1, 0)]);
        facet_labels.push(3);
        facets.extend_from_slice(&[vid(i, ny), vid(i + 1, ny)]);
        facet_labels.push(4);
    }

    let sm = SimplexMesh::new(2, points, cells, cell_labels, facets, facet_labels);
    Mesh::from_solver_mesh(&sm).expect("valid 2D multiscale mesh")
}

/// Builds a 3D tetrahedral block mesh on `[0, lx] x [0, ly] x [0, lz]` with an interior facet plane at `y = (ny / 2) * dy` (label 7)
/// and three y-domains: domain 1 (`0 <= j < j_slab0`), domain 2 (`j_slab0 <= j < j_slab1`), domain 1 (`j_slab1 <= j < ny`).
/// Boundary labels match `common::block`: `1` (`x = 0`), `4` (`x = lx`), `2` (`y = 0`), `5` (`y = ly`), `3` (`z = 0`), `6` (`z = lz`),
/// plus interior sheet plane `7` at `j = ny / 2`.
fn build_3d_waveguide_mesh(
    lx: f64,
    ly: f64,
    lz: f64,
    nx: usize,
    ny: usize,
    nz: usize,
    slab_j_range: (usize, usize),
) -> Mesh {
    let (dx, dy, dz) = (lx / nx as f64, ly / ny as f64, lz / nz as f64);
    let mut points = Vec::with_capacity((nx + 1) * (ny + 1) * (nz + 1) * 3);
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                points.extend_from_slice(&[i as f64 * dx, j as f64 * dy, k as f64 * dz]);
            }
        }
    }
    let vid = |i: usize, j: usize, k: usize| -> u32 { ((k * (ny + 1) + j) * (nx + 1) + i) as u32 };
    let mut cells = Vec::with_capacity(nx * ny * nz * 24);
    let mut cell_labels = Vec::with_capacity(nx * ny * nz * 6);

    for k in 0..nz {
        for j in 0..ny {
            let dom = if j >= slab_j_range.0 && j < slab_j_range.1 { 2 } else { 1 };
            for i in 0..nx {
                let v000 = vid(i, j, k);
                let v100 = vid(i + 1, j, k);
                let v010 = vid(i, j + 1, k);
                let v110 = vid(i + 1, j + 1, k);
                let v001 = vid(i, j, k + 1);
                let v101 = vid(i + 1, j, k + 1);
                let v011 = vid(i, j + 1, k + 1);
                let v111 = vid(i + 1, j + 1, k + 1);
                let tets = [
                    [v000, v100, v110, v111],
                    [v000, v110, v010, v111],
                    [v000, v010, v011, v111],
                    [v000, v011, v001, v111],
                    [v000, v001, v101, v111],
                    [v000, v101, v100, v111],
                ];
                for t in tets {
                    cells.extend_from_slice(&t);
                    cell_labels.push(dom);
                }
            }
        }
    }

    let mut facets = Vec::new();
    let mut facet_labels = Vec::new();
    // x = 0 (label 1) and x = lx (label 4)
    for k in 0..nz {
        for j in 0..ny {
            facets.extend_from_slice(&[vid(0, j, k), vid(0, j + 1, k), vid(0, j + 1, k + 1)]);
            facet_labels.push(1);
            facets.extend_from_slice(&[vid(0, j, k), vid(0, j + 1, k + 1), vid(0, j, k + 1)]);
            facet_labels.push(1);

            facets.extend_from_slice(&[vid(nx, j, k), vid(nx, j + 1, k + 1), vid(nx, j + 1, k)]);
            facet_labels.push(4);
            facets.extend_from_slice(&[vid(nx, j, k), vid(nx, j, k + 1), vid(nx, j + 1, k + 1)]);
            facet_labels.push(4);
        }
    }
    // y = 0 (label 2), y = ly (label 5), and interior plane j = ny / 2 (label 7)
    let j_mid = ny / 2;
    for k in 0..nz {
        for i in 0..nx {
            facets.extend_from_slice(&[vid(i, 0, k), vid(i + 1, 0, k + 1), vid(i + 1, 0, k)]);
            facet_labels.push(2);
            facets.extend_from_slice(&[vid(i, 0, k), vid(i, 0, k + 1), vid(i + 1, 0, k + 1)]);
            facet_labels.push(2);

            facets.extend_from_slice(&[vid(i, ny, k), vid(i + 1, ny, k), vid(i + 1, ny, k + 1)]);
            facet_labels.push(5);
            facets.extend_from_slice(&[vid(i, ny, k), vid(i + 1, ny, k + 1), vid(i, ny, k + 1)]);
            facet_labels.push(5);

            if j_mid > 0 && j_mid < ny {
                facets.extend_from_slice(&[vid(i, j_mid, k), vid(i + 1, j_mid, k), vid(i + 1, j_mid, k + 1)]);
                facet_labels.push(7);
                facets.extend_from_slice(&[vid(i, j_mid, k), vid(i + 1, j_mid, k + 1), vid(i, j_mid, k + 1)]);
                facet_labels.push(7);
            }
        }
    }
    // z = 0 (label 3) and z = lz (label 6)
    for j in 0..ny {
        for i in 0..nx {
            facets.extend_from_slice(&[vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0)]);
            facet_labels.push(3);
            facets.extend_from_slice(&[vid(i, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0)]);
            facet_labels.push(3);

            facets.extend_from_slice(&[vid(i, j, nz), vid(i + 1, j + 1, nz), vid(i + 1, j, nz)]);
            facet_labels.push(6);
            facets.extend_from_slice(&[vid(i, j, nz), vid(i, j + 1, nz), vid(i + 1, j + 1, nz)]);
            facet_labels.push(6);
        }
    }

    let sm = SimplexMesh::new(3, points, cells, cell_labels, facets, facet_labels);
    Mesh::from_solver_mesh(&sm).expect("valid 3D waveguide mesh")
}

#[test]
fn broadband_s_parameters_match_frequency_domain_driven_across_telecom_band() {
    // 3D TEM waveguide (`1.2 - 1.7 um`, 176.35 THz to 249.83 THz) with PEC on z = 0, lz (labels 3, 6),
    // PMC on x = 0, lx (labels 1, 4), and waveports on y = 0 (label 2) and y = ly (label 5).
    let (lx, ly, lz) = (0.4e-6, 1.5e-6, 0.4e-6);
    let mesh = build_3d_waveguide_mesh(lx, ly, lz, 1, 15, 1, (5, 10));

    let lambdas_m = [1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.55e-6, 1.6e-6, 1.7e-6];
    let freqs_hz: Vec<f64> = lambdas_m.iter().rev().map(|&lam| C0 / lam).collect();
    let f0 = 0.5 * (freqs_hz[0] + freqs_hz[freqs_hz.len() - 1]);
    let fwidth = 0.35 * f0;

    // 1. Frequency-domain reference solve via Study::Driven
    let mut p_driven = Problem::new(
        Study::Driven {
            f_min: freqs_hz[0],
            f_max: *freqs_hz.last().unwrap(),
            n_points: freqs_hz.len(),
            freqs_hz: Some(freqs_hz.clone()),
            prom_tol: None,
            max_hdm_solves: 12,
            strict_prom_tol: false,
            seed_freqs_hz: Vec::new(),
            solver: "direct".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p_driven.materials.set(1, "eps_r", 1.0);
    p_driven.materials.set(2, "eps_r", 1.0);
    p_driven.pec = vec![3, 6];
    p_driven.pmc = vec![1, 4];
    p_driven.waveports = vec![
        WavePortBc { name: "p1".into(), labels: vec![2], mode: 1, n_eff_guess: None, z0_ohm: None },
        WavePortBc { name: "p2".into(), labels: vec![5], mode: 1, n_eff_guess: None, z0_ohm: None },
    ];
    let res_driven =
        nyx_ems::run(&mesh, &p_driven, &femcore::progress::Progress::new()).expect("driven reference solve");
    let pd_driven = DrivenPortData::of(&res_driven).expect("driven port data");

    // 2. Time-domain solve via Study::Dgtd (3-stage kernel: assemble_dgtd -> solve_dgtd -> dgtd_fields)
    let dgtd_cfg = DgtdConfig {
        t_final: 50e-15,
        cfl: 0.75,
        flux: DgtdFlux::Upwind,
        time_stepper: DgtdTimeStepper::Lts,
        max_lts_tiers: 4,
        polarization: DgtdPolarization::Full,
        basis_order: Some(2),
        pulse: DgtdPulse::GaussianModulated {
            f0_hz: f0,
            fwidth_hz: fwidth,
            t0_s: None,
            amplitude: 1.0,
            phase_rad: 0.0,
        },
        freqs_hz: freqs_hz.clone(),
        excited_port: 1,
        absorbing_labels: Vec::new(),
        probes: vec![[0.5 * lx, 0.5 * ly, 0.5 * lz]],
        flux_monitors: vec![DgtdFluxMonitorSpec { name: "mid_flux".into(), labels: vec![7] }],
        dispersion_poles: Default::default(),
        initial_cavity_mode: None,
        ..Default::default()
    };
    let mut p_dgtd = Problem::new(Study::Dgtd(dgtd_cfg), 2);
    p_dgtd.materials = p_driven.materials.clone();
    p_dgtd.pec = p_driven.pec.clone();
    p_dgtd.pmc = p_driven.pmc.clone();
    p_dgtd.waveports = p_driven.waveports.clone();

    let asm = assemble_dgtd(&mesh, &p_dgtd).expect("assemble_dgtd");
    let sol = solve_dgtd(&asm).expect("solve_dgtd");
    let res_dgtd = dgtd_fields(&mesh, &asm, &sol).expect("dgtd_fields");
    let pd_dgtd = DrivenPortData::of(&res_dgtd).expect("dgtd DrivenPortData");
    let td_data = nyx_ems::DgtdData::of(&res_dgtd).expect("dgtd DgtdData");
    assert_eq!(td_data.flux_monitors.len(), 1);
    assert_eq!(td_data.time_monitors.len(), 1);

    for (fi, &freq) in freqs_hz.iter().enumerate() {
        let s11_dr = c(pd_driven.s_matrix[fi][0][0][0], pd_driven.s_matrix[fi][0][0][1]).norm_sqr();
        let s21_dr = c(pd_driven.s_matrix[fi][1][0][0], pd_driven.s_matrix[fi][1][0][1]).norm_sqr();
        let s11_dg = c(pd_dgtd.s_matrix[fi][0][0][0], pd_dgtd.s_matrix[fi][0][0][1]).norm_sqr();
        let s21_dg = c(pd_dgtd.s_matrix[fi][1][0][0], pd_dgtd.s_matrix[fi][1][0][1]).norm_sqr();
        let err_s11 = (s11_dg - s11_dr).abs();
        let err_s21 = (s21_dg - s21_dr).abs();
        assert!(
            err_s11 < 1e-3 && err_s21 < 1e-3,
            "at f = {:.2} THz: |S11|^2 dgtd={:.5} vs driven={:.5} (err={:.2e}), |S21|^2 dgtd={:.5} vs driven={:.5} (err={:.2e})",
            freq / 1e12,
            s11_dg,
            s11_dr,
            err_s11,
            s21_dg,
            s21_dr,
            err_s21
        );
    }
}

#[test]
fn mixed_quad_hex_prism_and_p3_bases_use_diagonal_gll_and_local_mass_inverses() {
    // 1. 2D pure Quad mesh at p = 3 -> diagonal GLL mass matrix via femcore::linalg::matfree::Gll
    let mut cm2 = CellMesh::<2>::new();
    cm2.points = vec![[0.0, 0.0], [0.5e-6, 0.0], [1.0e-6, 0.0], [0.0, 0.5e-6], [0.5e-6, 0.5e-6], [1.0e-6, 0.5e-6]];
    cm2.push_cell(CellType::Quad, &[0, 1, 4, 3], 1);
    cm2.push_cell(CellType::Quad, &[1, 2, 5, 4], 1);
    cm2.push_facet(&[0, 3], 1);
    cm2.push_facet(&[2, 5], 2);
    cm2.push_facet(&[0, 1], 3);
    cm2.push_facet(&[1, 2], 3);
    cm2.push_facet(&[3, 4], 4);
    cm2.push_facet(&[4, 5], 4);
    let mesh_quad = Mesh::D2(cm2);

    let mut cfg =
        DgtdConfig { t_final: 2.0e-15, basis_order: Some(3), polarization: DgtdPolarization::Tm, ..Default::default() };
    let p_quad = Problem::new(Study::Dgtd(cfg.clone()), 3);
    let asm_quad = assemble_dgtd(&mesh_quad, &p_quad).expect("assemble p=3 quad dgtd");
    assert_eq!(asm_quad.basis_order, 3);
    assert!(asm_quad.cells.iter().all(|c| c.is_diagonal_mass && c.n_p == 16));
    let res_quad = run_dgtd(&mesh_quad, &p_quad).expect("run p=3 quad dgtd");
    assert_eq!(res_quad.study, "dgtd");

    // 2. 3D mixed Hex + Prism mesh at p = 2 -> Hex uses diagonal GLL mass (27 nodes), Prism uses local dense inverse (18 nodes)
    let mut cm3 = CellMesh::<3>::new();
    cm3.points = vec![
        [0.0, 0.0, 0.0],          // 0
        [1.0e-6, 0.0, 0.0],       // 1
        [1.0e-6, 1.0e-6, 0.0],    // 2
        [0.0, 1.0e-6, 0.0],       // 3
        [2.0e-6, 0.5e-6, 0.0],    // 4
        [0.0, 0.0, 1.0e-6],       // 5
        [1.0e-6, 0.0, 1.0e-6],    // 6
        [1.0e-6, 1.0e-6, 1.0e-6], // 7
        [0.0, 1.0e-6, 1.0e-6],    // 8
        [2.0e-6, 0.5e-6, 1.0e-6], // 9
    ];
    cm3.push_cell(CellType::Hex, &[0, 1, 2, 3, 5, 6, 7, 8], 1);
    cm3.push_cell(CellType::Prism, &[1, 4, 2, 6, 9, 7], 1);
    cm3.push_facet(&[0, 3, 8, 5], 1);
    cm3.push_facet(&[1, 4, 9, 6], 2);
    let mesh_3d = Mesh::D3(cm3);

    cfg.basis_order = Some(2);
    cfg.polarization = DgtdPolarization::Full;
    let p_3d = Problem::new(Study::Dgtd(cfg), 2);
    let asm_3d = assemble_dgtd(&mesh_3d, &p_3d).expect("assemble mixed Hex+Prism dgtd");
    assert!(asm_3d.cells[0].is_diagonal_mass && asm_3d.cells[0].n_p == 27);
    assert!(!asm_3d.cells[1].is_diagonal_mass && asm_3d.cells[1].n_p == 18);
    let res_3d = run_dgtd(&mesh_3d, &p_3d).expect("run mixed Hex+Prism dgtd");
    assert_eq!(res_3d.dim, 3);
}

#[test]
fn multiscale_nanogap_resonator_lts_achieves_over_5x_speedup_and_waveform_parity() {
    // Multiscale 2D nanogap resonator:
    // 2 fine columns of width dx_fine = 10 nm (4 triangles) adjacent to 120 coarse columns of width dx_coarse = 800 nm (240 triangles),
    // with ly = 1.6 um -> h_max / h_min > 65 (> 50).
    let mesh = build_multiscale_2d_mesh(2, 10.0e-9, 120, 800.0e-9, 1, 1.6e-6);
    let lx_tot = 2.0 * 10.0e-9 + 120.0 * 800.0e-9;
    let probe_pos = [0.35 * lx_tot, 0.8e-6, 0.0];

    let mut base_cfg = DgtdConfig {
        t_final: 4.0e-15,
        cfl: 0.75,
        flux: DgtdFlux::Upwind,
        time_stepper: DgtdTimeStepper::Lts,
        max_lts_tiers: 7,
        polarization: DgtdPolarization::Tm,
        basis_order: Some(2),
        pulse: DgtdPulse::Ricker { f0_hz: 1.5e14, t0_s: Some(1.5e-15), amplitude: 1.0 },
        freqs_hz: vec![1.5e14],
        excited_port: 1,
        absorbing_labels: vec![1, 2],
        probes: vec![probe_pos],
        flux_monitors: Vec::new(),
        dispersion_poles: Default::default(),
        initial_cavity_mode: None,
        ..Default::default()
    };

    let mut p_lts = Problem::new(Study::Dgtd(base_cfg.clone()), 2);
    p_lts.materials.set(1, "eps_r", 1.0);
    p_lts.materials.set(2, "eps_r", 1.0);
    p_lts.pmc = vec![3, 4];
    p_lts.excitations.push(ExcitationSpec {
        name: "dipole_src".into(),
        domains: vec![],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: false,
        excitation: Excitation::Dipole {
            position: [0.30 * lx_tot, 0.8e-6, 0.0],
            moment: [c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)],
            magnetic: false,
        },
    });

    let asm_lts = assemble_dgtd(&mesh, &p_lts).expect("assemble LTS multiscale");
    let h_ratio = asm_lts.h_max / asm_lts.h_min;
    assert!(h_ratio > 50.0, "expected h_max / h_min > 50, got {h_ratio:.2}");
    assert!(asm_lts.n_lts_tiers >= 5, "expected >= 5 LTS tiers, got {}", asm_lts.n_lts_tiers);

    let stride = 1usize << (asm_lts.n_lts_tiers - 1);
    let sol_lts = solve_dgtd(&asm_lts).expect("solve LTS multiscale");
    let res_lts = dgtd_fields(&mesh, &asm_lts, &sol_lts).expect("fields LTS multiscale");
    let data_lts = nyx_ems::DgtdData::of(&res_lts).expect("DgtdData LTS");

    // Uniform global dt_0 baseline (Rk4) over the exact same total physical duration n_macro_steps * stride * dt_0
    base_cfg.time_stepper = DgtdTimeStepper::Rk4;
    base_cfg.t_final = (sol_lts.n_macro_steps as f64) * (stride as f64) * asm_lts.dt_0;
    let mut p_uni = p_lts.clone();
    p_uni.study = Study::Dgtd(base_cfg);
    let asm_uni = assemble_dgtd(&mesh, &p_uni).expect("assemble uniform multiscale");
    let sol_uni = solve_dgtd(&asm_uni).expect("solve uniform multiscale");

    let work_speedup = data_lts.lts_work_speedup;
    let runtime_speedup = sol_uni.solve_time_s / sol_lts.solve_time_s.max(1e-9);
    assert!(
        work_speedup >= 5.0,
        "expected LTS element-update work speedup >= 5x, got {work_speedup:.2}x (updates: {} vs {})",
        sol_lts.total_element_updates,
        sol_lts.equivalent_uniform_updates
    );
    assert!(
        runtime_speedup >= 5.0,
        "expected LTS measured wall-clock speedup >= 5x, got {runtime_speedup:.2}x ({:.4}s vs {:.4}s)",
        sol_uni.solve_time_s,
        sol_lts.solve_time_s
    );

    // Waveform parity at the resonator probe at every macro step timestamp (< 1e-4 relative L2/Linf discrepancy)
    let wf_lts = &sol_lts.probe_e[0];
    let wf_uni = &sol_uni.probe_e[0];
    let n_cmp = wf_lts.len().min(wf_uni.len().div_ceil(stride));
    assert!(n_cmp > 10);
    let peak = (0..n_cmp).map(|k| wf_uni[k * stride][2].abs()).fold(0.0_f64, f64::max).max(1e-30);
    let max_diff = (0..n_cmp).map(|k| (wf_lts[k][2] - wf_uni[k * stride][2]).abs()).fold(0.0_f64, f64::max);
    let rel_err = max_diff / peak;
    assert!(rel_err < 1e-4, "expected LTS vs uniform global dt waveform parity < 1e-4, got {rel_err:.2e}");
}

#[test]
fn central_flux_conserves_cavity_energy_and_sheet2d_matches_fresnel_transmission() {
    // 1. Energy conservation with Central flux (alpha = 0) on a closed 2D PEC cavity initialized to TM_11 mode
    let sm_cav = femcore::testkit::rect(1.0e-6, 1.0e-6, 5, 5);
    let mesh_cav = Mesh::from_solver_mesh(&sm_cav).expect("2D PEC cavity mesh");
    let cfg_cav = DgtdConfig {
        t_final: 15.0e-15,
        cfl: 0.4,
        flux: DgtdFlux::Central,
        time_stepper: DgtdTimeStepper::Rk4,
        max_lts_tiers: 1,
        polarization: DgtdPolarization::Tm,
        basis_order: Some(2),
        pulse: DgtdPulse::Ricker { f0_hz: 1.0e14, t0_s: Some(0.0), amplitude: 0.0 },
        freqs_hz: vec![2.12e14],
        excited_port: 1,
        absorbing_labels: Vec::new(),
        probes: vec![[0.5e-6, 0.5e-6, 0.0]],
        flux_monitors: Vec::new(),
        dispersion_poles: Default::default(),
        initial_cavity_mode: Some([1, 1, 0]),
        ..Default::default()
    };
    let mut p_cav = Problem::new(Study::Dgtd(cfg_cav), 2);
    p_cav.materials.set(1, "eps_r", 1.0);
    let res_cav = run_dgtd(&mesh_cav, &p_cav).expect("run central-flux PEC cavity");
    let data_cav = nyx_ems::DgtdData::of(&res_cav).expect("cavity DgtdData");
    let e0 = *data_cav.energy_j.first().expect("initial energy");
    let ef = *data_cav.energy_j.last().expect("final energy");
    assert!(e0 > 0.0);
    let rel_energy_drift = (ef - e0).abs() / e0;
    assert!(
        rel_energy_drift < 1e-6,
        "expected central-flux PEC cavity energy conservation < 1e-6, got {rel_energy_drift:.2e} (e0={e0:.4e}, ef={ef:.4e})"
    );

    // 2. Penetrable zero-thickness Sheet2D (FR-35) normal-incidence Fresnel transmission parity (< 1e-3)
    // Sheet with surface conductivity sigma_2D = 1 / ETA0 placed at y = ly / 2 (facet label 7):
    // Exact Fresnel transmission T = |2 / (2 + ETA0 * sigma_2D)|^2 = (2 / 3)^2 = 4 / 9 = 0.444444...
    // Exact Fresnel reflection   R = | -ETA0 * sigma_2D / (2 + ETA0 * sigma_2D) |^2 = (1 / 3)^2 = 1 / 9 = 0.111111...
    let (lx, ly, lz) = (0.4e-6, 1.6e-6, 0.4e-6);
    let mesh_sh = build_3d_waveguide_mesh(lx, ly, lz, 1, 16, 1, (8, 16));
    let freqs_sh = vec![1.85e14, 1.934e14, 2.05e14];
    let cfg_sh = DgtdConfig {
        t_final: 55.0e-15,
        cfl: 0.75,
        flux: DgtdFlux::Upwind,
        time_stepper: DgtdTimeStepper::Lts,
        max_lts_tiers: 4,
        polarization: DgtdPolarization::Full,
        basis_order: Some(2),
        pulse: DgtdPulse::GaussianModulated {
            f0_hz: 1.934e14,
            fwidth_hz: 5.5e13,
            t0_s: None,
            amplitude: 1.0,
            phase_rad: 0.0,
        },
        freqs_hz: freqs_sh.clone(),
        excited_port: 1,
        absorbing_labels: Vec::new(),
        probes: Vec::new(),
        flux_monitors: Vec::new(),
        dispersion_poles: Default::default(),
        initial_cavity_mode: None,
        ..Default::default()
    };
    let mut p_sh = Problem::new(Study::Dgtd(cfg_sh), 2);
    p_sh.materials.set(1, "eps_r", 1.0);
    p_sh.materials.set(2, "eps_r", 1.0);
    p_sh.pec = vec![3, 6];
    p_sh.pmc = vec![1, 4];
    p_sh.waveports = vec![
        WavePortBc { name: "p1".into(), labels: vec![2], mode: 1, n_eff_guess: None, z0_ohm: None },
        WavePortBc { name: "p2".into(), labels: vec![5], mode: 1, n_eff_guess: None, z0_ohm: None },
    ];
    let sigma_s = 1.0_f64 / ETA0;
    p_sh.sheet_2d.push(Sheet2dBc {
        labels: vec![7],
        interface_domains: Vec::new(),
        sheet: femcore::material::Sheet2d::kubo_graphene(0.3 * 1.602_176_634e-19, 1e-13, 300.0, 1e6, 0.0).unwrap(),
        inline_sigma_uv: Some([[sigma_s, 0.0], [sigma_s, 0.0]]),
        c_ox: 0.0,
        v_cnp: 0.0,
        v_gate: 0.0,
        f_ref: 0.0,
    });

    let res_sh = run_dgtd(&mesh_sh, &p_sh).expect("run Sheet2D DGTD");
    let pd_sh = DrivenPortData::of(&res_sh).expect("Sheet2D DrivenPortData");
    let t_fresnel = (2.0_f64 / (2.0 + ETA0 * sigma_s)).powi(2);
    let r_fresnel = ((ETA0 * sigma_s) / (2.0 + ETA0 * sigma_s)).powi(2);
    for (fi, &freq) in freqs_sh.iter().enumerate() {
        let r_dg = c(pd_sh.s_matrix[fi][0][0][0], pd_sh.s_matrix[fi][0][0][1]).norm_sqr();
        let t_dg = c(pd_sh.s_matrix[fi][1][0][0], pd_sh.s_matrix[fi][1][0][1]).norm_sqr();
        assert!(
            (t_dg - t_fresnel).abs() < 1e-3 && (r_dg - r_fresnel).abs() < 1e-3,
            "at f={:.2} THz: Sheet2D T={:.5} (exact {:.5}, err={:.2e}), R={:.5} (exact {:.5}, err={:.2e})",
            freq / 1e12,
            t_dg,
            t_fresnel,
            (t_dg - t_fresnel).abs(),
            r_dg,
            r_fresnel,
            (r_dg - r_fresnel).abs()
        );
    }

    // 3. Dispersive Drude + Lorentz + GyrotropicDrude ADE poles
    let mut cfg_ade = p_sh.study.clone();
    if let Study::Dgtd(ref mut c_ade) = cfg_ade {
        c_ade.t_final = 10.0e-15;
        c_ade.dispersion_poles.insert(
            2,
            vec![
                DgtdDispersionPole::Drude { omega_p_rad_s: 2.0e14, gamma_rad_s: 1.0e13 },
                DgtdDispersionPole::Lorentz { omega_0_rad_s: 1.2e15, gamma_rad_s: 2.0e13, delta_eps: 0.25 },
                DgtdDispersionPole::GyrotropicDrude {
                    omega_p_rad_s: 1.0e14,
                    gamma_rad_s: 1.0e13,
                    omega_c_rad_s: [0.0, 5.0e13, 0.0],
                },
            ],
        );
    }
    let mut p_ade = p_sh.clone();
    p_ade.study = cfg_ade;
    p_ade.sheet_2d.clear();
    let res_ade = run_dgtd(&mesh_sh, &p_ade).expect("run ADE Drude/Lorentz/Gyrotropic DGTD");
    assert!(res_ade.solutions[0].fields[0].values.iter().all(|v| v.is_finite()));
}

#[test]
fn curved_tet_scatterer_and_subprocess_multi_thread_bitwise_determinism() {
    // 3D tetrahedral mesh with curved quadratic edge midpoints (CurvedTet) illuminated by a TFSF PlaneWave pulse
    let mut mesh = build_3d_waveguide_mesh(1.0e-6, 1.0e-6, 1.0e-6, 2, 2, 2, (1, 2));
    if let Mesh::D3(ref mut cm) = mesh {
        let p0 = cm.points[0];
        let p1 = cm.points[1];
        let mid = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1]) + 0.015e-6, 0.5 * (p0[2] + p1[2]) + 0.01e-6];
        cm.set_edge_nodes([([0, 1], mid)]);
    }

    let cfg = DgtdConfig {
        t_final: 8.0e-15,
        cfl: 0.65,
        flux: DgtdFlux::Upwind,
        time_stepper: DgtdTimeStepper::Lts,
        max_lts_tiers: 4,
        polarization: DgtdPolarization::Full,
        basis_order: Some(2),
        pulse: DgtdPulse::GaussianModulated {
            f0_hz: 1.934e14,
            fwidth_hz: 5.0e13,
            t0_s: None,
            amplitude: 1.0,
            phase_rad: 0.0,
        },
        freqs_hz: vec![1.934e14],
        excited_port: 1,
        absorbing_labels: vec![1, 2, 3, 4, 5, 6],
        probes: vec![[0.5e-6, 0.5e-6, 0.5e-6]],
        flux_monitors: vec![DgtdFluxMonitorSpec { name: "scatter_plane".into(), labels: vec![7] }],
        dispersion_poles: Default::default(),
        initial_cavity_mode: None,
        ..Default::default()
    };
    let mut problem = Problem::new(Study::Dgtd(cfg), 2);
    problem.materials.set(1, "eps_r", 1.0);
    problem.materials.set(2, "eps_r", 2.25);
    problem.excitations.push(ExcitationSpec {
        name: "pw_src".into(),
        domains: vec![1, 2],
        eps_bg: c(1.0, 0.0),
        mu_bg: 1.0,
        keep_scattered: false,
        excitation: Excitation::PlaneWave {
            k_hat: [0.0, 1.0, 0.0],
            pol: [c(0.0, 0.0), c(0.0, 0.0), c(1.0, 0.0)],
            amplitude: 1.0,
            phase: 0.0,
            origin: [0.0, 0.0, 0.0],
        },
    });

    // Subprocess worker mode: when NYX_FR42_SUBPROC_DUMP is set, serialize result to stdout and return
    if std::env::var("NYX_FR42_SUBPROC_DUMP").is_ok() {
        let res = run_dgtd(&mesh, &problem).expect("subprocess run_dgtd");
        let json = serde_json::to_string(&res).expect("serialize SolveResult");
        println!("NYX_DGTD_RESULT_JSON:{json}");
        return;
    }

    // Spawn two isolated subprocesses with NYX_THREADS=1 and NYX_THREADS=8 and verify bit-identical results
    let exe = std::env::current_exe().expect("current test binary");
    let run_with_threads = |n_threads: &str| -> nyx_physics::SolveResult {
        let out = std::process::Command::new(&exe)
            .arg("--exact")
            .arg("curved_tet_scatterer_and_subprocess_multi_thread_bitwise_determinism")
            .arg("--nocapture")
            .env("NYX_THREADS", n_threads)
            .env("NYX_FR42_SUBPROC_DUMP", "1")
            .output()
            .expect("spawn subprocess");
        assert!(out.status.success(), "subprocess failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8(out.stdout).expect("utf8 stdout");
        let line = stdout
            .lines()
            .find_map(|l| l.strip_prefix("NYX_DGTD_RESULT_JSON:"))
            .expect("found NYX_DGTD_RESULT_JSON line");
        serde_json::from_str(line).expect("deserialize SolveResult")
    };

    let r1 = run_with_threads("1");
    let r8 = run_with_threads("8");

    for (s1, s8) in r1.solutions.iter().zip(&r8.solutions) {
        for (f1, f8) in s1.fields.iter().zip(&s8.fields) {
            let b1: Vec<u64> = f1.values.iter().map(|v| v.to_bits()).collect();
            let b8: Vec<u64> = f8.values.iter().map(|v| v.to_bits()).collect();
            assert_eq!(b1, b8, "field '{}' must be bitwise identical between NYX_THREADS=1 and 8", f1.name);
        }
    }
    let d1 = r1.scalars.iter().find(|(k, _)| k == DGTD_DATA_KEY).expect("DGTD_DATA_KEY on r1");
    let d8 = r8.scalars.iter().find(|(k, _)| k == DGTD_DATA_KEY).expect("DGTD_DATA_KEY on r8");
    let td1 = nyx_ems::DgtdData::of(&r1).expect("td1");
    let td8 = nyx_ems::DgtdData::of(&r8).expect("td8");
    let e1_bits: Vec<u64> = td1.energy_j.iter().map(|v| v.to_bits()).collect();
    let e8_bits: Vec<u64> = td8.energy_j.iter().map(|v| v.to_bits()).collect();
    assert_eq!(e1_bits, e8_bits, "energy history must be bitwise identical across NYX_THREADS=1 and 8");
    let p1_bits: Vec<u64> = td1.time_monitors[0].e_t.iter().flatten().map(|v| v.to_bits()).collect();
    let p8_bits: Vec<u64> = td8.time_monitors[0].e_t.iter().flatten().map(|v| v.to_bits()).collect();
    assert_eq!(p1_bits, p8_bits, "probe time-series must be bitwise identical across NYX_THREADS=1 and 8");
    let f1_bits: Vec<u64> = td1.flux_monitors[0].flux_w.iter().map(|v| v.to_bits()).collect();
    let f8_bits: Vec<u64> = td8.flux_monitors[0].flux_w.iter().map(|v| v.to_bits()).collect();
    assert_eq!(f1_bits, f8_bits, "flux monitor spectrum must be bitwise identical across NYX_THREADS=1 and 8");
    let _ = (d1, d8);
}
