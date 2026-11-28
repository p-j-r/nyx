// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Verification suite for FR-43 Part B (`Study::Fdtd` and `Study::Dgtd` opt-in `wgpu` compute
//! kernels, `backend = "auto" | "cpu" | "gpu"` selection, CPU fallback diagnostics, and
//! discrete GPU performance benchmarks, ADR 146).

use std::collections::{BTreeMap, HashMap};

use femcore::mesh::SimplexMesh;
use nyx_ems::dgtd::{
    assemble_dgtd, run_dgtd, solve_dgtd, DgtdConfig, DgtdFlux, DgtdFluxMonitorSpec, DgtdPolarization, DgtdPulse,
    DgtdTimeStepper,
};
use nyx_ems::fdtd::{
    assemble_fdtd, run_fdtd, solve_fdtd, FdtdBcFace, FdtdBoundaries, FdtdMonitorSpec, FdtdSourceSpec, FdtdStudyConfig,
    FdtdWaveform, SubpixelMode,
};
use nyx_ems::problem::{Problem, Study, TdBackend, WavePortBc};
use nyx_physics::{Mesh, Physics, ProblemSpec};

#[cfg(feature = "gpu")]
use femcore::material::{Chi2, Chi3, Nonlinear, QpmSpec};
#[cfg(feature = "gpu")]
use nyx_ems::dgtd::{dgtd_fields, DgtdDispersionPole};
#[cfg(feature = "gpu")]
use nyx_ems::fdtd::gpu::EmsGpuContext;
#[cfg(feature = "gpu")]
use nyx_ems::fdtd::{fdtd_fields, DispersiveDomainSpec, SusceptibilityPole};
#[cfg(feature = "gpu")]
use nyx_ems::problem::TdPrecision;
#[cfg(feature = "gpu")]
use nyx_physics::{check_result_finite, DrivenPortData};

/// Conforming 6-tet-per-hex Kuhn 3D `SimplexMesh` over `[lo, hi]` with exterior facet labels:
/// `1 = -x`, `2 = +x`, `3 = -y`, `4 = +y`, `5 = -z`, `6 = +z`.
fn kuhn_3d_mesh(lo: [f64; 3], hi: [f64; 3], n: [usize; 3], mut domain_fn: impl FnMut([f64; 3]) -> i32) -> Mesh {
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
    let sm = SimplexMesh::new(3, pts, cells, cell_labels, facets, flabels);
    Mesh::from_solver_mesh(&sm).expect("Mesh::from_solver_mesh")
}

fn build_fdtd_problem(backend: TdBackend, total_steps: usize) -> (Mesh, Problem) {
    let mesh = kuhn_3d_mesh([-0.45e-6, 0.0, -0.30e-6], [0.45e-6, 2.0e-6, 0.30e-6], [3, 5, 2], |_| 1);

    let f0 = 1.934e14;
    let fwidth = 4.0e13;
    let waveform = FdtdWaveform::GaussianPulse {
        f0_hz: f0,
        fwidth_hz: fwidth,
        amplitude: 1.0,
        phase_rad: 0.0,
        peak_offset_sigma: 4.5,
    };
    let freqs = vec![1.88e14, 1.934e14, 1.98e14];
    let cpml4 = FdtdBcFace::Cpml { layers: 4, order: 3.5, kappa_max: 2.0, alpha_max: 0.02, target_reflection: 1e-6 };

    let mut cfg = FdtdStudyConfig {
        max_steps: total_steps,
        shutoff_energy_ratio: 0.0,
        subpixel: SubpixelMode::Staircased,
        port_z0_ohm: 50.0,
        boundaries: FdtdBoundaries {
            faces: [FdtdBcFace::Pec, FdtdBcFace::Pec, cpml4, cpml4, FdtdBcFace::Pec, FdtdBcFace::Pec],
        },
        sources: vec![FdtdSourceSpec::ModeSource {
            name: "port1_src".into(),
            axis: 1,
            position_m: 0.45e-6,
            direction: 1.0,
            mode_index: 1,
            n_eff_guess: Some(1.5),
            window_lo: None,
            window_hi: None,
            waveform,
        }],
        monitors: vec![
            FdtdMonitorSpec::ModeMonitor {
                name: "port1".into(),
                axis: 1,
                position_m: 0.55e-6,
                forward_sign: 1.0,
                mode_index: 1,
                n_eff_guess: Some(1.5),
                window_lo: None,
                window_hi: None,
                freqs_hz: freqs.clone(),
            },
            FdtdMonitorSpec::ModeMonitor {
                name: "port2".into(),
                axis: 1,
                position_m: 1.45e-6,
                forward_sign: 1.0,
                mode_index: 1,
                n_eff_guess: Some(1.5),
                window_lo: None,
                window_hi: None,
                freqs_hz: freqs.clone(),
            },
            FdtdMonitorSpec::FluxMonitor {
                name: "mid_flux".into(),
                axis: 1,
                position_m: 1.0e-6,
                normal_sign: 1.0,
                window_lo: None,
                window_hi: None,
                freqs_hz: freqs,
            },
            FdtdMonitorSpec::FieldMonitor {
                name: "center_vol".into(),
                lo: Some([-0.24e-6, 0.85e-6, -0.18e-6]),
                hi: Some([0.24e-6, 1.15e-6, 0.18e-6]),
                freqs_hz: vec![1.934e14],
            },
            FdtdMonitorSpec::FieldTimeMonitor {
                name: "center_probe".into(),
                position: [0.0, 1.0e-6, 0.0],
                component: 2,
                sample_stride: 2,
                prony_freq_range_hz: None,
                prony_start_time_s: None,
            },
        ],
        backend,
        gpu_min_cells: 250_000,
        ..Default::default()
    };
    cfg.grid.uniform_dl = Some([0.15e-6, 0.15e-6, 0.15e-6]);

    let mut problem = Problem::new(Study::Fdtd(Box::new(cfg)), 1);
    problem.materials.set(1, "eps_r", 2.25);
    (mesh, problem)
}

fn build_dgtd_problem(backend: TdBackend, stepper: DgtdTimeStepper, flux: DgtdFlux) -> (Mesh, Problem) {
    let mesh = kuhn_3d_mesh([-0.5e-6, 0.0, -0.35e-6], [0.5e-6, 1.6e-6, 0.35e-6], [2, 4, 2], |_| 1);

    let f0 = 1.934e14;
    let fwidth = 4.5e13;
    let freqs = vec![1.88e14, 1.934e14, 1.98e14];

    let cfg = DgtdConfig {
        t_final: 1.8e-15,
        cfl: 0.75,
        flux,
        time_stepper: stepper,
        max_lts_tiers: 1,
        polarization: DgtdPolarization::Full,
        basis_order: Some(1),
        pulse: DgtdPulse::GaussianModulated {
            f0_hz: f0,
            fwidth_hz: fwidth,
            t0_s: Some(0.6e-15),
            amplitude: 1.0,
            phase_rad: 0.0,
        },
        freqs_hz: freqs,
        excited_port: 1,
        absorbing_labels: Vec::new(),
        probes: vec![[0.0, 0.8e-6, 0.0]],
        flux_monitors: vec![DgtdFluxMonitorSpec { name: "out_flux".into(), labels: vec![4] }],
        dispersion_poles: BTreeMap::new(),
        initial_cavity_mode: None,
        backend,
        precision: Default::default(),
        nonlinear: false,
        gpu_min_dofs: 50_000,
    };

    let mut problem = Problem::new(Study::Dgtd(cfg), 1);
    problem.materials.set(1, "eps_r", 2.25);
    problem.pec = vec![1, 2];
    problem.pmc = vec![5, 6];
    problem.waveports = vec![
        WavePortBc { name: "p1".into(), labels: vec![3], mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
        WavePortBc { name: "p2".into(), labels: vec![4], mode: 1, n_eff_guess: None, z0_ohm: Some(50.0) },
    ];
    (mesh, problem)
}

#[test]
fn td_backend_cpu_and_auto_preserve_zero_ulp_determinism() {
    // 1. Study::Fdtd: "cpu" vs "auto" (below gpu_min_cells threshold) across NYX_THREADS=1 vs 8
    let (mesh_f, prob_cpu_f) = build_fdtd_problem(TdBackend::Cpu, 20);
    let (_, prob_auto_f) = build_fdtd_problem(TdBackend::Auto, 20);

    std::env::set_var("NYX_THREADS", "1");
    let asm_f = assemble_fdtd(&mesh_f, &prob_cpu_f).expect("assemble_fdtd cpu t1");
    let sol_cpu_t1 = solve_fdtd(&asm_f).expect("solve_fdtd cpu t1");

    std::env::set_var("NYX_THREADS", "8");
    let sol_cpu_t8 = solve_fdtd(&asm_f).expect("solve_fdtd cpu t8");

    let asm_fa = assemble_fdtd(&mesh_f, &prob_auto_f).expect("assemble_fdtd auto");
    let sol_auto = solve_fdtd(&asm_fa).expect("solve_fdtd auto");

    assert_eq!(sol_cpu_t1.backend_used, "cpu");
    assert_eq!(sol_auto.backend_used, "cpu");
    assert_eq!(sol_cpu_t1.fields_re.ex.as_slice(), sol_cpu_t8.fields_re.ex.as_slice());
    assert_eq!(sol_cpu_t1.fields_re.ey.as_slice(), sol_cpu_t8.fields_re.ey.as_slice());
    assert_eq!(sol_cpu_t1.fields_re.ez.as_slice(), sol_cpu_t8.fields_re.ez.as_slice());
    assert_eq!(sol_cpu_t1.fields_re.ex.as_slice(), sol_auto.fields_re.ex.as_slice());
    assert_eq!(sol_cpu_t1.peak_energy, sol_cpu_t8.peak_energy);
    assert_eq!(sol_cpu_t1.final_energy, sol_cpu_t8.final_energy);
    assert_eq!(sol_cpu_t1.peak_energy, sol_auto.peak_energy);
    assert_eq!(sol_cpu_t1.final_energy, sol_auto.final_energy);

    // 2. Study::Dgtd: "cpu" vs "auto" (below gpu_min_dofs threshold) across NYX_THREADS=1 vs 8
    let (mesh_d, prob_cpu_d) = build_dgtd_problem(TdBackend::Cpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
    let (_, prob_auto_d) = build_dgtd_problem(TdBackend::Auto, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);

    std::env::set_var("NYX_THREADS", "1");
    let asm_d1 = assemble_dgtd(&mesh_d, &prob_cpu_d).expect("assemble_dgtd cpu t1");
    let dsol_cpu_t1 = solve_dgtd(&asm_d1).expect("solve_dgtd cpu t1");

    std::env::set_var("NYX_THREADS", "8");
    let asm_d8 = assemble_dgtd(&mesh_d, &prob_cpu_d).expect("assemble_dgtd cpu t8");
    let dsol_cpu_t8 = solve_dgtd(&asm_d8).expect("solve_dgtd cpu t8");

    let asm_da = assemble_dgtd(&mesh_d, &prob_auto_d).expect("assemble_dgtd auto");
    let dsol_auto = solve_dgtd(&asm_da).expect("solve_dgtd auto");
    std::env::remove_var("NYX_THREADS");

    assert_eq!(dsol_cpu_t1.backend_used, "cpu");
    assert_eq!(dsol_auto.backend_used, "cpu");
    assert_eq!(dsol_cpu_t1.final_state.e, dsol_cpu_t8.final_state.e);
    assert_eq!(dsol_cpu_t1.final_state.h, dsol_cpu_t8.final_state.h);
    assert_eq!(dsol_cpu_t1.final_state.e, dsol_auto.final_state.e);
    assert_eq!(dsol_cpu_t1.energy_j, dsol_cpu_t8.energy_j);
    assert_eq!(dsol_cpu_t1.energy_j, dsol_auto.energy_j);

    // 3. ProblemSpec JSON decoding of "backend" parameter
    let spec: ProblemSpec = serde_json::from_value(serde_json::json!({
        "format": 1,
        "physics": "ems",
        "order": 1,
        "study": {
            "kind": "fdtd",
            "params": {
                "run_time": 1.0e-15,
                "max_steps": 5,
                "backend": "cpu"
            }
        },
        "materials": { "1": { "eps_r": 2.25 } },
        "bcs": [],
        "params": {}
    }))
    .expect("valid ProblemSpec with backend=cpu");
    nyx_ems::Ems.validate(&spec, &mesh_f).expect("validate fdtd backend=cpu");

    // 4. Without `--features gpu`, requesting `backend = "gpu"` must return a clear error
    #[cfg(not(feature = "gpu"))]
    {
        let (mesh_g, prob_gpu_f) = build_fdtd_problem(TdBackend::Gpu, 10);
        let err_f = run_fdtd(&mesh_g, &prob_gpu_f).expect_err("fdtd gpu without feature must fail");
        assert!(err_f.contains("--features gpu"), "expected error mentioning --features gpu, got: {err_f}");

        let (mesh_gd, prob_gpu_d) = build_dgtd_problem(TdBackend::Gpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
        let err_d = run_dgtd(&mesh_gd, &prob_gpu_d).expect_err("dgtd gpu without feature must fail");
        assert!(err_d.contains("--features gpu"), "expected error mentioning --features gpu, got: {err_d}");
    }
}

#[cfg(feature = "gpu")]
#[test]
fn fdtd_gpu_vs_cpu_waveguide_mode_port_and_dft_parity() {
    let (mesh, prob_cpu) = build_fdtd_problem(TdBackend::Cpu, 65);
    let (_, prob_gpu) = build_fdtd_problem(TdBackend::Gpu, 65);

    let asm_cpu = assemble_fdtd(&mesh, &prob_cpu).expect("assemble_fdtd cpu");
    let sol_cpu = solve_fdtd(&asm_cpu).expect("solve_fdtd cpu");
    let res_cpu = fdtd_fields(&asm_cpu, &sol_cpu).expect("fdtd_fields cpu");
    check_result_finite(&res_cpu).expect("cpu result finite");

    let asm_gpu = assemble_fdtd(&mesh, &prob_gpu).expect("assemble_fdtd gpu");
    let sol_gpu = solve_fdtd(&asm_gpu).expect("solve_fdtd gpu");
    let res_gpu = fdtd_fields(&asm_gpu, &sol_gpu).expect("fdtd_fields gpu");
    check_result_finite(&res_gpu).expect("gpu result finite");

    assert_eq!(sol_gpu.backend_used, "gpu", "expected GPU execution on NVIDIA T400");
    assert!(sol_gpu.gpu_adapter.is_some());
    assert!(sol_gpu.fallback_warning.is_none());

    // 1. DrivenPortData S-parameter parity (< 1e-4 in |S|^2 and complex S)
    let s_cpu = DrivenPortData::of(&res_cpu).expect("cpu DrivenPortData");
    let s_gpu = DrivenPortData::of(&res_gpu).expect("gpu DrivenPortData");

    assert_eq!(s_cpu.freqs_hz.len(), s_gpu.freqs_hz.len());
    for fi in 0..s_cpu.freqs_hz.len() {
        for p in 0..s_cpu.ports.len() {
            for q in 0..s_cpu.ports.len() {
                let sc = s_cpu.s_matrix[fi][p][q];
                let sg = s_gpu.s_matrix[fi][p][q];
                let mag2_c = sc[0] * sc[0] + sc[1] * sc[1];
                let mag2_g = sg[0] * sg[0] + sg[1] * sg[1];
                let diff_mag2 = (mag2_c - mag2_g).abs();
                let diff_c = ((sc[0] - sg[0]).powi(2) + (sc[1] - sg[1]).powi(2)).sqrt();
                assert!(
                    diff_mag2 < 1e-4 && diff_c < 1e-4,
                    "FDTD S_{}{} at freq[{}] mismatch: cpu={:?} gpu={:?} (diff_mag2={:.3e}, diff_c={:.3e})",
                    p + 1,
                    q + 1,
                    fi,
                    sc,
                    sg,
                    diff_mag2,
                    diff_c
                );
            }
        }
    }

    // 2. Planar DFT accumulator parity (< 1e-4 relative in f32)
    let mut dft_norm2_cpu = 0.0f64;
    let mut dft_diff2 = 0.0f64;
    for (ac, ag) in sol_cpu.planar_dfts.iter().zip(&sol_gpu.planar_dfts) {
        for (eh_c, eh_g) in ac.dft_eh.iter().zip(&ag.dft_eh) {
            for c_idx in 0..4 {
                dft_norm2_cpu += eh_c[c_idx].norm_sqr();
                dft_diff2 += (eh_c[c_idx] - eh_g[c_idx]).norm_sqr();
            }
        }
    }
    let rel_planar_dft = (dft_diff2 / dft_norm2_cpu.max(1e-30)).sqrt();
    assert!(rel_planar_dft < 1e-4, "FDTD planar DFT relative error {:.3e} exceeds 1e-4", rel_planar_dft);

    // 3. Volumetric FieldMonitor DFT parity (< 1e-4 relative in f32)
    let vc = &sol_cpu.vol_dfts[0];
    let vg = &sol_gpu.vol_dfts[0];
    let mut vol_norm2 = 0.0f64;
    let mut vol_diff2 = 0.0f64;
    for (eh_c, eh_g) in vc.eh_dft.iter().zip(&vg.eh_dft) {
        for c_idx in 0..6 {
            vol_norm2 += eh_c[c_idx].norm_sqr();
            vol_diff2 += (eh_c[c_idx] - eh_g[c_idx]).norm_sqr();
        }
    }
    let rel_vol_dft = (vol_diff2 / vol_norm2.max(1e-30)).sqrt();
    assert!(rel_vol_dft < 1e-4, "FDTD volume DFT relative error {:.3e} exceeds 1e-4", rel_vol_dft);

    // 4. FieldTimeMonitor probe trace parity (< 1e-4 relative)
    let tc = &sol_cpu.time_traces[0];
    let tg = &sol_gpu.time_traces[0];
    assert_eq!(tc.0.len(), tg.0.len());
    let mut t_norm2 = 0.0f64;
    let mut t_diff2 = 0.0f64;
    for (&vc_t, &vg_t) in tc.1.iter().zip(&tg.1) {
        t_norm2 += vc_t * vc_t;
        t_diff2 += (vc_t - vg_t).powi(2);
    }
    let rel_probe = (t_diff2 / t_norm2.max(1e-30)).sqrt();
    assert!(rel_probe < 1e-4, "FDTD time probe relative error {:.3e} exceeds 1e-4", rel_probe);

    // 5. Repeated GPU execution must be 0-ULP bit-identical (atomic-free deterministic reductions)
    let sol_gpu_2 = solve_fdtd(&asm_gpu).expect("solve_fdtd gpu repeat");
    assert_eq!(sol_gpu.peak_energy, sol_gpu_2.peak_energy);
    assert_eq!(sol_gpu.final_energy, sol_gpu_2.final_energy);
    assert_eq!(sol_gpu.vol_dfts[0].eh_dft, sol_gpu_2.vol_dfts[0].eh_dft);
    assert_eq!(sol_gpu.time_traces[0].1, sol_gpu_2.time_traces[0].1);

    // 6. Fallback test: FDTD with dispersive pole and backend="gpu" falls back to CPU with DiagnosticWarning
    let (_, mut prob_disp) = build_fdtd_problem(TdBackend::Gpu, 15);
    if let Study::Fdtd(ref mut cfg) = prob_disp.study {
        cfg.dispersive_domains.insert(
            1,
            DispersiveDomainSpec {
                eps_inf: 2.25,
                poles: vec![SusceptibilityPole { omega_0: 1.0e14, gamma: 1.0e13, strength: 1.0e28 }],
            },
        );
    }
    let res_fb = run_fdtd(&mesh, &prob_disp).expect("fdtd fallback to cpu");
    let diag = res_fb.diagnostics.as_ref().expect("fdtd diagnostics present");
    assert!(
        diag.warnings.iter().any(|w| w.code == "gpu_fallback_to_cpu"),
        "expected gpu_fallback_to_cpu warning in diagnostics, got: {:?}",
        diag.warnings
    );
}

#[cfg(feature = "gpu")]
#[test]
fn dgtd_gpu_vs_cpu_lserk4_and_rk4_parity_and_fallbacks() {
    for (stepper, flux) in [(DgtdTimeStepper::Lserk4, DgtdFlux::Upwind), (DgtdTimeStepper::Rk4, DgtdFlux::Central)] {
        let (mesh, prob_cpu) = build_dgtd_problem(TdBackend::Cpu, stepper, flux);
        let (_, prob_gpu) = build_dgtd_problem(TdBackend::Gpu, stepper, flux);

        let asm_cpu = assemble_dgtd(&mesh, &prob_cpu).expect("assemble_dgtd cpu");
        let sol_cpu = solve_dgtd(&asm_cpu).expect("solve_dgtd cpu");
        let res_cpu = dgtd_fields(&mesh, &asm_cpu, &sol_cpu).expect("dgtd_fields cpu");
        check_result_finite(&res_cpu).expect("dgtd cpu finite");

        let asm_gpu = assemble_dgtd(&mesh, &prob_gpu).expect("assemble_dgtd gpu");
        let sol_gpu = solve_dgtd(&asm_gpu).expect("solve_dgtd gpu");
        let res_gpu = dgtd_fields(&mesh, &asm_gpu, &sol_gpu).expect("dgtd_fields gpu");
        check_result_finite(&res_gpu).expect("dgtd gpu finite");

        assert_eq!(sol_gpu.backend_used, "gpu");
        assert!(sol_gpu.gpu_adapter.is_some());
        assert!(sol_gpu.fallback_warning.is_none());

        // 1. Final nodal state parity (< 1e-4 relative in f32)
        let mut e_norm2 = 0.0f64;
        let mut e_diff2 = 0.0f64;
        for (ec, eg) in sol_cpu.final_state.e.iter().zip(&sol_gpu.final_state.e) {
            for k in 0..3 {
                e_norm2 += ec[k] * ec[k];
                e_diff2 += (ec[k] - eg[k]).powi(2);
            }
        }
        let rel_e = (e_diff2 / e_norm2.max(1e-30)).sqrt();
        assert!(
            rel_e < 1e-4,
            "DGTD ({:?}, {:?}) final E state relative error {:.3e} exceeds 1e-4",
            stepper,
            flux,
            rel_e
        );

        // 2. Nodal DFT parity (< 1e-4 relative in f32)
        let mut dft_norm2 = 0.0f64;
        let mut dft_diff2 = 0.0f64;
        for fi in 0..sol_cpu.dft_e.len() {
            for (ec, eg) in sol_cpu.dft_e[fi].iter().zip(&sol_gpu.dft_e[fi]) {
                for k in 0..3 {
                    dft_norm2 += ec[k].norm_sqr();
                    dft_diff2 += (ec[k] - eg[k]).norm_sqr();
                }
            }
        }
        let rel_dft = (dft_diff2 / dft_norm2.max(1e-30)).sqrt();
        assert!(
            rel_dft < 1e-4,
            "DGTD ({:?}, {:?}) nodal DFT relative error {:.3e} exceeds 1e-4",
            stepper,
            flux,
            rel_dft
        );

        // 3. DrivenPortData S-parameter parity (< 1e-4)
        let s_cpu = DrivenPortData::of(&res_cpu).expect("dgtd cpu DrivenPortData");
        let s_gpu = DrivenPortData::of(&res_gpu).expect("dgtd gpu DrivenPortData");

        for fi in 0..s_cpu.freqs_hz.len() {
            for p in 0..s_cpu.ports.len() {
                for q in 0..s_cpu.ports.len() {
                    let sc = s_cpu.s_matrix[fi][p][q];
                    let sg = s_gpu.s_matrix[fi][p][q];
                    let diff_c = ((sc[0] - sg[0]).powi(2) + (sc[1] - sg[1]).powi(2)).sqrt();
                    assert!(
                        diff_c < 1e-4,
                        "DGTD ({:?}, {:?}) S_{}{} freq[{}] mismatch: cpu={:?} gpu={:?} (diff={:.3e})",
                        stepper,
                        flux,
                        p + 1,
                        q + 1,
                        fi,
                        sc,
                        sg,
                        diff_c
                    );
                }
            }
        }

        // 4. Energy history & probe time series parity (< 1e-4 relative)
        let mut en_norm2 = 0.0f64;
        let mut en_diff2 = 0.0f64;
        for (&wc, &wg) in sol_cpu.energy_j.iter().zip(&sol_gpu.energy_j) {
            en_norm2 += wc * wc;
            en_diff2 += (wc - wg).powi(2);
        }
        let rel_en = (en_diff2 / en_norm2.max(1e-30)).sqrt();
        assert!(
            rel_en < 1e-4,
            "DGTD ({:?}, {:?}) energy history relative error {:.3e} exceeds 1e-4",
            stepper,
            flux,
            rel_en
        );

        // 5. Repeated GPU execution must be 0-ULP bit-identical
        let sol_gpu_2 = solve_dgtd(&asm_gpu).expect("solve_dgtd gpu repeat");
        assert_eq!(sol_gpu.final_state.e, sol_gpu_2.final_state.e);
        assert_eq!(sol_gpu.final_state.h, sol_gpu_2.final_state.h);
        assert_eq!(sol_gpu.energy_j, sol_gpu_2.energy_j);
        assert_eq!(sol_gpu.probe_e, sol_gpu_2.probe_e);
    }

    // 6. Fallback test: DGTD with dispersive pole and backend="gpu" falls back to CPU with DiagnosticWarning
    let (mesh, mut prob_disp) = build_dgtd_problem(TdBackend::Gpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
    if let Study::Dgtd(ref mut cfg) = prob_disp.study {
        cfg.dispersion_poles.insert(1, vec![DgtdDispersionPole::Drude { omega_p_rad_s: 1.0e14, gamma_rad_s: 1.0e13 }]);
    }
    let res_fb = run_dgtd(&mesh, &prob_disp).expect("dgtd fallback to cpu");
    let diag = res_fb.diagnostics.as_ref().expect("dgtd diagnostics present");
    assert!(
        diag.warnings.iter().any(|w| w.code == "gpu_fallback_to_cpu"),
        "expected gpu_fallback_to_cpu warning in DGTD diagnostics, got: {:?}",
        diag.warnings
    );
}

#[cfg(feature = "gpu")]
#[test]
fn fdtd_gpu_nvidia_t400_million_cell_benchmark() {
    // 100 x 100 x 100 = 1,000,000 Yee cells (>= 10^6 cells requirement)
    let mesh = kuhn_3d_mesh([0.0; 3], [10.0e-6; 3], [2, 2, 2], |_| 1);
    let n_steps = 30usize;
    let waveform = FdtdWaveform::GaussianPulse {
        f0_hz: 1.934e14,
        fwidth_hz: 4.0e13,
        amplitude: 1.0,
        phase_rad: 0.0,
        peak_offset_sigma: 4.5,
    };
    let cpml6 = FdtdBcFace::Cpml { layers: 6, order: 3.5, kappa_max: 2.0, alpha_max: 0.02, target_reflection: 1e-6 };

    let make_prob = |backend: TdBackend| -> Problem {
        let mut cfg = FdtdStudyConfig {
            max_steps: n_steps,
            shutoff_energy_ratio: 0.0,
            subpixel: SubpixelMode::Staircased,
            boundaries: FdtdBoundaries { faces: [cpml6; 6] },
            sources: vec![FdtdSourceSpec::PointDipole {
                name: "center_dip".into(),
                position: [5.0e-6, 5.0e-6, 5.0e-6],
                component: 2,
                magnetic: false,
                waveform: waveform.clone(),
            }],
            monitors: vec![FdtdMonitorSpec::FieldTimeMonitor {
                name: "probe".into(),
                position: [5.3e-6, 5.0e-6, 5.0e-6],
                component: 2,
                sample_stride: 1,
                prony_freq_range_hz: None,
                prony_start_time_s: None,
            }],
            backend,
            gpu_min_cells: 250_000,
            ..Default::default()
        };
        cfg.grid.uniform_dl = Some([0.10e-6, 0.10e-6, 0.10e-6]);
        let mut p = Problem::new(Study::Fdtd(Box::new(cfg)), 1);
        p.materials.set(1, "eps_r", 2.25);
        p
    };

    // With 1,000,000 cells >= gpu_min_cells (250,000), `TdBackend::Auto` automatically selects the discrete GPU!
    let prob_auto_gpu = make_prob(TdBackend::Auto);
    let asm = assemble_fdtd(&mesh, &prob_auto_gpu).expect("assemble 1M cell fdtd");
    let n_cells = asm.grid.n_cells();
    assert!(n_cells >= 1_000_000, "expected >= 1,000,000 FDTD cells, got {n_cells}");

    let sol_gpu = solve_fdtd(&asm).expect("solve 1M cell fdtd on GPU");
    assert_eq!(sol_gpu.backend_used, "gpu", "TdBackend::Auto on 1M-cell grid must select discrete GPU");
    let adapter_name = sol_gpu.gpu_adapter.as_deref().unwrap_or("unknown");

    let prob_cpu = make_prob(TdBackend::Cpu);
    let asm_cpu = assemble_fdtd(&mesh, &prob_cpu).expect("assemble 1M cell fdtd cpu");
    let sol_cpu = solve_fdtd(&asm_cpu).expect("solve 1M cell fdtd on CPU");

    let total_cell_updates = (n_cells as f64) * (sol_gpu.n_steps_run as f64);
    let gpu_mcells_s = (total_cell_updates / 1e6) / sol_gpu.solve_seconds.max(1e-9);
    let cpu_mcells_s = (total_cell_updates / 1e6) / sol_cpu.solve_seconds.max(1e-9);
    let speedup = sol_cpu.solve_seconds / sol_gpu.solve_seconds.max(1e-9);

    eprintln!(
        "[FR-43B BENCHMARK] 3D FDTD grid: {}x{}x{} = {} cells, {} steps | Adapter: {} | GPU: {:.3}s ({:.1} Mcells/s) | CPU: {:.3}s ({:.1} Mcells/s) | Speedup: {:.2}x",
        asm.grid.nx,
        asm.grid.ny,
        asm.grid.nz,
        n_cells,
        sol_gpu.n_steps_run,
        adapter_name,
        sol_gpu.solve_seconds,
        gpu_mcells_s,
        sol_cpu.solve_seconds,
        cpu_mcells_s,
        speedup
    );

    // Verify probe waveform parity on the 1M-cell grid (< 1e-4 relative)
    let tc = &sol_cpu.time_traces[0];
    let tg = &sol_gpu.time_traces[0];
    let mut p_norm2 = 0.0f64;
    let mut p_diff2 = 0.0f64;
    for (&vc_t, &vg_t) in tc.1.iter().zip(&tg.1) {
        p_norm2 += vc_t * vc_t;
        p_diff2 += (vc_t - vg_t).powi(2);
    }
    let rel_probe = (p_diff2 / p_norm2.max(1e-30)).sqrt();
    assert!(rel_probe < 1e-4, "1M-cell FDTD GPU vs CPU probe relative error {:.3e} exceeds 1e-4", rel_probe);
}

#[cfg(feature = "gpu")]
fn rel_l2_f64_slice(ref_vals: &[f64], test_vals: &[f64]) -> f64 {
    assert_eq!(ref_vals.len(), test_vals.len());
    let mut n2 = 0.0_f64;
    let mut d2 = 0.0_f64;
    for (&r, &t) in ref_vals.iter().zip(test_vals) {
        n2 += r * r;
        d2 += (r - t).powi(2);
    }
    (d2 / n2.max(1e-60)).sqrt()
}

#[cfg(feature = "gpu")]
fn make_rich_nonlinear_block() -> Nonlinear {
    let mut d = [[0.0_f64; 6]; 3];
    d[2][2] = 3.0e-11; // d33 = 30 pm/V
    d[2][0] = 5.0e-12;
    d[0][4] = 5.0e-12;
    Nonlinear {
        chi2: Some(Chi2 { source: Vec::new(), wavelength: None, kleinman: true, broken_centrosymmetry: false, d }),
        qpm: Some(QpmSpec {
            sign: 1.0,
            period: Some(0.6e-6),
            axis: [0.0, 1.0, 0.0],
            duty_cycle: 0.5,
            phase_offset: 0.0,
        }),
        chi3: Some(Chi3 {
            chi3_xxxx: Some(2.0e-19),
            beta_tpa: Some(1.0e-11),
            sigma_fca: Some(1.45e-21),
            k_fcd: Some(-1.0e-27),
            tau_c: Some(1.0e-9),
            ..Default::default()
        }),
        ..Default::default()
    }
}

#[cfg(feature = "gpu")]
#[test]
fn fdtd_gpu_nonlinear_kerr_shg_tpa_parity_and_autoskip() {
    let ctx = EmsGpuContext::acquire_default().expect("GPU context required");
    let nl_block = make_rich_nonlinear_block();

    // 1. Default `nonlinear = false` auto-skips `[nonlinear]` so CPU/GPU stay on the fast linear path
    let (mesh, prob_lin) = build_fdtd_problem(TdBackend::Cpu, 45);
    let sol_lin = solve_fdtd(&assemble_fdtd(&mesh, &prob_lin).expect("asm lin")).expect("solve lin");

    let (_, mut prob_skip_cpu) = build_fdtd_problem(TdBackend::Cpu, 45);
    prob_skip_cpu.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
    let sol_skip_cpu =
        solve_fdtd(&assemble_fdtd(&mesh, &prob_skip_cpu).expect("asm skip cpu")).expect("solve skip cpu");
    assert_eq!(sol_skip_cpu.time_traces[0].1, sol_lin.time_traces[0].1, "nonlinear=false must 0-ULP match linear");
    assert_eq!(sol_skip_cpu.peak_energy, sol_lin.peak_energy);

    let (_, mut prob_skip_gpu) = build_fdtd_problem(TdBackend::Gpu, 45);
    prob_skip_gpu.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
    let sol_skip_gpu =
        solve_fdtd(&assemble_fdtd(&mesh, &prob_skip_gpu).expect("asm skip gpu")).expect("solve skip gpu");
    assert_eq!(sol_skip_gpu.backend_used, "gpu");
    assert!(sol_skip_gpu.fallback_warning.is_none());
    assert!(
        sol_skip_gpu.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f32")),
        "linear auto-skipped run should use fast f32 GPU path, got {:?}",
        sol_skip_gpu.gpu_adapter
    );
    assert!(rel_l2_f64_slice(&sol_lin.time_traces[0].1, &sol_skip_gpu.time_traces[0].1) < 1e-4);

    // 2. Opt-in `nonlinear = true` at strong guided optical field (`amplitude = 10.0` on 1W-normalized mode -> E ~ 3e8 V/m)
    let prep_strong = |backend: TdBackend, precision: TdPrecision, nonlinear: bool, with_mat: bool| -> Problem {
        let (_, mut p) = build_fdtd_problem(backend, 45);
        if with_mat {
            p.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
        }
        if let Study::Fdtd(ref mut cfg) = p.study {
            cfg.precision = precision;
            cfg.nonlinear = nonlinear;
            if let Some(FdtdSourceSpec::ModeSource {
                waveform: FdtdWaveform::GaussianPulse { ref mut amplitude, .. },
                ..
            }) = cfg.sources.first_mut()
            {
                *amplitude = 10.0;
            }
        }
        p
    };

    let sol_strong_lin =
        solve_fdtd(&assemble_fdtd(&mesh, &prep_strong(TdBackend::Cpu, TdPrecision::Auto, false, false)).unwrap())
            .unwrap();
    if ctx.supports_f64() {
        let sol_strong_lin_gpu_f64 =
            solve_fdtd(&assemble_fdtd(&mesh, &prep_strong(TdBackend::Gpu, TdPrecision::F64, false, false)).unwrap())
                .unwrap();
        let err_lin_f64 = rel_l2_f64_slice(&sol_strong_lin.time_traces[0].1, &sol_strong_lin_gpu_f64.time_traces[0].1);
        assert!(err_lin_f64 < 1e-12, "FDTD linear SHADER_F64 probe rel error {err_lin_f64:.3e} >= 1e-12");
    }
    let sol_cpu_nl =
        solve_fdtd(&assemble_fdtd(&mesh, &prep_strong(TdBackend::Cpu, TdPrecision::Auto, true, true)).unwrap())
            .unwrap();
    assert!(
        rel_l2_f64_slice(&sol_strong_lin.time_traces[0].1, &sol_cpu_nl.time_traces[0].1) > 1e-4,
        "strong nonlinear FDTD run should diverge measurably from linear baseline"
    );

    if ctx.supports_f64() {
        let asm_f64 = assemble_fdtd(&mesh, &prep_strong(TdBackend::Gpu, TdPrecision::Auto, true, true)).unwrap();
        let sol_gpu_f64 = solve_fdtd(&asm_f64).unwrap();
        assert_eq!(sol_gpu_f64.backend_used, "gpu");
        assert!(sol_gpu_f64.fallback_warning.is_none());
        assert!(sol_gpu_f64.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f64")));
        let err_probe = rel_l2_f64_slice(&sol_cpu_nl.time_traces[0].1, &sol_gpu_f64.time_traces[0].1);
        assert!(err_probe < 1e-7, "FDTD GPU SHADER_F64 nonlinear probe rel error {err_probe:.3e} >= 1e-7");
        let sol_gpu_f64_2 = solve_fdtd(&asm_f64).unwrap();
        assert_eq!(sol_gpu_f64.time_traces[0].1, sol_gpu_f64_2.time_traces[0].1);
    }

    let asm_f32 = assemble_fdtd(&mesh, &prep_strong(TdBackend::Gpu, TdPrecision::F32, true, true)).unwrap();
    let sol_gpu_f32 = solve_fdtd(&asm_f32).unwrap();
    assert_eq!(sol_gpu_f32.backend_used, "gpu");
    assert!(sol_gpu_f32.fallback_warning.is_none());
    assert!(sol_gpu_f32.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f32")));
    let err_probe_f32 = rel_l2_f64_slice(&sol_cpu_nl.time_traces[0].1, &sol_gpu_f32.time_traces[0].1);
    assert!(err_probe_f32 < 1e-3, "FDTD GPU normalized f32 nonlinear probe rel error {err_probe_f32:.3e} >= 1e-3");
    let sol_gpu_f32_2 = solve_fdtd(&asm_f32).unwrap();
    assert_eq!(sol_gpu_f32.time_traces[0].1, sol_gpu_f32_2.time_traces[0].1);
}

#[cfg(feature = "gpu")]
#[test]
fn dgtd_gpu_nonlinear_kerr_shg_tpa_and_nanoscale_energy_parity() {
    let ctx = EmsGpuContext::acquire_default().expect("GPU context required");
    let nl_block = make_rich_nonlinear_block();

    // 1. Default `nonlinear = false` auto-skips `[nonlinear]`
    let (mesh, prob_lin) = build_dgtd_problem(TdBackend::Cpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
    let sol_lin = solve_dgtd(&assemble_dgtd(&mesh, &prob_lin).unwrap()).unwrap();

    let (_, mut prob_skip_cpu) = build_dgtd_problem(TdBackend::Cpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
    prob_skip_cpu.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
    let sol_skip_cpu = solve_dgtd(&assemble_dgtd(&mesh, &prob_skip_cpu).unwrap()).unwrap();
    assert_eq!(sol_skip_cpu.final_state.e, sol_lin.final_state.e, "DGTD nonlinear=false must 0-ULP match linear");

    let (_, mut prob_skip_gpu) = build_dgtd_problem(TdBackend::Gpu, DgtdTimeStepper::Lserk4, DgtdFlux::Upwind);
    prob_skip_gpu.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
    let sol_skip_gpu = solve_dgtd(&assemble_dgtd(&mesh, &prob_skip_gpu).unwrap()).unwrap();
    assert_eq!(sol_skip_gpu.backend_used, "gpu");
    assert!(sol_skip_gpu.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f32")));

    // 2. Opt-in `nonlinear = true` across Lserk4 & Rk4 on GPU (`f64` Auto & normalized `f32`)
    for stepper in [DgtdTimeStepper::Lserk4, DgtdTimeStepper::Rk4] {
        let prep_nl = |backend: TdBackend, precision: TdPrecision, nonlinear: bool, with_mat: bool| -> Problem {
            let (_, mut p) = build_dgtd_problem(backend, stepper, DgtdFlux::Upwind);
            if with_mat {
                p.materials.aniso.entry(1).or_default().nonlinear = Some(nl_block.clone());
            }
            if let Study::Dgtd(ref mut cfg) = p.study {
                cfg.precision = precision;
                cfg.nonlinear = nonlinear;
                if let DgtdPulse::GaussianModulated { ref mut amplitude, .. } = cfg.pulse {
                    *amplitude = 1.0e3;
                }
            }
            p
        };

        let sol_lin_s =
            solve_dgtd(&assemble_dgtd(&mesh, &prep_nl(TdBackend::Cpu, TdPrecision::Auto, false, false)).unwrap())
                .unwrap();
        let flat_e = |s: &[[f64; 3]]| -> Vec<f64> { s.iter().flat_map(|v| v.iter().copied()).collect() };
        let sol_cpu_nl =
            solve_dgtd(&assemble_dgtd(&mesh, &prep_nl(TdBackend::Cpu, TdPrecision::Auto, true, true)).unwrap())
                .unwrap();
        let rel_delta = rel_l2_f64_slice(&flat_e(&sol_lin_s.final_state.e), &flat_e(&sol_cpu_nl.final_state.e));
        if ctx.supports_f64() {
            let sol_lin_gpu_f64 =
                solve_dgtd(&assemble_dgtd(&mesh, &prep_nl(TdBackend::Gpu, TdPrecision::F64, false, false)).unwrap())
                    .unwrap();
            let err_lin_e =
                rel_l2_f64_slice(&flat_e(&sol_lin_s.final_state.e), &flat_e(&sol_lin_gpu_f64.final_state.e));
            let err_lin_en = rel_l2_f64_slice(&sol_lin_s.energy_j, &sol_lin_gpu_f64.energy_j);
            assert!(
                err_lin_e < 1e-12 && err_lin_en < 1e-12,
                "DGTD ({stepper:?}) linear SHADER_F64 errors: E={err_lin_e:.3e}, W={err_lin_en:.3e}"
            );
        }
        assert!(
            rel_delta > 1e-4,
            "strong nonlinear DGTD should differ measurably from linear baseline (got {rel_delta:.3e})"
        );

        if ctx.supports_f64() {
            let asm_f64 = assemble_dgtd(&mesh, &prep_nl(TdBackend::Gpu, TdPrecision::Auto, true, true)).unwrap();
            let sol_gpu_f64 = solve_dgtd(&asm_f64).unwrap();
            assert_eq!(sol_gpu_f64.backend_used, "gpu");
            assert!(sol_gpu_f64.fallback_warning.is_none());
            assert!(sol_gpu_f64.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f64")));
            let err_e = rel_l2_f64_slice(&flat_e(&sol_cpu_nl.final_state.e), &flat_e(&sol_gpu_f64.final_state.e));
            let err_nc = rel_l2_f64_slice(&sol_cpu_nl.final_state.n_carrier, &sol_gpu_f64.final_state.n_carrier);
            let err_en = rel_l2_f64_slice(&sol_cpu_nl.energy_j, &sol_gpu_f64.energy_j);
            assert!(
                err_e < 1e-7 && err_nc < 1e-7 && err_en < 1e-7,
                "DGTD GPU f64 nonlinear errors: E={err_e:.3e}, Nc={err_nc:.3e}, W={err_en:.3e}"
            );
            let sol_gpu_f64_2 = solve_dgtd(&asm_f64).unwrap();
            assert_eq!(sol_gpu_f64.final_state.e, sol_gpu_f64_2.final_state.e);
            assert_eq!(sol_gpu_f64.final_state.n_carrier, sol_gpu_f64_2.final_state.n_carrier);
        }

        let asm_f32 = assemble_dgtd(&mesh, &prep_nl(TdBackend::Gpu, TdPrecision::F32, true, true)).unwrap();
        let sol_gpu_f32 = solve_dgtd(&asm_f32).unwrap();
        assert_eq!(sol_gpu_f32.backend_used, "gpu");
        assert!(sol_gpu_f32.fallback_warning.is_none());
        assert!(sol_gpu_f32.gpu_adapter.as_deref().is_some_and(|s| s.contains("precision=f32")));
        let err_e32 = rel_l2_f64_slice(&flat_e(&sol_cpu_nl.final_state.e), &flat_e(&sol_gpu_f32.final_state.e));
        let err_nc32 = rel_l2_f64_slice(&sol_cpu_nl.final_state.n_carrier, &sol_gpu_f32.final_state.n_carrier);
        assert!(
            err_e32 < 1e-3 && err_nc32 < 1e-3,
            "DGTD GPU normalized f32 nonlinear errors: E={err_e32:.3e}, Nc={err_nc32:.3e}"
        );
    }

    // 3. Nanoscale 50 nm DGTD cavity energy check in explicit f32: zero subnormal underflow!
    let nano_mesh = kuhn_3d_mesh([0.0; 3], [0.05e-6; 3], [2, 2, 2], |_| 1);
    let make_nano = |backend: TdBackend| -> Problem {
        let mut p = Problem::new(
            Study::Dgtd(DgtdConfig {
                t_final: 1.5e-16,
                cfl: 0.75,
                flux: DgtdFlux::Upwind,
                time_stepper: DgtdTimeStepper::Lserk4,
                max_lts_tiers: 1,
                basis_order: Some(1),
                initial_cavity_mode: Some([1, 1, 0]),
                backend,
                precision: TdPrecision::F32,
                ..Default::default()
            }),
            1,
        );
        p.materials.set(1, "eps_r", 2.25);
        p.pec = vec![1, 2, 3, 4, 5, 6];
        p
    };
    let sol_nano_cpu = solve_dgtd(&assemble_dgtd(&nano_mesh, &make_nano(TdBackend::Cpu)).unwrap()).unwrap();
    let sol_nano_gpu = solve_dgtd(&assemble_dgtd(&nano_mesh, &make_nano(TdBackend::Gpu)).unwrap()).unwrap();
    assert_eq!(sol_nano_gpu.backend_used, "gpu");
    assert!(
        sol_nano_gpu.energy_j.iter().all(|&w| w.is_finite() && w > 0.0),
        "nanoscale energy underflowed: {:?}",
        sol_nano_gpu.energy_j
    );
    let err_nano = rel_l2_f64_slice(&sol_nano_cpu.energy_j, &sol_nano_gpu.energy_j);
    assert!(err_nano < 1e-4, "nanoscale DGTD f32 energy rel error {err_nano:.3e} >= 1e-4");
}
