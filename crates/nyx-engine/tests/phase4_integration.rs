// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 cross-FR integration verification suite (`nyx-engine`, ADR 146 / `FR-43`):
//! 1. Default CPU `Engine::call` verification:
//!    - `0-ULP` bit-exact parity of `Study::Fdtd` and `Study::Dgtd` between `backend = "cpu"` and
//!      `backend = "auto"` (below discrete GPU offload thresholds, or without `--features gpu`).
//!    - Clear diagnostic error when `backend = "gpu"` is requested without `--features gpu`.
//! 2. Sparse format and mixed-precision iterative FEM parity:
//!    - Real 3D `nyx-dfs` electrostatic finite-element stiffness system solved via `Csr64<f64>`,
//!      `SellCSigma<f64>`, parallel Chebyshev `Amg<f64>`, and mixed-precision `AmgF32` driving
//!      outer `f64` `Cg` and `fgmres_vec` (`0-ULP` SpMV parity and `< 1e-10` relative residual).
//!    - 3D `H(curl)` Nédélec curl-curl (`Ams<f64>`, `AmsF32`, and complex-symmetric `AmsFgmresSolver`
//!      with `SmootherKind::Chebyshev`).
//! 3. Subprocess `NYX_THREADS=1` vs `NYX_THREADS=8` bitwise `0-ULP` determinism across `Engine`
//!    `Study::Fdtd`, `Study::Dgtd`, `Csr64`, `SellCSigma`, and Chebyshev `AmgF32` `cg_vec`.
//! 4. Real discrete GPU verification (`#[cfg(feature = "gpu")]`):
//!    - Verifies `femcore::linalg::gpu::GpuContext::shared()` and `nyx_ems::fdtd::gpu::EmsGpuContext::acquire_default()`
//!      share a single `wgpu::Device` and `wgpu::Queue` on the discrete GPU.
//!    - End-to-end `Engine` solving `Study::Fdtd` and `Study::Dgtd` with `backend = "gpu"` (`< 1e-4`
//!      S-parameter, flux, DFT, and energy parity vs CPU, plus `0-ULP` repeated-run GPU determinism
//!      and graceful `gpu_fallback_to_cpu` warning on multi-tier LTS).
//!    - `GpuAmgF32` + `GpuSellCSigmaF32` + outer `f64` PCG / `fgmres_vec` solving the 3D `nyx-dfs`
//!      FEM stiffness system on the exact same shared `GpuContext::shared()` interleaved with `Engine`
//!      time-domain GPU runs.

use femcore::fem::free_dofs;
use femcore::linalg::amg::{Amg, AmgF32, AmgOptions, SmootherKind};
use femcore::linalg::ams::{Ams, AmsF32, AmsFgmresOptions, AmsFgmresSolver, AmsOptions};
use femcore::linalg::csr::{Csr, CsrError};
use femcore::linalg::device::{cg_vec, fgmres_vec, tree_reduce_sum, CpuVec, Csr64, SellCSigma, Vector};
use femcore::linalg::driver::DirectSolver;
use femcore::linalg::gpu::{gpu_available, probe_gpu_device};
use femcore::linalg::ksp::{Cg, CgNorm, CgOptions, GmresOptions, Orthogonalization, PcSide, Tolerances};
use femcore::linalg::maxwell::nedelec_cube;
use femcore::linalg::scalar::{Scalar, C64};
use femcore::linalg::vecops;
use nyx_dfs::problem::{Bc, BcKind};
use nyx_dfs::{assemble_electrostatics, Problem as DfsProblem, Study as DfsStudy};
use nyx_ems::{DgtdData, FdtdMonitorData, C0};
use nyx_engine::Engine;
use nyx_physics::{canonical_json, discretise, DrivenPortData, SolveResult};
use serde_json::{json, Value};
use std::process::Command;

fn call(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

fn assert_f64_bits_eq(a: &[f64], b: &[f64], ctx: &str) {
    assert_eq!(a.len(), b.len(), "{ctx}: slice length mismatch ({} vs {})", a.len(), b.len());
    for (i, (&x, &y)) in a.iter().zip(b).enumerate() {
        assert_eq!(x.to_bits(), y.to_bits(), "{ctx}: ULP mismatch at [{i}]: {x:e} vs {y:e}");
    }
}

fn true_rel_residual<T: Scalar>(a: &Csr<T>, x: &[T], b: &[T]) -> f64 {
    let mut r = vec![T::ZERO; b.len()];
    a.spmv(x, &mut r);
    vecops::aypx(-T::ONE, b, &mut r);
    vecops::norm2(&r) / vecops::norm2(b).max(1e-300)
}

fn build_engine_fdtd(backend: &str, max_steps: usize) -> Engine {
    let w_x = 0.24e-6;
    let h_z = 0.24e-6;
    let y_lo = -0.50e-6;
    let y_hi = 1.50e-6;
    let f0 = C0 / 1.55e-6;
    let fwidth = 0.22 * f0;
    let freqs_hz = [C0 / 1.60e-6, f0, C0 / 1.50e-6];

    let fdtd_cfg = json!({
        "courant": 0.95,
        "max_steps": max_steps,
        "shutoff_energy_ratio": 0.0,
        "subpixel": "staircased",
        "port_z0_ohm": 50.0,
        "backend": backend,
        "gpu_min_cells": 250_000,
        "grid": {
            "uniform_dl": [0.08e-6, 0.08e-6, 0.08e-6]
        },
        "boundaries": {
            "faces": [
                { "kind": "pmc" },
                { "kind": "pmc" },
                { "kind": "cpml", "layers": 6, "order": 3.5, "kappa_max": 2.0, "alpha_max": 0.02, "target_reflection": 1e-6 },
                { "kind": "cpml", "layers": 6, "order": 3.5, "kappa_max": 2.0, "alpha_max": 0.02, "target_reflection": 1e-6 },
                { "kind": "pec" },
                { "kind": "pec" }
            ]
        },
        "sources": [{
            "kind": "mode_source",
            "name": "src_p1",
            "axis": 1,
            "position_m": -0.05e-6,
            "direction": 1.0,
            "mode_index": 1,
            "waveform": {
                "kind": "gaussian_pulse",
                "f0_hz": f0,
                "fwidth_hz": fwidth,
                "amplitude": 1.0,
                "phase_rad": 0.0,
                "peak_offset_sigma": 4.2
            }
        }],
        "monitors": [
            {
                "kind": "mode_monitor",
                "name": "p1",
                "axis": 1,
                "position_m": 0.15e-6,
                "forward_sign": 1.0,
                "mode_index": 1,
                "freqs_hz": freqs_hz
            },
            {
                "kind": "mode_monitor",
                "name": "p2",
                "axis": 1,
                "position_m": 0.85e-6,
                "forward_sign": 1.0,
                "mode_index": 1,
                "freqs_hz": freqs_hz
            },
            {
                "kind": "flux_monitor",
                "name": "mid_flux",
                "axis": 1,
                "position_m": 0.50e-6,
                "normal_sign": 1.0,
                "freqs_hz": freqs_hz
            },
            {
                "kind": "field_monitor",
                "name": "vol_mon",
                "lo": [0.04e-6, 0.35e-6, 0.04e-6],
                "hi": [0.20e-6, 0.65e-6, 0.20e-6],
                "freqs_hz": [f0]
            }
        ]
    });

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p4-fdtd", "dim": 3, "h": 0.12e-6 }));
    call(&mut e, "add_box", json!({ "name": "b_wg", "min": [0.0, y_lo, 0.0], "max": [w_x, y_hi, h_z] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_wg" }] }));
    call(&mut e, "set_mesh", json!({ "h": 0.12e-6, "sweep": [{ "thickness": y_hi - y_lo, "n": 8 }] }));
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "fdtd",
            "order": 1,
            "study_params": {
                "config_json": fdtd_cfg.to_string(),
                "backend": backend
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": { "eps_r": 2.25 } }));
    e
}

fn build_engine_dgtd(backend: &str, stepper: &str, max_lts_tiers: usize, t_final: f64) -> Engine {
    let w_x = 0.24e-6;
    let h_z = 0.24e-6;
    let length_y = 1.00e-6;
    let f0 = C0 / 1.55e-6;
    let fwidth = 0.24 * f0;
    let freqs_hz = [C0 / 1.60e-6, f0, C0 / 1.50e-6];

    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p4-dgtd", "dim": 3, "h": 0.12e-6 }));
    call(&mut e, "add_box", json!({ "name": "b_wg", "min": [0.0, 0.0, 0.0], "max": [w_x, length_y, h_z] }));
    call(&mut e, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_wg" }] }));
    call(&mut e, "set_mesh", json!({ "h": 0.12e-6, "sweep": [{ "thickness": length_y, "n": 6 }] }));
    call(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "dgtd",
            "order": 2,
            "study_params": {
                "t_final": t_final,
                "cfl": 0.75,
                "flux": "upwind",
                "time_stepper": stepper,
                "max_lts_tiers": max_lts_tiers,
                "f0": f0,
                "fwidth": fwidth,
                "freqs_hz": freqs_hz,
                "excited_port": 1,
                "flux_monitors_on": "extreme(+y)",
                "backend": backend
            }
        }),
    );
    call(&mut e, "set_material", json!({ "physics": "ems", "region": "core", "material": { "eps_r": 2.25 } }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
    call(&mut e, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
    );
    e
}

/// Assembles a real 3D `nyx-dfs` P2 electrostatic capacitor stiffness matrix reduced to free DOFs
/// (`k_ff`, `rhs_free`) so AMG / SellCSigma / Csr64 / GpuAmgF32 tests exercise genuine unstructured FEM operators.
fn assemble_3d_dfs_capacitor_free_system(h: f64, n_layers: usize) -> (Csr<f64>, Vec<f64>) {
    let mut e = Engine::new();
    call(&mut e, "project_new", json!({ "title": "p4-dfs-fem", "dim": 3, "h": h }));
    call(&mut e, "add_box", json!({ "name": "b_lo", "min": [0.0, 0.0, 0.0], "max": [1.0, 0.5, 1.0] }));
    call(&mut e, "add_box", json!({ "name": "b_hi", "min": [0.0, 0.5, 0.0], "max": [1.0, 1.0, 1.0] }));
    call(
        &mut e,
        "add_overlay",
        json!({
            "regions": [
                { "name": "dielectric_a", "input": "b_lo" },
                { "name": "dielectric_b", "input": "b_hi" }
            ]
        }),
    );
    call(
        &mut e,
        "set_mesh",
        json!({
            "h": h,
            "sweep": [
                { "thickness": 0.5, "n": n_layers },
                { "thickness": 0.5, "n": n_layers }
            ]
        }),
    );
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 2 }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "dielectric_a", "material": { "eps_r": 3.9 } }));
    call(&mut e, "set_material", json!({ "physics": "dfs", "region": "dielectric_b", "material": { "eps_r": 11.7 } }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "extreme(-y)", "kind": "potential", "params": { "value": 0.0 } }),
    );
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "extreme(+y)", "kind": "potential", "params": { "value": 1.25 } }),
    );
    call(&mut e, "solve_run", json!({}));
    let mesh: nyx_physics::Mesh = (**e.result_model("dfs").expect("built 3D dfs model").mesh()).clone();
    let disc = discretise(&mesh).expect("discretise 3D mesh");

    // Identify boundary labels on extreme(-y) (y = 0) and extreme(+y) (y = 1)
    let mut lbl_bot = std::collections::BTreeSet::new();
    let mut lbl_top = std::collections::BTreeSet::new();
    for fi in 0..disc.n_facets() {
        let fv = disc.facet(fi);
        let yc = (disc.point(fv[0] as usize)[1] + disc.point(fv[1] as usize)[1] + disc.point(fv[2] as usize)[1]) / 3.0;
        if yc.abs() < 1e-9 {
            lbl_bot.insert(disc.facet_labels[fi]);
        } else if (yc - 1.0).abs() < 1e-9 {
            lbl_top.insert(disc.facet_labels[fi]);
        }
    }

    let mut prob = DfsProblem::new(DfsStudy::Electrostatic, 2);
    for &dom in &disc.cell_labels {
        prob.materials.set(dom, "eps_r", if dom == 1 { 3.9 } else { 11.7 });
    }
    prob.bcs.push(Bc { labels: lbl_bot.into_iter().collect(), kind: BcKind::Potential(0.0) });
    prob.bcs.push(Bc { labels: lbl_top.into_iter().collect(), kind: BcKind::Potential(1.25) });

    let asm = assemble_electrostatics(&disc, &prob).expect("assemble 3D electrostatics");
    let n = asm.k.nrows();
    let free = free_dofs(n, &asm.fixed);
    let k_ff = asm.k.submatrix(&free, &free);
    let k_fx = asm.k.submatrix(&free, &asm.fixed);
    let kx = k_fx.matvec(&asm.vals);
    let rhs_free: Vec<f64> = free.iter().enumerate().map(|(a, &i)| asm.f[i] - kx[a]).collect();
    (k_ff, rhs_free)
}

fn collect_solve_result_bits(res: &SolveResult) -> Vec<u64> {
    let mut bits = Vec::new();
    if let Some(dp) = DrivenPortData::of(res) {
        for f_mat in &dp.s_matrix {
            for row in f_mat {
                for s in row {
                    bits.push(s[0].to_bits());
                    bits.push(s[1].to_bits());
                }
            }
        }
    }
    if let Some(md) = FdtdMonitorData::of(res) {
        for fm in &md.flux_monitors {
            for &v in &fm.raw_flux_w_s2 {
                bits.push(v.to_bits());
            }
        }
    }
    if let Some(dd) = DgtdData::of(res) {
        for &w in &dd.energy_j {
            bits.push(w.to_bits());
        }
        for fm in &dd.flux_monitors {
            for &v in &fm.raw_flux_w_s2 {
                bits.push(v.to_bits());
            }
        }
    }
    for sol in &res.solutions {
        for f in &sol.fields {
            for &v in &f.values {
                bits.push(v.to_bits());
            }
        }
        for cf in &sol.cell_fields {
            for &v in &cf.re {
                bits.push(v.to_bits());
            }
            for &v in &cf.im {
                bits.push(v.to_bits());
            }
        }
    }
    bits
}

#[test]
fn engine_td_cpu_vs_auto_zero_ulp_and_gpu_feature_diagnostics() {
    // 1. Study::Fdtd: "cpu" vs "auto" (below gpu_min_cells = 250_000) must be 0-ULP bit-identical
    let mut e_fdtd_cpu = build_engine_fdtd("cpu", 30);
    call(&mut e_fdtd_cpu, "solve_run", json!({}));
    let res_f_cpu = e_fdtd_cpu.result("ems").expect("fdtd cpu result").0;
    assert_eq!(res_f_cpu.diagnostics.as_ref().expect("fdtd diag").solver_kind, "yee_subpixel_fdtd");

    let mut e_fdtd_auto = build_engine_fdtd("auto", 30);
    call(&mut e_fdtd_auto, "solve_run", json!({}));
    let res_f_auto = e_fdtd_auto.result("ems").expect("fdtd auto result").0;
    assert_eq!(res_f_auto.diagnostics.as_ref().expect("fdtd auto diag").solver_kind, "yee_subpixel_fdtd");

    assert_eq!(
        collect_solve_result_bits(res_f_cpu),
        collect_solve_result_bits(res_f_auto),
        "Engine Study::Fdtd backend='cpu' vs 'auto' must be 0-ULP bit-identical"
    );

    // 2. Study::Dgtd: "cpu" vs "auto" (below gpu_min_dofs = 50_000) must be 0-ULP bit-identical
    let mut e_dgtd_cpu = build_engine_dgtd("cpu", "lserk4", 1, 2.5e-15);
    call(&mut e_dgtd_cpu, "solve_run", json!({}));
    let res_d_cpu = e_dgtd_cpu.result("ems").expect("dgtd cpu result").0;

    let mut e_dgtd_auto = build_engine_dgtd("auto", "lserk4", 1, 2.5e-15);
    call(&mut e_dgtd_auto, "solve_run", json!({}));
    let res_d_auto = e_dgtd_auto.result("ems").expect("dgtd auto result").0;

    assert_eq!(
        collect_solve_result_bits(res_d_cpu),
        collect_solve_result_bits(res_d_auto),
        "Engine Study::Dgtd backend='cpu' vs 'auto' must be 0-ULP bit-identical"
    );

    // 3. Without `--features gpu`, GPU probes report unavailable and `backend = "gpu"` returns a clear error
    #[cfg(not(feature = "gpu"))]
    {
        assert!(!gpu_available(), "default CPU build must report gpu_available() == false");
        assert!(probe_gpu_device().is_none(), "default CPU build must report probe_gpu_device() == None");

        let mut e_fdtd_gpu = build_engine_fdtd("gpu", 10);
        let err_f = e_fdtd_gpu.call("solve_run", &json!({})).expect_err("FDTD backend='gpu' without feature must fail");
        assert!(
            err_f.to_string().contains("--features gpu"),
            "expected FDTD diagnostic mentioning --features gpu, got: {err_f}"
        );

        let mut e_dgtd_gpu = build_engine_dgtd("gpu", "lserk4", 1, 1.5e-15);
        let err_d = e_dgtd_gpu.call("solve_run", &json!({})).expect_err("DGTD backend='gpu' without feature must fail");
        assert!(
            err_d.to_string().contains("--features gpu"),
            "expected DGTD diagnostic mentioning --features gpu, got: {err_d}"
        );
    }
}

#[test]
fn fem_dfs_and_ems_csr64_sell_c_sigma_chebyshev_and_mixed_amg_parity() {
    // 1. Assemble 3D P2 inhomogeneous dielectric capacitor FEM matrix via `nyx-dfs`
    let (k_ff, rhs_free) = assemble_3d_dfs_capacitor_free_system(0.22, 6);
    let n = k_ff.nrows();
    assert!(n >= 2_000, "expected >= 2,000 P2 free DOFs, got {n}");

    // Direct reference solution
    let direct = DirectSolver::new(&k_ff, true).expect("DirectSolver Cholesky");
    let u_ref = direct.solve(&rhs_free);
    assert!(true_rel_residual(&k_ff, &u_ref, &rhs_free) < 1e-12);

    // 2. Csr64<f64> lossless conversion + > u32::MAX indexing + 0-ULP SpMV parity
    let k64 = Csr64::from(&k_ff);
    assert!(!k64.exceeds_u32_limits());
    let k_back = Csr::try_from(&k64).expect("Csr64 -> Csr roundtrip");
    assert_eq!(k_ff, k_back);

    let huge_col = (u32::MAX as u64) + 2_048;
    let huge_csr = Csr64::from_triplets64(2, huge_col + 8, &[(0, 3, 1.5), (1, huge_col, -2.25)]).expect("Csr64 > 4B");
    assert!(huge_csr.exceeds_u32_limits());
    assert_eq!(Csr::try_from(&huge_csr), Err(CsrError::DimensionTooLarge));

    // 3. SellCSigma<f64> (C = 32, sigma = 128) 0-ULP SpMV parity against Csr<f64> and Csr64<f64>
    let sell = SellCSigma::from_csr(&k_ff).expect("SellCSigma::from_csr");
    let (mut y_csr, mut y_csr64, mut y_sell) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    k_ff.spmv(&u_ref, &mut y_csr);
    k64.spmv_checked(&u_ref, &mut y_csr64).expect("Csr64 spmv_checked");
    sell.spmv(&u_ref, &mut y_sell);
    assert_f64_bits_eq(&y_csr, &y_csr64, "Csr64 vs Csr FEM SpMV");
    assert_f64_bits_eq(&y_csr, &y_sell, "SellCSigma vs Csr FEM SpMV");

    // 4. Parallel Chebyshev Amg<f64> + mixed-precision AmgF32 over SellCSigma<f64>
    let amg_cheb =
        Amg::new(&k_ff, AmgOptions { smoother: SmootherKind::Chebyshev { degree: 3 }, ..Default::default() })
            .expect("Chebyshev Amg<f64>");
    let amg_f32 = AmgF32::new_chebyshev(&k_ff, 3).expect("Chebyshev AmgF32");

    let cg_opts = CgOptions {
        tol: Tolerances { rtol: 1e-11, max_it: 150, ..Default::default() },
        norm: CgNorm::Unpreconditioned,
        ..Default::default()
    };
    let b_vec = CpuVec::new(rhs_free.clone());
    let mut u_cheb = CpuVec::zeros(n);
    let info_cheb = cg_vec(&sell, &amg_cheb, &b_vec, &mut u_cheb, &cg_opts);
    assert!(info_cheb.reason.is_converged() && info_cheb.iterations <= 40, "{info_cheb:?}");
    assert!(true_rel_residual(&k_ff, u_cheb.as_slice(), &rhs_free) < 1e-10);

    let mut u_mixed_pcg = vec![0.0; n];
    let info_mixed = Cg::new(cg_opts).solve(&k_ff, &amg_f32, &rhs_free, &mut u_mixed_pcg);
    assert!(info_mixed.reason.is_converged() && info_mixed.iterations <= 42, "{info_mixed:?}");
    assert!(true_rel_residual(&k_ff, &u_mixed_pcg, &rhs_free) < 1e-10);

    let mut u_mixed_fg = CpuVec::zeros(n);
    let fg_opts = GmresOptions {
        tol: Tolerances { rtol: 1e-11, max_it: 120, ..Default::default() },
        restart: 30,
        side: PcSide::Right,
        flexible: true,
        orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
        ..Default::default()
    };
    let info_fg = fgmres_vec(&sell, &amg_f32, &b_vec, &mut u_mixed_fg, &fg_opts);
    assert!(info_fg.reason.is_converged(), "{info_fg:?}");
    assert!(true_rel_residual(&k_ff, u_mixed_fg.as_slice(), &rhs_free) < 1e-10);

    let mut diff_sq = 0.0f64;
    for (&u_m, &u_r) in u_mixed_pcg.iter().zip(&u_ref) {
        diff_sq += (u_m - u_r).powi(2);
    }
    let rel_u = diff_sq.sqrt() / vecops::norm2(&u_ref);
    assert!(rel_u < 1e-9, "AmgF32 FEM potential relative error vs Cholesky: {rel_u:e}");

    // 5. 3D Nedelec H(curl) Chebyshev Ams, mixed-precision AmsF32, and complex-symmetric Chebyshev AmsFgmresSolver
    let cube = nedelec_cube(5);
    let a_curl = Csr::axpby(1.0, &cube.stiffness, 1.0, &cube.mass);
    let b_curl: Vec<f64> = (0..a_curl.nrows()).map(|i| ((i as f64) * 0.37).sin()).collect();
    let ams_opts = AmsOptions {
        amg: AmgOptions { smoother: SmootherKind::Chebyshev { degree: 3 }, ..Default::default() },
        smoother_sweeps: 1,
        smoother: SmootherKind::Chebyshev { degree: 3 },
    };
    let ams_cheb = Ams::new(&a_curl, &cube.gradient, &cube.edge_vectors, &ams_opts).expect("Chebyshev Ams");
    let ams_f32 = AmsF32::new(&a_curl, &cube.gradient, &cube.edge_vectors, &ams_opts).expect("Chebyshev AmsF32");

    let curl_cg = CgOptions {
        tol: Tolerances { rtol: 1e-10, max_it: 120, ..Default::default() },
        norm: CgNorm::Unpreconditioned,
        ..Default::default()
    };
    let (mut x_ams64, mut x_ams32) = (vec![0.0; a_curl.nrows()], vec![0.0; a_curl.nrows()]);
    assert!(Cg::new(curl_cg.clone()).solve(&a_curl, &ams_cheb, &b_curl, &mut x_ams64).reason.is_converged());
    assert!(Cg::new(curl_cg).solve(&a_curl, &ams_f32, &b_curl, &mut x_ams32).reason.is_converged());
    assert!(true_rel_residual(&a_curl, &x_ams64, &b_curl) < 1e-9);
    assert!(true_rel_residual(&a_curl, &x_ams32, &b_curl) < 1e-9);

    let k_c = cube.stiffness.to_complex();
    let m_c = cube.mass.to_complex();
    let c_damp = Csr::axpby(C64::new(0.35, 0.0), &m_c, C64::ZERO, &m_c);
    let mut ams_fg_opts = AmsFgmresOptions { rtol: 1e-9, max_it: 200, restart: 35, ..Default::default() };
    ams_fg_opts.ams = ams_opts;
    let ams_solver = AmsFgmresSolver::new_shifted_complex_p1(
        &k_c,
        Some(&c_damp),
        &m_c,
        2.8,
        &cube.gradient,
        &cube.edge_vectors,
        ams_fg_opts,
    )
    .expect("Chebyshev AmsFgmresSolver<C64>");
    let rhs_c: Vec<C64> = (0..ams_solver.dim()).map(|i| C64::new(((i as f64) * 0.21).cos(), -0.15)).collect();
    let (_x_c, st_c) = ams_solver.solve_one(&rhs_c);
    assert!(st_c.all_converged && st_c.max_relative_residual < 1e-8, "{st_c:?}");
}

fn run_phase4_determinism_payload() -> String {
    let mut bits: Vec<u64> = Vec::new();

    // 1. Engine Study::Fdtd (CPU)
    let mut e_fdtd = build_engine_fdtd("cpu", 25);
    call(&mut e_fdtd, "solve_run", json!({}));
    bits.extend(collect_solve_result_bits(e_fdtd.result("ems").expect("det fdtd").0));

    // 2. Engine Study::Dgtd (CPU)
    let mut e_dgtd = build_engine_dgtd("cpu", "lserk4", 1, 2.0e-15);
    call(&mut e_dgtd, "solve_run", json!({}));
    bits.extend(collect_solve_result_bits(e_dgtd.result("ems").expect("det dgtd").0));

    // 3. 3D FEM stiffness + SellCSigma + Csr64 + Chebyshev AmgF32 + cg_vec + tree_reduce_sum
    let (k_ff, rhs_free) = assemble_3d_dfs_capacitor_free_system(0.22, 6);
    let sell = SellCSigma::from_csr(&k_ff).expect("det sell");
    let k64 = Csr64::from(&k_ff);
    let amg_f32 = AmgF32::new_chebyshev(&k_ff, 3).expect("det amg_f32");
    let mut u = CpuVec::zeros(k_ff.nrows());
    let opts = CgOptions {
        tol: Tolerances { rtol: 1e-10, max_it: 120, ..Default::default() },
        norm: CgNorm::Unpreconditioned,
        ..Default::default()
    };
    let info = cg_vec(&sell, &amg_f32, &CpuVec::new(rhs_free), &mut u, &opts);
    bits.push(info.iterations as u64);
    bits.push(tree_reduce_sum(u.as_slice()).to_bits());
    bits.push(u.norm2_sq().to_bits());

    let mut y64 = vec![0.0; k_ff.nrows()];
    k64.spmv_checked(u.as_slice(), &mut y64).expect("det spmv64");
    bits.push(tree_reduce_sum(&y64).to_bits());

    canonical_json(&json!({ "bits": bits })).to_string()
}

#[test]
fn phase4_subprocess_nyx_threads_1_vs_8_determinism() {
    if std::env::var("NYX_P4_DET_CHILD").as_deref() == Ok("1") {
        println!("NYX_P4_DET_PAYLOAD={}", run_phase4_determinism_payload());
        return;
    }

    let exe = std::env::current_exe().expect("current test executable");
    let run_with = |nt: &str| -> String {
        let out = Command::new(&exe)
            .args(["--exact", "phase4_subprocess_nyx_threads_1_vs_8_determinism", "--nocapture"])
            .env("NYX_P4_DET_CHILD", "1")
            .env("NYX_THREADS", nt)
            .output()
            .expect("spawn determinism child");
        assert!(
            out.status.success(),
            "NYX_THREADS={nt} child failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .find_map(|l| l.strip_prefix("NYX_P4_DET_PAYLOAD="))
            .expect("missing NYX_P4_DET_PAYLOAD line")
            .to_string()
    };

    let p1 = run_with("1");
    let p8 = run_with("8");
    assert_eq!(p1, p8, "Phase 4 outputs differed between NYX_THREADS=1 and NYX_THREADS=8");
}

#[cfg(feature = "gpu")]
mod gpu_integration {
    use super::*;
    use femcore::linalg::gpu::{GpuAmgF32, GpuContext, GpuSellCSigmaF32, GpuVecF32};
    use nyx_ems::fdtd::gpu::EmsGpuContext;
    use nyx_physics::check_result_finite;

    #[test]
    fn shared_gpu_context_end_to_end_engine_fdtd_dgtd_and_gpu_amg_pcg() {
        // 1. Verify shared `wgpu::Device` & `wgpu::Queue` between `femcore` and `nyx-ems`
        assert!(gpu_available(), "expected discrete wgpu compute adapter to be available");
        let probed = probe_gpu_device().expect("probe_gpu_device");
        let shared_ctx = GpuContext::shared().expect("GpuContext::shared");
        let ems_ctx = EmsGpuContext::acquire_default().expect("EmsGpuContext::acquire_default");

        assert_eq!(shared_ctx.info(), &probed);
        assert!(shared_ctx.info().is_discrete, "expected discrete GPU adapter, got {:?}", shared_ctx.info());
        assert!(ems_ctx.is_discrete(), "expected EmsGpuContext discrete adapter");
        assert_eq!(shared_ctx.info().name, ems_ctx.adapter_info.name);
        assert!(
            ems_ctx.device == *shared_ctx.wgpu_device(),
            "EmsGpuContext and femcore::GpuContext::shared() must share the exact same wgpu::Device"
        );
        assert!(
            ems_ctx.queue == *shared_ctx.wgpu_queue(),
            "EmsGpuContext and femcore::GpuContext::shared() must share the exact same wgpu::Queue"
        );

        // 2. End-to-end `Engine` solving `Study::Fdtd` with `backend = "gpu"` vs `backend = "cpu"`
        let mut e_fdtd_cpu = build_engine_fdtd("cpu", 70);
        call(&mut e_fdtd_cpu, "solve_run", json!({}));
        let res_f_cpu = e_fdtd_cpu.result("ems").expect("fdtd cpu").0;
        check_result_finite(res_f_cpu).expect("fdtd cpu finite");

        let mut e_fdtd_gpu1 = build_engine_fdtd("gpu", 70);
        call(&mut e_fdtd_gpu1, "solve_run", json!({}));
        let res_f_gpu1 = e_fdtd_gpu1.result("ems").expect("fdtd gpu 1").0;
        check_result_finite(res_f_gpu1).expect("fdtd gpu 1 finite");

        let diag_f = res_f_gpu1.diagnostics.as_ref().expect("fdtd gpu diagnostics");
        assert_eq!(diag_f.solver_kind, "yee_subpixel_fdtd_gpu");
        assert!(diag_f.warnings.is_empty(), "unexpected FDTD GPU fallback warnings: {:?}", diag_f.warnings);
        assert!(
            res_f_gpu1.log.iter().any(|l| l.contains(&shared_ctx.info().name)),
            "FDTD GPU log should cite shared adapter '{}': {:?}",
            shared_ctx.info().name,
            res_f_gpu1.log
        );

        let dp_f_cpu = DrivenPortData::of(res_f_cpu).expect("fdtd cpu DrivenPortData");
        let dp_f_gpu = DrivenPortData::of(res_f_gpu1).expect("fdtd gpu DrivenPortData");
        for fi in 0..dp_f_cpu.freqs_hz.len() {
            for p in 0..dp_f_cpu.ports.len() {
                for q in 0..dp_f_cpu.ports.len() {
                    let sc = dp_f_cpu.s_matrix[fi][p][q];
                    let sg = dp_f_gpu.s_matrix[fi][p][q];
                    let diff = (sc[0] - sg[0]).hypot(sc[1] - sg[1]);
                    assert!(
                        diff < 1e-4,
                        "Engine FDTD GPU vs CPU S_{}{} at freq[{fi}] mismatch: cpu={sc:?} gpu={sg:?} (diff={diff:.3e})",
                        p + 1,
                        q + 1
                    );
                }
            }
        }
        let md_cpu = FdtdMonitorData::of(res_f_cpu).expect("fdtd cpu monitor data");
        let md_gpu = FdtdMonitorData::of(res_f_gpu1).expect("fdtd gpu monitor data");
        let (mut flx_n2, mut flx_d2) = (0.0f64, 0.0f64);
        for (&fc, &fg) in md_cpu.flux_monitors[0].raw_flux_w_s2.iter().zip(&md_gpu.flux_monitors[0].raw_flux_w_s2) {
            flx_n2 += fc * fc;
            flx_d2 += (fc - fg).powi(2);
        }
        assert!((flx_d2 / flx_n2.max(1e-30)).sqrt() < 1e-4, "FDTD FluxMonitor GPU vs CPU mismatch");

        // 3. Interleaved `GpuAmgF32` + `GpuSellCSigmaF32` + outer `f64` PCG on the SAME `shared_ctx`
        let (k_ff, rhs_free) = assemble_3d_dfs_capacitor_free_system(0.17, 7);
        let n_fem = k_ff.nrows();
        assert!(n_fem >= 4_000, "expected >= 4,000 3D P2 FEM free DOFs, got {n_fem}");

        let direct = DirectSolver::new(&k_ff, true).expect("DirectSolver FEM reference");
        let u_direct = direct.solve(&rhs_free);

        let gpu_amg = GpuAmgF32::new_chebyshev(&shared_ctx, &k_ff, 3).expect("GpuAmgF32 on shared_ctx");
        assert!(gpu_amg.context().wgpu_device() == &ems_ctx.device);
        let gpu_sell = GpuSellCSigmaF32::from_csr(&shared_ctx, &k_ff).expect("GpuSellCSigmaF32 on shared_ctx");

        let cg_opts = CgOptions {
            tol: Tolerances { rtol: 1e-10, max_it: 150, ..Default::default() },
            norm: CgNorm::Unpreconditioned,
            record_history: true,
            ..Default::default()
        };
        let mut u_gpu_pcg1 = vec![0.0; n_fem];
        let info_pcg1 = Cg::new(cg_opts.clone()).solve(&k_ff, &gpu_amg, &rhs_free, &mut u_gpu_pcg1);
        assert!(info_pcg1.reason.is_converged() && info_pcg1.iterations <= 45, "{info_pcg1:?}");
        let rel_res_pcg = true_rel_residual(&k_ff, &u_gpu_pcg1, &rhs_free);
        assert!(rel_res_pcg < 1e-9, "GpuAmgF32 + f64 PCG FEM residual {rel_res_pcg:e}");

        let err_vs_direct = u_gpu_pcg1.iter().zip(&u_direct).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
            / vecops::norm2(&u_direct);
        assert!(err_vs_direct < 1e-8, "GpuAmgF32 + f64 PCG vs Cholesky FEM error {err_vs_direct:e}");

        // Also run device-resident `cg_vec<GpuVecF32>` on the shared device
        let b_dev = GpuVecF32::from_f64_slice(&shared_ctx, &rhs_free);
        let mut x_dev = GpuVecF32::zeros(&shared_ctx, n_fem);
        let dev_info = cg_vec(
            &gpu_sell,
            &gpu_amg,
            &b_dev,
            &mut x_dev,
            &CgOptions { tol: Tolerances { rtol: 1e-5, max_it: 80, ..Default::default() }, ..Default::default() },
        );
        assert!(dev_info.reason.is_converged(), "{dev_info:?}");
        assert!(true_rel_residual(&k_ff, &x_dev.to_host(), &rhs_free) < 1e-4);

        // 4. End-to-end `Engine` solving `Study::Dgtd` with `backend = "gpu"` vs `backend = "cpu"` on shared context
        let mut e_dgtd_cpu = build_engine_dgtd("cpu", "lserk4", 1, 4.5e-15);
        call(&mut e_dgtd_cpu, "solve_run", json!({}));
        let res_d_cpu = e_dgtd_cpu.result("ems").expect("dgtd cpu").0;
        check_result_finite(res_d_cpu).expect("dgtd cpu finite");

        let mut e_dgtd_gpu1 = build_engine_dgtd("gpu", "lserk4", 1, 4.5e-15);
        call(&mut e_dgtd_gpu1, "solve_run", json!({}));
        let res_d_gpu1 = e_dgtd_gpu1.result("ems").expect("dgtd gpu 1").0;
        check_result_finite(res_d_gpu1).expect("dgtd gpu 1 finite");

        let diag_d = res_d_gpu1.diagnostics.as_ref().expect("dgtd gpu diagnostics");
        assert_eq!(diag_d.solver_kind, "dgtd_gpu");
        assert!(diag_d.warnings.is_empty(), "unexpected DGTD GPU fallback warnings: {:?}", diag_d.warnings);

        let dp_d_cpu = DrivenPortData::of(res_d_cpu).expect("dgtd cpu DrivenPortData");
        let dp_d_gpu = DrivenPortData::of(res_d_gpu1).expect("dgtd gpu DrivenPortData");
        for fi in 0..dp_d_cpu.freqs_hz.len() {
            for p in 0..dp_d_cpu.ports.len() {
                for q in 0..dp_d_cpu.ports.len() {
                    let sc = dp_d_cpu.s_matrix[fi][p][q];
                    let sg = dp_d_gpu.s_matrix[fi][p][q];
                    let diff = (sc[0] - sg[0]).hypot(sc[1] - sg[1]);
                    assert!(
                        diff < 1e-4,
                        "Engine DGTD GPU vs CPU S_{}{} at freq[{fi}] mismatch: cpu={sc:?} gpu={sg:?} (diff={diff:.3e})",
                        p + 1,
                        q + 1
                    );
                }
            }
        }
        let dd_cpu = DgtdData::of(res_d_cpu).expect("dgtd cpu data");
        let dd_gpu = DgtdData::of(res_d_gpu1).expect("dgtd gpu data");
        let (mut en_n2, mut en_d2) = (0.0f64, 0.0f64);
        for (&wc, &wg) in dd_cpu.energy_j.iter().zip(&dd_gpu.energy_j) {
            en_n2 += wc * wc;
            en_d2 += (wc - wg).powi(2);
        }
        assert!((en_d2 / en_n2.max(1e-30)).sqrt() < 1e-4, "Engine DGTD energy history GPU vs CPU mismatch");

        // 5. Repeat FDTD GPU, DGTD GPU, and GpuAmgF32 PCG after interleaving on the shared queue -> 0-ULP determinism!
        let mut e_fdtd_gpu2 = build_engine_fdtd("gpu", 70);
        call(&mut e_fdtd_gpu2, "solve_run", json!({}));
        assert_eq!(
            collect_solve_result_bits(res_f_gpu1),
            collect_solve_result_bits(e_fdtd_gpu2.result("ems").expect("fdtd gpu 2").0),
            "Repeated Engine FDTD GPU execution on shared GpuContext must be 0-ULP deterministic"
        );

        let mut e_dgtd_gpu2 = build_engine_dgtd("gpu", "lserk4", 1, 4.5e-15);
        call(&mut e_dgtd_gpu2, "solve_run", json!({}));
        assert_eq!(
            collect_solve_result_bits(res_d_gpu1),
            collect_solve_result_bits(e_dgtd_gpu2.result("ems").expect("dgtd gpu 2").0),
            "Repeated Engine DGTD GPU execution on shared GpuContext must be 0-ULP deterministic"
        );

        let mut u_gpu_pcg2 = vec![0.0; n_fem];
        let info_pcg2 = Cg::new(cg_opts).solve(&k_ff, &gpu_amg, &rhs_free, &mut u_gpu_pcg2);
        assert_eq!(info_pcg1.iterations, info_pcg2.iterations);
        assert_f64_bits_eq(&info_pcg1.history, &info_pcg2.history, "GpuAmgF32 PCG history determinism");
        assert_f64_bits_eq(&u_gpu_pcg1, &u_gpu_pcg2, "GpuAmgF32 PCG solution 0-ULP determinism");

        // 6. Auto-skip `[nonlinear]` by default (`nonlinear = false`) and GPU nonlinear execution (`nonlinear = true`) on silicon_nitride
        let mut e_sin_lin = build_engine_dgtd("gpu", "lserk4", 1, 2.0e-15);
        call(
            &mut e_sin_lin,
            "set_material",
            json!({ "physics": "ems", "region": "core", "material": "silicon_nitride" }),
        );
        call(&mut e_sin_lin, "solve_run", json!({}));
        let res_sin_lin = e_sin_lin.result("ems").expect("dgtd auto-skip result").0;
        let diag_sin_lin = res_sin_lin.diagnostics.as_ref().expect("dgtd auto-skip diagnostics");
        assert_eq!(diag_sin_lin.solver_kind, "dgtd_gpu");
        assert!(
            !diag_sin_lin.warnings.iter().any(|w| w.code == "gpu_fallback_to_cpu"),
            "nonlinear=false on silicon_nitride should run on GPU without gpu_fallback_to_cpu, got: {:?}",
            diag_sin_lin.warnings
        );

        let build_engine_dgtd_nl = |backend: &str, precision: &str| -> Engine {
            let f0 = C0 / 1.55e-6;
            let mut eng = Engine::new();
            call(&mut eng, "project_new", json!({ "title": "p4-dgtd-nl", "dim": 3, "h": 0.12e-6 }));
            call(
                &mut eng,
                "add_box",
                json!({ "name": "b_wg", "min": [0.0, 0.0, 0.0], "max": [0.24e-6, 1.00e-6, 0.24e-6] }),
            );
            call(&mut eng, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_wg" }] }));
            call(&mut eng, "set_mesh", json!({ "h": 0.12e-6, "sweep": [{ "thickness": 1.00e-6, "n": 6 }] }));
            call(
                &mut eng,
                "add_physics",
                json!({
                    "physics": "ems",
                    "study": "dgtd",
                    "order": 2,
                    "study_params": {
                        "t_final": 2.0e-15,
                        "cfl": 0.75,
                        "flux": "upwind",
                        "time_stepper": "lserk4",
                        "max_lts_tiers": 1,
                        "f0": f0,
                        "fwidth": 0.24 * f0,
                        "freqs_hz": [C0 / 1.60e-6, f0, C0 / 1.50e-6],
                        "excited_port": 1,
                        "flux_monitors_on": "extreme(+y)",
                        "backend": backend,
                        "precision": precision,
                        "nonlinear": true
                    }
                }),
            );
            call(
                &mut eng,
                "set_material",
                json!({ "physics": "ems", "region": "core", "material": "silicon_nitride" }),
            );
            call(&mut eng, "add_condition", json!({ "physics": "ems", "on": "extreme(-x)", "kind": "pmc" }));
            call(&mut eng, "add_condition", json!({ "physics": "ems", "on": "extreme(+x)", "kind": "pmc" }));
            call(
                &mut eng,
                "add_condition",
                json!({ "physics": "ems", "on": "extreme(-y)", "kind": "waveport", "params": { "name": "p1", "mode": 1 } }),
            );
            call(
                &mut eng,
                "add_condition",
                json!({ "physics": "ems", "on": "extreme(+y)", "kind": "waveport", "params": { "name": "p2", "mode": 1 } }),
            );
            eng
        };

        let mut e_sin_nl_cpu = build_engine_dgtd_nl("cpu", "auto");
        call(&mut e_sin_nl_cpu, "solve_run", json!({}));
        let res_sin_nl_cpu = e_sin_nl_cpu.result("ems").expect("dgtd nl cpu").0;

        let mut e_sin_nl_gpu = build_engine_dgtd_nl("gpu", "auto");
        call(&mut e_sin_nl_gpu, "solve_run", json!({}));
        let res_sin_nl_gpu = e_sin_nl_gpu.result("ems").expect("dgtd nl gpu").0;
        let diag_sin_nl_gpu = res_sin_nl_gpu.diagnostics.as_ref().expect("dgtd nl gpu diagnostics");
        assert_eq!(diag_sin_nl_gpu.solver_kind, "dgtd_gpu");
        assert!(
            !diag_sin_nl_gpu.warnings.iter().any(|w| w.code == "gpu_fallback_to_cpu"),
            "nonlinear=true on silicon_nitride should run on GPU without gpu_fallback_to_cpu, got: {:?}",
            diag_sin_nl_gpu.warnings
        );
        let d_nl_cpu = DgtdData::of(res_sin_nl_cpu).expect("cpu DgtdData");
        let d_nl_gpu = DgtdData::of(res_sin_nl_gpu).expect("gpu DgtdData");
        let (mut en_n2_nl, mut en_d2_nl) = (0.0f64, 0.0f64);
        for (&wc, &wg) in d_nl_cpu.energy_j.iter().zip(&d_nl_gpu.energy_j) {
            en_n2_nl += wc * wc;
            en_d2_nl += (wc - wg).powi(2);
        }
        assert!((en_d2_nl / en_n2_nl.max(1e-30)).sqrt() < 1e-4, "Engine nonlinear DGTD energy GPU vs CPU mismatch");

        // 7. Graceful Engine CPU fallback warning when dispersive ADE poles are active with `backend = "gpu"`
        let f0 = C0 / 1.55e-6;
        let fwidth = 0.24 * f0;
        let freqs_hz = [C0 / 1.60e-6, f0, C0 / 1.50e-6];
        let fdtd_disp_cfg = json!({
            "max_steps": 15,
            "shutoff_energy_ratio": 0.0,
            "subpixel": "staircased",
            "grid": { "uniform_dl": [0.06e-6, 0.06e-6, 0.06e-6] },
            "boundaries": { "faces": [
                { "kind": "pec" }, { "kind": "pec" },
                { "kind": "cpml", "layers": 4, "order": 3.5, "kappa_max": 2.0, "alpha_max": 0.02, "target_reflection": 1e-6 },
                { "kind": "cpml", "layers": 4, "order": 3.5, "kappa_max": 2.0, "alpha_max": 0.02, "target_reflection": 1e-6 },
                { "kind": "pec" }, { "kind": "pec" }
            ] },
            "dispersive_domains": {
                "1": { "eps_inf": 2.25, "poles": [{ "omega_0": 1.0e14, "gamma": 1.0e13, "strength": 1.0e28 }] }
            },
            "sources": [{
                "kind": "mode_source", "name": "s1", "axis": 1, "position_m": 0.30e-6, "direction": 1.0, "mode_index": 1, "n_eff_guess": 1.5,
                "waveform": { "kind": "gaussian_pulse", "f0_hz": f0, "fwidth_hz": fwidth, "amplitude": 1.0, "phase_rad": 0.0, "peak_offset_sigma": 4.2 }
            }],
            "monitors": [{ "kind": "mode_monitor", "name": "p1", "axis": 1, "position_m": 0.15e-6, "forward_sign": 1.0, "mode_index": 1, "freqs_hz": freqs_hz }]
        });
        let mut e_disp_fb = Engine::new();
        call(&mut e_disp_fb, "project_new", json!({ "title": "p4-fdtd-disp", "dim": 3, "h": 0.12e-6 }));
        call(
            &mut e_disp_fb,
            "add_box",
            json!({ "name": "b_wg", "min": [0.0, 0.0, 0.0], "max": [0.24e-6, 1.00e-6, 0.24e-6] }),
        );
        call(&mut e_disp_fb, "add_overlay", json!({ "regions": [{ "name": "core", "input": "b_wg" }] }));
        call(&mut e_disp_fb, "set_mesh", json!({ "h": 0.12e-6, "sweep": [{ "thickness": 1.00e-6, "n": 8 }] }));
        call(
            &mut e_disp_fb,
            "add_physics",
            json!({
                "physics": "ems",
                "study": "fdtd",
                "order": 1,
                "study_params": { "config_json": fdtd_disp_cfg.to_string(), "backend": "gpu" }
            }),
        );
        call(
            &mut e_disp_fb,
            "set_material",
            json!({ "physics": "ems", "region": "core", "material": { "eps_r": 2.25 } }),
        );
        call(&mut e_disp_fb, "solve_run", json!({}));
        let res_fb = e_disp_fb.result("ems").expect("fdtd dispersive fallback result").0;
        let diag_fb = res_fb.diagnostics.as_ref().expect("fdtd dispersive fallback diagnostics");
        assert_eq!(diag_fb.solver_kind, "yee_subpixel_fdtd");
        assert!(
            diag_fb.warnings.iter().any(|w| w.code == "gpu_fallback_to_cpu"),
            "expected gpu_fallback_to_cpu warning on dispersive FDTD, got: {:?}",
            diag_fb.warnings
        );
    }
}
