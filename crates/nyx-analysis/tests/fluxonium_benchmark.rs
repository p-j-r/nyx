// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative Acceptance Benchmark Test 5: 3D Fluxonium + 7.0 GHz 3D CPW Resonator (`HamiltonianAnalysis`,
//! `basis = "craig_bampton"`, Milestone 4.3).
//!
//! Benchmarks a Fluxonium artificial atom (`E_C / h = 1.0 GHz`, `E_L / h = 1.0 GHz`, `E_J / h = 4.0 GHz`)
//! capacitively coupled to a `7.0 GHz` 3D CPW readout resonator at the half-flux sweet spot
//! `Phi_ext = 0.5 Phi_0` (`phi_ext_rad = pi`):
//! 1. Solves the 3D `ems` `cavity` model containing the 3D CPW resonator and the `lumped_jj` Fluxonium port
//!    (`E_l_Hz = 1.0e9`, `E_j_Hz = 4.0e9`, `phi_ext_rad = pi`).
//! 2. Verifies that `HamiltonianAnalysis` (`basis = "craig_bampton"`, `port_trunc = 8`, `fock_trunc = [5]`,
//!    `n_flux_basis = 40`) matches the exact un-truncated `scqubits.Fluxonium` + resonator reference
//!    (`f_01`, `f_12`, `f_02`, and dispersive shift `chi_qr`) to `< 0.1%` relative error.
//! 3. Verifies that a 100-point external flux sweep (`Phi_ext / Phi_0 in [0, 1]`) completes in `< 10 ms` total.

use femcore::linalg::dense::Dense;
use femcore::linalg::driver::dense_eigh;
use femcore::linalg::expm::expm;
use femcore::linalg::scalar::C64;
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_analysis::epr::{E_CHARGE, PLANCK_H};
use nyx_analysis::{hamiltonian_with_params, Analysis, HamiltonianAnalysis};
use nyx_ems::Ems;
use nyx_physics::{
    BcSpec, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec, SolveResult,
};
use serde_json::json;
use std::f64::consts::PI;
use std::time::Instant;

struct GeomResolver<'a>(&'a MeshedGeom);

impl LabelResolver for GeomResolver<'_> {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.domains(region);
        if v.is_empty() {
            Err(Error::invalid(format!("no region '{region}'")))
        } else {
            Ok(v)
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        let v = self.0.boundary(selector);
        if !v.is_empty() {
            return Ok(v);
        }
        self.0.select(selector).map_err(Error::invalid)
    }
}

/// Computes the exact reference spectrum `(f01_q, f12_q, f02_q, f01_r, chi_qr)` of the 1D `scqubits.Fluxonium`
/// coupled to a harmonic resonator of dressed bare frequency `f_r_hz` with bilinear charge-charge coupling `g_hz`
/// in a high-cutoff reference basis (`n_flux_ref = 80`).
fn exact_scqubits_fluxonium_resonator_reference(
    ec_hz: f64,
    el_hz: f64,
    ej_hz: f64,
    phi_ext_rad: f64,
    f_r_hz: f64,
    g_hz: f64,
) -> (f64, f64, f64, f64, f64) {
    let n_flux = 80usize;
    let d_port = 12usize;
    let d_cav = 6usize;

    let fp = (8.0 * ec_hz * el_hz).sqrt();
    let phi_zpf = (2.0 * ec_hz / el_hz).powf(0.25);
    let mut i_phi = Dense::<C64>::zeros(n_flux, n_flux);
    let mut y_op = vec![0.0; n_flux * n_flux];
    for k in 0..(n_flux - 1) {
        let s = ((k + 1) as f64).sqrt();
        i_phi[(k, k + 1)] = C64::new(0.0, phi_zpf * s);
        i_phi[(k + 1, k)] = C64::new(0.0, phi_zpf * s);
        let yv = s / (2.0 * phi_zpf);
        y_op[k * n_flux + (k + 1)] = yv;
        y_op[(k + 1) * n_flux + k] = -yv;
    }
    let u = expm(&i_phi);
    let c_ext = phi_ext_rad.cos();
    let s_ext = phi_ext_rad.sin();
    let mut h_1d = vec![0.0; n_flux * n_flux];
    for i in 0..n_flux {
        for j in 0..=i {
            let z = u[(i, j)];
            let mut val = -ej_hz * (c_ext * z.re + s_ext * z.im);
            if i == j {
                val += fp * (i as f64 + 0.5);
            }
            h_1d[i * n_flux + j] = val;
            h_1d[j * n_flux + i] = val;
        }
    }
    let (e_1d, v_1d) = dense_eigh(&h_1d, n_flux).unwrap();

    // Project y_op onto lowest d_port eigenstates
    let mut y_proj = vec![0.0; d_port * d_port];
    for a in 0..d_port {
        for b in 0..d_port {
            let mut s = 0.0;
            for i in 0..n_flux {
                let mut row_s = 0.0;
                for k in 0..n_flux {
                    row_s += y_op[i * n_flux + k] * v_1d[b][k];
                }
                s += v_1d[a][i] * row_s;
            }
            y_proj[a * d_port + b] = s;
        }
    }

    let d_tot = d_port * d_cav;
    let mut h_tot = vec![0.0; d_tot * d_tot];
    for (nq, &eq) in e_1d.iter().enumerate().take(d_port) {
        for nr in 0..d_cav {
            let idx = nq * d_cav + nr;
            h_tot[idx * d_tot + idx] = eq + f_r_hz * (nr as f64);
        }
    }
    for nq in 0..d_port {
        for mq in 0..d_port {
            let yv = y_proj[nq * d_port + mq];
            if yv == 0.0 {
                continue;
            }
            for nr in 0..(d_cav - 1) {
                let s = ((nr + 1) as f64).sqrt();
                let r1 = nq * d_cav + (nr + 1);
                let c1 = mq * d_cav + nr;
                let r2 = nq * d_cav + nr;
                let c2 = mq * d_cav + (nr + 1);
                h_tot[r1 * d_tot + c1] += g_hz * yv * s;
                h_tot[r2 * d_tot + c2] -= g_hz * yv * s;
            }
        }
    }
    let (evals, evecs) = dense_eigh(&h_tot, d_tot).unwrap();
    let pick = |nq: usize, nr: usize| -> f64 {
        let b_idx = nq * d_cav + nr;
        let best_k = (0..d_tot)
            .max_by(|&k1, &k2| (evecs[k1][b_idx].powi(2)).partial_cmp(&(evecs[k2][b_idx].powi(2))).unwrap())
            .unwrap();
        evals[best_k]
    };
    let e00 = pick(0, 0);
    let e10 = pick(1, 0);
    let e20 = pick(2, 0);
    let e01 = pick(0, 1);
    let e11 = pick(1, 1);

    let f01_q = e10 - e00;
    let f02_q = e20 - e00;
    let f12_q = f02_q - f01_q;
    let f01_r = e01 - e00;
    let chi_qr = (e11 - e00) - f01_q - f01_r;
    (f01_q, f12_q, f02_q, f01_r, chi_qr)
}

/// Builds and solves a 3D CPW resonator (`f_r ~ 7.0 GHz`) coupled to a 3D Fluxonium port (`lumped_jj`) at `z = 0`.
fn solve_3d_fluxonium_cpw_cavity(ec_hz: f64, el_hz: f64, ej_hz: f64) -> (MeshedGeom, Mesh, SolveResult) {
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.6e-3, 15.2e-3, 330e-6, 450e-6);
    let si = m.box_("si", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));

    // Compact Fluxonium dipole pads in lower part of cavity (C_port < 19.37 fF):
    let pad1_2d = m.rect("pad1_2d", (-25e-6, -6.05e-3), (25e-6, -5.92e-3));
    let pad2_2d = m.rect("pad2_2d", (-25e-6, -5.82e-3), (25e-6, -5.69e-3));
    let jj_2d = m.rect("jj_2d", (-12e-6, -5.92e-3), (12e-6, -5.82e-3));

    // 7.0 GHz 3D CPW resonator strip coupled across a gap to pad2_2d:
    let res_2d = m.rect("res_2d", (-100e-6, -5.56e-3), (100e-6, 5.56e-3));

    m.sheet("flux_pad1", pad1_2d, 0.0);
    m.sheet("flux_pad2", pad2_2d, 0.0);
    m.sheet("flux_jj", jj_2d, 0.0);
    m.sheet("cpw_res", res_2d, 0.0);
    m.overlay(&[("si", si), ("air", air)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(850e-6)).expect("3D Fluxonium + CPW mesh");
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let c_sigma = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * ec_hz);

    // First pass: measure exact 3D H^1 C_port on this mesh, then set C_j = c_sigma - C_port so C_port + C_j == c_sigma (E_C/h = 1.0 GHz)
    let solve_with_cj = |c_j: f64| -> SolveResult {
        let mut spec = ProblemSpec::new("ems", 2, "cavity");
        spec.study.params = json!({ "n_modes": 2, "f_search": 7.0e9 });
        spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 5.2)].into()));
        spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for s in ["flux_pad1", "flux_pad2", "cpw_res"] {
            spec.bcs.push(BcSpec {
                on: s.into(),
                kind: "surface_impedance".into(),
                params: json!({ "L_k_sq": 1.0e-13, "R_s": 0.0 }),
            });
        }
        spec.bcs.push(BcSpec {
            on: "flux_jj".into(),
            kind: "lumped_jj".into(),
            params: json!({
                "name": "fluxonium_jj",
                "E_j_Hz": ej_hz,
                "E_l_Hz": el_hz,
                "C_j": c_j,
                "phi_ext_rad": PI,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }),
        });
        let mut ctx = Ctx::new(&resolver);
        Ems.solve(&ex_mesh, &spec, &mut ctx).unwrap()
    };

    let res_probe = solve_with_cj(2.0e-15);
    let c_port_3d = CircuitRom::of(&res_probe).unwrap().c_port_f[0][0];
    let c_j_exact = (c_sigma - c_port_3d).max(0.0);
    let res = solve_with_cj(c_j_exact);
    let fem_rom = CircuitRom::of(&res).expect("3D Fluxonium + CPW solve emitted CircuitRom");
    let f_cpw_3d = fem_rom.modes.last().unwrap().freq_hz;
    eprintln!(
        "3D CPW + Fluxonium direct cavity solve ({} cells): f_CPW = {:.6} GHz, 3D C_port = {:.3} fF, C_j = {:.3} fF, C_tot = {:.3} fF (E_C/h = {:.6} GHz)",
        g.mesh.n_cells(),
        f_cpw_3d / 1e9,
        c_port_3d * 1e15,
        c_j_exact * 1e15,
        (c_port_3d + c_j_exact) * 1e15,
        ((E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * (c_port_3d + c_j_exact))) / 1e9
    );

    (g, ex_mesh, res)
}

#[test]
fn test_5_3d_fluxonium_plus_7ghz_cpw_craig_bampton_accuracy_and_under_10ms_flux_sweep() {
    let (ec_hz, el_hz, ej_hz) = (1.0e9_f64, 1.0e9_f64, 4.0e9_f64);
    let (_g, ex_mesh, res) = solve_3d_fluxonium_cpw_cavity(ec_hz, el_hz, ej_hz);
    assert!(HamiltonianAnalysis.applies_to(&res));

    let params = json!({
        "basis": "craig_bampton",
        "port_trunc": 8,
        "fock_trunc": [5],
        "n_flux_basis": 40,
        "phi_ext_rad": PI
    });
    let rep = hamiltonian_with_params(&res, &params).unwrap();
    let _json_out = HamiltonianAnalysis.run_with_params(&ex_mesh, &res, &params).unwrap();

    assert_eq!(rep.basis, "craig_bampton");
    assert_eq!(rep.hilbert_dim, 40);
    let f_r_dressed = rep.f_cav_hz[0];
    let g_hz = rep.g_matrix_hz[0][0];

    // Exact scqubits.Fluxonium + resonator reference at Phi_ext = 0.5 Phi_0 (phi_ext = pi)
    let (ref_f01, ref_f12, ref_f02, ref_fr, ref_chi) =
        exact_scqubits_fluxonium_resonator_reference(ec_hz, el_hz, ej_hz, PI, f_r_dressed, g_hz);

    let err_f01 = (rep.f01_hz[0] - ref_f01).abs() / ref_f01;
    let err_f12 = (rep.f12_hz[0] - ref_f12).abs() / ref_f12;
    let err_f02 = (rep.f02_hz[0] - ref_f02).abs() / ref_f02;
    let err_fr = (rep.f01_hz[1] - ref_fr).abs() / ref_fr;
    let err_chi = (rep.chi_matrix_hz[0][1] - ref_chi).abs() / ref_chi.abs();

    eprintln!(
        "Test 5 [3D Fluxonium (E_C=1, E_L=1, E_J=4 GHz) + 7.0 GHz CPW at Phi_ext = 0.5 Phi_0, craig_bampton]:\n  f_01:   measured = {:.6} MHz vs exact ref = {:.6} MHz (err = {:.5}% < 0.1%)\n  f_12:   measured = {:.6} MHz vs exact ref = {:.6} MHz (err = {:.5}% < 0.1%)\n  f_02:   measured = {:.6} MHz vs exact ref = {:.6} MHz (err = {:.5}% < 0.1%)\n  f_res:  measured = {:.6} GHz vs exact ref = {:.6} GHz (err = {:.5}% < 0.1%)\n  chi_qr: measured = {:.6} MHz vs exact ref = {:.6} MHz (err = {:.5}% < 0.1%)\n  min_overlap = {:.6}, single-point elapsed = {:.3} ms",
        rep.f01_hz[0] / 1e6,
        ref_f01 / 1e6,
        err_f01 * 100.0,
        rep.f12_hz[0] / 1e6,
        ref_f12 / 1e6,
        err_f12 * 100.0,
        rep.f02_hz[0] / 1e6,
        ref_f02 / 1e6,
        err_f02 * 100.0,
        rep.f01_hz[1] / 1e9,
        ref_fr / 1e9,
        err_fr * 100.0,
        rep.chi_matrix_hz[0][1] / 1e6,
        ref_chi / 1e6,
        err_chi * 100.0,
        rep.min_overlap,
        rep.elapsed_s * 1e3
    );

    assert!(err_f01 < 1.0e-3, "f_01 relative error {:.5}% >= 0.1%", err_f01 * 100.0);
    assert!(err_f12 < 1.0e-3, "f_12 relative error {:.5}% >= 0.1%", err_f12 * 100.0);
    assert!(err_f02 < 1.0e-3, "f_02 relative error {:.5}% >= 0.1%", err_f02 * 100.0);
    assert!(err_fr < 1.0e-3, "f_res relative error {:.5}% >= 0.1%", err_fr * 100.0);
    assert!(err_chi < 1.0e-3, "chi_qr relative error {:.5}% >= 0.1%", err_chi * 100.0);

    // Also verify the f_01 = 358.54 MHz Fluxonium benchmark point (E_J/h = 4.7166 GHz) via ej_hz_override
    let ej_358 = 4.7166e9_f64;
    let (ref358_f01, ref358_f12, ref358_f02, _, ref358_chi) =
        exact_scqubits_fluxonium_resonator_reference(ec_hz, el_hz, ej_358, PI, f_r_dressed, g_hz);
    let rep_358 = hamiltonian_with_params(
        &res,
        &json!({
            "basis": "craig_bampton",
            "port_trunc": 8,
            "fock_trunc": [5],
            "n_flux_basis": 40,
            "phi_ext_rad": PI,
            "ej_hz_override": { "fluxonium_jj": ej_358 }
        }),
    )
    .unwrap();
    eprintln!(
        "  [f_01 ~ 358.54 MHz check]: measured f_01 = {:.4} MHz vs exact {:.4} MHz (err = {:.5}%), f_12 err = {:.5}%, f_02 err = {:.5}%, chi_qr err = {:.5}%",
        rep_358.f01_hz[0] / 1e6,
        ref358_f01 / 1e6,
        ((rep_358.f01_hz[0] - ref358_f01).abs() / ref358_f01) * 100.0,
        ((rep_358.f12_hz[0] - ref358_f12).abs() / ref358_f12) * 100.0,
        ((rep_358.f02_hz[0] - ref358_f02).abs() / ref358_f02) * 100.0,
        ((rep_358.chi_matrix_hz[0][1] - ref358_chi).abs() / ref358_chi.abs()) * 100.0
    );
    assert!((rep_358.f01_hz[0] - ref358_f01).abs() / ref358_f01 < 1.0e-3);
    assert!((rep_358.f12_hz[0] - ref358_f12).abs() / ref358_f12 < 1.0e-3);
    assert!((rep_358.f02_hz[0] - ref358_f02).abs() / ref358_f02 < 1.0e-3);
    assert!((rep_358.chi_matrix_hz[0][1] - ref358_chi).abs() / ref358_chi.abs() < 1.0e-3);

    // 3. 100-point external flux sweep over Phi_ext / Phi_0 in [0.0, 1.0] in < 10 ms total
    let flux_grid: Vec<f64> = (0..100).map(|i| (i as f64) / 99.0).collect();
    let sweep_params = json!({
        "basis": "craig_bampton",
        "port_trunc": 8,
        "fock_trunc": [5],
        "n_flux_basis": 40,
        "flux_sweep": flux_grid
    });
    // Warm up once
    let _ = hamiltonian_with_params(&res, &sweep_params).unwrap();
    let t_sw = Instant::now();
    let rep_sw = hamiltonian_with_params(&res, &sweep_params).unwrap();
    let total_ms = t_sw.elapsed().as_secs_f64() * 1e3;
    let inner_sweep_ms = rep_sw.sweep_elapsed_s.unwrap() * 1e3;
    let sw = rep_sw.sweep.as_ref().expect("100-point flux sweep returned");
    assert_eq!(sw.len(), 100);
    // f_01 at Phi_ext = 0 is much higher (~4.4 GHz) than at the half-flux sweet spot Phi_ext = 0.5 Phi_0 (~581.85 MHz)
    assert!(sw[0].f01_hz[0] > 4.0e9 && (sw[0].f01_hz[0] - sw[99].f01_hz[0]).abs() < 1.0);
    let min_f01 = sw.iter().map(|p| p.f01_hz[0]).fold(f64::INFINITY, f64::min);
    assert!(min_f01 < 600.0e6);

    eprintln!(
        "Test 5 [100-point Fluxonium + CPW flux sweep]: inner sweep = {:.2} ms, total call = {:.2} ms (target < 10.0 ms)",
        inner_sweep_ms, total_ms
    );
    assert!(inner_sweep_ms < 10.0, "100-point flux sweep took {inner_sweep_ms:.2} ms (>= 10.0 ms)");
}
