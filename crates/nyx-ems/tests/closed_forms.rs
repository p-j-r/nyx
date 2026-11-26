// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pillar (a) Verification: Exact Analytical Closed-Form Benchmarks (§4.2.3).
//!
//! Verifies all `nyx-ems` kernels (`modes`, `modes_complex`, `cavity`, `bloch`, `driven`, `ibc`,
//! `lumped_port`, `lumped_jj`) against independent closed-form electromagnetic solutions:
//! 1. WR-90 rectangular waveguide `TE10` cutoff and propagating/evanescent dispersion (`< 1%`)
//! 2. Circular waveguide `TE11` and `TM01` Bessel roots (`< 1%`)
//! 3. 3D rectangular and curved spherical PEC cavity eigenfrequencies (`< 1%`)
//! 4. Rectangular cavity IBC `Q` and rectangular/parallel-plate waveguide IBC attenuation `alpha` (`< 1%`)
//! 5. 1D periodic layered Bragg / Kronig-Penney analytic Bloch dispersion (`< 0.5%`)
//! 6. Lossy WR-90 driven waveguide `S21 = exp(-(alpha + j beta) L)` (`< 1%`)
//! 7. Lumped-port quarter-wave RLC cavity resonance `(f0, Q_L, S11(f0))` (`< 1%`)
//! 8. `lumped_jj` `distributed_sheet` vs `rank1` inductive frequency shift agreement (`< 1%`)

mod common;

use femcore::linalg::csr::c;
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::geomdrive3::sphere_cavity_with;
use femcore::meshing::{mesh_model, MeshPlan};
use femcore::optical::{EPS0, MU0};
use femcore::props::Periodic;
use geom::model::Model;
use nyx_ems::problem::{LumpedJjBc, LumpedPortBc, Problem, Study, SurfaceImpedanceBc, WavePortBc};
use nyx_ems::{Ems, C0};
use nyx_physics::{
    BcSpec, CircuitRom, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec,
    SolveResult,
};
use serde_json::json;
use std::f64::consts::PI;

struct GeomRes<'a>(&'a femcore::meshing::geomdrive::MeshedGeom);

impl LabelResolver for GeomRes<'_> {
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

fn solve_problem(mesh: &nyx_physics::Discretisation, p: &Problem) -> SolveResult {
    common::solve(mesh, p).expect("closed-form solve failed")
}

#[test]
fn wr90_te10_cutoff_and_dispersion_closed_form() {
    let (a, b) = (22.86e-3, 10.16e-3);
    let fc_exact = C0 / (2.0 * a);
    let mesh = common::rect(a, b, 18, 9);

    // 1. Just above cutoff: f = 1.05 * fc -> propagating TE10
    let f_above = 1.05 * fc_exact;
    let neff_above_exact = (1.0 - (fc_exact / f_above).powi(2)).sqrt();
    let mut p_ab = Problem::new(
        Study::Modes { freq: f_above, n_target: neff_above_exact * 1.002, evanescent: false, n_modes: 2 },
        2,
    );
    p_ab.materials.default.insert("eps_r".into(), 1.0);
    let r_ab = solve_problem(&mesh, &p_ab);
    let neff_ab = common::neff(&r_ab.solutions[0]).0;
    let fc_from_above = f_above * (1.0 - neff_ab * neff_ab).sqrt();
    let rel_ab = ((neff_ab - neff_above_exact) / neff_above_exact).abs();
    let rel_fc_ab = ((fc_from_above - fc_exact) / fc_exact).abs();
    assert!(rel_ab < 0.01 && rel_fc_ab < 0.01, "WR-90 above cutoff: neff rel {rel_ab:.3e}, fc rel {rel_fc_ab:.3e}");

    // 2. Just below cutoff: f = 0.95 * fc -> evanescent TE10
    let f_below = 0.95 * fc_exact;
    let neff_below_exact = ((fc_exact / f_below).powi(2) - 1.0).sqrt();
    let mut p_be = Problem::new(
        Study::Modes { freq: f_below, n_target: neff_below_exact * 1.002, evanescent: true, n_modes: 3 },
        2,
    );
    p_be.materials.default.insert("eps_r".into(), 1.0);
    let r_be = solve_problem(&mesh, &p_be);
    let best_be = r_be
        .solutions
        .iter()
        .map(|s| common::neff(s).1.abs())
        .min_by(|x, y| (x - neff_below_exact).abs().total_cmp(&(y - neff_below_exact).abs()))
        .unwrap();
    let fc_from_below = f_below * (1.0 + best_be * best_be).sqrt();
    let rel_be = ((best_be - neff_below_exact) / neff_below_exact).abs();
    let rel_fc_be = ((fc_from_below - fc_exact) / fc_exact).abs();
    assert!(rel_be < 0.01 && rel_fc_be < 0.01, "WR-90 below cutoff: imag rel {rel_be:.3e}, fc rel {rel_fc_be:.3e}");
}

#[test]
fn circular_waveguide_te11_tm01_bessel_roots() {
    let r = 10.0e-3;
    let mut m = Model::new_2d();
    let disk = m.circle("wg", (0.0, 0.0), r);
    m.overlay(&[("wg", disk)]);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(r / 14.0)).unwrap();

    let freq = 15.0e9;
    let k0 = 2.0 * PI * freq / C0;
    // Exact Bessel roots: J1'(p'_11) = 0 and J0(p_01) = 0
    let p_te11 = 1.8411837813406593;
    let p_tm01 = 2.404_825_557_695_773;
    let neff_te11 = (1.0 - (p_te11 / (k0 * r)).powi(2)).sqrt();
    let neff_tm01 = (1.0 - (p_tm01 / (k0 * r)).powi(2)).sqrt();

    let mut p = Problem::new(Study::Modes { freq, n_target: 0.90, evanescent: false, n_modes: 6 }, 2);
    p.materials.default.insert("eps_r".into(), 1.0);
    let res = solve_problem(&g.mesh, &p);
    let kc_rs: Vec<f64> = res
        .solutions
        .iter()
        .map(|s| {
            let ne = common::neff(s).0;
            r * k0 * (1.0 - ne * ne).sqrt()
        })
        .collect();

    let mut te11_matches = 0usize;
    let mut tm01_best = f64::INFINITY;
    for &kc_r in &kc_rs {
        if ((kc_r - p_te11) / p_te11).abs() < 0.01 {
            te11_matches += 1;
        }
        tm01_best = tm01_best.min(((kc_r - p_tm01) / p_tm01).abs());
    }
    assert!(
        te11_matches >= 2,
        "expected doubly-degenerate TE11 root {p_te11}, found kc*R = {kc_rs:?} (neff_te11={neff_te11:.4})"
    );
    assert!(
        tm01_best < 0.01,
        "expected TM01 root {p_tm01}, best rel err {tm01_best:.3e} (neff_tm01={neff_tm01:.4}) in {kc_rs:?}"
    );
}

#[test]
fn pec_rectangular_and_curved_spherical_cavity_eigenfrequencies() {
    // 1. 3D rectangular PEC box cavity
    let (a, b, d) = (0.04, 0.03, 0.02);
    let mut exact_rect = Vec::new();
    for m in 0..4usize {
        for n in 0..4usize {
            for p in 0..4usize {
                let nz = (m > 0) as usize + (n > 0) as usize + (p > 0) as usize;
                if nz < 2 {
                    continue;
                }
                let f = 0.5 * C0 * ((m as f64 / a).powi(2) + (n as f64 / b).powi(2) + (p as f64 / d).powi(2)).sqrt();
                exact_rect.push(f);
                if nz == 3 {
                    exact_rect.push(f);
                }
            }
        }
    }
    exact_rect.sort_by(|x, y| x.total_cmp(y));

    let box_mesh = common::block(a, b, d, 6, 5, 4);
    let mut pb = Problem::new(Study::Cavity { n_modes: 6, f_search: 0.9 * exact_rect[0] }, 2);
    pb.materials.default.insert("eps_r".into(), 1.0);
    pb.materials.default.insert("mu_r".into(), 1.0);
    let rb = solve_problem(&box_mesh, &pb);
    for (i, (sol, &ex)) in rb.solutions.iter().zip(&exact_rect).enumerate() {
        let got = sol.value.unwrap();
        let rel = ((got - ex) / ex).abs();
        assert!(rel < 0.01, "rect cavity mode {i}: got {got:.4e} vs {ex:.4e} (rel {rel:.3e})");
    }

    // 2. Curved 3D spherical PEC cavity: TM101 (kR = 2.743707...) and TE101 (kR = 4.493409...) triplets
    let r = 0.05;
    let mut opts = GeomMeshOpts::new(r / 4.2);
    opts.curved = true;
    let smesh = sphere_cavity_with(r, &opts).unwrap().mesh;
    assert!(smesh.has_curved());
    let kr_tm101 = 2.743707269992266;
    let kr_te101 = 4.493409457909064;
    let f_tm101 = C0 * kr_tm101 / (2.0 * PI * r);
    let f_te101 = C0 * kr_te101 / (2.0 * PI * r);

    let mut ps_tm = Problem::new(Study::Cavity { n_modes: 3, f_search: 0.95 * f_tm101 }, 2);
    ps_tm.materials.default.insert("eps_r".into(), 1.0);
    let r_tm = solve_problem(&smesh, &ps_tm);
    for (i, sol) in r_tm.solutions.iter().enumerate() {
        let rel = ((sol.value.unwrap() - f_tm101) / f_tm101).abs();
        assert!(rel < 0.01, "sphere TM101 mode {i}: rel {rel:.3e}");
    }

    let mut ps_te = Problem::new(Study::Cavity { n_modes: 3, f_search: 0.99 * f_te101 }, 2);
    ps_te.materials.default.insert("eps_r".into(), 1.0);
    let r_te = solve_problem(&smesh, &ps_te);
    for (i, sol) in r_te.solutions.iter().enumerate() {
        let rel = ((sol.value.unwrap() - f_te101) / f_te101).abs();
        assert!(rel < 0.01, "sphere TE101 mode {i}: rel {rel:.3e}");
    }
}

#[test]
fn ibc_cavity_q_and_waveguide_attenuation_closed_form() {
    let eta0 = (MU0 / EPS0).sqrt();
    let r_s = 0.05;

    // 1. Rectangular cavity TE110-type mode with resistive IBC (R_s = 0.05 Ohm) on z = 0, d (labels 3, 4)
    // E_z(x, y) = E0 sin(pi x / a) sin(pi y / b) => U = U_e + U_m = eps0 a b d E0^2 / 8, Q = k0 * eta0 * d / (2 * R_s).
    let (a, b, d): (f64, f64, f64) = (0.04, 0.03, 0.02);
    let k0_cav = PI * ((1.0 / a).powi(2) + (1.0 / b).powi(2)).sqrt();
    let f_cav = 0.5 * C0 * ((1.0 / a).powi(2) + (1.0 / b).powi(2)).sqrt();
    let q_exact = k0_cav * eta0 * d / (2.0 * r_s);
    let mesh3 = common::block(a, b, d, 6, 5, 3);
    for ibc_mode in ["direct", "perturbative"] {
        let mut p = Problem::new(Study::Cavity { n_modes: 1, f_search: 0.98 * f_cav }, 2);
        p.materials.default.insert("eps_r".into(), 1.0);
        p.surface_impedance.push(SurfaceImpedanceBc {
            labels: vec![3, 4],
            l_k_sq: 0.0,
            r_s,
            material: String::new(),
            thickness: 0.0,
            temperature: 0.0,
            f_ref: f_cav,
            film_model: "kerr".into(),
            lossless: false,
            ibc_mode: ibc_mode.into(),
        });
        let res = solve_problem(&mesh3, &p);
        let q_num = res.scalar("Q1").expect("IBC cavity Q1 recorded");
        let rel = ((q_num - q_exact) / q_exact).abs();
        assert!(rel < 0.01, "cavity IBC ({ibc_mode}) Q1 = {q_num:.2} vs exact {q_exact:.2} (rel {rel:.3e})");
    }

    // 2. WR-90 rectangular waveguide TE10 IBC attenuation alpha_c (Pozar Eq. 3.96)
    let (wa, wb) = (22.86e-3, 10.16e-3);
    let freq = 10.0e9;
    let k0 = 2.0 * PI * freq / C0;
    let fc = C0 / (2.0 * wa);
    let ratio = fc / freq;
    let alpha_wr90_exact = (r_s / (wb * eta0 * (1.0 - ratio * ratio).sqrt())) * (1.0 + 2.0 * (wb / wa) * ratio * ratio);
    let neff_wr90_re = (1.0 - ratio * ratio).sqrt();
    let mesh_wr = common::rect(wa, wb, 20, 10);
    let mut p_wr =
        Problem::new(Study::Modes { freq, n_target: neff_wr90_re * 1.001, evanescent: false, n_modes: 2 }, 2);
    p_wr.materials.default.insert("eps_r".into(), 1.0);
    p_wr.surface_impedance.push(SurfaceImpedanceBc {
        labels: vec![1, 2, 3, 4],
        l_k_sq: 0.0,
        r_s,
        material: String::new(),
        thickness: 0.0,
        temperature: 0.0,
        f_ref: freq,
        film_model: "kerr".into(),
        lossless: false,
        ibc_mode: "direct".into(),
    });
    let r_wr = solve_problem(&mesh_wr, &p_wr);
    let ne_wr = common::neff(&r_wr.solutions[0]);
    let alpha_wr_num = k0 * (-ne_wr.1);
    let rel_wr = ((alpha_wr_num - alpha_wr90_exact) / alpha_wr90_exact).abs();
    assert!(
        rel_wr < 0.01,
        "WR-90 IBC alpha = {alpha_wr_num:.6} Np/m vs exact {alpha_wr90_exact:.6} (rel {rel_wr:.3e})"
    );

    // 3. Parallel-plate TEM waveguide attenuation alpha_PP = R_s / (eta0 * h_gap)
    let (pp_w, pp_h) = (10.0e-3, 5.0e-3);
    let alpha_pp_exact = r_s / (eta0 * pp_h);
    let mesh_pp = common::rect(pp_w, pp_h, 8, 12);
    let mut p_pp = Problem::new(Study::Modes { freq, n_target: 1.001, evanescent: false, n_modes: 2 }, 2);
    p_pp.materials.default.insert("eps_r".into(), 1.0);
    p_pp.pmc = vec![1, 4];
    p_pp.surface_impedance.push(SurfaceImpedanceBc {
        labels: vec![2, 3],
        l_k_sq: 0.0,
        r_s,
        material: String::new(),
        thickness: 0.0,
        temperature: 0.0,
        f_ref: freq,
        film_model: "kerr".into(),
        lossless: false,
        ibc_mode: "direct".into(),
    });
    let r_pp = solve_problem(&mesh_pp, &p_pp);
    let ne_pp = common::neff(&r_pp.solutions[0]);
    let alpha_pp_num = k0 * (-ne_pp.1);
    let rel_pp = ((alpha_pp_num - alpha_pp_exact) / alpha_pp_exact).abs();
    assert!(
        rel_pp < 0.01,
        "parallel-plate TEM IBC alpha = {alpha_pp_num:.6} vs exact {alpha_pp_exact:.6} (rel {rel_pp:.3e})"
    );
}

#[test]
fn kronig_penney_1d_periodic_bragg_bloch_dispersion() {
    // 1D periodic bilayer stack along y: (n1, d1) = (3.4, 5 mm), (n2, d2) = (1.5, 8 mm), period L = 13 mm.
    let (n1, d1) = (3.4_f64, 5.0e-3);
    let (n2, d2) = (1.5_f64, 8.0e-3);
    let l = d1 + d2;
    let a = 4.0e-3;
    let mesh = common::stack_y(a, &[(d1, 10), (d2, 16)], a / 2.0);

    for freq in [2.0e9, 3.8e9] {
        let k0 = 2.0 * PI * freq / C0;
        let (p1, p2) = (n1 * k0 * d1, n2 * k0 * d2);
        let rhs = p1.cos() * p2.cos() - 0.5 * (n1 / n2 + n2 / n1) * p1.sin() * p2.sin();
        assert!(rhs.abs() < 0.99, "frequency {freq} Hz must be inside passband, cos(KL)={rhs}");
        let nb_exact = rhs.acos() / (k0 * l);

        let mut pb = Problem::new(Study::Bloch { freq, n_target: nb_exact * 1.005, n_modes: 2, brillouin_zone: 1 }, 2);
        pb.materials.set(1, "eps_r", n1 * n1);
        pb.materials.set(2, "eps_r", n2 * n2);
        pb.periodic = vec![
            Periodic { source: vec![1], dest: vec![6], phase: 0.0 },
            Periodic { source: vec![3], dest: vec![4], phase: 0.0 },
            Periodic { source: vec![2], dest: vec![5], phase: 0.0 },
        ];
        let res = solve_problem(&mesh, &pb);
        let nb_num = common::neff(&res.solutions[0]).0;
        let rel = ((nb_num - nb_exact) / nb_exact).abs();
        assert!(
            rel < 0.005,
            "Kronig-Penney at {freq:.2e} Hz: Bloch n_B = {nb_num:.6} vs exact {nb_exact:.6} (rel {rel:.3e})"
        );
    }
}

#[test]
fn lossy_wr90_driven_waveport_s21_attenuation_and_phase() {
    let (a, l, b) = (22.86e-3, 25.0e-3, 10.16e-3);
    let mesh = common::block(a, l, b, 10, 10, 2);
    let (eps_r, eps_r_im) = (2.2_f64, 0.11_f64);
    let freq = 10.0e9;

    let mut p = Problem::new(
        Study::Driven {
            f_min: freq,
            f_max: freq,
            n_points: 1,
            freqs_hz: Some(vec![freq]),
            prom_tol: None,
            max_hdm_solves: 10,
            strict_prom_tol: false,
            seed_freqs_hz: Vec::new(),
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p.materials.default.insert("eps_r".into(), eps_r);
    p.materials.default.insert("eps_r_im".into(), eps_r_im);
    p.waveports.push(WavePortBc { labels: vec![2], name: "p1".into(), mode: 1, n_eff_guess: None, z0_ohm: None });
    p.waveports.push(WavePortBc { labels: vec![5], name: "p2".into(), mode: 1, n_eff_guess: None, z0_ohm: None });

    let res = solve_problem(&mesh, &p);
    let data = DrivenPortData::of(&res).expect("DrivenPortData emitted");
    let s21 = c(data.s_matrix[0][1][0][0], data.s_matrix[0][1][0][1]);

    // Exact lossy WR-90 TE10 propagation constant: gamma = sqrt((pi/a)^2 - k0^2 * (eps_r - j eps_r_im))
    let k0 = 2.0 * PI * freq / C0;
    let eps_c = c(eps_r, -eps_r_im);
    let g2 = c((PI / a).powi(2), 0.0) - eps_c * (k0 * k0);
    let mut gamma = g2.sqrt();
    if gamma.re < 0.0 || gamma.im < 0.0 {
        gamma = -gamma;
    }
    let s21_exact = (-gamma * l).exp();
    let mag_rel = ((s21.norm() - s21_exact.norm()) / s21_exact.norm()).abs();
    let complex_rel = (s21 - s21_exact).norm() / s21_exact.norm();
    assert!(
        mag_rel < 0.01 && complex_rel < 0.01,
        "lossy WR-90 S21 = {s21} vs exact {s21_exact} (mag rel {mag_rel:.3e}, complex rel {complex_rel:.3e})"
    );
}

#[test]
fn lumped_port_rlc_cavity_resonance_f0_ql_and_s11() {
    // Parallel-plate TEM quarter-wave stub resonator [0, a] x [0, L] x [0, b]:
    // PMC at x=0,a (labels 1, 6), PEC at z=0,b (labels 3,4) and y=L (label 5), lumped_port at y=0 (label 2).
    // Uniform Ez(y) = E0 cos(pi y / (2 L)) => U = eps0 a b L E0^2 / 4,
    // exact f0 = c0/(4 L), Q0 = 1/tan_delta, Q_ext = pi a Z0 / (4 b eta0).
    let (a, l, b) = (10.0e-3, 12.0e-3, 5.0e-3);
    let mesh = common::block(a, l, b, 2, 12, 2);
    let tan_delta = 2.0e-3;
    let z0 = 120_000.0;
    let eta0 = (MU0 / EPS0).sqrt();

    let f0_exact = C0 / (4.0 * l);
    let q0_exact = 1.0 / tan_delta;
    let q_ext_exact = PI * a * z0 / (4.0 * b * eta0);
    let ql_exact = 1.0 / (1.0 / q0_exact + 1.0 / q_ext_exact);
    let s11_res_exact = ((q_ext_exact - q0_exact) / (q_ext_exact + q0_exact)).abs();

    // 1. Check Cavity eigenmode (f0, Q_L) with complex_damping lumped_port
    let mut pc = Problem::new(Study::Cavity { n_modes: 1, f_search: 0.95 * f0_exact }, 2);
    pc.materials.default.insert("eps_r".into(), 1.0);
    pc.materials.default.insert("eps_r_im".into(), tan_delta);
    pc.pmc = vec![1, 6];
    pc.lumped_ports.push(LumpedPortBc {
        labels: vec![2],
        name: "p1".into(),
        z0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: true,
        mode: "single".into(),
    });
    let rc = solve_problem(&mesh, &pc);
    let f0_cav = rc.solutions[0].value.unwrap();
    let ql_cav = rc.scalar("Q1").expect("RLC cavity Q1 recorded");
    assert!(((f0_cav - f0_exact) / f0_exact).abs() < 0.01, "RLC cavity f0 {f0_cav:.4e} vs {f0_exact:.4e}");
    assert!(((ql_cav - ql_exact) / ql_exact).abs() < 0.01, "RLC cavity Q_L {ql_cav:.2} vs {ql_exact:.2}");

    // 2. Driven sweep around resonance: extract f0 (at min |S11|), S11(f0), and Q_L from FWHM of (1 - |S11|^2)
    let span = 0.05 * f0_exact;
    let mut pd = Problem::new(
        Study::Driven {
            f_min: f0_exact - span,
            f_max: f0_exact + span,
            n_points: 401,
            freqs_hz: None,
            prom_tol: Some(1e-5),
            max_hdm_solves: 12,
            strict_prom_tol: true,
            seed_freqs_hz: vec![f0_exact],
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    pd.materials = pc.materials.clone();
    pd.pmc = vec![1, 6];
    pd.lumped_ports.push(LumpedPortBc {
        labels: vec![2],
        name: "p1".into(),
        z0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: false,
        mode: "single".into(),
    });
    let rd = solve_problem(&mesh, &pd);
    let data = DrivenPortData::of(&rd).expect("DrivenPortData emitted");
    let abs_loss: Vec<f64> = data.s_matrix.iter().map(|s| 1.0 - c(s[0][0][0], s[0][0][1]).norm_sqr()).collect();
    let (k_pk, &pk_val) = abs_loss.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap();
    let f0_drv = data.freqs_hz[k_pk];
    let s11_drv = c(data.s_matrix[k_pk][0][0][0], data.s_matrix[k_pk][0][0][1]).norm();
    let half = 0.5 * pk_val;
    let f_lo = (0..k_pk)
        .rev()
        .find_map(|k| {
            (abs_loss[k] <= half && abs_loss[k + 1] >= half).then(|| {
                let t = (half - abs_loss[k]) / (abs_loss[k + 1] - abs_loss[k]);
                data.freqs_hz[k] + t * (data.freqs_hz[k + 1] - data.freqs_hz[k])
            })
        })
        .expect("lower half-max");
    let f_hi = (k_pk..(abs_loss.len() - 1))
        .find_map(|k| {
            (abs_loss[k] >= half && abs_loss[k + 1] <= half).then(|| {
                let t = (abs_loss[k] - half) / (abs_loss[k] - abs_loss[k + 1]);
                data.freqs_hz[k] + t * (data.freqs_hz[k + 1] - data.freqs_hz[k])
            })
        })
        .expect("upper half-max");
    let ql_drv = f0_drv / (f_hi - f_lo);

    assert!(((f0_drv - f0_exact) / f0_exact).abs() < 0.01, "driven f0 {f0_drv:.4e} vs {f0_exact:.4e}");
    assert!(((ql_drv - ql_exact) / ql_exact).abs() < 0.01, "driven Q_L {ql_drv:.2} vs {ql_exact:.2}");
    assert!((s11_drv - s11_res_exact).abs() < 0.01, "driven |S11(f0)| = {s11_drv:.4} vs {s11_res_exact:.4}");
}

#[test]
fn lumped_jj_distributed_sheet_vs_rank1_frequency_shift_agreement() {
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.2e-3, 1.6e-3, 300e-6, 300e-6);
    let si = m.box_("si_box", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));
    let pad1_2d = m.rect("pad1_2d", (-160e-6, -500e-6), (160e-6, -40e-6));
    m.sheet("pad1", pad1_2d, 0.0);
    let pad2_2d = m.rect("pad2_2d", (-160e-6, 40e-6), (160e-6, 500e-6));
    m.sheet("pad2", pad2_2d, 0.0);
    let jj_2d = m.rect("jj_2d", (-30e-6, -40e-6), (30e-6, 40e-6));
    m.sheet("jj_sheet", jj_2d, 0.0);
    m.overlay(&[("si", si), ("air", air)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(220e-6)).unwrap();
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomRes(&g);

    let solve_with_jj = |e_j_hz: f64, f_search: f64, formulation: &str| -> f64 {
        let mut spec = ProblemSpec::new("ems", 2, "cavity");
        spec.study.params = json!({ "n_modes": 1, "f_search": f_search });
        spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".into(), 11.45)].into()));
        spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".into(), 1.0)].into()));
        spec.bcs.push(BcSpec { on: "pad1".into(), kind: "pec".into(), params: json!({}) });
        spec.bcs.push(BcSpec { on: "pad2".into(), kind: "pec".into(), params: json!({}) });
        spec.bcs.push(BcSpec {
            on: "jj_sheet".into(),
            kind: "lumped_jj".into(),
            params: json!({
                "name": "Q1_JJ",
                "E_j_Hz": e_j_hz,
                "C_j": 35.0e-15,
                "direction": [0.0, 1.0, 0.0],
                "formulation": formulation
            }),
        });
        let mut ctx = Ctx::new(&resolver);
        let res = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
        let rom = CircuitRom::of(&res).unwrap();
        rom.modes[0].freq_hz
    };

    // Compare inductive frequency and JJ tuning shift Delta f = f(18 GHz) - f(9 GHz) between rank1 and distributed_sheet
    let f_r_hi = solve_with_jj(18.0e9, 4.5e9, "rank1");
    let f_d_hi = solve_with_jj(18.0e9, 4.5e9, "distributed_sheet");
    let f_r_lo = solve_with_jj(9.0e9, 3.0e9, "rank1");
    let f_d_lo = solve_with_jj(9.0e9, 3.0e9, "distributed_sheet");
    let df_rank1 = (f_r_hi - f_r_lo).abs();
    let df_dist = (f_d_hi - f_d_lo).abs();
    let rel_f = ((f_d_hi - f_r_hi) / f_r_hi).abs();
    let rel_df = ((df_dist - df_rank1) / df_rank1).abs();
    assert!(
        rel_f < 0.01 && rel_df < 0.01,
        "lumped_jj distributed_sheet ({f_d_hi:.6e}, df={df_dist:.4e}) vs rank1 ({f_r_hi:.6e}, df={df_rank1:.4e}): rel_f={rel_f:.3e}, rel_df={rel_df:.3e}"
    );
    let _ = (LumpedJjBc {
        labels: vec![],
        name: String::new(),
        l_j: 0.0,
        e_j_hz: 0.0,
        c_j: 0.0,
        e_l_hz: 0.0,
        phi_ext_rad: 0.0,
        n_g: 0.0,
        squid_asymmetry: None,
        direction: [0.0, 1.0, 0.0],
        formulation: "rank1".into(),
    },);
}
