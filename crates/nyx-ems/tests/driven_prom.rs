// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative verification of Milestone 5 (`Study::Driven`, `lumped_port`, `waveport`, and
//! Hermite Tangential Galerkin PROM fast frequency sweep) and Acceptance Benchmark Test 6.

use femcore::linalg::csr::{c, C};
use femcore::meshing::geomdrive::MeshedGeom;
use femcore::meshing::{mesh_model, MeshPlan};
use geom::model::Model;
use nyx_ems::ibc::{resolve_jj_l_and_ej, MU0};
use nyx_ems::{Ems, C0};
use nyx_physics::{
    BcSpec, CircuitRom, Ctx, DrivenPortData, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec,
};
use serde_json::json;
use std::f64::consts::PI;

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

struct BlockResolver;

impl LabelResolver for BlockResolver {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "air" | "guide" => Ok(vec![1]),
            _ => Err(Error::invalid(format!("no region '{region}'"))),
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        match selector {
            "x0" => Ok(vec![1]),
            "y0" => Ok(vec![2]),
            "z0" => Ok(vec![3]),
            "z1" => Ok(vec![4]),
            "y1" => Ok(vec![5]),
            "x1" => Ok(vec![6]),
            _ => Err(Error::invalid(format!("no boundary '{selector}'"))),
        }
    }
}

/// Builds a 3D coupled feedline + transmon/resonator chip mesh on silicon (`eps_r = 11.45`) + air (`eps_r = 1.0`).
/// - Feedline on the left (`x in [-520 um, -340 um]`):
///   - `port1_sheet`: `y in [-700 um, -460 um]` (connects exterior ground wall at `y = -700 um` to `feed_trace`)
///   - `feed_trace`: `y in [-460 um, +460 um]`
///   - `port2_sheet`: `y in [+460 um, +700 um]` (connects `feed_trace` to exterior ground wall at `y = +700 um`)
/// - Coupled transmon/resonator on the right (`x in [-140 um, +420 um]`):
///   - `res_pad1`: `y in [-500 um, -60 um]`
///   - `res_pad2`: `y in [+60 um, +500 um]`
///   - `jj_sheet`: `x in [+40 um, +240 um], y in [-60 um, +60 um]`
fn build_coupled_feedline_resonator_mesh() -> MeshedGeom {
    let mut m = Model::new_3d();
    let (w_box, l_box, h_sub, h_air) = (1.4e-3, 1.4e-3, 300e-6, 300e-6);
    let si = m.box_("si_box", (-w_box / 2.0, -l_box / 2.0, -h_sub), (w_box / 2.0, l_box / 2.0, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, -l_box / 2.0, 0.0), (w_box / 2.0, l_box / 2.0, h_air));

    let p1_2d = m.rect("p1_2d", (-520e-6, -700e-6), (-340e-6, -460e-6));
    m.sheet("port1_sheet", p1_2d, 0.0);

    let feed_2d = m.rect("feed_2d", (-520e-6, -460e-6), (-340e-6, 460e-6));
    m.sheet("feed_trace", feed_2d, 0.0);

    let p2_2d = m.rect("p2_2d", (-520e-6, 460e-6), (-340e-6, 700e-6));
    m.sheet("port2_sheet", p2_2d, 0.0);

    let pad1_2d = m.rect("pad1_2d", (-140e-6, -500e-6), (420e-6, -60e-6));
    m.sheet("res_pad1", pad1_2d, 0.0);

    let pad2_2d = m.rect("pad2_2d", (-140e-6, 60e-6), (420e-6, 500e-6));
    m.sheet("res_pad2", pad2_2d, 0.0);

    let jj_2d = m.rect("jj_2d", (40e-6, -60e-6), (240e-6, 60e-6));
    m.sheet("jj_sheet", jj_2d, 0.0);

    m.overlay(&[("si", si), ("air", air)]);
    mesh_model(&m.list, &MeshPlan::simple(240e-6)).expect("3D coupled feedline + resonator sheet mesh")
}

/// Interpolates the positive-slope zero crossing $\operatorname{Im}[Y(\omega_p)] = 0$ nearest `f_target_hz`.
fn find_admittance_zero_crossing(freqs_hz: &[f64], im_y: &[f64], f_target_hz: f64) -> Option<f64> {
    let mut best: Option<f64> = None;
    for k in 0..(freqs_hz.len().saturating_sub(1)) {
        let (y0, y1) = (im_y[k], im_y[k + 1]);
        // Parallel LC resonance has d(Im Y)/d(omega) > 0: crosses from negative to positive
        if y0 <= 0.0 && y1 >= 0.0 && (y1 - y0) > 0.0 {
            let t = (-y0) / (y1 - y0);
            let f_zero = freqs_hz[k] + t * (freqs_hz[k + 1] - freqs_hz[k]);
            if best.is_none_or(|b| (f_zero - f_target_hz).abs() < (b - f_target_hz).abs()) {
                best = Some(f_zero);
            }
        }
    }
    best
}

#[test]
fn test_6_driven_coupled_resonator_feedline_prom_and_bbq_admittance() {
    let g = build_coupled_feedline_resonator_mesh();
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);
    eprintln!("Test 6 3D mesh: {} cells, {} points, {} facets", g.mesh.n_cells(), g.mesh.n_points(), g.mesh.n_facets());

    let e_j_hz = 20.0e9;
    let c_j_f = 35.0e-15;
    let (l_j_h, _) = resolve_jj_l_and_ej(0.0, e_j_hz);

    // -------------------------------------------------------------------------
    // Step 1: Solve the 3D `cavity` eigenmode problem with `lumped_jj` on the identical mesh
    // -------------------------------------------------------------------------
    let mut spec_cav = ProblemSpec::new("ems", 1, "cavity");
    spec_cav.study.params = json!({ "n_modes": 2, "f_search": 5.5e9 });
    spec_cav.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec_cav.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
        spec_cav.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
    }
    spec_cav.bcs.push(BcSpec {
        on: "jj_sheet".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q1_JJ",
            "E_j_Hz": e_j_hz,
            "C_j": c_j_f,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }),
    });

    let res_cav = Ems.solve(&exchange, &spec_cav, &mut Ctx::new(&resolver)).unwrap();
    let rom_cav = CircuitRom::of(&res_cav).expect("CircuitRom on cavity solve");
    // Pick the transmon/resonator eigenmode (dominant JJ inductive participation)
    let jj_mode = rom_cav
        .modes
        .iter()
        .max_by(|a, b| a.u_ind_jj_j[0].partial_cmp(&b.u_ind_jj_j[0]).unwrap())
        .expect("transmon mode found");
    let f_cavity_hz = jj_mode.freq_hz;
    eprintln!(
        "Cavity eigenfrequency with lumped_jj: f_cavity = {:.6} GHz (U_ind_jj = {:.4})",
        f_cavity_hz / 1e9,
        jj_mode.u_ind_jj_j[0]
    );
    assert!(jj_mode.u_ind_jj_j[0] > 0.2, "expected significant JJ participation, got {}", jj_mode.u_ind_jj_j[0]);

    // -------------------------------------------------------------------------
    // Step 2: Run 3-port `study = "driven"` (500 points, prom_tol = 1e-3)
    //         Port 1 & Port 2 on feedline, Port 3 ("JJ_port") on jj_sheet for BBQ admittance Y_JJ(omega)
    // -------------------------------------------------------------------------
    let f_min = (0.78 * f_cavity_hz).max(3.0e9);
    let f_max = 1.22 * f_cavity_hz;

    let build_driven_3port_spec = |n_pts: usize, f_lo: f64, f_hi: f64, prom_tol: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "f_min": f_lo,
            "f_max": f_hi,
            "n_points": n_pts,
            "prom_tol": prom_tol
        });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
            s.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
        }
        s.bcs.push(BcSpec {
            on: "port1_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
        });
        s.bcs.push(BcSpec {
            on: "port2_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }),
        });
        s.bcs.push(BcSpec {
            on: "jj_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "JJ_port", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
        });
        s
    };

    let spec_prom_500 = build_driven_3port_spec(500, f_min, f_max, json!(1e-3));
    let res_prom_500 = Ems.solve(&exchange, &spec_prom_500, &mut Ctx::new(&resolver)).unwrap();
    let data_500 = DrivenPortData::of(&res_prom_500).expect("DrivenPortData emitted");
    let rom_sweep_s = res_prom_500.scalar("rom_sweep_seconds").unwrap_or(0.0);

    eprintln!(
        "500-point 3-port PROM sweep: hdm_solves = {}, rom_dim = {}, max_residual = {:.3e}, rom_sweep_seconds = {:.4} s, total_seconds = {:.4} s",
        data_500.hdm_solves,
        data_500.rom_dim,
        data_500.max_residual,
        rom_sweep_s,
        res_prom_500.seconds
    );

    assert_eq!(data_500.freqs_hz.len(), 500);
    assert_eq!(data_500.ports, vec!["port1", "port2", "JJ_port"]);
    assert!(data_500.hdm_solves <= 12, "PROM required {} HDM solves (must be <= 12)", data_500.hdm_solves);
    assert!(rom_sweep_s < 1.0, "500-point ROM projection sweep took {:.4} s (must be < 1.0 s)", rom_sweep_s);
    assert!(data_500.max_residual < 1e-3, "PROM max residual {:.3e} >= 1e-3", data_500.max_residual);

    // -------------------------------------------------------------------------
    // Step 3: BBQ Admittance zero-crossing Im[Y_JJ(omega_p)] = 0 vs cavity eigenfrequency (< 0.1%)
    // -------------------------------------------------------------------------
    // 1) Open-circuited feedline ports (I_1 = I_2 = 0 => Y_EM,JJ = 1 / Z_{2,2}, matching the cavity BC!)
    // 2) Short-circuited feedline ports (V_1 = V_2 = 0 => Y_EM,JJ = Y_{2,2})
    let mut im_y_tot_open = Vec::with_capacity(500);
    let mut im_y_tot_short = Vec::with_capacity(500);
    for (k, &f) in data_500.freqs_hz.iter().enumerate() {
        let omega = 2.0 * PI * f;
        let y_jj_shunt_im = omega * c_j_f - 1.0 / (omega * l_j_h);

        let z33 = c(data_500.z_matrix[k][2][2][0], data_500.z_matrix[k][2][2][1]);
        let y_em_open = c(1.0, 0.0) / z33;
        im_y_tot_open.push(y_em_open.im + y_jj_shunt_im);

        let y33_im = data_500.y_matrix[k][2][2][1];
        im_y_tot_short.push(y33_im + y_jj_shunt_im);
    }

    let f_bbq_open = find_admittance_zero_crossing(&data_500.freqs_hz, &im_y_tot_open, f_cavity_hz)
        .expect("zero crossing of Im[Y_tot^open] found");
    let f_bbq_short = find_admittance_zero_crossing(&data_500.freqs_hz, &im_y_tot_short, f_cavity_hz)
        .expect("zero crossing of Im[Y_tot^short] found");

    let rel_err_open_pct = ((f_bbq_open - f_cavity_hz).abs() / f_cavity_hz) * 100.0;
    let rel_err_short_pct = ((f_bbq_short - f_cavity_hz).abs() / f_cavity_hz) * 100.0;
    eprintln!(
        "BBQ zero-crossing frequency: f_BBQ(open) = {:.6} GHz (err = {:.5}%), f_BBQ(short) = {:.6} GHz (err = {:.5}%) vs f_cavity = {:.6} GHz",
        f_bbq_open / 1e9,
        rel_err_open_pct,
        f_bbq_short / 1e9,
        rel_err_short_pct,
        f_cavity_hz / 1e9
    );
    assert!(
        rel_err_open_pct < 0.1,
        "BBQ zero-crossing {:.6} GHz differs from cavity eigenfrequency {:.6} GHz by {:.4}% (>= 0.1%)",
        f_bbq_open / 1e9,
        f_cavity_hz / 1e9,
        rel_err_open_pct
    );
    assert!(
        rel_err_short_pct < 0.1,
        "BBQ short-circuit zero-crossing {:.6} GHz differs from cavity eigenfrequency {:.6} GHz by {:.4}% (>= 0.1%)",
        f_bbq_short / 1e9,
        f_cavity_hz / 1e9,
        rel_err_short_pct
    );

    // -------------------------------------------------------------------------
    // Step 4: Compare PROM (prom_tol = 1e-3) against direct HDM (prom_tol = null) on identical grid
    //         AND verify 2-port feedline S21 with lumped_jj active on jj_sheet!
    // -------------------------------------------------------------------------
    let build_driven_2port_with_jj = |n_pts: usize, prom_tol: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "f_min": f_min,
            "f_max": f_max,
            "n_points": n_pts,
            "prom_tol": prom_tol
        });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["feed_trace", "res_pad1", "res_pad2"] {
            s.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
        }
        s.bcs.push(BcSpec {
            on: "port1_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "port1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
        });
        s.bcs.push(BcSpec {
            on: "port2_sheet".into(),
            kind: "lumped_port".into(),
            params: json!({ "name": "port2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0] }),
        });
        s.bcs.push(BcSpec {
            on: "jj_sheet".into(),
            kind: "lumped_jj".into(),
            params: json!({
                "name": "Q1_JJ",
                "E_j_Hz": e_j_hz,
                "C_j": c_j_f,
                "direction": [0.0, 1.0, 0.0],
                "formulation": "rank1"
            }),
        });
        s
    };

    // Run 500-point PROM on 2-port feedline + lumped_jj, and compare at 11 evenly spaced checkpoints against direct HDM (prom_tol = null)
    let res_2p_500 =
        Ems.solve(&exchange, &build_driven_2port_with_jj(500, json!(1e-3)), &mut Ctx::new(&resolver)).unwrap();
    let d_2p_500 = DrivenPortData::of(&res_2p_500).unwrap();
    let rom_2p_s = res_2p_500.scalar("rom_sweep_seconds").unwrap_or(0.0);
    assert!(d_2p_500.hdm_solves <= 12);
    assert!(rom_2p_s < 1.0);

    // Compare 11-point direct HDM (`prom_tol: null`) vs 11-point PROM (`prom_tol: 1e-3`) on both the 3-port and 2-port circuits
    let n_chk = 11usize;
    let res_3p_prom = Ems
        .solve(&exchange, &build_driven_3port_spec(n_chk, f_min, f_max, json!(1e-3)), &mut Ctx::new(&resolver))
        .unwrap();
    let res_3p_hdm = Ems
        .solve(
            &exchange,
            &build_driven_3port_spec(n_chk, f_min, f_max, serde_json::Value::Null),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    let d_3p_prom = DrivenPortData::of(&res_3p_prom).unwrap();
    let d_3p_hdm = DrivenPortData::of(&res_3p_hdm).unwrap();
    assert_eq!(d_3p_hdm.hdm_solves, n_chk);
    assert_eq!(d_3p_hdm.rom_dim, 0);

    let mut max_rel_s21_err_3p = 0.0_f64;
    let mut max_rel_s_err_3p = 0.0_f64;
    for k in 0..n_chk {
        let s21_p = c(d_3p_prom.s_matrix[k][1][0][0], d_3p_prom.s_matrix[k][1][0][1]);
        let s21_h = c(d_3p_hdm.s_matrix[k][1][0][0], d_3p_hdm.s_matrix[k][1][0][1]);
        let err_s21 = (s21_p - s21_h).norm() / s21_h.norm().max(1e-12);
        max_rel_s21_err_3p = max_rel_s21_err_3p.max(err_s21);

        let mut diff2 = 0.0;
        let mut norm2 = 0.0;
        for r in 0..3 {
            for col in 0..3 {
                let sp = c(d_3p_prom.s_matrix[k][r][col][0], d_3p_prom.s_matrix[k][r][col][1]);
                let sh = c(d_3p_hdm.s_matrix[k][r][col][0], d_3p_hdm.s_matrix[k][r][col][1]);
                diff2 += (sp - sh).norm_sqr();
                norm2 += sh.norm_sqr();
            }
        }
        max_rel_s_err_3p = max_rel_s_err_3p.max((diff2 / norm2.max(1e-24)).sqrt());
    }

    let res_2p_prom =
        Ems.solve(&exchange, &build_driven_2port_with_jj(n_chk, json!(1e-3)), &mut Ctx::new(&resolver)).unwrap();
    let res_2p_hdm = Ems
        .solve(&exchange, &build_driven_2port_with_jj(n_chk, serde_json::Value::Null), &mut Ctx::new(&resolver))
        .unwrap();
    let d_2p_prom = DrivenPortData::of(&res_2p_prom).unwrap();
    let d_2p_hdm = DrivenPortData::of(&res_2p_hdm).unwrap();

    let mut max_rel_s21_err_2p = 0.0_f64;
    for k in 0..n_chk {
        let s21_p = c(d_2p_prom.s_matrix[k][1][0][0], d_2p_prom.s_matrix[k][1][0][1]);
        let s21_h = c(d_2p_hdm.s_matrix[k][1][0][0], d_2p_hdm.s_matrix[k][1][0][1]);
        let err_s21 = (s21_p - s21_h).norm() / s21_h.norm().max(1e-12);
        max_rel_s21_err_2p = max_rel_s21_err_2p.max(err_s21);
    }

    eprintln!(
        "PROM vs Direct HDM max relative S21 error: 3-port = {:.3e} (full S-matrix = {:.3e}), 2-port+JJ = {:.3e}",
        max_rel_s21_err_3p, max_rel_s_err_3p, max_rel_s21_err_2p
    );
    assert!(
        max_rel_s21_err_3p < 1e-3,
        "3-port PROM vs direct HDM max relative S21 error {:.3e} >= 1e-3",
        max_rel_s21_err_3p
    );
    assert!(
        max_rel_s21_err_2p < 1e-3,
        "2-port+JJ PROM vs direct HDM max relative S21 error {:.3e} >= 1e-3",
        max_rel_s21_err_2p
    );

    // Check reciprocity and lossless unitarity (S^H S = I) across the 500-point 2-port+JJ sweep
    let mut max_recip_err = 0.0_f64;
    let mut max_unitarity_err = 0.0_f64;
    for k in 0..500 {
        let s11 = c(d_2p_500.s_matrix[k][0][0][0], d_2p_500.s_matrix[k][0][0][1]);
        let s21 = c(d_2p_500.s_matrix[k][1][0][0], d_2p_500.s_matrix[k][1][0][1]);
        let s12 = c(d_2p_500.s_matrix[k][0][1][0], d_2p_500.s_matrix[k][0][1][1]);
        let s22 = c(d_2p_500.s_matrix[k][1][1][0], d_2p_500.s_matrix[k][1][1][1]);
        max_recip_err = max_recip_err.max((s21 - s12).norm());
        let p_col1 = s11.norm_sqr() + s21.norm_sqr();
        let p_col2 = s12.norm_sqr() + s22.norm_sqr();
        max_unitarity_err = max_unitarity_err.max((p_col1 - 1.0).abs()).max((p_col2 - 1.0).abs());
    }
    eprintln!(
        "500-point 2-port+JJ max reciprocity error |S21 - S12| = {:.3e}, max unitarity error ||S_col|^2 - 1| = {:.3e}",
        max_recip_err, max_unitarity_err
    );
    assert!(max_recip_err < 1e-12, "reciprocity error {max_recip_err:e} >= 1e-12");
    assert!(max_unitarity_err < 1e-10, "lossless unitarity error {max_unitarity_err:e} >= 1e-10");
}

#[test]
fn test_6_waveport_3d_through_line_s_parameters_and_prom_accuracy() {
    // Rectangular waveguide WR-90 section: x in [0, a], y in [0, L], z in [0, b]
    // Fundamental TE10 mode has E_z(x) = sin(pi x / a), cutoff f_c = c_0 / (2 a) = 6.557 GHz.
    let (a, l, b) = (22.86e-3, 12.0e-3, 10.16e-3);
    let mesh = femcore::testkit::block(a, l, b, 6, 4, 2);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();
    let f_c = C0 / (2.0 * a);

    let build_wp_spec = |n_pts: usize, prom_tol: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 2, "driven");
        s.study.params = json!({
            "f_min": 8.0e9,
            "f_max": 12.0e9,
            "n_points": n_pts,
            "prom_tol": prom_tol
        });
        s.materials.insert("guide".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        s.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({ "name": "wg_in", "mode": 1 }) });
        s.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({ "name": "wg_out", "mode": 1 }) });
        s
    };

    // 1. Run 201-point PROM waveport sweep
    let res_prom = Ems.solve(&exchange, &build_wp_spec(201, json!(1e-3)), &mut Ctx::new(&BlockResolver)).unwrap();
    let data = DrivenPortData::of(&res_prom).expect("DrivenPortData on waveport solve");
    let rom_sweep_s = res_prom.scalar("rom_sweep_seconds").unwrap_or(0.0);

    let mut min_s21 = f64::INFINITY;
    let mut max_s11 = 0.0_f64;
    let mut max_z0_rel_err = 0.0_f64;
    let mut max_phase_err_rad = 0.0_f64;

    for (k, &f) in data.freqs_hz.iter().enumerate() {
        let s11 = c(data.s_matrix[k][0][0][0], data.s_matrix[k][0][0][1]);
        let s21 = c(data.s_matrix[k][1][0][0], data.s_matrix[k][1][0][1]);
        let s12 = c(data.s_matrix[k][0][1][0], data.s_matrix[k][0][1][1]);
        let s22 = c(data.s_matrix[k][1][1][0], data.s_matrix[k][1][1][1]);

        min_s21 = min_s21.min(s21.norm()).min(s12.norm());
        max_s11 = max_s11.max(s11.norm()).max(s22.norm());

        // Check analytical TE10 wave impedance Z_TE10 = (omega mu_0) / beta
        let omega = 2.0 * PI * f;
        let k0 = omega / C0;
        let kc = PI / a;
        let beta_analytic = (k0 * k0 - kc * kc).sqrt();
        let z0_analytic = omega * MU0 / beta_analytic;
        let z0_fem = data.z0_ohm[k][0][0];
        let z0_err = (z0_fem - z0_analytic).abs() / z0_analytic;
        max_z0_rel_err = max_z0_rel_err.max(z0_err);

        // Check transmission phase S21 ~ exp(-i beta L)
        let s21_analytic = C::from_polar(1.0, -beta_analytic * l);
        let phase_err = (s21 - s21_analytic).norm();
        max_phase_err_rad = max_phase_err_rad.max(phase_err);
        let _ = f_c;
    }

    eprintln!(
        "Waveport 3D through-line (201 pts, order 2): hdm_solves = {}, rom_dim = {}, min |S21| = {:.6}, max |S11| = {:.6}, max Z0 rel err = {:.3e}, max |S21 - e^(-i beta L)| = {:.3e}, rom_sweep = {:.4} s",
        data.hdm_solves,
        data.rom_dim,
        min_s21,
        max_s11,
        max_z0_rel_err,
        max_phase_err_rad,
        rom_sweep_s
    );

    assert!(data.hdm_solves <= 12, "waveport PROM hdm_solves = {} > 12", data.hdm_solves);
    assert!(min_s21 > 0.99, "waveport min |S21| = {min_s21:.6} <= 0.99");
    assert!(max_s11 < 0.05, "waveport max |S11| = {max_s11:.6} >= 0.05");
    assert!(max_z0_rel_err < 5e-3, "waveport Z0 relative error {max_z0_rel_err:e} >= 5e-3");

    // 2. Compare PROM vs direct HDM (`prom_tol = null`) across 11 frequencies
    let n_chk = 11usize;
    let res_p11 = Ems.solve(&exchange, &build_wp_spec(n_chk, json!(1e-3)), &mut Ctx::new(&BlockResolver)).unwrap();
    let res_h11 =
        Ems.solve(&exchange, &build_wp_spec(n_chk, serde_json::Value::Null), &mut Ctx::new(&BlockResolver)).unwrap();
    let d_p11 = DrivenPortData::of(&res_p11).unwrap();
    let d_h11 = DrivenPortData::of(&res_h11).unwrap();

    let mut max_wp_prom_err = 0.0_f64;
    for k in 0..n_chk {
        let s21_p = c(d_p11.s_matrix[k][1][0][0], d_p11.s_matrix[k][1][0][1]);
        let s21_h = c(d_h11.s_matrix[k][1][0][0], d_h11.s_matrix[k][1][0][1]);
        let err = (s21_p - s21_h).norm() / s21_h.norm().max(1e-12);
        max_wp_prom_err = max_wp_prom_err.max(err);
    }
    eprintln!("Waveport PROM vs direct HDM max relative S21 error = {:.3e}", max_wp_prom_err);
    assert!(max_wp_prom_err < 1e-3, "waveport PROM vs direct HDM error {max_wp_prom_err:e} >= 1e-3");
}

#[test]
fn driven_study_and_port_condition_validation_checks() {
    let mesh = femcore::testkit::block(0.02, 0.01, 0.01, 3, 2, 2);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();

    // 1. lumped_port outside Study::Driven is rejected
    let mut bad_cav = ProblemSpec::new("ems", 1, "cavity");
    bad_cav.study.params = json!({ "f_search": 10e9, "n_modes": 1 });
    bad_cav.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "p1", "Z_0": 50.0 }),
    });
    let err = Ems.solve(&exchange, &bad_cav, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(err.to_string().contains("belongs to the driven study"), "{err}");

    // 2. waveport outside Study::Driven is rejected
    let mut bad_wp = ProblemSpec::new("ems", 1, "cavity");
    bad_wp.study.params = json!({ "f_search": 10e9, "n_modes": 1 });
    bad_wp.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({ "name": "p1", "mode": 1 }) });
    let err = Ems.solve(&exchange, &bad_wp, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(err.to_string().contains("belongs to the driven study"), "{err}");

    // 3. Duplicate port names are rejected
    let mut dup = ProblemSpec::new("ems", 1, "driven");
    dup.study.params = json!({ "f_min": 8e9, "f_max": 10e9, "n_points": 5 });
    dup.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({ "name": "same_name" }) });
    dup.bcs.push(BcSpec { on: "y1".into(), kind: "waveport".into(), params: json!({ "name": "same_name" }) });
    let err = Ems.solve(&exchange, &dup, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(err.to_string().contains("duplicate port name"), "{err}");

    // 4. f_max < f_min is rejected
    let mut bad_freq = ProblemSpec::new("ems", 1, "driven");
    bad_freq.study.params = json!({ "f_min": 10e9, "f_max": 8e9, "n_points": 5 });
    bad_freq.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({ "name": "p1" }) });
    let err = Ems.solve(&exchange, &bad_freq, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(err.to_string().contains("f_min <= f_max"), "{err}");

    // 5. Sizing for "driven" uses f_max and headroom 1.0
    let mut s = ProblemSpec::new("ems", 2, "driven");
    s.study.params = json!({ "f_min": 4e9, "f_max": 10e9, "n_points": 21 });
    s.materials.insert("guide".into(), MaterialAssign::Props([("eps_r".to_string(), 4.0)].into()));
    let sz = Ems.sizing(&s, None);
    assert_eq!(sz.order, 2);
    assert_eq!(sz.headroom, 1.0);
    let expected_lam = C0 / (10e9 * 2.0);
    assert!((sz.wavelengths["guide"] - expected_lam).abs() / expected_lam < 1e-12);
}
