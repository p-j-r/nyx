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
    assert!(max_phase_err_rad < 5e-3, "waveport max |S21 - exp(-i beta L)| = {max_phase_err_rad:e} >= 5e-3");

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

    // 1. lumped_port outside Study::Cavity | Study::Driven (e.g. 2D modes) is rejected, and succeeds in cavity
    let mesh2 = femcore::testkit::rect_split_x(0.02, 0.01, 0.01, 0.005);
    let exchange2 = Mesh::from_solver_mesh(&mesh2).unwrap();
    let mut bad_modes = ProblemSpec::new("ems", 1, "modes");
    bad_modes.study.params = json!({ "freq": 10e9, "n_target": 1.0, "n_modes": 1 });
    bad_modes.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "p1", "Z_0": 50.0 }),
    });
    let err = Ems.solve(&exchange2, &bad_modes, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(err.to_string().contains("belongs to the cavity or driven study"), "{err}");

    let mut ok_cav = ProblemSpec::new("ems", 1, "cavity");
    ok_cav.study.params = json!({ "f_search": 10e9, "n_modes": 1 });
    ok_cav.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "p1", "Z_0": 50.0 }),
    });
    let cav_res = Ems.solve(&exchange, &ok_cav, &mut Ctx::new(&BlockResolver)).unwrap();
    let cav_rom = CircuitRom::of(&cav_res).expect("CircuitRom emitted for cavity with lumped_port");
    assert_eq!(cav_rom.ports.len(), 1);
    assert_eq!(cav_rom.ports[0].kind, "lumped_port");
    assert_eq!(cav_rom.ports[0].z0_ohm, Some(50.0));

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

    // 5. Sizing for "driven" uses f_max and headroom 1.0 (and derives f_max from freqs_hz when given)
    let mut s = ProblemSpec::new("ems", 2, "driven");
    s.study.params = json!({ "f_min": 4e9, "f_max": 10e9, "n_points": 21 });
    s.materials.insert("guide".into(), MaterialAssign::Props([("eps_r".to_string(), 4.0)].into()));
    let sz = Ems.sizing(&s, None);
    assert_eq!(sz.order, 2);
    assert_eq!(sz.headroom, 1.0);
    let expected_lam = C0 / (10e9 * 2.0);
    assert!((sz.wavelengths["guide"] - expected_lam).abs() / expected_lam < 1e-12);

    let mut s_list = ProblemSpec::new("ems", 2, "driven");
    s_list.study.params = json!({ "freqs_hz": [4e9, 6e9, 12e9] });
    s_list.materials.insert("guide".into(), MaterialAssign::Props([("eps_r".to_string(), 4.0)].into()));
    let sz_list = Ems.sizing(&s_list, None);
    let expected_lam_list = C0 / (12e9 * 2.0);
    assert!((sz_list.wavelengths["guide"] - expected_lam_list).abs() / expected_lam_list < 1e-12);

    // 6. Non-ascending freqs_hz is rejected with Error::Invalid
    let mut bad_list = ProblemSpec::new("ems", 1, "driven");
    bad_list.study.params = json!({ "freqs_hz": [6e9, 5e9, 7e9] });
    bad_list.bcs.push(BcSpec { on: "y0".into(), kind: "waveport".into(), params: json!({ "name": "p1" }) });
    let err = Ems.solve(&exchange, &bad_list, &mut Ctx::new(&BlockResolver)).unwrap_err();
    assert!(matches!(err, Error::Invalid(_)), "expected Error::Invalid, got {err:?}");
}

#[test]
fn fr02_and_fr20_prom_budget_strict_tol_zoom_symbolic_lu_and_freqs_hz() {
    let g = build_coupled_feedline_resonator_mesh();
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let e_j_hz = 20.0e9;
    let c_j_f = 35.0e-15;

    // Reference cavity frequency
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
    let rom_cav = CircuitRom::of(&res_cav).unwrap();
    let f_cav =
        rom_cav.modes.iter().max_by(|a, b| a.u_ind_jj_j[0].partial_cmp(&b.u_ind_jj_j[0]).unwrap()).unwrap().freq_hz;

    let make_3port = |study_params: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = study_params;
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

    // 1. 3-port PROM with max_hdm_solves = 24, strict_prom_tol = true, seed_freqs_hz = [f_cav]
    let res_seeded = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "f_min": 0.78 * f_cav,
                "f_max": 1.22 * f_cav,
                "n_points": 101,
                "prom_tol": 1e-4,
                "max_hdm_solves": 24,
                "strict_prom_tol": true,
                "seed_freqs_hz": [f_cav]
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    let d_seeded = DrivenPortData::of(&res_seeded).unwrap();
    assert!(d_seeded.max_residual <= 1e-4, "residual {}", d_seeded.max_residual);
    assert!(d_seeded.hdm_solves <= 24);
    assert_eq!(res_seeded.scalar("symbolic_count"), Some(1.0), "symbolic LU must be reused across all anchors");
    assert!(
        d_seeded.anchor_freqs_hz.iter().any(|&f| (f - f_cav).abs() < 1e-3 * f_cav),
        "seed frequency should appear in anchors: {:?}",
        d_seeded.anchor_freqs_hz
    );

    // 2. strict_prom_tol = true fails cleanly when budget is exhausted before reaching tol
    let err_strict = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "f_min": 0.78 * f_cav,
                "f_max": 1.22 * f_cav,
                "n_points": 51,
                "prom_tol": 1e-12,
                "max_hdm_solves": 1,
                "strict_prom_tol": true
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap_err();
    let msg = err_strict.to_string();
    assert!(msg.contains("max_prom_residual") && msg.contains("prom_tol"), "unexpected error: {msg}");

    // 3. 2 MHz narrowband zoom sweep (verifies unit-L2 y1/y2 normalization + symbolic LU reuse on both PROM & HDM)
    let f_zoom_lo = f_cav - 1.0e6;
    let f_zoom_hi = f_cav + 1.0e6;
    let res_zoom_prom = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "f_min": f_zoom_lo,
                "f_max": f_zoom_hi,
                "n_points": 5,
                "prom_tol": 1e-5,
                "strict_prom_tol": true
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    let res_zoom_hdm = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "f_min": f_zoom_lo,
                "f_max": f_zoom_hi,
                "n_points": 5,
                "prom_tol": null
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    assert_eq!(res_zoom_prom.scalar("symbolic_count"), Some(1.0));
    assert_eq!(res_zoom_hdm.scalar("symbolic_count"), Some(1.0), "direct HDM loop must reuse symbolic LU");
    let dz_p = DrivenPortData::of(&res_zoom_prom).unwrap();
    let dz_h = DrivenPortData::of(&res_zoom_hdm).unwrap();
    for k in 0..5 {
        let s21_p = c(dz_p.s_matrix[k][1][0][0], dz_p.s_matrix[k][1][0][1]);
        let s21_h = c(dz_h.s_matrix[k][1][0][0], dz_h.s_matrix[k][1][0][1]);
        let err = (s21_p - s21_h).norm() / s21_h.norm().max(1e-12);
        assert!(err < 1e-4, "2 MHz zoom PROM vs HDM error at k={k}: {err:e}");
    }

    // 4. Explicit non-uniform freqs_hz without f_min/f_max
    let custom_freqs =
        vec![0.80 * f_cav, 0.95 * f_cav, 0.995 * f_cav, f_cav, 1.005 * f_cav, 1.05 * f_cav, 1.20 * f_cav];
    let res_custom_prom = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "freqs_hz": custom_freqs,
                "prom_tol": 1e-4,
                "strict_prom_tol": true
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    let res_custom_hdm = Ems
        .solve(
            &exchange,
            &make_3port(json!({
                "freqs_hz": custom_freqs,
                "prom_tol": null
            })),
            &mut Ctx::new(&resolver),
        )
        .unwrap();
    let dc_p = DrivenPortData::of(&res_custom_prom).unwrap();
    let dc_h = DrivenPortData::of(&res_custom_hdm).unwrap();
    assert_eq!(dc_p.freqs_hz, custom_freqs);
    assert_eq!(res_custom_prom.scalar("symbolic_count"), Some(1.0));
    assert_eq!(res_custom_hdm.scalar("symbolic_count"), Some(1.0));
    for k in 0..custom_freqs.len() {
        let s21_p = c(dc_p.s_matrix[k][1][0][0], dc_p.s_matrix[k][1][0][1]);
        let s21_h = c(dc_h.s_matrix[k][1][0][0], dc_h.s_matrix[k][1][0][1]);
        let err = (s21_p - s21_h).norm() / s21_h.norm().max(1e-12);
        assert!(err < 1e-3, "non-uniform freqs_hz PROM vs HDM error at k={k}: {err:e}");
    }
}

#[test]
fn fr13_driven_dielectric_loss_tangent_matches_qi_10000() {
    // 3D PEC box cavity (20 mm x 20 mm x 10 mm) with uniform inline loss_tangent = 1e-4 -> exact Q_i = 10000
    let (a, l, b) = (20e-3, 20e-3, 10e-3);
    let mesh = femcore::testkit::block(a, l, b, 4, 4, 2);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();
    let f_analytic = (C0 / 2.0) * (1.0 / (a * a) + 1.0 / (l * l)).sqrt();

    // Cavity eigenmode with PMC on y0 (matching high-Z_0 lumped_port open boundary) and loss_tangent = 1e-4 -> exact Q_i = 10000
    let f_pmc_analytic = (C0 / 2.0) * (1.0 / (a * a) + 1.0 / (4.0 * l * l)).sqrt();
    let mut spec_cav = ProblemSpec::new("ems", 2, "cavity");
    spec_cav.study.params = json!({ "n_modes": 1, "f_search": f_pmc_analytic });
    spec_cav.materials.insert(
        "air".into(),
        MaterialAssign::Props([("eps_r".to_string(), 1.0), ("loss_tangent".to_string(), 1.0e-4)].into()),
    );
    spec_cav.bcs.push(BcSpec { on: "y0".into(), kind: "pmc".into(), params: json!({}) });
    let res_cav = Ems.solve(&exchange, &spec_cav, &mut Ctx::new(&BlockResolver)).unwrap();
    let f_cav = res_cav.solutions[0].value.unwrap();
    let q_cav = res_cav.scalar("Q1").unwrap();
    assert!((q_cav - 10_000.0).abs() / 10_000.0 < 1e-3, "cavity Q_i = {q_cav:.2} != 10000");
    let _ = f_analytic;

    // Weakly-coupled driven sweep (Z_0 = 1e9 Ohm on y0 so Q_ext >> Q_i) around f_cav
    let df_span = 0.001 * f_cav; // +/- 0.1% spans ~20 linewidths for Q = 10000
    let mut spec_drv = ProblemSpec::new("ems", 2, "driven");
    spec_drv.study.params = json!({
        "f_min": f_cav - df_span,
        "f_max": f_cav + df_span,
        "n_points": 401,
        "prom_tol": 1e-5,
        "strict_prom_tol": true
    });
    spec_drv.materials.insert(
        "air".into(),
        MaterialAssign::Props([("eps_r".to_string(), 1.0), ("loss_tangent".to_string(), 1.0e-4)].into()),
    );
    spec_drv.bcs.push(BcSpec {
        on: "y0".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "probe", "Z_0": 1.0e9, "direction": [0.0, 0.0, 1.0] }),
    });
    let res_drv = Ems.solve(&exchange, &spec_drv, &mut Ctx::new(&BlockResolver)).unwrap();
    let data = DrivenPortData::of(&res_drv).unwrap();

    // Find peak of Re[Z11(f)] and extract Q_i = f_peak / FWHM (half-maximum of Lorentzian Re[Z11])
    let re_z: Vec<f64> = data.z_matrix.iter().map(|z| z[0][0][0]).collect();
    let (k_pk, &z_max) = re_z.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).expect("non-empty sweep");
    assert!(z_max > 0.0, "Re[Z11] must be positive for passive lossy cavity");
    let half = 0.5 * z_max;
    let mut f_lo_hm = None;
    for k in (0..k_pk).rev() {
        if re_z[k] <= half && re_z[k + 1] >= half {
            let t = (half - re_z[k]) / (re_z[k + 1] - re_z[k]);
            f_lo_hm = Some(data.freqs_hz[k] + t * (data.freqs_hz[k + 1] - data.freqs_hz[k]));
            break;
        }
    }
    let mut f_hi_hm = None;
    for k in k_pk..(re_z.len() - 1) {
        if re_z[k] >= half && re_z[k + 1] <= half {
            let t = (re_z[k] - half) / (re_z[k] - re_z[k + 1]);
            f_hi_hm = Some(data.freqs_hz[k] + t * (data.freqs_hz[k + 1] - data.freqs_hz[k]));
            break;
        }
    }
    let f_lo_hm = f_lo_hm.expect("lower half-max found");
    let f_hi_hm = f_hi_hm.expect("upper half-max found");
    let q_driven = data.freqs_hz[k_pk] / (f_hi_hm - f_lo_hm);
    eprintln!(
        "FR-13 driven loss_tangent=1e-4: f_cav = {:.6} GHz, Q_cav = {:.1}, Q_driven = {:.1}",
        f_cav / 1e9,
        q_cav,
        q_driven
    );
    assert!((q_driven - 10_000.0).abs() / 10_000.0 < 0.02, "driven Q_i = {q_driven:.1} differs from 10000 by > 2%");
    // Passive: |S11| < 1 at resonance
    let s11_pk = c(data.s_matrix[k_pk][0][0][0], data.s_matrix[k_pk][0][0][1]).norm();
    assert!(s11_pk < 1.0 - 1e-6, "lossy cavity |S11| at resonance must be < 1, got {s11_pk}");
}

#[test]
fn fr14_cpw_waveport_n_eff_guess_and_z0_ohm() {
    // 3D CPW through-line on Si (eps_r = 11.45) + air (eps_r = 1.0): quasi-TEM n_eff ~ 2.5
    let mut m = Model::new_3d();
    let (w_box, l_line, h_sub, h_air) = (600e-6, 1.2e-3, 250e-6, 250e-6);
    let (w_trace, s_gap) = (80e-6, 50e-6);
    let si = m.box_("si_box", (-w_box / 2.0, 0.0, -h_sub), (w_box / 2.0, l_line, 0.0));
    let air = m.box_("air_box", (-w_box / 2.0, 0.0, 0.0), (w_box / 2.0, l_line, h_air));
    let tr = m.rect("tr", (-w_trace / 2.0, 0.0), (w_trace / 2.0, l_line));
    m.sheet("cpw_trace", tr, 0.0);
    let gl = m.rect("gl", (-w_box / 2.0, 0.0), (-w_trace / 2.0 - s_gap, l_line));
    m.sheet("cpw_gnd_l", gl, 0.0);
    let gr = m.rect("gr", (w_trace / 2.0 + s_gap, 0.0), (w_box / 2.0, l_line));
    m.sheet("cpw_gnd_r", gr, 0.0);
    m.overlay(&[("si", si), ("air", air)]);
    let g = mesh_model(&m.list, &MeshPlan::simple(180e-6)).expect("CPW waveport mesh");
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let mut spec = ProblemSpec::new("ems", 1, "driven");
    spec.study.params = json!({ "f_min": 5.0e9, "f_max": 7.0e9, "n_points": 21, "prom_tol": 1e-3 });
    spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
        spec.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
    }
    spec.bcs.push(BcSpec {
        on: "plane(y = 0)".into(),
        kind: "waveport".into(),
        params: json!({ "name": "wp_in", "mode": 1, "n_eff_guess": 2.5, "z0_ohm": 50.0 }),
    });
    spec.bcs.push(BcSpec {
        on: format!("plane(y = {l_line:e})"),
        kind: "waveport".into(),
        params: json!({ "name": "wp_out", "mode": 1, "n_eff_guess": 2.5, "z0_ohm": 50.0 }),
    });

    let res = Ems.solve(&exchange, &spec, &mut Ctx::new(&resolver)).unwrap();
    let data = DrivenPortData::of(&res).unwrap();
    for k in 0..data.freqs_hz.len() {
        assert!((data.z0_ohm[k][0][0] - 50.0).abs() < 1e-12);
        assert!((data.z0_ohm[k][1][0] - 50.0).abs() < 1e-12);
        let s11 = c(data.s_matrix[k][0][0][0], data.s_matrix[k][0][0][1]);
        let s21 = c(data.s_matrix[k][1][0][0], data.s_matrix[k][1][0][1]);
        let s12 = c(data.s_matrix[k][0][1][0], data.s_matrix[k][0][1][1]);
        let s22 = c(data.s_matrix[k][1][1][0], data.s_matrix[k][1][1][1]);
        assert!(s21.norm() > 0.98, "CPW waveport |S21| = {} at k={k}", s21.norm());
        let p_tot = s11.norm_sqr() + s21.norm_sqr();
        assert!(p_tot <= 1.0 + 1e-8 && p_tot > 0.98, "waveport passivity/match |S11|^2+|S21|^2 = {p_tot} at k={k}");
        assert!((s21 - s12).norm() < 1e-10, "reciprocity");
        // Verify Z matrix uses z0_ohm = 50.0: Z (I - S) = 50 (I + S)
        let z11 = c(data.z_matrix[k][0][0][0], data.z_matrix[k][0][0][1]);
        let z12 = c(data.z_matrix[k][0][1][0], data.z_matrix[k][0][1][1]);
        let lhs = z11 * (c(1.0, 0.0) - s11) - z12 * s21;
        let rhs = c(50.0, 0.0) * (c(1.0, 0.0) + s11);
        let rel = (lhs - rhs).norm() / rhs.norm().max(1.0);
        assert!(rel < 1e-10, "Z matrix must reference z0_ohm = 50 Ohm: rel err {rel:e}");
        let _ = s22;
    }
}

#[test]
fn fr24_cpw_even_lumped_port_matches_endcap_and_validates_contact() {
    // 1. Compare balanced `cpw_even` (two lateral slot sheets) vs single end-cap lumped port on a 3D CPW through-line
    let (w_box, l_tot, h_sub, h_air) = (640e-6, 1.4e-3, 250e-6, 250e-6);
    let (w_trace, s_gap, l_port) = (100e-6, 60e-6, 120e-6);

    // Model A: single-gap end-cap ports at y in [0, l_port] and [l_tot - l_port, l_tot]
    let mut m_end = Model::new_3d();
    let si_a = m_end.box_("si_box", (-w_box / 2.0, 0.0, -h_sub), (w_box / 2.0, l_tot, 0.0));
    let air_a = m_end.box_("air_box", (-w_box / 2.0, 0.0, 0.0), (w_box / 2.0, l_tot, h_air));
    let tr_a = m_end.rect("tr", (-w_trace / 2.0, l_port), (w_trace / 2.0, l_tot - l_port));
    m_end.sheet("cpw_trace", tr_a, 0.0);
    let p1_a = m_end.rect("p1", (-w_trace / 2.0, 0.0), (w_trace / 2.0, l_port));
    m_end.sheet("port1", p1_a, 0.0);
    let p2_a = m_end.rect("p2", (-w_trace / 2.0, l_tot - l_port), (w_trace / 2.0, l_tot));
    m_end.sheet("port2", p2_a, 0.0);
    let gl_a = m_end.rect("gl", (-w_box / 2.0, 0.0), (-w_trace / 2.0 - s_gap, l_tot));
    m_end.sheet("cpw_gnd_l", gl_a, 0.0);
    let gr_a = m_end.rect("gr", (w_trace / 2.0 + s_gap, 0.0), (w_box / 2.0, l_tot));
    m_end.sheet("cpw_gnd_r", gr_a, 0.0);
    m_end.overlay(&[("si", si_a), ("air", air_a)]);
    let g_end = mesh_model(&m_end.list, &MeshPlan::simple(160e-6)).expect("CPW end-cap mesh");

    // Model B: balanced cpw_even ports across left and right lateral slots near y=0 and y=l_tot
    let mut m_even = Model::new_3d();
    let si_b = m_even.box_("si_box", (-w_box / 2.0, 0.0, -h_sub), (w_box / 2.0, l_tot, 0.0));
    let air_b = m_even.box_("air_box", (-w_box / 2.0, 0.0, 0.0), (w_box / 2.0, l_tot, h_air));
    let tr_b = m_even.rect("tr", (-w_trace / 2.0, l_port), (w_trace / 2.0, l_tot - l_port));
    m_even.sheet("cpw_trace", tr_b, 0.0);
    let gl_b = m_even.rect("gl", (-w_box / 2.0, 0.0), (-w_trace / 2.0 - s_gap, l_tot));
    m_even.sheet("cpw_gnd_l", gl_b, 0.0);
    let gr_b = m_even.rect("gr", (w_trace / 2.0 + s_gap, 0.0), (w_box / 2.0, l_tot));
    m_even.sheet("cpw_gnd_r", gr_b, 0.0);
    // Two disjoint lateral slots for port 1 (y in [l_port, 2*l_port]) and port 2 (y in [l_tot - 2*l_port, l_tot - l_port])
    let p1_l = m_even.rect("p1_l", (-w_trace / 2.0 - s_gap, l_port), (-w_trace / 2.0, 2.0 * l_port));
    let p1_r = m_even.rect("p1_r", (w_trace / 2.0, l_port), (w_trace / 2.0 + s_gap, 2.0 * l_port));
    let p1_both = m_even.union("p1_both", &[p1_l, p1_r]);
    m_even.sheet("port1_even", p1_both, 0.0);
    let p2_l = m_even.rect("p2_l", (-w_trace / 2.0 - s_gap, l_tot - 2.0 * l_port), (-w_trace / 2.0, l_tot - l_port));
    let p2_r = m_even.rect("p2_r", (w_trace / 2.0, l_tot - 2.0 * l_port), (w_trace / 2.0 + s_gap, l_tot - l_port));
    let p2_both = m_even.union("p2_both", &[p2_l, p2_r]);
    m_even.sheet("port2_even", p2_both, 0.0);
    m_even.overlay(&[("si", si_b), ("air", air_b)]);
    let g_even = mesh_model(&m_even.list, &MeshPlan::simple(160e-6)).expect("CPW cpw_even mesh");

    let run_cpw = |g: &MeshedGeom, even: bool| -> DrivenPortData {
        let ex = Mesh::from_solver_mesh(&g.mesh).unwrap();
        let resv = GeomResolver(g);
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({ "f_min": 4.0e9, "f_max": 8.0e9, "n_points": 21, "prom_tol": 1e-3 });
        s.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
        s.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
        for sheet in ["cpw_trace", "cpw_gnd_l", "cpw_gnd_r"] {
            s.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
        }
        if even {
            s.bcs.push(BcSpec {
                on: "port1_even".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "p1", "Z_0": 50.0, "direction": [1.0, 0.0, 0.0], "mode": "cpw_even" }),
            });
            s.bcs.push(BcSpec {
                on: "port2_even".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "p2", "Z_0": 50.0, "direction": [1.0, 0.0, 0.0], "mode": "cpw_even" }),
            });
        } else {
            s.bcs.push(BcSpec {
                on: "port1".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "p1", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0], "mode": "single" }),
            });
            s.bcs.push(BcSpec {
                on: "port2".into(),
                kind: "lumped_port".into(),
                params: json!({ "name": "p2", "Z_0": 50.0, "direction": [0.0, -1.0, 0.0], "mode": "single" }),
            });
        }
        let r = Ems.solve(&ex, &s, &mut Ctx::new(&resv)).unwrap();
        DrivenPortData::of(&r).unwrap()
    };

    let d_end = run_cpw(&g_end, false);
    let d_even = run_cpw(&g_even, true);
    let mut max_db_diff = 0.0_f64;
    for k in 0..d_end.freqs_hz.len() {
        let s21_a = c(d_end.s_matrix[k][1][0][0], d_end.s_matrix[k][1][0][1]);
        let s21_b = c(d_even.s_matrix[k][1][0][0], d_even.s_matrix[k][1][0][1]);
        let s12_b = c(d_even.s_matrix[k][0][1][0], d_even.s_matrix[k][0][1][1]);
        let s11_b = c(d_even.s_matrix[k][0][0][0], d_even.s_matrix[k][0][0][1]);
        let db_a = 20.0 * s21_a.norm().max(1e-30).log10();
        let db_b = 20.0 * s21_b.norm().max(1e-30).log10();
        max_db_diff = max_db_diff.max((db_a - db_b).abs());
        assert!((s21_b - s12_b).norm() < 1e-12, "cpw_even reciprocity");
        assert!((s11_b.norm_sqr() + s21_b.norm_sqr() - 1.0).abs() < 1e-10, "cpw_even lossless unitarity");
    }
    eprintln!("FR-24 cpw_even vs end-cap max |S21| dB diff = {:.4} dB (< 0.2 dB)", max_db_diff);
    assert!(max_db_diff < 0.2, "cpw_even vs end-cap differed by {max_db_diff:.4} dB (>= 0.2 dB)");

    // 2. Disconnected lumped_port with a 10 nm gap from conductor pad -> Error::Invalid (NOT Error::Solver!)
    let mut m_gap = Model::new_3d();
    let si_g = m_gap.box_("si_box", (-500e-6, -500e-6, -200e-6), (500e-6, 500e-6, 0.0));
    let air_g = m_gap.box_("air_box", (-500e-6, -500e-6, 0.0), (500e-6, 500e-6, 200e-6));
    let p_a = m_gap.rect("pad1_r", (-150e-6, -300e-6), (150e-6, -50e-6));
    m_gap.sheet("pad1", p_a, 0.0);
    let p_b = m_gap.rect("pad2_r", (-150e-6, 50e-6 + 10e-9), (150e-6, 300e-6));
    m_gap.sheet("pad2", p_b, 0.0);
    // Sheet touches pad1 at y = -50 um, stops 10 nm short of pad2 at y = +50 um
    let br = m_gap.rect("bridge_r", (-40e-6, -50e-6), (40e-6, 50e-6));
    m_gap.sheet("bridge", br, 0.0);
    m_gap.overlay(&[("si", si_g), ("air", air_g)]);
    let g_gap = mesh_model(&m_gap.list, &MeshPlan::simple(180e-6)).expect("disconnected gap mesh");
    let ex_gap = Mesh::from_solver_mesh(&g_gap.mesh).unwrap();
    let resv_gap = GeomResolver(&g_gap);

    let mut s_gap_drv = ProblemSpec::new("ems", 1, "driven");
    s_gap_drv.study.params = json!({ "f_min": 4e9, "f_max": 6e9, "n_points": 5 });
    s_gap_drv.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    s_gap_drv.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    s_gap_drv.bcs.push(BcSpec { on: "pad1".into(), kind: "pec".into(), params: json!({}) });
    s_gap_drv.bcs.push(BcSpec { on: "pad2".into(), kind: "pec".into(), params: json!({}) });
    s_gap_drv.bcs.push(BcSpec {
        on: "bridge".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "bad_p", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
    });
    let err_drv = Ems.solve(&ex_gap, &s_gap_drv, &mut Ctx::new(&resv_gap)).unwrap_err();
    assert!(
        matches!(err_drv, Error::Invalid(_)),
        "disconnected lumped_port must return Error::Invalid, got {err_drv:?}"
    );

    let mut s_gap_cav = ProblemSpec::new("ems", 1, "cavity");
    s_gap_cav.study.params = json!({ "n_modes": 1, "f_search": 5e9 });
    s_gap_cav.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    s_gap_cav.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    s_gap_cav.bcs.push(BcSpec { on: "pad1".into(), kind: "pec".into(), params: json!({}) });
    s_gap_cav.bcs.push(BcSpec { on: "pad2".into(), kind: "pec".into(), params: json!({}) });
    s_gap_cav.bcs.push(BcSpec {
        on: "bridge".into(),
        kind: "lumped_jj".into(),
        params: json!({ "name": "bad_jj", "E_j_Hz": 20e9, "direction": [0.0, 1.0, 0.0] }),
    });
    let err_cav = Ems.solve(&ex_gap, &s_gap_cav, &mut Ctx::new(&resv_gap)).unwrap_err();
    assert!(matches!(err_cav, Error::Invalid(_)), "disconnected lumped_jj must return Error::Invalid, got {err_cav:?}");

    // 3. Port stamped inside a single continuous PEC sheet (shorted) -> Error::Invalid
    let mut m_short = Model::new_3d();
    let si_s = m_short.box_("si_box", (-500e-6, -500e-6, -200e-6), (500e-6, 500e-6, 0.0));
    let air_s = m_short.box_("air_box", (-500e-6, -500e-6, 0.0), (500e-6, 500e-6, 200e-6));
    let outer = m_short.rect("outer", (-300e-6, -300e-6), (300e-6, 300e-6));
    let inner = m_short.rect("inner", (-60e-6, -60e-6), (60e-6, 60e-6));
    let frame = m_short.difference("frame", outer, &[inner]);
    m_short.sheet("pec_frame", frame, 0.0);
    m_short.sheet("shorted_port", inner, 0.0);
    m_short.overlay(&[("si", si_s), ("air", air_s)]);
    let g_short = mesh_model(&m_short.list, &MeshPlan::simple(180e-6)).expect("shorted port mesh");
    let ex_short = Mesh::from_solver_mesh(&g_short.mesh).unwrap();
    let resv_short = GeomResolver(&g_short);

    let mut s_short = ProblemSpec::new("ems", 1, "driven");
    s_short.study.params = json!({ "f_min": 4e9, "f_max": 6e9, "n_points": 5 });
    s_short.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    s_short.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    s_short.bcs.push(BcSpec { on: "pec_frame".into(), kind: "pec".into(), params: json!({}) });
    s_short.bcs.push(BcSpec {
        on: "shorted_port".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "short_p", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
    });
    let err_short = Ems.solve(&ex_short, &s_short, &mut Ctx::new(&resv_short)).unwrap_err();
    assert!(
        matches!(err_short, Error::Invalid(_)),
        "shorted lumped_port must return Error::Invalid, got {err_short:?}"
    );
}

#[test]
fn fr33_driven_solver_backends_and_factor_backend_logging() {
    let g = build_coupled_feedline_resonator_mesh();
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let make_spec = |solver: &str, blr_tol: f64, prom_tol: serde_json::Value| -> ProblemSpec {
        let mut s = ProblemSpec::new("ems", 1, "driven");
        s.study.params = json!({
            "freqs_hz": [5.0e9, 6.0e9],
            "prom_tol": prom_tol,
            "solver": solver,
            "blr_tol": blr_tol
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
        s
    };

    // Baseline: auto (PROM)
    let res_auto = Ems.solve(&exchange, &make_spec("auto", 1e-6, json!(1e-3)), &mut Ctx::new(&resolver)).unwrap();
    let d_auto = DrivenPortData::of(&res_auto).unwrap();
    assert!(
        res_auto.log.iter().any(|l| l.contains("solver=auto") && l.contains("backend=")),
        "auto must log solver & backend: {:?}",
        res_auto.log
    );

    // Direct, LBLT, BLR, and AMS+FGMRES must match auto on S-parameters and log their backend
    for (solver_name, expected_tag, prom_tol) in [
        ("direct", "backend=", json!(1e-3)),
        ("lblt", "backend=", json!(1e-3)),
        ("blr", "backend=BlrNd", json!(1e-3)),
        ("ams_fgmres", "backend=AmsFgmres", json!(null)),
    ] {
        let res = Ems.solve(&exchange, &make_spec(solver_name, 1e-6, prom_tol), &mut Ctx::new(&resolver)).unwrap();
        let d = DrivenPortData::of(&res).unwrap();
        assert!(
            res.log.iter().any(|l| l.contains(&format!("solver={solver_name}")) && l.contains(expected_tag)),
            "{solver_name} log missing expected tag '{expected_tag}': {:?}",
            res.log
        );
        if matches!(solver_name, "blr" | "ams_fgmres") {
            assert_eq!(
                res.scalar("symbolic_count"),
                Some(1.0),
                "{solver_name} must reuse cached symbolic/topology across frequencies (hdm_solves = {})",
                d.hdm_solves
            );
        }
        let mut max_diff = 0.0_f64;
        for k in 0..d_auto.freqs_hz.len() {
            for p in 0..2 {
                for q in 0..2 {
                    let s0 = c(d_auto.s_matrix[k][p][q][0], d_auto.s_matrix[k][p][q][1]);
                    let s1 = c(d.s_matrix[k][p][q][0], d.s_matrix[k][p][q][1]);
                    max_diff = max_diff.max((s1 - s0).norm());
                }
            }
        }
        eprintln!("FR-33 driven solver={solver_name}: max |S - S_auto| = {max_diff:.3e}, log = {:?}", res.log);
        assert!(max_diff < 1e-3, "FR-33 solver={solver_name} S-matrix diff = {max_diff:.3e} >= 1e-3");
    }

    // Order-2 P2 AMS+FGMRES parity check (`AmsFgmresSolver::new_shifted_complex_p2` + `update_matrix`)
    let make_p2_spec = |solver: &str| -> ProblemSpec {
        let mut s = make_spec(solver, 1e-6, json!(null));
        s.order = 2;
        s
    };
    let res_p2_auto = Ems.solve(&exchange, &make_p2_spec("auto"), &mut Ctx::new(&resolver)).unwrap();
    let res_p2_ams = Ems.solve(&exchange, &make_p2_spec("ams_fgmres"), &mut Ctx::new(&resolver)).unwrap();
    assert_eq!(res_p2_ams.scalar("symbolic_count"), Some(1.0));
    assert!(
        res_p2_ams.log.iter().any(|l| l.contains("solver=ams_fgmres") && l.contains("backend=AmsFgmres")),
        "P2 ams_fgmres log: {:?}",
        res_p2_ams.log
    );
    let d_p2_auto = DrivenPortData::of(&res_p2_auto).unwrap();
    let d_p2_ams = DrivenPortData::of(&res_p2_ams).unwrap();
    let mut max_p2_diff = 0.0_f64;
    for k in 0..d_p2_auto.freqs_hz.len() {
        for p in 0..2 {
            for q in 0..2 {
                let s0 = c(d_p2_auto.s_matrix[k][p][q][0], d_p2_auto.s_matrix[k][p][q][1]);
                let s1 = c(d_p2_ams.s_matrix[k][p][q][0], d_p2_ams.s_matrix[k][p][q][1]);
                max_p2_diff = max_p2_diff.max((s1 - s0).norm());
            }
        }
    }
    eprintln!("FR-33 P2 driven solver=ams_fgmres: max |S - S_auto| = {max_p2_diff:.3e}, log = {:?}", res_p2_ams.log);
    assert!(max_p2_diff < 1e-3, "FR-33 P2 solver=ams_fgmres S-matrix diff = {max_p2_diff:.3e} >= 1e-3");

    // Validation: unknown solver name is rejected
    let err_bad =
        Ems.solve(&exchange, &make_spec("bogus_solver", 1e-6, json!(null)), &mut Ctx::new(&resolver)).unwrap_err();
    assert!(matches!(err_bad, Error::Invalid(_)), "unknown solver must return Error::Invalid, got {err_bad:?}");
}
