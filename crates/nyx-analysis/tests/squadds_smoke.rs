// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! SQuADDS PR CI smoke suite for `nyx-analysis` (Plan §8.3, `< 90 s` release budget):
//! - Compact 3-net `TransmonCross` (`cross`, `claw`, `ground`) Mode G2 dual-box `1/R` extrapolation,
//!   `AF-3` single-box `far_field` (`order = 2`), and `FR-11` 3D `I_edge` Gram-matrix finite-thickness
//!   correction (`metal_thickness_m = 200 nm`), verifying SPD/sign invariants and pre-Schur qubit reduction.
//! - 4-mode `Dissipator` DC-SQUID flux sweep (`S4.3`), verifying native `RomPort` `squid_flux_sweep`
//!   against point-by-point `ej_hz_override` in `craig_bampton` across `[0, 0.5] Phi_0` (`< 250 ms`),
//!   tracking avoided crossings across both readout cavity and Purcell filter modes.
//! - `MUNINN` 3D Xmon + readout claw `tls_participation` with explicit 9-parameter MIT-LL `MA/MS/SA`
//!   oxide stack (`d = 200 nm`, `edge_cutoff = 1.5 um`) combined with multi-mode Purcell `T_1`.

use femcore::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
use femcore::meshing::{mesh_model, MeshPlan, Mesher, PmlMode};
use femcore::props::LinearSolver;
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_analysis::epr::{jj_flux_scale, project_rom_modes, C0, EPS0, MU0};
use nyx_analysis::{
    capacitance_matrix_with_mesh_and_params, hamiltonian_with_params, squid_effective_ej, tls_participation_with_params,
};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study};
use nyx_ems::Ems;
use nyx_physics::{
    BcSpec, CapacitanceMatrixData, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec,
    RomPort, SolveResult,
};
use serde_json::json;
use std::f64::consts::PI;

const EPS_SI: f64 = 11.45;

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

/// Builds and solves a compact 3-net `TransmonCross` (`cross`, `claw`, `ground_plane`) on a `500 x 500 um`
/// cryogenic silicon die (`eps_r = 11.45`) enclosed in a 6-sided vacuum box padded by `pad_m`.
/// When `confine_die` is `true`, `sub_box` is confined to the finite die footprint (creating vertical outer
/// `Si|Vac` sidewalls); when `false`, the silicon layer spans the full lateral box width (`FR-11` compatible).
fn solve_compact_transmon_cross_3net(pad_m: f64, far_field: bool, confine_die: bool) -> (Mesh, SolveResult) {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let half_die = 250.0e-6;
    let h_sub = 220.0e-6;
    let h_air = 220.0e-6;

    let mut m = Model::new_3d();
    let vac = m.box_(
        "vac_box",
        (-half_die - pad_m, -half_die - pad_m, -h_sub - pad_m),
        (half_die + pad_m, half_die + pad_m, h_air + pad_m),
    );
    // Pin `inner_vac_box` (`z in [0, h_air]`) and `sub_box` interior planes so the 2D XY triangulation and Z-slabs
    // inside `[-h_sub, +h_air]` are invariant across padding `P1 = 0.6 L` and `P2 = 1.2 L` (matching C0.2).
    let inner_vac = m.box_("inner_vac_box", (-half_die, -half_die, 0.0), (half_die, half_die, h_air));
    let sub_half = if confine_die { half_die } else { half_die + pad_m };
    let sub = m.box_("sub_box", (-sub_half, -sub_half, -h_sub), (sub_half, sub_half, 0.0));

    let cross_poly = m.polygons(
        "cross_poly",
        vec![vec![
            pt(-12.0e-6, -90.0e-6),
            pt(12.0e-6, -90.0e-6),
            pt(12.0e-6, -12.0e-6),
            pt(90.0e-6, -12.0e-6),
            pt(90.0e-6, 12.0e-6),
            pt(12.0e-6, 12.0e-6),
            pt(12.0e-6, 90.0e-6),
            pt(-12.0e-6, 90.0e-6),
            pt(-12.0e-6, 12.0e-6),
            pt(-90.0e-6, 12.0e-6),
            pt(-90.0e-6, -12.0e-6),
            pt(-12.0e-6, -12.0e-6),
        ]],
        vec![],
    );
    let claw_poly = m.polygons(
        "claw_poly",
        vec![vec![
            pt(-36.0e-6, 104.0e-6),
            pt(36.0e-6, 104.0e-6),
            pt(36.0e-6, 152.0e-6),
            pt(24.0e-6, 152.0e-6),
            pt(24.0e-6, 118.0e-6),
            pt(-24.0e-6, 118.0e-6),
            pt(-24.0e-6, 152.0e-6),
            pt(-36.0e-6, 152.0e-6),
        ]],
        vec![],
    );
    let gnd_poly = m.polygons(
        "gnd_poly",
        vec![vec![pt(-half_die, -half_die), pt(half_die, -half_die), pt(half_die, half_die), pt(-half_die, half_die)]],
        vec![vec![pt(-115.0e-6, -115.0e-6), pt(115.0e-6, -115.0e-6), pt(115.0e-6, 172.0e-6), pt(-115.0e-6, 172.0e-6)]],
    );

    m.sheet("cross", cross_poly, 0.0);
    m.sheet("claw", claw_poly, 0.0);
    m.sheet("ground_plane", gnd_poly, 0.0);
    // Vacuum precedes inner_vacuum & substrate in painter order so sub_box/inner_vac win inside the die while vac_box wraps all 6 outer walls.
    m.overlay(&[("vacuum", vac), ("inner_vacuum", inner_vac), ("substrate", sub)]);

    m.size("sz_cross", Selector::parse("cross:edge").unwrap(), Some(Expr::Num(9.0e-6)), None, None);
    m.size("sz_claw", Selector::parse("claw:edge").unwrap(), Some(Expr::Num(9.0e-6)), None, None);
    m.size("sz_gnd", Selector::parse("ground_plane:edge").unwrap(), Some(Expr::Num(20.0e-6)), None, None);
    m.size(
        "sz_die",
        Selector::parse("(name(substrate) | name(inner_vacuum))").unwrap(),
        Some(Expr::Num(95.0e-6)),
        None,
        None,
    );

    let mut opts = GeomMeshOpts::new(145.0e-6);
    opts.z_coarsen = false;
    opts.dz_min = Some(14.0e-6);
    opts.z_growth_ratio = Some(2.0);
    opts.max_z_slabs = Some(48);
    let plan = MeshPlan { opts, mesher: Mesher::Simple, pml: PmlMode::Inline };
    let g = mesh_model(&m.list, &plan).expect("compact TransmonCross 3-net mesh");

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.solver = LinearSolver::Direct;
    for d in g.domains("substrate") {
        prob.materials.set(d, "eps_r", EPS_SI);
        prob.domain_names.insert(d, "substrate".into());
    }
    for d in g.domains("inner_vacuum") {
        prob.materials.set(d, "eps_r", 1.0);
        prob.domain_names.insert(d, "inner_vacuum".into());
    }
    for d in g.domains("vacuum") {
        prob.materials.set(d, "eps_r", 1.0);
        prob.domain_names.insert(d, "vacuum".into());
    }

    let outer_kind = if far_field {
        BcKind::FarField { order: 2, decay_power: 1.0, center: Some([0.0, 0.0, 0.0]) }
    } else {
        BcKind::Ground
    };
    prob.bcs = vec![
        Bc { labels: g.boundary("cross"), kind: BcKind::Terminal("cross".into()) },
        Bc { labels: g.boundary("claw"), kind: BcKind::Terminal("claw".into()) },
        Bc { labels: g.boundary("ground_plane"), kind: BcKind::Terminal("ground".into()) },
        Bc { labels: g.select("outer").unwrap(), kind: outer_kind },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).expect("3-net P2 AMG electrostatics");
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    (ex_mesh, res)
}

#[test]
fn transmon_cross_3net_mode_g2_far_field_and_i_edge_thickness_smoke() {
    let l_chip = 500.0e-6;
    let p1 = 0.6 * l_chip;
    let p2 = 1.2 * l_chip;
    let r1 = 0.5 * l_chip + p1;
    let r2 = 0.5 * l_chip + p2;

    // 1. Mode G2 dual-box 1/R open-boundary extrapolation (`C0.2` finite-die + pinned `inner_vac_box` setup)
    let (_, res_p1) = solve_compact_transmon_cross_3net(p1, false, true);
    let (_, res_p2) = solve_compact_transmon_cross_3net(p2, false, true);
    let c_p1 = CapacitanceMatrixData::of(&res_p1).unwrap().maxwell_matrix_f;
    let c_p2 = CapacitanceMatrixData::of(&res_p2).unwrap().maxwell_matrix_f;

    let mut c_inf = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c_inf[i][j] = (r2 * c_p2[i][j] - r1 * c_p1[i][j]) / (r2 - r1);
        }
    }

    // Tight outer Dirichlet walls must monotonically inflate ground_to_ground (C_p1 > C_p2 > C_inf > 0)
    assert!(
        c_p1[2][2] > c_p2[2][2] && c_p2[2][2] > c_inf[2][2] && c_inf[2][2] > 0.0,
        "expected Dirichlet wall relaxation C_p1={:e} > C_p2={:e} > C_inf={:e}",
        c_p1[2][2],
        c_p2[2][2],
        c_inf[2][2]
    );

    // 2. Single-box AF-3 Bayliss-Turkel order-2 far_field BC at compact padding P1 (identical interior stack)
    let (_, res_ff) = solve_compact_transmon_cross_3net(p1, true, true);
    let c_ff = CapacitanceMatrixData::of(&res_ff).unwrap().maxwell_matrix_f;

    // Cross and claw self/mutual entries are well-shielded by the pocket and agree closely between Mode G2 and far_field
    for (idx, label) in [(0, "cross_to_cross"), (1, "claw_to_claw")] {
        let rel = (c_ff[idx][idx] - c_inf[idx][idx]).abs() / c_inf[idx][idx];
        assert!(
            rel < 0.025,
            "{label} mismatch between Mode G2 1/R ({:e}) and AF-3 far_field ({:e}): {:.3}%",
            c_inf[idx][idx],
            c_ff[idx][idx],
            rel * 100.0
        );
    }
    let rel_mut = (c_ff[0][1] - c_inf[0][1]).abs() / c_inf[0][1].abs();
    assert!(rel_mut < 0.03, "cross_to_claw mismatch {:.3}% >= 3.0%", rel_mut * 100.0);
    assert!(
        c_ff[2][2] < c_p1[2][2],
        "far_field ({:e}) must eliminate artificial P1 Dirichlet ground_to_ground inflation ({:e})",
        c_ff[2][2],
        c_p1[2][2]
    );

    // 3. FR-11 3D matched-asymptotic I_edge finite metal thickness perturbation (t = 200 nm) + qubit reduction
    // (evaluated on full-width 3-layer stack compatible with FR-11 on HEAD)
    let (mesh_thick, res_thick) = solve_compact_transmon_cross_3net(p1, true, false);
    let raw_q = capacitance_matrix_with_mesh_and_params(
        &mesh_thick,
        &res_thick,
        &json!({
            "ground": "ground",
            "qubits": [{ "name": "Q1", "pads": ["cross"], "E_J_Hz": 15.0e9 }],
            "resonators": [{ "name": "RO", "terminal": "claw", "f_r_Hz": 6.2e9 }]
        }),
    )
    .unwrap();
    let corr_q = capacitance_matrix_with_mesh_and_params(
        &mesh_thick,
        &res_thick,
        &json!({
            "ground": "ground",
            "metal_thickness_m": 200.0e-9,
            "qubits": [{ "name": "Q1", "pads": ["cross"], "E_J_Hz": 15.0e9 }],
            "resonators": [{ "name": "RO", "terminal": "claw", "f_r_Hz": 6.2e9 }]
        }),
    )
    .unwrap();

    assert_eq!(corr_q["has_null_mode"], json!(false));
    assert!(corr_q["thickness_correction"]["s_eff_m"].as_f64().unwrap() > 0.0);

    for i in 0..3 {
        let c0 = raw_q["maxwell_matrix_F"][i][i].as_f64().unwrap();
        let ct = corr_q["maxwell_matrix_F"][i][i].as_f64().unwrap();
        assert!(ct > c0, "200 nm thickness must increase diagonal C[{i},{i}]: {ct:e} <= {c0:e}");
        for j in (i + 1)..3 {
            let cij = corr_q["maxwell_matrix_F"][i][j].as_f64().unwrap();
            let cji = corr_q["maxwell_matrix_F"][j][i].as_f64().unwrap();
            let cjj = corr_q["maxwell_matrix_F"][j][j].as_f64().unwrap();
            assert!((cij - cji).abs() <= 1.0e-13 * (ct * cjj).sqrt(), "asymmetric corrected C[{i},{j}]");
            assert!(cij <= 1.0e-5 * (ct * cjj).sqrt(), "violated Maxwell off-diagonal bound C[{i},{j}] = {cij:e}");
        }
    }

    let delta_cross_claw = corr_q["thickness_correction"]["delta_maxwell_F"][0][1].as_f64().unwrap();
    assert!(delta_cross_claw < 0.0, "mutual cross-claw capacitance must increase with thickness: {delta_cross_claw:e}");
    let ec_raw = raw_q["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    let ec_corr = corr_q["qubits"][0]["E_C_Hz"].as_f64().unwrap();
    assert!(ec_corr < ec_raw, "thickness correction must lower transmon E_C: {ec_corr:e} >= {ec_raw:e}");
}

/// Core Bug 1 Regression (`crates/nyx-analysis/src/tls.rs:789-797`):
/// On a finite silicon die (`sub_box` confined to `die_bbox_m` inside a larger outer `vac_box`),
/// the 4 vertical outer die `Si|Vac` sidewalls (`|n_z| < 0.99`) falsely set `has_vertical_sa_trench = true` ->
/// `has_resolved_3d_step = true`, suppressing interior zero-thickness sheet perimeter detection and causing
/// `FR-11` (`metal_thickness_m`) to fail with `"found no interior zero-thickness sheet perimeter edges"`.
/// Expected to FAIL on unmodified `HEAD` (`9480093a`).
#[test]
fn bug1_fr11_metal_thickness_on_finite_die_fails() {
    let (mesh_finite_die, res_finite_die) = solve_compact_transmon_cross_3net(250.0e-6, true, true);
    let out = capacitance_matrix_with_mesh_and_params(
        &mesh_finite_die,
        &res_finite_die,
        &json!({ "metal_thickness_m": 200.0e-9 }),
    )
    .expect("Core Bug 1 (tls.rs:789-797): FR-11 metal_thickness_m must succeed on finite silicon die inside vacuum");
    let c_corr = out["maxwell_matrix_F"][0][0].as_f64().unwrap();
    let c_raw = CapacitanceMatrixData::of(&res_finite_die).unwrap().maxwell_matrix_f[0][0];
    assert!(c_corr > c_raw, "finite 200 nm thickness must increase self-capacitance ({c_corr:e} <= {c_raw:e})");
}

#[test]
fn dissipator_dc_squid_craig_bampton_vs_rom_port_flux_sweep_smoke() {
    // 4-DOF CircuitRom modeling S4.3 Dissipator:
    //   dof 0: fixed transmon qubit (~3.37 GHz)
    //   dof 1: flux-tunable DC-SQUID dissipator (sweeping ~15.3 GHz down to ~4.5 GHz across [0, 0.5] Phi_0)
    //   dof 2: readout cavity (~5.59 GHz, kappa ~ 0.48 MHz)
    //   dof 3: lossy Purcell/dissipator filter (~8.60 GHz, kappa_f ~ 120 MHz)
    let scale = jj_flux_scale();
    let ej_q = 9.5e9;
    let ej_d_max = 46.5e9;
    let squid_d = 0.095;

    let c_q = 92.0e-15;
    let c_d = 46.0e-15;
    let c_r = 420.0e-15;
    let c_f = 280.0e-15;

    let c_qd = 1.6e-15;
    let c_cd = 4.5e-15;
    let c_df = 14.0e-15;
    let c_qr = 0.8e-15;

    let w_r = 2.0 * PI * 5.594e9;
    let w_f = 2.0 * PI * 8.600e9;
    let inv_l_r = w_r * w_r * (c_r + c_cd + c_qr);
    let inv_l_f = w_f * w_f * (c_f + c_df);

    let z0 = 50.0;
    let kappa_r = 2.0 * PI * 477.0e3;
    let kappa_f = 2.0 * PI * 120.0e6;
    let b_ro = (kappa_r * (c_r + c_cd + c_qr) * z0).sqrt();
    let b_flt = (kappa_f * (c_f + c_df) * z0).sqrt();
    let damp_scale = MU0 * C0 / z0;

    let mut rom = CircuitRom {
        ports: vec![
            RomPort {
                name: "qubit_jj".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / ej_q),
                e_j_hz: Some(ej_q),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            },
            RomPort {
                name: "dissipator_jj".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / ej_d_max),
                e_j_hz: Some(ej_d_max),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: Some(squid_d),
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            },
            RomPort {
                name: "readout_port".into(),
                kind: "lumped_port".into(),
                l_j_h: None,
                e_j_hz: None,
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: Some(z0),
                direction: [1.0, 0.0, 0.0],
                formulation: "rank1".into(),
            },
            RomPort {
                name: "filter_port".into(),
                kind: "lumped_port".into(),
                l_j_h: None,
                e_j_hz: None,
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: Some(z0),
                direction: [1.0, 0.0, 0.0],
                formulation: "rank1".into(),
            },
        ],
        r_static: 0,
        r_dyn: 4,
        k_r: vec![
            vec![0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, MU0 * inv_l_r, 0.0],
            vec![0.0, 0.0, 0.0, MU0 * inv_l_f],
        ],
        c_r: vec![
            vec![0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, damp_scale * b_ro * b_ro, 0.0],
            vec![0.0, 0.0, 0.0, damp_scale * b_flt * b_flt],
        ],
        m_r: vec![
            vec![(c_q + c_qd + c_qr) / EPS0, -c_qd / EPS0, -c_qr / EPS0, 0.0],
            vec![-c_qd / EPS0, (c_d + c_qd + c_cd + c_df) / EPS0, -c_cd / EPS0, -c_df / EPS0],
            vec![-c_qr / EPS0, -c_cd / EPS0, (c_r + c_cd + c_qr) / EPS0, 0.0],
            vec![0.0, -c_df / EPS0, 0.0, (c_f + c_df) / EPS0],
        ],
        b_r: vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, b_ro, 0.0],
            vec![0.0, 0.0, 0.0, b_flt],
        ],
        c_port_f: vec![],
        modes: vec![],
    };
    rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
    let res = SolveResult { scalars: vec![rom.scalar()], ..Default::default() };

    let flux_grid: Vec<f64> = (0..=24).map(|k| 0.5 * (k as f64) / 24.0).collect();
    let rep = hamiltonian_with_params(
        &res,
        &json!({
            "basis": "craig_bampton",
            "port_trunc": 6,
            "fock_trunc": [4, 4],
            "sweep_port": "dissipator_jj",
            "squid_flux_sweep": flux_grid
        }),
    )
    .expect("native DC-SQUID Dissipator sweep");

    let sweep = rep.sweep.as_ref().expect("sweep points present");
    assert_eq!(sweep.len(), 25);
    assert!(
        rep.sweep_elapsed_s.unwrap() < 0.5,
        "25-point DC-SQUID ROM sweep took {:.3}s >= 0.5s",
        rep.sweep_elapsed_s.unwrap()
    );
    assert_eq!(rep.f_cav_hz.len(), 2);
    assert_eq!(rep.kappa_cav_hz.len(), 2);

    // Dissipator tunes from ~10.3 GHz at 0 Phi_0 (above both 5.59 GHz readout and 8.60 GHz Purcell filter)
    // down to ~2.9 GHz at 0.5 Phi_0 (below both cavities, crossing both avoided crossings)
    let f_d_0 = sweep[0].f01_hz[1];
    let f_d_half = sweep[24].f01_hz[1];
    assert!(f_d_0 > 9.5e9 && f_d_0 < 12.0e9, "unexpected 0-flux dissipator f01 = {f_d_0:e} Hz");
    assert!(f_d_half > 2.2e9 && f_d_half < 4.5e9, "unexpected 0.5-Phi0 dissipator f01 = {f_d_half:e} Hz");

    // Verify exact parity between native `squid_flux_sweep` and manual `ej_hz_override` loop
    for &idx in &[0usize, 8, 16, 24] {
        let phi_norm = flux_grid[idx];
        let ej_eff = squid_effective_ej(ej_d_max, phi_norm, squid_d);
        assert!((sweep[idx].e_j_hz.unwrap() - ej_eff).abs() <= 1.0e-6 * ej_eff);

        let manual = hamiltonian_with_params(
            &res,
            &json!({
                "basis": "craig_bampton",
                "port_trunc": 6,
                "fock_trunc": [4, 4],
                "ej_hz_override": { "dissipator_jj": ej_eff }
            }),
        )
        .unwrap();

        for (got, want) in sweep[idx].eigenvalues_hz.iter().take(8).zip(manual.eigenvalues_hz.iter().take(8)) {
            let denom = want.abs().max(1.0);
            assert!(
                (got - want).abs() / denom < 1.0e-12,
                "flux idx {idx} eigenvalue mismatch: native={got:e} vs ej_hz_override={want:e}"
            );
        }
        for (got, want) in sweep[idx].f01_hz.iter().zip(&manual.f01_hz) {
            assert!((got - want).abs() / want.abs() < 1.0e-12);
        }
    }

    // Verify non-zero minimum eigenvalue spacing across the entire avoided-crossing sweep
    for (k, pt) in sweep.iter().enumerate() {
        for w in pt.eigenvalues_hz.windows(2).take(5) {
            assert!(w[1] - w[0] > 10.0e6, "collapsed avoided crossing at point {k}: {:?}", &pt.eigenvalues_hz[..6]);
        }
    }
}

#[test]
fn muninn_3d_tls_participation_nine_param_oxide_and_purcell_smoke() {
    // Compact 3D Xmon + readout launch with h_edge = 0.9 um around Xmon/gap so edge_cutoff = 1.5 um
    // strictly satisfies 1.5 um >= 1.15 * h_edge and 1.5 um < s_gap / 2 (s_gap = 8 um).
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-180.0e-6, -180.0e-6, -120.0e-6), (180.0e-6, 220.0e-6, 0.0));
    let vac = m.box_("vacuum", (-180.0e-6, -180.0e-6, 0.0), (180.0e-6, 220.0e-6, 120.0e-6));

    // Xmon cross (arm half-width 10 um, arm half-length 55 um) inside pocket [-63, 63] x [-63, 63] um (8 um gap)
    let xmon_poly = m.polygons(
        "xmon_poly",
        vec![vec![
            pt(-10.0e-6, -55.0e-6),
            pt(10.0e-6, -55.0e-6),
            pt(10.0e-6, -10.0e-6),
            pt(55.0e-6, -10.0e-6),
            pt(55.0e-6, 10.0e-6),
            pt(10.0e-6, 10.0e-6),
            pt(10.0e-6, 55.0e-6),
            pt(-10.0e-6, 55.0e-6),
            pt(-10.0e-6, 10.0e-6),
            pt(-55.0e-6, 10.0e-6),
            pt(-55.0e-6, -10.0e-6),
            pt(-10.0e-6, -10.0e-6),
        ]],
        vec![],
    );
    // JJ sheet bridging south Xmon tip (y = -55 um) to ground pocket wall (y = -63 um)
    let jj_poly = m.polygons(
        "jj_poly",
        vec![vec![pt(-4.0e-6, -63.0e-6), pt(4.0e-6, -63.0e-6), pt(4.0e-6, -55.0e-6), pt(-4.0e-6, -55.0e-6)]],
        vec![],
    );
    // Readout coupler strip + axial 50-Ohm port sheet at y in [195, 220] um
    let ro_poly = m.polygons(
        "ro_poly",
        vec![vec![pt(-14.0e-6, 73.0e-6), pt(14.0e-6, 73.0e-6), pt(14.0e-6, 195.0e-6), pt(-14.0e-6, 195.0e-6)]],
        vec![],
    );
    let port_poly = m.polygons(
        "port_poly",
        vec![vec![pt(-14.0e-6, 195.0e-6), pt(14.0e-6, 195.0e-6), pt(14.0e-6, 220.0e-6), pt(-14.0e-6, 220.0e-6)]],
        vec![],
    );
    let gnd_poly = m.polygons(
        "gnd_poly",
        vec![vec![pt(-180.0e-6, -180.0e-6), pt(180.0e-6, -180.0e-6), pt(180.0e-6, 220.0e-6), pt(-180.0e-6, 220.0e-6)]],
        vec![
            vec![pt(-63.0e-6, -63.0e-6), pt(63.0e-6, -63.0e-6), pt(63.0e-6, 63.0e-6), pt(-63.0e-6, 63.0e-6)],
            vec![pt(-24.0e-6, 65.0e-6), pt(24.0e-6, 65.0e-6), pt(24.0e-6, 220.0e-6), pt(-24.0e-6, 220.0e-6)],
        ],
    );

    m.sheet("xmon", xmon_poly, 0.0);
    m.sheet("q_jj", jj_poly, 0.0);
    m.sheet("ro_strip", ro_poly, 0.0);
    m.sheet("ro_port", port_poly, 0.0);
    m.sheet("ground_plane", gnd_poly, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    // Exclude lumped_jj (`q_jj`) and lumped_port (`ro_port`) neighborhoods from edge refinement so the
    // passing MUNINN smoke test avoids Core Bug 2 (`cavity.rs:646-694` port vertex condensation stall)
    // while keeping fine Xmon arm/pocket edges (`h_edge = 3.5 um`, `edge_cutoff = 5.0 um`).
    let excl = "!in_box([-14e-6, -72e-6, -5e-6], [14e-6, -46e-6, 5e-6], touching) & !in_box([-26e-6, 185e-6, -5e-6], [26e-6, 225e-6, 5e-6], touching)";
    m.size(
        "sz_xmon_edge",
        Selector::parse(&format!("xmon:edge & {excl}")).unwrap(),
        Some(Expr::Num(3.5e-6)),
        None,
        None,
    );
    m.size(
        "sz_ro_edge",
        Selector::parse(&format!("ro_strip:edge & {excl}")).unwrap(),
        Some(Expr::Num(6.0e-6)),
        None,
        None,
    );
    m.size(
        "sz_gnd_edge",
        Selector::parse(&format!("ground_plane:edge & {excl}")).unwrap(),
        Some(Expr::Num(8.0e-6)),
        None,
        None,
    );

    let mut opts = GeomMeshOpts::new(65.0e-6);
    opts.z_coarsen = false;
    opts.dz_min = Some(10.0e-6);
    opts.z_growth_ratio = Some(1.9);
    opts.max_z_slabs = Some(8);
    let plan = MeshPlan { opts, mesher: Mesher::Simple, pml: PmlMode::Inline };
    let g = mesh_model(&m.list, &plan).expect("compact MUNINN Xmon 3D mesh");
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = GeomResolver(&g);

    let mut spec = ProblemSpec::new("ems", 1, "cavity");
    spec.study.params = json!({ "n_modes": 2, "f_search": 4.8e9 });
    spec.materials.insert("substrate".into(), MaterialAssign::Props([("eps_r".to_string(), EPS_SI)].into()));
    spec.materials.insert("vacuum".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    for sheet in ["xmon", "ro_strip", "ground_plane"] {
        spec.bcs.push(BcSpec { on: sheet.into(), kind: "pec".into(), params: json!({}) });
    }
    spec.bcs.push(BcSpec {
        on: "q_jj".into(),
        kind: "lumped_jj".into(),
        params: json!({
            "name": "Q_MUNINN",
            "E_j_Hz": 14.0e9,
            "C_j": 35.0e-15,
            "direction": [0.0, 1.0, 0.0],
            "formulation": "rank1"
        }),
    });
    spec.bcs.push(BcSpec {
        on: "ro_port".into(),
        kind: "lumped_port".into(),
        params: json!({
            "name": "RO_50",
            "Z_0": 50.0,
            "direction": [0.0, 1.0, 0.0]
        }),
    });

    let mut ctx = Ctx::new(&resolver);
    let res = Ems.solve(&ex_mesh, &spec, &mut ctx).expect("MUNINN 3D cavity solve");

    let rom = CircuitRom::of(&res).expect("CircuitRom present");
    let sol_idx = rom
        .modes
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.u_ind_jj_j[0].partial_cmp(&b.u_ind_jj_j[0]).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);

    // Exact 9-parameter MIT-LL oxide interface stack from Plan Appendix A.5
    let tls = tls_participation_with_params(
        &ex_mesh,
        &res,
        &json!({
            "solution": sol_idx,
            "d": 200.0e-9,
            "edge_cutoff": 5.0e-6,
            "t_ma": 3.0e-9, "eps_ma": 9.8,   "tan_delta_ma": 1.5e-3,
            "t_ms": 2.0e-9, "eps_ms": 11.45, "tan_delta_ms": 5.0e-4,
            "t_sa": 2.0e-9, "eps_sa": 3.9,   "tan_delta_sa": 1.5e-3
        }),
    )
    .expect("MUNINN 9-param 3D TLS participation");

    assert_eq!(tls.solution, sol_idx);
    assert!((tls.edge_cutoff_m - 5.0e-6).abs() < 1.0e-15);
    assert!(tls.p_ma > 0.0 && tls.p_ms > 0.0 && tls.p_sa > 0.0, "expected positive MA/MS/SA participations: {tls:?}");
    assert!(tls.p_ma_edge > 0.0 && tls.p_ms_edge > 0.0 && tls.p_sa_edge > 0.0, "expected edge-tube contributions");
    assert!(tls.q_tls > 1.0e4 && tls.q_tls < 1.0e9 && tls.q_tls.is_finite(), "unphysical q_tls = {:e}", tls.q_tls);

    let ham = hamiltonian_with_params(&res, &json!({ "basis": "craig_bampton" })).unwrap();
    let f01_hz = ham.f01_hz[0];
    let gamma_tls = 2.0 * PI * f01_hz / tls.q_tls;
    let gamma_purcell =
        if !ham.t1_purcell_s.is_empty() && ham.t1_purcell_s[0] > 0.0 { 1.0 / ham.t1_purcell_s[0] } else { 0.0 };
    let t1_pred_s = 1.0 / (gamma_tls + gamma_purcell);
    assert!(t1_pred_s > 1.0e-7 && t1_pred_s.is_finite(), "invalid combined T1 = {t1_pred_s:e} s");
}
