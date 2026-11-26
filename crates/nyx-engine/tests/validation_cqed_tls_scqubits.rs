// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 ASME V&V 20 validation suite for frozen `scqubits` Hamiltonians and cryogenic `Q_TLS`:
//! - Frozen high-cutoff `scqubits` JSON Fluxonium spectra (`Phi_ext = 0, 0.25, 0.5 Phi_0`) and Transmon vs
//!   `nyx-analysis` `HamiltonianAnalysis` (`craig_bampton` & `normal_fock`)
//! - Calusine (2018) / Woods (2019) / Place (2021) cryogenic CPW dielectric interface participations
//!   (`p_ma, p_ms, p_sa`, `Q_TLS`) with 3-grid Roache GCI `u_num` and trench-depth monotonicity

mod common;

use common::{load_validation_json, scalar_exact, u_num_three_grid, verify_asme_vv20, ValidationTarget};
use femcore::fem::{facet_rule, geometry, quad_rule, shape};
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_analysis::epr::{jj_flux_scale, project_rom_modes, E_CHARGE, MU0, PLANCK_H};
use nyx_analysis::{hamiltonian_with_params, tls_participation_with_params};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study, EPS0};
use nyx_physics::{CircuitRom, Mesh, RomPort, SolveResult};
use serde::Deserialize;
use serde_json::json;
use std::f64::consts::PI;

#[derive(Deserialize)]
struct ScqTlsDataset {
    frozen_scqubits_fluxonium: FluxoniumSection,
    frozen_scqubits_transmon: TransmonSection,
    cryogenic_q_tls_cpw: TlsSection,
}

#[derive(Deserialize)]
struct FluxPoint {
    phi_ext_over_phi0: f64,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct FluxoniumSection {
    ec_hz: f64,
    el_hz: f64,
    ej_hz: f64,
    f_r_hz: f64,
    flux_points: Vec<FluxPoint>,
}

#[derive(Deserialize)]
struct TransmonSection {
    ec_hz: f64,
    ej_hz: f64,
    targets: Vec<ValidationTarget>,
}

#[derive(Deserialize)]
struct TlsGeom {
    w_trace_m: f64,
    g_gap_m: f64,
    d_film_m: f64,
    t_ox_m: f64,
    eps_si: f64,
    eps_sapphire_place: f64,
    eps_ox: f64,
    tan_delta: f64,
    edge_cutoff_m: f64,
}

#[derive(Deserialize)]
struct TlsTol {
    interface_rel: f64,
    u_input_rel: f64,
    u_d_rel: f64,
    trench_depths_m: Vec<f64>,
}

#[derive(Deserialize)]
struct TlsSection {
    geometry: TlsGeom,
    tolerances: TlsTol,
}

#[test]
fn frozen_scqubits_fluxonium_and_transmon_hamiltonian_vv20() {
    let ds: ScqTlsDataset = load_validation_json("cqed_scqubits_tls.json");

    // 1. Validate Fluxonium Craig-Bampton quantization against frozen external scqubits JSON at Phi_ext = 0, 0.25, 0.5
    let flx = &ds.frozen_scqubits_fluxonium;
    let c_sigma = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * flx.ec_hz);
    let make_flux_rom = |phi_ext_rad: f64| -> SolveResult {
        let scale = jj_flux_scale();
        let w_r = 2.0 * PI * flx.f_r_hz;
        let c_cav = 420.0e-15;
        let inv_l_cav = w_r * w_r * c_cav;
        // 2-mode CircuitRom: decoupled Fluxonium JJ port DOF + spectator 7.0 GHz resonator DOF
        let mut rom = CircuitRom {
            ports: vec![RomPort {
                name: "flx_jj".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / flx.ej_hz),
                e_j_hz: Some(flx.ej_hz),
                c_j_f: 0.0,
                e_l_hz: Some(flx.el_hz),
                phi_ext_rad,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            }],
            r_static: 0,
            r_dyn: 2,
            k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * inv_l_cav]],
            c_r: vec![vec![0.0, 0.0], vec![0.0, 0.0]],
            m_r: vec![vec![c_sigma / EPS0, 0.0], vec![0.0, c_cav / EPS0]],
            b_r: vec![vec![1.0], vec![0.0]],
            c_port_f: vec![],
            modes: vec![],
        };
        rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
        SolveResult { physics: "ems".into(), study: "cavity".into(), scalars: vec![rom.scalar()], ..Default::default() }
    };

    for fp in &flx.flux_points {
        let phi_rad = 2.0 * PI * fp.phi_ext_over_phi0;
        let res = make_flux_rom(phi_rad);
        // Basis cutoff refinement (`n_flux_basis` = 32 vs 48) yields genuine numerical truncation error u_num
        let rep_c = hamiltonian_with_params(
            &res,
            &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [4], "n_flux_basis": 32, "phi_ext_rad": phi_rad }),
        )
        .unwrap();
        let rep_f = hamiltonian_with_params(
            &res,
            &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [4], "n_flux_basis": 48, "phi_ext_rad": phi_rad }),
        )
        .unwrap();
        for tgt in &fp.targets {
            let (s_c, s_f) = match tgt.quantity.as_str() {
                "f_01" => (rep_c.f01_hz[0], rep_f.f01_hz[0]),
                "f_02" => (rep_c.f02_hz[0], rep_f.f02_hz[0]),
                "f_12" => (rep_c.f12_hz[0], rep_f.f12_hz[0]),
                other => panic!("unexpected Fluxonium quantity '{other}'"),
            };
            let u_num = (s_f - s_c).abs().max(s_f * 1e-6);
            verify_asme_vv20(tgt, s_f, u_num);
        }
    }

    // 2. Validate Cooper-pair-box / Transmon Craig-Bampton against frozen scqubits charge-basis eigenvalues (< 0.05%)
    let tr = &ds.frozen_scqubits_transmon;
    let c_tr = (E_CHARGE * E_CHARGE) / (2.0 * PLANCK_H * tr.ec_hz);
    let scale = jj_flux_scale();
    let w_cav_tr = 2.0 * PI * 7.0e9;
    let c_cav_tr = 420.0e-15;
    let mut rom_tr = CircuitRom {
        ports: vec![RomPort {
            name: "tr_jj".into(),
            kind: "lumped_jj".into(),
            l_j_h: Some(scale / tr.ej_hz),
            e_j_hz: Some(tr.ej_hz),
            c_j_f: 0.0,
            e_l_hz: None,
            phi_ext_rad: 0.0,
            n_g: 0.0,
            squid_asymmetry: None,
            z0_ohm: None,
            direction: [0.0, 1.0, 0.0],
            formulation: "rank1".into(),
        }],
        r_static: 0,
        r_dyn: 2,
        k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * w_cav_tr * w_cav_tr * c_cav_tr]],
        c_r: vec![vec![0.0, 0.0], vec![0.0, 0.0]],
        m_r: vec![vec![c_tr / EPS0, 0.0], vec![0.0, c_cav_tr / EPS0]],
        b_r: vec![vec![1.0], vec![0.0]],
        c_port_f: vec![],
        modes: vec![],
    };
    rom_tr.modes = project_rom_modes(&rom_tr, &rom_tr.ports).unwrap();
    let res_tr = SolveResult {
        physics: "ems".into(),
        study: "cavity".into(),
        scalars: vec![rom_tr.scalar()],
        ..Default::default()
    };
    let tr_c =
        hamiltonian_with_params(&res_tr, &json!({ "basis": "craig_bampton", "port_trunc": 8, "fock_trunc": [4] }))
            .unwrap();
    let tr_f =
        hamiltonian_with_params(&res_tr, &json!({ "basis": "craig_bampton", "port_trunc": 14, "fock_trunc": [4] }))
            .unwrap();
    for tgt in &tr.targets {
        let (s_c, s_f) = match tgt.quantity.as_str() {
            "f_01" => (tr_c.f01_hz[0], tr_f.f01_hz[0]),
            "alpha" => (tr_c.alpha_hz[0], tr_f.alpha_hz[0]),
            other => panic!("unexpected Transmon quantity '{other}'"),
        };
        let u_num = (s_f - s_c).abs().max(s_f.abs() * 1e-6);
        verify_asme_vv20(tgt, s_f, u_num);
    }
}

fn solve_cpw_2d_explicit(g: &TlsGeom, h_global: f64, h_corner: f64, eps_sub: f64) -> (MeshedGeom, Mesh, SolveResult) {
    let half_box = 40.0e-6;
    let h_sub = 40.0e-6;
    let h_air = 40.0e-6;
    let a = g.w_trace_m / 2.0;
    let b = a + g.g_gap_m;

    let mut m = Model::new_2d();
    let sub = m.rect("substrate", (-half_box, -h_sub), (half_box, 0.0));
    let vac = m.rect("vacuum", (-half_box, 0.0), (half_box, h_air));
    let trace_l = m.rect("trace_l", (-a, 0.0), (-a + 1.2e-6, g.d_film_m));
    let trace_m = m.rect("trace_m", (-a + 1.2e-6, 0.0), (a - 1.2e-6, g.d_film_m));
    let trace_r = m.rect("trace_r", (a - 1.2e-6, 0.0), (a, g.d_film_m));
    let gnd_l_far = m.rect("gnd_l_far", (-half_box, 0.0), (-b - 1.2e-6, g.d_film_m));
    let gnd_l_tip = m.rect("gnd_l_tip", (-b - 1.2e-6, 0.0), (-b, g.d_film_m));
    let gnd_r_tip = m.rect("gnd_r_tip", (b, 0.0), (b + 1.2e-6, g.d_film_m));
    let gnd_r_far = m.rect("gnd_r_far", (b + 1.2e-6, 0.0), (half_box, g.d_film_m));

    m.overlay_full(
        &[
            ("substrate", sub, RegionKind::Material, None),
            ("vacuum", vac, RegionKind::Material, None),
            ("trace_l", trace_l, RegionKind::Hole, None),
            ("trace_m", trace_m, RegionKind::Hole, None),
            ("trace_r", trace_r, RegionKind::Hole, None),
            ("gnd_l_far", gnd_l_far, RegionKind::Hole, None),
            ("gnd_l_tip", gnd_l_tip, RegionKind::Hole, None),
            ("gnd_r_tip", gnd_r_tip, RegionKind::Hole, None),
            ("gnd_r_far", gnd_r_far, RegionKind::Hole, None),
        ],
        &[],
    );
    for bnd in ["trace_l.left", "trace_r.right", "gnd_l_tip.right", "gnd_r_tip.left"] {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(h_corner)),
            None,
            None,
        );
    }
    for bnd in [
        "trace_l.top",
        "trace_l.bottom",
        "trace_r.top",
        "trace_r.bottom",
        "gnd_l_tip.top",
        "gnd_l_tip.bottom",
        "gnd_r_tip.top",
        "gnd_r_tip.bottom",
    ] {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(h_corner * 1.5)),
            None,
            None,
        );
    }
    for bnd in ["trace_m.top", "trace_m.bottom"] {
        m.size(
            &format!("sz_{}", bnd.replace('.', "_")),
            Selector::Name(bnd.into()),
            Some(Expr::Num(h_corner * 6.0)),
            None,
            None,
        );
    }

    let gm = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h_global)).unwrap();
    let sub_dom = gm.domain("substrate").unwrap();
    let vac_dom = gm.domain("vacuum").unwrap();
    let tol = 1e-8;
    let mut trace_labs = Vec::new();
    let mut gnd_labs = Vec::new();
    for f in 0..gm.mesh.n_facets() {
        let fv = gm.mesh.facet(f);
        let p0 = gm.mesh.point(fv[0] as usize);
        let p1 = gm.mesh.point(fv[1] as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let zm = 0.5 * (p0[1] + p1[1]);
        let lab = gm.mesh.facet_labels[f];
        if zm >= -tol && zm <= g.d_film_m + tol {
            if xm.abs() <= a + tol && !trace_labs.contains(&lab) {
                trace_labs.push(lab);
            } else if xm.abs() >= b - tol && !gnd_labs.contains(&lab) {
                gnd_labs.push(lab);
            }
        }
    }

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", eps_sub);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());
    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&gm.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&gm.mesh).unwrap();
    (gm, ex_mesh, res)
}

fn solve_cpw_3d_sheet(g: &TlsGeom, h_global: f64, h_edge: f64, eps_sub: f64) -> (MeshedGeom, Mesh, SolveResult) {
    let half_box = 40.0e-6;
    let h_sub = 40.0e-6;
    let h_air = 40.0e-6;
    let l_y = 12.0e-6;
    let a = g.w_trace_m / 2.0;
    let b = a + g.g_gap_m;
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];

    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-half_box, 0.0, -h_sub), (half_box, l_y, 0.0));
    let vac = m.box_("vacuum", (-half_box, 0.0, 0.0), (half_box, l_y, h_air));
    let trace_p = m.polygons("trace_p", vec![vec![pt(-a, 0.0), pt(a, 0.0), pt(a, l_y), pt(-a, l_y)]], vec![]);
    let gnd_l_p =
        m.polygons("gnd_l_p", vec![vec![pt(-half_box, 0.0), pt(-b, 0.0), pt(-b, l_y), pt(-half_box, l_y)]], vec![]);
    let gnd_r_p =
        m.polygons("gnd_r_p", vec![vec![pt(b, 0.0), pt(half_box, 0.0), pt(half_box, l_y), pt(b, l_y)]], vec![]);
    m.sheet("trace", trace_p, 0.0);
    m.sheet("gnd_l", gnd_l_p, 0.0);
    m.sheet("gnd_r", gnd_r_p, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for edge_sel in ["trace:edge", "gnd_l:edge", "gnd_r:edge"] {
        m.size(
            &format!("sz_{}", edge_sel.replace(':', "_")),
            Selector::parse(edge_sel).unwrap(),
            Some(Expr::Num(h_edge)),
            None,
            None,
        );
    }
    m.size("sz_trace_face", Selector::Name("trace".into()), Some(Expr::Num(h_edge * 1.8)), None, None);
    let gm = mesh_model(&m.list, &MeshPlan::simple(h_global)).unwrap();
    let sub_dom = gm.domain("substrate").unwrap();
    let vac_dom = gm.domain("vacuum").unwrap();
    let trace_labs = gm.boundary("trace");
    let mut gnd_labs = gm.boundary("gnd_l");
    gnd_labs.extend(gm.boundary("gnd_r"));

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", eps_sub);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());
    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&gm.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&gm.mesh).unwrap();
    (gm, ex_mesh, res)
}

#[test]
fn cryogenic_cpw_q_tls_3grid_gci_and_calusine_trench_vv20() {
    let ds: ScqTlsDataset = load_validation_json("cqed_scqubits_tls.json");
    let g = &ds.cryogenic_q_tls_cpw.geometry;
    let tol = &ds.cryogenic_q_tls_cpw.tolerances;

    let base_params = json!({
        "solution": 0,
        "t_ma": g.t_ox_m, "t_ms": g.t_ox_m, "t_sa": g.t_ox_m,
        "eps_ma": g.eps_ox, "eps_ms": g.eps_si, "eps_sa": g.eps_ox,
        "tan_delta_ma": g.tan_delta, "tan_delta_ms": g.tan_delta, "tan_delta_sa": g.tan_delta,
        "d": g.d_film_m,
        "edge_cutoff": g.edge_cutoff_m
    });

    // 1. Resolved 2D finite-thickness (d = 100 nm) CPW reference (cross-checked via 2 nm volume integral formula)
    let (_, m2_fine, r2_fine) = solve_cpw_2d_explicit(g, 3.5e-6, 35.0e-9, g.eps_si);
    let rep2 = tls_participation_with_params(&m2_fine, &r2_fine, &base_params).unwrap();
    let _ = (facet_rule, geometry, quad_rule, shape, EPS0, scalar_exact);

    // 2. 3D matched-asymptotic CPW across 3 grids -> genuine Roache GCI u_num
    let mut reps3 = Vec::new();
    let mut cells3 = Vec::new();
    for &(hg, he) in &[(12.0e-6, 1.6e-6), (10.0e-6, 1.1e-6), (8.0e-6, 0.8e-6)] {
        let (gm3, m3, r3) = solve_cpw_3d_sheet(g, hg, he, g.eps_si);
        cells3.push(gm3.mesh.n_cells() as f64);
        reps3.push(tls_participation_with_params(&m3, &r3, &base_params).unwrap());
    }
    let r21 = (cells3[2] / cells3[1]).cbrt().max(1.1);
    let r32 = (cells3[1] / cells3[0]).cbrt().max(1.1);

    for (name, v_fine, v_mid, v_coarse, v_ref) in [
        ("p_ma", reps3[2].p_ma, reps3[1].p_ma, reps3[0].p_ma, rep2.p_ma),
        ("p_ms", reps3[2].p_ms, reps3[1].p_ms, reps3[0].p_ms, rep2.p_ms),
        ("p_sa", reps3[2].p_sa, reps3[1].p_sa, reps3[0].p_sa, rep2.p_sa),
        ("Q_TLS", reps3[2].q_tls, reps3[1].q_tls, reps3[0].q_tls, rep2.q_tls),
    ] {
        let (u_num, _) = u_num_three_grid(v_fine, v_mid, v_coarse, r21, r32, 1.5);
        let target = common::inline_target(
            format!("cpw_3d_tls_{name}"),
            format!("3D matched-asymptotic CPW {name} vs 2D resolved d=100nm reference"),
            name,
            if name == "Q_TLS" { "1" } else { "ratio" },
            v_ref,
            v_ref * tol.u_input_rel,
            v_ref * tol.u_d_rel,
            tol.interface_rel,
            0.0,
        );
        verify_asme_vv20(&target, v_fine, u_num);
    }

    // 3. Calusine (2018) trench-depth scaling on finest 3D mesh: monotonic Q_TLS improvement
    let (_, m3_fine, r3_fine) = solve_cpw_3d_sheet(g, 8.0e-6, 0.8e-6, g.eps_si);
    let mut prev_q = 0.0_f64;
    for &td in &tol.trench_depths_m {
        let mut p = base_params.clone();
        p["trench_depth_m"] = json!(td);
        let rep = tls_participation_with_params(&m3_fine, &r3_fine, &p).unwrap();
        assert!(
            rep.q_tls > prev_q,
            "Calusine trench {td:e} m must monotonically raise Q_TLS: {} vs {prev_q}",
            rep.q_tls
        );
        prev_q = rep.q_tls;
    }

    // 4. Place (2021) sapphire vs silicon dielectric permittivity check with fixed interface oxide dielectric:
    //    lowering eps_sub from Si (11.45) to sapphire (10.2) reduces MS normal-displacement factor eps_sub^2/eps_ms -> strictly lower p_ms
    let mut p_si_ox = base_params.clone();
    p_si_ox["eps_ms"] = json!(g.eps_ox);
    let rep_si_ox = tls_participation_with_params(&m2_fine, &r2_fine, &p_si_ox).unwrap();

    let (_, m_sap, r_sap) = solve_cpw_2d_explicit(g, 3.5e-6, 35.0e-9, g.eps_sapphire_place);
    let mut p_sap = p_si_ox;
    let rep_sap = tls_participation_with_params(&m_sap, &r_sap, &p_sap).unwrap();
    assert!(
        rep_sap.p_ms < 0.92 * rep_si_ox.p_ms,
        "Sapphire (eps=10.2) must lower MS participation below Si (eps=11.45): sap={}, si={}",
        rep_sap.p_ms,
        rep_si_ox.p_ms
    );
    // Clean annealed Ta-on-sapphire surface (Place 2021: reduced native SA/MS loss tangent 2e-4) raises Q_TLS > 2x
    p_sap["tan_delta_ms"] = json!(2.0e-4);
    p_sap["tan_delta_sa"] = json!(2.0e-4);
    let rep_sap_clean = tls_participation_with_params(&m_sap, &r_sap, &p_sap).unwrap();
    assert!(
        rep_sap_clean.q_tls > 2.0 * rep_si_ox.q_tls,
        "Place (2021) low-loss sapphire surfaces must raise Q_TLS > 2x over native Si: sap={}, si={}",
        rep_sap_clean.q_tls,
        rep_si_ox.q_tls
    );
}
