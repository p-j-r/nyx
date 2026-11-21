// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! AF-6 integration tests for 3D thick-sheet (`thickness > 0`) and etched-trench (`etch_depth > 0`)
//! TLS surface participation classification (`MA`, `MS`, `SA`) in `nyx-analysis`.

use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::Model;
use geom::select::Selector;
use nyx_analysis::{tls_participation_with_params, TlsReport};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study};
use nyx_physics::Mesh;
use serde_json::json;

const W_TRACE: f64 = 10.0e-6;
const G_GAP: f64 = 6.0e-6;
const EPS_SI: f64 = 11.45;
const EPS_OX: f64 = 10.0;
const T_OX: f64 = 2.0e-9;
const HALF_BOX: f64 = 35.0e-6;
const H_SUB: f64 = 35.0e-6;
const H_AIR: f64 = 35.0e-6;

fn solve_thick_cpw_3d(t_metal: f64, t_etch: f64, h_global: f64, h_edge: f64) -> TlsReport {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let l_y = 10.0e-6;
    let a = W_TRACE / 2.0;
    let b = a + G_GAP;

    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-HALF_BOX, 0.0, -H_SUB), (HALF_BOX, l_y, 0.0));
    let vac = m.box_("vacuum", (-HALF_BOX, 0.0, 0.0), (HALF_BOX, l_y, H_AIR));

    let trace_p = m.polygons("trace_p", vec![vec![pt(-a, 0.0), pt(a, 0.0), pt(a, l_y), pt(-a, l_y)]], vec![]);
    let gnd_l_p =
        m.polygons("gnd_l_p", vec![vec![pt(-HALF_BOX, 0.0), pt(-b, 0.0), pt(-b, l_y), pt(-HALF_BOX, l_y)]], vec![]);
    let gnd_r_p =
        m.polygons("gnd_r_p", vec![vec![pt(b, 0.0), pt(HALF_BOX, 0.0), pt(HALF_BOX, l_y), pt(b, l_y)]], vec![]);

    m.sheet_thick("trace", trace_p, 0.0, t_metal, t_etch);
    m.sheet_thick("gnd_l", gnd_l_p, 0.0, t_metal, t_etch);
    m.sheet_thick("gnd_r", gnd_r_p, 0.0, t_metal, t_etch);
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
    m.size("sz_trace_face", Selector::Name("trace".into()), Some(Expr::Num(h_edge * 2.0)), None, None);

    let g = mesh_model(&m.list, &MeshPlan::simple(h_global)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    if t_metal > 0.0 {
        for lab in g.boundary("trace") {
            assert_eq!(g.labels.sheet_thickness_m.get(&lab).copied(), Some(t_metal));
        }
    }
    if t_etch > 0.0 {
        for lab in g.boundary("trace") {
            assert_eq!(g.labels.sheet_etch_depth_m.get(&lab).copied(), Some(t_etch));
        }
    }

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", EPS_SI);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());

    let trace_labs = g.boundary("trace");
    let mut gnd_labs = g.boundary("gnd_l");
    gnd_labs.extend(g.boundary("gnd_r"));

    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let params = json!({
        "solution": 0,
        "t_ma": T_OX,
        "t_ms": T_OX,
        "t_sa": T_OX,
        "eps_ma": EPS_OX,
        "eps_ms": EPS_SI,
        "eps_sa": EPS_OX,
        "tan_delta_ma": 1.0e-3,
        "tan_delta_ms": 1.0e-3,
        "tan_delta_sa": 1.0e-3,
    });
    tls_participation_with_params(&ex_mesh, &res, &params).unwrap()
}

#[test]
fn thick_sheet_and_trench_3d_tls_classification_disables_knife_edge_double_counting() {
    let rep_base = solve_thick_cpw_3d(150.0e-9, 0.0, 9.0e-6, 1.5e-6);
    assert!(rep_base.p_ma > 0.0, "expected positive MA participation on thick 3D conductor");
    assert!(rep_base.p_ms > 0.0, "expected positive MS participation on thick 3D conductor");
    assert!(rep_base.p_sa > 0.0, "expected positive SA participation on uncovered substrate");
    assert!(
        (rep_base.p_total - (rep_base.p_ma + rep_base.p_ms + rep_base.p_sa)).abs() < 1e-15,
        "p_total must equal p_ma + p_ms + p_sa"
    );
    // Resolved 3D thick sheet must NOT apply the 1D singular knife-edge tube correction
    assert_eq!(rep_base.p_ma_edge, 0.0);
    assert_eq!(rep_base.p_ms_edge, 0.0);
    assert_eq!(rep_base.p_sa_edge, 0.0);
    assert_eq!(rep_base.edge_cutoff_m, 0.0);
    assert_eq!(rep_base.convergence_check["perimeter_segments"].as_u64(), Some(0));
    assert!((rep_base.convergence_check["d_film_m"].as_f64().unwrap() - 150.0e-9).abs() < 1e-15);

    // Adding an 80 nm substrate trench etch exposes vertical SA trench sidewalls and increases p_sa
    let rep_trench = solve_thick_cpw_3d(150.0e-9, 80.0e-9, 9.0e-6, 1.5e-6);
    assert_eq!(rep_trench.p_ma_edge, 0.0);
    assert_eq!(rep_trench.p_ms_edge, 0.0);
    assert_eq!(rep_trench.p_sa_edge, 0.0);
    assert!(
        rep_trench.p_sa > 0.0 && rep_trench.p_ma > 0.0 && rep_trench.p_ms > 0.0,
        "trench CPW must classify MA, MS, and SA (including vertical trench sidewalls)"
    );
}
