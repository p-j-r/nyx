// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative Acceptance Benchmark Test 1 (Multi-RHS Capacitance Matrix) and Milestone 2.2
//! Gauge-Stabilized London Magnetostatics verification tests (`cross_section`, `planar_stream`, `mixed_3d`).

use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::mesher::{extrude_y, Layer};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_dfs::{solve_electrostatics, solve_magnetostatics, Bc, BcKind, Dfs, Problem, Study, EPS0, MU0};
use nyx_physics::{
    CapacitanceMatrixData, Ctx, InductanceMatrixData, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec,
    SolveResult, TlsMeta,
};
use serde_json::json;
use std::f64::consts::PI;

/// Builds a 3D coaxial cylinder of inner radius `a`, outer radius `b`, and length `l` along `y`,
/// with curved quadratic (`P2`) mid-edge geometry nodes projected onto the exact circular surfaces `r = a` and `r = b`.
fn build_curved_3d_coax_mesh(
    a: f64,
    b: f64,
    l: f64,
    h2d: f64,
    n_layers: usize,
) -> (femcore::mesh::SimplexMesh, Vec<i32>, Vec<i32>) {
    let mut m = Model::new_2d();
    let outer = m.circle("diel", (0.0, 0.0), b);
    let inner = m.circle("core", (0.0, 0.0), a);
    m.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let mut opts = GeomMeshOpts::new(h2d);
    opts.max_dev = Some(h2d / 200.0);
    let g2 = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (a + b);
    let mut inner_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        if rm < r_mid {
            if !inner_labs.contains(&lab) {
                inner_labs.push(lab);
            }
        } else if !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }
    assert!(!inner_labs.is_empty() && !outer_labs.is_empty());

    let layers = vec![Layer { thickness: l, n: n_layers, domain_offset: 0, fill: None }];
    let mut m3 = extrude_y(&g2.mesh, &layers, 0.0).unwrap();

    // Curve all edges whose endpoints both lie on r = a or both lie on r = b
    let tol_a = 1e-3 * a;
    let tol_b = 1e-3 * b;
    let mut edge_nodes = Vec::new();
    for &(u, v) in m3.edges() {
        let pu = m3.point(u as usize);
        let pv = m3.point(v as usize);
        let ru = (pu[0] * pu[0] + pu[2] * pu[2]).sqrt();
        let rv = (pv[0] * pv[0] + pv[2] * pv[2]).sqrt();
        let target_r = if (ru - a).abs() < tol_a && (rv - a).abs() < tol_a {
            Some(a)
        } else if (ru - b).abs() < tol_b && (rv - b).abs() < tol_b {
            Some(b)
        } else {
            None
        };
        if let Some(r_exact) = target_r {
            let mx = 0.5 * (pu[0] + pv[0]);
            let my = 0.5 * (pu[1] + pv[1]);
            let mz = 0.5 * (pu[2] + pv[2]);
            let mr = (mx * mx + mz * mz).sqrt();
            if mr > 1e-12 && (mr - r_exact).abs() > 1e-12 * r_exact {
                edge_nodes.push(([u, v], [mx * (r_exact / mr), my, mz * (r_exact / mr)]));
            }
        }
    }
    m3.set_edge_nodes(&edge_nodes);
    assert!(m3.has_curved(), "3D coaxial mesh must carry curved mid-edge nodes");
    (m3, inner_labs, outer_labs)
}

#[test]
fn test_1a_3d_coaxial_cylinder_analytic_capacitance_check() {
    // Test 1a: a = 10 um, b = 30 um, L = 100 um, eps_r = 11.7
    // Exact analytic capacitance: C = 2 pi eps_0 eps_r L / ln(b / a)
    let a: f64 = 10.0e-6;
    let b: f64 = 30.0e-6;
    let l: f64 = 100.0e-6;
    let eps_r: f64 = 11.7;
    let c_exact = (2.0 * PI * EPS0 * eps_r * l) / (b / a).ln();

    let (mesh3, inner_labs, outer_labs) = build_curved_3d_coax_mesh(a, b, l, 4.0e-6, 2);

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.default.insert("eps_r".into(), eps_r);
    prob.bcs = vec![
        Bc { labels: inner_labs, kind: BcKind::Terminal("inner".into()) },
        Bc { labels: outer_labs, kind: BcKind::Ground },
    ];

    let res = solve_electrostatics(&mesh3, &prob).unwrap();
    assert!(
        res.log.iter().any(|line| line.contains("1 factorization")),
        "expected log confirming 1 factorization, got {:?}",
        res.log
    );

    let cap = CapacitanceMatrixData::of(&res).expect("CapacitanceMatrixData must be present");
    assert_eq!(cap.terminals, vec!["inner"]);
    let c_fem = cap.maxwell_matrix_f[0][0];
    let rel_err = (c_fem - c_exact).abs() / c_exact;
    eprintln!(
        "Test 1a (3D curved P2 coaxial cylinder): C_fem = {c_fem:.8e} F, C_exact = {c_exact:.8e} F, rel_err = {:.5}%",
        rel_err * 100.0
    );
    assert!(
        rel_err < 1.0e-3,
        "Test 1a failed: relative error {:.4}% exceeds 0.1% threshold (C_fem = {c_fem:e}, C_exact = {c_exact:e})",
        rel_err * 100.0
    );
    assert!(TlsMeta::of(&res).is_some(), "TlsMeta must be attached to electrostatic result");
}

fn build_planar_transmon_coupler_model(pad_edge_h: f64) -> Model {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    // 600 um x 600 um silicon + vacuum box around z = 0
    let sub = m.box_("substrate", (-300e-6, -300e-6, -150e-6), (300e-6, 300e-6, 0.0));
    let vac = m.box_("vacuum", (-300e-6, -300e-6, 0.0), (300e-6, 300e-6, 150e-6));

    // 2-pad transmon + coupler at z = 0:
    // pad_1: x in [-160, -60] um, y in [-50, 50] um
    let p1_poly = m.polygons(
        "p1_poly",
        vec![vec![pt(-160e-6, -50e-6), pt(-60e-6, -50e-6), pt(-60e-6, 50e-6), pt(-160e-6, 50e-6)]],
        vec![],
    );
    // pad_2: x in [-20, 80] um, y in [-50, 50] um (40 um gap to pad_1)
    let p2_poly = m.polygons(
        "p2_poly",
        vec![vec![pt(-20e-6, -50e-6), pt(80e-6, -50e-6), pt(80e-6, 50e-6), pt(-20e-6, 50e-6)]],
        vec![],
    );
    // coupler: x in [120, 180] um, y in [-50, 50] um (40 um gap to pad_2)
    let c_poly = m.polygons(
        "c_poly",
        vec![vec![pt(120e-6, -50e-6), pt(180e-6, -50e-6), pt(180e-6, 50e-6), pt(120e-6, 50e-6)]],
        vec![],
    );

    m.sheet("pad_1", p1_poly, 0.0);
    m.sheet("pad_2", p2_poly, 0.0);
    m.sheet("coupler", c_poly, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    for s in ["pad_1:edge", "pad_2:edge", "coupler:edge"] {
        m.size(
            &format!("sz_{}", s.replace(':', "_")),
            Selector::parse(s).unwrap(),
            Some(Expr::Num(pad_edge_h)),
            None,
            None,
        );
    }
    for s in ["pad_1", "pad_2", "coupler"] {
        m.size(&format!("sz_face_{s}"), Selector::Name(s.into()), Some(Expr::Num(pad_edge_h * 1.5)), None, None);
    }
    m
}

fn solve_transmon_3port(global_h: f64, pad_edge_h: f64) -> (CapacitanceMatrixData, SolveResult, usize) {
    let m = build_planar_transmon_coupler_model(pad_edge_h);
    let g = mesh_model(&m.list, &MeshPlan::simple(global_h)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", 11.45);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());

    let p1_labs = g.boundary("pad_1");
    let p2_labs = g.boundary("pad_2");
    let c_labs = g.boundary("coupler");
    assert!(!p1_labs.is_empty() && !p2_labs.is_empty() && !c_labs.is_empty());

    // Exterior box walls as ground reference
    let sheet_set: std::collections::BTreeSet<i32> = p1_labs.iter().chain(&p2_labs).chain(&c_labs).copied().collect();
    let ext = g.mesh.exterior_facets();
    let mut ground_labs = Vec::new();
    for (f, &is_ext) in ext.iter().enumerate() {
        let lab = g.mesh.facet_labels[f];
        if is_ext && !sheet_set.contains(&lab) && !ground_labs.contains(&lab) {
            ground_labs.push(lab);
        }
    }

    prob.bcs = vec![
        Bc { labels: p1_labs, kind: BcKind::Terminal("pad_1".into()) },
        Bc { labels: p2_labs, kind: BcKind::Terminal("pad_2".into()) },
        Bc { labels: c_labs, kind: BcKind::Terminal("coupler".into()) },
        Bc { labels: ground_labs, kind: BcKind::Ground },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let cap = CapacitanceMatrixData::of(&res).unwrap();
    (cap, res, g.mesh.n_cells())
}

#[test]
fn test_1b_3d_zero_thickness_planar_transmon_coupler_convergence_and_symmetry() {
    let (cap_coarse, _res_coarse, cells_coarse) = solve_transmon_3port(70e-6, 11e-6);
    let (cap_fine, res_fine, cells_fine) = solve_transmon_3port(60e-6, 9e-6);

    assert!(cells_fine > cells_coarse);
    assert!(
        res_fine
            .log
            .iter()
            .any(|l| (l.contains("1 factorization") || l.contains("1 hierarchy setup")) && l.contains("3 terminal")),
        "expected 1 factorization/setup for 3 terminals in log: {:?}",
        res_fine.log
    );
    assert_eq!(res_fine.solutions.len(), 3);
    assert_eq!(cap_fine.terminals, vec!["pad_1", "pad_2", "coupler"]);

    // Verify 3x3 symmetry (< 1e-12) and physical sign properties
    for i in 0..3 {
        assert!(cap_fine.maxwell_matrix_f[i][i] > 0.0);
        assert!(cap_fine.mutual_matrix_f[i][i] > 0.0);
        assert!(cap_fine.ec_matrix_hz[i][i] > 0.0);
        for j in 0..3 {
            let sym_c = (cap_fine.maxwell_matrix_f[i][j] - cap_fine.maxwell_matrix_f[j][i]).abs()
                / cap_fine.maxwell_matrix_f[i][i];
            let sym_ec =
                (cap_fine.ec_matrix_hz[i][j] - cap_fine.ec_matrix_hz[j][i]).abs() / cap_fine.ec_matrix_hz[i][i];
            assert!(sym_c < 1e-12, "C_M symmetry violation at ({i},{j}): {sym_c:e}");
            assert!(sym_ec < 1e-12, "E_C symmetry violation at ({i},{j}): {sym_ec:e}");
            if i != j {
                assert!(cap_fine.maxwell_matrix_f[i][j] < 0.0, "off-diagonal Maxwell C_{i}{j} must be negative");
                assert!(cap_fine.mutual_matrix_f[i][j] > 0.0, "off-diagonal mutual c_{i}{j} must be positive");
            }
        }
    }

    // Verify mesh-refinement convergence < 1.0% on the self/total capacitance norm
    let mut max_rel_diff = 0.0_f64;
    for i in 0..3 {
        let c_c = cap_coarse.maxwell_matrix_f[i][i];
        let c_f = cap_fine.maxwell_matrix_f[i][i];
        let rel = (c_f - c_c).abs() / c_f;
        max_rel_diff = max_rel_diff.max(rel);
    }
    eprintln!(
        "Test 1b (3D planar 2-pad transmon + coupler): coarse ({cells_coarse} cells) vs fine ({cells_fine} cells) max diagonal C_ii shift = {:.4}%",
        max_rel_diff * 100.0
    );
    assert!(
        max_rel_diff < 0.01,
        "Test 1b mesh-refinement convergence {:.4}% exceeds 1.0% threshold",
        max_rel_diff * 100.0
    );

    // Verify TlsMeta captures interior sheet labels and domain permittivities
    let tls = TlsMeta::of(&res_fine).expect("TlsMeta must be present");
    assert_eq!(tls.domain_eps_r.len(), 2);
    assert!(!tls.sheet_labels.is_empty(), "interior transmon sheets must be recorded in TlsMeta::sheet_labels");
}

#[test]
fn test_magnetostatics_2d_cross_section_coax_with_and_without_kinetic_inductance() {
    // 2D coaxial cross-section: a = 10 um, b = 30 um, mu_r = 1.0
    // Analytic geometric inductance per unit length: L'_geom = (mu_0 / 2 pi) ln(b / a)
    // Analytic kinetic inductance per unit length on inner conductor: L'_kin = L_k,sq / (2 pi a)
    let a: f64 = 10.0e-6;
    let b: f64 = 30.0e-6;
    let mut m = Model::new_2d();
    let outer = m.circle("diel", (0.0, 0.0), b);
    let inner = m.circle("core", (0.0, 0.0), a);
    m.overlay_full(&[("diel", outer, RegionKind::Material, None), ("core", inner, RegionKind::Hole, None)], &[]);
    let mut opts = GeomMeshOpts::new(2.5e-6);
    opts.max_dev = Some(2.5e-6 / 100.0);
    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &opts).unwrap();
    let r_mid = 0.5 * (a + b);
    let mut inner_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        let p0 = g.mesh.point(fv[0] as usize);
        let p1 = g.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g.mesh.facet_labels[f];
        if rm < r_mid {
            if !inner_labs.contains(&lab) {
                inner_labs.push(lab);
            }
        } else if !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }

    let l_geom_exact = (MU0 / (2.0 * PI)) * (b / a).ln();

    // 1. Pure geometric inductance (no london_sheet)
    let mut p_geom = Problem::new(Study::Magnetostatic { formulation: "cross_section".into() }, 2);
    p_geom.bcs = vec![
        Bc {
            labels: inner_labs.clone(),
            kind: BcKind::CurrentTerminal { name: "center".into(), current: 1.0, direction: [0.0, 1.0, 0.0] },
        },
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
    ];
    let r_geom = solve_magnetostatics(&g.mesh, &p_geom).unwrap();
    let ind_geom = InductanceMatrixData::of(&r_geom).unwrap();
    assert_eq!(ind_geom.formulation, "cross_section");
    let l_fem_geom = ind_geom.inductance_matrix_h[0][0];
    let err_geom = (l_fem_geom - l_geom_exact).abs() / l_geom_exact;
    eprintln!(
        "2D cross_section coax geometric L': fem = {l_fem_geom:.6e} H/m, exact = {l_geom_exact:.6e} H/m, err = {:.4}%",
        err_geom * 100.0
    );
    assert!(err_geom < 5e-3, "2D coax geometric L' error {err_geom:e} too large");

    // 2. Geometric + kinetic inductance on inner conductor (L_k,sq = 1.5 pH/sq)
    let l_k_sq = 1.5e-12;
    let l_kin_exact = l_k_sq / (2.0 * PI * a);
    let mut p_tot = p_geom.clone();
    p_tot.bcs.push(Bc { labels: inner_labs, kind: BcKind::LondonSheet { l_k_sq } });
    let r_tot = solve_magnetostatics(&g.mesh, &p_tot).unwrap();
    let l_fem_tot = InductanceMatrixData::of(&r_tot).unwrap().inductance_matrix_h[0][0];
    let err_kin = ((l_fem_tot - l_fem_geom) - l_kin_exact).abs() / l_kin_exact;
    eprintln!(
        "2D cross_section coax kinetic delta L': fem = {:.6e} H/m, exact = {l_kin_exact:.6e} H/m, err = {:.4}%",
        l_fem_tot - l_fem_geom,
        err_kin * 100.0
    );
    assert!(err_kin < 5e-3, "2D coax kinetic L' error {err_kin:e} too large");
}

#[test]
fn test_magnetostatics_2d_planar_stream_strip_kinetic_and_geometric_inductance() {
    // Superconducting strip of length L = 40 um and width W = 10 um (N_sq = L / W = 4 squares):
    // Boundaries of femcore::testkit::rect(L, W): 1 is x=0, 2 is z=0, 3 is z=W, 4 is x=L.
    // Stream function g = 0 on boundary 2 (z = 0) and g = I on boundary 3 (z = W), natural dg/dn = 0 on ports 1 & 4.
    let (l_strip, w_strip) = (40.0e-6, 10.0e-6);
    let mesh = femcore::testkit::rect(l_strip, w_strip, 12, 4);

    // 1. Pure geometric Biot-Savart inductance (L_k,sq = 0)
    let mut p_bs = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_bs.bcs = vec![
        Bc { labels: vec![2], kind: BcKind::Ground },
        Bc {
            labels: vec![3],
            kind: BcKind::CurrentTerminal { name: "strip".into(), current: 1.0, direction: [1.0, 0.0, 0.0] },
        },
    ];
    let r_bs = solve_magnetostatics(&mesh, &p_bs).unwrap();
    let l_bs = InductanceMatrixData::of(&r_bs).unwrap().inductance_matrix_h[0][0];
    assert!(l_bs > 5.0e-12 && l_bs < 50.0e-12, "unexpected 2D Biot-Savart strip geometric inductance {l_bs:e} H");

    // 2. Add London sheet kinetic inductance L_k,sq = 2.0 pH/sq => L_kin = L_k,sq * (L / W) = 8.0 pH
    let l_k_sq = 2.0e-12;
    let l_kin_exact = l_k_sq * (l_strip / w_strip);
    let mut p_tot = p_bs.clone();
    p_tot.bcs.push(Bc { labels: vec![1, 2, 3, 4], kind: BcKind::LondonSheet { l_k_sq } });
    let r_tot = solve_magnetostatics(&mesh, &p_tot).unwrap();
    let l_tot = InductanceMatrixData::of(&r_tot).unwrap().inductance_matrix_h[0][0];
    assert!(
        l_tot > l_bs + 0.95 * l_kin_exact && l_tot < l_bs + 1.05 * l_kin_exact,
        "expected L_tot ({l_tot:e}) == L_bs ({l_bs:e}) + L_kin ({l_kin_exact:e})"
    );
}

#[test]
fn test_magnetostatics_3d_mixed_3d_superconducting_sheet_trace_and_bcs_material_decode() {
    // 3D superconducting trace embedded at z = 0 inside a 3D box using Op::Sheet
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("substrate", (-60e-6, -40e-6, -30e-6), (60e-6, 40e-6, 0.0));
    let vac = m.box_("vacuum", (-60e-6, -40e-6, 0.0), (60e-6, 40e-6, 30e-6));
    let tr_poly = m.polygons(
        "trace_poly",
        vec![vec![pt(-40e-6, -10e-6), pt(40e-6, -10e-6), pt(40e-6, 10e-6), pt(-40e-6, 10e-6)]],
        vec![],
    );
    m.sheet("trace", tr_poly, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(25e-6)).unwrap();
    let trace_labs = g.boundary("trace");
    assert!(!trace_labs.is_empty());

    // Solve via ProblemSpec + Library BCS material resolution (`nb` at 10 mK, d = 100 nm, film_model = "kerr")
    struct GResolver<'a>(&'a femcore::meshing::geomdrive::MeshedGeom);
    impl LabelResolver for GResolver<'_> {
        fn domains(&self, region: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            self.0.domain(region).map(|d| vec![d]).ok_or_else(|| nyx_physics::Error::invalid(region))
        }
        fn boundaries(&self, selector: &str) -> Result<Vec<i32>, nyx_physics::Error> {
            let v = self.0.boundary(selector);
            if v.is_empty() {
                Err(nyx_physics::Error::invalid(selector))
            } else {
                Ok(v)
            }
        }
    }

    let mut spec = ProblemSpec::new("dfs", 1, "magnetostatic");
    spec.study.params = json!({ "formulation": "mixed_3d" });
    spec.materials.insert("substrate".into(), MaterialAssign::Props([("mu_r".to_string(), 1.0)].into()));
    spec.materials.insert("vacuum".into(), MaterialAssign::Props([("mu_r".to_string(), 1.0)].into()));
    spec.bcs.push(nyx_physics::BcSpec {
        on: "trace".into(),
        kind: "london_sheet".into(),
        params: json!({
            "L_k_sq": 0.0,
            "material": "niobium",
            "thickness": 100e-9,
            "temperature": 0.01,
            "film_model": "kerr"
        }),
    });
    spec.bcs.push(nyx_physics::BcSpec {
        on: "trace".into(),
        kind: "current_terminal".into(),
        params: json!({
            "name": "trace_port",
            "current": 1.0,
            "direction": [1.0, 0.0, 0.0]
        }),
    });

    let resolver = GResolver(&g);
    let mut ctx = Ctx::new(&resolver);
    let mat_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../materials");
    ctx.library = Some(std::sync::Arc::new(femcore::material::Library::load_dir(&mat_dir).unwrap()));

    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let res = Dfs.solve(&ex_mesh, &spec, &mut ctx).unwrap();
    let ind = InductanceMatrixData::of(&res).expect("3D mixed_3d result must carry InductanceMatrixData");
    assert_eq!(ind.formulation, "mixed_3d");
    assert_eq!(ind.terminals, vec!["trace_port"]);
    let l_3d = ind.inductance_matrix_h[0][0];
    eprintln!("3D mixed_3d superconducting trace inductance L = {l_3d:.6e} H");
    assert!(l_3d > 0.0 && l_3d.is_finite());
    assert!(
        res.log.iter().any(|l| l.contains("1 factorization") && l.contains("deduplicated triangles")),
        "unexpected mixed_3d log: {:?}",
        res.log
    );
}

#[test]
fn test_magnetostatics_2d_planar_stream_and_3d_mixed_flux_loop() {
    // 1. 2D planar_stream superconducting ring/washer with inner hole flux_loop
    let r_in: f64 = 15.0e-6;
    let r_out: f64 = 35.0e-6;
    let mut m2 = Model::new_2d();
    let outer = m2.circle("washer", (0.0, 0.0), r_out);
    let inner = m2.circle("hole", (0.0, 0.0), r_in);
    m2.overlay_full(&[("washer", outer, RegionKind::Material, None), ("hole", inner, RegionKind::Hole, None)], &[]);
    let g2 = mesh_geom_2d(evaluate(&m2.list).unwrap(), &GeomMeshOpts::new(6.0e-6)).unwrap();
    let r_mid = 0.5 * (r_in + r_out);
    let mut hole_labs = Vec::new();
    let mut outer_labs = Vec::new();
    for f in 0..g2.mesh.n_facets() {
        let fv = g2.mesh.facet(f);
        let p0 = g2.mesh.point(fv[0] as usize);
        let p1 = g2.mesh.point(fv[1] as usize);
        let rm = (0.5 * (p0[0] + p1[0])).hypot(0.5 * (p0[1] + p1[1]));
        let lab = g2.mesh.facet_labels[f];
        if rm < r_mid {
            if !hole_labs.contains(&lab) {
                hole_labs.push(lab);
            }
        } else if !outer_labs.contains(&lab) {
            outer_labs.push(lab);
        }
    }

    let mut p_loop2 = Problem::new(Study::Magnetostatic { formulation: "planar_stream".into() }, 1);
    p_loop2.bcs = vec![
        Bc { labels: outer_labs.clone(), kind: BcKind::Ground },
        Bc {
            labels: hole_labs.clone(),
            kind: BcKind::FluxLoop { name: "squid_loop".into(), flux_wb: 2.067833848e-15 },
        },
        Bc { labels: outer_labs, kind: BcKind::LondonSheet { l_k_sq: 1.0e-12 } },
    ];
    let r_loop2 = solve_magnetostatics(&g2.mesh, &p_loop2).unwrap();
    let ind_loop2 = InductanceMatrixData::of(&r_loop2).unwrap();
    assert_eq!(ind_loop2.terminals, vec!["squid_loop"]);
    assert!(ind_loop2.inductance_matrix_h[0][0] > 0.0);

    // 2. 3D mixed_3d square superconducting washer loop at z = 0 driven by flux_loop
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m3 = Model::new_3d();
    let sub = m3.box_("substrate", (-60e-6, -60e-6, -30e-6), (60e-6, 60e-6, 0.0));
    let vac = m3.box_("vacuum", (-60e-6, -60e-6, 0.0), (60e-6, 60e-6, 30e-6));
    let outer_ring = vec![pt(-35e-6, -35e-6), pt(35e-6, -35e-6), pt(35e-6, 35e-6), pt(-35e-6, 35e-6)];
    let inner_hole = vec![pt(-15e-6, -15e-6), pt(15e-6, -15e-6), pt(15e-6, 15e-6), pt(-15e-6, 15e-6)];
    let washer_poly = m3.polygons("washer_poly", vec![outer_ring], vec![inner_hole]);
    m3.sheet("washer", washer_poly, 0.0);
    m3.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let g3 = mesh_model(&m3.list, &MeshPlan::simple(25e-6)).unwrap();
    let washer_labs = g3.boundary("washer");
    let mut p_loop3 = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
    p_loop3.bcs = vec![
        Bc { labels: washer_labs.clone(), kind: BcKind::LondonSheet { l_k_sq: 1.5e-12 } },
        Bc { labels: washer_labs, kind: BcKind::FluxLoop { name: "washer_loop".into(), flux_wb: 1.0 } },
    ];
    let r_loop3 = solve_magnetostatics(&g3.mesh, &p_loop3).unwrap();
    let ind_loop3 = InductanceMatrixData::of(&r_loop3).unwrap();
    assert_eq!(ind_loop3.terminals, vec!["washer_loop"]);
    assert!(ind_loop3.inductance_matrix_h[0][0] > 0.0 && ind_loop3.inductance_matrix_h[0][0].is_finite());
}
