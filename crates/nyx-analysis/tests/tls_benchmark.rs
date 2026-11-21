// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Quantitative Acceptance Benchmark Test 2: Dielectric TLS Surface Participation (`TlsParticipationAnalysis`, Milestone 4.4).
//!
//! Benchmarks a superconducting coplanar waveguide (CPW) cross-section with center trace width `w = 10 um`,
//! gap `g = 6 um`, physical metal thickness `d = 100 nm`, native oxide thickness `t_ox = 2 nm` (`t_ma = t_ms = t_sa = 2 nm`,
//! `eps_ma = eps_sa = 10.0`, `eps_ms = 11.45`, `tan_delta = 1e-3`), on silicon (`eps_r = 11.45`):
//! 1. **2D finite-thickness (`d = 100 nm`) CPW cross-section**: compares the boundary perturbation formula
//!    (`TlsParticipationAnalysis`) against the direct volume integral over explicit `t_ox = 2 nm` oxide layers
//!    adjacent to the `MA`, `MS`, and `SA` interfaces, and verifies mesh-refinement convergence as `h -> 0`.
//! 2. **3D zero-thickness sheet CPW (`Op::Sheet`)**: runs `TlsParticipationAnalysis` with 3D matched asymptotic
//!    edge-tube regularization across a mesh-refinement sequence (`h -> 0`) and an `edge_cutoff` (`R_0`) stability
//!    sweep, asserting smooth convergence, `R_0`-invariance, and `< 2.0%` agreement with the 2D finite-thickness
//!    reference.

use femcore::fem::{facet_rule, geometry, quad_rule, shape};
use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_analysis::{tls_participation_with_params, Analysis, TlsParticipationAnalysis, TlsReport};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study, EPS0};
use nyx_physics::{Mesh, SolveResult};
use serde_json::json;

const W_TRACE: f64 = 10.0e-6;
const G_GAP: f64 = 6.0e-6;
const D_FILM: f64 = 100.0e-9;
const T_OX: f64 = 2.0e-9;
const EPS_SI: f64 = 11.45;
const EPS_OX: f64 = 10.0;
const HALF_BOX: f64 = 40.0e-6;
const H_SUB: f64 = 40.0e-6;
const H_AIR: f64 = 40.0e-6;

/// Builds and solves a 2D finite-thickness (`d = 100 nm`) CPW cross-section on `[-HALF_BOX, HALF_BOX] x [-H_SUB, H_AIR]`.
fn solve_cpw_2d_cross_section(h_global: f64, h_corner: f64) -> (MeshedGeom, Mesh, SolveResult) {
    let mut m = Model::new_2d();
    let sub = m.rect("substrate", (-HALF_BOX, -H_SUB), (HALF_BOX, 0.0));
    let vac = m.rect("vacuum", (-HALF_BOX, 0.0), (HALF_BOX, H_AIR));

    let a = W_TRACE / 2.0; // 5 um
    let b = a + G_GAP; // 11 um
    let trace_l = m.rect("trace_l", (-a, 0.0), (-a + 1.2e-6, D_FILM));
    let trace_m = m.rect("trace_m", (-a + 1.2e-6, 0.0), (a - 1.2e-6, D_FILM));
    let trace_r = m.rect("trace_r", (a - 1.2e-6, 0.0), (a, D_FILM));
    let gnd_l_far = m.rect("gnd_l_far", (-HALF_BOX, 0.0), (-b - 1.2e-6, D_FILM));
    let gnd_l_tip = m.rect("gnd_l_tip", (-b - 1.2e-6, 0.0), (-b, D_FILM));
    let gnd_r_tip = m.rect("gnd_r_tip", (b, 0.0), (b + 1.2e-6, D_FILM));
    let gnd_r_far = m.rect("gnd_r_far", (b + 1.2e-6, 0.0), (HALF_BOX, D_FILM));

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

    let g = mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(h_global)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

    // Identify facets on the center trace (`|x| <= a + tol, z in [0, d]`) and ground planes (`|x| >= b - tol, z in [0, d]`)
    let tol = 1e-8;
    let mut trace_labs = Vec::new();
    let mut gnd_labs = Vec::new();
    for f in 0..g.mesh.n_facets() {
        let fv = g.mesh.facet(f);
        let p0 = g.mesh.point(fv[0] as usize);
        let p1 = g.mesh.point(fv[1] as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let zm = 0.5 * (p0[1] + p1[1]);
        let lab = g.mesh.facet_labels[f];
        if zm >= -tol && zm <= D_FILM + tol {
            if xm.abs() <= a + tol {
                if !trace_labs.contains(&lab) {
                    trace_labs.push(lab);
                }
            } else if xm.abs() >= b - tol && !gnd_labs.contains(&lab) {
                gnd_labs.push(lab);
            }
        }
    }
    assert!(!trace_labs.is_empty() && !gnd_labs.is_empty());

    let mut prob = Problem::new(Study::Electrostatic, 2);
    prob.materials.set(sub_dom, "eps_r", EPS_SI);
    prob.materials.set(vac_dom, "eps_r", 1.0);
    prob.domain_names.insert(sub_dom, "substrate".into());
    prob.domain_names.insert(vac_dom, "vacuum".into());
    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];

    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let ex_mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    (g, ex_mesh, res)
}

/// Evaluates the direct 2D volume integral over explicit `t_ox = 2 nm` oxide layers (`Omega_MA`, `Omega_MS`, `Omega_SA`)
/// by 2D tensor-product Gauss quadrature across the 2 nm oxide thickness inside each adjacent P2 triangle.
fn explicit_2nm_oxide_volume_participations_2d(g: &MeshedGeom, res: &SolveResult) -> (f64, f64, f64) {
    let mesh = &g.mesh;
    let sub_dom = g.domain("substrate").unwrap();
    let sol = &res.solutions[0];
    let ex = sol.cell_fields.iter().find(|f| f.name == "Ex").unwrap();
    let ez = sol.cell_fields.iter().find(|f| f.name == "Ez").unwrap();
    let geo = geometry(mesh);
    let quad2 = quad_rule(2, 4);
    let nloc = ex.nloc;

    let eval_e = |ci: usize, lam: &[f64; 3]| -> [f64; 2] {
        let (sh, _) = shape(lam, 2);
        let base = ci * nloc;
        let (mut vx, mut vz) = (0.0, 0.0);
        for (a, &sh_a) in sh.iter().enumerate().take(nloc) {
            vx += sh_a * ex.re[base + a];
            vz += sh_a * ez.re[base + a];
        }
        [vx, vz]
    };

    // Total bulk energy U_elec
    let mut u_elec = 0.0;
    for ci in 0..mesh.n_cells() {
        let eps = if mesh.cell_labels[ci] == sub_dom { EPS_SI } else { 1.0 };
        let mut s = 0.0;
        for (q, l) in quad2.lam.iter().enumerate() {
            let e = eval_e(ci, &[l[0], l[1], l[2]]);
            s += quad2.w[q] * (e[0] * e[0] + e[1] * e[1]);
        }
        u_elec += 0.5 * EPS0 * eps * geo.vol[ci] * s;
    }

    // Build edge -> adjacent cells map
    type EdgeCellList = Vec<(usize, [usize; 2])>;
    let mut edge_cells: std::collections::HashMap<(u32, u32), EdgeCellList> = std::collections::HashMap::new();
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for [la, lb] in [[0, 1], [1, 2], [0, 2]] {
            let (va, vb) = (v[la], v[lb]);
            let (k, s) = if va < vb { ((va, vb), [la, lb]) } else { ((vb, va), [lb, la]) };
            edge_cells.entry(k).or_default().push((ci, s));
        }
    }

    let f_rule = facet_rule(1);
    // 2-point Gauss quadrature across the 2 nm oxide layer depth y_n in [0, T_OX]
    let g_depth = [0.5 * T_OX * (1.0 - 1.0 / 3.0_f64.sqrt()), 0.5 * T_OX * (1.0 + 1.0 / 3.0_f64.sqrt())];
    let w_depth = 0.5 * T_OX;

    // Sample inside the adjacent cell at normal offset y_n in [0, T_OX] by affine barycentric shift
    let eval_in_cell_at = |ci: usize, x_pt: [f64; 2]| -> [f64; 2] {
        let v = mesh.cell(ci);
        let p0 = mesh.point(v[0] as usize);
        let dx = [x_pt[0] - p0[0], x_pt[1] - p0[1]];
        let base = ci * 6;
        let mut lam = [0.0; 3];
        for (b, lam_b) in lam.iter_mut().enumerate().skip(1) {
            *lam_b = (geo.g[base + b * 2] * dx[0] + geo.g[base + b * 2 + 1] * dx[1]).clamp(0.0, 1.0);
        }
        lam[0] = (1.0 - lam[1] - lam[2]).clamp(0.0, 1.0);
        eval_e(ci, &lam)
    };

    let (mut u_ma, mut u_ms, mut u_sa) = (0.0, 0.0, 0.0);
    let a = W_TRACE / 2.0;
    let b = a + G_GAP;
    let tol = 1e-8;

    for (&(v0, v1), cells) in &edge_cells {
        let p0 = mesh.point(v0 as usize);
        let p1 = mesh.point(v1 as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let zm = 0.5 * (p0[1] + p1[1]);
        if (p0[0].abs() - HALF_BOX).abs() < tol && (p1[0].abs() - HALF_BOX).abs() < tol {
            continue;
        }
        let dx = p1[0] - p0[0];
        let dz = p1[1] - p0[1];
        let len = dx.hypot(dz);
        if len <= 0.0 {
            continue;
        }
        let tau = [dx / len, dz / len];

        let is_metal_edge = zm >= -tol && zm <= D_FILM + tol && (xm.abs() <= a + tol || xm.abs() >= b - tol);
        if is_metal_edge && cells.len() == 1 {
            let (ci, _) = cells[0];
            let vc = mesh.cell(ci);
            let cen = [
                (mesh.point(vc[0] as usize)[0] + mesh.point(vc[1] as usize)[0] + mesh.point(vc[2] as usize)[0]) / 3.0,
                (mesh.point(vc[0] as usize)[1] + mesh.point(vc[1] as usize)[1] + mesh.point(vc[2] as usize)[1]) / 3.0,
            ];
            let mut n_in = [-tau[1], tau[0]];
            if (cen[0] - xm) * n_in[0] + (cen[1] - zm) * n_in[1] < 0.0 {
                n_in = [-n_in[0], -n_in[1]];
            }
            let is_sub = mesh.cell_labels[ci] == sub_dom;
            let (eps_host, eps_ox) = if is_sub { (EPS_SI, EPS_SI) } else { (1.0, EPS_OX) };

            let mut strip_u = 0.0;
            for (q, l1) in f_rule.lam.iter().enumerate() {
                let xs = [l1[0] * p0[0] + l1[1] * p1[0], l1[0] * p0[1] + l1[1] * p1[1]];
                let ws = f_rule.w[q] * len;
                for &yn in &g_depth {
                    let x_pt = [xs[0] + yn * n_in[0], xs[1] + yn * n_in[1]];
                    let e_host = eval_in_cell_at(ci, x_pt);
                    let en_host = e_host[0] * n_in[0] + e_host[1] * n_in[1];
                    let et_host = e_host[0] * tau[0] + e_host[1] * tau[1];
                    // Inside the 2 nm oxide layer: D_n = eps_host E_{n,host} = eps_ox E_{n,ox}, E_{t,ox} = E_{t,host}
                    let en_ox = (eps_host / eps_ox) * en_host;
                    let e2_ox = en_ox * en_ox + et_host * et_host;
                    strip_u += 0.5 * EPS0 * eps_ox * ws * w_depth * e2_ox;
                }
            }
            if is_sub {
                u_ms += strip_u;
            } else {
                u_ma += strip_u;
            }
        } else if !is_metal_edge && cells.len() == 2 {
            let (c0, _) = cells[0];
            let (c1, _) = cells[1];
            let c_vac = if mesh.cell_labels[c0] != sub_dom && mesh.cell_labels[c1] == sub_dom {
                Some(c0)
            } else if mesh.cell_labels[c1] != sub_dom && mesh.cell_labels[c0] == sub_dom {
                Some(c1)
            } else {
                None
            };
            if let Some(ci) = c_vac {
                let n_in = [0.0, 1.0];
                let mut strip_u = 0.0;
                for (q, l1) in f_rule.lam.iter().enumerate() {
                    let xs = [l1[0] * p0[0] + l1[1] * p1[0], l1[0] * p0[1] + l1[1] * p1[1]];
                    let ws = f_rule.w[q] * len;
                    for &yn in &g_depth {
                        let x_pt = [xs[0] + yn * n_in[0], xs[1] + yn * n_in[1]];
                        let e_vac = eval_in_cell_at(ci, x_pt);
                        let en_vac = e_vac[0] * n_in[0] + e_vac[1] * n_in[1];
                        let et_vac = e_vac[0] * tau[0] + e_vac[1] * tau[1];
                        let en_ox = (1.0 / EPS_OX) * en_vac;
                        let e2_ox = en_ox * en_ox + et_vac * et_vac;
                        strip_u += 0.5 * EPS0 * EPS_OX * ws * w_depth * e2_ox;
                    }
                }
                u_sa += strip_u;
            }
        }
    }

    (u_ma / u_elec, u_ms / u_elec, u_sa / u_elec)
}

/// Builds and solves a 3D zero-thickness sheet CPW (`Op::Sheet`) with identical cross-section `(w = 10 um, g = 6 um)`.
fn solve_cpw_3d_sheet(h_global: f64, h_edge: f64) -> (MeshedGeom, Mesh, SolveResult) {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let l_y = 12.0e-6;
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

    let g = mesh_model(&m.list, &MeshPlan::simple(h_global)).unwrap();
    let sub_dom = g.domain("substrate").unwrap();
    let vac_dom = g.domain("vacuum").unwrap();

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
    (g, ex_mesh, res)
}

fn tls_params_with_cutoff(r0: Option<f64>) -> serde_json::Value {
    let mut v = json!({
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
        "d": D_FILM,
    });
    if let Some(r) = r0 {
        v["edge_cutoff"] = json!(r);
    }
    v
}

#[test]
fn test_2_cpw_2d_boundary_vs_explicit_oxide_and_3d_matched_asymptotic_convergence() {
    let params = tls_params_with_cutoff(None);

    // 1. 2D finite-thickness (d = 100 nm) CPW mesh convergence & explicit 2 nm oxide check
    let (_, mesh2_coarse, res2_coarse) = solve_cpw_2d_cross_section(4.0e-6, 45.0e-9);
    let (g2_fine, mesh2_fine, res2_fine) = solve_cpw_2d_cross_section(3.5e-6, 35.0e-9);
    assert!(TlsParticipationAnalysis.applies_to(&res2_fine));

    let rep2_coarse = tls_participation_with_params(&mesh2_coarse, &res2_coarse, &params).unwrap();
    let rep2_fine = tls_participation_with_params(&mesh2_fine, &res2_fine, &params).unwrap();
    let (p_ma_ox, p_ms_ox, p_sa_ox) = explicit_2nm_oxide_volume_participations_2d(&g2_fine, &res2_fine);
    let p_tot_ox = p_ma_ox + p_ms_ox + p_sa_ox;

    let diff_2d_h = (rep2_fine.p_total - rep2_coarse.p_total).abs() / rep2_fine.p_total;
    let diff_2d_ox = (rep2_fine.p_total - p_tot_ox).abs() / p_tot_ox;
    eprintln!(
        "Test 2 [2D finite-thickness d=100nm]:\n  2D boundary (coarse): p_ma={:.5e}, p_ms={:.5e}, p_sa={:.5e}, p_tot={:.5e}\n  2D boundary (fine):   p_ma={:.5e}, p_ms={:.5e}, p_sa={:.5e}, p_tot={:.5e}, Q_tls={:.1}\n  2D explicit 2nm ox:   p_ma={:.5e}, p_ms={:.5e}, p_sa={:.5e}, p_tot={:.5e}\n  2D h-refinement diff = {:.3}%, 2D boundary vs explicit 2nm oxide diff = {:.3}%",
        rep2_coarse.p_ma, rep2_coarse.p_ms, rep2_coarse.p_sa, rep2_coarse.p_total,
        rep2_fine.p_ma, rep2_fine.p_ms, rep2_fine.p_sa, rep2_fine.p_total, rep2_fine.q_tls,
        p_ma_ox, p_ms_ox, p_sa_ox, p_tot_ox,
        diff_2d_h * 100.0,
        diff_2d_ox * 100.0
    );
    assert!(diff_2d_h < 0.02, "2D mesh refinement shift {:.3}% >= 2.0%", diff_2d_h * 100.0);
    assert!(
        diff_2d_ox < 0.02,
        "2D boundary formula vs explicit 2nm oxide integral diff {:.3}% >= 2.0%",
        diff_2d_ox * 100.0
    );

    // 2. 3D zero-thickness sheet CPW with matched asymptotic edge-tube regularization:
    //    check convergence across 3 mesh refinement levels (h_edge = 1.6 um -> 1.1 um -> 0.8 um)
    let r0_fixed = 1.5e-6;
    let params_r0 = tls_params_with_cutoff(Some(r0_fixed));
    let mut reps_3d: Vec<(usize, TlsReport)> = Vec::new();
    for &(h_glob, h_ed) in &[(12.0e-6, 1.6e-6), (10.0e-6, 1.1e-6), (8.0e-6, 0.8e-6)] {
        let (g3, m3, r3) = solve_cpw_3d_sheet(h_glob, h_ed);
        let rep = tls_participation_with_params(&m3, &r3, &params_r0).unwrap();
        eprintln!(
            "Test 2 [3D sheet h_edge={:.2}um, {} cells, R0={:.2}um]: p_ma={:.5e}, p_ms={:.5e}, p_sa={:.5e}, p_tot={:.5e}",
            h_ed * 1e6,
            g3.mesh.n_cells(),
            r0_fixed * 1e6,
            rep.p_ma,
            rep.p_ms,
            rep.p_sa,
            rep.p_total
        );
        reps_3d.push((g3.mesh.n_cells(), rep));
    }

    let rep3_med = &reps_3d[1].1;
    let rep3_fine = &reps_3d[2].1;
    let diff_3d_h = (rep3_fine.p_total - rep3_med.p_total).abs() / rep3_fine.p_total;
    eprintln!("Test 2 [3D sheet h-refinement convergence (medium -> fine)]: {:.3}%", diff_3d_h * 100.0);
    assert!(diff_3d_h < 0.02, "3D matched-asymptotic h-refinement shift {:.3}% >= 2.0%", diff_3d_h * 100.0);

    // 3. Verify R_0 (`edge_cutoff`) invariance on the fine 3D mesh across R_0 in [1.1 um, 1.5 um, 1.9 um]
    let (_, m3_fine, r3_fine) = solve_cpw_3d_sheet(8.0e-6, 0.8e-6);
    for &r0_test in &[1.1e-6, 1.5e-6, 1.9e-6] {
        let rep_r = tls_participation_with_params(&m3_fine, &r3_fine, &tls_params_with_cutoff(Some(r0_test))).unwrap();
        let rel_r0 = (rep_r.p_total - rep3_fine.p_total).abs() / rep3_fine.p_total;
        eprintln!(
            "  R0 stability check (R0 = {:.2} um): p_tot = {:.5e} (rel diff vs R0=1.5um: {:.3}%)",
            r0_test * 1e6,
            rep_r.p_total,
            rel_r0 * 100.0
        );
        assert!(rel_r0 < 0.02, "R_0 sensitivity {:.3}% exceeds 2.0% at R0={r0_test:e}", rel_r0 * 100.0);
    }

    // 4. Verify < 2.0% agreement between 3D matched-asymptotic and 2D finite-thickness (d = 100 nm) reference
    let err_p_tot = (rep3_fine.p_total - rep2_fine.p_total).abs() / rep2_fine.p_total;
    let err_p_ms = (rep3_fine.p_ms - rep2_fine.p_ms).abs() / rep2_fine.p_ms;
    let err_p_sa = (rep3_fine.p_sa - rep2_fine.p_sa).abs() / rep2_fine.p_sa;
    let err_p_ma = (rep3_fine.p_ma - rep2_fine.p_ma).abs() / rep2_fine.p_ma;
    eprintln!(
        "Test 2 [3D matched-asymptotic vs 2D reference]:\n  p_total err = {:.3}%\n  p_ms err    = {:.3}%\n  p_sa err    = {:.3}%\n  p_ma err    = {:.3}%",
        err_p_tot * 100.0,
        err_p_ms * 100.0,
        err_p_sa * 100.0,
        err_p_ma * 100.0
    );
    assert!(
        err_p_tot < 0.02 && err_p_ms < 0.02 && err_p_sa < 0.02 && err_p_ma < 0.02,
        "3D matched-asymptotic vs 2D finite-thickness error exceeds 2.0%: p_tot={:.3}%, p_ms={:.3}%, p_sa={:.3}%, p_ma={:.3}%",
        err_p_tot * 100.0,
        err_p_ms * 100.0,
        err_p_sa * 100.0,
        err_p_ma * 100.0
    );
}
