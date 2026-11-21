// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Integration tests for `AF-1` (`femcore::adapt::electro`, interior-sheet `adapt3` / `Sig::Sheet`,
//! and `femcore::meshing::sheet3d_amr::OverlayOracle2`).

use femcore::adapt::electro::{
    combine_goal_oriented, electro_goal_indicators, electro_metrics_3d, electro_overlay_metrics_2d,
    electro_term_indicators, project_indicators_to_overlay_2d, ElectroEstimatorOptions, ElectroMetric2dOptions,
    ElectroMetric3dOptions, Permittivity, DEFAULT_EPS_FLOOR, EPS0,
};
use femcore::adapt::metric::sizes;
use femcore::fem::{
    assemble_boundary_mass, assemble_laplace, build_face_to_cell_3d, build_face_to_cells_3d,
    dirichlet_constraints, facet_measure, solve_spd, CurvedTet, TET_EDGES,
};
use femcore::props::LinearSolver;
use femcore::material::M3;
use femcore::mesh::SimplexMesh;
use femcore::meshing::sheet3d_amr::{
    adapt_overlay_2d, adapt_sheet_full_3d, build_overlay_mesh_2d, extract_overlay_sheet_features, extrude_overlay_2d,
    sheet_rim_vertices_3d,
};
use femcore::meshing::{mesh_model, Adaptive, MeshPlan};
use geom::expr::Expr;
use geom::model::Model;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

fn tri_area3(pts: &[f64], t: &[u32]) -> f64 {
    let p = |v: u32| [pts[3 * v as usize], pts[3 * v as usize + 1], pts[3 * v as usize + 2]];
    let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

/// Solves a multi-terminal P1/P2 Laplace capacitance problem directly on `mesh` using pure `femcore::fem`.
fn solve_terminals_capacitance(
    mesh: &SimplexMesh,
    order: usize,
    eps_r: &[f64],
    term_labels: &[Vec<i32>],
    ground_labels: &[i32],
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let coeff: Vec<f64> = eps_r.iter().map(|&e| EPS0 * e).collect();
    let k = assemble_laplace(mesh, order, &coeff);
    let n = k.nrows();
    let zero_rhs = vec![0.0; n];
    let n_term = term_labels.len();
    let mut sols = Vec::with_capacity(n_term);

    for ti in 0..n_term {
        let mut volt = BTreeMap::new();
        for &gl in ground_labels {
            volt.insert(gl, 0.0);
        }
        for (tj, labs) in term_labels.iter().enumerate() {
            let val = if tj == ti { 1.0 } else { 0.0 };
            for &l in labs {
                volt.entry(l).or_insert(val);
            }
        }
        let (fixed, vals) = dirichlet_constraints(mesh, &volt, order);
        let (u, _) = solve_spd(mesh, order, 1, LinearSolver::Auto, &k, &zero_rhs, &fixed, &vals).unwrap();
        sols.push(u);
    }

    let mut cmat = vec![vec![0.0; n_term]; n_term];
    for j in 0..n_term {
        let kvj = k.matvec(&sols[j]);
        for i in 0..n_term {
            cmat[i][j] = sols[i].iter().zip(&kvj).map(|(a, b)| a * b).sum();
        }
    }
    (sols, cmat)
}

#[test]
fn face_to_cell_3d_and_face_to_cells_3d_sorted_local_indices() {
    let m = femcore::meshing::geomdrive::block(1.0, 2.0, 3.0, 0.45, 3).unwrap().mesh;
    let one = build_face_to_cell_3d(&m);
    let all = build_face_to_cells_3d(&m);
    assert_eq!(one.len(), all.len());

    for (key, cs) in &all {
        assert!(key[0] < key[1] && key[1] < key[2], "key not sorted: {key:?}");
        assert!((1..=2).contains(&cs.len()));
        if cs.len() == 2 {
            assert!(cs[0].0 < cs[1].0, "cell indices must be ascending");
        }
        assert_eq!(one[key], cs[0]);
        for &(ci, lf) in cs {
            let v = m.cell(ci);
            assert_eq!([v[lf[0]], v[lf[1]], v[lf[2]]], *key, "local indices must map to sorted key");
        }
    }
}

#[test]
fn estimator_exact_linear_and_curved_isoparametric_zero_residual_and_tensor_identity() {
    let mut m = femcore::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.4, 3).unwrap().mesh;
    // Curve interior edges smoothly while keeping all elements valid (`CurvedTet`)
    let mut edge_nodes = Vec::new();
    for &(a, b) in m.edges() {
        let pa = m.point(a as usize);
        let pb = m.point(b as usize);
        let mid = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2) + (pb[2] - pa[2]).powi(2)).sqrt();
        let interior = mid.iter().all(|&c| c > 0.15 && c < 0.85);
        if interior {
            edge_nodes.push(([a, b], [mid[0] + 0.03 * len, mid[1] - 0.02 * len, mid[2] + 0.025 * len]));
        }
    }
    assert!(!edge_nodes.is_empty());
    m.set_edge_nodes(&edge_nodes);
    assert!(m.has_curved());
    for c in 0..m.n_cells() {
        if m.cell_is_curved(c) {
            let cv = m.cell(c);
            let v4: [[f64; 3]; 4] = std::array::from_fn(|i| {
                let p = m.point(cv[i] as usize);
                [p[0], p[1], p[2]]
            });
            let det0 = CurvedTet::straight(v4).jacobian(&[0.25, 0.25, 0.25, 0.25]).1;
            let ct = CurvedTet::new(m.cell_geometry_nodes(c));
            for &(i, j) in TET_EDGES.iter() {
                let mut lam = [0.0; 4];
                lam[i] = 0.5;
                lam[j] = 0.5;
                assert!(ct.jacobian(&lam).1 * det0 > 0.0);
            }
        }
    }

    // Exact physical linear potential V(x, y, z) = 1.7 x - 0.9 y + 2.3 z + 0.5 sampled at true P2 geometry nodes
    let nv = m.n_points();
    let eval_lin = |p: [f64; 3]| 1.7 * p[0] - 0.9 * p[1] + 2.3 * p[2] + 0.5;
    let mut u_p2 = Vec::with_capacity(m.n_dofs(2));
    for i in 0..nv {
        let p = m.point(i);
        u_p2.push(eval_lin([p[0], p[1], p[2]]));
    }
    for (ei, &(a, b)) in m.edges().iter().enumerate() {
        let pos = m.edge_node(ei).unwrap_or_else(|| {
            let (pa, pb) = (m.point(a as usize), m.point(b as usize));
            [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]
        });
        u_p2.push(eval_lin(pos));
    }

    let eps_scalar = vec![3.9; m.n_cells()];
    let all_dirichlet = vec![true; m.n_facets()];
    let opt_s = ElectroEstimatorOptions::new(2, Permittivity::Scalar(&eps_scalar), &all_dirichlet);
    let ind_s = electro_term_indicators(&m, std::slice::from_ref(&u_p2), &opt_s).unwrap();
    assert!(ind_s.vol_part[0] < 1e-34, "curved isoparametric linear field must have zero div residual, got {:e}", ind_s.vol_part[0]);
    assert!(ind_s.jump_part[0] < 1e-34, "curved isoparametric linear field must have zero jump residual, got {:e}", ind_s.jump_part[0]);
    assert_eq!(ind_s.robin_part[0], 0.0);

    // Scalar vs isotropic Tensor bit-identity on a non-linear field `V = x^2 + 0.5 y^2 + z^3`
    let u_nonlin: Vec<f64> = m
        .dof_coordinates(2)
        .chunks(3)
        .map(|p| p[0] * p[0] + 0.5 * p[1] * p[1] + p[2].powi(3))
        .collect();
    let no_dirichlet = vec![false; m.n_facets()];
    let robin_w = vec![2.5 * EPS0; m.n_facets()];
    let mut opt1 = ElectroEstimatorOptions::new(2, Permittivity::Scalar(&eps_scalar), &no_dirichlet);
    opt1.robin_alpha = Some(&robin_w);
    let res1 = electro_term_indicators(&m, std::slice::from_ref(&u_nonlin), &opt1).unwrap();

    let eps_tensor: Vec<M3> = vec![[[3.9, 0.0, 0.0], [0.0, 3.9, 0.0], [0.0, 0.0, 3.9]]; m.n_cells()];
    let mut opt2 = ElectroEstimatorOptions::new(2, Permittivity::Tensor(&eps_tensor), &no_dirichlet);
    opt2.robin_alpha = Some(&robin_w);
    let res2 = electro_term_indicators(&m, std::slice::from_ref(&u_nonlin), &opt2).unwrap();

    assert!(res1.vol_part[0] > 0.0 && res1.jump_part[0] > 0.0 && res1.robin_part[0] > 0.0);
    for c in 0..m.n_cells() {
        let a = res1.eta2_by_term[0][c];
        let b = res2.eta2_by_term[0][c];
        assert!(
            (a - b).abs() <= 1e-14 * a.abs().max(1e-300),
            "scalar vs isotropic tensor mismatch on cell {c}: {a:e} vs {b:e}"
        );
    }

    // Robin residual vanishes when `n . (eps0 eps_r grad V) + alpha V == 0` on a boundary face:
    // Take straight unit cube, `V(z) = 1.0 + z`, `eps_r = 2.0`. On top face `z = 1` (`n = +z`),
    // `n . (eps0 eps_r grad V) = +2 eps0`, while on bottom `z = 0` (`n = -z`), make `z = 0` and side walls Dirichlet,
    // and check Robin on `z = 1` (`V = 2.0`): `2 eps0 + alpha * 2 = 0` needs negative alpha, so test `V(z) = 2.0 - z`:
    // at `z = 1`, `V = 1.0`, `dV/dz = -1.0` => `n . D_grad = -2 eps0`, so `alpha = +2 eps0` cancels identically!
    let m_flat = femcore::meshing::geomdrive::block(1.0, 1.0, 1.0, 0.5, 3).unwrap().mesh;
    let u_rob: Vec<f64> = m_flat.dof_coordinates(1).chunks(3).map(|p| 2.0 - p[2]).collect();
    let mut dir_mask = vec![true; m_flat.n_facets()];
    let mut rob_alpha = vec![0.0; m_flat.n_facets()];
    let mut n_top = 0usize;
    for f in 0..m_flat.n_facets() {
        let fv = m_flat.facet(f);
        if fv.iter().all(|&vi| (m_flat.point(vi as usize)[2] - 1.0).abs() < 1e-12) {
            dir_mask[f] = false;
            rob_alpha[f] = 2.0 * EPS0;
            n_top += 1;
        }
    }
    assert!(n_top > 0);
    let eps2 = vec![2.0; m_flat.n_cells()];
    let mut opt_rob = ElectroEstimatorOptions::new(1, Permittivity::Scalar(&eps2), &dir_mask);
    opt_rob.robin_alpha = Some(&rob_alpha);
    let ind_rob = electro_term_indicators(&m_flat, std::slice::from_ref(&u_rob), &opt_rob).unwrap();
    assert!(ind_rob.robin_part[0] < 1e-36, "matched Robin BC must give zero residual, got {:e}", ind_rob.robin_part[0]);
    let _ = (assemble_boundary_mass, facet_measure);
}

#[test]
fn estimator_dielectric_interface_jump_and_dirichlet_sheet_exclusion() {
    // Two-layer box z in [-1, 1] with interface/sheet at z = 0
    let mut m = Model::new_3d();
    let bot = m.box_("bot", (0.0, 0.0, -1.0), (1.0, 1.0, 0.0));
    let top = m.box_("top", (0.0, 0.0, 0.0), (1.0, 1.0, 1.0));
    let pad_poly = m.polygons(
        "pad_poly",
        vec![vec![
            [Expr::Num(0.25), Expr::Num(0.25)],
            [Expr::Num(0.75), Expr::Num(0.25)],
            [Expr::Num(0.75), Expr::Num(0.75)],
            [Expr::Num(0.25), Expr::Num(0.75)],
        ]],
        vec![],
    );
    m.sheet("pad", pad_poly, 0.0);
    m.overlay(&[("sub", bot), ("vac", top)]);

    let g = mesh_model(&m.list, &MeshPlan::simple(0.35)).unwrap();
    let mesh = &g.mesh;
    let sub_dom = g.domain("sub").unwrap();
    let eps_r: Vec<f64> = mesh.cell_labels.iter().map(|&d| if d == sub_dom { 11.7 } else { 1.0 }).collect();

    // Potential `V(z) = |z|` has `dV/dz = +1` above z=0 (`eps_r = 1`) and `dV/dz = -1` below z=0 (`eps_r = 11.7`),
    // so normal flux jump across z=0 is `|1 * (+1) - 11.7 * (-1)| = 12.7` everywhere at z = 0.
    let u: Vec<f64> = mesh.dof_coordinates(1).chunks(3).map(|p| p[2].abs()).collect();

    let pad_labs: BTreeSet<i32> = g.boundary("pad").into_iter().collect();
    assert!(!pad_labs.is_empty());

    // Without marking `pad` as Dirichlet, jump runs over the full 1.0 m^2 z=0 plane;
    // with `pad` (area 0.25 m^2) marked Dirichlet, those interior sheet triangles are excluded!
    let dir_outer_only = mesh.exterior_facets();
    let mut dir_with_sheet = dir_outer_only.clone();
    for (d, &lab) in dir_with_sheet.iter_mut().zip(&mesh.facet_labels) {
        if pad_labs.contains(&lab) {
            *d = true;
        }
    }

    let opt_all = ElectroEstimatorOptions::new(1, Permittivity::Scalar(&eps_r), &dir_outer_only);
    let opt_excl = ElectroEstimatorOptions::new(1, Permittivity::Scalar(&eps_r), &dir_with_sheet);
    let ind_all = electro_term_indicators(mesh, std::slice::from_ref(&u), &opt_all).unwrap();
    let ind_excl = electro_term_indicators(mesh, std::slice::from_ref(&u), &opt_excl).unwrap();

    assert!(ind_all.jump_part[0] > 0.0);
    assert!(
        ind_excl.jump_part[0] < 0.85 * ind_all.jump_part[0] && ind_excl.jump_part[0] > 0.5 * ind_all.jump_part[0],
        "excluding 25% sheet area must reduce z=0 jump integral correspondingly: excl={:e} vs all={:e}",
        ind_excl.jump_part[0],
        ind_all.jump_part[0]
    );
}

#[test]
fn goal_oriented_cauchy_schwarz_weights_and_eps_floor() {
    let terms = femcore::adapt::electro::ElectroTermIndicators {
        eta2_by_term: vec![vec![4.0e-30, 0.0], vec![9.0e-30, 16.0e-30]],
        vol_part: vec![0.0, 0.0],
        jump_part: vec![4.0e-30, 25.0e-30],
        robin_part: vec![0.0, 0.0],
    };
    // C_00 = 100 fF, C_11 = 400 fF => sqrt(C_00 C_11) = 200 fF; C_01 = -1e-19 F < 1e-3 * 200 fF = 0.2 fF
    let cmat = vec![vec![100e-15, -1.0e-19], vec![-1.0e-19, 400e-15]];
    let goal = combine_goal_oriented(terms, &cmat, DEFAULT_EPS_FLOOR).unwrap();

    let w00 = 1.0 / 100e-15;
    let w11 = 1.0 / 400e-15;
    let w01 = 1.0 / (1e-3 * 200e-15);
    let want0 = w00 * 4.0e-30 + w01 * (2.0e-15 * 3.0e-15) + w11 * 9.0e-30;
    let want1 = w11 * 16.0e-30;
    assert!((goal.eta2[0] - want0).abs() < 1e-12 * want0);
    assert!((goal.eta2[1] - want1).abs() < 1e-12 * want1);
    assert!((goal.total - (want0 + want1).sqrt()).abs() < 1e-12 * goal.total);
}

#[test]
fn overlay_2d_metrics_corner_ball_prevents_anisotropic_collision_and_clamps_hmin() {
    let m2 = femcore::meshing::geomdrive::plate(1.0, 1.0, 0.05).unwrap().mesh;
    // L-shaped 90-degree corner at (0.5, 0.5) formed by two orthogonal rim segments
    let rim_segs = [([0.5, 0.5], [0.9, 0.5]), ([0.5, 0.5], [0.5, 0.9])];
    let corners = [[0.5, 0.5]];
    let eta2_2d: Vec<f64> = (0..m2.n_cells())
        .map(|t| {
            let v = m2.cell(t);
            let cx = (m2.point(v[0] as usize)[0] + m2.point(v[1] as usize)[0] + m2.point(v[2] as usize)[0]) / 3.0;
            let cy = (m2.point(v[0] as usize)[1] + m2.point(v[1] as usize)[1] + m2.point(v[2] as usize)[1]) / 3.0;
            let d = (cx - 0.5).hypot(cy - 0.5);
            (-12.0 * d).exp()
        })
        .collect();

    let mut opt = ElectroMetric2dOptions::new(2, 0.2, 0.05, 2.5 * m2.n_cells() as f64);
    opt.rim_anisotropy = 8.0;
    opt.max_aspect = 20.0;
    opt.rim_segments = &rim_segs;
    opt.corners = &corners;

    let met = electro_overlay_metrics_2d(&m2, &eta2_2d, &opt, None).unwrap();
    let h_min = 0.2 * 0.05;
    let mut checked_corner = false;
    let mut saw_aniso_away = false;
    for (vi, m_vi) in met.iter().enumerate().take(m2.n_points()) {
        let p = m2.point(vi);
        let (hs, _) = sizes(m_vi);
        assert!(hs[0] >= h_min * (1.0 - 1e-9) && hs[1] <= 0.2 * (1.0 + 1e-9), "clamps violated: {hs:?}");
        assert!(hs[1] / hs[0] <= 20.0 * (1.0 + 1e-9));
        let d_corner = (p[0] - 0.5).hypot(p[1] - 0.5);
        if d_corner < 0.03 {
            // Inside the inner corner ball, metric must be strictly isotropic!
            assert!((hs[1] / hs[0] - 1.0).abs() < 1e-6, "corner vertex {p:?} remained anisotropic: {hs:?}");
            checked_corner = true;
        }
        if (p[0] - 0.8).abs() < 0.05 && (p[1] - 0.5).abs() < 0.03 && hs[1] / hs[0] > 1.5 {
            saw_aniso_away = true;
        }
    }
    assert!(checked_corner, "expected vertices near (0.5, 0.5)");
    assert!(saw_aniso_away, "expected anisotropic stretching along straight rim away from corner");
}

#[test]
fn overlay_oracle_and_sheet_2d_z_multipass_amr_converges_and_preserves_circles() {
    // 3D stack with Si [-150, 0] um, Vac [0, 150] um, a circular sheet `disk` (radius 40 um at (-55, 0))
    // and a rectangular sheet `pad` ([15..95] x [-40..40] um) on z = 0, outer box grounded on top/bottom.
    let make_model = || -> Model {
        let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
        let mut m = Model::new_3d();
        let sub = m.box_("sub", (-180e-6, -150e-6, -150e-6), (180e-6, 150e-6, 0.0));
        let vac = m.box_("vac", (-180e-6, -150e-6, 0.0), (180e-6, 150e-6, 150e-6));
        let disk = m.circle("disk_prim", (-55e-6, 0.0), 40e-6);
        let pad = m.polygons(
            "pad_prim",
            vec![vec![pt(15e-6, -40e-6), pt(95e-6, -40e-6), pt(95e-6, 40e-6), pt(15e-6, 40e-6)]],
            vec![],
        );
        m.sheet("disk", disk, 0.0);
        m.sheet("pad", pad, 0.0);
        m.overlay(&[("substrate", sub), ("vacuum", vac)]);
        m
    };

    let model = make_model();
    let h_init = 45e-6;
    let plan = MeshPlan::simple(h_init);
    let ev0 = femcore::meshing::geomdrive::evaluate(&model.list).unwrap();
    let feats = extract_overlay_sheet_features(ev0.complex()).unwrap();
    assert!(!feats.rim_segments.is_empty());
    assert_eq!(feats.corners.len(), 4); // 4 sharp 90-degree rectangle corners (smooth circle has 0 corners)

    let mut m2 = build_overlay_mesh_2d(&ev0, &plan).unwrap();
    let mut totals = Vec::new();
    let mut c00_hist = Vec::new();
    let mut c01_hist = Vec::new();

    for pass in 0..4 {
        let ev_k = femcore::meshing::geomdrive::evaluate(&model.list).unwrap();
        let g = extrude_overlay_2d(ev_k, &m2, h_init, None).unwrap();
        let mesh = &g.mesh;

        // Check exact circular sheet perimeter projection on pass >= 1
        let disk_labs = g.boundary("disk");
        let rim3 = sheet_rim_vertices_3d(mesh, &disk_labs);
        let mut n_rim_disk = 0usize;
        for (vi, &on_r) in rim3.iter().enumerate() {
            if on_r {
                let p = mesh.point(vi);
                let r = (p[0] - (-55e-6)).hypot(p[1] - 0.0);
                assert!((r - 40e-6).abs() < 1e-10, "pass {pass}: disk rim vertex {p:?} off R=40um (r={r:e})");
                assert!(p[2].abs() < 1e-12);
                n_rim_disk += 1;
            }
        }
        assert!(n_rim_disk >= 6);

        let sub_dom = g.domain("substrate").unwrap();
        let eps_r: Vec<f64> = mesh.cell_labels.iter().map(|&d| if d == sub_dom { 11.7 } else { 1.0 }).collect();
        let pad_labs = g.boundary("pad");
        let mut gnd_labs = g.boundary("sub.bottom");
        gnd_labs.extend(g.boundary("vac.top"));
        let term_labs = vec![disk_labs.clone(), pad_labs.clone()];

        let (sols, cmat) = solve_terminals_capacitance(mesh, 2, &eps_r, &term_labs, &gnd_labs);
        let mut all_dir: BTreeSet<i32> = gnd_labs.iter().copied().collect();
        all_dir.extend(disk_labs.iter().copied());
        all_dir.extend(pad_labs.iter().copied());
        let dir_mask: Vec<bool> = mesh.facet_labels.iter().map(|l| all_dir.contains(l)).collect();

        let est_opt = ElectroEstimatorOptions::new(2, Permittivity::Scalar(&eps_r), &dir_mask);
        let goal = electro_goal_indicators(mesh, &sols, &cmat, DEFAULT_EPS_FLOOR, &est_opt).unwrap();
        totals.push(goal.total);
        c00_hist.push(cmat[0][0]);
        c01_hist.push(-cmat[0][1]);

        if pass < 3 {
            let eta2_2d = project_indicators_to_overlay_2d(mesh, &m2, &goal.eta2, feats.plane_axes).unwrap();
            let mut mopt = ElectroMetric2dOptions::new(2, h_init, 0.05, 1.5 * m2.n_cells() as f64);
            mopt.rim_segments = &feats.rim_segments;
            mopt.corners = &feats.corners;
            let met2 = electro_overlay_metrics_2d(&m2, &eta2_2d, &mopt, None).unwrap();
            let ev_ref = femcore::meshing::geomdrive::evaluate(&model.list).unwrap();
            let ov_ref = &ev_ref.complex().stack.as_ref().unwrap().overlay;
            let a = Adaptive { scale: 1.0, grade: 0.35, max_cells: Some(50_000) };
            let (m2_next, _) = adapt_overlay_2d(ov_ref, &m2, &met2, &a).unwrap();
            // Pinned overlay corners of `pad` must remain exact vertices of `m2_next`
            for &[cx, cy] in &[
                [15e-6, -40e-6],
                [95e-6, -40e-6],
                [95e-6, 40e-6],
                [15e-6, 40e-6],
            ] {
                let min_d = (0..m2_next.n_points())
                    .map(|i| (m2_next.point(i)[0] - cx).hypot(m2_next.point(i)[1] - cy))
                    .fold(f64::INFINITY, f64::min);
                assert!(min_d < 1e-12, "overlay corner ({cx:e}, {cy:e}) drifted by {min_d:e}");
            }
            m2 = m2_next;
        }
    }

    assert!(
        totals.windows(2).all(|w| w[1] < w[0]) && totals[3] < 0.72 * totals[0],
        "goal-oriented total indicator must decay strictly monotonically across sheet_2d_z AMR passes: {totals:?}"
    );
    assert!(c00_hist[3] > 0.0 && c01_hist[3] > 0.0);
    // Relative step on final pass is smaller than on first pass
    let d_first = (c00_hist[1] - c00_hist[0]).abs() / c00_hist[1];
    let d_last = (c00_hist[3] - c00_hist[2]).abs() / c00_hist[3];
    assert!(d_last < d_first, "capacitance changes should contract: d_first={d_first:e}, d_last={d_last:e}");
}

#[test]
fn full_3d_interior_sheet_adapt3_preserves_l0_eq_l1_and_l0_ne_l1_sheets_and_creases() {
    // 3D model with:
    // 1. `disk_iface`: circular sheet of radius 35 um at z = 0 on `substrate|vacuum` (`l0 != l1`),
    // 2. `pad_float`: rectangular sheet [10..70] x [-30..30] um at z = +50 um STRICTLY INSIDE `vacuum` (`l0 == l1`!).
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-120e-6, -100e-6, -80e-6), (120e-6, 100e-6, 0.0));
    let vac = m.box_("vac", (-120e-6, -100e-6, 0.0), (120e-6, 100e-6, 120e-6));
    let disk = m.circle("disk_prim", (-45e-6, 0.0), 35e-6);
    let pad = m.polygons(
        "pad_prim",
        vec![vec![pt(10e-6, -30e-6), pt(70e-6, -30e-6), pt(70e-6, 30e-6), pt(10e-6, 30e-6)]],
        vec![],
    );
    m.sheet("disk_iface", disk, 0.0);
    m.sheet("pad_float", pad, 50e-6);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let mut g = mesh_model(&m.list, &MeshPlan::simple(45e-6)).unwrap();
    let die_area = 240e-6 * 200e-6;
    let exact_pad_area = 60e-6 * 60e-6;
    let exact_disk_area = std::f64::consts::PI * (35e-6_f64).powi(2);

    let mut totals = Vec::new();
    for pass in 0..3 {
        let mesh = &g.mesh;
        let sub_dom = g.domain("substrate").unwrap();
        let vac_dom = g.domain("vacuum").unwrap();
        let disk_labs = g.boundary("disk_iface");
        let pad_labs = g.boundary("pad_float");
        let iface_labs = g.boundary("substrate|vacuum");

        // Verify single-facet emission, zero duplicates across ALL facets, exact domain sides, and exact areas!
        let mut tet_faces: HashMap<[u32; 3], Vec<i32>> = HashMap::new();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let dom = mesh.cell_labels[c];
            for skip in 0..4 {
                let mut f = [0u32; 3];
                let mut k = 0;
                for (i, &vi) in v.iter().enumerate() {
                    if i != skip {
                        f[k] = vi;
                        k += 1;
                    }
                }
                f.sort_unstable();
                tet_faces.entry(f).or_default().push(dom);
            }
        }
        let mut seen_facets: HashSet<[u32; 3]> = HashSet::new();
        for f in 0..mesh.n_facets() {
            let fv = mesh.facet(f);
            let mut key = [fv[0], fv[1], fv[2]];
            key.sort_unstable();
            assert!(seen_facets.insert(key), "pass {pass}: duplicate facet {key:?} in mesh.facets!");
        }

        // `pad_float` (`l0 == l1`): every facet shared by 2 `vac_dom` tets, exact area, exact z = 50 um
        let mut pad_area = 0.0;
        for &lab in &pad_labs {
            for tri in mesh.sheet_facets(lab) {
                pad_area += tri_area3(&mesh.points, tri);
                let mut key = [tri[0], tri[1], tri[2]];
                key.sort_unstable();
                assert_eq!(tet_faces[&key], vec![vac_dom, vac_dom], "pad_float must keep l0 == l1 == vac_dom");
                for &vi in tri {
                    assert!((mesh.point(vi as usize)[2] - 50e-6).abs() < 1e-12);
                }
            }
        }
        assert!(
            (pad_area - exact_pad_area).abs() < 1e-10 * exact_pad_area,
            "pass {pass}: pad_float area {pad_area:e} vs {exact_pad_area:e}"
        );

        // `disk_iface` (`l0 != l1`): shared by `sub_dom` and `vac_dom`, rim projected onto exact circle r = 35 um at z = 0
        let mut disk_area = 0.0;
        for &lab in &disk_labs {
            for tri in mesh.sheet_facets(lab) {
                disk_area += tri_area3(&mesh.points, tri);
                let mut key = [tri[0], tri[1], tri[2]];
                key.sort_unstable();
                let doms = &tet_faces[&key];
                assert_eq!(doms.len(), 2);
                assert!(doms.contains(&sub_dom) && doms.contains(&vac_dom));
            }
        }
        let mut iface_area = 0.0;
        for &lab in &iface_labs {
            for tri in mesh.sheet_facets(lab) {
                iface_area += tri_area3(&mesh.points, tri);
            }
        }
        assert!(
            (disk_area + iface_area - die_area).abs() < 1e-10 * die_area,
            "pass {pass}: z=0 complement sum {:e} vs die_area {die_area:e}",
            disk_area + iface_area
        );
        if pass == 2 {
            // Chord polygon area of disk_iface approaches exact pi * R^2 as rim creases split onto Complex::edges3
            assert!(
                (disk_area - exact_disk_area).abs() < 0.04 * exact_disk_area,
                "pass {pass}: disk_area={disk_area:e} vs exact={exact_disk_area:e}"
            );
        }
        let rim_disk = sheet_rim_vertices_3d(mesh, &disk_labs);
        for (vi, &is_r) in rim_disk.iter().enumerate() {
            if is_r {
                let p = mesh.point(vi);
                let r = (p[0] - (-45e-6)).hypot(p[1]);
                assert!((r - 35e-6).abs() < 1e-10, "pass {pass}: disk_iface rim vertex {p:?} off circle (r={r:e})");
                assert!(p[2].abs() < 1e-11, "pass {pass}: disk_iface rim z drifted: {:e}", p[2]);
            }
        }

        // Solve P1 electrostatics and adapt in full 3D
        let eps_r: Vec<f64> = mesh.cell_labels.iter().map(|&d| if d == sub_dom { 11.7 } else { 1.0 }).collect();
        let mut gnd_labs = g.boundary("sub.bottom");
        gnd_labs.extend(g.boundary("vac.top"));
        let term_labs = vec![disk_labs.clone(), pad_labs.clone()];
        let (sols, cmat) = solve_terminals_capacitance(mesh, 1, &eps_r, &term_labs, &gnd_labs);

        let mut all_dir: BTreeSet<i32> = gnd_labs.iter().copied().collect();
        all_dir.extend(disk_labs.iter().copied());
        all_dir.extend(pad_labs.iter().copied());
        let dir_mask: Vec<bool> = mesh.facet_labels.iter().map(|l| all_dir.contains(l)).collect();
        let est_opt = ElectroEstimatorOptions::new(1, Permittivity::Scalar(&eps_r), &dir_mask);
        let goal = electro_goal_indicators(mesh, &sols, &cmat, DEFAULT_EPS_FLOOR, &est_opt).unwrap();
        totals.push(goal.total);

        if pass < 2 {
            let mut sheet_all = disk_labs.clone();
            sheet_all.extend(pad_labs.iter().copied());
            let rim_all = sheet_rim_vertices_3d(mesh, &sheet_all);
            let mut opt3 = ElectroMetric3dOptions::new(1, 45e-6, 0.08, 1.8 * mesh.n_cells() as f64);
            opt3.rim_vertices = Some(&rim_all);
            let met3 = electro_metrics_3d(mesh, &goal.eta2, &opt3, None).unwrap();
            let a = Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(80_000) };
            let (g_next, _) = adapt_sheet_full_3d(g, &met3, &a).unwrap();
            g = g_next;
        }
    }

    assert!(totals[2] < totals[0], "full_3d AMR indicator must decrease across passes: {totals:?}");
}

fn compute_thread_determinism_digest() -> String {
    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-100e-6, -100e-6, -60e-6), (100e-6, 100e-6, 0.0));
    let vac = m.box_("vac", (-100e-6, -100e-6, 0.0), (100e-6, 100e-6, 60e-6));
    let p1 = m.polygons(
        "p1",
        vec![vec![pt(-60e-6, -30e-6), pt(-10e-6, -30e-6), pt(-10e-6, 30e-6), pt(-60e-6, 30e-6)]],
        vec![],
    );
    let p2 = m.polygons(
        "p2",
        vec![vec![pt(10e-6, -30e-6), pt(60e-6, -30e-6), pt(60e-6, 30e-6), pt(10e-6, 30e-6)]],
        vec![],
    );
    m.sheet("c0", p1, 0.0);
    m.sheet("c1", p2, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let plan = MeshPlan::simple(18e-6);
    let ev = femcore::meshing::geomdrive::evaluate(&m.list).unwrap();
    let feats = extract_overlay_sheet_features(ev.complex()).unwrap();
    let m2 = build_overlay_mesh_2d(&ev, &plan).unwrap();
    let g = extrude_overlay_2d(ev, &m2, 18e-6, None).unwrap();
    let mesh = &g.mesh;

    let sub_dom = g.domain("substrate").unwrap();
    let eps_r: Vec<f64> = mesh.cell_labels.iter().map(|&d| if d == sub_dom { 11.7 } else { 1.0 }).collect();
    let c0_labs = g.boundary("c0");
    let c1_labs = g.boundary("c1");
    let mut gnd_labs = g.boundary("sub.bottom");
    gnd_labs.extend(g.boundary("vac.top"));
    let (sols, cmat) = solve_terminals_capacitance(mesh, 2, &eps_r, &[c0_labs.clone(), c1_labs.clone()], &gnd_labs);

    let mut all_dir: BTreeSet<i32> = gnd_labs.into_iter().collect();
    all_dir.extend(c0_labs);
    all_dir.extend(c1_labs);
    let dir_mask: Vec<bool> = mesh.facet_labels.iter().map(|l| all_dir.contains(l)).collect();
    let est_opt = ElectroEstimatorOptions::new(2, Permittivity::Scalar(&eps_r), &dir_mask);
    let goal = electro_goal_indicators(mesh, &sols, &cmat, DEFAULT_EPS_FLOOR, &est_opt).unwrap();

    let eta2_2d = project_indicators_to_overlay_2d(mesh, &m2, &goal.eta2, feats.plane_axes).unwrap();
    let mut mopt = ElectroMetric2dOptions::new(2, 18e-6, 0.05, 1.5 * m2.n_cells() as f64);
    mopt.rim_segments = &feats.rim_segments;
    mopt.corners = &feats.corners;
    let met2 = electro_overlay_metrics_2d(&m2, &eta2_2d, &mopt, None).unwrap();
    let opt3 = ElectroMetric3dOptions::new(2, 18e-6, 0.05, 1.5 * mesh.n_cells() as f64);
    let met3 = electro_metrics_3d(mesh, &goal.eta2, &opt3, None).unwrap();

    let mut out = format!("total={:016x};", goal.total.to_bits());
    for &v in &goal.eta2 {
        out.push_str(&format!("{:016x},", v.to_bits()));
    }
    out.push(';');
    for m in &met2 {
        for r in m {
            for &x in r {
                out.push_str(&format!("{:016x},", x.to_bits()));
            }
        }
    }
    out.push(';');
    for m in &met3 {
        for r in m {
            for &x in r {
                out.push_str(&format!("{:016x},", x.to_bits()));
            }
        }
    }
    out
}

#[test]
fn thread_determinism_worker_entry() {
    if std::env::var_os("NYX_ELECTRO_AMR_WORKER").is_some() {
        println!("DIGEST:{}", compute_thread_determinism_digest());
    }
}

#[test]
fn thread_determinism_nyx_threads_1_vs_4_bit_identical() {
    let exe = std::env::current_exe().unwrap();
    let run_with_threads = |nt: usize| -> String {
        let out = std::process::Command::new(&exe)
            .env("NYX_THREADS", nt.to_string())
            .env("NYX_ELECTRO_AMR_WORKER", "1")
            .args(["--exact", "thread_determinism_worker_entry", "--nocapture"])
            .output()
            .expect("failed to spawn worker test process");
        assert!(out.status.success(), "worker failed: {}", String::from_utf8_lossy(&out.stderr));
        let stdout = String::from_utf8(out.stdout).unwrap();
        stdout
            .lines()
            .find_map(|l| l.strip_prefix("DIGEST:"))
            .expect("missing DIGEST line in worker output")
            .to_string()
    };

    let d1 = run_with_threads(1);
    let d4 = run_with_threads(4);
    assert_eq!(d1, d4, "NYX_THREADS=1 and NYX_THREADS=4 estimator/metric outputs must be bit-identical!");
}
