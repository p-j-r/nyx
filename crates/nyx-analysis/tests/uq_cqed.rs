// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 4 VVUQ verification tests for zero-solve cQED uncertainty quantification (`nyx-analysis::uq`):
//! - Zero-solve analytical/FD sensitivities + deterministic `SplitMix64` Monte Carlo over `CircuitRom` (`epr` +
//!   `hamiltonian` in `normal_fock` and `craig_bampton`), `CapacitanceMatrixData` (+ 3D thickness Gram cache),
//!   and 2D/3D dielectric TLS (`p_ma`, `p_ms`, `p_sa`, `loss_tan_total`, `Q_TLS`);
//! - Flux sweet-spot behavior (`linear_sigma == 0` while `mc_sigma > 0` captures quadratic dephasing cleanly);
//! - Strict ASME V&V 20 epistemic vs aleatoric separation: Prager–Synge `[low, up]` and `error_budget`
//!   discretization brackets NEVER collapse into Gaussian sigmas, remain invariant when fabrication sigmas change,
//!   do not perturb aleatoric sigmas when widened 10x, and compose into explicit `outer_envelope` p-box brackets.
#![allow(clippy::needless_range_loop)]

use femcore::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts};
use femcore::meshing::{mesh_model, MeshPlan};
use geom::expr::Expr;
use geom::model::{Model, RegionKind};
use geom::select::Selector;
use nyx_analysis::epr::{jj_flux_scale, project_rom_modes, C0, EPS0, MU0};
use nyx_analysis::{
    hamiltonian_with_params, tls_participation_with_params, uq_cqed_with_companion_and_params,
    uq_cqed_with_mesh_and_params, uq_cqed_with_params, Analysis, UqCqedAnalysis,
};
use nyx_dfs::{solve_electrostatics, Bc, BcKind, Problem, Study};
use nyx_physics::{
    CapacitanceMatrixData, CircuitRom, Mesh, NumericalErrorBudget, RomPort, SolveResult, SpatialErrorEstimate,
};
use serde_json::{json, Value};
use std::f64::consts::PI;

fn mutual_from_maxwell(c: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = c.len();
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        let mut s = 0.0;
        for (k, &v) in c[i].iter().enumerate() {
            s += v;
            if i != k {
                m[i][k] = -v;
            }
        }
        m[i][i] = s;
    }
    m
}

fn mat_add_rank1(base: &[Vec<f64>], v: &[f64], scale: f64) -> Vec<Vec<f64>> {
    let n = base.len();
    let mut out = base.to_vec();
    for i in 0..n {
        for j in 0..n {
            out[i][j] += scale * v[i] * v[j];
        }
    }
    out
}

fn make_transmon_resonator_rom(ej_hz: f64, squid_asymmetry: Option<f64>, c_q_f: f64) -> SolveResult {
    let scale = jj_flux_scale();
    let c_r = 50.0e-15_f64;
    let c_g = 1.8e-15_f64;
    let omega_r = 2.0 * PI * 6.8e9;
    let l_r_inv = omega_r * omega_r * (c_r + c_g);
    let b_feed = 0.025_f64;
    let d_port = (MU0 * C0 / 50.0) * b_feed * b_feed;

    let mut rom = CircuitRom {
        ports: vec![
            RomPort {
                name: "Q1_JJ".into(),
                kind: "lumped_jj".into(),
                l_j_h: Some(scale / ej_hz),
                e_j_hz: Some(ej_hz),
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry,
                z0_ohm: None,
                direction: [0.0, 1.0, 0.0],
                formulation: "rank1".into(),
            },
            RomPort {
                name: "ro_port".into(),
                kind: "lumped_port".into(),
                l_j_h: None,
                e_j_hz: None,
                c_j_f: 0.0,
                e_l_hz: None,
                phi_ext_rad: 0.0,
                n_g: 0.0,
                squid_asymmetry: None,
                z0_ohm: Some(50.0),
                direction: [1.0, 0.0, 0.0],
                formulation: "rank1".into(),
            },
        ],
        r_static: 0,
        r_dyn: 2,
        k_r: vec![vec![0.0, 0.0], vec![0.0, MU0 * l_r_inv]],
        c_r: vec![vec![0.0, 0.0], vec![0.0, d_port]],
        m_r: vec![vec![(c_q_f + c_g) / EPS0, -c_g / EPS0], vec![-c_g / EPS0, (c_r + c_g) / EPS0]],
        b_r: vec![vec![1.0, 0.0], vec![0.0, b_feed]],
        c_port_f: vec![],
        modes: vec![],
    };
    rom.modes = project_rom_modes(&rom, &rom.ports).unwrap();
    SolveResult {
        physics: "ems".into(),
        study: "cavity".into(),
        dim: 3,
        scalars: vec![rom.scalar()],
        ..Default::default()
    }
}

fn solve_fast_2d_cpw_tls() -> (Mesh, SolveResult) {
    let mut m = Model::new_2d();
    let sub = m.rect("substrate", (-30.0e-6, -25.0e-6), (30.0e-6, 0.0));
    let vac = m.rect("vacuum", (-30.0e-6, 0.0), (30.0e-6, 25.0e-6));
    let trace = m.rect("trace", (-4.0e-6, 0.0), (4.0e-6, 0.25e-6));
    let gnd_l = m.rect("gnd_l", (-30.0e-6, 0.0), (-9.0e-6, 0.25e-6));
    let gnd_r = m.rect("gnd_r", (9.0e-6, 0.0), (30.0e-6, 0.25e-6));

    m.overlay_full(
        &[
            ("substrate", sub, RegionKind::Material, None),
            ("vacuum", vac, RegionKind::Material, None),
            ("trace", trace, RegionKind::Hole, None),
            ("gnd_l", gnd_l, RegionKind::Hole, None),
            ("gnd_r", gnd_r, RegionKind::Hole, None),
        ],
        &[],
    );
    let eg = evaluate(&m.list).unwrap();
    let mut gm = mesh_geom_2d(eg, &GeomMeshOpts::new(2.5e-6)).unwrap();

    let sub_d = gm.domain("substrate").unwrap();
    let vac_d = gm.domain("vacuum").unwrap();
    let ext = gm.mesh.exterior_facets();
    for (f, &is_ext) in ext.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let fv = gm.mesh.facet(f);
        let p0 = gm.mesh.point(fv[0] as usize);
        let p1 = gm.mesh.point(fv[1] as usize);
        let xm = 0.5 * (p0[0] + p1[0]);
        let zm = 0.5 * (p0[1] + p1[1]);
        let is_outer = ((p0[0].abs() - 30.0e-6).abs() < 1e-11 && (p1[0].abs() - 30.0e-6).abs() < 1e-11)
            || ((p0[1].abs() - 25.0e-6).abs() < 1e-11 && (p1[1].abs() - 25.0e-6).abs() < 1e-11);
        if is_outer {
            gm.mesh.facet_labels[f] = 300;
        } else if (-1e-11..=0.25e-6 + 1e-11).contains(&zm) && xm.abs() <= 4.0e-6 + 1e-11 {
            gm.mesh.facet_labels[f] = 100;
        } else if (-1e-11..=0.25e-6 + 1e-11).contains(&zm) && xm.abs() >= 9.0e-6 - 1e-11 {
            gm.mesh.facet_labels[f] = 200;
        } else {
            gm.mesh.facet_labels[f] = 300;
        }
    }
    let trace_labs = vec![100];
    let gnd_labs = vec![200];

    let mut prob = Problem::new(Study::Electrostatic, 1);
    prob.materials.set(sub_d, "eps_r", 11.45);
    prob.domain_names.insert(sub_d, "substrate".into());
    prob.materials.set(vac_d, "eps_r", 1.0);
    prob.domain_names.insert(vac_d, "vacuum".into());
    prob.bcs = vec![
        Bc { labels: trace_labs, kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd_labs, kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&gm.mesh, &prob).unwrap();
    let mesh = Mesh::from_solver_mesh(&gm.mesh).unwrap();
    (mesh, res)
}

fn solve_fast_3d_cpw_sheet() -> (Mesh, SolveResult) {
    let mut m = Model::new_3d();
    let hb = 30.0e-6;
    let ly = 2.5e-6;
    let a = 4.0e-6;
    let b = 9.0e-6;
    let sub = m.box_("substrate", (-hb, 0.0, -30.0e-6), (hb, ly, 0.0));
    let vac = m.box_("vacuum", (-hb, 0.0, 0.0), (hb, ly, 30.0e-6));

    let pt = |x: f64, y: f64| [Expr::Num(x), Expr::Num(y)];
    let trace_p = m.polygons("trace_p", vec![vec![pt(-a, 0.0), pt(a, 0.0), pt(a, ly), pt(-a, ly)]], vec![]);
    let gnd_l_p = m.polygons("gnd_l_p", vec![vec![pt(-hb, 0.0), pt(-b, 0.0), pt(-b, ly), pt(-hb, ly)]], vec![]);
    let gnd_r_p = m.polygons("gnd_r_p", vec![vec![pt(b, 0.0), pt(hb, 0.0), pt(hb, ly), pt(b, ly)]], vec![]);
    m.sheet("trace", trace_p, 0.0);
    m.sheet("gnd_l", gnd_l_p, 0.0);
    m.sheet("gnd_r", gnd_r_p, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);
    for sel in ["trace:edge", "gnd_l:edge", "gnd_r:edge"] {
        m.size(
            &format!("sz_{}", sel.replace(':', "_")),
            Selector::parse(sel).unwrap(),
            Some(Expr::Num(1.2e-6)),
            None,
            None,
        );
    }
    let g = mesh_model(&m.list, &MeshPlan::simple(7.0e-6)).unwrap();

    let mut prob = Problem::new(Study::Electrostatic, 2);
    for d in g.domains("substrate") {
        prob.materials.set(d, "eps_r", 11.45);
        prob.domain_names.insert(d, "substrate".into());
    }
    for d in g.domains("vacuum") {
        prob.materials.set(d, "eps_r", 1.0);
        prob.domain_names.insert(d, "vacuum".into());
    }
    let mut gnd = g.boundary("gnd_l");
    gnd.extend(g.boundary("gnd_r"));
    prob.bcs = vec![
        Bc { labels: g.boundary("trace"), kind: BcKind::Terminal("trace".into()) },
        Bc { labels: gnd, kind: BcKind::Ground },
    ];
    let res = solve_electrostatics(&g.mesh, &prob).unwrap();
    let mesh = Mesh::from_solver_mesh(&g.mesh).unwrap();
    (mesh, res)
}

#[test]
fn uq_cqed_rom_linear_vs_mc_determinism_and_flux_sweet_spot() {
    let res = make_transmon_resonator_rom(16.0e9, Some(0.25), 75.0e-15);
    assert!(UqCqedAnalysis.applies_to(&res));

    let params = json!({
        "basis": "normal_fock",
        "fock_trunc": [6, 5],
        "seed": 20260401u64,
        "n_samples": 600,
        "uncertainties": {
            "E_J_Hz": { "rel_tol": 0.02 },
            "Z0_ohm:ro_port": { "rel_tol": 0.10 }
        }
    });
    let rep_a = uq_cqed_with_params(&res, &params).unwrap();
    let rep_b = uq_cqed_with_params(&res, &params).unwrap();
    assert!(rep_a.zero_fem_solves);
    assert!(rep_a.epistemic_aleatoric_separated);

    let f01 = rep_a.quantity("f01_Hz").unwrap();
    let f01_b = rep_b.quantity("f01_Hz").unwrap();
    let ham_nom = hamiltonian_with_params(&res, &json!({ "basis": "normal_fock", "fock_trunc": [6, 5] })).unwrap();
    assert!((f01.nominal - ham_nom.f01_hz[0]).abs() < 1.0);

    let a_a = f01.aleatoric.as_ref().unwrap();
    let a_b = f01_b.aleatoric.as_ref().unwrap();
    assert_eq!(a_a.linear_sigma.to_bits(), a_b.linear_sigma.to_bits());
    assert_eq!(a_a.mc_mean.to_bits(), a_b.mc_mean.to_bits());
    assert_eq!(a_a.mc_sigma.to_bits(), a_b.mc_sigma.to_bits());
    assert_eq!(a_a.ci_95[0].to_bits(), a_b.ci_95[0].to_bits());
    assert_eq!(a_a.ci_95[1].to_bits(), a_b.ci_95[1].to_bits());

    // Weak 2% E_J uncertainty -> linear_sigma and mc_sigma agree closely (~1% shift on f01 ~ sqrt(8 E_J E_C))
    let rel_diff = (a_a.linear_sigma - a_a.mc_sigma).abs() / a_a.linear_sigma;
    assert!(rel_diff < 0.12, "linear {} vs mc {}", a_a.linear_sigma, a_a.mc_sigma);
    assert!(a_a.ci_95[0] < f01.nominal && f01.nominal < a_a.ci_95[1]);
    let sum_vc: f64 = a_a.variance_contributions.values().sum();
    assert!((sum_vc - 1.0).abs() < 1e-12);
    assert!(a_a.variance_contributions["E_J_Hz"] > 0.99);

    // Cavity kappa is dominated by Z0_ohm:ro_port (~100% of variance); Purcell T1 reflects both Z0_ohm and detuning via E_J_Hz
    let kappa_a = rep_a.quantity("kappa_cav_Hz").unwrap().aleatoric.as_ref().unwrap();
    assert!(kappa_a.variance_contributions["Z0_ohm:ro_port"] > 0.98);
    let t1_p = rep_a.quantity("t1_purcell_s").unwrap().aleatoric.as_ref().unwrap();
    assert!(t1_p.linear_sigma > 0.0);
    assert!(t1_p.variance_contributions["Z0_ohm:ro_port"] > 0.20);
    assert!(t1_p.variance_contributions["E_J_Hz"] > 0.20);
    assert!((t1_p.variance_contributions.values().sum::<f64>() - 1.0).abs() < 1e-12);

    // Also test craig_bampton basis compatibility
    let cb_rep = uq_cqed_with_params(
        &res,
        &json!({
            "basis": "craig_bampton",
            "port_trunc": 7,
            "n_samples": 64,
            "seed": 7,
            "uncertainties": { "L_J_H": { "rel_tol": 0.02 } }
        }),
    )
    .unwrap();
    let f01_cb = cb_rep.quantity("f01_Hz").unwrap();
    assert!((f01_cb.nominal - f01.nominal).abs() / f01.nominal < 1e-3);
    assert!(f01_cb.aleatoric.as_ref().unwrap().linear_sigma > 0.0);

    // Flux sweet spot at Phi_ext = 0: first derivative is 0, quadratic curvature drives mc_sigma > 0 with no NaNs
    let sweet_rep = uq_cqed_with_params(
        &res,
        &json!({
            "basis": "normal_fock",
            "fock_trunc": [6, 4],
            "seed": 99,
            "n_samples": 400,
            "squid_flux_phi0": 0.0,
            "uncertainties": { "squid_flux_phi0": 0.03 }
        }),
    )
    .unwrap();
    let f_sweet = sweet_rep.quantity("f01_Hz").unwrap();
    let a_sweet = f_sweet.aleatoric.as_ref().unwrap();
    assert!(a_sweet.linear_sigma.abs() < 1.0e2, "linear_sigma at sweet spot = {}", a_sweet.linear_sigma);
    assert!(a_sweet.mc_sigma > 1.0e6, "quadratic flux dephasing mc_sigma = {}", a_sweet.mc_sigma);
    assert!(a_sweet.mc_mean < f_sweet.nominal);
    assert!(a_sweet.variance_contributions.is_empty());
}

#[test]
fn uq_cqed_asme_vv20_strict_separation_between_prager_synge_and_aleatoric() {
    let f = 1.0e-15;
    let c_exact = vec![
        vec![190.0 * f, -35.0 * f, -55.0 * f, -60.0 * f],
        vec![-35.0 * f, 48.0 * f, -8.5 * f, -0.8 * f],
        vec![-55.0 * f, -8.5 * f, 85.0 * f, -20.0 * f],
        vec![-60.0 * f, -0.8 * f, -20.0 * f, 83.0 * f],
    ];
    let make_cap_res = |scale: f64| -> SolveResult {
        // Construct exact symmetric PSD brackets C_mid +/- scale * Delta around fixed nominal c_exact
        let v1 = [1.0, -0.4, 0.8, -0.7];
        let v2 = [0.3, 0.5, 0.6, 0.5];
        let mut c_low = mat_add_rank1(&c_exact, &v1, -3.0 * scale * f);
        c_low = mat_add_rank1(&c_low, &v2, -2.0 * scale * f);
        let mut c_up = mat_add_rank1(&c_exact, &v1, 3.0 * scale * f);
        c_up = mat_add_rank1(&c_up, &v2, 2.0 * scale * f);
        let cap = CapacitanceMatrixData {
            terminals: vec!["c0".into(), "c1".into(), "c2".into(), "c3".into()],
            maxwell_matrix_f: c_exact.clone(),
            mutual_matrix_f: mutual_from_maxwell(&c_exact),
            ec_matrix_hz: vec![vec![0.0; 4]; 4],
            maxwell_lower_matrix_f: Some(c_low),
            maxwell_upper_matrix_f: Some(c_up),
            ..Default::default()
        };
        SolveResult {
            physics: "dfs".into(),
            study: "electrostatic".into(),
            dim: 3,
            scalars: vec![cap.scalar()],
            ..Default::default()
        }
    };

    let cap_narrow = make_cap_res(0.5);
    let cap_wide = make_cap_res(4.0); // 8x wider Prager-Synge epistemic interval around the exact same midpoint!

    let cap_cfg = json!({
        "ground": "c0",
        "eliminate": "drop",
        "qubits": [{ "name": "Q1", "pads": ["c2", "c3"], "E_J_Hz": 16.0e9 }],
        "resonators": [{ "name": "R1", "terminal": "c1", "f_r_Hz": 8.2e9, "Z0_ohm": 50.0 }],
        "n_samples": 400,
        "seed": 11
    });

    // 1. Epistemic-only run (`uncertainties` omitted): NO fake Gaussian sigma should ever be synthesized
    let ep_only = uq_cqed_with_params(&cap_narrow, &cap_cfg).unwrap();
    let cs_only = ep_only.quantity("C_Sigma_F").unwrap();
    assert!(cs_only.aleatoric.is_none());
    let ep_cs = cs_only.epistemic.as_ref().unwrap();
    assert!(ep_cs.low < cs_only.nominal && cs_only.nominal < ep_cs.up);
    assert!(ep_only.quantity("f01_Hz").unwrap().epistemic.as_ref().unwrap().half_width > 0.0);
    assert!(ep_only.quantity("chi_Hz").unwrap().epistemic.as_ref().unwrap().half_width > 0.0);

    // 2. Widening Prager-Synge 8x MUST leave aleatoric sigmas 0-ULP invariant while growing epistemic + outer p-box
    let mut joint_narrow = cap_cfg.clone();
    joint_narrow["uncertainties"] = json!({ "E_J_Hz": { "rel_tol": 0.025 }, "f_r_Hz": { "rel_tol": 0.005 } });
    let r_narrow = uq_cqed_with_params(&cap_narrow, &joint_narrow).unwrap();
    let r_wide = uq_cqed_with_params(&cap_wide, &joint_narrow).unwrap();

    for qname in ["f01_Hz", "chi_Hz", "g_Hz"] {
        let qn = r_narrow.quantity(qname).unwrap();
        let qw = r_wide.quantity(qname).unwrap();
        let an = qn.aleatoric.as_ref().unwrap();
        let aw = qw.aleatoric.as_ref().unwrap();
        assert_eq!(an.linear_sigma.to_bits(), aw.linear_sigma.to_bits(), "{qname} linear_sigma polluted by epistemic");
        assert_eq!(an.mc_sigma.to_bits(), aw.mc_sigma.to_bits(), "{qname} mc_sigma polluted by epistemic");

        let en = qn.epistemic.as_ref().unwrap();
        let ew = qw.epistemic.as_ref().unwrap();
        assert!(ew.half_width > 4.0 * en.half_width);

        let pn = qn.outer_envelope.as_ref().unwrap();
        let pw = qw.outer_envelope.as_ref().unwrap();
        assert!(pw.bracket_95[0] < pn.bracket_95[0] && pw.bracket_95[1] > pn.bracket_95[1]);
        assert!(pw.bracket_95[0] <= aw.ci_95[0] && pw.bracket_95[1] >= aw.ci_95[1]);
    }

    // 3. Changing aleatoric E_J_Hz sigma 4x MUST leave epistemic Prager-Synge brackets unchanged
    let mut larger_aleatoric = joint_narrow.clone();
    larger_aleatoric["uncertainties"] = json!({ "E_J_Hz": { "rel_tol": 0.10 }, "f_r_Hz": { "rel_tol": 0.02 } });
    let r_big_al = uq_cqed_with_params(&cap_narrow, &larger_aleatoric).unwrap();
    assert_eq!(
        r_narrow.quantity("C_Sigma_F").unwrap().epistemic.as_ref().unwrap().low.to_bits(),
        r_big_al.quantity("C_Sigma_F").unwrap().epistemic.as_ref().unwrap().low.to_bits()
    );
    assert_eq!(
        r_narrow.quantity("E_C_Hz").unwrap().epistemic.as_ref().unwrap().up.to_bits(),
        r_big_al.quantity("E_C_Hz").unwrap().epistemic.as_ref().unwrap().up.to_bits()
    );

    // 4. Companion ems CircuitRom + dfs Prager-Synge bracket propagation + Roache GCI error_budget
    let mut rom_res = make_transmon_resonator_rom(16.0e9, None, 75.0e-15);
    let dummy_mesh = Mesh::new(3, vec![], vec![], vec![], vec![], vec![]);
    let comp_rep = uq_cqed_with_companion_and_params(
        None,
        &rom_res,
        Some((&dummy_mesh, &cap_wide)),
        &json!({
            "basis": "normal_fock",
            "fock_trunc": [6, 4],
            "n_samples": 128,
            "seed": 5,
            "ground": "c0",
            "qubits": [{ "name": "Q1_JJ", "pads": ["c2", "c3"], "E_J_Hz": 16.0e9 }],
            "uncertainties": { "E_J_Hz": { "rel_tol": 0.02 } }
        }),
    )
    .unwrap();
    let f01_joint = comp_rep.quantity("f01_Hz").unwrap();
    let f01_ep = f01_joint.epistemic.as_ref().unwrap();
    assert_eq!(f01_ep.source, "prager_synge_capacitance");
    assert!(f01_ep.low < f01_joint.nominal && f01_joint.nominal < f01_ep.up);
    assert!(f01_joint.aleatoric.as_ref().unwrap().linear_sigma > 0.0);

    rom_res.error_budget = Some(NumericalErrorBudget {
        algebraic_rel: 1e-11,
        spatial_discretization: Some(SpatialErrorEstimate {
            method: "roache_gci".into(),
            rel_bound: 0.015,
            is_epistemic_bound: true,
            ..Default::default()
        }),
        ..Default::default()
    });
    let gci_rep = uq_cqed_with_params(
        &rom_res,
        &json!({
            "basis": "normal_fock",
            "fock_trunc": [6, 4],
            "n_samples": 128,
            "seed": 5,
            "uncertainties": { "E_J_Hz": { "rel_tol": 0.02 } }
        }),
    )
    .unwrap();
    let gci_f01 = gci_rep.quantity("f01_Hz").unwrap();
    assert_eq!(gci_f01.epistemic.as_ref().unwrap().source, "error_budget:roache_gci");
    assert!((gci_f01.epistemic.as_ref().unwrap().rel_half_width - 0.015).abs() < 1e-12);
    assert_eq!(
        gci_f01.aleatoric.as_ref().unwrap().linear_sigma.to_bits(),
        f01_joint.aleatoric.as_ref().unwrap().linear_sigma.to_bits()
    );
}

#[test]
fn uq_cqed_tls_2d_and_3d_thickness_gram_zero_solve() {
    // 1. 2D TLS solve once, then propagate interface & tan_delta UQ with zero extra FEM solves
    let (mesh2, mut res2) = solve_fast_2d_cpw_tls();
    res2.error_budget = Some(NumericalErrorBudget {
        algebraic_rel: 1e-10,
        spatial_discretization: Some(SpatialErrorEstimate {
            method: "prager_synge".into(),
            rel_bound: 0.04,
            is_epistemic_bound: true,
            ..Default::default()
        }),
        ..Default::default()
    });
    let tls_nom = tls_participation_with_params(&mesh2, &res2, &Value::Null).unwrap();
    let uq2 = uq_cqed_with_mesh_and_params(
        &mesh2,
        &res2,
        &json!({
            "seed": 77,
            "n_samples": 800,
            "uncertainties": {
                "tan_delta": { "rel_tol": 0.10 },
                "t_ms": { "rel_tol": 0.05 }
            }
        }),
    )
    .unwrap();
    let lt = uq2.quantity("loss_tan_total").unwrap();
    assert!((lt.nominal - tls_nom.loss_tan_total).abs() / tls_nom.loss_tan_total < 1e-12);
    let lt_al = lt.aleatoric.as_ref().unwrap();
    let lt_ep = lt.epistemic.as_ref().unwrap();
    assert_eq!(lt_ep.source, "error_budget:prager_synge");
    assert!((lt_ep.rel_half_width - 0.04).abs() < 1e-12);
    assert!((lt_al.linear_sigma - lt_al.mc_sigma).abs() / lt_al.linear_sigma < 0.12);
    let p_sum = lt_al.variance_contributions["tan_delta_scale"] + lt_al.variance_contributions["t_ms"];
    assert!((p_sum - 1.0).abs() < 1e-12);
    let qtls = uq2.quantity("Q_TLS").unwrap();
    assert!((qtls.nominal - tls_nom.q_tls).abs() / tls_nom.q_tls < 1e-12);
    assert!(qtls.aleatoric.as_ref().unwrap().rel_linear_sigma > 0.10);

    // 2. 3D CPW sheet: verify cached PreparedThicknessGram + PreparedTls3d propagate 1000 samples fast
    let (mesh3, res3) = solve_fast_3d_cpw_sheet();
    let uq3 = uq_cqed_with_mesh_and_params(
        &mesh3,
        &res3,
        &json!({
            "seed": 1234,
            "n_samples": 600,
            "metal_thickness_m": 150.0e-9,
            "qubits": [{ "name": "Q_trace", "pads": ["trace"], "E_J_Hz": 15.0e9 }],
            "d": 150.0e-9,
            "trench_depth_m": 40.0e-9,
            "uncertainties": {
                "metal_thickness_m": { "rel_tol": 0.08 },
                "d": { "rel_tol": 0.10 },
                "trench_depth_m": 10.0e-9,
                "tan_delta_ms": { "rel_tol": 0.12 }
            }
        }),
    )
    .unwrap();
    assert!(uq3.elapsed_s < 2.0, "3D zero-solve UQ took {} s", uq3.elapsed_s);
    assert!(uq3.quantity("C_Sigma_F").unwrap().aleatoric.as_ref().unwrap().linear_sigma > 0.0);
    assert!(uq3.quantity("f01_Hz").unwrap().aleatoric.as_ref().unwrap().linear_sigma > 0.0);
    assert!(uq3.quantity("Q_TLS").unwrap().aleatoric.as_ref().unwrap().linear_sigma > 0.0);
}

#[test]
fn uq_cqed_validation_and_refusal_cases() {
    let res = make_transmon_resonator_rom(16.0e9, None, 75.0e-15);
    assert!(uq_cqed_with_params(
        &res,
        &json!({ "uncertainties": { "E_J_Hz": { "rel_tol": 0.01 }, "L_J_H": { "rel_tol": 0.01 } } })
    )
    .is_err());

    let (mesh2, res2) = solve_fast_2d_cpw_tls();
    assert!(uq_cqed_with_mesh_and_params(
        &mesh2,
        &res2,
        &json!({ "uncertainties": { "tan_delta": { "rel_tol": 0.1 }, "tan_delta_ma": { "rel_tol": 0.1 } } })
    )
    .is_err());

    assert!(uq_cqed_with_params(&SolveResult::default(), &Value::Null).is_err());
}
