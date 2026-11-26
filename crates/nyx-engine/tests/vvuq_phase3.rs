// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Phase 3 VVUQ verification tests covering §5.1, §5.3, §5.4, and §5.5 across `nyx-engine`:
//! - 3-grid Roache GCI refinement and inverted/explicit coarsening ladders
//! - `r_eff < 1.3` conservative bracket fallback and warning
//! - `ValidationTier` strict `weakness()` ordering, `BoundCertified` elevation, `NotAffordable` gate, and `SetProject` check clearing
//! - `apply_extrapolated: true` setting `extrapolated_scalars_only` in `NumericalErrorBudget`
//! - R11a structured `ModeScalarMeta` pairing with legacy string fallback
//! - Worker-isolated `convergence_check` through `nyx-worker`
//! - `ems::cavity` Maxwell residual AMR, `sms` corner-refining Hessian AMR (R6f), `retain_amr_data`, and `quarter_point` + `flux_bounds` gate
//! - `reduce`: P2 interior extrema refinement (`interior: true` and candidate reporting) and non-polynomial / high-degree quadrature warning

use nyx_engine::convergence::{encode_mode_scalar_meta, ModeScalarMeta};
use nyx_engine::validated::{Validation, ValidationTier};
use nyx_engine::{Engine, FemcoreMesher};
use nyx_physics::{
    BcDesc, Ctx, Dims, ElectroAmrData, Error, FieldDesc, FieldLink, Mesh, ParamDesc, Physics, ProblemSpec, PropDesc,
    Registry, Scalar, Sizing, Solution, SolveResult, StudyDesc, ELECTRO_AMR_KEY,
};
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;

fn example(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples")
        .join(format!("{name}.nyx.json"))
        .to_string_lossy()
        .into_owned()
}

fn worker_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_nyx-worker"))
}

fn c(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|err| panic!("{cmd}({args}) failed: {err}"))
}

#[test]
fn validation_tier_weakness_ordering_and_set_project_clears_checks() {
    assert_eq!(ValidationTier::BoundCertified.weakness(), ValidationTier::GciVerified.weakness());
    let chain = [
        ValidationTier::GciVerified,
        ValidationTier::TwoGridChecked,
        ValidationTier::ResidualChecked,
        ValidationTier::Exploratory,
        ValidationTier::NotAffordable,
        ValidationTier::NonConvergent,
    ];
    for w in chain.windows(2) {
        assert!(w[0].weakness() < w[1].weakness(), "{:?} must be stronger than {:?}", w[0], w[1]);
    }

    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("02_cantilever") }));
    c(&mut e, "solve_run", json!({}));
    c(&mut e, "convergence_check", json!({ "physics": "sms", "factor": 0.6, "tolerance": 0.5 }));
    assert_eq!(e.validated("sms").state, Validation::Checked);
    assert_eq!(e.validated("sms").tier, ValidationTier::TwoGridChecked);

    c(&mut e, "project_new", json!({ "title": "Fresh", "dim": 2 }));
    assert_eq!(e.validated("sms").state, Validation::Exploratory);
    assert_eq!(e.validated("sms").tier, ValidationTier::Exploratory);
}

#[test]
fn three_grid_roache_gci_coarsening_ladder_and_warnings() {
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("04_thermal_bar") }));
    c(&mut e, "set_mesh", json!({ "adapt_passes": 2 }));
    c(&mut e, "solve_run", json!({ "physics": "dfs" }));
    let present_cells = e.model().unwrap().mesh().n_cells() as u64;

    // 3-grid refinement with loose tolerance (> 0.05) on adapt_passes > 1 triggers both warnings and Roache GCI.
    let r = c(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.6, "grids": 3, "tolerance": 0.1 }));
    assert_eq!(r["grids"], json!(3), "{r}");
    let cc = r["cells"]["coarse"].as_u64().unwrap();
    let cm = r["cells"]["medium"].as_u64().unwrap();
    let cf = r["cells"]["refined"].as_u64().unwrap();
    assert!(cc < cm && cm < cf, "expected coarse < medium < refined: {r}");
    assert_eq!(cc, present_cells);
    assert!(r["gci"].is_object(), "{r}");
    assert!(r["r_eff"].as_f64().unwrap() >= 1.3, "{r}");
    assert_eq!(r["spatial_error_estimate"]["method"], "roache_gci_3grid", "{r}");

    let codes: Vec<&str> = r["warnings"].as_array().unwrap().iter().filter_map(|w| w["code"].as_str()).collect();
    assert!(codes.contains(&"adapt_passes_downgraded"), "{codes:?}");
    assert!(codes.contains(&"loose_tolerance"), "{codes:?}");

    let val = e.validated("dfs");
    assert_eq!(val.state, Validation::Checked, "{val:?}");
    assert_eq!(val.tier, ValidationTier::GciVerified, "{val:?}");

    // Inverted coarsening ladder (coarsen: true with f < 1.0 -> [1, 1/f, (1/f)^2]).
    let rc = c(
        &mut e,
        "convergence_check",
        json!({ "physics": "dfs", "factor": 0.6, "grids": 3, "coarsen": true, "tolerance": 0.25 }),
    );
    assert_eq!(rc["cells"]["refined"].as_u64().unwrap(), present_cells, "{rc}");
    assert!(
        rc["cells"]["refined"].as_u64().unwrap() > rc["cells"]["medium"].as_u64().unwrap()
            && rc["cells"]["medium"].as_u64().unwrap() > rc["cells"]["coarse"].as_u64().unwrap(),
        "{rc}"
    );

    // Explicit coarsen_factors [1.5, 2.25] also keeps present as finest.
    let rexp = c(
        &mut e,
        "convergence_check",
        json!({ "physics": "dfs", "factor": 0.6, "grids": 3, "coarsen_factors": [1.5, 2.25], "tolerance": 0.25 }),
    );
    assert_eq!(rexp["cells"]["refined"].as_u64().unwrap(), present_cells, "{rexp}");
    assert!(rexp["cells"]["refined"].as_u64().unwrap() > rexp["cells"]["coarse"].as_u64().unwrap());
}

#[test]
fn low_refinement_ratio_falls_back_to_conservative_bracket_and_withholds_checked() {
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("04_thermal_bar") }));
    c(&mut e, "solve_run", json!({ "physics": "dfs" }));

    // 2-grid check with factor = 0.9 gives r_eff ~= 1.11 < 1.3.
    let r2 = c(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.9, "tolerance": 0.5 }));
    assert!(r2["r_eff"].as_f64().unwrap() < 1.3, "{r2}");
    let codes2: Vec<&str> = r2["warnings"].as_array().unwrap().iter().filter_map(|w| w["code"].as_str()).collect();
    assert!(codes2.contains(&"low_refinement_ratio"), "{r2}");
    assert_eq!(r2["spatial_error_estimate"]["method"], "conservative_bracket", "{r2}");
    assert_eq!(e.validated("dfs").state, Validation::Exploratory);
    assert_eq!(e.validated("dfs").tier, ValidationTier::Exploratory);

    // 3-grid check with factor = 0.9 falls back to conservative bracket and flags low_refinement_fallback.
    let r3 = c(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.9, "grids": 3, "tolerance": 0.5 }));
    assert_eq!(r3["gci"]["low_refinement_fallback"], json!(true), "{r3}");
    assert_eq!(r3["spatial_error_estimate"]["method"], "conservative_bracket", "{r3}");
    assert_eq!(e.validated("dfs").state, Validation::Exploratory);
}

#[test]
fn prager_synge_flux_bounds_elevate_to_bound_certified_even_when_refinement_unaffordable() {
    let build_lbracket = |target_tol: f64| -> Engine {
        let mut e = Engine::new();
        c(&mut e, "project_new", json!({ "title": "BoundCert", "dim": 2, "h": 2.5e-5 }));
        c(
            &mut e,
            "add_polygons",
            json!({
                "name": "l_vac",
                "polygons": [[[0.0, 0.0], [2e-4, 0.0], [2e-4, 1e-4], [1e-4, 1e-4], [1e-4, 2e-4], [0.0, 2e-4]]]
            }),
        );
        c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Vac", "input": "l_vac" }] }));
        c(&mut e, "set_mesh", json!({ "target_tol": target_tol }));
        c(
            &mut e,
            "add_physics",
            json!({ "physics": "dfs", "study": "electrostatic", "order": 1, "study_params": { "flux_bounds": "midpoint" } }),
        );
        c(&mut e, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }));
        c(&mut e, "add_condition", json!({ "physics": "dfs", "on": "outer & plane(x = 0)", "kind": "ground" }));
        c(
            &mut e,
            "add_condition",
            json!({ "physics": "dfs", "on": "outer & plane(x = 2e-4)", "kind": "terminal", "params": { "name": "sig" } }),
        );
        c(&mut e, "solve_run", json!({}));
        let cells = e.model().unwrap().mesh().n_cells() as u64;
        // Current solve fits 4 GB (2.0 GB < 3.6 GB), while 0.6x refinement (~2.8x cells -> ~5.6 GB) is NotAffordable.
        e.set_memory_probe(Arc::new(|| Some((4_000_000_000, 4_000_000_000))));
        let row = |cells: u64, bytes: u64| nyx_engine::estimate::Row {
            chain: "dfs".into(),
            dim: 2,
            cells,
            dofs: cells,
            peak_working_set: bytes,
            peak_commit: bytes,
            seconds: 0.5,
        };
        e.history().add(row(cells / 2, 1_000_000_000));
        e.history().add(row(cells, 2_000_000_000));
        e
    };

    let e_ok = build_lbracket(0.5);
    let v_ok = e_ok.validated("dfs");
    assert_eq!(v_ok.state, Validation::Checked, "{v_ok:?}");
    assert_eq!(v_ok.tier, ValidationTier::BoundCertified, "{v_ok:?}");

    // Tightening target_tol below the Prager-Synge gap drops certification back to NotAffordable.
    let e_tight = build_lbracket(1e-12);
    let v_tight = e_tight.validated("dfs");
    assert_eq!(v_tight.state, Validation::NotAffordable, "{v_tight:?}");
    assert_eq!(v_tight.tier, ValidationTier::NotAffordable, "{v_tight:?}");
}

#[test]
fn apply_extrapolated_marks_extrapolated_scalars_only_in_error_budget() {
    let mut e = Engine::new();
    c(&mut e, "project_new", json!({ "title": "Extrap", "dim": 2, "h": 2.5e-5 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [-1e-4, -1e-4], "max": [1e-4, 1e-4] }));
    c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Vac", "input": "box" }] }));
    c(&mut e, "add_physics", json!({ "physics": "dfs", "study": "electrostatic", "order": 1 }));
    c(&mut e, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }));
    c(&mut e, "add_condition", json!({ "physics": "dfs", "on": "outer & plane(z = -1e-4)", "kind": "ground" }));
    c(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer & plane(z = 1e-4)", "kind": "terminal", "params": { "name": "top" } }),
    );
    c(&mut e, "solve_run", json!({}));

    let r = c(
        &mut e,
        "convergence_check",
        json!({
            "physics": "dfs",
            "factor": 0.5,
            "tolerance": 0.25,
            "extrapolate": "richardson_log_spd",
            "apply_extrapolated": true
        }),
    );
    assert_eq!(r["spatial_error_estimate"]["extrapolated_scalars_only"], json!(true), "{r}");
    let info = c(&mut e, "result_info", json!({ "physics": "dfs" }));
    assert_eq!(info[0]["error_budget"]["spatial_discretization"]["extrapolated_scalars_only"], json!(true), "{info}");
}

#[test]
fn worker_isolated_convergence_check_runs_through_nyx_worker() {
    let mut e = Engine::new();
    c(&mut e, "project_open", json!({ "path": example("04_thermal_bar") }));
    e.set_worker(worker_bin(), true);
    let r = c(&mut e, "convergence_check", json!({ "physics": "dfs", "factor": 0.6, "tolerance": 0.1 }));
    assert_eq!(r["worker_isolated"], json!(true), "{r}");
    assert_eq!(r["converged"], json!(true), "{r}");
}

#[test]
fn r11a_mode_scalar_metadata_takes_precedence_and_keeps_legacy_fallback() {
    struct SwappingModesPhysics;
    impl Physics for SwappingModesPhysics {
        fn id(&self) -> &'static str {
            "swp"
        }
        fn name(&self) -> &'static str {
            "Mode swap test"
        }
        fn studies(&self) -> Vec<StudyDesc> {
            vec![StudyDesc {
                kind: "modes".into(),
                name: "Modes".into(),
                doc: "".into(),
                dims: Dims::Two,
                params: vec![],
            }]
        }
        fn bc_kinds(&self) -> Vec<BcDesc> {
            vec![]
        }
        fn properties(&self) -> Vec<PropDesc> {
            vec![]
        }
        fn params(&self) -> Vec<ParamDesc> {
            vec![]
        }
        fn fields(&self, _: &str) -> Vec<FieldDesc> {
            vec![]
        }
        fn exports(&self) -> Vec<FieldLink> {
            vec![]
        }
        fn imports(&self) -> Vec<FieldLink> {
            vec![]
        }
        fn validate(&self, _: &ProblemSpec, _: &Mesh) -> Result<(), Error> {
            Ok(())
        }
        fn sizing(&self, _: &ProblemSpec, _: Option<&Arc<femcore::material::Library>>) -> Sizing {
            Sizing::default()
        }
        fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, _: &mut Ctx) -> Result<SolveResult, Error> {
            let coarse = mesh.n_cells() <= 100;
            // On refined mesh, mode order swaps while physical branch values converge cleanly.
            let sols = if coarse {
                vec![
                    Solution {
                        label: "m0".into(),
                        value: Some(2.0),
                        eigenvalue: Some([2.0, 0.0]),
                        fields: vec![],
                        cell_fields: vec![],
                    },
                    Solution {
                        label: "m1".into(),
                        value: Some(3.0),
                        eigenvalue: Some([3.0, 0.0]),
                        fields: vec![],
                        cell_fields: vec![],
                    },
                ]
            } else {
                vec![
                    Solution {
                        label: "m0".into(),
                        value: Some(3.001),
                        eigenvalue: Some([3.001, 0.0]),
                        fields: vec![],
                        cell_fields: vec![],
                    },
                    Solution {
                        label: "m1".into(),
                        value: Some(2.001),
                        eigenvalue: Some([2.001, 0.0]),
                        fields: vec![],
                        cell_fields: vec![],
                    },
                ]
            };
            let meta = vec![
                ModeScalarMeta {
                    name: "branch gamma".into(),
                    base: "branch gamma".into(),
                    unit: "(a.u.)".into(),
                    mode_index: Some(if coarse { 0 } else { 1 }),
                },
                ModeScalarMeta {
                    name: "summary count 2".into(),
                    base: "summary count".into(),
                    unit: String::new(),
                    mode_index: None,
                },
            ];
            let scalars = vec![
                ("branch gamma".into(), Scalar::Num(if coarse { 10.0 } else { 10.001 })),
                // Legacy trailing-digit scalar still parses via fallback when no metadata exists.
                ("legacy metric 1".into(), Scalar::Num(if coarse { 5.0 } else { 99.0 })),
                ("legacy metric 2".into(), Scalar::Num(if coarse { 99.0 } else { 5.001 })),
                ("summary count 2".into(), Scalar::Num(4.0)),
                encode_mode_scalar_meta(&meta),
            ];
            Ok(SolveResult {
                physics: "swp".into(),
                study: spec.study.kind.clone(),
                order: spec.order,
                dim: 2,
                dofs: mesh.n_cells(),
                seconds: 0.01,
                solutions: sols,
                scalars,
                log: vec![],
                diagnostics: None,
                error_budget: None,
            })
        }
    }

    let mut reg = Registry::new();
    reg.register(Arc::new(SwappingModesPhysics)).unwrap();
    let mut e = Engine::with(reg, Arc::new(FemcoreMesher));
    c(&mut e, "project_new", json!({ "title": "R11a", "dim": 2, "h": 0.25 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "R", "input": "box" }] }));
    c(&mut e, "add_physics", json!({ "physics": "swp", "study": "modes", "order": 1 }));
    c(&mut e, "solve_run", json!({}));

    let r = c(
        &mut e,
        "convergence_check",
        json!({ "physics": "swp", "factor": 0.4, "tolerance": 0.02, "match": "nearest" }),
    );
    assert_eq!(r["converged"], json!(true), "{r}");
    let row_names: Vec<&str> = r["rows"].as_array().unwrap().iter().filter_map(|x| x["quantity"].as_str()).collect();
    assert!(row_names.iter().any(|n| n.contains("branch gamma")), "{row_names:?}");
    assert!(row_names.iter().any(|n| n.contains("legacy metric")), "{row_names:?}");
    // Explicit `mode_index: None` prevents `summary count 2` from matching as mode 2.
    assert!(!row_names.iter().any(|n| n.contains("summary count")), "{row_names:?}");
}

#[test]
fn ems_cavity_maxwell_amr_and_sms_corner_amr_and_retain_amr_data() {
    // 1. 3D PEC box cavity Maxwell residual AMR (`femcore::adapt::maxwell`).
    let mut em = Engine::new();
    c(&mut em, "project_new", json!({ "title": "BoxCavityAMR", "dim": 3 }));
    c(&mut em, "add_box", json!({ "name": "cav", "min": [0.0, 0.0, 0.0], "max": [0.04, 0.03, 0.02] }));
    c(&mut em, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Air", "input": "cav" }] }));
    c(&mut em, "set_mesh", json!({ "h": 0.015, "adapt_passes": 2, "adapt_growth": 1.6 }));
    c(
        &mut em,
        "add_physics",
        json!({ "physics": "ems", "study": "cavity", "order": 1, "study_params": { "f_search": 6.0e9, "n_modes": 1 } }),
    );
    c(&mut em, "set_material", json!({ "physics": "ems", "region": "Air", "material": { "eps_r": 1.0, "mu_r": 1.0 } }));
    c(&mut em, "mesh_build", json!({}));
    let initial_cells = em.model().unwrap().mesh().n_cells();
    c(&mut em, "solve_run", json!({}));
    let (res, stale) = em.result("ems").unwrap();
    assert!(!stale);
    assert!(res.log.iter().any(|l| l.contains("Maxwell AMR pass 2/2")), "{:?}", res.log);
    assert!(em.model().unwrap().mesh().n_cells() > initial_cells, "expected Maxwell AMR growth");

    // 2. 2D L-bracket structural corner-refining AMR (`R6f` using `femcore::adapt::estimate`).
    let mut sm = Engine::new();
    c(&mut sm, "project_new", json!({ "title": "LBracketAMR", "dim": 2 }));
    c(
        &mut sm,
        "add_polygons",
        json!({
            "name": "l_shape",
            "polygons": [[[0.0, 0.0], [1.0, 0.0], [1.0, 0.25], [0.5, 0.25], [0.5, 0.5], [0.0, 0.5]]]
        }),
    );
    c(&mut sm, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Steel", "input": "l_shape" }] }));
    c(&mut sm, "set_mesh", json!({ "h": 0.12, "adapt_passes": 3, "adapt_growth": 1.6, "adapt_rel_tol": 1e-6 }));
    c(&mut sm, "add_physics", json!({ "physics": "sms", "study": "static", "order": 1 }));
    c(&mut sm, "set_physics_params", json!({ "physics": "sms", "params": { "plane": "stress" } }));
    c(
        &mut sm,
        "set_material",
        json!({ "physics": "sms", "region": "Steel", "material": { "E": 2.0e11, "nu": 0.3, "rho": 7850.0, "alpha": 1.2e-5 } }),
    );
    c(
        &mut sm,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 0)", "kind": "fixed", "params": { "ux": 0.0, "uy": 0.0, "uz": 0.0 } }),
    );
    c(
        &mut sm,
        "add_condition",
        json!({ "physics": "sms", "on": "outer & plane(x = 1)", "kind": "load", "params": { "traction": [0.0, 0.0, -1.0e5] } }),
    );
    c(&mut sm, "mesh_build", json!({}));
    let init_sm_mesh = sm.model().unwrap().solver_mesh();
    let init_sm_cells = init_sm_mesh.n_cells();
    let min_corner_h = |m: &femcore::mesh::SimplexMesh| -> f64 {
        let mut best = f64::INFINITY;
        for c in 0..m.n_cells() {
            let v = m.cell(c);
            let pts: [[f64; 2]; 3] = std::array::from_fn(|i| {
                let p = m.point(v[i] as usize);
                [p[0], p[1]]
            });
            let near_corner = pts.iter().any(|p| ((p[0] - 0.5).powi(2) + (p[1] - 0.25).powi(2)).sqrt() < 0.08);
            if near_corner {
                for &(i, j) in &[(0, 1), (1, 2), (2, 0)] {
                    best = best.min(((pts[i][0] - pts[j][0]).powi(2) + (pts[i][1] - pts[j][1]).powi(2)).sqrt());
                }
            }
        }
        best
    };
    let h_corner_before = min_corner_h(&init_sm_mesh);
    c(&mut sm, "solve_run", json!({}));
    let (sm_res, _) = sm.result("sms").unwrap();
    assert!(sm_res.log.iter().any(|l| l.contains("SMS AMR pass")), "{:?}", sm_res.log);
    let final_sm_mesh = sm.model().unwrap().solver_mesh();
    assert!(final_sm_mesh.n_cells() > init_sm_cells, "expected SMS AMR cell refinement");
    let h_corner_after = min_corner_h(&final_sm_mesh);
    assert!(
        h_corner_after < 0.8 * h_corner_before,
        "expected corner refinement near (0.5, 0.25): {h_corner_before:.4} -> {h_corner_after:.4}"
    );

    // 3. `retain_amr_data` opt-in retains ElectroAmrData, while `quarter_point` + `flux_bounds` remains refused.
    let build_dfs_cap = |qp: bool, flux_bounds: &str, adapt_passes: u64| -> Engine {
        let mut d = Engine::new();
        c(&mut d, "project_new", json!({ "title": "RetainAMR", "dim": 2, "h": 2.5e-5 }));
        c(&mut d, "add_rect", json!({ "name": "box", "min": [-1e-4, -1e-4], "max": [1e-4, 1e-4] }));
        c(&mut d, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Vac", "input": "box" }] }));
        c(&mut d, "set_mesh", json!({ "singular_edges": if qp { "quarter_point" } else { "none" } }));
        c(
            &mut d,
            "add_physics",
            json!({
                "physics": "dfs",
                "study": "electrostatic",
                "order": 1,
                "study_params": { "flux_bounds": flux_bounds, "adapt_passes": adapt_passes }
            }),
        );
        c(&mut d, "set_material", json!({ "physics": "dfs", "region": "Vac", "material": { "eps_r": 1.0 } }));
        c(&mut d, "add_condition", json!({ "physics": "dfs", "on": "outer & plane(z = -1e-4)", "kind": "ground" }));
        c(
            &mut d,
            "add_condition",
            json!({ "physics": "dfs", "on": "outer & plane(z = 1e-4)", "kind": "terminal", "params": { "name": "top" } }),
        );
        d
    };
    let mut d = build_dfs_cap(false, "none", 2);
    c(&mut d, "solve_run", json!({}));
    assert!(!d.result("dfs").unwrap().0.scalars.iter().any(|(k, _)| k == ELECTRO_AMR_KEY));

    c(&mut d, "solve_run", json!({ "retain_amr_data": true }));
    assert!(
        ElectroAmrData::of(d.result("dfs").unwrap().0).is_some(),
        "retain_amr_data: true should keep ElectroAmrData"
    );

    c(&mut d, "set_mesh", json!({ "retain_amr_data": true }));
    c(&mut d, "solve_run", json!({}));
    assert!(
        ElectroAmrData::of(d.result("dfs").unwrap().0).is_some(),
        "mesh.retain_amr_data: true should keep ElectroAmrData"
    );

    let mut d_qp = build_dfs_cap(true, "report", 1);
    let err = d_qp.call("solve_run", &json!({})).unwrap_err().to_string();
    assert!(err.contains("quarter_point"), "unexpected error: {err}");
}

#[test]
fn reduce_p2_interior_extremum_and_quadrature_warning() {
    struct InteriorPeakPhysics;
    impl Physics for InteriorPeakPhysics {
        fn id(&self) -> &'static str {
            "pik"
        }
        fn name(&self) -> &'static str {
            "Interior quadratic peak"
        }
        fn studies(&self) -> Vec<StudyDesc> {
            vec![StudyDesc {
                kind: "poly".into(),
                name: "Poly".into(),
                doc: "".into(),
                dims: Dims::Two,
                params: vec![],
            }]
        }
        fn bc_kinds(&self) -> Vec<BcDesc> {
            vec![]
        }
        fn properties(&self) -> Vec<PropDesc> {
            vec![]
        }
        fn params(&self) -> Vec<ParamDesc> {
            vec![]
        }
        fn fields(&self, _: &str) -> Vec<FieldDesc> {
            vec![]
        }
        fn exports(&self) -> Vec<FieldLink> {
            vec![]
        }
        fn imports(&self) -> Vec<FieldLink> {
            vec![]
        }
        fn validate(&self, _: &ProblemSpec, _: &Mesh) -> Result<(), Error> {
            Ok(())
        }
        fn sizing(&self, _: &ProblemSpec, _: Option<&Arc<femcore::material::Library>>) -> Sizing {
            Sizing::default()
        }
        fn solve(&self, mesh: &Mesh, _: &ProblemSpec, _: &mut Ctx) -> Result<SolveResult, Error> {
            let sm = mesh.discretise().map_err(Error::Invalid)?;
            let nloc = sm.nloc(2);
            let mut vals = Vec::with_capacity(sm.n_cells() * nloc);
            let peak = [0.37, 0.19];
            let eval_f = |p: &[f64]| 2.0 - ((p[0] - peak[0]).powi(2) + (p[1] - peak[1]).powi(2));
            for c in 0..sm.n_cells() {
                let v = sm.cell(c);
                for a in 0..nloc {
                    let pt = if a < 3 {
                        sm.point(v[a] as usize).to_vec()
                    } else {
                        let (i, j) = match a {
                            3 => (0, 1),
                            4 => (0, 2),
                            _ => (1, 2),
                        };
                        let p0 = sm.point(v[i] as usize);
                        let p1 = sm.point(v[j] as usize);
                        vec![0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])]
                    };
                    vals.push(eval_f(&pt));
                }
            }
            Ok(SolveResult {
                physics: "pik".into(),
                study: "poly".into(),
                order: 2,
                dim: 2,
                dofs: vals.len(),
                seconds: 0.0,
                solutions: vec![Solution {
                    label: "poly".into(),
                    value: None,
                    eigenvalue: None,
                    fields: vec![],
                    cell_fields: vec![nyx_physics::CellField {
                        name: "p".into(),
                        unit: String::new(),
                        nloc,
                        re: vals,
                        im: vec![],
                    }],
                }],
                scalars: vec![],
                log: vec![],
                diagnostics: None,
                error_budget: None,
            })
        }
    }

    let mut reg = Registry::new();
    reg.register(Arc::new(InteriorPeakPhysics)).unwrap();
    let mut e = Engine::with(reg, Arc::new(FemcoreMesher));
    c(&mut e, "project_new", json!({ "title": "PeakP2", "dim": 2, "h": 0.45 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "All", "input": "box" }] }));
    c(&mut e, "add_physics", json!({ "physics": "pik", "study": "poly", "order": 2 }));
    c(&mut e, "solve_run", json!({}));

    // Default keeps nodal value intact while surfacing `interior_extremum` + warning.
    let nodal_ans = c(&mut e, "extremum", json!({ "expr": "p", "physics": "pik" }));
    let nodal_v = nodal_ans["value"].as_f64().unwrap();
    assert!(nodal_v < 2.0 - 1e-4, "coarse nodal peak must miss (0.37, 0.19): {nodal_ans}");
    assert!((nodal_ans["interior_extremum"]["value"].as_f64().unwrap() - 2.0).abs() < 1e-10, "{nodal_ans}");
    assert_eq!(nodal_ans["warnings"][0]["code"], "p2_interior_extremum_exceeds_nodal", "{nodal_ans}");

    // `interior: true` directly returns the P2 Newton-refined interior extremum and peak coordinates.
    let int_ans = c(&mut e, "extremum", json!({ "expr": "p", "physics": "pik", "interior": true }));
    assert!((int_ans["value"].as_f64().unwrap() - 2.0).abs() < 1e-10, "{int_ans}");
    let pt = int_ans["point"].as_array().unwrap();
    assert!(
        (pt[0].as_f64().unwrap() - 0.37).abs() < 1e-4 && (pt[1].as_f64().unwrap() - 0.19).abs() < 1e-4,
        "{int_ans}"
    );

    // Quadrature warning: degree-2 polynomial `p` emits no warnings; non-poly/deg-6 formulas warn without double pass.
    let exact_int = c(&mut e, "integrate", json!({ "expr": "p", "physics": "pik" }));
    assert!(exact_int.get("warnings").is_none(), "{exact_int}");
    for expr in ["sqrt(p)", "p * p * p", "p / (p + 2.0)"] {
        let non_poly = c(&mut e, "integrate", json!({ "expr": expr, "physics": "pik" }));
        assert_eq!(non_poly["warnings"][0]["code"], "low_quadrature_order_for_formula", "{expr}: {non_poly}");
    }
}
