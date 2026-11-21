// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Proportional local sizing (`scale_sizes`) and mode-signature tracking (`track` with `region`/`min_share`, `component`/`min_component`,
//! `max_layer_fraction`, and `f`/`n_b`/`n_eff` targets) in `convergence_check` (FR-6).

use geom::model::FeatureList;
use geom::sim::SizeRule;
use nyx_engine::convergence::Limits;
use nyx_engine::{Engine, FemcoreMesher, MeshSettings, MeshedModel, Mesher};
use nyx_physics::{
    BcDesc, Ctx, Dims, Error, FieldDesc, FieldLink, Mesh, ParamDesc, Physics, ProblemSpec, PropDesc, Registry, Scalar,
    Sizing, Solution, SolveResult, StudyDesc,
};
use serde_json::{json, Value};
use std::sync::{Arc, Mutex};

#[derive(Clone)]
struct RecordedMesh {
    rules: Vec<SizeRule>,
    total_cells: usize,
    core_cells: usize,
}

struct RecordingMesher {
    inner: FemcoreMesher,
    seen: Mutex<Vec<RecordedMesh>>,
}

impl RecordingMesher {
    fn new() -> Arc<Self> {
        Arc::new(Self { inner: FemcoreMesher, seen: Mutex::new(Vec::new()) })
    }

    fn take(&self) -> Vec<RecordedMesh> {
        std::mem::take(&mut *self.seen.lock().unwrap())
    }
}

impl Mesher for RecordingMesher {
    fn mesh(
        &self,
        geometry: &FeatureList,
        settings: &MeshSettings,
        sizing: Option<&Sizing>,
    ) -> Result<MeshedModel, Error> {
        let m = self.inner.mesh(geometry, settings, sizing)?;
        let si_id = m.regions().into_iter().find(|(_, name)| name == "Si").map(|(id, _)| id);
        let sm = m.solver_mesh();
        let core_cells = match si_id {
            Some(id) => sm.cell_labels.iter().filter(|&&l| l == id).count(),
            None => 0,
        };
        self.seen.lock().unwrap().push(RecordedMesh {
            rules: m.evaluated().size_rules().to_vec(),
            total_cells: m.mesh().n_cells(),
            core_cells,
        });
        Ok(m)
    }
}

fn c(e: &mut Engine, cmd: &str, args: Value) -> Value {
    e.call(cmd, &args).unwrap_or_else(|er| panic!("{cmd}({args}): {er}"))
}

#[test]
fn scale_sizes_refines_numeric_and_text_size_rules_preserves_project_and_sets_checked() {
    assert!(Limits::default().scale_sizes, "scale_sizes must default to true");

    // 1. Numeric `Expr::Num` rules on `h`, `max_deviation`, and `min_layers`, plus a suppressed `add_size`
    let probe = RecordingMesher::new();
    let mut e = Engine::new();
    e.set_mesher(probe.clone());
    c(&mut e, "project_new", json!({ "title": "soi_local_sizes_num", "dim": 2, "h": 8e-7 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [-2.0e-6, -2.0e-6], "max": [2.0e-6, 0.0] }));
    c(&mut e, "add_rect", json!({ "name": "clad", "min": [-2.0e-6, 0.0], "max": [2.0e-6, 2.0e-6] }));
    c(
        &mut e,
        "add_rounded_rect",
        json!({ "name": "core", "min": [-2.5e-7, 2.0e-8], "max": [2.5e-7, 2.7e-7], "r": 4e-8 }),
    );
    c(
        &mut e,
        "add_overlay",
        json!({
            "name": "regions",
            "regions": [
                { "name": "BOX", "input": "box" },
                { "name": "clad", "input": "clad" },
                { "name": "Si", "input": "core" }
            ]
        }),
    );
    // Active numeric size rule on the core (`h = 50 nm` << global `800 nm`, `max_deviation = 6 nm`, `min_layers = 3`)
    c(
        &mut e,
        "add_size",
        json!({ "name": "size_core", "sel": "region(Si)", "h": 5.0e-8, "max_deviation": 6.0e-9, "min_layers": 3.0 }),
    );
    // Suppressed size rule: must stay untouched
    c(&mut e, "add_size", json!({ "name": "size_ignored", "sel": "region(BOX)", "h": 2.0e-8 }));
    c(&mut e, "suppress_feature", json!({ "name": "size_ignored", "suppressed": true }));
    c(
        &mut e,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "modes",
            "order": 2,
            "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.5, "n_modes": 1 }
        }),
    );
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "BOX", "material": "silica" }));
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" }));
    c(&mut e, "set_material", json!({ "physics": "ems", "region": "Si", "material": "silicon" }));

    let project_before = e.project().to_json();

    // Verify `plan_with` for both `scale_sizes = false` and `scale_sizes = true` at factor 0.6
    let plan_unscaled = e.convergence_job().plan_with("ems", 0.6, false).expect("plan_with false");
    let meshes_unscaled = probe.take();
    assert_eq!(meshes_unscaled.len(), 2);
    let (coarse_m, fine_unscaled_m) = (&meshes_unscaled[0], &meshes_unscaled[1]);
    assert_eq!(plan_unscaled["scale_sizes"], json!(false));
    assert_eq!(coarse_m.rules.len(), 1, "suppressed rule excluded");
    assert_eq!(fine_unscaled_m.rules[0].h, Some(5.0e-8), "scale_sizes=false leaves h untouched");
    assert_eq!(fine_unscaled_m.rules[0].max_dev, Some(6.0e-9), "scale_sizes=false leaves max_deviation untouched");
    assert_eq!(fine_unscaled_m.rules[0].min_layers, Some(3.0), "scale_sizes=false leaves min_layers untouched");

    let plan_scaled = e.convergence_job().plan("ems", 0.6).expect("default plan scales sizes");
    let meshes_scaled = probe.take();
    assert_eq!(meshes_scaled.len(), 2);
    let fine_scaled_m = &meshes_scaled[1];
    assert_eq!(plan_scaled["scale_sizes"], json!(true));
    assert!((fine_scaled_m.rules[0].h.unwrap() - 5.0e-8 * 0.6).abs() < 1e-20);
    assert!((fine_scaled_m.rules[0].max_dev.unwrap() - 6.0e-9 * 0.6).abs() < 1e-21);
    // (3.0 / 0.6).ceil() == 5.0
    assert_eq!(fine_scaled_m.rules[0].min_layers, Some((3.0_f64 / 0.6).ceil()));
    assert!(
        fine_scaled_m.core_cells >= 2 * coarse_m.core_cells,
        "scale_sizes=true refines the local Si core ({} -> {} core cells)",
        coarse_m.core_cells,
        fine_scaled_m.core_cells
    );
    assert!(
        fine_scaled_m.core_cells >= 2 * fine_unscaled_m.core_cells,
        "scale_sizes=false did NOT refine the local Si core ({} unscaled vs {} scaled core cells, coarse={})",
        fine_unscaled_m.core_cells,
        fine_scaled_m.core_cells,
        coarse_m.core_cells
    );

    // Solve on the engine, then check with `scale_sizes: false` and `scale_sizes: true`
    c(&mut e, "solve_run", json!({}));
    probe.take();

    let check_false = c(&mut e, "convergence_check", json!({ "physics": "ems", "factor": 0.6, "scale_sizes": false }));
    assert_eq!(check_false["scale_sizes"], json!(false));
    assert_eq!(check_false["cells"], plan_unscaled["cells"], "plan_with(false) matches check(false) cells");
    probe.take();

    let check_true = c(&mut e, "convergence_check", json!({ "physics": "ems", "factor": 0.6, "tolerance": 4e-2 }));
    assert_eq!(check_true["scale_sizes"], json!(true), "default check has scale_sizes=true");
    assert_eq!(check_true["cells"], plan_scaled["cells"], "plan_with(true) matches check(true) cells");
    assert_eq!(check_true["converged"], json!(true), "{check_true}");
    assert_eq!(check_true["validated"]["state"], json!("checked"), "{check_true}");
    assert_eq!(c(&mut e, "result_info", json!({ "physics": "ems" }))[0]["validated"]["state"], json!("checked"));
    assert_eq!(e.project().to_json(), project_before, "active project must never be mutated");

    // 2. Parametric `Expr::Text` rules (`h`, `max_deviation`, `min_layers` with expressions needing parentheses)
    let probe_txt = RecordingMesher::new();
    let mut et = Engine::new();
    et.set_mesher(probe_txt.clone());
    c(&mut et, "project_new", json!({ "title": "soi_local_sizes_expr", "dim": 2, "h": 8e-7 }));
    c(&mut et, "set_param", json!({ "name": "h0", "value": 3.0e-8 }));
    c(&mut et, "set_param", json!({ "name": "d0", "value": 4.0e-9 }));
    c(&mut et, "set_param", json!({ "name": "nl", "value": 2.0 }));
    c(&mut et, "add_rect", json!({ "name": "box", "min": [-2.0e-6, -2.0e-6], "max": [2.0e-6, 0.0] }));
    c(&mut et, "add_rect", json!({ "name": "clad", "min": [-2.0e-6, 0.0], "max": [2.0e-6, 2.0e-6] }));
    c(
        &mut et,
        "add_rounded_rect",
        json!({ "name": "core", "min": [-2.5e-7, 2.0e-8], "max": [2.5e-7, 2.7e-7], "r": 4e-8 }),
    );
    c(
        &mut et,
        "add_overlay",
        json!({
            "name": "regions",
            "regions": [
                { "name": "BOX", "input": "box" },
                { "name": "clad", "input": "clad" },
                { "name": "Si", "input": "core" }
            ]
        }),
    );
    // Sum expressions (`h0 + 2e-8`) verify `({s}) * factor` parenthesizes properly instead of scaling only the second term!
    c(
        &mut et,
        "add_size",
        json!({
            "name": "size_core_expr",
            "sel": "region(Si)",
            "h": "h0 + 2e-8",
            "max_deviation": "d0 + 2e-9",
            "min_layers": "nl + 1"
        }),
    );
    c(
        &mut et,
        "add_physics",
        json!({
            "physics": "ems",
            "study": "modes",
            "order": 2,
            "study_params": { "freq": 299792458.0 / 1.55e-6, "n_target": 2.5, "n_modes": 1 }
        }),
    );
    c(&mut et, "set_material", json!({ "physics": "ems", "region": "BOX", "material": "silica" }));
    c(&mut et, "set_material", json!({ "physics": "ems", "region": "clad", "material": "silica" }));
    c(&mut et, "set_material", json!({ "physics": "ems", "region": "Si", "material": "silicon" }));

    let expr_before = et.project().to_json();
    c(&mut et, "solve_run", json!({}));
    probe_txt.take();

    let r_expr = c(
        &mut et,
        "convergence_check",
        json!({ "physics": "ems", "factor": 0.5, "scale_sizes": true, "tolerance": 4e-2 }),
    );
    let seen_expr = probe_txt.take();
    // Coarse solve was reused from `solve_run`, so only the refined mesh was built!
    assert_eq!(seen_expr.len(), 1);
    let fine_rule = &seen_expr[0].rules[0];
    assert!((fine_rule.h.unwrap() - 5.0e-8 * 0.5).abs() < 1e-20, "({:?})", fine_rule.h);
    assert!((fine_rule.max_dev.unwrap() - 6.0e-9 * 0.5).abs() < 1e-21, "({:?})", fine_rule.max_dev);
    assert!((fine_rule.min_layers.unwrap() - 3.0 / 0.5).abs() < 1e-12, "({:?})", fine_rule.min_layers);
    assert!(seen_expr[0].core_cells > fine_unscaled_m.core_cells * 2);
    assert_eq!(r_expr["converged"], json!(true), "{r_expr}");
    assert_eq!(r_expr["validated"]["state"], json!("checked"), "{r_expr}");
    assert_eq!(et.project().to_json(), expr_before, "active project with Expr::Text must never be mutated");
    let _ = seen_expr[0].total_cells;
}

#[test]
fn track_signature_excludes_cladding_and_pml_modes_on_guided_waveguide() {
    // 2D elastic waveguide (`sms`, `study: "modes"`) with slow-wave core `Si` embedded in `BOX` + `clad` + absorbing PML:
    // emits `share in <region> N`, `lateral/vertical/longitudinal share N`, and `layer fraction N`.
    let mut e = Engine::new();
    c(&mut e, "project_new", json!({ "title": "guided_acoustic_core", "dim": 2, "h": 2.5e-7 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [-1.8e-6, -1.5e-6], "max": [1.8e-6, 0.0] }));
    c(&mut e, "add_rect", json!({ "name": "clad", "min": [-1.8e-6, 0.0], "max": [1.8e-6, 1.5e-6] }));
    c(&mut e, "add_rect", json!({ "name": "core", "min": [-2.5e-7, -1.5e-7], "max": [2.5e-7, 1.5e-7] }));
    c(&mut e, "add_rect", json!({ "name": "win", "min": [-1.4e-6, -1.1e-6], "max": [1.4e-6, 1.1e-6] }));
    c(
        &mut e,
        "add_overlay",
        json!({
            "name": "regions",
            "regions": [
                { "name": "BOX", "input": "box" },
                { "name": "clad", "input": "clad" },
                { "name": "Si", "input": "core" }
            ]
        }),
    );
    c(
        &mut e,
        "add_pml",
        json!({
            "name": "pml",
            "window": "win",
            "sides": ["XMin", "XMax", "ZMin", "ZMax"],
            "thickness": 4.0e-7,
            "layers": 2.0,
            "grading": 2.0,
            "regions": ["BOX", "clad"]
        }),
    );
    c(&mut e, "add_size", json!({ "name": "core_size", "sel": "region(Si)", "h": 8.0e-8 }));
    c(
        &mut e,
        "add_physics",
        json!({
            "physics": "sms",
            "study": "modes",
            "order": 2,
            "study_params": { "wavevector": 6.0e6, "f_search": 3.5e9, "n_modes": 12 }
        }),
    );
    // Stiff fast-wave cladding (`v_s ~ 5786 m/s`) surrounding a slower guided core (`v_s ~ 2000 m/s`)
    let stiff = json!({ "E": 180e9, "nu": 0.22, "rho": 2200.0 });
    let slow_core = json!({ "E": 40e9, "nu": 0.25, "rho": 4000.0 });
    c(&mut e, "set_material", json!({ "physics": "sms", "region": "BOX*", "material": stiff }));
    c(&mut e, "set_material", json!({ "physics": "sms", "region": "clad*", "material": stiff }));
    c(&mut e, "set_material", json!({ "physics": "sms", "region": "Si", "material": slow_core }));
    c(&mut e, "add_pml_conditions", json!({ "physics": "sms", "strength": 2.5 }));

    let before = e.project().to_json();
    c(&mut e, "solve_run", json!({}));
    let solved = e.result("sms").expect("sms result").0.clone();
    let num = |name: &str| -> f64 {
        solved
            .scalars
            .iter()
            .find_map(|(k, v)| match v {
                Scalar::Num(x) if k == name => Some(*x),
                _ => None,
            })
            .unwrap_or_else(|| panic!("missing scalar {name}"))
    };

    // Locate a cladding/layer mode (`share in Si < 0.2`) and a guided core mode (`share in Si > 0.6`, `layer fraction < 0.05`)
    let clad_idx = (0..solved.solutions.len())
        .find(|&i| num(&format!("share in Si {}", i + 1)) < 0.2)
        .expect("at least one cladding/layer mode near f_search");
    let guided_idx = (0..solved.solutions.len())
        .find(|&i| num(&format!("share in Si {}", i + 1)) > 0.6 && num(&format!("layer fraction {}", i + 1)) < 0.05)
        .expect("at least one guided core mode");
    let (dom_comp, dom_share) = ["lateral", "vertical", "longitudinal"]
        .into_iter()
        .map(|comp| (comp, num(&format!("{comp} share {}", guided_idx + 1))))
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap();
    let min_comp = (dom_share * 0.8).min(0.95);
    let f_clad = solved.solutions[clad_idx].eigenvalue.unwrap()[0];

    // Without mode-signature filtering, targeting `f_clad` locks onto the cladding mode (`clad_idx`)
    let plain = c(&mut e, "convergence_check", json!({ "physics": "sms", "factor": 0.7, "track": { "f": f_clad } }));
    assert_eq!(plain["pairs"][0]["coarse"].as_u64(), Some(clad_idx as u64), "{plain}");
    assert_eq!(plain["track"], json!({ "f": f_clad }), "backwards-compatible track JSON for plain {{f}}");

    // With mode-signature filtering, even though `f` is set right on `f_clad` (cladding mode sits at distance 0!),
    // cladding/PML modes are excluded and the guided core mode is paired and converges!
    let tracked = c(
        &mut e,
        "convergence_check",
        json!({
            "physics": "sms",
            "factor": 0.7,
            "tolerance": 0.10,
            "track": {
                "f": f_clad,
                "region": "Si",
                "min_share": 0.50,
                "component": dom_comp,
                "min_component": min_comp,
                "max_layer_fraction": 0.05
            }
        }),
    );
    assert_eq!(tracked["match"], json!("track"));
    assert_eq!(tracked["matching"], json!("track"));
    assert_eq!(tracked["pairs"].as_array().unwrap().len(), 1);
    let paired_coarse = tracked["pairs"][0]["coarse"].as_u64().unwrap() as usize;
    assert!(
        num(&format!("share in Si {}", paired_coarse + 1)) >= 0.50
            && num(&format!("{dom_comp} share {}", paired_coarse + 1)) >= min_comp
            && num(&format!("layer fraction {}", paired_coarse + 1)) <= 0.05,
        "paired coarse mode {paired_coarse} must satisfy all signature filters (expected guided e.g. {guided_idx})"
    );
    assert_ne!(paired_coarse, clad_idx, "cladding mode at exact target frequency was excluded");
    let f_rel = tracked["rows"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["quantity"] == format!("eigenvalue[{paired_coarse}].re"))
        .unwrap()["relative"]
        .as_f64()
        .unwrap();
    assert!(f_rel < 2e-3, "guided mode real frequency converged to < 0.2% (got {f_rel:.3e})");
    assert_eq!(tracked["converged"], json!(true), "guided core mode converged cleanly: {tracked}");
    assert_eq!(tracked["validated"]["state"], json!("checked"), "{tracked}");
    assert_eq!(tracked["validated"]["check"]["matching"], json!("track"));
    assert!(tracked["validated"]["check"]["observable"].as_str().unwrap().contains("tracked mode"));
    assert_eq!(e.project().to_json(), before, "project untouched");
}

/// Controlled mock physics whose coarse vs fine solve differs by mesh cell count: proves candidate filtering runs on BOTH
/// coarse and fine `SolveResult`s before pairing, tests `n_b`/`n_eff`/target-less `pair_nearest`, and checks diagnostics.
struct MockModesPhysics;

impl Physics for MockModesPhysics {
    fn id(&self) -> &'static str {
        "tst"
    }
    fn name(&self) -> &'static str {
        "Mock mode-signature physics"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc { kind: "modes".into(), name: "Modes".into(), doc: "".into(), dims: Dims::Both, params: vec![] }]
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
    fn sizing(&self, spec: &ProblemSpec, _: Option<&Arc<femcore::material::Library>>) -> Sizing {
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, _ctx: &mut Ctx) -> Result<SolveResult, Error> {
        let is_fine = mesh.n_cells() > 40;
        // Tuple: (f_hz, si_share, box_share, ex_share, ez_share, layer_frac)
        let rows: Vec<(f64, f64, f64, f64, f64, f64)> = if !is_fine {
            vec![
                // Mode 0: cladding mode sitting right at 194.40 THz
                (194.400e12, 0.06, 0.94, 0.90, 0.05, 0.02),
                // Mode 1: PML radiation mode with high core overlap attempt but huge layer_fraction
                (194.410e12, 0.62, 0.38, 0.88, 0.06, 0.55),
                // Mode 2: TM0 guided mode in Si
                (194.420e12, 0.74, 0.26, 0.04, 0.91, 0.01),
                // Mode 3: TE0 guided mode in Si (share 0.850 on coarse, 0.848 on fine to test fine-only filter boundary too)
                (194.450e12, 0.850, 0.150, 0.92, 0.03, 0.01),
            ]
        } else {
            vec![
                // Mode 0 (fine): cladding distractor right next to coarse TE0 (194.45001 THz vs coarse 194.450 THz!)
                // If fine candidates were NOT filtered first, nearest-eigenvalue to coarse TE0 would pick this cladding mode!
                (194.450_01e12, 0.05, 0.95, 0.89, 0.06, 0.02),
                // Mode 1 (fine): true refined TE0 guided mode in Si
                (194.450_50e12, 0.848, 0.152, 0.92, 0.03, 0.01),
                // Mode 2 (fine): true refined TM0 guided mode in Si
                (194.420_30e12, 0.741, 0.259, 0.04, 0.91, 0.01),
                // Mode 3 (fine): PML mode
                (194.412_00e12, 0.61, 0.39, 0.87, 0.07, 0.52),
            ]
        };
        let mut solutions = Vec::new();
        let mut scalars = Vec::new();
        for (idx, (f, si, bx, ex, ez, lf)) in rows.into_iter().enumerate() {
            let n = idx + 1;
            solutions.push(Solution {
                label: format!("mode {n}"),
                eigenvalue: Some([f, 0.0]),
                value: Some(f),
                ..Default::default()
            });
            scalars.push((format!("f re {n} (Hz)"), Scalar::Num(f)));
            scalars.push((format!("share in Si {n}"), Scalar::Num(si)));
            scalars.push((format!("share in BOX {n}"), Scalar::Num(bx)));
            scalars.push((format!("Ex share {n}"), Scalar::Num(ex)));
            scalars.push((format!("Ez share {n}"), Scalar::Num(ez)));
            scalars.push((format!("layer fraction {n}"), Scalar::Num(lf)));
        }
        Ok(SolveResult {
            physics: "tst".into(),
            study: "modes".into(),
            order: spec.order,
            dim: mesh.dim() as u8,
            solutions,
            scalars,
            log: vec![],
            seconds: 0.001,
            dofs: mesh.n_cells() * 3,
        })
    }
}

#[test]
fn track_filters_both_meshes_supports_nb_neff_targetless_pairing_and_validates_inputs() {
    let mut reg = Registry::new();
    reg.register(Arc::new(MockModesPhysics)).unwrap();
    let mut e = Engine::with(reg, Arc::new(FemcoreMesher));
    c(&mut e, "project_new", json!({ "title": "mock_track", "dim": 2, "h": 0.35 }));
    c(&mut e, "add_rect", json!({ "name": "box", "min": [0.0, 0.0], "max": [1.0, 1.0] }));
    c(&mut e, "add_overlay", json!({ "name": "regions", "regions": [{ "name": "Si", "input": "box" }] }));
    c(&mut e, "add_physics", json!({ "physics": "tst", "study": "modes", "order": 1 }));

    // 1. Coarse cladding Mode 0 sits right at target 194.40 THz AND fine cladding Mode 0 sits right at 194.45001 THz
    // (closer to coarse TE0 194.45 THz than fine TE0 194.4505 THz!). Filtering both meshes pairs coarse 3 -> fine 1!
    let te = c(
        &mut e,
        "convergence_check",
        json!({
            "physics": "tst",
            "factor": 0.5,
            "tolerance": 5e-2,
            "track": {
                "f": 194.400e12,
                "region": "Si",
                "min_share": 0.50,
                "component": "Ex",
                "min_component": 0.60,
                "max_layer_fraction": 0.10
            }
        }),
    );
    assert_eq!(te["pairs"].as_array().unwrap().len(), 1);
    assert_eq!(te["pairs"][0]["coarse"], json!(3), "coarse TE0 selected over cladding, PML, and TM0: {te}");
    assert_eq!(te["pairs"][0]["refined"], json!(1), "fine TE0 selected over closer fine cladding distractor 0: {te}");
    assert_eq!(te["converged"], json!(true), "{te}");

    // 2. `n_b` and `n_eff` target keys
    let tm_nb = c(
        &mut e,
        "convergence_check",
        json!({
            "physics": "tst",
            "factor": 0.5,
            "track": {
                "n_b": 194.400e12,
                "region": "Si",
                "min_share": 0.50,
                "component": "Ez",
                "min_component": 0.60,
                "max_layer_fraction": 0.10
            }
        }),
    );
    assert_eq!((tm_nb["pairs"][0]["coarse"].as_u64(), tm_nb["pairs"][0]["refined"].as_u64()), (Some(2), Some(2)));

    let neff_only =
        c(&mut e, "convergence_check", json!({ "physics": "tst", "factor": 0.5, "track": { "n_eff": 194.400e12 } }));
    assert_eq!(neff_only["pairs"][0]["coarse"], json!(0));

    // 3. Target-less `track` (`target: None`): runs `pair_nearest` restricted to filtered candidates (pairs both TM0 and TE0)
    let both_guided = c(
        &mut e,
        "convergence_check",
        json!({
            "physics": "tst",
            "factor": 0.5,
            "tolerance": 5e-2,
            "track": { "region": "Si", "min_share": 0.50, "max_layer_fraction": 0.10 }
        }),
    );
    let pairs: Vec<(u64, u64)> = both_guided["pairs"]
        .as_array()
        .unwrap()
        .iter()
        .map(|p| (p["coarse"].as_u64().unwrap(), p["refined"].as_u64().unwrap()))
        .collect();
    assert_eq!(pairs, vec![(2, 2), (3, 1)], "both guided modes paired, cladding and PML excluded: {both_guided}");
    assert_eq!(both_guided["unpaired"]["coarse"], json!([]));
    assert_eq!(both_guided["unpaired"]["refined"], json!([]));
    assert_eq!(both_guided["converged"], json!(true));

    // 4. Diagnostic errors when zero coarse modes or zero fine modes match the filter
    let err_coarse = e
        .call(
            "convergence_check",
            &json!({ "physics": "tst", "factor": 0.5, "track": { "region": "Si", "min_share": 0.95 } }),
        )
        .unwrap_err()
        .to_string();
    assert!(err_coarse.contains("present mesh") && err_coarse.contains("share in Si >= 0.95"), "{err_coarse}");

    // Coarse TE0 has `share in Si = 0.850 >= 0.849`, while fine TE0 has `0.848 < 0.849` -> zero fine modes match!
    let err_fine = e
        .call(
            "convergence_check",
            &json!({ "physics": "tst", "factor": 0.5, "track": { "region": "Si", "min_share": 0.849 } }),
        )
        .unwrap_err()
        .to_string();
    assert!(err_fine.contains("refined mesh") && err_fine.contains("share in Si >= 0.849"), "{err_fine}");

    // 5. Input validation errors
    for bad in [
        json!({ "track": {} }),
        json!({ "track": 123 }),
        json!({ "track": { "unknown_key": 1.0 } }),
        json!({ "track": { "f": 1.0, "n_eff": 2.0 } }),
        json!({ "track": { "f": 1.0, "n_b": 2.0 } }),
        json!({ "track": { "region": "Si" } }),
        json!({ "track": { "min_share": 0.5 } }),
        json!({ "track": { "region": "", "min_share": 0.5 } }),
        json!({ "track": { "region": "Si", "min_share": 1.2 } }),
        json!({ "track": { "region": "Si", "min_share": -0.1 } }),
        json!({ "track": { "component": "Ex" } }),
        json!({ "track": { "min_component": 0.5 } }),
        json!({ "track": { "component": "Ex", "min_component": 1.5 } }),
        json!({ "track": { "max_layer_fraction": -0.05 } }),
        json!({ "track": { "max_layer_fraction": 1.05 } }),
        json!({ "match": "index", "track": { "region": "Si", "min_share": 0.5 } }),
        json!({ "scale_sizes": "yes" }),
    ] {
        let mut a = json!({ "physics": "tst", "factor": 0.5 });
        for (k, v) in bad.as_object().unwrap() {
            a[k] = v.clone();
        }
        assert!(e.call("convergence_check", &a).is_err(), "expected validation failure for {a}");
    }
}
