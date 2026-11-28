// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Contracts of `nyx-physics` (ADR 19, ADR 25):
//! * `problem_spec.v3.json`: the serde shape of the problem envelope, one spec per physics
//! * `solve_result.v2.json`: the serde shape of a result (the old `SolveResult` was not serialisable, R27)
//! * `physics_descriptors.v1.json`: the serde shape of the descriptors, one of each kind
//! * `physics_api` (tests only): the trait is object safe and has the listed methods; the registry rules; the error codes
//!
//! Rewrite an exemplar: `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-physics --test contracts <test name>`, with a `VERSIONS.md` line.
mod common;
use common::*;
use nyx_physics::*;
use serde_json::{json, Value};
use std::sync::Arc;

fn exemplar_specs() -> Vec<ProblemSpec> {
    let mut ems = ProblemSpec::new("ems", 2, "modes");
    ems.study.params = json!({ "freq": 1.934e14, "n_target": 2.4, "evanescent": false, "n_modes": 4 });
    ems.materials.insert("core".into(), MaterialAssign::Props([("eps_r".to_string(), 12.08)].into()));
    ems.materials.insert("clad".into(), MaterialAssign::Named("silica".into()));
    ems.bcs.push(BcSpec { on: "clad_pml".into(), kind: "pml".into(), params: json!({ "x": "plus", "strength": 3.0 }) });

    let mut sms = ProblemSpec::new("sms", 2, "eigenfrequency");
    sms.study.params = json!({ "n_modes": 6, "f_search": 1.0e9 });
    sms.materials.insert(
        "beam".into(),
        MaterialAssign::Props([("E".to_string(), 169e9), ("nu".to_string(), 0.22), ("rho".to_string(), 2329.0)].into()),
    );
    sms.bcs.push(BcSpec {
        on: "outer & plane(x = 0)".into(),
        kind: "fixed".into(),
        params: json!({ "ux": 0.0, "uy": 0.0, "uz": 0.0 }),
    });
    sms.params = json!({ "plane": "strain", "inactive": ["air"] });

    let mut dfs = ProblemSpec::new("dfs", 1, "heat_transient");
    dfs.study.params = json!({ "t_end": 2.0, "n_out": 10, "t0": "293.15 + 50*exp(-((x-0.5)^2)/0.01)" });
    dfs.bcs.push(BcSpec {
        on: "outer".into(),
        kind: "convection".into(),
        params: json!({ "h": 10.0, "t_ext": 293.15 }),
    });

    vec![dfs, sms, ems]
}

#[test]
fn problem_spec_v3_exemplar_is_stable() {
    let specs = exemplar_specs();
    check_exemplar("problem_spec.v3.json", &serde_json::to_value(&specs).unwrap());
    // and the checked-in text reads back to the same values
    let text = std::fs::read_to_string(contract_file("problem_spec.v3.json")).unwrap();
    let back: Vec<ProblemSpec> = serde_json::from_str(&text).unwrap();
    assert_eq!(back, specs);
}

#[test]
fn problem_spec_v3_shape_that_other_crates_read_without_types() {
    // the store, the MCP server and the window read specs as untyped JSON: these names are the contract
    let v = serde_json::to_value(&exemplar_specs()[2]).unwrap();
    for k in ["format", "physics", "order", "study", "materials", "bcs"] {
        assert!(v.get(k).is_some(), "ProblemSpec lost the key {k}");
    }
    assert_eq!(v["study"]["kind"], json!("modes"));
    assert_eq!(v["materials"]["clad"], json!("silica"));
    assert_eq!(v["materials"]["core"]["eps_r"], json!(12.08));
    assert_eq!(v["bcs"][0]["on"], json!("clad_pml"));
    // the key is written in every problem (ADR 63)
    assert_eq!(v["format"], json!(1));
}

#[test]
fn problem_spec_v3_reading_is_strict_about_the_physics_and_lenient_about_absent_options() {
    // Documented, not endorsed (R20): unknown top-level keys are ignored by serde; the descriptors (`validate_spec`) are what refuse
    // unknown parameters inside `study.params`, `bcs[].params` and `params`.
    let minimal: ProblemSpec =
        serde_json::from_value(json!({ "physics": "dfs", "order": 1, "study": { "kind": "electrostatic" } })).unwrap();
    assert!(minimal.materials.is_empty() && minimal.bcs.is_empty() && minimal.params.is_null());
    // a problem given to a command without the key is read as the current format; a project file must carry it (`Project::from_json`)
    assert_eq!(minimal.format, PROBLEM_FORMAT);
    assert!(serde_json::from_value::<ProblemSpec>(json!({ "order": 1 })).is_err());
}

#[test]
fn solve_result_v2_exemplar_is_stable() {
    let r = SolveResult {
        physics: "ems".into(),
        study: "modes".into(),
        order: 2,
        dim: 2,
        dofs: 12,
        seconds: 0.5,
        solutions: vec![Solution {
            label: "mode 1".into(),
            fields: vec![Field {
                name: "E".into(),
                unit: "a.u.".into(),
                components: vec!["x".into(), "z".into()],
                values: vec![0.0, 1.0, 0.5, 0.25],
            }],
            value: Some(1.0),
            eigenvalue: Some([2.4, -1e-6]),
            cell_fields: vec![CellField {
                name: "Ex".into(),
                unit: "a.u.".into(),
                nloc: 3,
                re: vec![1.0, 2.0, 3.0],
                im: vec![],
            }],
        }],
        scalars: vec![("neff 1".into(), Scalar::Num(2.4)), ("polarization 1".into(), Scalar::Text("TE".into()))],
        log: vec!["solving ems on a 2D mesh".into()],
        diagnostics: Some(nyx_physics::SolverDiagnostics {
            solver_kind: "shift_invert_complex".into(),
            converged: true,
            iterations: Some(12),
            restarts: Some(1),
            max_algebraic_residual: 1e-12,
            per_mode_residuals: vec![1e-12],
            rcond_estimate: None,
            cholesky_fallbacks: 0,
            warnings: vec![],
        }),
        error_budget: Some(nyx_physics::NumericalErrorBudget {
            algebraic_rel: 1e-12,
            spatial_discretization: Some(nyx_physics::SpatialErrorEstimate {
                method: "two_grid_delta".into(),
                rel_bound: 1e-4,
                is_epistemic_bound: false,
                observed_order: None,
                refinement_ratio: Some(1.6),
                safety_factor: None,
                extrapolated_scalars_only: None,
            }),
            temporal_rel: None,
            model_truncation_rel: None,
            combined_num_uncertainty_rel: Some(1e-4),
            alg_below_disc_ratio: Some(1e-8),
        }),
    };
    check_exemplar("solve_result.v2.json", &serde_json::to_value(&r).unwrap());
    let back: SolveResult = serde_json::from_value(serde_json::to_value(&r).unwrap()).unwrap();
    assert_eq!(back, r);
    // Legacy v2.0 JSON omitting diagnostics & error_budget still deserializes cleanly
    let mut legacy = serde_json::to_value(&r).unwrap();
    let obj = legacy.as_object_mut().unwrap();
    obj.remove("diagnostics");
    obj.remove("error_budget");
    let legacy_res: SolveResult = serde_json::from_value(legacy).unwrap();
    assert!(legacy_res.diagnostics.is_none() && legacy_res.error_budget.is_none());
    assert_eq!(r.scalar("neff"), Some(2.4));
    assert_eq!(r.scalar("polarization"), None);
    assert!(r.solutions[0].field("E").is_some());
}

#[test]
fn descriptors_exemplar_is_stable() {
    let m = Mock::new("tst");
    let v = json!({
        "id": m.id(),
        "studies": serde_json::to_value(m.studies()).unwrap(),
        "bc_kinds": serde_json::to_value(m.bc_kinds()).unwrap(),
        "properties": serde_json::to_value(m.properties()).unwrap(),
        "params": serde_json::to_value(m.params()).unwrap(),
        "fields": serde_json::to_value(m.fields("modes")).unwrap(),
        "link": serde_json::to_value(FieldLink::new("tst.T", "temperature")).unwrap(),
        "sizing": serde_json::to_value(Sizing { order: 2, headroom: 1.5, wavelengths: [("core".to_string(), 5e-7)].into() }).unwrap(),
    });
    check_exemplar("physics_descriptors.v1.json", &v);
    let studies: Vec<StudyDesc> = serde_json::from_value(v["studies"].clone()).unwrap();
    assert_eq!(studies, m.studies());
}

#[test]
fn the_physics_trait_is_object_safe_and_has_these_methods() {
    // a change to this list is a change of the interface: version it (ADR 19)
    let p: Box<dyn Physics> = Box::new(Mock::new("tst"));
    let _ = (
        p.id(),
        p.name(),
        p.studies(),
        p.bc_kinds(),
        p.properties(),
        p.params(),
        p.fields("modes"),
        p.exports(),
        p.imports(),
    );
    let spec = ProblemSpec::new("tst", 2, "modes");
    let _ = p.sizing(&spec, None);
    let sig: fn(&dyn Physics, &Mesh, &ProblemSpec) -> Result<(), Error> = |p, m, s| p.validate(s, m);
    let _ = sig;
    let solve: fn(&dyn Physics, &Mesh, &ProblemSpec, &mut Ctx) -> Result<SolveResult, Error> =
        |p, m, s, c| p.solve(m, s, c);
    let _ = solve;
    type SensFn =
        fn(&dyn Physics, &Mesh, &ProblemSpec, &mut Ctx, &SensitivityRequest) -> Result<SensitivityReport, Error>;
    let sens: SensFn = |p, m, s, c, r| p.sensitivity(m, s, c, r);
    let _ = sens;
    fn send_sync<T: Send + Sync + ?Sized>() {}
    send_sync::<dyn Physics>();
}

#[test]
fn error_codes_are_stable_and_cancel_is_recognised() {
    assert_eq!(Error::invalid("x").code(), "invalid");
    assert_eq!(Error::Unsupported("x".into()).code(), "unsupported");
    assert_eq!(Error::Solver("x".into()).code(), "solver");
    assert_eq!(Error::NonFinite("x".into()).code(), "non_finite");
    assert_eq!(Error::Cancelled.code(), "cancelled");
    assert_eq!(Error::NotYetMoved("x".into()).code(), "not_yet_moved");
    // the kernels of femcore signal a cancel with the text "cancelled"
    assert_eq!(Error::from_kernel("cancelled".into()), Error::Cancelled);
    assert!(matches!(Error::from_kernel("the result contains NaN or infinity".into()), Error::NonFinite(_)));
    assert!(matches!(Error::from_kernel("singular matrix".into()), Error::Solver(_)));
    assert_eq!(Error::Cancelled.to_string(), "cancelled");
}

#[test]
fn ctx_reports_a_cancel() {
    let table = Table;
    let mut ctx = Ctx::new(&table);
    assert!(ctx.check().is_ok());
    ctx.progress.cancel();
    assert_eq!(ctx.check(), Err(Error::Cancelled));
    let mock = Mock::new("tst");
    let mesh = Mesh::new(2, vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0], vec![0, 1, 2], vec![1], vec![], vec![]);
    assert_eq!(mock.solve(&mesh, &ProblemSpec::new("tst", 1, "modes"), &mut ctx), Err(Error::Cancelled));
}

#[test]
fn the_registry_takes_three_lowercase_letters_once() {
    let mut r = Registry::new();
    assert!(r.is_empty());
    r.register(Arc::new(Mock::new("tst"))).unwrap();
    assert!(r.register(Arc::new(Mock::new("tst"))).is_err());
    for bad in ["ab", "abcd", "AB1", "Abc", "a-b"] {
        assert!(r.register(Arc::new(Mock::new(Box::leak(bad.to_string().into_boxed_str())))).is_err(), "{bad}");
    }
    assert_eq!(r.ids(), vec!["tst"]);
    assert_eq!(r.len(), 1);
    assert_eq!(r.get("tst").unwrap().name(), "Mock physics");
    assert!(r.get("xyz").is_none());
}

#[test]
fn a_one_way_chain_is_solved_in_dependency_order_and_a_cycle_is_refused() {
    let mut heat = Mock::new("aaa");
    heat.exports = vec![FieldLink::new("aaa.T", "temperature")];
    let mut solid = Mock::new("zzz");
    solid.imports = vec![FieldLink::new("zzz.temperature", "temperature")];
    let optics = Mock::new("mmm");
    let mut r = Registry::new();
    for p in [heat, solid, optics] {
        r.register(Arc::new(p)).unwrap();
    }
    // the heat physics has the alphabetically first code but the dependency decides the order anyway: zzz after aaa
    assert_eq!(r.solve_order(&["zzz", "aaa"]).unwrap(), vec!["aaa", "zzz"]);
    // an import with no exporter among the chosen physics is not an error
    assert_eq!(r.solve_order(&["zzz"]).unwrap(), vec!["zzz"]);
    assert_eq!(r.solve_order(&["mmm", "aaa", "zzz"]).unwrap().len(), 3);
    assert!(r.solve_order(&["nope"]).is_err());
    // a cycle: two physics that import what the other exports
    let mut a = Mock::new("ccc");
    a.exports = vec![FieldLink::new("ccc.A", "alpha")];
    a.imports = vec![FieldLink::new("ccc.B", "beta")];
    let mut b = Mock::new("ddd");
    b.exports = vec![FieldLink::new("ddd.B", "beta")];
    b.imports = vec![FieldLink::new("ddd.A", "alpha")];
    let mut r2 = Registry::new();
    r2.register(Arc::new(a)).unwrap();
    r2.register(Arc::new(b)).unwrap();
    let e = r2.solve_order(&["ccc", "ddd"]).unwrap_err();
    assert!(e.to_string().contains("cycle"), "{e}");
}

#[test]
fn canonical_json_sorts_keys_at_every_depth() {
    let v = json!({ "b": { "z": 1, "a": [ { "y": 1, "x": 2 } ] }, "a": 0 });
    let c = canonical_json(&v);
    let keys: Vec<&String> = c.as_object().unwrap().keys().collect();
    assert_eq!(keys, vec!["a", "b"]);
    let inner: Vec<&String> = c["b"].as_object().unwrap().keys().collect();
    assert_eq!(inner, vec!["a", "z"]);
    let deep: Vec<&String> = c["b"]["a"][0].as_object().unwrap().keys().collect();
    assert_eq!(deep, vec!["x", "y"]);
    assert_eq!(c, v, "equal as values whatever the order");
    let _: Value = c;
}
