//! Behaviour of the engine core with stand-in physics (the real ones answer `NotYetMoved` until restructure phase 4): the project is built
//! through commands only (the strip, the rib and the ring of gate G5), meshed with the real mesher, solved through the `Physics` trait in
//! coupling order, cancelled, stored and reopened.
use nyx_engine::*;
use nyx_physics::*;
use serde_json::{json, Value};
use std::sync::Arc;

// ------------------------------------------------------------------------------------------------------ stand-in physics

/// `tst`: writes the label of the region of each material into a field `T`, resolves the selector of each condition, exports `T` as
/// `temperature`.
struct Tst;
impl Physics for Tst {
    fn id(&self) -> &'static str {
        "tst"
    }
    fn name(&self) -> &'static str {
        "Test physics"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc {
            kind: "const".into(),
            name: "Constant".into(),
            doc: "".into(),
            dims: Dims::Both,
            params: vec![],
        }]
    }
    fn bc_kinds(&self) -> Vec<BcDesc> {
        vec![BcDesc {
            kind: "fix".into(),
            name: "Fixed".into(),
            doc: "".into(),
            applies: Applies::Boundary,
            params: vec![ParamDesc::number("value", "Value", "", 0.0, "")],
        }]
    }
    fn properties(&self) -> Vec<PropDesc> {
        vec![PropDesc::new("k", "Conductivity", "W/m/K", PropKind::Scalar, &[], Some(1.0))]
    }
    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![FieldDesc::new("T", "K", &[], "")]
    }
    fn exports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("tst.T", "temperature")]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn validate(&self, _s: &ProblemSpec, _m: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, _s: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        Sizing::default()
    }
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        let mut t = vec![0.0; 3];
        for region in spec.materials.keys() {
            let labels = ctx.resolver.domains(region)?;
            t.push(labels[0] as f64);
        }
        let mut n_boundary = 0;
        for bc in &spec.bcs {
            n_boundary += ctx.resolver.boundaries(&bc.on)?.len();
        }
        Ok(SolveResult {
            physics: "tst".into(),
            study: spec.study.kind.clone(),
            order: spec.order,
            dim: mesh.dim() as u8,
            dofs: mesh.n_cells(),
            seconds: 0.01,
            solutions: vec![Solution {
                label: "only".into(),
                fields: vec![Field { name: "T".into(), unit: "K".into(), components: vec![], values: t }],
                value: Some(1.0),
                ..Default::default()
            }],
            scalars: vec![
                ("boundary labels".into(), Scalar::Num(n_boundary as f64)),
                ("cells".into(), Scalar::Num(mesh.n_cells() as f64)),
            ],
            log: vec!["solved the stand-in".into()],
        })
    }
}

/// `use`: imports `temperature` and reports the sum of what it got; proves the order of a chain and the value passing.
struct Use;
impl Physics for Use {
    fn id(&self) -> &'static str {
        "use"
    }
    fn name(&self) -> &'static str {
        "Importer"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc { kind: "sum".into(), name: "Sum".into(), doc: "".into(), dims: Dims::Both, params: vec![] }]
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
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![FieldDesc::new("S", "", &[], "")]
    }
    fn exports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![FieldLink::new("use.T", "temperature")]
    }
    fn validate(&self, _s: &ProblemSpec, _m: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, _s: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        Sizing::default()
    }
    fn solve(&self, _mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        let t = ctx.inputs.get("temperature").ok_or_else(|| Error::invalid("no temperature was imported"))?;
        Ok(SolveResult {
            physics: "use".into(),
            study: spec.study.kind.clone(),
            solutions: vec![Solution { label: "s".into(), ..Default::default() }],
            scalars: vec![("sum".into(), Scalar::Num(t.iter().sum()))],
            ..Default::default()
        })
    }
}

/// `slo`: runs until cancelled (at most 20 s).
struct Slo;
impl Physics for Slo {
    fn id(&self) -> &'static str {
        "slo"
    }
    fn name(&self) -> &'static str {
        "Slow"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc { kind: "wait".into(), name: "Wait".into(), doc: "".into(), dims: Dims::Both, params: vec![] }]
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
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![]
    }
    fn exports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn validate(&self, _s: &ProblemSpec, _m: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, _s: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        Sizing::default()
    }
    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        for _ in 0..2000 {
            ctx.check()?;
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        Ok(SolveResult::default())
    }
}

fn test_engine() -> Engine {
    let mut r = Registry::new();
    r.register(Arc::new(Tst)).unwrap();
    r.register(Arc::new(Use)).unwrap();
    r.register(Arc::new(Slo)).unwrap();
    Engine::with(r, Arc::new(FemcoreMesher))
}

fn temp(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("nyx-engine-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn call(e: &mut Engine, name: &str, args: Value) -> Value {
    e.call(name, &args).unwrap_or_else(|err| panic!("{name} failed: {err}"))
}

/// The strip of the SOI wire: oxide, silicon core, air above. Nothing but commands.
fn build_strip(e: &mut Engine) {
    call(e, "project_new", json!({ "title": "Strip", "h": 2e-7 }));
    call(e, "set_param", json!({ "name": "w", "value": 5e-7 }));
    call(e, "set_param", json!({ "name": "t", "value": 2.2e-7 }));
    call(e, "add_rect", json!({ "name": "window", "min": [-2e-6, -1.5e-6], "max": [2e-6, 1.5e-6] }));
    call(e, "add_rect", json!({ "name": "box", "min": [-2e-6, -1.5e-6], "max": [2e-6, 0.0] }));
    call(e, "add_rect", json!({ "name": "core", "min": ["-w/2", 0.0], "max": ["w/2", "t"] }));
    call(
        e,
        "add_overlay",
        json!({ "regions": [
        { "name": "air", "input": "window" }, { "name": "oxide", "input": "box" }, { "name": "silicon", "input": "core" } ] }),
    );
}

// ------------------------------------------------------------------------------------------------------ tests

#[test]
fn the_strip_is_built_meshed_and_solved_through_commands_alone() {
    let mut e = test_engine();
    build_strip(&mut e);
    let names = call(&mut e, "geometry_names", json!({}));
    let regions: Vec<&str> = names["regions"].as_array().unwrap().iter().map(|r| r["name"].as_str().unwrap()).collect();
    assert_eq!(regions, ["air", "oxide", "silicon"], "regions are numbered in overlay order");
    assert!(names["boundaries"].as_array().unwrap().len() >= 3);

    call(&mut e, "add_physics", json!({ "physics": "tst", "study": "const" }));
    call(&mut e, "set_material", json!({ "physics": "tst", "region": "silicon", "material": { "k": 150.0 } }));
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "tst", "on": "outer", "kind": "fix", "params": { "value": 293.0 } }),
    );
    let mesh = call(&mut e, "mesh_build", json!({}));
    assert!(mesh["mesh"]["cells"].as_u64().unwrap() > 50, "{mesh}");
    let solved = call(&mut e, "solve_run", json!({}));
    assert_eq!(solved["results"][0]["physics"], "tst");
    // the stand-in wrote the label of the region of the one material: silicon is region 3
    let field = call(&mut e, "result_field", json!({ "physics": "tst", "field": "T" }));
    assert_eq!(field["values"], json!([0.0, 0.0, 0.0, 3.0]));
    let info = call(&mut e, "result_info", json!({}));
    assert_eq!(info[0]["stale"], false);
    assert!(info[0]["solutions"][0]["fields"][0]["name"] == "T");
    // the condition `outer` resolved to boundary labels
    assert!(solved["results"][0]["scalars"].to_string().contains("boundary labels"));

    // an edit makes the result stale and the mesh not current, without any invalidation call
    call(&mut e, "set_param", json!({ "name": "w", "value": 6e-7 }));
    assert_eq!(call(&mut e, "result_info", json!({}))[0]["stale"], true);
    assert!(call(&mut e, "project_info", json!({}))["mesh"].is_null());
}

#[test]
fn the_rib_and_the_ring_are_built_and_meshed_through_commands_alone() {
    let mut e = test_engine();
    call(&mut e, "project_new", json!({ "title": "Rib", "h": 3e-7 }));
    call(&mut e, "add_rect", json!({ "name": "window", "min": [-2e-6, -1.5e-6], "max": [2e-6, 1.5e-6] }));
    call(&mut e, "add_rect", json!({ "name": "slab", "min": [-2e-6, 0.0], "max": [2e-6, 1e-7] }));
    call(&mut e, "add_rect", json!({ "name": "rib", "min": [-2.5e-7, 1e-7], "max": [2.5e-7, 3e-7] }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [
        { "name": "air", "input": "window" }, { "name": "slab", "input": "slab" }, { "name": "rib", "input": "rib" } ] }),
    );
    assert_eq!(call(&mut e, "mesh_build", json!({}))["mesh"]["regions"].as_array().unwrap().len(), 3);

    call(&mut e, "project_new", json!({ "title": "Ring", "h": 4e-7 }));
    call(&mut e, "add_rect", json!({ "name": "window", "min": [-8e-6, -8e-6], "max": [8e-6, 8e-6] }));
    call(&mut e, "add_circle", json!({ "name": "outer_disc", "center": [0.0, 0.0], "r": 5e-6 }));
    call(&mut e, "add_circle", json!({ "name": "inner_disc", "center": [0.0, 0.0], "r": 4.5e-6 }));
    call(&mut e, "add_difference", json!({ "name": "ring", "a": "outer_disc", "b": ["inner_disc"] }));
    call(
        &mut e,
        "add_overlay",
        json!({ "regions": [ { "name": "air", "input": "window" }, { "name": "ring", "input": "ring" } ] }),
    );
    let m = call(&mut e, "mesh_build", json!({}));
    assert_eq!(m["mesh"]["regions"].as_array().unwrap().len(), 2, "{m}");
    // editing by name: widen the ring by moving the inner disc
    call(&mut e, "update_feature", json!({ "name": "inner_disc", "fields": { "r": 4.0e-6 } }));
    assert_eq!(call(&mut e, "list_features", json!({}))["features"][2]["args"]["r"], 4.0e-6);
}

#[test]
fn refusals_leave_the_project_as_it_was_and_say_what_is_wrong() {
    let mut e = test_engine();
    build_strip(&mut e);
    let before = e.project().clone();
    let err = |e: &mut Engine, name: &str, a: Value| e.call(name, &a).unwrap_err().to_string();
    assert!(err(&mut e, "add_rect", json!({ "name": "box", "min": [0, 0], "max": [1, 1] })).contains("exists"));
    assert!(err(&mut e, "add_rect", json!({ "min": [0, 0] })).contains("needs `max`"));
    assert!(err(&mut e, "remove_feature", json!({ "name": "box" })).contains("is used by overlay1"));
    assert!(err(&mut e, "add_overlay", json!({ "regions": [{ "name": "x", "input": "core" }] }))
        .contains("only one Overlay"));
    assert!(err(&mut e, "add_union", json!({ "ids": ["nothing"] })).contains("no feature named 'nothing'"));
    assert!(err(&mut e, "rename_feature", json!({ "name": "core", "new_name": "a.b" })).contains("reserved"));
    assert!(err(&mut e, "frobnicate", json!({})).contains("list_commands"));
    assert!(err(&mut e, "add_physics", json!({ "physics": "zzz" })).contains("no physics 'zzz'"));
    assert!(err(&mut e, "set_material", json!({ "physics": "tst", "region": "silicon", "material": "x" }))
        .contains("the project has no physics 'tst'"));
    call(&mut e, "add_physics", json!({ "physics": "tst" }));
    assert!(err(&mut e, "set_material", json!({ "physics": "tst", "region": "nowhere", "material": "x" }))
        .contains("no region 'nowhere'"));
    assert_eq!(e.project().geometry, before.geometry, "no refused edit changed the geometry");
    // the project check names a property no physics reads
    call(&mut e, "set_material", json!({ "physics": "tst", "region": "silicon", "material": { "bogus": 1.0 } }));
    assert!(err(&mut e, "check_project", json!({})).contains("does not read a property 'bogus'"));
}

#[test]
fn a_coupled_chain_solves_in_order_and_passes_the_exported_field() {
    let mut e = test_engine();
    build_strip(&mut e);
    call(&mut e, "add_physics", json!({ "physics": "use", "study": "sum" }));
    call(&mut e, "add_physics", json!({ "physics": "tst" }));
    call(&mut e, "set_material", json!({ "physics": "tst", "region": "oxide", "material": { "k": 1.4 } }));
    let r = call(&mut e, "solve_run", json!({ "physics": ["use", "tst"] }));
    let order: Vec<&str> = r["results"].as_array().unwrap().iter().map(|x| x["physics"].as_str().unwrap()).collect();
    assert_eq!(order, ["tst", "use"], "the exporter is solved first");
    // tst wrote [0, 0, 0, 2] (oxide is region 2); use summed it
    let sum = call(&mut e, "result_info", json!({ "physics": "use" }))[0]["scalars"].clone();
    assert_eq!(sum[0][1], 2.0, "{sum}");
}

#[test]
fn a_running_solve_can_be_cancelled_and_the_cancelled_run_is_recorded() {
    let mut e = test_engine();
    build_strip(&mut e);
    call(&mut e, "add_physics", json!({ "physics": "slo", "study": "wait" }));
    let ws = temp("cancel");
    call(&mut e, "store_save", json!({ "workspace": ws.to_str().unwrap() }));
    let id = call(&mut e, "solve_run", json!({ "wait": false }))["job"].as_u64().unwrap();
    assert_eq!(call(&mut e, "job_status", json!({ "id": id }))["state"], "running");
    assert_eq!(call(&mut e, "job_cancel", json!({ "id": id }))["cancelled"], true);
    let err = e.call("job_wait", &json!({ "id": id })).unwrap_err();
    assert_eq!(err.code(), "cancelled");
    assert_eq!(call(&mut e, "job_status", json!({ "id": id }))["state"], "cancelled");
    let runs = call(&mut e, "store_runs", json!({}));
    assert_eq!(runs[0]["status"], "cancelled", "{runs}");
    // a second wait finds nothing to apply, and says so
    assert!(e.call("job_wait", &json!({ "id": id })).unwrap_err().to_string().contains("applied already"));
    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn polling_applies_a_finished_job_without_blocking() {
    let mut e = test_engine();
    build_strip(&mut e);
    let id = call(&mut e, "mesh_build", json!({ "wait": false }))["job"].as_u64().unwrap();
    let mut applied = vec![];
    for _ in 0..600 {
        applied = e.poll();
        if !applied.is_empty() {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(20));
    }
    assert_eq!(applied.len(), 1);
    assert_eq!(applied[0].0, id);
    assert!(applied[0].1.as_ref().unwrap()["mesh"]["cells"].as_u64().unwrap() > 0);
    assert!(e.model().is_some());
}

#[test]
fn a_project_is_stored_solved_into_the_run_table_and_reopened() {
    let mut e = test_engine();
    build_strip(&mut e);
    call(&mut e, "add_physics", json!({ "physics": "tst", "study": "const" }));
    call(&mut e, "set_material", json!({ "physics": "tst", "region": "silicon", "material": { "k": 150.0 } }));
    let ws = temp("store");
    let saved = call(&mut e, "store_save", json!({ "workspace": ws.to_str().unwrap() }));
    let folder = std::path::PathBuf::from(saved["folder"].as_str().unwrap());
    assert_eq!(saved["version"], 1);

    call(&mut e, "solve_run", json!({}));
    let runs = call(&mut e, "store_runs", json!({}));
    assert_eq!(runs.as_array().unwrap().len(), 1);
    assert_eq!(runs[0]["status"], "ok");
    assert!(runs[0]["metrics"]["cells"].as_f64().unwrap() > 50.0);
    assert_eq!(runs[0]["params"]["w"], 5e-7, "numeric model parameters are the run's parameters");
    let run_dir = folder.join("e1").join("runs").join(runs[0]["id"].as_str().unwrap());
    assert!(run_dir.join("result_tst.json").exists() && run_dir.join("log.txt").exists());

    // an edit and a save is a new version; the card is made by the engine
    call(&mut e, "set_param", json!({ "name": "w", "value": 6e-7 }));
    let again = call(&mut e, "store_save", json!({}));
    assert_eq!((again["version"].clone(), again["changed"].clone()), (json!(2), json!(true)));
    let card = nyx_store::WorkflowDir::open_with(&folder, nyx_engine::store_link::card_of).unwrap().wf.experiments[0]
        .card
        .clone();
    assert!(card.summary.starts_with("2D, 4 features, 3 regions, tst const"), "{}", card.summary);
    assert!(card.keywords.contains(&"silicon".to_string()));

    // reopen in a new engine: the same project
    let mut f = test_engine();
    call(&mut f, "store_open", json!({ "folder": folder.to_str().unwrap() }));
    assert_eq!(f.project(), e.project());
    assert_eq!(call(&mut f, "store_runs", json!({})).as_array().unwrap().len(), 1);
    let _ = std::fs::remove_dir_all(&ws);
}

#[test]
fn the_project_json_is_stable_and_round_trips() {
    let mut e = test_engine();
    build_strip(&mut e);
    call(&mut e, "add_physics", json!({ "physics": "tst" }));
    let a = call(&mut e, "project_json", json!({}))["json"].as_str().unwrap().to_string();
    assert_eq!(a, e.project().to_json());
    let mut f = test_engine();
    call(&mut f, "project_load", json!({ "json": a.clone() }));
    assert_eq!(f.project(), e.project());
    assert_eq!(f.project().to_json(), a);
    let path = temp("file").join("p.nyx.json");
    call(&mut e, "project_save", json!({ "path": path.to_str().unwrap() }));
    call(&mut f, "project_new", json!({}));
    call(&mut f, "project_open", json!({ "path": path.to_str().unwrap() }));
    assert_eq!(f.project(), e.project());
    assert!(f
        .call("project_load", &json!({ "json": "{\"title\": 1}" }))
        .unwrap_err()
        .to_string()
        .contains("not a nyx project"));
}

#[test]
fn dfs_solves_the_strip_through_the_engine() {
    // restructure phase 4: the real `dfs` of this build solves through `Physics::solve`, with regions and the boundary selector resolved by
    // the engine's label map. A heat source in the silicon core, a fixed temperature on the outer boundary.
    let mut e = Engine::new();
    build_strip(&mut e);
    call(&mut e, "add_physics", json!({ "physics": "dfs", "study": "heat_stationary" }));
    for (region, k, q) in [("air", 0.026, 0.0), ("oxide", 1.4, 0.0), ("silicon", 150.0, 1e13)] {
        call(&mut e, "set_material", json!({ "physics": "dfs", "region": region, "material": { "k": k, "Q": q } }));
    }
    call(
        &mut e,
        "add_condition",
        json!({ "physics": "dfs", "on": "outer", "kind": "temperature", "params": { "value": 300.0 } }),
    );
    call(&mut e, "mesh_build", json!({}));
    let solved = call(&mut e, "solve_run", json!({}));
    assert_eq!(solved["results"][0]["physics"], "dfs");
    let field = call(&mut e, "result_field", json!({ "physics": "dfs", "field": "T" }));
    let t: Vec<f64> = field["values"].as_array().unwrap().iter().map(|v| v.as_f64().unwrap()).collect();
    let (lo, hi) = t.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    assert!((lo - 300.0).abs() < 1e-6, "the boundary holds 300 K, the minimum is {lo}");
    assert!(hi > 300.0 + 1e-3, "the heated core is warmer than the boundary: {hi}");
    // a condition that belongs to another study of the same physics is refused with the reason
    call(&mut e, "add_condition", json!({ "physics": "dfs", "on": "outer", "kind": "ground" }));
    let err = e.call("solve_run", &json!({})).unwrap_err();
    assert_eq!(err.code(), "invalid", "{err}");
    assert!(err.to_string().contains("not a heat-transfer boundary condition"), "{err}");
    assert_eq!(e.registry().ids(), ["dfs", "ems", "sms"]);
    // the sibling physics are still skeletons until their phase 4 threads land; every physics is described from its descriptors
    for id in e.registry().ids() {
        let d = call(&mut e, "describe_physics", json!({ "physics": id }));
        assert!(!d["studies"].as_array().unwrap().is_empty());
    }
}

#[test]
fn every_geometry_operation_has_an_add_command_and_every_command_has_a_schema() {
    let e = Engine::new();
    let cat = e.catalog();
    let names: Vec<&str> = cat.iter().map(|c| c.name.as_str()).collect();
    let mut sorted = names.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), names.len(), "a command name is used twice");
    for d in nyx_engine::geomdesc::all() {
        assert!(names.contains(&nyx_engine::geomdesc::tool_name(d.kind).as_str()), "no command for {}", d.kind);
    }
    for c in &cat {
        assert_eq!(c.schema["type"], "object", "{}", c.name);
        assert!(!c.title.is_empty(), "{}", c.name);
    }
    // no command is in the catalog but unhandled: with no arguments each one succeeds or says what it needs
    for c in &cat {
        let mut e = Engine::new();
        if let Err(err) = e.call(&c.name, &json!({})) {
            assert!(!err.to_string().contains("there is no command"), "{} is in the catalog but not handled", c.name);
        }
    }
}
