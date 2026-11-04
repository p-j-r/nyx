//! Contract of `nyx-dfs` (ADR 19, ADR 25): its descriptors are pinned by `docs/architecture/contracts/physics_dfs.v1.json`, and an example problem
//! passes the structural validation; the envelope solves through `Physics::solve` and agrees with the typed problem; the stage functions are public.
//!
//! Rewrite the exemplar after a deliberate change (with a `VERSIONS.md` line):
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-dfs --test contracts`
use nyx_dfs::Dfs;
use nyx_physics::*;
use serde_json::json;

fn example() -> ProblemSpec {
    let mut s = ProblemSpec::new("dfs", 1, "heat_transient");
    s.study.params = json!({ "t_end": 2.0, "n_out": 10, "t0": "293.15" });
    s.materials.insert(
        "plate".into(),
        MaterialAssign::Props(
            [("k".to_string(), 150.0), ("rho".to_string(), 2329.0), ("Cp".to_string(), 700.0)].into(),
        ),
    );
    s.bcs.push(BcSpec { on: "outer".into(), kind: "convection".into(), params: json!({ "h": 10.0, "t_ext": 293.15 }) });
    s
}

#[test]
fn descriptors_exemplar_is_stable() {
    let v = contract::descriptor_snapshot(&Dfs);
    contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), "physics_dfs.v1.json", &v);
}

#[test]
fn the_code_is_three_lowercase_letters_and_the_registry_takes_it() {
    assert_eq!(Dfs.id(), "dfs");
    assert!(registry::is_valid_code(Dfs.id()));
    let mut r = Registry::new();
    r.register(std::sync::Arc::new(Dfs)).unwrap();
}

#[test]
fn the_example_problem_passes_validation_and_every_study_has_fields() {
    validate_spec(&Dfs, &example()).unwrap();
    for s in Dfs.studies() {
        assert!(!Dfs.fields(&s.kind).is_empty(), "study {} lists no fields", s.kind);
    }
    // a property the physics does not read is refused (R12)
    let mut bad = example();
    bad.materials.insert("zzz".into(), MaterialAssign::Props([("not_a_property".to_string(), 1.0)].into()));
    assert!(validate_spec(&Dfs, &bad).is_err());
    // an unknown study kind is refused with the list of the known ones
    let mut bad = example();
    bad.study.kind = "nope".into();
    assert!(validate_spec(&Dfs, &bad).unwrap_err().to_string().contains(&Dfs.studies()[0].kind));
}

#[test]
fn the_json_schema_accepts_the_example_shape() {
    let schema = schema::physics_schema(&Dfs);
    assert_eq!(schema["properties"]["physics"]["const"], json!("dfs"));
    let kinds: Vec<_> = schema["properties"]["study"]["oneOf"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["properties"]["kind"]["const"].clone())
        .collect();
    assert_eq!(kinds.len(), Dfs.studies().len());
}

/// A table of names for the tests: region `plate` is domain 1, the boundary selector `outer` is boundaries 1 to 4 and `left`/`right` are 1 and 4
/// (the numbering of `femcore::testkit::rect`).
struct Table;
impl LabelResolver for Table {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        match r {
            "plate" => Ok(vec![1]),
            _ => Err(Error::invalid(format!("no region {r}"))),
        }
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        match s {
            "outer" => Ok(vec![1, 2, 3, 4]),
            "left" => Ok(vec![1]),
            "right" => Ok(vec![4]),
            _ => Err(Error::invalid(format!("no selector {s}"))),
        }
    }
}

#[test]
fn the_envelope_solves_through_the_trait_and_agrees_with_the_typed_problem() {
    // electrostatics of a parallel-plate capacitor, once through `Physics::solve` (names resolved by the table) and once through the typed
    // problem with the numbers written out: the two must give the same bits
    let mesh = femcore::testkit::rect(1.0, 0.5, 6, 3);
    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.materials.insert("plate".into(), MaterialAssign::Props([("eps_r".to_string(), 4.2)].into()));
    spec.bcs.push(BcSpec { on: "left".into(), kind: "potential".into(), params: json!({ "value": 1.0 }) });
    spec.bcs.push(BcSpec { on: "right".into(), kind: "ground".into(), params: serde_json::Value::Null });
    let table = Table;
    let mut ctx = Ctx::new(&table);
    let r = Dfs.solve(&mesh, &spec, &mut ctx).unwrap();
    assert_eq!((r.physics.as_str(), r.study.as_str(), r.order, r.dim), ("dfs", "electrostatic", 2, 2));
    let mut p = nyx_dfs::Problem::new(nyx_dfs::Study::Electrostatic, 2);
    p.materials.set(1, "eps_r", 4.2);
    p.bcs = vec![
        nyx_dfs::Bc { labels: vec![1], kind: nyx_dfs::BcKind::Potential(1.0) },
        nyx_dfs::Bc { labels: vec![4], kind: nyx_dfs::BcKind::Ground },
    ];
    let t = nyx_dfs::solve_electrostatics(&mesh, &p).unwrap();
    assert_eq!(r.solutions[0].fields[0].values, t.solutions[0].fields[0].values);
    // C = eps0 eps_r (0.5 / 1) per metre
    let c = r.scalar("capacitance").unwrap();
    assert!((c / (nyx_dfs::EPS0 * 4.2) - 0.5).abs() < 1e-6);
    // the result is plain data: it survives a round trip through JSON
    let back: SolveResult = serde_json::from_str(&serde_json::to_string(&r).unwrap()).unwrap();
    assert_eq!(back.solutions[0].fields[0].values, r.solutions[0].fields[0].values);
    // and the sizing hint carries the order
    assert_eq!(Dfs.sizing(&spec).order, spec.order);
}

#[test]
fn the_envelope_is_refused_where_it_is_wrong() {
    let mesh = femcore::testkit::rect(1.0, 0.5, 4, 2);
    let table = Table;
    // an unknown region
    let mut spec = ProblemSpec::new("dfs", 2, "heat_stationary");
    spec.materials.insert("nowhere".into(), MaterialAssign::Props([("k".to_string(), 1.0)].into()));
    assert!(Dfs.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().to_string().contains("no region nowhere"));
    // a heat condition in an electrostatics study
    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.bcs.push(BcSpec { on: "left".into(), kind: "temperature".into(), params: json!({ "value": 300.0 }) });
    let e = Dfs.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err();
    assert_eq!(e.code(), "invalid");
    assert!(e.to_string().contains("not an electrostatics boundary condition"), "{e}");
    // a named material needs the material files of phase 6
    let mut spec = ProblemSpec::new("dfs", 2, "heat_stationary");
    spec.materials.insert("plate".into(), MaterialAssign::Named("silicon".into()));
    assert_eq!(Dfs.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().code(), "unsupported");
    // a cancelled run says so
    let mut spec = ProblemSpec::new("dfs", 2, "electrostatic");
    spec.bcs.push(BcSpec { on: "left".into(), kind: "ground".into(), params: serde_json::Value::Null });
    let mut ctx = Ctx::new(&table);
    ctx.progress.cancel();
    assert_eq!(Dfs.solve(&mesh, &spec, &mut ctx).unwrap_err().code(), "cancelled");
}

#[test]
fn kernel_stage_api_is_public() {
    // ADR 9: naming the functions makes removing or renaming one a compile error here
    let _ = (nyx_dfs::assemble_electrostatics, nyx_dfs::solve_potential, nyx_dfs::electrostatics_fields);
    let _ = (nyx_dfs::assemble_heat, nyx_dfs::solve_heat_states, nyx_dfs::heat_fields);
    let _ = (nyx_dfs::solve_electrostatics, nyx_dfs::solve_heat, nyx_dfs::run, nyx_dfs::decode);
}
