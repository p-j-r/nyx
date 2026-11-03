//! Contract of `nyx-dfs` (ADR 19, ADR 25): its descriptors are pinned by `docs/architecture/contracts/physics_dfs.v1.json`, and an example problem
//! passes the structural validation. The kernels are not in this crate yet (restructure phase 4), so `solve` answers `NotYetMoved`.
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

struct NoNames;
impl LabelResolver for NoNames {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        Err(Error::invalid(format!("no region {r}")))
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        Err(Error::invalid(format!("no selector {s}")))
    }
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

#[test]
fn solve_is_not_moved_in_yet() {
    let mesh = Mesh::new(2, vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0], vec![0, 1, 2], vec![1], vec![], vec![]);
    let table = NoNames;
    let mut ctx = Ctx::new(&table);
    let e = Dfs.solve(&mesh, &example(), &mut ctx).unwrap_err();
    assert_eq!(e.code(), "not_yet_moved");
    // and the sizing hint carries the order
    assert_eq!(Dfs.sizing(&example()).order, example().order);
}
