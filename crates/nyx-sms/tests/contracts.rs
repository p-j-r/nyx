//! Contract of `nyx-sms` (ADR 19, ADR 25): its descriptors are pinned by `docs/architecture/contracts/physics_sms.v1.json`, and an example problem
//! passes the structural validation. The kernels are not in this crate yet (restructure phase 4), so `solve` answers `NotYetMoved`.
//!
//! Rewrite the exemplar after a deliberate change (with a `VERSIONS.md` line):
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-sms --test contracts`
use nyx_physics::*;
use nyx_sms::Sms;
use serde_json::json;

fn example() -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, "eigenfrequency");
    s.study.params = json!({ "n_modes": 6, "f_search": 1e9, "kz": 3.0 });
    s.materials.insert(
        "beam".into(),
        MaterialAssign::Props([("E".to_string(), 169e9), ("nu".to_string(), 0.22), ("rho".to_string(), 2329.0)].into()),
    );
    s.bcs.push(BcSpec {
        on: "outer & plane(x = 0)".into(),
        kind: "fixed".into(),
        params: json!({ "ux": 0.0, "uz": 0.0 }),
    });
    s.bcs.push(BcSpec {
        on: "substrate_pml".into(),
        kind: "pml".into(),
        params: json!({ "z": "minus", "strength": 3.0 }),
    });
    s.params = json!({ "plane": "strain", "inactive": ["air"] });
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
    let v = contract::descriptor_snapshot(&Sms);
    contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), "physics_sms.v1.json", &v);
}

#[test]
fn the_code_is_three_lowercase_letters_and_the_registry_takes_it() {
    assert_eq!(Sms.id(), "sms");
    assert!(registry::is_valid_code(Sms.id()));
    let mut r = Registry::new();
    r.register(std::sync::Arc::new(Sms)).unwrap();
}

#[test]
fn the_example_problem_passes_validation_and_every_study_has_fields() {
    validate_spec(&Sms, &example()).unwrap();
    for s in Sms.studies() {
        assert!(!Sms.fields(&s.kind).is_empty(), "study {} lists no fields", s.kind);
    }
    // a property the physics does not read is refused (R12)
    let mut bad = example();
    bad.materials.insert("zzz".into(), MaterialAssign::Props([("not_a_property".to_string(), 1.0)].into()));
    assert!(validate_spec(&Sms, &bad).is_err());
    // an unknown study kind is refused with the list of the known ones
    let mut bad = example();
    bad.study.kind = "nope".into();
    assert!(validate_spec(&Sms, &bad).unwrap_err().to_string().contains(&Sms.studies()[0].kind));
}

#[test]
fn the_json_schema_accepts_the_example_shape() {
    let schema = schema::physics_schema(&Sms);
    assert_eq!(schema["properties"]["physics"]["const"], json!("sms"));
    let kinds: Vec<_> = schema["properties"]["study"]["oneOf"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["properties"]["kind"]["const"].clone())
        .collect();
    assert_eq!(kinds.len(), Sms.studies().len());
}

#[test]
fn solve_is_not_moved_in_yet() {
    let mesh = Mesh::new(2, vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0], vec![0, 1, 2], vec![1], vec![], vec![]);
    let table = NoNames;
    let mut ctx = Ctx::new(&table);
    let e = Sms.solve(&mesh, &example(), &mut ctx).unwrap_err();
    assert_eq!(e.code(), "not_yet_moved");
    // and the sizing hint carries the order
    assert_eq!(Sms.sizing(&example()).order, example().order);
}
