//! Contract of `nyx-ems` (ADR 19, ADR 25, ADR 45): its descriptors are pinned by `docs/architecture/contracts/physics_ems.v1.json`, and an
//! example problem passes the structural validation; the envelope solves through `Physics::solve` and agrees with the typed problem; the stage
//! functions are public; the shape of a result is pinned by `physics_ems_result.v1.json`.
//!
//! Rewrite the exemplars after a deliberate change (with a `VERSIONS.md` line):
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-ems --test contracts`
mod common;
use nyx_ems::Ems;
use nyx_physics::*;
use serde_json::json;

fn example() -> ProblemSpec {
    let mut s = ProblemSpec::new("ems", 2, "modes");
    s.study.params = json!({ "freq": 1.934e14, "n_target": 2.4, "evanescent": false, "n_modes": 4 });
    s.materials.insert("core".into(), MaterialAssign::Props([("eps_r".to_string(), 12.08)].into()));
    s.materials.insert("clad".into(), MaterialAssign::Props([("eps_r".to_string(), 2.085)].into()));
    s.bcs.push(BcSpec { on: "clad_pml".into(), kind: "pml".into(), params: json!({ "x": "plus" }) });
    s
}

#[test]
fn descriptors_exemplar_is_stable() {
    let v = contract::descriptor_snapshot(&Ems);
    contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), "physics_ems.v1.json", &v);
}

#[test]
fn the_code_is_three_lowercase_letters_and_the_registry_takes_it() {
    assert_eq!(Ems.id(), "ems");
    assert!(registry::is_valid_code(Ems.id()));
    let mut r = Registry::new();
    r.register(std::sync::Arc::new(Ems)).unwrap();
}

#[test]
fn the_example_problem_passes_validation_and_every_study_has_fields() {
    validate_spec(&Ems, &example()).unwrap();
    for s in Ems.studies() {
        assert!(!Ems.fields(&s.kind).is_empty(), "study {} lists no fields", s.kind);
    }
    // a property the physics does not read is refused (R12)
    let mut bad = example();
    bad.materials.insert("zzz".into(), MaterialAssign::Props([("not_a_property".to_string(), 1.0)].into()));
    assert!(validate_spec(&Ems, &bad).is_err());
    // an unknown study kind is refused with the list of the known ones
    let mut bad = example();
    bad.study.kind = "nope".into();
    assert!(validate_spec(&Ems, &bad).unwrap_err().to_string().contains(&Ems.studies()[0].kind));
}

#[test]
fn the_json_schema_accepts_the_example_shape() {
    let schema = schema::physics_schema(&Ems);
    assert_eq!(schema["properties"]["physics"]["const"], json!("ems"));
    let kinds: Vec<_> = schema["properties"]["study"]["oneOf"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| s["properties"]["kind"]["const"].clone())
        .collect();
    assert_eq!(kinds.len(), Ems.studies().len());
}

/// A table of names for the tests. Regions: `guide` and `air` are domain 1 and `layer` is domain 2 (the numbering of `rect_split_x`).
/// Boundaries: `y0` is boundary 2 and `y1` is boundary 5 (the Floquet faces of `common::block`).
struct Table;
impl LabelResolver for Table {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        match r {
            "guide" | "air" => Ok(vec![1]),
            "layer" => Ok(vec![2]),
            _ => Err(Error::invalid(format!("no region {r}"))),
        }
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        match s {
            "y0" => Ok(vec![2]),
            "y1" => Ok(vec![5]),
            _ => Err(Error::invalid(format!("no selector {s}"))),
        }
    }
}

fn bits(r: &SolveResult) -> Vec<Vec<u64>> {
    r.solutions.iter().flat_map(|s| s.fields.iter().map(|f| f.values.iter().map(|v| v.to_bits()).collect())).collect()
}

#[test]
fn the_envelope_solves_through_the_trait_and_agrees_with_the_typed_problem() {
    // modes of a PEC guide with an absorbing layer: once through `Physics::solve` (names resolved by the table), once through the typed problem
    // with the numbers written out; the two must give the same bits
    let mesh = femcore::testkit::rect_split_x(0.03, 0.015, 0.02, 0.0025);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();
    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({ "freq": 12e9, "n_target": 1.5, "n_modes": 3 });
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 4.0)].into()));
    spec.materials.insert("layer".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    spec.bcs.push(BcSpec { on: "layer".into(), kind: "pml".into(), params: json!({ "x": "plus" }) });
    let table = Table;
    let mut ctx = Ctx::new(&table);
    let r = Ems.solve(&exchange, &spec, &mut ctx).unwrap();
    assert_eq!((r.physics.as_str(), r.study.as_str(), r.order, r.dim), ("ems", "modes", 2, 2));
    let mut p = common::emw();
    p.study = nyx_ems::Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 3 };
    p.materials.set(1, "eps_r", 4.0);
    p.materials.set(2, "eps_r", 1.0);
    p.materials.pml.insert(2, femcore::props::Pml { x: femcore::props::PmlDir::Plus, ..Default::default() });
    let t = common::solve(&mesh, &p).unwrap();
    assert_eq!(bits(&r), bits(&t));
    assert_eq!(r.solutions[0].eigenvalue, t.solutions[0].eigenvalue);
    // the result is plain data: it survives a round trip through JSON
    let back: SolveResult = serde_json::from_str(&serde_json::to_string(&r).unwrap()).unwrap();
    assert_eq!(bits(&back), bits(&r));
    // and the sizing hint carries the order
    assert_eq!(Ems.sizing(&spec).order, spec.order);
}

#[test]
fn a_floquet_pair_in_the_envelope_is_the_same_cavity_as_the_typed_one() {
    let (a, l, b) = (0.03, 0.01, 0.015);
    let mesh = common::block(a, l, b, 4, 2, 2);
    let exchange = Mesh::from_solver_mesh(&mesh).unwrap();
    let phase = 0.8;
    let f_te10 = 0.5 * nyx_ems::C0 * ((1.0 / a).powi(2) + (phase / l / std::f64::consts::PI).powi(2)).sqrt();
    let mut spec = ProblemSpec::new("ems", 2, "cavity");
    spec.study.params = json!({ "n_modes": 3, "f_search": 0.95 * f_te10 });
    spec.bcs.push(BcSpec { on: "y0".into(), kind: "periodic".into(), params: json!({ "dest": "y1", "phase": phase }) });
    let table = Table;
    let r = Ems.solve(&exchange, &spec, &mut Ctx::new(&table)).unwrap();
    let mut p = nyx_ems::Problem::new(nyx_ems::Study::Cavity { n_modes: 3, f_search: 0.95 * f_te10 }, 2);
    p.periodic = vec![femcore::props::Periodic { source: vec![2], dest: vec![5], phase }];
    let t = common::solve(&mesh, &p).unwrap();
    assert_eq!(bits(&r), bits(&t));
    assert!(r.solutions.iter().any(|s| (s.value.unwrap() / f_te10 - 1.0).abs() < 5e-3));
    // the Floquet pair of the cavity study is refused in a modes study
    let mut bad = ProblemSpec::new("ems", 2, "modes");
    bad.study.params = json!({ "freq": 12e9 });
    bad.bcs.push(BcSpec { on: "y0".into(), kind: "periodic".into(), params: json!({ "dest": "y1" }) });
    let flat = Mesh::from_solver_mesh(&femcore::testkit::rect(0.03, 0.015, 4, 2)).unwrap();
    let e = Ems.solve(&flat, &bad, &mut Ctx::new(&table)).unwrap_err();
    assert_eq!(e.code(), "invalid");
    assert!(e.to_string().contains("belongs to the cavity study"), "{e}");
}

#[test]
fn the_envelope_is_refused_where_it_is_wrong() {
    let mesh = Mesh::from_solver_mesh(&femcore::testkit::rect(0.03, 0.015, 6, 3)).unwrap();
    let table = Table;
    // an unknown region
    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({ "freq": 12e9 });
    spec.materials.insert("nowhere".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    assert!(Ems.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().to_string().contains("no region nowhere"));
    // a named material needs the material files of phase 6
    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({ "freq": 12e9 });
    spec.materials.insert("guide".into(), MaterialAssign::Named("silicon".into()));
    assert_eq!(Ems.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().code(), "unsupported");
    // the study names the kernel: a cavity study on a 2D mesh, and modes on a 3D one, are refused by name, not guessed from the mesh
    let mut cavity = ProblemSpec::new("ems", 2, "cavity");
    cavity.study.params = json!({ "f_search": 10e9 });
    let e = Ems.solve(&mesh, &cavity, &mut Ctx::new(&table)).unwrap_err();
    assert_eq!(e.code(), "invalid");
    assert!(e.to_string().contains("'cavity' needs a 3D mesh"), "{e}");
    let mesh3 = Mesh::from_solver_mesh(&common::block(0.03, 0.01, 0.015, 3, 2, 2)).unwrap();
    let mut modes = ProblemSpec::new("ems", 2, "modes");
    modes.study.params = json!({ "freq": 12e9 });
    let e = Ems.solve(&mesh3, &modes, &mut Ctx::new(&table)).unwrap_err();
    assert!(e.to_string().contains("'modes' needs a 2D mesh"), "{e}");
    // a layer along y does not exist on a cross-section
    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({ "freq": 12e9 });
    spec.bcs.push(BcSpec { on: "guide".into(), kind: "pml".into(), params: json!({ "y": "plus" }) });
    let e = Ems.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err();
    assert!(e.to_string().contains("out of the plane of a 2D cross-section"), "{e}");
    // a NaN in a material is refused by name (the envelope cannot hold one, the typed problem can)
    let mut p = common::emw();
    p.materials.set(1, "eps_r", f64::NAN);
    let e = nyx_ems::run(&mesh, &p, &femcore::progress::Progress::new()).unwrap_err();
    assert_eq!(e.code(), "non_finite");
    assert!(e.to_string().contains("eps_r"), "{e}");
    // a cancelled run says so
    let mut spec = ProblemSpec::new("ems", 2, "modes");
    spec.study.params = json!({ "freq": 12e9 });
    let mut ctx = Ctx::new(&table);
    ctx.progress.cancel();
    assert_eq!(Ems.solve(&mesh, &spec, &mut ctx).unwrap_err().code(), "cancelled");
}

#[test]
fn kernel_stage_api_is_public() {
    // ADR 9: naming the functions makes removing or renaming one a compile error here
    let _ = (nyx_ems::assemble_modes, nyx_ems::solve_mode_pairs, nyx_ems::mode_fields_cells);
    let _ = (nyx_ems::assemble_modes_complex, nyx_ems::solve_mode_pairs_complex, nyx_ems::mode_fields_complex_cells);
    let _ = (nyx_ems::assemble_axisym, nyx_ems::solve_axisym_pairs, nyx_ems::axisym_fields_cells);
    let _ = (nyx_ems::assemble_emw_eigen, nyx_ems::solve_emw_modes, nyx_ems::emw_eigen_fields);
    // the one-call kernels of the studies, the run and the decoder
    let _ = (
        nyx_ems::solve_waveguide,
        nyx_ems::solve_waveguide_complex,
        nyx_ems::solve_axisymmetric,
        nyx_ems::solve_cavity,
    );
    let _ = (nyx_ems::run, nyx_ems::run_analysed, nyx_ems::decode, nyx_ems::solve_optical);
}

/// The shape of a result, without numbers: what a front-end, the store and the analyses may rely on (names, units, components, the per-cell
/// fields, the scalars, the number of solutions).
fn result_shape(r: &SolveResult) -> serde_json::Value {
    json!({
        "physics": r.physics,
        "study": r.study,
        "solutions": r.solutions.iter().map(|s| json!({
            "label_prefix": s.label.split(':').next(),
            "has_value": s.value.is_some(),
            "has_eigenvalue": s.eigenvalue.is_some(),
            "fields": s.fields.iter().map(|f| json!({ "name": f.name, "unit": f.unit, "components": f.components })).collect::<Vec<_>>(),
            "cell_fields": s.cell_fields.iter().map(|f| f.name.clone()).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
        "scalars": r.scalars.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>(),
    })
}

#[test]
fn result_shape_exemplar_is_stable() {
    // one small problem per kind of kernel; the exemplar `physics_ems_result.v1.json` lists their shapes
    use femcore::props::{Pml, PmlDir};
    use nyx_ems::Study;
    let run = |mesh: &femcore::mesh::SimplexMesh, p: &nyx_ems::Problem| {
        let exchange = Mesh::from_solver_mesh(mesh).unwrap();
        result_shape(&nyx_ems::run(&exchange, p, &femcore::progress::Progress::new()).unwrap())
    };
    let guide = femcore::testkit::rect(0.03, 0.015, 8, 4);
    let mut real = common::emw();
    real.study = Study::Modes { freq: 12e9, n_target: 0.9, evanescent: false, n_modes: 2 };
    let split = femcore::testkit::rect_split_x(0.03, 0.015, 0.02, 0.0025);
    let mut complex = common::emw();
    complex.study = Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 2 };
    complex.materials.set(1, "eps_r", 4.0);
    complex.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Pml::default() });
    let coax = femcore::testkit::parts_mesh(&[femcore::testkit::part(0.10, 0.0, 0.20, 0.15, 1, Some(0.03))], 0.03);
    let axi = nyx_ems::Problem::new(Study::AxisymmetricModes { f_search: 1.0e9, azimuthal_order: 0, n_modes: 1 }, 2);
    let boxm = common::block(0.04, 0.03, 0.02, 3, 2, 2);
    let cavity = nyx_ems::Problem::new(Study::Cavity { n_modes: 2, f_search: 6.0e9 }, 2);
    let got = json!({
        "modes_real": run(&guide, &real),
        "modes_complex": run(&split, &complex),
        "axisymmetric_modes": run(&coax, &axi),
        "cavity": run(&boxm, &cavity),
    });
    contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), "physics_ems_result.v1.json", &got);
}
