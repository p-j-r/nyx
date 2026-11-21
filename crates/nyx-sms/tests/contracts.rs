// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Contract of `nyx-sms` (ADR 19, ADR 25, ADR 42): its descriptors are pinned by `docs/architecture/contracts/physics_sms.v1.json`, an example problem
//! passes the structural validation; the envelope solves through `Physics::solve` and agrees with the typed problem; the shape of a result is pinned by
//! `physics_sms_result.v1.json`; the stage functions are public.
//!
//! Rewrite the exemplar after a deliberate change (with a `VERSIONS.md` line):
//! `NYX_WRITE_CONTRACT=1 cargo test --release -p nyx-sms --test contracts`
use nyx_physics::*;
use nyx_sms::Sms;
use serde_json::json;

fn example() -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, "eigenfrequency");
    s.study.params = json!({ "n_modes": 6, "f_search": 1e9 });
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

fn axisymmetric_example() -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, "axisymmetric_eigenfrequency");
    s.study.params = json!({ "n_modes": 6, "f_search": 715e6, "azimuthal_order": 0 });
    s.materials.insert(
        "disk".into(),
        MaterialAssign::Props([("E".to_string(), 160e9), ("nu".to_string(), 0.22), ("rho".to_string(), 2330.0)].into()),
    );
    s.bcs.push(BcSpec { on: "outer".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0, "uy": 0.0 }) });
    s.bcs.push(BcSpec {
        on: "substrate_pml".into(),
        kind: "pml".into(),
        params: json!({ "x": "plus", "z": "minus", "strength": 40.0, "order": 1.0 }),
    });
    s
}

fn modes_example() -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, "modes");
    s.study.params = json!({ "n_modes": 6, "f_search": 7.45e9, "wavevector": 2.3468e7 });
    s.materials.insert(
        "core".into(),
        MaterialAssign::Props([("E".to_string(), 169e9), ("nu".to_string(), 0.22), ("rho".to_string(), 2329.0)].into()),
    );
    s.bcs.push(BcSpec { on: "outer".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0, "uy": 0.0 }) });
    s.bcs.push(BcSpec {
        on: "substrate_pml".into(),
        kind: "pml".into(),
        params: json!({ "x": "plus", "z": "minus", "strength": 3.0, "order": 2.0 }),
    });
    s
}

#[test]
fn the_axisymmetric_study_example_passes_validation() {
    validate_spec(&Sms, &axisymmetric_example()).unwrap();
    validate_spec(&Sms, &modes_example()).unwrap();
    // its parameters are checked against the descriptors: a fractional azimuthal order and an unknown parameter are refused
    let mut bad = axisymmetric_example();
    bad.study.params = json!({ "n_modes": 6, "f_search": 715e6, "azimuthal_order": 0.5 });
    assert!(validate_spec(&Sms, &bad).is_err());
    let mut bad = axisymmetric_example();
    bad.study.params = json!({ "n_modes": 6, "wave_number": 1.0 });
    assert!(validate_spec(&Sms, &bad).is_err());
    let mut bad = modes_example();
    bad.study.params = json!({ "n_modes": 6, "azimuthal_order": 0 });
    assert!(validate_spec(&Sms, &bad).is_err());
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

/// A table of names for the tests (the numbering of `femcore::testkit::rect`): region `bar` is domain 1, `left` 1, `bottom` 2, `top` 3, `right` 4.
struct Table;
impl LabelResolver for Table {
    fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
        match r {
            "bar" => Ok(vec![1]),
            _ => Err(Error::invalid(format!("no region {r}"))),
        }
    }
    fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
        match s {
            "left" => Ok(vec![1]),
            "bottom" => Ok(vec![2]),
            "top" => Ok(vec![3]),
            "right" => Ok(vec![4]),
            _ => Err(Error::invalid(format!("no selector {s}"))),
        }
    }
}

fn exchange(mesh: &femcore::mesh::SimplexMesh) -> Mesh {
    Mesh::from_solver_mesh(mesh).unwrap()
}

fn bar_spec(study: &str) -> ProblemSpec {
    let mut s = ProblemSpec::new("sms", 2, study);
    s.materials.insert(
        "bar".into(),
        MaterialAssign::Props([("E".to_string(), 200e9), ("nu".to_string(), 0.3), ("rho".to_string(), 7850.0)].into()),
    );
    s.bcs.push(BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0, "uz": 0.0 }) });
    s
}

fn typed_bar(study: nyx_sms::Study) -> nyx_sms::Problem {
    let mut p = nyx_sms::Problem::new(study, 2);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        p.materials.set(1, k, v);
    }
    p.bcs = vec![nyx_sms::Bc {
        labels: vec![1],
        kind: nyx_sms::BcKind::Fixed { fixed: [true, false, true], values: [0.0; 3] },
    }];
    p
}

#[test]
fn the_envelope_solves_through_the_trait_and_agrees_with_the_typed_problem() {
    // a loaded cantilever, once through `Physics::solve` (names resolved by the table) and once through the typed problem with the numbers
    // written out: the two must give the same bits
    let mesh = femcore::testkit::rect(1.0, 0.2, 10, 3);
    let mut spec = bar_spec("static");
    spec.bcs.push(BcSpec { on: "right".into(), kind: "load".into(), params: json!({ "traction": [0.0, 0.0, -1e5] }) });
    spec.params = json!({ "plane": "stress" });
    let table = Table;
    let r = Sms.solve(&exchange(&mesh), &spec, &mut Ctx::new(&table)).unwrap();
    assert_eq!((r.physics.as_str(), r.study.as_str(), r.order, r.dim), ("sms", "static", 2, 2));
    let mut p = typed_bar(nyx_sms::Study::Static);
    p.plane = nyx_sms::Plane::Stress;
    p.bcs.push(nyx_sms::Bc { labels: vec![4], kind: nyx_sms::BcKind::Load([0.0, 0.0, -1e5]) });
    let t = nyx_sms::solve_solid(&mesh, &p).unwrap();
    assert_eq!(r.solutions[0].field("u").unwrap().values, t.solutions[0].field("u").unwrap().values);
    assert_eq!(r.solutions[0].flat("von Mises"), t.solutions[0].flat("von Mises"));
    assert!(r.scalar("max |u|").unwrap() > 0.0);
    // the result is plain data: it survives a round trip through JSON
    let back: SolveResult = serde_json::from_str(&serde_json::to_string(&r).unwrap()).unwrap();
    assert_eq!(back.solutions[0].field("u").unwrap().values, r.solutions[0].field("u").unwrap().values);
    // and the sizing hint carries the order
    assert_eq!(Sms.sizing(&spec, None).order, spec.order);
}

#[test]
fn absorbing_layers_and_periodic_pairs_come_through_the_envelope() {
    let mesh = femcore::testkit::rect(1.0, 0.1, 20, 2);
    let table = Table;
    // a perfectly matched layer on the region: complex frequencies, the same as the typed problem
    let mut spec = bar_spec("eigenfrequency");
    spec.study.params = json!({ "n_modes": 3 });
    spec.params = json!({ "plane": "stress" });
    spec.bcs.push(BcSpec {
        on: "bar".into(),
        kind: "pml".into(),
        params: json!({ "x": "plus", "strength": 0.1, "order": 0.0 }),
    });
    let r = Sms.solve(&exchange(&mesh), &spec, &mut Ctx::new(&table)).unwrap();
    let mut p = typed_bar(nyx_sms::Study::Eigenfrequency { n_modes: 3, f_search: 0.0 });
    p.plane = nyx_sms::Plane::Stress;
    p.materials.pml.insert(
        1,
        femcore::props::Pml { x: femcore::props::PmlDir::Plus, strength: 0.1, order: 0.0, ..Default::default() },
    );
    let t = nyx_sms::solve_solid(&mesh, &p).unwrap();
    let f = |r: &SolveResult| r.solutions.iter().map(|s| s.value.unwrap()).collect::<Vec<_>>();
    assert_eq!(f(&r), f(&t));
    assert!(r.scalar("Q1").is_some());
    // a Floquet pair: `on` is the source, `dest` is a selector of its own
    let mut spec = ProblemSpec::new("sms", 2, "eigenfrequency");
    spec.study.params = json!({ "n_modes": 3, "f_search": 1.0e6 });
    spec.materials.insert(
        "bar".into(),
        MaterialAssign::Props([("E".to_string(), 70e9), ("nu".to_string(), 0.3), ("rho".to_string(), 2700.0)].into()),
    );
    spec.bcs.push(BcSpec {
        on: "left".into(),
        kind: "periodic".into(),
        params: json!({ "dest": "right", "phase": 1.2 }),
    });
    let mesh = femcore::testkit::rect(1.0e-3, 0.3e-3, 8, 3);
    let r = Sms.solve(&exchange(&mesh), &spec, &mut Ctx::new(&table)).unwrap();
    let mut p = nyx_sms::Problem::new(nyx_sms::Study::Eigenfrequency { n_modes: 3, f_search: 1.0e6 }, 2);
    for (k, v) in [("E", 70e9), ("nu", 0.3), ("rho", 2700.0)] {
        p.materials.set(1, k, v);
    }
    p.periodic = vec![femcore::props::Periodic { source: vec![1], dest: vec![4], phase: 1.2 }];
    let t = nyx_sms::solve_solid(&mesh, &p).unwrap();
    assert_eq!(f(&r), f(&t));
}

#[test]
fn an_imported_temperature_is_the_thermal_load() {
    // the temperature another physics exported (`dfs.T`, the role `temperature`, one value per dof) loads the solid like the same uniform one
    let mesh = femcore::testkit::rect(1.0, 0.1, 10, 2);
    let table = Table;
    let mut spec = bar_spec("static");
    spec.materials.insert(
        "bar".into(),
        MaterialAssign::Props([("E".to_string(), 70e9), ("nu".to_string(), 0.33), ("alpha".to_string(), 23e-6)].into()),
    );
    spec.bcs = vec![
        BcSpec { on: "left".into(), kind: "fixed".into(), params: json!({ "ux": 0.0 }) },
        BcSpec { on: "bottom".into(), kind: "fixed".into(), params: json!({ "uz": 0.0 }) },
    ];
    spec.params = json!({ "plane": "stress", "thermal": "uniform", "temperature": 393.15, "t_ref": 293.15 });
    let uniform = Sms.solve(&exchange(&mesh), &spec, &mut Ctx::new(&table)).unwrap();
    spec.params = json!({ "plane": "stress", "thermal": "imported", "t_ref": 293.15 });
    // no temperature given
    let e = Sms.solve(&exchange(&mesh), &spec, &mut Ctx::new(&table)).unwrap_err();
    assert_eq!(e.code(), "invalid");
    assert!(e.to_string().contains("temperature"), "{e}");
    let mut ctx = Ctx::new(&table);
    ctx.inputs.insert("temperature".into(), vec![393.15; mesh.n_dofs(2)]);
    let imported = Sms.solve(&exchange(&mesh), &spec, &mut ctx).unwrap();
    assert_eq!(imported.solutions[0].field("u").unwrap().values, uniform.solutions[0].field("u").unwrap().values);
    // a temperature on another mesh is refused with the old message
    let mut ctx = Ctx::new(&table);
    ctx.inputs.insert("temperature".into(), vec![393.15; 3]);
    let e = Sms.solve(&exchange(&mesh), &spec, &mut ctx).unwrap_err();
    assert!(e.to_string().contains("no heat-transfer result on this mesh"), "{e}");
}

#[test]
fn the_envelope_is_refused_where_it_is_wrong() {
    let mesh = exchange(&femcore::testkit::rect(1.0, 0.5, 4, 2));
    let table = Table;
    // an unknown region
    let mut spec = bar_spec("static");
    spec.materials.insert("nowhere".into(), MaterialAssign::Props([("E".to_string(), 1.0)].into()));
    assert!(Sms.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().to_string().contains("no region nowhere"));
    // conditions that belong to the eigenfrequency study, in a static one
    for (kind, on, params) in [
        ("pml", "bar", json!({ "x": "plus" })),
        ("low_reflecting", "right", serde_json::Value::Null),
        ("periodic", "left", json!({ "dest": "right" })),
    ] {
        let mut spec = bar_spec("static");
        spec.bcs.push(BcSpec { on: on.into(), kind: kind.into(), params });
        let e = Sms.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err();
        assert_eq!(e.code(), "invalid", "{kind}");
        assert!(e.to_string().contains("eigenfrequency study"), "{kind}: {e}");
    }
    // a boundary load in an eigenfrequency study
    let mut spec = bar_spec("eigenfrequency");
    spec.bcs.push(BcSpec { on: "right".into(), kind: "load".into(), params: json!({ "traction": [1.0, 0.0, 0.0] }) });
    assert_eq!(Sms.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().code(), "invalid");
    // a named material needs the material files of phase 6
    let mut spec = bar_spec("static");
    spec.materials.insert("bar".into(), MaterialAssign::Named("silicon".into()));
    assert_eq!(Sms.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err().code(), "unsupported");
    // a static problem needs a fixed boundary
    let mut spec = bar_spec("static");
    spec.bcs.clear();
    let e = Sms.solve(&mesh, &spec, &mut Ctx::new(&table)).unwrap_err();
    assert!(e.to_string().contains("needs a fixed boundary"), "{e}");
    // a cancelled run says so
    let mut ctx = Ctx::new(&table);
    ctx.progress.cancel();
    assert_eq!(Sms.solve(&mesh, &bar_spec("static"), &mut ctx).unwrap_err().code(), "cancelled");
}

#[test]
fn kernel_stage_api_is_public() {
    // ADR 9: naming the functions makes removing or renaming one a compile error here (contract `kernel_stages`)
    use nyx_sms::{axisym, cross_section, floquet, lrb, pml, solid, waveguide};
    let _ = (solid::setup, solid::assemble_solid, solid::assemble_static, solid::solve_static, solid::solve_real_modes);
    let _ = (pml::assemble_solid_pml, pml::solve_pml_modes, pml::pml_fields);
    let _ = (floquet::assemble_floquet, floquet::solve_floquet_modes, floquet::floquet_fields);
    let _ = (lrb::assemble_lrb, lrb::solve_lrb_modes, lrb::lrb_fields);
    let _ = (
        cross_section::assemble_cross_section,
        cross_section::solve_cross_section_modes,
        cross_section::cross_section_fields,
        cross_section::solve_axisymmetric,
        cross_section::solve_waveguide_modes,
    );
    let _ = (axisym::assemble_stiffness_axisym, axisym::assemble_mass_axisym, axisym::average_about_axis);
    let _ = (waveguide::assemble_stiffness_waveguide, waveguide::mode_energy_waveguide);
    let _ = (
        nyx_sms::solve_solid,
        nyx_sms::solve_solid_axisym,
        nyx_sms::solve_solid_waveguide,
        nyx_sms::solve_waveguide_modes,
        nyx_sms::run,
        nyx_sms::decode,
    );
}

/// The shape of a result, without numbers: what a front-end, the store and the analyses may rely on (names, units, components, the per-cell
/// fields, the scalars, the number of solutions).
fn result_shape(r: &SolveResult) -> serde_json::Value {
    json!({
        "physics": r.physics,
        "study": r.study,
        "solutions": r.solutions.iter().map(|s| json!({
            "label": s.label,
            "has_value": s.value.is_some(),
            "fields": s.fields.iter().map(|f| json!({ "name": f.name, "unit": f.unit, "components": f.components })).collect::<Vec<_>>(),
            "flat": s.flat_fields().iter().map(|f| f.0.clone()).collect::<Vec<_>>(),
            "cell_fields": s.cell_fields.iter().map(|f| f.name.clone()).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
        "scalars": r.scalars.iter().map(|(n, _)| n.clone()).collect::<Vec<_>>(),
    })
}

#[test]
fn result_shape_exemplar_is_stable() {
    // one small problem per study and per variant of the eigenfrequency kernel (real, absorbing layer, low-reflecting boundary, Floquet pair);
    // the exemplar `physics_sms_result.v1.json` lists their shapes (the first solution of each: the others have the same shape)
    use femcore::props::{Periodic, Pml, PmlDir};
    use nyx_sms::{Bc, BcKind, Problem, Study};
    let mesh = femcore::testkit::rect(1.0, 0.1, 12, 2);
    let ex = exchange(&mesh);
    let run = |p: &Problem| {
        let mut r = nyx_sms::run(&ex, p, &femcore::progress::Progress::new()).unwrap();
        r.solutions.truncate(1);
        result_shape(&r)
    };
    let clamp = Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } };
    let mut stat = Problem::new(Study::Static, 2);
    for (k, v) in [("E", 200e9), ("nu", 0.3), ("rho", 7850.0)] {
        stat.materials.default.insert(k.into(), v);
    }
    stat.bcs = vec![clamp.clone(), Bc { labels: vec![4], kind: BcKind::Load([0.0, 0.0, -1e5]) }];
    let mut eigen = stat.clone();
    eigen.study = Study::Eigenfrequency { n_modes: 2, f_search: 0.0 };
    eigen.bcs = vec![clamp.clone()];
    let mut with_pml = eigen.clone();
    with_pml.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    let mut with_lrb = eigen.clone();
    with_lrb.study = Study::Eigenfrequency { n_modes: 2, f_search: 2.0e3 };
    with_lrb.bcs.push(Bc { labels: vec![4], kind: BcKind::LowReflecting });
    let mut floquet = eigen.clone();
    floquet.bcs.clear();
    floquet.study = Study::Eigenfrequency { n_modes: 2, f_search: 1.0e3 };
    floquet.periodic = vec![Periodic { source: vec![1], dest: vec![4], phase: 1.2 }];
    let mut axi = eigen.clone();
    axi.study = Study::AxisymmetricEigenfrequency { n_modes: 2, f_search: 0.0, azimuthal_order: 0 };
    axi.bcs = vec![clamp.clone()];
    let mut axi_pml = axi.clone();
    axi_pml.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    let mut wg = eigen.clone();
    wg.study = Study::Modes { n_modes: 2, f_search: 0.0, wavevector: 10.0 };
    wg.bcs = vec![clamp.clone()];
    let mut wg_pml = wg.clone();
    wg_pml.materials.pml.insert(1, Pml { x: PmlDir::Plus, strength: 0.1, order: 0.0, ..Pml::default() });
    // the half plane must lie in x >= 0: the same bar moved to r = 1 .. 2
    let ring = {
        let flat = femcore::testkit::rect(1.0, 0.1, 12, 2);
        let mut pts = flat.points.clone();
        for i in 0..flat.n_points() {
            pts[2 * i] += 1.0;
        }
        femcore::mesh::SimplexMesh::new(
            2,
            pts,
            flat.cells.clone(),
            flat.cell_labels.clone(),
            flat.facets.clone(),
            flat.facet_labels.clone(),
        )
    };
    let ring_ex = exchange(&ring);
    let run_ring = |p: &Problem| {
        let mut r = nyx_sms::run(&ring_ex, p, &femcore::progress::Progress::new()).unwrap();
        r.solutions.truncate(1);
        result_shape(&r)
    };
    let got = json!({
        "static": run(&stat),
        "eigenfrequency": run(&eigen),
        "eigenfrequency_pml": run(&with_pml),
        "eigenfrequency_low_reflecting": run(&with_lrb),
        "eigenfrequency_floquet": run(&floquet),
        "axisymmetric_eigenfrequency": run_ring(&axi),
        "axisymmetric_eigenfrequency_pml": run_ring(&axi_pml),
        "modes": run(&wg),
        "modes_pml": run(&wg_pml),
    });
    contract::check_exemplar(env!("CARGO_MANIFEST_DIR"), "physics_sms_result.v1.json", &got);
}
