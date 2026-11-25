// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Inline diagonal and symmetric tensor material shorthand tests (FR-7): `PropValue`, `MaterialAssign::PropsTensor`,
//! structural validation (tensor-only kind check, finiteness, symmetry tolerance), decoding into `Materials::aniso`,
//! and JSON Schema generation.

use femcore::material::Tensor;
use femcore::mesh::SimplexMesh;
use femcore::props::TensorKind;
use nyx_physics::decode::materials_by_label;
use nyx_physics::schema::physics_schema;
use nyx_physics::*;
use serde_json::json;

struct TensorMock;

impl Physics for TensorMock {
    fn id(&self) -> &'static str {
        "tst"
    }
    fn name(&self) -> &'static str {
        "Tensor mock physics"
    }
    fn studies(&self) -> Vec<StudyDesc> {
        vec![StudyDesc {
            kind: "modes".into(),
            name: "Modes".into(),
            doc: "a study".into(),
            dims: Dims::Both,
            params: vec![ParamDesc::number("freq", "Frequency", "Hz", 1e9, "").required().range(Some(1.0), Some(1e18))],
        }]
    }
    fn bc_kinds(&self) -> Vec<BcDesc> {
        vec![]
    }
    fn properties(&self) -> Vec<PropDesc> {
        vec![
            PropDesc::new("E", "Young's modulus", "Pa", PropKind::Scalar, &[], None),
            PropDesc::new("mu_r", "Relative permeability", "", PropKind::Scalar, &[], Some(1.0)),
            PropDesc::new("eps_r", "Permittivity", "", PropKind::Tensor, &["modes"], Some(1.0)),
            PropDesc::new("k", "Thermal conductivity", "W/(m K)", PropKind::Tensor, &[], None),
            PropDesc::new("alpha", "Thermal expansion", "1/K", PropKind::Tensor, &[], Some(0.0)),
            PropDesc::new("piezo_e", "Piezoelectric coupling", "C/m^2", PropKind::Matrix3x6, &[], None),
        ]
    }
    fn params(&self) -> Vec<ParamDesc> {
        vec![]
    }
    fn fields(&self, _study: &str) -> Vec<FieldDesc> {
        vec![FieldDesc::new("E", "V/m", &["x", "y", "z"], "field")]
    }
    fn exports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn imports(&self) -> Vec<FieldLink> {
        vec![]
    }
    fn validate(&self, _spec: &ProblemSpec, _mesh: &Mesh) -> Result<(), Error> {
        Ok(())
    }
    fn sizing(&self, spec: &ProblemSpec, _library: Option<&std::sync::Arc<femcore::material::Library>>) -> Sizing {
        Sizing { order: spec.order, headroom: 1.0, wavelengths: Default::default() }
    }
    fn solve(&self, _mesh: &Mesh, _spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error> {
        ctx.check()?;
        Ok(SolveResult { physics: "tst".into(), ..Default::default() })
    }
}

struct OneDomain;

impl LabelResolver for OneDomain {
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error> {
        match region {
            "core" => Ok(vec![1]),
            other => Err(Error::invalid(format!("no region '{other}'"))),
        }
    }
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error> {
        Err(Error::invalid(format!("no boundary '{selector}'")))
    }
}

fn base_spec() -> ProblemSpec {
    let mut s = ProblemSpec::new("tst", 2, "modes");
    s.study.params = json!({ "freq": 2e14 });
    s
}

fn one_tet_mesh() -> SimplexMesh {
    Mesh::new(
        3,
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        vec![0, 1, 2, 3],
        vec![1],
        vec![0, 1, 2],
        vec![1],
    )
    .discretise()
    .unwrap()
}

#[test]
fn scalar_only_map_stays_props_and_tensor_shorthand_deserializes_as_props_tensor() {
    let scalar_spec: ProblemSpec = serde_json::from_value(json!({
        "physics": "tst",
        "order": 2,
        "study": { "kind": "modes", "params": { "freq": 2e14 } },
        "materials": {
            "core": { "eps_r": 12.08, "mu_r": 1.0 }
        }
    }))
    .unwrap();
    assert!(matches!(scalar_spec.materials.get("core"), Some(MaterialAssign::Props(_))));

    let diag_spec: ProblemSpec = serde_json::from_value(json!({
        "physics": "tst",
        "order": 2,
        "study": { "kind": "modes", "params": { "freq": 2e14 } },
        "materials": {
            "core": { "eps_r": [9.4, 9.4, 11.6], "mu_r": 1.0 }
        }
    }))
    .unwrap();
    match diag_spec.materials.get("core").unwrap() {
        MaterialAssign::PropsTensor(map) => {
            assert_eq!(map.get("eps_r"), Some(&PropValue::Diag([9.4, 9.4, 11.6])));
            assert_eq!(map.get("mu_r"), Some(&PropValue::Scalar(1.0)));
        }
        other => panic!("expected PropsTensor, got {other:?}"),
    }

    let mat = [[9.4, 0.5, -0.2], [0.5, 10.1, 0.3], [-0.2, 0.3, 11.6]];
    let mat_spec: ProblemSpec = serde_json::from_value(json!({
        "physics": "tst",
        "order": 2,
        "study": { "kind": "modes", "params": { "freq": 2e14 } },
        "materials": {
            "core": { "eps_r": mat }
        }
    }))
    .unwrap();
    match mat_spec.materials.get("core").unwrap() {
        MaterialAssign::PropsTensor(map) => {
            assert_eq!(map.get("eps_r"), Some(&PropValue::Matrix(mat)));
        }
        other => panic!("expected PropsTensor, got {other:?}"),
    }
}

#[test]
fn validate_spec_enforces_tensor_kind_finiteness_and_symmetry() {
    // Valid diagonal and symmetric matrix on tensor properties pass
    let mut s = base_spec();
    s.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor(
            [
                ("eps_r".into(), PropValue::Diag([9.4, 9.4, 11.6])),
                ("k".into(), PropValue::Matrix([[150.0, 5.0, 0.0], [5.0, 150.0, 0.0], [0.0, 0.0, 120.0]])),
                ("mu_r".into(), PropValue::Scalar(1.0)),
            ]
            .into(),
        ),
    );
    validate_spec(&TensorMock, &s).unwrap();

    // Tensor value on a scalar property (mu_r or E) is rejected
    let mut bad_kind = base_spec();
    bad_kind
        .materials
        .insert("core".into(), MaterialAssign::PropsTensor([("mu_r".into(), PropValue::Diag([1.0, 1.0, 1.0]))].into()));
    let err = validate_spec(&TensorMock, &bad_kind).unwrap_err().to_string();
    assert!(err.contains("mu_r") && err.contains("scalar"), "{err}");

    // Non-finite component in Diag or Matrix is rejected
    let mut bad_nan = base_spec();
    bad_nan.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor([("eps_r".into(), PropValue::Diag([9.4, f64::NAN, 11.6]))].into()),
    );
    let err = validate_spec(&TensorMock, &bad_nan).unwrap_err().to_string();
    assert!(err.contains("finite"), "{err}");

    // Asymmetric 3x3 matrix is rejected
    let mut bad_asym = base_spec();
    bad_asym.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor(
            [("eps_r".into(), PropValue::Matrix([[9.4, 0.5, 0.0], [-0.5, 9.4, 0.0], [0.0, 0.0, 11.6]]))].into(),
        ),
    );
    let err = validate_spec(&TensorMock, &bad_asym).unwrap_err().to_string();
    assert!(err.contains("symmetric"), "{err}");
}

#[test]
fn decode_populates_aniso_tensors_and_evaluates_cell_tensor() {
    let mesh = one_tet_mesh();
    let ctx = Ctx::new(&OneDomain);

    // 1. Diagonal anisotropic tensor
    let mut s = base_spec();
    s.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor(
            [
                ("eps_r".into(), PropValue::Diag([9.4, 9.8, 11.6])),
                ("k".into(), PropValue::Diag([40.0, 40.0, 30.0])),
                ("alpha".into(), PropValue::Diag([5e-6, 5e-6, 7e-6])),
                ("mu_r".into(), PropValue::Scalar(1.25)),
            ]
            .into(),
        ),
    );
    let m = materials_by_label(&s, &ctx, &TensorMock.properties()).unwrap();
    assert_eq!(m.get(1, "mu_r"), Some(1.25));
    let an = &m.aniso[&1];
    assert_eq!(an.euler_deg, [0.0; 3]);
    assert_eq!(an.permittivity, Some(Tensor::Diagonal { x: 9.4, y: 9.8, z: 11.6 }));
    assert_eq!(an.conductivity, Some(Tensor::Diagonal { x: 40.0, y: 40.0, z: 30.0 }));
    assert_eq!(an.expansion, Some(Tensor::Diagonal { x: 5e-6, y: 5e-6, z: 7e-6 }));
    let ct = m.cell_tensor(&mesh, TensorKind::Permittivity, Some(1.0)).unwrap();
    assert_eq!(ct[0], [[9.4, 0.0, 0.0], [0.0, 9.8, 0.0], [0.0, 0.0, 11.6]]);

    // 2. Full 3x3 symmetric matrix
    let mat = [[9.4, 0.5, -0.2], [0.5, 10.1, 0.3], [-0.2, 0.3, 11.6]];
    let mut s2 = base_spec();
    s2.materials.insert("core".into(), MaterialAssign::PropsTensor([("eps_r".into(), PropValue::Matrix(mat))].into()));
    let m2 = materials_by_label(&s2, &ctx, &TensorMock.properties()).unwrap();
    assert_eq!(
        m2.aniso[&1].permittivity,
        Some(Tensor::Symmetric { xx: 9.4, yy: 10.1, zz: 11.6, xy: 0.5, yz: 0.3, xz: -0.2 })
    );
    let ct2 = m2.cell_tensor(&mesh, TensorKind::Permittivity, Some(1.0)).unwrap();
    assert_eq!(ct2[0], mat);

    // 2b. 6-component symmetric [xx, yy, zz, xy, xz, yz] shorthand
    let mut s2b = base_spec();
    s2b.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor([("eps_r".into(), PropValue::Sym6([9.4, 10.1, 11.6, 0.5, -0.2, 0.3]))].into()),
    );
    validate_spec(&TensorMock, &s2b).unwrap();
    let m2b = materials_by_label(&s2b, &ctx, &TensorMock.properties()).unwrap();
    let ct2b = m2b.cell_tensor(&mesh, TensorKind::Permittivity, Some(1.0)).unwrap();
    assert_eq!(ct2b[0], mat);

    // 3. Isotropic diagonal tensor reduces to scalar
    let mut s3 = base_spec();
    s3.materials.insert(
        "core".into(),
        MaterialAssign::PropsTensor([("eps_r".into(), PropValue::Diag([4.5, 4.5, 4.5]))].into()),
    );
    let m3 = materials_by_label(&s3, &ctx, &TensorMock.properties()).unwrap();
    assert_eq!(m3.get(1, "eps_r"), Some(4.5));
    assert!(!m3.has_tensor(TensorKind::Permittivity));
    let ct3 = m3.cell_tensor(&mesh, TensorKind::Permittivity, Some(1.0)).unwrap();
    assert_eq!(ct3[0], [[4.5, 0.0, 0.0], [0.0, 4.5, 0.0], [0.0, 0.0, 4.5]]);
}

#[test]
fn physics_schema_allows_scalar_diag_and_matrix_for_tensor_properties_only() {
    let schema = physics_schema(&TensorMock);
    let mat_props = &schema["properties"]["materials"]["additionalProperties"]["oneOf"][1]["properties"];
    assert_eq!(mat_props["E"]["type"], json!("number"));
    assert_eq!(mat_props["mu_r"]["type"], json!("number"));
    let eps_one_of = mat_props["eps_r"]["oneOf"].as_array().expect("eps_r should have oneOf schema");
    assert_eq!(eps_one_of.len(), 4);
    assert_eq!(eps_one_of[0]["type"], json!("number"));
    assert_eq!(eps_one_of[1]["minItems"], json!(3));
    assert_eq!(eps_one_of[2]["minItems"], json!(6));
    assert_eq!(eps_one_of[3]["items"]["minItems"], json!(3));
    assert_eq!(mat_props["piezo_e"]["minItems"], json!(3));
    assert_eq!(mat_props["piezo_e"]["items"]["minItems"], json!(6));
    let mat_one_of = schema["properties"]["materials"]["additionalProperties"]["oneOf"].as_array().unwrap();
    assert_eq!(mat_one_of.len(), 3);
    assert_eq!(mat_one_of[2]["required"], json!(["material", "euler_deg"]));
}

#[test]
fn inline_piezo_e_matrix3x6_validates_decodes_and_rejects_scalar() {
    let e_mat =
        [[0.0, 0.0, 0.0, 0.0, -0.48, 0.0], [0.0, 0.0, 0.0, -0.48, 0.0, 0.0], [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0]];
    let spec: ProblemSpec = serde_json::from_value(json!({
        "physics": "tst",
        "order": 2,
        "study": { "kind": "modes", "params": { "freq": 2e9 } },
        "materials": {
            "core": { "E": 330e9, "eps_r": 9.5, "piezo_e": e_mat }
        }
    }))
    .unwrap();
    validate_spec(&TensorMock, &spec).unwrap();
    let ctx = Ctx::new(&OneDomain);
    let decoded = materials_by_label(&spec, &ctx, &TensorMock.properties()).unwrap();
    assert_eq!(decoded.aniso[&1].piezo_e, Some(e_mat));

    // Scalar/Props or 3x3 for piezo_e must fail
    let mut bad_scalar = base_spec();
    bad_scalar.materials.insert("core".into(), MaterialAssign::Props([("piezo_e".into(), 1.55)].into()));
    assert!(validate_spec(&TensorMock, &bad_scalar).is_err());

    let mut bad_mat = base_spec();
    bad_mat
        .materials
        .insert("core".into(), MaterialAssign::PropsTensor([("eps_r".into(), PropValue::Matrix3x6(e_mat))].into()));
    assert!(validate_spec(&TensorMock, &bad_mat).is_err());
}
