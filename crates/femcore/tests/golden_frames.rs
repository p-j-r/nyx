//! Frozen reference output of the 2D kernels, taken through `solve::solve` before a kernel is ported to the wafer axes.
//!
//! Each case is a wafer-frame problem. The golden file `tests/golden_frames/<case>.json` holds the names of every field of every solution, the
//! scalars and effective indices, and for every field a few numbers (sum of |v|, max |v|, the values at fixed dofs). A ported kernel must
//! reproduce them to round-off. `NYX_WRITE_GOLDEN=1 cargo test --test golden_frames` writes the files (do that only on a kernel that is known
//! to be right, i.e. before porting it).
mod common;
use femcore::material::*;
use femcore::problem::*;
use serde_json::{json, Value};
use std::collections::BTreeMap;

const C0: f64 = 299_792_458.0;

fn record(r: &SolveResult) -> Value {
    let mut sols = vec![];
    for s in &r.solutions {
        let mut fields = BTreeMap::new();
        for f in &s.fields {
            let v = &f.values;
            let n = v.len();
            let l1: f64 = v.iter().map(|x| x.abs()).sum();
            let mx = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            fields.insert(f.name.clone(), json!({"n": n, "l1": l1, "max": mx, "probe": [v[0], v[n / 3], v[n / 2], v[n - 1]]}));
        }
        sols.push(json!({"neff": s.neff.map(|(a, b)| json!([a, b])), "value": s.value, "fields": fields}));
    }
    let mut scalars = BTreeMap::new();
    for (k, v) in &r.scalars {
        scalars.insert(k.clone(), match v {
            Scalar::Num(x) => json!(x),
            Scalar::Text(t) => json!(t),
        });
    }
    json!({"solutions": sols, "scalars": scalars})
}

fn close(a: f64, b: f64, scale: f64, tol: f64) -> bool {
    (a - b).abs() <= tol * scale.max(1e-300)
}

fn compare(path: &str, want: &Value, got: &Value, scale: f64, tol: f64, out: &mut Vec<String>) {
    match (want, got) {
        (Value::Number(a), Value::Number(b)) => {
            let (a, b) = (a.as_f64().unwrap(), b.as_f64().unwrap());
            if !close(a, b, scale.max(a.abs()), tol) {
                out.push(format!("{path}: {a:e} against {b:e}"));
            }
        }
        (Value::Array(a), Value::Array(b)) => {
            if a.len() != b.len() {
                out.push(format!("{path}: {} entries against {}", a.len(), b.len()));
            }
            for (i, (x, y)) in a.iter().zip(b).enumerate() {
                compare(&format!("{path}[{i}]"), x, y, scale, tol, out);
            }
        }
        (Value::Object(a), Value::Object(b)) => {
            let (ka, kb): (Vec<_>, Vec<_>) = (a.keys().collect(), b.keys().collect());
            if ka != kb {
                out.push(format!("{path}: names {ka:?} against {kb:?}"));
                return;
            }
            // the scale of a field's probe values is the field's own maximum
            let sc = a.get("max").and_then(|m| m.as_f64()).unwrap_or(scale);
            for k in a.keys() {
                compare(&format!("{path}.{k}"), &a[k], &b[k], sc, tol, out);
            }
        }
        _ => {
            if want != got {
                out.push(format!("{path}: {want} against {got}"));
            }
        }
    }
}

/// Checks (or, with NYX_WRITE_GOLDEN, writes) the record of `r` against the golden file `name`.
fn golden(name: &str, r: &SolveResult, tol: f64) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/golden_frames/{name}.json"));
    let got = record(r);
    if std::env::var("NYX_WRITE_GOLDEN").is_ok() {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, serde_json::to_string_pretty(&got).unwrap() + "\n").unwrap();
        return;
    }
    let want: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {path:?}"))).unwrap();
    // the text of an effective index and Q = Re/(2|Im|) carry the round-off of an imaginary part that is zero: the numbers are in `neff`
    let strip = |v: &Value| -> Value {
        let mut v = v.clone();
        if let Some(sc) = v["scalars"].as_object_mut() {
            sc.retain(|k, x| !(k.starts_with("Q ") || x.is_string() && k.starts_with("neff")));
        }
        v
    };
    let mut diffs = vec![];
    compare(name, &strip(&want), &strip(&got), 1.0, tol, &mut diffs);
    assert!(diffs.is_empty(), "{} differences from the golden file {name}:\n{}", diffs.len(), diffs.iter().take(12).cloned().collect::<Vec<_>>().join("\n"));
}

fn solve(mesh: &femcore::mesh::SimplexMesh, spec: &ProblemSpec) -> SolveResult {
    femcore::solve::solve(mesh, spec, &mut SolveState::default()).unwrap()
}

// ---------------------------------------------------------------------- scalar kernels

#[test]
fn golden_electrostatics_anisotropic_plates() {
    let mesh = common::rect(1.0, 0.5, 8, 4);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    // a generic orientation, so that every entry of the tensor (xz, xy, zy) is non-zero
    spec.materials.aniso.insert(1, Anisotropy { permittivity: Some(Tensor::Diagonal { x: 7.0, y: 3.0, z: 5.0 }), euler_deg: [30.0, 40.0, 10.0], ..Default::default() });
    spec.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![3], kind: BcKind::Ground }];
    golden("es_aniso_plates", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_electrostatics_series_with_charge() {
    let mesh = common::rect_split_x(1.0, 1.0, 0.5, 0.125);
    let mut spec = ProblemSpec::new(Physics::Electrostatics);
    spec.materials.set(1, "eps_r", 2.0);
    spec.materials.set(2, "eps_r", 8.0);
    spec.materials.set(2, "rho_v", 1e-9);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Ground }, Bc { on: None, labels: vec![3], kind: BcKind::SurfaceCharge(1e-10) }];
    golden("es_series_charge", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_anisotropic_stationary() {
    let mesh = common::rect(1.0, 0.4, 10, 4);
    let mut spec = ProblemSpec::new(Physics::Heat);
    spec.materials.aniso.insert(1, Anisotropy { conductivity: Some(Tensor::Diagonal { x: 30.0, y: 5.0, z: 2.0 }), euler_deg: [10.0, 30.0, 50.0], ..Default::default() });
    spec.materials.default.insert("Q".into(), 1e4);
    spec.bcs = vec![Bc { on: None, labels: vec![1], kind: BcKind::Temperature(300.0) }, Bc { on: None, labels: vec![4], kind: BcKind::Convection { h: 20.0, text: 280.0 } }, Bc { on: None, labels: vec![3], kind: BcKind::Flux(500.0) }];
    golden("heat_aniso_stationary", &solve(&mesh, &spec), 1e-9);
}

#[test]
fn golden_heat_transient_dense() {
    let mesh = common::rect(1.0, 0.2, 10, 2);
    let mut spec = ProblemSpec::new(Physics::Heat);
    spec.materials.default.insert("k".into(), 2.0);
    spec.materials.default.insert("rho".into(), 1.0);
    spec.materials.default.insert("Cp".into(), 1.0);
    spec.bcs = vec![Bc { on: None, labels: vec![1, 4], kind: BcKind::Temperature(0.0) }];
    spec.study = Study::Transient { t_end: 0.2, n_out: 4, t0: "sin(pi*x)".into() };
    golden("heat_transient", &solve(&mesh, &spec), 1e-9);
}

// ---------------------------------------------------------------------- optical modes

#[test]
fn golden_soi_wire_modes() {
    let wire = femcore::usecases::soi::SoiWire { n_modes: 2, ..Default::default() };
    golden("emw_soi_wire", &wire.solve().unwrap(), 1e-8);
}

#[test]
fn golden_soi_wire_with_sloped_wall_hybrid_modes() {
    let wire = femcore::usecases::soi::SoiWire { width: 300e-9, sidewall_left_deg: Some(75.0), n_modes: 2, ..Default::default() };
    golden("emw_soi_hybrid", &wire.solve().unwrap(), 1e-8);
}

fn crystal(tilt_deg: f64) -> Anisotropy {
    Anisotropy { permittivity: Some(Tensor::Diagonal { x: 4.0, y: 4.0, z: 6.5 }), euler_deg: [0.0, tilt_deg, 0.0], ..Default::default() }
}

#[test]
fn golden_tilted_crystal_in_a_pec_guide() {
    // the coupled (xz/yz) path of the complex solver, then the uncoupled real path
    for (name, tilt) in [("emw_pec_tilted", 40.0), ("emw_pec_straight", 0.0)] {
        let mesh = common::rect(0.03, 0.015, 12, 6);
        let mut s = ProblemSpec::new(Physics::Emw);
        s.order = 2;
        s.materials.aniso.insert(1, crystal(tilt));
        s.study = Study::Modes { freq: 12e9, n_target: 1.6, evanescent: false, n_modes: 3 };
        golden(name, &solve(&mesh, &s), 1e-7);
    }
}

#[test]
fn golden_pec_guide_with_a_pml() {
    let mesh = common::rect_split_x(0.03, 0.015, 0.02, 0.0025);
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.set(1, "eps_r", 4.0);
    s.materials.set(2, "eps_r", 1.0);
    // absorbing layer on the right (domain 2) along x
    s.materials.pml.insert(2, Pml { x: PmlDir::Plus, ..Pml::default() });
    s.study = Study::Modes { freq: 12e9, n_target: 1.5, evanescent: false, n_modes: 3 };
    let _ = C0;
    golden("emw_pml_x", &solve(&mesh, &s), 1e-7);
}
