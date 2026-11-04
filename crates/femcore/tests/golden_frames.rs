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
            // the order of the keys of an object is not part of a golden: `serde_json` keeps insertion order when any crate of the build
            // enables `preserve_order` (the app does), and sorts them otherwise, so `cargo test --workspace` and `-p femcore` differ
            let (mut ka, mut kb): (Vec<_>, Vec<_>) = (a.keys().collect(), b.keys().collect());
            ka.sort();
            kb.sort();
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
    // the kernels that have a fields stage hand out the cell fields
    if ["emw_", "es_", "heat_", "solid_"].iter().any(|p| name.starts_with(p)) {
        assert!(r.solutions.iter().all(|s| s.cell_fields.is_some()), "{name}: a solution without cell fields");
    }
    // NYX_WRITE_GOLDEN=1 writes every file, NYX_WRITE_GOLDEN=<prefix> only the files whose name starts with the prefix
    if std::env::var("NYX_WRITE_GOLDEN").is_ok_and(|w| w == "1" || name.starts_with(&w)) {
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
            // a Stokes parameter that is zero prints with either sign
            for x in sc.values_mut() {
                if let Some(t) = x.as_str() {
                    let t = t.replace("-0.000", "+0.000");
                    *x = Value::String(t);
                }
            }
        }
        v
    };
    let (want, got) = (strip(&want), strip(&got));
    let mut diffs = vec![];
    compare(name, &want, &got, 1.0, tol, &mut diffs);
    if !diffs.is_empty() {
        // an eigenvector is fixed only up to its sign when two of its entries tie in modulus: accept a solution whose signed fields
        // are all negated (the other fields, scalars and the rest of the record must still agree)
        let mut gauged = got.clone();
        if let (Some(ws), Some(gs)) = (want["solutions"].as_array(), gauged["solutions"].as_array_mut()) {
            for (w, g) in ws.iter().zip(gs.iter_mut()) {
                let mut d = vec![];
                compare("", w, g, 1.0, tol, &mut d);
                let mut flipped = g.clone();
                if !d.is_empty() {
                    if let Some(fs) = flipped["fields"].as_object_mut() {
                        for (fname, f) in fs.iter_mut() {
                            if !fname.starts_with('|') && fname != "von Mises" {
                                if let Some(pr) = f["probe"].as_array_mut() {
                                    pr.iter_mut().for_each(|x| *x = json!(-x.as_f64().unwrap()));
                                }
                            }
                        }
                    }
                    let mut d2 = vec![];
                    compare("", w, &flipped, 1.0, tol, &mut d2);
                    if d2.is_empty() {
                        *g = flipped;
                    }
                }
            }
        }
        diffs.clear();
        compare(name, &want, &gauged, 1.0, tol, &mut diffs);
    }
    assert!(diffs.is_empty(), "{} differences from the golden file {name}:\n{}", diffs.len(), diffs.iter().take(12).cloned().collect::<Vec<_>>().join("\n"));
}

fn solve(mesh: &femcore::mesh::SimplexMesh, spec: &ProblemSpec) -> SolveResult {
    let r = femcore::solve::solve(mesh, spec, &mut SolveState::default()).unwrap();
    cell_fields_reproduce_the_plot_fields(mesh, &r);
    r
}

/// The fields stage hands out the same quantities as the plot fields before the average over the cells at a node: every
/// cell field has the name of a plot field (or its `Re ` version) and its node average is that plot field.
fn cell_fields_reproduce_the_plot_fields(mesh: &femcore::mesh::SimplexMesh, r: &SolveResult) {
    for sol in &r.solutions {
        let Some(cf) = &sol.cell_fields else { continue };
        assert_eq!(cf.nloc, mesh.nloc(r.order));
        assert!(!cf.fields.is_empty());
        for f in &cf.fields {
            assert_eq!(f.re.len(), mesh.n_cells() * cf.nloc, "{}", f.name);
            assert!(f.im.is_empty() || f.im.len() == f.re.len(), "{}", f.name);
            let plot = sol
                .field(&f.name)
                .or_else(|| sol.field(&format!("Re {}", f.name)))
                .unwrap_or_else(|| panic!("cell field {} has no plot field in {:?}", f.name, sol.fields.iter().map(|p| &p.name).collect::<Vec<_>>()));
            let avg = femcore::fields::node_average(mesh, &f.re, 1, cf.order);
            let scale = plot.values.iter().fold(0.0f64, |m, x| m.max(x.abs())).max(1e-300);
            let worst = avg.iter().zip(&plot.values).fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(worst <= 1e-9 * scale, "node average of cell field {} differs from the plot field by {worst:e} (scale {scale:e})", f.name);
        }
    }
}

// ---------------------------------------------------------------------- scalar kernels

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

