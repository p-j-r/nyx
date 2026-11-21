// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Test kit for the physics crates: the golden-file harness (ADR 10). A golden file holds, for every solution of a case, the names of its fields
//! (flat names: a vector `E` is `|E|`, `Ex`, `Ez`), `n`, the sum of |v|, the max and a few probe values per field, the scalars and the
//! eigenvalue. A kernel that is ported must reproduce them. `NYX_WRITE_GOLDEN=1` writes every file, `NYX_WRITE_GOLDEN=<prefix>` the files whose
//! name starts with the prefix; do that only for a kernel known to be right, and say so in the commit.

use crate::result::{Scalar, SolveResult};
use serde_json::{json, Value};
use std::collections::BTreeMap;

/// The record of a result that goes into a golden file.
pub fn record(r: &SolveResult) -> Value {
    let mut sols = vec![];
    for s in &r.solutions {
        let mut fields = BTreeMap::new();
        for (name, _unit, v) in s.flat_fields() {
            let n = v.len();
            let l1: f64 = v.iter().map(|x| x.abs()).sum();
            let mx = v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            fields.insert(name, json!({"n": n, "l1": l1, "max": mx, "probe": [v[0], v[n / 3], v[n / 2], v[n - 1]]}));
        }
        sols.push(json!({"neff": s.eigenvalue.map(|e| json!([e[0], e[1]])), "value": s.value, "fields": fields}));
    }
    let mut scalars = BTreeMap::new();
    for (k, v) in &r.scalars {
        scalars.insert(
            k.clone(),
            match v {
                Scalar::Num(x) => json!(x),
                Scalar::Text(t) => json!(t),
            },
        );
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
pub fn golden(dir: &std::path::Path, name: &str, r: &SolveResult, tol: f64) {
    let path = dir.join(format!("{name}.json"));
    let got = record(r);
    // the kernels that have a fields stage hand out the cell fields
    if ["emw_", "es_", "heat_", "solid_"].iter().any(|p| name.starts_with(p)) {
        assert!(r.solutions.iter().all(|s| !s.cell_fields.is_empty()), "{name}: a solution without cell fields");
    }
    // NYX_WRITE_GOLDEN=1 writes every file, NYX_WRITE_GOLDEN=<prefix> only the files whose name starts with the prefix
    if std::env::var("NYX_WRITE_GOLDEN").is_ok_and(|w| w == "1" || name.starts_with(&w)) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(&path, serde_json::to_string_pretty(&got).unwrap() + "\n").unwrap();
        return;
    }
    let want: Value =
        serde_json::from_str(&std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {path:?}"))).unwrap();
    // the text of an effective index and Q = Re/(2|Im|) carry the round-off of an imaginary part that is zero: the numbers are in `neff`
    const ADDITIVE_FR3_SCALAR_PREFIXES: &[&str] = &[
        "n_g ",
        "v_g ",
        "n_B ",
        "Ex share ",
        "Ey share ",
        "Ez share ",
        "share in domain ",
        "layer fraction ",
        "f re ",
        "f im ",
    ];
    let want_scalars = want["scalars"].as_object().cloned().unwrap_or_default();
    let strip = |v: &Value, is_got: bool| -> Value {
        let mut v = v.clone();
        if let Some(sc) = v["scalars"].as_object_mut() {
            sc.retain(|k, x| {
                if k.starts_with("Q ") || (x.is_string() && k.starts_with("neff")) {
                    return false;
                }
                if is_got
                    && !want_scalars.contains_key(k)
                    && ADDITIVE_FR3_SCALAR_PREFIXES.iter().any(|p| k.starts_with(p))
                {
                    return false;
                }
                true
            });
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
    let (want, got) = (strip(&want, false), strip(&got, true));
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
    assert!(
        diffs.is_empty(),
        "{} differences from the golden file {name}:\n{}",
        diffs.len(),
        diffs.iter().take(12).cloned().collect::<Vec<_>>().join("\n")
    );
}

/// The fields stage hands out the same quantities as the plot fields before the average over the cells at a node: every cell field has the
/// name of a plot field (or its `Re ` version) and its node average is that plot field.
pub fn cell_fields_reproduce_the_plot_fields(mesh: &crate::Discretisation, r: &SolveResult) {
    for sol in &r.solutions {
        if sol.cell_fields.is_empty() {
            continue;
        }
        let flat = sol.flat_fields();
        let find = |n: &str| flat.iter().find(|(name, _, _)| name == n);
        for f in &sol.cell_fields {
            assert_eq!(f.nloc, mesh.nloc(r.order as usize));
            assert_eq!(f.re.len(), mesh.n_cells() * f.nloc, "{}", f.name);
            assert!(f.im.is_empty() || f.im.len() == f.re.len(), "{}", f.name);
            let plot = find(&f.name).or_else(|| find(&format!("Re {}", f.name))).unwrap_or_else(|| {
                panic!("cell field {} has no plot field in {:?}", f.name, flat.iter().map(|p| &p.0).collect::<Vec<_>>())
            });
            let avg = femcore::fields::node_average(mesh, &f.re, 1, r.order as usize);
            let scale = plot.2.iter().fold(0.0f64, |m, x| m.max(x.abs())).max(1e-300);
            let worst = avg.iter().zip(&plot.2).fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
            assert!(
                worst <= 1e-9 * scale,
                "node average of cell field {} differs from the plot field by {worst:e} (scale {scale:e})",
                f.name
            );
        }
    }
}
