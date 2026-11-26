// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME VVUQ Pillar (f): Parsed `serde_json::Value` golden test (ADR 26) for a 3D combined
//! LRB + Kelvin-Voigt viscoelastic eigenfrequency solve (`tests/golden/solid_3d_lrb_viscoelastic.json`).

mod common;

use femcore::material::{Anisotropy, Stiffness};
use nyx_physics::{testkit::cell_fields_reproduce_the_plot_fields, Scalar, SolveResult};
use nyx_sms::*;
use serde_json::{json, Value};
use std::path::PathBuf;

fn canonical_3d_lrb_viscoelastic_result() -> SolveResult {
    let mesh = common::block(1.0e-3, 0.3e-3, 0.25e-3, 4, 2, 2);
    let mut sp = Problem::new(Study::Eigenfrequency { n_modes: 3, f_search: 1.5e6 }, 2);
    sp.materials.aniso.insert(
        1,
        Anisotropy {
            stiffness: Some(Stiffness::Trigonal {
                c11: 496.8e9,
                c12: 163.6e9,
                c13: 110.9e9,
                c14: -23.5e9,
                c33: 498.1e9,
                c44: 147.4e9,
            }),
            euler_deg: [25.0, 40.0, 10.0],
            ..Default::default()
        },
    );
    sp.materials.default.insert("rho".into(), 3980.0);
    // Add proportional-scale Kelvin-Voigt shear/bulk viscosity alongside outer LRB radiation on right face (label 6).
    sp.materials.default.insert("eta".into(), 5.0e3);
    sp.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } },
        Bc { labels: vec![6], kind: BcKind::LowReflecting },
    ];
    let res = common::solve(&mesh, &sp);
    cell_fields_reproduce_the_plot_fields(&mesh, &res);
    res
}

fn parse_complex_hz(txt: &str) -> Option<[f64; 2]> {
    let trimmed = txt.trim().strip_suffix('i')?;
    let sep = trimmed.rfind(['+', '-'])?;
    if sep == 0 {
        return None;
    }
    let re: f64 = trimmed[..sep].trim().parse().ok()?;
    let im: f64 = trimmed[sep..].trim().parse().ok()?;
    Some([re, im])
}

fn summarize_value(r: &SolveResult) -> Value {
    let mut scalars = serde_json::Map::new();
    for (k, v) in &r.scalars {
        match v {
            Scalar::Num(q) => {
                scalars.insert(k.clone(), json!(q));
            }
            Scalar::Text(txt) => {
                if let Some([re, im]) = parse_complex_hz(txt) {
                    scalars.insert(k.clone(), json!([re, im]));
                }
            }
        }
    }
    let solutions: Vec<Value> = r
        .solutions
        .iter()
        .map(|s| {
            let mut fields = serde_json::Map::new();
            for name in ["|u|", "Re ux", "Re uy", "Re uz"] {
                let vals = s.flat(name).unwrap();
                let l1: f64 = vals.iter().map(|x| x.abs()).sum();
                let max_abs: f64 = vals.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
                let step = (vals.len() / 4).max(1);
                let probe: Vec<f64> = (0..4).map(|k| vals[(k * step).min(vals.len() - 1)]).collect();
                fields.insert(name.into(), json!({ "n": vals.len(), "l1": l1, "max": max_abs, "probe": probe }));
            }
            json!({
                "value": s.value,
                "fields": fields,
            })
        })
        .collect();
    json!({
        "physics": r.physics,
        "study": r.study,
        "dim": r.dim,
        "order": r.order,
        "dofs": r.dofs,
        "unstable_modes": r.scalar("unstable_modes").unwrap_or(0.0),
        "scalars": scalars,
        "solutions": solutions,
    })
}

fn field_tolerance(path: &str) -> f64 {
    if path.contains(".fields.Re u") {
        // Minor eigenvector components of far-from-shift overdamped modes vary up to ~6e-4 across glibc/UCRT libm.
        1e-3
    } else if path.starts_with("root.scalars.Q") || path.contains(".fields.|u|.") {
        2e-4
    } else if path.starts_with("root.scalars.f") || path.ends_with(".value") {
        1e-5
    } else {
        1e-7
    }
}

fn collect_json_diffs(got: &Value, exp: &Value, path: &str, diffs: &mut Vec<String>) {
    match (got, exp) {
        (Value::Number(ga), Value::Number(ea)) => {
            let g = ga.as_f64().unwrap();
            let e = ea.as_f64().unwrap();
            let scale = e.abs().max(1.0);
            let tol = field_tolerance(path);
            if (g - e).abs() > tol * scale {
                diffs.push(format!(
                    "{path}: got {g:.12e} vs expected {e:.12e} (rel={:.3e}, tol={tol:.1e})",
                    (g - e).abs() / scale
                ));
            }
        }
        (Value::String(gs), Value::String(es)) => {
            if gs != es {
                diffs.push(format!("string mismatch at {path}: got {gs} vs expected {es}"));
            }
        }
        (Value::Array(ga), Value::Array(ea)) => {
            if ga.len() != ea.len() {
                diffs.push(format!("array length mismatch at {path}: {} vs {}", ga.len(), ea.len()));
                return;
            }
            for (i, (g, e)) in ga.iter().zip(ea).enumerate() {
                collect_json_diffs(g, e, &format!("{path}[{i}]"), diffs);
            }
        }
        (Value::Object(go), Value::Object(eo)) => {
            if go.len() != eo.len() {
                diffs.push(format!("object keys mismatch at {path}: {:?} vs {:?}", go.keys(), eo.keys()));
                return;
            }
            for (k, ev) in eo {
                match go.get(k) {
                    Some(gv) => collect_json_diffs(gv, ev, &format!("{path}.{k}"), diffs),
                    None => diffs.push(format!("missing key {path}.{k}")),
                }
            }
        }
        _ => {
            if got != exp {
                diffs.push(format!("type/value mismatch at {path}: {got} vs {exp}"));
            }
        }
    }
}

#[test]
fn golden_solid_3d_lrb_viscoelastic_parsed_json_parity() {
    let res = canonical_3d_lrb_viscoelastic_result();
    assert_eq!(res.solutions.len(), 3);
    assert!(res.scalar("Q1").is_some_and(|q| q.is_finite() && q > 0.0));
    let got = summarize_value(&res);

    let path: PathBuf = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/solid_3d_lrb_viscoelastic.json");
    if std::env::var_os("NYX_WRITE_VVUQ_GOLDEN").is_some() {
        std::fs::write(&path, format!("{}\n", serde_json::to_string_pretty(&got).unwrap())).unwrap();
        return;
    }
    let text = std::fs::read_to_string(&path).expect("read tests/golden/solid_3d_lrb_viscoelastic.json");
    let exp: Value = serde_json::from_str(&text).expect("parse solid_3d_lrb_viscoelastic.json");
    let mut diffs = Vec::new();
    collect_json_diffs(&got, &exp, "root", &mut diffs);
    assert!(diffs.is_empty(), "golden mismatches:\n{}", diffs.join("\n"));
}
