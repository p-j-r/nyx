// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! ASME Layer 2 Pillar (f): Brand-new canonical JSON golden tests (ADR 10, ADR 26) for 3D curved-P2
//! electrostatics and `magneto::mixed_3d`, verified via parsed `serde_json::Value` comparisons without
//! altering any existing `tests/golden/*.json` files.

use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::SimplexMesh;
use nyx_dfs::{solve_electrostatics, solve_magnetostatics, Bc, BcKind, Problem, Study};
use nyx_physics::testkit::record;
use serde_json::Value;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::path::PathBuf;

fn golden_dir() -> PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn assert_json_value_close(path: &str, want: &Value, got: &Value, scale: f64, tol: f64, errs: &mut Vec<String>) {
    match (want, got) {
        (Value::Number(a), Value::Number(b)) => {
            let (fa, fb) = (a.as_f64().unwrap(), b.as_f64().unwrap());
            let sc = scale.max(fa.abs()).max(1e-300);
            if (fa - fb).abs() > tol * sc {
                errs.push(format!("{path}: wanted {fa:e}, got {fb:e}"));
            }
        }
        (Value::Array(a), Value::Array(b)) => {
            if a.len() != b.len() {
                errs.push(format!("{path}: array length {} != {}", a.len(), b.len()));
                return;
            }
            for (i, (x, y)) in a.iter().zip(b).enumerate() {
                assert_json_value_close(&format!("{path}[{i}]"), x, y, scale, tol, errs);
            }
        }
        (Value::Object(a), Value::Object(b)) => {
            let mut ka: Vec<_> = a.keys().collect();
            let mut kb: Vec<_> = b.keys().collect();
            ka.sort();
            kb.sort();
            if ka != kb {
                errs.push(format!("{path}: keys {ka:?} != {kb:?}"));
                return;
            }
            let sc = a.get("max").and_then(|v| v.as_f64()).unwrap_or(scale);
            for k in ka {
                assert_json_value_close(&format!("{path}.{k}"), &a[k], &b[k], sc, tol, errs);
            }
        }
        _ if want == got => {}
        _ => errs.push(format!("{path}: wanted {want}, got {got}")),
    }
}

fn verify_new_golden(name: &str, got: &Value) {
    let path = golden_dir().join(format!("{name}.json"));
    if std::env::var("NYX_WRITE_VVUQ_GOLDEN").is_ok() {
        let text = serde_json::to_string_pretty(got).unwrap() + "\n";
        std::fs::write(&path, text).unwrap();
    }
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("missing golden {}: {e}", path.display()));
    let want: Value = serde_json::from_str(&text).expect("golden JSON must parse cleanly");
    let mut errs = Vec::new();
    assert_json_value_close(name, &want, got, 1.0, 1e-9, &mut errs);
    assert!(errs.is_empty(), "golden mismatch for {name}:\n{}", errs.join("\n"));
}

fn build_curved_spherical_shell_golden(r_in: f64, r_out: f64, m: usize, n_r: usize) -> SimplexMesh {
    let mut grid_to_surf: BTreeMap<(usize, usize, usize), usize> = BTreeMap::new();
    let mut u_verts: Vec<[f64; 3]> = Vec::new();
    let mut get_surf = |ix: usize, iy: usize, iz: usize| -> usize {
        if let Some(&id) = grid_to_surf.get(&(ix, iy, iz)) {
            return id;
        }
        let id = u_verts.len();
        let q = [
            (0.25 * PI * (2.0 * (ix as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iy as f64) / (m as f64) - 1.0)).tan(),
            (0.25 * PI * (2.0 * (iz as f64) / (m as f64) - 1.0)).tan(),
        ];
        let len = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
        u_verts.push([q[0] / len, q[1] / len, q[2] / len]);
        grid_to_surf.insert((ix, iy, iz), id);
        id
    };
    let mut quads = Vec::new();
    for i in 0..m {
        for j in 0..m {
            quads.push([get_surf(m, i, j), get_surf(m, i + 1, j), get_surf(m, i + 1, j + 1), get_surf(m, i, j + 1)]);
            quads.push([get_surf(0, j, i), get_surf(0, j + 1, i), get_surf(0, j + 1, i + 1), get_surf(0, j, i + 1)]);
            quads.push([get_surf(j, m, i), get_surf(j + 1, m, i), get_surf(j + 1, m, i + 1), get_surf(j, m, i + 1)]);
            quads.push([get_surf(i, 0, j), get_surf(i + 1, 0, j), get_surf(i + 1, 0, j + 1), get_surf(i, 0, j + 1)]);
            quads.push([get_surf(i, j, m), get_surf(i + 1, j, m), get_surf(i + 1, j + 1, m), get_surf(i, j + 1, m)]);
            quads.push([get_surf(j, i, 0), get_surf(j + 1, i, 0), get_surf(j + 1, i + 1, 0), get_surf(j, i + 1, 0)]);
        }
    }
    let n_s = u_verts.len();
    let mut cm = CellMesh::<3>::new();
    for k in 0..=n_r {
        let r = r_in + (r_out - r_in) * (k as f64) / (n_r as f64);
        for u in &u_verts {
            cm.points.push([r * u[0], r * u[1], r * u[2]]);
        }
    }
    let vid = |k: usize, s: usize| -> u32 { (k * n_s + s) as u32 };
    for k in 0..n_r {
        for &[s0, s1, s2, s3] in &quads {
            cm.push_cell(
                CellType::Hex,
                &[
                    vid(k, s0),
                    vid(k, s1),
                    vid(k, s2),
                    vid(k, s3),
                    vid(k + 1, s0),
                    vid(k + 1, s1),
                    vid(k + 1, s2),
                    vid(k + 1, s3),
                ],
                1,
            );
        }
    }
    for &[s0, s1, s2, s3] in &quads {
        cm.push_facet(&[vid(0, s0), vid(0, s3), vid(0, s2), vid(0, s1)], 1);
        cm.push_facet(&[vid(n_r, s0), vid(n_r, s1), vid(n_r, s2), vid(n_r, s3)], 2);
    }
    let mut sm = SimplexMesh::from_cells(&cm).unwrap();
    let mut en = Vec::new();
    for &(va, vb) in sm.edges() {
        let pa = sm.point(va as usize);
        let pb = sm.point(vb as usize);
        let ra = (pa[0] * pa[0] + pa[1] * pa[1] + pa[2] * pa[2]).sqrt();
        let rb = (pb[0] * pb[0] + pb[1] * pb[1] + pb[2] * pb[2]).sqrt();
        if (ra - rb).abs() < 1e-11 {
            let mid: [f64; 3] = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
            let rmid = (mid[0] * mid[0] + mid[1] * mid[1] + mid[2] * mid[2]).sqrt();
            if rmid > 1e-15 {
                let sc = ra / rmid;
                en.push(([va, vb], [mid[0] * sc, mid[1] * sc, mid[2] * sc]));
            }
        }
    }
    sm.set_edge_nodes(&en);
    sm
}

#[test]
fn vvuq_golden_es3d_curved_spheres() {
    let mesh = build_curved_spherical_shell_golden(0.4, 1.0, 4, 2);
    let mut spec = Problem::new(Study::Electrostatic, 2);
    spec.materials.default.insert("eps_r".into(), 2.0);
    spec.bcs = vec![Bc { labels: vec![1], kind: BcKind::Potential(1.0) }, Bc { labels: vec![2], kind: BcKind::Ground }];
    let mut res = solve_electrostatics(&mesh, &spec).unwrap();
    res.scalars.retain(|(k, _)| k != nyx_physics::TLS_META_KEY);
    let val = record(&res);
    verify_new_golden("es3d_curved_coax", &val);
}

#[test]
fn vvuq_golden_mag3d_mixed_3d_block() {
    let mesh = femcore::testkit::block(1.0e-3, 1.0e-3, 1.0e-3, 3, 3, 3);
    let mut spec = Problem::new(Study::Magnetostatic { formulation: "mixed_3d".into() }, 1);
    spec.materials.default.insert("mu_r".into(), 1.0);
    spec.bcs = vec![
        Bc { labels: vec![1], kind: BcKind::LondonSheet { l_k_sq: 1.5e-12 } },
        Bc {
            labels: vec![1],
            kind: BcKind::CurrentTerminal { name: "sheet".into(), current: 1.0, direction: [0.0, 0.0, 1.0] },
        },
        Bc { labels: vec![2, 3, 4, 5, 6], kind: BcKind::Ground },
    ];
    let res = solve_magnetostatics(&mesh, &spec).unwrap();
    let val = record(&res);
    verify_new_golden("mag3d_mixed_loop", &val);
}
