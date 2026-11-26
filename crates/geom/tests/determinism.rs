// SPDX-License-Identifier: BSD-2-Clause
// Copyright (C) 2026 p-j-r
//! Subprocess two-tier determinism verification for `geom` (`RAYON_NUM_THREADS=1` vs `4`, 0-ULP):
//! verifies bit-identical (`to_bits()`) 2D arrangements, 3D solid evaluations, and projection oracle
//! queries across child processes with `RAYON_NUM_THREADS=1` and `RAYON_NUM_THREADS=4` under `--features parallel`.

use geom::eval::Tier;
use geom::math::V2;
use geom::model::Model;
use geom::sim::EntityRef;
use std::f64::consts::PI;

fn push_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_opt_u32(buf: &mut Vec<u8>, v: Option<u32>) {
    match v {
        Some(x) => {
            buf.push(1);
            push_u64(buf, x as u64);
        }
        None => {
            buf.push(0);
            push_u64(buf, 0);
        }
    }
}

fn collect_2d_arrangement_bits() -> Vec<u8> {
    let pitch = 18e-6;
    let r_c = 0.5 * pitch;
    let y_span = 160e-6;
    let n_turns = 80usize;
    let pts_per_turn = 16usize;
    let mut centerline: Vec<V2> = Vec::with_capacity(n_turns * pts_per_turn);
    for k in 0..n_turns {
        let cx = (k as f64 + 0.5) * pitch;
        let cy = if k % 2 == 0 { y_span } else { -y_span };
        for j in 0..pts_per_turn {
            let frac = j as f64 / (pts_per_turn - 1) as f64;
            let phi = if k % 2 == 0 { PI * (1.0 - frac) } else { -PI * (1.0 - frac) };
            centerline.push([cx + r_c * phi.cos(), cy + r_c * phi.sin()]);
        }
    }
    let expanded = geom::gds::expand_path_polygons(&centerline, 3.5e-6, 0, 0.0, 0.0);
    let pts: Vec<(f64, f64)> = expanded[0].iter().map(|p| (p[0], p[1])).collect();

    let mut m = Model::new_2d();
    let die = m.rect("die", (-60e-6, -240e-6), (1560e-6, 240e-6));
    let meander = m.polygon("meander", &pts);
    let pad_a = m.rounded_rect("pad_a", (40e-6, -220e-6), (220e-6, -180e-6), 8e-6);
    let pad_b = m.circle("pad_b", (800e-6, 200e-6), 16e-6);
    m.overlay(&[("die", die), ("meander", meander), ("pad_a", pad_a), ("pad_b", pad_b)]);

    let mut cache = geom::scene::Cache::new();
    let scene = geom::scene::build_scene(&m.list, &mut cache).unwrap();
    let cx = geom::arrange::arrange(&scene, &mut cache).unwrap();
    cx.check().unwrap();

    let mut out = Vec::new();
    push_u64(&mut out, cx.regions.len() as u64);
    for r in &cx.regions {
        push_u64(&mut out, r.measure.to_bits());
        for x in r.bbox.lo.iter().chain(r.bbox.hi.iter()) {
            push_u64(&mut out, x.to_bits());
        }
    }
    push_u64(&mut out, cx.edges.len() as u64);
    for e in &cx.edges {
        push_u64(&mut out, e.length.to_bits());
        push_u64(&mut out, e.t0.to_bits());
        push_u64(&mut out, e.t1.to_bits());
        push_u64(&mut out, e.v0 as u64);
        push_u64(&mut out, e.v1 as u64);
        push_opt_u32(&mut out, e.left);
        push_opt_u32(&mut out, e.right);
    }
    push_u64(&mut out, cx.vertices.len() as u64);
    for v in &cx.vertices {
        for c in v.p {
            push_u64(&mut out, c.to_bits());
        }
    }
    for i in 0..200 {
        let x = (i as f64 / 200.0) * 1500e-6;
        let z = ((i % 20) as f64 - 10.0) * 18e-6;
        push_opt_u32(&mut out, cx.locate([x, z]));
    }
    out
}

fn collect_3d_eval_and_oracle_bits() -> Vec<u8> {
    let mut m = Model::new_3d();
    let clad = m.box_("clad", (-4.0, -2.0, -3.0), (4.0, 2.0, 3.0));
    let hole = m.cylinder("hole", (0.0, 0.0, -3.0), 1.1, 6.0);
    let dome = m.sphere("dome", (2.4, 0.0, 0.0), 0.8);
    let shell = m.difference("shell", clad, &[hole]);
    m.overlay(&[("shell", shell), ("dome", dome)]);

    let ev = m.eval(Tier::Complex).unwrap();
    let cx = ev.complex();
    cx.check().unwrap();
    let oracle = ev.oracle();

    let mut out = Vec::new();
    push_u64(&mut out, ev.regions().len() as u64);
    for r in ev.regions() {
        push_u64(&mut out, r.measure.to_bits());
    }
    push_u64(&mut out, cx.vertices.len() as u64);
    for v in &cx.vertices {
        for c in v.p {
            push_u64(&mut out, c.to_bits());
        }
    }
    push_u64(&mut out, cx.faces.len() as u64);
    for (fi, f) in cx.faces.iter().enumerate() {
        push_u64(&mut out, f.area.to_bits());
        push_opt_u32(&mut out, f.first);
        push_opt_u32(&mut out, f.second);
        for c in f.bbox.lo.into_iter().chain(f.bbox.hi) {
            push_u64(&mut out, c.to_bits());
        }
        if let Some(q) = oracle.project([0.7, -0.5, 0.8], EntityRef::Face(fi as u32)) {
            out.push(1);
            for coord in q {
                push_u64(&mut out, coord.to_bits());
            }
        } else {
            out.push(0);
        }
    }
    out
}

#[test]
fn geom_2d_and_3d_eval_are_bit_exact_across_rayon_thread_counts() {
    if let Ok(out_path) = std::env::var("NYX_GEOM_DET_OUT") {
        let mut bytes = collect_2d_arrangement_bits();
        bytes.extend_from_slice(&collect_3d_eval_and_oracle_bits());
        std::fs::write(&out_path, &bytes).unwrap();
        return;
    }

    let exe = std::env::current_exe().unwrap();
    let run_child = |nt: &str| -> Vec<u8> {
        let out_path = std::env::temp_dir().join(format!("nyx-geom-det-{}-t{nt}.bin", std::process::id()));
        let out = std::process::Command::new(&exe)
            .args(["--exact", "geom_2d_and_3d_eval_are_bit_exact_across_rayon_thread_counts", "--nocapture"])
            .env("RAYON_NUM_THREADS", nt)
            .env("NYX_THREADS", nt)
            .env("NYX_GEOM_DET_OUT", &out_path)
            .output()
            .expect("spawn geom determinism subprocess");
        assert!(
            out.status.success(),
            "child RAYON_NUM_THREADS={nt} failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
        let bytes = std::fs::read(&out_path).unwrap();
        let _ = std::fs::remove_file(&out_path);
        bytes
    };

    let b1 = run_child("1");
    let b4 = run_child("4");
    assert!(!b1.is_empty());
    assert_eq!(b1, b4, "geom 2D/3D outputs must be 0-ULP bit-identical between RAYON_NUM_THREADS=1 and 4");
}
