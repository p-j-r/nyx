//! `fs-contour <bundle-dir> [options]`: all eigenvalues of `A x = λ B x` inside an ellipse, by the
//! contour-integral method; prints one JSON line (same fields as the SLEPc CISS mode of
//! `harness/c/slepc_eps.c`).
//!
//!   --center-re x --center-im y --a r --b s   ellipse semi-axes (b defaults to a)
//!   --nodes 24 --block 8 --moments 4 --workers 1 --reps 3 --rank-tol 1e-12

use femcore::linalg::bundle::{Bundle, read_meta};
use femcore::linalg::contour::{ContourOptions, Ellipse, contour_eigs};
use femcore::linalg::{C64, Csr};
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

fn run<T: femcore::linalg::faer_glue::FaerScalar>(dir: &str, kv: &HashMap<String, String>) {
    let bun = Bundle::<T>::read(Path::new(dir)).expect("read bundle");
    let f = |k: &str, d: f64| kv.get(k).map_or(d, |s| s.parse().expect(k));
    let u = |k: &str, d: usize| kv.get(k).map_or(d, |s| s.parse().expect(k));
    let a = f("a", 1.0);
    let region = Ellipse { center: C64::new(f("center-re", 0.0), f("center-im", 0.0)), a, b: f("b", a) };
    let opts = ContourOptions {
        nodes: u("nodes", 24),
        block: u("block", 8),
        moments: u("moments", 4),
        workers: u("workers", 1),
        rank_tol: f("rank-tol", 1e-12),
        ..ContourOptions::new(region)
    };
    let reps = u("reps", 3).max(1);
    let mut res = contour_eigs(&bun.a, bun.b.as_ref(), &opts).expect("contour");
    let mut best = f64::MAX;
    for _ in 0..reps {
        let t = Instant::now();
        res = contour_eigs(&bun.a, bun.b.as_ref(), &opts).expect("contour");
        best = best.min(t.elapsed().as_secs_f64());
    }
    let vals: Vec<String> = res.values.iter().map(|z| format!("[{:.17e},{:.17e}]", z.re, z.im)).collect();
    let resid: Vec<String> = res.residuals.iter().map(|x| format!("{x:.17e}")).collect();
    println!(
        "{{\"impl\":\"femsolve\",\"kind\":\"contour\",\"n\":{},\"nnz\":{},\"count\":{},\"estimated_count\":{:.6},\"count_stderr\":{:.6},\"rank\":{},\"saturated\":{},\"rejected\":{},\"total_min_s\":{best:.6e},\"t_factor\":{:.6e},\"t_solve\":{:.6e},\"t_svd\":{:.6e},\"t_project\":{:.6e},\"values\":[{}],\"true_rel_resid\":[{}]}}",
        bun.a.nrows(),
        bun.a.nnz(),
        res.values.len(),
        res.estimated_count,
        res.count_stderr,
        res.rank,
        res.saturated,
        res.rejected,
        res.timings.factor,
        res.timings.solve,
        res.timings.svd,
        res.timings.project,
        vals.join(","),
        resid.join(",")
    );
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    if raw.is_empty() {
        eprintln!("usage: fs-contour <bundle-dir> [options] (see the module docs of src/bin/fs_contour.rs)");
        std::process::exit(2);
    }
    let mut kv = HashMap::new();
    let mut i = 1;
    while i < raw.len() {
        kv.insert(raw[i].trim_start_matches("--").to_string(), raw.get(i + 1).cloned().unwrap_or_default());
        i += 2;
    }
    let meta = read_meta(Path::new(&raw[0])).expect("meta.txt");
    match meta.get("scalar").map(String::as_str) {
        Some("complex") => run::<C64>(&raw[0], &kv),
        _ => run::<f64>(&raw[0], &kv),
    }
    let _: Option<Csr<f64>> = None;
}
