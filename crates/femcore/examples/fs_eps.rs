// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! `fs-eps <bundle-dir> [options]`: solves a bundle with Krylov–Schur and prints one JSON line.
//!
//! The JSON has the same fields as the SLEPc driver (`harness/c/slepc_eps.c`), so the two can
//! be compared by `harness/run.py`.
//!
//! Options (defaults in brackets):
//!   --mode shift|direct   [shift if the bundle suggests a shift, else direct]
//!   --sigma-re x --sigma-im y   shift [bundle's sigma_re / sigma_im]
//!   --nev k [bundle's nev]   --ncv c [SLEPc default]   --tol t [1e-8]   --maxit m
//!   --which largest|smallest|largest_real|smallest_real   (direct mode only) [largest]
//!   --block b [1]   --keep f [0.5]   --reps r [3]   --definite (try Cholesky first)
//!   --conv theta|original   convergence test: SLEPc's ‖r‖/|θ| [theta, for parity] or the relative
//!                           residual of the original problem [original]
//!   --factor-threads global|seq|<n> [global]   --solve-threads global|seq|<n> [seq]
//!   --history (log the error estimates of every restart)

use femcore::linalg::C64;
use femcore::linalg::bundle::{Bundle, read_meta};
use femcore::linalg::direct::{Definiteness, ThreadPolicy, Threads};
use femcore::linalg::eps::{EpsOptions, EpsReason, Which, krylov_schur};
use femcore::linalg::faer_glue::FaerScalar;
use femcore::linalg::problem::{EigenProblem, KrylovSchurShiftInvert, ShiftInvertEigenSolver, Symmetry};
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

struct Args {
    dir: String,
    kv: HashMap<String, String>,
    flags: Vec<String>,
}

impl Args {
    fn parse() -> Self {
        let raw: Vec<String> = std::env::args().skip(1).collect();
        if raw.is_empty() {
            eprintln!("usage: fs-eps <bundle-dir> [options]   (see the module docs of src/bin/fs_eps.rs)");
            std::process::exit(2);
        }
        let (mut kv, mut flags) = (HashMap::new(), Vec::new());
        let mut i = 1;
        while i < raw.len() {
            let k = raw[i].trim_start_matches("--").to_string();
            if matches!(k.as_str(), "definite" | "history") {
                flags.push(k);
                i += 1;
            } else {
                kv.insert(k, raw.get(i + 1).cloned().unwrap_or_default());
                i += 2;
            }
        }
        Self { dir: raw[0].clone(), kv, flags }
    }
    fn f64(&self, k: &str) -> Option<f64> {
        self.kv.get(k).and_then(|s| s.parse().ok())
    }
    fn usize(&self, k: &str) -> Option<usize> {
        self.kv.get(k).and_then(|s| s.parse().ok())
    }
    fn flag(&self, k: &str) -> bool {
        self.flags.iter().any(|f| f == k)
    }
}

fn json_values(v: &[C64]) -> String {
    v.iter().map(|z| format!("[{:.17e},{:.17e}]", z.re, z.im)).collect::<Vec<_>>().join(",")
}

fn history_json(h: &[femcore::linalg::eps::HistoryEntry]) -> String {
    if h.is_empty() {
        return String::new();
    }
    let items: Vec<String> = h
        .iter()
        .map(|e| {
            format!(
                "{{\"cycle\":{},\"nconv\":{},\"nv\":{},\"matvecs\":{},\"errors\":[{}],\"ritz\":[{}]}}",
                e.cycle,
                e.nconv,
                e.nv,
                e.matvecs,
                json_list(&e.errors),
                json_values(&e.ritz)
            )
        })
        .collect();
    format!(",\"history\":[{}]", items.join(","))
}

fn json_list(v: &[f64]) -> String {
    v.iter().map(|x| format!("{x:.17e}")).collect::<Vec<_>>().join(",")
}

fn run<T: FaerScalar>(args: &Args) {
    let bun = Bundle::<T>::read(Path::new(&args.dir)).expect("read bundle");
    let n = bun.a.nrows();
    let nev = args.usize("nev").or_else(|| bun.get_f64("nev").map(|x| x as usize)).unwrap_or(4);
    let (sre, sim) = (
        args.f64("sigma-re").or(bun.get_f64("sigma_re")).unwrap_or(0.0),
        args.f64("sigma-im").or(bun.get_f64("sigma_im")).unwrap_or(0.0),
    );
    let mode = args.kv.get("mode").cloned().unwrap_or_else(|| {
        if bun.b.is_some() || bun.meta.contains_key("sigma_re") { "shift".into() } else { "direct".into() }
    });
    let reps = args.usize("reps").unwrap_or(3).max(1);
    let hermitian = bun.meta.get("symmetry").map(|s| s == "hermitian").unwrap_or(false);
    let symmetry = if hermitian { Symmetry::Hermitian } else { Symmetry::General };
    let mut opts = if hermitian { EpsOptions::hermitian(nev) } else { EpsOptions::general(nev) };
    opts.ncv = args.usize("ncv");
    opts.tol = args.f64("tol").unwrap_or(1e-8);
    opts.max_it = args.usize("maxit");
    opts.block_size = args.usize("block").unwrap_or(1);
    opts.keep = args.f64("keep").unwrap_or(0.5);
    opts.record_history = args.flag("history");
    let nnz = bun.a.nnz();
    let scalar = if T::IS_COMPLEX { "complex" } else { "real" };

    let (values, resid, restarts, matvecs, nconv, reason, total_min, eps_min, factor_s, extra, ncv_used);
    let mut cold_s = 0.0;
    if mode == "shift" {
        let problem = EigenProblem { a: bun.a.clone(), b: bun.b.clone(), symmetry };
        let def = if args.flag("definite") { Definiteness::PositiveDefinite } else { Definiteness::General };
        let mut s = KrylovSchurShiftInvert::new(&problem, def).expect("setup");
        s.opts = EpsOptions { which: Which::LargestMagnitude, ..opts.clone() };
        s.start = bun.v0.clone();
        let parse_threads = |k: &str, default: Threads| match args.kv.get(k).map(String::as_str) {
            None => default,
            Some("global") => Threads::Global,
            Some("seq") => Threads::Sequential,
            Some(n) => Threads::Fixed(n.parse().expect("--factor-threads / --solve-threads: global|seq|<n>")),
        };
        let dflt = ThreadPolicy::default();
        s.set_thread_policy(ThreadPolicy {
            factor: parse_threads("factor-threads", dflt.factor),
            solve: parse_threads("solve-threads", dflt.solve),
        });
        s.theta_convergence = args.kv.get("conv").map(String::as_str) != Some("original");
        let sigma = T::from_parts(sre, sim);
        let t_cold = Instant::now();
        let first = s.solve_near(sigma, nev).expect("solve"); // warm-up (includes the symbolic analysis)
        cold_s = t_cold.elapsed().as_secs_f64() + s.factor_stats().analyse_secs * 0.0;
        let (mut best_total, mut best_eps, mut best_fac) = (f64::MAX, f64::MAX, f64::MAX);
        let mut last = first;
        for _ in 0..reps {
            let t0 = Instant::now();
            last = s.solve_near(sigma, nev).expect("solve");
            let tt = t0.elapsed().as_secs_f64();
            if tt < best_total {
                best_total = tt;
                best_eps = last.eps_timings.total;
                best_fac = last.factor_stats.numeric_secs_last;
            }
        }
        values = last.values.iter().take(nev).copied().collect::<Vec<_>>();
        resid = last.residuals.iter().take(nev).copied().collect::<Vec<_>>();
        restarts = last.restarts;
        ncv_used = last.ncv;
        matvecs = last.matvecs;
        nconv = last.values.len();
        reason = if last.converged { EpsReason::Converged } else { EpsReason::DivergedIts };
        total_min = best_total;
        eps_min = best_eps;
        factor_s = best_fac;
        let st = last.factor_stats;
        let hist = history_json(&last.history);
        extra = format!(
            "{hist},\"factor\":\"{:?}\",\"symbolic_s\":{:.6e},\"symbolic_count\":{},\"breakdowns\":{},\"orth_calls\":{},\"orth_refinements\":{},\"t_operator\":{:.6e},\"t_weight\":{:.6e},\"t_orth\":{:.6e},\"t_dense\":{:.6e},\"t_restart\":{:.6e}",
            st.used,
            st.analyse_secs,
            st.symbolic_count,
            last.breakdowns,
            last.orth_calls,
            last.orth_refinements,
            last.eps_timings.operator,
            last.eps_timings.weight,
            last.eps_timings.orthogonalize,
            last.eps_timings.dense,
            last.eps_timings.restart
        );
    } else {
        let which = match args.kv.get("which").map(String::as_str).unwrap_or("largest") {
            "largest" => Which::LargestMagnitude,
            "smallest" => Which::SmallestMagnitude,
            "largest_real" => Which::LargestReal,
            "smallest_real" => Which::SmallestReal,
            other => panic!("unknown --which {other}"),
        };
        opts.which = which;
        assert!(bun.b.is_none(), "direct mode needs a standard problem (no B)");
        let _ = krylov_schur(&bun.a, None, bun.v0.as_deref(), &opts).expect("solve");
        let (mut best, mut last) = (f64::MAX, None);
        for _ in 0..reps {
            let t0 = Instant::now();
            let r = krylov_schur(&bun.a, None, bun.v0.as_deref(), &opts).expect("solve");
            best = best.min(t0.elapsed().as_secs_f64());
            last = Some(r);
        }
        let r = last.unwrap();
        values = r.values.iter().take(nev).copied().collect::<Vec<_>>();
        let mut ax = vec![T::ZERO; n];
        resid = (0..values.len())
            .map(|i| {
                bun.a.spmv(r.vector(i), &mut ax);
                let s: f64 =
                    ax.iter().zip(r.vector(i)).map(|(a, x)| (a.to_c64() - values[i] * x.to_c64()).norm_sqr()).sum();
                s.sqrt() / values[i].norm()
            })
            .collect();
        restarts = r.iterations;
        ncv_used = r.ncv;
        matvecs = r.matvecs;
        nconv = r.nconv;
        reason = r.reason;
        total_min = best;
        eps_min = best;
        factor_s = 0.0;
        extra = format!(
            ",\"breakdowns\":{},\"t_operator\":{:.6e},\"t_orth\":{:.6e},\"t_dense\":{:.6e},\"t_restart\":{:.6e},\"orth_refinements\":{},\"orth_calls\":{}",
            r.breakdowns,
            r.timings.operator,
            r.timings.orthogonalize,
            r.timings.dense,
            r.timings.restart,
            r.orth_refinements,
            r.orth_calls
        );
    }
    println!(
        "{{\"impl\":\"femsolve\",\"kind\":\"eps\",\"n\":{n},\"nnz\":{nnz},\"scalar\":\"{scalar}\",\"mode\":\"{mode}\",\"nev\":{nev},\"ncv\":{ncv_used},\"block\":{},\"tol\":{:e},\"restarts\":{restarts},\"matvecs\":{matvecs},\"nconv\":{nconv},\"reason\":\"{:?}\",\"cold_s\":{cold_s:.6e},\"total_min_s\":{total_min:.6e},\"eps_min_s\":{eps_min:.6e},\"factor_s\":{factor_s:.6e},\"values\":[{}],\"true_rel_resid\":[{}]{extra}}}",
        opts.block_size,
        opts.tol,
        reason,
        json_values(&values),
        json_list(&resid)
    );
}

fn main() {
    let args = Args::parse();
    let meta = read_meta(Path::new(&args.dir)).expect("read meta.txt");
    match meta.get("scalar").map(String::as_str) {
        Some("complex") => run::<C64>(&args),
        _ => run::<f64>(&args),
    }
}
