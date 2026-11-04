//! `fs-ksp <bundle-dir> [options]`: solves `A x = b` (`b` = the bundle's `v0.mtx`) with a Krylov
//! method and prints one JSON line in the same schema as the PETSc driver
//! (`harness/c/petsc_ksp.c`).
//!
//!   --solver cg|cocg|gmres|fgmres|minres|bcgs [cocg]   --pc none|jacobi|amg|ilu|icc [none]   (amg: --hermitian 1 for R = Pᴴ; ilu/icc: --levels k [0]; asm: --blocks p --overlap k --asm-type restrict|basic)   --rtol 1e-5   --maxit 10000
//!   --restart 30   --side left|right [left]   --orth cgs|mgs|cgs-refine [cgs]   --norm pc|unpc|natural [pc]
//!   --reps 3   --petsc-sign-test 1 (reproduce PETSc's COCG breakdown test, for parity runs)

use femcore::linalg::amg::{Amg, AmgOptions};
use femcore::linalg::bundle::{Bundle, read_meta};
use femcore::linalg::ilu::{Ic, Ilu};
use femcore::linalg::ksp::{
    Bicgstab, BicgstabOptions, Cg, CgKind, CgNorm, CgOptions, Gmres, GmresOptions, Minres, MinresOptions,
    Orthogonalization, PcSide, SolveInfo, Tolerances,
};
use femcore::linalg::linop::{Identity, Preconditioner};
use femcore::linalg::pc::Jacobi;
use femcore::linalg::schwarz::{Schwarz, SchwarzKind, SchwarzOptions, contiguous_partition};
use femcore::linalg::{C64, Csr};
use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

fn run<T: femcore::linalg::faer_glue::FaerScalar>(dir: &str, kv: &HashMap<String, String>) {
    let bun = Bundle::<T>::read(Path::new(dir)).expect("read bundle");
    let a = &bun.a;
    let b = bun.v0.clone().expect("bundle needs v0.mtx as the right-hand side");
    let get = |k: &str, d: &str| kv.get(k).cloned().unwrap_or_else(|| d.to_string());
    let tol = Tolerances {
        rtol: get("rtol", "1e-5").parse().unwrap(),
        max_it: get("maxit", "10000").parse().unwrap(),
        ..Default::default()
    };
    let reps: usize = get("reps", "3").parse().unwrap();
    let solver = get("solver", "cocg");
    let jac = Jacobi::new(a);
    let id = Identity;
    let t_setup = Instant::now();
    let amg = if get("pc", "none") == "amg" {
        Some(
            Amg::new(a, AmgOptions { hermitian: kv.contains_key("hermitian"), ..Default::default() })
                .expect("amg setup"),
        )
    } else {
        None
    };
    let levels: u32 = get("levels", "0").parse().unwrap();
    let ilu = (get("pc", "none") == "ilu").then(|| Ilu::ilu_k(a, levels));
    let icc = (get("pc", "none") == "icc").then(|| Ic::ic_k(a, levels, true));
    let asm = (get("pc", "none") == "asm").then(|| {
        let part = contiguous_partition(a.nrows(), get("blocks", "4").parse().unwrap());
        let kind = if get("asm-type", "restrict") == "basic" { SchwarzKind::Basic } else { SchwarzKind::Restricted };
        Schwarz::new(
            a,
            &part,
            &SchwarzOptions { overlap: get("overlap", "1").parse().unwrap(), kind, ..Default::default() },
        )
        .expect("schwarz setup")
    });
    let setup_s = t_setup.elapsed().as_secs_f64();
    let pc: &dyn Preconditioner<T> = match (get("pc", "none").as_str(), &amg) {
        ("jacobi", _) => &jac,
        ("amg", Some(m)) => m,
        ("ilu", _) => ilu.as_ref().unwrap(),
        ("icc", _) => icc.as_ref().unwrap(),
        ("asm", _) => asm.as_ref().unwrap(),
        _ => &id,
    };
    let norm = match get("norm", "pc").as_str() {
        "unpc" => CgNorm::Unpreconditioned,
        "natural" => CgNorm::Natural,
        _ => CgNorm::Preconditioned,
    };
    let orthog = match get("orth", "cgs").as_str() {
        "mgs" => Orthogonalization::ModifiedGramSchmidt,
        "cgs-refine" => Orthogonalization::ClassicalGramSchmidtIfNeeded,
        _ => Orthogonalization::ClassicalGramSchmidt,
    };
    let side = if get("side", "left") == "right" { PcSide::Right } else { PcSide::Left };
    let restart: usize = get("restart", "30").parse().unwrap();
    let n = a.nrows();
    let solve = |x: &mut [T]| -> SolveInfo {
        match solver.as_str() {
            "cg" | "cocg" => {
                let kind = if solver == "cg" { CgKind::Hermitian } else { CgKind::Symmetric };
                Cg::new(CgOptions {
                    tol,
                    kind,
                    norm,
                    record_history: true,
                    petsc_sign_test: kv.contains_key("petsc-sign-test"),
                    ..Default::default()
                })
                .solve(a, pc, &b, x)
            }
            "gmres" | "fgmres" => Gmres::new(GmresOptions {
                tol,
                restart,
                side,
                orthog,
                flexible: solver == "fgmres",
                record_history: true,
                ..Default::default()
            })
            .solve(a, pc, &b, x),
            "minres" => Minres::new(MinresOptions { tol, record_history: true }).solve(a, pc, &b, x),
            "bcgs" => Bicgstab::new(BicgstabOptions { tol, record_history: true }).solve(a, pc, &b, x),
            other => panic!("unknown solver {other}"),
        }
    };
    let mut x = vec![T::ZERO; n];
    let mut info = solve(&mut x); // warm-up
    let mut best = f64::MAX;
    for _ in 0..reps {
        let t = Instant::now();
        info = solve(&mut x);
        best = best.min(t.elapsed().as_secs_f64());
    }
    let mut ax = vec![T::ZERO; n];
    a.spmv(&x, &mut ax);
    let rr: f64 = ax.iter().zip(&b).map(|(p, q)| (*q - *p).modulus_sq()).sum::<f64>().sqrt();
    let bn: f64 = b.iter().map(|v| v.modulus_sq()).sum::<f64>().sqrt();
    let hist: Vec<String> = info.history.iter().map(|v| format!("{v:.17e}")).collect();
    println!(
        "{{\"impl\":\"femsolve\",\"kind\":\"ksp\",\"n\":{n},\"nnz\":{},\"scalar\":\"{}\",\"solver\":\"{solver}\",\"iterations\":{},\"reason\":\"{:?}\",\"setup_s\":{setup_s:.6e},\"solve_min_s\":{best:.6e},\"true_rel_resid\":{:.17e},\"history\":[{}]}}",
        a.nnz(),
        if T::IS_COMPLEX { "complex" } else { "real" },
        info.iterations,
        info.reason,
        rr / bn,
        hist.join(",")
    );
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    if raw.is_empty() {
        eprintln!("usage: fs-ksp <bundle-dir> [options] (see the module docs of src/bin/fs_ksp.rs)");
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
