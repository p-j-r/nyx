//! `fs-snes --problem bratu|atan [--nx 80] [--mu 6] [--n 10] [--x0 2] [--rtol 1e-8] [--ls bt|none] [--lag 1] [--maxit 50]`
//!
//! Newton on the same two problems as `harness/c/petsc_snes.c` (Bratu `−Δu − μ e^u` on an `nx × nx` grid with
//! spacing `1/(nx+1)`, and the diagonal arctan system) and one JSON line with the residual-norm history.

use femcore::solvers::Csr;
use femcore::solvers::csr::laplacian_2d;
use femcore::solvers::direct::Definiteness;
use femcore::solvers::snes::{LineSearch, NewtonOptions, newton};
use std::collections::HashMap;

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut kv: HashMap<String, String> = HashMap::new();
    let mut i = 0;
    while i + 1 < raw.len() + 1 && i < raw.len() {
        kv.insert(raw[i].trim_start_matches("--").to_string(), raw.get(i + 1).cloned().unwrap_or_default());
        i += 2;
    }
    let get = |k: &str, d: &str| kv.get(k).cloned().unwrap_or_else(|| d.to_string());
    let problem = get("problem", "bratu");
    let rtol: f64 = get("rtol", "1e-8").parse().unwrap();
    let lag: usize = get("lag", "1").parse().unwrap();
    let maxit: usize = get("maxit", "50").parse().unwrap();
    let ls = if get("ls", "bt") == "none" { LineSearch::None } else { LineSearch::default() };
    let opts = NewtonOptions { rtol, atol: 1e-50, stol: 1e-14, max_it: maxit, line_search: ls, jacobian_lag: lag, definiteness: Definiteness::General };
    let (n, res) = if problem == "atan" {
        let n: usize = get("n", "10").parse().unwrap();
        let x0: f64 = get("x0", "2").parse().unwrap();
        let mut f = |x: &[f64], o: &mut [f64]| x.iter().zip(o.iter_mut()).for_each(|(a, b)| *b = a.atan());
        let mut jac = |x: &[f64]| {
            let t: Vec<_> = x.iter().enumerate().map(|(i, a)| (i, i, 1.0 / (1.0 + a * a))).collect();
            Csr::from_triplets(n, n, &t).unwrap()
        };
        (n, newton(&mut f, &mut jac, &vec![x0; n], &opts).expect("newton"))
    } else {
        let nx: usize = get("nx", "80").parse().unwrap();
        let mu: f64 = get("mu", "6").parse().unwrap();
        let n = nx * nx;
        let h2inv = (nx as f64 + 1.0).powi(2);
        let l = laplacian_2d(nx, nx).map(|v| v * h2inv);
        let mut f = |u: &[f64], o: &mut [f64]| {
            l.spmv(u, o);
            for (oi, ui) in o.iter_mut().zip(u) {
                *oi -= mu * ui.exp();
            }
        };
        let mut jac = |u: &[f64]| {
            let d: Vec<_> = u.iter().enumerate().map(|(i, v)| (i, i, -mu * v.exp())).collect();
            Csr::axpby(1.0, &l, 1.0, &Csr::from_triplets(n, n, &d).unwrap())
        };
        (n, newton(&mut f, &mut jac, &vec![0.0; n], &opts).expect("newton"))
    };
    let norms: Vec<String> = res.residual_norms.iter().map(|v| format!("{v:.17e}")).collect();
    println!(
        "{{\"impl\":\"femsolve\",\"kind\":\"snes\",\"problem\":\"{problem}\",\"n\":{n},\"iterations\":{},\"reason\":\"{:?}\",\"function_evals\":{},\"jacobian_evals\":{},\"backtracks\":{},\"residual_norms\":[{}]}}",
        res.iterations,
        res.reason,
        res.function_evals,
        res.jacobian_evals,
        res.backtracks,
        norms.join(",")
    );
}
