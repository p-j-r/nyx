//! M2 study: matrix-free sum factorisation against the assembled CSR matrix for the 3-D spectral-element
//! Helmholtz operator on a warped mesh (on an axis-aligned box the assembled GLL matrices are very sparse, ~3(p+1) per
//! row, and assembly wins: see docs/LEARNINGS.md) at (almost) fixed unknown count and growing polynomial degree.
//! `cargo run --release --example matfree_bench [axis_nodes_target]`  (build with `--features parallel` for threads)
use femcore::solvers::linop::LinearOperator;
use femcore::solvers::matfree::{Boundary, SpectralBox};
use femcore::solvers::rng::random_vector;
use std::time::Instant;

fn time_apply(op: &dyn LinearOperator<f64>, x: &[f64], y: &mut [f64]) -> f64 {
    op.apply(x, y); // warm up
    let reps = 5;
    let t = Instant::now();
    for _ in 0..reps {
        op.apply(x, y);
    }
    t.elapsed().as_secs_f64() / reps as f64
}

fn main() {
    let axis = std::env::args().nth(1).and_then(|s| s.parse::<usize>().ok()).unwrap_or(24);
    println!(
        "3-D Helmholtz-type SEM operator on a WARPED (non-affine) mesh, ~{}^3 nodes; assembled element blocks are dense",
        axis + 1
    );
    println!(
        "{:>3} {:>6} {:>9} {:>12} {:>12} {:>9} {:>14} {:>12}",
        "p", "elems", "dofs", "matfree [ms]", "csr [ms]", "speed-up", "csr bytes/dof", "assemble [s]"
    );
    for p in [1usize, 2, 3, 4, 6, 8, 12] {
        if axis % p != 0 {
            continue;
        }
        let nel = axis / p;
        let h = 1.0 / nel as f64;
        let base = SpectralBox::<f64>::uniform(p, &[nel, nel, nel], &[h, h, h], 1.0, 0.0, Boundary::Dirichlet);
        let nn = base.nodes_per_axis()[0];
        let coords: Vec<Vec<f64>> = (0..3)
            .map(|d| {
                (0..base.len())
                    .map(|g| {
                        let x = [base.coord(0, g % nn), base.coord(1, (g / nn) % nn), base.coord(2, g / (nn * nn))];
                        let bump: f64 = x.iter().map(|v| (std::f64::consts::PI * v).sin()).product();
                        x[d] + 0.08 * (1.0 + 0.5 * d as f64) * bump
                    })
                    .collect()
            })
            .collect();
        let ne = nel.pow(3);
        let op = SpectralBox::<f64>::mapped(
            p,
            &[nel, nel, nel],
            &coords,
            vec![1.0; ne],
            vec![-20.0; ne],
            Boundary::Dirichlet,
        );
        let n = op.len();
        let x = random_vector::<f64>(n, 1);
        let mut y = vec![0.0; n];
        let t_free = time_apply(&op, &x, &mut y);
        let (csr_ms, bytes, asm) = if p <= 6 {
            let t = Instant::now();
            let a = op.assemble();
            let asm = t.elapsed().as_secs_f64();
            let mut y2 = vec![0.0; n];
            let ta = time_apply(&a, &x, &mut y2);
            let err = y.iter().zip(&y2).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
                / y2.iter().map(|v| v * v).sum::<f64>().sqrt();
            assert!(err < 1e-12, "mismatch {err}");
            (
                format!("{:.2}", ta * 1e3),
                format!("{:.0}", (a.nnz() * 12 + (n + 1) * 8) as f64 / n as f64),
                format!("{asm:.2}"),
            )
        } else {
            ("-".into(), "-".into(), "-".into())
        };
        let speed = if p <= 6 {
            let a = csr_ms.parse::<f64>().unwrap() / (t_free * 1e3);
            format!("{a:.1}x")
        } else {
            "-".into()
        };
        println!(
            "{p:>3} {:>6} {n:>9} {:>12.2} {csr_ms:>12} {speed:>9} {bytes:>14} {asm:>12}",
            nel.pow(3),
            t_free * 1e3
        );
    }
}
