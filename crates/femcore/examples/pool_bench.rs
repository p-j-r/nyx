//! Cost of the parallel kernels of the solvers: sparse matrix-vector product at several sizes and the time of one
//! multigrid-preconditioned CG iteration. Used to compare thread strategies.
//!
//!   cargo run --release -p femcore --example pool_bench            (NYX_THREADS=N limits the pool)
use femcore::amg::{pcg, Amg, AmgOptions};
use femcore::fem::assemble_laplace;
use femcore::mesh::SimplexMesh;
use femcore::sparse::{Csr, Csr32};
use std::time::Instant;

/// median and minimum of `reps` timings in microseconds
fn bench<F: FnMut()>(mut f: F, reps: usize) -> (f64, f64) {
    for _ in 0..3 {
        f();
    }
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let s = Instant::now();
            f();
            s.elapsed().as_secs_f64() * 1e6
        })
        .collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (t[t.len() / 2], t[0])
}

fn main() {
    println!("== Csr32::matvec_into on P1 Laplace matrices of a block (median / min microseconds) ==");
    for n in [14usize, 24, 32, 44, 60, 80, 110] {
        let m = SimplexMesh::block(1.0, 1.0, 1.0, n, n, n);
        let a: Csr = assemble_laplace(&m, 1, &vec![1.0; m.n_cells()]);
        let a32 = Csr32::from_csr(&a);
        let x: Vec<f64> = (0..a.ncols).map(|i| ((i * 7919) % 101) as f64 / 101.0).collect();
        let mut y = vec![0.0; a.nrows];
        let (med, min) = bench(|| a32.matvec_into(&x, &mut y), if a32.nnz() > 3_000_000 { 60 } else { 300 });
        println!("{:>9} rows {:>10} nnz: {med:9.1} / {min:9.1}", a32.nrows, a32.nnz());
    }
    for (label, order, n) in [("P1 Laplace n=52", 1usize, 52usize), ("P2 Laplace n=22", 2, 22)] {
        let m = SimplexMesh::block(1.0, 1.0, 1.0, n, n, n);
        let a: Csr = assemble_laplace(&m, order, &vec![1.0; m.n_cells()]);
        println!("\n== {label}: {} rows, {} nnz ==", a.nrows, a.nnz());
        let t = Instant::now();
        let amg = Amg::new(&a, vec![vec![1.0; a.nrows]], 1, &AmgOptions::default()).unwrap();
        println!("AMG setup {:.2} s, levels {:?}", t.elapsed().as_secs_f64(), amg.sizes());
        let b: Vec<f64> = (0..a.nrows).map(|i| ((i * 31) % 17) as f64 / 17.0).collect();
        let mut best = f64::MAX;
        for _ in 0..5 {
            let t = Instant::now();
            let r = pcg(&a, &b, &amg, 1e-30, 25);
            best = best.min(t.elapsed().as_secs_f64() / r.iterations as f64);
        }
        let t = Instant::now();
        let r = pcg(&a, &b, &amg, 1e-8, 200);
        println!("AMG-PCG: {:.0} us per iteration (best of 5 x 25); to 1e-8: {} iterations, {:.3} s", best * 1e6, r.iterations, t.elapsed().as_secs_f64());
    }
}
