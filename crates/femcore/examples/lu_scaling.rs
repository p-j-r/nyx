//! How the sparse factorisations of faer scale with the number of threads (NYX_THREADS=N): the complex LU used
//! by the mode solvers on a 2D matrix of 160,000 unknowns with a 25-point stencil (like a Nedelec P2 cross-section),
//! and the real Cholesky of a 3D 7-point matrix of 125,000 unknowns.
//!
//!   for %n in (1 2 4 8 20) do set NYX_THREADS=%n && cargo run --release -p femcore --example lu_scaling
use femcore::csparse::{c, CCoo, CSolver};
use femcore::linalg::SparseSolver;
use femcore::sparse::Coo;
use std::time::Instant;

fn main() {
    let threads = femcore::par::threads();
    // complex 2D: 400 x 400 grid, entries to all neighbours within +-2 in x and y
    let n = 400usize;
    let mut a = CCoo::new(n * n, n * n);
    for y in 0..n {
        for x in 0..n {
            let i = y * n + x;
            for dy in -2i64..=2 {
                for dx in -2i64..=2 {
                    let (xx, yy) = (x as i64 + dx, y as i64 + dy);
                    if xx < 0 || yy < 0 || xx >= n as i64 || yy >= n as i64 {
                        continue;
                    }
                    let j = yy as usize * n + xx as usize;
                    let v = if i == j { c(30.0, 0.5) } else { c(-1.0 / (1.0 + (dx * dx + dy * dy) as f64), 0.01) };
                    a.push_block(&[i], &[j], &[v]);
                }
            }
        }
    }
    let a = a.to_csr();
    let t = Instant::now();
    let _s = CSolver::new(&a).unwrap();
    println!("threads {threads:2}: complex LU 2D, {} unknowns, {} nnz: {:.2} s", a.nrows, a.data.len(), t.elapsed().as_secs_f64());
    // real 3D Cholesky: 50^3 grid, 7-point
    let m = 50usize;
    let mut k = Coo::new(m * m * m, m * m * m);
    let id = |x: usize, y: usize, z: usize| (z * m + y) * m + x;
    for z in 0..m {
        for y in 0..m {
            for x in 0..m {
                k.push(id(x, y, z), id(x, y, z), 6.1);
                for (dx, dy, dz) in [(1i64, 0i64, 0i64), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)] {
                    let (xx, yy, zz) = (x as i64 + dx, y as i64 + dy, z as i64 + dz);
                    if xx >= 0 && yy >= 0 && zz >= 0 && (xx as usize) < m && (yy as usize) < m && (zz as usize) < m {
                        k.push(id(x, y, z), id(xx as usize, yy as usize, zz as usize), -1.0);
                    }
                }
            }
        }
    }
    let k = k.to_csr();
    let t = Instant::now();
    let _s = SparseSolver::new(&k, true).unwrap();
    println!("threads {threads:2}: real Cholesky 3D, {} unknowns, {} nnz: {:.2} s", k.nrows, k.nnz(), t.elapsed().as_secs_f64());
}
