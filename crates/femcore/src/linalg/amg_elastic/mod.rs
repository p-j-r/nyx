//! Smoothed-aggregation algebraic multigrid and a preconditioned conjugate gradient solver for
//! large symmetric positive definite systems (3D elasticity and diffusion with millions of
//! unknowns, where a sparse direct factorisation needs too much memory).
//!
//! The method is that of Vanek, Mandel and Brezina (1996): nodes are grouped into aggregates by
//! a greedy pass over the strength-of-connection graph, a tentative prolongator is built from the
//! near-nullspace (constants for diffusion, rigid-body modes for elasticity) by a local QR
//! factorisation per aggregate, and smoothed by one damped-Jacobi step. Coarse matrices are
//! Galerkin products; the smoother is symmetric Gauss-Seidel, so one V-cycle is a symmetric
//! positive definite preconditioner for conjugate gradients.

mod single;
use crate::linalg::csr::RCsr as Csr;
use crate::linalg::driver::DirectSolver;
use single::{dot, Csr32};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct AmgOptions {
    pub theta: f64,
    /// pre- and post-smoothing sweeps (symmetric Gauss-Seidel)
    pub sweeps: usize,
    /// Chebyshev polynomial smoother of this degree (0: symmetric Gauss-Seidel). The Chebyshev
    /// smoother is built from matrix-vector products only and runs on all cores.
    pub chebyshev: usize,
    pub max_coarse: usize,
    pub max_levels: usize,
    /// coarse-level visits per cycle: 1 is a V-cycle, 2 a W-cycle
    pub gamma: usize,
    /// steps of prolongator smoothing (1 is the classical smoothed aggregation)
    pub smooth_p: usize,
    /// use the inverse of the node blocks (3x3 for elasticity) instead of the scalar diagonal in the
    /// smoother and the prolongator smoothing
    pub block: bool,
    /// for P2 systems do not use the P1 subspace as the first coarse space (experiments)
    pub no_p1_coarse: bool,
    /// steps of prolongator smoothing applied to the P2 -> P1 interpolation
    pub smooth_lin: usize,
}

impl Default for AmgOptions {
    fn default() -> Self {
        AmgOptions { theta: 0.04, sweeps: 1, chebyshev: 3, max_coarse: 1500, max_levels: 12, gamma: 1, smooth_p: 1, block: true, no_p1_coarse: false, smooth_lin: 0 }.with_env()
    }
}

impl AmgOptions {
    /// Overrides from the environment (NYX_AMG_THETA, NYX_AMG_SWEEPS), for experiments.
    pub fn with_env(mut self) -> Self {
        if let Some(v) = std::env::var("NYX_AMG_THETA").ok().and_then(|v| v.parse().ok()) {
            self.theta = v;
        }
        if let Some(v) = std::env::var("NYX_AMG_CHEB").ok().and_then(|v| v.parse().ok()) {
            self.chebyshev = v;
        }
        if let Some(v) = std::env::var("NYX_AMG_GAMMA").ok().and_then(|v| v.parse().ok()) {
            self.gamma = v;
        }
        if let Some(v) = std::env::var("NYX_AMG_SMOOTHP").ok().and_then(|v| v.parse().ok()) {
            self.smooth_p = v;
        }
        if let Some(v) = std::env::var("NYX_AMG_SMOOTHLIN").ok().and_then(|v| v.parse().ok()) {
            self.smooth_lin = v;
        }
        if std::env::var("NYX_AMG_NOP1").is_ok() {
            self.no_p1_coarse = true;
        }
        if let Ok(v) = std::env::var("NYX_AMG_BLOCK") {
            self.block = v != "0";
        }
        if let Some(v) = std::env::var("NYX_AMG_SWEEPS").ok().and_then(|v| v.parse().ok()) {
            self.sweeps = v;
        }
        self
    }
}

struct Level {
    // single-precision copies: the cycle is only a preconditioner and is limited by memory bandwidth
    a: Csr32,
    p: Csr32,
    r: Csr32,
    diag: Vec<f64>,
    /// inverse of the diagonal blocks (row-major, `bs` x `bs` each) and the block size
    binv: Vec<f64>,
    bs: usize,
    /// estimate of the largest eigenvalue of D^-1 A
    lmax: f64,
}

pub struct Amg {
    gamma: usize,
    sweeps: usize,
    cheb: usize,
    levels: Vec<Level>,
    coarse_rows: usize,
    coarse_a32: Csr32,
    coarse: DirectSolver<f64>,
}

/// Sparse matrix product, on all cores by blocks of rows of `a` for large matrices.
/// Sparse matrix product, on all cores by blocks of rows of `a` for large matrices.
pub fn spmm(a: &Csr, b: &Csr) -> Csr {
    assert_eq!(a.ncols(), b.nrows());
    let nt = if a.nrows() < 20_000 { 1 } else { crate::par::threads() };
    let chunk = a.nrows().div_ceil(nt).max(1);
    let block = |r0: usize, r1: usize| -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let mut counts = Vec::with_capacity(r1 - r0);
        let (mut indices, mut data): (Vec<usize>, Vec<f64>) = (vec![], vec![]);
        let mut mark = vec![usize::MAX; b.ncols()];
        let mut acc = vec![0.0; b.ncols()];
        let mut touched: Vec<usize> = vec![];
        for i in r0..r1 {
            touched.clear();
            for ka in a.indptr()[i]..a.indptr()[i + 1] {
                let (k, av) = (a.indices()[ka] as usize, a.values()[ka]);
                for kb in b.indptr()[k]..b.indptr()[k + 1] {
                    let j = b.indices()[kb] as usize;
                    if mark[j] != i {
                        mark[j] = i;
                        acc[j] = 0.0;
                        touched.push(j);
                    }
                    acc[j] += av * b.values()[kb];
                }
            }
            touched.sort_unstable();
            for &j in &touched {
                indices.push(j);
                data.push(acc[j]);
            }
            counts.push(touched.len());
        }
        (counts, indices, data)
    };
    let blocks: Vec<(Vec<usize>, Vec<usize>, Vec<f64>)> = if nt == 1 {
        vec![block(0, a.nrows())]
    } else {
        (0..a.nrows().div_ceil(chunk)).into_par_iter().map(|t| block(t * chunk, ((t + 1) * chunk).min(a.nrows()))).collect()
    };
    let mut indptr = Vec::with_capacity(a.nrows() + 1);
    indptr.push(0usize);
    let (mut indices, mut data) = (Vec::new(), Vec::new());
    for (counts, idx, val) in blocks {
        for c in counts {
            indptr.push(indptr.last().unwrap() + c);
        }
        indices.extend(idx);
        data.extend(val);
    }
    Csr::from_raw(a.nrows(), b.ncols(), indptr, indices, data)
}

/// Inverse of a small dense matrix (Gauss-Jordan, partial pivoting); `None` when singular.
fn invert_small(m: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut a = m.to_vec();
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    let scale = (0..n).map(|i| a[i * n + i].abs()).fold(0.0, f64::max);
    for c in 0..n {
        let piv = (c..n).max_by(|&x, &y| a[x * n + c].abs().partial_cmp(&a[y * n + c].abs()).unwrap())?;
        if a[piv * n + c].abs() <= 1e-13 * scale {
            return None;
        }
        if piv != c {
            for j in 0..n {
                a.swap(c * n + j, piv * n + j);
                inv.swap(c * n + j, piv * n + j);
            }
        }
        let d = 1.0 / a[c * n + c];
        for j in 0..n {
            a[c * n + j] *= d;
            inv[c * n + j] *= d;
        }
        for r in 0..n {
            if r != c {
                let f = a[r * n + c];
                if f != 0.0 {
                    for j in 0..n {
                        a[r * n + j] -= f * a[c * n + j];
                        inv[r * n + j] -= f * inv[c * n + j];
                    }
                }
            }
        }
    }
    Some(inv)
}

/// Inverses of the `bs` x `bs` diagonal blocks of `a` (row-major, block after block); with
/// `bs == 1` or a singular block the reciprocal of the scalar diagonal is used.
fn block_inverses(a: &Csr, bs: usize) -> Vec<f64> {
    let nb = a.nrows() / bs;
    let mut out = vec![0.0; nb * bs * bs];
    for i in 0..nb {
        let mut m = vec![0.0; bs * bs];
        for r in 0..bs {
            let row = i * bs + r;
            for k in a.indptr()[row]..a.indptr()[row + 1] {
                let c = a.indices()[k] as usize;
                if c / bs == i {
                    m[r * bs + c % bs] = a.values()[k];
                }
            }
        }
        let inv = if bs > 1 { invert_small(&m, bs) } else { None };
        let blk = &mut out[i * bs * bs..(i + 1) * bs * bs];
        match inv {
            Some(v) => blk.copy_from_slice(&v),
            None => {
                for r in 0..bs {
                    let d = m[r * bs + r];
                    blk[r * bs + r] = if d.abs() > 1e-300 { 1.0 / d } else { 1.0 };
                }
            }
        }
    }
    out
}

/// `z = D^-1 x` with the inverse diagonal blocks of `block_inverses`.
fn apply_binv(binv: &[f64], bs: usize, x: &[f64]) -> Vec<f64> {
    let mut z = vec![0.0; x.len()];
    for (i, blk) in binv.chunks(bs * bs).enumerate() {
        for r in 0..bs {
            let mut s = 0.0;
            for c in 0..bs {
                s += blk[r * bs + c] * x[i * bs + c];
            }
            z[i * bs + r] = s;
        }
    }
    z
}

/// The block diagonal `D^-1` as a sparse matrix, to scale rows of products.
fn binv_csr(binv: &[f64], bs: usize) -> Csr {
    let n = binv.len() / bs;
    let mut indptr = Vec::with_capacity(n + 1);
    let mut indices = Vec::with_capacity(binv.len());
    indptr.push(0);
    for (i, blk) in binv.chunks(bs * bs).enumerate() {
        for r in 0..bs {
            for c in 0..bs {
                indices.push(i * bs + c);
                let _ = (r, c);
            }
            indptr.push(indices.len());
        }
        let _ = blk;
    }
    Csr::from_raw(n, n, indptr, indices, binv.to_vec())
}

/// Spectral radius estimate of `D^-1 A` by the power method.
fn rho_dinv_a(a: &Csr, binv: &[f64], bs: usize) -> f64 {
    let n = a.nrows();
    let mut v: Vec<f64> = (0..n).map(|i| ((i * 2654435761usize) % 1000) as f64 / 1000.0 + 0.1).collect();
    let mut rho = 1.0;
    for _ in 0..15 {
        let av = a.matvec(&v);
        let w = apply_binv(binv, bs, &av);
        let nrm = dot(&w, &w).sqrt().max(1e-300);
        rho = nrm / dot(&v, &v).sqrt().max(1e-300);
        v = w.iter().map(|x| x / nrm).collect();
    }
    rho * 1.05
}

/// Greedy aggregation of the nodes of `a` (block size `bs`); returns the aggregate of each node.
fn aggregate(a: &Csr, bs: usize, theta: f64) -> (Vec<usize>, usize) {
    let nn = a.nrows() / bs;
    // block strengths: Frobenius norm of the block (i, j)
    let mut diag_norm = vec![0.0; nn];
    let mut nbrs: Vec<Vec<(usize, f64)>> = vec![vec![]; nn];
    for i in 0..nn {
        let mut acc: std::collections::BTreeMap<usize, f64> = Default::default();
        for r in i * bs..(i + 1) * bs {
            for k in a.indptr()[r]..a.indptr()[r + 1] {
                *acc.entry(a.indices()[k] as usize / bs).or_insert(0.0) += a.values()[k] * a.values()[k];
            }
        }
        for (j, s) in acc {
            if j == i {
                diag_norm[i] = s.sqrt();
            } else {
                nbrs[i].push((j, s.sqrt()));
            }
        }
    }
    let strong: Vec<Vec<usize>> = (0..nn)
        .map(|i| nbrs[i].iter().filter(|(j, s)| *s >= theta * (diag_norm[i] * diag_norm[*j]).sqrt()).map(|(j, _)| *j).collect())
        .collect();
    let mut agg = vec![usize::MAX; nn];
    let mut na = 0;
    // pass 1: roots whose strong neighbours are all free
    for i in 0..nn {
        if agg[i] == usize::MAX && strong[i].iter().all(|&j| agg[j] == usize::MAX) {
            agg[i] = na;
            for &j in &strong[i] {
                agg[j] = na;
            }
            na += 1;
        }
    }
    // pass 2: attach the rest to the neighbouring aggregate with the strongest link
    let snapshot = agg.clone();
    for i in 0..nn {
        if agg[i] == usize::MAX {
            let mut best: Option<(f64, usize)> = None;
            for &(j, s) in &nbrs[i] {
                if snapshot[j] != usize::MAX && strong[i].contains(&j) && best.is_none_or(|(bs_, _)| s > bs_) {
                    best = Some((s, snapshot[j]));
                }
            }
            if let Some((_, g)) = best {
                agg[i] = g;
            }
        }
    }
    // pass 3: whatever is left forms its own aggregates (with free neighbours)
    for i in 0..nn {
        if agg[i] == usize::MAX {
            agg[i] = na;
            for &j in &strong[i] {
                if agg[j] == usize::MAX {
                    agg[j] = na;
                }
            }
            na += 1;
        }
    }
    (agg, na)
}

/// Tentative prolongator and coarse nullspace from aggregates and the nullspace `b` (k columns,
/// each of length n = nn*bs).
fn tentative(agg: &[usize], na: usize, bs: usize, b: &[Vec<f64>]) -> (Csr, Vec<Vec<f64>>) {
    let k = b.len();
    let nn = agg.len();
    let mut members: Vec<Vec<usize>> = vec![vec![]; na];
    for (i, &g) in agg.iter().enumerate() {
        members[g].push(i);
    }
    let mut indptr = vec![0usize; nn * bs + 1];
    let mut indices: Vec<usize> = vec![];
    let mut data: Vec<f64> = vec![];
    let mut coarse_b = vec![vec![0.0; na * k]; k];
    // rows of P are node-major: build per row after computing the aggregate blocks
    let mut qblocks: Vec<Vec<f64>> = Vec::with_capacity(na); // rows x k, row-major per aggregate
    for (g, mem) in members.iter().enumerate() {
        let rows = mem.len() * bs;
        let mut q = vec![0.0; rows * k]; // column-major copy of B restricted
        for (mi, &node) in mem.iter().enumerate() {
            for c in 0..bs {
                for j in 0..k {
                    q[j * rows + mi * bs + c] = b[j][node * bs + c];
                }
            }
        }
        // modified Gram-Schmidt: q <- Q, r upper triangular k x k
        let mut r = vec![0.0; k * k];
        for j in 0..k {
            let col0 = (0..rows).map(|t| q[j * rows + t].powi(2)).sum::<f64>().sqrt();
            for i in 0..j {
                let mut s = 0.0;
                for t in 0..rows {
                    s += q[i * rows + t] * q[j * rows + t];
                }
                r[i * k + j] = s;
                for t in 0..rows {
                    q[j * rows + t] -= s * q[i * rows + t];
                }
            }
            let nrm = (0..rows).map(|t| q[j * rows + t].powi(2)).sum::<f64>().sqrt();
            if nrm > 1e-8 * col0 && nrm > 0.0 {
                r[j * k + j] = nrm;
                for t in 0..rows {
                    q[j * rows + t] /= nrm;
                }
            } else {
                for t in 0..rows {
                    q[j * rows + t] = 0.0;
                }
            }
        }
        for i in 0..k {
            for j in 0..k {
                coarse_b[j][g * k + i] = r[i * k + j];
            }
        }
        // store Q row-major
        let mut qrow = vec![0.0; rows * k];
        for t in 0..rows {
            for j in 0..k {
                qrow[t * k + j] = q[j * rows + t];
            }
        }
        qblocks.push(qrow);
    }
    // assemble P row by row
    let mut pos_in_agg = vec![0usize; nn];
    for mem in &members {
        for (mi, &node) in mem.iter().enumerate() {
            pos_in_agg[node] = mi;
        }
    }
    let mut count = 0;
    for node in 0..nn {
        let g = agg[node];
        for c in 0..bs {
            let row = node * bs + c;
            let qrow = &qblocks[g];
            let t = pos_in_agg[node] * bs + c;
            for j in 0..k {
                let v = qrow[t * k + j];
                if v != 0.0 {
                    indices.push(g * k + j);
                    data.push(v);
                    count += 1;
                }
            }
            indptr[row + 1] = count;
        }
    }
    (Csr::from_raw(nn * bs, na * k, indptr, indices, data), coarse_b)
}

impl Amg {
    /// `nullspace`: the near-nullspace vectors in dof order (node-major, block size `bs`).
    pub fn new(a: &Csr, nullspace: Vec<Vec<f64>>, bs: usize, opts: &AmgOptions) -> Result<Amg, String> {
        Self::build(a, nullspace, bs, opts, None)
    }

    /// Like `new` for a quadratic (P2) discretisation: the first coarse space is the linear (P1)
    /// subspace, given by `lin`, the node-level interpolation (all nodes x vertex nodes, vertex
    /// nodes first) of linear finite elements into the quadratic ones. The aggregation hierarchy
    /// then starts from the P1 matrix, which is much better behaved than the P2 matrix itself.
    pub fn new_p2(a: &Csr, nullspace: Vec<Vec<f64>>, bs: usize, opts: &AmgOptions, lin: &Csr) -> Result<Amg, String> {
        Self::build(a, nullspace, bs, opts, Some(lin))
    }

    fn build(a: &Csr, nullspace: Vec<Vec<f64>>, bs: usize, opts: &AmgOptions, lin: Option<&Csr>) -> Result<Amg, String> {
        let mut levels: Vec<Level> = vec![];
        let mut tl = std::time::Instant::now();
        let _stage = crate::progress::stage(format!("building the multigrid hierarchy ({} unknowns)", a.nrows()));
        // the finest matrix is only borrowed (a copy would be as large as the matrix); coarser ones are owned
        let mut cur: std::borrow::Cow<Csr> = std::borrow::Cow::Borrowed(a);
        let mut b = nullspace;
        let mut bs = bs;
        loop {
            crate::progress::check()?;
            let k = b.len();
            if cur.nrows() <= opts.max_coarse || levels.len() + 1 >= opts.max_levels {
                break;
            }
            let first_linear = levels.is_empty() && lin.is_some();
            let (pt, bc, next_bs) = if first_linear {
                let j = lin.unwrap();
                // block expansion: node row r, component c  ->  vertex node column, component c
                let mut indptr = vec![0usize];
                let (mut indices, mut data) = (vec![], vec![]);
                for node in 0..j.nrows() {
                    for c in 0..bs {
                        for q in j.indptr()[node]..j.indptr()[node + 1] {
                            indices.push(j.indices()[q] as usize * bs + c);
                            data.push(j.values()[q]);
                        }
                        indptr.push(indices.len());
                    }
                }
                let p0 = Csr::from_raw(j.nrows() * bs, j.ncols() * bs, indptr, indices, data);
                // linear fields are reproduced exactly, so the coarse near-nullspace is the restriction to the vertices
                let bc: Vec<Vec<f64>> = b.iter().map(|v| v[..j.ncols() * bs].to_vec()).collect();
                (p0, bc, bs)
            } else {
                let (agg, na) = aggregate(&cur, bs, opts.theta * 0.5f64.powi(levels.len() as i32));
                crate::progress::lap_since(&mut tl, "      amg: aggregate");
                if na * k >= cur.nrows() * 9 / 10 {
                    break; // no useful coarsening
                }
                let (pt, bc) = tentative(&agg, na, bs, &b);
                crate::progress::lap_since(&mut tl, "      amg: tentative");
                (pt, bc, k)
            };
            let diag: Vec<f64> = cur.diagonal().iter().map(|d| if d.abs() > 1e-300 { *d } else { 1.0 }).collect();
            let blk = if opts.block { bs } else { 1 };
            let binv = block_inverses(&cur, blk);
            let lmax = rho_dinv_a(&cur, &binv, blk);
            let omega = 4.0 / 3.0 / lmax;
            crate::progress::lap_since(&mut tl, "      amg: block inverses, rho");
            let dinv = binv_csr(&binv, blk);
            // P = (I - omega D^-1 A)^s Pt
            let mut p = pt;
            for _ in 0..if first_linear { opts.smooth_lin } else { opts.smooth_p.max(1) } {
                let mut scaled = spmm(&dinv, &spmm(&cur, &p));
                scaled.values_mut().iter_mut().for_each(|v| *v *= omega);
                p = p.add_scaled(&scaled, -1.0);
            }
            crate::progress::lap_since(&mut tl, "      amg: smooth P");
            let r = p.transpose();
            let ap = spmm(&cur, &p);
            let mut ac = spmm(&r, &ap);
            crate::progress::lap_since(&mut tl, "      amg: Galerkin product");
            // empty coarse dofs (rank-deficient aggregates) get an identity row
            let empty: Vec<usize> = (0..ac.nrows()).filter(|&i| ac.get(i, i).abs() < 1e-300).collect();
            if !empty.is_empty() {
                let mut d = vec![0.0; ac.nrows()];
                for &i in &empty {
                    d[i] = 1.0;
                }
                ac = ac.add_scaled(&Csr::from_diag(&d), 1.0);
            }
            levels.push(Level { a: Csr32::from_csr(&cur), p: Csr32::from_csr(&p), r: Csr32::from_csr(&r), diag, binv, bs: blk, lmax });
            crate::progress::lap_since(&mut tl, "      amg: convert level");
            cur = std::borrow::Cow::Owned(ac);
            b = bc;
            bs = next_bs;
        }
        let coarse = DirectSolver::new(&cur, true)?;
        crate::progress::note(format!("multigrid levels (unknowns): {:?} + coarsest {}", levels.iter().map(|l| l.a.nrows).collect::<Vec<_>>(), cur.nrows()));
        Ok(Amg { gamma: opts.gamma.max(1), cheb: opts.chebyshev, sweeps: opts.sweeps.max(1), levels, coarse_rows: cur.nrows(), coarse_a32: Csr32::from_csr(&cur), coarse })
    }

    pub fn n_levels(&self) -> usize {
        self.levels.len() + 1
    }

    /// Sizes of the matrices on each level, finest first.
    pub fn sizes(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.levels.iter().map(|l| l.a.nrows).collect();
        v.push(self.coarse_rows);
        v
    }

    fn sgs(a: &Csr32, diag: &[f64], b: &[f64], x: &mut [f64], forward: bool) {
        let n = a.nrows;
        let order: Box<dyn Iterator<Item = usize>> = if forward { Box::new(0..n) } else { Box::new((0..n).rev()) };
        for i in order {
            let mut s = b[i];
            for k in a.indptr[i]..a.indptr[i + 1] {
                let j = a.indices[k] as usize;
                if j != i {
                    s -= a.data[k] as f64 * x[j];
                }
            }
            x[i] = s / diag[i];
        }
    }

    /// `zero`: the guess `x` is zero, which saves the first matrix-vector product.
    fn smooth(&self, lev: &Level, b: &[f64], x: &mut [f64], forward: bool, zero: bool) {
        if self.cheb == 0 {
            for _ in 0..self.sweeps {
                Self::sgs(&lev.a, &lev.diag, b, x, forward);
            }
            return;
        }
        // Chebyshev iteration on D^-1 A over [lmax/30, 1.1 lmax]
        let n = x.len();
        let (hi, lo) = (1.1 * lev.lmax, 1.1 * lev.lmax / 30.0);
        let (theta, delta) = (0.5 * (hi + lo), 0.5 * (hi - lo));
        let sigma = theta / delta;
        let mut rho_old = 1.0 / sigma;
        let resid = |x: &[f64]| -> Vec<f64> {
            let ax = lev.a.matvec(x);
            let rr: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
            apply_binv(&lev.binv, lev.bs, &rr)
        };
        let first = if zero { apply_binv(&lev.binv, lev.bs, b) } else { resid(x) };
        let mut d: Vec<f64> = first.iter().map(|r| r / theta).collect();
        for i in 0..n {
            x[i] += d[i];
        }
        for _ in 1..self.cheb {
            let rho_new = 1.0 / (2.0 * sigma - rho_old);
            let r = resid(x);
            for i in 0..n {
                d[i] = rho_new * rho_old * d[i] + 2.0 * rho_new / delta * r[i];
                x[i] += d[i];
            }
            rho_old = rho_new;
        }
    }

    fn cycle(&self, l: usize, b: &[f64], x: &mut [f64]) {
        if l == self.levels.len() {
            let s = self.coarse.solve(b);
            x.copy_from_slice(&s);
            return;
        }
        let lev = &self.levels[l];
        x.iter_mut().for_each(|v| *v = 0.0);
        self.smooth(lev, b, x, true, true);
        let ax = lev.a.matvec(x);
        let res: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let rc = lev.r.matvec(&res);
        let mut xc = vec![0.0; rc.len()];
        self.cycle(l + 1, &rc, &mut xc);
        if self.gamma > 1 {
            let ac = if l + 1 == self.levels.len() { &self.coarse_a32 } else { &self.levels[l + 1].a };
            for _ in 1..self.gamma {
                let axc = ac.matvec(&xc);
                let r2: Vec<f64> = rc.iter().zip(&axc).map(|(p, q)| p - q).collect();
                let mut dx = vec![0.0; rc.len()];
                self.cycle(l + 1, &r2, &mut dx);
                xc.iter_mut().zip(&dx).for_each(|(v, c)| *v += c);
            }
        }
        let corr = lev.p.matvec(&xc);
        x.iter_mut().zip(&corr).for_each(|(v, c)| *v += c);
        self.smooth(lev, b, x, false, false);
    }

    /// One V-cycle applied to `b` from a zero guess.
    pub fn apply(&self, b: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; b.len()];
        self.cycle(0, b, &mut x);
        x
    }
}

pub struct PcgResult {
    pub x: Vec<f64>,
    pub iterations: usize,
    pub relative_residual: f64,
}

/// Conjugate gradients preconditioned by one AMG V-cycle.
pub fn pcg(a: &Csr, b: &[f64], amg: &Amg, tol: f64, max_iter: usize) -> PcgResult {
    let n = a.nrows();
    let trace = std::env::var("NYX_AMG_TRACE").is_ok();
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let bn = dot(b, b).sqrt().max(1e-300);
    let mut z = amg.apply(&r);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    for it in 0..max_iter {
        let ap = a.matvec(&p);
        let alpha = rz / dot(&p, &ap);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rel = dot(&r, &r).sqrt() / bn;
        if trace && it % 10 == 0 {
            eprintln!("  pcg it {it}: {rel:.3e}");
        }
        if it % 10 == 9 {
            crate::progress::note(format!("conjugate gradients iteration {}: relative residual {rel:.1e} (goal {tol:.0e})", it + 1));
        }
        if crate::progress::check().is_err() {
            // cancelled: report "not converged"; callers that fall back to a factorisation hit the cancel check there
            return PcgResult { x, iterations: it + 1, relative_residual: f64::INFINITY };
        }
        if rel < tol {
            return PcgResult { x, iterations: it + 1, relative_residual: rel };
        }
        z = amg.apply(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    let rel = dot(&r, &r).sqrt() / bn;
    PcgResult { x, iterations: max_iter, relative_residual: rel }
}

