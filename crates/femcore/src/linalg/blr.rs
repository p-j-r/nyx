//! **Low-rank (block low-rank, BLR) multifrontal sparse LU** (port item S4, the sparse-solver half).
//!
//! A multifrontal factorisation eliminates the unknowns along a *nested-dissection tree*: each tree node (a
//! separator, or a small leaf subdomain) owns a dense **front** `F = [F₁₁ F₁₂; F₂₁ F₂₂]` whose pivot block `F₁₁`
//! (the node's own unknowns) is factorised, whose panels `L₂₁ = F₂₁U₁₁⁻¹` and `U₁₂ = L₁₁⁻¹F₁₂` couple the node to the
//! ancestor separators, and whose Schur complement `F₂₂ − L₂₁U₁₂` (the *update*) is added into the parent's front.
//! The big fronts are the top separators: in 3-D a separator is a surface of `n^{2/3}` unknowns and its front costs
//! `O(n²)` flops and `O(n^{4/3})` memory.
//!
//! **Why low rank.** The panels `L₂₁`, `U₁₂` describe the coupling between two *pieces of surfaces*. Blocks of them
//! that join clusters far from each other are numerically low-rank (the same smoothness that makes the Schur complement
//! of [`crate::linalg::lowrank`]'s study compress), so cutting the panels into `b × b` blocks and replacing every block that is
//! compressible by `U Vᴴ` (found by full-pivoting adaptive cross approximation, [`aca_full`]) reduces the *stored* factor,
//! the work of the triangular solves, and — if the Schur update were also done from the compressed blocks — the flops
//! (see below). This is the BLR format of Amestoy, Ashcraft, Boiteau, Buttari, L'Excellent & Weisbecker (2015),
//! written here from the paper and from the building blocks already in this crate.
//!
//! # What is implemented
//! * Ordering, tree and symbolic structure: [`nested_dissection_tree`] on the
//!   pattern of `A + Aᵀ`; the front of a node has its own unknowns first, then the sorted ancestors it touches. The
//!   symbolic object ([`BlrSymbolic`]) depends only on the pattern, so it is built once and reused for every shift.
//! * Numeric phase: dense fronts assembled from the permuted matrix and the children's updates (extend-add), pivot
//!   block `LU` with **partial pivoting inside the pivot block** (`faer`'s dense kernels, `faer_glue.rs`), panels by
//!   triangular solves, update by a matrix product (BLAS-3). Real or complex scalars. The matrix must be
//!   *structurally* symmetric (values may be unsymmetric).
//! * Compression: only fronts with pivot size `p ≥ min_front` are compressed; the packed `L₁₁U₁₁` of the pivot block and the two
//!   panels are cut in `block × block` blocks (the diagonal blocks stay dense); a block is stored as `U Vᴴ` if the cross approximation reaches `‖E‖_F ≤ tol ‖block‖_F` with
//!   `rank · (rows + cols) ≤ accept_ratio · rows · cols`, and densely otherwise. `tol = 0` is the exact solver.
//! * Solves ([`BlrLu::solve_in_place`]): forward over the tree, backward over the reversed tree, the panels
//!   applied block by block. [`BlrLu`] is a [`Preconditioner`]; [`ShiftedBlrLu`] a [`ShiftedSolver`] (symbolic data reused
//!   across shifts); [`BlrLu::solve_refined`] adds iterative refinement.
//!
//! # What this does *not* do (and why it matters)
//! The compression is applied to the panels **after** they were computed densely, and the Schur updates are computed
//! from the **exact** dense panels. Consequences: (1) the stored factor and the solves benefit fully, (2) the
//! factorisation flops are *not* reduced and the peak memory (one dense front at a time plus stacked updates) is not
//! reduced, (3) errors do not compound from level to level (every stored panel is within `tol` of the exact one), so
//! the solve error is `O(tol · levels)`. The flop-saving variants (compress the assembled panels first and do the
//! triangular solves and the update on the thin factors, MUMPS "FCSU/UCFS") need left-looking panel solves on block
//! columns and compressed contribution blocks; they are the natural next step. No look-ahead pivoting across fronts:
//! a pivot block that is exactly singular is an error (`Singular`), a small one gives a poor factor (use refinement or
//! GMRES with this as the preconditioner).

use crate::linalg::csr::Csr;
use crate::linalg::direct::{FactorBackend, FactorStats, ShiftedSolver, SolverError, Used};
use crate::linalg::faer_glue::{FaerScalar, dense_lu, gemm_sub, solve_unit_lower, solve_upper_right};
use crate::linalg::linop::Preconditioner;
use crate::linalg::lowrank::{LowRank, aca_full};
use crate::linalg::ordering::nested_dissection_tree;
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{axpy, dotc, norm2};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

#[derive(Debug, Clone)]
pub struct BlrOptions {
    /// Relative accuracy of every compressed block; `0` disables compression (exact multifrontal LU).
    pub tol: f64,
    /// Nested-dissection subdomains of at most this many unknowns become leaf fronts.
    pub leaf: usize,
    /// Block size `b` of the BLR partition of the panels.
    pub block: usize,
    /// Compress the factors of a front only if its pivot size `p` is at least this large.
    pub min_front: usize,
    /// A block is kept low-rank only if `rank·(rows+cols) ≤ accept_ratio·rows·cols`.
    pub accept_ratio: f64,
}

impl Default for BlrOptions {
    fn default() -> Self {
        Self { tol: 1e-8, leaf: 64, block: 64, min_front: 128, accept_ratio: 0.7 }
    }
}

/// Memory and cost figures of a factorisation.
#[derive(Debug, Clone, Copy, Default)]
pub struct BlrStats {
    /// Entries an exact multifrontal factor stores (`Σ p² + 2pu`).
    pub exact_entries: usize,
    /// Entries stored by this factorisation (dense pivot blocks + panels, compressed where possible).
    pub stored_entries: usize,
    pub fronts: usize,
    /// Pivot and update size of the largest front.
    pub largest_front: (usize, usize),
    /// Panel blocks examined for compression, how many were stored low-rank, and the largest rank.
    pub blocks: usize,
    pub low_rank_blocks: usize,
    pub max_rank: usize,
    pub analyse_secs: f64,
    pub factor_secs: f64,
}

impl BlrStats {
    /// `stored_entries / exact_entries`.
    pub fn compression(&self) -> f64 {
        self.stored_entries as f64 / self.exact_entries.max(1) as f64
    }
}

struct SymNode {
    start: usize,
    end: usize,
    children: Vec<usize>,
    /// Sorted positions (in the elimination order) of the ancestors this node's front touches.
    strct: Vec<usize>,
}

/// Pattern-only data of the factorisation: ordering, tree and front structures. Reusable for every matrix with the
/// same (or a smaller) pattern.
pub struct BlrSymbolic {
    n: usize,
    /// `order[k]` = original unknown eliminated `k`-th.
    order: Vec<usize>,
    /// `pinv[i]` = position of original unknown `i`.
    pinv: Vec<usize>,
    nodes: Vec<SymNode>,
    exact_entries: usize,
}

fn add_pattern<T: Scalar>(adj: &mut [Vec<u32>], a: &Csr<T>) {
    for r in 0..a.nrows() {
        for &c in a.row(r).0 {
            if c as usize != r {
                adj[r].push(c);
                adj[c as usize].push(r as u32);
            }
        }
    }
}

impl BlrSymbolic {
    /// Analyses the union of the patterns of `matrices` (all `n × n`).
    pub fn analyse<T: Scalar>(matrices: &[&Csr<T>], leaf: usize) -> Result<Self, SolverError> {
        let n = matrices.first().ok_or(SolverError::InvalidArgument("at least one matrix"))?.nrows();
        if matrices.iter().any(|m| m.nrows() != n || m.ncols() != n) {
            return Err(SolverError::InvalidArgument("matrices must be square of the same size"));
        }
        let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
        for m in matrices {
            add_pattern(&mut adj, m);
        }
        for l in &mut adj {
            l.sort_unstable();
            l.dedup();
        }
        let tree = nested_dissection_tree(&adj, leaf);
        let order = tree.order;
        let mut pinv = vec![0usize; n];
        for (k, &o) in order.iter().enumerate() {
            pinv[o] = k;
        }
        let mut nodes: Vec<SymNode> = tree
            .nodes
            .iter()
            .map(|t| SymNode { start: t.start, end: t.end, children: t.children.clone(), strct: Vec::new() })
            .collect();
        let mut mark = vec![usize::MAX; n];
        let mut exact = 0usize;
        for t in 0..nodes.len() {
            let (start, end) = (nodes[t].start, nodes[t].end);
            let mut strct: Vec<usize> = Vec::new();
            for i in start..end {
                for &v in &adj[order[i]] {
                    let j = pinv[v as usize];
                    if j >= end && mark[j] != t {
                        mark[j] = t;
                        strct.push(j);
                    }
                }
            }
            for &c in &nodes[t].children {
                for &j in &nodes[c].strct {
                    if j >= end && mark[j] != t {
                        mark[j] = t;
                        strct.push(j);
                    }
                }
            }
            strct.sort_unstable();
            let (p, u) = (end - start, strct.len());
            exact += p * p + 2 * p * u;
            nodes[t].strct = strct;
        }
        Ok(Self { n, order, pinv, nodes, exact_entries: exact })
    }

    pub fn dim(&self) -> usize {
        self.n
    }

    /// Number of fronts.
    pub fn num_fronts(&self) -> usize {
        self.nodes.len()
    }

    /// `(pivot size, update size)` of every front.
    pub fn front_sizes(&self) -> Vec<(usize, usize)> {
        self.nodes.iter().map(|t| (t.end - t.start, t.strct.len())).collect()
    }

    /// Entries of the exact multifrontal factor.
    pub fn exact_entries(&self) -> usize {
        self.exact_entries
    }
}

enum Block<T: Scalar> {
    Dense(Vec<T>), // column-major
    Low(LowRank<T>),
}

/// A matrix cut into blocks, each dense or low-rank; used for the panels `L₂₁` (u × p) and `U₁₂` (p × u).
struct Panel<T: Scalar> {
    rb: Vec<usize>,
    cb: Vec<usize>,
    blocks: Vec<Block<T>>,
}

fn partition(n: usize, b: usize) -> Vec<usize> {
    let mut v = vec![0];
    let mut s = 0;
    while s < n {
        let mut e = (s + b).min(n);
        if n - e < b / 2 {
            e = n; // merge a short last block
        }
        v.push(e);
        s = e;
    }
    v
}

impl<T: Scalar> Panel<T> {
    fn dense(data: Vec<T>, rows: usize, cols: usize) -> Self {
        if rows == 0 || cols == 0 {
            return Self { rb: vec![0], cb: vec![0], blocks: Vec::new() };
        }
        Self { rb: vec![0, rows], cb: vec![0, cols], blocks: vec![Block::Dense(data)] }
    }

    fn stored(&self) -> usize {
        self.blocks
            .iter()
            .map(|b| match b {
                Block::Dense(d) => d.len(),
                Block::Low(l) => l.storage(),
            })
            .sum()
    }

    /// `y -= P x`.
    fn apply_sub(&self, x: &[T], y: &mut [T]) {
        let nbc = self.cb.len() - 1;
        for (idx, blk) in self.blocks.iter().enumerate() {
            let (bi, bj) = (idx / nbc, idx % nbc);
            let (r0, r1, c0, c1) = (self.rb[bi], self.rb[bi + 1], self.cb[bj], self.cb[bj + 1]);
            apply_block_sub(blk, &x[c0..c1], &mut y[r0..r1]);
        }
    }
}

/// Stores the `mr × nc` block `blk` as `U Vᴴ` if a cross approximation reaches `tol` with an acceptable rank, else densely.
fn compress_block<T: Scalar>(blk: Vec<T>, mr: usize, nc: usize, o: &BlrOptions, st: &mut BlrStats) -> Block<T> {
    st.blocks += 1;
    let norm0 = norm2(&blk);
    if norm0 == 0.0 {
        st.low_rank_blocks += 1;
        return Block::Low(LowRank::zero(mr, nc));
    }
    let maxr = ((o.accept_ratio * (mr * nc) as f64) / (mr + nc) as f64) as usize;
    if maxr >= 1 {
        let mut work = blk.clone();
        let lr = aca_full(&mut work, mr, nc, o.tol, maxr);
        if lr.r <= maxr && norm2(&work) <= o.tol * norm0 * (1.0 + 1e-9) {
            st.low_rank_blocks += 1;
            st.max_rank = st.max_rank.max(lr.r);
            return Block::Low(lr);
        }
    }
    Block::Dense(blk)
}

fn extract<T: Scalar>(data: &[T], ld: usize, r0: usize, r1: usize, c0: usize, c1: usize) -> Vec<T> {
    let mut blk = Vec::with_capacity((r1 - r0) * (c1 - c0));
    for c in c0..c1 {
        blk.extend_from_slice(&data[c * ld + r0..c * ld + r1]);
    }
    blk
}

/// Cuts `data` (`rows × cols`, column-major) into blocks and compresses what can be compressed.
fn compress_panel<T: Scalar>(data: &[T], rows: usize, cols: usize, o: &BlrOptions, st: &mut BlrStats) -> Panel<T> {
    let (rb, cb) = (partition(rows, o.block), partition(cols, o.block));
    let mut blocks = Vec::with_capacity((rb.len() - 1) * (cb.len() - 1));
    for bi in 0..rb.len() - 1 {
        for bj in 0..cb.len() - 1 {
            let (r0, r1, c0, c1) = (rb[bi], rb[bi + 1], cb[bj], cb[bj + 1]);
            blocks.push(compress_block(extract(data, rows, r0, r1, c0, c1), r1 - r0, c1 - c0, o, st));
        }
    }
    Panel { rb, cb, blocks }
}

/// `y ← L⁻¹ y` for the unit lower triangle of the packed `n × n` matrix `d`.
fn unit_lower_solve<T: Scalar>(d: &[T], n: usize, y: &mut [T]) {
    for j in 0..n {
        let yj = y[j];
        if yj != T::ZERO {
            axpy(-yj, &d[j * n + j + 1..(j + 1) * n], &mut y[j + 1..]);
        }
    }
}

/// `y ← U⁻¹ y` for the upper triangle (with diagonal) of the packed `n × n` matrix `d`.
fn upper_solve<T: Scalar>(d: &[T], n: usize, y: &mut [T]) {
    for j in (0..n).rev() {
        y[j] = y[j] / d[j * n + j];
        let yj = y[j];
        if yj != T::ZERO {
            axpy(-yj, &d[j * n..j * n + j], &mut y[..j]);
        }
    }
}

/// The packed `LU` of a pivot block: dense, or cut in blocks with the off-diagonal blocks of `L` and `U` compressed.
enum PivotLu<T: Scalar> {
    Dense(Vec<T>),
    Blocked {
        bb: Vec<usize>,
        diag: Vec<Vec<T>>,
        /// `(i, j)` with `i > j` at `i(i−1)/2 + j`
        lower: Vec<Block<T>>,
        /// `(i, j)` with `i < j` at `j(j−1)/2 + i`
        upper: Vec<Block<T>>,
    },
}

impl<T: Scalar> PivotLu<T> {
    fn blocked(lu: &[T], p: usize, o: &BlrOptions, st: &mut BlrStats) -> Self {
        let bb = partition(p, o.block);
        let nb = bb.len() - 1;
        let mut diag = Vec::with_capacity(nb);
        let mut lower = Vec::new();
        let mut upper: Vec<Block<T>> = Vec::new();
        for i in 0..nb {
            diag.push(extract(lu, p, bb[i], bb[i + 1], bb[i], bb[i + 1]));
            for j in 0..i {
                lower.push(compress_block(
                    extract(lu, p, bb[i], bb[i + 1], bb[j], bb[j + 1]),
                    bb[i + 1] - bb[i],
                    bb[j + 1] - bb[j],
                    o,
                    st,
                ));
            }
        }
        for j in 0..nb {
            for i in 0..j {
                upper.push(compress_block(
                    extract(lu, p, bb[i], bb[i + 1], bb[j], bb[j + 1]),
                    bb[i + 1] - bb[i],
                    bb[j + 1] - bb[j],
                    o,
                    st,
                ));
            }
        }
        Self::Blocked { bb, diag, lower, upper }
    }

    fn stored(&self) -> usize {
        let blk = |b: &Block<T>| match b {
            Block::Dense(d) => d.len(),
            Block::Low(l) => l.storage(),
        };
        match self {
            Self::Dense(d) => d.len(),
            Self::Blocked { diag, lower, upper, .. } => {
                diag.iter().map(Vec::len).sum::<usize>()
                    + lower.iter().map(blk).sum::<usize>()
                    + upper.iter().map(blk).sum::<usize>()
            }
        }
    }

    /// `y ← L⁻¹ y`.
    fn forward(&self, p: usize, y: &mut [T]) {
        match self {
            Self::Dense(d) => unit_lower_solve(d, p, y),
            Self::Blocked { bb, diag, lower, .. } => {
                for i in 0..bb.len() - 1 {
                    let (head, tail) = y.split_at_mut(bb[i]);
                    let yi = &mut tail[..bb[i + 1] - bb[i]];
                    for j in 0..i {
                        apply_block_sub(&lower[i * (i - 1) / 2 + j], &head[bb[j]..bb[j + 1]], yi);
                    }
                    unit_lower_solve(&diag[i], bb[i + 1] - bb[i], yi);
                }
            }
        }
    }

    /// `y ← U⁻¹ y`.
    fn backward(&self, p: usize, y: &mut [T]) {
        match self {
            Self::Dense(d) => upper_solve(d, p, y),
            Self::Blocked { bb, diag, upper, .. } => {
                let nb = bb.len() - 1;
                for i in (0..nb).rev() {
                    let (head, tail) = y.split_at_mut(bb[i + 1]);
                    let yi = &mut head[bb[i]..];
                    for j in i + 1..nb {
                        apply_block_sub(
                            &upper[j * (j - 1) / 2 + i],
                            &tail[bb[j] - bb[i + 1]..bb[j + 1] - bb[i + 1]],
                            yi,
                        );
                    }
                    upper_solve(&diag[i], bb[i + 1] - bb[i], yi);
                }
            }
        }
    }
}

/// `y -= B x` for one block.
fn apply_block_sub<T: Scalar>(blk: &Block<T>, x: &[T], y: &mut [T]) {
    match blk {
        Block::Dense(d) => {
            let mr = y.len();
            for (c, &xc) in x.iter().enumerate() {
                axpy(-xc, &d[c * mr..(c + 1) * mr], y);
            }
        }
        Block::Low(l) => {
            for k in 0..l.r {
                let c = dotc(l.vcol(k), x);
                axpy(-c, l.ucol(k), y);
            }
        }
    }
}

struct Front<T: Scalar> {
    p: usize,
    /// Packed `L₁₁U₁₁` of the pivot block and its row permutation `perm`.
    piv: PivotLu<T>,
    perm: Vec<usize>,
    l21: Panel<T>,
    u12: Panel<T>,
}

/// The numeric factors of a [`BlrSymbolic`] for one matrix.
pub struct BlrFactors<T: Scalar> {
    fronts: Vec<Front<T>>,
    stats: BlrStats,
}

impl<T: Scalar> BlrFactors<T> {
    pub fn stats(&self) -> &BlrStats {
        &self.stats
    }
}

fn permuted<T: Scalar>(sym: &BlrSymbolic, a: &Csr<T>) -> Csr<T> {
    let mut trip = Vec::with_capacity(a.nnz());
    for r in 0..a.nrows() {
        let (cols, vals) = a.row(r);
        for (&c, &v) in cols.iter().zip(vals) {
            trip.push((sym.pinv[r], sym.pinv[c as usize], v));
        }
    }
    Csr::from_triplets(a.nrows(), a.ncols(), &trip).expect("valid permuted triplets")
}

/// Numeric factorisation of `a` on the structure `sym`.
pub fn factor_numeric<T: FaerScalar>(
    sym: &BlrSymbolic,
    a: &Csr<T>,
    opts: &BlrOptions,
) -> Result<BlrFactors<T>, SolverError> {
    let t0 = Instant::now();
    let n = sym.n;
    if a.nrows() != n || a.ncols() != n {
        return Err(SolverError::InvalidArgument("matrix does not match the symbolic analysis"));
    }
    let ap = permuted(sym, a);
    let apt = ap.transpose();
    let mut updates: Vec<Option<Vec<T>>> = (0..sym.nodes.len()).map(|_| None).collect();
    let mut loc = vec![usize::MAX; n];
    let mut fronts: Vec<Front<T>> = Vec::with_capacity(sym.nodes.len());
    let mut st = BlrStats { exact_entries: sym.exact_entries, fronts: sym.nodes.len(), ..Default::default() };
    for (t, node) in sym.nodes.iter().enumerate() {
        let (start, end) = (node.start, node.end);
        let (p, u) = (end - start, node.strct.len());
        let m = p + u;
        for i in 0..p {
            loc[start + i] = i;
        }
        for (k, &g) in node.strct.iter().enumerate() {
            loc[g] = p + k;
        }
        let mut f = vec![T::ZERO; m * m];
        for i in start..end {
            let r = i - start;
            let (cols, vals) = ap.row(i);
            for (&j, &v) in cols.iter().zip(vals) {
                let j = j as usize;
                if j >= start {
                    let c = loc[j];
                    if c == usize::MAX {
                        return Err(SolverError::InvalidArgument(
                            "matrix pattern is not contained in the analysed pattern",
                        ));
                    }
                    f[c * m + r] += v;
                }
            }
            let (rows_, vals_) = apt.row(i);
            for (&j, &v) in rows_.iter().zip(vals_) {
                let j = j as usize;
                if j >= end {
                    let c = loc[j];
                    if c == usize::MAX {
                        return Err(SolverError::InvalidArgument(
                            "matrix pattern is not contained in the analysed pattern",
                        ));
                    }
                    f[r * m + c] += v;
                }
            }
        }
        // extend-add of the children's updates
        for &c in &node.children {
            if let Some(upd) = updates[c].take() {
                let sc = &sym.nodes[c].strct;
                let uc = sc.len();
                for b in 0..uc {
                    let lb = loc[sc[b]];
                    for a_ in 0..uc {
                        f[lb * m + loc[sc[a_]]] += upd[b * uc + a_];
                    }
                }
            }
        }
        for i in 0..p {
            loc[start + i] = usize::MAX;
        }
        for &g in &node.strct {
            loc[g] = usize::MAX;
        }
        // blocks of the front
        let mut f11 = Vec::with_capacity(p * p);
        let mut f21 = Vec::with_capacity(u * p);
        for c in 0..p {
            f11.extend_from_slice(&f[c * m..c * m + p]);
            f21.extend_from_slice(&f[c * m + p..(c + 1) * m]);
        }
        let mut f12 = Vec::with_capacity(p * u);
        let mut f22 = Vec::with_capacity(u * u);
        for c in 0..u {
            f12.extend_from_slice(&f[(p + c) * m..(p + c) * m + p]);
            f22.extend_from_slice(&f[(p + c) * m + p..(p + c + 1) * m]);
        }
        drop(f);
        let perm = dense_lu(&mut f11, p).ok_or(SolverError::Singular)?;
        if u > 0 {
            // rows of F₁₂ in the pivoted order (one column scratch at a time, avoiding a full p×u clone),
            // then U₁₂ = L₁₁⁻¹ P F₁₂, L₂₁ = F₂₁ U₁₁⁻¹
            let mut col_tmp = vec![T::ZERO; p];
            for c in 0..u {
                let col = &mut f12[c * p..(c + 1) * p];
                col_tmp.copy_from_slice(col);
                for i in 0..p {
                    col[i] = col_tmp[perm[i]];
                }
            }
            drop(col_tmp);
            solve_unit_lower(&f11, p, &mut f12, u);
            solve_upper_right(&f11, p, &mut f21, u);
            gemm_sub(&mut f22, &f21, &f12, u, u, p);
            updates[t] = Some(f22);
        }
        let compress = opts.tol > 0.0 && p >= opts.min_front;
        let (piv, l21, u12) = if compress {
            let piv = PivotLu::blocked(&f11, p, opts, &mut st);
            drop(f11);
            let l21 = compress_panel(&f21, u, p, opts, &mut st);
            drop(f21);
            let u12 = compress_panel(&f12, p, u, opts, &mut st);
            drop(f12);
            (piv, l21, u12)
        } else {
            (PivotLu::Dense(f11), Panel::dense(f21, u, p), Panel::dense(f12, p, u))
        };
        st.stored_entries += piv.stored() + l21.stored() + u12.stored();
        if p * p + 2 * p * u > st.largest_front.0 * st.largest_front.0 + 2 * st.largest_front.0 * st.largest_front.1 {
            st.largest_front = (p, u);
        }
        fronts.push(Front { p, piv, perm, l21, u12 });
    }
    st.factor_secs = t0.elapsed().as_secs_f64();
    Ok(BlrFactors { fronts, stats: st })
}

/// A low-rank multifrontal LU factorisation of one sparse matrix.
pub struct BlrLu<T: Scalar> {
    sym: BlrSymbolic,
    f: BlrFactors<T>,
}

impl<T: FaerScalar> BlrLu<T> {
    /// Analyses and factorises `a` (structurally symmetric).
    pub fn factor(a: &Csr<T>, opts: &BlrOptions) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let sym = BlrSymbolic::analyse(&[a], opts.leaf)?;
        let analyse = t0.elapsed().as_secs_f64();
        let mut f = factor_numeric(&sym, a, opts)?;
        f.stats.analyse_secs = analyse;
        Ok(Self { sym, f })
    }

    pub fn stats(&self) -> &BlrStats {
        &self.f.stats
    }

    pub fn symbolic(&self) -> &BlrSymbolic {
        &self.sym
    }

    pub fn dim(&self) -> usize {
        self.sym.n
    }

    /// Solves `A x = b` in place.
    pub fn solve_in_place(&self, b: &mut [T]) {
        solve_with(&self.sym, &self.f, b)
    }

    pub fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        let n = self.sym.n;
        for c in 0..nrhs {
            self.solve_in_place(&mut rhs[c * n..(c + 1) * n]);
        }
    }

    /// Solves `A x = b` to `rtol` by iterative refinement with this (approximate) factorisation: `x += M⁻¹(b − A x)`.
    /// Returns the iterations used and the final relative residual. `a` must be the matrix that was factorised.
    pub fn solve_refined(&self, a: &Csr<T>, b: &[T], x: &mut [T], max_it: usize, rtol: f64) -> (usize, f64) {
        let n = self.sym.n;
        x.copy_from_slice(b);
        self.solve_in_place(x);
        let bn = norm2(b).max(f64::MIN_POSITIVE);
        let mut r = vec![T::ZERO; n];
        let mut best = f64::INFINITY;
        for it in 0..=max_it {
            a.spmv(x, &mut r);
            for (ri, bi) in r.iter_mut().zip(b) {
                *ri = *bi - *ri;
            }
            let rel = norm2(&r) / bn;
            if rel <= rtol || it == max_it || rel > 0.9 * best {
                return (it, rel);
            }
            best = rel;
            self.solve_in_place(&mut r);
            axpy(T::ONE, &r, x);
        }
        unreachable!()
    }
}

fn solve_with<T: Scalar>(sym: &BlrSymbolic, f: &BlrFactors<T>, b: &mut [T]) {
    let n = sym.n;
    assert_eq!(b.len(), n);
    let mut w: Vec<T> = sym.order.iter().map(|&o| b[o]).collect();
    // forward
    for (node, fr) in sym.nodes.iter().zip(&f.fronts) {
        let (start, p) = (node.start, fr.p);
        let mut y: Vec<T> = (0..p).map(|i| w[start + fr.perm[i]]).collect();
        fr.piv.forward(p, &mut y);
        if !node.strct.is_empty() {
            let mut upd = vec![T::ZERO; node.strct.len()];
            fr.l21.apply_sub(&y, &mut upd); // upd = −L₂₁ y
            for (&g, &d) in node.strct.iter().zip(&upd) {
                w[g] += d;
            }
        }
        w[start..start + p].copy_from_slice(&y);
    }
    // backward
    for (node, fr) in sym.nodes.iter().zip(&f.fronts).rev() {
        let (start, p) = (node.start, fr.p);
        let mut y: Vec<T> = w[start..start + p].to_vec();
        if !node.strct.is_empty() {
            let xs: Vec<T> = node.strct.iter().map(|&g| w[g]).collect();
            fr.u12.apply_sub(&xs, &mut y);
        }
        fr.piv.backward(p, &mut y);
        w[start..start + p].copy_from_slice(&y);
    }
    for (k, &o) in sym.order.iter().enumerate() {
        b[o] = w[k];
    }
}

impl<T: FaerScalar> Preconditioner<T> for BlrLu<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        self.solve_in_place(z);
    }
}

/// [`ShiftedSolver`] on top of the BLR factorisation: `(A − σB)⁻¹`, the ordering, tree and front structures computed
/// once from the union pattern of `A` and `B`; every shift repeats only the numeric phase.
pub struct ShiftedBlrLu<T: FaerScalar> {
    a: Csr<T>,
    b: Option<Csr<T>>,
    sym: BlrSymbolic,
    opts: BlrOptions,
    factors: Option<BlrFactors<T>>,
    stats: FactorStats,
    solves: AtomicUsize,
}

impl<T: FaerScalar> ShiftedBlrLu<T> {
    pub fn new(a: &Csr<T>, b: Option<&Csr<T>>, opts: &BlrOptions) -> Result<Self, SolverError> {
        let t0 = Instant::now();
        let mut pats = vec![a];
        if let Some(b) = b {
            pats.push(b);
        }
        let sym = BlrSymbolic::analyse(&pats, opts.leaf)?;
        Ok(Self {
            a: a.clone(),
            b: b.cloned(),
            sym,
            opts: opts.clone(),
            factors: None,
            stats: FactorStats {
                analyse_secs: t0.elapsed().as_secs_f64(),
                symbolic_count: 1,
                numeric_secs_total: 0.0,
                numeric_secs_last: 0.0,
                numeric_count: 0,
                assemble_secs_last: 0.0,
                used: Used::None,
                backend: FactorBackend::None,
                cholesky_fallbacks: 0,
                solves: 0,
                factor_nnz: 0,
                solve_secs_total: 0.0,
            },
            solves: AtomicUsize::new(0),
        })
    }

    /// Replaces the numeric matrices while reusing the cached symbolic nested-dissection tree
    /// (`symbolic_count` stays 1). The new sparsity pattern must be contained in the analysed pattern;
    /// call [`ShiftedSolver::set_shift`] afterwards before solving.
    pub fn replace_matrices(&mut self, a: &Csr<T>, b: Option<&Csr<T>>) -> Result<(), SolverError> {
        let n = self.sym.n;
        if a.nrows() != n || a.ncols() != n {
            return Err(SolverError::InvalidArgument("matrix does not match the symbolic analysis"));
        }
        if b.is_some_and(|m| m.nrows() != n || m.ncols() != n) {
            return Err(SolverError::InvalidArgument("B matrix does not match the symbolic analysis"));
        }
        self.factors = None;
        self.a = a.clone();
        self.b = b.cloned();
        Ok(())
    }

    /// Statistics of the current BLR factorisation (compression ratio, ranks, stored/exact entries).
    pub fn blr_stats(&self) -> Option<&BlrStats> {
        self.factors.as_ref().map(|f| &f.stats)
    }

    pub fn symbolic(&self) -> &BlrSymbolic {
        &self.sym
    }

    pub fn options(&self) -> &BlrOptions {
        &self.opts
    }

    /// Solves multiple right-hand-side columns against the current BLR factorisation.
    pub fn solve_many(&self, bs: &[Vec<T>]) -> Vec<Vec<T>> {
        let f = self.factors.as_ref().expect("set_shift must be called before solving");
        let out: Vec<Vec<T>> = bs
            .iter()
            .map(|b| {
                let mut x = b.clone();
                solve_with(&self.sym, f, &mut x);
                x
            })
            .collect();
        self.solves.fetch_add(bs.len(), Ordering::Relaxed);
        out
    }

    /// Solves `a x = b` to `rtol` by deterministic iterative refinement using the current BLR factor.
    pub fn solve_refined(&self, a: &Csr<T>, b: &[T], x: &mut [T], max_it: usize, rtol: f64) -> (usize, f64) {
        let f = self.factors.as_ref().expect("set_shift must be called before solving");
        let n = self.sym.n;
        assert!(b.len() == n && x.len() == n);
        x.copy_from_slice(b);
        solve_with(&self.sym, f, x);
        self.solves.fetch_add(1, Ordering::Relaxed);
        let bn = norm2(b).max(f64::MIN_POSITIVE);
        let mut r = vec![T::ZERO; n];
        let mut best = f64::INFINITY;
        for it in 0..=max_it {
            a.spmv(x, &mut r);
            for (ri, bi) in r.iter_mut().zip(b) {
                *ri = *bi - *ri;
            }
            let rel = norm2(&r) / bn;
            if rel <= rtol || it == max_it || rel > 0.9 * best {
                return (it, rel);
            }
            best = rel;
            solve_with(&self.sym, f, &mut r);
            self.solves.fetch_add(1, Ordering::Relaxed);
            axpy(T::ONE, &r, x);
        }
        unreachable!()
    }
}

impl<T: FaerScalar> ShiftedSolver<T> for ShiftedBlrLu<T> {
    fn dim(&self) -> usize {
        self.sym.n
    }

    fn set_shift(&mut self, sigma: T) -> Result<(), SolverError> {
        // Release previous numeric factors before building the new shift so peak RAM never holds two shifts.
        self.factors = None;
        let t0 = Instant::now();
        let owned_shifted;
        let mat = if sigma == T::ZERO {
            &self.a
        } else {
            owned_shifted = match &self.b {
                Some(b) => Csr::axpby(T::ONE, &self.a, -sigma, b),
                None => Csr::axpby(T::ONE, &self.a, -sigma, &Csr::identity(self.sym.n)),
            };
            &owned_shifted
        };
        self.stats.assemble_secs_last = t0.elapsed().as_secs_f64();
        let t1 = Instant::now();
        let f = factor_numeric(&self.sym, mat, &self.opts)?;
        self.stats.numeric_secs_last = t1.elapsed().as_secs_f64();
        self.stats.numeric_secs_total += self.stats.numeric_secs_last;
        self.stats.numeric_count += 1;
        self.stats.used = Used::Lu;
        self.stats.backend = FactorBackend::BlrNd;
        self.stats.factor_nnz = f.stats.stored_entries;
        self.factors = Some(f);
        Ok(())
    }

    fn solve_block_in_place(&self, rhs: &mut [T], nrhs: usize) {
        let f = self.factors.as_ref().expect("set_shift must be called before solving");
        let n = self.sym.n;
        for c in 0..nrhs {
            solve_with(&self.sym, f, &mut rhs[c * n..(c + 1) * n]);
        }
        self.solves.fetch_add(nrhs, Ordering::Relaxed);
    }

    fn stats(&self) -> FactorStats {
        FactorStats { solves: self.solves.load(Ordering::Relaxed), ..self.stats }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::maxwell::nedelec_cube;
    use crate::linalg::rng::random_vector;
    use crate::linalg::scalar::C64;

    #[test]
    fn shifted_blr_lu_complex_symmetric_nedelec_3d() {
        let cube = nedelec_cube(4);
        let k = cube.stiffness.to_complex();
        let m = cube.mass.to_complex();
        // Complex-symmetric driven operator A = K + i k0 C - k0^2 M (A = A^T != A^H)
        let k0 = 3.5f64;
        let a_loss = Csr::axpby(C64::ONE, &k, C64::new(0.0, 0.15 * k0), &m);
        let sigma = C64::new(k0 * k0, 0.0);
        let a_exact = Csr::axpby(C64::ONE, &a_loss, -sigma, &m);
        assert!(a_exact.is_symmetric(1e-12));
        assert!(!a_exact.is_hermitian(1e-6));

        let opts = BlrOptions { tol: 1e-6, leaf: 32, block: 32, min_front: 48, accept_ratio: 0.7 };
        let mut blr = ShiftedBlrLu::<C64>::new(&a_loss, Some(&m), &opts).unwrap();
        blr.set_shift(sigma).unwrap();

        assert_eq!(blr.stats().backend, FactorBackend::BlrNd);
        assert_eq!(blr.stats().symbolic_count, 1);
        assert_eq!(blr.stats().numeric_count, 1);
        let bs = blr.blr_stats().expect("blr_stats must be populated after set_shift");
        assert!(bs.exact_entries > 0 && bs.stored_entries > 0);
        assert!(bs.compression() <= 1.0001);

        let n = a_exact.nrows();
        let rhs0 = random_vector::<C64>(n, 101);
        let rhs1 = random_vector::<C64>(n, 202);
        let sols = blr.solve_many(&[rhs0.clone(), rhs1.clone()]);
        for (x, b) in sols.iter().zip([&rhs0, &rhs1]) {
            let mut ax = vec![C64::ZERO; n];
            a_exact.spmv(x, &mut ax);
            let err = ax.iter().zip(b).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
            let rel = err / norm2(b).max(1e-300);
            assert!(rel < 1e-4, "BLR direct solve relative residual too large: {rel:e}");
        }

        let mut x_ref = vec![C64::ZERO; n];
        let (_it, rel_ref) = blr.solve_refined(&a_exact, &rhs0, &mut x_ref, 5, 1e-12);
        assert!(rel_ref < 1e-11, "BLR refined residual {rel_ref:e}");

        // Frequency update via replace_matrices reuses symbolic tree (`symbolic_count == 1`)
        let k0_2 = 4.2f64;
        let a_2 = Csr::axpby(C64::ONE, &k, C64::new(-k0_2 * k0_2, 0.15 * k0_2), &m);
        blr.replace_matrices(&a_2, None).unwrap();
        blr.set_shift(C64::ZERO).unwrap();
        assert_eq!(blr.stats().symbolic_count, 1);
        assert_eq!(blr.stats().numeric_count, 2);
        let mut x2 = rhs0.clone();
        blr.solve_block_in_place(&mut x2, 1);
        let mut ax2 = vec![C64::ZERO; n];
        a_2.spmv(&x2, &mut ax2);
        let rel2 = ax2.iter().zip(&rhs0).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt() / norm2(&rhs0);
        assert!(rel2 < 1e-4, "BLR updated frequency residual {rel2:e}");
    }
}
