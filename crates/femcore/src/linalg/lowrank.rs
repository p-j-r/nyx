//! **Low-rank (rank-structured) matrices and a fast direct solver** (port item S4): adaptive cross
//! approximation, the HODLR format (hierarchically off-diagonal low-rank) and its `O(n log² n)` direct solve.
//!
//! # Why this exists
//! After nested dissection the large dense blocks of a sparse factorisation sit on the *separators*, and
//! the dense Schur complements there have off-diagonal blocks of very low numerical rank (the interaction
//! of two parts of a smooth interface is smooth). Block low-rank (BLR) and hierarchical (HSS/HODLR) solvers
//! exploit that to cut the memory and the flops of the dense fronts. The same structure appears in boundary
//! element matrices, covariance matrices, integral operators and in the inverses of banded operators.
//!
//! # What is here
//! * [`LowRank`]: `A ≈ U Vᴴ`; [`aca_partial`] builds it from a few rows and columns of the block (it never
//!   forms the block: `O(r(m+n))` entry evaluations), [`aca_full`] from an explicit dense block (pivot =
//!   largest remaining entry; the robust choice when the block exists anyway). Stopping rule: the newest
//!   rank-1 term is below `tol` times the running estimate of `‖A‖_F` (partial pivoting additionally requires
//!   two consecutive small terms, because a single small term can be a pivot-row accident).
//! * [`Hodlr`]: balanced binary tree over the index range; diagonal blocks dense at the leaves, every
//!   off-diagonal block `A₁₂`, `A₂₁` of every node compressed. Built from an entry function
//!   ([`Hodlr::from_entries`]), a dense matrix or a sparse one. Fast `matvec` (it is a [`LinearOperator`]).
//! * [`HodlrLu`] (from [`Hodlr::factor`]): the recursive *Woodbury* direct solver of Ambikasaran & Darve (2013;
//!   Kressner–Massei–Robol 2019 for the same idea in matrix-function form). With `A = D + U_c V_cᴴ`,
//!   `D = diag(A₁₁, A₂₂)`, `U_c = diag(P₁, P₂)`, `V_c = [0 Q₁; Q₂ 0]` (so `A₁₂ = P₁Q₂ᴴ`, `A₂₁ = P₂Q₁ᴴ`),
//!   `A⁻¹ = D⁻¹ − D⁻¹U_c (I + V_cᴴD⁻¹U_c)⁻¹ V_cᴴD⁻¹`; `D⁻¹` is the same solver one level down, applied to the
//!   `r` columns of `U_c` once at factorisation time. No pivoting: it is intended for (block) diagonally
//!   dominant, SPD or otherwise well-conditioned-leaf matrices, and as a *preconditioner* when loose tolerances are used
//!   ([`HodlrLu`] is a [`Preconditioner`]).
//!
//! # What is not here
//! The sparse use case itself — a multifrontal factorisation with low-rank fronts (MUMPS-BLR, STRUMPACK) —
//! needs access to the elimination tree and the frontal matrices, which neither `faer`'s LU nor
//! [`crate::linalg::sparse_lu`] exposes. [`Hodlr::from_csr`] with a dense Schur complement (see the tests) shows
//! the compression that such a solver would exploit; wiring it into a supernodal sparse factorisation
//! is left open (`docs/STATUS.md`). HSS (nested bases, `O(n)`) and H² formats are not implemented: HODLR is the
//! simplest member of the family and its cost is only a log factor worse.

use crate::linalg::csr::Csr;
use crate::linalg::linop::{LinearOperator, Preconditioner};
use crate::linalg::scalar::Scalar;
use crate::linalg::vecops::{axpy, dotc, norm2};

/// `A ≈ U Vᴴ` with `U: m × r`, `V: n × r` (column-major).
#[derive(Debug, Clone)]
pub struct LowRank<T: Scalar> {
    pub m: usize,
    pub n: usize,
    pub r: usize,
    pub u: Vec<T>,
    pub v: Vec<T>,
}

impl<T: Scalar> LowRank<T> {
    pub fn zero(m: usize, n: usize) -> Self {
        Self { m, n, r: 0, u: Vec::new(), v: Vec::new() }
    }

    fn from_terms(m: usize, n: usize, us: Vec<Vec<T>>, vs: Vec<Vec<T>>) -> Self {
        let r = us.len();
        Self { m, n, r, u: us.concat(), v: vs.concat() }
    }

    pub fn storage(&self) -> usize {
        self.r * (self.m + self.n)
    }

    pub fn ucol(&self, k: usize) -> &[T] {
        &self.u[k * self.m..(k + 1) * self.m]
    }
    pub fn vcol(&self, k: usize) -> &[T] {
        &self.v[k * self.n..(k + 1) * self.n]
    }

    /// `y += U (Vᴴ x)`.
    pub fn apply_acc(&self, x: &[T], y: &mut [T]) {
        for k in 0..self.r {
            let c = dotc(self.vcol(k), x);
            axpy(c, self.ucol(k), y);
        }
    }

    /// `Vᴴ x`.
    pub fn vh_apply(&self, x: &[T]) -> Vec<T> {
        (0..self.r).map(|k| dotc(self.vcol(k), x)).collect()
    }

    /// The block as a dense column-major array (testing and error measurement).
    pub fn to_dense(&self) -> Vec<T> {
        let mut d = vec![T::ZERO; self.m * self.n];
        for k in 0..self.r {
            for j in 0..self.n {
                let c = self.vcol(k)[j].conj();
                for i in 0..self.m {
                    d[j * self.m + i] += self.ucol(k)[i] * c;
                }
            }
        }
        d
    }
}

/// Adaptive cross approximation with **partial pivoting** (Bebendorf): the block `A[0..m, 0..n]` is given
/// only through `entry(i, j)`. Returns the factors with relative Frobenius accuracy about `tol` (heuristic, as
/// for every ACA), at most `max_rank` terms.
pub fn aca_partial<T: Scalar>(
    m: usize,
    n: usize,
    entry: &dyn Fn(usize, usize) -> T,
    tol: f64,
    max_rank: usize,
) -> LowRank<T> {
    if m == 0 || n == 0 {
        return LowRank::zero(m, n);
    }
    let max_rank = max_rank.min(m.min(n));
    let mut us: Vec<Vec<T>> = Vec::new();
    let mut vs: Vec<Vec<T>> = Vec::new(); // plain (not conjugated) rows: A ≈ Σ u_k v_kᵀ
    let mut row_used = vec![false; m];
    let mut frob2 = 0.0f64;
    let mut small_in_a_row = 0usize;
    let mut i = 0usize;
    let mut zero_rows = 0usize;
    while us.len() < max_rank {
        row_used[i] = true;
        // residual row i
        let mut row: Vec<T> = (0..n).map(|j| entry(i, j)).collect();
        for (u, v) in us.iter().zip(&vs) {
            let ui = u[i];
            if ui != T::ZERO {
                for (rj, vj) in row.iter_mut().zip(v) {
                    *rj -= ui * *vj;
                }
            }
        }
        let (jmax, amax) =
            row.iter().enumerate().map(|(j, v)| (j, v.modulus())).fold((0, 0.0), |b, c| if c.1 > b.1 { c } else { b });
        if amax <= 1e-300 {
            // a zero residual row: try another one (give up when every row was tried and nothing was found)
            zero_rows += 1;
            match (0..m).find(|&q| !row_used[q]) {
                Some(next) if zero_rows <= m => {
                    i = next;
                    continue;
                }
                _ => break,
            }
        }
        let piv = row[jmax];
        let v: Vec<T> = row.iter().map(|x| *x / piv).collect();
        let mut u: Vec<T> = (0..m).map(|q| entry(q, jmax)).collect();
        for (uk, vk) in us.iter().zip(&vs) {
            let c = vk[jmax];
            if c != T::ZERO {
                for (ui, uki) in u.iter_mut().zip(uk) {
                    *ui -= *uki * c;
                }
            }
        }
        // ‖Σ u_k v_kᵀ‖_F² update:  ‖u‖²‖v‖² + 2 Re Σ_l (u_lᴴ u)(v_lᴴ v)‾ …
        let (nu, nv) = (norm2(&u), norm2(&v));
        let mut cross = 0.0;
        for (ul, vl) in us.iter().zip(&vs) {
            let a = dotc(ul, &u);
            cross += (a * dotc(vl, &v)).re();
        }
        frob2 += nu * nu * nv * nv + 2.0 * cross;
        let term = nu * nv;
        us.push(u.clone());
        vs.push(v);
        if term <= tol * frob2.max(0.0).sqrt() {
            small_in_a_row += 1;
            if small_in_a_row >= 2 {
                break;
            }
        } else {
            small_in_a_row = 0;
        }
        // next pivot row: the largest entry of the new column among unused rows
        let next = (0..m).filter(|&q| !row_used[q]).max_by(|&a, &b| u[a].modulus().total_cmp(&u[b].modulus()));
        match next {
            Some(q) => i = q,
            None => break,
        }
    }
    let vconj: Vec<Vec<T>> = vs.into_iter().map(|v| v.into_iter().map(|x| x.conj()).collect()).collect();
    LowRank::from_terms(m, n, us, vconj)
}

/// ACA with **full pivoting** on an explicit column-major `m × n` block (destroyed; on return it holds the residual):
/// at each step the largest remaining entry is the pivot, and the rank-1 update, the remaining Frobenius norm and the
/// next pivot are computed in one pass. Stops when the remaining norm is below `tol · ‖A‖_F` or at `max_rank`.
pub fn aca_full<T: Scalar>(a: &mut [T], m: usize, n: usize, tol: f64, max_rank: usize) -> LowRank<T> {
    assert_eq!(a.len(), m * n);
    if m == 0 || n == 0 {
        return LowRank::zero(m, n);
    }
    let max_rank = max_rank.min(m.min(n));
    // one pass gives the Frobenius norm and the position of the largest entry; every later pass fuses the rank-1
    // update with the recomputation of both
    let scan = |a: &[T]| -> (f64, usize, f64) {
        let (mut sum, mut best, mut at) = (0.0f64, 0.0f64, 0usize);
        for (k, v) in a.iter().enumerate() {
            let q = v.modulus_sq();
            sum += q;
            if q > best {
                best = q;
                at = k;
            }
        }
        (sum, at, best)
    };
    let (sum0, mut at, mut best) = scan(a);
    let norm0 = sum0.sqrt();
    if norm0 == 0.0 {
        return LowRank::zero(m, n);
    }
    let mut rem2 = sum0;
    let mut us: Vec<Vec<T>> = Vec::new();
    let mut vs: Vec<Vec<T>> = Vec::new();
    while us.len() < max_rank {
        if rem2.max(0.0).sqrt() <= tol * norm0 || best == 0.0 {
            break;
        }
        let (bi, bj) = (at % m, at / m);
        let piv = a[bj * m + bi];
        let u: Vec<T> = (0..m).map(|i| a[bj * m + i]).collect();
        let v: Vec<T> = (0..n).map(|j| a[j * m + bi] / piv).collect();
        let (mut sum, mut nbest, mut nat) = (0.0f64, 0.0f64, 0usize);
        for j in 0..n {
            let vj = v[j];
            let col = &mut a[j * m..(j + 1) * m];
            for (i, (x, ui)) in col.iter_mut().zip(&u).enumerate() {
                *x -= *ui * vj;
                let q = x.modulus_sq();
                sum += q;
                if q > nbest {
                    nbest = q;
                    nat = j * m + i;
                }
            }
        }
        rem2 = sum;
        best = nbest;
        at = nat;
        us.push(u);
        vs.push(v.into_iter().map(|x| x.conj()).collect());
    }
    LowRank::from_terms(m, n, us, vs)
}

// ---------------------------------------------------------------------------------------------------------

/// Dense LU with partial pivoting (the leaves and the small Woodbury systems).
#[derive(Debug, Clone)]
struct DenseLu<T: Scalar> {
    n: usize,
    lu: Vec<T>,
    piv: Vec<usize>,
}

impl<T: Scalar> DenseLu<T> {
    fn new(mut a: Vec<T>, n: usize) -> Option<Self> {
        let mut piv = vec![0usize; n];
        for k in 0..n {
            let mut p = k;
            let mut best = a[k * n + k].modulus();
            for i in k + 1..n {
                let v = a[k * n + i].modulus();
                if v > best {
                    best = v;
                    p = i;
                }
            }
            if best == 0.0 || !best.is_finite() {
                return None;
            }
            piv[k] = p;
            if p != k {
                for j in 0..n {
                    a.swap(j * n + k, j * n + p);
                }
            }
            let d = a[k * n + k];
            for i in k + 1..n {
                a[k * n + i] = a[k * n + i] / d;
            }
            for j in k + 1..n {
                let c = a[j * n + k];
                if c != T::ZERO {
                    for i in k + 1..n {
                        let l = a[k * n + i];
                        a[j * n + i] -= l * c;
                    }
                }
            }
        }
        Some(Self { n, lu: a, piv })
    }

    fn solve(&self, b: &mut [T]) {
        let n = self.n;
        for k in 0..n {
            b.swap(k, self.piv[k]);
        }
        for k in 0..n {
            let bk = b[k];
            for (bi, l) in b[k + 1..n].iter_mut().zip(&self.lu[k * n + k + 1..(k + 1) * n]) {
                *bi -= *l * bk;
            }
        }
        for k in (0..n).rev() {
            b[k] = b[k] / self.lu[k * n + k];
            let bk = b[k];
            for (bi, u) in b[..k].iter_mut().zip(&self.lu[k * n..k * n + k]) {
                *bi -= *u * bk;
            }
        }
    }
}

#[derive(Debug, Clone)]
enum Node<T: Scalar> {
    Leaf(Vec<T>), // column-major n×n
    Inner { n1: usize, left: Box<Node<T>>, right: Box<Node<T>>, a12: LowRank<T>, a21: LowRank<T> },
}

/// Options of the HODLR construction.
#[derive(Debug, Clone, Copy)]
pub struct HodlrOptions {
    /// Diagonal blocks of at most this size are stored densely.
    pub leaf: usize,
    /// Relative Frobenius tolerance of each off-diagonal block.
    pub tol: f64,
    pub max_rank: usize,
}

impl Default for HodlrOptions {
    fn default() -> Self {
        Self { leaf: 64, tol: 1e-8, max_rank: 200 }
    }
}

/// A square matrix in HODLR format.
#[derive(Debug, Clone)]
pub struct Hodlr<T: Scalar> {
    n: usize,
    root: Node<T>,
}

impl<T: Scalar> Hodlr<T> {
    /// Builds the representation from `entry(i, j)` without forming the matrix (off-diagonal blocks by
    /// [`aca_partial`], diagonal leaves evaluated explicitly).
    pub fn from_entries(n: usize, entry: &dyn Fn(usize, usize) -> T, opts: HodlrOptions) -> Self {
        fn build<T: Scalar>(lo: usize, hi: usize, entry: &dyn Fn(usize, usize) -> T, o: HodlrOptions) -> Node<T> {
            let n = hi - lo;
            if n <= o.leaf.max(2) {
                let mut a = vec![T::ZERO; n * n];
                for j in 0..n {
                    for i in 0..n {
                        a[j * n + i] = entry(lo + i, lo + j);
                    }
                }
                return Node::Leaf(a);
            }
            let n1 = n / 2;
            let mid = lo + n1;
            let a12 = aca_partial(n1, n - n1, &|i, j| entry(lo + i, mid + j), o.tol, o.max_rank);
            let a21 = aca_partial(n - n1, n1, &|i, j| entry(mid + i, lo + j), o.tol, o.max_rank);
            Node::Inner {
                n1,
                left: Box::new(build(lo, mid, entry, o)),
                right: Box::new(build(mid, hi, entry, o)),
                a12,
                a21,
            }
        }
        Self { n, root: build(0, n, entry, opts) }
    }

    /// From a column-major dense matrix; off-diagonal blocks by [`aca_full`].
    pub fn from_dense(a: &[T], n: usize, opts: HodlrOptions) -> Self {
        assert_eq!(a.len(), n * n);
        fn build<T: Scalar>(a: &[T], ld: usize, lo: usize, hi: usize, o: HodlrOptions) -> Node<T> {
            let n = hi - lo;
            if n <= o.leaf.max(2) {
                let mut d = vec![T::ZERO; n * n];
                for j in 0..n {
                    for i in 0..n {
                        d[j * n + i] = a[(lo + j) * ld + lo + i];
                    }
                }
                return Node::Leaf(d);
            }
            let n1 = n / 2;
            let mid = lo + n1;
            let block = |r0: usize, r1: usize, c0: usize, c1: usize| -> Vec<T> {
                let mut b = Vec::with_capacity((r1 - r0) * (c1 - c0));
                for j in c0..c1 {
                    b.extend_from_slice(&a[j * ld + r0..j * ld + r1]);
                }
                b
            };
            let mut b12 = block(lo, mid, mid, hi);
            let mut b21 = block(mid, hi, lo, mid);
            Node::Inner {
                n1,
                left: Box::new(build(a, ld, lo, mid, o)),
                right: Box::new(build(a, ld, mid, hi, o)),
                a12: aca_full(&mut b12, n1, n - n1, o.tol, o.max_rank),
                a21: aca_full(&mut b21, n - n1, n1, o.tol, o.max_rank),
            }
        }
        Self { n, root: build(a, n, 0, n, opts) }
    }

    /// From a sparse matrix (off-diagonal blocks are extracted densely, so memory is `O(n²/4)` transiently:
    /// meant for matrices up to a few thousand rows, such as dense Schur complements stored sparsely).
    pub fn from_csr(a: &Csr<T>, opts: HodlrOptions) -> Self {
        let n = a.nrows();
        assert_eq!(a.ncols(), n);
        let dense = a.to_dense(); // row-major
        let mut cm = vec![T::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                cm[j * n + i] = dense[i * n + j];
            }
        }
        Self::from_dense(&cm, n, opts)
    }

    pub fn dim(&self) -> usize {
        self.n
    }

    /// Scalars stored (dense leaves + low-rank factors); compare with `n²`.
    pub fn storage(&self) -> usize {
        fn go<T: Scalar>(n: &Node<T>) -> usize {
            match n {
                Node::Leaf(a) => a.len(),
                Node::Inner { left, right, a12, a21, .. } => go(left) + go(right) + a12.storage() + a21.storage(),
            }
        }
        go(&self.root)
    }

    /// Largest numerical rank of the off-diagonal blocks on each level (root first).
    pub fn rank_profile(&self) -> Vec<usize> {
        fn go<T: Scalar>(n: &Node<T>, level: usize, out: &mut Vec<usize>) {
            if let Node::Inner { left, right, a12, a21, .. } = n {
                if out.len() <= level {
                    out.resize(level + 1, 0);
                }
                out[level] = out[level].max(a12.r).max(a21.r);
                go(left, level + 1, out);
                go(right, level + 1, out);
            }
        }
        let mut out = Vec::new();
        go(&self.root, 0, &mut out);
        out
    }

    /// `y = A x`.
    pub fn matvec(&self, x: &[T], y: &mut [T]) {
        assert!(x.len() == self.n && y.len() == self.n);
        fn go<T: Scalar>(node: &Node<T>, x: &[T], y: &mut [T]) {
            match node {
                Node::Leaf(a) => {
                    let n = x.len();
                    y.fill(T::ZERO);
                    for j in 0..n {
                        axpy(x[j], &a[j * n..(j + 1) * n], y);
                    }
                }
                Node::Inner { n1, left, right, a12, a21 } => {
                    let (x1, x2) = x.split_at(*n1);
                    let (y1, y2) = y.split_at_mut(*n1);
                    go(left, x1, y1);
                    go(right, x2, y2);
                    a12.apply_acc(x2, y1);
                    a21.apply_acc(x1, y2);
                }
            }
        }
        go(&self.root, x, y)
    }

    /// Factorises the matrix (recursive Woodbury, see the module documentation). `None` if a leaf block or a
    /// small capacitance system is singular.
    pub fn factor(&self) -> Option<HodlrLu<T>> {
        fn go<T: Scalar>(node: &Node<T>) -> Option<FNode<T>> {
            match node {
                Node::Leaf(a) => {
                    let n = (a.len() as f64).sqrt().round() as usize;
                    Some(FNode::Leaf(DenseLu::new(a.clone(), n)?))
                }
                Node::Inner { n1, left, right, a12, a21 } => {
                    let (l, r) = (go(left)?, go(right)?);
                    let n2 = a12.n;
                    let (r1, r2) = (a12.r, a21.r);
                    // W1 = A11⁻¹ P1  (n1 × r1),  W2 = A22⁻¹ P2  (n2 × r2)
                    let mut w1 = a12.u.clone();
                    for c in 0..r1 {
                        l.solve(&mut w1[c * n1..(c + 1) * n1]);
                    }
                    let mut w2 = a21.u.clone();
                    for c in 0..r2 {
                        r.solve(&mut w2[c * n2..(c + 1) * n2]);
                    }
                    // S = [[I, Q2ᴴ W2],[Q1ᴴ W1, I]]
                    let rr = r1 + r2;
                    let mut s = vec![T::ZERO; rr * rr];
                    for k in 0..rr {
                        s[k * rr + k] = T::ONE;
                    }
                    for c in 0..r2 {
                        for a in 0..r1 {
                            s[(r1 + c) * rr + a] = dotc(a12.vcol(a), &w2[c * n2..(c + 1) * n2]);
                        }
                    }
                    for c in 0..r1 {
                        for a in 0..r2 {
                            s[c * rr + r1 + a] = dotc(a21.vcol(a), &w1[c * n1..(c + 1) * n1]);
                        }
                    }
                    let s_lu = if rr > 0 { Some(DenseLu::new(s, rr)?) } else { None };
                    Some(FNode::Inner {
                        n1: *n1,
                        left: Box::new(l),
                        right: Box::new(r),
                        w1,
                        w2,
                        q2: a12.v.clone(),
                        q1: a21.v.clone(),
                        r1,
                        r2,
                        s: s_lu,
                    })
                }
            }
        }
        Some(HodlrLu { n: self.n, root: go(&self.root)? })
    }
}

impl<T: Scalar> LinearOperator<T> for Hodlr<T> {
    fn dim(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[T], y: &mut [T]) {
        self.matvec(x, y)
    }
}

enum FNode<T: Scalar> {
    Leaf(DenseLu<T>),
    Inner {
        n1: usize,
        left: Box<FNode<T>>,
        right: Box<FNode<T>>,
        w1: Vec<T>,
        w2: Vec<T>,
        q2: Vec<T>,
        q1: Vec<T>,
        r1: usize,
        r2: usize,
        s: Option<DenseLu<T>>,
    },
}

impl<T: Scalar> FNode<T> {
    fn solve(&self, b: &mut [T]) {
        match self {
            FNode::Leaf(lu) => lu.solve(b),
            FNode::Inner { n1, left, right, w1, w2, q2, q1, r1, r2, s } => {
                let n = b.len();
                let n2 = n - n1;
                let (b1, b2) = b.split_at_mut(*n1);
                left.solve(b1);
                right.solve(b2);
                let Some(s) = s else { return };
                // t = Vcᴴ z = [Q2ᴴ z2 ; Q1ᴴ z1]
                let mut t: Vec<T> = Vec::with_capacity(r1 + r2);
                for a in 0..*r1 {
                    t.push(dotc(&q2[a * n2..(a + 1) * n2], b2));
                }
                for a in 0..*r2 {
                    t.push(dotc(&q1[a * n1..(a + 1) * n1], b1));
                }
                s.solve(&mut t);
                for a in 0..*r1 {
                    axpy(-t[a], &w1[a * n1..(a + 1) * n1], b1);
                }
                for a in 0..*r2 {
                    axpy(-t[r1 + a], &w2[a * n2..(a + 1) * n2], b2);
                }
            }
        }
    }
}

/// Factorisation of a [`Hodlr`] matrix.
pub struct HodlrLu<T: Scalar> {
    n: usize,
    root: FNode<T>,
}

impl<T: Scalar> HodlrLu<T> {
    pub fn dim(&self) -> usize {
        self.n
    }
    pub fn solve_in_place(&self, b: &mut [T]) {
        assert_eq!(b.len(), self.n);
        self.root.solve(b)
    }
}

impl<T: Scalar> Preconditioner<T> for HodlrLu<T> {
    fn apply(&self, r: &[T], z: &mut [T]) {
        z.copy_from_slice(r);
        self.solve_in_place(z);
    }
}
