//! The one matrix type of the finite-element code: assembly (`Triplets`), the operations the kernels use, and the aliases that keep
//! kernel code readable (`RCsr` real, `CCsr` complex, `Coo`, `CCoo`). Before phase 3 the kernels had a real and a complex copy of all
//! of this (`sparse`, `csparse`). A child module of `csr` so that it can build `Csr` values from their private fields.

use super::*;

/// Real matrix.
pub type RCsr = Csr<f64>;
/// Complex matrix.
pub type CCsr = Csr<C64>;
/// Complex number as the kernels write it.
pub type C = C64;

/// `re + i im`.
pub fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

/// Coordinate-format builder; duplicate entries are summed when converted ([`Triplets::to_csr`]).
#[derive(Clone, Debug)]
pub struct Triplets<T: Scalar> {
    nrows: usize,
    ncols: usize,
    rows: Vec<u32>,
    cols: Vec<u32>,
    vals: Vec<T>,
}

pub type Coo = Triplets<f64>;
pub type CCoo = Triplets<C64>;

impl<T: Scalar> Triplets<T> {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        assert!(nrows <= u32::MAX as usize && ncols <= u32::MAX as usize, "Triplets dimensions ({nrows}, {ncols}) overflow u32");
        Triplets { nrows, ncols, rows: vec![], cols: vec![], vals: vec![] }
    }
    pub fn with_capacity(nrows: usize, ncols: usize, cap: usize) -> Self {
        assert!(nrows <= u32::MAX as usize && ncols <= u32::MAX as usize, "Triplets dimensions ({nrows}, {ncols}) overflow u32");
        Triplets { nrows, ncols, rows: Vec::with_capacity(cap), cols: Vec::with_capacity(cap), vals: Vec::with_capacity(cap) }
    }
    pub fn nrows(&self) -> usize {
        self.nrows
    }
    pub fn ncols(&self) -> usize {
        self.ncols
    }
    pub fn rows(&self) -> &[u32] {
        &self.rows
    }
    pub fn cols(&self) -> &[u32] {
        &self.cols
    }
    pub fn vals(&self) -> &[T] {
        &self.vals
    }
    #[inline]
    pub fn push(&mut self, r: usize, c: usize, v: T) {
        assert!(r < self.nrows && c < self.ncols, "Triplets::push ({r}, {c}) out of bounds ({}, {})", self.nrows, self.ncols);
        let r_u32 = u32::try_from(r).expect("Triplets row index overflows u32");
        let c_u32 = u32::try_from(c).expect("Triplets column index overflows u32");
        self.rows.push(r_u32);
        self.cols.push(c_u32);
        self.vals.push(v);
    }
    /// Appends all entries from `other` after checking matching `(nrows, ncols)` dimensions.
    pub fn extend(&mut self, other: &Self) {
        assert!(
            self.nrows == other.nrows && self.ncols == other.ncols,
            "Triplets::extend shape mismatch: ({}, {}) != ({}, {})",
            self.nrows,
            self.ncols,
            other.nrows,
            other.ncols
        );
        assert!(
            other.rows.len() == other.cols.len() && other.rows.len() == other.vals.len(),
            "Triplets::extend buffer length mismatch"
        );
        self.rows.extend_from_slice(&other.rows);
        self.cols.extend_from_slice(&other.cols);
        self.vals.extend_from_slice(&other.vals);
    }
    /// Adds the dense block `local[a*cols.len() + b]` at `(rows[a], cols[b])`; exact zeros are skipped.
    pub fn push_block(&mut self, rows: &[usize], cols: &[usize], local: &[T]) {
        let expected = rows.len().checked_mul(cols.len()).expect("push_block block size overflow");
        assert_eq!(local.len(), expected);
        for (a, &r) in rows.iter().enumerate() {
            for (b, &cc) in cols.iter().enumerate() {
                let v = local[a * cols.len() + b];
                if v != T::ZERO {
                    self.push(r, cc, v);
                }
            }
        }
    }
    /// Counting sort by row, then each row sorted by column with duplicates summed.
    pub fn to_csr(&self) -> Csr<T> {
        let n = self.nrows;
        let nnz = self.rows.len();
        assert!(nnz == self.cols.len() && nnz == self.vals.len(), "Triplets buffer length mismatch");
        let mut count = vec![0usize; n.checked_add(1).expect("Triplets nrows overflow")];
        for &r in &self.rows {
            let ru = r as usize;
            assert!(ru < n, "Triplets row index {ru} >= {n}");
            count[ru + 1] += 1;
        }
        for i in 0..n {
            count[i + 1] += count[i];
        }
        let mut order = vec![0usize; nnz];
        let mut next = count.clone();
        for k in 0..nnz {
            let r = self.rows[k] as usize;
            order[next[r]] = k;
            next[r] += 1;
        }
        let mut indptr = Vec::with_capacity(n + 1);
        let mut indices: Vec<u32> = Vec::with_capacity(nnz);
        let mut values: Vec<T> = Vec::with_capacity(nnz);
        indptr.push(0);
        let mut row: Vec<(u32, T)> = Vec::new();
        for r in 0..n {
            row.clear();
            for &k in &order[count[r]..count[r + 1]] {
                let col = self.cols[k];
                assert!((col as usize) < self.ncols, "Triplets col index {col} >= {}", self.ncols);
                row.push((col, self.vals[k]));
            }
            row.sort_unstable_by_key(|e| e.0);
            let mut i = 0;
            while i < row.len() {
                let col = row[i].0;
                let mut s = T::ZERO;
                while i < row.len() && row[i].0 == col {
                    s += row[i].1;
                    i += 1;
                }
                indices.push(col);
                values.push(s);
            }
            indptr.push(indices.len());
        }
        Csr::try_new(n, self.ncols, indptr, indices, values).expect("valid CSR from Triplets")
    }
}

pub(super) type SumBlock<T> = (Vec<usize>, Vec<u32>, Vec<T>);

pub(super) fn validate_sum_parts<T: Scalar>(parts: &[Csr<T>]) -> (usize, usize) {
    let (nrows, ncols) = (parts[0].nrows, parts[0].ncols);
    for (idx, p) in parts.iter().enumerate().skip(1) {
        assert!(
            p.nrows == nrows && p.ncols == ncols,
            "Csr::sum_all shape mismatch at part {idx}: ({}, {}) != ({nrows}, {ncols})",
            p.nrows,
            p.ncols
        );
    }
    (nrows, ncols)
}

pub(super) fn sum_block_rows<T: Scalar>(parts: &[Csr<T>], r0: usize, r1: usize) -> SumBlock<T> {
    let mut counts = Vec::with_capacity(r1 - r0);
    let (mut idx, mut val): (Vec<u32>, Vec<T>) = (Vec::new(), Vec::new());
    let mut row: Vec<(u32, T)> = Vec::new();
    for r in r0..r1 {
        row.clear();
        for p in parts {
            let (cs, vs) = p.row(r);
            row.extend(cs.iter().copied().zip(vs.iter().copied()));
        }
        row.sort_unstable_by_key(|e| e.0);
        let start = idx.len();
        for &(cc, v) in &row {
            if idx.len() > start && *idx.last().unwrap() == cc {
                *val.last_mut().unwrap() += v;
            } else {
                idx.push(cc);
                val.push(v);
            }
        }
        counts.push(idx.len() - start);
    }
    (counts, idx, val)
}

pub(super) fn assemble_sum_blocks<T: Scalar>(
    nrows: usize,
    ncols: usize,
    blocks: Vec<SumBlock<T>>,
) -> Csr<T> {
    let mut indptr = Vec::with_capacity(nrows + 1);
    indptr.push(0usize);
    let total: usize = blocks.iter().map(|b| b.1.len()).sum();
    let (mut indices, mut values) = (Vec::with_capacity(total), Vec::with_capacity(total));
    for (counts, idx, val) in blocks {
        for cnt in counts {
            indptr.push(indptr.last().unwrap() + cnt);
        }
        indices.extend(idx);
        values.extend(val);
    }
    Csr::try_new(nrows, ncols, indptr, indices, values).expect("valid summed CSR")
}

impl<T: Scalar> Csr<T> {
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        let len = nrows.checked_add(1).expect("Csr::zeros nrows overflow");
        Self::try_new(nrows, ncols, vec![0; len], vec![], vec![]).expect("valid zero CSR")
    }
    pub fn from_diag(d: &[T]) -> Self {
        let n = d.len();
        assert!(n <= u32::MAX as usize, "Csr::from_diag dimension {n} overflows u32");
        let indices: Vec<u32> = (0..n).map(|i| u32::try_from(i).expect("u32 index")).collect();
        Self::try_new(n, n, (0..=n).collect(), indices, d.to_vec()).expect("valid diagonal CSR")
    }
    /// From `usize` column indices (rows sorted by the caller); panics on an inconsistent matrix, which is a bug of the caller.
    pub fn from_raw(nrows: usize, ncols: usize, indptr: Vec<usize>, indices: Vec<usize>, values: Vec<T>) -> Self {
        let u32_indices: Vec<u32> = indices
            .into_iter()
            .map(|c| u32::try_from(c).expect("Csr::from_raw column index overflows u32"))
            .collect();
        Self::try_new(nrows, ncols, indptr, u32_indices, values).expect("a consistent CSR matrix")
    }
    #[inline]
    pub fn row_range(&self, r: usize) -> std::ops::Range<usize> {
        self.indptr[r]..self.indptr[r + 1]
    }
    /// Entry `(r, c)` (zero where not stored); rows must be sorted.
    pub fn get(&self, r: usize, c: usize) -> T {
        let (cols, vals) = self.row(r);
        if c >= self.ncols {
            return T::ZERO;
        }
        let Ok(c_u32) = u32::try_from(c) else {
            return T::ZERO;
        };
        match cols.binary_search(&c_u32) {
            Ok(p) => vals[p],
            Err(_) => T::ZERO,
        }
    }
    /// `A x` as a new vector.
    pub fn matvec(&self, x: &[T]) -> Vec<T> {
        let mut y = vec![T::ZERO; self.nrows];
        self.spmv(x, &mut y);
        y
    }
    pub fn matvec_into(&self, x: &[T], y: &mut [T]) {
        self.spmv(x, y)
    }
    /// `self + alpha * other` (same shape).
    pub fn add_scaled(&self, other: &Csr<T>, alpha: T) -> Csr<T> {
        Csr::axpby(T::ONE, self, alpha, other)
    }
    pub fn scaled(&self, a: T) -> Csr<T> {
        let mut m = self.clone();
        m.values.iter_mut().for_each(|v| *v *= a);
        m
    }
    /// The stored values, mutable (same pattern): the row and column scaling of a generalised problem.
    pub fn values_mut(&mut self) -> &mut [T] {
        &mut self.values
    }
    /// Rows `rows` and columns `cols` (index lists in the given order); entries whose column is not in `cols` are dropped.
    pub fn submatrix(&self, rows: &[usize], cols: &[usize]) -> Csr<T> {
        assert!(cols.len() <= u32::MAX as usize, "Csr::submatrix column count overflows u32");
        let mut map = vec![u32::MAX; self.ncols];
        for (j, &c) in cols.iter().enumerate() {
            map[c] = u32::try_from(j).expect("u32 index");
        }
        let mut indptr = vec![0usize];
        let mut indices: Vec<u32> = Vec::new();
        let mut values: Vec<T> = Vec::new();
        let mut row: Vec<(u32, T)> = Vec::new();
        for &r in rows {
            row.clear();
            let (cs, vs) = self.row(r);
            for (&cc, &v) in cs.iter().zip(vs) {
                let j = map[cc as usize];
                if j != u32::MAX {
                    row.push((j, v));
                }
            }
            row.sort_unstable_by_key(|e| e.0);
            for &(j, v) in &row {
                indices.push(j);
                values.push(v);
            }
            indptr.push(indices.len());
        }
        Self::try_new(rows.len(), cols.len(), indptr, indices, values).expect("valid submatrix CSR")
    }
    pub fn block_diag(blocks: &[&Csr<T>]) -> Csr<T> {
        let nrows = blocks
            .iter()
            .try_fold(0usize, |acc, b| acc.checked_add(b.nrows))
            .expect("Csr::block_diag nrows overflow");
        let ncols = blocks
            .iter()
            .try_fold(0usize, |acc, b| acc.checked_add(b.ncols))
            .expect("Csr::block_diag ncols overflow");
        assert!(ncols <= u32::MAX as usize, "Csr::block_diag ncols {ncols} overflows u32");
        let mut indptr = vec![0usize];
        let (mut indices, mut values) = (Vec::new(), Vec::new());
        let mut co = 0usize;
        for b in blocks {
            let co_u32 = u32::try_from(co).expect("u32 column offset");
            for r in 0..b.nrows {
                let (cs, vs) = b.row(r);
                for (&cc, &v) in cs.iter().zip(vs) {
                    indices.push(cc.checked_add(co_u32).expect("Csr::block_diag index overflow"));
                    values.push(v);
                }
                indptr.push(indices.len());
            }
            co += b.ncols;
        }
        Self::try_new(nrows, ncols, indptr, indices, values).expect("valid block_diag CSR")
    }
    /// Sum of matrices of one shape. Blocks of rows are merged on all cores; a row of the sum is a small merge of at most `parts.len()`
    /// sorted lists.
    pub fn sum_all(mut parts: Vec<Csr<T>>) -> Csr<T> {
        use rayon::prelude::*;
        assert!(!parts.is_empty(), "Csr::sum_all requires at least one part");
        if parts.len() == 1 {
            return parts.pop().unwrap();
        }
        let (nrows, ncols) = validate_sum_parts(&parts);
        let chunk = nrows.div_ceil(crate::par::threads()).max(1);
        let parts_ref = &parts;
        let blocks: Vec<SumBlock<T>> = (0..nrows.div_ceil(chunk))
            .into_par_iter()
            .map(|b| {
                let (r0, r1) = (b * chunk, ((b + 1) * chunk).min(nrows));
                sum_block_rows(parts_ref, r0, r1)
            })
            .collect();
        // the parts are as large as the sum: free them before the blocks are copied into the final arrays
        drop(parts);
        assemble_sum_blocks(nrows, ncols, blocks)
    }
    /// Symmetrised copy `(A + A^T)/2`.
    pub fn symmetrized(&self) -> Csr<T> {
        self.add_scaled(&self.transpose(), T::ONE).scaled(T::from_f64(0.5))
    }
    /// Entrywise conjugate (identity for real matrices).
    pub fn conj(&self) -> Csr<T> {
        self.map(|v| v.conj())
    }
}

impl Csr<C64> {
    /// The real part of every entry (same pattern).
    pub fn real_part(&self) -> Csr<f64> {
        self.map(|z| z.re)
    }
    /// All imaginary parts are exactly zero.
    pub fn is_real(&self) -> bool {
        self.values.iter().all(|z| z.im == 0.0)
    }
}
