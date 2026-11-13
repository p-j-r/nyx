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
    pub nrows: usize,
    pub ncols: usize,
    pub rows: Vec<u32>,
    pub cols: Vec<u32>,
    pub vals: Vec<T>,
}

pub type Coo = Triplets<f64>;
pub type CCoo = Triplets<C64>;

impl<T: Scalar> Triplets<T> {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Triplets { nrows, ncols, rows: vec![], cols: vec![], vals: vec![] }
    }
    pub fn with_capacity(nrows: usize, ncols: usize, cap: usize) -> Self {
        Triplets { nrows, ncols, rows: Vec::with_capacity(cap), cols: Vec::with_capacity(cap), vals: Vec::with_capacity(cap) }
    }
    #[inline]
    pub fn push(&mut self, r: usize, c: usize, v: T) {
        debug_assert!(r < self.nrows && c < self.ncols);
        self.rows.push(r as u32);
        self.cols.push(c as u32);
        self.vals.push(v);
    }
    /// Adds the dense block `local[a*cols.len() + b]` at `(rows[a], cols[b])`; exact zeros are skipped.
    pub fn push_block(&mut self, rows: &[usize], cols: &[usize], local: &[T]) {
        debug_assert_eq!(local.len(), rows.len() * cols.len());
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
        let mut count = vec![0usize; n + 1];
        for &r in &self.rows {
            count[r as usize + 1] += 1;
        }
        for i in 0..n {
            count[i + 1] += count[i];
        }
        let nnz = self.rows.len();
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
                row.push((self.cols[k], self.vals[k]));
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
        Csr { nrows: n, ncols: self.ncols, indptr, indices, values }
    }
}

impl<T: Scalar> Csr<T> {
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Csr { nrows, ncols, indptr: vec![0; nrows + 1], indices: vec![], values: vec![] }
    }
    pub fn from_diag(d: &[T]) -> Self {
        let n = d.len();
        Csr { nrows: n, ncols: n, indptr: (0..=n).collect(), indices: (0..n as u32).collect(), values: d.to_vec() }
    }
    /// From `usize` column indices (rows sorted by the caller); panics on an inconsistent matrix, which is a bug of the caller.
    pub fn from_raw(nrows: usize, ncols: usize, indptr: Vec<usize>, indices: Vec<usize>, values: Vec<T>) -> Self {
        Self::try_new(nrows, ncols, indptr, indices.into_iter().map(|c| c as u32).collect(), values).expect("a consistent CSR matrix")
    }
    #[inline]
    pub fn row_range(&self, r: usize) -> std::ops::Range<usize> {
        self.indptr[r]..self.indptr[r + 1]
    }
    /// Entry `(r, c)` (zero where not stored); rows must be sorted.
    pub fn get(&self, r: usize, c: usize) -> T {
        let (cols, vals) = self.row(r);
        match cols.binary_search(&(c as u32)) {
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
        let mut map = vec![u32::MAX; self.ncols];
        for (j, &c) in cols.iter().enumerate() {
            map[c] = j as u32;
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
        Csr { nrows: rows.len(), ncols: cols.len(), indptr, indices, values }
    }
    pub fn block_diag(blocks: &[&Csr<T>]) -> Csr<T> {
        let nrows = blocks.iter().map(|b| b.nrows).sum();
        let ncols = blocks.iter().map(|b| b.ncols).sum();
        let mut indptr = vec![0usize];
        let (mut indices, mut values) = (Vec::new(), Vec::new());
        let mut co = 0u32;
        for b in blocks {
            for r in 0..b.nrows {
                let (cs, vs) = b.row(r);
                for (&cc, &v) in cs.iter().zip(vs) {
                    indices.push(cc + co);
                    values.push(v);
                }
                indptr.push(indices.len());
            }
            co += b.ncols as u32;
        }
        Csr { nrows, ncols, indptr, indices, values }
    }
    /// Sum of matrices of one shape. Blocks of rows are merged on all cores; a row of the sum is a small merge of at most `parts.len()`
    /// sorted lists.
    pub fn sum_all(mut parts: Vec<Csr<T>>) -> Csr<T> {
        use rayon::prelude::*;
        if parts.len() == 1 {
            return parts.pop().unwrap();
        }
        let (nrows, ncols) = (parts[0].nrows, parts[0].ncols);
        let chunk = nrows.div_ceil(crate::par::threads()).max(1);
        let parts_ref = &parts;
        type Block<T> = (Vec<usize>, Vec<u32>, Vec<T>);
        let blocks: Vec<Block<T>> = (0..nrows.div_ceil(chunk))
            .into_par_iter()
            .map(|b| {
                let (r0, r1) = (b * chunk, ((b + 1) * chunk).min(nrows));
                let mut counts = Vec::with_capacity(r1 - r0);
                let (mut idx, mut val): (Vec<u32>, Vec<T>) = (Vec::new(), Vec::new());
                let mut row: Vec<(u32, T)> = Vec::new();
                for r in r0..r1 {
                    row.clear();
                    for p in parts_ref {
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
            })
            .collect();
        // the parts are as large as the sum: free them before the blocks are copied into the final arrays
        drop(parts);
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
        Csr { nrows, ncols, indptr, indices, values }
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
