//! Sparse matrices: coordinate (triplet) assembly `Coo`, compressed rows `Csr` and the single-precision `Csr32`.
//!
//! `Coo` collects (row, column, value) triplets and turns them into a `Csr` (columns sorted, duplicates summed, which is
//! how element matrices are added up into the global matrix). `Csr` stores 8-byte column indices and `f64` values;
//! `Csr32` stores 4-byte indices and `f32` values for the multigrid cycles, whose products are limited by memory
//! bandwidth and only need single precision (the sums are still accumulated in double precision). Matrix-vector
//! products and sums of matrices run on the shared thread pool (`par`) above a size at which that pays off.

use rayon::prelude::*;
use std::ops::Range;

/// Matrix-vector products with fewer non-zeros than this run on the calling thread (the pool needs ~12 microseconds
/// to start a product, which pays off from about 30,000 non-zeros: 105,000 non-zeros take 68 us serial, 29 us on 20 threads).
const PAR_MIN_NNZ32: usize = 40_000;
const PAR_MIN_NNZ: usize = 80_000;

/// Compressed rows with 32-bit column indices and single-precision values, for matrix-vector
/// products that only need to be accurate to about 1e-7 (multigrid cycles used as a
/// preconditioner): 8 bytes per entry instead of 16, and these products are limited by memory
/// bandwidth. Products accumulate in double precision.
#[derive(Clone, Debug, Default)]
pub struct Csr32 {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<u32>,
    pub data: Vec<f32>,
}

impl Csr32 {
    pub fn from_csr(a: &Csr) -> Csr32 {
        Csr32 { nrows: a.nrows, ncols: a.ncols, indptr: a.indptr.clone(), indices: a.indices.iter().map(|&c| c as u32).collect(), data: a.data.iter().map(|&v| v as f32).collect() }
    }
    pub fn nnz(&self) -> usize {
        self.indices.len()
    }
    #[inline]
    fn row_dot(&self, r: usize, x: &[f64]) -> f64 {
        let mut s = 0.0;
        for k in self.indptr[r]..self.indptr[r + 1] {
            s += self.data[k] as f64 * x[self.indices[k] as usize];
        }
        s
    }
    pub fn matvec_into(&self, x: &[f64], y: &mut [f64]) {
        if self.nnz() < PAR_MIN_NNZ32 {
            for r in 0..self.nrows {
                y[r] = self.row_dot(r, x);
            }
            return;
        }
        let chunk = self.nrows.div_ceil(crate::par::threads()).max(1);
        y.par_chunks_mut(chunk).enumerate().for_each(|(ci, ys)| {
            for (i, v) in ys.iter_mut().enumerate() {
                *v = self.row_dot(ci * chunk + i, x);
            }
        });
    }
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.nrows];
        self.matvec_into(x, &mut y);
        y
    }
}

/// Coordinate-format builder; duplicate entries are summed when converted.
#[derive(Clone, Debug, Default)]
pub struct Coo {
    pub nrows: usize,
    pub ncols: usize,
    pub rows: Vec<u32>,
    pub cols: Vec<u32>,
    pub vals: Vec<f64>,
}

impl Coo {
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Coo { nrows, ncols, ..Default::default() }
    }
    pub fn with_capacity(nrows: usize, ncols: usize, cap: usize) -> Self {
        Coo {
            nrows,
            ncols,
            rows: Vec::with_capacity(cap),
            cols: Vec::with_capacity(cap),
            vals: Vec::with_capacity(cap),
        }
    }
    #[inline]
    pub fn push(&mut self, r: usize, c: usize, v: f64) {
        debug_assert!(r < self.nrows && c < self.ncols);
        self.rows.push(r as u32);
        self.cols.push(c as u32);
        self.vals.push(v);
    }
    /// Adds the dense block `local[a*cols.len() + b]` at (rows[a], cols[b]).
    pub fn push_block(&mut self, rows: &[usize], cols: &[usize], local: &[f64]) {
        debug_assert_eq!(local.len(), rows.len() * cols.len());
        for (a, &r) in rows.iter().enumerate() {
            for (b, &c) in cols.iter().enumerate() {
                let v = local[a * cols.len() + b];
                if v != 0.0 {
                    self.push(r, c, v);
                }
            }
        }
    }
    pub fn to_csr(&self) -> Csr {
        let n = self.nrows;
        // counting sort by row, then sort each row by column and merge duplicates
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
        let mut indices = Vec::with_capacity(nnz);
        let mut data = Vec::with_capacity(nnz);
        indptr.push(0);
        let mut row: Vec<(u32, f64)> = Vec::new();
        for r in 0..n {
            row.clear();
            for &k in &order[count[r]..count[r + 1]] {
                row.push((self.cols[k], self.vals[k]));
            }
            row.sort_unstable_by_key(|e| e.0);
            let mut i = 0;
            while i < row.len() {
                let c = row[i].0;
                let mut s = 0.0;
                while i < row.len() && row[i].0 == c {
                    s += row[i].1;
                    i += 1;
                }
                indices.push(c as usize);
                data.push(s);
            }
            indptr.push(indices.len());
        }
        Csr { nrows: n, ncols: self.ncols, indptr, indices, data }
    }
}

/// Compressed sparse row matrix with sorted column indices and no duplicates.
#[derive(Clone, Debug)]
pub struct Csr {
    pub nrows: usize,
    pub ncols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<f64>,
}

impl Csr {
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Csr { nrows, ncols, indptr: vec![0; nrows + 1], indices: vec![], data: vec![] }
    }
    pub fn identity(n: usize) -> Self {
        Csr { nrows: n, ncols: n, indptr: (0..=n).collect(), indices: (0..n).collect(), data: vec![1.0; n] }
    }
    pub fn from_diag(d: &[f64]) -> Self {
        let n = d.len();
        Csr { nrows: n, ncols: n, indptr: (0..=n).collect(), indices: (0..n).collect(), data: d.to_vec() }
    }
    pub fn nnz(&self) -> usize {
        self.data.len()
    }
    #[inline]
    pub fn row_range(&self, r: usize) -> Range<usize> {
        self.indptr[r]..self.indptr[r + 1]
    }
    pub fn get(&self, r: usize, c: usize) -> f64 {
        let rg = self.row_range(r);
        match self.indices[rg.clone()].binary_search(&c) {
            Ok(p) => self.data[rg.start + p],
            Err(_) => 0.0,
        }
    }
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.nrows];
        self.matvec_into(x, &mut y);
        y
    }
    pub fn matvec_into(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.ncols);
        let row = |r: usize| -> f64 {
            let mut s = 0.0;
            for k in self.indptr[r]..self.indptr[r + 1] {
                s += self.data[k] * x[self.indices[k]];
            }
            s
        };
        // large products run on all cores
        if self.nnz() < PAR_MIN_NNZ {
            for r in 0..self.nrows {
                y[r] = row(r);
            }
            return;
        }
        let chunk = self.nrows.div_ceil(crate::par::threads_up_to(16)).max(1);
        y.par_chunks_mut(chunk).enumerate().for_each(|(ci, ys)| {
            for (i, v) in ys.iter_mut().enumerate() {
                *v = row(ci * chunk + i);
            }
        });
    }
    /// Sum of matrices of one shape. Blocks of rows are merged on all cores; the rows of the parts
    /// are sorted, so a row of the sum is a small merge of at most `parts.len()` sorted lists.
    pub fn sum_all(mut parts: Vec<Csr>) -> Csr {
        if parts.len() == 1 {
            return parts.pop().unwrap();
        }
        let (nrows, ncols) = (parts[0].nrows, parts[0].ncols);
        let chunk = nrows.div_ceil(crate::par::threads()).max(1);
        let parts_ref = &parts;
        let blocks: Vec<(Vec<usize>, Vec<usize>, Vec<f64>)> = (0..nrows.div_ceil(chunk))
            .into_par_iter()
            .map(|b| {
                let (r0, r1) = (b * chunk, ((b + 1) * chunk).min(nrows));
                let mut counts = Vec::with_capacity(r1 - r0);
                let (mut idx, mut val) = (Vec::new(), Vec::new());
                let mut row: Vec<(usize, f64)> = Vec::new();
                for r in r0..r1 {
                    row.clear();
                    for p in parts_ref {
                        for k in p.indptr[r]..p.indptr[r + 1] {
                            row.push((p.indices[k], p.data[k]));
                        }
                    }
                    row.sort_unstable_by_key(|e| e.0);
                    let start = idx.len();
                    for &(c, v) in &row {
                        if idx.len() > start && *idx.last().unwrap() == c {
                            *val.last_mut().unwrap() += v;
                        } else {
                            idx.push(c);
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
        let (mut indices, mut data) = (Vec::with_capacity(total), Vec::with_capacity(total));
        for (counts, idx, val) in blocks {
            for c in counts {
                indptr.push(indptr.last().unwrap() + c);
            }
            indices.extend(idx);
            data.extend(val);
        }
        Csr { nrows, ncols, indptr, indices, data }
    }
    pub fn transpose(&self) -> Csr {
        let mut count = vec![0usize; self.ncols + 1];
        for &c in &self.indices {
            count[c + 1] += 1;
        }
        for i in 0..self.ncols {
            count[i + 1] += count[i];
        }
        let mut next = count.clone();
        let mut indices = vec![0usize; self.nnz()];
        let mut data = vec![0.0; self.nnz()];
        for r in 0..self.nrows {
            for k in self.row_range(r) {
                let c = self.indices[k];
                let p = next[c];
                indices[p] = r;
                data[p] = self.data[k];
                next[c] += 1;
            }
        }
        Csr { nrows: self.ncols, ncols: self.nrows, indptr: count, indices, data }
    }
    /// `self + alpha * other` (same shape).
    pub fn add_scaled(&self, other: &Csr, alpha: f64) -> Csr {
        assert_eq!((self.nrows, self.ncols), (other.nrows, other.ncols));
        let mut indptr = vec![0usize];
        let mut indices = Vec::with_capacity(self.nnz() + other.nnz());
        let mut data = Vec::with_capacity(self.nnz() + other.nnz());
        for r in 0..self.nrows {
            let (mut i, ie) = (self.indptr[r], self.indptr[r + 1]);
            let (mut j, je) = (other.indptr[r], other.indptr[r + 1]);
            while i < ie || j < je {
                let ci = if i < ie { self.indices[i] } else { usize::MAX };
                let cj = if j < je { other.indices[j] } else { usize::MAX };
                if ci == cj {
                    indices.push(ci);
                    data.push(self.data[i] + alpha * other.data[j]);
                    i += 1;
                    j += 1;
                } else if ci < cj {
                    indices.push(ci);
                    data.push(self.data[i]);
                    i += 1;
                } else {
                    indices.push(cj);
                    data.push(alpha * other.data[j]);
                    j += 1;
                }
            }
            indptr.push(indices.len());
        }
        Csr { nrows: self.nrows, ncols: self.ncols, indptr, indices, data }
    }
    pub fn scaled(&self, a: f64) -> Csr {
        let mut m = self.clone();
        m.data.iter_mut().for_each(|v| *v *= a);
        m
    }
    /// Rows `rows` and columns `cols` (both index lists, in the given order); entries whose
    /// column is not in `cols` are dropped.
    pub fn submatrix(&self, rows: &[usize], cols: &[usize]) -> Csr {
        let mut map = vec![usize::MAX; self.ncols];
        for (j, &c) in cols.iter().enumerate() {
            map[c] = j;
        }
        let mut indptr = vec![0usize];
        let mut indices = Vec::new();
        let mut data = Vec::new();
        let mut row: Vec<(usize, f64)> = Vec::new();
        for &r in rows {
            row.clear();
            for k in self.row_range(r) {
                let j = map[self.indices[k]];
                if j != usize::MAX {
                    row.push((j, self.data[k]));
                }
            }
            row.sort_unstable_by_key(|e| e.0);
            for &(j, v) in &row {
                indices.push(j);
                data.push(v);
            }
            indptr.push(indices.len());
        }
        Csr { nrows: rows.len(), ncols: cols.len(), indptr, indices, data }
    }
    pub fn block_diag(blocks: &[&Csr]) -> Csr {
        let nrows = blocks.iter().map(|b| b.nrows).sum();
        let ncols = blocks.iter().map(|b| b.ncols).sum();
        let mut indptr = vec![0usize];
        let mut indices = Vec::new();
        let mut data = Vec::new();
        let mut co = 0;
        for b in blocks {
            for r in 0..b.nrows {
                for k in b.row_range(r) {
                    indices.push(b.indices[k] + co);
                    data.push(b.data[k]);
                }
                indptr.push(indices.len());
            }
            co += b.ncols;
        }
        Csr { nrows, ncols, indptr, indices, data }
    }
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.nrows.min(self.ncols)).map(|i| self.get(i, i)).collect()
    }
    pub fn to_dense(&self) -> Vec<f64> {
        let mut d = vec![0.0; self.nrows * self.ncols];
        for r in 0..self.nrows {
            for k in self.row_range(r) {
                d[r * self.ncols + self.indices[k]] = self.data[k];
            }
        }
        d
    }
    /// Symmetrised copy `(A + A^T)/2`.
    pub fn symmetrized(&self) -> Csr {
        self.add_scaled(&self.transpose(), 1.0).scaled(0.5)
    }
}

/// Dot product.
#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assemble_and_ops() {
        let mut c = Coo::new(3, 3);
        c.push(0, 0, 1.0);
        c.push(0, 0, 2.0);
        c.push(1, 2, 4.0);
        c.push(2, 1, 5.0);
        c.push(1, 1, -1.0);
        let a = c.to_csr();
        assert_eq!(a.get(0, 0), 3.0);
        assert_eq!(a.nnz(), 4);
        assert_eq!(a.matvec(&[1.0, 2.0, 3.0]), vec![3.0, 10.0, 10.0]);
        let t = a.transpose();
        assert_eq!(t.get(2, 1), 4.0);
        let s = a.add_scaled(&t, 2.0);
        assert_eq!(s.get(1, 2), 4.0 + 2.0 * 5.0);
        let sub = a.submatrix(&[2, 1], &[1, 2]);
        assert_eq!(sub.to_dense(), vec![5.0, 0.0, -1.0, 4.0]);
    }
}
