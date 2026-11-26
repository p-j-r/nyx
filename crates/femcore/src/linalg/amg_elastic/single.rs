//! The single-precision matrix copy used by the multigrid cycles, and a dot product.
//!
//! `Csr32` stores 4-byte indices and `f32` values: the cycle is only a preconditioner and its products are limited by memory bandwidth, so
//! single precision halves the traffic (the sums are still accumulated in double precision). The hierarchy itself is built on the library
//! matrix type `Csr<f64>`.

use crate::linalg::csr::RCsr;

const PAR_MIN_NNZ32: usize = 40_000;
const ROW_CHUNK32: usize = 512;


/// Compressed rows with 32-bit column indices and single-precision values, for matrix-vector
/// products that only need to be accurate to about 1e-7 (multigrid cycles used as a
/// preconditioner): 8 bytes per entry instead of 16, and these products are limited by memory
/// bandwidth. Products accumulate in double precision.
#[derive(Clone, Debug, Default)]
pub struct Csr32 {
    pub nrows: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<u32>,
    pub data: Vec<f32>,
}


impl Csr32 {
    pub fn from_csr(a: &RCsr) -> Csr32 {
        Csr32 { nrows: a.nrows(), indptr: a.indptr().to_vec(), indices: a.indices().to_vec(), data: a.values().iter().map(|&v| v as f32).collect() }
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
        let par = self.nnz() >= PAR_MIN_NNZ32 && crate::par::threads() > 1;
        crate::par::for_chunks_mut(y, ROW_CHUNK32, par, |ci, ys| {
            for (i, v) in ys.iter_mut().enumerate() {
                *v = self.row_dot(ci * ROW_CHUNK32 + i, x);
            }
        });
    }
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.nrows];
        self.matvec_into(x, &mut y);
        y
    }
}


/// Dot product.
#[inline]
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
