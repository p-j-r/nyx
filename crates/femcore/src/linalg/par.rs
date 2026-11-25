//! Deterministic chunk parallelism.
//!
//! All parallel kernels split the rows into fixed-size chunks (independent of the number of
//! threads), compute per-chunk partial results, and fold the partials **in chunk order** on the
//! calling thread. The result is therefore bit-identical for every thread count, and identical
//! to the serial run (`NYX_THREADS=1`). That is the determinism rule nyx already uses
//! ("fixed pieces"), kept as the contract for every kernel in `solvers`.
//!
//! The work runs on nyx's one shared rayon pool (`crate::par`, sized by `NYX_THREADS`); `ensure_pool` makes sure
//! that pool exists before the first rayon call, so the solvers never fall back to rayon's own default
//! pool of a different size. Waking a sleeping worker cost 50-150 µs per region on the Linux VM the port was
//! written on; the thresholds below keep small problems serial.

/// Minimum number of scalar elements touched by one kernel call before it is run in parallel.
pub const PAR_MIN_ELEMS: usize = 1 << 17;

/// Creates the shared nyx pool (honouring `NYX_THREADS`) if no kernel has done so yet.
#[inline]
fn ensure_pool() {
    let _ = crate::par::threads();
}

/// `f(0), f(1), …, f(nchunks-1)` collected in order, computed in parallel when `parallel`.
pub fn chunk_map<R, F>(nchunks: usize, parallel: bool, f: F) -> Vec<R>
where
    R: Send,
    F: Fn(usize) -> R + Sync + Send,
{
    if parallel && nchunks > 1 {
        use rayon::prelude::*;
        ensure_pool();
        return (0..nchunks).into_par_iter().map(f).collect();
    }
    (0..nchunks).map(f).collect()
}

/// Applies `f(chunk_index, chunk)` to consecutive `chunk`-sized pieces of `data`.
pub fn for_chunks_mut<T, F>(data: &mut [T], chunk: usize, parallel: bool, f: F)
where
    T: Send,
    F: Fn(usize, &mut [T]) + Sync + Send,
{
    if parallel && data.len() > chunk {
        use rayon::prelude::*;
        ensure_pool();
        data.par_chunks_mut(chunk).enumerate().for_each(|(i, c)| f(i, c));
        return;
    }
    for (i, c) in data.chunks_mut(chunk).enumerate() {
        f(i, c);
    }
}

/// Applies `f` to every item, in parallel when `parallel`. The items are independent (typically disjoint
/// `&mut` slices of one output), so the result does not depend on the thread count.
pub fn for_each_mut<X, F>(items: &mut [X], parallel: bool, f: F)
where
    X: Send,
    F: Fn(&mut X) + Sync + Send,
{
    if parallel && items.len() > 1 {
        use rayon::prelude::*;
        ensure_pool();
        items.par_iter_mut().for_each(f);
        return;
    }
    items.iter_mut().for_each(f);
}

/// Whether a kernel touching `elems` scalars should run in parallel (never with `NYX_THREADS=1`).
#[inline]
pub fn worth_parallel(elems: usize) -> bool {
    elems >= PAR_MIN_ELEMS && crate::par::threads() > 1
}

/// Shared mutable access to disjoint row ranges of a column-major block (`ncols` columns of
/// length `n`, stored contiguously).
///
/// This is the one `unsafe` building block of the parallel basis kernels. Chunks of rows are
/// processed concurrently; each chunk only touches its own rows `s..e` of every column, so the
/// ranges handed out never overlap.
/// Pointer to `ncols` contiguous columns of length `n` inside `&'a mut [T]`, shared across
/// threads when each worker claims disjoint tiles.
pub struct RawCols<'a, T> {
    ptr: *mut T,
    n: usize,
    ncols: usize,
    _marker: std::marker::PhantomData<&'a mut [T]>,
}

// SAFETY: `RawCols` carries `'a` provenance from `&'a mut [T]`; sharing between threads is sound
// when concurrent callers claim disjoint `(c, start..end)` tiles.
unsafe impl<T: Send> Send for RawCols<'_, T> {}
unsafe impl<T: Send> Sync for RawCols<'_, T> {}

impl<'a, T> RawCols<'a, T> {
    pub fn new(data: &'a mut [T], n: usize, ncols: usize) -> Self {
        let total = n.checked_mul(ncols).expect("RawCols size overflow");
        assert!(data.len() >= total, "RawCols buffer too short: {} < {total}", data.len());
        Self { ptr: data.as_mut_ptr(), n, ncols, _marker: std::marker::PhantomData }
    }

    /// Mutable tile `start..end` of column `c` with release bounds checks.
    #[allow(clippy::mut_from_ref)]
    #[inline]
    pub fn chunk_mut(&self, c: usize, start: usize, end: usize) -> &'a mut [T] {
        assert!(
            c < self.ncols && start <= end && end <= self.n,
            "RawCols::chunk_mut out of bounds: c={c}/{} range={start}..{end}/{}",
            self.ncols,
            self.n
        );
        // SAFETY: `c < self.ncols`, `start <= end <= self.n`, and `self.n * self.ncols <= data.len()`
        // (verified with checked multiplication in `new`), so the slice lies inside the allocation.
        unsafe { std::slice::from_raw_parts_mut(self.ptr.add(c * self.n + start), end - start) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_cols_disjoint_tiles_work() {
        let mut buf = vec![0usize; 12];
        let cols = RawCols::new(&mut buf, 4, 3);
        cols.chunk_mut(0, 0, 2).copy_from_slice(&[1, 2]);
        cols.chunk_mut(2, 1, 4).copy_from_slice(&[7, 8, 9]);
        assert_eq!(&buf[0..4], &[1, 2, 0, 0]);
        assert_eq!(&buf[8..12], &[0, 7, 8, 9]);
    }

    #[test]
    #[should_panic(expected = "RawCols size overflow")]
    fn raw_cols_rejects_size_overflow() {
        let mut buf = [0u8; 4];
        let _ = RawCols::new(&mut buf, usize::MAX, 2);
    }

    #[test]
    #[should_panic]
    fn raw_cols_rejects_short_buffer() {
        let mut buf = [0u8; 7];
        let _ = RawCols::new(&mut buf, 4, 2);
    }

    #[test]
    #[should_panic]
    fn raw_cols_rejects_oob_col() {
        let mut buf = [0u8; 8];
        let cols = RawCols::new(&mut buf, 4, 2);
        let _ = cols.chunk_mut(2, 0, 1);
    }

    #[test]
    #[should_panic]
    fn raw_cols_rejects_oob_rows() {
        let mut buf = [0u8; 8];
        let cols = RawCols::new(&mut buf, 4, 2);
        let _ = cols.chunk_mut(0, 2, 5);
    }
}
