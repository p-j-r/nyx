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
pub struct RawCols<T> {
    ptr: *mut T,
    n: usize,
    ncols: usize,
}

// SAFETY: `RawCols` is only a pointer + shape; the safety contract of `rows` (disjoint ranges)
// is what makes sharing it between threads sound.
unsafe impl<T: Send> Send for RawCols<T> {}
unsafe impl<T: Send> Sync for RawCols<T> {}

impl<T> RawCols<T> {
    pub fn new(data: &mut [T], n: usize, ncols: usize) -> Self {
        assert!(data.len() >= n * ncols);
        Self { ptr: data.as_mut_ptr(), n, ncols }
    }

    /// Mutable rows `s..e` of column `c`.
    ///
    /// # Safety
    /// No other live reference may overlap rows `s..e` of column `c`, and the originating slice
    /// must outlive the returned reference.
    #[allow(clippy::mut_from_ref)]
    #[inline]
    pub unsafe fn rows(&self, c: usize, s: usize, e: usize) -> &mut [T] {
        debug_assert!(c < self.ncols && s <= e && e <= self.n);
        // SAFETY: in bounds by the assertions above; exclusivity is the caller's obligation.
        unsafe { std::slice::from_raw_parts_mut(self.ptr.add(c * self.n + s), e - s) }
    }

    /// Raw pointer to element `(row, col)`, for gather/scatter copies.
    #[inline]
    pub fn at(&self, row: usize, c: usize) -> *mut T {
        debug_assert!(c < self.ncols && row <= self.n);
        // SAFETY: only pointer arithmetic inside the allocation (one-past-the-end allowed).
        unsafe { self.ptr.add(c * self.n + row) }
    }
}
