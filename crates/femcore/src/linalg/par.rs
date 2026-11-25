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
#[cfg(not(miri))]
pub const PAR_MIN_ELEMS: usize = 1 << 17;
/// Lowered parallel threshold under Miri so small `n <= 32` tests exercise parallel `RawCols` tile paths.
#[cfg(miri)]
pub const PAR_MIN_ELEMS: usize = 8;

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
    fn raw_cols_parallel_tiles_small_n_disjoint() {
        let n = 16;
        let ncols = 3;
        let chunk = 4;
        let mut buf = vec![0usize; n * ncols];
        let cols = RawCols::new(&mut buf, n, ncols);
        let nchunks = n.div_ceil(chunk);
        // Force parallel path when `crate::par::threads() > 1` (and under Miri with `PAR_MIN_ELEMS = 8`).
        let run_par = worth_parallel(n * ncols) || nchunks > 1;
        let partials = chunk_map(nchunks, run_par, |ci| {
            let s = ci * chunk;
            let e = (s + chunk).min(n);
            let mut sum = 0usize;
            for c in 0..ncols {
                let tile = cols.chunk_mut(c, s, e);
                for (offset, slot) in tile.iter_mut().enumerate() {
                    let v = (c + 1) * 100 + s + offset;
                    *slot = v;
                    sum += v;
                }
            }
            sum
        });
        let total: usize = partials.into_iter().sum();
        assert_eq!(total, buf.iter().sum::<usize>());
        for c in 0..ncols {
            for r in 0..n {
                assert_eq!(buf[c * n + r], (c + 1) * 100 + r);
            }
        }
    }

    #[test]
    #[should_panic]
    fn raw_cols_rejects_oob_rows() {
        let mut buf = [0u8; 8];
        let cols = RawCols::new(&mut buf, 4, 2);
        let _ = cols.chunk_mut(0, 2, 5);
    }
}

#[cfg(kani)]
mod kani_proofs {
    use super::*;

    #[kani::proof]
    #[kani::unwind(16)]
    fn verify_raw_cols_disjoint_chunks() {
        let n: usize = kani::any();
        let ncols: usize = kani::any();
        kani::assume(n >= 1 && n <= 4);
        kani::assume(ncols >= 1 && ncols <= 4);

        let mut buf = [0u8; 16];
        let cols = RawCols::new(&mut buf[..n * ncols], n, ncols);

        let c1: usize = kani::any();
        let s1: usize = kani::any();
        let e1: usize = kani::any();
        kani::assume(c1 < ncols && s1 <= e1 && e1 <= n);

        let c2: usize = kani::any();
        let s2: usize = kani::any();
        let e2: usize = kani::any();
        kani::assume(c2 < ncols && s2 <= e2 && e2 <= n);

        // Disjoint tiles: either different columns or non-overlapping row ranges.
        kani::assume(c1 != c2 || e1 <= s2 || e2 <= s1);

        let slice1 = cols.chunk_mut(c1, s1, e1);
        for v in slice1.iter_mut() {
            *v = 1;
        }
        let slice2 = cols.chunk_mut(c2, s2, e2);
        for v in slice2.iter_mut() {
            *v = 2;
        }

        for r in s1..e1 {
            assert_eq!(buf[c1 * n + r], 1);
        }
        for r in s2..e2 {
            assert_eq!(buf[c2 * n + r], 2);
        }
    }

    #[kani::proof]
    #[kani::unwind(4)]
    #[kani::should_panic]
    fn verify_raw_cols_rejects_overflow() {
        let n: usize = kani::any();
        let ncols: usize = kani::any();
        kani::assume(n.checked_mul(ncols).is_none());
        let mut buf = [0u8; 4];
        let _ = RawCols::new(&mut buf, n, ncols);
    }
}
