//! Deterministic chunked parallelism on `std::thread::scope`.
//!
//! Work is cut into fixed-size chunks that do not depend on the thread count, and chunk results
//! are combined in chunk order. A result is therefore bit-identical for any number of threads
//! (the same contract as the fixed-chunk reductions in the solver port). The thread count comes
//! from `FEMMESH_THREADS` (default: available parallelism). No pool is kept here: the calls are
//! coarse enough that spawning scoped threads is negligible, and nyx's own pool replaces this
//! module when the code is moved over.

/// Chunk length used by default.
pub const CHUNK: usize = 512;

/// Number of worker threads.
pub fn threads() -> usize {
    std::env::var("FEMMESH_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&t| t >= 1)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1))
}

/// Apply `f(start, chunk)` to consecutive fixed-size chunks of `data`, in parallel, and return the
/// sum of the returned values taken in chunk order.
///
/// `f` must only touch its own chunk. With `nthreads == 1`, or fewer than two chunks, everything
/// runs on the calling thread with the same chunking.
pub fn fold_chunks_mut<T: Send, F>(data: &mut [T], chunk: usize, nthreads: usize, f: F) -> f64
where
    F: Fn(usize, &mut [T]) -> f64 + Sync,
{
    let n = data.len();
    let nchunks = n.div_ceil(chunk.max(1));
    let mut partial = vec![0.0f64; nchunks];
    if nthreads <= 1 || nchunks < 2 {
        for (ci, (d, p)) in data.chunks_mut(chunk).zip(partial.iter_mut()).enumerate() {
            *p = f(ci * chunk, d);
        }
    } else {
        let nt = nthreads.min(nchunks);
        // static round-robin assignment of (chunk index, chunk, partial slot) triples
        let mut work: Vec<Vec<(usize, &mut [T], &mut f64)>> = (0..nt).map(|_| Vec::new()).collect();
        for (ci, (d, p)) in data.chunks_mut(chunk).zip(partial.iter_mut()).enumerate() {
            work[ci % nt].push((ci, d, p));
        }
        std::thread::scope(|s| {
            let f = &f;
            for w in work {
                s.spawn(move || {
                    for (ci, d, p) in w {
                        *p = f(ci * chunk, d);
                    }
                });
            }
        });
    }
    partial.iter().sum()
}

/// Sum of `f(i)` over `0..n`, evaluated in fixed chunks and combined in chunk order.
pub fn sum_range<F>(n: usize, chunk: usize, nthreads: usize, f: F) -> f64
where
    F: Fn(usize) -> f64 + Sync,
{
    let mut idx = vec![(); n];
    fold_chunks_mut(&mut idx, chunk, nthreads, |start, c| {
        let mut s = 0.0;
        for k in 0..c.len() {
            s += f(start + k);
        }
        s
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sums_are_bit_identical_across_thread_counts() {
        let n = 100_003;
        let f = |i: usize| ((i as f64) * 0.618).sin() * 1e-3 + 1.0 / (1.0 + i as f64);
        let r1 = sum_range(n, CHUNK, 1, f);
        for t in [2, 3, 4, 7, 16] {
            assert_eq!(r1.to_bits(), sum_range(n, CHUNK, t, f).to_bits(), "threads {t}");
        }
    }

    #[test]
    fn chunk_start_offsets_are_right() {
        let mut v = vec![0usize; 5000];
        fold_chunks_mut(&mut v, 128, 4, |start, c| {
            for (k, x) in c.iter_mut().enumerate() {
                *x = start + k;
            }
            0.0
        });
        assert!(v.iter().enumerate().all(|(i, &x)| i == x));
    }
}
