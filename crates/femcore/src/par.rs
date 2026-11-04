//! The one shared pool of worker threads used by every parallel loop of the solvers.
//!
//! Before, each parallel kernel (sparse matrix-vector product, orthogonalisation, assembly, ...) started its own
//! operating-system threads with `std::thread::scope` and let them die again; starting 20 threads costs about
//! 0.5 ms, as much as a whole matrix-vector product of a medium problem, and a multigrid cycle does dozens of them.
//! Now all kernels run on a rayon pool that is created once and kept (fork-join of 20 tasks: ~12 microseconds).
//!
//! The number of threads is the number of logical processors (at most 32) unless the environment variable
//! `NYX_THREADS` says otherwise, for example `NYX_THREADS=4` to leave cores to other programs or
//! `NYX_THREADS=1` for a fully serial run. It is read once, when the first parallel kernel runs.
//!
//! Work is always split into the same pieces for a given thread count and the pieces are combined in a fixed
//! order, so a run gives the same numbers every time (different thread counts can differ in the last bits of sums).

use std::sync::OnceLock;

/// Largest number of threads ever used.
const MAX_THREADS: usize = 32;

/// Number of worker threads of the shared pool (and the number of pieces large loops are split into).
pub fn threads() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| {
        let n = std::env::var("NYX_THREADS").ok().and_then(|v| v.trim().parse::<usize>().ok()).filter(|&n| n >= 1).unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1));
        let n = n.min(MAX_THREADS);
        // the global pool can only be built once; if something else built it first, rayon's own size applies
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n).thread_name(|i| format!("nyx-pool-{i}")).build_global();
        n
    })
}

/// `threads()` capped at `cap` (some kernels stop scaling earlier).
pub fn threads_up_to(cap: usize) -> usize {
    threads().min(cap).max(1)
}

/// Threads used for sparse factorisations (supernodal Cholesky/LU of faer). They scale only to about 8 threads and
/// get slower again beyond: measured on a 6 P-core + 8 E-core (20 logical) processor, 160,000-unknown complex LU of a
/// 25-point 2D stencil: 1 thread 200 s, 4: 65 s, 6: 46 s, 8: 38 s, 12: 47 s, 20: 54 s; 3D Cholesky of 125,000 unknowns:
/// 9.6 s, 3.4 s, 2.8 s, 3.1 s, 3.7 s, 5.3 s. So factorisations run on their own smaller pool; `NYX_FACTOR_THREADS=N`
/// overrides the default of `min(threads, 8)`.
pub fn factor_threads() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| std::env::var("NYX_FACTOR_THREADS").ok().and_then(|v| v.trim().parse::<usize>().ok()).filter(|&n| n >= 1).unwrap_or(8).min(threads()))
}

/// Runs `f` (a sparse factorisation) on the factorisation pool: faer takes the number of threads of the pool it is called in.
pub fn on_factor_pool<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    static POOL: OnceLock<Option<rayon::ThreadPool>> = OnceLock::new();
    let n = factor_threads();
    if n >= threads() {
        return f();
    }
    match POOL.get_or_init(|| rayon::ThreadPoolBuilder::new().num_threads(n).thread_name(|i| format!("nyx-factor-{i}")).build().ok()) {
        Some(pool) => pool.install(f),
        None => f(),
    }
}


/// Chunk length used by default.
pub const CHUNK: usize = 512;

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
    use rayon::prelude::*;

    #[test]
    fn the_pool_has_the_requested_size_and_runs_work() {
        let n = threads();
        assert!((1..=MAX_THREADS).contains(&n));
        let sum: usize = (0..1000usize).into_par_iter().map(|i| i * 2).sum();
        assert_eq!(sum, 999 * 1000);
        assert_eq!(threads_up_to(1), 1);
    }

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
