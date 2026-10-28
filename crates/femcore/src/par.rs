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
}
