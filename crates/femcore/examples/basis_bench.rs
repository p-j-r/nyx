//! Memory-bandwidth check of the orthogonalisation kernels (the `BV` role).
//!
//!   cargo run --release --example basis_bench -- [n]
//!
//! For each scalar type and basis width it reports the effective bandwidth of one pass over the
//! basis (`bytes of V / time`) for each kernel, next to a plain streaming sum of the same data.
//! A kernel that is memory-bound should run at or near the streaming figure; if not, the code
//! generation (not the algorithm) is the problem. Pitfall 1 of the handoff brief: always look.

use femcore::linalg::basis::Basis;
use femcore::linalg::rng::random_vector;
use femcore::linalg::{C64, Scalar};
use std::hint::black_box;
use std::time::Instant;

fn best<F: FnMut()>(mut f: F, reps: usize) -> f64 {
    f();
    (0..reps)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::MAX, f64::min)
}

fn stream_sum(data: &[f64]) -> f64 {
    let mut acc = [0.0f64; 8];
    for c in data.chunks_exact(8) {
        for l in 0..8 {
            acc[l] += c[l];
        }
    }
    acc.iter().sum()
}

fn run<T: Scalar>(name: &str, n: usize, k: usize, nb: usize) {
    let mut v = Basis::<T>::new(n, k + nb + 1);
    for j in 0..k + nb {
        v.col_mut(j).copy_from_slice(&random_vector::<T>(n, j as u64 + 1));
    }
    let bytes = (n * k * std::mem::size_of::<T>()) as f64;
    // plain streaming reference over the same amount of memory, viewed as f64
    let flat: Vec<f64> = (0..(n * k * std::mem::size_of::<T>() / 8)).map(|i| (i % 7) as f64).collect();
    let t_stream = best(
        || {
            black_box(stream_sum(black_box(&flat)));
        },
        8,
    );

    let mut h = vec![T::ZERO; k * nb];
    let mut h2 = vec![T::ZERO; k * nb];
    let mut nrm = vec![0.0; nb];
    let mut w = random_vector::<T>(n * nb, 99);
    let (head, _tail) = v.split_cols_mut(k);
    let t_dots = best(|| head.dots_block(k, black_box(&w), nb, &mut h, Some(&mut nrm)), 8);
    let t_sub = best(|| head.sub_block(k, &h, &mut w, nb), 8);
    let t_fused = best(|| head.sub_block_then_dots(k, &h, &mut w, nb, &mut h2, &mut nrm), 8);
    // restart multiply: V[:, 0..k/2] <- V[:, 0..k] * Q  (reads k columns, writes k/2)
    let q = random_vector::<T>(k * (k / 2), 5);
    let t_mult = best(|| v.mult_in_place(0, k, k / 2, &q, k), 4);
    let gb = |t: f64| bytes / t / 1e9;
    println!(
        "{name:<8} n={n:>8} k={k:>3} b={nb}: stream {:5.1} GB/s | dots {:5.1} | sub {:5.1} | sub+dots(fused) {:5.1} (vs 2 passes) | V*Q {:5.1} GB/s (read side)",
        gb(t_stream),
        gb(t_dots),
        gb(t_sub),
        gb(t_fused),
        gb(t_mult)
    );
}

fn main() {
    let n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(400_000);
    for &(k, nb) in &[(20usize, 1usize), (60, 1), (60, 2), (60, 4)] {
        run::<f64>("f64", n, k, nb);
        run::<C64>("c64", n / 2, k, nb);
    }
}
