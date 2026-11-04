#![allow(clippy::needless_range_loop)]
//! M1 benchmark: fused derivative-to-metric kernel vs the straightforward multi-pass version.
//!
//! Minimum of repeated runs after warm-up, variants interleaved in one session.
//! Usage: `cargo run --release --example bench_metric -- [n_vertices=200000] [repeats=7]`

use femcore::adapt::metric::*;
use std::time::Instant;

thread_local! {
    static NOTES: std::cell::RefCell<Vec<String>> = const { std::cell::RefCell::new(Vec::new()) };
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

fn field<const D: usize>(n: usize) -> (Vec<Sym<D>>, Vec<f64>) {
    let mut g = Lcg(42);
    let hs = (0..n)
        .map(|_| {
            let mut a = [[0.0; D]; D];
            for i in 0..D {
                for j in i..D {
                    let v = (g.next() - 0.5) * 10f64.powf(g.next() * 4.0 - 1.0);
                    a[i][j] = v;
                    a[j][i] = v;
                }
            }
            a
        })
        .collect();
    (hs, vec![1.0 / n as f64; n])
}

fn min_time<R>(reps: usize, mut f: impl FnMut() -> R) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..reps {
        let t = Instant::now();
        std::hint::black_box(f());
        best = best.min(t.elapsed().as_secs_f64());
    }
    best
}

fn bench<const D: usize>(n: usize, reps: usize, nt: usize) {
    let (hs, vol) = field::<D>(n);
    let opt = LpOptions { p: 2.0, target_complexity: n as f64 / 50.0, hmin: 1e-3, hmax: 10.0, ..Default::default() };
    // warm-up
    std::hint::black_box(lp_metrics_with_threads(&hs, &vol, &opt, 1));
    std::hint::black_box(lp_metrics_unfused(&hs, &vol, &opt));
    let (mut tf1, mut tfn, mut tu) = (f64::INFINITY, f64::INFINITY, f64::INFINITY);
    for _ in 0..reps {
        tu = tu.min(min_time(1, || lp_metrics_unfused(&hs, &vol, &opt)));
        tf1 = tf1.min(min_time(1, || lp_metrics_with_threads(&hs, &vol, &opt, 1)));
        tfn = tfn.min(min_time(1, || lp_metrics_with_threads(&hs, &vol, &opt, nt)));
    }
    let r1 = lp_metrics_with_threads(&hs, &vol, &opt, 1);
    let rn = lp_metrics_with_threads(&hs, &vol, &opt, nt);
    let same = r1
        .metrics
        .iter()
        .zip(&rn.metrics)
        .all(|(a, b)| (0..D).all(|i| (0..D).all(|j| a[i][j].to_bits() == b[i][j].to_bits())));
    let un = lp_metrics_unfused(&hs, &vol, &opt);
    let t_eig = min_time(reps, || hs.iter().map(|h| sym_eig(h).0[0]).sum::<f64>());
    NOTES.with(|v| {
        v.borrow_mut().push(format!(
            "{D}D: {n} eigen-decompositions alone take {:.1} ms = {:.0}% of the fused 1-thread time",
            t_eig * 1e3,
            100.0 * t_eig / tf1
        ))
    });
    println!(
        "| {D}D | {n} | {:.1} | {:.1} | {:.1} | {:.2}x | {:.2}x | {} | {:.1e} | {:.1e} |",
        tu * 1e3,
        tf1 * 1e3,
        tfn * 1e3,
        tu / tf1,
        tf1 / tfn,
        if same { "yes" } else { "NO" },
        (r1.complexity / opt.target_complexity - 1.0).abs(),
        (complexity(&un, &vol) / opt.target_complexity - 1.0).abs(),
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(200_000);
    let reps: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(7);
    let nt = femcore::adapt::par::threads();
    println!("# lp_metrics: fused vs unfused ({} threads for the parallel column)\n", nt);
    println!(
        "| dim | vertices | unfused ms | fused 1T ms | fused {nt}T ms | fused/unfused speedup | thread speedup | bit-identical across threads | fused complexity miss | unfused complexity miss |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    bench::<2>(n, reps, nt);
    bench::<3>(n, reps, nt);
    println!();
    NOTES.with(|v| v.borrow().iter().for_each(|l| println!("* {l}")));
}
