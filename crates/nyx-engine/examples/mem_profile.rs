//! Heap use of a solve, stage by stage: a counting allocator records the bytes in use and a watcher thread notes the
//! largest value seen in every stage that the solver reports (see `femcore::progress`).
//!
//!   cargo run --release -p nyx-engine --example mem_profile -- NX NY NZ [p1|p2] [static|eigen] [auto|direct|iter]
#[path = "../../femcore/examples/benchmesh/mod.rs"]
mod benchmesh;
mod common;
use common::*;
use serde_json::Value;
use nyx_engine::facade::Progress;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

struct Counting;
static CURRENT: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = System.alloc(l);
        if !p.is_null() {
            let now = CURRENT.fetch_add(l.size(), Ordering::Relaxed) + l.size();
            PEAK.fetch_max(now, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l);
        CURRENT.fetch_sub(l.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        let q = System.realloc(p, l, new);
        if !q.is_null() {
            if new >= l.size() {
                let now = CURRENT.fetch_add(new - l.size(), Ordering::Relaxed) + new - l.size();
                PEAK.fetch_max(now, Ordering::Relaxed);
            } else {
                CURRENT.fetch_sub(l.size() - new, Ordering::Relaxed);
            }
        }
        q
    }
}

#[global_allocator]
static A: Counting = Counting;

fn gb(b: usize) -> f64 {
    b as f64 / 1e9
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: Vec<usize> = a.iter().take(3).map(|x| x.parse().unwrap()).collect();
    let order: u8 = if a.get(3).map(|x| x.as_str()) == Some("p1") { 1 } else { 2 };
    let eigen = a.get(4).map(|x| x.as_str()) == Some("eigen");
    let linear = match a.get(5).map(|x| x.as_str()) {
        Some("direct") => "direct",
        Some("iter") => "iterative",
        _ => "auto",
    };
    let mesh = benchmesh::block(n[0] as f64, n[1] as f64, n[2] as f64, n[0], n[1], n[2]);
    println!("mesh: {} cells, {:.2} GB in use after meshing", mesh.n_cells(), gb(CURRENT.load(Ordering::Relaxed)));
    let mut s = spec("sms", order, "static", Value::Null);
    steel(&mut s);
    solver(&mut s, linear);
    fixed(&mut s, "xmin");
    load(&mut s, "xmax", [0.0, 0.0, -1e5]);
    if eigen {
        study(&mut s, "eigenfrequency", serde_json::json!({ "n_modes": 6, "f_search": 0.0 }));
    }
    let base = CURRENT.load(Ordering::Relaxed);
    PEAK.store(base, Ordering::Relaxed);
    let p = Progress::new();
    let p2 = p.clone();
    let t0 = Instant::now();
    let h = std::thread::spawn(move || {
        solve(&mesh, &s, Some(p2)).map(|r| (r.dofs, r.log))
    });
    // largest heap use seen per stage, in order of appearance
    let mut stages: Vec<(String, usize)> = vec![];
    while !h.is_finished() {
        let st = p.status();
        let (cur, peak) =
            (CURRENT.load(Ordering::Relaxed), PEAK.swap(CURRENT.load(Ordering::Relaxed), Ordering::Relaxed));
        let hi = cur.max(peak);
        match stages.last_mut() {
            Some((name, m)) if *name == st.stage => *m = (*m).max(hi),
            _ => stages.push((st.stage.clone(), hi)),
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    let r = h.join().unwrap();
    println!(
        "solve: {:.1} s, {}",
        t0.elapsed().as_secs_f64(),
        r.as_ref().map(|r| format!("{} dofs", r.0)).unwrap_or_else(|e| e.clone())
    );
    println!("heap in use before the solve: {:.2} GB (mesh and spec)", gb(base));
    for (name, m) in &stages {
        println!("  peak {:6.2} GB  while: {}", gb(*m), if name.is_empty() { "(between stages)" } else { name });
    }
    if let Ok((dofs, _)) = &r {
        let top = stages.iter().map(|s| s.1).max().unwrap_or(0);
        println!("overall peak {:.2} GB = {:.0} bytes per unknown", gb(top), top as f64 / *dofs as f64);
    }
}
