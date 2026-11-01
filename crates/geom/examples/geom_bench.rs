//! Timings for the targets of the brief (section 10): a ~50-feature model, selector resolution over ~500 facets, a parameter sweep
//! showing cache hits, and a 200-hole pattern. `cargo run --release --example geom_bench`.

use geom::eval::{Evaluator, Tier};
use geom::model::*;
use geom::select::Selector;
use std::time::Instant;

fn time<T>(label: &str, reps: usize, mut f: impl FnMut() -> T) -> f64 {
    let t = Instant::now();
    for _ in 0..reps {
        std::hint::black_box(f());
    }
    let ms = t.elapsed().as_secs_f64() * 1e3 / reps as f64;
    println!("{label:<58} {ms:>9.3} ms");
    ms
}

/// 48 features: a clad/box stack, a 5x5 photonic-crystal-like grid of circles and rounded rectangles, a rib with sloped walls.
fn model_50(w: f64) -> FeatureList {
    let mut m = Model::new_2d();
    m.param("w", w).param("t", 300e-9).param("r", 80e-9);
    let clad = m.rect("clad", ("-4e-6", "-2e-6"), ("4e-6", "t+1e-6"));
    let bx = m.rect("box", ("-4e-6", "-2e-6"), ("4e-6", 0.0));
    let core = m.rect("core", ("-w/2", 0.0), ("w/2", "t"));
    let mut regs = vec![("clad".to_string(), clad), ("BOX".to_string(), bx), ("Si".to_string(), core)];
    for i in 0..5 {
        for j in 0..5 {
            let c = m.circle(
                &format!("c{i}_{j}"),
                (format!("{}e-6", -1.6 + 0.8 * i as f64), format!("{}e-7", 0.5 + 0.9 * j as f64)).into_exprs(),
                "r",
            );
            regs.push((format!("hole{i}_{j}"), c));
        }
    }
    for k in 0..20 {
        let rr = m.rounded_rect(
            &format!("rr{k}"),
            (format!("{}e-6", -3.8 + 0.37 * k as f64), "-1.8e-6".to_string()).into_exprs(),
            (format!("{}e-6", -3.7 + 0.37 * k as f64), "-1.7e-6".to_string()).into_exprs(),
            "10e-9",
        );
        regs.push((format!("rr{k}"), rr));
    }
    let regs: Vec<(&str, FeatureId)> = regs.iter().map(|(n, i)| (n.as_str(), *i)).collect();
    m.overlay(&regs);
    m.list
}

trait IntoExprs {
    fn into_exprs(self) -> (geom::expr::Expr, geom::expr::Expr);
}
impl IntoExprs for (String, String) {
    fn into_exprs(self) -> (geom::expr::Expr, geom::expr::Expr) {
        (self.0.into(), self.1.into())
    }
}

fn main() {
    let l = model_50(500e-9);
    println!("features: {}", l.features.len());
    let mut ev = Evaluator::new();
    let cold =
        time("evaluate ~50 features (cold cache each time)", 20, || Evaluator::new().eval(&l, Tier::Complex).unwrap());
    let e = ev.eval(&l, Tier::Complex).unwrap();
    println!("  -> {} regions, {} edges", e.regions().len(), e.edges().len());
    time("evaluate ~50 features (warm cache, same input)", 20, || ev.eval(&l, Tier::Complex).unwrap());
    time("Classify tier only", 200, || Evaluator::new().eval(&l, Tier::Classify).unwrap());
    let before = ev.stats();
    time("sweep of 100 widths (shared cache)", 1, || {
        for k in 0..100 {
            ev.eval(&l.with(&[("w", 400e-9 + 2e-9 * k as f64)]), Tier::Complex).unwrap();
        }
    });
    let s = ev.stats();
    println!(
        "  -> prims built {} cached {}; pairs computed {} cached {}",
        s.prims_built - before.prims_built,
        s.prims_cached - before.prims_cached,
        s.pairs_computed - before.pairs_computed,
        s.pairs_cached - before.pairs_cached
    );
    let sel =
        Selector::parse("boundary(region(Si)) | interface(BOX, *) | (outer & extreme(+y)) | kind(circle)").unwrap();
    time(&format!("resolve a 4-term selector over {} edges", e.edges().len()), 200, || e.resolve(&sel).unwrap());
    let sel2 = Selector::parse("name(BOX|Si)").unwrap();
    time("resolve name(...)", 2000, || e.resolve(&sel2).unwrap());
    time("label map", 200, || e.label_map().unwrap());
    time("tessellate (max_dev 1 nm)", 20, || e.tessellate(1e-9));

    let mut m = Model::new_2d();
    let plate = m.rect("plate", (0.0, 0.0), (210.0, 10.0));
    let h = m.circle("hole", (1.0, 5.0), 0.3);
    let hs = m.linear_pattern("holes", h, (1.0, 0.0, 0.0), 200.0);
    let d = m.difference("pd", plate, &[hs]);
    m.overlay(&[("plate", d)]);
    let e = Evaluator::new().eval(&m.list, Tier::Complex).unwrap();
    println!("200-hole pattern: {} edges", e.edges().len());
    let _ = cold;

    // 3D: one Bragg period (two layers), and a slab with a 10 x 10 grid of cylindrical holes (a photonic-crystal slab)
    let mut m = Model::new_3d();
    m.param("w", 500e-9).param("t", 300e-9).param("p", 280e-9).param("tooth", 80e-9);
    let clad = m.box_("clad", ("-1.2e-6", "-2e-6", 0.0), ("1.2e-6", "1.3e-6", "p"));
    let core = m.box_("core", ("-w/2", 0.0, 0.0), ("w/2", "t", "p"));
    let lt = m.box_("lt", ("-w/2-tooth", 0.0, 0.0), ("-w/2", "t", "p/2"));
    let rt = m.box_("rt", ("w/2", 0.0, 0.0), ("w/2+tooth", "t", "p/2"));
    m.overlay(&[("clad", clad), ("Si", core), ("lt", lt), ("rt", rt)]);
    let l3 = m.list.clone();
    time("evaluate one Bragg period (3D, 2 layers)", 50, || Evaluator::new().eval(&l3, Tier::Complex).unwrap());
    let mut m = Model::new_3d();
    let slab = m.box_("slab", (0.0, 0.0, 0.0), (11.0, 11.0, 0.22));
    let hole = m.cylinder("hole", (0.5, 0.5, -0.1), 0.2, 0.42);
    let row = m.linear_pattern("row", hole, (1.0, 0.0, 0.0), 10.0);
    let grid = m.linear_pattern("grid", row, (0.0, 1.0, 0.0), 10.0);
    let d = m.difference("pc", slab, &[grid]);
    m.overlay(&[("slab", d)]);
    let pc = m.list.clone();
    let e = Evaluator::new().eval(&pc, Tier::Complex).unwrap();
    println!(
        "photonic-crystal slab: {} faces, {} overlay edges",
        e.faces().len(),
        e.layer_stack().unwrap().overlay.edges.len()
    );
    time("evaluate the 100-hole photonic-crystal slab (3D)", 5, || Evaluator::new().eval(&pc, Tier::Complex).unwrap());
    time("tessellate it (max_dev 1e-3)", 5, || e.tessellate(1e-3));
    time("evaluate the 200-hole plate (cold)", 10, || Evaluator::new().eval(&m.list, Tier::Complex).unwrap());
}
