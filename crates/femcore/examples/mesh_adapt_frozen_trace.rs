use femcore::meshers::adapt::drive::AdaptOptions;
use femcore::meshers::adapt::geom::PolylineGeometry;
use femcore::meshers::adapt::metric::*;
use femcore::meshers::adapt::ops::*;
use femcore::meshers::adapt::topo::{TopoOptions, Work2};
use femcore::meshers::fixtures::{box_tri, side_label};

fn main() {
    let m = box_tri([[0.0, 1.0], [0.0, 1.0]], 5, 5);
    let f = |_: &[f64; 2]| isotropic::<2>(0.04);
    let ao = AdaptOptions {
        topo: TopoOptions { frozen_labels: vec![side_label(0, false)], ..Default::default() },
        ..Default::default()
    };
    let mut w = Work2::from_mesh(&m, &f, &ao.topo).unwrap();
    let geom = PolylineGeometry::from_work(&w);
    let mut step = 0;
    for it in 0..40 {
        for _ in 0..30 {
            let n = split_pass(&mut w, &f, &geom, &ao, it < 2);
            step += 1;
            if let Err(e) = w.check() {
                println!("step {step} split n={n}: {e}");
                return;
            }
            if n == 0 {
                break;
            }
            let n = flip_pass(&mut w);
            step += 1;
            if let Err(e) = w.check() {
                println!("step {step} flip n={n}: {e}");
                return;
            }
        }
        for _ in 0..10 {
            let n = collapse_pass(&mut w, &ao);
            step += 1;
            if let Err(e) = w.check() {
                println!("step {step} collapse n={n}: {e}");
                return;
            }
            if n == 0 {
                break;
            }
        }
        for _ in 0..3 {
            let n = smooth_pass(&mut w, &f, &geom, &ao);
            step += 1;
            if let Err(e) = w.check() {
                println!("step {step} smooth n={n}: {e}");
                return;
            }
        }
        println!("iteration {it} ok, {} alive tris", w.n_alive());
    }
}
