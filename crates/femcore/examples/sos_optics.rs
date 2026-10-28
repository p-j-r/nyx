//! Optical modes of a silicon strip on sapphire at 1550 nm (450 nm x 230 nm silicon strip, sapphire n = 1.746, silicon n = 3.45, air, PML).
//!   cargo run --release -p femcore --example sos_optics -- [h_fine_nm] [strip_width_nm] [strip_height_nm]
use femcore::mesher::*;
use femcore::problem::*;
use std::time::Instant;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let scale: f64 = std::env::var("SCALE").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
    let h_fine = a.first().copied().unwrap_or(40.0) * 1e-9 * scale;
    let (w_si, h_si) = (a.get(1).copied().unwrap_or(450.0) * 1e-9 * scale, a.get(2).copied().unwrap_or(230.0) * 1e-9 * scale);
    let scale: f64 = std::env::var("SCALE").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0);
    let um = 1e-6 * scale;
    let (w_sapph, h_sapph, h_air, w_pml) = (12.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let lambda = 1550e-9 * scale;
    let c0 = 299792458.0;
    // the physical domain is x in [-6, 6], y in [-5, 8]; PML layers of 1 um around it
    let (x0, x1, y0, y1) = (-w_sapph / 2.0, w_sapph / 2.0, -h_sapph / 2.0, h_sapph / 2.0 + h_air);
    let ys = h_sapph / 2.0; // sapphire surface
    let rect = |xa: f64, ya: f64, xb: f64, yb: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: xa, y0: ya, x1: xb, y1: yb }, domain: d, h };
    let envf = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d);
    let hc = envf("H_COARSE", 0.5) * um; // coarse size in the cladding
    let h_zone = envf("H_ZONE", 0.12) * um;
    let zone = 1.5 * um;
    let parts = vec![
        // domain 1 sapphire, 3 air
        rect(x0 - w_pml, y0 - w_pml, x1 + w_pml, y1 + w_pml, 3, None),
        rect(x0, y0, x1, ys, 1, None),
        // refined zones around the strip
        rect(-zone, ys - zone, zone, ys, 1, Some(h_zone)),
        rect(-zone, ys, zone, ys + zone, 3, Some(h_zone)),
        // the strip
        rect(-w_si / 2.0, ys, w_si / 2.0, ys + h_si, 2, Some(h_fine)),
        // PML layers: 4 left, 5 right, 6 bottom, 7 top, corners 8..11
        rect(x0 - w_pml, y0, x0, y1, 4, None),
        rect(x1, y0, x1 + w_pml, y1, 5, None),
        rect(x0, y0 - w_pml, x1, y0, 6, None),
        rect(x0, y1, x1, y1 + w_pml, 7, None),
        rect(x0 - w_pml, y0 - w_pml, x0, y0, 8, None),
        rect(x1, y0 - w_pml, x1 + w_pml, y0, 9, None),
        rect(x0 - w_pml, y1, x0, y1 + w_pml, 10, None),
        rect(x1, y1, x1 + w_pml, y1 + w_pml, 11, None),
    ];
    // the PML layers sit in sapphire below the surface and air above it: split them at the surface
    let spec2 = Spec2d { parts, h: hc, min_angle_deg: 25.0 };
    let t = Instant::now();
    let mesh = mesh_2d(&spec2).unwrap();
    println!("mesh: {} cells, {} points in {:.2} s", mesh.n_cells(), mesh.n_points(), t.elapsed().as_secs_f64());
    {
        let mut amin = f64::MAX;
        let mut amax = 0.0f64;
        let mut used = vec![false; mesh.n_points()];
        let mut neg = 0;
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let (p0, p1, p2) = (mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize));
            let ar = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
            if ar < 0.0 { neg += 1; }
            if ar.abs() < 1e-20 { println!("degenerate cell {c} label {}: {:?} {:?} {:?} area {ar:e}", mesh.cell_labels[c], p0, p1, p2); }
            amin = amin.min(ar.abs());
            amax = amax.max(ar.abs());
            for &i in v { used[i as usize] = true; }
        }
        println!("area min {:.3e} max {:.3e}, clockwise {}, unused points {}", amin, amax, neg, used.iter().filter(|&&u| !u).count());
    }
    let mut spec = ProblemSpec::new(Physics::Emw);
    spec.order = 2;
    let eps = |n: f64| n * n;
    for d in 1..=11 {
        // PML layers beside the substrate are sapphire up to the surface and air above; keep it simple
        // and use one material per layer (substrate side: sapphire for left/right/bottom, air for top)
        let n = match d {
            1 | 4 | 5 | 6 | 8 | 9 => 1.746,
            2 => 3.45,
            _ => 1.0,
        };
        let n = if std::env::var_os("ALLAIR").is_some() { 1.0 } else if std::env::var_os("NOSI").is_some() && d == 2 { 1.746 } else { n };
        spec.materials.set(d, "eps_r", eps(n));
        spec.materials.set(d, "mu_r", 1.0);
    }
    let pml = |x, y| Pml { x, y, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 };
    spec.materials.pml.insert(4, pml(PmlDir::Minus, PmlDir::Off));
    spec.materials.pml.insert(5, pml(PmlDir::Plus, PmlDir::Off));
    spec.materials.pml.insert(6, pml(PmlDir::Off, PmlDir::Minus));
    spec.materials.pml.insert(7, pml(PmlDir::Off, PmlDir::Plus));
    spec.materials.pml.insert(8, pml(PmlDir::Minus, PmlDir::Minus));
    spec.materials.pml.insert(9, pml(PmlDir::Plus, PmlDir::Minus));
    spec.materials.pml.insert(10, pml(PmlDir::Minus, PmlDir::Plus));
    spec.materials.pml.insert(11, pml(PmlDir::Plus, PmlDir::Plus));
    if std::env::var_os("NOPML").is_some() {
        spec.materials.pml.clear();
    }
    {
        let ext = mesh.exterior_facets();
        println!("facets {} (exterior {}), boundary labels {:?}", mesh.n_facets(), ext.iter().filter(|&&e| e).count(), mesh.boundary_labels());
    }
    spec.study = Study::Modes { freq: c0 / lambda, n_target: std::env::var("NT").ok().and_then(|v| v.parse().ok()).unwrap_or(3.0), evanescent: false, n_modes: 8 };
    let t = Instant::now();
    match femcore::solve::solve(&mesh, &spec, &mut SolveState::default()) {
        Ok(r) => {
            println!("solved {} dofs in {:.2} s", r.dofs, t.elapsed().as_secs_f64());
            for s in &r.solutions {
                println!("  {}", s.label);
            }
        }
        Err(e) => println!("failed: {e}"),
    }
}
