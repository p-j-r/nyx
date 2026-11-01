//! Acoustic modes of the silicon strip on sapphire with
//! an out-of-plane wave number: Si strip 450 nm x 230 nm on a sapphire block with fixed outer
//! sides, rotated sapphire (Euler 60, 57.6, 0 degrees), forward Brillouin scattering
//! kz = neff 2 pi f_ac / c (about 460 rad/m); the search is around 9.47 GHz.
//!   cargo run --release -p femcore --example sos_acoustics -- [kz rad/m] [f_search GHz] [n_modes]
use femcore::material::*;
use femcore::mesher::*;
use femcore::problem::*;
use std::time::Instant;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let kz = a.first().copied().unwrap_or(2.3172 * 2.0 * std::f64::consts::PI * 9.47e9 / 299792458.0);
    let f_search = a.get(1).copied().unwrap_or(9.47) * 1e9;
    let n_modes = a.get(2).copied().unwrap_or(12.0) as usize;
    let envf = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(d);
    let um = 1e-6;
    let (w_sapph, h_sapph, h_air, w_pml) = (12.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let (w_si, h_si) = (450e-9, 230e-9);
    let (x0, x1, y0, y1) = (-w_sapph / 2.0, w_sapph / 2.0, -h_sapph / 2.0, h_sapph / 2.0 + h_air);
    let ys = h_sapph / 2.0;
    let rect = |xa: f64, ya: f64, xb: f64, yb: f64, d: i32, h: Option<f64>| Part { shape: Shape::Rect { x0: xa, y0: ya, x1: xb, y1: yb }, domain: d, h };
    let hc = envf("H_COARSE", 0.4) * um;
    let zone = 1.5 * um;
    let parts = vec![
        rect(x0 - w_pml, y0 - w_pml, x1 + w_pml, y1 + w_pml, 3, None),
        rect(x0, y0, x1, ys, 1, None),
        rect(-zone, ys - zone, zone, ys, 1, Some(envf("H_ZONE", 0.08) * um)),
        rect(-zone, ys, zone, ys + zone, 3, Some(envf("H_ZONE", 0.08) * um)),
        rect(-w_si / 2.0, ys, w_si / 2.0, ys + h_si, 2, Some(envf("H_FINE", 25.0) * 1e-9)),
        rect(x0 - w_pml, y0, x0, ys, 4, None),
        rect(x1, y0, x1 + w_pml, ys, 5, None),
        rect(x0, y0 - w_pml, x1, y0, 6, None),
        rect(x0 - w_pml, y0 - w_pml, x0, y0, 8, None),
        rect(x1, y0 - w_pml, x1 + w_pml, y0, 9, None),
    ];
    let mesh = mesh_2d(&Spec2d { parts, h: hc, min_angle_deg: 25.0 }).unwrap();
    println!("mesh: {} cells, {} points", mesh.n_cells(), mesh.n_points());
    // the solid occupies the sapphire (domains 1, 4, 5, 6, 8, 9) and the silicon strip (2); air is inactive.
    // The mesh also contains the air-side PML domains of the optical model (3 above); every domain with no
    // solid is listed as inactive.
    let mut spec = ProblemSpec::new(Physics::Solid);
    spec.order = 2;
    let domains = mesh.domain_labels();
    spec.inactive = domains.iter().copied().filter(|d| ![1, 2, 4, 5, 6, 8, 9].contains(d)).collect();
    for &d in &domains {
        spec.materials.set(d, "E", 100e9);
        spec.materials.set(d, "nu", 0.25);
        spec.materials.set(d, "rho", 3980.0);
    }
    let eul = std::env::var("EULER").ok().map(|s| { let v: Vec<f64> = s.split(',').filter_map(|x| x.parse().ok()).collect(); [v[0], v[1], v[2]] }).unwrap_or([60.0, 57.6, 0.0]);
    let saph = Anisotropy { stiffness: Some(Stiffness::Trigonal { c11: 497.73e9, c12: 162.8e9, c13: 116.0e9, c14: -21.9e9, c33: 500.9e9, c44: 146.8e9 }), euler_deg: eul, ..Default::default() };
    for d in [1, 4, 5, 6, 8, 9] {
        spec.materials.aniso.insert(d, saph.clone());
        spec.materials.set(d, "rho", 3980.0);
    }
    spec.materials.aniso.insert(2, Anisotropy { stiffness: Some(Stiffness::Cubic { c11: 166e9, c12: 64e9, c44: 79e9 }), euler_deg: [0.0, 0.0, 0.0], ..Default::default() });
    spec.materials.set(2, "rho", 2329.0);
    // fixed: the outer sides and the bottom of the sapphire block
    let mut labels = vec![];
    for f in 0..mesh.n_facets() {
        let (pa, pb) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
        let (mx, my) = (0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]));
        let outer = (mx - (x0 - w_pml)).abs() < 1e-9 || (mx - (x1 + w_pml)).abs() < 1e-9 || (my - (y0 - w_pml)).abs() < 1e-9;
        if outer && my < ys + 1e-9 {
            labels.push(mesh.facet_labels[f]);
        }
    }
    labels.sort();
    labels.dedup();
    spec.bcs = vec![Bc { on: None, labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }];
    if std::env::var_os("PML").is_some() {
        let st = envf("PML_STRENGTH", 3.0);
        let pml = |x, y| Pml { x, y, z: PmlDir::Off, strength: st, order: 2.0, re: 1.0 };
        spec.materials.pml.insert(4, pml(PmlDir::Minus, PmlDir::Off));
        spec.materials.pml.insert(5, pml(PmlDir::Plus, PmlDir::Off));
        spec.materials.pml.insert(6, pml(PmlDir::Off, PmlDir::Minus));
        spec.materials.pml.insert(8, pml(PmlDir::Minus, PmlDir::Minus));
        spec.materials.pml.insert(9, pml(PmlDir::Plus, PmlDir::Minus));
    }
    spec.study = Study::Eigenfrequency { n_modes, f_search };
    spec.out_of_plane = Some(kz);
    println!("kz = {kz:.3} rad/m, searching {:.3} GHz", f_search / 1e9);
    let t = Instant::now();
    match femcore::solve::solve(&mesh, &spec, &mut SolveState::default()) {
        Ok(r) => {
            println!("solved {} dofs in {:.1} s; {}", r.dofs, t.elapsed().as_secs_f64(), r.log.join(" | "));
            for s in &r.solutions {
                println!("  {}", s.label);
            }
        }
        Err(e) => println!("failed: {e}"),
    }
}
