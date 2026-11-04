//! Electromagnetic eigenmodes of 3D cavities.
mod common;
use femcore::props::*;
use nyx_ems::{Problem, Study};

const C0: f64 = 299_792_458.0;

fn cavity_modes(a: f64, b: f64, d: f64, fmax: f64) -> Vec<f64> {
    // PEC box: f = c/2 sqrt((m/a)^2 + (n/b)^2 + (p/d)^2), at least two indices non-zero, none negative
    let mut out = vec![];
    for m in 0..5 {
        for n in 0..5 {
            for p in 0..5 {
                let nonzero = (m > 0) as usize + (n > 0) as usize + (p > 0) as usize;
                if nonzero < 2 {
                    continue;
                }
                let f = 0.5 * C0 * ((m as f64 / a).powi(2) + (n as f64 / b).powi(2) + (p as f64 / d).powi(2)).sqrt();
                if f <= fmax {
                    // TE/TM degeneracy: modes with all three indices non-zero appear twice
                    out.push(f);
                    if nonzero == 3 {
                        out.push(f);
                    }
                }
            }
        }
    }
    out.sort_by(|x, y| x.partial_cmp(y).unwrap());
    out
}

fn spec(order: usize, f_search: f64, n: usize) -> Problem {
    let mut s = common::emw();
    s.order = order;
    s.materials.default.insert("eps_r".into(), 1.0);
    s.materials.default.insert("mu_r".into(), 1.0);
    s.study = Study::Cavity { n_modes: n, f_search };
    s
}

#[test]
fn box_cavity_modes_match_the_analytic_frequencies() {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let mesh = common::block(a, b, d, 5, 4, 3);
    let exact = cavity_modes(a, b, d, 20e9);
    // quadratic Nedelec elements on 60 cells: all of the first six modes within about 1e-3, and no
    // spurious (gradient) modes below the first physical one
    let r = common::solve(&mesh, &spec(2, 0.9 * exact[0], 6)).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|s| s.value.unwrap()).collect();
    for (g, e) in got.iter().zip(&exact) {
        assert!((g / e - 1.0).abs() < 2e-3, "{g} vs {e}: all {got:?} vs {exact:?}");
    }
    assert!(got[0] > 0.95 * exact[0]);
    // a finer mesh improves them
    let fine = common::block(a, b, d, 8, 6, 4);
    let r = common::solve(&fine, &spec(2, 0.9 * exact[0], 6)).unwrap();
    for (s, e) in r.solutions.iter().zip(&exact) {
        assert!((s.value.unwrap() / e - 1.0).abs() < 6e-4, "{} vs {e}", s.value.unwrap());
    }
}

#[test]
fn first_order_elements_converge_more_slowly() {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, d, 20e9);
    let mut errs = vec![];
    for (nx, ny, nz) in [(4, 3, 2), (8, 6, 4)] {
        let mesh = common::block(a, b, d, nx, ny, nz);
        let r = common::solve(&mesh, &spec(1, 0.9 * exact[0], 1)).unwrap();
        errs.push((r.solutions[0].value.unwrap() / exact[0] - 1.0).abs());
    }
    assert!(errs[1] < 0.4 * errs[0], "second-order convergence expected: {errs:?}");
}

#[test]
fn spherical_cavity_lowest_mode_on_a_stuffed_mesh() {
    // TM mode of a PEC sphere: k R = 2.7437 (three-fold degenerate); the mesh follows the sphere to O(h^2). The sphere is a geom model
    // (a sphere in a hole box); the mesh is made by the 3D adapter (ADR 40)
    let r = 0.05;
    let mesh = femcore::meshing::geomdrive3::sphere_cavity(r, 0.012).unwrap().mesh;
    let exact = C0 * 2.7437 / (2.0 * std::f64::consts::PI * r);
    let res = common::solve(&mesh, &spec(2, exact, 4)).unwrap();
    let f: Vec<f64> = res.solutions.iter().map(|s| s.value.unwrap()).collect();
    assert!(f.iter().take(3).all(|v| (v / exact - 1.0).abs() < 0.01), "{f:?} vs {exact}");
    // three modes of the degenerate triplet, then a clearly higher one
    assert!(f[3] > 1.2 * exact || f.len() == 3, "{f:?}");
}

#[test]
fn a_real_stretch_makes_the_cavity_longer() {
    // PML machinery with a real stretch of 2 over the upper half in z equals a cavity 1.5 times as long
    let (a, b, d) = (0.04, 0.03, 0.02);
    let mesh = common::block_split_z(a, b, d, 0.5 * d, d / 6.0, 4);
    let mut s = spec(2, 0.0, 4);
    s.study =
        Study::Cavity { n_modes: 4, f_search: 0.9 * C0 / 2.0 * ((1.0 / a).powi(2) + (1.0 / (1.5 * d)).powi(2)).sqrt() };
    s.materials.set(1, "eps_r", 1.0);
    s.materials.set(2, "eps_r", 1.0);
    s.materials.set(1, "mu_r", 1.0);
    s.materials.set(2, "mu_r", 1.0);
    s.materials
        .pml
        .insert(2, Pml { x: PmlDir::Off, y: PmlDir::Off, z: PmlDir::Plus, strength: 0.0, order: 0.0, re: 2.0 });
    let res = common::solve(&mesh, &s).unwrap();
    let exact = cavity_modes(a, b, 1.5 * d, 20e9);
    let got: Vec<f64> = res.solutions.iter().map(|x| x.value.unwrap()).collect();
    for (g, e) in got.iter().zip(&exact) {
        assert!((g / e - 1.0).abs() < 3e-3, "{got:?} vs {exact:?}");
    }
}

/// Periodic cell of a rectangular PEC waveguide (cross-section a x b in the x-z plane, length l along y, the propagation axis): boundary
/// numbers of `common::block`: 1 x=0, 2 y=0, 3 z=0, 4 z=lz, 5 y=ly, 6 x=lx. Periodic y faces 2 -> 5; `nz` counts the elements along y.
fn guide_cell(order: usize, phase: f64, nx: usize, ny: usize, nz: usize, l: f64) -> Vec<f64> {
    let (a, b) = (0.03, 0.015);
    let mesh = common::block(a, l, b, nx, nz, ny);
    let mut s = spec(order, 0.0, 3);
    s.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase }];
    // TE10: f = c / 2 sqrt((1/a)^2 + (beta/pi)^2), beta = (phase + 2 pi m) / l
    let beta = phase / l;
    let f_te10 = 0.5 * C0 * ((1.0 / a).powi(2) + (beta / std::f64::consts::PI).powi(2)).sqrt();
    s.study = Study::Cavity { n_modes: 3, f_search: 0.95 * f_te10 };
    let r = common::solve(&mesh, &s).unwrap();
    r.solutions.iter().map(|x| x.value.unwrap()).collect()
}

#[test]
fn floquet_periodic_waveguide_cell_gives_the_guided_dispersion() {
    let l = 0.01;
    for phase in [0.8f64, 2.0] {
        let beta = phase / l;
        let a: f64 = 0.03;
        let f_te10 = 0.5 * C0 * ((1.0 / a).powi(2) + (beta / std::f64::consts::PI).powi(2)).sqrt();
        let got = guide_cell(2, phase, 4, 2, 2, l);
        assert!(got.iter().any(|f| (f / f_te10 - 1.0).abs() < 3e-3), "phase {phase}: {got:?} vs TE10 {f_te10}");
    }
}

#[test]
fn cavity_modes_carry_a_polarization_analysis() {
    // the lowest mode of the 40 x 30 x 20 mm box has E along z only (the surface normal of the wafer frame: quasi-TM), the second one (indices
    // 1, 0, 1) has E along y only (the propagation axis: all of its field is "longitudinal"); the lateral axis is x
    let (a, b, d) = (0.04, 0.03, 0.02);
    let mesh = common::block(a, b, d, 5, 4, 3);
    let exact = cavity_modes(a, b, d, 20e9);
    let solved = common::try_solve_analysed(&mesh, &spec(2, 0.9 * exact[0], 3)).unwrap();
    let (r, an) = (&solved.result, solved.analysis.as_ref().expect("a polarization analysis"));
    assert_eq!(an.value_label, "f (THz)");
    let zero = an.rows[0].report.as_ref().unwrap();
    assert!(zero.sphere[0] < -0.99 && zero.te_fraction < 0.01 && zero.long_fraction < 0.01, "{zero:?}");
    assert_eq!(zero.class.label(), "quasi-TM (z)");
    let first = an.rows[1].report.as_ref().unwrap();
    assert!(first.long_fraction > 0.99, "{first:?}");
    // a standing wave of a closed cavity has real fields: no handedness
    assert!(first.sphere[2].abs() < 1e-6);
    assert!(r.log.iter().any(|l| l.contains("polarization axes: lateral x, normal z, propagation y")));
}
