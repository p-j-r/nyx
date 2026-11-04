//! Nanophotonic scale: a 450 nm x 230 nm silicon strip (n = 3.45) on sapphire (n = 1.746) under air at 1550 nm with PML on
//! all sides. The guided mode must lie between the cladding and core indices, be lossless, and be converged under mesh
//! refinement (the finest meshes of the convergence study give neff = 2.3175).
mod common;
use femcore::props::*;
use nyx_ems::{Problem, Study};

fn sos_mesh(h_coarse: f64, h_zone: f64, h_fine: f64) -> femcore::mesh::SimplexMesh {
    let um = 1e-6;
    let (w_sapph, h_sapph, h_air, w_pml) = (12.0 * um, 10.0 * um, 3.0 * um, 1.0 * um);
    let (w_si, h_si) = (450e-9, 230e-9);
    let (x0, x1, y0, y1) = (-w_sapph / 2.0, w_sapph / 2.0, -h_sapph / 2.0, h_sapph / 2.0 + h_air);
    let ys = h_sapph / 2.0;
    let rect = common::part;
    let zone = 1.5 * um;
    let parts = [
        rect(x0 - w_pml, y0 - w_pml, x1 + w_pml, y1 + w_pml, 3, None),
        rect(x0, y0, x1, ys, 1, None),
        rect(-zone, ys - zone, zone, ys, 1, Some(h_zone)),
        rect(-zone, ys, zone, ys + zone, 3, Some(h_zone)),
        rect(-w_si / 2.0, ys, w_si / 2.0, ys + h_si, 2, Some(h_fine)),
        rect(x0 - w_pml, y0, x0, y1, 4, None),
        rect(x1, y0, x1 + w_pml, y1, 5, None),
        rect(x0, y0 - w_pml, x1, y0, 6, None),
        rect(x0, y1, x1, y1 + w_pml, 7, None),
        rect(x0 - w_pml, y0 - w_pml, x0, y0, 8, None),
        rect(x1, y0 - w_pml, x1 + w_pml, y0, 9, None),
        rect(x0 - w_pml, y1, x0, y1 + w_pml, 10, None),
        rect(x1, y1, x1 + w_pml, y1 + w_pml, 11, None),
    ];
    common::parts_mesh(&parts, h_coarse)
}

fn optical_spec() -> Problem {
    let mut spec = common::emw();
    spec.order = 2;
    for d in 1..=11 {
        let n: f64 = match d {
            1 | 4 | 5 | 6 | 8 | 9 => 1.746,
            2 => 3.45,
            _ => 1.0,
        };
        spec.materials.set(d, "eps_r", n * n);
        spec.materials.set(d, "mu_r", 1.0);
    }
    let pml = |x, z| Pml { x, y: PmlDir::Off, z, strength: 3.0, order: 2.0, re: 1.0 };
    for (d, p) in [
        (4, pml(PmlDir::Minus, PmlDir::Off)),
        (5, pml(PmlDir::Plus, PmlDir::Off)),
        (6, pml(PmlDir::Off, PmlDir::Minus)),
        (7, pml(PmlDir::Off, PmlDir::Plus)),
        (8, pml(PmlDir::Minus, PmlDir::Minus)),
        (9, pml(PmlDir::Plus, PmlDir::Minus)),
        (10, pml(PmlDir::Minus, PmlDir::Plus)),
        (11, pml(PmlDir::Plus, PmlDir::Plus)),
    ] {
        spec.materials.pml.insert(d, p);
    }
    spec.study = Study::Modes { freq: 299792458.0 / 1550e-9, n_target: 2.4, evanescent: false, n_modes: 4 };
    spec
}

#[test]
fn a_mesh_of_overlapping_rectangles_has_no_degenerate_cells() {
    // the layers touch the substrate and the strip at coincident edges; rounding must not leave slivers
    let mesh = sos_mesh(0.5e-6, 0.12e-6, 40e-9);
    for c in 0..mesh.n_cells() {
        let v = mesh.cell(c);
        let (p0, p1, p2) = (mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize));
        let area = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])).abs();
        let longest2 = (0..3)
            .map(|i| {
                let (a, b) = (mesh.point(v[i] as usize), mesh.point(v[(i + 1) % 3] as usize));
                (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)
            })
            .fold(0.0, f64::max);
        assert!(
            area > 1e-3 * longest2.max(1e-300) * 1e-3,
            "cell {c} is degenerate: area {area:e}, longest edge^2 {longest2:e}"
        );
    }
}

#[test]
fn silicon_strip_on_sapphire_guided_mode_is_converged() {
    let neff = |h_coarse, h_zone, h_fine| {
        let mesh = sos_mesh(h_coarse, h_zone, h_fine);
        common::neff(&common::solve(&mesh, &optical_spec()).unwrap().solutions[0])
    };
    let coarse = neff(0.5e-6, 0.12e-6, 40e-9);
    let fine = neff(0.5e-6, 0.06e-6, 20e-9);
    assert!(coarse.0 > 1.746 && coarse.0 < 3.45, "neff = {coarse:?}");
    assert!(coarse.1.abs() < 1e-5 && fine.1.abs() < 1e-5, "the guided mode is lossless: {coarse:?} {fine:?}");
    assert!((coarse.0 - fine.0).abs() < 1e-3, "coarse {coarse:?} against fine {fine:?}");
    assert!((fine.0 - 2.3175).abs() < 1e-3, "fine neff = {fine:?}");
}

#[test]
fn the_complex_mode_solver_with_pml_returns_the_polarization_of_the_guided_mode() {
    // the strip is x polarized (quasi-TE): the analysis comes from the complex field and leaves the absorbing layers out
    let (r, analysis) = {
        let s = common::try_solve_analysed(&sos_mesh(0.5e-6, 0.12e-6, 40e-9), &optical_spec()).unwrap();
        (s.result, s.analysis)
    };
    let a = analysis.as_ref().expect("a polarization analysis");
    let m = a.rows[0].report.as_ref().unwrap();
    assert!(m.te_fraction > 0.9 && m.dop > 0.9 && m.sphere[0] > 0.8, "{m:?}");
    assert!(m.sphere[2].abs() < 0.05, "a lossless guided mode of a mirror-symmetric strip is not circular: {m:?}");
    assert_eq!(m.class.label(), "quasi-TE (x)");
    assert!(r.scalars.iter().any(|(n, _)| n == "polarization 1"));
    assert!(r.log.iter().any(|l| l.contains("Stokes")));
}
