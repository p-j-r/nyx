//! The axis conditions of the axisymmetric solid kernel (K1, ADR 69) against closed forms of bodies that touch the axis `r = 0`:
//! torsion and radial modes of a solid cylinder (exact for any length) and the radial contour modes of a thin solid disk (the mode
//! family of benchmark A1). Printed relative errors are in the test output (`-- --nocapture`).
mod common;
use femcore::linalg::csr::C;
use nyx_physics::Discretisation;
use nyx_sms::*;

const PI: f64 = std::f64::consts::PI;

/// Bessel function of the first kind by the periodic trapezoid rule on `J_n(x) = (1 / 2 pi) int cos(n t - x sin t) dt` (spectrally accurate).
fn bessel(n: i32, x: f64) -> f64 {
    let steps = 400;
    (0..steps)
        .map(|k| {
            let t = 2.0 * PI * k as f64 / steps as f64;
            (n as f64 * t - x * t.sin()).cos()
        })
        .sum::<f64>()
        / steps as f64
}

/// Roots of `g` in `(lo, hi)` found by scanning for sign changes and bisecting; poles of `g` are rejected by a residual check.
fn roots(g: impl Fn(f64) -> f64, lo: f64, hi: f64) -> Vec<f64> {
    let mut out = vec![];
    let step = 0.005;
    let mut x = lo;
    while x + step < hi {
        let (a, b) = (x, x + step);
        if g(a) * g(b) < 0.0 {
            let (mut a, mut b) = (a, b);
            for _ in 0..60 {
                let mid = 0.5 * (a + b);
                if g(a) * g(mid) <= 0.0 {
                    b = mid;
                } else {
                    a = mid;
                }
            }
            let r = 0.5 * (a + b);
            if g(r).abs() < 1e-5 {
                out.push(r);
            }
        }
        x += step;
    }
    out
}

fn problem(e: f64, nu: f64, rho: f64, n_modes: usize, f_search: f64) -> Problem {
    let mut sp = Problem::new(Study::Static, 2);
    sp.order = 2;
    sp.materials.set(1, "E", e);
    sp.materials.set(1, "nu", nu);
    sp.materials.set(1, "rho", rho);
    sp.study = Study::Eigenfrequency { n_modes, f_search };
    sp
}

/// A solid cylinder of radius `rr` and length `len` (mesh x = r from the axis, y = z from `-len / 2` to `len / 2`).
fn cylinder(rr: f64, len: f64, h: f64) -> Discretisation {
    common::parts_mesh(&[common::part(0.0, -len / 2.0, rr, len / 2.0, 1, Some(h))], h)
}

/// The boundary numbers of the two end faces (`|z| = len / 2`).
fn end_labels(mesh: &Discretisation, len: f64) -> Vec<i32> {
    let mut labels = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let fv = mesh.facet(f);
        if fv.iter().all(|&v| (mesh.point(v as usize)[1].abs() - len / 2.0).abs() < 1e-9 * len) {
            labels.insert(mesh.facet_labels[f]);
        }
    }
    labels.into_iter().collect()
}

/// Share of the norm `sum |u_k|^2` that sits in displacement block `k` (0 u_r, 1 u_z, 2 u_theta).
fn block_share(v: &[C], n1: usize, k: usize) -> f64 {
    let part = |b: usize| -> f64 { v[b * n1..(b + 1) * n1].iter().map(|z| z.norm_sqr()).sum() };
    part(k) / (part(0) + part(1) + part(2))
}

/// The mode of `modes` whose frequency is nearest `target`: its index and frequency.
fn nearest_mode(modes: &CrossModes, target: f64) -> (usize, f64) {
    modes
        .freqs
        .iter()
        .enumerate()
        .min_by(|a, b| (a.1.re - target).abs().total_cmp(&(b.1.re - target).abs()))
        .map(|(i, z)| (i, z.re))
        .unwrap()
}

/// T1: torsion of a free solid cylinder, `m = 0`: `u_theta = J1(k r)` with `J2(k R) = 0`, `f = x c_s / (2 pi R)`; exact for any length.
/// It puts `u_theta` on the axis.
#[test]
fn a_solid_cylinder_has_the_torsional_frequencies_of_the_j1_solution() {
    let (e, nu, rho, rr, len) = (169e9_f64, 0.28_f64, 2329.0_f64, 1.0_f64, 1.0_f64);
    let cs = (e / (2.0 * (1.0 + nu)) / rho).sqrt();
    let zeros = roots(|x| bessel(2, x), 3.0, 12.0);
    assert!(zeros.len() >= 2 && (zeros[0] - 5.13562).abs() < 1e-4 && (zeros[1] - 8.41724).abs() < 1e-4, "{zeros:?}");
    let mesh = cylinder(rr, len, 0.1);
    for (i, &x) in zeros.iter().take(2).enumerate() {
        let want = x * cs / (2.0 * PI * rr);
        let sp = problem(e, nu, rho, 6, want);
        let (_, modes) = solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
        let (k, got) = nearest_mode(&modes, want);
        let share = block_share(&modes.vectors[k], mesh.n_dofs(2), 2);
        let rel = (got - want).abs() / want;
        println!(
            "T1 torsion mode {}: {got:.4} Hz against {want:.4} Hz, relative error {rel:.2e}, u_theta share {share:.3}",
            i + 1
        );
        assert!(share > 0.99, "the nearest mode is not torsional (share {share})");
        assert!(rel < 5e-4, "{got} against {want}");
    }
}

/// T2: radial plane-strain modes of a solid cylinder with `u_z = 0` on both ends, `m = 0`: `u_r = J1(k r)`,
/// `k R J0(k R) / J1(k R) = (1 - 2 nu) / (1 - nu)`, `f = k c_p / (2 pi)` with `c_p^2 = (lambda + 2 mu) / rho`. It puts `u_r` on the axis.
#[test]
fn a_solid_cylinder_has_the_radial_frequencies_of_the_j1_solution() {
    let (e, nu, rho, rr, len) = (169e9_f64, 0.28_f64, 2329.0_f64, 1.0_f64, 0.5_f64);
    let cp = (e * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) / rho).sqrt();
    let rhs = (1.0 - 2.0 * nu) / (1.0 - nu);
    let kr = roots(|x| x * bessel(0, x) / bessel(1, x) - rhs, 0.5, 12.0);
    assert!(kr.len() >= 2, "{kr:?}");
    let mesh = cylinder(rr, len, 0.1);
    let mut sp = problem(e, nu, rho, 6, 0.0);
    let labels = end_labels(&mesh, len);
    assert!(!labels.is_empty());
    sp.bcs.push(Bc { labels, kind: BcKind::Fixed { fixed: [false, false, true], values: [0.0; 3] } });
    for (i, &x) in kr.iter().take(2).enumerate() {
        let want = x * cp / (2.0 * PI * rr);
        sp.study = Study::Eigenfrequency { n_modes: 6, f_search: want };
        let (_, modes) = solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
        let (k, got) = nearest_mode(&modes, want);
        let share = block_share(&modes.vectors[k], mesh.n_dofs(2), 0);
        let rel = (got - want).abs() / want;
        println!(
            "T2 radial mode {}: {got:.4} Hz against {want:.4} Hz, relative error {rel:.2e}, u_r share {share:.3}",
            i + 1
        );
        assert!(share > 0.99, "the nearest mode is not radial (share {share})");
        assert!(rel < 1e-4, "{got} against {want}");
    }
}

/// T3: thin solid disk (`t / R = 0.02`), `m = 0`, radial contour modes: `lambda J0(lambda) / J1(lambda) = 1 - nu`,
/// `f = lambda / (2 pi R) sqrt(E / (rho (1 - nu^2)))`, first two roots; the thin-plate error is of order `(t / R)^2`. This is the mode family of A1.
#[test]
fn a_thin_solid_disk_has_the_radial_contour_frequencies() {
    let (e, nu, rho, rr, t) = (169e9_f64, 0.28_f64, 2329.0_f64, 1.0_f64, 0.02_f64);
    let c = (e / (rho * (1.0 - nu * nu))).sqrt();
    let lam = roots(|x| x * bessel(0, x) / bessel(1, x) - (1.0 - nu), 0.5, 10.0);
    assert!(lam.len() >= 2, "{lam:?}");
    let mesh = cylinder(rr, t, t / 2.0);
    for (i, &x) in lam.iter().take(2).enumerate() {
        let want = x * c / (2.0 * PI * rr);
        let sp = problem(e, nu, rho, 6, want);
        let (_, modes) = solve_solid_axisym(&mesh, &sp, 0.0).unwrap();
        let (k, got) = nearest_mode(&modes, want);
        let share = block_share(&modes.vectors[k], mesh.n_dofs(2), 0);
        let rel = (got - want).abs() / want;
        println!(
            "T3 contour mode {}: {got:.4} Hz against {want:.4} Hz, relative error {rel:.2e}, u_r share {share:.3}",
            i + 1
        );
        assert!(share > 0.9, "the nearest mode is not radial (share {share})");
        assert!(rel < 5e-4, "{got} against {want}");
    }
}

/// T3, `m = 2` on a disk (no closed form claimed): the lowest frequency is positive, finite and stable under one refinement to 1 %.
#[test]
fn the_lowest_m2_mode_of_a_disk_is_finite_and_stable_under_refinement() {
    let (e, nu, rho, rr, t) = (169e9_f64, 0.28_f64, 2329.0_f64, 1.0_f64, 0.1_f64);
    let lowest = |h: f64| -> f64 {
        let mesh = cylinder(rr, t, h);
        let sp = problem(e, nu, rho, 3, 0.0);
        let (_, modes) = solve_solid_axisym(&mesh, &sp, 2.0).unwrap();
        modes.freqs.iter().map(|z| z.re).fold(f64::INFINITY, f64::min)
    };
    let (coarse, fine) = (lowest(0.05), lowest(0.025));
    println!(
        "T3 m = 2 lowest frequency: {coarse:.4} Hz (h = 0.05), {fine:.4} Hz (h = 0.025), change {:.2e}",
        (coarse - fine).abs() / fine
    );
    assert!(coarse.is_finite() && coarse > 0.0 && fine.is_finite() && fine > 0.0);
    assert!((coarse - fine).abs() < 0.01 * fine, "{coarse} against {fine}");
}

/// `m = 1` on a free solid cylinder: the rigid translation across the axis and the rigid tilt have no strain, so two modes sit at zero
/// frequency (the tie `u_r = i m u_t` holds for them exactly). Checks that the `m = +-1` rule holds `u_z` on the axis and nothing else.
#[test]
fn m1_of_a_free_cylinder_keeps_its_two_rigid_modes() {
    let (e, nu, rho) = (169e9_f64, 0.28_f64, 2329.0_f64);
    let mesh = cylinder(1.0, 1.0, 0.1);
    let sp = problem(e, nu, rho, 4, 0.0);
    let (_, modes) = solve_solid_axisym(&mesh, &sp, 1.0).unwrap();
    let mut f: Vec<f64> = modes.freqs.iter().map(|z| z.re.abs()).collect();
    f.sort_by(|a, b| a.total_cmp(b));
    println!("m = 1 free cylinder, lowest |f|: {f:?}");
    assert!(f[1] < 1e-3 * f[2], "{f:?}");
}

/// A cell that crosses the axis is refused by name.
#[test]
fn a_mesh_that_crosses_the_axis_is_refused() {
    let mesh = common::parts_mesh(&[common::part(-0.1, 0.0, 0.5, 0.2, 1, Some(0.1))], 0.1);
    let sp = problem(169e9, 0.28, 2329.0, 3, 0.0);
    let err = solve_solid_axisym(&mesh, &sp, 0.0).err().expect("must be refused");
    assert!(err.contains("r >= 0"), "{err}");
}
