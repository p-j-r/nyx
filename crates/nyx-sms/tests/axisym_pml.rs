// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The cylindrical PML stretch of the axisymmetric solid kernel (K2 of the A1 plan, ADR 69): inside an absorbing layer with a radial (`x`)
//! stretch the radius becomes `r~` and `s_theta = r~ / r` enters the `1/r` terms and the weight `r`, as in `nyx-ems` (ADR 46).
//!
//! * T4 the large-radius limit against the plane-strain PML kernel, plus two exact consistency checks of the stretch (unit stretch changes
//!   nothing; a real radial stretch equals the longer ring).
//! * T5 a wheel on a post on a substrate with layers: the quality factor does not depend on the thickness or the strength of the layers.
//! * T6 the same model with `s_theta` forced to 1 (the old, untransformed form), the difference of the quality factors is printed.
//!
//! T5 and T6 use an annular wheel off the axis; T5b and T6b repeat them on the topology of benchmark A1, a solid disk on a post on a substrate that all touch the axis (K1, ADR 69).
//! Printed numbers are in the test output (`-- --nocapture`).
mod common;
use femcore::linalg::csr::C;
use femcore::props::*;
use nyx_sms::axisym::{assemble_mass_axisym, assemble_stiffness_axisym, cylindrical_stretches};
use nyx_sms::cross_section::{assemble_cross_section, solve_cross_section_modes, CrossAssembly, CrossSection};
use nyx_sms::pml::{assemble_solid_pml, solve_pml_modes};
use nyx_sms::*;

fn fixed_all(labels: Vec<i32>) -> Bc {
    Bc { labels, kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }
}

fn problem(n_modes: usize, f_search: f64) -> Problem {
    let mut s = Problem::new(Study::Static, 2);
    s.order = 2;
    s.plane = Plane::Strain;
    s.study = Study::Eigenfrequency { n_modes, f_search };
    for (k, v) in [("E", 1.0), ("nu", 0.3), ("rho", 1.0)] {
        s.materials.default.insert(k.into(), v);
    }
    s
}

fn shifted(mesh: &femcore::mesh::SimplexMesh, dx: f64) -> femcore::mesh::SimplexMesh {
    let mut m = mesh.clone();
    m.points.chunks_mut(2).for_each(|p| p[0] += dx);
    femcore::mesh::SimplexMesh::new(2, m.points, m.cells, m.cell_labels, m.facets, m.facet_labels)
}

fn modes(mesh: &femcore::mesh::SimplexMesh, s: &Problem, m: f64) -> (Vec<C>, Vec<f64>) {
    let (_, cm) = solve_solid_axisym(mesh, s, m).unwrap();
    let n1 = mesh.n_dofs(s.order);
    // the share of the squared displacement that is azimuthal (the torsional modes decouple at m = 0)
    let theta = cm
        .vectors
        .iter()
        .map(|u| {
            let tot: f64 = u.iter().map(|z| z.norm_sqr()).sum();
            let t: f64 = u[2 * n1..].iter().map(|z| z.norm_sqr()).sum();
            t / tot
        })
        .collect();
    (cm.freqs, theta)
}

/// T4a: a layer with `s = 1` (strength 0, `re` 1) is the real solver, so `s_theta = 1` exactly and the model is the unstretched ring.
#[test]
fn a_unit_stretch_changes_nothing_on_a_ring() {
    let mesh = shifted(&common::rect(0.3, 0.1, 12, 4), 1.0);
    let mut plain = problem(6, 0.0);
    plain.bcs = vec![fixed_all(vec![1])];
    let mut layered = plain.clone();
    layered
        .materials
        .pml
        .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.0, order: 0.0, re: 1.0 });
    for m in [0.0, 2.0] {
        let (a, _) = modes(&mesh, &plain, m);
        let (b, _) = modes(&mesh, &layered, m);
        for (x, y) in a.iter().zip(&b) {
            assert!(((x - y).norm() / x.norm()) < 1e-9, "m = {m}: {x} against {y}");
        }
    }
}

/// T4b: a real radial stretch `re` maps the layer `[lo, hi]` onto `[lo, lo + re (hi - lo)]`: the model is the longer ring. `s_theta` is the
/// value at the cell centre, so `r~` is exact only at the centres (error of order `h / r`, ADR 69); the error is printed for two meshes.
#[test]
fn a_real_radial_stretch_is_the_longer_ring() {
    let (r0, w, re) = (1.0, 0.3, 1.5);
    let mut err = vec![];
    for (nx, ny) in [(6, 2), (12, 4), (24, 8)] {
        let mesh = shifted(&common::rect(w, 0.1, nx, ny), r0);
        let long = shifted(&common::rect(w * re, 0.1, nx, ny), r0);
        let mut ref_p = problem(5, 0.0);
        ref_p.bcs = vec![fixed_all(vec![1])];
        let mut s = ref_p.clone();
        s.materials
            .pml
            .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.0, order: 0.0, re });
        let mut worst: f64 = 0.0;
        for m in [0.0, 3.0] {
            let (a, _) = modes(&long, &ref_p, m);
            let (b, _) = modes(&mesh, &s, m);
            for (x, y) in a.iter().zip(&b) {
                worst = worst.max((x - y).norm() / x.norm());
            }
        }
        println!("real stretch {re} on a ring of {nx} x {ny} cells: worst relative frequency error {worst:.3e}");
        err.push(worst);
    }
    // the stretch is evaluated per cell, so the error falls with h, to first order
    assert!(err[2] < 0.5 * err[0] && err[2] < 5e-3, "errors {err:?}");
}

/// T4: the large-radius limit. A clamped strip `[R0, R0 + L]` with an absorbing layer along `r`, `R0 = 1000 L`, at `m = 0`, has the
/// in-plane modes of the plane-strain PML kernel on the same strip (the cylindrical terms are of order `L / R0`); the torsional modes
/// `u_theta` are left out.
#[test]
fn at_a_large_radius_the_cylindrical_pml_is_the_plane_strain_pml() {
    let (l, r0) = (1.0, 1000.0);
    let strip = common::rect(l, 0.1, 40, 2);
    let far = shifted(&strip, r0);
    let mut s = problem(6, 0.0);
    s.bcs = vec![fixed_all(vec![1])];
    s.materials
        .pml
        .insert(1, Pml { x: PmlDir::Plus, y: PmlDir::Off, z: PmlDir::Off, strength: 0.2, order: 1.0, re: 1.0 });
    // plane strain on the strip at the origin
    let asm = assemble_solid_pml(&strip, &s).unwrap();
    let plane = solve_pml_modes(&asm).unwrap();
    let plane_f: Vec<C> = plane.freqs.iter().map(|(_, f)| *f).collect();
    // the axisymmetric strip: ask for many modes, keep the in-plane ones
    let mut sa = s.clone();
    sa.study = Study::Eigenfrequency { n_modes: 14, f_search: 0.0 };
    let (fa, theta) = modes(&far, &sa, 0.0);
    let mut inplane: Vec<C> = fa.iter().zip(&theta).filter(|(_, t)| **t < 0.5).map(|(f, _)| *f).collect();
    inplane.sort_by(|a, b| a.re.total_cmp(&b.re));
    let k = 4.min(inplane.len()).min(plane_f.len());
    assert!(k >= 3, "only {} in-plane modes found among {}", inplane.len(), fa.len());
    let mut worst: f64 = 0.0;
    for j in 0..k {
        let (a, b) = (plane_f[j], inplane[j]);
        let e = (a - b).norm() / a.norm();
        println!("large radius, mode {}: plane strain {a:.6e}, axisymmetric at R0 = {r0} L {b:.6e}, relative difference {e:.2e}", j + 1);
        worst = worst.max(e);
    }
    assert!(worst < 5e-3, "worst relative difference {worst:.3e}");
}

// ---- the wheel on a post on a substrate (T5, T6) -------------------------------------------------------------------------------------------

struct Wheel {
    mesh: femcore::mesh::SimplexMesh,
    prob: Problem,
}

/// A thin annular wheel on a short ring-shaped post on a substrate (mesh x = r, y = z), absorbing layers of thickness `t` to the left, the right
/// and below (domains 3 to 7), the outer edges clamped. Off the axis: the substrate starts at `r = 8 - 6 = 2` at the thickest layer.
fn wheel(t: f64, strength: f64, n_modes: usize, f_search: f64, dr: f64) -> Wheel {
    let rect = |x0, z0, x1, z1, d, h| common::part(x0, z0, x1, z1, d, h);
    let (a, b, d) = (8.0 + dr, 20.0 + dr, 3.0); // substrate core r in [a, b], depth d
    let parts = [
        rect(a, -d, b, 0.0, 2, None),
        rect(a - t, -d, a, 0.0, 3, None),
        rect(b, -d, b + t, 0.0, 4, None),
        rect(a, -d - t, b, -d, 5, None),
        rect(a - t, -d - t, a, -d, 6, None),
        rect(b, -d - t, b + t, -d, 7, None),
        rect(13.5 + dr, 0.0, 14.5 + dr, 1.0, 1, Some(0.25)),
        rect(12.0 + dr, 1.0, 16.0 + dr, 1.4, 1, Some(0.2)),
    ];
    let mesh = common::parts_mesh(&parts, 0.8);
    let (xl, xr, zb) = (a - t, b + t, -d - t);
    let on = |v: f64, w: f64| (v - w).abs() < 1e-9;
    let mut outer: Vec<i32> = (0..mesh.n_facets())
        .filter(|&f| {
            let (p, q) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
            (on(p[0], xl) && on(q[0], xl)) || (on(p[0], xr) && on(q[0], xr)) || (on(p[1], zb) && on(q[1], zb))
        })
        .map(|f| mesh.facet_labels[f])
        .collect();
    outer.sort();
    outer.dedup();
    let mut prob = problem(n_modes, f_search);
    prob.bcs = vec![fixed_all(outer)];
    let st = |x, z| Pml { x, y: PmlDir::Off, z, strength, order: 2.0, re: 1.0 };
    prob.materials.pml.insert(3, st(PmlDir::Minus, PmlDir::Off));
    prob.materials.pml.insert(4, st(PmlDir::Plus, PmlDir::Off));
    prob.materials.pml.insert(5, st(PmlDir::Off, PmlDir::Minus));
    prob.materials.pml.insert(6, st(PmlDir::Minus, PmlDir::Minus));
    prob.materials.pml.insert(7, st(PmlDir::Plus, PmlDir::Minus));
    Wheel { mesh, prob }
}

/// Complex frequency, quality factor and the share of the mass (`rho |u|^2 r`) that sits in the wheel and post (domain 1), per mode.
struct Mode {
    f: C,
    q: f64,
    share: f64,
}

fn wheel_modes(w: &Wheel, m: f64, theta_one: bool) -> Vec<Mode> {
    let mut asm = assemble_cross_section(&w.mesh, &w.prob, CrossSection::Axisymmetric(m)).unwrap();
    if theta_one {
        old_form(&w.mesh, &w.prob, &mut asm, m);
    }
    let solved = solve_cross_section_modes(&asm).unwrap();
    let order = w.prob.order;
    let physical = |rho: &[f64]| assemble_mass_axisym(&w.mesh, rho, order, None);
    let rho_res: Vec<f64> =
        (0..w.mesh.n_cells()).map(|c| if w.mesh.cell_labels[c] == 1 { asm.rho[c] } else { 0.0 }).collect();
    let (m_all, m_res) = (physical(&asm.rho), physical(&rho_res));
    let quad = |mm: &femcore::linalg::csr::CCsr, u: &[C]| -> f64 {
        mm.matvec(u).iter().zip(u).map(|(a, b)| (b.conj() * a).re).sum()
    };
    solved
        .modes
        .freqs
        .iter()
        .zip(&solved.modes.vectors)
        .map(|(f, u)| Mode { f: *f, q: f.re / (2.0 * f.im.abs()), share: quad(&m_res, u) / quad(&m_all, u) })
        .collect()
}

/// The old form: the same assembly with `s_theta = 1` (`1/r` and the weight `r` untransformed), through the public stage pieces.
fn old_form(mesh: &femcore::mesh::SimplexMesh, p: &Problem, asm: &mut CrossAssembly, m: f64) {
    let mut st = cylindrical_stretches(mesh, &p.materials.pml);
    st.iter_mut().for_each(|s| s[2] = C::new(1.0, 0.0));
    let k = assemble_stiffness_axisym(mesh, &asm.c6, p.order, m, Some(&st)).unwrap();
    let mm = assemble_mass_axisym(mesh, &asm.rho, p.order, Some(&st));
    asm.kf = k.submatrix(&asm.free, &asm.free);
    asm.mf = mm.submatrix(&asm.free, &asm.free);
}

/// The mode of the wheel nearest `f0` among those of `modes` that sit mostly in the wheel.
fn pick(modes: &[Mode], f0: f64) -> &Mode {
    modes
        .iter()
        .filter(|x| x.share > 0.3)
        .min_by(|a, b| (a.f.re - f0).abs().total_cmp(&(b.f.re - f0).abs()))
        .expect("no mode in the wheel")
}

const F_WHEEL: f64 = 0.086;

/// T5: the quality factor of the wheel mode near `F_WHEEL` does not depend on the layers: it changes by less than 5 % when the thickness
/// of the layers doubles (3 to 6) and when their strength doubles (3 to 6).
#[test]
fn the_quality_factor_does_not_depend_on_the_layers() {
    let run = |t: f64, strength: f64| {
        let w = wheel(t, strength, 14, F_WHEEL, 0.0);
        let all = wheel_modes(&w, 0.0, false);
        let x = pick(&all, F_WHEEL);
        println!(
            "layers t = {t}, strength {strength} ({} cells): f = {:.6e} {:+.3e}i, Q = {:.4}, share in the wheel {:.3}",
            w.mesh.n_cells(),
            x.f.re,
            x.f.im,
            x.q,
            x.share
        );
        (x.f.re, x.q)
    };
    let (f0, q0) = run(3.0, 3.0);
    let (f1, q1) = run(6.0, 3.0);
    let (f2, q2) = run(3.0, 6.0);
    for (f, q, what) in [(f1, q1, "thickness 3 -> 6"), (f2, q2, "strength 3 -> 6")] {
        let (df, dq) = ((f / f0 - 1.0).abs(), (q / q0 - 1.0).abs());
        println!("{what}: Re f changes by {df:.2e}, Q by {dq:.2e}");
        assert!(df < 5e-3 && dq < 5e-2, "{what}: Re f {f0} -> {f}, Q {q0} -> {q}");
    }
}

/// T6: the same model with the old form (`s_theta` forced to 1, `r` untransformed in the layers, the choice of Bindel and Govindjee). The
/// difference of the quality factors is printed, for the wheel at `r = 14` and moved towards the axis to `r = 11` where the layers see more
/// curvature; only finiteness is asserted (the numbers go into ADR 69).
#[test]
fn the_untransformed_form_gives_a_finite_quality_factor_and_its_difference_is_printed() {
    for (t, strength, dr) in [(3.0, 3.0, 0.0), (6.0, 3.0, 0.0), (3.0, 3.0, -3.0), (3.0, 6.0, -3.0)] {
        let w = wheel(t, strength, 14, F_WHEEL, dr);
        let new = wheel_modes(&w, 0.0, false);
        let old = wheel_modes(&w, 0.0, true);
        let (new, old) = (pick(&new, F_WHEEL), pick(&old, F_WHEEL));
        println!(
            "wheel moved by {dr}, layers t = {t}, strength {strength}: f = {:.5e}, Q with s_theta = r~/r {:.4}, with s_theta = 1 {:.4}, relative difference {:.3e}",
            new.f.re,
            new.q,
            old.q,
            (old.q / new.q - 1.0).abs()
        );
        assert!(old.q.is_finite() && old.q > 0.0 && new.q.is_finite() && new.q > 0.0);
    }
}

// ---- the same on a solid disk on a post on a substrate, all touching the axis (T5b, T6b) ----------------------------------------------------

/// A solid disk (`r` 0 to 8, 0.4 thick) on a post (`r` 0 to 1, 1 tall) on a substrate (`r` 0 to 20, depth 3) with absorbing layers of thickness `t`
/// to the right and below (domains 4, 5, 7), the outer edges clamped. Domain 1 is the disk and post; every part touches the axis.
fn disk_on_post(t: f64, strength: f64, n_modes: usize, f_search: f64) -> Wheel {
    let rect = |x0, z0, x1, z1, d, h| common::part(x0, z0, x1, z1, d, h);
    let (b, d) = (20.0, 3.0);
    let parts = [
        rect(0.0, -d, b, 0.0, 2, None),
        rect(b, -d, b + t, 0.0, 4, None),
        rect(0.0, -d - t, b, -d, 5, None),
        rect(b, -d - t, b + t, -d, 7, None),
        rect(0.0, 0.0, 1.0, 1.0, 1, Some(0.25)),
        rect(0.0, 1.0, 8.0, 1.4, 1, Some(0.2)),
    ];
    let mesh = common::parts_mesh(&parts, 0.8);
    let (xr, zb) = (b + t, -d - t);
    let on = |v: f64, w: f64| (v - w).abs() < 1e-9;
    let mut outer: Vec<i32> = (0..mesh.n_facets())
        .filter(|&f| {
            let (p, q) = (mesh.point(mesh.facet(f)[0] as usize), mesh.point(mesh.facet(f)[1] as usize));
            (on(p[0], xr) && on(q[0], xr)) || (on(p[1], zb) && on(q[1], zb))
        })
        .map(|f| mesh.facet_labels[f])
        .collect();
    outer.sort();
    outer.dedup();
    let mut prob = problem(n_modes, f_search);
    prob.bcs = vec![fixed_all(outer)];
    let st = |x, z| Pml { x, y: PmlDir::Off, z, strength, order: 2.0, re: 1.0 };
    prob.materials.pml.insert(4, st(PmlDir::Plus, PmlDir::Off));
    prob.materials.pml.insert(5, st(PmlDir::Off, PmlDir::Minus));
    prob.materials.pml.insert(7, st(PmlDir::Plus, PmlDir::Minus));
    Wheel { mesh, prob }
}

/// The first radial contour mode of the free disk, `lambda = 2.05` for `nu = 0.3` (plane stress): `f = lambda c / (2 pi R)`, `c = sqrt(E / (rho (1 - nu^2)))`.
const F_DISK: f64 = 2.05 * 1.0482 / (2.0 * std::f64::consts::PI * 8.0);

/// T5b: on a disk that touches the axis, the quality factor of the contour mode does not depend on the layers (5 % for a doubled thickness or
/// strength), and the frequency stays within 10 % of the free disk's closed form (the post and the substrate load it).
#[test]
fn the_quality_factor_of_a_disk_on_the_axis_does_not_depend_on_the_layers() {
    let run = |t: f64, strength: f64| {
        let w = disk_on_post(t, strength, 14, F_DISK);
        let all = wheel_modes(&w, 0.0, false);
        let x = pick(&all, F_DISK);
        println!(
            "disk on the axis, layers t = {t}, strength {strength} ({} cells): f = {:.6e} {:+.3e}i, Q = {:.4}, share in the disk {:.3}",
            w.mesh.n_cells(),
            x.f.re,
            x.f.im,
            x.q,
            x.share
        );
        (x.f.re, x.q)
    };
    let (f0, q0) = run(3.0, 3.0);
    assert!((f0 / F_DISK - 1.0).abs() < 0.1, "{f0} against the free disk {F_DISK}");
    assert!(q0.is_finite() && q0 > 0.0);
    let (f1, q1) = run(6.0, 3.0);
    let (f2, q2) = run(3.0, 6.0);
    for (f, q, what) in [(f1, q1, "thickness 3 -> 6"), (f2, q2, "strength 3 -> 6")] {
        let (df, dq) = ((f / f0 - 1.0).abs(), (q / q0 - 1.0).abs());
        println!("{what}: Re f changes by {df:.2e}, Q by {dq:.2e}");
        assert!(df < 5e-3 && dq < 5e-2, "{what}: Re f {f0} -> {f}, Q {q0} -> {q}");
    }
}

/// T6b: the old (untransformed) form on the disk that touches the axis: the difference of the quality factors is printed, only finiteness asserted.
#[test]
fn the_untransformed_form_on_a_disk_on_the_axis_gives_a_finite_quality_factor() {
    for (t, strength) in [(3.0, 3.0), (6.0, 3.0), (3.0, 6.0)] {
        let w = disk_on_post(t, strength, 14, F_DISK);
        let new = wheel_modes(&w, 0.0, false);
        let old = wheel_modes(&w, 0.0, true);
        let (new, old) = (pick(&new, F_DISK), pick(&old, F_DISK));
        println!(
            "disk on the axis, layers t = {t}, strength {strength}: f = {:.5e}, Q with s_theta = r~/r {:.4}, with s_theta = 1 {:.4}, relative difference {:.3e}",
            new.f.re,
            new.q,
            old.q,
            (old.q / new.q - 1.0).abs()
        );
        assert!(old.q.is_finite() && old.q > 0.0 && new.q.is_finite() && new.q > 0.0);
    }
}
