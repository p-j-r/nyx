//! The 3D mesher (isosurface stuffing on a body-centred cubic lattice): exact boxes, converging
//! curved volumes, watertight conforming meshes with bounded dihedral angles, labels, and
//! solutions on the meshes.
mod common;
use femcore::mesh::SimplexMesh;
use femcore::mesher3d::*;
use femcore::problem::*;
use std::collections::HashMap;

fn spec(bodies: Vec<(Solid, i32)>, h: f64) -> Spec3d {
    Spec3d { bodies: bodies.into_iter().map(|(solid, domain)| Body3d { solid, domain, h: None }).collect(), h, alpha: 0.3, snap_features: true, min_layers: 2.0 }
}

fn pt(m: &SimplexMesh, v: u32) -> [f64; 3] {
    [m.points[3 * v as usize], m.points[3 * v as usize + 1], m.points[3 * v as usize + 2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn volumes(m: &SimplexMesh) -> Vec<f64> {
    (0..m.n_cells())
        .map(|c| {
            let v = m.cell(c);
            let p0 = pt(m, v[0]);
            dot(sub(pt(m, v[1]), p0), cross(sub(pt(m, v[2]), p0), sub(pt(m, v[3]), p0))) / 6.0
        })
        .collect()
}

fn domain_volume(m: &SimplexMesh, d: i32) -> f64 {
    volumes(m).iter().zip(&m.cell_labels).filter(|(_, &l)| l == d).map(|(v, _)| v).sum()
}

/// Every interior face is shared by two cells, the faces used once are exactly the boundary
/// facets, and every cell has positive volume.
fn check_conforming(m: &SimplexMesh) {
    assert!(volumes(m).iter().all(|&v| v > 0.0), "all cells positively oriented and non-degenerate");
    let mut faces: HashMap<[u32; 3], usize> = HashMap::new();
    for c in 0..m.n_cells() {
        let v = m.cell(c);
        for skip in 0..4 {
            let mut f: Vec<u32> = (0..4).filter(|&i| i != skip).map(|i| v[i]).collect();
            f.sort_unstable();
            *faces.entry([f[0], f[1], f[2]]).or_insert(0) += 1;
        }
    }
    assert!(faces.values().all(|&n| n <= 2), "no face in more than two cells");
    let once: std::collections::HashSet<[u32; 3]> = faces.iter().filter(|(_, &n)| n == 1).map(|(k, _)| *k).collect();
    let given: std::collections::HashSet<[u32; 3]> = (0..m.n_facets())
        .map(|f| {
            let mut t = [m.facets[3 * f], m.facets[3 * f + 1], m.facets[3 * f + 2]];
            t.sort_unstable();
            t
        })
        .collect();
    assert_eq!(given.len(), m.n_facets(), "no duplicate facets");
    assert_eq!(once, given, "boundary facets are exactly the faces used by one cell");
}

fn min_dihedral_deg(m: &SimplexMesh) -> f64 {
    let mut best = 180.0f64;
    for c in 0..m.n_cells() {
        let v: Vec<[f64; 3]> = m.cell(c).iter().map(|&i| pt(m, i)).collect();
        for (a, b) in [(0usize, 1usize), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let others: Vec<usize> = (0..4).filter(|&i| i != a && i != b).collect();
            let e = sub(v[b], v[a]);
            let en = dot(e, e).sqrt();
            let proj = |p: [f64; 3]| {
                let d = sub(p, v[a]);
                let t = dot(d, e) / (en * en);
                [d[0] - t * e[0], d[1] - t * e[1], d[2] - t * e[2]]
            };
            let (u, w) = (proj(v[others[0]]), proj(v[others[1]]));
            let cosang = dot(u, w) / (dot(u, u).sqrt() * dot(w, w).sqrt());
            best = best.min(cosang.clamp(-1.0, 1.0).acos().to_degrees());
        }
    }
    best
}

fn facet_area(m: &SimplexMesh, f: usize) -> f64 {
    let (a, b, c) = (pt(m, m.facets[3 * f]), pt(m, m.facets[3 * f + 1]), pt(m, m.facets[3 * f + 2]));
    let x = cross(sub(b, a), sub(c, a));
    0.5 * dot(x, x).sqrt()
}

#[test]
fn a_box_is_reproduced_exactly() {
    let s = spec(vec![(Solid::Box { min: [0.0, 0.0, 0.0], max: [1.0, 0.5, 0.25] }, 1)], 0.1);
    let (m, rep) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert!((domain_volume(&m, 1) - 0.125).abs() < 1e-12, "{}", domain_volume(&m, 1));
    assert_eq!(m.boundary_labels(), vec![1, 2, 3, 4, 5, 6]);
    let mut area = [0.0; 7];
    for f in 0..m.n_facets() {
        area[m.facet_labels[f] as usize] += facet_area(&m, f);
    }
    for (l, want) in [(1, 0.125), (2, 0.125), (3, 0.25), (4, 0.25), (5, 0.5), (6, 0.5)] {
        assert!((area[l] - want).abs() < 1e-12, "face {l}: {} vs {want}", area[l]);
    }
    assert_eq!(rep.warped, 0);
}

#[test]
fn a_sphere_converges_and_is_watertight() {
    let exact = 4.0 / 3.0 * std::f64::consts::PI;
    let mut errs = vec![];
    for h in [0.2, 0.1, 0.05] {
        let s = spec(vec![(Solid::Sphere { c: [0.0; 3], r: 1.0 }, 1)], h);
        let (m, _) = mesh_3d(&s).unwrap();
        check_conforming(&m);
        let v = domain_volume(&m, 1);
        errs.push((v - exact).abs() / exact);
        // a closed surface: every edge of the boundary facets is in exactly two of them
        let mut edges: HashMap<(u32, u32), usize> = HashMap::new();
        for f in 0..m.n_facets() {
            let t = [m.facets[3 * f], m.facets[3 * f + 1], m.facets[3 * f + 2]];
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        assert!(edges.values().all(|&n| n == 2), "h = {h}: the boundary is not a closed surface");
        // surface points lie on the sphere
        for f in 0..m.n_facets() {
            for k in 0..3 {
                let p = pt(&m, m.facets[3 * f + k]);
                assert!((dot(p, p).sqrt() - 1.0).abs() < 2e-3, "h = {h}: surface vertex off the sphere");
            }
        }
        let ang = min_dihedral_deg(&m);
        assert!(ang > 8.0, "h = {h}: smallest dihedral angle {ang}");
    }
    assert!(errs[2] < errs[1] && errs[1] < errs[0] && errs[2] < 5e-3, "volume errors {errs:?}");
}

#[test]
fn a_hole_and_boundary_labels() {
    // a plate with a cylindrical hole through it: the hole wall is face 1 of the second body
    let s = spec(vec![(Solid::Box { min: [0.0, 0.0, 0.0], max: [2.0, 1.0, 0.4] }, 1), (Solid::Cylinder { base: [1.0, 0.5, -1.0], axis: 2, r: 0.25, len: 3.0 }, 0)], 0.05);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    let want = 2.0 * 1.0 * 0.4 - std::f64::consts::PI * 0.25 * 0.25 * 0.4;
    assert!((domain_volume(&m, 1) - want).abs() / want < 5e-3, "{} vs {want}", domain_volume(&m, 1));
    let labels = m.boundary_labels();
    assert_eq!(labels, vec![1, 2, 3, 4, 5, 6, 7], "box faces 1..6 and the hole wall 7: {labels:?}");
    for f in 0..m.n_facets() {
        if m.facet_labels[f] == 7 {
            for k in 0..3 {
                let p = pt(&m, m.facets[3 * f + k]);
                // the rim at the plate faces is a sharp edge: vertices there are within a fraction of h
                let off = (((p[0] - 1.0).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - 0.25).abs();
                assert!(off < 0.5 * 0.05, "hole wall vertex {off} off the cylinder");
            }
        }
    }
}

#[test]
fn two_bodies_get_their_domain_numbers() {
    let s = spec(vec![(Solid::Box { min: [0.0, 0.0, 0.0], max: [1.0, 1.0, 0.5] }, 1), (Solid::Box { min: [0.5, 0.0, 0.0], max: [1.0, 1.0, 0.5] }, 2)], 0.1);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert_eq!(m.domain_labels(), vec![1, 2]);
    let (v1, v2) = (domain_volume(&m, 1), domain_volume(&m, 2));
    assert!((v1 + v2 - 0.5).abs() < 1e-12);
    assert!((v2 - 0.25).abs() < 0.05 * 0.25, "interface within an element size: {v2}");
}

#[test]
fn uniaxial_tension_of_a_block_is_exact() {
    // with nu = 0 the displacement uz = sigma z / E is linear, so any conforming mesh of the bar
    // must reproduce it to rounding: this checks cells, orientation and the boundary labels together
    let s = spec(vec![(Solid::Box { min: [0.0, 0.0, 0.0], max: [0.3, 0.2, 1.0] }, 1)], 0.05);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.materials.default.insert("E".into(), 1e9);
    sp.materials.default.insert("nu".into(), 0.0);
    sp.bcs = vec![Bc { on: None, labels: vec![5], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![6], kind: BcKind::Load([0.0, 0.0, 1e6]) }];
    let r = femcore::solid::solve_solid(&m, &sp, &SolveState::default()).unwrap();
    let uz = &r.solutions[0].field("uz").unwrap().values;
    let x = m.dof_coordinates(2);
    let exact: Vec<f64> = (0..uz.len()).map(|i| 1e6 * x[3 * i + 2] / 1e9).collect();
    assert!(common::maxdiff(uz, &exact) < 1e-12, "{}", common::maxdiff(uz, &exact));
}

#[test]
fn tension_of_a_bar_with_a_locally_refined_slab_is_exact() {
    // the same bar with a slab (domain 2, off the lattice) that has a finer mesh: still exact, so the
    // refined cells, their orientation, the split facets and the labels are consistent
    let mut s = spec(vec![(Solid::Box { min: [0.0, 0.0, 0.0], max: [0.3, 0.2, 1.0] }, 1), (Solid::Box { min: [0.0, 0.0, 0.33], max: [0.3, 0.2, 0.57] }, 2)], 0.1);
    s.bodies[1].h = Some(0.03);
    let (m, rep) = mesh_3d(&s).unwrap();
    assert!(rep.refined > 0);
    check_conforming(&m);
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.materials.default.insert("E".into(), 1e9);
    sp.materials.default.insert("nu".into(), 0.0);
    sp.bcs = vec![Bc { on: None, labels: vec![5], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![6], kind: BcKind::Load([0.0, 0.0, 1e6]) }];
    let r = femcore::solid::solve_solid(&m, &sp, &SolveState::default()).unwrap();
    let uz = &r.solutions[0].field("uz").unwrap().values;
    let x = m.dof_coordinates(2);
    let exact: Vec<f64> = (0..uz.len()).map(|i| 1e6 * x[3 * i + 2] / 1e9).collect();
    assert!(common::maxdiff(uz, &exact) < 1e-11, "{}", common::maxdiff(uz, &exact));
}

#[test]
fn tension_of_a_cylinder_is_close_to_the_bar_solution() {
    // the polygonal approximation of the circle changes the cross-section by a few percent
    let s = spec(vec![(Solid::Cylinder { base: [0.0, 0.0, 0.0], axis: 2, r: 0.1, len: 1.0 }, 1)], 0.03);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert_eq!(m.boundary_labels(), vec![1, 2, 3]);
    let mut sp = ProblemSpec::new(Physics::Solid);
    sp.materials.default.insert("E".into(), 1e9);
    sp.materials.default.insert("nu".into(), 0.0);
    sp.bcs = vec![Bc { on: None, labels: vec![2], kind: BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] } }, Bc { on: None, labels: vec![3], kind: BcKind::Load([0.0, 0.0, 1e6]) }];
    let r = femcore::solid::solve_solid(&m, &sp, &SolveState::default()).unwrap();
    let uz = &r.solutions[0].field("uz").unwrap().values;
    let x = m.dof_coordinates(2);
    let exact: Vec<f64> = (0..uz.len()).map(|i| 1e6 * x[3 * i + 2] / 1e9).collect();
    let scale = exact.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    assert!(common::maxdiff(uz, &exact) < 0.05 * scale, "{} of {scale}", common::maxdiff(uz, &exact));
}

#[test]
fn rejects_bad_input() {
    assert!(mesh_3d(&Spec3d::default()).is_err());
    let hole_only = spec(vec![(Solid::Sphere { c: [0.0; 3], r: 1.0 }, 0)], 0.1);
    assert!(mesh_3d(&hole_only).is_err());
    let bad_h = spec(vec![(Solid::Sphere { c: [0.0; 3], r: 1.0 }, 1)], -1.0);
    assert!(mesh_3d(&bad_h).is_err());
    let too_fine = spec(vec![(Solid::Sphere { c: [0.0; 3], r: 1.0 }, 1)], 0.001);
    assert!(mesh_3d(&too_fine).is_err());
}

fn box_body(min: [f64; 3], max: [f64; 3], d: i32) -> (Solid, i32) {
    (Solid::Box { min, max }, d)
}

#[test]
fn a_planar_material_interface_off_the_lattice_is_fitted_exactly() {
    // the interface x = 0.37 is not on the lattice (h = 0.1): the cells are cut along it, so the volumes are exact
    let s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 0.5], 1), box_body([0.37, 0.0, 0.0], [1.0, 1.0, 0.5], 2)], 0.1);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert!((domain_volume(&m, 2) - 0.63 * 0.5).abs() < 1e-12, "{}", domain_volume(&m, 2));
    assert!((domain_volume(&m, 1) - 0.37 * 0.5).abs() < 1e-12, "{}", domain_volume(&m, 1));
    // no cell has vertices strictly on both sides of the interface
    for c in 0..m.n_cells() {
        let xs: Vec<f64> = m.cell(c).iter().map(|&v| pt(&m, v)[0]).collect();
        let (lo, hi) = (xs.iter().cloned().fold(f64::INFINITY, f64::min), xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
        assert!(hi <= 0.37 + 1e-9 || lo >= 0.37 - 1e-9, "cell {c} straddles the interface: {lo}..{hi}");
    }
}

#[test]
fn a_curved_material_interface_converges_faster_than_the_element_size() {
    // a ball of domain 2 inside a block of domain 1: no cell straddles the sphere, and its volume converges
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 0.25f64.powi(3);
    let mut errs = vec![];
    for h in [0.1, 0.05] {
        let s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1), (Solid::Sphere { c: [0.5, 0.5, 0.5], r: 0.25 }, 2)], h);
        let (m, _) = mesh_3d(&s).unwrap();
        check_conforming(&m);
        for c in 0..m.n_cells() {
            let d: Vec<f64> = m.cell(c).iter().map(|&v| { let p = pt(&m, v); ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2)).sqrt() - 0.25 }).collect();
            let (lo, hi) = (d.iter().cloned().fold(f64::INFINITY, f64::min), d.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
            assert!(lo >= -1e-8 || hi <= 1e-8, "h = {h}: cell {c} straddles the sphere: {lo}..{hi}");
        }
        errs.push((domain_volume(&m, 2) - exact).abs() / exact);
        assert!((domain_volume(&m, 1) + domain_volume(&m, 2) - 1.0).abs() < 1e-12);
    }
    assert!(errs[1] < 0.5 * errs[0] && errs[1] < 0.02, "{errs:?}");
}

#[test]
fn sharp_edges_off_the_lattice_are_fitted() {
    // an L-shaped solid whose re-entrant edge (x = 0.43, z = 0.5) is not on the lattice
    let make = |snap: bool| {
        let mut s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 0.5], 1), box_body([0.0, 0.0, 0.5], [0.42, 1.0, 1.0], 1)], 0.1);
        s.snap_features = snap;
        s
    };
    let exact = 0.5 + 0.42 * 0.5;
    let (fitted, _) = mesh_3d(&make(true)).unwrap();
    let (plain, _) = mesh_3d(&make(false)).unwrap();
    check_conforming(&fitted);
    check_conforming(&plain);
    let (ef, ep) = ((domain_volume(&fitted, 1) - exact).abs() / exact, (domain_volume(&plain, 1) - exact).abs() / exact);
    assert!(ef < 1e-9 || ef < 0.3 * ep, "fitted error {ef}, plain error {ep}");
    // the fitted mesh has vertices on the crease
    let on_crease = (0..fitted.n_points()).filter(|&v| { let p = pt(&fitted, v as u32); (p[0] - 0.42).abs() < 1e-9 && (p[2] - 0.5).abs() < 1e-9 }).count();
    assert!(on_crease >= 5, "vertices on the re-entrant edge: {on_crease}");
}

fn max_edge_in(m: &SimplexMesh, domain: i32) -> f64 {
    let mut worst = 0.0f64;
    for c in 0..m.n_cells() {
        if m.cell_labels[c] != domain {
            continue;
        }
        let v = m.cell(c);
        for i in 0..4 {
            for j in i + 1..4 {
                let d = sub(pt(m, v[i]), pt(m, v[j]));
                worst = worst.max(dot(d, d).sqrt());
            }
        }
    }
    worst
}

#[test]
fn a_body_with_its_own_element_size_is_refined_locally() {
    // domain 2 wants elements of 0.05 in a mesh of 0.2: a graded, conforming mesh with exact volumes
    let mut s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1), box_body([0.3, 0.3, 0.3], [0.7, 0.7, 0.7], 2)], 0.2);
    let (coarse, _) = mesh_3d(&s).unwrap();
    s.bodies[1].h = Some(0.05);
    let (m, rep) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert!(rep.refined > 0);
    assert!((domain_volume(&m, 2) - 0.064).abs() < 1e-12 && (domain_volume(&m, 1) - 0.936).abs() < 1e-12);
    let count = |mm: &SimplexMesh, d: i32| mm.cell_labels.iter().filter(|&&l| l == d).count();
    assert!(count(&m, 2) > 20 * count(&coarse, 2), "refined body: {} cells, was {}", count(&m, 2), count(&coarse, 2));
    assert!(max_edge_in(&m, 2) <= 1.3 * 0.05 * 1.0001, "longest edge in the refined body {}", max_edge_in(&m, 2));
    // the other body stays coarse away from the refined one
    assert!(count(&m, 1) < 6 * count(&coarse, 1), "{} against {}", count(&m, 1), count(&coarse, 1));
    // far from the refined body the cells are as large as before
    let far_coarse = (0..m.n_cells()).filter(|&c| {
        let v = m.cell(c);
        let cen = [0, 1, 2].map(|k| v.iter().map(|&i| pt(&m, i)[k]).sum::<f64>() / 4.0);
        cen[0] < 0.1 && m.cell_labels[c] == 1 && (0..4).all(|i| (i + 1..4).all(|j| { let d = sub(pt(&m, v[i]), pt(&m, v[j])); dot(d, d).sqrt() > 0.1 }))
    }).count();
    assert!(far_coarse > 20, "cells of the original size far from the refined body: {far_coarse}");
}

#[test]
fn local_refinement_keeps_curved_surfaces_on_the_surface() {
    // a refined ball in a coarse block: the new vertices on the ball move onto the sphere
    let r = 0.25;
    let exact = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
    let mut s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 1), (Solid::Sphere { c: [0.5, 0.5, 0.5], r }, 2)], 0.2);
    let (coarse, _) = mesh_3d(&s).unwrap();
    s.bodies[1].h = Some(0.06);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    let err = |mm: &SimplexMesh| (domain_volume(mm, 2) - exact).abs() / exact;
    assert!(err(&m) < 0.3 * err(&coarse) && err(&m) < 0.02, "volume error {} against {}", err(&m), err(&coarse));
    assert!((domain_volume(&m, 1) + domain_volume(&m, 2) - 1.0).abs() < 1e-12);
    // every vertex shared by both domains lies on the sphere
    let mut in1 = vec![false; m.n_points()];
    let mut in2 = vec![false; m.n_points()];
    for c in 0..m.n_cells() {
        for &v in m.cell(c) {
            if m.cell_labels[c] == 1 { in1[v as usize] = true } else { in2[v as usize] = true }
        }
    }
    let mut worst = 0.0f64;
    let mut on_sphere = 0;
    for v in 0..m.n_points() {
        if in1[v] && in2[v] {
            let p = pt(&m, v as u32);
            worst = worst.max((dot(sub(p, [0.5; 3]), sub(p, [0.5; 3])).sqrt() - r).abs());
            on_sphere += 1;
        }
    }
    assert!(on_sphere > 100 && worst < 1e-8, "{on_sphere} interface vertices, worst distance from the sphere {worst}");
}

#[test]
fn a_thin_film_gets_layers_across_its_thickness() {
    // a film of 0.05 on a substrate with a lattice size of 0.2: at least two cell layers across the film
    let s = spec(vec![box_body([0.0, 0.0, 0.0], [1.0, 1.0, 0.6], 1), box_body([0.0, 0.0, 0.6], [1.0, 1.0, 0.65], 2)], 0.2);
    let (m, _) = mesh_3d(&s).unwrap();
    check_conforming(&m);
    assert!((domain_volume(&m, 2) - 0.05).abs() < 1e-12 && (domain_volume(&m, 1) - 0.6).abs() < 1e-12);
    let mut zs: Vec<f64> = vec![];
    for c in 0..m.n_cells() {
        if m.cell_labels[c] == 2 {
            for &v in m.cell(c) {
                zs.push((pt(&m, v)[2] * 1e9).round() / 1e9);
            }
        }
    }
    zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    zs.dedup();
    assert!(zs.len() >= 3, "vertex layers in the film: {zs:?}");
    let mut off = s.clone();
    off.min_layers = 0.0;
    let (plain, _) = mesh_3d(&off).unwrap();
    assert!(m.n_cells() > plain.n_cells());
}
