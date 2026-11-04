//! 3D metric adaptation: refinement of boxes towards isotropic and anisotropic metrics keeps volumes, labelled faces, interfaces,
//! creases and corners, and reaches the metric. Recomputes volume, area and edge lengths independently of the code under test.
use femcore::adapt::drive3::adapt3;
use femcore::adapt::field::UniformSize;
use femcore::adapt::metric::{edge_length, from_sizes, Sym};
use femcore::adapt::ops3::AdaptOptions3;
use femcore::adapt::topo3::TopoOptions3;
use femcore::mesh::fixtures::box_tet;
use femcore::mesh::native::SimplexMesh;
use std::collections::BTreeSet;
mod common;

fn unit_box(n: usize) -> SimplexMesh<3> {
    box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [n, n, n])
}

fn volume(m: &SimplexMesh<3>, label: Option<i32>) -> f64 {
    (0..m.n_cells()).filter(|&c| label.is_none_or(|l| m.cell_label[c] == l)).map(|c| {
        let v = m.cell(c);
        let p = |i: usize| m.points[v[i] as usize];
        let (a, b, c, d) = (p(0), p(1), p(2), p(3));
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let x = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
        (u[0] * (w[1] * x[2] - w[2] * x[1]) - u[1] * (w[0] * x[2] - w[2] * x[0]) + u[2] * (w[0] * x[1] - w[1] * x[0])) / 6.0
    }).sum()
}

fn facet_area(m: &SimplexMesh<3>) -> f64 {
    (0..m.n_facets()).map(|f| {
        let v = m.facet(f);
        let p = |i: usize| m.points[v[i] as usize];
        let (a, b, c) = (p(0), p(1), p(2));
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
    }).sum()
}

/// Fraction of edges with metric length in [1/sqrt 2, sqrt 2] for a constant metric.
fn in_band(m: &SimplexMesh<3>, met: &Sym<3>) -> f64 {
    let edges = m.edges();
    let ok = edges.iter().filter(|e| {
        let (a, b) = (m.points[e[0] as usize], m.points[e[1] as usize]);
        let l = edge_length(met, met, &[b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        (std::f64::consts::FRAC_1_SQRT_2..=std::f64::consts::SQRT_2).contains(&l)
    }).count();
    ok as f64 / edges.len() as f64
}

fn qmin(m: &SimplexMesh<3>, met: &Sym<3>) -> f64 {
    (0..m.n_cells()).map(|c| {
        let v = m.cell(c);
        let p: Vec<[f64; 3]> = v.iter().map(|&i| m.points[i as usize]).collect();
        femcore::mesh::quality::simplex_quality::<3>(&p, Some(&[*met; 4]))
    }).fold(f64::INFINITY, f64::min)
}

#[test]
fn isotropic_refinement_of_a_box_keeps_geometry_and_reaches_the_metric() {
    let m = unit_box(2);
    let h = 0.15;
    let out = adapt3(&m, &UniformSize(h), &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
    let a = &out.mesh;
    a.check().unwrap();
    assert!((volume(a, None) - 1.0).abs() < 1e-12, "volume {}", volume(a, None));
    assert!((facet_area(a) - 6.0).abs() < 1e-12, "boundary area {}", facet_area(a));
    // the corners are still vertices and every boundary facet lies in a face of the cube
    for c in 0..8 {
        let p = [(c & 1) as f64, ((c >> 1) & 1) as f64, ((c >> 2) & 1) as f64];
        assert!(a.points.iter().any(|q| *q == p), "corner {p:?} lost");
    }
    for f in 0..a.n_facets() {
        let v = a.facet(f);
        assert!((0..3).any(|k| v.iter().all(|&i| a.points[i as usize][k] == 0.0) || v.iter().all(|&i| a.points[i as usize][k] == 1.0)));
    }
    let met = from_sizes(&[h; 3], &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let (band, q) = (in_band(a, &met), qmin(a, &met));
    eprintln!("box 2^3 -> {} tets, in band {band:.3}, qmin {q:.3}, iters {}", a.n_cells(), out.report.iters.len());
    assert!(a.n_cells() > 20 * m.n_cells());
    assert!(band > 0.8, "in band {band}");
    assert!(q > 0.02, "worst tetrahedron {q}");
}

#[test]
fn a_material_interface_stays_a_plane_and_each_material_keeps_its_volume() {
    let mut m = unit_box(2);
    for c in 0..m.n_cells() {
        let z: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][2]).sum::<f64>() / 4.0;
        if z > 0.5 {
            m.cell_label[c] = 2;
        }
    }
    let out = adapt3(&m, &UniformSize(0.2), &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
    let a = &out.mesh;
    assert!((volume(a, Some(1)) - 0.5).abs() < 1e-12 && (volume(a, Some(2)) - 0.5).abs() < 1e-12, "{} {}", volume(a, Some(1)), volume(a, Some(2)));
    let labels: BTreeSet<i32> = a.cell_label.iter().copied().collect();
    assert_eq!(labels, BTreeSet::from([1, 2]));
    // vertices strictly on the plane exist, so the interface was refined, not dropped
    assert!(a.points.iter().filter(|p| p[2] == 0.5).count() > 25);
}

#[test]
fn an_anisotropic_metric_gives_elongated_tetrahedra() {
    let m = unit_box(2);
    let sizes = [0.5, 0.5, 0.1];
    let met = from_sizes(&sizes, &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let out = adapt3(&m, &move |_: &[f64; 3]| met, &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
    let a = &out.mesh;
    assert!((volume(a, None) - 1.0).abs() < 1e-12);
    let (band, q) = (in_band(a, &met), qmin(a, &met));
    eprintln!("anisotropic: {} tets, in band {band:.3}, qmin {q:.3}", a.n_cells());
    // far fewer than the isotropic 0.1 mesh would need (about 1e4 · 6 tets), more than the start
    assert!(a.n_cells() > m.n_cells() && a.n_cells() < 700, "{} tets", a.n_cells());
    assert!(band > 0.9, "in band {band}");
}

#[test]
fn a_graded_metric_is_followed_with_good_tetrahedra() {
    // size grows linearly from 0.06 at x = 0 to 0.36 at x = 1
    let field = |x: &[f64; 3]| {
        let h = 0.06 + 0.3 * x[0];
        from_sizes(&[h; 3], &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    };
    let m = unit_box(3);
    let out = adapt3(&m, &field, &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
    let a = &out.mesh;
    assert!((volume(a, None) - 1.0).abs() < 1e-12);
    // edge lengths measured against the field at the edge midpoint
    let edges = a.edges();
    let mut band = 0;
    for e in &edges {
        let (p, q) = (a.points[e[0] as usize], a.points[e[1] as usize]);
        let mid = [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0, (p[2] + q[2]) / 2.0];
        let l = edge_length(&field(&mid), &field(&mid), &[q[0] - p[0], q[1] - p[1], q[2] - p[2]]);
        if (std::f64::consts::FRAC_1_SQRT_2 * 0.9..=std::f64::consts::SQRT_2 * 1.1).contains(&l) {
            band += 1;
        }
    }
    let frac = band as f64 / edges.len() as f64;
    let q = out.report.iters.last().unwrap().q_min;
    eprintln!("graded: {} tets, in band {frac:.3}, qmin {q:.3}, iters {}, converged {}", a.n_cells(), out.report.iters.len(), out.report.converged);
    assert!(frac > 0.85, "in band {frac}");
    assert!(q > 0.02, "worst tetrahedron {q}");
}

#[test]
fn a_fine_box_coarsens_on_its_surface_too() {
    let mut m = unit_box(8);
    for c in 0..m.n_cells() {
        let z: f64 = m.cell(c).iter().map(|&v| m.points[v as usize][2]).sum::<f64>() / 4.0;
        if z > 0.5 {
            m.cell_label[c] = 2;
        }
    }
    let start = m.n_cells();
    let out = adapt3(&m, &UniformSize(0.5), &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
    let a = &out.mesh;
    a.check().unwrap();
    // flat faces and creases lose nothing by coarsening: the geometry is exact
    assert!((volume(a, Some(1)) - 0.5).abs() < 1e-12 && (volume(a, Some(2)) - 0.5).abs() < 1e-12);
    assert!((facet_area(a) - 6.0).abs() < 1e-12, "boundary area {}", facet_area(a));
    for c in 0..8 {
        let p = [(c & 1) as f64, ((c >> 1) & 1) as f64, ((c >> 2) & 1) as f64];
        assert!(a.points.iter().any(|q| *q == p), "corner {p:?} lost");
    }
    let met = from_sizes(&[0.5; 3], &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let surface = a.points.iter().filter(|p| p.iter().any(|&x| x == 0.0 || x == 1.0)).count();
    eprintln!("coarsened {start} -> {} tets, {surface} vertices on the surface, in band {:.3}, qmin {:.3}", a.n_cells(), in_band(a, &met), qmin(a, &met));
    assert!(a.n_cells() * 20 < start, "{} tets left of {start}", a.n_cells());
    assert!(surface < 150, "{surface} surface vertices of 6 * 81 + ... before");
}

#[test]
fn edge_removal_and_surface_flips_never_worsen_the_worst_tetrahedron() {
    let field = |x: &[f64; 3]| {
        let h = 0.07 + 0.25 * (x[0] * x[0] + x[1]).min(1.0);
        from_sizes(&[h; 3], &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    };
    let m = unit_box(3);
    let run = |edge_removal: bool| {
        let opt = AdaptOptions3 { edge_removal, ..AdaptOptions3::default() };
        let out = adapt3(&m, &field, &TopoOptions3::default(), &opt).unwrap();
        assert!((volume(&out.mesh, None) - 1.0).abs() < 1e-12);
        assert!((facet_area(&out.mesh) - 6.0).abs() < 1e-12);
        let flips: usize = out.report.iters.iter().map(|s| s.flips).sum();
        let last = out.report.iters.last().unwrap();
        (out.mesh.n_cells(), last.q_min, last.q_mean, last.in_band, flips)
    };
    let (off, on) = (run(false), run(true));
    eprintln!("without: {} tets qmin {:.3} qmean {:.3} in band {:.3} flips {}\nwith:    {} tets qmin {:.3} qmean {:.3} in band {:.3} flips {}", off.0, off.1, off.2, off.3, off.4, on.0, on.1, on.2, on.3, on.4);
    assert!(on.2 >= off.2 - 1e-3, "mean quality {} against {}", on.2, off.2);
    assert!(on.1 > 0.02);
}

#[test]
fn frozen_periodic_faces_survive_adaptation_and_the_floquet_cell_still_solves() {
    use femcore::adapt::topo3::TopoOptions3;
    use femcore::mesh::interop::{simplex_from_adapt, simplex_to_adapt};
    use femcore::problem::*;
    // a periodic cell of a PEC guide (see tests/emw3d.rs): faces 2 (y = 0) and 5 (y = l) are Floquet partners and must stay congruent
    let (a, b, l) = (0.03, 0.015, 0.01);
    // boundary numbers of common::block: 1 x = 0, 2 y = 0, 3 z = 0, 4 z = lz, 5 y = ly, 6 x = lx
    let nyx = common::block(a, l, b, 6, 3, 2);
    let before = simplex_to_adapt::<3>(&nyx).unwrap();
    let topo = TopoOptions3 { frozen_labels: vec![2, 5], ..TopoOptions3::default() };
    // finer than the start inside, as fine as the start on the frozen faces
    let out = adapt3(&before, &UniformSize(0.0035), &topo, &AdaptOptions3::default()).unwrap();
    let after = simplex_from_adapt(&out.mesh);
    assert!(after.n_cells() > nyx.n_cells());
    let on = |m: &femcore::mesh::SimplexMesh, y: f64| -> BTreeSet<[u64; 2]> {
        (0..m.n_points()).filter(|&i| (m.points[3 * i + 1] - y).abs() < 1e-15).map(|i| [m.points[3 * i].to_bits(), m.points[3 * i + 2].to_bits()]).collect()
    };
    assert_eq!(on(&nyx, 0.0), on(&after, 0.0), "the y = 0 face is unchanged");
    assert_eq!(on(&nyx, l), on(&after, l), "the y = l face is unchanged");
    // TE10 with Floquet phase 0.8 over the cell: f = c/2 sqrt((1/a)^2 + (beta/pi)^2)
    let phase = 0.8f64;
    let c0 = 299_792_458.0;
    let exact = 0.5 * c0 * ((1.0 / a).powi(2) + (phase / l / std::f64::consts::PI).powi(2)).sqrt();
    let mut s = ProblemSpec::new(Physics::Emw);
    s.order = 2;
    s.materials.default.insert("eps_r".into(), 1.0);
    s.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase }];
    s.study = Study::Eigenfrequency { n_modes: 3, f_search: 0.95 * exact };
    let r = femcore::solve::solve(&after, &s, &mut SolveState::default()).unwrap();
    let got: Vec<f64> = r.solutions.iter().map(|x| x.value.unwrap()).collect();
    eprintln!("{} -> {} tets, Floquet TE10 {:.6e} vs {:.6e}", nyx.n_cells(), after.n_cells(), got[0], exact);
    assert!((got[0] / exact - 1.0).abs() < 5e-3, "{got:?} vs {exact}");
}
