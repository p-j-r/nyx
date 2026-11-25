// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Curved (isoparametric, quadratic geometry) tetrahedra, ADR 82: the map of `fem::CurvedTet`, the quadrature rule built for it and the
//! storage of curved edges in the meshes. Closed forms only.
#![allow(clippy::needless_range_loop)]
use femcore::fem::{geometry, quad_rule_collapsed, CurvedTet};
use femcore::mesh::cells::{CellMesh, CellType};
use femcore::mesh::{Mesh, SimplexMesh};

fn tet() -> [[f64; 3]; 4] {
    [[0.1, -0.2, 0.05], [1.3, 0.1, -0.1], [0.2, 1.1, 0.3], [-0.1, 0.2, 0.9]]
}

#[test]
fn a_curved_tet_with_midpoint_nodes_is_the_straight_tet() {
    let v = tet();
    let ct = CurvedTet::straight(v);
    let mesh = SimplexMesh::new(3, v.iter().flatten().copied().collect(), vec![0, 1, 2, 3], vec![1], vec![1, 2, 3], vec![1]);
    let geo = geometry(&mesh);
    for lam in [[0.25; 4], [0.7, 0.1, 0.1, 0.1], [0.0, 0.5, 0.3, 0.2], [0.1, 0.2, 0.3, 0.4]] {
        let (g, det) = ct.jacobian(&lam);
        for a in 0..4 {
            for k in 0..3 {
                assert!((g[a][k] - geo.g[a * 3 + k]).abs() < 1e-13, "gradient {a} {k} at {lam:?}: {} vs {}", g[a][k], geo.g[a * 3 + k]);
            }
        }
        assert!((det.abs() / 6.0 - geo.vol[0]).abs() < 1e-14);
        // the position is the affine one
        let p = ct.position(&lam);
        for k in 0..3 {
            let want: f64 = (0..4).map(|a| lam[a] * v[a][k]).sum();
            assert!((p[k] - want).abs() < 1e-14);
        }
    }
    assert!((ct.volume() - geo.vol[0]).abs() < 1e-14);
}

#[test]
fn the_collapsed_rule_integrates_monomials_exactly() {
    // int over the reference tetrahedron of l0^a l1^b l2^c l3^d = a! b! c! d! 3! / (a + b + c + d + 3)!  (times the volume 1/6 / (1/6) = 1 with weights summing to 1)
    fn fact(n: usize) -> f64 {
        (1..=n).map(|k| k as f64).product()
    }
    for degree in [2usize, 4, 6, 8] {
        let q = quad_rule_collapsed(degree);
        assert!((q.w.iter().sum::<f64>() - 1.0).abs() < 1e-14, "weights sum, degree {degree}");
        assert!(q.w.iter().all(|&w| w > 0.0));
        for a in 0..=degree {
            for b in 0..=degree - a {
                for c in 0..=degree - a - b {
                    let d = degree - a - b - c;
                    let exact = fact(a) * fact(b) * fact(c) * fact(d) * 6.0 / fact(a + b + c + d + 3);
                    let got: f64 = q.lam.iter().zip(&q.w).map(|(l, w)| w * l[0].powi(a as i32) * l[1].powi(b as i32) * l[2].powi(c as i32) * l[3].powi(d as i32)).sum();
                    assert!((got - exact).abs() < 1e-12 * exact.max(1e-3), "degree {degree} monomial {a}{b}{c}{d}: {got} vs {exact}");
                }
            }
        }
    }
}

/// A tetrahedron of the unit simplex whose edge nodes sit on a sphere-like quadratic map: the jacobian of the quadratic geometry is
/// checked against finite differences of `position`.
#[test]
fn the_jacobian_matches_finite_differences_of_the_position() {
    let v = tet();
    let mut x = CurvedTet::straight(v).x;
    // displace the edge nodes by hand
    let d = [[0.03, -0.02, 0.01], [-0.01, 0.04, 0.02], [0.02, 0.01, -0.03], [0.0, -0.03, 0.02], [0.03, 0.03, 0.0], [-0.02, 0.0, 0.04]];
    for e in 0..6 {
        for k in 0..3 {
            x[4 + e][k] += d[e][k];
        }
    }
    let ct = CurvedTet::new(x);
    let lam0 = [0.2, 0.3, 0.1, 0.4];
    let jm = ct.jacobian_matrix(&lam0);
    let h = 1e-6;
    for k in 0..3 {
        // reference coordinate k+1 is lambda_{k+1}; lambda_0 = 1 - sum
        let mut lp = lam0;
        let mut lm = lam0;
        lp[k + 1] += h;
        lp[0] -= h;
        lm[k + 1] -= h;
        lm[0] += h;
        let (pp, pm) = (ct.position(&lp), ct.position(&lm));
        for i in 0..3 {
            let fd = (pp[i] - pm[i]) / (2.0 * h);
            assert!((jm[i][k] - fd).abs() < 1e-8, "J[{i}][{k}] = {} vs {fd}", jm[i][k]);
        }
    }
    // gradients of the barycentric functions invert the map: sum_a grad lambda_a (x_a . e) ... checked through g_k . dx/dxi_j = delta_kj
    let (g, det) = ct.jacobian(&lam0);
    assert!(det > 0.0);
    for kk in 1..4 {
        for j in 0..3 {
            let s: f64 = (0..3).map(|i| g[kk][i] * jm[i][j]).sum();
            assert!((s - if kk == j + 1 { 1.0 } else { 0.0 }).abs() < 1e-12);
        }
    }
    for i in 0..3 {
        assert!((g[0][i] + g[1][i] + g[2][i] + g[3][i]).abs() < 1e-14);
    }
}

#[test]
fn the_volume_of_a_tet_with_a_quadratic_edge_follows_the_closed_form() {
    let v = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    // unit tetrahedron with the node of edge 0-1 moved by 0.05: the map has x_i += 4 s l0 l1, so J gets the entry 4 s (l0 - l1) in column 0
    // of one row. In rows 1 and 2 the matrix stays triangular (det J = 1 exactly); in row 0 the determinant is 1 + 4 s (l0 - l1), whose
    // integral over the cell is zero because l0 and l1 have the same mean. The volume is 1/6 in all three directions.
    // Correction of the comment above: for y or z the determinant is 1 - 4 s l1, whose integral is (1 - s) / 6; for x it is 1 / 6.
    let s = 0.05;
    for (k, want) in [(0usize, 1.0 / 6.0), (1, (1.0 - s) / 6.0), (2, (1.0 - s) / 6.0)] {
        let mut x = CurvedTet::straight(v).x;
        x[4][k] += s;
        let vol = CurvedTet::new(x).volume();
        assert!((vol - want).abs() < 1e-14, "axis {k}: {vol} vs {want}");
    }
}

#[test]
fn edge_nodes_survive_the_exchange_mesh_and_an_empty_list_is_the_straight_mesh() {
    // two tetrahedra sharing a face; the edge 0-1 is curved
    let pts = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let mut cm = CellMesh::<3>::new();
    cm.points = pts.chunks(3).map(|p| [p[0], p[1], p[2]]).collect();
    cm.push_cell(CellType::Tet, &[0, 1, 2, 3], 1);
    cm.push_cell(CellType::Tet, &[1, 4, 2, 3], 1);
    for f in [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1], [1, 4, 2], [4, 3, 2], [1, 3, 4]] {
        cm.push_facet(&f, 1);
    }
    let plain = SimplexMesh::from_cells(&cm).unwrap();
    assert!(!plain.has_curved() && !plain.cell_is_curved(0));
    cm.set_edge_nodes([([1u32, 0], [0.5, 0.0, 0.03]), ([3, 4], [0.5, 0.5, 1.0])]);
    assert_eq!(cm.edge_nodes.len(), 2);
    assert_eq!(cm.edge_node(0, 1), Some([0.5, 0.0, 0.03]));
    assert_eq!(cm.edge_node(0, 2), None);
    let s = SimplexMesh::from_cells(&cm).unwrap();
    assert!(s.has_curved());
    assert!(s.cell_is_curved(0), "cell 0 has the edge 0-1");
    assert!(s.cell_is_curved(1), "cell 1 has the edge 3-4");
    let nodes = s.cell_geometry_nodes(0);
    assert_eq!(nodes[4], [0.5, 0.0, 0.03]);
    assert_eq!(nodes[5], [0.0, 0.5, 0.0], "a straight edge gives its midpoint");
    // round trip through the exchange mesh keeps the list
    let back = Mesh::from_solver_mesh(&s).unwrap();
    let Mesh::D3(b) = &back else { panic!("3D") };
    assert_eq!(b.edge_nodes, cm.edge_nodes);
    let again = back.discretise().unwrap();
    assert_eq!(again.edge_node_list(), s.edge_node_list());
    // the straight mesh does not change
    let straight = Mesh::from_solver_mesh(&plain).unwrap();
    let Mesh::D3(b) = &straight else { panic!("3D") };
    assert!(b.edge_nodes.is_empty());
    // an entry for an edge that does not exist is dropped by the solver mesh and by `retain_existing_edge_nodes`
    cm.set_edge_nodes([([0u32, 4], [0.5, 0.5, 0.5])]);
    assert!(SimplexMesh::from_cells(&cm).unwrap().edge_node_list().is_empty());
    cm.retain_existing_edge_nodes();
    assert!(cm.edge_nodes.is_empty());
}

// ------------------------------------------------------------------ meshed sphere

fn sphere_volume_error(h: f64, curved: bool) -> f64 {
    let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(h);
    opts.curved = curved;
    let g = femcore::meshing::geomdrive3::sphere_cavity_with(1.0, &opts).unwrap();
    let m = &g.mesh;
    assert_eq!(m.has_curved(), curved);
    let vol: f64 = (0..m.n_cells()).map(|c| CurvedTet::new(m.cell_geometry_nodes(c)).volume()).sum();
    (vol / (4.0 / 3.0 * std::f64::consts::PI) - 1.0).abs()
}

#[test]
fn the_volume_of_a_meshed_sphere_converges_at_third_order_with_curved_edges() {
    let hs = [0.5, 0.35, 0.25];
    let straight: Vec<f64> = hs.iter().map(|&h| sphere_volume_error(h, false)).collect();
    let curved: Vec<f64> = hs.iter().map(|&h| sphere_volume_error(h, true)).collect();
    eprintln!("sphere volume error, h = {hs:?}: straight {straight:?}, curved {curved:?}");
    // the rates between the finest two sizes
    let rate = |e: &[f64]| (e[1] / e[2]).ln() / (hs[1] / hs[2]).ln();
    assert!(curved[2] < 1e-3, "curved error at h = r/4: {}", curved[2]);
    assert!(curved[2] < straight[2] / 5.0, "curved {curved:?} vs straight {straight:?}");
    assert!(rate(&curved) > 2.5, "curved rate {}: {curved:?}", rate(&curved));
    assert!(rate(&straight) < 2.6, "straight rate {}: {straight:?}", rate(&straight));
}

#[test]
fn the_curved_edges_do_not_depend_on_the_number_of_threads() {
    let make = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        pool.install(|| {
            let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(0.3);
            opts.curved = true;
            femcore::meshing::geomdrive3::sphere_cavity_with(1.0, &opts).unwrap().mesh
        })
    };
    let (a, b) = (make(1), make(5));
    assert!(a.has_curved());
    assert_eq!(a.edge_node_list(), b.edge_node_list());
    assert_eq!(a.points, b.points);
}

#[test]
fn sheet3d_circular_disk_curved_mid_edge_nodes_converge_area_and_stay_positive_jacobian() {
    use femcore::meshing::request::{mesh_request, MeshRequest, MesherRequest};
    use geom::model::Model;

    let r = 100e-6;
    let mut m = Model::new_3d();
    let sub = m.box_("sub", (-220e-6, -220e-6, -80e-6), (220e-6, 220e-6, 0.0));
    let vac = m.box_("vac", (-220e-6, -220e-6, 0.0), (220e-6, 220e-6, 80e-6));
    let disk_2d = m.circle("disk_2d", (0.0, 0.0), r);
    m.sheet("disk", disk_2d, 0.0);
    m.overlay(&[("substrate", sub), ("vacuum", vac)]);

    let req = MeshRequest {
        h: 55e-6,
        max_dev: Some(2.5e-6),
        mesher: MesherRequest::Simple,
        pml_mode: femcore::meshing::request::PmlModeRequest::Inline,
        structured_pml: None,
        waves: None,
        curved: true,
        gap_cells: None,
        z_coarsen: false,
        max_cells: None,
        dz_min: None,
        max_z_slabs: None,
        z_growth_ratio: None,
        max_steiner_multiplier: None,
        strict_mesh: true,
        singular_edges: None,
        quarter_point_xi: None,
        corner_grading: false,
    };
    let (g, rep) = mesh_request(&m.list, &req).unwrap();
    assert!(g.mesh.has_curved());
    let cs = rep.curved.unwrap();
    assert!(cs.edges > 0 && cs.cells > 0);
    assert!(cs.min_scaled_jacobian > 0.1, "min_scaled_jacobian = {}", cs.min_scaled_jacobian);
}
