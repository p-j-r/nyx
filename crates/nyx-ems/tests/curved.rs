//! The `ems` cavity kernel on curved (isoparametric) tetrahedra, ADR 82: closed forms only. A mesh with every edge node at the midpoint
//! gives the matrices of the straight mesh; the curl of the basis is the physical curl of its field; a constant field is represented
//! exactly on a mesh that is a quadratic image of a box.
mod common;
use femcore::fem::{quad_rule_collapsed, CurvedTet};
use femcore::linalg::csr::{c, C};
use femcore::mesh::SimplexMesh;
use nyx_ems::cavity::{assemble_maxwell, basis3, Ned3, CURVED_QUAD_DEGREE};

type C3 = [[C; 3]; 3];

fn unit_media(n: usize) -> (Vec<C3>, Vec<[C; 3]>) {
    let mut e = [[c(0.0, 0.0); 3]; 3];
    for (a, row) in e.iter_mut().enumerate() {
        row[a] = c(1.0, 0.0);
    }
    (vec![e; n], vec![[c(1.0, 0.0); 3]; n])
}

fn ned(mesh: &SimplexMesh, order: usize) -> Ned3 {
    Ned3::new(mesh, order, (0..mesh.n_points() as u32).collect())
}

/// Puts every edge node at the midpoint of its edge: the curved path of the assembly on straight geometry.
fn with_midpoint_nodes(m: &SimplexMesh) -> SimplexMesh {
    let mut out = SimplexMesh::new(
        3,
        m.points.clone(),
        m.cells.clone(),
        m.cell_labels.clone(),
        m.facets.clone(),
        m.facet_labels.clone(),
    );
    let nodes: Vec<([u32; 2], [f64; 3])> = m
        .edges()
        .iter()
        .map(|&(a, b)| {
            let (pa, pb) = (m.point(a as usize), m.point(b as usize));
            ([a, b], [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0])
        })
        .collect();
    out.set_edge_nodes(&nodes);
    out
}

fn rel_diff(a: &femcore::linalg::csr::CCsr, b: &femcore::linalg::csr::CCsr) -> f64 {
    let n = a.nrows();
    let x: Vec<C> =
        (0..n).map(|i| c(((i * 7919) % 101) as f64 / 101.0 - 0.5, ((i * 104729) % 53) as f64 / 53.0 - 0.5)).collect();
    let (ya, yb) = (a.matvec(&x), b.matvec(&x));
    let num: f64 = ya.iter().zip(&yb).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
    let den: f64 = ya.iter().map(|p| p.norm_sqr()).sum::<f64>().sqrt();
    num / den
}

#[test]
fn midpoint_edge_nodes_give_the_matrices_of_the_straight_mesh() {
    let straight = common::block(0.04, 0.03, 0.02, 2, 2, 2);
    let curved = with_midpoint_nodes(&straight);
    assert!(curved.has_curved() && (0..curved.n_cells()).all(|c| curved.cell_is_curved(c)));
    let (eps, mu) = unit_media(straight.n_cells());
    for order in [1, 2] {
        let (k0, m0) = assemble_maxwell(&straight, &ned(&straight, order), &eps, &mu);
        let (k1, m1) = assemble_maxwell(&curved, &ned(&curved, order), &eps, &mu);
        assert!(rel_diff(&k0, &k1) < 1e-11, "order {order}: K differs by {}", rel_diff(&k0, &k1));
        assert!(rel_diff(&m0, &m1) < 1e-11, "order {order}: M differs by {}", rel_diff(&m0, &m1));
    }
}

/// Inverse of the quadratic map by Newton iteration: the barycentric point whose image is `p`.
fn inverse_map(ct: &CurvedTet, p: [f64; 3], start: [f64; 4]) -> [f64; 4] {
    let mut lam = start;
    for _ in 0..30 {
        let q = ct.position(&lam);
        let r = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let (g, _) = ct.jacobian(&lam);
        let mut d = [0.0; 3];
        for k in 0..3 {
            d[k] = (0..3).map(|i| g[k + 1][i] * r[i]).sum();
        }
        for k in 0..3 {
            lam[k + 1] -= d[k];
        }
        lam[0] = 1.0 - lam[1] - lam[2] - lam[3];
        if d.iter().all(|x| x.abs() < 1e-15) {
            break;
        }
    }
    lam
}

#[test]
fn the_curl_of_the_basis_is_the_curl_of_its_field_on_a_curved_cell() {
    let v = [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.1, 1.0, 0.2], [0.0, 0.2, 1.0]];
    let mut x = CurvedTet::straight(v).x;
    let d = [
        [0.03, -0.02, 0.01],
        [-0.01, 0.04, 0.02],
        [0.02, 0.01, -0.03],
        [0.0, -0.03, 0.02],
        [0.03, 0.03, 0.0],
        [-0.02, 0.0, 0.04],
    ];
    for e in 0..6 {
        for k in 0..3 {
            x[4 + e][k] += d[e][k];
        }
    }
    let ct = CurvedTet::new(x);
    let gid = [3u32, 1, 4, 2];
    for order in [1usize, 2] {
        let nb = if order == 1 { 6 } else { 20 };
        let coef: Vec<f64> = (0..nb).map(|s| ((s * 37 + 11) % 17) as f64 / 17.0 - 0.5).collect();
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let field = |lam: [f64; 4], w: &mut Vec<[f64; 3]>, cu: &mut Vec<[f64; 3]>| -> ([f64; 3], [f64; 3]) {
            let (g, _) = ct.jacobian(&lam);
            basis3(order, &lam, &g, &gid, w, cu);
            let (mut e, mut cr) = ([0.0; 3], [0.0; 3]);
            for s in 0..nb {
                for k in 0..3 {
                    e[k] += coef[s] * w[s][k];
                    cr[k] += coef[s] * cu[s][k];
                }
            }
            (e, cr)
        };
        let lam0 = [0.3, 0.2, 0.25, 0.25];
        let p0 = ct.position(&lam0);
        let (_, exact) = field(lam0, &mut w, &mut cu);
        let h = 1e-5;
        let mut de = [[0.0; 3]; 3]; // de[j][i] = d E_i / d x_j
        for j in 0..3 {
            let (mut pp, mut pm) = (p0, p0);
            pp[j] += h;
            pm[j] -= h;
            let (lp, lm) = (inverse_map(&ct, pp, lam0), inverse_map(&ct, pm, lam0));
            let (ep, _) = field(lp, &mut w, &mut cu);
            let (em, _) = field(lm, &mut w, &mut cu);
            for i in 0..3 {
                de[j][i] = (ep[i] - em[i]) / (2.0 * h);
            }
        }
        let fd = [de[1][2] - de[2][1], de[2][0] - de[0][2], de[0][1] - de[1][0]];
        for k in 0..3 {
            assert!(
                (fd[k] - exact[k]).abs() < 1e-6 * (1.0 + exact[k].abs()),
                "order {order} component {k}: finite difference {} vs basis curl {}",
                fd[k],
                exact[k]
            );
        }
    }
}

/// A smooth quadratic map of the unit cube and its closed-form volume.
fn phi(p: [f64; 3]) -> [f64; 3] {
    let (a, a2, b, cc) = (0.10, 0.08, 0.06, 0.05);
    [p[0] + a * p[0] * p[0] + b * p[1] * p[2], p[1] + a2 * p[2] * p[2], p[2] + cc * p[0] * p[1]]
}
fn phi_volume() -> f64 {
    let (a, a2, b, cc) = (0.10, 0.08, 0.06, 0.05);
    1.0 + a - a2 * cc / 2.0 - (2.0 / 3.0) * a * a2 * cc + a2 * b * cc / 3.0 - b * cc / 3.0
}

#[test]
fn a_constant_field_is_in_the_second_order_space_on_a_quadratic_image_of_a_box() {
    let box_mesh = common::block(1.0, 1.0, 1.0, 2, 2, 2);
    let pts: Vec<f64> = box_mesh.points.chunks(3).flat_map(|p| phi([p[0], p[1], p[2]])).collect();
    let mut mesh = SimplexMesh::new(
        3,
        pts,
        box_mesh.cells.clone(),
        box_mesh.cell_labels.clone(),
        box_mesh.facets.clone(),
        box_mesh.facet_labels.clone(),
    );
    let nodes: Vec<([u32; 2], [f64; 3])> = box_mesh
        .edges()
        .iter()
        .map(|&(a, b)| {
            let (pa, pb) = (box_mesh.point(a as usize), box_mesh.point(b as usize));
            ([a, b], phi([(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0]))
        })
        .collect();
    mesh.set_edge_nodes(&nodes);
    // the cells fill the image of the box exactly
    let vol: f64 = (0..mesh.n_cells()).map(|ci| CurvedTet::new(mesh.cell_geometry_nodes(ci)).volume()).sum();
    assert!((vol - phi_volume()).abs() < 1e-12, "volume {vol} vs {}", phi_volume());
    let order = 2;
    let nd = ned(&mesh, order);
    let (eps, mu) = unit_media(mesh.n_cells());
    let (k, m) = assemble_maxwell(&mesh, &nd, &eps, &mu);
    // projection of the constant field cvec on the space: x = M^-1 b; it is exact when the constant is in the space
    let cvec = [0.3, -0.5, 0.8];
    let mut b = vec![0.0; nd.n_dofs()];
    let rule = quad_rule_collapsed(CURVED_QUAD_DEGREE + 2);
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    for ci in 0..mesh.n_cells() {
        let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
        let gid = nd.keys(&mesh, ci);
        for (l, wq) in rule.lam.iter().zip(&rule.w) {
            let lam = [l[0], l[1], l[2], l[3]];
            let (g, det) = ct.jacobian(&lam);
            basis3(order, &lam, &g, &gid, &mut w, &mut cu);
            for (s, ws) in w.iter().enumerate() {
                b[nd.cell_dofs[ci][s]] += wq * det.abs() / 6.0 * (ws[0] * cvec[0] + ws[1] * cvec[1] + ws[2] * cvec[2]);
            }
        }
    }
    let mr = m.real_part();
    let solver = femcore::linalg::driver::DirectSolver::new(&mr, true).unwrap();
    let xs = solver.solve(&b);
    let energy: f64 = xs.iter().zip(&b).map(|(p, q)| p * q).sum();
    let want = vol * (cvec[0] * cvec[0] + cvec[1] * cvec[1] + cvec[2] * cvec[2]);
    assert!((energy / want - 1.0).abs() < 1e-10, "projected norm {energy} vs {want}: the constant is not in the space");
    // its curl is zero: K x = 0 relative to M x
    let xc: Vec<C> = xs.iter().map(|&v| c(v, 0.0)).collect();
    let kx = k.matvec(&xc);
    let mx = m.matvec(&xc);
    let (nk, nm): (f64, f64) =
        (kx.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt(), mx.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt());
    assert!(nk < 1e-9 * nm, "curl of a constant field: |Kx| = {nk}, |Mx| = {nm}");
}

// ------------------------------------------------------------------ PEC sphere cavity

const C0: f64 = 299_792_458.0;

/// Mean relative error of the three modes of the lowest triplet (k r = 2.7437) of a PEC sphere of radius `r`, order 2.
fn sphere_error(r: f64, h: f64, curved: bool) -> (f64, usize) {
    let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(h);
    opts.curved = curved;
    let mesh = femcore::meshing::geomdrive3::sphere_cavity_with(r, &opts).unwrap().mesh;
    assert_eq!(mesh.has_curved(), curved);
    let exact = C0 * 2.7437 / (2.0 * std::f64::consts::PI * r);
    let mut p = common::emw();
    p.order = 2;
    p.materials.default.insert("eps_r".into(), 1.0);
    p.materials.default.insert("mu_r".into(), 1.0);
    p.study = nyx_ems::Study::Cavity { n_modes: 3, f_search: exact };
    let res = common::solve(&mesh, &p).unwrap();
    let f: Vec<f64> = res.solutions.iter().map(|s| s.value.unwrap()).collect();
    assert_eq!(f.len(), 3, "{f:?}");
    (f.iter().map(|v| (v / exact - 1.0).abs()).sum::<f64>() / 3.0, mesh.n_cells())
}

#[test]
#[ignore = "slow: four sphere meshes, straight and curved"]
fn curved_edges_raise_the_convergence_rate_of_the_pec_sphere_cavity() {
    let r = 0.05;
    let hs = [r / 2.5, r / 3.5, r / 5.0, r / 7.0];
    let mut rows = vec![];
    for &h in &hs {
        let (es, ns) = sphere_error(r, h, false);
        let (ec, nc) = sphere_error(r, h, true);
        eprintln!("h = r/{:.1}: straight {es:.3e} ({ns} cells), curved {ec:.3e} ({nc} cells)", r / h);
        rows.push((es, ec));
    }
    let rate = |a: f64, b: f64, h0: f64, h1: f64| (a / b).ln() / (h0 / h1).ln();
    let (rs, rc) = (rate(rows[1].0, rows[3].0, hs[1], hs[3]), rate(rows[1].1, rows[3].1, hs[1], hs[3]));
    eprintln!("rate over the finer sizes: straight {rs:.2}, curved {rc:.2}");
    assert!(rc > rs + 0.5, "curved rate {rc} vs straight {rs}");
    assert!(rows[1].1 * 5.0 < rows[1].0, "at the middle mesh curved {} vs straight {}", rows[1].1, rows[1].0);
}
