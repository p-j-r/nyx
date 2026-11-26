// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `ems` cavity kernel on curved (isoparametric) tetrahedra, ADR 82: closed forms only. A mesh with every edge node at the midpoint
//! gives the matrices of the straight mesh; the curl of the basis is the physical curl of its field; a constant field is represented
//! exactly on a mesh that is a quadratic image of a box.
mod common;
use femcore::fem::{quad_rule_collapsed, CurvedTet};
use femcore::linalg::csr::{c, C};
use femcore::mesh::SimplexMesh;
use nyx_ems::cavity::{assemble_maxwell, basis3, Ned3, CURVED_QUAD_DEGREE};

type C3 = [[C; 3]; 3];

fn unit_media(n: usize) -> (Vec<C3>, Vec<C3>) {
    let mut e = [[c(0.0, 0.0); 3]; 3];
    for (a, row) in e.iter_mut().enumerate() {
        row[a] = c(1.0, 0.0);
    }
    (vec![e; n], vec![e; n])
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

/// Mean relative error of the three modes of the lowest triplet (exact TM101 root k R = 2.743707269992266) of a PEC sphere of radius `r`, order 2.
fn sphere_error(r: f64, h: f64, curved: bool) -> (f64, usize) {
    let mut opts = femcore::meshing::geomdrive::GeomMeshOpts::new(h);
    opts.curved = curved;
    let mesh = femcore::meshing::geomdrive3::sphere_cavity_with(r, &opts).unwrap().mesh;
    assert_eq!(mesh.has_curved(), curved);
    const KR_TM101: f64 = 2.743707269992266;
    let exact = C0 * KR_TM101 / (2.0 * std::f64::consts::PI * r);
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
fn curved_edges_raise_the_convergence_rate_of_the_pec_sphere_cavity() {
    let r = 0.05;
    let hs = [r / 2.5, r / 3.5, r / 5.0, r / 7.0];
    let mut straight_errs = vec![];
    let mut curved_errs = vec![];
    for &h in &hs {
        let (es, ns) = sphere_error(r, h, false);
        let (ec, nc) = sphere_error(r, h, true);
        eprintln!("h = r/{:.1}: straight {es:.3e} ({ns} cells), curved {ec:.3e} ({nc} cells)", r / h);
        straight_errs.push(es);
        curved_errs.push(ec);
    }
    let rep_s = femcore::vvuq::convergence_rate(&hs, &straight_errs).unwrap();
    let rep_c = femcore::vvuq::convergence_rate(&hs, &curved_errs).unwrap();
    let p_straight = rep_s.ols_slope;
    let p_curved = rep_c.ols_slope;
    eprintln!(
        "OLS convergence rate: straight {p_straight:.2} (asym {:.2}), curved {p_curved:.2} (asym {:.2})",
        rep_s.asymptotic_slope, rep_c.asymptotic_slope
    );
    assert!(p_curved >= 3.5, "curved convergence rate {p_curved:.3} < 3.5 (errs {curved_errs:?})");
    assert!(p_curved > p_straight + 1.0, "curved rate {p_curved:.3} must exceed straight {p_straight:.3} + 1.0");
    assert!(
        curved_errs[1] * 5.0 < straight_errs[1],
        "at middle mesh curved {} vs straight {}",
        curved_errs[1],
        straight_errs[1]
    );
}

#[test]
fn r34_curved_spherical_pml_cavity_and_postsolve_metrics() {
    use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
    use femcore::meshing::geomdrive3::mesh_geom_3d;
    use femcore::props::{Pml, PmlCoord};
    use geom::model::{Model, RegionKind};
    use nyx_physics::Scalar;

    let (r_core, r_air, r_pml) = (0.03, 0.045, 0.065);
    let mut m = Model::new_3d();
    let pml_s = m.sphere("pml", (0.0, 0.0, 0.0), r_pml);
    let air_s = m.sphere("air", (0.0, 0.0, 0.0), r_air);
    let core_s = m.sphere("core", (0.0, 0.0, 0.0), r_core);
    m.overlay_full(
        &[
            ("pml", pml_s, RegionKind::Material, None),
            ("air", air_s, RegionKind::Material, None),
            ("core", core_s, RegionKind::Material, None),
        ],
        &[],
    );
    let mut opts = GeomMeshOpts::new(0.020);
    opts.curved = true;
    let g = mesh_geom_3d(evaluate(&m.list).unwrap(), &opts).unwrap();
    assert!(g.mesh.has_curved(), "3-shell spherical cavity+PML must have CurvedTet elements");
    let (d_core, d_air, d_pml) = (g.domain("core").unwrap(), g.domain("air").unwrap(), g.domain("pml").unwrap());

    let solve_pml = |strength: f64| -> (f64, f64, f64) {
        let mut p = common::emw();
        p.order = 2;
        p.materials.set(d_core, "eps_r", 16.0);
        p.materials.set(d_air, "eps_r", 1.0);
        p.materials.set(d_pml, "eps_r", 1.0);
        p.materials.pml.insert(
            d_pml,
            Pml {
                coord: PmlCoord::Spherical { center: [0.0; 3], r_in: r_air, thickness: r_pml - r_air },
                strength,
                order: 2.0,
                re: 1.0,
                ..Pml::default()
            },
        );
        // Dielectric sphere TE101-like Mie resonance: k0 * r_core * sqrt(16) ~ pi => f ~ 1.25 GHz
        p.study = nyx_ems::Study::Cavity { n_modes: 4, f_search: 1.22e9 };
        let res = common::solve(&g.mesh, &p).unwrap();
        let scalar = |name: &str| -> Option<f64> {
            res.scalars.iter().find_map(|(k, v)| match (k == name, v) {
                (true, Scalar::Num(x)) => Some(*x),
                _ => None,
            })
        };
        let mut best: Option<(f64, f64, f64)> = None;
        for (i, sol) in res.solutions.iter().enumerate() {
            let idx = i + 1;
            let lf = scalar(&format!("layer fraction {idx}")).expect("layer fraction emitted");
            assert!(lf.is_finite() && (0.0..=1.0).contains(&lf), "layer_fraction {idx} = {lf}");
            let f_re = sol.value.unwrap();
            let q = scalar(&format!("Q{idx}")).unwrap_or(f64::INFINITY);
            if lf < 0.25 && f_re > 0.5e9 && best.is_none() {
                best = Some((f_re, q, lf));
            }
        }
        best.expect("at least one guided Mie core resonance with layer_fraction < 0.25")
    };

    let (f3, q3, lf3) = solve_pml(3.0);
    let (f5, q5, lf5) = solve_pml(5.0);
    let rel_f = ((f3 - f5) / f3).abs();
    eprintln!("R34 curved spherical PML: s=3 f={f3:.4e} Q={q3:.2} lf={lf3:.3}, s=5 f={f5:.4e} Q={q5:.2} lf={lf5:.3}");
    assert!(q3.is_finite() && q3 > 1.0 && q5.is_finite() && q5 > 1.0, "expected passive radiative Q: {q3}, {q5}");
    assert!(rel_f < 0.02, "curved spherical PML resonant frequency drifted by {rel_f:.3e} across strengths");
}

#[test]
fn fr32_p2_quarter_point_cpw_cavity_equipartition_and_zero_spurious_modes() {
    use femcore::meshing::singularity::{apply_quarter_point_nodes, sheet_rim_vertices, DEFAULT_QUARTER_POINT_XI};
    use femcore::meshing::{mesh_model, MeshPlan};
    use geom::model::Model;
    use nyx_ems::Ems;
    use nyx_physics::{BcSpec, CircuitRom, Ctx, Error, LabelResolver, MaterialAssign, Mesh, Physics, ProblemSpec};
    use serde_json::json;
    use std::collections::BTreeSet;

    // 3D CPW-like strip resonator on Si/air with a lumped_port so CircuitRom & equipartition_error are evaluated.
    let mut m = Model::new_3d();
    let si = m.box_("si_box", (-400e-6, -700e-6, -240e-6), (400e-6, 700e-6, 0.0));
    let air = m.box_("air_box", (-400e-6, -700e-6, 0.0), (400e-6, 700e-6, 240e-6));
    let strip = m.rect("strip_r", (-80e-6, -480e-6), (80e-6, 480e-6));
    m.sheet("strip", strip, 0.0);
    let port_r = m.rect("port_r", (-80e-6, -700e-6), (80e-6, -480e-6));
    m.sheet("port_in", port_r, 0.0);
    m.overlay(&[("si", si), ("air", air)]);

    let mut g = mesh_model(&m.list, &MeshPlan::simple(340e-6)).expect("FR-32 CPW strip mesh");
    let strip_labels = g.boundary("strip");
    let rim_verts = sheet_rim_vertices(&g.mesh, &strip_labels);
    assert!(rim_verts.iter().any(|&r| r), "strip must have singular rim vertices");
    let rep = apply_quarter_point_nodes(&mut g.mesh, &rim_verts, DEFAULT_QUARTER_POINT_XI, &BTreeSet::new());
    assert!(rep.shifted_edges > 0, "expected quarter-point shifted edges > 0");
    assert!(g.mesh.has_qp() && g.mesh.has_curved(), "mesh must carry quarter-point curved tets");

    struct Res<'a>(&'a femcore::meshing::geomdrive::MeshedGeom);
    impl LabelResolver for Res<'_> {
        fn domains(&self, r: &str) -> Result<Vec<i32>, Error> {
            let v = self.0.domains(r);
            if v.is_empty() {
                Err(Error::invalid(format!("no region {r}")))
            } else {
                Ok(v)
            }
        }
        fn boundaries(&self, s: &str) -> Result<Vec<i32>, Error> {
            let v = self.0.boundary(s);
            if !v.is_empty() {
                return Ok(v);
            }
            self.0.select(s).map_err(Error::invalid)
        }
    }
    let exchange = Mesh::from_solver_mesh(&g.mesh).unwrap();
    let resolver = Res(&g);
    let mut spec = ProblemSpec::new("ems", 2, "cavity");
    spec.study.params = json!({ "n_modes": 2, "f_search": 25.0e9 });
    spec.materials.insert("si".into(), MaterialAssign::Props([("eps_r".to_string(), 11.45)].into()));
    spec.materials.insert("air".into(), MaterialAssign::Props([("eps_r".to_string(), 1.0)].into()));
    spec.bcs.push(BcSpec { on: "strip".into(), kind: "pec".into(), params: json!({}) });
    spec.bcs.push(BcSpec {
        on: "port_in".into(),
        kind: "lumped_port".into(),
        params: json!({ "name": "port_in", "Z_0": 50.0, "direction": [0.0, 1.0, 0.0] }),
    });

    let res = Ems.solve(&exchange, &spec, &mut Ctx::new(&resolver)).expect("FR-32 P2 quarter-point cavity solve");
    assert!(!res.solutions.is_empty());
    for sol in &res.solutions {
        let f_re = sol.eigenvalue.expect("eigenvalue present")[0];
        assert!(f_re > 1.0e8, "FR-32 quarter-point P2 cavity emitted spurious mode {f_re:e} Hz");
        for fld in &sol.fields {
            assert!(fld.values.iter().all(|v| v.is_finite()), "non-finite nodal field {}", fld.name);
        }
        for cf in &sol.cell_fields {
            assert!(
                cf.re.iter().all(|v| v.is_finite()) && cf.im.iter().all(|v| v.is_finite()),
                "non-finite cell field {}",
                cf.name
            );
        }
    }
    let eq = res.scalar("equipartition_error").expect("equipartition_error emitted");
    assert!(eq < 1e-9, "FR-32 quarter-point P2 cavity equipartition_error = {eq:e} >= 1e-9");
    let rom = CircuitRom::of(&res).expect("CircuitRom emitted");
    assert_eq!(rom.modes.len(), res.solutions.len());
}
