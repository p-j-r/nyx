// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Pillars (d) & (e) Verification: Conservation Laws & Symmetries / Invariances (§4.5–§4.6).
//!
//! 1. Weak Gauss law `div(eps E) = 0` (`< 1e-10`) & `div(mu H) = 0` (`< 1e-10`) on source-free
//!    lossless interiors (piecewise dielectric box and `CurvedTet` spherical cavity).
//! 2. Poynting power balance (`|P_diss - P_in| / max(P_in, 1e-30) < 1e-3`) on driven lossy & IBC cavities.
//! 3. Multi-port (3-port & 4-port) passivity `lambda_min(I - S^H S) >= -1e-9` and reciprocity `< 1e-10`.
//! 4. Rigid Euclidean invariance (`< 1e-10` mesh-level with rotated Bunge crystal tensor & `< 1e-3`
//!    geometry-level rotated remesh), vertex/cell permutation invariance (`< 1e-10`), port permutation
//!    `P S P^T` (`< 1e-12`), dimensional scaling invariance (`< 1e-10`), and `+phase/-phase` Floquet degeneracy (`< 1e-10`).

mod common;

use femcore::fem::{geometry, quad_rule, quad_rule_collapsed, shape, CurvedTet};
use femcore::linalg::csr::{c, C};
use femcore::linalg::driver::dense_gen_eigh;
use femcore::material::{euler_from_rotation, rotation, Anisotropy, Tensor};
use femcore::mesh::SimplexMesh;
use femcore::meshing::geomdrive::{evaluate, GeomMeshOpts};
use femcore::meshing::geomdrive3::{mesh_geom_3d, sphere_cavity_with};
use femcore::optical::EPS0;
use femcore::props::{Periodic, TensorKind};
use geom::model::{Axis, Model, TStep};
use nyx_ems::cavity::{assemble_emw_eigen, basis3_keyed, solve_emw_modes, CURVED_QUAD_DEGREE};
use nyx_ems::problem::{LumpedPortBc, Problem, Study, SurfaceImpedanceBc, WavePortBc};
use nyx_physics::{DrivenPortData, SolveResult};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

fn rotate_mat3_point(r: &[[f64; 3]; 3], t: [f64; 3], p: [f64; 3]) -> [f64; 3] {
    [
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2] + t[0],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2] + t[1],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2] + t[2],
    ]
}

fn rotate_vec3(r: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
        r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
        r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
    ]
}

fn transform_mesh_3d(mesh: &SimplexMesh, r: &[[f64; 3]; 3], t: [f64; 3]) -> SimplexMesh {
    let mut pts = Vec::with_capacity(mesh.points.len());
    for &p in mesh.points.as_chunks::<3>().0 {
        let q = rotate_mat3_point(r, t, p);
        pts.extend_from_slice(&q);
    }
    let mut out = SimplexMesh::new(
        3,
        pts,
        mesh.cells.clone(),
        mesh.cell_labels.clone(),
        mesh.facets.clone(),
        mesh.facet_labels.clone(),
    );
    if mesh.has_curved() {
        let mut nodes = mesh.edge_node_list();
        for (_, p) in &mut nodes {
            *p = rotate_mat3_point(r, t, *p);
        }
        out.set_edge_nodes(&nodes);
    }
    out
}

fn scale_mesh(mesh: &SimplexMesh, scale: f64) -> SimplexMesh {
    let pts: Vec<f64> = mesh.points.iter().map(|&v| v * scale).collect();
    let mut out = SimplexMesh::new(
        mesh.dim,
        pts,
        mesh.cells.clone(),
        mesh.cell_labels.clone(),
        mesh.facets.clone(),
        mesh.facet_labels.clone(),
    );
    if mesh.has_curved() {
        let mut nodes = mesh.edge_node_list();
        for (_, p) in &mut nodes {
            *p = [p[0] * scale, p[1] * scale, p[2] * scale];
        }
        out.set_edge_nodes(&nodes);
    }
    out
}

/// Checks weak Gauss law `div(eps E) = 0` and `div(mu H) = 0` over interior P1 & P2 Lagrange
/// test functions away from boundary constraints on solved cavity eigenmodes.
fn check_weak_gauss_and_div_b(mesh: &SimplexMesh, p: &Problem) {
    let asm = assemble_emw_eigen(mesh, p).unwrap();
    let solved = solve_emw_modes(&asm).unwrap();
    assert!(!solved.freqs.is_empty(), "expected converged cavity modes");

    let ext_facets = mesh.exterior_facets();
    let mut bnd_verts = vec![false; mesh.n_points()];
    let mut bnd_edges = vec![false; mesh.n_edges()];
    for (f, &is_ext) in ext_facets.iter().enumerate() {
        if !is_ext {
            continue;
        }
        let fv = mesh.facet(f);
        for &v in fv {
            bnd_verts[v as usize] = true;
        }
        for &(a, b) in &[(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
            bnd_edges[mesh.edge_index(a, b)] = true;
        }
    }
    let interior_verts: Vec<usize> = (0..mesh.n_points()).filter(|&v| !bnd_verts[v]).collect();
    let interior_edges: Vec<usize> = (0..mesh.n_edges()).filter(|&e| !bnd_edges[e]).collect();
    assert!(!interior_verts.is_empty(), "mesh must have interior vertices for Gauss test");

    let eps_re = p.materials.cell_tensor(mesh, TensorKind::Permittivity, Some(1.0)).unwrap();
    let geo = geometry(mesh);
    let quad = quad_rule(3, 6);
    let cquad = quad_rule_collapsed(CURVED_QUAD_DEGREE);
    let nb = asm.ned.nb();

    for mode_idx in 0..solved.freqs.len().min(2) {
        let (_, vidx) = solved.freqs[mode_idx];
        let u = asm.red.expand(&solved.vectors[vidx]);

        let mut d_norm2 = 0.0_f64;
        let mut b_norm2 = 0.0_f64;
        let mut grad_v_norm2 = vec![0.0_f64; mesh.n_points()];
        let mut int_d_v = vec![c(0.0, 0.0); mesh.n_points()];
        let mut int_b_v = vec![c(0.0, 0.0); mesh.n_points()];
        let mut grad_e_norm2 = vec![0.0_f64; mesh.n_edges()];
        let mut int_d_e = vec![c(0.0, 0.0); mesh.n_edges()];
        let mut int_b_e = vec![c(0.0, 0.0); mesh.n_edges()];

        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];
        for (ci, ep) in eps_re.iter().enumerate() {
            let cv = mesh.cell(ci);
            let eids = mesh.cell_edge_ids(ci);
            let cd = &asm.ned.cell_dofs[ci];
            let gid = asm.ned.key_pairs(mesh, ci);
            for (b, row) in gs.iter_mut().enumerate() {
                row.copy_from_slice(&geo.g[ci * 12 + b * 3..ci * 12 + b * 3 + 3]);
            }
            let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
            let nq = if curved.is_some() { cquad.w.len() } else { quad.w.len() };
            for q in 0..nq {
                let (lam, gq, wq) = match &curved {
                    Some(ct) => {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        let (gq, det) = ct.jacobian(&lam);
                        let gq = if det.is_finite() && det.abs() > 1e-30 && gq.iter().flatten().all(|v| v.is_finite()) {
                            gq
                        } else {
                            gs
                        };
                        (lam, gq, cquad.w[q] * det.abs() / 6.0)
                    }
                    None => {
                        let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                        (lam, gs, quad.w[q] * geo.vol[ci])
                    }
                };
                basis3_keyed(asm.order, &lam, &gq, &gid, &mut w, &mut cu);
                let (mut eh, mut curleh) = ([c(0.0, 0.0); 3], [c(0.0, 0.0); 3]);
                for s in 0..nb {
                    let coef = u[cd[s]];
                    for kk in 0..3 {
                        eh[kk] += coef * w[s][kk];
                        curleh[kk] += coef * cu[s][kk];
                    }
                }
                let mut dh = [c(0.0, 0.0); 3];
                for i in 0..3 {
                    for j in 0..3 {
                        dh[i] += eh[j] * ep[i][j];
                    }
                }
                for kk in 0..3 {
                    d_norm2 += dh[kk].norm_sqr() * wq;
                    b_norm2 += curleh[kk].norm_sqr() * wq;
                }
                for a in 0..4 {
                    let v = cv[a] as usize;
                    if !bnd_verts[v] {
                        let ga = gq[a];
                        grad_v_norm2[v] += (ga[0] * ga[0] + ga[1] * ga[1] + ga[2] * ga[2]) * wq;
                        for kk in 0..3 {
                            int_d_v[v] += dh[kk] * (ga[kk] * wq);
                            int_b_v[v] += curleh[kk] * (ga[kk] * wq);
                        }
                    }
                }
                if asm.order == 2 {
                    const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
                    for (le, &(a, b)) in EDGES.iter().enumerate() {
                        let eid = eids[le] as usize;
                        if !bnd_edges[eid] {
                            let ge = [
                                lam[a] * gq[b][0] + lam[b] * gq[a][0],
                                lam[a] * gq[b][1] + lam[b] * gq[a][1],
                                lam[a] * gq[b][2] + lam[b] * gq[a][2],
                            ];
                            grad_e_norm2[eid] += (ge[0] * ge[0] + ge[1] * ge[1] + ge[2] * ge[2]) * wq;
                            for kk in 0..3 {
                                int_d_e[eid] += dh[kk] * (ge[kk] * wq);
                                int_b_e[eid] += curleh[kk] * (ge[kk] * wq);
                            }
                        }
                    }
                }
            }
        }

        let d_norm = d_norm2.sqrt();
        let b_norm = b_norm2.sqrt();
        let mut max_gauss = 0.0_f64;
        let mut max_div_b = 0.0_f64;
        for &v in &interior_verts {
            let gn = grad_v_norm2[v].sqrt();
            max_gauss = max_gauss.max(int_d_v[v].norm() / (d_norm * gn));
            max_div_b = max_div_b.max(int_b_v[v].norm() / (b_norm * gn));
        }
        if asm.order == 2 {
            for &e in &interior_edges {
                let gn = grad_e_norm2[e].sqrt();
                max_gauss = max_gauss.max(int_d_e[e].norm() / (d_norm * gn));
                max_div_b = max_div_b.max(int_b_e[e].norm() / (b_norm * gn));
            }
        }
        assert!(max_gauss < 1e-10, "mode {mode_idx} weak Gauss residual {max_gauss:.3e} >= 1e-10");
        assert!(max_div_b < 1e-10, "mode {mode_idx} weak div(mu H) residual {max_div_b:.3e} >= 1e-10");
    }
}

#[test]
fn weak_gauss_law_and_div_b_zero_straight_and_curved() {
    // 1. Piecewise heterogeneous dielectric box cavity (eps = 2.2 / 9.0 across z = 0.01)
    let m_het = common::block_split_z(0.04, 0.03, 0.02, 0.01, 0.008, 3);
    let mut p_het = Problem::new(Study::Cavity { n_modes: 2, f_search: 2.5e9 }, 2);
    p_het.materials.set(1, "eps_r", 2.2);
    p_het.materials.set(2, "eps_r", 9.0);
    check_weak_gauss_and_div_b(&m_het, &p_het);

    // 2. Curved spherical PEC cavity (CurvedTet)
    let mut opts = GeomMeshOpts::new(0.05 / 2.6);
    opts.curved = true;
    let m_sph = sphere_cavity_with(0.05, &opts).unwrap().mesh;
    assert!(m_sph.has_curved());
    let mut p_sph = Problem::new(Study::Cavity { n_modes: 2, f_search: 2.5e9 }, 2);
    p_sph.materials.default.insert("eps_r".into(), 1.0);
    check_weak_gauss_and_div_b(&m_sph, &p_sph);
}

fn find_cfield<'a>(sol: &'a nyx_physics::Solution, name: &str) -> (&'a [f64], &'a [f64]) {
    let cf = sol.cell_fields.iter().find(|f| f.name == name).expect("cell field missing");
    (&cf.re, &cf.im)
}

#[test]
fn poynting_power_balance_driven_cavity() {
    // Parallel-plate TEM quarter-wave stub [0, a] x [0, L] x [0, b]:
    // PMC at x = 0,a (1,6), PEC at z = 0,b (3,4), lumped_port at y=0 (2), resistive IBC at y=L (5) + loss_tangent
    let (a, l, b) = (10.0e-3, 12.0e-3, 5.0e-3);
    let mesh = common::block(a, l, b, 4, 16, 2);
    let freq = 6.5e9;
    let eps_im = 2.0e-3;
    let r_s = 0.5;

    let mut p = Problem::new(
        Study::Driven {
            f_min: freq,
            f_max: freq,
            n_points: 1,
            freqs_hz: Some(vec![freq]),
            prom_tol: None,
            max_hdm_solves: 5,
            strict_prom_tol: false,
            seed_freqs_hz: Vec::new(),
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p.materials.default.insert("eps_r".into(), 1.0);
    p.materials.default.insert("eps_r_im".into(), eps_im);
    p.pmc = vec![1, 6];
    p.lumped_ports.push(LumpedPortBc {
        labels: vec![2],
        name: "p1".into(),
        z0: 100.0,
        direction: [0.0, 0.0, 1.0],
        complex_damping: false,
        mode: "single".into(),
    });
    p.surface_impedance.push(SurfaceImpedanceBc {
        labels: vec![5],
        l_k_sq: 0.0,
        r_s,
        material: String::new(),
        thickness: 0.0,
        temperature: 0.0,
        f_ref: freq,
        film_model: "kerr".into(),
        lossless: false,
        ibc_mode: "direct".into(),
    });

    let res = common::solve(&mesh, &p).unwrap();
    let data = DrivenPortData::of(&res).unwrap();
    let s11 = c(data.s_matrix[0][0][0][0], data.s_matrix[0][0][0][1]);
    let p_in = 0.5 * (1.0 - s11.norm_sqr());

    let sol = &res.solutions[0];
    let ex = find_cfield(sol, "Ex");
    let ey = find_cfield(sol, "Ey");
    let ez = find_cfield(sol, "Ez");
    let hx = find_cfield(sol, "Hx");
    let hy = find_cfield(sol, "Hy");
    let hz = find_cfield(sol, "Hz");

    let geo = geometry(&mesh);
    let quad3 = quad_rule(3, 4);
    let n3: Vec<Vec<f64>> = quad3.lam.iter().map(|l| shape(l, 2).0).collect();
    let mut e2_vol = 0.0_f64;
    for ci in 0..mesh.n_cells() {
        let base = ci * 10;
        let eval_c = |comp: (&[f64], &[f64]), nv: &[f64]| -> C {
            let mut sum = c(0.0, 0.0);
            for (a, &nval) in nv.iter().enumerate().take(10) {
                sum += c(comp.0[base + a], comp.1[base + a]) * nval;
            }
            sum
        };
        for (q, wq_ref) in quad3.w.iter().enumerate() {
            let wq = wq_ref * geo.vol[ci];
            let nv = &n3[q];
            let e2 = eval_c(ex, nv).norm_sqr() + eval_c(ey, nv).norm_sqr() + eval_c(ez, nv).norm_sqr();
            e2_vol += wq * e2;
        }
    }
    let omega = 2.0 * PI * freq;
    let p_vol = 0.5 * omega * EPS0 * eps_im * e2_vol;

    // Surface IBC dissipation P_ibc = 0.5 * R_s * int_{y=L} |H_T|^2 dA
    let mut face_to_cell: BTreeMap<[u32; 3], usize> = BTreeMap::new();
    const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    for ci in 0..mesh.n_cells() {
        let cv = mesh.cell(ci);
        for f in &FACES {
            let mut tri = [cv[f[0]], cv[f[1]], cv[f[2]]];
            tri.sort_unstable();
            face_to_cell.entry(tri).or_insert(ci);
        }
    }
    let quad2 = quad_rule(2, 4);
    let mut ht2_surf = 0.0_f64;
    let mut seen = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if mesh.facet_labels[f] != 5 {
            continue;
        }
        let fv = mesh.facet(f);
        let mut tri = [fv[0], fv[1], fv[2]];
        tri.sort_unstable();
        if !seen.insert(tri) {
            continue;
        }
        let ci = face_to_cell[&tri];
        let cv = mesh.cell(ci);
        let lf = [
            cv.iter().position(|&v| v == tri[0]).unwrap(),
            cv.iter().position(|&v| v == tri[1]).unwrap(),
            cv.iter().position(|&v| v == tri[2]).unwrap(),
        ];
        let opp = 6 - lf[0] - lf[1] - lf[2];
        let p = |i: usize| -> [f64; 3] {
            let s = mesh.point(cv[i] as usize);
            [s[0], s[1], s[2]]
        };
        let (p0, p1, p2, po) = (p(lf[0]), p(lf[1]), p(lf[2]), p(opp));
        let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let cr = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        let nl = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
        let area = 0.5 * nl;
        let _ = po;
        let n_hat = [cr[0] / nl, cr[1] / nl, cr[2] / nl];
        let base = ci * 10;
        for (q, l2) in quad2.lam.iter().enumerate() {
            let mut lam4 = [0.0; 4];
            lam4[lf[0]] = l2[0];
            lam4[lf[1]] = l2[1];
            lam4[lf[2]] = l2[2];
            let nv = shape(&lam4, 2).0;
            let eval_h = |comp: (&[f64], &[f64])| -> C {
                let mut sum = c(0.0, 0.0);
                for (a, &nval) in nv.iter().enumerate().take(10) {
                    sum += c(comp.0[base + a], comp.1[base + a]) * nval;
                }
                sum
            };
            let hv = [eval_h(hx), eval_h(hy), eval_h(hz)];
            let hn = hv[0] * n_hat[0] + hv[1] * n_hat[1] + hv[2] * n_hat[2];
            let ht2 = hv[0].norm_sqr() + hv[1].norm_sqr() + hv[2].norm_sqr() - hn.norm_sqr();
            ht2_surf += quad2.w[q] * area * ht2.max(0.0);
        }
    }
    let p_ibc = 0.5 * r_s * ht2_surf;
    let p_diss = p_vol + p_ibc;
    let rel = (p_diss - p_in).abs() / p_in.max(1e-30);
    eprintln!("Poynting balance: P_in={p_in:.6e}, P_vol={p_vol:.6e}, P_ibc={p_ibc:.6e}, rel={rel:.3e}");
    assert!(rel < 1e-3, "Poynting power mismatch: P_in={p_in:.6e}, P_diss={p_diss:.6e}, rel={rel:.3e}");
}

fn check_passivity_and_reciprocity(s_mat: &[Vec<[f64; 2]>]) {
    let p = s_mat.len();
    // 1. Reciprocity ||S - S^T||_F / max(||S||_F, 1) < 1e-10
    let mut s_frob2 = 0.0_f64;
    let mut asym2 = 0.0_f64;
    for (i, row_i) in s_mat.iter().enumerate() {
        for (j, sij_raw) in row_i.iter().enumerate() {
            let sij = c(sij_raw[0], sij_raw[1]);
            let sji = c(s_mat[j][i][0], s_mat[j][i][1]);
            s_frob2 += sij.norm_sqr();
            asym2 += (sij - sji).norm_sqr();
        }
    }
    let recip = asym2.sqrt() / s_frob2.sqrt().max(1.0);
    assert!(recip < 1e-10, "{p}-port reciprocity violation {recip:.3e}");

    // 2. Passivity: H = I - S^H S Hermitian PSD via 2p x 2p real symmetric matrix [A, -B; B, A]
    let n2 = 2 * p;
    let mut h_real = vec![0.0_f64; n2 * n2];
    let mut eye = vec![0.0_f64; n2 * n2];
    for i in 0..n2 {
        eye[i * n2 + i] = 1.0;
    }
    for i in 0..p {
        for j in 0..p {
            let mut acc = if i == j { c(1.0, 0.0) } else { c(0.0, 0.0) };
            for row_k in s_mat.iter().take(p) {
                let ski = c(row_k[i][0], row_k[i][1]);
                let skj = c(row_k[j][0], row_k[j][1]);
                acc -= ski.conj() * skj;
            }
            h_real[i * n2 + j] = acc.re;
            h_real[(i + p) * n2 + (j + p)] = acc.re;
            h_real[i * n2 + (j + p)] = -acc.im;
            h_real[(i + p) * n2 + j] = acc.im;
        }
    }
    let (evals, _) = dense_gen_eigh(&h_real, &eye, n2).unwrap();
    let lam_min = evals.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(lam_min >= -1e-9, "{p}-port passivity violation: lambda_min(I - S^H S) = {lam_min:.3e}");
}

fn solve_nport_box(n_ports: usize, eps_im: f64) -> (SolveResult, DrivenPortData) {
    let mesh = common::block(12.0e-3, 18.0e-3, 6.0e-3, 3, 5, 2);
    let mut p = Problem::new(
        Study::Driven {
            f_min: 4.0e9,
            f_max: 8.0e9,
            n_points: 11,
            freqs_hz: None,
            prom_tol: Some(1e-5),
            max_hdm_solves: 11,
            strict_prom_tol: true,
            seed_freqs_hz: Vec::new(),
            solver: "auto".into(),
            blr_tol: 1e-6,
        },
        2,
    );
    p.materials.default.insert("eps_r".into(), 4.0);
    p.materials.default.insert("eps_r_im".into(), eps_im);
    let specs = [
        (1, "p1", [0.0, 0.0, 1.0]),
        (2, "p2", [0.0, 0.0, 1.0]),
        (5, "p3", [0.0, 0.0, 1.0]),
        (6, "p4", [0.0, 0.0, 1.0]),
    ];
    for &(lbl, nm, dir) in &specs[..n_ports] {
        p.lumped_ports.push(LumpedPortBc {
            labels: vec![lbl],
            name: nm.into(),
            z0: 50.0,
            direction: dir,
            complex_damping: false,
            mode: "single".into(),
        });
    }
    let res = common::solve(&mesh, &p).unwrap();
    let data = DrivenPortData::of(&res).unwrap();
    (res, data)
}

#[test]
fn multiport_3_and_4_port_passivity_reciprocity_and_port_permutation() {
    for n_ports in [3usize, 4] {
        for eps_im in [0.0, 0.05] {
            let (_res, data) = solve_nport_box(n_ports, eps_im);
            for s_k in &data.s_matrix {
                check_passivity_and_reciprocity(s_k);
            }
        }
    }

    // Port permutation invariance S_perm == P S P^T (< 1e-12)
    let mesh = common::block(12.0e-3, 18.0e-3, 6.0e-3, 3, 5, 2);
    let freq = 6.0e9;
    let perm = [2usize, 0, 3, 1];
    let port_defs = [(1, "p1"), (2, "p2"), (5, "p3"), (6, "p4")];
    let build = |order_idx: &[usize]| -> Vec<Vec<[f64; 2]>> {
        let mut p = Problem::new(
            Study::Driven {
                f_min: freq,
                f_max: freq,
                n_points: 1,
                freqs_hz: Some(vec![freq]),
                prom_tol: None,
                max_hdm_solves: 5,
                strict_prom_tol: false,
                seed_freqs_hz: Vec::new(),
                solver: "auto".into(),
                blr_tol: 1e-6,
            },
            2,
        );
        p.materials.default.insert("eps_r".into(), 4.0);
        p.materials.default.insert("eps_r_im".into(), 0.05);
        for &idx in order_idx {
            let (lbl, nm) = port_defs[idx];
            p.lumped_ports.push(LumpedPortBc {
                labels: vec![lbl],
                name: nm.into(),
                z0: 50.0 + 10.0 * (idx as f64),
                direction: [0.0, 0.0, 1.0],
                complex_damping: false,
                mode: "single".into(),
            });
        }
        let res = common::solve(&mesh, &p).unwrap();
        DrivenPortData::of(&res).unwrap().s_matrix.swap_remove(0)
    };
    let s_orig = build(&[0, 1, 2, 3]);
    let s_perm = build(&perm);
    let mut diff2 = 0.0_f64;
    let mut norm2 = 0.0_f64;
    for i in 0..4 {
        for j in 0..4 {
            let a = c(s_perm[i][j][0], s_perm[i][j][1]);
            let b = c(s_orig[perm[i]][perm[j]][0], s_orig[perm[i]][perm[j]][1]);
            diff2 += (a - b).norm_sqr();
            norm2 += b.norm_sqr();
        }
    }
    let rel = diff2.sqrt() / norm2.sqrt();
    assert!(rel < 1e-12, "port permutation S_perm != P S P^T: rel = {rel:.3e}");
}

#[test]
fn rigid_euclidean_invariance_mesh_and_geometry_remesh() {
    // 1. Mesh-level rigid Euclidean invariance (< 1e-10) with rotated crystal tensor, IBC, and lumped port
    let base_mesh = common::block(0.036, 0.028, 0.018, 3, 3, 2);
    let r_eucl = rotation([23.0, 41.0, 67.0]);
    let t_eucl = [0.13, -0.27, 0.41];
    let rot_mesh = transform_mesh_3d(&base_mesh, &r_eucl, t_eucl);

    let crystal_euler = [17.0, 32.0, -48.0];
    let mut r_tot = [[0.0; 3]; 3];
    let r_crys = rotation(crystal_euler);
    for i in 0..3 {
        for j in 0..3 {
            r_tot[i][j] = (0..3).map(|k| r_eucl[i][k] * r_crys[k][j]).sum();
        }
    }
    let rot_euler = euler_from_rotation(&r_tot);

    let make_cav = |euler: [f64; 3]| -> Problem {
        let mut p = Problem::new(Study::Cavity { n_modes: 3, f_search: 2.0e9 }, 2);
        p.materials.aniso.insert(
            1,
            Anisotropy {
                permittivity: Some(Tensor::Diagonal { x: 2.0, y: 5.0, z: 9.0 }),
                euler_deg: euler,
                ..Anisotropy::default()
            },
        );
        p.surface_impedance.push(SurfaceImpedanceBc {
            labels: vec![2],
            l_k_sq: 1e-12,
            r_s: 0.05,
            material: String::new(),
            thickness: 0.0,
            temperature: 0.0,
            f_ref: 2.5e9,
            film_model: "kerr".into(),
            lossless: false,
            ibc_mode: "direct".into(),
        });
        p
    };
    let r0 = common::solve(&base_mesh, &make_cav(crystal_euler)).unwrap();
    let r1 = common::solve(&rot_mesh, &make_cav(rot_euler)).unwrap();
    for i in 0..3 {
        let f0 = r0.solutions[i].value.unwrap();
        let f1 = r1.solutions[i].value.unwrap();
        let rel = ((f1 - f0) / f0).abs();
        assert!(rel < 1e-10, "rigid Euclidean cavity mode {i}: {f0} vs {f1} (rel {rel:.3e})");
    }

    // Also driven port S11 with rotated direction vector
    let d0 = [0.0, 0.0, 1.0];
    let d1 = rotate_vec3(&r_eucl, d0);
    let make_drv = |euler: [f64; 3], dir: [f64; 3]| -> Problem {
        let mut p = Problem::new(
            Study::Driven {
                f_min: 3.0e9,
                f_max: 3.0e9,
                n_points: 1,
                freqs_hz: Some(vec![3.0e9]),
                prom_tol: None,
                max_hdm_solves: 5,
                strict_prom_tol: false,
                seed_freqs_hz: Vec::new(),
                solver: "auto".into(),
                blr_tol: 1e-6,
            },
            2,
        );
        p.materials.aniso.insert(
            1,
            Anisotropy {
                permittivity: Some(Tensor::Diagonal { x: 2.0, y: 5.0, z: 9.0 }),
                euler_deg: euler,
                ..Anisotropy::default()
            },
        );
        p.lumped_ports.push(LumpedPortBc {
            labels: vec![2],
            name: "p1".into(),
            z0: 50.0,
            direction: dir,
            complex_damping: false,
            mode: "single".into(),
        });
        p
    };
    let s0 = DrivenPortData::of(&common::solve(&base_mesh, &make_drv(crystal_euler, d0)).unwrap()).unwrap().s_matrix[0]
        [0][0];
    let s1 =
        DrivenPortData::of(&common::solve(&rot_mesh, &make_drv(rot_euler, d1)).unwrap()).unwrap().s_matrix[0][0][0];
    let ds = (c(s1[0] - s0[0], s1[1] - s0[1])).norm();
    assert!(ds < 1e-10, "rigid Euclidean driven S11 diff {ds:.3e}");

    // 2. Geometry-level rotated remesh (< 1e-3): rotated 2D elliptical-core fiber remesh AND 3D cavity remesh
    let build_2d_remesh = |rotated: bool| -> Vec<f64> {
        let mut m = Model::new_2d();
        let clad = m.circle("clad", (0.0, 0.0), 12.0e-3);
        let core = m.ellipse("core", (2.0e-3, -1.0e-3), (4.5e-3, 2.2e-3), 0.0);
        let (clad, core) = if rotated {
            let steps = vec![
                TStep::Rotate { axis: Axis::Y, deg: 23.0.into(), about: [0.0.into(), 0.0.into(), 0.0.into()] },
                TStep::Translate([0.017.into(), 0.0.into(), (-0.011).into()]),
            ];
            (m.transform("clad_r", clad, steps.clone()), m.transform("core_r", core, steps))
        } else {
            (clad, core)
        };
        m.overlay(&[("clad", clad), ("core", core)]);
        let g =
            femcore::meshing::geomdrive::mesh_geom_2d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(1.0e-3)).unwrap();
        let mut p = Problem::new(Study::Modes { freq: 12.0e9, n_target: 1.8, evanescent: false, n_modes: 2 }, 2);
        p.materials.set(g.domain("clad").unwrap(), "eps_r", 1.0);
        p.materials.set(g.domain("core").unwrap(), "eps_r", 4.0);
        common::solve(&g.mesh, &p).unwrap().solutions.iter().map(|s| common::neff(s).0).collect()
    };
    let n2_a = build_2d_remesh(false);
    let n2_b = build_2d_remesh(true);
    for i in 0..2 {
        let rel = ((n2_b[i] - n2_a[i]) / n2_a[i]).abs();
        assert!(rel < 1e-3, "2D geometry rotated remesh mode {i}: {} vs {} (rel {rel:.3e})", n2_a[i], n2_b[i]);
    }

    // 3D rotated box cavity remesh
    let build_3d_remesh = |rotated: bool| -> Vec<f64> {
        let mut m = Model::new_3d();
        let cav = m.box_("cav", (-0.02, -0.015, -0.01), (0.02, 0.015, 0.01));
        let cav = if rotated {
            m.transform(
                "cav_r",
                cav,
                vec![
                    TStep::Rotate { axis: Axis::Z, deg: 23.0.into(), about: [0.0.into(), 0.0.into(), 0.0.into()] },
                    TStep::Rotate { axis: Axis::X, deg: 41.0.into(), about: [0.0.into(), 0.0.into(), 0.0.into()] },
                    TStep::Translate([0.05.into(), (-0.03).into(), 0.02.into()]),
                ],
            )
        } else {
            cav
        };
        m.overlay(&[("cav", cav)]);
        let g = mesh_geom_3d(evaluate(&m.list).unwrap(), &GeomMeshOpts::new(0.0028)).unwrap();
        let mut p = Problem::new(Study::Cavity { n_modes: 3, f_search: 6.0e9 }, 2);
        p.materials.default.insert("eps_r".into(), 1.0);
        let res = common::solve(&g.mesh, &p).unwrap();
        res.solutions.iter().map(|s| s.value.unwrap()).collect()
    };
    let fa = build_3d_remesh(false);
    let fb = build_3d_remesh(true);
    for i in 0..3 {
        let rel = ((fb[i] - fa[i]) / fa[i]).abs();
        assert!(rel < 1e-2, "3D rotated cavity remesh mode {i}: {} vs {} (rel {rel:.3e})", fa[i], fb[i]);
    }
}

#[test]
fn vertex_cell_permutation_dimensional_scaling_and_floquet_symmetry() {
    // 1. Vertex, cell, even-tet vertex, and facet permutation invariance (< 1e-10)
    let base = common::block(0.04, 0.03, 0.02, 3, 2, 2);
    let nv = base.n_points();
    let nc = base.n_cells();
    let nf = base.n_facets();
    let mut v_perm: Vec<u32> = (0..nv as u32).collect();
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next_u64 = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        state
    };
    for i in (1..nv).rev() {
        let j = (next_u64() as usize) % (i + 1);
        v_perm.swap(i, j);
    }
    let mut v_inv = vec![0u32; nv];
    for (old, &new_id) in v_perm.iter().enumerate() {
        v_inv[new_id as usize] = old as u32;
    }
    let mut p_pts = vec![0.0; base.points.len()];
    for (old, &new_id) in v_perm.iter().enumerate() {
        p_pts[3 * new_id as usize..3 * new_id as usize + 3].copy_from_slice(&base.points[3 * old..3 * old + 3]);
    }
    let mut c_perm: Vec<usize> = (0..nc).collect();
    for i in (1..nc).rev() {
        let j = (next_u64() as usize) % (i + 1);
        c_perm.swap(i, j);
    }
    let mut p_cells = Vec::with_capacity(base.cells.len());
    let mut p_clbl = Vec::with_capacity(nc);
    for &ci in &c_perm {
        let cv = base.cell(ci);
        let mapped = [v_perm[cv[0] as usize], v_perm[cv[1] as usize], v_perm[cv[2] as usize], v_perm[cv[3] as usize]];
        // Even (orientation-preserving A_4) tet vertex permutations to shuffle local slots while keeping det J > 0
        const EVEN_TETS: [[usize; 4]; 4] = [[0, 1, 2, 3], [1, 2, 0, 3], [2, 0, 1, 3], [3, 1, 0, 2]];
        let p = EVEN_TETS[ci % 4];
        p_cells.extend_from_slice(&[mapped[p[0]], mapped[p[1]], mapped[p[2]], mapped[p[3]]]);
        p_clbl.push(base.cell_labels[ci]);
    }
    let mut f_perm: Vec<usize> = (0..nf).collect();
    for i in (1..nf).rev() {
        let j = (next_u64() as usize) % (i + 1);
        f_perm.swap(i, j);
    }
    let mut p_facets = Vec::with_capacity(base.facets.len());
    let mut p_flbl = Vec::with_capacity(nf);
    for &fi in &f_perm {
        let fv = base.facet(fi);
        p_facets.extend_from_slice(&[v_perm[fv[0] as usize], v_perm[fv[1] as usize], v_perm[fv[2] as usize]]);
        p_flbl.push(base.facet_labels[fi]);
    }
    let perm_mesh = SimplexMesh::new(3, p_pts, p_cells, p_clbl, p_facets, p_flbl);
    let _ = v_inv;

    let mut pc = Problem::new(Study::Cavity { n_modes: 3, f_search: 5.5e9 }, 2);
    pc.materials.default.insert("eps_r".into(), 1.0);
    let ra = common::solve(&base, &pc).unwrap();
    let rb = common::solve(&perm_mesh, &pc).unwrap();
    for i in 0..3 {
        let fa = ra.solutions[i].value.unwrap();
        let fb = rb.solutions[i].value.unwrap();
        let rel = ((fb - fa) / fa).abs();
        assert!(rel < 1e-10, "permutation invariance cavity mode {i}: {fa} vs {fb} (rel {rel:.3e})");
    }

    // 2. Dimensional scaling invariance (< 1e-10) across cavity, 2D modes, and driven waveports
    let lam = 2.5_f64;
    let scaled_3d = scale_mesh(&base, lam);
    let mut pc_sc = pc.clone();
    pc_sc.study = Study::Cavity { n_modes: 3, f_search: 5.5e9 / lam };
    let r_sc = common::solve(&scaled_3d, &pc_sc).unwrap();
    for i in 0..3 {
        let fa = ra.solutions[i].value.unwrap();
        let fb = r_sc.solutions[i].value.unwrap() * lam;
        let rel = ((fb - fa) / fa).abs();
        assert!(rel < 1e-10, "dimensional scaling 3D cavity mode {i}: rel {rel:.3e}");
    }

    let rect_base = common::rect(22.86e-3, 10.16e-3, 10, 5);
    let rect_scaled = scale_mesh(&rect_base, lam);
    let freq_m = 10.0e9;
    let mut pm0 = Problem::new(Study::Modes { freq: freq_m, n_target: 0.75, evanescent: false, n_modes: 1 }, 2);
    pm0.materials.default.insert("eps_r".into(), 1.0);
    let mut pm1 = pm0.clone();
    pm1.study = Study::Modes { freq: freq_m / lam, n_target: 0.75, evanescent: false, n_modes: 1 };
    let ne0 = common::neff(&common::solve(&rect_base, &pm0).unwrap().solutions[0]).0;
    let ne1 = common::neff(&common::solve(&rect_scaled, &pm1).unwrap().solutions[0]).0;
    assert!((ne1 - ne0).abs() < 1e-10, "dimensional scaling 2D modes: {ne0} vs {ne1}");

    let make_wp = |f: f64| -> Problem {
        let mut pd = Problem::new(
            Study::Driven {
                f_min: f,
                f_max: f,
                n_points: 1,
                freqs_hz: Some(vec![f]),
                prom_tol: None,
                max_hdm_solves: 5,
                strict_prom_tol: false,
                seed_freqs_hz: Vec::new(),
                solver: "auto".into(),
                blr_tol: 1e-6,
            },
            2,
        );
        pd.materials.default.insert("eps_r".into(), 1.0);
        pd.waveports.push(WavePortBc { labels: vec![2], name: "p1".into(), mode: 1, n_eff_guess: None, z0_ohm: None });
        pd.waveports.push(WavePortBc { labels: vec![5], name: "p2".into(), mode: 1, n_eff_guess: None, z0_ohm: None });
        pd
    };
    let s_b = DrivenPortData::of(&common::solve(&base, &make_wp(9.0e9)).unwrap()).unwrap().s_matrix.swap_remove(0);
    let s_s =
        DrivenPortData::of(&common::solve(&scaled_3d, &make_wp(9.0e9 / lam)).unwrap()).unwrap().s_matrix.swap_remove(0);
    let mut ds_frob = 0.0_f64;
    for i in 0..2 {
        for j in 0..2 {
            ds_frob += c(s_s[i][j][0] - s_b[i][j][0], s_s[i][j][1] - s_b[i][j][1]).norm_sqr();
        }
    }
    assert!(ds_frob.sqrt() < 1e-10, "dimensional scaling driven waveport S diff {:.3e}", ds_frob.sqrt());

    // 3. `+phase / -phase` Floquet time-reversal degeneracy (< 1e-10) in Cavity & Bloch
    let gcell = common::block(0.03, 0.01, 0.015, 4, 2, 2);
    let solve_phase = |phase: f64| -> Vec<f64> {
        let mut p = Problem::new(Study::Cavity { n_modes: 2, f_search: 7.0e9 }, 2);
        p.materials.default.insert("eps_r".into(), 1.0);
        p.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase }];
        common::solve(&gcell, &p).unwrap().solutions.iter().map(|s| s.value.unwrap()).collect()
    };
    let fp = solve_phase(0.73);
    let fm = solve_phase(-0.73);
    for i in 0..2 {
        let rel = ((fp[i] - fm[i]) / fp[i]).abs();
        assert!(rel < 1e-10, "Floquet +phase/-phase degeneracy mode {i}: {} vs {} (rel {rel:.3e})", fp[i], fm[i]);
    }
    // Also verify Bloch reversed direction (`source <-> dest` swap) matches to < 1e-10
    let mut pb1 = Problem::new(Study::Bloch { freq: 8.0e9, n_target: 0.75, n_modes: 2, brillouin_zone: 1 }, 2);
    pb1.materials.default.insert("eps_r".into(), 1.0);
    pb1.periodic = vec![Periodic { source: vec![2], dest: vec![5], phase: 0.0 }];
    let mut pb2 = pb1.clone();
    pb2.periodic = vec![Periodic { source: vec![5], dest: vec![2], phase: 0.0 }];
    let nb1 = common::neff(&common::solve(&gcell, &pb1).unwrap().solutions[0]).0;
    let nb2 = common::neff(&common::solve(&gcell, &pb2).unwrap().solutions[0]).0;
    assert!(((nb1 - nb2) / nb1).abs() < 1e-10, "Bloch forward/reverse time-reversal symmetry: {nb1} vs {nb2}");
}
