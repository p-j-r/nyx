// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Elastic modes of a body of revolution (a ring, a disk, a wheel, a hemisphere) on its half-plane
//! cross-section: the displacement is `u(r, z) exp(-i m theta)` with azimuthal order `m`, and the
//! mesh has `x = r`, `y = z`. Supports both `SimplexMesh` and native `MixedDiscretisation` (`Tri` +
//! `Quad`, straight and curved), with full `2x2` `(r, z)` curvilinear Jacobians and `s_theta = r_tilde / r`
//! for Cylindrical and Spherical PMLs (including shells touching `r = 0`).

use crate::pml::ResolvedPmlDomains;
use femcore::fem::*;
use femcore::fem_mixed::{lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::material::M6;
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::pml::{cell_stretches, radial_coordinate};
use femcore::props::{Pml, PmlDir};
use nyx_physics::Discretisation;

/// The stiffness (solver Voigt order xx, yy, zz, xy, yz, xz in the frame of the cross-section)
/// averaged over rotations of the crystal about the axis `axis` (0, 1 or 2): what a ring sees when the
/// crystal axes stay fixed while the local frame turns around the ring.
pub fn average_about_axis(c6: &M6, axis: usize, samples: usize) -> M6 {
    let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            for (p, q) in [(a, b), (b, a)] {
                for (s, u) in [(cc, d), (d, cc)] {
                    t[p][q][s][u] = c6[i][j];
                }
            }
        }
    }
    let mut avg = [[[[0.0; 3]; 3]; 3]; 3];
    for sidx in 0..samples {
        let th = 2.0 * std::f64::consts::PI * sidx as f64 / samples as f64;
        let (sn, cs) = th.sin_cos();
        let mut r = [[0.0; 3]; 3];
        let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
        r[axis][axis] = 1.0;
        r[a1][a1] = cs;
        r[a1][a2] = -sn;
        r[a2][a1] = sn;
        r[a2][a2] = cs;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                for cc in 0..3 {
                                    for d in 0..3 {
                                        s += r[i][a] * r[j][b] * r[k][cc] * r[l][d] * t[a][b][cc][d];
                                    }
                                }
                            }
                        }
                        avg[i][j][k][l] += s / samples as f64;
                    }
                }
            }
        }
    }
    let mut out = [[0.0; 6]; 6];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            out[i][j] = avg[a][b][cc][d];
        }
    }
    out
}

/// Per-cell stretch `[s_r, s_z, s_theta]` of the absorbing layers of an axisymmetric model (mesh x = r, mesh y = z).
/// Retained for legacy per-cell callers (`tests/axisym_pml.rs`).
pub fn cylindrical_stretches(mesh: &Discretisation, pml: &std::collections::BTreeMap<i32, Pml>) -> Vec<[C; 3]> {
    let mut out = cell_stretches(mesh, pml);
    let mut boxes: std::collections::BTreeMap<i32, (f64, f64)> = Default::default();
    for ci in 0..mesh.n_cells() {
        let l = mesh.cell_labels[ci];
        if pml.get(&l).is_some_and(|p| p.x != PmlDir::Off) {
            let e = boxes.entry(l).or_insert((f64::INFINITY, f64::NEG_INFINITY));
            for &v in mesh.cell(ci) {
                let r = mesh.point(v as usize)[0];
                e.0 = e.0.min(r);
                e.1 = e.1.max(r);
            }
        }
    }
    for (ci, s) in out.iter_mut().enumerate() {
        let l = mesh.cell_labels[ci];
        let (Some(p), Some(&(lo, hi))) = (pml.get(&l), boxes.get(&l)) else { continue };
        let r = (0..3).map(|k| mesh.point(mesh.cell(ci)[k] as usize)[0]).sum::<f64>() / 3.0;
        s[2] = radial_coordinate(p, r, lo, hi) / r;
    }
    out
}

#[allow(clippy::erasing_op)]
fn accumulate_axisym_qpoint(
    nloc: usize,
    dm: &M6,
    m: f64,
    r: f64,
    ist: C,
    w_k: C,
    vals: &[f64],
    gr_gz: &[(C, C)],
    b0: &mut [C],
    b1: &mut [f64],
    b2: &mut [f64],
    db: &mut [C],
    ke: &mut [C],
) {
    let nd = 3 * nloc;
    let beta = c(0.0, -m / r) * ist;
    let beta_h = c(0.0, m / r) * ist;
    let inv_r = ist / r;
    b0[..6 * nd].iter_mut().for_each(|x| *x = c(0.0, 0.0));
    b1[..6 * nd].iter_mut().for_each(|x| *x = 0.0);
    b2[..6 * nd].iter_mut().for_each(|x| *x = 0.0);
    for a in 0..nloc {
        let (gx, gy) = gr_gz[a];
        let nv = vals[a];
        // rows: 0 rr, 1 zz, 2 tt, 3 rz, 4 zt, 5 rt ; components 0 u_r, 1 u_z, 2 u_t
        b0[0 * nd + 0 * nloc + a] = gx;
        b0[nd + nloc + a] = gy;
        b2[2 * nd + 0 * nloc + a] = nv;
        b1[2 * nd + 2 * nloc + a] = nv;
        b0[3 * nd + 0 * nloc + a] = gy;
        b0[3 * nd + nloc + a] = gx;
        b0[4 * nd + 2 * nloc + a] = gy;
        b1[4 * nd + nloc + a] = nv;
        b0[5 * nd + 2 * nloc + a] = gx;
        b2[5 * nd + 2 * nloc + a] = -nv;
        b1[5 * nd + 0 * nloc + a] = nv;
    }
    for s in 0..6 {
        for j in 0..nd {
            let mut acc = c(0.0, 0.0);
            for t in 0..6 {
                acc += (b0[t * nd + j] + beta * b1[t * nd + j] + inv_r * b2[t * nd + j]) * dm[s][t];
            }
            db[s * nd + j] = acc;
        }
    }
    let w_kr = w_k * r;
    for i in 0..nd {
        for j in 0..nd {
            let mut acc = c(0.0, 0.0);
            for s in 0..6 {
                let bh = b0[s * nd + i] + beta_h * b1[s * nd + i] + inv_r * b2[s * nd + i];
                acc += bh * db[s * nd + j];
            }
            ke[i * nd + j] += acc * w_kr;
        }
    }
}

/// Stiffness matrix over `3 * n_dofs(order)` unknowns (component blocks u_r, u_z, u_theta). `stretch` holds `[s_r, s_z, s_theta]` per cell
/// (`cylindrical_stretches`).
#[allow(clippy::erasing_op)]
pub fn assemble_stiffness_axisym(
    mesh: &Discretisation,
    c6: &[M6],
    order: usize,
    m: f64,
    stretch: Option<&[[C; 3]]>,
) -> Result<CCsr, String> {
    let re = RefElement::new(2, order, 2 * order + 2);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let tol = axis_tolerance(mesh);
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut gr_gz = vec![(c(0.0, 0.0), c(0.0, 0.0)); nloc];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut b2 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    for ci in 0..mesh.n_cells() {
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let dm = &c6[ci];
        let (sx, sy, st) = stretch.map_or((c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1], s[ci][2]));
        let ist = c(1.0, 0.0) / st;
        let det = sx * sy * st;
        let v = mesh.cell(ci);
        let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
        if xs.iter().any(|&x| x < -tol) {
            return Err(format!("the mesh of an axisymmetric problem must lie in r >= 0 (cell {ci} crosses the axis)"));
        }
        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            if r <= 0.0 {
                return Err(format!(
                    "cell {ci} of an axisymmetric problem has no extent in r (all vertices on the axis)"
                ));
            }
            for a in 0..nloc {
                gr_gz[a] = (c(grad[a * 2], 0.0) / sx, c(grad[a * 2 + 1], 0.0) / sy);
            }
            accumulate_axisym_qpoint(
                nloc,
                dm,
                m,
                r,
                ist,
                c(re.w[q], 0.0),
                &re.n[q],
                &gr_gz,
                &mut b0,
                &mut b1,
                &mut b2,
                &mut db,
                &mut ke,
            );
        }
        let vol = c(geo.vol[ci], 0.0) * det;
        ke.iter_mut().for_each(|x| *x *= vol);
        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        coo.push_block(&full, &full, &ke);
    }
    Ok(coo.to_csr())
}

/// Quadrature-point curvilinear `(K, M)` assembly for 2D axisymmetric models on [`SimplexMesh`]
/// (supports Cartesian, Cylindrical, and Spherical PMLs touching `r = 0` with full `2x2` `(r, z)` coupling).
pub(crate) fn assemble_axisym_pml_simplex(
    mesh: &Discretisation,
    c6: &[M6],
    rho: &[f64],
    order: usize,
    m: f64,
    resolved: &ResolvedPmlDomains,
    inactive: &[i32],
) -> Result<(CCsr, CCsr), String> {
    let re = RefElement::new(2, order, 2 * order + 4);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let tol = axis_tolerance(mesh);
    let mut kc = CCoo::new(n, n);
    let mut mc = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut gr_gz = vec![(c(0.0, 0.0), c(0.0, 0.0)); nloc];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut b2 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    let mut me = vec![c(0.0, 0.0); nloc * nloc];
    let mut mblock = vec![c(0.0, 0.0); 9 * nloc * nloc];

    for ci in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[ci];
        if inactive.contains(&dom) {
            continue;
        }
        let v = mesh.cell(ci);
        let p0 = mesh.point(v[0] as usize);
        let p1 = mesh.point(v[1] as usize);
        let p2 = mesh.point(v[2] as usize);
        let xs = [p0[0], p1[0], p2[0]];
        let zs = [p0[1], p1[1], p2[1]];
        if xs.iter().any(|&x| x < -tol) {
            return Err(format!("the mesh of an axisymmetric problem must lie in r >= 0 (cell {ci} crosses the axis)"));
        }

        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        me.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let dm = &c6[ci];
        let vol = geo.vol[ci];

        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            let z = lam[0] * zs[0] + lam[1] * zs[1] + lam[2] * zs[2];
            if r <= 0.0 {
                return Err(format!(
                    "cell {ci} of an axisymmetric problem has no extent in r (all vertices on the axis)"
                ));
            }
            let (jac, loss) = resolved.eval_jacobian(dom, [r, 0.0, z], 2, true);
            let ist = c(1.0, 0.0) / jac.s_theta;
            for a in 0..nloc {
                let (dr, dz) = (grad[a * 2], grad[a * 2 + 1]);
                let gr = jac.j_inv[0][0] * dr + jac.j_inv[1][0] * dz;
                let gz = jac.j_inv[0][1] * dr + jac.j_inv[1][1] * dz;
                gr_gz[a] = (gr, gz);
            }
            let w_k = jac.det * loss * (re.w[q] * vol);
            accumulate_axisym_qpoint(
                nloc, dm, m, r, ist, w_k, &re.n[q], &gr_gz, &mut b0, &mut b1, &mut b2, &mut db, &mut ke,
            );
            let w_mr = jac.det * (re.w[q] * vol * rho[ci] * r);
            let nv = &re.n[q];
            for a in 0..nloc {
                let wa = w_mr * nv[a];
                for bb in 0..nloc {
                    me[a * nloc + bb] += wa * nv[bb];
                }
            }
        }

        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        kc.push_block(&full, &full, &ke);

        mblock.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        for k in 0..3 {
            for a in 0..nloc {
                for bb in 0..nloc {
                    mblock[(k * nloc + a) * nd + k * nloc + bb] = me[a * nloc + bb];
                }
            }
        }
        mc.push_block(&full, &full, &mblock[..nd * nd]);
    }
    Ok((kc.to_csr(), mc.to_csr()))
}

/// Native mixed-cell (`Tri` + `Quad`, straight and curved) `(K, M)` assembly for 2D axisymmetric models.
pub(crate) fn assemble_axisym_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    rho: &[f64],
    m: f64,
    resolved: Option<&ResolvedPmlDomains>,
    inactive: &[i32],
) -> Result<(CCsr, CCsr), String> {
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = 3 * n1;
    let tol = axis_tolerance_mixed(mixed);
    let mut kc = CCoo::new(n, n);
    let mut mc = CCoo::new(n, n);

    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    let mut gr_gz = Vec::new();
    let mut b0 = Vec::new();
    let mut b1 = Vec::new();
    let mut b2 = Vec::new();
    let mut db = Vec::new();
    let mut ke = Vec::new();
    let mut me = Vec::new();
    let mut full = Vec::new();

    for ci in 0..mixed.n_cells() {
        let dom = mixed.cell_labels[ci];
        if inactive.contains(&dom) {
            continue;
        }
        for &v in mixed.cell(ci) {
            if mixed.point(v as usize)[0] < -tol {
                return Err(format!(
                    "the mesh of an axisymmetric problem must lie in r >= 0 (cell {ci} crosses the axis)"
                ));
            }
        }
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let in_pml = resolved.is_some_and(|r| r.contains(dom));
        let deg = if in_pml || iso.curved { 2 * order + 4 } else { 2 * order + 2 };
        let rule = quad_rule_mixed(ct, deg);
        let cdofs = dof_map.cell_dofs(ci);
        let nloc = cdofs.len();
        let nd = 3 * nloc;
        let dm = &c6[ci];

        gr_gz.resize(nloc, (c(0.0, 0.0), c(0.0, 0.0)));
        b0.resize(6 * nd, c(0.0, 0.0));
        b1.resize(6 * nd, 0.0);
        b2.resize(6 * nd, 0.0);
        db.resize(6 * nd, c(0.0, 0.0));
        ke.clear();
        ke.resize(nd * nd, c(0.0, 0.0));
        me.clear();
        me.resize(nloc * nloc, c(0.0, 0.0));

        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            let r = jg.x[0];
            let z = jg.x[1];
            if r <= 0.0 {
                return Err(format!(
                    "cell {ci} of an axisymmetric problem has non-positive radius r = {r:.3e} at a quadrature point"
                ));
            }
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let (jac, loss) = match resolved {
                Some(res) => res.eval_jacobian(dom, [r, 0.0, z], 2, true),
                None => (femcore::pml::PmlJacobian::identity(), c(1.0, 0.0)),
            };
            let ist = c(1.0, 0.0) / jac.s_theta;
            for a in 0..nloc {
                let gp = jg.map_grad(dxi[a]);
                let gr = jac.j_inv[0][0] * gp[0] + jac.j_inv[1][0] * gp[1];
                let gz = jac.j_inv[0][1] * gp[0] + jac.j_inv[1][1] * gp[1];
                gr_gz[a] = (gr, gz);
            }
            let w_da = wq * jg.det.abs();
            let w_k = jac.det * loss * w_da;
            accumulate_axisym_qpoint(
                nloc, dm, m, r, ist, w_k, &vals, &gr_gz, &mut b0, &mut b1, &mut b2, &mut db, &mut ke,
            );
            let w_mr = jac.det * (w_da * rho[ci] * r);
            for a in 0..nloc {
                let wa = w_mr * vals[a];
                for bb in 0..nloc {
                    me[a * nloc + bb] += wa * vals[bb];
                }
            }
        }

        full.clear();
        for k in 0..3 {
            full.extend(cdofs.iter().map(|&gd| gd + k * n1));
        }
        kc.push_block(&full, &full, &ke);
        for k in 0..3 {
            let comp_dofs: Vec<usize> = cdofs.iter().map(|&gd| gd + k * n1).collect();
            mc.push_block(&comp_dofs, &comp_dofs, &me);
        }
    }
    Ok((kc.to_csr(), mc.to_csr()))
}

pub(crate) fn assemble_mass_axisym_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    rho: &[f64],
) -> CCsr {
    let space = dof_map.space;
    let deg_base = 2 * space.order() + 2;
    let n1 = dof_map.n_dofs;
    let n = 3 * n1;
    let mut mc = CCoo::new(n, n);
    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    let mut me = Vec::new();

    for ci in 0..mixed.n_cells() {
        if rho[ci] == 0.0 {
            continue;
        }
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let rule = quad_rule_mixed(ct, if iso.curved { deg_base.max(6) } else { deg_base });
        let cdofs = dof_map.cell_dofs(ci);
        let nloc = cdofs.len();
        me.clear();
        me.resize(nloc * nloc, c(0.0, 0.0));
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let w_mr = wq * jg.det.abs() * rho[ci] * jg.x[0].max(0.0);
            for a in 0..nloc {
                let wa = w_mr * vals[a];
                for bb in 0..nloc {
                    me[a * nloc + bb] += c(wa * vals[bb], 0.0);
                }
            }
        }
        for k in 0..3 {
            let comp_dofs: Vec<usize> = cdofs.iter().map(|&gd| gd + k * n1).collect();
            mc.push_block(&comp_dofs, &comp_dofs, &me);
        }
    }
    mc.to_csr()
}

pub(crate) fn mode_energy_axisym_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    rho: &[f64],
    class: &[usize],
    n_classes: usize,
    modes: &[Vec<C>],
) -> Vec<Vec<[f64; 3]>> {
    let space = dof_map.space;
    let deg_base = 2 * space.order() + 2;
    let n1 = dof_map.n_dofs;
    let mut out = vec![vec![[0.0f64; 3]; n_classes]; modes.len()];
    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    let mut me = Vec::new();

    for ci in 0..mixed.n_cells() {
        let cl = class[ci];
        if cl == usize::MAX || rho[ci] == 0.0 {
            continue;
        }
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let rule = quad_rule_mixed(ct, if iso.curved { deg_base.max(6) } else { deg_base });
        let cdofs = dof_map.cell_dofs(ci);
        let nloc = cdofs.len();
        me.clear();
        me.resize(nloc * nloc, 0.0);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let w_mr = wq * jg.det.abs() * rho[ci] * jg.x[0].max(0.0);
            for a in 0..nloc {
                let wa = w_mr * vals[a];
                for bb in 0..nloc {
                    me[a * nloc + bb] += wa * vals[bb];
                }
            }
        }
        for (j, u) in modes.iter().enumerate() {
            for k in 0..3 {
                let mut acc = 0.0;
                for a in 0..nloc {
                    let ua = u[k * n1 + cdofs[a]];
                    for bb in 0..nloc {
                        acc += me[a * nloc + bb] * (ua.conj() * u[k * n1 + cdofs[bb]]).re;
                    }
                }
                out[j][cl][k] += acc;
            }
        }
    }
    out
}

/// Barycentric coordinates of quadrature point `q` of the reference element rule.
fn barycentric(re: &RefElement, q: usize, order: usize) -> [f64; 3] {
    if order == 1 {
        return [re.n[q][0], re.n[q][1], re.n[q][2]];
    }
    let nv = &re.n[q];
    [nv[0] + (nv[3] + nv[4]) / 2.0, nv[1] + (nv[3] + nv[5]) / 2.0, nv[2] + (nv[4] + nv[5]) / 2.0]
}

/// Mass matrix `rho r` over `3 * n_dofs` unknowns, with the PML stretch determinant `s_r s_z s_theta` as a complex weight.
pub fn assemble_mass_axisym(mesh: &Discretisation, rho: &[f64], order: usize, stretch: Option<&[[C; 3]]>) -> CCsr {
    let re = RefElement::new(2, order, 2 * order + 1);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut me = vec![c(0.0, 0.0); nloc * nloc];
    let mut block = vec![c(0.0, 0.0); 9 * nloc * nloc];
    for ci in 0..mesh.n_cells() {
        me.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let v = mesh.cell(ci);
        let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
        for q in 0..re.w.len() {
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            for a in 0..nloc {
                for b in 0..nloc {
                    me[a * nloc + b] += c(re.w[q] * re.n[q][a] * re.n[q][b] * r, 0.0);
                }
            }
        }
        let det = stretch.map_or(c(1.0, 0.0), |s| s[ci][0] * s[ci][1] * s[ci][2]);
        let sc = det * (rho[ci] * geo.vol[ci]);
        block.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let nd = 3 * nloc;
        for k in 0..3 {
            for a in 0..nloc {
                for b in 0..nloc {
                    block[(k * nloc + a) * nd + k * nloc + b] = me[a * nloc + b] * sc;
                }
            }
        }
        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        coo.push_block(&full, &full, &block[..nd * nd]);
    }
    coo.to_csr()
}

/// The kinetic-energy weight `int rho |u_k|^2 r dr dz` of every mode, split over the three displacement components (radial, axial,
/// azimuthal) and over classes of cells.
pub fn mode_energy_axisym(
    mesh: &Discretisation,
    rho: &[f64],
    order: usize,
    class: &[usize],
    n_classes: usize,
    modes: &[Vec<C>],
) -> Vec<Vec<[f64; 3]>> {
    let re = RefElement::new(2, order, 2 * order + 1);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n1 = mesh.n_dofs(order);
    let mut out = vec![vec![[0.0f64; 3]; n_classes]; modes.len()];
    let mut dofs = Vec::new();
    let mut me = vec![0.0f64; nloc * nloc];
    for ci in 0..mesh.n_cells() {
        let cl = class[ci];
        if cl == usize::MAX || rho[ci] == 0.0 {
            continue;
        }
        me.iter_mut().for_each(|x| *x = 0.0);
        let v = mesh.cell(ci);
        let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
        for q in 0..re.w.len() {
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            for a in 0..nloc {
                for b in 0..nloc {
                    me[a * nloc + b] += re.w[q] * re.n[q][a] * re.n[q][b] * r;
                }
            }
        }
        let sc = rho[ci] * geo.vol[ci];
        mesh.cell_dofs(ci, order, &mut dofs);
        for (j, u) in modes.iter().enumerate() {
            for k in 0..3 {
                let mut acc = 0.0;
                for a in 0..nloc {
                    let ua = u[k * n1 + dofs[a]];
                    for b in 0..nloc {
                        acc += me[a * nloc + b] * (ua.conj() * u[k * n1 + dofs[b]]).re;
                    }
                }
                out[j][cl][k] += acc * sc;
            }
        }
    }
    out
}

/// Distance from the axis below which a point counts as on it: 1e-9 of the extent of the mesh.
pub fn axis_tolerance(mesh: &Discretisation) -> f64 {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..mesh.n_points() {
        let p = mesh.point(i);
        lo = lo.min(p[0].min(p[1]));
        hi = hi.max(p[0].max(p[1]));
    }
    1e-9 * (hi - lo).max(f64::MIN_POSITIVE)
}

pub(crate) fn axis_tolerance_mixed(mixed: &MixedDiscretisation) -> f64 {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..mixed.n_points() {
        let p = mixed.point(i);
        lo = lo.min(p[0].min(p[1]));
        hi = hi.max(p[0].max(p[1]));
    }
    1e-9 * (hi - lo).max(f64::MIN_POSITIVE)
}

/// The essential conditions on the axis `r = 0` (`0 u_r`, `1 u_z`, `2 u_theta`), and the number of dofs found on the axis.
pub fn axis_constraints(mesh: &Discretisation, order: usize, m: f64) -> Result<(Vec<(usize, usize)>, usize), String> {
    let tol = axis_tolerance(mesh);
    let coords = mesh.dof_coordinates(order);
    let on_axis: Vec<usize> = (0..coords.len() / 2).filter(|&d| coords[2 * d] <= tol).collect();
    if on_axis.is_empty() {
        return Ok((vec![], 0));
    }
    if (m - m.round()).abs() > 1e-9 {
        return Err(format!("the mesh touches the axis, so the azimuthal order must be an integer (got {m})"));
    }
    let blocks: &[usize] = match m.round().abs() as u64 {
        0 => &[0, 2],
        1 => &[1],
        _ => &[0, 1, 2],
    };
    let pairs = on_axis.iter().flat_map(|&d| blocks.iter().map(move |&k| (k, d))).collect();
    Ok((pairs, on_axis.len()))
}

pub(crate) fn axis_constraints_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    m: f64,
) -> Result<(Vec<(usize, usize)>, usize), String> {
    let tol = axis_tolerance_mixed(mixed);
    let coords = dof_map.dof_coordinates(mixed);
    let on_axis: Vec<usize> = (0..dof_map.n_dofs).filter(|&d| coords[2 * d] <= tol).collect();
    if on_axis.is_empty() {
        return Ok((vec![], 0));
    }
    if (m - m.round()).abs() > 1e-9 {
        return Err(format!("the mesh touches the axis, so the azimuthal order must be an integer (got {m})"));
    }
    let blocks: &[usize] = match m.round().abs() as u64 {
        0 => &[0, 2],
        1 => &[1],
        _ => &[0, 1, 2],
    };
    let pairs = on_axis.iter().flat_map(|&d| blocks.iter().map(move |&k| (k, d))).collect();
    Ok((pairs, on_axis.len()))
}

#[allow(dead_code)]
fn _uses_lagrange_space(_: LagrangeSpace) {}
