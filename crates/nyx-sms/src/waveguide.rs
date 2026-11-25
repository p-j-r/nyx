// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Elastic waveguide modes on a 2D `(x, z)` cross-section with longitudinal wavevector `q` along `+y`:
//! the displacement is `u(x, z) exp(-i q y)` with three components stored in component blocks
//! `(0: u_x, 1: u_z, 2: u_y)` matching the cross-section Voigt order `(xx, zz, yy, xz, zy, xy)` of
//! `femcore::section::plane_stiffness(c, 2)`. Supports both `SimplexMesh` and native
//! `MixedDiscretisation` (`Tri` + `Quad`, straight and curved) with full `2x2` curvilinear PML Jacobians.

use crate::pml::ResolvedPmlDomains;
use femcore::fem::*;
use femcore::fem_mixed::{lagrange_shape_into, quad_rule_mixed, IsoparametricCell, MixedLagrangeDofMap};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::material::M6;
use femcore::mesh::mixed::MixedDiscretisation;
use nyx_physics::Discretisation;

#[allow(clippy::erasing_op)]
fn accumulate_waveguide_qpoint(
    nloc: usize,
    dm: &M6,
    beta: C,
    beta_h: C,
    w_k: C,
    vals: &[f64],
    gx_gz: &[(C, C)],
    b0: &mut [C],
    b1: &mut [f64],
    db: &mut [C],
    ke: &mut [C],
) {
    let nd = 3 * nloc;
    b0[..6 * nd].iter_mut().for_each(|x| *x = c(0.0, 0.0));
    b1[..6 * nd].iter_mut().for_each(|x| *x = 0.0);
    for a in 0..nloc {
        let (gx, gy) = gx_gz[a];
        let nv = vals[a];
        // rows: 0 xx, 1 zz, 2 yy, 3 xz, 4 zy, 5 xy ; components: 0 u_x, 1 u_z, 2 u_y
        b0[0 * nd + 0 * nloc + a] = gx;
        b0[nd + nloc + a] = gy;
        b1[2 * nd + 2 * nloc + a] = nv;
        b0[3 * nd + 0 * nloc + a] = gy;
        b0[3 * nd + nloc + a] = gx;
        b0[4 * nd + 2 * nloc + a] = gy;
        b1[4 * nd + nloc + a] = nv;
        b0[5 * nd + 2 * nloc + a] = gx;
        b1[5 * nd + 0 * nloc + a] = nv;
    }
    for s in 0..6 {
        for j in 0..nd {
            let mut acc = c(0.0, 0.0);
            for t in 0..6 {
                acc += (b0[t * nd + j] + beta * b1[t * nd + j]) * dm[s][t];
            }
            db[s * nd + j] = acc;
        }
    }
    for i in 0..nd {
        for j in 0..nd {
            let mut acc = c(0.0, 0.0);
            for s in 0..6 {
                let bh = b0[s * nd + i] + beta_h * b1[s * nd + i];
                acc += bh * db[s * nd + j];
            }
            ke[i * nd + j] += acc * w_k;
        }
    }
}

/// Stiffness matrix over `3 * n_dofs(order)` unknowns (component blocks `u_x, u_z, u_y`) for a 2D
/// `(x, z)` waveguide cross-section at longitudinal wavevector `q` (rad/m).
#[allow(clippy::erasing_op)]
pub fn assemble_stiffness_waveguide(
    mesh: &Discretisation,
    c6: &[M6],
    order: usize,
    q: f64,
    stretch: Option<&[[C; 3]]>,
) -> Result<CCsr, String> {
    if mesh.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let re = RefElement::new(2, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut gx_gz = vec![(c(0.0, 0.0), c(0.0, 0.0)); nloc];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    let beta = c(0.0, -q);
    let beta_h = c(0.0, q);
    for ci in 0..mesh.n_cells() {
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let dm = &c6[ci];
        let (sx, sz) = stretch.map_or((c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1]));
        let det = sx * sz;
        for qp in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, qp, &mut grad);
            for a in 0..nloc {
                gx_gz[a] = (c(grad[a * 2], 0.0) / sx, c(grad[a * 2 + 1], 0.0) / sz);
            }
            accumulate_waveguide_qpoint(
                nloc,
                dm,
                beta,
                beta_h,
                c(re.w[qp], 0.0),
                &re.n[qp],
                &gx_gz,
                &mut b0,
                &mut b1,
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

fn tri_barycentric(re: &RefElement, q: usize, order: usize) -> [f64; 3] {
    if order == 1 {
        return [re.n[q][0], re.n[q][1], re.n[q][2]];
    }
    let nv = &re.n[q];
    [nv[0] + 0.5 * (nv[3] + nv[4]), nv[1] + 0.5 * (nv[3] + nv[5]), nv[2] + 0.5 * (nv[4] + nv[5])]
}

/// Quadrature-point curvilinear `(K, M)` assembly for 2.5D waveguide models on [`SimplexMesh`].
pub(crate) fn assemble_waveguide_pml_simplex(
    mesh: &Discretisation,
    c6: &[M6],
    rho: &[f64],
    order: usize,
    q: f64,
    resolved: &ResolvedPmlDomains,
    inactive: &[i32],
) -> Result<(CCsr, CCsr), String> {
    if mesh.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let re = RefElement::new(2, order, 2 * order + 4);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut kc = CCoo::new(n, n);
    let mut mc = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut gx_gz = vec![(c(0.0, 0.0), c(0.0, 0.0)); nloc];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    let mut me = vec![c(0.0, 0.0); nloc * nloc];
    let beta = c(0.0, -q);
    let beta_h = c(0.0, q);

    for ci in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[ci];
        if inactive.contains(&dom) {
            continue;
        }
        let v = mesh.cell(ci);
        let p0 = mesh.point(v[0] as usize);
        let p1 = mesh.point(v[1] as usize);
        let p2 = mesh.point(v[2] as usize);
        let dm = &c6[ci];
        let vol = geo.vol[ci];
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        me.iter_mut().for_each(|x| *x = c(0.0, 0.0));

        for qp in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, qp, &mut grad);
            let lam = tri_barycentric(&re, qp, order);
            let x = lam[0] * p0[0] + lam[1] * p1[0] + lam[2] * p2[0];
            let z = lam[0] * p0[1] + lam[1] * p1[1] + lam[2] * p2[1];
            let (jac, loss) = resolved.eval_jacobian(dom, [x, 0.0, z], 2, false);
            for a in 0..nloc {
                let (dx, dz) = (grad[a * 2], grad[a * 2 + 1]);
                let gx = jac.j_inv[0][0] * dx + jac.j_inv[1][0] * dz;
                let gz = jac.j_inv[0][1] * dx + jac.j_inv[1][1] * dz;
                gx_gz[a] = (gx, gz);
            }
            let w_k = jac.det * loss * (re.w[qp] * vol);
            accumulate_waveguide_qpoint(
                nloc, dm, beta, beta_h, w_k, &re.n[qp], &gx_gz, &mut b0, &mut b1, &mut db, &mut ke,
            );
            let w_m = jac.det * (re.w[qp] * vol * rho[ci]);
            let nv = &re.n[qp];
            for a in 0..nloc {
                let wa = w_m * nv[a];
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
        for k in 0..3 {
            let comp_dofs: Vec<usize> = scratch.iter().map(|&d| d + k * n1).collect();
            mc.push_block(&comp_dofs, &comp_dofs, &me);
        }
    }
    Ok((kc.to_csr(), mc.to_csr()))
}

/// Native mixed-cell (`Tri` + `Quad`, straight and curved) `(K, M)` assembly for 2.5D waveguide models.
pub(crate) fn assemble_waveguide_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    rho: &[f64],
    q: f64,
    resolved: Option<&ResolvedPmlDomains>,
    inactive: &[i32],
) -> Result<(CCsr, CCsr), String> {
    if mixed.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = 3 * n1;
    let beta = c(0.0, -q);
    let beta_h = c(0.0, q);
    let mut kc = CCoo::new(n, n);
    let mut mc = CCoo::new(n, n);

    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    let mut gx_gz = Vec::new();
    let mut b0 = Vec::new();
    let mut b1 = Vec::new();
    let mut db = Vec::new();
    let mut ke = Vec::new();
    let mut me = Vec::new();
    let mut full = Vec::new();

    for ci in 0..mixed.n_cells() {
        let dom = mixed.cell_labels[ci];
        if inactive.contains(&dom) {
            continue;
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

        gx_gz.resize(nloc, (c(0.0, 0.0), c(0.0, 0.0)));
        b0.resize(6 * nd, c(0.0, 0.0));
        b1.resize(6 * nd, 0.0);
        db.resize(6 * nd, c(0.0, 0.0));
        ke.clear();
        ke.resize(nd * nd, c(0.0, 0.0));
        me.clear();
        me.resize(nloc * nloc, c(0.0, 0.0));

        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let (jac, loss) = match resolved {
                Some(res) => res.eval_jacobian(dom, [jg.x[0], 0.0, jg.x[1]], 2, false),
                None => (femcore::pml::PmlJacobian::identity(), c(1.0, 0.0)),
            };
            for a in 0..nloc {
                let gp = jg.map_grad(dxi[a]);
                let gx = jac.j_inv[0][0] * gp[0] + jac.j_inv[1][0] * gp[1];
                let gz = jac.j_inv[0][1] * gp[0] + jac.j_inv[1][1] * gp[1];
                gx_gz[a] = (gx, gz);
            }
            let w_da = wq * jg.det.abs();
            let w_k = jac.det * loss * w_da;
            accumulate_waveguide_qpoint(nloc, dm, beta, beta_h, w_k, &vals, &gx_gz, &mut b0, &mut b1, &mut db, &mut ke);
            let w_m = jac.det * (w_da * rho[ci]);
            for a in 0..nloc {
                let wa = w_m * vals[a];
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

pub(crate) fn mode_energy_waveguide_mixed(
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
            let w_m = wq * jg.det.abs() * rho[ci];
            for a in 0..nloc {
                let wa = w_m * vals[a];
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

/// The planar kinetic-energy weight `int rho |u_k|^2 dx dz` of every mode, split over the three
/// displacement components (`0: u_x` lateral, `1: u_z` vertical, `2: u_y` longitudinal) and over
/// classes of cells.
pub fn mode_energy_waveguide(
    mesh: &Discretisation,
    rho: &[f64],
    order: usize,
    class: &[usize],
    n_classes: usize,
    modes: &[Vec<C>],
) -> Vec<Vec<[f64; 3]>> {
    let re = RefElement::new(2, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n1 = mesh.n_dofs(order);
    let mut out = vec![vec![[0.0f64; 3]; n_classes]; modes.len()];
    let mut dofs = Vec::new();
    let mut me = vec![0.0f64; nloc * nloc];
    for qp in 0..re.w.len() {
        for a in 0..nloc {
            for b in 0..nloc {
                me[a * nloc + b] += re.w[qp] * re.n[qp][a] * re.n[qp][b];
            }
        }
    }
    for ci in 0..mesh.n_cells() {
        let cl = class[ci];
        if cl == usize::MAX || rho[ci] == 0.0 {
            continue;
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
