// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Elastic waveguide modes on a 2D `(x, z)` cross-section with longitudinal wavevector `q` along `+y`:
//! the displacement is `u(x, z) exp(-i q y)` with three components stored in component blocks
//! `(0: u_x, 1: u_z, 2: u_y)` matching the cross-section Voigt order `(xx, zz, yy, xz, zy, xy)` of
//! `femcore::section::plane_stiffness(c, 2)`. Supports both `SimplexMesh` and native
//! `MixedDiscretisation` (`Tri` + `Quad`, straight and curved) with full `2x2` curvilinear PML Jacobians.

use crate::piezo_mat::PiezoCellMat;
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

#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
fn accumulate_waveguide_piezo_qpoint(
    nloc: usize,
    mat: &PiezoCellMat,
    beta: C,
    beta_h: C,
    w_k: C,
    w_coup: C,
    vals: &[f64],
    gx_gz: &[(C, C)],
    b0: &mut [C],
    b1: &mut [f64],
    db: &mut [C],
    ke_uu: &mut [C],
    ce_uu: Option<&mut [C]>,
    ke_up: &mut [C],
    ke_pu: &mut [C],
    ke_pp: &mut [C],
) {
    let nd = 3 * nloc;
    if mat.mech_active {
        accumulate_waveguide_qpoint(nloc, &mat.c6, beta, beta_h, w_k, vals, gx_gz, b0, b1, db, ke_uu);
        if let Some(ce) = ce_uu {
            if mat.eta.iter().flatten().any(|&v| v != 0.0) {
                for s in 0..6 {
                    for j in 0..nd {
                        let mut acc = c(0.0, 0.0);
                        for t in 0..6 {
                            acc += (b0[t * nd + j] + beta * b1[t * nd + j]) * mat.eta[s][t];
                        }
                        db[s * nd + j] = acc;
                    }
                }
                for i in 0..nd {
                    for j in 0..nd {
                        let mut acc = c(0.0, 0.0);
                        for s in 0..6 {
                            acc += (b0[s * nd + i] + beta_h * b1[s * nd + i]) * db[s * nd + j];
                        }
                        ce[i * nd + j] += acc * w_coup;
                    }
                }
            }
        }
    }
    if mat.elec_active {
        let mut bphi_t = Vec::with_capacity(nloc);
        let mut bphi_h = Vec::with_capacity(nloc);
        for a in 0..nloc {
            let (gx, gz) = gx_gz[a];
            let nv = vals[a];
            bphi_t.push([gx, gz, beta * nv]);
            bphi_h.push([gx, gz, beta_h * nv]);
        }
        if mat.mech_active && mat.e.iter().flatten().any(|&v| v != 0.0) {
            for i in 0..nd {
                let mut be_h = [c(0.0, 0.0); 3];
                for s in 0..6 {
                    let b_si = b0[s * nd + i] + beta_h * b1[s * nd + i];
                    for p in 0..3 {
                        be_h[p] += b_si * mat.e[p][s];
                    }
                }
                for a in 0..nloc {
                    let mut acc = c(0.0, 0.0);
                    for p in 0..3 {
                        acc += be_h[p] * bphi_t[a][p];
                    }
                    ke_up[i * nloc + a] += acc * w_coup;
                }
            }
            for j in 0..nd {
                let mut e_b = [c(0.0, 0.0); 3];
                for s in 0..6 {
                    let b_sj = b0[s * nd + j] + beta * b1[s * nd + j];
                    for p in 0..3 {
                        e_b[p] += mat.e[p][s] * b_sj;
                    }
                }
                for a in 0..nloc {
                    let mut acc = c(0.0, 0.0);
                    for p in 0..3 {
                        acc += bphi_h[a][p] * e_b[p];
                    }
                    ke_pu[a * nd + j] += acc * w_coup;
                }
            }
        }
        let mut eps_c = [[c(0.0, 0.0); 3]; 3];
        for p in 0..3 {
            for q in 0..3 {
                eps_c[p][q] = c(mat.eps_re[p][q], -mat.eps_im[p][q]);
            }
        }
        for b in 0..nloc {
            let mut eb = [c(0.0, 0.0); 3];
            for p in 0..3 {
                for q in 0..3 {
                    eb[p] += eps_c[p][q] * bphi_t[b][q];
                }
            }
            for a in 0..nloc {
                let mut acc = c(0.0, 0.0);
                for p in 0..3 {
                    acc += bphi_h[a][p] * eb[p];
                }
                ke_pp[a * nloc + b] += acc * w_coup;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn push_waveguide_piezo_cell(
    kc: &mut CCoo,
    mc: &mut CCoo,
    cc: &mut Option<CCoo>,
    n1: usize,
    cell_dofs: &[usize],
    mat: &PiezoCellMat,
    has_eta: bool,
    ke_uu: &[C],
    me: &[C],
    ce_uu: &[C],
    ke_up: &[C],
    ke_pu: &[C],
    ke_pp: &mut [C],
    u_dofs: &mut Vec<usize>,
    p_dofs: &mut Vec<usize>,
) {
    u_dofs.clear();
    for k in 0..3 {
        u_dofs.extend(cell_dofs.iter().map(|&d| d + k * n1));
    }
    p_dofs.clear();
    p_dofs.extend(cell_dofs.iter().map(|&d| d + 3 * n1));
    let nloc = cell_dofs.len();

    if mat.mech_active {
        kc.push_block(u_dofs, u_dofs, ke_uu);
        for k in 0..3 {
            let comp = &u_dofs[k * nloc..(k + 1) * nloc];
            mc.push_block(comp, comp, me);
        }
        if has_eta && mat.eta.iter().flatten().any(|&v| v != 0.0) {
            if let Some(c_coo) = cc.as_mut() {
                c_coo.push_block(u_dofs, u_dofs, ce_uu);
            }
        }
        if mat.elec_active && mat.e.iter().flatten().any(|&v| v != 0.0) {
            kc.push_block(u_dofs, p_dofs, ke_up);
            kc.push_block(p_dofs, u_dofs, ke_pu);
        }
    }
    if mat.elec_active {
        for x in ke_pp.iter_mut() {
            *x = -*x;
        }
        kc.push_block(p_dofs, p_dofs, ke_pp);
    }
}

/// Assembles 4-block `(u_x, u_z, u_y, phi)` coupled piezoelectric waveguide matrices
/// `(K_full, M_uu, C_eta_opt)` on a 2D [`Discretisation`] (`SimplexMesh`).
pub(crate) fn assemble_waveguide_piezo_simplex(
    mesh: &Discretisation,
    mats: &[PiezoCellMat],
    order: usize,
    q: f64,
    resolved: Option<&ResolvedPmlDomains>,
) -> Result<(CCsr, CCsr, Option<CCsr>), String> {
    if mesh.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let has_pml = resolved.is_some();
    let has_eta = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let deg = if has_pml { 2 * order + 4 } else { 2 * order + 2 };
    let re = RefElement::new(2, order, deg);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let mut kc = CCoo::new(4 * n1, 4 * n1);
    let mut mc = CCoo::new(3 * n1, 3 * n1);
    let mut cc = if has_eta { Some(CCoo::new(3 * n1, 3 * n1)) } else { None };

    let (mut scratch, mut u_dofs, mut p_dofs) = (Vec::new(), Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut gx_gz = vec![(c(0.0, 0.0), c(0.0, 0.0)); nloc];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke_uu = vec![c(0.0, 0.0); nd * nd];
    let mut ce_uu = vec![c(0.0, 0.0); nd * nd];
    let mut me = vec![c(0.0, 0.0); nloc * nloc];
    let mut ke_up = vec![c(0.0, 0.0); nd * nloc];
    let mut ke_pu = vec![c(0.0, 0.0); nloc * nd];
    let mut ke_pp = vec![c(0.0, 0.0); nloc * nloc];
    let beta = c(0.0, -q);
    let beta_h = c(0.0, q);

    for ci in 0..mesh.n_cells() {
        let mat = &mats[ci];
        if !mat.mech_active && !mat.elec_active {
            continue;
        }
        let dom = mesh.cell_labels[ci];
        let v = mesh.cell(ci);
        let p0 = mesh.point(v[0] as usize);
        let p1 = mesh.point(v[1] as usize);
        let p2 = mesh.point(v[2] as usize);
        let vol = geo.vol[ci];
        ke_uu.fill(c(0.0, 0.0));
        ce_uu.fill(c(0.0, 0.0));
        me.fill(c(0.0, 0.0));
        ke_up.fill(c(0.0, 0.0));
        ke_pu.fill(c(0.0, 0.0));
        ke_pp.fill(c(0.0, 0.0));

        for qp in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, qp, &mut grad);
            let (jac, loss) = match resolved {
                Some(res) => {
                    let lam = tri_barycentric(&re, qp, order);
                    let x = lam[0] * p0[0] + lam[1] * p1[0] + lam[2] * p2[0];
                    let z = lam[0] * p0[1] + lam[1] * p1[1] + lam[2] * p2[1];
                    res.eval_jacobian(dom, [x, 0.0, z], 2, false)
                }
                None => (femcore::pml::PmlJacobian::identity(), c(1.0, 0.0)),
            };
            for a in 0..nloc {
                let (dx, dz) = (grad[a * 2], grad[a * 2 + 1]);
                let gx = jac.j_inv[0][0] * dx + jac.j_inv[1][0] * dz;
                let gz = jac.j_inv[0][1] * dx + jac.j_inv[1][1] * dz;
                gx_gz[a] = (gx, gz);
            }
            let w_coup = jac.det * (re.w[qp] * vol);
            let w_k = w_coup * loss;
            accumulate_waveguide_piezo_qpoint(
                nloc,
                mat,
                beta,
                beta_h,
                w_k,
                w_coup,
                &re.n[qp],
                &gx_gz,
                &mut b0,
                &mut b1,
                &mut db,
                &mut ke_uu,
                if has_eta { Some(&mut ce_uu) } else { None },
                &mut ke_up,
                &mut ke_pu,
                &mut ke_pp,
            );
            if mat.mech_active && mat.rho != 0.0 {
                let w_m = jac.det * (re.w[qp] * vol * mat.rho);
                let nv = &re.n[qp];
                for a in 0..nloc {
                    let wa = w_m * nv[a];
                    for bb in 0..nloc {
                        me[a * nloc + bb] += wa * nv[bb];
                    }
                }
            }
        }
        mesh.cell_dofs(ci, order, &mut scratch);
        push_waveguide_piezo_cell(
            &mut kc,
            &mut mc,
            &mut cc,
            n1,
            &scratch,
            mat,
            has_eta,
            &ke_uu,
            &me,
            &ce_uu,
            &ke_up,
            &ke_pu,
            &mut ke_pp,
            &mut u_dofs,
            &mut p_dofs,
        );
    }
    Ok((kc.to_csr(), mc.to_csr(), cc.map(|x| x.to_csr())))
}

/// Native mixed-cell (`Tri` + `Quad`, straight and curved) 4-block `(u_x, u_z, u_y, phi)`
/// coupled piezoelectric waveguide assembly `(K_full, M_uu, C_eta_opt)`.
pub(crate) fn assemble_waveguide_piezo_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    mats: &[PiezoCellMat],
    q: f64,
    resolved: Option<&ResolvedPmlDomains>,
) -> Result<(CCsr, CCsr, Option<CCsr>), String> {
    if mixed.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let has_eta = mats.iter().any(|m| m.mech_active && m.eta.iter().flatten().any(|&v| v != 0.0));
    let beta = c(0.0, -q);
    let beta_h = c(0.0, q);
    let mut kc = CCoo::new(4 * n1, 4 * n1);
    let mut mc = CCoo::new(3 * n1, 3 * n1);
    let mut cc = if has_eta { Some(CCoo::new(3 * n1, 3 * n1)) } else { None };

    let (mut vals, mut dxi, mut sn, mut sdn) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut gx_gz, mut b0, mut b1, mut db) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut ke_uu, mut ce_uu, mut me) = (Vec::new(), Vec::new(), Vec::new());
    let (mut ke_up, mut ke_pu, mut ke_pp) = (Vec::new(), Vec::new(), Vec::new());
    let (mut u_dofs, mut p_dofs) = (Vec::new(), Vec::new());

    for ci in 0..mixed.n_cells() {
        let mat = &mats[ci];
        if !mat.mech_active && !mat.elec_active {
            continue;
        }
        let dom = mixed.cell_labels[ci];
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let in_pml = resolved.is_some_and(|r| r.contains(dom));
        let deg = if in_pml || iso.curved { 2 * order + 4 } else { 2 * order + 2 };
        let rule = quad_rule_mixed(ct, deg);
        let cdofs = dof_map.cell_dofs(ci);
        let nloc = cdofs.len();
        let nd = 3 * nloc;

        gx_gz.resize(nloc, (c(0.0, 0.0), c(0.0, 0.0)));
        b0.resize(6 * nd, c(0.0, 0.0));
        b1.resize(6 * nd, 0.0);
        db.resize(6 * nd, c(0.0, 0.0));
        ke_uu.clear();
        ke_uu.resize(nd * nd, c(0.0, 0.0));
        ce_uu.clear();
        ce_uu.resize(if has_eta { nd * nd } else { 0 }, c(0.0, 0.0));
        me.clear();
        me.resize(nloc * nloc, c(0.0, 0.0));
        ke_up.clear();
        ke_up.resize(nd * nloc, c(0.0, 0.0));
        ke_pu.clear();
        ke_pu.resize(nloc * nd, c(0.0, 0.0));
        ke_pp.clear();
        ke_pp.resize(nloc * nloc, c(0.0, 0.0));

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
            let w_coup = jac.det * w_da;
            let w_k = w_coup * loss;
            accumulate_waveguide_piezo_qpoint(
                nloc,
                mat,
                beta,
                beta_h,
                w_k,
                w_coup,
                &vals,
                &gx_gz,
                &mut b0,
                &mut b1,
                &mut db,
                &mut ke_uu,
                if has_eta { Some(&mut ce_uu) } else { None },
                &mut ke_up,
                &mut ke_pu,
                &mut ke_pp,
            );
            if mat.mech_active && mat.rho != 0.0 {
                let w_m = jac.det * (w_da * mat.rho);
                for a in 0..nloc {
                    let wa = w_m * vals[a];
                    for bb in 0..nloc {
                        me[a * nloc + bb] += wa * vals[bb];
                    }
                }
            }
        }
        push_waveguide_piezo_cell(
            &mut kc,
            &mut mc,
            &mut cc,
            n1,
            cdofs,
            mat,
            has_eta,
            &ke_uu,
            &me,
            &ce_uu,
            &ke_up,
            &ke_pu,
            &mut ke_pp,
            &mut u_dofs,
            &mut p_dofs,
        );
    }
    Ok((kc.to_csr(), mc.to_csr(), cc.map(|x| x.to_csr())))
}
