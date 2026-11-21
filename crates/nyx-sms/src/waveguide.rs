// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Elastic waveguide modes on a 2D `(x, z)` cross-section with longitudinal wavevector `q` along `+y`:
//! the displacement is `u(x, z) exp(-i q y)` with three components stored in component blocks
//! `(0: u_x, 1: u_z, 2: u_y)` matching the cross-section Voigt order `(xx, zz, yy, xz, zy, xy)` of
//! `femcore::section::plane_stiffness(c, 2)`.
//!
//! The strains are
//!
//! ```text
//!   e_xx = du_x/dx         e_zz = du_z/dz         e_yy = -i q u_y
//!   g_xz = du_x/dz + du_z/dx   g_zy = du_y/dz - i q u_z   g_xy = du_y/dx - i q u_x
//! ```
//!
//! Inside an absorbing layer (`femcore::pml::cell_stretches`), the in-plane derivatives become
//! `(1 / s_x) d/dx` and `(1 / s_z) d/dz` and the area element becomes `s_x s_z dx dz`. The test
//! function carries the conjugate longitudinal phase `exp(+i q y)` while the complex coordinate
//! stretches `s_x, s_z` are not conjugated on the test side, so the layer blocks remain complex
//! symmetric.

use femcore::fem::*;
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::material::M6;
use nyx_physics::Discretisation;

/// Stiffness matrix over `3 * n_dofs(order)` unknowns (component blocks `u_x, u_z, u_y`) for a 2D
/// `(x, z)` waveguide cross-section at longitudinal wavevector `q` (rad/m). `c6` is per-cell
/// stiffness already permuted by `femcore::section::plane_stiffness(c, 2)` into
/// `(xx, zz, yy, xz, zy, xy)`. `stretch` holds `[s_x, s_z, 1]` per cell (`femcore::pml::cell_stretches`).
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
            b0.iter_mut().for_each(|x| *x = c(0.0, 0.0));
            b1.iter_mut().for_each(|x| *x = 0.0);
            for a in 0..nloc {
                let (gx, gy, nv) = (c(grad[a * 2], 0.0) / sx, c(grad[a * 2 + 1], 0.0) / sz, re.n[qp][a]);
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
                    ke[i * nd + j] += acc * re.w[qp];
                }
            }
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

/// The planar kinetic-energy weight `int rho |u_k|^2 dx dz` of every mode, split over the three
/// displacement components (`0: u_x` lateral, `1: u_z` vertical, `2: u_y` longitudinal) and over
/// classes of cells: `class[cell]` is the class of a cell or `usize::MAX` for a cell that is left
/// out (an absorbing layer). No PML stretch enters; the result is `[mode][class][component]`.
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
