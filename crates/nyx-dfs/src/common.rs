// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Helpers shared by the electrostatics and heat kernels: per-cell gradients and values (not averaged over the cells at a node), the
//! component and norm arrays of a vector field, the minimum and maximum of an array.

use femcore::fields::gradients_at_nodes;
use femcore::material::M3;
use nyx_physics::Discretisation;

/// The wafer tensors of the cells in the axes of the mesh (on a cross-section the x-z block comes first, see `section::plane_tensor`).
pub(crate) fn plane_tensors(mesh: &Discretisation, t: Vec<M3>) -> Vec<M3> {
    t.iter().map(|m| femcore::section::plane_tensor(m, mesh.dim)).collect()
}

/// `scale[c] * grad u` at the local dofs of every cell, not averaged: `g[(c*nloc + a)*dim + k]`.
pub(crate) fn cell_gradient(mesh: &Discretisation, u: &[f64], order: usize, scale: &[f64]) -> Vec<f64> {
    let d = mesh.dim;
    let nloc = mesh.nloc(order);
    let mut g = gradients_at_nodes(mesh, u, 1, order);
    for c in 0..mesh.n_cells() {
        for a in 0..nloc {
            for k in 0..d {
                g[(c * nloc + a) * d + k] *= scale[c];
            }
        }
    }
    g
}

/// `sign * K grad(u)` at the local dofs of every cell, not averaged, `dim` values per local dof.
pub(crate) fn cell_tensor_gradient(mesh: &Discretisation, u: &[f64], order: usize, k: &[M3], sign: f64) -> Vec<f64> {
    let d = mesh.dim;
    let nloc = mesh.nloc(order);
    let mut g = gradients_at_nodes(mesh, u, 1, order);
    for c in 0..mesh.n_cells() {
        for a in 0..nloc {
            let base = (c * nloc + a) * d;
            let gv: Vec<f64> = g[base..base + d].to_vec();
            for i in 0..d {
                g[base + i] = sign * (0..d).map(|j| k[c][i][j] * gv[j]).sum::<f64>();
            }
        }
    }
    g
}

pub(crate) fn components(g: &[f64], d: usize) -> Vec<Vec<f64>> {
    (0..d).map(|k| g.chunks(d).map(|p| p[k]).collect()).collect()
}

pub(crate) fn norm_of(comps: &[Vec<f64>]) -> Vec<f64> {
    (0..comps[0].len()).map(|i| comps.iter().map(|c| c[i] * c[i]).sum::<f64>().sqrt()).collect()
}

pub(crate) fn minmax(v: &[f64]) -> (f64, f64) {
    v.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| (a.min(x), b.max(x)))
}
