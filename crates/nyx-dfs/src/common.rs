// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Helpers shared by the electrostatics, heat, and magnetostatics kernels: per-cell gradients and
//! values, native mixed-cell (`Quad`, `Prism`, `Hex`) assembly, mixed boundary/far-field Robin
//! operators, and static Infinite Element (`Rational` and `Exponential`) coordinate transformations.

use crate::problem::{BcKind, IeCoordKind, IeProfile, InfiniteElementSpec, Problem};
use femcore::farfield::FarFieldPermittivity;
use femcore::fem::{assemble_parallel, free_dofs, geometry, quad_rule_collapsed, shape, CurvedTet, Quad};
use femcore::fem_mixed::{
    eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, lagrange_shape_into, quad_rule_mixed,
    IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::gradients_at_nodes;
use femcore::linalg::csr::{c, Coo, RCsr as Csr, C};
use femcore::linalg::driver::DirectSolver;
use femcore::material::{self, M3};
use femcore::mesh::mixed::{FacetType, MixedDiscretisation, QUAD_EDGES};
use femcore::pml::{ie_exponential_jacobian, ie_rational_jacobian, PmlJacobian};
use femcore::props::{LinearSolver, Materials, Pml, PmlCoord, PmlDir, TensorKind};
use femcore::vvuq::SourceField;
use nyx_physics::{Discretisation, Mesh, MixedToSimplexProjector};
use std::collections::{BTreeMap, BTreeSet, HashMap};

const EPS_GEOM: f64 = 1.0e-3;
const RADIAL_PROJ_REG: f64 = 0.05;

/// The wafer tensors of the cells in the axes of the mesh (on a cross-section the x-z block comes first, see `section::plane_tensor`).
pub(crate) fn plane_tensors(mesh: &Discretisation, t: Vec<M3>) -> Vec<M3> {
    t.iter().map(|m| femcore::section::plane_tensor(m, mesh.dim)).collect()
}

/// The wafer tensors of the cells of a [`MixedDiscretisation`] in the axes of the mesh.
pub(crate) fn plane_tensors_mixed(mesh: &MixedDiscretisation, t: Vec<M3>) -> Vec<M3> {
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

// ------------------------------------------------------------------ mixed-mesh material & label helpers

pub(crate) fn check_boundary_labels_mixed(
    mesh: &MixedDiscretisation,
    labels: impl IntoIterator<Item = i32>,
) -> Result<(), String> {
    let known: BTreeSet<i32> = mesh.facet_labels.iter().copied().collect();
    for l in labels {
        if !known.contains(&l) {
            return Err(format!("the mesh has no boundary with number {l}"));
        }
    }
    Ok(())
}

pub(crate) fn check_domain_labels_mixed(
    mesh: &MixedDiscretisation,
    labels: impl IntoIterator<Item = i32>,
) -> Result<(), String> {
    let known: BTreeSet<i32> = mesh.cell_labels.iter().copied().collect();
    for l in labels {
        if !known.contains(&l) {
            return Err(format!("the mesh has no domain with number {l}"));
        }
    }
    Ok(())
}

pub(crate) fn cell_values_mixed(
    materials: &Materials,
    mesh: &MixedDiscretisation,
    name: &str,
    fallback: Option<f64>,
) -> Result<Vec<f64>, String> {
    let mut by_domain: BTreeMap<i32, f64> = BTreeMap::new();
    for d in mesh.domain_labels() {
        let v = materials.get(d, name).or(fallback).ok_or_else(|| format!("domain {d} has no value for '{name}'"))?;
        by_domain.insert(d, v);
    }
    Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
}

pub(crate) fn cell_tensor_mixed(
    materials: &Materials,
    mesh: &MixedDiscretisation,
    kind: TensorKind,
    fallback: Option<f64>,
) -> Result<Vec<M3>, String> {
    let (scalar_name, pick): (&str, fn(&femcore::material::Anisotropy) -> Option<&femcore::material::Tensor>) =
        match kind {
            TensorKind::Conductivity => ("k", |a| a.conductivity.as_ref()),
            TensorKind::Permittivity => ("eps_r", |a| a.permittivity.as_ref()),
            TensorKind::PermittivityIm => ("eps_r_im", |a| a.permittivity_im.as_ref()),
            TensorKind::Expansion => ("alpha", |a| a.expansion.as_ref()),
        };
    let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
    for d in mesh.domain_labels() {
        let t = match materials.aniso.get(&d).and_then(pick) {
            Some(t) => material::tensor_global(t, materials.aniso[&d].euler_deg),
            None => {
                let v = materials
                    .get(d, scalar_name)
                    .or(fallback)
                    .ok_or_else(|| format!("domain {d} has no value for '{scalar_name}'"))?;
                material::iso_tensor(v)
            }
        };
        by_domain.insert(d, t);
    }
    Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
}

// ------------------------------------------------------------------ Infinite Element resolution & Jacobians

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[inline]
fn permute_wafer_to_2d(m: [[C; 3]; 3]) -> [[C; 3]; 3] {
    const P: [usize; 3] = [0, 2, 1];
    let mut out = [[c(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = m[P[i]][P[j]];
        }
    }
    out
}

/// Resolved per-domain Infinite Element geometry and stretching parameters.
#[derive(Clone, Debug)]
pub(crate) struct ResolvedIeDomain {
    pub profile: IeProfile,
    pub pml: Pml,
    pub d_pole: f64,
    pub decay: f64,
    pub center: [f64; 3],
}

impl ResolvedIeDomain {
    /// Whether the outer boundary (`t = 1`, mapped to `\tilde{R} = \infty`) should be clamped to the
    /// asymptotic value at infinity (`V = 0`, `T = T_inf`, `n x A = 0`).
    pub fn locks_outer_dirichlet(&self) -> bool {
        matches!(
            (self.profile, self.pml.coord),
            (IeProfile::Rational, PmlCoord::Spherical { .. })
                | (IeProfile::Rational, PmlCoord::Cylindrical { radial_on: true, .. })
        )
    }

    /// Evaluates the real coordinate-stretching Jacobian `J_IE(x_q)` at mesh point `x_mesh` (`[x, z, 0]` in 2D, `[x, y, z]` in 3D).
    pub fn eval_jacobian(&self, x_mesh: [f64; 3], dim: usize) -> PmlJacobian {
        let x_wafer = femcore::section::wafer_point(&x_mesh[..dim], dim);
        match self.profile {
            IeProfile::Exponential => {
                ie_exponential_jacobian(&self.pml, self.pml.bounds, self.decay, x_wafer, dim, false)
            }
            IeProfile::Rational => {
                // Touch femcore::pml::ie_rational_jacobian for penetration depth `s_dist` and Cartesian fallback,
                // and evaluate the exact pole-at-t=1 `[R_in, R_out) -> [R_in, \infty)` rational mapping
                // `\tilde{u}(t) - u_0 = d_pole * t / (1 - t)`, `s_u(t) = d_pole / (d * (1 - t)^2)`.
                let base_jac = ie_rational_jacobian(&self.pml, self.pml.bounds, self.d_pole, x_wafer, dim, false);
                eval_exact_rational_ie_jacobian(&self.pml, self.d_pole, x_wafer, dim, None, base_jac)
            }
        }
    }

    /// Evaluates the coordinate-stretching Jacobian `J_IE` at quadrature point `jac_k` of cell `iso`,
    /// using the element-interpolated radial coordinate `r_eff = \sum_a N_a R_a` and discrete radial normal
    /// `\hat{n}_r = \nabla_x r_eff / |\nabla_x r_eff|` so that constant-`r_eff` cell layers are 100% free of
    /// radial/tangential cross-contamination on both straight (`P1`) and curved (`P2`) cells.
    pub fn eval_cell_jacobian(
        &self,
        iso: &IsoparametricCell,
        jac_k: &femcore::fem_mixed::JacobianEval,
        sn: &[f64],
        sdn: &[[f64; 3]],
        dim: usize,
    ) -> PmlJacobian {
        if self.profile == IeProfile::Rational {
            match self.pml.coord {
                PmlCoord::Spherical { center, .. } => {
                    let mut r_eff = 0.0;
                    let mut dr_dxi = [0.0; 3];
                    for (a, node) in iso.nodes.iter().enumerate() {
                        let pw = femcore::section::wafer_point(&node[..dim], dim);
                        let ra = norm3(sub3(pw, center));
                        r_eff += sn[a] * ra;
                        dr_dxi[0] += sdn[a][0] * ra;
                        dr_dxi[1] += sdn[a][1] * ra;
                        dr_dxi[2] += sdn[a][2] * ra;
                    }
                    let g_mesh = jac_k.map_grad(dr_dxi);
                    let g_wafer = if dim == 2 { [g_mesh[0], 0.0, g_mesh[1]] } else { g_mesh };
                    let gn = norm3(g_wafer);
                    let r_hat_override = (gn > 1e-30).then(|| [g_wafer[0] / gn, g_wafer[1] / gn, g_wafer[2] / gn]);
                    let x_wafer = femcore::section::wafer_point(&jac_k.x[..dim], dim);
                    let base_jac = ie_rational_jacobian(&self.pml, self.pml.bounds, self.d_pole, x_wafer, dim, false);
                    return eval_exact_rational_ie_jacobian(
                        &self.pml,
                        self.d_pole,
                        x_wafer,
                        dim,
                        Some((r_eff, r_hat_override)),
                        base_jac,
                    );
                }
                PmlCoord::Cylindrical { center, axis, radial_on: true, .. } => {
                    let ax_norm = norm3(axis);
                    let a_hat = if ax_norm > 1e-30 {
                        [axis[0] / ax_norm, axis[1] / ax_norm, axis[2] / ax_norm]
                    } else if dim == 2 {
                        [0.0, 1.0, 0.0]
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    let mut r_eff = 0.0;
                    let mut dr_dxi = [0.0; 3];
                    for (a, node) in iso.nodes.iter().enumerate() {
                        let pw = femcore::section::wafer_point(&node[..dim], dim);
                        let d_vec = sub3(pw, center);
                        let zeta = dot3(d_vec, a_hat);
                        let r_vec =
                            [d_vec[0] - zeta * a_hat[0], d_vec[1] - zeta * a_hat[1], d_vec[2] - zeta * a_hat[2]];
                        let ra = norm3(r_vec);
                        r_eff += sn[a] * ra;
                        dr_dxi[0] += sdn[a][0] * ra;
                        dr_dxi[1] += sdn[a][1] * ra;
                        dr_dxi[2] += sdn[a][2] * ra;
                    }
                    let g_mesh = jac_k.map_grad(dr_dxi);
                    let g_wafer = if dim == 2 { [g_mesh[0], 0.0, g_mesh[1]] } else { g_mesh };
                    let g_axial = dot3(g_wafer, a_hat);
                    let g_perp = [
                        g_wafer[0] - g_axial * a_hat[0],
                        g_wafer[1] - g_axial * a_hat[1],
                        g_wafer[2] - g_axial * a_hat[2],
                    ];
                    let gn = norm3(g_perp);
                    let r_hat_override = (gn > 1e-30).then(|| [g_perp[0] / gn, g_perp[1] / gn, g_perp[2] / gn]);
                    let x_wafer = femcore::section::wafer_point(&jac_k.x[..dim], dim);
                    let base_jac = ie_rational_jacobian(&self.pml, self.pml.bounds, self.d_pole, x_wafer, dim, false);
                    return eval_exact_rational_ie_jacobian(
                        &self.pml,
                        self.d_pole,
                        x_wafer,
                        dim,
                        Some((r_eff, r_hat_override)),
                        base_jac,
                    );
                }
                _ => {}
            }
        }
        self.eval_jacobian(jac_k.x, dim)
    }
}

fn eval_exact_rational_ie_jacobian(
    spec: &Pml,
    d_pole: f64,
    x: [f64; 3],
    dim: usize,
    radial_override: Option<(f64, Option<[f64; 3]>)>,
    fallback: PmlJacobian,
) -> PmlJacobian {
    let one = c(1.0, 0.0);
    let zero = c(0.0, 0.0);
    let dp = d_pole.max(1e-30);

    match spec.coord {
        PmlCoord::Spherical { center, r_in, thickness } => {
            if thickness <= 0.0 {
                return fallback;
            }
            let r_vec = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
            let r_euclid = norm3(r_vec);
            let big_r = radial_override.map(|(re, _)| re).unwrap_or(r_euclid);
            if big_r < r_in - 1e-12 * thickness {
                return PmlJacobian::identity();
            }
            let t = ((big_r - r_in) / thickness).clamp(0.0, 1.0 - 1e-12);
            let om_t = (1.0 - t).max(1e-12);
            let delta_r = dp * t / om_t;
            let s_r_re = dp / (thickness * om_t * om_t);
            let s_big_r = c(s_r_re, 0.0);
            let r_scale = thickness.max(r_in.abs()).max(1e-300);
            let s_omega = if big_r > 1e-15 * r_scale { c((r_in + delta_r) / big_r, 0.0) } else { s_big_r };
            let r_hat = radial_override.and_then(|(_, rh)| rh).unwrap_or_else(|| {
                if r_euclid > 1e-15 * r_scale {
                    [r_vec[0] / r_euclid, r_vec[1] / r_euclid, r_vec[2] / r_euclid]
                } else {
                    [0.0, 0.0, 1.0]
                }
            });

            let inv_sr = one / s_big_r;
            let inv_so = one / s_omega;
            let diff_s = s_big_r - s_omega;
            let diff_inv = inv_sr - inv_so;

            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    let rab = r_hat[a] * r_hat[b];
                    let kron = if a == b { 1.0 } else { 0.0 };
                    j_w[a][b] = s_omega * kron + diff_s * rab;
                    j_inv_w[a][b] = inv_so * kron + diff_inv * rab;
                }
            }

            if dim == 2 {
                j_w[1][1] = one;
                j_inv_w[1][1] = one;
                let (j, j_inv) = (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w));
                return PmlJacobian { j, j_inv, det: s_big_r * s_omega, s_theta: s_omega, s_dist: t };
            }

            PmlJacobian { j: j_w, j_inv: j_inv_w, det: s_big_r * s_omega * s_omega, s_theta: s_omega, s_dist: t }
        }
        PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness } => {
            let ax_norm = norm3(axis);
            let a_hat = if ax_norm > 1e-30 {
                [axis[0] / ax_norm, axis[1] / ax_norm, axis[2] / ax_norm]
            } else if dim == 2 {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            let d_vec = sub3(x, center);
            let zeta = dot3(d_vec, a_hat);
            let r_vec = [d_vec[0] - zeta * a_hat[0], d_vec[1] - zeta * a_hat[1], d_vec[2] - zeta * a_hat[2]];
            let r_euclid = norm3(r_vec);
            let r = radial_override.map(|(re, _)| re).unwrap_or(r_euclid);

            let (mut s_r, mut s_phi, mut t_r) = (one, one, 0.0_f64);
            if radial_on && r_thickness > 0.0 && r >= r_in - 1e-12 * r_thickness {
                t_r = ((r - r_in) / r_thickness).clamp(0.0, 1.0 - 1e-12);
                let om_t = (1.0 - t_r).max(1e-12);
                let delta_r = dp * t_r / om_t;
                let sr_re = dp / (r_thickness * om_t * om_t);
                s_r = c(sr_re, 0.0);
                let r_scale = r_thickness.max(r_in.abs()).max(1e-300);
                s_phi = if r > 1e-15 * r_scale { c((r_in + delta_r) / r, 0.0) } else { s_r };
            }

            let (mut s_zeta, mut t_z) = (one, 0.0_f64);
            if axial != PmlDir::Off && z_thickness > 0.0 {
                let raw_tz = match axial {
                    PmlDir::Plus => (zeta - z_in) / z_thickness,
                    PmlDir::Minus => (z_in - zeta) / z_thickness,
                    PmlDir::Off => 0.0,
                };
                if raw_tz >= 0.0 {
                    t_z = raw_tz.clamp(0.0, 1.0 - 1e-12);
                    let c_z = if radial_on { 1.0 } else { 0.95 };
                    let om_tz = (1.0 - c_z * t_z).max(1e-12);
                    let dp_z = z_in.abs().max(0.5 * z_thickness).max(1e-30);
                    s_zeta = c(dp_z / (z_thickness * om_tz * om_tz), 0.0);
                }
            }

            let r_scale = r_thickness.max(r_in.abs()).max(1e-300);
            let r_hat = radial_override.and_then(|(_, rh)| rh).unwrap_or_else(|| {
                if r_euclid > 1e-15 * r_scale {
                    [r_vec[0] / r_euclid, r_vec[1] / r_euclid, r_vec[2] / r_euclid]
                } else {
                    [0.0, 0.0, 0.0]
                }
            });

            let inv_sr = one / s_r;
            let inv_sphi = one / s_phi;
            let inv_sz = one / s_zeta;
            let diff_r = s_r - s_phi;
            let diff_inv_r = inv_sr - inv_sphi;
            let diff_z = s_zeta - s_phi;
            let diff_inv_z = inv_sz - inv_sphi;

            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            for a in 0..3 {
                for b in 0..3 {
                    let rab = r_hat[a] * r_hat[b];
                    let aab = a_hat[a] * a_hat[b];
                    let kron = if a == b { 1.0 } else { 0.0 };
                    j_w[a][b] = s_phi * kron + diff_r * rab + diff_z * aab;
                    j_inv_w[a][b] = inv_sphi * kron + diff_inv_r * rab + diff_inv_z * aab;
                }
            }

            let s_dist = t_r.max(t_z);
            let det = if dim == 2 { s_r * s_phi } else { s_r * s_phi * s_zeta };
            let (j, j_inv) =
                if dim == 2 { (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w)) } else { (j_w, j_inv_w) };
            PmlJacobian { j, j_inv, det, s_theta: s_phi, s_dist }
        }
        PmlCoord::Cartesian => {
            let Some((lo, hi)) = spec.bounds else {
                return fallback;
            };
            let wafer_dirs = if dim == 2 { [spec.x, PmlDir::Off, spec.z] } else { [spec.x, spec.y, spec.z] };
            let mut s_diag = [one; 3];
            let mut t_max = 0.0_f64;
            const C_CART: f64 = 0.95;
            for a in 0..3 {
                let len = (hi[a] - lo[a]).max(1e-300);
                let t = match wafer_dirs[a] {
                    PmlDir::Off => continue,
                    PmlDir::Plus => ((x[a] - lo[a]) / len).clamp(0.0, 1.0),
                    PmlDir::Minus => ((hi[a] - x[a]) / len).clamp(0.0, 1.0),
                };
                t_max = t_max.max(t);
                let om_t = (1.0 - C_CART * t).max(1e-6);
                let s_re = 1.0 / (om_t * om_t);
                s_diag[a] = c(s_re, 0.0);
            }
            let mut j_w = [[zero; 3]; 3];
            let mut j_inv_w = [[zero; 3]; 3];
            for a in 0..3 {
                j_w[a][a] = s_diag[a];
                j_inv_w[a][a] = one / s_diag[a];
            }
            let det = if dim == 2 { s_diag[0] * s_diag[2] } else { s_diag[0] * s_diag[1] * s_diag[2] };
            let (j, j_inv) =
                if dim == 2 { (permute_wafer_to_2d(j_w), permute_wafer_to_2d(j_inv_w)) } else { (j_w, j_inv_w) };
            PmlJacobian { j, j_inv, det, s_theta: one, s_dist: t_max }
        }
    }
}

/// Resolves all domain-applied `BcKind::InfiniteElement` (and `spec.materials.pml`) entries against `mesh`.
pub(crate) fn resolve_ie_domains_mixed(
    mesh: &MixedDiscretisation,
    spec: &Problem,
) -> Result<BTreeMap<i32, ResolvedIeDomain>, String> {
    let mut raw_specs: BTreeMap<i32, InfiniteElementSpec> = BTreeMap::new();
    for bc in &spec.bcs {
        if let BcKind::InfiniteElement(ie) = &bc.kind {
            check_domain_labels_mixed(mesh, bc.labels.iter().copied())?;
            for &dom in &bc.labels {
                raw_specs.insert(dom, ie.clone());
            }
        }
    }
    for (&dom, pml) in &spec.materials.pml {
        if pml.is_active() && !raw_specs.contains_key(&dom) {
            let (coord, center, axis, r_in, thickness) = match pml.coord {
                PmlCoord::Cartesian => (IeCoordKind::Cartesian, None, [0.0, 0.0, 1.0], None, None),
                PmlCoord::Cylindrical { center, axis, r_in, r_thickness, .. } => (
                    IeCoordKind::Cylindrical,
                    Some(center),
                    axis,
                    (r_in > 0.0).then_some(r_in),
                    (r_thickness > 0.0).then_some(r_thickness),
                ),
                PmlCoord::Spherical { center, r_in, thickness } => (
                    IeCoordKind::Spherical,
                    Some(center),
                    [0.0, 0.0, 1.0],
                    (r_in > 0.0).then_some(r_in),
                    (thickness > 0.0).then_some(thickness),
                ),
            };
            raw_specs.insert(
                dom,
                InfiniteElementSpec {
                    coord,
                    profile: IeProfile::Rational,
                    x: pml.x,
                    y: pml.y,
                    z: pml.z,
                    center,
                    axis,
                    r_in,
                    thickness,
                    pole_dist: None,
                    decay: 30.0,
                    bounds_lo: pml.bounds.map(|b| b.0),
                    bounds_hi: pml.bounds.map(|b| b.1),
                },
            );
        }
    }

    if raw_specs.is_empty() {
        return Ok(BTreeMap::new());
    }

    let d = mesh.dim;
    let np = mesh.n_points();
    let ie_dom_set: BTreeSet<i32> = raw_specs.keys().copied().collect();

    // Mark vertices belonging to non-IE interior cells vs IE cells
    let mut in_phys = vec![false; np];
    let mut in_ie = vec![false; np];
    let mut phys_lo = [f64::INFINITY; 3];
    let mut phys_hi = [f64::NEG_INFINITY; 3];
    let mut has_phys = false;

    for c in 0..mesh.n_cells() {
        let is_ie = ie_dom_set.contains(&mesh.cell_labels[c]);
        for &v in mesh.cell(c) {
            let vu = v as usize;
            let pw = femcore::section::wafer_point(mesh.point(vu), d);
            if is_ie {
                in_ie[vu] = true;
            } else {
                in_phys[vu] = true;
                has_phys = true;
                for k in 0..3 {
                    phys_lo[k] = phys_lo[k].min(pw[k]);
                    phys_hi[k] = phys_hi[k].max(pw[k]);
                }
            }
        }
    }

    let (m_lo_raw, m_hi_raw) = mesh.bbox();
    let mesh_lo = femcore::section::wafer_point(&m_lo_raw[..d], d);
    let mesh_hi = femcore::section::wafer_point(&m_hi_raw[..d], d);
    let bbox_mid = if has_phys {
        [0.5 * (phys_lo[0] + phys_hi[0]), 0.5 * (phys_lo[1] + phys_hi[1]), 0.5 * (phys_lo[2] + phys_hi[2])]
    } else {
        [0.5 * (mesh_lo[0] + mesh_hi[0]), 0.5 * (mesh_lo[1] + mesh_hi[1]), 0.5 * (mesh_lo[2] + mesh_hi[2])]
    };

    // Interface vertices (shared between physical interior and IE, or inner boundary vertices)
    let mut interface_verts: Vec<[f64; 3]> =
        (0..np).filter(|&v| in_phys[v] && in_ie[v]).map(|v| femcore::section::wafer_point(mesh.point(v), d)).collect();
    if interface_verts.is_empty() {
        // Standalone IE shell mesh: use all exterior boundary vertices of the mesh
        let mut bverts = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if mesh.exterior_facets[f] {
                bverts.extend(mesh.facet(f).iter().copied());
            }
        }
        for v in bverts {
            interface_verts.push(femcore::section::wafer_point(mesh.point(v as usize), d));
        }
    }

    // Auto-select `[0, 0, 0]` vs `bbox_mid` when `ie.center` is `None` by checking radial uniformity on inner interface vertices
    let auto_center = {
        let radial_rel_spread = |cand: [f64; 3]| -> f64 {
            if interface_verts.is_empty() {
                return 1.0;
            }
            let mut dists: Vec<f64> = interface_verts.iter().map(|&p| norm3(sub3(p, cand))).collect();
            dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // Use the inner half of boundary vertices (in case standalone shell includes both R_in and R_out)
            let n_inner = (dists.len() / 2).max(1);
            let slice = &dists[..n_inner];
            let mean = slice.iter().sum::<f64>() / (n_inner as f64);
            if mean <= 1e-30 {
                return 1.0;
            }
            let var = slice.iter().map(|&r| (r - mean) * (r - mean)).sum::<f64>() / (n_inner as f64);
            var.sqrt() / mean
        };
        let sp_origin = radial_rel_spread([0.0, 0.0, 0.0]);
        let sp_mid = radial_rel_spread(bbox_mid);
        if sp_origin <= sp_mid + 1e-6 {
            [0.0, 0.0, 0.0]
        } else {
            bbox_mid
        }
    };

    let mut resolved = BTreeMap::new();
    for (dom, ie) in raw_specs {
        let mut dom_lo = [f64::INFINITY; 3];
        let mut dom_hi = [f64::NEG_INFINITY; 3];
        let mut dom_verts = BTreeSet::new();
        for c in 0..mesh.n_cells() {
            if mesh.cell_labels[c] == dom {
                for &v in mesh.cell(c) {
                    dom_verts.insert(v as usize);
                    let pw = femcore::section::wafer_point(mesh.point(v as usize), d);
                    for k in 0..3 {
                        dom_lo[k] = dom_lo[k].min(pw[k]);
                        dom_hi[k] = dom_hi[k].max(pw[k]);
                    }
                }
            }
        }
        let center = ie.center.unwrap_or(auto_center);

        let has_cart_dirs = ie.x != PmlDir::Off
            || ie.y != PmlDir::Off
            || ie.z != PmlDir::Off
            || ie.bounds_lo.is_some()
            || ie.bounds_hi.is_some();

        let eff_coord = match ie.coord {
            IeCoordKind::Auto => {
                if has_cart_dirs {
                    IeCoordKind::Cartesian
                } else if d == 2 {
                    if dom_lo[0] >= -1e-12 && dom_lo[0] <= 1e-6 * (dom_hi[0] - dom_lo[0]).max(1e-12) {
                        IeCoordKind::Spherical
                    } else {
                        IeCoordKind::Cylindrical
                    }
                } else {
                    IeCoordKind::Spherical
                }
            }
            other => other,
        };

        let (pml, d_pole) = match eff_coord {
            IeCoordKind::Cartesian | IeCoordKind::Auto => {
                let blo = ie.bounds_lo.unwrap_or(dom_lo);
                let bhi = ie.bounds_hi.unwrap_or(dom_hi);
                let (mut dir_x, mut dir_y, mut dir_z) = (ie.x, ie.y, ie.z);
                if dir_x == PmlDir::Off && dir_y == PmlDir::Off && dir_z == PmlDir::Off {
                    let ref_lo = if has_phys { phys_lo } else { mesh_lo };
                    let ref_hi = if has_phys { phys_hi } else { mesh_hi };
                    let tol = 1e-6
                        * ((mesh_hi[0] - mesh_lo[0])
                            .max(mesh_hi[1] - mesh_lo[1])
                            .max(mesh_hi[2] - mesh_lo[2])
                            .max(1e-12));
                    if blo[0] >= ref_hi[0] - tol {
                        dir_x = PmlDir::Plus;
                    } else if bhi[0] <= ref_lo[0] + tol {
                        dir_x = PmlDir::Minus;
                    }
                    if d == 3 {
                        if blo[1] >= ref_hi[1] - tol {
                            dir_y = PmlDir::Plus;
                        } else if bhi[1] <= ref_lo[1] + tol {
                            dir_y = PmlDir::Minus;
                        }
                    }
                    if blo[2] >= ref_hi[2] - tol {
                        dir_z = PmlDir::Plus;
                    } else if bhi[2] <= ref_lo[2] + tol {
                        dir_z = PmlDir::Minus;
                    }
                }
                let l_char = (bhi[0] - blo[0]).max(bhi[1] - blo[1]).max(bhi[2] - blo[2]).max(1e-12);
                let dp = ie.pole_dist.unwrap_or(0.5 * l_char);
                let pml = Pml {
                    x: dir_x,
                    y: dir_y,
                    z: dir_z,
                    coord: PmlCoord::Cartesian,
                    bounds: Some((blo, bhi)),
                    ..Pml::default()
                };
                (pml, dp)
            }
            IeCoordKind::Spherical => {
                let mut r_min = f64::INFINITY;
                let mut r_max = 0.0_f64;
                for &v in &dom_verts {
                    let pw = femcore::section::wafer_point(mesh.point(v), d);
                    let r = norm3(sub3(pw, center));
                    r_min = r_min.min(r);
                    r_max = r_max.max(r);
                }
                let r_in = ie.r_in.unwrap_or(r_min).max(1e-15);
                let thickness = ie.thickness.unwrap_or((r_max - r_in).max(1e-12));
                let dp = ie.pole_dist.unwrap_or(r_in);
                let pml = Pml { coord: PmlCoord::Spherical { center, r_in, thickness }, ..Pml::default() };
                (pml, dp)
            }
            IeCoordKind::Cylindrical => {
                let a_hat = if d == 2 && (ie.axis[0].abs() < 1e-12 && ie.axis[1].abs() < 1e-12) {
                    // In 2D planar `(x, z)`, the cylinder axis perpendicular to the cross-section is wafer `+y`.
                    [0.0, 1.0, 0.0]
                } else {
                    let an = norm3(ie.axis).max(1e-30);
                    [ie.axis[0] / an, ie.axis[1] / an, ie.axis[2] / an]
                };
                let mut r_min = f64::INFINITY;
                let mut r_max = 0.0_f64;
                for &v in &dom_verts {
                    let pw = femcore::section::wafer_point(mesh.point(v), d);
                    let dv = sub3(pw, center);
                    let zeta = dot3(dv, a_hat);
                    let rv = [dv[0] - zeta * a_hat[0], dv[1] - zeta * a_hat[1], dv[2] - zeta * a_hat[2]];
                    let r = norm3(rv);
                    r_min = r_min.min(r);
                    r_max = r_max.max(r);
                }
                let r_in = ie.r_in.unwrap_or(r_min).max(1e-15);
                let r_thickness = ie.thickness.unwrap_or((r_max - r_in).max(1e-12));
                let dp = ie.pole_dist.unwrap_or(r_in);
                let pml = Pml {
                    coord: PmlCoord::Cylindrical {
                        center,
                        axis: a_hat,
                        r_in,
                        r_thickness,
                        radial_on: true,
                        axial: ie.z,
                        z_in: 0.0,
                        z_thickness: r_thickness,
                    },
                    ..Pml::default()
                };
                (pml, dp)
            }
        };

        resolved.insert(dom, ResolvedIeDomain { profile: ie.profile, pml, d_pole, decay: ie.decay, center });
    }

    Ok(resolved)
}

/// Maps each exterior facet `f` of `mesh` to its owning cell index `c_f`.
pub(crate) fn build_facet_to_cell_mixed(mesh: &MixedDiscretisation) -> Vec<Option<usize>> {
    let nf = mesh.n_facets();
    let nc = mesh.n_cells();
    let mut out = vec![None; nf];
    if nf == 0 {
        return out;
    }
    if mesh.dim == 2 {
        let mut edge_to_cell = vec![usize::MAX; mesh.n_edges()];
        for c in 0..nc {
            for &eid in mesh.cell_edge_ids(c) {
                if edge_to_cell[eid as usize] == usize::MAX {
                    edge_to_cell[eid as usize] = c;
                }
            }
        }
        for f in 0..nf {
            let fv = mesh.facet(f);
            if let Some(eid) = mesh.edge_index(fv[0], fv[1]) {
                let c = edge_to_cell[eid];
                if c != usize::MAX {
                    out[f] = Some(c);
                }
            }
        }
    } else {
        let mut tri_to_cell = vec![usize::MAX; mesh.n_tri_faces()];
        let mut quad_to_cell = vec![usize::MAX; mesh.n_quad_faces()];
        for c in 0..nc {
            for &tf in mesh.cell_tri_face_ids(c) {
                if tri_to_cell[tf as usize] == usize::MAX {
                    tri_to_cell[tf as usize] = c;
                }
            }
            for &qf in mesh.cell_quad_face_ids(c) {
                if quad_to_cell[qf as usize] == usize::MAX {
                    quad_to_cell[qf as usize] = c;
                }
            }
        }
        for f in 0..nf {
            let fv = mesh.facet(f);
            match mesh.facet_type(f) {
                FacetType::Tri3 => {
                    if let Some(tf) = mesh.tri_face_index([fv[0], fv[1], fv[2]]) {
                        let c = tri_to_cell[tf];
                        if c != usize::MAX {
                            out[f] = Some(c);
                        }
                    }
                }
                FacetType::Quad4 => {
                    if let Some(qf) = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]) {
                        let c = quad_to_cell[qf];
                        if c != usize::MAX {
                            out[f] = Some(c);
                        }
                    }
                }
                FacetType::Edge2 => {}
            }
        }
    }
    out
}

/// Classifies outer exterior facets of Infinite Element domains into:
/// - `rational_dirichlet_facets`: outer facets (`s_dist > 0.85`) of `Rational` IE domains where `V = 0` / `T = T_inf` / `n x A = 0` holds at `\tilde{R} = \infty`,
/// - `exponential_robin_facets`: `(facet_idx, cell_idx, s_omega_factor)` of `Exponential` (or Cartesian) IE outer facets where Order-1 `far_field` Robin applies automatically.
pub(crate) fn classify_ie_outer_facets_mixed(
    mesh: &MixedDiscretisation,
    ie_map: &BTreeMap<i32, ResolvedIeDomain>,
) -> (Vec<usize>, Vec<(usize, usize, f64)>) {
    if ie_map.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let f2c = build_facet_to_cell_mixed(mesh);
    let mut rational_facets = Vec::new();
    let mut exp_robin_facets = Vec::new();

    for f in 0..mesh.n_facets() {
        if !mesh.exterior_facets[f] {
            continue;
        }
        let Some(c_f) = f2c[f] else {
            continue;
        };
        let dom = mesh.cell_labels[c_f];
        let Some(ie) = ie_map.get(&dom) else {
            continue;
        };
        let fv = mesh.facet(f);
        let mut cen = [0.0; 3];
        for &v in fv {
            let p = mesh.point3(v as usize);
            for k in 0..3 {
                cen[k] += p[k] / (fv.len() as f64);
            }
        }
        let jac = ie.eval_jacobian(cen, mesh.dim);
        if jac.s_dist > 0.85 {
            if ie.locks_outer_dirichlet() {
                rational_facets.push(f);
            } else {
                let s_omega = if mesh.dim == 3 { jac.s_theta.re.abs().max(1.0) } else { 1.0 };
                exp_robin_facets.push((f, c_f, s_omega));
            }
        }
    }
    (rational_facets, exp_robin_facets)
}

// ------------------------------------------------------------------ Mixed-cell Laplace, Source & Boundary Assemblers

/// Assembles scalar Laplace stiffness `K_{ab} = \int (\mathbf{J}_{IE}^{-T}\nabla\phi_a)^T \mathbf{K}_c (\mathbf{J}_{IE}^{-T}\nabla\phi_b) |\det\mathbf{J}_{IE}| dV`
/// over a [`MixedDiscretisation`].
pub(crate) fn assemble_laplace_ie_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    coeff: &[f64],
    tensor: Option<&[M3]>,
    ie_map: &BTreeMap<i32, ResolvedIeDomain>,
) -> Csr {
    let d = mesh.dim;
    let n = dof_map.n_dofs;
    let base_deg = 2 * dof_map.space.order();

    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grads: Vec<[f64; 3]> = Vec::new();
        let mut ke = Vec::new();

        for c in range {
            let ct = mesh.cell_type(c);
            let iso = IsoparametricCell::from_mesh(mesh, c);
            let ie_opt = ie_map.get(&mesh.cell_labels[c]);
            let deg = if ie_opt.is_some() || iso.curved { base_deg.max(6) } else { base_deg };
            let rule = quad_rule_mixed(ct, deg);
            let dofs = dof_map.cell_dofs(c);
            let nloc = dofs.len();
            ke.clear();
            ke.resize(nloc * nloc, 0.0);

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jac_k = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
                grads.clear();

                let mut weight = wq * jac_k.det.abs();
                if let Some(ie) = ie_opt {
                    let jac_ie = ie.eval_cell_jacobian(&iso, &jac_k, &sn, &sdn, d);
                    weight *= jac_ie.det.re.abs();
                    for a in 0..nloc {
                        let gp = jac_k.map_grad(dxi[a]);
                        let mut gt = [0.0; 3];
                        for i in 0..d {
                            let mut s = 0.0;
                            for m in 0..d {
                                s += jac_ie.j_inv[m][i].re * gp[m];
                            }
                            gt[i] = s;
                        }
                        grads.push(gt);
                    }
                } else {
                    for a in 0..nloc {
                        grads.push(jac_k.map_grad(dxi[a]));
                    }
                }

                if let Some(kt) = tensor {
                    let kc = &kt[c];
                    for a in 0..nloc {
                        let mut kg = [0.0; 3];
                        for i in 0..d {
                            for j in 0..d {
                                kg[i] += kc[i][j] * grads[a][j];
                            }
                        }
                        for b in 0..nloc {
                            let mut dot = 0.0;
                            for i in 0..d {
                                dot += kg[i] * grads[b][i];
                            }
                            ke[a * nloc + b] += weight * dot;
                        }
                    }
                } else {
                    let w_cf = weight * coeff[c];
                    for a in 0..nloc {
                        for b in 0..nloc {
                            let mut dot = 0.0;
                            for k in 0..d {
                                dot += grads[a][k] * grads[b][k];
                            }
                            ke[a * nloc + b] += w_cf * dot;
                        }
                    }
                }
            }
            coo.push_block(dofs, dofs, &ke);
        }
    })
}

/// Assembles volumetric source vector `f_a = \int_\Omega source(c) \phi_a dV` over a [`MixedDiscretisation`].
pub(crate) fn assemble_source_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    source: &[f64],
) -> Vec<f64> {
    let n = dof_map.n_dofs;
    let base_deg = 2 * dof_map.space.order();
    let mut f = vec![0.0; n];
    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();

    for c in 0..mesh.n_cells() {
        let sc = source[c];
        if sc == 0.0 {
            continue;
        }
        let ct = mesh.cell_type(c);
        let iso = IsoparametricCell::from_mesh(mesh, c);
        let deg = if iso.curved { base_deg.max(6) } else { base_deg };
        let rule = quad_rule_mixed(ct, deg);
        let dofs = dof_map.cell_dofs(c);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jac_k = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
            let w = wq * jac_k.det.abs() * sc;
            for (a, &gd) in dofs.iter().enumerate() {
                f[gd] += w * vals[a];
            }
        }
    }
    f
}

/// Assembles boundary load `f_a = \int_{\Gamma(labels)} value * \phi_a dS` over facets of a [`MixedDiscretisation`].
pub(crate) fn assemble_boundary_load_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    labels: &[i32],
    value: f64,
) -> Vec<f64> {
    let mut out = vec![0.0; dof_map.n_dofs];
    if value == 0.0 || labels.is_empty() {
        return out;
    }
    let deg = 2 * dof_map.space.order() + 2;
    let mut dofs = Vec::new();
    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, deg);
        dof_map.facet_dofs(mesh, f, &mut dofs);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (_x_q, ds_factor, _n_hat) = eval_facet_jacobian(mesh, f, xi);
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            let w = wq * ds_factor * value;
            for (a, &gd) in dofs.iter().enumerate() {
                out[gd] += w * nv[a];
            }
        }
    }
    out
}

/// Assembles boundary mass matrix `M_{ab} = \int_{\Gamma(labels)} coeff * \phi_a \phi_b dS` over facets of a [`MixedDiscretisation`].
pub(crate) fn assemble_boundary_mass_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    labels: &[i32],
    coeff: f64,
) -> Csr {
    let n = dof_map.n_dofs;
    if coeff == 0.0 || labels.is_empty() {
        return Csr::zeros(n, n);
    }
    let deg = 2 * dof_map.space.order() + 2;
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut ke = Vec::new();
    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, deg);
        dof_map.facet_dofs(mesh, f, &mut dofs);
        let nloc = dofs.len();
        ke.clear();
        ke.resize(nloc * nloc, 0.0);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (_x_q, ds_factor, _n_hat) = eval_facet_jacobian(mesh, f, xi);
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            let w = wq * ds_factor * coeff;
            for a in 0..nloc {
                for b in 0..nloc {
                    ke[a * nloc + b] += w * nv[a] * nv[b];
                }
            }
        }
        coo.push_block(&dofs, &dofs, &ke);
    }
    coo.to_csr()
}

/// Total boundary measure (length in 2D, area in 3D) of facets matching `labels`.
pub(crate) fn boundary_measure_mixed(mesh: &MixedDiscretisation, labels: &[i32]) -> f64 {
    let mut total = 0.0;
    for f in 0..mesh.n_facets() {
        if labels.contains(&mesh.facet_labels[f]) {
            let ft = mesh.facet_type(f);
            let rule = facet_rule_mixed(ft, 4);
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let (_, ds_factor, _) = eval_facet_jacobian(mesh, f, xi);
                total += wq * ds_factor;
            }
        }
    }
    total
}

fn cell_centroid_3d(mesh: &MixedDiscretisation, c: usize) -> [f64; 3] {
    let vs = mesh.cell(c);
    let mut cen = [0.0; 3];
    let inv = 1.0 / (vs.len().max(1) as f64);
    for &v in vs {
        let p = mesh.point3(v as usize);
        cen[0] += p[0] * inv;
        cen[1] += p[1] * inv;
        cen[2] += p[2] * inv;
    }
    cen
}

#[inline]
fn normal_perm_at(perm: &FarFieldPermittivity<'_>, dim: usize, cell: usize, n_hat: [f64; 3]) -> f64 {
    let val = match perm {
        FarFieldPermittivity::Scalar(s) => s.get(cell).copied().unwrap_or(1.0),
        FarFieldPermittivity::Tensor(t) => {
            if let Some(m) = t.get(cell) {
                let mut acc = 0.0;
                for i in 0..dim.min(3) {
                    for j in 0..dim.min(3) {
                        acc += n_hat[i] * m[i][j] * n_hat[j];
                    }
                }
                acc
            } else {
                1.0
            }
        }
    };
    if val.is_finite() {
        val.max(0.0)
    } else {
        0.0
    }
}

/// Deterministic area/length-weighted centroid of facets matching `labels` on a [`MixedDiscretisation`].
pub(crate) fn farfield_centroid_mixed(mesh: &MixedDiscretisation, labels: &[i32]) -> [f64; 3] {
    let mut sum_m = 0.0;
    let mut sum_c = [0.0; 3];
    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, 2);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (x_q, ds_factor, _) = eval_facet_jacobian(mesh, f, xi);
            let m = wq * ds_factor;
            if m.is_finite() && m > 0.0 && x_q.iter().all(|v| v.is_finite()) {
                sum_m += m;
                sum_c[0] += m * x_q[0];
                sum_c[1] += m * x_q[1];
                sum_c[2] += m * x_q[2];
            }
        }
    }
    if sum_m > 1e-30 && sum_m.is_finite() {
        return [sum_c[0] / sum_m, sum_c[1] / sum_m, sum_c[2] / sum_m];
    }
    let (lo, hi) = mesh.bbox();
    [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])]
}

/// Outward-facing positive measure of facets matching `labels` with respect to `center` on a [`MixedDiscretisation`].
pub(crate) fn farfield_positive_measure_mixed(mesh: &MixedDiscretisation, labels: &[i32], center: [f64; 3]) -> f64 {
    if labels.is_empty() {
        return 0.0;
    }
    let f2c = build_facet_to_cell_mixed(mesh);
    let mut total = 0.0;
    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let Some(c_f) = f2c[f] else {
            continue;
        };
        let p_cell = cell_centroid_3d(mesh, c_f);
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, 2);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (x_q, ds_factor, mut n_hat) = eval_facet_jacobian(mesh, f, xi);
            if dot3(n_hat, sub3(x_q, p_cell)) < 0.0 {
                n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
            }
            if dot3(sub3(x_q, center), n_hat) > 0.0 {
                total += wq * ds_factor;
            }
        }
    }
    total
}

///Evaluates position, area factor, outward unit normal, and physical tangent vectors `(t0, t1)` on 3D facet `f` (`Tri3` or `Quad4`).
fn eval_facet_tangents_3d(
    mesh: &MixedDiscretisation,
    f: usize,
    p_cell: [f64; 3],
    xi: &[f64; 2],
) -> ([f64; 3], f64, [f64; 3], [f64; 3], [f64; 3]) {
    let ft = mesh.facet_type(f);
    let fv = mesh.facet(f);
    let corners: Vec<[f64; 3]> = fv.iter().map(|&v| mesh.point3(v as usize)).collect();
    let (space, nodes) = match ft {
        FacetType::Tri3 => {
            let e01 = mesh.edge_index(fv[0], fv[1]).unwrap();
            let e02 = mesh.edge_index(fv[0], fv[2]).unwrap();
            let e12 = mesh.edge_index(fv[1], fv[2]).unwrap();
            let any_curved =
                mesh.edge_node(e01).is_some() || mesh.edge_node(e02).is_some() || mesh.edge_node(e12).is_some();
            if any_curved {
                let mid = |eid: usize, i: usize, j: usize| {
                    mesh.edge_node(eid).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ])
                };
                (
                    LagrangeSpace::Order2Serendipity,
                    vec![corners[0], corners[1], corners[2], mid(e01, 0, 1), mid(e02, 0, 2), mid(e12, 1, 2)],
                )
            } else {
                (LagrangeSpace::Order1, corners)
            }
        }
        FacetType::Quad4 => {
            let qf_id = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
            let eids = [
                mesh.edge_index(fv[0], fv[1]).unwrap(),
                mesh.edge_index(fv[1], fv[2]).unwrap(),
                mesh.edge_index(fv[2], fv[3]).unwrap(),
                mesh.edge_index(fv[3], fv[0]).unwrap(),
            ];
            let has_qfn = mesh.quad_face_node(qf_id).is_some();
            let any_curved = has_qfn || eids.iter().any(|&e| mesh.edge_node(e).is_some());
            if any_curved {
                let mut ns = corners.clone();
                for (k, &(i, j)) in QUAD_EDGES.iter().enumerate() {
                    let m = mesh.edge_node(eids[k]).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ]);
                    ns.push(m);
                }
                if has_qfn {
                    ns.push(mesh.quad_face_node(qf_id).unwrap());
                    (LagrangeSpace::Order2Tensor, ns)
                } else {
                    (LagrangeSpace::Order2Serendipity, ns)
                }
            } else {
                (LagrangeSpace::Order1, corners)
            }
        }
        FacetType::Edge2 => unreachable!("eval_facet_tangents_3d requires 3D facet"),
    };
    let (n, dn) = facet_lagrange_shape(ft, space, xi);
    let (mut x, mut t0, mut t1) = ([0.0; 3], [0.0; 3], [0.0; 3]);
    for (a, p) in nodes.iter().enumerate() {
        for k in 0..3 {
            x[k] += n[a] * p[k];
            t0[k] += dn[a][0] * p[k];
            t1[k] += dn[a][1] * p[k];
        }
    }
    let cx = cross3(t0, t1);
    let ds = norm3(cx);
    let mut n_hat = if ds > 0.0 { [cx[0] / ds, cx[1] / ds, cx[2] / ds] } else { [0.0; 3] };
    if dot3(n_hat, sub3(x, p_cell)) < 0.0 {
        n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
    }
    (x, ds, n_hat, t0, t1)
}

/// Assembles Order-1 or 3D Order-2 Laplace–Beltrami `far_field` Robin boundary stiffness on a [`MixedDiscretisation`]
/// (supporting 2D `Edge2` and 3D `Tri3` + `Quad4` facets).
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_farfield_robin_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    labels: &[i32],
    bc_order: usize,
    decay_power: f64,
    center: [f64; 3],
    scale0: f64,
    perm: FarFieldPermittivity<'_>,
) -> Csr {
    let n = dof_map.n_dofs;
    if labels.is_empty() {
        return Csr::zeros(n, n);
    }
    let gamma = if decay_power.is_finite() && decay_power > 0.0 { decay_power } else { 1.0 };
    let f2c = build_facet_to_cell_mixed(mesh);
    let deg = 2 * dof_map.space.order() + 2;
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut ke = Vec::new();
    let mut grad_surf: Vec<[f64; 3]> = Vec::new();

    let (lo, hi) = mesh.bbox();
    let r_inf_2d = ((hi[0] - lo[0]).hypot(hi[1] - lo[1])).max(1e-12);

    for f in 0..mesh.n_facets() {
        if !labels.contains(&mesh.facet_labels[f]) {
            continue;
        }
        let Some(c_f) = f2c[f] else {
            continue;
        };
        let p_cell = cell_centroid_3d(mesh, c_f);
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, deg);
        dof_map.facet_dofs(mesh, f, &mut dofs);
        let nloc = dofs.len();
        ke.clear();
        ke.resize(nloc * nloc, 0.0);

        if mesh.dim == 2 {
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let (x_q, ds_factor, mut n_hat) = eval_facet_jacobian(mesh, f, xi);
                if dot3(n_hat, sub3(x_q, p_cell)) < 0.0 {
                    n_hat = [-n_hat[0], -n_hat[1], 0.0];
                }
                let r_vec = [x_q[0] - center[0], x_q[1] - center[1], 0.0];
                let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1];
                let rn = dot3(r_vec, n_hat);
                if r2 <= 1e-30 || rn <= 0.0 {
                    continue;
                }
                let r = r2.sqrt();
                let r_scale = r_inf_2d.max(r);
                let denom_ln = (2.0 * r_scale / r).ln().max(1.0);
                let inv_r_eff = rn / r2;
                let eps_nf = normal_perm_at(&perm, 2, c_f, n_hat);
                let alpha = scale0 * eps_nf * (gamma * inv_r_eff) / (2.0 * denom_ln);
                let ds = wq * ds_factor * alpha;
                let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
                for a in 0..nloc {
                    for b in 0..nloc {
                        ke[a * nloc + b] += ds * nv[a] * nv[b];
                    }
                }
            }
        } else {
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let (x_q, ds_factor, n_hat, t0, t1) = eval_facet_tangents_3d(mesh, f, p_cell, xi);
                if ds_factor <= 1e-30 {
                    continue;
                }
                let r_vec = sub3(x_q, center);
                let r2 = dot3(r_vec, r_vec);
                let rn = dot3(r_vec, n_hat);
                if r2 <= 1e-30 || rn <= 0.0 {
                    continue;
                }
                let r = r2.sqrt();
                let eps_nf = normal_perm_at(&perm, 3, c_f, n_hat);
                if eps_nf <= 0.0 {
                    continue;
                }
                let ds = wq * ds_factor;
                let (nv, dn_fe) = facet_lagrange_shape(ft, dof_map.space, xi);

                // Order-1 Robin contribution
                let alpha = scale0 * eps_nf * (gamma * rn / r2);
                let w_mass = ds * alpha;
                for a in 0..nloc {
                    for b in 0..nloc {
                        ke[a * nloc + b] += w_mass * nv[a] * nv[b];
                    }
                }

                // Order-2 symmetric Laplace-Beltrami contribution in 3D
                if bc_order >= 2 && rn > EPS_GEOM * r {
                    let g00 = dot3(t0, t0);
                    let g01 = dot3(t0, t1);
                    let g11 = dot3(t1, t1);
                    let det_g = g00 * g11 - g01 * g01;
                    if det_g > 1e-60 {
                        let inv_det = 1.0 / det_g;
                        let t0_contra = [
                            (g11 * t0[0] - g01 * t1[0]) * inv_det,
                            (g11 * t0[1] - g01 * t1[1]) * inv_det,
                            (g11 * t0[2] - g01 * t1[2]) * inv_det,
                        ];
                        let t1_contra = [
                            (-g01 * t0[0] + g00 * t1[0]) * inv_det,
                            (-g01 * t0[1] + g00 * t1[1]) * inv_det,
                            (-g01 * t0[2] + g00 * t1[2]) * inv_det,
                        ];
                        let r_hat = [r_vec[0] / r, r_vec[1] / r, r_vec[2] / r];
                        let cos_th = rn / r;
                        let r_gamma =
                            [r_hat[0] - cos_th * n_hat[0], r_hat[1] - cos_th * n_hat[1], r_hat[2] - cos_th * n_hat[2]];
                        let rg2 = dot3(r_gamma, r_gamma);
                        let proj_scale = 1.0 / (rg2 + RADIAL_PROJ_REG * RADIAL_PROJ_REG);

                        grad_surf.clear();
                        for a in 0..nloc {
                            let d0 = dn_fe[a][0];
                            let d1 = dn_fe[a][1];
                            let g_raw = [
                                d0 * t0_contra[0] + d1 * t1_contra[0],
                                d0 * t0_contra[1] + d1 * t1_contra[1],
                                d0 * t0_contra[2] + d1 * t1_contra[2],
                            ];
                            let gr = dot3(g_raw, r_gamma) * proj_scale;
                            grad_surf.push([
                                g_raw[0] - gr * r_gamma[0],
                                g_raw[1] - gr * r_gamma[1],
                                g_raw[2] - gr * r_gamma[2],
                            ]);
                        }
                        let beta = scale0 * eps_nf * r2 / (2.0 * rn.max(EPS_GEOM * r));
                        let w_lb = ds * beta;
                        for a in 0..nloc {
                            for b in 0..nloc {
                                ke[a * nloc + b] += w_lb * dot3(grad_surf[a], grad_surf[b]);
                            }
                        }
                    }
                }
            }
        }
        coo.push_block(&dofs, &dofs, &ke);
    }
    coo.to_csr()
}

/// Assembles automatic Order-1 `far_field` Robin on the outer facets of `Exponential` (or Cartesian) Infinite Element domains.
pub(crate) fn assemble_ie_outer_robin_mixed(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    exp_facets: &[(usize, usize, f64)],
    ie_map: &BTreeMap<i32, ResolvedIeDomain>,
    scale0: f64,
    perm: FarFieldPermittivity<'_>,
) -> Csr {
    let n = dof_map.n_dofs;
    if exp_facets.is_empty() {
        return Csr::zeros(n, n);
    }
    let deg = 2 * dof_map.space.order() + 2;
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut ke = Vec::new();
    let (lo, hi) = mesh.bbox();
    let r_inf_2d = ((hi[0] - lo[0]).hypot(hi[1] - lo[1])).max(1e-12);

    for &(f, c_f, s_omega) in exp_facets {
        let dom = mesh.cell_labels[c_f];
        let center = ie_map.get(&dom).map_or([0.0; 3], |ie| ie.center);
        let p_cell = cell_centroid_3d(mesh, c_f);
        let ft = mesh.facet_type(f);
        let rule = facet_rule_mixed(ft, deg);
        dof_map.facet_dofs(mesh, f, &mut dofs);
        let nloc = dofs.len();
        ke.clear();
        ke.resize(nloc * nloc, 0.0);

        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (x_q, ds_factor, mut n_hat) = eval_facet_jacobian(mesh, f, xi);
            if dot3(n_hat, sub3(x_q, p_cell)) < 0.0 {
                n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
            }
            let r_vec = sub3(x_q, center);
            let r2 = if mesh.dim == 2 { r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] } else { dot3(r_vec, r_vec) };
            let rn = dot3(r_vec, n_hat);
            if r2 <= 1e-30 || rn <= 0.0 {
                continue;
            }
            let eps_nf = normal_perm_at(&perm, mesh.dim, c_f, n_hat);
            let alpha = if mesh.dim == 2 {
                let r = r2.sqrt();
                let r_scale = (r_inf_2d * s_omega.max(1.0)).max(r);
                let denom_ln = (2.0 * r_scale / r).ln().max(1.0);
                scale0 * eps_nf * (rn / r2) / (2.0 * denom_ln)
            } else {
                scale0 * eps_nf * s_omega * (rn / r2)
            };
            let ds = wq * ds_factor * alpha;
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            for a in 0..nloc {
                for b in 0..nloc {
                    ke[a * nloc + b] += ds * nv[a] * nv[b];
                }
            }
        }
        coo.push_block(&dofs, &dofs, &ke);
    }
    coo.to_csr()
}

/// Solves a real SPD system `K u = f` with Dirichlet constraints `(fixed, vals)` on a [`MixedDiscretisation`].
pub(crate) fn solve_spd_mixed(
    _solver_pref: LinearSolver,
    k: &Csr,
    f: &[f64],
    fixed: &[usize],
    vals: &[f64],
) -> Result<(Vec<f64>, Option<String>), String> {
    let n = k.nrows();
    let mut u = vec![0.0; n];
    for (&d, &v) in fixed.iter().zip(vals) {
        u[d] = v;
    }
    let free = free_dofs(n, fixed);
    let n_free = free.len();
    if n_free == 0 {
        return Ok((u, Some("0 free dofs: determined by boundary constraints".into())));
    }
    let kff = k.submatrix(&free, &free);
    let rhs: Vec<f64> = if fixed.is_empty() || vals.iter().all(|&v| v == 0.0) {
        free.iter().map(|&i| f[i]).collect()
    } else {
        let kfx = k.submatrix(&free, fixed);
        let kx = kfx.matvec(vals);
        free.iter().enumerate().map(|(a, &i)| f[i] - kx[a]).collect()
    };
    let solver = DirectSolver::new(&kff, true)?;
    let u_free = solver.solve(&rhs);
    for (a, &i) in free.iter().enumerate() {
        u[i] = u_free[a];
    }
    Ok((u, Some(format!("direct Cholesky ({n_free} free mixed-cell dofs)"))))
}

/// Evaluates `scale[mc] * (tensor[mc] * \nabla_x u)` on the native parent cells of `mesh` at every
/// sub-simplex node of `proj` (`n_simplex_cells * nloc(target_order) * dim` values).
pub(crate) fn sample_mixed_gradients_on_simplex(
    mesh: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    proj: &MixedToSimplexProjector,
    u_mixed: &[f64],
    tensor: Option<&[M3]>,
    scale: &[f64],
    target_order: usize,
) -> Vec<f64> {
    let d = mesh.dim;
    let iso_cells: Vec<IsoparametricCell> =
        (0..mesh.n_cells()).map(|c| IsoparametricCell::from_mesh(mesh, c)).collect();
    let mut sn = Vec::new();
    let mut sdn = Vec::new();
    let mut vals = Vec::new();
    let mut dxi = Vec::new();

    let samples: Vec<[f64; 3]> = proj.sample_simplex_cell_nodes(mesh, target_order, |mc, xi| {
        let ct = mesh.cell_type(mc);
        let jac = iso_cells[mc].eval_into(xi, &mut sn, &mut sdn);
        lagrange_shape_into(ct, dof_map.space, xi, &mut vals, &mut dxi);
        let cdofs = dof_map.cell_dofs(mc);
        let mut g_ref = [0.0; 3];
        for (a, &gd) in cdofs.iter().enumerate() {
            let ua = u_mixed[gd];
            for k in 0..d {
                g_ref[k] += ua * dxi[a][k];
            }
        }
        let g_phys = jac.map_grad(g_ref);
        let sc = scale[mc];
        let mut out = [0.0; 3];
        if let Some(kt) = tensor {
            let kc = &kt[mc];
            for i in 0..d {
                let mut s = 0.0;
                for j in 0..d {
                    s += kc[i][j] * g_phys[j];
                }
                out[i] = sc * s;
            }
        } else {
            for i in 0..d {
                out[i] = sc * g_phys[i];
            }
        }
        out
    });

    let mut flat = Vec::with_capacity(samples.len() * d);
    for s in samples {
        flat.extend_from_slice(&s[..d]);
    }
    flat
}

/// Returns `true` if `mesh` carries any quadratic mid-edge node offsets (`edge_nodes` non-empty).
pub(crate) fn mesh_has_curved_edges(mesh: &Mesh) -> bool {
    match mesh {
        Mesh::D2(m) => !m.edge_nodes.is_empty(),
        Mesh::D3(m) => !m.edge_nodes.is_empty(),
    }
}

/// 5-point Gauss-Legendre rule on the reference segment (exact degree 9; barycentric `[lam0, lam1]`).
pub(crate) fn quad_seg_high() -> Quad {
    const T5: [f64; 5] =
        [-0.906_179_845_938_664, -0.538_469_310_105_683_1, 0.0, 0.538_469_310_105_683_1, 0.906_179_845_938_664];
    const W5: [f64; 5] = [
        0.236_926_885_056_189_1,
        0.478_628_670_499_366_5,
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_5,
        0.236_926_885_056_189_1,
    ];
    let mut lam = Vec::with_capacity(5);
    let mut w = Vec::with_capacity(5);
    for (&t, &wi) in T5.iter().zip(&W5) {
        let s = 0.5 * (t + 1.0);
        lam.push(vec![1.0 - s, s]);
        w.push(0.5 * wi);
    }
    Quad { lam, w }
}

/// 4x4 collapsed Gauss-Legendre rule on the reference triangle (exact degree 7; weights sum to 1).
pub(crate) fn quad_tri_high() -> Quad {
    const T4: [f64; 4] =
        [-0.861_136_311_594_052_6, -0.339_981_043_584_856_3, 0.339_981_043_584_856_3, 0.861_136_311_594_052_6];
    const W4: [f64; 4] =
        [0.347_854_845_137_453_8, 0.652_145_154_862_546_1, 0.652_145_154_862_546_1, 0.347_854_845_137_453_8];
    let mut lam = Vec::with_capacity(16);
    let mut w = Vec::with_capacity(16);
    for (&tu, &wu) in T4.iter().zip(&W4) {
        let u = 0.5 * (1.0 + tu);
        let du = 0.5 * wu;
        for (&tv, &wv) in T4.iter().zip(&W4) {
            let v = 0.5 * (1.0 + tv);
            let dv = 0.5 * wv;
            let l1 = u;
            let l2 = (1.0 - u) * v;
            let l0 = 1.0 - l1 - l2;
            lam.push(vec![l0, l1, l2]);
            w.push(2.0 * (1.0 - u) * du * dv);
        }
    }
    Quad { lam, w }
}

/// Physical coordinate of scalar DOF `d` (`0..n_points` vertices, then curved or straight edge midpoints).
pub(crate) fn dof_physical_coord<const D: usize>(mesh: &Discretisation, d: usize, order: usize) -> [f64; D] {
    let nv = mesh.n_points();
    if d < nv || order < 2 {
        let p = mesh.point(d);
        let mut out = [0.0; D];
        out.copy_from_slice(&p[..D]);
        return out;
    }
    let e = d - nv;
    if let Some(en) = mesh.edge_node(e) {
        let mut out = [0.0; D];
        out.copy_from_slice(&en[..D]);
        return out;
    }
    let (a, b) = mesh.edges()[e];
    let pa = mesh.point(a as usize);
    let pb = mesh.point(b as usize);
    let mut out = [0.0; D];
    for k in 0..D {
        out[k] = 0.5 * (pa[k] + pb[k]);
    }
    out
}

/// Assembles `int_Omega f(x, t, dom) * phi_i(x) dV` via monomorphized [`SourceField<D>`].
pub(crate) fn assemble_source_field<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    order: usize,
    source: &S,
    t: f64,
    base_uniform: Option<&[f64]>,
) -> Result<Vec<f64>, String> {
    if mesh.dim != D {
        return Err(format!("SourceField<{D}> called on {}D mesh", mesh.dim));
    }
    let n = mesh.n_dofs(order);
    let mut f = vec![0.0; n];
    let nloc = mesh.nloc(order);
    let geo = geometry(mesh);
    let mut dofs = Vec::new();
    if D == 2 {
        let qrule = quad_tri_high();
        let shapes: Vec<Vec<f64>> = qrule.lam.iter().map(|lam| shape(lam, order).0).collect();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let dom = mesh.cell_labels[c] as u32;
            let extra = base_uniform.map_or(0.0, |u| u[c]);
            let pts = [mesh.point(v[0] as usize), mesh.point(v[1] as usize), mesh.point(v[2] as usize)];
            mesh.cell_dofs(c, order, &mut dofs);
            let vol = geo.vol[c];
            for q in 0..qrule.w.len() {
                let lam = &qrule.lam[q];
                let mut xq = [0.0; D];
                for (b, pb) in pts.iter().enumerate() {
                    for k in 0..D {
                        xq[k] += lam[b] * pb[k];
                    }
                }
                let val = source.volume_rhs_re(xq, t, dom)[0] + extra;
                let wq = vol * qrule.w[q] * val;
                let nv = &shapes[q];
                for a in 0..nloc {
                    f[dofs[a]] += wq * nv[a];
                }
            }
        }
    } else {
        let qstraight = quad_rule_collapsed(4);
        let qcurved = quad_rule_collapsed(8);
        let s_straight: Vec<Vec<f64>> = qstraight.lam.iter().map(|lam| shape(lam, order).0).collect();
        let s_curved: Vec<Vec<f64>> = qcurved.lam.iter().map(|lam| shape(lam, order).0).collect();
        for c in 0..mesh.n_cells() {
            let v = mesh.cell(c);
            let dom = mesh.cell_labels[c] as u32;
            let extra = base_uniform.map_or(0.0, |u| u[c]);
            mesh.cell_dofs(c, order, &mut dofs);
            if mesh.cell_is_curved(c) {
                let ct = CurvedTet::new(mesh.cell_geometry_nodes(c));
                for q in 0..qcurved.w.len() {
                    let lam = [qcurved.lam[q][0], qcurved.lam[q][1], qcurved.lam[q][2], qcurved.lam[q][3]];
                    let x = ct.position(&lam);
                    let mut xq = [0.0; D];
                    xq.copy_from_slice(&x[..D]);
                    let (_, det) = ct.jacobian(&lam);
                    let val = source.volume_rhs_re(xq, t, dom)[0] + extra;
                    let wq = (qcurved.w[q] * det.abs() / 6.0) * val;
                    let nv = &s_curved[q];
                    for a in 0..nloc {
                        f[dofs[a]] += wq * nv[a];
                    }
                }
            } else {
                let pts = [
                    mesh.point(v[0] as usize),
                    mesh.point(v[1] as usize),
                    mesh.point(v[2] as usize),
                    mesh.point(v[3] as usize),
                ];
                let vol = geo.vol[c];
                for q in 0..qstraight.w.len() {
                    let lam = &qstraight.lam[q];
                    let mut xq = [0.0; D];
                    for (b, pb) in pts.iter().enumerate() {
                        for k in 0..D {
                            xq[k] += lam[b] * pb[k];
                        }
                    }
                    let val = source.volume_rhs_re(xq, t, dom)[0] + extra;
                    let wq = vol * qstraight.w[q] * val;
                    let nv = &s_straight[q];
                    for a in 0..nloc {
                        f[dofs[a]] += wq * nv[a];
                    }
                }
            }
        }
    }
    Ok(f)
}

/// Collects Dirichlet `(fixed_dofs, values)` from [`SourceField<D>`] on facets whose label is in `labels`.
pub(crate) fn dirichlet_constraints_field<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    labels: &BTreeSet<i32>,
    order: usize,
    source: &S,
    t: f64,
) -> (Vec<usize>, Vec<f64>) {
    let mut acc: BTreeMap<usize, (f64, usize)> = BTreeMap::new();
    let mut dofs = Vec::new();
    for f in 0..mesh.n_facets() {
        let lab = mesh.facet_labels[f];
        if labels.contains(&lab) {
            mesh.facet_dofs(f, order, &mut dofs);
            for &d in &dofs {
                let x = dof_physical_coord::<D>(mesh, d, order);
                let val = source.dirichlet_re(x, t, lab as u32)[0];
                let e = acc.entry(d).or_insert((0.0, 0));
                e.0 += val;
                e.1 += 1;
            }
        }
    }
    let fixed = acc.keys().copied().collect();
    let vals = acc.values().map(|(s, cnt)| s / *cnt as f64).collect();
    (fixed, vals)
}

fn facet_opp_vertices_2d(mesh: &Discretisation) -> HashMap<[u32; 2], u32> {
    let mut map = HashMap::new();
    for c in 0..mesh.n_cells() {
        let cv = mesh.cell(c);
        for (a, b, opp) in [(cv[0], cv[1], cv[2]), (cv[1], cv[2], cv[0]), (cv[2], cv[0], cv[1])] {
            map.entry([a.min(b), a.max(b)]).or_insert(opp);
        }
    }
    map
}

fn facet_opp_vertices_3d(mesh: &Discretisation) -> HashMap<[u32; 3], u32> {
    let mut map = HashMap::new();
    for c in 0..mesh.n_cells() {
        let cv = mesh.cell(c);
        for skip in 0..4 {
            let [i0, i1, i2] = match skip {
                0 => [1, 2, 3],
                1 => [0, 2, 3],
                2 => [0, 1, 3],
                _ => [0, 1, 2],
            };
            let mut key = [cv[i0], cv[i1], cv[i2]];
            key.sort_unstable();
            map.entry(key).or_insert(cv[skip]);
        }
    }
    map
}

fn facet_3d_geom_nodes(mesh: &Discretisation, fv: &[u32]) -> [[f64; 3]; 6] {
    let pt = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], p[2]]
    };
    let mid = |a: u32, b: u32| -> [f64; 3] {
        if let Some(en) = mesh.edge_node(mesh.edge_index(a, b)) {
            return en;
        }
        let (pa, pb) = (pt(a), pt(b));
        [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]
    };
    [pt(fv[0]), pt(fv[1]), pt(fv[2]), mid(fv[0], fv[1]), mid(fv[0], fv[2]), mid(fv[1], fv[2])]
}

/// Assembles `int_{Gamma} g(x, n, t, bnd) phi_i dS` via monomorphized [`SourceField<D>`], supporting curved 3D facets.
pub(crate) fn assemble_boundary_load_field<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    labels: &[i32],
    order: usize,
    source: &S,
    t: f64,
) -> Result<Vec<f64>, String> {
    if mesh.dim != D {
        return Err(format!("SourceField<{D}> called on {}D mesh", mesh.dim));
    }
    let n = mesh.n_dofs(order);
    let mut f_vec = vec![0.0; n];
    let label_set: BTreeSet<i32> = labels.iter().copied().collect();
    if label_set.is_empty() {
        return Ok(f_vec);
    }
    let mut dofs = Vec::new();
    if D == 2 {
        let opp_map = facet_opp_vertices_2d(mesh);
        let qrule = quad_seg_high();
        let shapes: Vec<Vec<f64>> = qrule.lam.iter().map(|lam| shape(lam, order).0).collect();
        for f_idx in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f_idx];
            if !label_set.contains(&lab) {
                continue;
            }
            let fv = mesh.facet(f_idx);
            let p0 = mesh.point(fv[0] as usize);
            let p1 = mesh.point(fv[1] as usize);
            let dx = p1[0] - p0[0];
            let dy = p1[1] - p0[1];
            let len = dx.hypot(dy);
            if len <= 1e-30 {
                continue;
            }
            let mut n_hat = [-dy / len, dx / len];
            if let Some(&opp) = opp_map.get(&[fv[0].min(fv[1]), fv[0].max(fv[1])]) {
                let po = mesh.point(opp as usize);
                let cen = [0.5 * (p0[0] + p1[0]), 0.5 * (p0[1] + p1[1])];
                if (cen[0] - po[0]) * n_hat[0] + (cen[1] - po[1]) * n_hat[1] < 0.0 {
                    n_hat = [-n_hat[0], -n_hat[1]];
                }
            }
            mesh.facet_dofs(f_idx, order, &mut dofs);
            for (q, lam) in qrule.lam.iter().enumerate() {
                let xq = [lam[0] * p0[0] + lam[1] * p1[0], lam[0] * p0[1] + lam[1] * p1[1]];
                let mut x_arr = [0.0; D];
                x_arr.copy_from_slice(&xq[..D]);
                let mut n_arr = [0.0; D];
                n_arr.copy_from_slice(&n_hat[..D]);
                let gval = source.neumann_robin_re(x_arr, n_arr, t, lab as u32)[0];
                let wq = len * qrule.w[q] * gval;
                for (a, &d) in dofs.iter().enumerate() {
                    f_vec[d] += wq * shapes[q][a];
                }
            }
        }
    } else {
        let opp_map = facet_opp_vertices_3d(mesh);
        let qrule = quad_tri_high();
        let shapes_fe: Vec<Vec<f64>> = qrule.lam.iter().map(|lam| shape(lam, order).0).collect();
        let shapes_g2: Vec<(Vec<f64>, Vec<f64>)> = qrule.lam.iter().map(|lam| shape(lam, 2)).collect();
        for f_idx in 0..mesh.n_facets() {
            let lab = mesh.facet_labels[f_idx];
            if !label_set.contains(&lab) {
                continue;
            }
            let fv = mesh.facet(f_idx);
            let mut key = [fv[0], fv[1], fv[2]];
            key.sort_unstable();
            let p_opp = opp_map.get(&key).map(|&opp| {
                let p = mesh.point(opp as usize);
                [p[0], p[1], p[2]]
            });
            mesh.facet_dofs(f_idx, order, &mut dofs);
            if mesh.has_curved() {
                let x6 = facet_3d_geom_nodes(mesh, fv);
                for (q, _) in qrule.lam.iter().enumerate() {
                    let (ref nv2, ref dn2) = shapes_g2[q];
                    let mut xq = [0.0; 3];
                    let mut t1 = [0.0; 3];
                    let mut t2 = [0.0; 3];
                    for a in 0..6 {
                        let d1 = dn2[a * 3 + 1] - dn2[a * 3];
                        let d2 = dn2[a * 3 + 2] - dn2[a * 3];
                        for k in 0..3 {
                            xq[k] += nv2[a] * x6[a][k];
                            t1[k] += d1 * x6[a][k];
                            t2[k] += d2 * x6[a][k];
                        }
                    }
                    let cx =
                        [t1[1] * t2[2] - t1[2] * t2[1], t1[2] * t2[0] - t1[0] * t2[2], t1[0] * t2[1] - t1[1] * t2[0]];
                    let n_len = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
                    if n_len <= 1e-30 {
                        continue;
                    }
                    let mut n_hat = [cx[0] / n_len, cx[1] / n_len, cx[2] / n_len];
                    if let Some(po) = p_opp {
                        let d_out =
                            (xq[0] - po[0]) * n_hat[0] + (xq[1] - po[1]) * n_hat[1] + (xq[2] - po[2]) * n_hat[2];
                        if d_out < 0.0 {
                            n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                        }
                    }
                    let mut x_arr = [0.0; D];
                    x_arr.copy_from_slice(&xq[..D]);
                    let mut n_arr = [0.0; D];
                    n_arr.copy_from_slice(&n_hat[..D]);
                    let gval = source.neumann_robin_re(x_arr, n_arr, t, lab as u32)[0];
                    let wq = (0.5 * n_len) * qrule.w[q] * gval;
                    for (a, &d) in dofs.iter().enumerate() {
                        f_vec[d] += wq * shapes_fe[q][a];
                    }
                }
            } else {
                let p0 = mesh.point(fv[0] as usize);
                let p1 = mesh.point(fv[1] as usize);
                let p2 = mesh.point(fv[2] as usize);
                let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                let cx = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
                let n_len = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
                if n_len <= 1e-30 {
                    continue;
                }
                let area = 0.5 * n_len;
                let mut n_hat = [cx[0] / n_len, cx[1] / n_len, cx[2] / n_len];
                if let Some(po) = p_opp {
                    let d_out = (p0[0] - po[0]) * n_hat[0] + (p0[1] - po[1]) * n_hat[1] + (p0[2] - po[2]) * n_hat[2];
                    if d_out < 0.0 {
                        n_hat = [-n_hat[0], -n_hat[1], -n_hat[2]];
                    }
                }
                for (q, lam) in qrule.lam.iter().enumerate() {
                    let xq = [
                        lam[0] * p0[0] + lam[1] * p1[0] + lam[2] * p2[0],
                        lam[0] * p0[1] + lam[1] * p1[1] + lam[2] * p2[1],
                        lam[0] * p0[2] + lam[1] * p1[2] + lam[2] * p2[2],
                    ];
                    let mut x_arr = [0.0; D];
                    x_arr.copy_from_slice(&xq[..D]);
                    let mut n_arr = [0.0; D];
                    n_arr.copy_from_slice(&n_hat[..D]);
                    let gval = source.neumann_robin_re(x_arr, n_arr, t, lab as u32)[0];
                    let wq = area * qrule.w[q] * gval;
                    for (ia, &d) in dofs.iter().enumerate() {
                        f_vec[d] += wq * shapes_fe[q][ia];
                    }
                }
            }
        }
    }
    Ok(f_vec)
}

/// Assembles `int_{Gamma} alpha * phi_i * phi_j dS` including curved 3D quadratic boundary triangles.
pub(crate) fn assemble_boundary_mass_curved_aware(
    mesh: &Discretisation,
    labels: &[i32],
    order: usize,
    alpha: f64,
) -> Csr {
    let n = mesh.n_dofs(order);
    let label_set: BTreeSet<i32> = labels.iter().copied().collect();
    if label_set.is_empty() || alpha == 0.0 {
        return Csr::zeros(n, n);
    }
    if mesh.dim != 3 || !mesh.has_curved() {
        return femcore::fem::assemble_boundary_mass(mesh, labels, order, alpha);
    }
    let qrule = quad_tri_high();
    let shapes_fe: Vec<Vec<f64>> = qrule.lam.iter().map(|lam| shape(lam, order).0).collect();
    let shapes_g2: Vec<(Vec<f64>, Vec<f64>)> = qrule.lam.iter().map(|lam| shape(lam, 2)).collect();
    let nloc = if order == 1 { 3 } else { 6 };
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut ke = vec![0.0; nloc * nloc];
    for f_idx in 0..mesh.n_facets() {
        if !label_set.contains(&mesh.facet_labels[f_idx]) {
            continue;
        }
        let fv = mesh.facet(f_idx);
        let x6 = facet_3d_geom_nodes(mesh, fv);
        ke.fill(0.0);
        for (q, _) in qrule.lam.iter().enumerate() {
            let (_, ref dn2) = shapes_g2[q];
            let mut t1 = [0.0; 3];
            let mut t2 = [0.0; 3];
            for a in 0..6 {
                let d1 = dn2[a * 3 + 1] - dn2[a * 3];
                let d2 = dn2[a * 3 + 2] - dn2[a * 3];
                for k in 0..3 {
                    t1[k] += d1 * x6[a][k];
                    t2[k] += d2 * x6[a][k];
                }
            }
            let cx = [t1[1] * t2[2] - t1[2] * t2[1], t1[2] * t2[0] - t1[0] * t2[2], t1[0] * t2[1] - t1[1] * t2[0]];
            let n_len = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            if n_len <= 1e-30 {
                continue;
            }
            let wq = alpha * (0.5 * n_len) * qrule.w[q];
            let nv = &shapes_fe[q];
            for ia in 0..nloc {
                for ib in 0..nloc {
                    ke[ia * nloc + ib] += wq * nv[ia] * nv[ib];
                }
            }
        }
        mesh.facet_dofs(f_idx, order, &mut dofs);
        coo.push_block(&dofs, &dofs, &ke);
    }
    coo.to_csr()
}
