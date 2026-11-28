// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Topology and shape optimization primitives on unstructured 2D and 3D meshes (ADR 142):
//!
//! * [`HelmholtzFilter`]: unstructured PDE density filter `(-l_f^2 div grad + 1) u = rho` assembled with P1
//!   Lagrange elements and solved via [`DirectSolver`], with exact partition of unity, volume conservation,
//!   and discrete transpose pullback `F^T g` for adjoint chain-rule sensitivities.
//! * [`heaviside_project`], [`simp_interpolate`], [`ramp_interpolate`], [`optical_eps_interpolate`]:
//!   smooth threshold projection, mechanical/thermal SIMP & RAMP interpolations, plasmon-safe complex
//!   refractive-index `(n, kappa)` optical permittivity interpolation, and `beta` continuation schedules.
//! * [`ShapeVelocityField`] and [`compute_shape_velocity`]: fixed-topology CAD boundary velocity extraction
//!   from `geom` parametric evaluations extended smoothly to interior vertices via inverse-volume-stiffened
//!   harmonic Laplacian solves, plus mesh displacement and exact affine simplex shape derivatives.

use crate::fem::{geometry, Geometry};
use crate::linalg::csr::{Coo, RCsr as Csr};
use crate::linalg::driver::DirectSolver;
use crate::linalg::scalar::C64;
use crate::mesh::cells::CellMesh;
use crate::mesh::{Mesh, SimplexMesh};
use crate::section::{wafer_axes, wafer_point};
use geom::complex::FaceGeom;
use geom::eval::{Evaluated, Evaluator, Tier};
use geom::labels;
use geom::model::FeatureList;
use geom::sim::EntityRef;
use std::collections::{BTreeMap, BTreeSet};

// ---------------------------------------------------------------------------
// Unstructured Helmholtz PDE density filter
// ---------------------------------------------------------------------------

/// Unstructured P1 Helmholtz PDE density filter `(-l_f^2 div grad + 1) u = rho` restricted to
/// an active set of design cells (all cells by default, or specific domain labels).
///
/// Forward filtering maps cellwise design variables `rho` (`len == n_design_cells`) to nodal field
/// `u = H^-1 B rho`, where `H = l_f^2 K_P1 + M_P1` (homogeneous Neumann on the design-subdomain
/// boundary) and `B_{i,e} = |Omega_e| / (d + 1)`, followed by cell averaging
/// `(A_avg u)_e = sum_{a=0}^d u_{i_a(e)} / (d + 1)`.
///
/// Consequently:
/// * Constant fields are preserved (`F 1 = 1`) and total integrated volume `sum_e |Omega_e| rho_e`
///   is conserved to machine precision on arbitrary 2D/3D unstructured meshes.
/// * `F` is self-adjoint in the volume-weighted inner product `<x, y>_V = sum_e |Omega_e| x_e y_e`.
/// * [`Self::pullback_gradient`] applies the exact Euclidean transpose `F^T g = B^T H^-1 A_avg^T g`
///   reusing the single cached Cholesky factorisation of `H`.
pub struct HelmholtzFilter {
    dim: usize,
    n_mesh_cells: usize,
    filter_radius: f64,
    design_cells: Vec<usize>,
    cell_volumes: Vec<f64>,
    cell_local_verts: Vec<usize>,
    n_active_verts: usize,
    solver: Option<DirectSolver<f64>>,
}

impl HelmholtzFilter {
    /// Builds and factorises the Helmholtz PDE density filter on `mesh`.
    ///
    /// When `design_domains` is `Some(labels)`, only cells carrying one of `labels` participate in
    /// the design subdomain; when `None`, all mesh cells participate. Setting `filter_radius == 0.0`
    /// selects an exact identity filter without factorising a matrix.
    pub fn new(
        mesh: &SimplexMesh,
        filter_radius: f64,
        design_domains: Option<&[i32]>,
    ) -> Result<Self, String> {
        if mesh.dim != 2 && mesh.dim != 3 {
            return Err(format!("HelmholtzFilter requires a 2D or 3D mesh, got dim = {}", mesh.dim));
        }
        if !filter_radius.is_finite() || filter_radius < 0.0 {
            return Err(format!(
                "HelmholtzFilter filter_radius must be finite and non-negative, got {filter_radius}"
            ));
        }
        let n_mesh_cells = mesh.n_cells();
        let design_cells: Vec<usize> = match design_domains {
            Some(doms) => (0..n_mesh_cells)
                .filter(|&c| doms.contains(&mesh.cell_labels[c]))
                .collect(),
            None => (0..n_mesh_cells).collect(),
        };
        if design_cells.is_empty() {
            return Err("HelmholtzFilter design domain contains zero mesh cells".into());
        }

        let geo = geometry(mesh);
        let d = mesh.dim;
        let nloc = d + 1;
        let mut cell_volumes = Vec::with_capacity(design_cells.len());
        let mut vert_set = BTreeSet::new();
        for &c in &design_cells {
            let v_vol = geo.vol[c];
            if !v_vol.is_finite() || v_vol <= 0.0 {
                return Err(format!("cell {c} has non-positive or non-finite volume {v_vol}"));
            }
            cell_volumes.push(v_vol);
            for &vid in mesh.cell(c) {
                vert_set.insert(vid as usize);
            }
        }

        let active_verts: Vec<usize> = vert_set.into_iter().collect();
        let n_active_verts = active_verts.len();
        let mut global_to_local = vec![usize::MAX; mesh.n_points()];
        for (li, &gv) in active_verts.iter().enumerate() {
            if gv >= global_to_local.len() {
                return Err(format!("vertex {gv} out of bounds for {} points", mesh.n_points()));
            }
            global_to_local[gv] = li;
        }

        let mut cell_local_verts = Vec::with_capacity(design_cells.len() * nloc);
        for &c in &design_cells {
            for &vid in mesh.cell(c) {
                cell_local_verts.push(global_to_local[vid as usize]);
            }
        }

        if filter_radius == 0.0 {
            return Ok(Self {
                dim: d,
                n_mesh_cells,
                filter_radius,
                design_cells,
                cell_volumes,
                cell_local_verts,
                n_active_verts,
                solver: None,
            });
        }

        let lf2 = filter_radius * filter_radius;
        let mass_denom = ((d + 1) * (d + 2)) as f64;
        let mut coo = Coo::with_capacity(
            n_active_verts,
            n_active_verts,
            design_cells.len() * nloc * nloc,
        );
        let mut lv = [0usize; 4];
        let mut he = [0.0f64; 16];
        for (idx, &c) in design_cells.iter().enumerate() {
            let vol = cell_volumes[idx];
            let g_base = c * nloc * d;
            for a in 0..nloc {
                lv[a] = cell_local_verts[idx * nloc + a];
            }
            for a in 0..nloc {
                for b in 0..nloc {
                    let mut grad_dot = 0.0;
                    for k in 0..d {
                        grad_dot += geo.g[g_base + a * d + k] * geo.g[g_base + b * d + k];
                    }
                    let m_ab = vol * (if a == b { 2.0 } else { 1.0 }) / mass_denom;
                    let k_ab = vol * grad_dot;
                    he[a * nloc + b] = lf2 * k_ab + m_ab;
                }
            }
            coo.push_block(&lv[..nloc], &lv[..nloc], &he[..nloc * nloc]);
        }
        let h_mat = coo.to_csr();
        let solver = DirectSolver::new(&h_mat, true)?;

        Ok(Self {
            dim: d,
            n_mesh_cells,
            filter_radius,
            design_cells,
            cell_volumes,
            cell_local_verts,
            n_active_verts,
            solver: Some(solver),
        })
    }

    /// Builds a [`HelmholtzFilter`] directly from an exchange [`Mesh`] via [`Mesh::discretise`].
    pub fn from_exchange_mesh(
        mesh: &Mesh,
        filter_radius: f64,
        design_domains: Option<&[i32]>,
    ) -> Result<Self, String> {
        let simplex = mesh.discretise()?;
        Self::new(&simplex, filter_radius, design_domains)
    }

    /// Filter length scale `l_f` (m).
    pub fn filter_radius(&self) -> f64 {
        self.filter_radius
    }

    /// Indices of the active design cells in `0..mesh.n_cells()`.
    pub fn design_cells(&self) -> &[usize] {
        &self.design_cells
    }

    /// Number of active design cells (`rho.len()`).
    pub fn n_design_cells(&self) -> usize {
        self.design_cells.len()
    }

    /// Volumes (2D areas, 3D volumes) of the active design cells in the order of [`Self::design_cells`].
    pub fn cell_volumes(&self) -> &[f64] {
        &self.cell_volumes
    }

    /// Total measure `sum_e |Omega_e|` of the active design region.
    pub fn total_volume(&self) -> f64 {
        self.cell_volumes.iter().sum()
    }

    /// Applies the forward Helmholtz PDE filter `rho_tilde = F rho` to design-cell densities `rho`.
    pub fn filter_densities(&self, rho: &[f64]) -> Result<Vec<f64>, String> {
        let nc = self.design_cells.len();
        if rho.len() != nc {
            return Err(format!(
                "HelmholtzFilter::filter_densities expected {nc} design densities, got {}",
                rho.len()
            ));
        }
        if let Some(pos) = rho.iter().position(|x| !x.is_finite()) {
            return Err(format!("non-finite design density at index {pos}"));
        }
        let Some(solver) = &self.solver else {
            return Ok(rho.to_vec());
        };

        let nloc = self.dim + 1;
        let inv_nloc = 1.0 / nloc as f64;
        let mut rhs = vec![0.0; self.n_active_verts];
        for idx in 0..nc {
            let contrib = self.cell_volumes[idx] * inv_nloc * rho[idx];
            let base = idx * nloc;
            for a in 0..nloc {
                rhs[self.cell_local_verts[base + a]] += contrib;
            }
        }

        let u_nodal = solver.solve(&rhs);
        let mut out = vec![0.0; nc];
        for idx in 0..nc {
            let base = idx * nloc;
            let mut sum = 0.0;
            for a in 0..nloc {
                sum += u_nodal[self.cell_local_verts[base + a]];
            }
            out[idx] = sum * inv_nloc;
        }
        Ok(out)
    }

    /// Applies the exact discrete transpose `dJ/d(rho) = F^T (dJ/d(rho_tilde))` so that
    /// `<g, F delta_rho> == <F^T g, delta_rho>` holds to machine precision for all `g`, `delta_rho`.
    pub fn pullback_gradient(&self, d_j_d_rho_tilde: &[f64]) -> Result<Vec<f64>, String> {
        let nc = self.design_cells.len();
        if d_j_d_rho_tilde.len() != nc {
            return Err(format!(
                "HelmholtzFilter::pullback_gradient expected {nc} entries, got {}",
                d_j_d_rho_tilde.len()
            ));
        }
        if let Some(pos) = d_j_d_rho_tilde.iter().position(|x| !x.is_finite()) {
            return Err(format!("non-finite sensitivity entry at index {pos}"));
        }
        let Some(solver) = &self.solver else {
            return Ok(d_j_d_rho_tilde.to_vec());
        };

        let nloc = self.dim + 1;
        let inv_nloc = 1.0 / nloc as f64;
        let mut adj_rhs = vec![0.0; self.n_active_verts];
        for idx in 0..nc {
            let contrib = inv_nloc * d_j_d_rho_tilde[idx];
            let base = idx * nloc;
            for a in 0..nloc {
                adj_rhs[self.cell_local_verts[base + a]] += contrib;
            }
        }

        let w_nodal = solver.solve(&adj_rhs);
        let mut grad_rho = vec![0.0; nc];
        for idx in 0..nc {
            let base = idx * nloc;
            let mut sum = 0.0;
            for a in 0..nloc {
                sum += w_nodal[self.cell_local_verts[base + a]];
            }
            grad_rho[idx] = self.cell_volumes[idx] * inv_nloc * sum;
        }
        Ok(grad_rho)
    }

    /// Scatters a design-cell vector (`len == n_design_cells`) onto all `n_mesh_cells` cells,
    /// filling non-design cells with `default_val`.
    pub fn scatter_to_mesh(&self, design_vals: &[f64], default_val: f64) -> Result<Vec<f64>, String> {
        if design_vals.len() != self.design_cells.len() {
            return Err(format!(
                "scatter_to_mesh expected {} values, got {}",
                self.design_cells.len(),
                design_vals.len()
            ));
        }
        let mut out = vec![default_val; self.n_mesh_cells];
        for (&c, &v) in self.design_cells.iter().zip(design_vals) {
            out[c] = v;
        }
        Ok(out)
    }
}

// ---------------------------------------------------------------------------
// Heaviside projection, SIMP / RAMP / plasmon-safe optical interpolation
// ---------------------------------------------------------------------------

/// Smooth threshold Heaviside projection (Wang et al. 2011) at a single point:
/// returns `(rho_bar, d(rho_bar) / d(rho_tilde))`.
///
/// Smoothly reduces to the identity `(rho_tilde, 1.0)` as `beta -> 0^+`.
#[inline]
pub fn heaviside_project_scalar(rho_tilde: f64, beta: f64, eta: f64) -> (f64, f64) {
    let eta_c = eta.clamp(1e-6, 1.0 - 1e-6);
    if beta <= 1e-8 {
        return (rho_tilde, 1.0);
    }
    let th_eta = (beta * eta_c).tanh();
    let denom = th_eta + (beta * (1.0 - eta_c)).tanh();
    if denom.abs() < 1e-30 {
        return (rho_tilde, 1.0);
    }
    let arg = beta * (rho_tilde - eta_c);
    let rho_bar = (th_eta + arg.tanh()) / denom;
    let sech2 = if arg.abs() > 50.0 {
        0.0
    } else {
        let ch = arg.cosh();
        1.0 / (ch * ch)
    };
    let deriv = beta * sech2 / denom;
    (rho_bar, deriv)
}

/// Applies [`heaviside_project_scalar`] elementwise, returning `(rho_bar, d_rho_bar_d_rho_tilde)`.
pub fn heaviside_project(rho_tilde: &[f64], beta: f64, eta: f64) -> (Vec<f64>, Vec<f64>) {
    let mut bar = Vec::with_capacity(rho_tilde.len());
    let mut deriv = Vec::with_capacity(rho_tilde.len());
    for &r in rho_tilde {
        let (b, d) = heaviside_project_scalar(r, beta, eta);
        bar.push(b);
        deriv.push(d);
    }
    (bar, deriv)
}

/// Generates a geometric `beta`-continuation schedule `(beta_stage, iters_stage)` from `beta_init`
/// to `beta_max` across `stages` stages.
pub fn continuation_schedule(
    beta_init: f64,
    beta_max: f64,
    stages: usize,
    iters_per_stage: usize,
) -> Vec<(f64, usize)> {
    let n = stages.max(1);
    let it = iters_per_stage.max(1);
    let b0 = beta_init.max(1e-3);
    let b1 = beta_max.max(b0);
    if n == 1 {
        return vec![(b1, it)];
    }
    let ratio = (b1 / b0).powf(1.0 / (n - 1) as f64);
    (0..n).map(|k| (b0 * ratio.powi(k as i32), it)).collect()
}

/// Standard SIMP interpolation `val(rho) = val_min + rho^p * (val_max - val_min)` and analytic
/// derivative `d(val) / d(rho)`.
#[inline]
pub fn simp_interpolate(rho_bar: f64, val_min: f64, val_max: f64, penalty: f64) -> (f64, f64) {
    let p = penalty.max(1.0);
    let dv = val_max - val_min;
    if (p - 1.0).abs() < 1e-14 {
        return (val_min + rho_bar * dv, dv);
    }
    let r = rho_bar.max(0.0);
    let w = r.powf(p);
    let dw = if r <= 0.0 { 0.0 } else { p * r.powf(p - 1.0) };
    (val_min + w * dv, dw * dv)
}

/// Rational Approximation of Material Properties (RAMP, Stolpe & Svanberg 2001):
/// `val(rho) = val_min + (rho / (1 + q (1 - rho))) * (val_max - val_min)` and `d(val)/d(rho)`.
/// Has strictly positive slope at `rho = 0` for all `q >= 0`.
#[inline]
pub fn ramp_interpolate(rho_bar: f64, val_min: f64, val_max: f64, q: f64) -> (f64, f64) {
    let q_pos = q.max(0.0);
    let denom = 1.0 + q_pos * (1.0 - rho_bar);
    let dv = val_max - val_min;
    let w = rho_bar / denom;
    let dw = (1.0 + q_pos) / (denom * denom);
    (val_min + w * dv, dw * dv)
}

/// Interpolation law for complex relative permittivity `eps_r(rho_bar)` in optical / RF topology
/// optimization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum OpticalInterpMode {
    /// Automatically selects `DielectricSimp` when both endpoints have `Re(eps) > 0` with small
    /// extinction, and switches to plasmon-safe `RefractiveIndex` whenever `Re(eps) <= 0` (metals)
    /// or `|Im(eps)| > 1e-3 * |eps|`.
    Auto { penalty: f64, damping_floor: f64 },
    /// Power-law dielectric interpolation `eps_void + rho^p * (eps_solid - eps_void)`.
    DielectricSimp { penalty: f64 },
    /// RAMP dielectric interpolation `eps_void + (rho / (1 + q*(1 - rho))) * (eps_solid - eps_void)`.
    DielectricRamp { q: f64 },
    /// Refractive-index `(n, kappa)` interpolation `eps(rho) = N(rho)^2` with optional intermediate
    /// attenuation `damping_floor * 4 * rho * (1 - rho)` so `Re(eps)` never crosses zero undamped
    /// (preventing artificial localized plasmon resonances at grey densities).
    RefractiveIndex { damping_floor: f64 },
}

#[inline]
fn passive_sqrt_eps(eps: C64, loss_im_sign: f64) -> C64 {
    let r = (eps.re * eps.re + eps.im * eps.im).sqrt();
    if r == 0.0 {
        return C64::ZERO;
    }
    let n_re = (0.5 * (r + eps.re).max(0.0)).sqrt();
    let mut n_im = (0.5 * (r - eps.re).max(0.0)).sqrt();
    if eps.im < 0.0 || (eps.im == 0.0 && loss_im_sign < 0.0) {
        n_im = -n_im;
    }
    C64::new(n_re, n_im)
}

/// Evaluates `(eps_r(rho_bar), d(eps_r) / d(rho_bar))` under [`OpticalInterpMode`].
#[inline]
pub fn optical_eps_interpolate(
    rho_bar: f64,
    eps_void: C64,
    eps_solid: C64,
    mode: OpticalInterpMode,
) -> (C64, C64) {
    let eval_refractive = |damping_floor: f64| -> (C64, C64) {
        let loss_sign = if eps_void.im > 0.0 || eps_solid.im > 0.0 { 1.0 } else { -1.0 };
        let n0 = passive_sqrt_eps(eps_void, loss_sign);
        let n1 = passive_sqrt_eps(eps_solid, loss_sign);
        let dn_lin = n1 - n0;
        let damp = loss_sign * damping_floor.max(0.0) * 4.0 * rho_bar * (1.0 - rho_bar);
        let ddamp = loss_sign * damping_floor.max(0.0) * 4.0 * (1.0 - 2.0 * rho_bar);
        let n_eff = C64::new(n0.re + rho_bar * dn_lin.re, n0.im + rho_bar * dn_lin.im + damp);
        let dn_eff = C64::new(dn_lin.re, dn_lin.im + ddamp);
        let eps = n_eff * n_eff;
        let deps = C64::new(2.0, 0.0) * n_eff * dn_eff;
        (eps, deps)
    };
    let eval_simp = |penalty: f64| -> (C64, C64) {
        let (re, dre) = simp_interpolate(rho_bar, eps_void.re, eps_solid.re, penalty);
        let (im, dim) = simp_interpolate(rho_bar, eps_void.im, eps_solid.im, penalty);
        (C64::new(re, im), C64::new(dre, dim))
    };

    match mode {
        OpticalInterpMode::DielectricSimp { penalty } => eval_simp(penalty),
        OpticalInterpMode::DielectricRamp { q } => {
            let (re, dre) = ramp_interpolate(rho_bar, eps_void.re, eps_solid.re, q);
            let (im, dim) = ramp_interpolate(rho_bar, eps_void.im, eps_solid.im, q);
            (C64::new(re, im), C64::new(dre, dim))
        }
        OpticalInterpMode::RefractiveIndex { damping_floor } => eval_refractive(damping_floor),
        OpticalInterpMode::Auto { penalty, damping_floor } => {
            let v_abs = (eps_void.re * eps_void.re + eps_void.im * eps_void.im).sqrt();
            let s_abs = (eps_solid.re * eps_solid.re + eps_solid.im * eps_solid.im).sqrt();
            let is_metallic_or_lossy = eps_void.re <= 0.0
                || eps_solid.re <= 0.0
                || eps_void.im.abs() > 1e-3 * v_abs.max(1.0)
                || eps_solid.im.abs() > 1e-3 * s_abs.max(1.0);
            if is_metallic_or_lossy {
                let floor = if damping_floor > 0.0 {
                    damping_floor
                } else if eps_void.re * eps_solid.re <= 0.0 && eps_void.im == 0.0 && eps_solid.im == 0.0 {
                    0.25
                } else {
                    0.0
                };
                eval_refractive(floor)
            } else {
                eval_simp(penalty)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Parametric CAD / Harmonic Mesh-Velocity Shape Sensitivity
// ---------------------------------------------------------------------------

/// Smooth fixed-topology mesh design-velocity field `V_mesh = dx_mesh / dp` in mesh coordinates
/// (`[vx, vz, 0]` in 2D cross-sections, `[vx, vy, vz]` in 3D), with curved mid-edge velocities
/// for quadratic isoparametric edges (`ADR 82`).
#[derive(Clone, Debug, PartialEq)]
pub struct ShapeVelocityField {
    pub dim: usize,
    /// Design velocity per mesh vertex (`len == n_points`), with leading `dim` entries in mesh axes.
    pub vertex_vel: Vec<[f64; 3]>,
    /// Design velocity of curved edge midpoints `(sorted_vertex_pair, vel_in_mesh_axes)`.
    pub edge_mid_vel: Vec<([u32; 2], [f64; 3])>,
}

impl ShapeVelocityField {
    /// Supremum norm `max_i ||v_i||_2` across all vertices.
    pub fn max_norm(&self) -> f64 {
        self.vertex_vel
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0f64, f64::max)
    }

    /// Chooses an optimal fixed-topology step `h` around `base_rel * char_len / ||V||_inf` for
    /// directional operator assembly derivatives `(K(x + h V) - K(x - h V)) / (2 h)`.
    pub fn optimal_step(&self, char_len: f64, base_rel: f64) -> f64 {
        let vmax = self.max_norm();
        if vmax <= 1e-30 {
            base_rel.max(1e-8)
        } else {
            (base_rel * char_len.max(1e-12) / vmax).max(1e-14)
        }
    }
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn wafer_vel_to_mesh(v_wafer: [f64; 3], dim: usize) -> [f64; 3] {
    let axes = wafer_axes(dim);
    let mut out = [0.0; 3];
    for (m_ax, &w_ax) in axes.iter().enumerate() {
        out[m_ax] = v_wafer[w_ax];
    }
    out
}

fn mesh_vel_to_wafer(v_mesh: [f64; 3], dim: usize) -> [f64; 3] {
    let axes = wafer_axes(dim);
    let mut out = [0.0; 3];
    for (m_ax, &w_ax) in axes.iter().enumerate() {
        out[w_ax] = v_mesh[m_ax];
    }
    out
}

/// Computes a smooth fixed-topology mesh velocity field `V_mesh = dx_mesh / dp` with respect to
/// named geometry parameter `param_name` in `list`:
///
/// 1. Perturbs `param_name` by `+-dp` around its nominal value in `eval0` (`Tier::Complex`) and
///    verifies CAD topological invariance (`vertices`, `edges`, `edges3`, `faces` counts preserved).
/// 2. Extracts exact CAD boundary velocities at constrained vertices (0D CAD vertices, 1D curves
///    and 3D crease/loop edges `edges3`, 2D curved/extruded lateral faces, and normal components
///    on 3D planar/cap faces so cap meshes slide tangentially without element squeezing).
/// 3. Solves an inverse-volume-stiffened harmonic Laplacian extension `-div(mu_e grad V_alpha) = 0`
///    (`mu_e = 1 / |Omega_e|`) per mesh coordinate via [`DirectSolver`] to propagate boundary
///    motion smoothly to interior vertices.
/// 4. Projects curved mid-edge geometry nodes (`CellMesh::edge_nodes`) onto the perturbed CAD
///    entities at `p0 +- dp`.
pub fn compute_shape_velocity(
    list: &FeatureList,
    eval0: &Evaluated,
    mesh: &Mesh,
    param_name: &str,
    step_rel: Option<f64>,
) -> Result<ShapeVelocityField, String> {
    let simplex = mesh.discretise()?;
    let dim = simplex.dim;
    let Some(p0) = eval0.scene().env.get(param_name) else {
        return Err(format!("geometry parameter '{param_name}' was not found in model parameters"));
    };
    if !p0.is_finite() {
        return Err(format!("geometry parameter '{param_name}' has non-finite value {p0}"));
    }

    let diag = eval0.bounds().diagonal().max(1e-9);
    let rel = step_rel.unwrap_or(1e-6).clamp(1e-9, 1e-3);
    // Keep perturbation well below minimum element altitude while well above f64 roundoff.
    let geo = geometry(&simplex);
    let min_vol = geo
        .vol
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
        .max(1e-36);
    let h_min = min_vol.powf(1.0 / dim as f64);
    let scale = p0.abs().max(h_min).max(diag * 1e-3).max(1e-9);
    let dp = (rel * scale).min(0.05 * h_min.max(1e-12)).max(1e-12 * scale.max(1.0));

    let mut ev_cache = Evaluator::new();
    let list_plus = list.with(&[(param_name, p0 + dp)]);
    let list_minus = list.with(&[(param_name, p0 - dp)]);
    let ev_plus = ev_cache
        .eval(&list_plus, Tier::Complex)
        .map_err(|e| format!("perturbed (+dp) geometry evaluation for '{param_name}' failed: {e}"))?;
    let ev_minus = ev_cache
        .eval(&list_minus, Tier::Complex)
        .map_err(|e| format!("perturbed (-dp) geometry evaluation for '{param_name}' failed: {e}"))?;

    if ev_plus.vertices().len() != eval0.vertices().len()
        || ev_minus.vertices().len() != eval0.vertices().len()
        || ev_plus.edges().len() != eval0.edges().len()
        || ev_minus.edges().len() != eval0.edges().len()
        || ev_plus.edges3().len() != eval0.edges3().len()
        || ev_minus.edges3().len() != eval0.edges3().len()
        || ev_plus.faces().len() != eval0.faces().len()
        || ev_minus.faces().len() != eval0.faces().len()
    {
        return Err(format!(
            "perturbing parameter '{param_name}' by +- {dp:e} changed CAD topology; shape sensitivity requires fixed topology"
        ));
    }

    let label_map = labels::compile(eval0)
        .map_err(|e| format!("failed to compile CAD labels for shape velocity: {e}"))?;
    let n_pts = simplex.n_points();
    let mean_vol = (geo.vol.iter().sum::<f64>() / simplex.n_cells().max(1) as f64).max(1e-36);
    let h_char = mean_vol.powf(1.0 / dim as f64);
    let tol = (1e-4 * h_char).max(1e-5 * diag).max(1e-11);
    let inv_2dp = 1.0 / (2.0 * dp);

    // Collect incident facet labels per vertex and detect interface/exterior vertices.
    let mut vert_facet_labels: Vec<BTreeSet<i32>> = vec![BTreeSet::new(); n_pts];
    for f in 0..simplex.n_facets() {
        let lbl = simplex.facet_labels[f];
        for &vid in simplex.facet(f) {
            vert_facet_labels[vid as usize].insert(lbl);
        }
    }
    let mut vert_domains: Vec<Option<i32>> = vec![None; n_pts];
    let mut vert_multi_domain = vec![false; n_pts];
    for c in 0..simplex.n_cells() {
        let dom = simplex.cell_labels[c];
        for &vid in simplex.cell(c) {
            let vi = vid as usize;
            match vert_domains[vi] {
                None => vert_domains[vi] = Some(dom),
                Some(d0) if d0 != dom => vert_multi_domain[vi] = true,
                _ => {}
            }
        }
    }
    // Also detect unlabelled exterior boundary vertices.
    let mut exterior_verts = vec![false; n_pts];
    {
        use std::collections::HashMap;
        let nloc = dim + 1;
        let mut fcount: HashMap<[u32; 3], (u32, [u32; 3])> = HashMap::with_capacity(simplex.n_cells() * 2);
        for c in 0..simplex.n_cells() {
            let cv = simplex.cell(c);
            for skip in 0..nloc {
                let mut key = [u32::MAX; 3];
                let mut raw = [u32::MAX; 3];
                let mut p = 0;
                for i in 0..nloc {
                    if i != skip {
                        key[p] = cv[i];
                        raw[p] = cv[i];
                        p += 1;
                    }
                }
                key[..dim].sort_unstable();
                let entry = fcount.entry(key).or_insert((0, raw));
                entry.0 += 1;
            }
        }
        for (_, (cnt, raw)) in fcount {
            if cnt == 1 {
                for &vid in &raw[..dim] {
                    exterior_verts[vid as usize] = true;
                }
            }
        }
    }

    // Per mesh axis `alpha in 0..dim`, record Dirichlet boundary values `fixed_comp[alpha]: BTreeMap<usize, f64>`.
    let mut fixed_comp: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); dim];
    let oracle0 = eval0.oracle();
    let n_cad_verts = eval0.vertices().len();
    let n_cad_edges2 = eval0.edges().len();
    let n_cad_edges3 = eval0.edges3().len();
    let n_cad_faces = eval0.faces().len();

    for vi in 0..n_pts {
        let is_cand = !vert_facet_labels[vi].is_empty() || vert_multi_domain[vi] || exterior_verts[vi];
        if !is_cand {
            continue;
        }
        let pt0 = wafer_point(simplex.point(vi), dim);

        // 1. 0D CAD Vertex match (2D complex vertices are in world frame; 3D corners are endpoints of `edges3`)
        if dim == 2 {
            let mut best_v: Option<(usize, f64)> = None;
            for cvi in 0..n_cad_verts {
                let d = dist3(pt0, eval0.vertices()[cvi].p);
                if d <= tol && best_v.is_none_or(|b| d < b.1) {
                    best_v = Some((cvi, d));
                }
            }
            if let Some((cvi, _)) = best_v {
                let pp = ev_plus.vertices()[cvi].p;
                let pm = ev_minus.vertices()[cvi].p;
                let vw = [
                    (pp[0] - pm[0]) * inv_2dp,
                    (pp[1] - pm[1]) * inv_2dp,
                    (pp[2] - pm[2]) * inv_2dp,
                ];
                let vm = wafer_vel_to_mesh(vw, dim);
                for a in 0..dim {
                    fixed_comp[a].insert(vi, vm[a]);
                }
                continue;
            }
        }

        if dim == 2 {
            // 2D CAD Edge match (prefer edges belonging to incident facet labels within surface tolerance)
            let mut cand_edges: Vec<u32> = vert_facet_labels[vi]
                .iter()
                .flat_map(|&lbl| label_map.edges_of_label(lbl))
                .filter(|&ei| (ei as usize) < n_cad_edges2)
                .collect();
            cand_edges.sort_unstable();
            cand_edges.dedup();
            let find_best_edge = |list_e: &[u32], max_d: f64| -> Option<(u32, f64)> {
                let mut best: Option<(u32, f64)> = None;
                for &ei in list_e {
                    if let Some(proj) = oracle0.project(pt0, EntityRef::Edge(ei)) {
                        let d = dist3(pt0, proj);
                        if d <= max_d && best.is_none_or(|b| d < b.1) {
                            best = Some((ei, d));
                        }
                    }
                }
                best
            };
            let matched_edge = find_best_edge(&cand_edges, 10.0 * tol)
                .or_else(|| {
                    let all_e: Vec<u32> = (0..n_cad_edges2 as u32).collect();
                    find_best_edge(&all_e, tol)
                })
                .map(|(ei, _)| ei);
            if let Some(ei) = matched_edge {
                if let Some(fr) = eval0.fractions_of(EntityRef::Edge(ei), pt0) {
                    if let (Some(pp), Some(pm)) = (
                        ev_plus.position_at(EntityRef::Edge(ei), fr),
                        ev_minus.position_at(EntityRef::Edge(ei), fr),
                    ) {
                        let vw = [
                            (pp[0] - pm[0]) * inv_2dp,
                            (pp[1] - pm[1]) * inv_2dp,
                            (pp[2] - pm[2]) * inv_2dp,
                        ];
                        let vm = wafer_vel_to_mesh(vw, dim);
                        for a in 0..dim {
                            fixed_comp[a].insert(vi, vm[a]);
                        }
                    }
                }
            }
        } else {
            // 3D: first check 1D CAD crease/loop edges `edges3` within surface coincidence tolerance `tol`
            // (in cut-cell 3D meshes a bevel triangle near a crease has vertices on distinct faces, so
            // surface coincidence `d <= tol` ensures only vertices truly on the entity are constrained).
            let mut best_e3: Option<(u32, f64)> = None;
            for ei in 0..n_cad_edges3 as u32 {
                if let Some(proj) = oracle0.project(pt0, EntityRef::Edge(ei)) {
                    let d = dist3(pt0, proj);
                    if d <= tol && best_e3.is_none_or(|b| d < b.1) {
                        best_e3 = Some((ei, d));
                    }
                }
            }
            if let Some((ei, _)) = best_e3 {
                if let Some(fr) = eval0.fractions_of(EntityRef::Edge(ei), pt0) {
                    if let (Some(pp), Some(pm)) = (
                        ev_plus.position_at(EntityRef::Edge(ei), fr),
                        ev_minus.position_at(EntityRef::Edge(ei), fr),
                    ) {
                        let vw = [
                            (pp[0] - pm[0]) * inv_2dp,
                            (pp[1] - pm[1]) * inv_2dp,
                            (pp[2] - pm[2]) * inv_2dp,
                        ];
                        for a in 0..3 {
                            fixed_comp[a].insert(vi, vw[a]);
                        }
                        continue;
                    }
                }
            }

            // 3D: check 2D CAD Faces within coincidence tolerance `tol` (preferring incident facet labels).
            let mut cand_faces: Vec<u32> = vert_facet_labels[vi]
                .iter()
                .flat_map(|&lbl| label_map.edges_of_label(lbl))
                .filter(|&fi| (fi as usize) < n_cad_faces)
                .collect();
            cand_faces.sort_unstable();
            cand_faces.dedup();
            let collect_matching_faces = |list_f: &[u32], max_d: f64| -> Vec<u32> {
                let mut hits = Vec::new();
                for &fi in list_f {
                    if let Some(proj) = oracle0.project(pt0, EntityRef::Face(fi)) {
                        if dist3(pt0, proj) <= max_d {
                            hits.push(fi);
                        }
                    }
                }
                hits
            };
            let mut matched_faces = collect_matching_faces(&cand_faces, tol);
            if matched_faces.is_empty() {
                let all_f: Vec<u32> = (0..n_cad_faces as u32).collect();
                matched_faces = collect_matching_faces(&all_f, tol);
            }

            // Process planar/cap faces first (constraining their normal components), then lateral/quadric faces.
            matched_faces.sort_by_key(|&fi| {
                matches!(
                    eval0.faces()[fi as usize].geom,
                    FaceGeom::Lateral { .. } | FaceGeom::Quadric(_)
                )
            });

            for fi in matched_faces {
                let face_geom = &eval0.faces()[fi as usize].geom;
                let fr_opt = match face_geom {
                    FaceGeom::Lateral { edge, z0, z1 } => {
                        let cx0 = eval0.complex();
                        cx0.stack.as_ref().and_then(|st| {
                            let ed = st.overlay.edges.get(*edge as usize)?;
                            let q0 = cx0.iv(pt0);
                            let (t, _) = st.overlay.edge_closest(*edge, [q0[0], q0[1]]);
                            let u = if (ed.t1 - ed.t0).abs() > 1e-30 {
                                ((t - ed.t0) / (ed.t1 - ed.t0)).clamp(0.0, 1.0)
                            } else {
                                0.0
                            };
                            let v = if (z1 - z0).abs() > 1e-30 {
                                ((q0[2] - z0) / (z1 - z0)).clamp(0.0, 1.0)
                            } else {
                                0.0
                            };
                            Some(geom::query::Fractions { u, v })
                        })
                    }
                    _ => eval0.fractions_of(EntityRef::Face(fi), pt0),
                };
                let Some(fr) = fr_opt else {
                    continue;
                };
                let eval_face_exact = |ev_cur: &Evaluated| -> Option<[f64; 3]> {
                    let cx_cur = ev_cur.complex();
                    match &cx_cur.faces.get(fi as usize)?.geom {
                        FaceGeom::Lateral { edge, z0, z1 } => {
                            let ov = &cx_cur.stack.as_ref()?.overlay;
                            let ed = ov.edges.get(*edge as usize)?;
                            let u = fr.u.clamp(0.0, 1.0);
                            let v = fr.v.clamp(0.0, 1.0);
                            let q = ov.edge_point(*edge, ed.t0 + u * (ed.t1 - ed.t0));
                            Some(cx_cur.wv([q[0], q[1], z0 + v * (z1 - z0)]))
                        }
                        _ => ev_cur.position_at(EntityRef::Face(fi), fr),
                    }
                };
                let (Some(pp), Some(pm)) = (eval_face_exact(&ev_plus), eval_face_exact(&ev_minus))
                else {
                    continue;
                };
                let vw = [
                    (pp[0] - pm[0]) * inv_2dp,
                    (pp[1] - pm[1]) * inv_2dp,
                    (pp[2] - pm[2]) * inv_2dp,
                ];
                match face_geom {
                    FaceGeom::Lateral { .. } | FaceGeom::Quadric(_) => {
                        for a in 0..3 {
                            fixed_comp[a].insert(vi, vw[a]);
                        }
                    }
                    FaceGeom::Cap { .. } | FaceGeom::Planar { .. } => {
                        let normals = eval0.facet_normals(fi);
                        if let Some(&n) = normals.first() {
                            let axis_aligned = (0..3).find(|&a| n[a].abs() > 0.9999);
                            if let Some(ax) = axis_aligned {
                                fixed_comp[ax].entry(vi).or_insert(vw[ax]);
                            } else {
                                for a in 0..3 {
                                    fixed_comp[a].entry(vi).or_insert(vw[a]);
                                }
                            }
                        } else {
                            for a in 0..3 {
                                fixed_comp[a].entry(vi).or_insert(vw[a]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Assemble inverse-volume-stiffened P1 harmonic Laplacian `K_harm` (`-div(mu_e grad V_alpha) = 0`),
    // protecting small cut boundary cells from inversion during finite shape steps.
    let stiff_coeff: Vec<f64> = geo
        .vol
        .iter()
        .map(|&v| (mean_vol / v.max(1e-36)).clamp(1e-3, 1e4))
        .collect();
    let k_harm = crate::fem::assemble_laplace(&simplex, 1, &stiff_coeff);
    let zero_rhs = vec![0.0; n_pts];

    let mut vertex_vel = vec![[0.0f64; 3]; n_pts];
    for a in 0..dim {
        if fixed_comp[a].is_empty() || fixed_comp[a].values().all(|&v| v.abs() == 0.0) {
            continue;
        }
        let fixed_idx: Vec<usize> = fixed_comp[a].keys().copied().collect();
        let fixed_vals: Vec<f64> = fixed_comp[a].values().copied().collect();
        // Add tiny mass regularization only if a disconnected component had zero Dirichlet constraints.
        let mut is_fixed = vec![false; n_pts];
        for &i in &fixed_idx {
            is_fixed[i] = true;
        }
        let free: Vec<usize> = (0..n_pts).filter(|&i| !is_fixed[i]).collect();
        for (&i, &val) in fixed_idx.iter().zip(&fixed_vals) {
            vertex_vel[i][a] = val;
        }
        if !free.is_empty() {
            let mut kff = k_harm.submatrix(&free, &free);
            // Tiny diagonal safeguard `1e-14 * trace_mean` in case an isolated floating island has no fixed nodes.
            let diag_mean = (kff.diagonal().iter().sum::<f64>() / free.len() as f64).abs().max(1.0);
            let reg = Csr::axpby(1.0, &kff, 1e-14 * diag_mean, &Csr::identity(free.len()));
            kff = reg;
            let kfx = k_harm.submatrix(&free, &fixed_idx);
            let kx = kfx.matvec(&fixed_vals);
            let rhs: Vec<f64> = free.iter().enumerate().map(|(k, &gi)| zero_rhs[gi] - kx[k]).collect();
            let sol = DirectSolver::new(&kff, true)?.solve(&rhs);
            for (k, &gi) in free.iter().enumerate() {
                if !sol[k].is_finite() {
                    return Err(format!("non-finite harmonic mesh velocity at vertex {gi}, axis {a}"));
                }
                vertex_vel[gi][a] = sol[k];
            }
        }
    }

    // Compute curved mid-edge velocities when `mesh` carries quadratic edge nodes (`ADR 82`).
    let existing_edge_nodes: Vec<([u32; 2], [f64; 3])> = match mesh {
        Mesh::D2(cm) => cm
            .edge_nodes
            .iter()
            .map(|&(pr, p)| (pr, [p[0], p[1], 0.0]))
            .collect(),
        Mesh::D3(cm) => cm.edge_nodes.clone(),
    };

    let mut edge_mid_vel = Vec::with_capacity(existing_edge_nodes.len());
    if !existing_edge_nodes.is_empty() {
        let oracle_plus = ev_plus.oracle();
        let oracle_minus = ev_minus.oracle();
        for (pair, mid_m) in existing_edge_nodes {
            let (va, vb) = (pair[0] as usize, pair[1] as usize);
            let mid0 = wafer_point(&mid_m[..dim], dim);
            let pa0 = wafer_point(simplex.point(va), dim);
            let pb0 = wafer_point(simplex.point(vb), dim);
            let vwa = mesh_vel_to_wafer(vertex_vel[va], dim);
            let vwb = mesh_vel_to_wafer(vertex_vel[vb], dim);
            let chord_mid_plus = [
                0.5 * ((pa0[0] + dp * vwa[0]) + (pb0[0] + dp * vwb[0])),
                0.5 * ((pa0[1] + dp * vwa[1]) + (pb0[1] + dp * vwb[1])),
                0.5 * ((pa0[2] + dp * vwa[2]) + (pb0[2] + dp * vwb[2])),
            ];
            let chord_mid_minus = [
                0.5 * ((pa0[0] - dp * vwa[0]) + (pb0[0] - dp * vwb[0])),
                0.5 * ((pa0[1] - dp * vwa[1]) + (pb0[1] - dp * vwb[1])),
                0.5 * ((pa0[2] - dp * vwa[2]) + (pb0[2] - dp * vwb[2])),
            ];

            let mut v_mid_mesh = [
                0.5 * (vertex_vel[va][0] + vertex_vel[vb][0]),
                0.5 * (vertex_vel[va][1] + vertex_vel[vb][1]),
                0.5 * (vertex_vel[va][2] + vertex_vel[vb][2]),
            ];

            if dim == 2 {
                let mut best_e: Option<(u32, f64)> = None;
                for ei in 0..n_cad_edges2 as u32 {
                    if let Some(proj) = oracle0.project(mid0, EntityRef::Edge(ei)) {
                        let d = dist3(mid0, proj);
                        if d <= tol && best_e.is_none_or(|b| d < b.1) {
                            best_e = Some((ei, d));
                        }
                    }
                }
                if let Some((ei, _)) = best_e {
                    if let (Some(mp), Some(mm)) = (
                        oracle_plus.project(chord_mid_plus, EntityRef::Edge(ei)),
                        oracle_minus.project(chord_mid_minus, EntityRef::Edge(ei)),
                    ) {
                        let vw = [
                            (mp[0] - mm[0]) * inv_2dp,
                            (mp[1] - mm[1]) * inv_2dp,
                            (mp[2] - mm[2]) * inv_2dp,
                        ];
                        v_mid_mesh = wafer_vel_to_mesh(vw, dim);
                    }
                }
            } else {
                let mut projected = false;
                let mut best_e3: Option<(u32, f64)> = None;
                for ei in 0..n_cad_edges3 as u32 {
                    if let Some(proj) = oracle0.project(mid0, EntityRef::Edge(ei)) {
                        let d = dist3(mid0, proj);
                        if d <= tol && best_e3.is_none_or(|b| d < b.1) {
                            best_e3 = Some((ei, d));
                        }
                    }
                }
                if let Some((ei, _)) = best_e3 {
                    if let (Some(mp), Some(mm)) = (
                        oracle_plus.project(chord_mid_plus, EntityRef::Edge(ei)),
                        oracle_minus.project(chord_mid_minus, EntityRef::Edge(ei)),
                    ) {
                        v_mid_mesh = [
                            (mp[0] - mm[0]) * inv_2dp,
                            (mp[1] - mm[1]) * inv_2dp,
                            (mp[2] - mm[2]) * inv_2dp,
                        ];
                        projected = true;
                    }
                }
                if !projected {
                    let mut best_f: Option<(u32, f64)> = None;
                    for fi in 0..n_cad_faces as u32 {
                        if matches!(
                            eval0.faces()[fi as usize].geom,
                            FaceGeom::Lateral { .. } | FaceGeom::Quadric(_)
                        ) {
                            if let Some(proj) = oracle0.project(mid0, EntityRef::Face(fi)) {
                                let d = dist3(mid0, proj);
                                if d <= tol && best_f.is_none_or(|b| d < b.1) {
                                    best_f = Some((fi, d));
                                }
                            }
                        }
                    }
                    if let Some((fi, _)) = best_f {
                        if let (Some(mp), Some(mm)) = (
                            oracle_plus.project(chord_mid_plus, EntityRef::Face(fi)),
                            oracle_minus.project(chord_mid_minus, EntityRef::Face(fi)),
                        ) {
                            v_mid_mesh = [
                                (mp[0] - mm[0]) * inv_2dp,
                                (mp[1] - mm[1]) * inv_2dp,
                                (mp[2] - mm[2]) * inv_2dp,
                            ];
                        }
                    }
                }
            }
            edge_mid_vel.push((pair, v_mid_mesh));
        }
    }

    Ok(ShapeVelocityField { dim, vertex_vel, edge_mid_vel })
}

fn displace_cell_mesh<const D: usize>(
    cm: &CellMesh<D>,
    step: f64,
    vel: &ShapeVelocityField,
) -> Result<CellMesh<D>, String> {
    if vel.dim != D {
        return Err(format!("ShapeVelocityField dim {} does not match mesh dim {D}", vel.dim));
    }
    if cm.n_points() != vel.vertex_vel.len() {
        return Err(format!(
            "ShapeVelocityField vertex count {} does not match mesh point count {}",
            vel.vertex_vel.len(),
            cm.n_points()
        ));
    }
    let mut out = cm.clone();
    for (i, p) in out.points.iter_mut().enumerate() {
        for a in 0..D {
            p[a] += step * vel.vertex_vel[i][a];
        }
    }
    if !out.edge_nodes.is_empty() {
        for (pair, mid) in &mut out.edge_nodes {
            let v_mid = vel
                .edge_mid_vel
                .binary_search_by_key(pair, |x| x.0)
                .ok()
                .map(|idx| vel.edge_mid_vel[idx].1)
                .unwrap_or_else(|| {
                    let (va, vb) = (pair[0] as usize, pair[1] as usize);
                    [
                        0.5 * (vel.vertex_vel[va][0] + vel.vertex_vel[vb][0]),
                        0.5 * (vel.vertex_vel[va][1] + vel.vertex_vel[vb][1]),
                        0.5 * (vel.vertex_vel[va][2] + vel.vertex_vel[vb][2]),
                    ]
                });
            for a in 0..D {
                mid[a] += step * v_mid[a];
            }
        }
    }
    Ok(out)
}

/// Displaces an exchange [`Mesh`] along `vel` by `x -> x + step * V_mesh` while preserving
/// exact cell, vertex, edge, and facet topology and numbering.
pub fn displace_mesh(mesh: &Mesh, step: f64, vel: &ShapeVelocityField) -> Result<Mesh, String> {
    match mesh {
        Mesh::D2(cm) => Ok(Mesh::D2(displace_cell_mesh(cm, step, vel)?)),
        Mesh::D3(cm) => Ok(Mesh::D3(displace_cell_mesh(cm, step, vel)?)),
    }
}

/// Displaces a [`SimplexMesh`] along `vel` by `x -> x + step * V_mesh` while preserving
/// identical connectivity, edge ordering, and P1/P2/Nedelec DOF numbering.
pub fn displace_simplex_mesh(
    mesh: &SimplexMesh,
    step: f64,
    vel: &ShapeVelocityField,
) -> Result<SimplexMesh, String> {
    let d = mesh.dim;
    if vel.dim != d {
        return Err(format!("ShapeVelocityField dim {} does not match SimplexMesh dim {d}", vel.dim));
    }
    if mesh.n_points() != vel.vertex_vel.len() {
        return Err(format!(
            "ShapeVelocityField vertex count {} does not match SimplexMesh point count {}",
            vel.vertex_vel.len(),
            mesh.n_points()
        ));
    }
    let mut out = mesh.clone();
    for i in 0..mesh.n_points() {
        for a in 0..d {
            out.points[i * d + a] += step * vel.vertex_vel[i][a];
        }
    }
    if mesh.has_curved() {
        let mut nodes = mesh.edge_node_list();
        for (pair, mid) in &mut nodes {
            let v_mid = vel
                .edge_mid_vel
                .binary_search_by_key(pair, |x| x.0)
                .ok()
                .map(|idx| vel.edge_mid_vel[idx].1)
                .unwrap_or_else(|| {
                    let (va, vb) = (pair[0] as usize, pair[1] as usize);
                    [
                        0.5 * (vel.vertex_vel[va][0] + vel.vertex_vel[vb][0]),
                        0.5 * (vel.vertex_vel[va][1] + vel.vertex_vel[vb][1]),
                        0.5 * (vel.vertex_vel[va][2] + vel.vertex_vel[vb][2]),
                    ]
                });
            for a in 0..d {
                mid[a] += step * v_mid[a];
            }
        }
        let qp = mesh.qp_edge_list();
        out.set_edge_nodes(&nodes);
        if !qp.is_empty() {
            out.set_qp_edges(&qp);
        }
    }
    Ok(out)
}

/// Exact analytic directional derivatives `(d_vol, d_g)` of affine simplex cell volumes and
/// barycentric coordinate gradients `grad lambda_b` along `vel`:
///
/// * `d(|Omega_c|) / dp = |Omega_c| * sum_{a=0}^d (grad lambda_a . v_a)`
/// * `d(grad lambda_b) / dp = -sum_{a=0}^d (grad lambda_b . v_a) grad lambda_a`
pub fn simplex_geometry_shape_deriv(
    mesh: &SimplexMesh,
    geo: &Geometry,
    vel: &ShapeVelocityField,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let d = mesh.dim;
    if vel.dim != d || vel.vertex_vel.len() != mesh.n_points() {
        return Err("ShapeVelocityField dimensions do not match SimplexMesh".into());
    }
    let nc = mesh.n_cells();
    let nloc = d + 1;
    let mut d_vol = vec![0.0; nc];
    let mut d_g = vec![0.0; nc * nloc * d];

    for c in 0..nc {
        let vs = mesh.cell(c);
        let base = c * nloc * d;
        let mut div_v = 0.0;
        let mut gv = [[0.0f64; 4]; 4]; // gv[b][a] = grad lambda_b . v_a
        for b in 0..nloc {
            for a in 0..nloc {
                let va = vel.vertex_vel[vs[a] as usize];
                let mut dot = 0.0;
                for k in 0..d {
                    dot += geo.g[base + b * d + k] * va[k];
                }
                gv[b][a] = dot;
            }
            div_v += gv[b][b];
        }
        d_vol[c] = geo.vol[c] * div_v;
        for b in 0..nloc {
            for k in 0..d {
                let mut acc = 0.0;
                for a in 0..nloc {
                    acc -= gv[b][a] * geo.g[base + a * d + k];
                }
                d_g[base + b * d + k] = acc;
            }
        }
    }
    Ok((d_vol, d_g))
}

/// Exact analytic shape derivative `dK / dp` of the scalar stiffness matrix
/// `K_ij = integral coeff_c grad(phi_i) . grad(phi_j) dV` on affine simplices (`order == 1` or `2`).
pub fn assemble_laplace_shape_deriv(
    mesh: &SimplexMesh,
    order: usize,
    coeff: &[f64],
    vel: &ShapeVelocityField,
) -> Result<Csr, String> {
    let d = mesh.dim;
    let re = crate::fem::RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let (d_vol, d_g) = simplex_geometry_shape_deriv(mesh, &geo, vel)?;
    let nloc = re.nloc;
    let nv_cell = d + 1;
    let n = mesh.n_dofs(order);
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut grad = vec![0.0; nloc * d];
    let mut dgrad = vec![0.0; nloc * d];
    let mut dke = vec![0.0; nloc * nloc];

    for c in 0..mesh.n_cells() {
        dke.fill(0.0);
        let base = c * nv_cell * d;
        let vol_c = geo.vol[c];
        let dvol_c = d_vol[c];
        for q in 0..re.w.len() {
            let dn = &re.dn[q];
            for a in 0..nloc {
                for k in 0..d {
                    let mut g_ak = 0.0;
                    let mut dg_ak = 0.0;
                    for b in 0..nv_cell {
                        let dnb = dn[a * nv_cell + b];
                        g_ak += dnb * geo.g[base + b * d + k];
                        dg_ak += dnb * d_g[base + b * d + k];
                    }
                    grad[a * d + k] = g_ak;
                    dgrad[a * d + k] = dg_ak;
                }
            }
            for a in 0..nloc {
                for b in 0..nloc {
                    let mut gg = 0.0;
                    let mut dgg = 0.0;
                    for k in 0..d {
                        gg += grad[a * d + k] * grad[b * d + k];
                        dgg += dgrad[a * d + k] * grad[b * d + k] + grad[a * d + k] * dgrad[b * d + k];
                    }
                    dke[a * nloc + b] += re.w[q] * coeff[c] * (dvol_c * gg + vol_c * dgg);
                }
            }
        }
        mesh.cell_dofs(c, order, &mut dofs);
        coo.push_block(&dofs, &dofs, &dke);
    }
    Ok(coo.to_csr())
}

/// Exact analytic shape derivative `dM / dp` of the scalar mass matrix
/// `M_ij = integral coeff_c phi_i phi_j dV` on affine simplices (`order == 1` or `2`).
pub fn assemble_mass_shape_deriv(
    mesh: &SimplexMesh,
    order: usize,
    coeff: &[f64],
    vel: &ShapeVelocityField,
) -> Result<Csr, String> {
    let d = mesh.dim;
    let re = crate::fem::RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let (d_vol, _) = simplex_geometry_shape_deriv(mesh, &geo, vel)?;
    let nloc = re.nloc;
    let n = mesh.n_dofs(order);
    let mut me = vec![0.0; nloc * nloc];
    for q in 0..re.w.len() {
        for a in 0..nloc {
            for b in 0..nloc {
                me[a * nloc + b] += re.w[q] * re.n[q][a] * re.n[q][b];
            }
        }
    }
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut dme = vec![0.0; nloc * nloc];
    for c in 0..mesh.n_cells() {
        let sc = d_vol[c] * coeff[c];
        for (o, &m) in dme.iter_mut().zip(&me) {
            *o = m * sc;
        }
        mesh.cell_dofs(c, order, &mut dofs);
        coo.push_block(&dofs, &dofs, &dme);
    }
    Ok(coo.to_csr())
}
