// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Eigenfrequencies of elastic bodies with perfectly matched layers (curvilinear complex coordinate
//! stretching across Cartesian, Cylindrical, and Spherical topologies), pre-solve Christoffel
//! backward-wave stabilization, outer PML clamping, signed quality factors `Q = Re f / (2 Im f)`,
//! and native mixed-cell (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) assembly.

use crate::floquet::mixed_nodal_relations;
use crate::lrb::AMPLIFYING_IM_TOL;
use crate::problem::*;
use crate::solid::{d_matrix, dof_layout, setup};
use femcore::fem::*;
use femcore::fem_mixed::{
    assemble_mass_mixed, lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::material::{
    coupled::EPSILON_0, isotropic_stiffness, piezo_e_global, stiffness_global, tensor_global, M3, M6,
};
use femcore::mesh::cells::canonical_quad_face;
use femcore::mesh::mixed::{
    cell_local_edges, cell_local_quad_faces, cell_local_tri_faces, sorted_tri_face, FacetType, MixedDiscretisation,
    QUAD_EDGES, TRI_EDGES,
};
use femcore::meshing::physics_size::{
    backward_wave_mpml_ratio, piezo_backward_wave_mpml_ratio, piezo_wave_speed_range, wave_speed_range,
    BackwardWaveCheck,
};
use femcore::periodic::{congruence, match_vertices, nodal_relations, Reduction};
use femcore::pml::{adiabatic_loss_factor, auto_scale_params, pml_jacobian, PmlJacobian};
use femcore::props::{Pml, PmlCoord, PmlDir, PmlProfile, PmlStabilization};
use nyx_physics::{cell_fields_of, Discretisation, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::time::Instant;

const OUTER_T_TOL: f64 = 1.0e-5;

/// One resolved active PML domain with wafer-frame bounds, optional auto-scaling, and Christoffel stabilization.
#[derive(Clone, Debug)]
pub(crate) struct ResolvedPmlDomain {
    pub spec: Pml,
    pub bounds: ([f64; 3], [f64; 3]),
}

/// Pre-solve resolved PML configuration across all active PML domains of a mesh.
#[derive(Clone, Debug, Default)]
pub(crate) struct ResolvedPmlDomains {
    pub domains: BTreeMap<i32, ResolvedPmlDomain>,
    pub diagnostics: Vec<String>,
}

fn domain_wafer_stiffness(materials: &femcore::props::Materials, dom: i32) -> Result<M6, String> {
    if let Some(st) = materials.aniso.get(&dom).and_then(|a| a.stiffness.as_ref()) {
        Ok(stiffness_global(st, materials.aniso[&dom].euler_deg))
    } else {
        let e = materials.get(dom, "E").ok_or_else(|| format!("domain {dom} has no value for 'E'"))?;
        let nu = materials.get(dom, "nu").ok_or_else(|| format!("domain {dom} has no value for 'nu'"))?;
        if nu >= 0.5 || nu <= -1.0 {
            return Err("Poisson's ratio must be between -1 and 0.5".into());
        }
        if e <= 0.0 {
            return Err("Young's modulus must be positive".into());
        }
        Ok(isotropic_stiffness(e, nu))
    }
}

fn normalize_axis(axis: [f64; 3], dim: usize, axisym: bool) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if n > 1e-30 {
        [axis[0] / n, axis[1] / n, axis[2] / n]
    } else if dim == 2 && !axisym {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn pml_outward_normals(spec: &Pml, dim: usize, axisym: bool) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    match spec.coord {
        PmlCoord::Cartesian => {
            match spec.x {
                PmlDir::Plus => out.push([1.0, 0.0, 0.0]),
                PmlDir::Minus => out.push([-1.0, 0.0, 0.0]),
                PmlDir::Off => {}
            }
            if dim == 3 {
                match spec.y {
                    PmlDir::Plus => out.push([0.0, 1.0, 0.0]),
                    PmlDir::Minus => out.push([0.0, -1.0, 0.0]),
                    PmlDir::Off => {}
                }
            }
            match spec.z {
                PmlDir::Plus => out.push([0.0, 0.0, 1.0]),
                PmlDir::Minus => out.push([0.0, 0.0, -1.0]),
                PmlDir::Off => {}
            }
        }
        PmlCoord::Cylindrical { axis, radial_on, axial, .. } => {
            let a_hat = normalize_axis(axis, dim, axisym);
            match axial {
                PmlDir::Plus => out.push(a_hat),
                PmlDir::Minus => out.push([-a_hat[0], -a_hat[1], -a_hat[2]]),
                PmlDir::Off => {}
            }
            if radial_on {
                if axisym {
                    out.push([1.0, 0.0, 0.0]);
                } else if dim == 2 {
                    for k in 0..8 {
                        let th = 2.0 * std::f64::consts::PI * (k as f64) / 8.0;
                        out.push([th.cos(), 0.0, th.sin()]);
                    }
                } else {
                    let ref_v = if a_hat[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
                    let dot = ref_v[0] * a_hat[0] + ref_v[1] * a_hat[1] + ref_v[2] * a_hat[2];
                    let mut u1 = [ref_v[0] - dot * a_hat[0], ref_v[1] - dot * a_hat[1], ref_v[2] - dot * a_hat[2]];
                    let l1 = (u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]).sqrt();
                    u1 = [u1[0] / l1, u1[1] / l1, u1[2] / l1];
                    let u2 = [
                        a_hat[1] * u1[2] - a_hat[2] * u1[1],
                        a_hat[2] * u1[0] - a_hat[0] * u1[2],
                        a_hat[0] * u1[1] - a_hat[1] * u1[0],
                    ];
                    for k in 0..8 {
                        let th = 2.0 * std::f64::consts::PI * (k as f64) / 8.0;
                        let (s, c0) = th.sin_cos();
                        out.push([c0 * u1[0] + s * u2[0], c0 * u1[1] + s * u2[1], c0 * u1[2] + s * u2[2]]);
                    }
                }
            }
        }
        PmlCoord::Spherical { .. } => {
            if axisym {
                for k in 0..=4 {
                    let th = std::f64::consts::PI * (k as f64) / 4.0;
                    out.push([th.sin(), 0.0, th.cos()]);
                }
            } else {
                out.extend_from_slice(&[
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                    [1.0, 1.0, 1.0],
                    [-1.0, 1.0, 1.0],
                    [1.0, -1.0, 1.0],
                    [1.0, 1.0, -1.0],
                ]);
            }
        }
    }
    out
}

fn infer_curvilinear_extents(spec: &mut Pml, pts_wafer: &[[f64; 3]], dim: usize, axisym: bool) {
    if pts_wafer.is_empty() {
        return;
    }
    match &mut spec.coord {
        PmlCoord::Cartesian => {}
        PmlCoord::Spherical { center, r_in, thickness } => {
            if *thickness <= 0.0 {
                let (mut r_lo, mut r_hi) = (f64::INFINITY, 0.0_f64);
                for p in pts_wafer {
                    let dx = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                    let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                    r_lo = r_lo.min(r);
                    r_hi = r_hi.max(r);
                }
                if r_hi > r_lo {
                    *r_in = r_lo;
                    *thickness = r_hi - r_lo;
                }
            }
        }
        PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness } => {
            let a_hat = normalize_axis(*axis, dim, axisym);
            let (mut r_lo, mut r_hi) = (f64::INFINITY, 0.0_f64);
            let (mut z_lo, mut z_hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for p in pts_wafer {
                let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                let zeta = d[0] * a_hat[0] + d[1] * a_hat[1] + d[2] * a_hat[2];
                let rv = [d[0] - zeta * a_hat[0], d[1] - zeta * a_hat[1], d[2] - zeta * a_hat[2]];
                let r = (rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]).sqrt();
                r_lo = r_lo.min(r);
                r_hi = r_hi.max(r);
                z_lo = z_lo.min(zeta);
                z_hi = z_hi.max(zeta);
            }
            if *radial_on && *r_thickness <= 0.0 && r_hi > r_lo {
                *r_in = r_lo;
                *r_thickness = r_hi - r_lo;
            }
            if *axial != PmlDir::Off && *z_thickness <= 0.0 && z_hi > z_lo {
                *z_in = if *axial == PmlDir::Plus { z_lo } else { z_hi };
                *z_thickness = z_hi - z_lo;
            }
        }
    }
}

fn effective_pml_thickness(spec: &Pml, bounds: ([f64; 3], [f64; 3]), dim: usize) -> Option<f64> {
    let (lo, hi) = bounds;
    match spec.coord {
        PmlCoord::Cartesian => {
            let mut d_min = f64::INFINITY;
            if spec.x != PmlDir::Off {
                d_min = d_min.min(hi[0] - lo[0]);
            }
            if dim == 3 && spec.y != PmlDir::Off {
                d_min = d_min.min(hi[1] - lo[1]);
            }
            if spec.z != PmlDir::Off {
                d_min = d_min.min(hi[2] - lo[2]);
            }
            (d_min.is_finite() && d_min > 0.0).then_some(d_min)
        }
        PmlCoord::Cylindrical { r_thickness, radial_on, axial, z_thickness, .. } => {
            let mut d_min = f64::INFINITY;
            if radial_on && r_thickness > 0.0 {
                d_min = d_min.min(r_thickness);
            }
            if axial != PmlDir::Off && z_thickness > 0.0 {
                d_min = d_min.min(z_thickness);
            }
            (d_min.is_finite() && d_min > 0.0).then_some(d_min)
        }
        PmlCoord::Spherical { thickness, .. } => (thickness > 0.0).then_some(thickness),
    }
}

impl ResolvedPmlDomains {
    fn build_from_domain_points(
        dim: usize,
        domain_pts: BTreeMap<i32, Vec<[f64; 3]>>,
        spec: &Problem,
        axisym: bool,
        f_search: f64,
    ) -> Result<Self, String> {
        let mut domains = BTreeMap::new();
        let mut diagnostics = Vec::new();

        for (&dom, raw_pml) in &spec.materials.pml {
            if spec.inactive.contains(&dom) {
                continue;
            }
            let Some(pts) = domain_pts.get(&dom) else { continue };
            let mut eff = raw_pml.clone();
            infer_curvilinear_extents(&mut eff, pts, dim, axisym);
            if !eff.is_active() {
                continue;
            }

            let bounds = if let Some(b) = eff.bounds {
                b
            } else {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for p in pts {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                        hi[k] = hi[k].max(p[k]);
                    }
                }
                (lo, hi)
            };

            let mut c6_w = domain_wafer_stiffness(&spec.materials, dom)?;
            let piezo_info: Option<([[f64; 6]; 3], M3)> = if spec.uses_piezoelectric() {
                spec.materials.aniso.get(&dom).and_then(|an| an.piezo_e.as_ref()).map(|e_tb| {
                    let e_sol = piezo_e_global(e_tb, spec.materials.aniso[&dom].euler_deg);
                    let eps_r = spec
                        .materials
                        .aniso
                        .get(&dom)
                        .and_then(|a| a.permittivity.as_ref())
                        .map(|t| tensor_global(t, spec.materials.aniso[&dom].euler_deg))
                        .unwrap_or_else(|| {
                            let s = spec.materials.get(dom, "eps_r").unwrap_or(1.0);
                            [[s, 0.0, 0.0], [0.0, s, 0.0], [0.0, 0.0, s]]
                        });
                    let mut eps_abs = [[0.0; 3]; 3];
                    for i in 0..3 {
                        for j in 0..3 {
                            eps_abs[i][j] = EPSILON_0 * eps_r[i][j];
                        }
                    }
                    (e_sol, eps_abs)
                })
            } else {
                None
            };
            if axisym {
                // In wafer `(x, y, z)`, the symmetry axis is `z` (index 2).
                c6_w = crate::axisym::average_about_axis(&c6_w, 2, 16);
            }
            let rho = spec.materials.get(dom, "rho").ok_or_else(|| format!("domain {dom} has no value for 'rho'"))?;

            if eff.auto_scale && f_search > 0.0 {
                let (c_min, c_max) = if let Some((ref e_sol, ref eps_abs)) = piezo_info {
                    piezo_wave_speed_range(&c6_w, e_sol, eps_abs, rho)
                } else {
                    wave_speed_range(&c6_w, rho)
                };
                let lam_max = (c_max / f_search).max(1e-30);
                let lam_min = (c_min / (1.5 * f_search)).clamp(1e-30, lam_max);
                let d_act = effective_pml_thickness(&eff, bounds, dim);
                let is_rat = matches!(eff.profile, PmlProfile::Rational { .. });
                let cal = auto_scale_params(lam_min, lam_max, d_act, eff.order, eff.target_reflection, true, is_rat);
                eff.strength = cal.psi_max;
                eff.profile = if is_rat {
                    PmlProfile::Rational { curvature: cal.rational_curvature() }
                } else {
                    PmlProfile::KappaPolynomial { kappa_max: cal.kappa_max }
                };
            }

            // Pre-solve Christoffel backward-wave check for anisotropic (or piezoelectric) domains
            if spec.materials.aniso.get(&dom).and_then(|a| a.stiffness.as_ref()).is_some() || piezo_info.is_some() {
                let mut worst =
                    BackwardWaveCheck { min_projection: 1.0, has_backward_wave: false, recommended_mpml_ratio: 0.0 };
                for n_out in pml_outward_normals(&eff, dim, axisym) {
                    let chk = if let Some((ref e_sol, ref eps_abs)) = piezo_info {
                        piezo_backward_wave_mpml_ratio(&c6_w, e_sol, eps_abs, rho, n_out)
                    } else {
                        backward_wave_mpml_ratio(&c6_w, rho, n_out)
                    };
                    if chk.min_projection < worst.min_projection {
                        worst.min_projection = chk.min_projection;
                    }
                    if chk.has_backward_wave {
                        worst.has_backward_wave = true;
                        worst.recommended_mpml_ratio = worst.recommended_mpml_ratio.max(chk.recommended_mpml_ratio);
                    }
                }
                if worst.has_backward_wave {
                    let p_min = worst.recommended_mpml_ratio / 1.25;
                    match eff.stabilization {
                        PmlStabilization::Auto => {
                            let ratio = worst.recommended_mpml_ratio;
                            eff.stabilization = PmlStabilization::Mpml { ratio };
                            diagnostics.push(format!(
                                "PML domain {dom}: Christoffel check detected backward-wave branch (min v_g.n/|v_g| = {:.4}, p_min = {:.4}); auto-enabled M-PML tangential damping (ratio = {:.4})",
                                worst.min_projection, p_min, ratio
                            ));
                        }
                        PmlStabilization::Off => {
                            diagnostics.push(format!(
                                "PML domain {dom}: Christoffel check detected backward-wave branch (min v_g.n/|v_g| = {:.4}, p_min = {:.4}) with stabilization = Off",
                                worst.min_projection, p_min
                            ));
                        }
                        PmlStabilization::Mpml { ratio } => {
                            diagnostics.push(format!(
                                "PML domain {dom}: Christoffel check detected backward-wave branch (p_min = {:.4}); using explicit M-PML (ratio = {ratio:.4})",
                                p_min
                            ));
                        }
                        PmlStabilization::Adiabatic { eta_max } => {
                            diagnostics.push(format!(
                                "PML domain {dom}: Christoffel check detected backward-wave branch (p_min = {:.4}); using adiabatic stiffness damping (eta_max = {eta_max:.4})",
                                p_min
                            ));
                        }
                    }
                }
            }

            domains.insert(dom, ResolvedPmlDomain { spec: eff, bounds });
        }

        Ok(Self { domains, diagnostics })
    }

    pub fn from_simplex(mesh: &Discretisation, spec: &Problem, axisym: bool, f_search: f64) -> Result<Self, String> {
        let d = mesh.dim;
        let mut dom_verts: BTreeMap<i32, BTreeSet<u32>> = BTreeMap::new();
        for ci in 0..mesh.n_cells() {
            let dom = mesh.cell_labels[ci];
            if spec.materials.pml.contains_key(&dom) && !spec.inactive.contains(&dom) {
                dom_verts.entry(dom).or_default().extend(mesh.cell(ci).iter().copied());
            }
        }
        let mut dom_pts = BTreeMap::new();
        for (dom, vs) in dom_verts {
            let pts: Vec<[f64; 3]> =
                vs.into_iter().map(|v| femcore::section::wafer_point(mesh.point(v as usize), d)).collect();
            dom_pts.insert(dom, pts);
        }
        Self::build_from_domain_points(d, dom_pts, spec, axisym, f_search)
    }

    pub fn from_mixed(
        mixed: &MixedDiscretisation,
        spec: &Problem,
        axisym: bool,
        f_search: f64,
    ) -> Result<Self, String> {
        let d = mixed.dim;
        let mut dom_verts: BTreeMap<i32, BTreeSet<u32>> = BTreeMap::new();
        for ci in 0..mixed.n_cells() {
            let dom = mixed.cell_labels[ci];
            if spec.materials.pml.contains_key(&dom) && !spec.inactive.contains(&dom) {
                dom_verts.entry(dom).or_default().extend(mixed.cell(ci).iter().copied());
            }
        }
        let mut dom_pts = BTreeMap::new();
        for (dom, vs) in dom_verts {
            let pts: Vec<[f64; 3]> =
                vs.into_iter().map(|v| femcore::section::wafer_point(mixed.point(v as usize), d)).collect();
            dom_pts.insert(dom, pts);
        }
        Self::build_from_domain_points(d, dom_pts, spec, axisym, f_search)
    }

    #[inline]
    pub fn contains(&self, dom: i32) -> bool {
        self.domains.contains_key(&dom)
    }

    #[inline]
    pub fn eval_jacobian(&self, dom: i32, x_wafer: [f64; 3], dim: usize, axisym: bool) -> (PmlJacobian, C) {
        if let Some(rd) = self.domains.get(&dom) {
            let jac = pml_jacobian(&rd.spec, Some(rd.bounds), x_wafer, dim, axisym);
            let loss = adiabatic_loss_factor(&rd.spec, jac.s_dist);
            (jac, loss)
        } else {
            (PmlJacobian::identity(), c(1.0, 0.0))
        }
    }

    /// True if an exterior face with corner coordinates `pts_wafer` belonging to PML domain `dom`
    /// is an outer PML termination wall (`clamp_outer == true` and all face corners sit at `t >= 1 - tol`
    /// along at least one active stretching coordinate).
    pub fn is_outer_clamped_facet(
        &self,
        dom: i32,
        pts_wafer: &[[f64; 3]],
        dim: usize,
        axisym: bool,
        axis_tol: f64,
    ) -> bool {
        let Some(rd) = self.domains.get(&dom) else {
            return false;
        };
        if !rd.spec.clamp_outer || pts_wafer.is_empty() {
            return false;
        }
        if axisym && pts_wafer.iter().all(|p| p[0] <= axis_tol) {
            return false;
        }
        is_outer_termination_face(&rd.spec, rd.bounds, pts_wafer, dim, axisym)
    }
}

/// Checks whether all corner points `pts_wafer` of an exterior face sit at the outer end (`t >= 1 - tol`)
/// of at least one active coordinate of `spec` (shared by both PMLs and Static Infinite Elements).
pub(crate) fn is_outer_termination_face(
    spec: &Pml,
    bounds: ([f64; 3], [f64; 3]),
    pts_wafer: &[[f64; 3]],
    dim: usize,
    axisym: bool,
) -> bool {
    if pts_wafer.is_empty() {
        return false;
    }
    let (lo, hi) = bounds;
    match spec.coord {
        PmlCoord::Cartesian => {
            let dirs = if dim == 2 { [spec.x, PmlDir::Off, spec.z] } else { [spec.x, spec.y, spec.z] };
            for a in 0..3 {
                let len = hi[a] - lo[a];
                if len <= 1e-300 || dirs[a] == PmlDir::Off {
                    continue;
                }
                let all_at_outer = pts_wafer.iter().all(|p| {
                    let t = match dirs[a] {
                        PmlDir::Plus => (p[a] - lo[a]) / len,
                        PmlDir::Minus => (hi[a] - p[a]) / len,
                        PmlDir::Off => 0.0,
                    };
                    t >= 1.0 - OUTER_T_TOL
                });
                if all_at_outer {
                    return true;
                }
            }
            false
        }
        PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness } => {
            let a_hat = normalize_axis(axis, dim, axisym);
            if radial_on && r_thickness > 0.0 {
                let all_radial_outer = pts_wafer.iter().all(|p| {
                    let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                    let zeta = d[0] * a_hat[0] + d[1] * a_hat[1] + d[2] * a_hat[2];
                    let rv = [d[0] - zeta * a_hat[0], d[1] - zeta * a_hat[1], d[2] - zeta * a_hat[2]];
                    let r = (rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]).sqrt();
                    (r - r_in) / r_thickness >= 1.0 - OUTER_T_TOL
                });
                if all_radial_outer {
                    return true;
                }
            }
            if axial != PmlDir::Off && z_thickness > 0.0 {
                let all_axial_outer = pts_wafer.iter().all(|p| {
                    let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                    let zeta = d[0] * a_hat[0] + d[1] * a_hat[1] + d[2] * a_hat[2];
                    let tz = match axial {
                        PmlDir::Plus => (zeta - z_in) / z_thickness,
                        PmlDir::Minus => (z_in - zeta) / z_thickness,
                        PmlDir::Off => 0.0,
                    };
                    tz >= 1.0 - OUTER_T_TOL
                });
                if all_axial_outer {
                    return true;
                }
            }
            false
        }
        PmlCoord::Spherical { center, r_in, thickness } => {
            if thickness <= 0.0 {
                return false;
            }
            pts_wafer.iter().all(|p| {
                let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                let big_r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                (big_r - r_in) / thickness >= 1.0 - OUTER_T_TOL
            })
        }
    }
}

fn explicit_protected_facet_labels(spec: &Problem) -> BTreeSet<i32> {
    let mut out = BTreeSet::new();
    for bc in &spec.bcs {
        match bc.kind {
            BcKind::Fixed { .. }
            | BcKind::LowReflecting
            | BcKind::Ground
            | BcKind::Potential(_)
            | BcKind::SurfaceCharge(_)
            | BcKind::FloatingPotential { .. }
            | BcKind::Terminal { .. } => {
                out.extend(bc.labels.iter().copied());
            }
            BcKind::Load(_) | BcKind::Conductor { .. } => {}
        }
    }
    for p in &spec.periodic {
        out.extend(p.source.iter().copied());
        out.extend(p.dest.iter().copied());
    }
    out
}

/// Returns sorted unique scalar DOFs (`0 .. mesh.n_dofs(order)`) lying on outer clamped PML termination
/// faces of a [`SimplexMesh`].
pub(crate) fn outer_clamped_dofs_simplex(
    mesh: &Discretisation,
    order: usize,
    spec: &Problem,
    resolved: &ResolvedPmlDomains,
    axisym: bool,
    axis_tol: f64,
) -> Vec<usize> {
    if resolved.domains.values().all(|d| !d.spec.clamp_outer) {
        return Vec::new();
    }
    let d = mesh.dim;
    let nv = mesh.n_points();
    let protected = explicit_protected_facet_labels(spec);
    let mut clamped = BTreeSet::new();

    if d == 2 {
        let mut facet_label_of: HashMap<[u32; 2], i32> = HashMap::new();
        for f in 0..mesh.n_facets() {
            let fv = mesh.facet(f);
            let k = if fv[0] <= fv[1] { [fv[0], fv[1]] } else { [fv[1], fv[0]] };
            facet_label_of.insert(k, mesh.facet_labels[f]);
        }
        let mut adj: HashMap<[u32; 2], (u32, usize)> = HashMap::new();
        for ci in 0..mesh.n_cells() {
            if spec.inactive.contains(&mesh.cell_labels[ci]) {
                continue;
            }
            let cv = mesh.cell(ci);
            for &(i, j) in &[(0usize, 1usize), (1, 2), (2, 0)] {
                let (a, b) = (cv[i], cv[j]);
                let k = if a <= b { [a, b] } else { [b, a] };
                let e = adj.entry(k).or_insert((0, ci));
                e.0 += 1;
            }
        }
        for (edge, (cnt, ci)) in adj {
            if cnt != 1 {
                continue;
            }
            if facet_label_of.get(&edge).is_some_and(|lab| protected.contains(lab)) {
                continue;
            }
            let dom = mesh.cell_labels[ci];
            let pts_w = [
                femcore::section::wafer_point(mesh.point(edge[0] as usize), 2),
                femcore::section::wafer_point(mesh.point(edge[1] as usize), 2),
            ];
            if resolved.is_outer_clamped_facet(dom, &pts_w, 2, axisym, axis_tol) {
                clamped.insert(edge[0] as usize);
                clamped.insert(edge[1] as usize);
                if order == 2 {
                    clamped.insert(nv + mesh.edge_index(edge[0], edge[1]));
                }
            }
        }
    } else {
        let mut facet_label_of: HashMap<[u32; 3], i32> = HashMap::new();
        for f in 0..mesh.n_facets() {
            let fv = mesh.facet(f);
            facet_label_of.insert(sorted_tri_face([fv[0], fv[1], fv[2]]), mesh.facet_labels[f]);
        }
        let mut adj: HashMap<[u32; 3], (u32, usize)> = HashMap::new();
        for ci in 0..mesh.n_cells() {
            if spec.inactive.contains(&mesh.cell_labels[ci]) {
                continue;
            }
            let cv = mesh.cell(ci);
            for &[i, j, k] in &[[1usize, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]] {
                let key = sorted_tri_face([cv[i], cv[j], cv[k]]);
                let e = adj.entry(key).or_insert((0, ci));
                e.0 += 1;
            }
        }
        for (tri, (cnt, ci)) in adj {
            if cnt != 1 {
                continue;
            }
            if facet_label_of.get(&tri).is_some_and(|lab| protected.contains(lab)) {
                continue;
            }
            let dom = mesh.cell_labels[ci];
            let pts_w = [
                femcore::section::wafer_point(mesh.point(tri[0] as usize), 3),
                femcore::section::wafer_point(mesh.point(tri[1] as usize), 3),
                femcore::section::wafer_point(mesh.point(tri[2] as usize), 3),
            ];
            if resolved.is_outer_clamped_facet(dom, &pts_w, 3, false, 0.0) {
                for &v in &tri {
                    clamped.insert(v as usize);
                }
                if order == 2 {
                    clamped.insert(nv + mesh.edge_index(tri[0], tri[1]));
                    clamped.insert(nv + mesh.edge_index(tri[0], tri[2]));
                    clamped.insert(nv + mesh.edge_index(tri[1], tri[2]));
                }
            }
        }
    }

    clamped.into_iter().collect()
}

/// Returns sorted unique scalar DOFs (`0 .. dof_map.n_dofs`) lying on outer clamped termination faces
/// of a [`MixedDiscretisation`] according to `predicate(dom, &pts_wafer)`.
pub(crate) fn outer_faces_mixed_where<F>(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    inactive: &[i32],
    protected_labels: &BTreeSet<i32>,
    mut predicate: F,
) -> Vec<usize>
where
    F: FnMut(i32, &[[f64; 3]]) -> bool,
{
    let d = mixed.dim;
    let nv = dof_map.n_vertices;
    let mut clamped = BTreeSet::new();

    if d == 2 {
        let mut facet_label_of: HashMap<(u32, u32), i32> = HashMap::new();
        for f in 0..mixed.n_facets() {
            if mixed.facet_type(f) == FacetType::Edge2 {
                let fv = mixed.facet(f);
                facet_label_of.insert((fv[0].min(fv[1]), fv[0].max(fv[1])), mixed.facet_labels[f]);
            }
        }
        let mut adj: HashMap<(u32, u32), (u32, usize)> = HashMap::new();
        for ci in 0..mixed.n_cells() {
            if inactive.contains(&mixed.cell_labels[ci]) {
                continue;
            }
            let vs = mixed.cell(ci);
            for &(i, j) in cell_local_edges(mixed.cell_type(ci)) {
                let (a, b) = (vs[i], vs[j]);
                let e = adj.entry((a.min(b), a.max(b))).or_insert((0, ci));
                e.0 += 1;
            }
        }
        for (edge, (cnt, ci)) in adj {
            if cnt != 1 {
                continue;
            }
            if facet_label_of.get(&edge).is_some_and(|lab| protected_labels.contains(lab)) {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let pts_w = [
                femcore::section::wafer_point(mixed.point(edge.0 as usize), 2),
                femcore::section::wafer_point(mixed.point(edge.1 as usize), 2),
            ];
            if predicate(dom, &pts_w) {
                clamped.insert(edge.0 as usize);
                clamped.insert(edge.1 as usize);
                if dof_map.space != LagrangeSpace::Order1 {
                    if let Some(eid) = mixed.edge_index(edge.0, edge.1) {
                        clamped.insert(nv + eid);
                    }
                }
            }
        }
    } else {
        let mut tri_label_of: HashMap<[u32; 3], i32> = HashMap::new();
        let mut quad_label_of: HashMap<[u32; 4], i32> = HashMap::new();
        for f in 0..mixed.n_facets() {
            let fv = mixed.facet(f);
            match mixed.facet_type(f) {
                FacetType::Tri3 => {
                    tri_label_of.insert(sorted_tri_face([fv[0], fv[1], fv[2]]), mixed.facet_labels[f]);
                }
                FacetType::Quad4 => {
                    quad_label_of.insert(canonical_quad_face([fv[0], fv[1], fv[2], fv[3]]), mixed.facet_labels[f]);
                }
                FacetType::Edge2 => {}
            }
        }
        let mut tri_adj: HashMap<[u32; 3], (u32, usize)> = HashMap::new();
        let mut quad_adj: HashMap<[u32; 4], (u32, usize)> = HashMap::new();
        for ci in 0..mixed.n_cells() {
            if inactive.contains(&mixed.cell_labels[ci]) {
                continue;
            }
            let vs = mixed.cell(ci);
            let ct = mixed.cell_type(ci);
            for lf in cell_local_tri_faces(ct) {
                let e = tri_adj.entry(sorted_tri_face([vs[lf[0]], vs[lf[1]], vs[lf[2]]])).or_insert((0, ci));
                e.0 += 1;
            }
            for lq in cell_local_quad_faces(ct) {
                let e = quad_adj
                    .entry(canonical_quad_face([vs[lq[0]], vs[lq[1]], vs[lq[2]], vs[lq[3]]]))
                    .or_insert((0, ci));
                e.0 += 1;
            }
        }
        for (tri, (cnt, ci)) in tri_adj {
            if cnt != 1 {
                continue;
            }
            if tri_label_of.get(&tri).is_some_and(|lab| protected_labels.contains(lab)) {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let pts_w = tri.map(|v| femcore::section::wafer_point(mixed.point(v as usize), 3));
            if predicate(dom, &pts_w) {
                for &v in &tri {
                    clamped.insert(v as usize);
                }
                if dof_map.space != LagrangeSpace::Order1 {
                    for &(i, j) in &TRI_EDGES {
                        if let Some(eid) = mixed.edge_index(tri[i], tri[j]) {
                            clamped.insert(nv + eid);
                        }
                    }
                }
            }
        }
        for (quad, (cnt, ci)) in quad_adj {
            if cnt != 1 {
                continue;
            }
            if quad_label_of.get(&quad).is_some_and(|lab| protected_labels.contains(lab)) {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let pts_w = quad.map(|v| femcore::section::wafer_point(mixed.point(v as usize), 3));
            if predicate(dom, &pts_w) {
                for &v in &quad {
                    clamped.insert(v as usize);
                }
                if dof_map.space != LagrangeSpace::Order1 {
                    for &(i, j) in &QUAD_EDGES {
                        if let Some(eid) = mixed.edge_index(quad[i], quad[j]) {
                            clamped.insert(nv + eid);
                        }
                    }
                    if dof_map.space == LagrangeSpace::Order2Tensor {
                        if let Some(qfid) = mixed.quad_face_index(quad) {
                            clamped.insert(nv + dof_map.n_edges + qfid);
                        }
                    }
                }
            }
        }
    }

    clamped.into_iter().collect()
}

pub(crate) fn outer_clamped_dofs_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    resolved: &ResolvedPmlDomains,
    axisym: bool,
    axis_tol: f64,
) -> Vec<usize> {
    if resolved.domains.values().all(|d| !d.spec.clamp_outer) {
        return Vec::new();
    }
    let protected = explicit_protected_facet_labels(spec);
    let d = mixed.dim;
    outer_faces_mixed_where(mixed, dof_map, &spec.inactive, &protected, |dom, pts_w| {
        resolved.is_outer_clamped_facet(dom, pts_w, d, axisym, axis_tol)
    })
}

/// Complex strain-displacement matrix from complex shape gradients `g[a*d + k]`.
pub(crate) fn b_matrix_c(g: &[C], nloc: usize, d: usize, b: &mut [C]) {
    let ns = if d == 2 { 3 } else { 6 };
    b[..ns * d * nloc].iter_mut().for_each(|x| *x = c(0.0, 0.0));
    let nd = d * nloc;
    for k in 0..d {
        for a in 0..nloc {
            b[k * nd + k * nloc + a] = g[a * d + k];
        }
    }
    let shear: &[(usize, usize)] = if d == 2 { &[(0, 1)] } else { &[(0, 1), (1, 2), (0, 2)] };
    for (r, &(p, q)) in shear.iter().enumerate() {
        let row = (d + r) * nd;
        for a in 0..nloc {
            b[row + p * nloc + a] = g[a * d + q];
            b[row + q * nloc + a] = g[a * d + p];
        }
    }
}

pub(crate) fn simplex_barycentric(re: &RefElement, q: usize, d: usize, order: usize) -> [f64; 4] {
    let nv = &re.n[q];
    if order == 1 {
        let mut out = [0.0; 4];
        out[..=d].copy_from_slice(&nv[..=d]);
        return out;
    }
    if d == 2 {
        [nv[0] + 0.5 * (nv[3] + nv[4]), nv[1] + 0.5 * (nv[3] + nv[5]), nv[2] + 0.5 * (nv[4] + nv[5]), 0.0]
    } else {
        // 3D P2 Tet edges: (0,1)=4, (0,2)=5, (0,3)=6, (1,2)=7, (1,3)=8, (2,3)=9
        [
            nv[0] + 0.5 * (nv[4] + nv[5] + nv[6]),
            nv[1] + 0.5 * (nv[4] + nv[7] + nv[8]),
            nv[2] + 0.5 * (nv[5] + nv[7] + nv[9]),
            nv[3] + 0.5 * (nv[6] + nv[8] + nv[9]),
        ]
    }
}

pub(crate) fn real_quad_norm(m1: &Csr, u: &[C], n1: usize, d: usize) -> f64 {
    let mut acc = 0.0;
    for comp in 0..d {
        let off = comp * n1;
        for r in 0..n1 {
            let mut row = c(0.0, 0.0);
            for k in m1.indptr()[r]..m1.indptr()[r + 1] {
                row += u[off + m1.indices()[k] as usize] * m1.values()[k];
            }
            acc += (u[off + r].conj() * row).re;
        }
    }
    acc
}

pub fn solve_solid_pml(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_solid_pml(mesh, spec)?;
    let solved = solve_pml_modes(&asm)?;
    let r = pml_fields(mesh, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

/// Native mixed-cell (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) PML eigenfrequency solve.
pub fn solve_solid_pml_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let order = spec.order;
    let d = mixed.dim;
    let Study::Eigenfrequency { n_modes, f_search } = spec.study else {
        return Err("PML needs the Eigenfrequency study".into());
    };
    crate::solid::mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    let thk = if d == 2 { spec.thickness } else { 1.0 };

    let c6: Vec<M6> = crate::solid::mixed_cell_stiffness_skipping(&spec.materials, mixed, &spec.inactive)?
        .iter()
        .map(|c| femcore::section::plane_stiffness(c, d))
        .collect();
    let rho = crate::solid::mixed_cell_rho(&spec.materials, mixed, &spec.inactive)?;
    let resolved = ResolvedPmlDomains::from_mixed(mixed, spec, false, f_search)?;

    let mut fixed_set: BTreeSet<usize> =
        crate::solid::mixed_untouched_dofs(mixed, &dof_map, d, &spec.inactive).into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = dof_map.boundary_dofs(mixed, &bc.labels);
            for (k, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                if fixed[w] {
                    for &dof in &bd {
                        fixed_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    for dof in outer_clamped_dofs_mixed(mixed, &dof_map, spec, &resolved, false, 0.0) {
        for k in 0..d {
            fixed_set.insert(dof + k * n1);
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();

    let nc = mixed.n_cells();
    let ns = if d == 2 { 3 } else { 6 };
    let nt = if nc < 256 { 1 } else { femcore::par::threads_up_to(8) };
    let chunk = nc.div_ceil(nt).max(1);

    let chunks: Vec<(CCoo, CCoo)> = (0..nt)
        .into_par_iter()
        .map(|tid| {
            let mut kc = CCoo::new(n, n);
            let mut mc = CCoo::new(n, n);
            let mut vals = Vec::new();
            let mut dxi = Vec::new();
            let mut sn = Vec::new();
            let mut sdn = Vec::new();
            let mut gc = Vec::new();
            let mut b = Vec::new();
            let mut db = Vec::new();
            let mut ke = Vec::new();
            let mut me = Vec::new();
            let mut full = Vec::new();

            for ci in (tid * chunk)..((tid + 1) * chunk).min(nc) {
                let dom = mixed.cell_labels[ci];
                if spec.inactive.contains(&dom) {
                    continue;
                }
                let ct = mixed.cell_type(ci);
                let iso = IsoparametricCell::from_mesh(mixed, ci);
                let in_pml = resolved.contains(dom);
                let deg = if in_pml || iso.curved { 2 * order + 4 } else { 2 * order + 2 };
                let rule = quad_rule_mixed(ct, deg);
                let cdofs = dof_map.cell_dofs(ci);
                let nloc = cdofs.len();
                let nd = d * nloc;
                let dm = d_matrix(&c6[ci], d, spec.plane);

                gc.resize(nloc * d, c(0.0, 0.0));
                b.resize(ns * nd, c(0.0, 0.0));
                db.resize(ns * nd, c(0.0, 0.0));
                ke.clear();
                ke.resize(nd * nd, c(0.0, 0.0));
                me.clear();
                me.resize(nloc * nloc, c(0.0, 0.0));

                for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                    let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                    lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                    let x_w = femcore::section::wafer_point(&jg.x[..d], d);
                    let (jp, loss) = resolved.eval_jacobian(dom, x_w, d, false);

                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            let mut acc = c(0.0, 0.0);
                            for k in 0..d {
                                acc += jp.j_inv[k][i] * gp[k];
                            }
                            gc[a * d + i] = acc;
                        }
                    }
                    b_matrix_c(&gc, nloc, d, &mut b);

                    let w_phys = wq * jg.det.abs() * thk;
                    let w_k = jp.det * loss * w_phys;
                    let w_m = jp.det * (w_phys * rho[ci]);

                    for r in 0..ns {
                        for j in 0..nd {
                            let mut acc = c(0.0, 0.0);
                            for t in 0..ns {
                                acc += b[t * nd + j] * dm[r * ns + t];
                            }
                            db[r * nd + j] = acc;
                        }
                    }
                    for i in 0..nd {
                        for j in 0..nd {
                            let mut acc = c(0.0, 0.0);
                            for r in 0..ns {
                                acc += b[r * nd + i] * db[r * nd + j];
                            }
                            ke[i * nd + j] += acc * w_k;
                        }
                    }
                    for a in 0..nloc {
                        let wa = w_m * vals[a];
                        for bb in 0..nloc {
                            me[a * nloc + bb] += wa * vals[bb];
                        }
                    }
                }

                full.clear();
                for k in 0..d {
                    full.extend(cdofs.iter().map(|&gd| gd + k * n1));
                }
                kc.push_block(&full, &full, &ke);
                for k in 0..d {
                    let comp_dofs: Vec<usize> = cdofs.iter().map(|&gd| gd + k * n1).collect();
                    mc.push_block(&comp_dofs, &comp_dofs, &me);
                }
            }
            (kc, mc)
        })
        .collect();

    let mut kc_all = CCoo::new(n, n);
    let mut mc_all = CCoo::new(n, n);
    for (kc, mc) in chunks {
        kc_all.extend(&kc);
        mc_all.extend(&mc);
    }
    let (kfull, mfull) = (kc_all.to_csr(), mc_all.to_csr());

    let (kff, mff, free, red, n_pairs) = if spec.periodic.is_empty() {
        let free = free_dofs(n, &fixed);
        if n_modes < 1 || n_modes + 1 >= free.len() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        (kfull.submatrix(&free, &free), mfull.submatrix(&free, &free), free, None, 0)
    } else {
        let mut relations = Vec::new();
        for p in &spec.periodic {
            relations.extend(mixed_nodal_relations(mixed, &dof_map, d, p)?);
        }
        let red = Reduction::build(n, &relations, &fixed);
        let pm = red.prolong();
        let nred = red.n_reduced();
        if n_modes < 1 || n_modes + 1 >= nred {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        (congruence(&pm, &kfull), congruence(&pm, &mfull), Vec::new(), Some(red), spec.periodic.len())
    };

    let pml_cells: Vec<bool> = (0..nc).map(|ci| resolved.contains(mixed.cell_labels[ci])).collect();
    let any_pml = pml_cells.iter().any(|&p| p);
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
    let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
    let m1_all = assemble_mass_mixed(mixed, &dof_map, &rho_all_t);
    let m1_pml = assemble_mass_mixed(mixed, &dof_map, &rho_pml_t);

    let asm = PmlAssembly {
        order,
        d,
        n,
        kff,
        mff,
        free,
        red,
        n_pairs,
        n_modes,
        f_search,
        has_fixed: !fixed.is_empty(),
        diagnostics: resolved.diagnostics,
        m1_phys_all: Some(m1_all),
        m1_phys_pml: Some(m1_pml),
        any_pml,
        has_interior,
    };
    let solved = solve_pml_modes(&asm)?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let r = pml_fields_mixed(mixed, &dof_map, &proj, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

/// Stage 1 of the solid kernel with perfectly matched layers: the complex stiffness and the mass (weighted by the
/// curvilinear stretch determinant) on the free dofs.
pub struct PmlAssembly {
    pub order: usize,
    pub d: usize,
    pub n: usize,
    pub kff: CCsr,
    pub mff: CCsr,
    pub free: Vec<usize>,
    pub red: Option<Reduction>,
    pub n_pairs: usize,
    pub n_modes: usize,
    pub f_search: f64,
    pub has_fixed: bool,
    pub(crate) diagnostics: Vec<String>,
    pub(crate) m1_phys_all: Option<Csr>,
    pub(crate) m1_phys_pml: Option<Csr>,
    pub(crate) any_pml: bool,
    pub(crate) has_interior: bool,
}

/// Stage 2 output: eigenvalue index and complex frequency (Hz) of every mode sorted by real part, the free-dof vectors.
pub struct PmlSolve {
    pub freqs: Vec<(usize, C)>,
    pub vectors: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
}

pub fn assemble_solid_pml(mesh: &Discretisation, spec: &Problem) -> Result<PmlAssembly, String> {
    check_simplex_straight(mesh)?;
    let order = spec.order;
    let mut st = setup(mesh, spec)?;
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let Study::Eigenfrequency { n_modes, f_search } = spec.study else {
        return Err("PML needs the Eigenfrequency study".into());
    };
    let thk = if d == 2 { spec.thickness } else { 1.0 };
    let rho = crate::solid::cell_rho(&spec.materials, mesh, &spec.inactive)?;
    let resolved = ResolvedPmlDomains::from_simplex(mesh, spec, false, f_search)?;

    let outer_clamped = outer_clamped_dofs_simplex(mesh, order, spec, &resolved, false, 0.0);
    if !outer_clamped.is_empty() {
        let mut fixed_set: BTreeSet<usize> = st.fixed.iter().copied().collect();
        for dof in outer_clamped {
            for k in 0..d {
                fixed_set.insert(dof + k * n1);
            }
        }
        st.fixed = fixed_set.into_iter().collect();
    }

    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let rm = RefElement::new(d, order, 2 * order);
    let re_pml = RefElement::new(d, order, 2 * order + 2);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let mut mloc = vec![0.0; nloc * nloc];
    for q in 0..rm.w.len() {
        for a in 0..nloc {
            for bb in 0..nloc {
                mloc[a * nloc + bb] += rm.w[q] * rm.n[q][a] * rm.n[q][bb];
            }
        }
    }

    let nc = mesh.n_cells();
    let nt = if nc < 512 { 1 } else { femcore::par::threads_up_to(8) };
    let chunk = nc.div_ceil(nt).max(1);

    let chunks: Vec<(CCoo, CCoo)> = (0..nt)
        .into_par_iter()
        .map(|tid| {
            let mut kc = CCoo::new(n, n);
            let mut mc = CCoo::new(n, n);
            let (mut scratch, mut full) = (Vec::new(), Vec::new());
            let mut grad = vec![0.0; nloc * d];
            let mut gc = vec![c(0.0, 0.0); nloc * d];
            let mut b = vec![c(0.0, 0.0); ns * nd];
            let mut db = vec![c(0.0, 0.0); ns * nd];
            let mut ke = vec![c(0.0, 0.0); nd * nd];
            let mut me = vec![c(0.0, 0.0); nloc * nloc];

            for ci in (tid * chunk)..((tid + 1) * chunk).min(nc) {
                let dom = mesh.cell_labels[ci];
                if spec.inactive.contains(&dom) {
                    continue;
                }
                let dm = d_matrix(&st.c6[ci], d, spec.plane);
                ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
                let vol_thk = geo.vol[ci] * thk;

                if resolved.contains(dom) {
                    me.iter_mut().for_each(|x| *x = c(0.0, 0.0));
                    let cv = mesh.cell(ci);
                    for q in 0..re_pml.w.len() {
                        shape_gradients(&re_pml, &geo, ci, q, &mut grad);
                        let lam = simplex_barycentric(&re_pml, q, d, order);
                        let mut x_phys = [0.0; 3];
                        for corner in 0..=d {
                            let pt = mesh.point(cv[corner] as usize);
                            for k in 0..d {
                                x_phys[k] += lam[corner] * pt[k];
                            }
                        }
                        let x_w = femcore::section::wafer_point(&x_phys[..d], d);
                        let (jac, loss) = resolved.eval_jacobian(dom, x_w, d, false);

                        for a in 0..nloc {
                            for i in 0..d {
                                let mut acc = c(0.0, 0.0);
                                for k in 0..d {
                                    acc += jac.j_inv[k][i] * grad[a * d + k];
                                }
                                gc[a * d + i] = acc;
                            }
                        }
                        b_matrix_c(&gc, nloc, d, &mut b);
                        let w_k = jac.det * loss * (re_pml.w[q] * vol_thk);
                        let w_m = jac.det * (re_pml.w[q] * vol_thk * rho[ci]);

                        for r in 0..ns {
                            for j in 0..nd {
                                let mut acc = c(0.0, 0.0);
                                for t in 0..ns {
                                    acc += b[t * nd + j] * dm[r * ns + t];
                                }
                                db[r * nd + j] = acc;
                            }
                        }
                        for i in 0..nd {
                            for j in 0..nd {
                                let mut acc = c(0.0, 0.0);
                                for r in 0..ns {
                                    acc += b[r * nd + i] * db[r * nd + j];
                                }
                                ke[i * nd + j] += acc * w_k;
                            }
                        }
                        let nv = &re_pml.n[q];
                        for a in 0..nloc {
                            let wa = w_m * nv[a];
                            for bb in 0..nloc {
                                me[a * nloc + bb] += wa * nv[bb];
                            }
                        }
                    }
                    dof_layout(mesh, order, ci, &mut scratch, &mut full);
                    kc.push_block(&full, &full, &ke);
                    mesh.cell_dofs(ci, order, &mut scratch);
                    for comp in 0..d {
                        let dofs: Vec<usize> = scratch.iter().map(|&x| x + comp * n1).collect();
                        mc.push_block(&dofs, &dofs, &me);
                    }
                } else {
                    for q in 0..re.w.len() {
                        shape_gradients(&re, &geo, ci, q, &mut grad);
                        for a in 0..nloc {
                            for k in 0..d {
                                gc[a * d + k] = c(grad[a * d + k], 0.0);
                            }
                        }
                        b_matrix_c(&gc, nloc, d, &mut b);
                        for r in 0..ns {
                            for j in 0..nd {
                                let mut acc = c(0.0, 0.0);
                                for t in 0..ns {
                                    acc += b[t * nd + j] * dm[r * ns + t];
                                }
                                db[r * nd + j] = acc;
                            }
                        }
                        for i in 0..nd {
                            for j in 0..nd {
                                let mut acc = c(0.0, 0.0);
                                for r in 0..ns {
                                    acc += b[r * nd + i] * db[r * nd + j];
                                }
                                ke[i * nd + j] += acc * re.w[q];
                            }
                        }
                    }
                    ke.iter_mut().for_each(|x| *x *= vol_thk);
                    dof_layout(mesh, order, ci, &mut scratch, &mut full);
                    kc.push_block(&full, &full, &ke);
                    mesh.cell_dofs(ci, order, &mut scratch);
                    let mrow: Vec<C> = mloc.iter().map(|m| c(rho[ci] * vol_thk * m, 0.0)).collect();
                    for comp in 0..d {
                        let dofs: Vec<usize> = scratch.iter().map(|&x| x + comp * n1).collect();
                        mc.push_block(&dofs, &dofs, &mrow);
                    }
                }
            }
            (kc, mc)
        })
        .collect();

    let mut kc_all = CCoo::new(n, n);
    let mut mc_all = CCoo::new(n, n);
    for (kc, mc) in chunks {
        kc_all.extend(&kc);
        mc_all.extend(&mc);
    }
    let (kfull, mfull) = (kc_all.to_csr(), mc_all.to_csr());

    let (kff, mff, free, red, n_pairs) = if spec.periodic.is_empty() {
        let free = free_dofs(n, &st.fixed);
        if n_modes < 1 || n_modes + 1 >= free.len() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let kff = kfull.submatrix(&free, &free);
        let mff = mfull.submatrix(&free, &free);
        (kff, mff, free, None, 0)
    } else {
        let mut relations = vec![];
        for p in &spec.periodic {
            let mt = match_vertices(mesh, p)?;
            relations.extend(nodal_relations(mesh, order, d, &mt, p.phase)?);
        }
        let red = Reduction::build(n, &relations, &st.fixed);
        let pm = red.prolong();
        let nred = red.n_reduced();
        if n_modes < 1 || n_modes + 1 >= nred {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let kff = congruence(&pm, &kfull);
        let mff = congruence(&pm, &mfull);
        (kff, mff, vec![], Some(red), spec.periodic.len())
    };

    let pml_cells: Vec<bool> = (0..nc).map(|ci| resolved.contains(mesh.cell_labels[ci])).collect();
    let any_pml = pml_cells.iter().any(|&p| p);
    let has_interior = pml_cells.iter().zip(&rho).any(|(&p, &r)| !p && r > 0.0);
    let rho_all_t: Vec<f64> = rho.iter().map(|&r| r * thk).collect();
    let rho_pml_t: Vec<f64> = rho.iter().zip(&pml_cells).map(|(&r, &p)| if p { r * thk } else { 0.0 }).collect();
    let m1_all = assemble_mass(mesh, order, &rho_all_t);
    let m1_pml = assemble_mass(mesh, order, &rho_pml_t);

    Ok(PmlAssembly {
        order,
        d,
        n,
        kff,
        mff,
        free,
        red,
        n_pairs,
        n_modes,
        f_search,
        has_fixed: !st.fixed.is_empty(),
        diagnostics: resolved.diagnostics,
        m1_phys_all: Some(m1_all),
        m1_phys_pml: Some(m1_pml),
        any_pml,
        has_interior,
    })
}

/// Stage 2: the complex eigenpairs nearest the search frequency.
pub fn solve_pml_modes(asm: &PmlAssembly) -> Result<PmlSolve, String> {
    let sigma = if asm.f_search > 0.0 { (2.0 * std::f64::consts::PI * asm.f_search).powi(2) } else { -1.0 };
    let pairs = eigs_complex(&asm.kff, &asm.mff, c(sigma, 0.0), asm.n_modes, 1e-9, 600)?;
    let mut freqs: Vec<(usize, C)> = pairs
        .values
        .iter()
        .enumerate()
        .filter(|(i, _)| asm.red.is_none() || pairs.residuals.get(*i).is_none_or(|r| *r < 1e-6))
        .map(|(i, l)| {
            let mut f = l.sqrt() / (2.0 * std::f64::consts::PI);
            if f.re < 0.0 {
                f = -f;
            }
            (i, f)
        })
        .collect();
    freqs.sort_by(|a, b| a.1.re.partial_cmp(&b.1.re).unwrap());
    if asm.red.is_some() && freqs.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    Ok(PmlSolve { freqs, vectors: pairs.vectors.clone(), dim: pairs.dim, converged: pairs.converged })
}

pub(crate) fn push_pml_mode_scalars_and_label(
    j: usize,
    f: C,
    layer_fraction: f64,
    phys_mass_norm: f64,
    any_pml: bool,
    has_interior: bool,
    scalars: &mut Vec<(String, Scalar)>,
    unstable_count: &mut usize,
) -> (String, f64) {
    let is_amplifying = f.im < -AMPLIFYING_IM_TOL * f.re.abs().max(1e-30);
    if is_amplifying {
        *unstable_count += 1;
    }
    let q = if f.im.abs() > 0.0 { f.re / (2.0 * f.im) } else { f64::INFINITY };
    let label = if is_amplifying {
        format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]", j + 1, f.re, f.im, q)
    } else {
        format!("mode {}: {:.6e} {:+.3e}i Hz", j + 1, f.re, f.im)
    };
    scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
    if q.is_finite() {
        scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
    }
    scalars.push((format!("layer fraction {}", j + 1), Scalar::Num(layer_fraction)));
    scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(f.re)));
    scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(f.im)));
    if has_interior && any_pml {
        let two_pi = 2.0 * std::f64::consts::PI;
        let p_pml = (two_pi * f.re) * (two_pi * f.im) * phys_mass_norm;
        scalars.push((format!("P_pml {} (W)", j + 1), Scalar::Num(p_pml)));
    }
    (label, f.re)
}

fn finish_pml_log(
    asm: &PmlAssembly,
    solved: &PmlSolve,
    unstable_count: usize,
    scalars: &mut Vec<(String, Scalar)>,
) -> Vec<String> {
    let mut log = asm.diagnostics.clone();
    log.push(if let Some(red) = &asm.red {
        format!(
            "{} Floquet-periodic pair(s) with PML: {} dofs, {} after constraints; complex eigenproblem, Krylov dimension {}{}",
            asm.n_pairs,
            asm.n,
            red.n_reduced(),
            solved.dim,
            if solved.converged { "" } else { " (not every pair converged)" }
        )
    } else {
        format!(
            "complex eigenproblem with PML, Krylov dimension {}{}",
            solved.dim,
            if solved.converged { "" } else { " (not fully converged)" }
        )
    });
    if !asm.has_fixed {
        log.push("no fixed boundary: the lowest modes are rigid-body motions".into());
    }
    if unstable_count > 0 || (asm.has_interior && asm.any_pml) {
        scalars.push(("unstable_modes".into(), Scalar::Num(unstable_count as f64)));
    }
    if unstable_count > 0 {
        log.push(format!(
            "WARNING: {unstable_count} amplifying backward-wave mode(s) detected (Im f < -1e-9 Re f); enable stabilization = \"auto\", \"mpml\", or \"adiabatic\""
        ));
    }
    log
}

/// Stage 3: plot fields and per-cell complex displacement of every mode, frequencies, signed Q,
/// `layer fraction N`, and radiated elastic power `P_pml N (W)`.
pub fn pml_fields(mesh: &Discretisation, asm: &PmlAssembly, solved: &PmlSolve) -> Result<SolveResult, String> {
    let (order, d, n) = (asm.order, asm.d, asm.n);
    let n1 = n / d;
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut unstable_count = 0usize;

    for (j, &(i, raw_f)) in solved.freqs.iter().enumerate() {
        let u = if let Some(red) = &asm.red {
            red.expand(&solved.vectors[i])
        } else {
            let mut u = vec![c(0.0, 0.0); n];
            for (a, &dof) in asm.free.iter().enumerate() {
                u[dof] = solved.vectors[i][a];
            }
            u
        };
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let scale2 = 1.0 / big.norm_sqr().max(1e-300);

        let (layer_fraction, phys_mass_norm) = match (&asm.m1_phys_all, &asm.m1_phys_pml) {
            (Some(ma), Some(mp)) => {
                let tot = real_quad_norm(ma, &u, n1, d) * scale2;
                let in_pml = real_quad_norm(mp, &u, n1, d) * scale2;
                let lf = if tot > 0.0 { (in_pml / tot).clamp(0.0, 1.0) } else { 0.0 };
                (lf, (tot - in_pml).max(0.0))
            }
            _ => (0.0, 0.0),
        };

        let (fields, cells) = crate::solid::complex_mode_fields(mesh, order, &u);
        let (label, f_re) = push_pml_mode_scalars_and_label(
            j,
            raw_f,
            layer_fraction,
            phys_mass_norm,
            asm.any_pml,
            asm.has_interior,
            &mut scalars,
            &mut unstable_count,
        );
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(f_re),
            ..Default::default()
        });
    }
    let log = finish_pml_log(asm, solved, unstable_count, &mut scalars);
    Ok(SolveResult { order: order as u8, dim: d as u8, solutions, scalars, log, dofs: n, ..Default::default() })
}

fn pml_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    proj: &MixedToSimplexProjector,
    asm: &PmlAssembly,
    solved: &PmlSolve,
) -> Result<SolveResult, String> {
    let (order, d, n) = (asm.order, asm.d, asm.n);
    let n1 = dof_map.n_dofs;
    let mut solutions = Vec::new();
    let mut scalars = Vec::new();
    let mut unstable_count = 0usize;

    for (j, &(i, raw_f)) in solved.freqs.iter().enumerate() {
        let u = if let Some(red) = &asm.red {
            red.expand(&solved.vectors[i])
        } else {
            let mut u = vec![c(0.0, 0.0); n];
            for (a, &dof) in asm.free.iter().enumerate() {
                u[dof] = solved.vectors[i][a];
            }
            u
        };
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let scale2 = 1.0 / big.norm_sqr().max(1e-300);

        let (layer_fraction, phys_mass_norm) = match (&asm.m1_phys_all, &asm.m1_phys_pml) {
            (Some(ma), Some(mp)) => {
                let tot = real_quad_norm(ma, &u, n1, d) * scale2;
                let in_pml = real_quad_norm(mp, &u, n1, d) * scale2;
                let lf = if tot > 0.0 { (in_pml / tot).clamp(0.0, 1.0) } else { 0.0 };
                (lf, (tot - in_pml).max(0.0))
            }
            _ => (0.0, 0.0),
        };

        let (fields, cell_fields) = crate::solid::complex_mode_fields_mixed(mixed, dof_map, proj, order, &u);
        let (label, f_re) = push_pml_mode_scalars_and_label(
            j,
            raw_f,
            layer_fraction,
            phys_mass_norm,
            asm.any_pml,
            asm.has_interior,
            &mut scalars,
            &mut unstable_count,
        );
        solutions.push(Solution { label, fields, cell_fields, value: Some(f_re), ..Default::default() });
    }
    let log = finish_pml_log(asm, solved, unstable_count, &mut scalars);
    Ok(SolveResult { order: order as u8, dim: d as u8, solutions, scalars, log, dofs: n, ..Default::default() })
}
