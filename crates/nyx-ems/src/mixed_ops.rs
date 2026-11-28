// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared quadrature-point curvilinear PML pullback (`EmPmlTable`), symmetric `3x3` `C3` permeability
//! contractions, native mixed-cell (`MixedDiscretisation`: `Tri`+`Quad` in 2D, `Tet`+`Prism`+`Hex` in 3D)
//! Maxwell / Floquet-Bloch / driven / sheet / port assemblers, and `Ams` (`P1`) / `AmsP2` (`P2`)
//! auxiliary-space iterative curl-curl solvers.

use crate::cavity::{assemble_complex, CURVED_QUAD_DEGREE};
use crate::ibc::{
    ibc_parallel_params, jj_inv_l_tot, sheet2d_tangent_frame_3d, sorted3, AssembledJjPort, IbcJjAssembly3D,
    PerturbativeIbcTri, SheetCellKind, ETA0, MU0,
};
use crate::modes::{assemble_modes_with_bcs, solve_mode_pairs, EpsCell, C0};
use crate::ports::{AssembledLumpedPort, AssembledWavePort};
pub use crate::problem::symmetry_facet_flags_mixed;
use crate::problem::{LumpedJjBc, LumpedPortBc, Problem, Sheet2dBc, SurfaceImpedanceBc};
use femcore::edge::{cell_data, edge_basis, tri_g};
use femcore::edge_mixed::{
    build_ams_p1_from_masks, build_ams_p2_from_masks, eval_cell_basis_3d, gradient_p1_full, gradient_p2_full,
    nedelec_nloc, scalar_p2_tensor_count, NedEvalScratch, NedMixed,
};
use femcore::fem::{geometry, quad_rule};
use femcore::fem_mixed::{
    cell_corner_ref_coords, eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, quad_rule_mixed,
    IsoparametricCell, LagrangeSpace,
};
use femcore::linalg::ams::{Ams, AmsOptions};
use femcore::linalg::csr::{c, CCsr, Coo, RCsr, C};
use femcore::linalg::driver::{
    eigs_curl_curl_p2_jd, solve_curl_curl_p2_iterative, CEigPairs, IterativeCurlCurlOptions, IterativeSolveReport,
};
use femcore::linalg::ksp::{Gmres, GmresOptions, Orthogonalization, PcSide, Tolerances};
use femcore::material::{self, M3};
use femcore::mesh::cells::CellType;
use femcore::mesh::mixed::{
    FacetType, MixedDiscretisation, HEX_QUAD_FACES, PRISM_QUAD_FACES, PRISM_TRI_FACES, QUAD_EDGES, TET_TRI_FACES,
    TRI_EDGES,
};
use femcore::mesh::SimplexMesh;
use femcore::periodic::phase_factor;
use femcore::pml::{auto_scale_params, pml_jacobian, pullback_maxwell, PmlJacobian};
use femcore::props::{Materials, Periodic, Pml, PmlCoord, PmlDir, PmlProfile};
use femcore::section::wafer_point;
use nyx_physics::Discretisation;
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::f64::consts::PI;

/// Symmetric `3x3` complex material tensor (`eps` or `mu_inv`).
pub type C3 = [[C; 3]; 3];

/// Converts an isotropic scalar permeability `mu_r` into a diagonal `3x3` `C3` inverse permeability `diag(1/mu_r)`.
#[inline]
pub fn iso_mu_inv_c3(mur: f64) -> C3 {
    let v = c(1.0 / mur, 0.0);
    let z = c(0.0, 0.0);
    [[v, z, z], [z, v, z], [z, z, v]]
}

/// Converts a real+imaginary loss pair `(re, im_loss)` into complex `eps = re - i * im_loss`.
#[inline]
pub fn complex_eps_c3(re: &M3, im_loss: &M3) -> C3 {
    let mut out = [[c(0.0, 0.0); 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            out[a][b] = c(re[a][b], -im_loss[a][b]);
        }
    }
    out
}

/// Bilinear contraction `\sum_{a,b=0..3} m[a][b] * u[a] * v[b]` for real basis vectors `u, v`.
#[inline]
pub fn c3_dot_rr(m: &C3, u: [f64; 3], v: [f64; 3]) -> C {
    let mut acc = c(0.0, 0.0);
    for a in 0..3 {
        for b in 0..3 {
            acc += m[a][b] * (u[a] * v[b]);
        }
    }
    acc
}

/// Resolved per-domain PML specification with inferred geometric bounds and optional wavelength auto-scaling.
#[derive(Clone, Debug)]
pub struct ResolvedDomainPml {
    pub spec: Pml,
    pub fallback_bounds: ([f64; 3], [f64; 3]),
}

/// Per-cell lookup table for quadrature-point curvilinear PML Jacobian evaluation and Maxwell tensor pullback.
#[derive(Clone, Debug)]
pub struct EmPmlTable {
    domains: BTreeMap<i32, ResolvedDomainPml>,
    cell_pml: Vec<Option<ResolvedDomainPml>>,
}

impl EmPmlTable {
    /// Builds a resolved PML table for a `SimplexMesh` (`Discretisation`).
    pub fn from_simplex(
        mesh: &Discretisation,
        materials: &Materials,
        f_ref: f64,
        is_eigen: bool,
        is_iterative: bool,
    ) -> Self {
        let mut dom_pts: BTreeMap<i32, Vec<[f64; 3]>> = BTreeMap::new();
        for ci in 0..mesh.n_cells() {
            let dom = mesh.cell_labels[ci];
            if materials.pml.contains_key(&dom) {
                let entry = dom_pts.entry(dom).or_default();
                for &v in mesh.cell(ci) {
                    entry.push(wafer_point(mesh.point(v as usize), mesh.dim));
                }
            }
        }
        Self::build_from_domain_points(
            mesh.dim,
            &mesh.cell_labels,
            &mesh.domain_labels(),
            dom_pts,
            materials,
            f_ref,
            is_eigen,
            is_iterative,
        )
    }

    /// Builds a resolved PML table for a `MixedDiscretisation`.
    pub fn from_mixed(
        mesh: &MixedDiscretisation,
        materials: &Materials,
        f_ref: f64,
        is_eigen: bool,
        is_iterative: bool,
    ) -> Self {
        let mut dom_pts: BTreeMap<i32, Vec<[f64; 3]>> = BTreeMap::new();
        for ci in 0..mesh.n_cells() {
            let dom = mesh.cell_labels[ci];
            if materials.pml.contains_key(&dom) {
                let entry = dom_pts.entry(dom).or_default();
                for &v in mesh.cell(ci) {
                    entry.push(wafer_point(mesh.point(v as usize), mesh.dim));
                }
            }
        }
        Self::build_from_domain_points(
            mesh.dim,
            &mesh.cell_labels,
            &mesh.domain_labels(),
            dom_pts,
            materials,
            f_ref,
            is_eigen,
            is_iterative,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build_from_domain_points(
        dim: usize,
        cell_labels: &[i32],
        all_domains: &[i32],
        dom_pts: BTreeMap<i32, Vec<[f64; 3]>>,
        materials: &Materials,
        f_ref: f64,
        is_eigen: bool,
        is_iterative: bool,
    ) -> Self {
        let (lam_min_all, lam_max_all) = estimate_wavelength_range(all_domains, materials, f_ref);
        let mut domains = BTreeMap::new();

        for (&dom, raw_spec) in &materials.pml {
            let pts = dom_pts.get(&dom).map(Vec::as_slice).unwrap_or(&[]);
            if pts.is_empty() {
                continue;
            }
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for p in pts {
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            let mut spec = raw_spec.clone();

            // Infer missing radial/axial/spherical shell extents from the PML domain's actual vertices when 0.
            match &mut spec.coord {
                PmlCoord::Cartesian => {}
                PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness } => {
                    let an = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
                    let a_hat = if an > 1e-30 {
                        [axis[0] / an, axis[1] / an, axis[2] / an]
                    } else if dim == 2 {
                        [0.0, 1.0, 0.0]
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    let mut r_lo = f64::INFINITY;
                    let mut r_hi = 0.0_f64;
                    let mut z_lo = f64::INFINITY;
                    let mut z_hi = f64::NEG_INFINITY;
                    for p in pts {
                        let dv = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                        let zeta = dv[0] * a_hat[0] + dv[1] * a_hat[1] + dv[2] * a_hat[2];
                        let rv = [dv[0] - zeta * a_hat[0], dv[1] - zeta * a_hat[1], dv[2] - zeta * a_hat[2]];
                        let r = (rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]).sqrt();
                        r_lo = r_lo.min(r);
                        r_hi = r_hi.max(r);
                        z_lo = z_lo.min(zeta);
                        z_hi = z_hi.max(zeta);
                    }
                    if *radial_on && *r_thickness <= 0.0 && r_hi > r_lo {
                        if *r_in <= 0.0 {
                            *r_in = r_lo;
                        }
                        *r_thickness = (r_hi - *r_in).max(0.0);
                    }
                    if *axial != PmlDir::Off && *z_thickness <= 0.0 && z_hi > z_lo {
                        if *z_in == 0.0 {
                            *z_in = match *axial {
                                PmlDir::Plus => z_lo,
                                PmlDir::Minus => z_hi,
                                PmlDir::Off => 0.0,
                            };
                        }
                        *z_thickness = match *axial {
                            PmlDir::Plus => (z_hi - *z_in).max(0.0),
                            PmlDir::Minus => (*z_in - z_lo).max(0.0),
                            PmlDir::Off => 0.0,
                        };
                    }
                }
                PmlCoord::Spherical { center, r_in, thickness } => {
                    if *thickness <= 0.0 {
                        let mut r_lo = f64::INFINITY;
                        let mut r_hi = 0.0_f64;
                        for p in pts {
                            let dv = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
                            let r = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
                            r_lo = r_lo.min(r);
                            r_hi = r_hi.max(r);
                        }
                        if r_hi > r_lo {
                            if *r_in <= 0.0 {
                                *r_in = r_lo;
                            }
                            *thickness = (r_hi - *r_in).max(0.0);
                        }
                    }
                }
            }

            if !spec.is_active() {
                continue;
            }

            if is_iterative && matches!(spec.profile, PmlProfile::Rational { .. }) {
                femcore::progress::note(format!(
                    "domain {dom}: rational PML profile replaced with kappa_polynomial (kappa_max = 2.0) for iterative AMS preconditioner conditioning"
                ));
                spec.profile = PmlProfile::KappaPolynomial { kappa_max: 2.0 };
            }

            if spec.auto_scale && f_ref > 0.0 && f_ref.is_finite() {
                let d_actual = active_pml_thickness(&spec, (lo, hi));
                let is_rat = matches!(spec.profile, PmlProfile::Rational { .. }) && !is_eigen && !is_iterative;
                let auto = auto_scale_params(
                    lam_min_all,
                    lam_max_all,
                    d_actual,
                    spec.order,
                    spec.target_reflection,
                    is_eigen,
                    is_rat,
                );
                spec.strength = auto.psi_max;
                spec.order = auto.order;
                spec.profile = if is_rat {
                    PmlProfile::Rational { curvature: auto.rational_curvature() }
                } else {
                    PmlProfile::KappaPolynomial { kappa_max: auto.kappa_max }
                };
            }

            domains.insert(dom, ResolvedDomainPml { spec, fallback_bounds: (lo, hi) });
        }

        let cell_pml = cell_labels.iter().map(|l| domains.get(l).cloned()).collect();
        Self { domains, cell_pml }
    }

    #[inline]
    pub fn has_any(&self) -> bool {
        !self.domains.is_empty()
    }

    #[inline]
    pub fn is_pml(&self, ci: usize) -> bool {
        self.cell_pml[ci].is_some()
    }

    pub fn pml_cells_mask(&self) -> Vec<bool> {
        self.cell_pml.iter().map(Option::is_some).collect()
    }

    /// Evaluates pulled-back `(eps_tilde, mu_inv_tilde)` in 3D wafer axes `(x, y, z)` at `x_wafer`.
    #[inline]
    pub fn eval_3d(&self, ci: usize, eps_0: C3, mu_inv_0: C3, x_wafer: [f64; 3]) -> (C3, C3) {
        match &self.cell_pml[ci] {
            Some(rp) => {
                let jac = pml_jacobian(&rp.spec, Some(rp.fallback_bounds), x_wafer, 3, false);
                pullback_maxwell(eps_0, mu_inv_0, &jac)
            }
            None => (eps_0, mu_inv_0),
        }
    }

    /// Evaluates `(eps_tilde, mu_inv_tilde, jac)` in 2D cross-section/axisymmetric local axes `(x/r, z, y/theta)`
    /// at 2D mesh point `x_2d = [x0, x1]`.
    #[inline]
    pub fn eval_2d(
        &self,
        ci: usize,
        eps_local_0: C3,
        mu_inv_local_0: C3,
        x_2d: [f64; 2],
        axisym: bool,
    ) -> (C3, C3, PmlJacobian) {
        match &self.cell_pml[ci] {
            Some(rp) => {
                let jac = pml_jacobian(&rp.spec, Some(rp.fallback_bounds), [x_2d[0], 0.0, x_2d[1]], 2, axisym);
                let (e, m) = pullback_maxwell(eps_local_0, mu_inv_local_0, &jac);
                (e, m, jac)
            }
            None => (eps_local_0, mu_inv_local_0, PmlJacobian::identity()),
        }
    }
}

fn active_pml_thickness(spec: &Pml, bounds: ([f64; 3], [f64; 3])) -> Option<f64> {
    let (lo, hi) = spec.bounds.unwrap_or(bounds);
    match spec.coord {
        PmlCoord::Cartesian => {
            let mut d = f64::INFINITY;
            for (a, dir) in [spec.x, spec.y, spec.z].into_iter().enumerate() {
                if dir != PmlDir::Off {
                    let span = hi[a] - lo[a];
                    if span > 0.0 {
                        d = d.min(span);
                    }
                }
            }
            d.is_finite().then_some(d)
        }
        PmlCoord::Cylindrical { r_thickness, radial_on, axial, z_thickness, .. } => {
            let mut d = f64::INFINITY;
            if radial_on && r_thickness > 0.0 {
                d = d.min(r_thickness);
            }
            if axial != PmlDir::Off && z_thickness > 0.0 {
                d = d.min(z_thickness);
            }
            d.is_finite().then_some(d)
        }
        PmlCoord::Spherical { thickness, .. } => (thickness > 0.0).then_some(thickness),
    }
}

fn estimate_wavelength_range(domains: &[i32], materials: &Materials, f_ref: f64) -> (f64, f64) {
    if !(f_ref > 0.0 && f_ref.is_finite()) {
        return (1.0, 1.0);
    }
    let lam0 = C0 / f_ref;
    let mut lam_min = f64::INFINITY;
    let mut lam_max = 0.0_f64;
    for &d in domains {
        let mur = materials.get(d, "mu_r").unwrap_or(1.0).max(1e-6);
        let (eps_min, eps_max) = match materials.aniso.get(&d).and_then(|a| a.permittivity.as_ref()) {
            Some(t) => {
                let g = material::tensor_global(t, materials.aniso[&d].euler_deg);
                let d0 = g[0][0].abs().max(1e-6);
                let d1 = g[1][1].abs().max(1e-6);
                let d2 = g[2][2].abs().max(1e-6);
                (d0.min(d1).min(d2), d0.max(d1).max(d2))
            }
            None => {
                let e = materials.get(d, "eps_r").unwrap_or(1.0).abs().max(1e-6);
                (e, e)
            }
        };
        let l_short = lam0 / (eps_max * mur).sqrt();
        let l_long = lam0 / (eps_min * mur).sqrt();
        lam_min = lam_min.min(l_short);
        lam_max = lam_max.max(l_long);
    }
    if !lam_min.is_finite() || lam_max <= 0.0 {
        (lam0, lam0)
    } else {
        (lam_min, lam_max)
    }
}

// ------------------------------------------------------------------ MixedDiscretisation material lookup

/// Returns `true` if a `3x3` complex tensor `m` violates transpose symmetry (`m[a][b] != m[b][a]`).
#[inline]
pub fn is_c3_non_symmetric(m: &C3) -> bool {
    let max_v = m.iter().flatten().fold(0.0_f64, |acc, z| acc.max(z.norm()));
    let tol = 1e-14 * max_v.max(1e-30);
    (m[0][1] - m[1][0]).norm() > tol || (m[0][2] - m[2][0]).norm() > tol || (m[1][2] - m[2][1]).norm() > tol
}

/// Returns `true` if any cell in `eps` or `mu_inv` is non-symmetric.
pub fn is_c3_slice_non_symmetric(eps: &[C3], mu_inv: &[C3]) -> bool {
    eps.iter().any(is_c3_non_symmetric) || mu_inv.iter().any(is_c3_non_symmetric)
}

/// Multiplies a `3x3` complex tensor `m` by a 3-vector `v`.
#[inline]
pub fn c3_mul_vec(m: &C3, v: [C; 3]) -> [C; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

/// Unstretched per-cell electromagnetic tensors on a `MixedDiscretisation`.
pub struct MixedCellMaterials {
    pub eps_0: Vec<C3>,
    pub mu_inv_0: Vec<C3>,
    pub eps_phys_re: Vec<M3>,
    pub eps_phys_im: Vec<M3>,
    pub d_omega_eps: Vec<M3>,
    pub mur_phys: Vec<f64>,
    pub is_non_symmetric: bool,
}

impl MixedCellMaterials {
    pub fn from_problem(mesh: &MixedDiscretisation, problem: &Problem) -> Result<Self, String> {
        let doms = mesh.domain_labels();
        let mut re_by_dom = BTreeMap::new();
        let mut im_by_dom = BTreeMap::new();
        let mut mu_by_dom = BTreeMap::new();
        let mut eps_c_by_dom = BTreeMap::new();
        let mut mui_c_by_dom = BTreeMap::new();
        let has_tensor_mu = problem.materials.has_tensor(femcore::props::TensorKind::Permeability)
            || problem.materials.has_tensor(femcore::props::TensorKind::PermeabilityIm)
            || problem.materials.aniso.values().any(|a| a.mu_gyrotropy.is_some() || a.mu_gyrotropy_im.is_some());

        for &d in &doms {
            let re = match problem.materials.aniso.get(&d).and_then(|a| a.permittivity.as_ref()) {
                Some(t) => material::tensor_global(t, problem.materials.aniso[&d].euler_deg),
                None => {
                    let v = problem.materials.get(d, "eps_r").unwrap_or(1.0);
                    material::iso_tensor(v)
                }
            };
            let im = match problem.materials.aniso.get(&d).and_then(|a| a.permittivity_im.as_ref()) {
                Some(t) => material::tensor_global(t, problem.materials.aniso[&d].euler_deg),
                None => {
                    let v = problem.materials.get(d, "eps_r_im").unwrap_or(0.0);
                    material::iso_tensor(v)
                }
            };
            let mur = problem.materials.get(d, "mu_r").unwrap_or(1.0);
            let eps_c = problem.materials.domain_full_complex_permittivity(d, Some(1.0))?;
            let mui_c = if has_tensor_mu {
                let mu_c = problem.materials.domain_full_complex_permeability(d, Some(1.0))?;
                femcore::props::inv3_c(&mu_c).map_err(|e| format!("domain {d}: {e}"))?
            } else {
                iso_mu_inv_c3(mur)
            };
            re_by_dom.insert(d, re);
            im_by_dom.insert(d, im);
            mu_by_dom.insert(d, mur);
            eps_c_by_dom.insert(d, eps_c);
            mui_c_by_dom.insert(d, mui_c);
        }

        let nc = mesh.n_cells();
        let mut eps_0 = Vec::with_capacity(nc);
        let mut mu_inv_0 = Vec::with_capacity(nc);
        let mut eps_phys_re = Vec::with_capacity(nc);
        let mut eps_phys_im = Vec::with_capacity(nc);
        let mut d_omega_eps = Vec::with_capacity(nc);
        let mut mur_phys = Vec::with_capacity(nc);

        let nv = mesh.n_points();
        let ct_opt = problem.cell_temperatures.as_deref();
        let tf_opt = problem.temperature_field.as_deref();
        let eo_pe = if problem.has_eo_or_pe_perturbation() { Some(problem.domain_eo_pe_global()) } else { None };
        let cell_e = if eo_pe.is_some() {
            problem.eval_cell_electric_field_with(nc, nv, mesh.dim, |ci| mesh.cell(ci))
        } else {
            None
        };
        let cell_s = if eo_pe.is_some() { problem.eval_cell_strain_with(nc, nv, |ci| mesh.cell(ci)) } else { None };
        for ci in 0..nc {
            let dom = mesh.cell_labels[ci];
            let (mut re, mut im) = (re_by_dom[&dom], im_by_dom[&dom]);
            let mut eps_c = eps_c_by_dom[&dom];
            if let Some(td) = problem.thermo_optic_domains.get(&dom) {
                let tc_opt = if let Some(ct) = ct_opt.filter(|v| v.len() == nc) {
                    Some(ct[ci])
                } else if let Some(tf) = tf_opt {
                    if tf.len() == nc {
                        Some(tf[ci])
                    } else if tf.len() >= nv {
                        let cv = mesh.cell(ci);
                        Some(cv.iter().map(|&v| tf[v as usize]).sum::<f64>() / (cv.len() as f64))
                    } else {
                        None
                    }
                } else {
                    None
                };
                if let Some(tc) = tc_opt.filter(|&t| t > 0.0 && t.is_finite()) {
                    if let Some(pair) = td.eval_at(tc)? {
                        re = pair.0;
                        im = pair.1;
                        eps_c = apply_thermo_optic_to_c3(eps_c, &re, &im);
                    }
                }
            }
            if let Some(ref map) = eo_pe {
                if let Some((r_g, p_g)) = map.get(&dom) {
                    let e_dc = cell_e.as_ref().map(|v| &v[ci]);
                    let s_dc = cell_s.as_ref().map(|v| &v[ci]);
                    let d_eps =
                        femcore::material::coupled::delta_eps_eo_pe(&re, r_g.as_ref(), p_g.as_ref(), e_dc, s_dc);
                    for a in 0..3 {
                        for b in 0..3 {
                            re[a][b] += d_eps[a][b];
                        }
                    }
                    eps_c = apply_thermo_optic_to_c3(eps_c, &re, &im);
                }
            }
            let mur = mu_by_dom[&dom];
            eps_0.push(eps_c);
            mu_inv_0.push(mui_c_by_dom[&dom]);
            eps_phys_re.push(re);
            eps_phys_im.push(im);
            d_omega_eps.push(problem.d_omega_eps.get(&dom).copied().unwrap_or(re));
            mur_phys.push(mur);
        }
        let is_non_symmetric = problem.materials.is_non_symmetric() || is_c3_slice_non_symmetric(&eps_0, &mu_inv_0);

        Ok(Self { eps_0, mu_inv_0, eps_phys_re, eps_phys_im, d_omega_eps, mur_phys, is_non_symmetric })
    }
}

/// Replaces the symmetric part of `base_eps` with `complex_eps_c3(re, im)` while preserving any
/// skew-symmetric gyrotropic off-diagonal entries `0.5 * (base_eps[a][b] - base_eps[b][a])`.
#[inline]
pub(crate) fn apply_thermo_optic_to_c3(base_eps: C3, re: &femcore::material::M3, im: &femcore::material::M3) -> C3 {
    let mut out = complex_eps_c3(re, im);
    let mut out_has_skew = false;
    for a in 0..3 {
        for b in (a + 1)..3 {
            if (out[a][b] - out[b][a]).norm() > 1e-15 {
                out_has_skew = true;
            }
        }
    }
    if !out_has_skew {
        for a in 0..3 {
            for b in 0..3 {
                out[a][b] += c(0.5, 0.0) * (base_eps[a][b] - base_eps[b][a]);
            }
        }
    }
    out
}

/// Evaluates unstretched `(eps_0, mu_inv_0)` on a `SimplexMesh` (`Discretisation`) in 3D wafer axes.
pub fn simplex_cell_materials_3d(mesh: &Discretisation, problem: &Problem) -> Result<(Vec<C3>, Vec<C3>), String> {
    let mut eps_0 = problem.materials.cell_full_complex_permittivity(mesh, Some(1.0))?;
    if problem.has_cell_material_perturbation() {
        let tens = problem.cell_complex_permittivity(mesh)?;
        let eo_pe = problem.domain_eo_pe_global();
        for ci in 0..mesh.n_cells() {
            let dom = mesh.cell_labels[ci];
            if problem.thermo_optic_domains.contains_key(&dom) || eo_pe.contains_key(&dom) {
                eps_0[ci] = apply_thermo_optic_to_c3(eps_0[ci], &tens[ci].0, &tens[ci].1);
            }
        }
    }
    let has_tensor_mu = problem.materials.has_tensor(femcore::props::TensorKind::Permeability)
        || problem.materials.has_tensor(femcore::props::TensorKind::PermeabilityIm)
        || problem.materials.aniso.values().any(|a| a.mu_gyrotropy.is_some() || a.mu_gyrotropy_im.is_some());
    let mu_inv_0 = if has_tensor_mu {
        problem.materials.cell_full_complex_permeability_inv(mesh, Some(1.0))?
    } else {
        let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
        mur.into_iter().map(iso_mu_inv_c3).collect()
    };
    Ok((eps_0, mu_inv_0))
}

// ------------------------------------------------------------------ 3D Mixed-cell Maxwell & Bloch assembly

#[inline]
pub fn add_maxwell_point_c3(
    which: usize,
    nb: usize,
    ke: &mut [C],
    w: &[[f64; 3]],
    cu: &[[f64; 3]],
    wq: f64,
    eps: &C3,
    mu_inv: &C3,
) {
    if which == 0 {
        for i in 0..nb {
            let ci_ = cu[i];
            for j in 0..nb {
                ke[i * nb + j] += c3_dot_rr(mu_inv, ci_, cu[j]) * wq;
            }
        }
    } else {
        for i in 0..nb {
            let wi = w[i];
            for j in 0..nb {
                ke[i * nb + j] += c3_dot_rr(eps, wi, w[j]) * wq;
            }
        }
    }
}

/// Assembles `(K, M)` on a 3D `MixedDiscretisation` (`Tet + Prism + Hex`) with quadrature-point PML pullback.
pub fn assemble_maxwell_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    eps_0: &[C3],
    mu_inv_0: &[C3],
    pml: Option<&EmPmlTable>,
) -> (CCsr, CCsr) {
    let nc = mesh.n_cells();
    let n = ned.n_dofs;
    let run = |which: usize| -> CCsr {
        assemble_complex(nc, n, |range, coo| {
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            let mut ke = Vec::new();
            for ci in range {
                let iso = IsoparametricCell::from_mesh(mesh, ci);
                let is_pml = pml.is_some_and(|p| p.is_pml(ci));
                let deg = if iso.curved || is_pml { CURVED_QUAD_DEGREE } else { 4 };
                let rule = quad_rule_mixed(iso.cell_type, deg);
                let nb = nedelec_nloc(iso.cell_type, ned.order);
                ke.clear();
                ke.resize(nb * nb, c(0.0, 0.0));
                for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_3d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let wq = w_ref * jac.det.abs();
                    let (ep, mui) = match pml {
                        Some(p) => p.eval_3d(ci, eps_0[ci], mu_inv_0[ci], jac.x),
                        None => (eps_0[ci], mu_inv_0[ci]),
                    };
                    add_maxwell_point_c3(which, nb, &mut ke, &w, &cu, wq, &ep, &mui);
                }
                let cdofs = ned.cell_dofs(ci);
                coo.push_block(cdofs, cdofs, &ke);
            }
        })
    };
    (run(0), run(1))
}

#[inline]
#[allow(clippy::too_many_arguments)]
pub fn add_bloch_point_c3(
    which: usize,
    nb: usize,
    ke: &mut [C],
    w: &[[f64; 3]],
    cu: &[[f64; 3]],
    wq: f64,
    inv_k0: f64,
    inv_k0_2: f64,
    eps: &C3,
    mu_inv: &C3,
) {
    match which {
        0 => {
            for i in 0..nb {
                let ci_ = cu[i];
                let wi = w[i];
                for j in 0..nb {
                    let k_curl = c3_dot_rr(mu_inv, ci_, cu[j]);
                    let m_eps = c3_dot_rr(eps, wi, w[j]);
                    ke[i * nb + j] += (k_curl * inv_k0_2 - m_eps) * wq;
                }
            }
        }
        1 => {
            // u_i = y_hat x w_i = [w_z, 0, -w_x]
            let pref = c(0.0, inv_k0 * wq);
            for i in 0..nb {
                let ui = [w[i][2], 0.0, -w[i][0]];
                let ci_ = cu[i];
                for j in 0..nb {
                    let uj = [w[j][2], 0.0, -w[j][0]];
                    let cj = cu[j];
                    let term = c3_dot_rr(mu_inv, ui, cj) - c3_dot_rr(mu_inv, ci_, uj);
                    ke[i * nb + j] += pref * term;
                }
            }
        }
        _ => {
            for i in 0..nb {
                let ui = [w[i][2], 0.0, -w[i][0]];
                for j in 0..nb {
                    let uj = [w[j][2], 0.0, -w[j][0]];
                    ke[i * nb + j] += c3_dot_rr(mu_inv, ui, uj) * wq;
                }
            }
        }
    }
}

/// Assembles `(A0, A1, A2)` for 3D Floquet-Bloch modes on a `MixedDiscretisation`.
pub fn assemble_bloch_matrices_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    k0: f64,
    eps_0: &[C3],
    mu_inv_0: &[C3],
    pml: Option<&EmPmlTable>,
) -> (CCsr, CCsr, CCsr) {
    let nc = mesh.n_cells();
    let n = ned.n_dofs;
    let inv_k0 = 1.0 / k0;
    let inv_k0_2 = inv_k0 * inv_k0;
    let run = |which: usize| -> CCsr {
        assemble_complex(nc, n, |range, coo| {
            let mut scratch = NedEvalScratch::new();
            let mut w = Vec::new();
            let mut cu = Vec::new();
            let mut ke = Vec::new();
            for ci in range {
                let iso = IsoparametricCell::from_mesh(mesh, ci);
                let is_pml = pml.is_some_and(|p| p.is_pml(ci));
                let deg = if iso.curved || is_pml { CURVED_QUAD_DEGREE } else { 4 };
                let rule = quad_rule_mixed(iso.cell_type, deg);
                let nb = nedelec_nloc(iso.cell_type, ned.order);
                ke.clear();
                ke.resize(nb * nb, c(0.0, 0.0));
                for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
                    let jac = eval_cell_basis_3d(mesh, ned, ci, &iso, xi, &mut scratch, &mut w, &mut cu);
                    let wq = w_ref * jac.det.abs();
                    let (ep, mui) = match pml {
                        Some(p) => p.eval_3d(ci, eps_0[ci], mu_inv_0[ci], jac.x),
                        None => (eps_0[ci], mu_inv_0[ci]),
                    };
                    add_bloch_point_c3(which, nb, &mut ke, &w, &cu, wq, inv_k0, inv_k0_2, &ep, &mui);
                }
                let cdofs = ned.cell_dofs(ci);
                coo.push_block(cdofs, cdofs, &ke);
            }
        })
    };
    (run(0), run(1), run(2))
}

// ------------------------------------------------------------------ Floquet periodicity on MixedDiscretisation

/// Matched periodic boundary pair on a `MixedDiscretisation`.
#[derive(Clone, Debug)]
pub struct MixedPeriodicMatch {
    pub map: HashMap<u32, u32>,
    pub src_facets: Vec<usize>,
    pub shift: [f64; 3],
}

pub fn match_vertices_mixed(mesh: &MixedDiscretisation, p: &Periodic) -> Result<MixedPeriodicMatch, String> {
    let collect = |labels: &[i32]| -> (Vec<usize>, Vec<u32>) {
        let mut fcs = Vec::new();
        let mut vs = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if labels.contains(&mesh.facet_labels[f]) {
                fcs.push(f);
                for &v in mesh.facet(f) {
                    vs.insert(v);
                }
            }
        }
        (fcs, vs.into_iter().collect())
    };
    let (src_facets, sv) = collect(&p.source);
    let (_, dv) = collect(&p.dest);
    if sv.is_empty() || dv.is_empty() {
        return Err(format!("periodic boundary: labels {:?} or {:?} select no facets", p.source, p.dest));
    }
    if sv.len() != dv.len() {
        return Err(format!(
            "periodic boundary: source ({} vertices) and destination ({} vertices) do not match",
            sv.len(),
            dv.len()
        ));
    }
    let centroid = |vs: &[u32]| -> [f64; 3] {
        let mut c = [0.0; 3];
        for &v in vs {
            let pt = mesh.point3(v as usize);
            for k in 0..mesh.dim {
                c[k] += pt[k];
            }
        }
        let inv = 1.0 / vs.len() as f64;
        [c[0] * inv, c[1] * inv, c[2] * inv]
    };
    let (cs, cd) = (centroid(&sv), centroid(&dv));
    let shift = [cd[0] - cs[0], cd[1] - cs[1], cd[2] - cs[2]];
    let (lo, hi) = mesh.bbox();
    let span = (0..mesh.dim).map(|k| hi[k] - lo[k]).fold(0.0_f64, f64::max).max(1e-30);
    let tol2 = (1e-6 * span).powi(2);

    let mut map = HashMap::with_capacity(sv.len());
    let mut used = BTreeSet::new();
    for &s in &sv {
        let ps = mesh.point3(s as usize);
        let want = [ps[0] + shift[0], ps[1] + shift[1], ps[2] + shift[2]];
        let mut best: Option<(u32, f64)> = None;
        for &t in &dv {
            let pt = mesh.point3(t as usize);
            let d2 = (0..mesh.dim).map(|k| (pt[k] - want[k]).powi(2)).sum::<f64>();
            if best.is_none_or(|(_, bd)| d2 < bd) {
                best = Some((t, d2));
            }
        }
        let (t, d2) = best.unwrap();
        if d2 > tol2 || !used.insert(t) {
            return Err("periodic boundary: source vertex has no unique translated counterpart on destination".into());
        }
        map.insert(s, t);
    }
    Ok(MixedPeriodicMatch { map, src_facets, shift })
}

pub fn vertex_keys_mixed(mesh: &MixedDiscretisation, matches: &[MixedPeriodicMatch]) -> Vec<u32> {
    let n = mesh.n_points();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(p: &mut [u32], mut i: u32) -> u32 {
        while p[i as usize] != i {
            p[i as usize] = p[p[i as usize] as usize];
            i = p[i as usize];
        }
        i
    }
    for m in matches {
        for (&s, &t) in &m.map {
            let (rs, rt) = (find(&mut parent, s), find(&mut parent, t));
            if rs != rt {
                parent[rs.max(rt) as usize] = rs.min(rt);
            }
        }
    }
    (0..n as u32).map(|i| find(&mut parent, i)).collect()
}

/// Periodic scalar nullspace relations `(dest_col, src_col, exp(-i phase))` matching `gradient_p1_full` (`order == 1`)
/// or `gradient_p2_full` (`order == 2`: vertices `0..nv`, edges `nv..nv+ne`, 3D quad faces `nv+ne..nv+ne+nqf`).
pub fn nodal_relations_mixed_3d(
    mesh: &MixedDiscretisation,
    order: usize,
    m: &MixedPeriodicMatch,
    phase: f64,
) -> Result<Vec<(usize, usize, C)>, String> {
    let f = phase_factor(phase);
    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    let mut out = Vec::new();
    for (&s, &t) in &m.map {
        out.push((t as usize, s as usize, f));
    }
    if order == 2 {
        let mut done_edges = BTreeSet::new();
        let mut done_quads = BTreeSet::new();
        for &fc in &m.src_facets {
            let fv = mesh.facet(fc);
            let ft = mesh.facet_type(fc);
            let le: &[(usize, usize)] = match ft {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in le {
                let (asrc, bsrc) = (fv[i], fv[j]);
                let key = (asrc.min(bsrc), asrc.max(bsrc));
                if !done_edges.insert(key) {
                    continue;
                }
                let es = mesh.edge_index(asrc, bsrc).unwrap();
                let et = mesh
                    .edge_index(m.map[&asrc], m.map[&bsrc])
                    .ok_or_else(|| "periodic source edge has no destination counterpart".to_string())?;
                out.push((nv + et, nv + es, f));
            }
            if ft == FacetType::Quad4 {
                let qfs = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
                if done_quads.insert(qfs) {
                    let qft = mesh
                        .quad_face_index([m.map[&fv[0]], m.map[&fv[1]], m.map[&fv[2]], m.map[&fv[3]]])
                        .ok_or_else(|| "periodic source Quad4 face has no destination counterpart".to_string())?;
                    out.push((nv + ne + qft, nv + ne + qfs, f));
                }
            }
        }
    }
    Ok(out)
}

/// Infers `(lateral, normal, propagation)` polarization axes on a 3D `MixedDiscretisation`.
pub fn axes_of_3d_periodic_mixed(periodic: &[Periodic], mesh: &MixedDiscretisation) -> ([usize; 3], bool) {
    for p in periodic {
        if let Ok(mt) = match_vertices_mixed(mesh, p) {
            let ax = (0..3).max_by(|&a, &b| mt.shift[a].abs().total_cmp(&mt.shift[b].abs())).unwrap_or(1);
            return match ax {
                0 => ([1, 2, 0], false),
                2 => ([0, 1, 2], false),
                _ => ([0, 2, 1], true),
            };
        }
    }
    ([0, 2, 1], true)
}

// ------------------------------------------------------------------ 3D Mixed-cell facets (IBC, JJ, ports, waveports)

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MixedFaceKey {
    Tri([u32; 3]),
    Quad([u32; 4]),
}

impl MixedFaceKey {
    pub fn from_slice(fv: &[u32]) -> Self {
        if fv.len() == 3 {
            MixedFaceKey::Tri(sorted3(fv[0], fv[1], fv[2]))
        } else {
            let mut q = [fv[0], fv[1], fv[2], fv[3]];
            q.sort_unstable();
            MixedFaceKey::Quad(q)
        }
    }

    pub fn vertices(&self) -> Vec<u32> {
        match *self {
            MixedFaceKey::Tri(t) => t.to_vec(),
            MixedFaceKey::Quad(q) => q.to_vec(),
        }
    }
}

/// Maps every canonical 3D cell face (`Tri3` or `Quad4`) to the smallest-index adjacent cell `ci`.
pub fn build_face_to_cell_3d_mixed(mesh: &MixedDiscretisation) -> HashMap<MixedFaceKey, usize> {
    let mut map = HashMap::new();
    for ci in 0..mesh.n_cells() {
        let cv = mesh.cell(ci);
        match mesh.cell_type(ci) {
            CellType::Tet => {
                for lf in &TET_TRI_FACES {
                    map.entry(MixedFaceKey::from_slice(&[cv[lf[0]], cv[lf[1]], cv[lf[2]]])).or_insert(ci);
                }
            }
            CellType::Prism => {
                for lf in &PRISM_TRI_FACES {
                    map.entry(MixedFaceKey::from_slice(&[cv[lf[0]], cv[lf[1]], cv[lf[2]]])).or_insert(ci);
                }
                for lq in &PRISM_QUAD_FACES {
                    map.entry(MixedFaceKey::from_slice(&[cv[lq[0]], cv[lq[1]], cv[lq[2]], cv[lq[3]]])).or_insert(ci);
                }
            }
            CellType::Hex => {
                for lq in &HEX_QUAD_FACES {
                    map.entry(MixedFaceKey::from_slice(&[cv[lq[0]], cv[lq[1]], cv[lq[2]], cv[lq[3]]])).or_insert(ci);
                }
            }
            _ => {}
        }
    }
    map
}

/// Evaluates tangential Nedelec basis vectors `psi_{s,||} = w_s - (w_s . n) n` at quadrature points on 3D facet `f`
/// (`Tri3` or `Quad4`, straight or curved). Returns `(ci, area, n_out_avg, bnd_dofs, quad_points_x, quad_weights, psi_q)`.
pub fn eval_facet_tangential_basis_3d_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    face_to_cell: &HashMap<MixedFaceKey, usize>,
    f: usize,
) -> Result<(usize, f64, [f64; 3], Vec<usize>, Vec<[f64; 3]>, Vec<f64>, Vec<Vec<[f64; 3]>>), String> {
    let fv = mesh.facet(f);
    let ft = mesh.facet_type(f);
    let key = MixedFaceKey::from_slice(fv);
    let &ci = face_to_cell
        .get(&key)
        .ok_or_else(|| format!("3D facet {f} ({fv:?}) is not a face of any cell in MixedDiscretisation"))?;
    let cv = mesh.cell(ci);
    let ct = mesh.cell_type(ci);
    let corners_ref = cell_corner_ref_coords(ct);
    let local_idx: Vec<usize> =
        fv.iter().map(|&v| cv.iter().position(|&x| x == v).expect("facet vertex belongs to adjacent cell")).collect();

    // Cell centroid to orient outward normal `n_out` from cell `ci`.
    let mut cell_cent = [0.0; 3];
    for &v in cv {
        let p = mesh.point3(v as usize);
        for k in 0..3 {
            cell_cent[k] += p[k] / cv.len() as f64;
        }
    }

    let iso = IsoparametricCell::from_mesh(mesh, ci);
    let deg = if iso.curved { CURVED_QUAD_DEGREE } else { 4 };
    let frule = facet_rule_mixed(ft, deg);

    let bnd_dofs = ned.facet_dofs(mesh, f);
    let cdofs = ned.cell_dofs(ci);
    let local_slots: Vec<usize> = bnd_dofs
        .iter()
        .map(|&gd| cdofs.iter().position(|&cd| cd == gd).expect("facet DOF belongs to adjacent cell"))
        .collect();

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut cu = Vec::new();
    let mut area = 0.0;
    let mut n_acc = [0.0; 3];
    let mut pts_x = Vec::with_capacity(frule.w.len());
    let mut weights = Vec::with_capacity(frule.w.len());
    let mut psi_q = Vec::with_capacity(frule.w.len());

    for (xi_f, &wf) in frule.pts.iter().zip(&frule.w) {
        let (fn1, _) = facet_lagrange_shape(ft, LagrangeSpace::Order1, xi_f);
        let mut xi_cell = [0.0; 3];
        for (k, &li) in local_idx.iter().enumerate() {
            let cr = corners_ref[li];
            for a in 0..3 {
                xi_cell[a] += fn1[k] * cr[a];
            }
        }
        let (xq, dsq, mut n_q) = eval_facet_jacobian(mesh, f, xi_f);
        if dot3(n_q, sub3(xq, cell_cent)) < 0.0 {
            n_q = [-n_q[0], -n_q[1], -n_q[2]];
        }
        eval_cell_basis_3d(mesh, ned, ci, &iso, &xi_cell, &mut scratch, &mut w, &mut cu);
        let mut row = Vec::with_capacity(bnd_dofs.len());
        for &s in &local_slots {
            let ws = w[s];
            let wn = dot3(ws, n_q);
            row.push([ws[0] - wn * n_q[0], ws[1] - wn * n_q[1], ws[2] - wn * n_q[2]]);
        }
        let wq = wf * dsq;
        area += wq;
        for a in 0..3 {
            n_acc[a] += wq * n_q[a];
        }
        pts_x.push(xq);
        weights.push(wq);
        psi_q.push(row);
    }

    if area <= 0.0 {
        return Err(format!("degenerate zero-area 3D boundary facet {f}"));
    }
    let nn = norm3(n_acc).max(1e-30);
    let n_out = [n_acc[0] / nn, n_acc[1] / nn, n_acc[2] / nn];
    Ok((ci, area, n_out, bnd_dofs, pts_x, weights, psi_q))
}

fn facet_is_curved_mixed(mesh: &MixedDiscretisation, f: usize) -> bool {
    if !mesh.has_curved() {
        return false;
    }
    let fv = mesh.facet(f);
    let le: &[(usize, usize)] = match mesh.facet_type(f) {
        FacetType::Edge2 => &[(0, 1)],
        FacetType::Tri3 => &TRI_EDGES,
        FacetType::Quad4 => &QUAD_EDGES,
    };
    if le.iter().any(|&(i, j)| mesh.edge_index(fv[i], fv[j]).is_some_and(|e| mesh.edge_node(e).is_some())) {
        return true;
    }
    mesh.facet_type(f) == FacetType::Quad4
        && mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).is_some_and(|qf| mesh.quad_face_node(qf).is_some())
}

fn facet_tri_decomposition(fv: &[u32]) -> Vec<[u32; 3]> {
    if fv.len() == 3 {
        vec![sorted3(fv[0], fv[1], fv[2])]
    } else {
        let m = (0..4).min_by_key(|&i| fv[i]).unwrap();
        let c = |k: usize| fv[(m + k) % 4];
        vec![sorted3(c(0), c(1), c(2)), sorted3(c(0), c(2), c(3))]
    }
}

/// Validates that every `lumped_jj` and `lumped_port` sheet on a 3D `MixedDiscretisation` contacts two conductors.
pub fn validate_port_contacts_3d_mixed(mesh: &MixedDiscretisation, problem: &Problem) -> Result<(), String> {
    if problem.lumped_jj.is_empty() && problem.lumped_ports.is_empty() {
        return Ok(());
    }
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let wp_labels: Vec<i32> = problem.waveports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let (sym_pec, sym_pmc) = symmetry_facet_flags_mixed(mesh, problem.symmetry)?;

    let mut cond_verts = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let is_cond = problem.pec.contains(&l)
            || ibc_labels.contains(&l)
            || sym_pec[f]
            || (mesh.exterior_facets[f]
                && !periodic_labels.contains(&l)
                && !problem.pmc.contains(&l)
                && !sym_pmc[f]
                && !jj_labels.contains(&l)
                && !lp_labels.contains(&l)
                && !wp_labels.contains(&l));
        if is_cond {
            for &v in mesh.facet(f) {
                cond_verts.insert(v);
            }
        }
    }

    let check_one = |kind: &str, name: &str, labels: &[i32], dir: [f64; 3]| -> Result<(), String> {
        let mut port_verts = BTreeSet::new();
        for f in 0..mesh.n_facets() {
            if labels.contains(&mesh.facet_labels[f]) {
                for &v in mesh.facet(f) {
                    port_verts.insert(v);
                }
            }
        }
        if port_verts.is_empty() {
            return Err(format!("{kind} '{name}' selects no facets on the 3D mesh"));
        }
        let (mut u_min, mut u_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &port_verts {
            let u = dot3(mesh.point3(v as usize), dir);
            u_min = u_min.min(u);
            u_max = u_max.max(u);
        }
        let span = (u_max - u_min).max(1e-30);
        let touches_lo = port_verts
            .iter()
            .any(|v| cond_verts.contains(v) && dot3(mesh.point3(*v as usize), dir) <= u_min + 0.25 * span);
        let touches_hi = port_verts
            .iter()
            .any(|v| cond_verts.contains(v) && dot3(mesh.point3(*v as usize), dir) >= u_max - 0.25 * span);
        if !(touches_lo && touches_hi) {
            return Err(format!(
                "{kind} '{name}' must contact conductor boundaries on both ends along direction {:?}",
                dir
            ));
        }
        Ok(())
    };

    for jj in &problem.lumped_jj {
        check_one("lumped_jj", &jj.name, &jj.labels, jj.direction)?;
    }
    for lp in &problem.lumped_ports {
        check_one("lumped_port", &lp.name, &lp.labels, lp.direction)?;
    }
    Ok(())
}

/// Assembles `pec`, `surface_impedance`, and `lumped_jj` boundary conditions on a 3D `MixedDiscretisation`,
/// returning `(IbcJjAssembly3D, pec_facets, ibc_facets, jj_facets_by_port)` (direct IBC mode for `driven`/`bloch`).
pub fn assemble_ibc_and_jj_3d_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f_ref: f64,
) -> Result<(IbcJjAssembly3D, Vec<usize>, Vec<usize>, Vec<Vec<usize>>), String> {
    assemble_ibc_and_jj_3d_mixed_with_mode(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f_ref, false)
}

/// Like [`assemble_ibc_and_jj_3d_mixed`], supporting `ibc_mode` (`"direct" | "perturbative" | "auto"`) for `cavity`.
pub fn assemble_ibc_and_jj_3d_mixed_cavity(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f_ref: f64,
) -> Result<(IbcJjAssembly3D, Vec<usize>, Vec<usize>, Vec<Vec<usize>>), String> {
    assemble_ibc_and_jj_3d_mixed_with_mode(mesh, ned, pec_labels, ibc_bcs, jj_bcs, f_ref, true)
}

fn assemble_ibc_and_jj_3d_mixed_with_mode(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    jj_bcs: &[LumpedJjBc],
    f_ref: f64,
    allow_perturbative: bool,
) -> Result<(IbcJjAssembly3D, Vec<usize>, Vec<usize>, Vec<Vec<usize>>), String> {
    let n = ned.n_dofs;
    let face_to_cell = build_face_to_cell_3d_mixed(mesh);
    let omega_ref = (2.0 * PI * f_ref).max(1e-30);
    let k0_ref = omega_ref / C0;

    let mut pec_seen = BTreeSet::new();
    let mut pec_facets = Vec::new();
    let mut pec_tris = Vec::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            let fv = mesh.facet(f);
            if pec_seen.insert(MixedFaceKey::from_slice(fv)) {
                pec_facets.push(f);
                pec_tris.extend(facet_tri_decomposition(fv));
            }
        }
    }

    let mut k_ibc_coo = Coo::new(n, n);
    let mut c_ibc_coo = Coo::new(n, n);
    let mut ibc_seen = BTreeSet::new();
    let mut ibc_facets = Vec::new();
    let mut ibc_tris = Vec::new();
    let mut pert_ibc_tris = Vec::new();
    let mut direct_ibc_labels = Vec::new();
    let mut warnings = Vec::new();
    let mut cell_sheet = vec![None; mesh.n_cells()];

    for &f in &pec_facets {
        let (ci, _, n_out, _, _, _, _) = eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
        cell_sheet[ci] = Some(SheetCellKind::Pec { n_hat: n_out });
    }

    for bc in ibc_bcs {
        let omega_bc = if bc.f_ref > 0.0 { 2.0 * PI * bc.f_ref } else { omega_ref };
        let (inv_l_eff, g_eff) = ibc_parallel_params(bc.r_s, bc.l_k_sq, omega_bc);
        let k_scale = MU0 * inv_l_eff;
        let c_scale = MU0 * C0 * g_eff;
        let mut bc_facets = Vec::new();
        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = MixedFaceKey::from_slice(fv);
            if pec_seen.contains(&key) || !ibc_seen.insert(key) {
                continue;
            }
            bc_facets.push(f);
        }
        if bc_facets.is_empty() {
            continue;
        }

        let all_tri = bc_facets.iter().all(|&f| mesh.facet_type(f) == FacetType::Tri3);
        let mut edge_len_sum = 0.0_f64;
        let mut area_tot = 0.0_f64;
        let mut edge_counts: BTreeMap<(u32, u32), usize> = BTreeMap::new();
        let mut n_tris_bc = 0usize;
        for &f in &bc_facets {
            for tri in facet_tri_decomposition(mesh.facet(f)) {
                let (p0, p1, p2) =
                    (mesh.point3(tri[0] as usize), mesh.point3(tri[1] as usize), mesh.point3(tri[2] as usize));
                let l01 = norm3(sub3(p1, p0));
                let l12 = norm3(sub3(p2, p1));
                let l20 = norm3(sub3(p0, p2));
                edge_len_sum += (l01 + l12 + l20) / 3.0;
                area_tot += 0.5 * norm3(cross3(sub3(p1, p0), sub3(p2, p0)));
                n_tris_bc += 1;
            }
            let fv = mesh.facet(f);
            let le: &[(usize, usize)] = match mesh.facet_type(f) {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in le {
                let a = fv[i];
                let b = fv[j];
                *edge_counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let h_sheet_mean = edge_len_sum / (n_tris_bc.max(1) as f64);
        let perim_tot: f64 = edge_counts
            .iter()
            .filter(|(_, &cnt)| cnt == 1)
            .map(|(&(a, b), _)| norm3(sub3(mesh.point3(b as usize), mesh.point3(a as usize))))
            .sum();
        let est_tris_across_width = if perim_tot > 0.0 && h_sheet_mean > 0.0 {
            (2.0 * area_tot / perim_tot) / h_sheet_mean
        } else {
            f64::INFINITY
        };
        let lambda_k = bc.l_k_sq / MU0;
        let z_s_mag = bc.r_s.hypot(omega_ref * bc.l_k_sq);
        let use_perturbative = allow_perturbative
            && match bc.ibc_mode.as_str() {
                "perturbative" => {
                    if !all_tri {
                        return Err(
                            "surface_impedance 'perturbative' mode on mixed meshes requires triangular physical facets"
                                .into(),
                        );
                    }
                    true
                }
                "auto" => all_tri && ned.order == 1 && h_sheet_mean > 10.0 * lambda_k && z_s_mag < 0.05 * ETA0,
                _ => false,
            };

        if use_perturbative {
            pec_facets.extend(bc_facets.iter().copied());
            for f in bc_facets {
                let fv = mesh.facet(f);
                let tri = sorted3(fv[0], fv[1], fv[2]);
                pec_tris.push(tri);
                pert_ibc_tris.push(PerturbativeIbcTri { tri, l_k_sq: bc.l_k_sq, r_s: bc.r_s, f_ref });
                let (ci, _, n_out, _, _, _, _) = eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
                cell_sheet[ci] = Some(SheetCellKind::Pec { n_hat: n_out });
            }
            continue;
        }

        direct_ibc_labels.extend(bc.labels.iter().copied());
        if allow_perturbative
            && ned.order == 1
            && ((bc.l_k_sq > 0.0 && h_sheet_mean > 20.0 * lambda_k) || est_tris_across_width < 3.0)
        {
            warnings.push(format!(
                "warning: surface_impedance in 'direct' mode at order 1 has coarse strip resolution (h_tri = {h_sheet_mean:.3e} m, L_k_sq/mu_0 = {lambda_k:.3e} m, est_tris_across_width = {est_tris_across_width:.2}); consider ibc_mode = 'perturbative' or 'auto' or order = 2"
            ));
        }

        for f in bc_facets {
            ibc_facets.push(f);
            ibc_tris.extend(facet_tri_decomposition(mesh.facet(f)));
            let (ci, _, n_out, bnd_dofs, _, weights, psi_q) =
                eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
            cell_sheet[ci] = Some(SheetCellKind::Ibc { inv_l_eff, g_eff, n_hat: n_out });
            let nb_f = bnd_dofs.len();
            let mut m_loc = vec![0.0; nb_f * nb_f];
            for (q, &wq) in weights.iter().enumerate() {
                for a in 0..nb_f {
                    for b in 0..nb_f {
                        m_loc[a * nb_f + b] += wq * dot3(psi_q[q][a], psi_q[q][b]);
                    }
                }
            }
            if k_scale > 0.0 {
                let kb: Vec<f64> = m_loc.iter().map(|&v| v * k_scale).collect();
                k_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &kb);
            }
            if c_scale > 0.0 {
                let cb: Vec<f64> = m_loc.iter().map(|&v| v * c_scale).collect();
                c_ibc_coo.push_block(&bnd_dofs, &bnd_dofs, &cb);
            }
        }
    }

    let k_ibc_full = k_ibc_coo.to_csr();
    let c_ibc_full = c_ibc_coo.to_csr();
    let kc_ibc_full = if c_ibc_full.nnz() > 0 {
        k_ibc_full.to_complex().add_scaled(&c_ibc_full.to_complex(), c(0.0, k0_ref))
    } else {
        k_ibc_full.to_complex()
    };

    let mut k_jj_full = RCsr::zeros(n, n);
    let mut m_jj_full = RCsr::zeros(n, n);
    let mut jj_tris = Vec::new();
    let mut jj_ports = Vec::with_capacity(jj_bcs.len());
    let mut jj_facets_by_port = Vec::with_capacity(jj_bcs.len());

    for (port_index, bc) in jj_bcs.iter().enumerate() {
        let mut seen = BTreeSet::new();
        let mut p_facets = Vec::new();
        let mut p_tris = Vec::new();
        for f in 0..mesh.n_facets() {
            if bc.labels.contains(&mesh.facet_labels[f]) {
                let fv = mesh.facet(f);
                let key = MixedFaceKey::from_slice(fv);
                if !pec_seen.contains(&key) && seen.insert(key) {
                    p_facets.push(f);
                    p_tris.extend(facet_tri_decomposition(fv));
                }
            }
        }
        if p_facets.is_empty() {
            return Err(format!("lumped_jj '{}' selects no facets", bc.name));
        }
        let (n_hat, d_hat, w_hat, u_min, u_max, length_m, width_m) =
            sheet_geometry_mixed(mesh, &face_to_cell, ned, &p_facets, bc.direction, "lumped_jj", &bc.name)?;
        let inv_l_tot = jj_inv_l_tot(bc.l_j, bc.e_l_hz);

        let mut b_full = vec![0.0; n];
        let mut dist_coo = Coo::new(n, n);
        let scale_b = 1.0 / width_m.max(1e-30);
        let scale_dist = length_m / width_m.max(1e-30);

        for &f in &p_facets {
            let (ci, _, _, bnd_dofs, _, weights, psi_q) =
                eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
            cell_sheet[ci] = Some(SheetCellKind::Jj { port_index, d_hat, width_m });
            let nb_f = bnd_dofs.len();
            for (q, &wq) in weights.iter().enumerate() {
                for a in 0..nb_f {
                    b_full[bnd_dofs[a]] += wq * dot3(psi_q[q][a], d_hat) * scale_b;
                }
                if bc.formulation == "distributed_sheet" {
                    for a in 0..nb_f {
                        let pa = dot3(psi_q[q][a], d_hat);
                        for b in 0..nb_f {
                            let pb = dot3(psi_q[q][b], d_hat);
                            dist_coo.push(bnd_dofs[a], bnd_dofs[b], wq * scale_dist * pa * pb);
                        }
                    }
                }
            }
        }

        let rank1_or_dist =
            if bc.formulation == "distributed_sheet" { dist_coo.to_csr() } else { outer_product_sparse(&b_full) };
        let k_port_full = rank1_or_dist.scaled(MU0 * inv_l_tot);
        let m_port_full = rank1_or_dist.scaled(MU0 * C0 * C0 * bc.c_j);
        if k_port_full.nnz() > 0 {
            k_jj_full = k_jj_full.add_scaled(&k_port_full, 1.0);
        }
        if m_port_full.nnz() > 0 {
            m_jj_full = m_jj_full.add_scaled(&m_port_full, 1.0);
        }

        jj_tris.extend(p_tris.iter().copied());
        jj_facets_by_port.push(p_facets);
        jj_ports.push(AssembledJjPort {
            bc: bc.clone(),
            d_hat,
            w_hat,
            n_hat,
            length_m,
            width_m,
            u_min,
            u_max,
            inv_l_tot,
            triangles: p_tris,
            b_full,
            b_red: Vec::new(),
            k_jj_red: k_port_full,
            m_jj_red: m_port_full,
        });
    }

    Ok((
        IbcJjAssembly3D {
            pec_tris,
            ibc_tris,
            pert_ibc_tris,
            direct_ibc_labels,
            jj_tris,
            k_ibc_full,
            c_ibc_full,
            kc_ibc_full,
            k_jj_full,
            m_jj_full,
            jj_ports,
            cell_sheet,
            warnings,
        },
        pec_facets,
        ibc_facets,
        jj_facets_by_port,
    ))
}

/// Appends penetrable 2D tensor surface conductivity (`sheet_2d`) jump integrals to an [`IbcJjAssembly3D`]
/// on a 3D [`MixedDiscretisation`] (`Tet + Prism + Hex`). Returns `Ok(true)` if any sheet has $\sigma_{uv} \neq \sigma_{vu}$.
pub fn assemble_sheet2d_3d_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    pec_labels: &[i32],
    sheet_bcs: &[Sheet2dBc],
    imported_potential: Option<&[f64]>,
    f_ref: f64,
    ij: &mut IbcJjAssembly3D,
) -> Result<bool, String> {
    if sheet_bcs.is_empty() {
        return Ok(false);
    }
    let n = ned.n_dofs;
    let face_to_cell = build_face_to_cell_3d_mixed(mesh);
    let omega_ref = (2.0 * PI * f_ref).max(1e-30);
    let k0_ref = omega_ref / C0;

    let mut seen = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        if pec_labels.contains(&mesh.facet_labels[f]) {
            seen.insert(MixedFaceKey::from_slice(mesh.facet(f)));
        }
    }

    let mut k_2d_coo = Coo::new(n, n);
    let mut c_2d_coo = Coo::new(n, n);
    let mut has_non_symmetric = false;

    for bc in sheet_bcs {
        if bc.is_non_symmetric() {
            has_non_symmetric = true;
        }
        let f_eval = if bc.f_ref > 0.0 { bc.f_ref } else { f_ref };
        ij.direct_ibc_labels.extend(bc.labels.iter().copied());

        for f in 0..mesh.n_facets() {
            if !bc.labels.contains(&mesh.facet_labels[f]) {
                continue;
            }
            let fv = mesh.facet(f);
            let ft = mesh.facet_type(f);
            let key = MixedFaceKey::from_slice(fv);
            if !seen.insert(key) {
                continue;
            }
            let (ci, _, n_out, bnd_dofs, _, weights, psi_q) =
                eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
            let (u_hat, v_hat, _) = sheet2d_tangent_frame_3d(n_out);
            let iso = IsoparametricCell::from_mesh(mesh, ci);
            let deg = if iso.curved { CURVED_QUAD_DEGREE } else { 4 };
            let frule = facet_rule_mixed(ft, deg);

            let nb_f = bnd_dofs.len();
            let mut k_loc = vec![0.0; nb_f * nb_f];
            let mut c_loc = vec![0.0; nb_f * nb_f];
            let mut sig_rep = [[c(0.0, 0.0); 2]; 2];

            for (q, (&wq, xi_f)) in weights.iter().zip(&frule.pts).enumerate() {
                let local_v = imported_potential.and_then(|pot| {
                    if pot.len() >= mesh.n_points() {
                        let (fn1, _) = facet_lagrange_shape(ft, LagrangeSpace::Order1, xi_f);
                        Some(fn1.iter().zip(fv).map(|(&sh, &vid)| sh * pot[vid as usize]).sum::<f64>())
                    } else {
                        None
                    }
                });
                let sig = bc.eval_tensor_at(f_eval, None, local_v)?;
                if q == 0 {
                    sig_rep = sig;
                }
                if (sig[0][1] - sig[1][0]).norm() > 1e-18 * (sig[0][0].norm() + sig[1][1].norm()).max(1e-30) {
                    has_non_symmetric = true;
                }
                let k_tensor = [
                    [-omega_ref * MU0 * sig[0][0].im, -omega_ref * MU0 * sig[0][1].im],
                    [-omega_ref * MU0 * sig[1][0].im, -omega_ref * MU0 * sig[1][1].im],
                ];
                let c_tensor = [
                    [MU0 * C0 * sig[0][0].re, MU0 * C0 * sig[0][1].re],
                    [MU0 * C0 * sig[1][0].re, MU0 * C0 * sig[1][1].re],
                ];
                for a in 0..nb_f {
                    let wa = psi_q[q][a];
                    let wa_u = dot3(wa, u_hat);
                    let wa_v = dot3(wa, v_hat);
                    for b in 0..nb_f {
                        let wb = psi_q[q][b];
                        let wb_u = dot3(wb, u_hat);
                        let wb_v = dot3(wb, v_hat);
                        let k_ab = wa_u * (k_tensor[0][0] * wb_u + k_tensor[0][1] * wb_v)
                            + wa_v * (k_tensor[1][0] * wb_u + k_tensor[1][1] * wb_v);
                        let c_ab = wa_u * (c_tensor[0][0] * wb_u + c_tensor[0][1] * wb_v)
                            + wa_v * (c_tensor[1][0] * wb_u + c_tensor[1][1] * wb_v);
                        k_loc[a * nb_f + b] += wq * k_ab;
                        c_loc[a * nb_f + b] += wq * c_ab;
                    }
                }
            }

            if k_loc.iter().any(|&v| v != 0.0) {
                k_2d_coo.push_block(&bnd_dofs, &bnd_dofs, &k_loc);
            }
            if c_loc.iter().any(|&v| v != 0.0) {
                c_2d_coo.push_block(&bnd_dofs, &bnd_dofs, &c_loc);
            }
            ij.cell_sheet[ci] = Some(SheetCellKind::Sheet2D { sigma_2d: sig_rep, u_hat, v_hat, n_hat: n_out });
        }
    }

    let k_2d = k_2d_coo.to_csr();
    let c_2d = c_2d_coo.to_csr();
    if k_2d.nnz() > 0 {
        ij.k_ibc_full = ij.k_ibc_full.add_scaled(&k_2d, 1.0);
        ij.kc_ibc_full = ij.kc_ibc_full.add_scaled(&k_2d.to_complex(), c(1.0, 0.0));
    }
    if c_2d.nnz() > 0 {
        ij.c_ibc_full = ij.c_ibc_full.add_scaled(&c_2d, 1.0);
        ij.kc_ibc_full = ij.kc_ibc_full.add_scaled(&c_2d.to_complex(), c(0.0, k0_ref));
    }
    Ok(has_non_symmetric)
}

/// Assembles `lumped_port` conditions on a 3D `MixedDiscretisation`.
pub fn assemble_lumped_ports_3d_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    pec_labels: &[i32],
    lumped_ports: &[LumpedPortBc],
) -> Result<(Vec<AssembledLumpedPort>, Vec<Vec<usize>>), String> {
    let n = ned.n_dofs;
    let face_to_cell = build_face_to_cell_3d_mixed(mesh);
    let pec_keys: BTreeSet<MixedFaceKey> = (0..mesh.n_facets())
        .filter(|&f| pec_labels.contains(&mesh.facet_labels[f]))
        .map(|f| MixedFaceKey::from_slice(mesh.facet(f)))
        .collect();

    let mut out = Vec::with_capacity(lumped_ports.len());
    let mut facets_by_port = Vec::with_capacity(lumped_ports.len());

    for bc in lumped_ports {
        let mut seen = BTreeSet::new();
        let mut p_facets = Vec::new();
        let mut p_tris = Vec::new();
        for f in 0..mesh.n_facets() {
            if bc.labels.contains(&mesh.facet_labels[f]) {
                let fv = mesh.facet(f);
                let key = MixedFaceKey::from_slice(fv);
                if !pec_keys.contains(&key) && seen.insert(key) {
                    p_facets.push(f);
                    p_tris.extend(facet_tri_decomposition(fv));
                }
            }
        }
        if p_facets.is_empty() {
            return Err(format!("lumped_port '{}' selects no facets", bc.name));
        }

        let (n_hat, d_hat, w_hat, u_min, u_max, length_m, width_m) =
            sheet_geometry_mixed(mesh, &face_to_cell, ned, &p_facets, bc.direction, "lumped_port", &bc.name)?;

        let mut b_full = vec![0.0; n];
        let cpw_slots = if bc.mode == "cpw_even" {
            let comps = partition_facet_components(mesh, &p_facets);
            if comps.len() != 2 {
                return Err(format!(
                    "lumped_port '{}' with mode = 'cpw_even' requires 2 disjoint slot components (found {})",
                    bc.name,
                    comps.len()
                ));
            }
            let cent_u = |fcs: &[usize]| -> f64 {
                let mut vs = BTreeSet::new();
                for &f in fcs {
                    for &v in mesh.facet(f) {
                        vs.insert(v);
                    }
                }
                vs.iter().map(|&v| dot3(mesh.point3(v as usize), d_hat)).sum::<f64>() / vs.len() as f64
            };
            let (left_fcs, right_fcs) = if cent_u(&comps[0]) <= cent_u(&comps[1]) {
                (comps[0].clone(), comps[1].clone())
            } else {
                (comps[1].clone(), comps[0].clone())
            };
            for (sign, slot_fcs) in [(0.5_f64, &left_fcs), (-0.5_f64, &right_fcs)] {
                let (_, _, _, _, _, _, slot_w) =
                    sheet_geometry_mixed(mesh, &face_to_cell, ned, slot_fcs, bc.direction, "lumped_port", &bc.name)?;
                let scale = sign / slot_w.max(1e-30);
                for &f in slot_fcs {
                    let (_, _, _, bnd_dofs, _, weights, psi_q) =
                        eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
                    for (q, &wq) in weights.iter().enumerate() {
                        for (a, &gd) in bnd_dofs.iter().enumerate() {
                            b_full[gd] += wq * dot3(psi_q[q][a], d_hat) * scale;
                        }
                    }
                }
            }
            let tris_l: Vec<[u32; 3]> = left_fcs.iter().flat_map(|&f| facet_tri_decomposition(mesh.facet(f))).collect();
            let tris_r: Vec<[u32; 3]> =
                right_fcs.iter().flat_map(|&f| facet_tri_decomposition(mesh.facet(f))).collect();
            Some((tris_l, tris_r))
        } else {
            let scale = 1.0 / width_m.max(1e-30);
            for &f in &p_facets {
                let (_, _, _, bnd_dofs, _, weights, psi_q) =
                    eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
                for (q, &wq) in weights.iter().enumerate() {
                    for (a, &gd) in bnd_dofs.iter().enumerate() {
                        b_full[gd] += wq * dot3(psi_q[q][a], d_hat) * scale;
                    }
                }
            }
            None
        };

        let c_port_full = outer_product_sparse(&b_full).scaled(MU0 * C0 / bc.z0.max(1e-30));
        facets_by_port.push(p_facets);
        out.push(AssembledLumpedPort {
            bc: bc.clone(),
            d_hat,
            w_hat,
            n_hat,
            length_m,
            width_m,
            u_min,
            u_max,
            triangles: p_tris,
            cpw_slots,
            b_full,
            b_red: Vec::new(),
            c_port_full,
        });
    }
    Ok((out, facets_by_port))
}

/// Assembles 2D boundary-mode `waveport` conditions on a 3D `MixedDiscretisation` (supporting both `Tri3` and `Quad4` port facets).
pub fn assemble_waveports_3d_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    problem: &Problem,
    f_min: f64,
    f_max: f64,
) -> Result<(Vec<AssembledWavePort>, Vec<Vec<usize>>), String> {
    if problem.waveports.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }
    let n = ned.n_dofs;
    let order = ned.order;
    let face_to_cell = build_face_to_cell_3d_mixed(mesh);
    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let (sym_pec, sym_pmc) = symmetry_facet_flags_mixed(mesh, problem.symmetry)?;
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let wp_all_labels: Vec<i32> = problem.waveports.iter().flat_map(|b| b.labels.iter().copied()).collect();

    // Collect conductor edges in 3D so the 2D waveport boundary mesh inherits exact PEC/PMC/IBC 1D edges.
    let mut pec_edges_3d = BTreeSet::new();
    let mut pmc_edges_3d = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        let le: &[(usize, usize)] = match mesh.facet_type(f) {
            FacetType::Edge2 => &[(0, 1)],
            FacetType::Tri3 => &TRI_EDGES,
            FacetType::Quad4 => &QUAD_EDGES,
        };
        let is_pec = problem.pec.contains(&l)
            || sym_pec[f]
            || (mesh.exterior_facets[f]
                && !periodic_labels.contains(&l)
                && !problem.pmc.contains(&l)
                && !sym_pmc[f]
                && !ibc_labels.contains(&l)
                && !jj_labels.contains(&l)
                && !lp_labels.contains(&l)
                && !wp_all_labels.contains(&l));
        let is_pmc = problem.pmc.contains(&l) || sym_pmc[f];
        for &(i, j) in le {
            let key = (fv[i].min(fv[j]), fv[i].max(fv[j]));
            if is_pec {
                pec_edges_3d.insert(key);
            } else if is_pmc {
                pmc_edges_3d.insert(key);
            }
        }
    }

    let mut out = Vec::with_capacity(problem.waveports.len());
    let mut facets_by_port = Vec::with_capacity(problem.waveports.len());

    for bc in &problem.waveports {
        let mut seen = BTreeSet::new();
        let mut p_facets = Vec::new();
        let mut p_tris = Vec::new();
        for f in 0..mesh.n_facets() {
            if bc.labels.contains(&mesh.facet_labels[f]) {
                if !mesh.exterior_facets[f] {
                    return Err(format!("waveport '{}' must sit on an exterior boundary face", bc.name));
                }
                let fv = mesh.facet(f);
                if seen.insert(MixedFaceKey::from_slice(fv)) {
                    p_facets.push(f);
                    p_tris.extend(facet_tri_decomposition(fv));
                }
            }
        }
        if p_facets.is_empty() {
            return Err(format!("waveport '{}' selects no exterior facets", bc.name));
        }
        for &f in &p_facets {
            let fv = mesh.facet(f);
            let ci = face_to_cell[&MixedFaceKey::from_slice(fv)];
            if facet_is_curved_mixed(mesh, f) || mesh.cell_is_curved(ci) {
                return Err(format!(
                    "waveport '{}' does not support curved port facets or adjacent curved cells",
                    bc.name
                ));
            }
        }

        // Compute port plane frame `(u_hat, v_hat, n_out)`
        let (_, _, n_out, _, _, _, _) = eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, p_facets[0])?;
        let ref_ax = if n_out[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u_raw = cross3(ref_ax, n_out);
        let un = norm3(u_raw).max(1e-30);
        let u_hat = [u_raw[0] / un, u_raw[1] / un, u_raw[2] / un];
        let v_hat = cross3(n_out, u_hat);

        // Build 2D SimplexMesh on the port plane from triangulated facets (`Tri3` -> 1 tri, `Quad4` -> 2 tris).
        let mut v3_to_v2: BTreeMap<u32, u32> = BTreeMap::new();
        let mut pts2 = Vec::new();
        let p_origin = mesh.point3(mesh.facet(p_facets[0])[0] as usize);
        for &f in &p_facets {
            for &v in mesh.facet(f) {
                v3_to_v2.entry(v).or_insert_with(|| {
                    let id = (pts2.len() / 2) as u32;
                    let d = sub3(mesh.point3(v as usize), p_origin);
                    pts2.push(dot3(d, u_hat));
                    pts2.push(dot3(d, v_hat));
                    id
                });
            }
        }

        let mut cells2 = Vec::new();
        let mut cell_labels2 = Vec::new();
        let mut eps2 = Vec::new();
        let mut mur2 = Vec::new();
        let mut edge_count2: BTreeMap<(u32, u32), (usize, u32, u32)> = BTreeMap::new();

        for &f in &p_facets {
            let fv = mesh.facet(f);
            let ci = face_to_cell[&MixedFaceKey::from_slice(fv)];
            let dom = mesh.cell_labels[ci];
            let er = &mats.eps_phys_re[ci];
            let euu = dot3(u_hat, [dot3(er[0], u_hat), dot3(er[1], u_hat), dot3(er[2], u_hat)]);
            let evv = dot3(v_hat, [dot3(er[0], v_hat), dot3(er[1], v_hat), dot3(er[2], v_hat)]);
            let euv = dot3(u_hat, [dot3(er[0], v_hat), dot3(er[1], v_hat), dot3(er[2], v_hat)]);
            let enn = dot3(n_out, [dot3(er[0], n_out), dot3(er[1], n_out), dot3(er[2], n_out)]);
            let ec = EpsCell { xx: euu, xy: euv, yy: evv, zz: enn };
            let mr = mats.mur_phys[ci];

            let le: &[(usize, usize)] = match mesh.facet_type(f) {
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
                _ => &[],
            };
            for &(i, j) in le {
                let (a3, b3) = (fv[i], fv[j]);
                let (a2, b2) = (v3_to_v2[&a3], v3_to_v2[&b3]);
                let k2 = (a2.min(b2), a2.max(b2));
                let entry = edge_count2.entry(k2).or_insert((0, a3.min(b3), a3.max(b3)));
                entry.0 += 1;
            }

            for tri3 in facet_tri_decomposition(fv) {
                let (i0, mut i1, mut i2) = (v3_to_v2[&tri3[0]], v3_to_v2[&tri3[1]], v3_to_v2[&tri3[2]]);
                let (p0, p1, p2) = (
                    [pts2[2 * i0 as usize], pts2[2 * i0 as usize + 1]],
                    [pts2[2 * i1 as usize], pts2[2 * i1 as usize + 1]],
                    [pts2[2 * i2 as usize], pts2[2 * i2 as usize + 1]],
                );
                if (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]) < 0.0 {
                    std::mem::swap(&mut i1, &mut i2);
                }
                cells2.extend_from_slice(&[i0, i1, i2]);
                cell_labels2.push(dom);
                eps2.push(ec);
                mur2.push(mr);
            }
        }

        let mut facets2 = Vec::new();
        let mut facet_labels2 = Vec::new();
        for (&(a2, b2), &(cnt, a3, b3)) in &edge_count2 {
            let k3 = (a3, b3);
            if pec_edges_3d.contains(&k3) {
                facets2.extend_from_slice(&[a2, b2]);
                facet_labels2.push(1); // PEC
            } else if pmc_edges_3d.contains(&k3) || cnt == 1 {
                facets2.extend_from_slice(&[a2, b2]);
                facet_labels2.push(if pmc_edges_3d.contains(&k3) { 2 } else { 1 });
            }
        }
        let mesh2 = SimplexMesh::new(2, pts2, cells2, cell_labels2, facets2, facet_labels2);
        let geo2 = geometry(&mesh2);
        let cd2 = cell_data(&mesh2, order);

        // Solve 2D boundary mode at 3 anchor frequencies across [f_min, f_max] to fit quadratic beta^2(k0^2).
        let f_mid = 0.5 * (f_min + f_max);
        let anchor_fs =
            if (f_max - f_min).abs() <= 1e-9 * f_mid.max(1.0) { [f_mid, f_mid, f_mid] } else { [f_min, f_mid, f_max] };
        let n_max =
            eps2.iter().zip(&mur2).map(|(e, &m)| ((e.xx.max(e.yy)) * m).max(1.0).sqrt()).fold(1.0_f64, f64::max);
        let (n_search, k_ask) = if let Some(ng) = bc.n_eff_guess {
            (ng, (bc.mode + 6).max(8))
        } else {
            (1.02 * n_max, (bc.mode + 3).max(4))
        };

        let solve_2d_sorted = |fa: f64| -> Result<Vec<(f64, Vec<f64>, f64)>, String> {
            let asm2 = assemble_modes_with_bcs(&mesh2, fa, &eps2, &mur2, order, &[1], &[], &[2])?;
            let k_actual = k_ask.min(asm2.free.len().saturating_sub(2)).max(1);
            let m2 = solve_mode_pairs(&asm2, n_search, false, k_actual)?;
            if m2.neff.is_empty() {
                return Err(format!("waveport '{}': no 2D guided modes converged at {:.4} GHz", bc.name, fa / 1e9));
            }
            let k0_a = 2.0 * PI * fa / C0;
            let mut cands: Vec<(f64, Vec<f64>, f64)> = m2
                .neff
                .iter()
                .enumerate()
                .map(|(i, &(nr, ni))| {
                    let b2 = if nr > 0.0 { (nr * k0_a).powi(2) } else { -(ni * k0_a).powi(2) };
                    let r = m2.residuals.get(i).copied().unwrap_or(0.0);
                    (b2, m2.et[i].clone(), r)
                })
                .collect();
            if let Some(ng) = bc.n_eff_guess {
                cands.sort_by(|a, b| match (a.0 > 0.0, b.0 > 0.0) {
                    (true, false) => std::cmp::Ordering::Less,
                    (false, true) => std::cmp::Ordering::Greater,
                    (true, true) => {
                        let da = (a.0.sqrt() / k0_a - ng).abs();
                        let db = (b.0.sqrt() / k0_a - ng).abs();
                        da.total_cmp(&db)
                    }
                    (false, false) => b.0.total_cmp(&a.0),
                });
            } else {
                cands.sort_by(|a, b| b.0.total_cmp(&a.0));
            }
            Ok(cands)
        };

        let mid_cands = solve_2d_sorted(f_mid)?;
        if bc.mode > mid_cands.len() {
            return Err(format!(
                "waveport '{}': requested mode {} but only {} mode(s) converged on 2D port cross-section",
                bc.name,
                bc.mode,
                mid_cands.len()
            ));
        }
        let (beta2_mid, mut et_mid, mid_res) = mid_cands[bc.mode - 1].clone();
        if beta2_mid <= 0.0 {
            return Err(format!(
                "waveport '{}' mode {} is evanescent (below cutoff) at midband {:.4} GHz",
                bc.name,
                bc.mode,
                f_mid / 1e9
            ));
        }

        let q2 = quad_rule(2, 4);
        let nloc2 = if order == 1 { 3 } else { 8 };
        let overlap_2d = |et_a: &[f64], et_b: &[f64]| -> f64 {
            let mut wa = vec![[0.0; 2]; nloc2];
            let mut ca = vec![0.0; nloc2];
            let mut acc = 0.0;
            for tc in 0..mesh2.n_cells() {
                let g2 = tri_g(&geo2, tc);
                let inv_mur = 1.0 / mur2[tc];
                for (q, l2) in q2.lam.iter().enumerate() {
                    edge_basis(l2, &g2, &cd2.sign[tc], order, &mut wa, &mut ca);
                    let (mut ea, mut eb) = ([0.0; 2], [0.0; 2]);
                    for s in 0..nloc2 {
                        let gd = cd2.gdofs[tc][s];
                        ea[0] += et_a[gd] * wa[s][0];
                        ea[1] += et_a[gd] * wa[s][1];
                        eb[0] += et_b[gd] * wa[s][0];
                        eb[1] += et_b[gd] * wa[s][1];
                    }
                    acc += q2.w[q] * geo2.vol[tc] * inv_mur * (ea[0] * eb[0] + ea[1] * eb[1]);
                }
            }
            acc
        };

        let norm_sq = overlap_2d(&et_mid, &et_mid);
        let e_scale = 1.0 / norm_sq.max(1e-30).sqrt();
        for v in &mut et_mid {
            *v *= e_scale;
        }

        let mut min_mode_overlap = 1.0_f64;
        let mut eigen_residuals = vec![mid_res];
        let beta2_samples = if (f_max - f_min).abs() <= 1e-9 * f_mid.max(1.0) {
            [beta2_mid, beta2_mid, beta2_mid]
        } else {
            let track_at = |fa: f64| -> Result<(f64, f64, f64), String> {
                let cands = solve_2d_sorted(fa)?;
                let mut best_ov = -1.0_f64;
                let mut best_b2 = beta2_mid;
                let mut best_res = mid_res;
                for (b2, et_c, r) in cands {
                    let ov = overlap_2d(&et_c, &et_mid).abs();
                    if ov > best_ov {
                        best_ov = ov;
                        best_b2 = b2;
                        best_res = r;
                    }
                }
                Ok((best_b2, best_ov.max(0.0), best_res))
            };
            let (b2_min, ov_min, res_min) = track_at(f_min)?;
            let (b2_max, ov_max, res_max) = track_at(f_max)?;
            min_mode_overlap = ov_min.min(ov_max);
            eigen_residuals.push(res_min);
            eigen_residuals.push(res_max);
            [b2_min, beta2_mid, b2_max]
        };
        let beta2_poly = fit_beta2_poly(anchor_fs, beta2_samples);
        let max_eigen_residual = eigen_residuals.iter().copied().fold(0.0_f64, f64::max);
        let mut warnings = Vec::new();
        if min_mode_overlap < 0.80 {
            warnings.push(nyx_physics::DiagnosticWarning {
                code: "waveport_mode_overlap_low".into(),
                severity: nyx_physics::DiagnosticSeverity::Error,
                message: format!(
                    "waveport '{}' mode {} modal overlap dropped to {:.3} (< 0.80) across sweep endpoints",
                    bc.name, bc.mode, min_mode_overlap
                ),
            });
        } else if min_mode_overlap < 0.90 {
            warnings.push(nyx_physics::DiagnosticWarning {
                code: "waveport_mode_overlap_low".into(),
                severity: nyx_physics::DiagnosticSeverity::Warning,
                message: format!(
                    "waveport '{}' mode {} modal overlap dropped to {:.3} (< 0.90) across sweep endpoints",
                    bc.name, bc.mode, min_mode_overlap
                ),
            });
        }

        // Assemble 3D `W_p` and `b_p` directly over the 3D `Tri3` / `Quad4` port facets!
        let mut w_coo = Coo::new(n, n);
        let mut b_full = vec![0.0; n];
        for &f in &p_facets {
            let (ci, _, _, bnd_dofs, pts_x, weights, psi_q) =
                eval_facet_tangential_basis_3d_mixed(mesh, ned, &face_to_cell, f)?;
            let inv_mur = 1.0 / mats.mur_phys[ci];
            let nb_f = bnd_dofs.len();
            let mut w_loc = vec![0.0; nb_f * nb_f];
            for (q, (&xq, &wq)) in pts_x.iter().zip(&weights).enumerate() {
                let d = sub3(xq, p_origin);
                let uv = [dot3(d, u_hat), dot3(d, v_hat)];
                let e2 = eval_2d_port_mode_at(&mesh2, &geo2, &cd2, order, &et_mid, uv);
                let e3 = [
                    e2[0] * u_hat[0] + e2[1] * v_hat[0],
                    e2[0] * u_hat[1] + e2[1] * v_hat[1],
                    e2[0] * u_hat[2] + e2[1] * v_hat[2],
                ];
                for a in 0..nb_f {
                    b_full[bnd_dofs[a]] += wq * inv_mur * dot3(psi_q[q][a], e3);
                    for b in 0..nb_f {
                        w_loc[a * nb_f + b] += wq * inv_mur * dot3(psi_q[q][a], psi_q[q][b]);
                    }
                }
            }
            w_coo.push_block(&bnd_dofs, &bnd_dofs, &w_loc);
        }

        let w_full = w_coo.to_csr();
        facets_by_port.push(p_facets);
        out.push(AssembledWavePort {
            bc: bc.clone(),
            n_out,
            u_hat,
            v_hat,
            triangles: p_tris,
            b_full,
            b_red: Vec::new(),
            w_full,
            w_red: RCsr::zeros(0, 0),
            beta2_poly,
            min_mode_overlap,
            max_eigen_residual,
            eigen_residuals,
            warnings,
        });
    }

    Ok((out, facets_by_port))
}

fn eval_2d_port_mode_at(
    mesh2: &SimplexMesh,
    geo2: &femcore::fem::Geometry,
    cd2: &femcore::edge::CellData,
    order: usize,
    et: &[f64],
    uv: [f64; 2],
) -> [f64; 2] {
    let nloc2 = if order == 1 { 3 } else { 8 };
    let mut w2 = vec![[0.0; 2]; nloc2];
    let mut curl2 = vec![0.0; nloc2];
    let mut best_tc = 0usize;
    let mut best_min_lam = f64::NEG_INFINITY;
    let mut best_lam = [1.0 / 3.0; 3];

    for tc in 0..mesh2.n_cells() {
        let cv = mesh2.cell(tc);
        let p0 = mesh2.point(cv[0] as usize);
        let p1 = mesh2.point(cv[1] as usize);
        let p2 = mesh2.point(cv[2] as usize);
        let det = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);
        if det.abs() <= 1e-30 {
            continue;
        }
        let inv = 1.0 / det;
        let l1 = ((uv[0] - p0[0]) * (p2[1] - p0[1]) - (uv[1] - p0[1]) * (p2[0] - p0[0])) * inv;
        let l2 = ((p1[0] - p0[0]) * (uv[1] - p0[1]) - (p1[1] - p0[1]) * (uv[0] - p0[0])) * inv;
        let l0 = 1.0 - l1 - l2;
        let min_l = l0.min(l1).min(l2);
        if min_l > best_min_lam {
            best_min_lam = min_l;
            best_tc = tc;
            best_lam = [l0, l1, l2];
            if min_l >= -1e-10 {
                break;
            }
        }
    }
    let g2 = tri_g(geo2, best_tc);
    edge_basis(&best_lam, &g2, &cd2.sign[best_tc], order, &mut w2, &mut curl2);
    let mut e2 = [0.0; 2];
    for s in 0..nloc2 {
        let cf = et[cd2.gdofs[best_tc][s]];
        e2[0] += cf * w2[s][0];
        e2[1] += cf * w2[s][1];
    }
    e2
}

fn fit_beta2_poly(fs: [f64; 3], b2: [f64; 3]) -> [f64; 3] {
    let s = fs.map(|f| (2.0 * PI * f / C0).powi(2));
    if (s[2] - s[0]).abs() <= 1e-18 * s[1].max(1.0) {
        return [0.0, b2[1] / s[1].max(1e-30), 0.0];
    }
    let l01 = (b2[1] - b2[0]) / (s[1] - s[0]);
    let l12 = (b2[2] - b2[1]) / (s[2] - s[1]);
    let c2 = (l12 - l01) / (s[2] - s[0]);
    let c1 = l01 - c2 * (s[0] + s[1]);
    let c0 = b2[0] - c1 * s[0] - c2 * s[0] * s[0];
    [c0, c1, c2]
}

fn sheet_geometry_mixed(
    mesh: &MixedDiscretisation,
    face_to_cell: &HashMap<MixedFaceKey, usize>,
    ned: &NedMixed,
    facets: &[usize],
    direction: [f64; 3],
    kind: &str,
    name: &str,
) -> Result<([f64; 3], [f64; 3], [f64; 3], f64, f64, f64, f64), String> {
    let mut n_acc = [0.0; 3];
    let mut n_ref: Option<[f64; 3]> = None;
    let mut verts = BTreeSet::new();
    for &f in facets {
        let (_, area, mut nr, _, _, _, _) = eval_facet_tangential_basis_3d_mixed(mesh, ned, face_to_cell, f)?;
        if let Some(r) = n_ref {
            if dot3(nr, r) < 0.0 {
                nr = [-nr[0], -nr[1], -nr[2]];
            }
        } else {
            n_ref = Some(nr);
        }
        for k in 0..3 {
            n_acc[k] += area * nr[k];
        }
        for &v in mesh.facet(f) {
            verts.insert(v);
        }
    }
    let nn = norm3(n_acc);
    if nn <= 0.0 {
        return Err(format!("{kind} '{name}': degenerate sheet normal"));
    }
    let n_hat = [n_acc[0] / nn, n_acc[1] / nn, n_acc[2] / nn];
    let dn = dot3(direction, n_hat);
    let d_proj = [direction[0] - dn * n_hat[0], direction[1] - dn * n_hat[1], direction[2] - dn * n_hat[2]];
    let dp_norm = norm3(d_proj);
    if dp_norm <= 1e-12 {
        return Err(format!("{kind} '{name}': direction vector is normal to sheet"));
    }
    let d_hat = [d_proj[0] / dp_norm, d_proj[1] / dp_norm, d_proj[2] / dp_norm];
    let w_hat = cross3(n_hat, d_hat);

    let (mut u_min, mut u_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut w_min, mut w_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for &v in &verts {
        let p = mesh.point3(v as usize);
        let u = dot3(p, d_hat);
        let w = dot3(p, w_hat);
        u_min = u_min.min(u);
        u_max = u_max.max(u);
        w_min = w_min.min(w);
        w_max = w_max.max(w);
    }
    let length_m = (u_max - u_min).max(1e-15);
    let width_m = (w_max - w_min).max(1e-15);
    Ok((n_hat, d_hat, w_hat, u_min, u_max, length_m, width_m))
}

fn partition_facet_components(mesh: &MixedDiscretisation, facets: &[usize]) -> Vec<Vec<usize>> {
    let m = facets.len();
    if m == 0 {
        return Vec::new();
    }
    let mut parent: Vec<usize> = (0..m).collect();
    fn find(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let mut vert_owner: HashMap<u32, usize> = HashMap::new();
    for (fi, &f) in facets.iter().enumerate() {
        for &v in mesh.facet(f) {
            if let Some(&prev) = vert_owner.get(&v) {
                let (ra, rb) = (find(&mut parent, prev), find(&mut parent, fi));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            } else {
                vert_owner.insert(v, fi);
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (fi, &f) in facets.iter().enumerate() {
        let r = find(&mut parent, fi);
        groups.entry(r).or_default().push(f);
    }
    groups.into_values().collect()
}

fn outer_product_sparse(b: &[f64]) -> RCsr {
    let n = b.len();
    let nz: Vec<(usize, f64)> = b.iter().copied().enumerate().filter(|&(_, v)| v.abs() > 0.0).collect();
    let mut coo = Coo::new(n, n);
    for &(i, vi) in &nz {
        for &(j, vj) in &nz {
            coo.push(i, j, vi * vj);
        }
    }
    coo.to_csr()
}

// ------------------------------------------------------------------ Iterative AMS / AmsP2 on MixedDiscretisation

/// Solves a 3D mixed-cell curl-curl linear system `A x = rhs` on free DOFs (order 1 `Ams` or order 2 `AmsP2`)
/// using right-preconditioned FGMRES, building the auxiliary-space hierarchy from `pc_matrix.unwrap_or(a)`.
/// Callers with PMLs should resolve their `EmPmlTable` with `is_iterative = true` so any `Rational` profile
/// safely falls back to `KappaPolynomial { kappa_max: 2.0 }`.
#[allow(clippy::too_many_arguments)]
pub fn solve_mixed_curl_curl_iterative(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    fixed_dofs: &[bool],
    fixed_edges: &[bool],
    fixed_verts: &[bool],
    fixed_quad_faces: &[bool],
    a: &CCsr,
    pc_matrix: Option<&CCsr>,
    rhs: &[C],
    opts: &IterativeCurlCurlOptions,
) -> Result<IterativeSolveReport<C>, String> {
    if ned.order == 2 {
        let ams_p2 = build_ams_p2_from_masks(mesh, ned, fixed_dofs, fixed_edges, fixed_verts, fixed_quad_faces);
        return solve_curl_curl_p2_iterative(a, pc_matrix, &ams_p2.topology, rhs, opts);
    }

    let n = a.nrows();
    if n == 0 || a.ncols() != n || rhs.len() != n {
        return Err("solve_mixed_curl_curl_iterative: matrix must be square and match rhs length".into());
    }
    let ams_p1 = build_ams_p1_from_masks(mesh, &ned.vkey, fixed_edges, fixed_verts);
    let pc_mat = pc_matrix.unwrap_or(a);
    let ams = Ams::new(pc_mat, &ams_p1.gradient_p1, &ams_p1.edge_vectors, &AmsOptions::default())
        .map_err(|e| format!("Ams (P1) setup failed: {e}"))?;
    let preconditioner_nnz = ams.nnz_storage();

    let mut gmres = Gmres::<C>::new(GmresOptions {
        tol: Tolerances { rtol: opts.rtol, atol: opts.atol, max_it: opts.max_it, ..Default::default() },
        restart: opts.restart.max(2),
        side: PcSide::Right,
        orthog: Orthogonalization::ClassicalGramSchmidtIfNeeded,
        flexible: opts.flexible,
        ..Default::default()
    });
    let mut x = vec![c(0.0, 0.0); n];
    let info = gmres.solve(a, &ams, rhs, &mut x);
    let ax = a.matvec(&x);
    let b_norm = rhs.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt().max(1e-300);
    let r_norm = ax.iter().zip(rhs).map(|(p, q)| (*p - *q).norm_sqr()).sum::<f64>().sqrt();
    let relative_residual = r_norm / b_norm;
    Ok(IterativeSolveReport {
        x,
        iterations: info.iterations,
        relative_residual,
        converged: info.reason.is_converged() && relative_residual <= (opts.rtol * 10.0).max(opts.atol),
        preconditioner_nnz,
    })
}

/// Computes `nev` interior eigenpairs of a 3D mixed-cell `Tet + Prism + Hex` cavity `K x = lambda M x` near `target`
/// iteratively via Jacobi-Davidson preconditioned by `Ams` (`order == 1`) or `AmsP2` (`order == 2`).
#[allow(clippy::too_many_arguments)]
pub fn eigs_mixed_curl_curl_ams(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    fixed_dofs: &[bool],
    fixed_edges: &[bool],
    fixed_verts: &[bool],
    fixed_quad_faces: &[bool],
    k: &CCsr,
    m: &CCsr,
    target: C,
    nev: usize,
    tol: f64,
    max_it: usize,
) -> Result<CEigPairs, String> {
    if ned.order == 2 {
        let ams_p2 = build_ams_p2_from_masks(mesh, ned, fixed_dofs, fixed_edges, fixed_verts, fixed_quad_faces);
        return eigs_curl_curl_p2_jd(k, m, &ams_p2.topology, target, nev, tol, max_it);
    }
    use femcore::linalg::jd::{jacobi_davidson, JdOptions};
    let n = k.nrows();
    if n == 0 || k.ncols() != n || m.nrows() != n || m.ncols() != n {
        return Err("eigs_mixed_curl_curl_ams: K and M must be square of identical dimension".into());
    }
    let ams_p1 = build_ams_p1_from_masks(mesh, &ned.vkey, fixed_edges, fixed_verts);
    let beta = target.norm().max(1.0);
    let spd = CCsr::axpby(c(1.0, 0.0), k, c(beta, 0.0), m);
    let ams = Ams::new(&spd, &ams_p1.gradient_p1, &ams_p1.edge_vectors, &AmsOptions::default())
        .map_err(|e| format!("Ams (P1) setup failed: {e}"))?;
    let mut opts = JdOptions::new(nev, target, true);
    opts.tol = tol;
    opts.max_it = max_it;
    opts.harmonic = false;
    opts.inner_its = 25;
    opts.m_min = (nev + 6).max(8);
    opts.m_max = (opts.m_min + 16).min(n.saturating_sub(2));
    let res = jacobi_davidson(k, Some(m), &ams, None, &opts).map_err(|e| format!("Jacobi-Davidson failed: {e}"))?;
    let vectors = (0..res.values.len()).map(|j| res.vectors.col(j).to_vec()).collect();
    Ok(CEigPairs {
        values: res.values,
        vectors,
        residuals: res.residuals,
        dim: opts.m_max,
        restarts: res.iterations,
        matvecs: res.matvecs,
        converged: res.converged,
    })
}

/// Used by `MixedLagrange` scalar nullspace reduction on `MixedDiscretisation`.
pub fn gradient_mixed_full(mesh: &MixedDiscretisation, ned: &NedMixed) -> RCsr {
    if ned.order == 1 {
        gradient_p1_full(mesh, &ned.vkey)
    } else {
        gradient_p2_full(mesh, ned)
    }
}

pub fn scalar_nullspace_count_mixed(mesh: &MixedDiscretisation, order: usize) -> usize {
    if order == 1 {
        mesh.n_points()
    } else {
        scalar_p2_tensor_count(mesh)
    }
}

#[allow(dead_code)]
fn _bfs_touch(adj: &[Vec<usize>], start: usize) -> Vec<bool> {
    let mut vis = vec![false; adj.len()];
    let mut q = VecDeque::new();
    vis[start] = true;
    q.push_back(start);
    while let Some(u) = q.pop_front() {
        for &v in &adj[u] {
            if !vis[v] {
                vis[v] = true;
                q.push_back(v);
            }
        }
    }
    vis
}
