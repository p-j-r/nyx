// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Linear elastic solid mechanics on P1/P2 triangles/quadrilaterals (plane stress or plane strain) and
//! tetrahedra/prisms/hexahedra: static solve with boundary loads, prescribed displacements, thermal
//! expansion, and locking-free Static Infinite Elements (`Rational` and `Exponential`), plus real/complex
//! eigenfrequencies, stress, and von Mises.
//! Materials may be anisotropic (a 6x6 stiffness per domain, see `femcore::material`).
//!
//! Everything is in the wafer axes (`femcore::section`): a 2D mesh is the x-z cross-section, y the axis out of it. Dof layout: all displacements
//! along the first mesh axis, then along the second (then the third): (x, z) on a cross-section, (x, y, z) in 3D. Voigt order of stresses
//! and strains with engineering shear strains: (xx, yy, zz, xy, yz, xz) in 3D; on a cross-section the Voigt order of the axes (x, z, y),
//! (xx, zz, yy, xz, zy, xy), that `section::plane_stiffness` gives the stiffness in.
//!
//! 2D: plane strain uses the (xx, zz, xz) block of the stiffness with eps_yy = 0; plane stress
//! condenses eps_yy out (sigma_yy = 0). Out-of-plane shear strains are taken as zero. The
//! thermal strain is `alpha (T - Tref)` (a tensor for anisotropic expansion); in plane strain
//! its out-of-plane part loads the in-plane stress through the stiffness, which for an isotropic
//! material is the (1 + nu) alpha dT in-plane free expansion that plane strain gives.

use crate::problem::*;
use femcore::fem::check_boundary_labels;
use femcore::fem::*;
use femcore::fem_mixed::{
    assemble_mass_mixed, eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, lagrange_ref_nodes,
    lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap,
};
use femcore::fields::{gradients_at_nodes, node_average, CellFields};
use femcore::linalg::csr::{c, Coo, RCsr as Csr, C};
use femcore::linalg::driver::{eigs_sym, EigOpts};
use femcore::material::{self, M3, M6};
use femcore::mesh::cells::{canonical_quad_face, CellType};
use femcore::mesh::mixed::{
    cell_local_edges, cell_local_quad_faces, cell_local_tri_faces, sorted_tri_face, FacetType, MixedDiscretisation,
    QUAD_EDGES, TRI_EDGES,
};
use femcore::pml::ie_exponential_stretch_1d;
use femcore::props::{LinearSolver, Materials, Pml, PmlCoord, PmlDir, TensorKind};
use femcore::vvuq::{SourceField, UniformSource};
use nyx_physics::{
    cell_fields_of, CellField, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::time::Instant;

fn check_labels(mesh: &Discretisation, bcs: &[Bc]) -> Result<(), String> {
    check_boundary_labels(
        mesh,
        bcs.iter().filter(|b| !b.kind.is_domain_condition()).flat_map(|b| b.labels.iter().copied()),
    )?;
    let doms = mesh.domain_labels();
    let bad_dom: Vec<i32> = bcs
        .iter()
        .filter(|b| b.kind.is_domain_condition())
        .flat_map(|b| b.labels.iter().copied())
        .filter(|l| !doms.contains(l))
        .collect();
    if !bad_dom.is_empty() {
        return Err(format!("domain number(s) {bad_dom:?} do not exist in this mesh (it has {doms:?})"));
    }
    Ok(())
}

pub(crate) fn mixed_check_boundary_labels(mixed: &MixedDiscretisation, bcs: &[Bc]) -> Result<(), String> {
    let mut known: Vec<i32> = mixed.facet_labels.clone();
    known.sort_unstable();
    known.dedup();
    let mut bad: Vec<i32> = bcs
        .iter()
        .filter(|b| !b.kind.is_domain_condition())
        .flat_map(|b| b.labels.iter().copied())
        .filter(|l| known.binary_search(l).is_err())
        .collect();
    bad.sort_unstable();
    bad.dedup();
    if !bad.is_empty() {
        return Err(format!("boundary number(s) {bad:?} do not exist in this mesh (it has {known:?})"));
    }
    let doms = mixed_domain_labels(mixed);
    let mut bad_dom: Vec<i32> = bcs
        .iter()
        .filter(|b| b.kind.is_domain_condition())
        .flat_map(|b| b.labels.iter().copied())
        .filter(|l| doms.binary_search(l).is_err())
        .collect();
    bad_dom.sort_unstable();
    bad_dom.dedup();
    if !bad_dom.is_empty() {
        return Err(format!("domain number(s) {bad_dom:?} do not exist in this mesh (it has {doms:?})"));
    }
    Ok(())
}

fn mixed_domain_labels(mixed: &MixedDiscretisation) -> Vec<i32> {
    let mut labels = mixed.cell_labels.clone();
    labels.sort_unstable();
    labels.dedup();
    labels
}

pub(crate) fn mixed_cell_stiffness_skipping(
    materials: &Materials,
    mixed: &MixedDiscretisation,
    skip: &[i32],
) -> Result<Vec<M6>, String> {
    let mut by_domain: BTreeMap<i32, M6> = BTreeMap::new();
    for d in mixed_domain_labels(mixed) {
        if skip.contains(&d) {
            by_domain.insert(d, [[0.0; 6]; 6]);
            continue;
        }
        let c = if let Some(st) = materials.aniso.get(&d).and_then(|a| a.stiffness.as_ref()) {
            material::stiffness_global(st, materials.aniso[&d].euler_deg)
        } else {
            let e = materials.get(d, "E").ok_or_else(|| format!("domain {d} has no value for 'E'"))?;
            let nu = materials.get(d, "nu").ok_or_else(|| format!("domain {d} has no value for 'nu'"))?;
            if nu >= 0.5 || nu <= -1.0 {
                return Err("Poisson's ratio must be between -1 and 0.5".into());
            }
            if e <= 0.0 {
                return Err("Young's modulus must be positive".into());
            }
            material::isotropic_stiffness(e, nu)
        };
        if !material::is_positive_definite(&c) {
            return Err(format!("the stiffness of domain {d} is not positive definite (check the constants)"));
        }
        by_domain.insert(d, c);
    }
    Ok(mixed.cell_labels.iter().map(|l| by_domain[l]).collect())
}

pub(crate) fn mixed_cell_viscosity_skipping(
    materials: &Materials,
    mixed: &MixedDiscretisation,
    skip: &[i32],
) -> Result<Vec<M6>, String> {
    let mut by_domain: BTreeMap<i32, M6> = BTreeMap::new();
    for d in mixed_domain_labels(mixed) {
        if skip.contains(&d) {
            by_domain.insert(d, [[0.0; 6]; 6]);
            continue;
        }
        let eta = if let Some(v) = materials.aniso.get(&d).and_then(|a| a.viscosity.as_ref()) {
            material::rotate_m6_global(v, materials.aniso[&d].euler_deg)
        } else if let Some(eta_s) = materials.get(d, "eta") {
            if eta_s < 0.0 {
                return Err(format!("viscosity 'eta' of domain {d} must be non-negative (found {eta_s})"));
            }
            if eta_s > 0.0 {
                material::isotropic_stiffness(eta_s, 0.0)
            } else {
                [[0.0; 6]; 6]
            }
        } else {
            [[0.0; 6]; 6]
        };
        by_domain.insert(d, eta);
    }
    Ok(mixed.cell_labels.iter().map(|l| by_domain[l]).collect())
}

pub(crate) fn mixed_cell_rho(
    materials: &Materials,
    mixed: &MixedDiscretisation,
    inactive: &[i32],
) -> Result<Vec<f64>, String> {
    let mut by_domain = BTreeMap::new();
    for d in mixed_domain_labels(mixed) {
        let v = if inactive.contains(&d) {
            0.0
        } else {
            materials.get(d, "rho").ok_or_else(|| format!("domain {d} has no value for 'rho'"))?
        };
        by_domain.insert(d, v);
    }
    Ok(mixed.cell_labels.iter().map(|l| by_domain[l]).collect())
}

pub(crate) fn mixed_cell_tensor(
    materials: &Materials,
    mixed: &MixedDiscretisation,
    kind: TensorKind,
    fallback: Option<f64>,
) -> Result<Vec<M3>, String> {
    let (scalar_name, pick): (&str, fn(&material::Anisotropy) -> Option<&material::Tensor>) = match kind {
        TensorKind::Conductivity => ("k", |a| a.conductivity.as_ref()),
        TensorKind::Permittivity => ("eps_r", |a| a.permittivity.as_ref()),
        TensorKind::PermittivityIm => ("eps_r_im", |a| a.permittivity_im.as_ref()),
        TensorKind::Expansion => ("alpha", |a| a.expansion.as_ref()),
    };
    let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
    for d in mixed_domain_labels(mixed) {
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
    Ok(mixed.cell_labels.iter().map(|l| by_domain[l]).collect())
}

pub(crate) fn mixed_untouched_dofs(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    ncomp: usize,
    inactive: &[i32],
) -> Vec<usize> {
    if inactive.is_empty() {
        return Vec::new();
    }
    let n1 = dof_map.n_dofs;
    let mut used = vec![false; n1];
    for c in 0..mixed.n_cells() {
        if !inactive.contains(&mixed.cell_labels[c]) {
            for &d in dof_map.cell_dofs(c) {
                used[d] = true;
            }
        }
    }
    let mut out = Vec::new();
    for k in 0..ncomp {
        out.extend((0..n1).filter(|&d| !used[d]).map(|d| d + k * n1));
    }
    out
}

/// Indices of the in-plane components (xx, zz, xz) in the solver's Voigt order on a cross-section.
const INPLANE: [usize; 3] = [0, 1, 3];

/// Constitutive matrix, row-major `ns x ns` (ns = 3 in 2D, 6 in 3D).
pub(crate) fn d_matrix(c: &M6, dim: usize, plane: Plane) -> Vec<f64> {
    if dim == 3 {
        return c.iter().flat_map(|r| r.iter().copied()).collect();
    }
    let mut d = vec![0.0; 9];
    for (a, &i) in INPLANE.iter().enumerate() {
        for (b, &j) in INPLANE.iter().enumerate() {
            d[a * 3 + b] = match plane {
                Plane::Strain => c[i][j],
                Plane::Stress => {
                    if c[2][2].abs() <= 1e-30 {
                        0.0
                    } else {
                        c[i][j] - c[i][2] * c[2][j] / c[2][2]
                    }
                }
            };
        }
    }
    d
}

/// Thermal strain vector (solver Voigt order, engineering shear) for expansion tensor `a` and
/// temperature rise `dt`.
pub(crate) fn thermal_strain(a: &M3, dt: f64) -> [f64; 6] {
    [a[0][0] * dt, a[1][1] * dt, a[2][2] * dt, 2.0 * a[0][1] * dt, 2.0 * a[1][2] * dt, 2.0 * a[0][2] * dt]
}

/// The stress (in the components that load the model: all six in 3D, xx, yy, xy in 2D) that the
/// thermal strain `eth` would produce if the body were fully constrained.
fn thermal_stress(c: &M6, dim: usize, plane: Plane, eth: &[f64; 6]) -> Vec<f64> {
    if dim == 3 {
        return (0..6).map(|i| (0..6).map(|j| c[i][j] * eth[j]).sum()).collect();
    }
    match plane {
        Plane::Strain => INPLANE.iter().map(|&i| (0..6).map(|j| c[i][j] * eth[j]).sum()).collect(),
        Plane::Stress => {
            let d = d_matrix(c, 2, Plane::Stress);
            (0..3).map(|a| (0..3).map(|b| d[a * 3 + b] * eth[INPLANE[b]]).sum()).collect()
        }
    }
}

/// Strain-displacement matrix `B[s * (d*nloc) + comp*nloc + a]` from shape gradients `grad[a*d + k]`.
pub(crate) fn b_matrix(grad: &[f64], nloc: usize, d: usize, b: &mut [f64]) {
    let ns = if d == 2 { 3 } else { 6 };
    b[..ns * d * nloc].iter_mut().for_each(|x| *x = 0.0);
    let nd = d * nloc;
    for k in 0..d {
        for a in 0..nloc {
            b[k * nd + k * nloc + a] = grad[a * d + k];
        }
    }
    let shear: &[(usize, usize)] = if d == 2 { &[(0, 1)] } else { &[(0, 1), (1, 2), (0, 2)] };
    for (r, &(p, q)) in shear.iter().enumerate() {
        let row = (d + r) * nd;
        for a in 0..nloc {
            b[row + p * nloc + a] = grad[a * d + q];
            b[row + q * nloc + a] = grad[a * d + p];
        }
    }
}

pub(crate) fn dof_layout(
    mesh: &Discretisation,
    order: usize,
    c: usize,
    scratch: &mut Vec<usize>,
    full: &mut Vec<usize>,
) {
    let n1 = mesh.n_dofs(order);
    mesh.cell_dofs(c, order, scratch);
    full.clear();
    for k in 0..mesh.dim {
        full.extend(scratch.iter().map(|&d| d + k * n1));
    }
}

/// Stiffness matrix with a (possibly anisotropic) stiffness per cell.
pub fn assemble_solid(mesh: &Discretisation, c6: &[M6], order: usize, plane: Plane, thickness: f64) -> Csr {
    let d = mesh.dim;
    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let n = d * mesh.n_dofs(order);
    assemble_parallel(mesh.n_cells(), n, n, |range, coo| {
        let (mut scratch, mut full) = (Vec::new(), Vec::new());
        let mut grad = vec![0.0; nloc * d];
        let mut b = vec![0.0; ns * nd];
        let mut db = vec![0.0; ns * nd];
        let mut ke = vec![0.0; nd * nd];
        for c in range {
            let dm = d_matrix(&c6[c], d, plane);
            ke.iter_mut().for_each(|x| *x = 0.0);
            for q in 0..re.w.len() {
                shape_gradients(&re, &geo, c, q, &mut grad);
                b_matrix(&grad, nloc, d, &mut b);
                for s in 0..ns {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for t in 0..ns {
                            acc += dm[s * ns + t] * b[t * nd + j];
                        }
                        db[s * nd + j] = acc;
                    }
                }
                for i in 0..nd {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * db[s * nd + j];
                        }
                        ke[i * nd + j] += re.w[q] * acc;
                    }
                }
            }
            let sc = geo.vol[c] * if d == 2 { thickness } else { 1.0 };
            ke.iter_mut().for_each(|x| *x *= sc);
            dof_layout(mesh, order, c, &mut scratch, &mut full);
            coo.push_block(&full, &full, &ke);
        }
    })
}

// ------------------------------------------------------------------ Static Infinite Elements (`ResolvedIeDomains`)

#[derive(Clone, Copy, Debug)]
pub(crate) struct IeJacobian {
    pub j_inv: [[f64; 3]; 3],
    pub det: f64,
    pub r_tilde: f64,
}

impl IeJacobian {
    #[inline]
    fn identity() -> Self {
        Self { j_inv: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], det: 1.0, r_tilde: 1.0 }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ResolvedIeDomain {
    pub ie: InfiniteElement,
    pub pml: Pml,
    pub bounds: ([f64; 3], [f64; 3]),
}

#[derive(Clone, Debug, Default)]
pub(crate) struct ResolvedIeDomains {
    pub domains: BTreeMap<i32, ResolvedIeDomain>,
    pub diagnostics: Vec<String>,
}

fn permute_real_wafer_to_2d(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    const P: [usize; 3] = [0, 2, 1];
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = m[P[i]][P[j]];
        }
    }
    out
}

fn normalize_ie_axis(axis: [f64; 3], dim: usize) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if n > 1e-30 {
        [axis[0] / n, axis[1] / n, axis[2] / n]
    } else if dim == 2 {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn infer_ie_curvilinear_extents(pml: &mut Pml, pts_wafer: &[[f64; 3]], dim: usize) {
    if pts_wafer.is_empty() {
        return;
    }
    match &mut pml.coord {
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
            let a_hat = normalize_ie_axis(*axis, dim);
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

/// Evaluates `(s_u, inv_su, du_tilde)` along a 1D IE ray `t in [0, 1]` of physical thickness `length`.
/// Computes `inv_su` directly without dividing by `s_u`, so at `t = 1` (`q = 0`) for `Rational`
/// `inv_su == 0.0` with zero `Inf`/`NaN`.
#[inline]
fn eval_ie_ray_1d(ie: &InfiniteElement, t: f64, length: f64, default_pole: f64) -> (f64, f64, f64, f64) {
    let t = t.clamp(0.0, 1.0);
    let len = length.max(1e-300);
    let d_pole = ie.pole_dist.unwrap_or(default_pole).max(1e-30);
    match ie.profile {
        IeProfile::Rational => {
            let q = (1.0 - t).clamp(0.0, 1.0);
            let q_safe = q.max(1e-12);
            let s_u = d_pole / (len * q_safe * q_safe);
            let inv_su = (len * q * q) / d_pole;
            let du = (d_pole * t) / q_safe;
            (s_u, inv_su, du, d_pole)
        }
        IeProfile::Exponential => {
            let (s_u, du) = ie_exponential_stretch_1d(t, len, ie.decay);
            (s_u, 1.0 / s_u, du, d_pole)
        }
    }
}

impl ResolvedIeDomains {
    pub(crate) fn from_mixed(mixed: &MixedDiscretisation, spec: &Problem) -> Self {
        if spec.infinite_elements.is_empty() {
            return Self::default();
        }
        let d = mixed.dim;
        let mut dom_verts: BTreeMap<i32, BTreeSet<u32>> = BTreeMap::new();
        for ci in 0..mixed.n_cells() {
            let dom = mixed.cell_labels[ci];
            if spec.infinite_elements.contains_key(&dom) && !spec.inactive.contains(&dom) {
                dom_verts.entry(dom).or_default().extend(mixed.cell(ci).iter().copied());
            }
        }
        let mut domains = BTreeMap::new();
        let mut diagnostics = Vec::new();
        for (&dom, ie) in &spec.infinite_elements {
            if spec.inactive.contains(&dom) || !ie.is_active() {
                continue;
            }
            let Some(vs) = dom_verts.get(&dom) else { continue };
            let pts: Vec<[f64; 3]> =
                vs.iter().map(|&v| femcore::section::wafer_point(mixed.point(v as usize), d)).collect();
            let mut pml = ie.pml.clone();
            infer_ie_curvilinear_extents(&mut pml, &pts, d);
            if !pml.is_active() {
                continue;
            }
            let bounds = if let Some(b) = pml.bounds {
                b
            } else {
                let mut lo = [f64::INFINITY; 3];
                let mut hi = [f64::NEG_INFINITY; 3];
                for p in &pts {
                    for k in 0..3 {
                        lo[k] = lo[k].min(p[k]);
                        hi[k] = hi[k].max(p[k]);
                    }
                }
                (lo, hi)
            };
            let prof_name = match ie.profile {
                IeProfile::Rational => "rational [R_in, R_out) -> [R_in, inf)",
                IeProfile::Exponential => "exponential + outer 1st-order spring Robin",
            };
            diagnostics.push(format!("static infinite element on domain {dom}: {prof_name}"));
            domains.insert(dom, ResolvedIeDomain { ie: ie.clone(), pml, bounds });
        }
        Self { domains, diagnostics }
    }

    #[inline]
    pub(crate) fn contains(&self, dom: i32) -> bool {
        self.domains.contains_key(&dom)
    }

    pub(crate) fn has_rational(&self) -> bool {
        self.domains.values().any(|rd| rd.ie.profile == IeProfile::Rational)
    }

    pub(crate) fn has_exponential_robin(&self) -> bool {
        self.domains.values().any(|rd| rd.ie.profile == IeProfile::Exponential && !rd.ie.clamp_outer)
    }

    /// Evaluates the Infinite Element Jacobian inside cell `iso`, using nodal-radius isoparametric
    /// interpolation on uncurved (`!iso.curved`) Spherical/Cylindrical shells so flat polyhedral
    /// outer facets map conformingly to `t = 1` (`R_tilde -> infinity`) without secant gap error.
    pub(crate) fn eval_cell_jacobian(
        &self,
        dom: i32,
        iso: &IsoparametricCell,
        jg: &femcore::fem_mixed::JacobianEval,
        sn: &[f64],
        sdn: &[[f64; 3]],
        dim: usize,
    ) -> IeJacobian {
        let x_w = femcore::section::wafer_point(&jg.x[..dim], dim);
        let Some(rd) = self.domains.get(&dom) else {
            return IeJacobian::identity();
        };
        if !iso.curved {
            if let PmlCoord::Spherical { center, r_in, thickness } = rd.pml.coord {
                if thickness > 0.0 {
                    let mut u = 0.0;
                    let mut g_m = [0.0; 3];
                    for (b, node) in iso.nodes.iter().enumerate() {
                        let pw = femcore::section::wafer_point(&node[..dim], dim);
                        let dx = [pw[0] - center[0], pw[1] - center[1], pw[2] - center[2]];
                        let rb = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                        u += sn[b] * rb;
                        let gp = jg.map_grad(sdn[b]);
                        for k in 0..dim {
                            g_m[k] += rb * gp[k];
                        }
                    }
                    if u < r_in {
                        return IeJacobian::identity();
                    }
                    let g_w = if dim == 2 { [g_m[0], 0.0, g_m[1]] } else { g_m };
                    let r_vec = [x_w[0] - center[0], x_w[1] - center[1], x_w[2] - center[2]];
                    let c1 = r_vec[0] * g_w[0] + r_vec[1] * g_w[1] + r_vec[2] * g_w[2];
                    let r_scale = thickness.max(r_in.abs()).max(1e-30);
                    if u > 1e-15 * r_scale && c1 > 1e-15 * u {
                        let t = ((u - r_in) / thickness).clamp(0.0, 1.0);
                        let default_pole = if r_in > 0.0 { r_in } else { thickness };
                        let (s_r, inv_sr, du_r, d_pole) = eval_ie_ray_1d(&rd.ie, t, thickness, default_pole);
                        let r_tilde = (r_in + du_r).max(1e-12 * r_scale);
                        let q_rem = (1.0 - t).clamp(0.0, 1.0);
                        let (gamma, inv_gamma) = match rd.ie.profile {
                            IeProfile::Rational => {
                                let q_safe = q_rem.max(1e-12);
                                let num = r_in * q_safe + d_pole * t;
                                (num / (u * q_safe), (u * q_rem) / (r_in * q_rem + d_pole * t).max(1e-300))
                            }
                            IeProfile::Exponential => (r_tilde / u, u / r_tilde.max(1e-300)),
                        };
                        let ratio = c1 / u;
                        let sigma = (1.0 - ratio) * gamma + ratio * s_r;
                        let inv_sigma = if rd.ie.profile == IeProfile::Rational && q_rem <= 1e-6 {
                            inv_sr / ratio.max(1e-300)
                        } else {
                            1.0 / sigma.max(1e-300)
                        };
                        let coeff = (inv_gamma - inv_sigma) / c1;
                        let mut j_inv_w = [[0.0; 3]; 3];
                        for p in 0..3 {
                            for q_idx in 0..3 {
                                let diag = if p == q_idx { inv_gamma } else { 0.0 };
                                j_inv_w[p][q_idx] = diag - coeff * r_vec[p] * g_w[q_idx];
                            }
                        }
                        if dim == 2 {
                            j_inv_w[1][1] = 1.0;
                            return IeJacobian {
                                j_inv: permute_real_wafer_to_2d(j_inv_w),
                                det: gamma * sigma,
                                r_tilde,
                            };
                        }
                        return IeJacobian { j_inv: j_inv_w, det: gamma * gamma * sigma, r_tilde };
                    }
                }
            }
        }
        self.eval_jacobian(dom, x_w, dim)
    }

    pub(crate) fn eval_jacobian(&self, dom: i32, x_w: [f64; 3], dim: usize) -> IeJacobian {
        let Some(rd) = self.domains.get(&dom) else {
            return IeJacobian::identity();
        };
        let (lo, hi) = rd.bounds;
        match rd.pml.coord {
            PmlCoord::Cartesian => {
                let dirs = if dim == 2 { [rd.pml.x, PmlDir::Off, rd.pml.z] } else { [rd.pml.x, rd.pml.y, rd.pml.z] };
                let mut s_w = [1.0; 3];
                let mut inv_w = [1.0; 3];
                let mut x_tilde = x_w;
                let mut r_scale = 0.0_f64;
                for a in 0..3 {
                    let len = (hi[a] - lo[a]).max(1e-300);
                    r_scale = r_scale.max(len).max(lo[a].abs()).max(hi[a].abs());
                    let (t, default_pole) = match dirs[a] {
                        PmlDir::Off => continue,
                        PmlDir::Plus => {
                            let pole = if lo[a] > 1e-12 * len { lo[a] } else { len };
                            (((x_w[a] - lo[a]) / len).clamp(0.0, 1.0), pole)
                        }
                        PmlDir::Minus => {
                            let pole = if hi[a] < -1e-12 * len { -hi[a] } else { len };
                            (((hi[a] - x_w[a]) / len).clamp(0.0, 1.0), pole)
                        }
                    };
                    let (su, inv_su, du, _) = eval_ie_ray_1d(&rd.ie, t, len, default_pole);
                    s_w[a] = su;
                    inv_w[a] = inv_su;
                    x_tilde[a] = match dirs[a] {
                        PmlDir::Plus => lo[a] + du,
                        PmlDir::Minus => hi[a] - du,
                        PmlDir::Off => x_w[a],
                    };
                }
                let mut j_inv_w = [[0.0; 3]; 3];
                for a in 0..3 {
                    j_inv_w[a][a] = inv_w[a];
                }
                let det = if dim == 2 { s_w[0] * s_w[2] } else { s_w[0] * s_w[1] * s_w[2] };
                let r_tilde = if dim == 2 {
                    x_tilde[0].hypot(x_tilde[2])
                } else {
                    (x_tilde[0] * x_tilde[0] + x_tilde[1] * x_tilde[1] + x_tilde[2] * x_tilde[2]).sqrt()
                }
                .max(1e-12 * r_scale.max(1e-30));
                let j_inv = if dim == 2 { permute_real_wafer_to_2d(j_inv_w) } else { j_inv_w };
                IeJacobian { j_inv, det, r_tilde }
            }
            PmlCoord::Cylindrical { center, axis, r_in, r_thickness, radial_on, axial, z_in, z_thickness } => {
                let a_hat = normalize_ie_axis(axis, dim);
                let d_vec = [x_w[0] - center[0], x_w[1] - center[1], x_w[2] - center[2]];
                let zeta = d_vec[0] * a_hat[0] + d_vec[1] * a_hat[1] + d_vec[2] * a_hat[2];
                let r_vec = [d_vec[0] - zeta * a_hat[0], d_vec[1] - zeta * a_hat[1], d_vec[2] - zeta * a_hat[2]];
                let r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
                let r_scale = r_thickness.max(r_in.abs()).max(z_thickness).max(1e-30);

                let (mut s_r, mut inv_sr, mut s_phi, mut inv_sphi, mut r_tilde_cyl) = (1.0, 1.0, 1.0, 1.0, r);
                if radial_on && r_thickness > 0.0 && r >= r_in {
                    let t_r = ((r - r_in) / r_thickness).clamp(0.0, 1.0);
                    let default_pole = if r_in > 0.0 { r_in } else { r_thickness };
                    let (sr_v, inv_sr_v, du_r, d_pole) = eval_ie_ray_1d(&rd.ie, t_r, r_thickness, default_pole);
                    s_r = sr_v;
                    inv_sr = inv_sr_v;
                    r_tilde_cyl = r_in + du_r;
                    if r > 1e-15 * r_scale {
                        match rd.ie.profile {
                            IeProfile::Rational => {
                                let q = (1.0 - t_r).clamp(0.0, 1.0);
                                let q_safe = q.max(1e-12);
                                let num = r_in * q_safe + d_pole * t_r;
                                s_phi = num / (r * q_safe);
                                inv_sphi = (r * q) / (r_in * q + d_pole * t_r).max(1e-300);
                            }
                            IeProfile::Exponential => {
                                s_phi = r_tilde_cyl / r;
                                inv_sphi = r / r_tilde_cyl.max(1e-300);
                            }
                        }
                    } else {
                        s_phi = s_r;
                        inv_sphi = inv_sr;
                    }
                }

                let (mut s_zeta, mut inv_szeta, mut zeta_tilde) = (1.0, 1.0, zeta);
                if dim == 3 && axial != PmlDir::Off && z_thickness > 0.0 {
                    let raw_tz = match axial {
                        PmlDir::Plus => (zeta - z_in) / z_thickness,
                        PmlDir::Minus => (z_in - zeta) / z_thickness,
                        PmlDir::Off => 0.0,
                    };
                    if raw_tz >= 0.0 {
                        let t_z = raw_tz.clamp(0.0, 1.0);
                        let default_pole = if z_in.abs() > 1e-12 * z_thickness { z_in.abs() } else { z_thickness };
                        let (sz_v, inv_sz_v, du_z, _) = eval_ie_ray_1d(&rd.ie, t_z, z_thickness, default_pole);
                        s_zeta = sz_v;
                        inv_szeta = inv_sz_v;
                        zeta_tilde = if axial == PmlDir::Plus { z_in + du_z } else { z_in - du_z };
                    }
                }

                let r_hat = if r > 1e-30 { [r_vec[0] / r, r_vec[1] / r, r_vec[2] / r] } else { [0.0; 3] };
                let mut j_inv_w = [[0.0; 3]; 3];
                let dinv_r = inv_sr - inv_sphi;
                for p in 0..3 {
                    for q in 0..3 {
                        let id_pq = if p == q { 1.0 } else { 0.0 };
                        let pa = a_hat[p] * a_hat[q];
                        let pr = r_hat[p] * r_hat[q];
                        j_inv_w[p][q] = inv_sphi * (id_pq - pa) + dinv_r * pr + inv_szeta * pa;
                    }
                }
                let det = if dim == 2 { s_r * s_phi } else { s_r * s_phi * s_zeta };
                let r_tilde = if dim == 2 { r_tilde_cyl } else { r_tilde_cyl.hypot(zeta_tilde) }.max(1e-12 * r_scale);
                let j_inv = if dim == 2 { permute_real_wafer_to_2d(j_inv_w) } else { j_inv_w };
                IeJacobian { j_inv, det, r_tilde }
            }
            PmlCoord::Spherical { center, r_in, thickness } => {
                if thickness <= 0.0 {
                    return IeJacobian::identity();
                }
                let r_vec = [x_w[0] - center[0], x_w[1] - center[1], x_w[2] - center[2]];
                let big_r = (r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2]).sqrt();
                if big_r < r_in {
                    return IeJacobian::identity();
                }
                let r_scale = thickness.max(r_in.abs()).max(1e-30);
                let t = ((big_r - r_in) / thickness).clamp(0.0, 1.0);
                let default_pole = if r_in > 0.0 { r_in } else { thickness };
                let (s_r, inv_sr, du_r, d_pole) = eval_ie_ray_1d(&rd.ie, t, thickness, default_pole);
                let r_tilde = (r_in + du_r).max(1e-12 * r_scale);
                let (s_om, inv_som) = if big_r > 1e-15 * r_scale {
                    match rd.ie.profile {
                        IeProfile::Rational => {
                            let q = (1.0 - t).clamp(0.0, 1.0);
                            let q_safe = q.max(1e-12);
                            let num = r_in * q_safe + d_pole * t;
                            (num / (big_r * q_safe), (big_r * q) / (r_in * q + d_pole * t).max(1e-300))
                        }
                        IeProfile::Exponential => (r_tilde / big_r, big_r / r_tilde.max(1e-300)),
                    }
                } else {
                    (s_r, inv_sr)
                };
                let r_hat =
                    if big_r > 1e-30 { [r_vec[0] / big_r, r_vec[1] / big_r, r_vec[2] / big_r] } else { [0.0; 3] };
                let mut j_inv_w = [[0.0; 3]; 3];
                let dinv = inv_sr - inv_som;
                for p in 0..3 {
                    for q in 0..3 {
                        let id_pq = if p == q { 1.0 } else { 0.0 };
                        let pr = r_hat[p] * r_hat[q];
                        j_inv_w[p][q] = inv_som * id_pq + dinv * pr;
                    }
                }
                if dim == 2 {
                    j_inv_w[1][1] = 1.0;
                    return IeJacobian { j_inv: permute_real_wafer_to_2d(j_inv_w), det: s_r * s_om, r_tilde };
                }
                IeJacobian { j_inv: j_inv_w, det: s_r * s_om * s_om, r_tilde }
            }
        }
    }

    pub(crate) fn outer_clamped_dofs(
        &self,
        mixed: &MixedDiscretisation,
        dof_map: &MixedLagrangeDofMap,
        spec: &Problem,
    ) -> Vec<usize> {
        if self.domains.values().all(|rd| !rd.ie.clamp_outer) {
            return Vec::new();
        }
        let mut protected = BTreeSet::new();
        for bc in &spec.bcs {
            if matches!(bc.kind, BcKind::Fixed { .. }) {
                protected.extend(bc.labels.iter().copied());
            }
        }
        let d = mixed.dim;
        crate::pml::outer_faces_mixed_where(mixed, dof_map, &spec.inactive, &protected, |dom, pts_w| {
            self.domains.get(&dom).is_some_and(|rd| {
                rd.ie.clamp_outer && crate::pml::is_outer_termination_face(&rd.pml, rd.bounds, pts_w, d, false)
            })
        })
    }
}

/// Evaluates `(x, ds, n_hat)` on an arbitrary face defined by corner vertex indices `fv` (`Edge2`, `Tri3`, `Quad4`)
/// regardless of whether the face appears in `mixed.facets`.
pub(crate) fn eval_raw_face_jacobian(
    mixed: &MixedDiscretisation,
    ft: FacetType,
    fv: &[u32],
    xi: &[f64; 2],
) -> ([f64; 3], f64, [f64; 3]) {
    let corners: Vec<[f64; 3]> = fv.iter().map(|&v| mixed.point3(v as usize)).collect();
    match ft {
        FacetType::Edge2 => {
            let eid = mixed.edge_index(fv[0], fv[1]).expect("face edge");
            let (space, nodes) = if let Some(mid) = mixed.edge_node(eid) {
                (LagrangeSpace::Order2Serendipity, vec![corners[0], corners[1], mid])
            } else {
                (LagrangeSpace::Order1, corners)
            };
            let (n, dn) = facet_lagrange_shape(ft, space, xi);
            let mut x = [0.0; 3];
            let mut t = [0.0; 2];
            for (a, p) in nodes.iter().enumerate() {
                x[0] += n[a] * p[0];
                x[1] += n[a] * p[1];
                t[0] += dn[a][0] * p[0];
                t[1] += dn[a][0] * p[1];
            }
            let ds = t[0].hypot(t[1]);
            let n_hat = if ds > 0.0 { [t[1] / ds, -t[0] / ds, 0.0] } else { [0.0; 3] };
            (x, ds, n_hat)
        }
        FacetType::Tri3 => {
            let e01 = mixed.edge_index(fv[0], fv[1]).unwrap();
            let e02 = mixed.edge_index(fv[0], fv[2]).unwrap();
            let e12 = mixed.edge_index(fv[1], fv[2]).unwrap();
            let any_curved =
                mixed.edge_node(e01).is_some() || mixed.edge_node(e02).is_some() || mixed.edge_node(e12).is_some();
            let (space, nodes) = if any_curved {
                let mid = |eid: usize, i: usize, j: usize| {
                    mixed.edge_node(eid).unwrap_or([
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
            let cx = [t0[1] * t1[2] - t0[2] * t1[1], t0[2] * t1[0] - t0[0] * t1[2], t0[0] * t1[1] - t0[1] * t1[0]];
            let ds = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            let n_hat = if ds > 0.0 { [cx[0] / ds, cx[1] / ds, cx[2] / ds] } else { [0.0; 3] };
            (x, ds, n_hat)
        }
        FacetType::Quad4 => {
            let qf_id = mixed.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
            let eids = [
                mixed.edge_index(fv[0], fv[1]).unwrap(),
                mixed.edge_index(fv[1], fv[2]).unwrap(),
                mixed.edge_index(fv[2], fv[3]).unwrap(),
                mixed.edge_index(fv[3], fv[0]).unwrap(),
            ];
            let has_qfn = mixed.quad_face_node(qf_id).is_some();
            let any_curved = has_qfn || eids.iter().any(|&e| mixed.edge_node(e).is_some());
            let (space, nodes) = if any_curved {
                let mut ns = corners.clone();
                for (k, &(i, j)) in QUAD_EDGES.iter().enumerate() {
                    ns.push(mixed.edge_node(eids[k]).unwrap_or([
                        0.5 * (corners[i][0] + corners[j][0]),
                        0.5 * (corners[i][1] + corners[j][1]),
                        0.5 * (corners[i][2] + corners[j][2]),
                    ]));
                }
                if has_qfn {
                    ns.push(mixed.quad_face_node(qf_id).unwrap());
                    (LagrangeSpace::Order2Tensor, ns)
                } else {
                    (LagrangeSpace::Order2Serendipity, ns)
                }
            } else {
                (LagrangeSpace::Order1, corners)
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
            let cx = [t0[1] * t1[2] - t0[2] * t1[1], t0[2] * t1[0] - t0[0] * t1[2], t0[0] * t1[1] - t0[1] * t1[0]];
            let ds = (cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]).sqrt();
            let n_hat = if ds > 0.0 { [cx[0] / ds, cx[1] / ds, cx[2] / ds] } else { [0.0; 3] };
            (x, ds, n_hat)
        }
    }
}

pub(crate) fn raw_face_dofs(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    ft: FacetType,
    fv: &[u32],
    out: &mut Vec<usize>,
) {
    out.clear();
    for &v in fv {
        out.push(v as usize);
    }
    if dof_map.space == LagrangeSpace::Order1 {
        return;
    }
    let nv = dof_map.n_vertices;
    match ft {
        FacetType::Edge2 => out.push(nv + mixed.edge_index(fv[0], fv[1]).unwrap()),
        FacetType::Tri3 => {
            for &(i, j) in &TRI_EDGES {
                out.push(nv + mixed.edge_index(fv[i], fv[j]).unwrap());
            }
        }
        FacetType::Quad4 => {
            for &(i, j) in &QUAD_EDGES {
                out.push(nv + mixed.edge_index(fv[i], fv[j]).unwrap());
            }
            if dof_map.space == LagrangeSpace::Order2Tensor {
                let qf = mixed.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap();
                out.push(nv + dof_map.n_edges + qf);
            }
        }
    }
}

/// Collects exterior termination faces of active exponential Infinite Element domains (when `clamp_outer == false`).
pub(crate) fn collect_exponential_ie_faces(
    mixed: &MixedDiscretisation,
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
    skip_domains: &[i32],
) -> Vec<(usize, FacetType, Vec<u32>)> {
    let d = mixed.dim;
    let mut faces: Vec<(usize, FacetType, Vec<u32>)> = Vec::new();
    let skip_cell = |ci: usize| {
        let dom = mixed.cell_labels[ci];
        inactive.contains(&dom) || skip_domains.contains(&dom)
    };
    if d == 2 {
        let mut adj: HashMap<(u32, u32), (u32, usize, [u32; 2])> = HashMap::new();
        for ci in 0..mixed.n_cells() {
            if skip_cell(ci) {
                continue;
            }
            let vs = mixed.cell(ci);
            for &(i, j) in cell_local_edges(mixed.cell_type(ci)) {
                let (a, b) = (vs[i], vs[j]);
                let e = adj.entry((a.min(b), a.max(b))).or_insert((0, ci, [a, b]));
                e.0 += 1;
            }
        }
        for (_, (cnt, ci, fv)) in adj {
            if cnt != 1 {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let Some(rd) = ie_domains.domains.get(&dom) else { continue };
            if rd.ie.profile != IeProfile::Exponential || rd.ie.clamp_outer {
                continue;
            }
            let pts_w = fv.map(|v| femcore::section::wafer_point(mixed.point(v as usize), 2));
            if crate::pml::is_outer_termination_face(&rd.pml, rd.bounds, &pts_w, 2, false) {
                faces.push((ci, FacetType::Edge2, fv.to_vec()));
            }
        }
    } else {
        let mut tri_adj: HashMap<[u32; 3], (u32, usize, [u32; 3])> = HashMap::new();
        let mut quad_adj: HashMap<[u32; 4], (u32, usize, [u32; 4])> = HashMap::new();
        for ci in 0..mixed.n_cells() {
            if skip_cell(ci) {
                continue;
            }
            let vs = mixed.cell(ci);
            let ct = mixed.cell_type(ci);
            for lf in cell_local_tri_faces(ct) {
                let raw = [vs[lf[0]], vs[lf[1]], vs[lf[2]]];
                let e = tri_adj.entry(sorted_tri_face(raw)).or_insert((0, ci, raw));
                e.0 += 1;
            }
            for lq in cell_local_quad_faces(ct) {
                let raw = [vs[lq[0]], vs[lq[1]], vs[lq[2]], vs[lq[3]]];
                let e = quad_adj.entry(canonical_quad_face(raw)).or_insert((0, ci, raw));
                e.0 += 1;
            }
        }
        for (_, (cnt, ci, fv)) in tri_adj {
            if cnt != 1 {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let Some(rd) = ie_domains.domains.get(&dom) else { continue };
            if rd.ie.profile != IeProfile::Exponential || rd.ie.clamp_outer {
                continue;
            }
            let pts_w = fv.map(|v| femcore::section::wafer_point(mixed.point(v as usize), 3));
            if crate::pml::is_outer_termination_face(&rd.pml, rd.bounds, &pts_w, 3, false) {
                faces.push((ci, FacetType::Tri3, fv.to_vec()));
            }
        }
        for (_, (cnt, ci, fv)) in quad_adj {
            if cnt != 1 {
                continue;
            }
            let dom = mixed.cell_labels[ci];
            let Some(rd) = ie_domains.domains.get(&dom) else { continue };
            if rd.ie.profile != IeProfile::Exponential || rd.ie.clamp_outer {
                continue;
            }
            let pts_w = fv.map(|v| femcore::section::wafer_point(mixed.point(v as usize), 3));
            if crate::pml::is_outer_termination_face(&rd.pml, rd.bounds, &pts_w, 3, false) {
                faces.push((ci, FacetType::Quad4, fv.to_vec()));
            }
        }
    }
    faces
}

/// Assembles the outer 1st-order asymptotic elastic spring Robin matrix on exterior termination
/// facets of `IeProfile::Exponential` domains (`t_tilde = -(1 / R_tilde) [G I + (M_p - G) n_tilde n_tilde^T] u`).
pub(crate) fn assemble_exponential_ie_robin(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
    thickness_2d: f64,
) -> Csr {
    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    if !ie_domains.has_exponential_robin() {
        return Csr::zeros(n, n);
    }
    let faces = collect_exponential_ie_faces(mixed, ie_domains, inactive, &[]);
    if faces.is_empty() {
        return Csr::zeros(n, n);
    }

    let deg = 2 * dof_map.space.order() + 4;
    let thk = if d == 2 { thickness_2d } else { 1.0 };
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();

    for (ci, ft, fv) in faces {
        let dom = mixed.cell_labels[ci];
        let (bulk, shear) = crate::lrb::voigt_moduli(&c6[ci]);
        let p_mod = bulk + 4.0 * shear / 3.0;
        raw_face_dofs(mixed, dof_map, ft, &fv, &mut dofs);
        let nloc = dofs.len();
        let rule = facet_rule_mixed(ft, deg);

        let mut blocks = vec![vec![0.0; nloc * nloc]; d * d];
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (x_m, ds, n_hat) = eval_raw_face_jacobian(mixed, ft, &fv, xi);
            if ds <= 0.0 {
                continue;
            }
            let x_w = femcore::section::wafer_point(&x_m[..d], d);
            let ie_j = ie_domains.eval_jacobian(dom, x_w, d);
            // Nanson's surface area/normal relation: d\tilde{\Gamma} \hat{\tilde{n}} = |det J_s| J_s^{-T} \hat{n}_phys d\Gamma_phys
            let mut n_tilde = [0.0; 3];
            for i in 0..d {
                let mut s = 0.0;
                for k in 0..d {
                    s += ie_j.j_inv[k][i] * n_hat[k];
                }
                n_tilde[i] = ie_j.det.abs() * s;
            }
            let area_scale = (0..d).map(|i| n_tilde[i] * n_tilde[i]).sum::<f64>().sqrt();
            if area_scale <= 0.0 {
                continue;
            }
            let nn: [f64; 3] = std::array::from_fn(|i| if i < d { n_tilde[i] / area_scale } else { 0.0 });
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            let w_ds = wq * ds * area_scale * thk / ie_j.r_tilde.max(1e-30);
            for i in 0..d {
                for j in 0..d {
                    let dij = (if i == j { shear } else { 0.0 }) + (p_mod - shear) * nn[i] * nn[j];
                    if dij == 0.0 {
                        continue;
                    }
                    let blk = &mut blocks[i * d + j];
                    for a in 0..nloc {
                        let wa = w_ds * dij * nv[a];
                        for b in 0..nloc {
                            blk[a * nloc + b] += wa * nv[b];
                        }
                    }
                }
            }
        }

        for i in 0..d {
            let ri: Vec<usize> = dofs.iter().map(|&x| x + i * n1).collect();
            for j in 0..d {
                let blk = &blocks[i * d + j];
                if blk.iter().any(|&v| v != 0.0) {
                    let cj: Vec<usize> = dofs.iter().map(|&x| x + j * n1).collect();
                    coo.push_block(&ri, &cj, blk);
                }
            }
        }
    }
    coo.to_csr()
}

/// Native mixed-cell elastic stiffness matrix (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`).
pub(crate) fn assemble_solid_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    plane: Plane,
    thickness: f64,
) -> Csr {
    assemble_solid_mixed_with_ie(mixed, dof_map, c6, plane, thickness, &ResolvedIeDomains::default(), &[])
}

fn assemble_solid_mixed_with_ie(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    plane: Plane,
    thickness: f64,
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
) -> Csr {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    let ns = if d == 2 { 3 } else { 6 };
    let thk = if d == 2 { thickness } else { 1.0 };

    let k_vol = assemble_parallel(mixed.n_cells(), n, n, |range, coo| {
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();
        let mut grad = Vec::new();
        let mut b = Vec::new();
        let mut db = Vec::new();
        let mut ke = Vec::new();
        let mut full = Vec::new();

        for ci in range {
            let dom = mixed.cell_labels[ci];
            if inactive.contains(&dom) {
                continue;
            }
            let ct = mixed.cell_type(ci);
            let iso = IsoparametricCell::from_mesh(mixed, ci);
            let in_ie = ie_domains.contains(dom);
            let deg = if in_ie {
                10
            } else if iso.curved {
                2 * order + 4
            } else {
                2 * order + 2
            };
            let rule = quad_rule_mixed(ct, deg);
            let cdofs = dof_map.cell_dofs(ci);
            let nloc = cdofs.len();
            let nd = d * nloc;
            let dm = d_matrix(&c6[ci], d, plane);

            grad.resize(nloc * d, 0.0);
            b.resize(ns * nd, 0.0);
            db.resize(ns * nd, 0.0);
            ke.clear();
            ke.resize(nd * nd, 0.0);

            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);

                let w_eff = if in_ie {
                    let ie_j = ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        for i in 0..d {
                            let mut acc = 0.0;
                            for k in 0..d {
                                acc += ie_j.j_inv[k][i] * gp[k];
                            }
                            grad[a * d + i] = acc;
                        }
                    }
                    wq * jg.det.abs() * ie_j.det * thk
                } else {
                    for a in 0..nloc {
                        let gp = jg.map_grad(dxi[a]);
                        grad[a * d..a * d + d].copy_from_slice(&gp[..d]);
                    }
                    wq * jg.det.abs() * thk
                };

                b_matrix(&grad, nloc, d, &mut b);
                for s in 0..ns {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for t in 0..ns {
                            acc += dm[s * ns + t] * b[t * nd + j];
                        }
                        db[s * nd + j] = acc;
                    }
                }
                for i in 0..nd {
                    for j in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * db[s * nd + j];
                        }
                        ke[i * nd + j] += w_eff * acc;
                    }
                }
            }

            full.clear();
            for k in 0..d {
                full.extend(cdofs.iter().map(|&gd| gd + k * n1));
            }
            coo.push_block(&full, &full, &ke);
        }
    });

    if ie_domains.has_exponential_robin() {
        let k_rob = assemble_exponential_ie_robin(mixed, dof_map, c6, ie_domains, inactive, thickness);
        if k_rob.nnz() > 0 {
            return k_vol.add_scaled(&k_rob, 1.0);
        }
    }
    k_vol
}

/// Load vector of the thermal strain `alpha (T - Tref)`; `t` is the temperature at the scalar
/// dofs of the mesh and `alpha` the expansion tensor of every cell.
#[allow(clippy::too_many_arguments)]
pub fn thermal_load(
    mesh: &Discretisation,
    c6: &[M6],
    alpha: &[M3],
    t_field: &[f64],
    tref: f64,
    order: usize,
    plane: Plane,
    thickness: f64,
) -> Vec<f64> {
    let d = mesh.dim;
    let re = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let n1 = mesh.n_dofs(order);
    let nc = mesh.n_cells();
    // cells in blocks on several threads, each with its own load vector, added up at the end
    let nt = if nc < 4000 { 1 } else { femcore::par::threads_up_to(8) };
    let chunk = nc.div_ceil(nt).max(1);
    let (re, geo) = (&re, &geo);
    let parts: Vec<Vec<f64>> = (0..nt)
        .into_par_iter()
        .map(|t| {
            let mut f = vec![0.0; d * n1];
            let (mut scratch, mut full) = (Vec::new(), Vec::new());
            let mut grad = vec![0.0; nloc * d];
            let mut b = vec![0.0; ns * nd];
            for c in t * chunk..((t + 1) * chunk).min(nc) {
                mesh.cell_dofs(c, order, &mut scratch);
                let mut fe = vec![0.0; nd];
                for q in 0..re.w.len() {
                    let tq: f64 = (0..nloc).map(|a| re.n[q][a] * t_field[scratch[a]]).sum();
                    let sig = thermal_stress(&c6[c], d, plane, &thermal_strain(&alpha[c], tq - tref));
                    shape_gradients(re, geo, c, q, &mut grad);
                    b_matrix(&grad, nloc, d, &mut b);
                    for i in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * sig[s];
                        }
                        fe[i] += re.w[q] * acc;
                    }
                }
                let sc = geo.vol[c] * if d == 2 { thickness } else { 1.0 };
                dof_layout(mesh, order, c, &mut scratch, &mut full);
                for (i, &dof) in full.iter().enumerate() {
                    f[dof] += fe[i] * sc;
                }
            }
            f
        })
        .collect();
    let mut parts = parts.into_iter();
    let mut f = parts.next().unwrap_or_else(|| vec![0.0; d * n1]);
    for p in parts {
        f.iter_mut().zip(&p).for_each(|(a, b)| *a += b);
    }
    f
}

#[allow(clippy::too_many_arguments)]
fn thermal_load_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    alpha: &[M3],
    t_field: &[f64],
    tref: f64,
    plane: Plane,
    thickness: f64,
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
) -> Vec<f64> {
    let d = mixed.dim;
    let space = dof_map.space;
    let order = space.order();
    let n1 = dof_map.n_dofs;
    let nc = mixed.n_cells();
    let ns = if d == 2 { 3 } else { 6 };
    let thk = if d == 2 { thickness } else { 1.0 };
    let nt = if nc < 2000 { 1 } else { femcore::par::threads_up_to(8) };
    let chunk = nc.div_ceil(nt).max(1);

    let parts: Vec<Vec<f64>> = (0..nt)
        .into_par_iter()
        .map(|tid| {
            let mut f = vec![0.0; d * n1];
            let mut vals = Vec::new();
            let mut dxi = Vec::new();
            let mut sn = Vec::new();
            let mut sdn = Vec::new();
            let mut grad = Vec::new();
            let mut b = Vec::new();
            for ci in (tid * chunk)..((tid + 1) * chunk).min(nc) {
                let dom = mixed.cell_labels[ci];
                if inactive.contains(&dom) {
                    continue;
                }
                let ct = mixed.cell_type(ci);
                let iso = IsoparametricCell::from_mesh(mixed, ci);
                let in_ie = ie_domains.contains(dom);
                let deg = if in_ie {
                    8
                } else if iso.curved {
                    2 * order + 4
                } else {
                    2 * order + 2
                };
                let rule = quad_rule_mixed(ct, deg);
                let cdofs = dof_map.cell_dofs(ci);
                let nloc = cdofs.len();
                let nd = d * nloc;
                grad.resize(nloc * d, 0.0);
                b.resize(ns * nd, 0.0);
                let mut fe = vec![0.0; nd];

                for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                    let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                    lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                    let tq: f64 = (0..nloc).map(|a| vals[a] * t_field[cdofs[a]]).sum();
                    if (tq - tref).abs() == 0.0 {
                        continue;
                    }
                    let sig = thermal_stress(&c6[ci], d, plane, &thermal_strain(&alpha[ci], tq - tref));
                    let w_eff = if in_ie {
                        let ie_j = ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                        for a in 0..nloc {
                            let gp = jg.map_grad(dxi[a]);
                            for i in 0..d {
                                let mut acc = 0.0;
                                for k in 0..d {
                                    acc += ie_j.j_inv[k][i] * gp[k];
                                }
                                grad[a * d + i] = acc;
                            }
                        }
                        wq * jg.det.abs() * ie_j.det * thk
                    } else {
                        for a in 0..nloc {
                            let gp = jg.map_grad(dxi[a]);
                            grad[a * d..a * d + d].copy_from_slice(&gp[..d]);
                        }
                        wq * jg.det.abs() * thk
                    };
                    b_matrix(&grad, nloc, d, &mut b);
                    for i in 0..nd {
                        let mut acc = 0.0;
                        for s in 0..ns {
                            acc += b[s * nd + i] * sig[s];
                        }
                        fe[i] += w_eff * acc;
                    }
                }
                for k in 0..d {
                    for (a, &gd) in cdofs.iter().enumerate() {
                        f[k * n1 + gd] += fe[k * nloc + a];
                    }
                }
            }
            f
        })
        .collect();

    let mut parts = parts.into_iter();
    let mut f = parts.next().unwrap_or_else(|| vec![0.0; d * n1]);
    for p in parts {
        f.iter_mut().zip(&p).for_each(|(a, b)| *a += b);
    }
    f
}

/// Cauchy stress (solver Voigt order, six components, also in 2D) at every local node of every
/// cell from the displacement gradients `grad[p*d*d + i*d + k] = d u_i / d x_k`. `eth` holds the
/// thermal strain vector per cell node, or is empty.
pub fn stress_at_nodes(
    grad: &[f64],
    nc: usize,
    nloc: usize,
    d: usize,
    c6: &[M6],
    plane: Plane,
    eth: &[[f64; 6]],
) -> Vec<[f64; 6]> {
    let mut out = vec![[0.0; 6]; nc * nloc];
    for c in 0..nc {
        let dm = if d == 2 && plane == Plane::Stress { Some(d_matrix(&c6[c], 2, Plane::Stress)) } else { None };
        for a in 0..nloc {
            let p = c * nloc + a;
            let g = &grad[p * d * d..(p + 1) * d * d];
            let mut e = [0.0; 6];
            for i in 0..d {
                e[i] = g[i * d + i];
            }
            e[3] = g[1] + g[d]; // gamma_xy = du/dy + dv/dx
            if d == 3 {
                e[4] = g[3 + 2] + g[3 * 2 + 1]; // gamma_yz
                e[5] = g[2] + g[3 * 2]; // gamma_xz
            }
            let th = if eth.is_empty() { [0.0; 6] } else { eth[p] };
            let mut m = [0.0; 6];
            for i in 0..6 {
                m[i] = e[i] - th[i];
            }
            match &dm {
                None => {
                    for i in 0..6 {
                        out[p][i] = (0..6).map(|j| c6[c][i][j] * m[j]).sum();
                    }
                }
                Some(dm) => {
                    // plane stress: only the in-plane components carry stress
                    for (ai, &i) in INPLANE.iter().enumerate() {
                        out[p][i] = (0..3).map(|bj| dm[ai * 3 + bj] * m[INPLANE[bj]]).sum();
                    }
                }
            }
        }
    }
    out
}

/// von Mises stress of a six-component stress.
pub fn von_mises(s: &[f64; 6]) -> f64 {
    let (sx, sy, sz, sxy, syz, sxz) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2)) + 3.0 * (sxy * sxy + syz * syz + sxz * sxz))
        .sqrt()
}

pub struct Setup {
    pub d: usize,
    pub n1: usize,
    pub c6: Vec<M6>,
    pub fixed: Vec<usize>,
    pub vals: Vec<f64>,
}

/// Setup stage (materials and essential Dirichlet constraints), supporting monomorphized [`SourceField`].
pub fn setup_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<Setup, String> {
    crate::problem::check_simplex_straight(mesh)?;
    let d = mesh.dim;
    if d != D {
        return Err(format!("SourceField<{D}> dimension does not match {d}D mesh"));
    }
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let c6: Vec<M6> = spec
        .materials
        .cell_stiffness_skipping(mesh, &spec.inactive)?
        .iter()
        .map(|c| femcore::section::plane_stiffness(c, d))
        .collect();
    let n1 = mesh.n_dofs(order);
    let mut seen = BTreeMap::new();
    for dof in untouched_dofs(mesh, order, d, &spec.inactive) {
        seen.insert(dof, 0.0);
    }
    let coords = if source.is_uniform() { Vec::new() } else { mesh.dof_coordinates(order) };
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, values } = &bc.kind {
            if source.is_uniform() {
                let bd = boundary_dofs(mesh, &bc.labels, order);
                for (k, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                    if fixed[w] {
                        for &dof in &bd {
                            seen.entry(dof + k * n1).or_insert(values[w]);
                        }
                    }
                }
            } else {
                for &lbl in &bc.labels {
                    let bd = boundary_dofs(mesh, &[lbl], order);
                    for &dof in &bd {
                        let mut x = [0.0; D];
                        x.copy_from_slice(&coords[dof * D..dof * D + D]);
                        let u_bc = source.dirichlet_re(x, 0.0, lbl as u32);
                        for (k, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                            if fixed[w] {
                                seen.entry(dof + k * n1).or_insert(u_bc[k]);
                            }
                        }
                    }
                }
            }
        }
    }
    let fixed: Vec<usize> = seen.keys().copied().collect();
    let vals: Vec<f64> = seen.values().copied().collect();
    Ok(Setup { d, n1, c6, fixed, vals })
}

pub fn setup(mesh: &Discretisation, spec: &Problem) -> Result<Setup, String> {
    match mesh.dim {
        2 => setup_with_source::<2, UniformSource>(mesh, spec, &UniformSource),
        3 => setup_with_source::<3, UniformSource>(mesh, spec, &UniformSource),
        d => Err(format!("unsupported mesh dimension {d}")),
    }
}

fn setup_mixed_with_source<const D: usize, S: SourceField<D>>(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    ie_domains: &ResolvedIeDomains,
    source: &S,
) -> Result<Setup, String> {
    let d = mixed.dim;
    if d != D {
        return Err(format!("SourceField<{D}> dimension does not match {d}D mesh"));
    }
    mixed_check_boundary_labels(mixed, &spec.bcs)?;
    let c6: Vec<M6> = mixed_cell_stiffness_skipping(&spec.materials, mixed, &spec.inactive)?
        .iter()
        .map(|c| femcore::section::plane_stiffness(c, d))
        .collect();
    let n1 = dof_map.n_dofs;
    let mut seen = BTreeMap::new();
    for dof in mixed_untouched_dofs(mixed, dof_map, d, &spec.inactive) {
        seen.insert(dof, 0.0);
    }
    let coords = if source.is_uniform() { Vec::new() } else { dof_map.dof_coordinates(mixed) };
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, values } = &bc.kind {
            if source.is_uniform() {
                let bd = dof_map.boundary_dofs(mixed, &bc.labels);
                for (k, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                    if fixed[w] {
                        for &dof in &bd {
                            seen.entry(dof + k * n1).or_insert(values[w]);
                        }
                    }
                }
            } else {
                for &lbl in &bc.labels {
                    let bd = dof_map.boundary_dofs(mixed, &[lbl]);
                    for &dof in &bd {
                        let mut x = [0.0; D];
                        x.copy_from_slice(&coords[dof * D..dof * D + D]);
                        let u_bc = source.dirichlet_re(x, 0.0, lbl as u32);
                        for (k, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                            if fixed[w] {
                                seen.entry(dof + k * n1).or_insert(u_bc[k]);
                            }
                        }
                    }
                }
            }
        }
    }
    for dof in ie_domains.outer_clamped_dofs(mixed, dof_map, spec) {
        for k in 0..d {
            seen.entry(dof + k * n1).or_insert(0.0);
        }
    }
    let fixed: Vec<usize> = seen.keys().copied().collect();
    let vals: Vec<f64> = seen.values().copied().collect();
    Ok(Setup { d, n1, c6, fixed, vals })
}

pub struct Thermal {
    pub alpha: Vec<M3>,
    pub t: Vec<f64>,
    pub tref: f64,
}

/// Stage 3 of the real solid kernels: the displacement as one vector field `u` (flat names `|u|`, `ux`, `uz`; its values are the
/// displacement per dof), `von Mises` and the stresses under wafer names, and the same components per cell without averaging.
fn fields_of(
    mesh: &Discretisation,
    u: &[f64],
    order: usize,
    st: &Setup,
    plane: Plane,
    thermal: Option<&Thermal>,
    with_stress: bool,
) -> (Vec<Field>, CellFields) {
    let d = st.d;
    let n1 = st.n1;
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();
    let names = femcore::section::component_names(d);
    let comps: Vec<Vec<f64>> = (0..d).map(|k| u[k * n1..(k + 1) * n1].to_vec()).collect();
    let mut fields = vec![Field::vector("u", "m", names, &comps)];
    let mut cells = CellFields::new(mesh, order);
    if with_stress {
        let mut tl = Instant::now();
        let grad = gradients_at_nodes(mesh, u, d, order);
        lap(&mut tl, "  gradients");
        let mut eth = vec![];
        if let Some(th) = thermal {
            let mut dofs = Vec::new();
            for c in 0..nc {
                mesh.cell_dofs(c, order, &mut dofs);
                for a in 0..nloc {
                    eth.push(thermal_strain(&th.alpha[c], th.t[dofs[a]] - th.tref));
                }
            }
        }
        let s = stress_at_nodes(&grad, nc, nloc, d, &st.c6, plane, &eth);
        lap(&mut tl, "  stress at nodes");
        let vm: Vec<f64> = s.iter().map(von_mises).collect();
        fields.push(Field::scalar("von Mises", "Pa", node_average(mesh, &vm, 1, order)));
        cells.push_real("von Mises", "Pa", vm);
        for (k, ax) in names.iter().enumerate() {
            cells.push_real(
                format!("u{ax}"),
                "m",
                femcore::fem::cell_values_at_nodes(mesh, &u[k * n1..(k + 1) * n1], order),
            );
        }
        let idxs: &[usize] = if d == 2 { &[0, 1, 3, 2] } else { &[0, 1, 2, 3, 4, 5] };
        for &idx in idxs {
            let v: Vec<f64> = s.iter().map(|x| x[idx]).collect();
            fields.push(Field::scalar(femcore::section::stress_name(d, idx), "Pa", node_average(mesh, &v, 1, order)));
            cells.push_real(femcore::section::stress_name(d, idx), "Pa", v);
        }
    } else {
        for (k, ax) in names.iter().enumerate() {
            cells.push_real(
                format!("u{ax}"),
                "m",
                femcore::fem::cell_values_at_nodes(mesh, &u[k * n1..(k + 1) * n1], order),
            );
        }
    }
    (fields, cells)
}

#[allow(clippy::too_many_arguments)]
fn fields_of_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    u: &[f64],
    u_unit: &str,
    order: usize,
    st: &Setup,
    plane: Plane,
    thermal: Option<&Thermal>,
    ie_domains: &ResolvedIeDomains,
    inactive: &[i32],
    with_stress: bool,
) -> (Vec<Field>, Vec<CellField>) {
    let d = st.d;
    let n1 = st.n1;
    let nc = mixed.n_cells();
    let space = dof_map.space;
    let names = femcore::section::component_names(d);

    let mut fields = vec![proj.project_vector_field(mixed, dof_map, "u", u_unit, names, u, order)];
    let mut cell_fields = Vec::new();

    let mut m_offsets = Vec::with_capacity(nc + 1);
    m_offsets.push(0);
    for mc in 0..nc {
        m_offsets.push(m_offsets[mc] + dof_map.cell_dofs(mc).len());
    }
    let total_cell_nodes = *m_offsets.last().unwrap();

    if with_stress {
        let mut s_all = Vec::with_capacity(total_cell_nodes);
        let mut vals = Vec::new();
        let mut dxi = Vec::new();
        let mut sn = Vec::new();
        let mut sdn = Vec::new();

        for mc in 0..nc {
            let dom = mixed.cell_labels[mc];
            let cdofs = dof_map.cell_dofs(mc);
            let nloc = cdofs.len();
            if inactive.contains(&dom) {
                s_all.extend(std::iter::repeat_n([0.0; 6], nloc));
                continue;
            }
            let ct = mixed.cell_type(mc);
            let iso = IsoparametricCell::from_mesh(mixed, mc);
            let in_ie = ie_domains.contains(dom);
            let ref_nodes = lagrange_ref_nodes(ct, space);
            let mut grad_cell = vec![0.0; nloc * d * d];
            let mut eth_cell = if thermal.is_some() { Vec::with_capacity(nloc) } else { Vec::new() };

            for (a, xi) in ref_nodes.iter().enumerate() {
                let jg = iso.eval_into(xi, &mut sn, &mut sdn);
                lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
                let mut g_phys = [[0.0; 3]; 3];
                for (b_idx, &gd) in cdofs.iter().enumerate() {
                    let gp = jg.map_grad(dxi[b_idx]);
                    for i in 0..d {
                        let ui = u[i * n1 + gd];
                        for k in 0..d {
                            g_phys[i][k] += ui * gp[k];
                        }
                    }
                }
                if in_ie {
                    let ie_j = ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, d);
                    for i in 0..d {
                        for m in 0..d {
                            let mut acc = 0.0;
                            for k in 0..d {
                                acc += g_phys[i][k] * ie_j.j_inv[k][m];
                            }
                            grad_cell[a * d * d + i * d + m] = acc;
                        }
                    }
                } else {
                    for i in 0..d {
                        for m in 0..d {
                            grad_cell[a * d * d + i * d + m] = g_phys[i][m];
                        }
                    }
                }
                if let Some(th) = thermal {
                    eth_cell.push(thermal_strain(&th.alpha[mc], th.t[cdofs[a]] - th.tref));
                }
            }
            let cell_stress = stress_at_nodes(&grad_cell, 1, nloc, d, &st.c6[mc..mc + 1], plane, &eth_cell);
            s_all.extend_from_slice(&cell_stress);
        }

        let vm_mixed: Vec<f64> = s_all.iter().map(von_mises).collect();
        let vm_cf = proj.project_cell_field(mixed, space, "von Mises", "Pa", &vm_mixed, &[], &m_offsets, order);
        fields.push(Field::scalar("von Mises", "Pa", node_average(simplex, &vm_cf.re, 1, order)));
        cell_fields.push(vm_cf);

        for (k, ax) in names.iter().enumerate() {
            let mut uk_mixed = Vec::with_capacity(total_cell_nodes);
            for mc in 0..nc {
                for &gd in dof_map.cell_dofs(mc) {
                    uk_mixed.push(u[k * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                space,
                format!("u{ax}"),
                u_unit,
                &uk_mixed,
                &[],
                &m_offsets,
                order,
            ));
        }

        let idxs: &[usize] = if d == 2 { &[0, 1, 3, 2] } else { &[0, 1, 2, 3, 4, 5] };
        for &idx in idxs {
            let sname = femcore::section::stress_name(d, idx);
            let comp_mixed: Vec<f64> = s_all.iter().map(|x| x[idx]).collect();
            let s_cf = proj.project_cell_field(mixed, space, sname, "Pa", &comp_mixed, &[], &m_offsets, order);
            fields.push(Field::scalar(sname, "Pa", node_average(simplex, &s_cf.re, 1, order)));
            cell_fields.push(s_cf);
        }
    } else {
        for (k, ax) in names.iter().enumerate() {
            let mut uk_mixed = Vec::with_capacity(total_cell_nodes);
            for mc in 0..nc {
                for &gd in dof_map.cell_dofs(mc) {
                    uk_mixed.push(u[k * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                space,
                format!("u{ax}"),
                u_unit,
                &uk_mixed,
                &[],
                &m_offsets,
                order,
            ));
        }
    }
    (fields, cell_fields)
}

/// Stage 3 of the complex-frequency solid kernels (PML, Floquet, low-reflecting): the mode shape `u` over all dofs
/// (component blocks) with its largest component made real and of unit modulus, as the plot fields `|u|` and `Re u{x,z}` (scalars: the
/// vector `Re u` would add a norm `|Re u|` that the pre-restructure results do not have), and the complex components per cell.
pub(crate) fn complex_mode_fields(mesh: &Discretisation, order: usize, u: &[C]) -> (Vec<Field>, CellFields) {
    let (d, n1) = (mesh.dim, mesh.n_dofs(order));
    let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let scale = 1.0 / big.norm().max(1e-300);
    let mag: Vec<f64> = (0..n1).map(|a| (0..d).map(|k| u[k * n1 + a].norm_sqr()).sum::<f64>().sqrt() * scale).collect();
    let mut fields = vec![Field::scalar("|u|", "a.u.", mag)];
    let mut cells = CellFields::new(mesh, order);
    let nloc = mesh.nloc(order);
    for (comp, ax) in femcore::section::component_names(d).iter().enumerate() {
        let v: Vec<f64> = (0..n1).map(|a| (u[comp * n1 + a] * rot * scale).re).collect();
        fields.push(Field::scalar(format!("Re u{ax}"), "a.u.", v));
        let (mut re_c, mut im_c) =
            (Vec::with_capacity(mesh.n_cells() * nloc), Vec::with_capacity(mesh.n_cells() * nloc));
        let mut dofs = Vec::new();
        for cell in 0..mesh.n_cells() {
            mesh.cell_dofs(cell, order, &mut dofs);
            for &dof in dofs.iter().take(nloc) {
                let z = u[comp * n1 + dof] * rot * scale;
                re_c.push(z.re);
                im_c.push(z.im);
            }
        }
        cells.push_complex(format!("u{ax}"), "a.u.", re_c, im_c);
    }
    (fields, cells)
}

/// Stage-3 complex-frequency mode fields projected from [`MixedDiscretisation`] onto [`Discretisation`] (`SimplexMesh`).
pub(crate) fn complex_mode_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    proj: &MixedToSimplexProjector,
    order: usize,
    u: &[C],
) -> (Vec<Field>, Vec<CellField>) {
    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let scale = 1.0 / big.norm().max(1e-300);

    let mut re_blocks = Vec::with_capacity(d * n1);
    let mut im_blocks = Vec::with_capacity(d * n1);
    for comp in 0..d {
        for a in 0..n1 {
            let z = u[comp * n1 + a] * rot * scale;
            re_blocks.push(z.re);
            im_blocks.push(z.im);
        }
    }
    let re_proj = proj.project_lagrange_blocks(mixed, dof_map, &re_blocks, d, order);
    let im_proj = proj.project_lagrange_blocks(mixed, dof_map, &im_blocks, d, order);
    let n_out = re_proj.len() / d;
    let mag: Vec<f64> = (0..n_out)
        .map(|a| {
            (0..d)
                .map(|k| {
                    let r0 = re_proj[k * n_out + a];
                    let i0 = im_proj[k * n_out + a];
                    r0 * r0 + i0 * i0
                })
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    let mut fields = vec![Field::scalar("|u|", "a.u.", mag)];

    let mut m_offsets = Vec::with_capacity(mixed.n_cells() + 1);
    m_offsets.push(0);
    for mc in 0..mixed.n_cells() {
        m_offsets.push(m_offsets[mc] + dof_map.cell_dofs(mc).len());
    }
    let total_nodes = *m_offsets.last().unwrap();
    let mut cell_fields = Vec::with_capacity(d);

    for (comp, ax) in femcore::section::component_names(d).iter().enumerate() {
        fields.push(Field::scalar(format!("Re u{ax}"), "a.u.", re_proj[comp * n_out..(comp + 1) * n_out].to_vec()));
        let mut re_c = Vec::with_capacity(total_nodes);
        let mut im_c = Vec::with_capacity(total_nodes);
        for mc in 0..mixed.n_cells() {
            for &gd in dof_map.cell_dofs(mc) {
                re_c.push(re_blocks[comp * n1 + gd]);
                im_c.push(im_blocks[comp * n1 + gd]);
            }
        }
        cell_fields.push(proj.project_cell_field(
            mixed,
            dof_map.space,
            format!("u{ax}"),
            "a.u.",
            &re_c,
            &im_c,
            &m_offsets,
            order,
        ));
    }
    (fields, cell_fields)
}

/// Stage 2 of the real solid eigenfrequency kernel: the `nm` eigenpairs of `K u = w^2 M u` on the free dofs nearest the
/// search frequency (shift-invert Lanczos; for the lowest modes of a large problem LOBPCG with a multigrid
/// preconditioner, with a factorisation as the fallback), and a note on which route was taken.
#[allow(clippy::too_many_arguments)]
pub fn solve_real_modes(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    k: &Csr,
    m: &Csr,
    free: &[usize],
    nm: usize,
    f_search: f64,
) -> Result<(femcore::linalg::driver::EigPairs, Option<String>), String> {
    let (d, order) = (st.d, spec.order);
    let kff = k.submatrix(free, free);
    let mff = m.submatrix(free, free);
    let mut eig_note: Option<String> = None;
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    // for the lowest modes (shift below the spectrum) K - sigma M is positive definite and large
    // 3D problems are solved with multigrid-preconditioned CG instead of a sparse factorisation
    let big = if mesh.dim == 3 {
        free.len() > femcore::fem::LOBPCG_MIN_FREE
    } else {
        free.len() > femcore::fem::AMG_THRESHOLD
    };
    let want_amg = match spec.solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => sigma <= 0.0,
        LinearSolver::Auto => sigma <= 0.0 && big,
    };
    let pairs = if want_amg {
        // no factorisation: block preconditioned CG (LOBPCG) with one multigrid cycle as preconditioner
        let shifted_full = k.add_scaled(m, -sigma);
        let sys = femcore::fem::AmgSystem::new(mesh, order, d, &shifted_full, &st.fixed)?;
        let factor = sys.convergence_factor();
        // multigrid that converges slowly (very thin or stretched parts) loses to a factorisation
        // unless the factorisation would not fit
        let use_amg = spec.solver == LinearSolver::Iterative
            || factor <= femcore::fem::LOBPCG_MAX_FACTOR
            || free.len() > femcore::fem::FACTORISE_MAX_FREE;
        if use_amg {
            // in Auto mode a stalled iteration gives up early and the factorisation takes over
            let guard = spec.solver != LinearSolver::Iterative;
            let res = femcore::linalg::driver::lobpcg_guarded(
                &kff,
                &mff,
                nm,
                nm.clamp(4, 8),
                &|r: &[f64]| sys.precondition_free(r),
                1e-8,
                400,
                guard,
            );
            match res {
                Ok(r) if r.converged || spec.solver == LinearSolver::Iterative => {
                    eig_note = Some(format!(
                        "LOBPCG with AMG preconditioner: {} levels {:?}, convergence factor {:.2}, setup {:.1} s{}",
                        sys.n_levels(),
                        sys.sizes(),
                        factor,
                        sys.setup_seconds,
                        if r.converged { "" } else { " (not fully converged)" }
                    ));
                    r
                }
                Err(e) if e.starts_with(femcore::linalg::driver::STALLED) => {
                    eig_note = Some(format!("{e}; used a sparse factorisation instead"));
                    eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
                }
                _ => eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?,
            }
        } else {
            eig_note =
                Some(format!("multigrid converges slowly here (factor {factor:.2}), used a sparse factorisation"));
            eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
        }
    } else {
        eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?
    };
    Ok((pairs, eig_note))
}

/// Stage 1 of the static solid kernel: the load vector (boundary tractions in wafer axes, thermal expansion) next to the
/// stiffness matrix. A static problem needs a fixed boundary.
pub struct StaticAssembly {
    pub k: Csr,
    pub f: Vec<f64>,
    pub thermal: Option<Thermal>,
}

fn simplex_outward_normal<const D: usize>(mesh: &Discretisation, fv: &[u32], ci: usize) -> [f64; D] {
    let pt = |v: u32| -> [f64; 3] {
        let p = mesh.point(v as usize);
        [p[0], p[1], if D == 3 { p[2] } else { 0.0 }]
    };
    let n3 = if D == 2 {
        let (a, b) = (pt(fv[0]), pt(fv[1]));
        [b[1] - a[1], -(b[0] - a[0]), 0.0]
    } else {
        let (a, b, c) = (pt(fv[0]), pt(fv[1]), pt(fv[2]));
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    };
    let len = (n3[0] * n3[0] + n3[1] * n3[1] + n3[2] * n3[2]).sqrt().max(1e-300);
    let mut n_out = [0.0; D];
    for k in 0..D {
        n_out[k] = n3[k] / len;
    }
    let cv = mesh.cell(ci);
    if let Some(&v_opp) = cv.iter().find(|v| !fv.contains(v)) {
        let opp = pt(v_opp);
        let mut fc = [0.0; D];
        for &v in fv {
            let pv = mesh.point(v as usize);
            for k in 0..D {
                fc[k] += pv[k] / fv.len() as f64;
            }
        }
        let dot: f64 = (0..D).map(|k| (fc[k] - opp[k]) * n_out[k]).sum();
        if dot < 0.0 {
            for nk in &mut n_out {
                *nk = -*nk;
            }
        }
    }
    n_out
}

fn add_simplex_mms_loads<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    thk_scale: f64,
    source: &S,
    f: &mut [f64],
) {
    let order = spec.order;
    let n1 = mesh.n_dofs(order);
    let geo = geometry(mesh);

    // Volumetric MMS body force: degree-6 simplex quadrature handles high-order smooth MMS fields.
    let (ct, ref_vol_inv) = if D == 2 { (CellType::Tri, 2.0) } else { (CellType::Tet, 6.0) };
    let vrule = quad_rule_mixed(ct, 6);
    let mut dofs = Vec::new();
    for ci in 0..mesh.n_cells() {
        let dom = mesh.cell_labels[ci];
        if spec.inactive.contains(&dom) {
            continue;
        }
        let cv = mesh.cell(ci);
        mesh.cell_dofs(ci, order, &mut dofs);
        let base_w = geo.vol[ci] * ref_vol_inv * thk_scale;
        for (xi, &wq) in vrule.pts.iter().zip(&vrule.w) {
            let lam = if D == 2 {
                vec![1.0 - xi[0] - xi[1], xi[0], xi[1]]
            } else {
                vec![1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]]
            };
            let mut x_q = [0.0; D];
            for (i, &v) in cv.iter().enumerate() {
                let pv = mesh.point(v as usize);
                for k in 0..D {
                    x_q[k] += lam[i] * pv[k];
                }
            }
            let (nv, _) = shape(&lam, order);
            let rhs = source.volume_rhs_re(x_q, 0.0, dom as u32);
            let w_phys = wq * base_w;
            for (a, &gd) in dofs.iter().enumerate() {
                let wa = w_phys * nv[a];
                for comp in 0..D {
                    f[comp * n1 + gd] += wa * rhs[comp];
                }
            }
        }
    }

    // Spatially varying Neumann traction on `BcKind::Load` boundaries.
    let has_load = spec.bcs.iter().any(|b| matches!(b.kind, BcKind::Load(_)));
    if !has_load {
        return;
    }
    let fkey = |v: &[u32]| -> [u32; 3] {
        let mut k = [u32::MAX; 3];
        k[..v.len()].copy_from_slice(v);
        k[..v.len()].sort_unstable();
        k
    };
    let mut cell_of: HashMap<[u32; 3], usize> = HashMap::new();
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for skip in 0..=D {
            let fv_buf: Vec<u32> = (0..=D).filter(|&i| i != skip).map(|i| v[i]).collect();
            cell_of.entry(fkey(&fv_buf)).or_insert(ci);
        }
    }
    let (ft, ref_fac_inv) = if D == 2 { (FacetType::Edge2, 0.5) } else { (FacetType::Tri3, 2.0) };
    let frule = facet_rule_mixed(ft, 6);
    for bc in &spec.bcs {
        let BcKind::Load(t) = &bc.kind else { continue };
        for fc in facets_with(mesh, &bc.labels) {
            let fv = mesh.facet(fc);
            let lbl = mesh.facet_labels[fc] as u32;
            let ci = cell_of[&fkey(fv)];
            let n_out = simplex_outward_normal::<D>(mesh, fv, ci);
            let area = facet_measure(mesh, fc);
            mesh.facet_dofs(fc, order, &mut dofs);
            for (xi, &wq) in frule.pts.iter().zip(&frule.w) {
                let lam = if D == 2 {
                    vec![0.5 * (1.0 - xi[0]), 0.5 * (1.0 + xi[0])]
                } else {
                    vec![1.0 - xi[0] - xi[1], xi[0], xi[1]]
                };
                let mut x_q = [0.0; D];
                for (i, &v) in fv.iter().enumerate() {
                    let pv = mesh.point(v as usize);
                    for k in 0..D {
                        x_q[k] += lam[i] * pv[k];
                    }
                }
                let (nv, _) = shape(&lam, order);
                let tr = source.neumann_robin_re(x_q, n_out, 0.0, lbl);
                let w_fac = wq * ref_fac_inv * area * thk_scale;
                for (a, &gd) in dofs.iter().enumerate() {
                    let wa = w_fac * nv[a];
                    for (comp, &w) in femcore::section::wafer_axes(D).iter().enumerate() {
                        f[comp * n1 + gd] += wa * (t[w] + tr[comp]);
                    }
                }
            }
        }
    }
}

/// Static assembly supporting monomorphized [`SourceField`].
pub fn assemble_static_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    k: Csr,
    source: &S,
) -> Result<StaticAssembly, String> {
    crate::problem::check_simplex_straight(mesh)?;
    if st.d != D || mesh.dim != D {
        return Err(format!("SourceField<{D}> dimension mismatch (mesh is {}D)", mesh.dim));
    }
    if st.fixed.is_empty() {
        return Err("a static solid problem needs a fixed boundary (otherwise it floats freely)".into());
    }
    let (d, n1, order, thk) = (st.d, st.n1, spec.order, spec.thickness);
    let n = d * n1;
    let mut f = vec![0.0; n];
    if source.is_uniform() {
        for bc in &spec.bcs {
            if let BcKind::Load(t) = &bc.kind {
                for (comp, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                    if t[w] != 0.0 {
                        let tc = t[w] * if d == 2 { thk } else { 1.0 };
                        add_boundary_load(mesh, &bc.labels, &|_| tc, order, &mut f, d, comp);
                    }
                }
            }
        }
    } else {
        let thk_scale = if d == 2 { thk } else { 1.0 };
        add_simplex_mms_loads::<D, S>(mesh, spec, thk_scale, source, &mut f);
    }
    let mut thermal: Option<Thermal> = None;
    if let Some(th) = &spec.thermal {
        let alpha: Vec<M3> = spec
            .materials
            .cell_tensor(mesh, TensorKind::Expansion, Some(0.0))?
            .iter()
            .map(|a| femcore::section::plane_tensor(a, d))
            .collect();
        let t = match &th.source {
            ThermalSource::Uniform(tv) => vec![*tv; n1],
            ThermalSource::Imported(t) if t.len() == n1 => t.clone(),
            ThermalSource::Imported(_) => {
                return Err("no heat-transfer result on this mesh and element order; solve Heat transfer first".into())
            }
        };
        let fth = thermal_load(mesh, &st.c6, &alpha, &t, th.t_ref, order, spec.plane, thk);
        f.iter_mut().zip(&fth).for_each(|(a, b)| *a += b);
        thermal = Some(Thermal { alpha, t, tref: th.t_ref });
    }
    Ok(StaticAssembly { k, f, thermal })
}

pub fn assemble_static(mesh: &Discretisation, spec: &Problem, st: &Setup, k: Csr) -> Result<StaticAssembly, String> {
    match mesh.dim {
        2 => assemble_static_with_source::<2, UniformSource>(mesh, spec, st, k, &UniformSource),
        3 => assemble_static_with_source::<3, UniformSource>(mesh, spec, st, k, &UniformSource),
        d => Err(format!("unsupported mesh dimension {d}")),
    }
}

pub(crate) fn add_boundary_load_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    labels: &[i32],
    tc: f64,
    out: &mut [f64],
    comp: usize,
) {
    let n1 = dof_map.n_dofs;
    let deg = 2 * dof_map.space.order() + 2;
    let mut dofs = Vec::new();
    for f in 0..mixed.n_facets() {
        if !labels.contains(&mixed.facet_labels[f]) {
            continue;
        }
        let ft = mixed.facet_type(f);
        dof_map.facet_dofs(mixed, f, &mut dofs);
        let rule = facet_rule_mixed(ft, deg);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (_x, ds, _n) = eval_facet_jacobian(mixed, f, xi);
            if ds <= 0.0 {
                continue;
            }
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            let w_ds = tc * wq * ds;
            for (a, &gd) in dofs.iter().enumerate() {
                out[comp * n1 + gd] += w_ds * nv[a];
            }
        }
    }
}

fn add_mixed_mms_loads<const D: usize, S: SourceField<D>>(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    ie_domains: &ResolvedIeDomains,
    thk_scale: f64,
    source: &S,
    f: &mut [f64],
) {
    let n1 = dof_map.n_dofs;
    let space = dof_map.space;
    let (mut vals, mut dxi, mut sn, mut sdn) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    for ci in 0..mixed.n_cells() {
        let dom = mixed.cell_labels[ci];
        if spec.inactive.contains(&dom) {
            continue;
        }
        let ct = mixed.cell_type(ci);
        let iso = IsoparametricCell::from_mesh(mixed, ci);
        let in_ie = ie_domains.contains(dom);
        let rule = quad_rule_mixed(ct, if in_ie { 10 } else { 6 });
        let cdofs = dof_map.cell_dofs(ci);
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let jg = iso.eval_into(xi, &mut sn, &mut sdn);
            lagrange_shape_into(ct, space, xi, &mut vals, &mut dxi);
            let ie_det = if in_ie { ie_domains.eval_cell_jacobian(dom, &iso, &jg, &sn, &sdn, D).det } else { 1.0 };
            let mut x_q = [0.0; D];
            x_q.copy_from_slice(&jg.x[..D]);
            let rhs = source.volume_rhs_re(x_q, 0.0, dom as u32);
            let w_eff = wq * jg.det.abs() * ie_det * thk_scale;
            for (a, &gd) in cdofs.iter().enumerate() {
                let wa = w_eff * vals[a];
                for comp in 0..D {
                    f[comp * n1 + gd] += wa * rhs[comp];
                }
            }
        }
    }

    let has_load = spec.bcs.iter().any(|b| matches!(b.kind, BcKind::Load(_)));
    if !has_load {
        return;
    }
    let mut edge_cell: HashMap<(u32, u32), usize> = HashMap::new();
    let mut tri_cell: HashMap<[u32; 3], usize> = HashMap::new();
    let mut quad_cell: HashMap<[u32; 4], usize> = HashMap::new();
    for ci in 0..mixed.n_cells() {
        if spec.inactive.contains(&mixed.cell_labels[ci]) {
            continue;
        }
        let vs = mixed.cell(ci);
        let ct = mixed.cell_type(ci);
        if D == 2 {
            for &(i, j) in cell_local_edges(ct) {
                let (a, b) = (vs[i], vs[j]);
                edge_cell.entry((a.min(b), a.max(b))).or_insert(ci);
            }
        } else {
            for lf in cell_local_tri_faces(ct) {
                tri_cell.entry(sorted_tri_face([vs[lf[0]], vs[lf[1]], vs[lf[2]]])).or_insert(ci);
            }
            for lq in cell_local_quad_faces(ct) {
                quad_cell.entry(canonical_quad_face([vs[lq[0]], vs[lq[1]], vs[lq[2]], vs[lq[3]]])).or_insert(ci);
            }
        }
    }
    let cell_centroid = |ci: usize| -> [f64; D] {
        let vs = mixed.cell(ci);
        let mut c0 = [0.0; D];
        for &v in vs {
            let pt = mixed.point(v as usize);
            for k in 0..D {
                c0[k] += pt[k] / vs.len() as f64;
            }
        }
        c0
    };
    let mut dofs = Vec::new();
    for bc in &spec.bcs {
        let BcKind::Load(t) = &bc.kind else { continue };
        for fc in 0..mixed.n_facets() {
            if !bc.labels.contains(&mixed.facet_labels[fc]) {
                continue;
            }
            let lbl = mixed.facet_labels[fc] as u32;
            let ft = mixed.facet_type(fc);
            let fv = mixed.facet(fc);
            let ci_opt = match ft {
                FacetType::Edge2 => edge_cell.get(&(fv[0].min(fv[1]), fv[0].max(fv[1]))).copied(),
                FacetType::Tri3 => tri_cell.get(&sorted_tri_face([fv[0], fv[1], fv[2]])).copied(),
                FacetType::Quad4 => quad_cell.get(&canonical_quad_face([fv[0], fv[1], fv[2], fv[3]])).copied(),
            };
            let xc = ci_opt.map(cell_centroid);
            dof_map.facet_dofs(mixed, fc, &mut dofs);
            let rule = facet_rule_mixed(ft, 6);
            for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
                let (x3, ds, n3) = eval_facet_jacobian(mixed, fc, xi);
                if ds <= 0.0 {
                    continue;
                }
                let mut x_q = [0.0; D];
                x_q.copy_from_slice(&x3[..D]);
                let mut n_out = [0.0; D];
                n_out.copy_from_slice(&n3[..D]);
                if let Some(c0) = xc {
                    let dot: f64 = (0..D).map(|k| (x_q[k] - c0[k]) * n_out[k]).sum();
                    if dot < 0.0 {
                        for nk in &mut n_out {
                            *nk = -*nk;
                        }
                    }
                }
                let (nv, _) = facet_lagrange_shape(ft, space, xi);
                let tr = source.neumann_robin_re(x_q, n_out, 0.0, lbl);
                let w_fac = wq * ds * thk_scale;
                for (a, &gd) in dofs.iter().enumerate() {
                    let wa = w_fac * nv[a];
                    for (comp, &w) in femcore::section::wafer_axes(D).iter().enumerate() {
                        f[comp * n1 + gd] += wa * (t[w] + tr[comp]);
                    }
                }
            }
        }
    }
}

fn assemble_static_mixed_with_source<const D: usize, S: SourceField<D>>(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    spec: &Problem,
    st: &Setup,
    ie_domains: &ResolvedIeDomains,
    k: Csr,
    source: &S,
) -> Result<StaticAssembly, String> {
    if st.d != D || mixed.dim != D {
        return Err(format!("SourceField<{D}> dimension mismatch (mesh is {}D)", mixed.dim));
    }
    if st.fixed.is_empty() && !ie_domains.has_exponential_robin() {
        return Err("a static solid problem needs a fixed boundary (otherwise it floats freely)".into());
    }
    let (d, n1, order, thk) = (st.d, st.n1, spec.order, spec.thickness);
    let n = d * n1;
    let mut f = vec![0.0; n];
    if source.is_uniform() {
        for bc in &spec.bcs {
            if let BcKind::Load(t) = &bc.kind {
                for (comp, &w) in femcore::section::wafer_axes(d).iter().enumerate() {
                    if t[w] != 0.0 {
                        let tc = t[w] * if d == 2 { thk } else { 1.0 };
                        add_boundary_load_mixed(mixed, dof_map, &bc.labels, tc, &mut f, comp);
                    }
                }
            }
        }
    } else {
        let thk_scale = if d == 2 { thk } else { 1.0 };
        add_mixed_mms_loads::<D, S>(mixed, dof_map, spec, ie_domains, thk_scale, source, &mut f);
    }
    let mut thermal: Option<Thermal> = None;
    if let Some(th) = &spec.thermal {
        let alpha: Vec<M3> = mixed_cell_tensor(&spec.materials, mixed, TensorKind::Expansion, Some(0.0))?
            .iter()
            .map(|a| femcore::section::plane_tensor(a, d))
            .collect();
        let t = match &th.source {
            ThermalSource::Uniform(tv) => vec![*tv; n1],
            ThermalSource::Imported(t) if t.len() == n1 => t.clone(),
            ThermalSource::Imported(t) if t.len() == simplex.n_dofs(order) => {
                let nv = dof_map.n_vertices;
                let mut mapped = Vec::with_capacity(n1);
                mapped.extend_from_slice(&t[..nv]);
                if dof_map.space != LagrangeSpace::Order1 {
                    for &(a, b) in &mixed.edges {
                        let seid = simplex.edge_index(a, b);
                        mapped.push(t[nv + seid]);
                    }
                }
                mapped
            }
            ThermalSource::Imported(_) => {
                return Err("no heat-transfer result on this mesh and element order; solve Heat transfer first".into())
            }
        };
        let fth = thermal_load_mixed(
            mixed,
            dof_map,
            &st.c6,
            &alpha,
            &t,
            th.t_ref,
            spec.plane,
            thk,
            ie_domains,
            &spec.inactive,
        );
        f.iter_mut().zip(&fth).for_each(|(a, b)| *a += b);
        thermal = Some(Thermal { alpha, t, tref: th.t_ref });
    }
    Ok(StaticAssembly { k, f, thermal })
}

/// Stage 2: the displacement (component blocks) by the linear solver `spec.solver` names, and its note.
pub fn solve_static(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    asm: &StaticAssembly,
) -> Result<(Vec<f64>, Option<String>), String> {
    let (u, note, _) = solve_static_checked(mesh, spec, st, asm)?;
    Ok((u, note))
}

pub(crate) fn solve_static_checked(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    asm: &StaticAssembly,
) -> Result<(Vec<f64>, Option<String>, crate::diag::SmsDiag), String> {
    use femcore::linalg::driver::{backward_residual_many, DirectSolver};
    let n = asm.k.nrows();
    let free = free_dofs(n, &st.fixed);
    let n_free = free.len();
    let iterative = match spec.solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => true,
        LinearSolver::Auto => femcore::fem::amg_wanted(mesh.dim, n_free),
    } && !st.fixed.is_empty();

    let solve_direct = |note_msg: Option<String>,
                        warn_fallback: bool|
     -> Result<(Vec<f64>, Option<String>, crate::diag::SmsDiag), String> {
        let mut u = vec![0.0; n];
        for (&d, &v) in st.fixed.iter().zip(&st.vals) {
            u[d] = v;
        }
        let mut diag = crate::diag::SmsDiag::new("direct_chol");
        if n_free > 0 {
            let kff = asm.k.submatrix(&free, &free);
            let rhs: Vec<f64> = if st.fixed.is_empty() || st.vals.iter().all(|&v| v == 0.0) {
                free.iter().map(|&i| asm.f[i]).collect()
            } else {
                let kfx = asm.k.submatrix(&free, &st.fixed);
                let kx = kfx.matvec(&st.vals);
                free.iter().enumerate().map(|(a, &i)| asm.f[i] - kx[a]).collect()
            };
            let solver = DirectSolver::new(&kff, true)?;
            let (uf, chk) = solver.solve_checked(&kff, &rhs);
            diag.absorb_direct(&chk);
            for (a, &i) in free.iter().enumerate() {
                u[i] = uf[a];
            }
        }
        if warn_fallback {
            diag.push_warning(
                "amg_fallback_to_direct",
                nyx_physics::DiagnosticSeverity::Warning,
                note_msg.as_deref().unwrap_or("AMG-CG fell back to direct Cholesky"),
            );
        }
        Ok((u, note_msg, diag))
    };

    if !iterative {
        return solve_direct(None, false);
    }

    let mut uvec = vec![0.0; n];
    for (&i, &v) in st.fixed.iter().zip(&st.vals) {
        uvec[i] = v;
    }
    let ku = asm.k.matvec(&uvec);
    let rhs: Vec<f64> = (0..n).map(|i| asm.f[i] - ku[i]).collect();
    let opts = femcore::linalg::amg_elastic::AmgOptions { line_smooth: true, ..Default::default() };
    let sys = AmgSystem::new_with(mesh, spec.order, st.d, &asm.k, &st.fixed, opts)?;
    let fv: Vec<(usize, f64)> = st.fixed.iter().copied().zip(st.vals.iter().copied()).collect();
    let tol = 1e-10;
    let t_solve = Instant::now();
    let (mut u, its, rel) = sys.solve_with(None, &rhs, &fv, tol, 600);
    let solve_s = t_solve.elapsed().as_secs_f64();
    femcore::progress::check()?;
    for (&i, &v) in st.fixed.iter().zip(&st.vals) {
        u[i] = v;
    }
    let kff = asm.k.submatrix(&free, &free);
    let rhs_free: Vec<f64> = free.iter().map(|&i| rhs[i]).collect();
    let u_free: Vec<f64> = free.iter().map(|&i| u[i]).collect();
    let bw_rel = if n_free > 0 { backward_residual_many(&kff, &[rhs_free], &[u_free]) } else { 0.0 };
    let worst_rel = rel.max(bw_rel);
    if !rel.is_finite() || rel > (tol * 1e3).max(1e-7) {
        if spec.solver == LinearSolver::Auto {
            return solve_direct(
                Some(format!(
                    "multigrid did not converge ({its} iterations, residual {rel:.1e}); used the direct solver"
                )),
                true,
            );
        }
        return Err(format!(
            "the iterative solver did not converge ({its} iterations, residual {rel:.2e}); choose the direct solver"
        ));
    }
    let mut diag = crate::diag::SmsDiag::new("amg_pcg");
    diag.iterations = Some(its);
    diag.max_algebraic_residual = worst_rel;
    let note = format!(
        "AMG-CG: {} levels {:?}, {} iterations, residual {:.1e} (setup {:.1} s, solve {:.1} s)",
        sys.n_levels(),
        sys.sizes(),
        its,
        rel,
        sys.setup_seconds,
        solve_s
    );
    Ok((u, Some(note), diag))
}

fn near_nullspace_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    ncomp: usize,
    is_fixed: &[bool],
    perm: &[usize],
) -> Vec<Vec<f64>> {
    let dim = mixed.dim;
    let n = is_fixed.len();
    let n1 = dof_map.n_dofs;
    let coords = dof_map.dof_coordinates(mixed);
    let mut centre = vec![0.0; dim];
    for p in coords.chunks(dim) {
        for a in 0..dim {
            centre[a] += p[a] / n1 as f64;
        }
    }
    let mut ns: Vec<Vec<f64>> = Vec::new();
    for t in 0..dim {
        let mut v = vec![0.0; n];
        for node in 0..n1 {
            v[node * ncomp + t] = 1.0;
        }
        ns.push(v);
    }
    let nrot = if dim == 2 { 1 } else { 3 };
    for r in 0..nrot {
        let mut v = vec![0.0; n];
        for node in 0..n1 {
            let p: Vec<f64> = (0..dim).map(|a| coords[node * dim + a] - centre[a]).collect();
            let (x, y, z) = (p[0], p[1], if dim == 3 { p[2] } else { 0.0 });
            let rot: [f64; 3] = if dim == 2 {
                [-y, x, 0.0]
            } else {
                match r {
                    0 => [0.0, -z, y],
                    1 => [z, 0.0, -x],
                    _ => [-y, x, 0.0],
                }
            };
            for c in 0..dim {
                v[node * ncomp + c] = rot[c];
            }
        }
        ns.push(v);
    }
    for v in ns.iter_mut() {
        for (r, x) in v.iter_mut().enumerate() {
            if is_fixed[perm[r]] {
                *x = 0.0;
            }
        }
    }
    ns
}

fn solve_static_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    st: &Setup,
    ie_domains: &ResolvedIeDomains,
    asm: &StaticAssembly,
) -> Result<(Vec<f64>, Option<String>, crate::diag::SmsDiag), String> {
    use femcore::linalg::driver::{backward_residual_many, DirectSolver};
    let n = asm.k.nrows();
    let free = free_dofs(n, &st.fixed);
    let n_free = free.len();
    let iterative = match spec.solver {
        LinearSolver::Direct => false,
        LinearSolver::Iterative => true,
        LinearSolver::Auto => femcore::fem::amg_wanted(mixed.dim, n_free) && !ie_domains.has_rational(),
    };
    let solve_direct = |note_msg: Option<String>,
                        warn_fallback: bool|
     -> Result<(Vec<f64>, Option<String>, crate::diag::SmsDiag), String> {
        let mut u = vec![0.0; n];
        for (&d, &v) in st.fixed.iter().zip(&st.vals) {
            u[d] = v;
        }
        let mut diag = crate::diag::SmsDiag::new("direct_chol");
        if n_free > 0 {
            let kff = asm.k.submatrix(&free, &free);
            let rhs: Vec<f64> = if st.fixed.is_empty() || st.vals.iter().all(|&v| v == 0.0) {
                free.iter().map(|&i| asm.f[i]).collect()
            } else {
                let kfx = asm.k.submatrix(&free, &st.fixed);
                let kx = kfx.matvec(&st.vals);
                free.iter().enumerate().map(|(a, &i)| asm.f[i] - kx[a]).collect()
            };
            let solver = DirectSolver::new(&kff, true)?;
            let (uf, chk) = solver.solve_checked(&kff, &rhs);
            diag.absorb_direct(&chk);
            for (a, &i) in free.iter().enumerate() {
                u[i] = uf[a];
            }
        }
        if warn_fallback {
            diag.push_warning(
                "amg_fallback_to_direct",
                nyx_physics::DiagnosticSeverity::Warning,
                note_msg.as_deref().unwrap_or("AMG-CG fell back to direct Cholesky"),
            );
        }
        Ok((u, note_msg, diag))
    };
    if !iterative || dof_map.space == LagrangeSpace::Order2Tensor {
        return solve_direct(None, false);
    }

    let t_setup = Instant::now();
    let n1 = dof_map.n_dofs;
    let ncomp = st.d;
    let mut is_fixed = vec![false; n];
    for &i in &st.fixed {
        is_fixed[i] = true;
    }
    let diag = asm.k.diagonal();
    let dscale = {
        let (mut s, mut cnt) = (0.0, 0usize);
        for i in 0..n {
            if !is_fixed[i] {
                s += diag[i].abs();
                cnt += 1;
            }
        }
        if cnt > 0 && s > 0.0 {
            s / cnt as f64
        } else {
            1.0
        }
    };

    let mut uvec = vec![0.0; n];
    for (&i, &v) in st.fixed.iter().zip(&st.vals) {
        uvec[i] = v;
    }
    let ku = asm.k.matvec(&uvec);
    let mut rhs: Vec<f64> = (0..n).map(|i| if is_fixed[i] { 0.0 } else { asm.f[i] - ku[i] }).collect();
    for (&i, &v) in st.fixed.iter().zip(&st.vals) {
        rhs[i] = dscale * v;
    }

    let perm: Vec<usize> = (0..n).map(|r| (r % ncomp) * n1 + r / ncomp).collect();
    let mp = {
        let (mut indptr, mut indices, mut data) = (vec![0usize], Vec::<usize>::new(), Vec::<f64>::new());
        for r in 0..n {
            if is_fixed[r] {
                indices.push(r);
                data.push(dscale);
            } else {
                for p in asm.k.indptr()[r]..asm.k.indptr()[r + 1] {
                    let col = asm.k.indices()[p] as usize;
                    if !is_fixed[col] {
                        indices.push(col);
                        data.push(asm.k.values()[p]);
                    }
                }
            }
            indptr.push(indices.len());
        }
        Csr::from_raw(n, n, indptr, indices, data).submatrix(&perm, &perm)
    };

    let ns = near_nullspace_mixed(mixed, dof_map, ncomp, &is_fixed, &perm);
    let opts = femcore::linalg::amg_elastic::AmgOptions::default();
    let amg_res = if dof_map.space == LagrangeSpace::Order2Serendipity {
        let nv = dof_map.n_vertices;
        let mut indptr = vec![0usize];
        let (mut indices, mut data) = (Vec::new(), Vec::new());
        for node in 0..n1 {
            if node < nv {
                indices.push(node);
                data.push(1.0);
            } else {
                let (p, q) = mixed.edges[node - nv];
                indices.extend_from_slice(&[p as usize, q as usize]);
                data.extend_from_slice(&[0.5, 0.5]);
            }
            indptr.push(indices.len());
        }
        let lin = Csr::from_raw(n1, nv, indptr, indices, data);
        femcore::linalg::amg_elastic::Amg::new_p2(&mp, ns, ncomp, &opts, &lin)
    } else {
        femcore::linalg::amg_elastic::Amg::new(&mp, ns, ncomp, &opts)
    };
    let amg = match amg_res {
        Ok(a) => a,
        Err(_) if spec.solver == LinearSolver::Auto => {
            return solve_direct(None, false);
        }
        Err(e) => return Err(e),
    };
    let setup_s = t_setup.elapsed().as_secs_f64();

    let rp: Vec<f64> = perm.iter().map(|&i| rhs[i]).collect();
    let t_solve = Instant::now();
    let res = femcore::linalg::amg_elastic::pcg(&mp, &rp, &amg, 1e-10, 800);
    let solve_s = t_solve.elapsed().as_secs_f64();
    femcore::progress::check()?;

    if !res.relative_residual.is_finite() || res.relative_residual > 1e-7 {
        if spec.solver == LinearSolver::Auto {
            return solve_direct(
                Some(format!(
                    "multigrid did not converge ({} iterations, residual {:.1e}); used the direct solver",
                    res.iterations, res.relative_residual
                )),
                true,
            );
        }
        return Err(format!(
            "the iterative solver did not converge ({} iterations, residual {:.2e}); choose the direct solver",
            res.iterations, res.relative_residual
        ));
    }

    let mut u = vec![0.0; n];
    for r in 0..n {
        u[perm[r]] = res.x[r];
    }
    for (&i, &v) in st.fixed.iter().zip(&st.vals) {
        u[i] = v;
    }
    let kff = asm.k.submatrix(&free, &free);
    let rhs_free: Vec<f64> = free.iter().map(|&i| rhs[i]).collect();
    let u_free: Vec<f64> = free.iter().map(|&i| u[i]).collect();
    let bw_rel = if n_free > 0 { backward_residual_many(&kff, &[rhs_free], &[u_free]) } else { 0.0 };
    let mut diag = crate::diag::SmsDiag::new("amg_pcg");
    diag.iterations = Some(res.iterations);
    diag.max_algebraic_residual = res.relative_residual.max(bw_rel);
    let note = format!(
        "AMG-CG: {} levels {:?}, {} iterations, residual {:.1e} (setup {:.1} s, solve {:.1} s)",
        amg.n_levels(),
        amg.sizes(),
        res.iterations,
        res.relative_residual,
        setup_s,
        solve_s
    );
    Ok((u, Some(note), diag))
}

pub(crate) use femcore::progress::lap_since as lap;

/// Density per cell; domains in `inactive` get `0.0` without requiring `"rho"`, while active domains without `"rho"` return an error.
pub(crate) fn cell_rho(
    materials: &femcore::props::Materials,
    mesh: &Discretisation,
    inactive: &[i32],
) -> Result<Vec<f64>, String> {
    let mut by_domain = BTreeMap::new();
    for d in mesh.domain_labels() {
        let v = if inactive.contains(&d) {
            0.0
        } else {
            materials.get(d, "rho").ok_or_else(|| format!("domain {d} has no value for 'rho'"))?
        };
        by_domain.insert(d, v);
    }
    Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
}

/// Solves the solid mechanics problem with a monomorphized [`SourceField<D>`] (VVUQ Pillar b MMS).
pub fn solve_solid_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    source: &S,
) -> Result<SolveResult, String> {
    crate::problem::check_simplex_straight(mesh)?;
    if mesh.dim != D {
        return Err(format!("SourceField<{D}> dimension mismatch (mesh is {}D)", mesh.dim));
    }
    spec.check_study_open_domains().map_err(|e| e.to_string())?;
    if matches!(spec.study, Study::AxisymmetricEigenfrequency { .. }) {
        return crate::cross_section::solve_axisymmetric(mesh, spec);
    }
    if matches!(spec.study, Study::Modes { .. }) {
        return crate::cross_section::solve_waveguide_modes(mesh, spec);
    }
    if spec.has_infinite_elements() {
        let raw = Mesh::from_solver_mesh(mesh)?;
        let mixed = nyx_physics::discretise_mixed(&raw).map_err(|e| e.to_string())?;
        return solve_solid_mixed_with_source::<D, S>(&raw, &mixed, spec, source);
    }
    let t0 = Instant::now();
    let mut tl = Instant::now();
    let mut result = solve_variant_with_source::<D, S>(mesh, spec, t0, &mut tl, source)?;
    result.physics = "sms".into();
    result.study = spec.study.kind().into();
    Ok(result)
}

/// Solves the problem: the three stages in order, with the variants for perfectly matched layers, Floquet pairs and low-reflecting
/// boundaries in their own modules. The result carries `physics`, `study`, `order`, `dim`, the solutions, scalars, log, time and dofs.
pub fn solve_solid(mesh: &Discretisation, spec: &Problem) -> Result<SolveResult, String> {
    match mesh.dim {
        2 => solve_solid_with_source::<2, UniformSource>(mesh, spec, &UniformSource),
        3 => solve_solid_with_source::<3, UniformSource>(mesh, spec, &UniformSource),
        d => Err(format!("unsupported mesh dimension {d}")),
    }
}

/// Native mixed-cell (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) solid solver with monomorphized [`SourceField<D>`].
pub fn solve_solid_mixed_with_source<const D: usize, S: SourceField<D>>(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
    source: &S,
) -> Result<SolveResult, String> {
    if mixed.dim != D {
        return Err(format!("SourceField<{D}> dimension mismatch (mesh is {}D)", mixed.dim));
    }
    spec.check_study_open_domains().map_err(|e| e.to_string())?;
    if matches!(spec.study, Study::AxisymmetricEigenfrequency { .. }) {
        return crate::cross_section::solve_axisymmetric_mixed(raw_mesh, mixed, spec);
    }
    if matches!(spec.study, Study::Modes { .. }) {
        return crate::cross_section::solve_waveguide_modes_mixed(raw_mesh, mixed, spec);
    }
    let t0 = Instant::now();
    let mut tl = Instant::now();
    let mut result = solve_variant_mixed_with_source::<D, S>(raw_mesh, mixed, spec, t0, &mut tl, source)?;
    result.physics = "sms".into();
    result.study = spec.study.kind().into();
    Ok(result)
}

/// Native mixed-cell (`Tri`/`Quad` and `Tet`/`Prism`/`Hex`) and Static Infinite Element solid solver.
pub fn solve_solid_mixed(raw_mesh: &Mesh, mixed: &MixedDiscretisation, spec: &Problem) -> Result<SolveResult, String> {
    match mixed.dim {
        2 => solve_solid_mixed_with_source::<2, UniformSource>(raw_mesh, mixed, spec, &UniformSource),
        3 => solve_solid_mixed_with_source::<3, UniformSource>(raw_mesh, mixed, spec, &UniformSource),
        d => Err(format!("unsupported mesh dimension {d}")),
    }
}

fn solve_variant_mixed_with_source<const D: usize, S: SourceField<D>>(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
    t0: Instant,
    tl: &mut Instant,
    source: &S,
) -> Result<SolveResult, String> {
    let order = spec.order;
    if spec.materials.has_pml() {
        return match spec.study {
            Study::Eigenfrequency { .. } => crate::pml::solve_solid_pml_mixed(raw_mesh, mixed, spec),
            Study::AxisymmetricEigenfrequency { .. } => {
                crate::cross_section::solve_axisymmetric_mixed(raw_mesh, mixed, spec)
            }
            Study::Modes { .. } => crate::cross_section::solve_waveguide_modes_mixed(raw_mesh, mixed, spec),
            _ => Err(
                "perfectly matched layers are used with the Eigenfrequency study (they give complex frequencies and Q)"
                    .into(),
            ),
        };
    }

    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let ie_domains = if matches!(spec.study, Study::Static) {
        ResolvedIeDomains::from_mixed(mixed, spec)
    } else {
        ResolvedIeDomains::default()
    };
    let st = setup_mixed_with_source::<D, S>(mixed, &dof_map, spec, &ie_domains, source)?;
    lap(tl, "setup (materials, constraints)");
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let thk = spec.thickness;
    let k = assemble_solid_mixed_with_ie(mixed, &dof_map, &st.c6, spec.plane, thk, &ie_domains, &spec.inactive);
    lap(tl, "stiffness assembly");

    match &spec.study {
        Study::AxisymmetricEigenfrequency { .. } => {
            Err("the axisymmetric study is solved by cross_section::solve_axisymmetric_mixed".into())
        }
        Study::Modes { .. } => {
            Err("the waveguide modes study is solved by cross_section::solve_waveguide_modes_mixed".into())
        }
        Study::Static => {
            let (simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
            let asm = assemble_static_mixed_with_source::<D, S>(
                mixed,
                &dof_map,
                &simplex,
                spec,
                &st,
                &ie_domains,
                k,
                source,
            )?;
            lap(tl, "loads");
            let (u, note, diag) = solve_static_mixed(mixed, &dof_map, spec, &st, &ie_domains, &asm)?;
            lap(tl, "linear solve");
            let (fields, cell_fields) = fields_of_mixed(
                mixed,
                &dof_map,
                &simplex,
                &proj,
                &u,
                "m",
                order,
                &st,
                spec.plane,
                asm.thermal.as_ref(),
                &ie_domains,
                &spec.inactive,
                true,
            );
            lap(tl, "fields (gradients, stress, projection)");
            let unorm = |i: usize| (0..d).map(|k| u[k * n1 + i] * u[k * n1 + i]).sum::<f64>().sqrt();
            let umax = (0..n1).map(unorm).fold(0.0, f64::max);
            let vmax =
                fields.iter().find(|f| f.name == "von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
            let mut log = ie_domains.diagnostics;
            log.extend(note);
            let (diagnostics, error_budget) = diag.finish();
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions: vec![Solution { label: "static".into(), fields, cell_fields, ..Default::default() }],
                scalars: vec![
                    ("max |u| (m)".into(), Scalar::Num(umax)),
                    ("max von Mises (Pa)".into(), Scalar::Num(vmax)),
                ],
                log,
                seconds: t0.elapsed().as_secs_f64(),
                dofs: n,
                diagnostics: Some(diagnostics),
                error_budget: Some(error_budget),
                ..Default::default()
            })
        }
        Study::Eigenfrequency { n_modes, f_search } => {
            let rho = mixed_cell_rho(&spec.materials, mixed, &spec.inactive)?;
            let rho_t: Vec<f64> = rho.iter().map(|r| r * if d == 2 { thk } else { 1.0 }).collect();
            let m1 = assemble_mass_mixed(mixed, &dof_map, &rho_t);
            let blocks: Vec<&Csr> = (0..d).map(|_| &m1).collect();
            let m = Csr::block_diag(&blocks);
            if !spec.periodic.is_empty() {
                return crate::floquet::solve_floquet_mixed(
                    raw_mesh, mixed, &dof_map, spec, &st.fixed, &k, &m, *n_modes, *f_search, t0,
                );
            }
            let lrb_labels: Vec<i32> = spec
                .bcs
                .iter()
                .filter(|b| matches!(b.kind, BcKind::LowReflecting))
                .flat_map(|b| b.labels.iter().copied())
                .collect();
            if !lrb_labels.is_empty() || (spec.materials.has_viscosity() && *f_search > 0.0) {
                return crate::lrb::solve_lrb_mixed(
                    raw_mesh,
                    mixed,
                    &dof_map,
                    spec,
                    &st.c6,
                    &st.fixed,
                    &k,
                    &m,
                    &rho,
                    &lrb_labels,
                    *n_modes,
                    *f_search,
                    t0,
                );
            }
            let free = free_dofs(n, &st.fixed);
            let nm = *n_modes;
            if nm < 1 || nm + 1 >= free.len() {
                return Err("number of modes must be between 1 and the number of free dofs".into());
            }
            let kff = k.submatrix(&free, &free);
            let mff = m.submatrix(&free, &free);
            let sigma = if *f_search > 0.0 { (2.0 * std::f64::consts::PI * *f_search).powi(2) } else { -1.0 };
            let pairs = eigs_sym(&kff, &mff, sigma, nm, &EigOpts::default(), None)?;

            let (simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
            let mut order_idx: Vec<usize> = (0..pairs.values.len()).collect();
            order_idx.sort_by(|&a, &b| pairs.values[a].partial_cmp(&pairs.values[b]).unwrap());
            let mut solutions = Vec::new();
            let mut scalars = Vec::new();
            for (j, &i) in order_idx.iter().enumerate() {
                let freq = pairs.values[i].max(0.0).sqrt() / (2.0 * std::f64::consts::PI);
                let mut u = vec![0.0; n];
                for (a, &dof) in free.iter().enumerate() {
                    u[dof] = pairs.vectors[i][a];
                }
                let mx = u.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
                if mx > 0.0 {
                    u.iter_mut().for_each(|x| *x /= mx);
                }
                let (fields, cell_fields) = fields_of_mixed(
                    mixed,
                    &dof_map,
                    &simplex,
                    &proj,
                    &u,
                    "a.u.",
                    order,
                    &st,
                    spec.plane,
                    None,
                    &ie_domains,
                    &spec.inactive,
                    false,
                );
                solutions.push(Solution {
                    label: format!("mode {}: {:.6} Hz", j + 1, freq),
                    fields,
                    cell_fields,
                    value: Some(freq),
                    ..Default::default()
                });
                scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
            }
            let mut log = vec![format!(
                "Krylov dimension {}{}",
                pairs.dim,
                if pairs.converged { "" } else { " (not fully converged)" }
            )];
            if st.fixed.is_empty() {
                log.push(
                    "no fixed boundary: the first six (3D) or three (2D) modes are rigid-body motions at ~0 Hz".into(),
                );
            }
            let (diagnostics, error_budget) =
                crate::diag::SmsDiag::from_eig_pairs("shift_invert_lanczos", &pairs).finish();
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions,
                scalars,
                log,
                seconds: t0.elapsed().as_secs_f64(),
                dofs: n,
                diagnostics: Some(diagnostics),
                error_budget: Some(error_budget),
                ..Default::default()
            })
        }
        Study::Driven { .. } => Err("driven study requires piezoelectric coupling".into()),
    }
}

fn solve_variant_with_source<const D: usize, S: SourceField<D>>(
    mesh: &Discretisation,
    spec: &Problem,
    t0: Instant,
    tl: &mut Instant,
    source: &S,
) -> Result<SolveResult, String> {
    let order = spec.order;
    if spec.materials.has_pml() {
        return match spec.study {
            Study::Eigenfrequency { .. } => crate::pml::solve_solid_pml(mesh, spec),
            Study::AxisymmetricEigenfrequency { .. } => crate::cross_section::solve_axisymmetric(mesh, spec),
            Study::Modes { .. } => crate::cross_section::solve_waveguide_modes(mesh, spec),
            _ => Err(
                "perfectly matched layers are used with the Eigenfrequency study (they give complex frequencies and Q)"
                    .into(),
            ),
        };
    }
    let st = setup_with_source::<D, S>(mesh, spec, source)?;
    lap(tl, "setup (materials, constraints)");
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let thk = spec.thickness;
    let k = assemble_solid(mesh, &st.c6, order, spec.plane, thk);
    lap(tl, "stiffness assembly");
    match &spec.study {
        Study::AxisymmetricEigenfrequency { .. } => {
            Err("the axisymmetric study is solved by cross_section::solve_axisymmetric".into())
        }
        Study::Modes { .. } => {
            Err("the waveguide modes study is solved by cross_section::solve_waveguide_modes".into())
        }
        Study::Static => {
            let asm = assemble_static_with_source::<D, S>(mesh, spec, &st, k, source)?;
            lap(tl, "loads");
            let (u, note, diag) = solve_static_checked(mesh, spec, &st, &asm)?;
            lap(tl, "linear solve");
            let (fields, cells) = fields_of(mesh, &u, order, &st, spec.plane, asm.thermal.as_ref(), true);
            lap(tl, "fields (gradients, stress, averaging)");
            let unorm = |i: usize| (0..d).map(|k| u[k * n1 + i] * u[k * n1 + i]).sum::<f64>().sqrt();
            let umax = (0..n1).map(unorm).fold(0.0, f64::max);
            let vmax =
                fields.iter().find(|f| f.name == "von Mises").unwrap().values.iter().cloned().fold(0.0, f64::max);
            let (diagnostics, error_budget) = diag.finish();
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions: vec![Solution {
                    label: "static".into(),
                    fields,
                    cell_fields: cell_fields_of(&cells),
                    ..Default::default()
                }],
                scalars: vec![
                    ("max |u| (m)".into(), Scalar::Num(umax)),
                    ("max von Mises (Pa)".into(), Scalar::Num(vmax)),
                ],
                log: note.into_iter().collect(),
                seconds: t0.elapsed().as_secs_f64(),
                dofs: n,
                diagnostics: Some(diagnostics),
                error_budget: Some(error_budget),
                ..Default::default()
            })
        }
        Study::Eigenfrequency { n_modes, f_search } => {
            let rho = cell_rho(&spec.materials, mesh, &spec.inactive)?;
            let rho_t: Vec<f64> = rho.iter().map(|r| r * if d == 2 { thk } else { 1.0 }).collect();
            let m1 = assemble_mass(mesh, order, &rho_t);
            let blocks: Vec<&Csr> = (0..d).map(|_| &m1).collect();
            let m = Csr::block_diag(&blocks);
            if !spec.periodic.is_empty() {
                return crate::floquet::solve_floquet(mesh, spec, &st, &k, &m, *n_modes, *f_search, t0);
            }
            let lrb_labels: Vec<i32> = spec
                .bcs
                .iter()
                .filter(|b| matches!(b.kind, BcKind::LowReflecting))
                .flat_map(|b| b.labels.iter().copied())
                .collect();
            if !lrb_labels.is_empty() || (spec.materials.has_viscosity() && *f_search > 0.0) {
                return crate::lrb::solve_lrb(mesh, spec, &st, &k, &m, &rho, &lrb_labels, *n_modes, *f_search, t0);
            }
            let free = free_dofs(n, &st.fixed);
            let nm = *n_modes;
            if nm < 1 || nm + 1 >= free.len() {
                return Err("number of modes must be between 1 and the number of free dofs".into());
            }
            let (pairs, eig_note) = solve_real_modes(mesh, spec, &st, &k, &m, &free, nm, *f_search)?;
            let mut order_idx: Vec<usize> = (0..pairs.values.len()).collect();
            order_idx.sort_by(|&a, &b| pairs.values[a].partial_cmp(&pairs.values[b]).unwrap());
            let mut solutions = vec![];
            let mut scalars = vec![];
            for (j, &i) in order_idx.iter().enumerate() {
                let freq = pairs.values[i].max(0.0).sqrt() / (2.0 * std::f64::consts::PI);
                let mut u = vec![0.0; n];
                for (a, &dof) in free.iter().enumerate() {
                    u[dof] = pairs.vectors[i][a];
                }
                let mx = u.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                if mx > 0.0 {
                    u.iter_mut().for_each(|x| *x /= mx);
                }
                let (fields, cells) = fields_of(mesh, &u, order, &st, spec.plane, None, false);
                solutions.push(Solution {
                    label: format!("mode {}: {:.6} Hz", j + 1, freq),
                    fields,
                    cell_fields: cell_fields_of(&cells),
                    value: Some(freq),
                    ..Default::default()
                });
                scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
            }
            let mut log = vec![format!(
                "Krylov dimension {}{}",
                pairs.dim,
                if pairs.converged { "" } else { " (not fully converged)" }
            )];
            let solver_kind = if eig_note.as_deref().is_some_and(|s| s.contains("LOBPCG")) {
                "lobpcg_amg"
            } else {
                "shift_invert_lanczos"
            };
            log.extend(eig_note);
            if st.fixed.is_empty() {
                log.push(
                    "no fixed boundary: the first six (3D) or three (2D) modes are rigid-body motions at ~0 Hz".into(),
                );
            }
            let (diagnostics, error_budget) = crate::diag::SmsDiag::from_eig_pairs(solver_kind, &pairs).finish();
            Ok(SolveResult {
                order: order as u8,
                dim: d as u8,
                solutions,
                scalars,
                log,
                seconds: t0.elapsed().as_secs_f64(),
                dofs: n,
                diagnostics: Some(diagnostics),
                error_budget: Some(error_budget),
                ..Default::default()
            })
        }
        Study::Driven { .. } => Err("driven study requires piezoelectric coupling".into()),
    }
}

/// Element-wise stress at the centroid of every cell, one row per cell:
/// `[sxx, sxz, szx, szz, syy, von Mises]` in 2D (the in-plane axes are x and z), `[s (3x3 row-major), von Mises]` in 3D.
pub fn centroid_stress(mesh: &Discretisation, u: &[f64], c6: &[M6], plane: Plane, order: usize) -> Vec<Vec<f64>> {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let geo = geometry(mesh);
    let lam: Vec<f64> = vec![1.0 / (d as f64 + 1.0); d + 1];
    let (_, dn) = shape(&lam, order);
    let nloc = mesh.nloc(order);
    let mut out = vec![];
    let mut dofs = Vec::new();
    for c in 0..mesh.n_cells() {
        mesh.cell_dofs(c, order, &mut dofs);
        let base = c * (d + 1) * d;
        let mut grad = vec![0.0; d * d];
        for i in 0..d {
            for k in 0..d {
                let mut s = 0.0;
                for a in 0..nloc {
                    let mut ga = 0.0;
                    for b in 0..=d {
                        ga += dn[a * (d + 1) + b] * geo.g[base + b * d + k];
                    }
                    s += u[i * n1 + dofs[a]] * ga;
                }
                grad[i * d + k] = s;
            }
        }
        let s = stress_at_nodes(&grad, 1, 1, d, &c6[c..c + 1], plane, &[])[0];
        let vm = von_mises(&s);
        out.push(if d == 2 {
            vec![s[0], s[3], s[3], s[1], s[2], vm]
        } else {
            vec![s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2], vm]
        });
    }
    out
}
