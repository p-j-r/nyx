// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electromagnetic eigenmodes of a body of revolution (ring resonators, disks, coaxial and spherical cavities)
//! with azimuthal dependence `exp(-i m theta)` on the `(r >= 0, z)` meridian half-plane (`mesh x = r`, `mesh z = z`;
//! `theta` is wafer `y`, so `m > 0` runs along `+y`).
//!
//! Unknowns: meridian `(E_r, E_z)` in Nedelec `H(curl)` (`order = 1` or `2`) paired with the covariant scaled
//! azimuthal potential `\Psi = i r E_theta` (`E_theta = -i \Psi / r`) in exact de Rham partner Lagrange `H^1`
//! (`P_1`/`P_2` on triangles, `Q_1`/`Q_2` `Order2Tensor` on quadrilaterals). Using `\Psi = i r E_theta`:
//! - eliminates `1/r` dipole locking at `|m| = 1` and keeps every discrete gradient `(grad phi, m phi)` strictly
//!   in `ker(curl_m)` at `\lambda = 0`,
//! - enforces exact `m`-dependent regularity on exterior `r = 0` symmetry-axis facets (`\Psi = 0` for all `m`;
//!   meridian tangential `E_z` free for `m = 0` and fixed to `0` for `|m| >= 1`),
//! - makes `A` and `B` real symmetric without PML and complex symmetric under Cylindrical / Spherical / Cartesian PMLs
//!   evaluated point-by-point via `femcore::pml::{pml_jacobian, pullback_maxwell}`.

use crate::ibc::{add_ibc_2d_axisym, add_ibc_2d_axisym_mixed, pec_edges_and_nodes_2d};
use crate::mixed_ops::{complex_eps_c3, iso_mu_inv_c3, EmPmlTable, MixedCellMaterials, C3};
use crate::modes::C0;
use crate::problem::{Problem, Study, SurfaceImpedanceBc};
use femcore::edge::{cell_data, edge_basis, tri_g, CellData};
use femcore::edge_mixed::{eval_cell_basis_2d, nedelec_nloc, NedEvalScratch, NedMixed};
use femcore::fem::*;
use femcore::fem_mixed::{lagrange_shape_into, quad_rule_mixed, IsoparametricCell, LagrangeSpace, MixedLagrangeDofMap};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::mesh::mixed::MixedDiscretisation;
use femcore::optical::{fine_rule, EPS0, MU0};
use femcore::pml::{radial_coordinate, stretch};
use femcore::props::{Materials, Pml, PmlDir, TensorKind};
use femcore::section::plane_tensor;
use nyx_physics::{
    cell_fields_of, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use std::time::Instant;

/// Material data of every cell of an axisymmetric problem in `(0 = r, 1 = z, 2 = theta)` order.
#[derive(Clone, Debug)]
pub struct AxiMaterial {
    pub eps: Vec<[C; 3]>,
    pub mu: Vec<[C; 3]>,
    pub eps_c3: Vec<C3>,
    pub mu_inv_c3: Vec<C3>,
    pub pml: Option<EmPmlTable>,
}

fn diag_to_c3(d: [C; 3]) -> C3 {
    let z = c(0.0, 0.0);
    [[d[0], z, z], [z, d[1], z], [z, z, d[2]]]
}

fn inv_diag_to_c3(d: [C; 3]) -> C3 {
    let one = c(1.0, 0.0);
    let z = c(0.0, 0.0);
    [[one / d[0], z, z], [z, one / d[1], z], [z, z, one / d[2]]]
}

impl AxiMaterial {
    /// Isotropic `eps_r = n^2` and `mu_r = 1` per domain (`index`: refractive index per domain label, default 1).
    pub fn from_indices(mesh: &Discretisation, index: &std::collections::BTreeMap<i32, f64>) -> AxiMaterial {
        let eps: Vec<[C; 3]> = (0..mesh.n_cells())
            .map(|ci| {
                let n = index.get(&mesh.cell_labels[ci]).copied().unwrap_or(1.0);
                [c(n * n, 0.0); 3]
            })
            .collect();
        let mu = vec![[c(1.0, 0.0); 3]; eps.len()];
        let eps_c3 = eps.iter().copied().map(diag_to_c3).collect();
        let mu_inv_c3 = mu.iter().copied().map(inv_diag_to_c3).collect();
        AxiMaterial { eps, mu, eps_c3, mu_inv_c3, pml: None }
    }

    /// `eps_r` and `mu_r` per domain from the problem's materials (default 1), absorbing layers (`Materials::pml`) included.
    pub fn from_materials(mesh: &Discretisation, m: &Materials) -> Result<AxiMaterial, String> {
        Self::from_materials_at_freq(mesh, m, 0.0)
    }

    /// Like [`Self::from_materials`], passing `f_ref` for `auto_scale` PML calibration.
    pub fn from_materials_at_freq(mesh: &Discretisation, m: &Materials, f_ref: f64) -> Result<AxiMaterial, String> {
        let nc = mesh.n_cells();
        let tens = m.cell_complex_permittivity(mesh, Some(1.0))?;
        let mur = m.cell_values(mesh, "mu_r", Some(1.0))?;
        let mut eps = Vec::with_capacity(nc);
        let mut mu = Vec::with_capacity(nc);
        let mut eps_c3 = Vec::with_capacity(nc);
        let mut mu_inv_c3 = Vec::with_capacity(nc);
        for ci in 0..nc {
            let re_rz = plane_tensor(&tens[ci].0, 2);
            let im_rz = plane_tensor(&tens[ci].1, 2);
            let ec = complex_eps_c3(&re_rz, &im_rz);
            eps.push([ec[0][0], ec[1][1], ec[2][2]]);
            mu.push([c(mur[ci], 0.0); 3]);
            eps_c3.push(ec);
            mu_inv_c3.push(iso_mu_inv_c3(mur[ci]));
        }
        let pml_table = EmPmlTable::from_simplex(mesh, m, f_ref, true, false);
        let pml = pml_table.has_any().then_some(pml_table);
        Ok(AxiMaterial { eps, mu, eps_c3, mu_inv_c3, pml })
    }

    /// Legacy cell-centroid Cartesian `r`/`z` stretch helper retained for callers constructing `AxiMaterial` manually.
    pub fn with_pml(mut self, mesh: &Discretisation, pml: &std::collections::BTreeMap<i32, Pml>) -> AxiMaterial {
        let mut boxes: std::collections::BTreeMap<i32, ([f64; 2], [f64; 2])> = Default::default();
        for ci in 0..mesh.n_cells() {
            let l = mesh.cell_labels[ci];
            if pml.get(&l).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off) {
                let e = boxes.entry(l).or_insert(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]));
                for &v in mesh.cell(ci) {
                    let pt = mesh.point(v as usize);
                    for k in 0..2 {
                        e.0[k] = e.0[k].min(pt[k]);
                        e.1[k] = e.1[k].max(pt[k]);
                    }
                }
            }
        }
        for ci in 0..mesh.n_cells() {
            let l = mesh.cell_labels[ci];
            let (Some(p), Some((lo, hi))) = (pml.get(&l), boxes.get(&l)) else { continue };
            let mut cen = [0.0; 2];
            for &v in mesh.cell(ci) {
                for k in 0..2 {
                    cen[k] += mesh.point(v as usize)[k] / 3.0;
                }
            }
            let (sr, sz) = stretch(p, [p.x, p.z], cen, *lo, *hi);
            let r = cen[0].max(1e-30);
            let rt = radial_coordinate(p, r, lo[0], hi[0]);
            let st = rt / r;
            let lam = [st * sz / sr, sr * st / sz, sr * sz / st];
            for k in 0..3 {
                self.eps[ci][k] *= lam[k];
                self.mu[ci][k] *= lam[k];
            }
            self.eps_c3[ci] = diag_to_c3(self.eps[ci]);
            self.mu_inv_c3[ci] = inv_diag_to_c3(self.mu[ci]);
        }
        self.pml = None;
        self
    }
}

/// Modes of one azimuthal order.
pub struct AxiModes {
    pub m: f64,
    /// real parts of the frequencies (Hz), ascending in distance from the shift
    pub freqs: Vec<f64>,
    /// complex frequencies (Hz): the imaginary part is the radiation loss (with absorbing layers), `Q = Re / (2 |Im|)`
    pub freqs_c: Vec<C>,
    /// coefficient vectors over all dofs (`net` edge dofs for `(E_r, E_z)` then Lagrange dofs of `\Psi = i r E_theta`)
    pub vectors: Vec<Vec<C>>,
    pub order: usize,
    pub converged: bool,
}

/// Fields of one mode at a point in `(0 = r, 1 = z, 2 = theta)` order.
#[derive(Clone, Copy, Debug, Default)]
pub struct Fields3 {
    pub e: [C; 3],
    pub d: [C; 3],
    pub h: [C; 3],
}

pub(crate) struct AxiCtx<'a> {
    pub mesh: &'a Discretisation,
    pub geo: Geometry,
    pub cd: CellData,
    pub order: usize,
    pub r_tol: f64,
}

impl<'a> AxiCtx<'a> {
    pub fn new(mesh: &'a Discretisation, order: usize) -> Self {
        let (lo, hi) = mesh.bbox();
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-30);
        AxiCtx { mesh, geo: geometry(mesh), cd: cell_data(mesh, order), order, r_tol: 1e-11 * extent }
    }

    pub fn radius_and_z(&self, ci: usize, l: &[f64; 3]) -> [f64; 2] {
        let v = self.mesh.cell(ci);
        let r = (0..3).map(|k| l[k] * self.mesh.point(v[k] as usize)[0]).sum();
        let z = (0..3).map(|k| l[k] * self.mesh.point(v[k] as usize)[1]).sum();
        [r, z]
    }

    /// Evaluates `E`, `D`, and `H` in `(r, z, theta)` order from `(E_rz, \Psi = i r E_theta)` at barycentric point `l`.
    pub fn fields(&self, x: &[C], m: f64, omega: f64, mat: &AxiMaterial, ci: usize, l: &[f64; 3]) -> Fields3 {
        let nloc_e = if self.order == 1 { 3 } else { 8 };
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        let g = tri_g(&self.geo, ci);
        edge_basis(l, &g, &self.cd.sign[ci], self.order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (mut er, mut ez, mut cth) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        for s in 0..nloc_e {
            let cf = x[self.cd.gdofs[ci][s]];
            er += cf * w[s][0];
            ez += cf * w[s][1];
            cth += cf * curl[s];
        }
        let (vals, dn) = shape(l, self.order);
        let mut pd = Vec::with_capacity(6);
        self.mesh.cell_dofs(ci, self.order, &mut pd);
        let (mut psi, mut dpsi_dr, mut dpsi_dz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        let gb = ci * 6;
        for a in 0..vals.len() {
            let cf = x[self.cd.net + pd[a]];
            psi += cf * vals[a];
            let (mut gx, mut gy) = (0.0, 0.0);
            for b in 0..3 {
                gx += dn[a * 3 + b] * self.geo.g[gb + b * 2];
                gy += dn[a * 3 + b] * self.geo.g[gb + b * 2 + 1];
            }
            dpsi_dr += cf * gx;
            dpsi_dz += cf * gy;
        }
        let rz = self.radius_and_z(ci, l);
        let r = rz[0];
        let m_abs_round = m.abs().round() as i64;
        let (et, cr, cz) = if r > self.r_tol {
            let et_v = c(0.0, -1.0) * psi / r;
            let cr_v = c(0.0, 1.0) * (dpsi_dz - ez * m) / r;
            let cz_v = c(0.0, 1.0) * (-dpsi_dr + er * m) / r;
            (et_v, cr_v, cz_v)
        } else {
            // Exact L'Hopital regularity limit on the symmetry axis `r = 0`:
            let et_axis = if m_abs_round == 1 { c(0.0, -1.0) * dpsi_dr } else { c(0.0, 0.0) };
            if m_abs_round != 1 {
                er = c(0.0, 0.0);
            }
            if m_abs_round != 0 {
                ez = c(0.0, 0.0);
            }
            (et_axis, c(0.0, 0.0), c(0.0, 0.0))
        };
        let ct = -cth;
        let (ep, mui) = match &mat.pml {
            Some(p) => {
                let (e_t, m_t, _) = p.eval_2d(ci, mat.eps_c3[ci], mat.mu_inv_c3[ci], rz, true);
                (e_t, m_t)
            }
            None => (mat.eps_c3[ci], mat.mu_inv_c3[ci]),
        };
        let q = c(0.0, 1.0) / (omega.max(1e-300) * MU0);
        let curl_vec = [cr, cz, ct];
        let mut h = [c(0.0, 0.0); 3];
        let e = [er, ez, et];
        let mut d = [c(0.0, 0.0); 3];
        for a in 0..3 {
            for b in 0..3 {
                h[a] += mui[a][b] * curl_vec[b] * q;
                d[a] += ep[a][b] * e[b] * EPS0;
            }
        }
        Fields3 { e, d, h }
    }
}

/// Modes with frequency nearest `f_target` at azimuthal order `m` (`nev` of them).
pub fn solve_axisym_modes(
    mesh: &Discretisation,
    mat: &AxiMaterial,
    m: f64,
    f_target: f64,
    nev: usize,
    order: usize,
) -> Result<AxiModes, String> {
    let asm = assemble_axisym(mesh, mat, m, order)?;
    solve_axisym_pairs(&asm, f_target, nev)
}

/// Stage 1: curl-curl matrix `A` and permittivity mass `B` of `A x = (omega/c)^2 B x` on free DOFs.
pub struct AxiAssembly {
    pub m: f64,
    pub order: usize,
    pub a: CCsr,
    pub b: CCsr,
    pub free: Vec<usize>,
    pub n_tot: usize,
}

pub fn assemble_axisym(mesh: &Discretisation, mat: &AxiMaterial, m: f64, order: usize) -> Result<AxiAssembly, String> {
    assemble_axisym_with_bcs(mesh, mat, m, order, 0.0, &[], &[], &[])
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn add_axisym_point(
    nloc_e: usize,
    nl: usize,
    m: f64,
    r: f64,
    wr: f64,
    w: &[[f64; 2]],
    curl: &[f64],
    vals: &[f64],
    gn: &[[f64; 2]],
    ep: &C3,
    mui: &C3,
    cb: &mut [[C; 3]],
    eb: &mut [[C; 3]],
    kl: &mut [C],
    ml: &mut [C],
) {
    let nb = nloc_e + nl;
    let inv_r = 1.0 / r;
    let m_inv_r = m * inv_r;
    for s in 0..nloc_e {
        cb[s] = [c(0.0, -m_inv_r * w[s][1]), c(0.0, m_inv_r * w[s][0]), c(-curl[s], 0.0)];
        eb[s] = [c(w[s][0], 0.0), c(w[s][1], 0.0), c(0.0, 0.0)];
    }
    for a in 0..nl {
        cb[nloc_e + a] = [c(0.0, gn[a][1] * inv_r), c(0.0, -gn[a][0] * inv_r), c(0.0, 0.0)];
        eb[nloc_e + a] = [c(0.0, 0.0), c(0.0, 0.0), c(0.0, -vals[a] * inv_r)];
    }
    for i in 0..nb {
        let ci_ = cb[i];
        let ei = eb[i];
        for j in 0..nb {
            let cj = cb[j];
            let ej = eb[j];
            let mut kk = c(0.0, 0.0);
            let mut mm = c(0.0, 0.0);
            for p in 0..3 {
                let cip = ci_[p].conj();
                let eip = ei[p].conj();
                for q in 0..3 {
                    kk += cip * mui[p][q] * cj[q];
                    mm += eip * ep[p][q] * ej[q];
                }
            }
            kl[i * nb + j] += kk * wr;
            ml[i * nb + j] += mm * wr;
        }
    }
}

pub fn assemble_axisym_with_bcs(
    mesh: &Discretisation,
    mat: &AxiMaterial,
    m: f64,
    order: usize,
    f_search: f64,
    pec_labels: &[i32],
    ibc_bcs: &[SurfaceImpedanceBc],
    pmc_labels: &[i32],
) -> Result<AxiAssembly, String> {
    if mesh.dim != 2 {
        return Err("an axisymmetric problem needs a 2D mesh in (r, z)".into());
    }
    if order != 1 && order != 2 {
        return Err("element order must be 1 or 2".into());
    }
    let (lo, hi) = mesh.bbox();
    let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-30);
    let r_axis_tol = 1e-9 * extent;
    if lo[0] < -r_axis_tol {
        return Err(format!(
            "the mesh of an axisymmetric problem must lie in r >= 0 (found minimum r = {:.3e})",
            lo[0]
        ));
    }

    let ctx = AxiCtx::new(mesh, order);
    let (nc, nv, ne) = (mesh.n_cells(), mesh.n_points(), mesh.n_edges());
    let net = ctx.cd.net;
    let nth = mesh.n_dofs(order);
    let n_tot = net + nth;
    let nloc_e = if order == 1 { 3 } else { 8 };
    let nl = mesh.nloc(order);
    let nb = nloc_e + nl;
    let (pts, wts) = fine_rule();
    let tab: Vec<(Vec<f64>, Vec<f64>)> = pts.iter().map(|l| shape(l, order)).collect();

    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut gn = vec![[0.0; 2]; nl];
    let mut kl = vec![c(0.0, 0.0); nb * nb];
    let mut ml = vec![c(0.0, 0.0); nb * nb];
    let mut cb = vec![[c(0.0, 0.0); 3]; nb];
    let mut eb = vec![[c(0.0, 0.0); 3]; nb];

    for ci in 0..nc {
        kl.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        ml.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let g = tri_g(&ctx.geo, ci);
        let area = ctx.geo.vol[ci];
        let eps_0 = mat.eps_c3[ci];
        let mui_0 = mat.mu_inv_c3[ci];
        for (q, (l, &wq)) in pts.iter().zip(&wts).enumerate() {
            let rz = ctx.radius_and_z(ci, l);
            let r = rz[0];
            if r <= 0.0 {
                return Err("degenerate axisymmetric triangle with non-positive interior radius".into());
            }
            edge_basis(l, &g, &ctx.cd.sign[ci], order, &mut w, &mut curl);
            let (vals, dn) = &tab[q];
            for a in 0..nl {
                gn[a][0] = (0..3).map(|b| dn[a * 3 + b] * g[b][0]).sum();
                gn[a][1] = (0..3).map(|b| dn[a * 3 + b] * g[b][1]).sum();
            }
            let (ep, mui) = match &mat.pml {
                Some(p) => {
                    let (e_t, m_t, _) = p.eval_2d(ci, eps_0, mui_0, rz, true);
                    (e_t, m_t)
                }
                None => (eps_0, mui_0),
            };
            add_axisym_point(
                nloc_e,
                nl,
                m,
                r,
                wq * area * r,
                &w,
                &curl,
                vals,
                &gn,
                &ep,
                &mui,
                &mut cb,
                &mut eb,
                &mut kl,
                &mut ml,
            );
        }
        mesh.cell_dofs(ci, order, &mut pd);
        let mut dofs: Vec<usize> = ctx.cd.gdofs[ci].clone();
        dofs.extend(pd.iter().map(|&p| net + p));
        a_coo.push_block(&dofs, &dofs, &kl);
        b_coo.push_block(&dofs, &dofs, &ml);
    }

    add_ibc_2d_axisym(mesh, order, f_search, pec_labels, ibc_bcs, &mut a_coo)?;
    let a_full: CCsr = a_coo.to_csr();
    let b_full: CCsr = b_coo.to_csr();

    // Separate exterior symmetry-axis facets (`r = 0`) from exterior PEC walls!
    let exterior = mesh.exterior_facets();
    let ibc_labels: Vec<i32> = ibc_bcs.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let mut pec_edges = Vec::new();
    let mut pec_nodes = Vec::new();
    let mut axis_edges = Vec::new();
    let mut axis_nodes = Vec::new();

    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        let on_axis =
            exterior[f] && !pec_labels.contains(&l) && fv.iter().all(|&v| mesh.point(v as usize)[0] <= r_axis_tol);
        if on_axis {
            axis_edges.push(mesh.edge_index(fv[0], fv[1]));
            axis_nodes.push(fv[0] as usize);
            axis_nodes.push(fv[1] as usize);
            continue;
        }
        let is_pec = pec_labels.contains(&l) || (exterior[f] && !ibc_labels.contains(&l) && !pmc_labels.contains(&l));
        if is_pec {
            pec_edges.push(mesh.edge_index(fv[0], fv[1]));
            pec_nodes.push(fv[0] as usize);
            pec_nodes.push(fv[1] as usize);
        }
    }
    let _ = pec_edges_and_nodes_2d;

    let mut fixed_set = std::collections::BTreeSet::new();
    // 1. PEC walls: tangential E_rz = 0 and Psi = i r E_theta = 0
    for &e in &pec_edges {
        fixed_set.insert(e);
        if order == 2 {
            fixed_set.insert(ne + e);
            fixed_set.insert(net + nv + e);
        }
    }
    for &v in &pec_nodes {
        fixed_set.insert(net + v);
    }
    // 2. Symmetry axis `r = 0`:
    //    - `\Psi = i r E_theta` vanishes on `r = 0` for EVERY `m`!
    //    - Meridian tangential `E_z` is FREE at `m = 0` and ZERO for `|m| >= 1`!
    for &v in &axis_nodes {
        fixed_set.insert(net + v);
    }
    for &e in &axis_edges {
        if order == 2 {
            fixed_set.insert(net + nv + e);
        }
        if m.abs() >= 0.5 {
            fixed_set.insert(e);
            if order == 2 {
                fixed_set.insert(ne + e);
            }
        }
    }
    // Also catch any isolated vertex touching `r = 0` (e.g. a single triangle apex on the axis):
    for v in 0..nv {
        if mesh.point(v)[0] <= r_axis_tol {
            fixed_set.insert(net + v);
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let free = free_dofs(n_tot, &fixed);
    Ok(AxiAssembly { m, order, a: a_full.submatrix(&free, &free), b: b_full.submatrix(&free, &free), free, n_tot })
}

/// Stage 2: the `nev` eigenpairs nearest `f_target` (Hz), sorted by closeness to it.
pub fn solve_axisym_pairs(asm: &AxiAssembly, f_target: f64, nev: usize) -> Result<AxiModes, String> {
    let (m, order, free, n_tot) = (asm.m, asm.order, &asm.free, asm.n_tot);
    let k0 = 2.0 * std::f64::consts::PI * f_target / C0;
    let sigma = k0 * k0;
    let pairs = eigs_complex(&asm.a, &asm.b, c(sigma, 0.0), nev, 1e-9, 800)?;
    let mut modes: Vec<(C, Vec<C>)> = vec![];
    for (i, lam) in pairs.values.iter().enumerate() {
        if lam.re <= 1e-6 * sigma.max(1e-30) || pairs.residuals.get(i).is_some_and(|r| *r > 1e-5) {
            continue;
        }
        let f = lam.sqrt() * (C0 / (2.0 * std::f64::consts::PI));
        let f = if f.re < 0.0 { -f } else { f };
        let mut x = vec![c(0.0, 0.0); n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = pairs.vectors[i][a];
        }
        modes.push((f, x));
    }
    modes.sort_by(|a, b| (a.0.re - f_target).abs().partial_cmp(&(b.0.re - f_target).abs()).unwrap());
    modes.truncate(nev);
    if modes.is_empty() {
        return Err("no mode converged near the target frequency".into());
    }
    Ok(AxiModes {
        m,
        freqs: modes.iter().map(|x| x.0.re).collect(),
        freqs_c: modes.iter().map(|x| x.0).collect(),
        vectors: modes.into_iter().map(|x| x.1).collect(),
        order,
        converged: pairs.converged,
    })
}

/// Power `P = 2 Int Re(conj(E_z) H_r - conj(E_r) H_z) dA` of mode `j` along `+theta`.
pub fn power_along_theta(mesh: &Discretisation, mat: &AxiMaterial, modes: &AxiModes, j: usize) -> f64 {
    let ctx = AxiCtx::new(mesh, modes.order);
    let (pts, wts) = fine_rule();
    let omega = 2.0 * std::f64::consts::PI * modes.freqs[j];
    let mut p = 0.0;
    for ci in 0..mesh.n_cells() {
        for (l, w) in pts.iter().zip(&wts) {
            let f = ctx.fields(&modes.vectors[j], modes.m, omega, mat, ci, l);
            p += 2.0 * (f.e[1].conj() * f.h[0] - f.e[0].conj() * f.h[1]).re * w * ctx.geo.vol[ci];
        }
    }
    p
}

/// Stage 3 on `SimplexMesh`: plot fields (`|E|`, `Re Ex`, `Re Ey`, `Re Ez`) and complex cell fields (`Ex`, `Ey`, `Ez`).
pub fn axisym_fields_cells(
    mesh: &Discretisation,
    mat: &AxiMaterial,
    modes: &AxiModes,
    j: usize,
) -> (Vec<Field>, CellFields) {
    let order = modes.order;
    let ctx = AxiCtx::new(mesh, order);
    let (nc, nloc) = (mesh.n_cells(), mesh.nloc(order));
    let lam = local_bary(2, order);
    let omega = 2.0 * std::f64::consts::PI * modes.freqs[j];
    let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    for ci in 0..nc {
        for (a, l) in lam.iter().enumerate() {
            e[ci * nloc + a] = ctx.fields(&modes.vectors[j], modes.m, omega, mat, ci, &[l[0], l[1], l[2]]).e;
        }
    }
    let big =
        e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, order);
    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
    ];
    let mut cells = CellFields::new(mesh, order);
    for (k, name) in [(0, "Ex"), (2, "Ey"), (1, "Ez")] {
        let z: Vec<C> = e.iter().map(|v| v[k] * rot * s).collect();
        cells.push_complex(name, "a.u.", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }
    (fields, cells)
}

/// Solves an `AxisymmetricModes` study on a `SimplexMesh` (`Discretisation`).
pub fn solve_axisymmetric(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let Study::AxisymmetricModes { f_search, azimuthal_order, n_modes } = problem.study.clone() else {
        return Err("an axisymmetric problem needs an AxisymmetricModes study".into());
    };
    if f_search <= 0.0 {
        return Err("the search frequency must be positive".into());
    }
    if n_modes < 1 {
        return Err("at least one mode is needed".into());
    }
    let mat = AxiMaterial::from_materials_at_freq(mesh, &problem.materials, f_search)?;
    let asm = assemble_axisym_with_bcs(
        mesh,
        &mat,
        azimuthal_order as f64,
        problem.order,
        f_search,
        &problem.pec,
        &problem.surface_impedance,
        &problem.pmc,
    )?;
    let modes = solve_axisym_pairs(&asm, f_search, n_modes)?;
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, f) in modes.freqs_c.iter().enumerate() {
        let (fields, cells) = axisym_fields_cells(mesh, &mat, &modes, j);
        let lossy = f.im.abs() > 1e-12 * f.re.abs().max(1.0);
        let label = if lossy {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, f.re / (2.0 * f.im.abs()))
        } else {
            format!("mode {}: {:.6} GHz", j + 1, f.re / 1e9)
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(f.re),
            eigenvalue: lossy.then_some([f.re, f.im]),
        });
        if lossy {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(f.re / (2.0 * f.im.abs()))));
        } else {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(f.re)));
        }
    }
    let mut log = vec![format!("azimuthal order {azimuthal_order}, {} dofs, shift {:.4e} Hz", asm.n_tot, f_search)];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    let _ = TensorKind::Permittivity;
    Ok(SolveResult {
        physics: "ems".into(),
        study: "axisymmetric_modes".into(),
        order: problem.order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.n_tot,
    })
}

// ------------------------------------------------------------------ Native 2D MixedDiscretisation (`Tri + Quad` + curved) axisymmetric solver

/// Solves a 2D `AxisymmetricModes` study natively on `Tri + Quad` (including curved 2D shells and `r = 0` Spherical PMLs).
pub fn solve_axisymmetric_mixed(
    mesh: &Mesh,
    mixed: &MixedDiscretisation,
    problem: &Problem,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let Study::AxisymmetricModes { f_search, azimuthal_order, n_modes } = problem.study.clone() else {
        return Err("an axisymmetric problem needs an AxisymmetricModes study".into());
    };
    if f_search <= 0.0 {
        return Err("the search frequency must be positive".into());
    }
    if n_modes < 1 {
        return Err("at least one mode is needed".into());
    }
    let m = azimuthal_order as f64;
    let (asm, ned, phi_map) = assemble_axisym_mixed(mixed, problem, m, f_search)?;
    let modes = solve_axisym_pairs(&asm, f_search, n_modes)?;

    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, mixed).map_err(|e| e.to_string())?;
    let mut solutions = Vec::with_capacity(modes.freqs_c.len());
    let mut scalars = Vec::new();
    for (j, f) in modes.freqs_c.iter().enumerate() {
        let (fields, cells) = axisym_fields_cells_mixed(mixed, &ned, &phi_map, &simplex, &proj, &modes, j);
        let lossy = f.im.abs() > 1e-12 * f.re.abs().max(1.0);
        let label = if lossy {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, f.re / (2.0 * f.im.abs()))
        } else {
            format!("mode {}: {:.6} GHz", j + 1, f.re / 1e9)
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(f.re),
            eigenvalue: lossy.then_some([f.re, f.im]),
        });
        if lossy {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(f.re / (2.0 * f.im.abs()))));
        } else {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(f.re)));
        }
    }
    let mut log = vec![format!(
        "axisymmetric mixed order {}, azimuthal order {azimuthal_order}, {} dofs, shift {:.4e} Hz",
        problem.order, asm.n_tot, f_search
    )];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    Ok(SolveResult {
        physics: "ems".into(),
        study: "axisymmetric_modes".into(),
        order: problem.order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.n_tot,
    })
}

pub fn assemble_axisym_mixed(
    mesh: &MixedDiscretisation,
    problem: &Problem,
    m: f64,
    f_search: f64,
) -> Result<(AxiAssembly, NedMixed, MixedLagrangeDofMap), String> {
    if mesh.dim != 2 {
        return Err("an axisymmetric problem needs a 2D mesh in (r, z)".into());
    }
    let order = problem.order;
    if order != 1 && order != 2 {
        return Err("element order must be 1 or 2".into());
    }
    let (lo, hi) = mesh.bbox();
    let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-30);
    let r_axis_tol = 1e-9 * extent;
    if lo[0] < -r_axis_tol {
        return Err(format!(
            "the mesh of an axisymmetric problem must lie in r >= 0 (found minimum r = {:.3e})",
            lo[0]
        ));
    }

    let ned = NedMixed::new(mesh, order);
    let phi_space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Tensor };
    let phi_map = MixedLagrangeDofMap::new(mesh, phi_space);
    let net = ned.n_dofs;
    let nth = phi_map.n_dofs;
    let n_tot = net + nth;

    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let pml_table = EmPmlTable::from_mixed(mesh, &problem.materials, f_search, true, false);
    let nc = mesh.n_cells();
    let mut eps_rz = Vec::with_capacity(nc);
    let mut mui_rz = Vec::with_capacity(nc);
    for ci in 0..nc {
        let re_rz = plane_tensor(&mats.eps_phys_re[ci], 2);
        let im_rz = plane_tensor(&mats.eps_phys_im[ci], 2);
        eps_rz.push(complex_eps_c3(&re_rz, &im_rz));
        mui_rz.push(iso_mu_inv_c3(mats.mur_phys[ci]));
    }

    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut vals = Vec::new();
    let mut dxi = Vec::new();
    let mut gn = Vec::new();
    let mut cb = Vec::new();
    let mut eb = Vec::new();
    let mut kl = Vec::new();
    let mut ml = Vec::new();

    for ci in 0..nc {
        let iso = IsoparametricCell::from_mesh(mesh, ci);
        let rule = quad_rule_mixed(iso.cell_type, 6);
        let nloc_e = nedelec_nloc(iso.cell_type, order);
        let pd = phi_map.cell_dofs(ci);
        let nl = pd.len();
        let nb = nloc_e + nl;
        gn.resize(nl, [0.0; 2]);
        cb.resize(nb, [c(0.0, 0.0); 3]);
        eb.resize(nb, [c(0.0, 0.0); 3]);
        kl.clear();
        kl.resize(nb * nb, c(0.0, 0.0));
        ml.clear();
        ml.resize(nb * nb, c(0.0, 0.0));

        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_2d(mesh, &ned, ci, &iso, xi, &mut scratch, &mut w, &mut curl);
            let r = jac.x[0];
            if r <= 0.0 {
                return Err("degenerate axisymmetric cell with non-positive interior radius".into());
            }
            lagrange_shape_into(iso.cell_type, phi_space, xi, &mut vals, &mut dxi);
            for a in 0..nl {
                let g3 = jac.map_grad(dxi[a]);
                gn[a] = [g3[0], g3[1]];
            }
            let (ep, mui, _) = pml_table.eval_2d(ci, eps_rz[ci], mui_rz[ci], [jac.x[0], jac.x[1]], true);
            add_axisym_point(
                nloc_e,
                nl,
                m,
                r,
                w_ref * jac.det.abs() * r,
                &w,
                &curl,
                &vals,
                &gn,
                &ep,
                &mui,
                &mut cb,
                &mut eb,
                &mut kl,
                &mut ml,
            );
        }

        let mut dofs = ned.cell_dofs(ci).to_vec();
        dofs.extend(pd.iter().map(|&p| net + p));
        a_coo.push_block(&dofs, &dofs, &kl);
        b_coo.push_block(&dofs, &dofs, &ml);
    }

    add_ibc_2d_axisym_mixed(mesh, &ned, &phi_map, f_search, &problem.pec, &problem.surface_impedance, &mut a_coo)?;
    let a_full = a_coo.to_csr();
    let b_full = b_coo.to_csr();

    let nv = mesh.n_points();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let mut fixed_set = std::collections::BTreeSet::new();

    for f in 0..mesh.n_facets() {
        let l = mesh.facet_labels[f];
        let fv = mesh.facet(f);
        let e = mesh.edge_index(fv[0], fv[1]).expect("2D facet edge");
        let on_axis = mesh.exterior_facets[f]
            && !problem.pec.contains(&l)
            && fv.iter().all(|&v| mesh.point(v as usize)[0] <= r_axis_tol);
        if on_axis {
            fixed_set.insert(net + fv[0] as usize);
            fixed_set.insert(net + fv[1] as usize);
            if order == 2 {
                fixed_set.insert(net + nv + e);
            }
            if m.abs() >= 0.5 {
                if order == 1 {
                    fixed_set.insert(e);
                } else {
                    fixed_set.insert(2 * e);
                    fixed_set.insert(2 * e + 1);
                }
            }
            continue;
        }
        let is_pec = problem.pec.contains(&l)
            || (mesh.exterior_facets[f] && !ibc_labels.contains(&l) && !problem.pmc.contains(&l));
        if is_pec {
            if order == 1 {
                fixed_set.insert(e);
            } else {
                fixed_set.insert(2 * e);
                fixed_set.insert(2 * e + 1);
                fixed_set.insert(net + nv + e);
            }
            fixed_set.insert(net + fv[0] as usize);
            fixed_set.insert(net + fv[1] as usize);
        }
    }
    for v in 0..nv {
        if mesh.point(v)[0] <= r_axis_tol {
            fixed_set.insert(net + v);
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let free = free_dofs(n_tot, &fixed);
    Ok((
        AxiAssembly { m, order, a: a_full.submatrix(&free, &free), b: b_full.submatrix(&free, &free), free, n_tot },
        ned,
        phi_map,
    ))
}

#[allow(clippy::too_many_arguments)]
fn axisym_fields_cells_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    phi_map: &MixedLagrangeDofMap,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    modes: &AxiModes,
    j: usize,
) -> (Vec<Field>, CellFields) {
    let order = modes.order;
    let m_abs_round = modes.m.abs().round() as i64;
    let (lo, hi) = mixed.bbox();
    let r_tol = 1e-11 * (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-30);
    let net = ned.n_dofs;
    let x = &modes.vectors[j];
    let isos: Vec<IsoparametricCell> = (0..mixed.n_cells()).map(|ci| IsoparametricCell::from_mesh(mixed, ci)).collect();

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut curl = Vec::new();
    let mut vals = Vec::new();
    let mut dxi = Vec::new();

    let e = proj.sample_simplex_cell_nodes(mixed, order, |mc, xi| {
        let jac = eval_cell_basis_2d(mixed, ned, mc, &isos[mc], xi, &mut scratch, &mut w, &mut curl);
        lagrange_shape_into(isos[mc].cell_type, phi_map.space, xi, &mut vals, &mut dxi);
        let cdofs = ned.cell_dofs(mc);
        let pd = phi_map.cell_dofs(mc);
        let (mut er, mut ez, mut psi, mut dpsi_dr) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        for (s, &gd) in cdofs.iter().enumerate() {
            let cf = x[gd];
            er += cf * w[s][0];
            ez += cf * w[s][1];
        }
        for (a, &gp) in pd.iter().enumerate() {
            let cf = x[net + gp];
            psi += cf * vals[a];
            let g3 = jac.map_grad(dxi[a]);
            dpsi_dr += cf * g3[0];
        }
        let r = jac.x[0];
        let et = if r > r_tol {
            c(0.0, -1.0) * psi / r
        } else {
            if m_abs_round != 1 {
                er = c(0.0, 0.0);
            }
            if m_abs_round != 0 {
                ez = c(0.0, 0.0);
            }
            if m_abs_round == 1 {
                c(0.0, -1.0) * dpsi_dr
            } else {
                c(0.0, 0.0)
            }
        };
        [er, ez, et]
    });

    let big = e
        .iter()
        .flat_map(|v| v.iter())
        .cloned()
        .fold(c(0.0, 0.0), |acc, z| if z.norm() > acc.norm() { z } else { acc });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |v: Vec<f64>| node_average(simplex, &v.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, order);
    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
    ];
    let mut cells = CellFields::new(simplex, order);
    for (k, name) in [(0, "Ex"), (2, "Ey"), (1, "Ez")] {
        let z: Vec<C> = e.iter().map(|v| v[k] * rot * s).collect();
        cells.push_complex(name, "a.u.", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }
    (fields, cells)
}
