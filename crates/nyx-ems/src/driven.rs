// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Driven frequency-domain multi-port solver (`Study::Driven`) with lumped sheet ports (`lumped_port`),
//! 2D boundary-mode waveguide ports (`waveport`), kinetic-inductance surface impedance (`surface_impedance`),
//! linear Josephson junctions (`lumped_jj`), and greedy Hermite Tangential Galerkin projection-based
//! reduced-order modeling (PROM).
//!
//! At each frequency $\omega = 2\pi f$ ($k_0 = \omega / c_0$), the port-terminated 3D curl-curl system is:
//! $$\tilde{\mathbf{A}}(\omega)\,\mathbf{X}(\omega) = \text{RHS}(\omega),$$
//! $$\tilde{\mathbf{A}}(\omega) = \mathbf{K}_{\text{tot}} + i k_0\big(\mathbf{C}_{\text{IBC}} + \mu_0 c_0 \mathbf{B}\mathbf{Z}_0^{-1}\mathbf{B}^T\big) + \sum_{p\in\text{wave}} i\beta_p(\omega)\mathbf{W}_p - k_0^2\mathbf{M}_{\text{tot}},$$
//! where column $q \in \{0..P-1\}$ of $\text{RHS}(\omega)$ corresponds to a unit $1\,\text{W}^{1/2}$ incident
//! power wave on port $q$: $\text{RHS}_q(\omega) = 2 i \omega \mu_0 Z_{0,q}(\omega)^{-1/2} \mathbf{b}_q$.
//!
//! When `prom_tol = Some(tol)`, each greedy anchor frequency $\omega_a$ factorizes $\tilde{\mathbf{A}}(\omega_a)$
//! **once** via `DirectSolver::new(&a_anchor, false)` and solves for both the field $\mathbf{X}(\omega_a)$ and its
//! frequency derivative $\mathbf{X}'(\omega_a)$ (plus 2nd derivative $\mathbf{X}''(\omega_a)$) on the same LU
//! factorization, accumulating an orthonormal real projection basis $\mathbf{Q} \in \mathbb{R}^{n_{\text{red}}\times r}$
//! and an isometric QR image-space basis $\mathbf{U} \in \mathbb{R}^{n_{\text{red}}\times s}$ for microsecond
//! full-space residual certification across all sweep frequencies.

use crate::cavity::{assemble_maxwell, basis3, Ned3};
use crate::ibc::{assemble_ibc_and_jj_3d, reduce_real_csr, reduce_real_vec, AssembledJjPort, SheetCellKind, MU0};
use crate::modes::C0;
use crate::ports::{
    assemble_lumped_ports_3d, assemble_waveports_3d, validate_port_contacts_3d, AssembledLumpedPort, AssembledWavePort,
};
use crate::problem::{Problem, Study};
use crate::rom::build_tls_meta;
use femcore::fem::{geometry, CurvedTet};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCsr, RCsr, C};
use femcore::linalg::driver::{thin_svd_c64, DirectSolver};
use femcore::pml::cell_stretches;
use nyx_physics::{cell_fields_of, Discretisation, DrivenPortData, Field, FieldFrame, Scalar, Solution, SolveResult};
use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::time::Instant;

type C3 = [[C; 3]; 3];

#[inline]
fn dot_real(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}

/// Equivalence classes of vertices under Floquet-periodic boundaries.
fn vertex_keys(mesh: &Discretisation, matches: &[femcore::periodic::Match]) -> Vec<u32> {
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
                let (lo, hi) = (rs.min(rt), rs.max(rt));
                parent[hi as usize] = lo;
            }
        }
    }
    (0..n as u32).map(|i| find(&mut parent, i)).collect()
}

/// Periodic Nedelec DOF relations between the two sides of a Floquet pair.
fn ned_relations(
    mesh: &Discretisation,
    ned: &Ned3,
    m: &femcore::periodic::Match,
    phase: f64,
) -> Result<Vec<(usize, usize, C)>, String> {
    let f = femcore::periodic::phase_factor(phase);
    let mut out = vec![];
    let mut done_edges = BTreeSet::new();
    for &fc in &m.src_facets {
        let v = mesh.facet(fc);
        let vs = [v[0], v[1], v[2]];
        let vt = [m.map[&vs[0]], m.map[&vs[1]], m.map[&vs[2]]];
        for (i, j) in [(0, 1), (1, 2), (0, 2)] {
            let key = (vs[i].min(vs[j]), vs[i].max(vs[j]));
            if !done_edges.insert(key) {
                continue;
            }
            let es = mesh.edge_index(vs[i], vs[j]);
            let et = mesh.edges().binary_search(&(vt[i].min(vt[j]), vt[i].max(vt[j]))).map_err(|_| {
                "an edge of the source boundary has no counterpart on the destination boundary".to_string()
            })?;
            if ned.order == 1 {
                out.push((et, es, f));
            } else {
                out.push((2 * et, 2 * es, f));
                out.push((2 * et + 1, 2 * es + 1, f));
            }
        }
        if ned.order == 2 {
            let ds = ned.face_dofs(vs);
            let dt = ned.face_dofs(vt);
            out.push((dt[0], ds[0], f));
            out.push((dt[1], ds[1], f));
        }
    }
    Ok(out)
}

/// Kind of driven port in the unified port list `0..P`.
#[derive(Clone, Debug)]
pub enum DrivenPortKind {
    Lumped { z0: f64 },
    Wave { wave_index: usize },
}

/// Unified metadata for one driven port $p \in \{0..P-1\}$.
#[derive(Clone, Debug)]
pub struct DrivenPortEntry {
    pub name: String,
    pub kind: DrivenPortKind,
    /// Equilibrated reduced port functional $\hat{\mathbf{b}}_p = \mathbf{D} \mathbf{b}_{p,\text{red}}$ (`nred`).
    pub b_hat: Vec<f64>,
}

impl DrivenPortEntry {
    /// Returns `(z0_p(omega), gamma_p(omega), d(gamma_p)/d(omega), d^2(gamma_p)/d(omega)^2)` at angular frequency `omega` (rad/s),
    /// where $\text{RHS}_p(\omega) = \gamma_p(\omega) \hat{\mathbf{b}}_p$ and $\gamma_p(\omega) = 2 i \omega \mu_0 Z_{0,p}(\omega)^{-1/2}$.
    #[inline]
    pub fn excitation_coeffs(&self, omega: f64, waveports: &[AssembledWavePort]) -> (f64, C, C, C) {
        match self.kind {
            DrivenPortKind::Lumped { z0 } => {
                let c_scale = 2.0 * MU0 * z0.max(1e-30).sqrt().recip();
                (z0, c(0.0, omega * c_scale), c(0.0, c_scale), c(0.0, 0.0))
            }
            DrivenPortKind::Wave { wave_index } => {
                let (beta, dbeta, d2beta, z0) = waveports[wave_index].beta_and_derivs(omega);
                // omega * mu_0 * Z_{0,p}^(-1/2) = sqrt(omega * mu_0 * beta)
                let g = (omega * MU0 * beta).max(1e-30);
                let dg = MU0 * (beta + omega * dbeta);
                let d2g = MU0 * (2.0 * dbeta + omega * d2beta);
                let sqrt_g = g.sqrt();
                let d_sqrt_g = 0.5 * dg / sqrt_g;
                let d2_sqrt_g = 0.5 * d2g / sqrt_g - 0.25 * dg * dg / (g * sqrt_g);
                (z0, c(0.0, 2.0 * sqrt_g), c(0.0, 2.0 * d_sqrt_g), c(0.0, 2.0 * d2_sqrt_g))
            }
        }
    }
    /// Returns the circuit reference impedance `Z_0` used for `DrivenPortData::z0_ohm` and `S -> Z/Y` conversion.
    /// For a `waveport` with an explicit `z0_ohm`, returns that circuit impedance while power-wave `excitation_coeffs`
    /// and `s_k` continue to use the modal wave impedance `z_wave = omega mu_0 / beta`.
    #[inline]
    pub fn circuit_z0(&self, z_wave: f64, waveports: &[AssembledWavePort]) -> f64 {
        match self.kind {
            DrivenPortKind::Lumped { z0 } => z0,
            DrivenPortKind::Wave { wave_index } => waveports[wave_index].bc.z0_ohm.unwrap_or(z_wave),
        }
    }
}

/// Stage 1 of the 3D driven frequency-domain solver: Nedelec numbering, PEC/periodic reduction,
/// Ruiz-equilibrated affine system matrices $(\hat{\mathbf{K}}_{\text{tot}}, \hat{\mathbf{C}}_{\text{tot}}, \hat{\mathbf{M}}_{\text{tot}}, \hat{\mathbf{W}}_w)$,
/// and assembled lumped/waveguide ports.
pub struct DrivenAssembly {
    pub order: usize,
    pub n: usize,
    pub nred: usize,
    pub ned: Ned3,
    pub red: femcore::periodic::Reduction,
    pub f_min: f64,
    pub f_max: f64,
    pub n_points: usize,
    pub freqs_hz: Option<Vec<f64>>,
    pub prom_tol: Option<f64>,
    pub max_hdm_solves: usize,
    pub strict_prom_tol: bool,
    pub seed_freqs_hz: Vec<f64>,
    pub d_ruiz: Vec<f64>,
    /// Affine operator terms $\hat{\mathbf{A}}_m = \mathbf{D} \mathbf{A}_m \mathbf{D}$ ($m = 0 \dots 3 + P_{\text{wave}} - 1$):
    /// `[K_tot, C_tot, M_tot, W_0, ..., W_{P_wave - 1}]`.
    pub terms_hat: Vec<CCsr>,
    pub ports: Vec<DrivenPortEntry>,
    pub lumped_ports: Vec<AssembledLumpedPort>,
    pub waveports: Vec<AssembledWavePort>,
    pub jj_ports: Vec<AssembledJjPort>,
    pub cell_sheet: Vec<Option<SheetCellKind>>,
    pub mu_inv: Vec<[C; 3]>,
    pub has_sheets: bool,
}

impl DrivenAssembly {
    /// Evaluates the complex scalar coefficients $(\alpha_m(\omega), \alpha_m'(\omega), \alpha_m''(\omega))$ for each term $\hat{\mathbf{A}}_m$.
    pub fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)> {
        let k0 = omega / C0;
        let dk0 = 1.0 / C0;
        let mut out = Vec::with_capacity(3 + self.waveports.len());
        // Term 0: K_tot * 1
        out.push((c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)));
        // Term 1: C_tot * (i k_0)
        out.push((c(0.0, k0), c(0.0, dk0), c(0.0, 0.0)));
        // Term 2: M_tot * (-k_0^2)
        out.push((c(-k0 * k0, 0.0), c(-2.0 * k0 * dk0, 0.0), c(-2.0 * dk0 * dk0, 0.0)));
        // Terms 3..: W_w * (i beta_w(omega))
        for wp in &self.waveports {
            let (beta, dbeta, d2beta, _) = wp.beta_and_derivs(omega);
            out.push((c(0.0, beta), c(0.0, dbeta), c(0.0, d2beta)));
        }
        out
    }

    /// Assembles the full equilibrated matrix $\hat{\mathbf{A}}(\omega) = \sum_m \alpha_m(\omega) \hat{\mathbf{A}}_m$
    /// on the fixed union CSR sparsity pattern of all non-empty affine terms (preserving explicit zeros so
    /// `DirectSolver::update_matrix` can reuse its symbolic LU factorization across every anchor frequency).
    pub fn assemble_matrix_at(&self, omega: f64) -> CCsr {
        let coeffs = self.term_coeffs(omega);
        let mut mat = self.terms_hat[0].scaled(coeffs[0].0);
        for (m, &(alpha, _, _)) in coeffs.iter().enumerate().skip(1) {
            if self.terms_hat[m].nnz() > 0 {
                mat = mat.add_scaled(&self.terms_hat[m], alpha);
            }
        }
        mat
    }

    /// Computes $\mathbf{y} = \hat{\mathbf{A}}'(\omega) \mathbf{x}$ or $\hat{\mathbf{A}}''(\omega) \mathbf{x}$.
    pub fn apply_matrix_deriv(&self, omega: f64, order_deriv: usize, x: &[C]) -> Vec<C> {
        let coeffs = self.term_coeffs(omega);
        let mut out = vec![c(0.0, 0.0); self.nred];
        for (m, &(_, d1, d2)) in coeffs.iter().enumerate() {
            let coeff = if order_deriv == 1 { d1 } else { d2 };
            if (coeff.re == 0.0 && coeff.im == 0.0) || self.terms_hat[m].nnz() == 0 {
                continue;
            }
            let mx = self.terms_hat[m].matvec(x);
            for (yi, &mxi) in out.iter_mut().zip(&mx) {
                *yi += coeff * mxi;
            }
        }
        out
    }
}

/// Assembles Stage 1 (`DrivenAssembly`) for a 3D `Study::Driven` problem.
pub fn assemble_emw_driven(mesh: &Discretisation, problem: &Problem) -> Result<DrivenAssembly, String> {
    if mesh.dim != 3 {
        return Err("the driven study of EM waves needs a 3D mesh".into());
    }
    let Study::Driven {
        f_min: study_f_min,
        f_max: study_f_max,
        n_points: study_n_points,
        freqs_hz,
        prom_tol,
        max_hdm_solves,
        strict_prom_tol,
        seed_freqs_hz,
    } = problem.study.clone()
    else {
        return Err("assemble_emw_driven requires a Driven study".into());
    };
    let (f_min, f_max, n_points) = if let Some(fs) = freqs_hz.as_ref().filter(|v| !v.is_empty()) {
        (fs[0], fs[fs.len() - 1], fs.len())
    } else {
        (study_f_min, study_f_max, study_n_points)
    };
    if problem.lumped_ports.is_empty() && problem.waveports.is_empty() {
        return Err("driven study requires at least one 'lumped_port' or 'waveport' boundary condition".into());
    }
    validate_port_contacts_3d(mesh, problem)?;
    let order = if problem.order == 1 { 1 } else { 2 };
    let tens = problem.materials.cell_complex_permittivity(mesh, Some(1.0))?;
    let eps_phys_re: Vec<[[f64; 3]; 3]> = tens.iter().map(|(re, _)| *re).collect();
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let has_pml = problem.materials.has_pml();
    let stretch = if has_pml { Some(cell_stretches(mesh, &problem.materials.pml)) } else { None };
    let mut eps: Vec<C3> = Vec::with_capacity(mesh.n_cells());
    let mut mu_inv: Vec<[C; 3]> = Vec::with_capacity(mesh.n_cells());
    for ci in 0..mesh.n_cells() {
        let (sx, sy, sz) =
            stretch.as_ref().map_or((c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1], s[ci][2]));
        let lam = [sy * sz / sx, sx * sz / sy, sx * sy / sz];
        let pml_here = stretch.is_some()
            && (sx.re != 1.0 || sx.im != 0.0 || sy.re != 1.0 || sy.im != 0.0 || sz.re != 1.0 || sz.im != 0.0);
        let (t_re, t_im) = tens[ci];
        let (sv, det) = ([sx, sy, sz], sx * sy * sz);
        let mut e = [[c(0.0, 0.0); 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let ec = c(t_re[a][b], -t_im[a][b]);
                e[a][b] = if pml_here { det * ec / (sv[a] * sv[b]) } else { ec };
            }
        }
        eps.push(e);
        let one = c(1.0, 0.0);
        let mu0 = c(1.0 / mur[ci], 0.0);
        mu_inv.push([
            mu0 / if pml_here { lam[0] } else { one },
            mu0 / if pml_here { lam[1] } else { one },
            mu0 / if pml_here { lam[2] } else { one },
        ]);
    }

    let mut matches = vec![];
    for p in &problem.periodic {
        matches.push(femcore::periodic::match_vertices(mesh, p)?);
    }
    let vkey = vertex_keys(mesh, &matches);
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let (k, m) = assemble_maxwell(mesh, &ned, &eps, &mu_inv);

    let f0 = 0.5 * (f_min + f_max);
    let has_sheets = !problem.pec.is_empty() || !problem.surface_impedance.is_empty() || !problem.lumped_jj.is_empty();
    let mut ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d(mesh, &ned, &problem.pec, &problem.surface_impedance, &problem.lumped_jj, f0)?)
    } else {
        None
    };

    let mut lumped_ports = assemble_lumped_ports_3d(mesh, &ned, &problem.pec, &problem.lumped_ports)?;
    let mut waveports = assemble_waveports_3d(
        mesh,
        &ned,
        &problem.pec,
        &problem.pmc,
        &problem.surface_impedance,
        &problem.waveports,
        &eps_phys_re,
        &mur,
        f_min,
        f_max,
    )?;

    // Exterior PEC walls exclude periodic, pmc, surface_impedance, lumped_jj, lumped_port, and waveport boundaries
    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let wp_labels: Vec<i32> = problem.waveports.iter().flat_map(|b| b.labels.iter().copied()).collect();

    let mut fixed_set = BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        if exterior[f]
            && !periodic_labels.contains(&label)
            && !problem.pmc.contains(&label)
            && !ibc_labels.contains(&label)
            && !jj_labels.contains(&label)
            && !lp_labels.contains(&label)
            && !wp_labels.contains(&label)
        {
            let fv = mesh.facet(f);
            for d in ned.boundary_dofs(mesh, [fv[0], fv[1], fv[2]]) {
                fixed_set.insert(d);
            }
        }
    }
    if let Some(ref ij) = ibc_jj {
        for &fv in &ij.pec_tris {
            for d in ned.boundary_dofs(mesh, fv) {
                fixed_set.insert(d);
            }
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let mut relations = vec![];
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        relations.extend(ned_relations(mesh, &ned, mt, p.phase)?);
    }
    let red = femcore::periodic::Reduction::build(n, &relations, &fixed);
    let pmat = red.prolong();
    let nred = red.n_reduced();
    if nred == 0 {
        return Err("driven study has no free degrees of freedom after PEC constraints".into());
    }

    let mut k_tot_full = k;
    let mut m_tot_full = m;
    let mut c_tot_r_full = RCsr::zeros(n, n);

    let (jj_ports, cell_sheet) = if let Some(ref mut ij) = ibc_jj {
        if ij.k_ibc_full.nnz() > 0 {
            k_tot_full = k_tot_full.add_scaled(&ij.k_ibc_full.to_complex(), c(1.0, 0.0));
        }
        if ij.k_jj_full.nnz() > 0 {
            k_tot_full = k_tot_full.add_scaled(&ij.k_jj_full.to_complex(), c(1.0, 0.0));
        }
        if ij.m_jj_full.nnz() > 0 {
            m_tot_full = m_tot_full.add_scaled(&ij.m_jj_full.to_complex(), c(1.0, 0.0));
        }
        if ij.c_ibc_full.nnz() > 0 {
            c_tot_r_full = c_tot_r_full.add_scaled(&ij.c_ibc_full, 1.0);
        }
        for port in &mut ij.jj_ports {
            port.b_red = reduce_real_vec(&red, &port.b_full);
            port.k_jj_red = reduce_real_csr(&red, &pmat, &port.k_jj_red);
            port.m_jj_red = reduce_real_csr(&red, &pmat, &port.m_jj_red);
        }
        (std::mem::take(&mut ij.jj_ports), std::mem::take(&mut ij.cell_sheet))
    } else {
        (Vec::new(), vec![None; mesh.n_cells()])
    };

    for lp in &mut lumped_ports {
        lp.b_red = reduce_real_vec(&red, &lp.b_full);
        if lp.c_port_full.nnz() > 0 {
            c_tot_r_full = c_tot_r_full.add_scaled(&lp.c_port_full, 1.0);
        }
    }
    for wp in &mut waveports {
        wp.b_red = reduce_real_vec(&red, &wp.b_full);
        wp.w_red = reduce_real_csr(&red, &pmat, &wp.w_full);
    }

    let k_tot_red = femcore::periodic::congruence(&pmat, &k_tot_full);
    let m_tot_red = femcore::periodic::congruence(&pmat, &m_tot_full);
    let c_tot_red = if c_tot_r_full.nnz() > 0 {
        femcore::periodic::congruence(&pmat, &c_tot_r_full.to_complex())
    } else {
        CCsr::zeros(nred, nred)
    };

    // Compute iterative symmetric Ruiz diagonal equilibration at midband omega_mid
    let omega_mid = 2.0 * PI * f0;
    let k0_mid = (omega_mid / C0).max(1.0);
    let k0_mid2 = k0_mid * k0_mid;
    let wp_betas_mid: Vec<f64> = waveports.iter().map(|wp| wp.beta_and_derivs(omega_mid).0).collect();

    let mut d_ruiz = vec![1.0_f64; nred];
    for _ in 0..8 {
        let mut rmax = vec![0.0_f64; nred];
        for i in 0..nred {
            let di = d_ruiz[i];
            let mut m_val = 0.0_f64;
            let (ck, vk) = k_tot_red.row(i);
            for (&j, &v) in ck.iter().zip(vk) {
                m_val = m_val.max(di * v.norm() * d_ruiz[j as usize]);
            }
            if c_tot_red.nnz() > 0 {
                let (cc, vc) = c_tot_red.row(i);
                for (&j, &v) in cc.iter().zip(vc) {
                    m_val = m_val.max(di * k0_mid * v.norm() * d_ruiz[j as usize]);
                }
            }
            let (cm, vm) = m_tot_red.row(i);
            for (&j, &v) in cm.iter().zip(vm) {
                m_val = m_val.max(di * k0_mid2 * v.norm() * d_ruiz[j as usize]);
            }
            for (w_idx, wp) in waveports.iter().enumerate() {
                if wp.w_red.nnz() > 0 {
                    let bw = wp_betas_mid[w_idx];
                    let (cw, vw) = wp.w_red.row(i);
                    for (&j, &v) in cw.iter().zip(vw) {
                        m_val = m_val.max(di * bw * v.abs() * d_ruiz[j as usize]);
                    }
                }
            }
            rmax[i] = m_val;
        }
        for i in 0..nred {
            if rmax[i] > 1e-300 {
                d_ruiz[i] = (d_ruiz[i] / rmax[i].sqrt()).clamp(1e-14, 1e14);
            }
        }
    }

    let scale_sym_c = |mat: &CCsr| -> CCsr {
        if mat.nnz() == 0 {
            return CCsr::zeros(nred, nred);
        }
        let mut out = mat.clone();
        for i in 0..nred {
            let di = d_ruiz[i];
            let range = out.row_range(i);
            let cols: Vec<u32> = out.indices()[range.clone()].to_vec();
            let vals = &mut out.values_mut()[range];
            for (v, &col) in vals.iter_mut().zip(&cols) {
                *v *= di * d_ruiz[col as usize];
            }
        }
        out
    };

    let mut terms_hat = Vec::with_capacity(3 + waveports.len());
    terms_hat.push(scale_sym_c(&k_tot_red));
    terms_hat.push(scale_sym_c(&c_tot_red));
    terms_hat.push(scale_sym_c(&m_tot_red));
    for wp in &waveports {
        terms_hat.push(scale_sym_c(&wp.w_red.to_complex()));
    }

    let mut ports = Vec::with_capacity(lumped_ports.len() + waveports.len());
    for lp in &lumped_ports {
        let b_hat: Vec<f64> = lp.b_red.iter().zip(&d_ruiz).map(|(&b, &d)| b * d).collect();
        ports.push(DrivenPortEntry { name: lp.bc.name.clone(), kind: DrivenPortKind::Lumped { z0: lp.bc.z0 }, b_hat });
    }
    for (wave_index, wp) in waveports.iter().enumerate() {
        let b_hat: Vec<f64> = wp.b_red.iter().zip(&d_ruiz).map(|(&b, &d)| b * d).collect();
        ports.push(DrivenPortEntry { name: wp.bc.name.clone(), kind: DrivenPortKind::Wave { wave_index }, b_hat });
    }

    Ok(DrivenAssembly {
        order,
        n,
        nred,
        ned,
        red,
        f_min,
        f_max,
        n_points,
        freqs_hz,
        prom_tol,
        max_hdm_solves,
        strict_prom_tol,
        seed_freqs_hz,
        d_ruiz,
        terms_hat,
        ports,
        lumped_ports,
        waveports,
        jj_ports,
        cell_sheet,
        mu_inv,
        has_sheets,
    })
}

/// Two-pass Modified Gram-Schmidt $l_2$-orthonormalization of `cand` into `basis`.
fn append_l2_orthonormal(basis: &mut Vec<Vec<f64>>, cand: &[f64]) -> bool {
    let norm0 = dot_real(cand, cand).sqrt();
    if norm0 <= 1e-30 {
        return false;
    }
    let mut v = cand.to_vec();
    for _pass in 0..2 {
        for q in basis.iter() {
            let proj = dot_real(&v, q);
            for (vi, &qi) in v.iter_mut().zip(q) {
                *vi -= proj * qi;
            }
        }
    }
    let norm = dot_real(&v, &v).sqrt();
    if norm <= 1e-12 * norm0 {
        return false;
    }
    let inv = 1.0 / norm;
    for vi in &mut v {
        *vi *= inv;
    }
    basis.push(v);
    true
}

/// Solves an $r \times r$ complex dense linear system $A X = B$ for $P$ right-hand sides (`b_cols[q]` of length $r$)
/// using Gaussian elimination with partial pivoting.
fn solve_dense_complex_multi(a_mat: &[C], b_cols: &[Vec<C>], r: usize) -> Option<Vec<Vec<C>>> {
    let p = b_cols.len();
    if r == 0 {
        return Some(vec![Vec::new(); p]);
    }
    let mut a = a_mat.to_vec();
    let mut rhs = vec![c(0.0, 0.0); r * p];
    for q in 0..p {
        for i in 0..r {
            rhs[i * p + q] = b_cols[q][i];
        }
    }

    for k in 0..r {
        let mut piv_row = k;
        let mut piv_max = a[k * r + k].norm();
        for i in (k + 1)..r {
            let val = a[i * r + k].norm();
            if val > piv_max {
                piv_max = val;
                piv_row = i;
            }
        }
        if piv_max <= 1e-300 {
            return None;
        }
        if piv_row != k {
            for j in k..r {
                a.swap(k * r + j, piv_row * r + j);
            }
            for q in 0..p {
                rhs.swap(k * p + q, piv_row * p + q);
            }
        }
        let inv_piv = c(1.0, 0.0) / a[k * r + k];
        for i in (k + 1)..r {
            let factor = a[i * r + k] * inv_piv;
            if factor.re == 0.0 && factor.im == 0.0 {
                continue;
            }
            for j in (k + 1)..r {
                let akj = a[k * r + j];
                a[i * r + j] -= factor * akj;
            }
            for q in 0..p {
                let rkq = rhs[k * p + q];
                rhs[i * p + q] -= factor * rkq;
            }
        }
    }

    let mut x_cols = vec![vec![c(0.0, 0.0); r]; p];
    for i in (0..r).rev() {
        let inv_diag = c(1.0, 0.0) / a[i * r + i];
        for q in 0..p {
            let mut sum = rhs[i * p + q];
            for j in (i + 1)..r {
                sum -= a[i * r + j] * x_cols[q][j];
            }
            x_cols[q][i] = sum * inv_diag;
        }
    }
    Some(x_cols)
}

/// Converts a $P \times P$ scattering matrix $\mathbf{S}$ and reference impedances $Z_{0,p}$ into
/// impedance $\mathbf{Z} = \mathbf{Z}_0^{1/2}(\mathbf{I} + \mathbf{S})(\mathbf{I} - \mathbf{S})^{-1}\mathbf{Z}_0^{1/2}$
/// and admittance $\mathbf{Y} = \mathbf{Z}_0^{-1/2}(\mathbf{I} - \mathbf{S})(\mathbf{I} + \mathbf{S})^{-1}\mathbf{Z}_0^{-1/2}$.
fn s_to_z_and_y(s_mat: &[Vec<C>], z0_ports: &[f64]) -> (Vec<Vec<[f64; 2]>>, Vec<Vec<[f64; 2]>>, Vec<Vec<[f64; 2]>>) {
    let p = s_mat.len();
    let sqrt_z0: Vec<f64> = z0_ports.iter().map(|&z| z.max(1e-30).sqrt()).collect();
    let inv_sqrt_z0: Vec<f64> = sqrt_z0.iter().map(|&sz| 1.0 / sz).collect();

    // Form V_mat = Z_0^{1/2} (I + S) and I_mat = Z_0^{-1/2} (I - S)
    // Since V = Z * I_mat => Z^T = (I_mat^T)^{-1} V^T, and I_mat = Y * V => Y^T = (V^T)^{-1} I_mat^T
    let mut v_t = vec![c(0.0, 0.0); p * p];
    let mut i_t = vec![c(0.0, 0.0); p * p];
    let mut v_cols_t = vec![vec![c(0.0, 0.0); p]; p];
    let mut i_cols_t = vec![vec![c(0.0, 0.0); p]; p];

    for r in 0..p {
        for col in 0..p {
            let delta = if r == col { c(1.0, 0.0) } else { c(0.0, 0.0) };
            let v_rc = (delta + s_mat[r][col]) * sqrt_z0[r];
            let i_rc = (delta - s_mat[r][col]) * inv_sqrt_z0[r];
            // Transpose: entry (col, r)
            v_t[col * p + r] = v_rc;
            i_t[col * p + r] = i_rc;
            v_cols_t[r][col] = v_rc;
            i_cols_t[r][col] = i_rc;
        }
    }
    // Slight diagonal regularization only if singular
    let solve_with_reg = |mat_t: &[C], rhs_cols_t: &[Vec<C>]| -> Vec<Vec<C>> {
        if let Some(sol) = solve_dense_complex_multi(mat_t, rhs_cols_t, p) {
            return sol;
        }
        let mut reg = mat_t.to_vec();
        for d in 0..p {
            reg[d * p + d] += c(1e-14, 0.0);
        }
        solve_dense_complex_multi(&reg, rhs_cols_t, p).unwrap_or_else(|| vec![vec![c(0.0, 0.0); p]; p])
    };

    let z_rows = solve_with_reg(&i_t, &v_cols_t);
    let y_rows = solve_with_reg(&v_t, &i_cols_t);

    let mut s_out = vec![vec![[0.0; 2]; p]; p];
    let mut z_out = vec![vec![[0.0; 2]; p]; p];
    let mut y_out = vec![vec![[0.0; 2]; p]; p];
    for r in 0..p {
        for col in 0..p {
            s_out[r][col] = [s_mat[r][col].re, s_mat[r][col].im];
            z_out[r][col] = [z_rows[r][col].re, z_rows[r][col].im];
            y_out[r][col] = [y_rows[r][col].re, y_rows[r][col].im];
        }
    }
    (s_out, z_out, y_out)
}

/// Real Hermite Tangential Galerkin PROM with isometric QR image-space basis $\mathbf{U}$ for exact
/// microsecond full-space residual certification across all sweep frequencies.
struct DrivenProm<'a> {
    asm: &'a DrivenAssembly,
    /// Persistent sparse LU solver reusing its cached symbolic factorization across all anchor frequencies.
    solver: Option<DirectSolver<C>>,
    /// Orthonormal real basis columns $\mathbf{q}_i \in \mathbb{R}^{n_{\text{red}}}$ (`r` columns).
    q_cols: Vec<Vec<f64>>,
    /// Reduced operator matrices $\mathbf{A}_{m,r} = \mathbf{Q}^T \hat{\mathbf{A}}_m \mathbf{Q}$ (`M` matrices, each `r x r` row-major).
    a_r: Vec<Vec<C>>,
    /// Reduced port vectors $\mathbf{B}_r = \mathbf{Q}^T \hat{\mathbf{B}}$ (`r x P`, indexed as `b_r[i][q]`).
    b_r: Vec<Vec<f64>>,
    /// Orthonormal image-space basis $\mathbf{U} = [\mathbf{u}_0, \dots, \mathbf{u}_{s-1}] \in \mathbb{R}^{n_{\text{red}} \times s}$.
    u_cols: Vec<Vec<f64>>,
    /// Coordinates of port vectors $\hat{\mathbf{b}}_q$ in $\mathbf{U}$: `u_b[q]` has length `s`.
    u_b: Vec<Vec<f64>>,
    /// Squared norms $\|\hat{\mathbf{b}}_q\|_2^2$.
    b_norm2: Vec<f64>,
    /// Coordinates of image vectors $\hat{\mathbf{A}}_m \mathbf{q}_i$ in $\mathbf{U}$: `u_aq[m][i]` has length `s`.
    u_aq: Vec<Vec<Vec<C>>>,
}

impl<'a> DrivenProm<'a> {
    fn new(asm: &'a DrivenAssembly) -> Self {
        let m_terms = asm.terms_hat.len();
        let p_ports = asm.ports.len();
        let mut u_cols: Vec<Vec<f64>> = Vec::new();
        for port in &asm.ports {
            append_l2_orthonormal(&mut u_cols, &port.b_hat);
        }
        let s = u_cols.len();
        let mut u_b = Vec::with_capacity(p_ports);
        let mut b_norm2 = Vec::with_capacity(p_ports);
        for port in &asm.ports {
            b_norm2.push(dot_real(&port.b_hat, &port.b_hat));
            let coords: Vec<f64> = u_cols.iter().map(|u| dot_real(u, &port.b_hat)).collect();
            debug_assert_eq!(coords.len(), s);
            u_b.push(coords);
        }
        Self {
            asm,
            solver: None,
            q_cols: Vec::new(),
            a_r: vec![Vec::new(); m_terms],
            b_r: Vec::new(),
            u_cols,
            u_b,
            b_norm2,
            u_aq: vec![Vec::new(); m_terms],
        }
    }

    fn dim(&self) -> usize {
        self.q_cols.len()
    }

    fn symbolic_count(&self) -> usize {
        self.solver.as_ref().map_or(0, |s| s.stats().symbolic_count)
    }

    /// Appends a real candidate vector to $\mathbf{Q}$ (if not already in the span) and updates
    /// $\mathbf{A}_{m,r}$, $\mathbf{B}_r$, and the isometric image-space coordinates in $\mathbf{U}$.
    fn add_real_vector(&mut self, cand: &[f64]) -> bool {
        if !append_l2_orthonormal(&mut self.q_cols, cand) {
            return false;
        }
        let r_old = self.q_cols.len() - 1;
        let r_new = r_old + 1;
        let q_new = &self.q_cols[r_old];
        let m_terms = self.asm.terms_hat.len();
        let nred = self.asm.nred;

        // 1. Update B_r = Q^T B_hat
        let b_row: Vec<f64> = self.asm.ports.iter().map(|p| dot_real(q_new, &p.b_hat)).collect();
        self.b_r.push(b_row);

        // 2. Compute images y_m = A_m q_new, update A_{m,r} and image basis U
        let q_new_c: Vec<C> = q_new.iter().map(|&x| c(x, 0.0)).collect();
        for m in 0..m_terms {
            let ym = if self.asm.terms_hat[m].nnz() > 0 {
                self.asm.terms_hat[m].matvec(&q_new_c)
            } else {
                vec![c(0.0, 0.0); nred]
            };

            // Grow A_{m,r} from (r_old x r_old) to (r_new x r_new)
            let mut grown = vec![c(0.0, 0.0); r_new * r_new];
            for i in 0..r_old {
                for j in 0..r_old {
                    grown[i * r_new + j] = self.a_r[m][i * r_old + j];
                }
            }
            for i in 0..r_new {
                let qi = &self.q_cols[i];
                let mut dot_val = c(0.0, 0.0);
                for (&qik, &ymk) in qi.iter().zip(&ym) {
                    dot_val += ymk * qik;
                }
                // Since A_m is complex symmetric (A_m^T = A_m) and Q is real, A_{m,r} is complex symmetric
                grown[i * r_new + r_old] = dot_val;
                grown[r_old * r_new + i] = dot_val;
            }
            self.a_r[m] = grown;

            // Expand image basis U with Re(ym) and Im(ym)
            let re_y: Vec<f64> = ym.iter().map(|z| z.re).collect();
            let im_y: Vec<f64> = ym.iter().map(|z| z.im).collect();
            append_l2_orthonormal(&mut self.u_cols, &re_y);
            append_l2_orthonormal(&mut self.u_cols, &im_y);

            let s_now = self.u_cols.len();
            let mut coords = Vec::with_capacity(s_now);
            for u in &self.u_cols {
                let mut d_re = 0.0;
                let mut d_im = 0.0;
                for (k, &uk) in u.iter().enumerate() {
                    d_re += uk * ym[k].re;
                    d_im += uk * ym[k].im;
                }
                coords.push(c(d_re, d_im));
            }
            self.u_aq[m].push(coords);
        }

        // Pad earlier u_b and u_aq vectors with zeros up to the current image-space dimension s_now
        let s_now = self.u_cols.len();
        for ub in &mut self.u_b {
            ub.resize(s_now, 0.0);
        }
        for m in 0..m_terms {
            for col in &mut self.u_aq[m] {
                col.resize(s_now, c(0.0, 0.0));
            }
        }
        true
    }

    /// Splits a complex snapshot vector `z` into `Re(z)` and `Im(z)` (after phase rotation to align its dominant part)
    /// and appends both to $\mathbf{Q}$.
    fn add_complex_snapshot(&mut self, z: &[C]) {
        let big = z.iter().copied().fold(c(0.0, 0.0), |m, v| if v.norm_sqr() > m.norm_sqr() { v } else { m });
        if big.norm_sqr() <= 1e-300 {
            return;
        }
        let rot = big.conj() / big.norm();
        let re_v: Vec<f64> = z.iter().map(|&v| (v * rot).re).collect();
        let im_v: Vec<f64> = z.iter().map(|&v| (v * rot).im).collect();
        self.add_real_vector(&re_v);
        self.add_real_vector(&im_v);
    }

    /// Performs one full HDM numeric factorization at angular frequency `omega_a` (rad/s), reusing the
    /// cached symbolic LU analysis, normalizes each of the `3P` Hermite snapshots (`Y`, `Y'`, `Y''`) directly
    /// to unit $L_2$ norm, compresses the anchor block via complex thin SVD (`s_k >= 1e-6 * s_0`), and adds
    /// the surviving left singular vectors to $\mathbf{Q}$.
    fn add_hermite_anchor(&mut self, omega_a: f64) -> Result<(), String> {
        let a_mat = self.asm.assemble_matrix_at(omega_a);
        match self.solver.as_mut() {
            Some(solver) => solver.update_matrix(&a_mat)?,
            None => self.solver = Some(DirectSolver::new(&a_mat, false)?),
        }
        let p_ports = self.asm.ports.len();
        let nred = self.asm.nred;

        let mut rhs_0 = Vec::with_capacity(p_ports);
        let mut gammas = Vec::with_capacity(p_ports);
        for port in &self.asm.ports {
            let (_, g0, g1, g2) = port.excitation_coeffs(omega_a, &self.asm.waveports);
            gammas.push((g0, g1, g2));
            rhs_0.push(port.b_hat.iter().map(|&b| g0 * b).collect::<Vec<C>>());
        }

        // 1. Value solve: A(omega_a) Y = RHS(omega_a)
        let y0_cols = self.solver.as_ref().unwrap().solve_many(&rhs_0);

        // 2. First derivative solve on the SAME LU factorization:
        //    A(omega_a) Y' = RHS'(omega_a) - A'(omega_a) Y
        let mut rhs_1 = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let a1_y0 = self.asm.apply_matrix_deriv(omega_a, 1, &y0_cols[q]);
            let g1 = gammas[q].1;
            let b_hat = &self.asm.ports[q].b_hat;
            let col: Vec<C> = (0..nred).map(|i| g1 * b_hat[i] - a1_y0[i]).collect();
            rhs_1.push(col);
        }
        let y1_cols = self.solver.as_ref().unwrap().solve_many(&rhs_1);

        // 3. Second derivative solve on the SAME LU factorization:
        //    A(omega_a) Y'' = RHS''(omega_a) - 2 A'(omega_a) Y' - A''(omega_a) Y
        let mut rhs_2 = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let a1_y1 = self.asm.apply_matrix_deriv(omega_a, 1, &y1_cols[q]);
            let a2_y0 = self.asm.apply_matrix_deriv(omega_a, 2, &y0_cols[q]);
            let g2 = gammas[q].2;
            let b_hat = &self.asm.ports[q].b_hat;
            let col: Vec<C> = (0..nred).map(|i| g2 * b_hat[i] - c(2.0, 0.0) * a1_y1[i] - a2_y0[i]).collect();
            rhs_2.push(col);
        }
        let y2_cols = self.solver.as_ref().unwrap().solve_many(&rhs_2);

        // Normalize every Y, Y', Y'' snapshot directly to unit L2 norm (bandwidth-independent),
        // then compress the 3P complex block via thin SVD keeping s_k >= 1e-6 * s_0.
        let mut anchor_cols: Vec<Vec<C>> = Vec::with_capacity(3 * p_ports);
        for col in y0_cols.iter().chain(y1_cols.iter()).chain(y2_cols.iter()) {
            let nrm = col.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
            if nrm > 0.0 && nrm.is_finite() {
                let inv = 1.0 / nrm;
                anchor_cols.push(col.iter().map(|&v| v * inv).collect());
            }
        }
        if !anchor_cols.is_empty() {
            let (s_vals, u_cols) = thin_svd_c64(&anchor_cols)?;
            if let Some(&s0) = s_vals.first() {
                if s0 > 0.0 {
                    let cutoff = 1e-6 * s0;
                    for (&sk, uk) in s_vals.iter().zip(&u_cols) {
                        if sk >= cutoff {
                            self.add_complex_snapshot(uk);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Solves the $r \times r$ reduced system at `omega` (rad/s) for all $P$ ports and returns
    /// `(x_r_cols, exact_relative_residual)`.
    fn solve_reduced_and_residual(&self, omega: f64) -> Option<(Vec<Vec<C>>, f64)> {
        let r = self.dim();
        let p_ports = self.asm.ports.len();
        let m_terms = self.asm.terms_hat.len();
        let coeffs = self.asm.term_coeffs(omega);

        let mut t_r = vec![c(0.0, 0.0); r * r];
        for m in 0..m_terms {
            let alpha = coeffs[m].0;
            if alpha.re == 0.0 && alpha.im == 0.0 {
                continue;
            }
            let am = &self.a_r[m];
            for idx in 0..(r * r) {
                t_r[idx] += alpha * am[idx];
            }
        }

        let mut gammas = Vec::with_capacity(p_ports);
        let mut rhs_r = Vec::with_capacity(p_ports);
        for q in 0..p_ports {
            let (_, g0, _, _) = self.asm.ports[q].excitation_coeffs(omega, &self.asm.waveports);
            gammas.push(g0);
            let col: Vec<C> = (0..r).map(|i| g0 * self.b_r[i][q]).collect();
            rhs_r.push(col);
        }

        let xr_cols = solve_dense_complex_multi(&t_r, &rhs_r, r)?;

        // Compute exact full-space relative residual in the isometric image basis U (dim s)
        let s = self.u_cols.len();
        let mut res_num2 = 0.0;
        let mut res_den2 = 0.0;
        let mut r_vec = vec![c(0.0, 0.0); s];

        for q in 0..p_ports {
            let g0 = gammas[q];
            res_den2 += g0.norm_sqr() * self.b_norm2[q];
            for (k, &ub_k) in self.u_b[q].iter().enumerate() {
                r_vec[k] = -g0 * ub_k;
            }
            let xrq = &xr_cols[q];
            for m in 0..m_terms {
                let alpha = coeffs[m].0;
                if alpha.re == 0.0 && alpha.im == 0.0 {
                    continue;
                }
                let u_aq_m = &self.u_aq[m];
                for i in 0..r {
                    let weight = alpha * xrq[i];
                    if weight.re == 0.0 && weight.im == 0.0 {
                        continue;
                    }
                    let col_u = &u_aq_m[i];
                    for k in 0..s {
                        r_vec[k] += weight * col_u[k];
                    }
                }
            }
            for &rv in &r_vec {
                res_num2 += rv.norm_sqr();
            }
        }

        let rel_res = (res_num2 / res_den2.max(1e-300)).sqrt();
        Some((xr_cols, rel_res))
    }

    /// Reconstructs the unscaled reduced-DOF physical electric field vector $\mathbf{x} = \mathbf{D} \mathbf{Q} \mathbf{x}_r$.
    fn reconstruct_unscaled(&self, xr: &[C]) -> Vec<C> {
        let nred = self.asm.nred;
        let mut y = vec![c(0.0, 0.0); nred];
        for (i, &c_val) in xr.iter().enumerate() {
            let qi = &self.q_cols[i];
            for k in 0..nred {
                y[k] += c_val * qi[k];
            }
        }
        for (yk, &dk) in y.iter_mut().zip(&self.asm.d_ruiz) {
            *yk *= dk;
        }
        y
    }
}

/// Builds a `Solution` (`|E|, Re Ex, Re Ey, Re Ez` nodal fields + `Ex, Ey, Ez` `"V/m"` and `Hx, Hy, Hz` `"A/m"` cell fields)
/// from an unscaled reduced-DOF solution vector `x_red` at frequency `freq_hz`.
fn build_driven_solution(
    mesh: &Discretisation,
    asm: &DrivenAssembly,
    freq_hz: f64,
    x_red: &[C],
    label: String,
) -> Solution {
    let x = asm.red.expand(x_red);
    let geo = geometry(mesh);
    let lamn = local_bary(3, 2);
    let nloc = lamn.len();
    let nc = mesh.n_cells();
    let nb = asm.ned.nb();

    let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    let mut curl_e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];

    for ci in 0..nc {
        let gid = asm.ned.keys(mesh, ci);
        for b in 0..4 {
            for kk in 0..3 {
                gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
            }
        }
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        for (a, l) in lamn.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            match &curved {
                Some(ct) => basis3(asm.order, &lam, &ct.jacobian(&lam).0, &gid, &mut w, &mut cu),
                None => basis3(asm.order, &lam, &gs, &gid, &mut w, &mut cu),
            }
            let mut val = [c(0.0, 0.0); 3];
            let mut cval = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let coef = x[asm.ned.cell_dofs[ci][s]];
                for kk in 0..3 {
                    val[kk] += coef * w[s][kk];
                    cval[kk] += coef * cu[s][kk];
                }
            }
            e[ci * nloc + a] = val;
            curl_e[ci * nloc + a] = cval;
        }
    }

    let big =
        e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, 2);

    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
    ];

    let mut cells = CellFields::new(mesh, 2);
    for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
        let z: Vec<C> = e.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "V/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }

    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let h_prefactor = c(0.0, 1.0 / (omega * MU0));
    let mut h_field = vec![[c(0.0, 0.0); 3]; nc * nloc];
    for ci in 0..nc {
        let mui = asm.mu_inv[ci];
        for a in 0..nloc {
            let idx = ci * nloc + a;
            let ce = curl_e[idx];
            h_field[idx] = [ce[0] * mui[0] * h_prefactor, ce[1] * mui[1] * h_prefactor, ce[2] * mui[2] * h_prefactor];
        }
    }
    for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
        let z: Vec<C> = h_field.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }

    Solution { label, fields, cell_fields: cell_fields_of(&cells), value: Some(freq_hz), eigenvalue: None }
}

/// Solves a 3D `Study::Driven` frequency sweep via Hermite Tangential Galerkin PROM (`prom_tol = Some(tol)`)
/// or direct point-by-point HDM factorizations (`prom_tol = None`).
pub fn solve_driven(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_emw_driven(mesh, problem)?;
    let freqs_hz: Vec<f64> = if let Some(fs) = asm.freqs_hz.as_ref().filter(|v| !v.is_empty()) {
        fs.to_vec()
    } else {
        let n_pts = asm.n_points.max(1);
        if n_pts == 1 {
            vec![asm.f_min]
        } else {
            (0..n_pts).map(|k| asm.f_min + (asm.f_max - asm.f_min) * (k as f64) / ((n_pts - 1) as f64)).collect()
        }
    };
    let n_points = freqs_hz.len();
    let omegas: Vec<f64> = freqs_hz.iter().map(|&f| 2.0 * PI * f).collect();
    let p_ports = asm.ports.len();

    let mut z0_ohm = Vec::with_capacity(n_points);
    let mut s_matrix = Vec::with_capacity(n_points);
    let mut z_matrix = Vec::with_capacity(n_points);
    let mut y_matrix = Vec::with_capacity(n_points);

    let mut hdm_solves = 0usize;
    let symbolic_count;
    let mut rom_dim = 0usize;
    let mut anchor_freqs_hz = Vec::new();
    let mut max_residual = 0.0_f64;
    let mut rom_sweep_seconds = 0.0_f64;
    let mut rep_freq_idx = n_points / 2;
    let mut rep_x_red: Vec<C> = vec![c(0.0, 0.0); asm.nred];

    if let Some(tol) = asm.prom_tol {
        let mut prom = DrivenProm::new(&asm);
        let max_solves = asm.max_hdm_solves.max(1);

        // Initial greedy anchors: [f_min, f_mid, f_max] merged with optional seed_freqs_hz
        let mut base_indices: Vec<usize> = if n_points == 1 || (asm.f_max - asm.f_min).abs() <= 1e-12 * asm.f_min {
            vec![0]
        } else if n_points == 2 {
            vec![0, 1]
        } else {
            vec![0, (n_points - 1) / 2, n_points - 1]
        };
        base_indices.dedup();

        let mut initial_anchors: Vec<(Option<usize>, f64)> =
            base_indices.iter().map(|&idx| (Some(idx), freqs_hz[idx])).collect();
        for &sf in &asm.seed_freqs_hz {
            let mut best_k = 0usize;
            let mut best_diff = f64::INFINITY;
            for (k, &fk) in freqs_hz.iter().enumerate() {
                let diff = (fk - sf).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_k = k;
                }
            }
            if best_diff <= 1e-9 * sf.max(1.0) {
                if !initial_anchors.iter().any(|&(oi, _)| oi == Some(best_k)) {
                    initial_anchors.push((Some(best_k), freqs_hz[best_k]));
                }
            } else if !initial_anchors.iter().any(|&(_, f)| (f - sf).abs() <= 1e-9 * sf.max(1.0)) {
                initial_anchors.push((None, sf));
            }
        }

        let mut anchor_indices: Vec<usize> = Vec::new();
        for (opt_idx, f_anchor) in initial_anchors.into_iter().take(max_solves) {
            femcore::progress::check()?;
            prom.add_hermite_anchor(2.0 * PI * f_anchor)?;
            hdm_solves += 1;
            anchor_freqs_hz.push(f_anchor);
            if let Some(idx) = opt_idx {
                anchor_indices.push(idx);
            }
        }

        const SWEEP_CHUNK: usize = 32;
        let mut xr_all: Vec<Vec<Vec<C>>>;

        loop {
            femcore::progress::check()?;
            let t_sweep = Instant::now();
            let dim_now = prom.dim();
            let prom_ref = &prom;
            let sweep_chunks = femcore::par::map_chunks(&omegas, SWEEP_CHUNK, |chunk| {
                chunk
                    .iter()
                    .map(|&w| match prom_ref.solve_reduced_and_residual(w) {
                        Some(pair) => pair,
                        None => (vec![vec![c(0.0, 0.0); dim_now]; p_ports], f64::INFINITY),
                    })
                    .collect::<Vec<_>>()
            });
            xr_all = Vec::with_capacity(n_points);
            let mut worst_res = 0.0_f64;
            let mut worst_cand_res = -1.0_f64;
            let mut worst_cand_idx = None;
            for (k, (xr_k, res_k)) in sweep_chunks.into_iter().flatten().enumerate() {
                if res_k > worst_res {
                    worst_res = res_k;
                }
                if !anchor_indices.contains(&k) && res_k > worst_cand_res {
                    worst_cand_res = res_k;
                    worst_cand_idx = Some(k);
                }
                xr_all.push(xr_k);
            }
            rom_sweep_seconds = t_sweep.elapsed().as_secs_f64();
            max_residual = worst_res;

            if worst_res <= tol || hdm_solves >= max_solves {
                break;
            }
            let Some(next_idx) = worst_cand_idx else {
                break;
            };
            anchor_indices.push(next_idx);
            prom.add_hermite_anchor(omegas[next_idx])?;
            hdm_solves += 1;
            anchor_freqs_hz.push(freqs_hz[next_idx]);
        }

        if asm.strict_prom_tol && max_residual > tol {
            return Err(format!(
                "PROM did not converge to prom_tol = {tol:.3e} within max_hdm_solves = {max_solves} (max_prom_residual = {max_residual:.3e})"
            ));
        }

        rom_dim = prom.dim();
        symbolic_count = prom.symbolic_count();

        // Extract S, Z, Y across all n_points from reduced coordinates X_r(omega_k) using deterministic fixed-chunk Rayon
        let t_post = Instant::now();
        let indices: Vec<usize> = (0..n_points).collect();
        let asm_ref = &asm;
        let prom_ref = &prom;
        let xr_ref = &xr_all;
        let omegas_ref = &omegas;
        let post_chunks = femcore::par::map_chunks(&indices, SWEEP_CHUNK, |chunk| {
            chunk
                .iter()
                .map(|&k| {
                    let w = omegas_ref[k];
                    let xr_k = &xr_ref[k];
                    let mut z_wave_k = Vec::with_capacity(p_ports);
                    let mut z_circ_k = Vec::with_capacity(p_ports);
                    let mut z0_pairs_k = Vec::with_capacity(p_ports);
                    for port in &asm_ref.ports {
                        let (z_wave_p, _, _, _) = port.excitation_coeffs(w, &asm_ref.waveports);
                        let z_circ_p = port.circuit_z0(z_wave_p, &asm_ref.waveports);
                        z_wave_k.push(z_wave_p);
                        z_circ_k.push(z_circ_p);
                        z0_pairs_k.push([z_circ_p, 0.0]);
                    }

                    let mut s_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
                    for p in 0..p_ports {
                        let inv_sqrt_z_wave_p = 1.0 / z_wave_k[p].max(1e-30).sqrt();
                        for q in 0..p_ports {
                            let mut v_pq = c(0.0, 0.0);
                            for i in 0..rom_dim {
                                v_pq += xr_k[q][i] * prom_ref.b_r[i][p];
                            }
                            let delta = if p == q { c(1.0, 0.0) } else { c(0.0, 0.0) };
                            s_k[p][q] = v_pq * inv_sqrt_z_wave_p - delta;
                        }
                    }
                    let metric = if p_ports >= 2 { s_k[1][0].norm() } else { s_k[0][0].norm() };
                    let (s_out, z_out, y_out) = s_to_z_and_y(&s_k, &z_circ_k);
                    (z0_pairs_k, s_out, z_out, y_out, metric)
                })
                .collect::<Vec<_>>()
        });
        let mut min_s_metric = f64::INFINITY;
        for (k, (z0_pairs_k, s_out, z_out, y_out, metric)) in post_chunks.into_iter().flatten().enumerate() {
            if metric < min_s_metric {
                min_s_metric = metric;
                rep_freq_idx = k;
            }
            z0_ohm.push(z0_pairs_k);
            s_matrix.push(s_out);
            z_matrix.push(z_out);
            y_matrix.push(y_out);
        }
        rom_sweep_seconds += t_post.elapsed().as_secs_f64();
        rep_x_red = prom.reconstruct_unscaled(&xr_all[rep_freq_idx][0]);
    } else {
        // Direct point-by-point HDM solve (prom_tol = None), reusing symbolic LU factorization across all frequencies
        let mut solver: Option<DirectSolver<C>> = None;
        let mut min_s_metric = f64::INFINITY;
        for (k, &w) in omegas.iter().enumerate() {
            femcore::progress::check()?;
            let a_mat = asm.assemble_matrix_at(w);
            match solver.as_mut() {
                Some(s) => s.update_matrix(&a_mat)?,
                None => solver = Some(DirectSolver::new(&a_mat, false)?),
            }
            let s_ref = solver.as_ref().unwrap();
            hdm_solves += 1;
            anchor_freqs_hz.push(freqs_hz[k]);

            let mut z_wave_k = Vec::with_capacity(p_ports);
            let mut z_circ_k = Vec::with_capacity(p_ports);
            let mut z0_pairs_k = Vec::with_capacity(p_ports);
            let mut rhs_cols = Vec::with_capacity(p_ports);
            for port in &asm.ports {
                let (z_wave_p, g0, _, _) = port.excitation_coeffs(w, &asm.waveports);
                let z_circ_p = port.circuit_z0(z_wave_p, &asm.waveports);
                z_wave_k.push(z_wave_p);
                z_circ_k.push(z_circ_p);
                z0_pairs_k.push([z_circ_p, 0.0]);
                rhs_cols.push(port.b_hat.iter().map(|&b| g0 * b).collect::<Vec<C>>());
            }
            z0_ohm.push(z0_pairs_k);

            let y_cols = s_ref.solve_many(&rhs_cols);
            let mut s_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
            for p in 0..p_ports {
                let inv_sqrt_z_wave_p = 1.0 / z_wave_k[p].max(1e-30).sqrt();
                let b_hat_p = &asm.ports[p].b_hat;
                for q in 0..p_ports {
                    let mut v_pq = c(0.0, 0.0);
                    for (&bp, &yq) in b_hat_p.iter().zip(&y_cols[q]) {
                        v_pq += yq * bp;
                    }
                    let delta = if p == q { c(1.0, 0.0) } else { c(0.0, 0.0) };
                    s_k[p][q] = v_pq * inv_sqrt_z_wave_p - delta;
                }
            }
            let metric = if p_ports >= 2 { s_k[1][0].norm() } else { s_k[0][0].norm() };
            if metric < min_s_metric {
                min_s_metric = metric;
                rep_freq_idx = k;
                rep_x_red = y_cols[0].iter().zip(&asm.d_ruiz).map(|(&yi, &di)| yi * di).collect();
            }
            let (s_out, z_out, y_out) = s_to_z_and_y(&s_k, &z_circ_k);
            s_matrix.push(s_out);
            z_matrix.push(z_out);
            y_matrix.push(y_out);
        }
        symbolic_count = solver.as_ref().map_or(0, |s| s.stats().symbolic_count);
    }

    let f_res_hz = freqs_hz[rep_freq_idx];
    let s_dip_mag = if p_ports >= 2 {
        let s21 = s_matrix[rep_freq_idx][1][0];
        s21[0].hypot(s21[1])
    } else {
        let s11 = s_matrix[rep_freq_idx][0][0];
        s11[0].hypot(s11[1])
    };
    let s21_min_db = 20.0 * s_dip_mag.max(1e-30).log10();

    let port_names: Vec<String> = asm.ports.iter().map(|p| p.name.clone()).collect();
    let driven_data = DrivenPortData {
        freqs_hz,
        ports: port_names,
        z0_ohm,
        s_matrix,
        z_matrix,
        y_matrix,
        hdm_solves,
        rom_dim,
        anchor_freqs_hz,
        max_residual,
    };

    let rep_sol = build_driven_solution(
        mesh,
        &asm,
        f_res_hz,
        &rep_x_red,
        format!("driven ({}): {:.6} GHz", asm.ports[0].name, f_res_hz / 1e9),
    );

    let (axes, left_handed) = femcore::section::axes_of_3d_periodic(&problem.periodic, mesh);
    let mut scalars = vec![
        driven_data.scalar(),
        ("hdm_solves".into(), Scalar::Num(hdm_solves as f64)),
        ("symbolic_count".into(), Scalar::Num(symbolic_count as f64)),
        ("rom_dim".into(), Scalar::Num(rom_dim as f64)),
        ("max_prom_residual".into(), Scalar::Num(max_residual)),
        ("rom_sweep_seconds".into(), Scalar::Num(rom_sweep_seconds)),
        ("S21_min_dB".into(), Scalar::Num(s21_min_db)),
        ("f_res_Hz".into(), Scalar::Num(f_res_hz)),
        FieldFrame { lateral: axes[0], normal: axes[1], propagation: axes[2], left_handed, outside: vec![] }.scalar(),
    ];
    if asm.has_sheets {
        scalars.push(build_tls_meta(mesh, problem).scalar());
    }

    let log = vec![format!(
        "driven sweep [{:.4}..{:.4} GHz, {} points, {} port(s)]: {} HDM solve(s), ROM dim {}, max residual {:.2e}, ROM sweep {:.4} s",
        asm.f_min / 1e9,
        asm.f_max / 1e9,
        n_points,
        p_ports,
        hdm_solves,
        rom_dim,
        max_residual,
        rom_sweep_seconds
    )];

    Ok(SolveResult {
        physics: "ems".into(),
        study: "driven".into(),
        order: asm.order as u8,
        dim: 3,
        solutions: vec![rep_sol],
        scalars,
        log,
        dofs: asm.n,
        seconds: t0.elapsed().as_secs_f64(),
    })
}
