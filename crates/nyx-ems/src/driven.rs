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

use crate::cavity::{assemble_maxwell_with_pml, basis3_keyed, gradient_full, safe_curved_jacobian, Ned3};
use crate::ibc::{assemble_ibc_and_jj_3d, reduce_real_csr, reduce_real_vec, AssembledJjPort, SheetCellKind, MU0};
use crate::mixed_ops::{
    assemble_ibc_and_jj_3d_mixed, assemble_lumped_ports_3d_mixed, assemble_maxwell_mixed, assemble_waveports_3d_mixed,
    axes_of_3d_periodic_mixed, iso_mu_inv_c3, match_vertices_mixed, simplex_cell_materials_3d,
    symmetry_facet_flags_mixed, validate_port_contacts_3d_mixed, vertex_keys_mixed, EmPmlTable, MixedCellMaterials,
};
use crate::modes::C0;
use crate::ports::{
    assemble_lumped_ports_3d, assemble_waveports_3d, validate_port_contacts_3d, AssembledLumpedPort, AssembledWavePort,
};
use crate::problem::{symmetry_facet_flags, Problem, Study};
use crate::rom::build_tls_meta;
use femcore::edge_mixed::{
    build_ams_p1_from_masks, build_ams_p2_from_masks, eval_cell_basis_3d, ned_mixed_periodic_relations, NedEvalScratch,
    NedMixed,
};
use femcore::fem::{geometry, CurvedTet};
use femcore::fem_mixed::IsoparametricCell;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::ams::{assemble_shifted_curl_curl, AmsFgmresOptions, AmsFgmresSolver, AmsP2Topology};
use femcore::linalg::csr::{c, CCsr, Coo, RCsr, C};
use femcore::linalg::equilibrate::{ruiz_symmetric_terms, scale_symmetric, scale_vec, RuizTerm, DEFAULT_RUIZ_PASSES};
use femcore::linalg::prom::{
    run_greedy_prom_sweep, s_to_z_and_y, ImageProm, PromLinearSolver, PromOperator, SparseDrivenSolver,
};
use femcore::mesh::mixed::{FacetType, MixedDiscretisation, QUAD_EDGES, TRI_EDGES};
use nyx_physics::{
    cell_fields_of, Discretisation, DrivenPortData, Field, FieldFrame, Mesh, MixedToSimplexProjector, Scalar, Solution,
    SolveResult,
};
use std::collections::BTreeSet;
use std::f64::consts::PI;
use std::time::Instant;

type C3 = [[C; 3]; 3];

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

/// Auxiliary-space topology for `solver == "ams_fgmres"` on $P_1$ or $P_2$ Nedelec spaces (`FR-33`).
pub enum DrivenAmsTopology {
    P1 { g_red: RCsr, edge_vecs: Vec<[f64; 3]> },
    P2(AmsP2Topology),
}

fn equilibrate_driven_system(
    nred: usize,
    f0: f64,
    k_tot_red: &CCsr,
    c_tot_red: &CCsr,
    m_tot_red: &CCsr,
    lumped_ports: &[AssembledLumpedPort],
    waveports: &[AssembledWavePort],
) -> (Vec<f64>, Vec<CCsr>, Vec<DrivenPortEntry>) {
    let omega_mid = 2.0 * PI * f0;
    let k0_mid = (omega_mid / C0).max(1.0);
    let k0_mid2 = k0_mid * k0_mid;
    let wp_betas_mid: Vec<f64> = waveports.iter().map(|wp| wp.beta_and_derivs(omega_mid).0).collect();

    let mut ruiz_terms = vec![
        RuizTerm::Complex(k_tot_red, 1.0),
        RuizTerm::Complex(c_tot_red, k0_mid),
        RuizTerm::Complex(m_tot_red, k0_mid2),
    ];
    for (w_idx, wp) in waveports.iter().enumerate() {
        if wp.w_red.nnz() > 0 {
            ruiz_terms.push(RuizTerm::Real(&wp.w_red, wp_betas_mid[w_idx]));
        }
    }
    let d_ruiz = ruiz_symmetric_terms(nred, &ruiz_terms, DEFAULT_RUIZ_PASSES);

    let mut terms_hat = Vec::with_capacity(3 + waveports.len());
    terms_hat.push(scale_symmetric(k_tot_red, &d_ruiz));
    terms_hat.push(scale_symmetric(c_tot_red, &d_ruiz));
    terms_hat.push(scale_symmetric(m_tot_red, &d_ruiz));
    for wp in waveports {
        terms_hat.push(scale_symmetric(&wp.w_red.to_complex(), &d_ruiz));
    }

    let mut ports = Vec::with_capacity(lumped_ports.len() + waveports.len());
    for lp in lumped_ports {
        ports.push(DrivenPortEntry {
            name: lp.bc.name.clone(),
            kind: DrivenPortKind::Lumped { z0: lp.bc.z0 },
            b_hat: scale_vec(&lp.b_red, &d_ruiz),
        });
    }
    for (wave_index, wp) in waveports.iter().enumerate() {
        ports.push(DrivenPortEntry {
            name: wp.bc.name.clone(),
            kind: DrivenPortKind::Wave { wave_index },
            b_hat: scale_vec(&wp.b_red, &d_ruiz),
        });
    }
    (d_ruiz, terms_hat, ports)
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
    pub solver: String,
    pub blr_tol: f64,
    pub d_ruiz: Vec<f64>,
    /// Affine operator terms $\hat{\mathbf{A}}_m = \mathbf{D} \mathbf{A}_m \mathbf{D}$ ($m = 0 \dots 3 + P_{\text{wave}} - 1$):
    /// `[K_tot, C_tot, M_tot, W_0, ..., W_{P_wave - 1}]`.
    pub terms_hat: Vec<CCsr>,
    /// Optional reduced interior gradient / prolongation topology for `solver == "ams_fgmres"` (`order == 1` or `2`).
    pub ams_data: Option<DrivenAmsTopology>,
    pub ports: Vec<DrivenPortEntry>,
    pub lumped_ports: Vec<AssembledLumpedPort>,
    pub waveports: Vec<AssembledWavePort>,
    pub jj_ports: Vec<AssembledJjPort>,
    pub cell_sheet: Vec<Option<SheetCellKind>>,
    pub mu_inv: Vec<C3>,
    pub mur_phys: Vec<f64>,
    pub has_sheets: bool,
    pub has_pml: bool,
    pub k_pml_red: CCsr,
    pub m_pml_red: CCsr,
    pub m_phys_red: RCsr,
    pub m_pml_phys_red: RCsr,
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

impl PromOperator for DrivenAssembly {
    fn n_dofs(&self) -> usize {
        self.nred
    }

    fn n_ports(&self) -> usize {
        self.ports.len()
    }

    fn terms_hat(&self) -> &[CCsr] {
        &self.terms_hat
    }

    fn term_coeffs(&self, omega: f64) -> Vec<(C, C, C)> {
        DrivenAssembly::term_coeffs(self, omega)
    }

    fn d_ruiz(&self) -> &[f64] {
        &self.d_ruiz
    }

    fn port_b_hat(&self, port: usize) -> Vec<C> {
        self.ports[port].b_hat.iter().map(|&b| c(b, 0.0)).collect()
    }

    fn port_gamma(&self, port: usize, omega: f64) -> (C, C, C) {
        let (_, g0, g1, g2) = self.ports[port].excitation_coeffs(omega, &self.waveports);
        (g0, g1, g2)
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
        solver,
        blr_tol,
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
    let (eps, mu_inv) = simplex_cell_materials_3d(mesh, problem)?;
    let f0 = 0.5 * (f_min + f_max);
    let pml_table = EmPmlTable::from_simplex(mesh, &problem.materials, f0, false, solver == "ams_fgmres");
    let pml_cells = pml_table.pml_cells_mask();
    let has_pml = pml_table.has_any();

    let mut matches = vec![];
    for p in &problem.periodic {
        matches.push(femcore::periodic::match_vertices(mesh, p)?);
    }
    let vkey = vertex_keys(mesh, &matches);
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let (k, m) = assemble_maxwell_with_pml(mesh, &ned, &eps, &mu_inv, has_pml.then_some(&pml_table));

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

    let (sym_pec, sym_pmc) = symmetry_facet_flags(mesh, problem.symmetry)?;
    let mut fixed_set = BTreeSet::new();
    let mut ext_pec_verts = BTreeSet::new();
    let mut cond_verts = BTreeSet::new();
    let mut cond_dsu_edges: Vec<(usize, usize)> = Vec::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        if exterior[f]
            && !periodic_labels.contains(&label)
            && (sym_pec[f]
                || (!problem.pmc.contains(&label)
                    && !sym_pmc[f]
                    && !ibc_labels.contains(&label)
                    && !jj_labels.contains(&label)
                    && !lp_labels.contains(&label)
                    && !wp_labels.contains(&label)))
        {
            let fv = mesh.facet(f);
            for d in ned.boundary_dofs(mesh, [fv[0], fv[1], fv[2]]) {
                fixed_set.insert(d);
            }
            for &v in fv {
                ext_pec_verts.insert(v as usize);
            }
        }
    }
    if let Some(ref ij) = ibc_jj {
        for &fv in &ij.pec_tris {
            for d in ned.boundary_dofs(mesh, fv) {
                fixed_set.insert(d);
            }
            for v in fv {
                cond_verts.insert(v as usize);
            }
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                cond_dsu_edges.push((a as usize, b as usize));
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
    let (d_ruiz, terms_hat, ports) =
        equilibrate_driven_system(nred, f0, &k_tot_red, &c_tot_red, &m_tot_red, &lumped_ports, &waveports);

    let ams_data = if solver == "ams_fgmres" && problem.periodic.is_empty() {
        let nv = mesh.n_points();
        let mut pec_vertex = ext_pec_verts.clone();
        let mut node_rel = Vec::new();
        if !cond_verts.is_empty() {
            let mut parent: Vec<usize> = (0..nv).collect();
            let find = |p: &mut [usize], mut i: usize| -> usize {
                while p[i] != i {
                    p[i] = p[p[i]];
                    i = p[i];
                }
                i
            };
            for &(a, b) in &cond_dsu_edges {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    let (lo, hi) = (ra.min(rb), ra.max(rb));
                    parent[hi] = lo;
                }
            }
            let mut comp_touches_ext = BTreeSet::new();
            for &v in &cond_verts {
                if ext_pec_verts.contains(&v) {
                    comp_touches_ext.insert(find(&mut parent, v));
                }
            }
            for &v in &cond_verts {
                let root = find(&mut parent, v);
                if comp_touches_ext.contains(&root) {
                    pec_vertex.insert(v);
                } else if v != root {
                    node_rel.push((v, root, c(1.0, 0.0)));
                }
            }
        }
        let mut node_fixed: Vec<usize> = pec_vertex.into_iter().collect();
        if node_fixed.is_empty() && nv > 0 {
            node_fixed.push(0);
        }
        let nred_nodes = femcore::periodic::Reduction::build(nv, &node_rel, &node_fixed);
        let nn_p1 = nred_nodes.n_reduced();
        if nn_p1 == 0 {
            None
        } else if order == 1 {
            let g_full = gradient_full(mesh, &ned);
            let g_nodes = femcore::linalg::amg_elastic::spmm(&g_full, &nred_nodes.prolong().real_part());
            let (mut indptr, mut indices, mut data) = (vec![0usize], Vec::new(), Vec::new());
            let mut edge_vecs = Vec::with_capacity(nred);
            let edges = mesh.edges();
            for &r in &red.roots {
                for q in g_nodes.indptr()[r]..g_nodes.indptr()[r + 1] {
                    indices.push(g_nodes.indices()[q] as usize);
                    data.push(g_nodes.values()[q]);
                }
                indptr.push(indices.len());
                let (a, b) = edges[r];
                let (lo, hi) = if (ned.vkey[a as usize], a) < (ned.vkey[b as usize], b) { (a, b) } else { (b, a) };
                let plo = mesh.point(lo as usize);
                let phi = mesh.point(hi as usize);
                edge_vecs.push([phi[0] - plo[0], phi[1] - plo[1], phi[2] - plo[2]]);
            }
            let g_red = RCsr::from_raw(nred, nn_p1, indptr, indices, data);
            Some(DrivenAmsTopology::P1 { g_red, edge_vecs })
        } else {
            let mut red_of_full = vec![None; n];
            for (red_idx, &full_dof) in red.roots.iter().enumerate() {
                red_of_full[full_dof] = Some(red_idx);
            }
            let mut free_edges = Vec::new();
            for e in 0..ned.n_edges {
                if let Some(row_p2) = red_of_full[2 * e] {
                    free_edges.push((e, row_p2));
                }
            }
            if free_edges.is_empty() {
                None
            } else {
                let ne_p1 = free_edges.len();
                let p_nodes_p1 = nred_nodes.prolong().real_part();
                let vert_col: Vec<Option<usize>> =
                    (0..nv).map(|v| p_nodes_p1.row(v).0.first().map(|&col| col as usize)).collect();
                let mut coo_p = Coo::new(nred, ne_p1);
                let mut coo_g1 = Coo::new(ne_p1, nn_p1);
                let mut free_bubbles = Vec::new();
                for e in 0..ned.n_edges {
                    if let Some(row_bub) = red_of_full[2 * e + 1] {
                        free_bubbles.push(row_bub);
                    }
                }
                let nn_p2 = nn_p1 + free_bubbles.len();
                let mut coo_g2 = Coo::new(nred, nn_p2);
                let mut edge_vectors = Vec::with_capacity(ne_p1);
                let edges = mesh.edges();
                for (col_e, &(e, row_p2)) in free_edges.iter().enumerate() {
                    coo_p.push(row_p2, col_e, 1.0);
                    let (a, b) = edges[e];
                    let (lo, hi) = if (ned.vkey[a as usize], a) < (ned.vkey[b as usize], b) { (a, b) } else { (b, a) };
                    let (c_hi_opt, c_lo_opt) = (vert_col[hi as usize], vert_col[lo as usize]);
                    if c_hi_opt != c_lo_opt {
                        if let Some(c_hi) = c_hi_opt {
                            coo_g1.push(col_e, c_hi, 1.0);
                            coo_g2.push(row_p2, c_hi, 1.0);
                        }
                        if let Some(c_lo) = c_lo_opt {
                            coo_g1.push(col_e, c_lo, -1.0);
                            coo_g2.push(row_p2, c_lo, -1.0);
                        }
                    }
                    let plo = mesh.point(lo as usize);
                    let phi = mesh.point(hi as usize);
                    edge_vectors.push([phi[0] - plo[0], phi[1] - plo[1], phi[2] - plo[2]]);
                }
                for (bub_idx, &row_bub) in free_bubbles.iter().enumerate() {
                    coo_g2.push(row_bub, nn_p1 + bub_idx, 1.0);
                }
                let mut blocks = Vec::with_capacity(ned.n_edges + ned.n_faces);
                for e in 0..ned.n_edges {
                    let blk: Vec<usize> = [red_of_full[2 * e], red_of_full[2 * e + 1]].into_iter().flatten().collect();
                    if !blk.is_empty() {
                        blocks.push(blk);
                    }
                }
                for f in 0..ned.n_faces {
                    let base = 2 * ned.n_edges + 2 * f;
                    let blk: Vec<usize> = [red_of_full[base], red_of_full[base + 1]].into_iter().flatten().collect();
                    if !blk.is_empty() {
                        blocks.push(blk);
                    }
                }
                Some(DrivenAmsTopology::P2(AmsP2Topology {
                    p1_to_p2: coo_p.to_csr(),
                    gradient_p1: coo_g1.to_csr(),
                    gradient_p2: Some(coo_g2.to_csr()),
                    edge_vectors,
                    blocks,
                }))
            }
        }
    } else {
        None
    };

    let (k_pml_red, m_pml_red, m_phys_red, m_pml_phys_red) = if has_pml {
        let zero_c3 = [[c(0.0, 0.0); 3]; 3];
        let mut eps_pml_only = Vec::with_capacity(mesh.n_cells());
        let mut mui_pml_only = Vec::with_capacity(mesh.n_cells());
        let mut eps_phys_only = Vec::with_capacity(mesh.n_cells());
        let mut eps_pml_unstr = Vec::with_capacity(mesh.n_cells());
        let mut mui_zero = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            let mut er_c = zero_c3;
            for a in 0..3 {
                for b in 0..3 {
                    er_c[a][b] = c(eps_phys_re[ci][a][b], 0.0);
                }
            }
            mui_zero.push(zero_c3);
            if pml_cells[ci] {
                eps_pml_only.push(eps[ci]);
                mui_pml_only.push(mu_inv[ci]);
                eps_phys_only.push(zero_c3);
                eps_pml_unstr.push(er_c);
            } else {
                eps_pml_only.push(zero_c3);
                mui_pml_only.push(zero_c3);
                eps_phys_only.push(er_c);
                eps_pml_unstr.push(zero_c3);
            }
        }
        let (kp, mp) = assemble_maxwell_with_pml(mesh, &ned, &eps_pml_only, &mui_pml_only, Some(&pml_table));
        let (_, m_ph) = assemble_maxwell_with_pml(mesh, &ned, &eps_phys_only, &mui_zero, None);
        let (_, m_pm_u) = assemble_maxwell_with_pml(mesh, &ned, &eps_pml_unstr, &mui_zero, None);
        (
            femcore::periodic::congruence(&pmat, &kp),
            femcore::periodic::congruence(&pmat, &mp),
            femcore::periodic::congruence(&pmat, &m_ph).real_part(),
            femcore::periodic::congruence(&pmat, &m_pm_u).real_part(),
        )
    } else {
        (CCsr::zeros(0, 0), CCsr::zeros(0, 0), RCsr::zeros(0, 0), RCsr::zeros(0, 0))
    };
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
        solver,
        blr_tol,
        d_ruiz,
        terms_hat,
        ams_data,
        ports,
        lumped_ports,
        waveports,
        jj_ports,
        cell_sheet,
        mu_inv,
        mur_phys: mur,
        has_sheets,
        has_pml,
        k_pml_red,
        m_pml_red,
        m_phys_red,
        m_pml_phys_red,
    })
}

/// Unified linear solver backend for 3D driven frequency-domain solves (FR-33).
enum DrivenLinearSolver {
    /// Sparse direct / LBLT / BLR multifrontal LU via shared `femcore::linalg::prom::SparseDrivenSolver`.
    Sparse(Box<SparseDrivenSolver>),
    /// Auxiliary-space Maxwell (`Ams<C>` / `AmsP2<C>`) preconditioned deterministic FGMRES (`ams_fgmres`)
    /// on the unscaled complex curl-curl operator (`order == 1` and `order == 2`).
    AmsFgmres { solver: Box<AmsFgmresSolver<C>>, d_ruiz: Vec<f64>, symbolic_count: usize },
}

fn unscale_sym_c(mat: &CCsr, d_ruiz: &[f64]) -> CCsr {
    let mut out = mat.clone();
    let n = out.nrows();
    for i in 0..n {
        let inv_di = 1.0 / d_ruiz[i];
        let range = out.row_range(i);
        let cols: Vec<u32> = out.indices()[range.clone()].to_vec();
        let vals = &mut out.values_mut()[range];
        for (v, &col) in vals.iter_mut().zip(&cols) {
            *v *= inv_di / d_ruiz[col as usize];
        }
    }
    out
}

impl DrivenLinearSolver {
    fn build_or_update(slot: &mut Option<Self>, asm: &DrivenAssembly, omega: f64, a_mat: &CCsr) -> Result<(), String> {
        match asm.solver.as_str() {
            "ams_fgmres" => {
                let Some(ref ams_topo) = asm.ams_data else {
                    return Err(
                        "driven solver 'ams_fgmres' requires a 3D mesh with interior non-PEC vertices and no periodic boundaries (use 'blr', 'lblt', or 'auto' otherwise)"
                            .into(),
                    );
                };
                let k0 = omega / C0;
                let coeffs = asm.term_coeffs(omega);
                let mut k_eff = unscale_sym_c(&asm.terms_hat[0], &asm.d_ruiz);
                for (m_idx, &(alpha, _, _)) in coeffs.iter().enumerate().skip(3) {
                    if asm.terms_hat[m_idx].nnz() > 0 && (alpha.re != 0.0 || alpha.im != 0.0) {
                        let w_unscaled = unscale_sym_c(&asm.terms_hat[m_idx], &asm.d_ruiz);
                        k_eff = k_eff.add_scaled(&w_unscaled, alpha);
                    }
                }
                let c_unscaled = (asm.terms_hat[1].nnz() > 0).then(|| unscale_sym_c(&asm.terms_hat[1], &asm.d_ruiz));
                let m_unscaled = unscale_sym_c(&asm.terms_hat[2], &asm.d_ruiz);

                match slot.as_mut() {
                    Some(Self::AmsFgmres { solver, .. }) => {
                        let (a_sys, a_pc) = assemble_shifted_curl_curl(&k_eff, c_unscaled.as_ref(), &m_unscaled, k0);
                        solver
                            .update_matrix(&a_sys, Some(&a_pc))
                            .map_err(|e| format!("AMS+FGMRES preconditioner update failed: {e}"))?;
                        Ok(())
                    }
                    _ => {
                        *slot = None;
                        let mut opts = AmsFgmresOptions {
                            rtol: 1e-8,
                            atol: 1e-30,
                            max_it: 800,
                            restart: 100,
                            flexible: true,
                            ..AmsFgmresOptions::default()
                        };
                        opts.ams.smoother_sweeps = 2;
                        opts.ams_p2.smoother_sweeps = 2;
                        opts.ams_p2.coarse_sweeps = 2;
                        opts.ams_p2.coarse_ams.smoother_sweeps = 2;
                        let ams_solver = match ams_topo {
                            DrivenAmsTopology::P1 { g_red, edge_vecs } => AmsFgmresSolver::new_shifted_complex_p1(
                                &k_eff,
                                c_unscaled.as_ref(),
                                &m_unscaled,
                                k0,
                                g_red,
                                edge_vecs,
                                opts,
                            ),
                            DrivenAmsTopology::P2(topo) => AmsFgmresSolver::new_shifted_complex_p2(
                                &k_eff,
                                c_unscaled.as_ref(),
                                &m_unscaled,
                                k0,
                                topo,
                                opts,
                            ),
                        }
                        .map_err(|e| format!("AMS+FGMRES setup failed: {e}"))?;
                        *slot = Some(Self::AmsFgmres {
                            solver: Box::new(ams_solver),
                            d_ruiz: asm.d_ruiz.clone(),
                            symbolic_count: 1,
                        });
                        Ok(())
                    }
                }
            }
            _ => match slot.as_mut() {
                Some(Self::Sparse(sp)) => sp.update_matrix(a_mat),
                _ => {
                    *slot = None;
                    let mut sp = SparseDrivenSolver::new(&asm.solver, asm.blr_tol);
                    sp.update_matrix(a_mat)?;
                    *slot = Some(Self::Sparse(Box::new(sp)));
                    Ok(())
                }
            },
        }
    }

    fn symbolic_count(&self) -> usize {
        match self {
            Self::Sparse(sp) => sp.symbolic_count(),
            Self::AmsFgmres { symbolic_count, .. } => *symbolic_count,
        }
    }

    fn backend_name(&self) -> String {
        match self {
            Self::Sparse(sp) => sp.backend_name(),
            Self::AmsFgmres { .. } => "AmsFgmres".to_string(),
        }
    }

    fn solve_many(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        match self {
            Self::Sparse(sp) => sp.solve_many(rhs_cols),
            Self::AmsFgmres { solver, d_ruiz, .. } => {
                let rhs_unscaled: Vec<Vec<C>> =
                    rhs_cols.iter().map(|b_hat| b_hat.iter().zip(d_ruiz).map(|(&bi, &di)| bi / di).collect()).collect();
                let (xs_unscaled, st) = solver.solve_many(&rhs_unscaled);
                if !st.all_converged {
                    return Err(format!(
                        "driven solver 'ams_fgmres' failed to converge after {} total iterations (max relative residual {:.3e})",
                        st.total_iterations, st.max_relative_residual
                    ));
                }
                let out: Vec<Vec<C>> = xs_unscaled
                    .into_iter()
                    .map(|x_unscaled| x_unscaled.into_iter().zip(d_ruiz).map(|(xi, &di)| xi / di).collect())
                    .collect();
                Ok(out)
            }
        }
    }
}

struct EmsPromSolver<'a> {
    asm: &'a DrivenAssembly,
    slot: Option<DrivenLinearSolver>,
}

impl PromLinearSolver for EmsPromSolver<'_> {
    fn update_and_solve(&mut self, omega: f64, a_mat: &CCsr, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        DrivenLinearSolver::build_or_update(&mut self.slot, self.asm, omega, a_mat)?;
        self.slot
            .as_ref()
            .ok_or_else(|| "driven linear solver missing after build_or_update".to_string())?
            .solve_many(rhs_cols)
    }

    fn solve_same_factor(&self, rhs_cols: &[Vec<C>]) -> Result<Vec<Vec<C>>, String> {
        self.slot
            .as_ref()
            .ok_or_else(|| "driven linear solver not factored before solve_same_factor".to_string())?
            .solve_many(rhs_cols)
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
        let gid = asm.ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for kk in 0..3 {
                gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
            }
        }
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        for (a, l) in lamn.iter().enumerate() {
            let lam = [l[0], l[1], l[2], l[3]];
            match &curved {
                Some(ct) => {
                    let gq = safe_curved_jacobian(ct, &lam, &gs);
                    basis3_keyed(asm.order, &lam, &gq, &gid, &mut w, &mut cu);
                }
                None => basis3_keyed(asm.order, &lam, &gs, &gid, &mut w, &mut cu),
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
        let mui = c(1.0 / asm.mur_phys[ci], 0.0);
        for a in 0..nloc {
            let idx = ci * nloc + a;
            let ce = curl_e[idx];
            h_field[idx] = [ce[0] * mui * h_prefactor, ce[1] * mui * h_prefactor, ce[2] * mui * h_prefactor];
        }
    }
    for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
        let z: Vec<C> = h_field.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }

    Solution { label, fields, cell_fields: cell_fields_of(&cells), value: Some(freq_hz), eigenvalue: None }
}

/// Evaluates `(P_pml_watts, layer_fraction)` at `freq_hz` for unscaled reduced field `x_red`.
fn compute_driven_pml_metrics(asm: &DrivenAssembly, freq_hz: f64, x_red: &[C]) -> (f64, f64) {
    if !asm.has_pml || x_red.is_empty() {
        return (0.0, 0.0);
    }
    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let k0_2 = (omega / C0).powi(2);
    let kx = asm.k_pml_red.matvec(x_red);
    let mx = asm.m_pml_red.matvec(x_red);
    let mut quad_im = 0.0_f64;
    for i in 0..x_red.len() {
        quad_im += (x_red[i].conj() * (kx[i] - mx[i] * k0_2)).im;
    }
    let p_pml_w = (quad_im / (omega * MU0)).max(0.0);

    let quad_real = |mat: &RCsr| -> f64 {
        let mut acc = 0.0;
        for i in 0..mat.nrows() {
            let xi = x_red[i];
            let (cols, vals) = mat.row(i);
            let mut rs = c(0.0, 0.0);
            for (&col, &v) in cols.iter().zip(vals) {
                rs += x_red[col as usize] * v;
            }
            acc += (xi.conj() * rs).re;
        }
        acc.max(0.0)
    };
    let u_phys = quad_real(&asm.m_phys_red);
    let u_pml = quad_real(&asm.m_pml_phys_red);
    let u_tot = u_phys + u_pml;
    let layer_fraction = if u_tot > 0.0 { (u_pml / u_tot).clamp(0.0, 1.0) } else { 0.0 };
    (p_pml_w, layer_fraction)
}

struct DrivenSweepCoreOut {
    driven_data: DrivenPortData,
    symbolic_count: usize,
    backend_name: String,
    rom_sweep_seconds: f64,
    s21_min_db: f64,
    f_res_hz: f64,
    rep_x_red: Vec<C>,
}

fn run_driven_sweep_core(asm: &DrivenAssembly) -> Result<DrivenSweepCoreOut, String> {
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
    let backend_name;
    let mut rom_dim = 0usize;
    let mut anchor_freqs_hz = Vec::new();
    let mut max_residual = 0.0_f64;
    let mut rom_sweep_seconds = 0.0_f64;
    let mut rep_freq_idx = n_points / 2;
    let mut rep_x_red: Vec<C> = vec![c(0.0, 0.0); asm.nred];

    if let Some(tol) = asm.prom_tol {
        let mut solver = EmsPromSolver { asm, slot: None };
        let mut min_s_metric = f64::INFINITY;
        const SWEEP_CHUNK: usize = 32;
        let summary = run_greedy_prom_sweep(
            asm,
            &mut solver,
            &freqs_hz,
            tol,
            asm.max_hdm_solves,
            asm.strict_prom_tol,
            &asm.seed_freqs_hz,
            |band, b_omegas, prom: &ImageProm<'_, DrivenAssembly>, xr_band| {
                let nb_pts = b_omegas.len();
                let t_post = Instant::now();
                let indices: Vec<usize> = (0..nb_pts).collect();
                let b_dim = prom.dim();
                let post_chunks = femcore::par::map_chunks(&indices, SWEEP_CHUNK, |chunk| {
                    chunk
                        .iter()
                        .map(|&k| {
                            let w = b_omegas[k];
                            let xr_k = &xr_band[k];
                            let mut z_wave_k = Vec::with_capacity(p_ports);
                            let mut z_circ_k = Vec::with_capacity(p_ports);
                            let mut z0_pairs_k = Vec::with_capacity(p_ports);
                            for port in &asm.ports {
                                let (z_wave_p, _, _, _) = port.excitation_coeffs(w, &asm.waveports);
                                let z_circ_p = port.circuit_z0(z_wave_p, &asm.waveports);
                                z_wave_k.push(z_wave_p);
                                z_circ_k.push(z_circ_p);
                                z0_pairs_k.push([z_circ_p, 0.0]);
                            }
                            let mut s_k = vec![vec![c(0.0, 0.0); p_ports]; p_ports];
                            for p in 0..p_ports {
                                let inv_sqrt_z_wave_p = 1.0 / z_wave_k[p].max(1e-30).sqrt();
                                for q in 0..p_ports {
                                    let mut v_pq = c(0.0, 0.0);
                                    for i in 0..b_dim {
                                        v_pq += xr_k[q][i] * prom.b_r_drive[i][p];
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
                let mut band_best_local = None;
                for (k_loc, (z0_pairs_k, s_out, z_out, y_out, metric)) in post_chunks.into_iter().flatten().enumerate()
                {
                    if metric < min_s_metric {
                        min_s_metric = metric;
                        rep_freq_idx = band.start + k_loc;
                        band_best_local = Some(k_loc);
                    }
                    z0_ohm.push(z0_pairs_k);
                    s_matrix.push(s_out);
                    z_matrix.push(z_out);
                    y_matrix.push(y_out);
                }
                if let Some(kl) = band_best_local {
                    rep_x_red = prom.reconstruct_unscaled(&xr_band[kl][0]);
                }
                rom_sweep_seconds += t_post.elapsed().as_secs_f64();
                Ok(())
            },
        )?;
        symbolic_count = solver.slot.as_ref().map_or(0, DrivenLinearSolver::symbolic_count);
        backend_name = solver.slot.as_ref().map_or_else(|| "None".to_string(), DrivenLinearSolver::backend_name);
        // FR-33: Drop sparse factor buffers immediately after anchor solves finish, before post-sweep & 3D CellFields.
        drop(solver.slot.take());
        hdm_solves = summary.hdm_solves;
        rom_dim = summary.rom_dim;
        anchor_freqs_hz = summary.anchor_freqs_hz;
        max_residual = summary.max_residual;
        rom_sweep_seconds += summary.rom_sweep_seconds;
    } else {
        // Direct point-by-point HDM solve (prom_tol = None), reusing symbolic factorization across all frequencies
        let mut solver: Option<DrivenLinearSolver> = None;
        let mut min_s_metric = f64::INFINITY;
        for (k, &w) in omegas.iter().enumerate() {
            femcore::progress::check()?;
            let a_mat = asm.assemble_matrix_at(w);
            DrivenLinearSolver::build_or_update(&mut solver, asm, w, &a_mat)?;
            let s_ref =
                solver.as_ref().ok_or_else(|| "driven linear solver missing after build_or_update".to_string())?;
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

            let y_cols = s_ref.solve_many(&rhs_cols)?;
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
        symbolic_count = solver.as_ref().map_or(0, DrivenLinearSolver::symbolic_count);
        backend_name = solver.as_ref().map_or_else(|| "None".to_string(), DrivenLinearSolver::backend_name);
        // FR-33: Drop sparse solver before `build_driven_solution` allocates 3D `CellFields`.
        drop(solver);
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

    Ok(DrivenSweepCoreOut {
        driven_data,
        symbolic_count,
        backend_name,
        rom_sweep_seconds,
        s21_min_db,
        f_res_hz,
        rep_x_red,
    })
}

/// Solves a 3D `Study::Driven` frequency sweep via Hermite Tangential Galerkin PROM (`prom_tol = Some(tol)`)
/// or direct point-by-point HDM factorizations (`prom_tol = None`).
pub fn solve_driven(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_emw_driven(mesh, problem)?;
    let out = run_driven_sweep_core(&asm)?;
    let (p_pml_w, layer_fraction) = compute_driven_pml_metrics(&asm, out.f_res_hz, &out.rep_x_red);

    let rep_sol = build_driven_solution(
        mesh,
        &asm,
        out.f_res_hz,
        &out.rep_x_red,
        format!("driven ({}): {:.6} GHz", asm.ports[0].name, out.f_res_hz / 1e9),
    );

    let n_points = out.driven_data.freqs_hz.len();
    let p_ports = asm.ports.len();
    let hdm_solves = out.driven_data.hdm_solves;
    let rom_dim = out.driven_data.rom_dim;
    let max_residual = out.driven_data.max_residual;

    let (axes, left_handed) = femcore::section::axes_of_3d_periodic(&problem.periodic, mesh);
    let mut scalars = vec![
        out.driven_data.scalar(),
        ("hdm_solves".into(), Scalar::Num(hdm_solves as f64)),
        ("symbolic_count".into(), Scalar::Num(out.symbolic_count as f64)),
        ("rom_dim".into(), Scalar::Num(rom_dim as f64)),
        ("max_prom_residual".into(), Scalar::Num(max_residual)),
        ("rom_sweep_seconds".into(), Scalar::Num(out.rom_sweep_seconds)),
        ("S21_min_dB".into(), Scalar::Num(out.s21_min_db)),
        ("f_res_Hz".into(), Scalar::Num(out.f_res_hz)),
        ("P_pml (W)".into(), Scalar::Num(p_pml_w)),
        ("layer fraction 1".into(), Scalar::Num(layer_fraction)),
        FieldFrame {
            lateral: axes[0],
            normal: axes[1],
            propagation: axes[2],
            left_handed,
            outside: crate::optical::absorbing_domains(problem),
        }
        .scalar(),
    ];
    if asm.has_sheets {
        scalars.push(build_tls_meta(mesh, problem).scalar());
    }

    let log = vec![format!(
        "driven sweep [{:.4}..{:.4} GHz, {} points, {} port(s), solver={}, backend={}]: {} HDM solve(s), ROM dim {}, max residual {:.2e}, ROM sweep {:.4} s",
        asm.f_min / 1e9,
        asm.f_max / 1e9,
        n_points,
        p_ports,
        asm.solver,
        out.backend_name,
        hdm_solves,
        rom_dim,
        max_residual,
        out.rom_sweep_seconds
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
        ..Default::default()
    })
}

// ------------------------------------------------------------------ Native 3D MixedDiscretisation (`Tet + Prism + Hex`) driven solver

/// Solves a 3D `Study::Driven` sweep natively on `Tet + Prism + Hex` (`MixedDiscretisation`).
pub fn solve_driven_mixed(mesh: &Mesh, mixed: &MixedDiscretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let (asm, ned_mixed) = assemble_emw_driven_mixed(mixed, problem)?;
    let out = run_driven_sweep_core(&asm)?;
    let (p_pml_w, layer_fraction) = compute_driven_pml_metrics(&asm, out.f_res_hz, &out.rep_x_red);

    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, mixed).map_err(|e| e.to_string())?;
    let rep_sol = build_driven_solution_mixed(
        mixed,
        &ned_mixed,
        &simplex,
        &proj,
        &asm,
        out.f_res_hz,
        &out.rep_x_red,
        format!("driven ({}): {:.6} GHz", asm.ports[0].name, out.f_res_hz / 1e9),
    );

    let n_points = out.driven_data.freqs_hz.len();
    let p_ports = asm.ports.len();
    let hdm_solves = out.driven_data.hdm_solves;
    let rom_dim = out.driven_data.rom_dim;
    let max_residual = out.driven_data.max_residual;

    let (axes, left_handed) = axes_of_3d_periodic_mixed(&problem.periodic, mixed);
    let mut scalars = vec![
        out.driven_data.scalar(),
        ("hdm_solves".into(), Scalar::Num(hdm_solves as f64)),
        ("symbolic_count".into(), Scalar::Num(out.symbolic_count as f64)),
        ("rom_dim".into(), Scalar::Num(rom_dim as f64)),
        ("max_prom_residual".into(), Scalar::Num(max_residual)),
        ("rom_sweep_seconds".into(), Scalar::Num(out.rom_sweep_seconds)),
        ("S21_min_dB".into(), Scalar::Num(out.s21_min_db)),
        ("f_res_Hz".into(), Scalar::Num(out.f_res_hz)),
        ("P_pml (W)".into(), Scalar::Num(p_pml_w)),
        ("layer fraction 1".into(), Scalar::Num(layer_fraction)),
        FieldFrame {
            lateral: axes[0],
            normal: axes[1],
            propagation: axes[2],
            left_handed,
            outside: crate::optical::absorbing_domains(problem),
        }
        .scalar(),
    ];
    if asm.has_sheets {
        scalars.push(build_tls_meta(&simplex, problem).scalar());
    }

    let log = vec![format!(
        "driven mixed sweep [{:.4}..{:.4} GHz, {} points, {} port(s), solver={}, backend={}]: {} HDM solve(s), ROM dim {}, max residual {:.2e}, ROM sweep {:.4} s",
        asm.f_min / 1e9,
        asm.f_max / 1e9,
        n_points,
        p_ports,
        asm.solver,
        out.backend_name,
        hdm_solves,
        rom_dim,
        max_residual,
        out.rom_sweep_seconds
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
        ..Default::default()
    })
}

pub fn assemble_emw_driven_mixed(
    mesh: &MixedDiscretisation,
    problem: &Problem,
) -> Result<(DrivenAssembly, NedMixed), String> {
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
        solver,
        blr_tol,
    } = problem.study.clone()
    else {
        return Err("assemble_emw_driven_mixed requires a Driven study".into());
    };
    let (f_min, f_max, n_points) = if let Some(fs) = freqs_hz.as_ref().filter(|v| !v.is_empty()) {
        (fs[0], fs[fs.len() - 1], fs.len())
    } else {
        (study_f_min, study_f_max, study_n_points)
    };
    if problem.lumped_ports.is_empty() && problem.waveports.is_empty() {
        return Err("driven study requires at least one 'lumped_port' or 'waveport' boundary condition".into());
    }
    validate_port_contacts_3d_mixed(mesh, problem)?;
    let order = if problem.order == 1 { 1 } else { 2 };
    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let f0 = 0.5 * (f_min + f_max);
    let pml_table = EmPmlTable::from_mixed(mesh, &problem.materials, f0, false, solver == "ams_fgmres");
    let pml_cells = pml_table.pml_cells_mask();
    let has_pml = pml_table.has_any();

    let mut matches = Vec::with_capacity(problem.periodic.len());
    for p in &problem.periodic {
        matches.push(match_vertices_mixed(mesh, p)?);
    }
    let vkey = vertex_keys_mixed(mesh, &matches);
    let ned = NedMixed::with_vkey(mesh, order, vkey);
    let n = ned.n_dofs;
    let (k, m) = assemble_maxwell_mixed(mesh, &ned, &mats.eps_0, &mats.mu_inv_0, has_pml.then_some(&pml_table));

    let has_sheets = !problem.pec.is_empty() || !problem.surface_impedance.is_empty() || !problem.lumped_jj.is_empty();
    let mut ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d_mixed(
            mesh,
            &ned,
            &problem.pec,
            &problem.surface_impedance,
            &problem.lumped_jj,
            f0,
        )?)
    } else {
        None
    };

    let (mut lumped_ports, _) = assemble_lumped_ports_3d_mixed(mesh, &ned, &problem.pec, &problem.lumped_ports)?;
    let (mut waveports, _) = assemble_waveports_3d_mixed(mesh, &ned, problem, f_min, f_max)?;

    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let wp_labels: Vec<i32> = problem.waveports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let (sym_pec, sym_pmc) = symmetry_facet_flags_mixed(mesh, problem.symmetry)?;

    let nv = mesh.n_points();
    let ne = mesh.n_edges();
    let mut fixed_set = BTreeSet::new();
    let mut fixed_verts = vec![false; nv];
    let mut fixed_edges = vec![false; ne];
    let mut fixed_quad_faces = vec![false; mesh.n_quad_faces()];
    let mut mark_facet = |f: usize, fv: &[u32], ft: FacetType| {
        for &v in fv {
            fixed_verts[v as usize] = true;
        }
        let le: &[(usize, usize)] = match ft {
            FacetType::Edge2 => &[(0, 1)],
            FacetType::Tri3 => &TRI_EDGES,
            FacetType::Quad4 => &QUAD_EDGES,
        };
        for &(i, j) in le {
            if let Some(eid) = mesh.edge_index(fv[i], fv[j]) {
                fixed_edges[eid] = true;
            }
        }
        if ft == FacetType::Quad4 {
            if let Some(qf) = mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]) {
                fixed_quad_faces[qf] = true;
            }
        }
        let _ = f;
    };
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        if mesh.exterior_facets[f]
            && !periodic_labels.contains(&label)
            && (sym_pec[f]
                || (!problem.pmc.contains(&label)
                    && !sym_pmc[f]
                    && !ibc_labels.contains(&label)
                    && !jj_labels.contains(&label)
                    && !lp_labels.contains(&label)
                    && !wp_labels.contains(&label)))
        {
            for d in ned.facet_dofs(mesh, f) {
                fixed_set.insert(d);
            }
            mark_facet(f, mesh.facet(f), mesh.facet_type(f));
        }
    }
    if let Some((_, ref pec_fcs, _, _)) = ibc_jj {
        for &f in pec_fcs {
            for d in ned.facet_dofs(mesh, f) {
                fixed_set.insert(d);
            }
            mark_facet(f, mesh.facet(f), mesh.facet_type(f));
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let mut relations = Vec::new();
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        relations.extend(ned_mixed_periodic_relations(mesh, &ned, &mt.src_facets, &mt.map, p.phase)?);
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

    let (jj_ports, cell_sheet) = if let Some((ref mut ij, _, _, _)) = ibc_jj {
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

    let (k_pml_red, m_pml_red, m_phys_red, m_pml_phys_red) = if has_pml {
        let zero_c3 = [[c(0.0, 0.0); 3]; 3];
        let mut eps_pml_only = Vec::with_capacity(mesh.n_cells());
        let mut mui_pml_only = Vec::with_capacity(mesh.n_cells());
        let mut eps_phys_only = Vec::with_capacity(mesh.n_cells());
        let mut eps_pml_unstr = Vec::with_capacity(mesh.n_cells());
        let mut mui_zero = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            let mut er_c = zero_c3;
            for a in 0..3 {
                for b in 0..3 {
                    er_c[a][b] = c(mats.eps_phys_re[ci][a][b], 0.0);
                }
            }
            mui_zero.push(zero_c3);
            if pml_cells[ci] {
                eps_pml_only.push(mats.eps_0[ci]);
                mui_pml_only.push(mats.mu_inv_0[ci]);
                eps_phys_only.push(zero_c3);
                eps_pml_unstr.push(er_c);
            } else {
                eps_pml_only.push(zero_c3);
                mui_pml_only.push(zero_c3);
                eps_phys_only.push(er_c);
                eps_pml_unstr.push(zero_c3);
            }
        }
        let (kp, mp) = assemble_maxwell_mixed(mesh, &ned, &eps_pml_only, &mui_pml_only, Some(&pml_table));
        let (_, m_ph) = assemble_maxwell_mixed(mesh, &ned, &eps_phys_only, &mui_zero, None);
        let (_, m_pm_u) = assemble_maxwell_mixed(mesh, &ned, &eps_pml_unstr, &mui_zero, None);
        (
            femcore::periodic::congruence(&pmat, &kp),
            femcore::periodic::congruence(&pmat, &mp),
            femcore::periodic::congruence(&pmat, &m_ph).real_part(),
            femcore::periodic::congruence(&pmat, &m_pm_u).real_part(),
        )
    } else {
        (CCsr::zeros(0, 0), CCsr::zeros(0, 0), RCsr::zeros(0, 0), RCsr::zeros(0, 0))
    };

    let (d_ruiz, terms_hat, ports) =
        equilibrate_driven_system(nred, f0, &k_tot_red, &c_tot_red, &m_tot_red, &lumped_ports, &waveports);

    let ams_data = if solver == "ams_fgmres" && problem.periodic.is_empty() {
        if !fixed_verts.iter().any(|&v| v) && nv > 0 {
            fixed_verts[0] = true;
        }
        if !fixed_verts.iter().all(|&v| v) {
            if order == 1 {
                let p1 = build_ams_p1_from_masks(mesh, &ned.vkey, &fixed_edges, &fixed_verts);
                (p1.free_edges == red.roots)
                    .then_some(DrivenAmsTopology::P1 { g_red: p1.gradient_p1, edge_vecs: p1.edge_vectors })
            } else {
                let mut fixed_dofs = vec![false; n];
                for &d in &fixed {
                    fixed_dofs[d] = true;
                }
                let p2 =
                    build_ams_p2_from_masks(mesh, &ned, &fixed_dofs, &fixed_edges, &fixed_verts, &fixed_quad_faces);
                (p2.free_dofs_p2 == red.roots).then_some(DrivenAmsTopology::P2(p2.topology))
            }
        } else {
            None
        }
    } else {
        None
    };

    let _ = iso_mu_inv_c3;
    Ok((
        DrivenAssembly {
            order,
            n,
            nred,
            ned: Ned3::empty(order),
            red,
            f_min,
            f_max,
            n_points,
            freqs_hz,
            prom_tol,
            max_hdm_solves,
            strict_prom_tol,
            seed_freqs_hz,
            solver,
            blr_tol,
            d_ruiz,
            terms_hat,
            ams_data,
            ports,
            lumped_ports,
            waveports,
            jj_ports,
            cell_sheet,
            mu_inv: mats.mu_inv_0,
            mur_phys: mats.mur_phys,
            has_sheets,
            has_pml,
            k_pml_red,
            m_pml_red,
            m_phys_red,
            m_pml_phys_red,
        },
        ned,
    ))
}

#[allow(clippy::too_many_arguments)]
fn build_driven_solution_mixed(
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    simplex: &Discretisation,
    proj: &MixedToSimplexProjector,
    asm: &DrivenAssembly,
    freq_hz: f64,
    x_red: &[C],
    label: String,
) -> Solution {
    let x = asm.red.expand(x_red);
    let isos: Vec<IsoparametricCell> = (0..mixed.n_cells()).map(|ci| IsoparametricCell::from_mesh(mixed, ci)).collect();
    let nsc = simplex.n_cells();
    let nloc = 10usize;
    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut cu = Vec::new();

    let sampled = proj.sample_simplex_cell_nodes(mixed, 2, |mc, xi| {
        eval_cell_basis_3d(mixed, ned, mc, &isos[mc], xi, &mut scratch, &mut w, &mut cu);
        let cdofs = ned.cell_dofs(mc);
        let mut val = [c(0.0, 0.0); 3];
        let mut cval = [c(0.0, 0.0); 3];
        for (s, &gd) in cdofs.iter().enumerate() {
            let coef = x[gd];
            for kk in 0..3 {
                val[kk] += coef * w[s][kk];
                cval[kk] += coef * cu[s][kk];
            }
        }
        (val, cval)
    });
    let e: Vec<[C; 3]> = sampled.iter().map(|p| p.0).collect();
    let curl_e: Vec<[C; 3]> = sampled.iter().map(|p| p.1).collect();

    let big =
        e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |vals: Vec<f64>| node_average(simplex, &vals.iter().map(|v| v * s).collect::<Vec<f64>>(), 1, 2);

    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
    ];

    let mut cells = CellFields::new(simplex, 2);
    for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
        let z: Vec<C> = e.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "V/m", z.iter().map(|v| v.re).collect(), z.iter().map(|v| v.im).collect());
    }

    let omega = (2.0 * PI * freq_hz).max(1e-300);
    let h_prefactor = c(0.0, 1.0 / (omega * MU0));
    let mut h_field = vec![[c(0.0, 0.0); 3]; nsc * nloc];
    for sc in 0..nsc {
        let mc = proj.parent_cell(sc);
        let mui = c(1.0 / asm.mur_phys[mc], 0.0);
        for a in 0..nloc {
            let idx = sc * nloc + a;
            let ce = curl_e[idx];
            h_field[idx] = [ce[0] * mui * h_prefactor, ce[1] * mui * h_prefactor, ce[2] * mui * h_prefactor];
        }
    }
    for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
        let z: Vec<C> = h_field.iter().map(|v| v[k]).collect();
        cells.push_complex(*name, "A/m", z.iter().map(|v| v.re).collect(), z.iter().map(|v| v.im).collect());
    }

    Solution { label, fields, cell_fields: cell_fields_of(&cells), value: Some(freq_hz), eigenvalue: None }
}
