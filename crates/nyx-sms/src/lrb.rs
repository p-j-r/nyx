// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Low-reflecting boundaries for solid mechanics: a boundary traction that absorbs outgoing
//! pressure and shear waves,
//!
//! ```text
//!   t = -rho c_p (v . n) n - rho c_s (v - (v . n) n),      v = i omega u,
//! ```
//!
//! with the wave speeds taken from the (Voigt-averaged) stiffness of the adjacent cell. The
//! eigenfrequency problem becomes quadratic in `omega`: `(K + i omega C - omega^2 M) u = 0` with
//! `C = int N^T D N dS`, `D = rho (c_s I + (c_p - c_s) n n^T)`. It is solved as the linear pencil
//! `[0 I; -K -C] z = lambda [I 0; 0 M] z` with `lambda = i omega` and `z = (u, lambda u)` by complex
//! shift-and-invert Arnoldi; the frequencies are complex and give the signed quality factor
//! `Q = Re f / (2 Im f)`.

use crate::problem::Problem;
use femcore::fem::*;
use femcore::fem_mixed::{eval_facet_jacobian, facet_lagrange_shape, facet_rule_mixed, MixedLagrangeDofMap};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::csr::{Coo, RCsr as Csr};
use femcore::linalg::driver::eigs_complex;
use femcore::material::M6;
use femcore::mesh::cells::canonical_quad_face;
use femcore::mesh::mixed::{
    cell_local_edges, cell_local_quad_faces, cell_local_tri_faces, sorted_tri_face, FacetType, MixedDiscretisation,
};
use nyx_physics::{cell_fields_of, Discretisation, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult};
use std::collections::HashMap;

/// Relative imaginary-frequency threshold below which `f.im < 0` flags an unphysical amplifying mode.
pub(crate) const AMPLIFYING_IM_TOL: f64 = 1e-9;

/// Pressure and shear speeds of an anisotropic stiffness (Voigt averages) and density.
pub fn voigt_speeds(c6: &M6, rho: f64) -> (f64, f64) {
    let bulk = (c6[0][0] + c6[1][1] + c6[2][2] + 2.0 * (c6[0][1] + c6[0][2] + c6[1][2])) / 9.0;
    let shear = (c6[0][0] + c6[1][1] + c6[2][2] - c6[0][1] - c6[0][2] - c6[1][2]
        + 3.0 * (c6[3][3] + c6[4][4] + c6[5][5]))
        / 15.0;
    (((bulk + 4.0 * shear / 3.0) / rho).sqrt(), (shear / rho).sqrt())
}

/// Voigt-averaged `(bulk_K, shear_G)` moduli of a `6x6` stiffness matrix `c6` without dividing by `rho`.
pub(crate) fn voigt_moduli(c6: &M6) -> (f64, f64) {
    let bulk = (c6[0][0] + c6[1][1] + c6[2][2] + 2.0 * (c6[0][1] + c6[0][2] + c6[1][2])) / 9.0;
    let shear = (c6[0][0] + c6[1][1] + c6[2][2] - c6[0][1] - c6[0][2] - c6[1][2]
        + 3.0 * (c6[3][3] + c6[4][4] + c6[5][5]))
        / 15.0;
    (bulk.max(0.0), shear.max(0.0))
}

/// The real symmetric damping matrix `C` over `ncomp` component blocks of the boundaries `labels`.
/// For a 2D mesh with `ncomp = 3` (out-of-plane wave number) the normal lies in the plane.
pub fn damping_matrix(
    mesh: &Discretisation,
    c6: &[M6],
    rho: &[f64],
    labels: &[i32],
    order: usize,
    ncomp: usize,
) -> Csr {
    damping_matrix_impl(mesh, c6, rho, labels, order, ncomp, false)
}

/// Axisymmetric real symmetric damping matrix `C` with the radial Jacobian weight `r`
/// (quadrature-point `int N^T D N r ds`, matching `axisym`'s `r dr dz` weak form, `R6b`).
pub fn damping_matrix_axisym(
    mesh: &Discretisation,
    c6: &[M6],
    rho: &[f64],
    labels: &[i32],
    order: usize,
    ncomp: usize,
) -> Csr {
    damping_matrix_impl(mesh, c6, rho, labels, order, ncomp, true)
}

fn damping_matrix_impl(
    mesh: &Discretisation,
    c6: &[M6],
    rho: &[f64],
    labels: &[i32],
    order: usize,
    ncomp: usize,
    radial_weight: bool,
) -> Csr {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let fd = d - 1;
    let rule = facet_rule(fd);
    let nloc = if order == 1 { d } else { d + d * (d - 1) / 2 };
    let shapes_q: Vec<Vec<f64>> = rule.lam.iter().map(|lam| shape(lam, order).0).collect();
    let mut me_const = vec![0.0; nloc * nloc];
    if !radial_weight {
        for (q, nv) in shapes_q.iter().enumerate() {
            for a in 0..nloc {
                for b in 0..nloc {
                    me_const[a * nloc + b] += rule.w[q] * nv[a] * nv[b];
                }
            }
        }
    }
    let key = |v: &[u32]| -> [u32; 3] {
        let mut k = [u32::MAX; 3];
        k[..v.len()].copy_from_slice(v);
        k[..v.len()].sort_unstable();
        k
    };
    let wanted: std::collections::HashSet<[u32; 3]> = facets_with(mesh, labels).map(|f| key(mesh.facet(f))).collect();
    let mut cell_of: HashMap<[u32; 3], usize> = HashMap::new();
    for ci in 0..mesh.n_cells() {
        let v = mesh.cell(ci);
        for skip in 0..=d {
            let f: Vec<u32> = (0..=d).filter(|&i| i != skip).map(|i| v[i]).collect();
            let k = key(&f);
            if wanted.contains(&k) {
                cell_of.entry(k).or_insert(ci);
            }
        }
    }
    let n = ncomp * n1;
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    let mut me_r = vec![0.0; nloc * nloc];
    for f in facets_with(mesh, labels) {
        let fv = mesh.facet(f);
        let ci = cell_of[&key(fv)];
        if rho[ci] <= 0.0 {
            continue;
        }
        let (cp, cs) = voigt_speeds(&c6[ci], rho[ci]);
        let p = |i: usize| -> [f64; 3] {
            let q = mesh.point(fv[i] as usize);
            [q[0], q[1], if d == 3 { q[2] } else { 0.0 }]
        };
        let nrm = if d == 2 {
            let (a, b) = (p(0), p(1));
            [b[1] - a[1], -(b[0] - a[0]), 0.0]
        } else {
            let (a, b, cc) = (p(0), p(1), p(2));
            let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [cc[0] - a[0], cc[1] - a[1], cc[2] - a[2]]);
            [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        };
        let len = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt().max(1e-300);
        let nn = [nrm[0] / len, nrm[1] / len, nrm[2] / len];
        let area = facet_measure(mesh, f);
        let me: &[f64] = if radial_weight {
            me_r.fill(0.0);
            for (q, lam) in rule.lam.iter().enumerate() {
                let rq: f64 = (0..fv.len()).map(|i| lam[i] * mesh.point(fv[i] as usize)[0]).sum();
                let w = rule.w[q] * rq;
                let nv = &shapes_q[q];
                for a in 0..nloc {
                    for b in 0..nloc {
                        me_r[a * nloc + b] += w * nv[a] * nv[b];
                    }
                }
            }
            &me_r
        } else {
            &me_const
        };
        mesh.facet_dofs(f, order, &mut dofs);
        for i in 0..ncomp {
            for j in 0..ncomp {
                let dij = rho[ci] * (if i == j { cs } else { 0.0 } + (cp - cs) * nn[i] * nn[j]);
                if dij == 0.0 {
                    continue;
                }
                let ke: Vec<f64> = me.iter().map(|m| m * area * dij).collect();
                let (ri, cj): (Vec<usize>, Vec<usize>) =
                    (dofs.iter().map(|&x| x + i * n1).collect(), dofs.iter().map(|&x| x + j * n1).collect());
                coo.push_block(&ri, &cj, &ke);
            }
        }
    }
    coo.to_csr()
}

/// Real symmetric low-reflecting damping matrix `C` assembled natively on a [`MixedDiscretisation`]
/// (`Edge2`, `Tri3`, and `Quad4` facets, straight or curved).
pub(crate) fn damping_matrix_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    rho: &[f64],
    labels: &[i32],
    ncomp: usize,
    thickness_2d: f64,
) -> Csr {
    damping_matrix_mixed_weighted(mixed, dof_map, c6, rho, labels, ncomp, thickness_2d, false)
}

/// Low-reflecting damping matrix on a [`MixedDiscretisation`] with optional radial Jacobian weight `r` (`R6b`).
pub(crate) fn damping_matrix_mixed_weighted(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    c6: &[M6],
    rho: &[f64],
    labels: &[i32],
    ncomp: usize,
    thickness_2d: f64,
    radial_weight: bool,
) -> Csr {
    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let n = ncomp * n1;
    if labels.is_empty() {
        return Csr::zeros(n, n);
    }

    // Build lookup from canonical facet vertices to an active adjacent cell (`rho[ci] > 0`).
    let mut edge_cell: HashMap<(u32, u32), usize> = HashMap::new();
    let mut tri_cell: HashMap<[u32; 3], usize> = HashMap::new();
    let mut quad_cell: HashMap<[u32; 4], usize> = HashMap::new();
    for ci in 0..mixed.n_cells() {
        if rho[ci] <= 0.0 {
            continue;
        }
        let vs = mixed.cell(ci);
        let ct = mixed.cell_type(ci);
        if d == 2 {
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

    let deg = 2 * dof_map.space.order() + if radial_weight { 4 } else { 2 };
    let mut coo = Coo::new(n, n);
    let mut dofs = Vec::new();
    for f in 0..mixed.n_facets() {
        if !labels.contains(&mixed.facet_labels[f]) {
            continue;
        }
        let ft = mixed.facet_type(f);
        let fv = mixed.facet(f);
        let ci_opt = match ft {
            FacetType::Edge2 => edge_cell.get(&(fv[0].min(fv[1]), fv[0].max(fv[1]))).copied(),
            FacetType::Tri3 => tri_cell.get(&sorted_tri_face([fv[0], fv[1], fv[2]])).copied(),
            FacetType::Quad4 => quad_cell.get(&canonical_quad_face([fv[0], fv[1], fv[2], fv[3]])).copied(),
        };
        let Some(ci) = ci_opt else { continue };
        let (cp, cs) = voigt_speeds(&c6[ci], rho[ci]);
        dof_map.facet_dofs(mixed, f, &mut dofs);
        let nloc = dofs.len();
        let rule = facet_rule_mixed(ft, deg);
        let thk = if d == 2 && ncomp == 2 { thickness_2d } else { 1.0 };

        let mut blocks = vec![vec![0.0; nloc * nloc]; ncomp * ncomp];
        for (xi, &wq) in rule.pts.iter().zip(&rule.w) {
            let (xq, ds, n_hat) = eval_facet_jacobian(mixed, f, xi);
            if ds <= 0.0 {
                continue;
            }
            let (nv, _) = facet_lagrange_shape(ft, dof_map.space, xi);
            let rq = if radial_weight { xq[0] } else { 1.0 };
            let w_ds = wq * ds * thk * rq;
            for i in 0..ncomp {
                for j in 0..ncomp {
                    let ni = if i < d { n_hat[i] } else { 0.0 };
                    let nj = if j < d { n_hat[j] } else { 0.0 };
                    let dij = rho[ci] * (if i == j { cs } else { 0.0 } + (cp - cs) * ni * nj);
                    if dij == 0.0 {
                        continue;
                    }
                    let blk = &mut blocks[i * ncomp + j];
                    for a in 0..nloc {
                        let wa = w_ds * dij * nv[a];
                        for b in 0..nloc {
                            blk[a * nloc + b] += wa * nv[b];
                        }
                    }
                }
            }
        }
        for i in 0..ncomp {
            let ri: Vec<usize> = dofs.iter().map(|&x| x + i * n1).collect();
            for j in 0..ncomp {
                let blk = &blocks[i * ncomp + j];
                if blk.iter().any(|&v| v != 0.0) {
                    let cj: Vec<usize> = dofs.iter().map(|&x| x + j * n1).collect();
                    coo.push_block(&ri, &cj, blk);
                }
            }
        }
    }
    coo.to_csr()
}

fn to_complex(a: &Csr) -> CCsr {
    a.to_complex()
}

/// Eigenvalues, eigenvectors, whether all converged, and the Krylov dimension used.
pub type QuadraticModes = (Vec<C>, Vec<Vec<C>>, bool, usize);

/// Eigenpairs `omega` (complex, rad/s) near `omega_s` of `(K + i omega C - omega^2 M) u = 0` over the free
/// dofs; also returns the displacement vectors, whether all converged and the Krylov dimension.
pub fn quadratic_modes(kf: &CCsr, mf: &CCsr, cf: &Csr, omega_s: f64, nev: usize) -> Result<QuadraticModes, String> {
    quadratic_modes_complex(kf, mf, &to_complex(cf), omega_s, nev)
}

/// Complex damping variant of [`quadratic_modes`] (supports congruence-reduced complex LRB matrices).
pub(crate) fn quadratic_modes_complex(
    kf: &CCsr,
    mf: &CCsr,
    cc: &CCsr,
    omega_s: f64,
    nev: usize,
) -> Result<QuadraticModes, String> {
    let n = kf.nrows();
    // A = [[0, I], [-K, -C]], B = [[I, 0], [0, M]]
    let mut a = CCoo::new(2 * n, 2 * n);
    let mut b = CCoo::new(2 * n, 2 * n);
    for i in 0..n {
        a.push(i, n + i, c(1.0, 0.0));
        b.push(i, i, c(1.0, 0.0));
    }
    for r in 0..n {
        for p in kf.indptr()[r]..kf.indptr()[r + 1] {
            a.push(n + r, kf.indices()[p] as usize, -kf.values()[p]);
        }
        for p in cc.indptr()[r]..cc.indptr()[r + 1] {
            a.push(n + r, n + cc.indices()[p] as usize, -cc.values()[p]);
        }
        for p in mf.indptr()[r]..mf.indptr()[r + 1] {
            b.push(n + r, n + mf.indices()[p] as usize, mf.values()[p]);
        }
    }
    let (ac, bc) = (a.to_csr(), b.to_csr());
    // lambda = i omega; the shift sits on the positive imaginary axis
    let pairs = eigs_complex(&ac, &bc, c(0.0, omega_s), nev, 1e-9, 600)?;
    let mut omegas = vec![];
    let mut vecs = vec![];
    for (i, l) in pairs.values.iter().enumerate() {
        if pairs.residuals.get(i).is_some_and(|r| *r > 1e-6) {
            continue;
        }
        let omega = c(0.0, -1.0) * *l; // omega = -i lambda
        if omega.re < 0.0 {
            continue;
        }
        omegas.push(omega);
        vecs.push(pairs.vectors[i][..n].to_vec());
    }
    Ok((omegas, vecs, pairs.converged, pairs.dim))
}

/// Eigenfrequencies with low-reflecting boundaries for the real stiffness and mass matrices of
/// `solid.rs` (2D plane problems and 3D).
#[allow(clippy::too_many_arguments)]
pub fn solve_lrb(
    mesh: &Discretisation,
    spec: &Problem,
    st: &crate::solid::Setup,
    k: &Csr,
    m: &Csr,
    rho_cells: &[f64],
    labels: &[i32],
    n_modes: usize,
    f_search: f64,
    t0: std::time::Instant,
) -> Result<SolveResult, String> {
    let asm = assemble_lrb(mesh, spec, st, k, m, rho_cells, labels, n_modes, f_search)?;
    let solved = solve_lrb_modes(&asm)?;
    let r = lrb_fields(mesh, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_lrb_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    c6: &[M6],
    fixed: &[usize],
    k: &Csr,
    m: &Csr,
    rho_cells: &[f64],
    labels: &[i32],
    n_modes: usize,
    f_search: f64,
    t0: std::time::Instant,
) -> Result<SolveResult, String> {
    if f_search <= 0.0 {
        return Err(
            "low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into()
        );
    }
    if !spec.periodic.is_empty() {
        return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
    }
    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    let free = free_dofs(n, fixed);
    if n_modes < 1 || n_modes + 1 >= free.len() {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kf = to_complex(&k.submatrix(&free, &free));
    let mf = to_complex(&m.submatrix(&free, &free));
    let thk = if d == 2 { spec.thickness } else { 1.0 };
    let mut cfull = damping_matrix_mixed(mixed, dof_map, c6, rho_cells, labels, d, thk);
    if spec.materials.has_viscosity() {
        let eta6: Vec<M6> = crate::solid::mixed_cell_viscosity_skipping(&spec.materials, mixed, &spec.inactive)?
            .iter()
            .map(|c| femcore::section::plane_stiffness(c, d))
            .collect();
        let c_visc = crate::solid::assemble_solid_mixed(mixed, dof_map, &eta6, spec.plane, spec.thickness);
        cfull = if cfull.nnz() == 0 { c_visc } else { cfull.add_scaled(&c_visc, 1.0) };
    }
    let cf = cfull.submatrix(&free, &free);
    let asm = LrbAssembly { order: spec.order, d, n, kf, mf, cf, free, n_modes, f_search };
    let solved = solve_lrb_modes(&asm)?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;

    let mut order_idx: Vec<usize> = (0..solved.omegas.len()).collect();
    order_idx.sort_by(|&a, &b| solved.omegas[a].re.partial_cmp(&solved.omegas[b].re).unwrap());
    let mut solutions = Vec::new();
    let mut scalars = Vec::new();
    let mut unstable_count = 0usize;
    for (j, &i) in order_idx.iter().enumerate() {
        let f = solved.omegas[i] / (2.0 * std::f64::consts::PI);
        let is_amplifying = f.im < -AMPLIFYING_IM_TOL * f.re.abs().max(1e-30);
        if is_amplifying {
            unstable_count += 1;
        }
        let mut u = vec![c(0.0, 0.0); n];
        for (a, &dof) in asm.free.iter().enumerate() {
            u[dof] = solved.vectors[i][a];
        }
        let (fields, cell_fields) = crate::solid::complex_mode_fields_mixed(mixed, dof_map, &proj, asm.order, &u);
        let q = if f.im.abs() > 0.0 { f.re / (2.0 * f.im) } else { f64::INFINITY };
        let label = if is_amplifying {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]", j + 1, f.re, f.im, q)
        } else {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, q)
        };
        solutions.push(Solution { label, fields, cell_fields, value: Some(f.re), ..Default::default() });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
        if q.is_finite() {
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        }
    }
    let mut log = vec![format!(
        "low-reflecting boundaries: quadratic eigenproblem solved as a pencil of size {}, Krylov dimension {}{}",
        2 * asm.free.len(),
        solved.dim,
        if solved.converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    if unstable_count > 0 {
        scalars.push(("unstable_modes".into(), Scalar::Num(unstable_count as f64)));
        log.push(format!("WARNING: {unstable_count} amplifying mode(s) detected (Im f < -1e-9 Re f)"));
    }
    Ok(SolveResult {
        order: asm.order as u8,
        dim: d as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        ..Default::default()
    })
}

/// Stage 1 of the solid kernel with low-reflecting boundaries: stiffness, mass and the boundary damping on the free dofs
/// (the eigenproblem is quadratic in the frequency).
pub struct LrbAssembly {
    pub order: usize,
    pub d: usize,
    pub n: usize,
    pub kf: CCsr,
    pub mf: CCsr,
    pub cf: Csr,
    pub free: Vec<usize>,
    pub n_modes: usize,
    pub f_search: f64,
}

/// Stage 2 output: complex angular frequencies and the free-dof vectors of the converged modes.
pub struct LrbSolve {
    pub omegas: Vec<C>,
    pub vectors: Vec<Vec<C>>,
    pub converged: bool,
    pub dim: usize,
}

#[allow(clippy::too_many_arguments)]
pub fn assemble_lrb(
    mesh: &Discretisation,
    spec: &Problem,
    st: &crate::solid::Setup,
    k: &Csr,
    m: &Csr,
    rho_cells: &[f64],
    labels: &[i32],
    n_modes: usize,
    f_search: f64,
) -> Result<LrbAssembly, String> {
    crate::problem::check_simplex_straight(mesh)?;
    if f_search <= 0.0 {
        return Err(
            "low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into()
        );
    }
    if !spec.periodic.is_empty() {
        return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
    }
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let free = free_dofs(n, &st.fixed);
    if n_modes < 1 || n_modes + 1 >= free.len() {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kf = to_complex(&k.submatrix(&free, &free));
    let mf = to_complex(&m.submatrix(&free, &free));
    let mut cfull = if labels.is_empty() {
        Csr::zeros(n, n)
    } else {
        damping_matrix(mesh, &st.c6, rho_cells, labels, spec.order, d)
    };
    if spec.materials.has_viscosity() {
        let eta6: Vec<M6> = spec
            .materials
            .cell_viscosity_skipping(mesh, &spec.inactive)?
            .iter()
            .map(|c| femcore::section::plane_stiffness(c, d))
            .collect();
        let c_visc = crate::solid::assemble_solid(mesh, &eta6, spec.order, spec.plane, spec.thickness);
        cfull = if cfull.nnz() == 0 { c_visc } else { cfull.add_scaled(&c_visc, 1.0) };
    }
    let cf = cfull.submatrix(&free, &free);
    Ok(LrbAssembly { order: spec.order, d, n, kf, mf, cf, free, n_modes, f_search })
}

/// Stage 2: the quadratic eigenpairs nearest the search frequency.
pub fn solve_lrb_modes(asm: &LrbAssembly) -> Result<LrbSolve, String> {
    let (omegas, vectors, converged, dim) =
        quadratic_modes(&asm.kf, &asm.mf, &asm.cf, 2.0 * std::f64::consts::PI * asm.f_search, asm.n_modes)?;
    if omegas.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    Ok(LrbSolve { omegas, vectors, converged, dim })
}

/// Stage 3: plot fields and per-cell complex displacement of every mode, the frequencies and signed Q.
pub fn lrb_fields(mesh: &Discretisation, asm: &LrbAssembly, solved: &LrbSolve) -> Result<SolveResult, String> {
    let (n, d) = (asm.n, asm.d);
    let (omegas, vecs, free) = (&solved.omegas, &solved.vectors, &asm.free);
    let (converged, dim) = (solved.converged, solved.dim);
    let mut order_idx: Vec<usize> = (0..omegas.len()).collect();
    order_idx.sort_by(|&a, &b| omegas[a].re.partial_cmp(&omegas[b].re).unwrap());
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut unstable_count = 0usize;
    for (j, &i) in order_idx.iter().enumerate() {
        let f = omegas[i] / (2.0 * std::f64::consts::PI);
        let is_amplifying = f.im < -AMPLIFYING_IM_TOL * f.re.abs().max(1e-30);
        if is_amplifying {
            unstable_count += 1;
        }
        let mut u = vec![c(0.0, 0.0); n];
        for (a, &dof) in free.iter().enumerate() {
            u[dof] = vecs[i][a];
        }
        let (fields, cells) = crate::solid::complex_mode_fields(mesh, asm.order, &u);
        let q = if f.im.abs() > 0.0 { f.re / (2.0 * f.im) } else { f64::INFINITY };
        let label = if is_amplifying {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]", j + 1, f.re, f.im, q)
        } else {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, q)
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(f.re),
            ..Default::default()
        });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
        if q.is_finite() {
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        }
    }
    let mut log = vec![format!(
        "low-reflecting boundaries: quadratic eigenproblem solved as a pencil of size {}, Krylov dimension {dim}{}",
        2 * free.len(),
        if converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    if unstable_count > 0 {
        scalars.push(("unstable_modes".into(), Scalar::Num(unstable_count as f64)));
        log.push(format!("WARNING: {unstable_count} amplifying mode(s) detected (Im f < -1e-9 Re f)"));
    }
    Ok(SolveResult { order: asm.order as u8, dim: d as u8, solutions, scalars, log, dofs: n, ..Default::default() })
}
