//! Low-reflecting boundaries for solid mechanics (low-reflecting boundary): a boundary
//! traction that absorbs outgoing pressure and shear waves,
//!
//! ```text
//!   t = -rho c_p (v . n) n - rho c_s (v - (v . n) n),      v = i omega u,
//! ```
//!
//! with the wave speeds taken from the (Voigt-averaged) stiffness of the adjacent cell. The
//! eigenfrequency problem becomes quadratic in `omega`: `(K + i omega C - omega^2 M) u = 0` with
//! `C = int N^T D N dS`, `D = rho (c_s I + (c_p - c_s) n n^T)`. It is solved as the linear pencil
//! `[0 I; -K -C] z = lambda [I 0; 0 M] z` with `lambda = i omega` and `z = (u, lambda u)` by complex
//! shift-and-invert Arnoldi; the frequencies are complex and give the quality factor.

use crate::csparse::{c, eigs_complex, CCoo, CCsr, C};
use crate::fem::*;
use crate::material::M6;
use crate::mesh::SimplexMesh;
use crate::sparse::{Coo, Csr};
use std::collections::HashMap;

/// Pressure and shear speeds of an anisotropic stiffness (Voigt averages) and density.
pub fn voigt_speeds(c6: &M6, rho: f64) -> (f64, f64) {
    let bulk = (c6[0][0] + c6[1][1] + c6[2][2] + 2.0 * (c6[0][1] + c6[0][2] + c6[1][2])) / 9.0;
    let shear = (c6[0][0] + c6[1][1] + c6[2][2] - c6[0][1] - c6[0][2] - c6[1][2] + 3.0 * (c6[3][3] + c6[4][4] + c6[5][5])) / 15.0;
    (((bulk + 4.0 * shear / 3.0) / rho).sqrt(), (shear / rho).sqrt())
}

/// The real symmetric damping matrix `C` over `ncomp` component blocks of the boundaries `labels`.
/// For a 2D mesh with `ncomp = 3` (out-of-plane wave number) the normal lies in the plane.
pub fn damping_matrix(mesh: &SimplexMesh, c6: &[M6], rho: &[f64], labels: &[i32], order: usize, ncomp: usize) -> Csr {
    let d = mesh.dim;
    let n1 = mesh.n_dofs(order);
    let fd = d - 1;
    let rule = facet_rule(fd);
    let nloc = if order == 1 { d } else { d + d * (d - 1) / 2 };
    let mut me = vec![0.0; nloc * nloc];
    for (q, lam) in rule.lam.iter().enumerate() {
        let (nv, _) = shape(lam, order);
        for a in 0..nloc {
            for b in 0..nloc {
                me[a * nloc + b] += rule.w[q] * nv[a] * nv[b];
            }
        }
    }
    // the cell on each boundary facet
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
    for f in facets_with(mesh, labels) {
        let fv = mesh.facet(f);
        let ci = cell_of[&key(fv)];
        let (cp, cs) = voigt_speeds(&c6[ci], rho[ci]);
        // unit normal
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
        mesh.facet_dofs(f, order, &mut dofs);
        for i in 0..ncomp {
            for j in 0..ncomp {
                let dij = rho[ci] * (if i == j { cs } else { 0.0 } + (cp - cs) * nn[i] * nn[j]);
                if dij == 0.0 {
                    continue;
                }
                let ke: Vec<f64> = me.iter().map(|m| m * area * dij).collect();
                let (ri, cj): (Vec<usize>, Vec<usize>) = (dofs.iter().map(|&x| x + i * n1).collect(), dofs.iter().map(|&x| x + j * n1).collect());
                coo.push_block(&ri, &cj, &ke);
            }
        }
    }
    coo.to_csr()
}

fn to_complex(a: &Csr) -> CCsr {
    CCsr { nrows: a.nrows, ncols: a.ncols, indptr: a.indptr.clone(), indices: a.indices.clone(), data: a.data.iter().map(|&v| c(v, 0.0)).collect() }
}

/// Eigenvalues, eigenvectors, whether all converged, and the Krylov dimension used.
pub type QuadraticModes = (Vec<C>, Vec<Vec<C>>, bool, usize);

/// Eigenpairs `omega` (complex, rad/s) near `omega_s` of `(K + i omega C - omega^2 M) u = 0` over the free
/// dofs; also returns the displacement vectors, whether all converged and the Krylov dimension.
pub fn quadratic_modes(kf: &CCsr, mf: &CCsr, cf: &Csr, omega_s: f64, nev: usize) -> Result<QuadraticModes, String> {
    let n = kf.nrows;
    let cc = to_complex(cf);
    // A = [[0, I], [-K, -C]], B = [[I, 0], [0, M]]
    let mut a = CCoo::new(2 * n, 2 * n);
    let mut b = CCoo::new(2 * n, 2 * n);
    for i in 0..n {
        a.rows.push(i as u32);
        a.cols.push((n + i) as u32);
        a.vals.push(c(1.0, 0.0));
        b.rows.push(i as u32);
        b.cols.push(i as u32);
        b.vals.push(c(1.0, 0.0));
    }
    for r in 0..n {
        for p in kf.indptr[r]..kf.indptr[r + 1] {
            a.rows.push((n + r) as u32);
            a.cols.push(kf.indices[p] as u32);
            a.vals.push(-kf.data[p]);
        }
        for p in cc.indptr[r]..cc.indptr[r + 1] {
            a.rows.push((n + r) as u32);
            a.cols.push((n + cc.indices[p]) as u32);
            a.vals.push(-cc.data[p]);
        }
        for p in mf.indptr[r]..mf.indptr[r + 1] {
            b.rows.push((n + r) as u32);
            b.cols.push((n + mf.indices[p]) as u32);
            b.vals.push(mf.data[p]);
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
pub(crate) fn solve_lrb(mesh: &SimplexMesh, spec: &crate::problem::ProblemSpec, st: &crate::solid::Setup, k: &Csr, m: &Csr, rho_cells: &[f64], labels: &[i32], n_modes: usize, f_search: f64, t0: std::time::Instant) -> Result<crate::problem::SolveResult, String> {
    use crate::problem::*;
    if f_search <= 0.0 {
        return Err("low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into());
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
    let cfull = damping_matrix(mesh, &st.c6, rho_cells, labels, spec.order, d);
    let cf = cfull.submatrix(&free, &free);
    let (omegas, vecs, converged, dim) = quadratic_modes(&kf, &mf, &cf, 2.0 * std::f64::consts::PI * f_search, n_modes)?;
    if omegas.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let mut order_idx: Vec<usize> = (0..omegas.len()).collect();
    order_idx.sort_by(|&a, &b| omegas[a].re.partial_cmp(&omegas[b].re).unwrap());
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, &i) in order_idx.iter().enumerate() {
        let f = omegas[i] / (2.0 * std::f64::consts::PI);
        let mut u = vec![c(0.0, 0.0); n];
        for (a, &dof) in free.iter().enumerate() {
            u[dof] = vecs[i][a];
        }
        let big = u.iter().cloned().fold(c(0.0, 0.0), |mx, z| if z.norm() > mx.norm() { z } else { mx });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let mag: Vec<f64> = (0..n1).map(|a| (0..d).map(|kk| u[kk * n1 + a].norm_sqr()).sum::<f64>().sqrt() * scale).collect();
        let mut fields = vec![Field { name: "|u|".into(), unit: "a.u.".into(), values: mag }];
        let mut disp = vec![0.0; n1 * d];
        for (comp, ax) in ["x", "y", "z"][..d].iter().enumerate() {
            let v: Vec<f64> = (0..n1).map(|a| (u[comp * n1 + a] * rot * scale).re).collect();
            for a in 0..n1 {
                disp[a * d + comp] = v[a];
            }
            fields.push(Field { name: format!("Re u{ax}"), unit: "a.u.".into(), values: v });
        }
        let q = if f.im.abs() > 0.0 { f.re / (2.0 * f.im.abs()) } else { f64::INFINITY };
        solutions.push(Solution { label: format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, q), fields, disp: Some(disp), value: Some(f.re), ..Default::default() });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
        if q.is_finite() {
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        }
    }
    let log = vec![format!("low-reflecting boundaries: quadratic eigenproblem solved as a pencil of size {}, Krylov dimension {dim}{}", 2 * free.len(), if converged { "" } else { " (not every pair converged; those are omitted)" })];
    Ok(SolveResult { order: spec.order, dim: d, solutions, scalars, log, seconds: t0.elapsed().as_secs_f64(), dofs: n })
}
