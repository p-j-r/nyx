//! Eigenfrequencies of elastic bodies with perfectly matched layers (complex coordinate
//! stretching), giving complex frequencies and quality factors (anchor loss into a substrate).
//!
//! With stretch factors `s = (sx, sy, sz)` per cell, every derivative `d/dx_k` becomes
//! `(1/s_k) d/dx_k` and every volume integral gains the factor `det(diag(s))`:
//!
//!   K = int B(s)^T C B(s) det(s),    M = int rho N N det(s),    K u = omega^2 M u,
//!
//! with complex symmetric K and M. The frequency is `f = omega / (2 pi)` (complex) and the
//! quality factor `Q = Re f / (2 |Im f|)`.

use crate::solvers::csr::{c, CCoo, CCsr, C};
use crate::solvers::driver::eigs_complex;
use crate::ems_pml::cell_stretches;
use crate::fem::*;
use crate::mesh::SimplexMesh;
use crate::problem::*;
use crate::sms::{d_matrix, dof_layout, setup};
use std::time::Instant;

/// Complex strain-displacement matrix from complex shape gradients `g[a*d + k]`.
fn b_matrix_c(g: &[C], nloc: usize, d: usize, b: &mut [C]) {
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

pub fn solve_solid_pml(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_solid_pml(mesh, spec)?;
    let solved = solve_pml_modes(&asm)?;
    let r = pml_fields(mesh, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

/// Stage 1 of the solid kernel with perfectly matched layers: the complex stiffness and the mass (weighted by the
/// stretch determinant) on the free dofs.
pub struct PmlAssembly {
    pub order: usize,
    pub d: usize,
    pub n: usize,
    pub kff: CCsr,
    pub mff: CCsr,
    pub free: Vec<usize>,
    pub n_modes: usize,
    pub f_search: f64,
    pub has_fixed: bool,
}

/// Stage 2 output: eigenvalue index and complex frequency (Hz) of every mode sorted by real part, the free-dof vectors.
pub struct PmlSolve {
    pub freqs: Vec<(usize, C)>,
    pub vectors: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
}

pub fn assemble_solid_pml(mesh: &SimplexMesh, spec: &ProblemSpec) -> Result<PmlAssembly, String> {
    let order = spec.order;
    let st = setup(mesh, spec)?;
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let Study::Eigenfrequency { n_modes, f_search } = spec.study else {
        return Err("PML needs the Eigenfrequency study".into());
    };
    let thk = if d == 2 { spec.thickness } else { 1.0 };
    let rho = spec.materials.cell_values(mesh, "rho", None)?;
    let stretch = cell_stretches(mesh, &spec.materials.pml);

    let re = RefElement::new(d, order, if order > 1 { 2 * (order - 1) } else { 1 });
    let rm = RefElement::new(d, order, 2 * order);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = d * nloc;
    let ns = if d == 2 { 3 } else { 6 };
    let mut kc = CCoo::new(n, n);
    let mut mc = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * d];
    let mut gc = vec![c(0.0, 0.0); nloc * d];
    let mut b = vec![c(0.0, 0.0); ns * nd];
    let mut db = vec![c(0.0, 0.0); ns * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    let mut mloc = vec![0.0; nloc * nloc];
    for q in 0..rm.w.len() {
        for a in 0..nloc {
            for bb in 0..nloc {
                mloc[a * nloc + bb] += rm.w[q] * rm.n[q][a] * rm.n[q][bb];
            }
        }
    }
    for ci in 0..mesh.n_cells() {
        let s = stretch[ci];
        let det = if d == 2 { s[0] * s[1] } else { s[0] * s[1] * s[2] };
        let dm = d_matrix(&st.c6[ci], d, spec.plane);
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            for a in 0..nloc {
                for k in 0..d {
                    gc[a * d + k] = c(grad[a * d + k], 0.0) / s[k];
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
        let sc = det * (geo.vol[ci] * thk);
        ke.iter_mut().for_each(|x| *x *= sc);
        dof_layout(mesh, order, ci, &mut scratch, &mut full);
        kc.push_block(&full, &full, &ke);
        // mass: one scalar block per displacement component
        mesh.cell_dofs(ci, order, &mut scratch);
        let mrow: Vec<C> = mloc.iter().map(|m| det * (rho[ci] * thk * geo.vol[ci] * m)).collect();
        for comp in 0..d {
            let dofs: Vec<usize> = scratch.iter().map(|&x| x + comp * n1).collect();
            mc.push_block(&dofs, &dofs, &mrow);
        }
    }
    let (kfull, mfull) = (kc.to_csr(), mc.to_csr());
    let free = free_dofs(n, &st.fixed);
    if n_modes < 1 || n_modes + 1 >= free.len() {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kff = kfull.submatrix(&free, &free);
    let mff = mfull.submatrix(&free, &free);
    Ok(PmlAssembly { order, d, n, kff, mff, free, n_modes, f_search, has_fixed: !st.fixed.is_empty() })
}

/// Stage 2: the complex eigenpairs nearest the search frequency.
pub fn solve_pml_modes(asm: &PmlAssembly) -> Result<PmlSolve, String> {
    let sigma = if asm.f_search > 0.0 { (2.0 * std::f64::consts::PI * asm.f_search).powi(2) } else { -1.0 };
    let pairs = eigs_complex(&asm.kff, &asm.mff, c(sigma, 0.0), asm.n_modes, 1e-9, 600)?;
    // frequencies, sorted by real part
    let mut freqs: Vec<(usize, C)> = pairs.values.iter().enumerate().map(|(i, l)| {
        let mut f = l.sqrt() / (2.0 * std::f64::consts::PI);
        if f.re < 0.0 {
            f = -f;
        }
        (i, f)
    }).collect();
    freqs.sort_by(|a, b| a.1.re.partial_cmp(&b.1.re).unwrap());
    Ok(PmlSolve { freqs, vectors: pairs.vectors.clone(), dim: pairs.dim, converged: pairs.converged })
}

/// Stage 3: plot fields and per-cell complex displacement of every mode, the frequencies and Q.
pub fn pml_fields(mesh: &SimplexMesh, asm: &PmlAssembly, solved: &PmlSolve) -> Result<SolveResult, String> {
    let (order, d, n) = (asm.order, asm.d, asm.n);
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, (i, f)) in solved.freqs.iter().enumerate() {
        let mut u = vec![c(0.0, 0.0); n];
        for (a, &dof) in asm.free.iter().enumerate() {
            u[dof] = solved.vectors[*i][a];
        }
        let (fields, disp, cells) = crate::sms::complex_mode_fields(mesh, order, &u);
        let q = if f.im.abs() > 0.0 { f.re / (2.0 * f.im.abs()) } else { f64::INFINITY };
        solutions.push(Solution { label: format!("mode {}: {:.6e} {:+.3e}i Hz", j + 1, f.re, f.im), fields, cell_fields: Some(cells), disp: Some(disp), value: Some(f.re), ..Default::default() });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
        if q.is_finite() {
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        }
    }
    let mut log = vec![format!("complex eigenproblem with PML, Krylov dimension {}{}", solved.dim, if solved.converged { "" } else { " (not fully converged)" })];
    if !asm.has_fixed {
        log.push("no fixed boundary: the lowest modes are rigid-body motions".into());
    }
    Ok(SolveResult { order, dim: d, solutions, scalars, log, dofs: n, ..Default::default() })
}
