//! Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals, periodic
//! waveguide cells): the stiffness and mass matrices are reduced by the periodic constraints and the
//! complex Hermitian eigenproblem is solved with the complex shift-and-invert Arnoldi solver.

use crate::csparse::{c, eigs_complex, CCsr, C};
use crate::mesh::SimplexMesh;
use crate::periodic::{congruence, match_vertices, nodal_relations, Reduction};
use crate::problem::*;
use crate::sms::Setup;
use crate::sparse::Csr;
use std::time::Instant;

fn to_complex(a: &Csr) -> CCsr {
    CCsr { nrows: a.nrows, ncols: a.ncols, indptr: a.indptr.clone(), indices: a.indices.clone(), data: a.data.iter().map(|&v| c(v, 0.0)).collect() }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_floquet(mesh: &SimplexMesh, spec: &ProblemSpec, st: &Setup, k: &Csr, m: &Csr, n_modes: usize, f_search: f64, t0: Instant) -> Result<SolveResult, String> {
    let asm = assemble_floquet(mesh, spec, st, k, m, n_modes, f_search)?;
    let solved = solve_floquet_modes(&asm)?;
    let r = floquet_fields(mesh, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

/// Stage 1 of the Floquet-periodic solid kernel: the stiffness and mass reduced to the independent dofs of the
/// periodic pairs (complex, with the Bloch phase) and the reduction that expands a mode back to all dofs.
pub struct FloquetAssembly {
    pub order: usize,
    pub n: usize,
    pub nred: usize,
    pub n_pairs: usize,
    pub red: Reduction,
    pub kc: CCsr,
    pub mc: CCsr,
    pub n_modes: usize,
    pub f_search: f64,
}

/// Stage 2 output: frequency and index of each converged mode (by frequency), the reduced vectors.
pub struct FloquetSolve {
    pub modes: Vec<(f64, usize)>,
    pub vectors: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
}

pub fn assemble_floquet(mesh: &SimplexMesh, spec: &ProblemSpec, st: &Setup, k: &Csr, m: &Csr, n_modes: usize, f_search: f64) -> Result<FloquetAssembly, String> {
    let (d, n1) = (st.d, st.n1);
    let n = d * n1;
    let order = spec.order;
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
    let kc = congruence(&pm, &to_complex(k));
    let mc = congruence(&pm, &to_complex(m));
    Ok(FloquetAssembly { order, n, nred, n_pairs: spec.periodic.len(), red, kc, mc, n_modes, f_search })
}

/// Stage 2: the complex Hermitian eigenpairs of the reduced problem nearest the search frequency.
pub fn solve_floquet_modes(asm: &FloquetAssembly) -> Result<FloquetSolve, String> {
    let sigma = if asm.f_search > 0.0 { (2.0 * std::f64::consts::PI * asm.f_search).powi(2) } else { -1.0 };
    let pairs = eigs_complex(&asm.kc, &asm.mc, c(sigma, 0.0), asm.n_modes, 1e-9, 500)?;
    let mut modes: Vec<(f64, usize)> = pairs
        .values
        .iter()
        .enumerate()
        .filter(|(i, _)| pairs.residuals.get(*i).is_none_or(|r| *r < 1e-6))
        .map(|(i, l)| (l.re.max(0.0).sqrt() / (2.0 * std::f64::consts::PI), i))
        .collect();
    modes.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    if modes.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    Ok(FloquetSolve { modes, vectors: pairs.vectors.clone(), dim: pairs.dim, converged: pairs.converged })
}

/// Stage 3: plot fields and per-cell complex displacement of every mode.
pub fn floquet_fields(mesh: &SimplexMesh, asm: &FloquetAssembly, solved: &FloquetSolve) -> Result<SolveResult, String> {
    let (order, n, nred) = (asm.order, asm.n, asm.nred);
    let d = mesh.dim;
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, &(freq, i)) in solved.modes.iter().enumerate() {
        let u: Vec<C> = asm.red.expand(&solved.vectors[i]);
        let (fields, disp, cells) = crate::sms::complex_mode_fields(mesh, order, &u);
        solutions.push(Solution { label: format!("mode {}: {:.6} Hz", j + 1, freq), fields, cell_fields: Some(cells), disp: Some(disp), value: Some(freq), ..Default::default() });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
    }
    let log = vec![format!("{} Floquet-periodic pair(s): {n} dofs, {nred} after constraints; complex Hermitian eigenproblem, Krylov dimension {}{}", asm.n_pairs, solved.dim, if solved.converged { "" } else { " (not every pair converged; those are omitted)" })];
    Ok(SolveResult { order, dim: d, solutions, scalars, log, dofs: n, ..Default::default() })
}
