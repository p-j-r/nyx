//! Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals, periodic
//! waveguide cells): the stiffness and mass matrices are reduced by the periodic constraints and the
//! complex Hermitian eigenproblem is solved with the complex shift-and-invert Arnoldi solver.

use crate::csparse::{c, eigs_complex, CCsr, C};
use crate::fem::free_dofs;
use crate::fields::node_average;
use crate::mesh::SimplexMesh;
use crate::periodic::{congruence, match_vertices, nodal_relations, Reduction};
use crate::problem::*;
use crate::solid::Setup;
use crate::sparse::Csr;
use std::time::Instant;

fn to_complex(a: &Csr) -> CCsr {
    CCsr { nrows: a.nrows, ncols: a.ncols, indptr: a.indptr.clone(), indices: a.indices.clone(), data: a.data.iter().map(|&v| c(v, 0.0)).collect() }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_floquet(mesh: &SimplexMesh, spec: &ProblemSpec, st: &Setup, k: &Csr, m: &Csr, n_modes: usize, f_search: f64, t0: Instant) -> Result<SolveResult, String> {
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
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    let pairs = eigs_complex(&kc, &mc, c(sigma, 0.0), n_modes, 1e-9, 500)?;
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
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, &(freq, i)) in modes.iter().enumerate() {
        let u: Vec<C> = red.expand(&pairs.vectors[i]);
        let big = u.iter().cloned().fold(c(0.0, 0.0), |mx, z| if z.norm() > mx.norm() { z } else { mx });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let mag: Vec<f64> = (0..n1).map(|a| (0..d).map(|k| u[k * n1 + a].norm_sqr()).sum::<f64>().sqrt() * scale).collect();
        let mut fields = vec![Field { name: "|u|".into(), unit: "a.u.".into(), values: mag }];
        let mut disp = vec![0.0; n1 * d];
        for (comp, ax) in ["x", "y", "z"][..d].iter().enumerate() {
            let v: Vec<f64> = (0..n1).map(|a| (u[comp * n1 + a] * rot * scale).re).collect();
            for a in 0..n1 {
                disp[a * d + comp] = v[a];
            }
            fields.push(Field { name: format!("Re u{ax}"), unit: "a.u.".into(), values: v });
        }
        solutions.push(Solution { label: format!("mode {}: {:.6} Hz", j + 1, freq), fields, disp: Some(disp), value: Some(freq), ..Default::default() });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
    }
    let _ = (node_average, free_dofs);
    let log = vec![format!("{} Floquet-periodic pair(s): {n} dofs, {nred} after constraints; complex Hermitian eigenproblem, Krylov dimension {}{}", spec.periodic.len(), pairs.dim, if pairs.converged { "" } else { " (not every pair converged; those are omitted)" })];
    Ok(SolveResult { order, dim: d, solutions, scalars, log, seconds: t0.elapsed().as_secs_f64(), dofs: n })
}
