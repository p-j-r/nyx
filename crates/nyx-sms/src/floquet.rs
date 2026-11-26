// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Eigenfrequencies of a solid with Floquet-periodic boundaries (phononic crystals, periodic
//! waveguide cells): the stiffness and mass matrices are reduced by the periodic constraints and the
//! complex Hermitian eigenproblem is solved with the complex shift-and-invert Arnoldi solver.

use crate::problem::*;
use crate::solid::Setup;
use femcore::fem_mixed::{LagrangeSpace, MixedLagrangeDofMap};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::mesh::mixed::{FacetType, MixedDiscretisation, QUAD_EDGES, TRI_EDGES};
use femcore::periodic::{congruence, match_vertices, nodal_relations, phase_factor, Reduction};
use femcore::props::Periodic;
use nyx_physics::{cell_fields_of, Discretisation, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult};
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

fn to_complex(a: &Csr) -> CCsr {
    a.to_complex()
}

/// Builds Floquet-periodic DOF relations `(dest_dof, src_dof, exp(-i phase))` directly on a
/// [`MixedDiscretisation`] (`Edge2`, `Tri3`, and `Quad4` periodic boundary facets).
pub(crate) fn mixed_nodal_relations(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    ncomp: usize,
    p: &Periodic,
) -> Result<Vec<(usize, usize, C)>, String> {
    let d = mixed.dim;
    let group_vertices = |labels: &[i32]| -> (Vec<usize>, Vec<u32>) {
        let mut facets = Vec::new();
        let mut verts = BTreeSet::new();
        for f in 0..mixed.n_facets() {
            if labels.contains(&mixed.facet_labels[f]) {
                facets.push(f);
                for &v in mixed.facet(f) {
                    verts.insert(v);
                }
            }
        }
        (facets, verts.into_iter().collect())
    };

    let (src_facets, sv) = group_vertices(&p.source);
    let (_, dv) = group_vertices(&p.dest);
    if sv.is_empty() || dv.is_empty() {
        return Err("a periodic condition needs boundary numbers on both sides that exist in the mesh".into());
    }
    if sv.len() != dv.len() {
        return Err(format!(
            "the periodic boundaries have different numbers of vertices ({} and {}): their meshes must match",
            sv.len(),
            dv.len()
        ));
    }

    let centroid = |vs: &[u32]| -> Vec<f64> {
        let mut c0 = vec![0.0; d];
        for &v in vs {
            let pt = mixed.point(v as usize);
            for k in 0..d {
                c0[k] += pt[k] / vs.len() as f64;
            }
        }
        c0
    };
    let (cs, cd) = (centroid(&sv), centroid(&dv));
    let shift: Vec<f64> = (0..d).map(|k| cd[k] - cs[k]).collect();

    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for i in 0..mixed.n_points() {
        let pt = mixed.point(i);
        for k in 0..d {
            lo[k] = lo[k].min(pt[k]);
            hi[k] = hi[k].max(pt[k]);
        }
    }
    let extent = (0..d).map(|k| hi[k] - lo[k]).fold(0.0, f64::max).max(1e-300);
    let tol = 1e-6 * extent;
    let cell = 4.0 * tol;
    let key = |pt: &[f64]| -> Vec<i64> { pt.iter().map(|x| (x / cell).floor() as i64).collect() };
    let mut grid: HashMap<Vec<i64>, Vec<u32>> = HashMap::new();
    for &v in &dv {
        grid.entry(key(mixed.point(v as usize))).or_default().push(v);
    }

    let mut offsets: Vec<Vec<i64>> = vec![vec![]];
    for _ in 0..d {
        offsets = offsets
            .into_iter()
            .flat_map(|o| {
                (-1..=1).map(move |dk| {
                    let mut n = o.clone();
                    n.push(dk);
                    n
                })
            })
            .collect();
    }

    let mut map: HashMap<u32, u32> = HashMap::new();
    for &s in &sv {
        let ps = mixed.point(s as usize);
        let target: Vec<f64> = (0..d).map(|k| ps[k] + shift[k]).collect();
        let kt = key(&target);
        let mut found = None;
        for off in &offsets {
            let kk: Vec<i64> = kt.iter().zip(off).map(|(a, b)| a + b).collect();
            if let Some(list) = grid.get(&kk) {
                for &t in list {
                    let pt = mixed.point(t as usize);
                    if (0..d).all(|k| (pt[k] - target[k]).abs() <= tol) {
                        found = Some(t);
                    }
                }
            }
        }
        match found {
            Some(t) => {
                map.insert(s, t);
            }
            None => {
                return Err(format!(
                    "no vertex on the destination boundary matches source vertex {s} shifted by the translation: the meshes of the periodic boundaries must be identical"
                ))
            }
        }
    }

    let f = phase_factor(p.phase);
    let n1 = dof_map.n_dofs;
    let nv = dof_map.n_vertices;
    let mut out = Vec::new();
    for (&s, &t) in &map {
        for k in 0..ncomp {
            out.push((t as usize + k * n1, s as usize + k * n1, f));
        }
    }
    if dof_map.space != LagrangeSpace::Order1 {
        let mut done = BTreeSet::new();
        for &fc in &src_facets {
            let fv = mixed.facet(fc);
            let local_edges: &[(usize, usize)] = match mixed.facet_type(fc) {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in local_edges {
                let (a, b) = (fv[i], fv[j]);
                let ekey = (a.min(b), a.max(b));
                if !done.insert(ekey) {
                    continue;
                }
                let es = mixed.edge_index(a, b).ok_or("source boundary facet edge missing from MixedDiscretisation")?;
                let et = mixed
                    .edge_index(map[&a], map[&b])
                    .ok_or("an edge of the source boundary has no counterpart on the destination boundary")?;
                for k in 0..ncomp {
                    out.push((nv + et + k * n1, nv + es + k * n1, f));
                }
            }
        }
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_floquet(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    k: &Csr,
    m: &Csr,
    n_modes: usize,
    f_search: f64,
    t0: Instant,
) -> Result<SolveResult, String> {
    let asm = assemble_floquet(mesh, spec, st, k, m, n_modes, f_search)?;
    let solved = solve_floquet_modes(&asm)?;
    let r = floquet_fields(mesh, &asm, &solved)?;
    Ok(SolveResult { seconds: t0.elapsed().as_secs_f64(), ..r })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_floquet_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    spec: &Problem,
    fixed: &[usize],
    k: &Csr,
    m: &Csr,
    n_modes: usize,
    f_search: f64,
    t0: Instant,
) -> Result<SolveResult, String> {
    let d = mixed.dim;
    let n1 = dof_map.n_dofs;
    let n = d * n1;
    let mut relations = Vec::new();
    for p in &spec.periodic {
        relations.extend(mixed_nodal_relations(mixed, dof_map, d, p)?);
    }
    let red = Reduction::build(n, &relations, fixed);
    let pm = red.prolong();
    let nred = red.n_reduced();
    if n_modes < 1 || n_modes + 1 >= nred {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kc = congruence(&pm, &to_complex(k));
    let mc = congruence(&pm, &to_complex(m));
    let asm =
        FloquetAssembly { order: spec.order, n, nred, n_pairs: spec.periodic.len(), red, kc, mc, n_modes, f_search };
    let solved = solve_floquet_modes(&asm)?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let mut solutions = Vec::new();
    let mut scalars = Vec::new();
    for (j, &(freq, i)) in solved.modes.iter().enumerate() {
        let u: Vec<C> = asm.red.expand(&solved.vectors[i]);
        let (fields, cell_fields) = crate::solid::complex_mode_fields_mixed(mixed, dof_map, &proj, spec.order, &u);
        solutions.push(Solution {
            label: format!("mode {}: {:.6} Hz", j + 1, freq),
            fields,
            cell_fields,
            value: Some(freq),
            ..Default::default()
        });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
    }
    let log = vec![format!(
        "{} Floquet-periodic pair(s): {n} dofs, {nred} after constraints; complex Hermitian eigenproblem, Krylov dimension {}{}",
        asm.n_pairs,
        solved.dim,
        if solved.converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    Ok(SolveResult {
        order: spec.order as u8,
        dim: d as u8,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: n,
        ..Default::default()
    })
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

pub fn assemble_floquet(
    mesh: &Discretisation,
    spec: &Problem,
    st: &Setup,
    k: &Csr,
    m: &Csr,
    n_modes: usize,
    f_search: f64,
) -> Result<FloquetAssembly, String> {
    check_simplex_straight(mesh)?;
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
pub fn floquet_fields(
    mesh: &Discretisation,
    asm: &FloquetAssembly,
    solved: &FloquetSolve,
) -> Result<SolveResult, String> {
    let (order, n, nred) = (asm.order, asm.n, asm.nred);
    let d = mesh.dim;
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, &(freq, i)) in solved.modes.iter().enumerate() {
        let u: Vec<C> = asm.red.expand(&solved.vectors[i]);
        let (fields, cells) = crate::solid::complex_mode_fields(mesh, order, &u);
        solutions.push(Solution {
            label: format!("mode {}: {:.6} Hz", j + 1, freq),
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(freq),
            ..Default::default()
        });
        scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
    }
    let log = vec![format!(
        "{} Floquet-periodic pair(s): {n} dofs, {nred} after constraints; complex Hermitian eigenproblem, Krylov dimension {}{}",
        asm.n_pairs,
        solved.dim,
        if solved.converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    Ok(SolveResult { order: order as u8, dim: d as u8, solutions, scalars, log, dofs: n, ..Default::default() })
}
