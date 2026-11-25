// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Shared kernel of the 2D elastic cross-section models (`axisymmetric_eigenfrequency` and `modes`):
//! stiffness, mass, low-reflecting damping, curvilinear PML stretching (Cartesian, Cylindrical, and
//! Spherical touching `r = 0`), pre-solve Christoffel backward-wave stabilization, outer PML clamping,
//! signed quality factors `Q = Re f / (2 Im f)`, radiated power `P_pml`, and native `Tri` + `Quad`
//! mixed-cell assembly.

use crate::floquet::mixed_nodal_relations;
use crate::lrb::AMPLIFYING_IM_TOL;
use crate::pml::{outer_clamped_dofs_mixed, outer_clamped_dofs_simplex, ResolvedPmlDomains};
use crate::problem::*;
use femcore::fem::check_boundary_labels;
use femcore::fem::*;
use femcore::fem_mixed::{assemble_mass_mixed, LagrangeSpace, MixedLagrangeDofMap};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCsr, C};
use femcore::material::M6;
use femcore::mesh::mixed::MixedDiscretisation;
use nyx_physics::{
    cell_fields_of, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
use std::collections::BTreeSet;
use std::time::Instant;

fn check_labels(mesh: &Discretisation, bcs: &[Bc]) -> Result<(), String> {
    check_boundary_labels(mesh, bcs.iter().flat_map(|b| b.labels.iter().copied()))
}

fn to_complex(a: &Csr) -> CCsr {
    a.to_complex()
}

/// The wafer axis (0 = x, 1 = y, 2 = z) of each displacement block: x, z, y (for a ring r, z, theta).
const PLANE_AXES: [usize; 3] = [0, 2, 1];

/// The complex eigenvectors behind a cross-section solid result: displacement over all `3 * n_dofs` dofs (component blocks), the complex
/// frequencies in Hz, in the order of the solutions.
pub struct CrossModes {
    pub m_az: f64,
    pub freqs: Vec<C>,
    pub vectors: Vec<Vec<C>>,
    pub c6: Vec<M6>,
    pub rho: Vec<f64>,
    pub converged: bool,
}

#[derive(Clone, Copy, Debug)]
pub enum CrossSection {
    Axisymmetric(f64),
    Waveguide(f64),
}

pub fn solve_solid_axisym(mesh: &Discretisation, spec: &Problem, m: f64) -> Result<(SolveResult, CrossModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Axisymmetric(m))
}

pub fn solve_solid_waveguide(
    mesh: &Discretisation,
    spec: &Problem,
    q: f64,
) -> Result<(SolveResult, CrossModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Waveguide(q))
}

pub fn solve_axisymmetric(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let Study::AxisymmetricEigenfrequency { azimuthal_order, .. } = problem.study else {
        return Err("solve_axisymmetric needs the axisymmetric_eigenfrequency study".into());
    };
    problem.check_axisymmetric().map_err(|e| e.to_string())?;
    check_half_plane(mesh)?;
    let (mut result, _) = solve_cross_section(mesh, problem, CrossSection::Axisymmetric(azimuthal_order as f64))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("axisymmetric study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

pub fn solve_axisymmetric_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    problem: &Problem,
) -> Result<SolveResult, String> {
    let Study::AxisymmetricEigenfrequency { azimuthal_order, .. } = problem.study else {
        return Err("solve_axisymmetric needs the axisymmetric_eigenfrequency study".into());
    };
    problem.check_axisymmetric().map_err(|e| e.to_string())?;
    check_half_plane_mixed(mixed)?;
    let (mut result, _) =
        solve_cross_section_mixed(raw_mesh, mixed, problem, CrossSection::Axisymmetric(azimuthal_order as f64))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("axisymmetric study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

pub fn solve_waveguide_modes(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let Study::Modes { wavevector, .. } = problem.study else {
        return Err("solve_waveguide_modes needs the modes study".into());
    };
    problem.check_modes().map_err(|e| e.to_string())?;
    let (mut result, _) = solve_cross_section(mesh, problem, CrossSection::Waveguide(wavevector))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("waveguide modes study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

pub fn solve_waveguide_modes_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    problem: &Problem,
) -> Result<SolveResult, String> {
    let Study::Modes { wavevector, .. } = problem.study else {
        return Err("solve_waveguide_modes needs the modes study".into());
    };
    problem.check_modes().map_err(|e| e.to_string())?;
    let (mut result, _) = solve_cross_section_mixed(raw_mesh, mixed, problem, CrossSection::Waveguide(wavevector))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("waveguide modes study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

fn check_half_plane(mesh: &Discretisation) -> Result<(), String> {
    let xs = (0..mesh.n_points()).map(|i| mesh.point(i)[0]);
    let (lo, hi) = xs.fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| (l.min(x), h.max(x)));
    if lo < -1e-9 * (hi - lo).abs() {
        return Err(format!("the (r, z) half plane must lie in x >= 0 (x = r): the mesh reaches x = {lo:.6e}"));
    }
    Ok(())
}

fn check_half_plane_mixed(mixed: &MixedDiscretisation) -> Result<(), String> {
    let xs = (0..mixed.n_points()).map(|i| mixed.point(i)[0]);
    let (lo, hi) = xs.fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| (l.min(x), h.max(x)));
    if lo < -1e-9 * (hi - lo).abs() {
        return Err(format!("the (r, z) half plane must lie in x >= 0 (x = r): the mesh reaches x = {lo:.6e}"));
    }
    Ok(())
}

pub fn solve_cross_section(
    mesh: &Discretisation,
    spec: &Problem,
    kind: CrossSection,
) -> Result<(SolveResult, CrossModes), String> {
    let t0 = Instant::now();
    let asm = assemble_cross_section(mesh, spec, kind)?;
    let solved = solve_cross_section_modes(&asm)?;
    let result = cross_section_fields(mesh, &asm, &solved)?;
    let result = SolveResult { seconds: t0.elapsed().as_secs_f64(), ..result };
    Ok((result, solved.modes))
}

pub fn solve_cross_section_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
    kind: CrossSection,
) -> Result<(SolveResult, CrossModes), String> {
    if mixed.dim != 2 {
        return Err("a cross-section model needs a 2D mesh".into());
    }
    let t0 = Instant::now();
    let (Study::Eigenfrequency { n_modes, f_search }
    | Study::AxisymmetricEigenfrequency { n_modes, f_search, .. }
    | Study::Modes { n_modes, f_search, .. }) = spec.study.clone()
    else {
        return Err("the cross-section solid models need an eigenfrequency study".into());
    };
    let order = spec.order;
    crate::solid::mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let n1 = dof_map.n_dofs;
    let n = 3 * n1;

    let mut c6 = crate::solid::mixed_cell_stiffness_skipping(&spec.materials, mixed, &spec.inactive)?;
    c6 = c6.iter().map(|c| femcore::section::plane_stiffness(c, 2)).collect();
    if matches!(kind, CrossSection::Axisymmetric(_)) {
        for m in c6.iter_mut() {
            *m = crate::axisym::average_about_axis(m, 1, 16);
        }
    }
    let rho = crate::solid::mixed_cell_rho(&spec.materials, mixed, &spec.inactive)?;

    let is_axisym = matches!(kind, CrossSection::Axisymmetric(_));
    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_mixed(mixed, spec, is_axisym, f_search)?)
    } else {
        None
    };

    let mut fixed_set: BTreeSet<usize> =
        crate::solid::mixed_untouched_dofs(mixed, &dof_map, 3, &spec.inactive).into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = dof_map.boundary_dofs(mixed, &bc.labels);
            for k in 0..3 {
                if fixed[PLANE_AXES[k]] {
                    for &dof in &bd {
                        fixed_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    let axis_tol = if is_axisym { crate::axisym::axis_tolerance_mixed(mixed) } else { 0.0 };
    let axis_dofs = match kind {
        CrossSection::Axisymmetric(m_az) => {
            let (axis_pairs, axis_dofs) = crate::axisym::axis_constraints_mixed(mixed, &dof_map, m_az)?;
            fixed_set.extend(axis_pairs.into_iter().map(|(k, dof)| dof + k * n1));
            axis_dofs
        }
        CrossSection::Waveguide(_) => 0,
    };
    if let Some(res) = &resolved {
        for dof in outer_clamped_dofs_mixed(mixed, &dof_map, spec, res, is_axisym, axis_tol) {
            for k in 0..3 {
                fixed_set.insert(dof + k * n1);
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();

    let (k, mc) = match kind {
        CrossSection::Axisymmetric(m) => {
            crate::axisym::assemble_axisym_mixed(mixed, &dof_map, &c6, &rho, m, resolved.as_ref(), &spec.inactive)?
        }
        CrossSection::Waveguide(q) => crate::waveguide::assemble_waveguide_mixed(
            mixed,
            &dof_map,
            &c6,
            &rho,
            q,
            resolved.as_ref(),
            &spec.inactive,
        )?,
    };

    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    if !lrb_labels.is_empty() && f_search <= 0.0 {
        return Err(
            "low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into()
        );
    }

    let (free, kf, mf, cf, red) = if matches!(kind, CrossSection::Waveguide(_)) && !spec.periodic.is_empty() {
        if !lrb_labels.is_empty() {
            return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
        }
        let mut relations = Vec::new();
        for p in &spec.periodic {
            relations.extend(mixed_nodal_relations(mixed, &dof_map, 3, p)?);
        }
        let red = femcore::periodic::Reduction::build(n, &relations, &fixed);
        if n_modes < 1 || n_modes + 1 >= red.n_reduced() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let pm = red.prolong();
        let kf = femcore::periodic::congruence(&pm, &k);
        let mf = femcore::periodic::congruence(&pm, &mc);
        (vec![], kf, mf, None, Some(red))
    } else {
        let free = free_dofs(n, &fixed);
        if n_modes < 1 || n_modes + 1 >= free.len() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let kf = k.submatrix(&free, &free);
        let mf = mc.submatrix(&free, &free);
        let cf = if lrb_labels.is_empty() {
            None
        } else {
            Some(
                crate::lrb::damping_matrix_mixed(mixed, &dof_map, &c6, &rho, &lrb_labels, 3, 1.0)
                    .submatrix(&free, &free),
            )
        };
        (free, kf, mf, cf, None)
    };

    let pml_cells: Vec<bool> =
        (0..mixed.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mixed.cell_labels[ci]))).collect();
    let has_loss = pml_cells.iter().any(|&p| p) || !lrb_labels.is_empty();
    let diagnostics = resolved.map_or_else(Vec::new, |r| r.diagnostics);

    let asm = CrossAssembly {
        kind,
        order,
        n,
        n1,
        kf,
        mf,
        cf,
        free,
        fixed,
        c6,
        rho,
        n_modes,
        f_search,
        has_loss,
        axis_dofs,
        pml_cells,
        red,
        diagnostics,
    };
    let solved = solve_cross_section_modes(&asm)?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let result = cross_section_fields_mixed(mixed, &dof_map, &proj, &asm, &solved)?;
    let result = SolveResult { seconds: t0.elapsed().as_secs_f64(), ..result };
    Ok((result, solved.modes))
}

/// Stage 1 of the cross-section solid kernels: stiffness, mass and (low-reflecting) damping on the free dofs, the
/// constraints, and the material data the coupling integrals need.
pub struct CrossAssembly {
    pub kind: CrossSection,
    pub order: usize,
    pub n: usize,
    pub n1: usize,
    pub kf: CCsr,
    pub mf: CCsr,
    pub cf: Option<Csr>,
    pub free: Vec<usize>,
    pub fixed: Vec<usize>,
    pub c6: Vec<M6>,
    pub rho: Vec<f64>,
    pub n_modes: usize,
    pub f_search: f64,
    pub has_loss: bool,
    pub axis_dofs: usize,
    pub pml_cells: Vec<bool>,
    pub red: Option<femcore::periodic::Reduction>,
    pub(crate) diagnostics: Vec<String>,
}

pub fn assemble_cross_section(
    mesh: &Discretisation,
    spec: &Problem,
    kind: CrossSection,
) -> Result<CrossAssembly, String> {
    if mesh.dim != 2 {
        return Err("a cross-section model needs a 2D mesh".into());
    }
    let (Study::Eigenfrequency { n_modes, f_search }
    | Study::AxisymmetricEigenfrequency { n_modes, f_search, .. }
    | Study::Modes { n_modes, f_search, .. }) = spec.study.clone()
    else {
        return Err("the cross-section solid models need an eigenfrequency study".into());
    };
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let mut c6 = spec.materials.cell_stiffness_skipping(mesh, &spec.inactive)?;
    c6 = c6.iter().map(|c| femcore::section::plane_stiffness(c, 2)).collect();
    let is_axisym = matches!(kind, CrossSection::Axisymmetric(_));
    if is_axisym {
        for m in c6.iter_mut() {
            *m = crate::axisym::average_about_axis(m, 1, 16);
        }
    }
    let rho = crate::solid::cell_rho(&spec.materials, mesh, &spec.inactive)?;

    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_simplex(mesh, spec, is_axisym, f_search)?)
    } else {
        None
    };

    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut fixed_set: BTreeSet<usize> = untouched_dofs(mesh, order, 3, &spec.inactive).into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = boundary_dofs(mesh, &bc.labels, order);
            for k in 0..3 {
                if fixed[PLANE_AXES[k]] {
                    for &dof in &bd {
                        fixed_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    let axis_tol = if is_axisym { crate::axisym::axis_tolerance(mesh) } else { 0.0 };
    let axis_dofs = match kind {
        CrossSection::Axisymmetric(m_az) => {
            let (axis_pairs, axis_dofs) = crate::axisym::axis_constraints(mesh, order, m_az)?;
            fixed_set.extend(axis_pairs.into_iter().map(|(k, dof)| dof + k * n1));
            axis_dofs
        }
        CrossSection::Waveguide(_) => 0,
    };
    if let Some(res) = &resolved {
        for dof in outer_clamped_dofs_simplex(mesh, order, spec, res, is_axisym, axis_tol) {
            for k in 0..3 {
                fixed_set.insert(dof + k * n1);
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();

    let (k, mc) = match (kind, &resolved) {
        (CrossSection::Axisymmetric(m), Some(res)) => {
            crate::axisym::assemble_axisym_pml_simplex(mesh, &c6, &rho, order, m, res, &spec.inactive)?
        }
        (CrossSection::Axisymmetric(m), None) => (
            crate::axisym::assemble_stiffness_axisym(mesh, &c6, order, m, None)?,
            crate::axisym::assemble_mass_axisym(mesh, &rho, order, None),
        ),
        (CrossSection::Waveguide(q), Some(res)) => {
            crate::waveguide::assemble_waveguide_pml_simplex(mesh, &c6, &rho, order, q, res, &spec.inactive)?
        }
        (CrossSection::Waveguide(q), None) => {
            let k0 = crate::waveguide::assemble_stiffness_waveguide(mesh, &c6, order, q, None)?;
            let m1 = assemble_mass(mesh, order, &rho);
            let m0 = to_complex(&Csr::block_diag(&[&m1, &m1, &m1]));
            (k0, m0)
        }
    };

    let lrb_labels: Vec<i32> = spec
        .bcs
        .iter()
        .filter(|b| matches!(b.kind, BcKind::LowReflecting))
        .flat_map(|b| b.labels.iter().copied())
        .collect();
    if !lrb_labels.is_empty() && f_search <= 0.0 {
        return Err(
            "low-reflecting boundaries need a search frequency (the problem is quadratic in the frequency)".into()
        );
    }
    let (free, kf, mf, cf, red) = if matches!(kind, CrossSection::Waveguide(_)) && !spec.periodic.is_empty() {
        if !lrb_labels.is_empty() {
            return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
        }
        let mut relations = vec![];
        for p in &spec.periodic {
            let mt = femcore::periodic::match_vertices(mesh, p)?;
            relations.extend(femcore::periodic::nodal_relations(mesh, order, 3, &mt, p.phase)?);
        }
        let red = femcore::periodic::Reduction::build(n, &relations, &fixed);
        if n_modes < 1 || n_modes + 1 >= red.n_reduced() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let pm = red.prolong();
        let kf = femcore::periodic::congruence(&pm, &k);
        let mf = femcore::periodic::congruence(&pm, &mc);
        (vec![], kf, mf, None, Some(red))
    } else {
        let free = free_dofs(n, &fixed);
        if n_modes < 1 || n_modes + 1 >= free.len() {
            return Err("number of modes must be between 1 and the number of free dofs".into());
        }
        let kf = k.submatrix(&free, &free);
        let mf = mc.submatrix(&free, &free);
        let cf = if lrb_labels.is_empty() {
            None
        } else {
            Some(crate::lrb::damping_matrix(mesh, &c6, &rho, &lrb_labels, order, 3).submatrix(&free, &free))
        };
        (free, kf, mf, cf, None)
    };
    let pml_cells: Vec<bool> =
        (0..mesh.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mesh.cell_labels[ci]))).collect();
    let has_loss = pml_cells.iter().any(|&p| p) || !lrb_labels.is_empty();
    let diagnostics = resolved.map_or_else(Vec::new, |r| r.diagnostics);

    Ok(CrossAssembly {
        kind,
        order,
        n,
        n1,
        kf,
        mf,
        cf,
        free,
        fixed,
        c6,
        rho,
        n_modes,
        f_search,
        has_loss,
        axis_dofs,
        pml_cells,
        red,
        diagnostics,
    })
}

pub struct CrossSolve {
    pub modes: CrossModes,
    pub dim: usize,
}

pub fn solve_cross_section_modes(asm: &CrossAssembly) -> Result<CrossSolve, String> {
    let (CrossSection::Axisymmetric(m_az) | CrossSection::Waveguide(m_az)) = asm.kind;
    let (n_modes, f_search, free, n) = (asm.n_modes, asm.f_search, &asm.free, asm.n);
    let (kf, mf) = (&asm.kf, &asm.mf);
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    type Found = (Vec<(C, usize)>, Vec<Vec<C>>, bool, usize);
    let (mut modes, vectors, converged, kdim): Found = if let Some(cf) = asm.cf.as_ref() {
        let (omegas, vecs, conv, dim) =
            crate::lrb::quadratic_modes(kf, mf, cf, 2.0 * std::f64::consts::PI * f_search, n_modes)?;
        let modes = omegas.iter().enumerate().map(|(i, w)| (*w / (2.0 * std::f64::consts::PI), i)).collect();
        (modes, vecs, conv, dim)
    } else {
        let pairs = femcore::linalg::driver::eigs_complex(kf, mf, c(sigma, 0.0), n_modes, 1e-9, 400)?;
        let modes = pairs
            .values
            .iter()
            .enumerate()
            .filter(|(i, _)| pairs.residuals.get(*i).is_none_or(|r| *r < 1e-6))
            .map(|(i, l)| {
                let mut f = l.sqrt() / (2.0 * std::f64::consts::PI);
                if f.re < 0.0 {
                    f = -f;
                }
                (f, i)
            })
            .collect();
        (modes, pairs.vectors.clone(), pairs.converged, pairs.dim)
    };
    modes.sort_by(|a, b| a.0.re.partial_cmp(&b.0.re).unwrap());
    let mut out_freqs: Vec<C> = vec![];
    let mut out_vectors: Vec<Vec<C>> = vec![];
    for &(fz, i) in modes.iter() {
        let u = if let Some(red) = &asm.red {
            red.expand(&vectors[i])
        } else {
            let mut u = vec![C::new(0.0, 0.0); n];
            for (a, &dof) in free.iter().enumerate() {
                u[dof] = vectors[i][a];
            }
            u
        };
        out_freqs.push(fz);
        out_vectors.push(u);
    }
    Ok(CrossSolve {
        modes: CrossModes {
            m_az,
            freqs: out_freqs,
            vectors: out_vectors,
            c6: asm.c6.clone(),
            rho: asm.rho.clone(),
            converged,
        },
        dim: kdim,
    })
}

fn mass_norm(m: &CCsr, u: &[C]) -> f64 {
    let mut acc = 0.0;
    for r in 0..m.nrows() {
        let mut row = c(0.0, 0.0);
        for k in m.indptr()[r]..m.indptr()[r + 1] {
            row += m.values()[k] * u[m.indices()[k] as usize];
        }
        acc += (u[r].conj() * row).re;
    }
    acc
}

pub fn cross_section_fields(
    mesh: &Discretisation,
    asm: &CrossAssembly,
    solved: &CrossSolve,
) -> Result<SolveResult, String> {
    let (n1, order, has_loss) = (asm.n1, asm.order, asm.has_loss);
    let modes = &solved.modes;
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut unstable_count = 0usize;

    let any_pml = asm.pml_cells.iter().any(|&p| p);
    let has_interior = asm.pml_cells.iter().zip(&asm.rho).any(|(&p, &r)| !p && r > 0.0);
    let (m_all, m_pml) = if any_pml {
        let rho_pml: Vec<f64> = asm.rho.iter().zip(&asm.pml_cells).map(|(r, &p)| if p { *r } else { 0.0 }).collect();
        match asm.kind {
            CrossSection::Axisymmetric(_) => (
                Some(crate::axisym::assemble_mass_axisym(mesh, &asm.rho, order, None)),
                Some(crate::axisym::assemble_mass_axisym(mesh, &rho_pml, order, None)),
            ),
            CrossSection::Waveguide(_) => {
                let m1_all = assemble_mass(mesh, order, &asm.rho);
                let m1_pml = assemble_mass(mesh, order, &rho_pml);
                (
                    Some(to_complex(&Csr::block_diag(&[&m1_all, &m1_all, &m1_all]))),
                    Some(to_complex(&Csr::block_diag(&[&m1_pml, &m1_pml, &m1_pml]))),
                )
            }
        }
    } else {
        (None, None)
    };

    let mut domains: Vec<i32> = (0..mesh.n_cells())
        .filter(|&ci| !asm.pml_cells[ci] && asm.rho[ci] != 0.0)
        .map(|ci| mesh.cell_labels[ci])
        .collect();
    domains.sort_unstable();
    domains.dedup();
    let class: Vec<usize> = (0..mesh.n_cells())
        .map(|ci| {
            if asm.pml_cells[ci] {
                usize::MAX
            } else {
                domains.binary_search(&mesh.cell_labels[ci]).unwrap_or(usize::MAX)
            }
        })
        .collect();
    let energy = match asm.kind {
        CrossSection::Axisymmetric(_) => {
            crate::axisym::mode_energy_axisym(mesh, &asm.rho, order, &class, domains.len(), &modes.vectors)
        }
        CrossSection::Waveguide(_) => {
            crate::waveguide::mode_energy_waveguide(mesh, &asm.rho, order, &class, domains.len(), &modes.vectors)
        }
    };
    let share_names = match asm.kind {
        CrossSection::Axisymmetric(_) => ["radial", "axial", "azimuthal"],
        CrossSection::Waveguide(_) => ["lateral", "vertical", "longitudinal"],
    };
    let vol_factor = match asm.kind {
        CrossSection::Axisymmetric(_) => 2.0 * std::f64::consts::PI,
        CrossSection::Waveguide(_) => 1.0,
    };

    for (j, &fz) in modes.freqs.iter().enumerate() {
        let freq = fz.re;
        let is_amplifying = fz.im < -AMPLIFYING_IM_TOL * fz.re.abs().max(1e-30);
        if is_amplifying {
            unstable_count += 1;
        }

        let u = &modes.vectors[j];
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let comp = |k: usize| -> Vec<C> { (0..n1).map(|a| u[k * n1 + a] * rot * scale).collect() };
        let (ux, uy, uz) = (comp(0), comp(1), comp(2));
        let mag: Vec<f64> = (0..n1).map(|a| (ux[a].norm_sqr() + uy[a].norm_sqr() + uz[a].norm_sqr()).sqrt()).collect();
        let re = |v: &[C]| -> Vec<f64> { v.iter().map(|z| z.re).collect() };
        let names = ["ux", "uz", "uy"];
        let fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar(format!("Re {}", names[0]), "a.u.", re(&ux)),
            Field::scalar(format!("Re {}", names[1]), "a.u.", re(&uy)),
            Field::scalar(format!("Re {}", names[2]), "a.u.", re(&uz)),
            Field::scalar(format!("|{}|", names[2]), "a.u.", uz.iter().map(|z| z.norm()).collect()),
        ];
        let (layer_fraction, phys_mass_norm) = match (&m_all, &m_pml) {
            (Some(ma), Some(mp)) => {
                let (tot, in_layer) = (mass_norm(ma, u), mass_norm(mp, u));
                let lf = if tot > 0.0 { (in_layer / tot).clamp(0.0, 1.0) } else { 0.0 };
                (lf, (tot - in_layer).max(0.0) * scale * scale * vol_factor)
            }
            _ => (0.0, 0.0),
        };
        let nloc = mesh.nloc(order);
        let mut cells = femcore::fields::CellFields::new(mesh, order);
        for (k, name) in names.iter().enumerate() {
            let (mut re_c, mut im_c) =
                (Vec::with_capacity(mesh.n_cells() * nloc), Vec::with_capacity(mesh.n_cells() * nloc));
            let mut dofs = Vec::new();
            for cell in 0..mesh.n_cells() {
                mesh.cell_dofs(cell, order, &mut dofs);
                for &d in dofs.iter().take(nloc) {
                    let z = u[k * n1 + d] * rot * scale;
                    re_c.push(z.re);
                    im_c.push(z.im);
                }
            }
            cells.push_complex(*name, "a.u.", re_c, im_c);
        }
        if has_loss && fz.im.abs() > 0.0 {
            let q = fz.re / (2.0 * fz.im);
            let label = if is_amplifying {
                format!(
                    "mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]",
                    j + 1,
                    fz.re,
                    fz.im,
                    q
                )
            } else {
                format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, fz.re, fz.im, q)
            };
            solutions.push(Solution {
                label,
                fields,
                cell_fields: cell_fields_of(&cells),
                value: Some(freq),
                eigenvalue: Some([fz.re, fz.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", fz.re, fz.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        } else {
            solutions.push(Solution {
                label: format!("mode {}: {:.6} Hz", j + 1, freq),
                fields,
                cell_fields: cell_fields_of(&cells),
                value: Some(freq),
                eigenvalue: Some([fz.re, fz.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
        }
        scalars.push((format!("layer fraction {}", j + 1), Scalar::Num(layer_fraction)));
        scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(fz.re)));
        scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(fz.im)));
        if has_interior && any_pml {
            let two_pi = 2.0 * std::f64::consts::PI;
            let p_pml = (two_pi * fz.re) * (two_pi * fz.im) * phys_mass_norm;
            scalars.push((format!("P_pml {} (W)", j + 1), Scalar::Num(p_pml)));
        }
        let total: f64 = energy[j].iter().map(|e| e[0] + e[1] + e[2]).sum();
        let share = |x: f64| if total > 0.0 { x / total } else { 0.0 };
        for (k, name) in share_names.iter().enumerate() {
            let part: f64 = energy[j].iter().map(|e| e[k]).sum();
            scalars.push((format!("{name} share {}", j + 1), Scalar::Num(share(part))));
        }
        for (d, label) in domains.iter().enumerate() {
            let e = energy[j][d];
            scalars.push((format!("share in domain {label} {}", j + 1), Scalar::Num(share(e[0] + e[1] + e[2]))));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let log = finish_cross_section_log(asm, solved, unstable_count, any_pml, has_interior, &mut scalars);
    Ok(SolveResult { order: order as u8, dim: 2, solutions, scalars, log, dofs: asm.n, ..Default::default() })
}

fn cross_section_fields_mixed(
    mixed: &MixedDiscretisation,
    dof_map: &MixedLagrangeDofMap,
    proj: &MixedToSimplexProjector,
    asm: &CrossAssembly,
    solved: &CrossSolve,
) -> Result<SolveResult, String> {
    let (n1, order, has_loss) = (asm.n1, asm.order, asm.has_loss);
    let modes = &solved.modes;
    let mut solutions = Vec::new();
    let mut scalars = Vec::new();
    let mut unstable_count = 0usize;

    let any_pml = asm.pml_cells.iter().any(|&p| p);
    let has_interior = asm.pml_cells.iter().zip(&asm.rho).any(|(&p, &r)| !p && r > 0.0);
    let (m_all, m_pml) = if any_pml {
        let rho_pml: Vec<f64> = asm.rho.iter().zip(&asm.pml_cells).map(|(r, &p)| if p { *r } else { 0.0 }).collect();
        match asm.kind {
            CrossSection::Axisymmetric(_) => (
                Some(crate::axisym::assemble_mass_axisym_mixed(mixed, dof_map, &asm.rho)),
                Some(crate::axisym::assemble_mass_axisym_mixed(mixed, dof_map, &rho_pml)),
            ),
            CrossSection::Waveguide(_) => {
                let m1_all = assemble_mass_mixed(mixed, dof_map, &asm.rho);
                let m1_pml = assemble_mass_mixed(mixed, dof_map, &rho_pml);
                (
                    Some(to_complex(&Csr::block_diag(&[&m1_all, &m1_all, &m1_all]))),
                    Some(to_complex(&Csr::block_diag(&[&m1_pml, &m1_pml, &m1_pml]))),
                )
            }
        }
    } else {
        (None, None)
    };

    let mut domains: Vec<i32> = (0..mixed.n_cells())
        .filter(|&ci| !asm.pml_cells[ci] && asm.rho[ci] != 0.0)
        .map(|ci| mixed.cell_labels[ci])
        .collect();
    domains.sort_unstable();
    domains.dedup();
    let class: Vec<usize> = (0..mixed.n_cells())
        .map(|ci| {
            if asm.pml_cells[ci] {
                usize::MAX
            } else {
                domains.binary_search(&mixed.cell_labels[ci]).unwrap_or(usize::MAX)
            }
        })
        .collect();
    let energy = match asm.kind {
        CrossSection::Axisymmetric(_) => {
            crate::axisym::mode_energy_axisym_mixed(mixed, dof_map, &asm.rho, &class, domains.len(), &modes.vectors)
        }
        CrossSection::Waveguide(_) => crate::waveguide::mode_energy_waveguide_mixed(
            mixed,
            dof_map,
            &asm.rho,
            &class,
            domains.len(),
            &modes.vectors,
        ),
    };
    let share_names = match asm.kind {
        CrossSection::Axisymmetric(_) => ["radial", "axial", "azimuthal"],
        CrossSection::Waveguide(_) => ["lateral", "vertical", "longitudinal"],
    };
    let vol_factor = match asm.kind {
        CrossSection::Axisymmetric(_) => 2.0 * std::f64::consts::PI,
        CrossSection::Waveguide(_) => 1.0,
    };

    let mut m_offsets = Vec::with_capacity(mixed.n_cells() + 1);
    m_offsets.push(0);
    for mc in 0..mixed.n_cells() {
        m_offsets.push(m_offsets[mc] + dof_map.cell_dofs(mc).len());
    }

    for (j, &fz) in modes.freqs.iter().enumerate() {
        let freq = fz.re;
        let is_amplifying = fz.im < -AMPLIFYING_IM_TOL * fz.re.abs().max(1e-30);
        if is_amplifying {
            unstable_count += 1;
        }

        let u = &modes.vectors[j];
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);

        let mut re_blocks = Vec::with_capacity(3 * n1);
        let mut im_blocks = Vec::with_capacity(3 * n1);
        for k in 0..3 {
            for a in 0..n1 {
                let z = u[k * n1 + a] * rot * scale;
                re_blocks.push(z.re);
                im_blocks.push(z.im);
            }
        }
        let re_proj = proj.project_lagrange_blocks(mixed, dof_map, &re_blocks, 3, order);
        let im_proj = proj.project_lagrange_blocks(mixed, dof_map, &im_blocks, 3, order);
        let n_out = re_proj.len() / 3;
        let mag: Vec<f64> = (0..n_out)
            .map(|a| {
                (0..3)
                    .map(|k| {
                        let r0 = re_proj[k * n_out + a];
                        let i0 = im_proj[k * n_out + a];
                        r0 * r0 + i0 * i0
                    })
                    .sum::<f64>()
                    .sqrt()
            })
            .collect();
        let uy_abs: Vec<f64> = (0..n_out).map(|a| re_proj[2 * n_out + a].hypot(im_proj[2 * n_out + a])).collect();

        let names = ["ux", "uz", "uy"];
        let fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar("Re ux", "a.u.", re_proj[0..n_out].to_vec()),
            Field::scalar("Re uz", "a.u.", re_proj[n_out..2 * n_out].to_vec()),
            Field::scalar("Re uy", "a.u.", re_proj[2 * n_out..3 * n_out].to_vec()),
            Field::scalar("|uy|", "a.u.", uy_abs),
        ];

        let mut cell_fields = Vec::with_capacity(3);
        for (k, name) in names.iter().enumerate() {
            let mut re_c = Vec::with_capacity(*m_offsets.last().unwrap());
            let mut im_c = Vec::with_capacity(*m_offsets.last().unwrap());
            for mc in 0..mixed.n_cells() {
                for &gd in dof_map.cell_dofs(mc) {
                    re_c.push(re_blocks[k * n1 + gd]);
                    im_c.push(im_blocks[k * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                dof_map.space,
                *name,
                "a.u.",
                &re_c,
                &im_c,
                &m_offsets,
                order,
            ));
        }

        let (layer_fraction, phys_mass_norm) = match (&m_all, &m_pml) {
            (Some(ma), Some(mp)) => {
                let (tot, in_layer) = (mass_norm(ma, u), mass_norm(mp, u));
                let lf = if tot > 0.0 { (in_layer / tot).clamp(0.0, 1.0) } else { 0.0 };
                (lf, (tot - in_layer).max(0.0) * scale * scale * vol_factor)
            }
            _ => (0.0, 0.0),
        };

        if has_loss && fz.im.abs() > 0.0 {
            let q = fz.re / (2.0 * fz.im);
            let label = if is_amplifying {
                format!(
                    "mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e}) [AMPLIFYING backward-wave mode]",
                    j + 1,
                    fz.re,
                    fz.im,
                    q
                )
            } else {
                format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, fz.re, fz.im, q)
            };
            solutions.push(Solution {
                label,
                fields,
                cell_fields,
                value: Some(freq),
                eigenvalue: Some([fz.re, fz.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", fz.re, fz.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        } else {
            solutions.push(Solution {
                label: format!("mode {}: {:.6} Hz", j + 1, freq),
                fields,
                cell_fields,
                value: Some(freq),
                eigenvalue: Some([fz.re, fz.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
        }
        scalars.push((format!("layer fraction {}", j + 1), Scalar::Num(layer_fraction)));
        scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(fz.re)));
        scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(fz.im)));
        if has_interior && any_pml {
            let two_pi = 2.0 * std::f64::consts::PI;
            let p_pml = (two_pi * fz.re) * (two_pi * fz.im) * phys_mass_norm;
            scalars.push((format!("P_pml {} (W)", j + 1), Scalar::Num(p_pml)));
        }
        let total: f64 = energy[j].iter().map(|e| e[0] + e[1] + e[2]).sum();
        let share = |x: f64| if total > 0.0 { x / total } else { 0.0 };
        for (k, name) in share_names.iter().enumerate() {
            let part: f64 = energy[j].iter().map(|e| e[k]).sum();
            scalars.push((format!("{name} share {}", j + 1), Scalar::Num(share(part))));
        }
        for (d, label) in domains.iter().enumerate() {
            let e = energy[j][d];
            scalars.push((format!("share in domain {label} {}", j + 1), Scalar::Num(share(e[0] + e[1] + e[2]))));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let log = finish_cross_section_log(asm, solved, unstable_count, any_pml, has_interior, &mut scalars);
    Ok(SolveResult { order: order as u8, dim: 2, solutions, scalars, log, dofs: asm.n, ..Default::default() })
}

fn finish_cross_section_log(
    asm: &CrossAssembly,
    solved: &CrossSolve,
    unstable_count: usize,
    any_pml: bool,
    has_interior: bool,
    scalars: &mut Vec<(String, Scalar)>,
) -> Vec<String> {
    let m_az = solved.modes.m_az;
    let what = match asm.kind {
        CrossSection::Axisymmetric(_) => format!("axisymmetric, azimuthal order m = {m_az}"),
        CrossSection::Waveguide(q) => format!("waveguide modes, longitudinal wavevector q = {q} rad/m"),
    };
    let mut log = asm.diagnostics.clone();
    log.push(format!(
        "{what}; complex {} eigenproblem, Krylov dimension {}{}",
        if asm.has_loss { "symmetric" } else { "Hermitian" },
        solved.dim,
        if solved.modes.converged { "" } else { " (not every pair converged; those are omitted)" }
    ));
    if asm.axis_dofs > 0 {
        let held = match m_az.round().abs() as u64 {
            0 => "u_r and u_theta",
            1 => "u_z",
            _ => "u_r, u_z and u_theta",
        };
        log.push(format!("axis r = 0: {} dofs, {held} held at zero (regularity for m = {m_az})", asm.axis_dofs));
    }
    if asm.fixed.is_empty() {
        log.push("no fixed boundary: the lowest modes may be rigid-body-like".into());
    }
    if unstable_count > 0 || (has_interior && any_pml) {
        scalars.push(("unstable_modes".into(), Scalar::Num(unstable_count as f64)));
    }
    if unstable_count > 0 {
        log.push(format!(
            "WARNING: {unstable_count} amplifying backward-wave mode(s) detected (Im f < -1e-9 Re f); enable stabilization = \"auto\", \"mpml\", or \"adiabatic\""
        ));
    }
    log
}
