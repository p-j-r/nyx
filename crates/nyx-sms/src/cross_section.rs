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
use std::collections::{BTreeSet, VecDeque};
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
    let mut result = if problem.uses_piezoelectric() {
        solve_axisymmetric_piezo(mesh, problem, azimuthal_order as f64)?
    } else {
        solve_cross_section(mesh, problem, CrossSection::Axisymmetric(azimuthal_order as f64))?.0
    };
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
    let mut result = if problem.uses_piezoelectric() {
        solve_axisymmetric_piezo_mixed(raw_mesh, mixed, problem, azimuthal_order as f64)?
    } else {
        solve_cross_section_mixed(raw_mesh, mixed, problem, CrossSection::Axisymmetric(azimuthal_order as f64))?.0
    };
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
    let mut result = if problem.uses_piezoelectric() {
        solve_waveguide_piezo(mesh, problem, wavevector)?
    } else {
        solve_cross_section(mesh, problem, CrossSection::Waveguide(wavevector))?.0
    };
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
    let mut result = if problem.uses_piezoelectric() {
        solve_waveguide_piezo_mixed(raw_mesh, mixed, problem, wavevector)?
    } else {
        solve_cross_section_mixed(raw_mesh, mixed, problem, CrossSection::Waveguide(wavevector))?.0
    };
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
        is_piezo: false,
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
    pub is_piezo: bool,
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
        is_piezo: false,
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
        let u_mech = &u[..3 * n1];
        let big = u_mech.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let comp = |k: usize| -> Vec<C> { (0..n1).map(|a| u[k * n1 + a] * rot * scale).collect() };
        let (ux, uy, uz) = (comp(0), comp(1), comp(2));
        let mag: Vec<f64> = (0..n1).map(|a| (ux[a].norm_sqr() + uy[a].norm_sqr() + uz[a].norm_sqr()).sqrt()).collect();
        let re = |v: &[C]| -> Vec<f64> { v.iter().map(|z| z.re).collect() };
        let names = ["ux", "uz", "uy"];
        let mut fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar(format!("Re {}", names[0]), "a.u.", re(&ux)),
            Field::scalar(format!("Re {}", names[1]), "a.u.", re(&uy)),
            Field::scalar(format!("Re {}", names[2]), "a.u.", re(&uz)),
            Field::scalar(format!("|{}|", names[2]), "a.u.", uz.iter().map(|z| z.norm()).collect()),
        ];
        if asm.is_piezo {
            let phi = comp(3);
            fields.push(Field::scalar("|phi|", "a.u.", phi.iter().map(|z| z.norm()).collect()));
            fields.push(Field::scalar("Re phi", "a.u.", re(&phi)));
        }
        let (layer_fraction, phys_mass_norm) = match (&m_all, &m_pml) {
            (Some(ma), Some(mp)) => {
                let (tot, in_layer) = (mass_norm(ma, u_mech), mass_norm(mp, u_mech));
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
        if asm.is_piezo {
            let (mut re_p, mut im_p) =
                (Vec::with_capacity(mesh.n_cells() * nloc), Vec::with_capacity(mesh.n_cells() * nloc));
            let mut dofs = Vec::new();
            for cell in 0..mesh.n_cells() {
                mesh.cell_dofs(cell, order, &mut dofs);
                for &d in dofs.iter().take(nloc) {
                    let z = u[3 * n1 + d] * rot * scale;
                    re_p.push(z.re);
                    im_p.push(z.im);
                }
            }
            cells.push_complex("phi", "a.u.", re_p, im_p);
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
        let u_mech = &u[..3 * n1];
        let big = u_mech.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);

        let n_comp = if asm.is_piezo { 4 } else { 3 };
        let mut re_blocks = Vec::with_capacity(n_comp * n1);
        let mut im_blocks = Vec::with_capacity(n_comp * n1);
        for k in 0..n_comp {
            for a in 0..n1 {
                let z = u[k * n1 + a] * rot * scale;
                re_blocks.push(z.re);
                im_blocks.push(z.im);
            }
        }
        let re_proj = proj.project_lagrange_blocks(mixed, dof_map, &re_blocks, n_comp, order);
        let im_proj = proj.project_lagrange_blocks(mixed, dof_map, &im_blocks, n_comp, order);
        let n_out = re_proj.len() / n_comp;
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
        let mut fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar("Re ux", "a.u.", re_proj[0..n_out].to_vec()),
            Field::scalar("Re uz", "a.u.", re_proj[n_out..2 * n_out].to_vec()),
            Field::scalar("Re uy", "a.u.", re_proj[2 * n_out..3 * n_out].to_vec()),
            Field::scalar("|uy|", "a.u.", uy_abs),
        ];
        if asm.is_piezo {
            let phi_abs: Vec<f64> = (0..n_out).map(|a| re_proj[3 * n_out + a].hypot(im_proj[3 * n_out + a])).collect();
            fields.push(Field::scalar("|phi|", "a.u.", phi_abs));
            fields.push(Field::scalar("Re phi", "a.u.", re_proj[3 * n_out..4 * n_out].to_vec()));
        }

        let mut cell_fields = Vec::with_capacity(n_comp);
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
        if asm.is_piezo {
            let mut re_pc = Vec::with_capacity(*m_offsets.last().unwrap());
            let mut im_pc = Vec::with_capacity(*m_offsets.last().unwrap());
            for mc in 0..mixed.n_cells() {
                for &gd in dof_map.cell_dofs(mc) {
                    re_pc.push(re_blocks[3 * n1 + gd]);
                    im_pc.push(im_blocks[3 * n1 + gd]);
                }
            }
            cell_fields.push(proj.project_cell_field(
                mixed,
                dof_map.space,
                "phi",
                "a.u.",
                &re_pc,
                &im_pc,
                &m_offsets,
                order,
            ));
        }

        let (layer_fraction, phys_mass_norm) = match (&m_all, &m_pml) {
            (Some(ma), Some(mp)) => {
                let (tot, in_layer) = (mass_norm(ma, u_mech), mass_norm(mp, u_mech));
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
        "{what}; {}complex {} eigenproblem, Krylov dimension {}{}",
        if asm.is_piezo { "coupled piezoelectric " } else { "" },
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

struct CrossPiezoSolveInput<'a> {
    kind: CrossSection,
    spec: &'a Problem,
    n1: usize,
    n_cells: usize,
    cell_labels: &'a [i32],
    mats: &'a [crate::piezo_mat::PiezoCellMat],
    mech_untouched: Vec<usize>,
    elec_untouched: Vec<usize>,
    pml_clamped: Vec<usize>,
    axis_mech_pairs: Vec<(usize, usize)>,
    axis_phi_fixed: Vec<usize>,
    axis_dofs: usize,
    rel_u: Vec<(usize, usize, C)>,
    rel_p_local: Vec<(usize, usize, C)>,
    lrb_labels: &'a [i32],
    pml_cells: Vec<bool>,
    diagnostics: Vec<String>,
    has_pml: bool,
    k_full: CCsr,
    m_uu: CCsr,
    c_uu_full: Option<CCsr>,
    n_modes: usize,
    f_search: f64,
    get_cell_dofs: &'a dyn Fn(usize, &mut Vec<usize>),
    get_bd: &'a dyn Fn(&[i32]) -> Vec<usize>,
}

fn solve_cross_section_piezo_common(inp: CrossPiezoSolveInput<'_>) -> Result<(CrossAssembly, CrossSolve), String> {
    let n1 = inp.n1;
    let n_mech = 3 * n1;
    let n_tot = 4 * n1;
    let spec = inp.spec;
    if matches!(inp.kind, CrossSection::Axisymmetric(m) if m.abs() > 0.0)
        && spec.bcs.iter().any(|b| {
            matches!(b.kind, BcKind::FloatingPotential { .. } | BcKind::Conductor { role: ConductorRole::Floating, .. })
        })
    {
        return Err("floating equipotential electrodes are only defined for axisymmetric azimuthal order m = 0".into());
    }

    let mut fixed_u_set: BTreeSet<usize> = inp.mech_untouched.into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = (inp.get_bd)(&bc.labels);
            for k in 0..3 {
                if fixed[PLANE_AXES[k]] {
                    for &dof in &bd {
                        fixed_u_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    for (k, dof) in inp.axis_mech_pairs {
        fixed_u_set.insert(dof + k * n1);
    }
    for &dof in &inp.pml_clamped {
        for k in 0..3 {
            fixed_u_set.insert(dof + k * n1);
        }
    }
    let fixed_u: Vec<usize> = fixed_u_set.into_iter().collect();

    // Electrical constraints on local phi DOFs `0..n1`
    let mut elec_active_dofs = vec![false; n1];
    let mut cell_buf = Vec::new();
    for ci in 0..inp.n_cells {
        if inp.mats[ci].elec_active {
            (inp.get_cell_dofs)(ci, &mut cell_buf);
            for &d in &cell_buf {
                elec_active_dofs[d] = true;
            }
        }
    }
    let cond_interface_dofs = |labels: &[i32]| -> Vec<usize> {
        let mut in_cond = vec![false; n1];
        let mut buf = Vec::new();
        for ci in 0..inp.n_cells {
            if labels.contains(&inp.cell_labels[ci]) {
                (inp.get_cell_dofs)(ci, &mut buf);
                for &d in &buf {
                    in_cond[d] = true;
                }
            }
        }
        (0..n1).filter(|&d| in_cond[d] && elec_active_dofs[d]).collect()
    };

    let elec_untouched_set: BTreeSet<usize> = inp.elec_untouched.iter().copied().collect();
    let mut fixed_p_set: BTreeSet<usize> = inp.elec_untouched.into_iter().collect();
    for &dof in &inp.pml_clamped {
        fixed_p_set.insert(dof);
    }
    for &dof in &inp.axis_phi_fixed {
        fixed_p_set.insert(dof);
    }

    let mut floating_groups: Vec<(String, Vec<usize>)> = Vec::new();
    for bc in &spec.bcs {
        match &bc.kind {
            BcKind::Terminal { .. } | BcKind::Conductor { role: ConductorRole::Terminal, .. } => {
                return Err("terminal boundaries are only valid in the driven study".into());
            }
            BcKind::Ground | BcKind::Potential(_) => {
                fixed_p_set.extend((inp.get_bd)(&bc.labels));
            }
            BcKind::FloatingPotential { .. } => {
                let label = format!("floating_potential {:?}", bc.labels);
                let bd = (inp.get_bd)(&bc.labels);
                if bd.is_empty() {
                    return Err(format!("floating potential '{label}' has no boundary degrees of freedom"));
                }
                floating_groups.push((label, bd));
            }
            BcKind::Conductor { role, name, .. } => {
                let iface = cond_interface_dofs(&bc.labels);
                match role {
                    ConductorRole::Ground | ConductorRole::Fixed => {
                        fixed_p_set.extend(iface);
                    }
                    ConductorRole::Floating => {
                        let label = if name.is_empty() { format!("conductor {:?}", bc.labels) } else { name.clone() };
                        if iface.is_empty() {
                            return Err(format!("floating conductor '{label}' has no active dielectric interface"));
                        }
                        floating_groups.push((label, iface));
                    }
                    ConductorRole::Terminal => unreachable!(),
                }
            }
            BcKind::Fixed { .. } | BcKind::Load(_) | BcKind::LowReflecting | BcKind::SurfaceCharge(_) => {}
        }
    }

    let mut rel_all: Vec<(usize, usize, C)> = inp.rel_u.clone();
    let mut adj_p: Vec<(usize, usize)> = Vec::new();
    for &(dst, src, w) in &inp.rel_p_local {
        rel_all.push((n_mech + dst, n_mech + src, w));
        adj_p.push((dst, src));
    }
    for (name, group) in &floating_groups {
        let root = group[0];
        if group.iter().any(|d| fixed_p_set.contains(d)) {
            return Err(format!("floating electrode '{name}' is shorted to a fixed or grounded potential"));
        }
        for &d in &group[1..] {
            rel_all.push((n_mech + d, n_mech + root, c(1.0, 0.0)));
            adj_p.push((d, root));
        }
    }

    let skip_gauge = match inp.kind {
        CrossSection::Waveguide(q) => q.abs() > 1e-30,
        CrossSection::Axisymmetric(m) => m.abs() > 0.0,
    };
    if !skip_gauge {
        let mut graph: Vec<Vec<usize>> = vec![Vec::new(); n1];
        for ci in 0..inp.n_cells {
            if !inp.mats[ci].elec_active {
                continue;
            }
            (inp.get_cell_dofs)(ci, &mut cell_buf);
            if let Some(&first) = cell_buf.first() {
                for &other in &cell_buf[1..] {
                    graph[first].push(other);
                    graph[other].push(first);
                }
            }
        }
        for &(a, b) in &adj_p {
            graph[a].push(b);
            graph[b].push(a);
        }
        let mut visited = vec![false; n1];
        let mut queue = VecDeque::new();
        for &seed in &fixed_p_set {
            if !visited[seed] {
                visited[seed] = true;
                queue.push_back(seed);
            }
        }
        let run_bfs = |visited: &mut [bool], queue: &mut VecDeque<usize>| {
            while let Some(u) = queue.pop_front() {
                for &v in &graph[u] {
                    if !visited[v] {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
        };
        run_bfs(&mut visited, &mut queue);
        for d in 0..n1 {
            if elec_active_dofs[d] && !elec_untouched_set.contains(&d) && !visited[d] {
                fixed_p_set.insert(d);
                visited[d] = true;
                queue.push_back(d);
                run_bfs(&mut visited, &mut queue);
            }
        }
    }

    let mut fixed_all_vec: Vec<usize> = fixed_u.clone();
    fixed_all_vec.extend(fixed_p_set.iter().map(|&d| n_mech + d));
    fixed_all_vec.sort_unstable();
    fixed_all_vec.dedup();

    let mut c6_out: Vec<M6> = Vec::with_capacity(inp.mats.len());
    let mut rho_out: Vec<f64> = Vec::with_capacity(inp.mats.len());
    for m in inp.mats {
        c6_out.push(m.c6);
        rho_out.push(m.rho);
    }
    let has_dielectric_loss = inp.mats.iter().any(|m| m.eps_im.iter().flatten().any(|&v| v != 0.0));
    let has_damping = inp.c_uu_full.is_some();
    let has_loss = inp.has_pml || !inp.lrb_labels.is_empty() || has_dielectric_loss || has_damping;
    let is_hermitian = !inp.has_pml && !has_dielectric_loss && !has_damping;

    let (free_out, red_out, n_u_red, n_tot_red, k_red, m_uu_red, c_uu_red) = if rel_all.is_empty() {
        let free_all = free_dofs(n_tot, &fixed_all_vec);
        let free_u: Vec<usize> = free_all.iter().copied().filter(|&d| d < n_mech).collect();
        let n_u_r = free_u.len();
        let n_tot_r = free_all.len();
        let k_r = inp.k_full.submatrix(&free_all, &free_all);
        let m_r = inp.m_uu.submatrix(&free_u, &free_u);
        let c_r = inp.c_uu_full.as_ref().map(|cm| cm.submatrix(&free_u, &free_u));
        (free_all, None, n_u_r, n_tot_r, k_r, m_r, c_r)
    } else {
        let red_u = femcore::periodic::Reduction::try_build(n_mech, &inp.rel_u, &fixed_u).map_err(|e| e.to_string())?;
        let red_all =
            femcore::periodic::Reduction::try_build(n_tot, &rel_all, &fixed_all_vec).map_err(|e| e.to_string())?;
        let n_u_r = red_u.n_reduced();
        let n_tot_r = red_all.n_reduced();
        if red_all.roots[..n_u_r] != red_u.roots[..] {
            return Err("internal error: periodic/floating reduction broke (u, phi) partition".into());
        }
        let pm_u = red_u.prolong();
        let pm_all = red_all.prolong();
        let k_r = femcore::periodic::congruence(&pm_all, &inp.k_full);
        let m_r = femcore::periodic::congruence(&pm_u, &inp.m_uu);
        let c_r = inp.c_uu_full.as_ref().map(|cm| femcore::periodic::congruence(&pm_u, cm));
        (Vec::new(), Some(red_all), n_u_r, n_tot_r, k_r, m_r, c_r)
    };

    let k_pp_red = crate::piezo_asm::sub_csr(&k_red, n_u_red, n_tot_red, n_u_red, n_tot_red);
    let k_pu_red = crate::piezo_asm::sub_csr(&k_red, n_u_red, n_tot_red, 0, n_u_red);

    let solved = crate::piezo::solve_piezo_complex_schur(&crate::piezo::PiezoComplexEigenInput {
        n_u: n_u_red,
        n_tot: n_tot_red,
        n_modes: inp.n_modes,
        f_search: inp.f_search,
        is_hermitian,
        has_pml: inp.has_pml,
        has_dielectric_loss,
        n_pairs: spec.periodic.len(),
        k_red: &k_red,
        m_uu_red: &m_uu_red,
        c_uu_red: c_uu_red.as_ref(),
        k_pp_red: &k_pp_red,
        k_pu_red: &k_pu_red,
    })?;

    let full_vectors: Vec<Vec<C>> = solved
        .reduced_vecs
        .iter()
        .map(|xr| {
            if let Some(red) = &red_out {
                red.expand(xr)
            } else {
                let mut full = vec![c(0.0, 0.0); n_tot];
                for (a, &dof) in free_out.iter().enumerate() {
                    full[dof] = xr[a];
                }
                full
            }
        })
        .collect();

    let (CrossSection::Axisymmetric(m_az) | CrossSection::Waveguide(m_az)) = inp.kind;
    let asm = CrossAssembly {
        kind: inp.kind,
        order: spec.order,
        n: n_tot,
        n1,
        kf: k_red,
        mf: m_uu_red,
        cf: None,
        free: free_out,
        fixed: fixed_all_vec,
        c6: c6_out.clone(),
        rho: rho_out.clone(),
        n_modes: inp.n_modes,
        f_search: inp.f_search,
        has_loss,
        axis_dofs: inp.axis_dofs,
        pml_cells: inp.pml_cells,
        red: red_out,
        diagnostics: inp.diagnostics,
        is_piezo: true,
    };
    let csolve = CrossSolve {
        modes: CrossModes {
            m_az,
            freqs: solved.freqs,
            vectors: full_vectors,
            c6: c6_out,
            rho: rho_out,
            converged: solved.converged,
        },
        dim: solved.dim,
    };
    Ok((asm, csolve))
}

fn solve_waveguide_piezo(mesh: &Discretisation, spec: &Problem, q: f64) -> Result<SolveResult, String> {
    if mesh.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let t0 = Instant::now();
    let Study::Modes { n_modes, f_search, .. } = spec.study else {
        return Err("solve_waveguide_piezo needs the modes study".into());
    };
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let (_by_dom, mats) = crate::piezo_mat::resolve_piezo_mats_simplex(mesh, spec, true)?;
    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_simplex(mesh, spec, false, f_search)?)
    } else {
        None
    };

    let n1 = mesh.n_dofs(order);
    let cond_domains = spec.conductor_domains();
    let mech_untouched = untouched_dofs(mesh, order, 3, &spec.inactive);
    let elec_untouched = untouched_dofs(mesh, order, 1, &cond_domains);
    let pml_clamped =
        resolved.as_ref().map_or_else(Vec::new, |r| outer_clamped_dofs_simplex(mesh, order, spec, r, false, 0.0));

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
    if !lrb_labels.is_empty() && !spec.periodic.is_empty() {
        return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
    }

    let mut rel_u = Vec::new();
    let mut rel_p_local = Vec::new();
    for p in &spec.periodic {
        let mt = femcore::periodic::match_vertices(mesh, p)?;
        rel_u.extend(femcore::periodic::nodal_relations(mesh, order, 3, &mt, p.phase)?);
        rel_p_local.extend(femcore::periodic::nodal_relations(mesh, order, 1, &mt, p.phase)?);
    }

    let (k_full, m_uu, c_eta) =
        crate::waveguide::assemble_waveguide_piezo_simplex(mesh, &mats, order, q, resolved.as_ref())?;
    let c_lrb = if lrb_labels.is_empty() {
        None
    } else {
        let c6_raw: Vec<M6> = mats.iter().map(|m| m.c6).collect();
        let rho_raw: Vec<f64> = mats.iter().map(|m| m.rho).collect();
        Some(to_complex(&crate::lrb::damping_matrix(mesh, &c6_raw, &rho_raw, &lrb_labels, order, 3)))
    };
    let c_uu_full = match (c_eta, c_lrb) {
        (Some(ce), Some(cl)) => Some(femcore::linalg::csr::Csr::axpby(c(1.0, 0.0), &ce, c(1.0, 0.0), &cl)),
        (Some(ce), None) => Some(ce),
        (None, Some(cl)) => Some(cl),
        (None, None) => None,
    };

    let pml_cells: Vec<bool> =
        (0..mesh.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mesh.cell_labels[ci]))).collect();
    let has_pml = pml_cells.iter().any(|&p| p);
    let diagnostics = resolved.map_or_else(Vec::new, |r| r.diagnostics);

    let get_cell_dofs = |ci: usize, out: &mut Vec<usize>| mesh.cell_dofs(ci, order, out);
    let get_bd = |labels: &[i32]| boundary_dofs(mesh, labels, order);
    let (asm, solved) = solve_cross_section_piezo_common(CrossPiezoSolveInput {
        kind: CrossSection::Waveguide(q),
        spec,
        n1,
        n_cells: mesh.n_cells(),
        cell_labels: &mesh.cell_labels,
        mats: &mats,
        mech_untouched,
        elec_untouched,
        pml_clamped,
        axis_mech_pairs: Vec::new(),
        axis_phi_fixed: Vec::new(),
        axis_dofs: 0,
        rel_u,
        rel_p_local,
        lrb_labels: &lrb_labels,
        pml_cells,
        diagnostics,
        has_pml,
        k_full,
        m_uu,
        c_uu_full,
        n_modes,
        f_search,
        get_cell_dofs: &get_cell_dofs,
        get_bd: &get_bd,
    })?;
    let mut result = cross_section_fields(mesh, &asm, &solved)?;
    result.seconds = t0.elapsed().as_secs_f64();
    Ok(result)
}

fn solve_waveguide_piezo_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
    q: f64,
) -> Result<SolveResult, String> {
    if mixed.dim != 2 {
        return Err("a waveguide cross-section model needs a 2D mesh".into());
    }
    let t0 = Instant::now();
    let Study::Modes { n_modes, f_search, .. } = spec.study else {
        return Err("solve_waveguide_piezo_mixed needs the modes study".into());
    };
    let order = spec.order;
    crate::solid::mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let n1 = dof_map.n_dofs;
    let (_by_dom, mats) = crate::piezo_mat::resolve_piezo_mats_mixed(mixed, spec, true)?;
    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_mixed(mixed, spec, false, f_search)?)
    } else {
        None
    };

    let cond_domains = spec.conductor_domains();
    let mech_untouched = crate::solid::mixed_untouched_dofs(mixed, &dof_map, 3, &spec.inactive);
    let elec_untouched = crate::solid::mixed_untouched_dofs(mixed, &dof_map, 1, &cond_domains);
    let pml_clamped =
        resolved.as_ref().map_or_else(Vec::new, |r| outer_clamped_dofs_mixed(mixed, &dof_map, spec, r, false, 0.0));

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
    if !lrb_labels.is_empty() && !spec.periodic.is_empty() {
        return Err("low-reflecting boundaries cannot be combined with periodic conditions yet".into());
    }

    let mut rel_u = Vec::new();
    let mut rel_p_local = Vec::new();
    for p in &spec.periodic {
        rel_u.extend(mixed_nodal_relations(mixed, &dof_map, 3, p)?);
        rel_p_local.extend(mixed_nodal_relations(mixed, &dof_map, 1, p)?);
    }

    let (k_full, m_uu, c_eta) =
        crate::waveguide::assemble_waveguide_piezo_mixed(mixed, &dof_map, &mats, q, resolved.as_ref())?;
    let c_lrb = if lrb_labels.is_empty() {
        None
    } else {
        let c6_raw: Vec<M6> = mats.iter().map(|m| m.c6).collect();
        let rho_raw: Vec<f64> = mats.iter().map(|m| m.rho).collect();
        Some(to_complex(&crate::lrb::damping_matrix_mixed(mixed, &dof_map, &c6_raw, &rho_raw, &lrb_labels, 3, 1.0)))
    };
    let c_uu_full = match (c_eta, c_lrb) {
        (Some(ce), Some(cl)) => Some(femcore::linalg::csr::Csr::axpby(c(1.0, 0.0), &ce, c(1.0, 0.0), &cl)),
        (Some(ce), None) => Some(ce),
        (None, Some(cl)) => Some(cl),
        (None, None) => None,
    };

    let pml_cells: Vec<bool> =
        (0..mixed.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mixed.cell_labels[ci]))).collect();
    let has_pml = pml_cells.iter().any(|&p| p);
    let diagnostics = resolved.map_or_else(Vec::new, |r| r.diagnostics);

    let get_cell_dofs = |ci: usize, out: &mut Vec<usize>| {
        out.clear();
        out.extend_from_slice(dof_map.cell_dofs(ci));
    };
    let get_bd = |labels: &[i32]| dof_map.boundary_dofs(mixed, labels);
    let (asm, solved) = solve_cross_section_piezo_common(CrossPiezoSolveInput {
        kind: CrossSection::Waveguide(q),
        spec,
        n1,
        n_cells: mixed.n_cells(),
        cell_labels: &mixed.cell_labels,
        mats: &mats,
        mech_untouched,
        elec_untouched,
        pml_clamped,
        axis_mech_pairs: Vec::new(),
        axis_phi_fixed: Vec::new(),
        axis_dofs: 0,
        rel_u,
        rel_p_local,
        lrb_labels: &lrb_labels,
        pml_cells,
        diagnostics,
        has_pml,
        k_full,
        m_uu,
        c_uu_full,
        n_modes,
        f_search,
        get_cell_dofs: &get_cell_dofs,
        get_bd: &get_bd,
    })?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let mut result = cross_section_fields_mixed(mixed, &dof_map, &proj, &asm, &solved)?;
    result.seconds = t0.elapsed().as_secs_f64();
    Ok(result)
}

fn solve_axisymmetric_piezo(mesh: &Discretisation, spec: &Problem, m: f64) -> Result<SolveResult, String> {
    if mesh.dim != 2 {
        return Err("an axisymmetric model needs a 2D mesh".into());
    }
    let t0 = Instant::now();
    let Study::AxisymmetricEigenfrequency { n_modes, f_search, .. } = spec.study else {
        return Err("solve_axisymmetric_piezo needs the axisymmetric_eigenfrequency study".into());
    };
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let (_by_dom, mut mats) = crate::piezo_mat::resolve_piezo_mats_simplex(mesh, spec, true)?;
    let mut diagnostics = crate::axisym::average_axisym_piezo_mats(&mut mats, &mesh.cell_labels)?;
    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_simplex(mesh, spec, true, f_search)?)
    } else {
        None
    };

    let n1 = mesh.n_dofs(order);
    let cond_domains = spec.conductor_domains();
    let mech_untouched = untouched_dofs(mesh, order, 3, &spec.inactive);
    let elec_untouched = untouched_dofs(mesh, order, 1, &cond_domains);
    let axis_tol = crate::axisym::axis_tolerance(mesh);
    let pml_clamped =
        resolved.as_ref().map_or_else(Vec::new, |r| outer_clamped_dofs_simplex(mesh, order, spec, r, true, axis_tol));
    let (axis_mech_pairs, axis_phi_fixed, axis_dofs) = crate::axisym::axis_constraints_piezo(mesh, order, m)?;

    let (k_full, m_uu, c_uu_full) =
        crate::axisym::assemble_axisym_piezo_simplex(mesh, &mats, order, m, resolved.as_ref())?;
    let pml_cells: Vec<bool> =
        (0..mesh.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mesh.cell_labels[ci]))).collect();
    let has_pml = pml_cells.iter().any(|&p| p);
    if let Some(res) = resolved {
        diagnostics.extend(res.diagnostics);
    }

    let get_cell_dofs = |ci: usize, out: &mut Vec<usize>| mesh.cell_dofs(ci, order, out);
    let get_bd = |labels: &[i32]| boundary_dofs(mesh, labels, order);
    let (asm, solved) = solve_cross_section_piezo_common(CrossPiezoSolveInput {
        kind: CrossSection::Axisymmetric(m),
        spec,
        n1,
        n_cells: mesh.n_cells(),
        cell_labels: &mesh.cell_labels,
        mats: &mats,
        mech_untouched,
        elec_untouched,
        pml_clamped,
        axis_mech_pairs,
        axis_phi_fixed,
        axis_dofs,
        rel_u: Vec::new(),
        rel_p_local: Vec::new(),
        lrb_labels: &[],
        pml_cells,
        diagnostics,
        has_pml,
        k_full,
        m_uu,
        c_uu_full,
        n_modes,
        f_search,
        get_cell_dofs: &get_cell_dofs,
        get_bd: &get_bd,
    })?;
    let mut result = cross_section_fields(mesh, &asm, &solved)?;
    result.seconds = t0.elapsed().as_secs_f64();
    Ok(result)
}

fn solve_axisymmetric_piezo_mixed(
    raw_mesh: &Mesh,
    mixed: &MixedDiscretisation,
    spec: &Problem,
    m: f64,
) -> Result<SolveResult, String> {
    if mixed.dim != 2 {
        return Err("an axisymmetric model needs a 2D mesh".into());
    }
    let t0 = Instant::now();
    let Study::AxisymmetricEigenfrequency { n_modes, f_search, .. } = spec.study else {
        return Err("solve_axisymmetric_piezo_mixed needs the axisymmetric_eigenfrequency study".into());
    };
    let order = spec.order;
    crate::solid::mixed_check_boundary_labels(mixed, &spec.bcs)?;

    let space = if order == 1 { LagrangeSpace::Order1 } else { LagrangeSpace::Order2Serendipity };
    let dof_map = MixedLagrangeDofMap::new(mixed, space);
    let n1 = dof_map.n_dofs;
    let (_by_dom, mut mats) = crate::piezo_mat::resolve_piezo_mats_mixed(mixed, spec, true)?;
    let mut diagnostics = crate::axisym::average_axisym_piezo_mats(&mut mats, &mixed.cell_labels)?;
    let resolved = if spec.materials.has_pml() {
        Some(ResolvedPmlDomains::from_mixed(mixed, spec, true, f_search)?)
    } else {
        None
    };

    let cond_domains = spec.conductor_domains();
    let mech_untouched = crate::solid::mixed_untouched_dofs(mixed, &dof_map, 3, &spec.inactive);
    let elec_untouched = crate::solid::mixed_untouched_dofs(mixed, &dof_map, 1, &cond_domains);
    let axis_tol = crate::axisym::axis_tolerance_mixed(mixed);
    let pml_clamped =
        resolved.as_ref().map_or_else(Vec::new, |r| outer_clamped_dofs_mixed(mixed, &dof_map, spec, r, true, axis_tol));
    let (axis_mech_pairs, axis_phi_fixed, axis_dofs) = crate::axisym::axis_constraints_piezo_mixed(mixed, &dof_map, m)?;

    let (k_full, m_uu, c_uu_full) =
        crate::axisym::assemble_axisym_piezo_mixed(mixed, &dof_map, &mats, m, resolved.as_ref())?;
    let pml_cells: Vec<bool> =
        (0..mixed.n_cells()).map(|ci| resolved.as_ref().is_some_and(|r| r.contains(mixed.cell_labels[ci]))).collect();
    let has_pml = pml_cells.iter().any(|&p| p);
    if let Some(res) = resolved {
        diagnostics.extend(res.diagnostics);
    }

    let get_cell_dofs = |ci: usize, out: &mut Vec<usize>| {
        out.clear();
        out.extend_from_slice(dof_map.cell_dofs(ci));
    };
    let get_bd = |labels: &[i32]| dof_map.boundary_dofs(mixed, labels);
    let (asm, solved) = solve_cross_section_piezo_common(CrossPiezoSolveInput {
        kind: CrossSection::Axisymmetric(m),
        spec,
        n1,
        n_cells: mixed.n_cells(),
        cell_labels: &mixed.cell_labels,
        mats: &mats,
        mech_untouched,
        elec_untouched,
        pml_clamped,
        axis_mech_pairs,
        axis_phi_fixed,
        axis_dofs,
        rel_u: Vec::new(),
        rel_p_local: Vec::new(),
        lrb_labels: &[],
        pml_cells,
        diagnostics,
        has_pml,
        k_full,
        m_uu,
        c_uu_full,
        n_modes,
        f_search,
        get_cell_dofs: &get_cell_dofs,
        get_bd: &get_bd,
    })?;
    let (_simplex, proj) = MixedToSimplexProjector::from_mesh(raw_mesh, mixed).map_err(|e| e.to_string())?;
    let mut result = cross_section_fields_mixed(mixed, &dof_map, &proj, &asm, &solved)?;
    result.seconds = t0.elapsed().as_secs_f64();
    Ok(result)
}
