// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The shared kernel of the elastic cross-section models: stiffness and mass on a 2D mesh with three displacement components, the
//! constraints, absorbing layers and low-reflecting boundaries, the complex eigensolve and the fields. It serves the axisymmetric
//! model (`axisym`: bodies of revolution, azimuthal order `m`); it was `femcore::sms_kz`, the shared part of the removed out-of-plane
//! wave-number model (paul, 2026-11-04), and moved here in restructure phase 4 (ADR 42). A library entry (`solve_solid_axisym`), not
//! reachable through the envelope yet: the ring checks (`tests/axisym_closed_form.rs`, `tests/ring_frame.rs`) and the golden `ring_model_axes` use it.

use crate::problem::*;
use femcore::fem::check_boundary_labels;
use femcore::fem::*;
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCsr, C};
use femcore::material::M6;
use femcore::props::PmlDir;
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::time::Instant;

fn check_labels(mesh: &Discretisation, bcs: &[Bc]) -> Result<(), String> {
    check_boundary_labels(mesh, bcs.iter().flat_map(|b| b.labels.iter().copied()))
}

fn ccsr_block_diag(blocks: &[&CCsr]) -> CCsr {
    let n: usize = blocks.iter().map(|b| b.nrows()).sum();
    let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
    let mut off = 0;
    for b in blocks {
        for r in 0..b.nrows() {
            for k in b.indptr()[r]..b.indptr()[r + 1] {
                indices.push(b.indices()[k] as usize + off);
                data.push(b.values()[k]);
            }
            indptr.push(indices.len());
        }
        off += b.ncols();
    }
    CCsr::from_raw(n, off, indptr, indices, data)
}

fn to_complex(a: &Csr) -> CCsr {
    a.to_complex()
}

/// The wafer axis (0 = x, 1 = y, 2 = z) of each displacement block: x, z, y (for a ring r, z, theta).
const PLANE_AXES: [usize; 3] = [0, 2, 1];

/// The complex eigenvectors behind a cross-section solid result: displacement over all `3 * n_dofs(order)` dofs (component blocks), the complex
/// frequencies in Hz, in the order of the solutions.
pub struct CrossModes {
    pub m_az: f64,
    pub freqs: Vec<C>,
    pub vectors: Vec<Vec<C>>,
    pub c6: Vec<M6>,
    pub rho: Vec<f64>,
    pub converged: bool,
}

/// What the third displacement component and the third dimension are: the azimuth of an axisymmetric body with azimuthal order `m`
/// (mesh coordinates x = r, y = z), or the longitudinal propagation axis `+y` of a 2.5D waveguide cross-section with wavevector `q`
/// (rad/m, displacement `exp(-i q y)`).
#[derive(Clone, Copy, Debug)]
pub enum CrossSection {
    Axisymmetric(f64),
    Waveguide(f64),
}

/// Elastic modes of a body of revolution with azimuthal dependence `exp(-i m theta)` on the (r, z) half
/// plane (mesh x = r, y = z): the displacement components are (u_r, u_z, u_theta) in this order.
pub fn solve_solid_axisym(mesh: &Discretisation, spec: &Problem, m: f64) -> Result<(SolveResult, CrossModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Axisymmetric(m))
}

/// Elastic modes of a 2D `(x, z)` waveguide cross-section with longitudinal dependence `exp(-i q y)`:
/// the displacement components are `(u_x, u_z, u_y)` in this order.
pub fn solve_solid_waveguide(
    mesh: &Discretisation,
    spec: &Problem,
    q: f64,
) -> Result<(SolveResult, CrossModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Waveguide(q))
}

/// The axisymmetric eigenfrequency study (ADR 68): the azimuthal order is taken from the study, the result carries `physics` and `study`.
/// Refuses what the study cannot say (`Problem::check_axisymmetric`) and a mesh that leaves the half plane `x >= 0`; with absorbing layers the
/// frequencies are complex (`Q` and `layer fraction` scalars, `eigenvalue` on every mode).
pub fn solve_axisymmetric(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let Study::AxisymmetricEigenfrequency { azimuthal_order, .. } = problem.study else {
        return Err("solve_axisymmetric needs the axisymmetric_eigenfrequency study".into());
    };
    problem.check_axisymmetric().map_err(|e| e.to_string())?;
    check_half_plane(mesh)?;
    let (mut result, _) = solve_cross_section(mesh, problem, CrossSection::Axisymmetric(azimuthal_order as f64))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    // the study reads neither the plane assumption nor the thickness: say so when one was set
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("axisymmetric study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

/// The 2.5D elastic waveguide modes study: the longitudinal wavevector `q` is taken from the study, the result carries `physics` and `study`.
/// Refuses what the study cannot say (`Problem::check_modes`); with absorbing layers or low-reflecting boundaries the frequencies are complex
/// (`Q` and `layer fraction` scalars, `eigenvalue` on every mode).
pub fn solve_waveguide_modes(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let Study::Modes { wavevector, .. } = problem.study else {
        return Err("solve_waveguide_modes needs the modes study".into());
    };
    problem.check_modes().map_err(|e| e.to_string())?;
    let (mut result, _) = solve_cross_section(mesh, problem, CrossSection::Waveguide(wavevector))?;
    result.physics = "sms".into();
    result.study = problem.study.kind().into();
    // the study reads neither the plane assumption nor the thickness: say so when one was set
    if problem.plane != Plane::Strain || problem.thickness != 1.0 {
        result.log.push("waveguide modes study: the 2D parameters plane and thickness are not used".into());
    }
    Ok(result)
}

/// The (r, z) half plane lies in `x >= 0` (up to 1e-9 of the width of the mesh).
fn check_half_plane(mesh: &Discretisation) -> Result<(), String> {
    let xs = (0..mesh.n_points()).map(|i| mesh.point(i)[0]);
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

/// Stage 1 of the cross-section solid kernels: stiffness, mass and (low-reflecting) damping on the free dofs, the
/// constraints, and the material data the coupling integrals need.
pub struct CrossAssembly {
    pub kind: CrossSection,
    pub order: usize,
    /// all dofs (three displacement blocks of `n1`)
    pub n: usize,
    pub n1: usize,
    pub kf: CCsr,
    pub mf: CCsr,
    /// damping on the free dofs when a boundary is low-reflecting
    pub cf: Option<Csr>,
    pub free: Vec<usize>,
    pub fixed: Vec<usize>,
    pub c6: Vec<M6>,
    pub rho: Vec<f64>,
    pub n_modes: usize,
    pub f_search: f64,
    /// a PML or a low-reflecting boundary makes the frequencies complex
    pub has_loss: bool,
    /// dofs of the mesh on the axis `r = 0` (nodes, not components) whose regularity conditions are in `fixed`
    pub axis_dofs: usize,
    /// per cell: the cell belongs to a perfectly matched layer (for the `layer fraction` scalar)
    pub pml_cells: Vec<bool>,
    /// Floquet-periodic reduction when lateral periodic boundary pairs are active in a waveguide model
    pub red: Option<femcore::periodic::Reduction>,
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
    // the blocks of the displacement are (x, z, y) (for a ring: r, z, theta): the stiffness in the Voigt order of those axes
    c6 = c6.iter().map(|c| femcore::section::plane_stiffness(c, 2)).collect();
    if matches!(kind, CrossSection::Axisymmetric(_)) {
        // a crystal's axes stay fixed while the local frame turns around the ring: the stiffness of the axisymmetric
        // model is the average over rotations about the ring axis (the y direction of the (r, z) plane); exact for isotropic materials
        for m in c6.iter_mut() {
            *m = crate::axisym::average_about_axis(m, 1, 16);
        }
    }
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    // fixed boundaries (zero displacement components)
    let mut fixed_set: std::collections::BTreeSet<usize> =
        untouched_dofs(mesh, order, 3, &spec.inactive).into_iter().collect();
    for bc in &spec.bcs {
        if let BcKind::Fixed { fixed, .. } = &bc.kind {
            let bd = boundary_dofs(mesh, &bc.labels, order);
            for k in 0..3 {
                // block k is wafer axis `PLANE_AXES[k]`
                if fixed[PLANE_AXES[k]] {
                    for &dof in &bd {
                        fixed_set.insert(dof + k * n1);
                    }
                }
            }
        }
    }
    // the axis of a body of revolution: regularity conditions on the dofs at r = 0 (axisym header, ADR 69)
    let axis_dofs = match kind {
        CrossSection::Axisymmetric(m_az) => {
            let (axis_pairs, axis_dofs) = crate::axisym::axis_constraints(mesh, order, m_az)?;
            fixed_set.extend(axis_pairs.into_iter().map(|(k, dof)| dof + k * n1));
            axis_dofs
        }
        CrossSection::Waveguide(_) => 0,
    };
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let stretch = if spec.materials.has_pml() {
        Some(match kind {
            CrossSection::Axisymmetric(_) => crate::axisym::cylindrical_stretches(mesh, &spec.materials.pml),
            CrossSection::Waveguide(_) => femcore::pml::cell_stretches(mesh, &spec.materials.pml),
        })
    } else {
        None
    };
    let k = match kind {
        CrossSection::Axisymmetric(m) => {
            crate::axisym::assemble_stiffness_axisym(mesh, &c6, order, m, stretch.as_deref())?
        }
        CrossSection::Waveguide(q) => {
            crate::waveguide::assemble_stiffness_waveguide(mesh, &c6, order, q, stretch.as_deref())?
        }
    };
    let rho = crate::solid::cell_rho(&spec.materials, mesh, &spec.inactive)?;
    // mass: weighted by the determinant of the stretch in the layers (complex there)
    let mc = match &stretch {
        _ if matches!(kind, CrossSection::Axisymmetric(_)) => {
            crate::axisym::assemble_mass_axisym(mesh, &rho, order, stretch.as_deref())
        }
        None => to_complex(&Csr::block_diag(&[
            &assemble_mass(mesh, order, &rho),
            &assemble_mass(mesh, order, &rho),
            &assemble_mass(mesh, order, &rho),
        ])),
        Some(st) => {
            let w: Vec<C> = (0..mesh.n_cells()).map(|i| st[i][0] * st[i][1]).collect();
            let re_part: Vec<f64> = rho.iter().zip(&w).map(|(r, z)| r * z.re).collect();
            let im_part: Vec<f64> = rho.iter().zip(&w).map(|(r, z)| r * z.im).collect();
            let (mr, mi) = (assemble_mass(mesh, order, &re_part), assemble_mass(mesh, order, &im_part));
            let m1 = |a: &Csr, b: &Csr| -> CCsr { to_complex(a).add_scaled(&to_complex(b), c(0.0, 1.0)) };
            let one = m1(&mr, &mi);
            let blocks: Vec<&CCsr> = vec![&one, &one, &one];
            ccsr_block_diag(&blocks)
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
    let has_loss = stretch.is_some() || !lrb_labels.is_empty();
    let pml_cells: Vec<bool> = (0..mesh.n_cells())
        .map(|ci| {
            spec.materials
                .pml
                .get(&mesh.cell_labels[ci])
                .is_some_and(|p| p.x != PmlDir::Off || p.y != PmlDir::Off || p.z != PmlDir::Off)
        })
        .collect();
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
    })
}

/// Stage 2 output: the mode shapes over all dofs and what the eigensolver says about them.
pub struct CrossSolve {
    pub modes: CrossModes,
    /// Krylov dimension of the eigensolver
    pub dim: usize,
}

/// Stage 2: the complex eigenpairs (`K - w^2 M` or the quadratic problem with the low-reflecting damping) nearest the
/// search frequency, as complex frequencies in Hz with their displacement over all dofs, sorted by real part.
pub fn solve_cross_section_modes(asm: &CrossAssembly) -> Result<CrossSolve, String> {
    let (CrossSection::Axisymmetric(m_az) | CrossSection::Waveguide(m_az)) = asm.kind;
    let (n_modes, f_search, free, n) = (asm.n_modes, asm.f_search, &asm.free, asm.n);
    let (kf, mf) = (&asm.kf, &asm.mf);
    // the shift: the requested frequency, or just below the spectrum for the lowest modes
    let sigma = if f_search > 0.0 { (2.0 * std::f64::consts::PI * f_search).powi(2) } else { -1.0 };
    // (frequencies, displacement vectors over the free dofs, converged, Krylov dimension)
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

/// `u^H M u` of a real-weighted block-diagonal mass matrix (the complex matrix `assemble_mass_axisym` gives without a stretch).
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

/// Stage 3: the plot fields `|u|`, `Re ux` (radial/lateral), `Re uz` (axial/vertical), `Re uy` (azimuthal/longitudinal), `|uy|` of every mode (largest component
/// real, unit maximum), the displacement of the (x, z) blocks, and the un-averaged per-cell complex components `ux`, `uz`, `uy`.
/// Every solution carries `eigenvalue = [Re f, Im f]` (Hz); the scalar `layer fraction N` is the share of kinetic energy that lies in
/// absorbing layers (0 without layers), which tells a resonator/guided mode from a layer mode. What identifies a mode (ADR 80): `f re N (Hz)` and
/// `f im N (Hz)`, always numbers; component shares (`radial`/`axial`/`azimuthal` or `lateral`/`vertical`/`longitudinal`), the split of kinetic energy outside the layers
/// over the three components (sum 1); `share in domain <label> N`, the same integral per domain outside the layers (sum 1).
pub fn cross_section_fields(
    mesh: &Discretisation,
    asm: &CrossAssembly,
    solved: &CrossSolve,
) -> Result<SolveResult, String> {
    let (n1, order, has_loss) = (asm.n1, asm.order, asm.has_loss);
    let modes = &solved.modes;
    let mut solutions = vec![];
    let mut scalars = vec![];
    // mass weights over all cells and over the cells of the layers only
    let any_pml = asm.pml_cells.iter().any(|&p| p);
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
    // what identifies a mode (ADR 80): the split of kinetic energy outside the layers over the three components and over the domains
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
    for (j, &fz) in modes.freqs.iter().enumerate() {
        let freq = fz.re;
        let u = &modes.vectors[j];
        // rotate the phase so the largest component is real, scale to unit maximum modulus
        let big = u.iter().cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let scale = 1.0 / big.norm().max(1e-300);
        let comp = |k: usize| -> Vec<C> { (0..n1).map(|a| u[k * n1 + a] * rot * scale).collect() };
        let (ux, uy, uz) = (comp(0), comp(1), comp(2));
        let mag: Vec<f64> = (0..n1).map(|a| (ux[a].norm_sqr() + uy[a].norm_sqr() + uz[a].norm_sqr()).sqrt()).collect();
        let re = |v: &[C]| -> Vec<f64> { v.iter().map(|z| z.re).collect() };
        // the names of the three blocks: wafer (x, z, y); for a ring (r, z, theta) = (x, z, y)
        let names = ["ux", "uz", "uy"];
        // the node fields say `Re` because the modes are complex; the cell fields below keep the complex values under the plain names
        let fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar(format!("Re {}", names[0]), "a.u.", re(&ux)),
            Field::scalar(format!("Re {}", names[1]), "a.u.", re(&uy)),
            Field::scalar(format!("Re {}", names[2]), "a.u.", re(&uz)),
            Field::scalar(format!("|{}|", names[2]), "a.u.", uz.iter().map(|z| z.norm()).collect()),
        ];
        let layer_fraction = match (&m_all, &m_pml) {
            (Some(ma), Some(mp)) => {
                let (tot, in_layer) = (mass_norm(ma, u), mass_norm(mp, u));
                if tot > 0.0 {
                    (in_layer / tot).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            }
            _ => 0.0,
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
            let q = fz.re / (2.0 * fz.im.abs());
            solutions.push(Solution {
                label: format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, fz.re, fz.im, q),
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
        // numeric parts of the frequency, with or without loss (the text scalar `fN (Hz)` of a lossy mode stays, contract 1.1)
        scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(fz.re)));
        scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(fz.im)));
        let total: f64 = energy[j].iter().map(|e| e[0] + e[1] + e[2]).sum();
        let share = |x: f64| if total > 0.0 { x / total } else { 0.0 };
        for (k, name) in share_names.iter().enumerate() {
            let part: f64 = energy[j].iter().map(|e| e[k]).sum();
            scalars.push((format!("{name} share {}", j + 1), Scalar::Num(share(part))));
        }
        // the kernel knows domain numbers; the engine puts the region names in their place (`share in <region> N`)
        for (d, label) in domains.iter().enumerate() {
            let e = energy[j][d];
            scalars.push((format!("share in domain {label} {}", j + 1), Scalar::Num(share(e[0] + e[1] + e[2]))));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let m_az = modes.m_az;
    let what = match asm.kind {
        CrossSection::Axisymmetric(_) => format!("axisymmetric, azimuthal order m = {m_az}"),
        CrossSection::Waveguide(q) => format!("waveguide modes, longitudinal wavevector q = {q} rad/m"),
    };
    let mut log = vec![format!(
        "{what}; complex Hermitian eigenproblem, Krylov dimension {}{}",
        solved.dim,
        if modes.converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
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
    Ok(SolveResult { order: order as u8, dim: 2, solutions, scalars, log, dofs: asm.n, ..Default::default() })
}
