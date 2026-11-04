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
/// (mesh coordinates x = r, y = z).
#[derive(Clone, Copy, Debug)]
pub enum CrossSection {
    Axisymmetric(f64),
}

/// Elastic modes of a body of revolution with azimuthal dependence `exp(-i m theta)` on the (r, z) half
/// plane (mesh x = r, y = z): the displacement components are (u_r, u_z, u_theta) in this order.
pub fn solve_solid_axisym(mesh: &Discretisation, spec: &Problem, m: f64) -> Result<(SolveResult, CrossModes), String> {
    solve_cross_section(mesh, spec, CrossSection::Axisymmetric(m))
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
}

pub fn assemble_cross_section(
    mesh: &Discretisation,
    spec: &Problem,
    kind: CrossSection,
) -> Result<CrossAssembly, String> {
    if mesh.dim != 2 {
        return Err("a cross-section model needs a 2D mesh".into());
    }
    let Study::Eigenfrequency { n_modes, f_search } = spec.study.clone() else {
        return Err("the cross-section solid models need an Eigenfrequency study".into());
    };
    let order = spec.order;
    check_labels(mesh, &spec.bcs)?;
    let mut c6 = spec.materials.cell_stiffness_skipping(mesh, &spec.inactive)?;
    // the blocks of the displacement are (x, z, y) (for a ring: r, z, theta): the stiffness in the Voigt order of those axes
    c6 = c6.iter().map(|c| femcore::section::plane_stiffness(c, 2)).collect();
    {
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
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let stretch =
        if spec.materials.has_pml() { Some(femcore::ems_pml::cell_stretches(mesh, &spec.materials.pml)) } else { None };
    let k = match kind {
        CrossSection::Axisymmetric(m) => {
            crate::axisym::assemble_stiffness_axisym(mesh, &c6, order, m, stretch.as_deref())?
        }
    };
    let mut rho = spec.materials.cell_values(mesh, "rho", None)?;
    for (c, r) in rho.iter_mut().enumerate() {
        if spec.inactive.contains(&mesh.cell_labels[c]) {
            *r = 0.0;
        }
    }
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
    let free = free_dofs(n, &fixed);
    if n_modes < 1 || n_modes + 1 >= free.len() {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }
    let kf = k.submatrix(&free, &free);
    let mf = mc.submatrix(&free, &free);
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
    let cf = if lrb_labels.is_empty() {
        None
    } else {
        Some(crate::lrb::damping_matrix(mesh, &c6, &rho, &lrb_labels, order, 3).submatrix(&free, &free))
    };
    let has_loss = stretch.is_some() || !lrb_labels.is_empty();
    Ok(CrossAssembly { kind, order, n, n1, kf, mf, cf, free, fixed, c6, rho, n_modes, f_search, has_loss })
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
    let CrossSection::Axisymmetric(m_az) = asm.kind;
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
        let mut u = vec![C::new(0.0, 0.0); n];
        for (a, &dof) in free.iter().enumerate() {
            u[dof] = vectors[i][a];
        }
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

/// Stage 3: the plot fields `|u|`, `ux`, `uz`, `uy`, `|uy|` of every mode (largest component real, unit maximum), the
/// displacement of the (x, z) blocks, and the un-averaged per-cell complex components `ux`, `uz`, `uy`.
pub fn cross_section_fields(
    mesh: &Discretisation,
    asm: &CrossAssembly,
    solved: &CrossSolve,
) -> Result<SolveResult, String> {
    let (n1, order, has_loss) = (asm.n1, asm.order, asm.has_loss);
    let modes = &solved.modes;
    let mut solutions = vec![];
    let mut scalars = vec![];
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
        let fields = vec![
            Field::scalar("|u|", "a.u.", mag),
            Field::scalar(names[0], "a.u.", re(&ux)),
            Field::scalar(names[1], "a.u.", re(&uy)),
            Field::scalar(names[2], "a.u.", re(&uz)),
            Field::scalar(format!("|{}|", names[2]), "a.u.", uz.iter().map(|z| z.norm()).collect()),
        ];
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
                ..Default::default()
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", fz.re, fz.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        } else {
            solutions.push(Solution {
                label: format!("mode {}: {:.6} Hz", j + 1, freq),
                fields,
                cell_fields: cell_fields_of(&cells),
                value: Some(freq),
                ..Default::default()
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(freq)));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let m_az = modes.m_az;
    let what = match asm.kind {
        CrossSection::Axisymmetric(_) => format!("axisymmetric, azimuthal order m = {m_az}"),
    };
    let mut log = vec![format!(
        "{what}; complex Hermitian eigenproblem, Krylov dimension {}{}",
        solved.dim,
        if modes.converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    if asm.fixed.is_empty() {
        log.push("no fixed boundary: the lowest modes may be rigid-body-like".into());
    }
    Ok(SolveResult { order: order as u8, dim: 2, solutions, scalars, log, dofs: asm.n, ..Default::default() })
}
