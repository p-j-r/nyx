// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The mode solvers' results as an `OpticalSet` (`femcore::optical`): the modes together with the material data needed to evaluate their
//! fields, for tests of the fields (the polarization analysis reads the cell fields of the result, `nyx-analysis`, ADR 56).
//! `solve_optical` solves a `Modes` problem and returns the set (it was in `femcore::optical` before the kernels moved here, ADR 45).

use crate::modes::{solve_modes_tensor, EpsCell, Modes};
use crate::modes_complex::{solve_modes_complex, ModesC};
use crate::problem::{Problem, Study};
use femcore::linalg::csr::c;
use femcore::optical::{cell_materials_signed, OptMode, OpticalSet};
use femcore::props::{Pml, TensorKind};
use nyx_physics::Discretisation;

fn from_real(m: &Modes) -> Vec<OptMode> {
    (0..m.neff.len())
        .map(|j| OptMode {
            neff: c(m.neff[j].0, m.neff[j].1),
            et: m.et[j].iter().map(|&x| c(x, 0.0)).collect(),
            phi: m.phi[j].iter().map(|&x| c(x, 0.0)).collect(),
        })
        .collect()
}

fn from_complex(m: &ModesC) -> Vec<OptMode> {
    (0..m.neff.len()).map(|j| OptMode { neff: m.neff[j], et: m.et[j].clone(), phi: m.phi[j].clone() }).collect()
}

fn in_pml(mesh: &Discretisation, problem: &Problem) -> Vec<bool> {
    (0..mesh.n_cells()).map(|ci| problem.materials.pml.get(&mesh.cell_labels[ci]).is_some_and(Pml::is_active)).collect()
}

/// Solves the optical modes of `problem` (a `Modes` study in wafer axes), propagating along `+y`; complex materials and PMLs go through the
/// complex solver.
pub fn solve_optical(mesh: &Discretisation, problem: &Problem) -> Result<OpticalSet, String> {
    solve_optical_wafer(mesh, problem, femcore::section::PROP_SIGN)
}

/// Like `solve_optical` for a propagation axis along `prop_sign * y` (only the reciprocity test uses a sign other than `PROP_SIGN`).
pub fn solve_optical_wafer(mesh: &Discretisation, problem: &Problem, prop_sign: f64) -> Result<OpticalSet, String> {
    if mesh.has_curved() {
        return Err(
            "curved meshes are not supported by simplex optical mode kernels; use mixed/isoparametric path or a straight mesh"
                .into(),
        );
    }
    let Study::Modes { freq, n_target, evanescent, n_modes } = problem.study.clone() else {
        return Err("the optical problem needs a Modes study".into());
    };
    if freq <= 0.0 || n_modes < 1 {
        return Err("frequency must be positive and at least one mode is needed".into());
    }
    let (eps, mu) = cell_materials_signed(mesh, &problem.materials, prop_sign)?;
    let in_pml = in_pml(mesh, problem);
    let coupled = eps.iter().any(|e| e.xz.norm() > 0.0 || e.yz.norm() > 0.0);
    let (modes, converged) = if problem.materials.has_pml() || coupled {
        let target = if evanescent { c(0.0, n_target) } else { c(n_target, 0.0) };
        let m = solve_modes_complex(mesh, freq, &eps, &mu, target, n_modes, problem.order)?;
        (from_complex(&m), m.converged)
    } else {
        let tens: Vec<_> = problem
            .materials
            .cell_tensor(mesh, TensorKind::Permittivity, Some(1.0))?
            .iter()
            .map(|t| femcore::section::local_tensor_signed(t, prop_sign))
            .collect();
        let e: Vec<EpsCell> =
            tens.iter().map(|t| EpsCell { xx: t[0][0], xy: t[0][1], yy: t[1][1], zz: t[2][2] }).collect();
        let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
        let m = solve_modes_tensor(mesh, freq, &e, &mur, n_target, evanescent, n_modes, problem.order)?;
        (from_real(&m), m.converged)
    };
    Ok(OpticalSet {
        freq,
        order: problem.order,
        eps,
        mu,
        modes,
        converged,
        in_pml,
        left_handed: femcore::section::is_left_handed(prop_sign),
    })
}

/// Domains of the problem that are absorbing layers: the polarization analysis leaves their cells out.
pub(crate) fn absorbing_domains(problem: &Problem) -> Vec<i32> {
    problem.materials.pml.iter().filter(|(_, p)| p.is_active()).map(|(&d, _)| d).collect()
}

/// The frame of the modes of a cross-section study: lateral x, normal z, propagation y (the local triple of the mode solvers, left-handed),
/// without the absorbing layers. Written to every mode result as the scalar `field frame`.
pub(crate) fn mode_frame(problem: &Problem) -> nyx_physics::FieldFrame {
    nyx_physics::FieldFrame {
        lateral: 0,
        normal: 2,
        propagation: 1,
        left_handed: femcore::section::is_left_handed(femcore::section::PROP_SIGN),
        outside: absorbing_domains(problem),
    }
}
