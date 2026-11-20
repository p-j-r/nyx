// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The shared data of a problem that is not specific to one physics: material properties per domain (scalars, anisotropic tensors with a
//! crystal orientation, perfectly matched layers), Floquet-periodic boundary pairs and the choice of linear solver. The physics crates
//! (`nyx-dfs`, `nyx-sms`, `nyx-ems`) decode their problem into these; the old `problem` module re-exports them for the window until
//! phase 7. Property names: `eps_r`, `mu_r`, `rho_v`, `k`, `rho`, `Cp`, `Q`, `E`, `nu`, `alpha`.

use crate::material::{self, Anisotropy, M3, M6};
use crate::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Material properties per domain; `default` applies to domains without their own value.
/// Property names: `eps_r`, `mu_r`, `rho_v` (charge density), `k`, `rho`, `Cp`, `Q` (heat
/// source), `E`, `nu`, `alpha`.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Materials {
    pub domains: BTreeMap<i32, BTreeMap<String, f64>>,
    pub default: BTreeMap<String, f64>,
    /// Anisotropic data per domain (crystal stiffness, tensor conductivity, permittivity and
    /// expansion, orientation). Where present it replaces the scalar properties it covers.
    #[serde(default)]
    pub aniso: BTreeMap<i32, Anisotropy>,
    /// Perfectly matched layers per domain (waveguide mode analysis).
    #[serde(default)]
    pub pml: BTreeMap<i32, Pml>,
}

/// Which side of a domain a perfectly matched layer absorbs towards, per axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PmlDir {
    Off,
    /// the layer starts at the domain's minimum coordinate and ends at its maximum
    Plus,
    /// the layer starts at the maximum coordinate and ends at the minimum
    Minus,
}

/// A perfectly matched layer by coordinate stretching, `s = re - j strength t^order` with `t`
/// running from 0 at the interface to 1 at the outer edge of the domain along each active axis.
/// With `order` 0 the stretch is the constant `re - j strength`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Pml {
    pub x: PmlDir,
    /// 3D only: y is out of the plane of a 2D cross-section (x-z), where only `x` and `z` exist
    pub y: PmlDir,
    #[serde(default = "pml_off")]
    pub z: PmlDir,
    pub strength: f64,
    pub order: f64,
    pub re: f64,
}

fn pml_off() -> PmlDir {
    PmlDir::Off
}

impl Default for Pml {
    fn default() -> Self {
        Pml { x: PmlDir::Off, y: PmlDir::Off, z: PmlDir::Off, strength: 3.0, order: 2.0, re: 1.0 }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TensorKind {
    Conductivity,
    Permittivity,
    Expansion,
}

impl TensorKind {
    fn scalar_name(self) -> &'static str {
        match self {
            TensorKind::Conductivity => "k",
            TensorKind::Permittivity => "eps_r",
            TensorKind::Expansion => "alpha",
        }
    }
    fn pick(self, a: &Anisotropy) -> Option<&material::Tensor> {
        match self {
            TensorKind::Conductivity => a.conductivity.as_ref(),
            TensorKind::Permittivity => a.permittivity.as_ref(),
            TensorKind::Expansion => a.expansion.as_ref(),
        }
    }
}

impl Materials {
    /// True if any domain has an active perfectly matched layer.
    pub fn has_pml(&self) -> bool {
        self.pml.values().any(|p| p.x != PmlDir::Off || p.y != PmlDir::Off || p.z != PmlDir::Off)
    }

    pub fn get(&self, domain: i32, name: &str) -> Option<f64> {
        self.domains.get(&domain).and_then(|m| m.get(name)).or_else(|| self.default.get(name)).copied()
    }
    pub fn set(&mut self, domain: i32, name: &str, v: f64) {
        self.domains.entry(domain).or_default().insert(name.to_string(), v);
    }
    /// True if any domain defines the tensor `kind` anisotropically.
    pub fn has_tensor(&self, kind: TensorKind) -> bool {
        self.aniso.values().any(|a| kind.pick(a).is_some())
    }

    pub fn has_stiffness(&self) -> bool {
        self.aniso.values().any(|a| a.stiffness.is_some())
    }

    /// Stiffness matrix (solver Voigt order) of every cell: crystal data where defined, else
    /// isotropic from E and nu.
    pub fn cell_stiffness(&self, mesh: &SimplexMesh) -> Result<Vec<M6>, String> {
        self.cell_stiffness_skipping(mesh, &[])
    }

    /// Stiffness of every cell; domains in `skip` (not part of this physics) get a zero matrix.
    pub fn cell_stiffness_skipping(&self, mesh: &SimplexMesh, skip: &[i32]) -> Result<Vec<M6>, String> {
        let mut by_domain: BTreeMap<i32, M6> = BTreeMap::new();
        for d in mesh.domain_labels() {
            if skip.contains(&d) {
                by_domain.insert(d, [[0.0; 6]; 6]);
                continue;
            }
            let c = if let Some(st) = self.aniso.get(&d).and_then(|a| a.stiffness.as_ref()) {
                material::stiffness_global(st, self.aniso[&d].euler_deg)
            } else {
                let e = self.get(d, "E").ok_or_else(|| format!("domain {d} has no value for 'E'"))?;
                let nu = self.get(d, "nu").ok_or_else(|| format!("domain {d} has no value for 'nu'"))?;
                if nu >= 0.5 || nu <= -1.0 {
                    return Err("Poisson's ratio must be between -1 and 0.5".into());
                }
                if e <= 0.0 {
                    return Err("Young's modulus must be positive".into());
                }
                material::isotropic_stiffness(e, nu)
            };
            if !material::is_positive_definite(&c) {
                return Err(format!("the stiffness of domain {d} is not positive definite (check the constants)"));
            }
            by_domain.insert(d, c);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
    }

    /// Symmetric tensor property of every cell (global frame): tensor data where defined,
    /// else the scalar property times the identity.
    pub fn cell_tensor(&self, mesh: &SimplexMesh, kind: TensorKind, fallback: Option<f64>) -> Result<Vec<M3>, String> {
        let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
        for d in mesh.domain_labels() {
            let t = match self.aniso.get(&d).and_then(|a| kind.pick(a)) {
                Some(t) => material::tensor_global(t, self.aniso[&d].euler_deg),
                None => {
                    let name = kind.scalar_name();
                    let v = self.get(d, name).or(fallback).ok_or_else(|| format!("domain {d} has no value for '{name}'"))?;
                    material::iso_tensor(v)
                }
            };
            by_domain.insert(d, t);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
    }

    /// One value of `name` per cell. Without a `fallback`, a domain lacking the property is an error.
    pub fn cell_values(&self, mesh: &SimplexMesh, name: &str, fallback: Option<f64>) -> Result<Vec<f64>, String> {
        let mut by_domain: BTreeMap<i32, f64> = BTreeMap::new();
        for d in mesh.domain_labels() {
            let v = self
                .get(d, name)
                .or(fallback)
                .ok_or_else(|| format!("domain {d} has no value for '{name}'"))?;
            by_domain.insert(d, v);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
    }
}

/// A periodic condition with a Floquet phase between two boundaries of equal meshes: the dofs on the
/// `dest` boundaries are the dofs on the `source` boundaries shifted by the vector between the two
/// groups, times `exp(-i phase)` (the phase advance `k . d` over the period).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Periodic {
    pub source: Vec<i32>,
    pub dest: Vec<i32>,
    /// phase advance in radians over the period
    pub phase: f64,
}

/// How the large symmetric positive definite systems of static problems are solved.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum LinearSolver {
    /// sparse direct (Cholesky) for small systems, algebraic multigrid for large ones (3D from about 30 000 free dofs, 2D from 400 000)
    #[default]
    Auto,
    Direct,
    /// smoothed-aggregation algebraic multigrid with conjugate gradients
    Iterative,
}

impl LinearSolver {
    pub fn name(self) -> &'static str {
        match self {
            LinearSolver::Auto => "Automatic",
            LinearSolver::Direct => "Direct (Cholesky)",
            LinearSolver::Iterative => "Iterative (AMG-CG)",
        }
    }
}

