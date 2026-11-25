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

/// Coordinate system and reference geometry of a perfectly matched layer (or infinite element).
/// All 3D position and axis vectors (`center`, `axis`) are expressed in the global wafer frame `(x, y, z)`,
/// where a 2D mesh lies in the `(x, z)` cross-section (`y = 0`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum PmlCoord {
    /// Axis-aligned Cartesian slab, edge, or corner (`Pml::x`, `Pml::y`, `Pml::z`).
    #[default]
    Cartesian,
    /// Cylindrical shell (`radial_on`: `r_in .. r_in + r_thickness` perpendicular to unit `axis`),
    /// axial endcap (`axial`: `Plus` runs `z_in .. z_in + z_thickness`, `Minus` runs `z_in - z_thickness .. z_in`
    /// along `axis` measured from `center`), or corner annulus (`radial_on && axial != Off`).
    Cylindrical {
        center: [f64; 3],
        axis: [f64; 3],
        r_in: f64,
        r_thickness: f64,
        radial_on: bool,
        axial: PmlDir,
        z_in: f64,
        z_thickness: f64,
    },
    /// Concentric spherical shell (`r_in .. r_in + thickness`) around `center`.
    Spherical {
        center: [f64; 3],
        r_in: f64,
        thickness: f64,
    },
}

/// Radial/normal coordinate-stretching law inside a PML (`t` in `[0, 1]`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum PmlProfile {
    /// Bounded ABCS complex dilation `s(t) = 1 + (kappa_max - 1) t^m - i * strength * t^m`
    /// (with legacy constant `re` honored when `kappa_max == 1.0 && re != 1.0`).
    /// Default for eigenfrequency studies, waveguides, and iterative `AMS`/`GMRES` solves.
    KappaPolynomial { kappa_max: f64 },
    /// Locking-free rational stretching `s(t) = (1 - i * strength * t^m) / (1 - curvature * t)^2`
    /// (`0 <= curvature < 1`) for direct-solved open radiation/scattering problems.
    Rational { curvature: f64 },
}

impl Default for PmlProfile {
    fn default() -> Self {
        PmlProfile::KappaPolynomial { kappa_max: 1.0 }
    }
}

/// Stabilization mode against backward-propagating anisotropic elastic or guided Lamb waves.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum PmlStabilization {
    /// Pure analytic continuation along the normal ray only (`s_perp = 1`).
    Off,
    /// Pre-solve Christoffel check enables M-PML automatically when bulk backward waves exist.
    #[default]
    Auto,
    /// Multi-axial PML auxiliary tangential damping `s_parallel = 1 - i * ratio * psi(t)`.
    Mpml { ratio: f64 },
    /// Adiabatic stiffness loss ramp `C_tilde(t) = (1 + i * eta_max * t^2) C`.
    Adiabatic { eta_max: f64 },
}

fn pml_off() -> PmlDir {
    PmlDir::Off
}

fn default_target_reflection() -> f64 {
    1e-6
}

fn default_clamp_outer() -> bool {
    true
}

/// A perfectly matched layer by coordinate stretching. Along each active axis `t` runs from `0` at
/// the physical interface to `1` at the outer termination wall. Legacy fields (`x`, `y`, `z`,
/// `strength`, `order`, `re`) remain fully compatible; curvilinear coordinates, explicit slab
/// bounds, dual profiles, physics-based auto-scaling, and elastic stabilization extend it cleanly.
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
    /// Coordinate topology (`Cartesian`, `Cylindrical`, or `Spherical`).
    #[serde(default)]
    pub coord: PmlCoord,
    /// Explicit `(lo, hi)` bounding box in wafer coordinates `(x, y, z)` for Cartesian slabs so
    /// disconnected or shared PML domains never miscompute domain-local bounds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bounds: Option<([f64; 3], [f64; 3])>,
    /// Stretching law (`KappaPolynomial` or `Rational`).
    #[serde(default)]
    pub profile: PmlProfile,
    /// When `true`, `strength`, `profile`, and thickness/layer counts calibrate from wavelengths.
    #[serde(default)]
    pub auto_scale: bool,
    /// Target normal-incidence amplitude reflection `R_0` (default `1e-6`, i.e. `-120 dB` power).
    #[serde(default = "default_target_reflection")]
    pub target_reflection: f64,
    /// Backward-wave stabilization strategy for elastodynamics (`nyx-sms`).
    #[serde(default)]
    pub stabilization: PmlStabilization,
    /// Whether exterior PML termination walls are clamped (`u = 0` Dirichlet) in `nyx-sms`.
    #[serde(default = "default_clamp_outer")]
    pub clamp_outer: bool,
}

impl Default for Pml {
    fn default() -> Self {
        Pml {
            x: PmlDir::Off,
            y: PmlDir::Off,
            z: PmlDir::Off,
            strength: 3.0,
            order: 2.0,
            re: 1.0,
            coord: PmlCoord::Cartesian,
            bounds: None,
            profile: PmlProfile::default(),
            auto_scale: false,
            target_reflection: default_target_reflection(),
            stabilization: PmlStabilization::Auto,
            clamp_outer: default_clamp_outer(),
        }
    }
}

impl Pml {
    /// Creates an auto-scaled PML with the given coordinate topology (`auto_scale = true`).
    pub fn auto(coord: PmlCoord) -> Self {
        Pml { coord, auto_scale: true, ..Self::default() }
    }

    /// Creates an auto-scaled Cartesian PML along `(x, y, z)` in wafer axes (`auto_scale = true`).
    pub fn auto_cartesian(x: PmlDir, y: PmlDir, z: PmlDir) -> Self {
        Pml { x, y, z, auto_scale: true, ..Self::default() }
    }

    /// True if this PML actively stretches at least one coordinate direction.
    pub fn is_active(&self) -> bool {
        match self.coord {
            PmlCoord::Cartesian => self.x != PmlDir::Off || self.y != PmlDir::Off || self.z != PmlDir::Off,
            PmlCoord::Cylindrical { radial_on, axial, .. } => radial_on || axial != PmlDir::Off,
            PmlCoord::Spherical { thickness, .. } => thickness > 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TensorKind {
    Conductivity,
    Permittivity,
    PermittivityIm,
    Expansion,
}

impl TensorKind {
    fn scalar_name(self) -> &'static str {
        match self {
            TensorKind::Conductivity => "k",
            TensorKind::Permittivity => "eps_r",
            TensorKind::PermittivityIm => "eps_r_im",
            TensorKind::Expansion => "alpha",
        }
    }
    fn pick(self, a: &Anisotropy) -> Option<&material::Tensor> {
        match self {
            TensorKind::Conductivity => a.conductivity.as_ref(),
            TensorKind::Permittivity => a.permittivity.as_ref(),
            TensorKind::PermittivityIm => a.permittivity_im.as_ref(),
            TensorKind::Expansion => a.expansion.as_ref(),
        }
    }
}

impl Materials {
    /// True if any domain has an active perfectly matched layer.
    pub fn has_pml(&self) -> bool {
        self.pml.values().any(Pml::is_active)
    }

    /// True if any domain carries a non-zero imaginary permittivity (dielectric/optical loss).
    pub fn has_complex_permittivity(&self) -> bool {
        self.default.get("eps_r_im").is_some_and(|v| v.abs() > 0.0)
            || self.domains.values().any(|m| m.get("eps_r_im").is_some_and(|v| v.abs() > 0.0))
            || self.aniso.values().any(|a| a.permittivity_im.is_some())
    }

    /// True if any domain carries a non-zero phonon/structural viscosity tensor or scalar `eta`.
    pub fn has_viscosity(&self) -> bool {
        self.default.get("eta").is_some_and(|v| *v > 0.0)
            || self.domains.values().any(|m| m.get("eta").is_some_and(|v| *v > 0.0))
            || self.aniso.values().any(|a| a.viscosity.is_some())
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

    /// Global-frame `(eps_re, eps_im_loss)` permittivity pair of every cell.
    /// Under the `exp(+i omega t)` convention, passive dielectric loss corresponds to `eps = eps_re - i eps_im_loss`.
    pub fn cell_complex_permittivity(&self, mesh: &SimplexMesh, fallback_re: Option<f64>) -> Result<Vec<(M3, M3)>, String> {
        let re = self.cell_tensor(mesh, TensorKind::Permittivity, fallback_re)?;
        let im = self.cell_tensor(mesh, TensorKind::PermittivityIm, Some(0.0))?;
        Ok(re.into_iter().zip(im).collect())
    }

    /// Global-frame `6x6` dynamic viscosity tensor `eta` (solver Voigt order, `Pa*s`) of every cell.
    pub fn cell_viscosity(&self, mesh: &SimplexMesh) -> Result<Vec<M6>, String> {
        self.cell_viscosity_skipping(mesh, &[])
    }

    /// Global-frame `6x6` dynamic viscosity tensor `eta` (solver Voigt order, `Pa*s`) of every cell;
    /// domains in `skip` or without viscosity get a zero matrix.
    pub fn cell_viscosity_skipping(&self, mesh: &SimplexMesh, skip: &[i32]) -> Result<Vec<M6>, String> {
        let mut by_domain: BTreeMap<i32, M6> = BTreeMap::new();
        for d in mesh.domain_labels() {
            if skip.contains(&d) {
                by_domain.insert(d, [[0.0; 6]; 6]);
                continue;
            }
            let eta = if let Some(v) = self.aniso.get(&d).and_then(|a| a.viscosity.as_ref()) {
                material::rotate_m6_global(v, self.aniso[&d].euler_deg)
            } else if let Some(eta_s) = self.get(d, "eta") {
                if eta_s < 0.0 {
                    return Err(format!("viscosity 'eta' of domain {d} must be non-negative (found {eta_s})"));
                }
                if eta_s > 0.0 {
                    material::isotropic_stiffness(eta_s, 0.0)
                } else {
                    [[0.0; 6]; 6]
                }
            } else {
                [[0.0; 6]; 6]
            };
            by_domain.insert(d, eta);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_pml_json_deserializes_with_curvilinear_defaults() {
        let raw = r#"{"x":"Plus","y":"Off","z":"Minus","strength":4.5,"order":2.0,"re":1.0}"#;
        let p: Pml = serde_json::from_str(raw).unwrap();
        assert_eq!(p.x, PmlDir::Plus);
        assert_eq!(p.y, PmlDir::Off);
        assert_eq!(p.z, PmlDir::Minus);
        assert_eq!(p.coord, PmlCoord::Cartesian);
        assert_eq!(p.bounds, None);
        assert_eq!(p.profile, PmlProfile::KappaPolynomial { kappa_max: 1.0 });
        assert!(!p.auto_scale);
        assert_eq!(p.target_reflection, 1e-6);
        assert_eq!(p.stabilization, PmlStabilization::Auto);
        assert!(p.clamp_outer);
    }

    #[test]
    fn has_pml_recognizes_cartesian_cylindrical_and_spherical_layers() {
        let mut m = Materials::default();
        assert!(!m.has_pml());

        m.pml.insert(1, Pml::default());
        assert!(!m.has_pml(), "inactive default Pml must not trigger has_pml");

        m.pml.insert(1, Pml::auto(PmlCoord::Spherical { center: [0.0; 3], r_in: 1.0, thickness: 0.0 }));
        assert!(!m.has_pml(), "zero-thickness spherical Pml is inactive");

        m.pml.insert(1, Pml::auto(PmlCoord::Spherical { center: [0.0; 3], r_in: 1.0, thickness: 0.2 }));
        assert!(m.has_pml());
        assert!(m.pml[&1].auto_scale);

        m.pml.insert(
            1,
            Pml::auto(PmlCoord::Cylindrical {
                center: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                r_in: 1.0,
                r_thickness: 0.2,
                radial_on: false,
                axial: PmlDir::Plus,
                z_in: 0.5,
                z_thickness: 0.2,
            }),
        );
        assert!(m.has_pml(), "axial-only cylindrical endcap PML must be recognized as active");
    }
}


