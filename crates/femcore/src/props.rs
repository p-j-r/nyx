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
    Permeability,
    PermeabilityIm,
    Expansion,
}

impl TensorKind {
    fn scalar_name(self) -> &'static str {
        match self {
            TensorKind::Conductivity => "k",
            TensorKind::Permittivity => "eps_r",
            TensorKind::PermittivityIm => "eps_r_im",
            TensorKind::Permeability => "mu_r",
            TensorKind::PermeabilityIm => "mu_r_im",
            TensorKind::Expansion => "alpha",
        }
    }
    fn pick(self, a: &Anisotropy) -> Option<&material::Tensor> {
        match self {
            TensorKind::Conductivity => a.conductivity.as_ref(),
            TensorKind::Permittivity => a.permittivity.as_ref(),
            TensorKind::PermittivityIm => a.permittivity_im.as_ref(),
            TensorKind::Permeability => a.permeability.as_ref(),
            TensorKind::PermeabilityIm => a.permeability_im.as_ref(),
            TensorKind::Expansion => a.expansion.as_ref(),
        }
    }
}

/// Exact `3x3` complex matrix inverse.
pub fn inv3_c(a: &[[crate::linalg::csr::C; 3]; 3]) -> Result<[[crate::linalg::csr::C; 3]; 3], String> {
    use crate::linalg::csr::c;
    let c00 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    let c01 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    let c02 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    let det = a[0][0] * c00 + a[0][1] * c01 + a[0][2] * c02;
    let scale = a.iter().flatten().fold(0.0_f64, |acc, z| acc.max(z.norm())).powi(3);
    if !det.re.is_finite() || !det.im.is_finite() || det.norm() <= 1e-14 * scale.max(1e-300) {
        return Err("3x3 complex permeability/permittivity tensor is singular and cannot be inverted".into());
    }
    let inv = c(1.0, 0.0) / det;
    Ok([
        [
            c00 * inv,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv,
        ],
        [
            c01 * inv,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv,
        ],
        [
            c02 * inv,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv,
        ],
    ])
}

#[allow(clippy::type_complexity)]
fn decompose_c3_into_sym_and_gyro(
    mat: &[[crate::linalg::csr::C; 3]; 3],
) -> (material::Tensor, Option<material::Tensor>, Option<[f64; 3]>, Option<[f64; 3]>) {
    let s_re = [
        [mat[0][0].re, 0.5 * (mat[0][1].re + mat[1][0].re), 0.5 * (mat[0][2].re + mat[2][0].re)],
        [0.5 * (mat[0][1].re + mat[1][0].re), mat[1][1].re, 0.5 * (mat[1][2].re + mat[2][1].re)],
        [0.5 * (mat[0][2].re + mat[2][0].re), 0.5 * (mat[1][2].re + mat[2][1].re), mat[2][2].re],
    ];
    let s_im = [
        [-mat[0][0].im, -0.5 * (mat[0][1].im + mat[1][0].im), -0.5 * (mat[0][2].im + mat[2][0].im)],
        [-0.5 * (mat[0][1].im + mat[1][0].im), -mat[1][1].im, -0.5 * (mat[1][2].im + mat[2][1].im)],
        [-0.5 * (mat[0][2].im + mat[2][0].im), -0.5 * (mat[1][2].im + mat[2][1].im), -mat[2][2].im],
    ];
    let g_re = [
        0.5 * (mat[2][1].im - mat[1][2].im),
        0.5 * (mat[0][2].im - mat[2][0].im),
        0.5 * (mat[1][0].im - mat[0][1].im),
    ];
    let g_im = [
        0.5 * (mat[1][2].re - mat[2][1].re),
        0.5 * (mat[2][0].re - mat[0][2].re),
        0.5 * (mat[0][1].re - mat[1][0].re),
    ];
    let to_sym = |m: &M3| material::Tensor::Symmetric {
        xx: m[0][0],
        yy: m[1][1],
        zz: m[2][2],
        xy: m[0][1],
        yz: m[1][2],
        xz: m[0][2],
    };
    let t_re = to_sym(&s_re);
    let t_im = s_im.iter().flatten().any(|&v| v != 0.0).then(|| to_sym(&s_im));
    let gy_re = g_re.iter().any(|&v| v != 0.0).then_some(g_re);
    let gy_im = g_im.iter().any(|&v| v != 0.0).then_some(g_im);
    (t_re, t_im, gy_re, gy_im)
}

impl Materials {
    /// True if any domain has an active perfectly matched layer.
    pub fn has_pml(&self) -> bool {
        self.pml.values().any(Pml::is_active)
    }

    /// True if any domain has a non-symmetric (gyrotropic / magneto-optic) permittivity or permeability tensor.
    pub fn is_non_symmetric(&self) -> bool {
        self.aniso.values().any(Anisotropy::is_non_symmetric)
    }

    /// True if any domain carries a non-zero imaginary permittivity (dielectric/optical loss or gyrotropy).
    pub fn has_complex_permittivity(&self) -> bool {
        self.default.get("eps_r_im").is_some_and(|v| v.abs() > 0.0)
            || self.domains.values().any(|m| m.get("eps_r_im").is_some_and(|v| v.abs() > 0.0))
            || self.aniso.values().any(|a| a.permittivity_im.is_some() || a.is_non_symmetric())
    }

    /// Sets the real electric gyrotropy / Voigt axial vector $\mathbf{g} = (g_x, g_y, g_z)$ on `domain`
    /// ($\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}_{\text{sym}} + i \mathbf{G}_\times$, where $\mathbf{G}_\times \mathbf{E} = \mathbf{g} \times \mathbf{E}$).
    pub fn set_gyrotropy(&mut self, domain: i32, g: [f64; 3]) {
        self.aniso.entry(domain).or_default().gyrotropy = Some(g);
    }

    /// Sets the complex electric gyrotropy axial vector $\mathbf{g} = \mathbf{g}_{\text{re}} + i \mathbf{g}_{\text{im}}$ on `domain`.
    pub fn set_gyrotropy_complex(&mut self, domain: i32, g_re: [f64; 3], g_im: [f64; 3]) {
        let an = self.aniso.entry(domain).or_default();
        an.gyrotropy = Some(g_re);
        an.gyrotropy_im = Some(g_im);
    }

    /// Sets the real magnetic gyrotropy / Polder axial vector $\mathbf{g}_\mu = (g_{\mu,x}, g_{\mu,y}, g_{\mu,z})$ on `domain`
    /// ($\boldsymbol{\mu} = \boldsymbol{\mu}_{\text{sym}} + i \mathbf{G}_{\mu,\times}$).
    pub fn set_mu_gyrotropy(&mut self, domain: i32, g: [f64; 3]) {
        self.aniso.entry(domain).or_default().mu_gyrotropy = Some(g);
    }

    /// Sets the complex magnetic gyrotropy axial vector $\mathbf{g}_\mu = \mathbf{g}_{\mu,\text{re}} + i \mathbf{g}_{\mu,\text{im}}$ on `domain`.
    pub fn set_mu_gyrotropy_complex(&mut self, domain: i32, g_re: [f64; 3], g_im: [f64; 3]) {
        let an = self.aniso.entry(domain).or_default();
        an.mu_gyrotropy = Some(g_re);
        an.mu_gyrotropy_im = Some(g_im);
    }

    /// Sets an arbitrary 9-entry complex $3 \times 3$ permittivity tensor $\boldsymbol{\varepsilon} \in \mathbb{C}^{3\times 3}$ on `domain`
    /// (in global coordinates, resetting `euler_deg = [0.0; 3]`).
    pub fn set_permittivity_tensor_c(&mut self, domain: i32, eps: [[crate::linalg::csr::C; 3]; 3]) {
        let (t_re, t_im, g_re, g_im) = decompose_c3_into_sym_and_gyro(&eps);
        let an = self.aniso.entry(domain).or_default();
        an.euler_deg = [0.0; 3];
        an.permittivity = Some(t_re);
        an.permittivity_im = t_im;
        an.gyrotropy = g_re;
        an.gyrotropy_im = g_im;
    }

    /// Sets an arbitrary 9-entry complex $3 \times 3$ permeability tensor $\boldsymbol{\mu} \in \mathbb{C}^{3\times 3}$ on `domain`
    /// (in global coordinates, resetting `euler_deg = [0.0; 3]`).
    pub fn set_permeability_tensor_c(&mut self, domain: i32, mu: [[crate::linalg::csr::C; 3]; 3]) {
        let (t_re, t_im, g_re, g_im) = decompose_c3_into_sym_and_gyro(&mu);
        let an = self.aniso.entry(domain).or_default();
        an.euler_deg = [0.0; 3];
        an.permeability = Some(t_re);
        an.permeability_im = t_im;
        an.mu_gyrotropy = g_re;
        an.mu_gyrotropy_im = g_im;
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

    /// True if any domain carries a piezoelectric coupling tensor `e_iJ`.
    pub fn has_piezoelectric(&self) -> bool {
        self.aniso.values().any(|a| a.piezo_e.is_some())
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

    /// Global-frame `3x6` piezoelectric stress-charge matrix `e` (solver Voigt order, `C/m^2`) of every cell.
    pub fn cell_piezo_e(&self, mesh: &SimplexMesh) -> Vec<[[f64; 6]; 3]> {
        self.cell_piezo_e_skipping(mesh, &[])
    }

    /// Global-frame `3x6` piezoelectric stress-charge matrix `e` (solver Voigt order, `C/m^2`) of every cell;
    /// domains in `skip` or without `[piezoelectric]` get a zero matrix.
    pub fn cell_piezo_e_skipping(&self, mesh: &SimplexMesh, skip: &[i32]) -> Vec<[[f64; 6]; 3]> {
        let mut by_domain: BTreeMap<i32, [[f64; 6]; 3]> = BTreeMap::new();
        for d in mesh.domain_labels() {
            let e = if skip.contains(&d) {
                [[0.0; 6]; 3]
            } else if let Some(e_tb) = self.aniso.get(&d).and_then(|a| a.piezo_e.as_ref()) {
                material::piezo_e_global(e_tb, self.aniso[&d].euler_deg)
            } else {
                [[0.0; 6]; 3]
            };
            by_domain.insert(d, e);
        }
        mesh.cell_labels.iter().map(|l| by_domain[l]).collect()
    }

    /// Symmetric tensor property of every cell (global frame): tensor data where defined,
    /// else the scalar property times the identity.
    pub fn cell_tensor(&self, mesh: &SimplexMesh, kind: TensorKind, fallback: Option<f64>) -> Result<Vec<M3>, String> {
        self.cell_tensor_skipping(mesh, kind, fallback, &[])
    }

    /// Symmetric tensor property of every cell (global frame); domains in `skip` get a zero matrix
    /// without requiring the scalar or tensor property to exist on that domain.
    pub fn cell_tensor_skipping(
        &self,
        mesh: &SimplexMesh,
        kind: TensorKind,
        fallback: Option<f64>,
        skip: &[i32],
    ) -> Result<Vec<M3>, String> {
        let mut by_domain: BTreeMap<i32, M3> = BTreeMap::new();
        for d in mesh.domain_labels() {
            if skip.contains(&d) {
                by_domain.insert(d, [[0.0; 3]; 3]);
                continue;
            }
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
        self.cell_complex_permittivity_skipping(mesh, fallback_re, &[])
    }

    /// Global-frame `(eps_re, eps_im_loss)` permittivity pair of every cell; domains in `skip` get zero tensors.
    pub fn cell_complex_permittivity_skipping(
        &self,
        mesh: &SimplexMesh,
        fallback_re: Option<f64>,
        skip: &[i32],
    ) -> Result<Vec<(M3, M3)>, String> {
        let re = self.cell_tensor_skipping(mesh, TensorKind::Permittivity, fallback_re, skip)?;
        let im = self.cell_tensor_skipping(mesh, TensorKind::PermittivityIm, Some(0.0), skip)?;
        Ok(re.into_iter().zip(im).collect())
    }

    /// Evaluates the full 9-entry complex relative permittivity tensor $\boldsymbol{\varepsilon} = \boldsymbol{\varepsilon}_{\text{sym}} + i \mathbf{G}_\times$
    /// in the global frame for `domain`.
    pub fn domain_full_complex_permittivity(
        &self,
        domain: i32,
        fallback_re: Option<f64>,
    ) -> Result<[[crate::linalg::csr::C; 3]; 3], String> {
        use crate::linalg::csr::c;
        let an_opt = self.aniso.get(&domain);
        let euler = an_opt.map_or([0.0; 3], |a| a.euler_deg);
        let t_re = match an_opt.and_then(|a| a.permittivity.as_ref()) {
            Some(t) => material::tensor_global(t, euler),
            None => {
                let v = self
                    .get(domain, "eps_r")
                    .or(fallback_re)
                    .ok_or_else(|| format!("domain {domain} has no value for 'eps_r'"))?;
                material::iso_tensor(v)
            }
        };
        let t_im = match an_opt.and_then(|a| a.permittivity_im.as_ref()) {
            Some(t) => material::tensor_global(t, euler),
            None => {
                let v = self.get(domain, "eps_r_im").unwrap_or(0.0);
                material::iso_tensor(v)
            }
        };
        let mut out = [[c(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = c(t_re[i][j], -t_im[i][j]);
            }
        }
        if let Some((g_re, g_im)) = an_opt.and_then(Anisotropy::gyrotropy_global) {
            // +i * G_x where (G_x E) = g x E and i * (g_re + i g_im) = -g_im + i g_re
            let ig_x = c(-g_im[0], g_re[0]);
            let ig_y = c(-g_im[1], g_re[1]);
            let ig_z = c(-g_im[2], g_re[2]);
            out[0][1] -= ig_z;
            out[1][0] += ig_z;
            out[0][2] += ig_y;
            out[2][0] -= ig_y;
            out[1][2] -= ig_x;
            out[2][1] += ig_x;
        }
        Ok(out)
    }

    /// Evaluates the full 9-entry complex relative permeability tensor $\boldsymbol{\mu} = \boldsymbol{\mu}_{\text{sym}} + i \mathbf{G}_{\mu,\times}$
    /// in the global frame for `domain`.
    pub fn domain_full_complex_permeability(
        &self,
        domain: i32,
        fallback_re: Option<f64>,
    ) -> Result<[[crate::linalg::csr::C; 3]; 3], String> {
        use crate::linalg::csr::c;
        let an_opt = self.aniso.get(&domain);
        let euler = an_opt.map_or([0.0; 3], |a| a.euler_deg);
        let t_re = match an_opt.and_then(|a| a.permeability.as_ref()) {
            Some(t) => material::tensor_global(t, euler),
            None => {
                let v = self
                    .get(domain, "mu_r")
                    .or(fallback_re)
                    .ok_or_else(|| format!("domain {domain} has no value for 'mu_r'"))?;
                material::iso_tensor(v)
            }
        };
        let t_im = match an_opt.and_then(|a| a.permeability_im.as_ref()) {
            Some(t) => material::tensor_global(t, euler),
            None => {
                let v = self.get(domain, "mu_r_im").unwrap_or(0.0);
                material::iso_tensor(v)
            }
        };
        let mut out = [[c(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = c(t_re[i][j], -t_im[i][j]);
            }
        }
        if let Some((g_re, g_im)) = an_opt.and_then(Anisotropy::mu_gyrotropy_global) {
            let ig_x = c(-g_im[0], g_re[0]);
            let ig_y = c(-g_im[1], g_re[1]);
            let ig_z = c(-g_im[2], g_re[2]);
            out[0][1] -= ig_z;
            out[1][0] += ig_z;
            out[0][2] += ig_y;
            out[2][0] -= ig_y;
            out[1][2] -= ig_x;
            out[2][1] += ig_x;
        }
        Ok(out)
    }

    /// Global-frame full 9-entry complex permittivity tensor $\boldsymbol{\varepsilon} \in \mathbb{C}^{3\times 3}$ of every cell.
    pub fn cell_full_complex_permittivity(
        &self,
        mesh: &SimplexMesh,
        fallback_re: Option<f64>,
    ) -> Result<Vec<[[crate::linalg::csr::C; 3]; 3]>, String> {
        let mut by_domain = BTreeMap::new();
        for d in mesh.domain_labels() {
            by_domain.insert(d, self.domain_full_complex_permittivity(d, fallback_re)?);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
    }

    /// Global-frame full 9-entry complex permeability tensor $\boldsymbol{\mu} \in \mathbb{C}^{3\times 3}$ of every cell.
    pub fn cell_full_complex_permeability(
        &self,
        mesh: &SimplexMesh,
        fallback_re: Option<f64>,
    ) -> Result<Vec<[[crate::linalg::csr::C; 3]; 3]>, String> {
        let mut by_domain = BTreeMap::new();
        for d in mesh.domain_labels() {
            by_domain.insert(d, self.domain_full_complex_permeability(d, fallback_re)?);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
    }

    /// Global-frame inverse complex permeability tensor $\boldsymbol{\mu}^{-1} \in \mathbb{C}^{3\times 3}$ of every cell.
    pub fn cell_full_complex_permeability_inv(
        &self,
        mesh: &SimplexMesh,
        fallback_re: Option<f64>,
    ) -> Result<Vec<[[crate::linalg::csr::C; 3]; 3]>, String> {
        let mut by_domain = BTreeMap::new();
        for d in mesh.domain_labels() {
            let mu = self.domain_full_complex_permeability(d, fallback_re)?;
            let inv = inv3_c(&mu).map_err(|e| format!("domain {d}: {e}"))?;
            by_domain.insert(d, inv);
        }
        Ok(mesh.cell_labels.iter().map(|l| by_domain[l]).collect())
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

    #[test]
    fn skipping_helpers_bypass_excluded_domains() {
        let mut mesh = crate::testkit::rect(1.0, 1.0, 2, 2);
        // Mark right-half cells as domain 2 (e.g., conductor without eps_r)
        for ci in 0..mesh.n_cells() {
            let v = mesh.cell(ci);
            let cx = (0..3).map(|k| mesh.point(v[k] as usize)[0]).sum::<f64>() / 3.0;
            mesh.cell_labels[ci] = if cx > 0.5 { 2 } else { 1 };
        }
        let mut m = Materials::default();
        m.set(1, "eps_r", 9.5);
        m.aniso.insert(
            1,
            material::Anisotropy {
                piezo_e: Some([[0.0, 0.0, 0.0, 0.0, -0.48, 0.0], [0.0, 0.0, 0.0, -0.48, 0.0, 0.0], [-0.58, -0.58, 1.55, 0.0, 0.0, 0.0]]),
                ..Default::default()
            },
        );
        assert!(m.has_piezoelectric());
        // Domain 2 has no eps_r: without skipping it fails, with skipping [2] it succeeds
        assert!(m.cell_complex_permittivity(&mesh, None).is_err());
        let eps = m.cell_complex_permittivity_skipping(&mesh, None, &[2]).unwrap();
        let pz = m.cell_piezo_e_skipping(&mesh, &[2]);
        for ci in 0..mesh.n_cells() {
            if mesh.cell_labels[ci] == 2 {
                assert_eq!(eps[ci].0, [[0.0; 3]; 3]);
                assert_eq!(pz[ci], [[0.0; 6]; 3]);
            } else {
                assert!((eps[ci].0[0][0] - 9.5).abs() < 1e-14);
                assert!((pz[ci][2][2] - 1.55).abs() < 1e-14);
            }
        }
    }
}


