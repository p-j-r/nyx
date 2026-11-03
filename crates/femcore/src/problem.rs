//! The description of a problem and of its result: what the window edits and saves, what the solvers take in and give back.
//!
//! A `ProblemSpec` names one physics (electrostatics, heat transfer, solid mechanics, EM waves),
//! the element order (1 or 2), the material data per domain (`Materials`: scalar properties, anisotropic tensors with a
//! crystal orientation, perfectly matched layers), the boundary conditions (`Bc`: boundary numbers and a `BcKind`),
//! the study (`Study`: stationary/static, transient, eigenfrequencies, waveguide modes) and the options (linear solver,
//! plane stress/strain and thickness in 2D, thermal expansion, out-of-plane wave number, inactive domains, Floquet-periodic
//! pairs, results file). Everything derives serde, so it is saved as JSON or, through `xmlio`, as XML.
//!
//! A problem is written in the axes of the wafer: x lateral, y the propagation direction, z the surface normal (right-handed; a 2D
//! cross-section is the x-z plane, with y into the page). A wave along `+y` has a positive out-of-plane wave number `k`. Every kernel reads the problem in these axes and
//! names its result fields in them (`section` is the one place with axis numbers).
//!
//! A solve returns a `SolveResult`: one or more `Solution`s with named fields (one value per dof, P2: vertices then edge
//! midpoints) and optional displacement, scalar results, a text log (the stage log of `progress`), the time and the
//! number of unknowns. `SolveState` carries the last heat-transfer result so that a later solid study can use it.

use crate::material::{self, Anisotropy, M3, M6};
use crate::mesh::SimplexMesh;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Physics {
    Electrostatics,
    Heat,
    Solid,
    Emw,
}

impl Physics {
    pub const ALL: [Physics; 4] = [Physics::Electrostatics, Physics::Heat, Physics::Solid, Physics::Emw];
    pub fn name(self) -> &'static str {
        match self {
            Physics::Electrostatics => "Electrostatics",
            Physics::Heat => "Heat transfer",
            Physics::Solid => "Solid mechanics",
            Physics::Emw => "EM waves (mode analysis)",
        }
    }
}

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

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum BcKind {
    Potential(f64),
    Ground,
    /// surface charge density, C/m^2
    SurfaceCharge(f64),
    Temperature(f64),
    /// inward heat flux, W/m^2
    Flux(f64),
    Convection { h: f64, text: f64 },
    /// per-component constraint (x, y, z) to the given displacement
    Fixed { fixed: [bool; 3], values: [f64; 3] },
    /// boundary load (force per area)
    Load([f64; 3]),
    /// absorbing traction for outgoing pressure and shear waves (eigenfrequency studies: complex frequencies, Q)
    LowReflecting,
}

impl BcKind {
    pub fn name(&self) -> &'static str {
        match self {
            BcKind::Potential(_) => "Electric potential",
            BcKind::Ground => "Ground",
            BcKind::SurfaceCharge(_) => "Surface charge",
            BcKind::Temperature(_) => "Temperature",
            BcKind::Flux(_) => "Heat flux",
            BcKind::Convection { .. } => "Convective cooling",
            BcKind::Fixed { .. } => "Fixed constraint",
            BcKind::Load(_) => "Boundary load",
            BcKind::LowReflecting => "Low-reflecting boundary",
        }
    }
    /// Boundary conditions that make sense for a physics.
    pub fn available(p: Physics) -> Vec<BcKind> {
        match p {
            Physics::Electrostatics => vec![BcKind::Potential(1.0), BcKind::Ground, BcKind::SurfaceCharge(0.0)],
            Physics::Heat => vec![BcKind::Temperature(293.15), BcKind::Flux(0.0), BcKind::Convection { h: 10.0, text: 293.15 }],
            Physics::Solid => vec![BcKind::Fixed { fixed: [true; 3], values: [0.0; 3] }, BcKind::Load([0.0; 3]), BcKind::LowReflecting],
            Physics::Emw => vec![],
        }
    }
    pub fn belongs_to(&self, p: Physics) -> bool {
        match self {
            BcKind::Potential(_) | BcKind::Ground | BcKind::SurfaceCharge(_) => p == Physics::Electrostatics,
            BcKind::Temperature(_) | BcKind::Flux(_) | BcKind::Convection { .. } => p == Physics::Heat,
            BcKind::Fixed { .. } | BcKind::Load(_) | BcKind::LowReflecting => p == Physics::Solid,
        }
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

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Bc {
    /// a selector over the geometry model (`outer & plane(x = 0)`); when set the labels are what it resolves to, refreshed whenever the model or mesh changes
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub on: Option<String>,
    pub labels: Vec<i32>,
    pub kind: BcKind,
}

/// Where a simulation writes its results besides showing them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum OutputMode {
    #[default]
    None,
    /// HDF5 through the standard sequential driver
    Hdf5Serial,
    /// parallel HDF5 (MPI-IO); needs a build with the `phdf5` feature
    Hdf5Parallel,
}

impl OutputMode {
    pub fn name(self) -> &'static str {
        match self {
            OutputMode::None => "No file",
            OutputMode::Hdf5Serial => "HDF5 (serial)",
            OutputMode::Hdf5Parallel => "HDF5 (parallel, MPI-IO)",
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct OutputSpec {
    pub mode: OutputMode,
    /// path of the results file (written after every Compute when the mode is not `None`)
    pub path: String,
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Plane {
    Strain,
    Stress,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ThermalSource {
    Uniform(f64),
    /// the last heat-transfer solution (final time for a transient run)
    FromHeat,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ThermalExpansion {
    pub source: ThermalSource,
    pub t_ref: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Study {
    Stationary,
    /// `t0` is an expression in x, y, z for the initial temperature
    Transient { t_end: f64, n_out: usize, t0: String },
    Static,
    /// eigenfrequencies; `f_search` > 0 finds modes near that frequency (Hz)
    Eigenfrequency { n_modes: usize, f_search: f64 },
    /// mode analysis at `freq` (Hz) near effective index `n_target` (imaginary if `evanescent`)
    Modes { freq: f64, n_target: f64, evanescent: bool, n_modes: usize },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ProblemSpec {
    pub physics: Physics,
    /// element order 1 or 2
    pub order: usize,
    pub materials: Materials,
    pub bcs: Vec<Bc>,
    pub study: Study,
    pub plane: Plane,
    /// 2D out-of-plane thickness (m)
    pub thickness: f64,
    pub thermal: Option<ThermalExpansion>,
    #[serde(default)]
    pub solver: LinearSolver,
    #[serde(default)]
    pub output: OutputSpec,
    /// solid mechanics on a 2D cross-section: out-of-plane wave number k (rad/m), positive for a wave along +y; the three
    /// displacement components live on the cross-section (eigenfrequencies only)
    #[serde(default)]
    pub out_of_plane: Option<f64>,
    /// domains that do not take part in this physics (for example air in a solid mechanics study):
    /// they carry no stiffness or mass and the dofs they alone touch are held at zero
    #[serde(default)]
    pub inactive: Vec<i32>,
    /// Floquet-periodic boundary pairs (eigenfrequency studies of solid mechanics and 3D EM waves)
    #[serde(default)]
    pub periodic: Vec<Periodic>,
}

impl ProblemSpec {
    pub fn new(physics: Physics) -> Self {
        let study = match physics {
            Physics::Electrostatics => Study::Stationary,
            Physics::Heat => Study::Stationary,
            Physics::Solid => Study::Static,
            Physics::Emw => Study::Modes { freq: 10e9, n_target: 1.0, evanescent: false, n_modes: 4 },
        };
        ProblemSpec { physics, order: 2, materials: Materials::default(), bcs: vec![], study, plane: Plane::Strain, thickness: 1.0, thermal: None, solver: LinearSolver::Auto, output: OutputSpec::default(), out_of_plane: None, inactive: vec![], periodic: vec![] }
    }
}

// ------------------------------------------------------------------ results

#[derive(Clone, Debug)]
pub struct Field {
    pub name: String,
    pub unit: String,
    /// value at every dof (P1: vertices; P2: vertices then edge midpoints)
    pub values: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct Solution {
    pub label: String,
    pub fields: Vec<Field>,
    /// displacement per dof, `dim` values each (solid mechanics)
    pub disp: Option<Vec<f64>>,
    /// time (s), frequency (Hz) or neff (re, im) of this solution
    pub value: Option<f64>,
    pub neff: Option<(f64, f64)>,
    /// the same derived quantities as `fields` before averaging over the cells at a node (wafer names), when the
    /// kernel's fields stage provides them
    pub cell_fields: Option<crate::fields::CellFields>,
}

impl Solution {
    pub fn field(&self, name: &str) -> Option<&Field> {
        self.fields.iter().find(|f| f.name == name)
    }
}

#[derive(Clone, Debug)]
pub enum Scalar {
    Num(f64),
    Text(String),
}

#[derive(Clone, Debug, Default)]
pub struct SolveResult {
    pub order: usize,
    pub dim: usize,
    pub solutions: Vec<Solution>,
    pub scalars: Vec<(String, Scalar)>,
    pub log: Vec<String>,
    pub seconds: f64,
    pub dofs: usize,
    /// polarization of the modes of an optical mode study (Stokes parameters, Poincare sphere), see `analysis`
    pub analysis: Option<crate::analysis::Analysis>,
}

impl SolveResult {
    pub fn scalar(&self, name: &str) -> Option<f64> {
        self.scalars.iter().find(|(n, _)| n.starts_with(name)).and_then(|(_, v)| match v {
            Scalar::Num(x) => Some(*x),
            _ => None,
        })
    }
}

/// State carried between solves (a heat result feeds thermal expansion).
#[derive(Clone, Debug, Default)]
pub struct SolveState {
    pub heat_t: Option<Vec<f64>>,
}
