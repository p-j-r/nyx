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

pub use crate::props::*;
use serde::{Deserialize, Serialize};

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
        ProblemSpec { physics, order: 2, materials: Materials::default(), bcs: vec![], study, plane: Plane::Strain, thickness: 1.0, thermal: None, solver: LinearSolver::Auto, output: OutputSpec::default(), inactive: vec![], periodic: vec![] }
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
