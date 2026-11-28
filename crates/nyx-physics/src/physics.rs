// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The `Physics` trait: everything the rest of nyx may ask of a physics (restructure plan, section 2; ADR 25).
//!
//! Only plain data crosses the trait: descriptors and specs are serde values, meshes and fields are flat arrays. That keeps a
//! separate-process plug-in possible later without changing the interface. The assemble, solve and fields stages of a kernel
//! stay inside the physics crate (ADR 9) and are public API of that crate; the trait is the uniform entry the engine uses.

use crate::desc::{BcDesc, FieldDesc, FieldLink, ParamDesc, PropDesc, Sizing, StudyDesc};
use crate::error::Error;
use crate::result::SolveResult;
use crate::sensitivity::{SensitivityReport, SensitivityRequest};
use crate::spec::ProblemSpec;
use crate::Mesh;
use femcore::progress::Progress;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Turns the names of a problem (`geom` region names and boundary selectors) into the integer labels of the mesh being solved.
/// Domains and boundaries are numbered from 1 (restructure phase 2). The engine implements it over the meshed geometry; a test
/// implements it over a table.
pub trait LabelResolver {
    /// The domain labels of a region name or a pattern of names (`sapphire*`).
    fn domains(&self, region: &str) -> Result<Vec<i32>, Error>;
    /// The boundary labels a selector picks (`outer & plane(x = 0)`).
    fn boundaries(&self, selector: &str) -> Result<Vec<i32>, Error>;
    /// Optional 1-based boundary label -> physical metal film thickness `t_m` (m) for finite-thickness `Op::Sheet` facets.
    fn sheet_thickness_m(&self) -> BTreeMap<i32, f64> {
        BTreeMap::new()
    }
    /// Optional 1-based boundary label -> substrate over-etch trench depth `t_e` (m) for etched `Op::Sheet` facets.
    fn sheet_etch_depth_m(&self) -> BTreeMap<i32, f64> {
        BTreeMap::new()
    }
}

/// A field another physics exported, kept on the mesh it was solved on (per-physics meshing, ADR 135). The engine fills
/// `Ctx::transferred_inputs` with one per imported role whenever the exporting physics was solved on a different mesh, and still puts the
/// field sampled at the dofs of the importing mesh into `Ctx::inputs`, so a kernel that only reads `inputs` keeps working unchanged. A kernel
/// that wants the exact source (to sample it at its own quadrature points, or per cell across an interface) reads this.
#[derive(Clone, Debug)]
pub struct TransferredField {
    /// the physics that exported it and the field name (`dfs`, `T`)
    pub physics: String,
    pub field: String,
    /// the mesh the field was solved on
    pub mesh: Arc<Mesh>,
    /// the simplex discretisation of `mesh` the values number (the engine shares the cached one; `None` makes one on demand)
    pub discretisation: Option<Arc<crate::Discretisation>>,
    /// canonical domain and boundary names of `mesh` (the `geom` names of `LabelMap`), by label
    pub domain_names: BTreeMap<i32, String>,
    pub boundary_names: BTreeMap<i32, String>,
    /// the element order of the values (1 or 2)
    pub order: usize,
    /// components per dof (1 for a scalar)
    pub components: usize,
    /// one value per dof of `mesh` at `order`, `components` interleaved per dof
    pub values: Vec<f64>,
    /// the un-averaged per-cell fields of the same solution, when it has them
    pub cell_fields: Vec<crate::result::CellField>,
}

impl TransferredField {
    /// Whether `mesh` is the mesh the field lives on (same vertices, cells and labels): the zero-cost path where no transfer is needed.
    pub fn is_same_mesh(&self, mesh: &Mesh) -> bool {
        femcore::transfer::same_exchange_mesh(&self.mesh, mesh)
    }

    fn source(&self) -> Result<Arc<crate::Discretisation>, Error> {
        match &self.discretisation {
            Some(d) => Ok(d.clone()),
            None => Ok(Arc::new(crate::discretise(&self.mesh)?)),
        }
    }

    /// The field at the local nodes of every cell of `target` (`nloc(order)` per cell, `components` interleaved), each node read only in
    /// source cells of the same canonical domain as its target cell (`target_domains`). `found` marks the nodes a source cell was found for.
    pub fn sample_cell_nodes(
        &self,
        target: &crate::Discretisation,
        target_domains: &BTreeMap<i32, String>,
        order: usize,
    ) -> Result<femcore::transfer::Sampled, Error> {
        use femcore::transfer::{SourceField, SourceMesh, TransferIndex, TransferOptions};
        let src = self.source()?;
        let index = TransferIndex::new(SourceMesh::Simplex(&src), &self.domain_names, TransferOptions::default())
            .map_err(|e| Error::invalid(e.to_string()))?;
        let field = SourceField::Nodal { values: &self.values, ncomp: self.components.max(1), order: self.order };
        index.sample_cell_nodes(&field, target, target_domains, order).map_err(|e| Error::invalid(e.to_string()))
    }

    /// The field at the dofs of `target` (`order`, `components` interleaved per dof): per-cell domain-filtered samples averaged at each dof.
    pub fn sample_dofs(
        &self,
        target: &crate::Discretisation,
        target_domains: &BTreeMap<i32, String>,
        order: usize,
    ) -> Result<femcore::transfer::Sampled, Error> {
        use femcore::transfer::{SourceField, SourceMesh, TransferIndex, TransferOptions};
        let src = self.source()?;
        let index = TransferIndex::new(SourceMesh::Simplex(&src), &self.domain_names, TransferOptions::default())
            .map_err(|e| Error::invalid(e.to_string()))?;
        let field = SourceField::Nodal { values: &self.values, ncomp: self.components.max(1), order: self.order };
        index.sample_dofs(&field, target, target_domains, order).map_err(|e| Error::invalid(e.to_string()))
    }
}

/// What a solve may use besides the mesh and the problem: the label resolver, progress and cancellation, and the fields that other
/// physics exported earlier in a one-way chain.
pub struct Ctx<'a> {
    pub resolver: &'a dyn LabelResolver,
    pub progress: Arc<Progress>,
    /// the material files a problem can name (`MaterialAssign::Named`); `None` where there are none (a unit test)
    pub library: Option<Arc<femcore::material::Library>>,
    /// the free-space wavelength (m) of the study, for materials given as optical models (the engine derives it from the frequency of the study)
    pub wavelength: Option<f64>,
    /// the reference frequency (Hz) of the study, for microwave/RF dielectric and superconducting models
    pub frequency: Option<f64>,
    /// the reference temperature (K) of the study, for temperature-dependent thermal, dielectric and thermo-optic models
    pub temperature: Option<f64>,
    /// imported fields by role (`temperature`), values per dof, from physics that ran before this one
    pub inputs: BTreeMap<String, Vec<f64>>,
    /// imported fields by role that were solved on another mesh (per-physics meshing, ADR 135), with that mesh; empty when every physics of
    /// the chain shares one mesh. `inputs` still holds the same roles sampled on this solve's mesh
    pub transferred_inputs: BTreeMap<String, TransferredField>,
}

impl<'a> Ctx<'a> {
    pub fn new(resolver: &'a dyn LabelResolver) -> Ctx<'a> {
        Ctx {
            resolver,
            progress: Progress::new(),
            library: None,
            wavelength: None,
            frequency: None,
            temperature: None,
            inputs: BTreeMap::new(),
            transferred_inputs: BTreeMap::new(),
        }
    }

    /// `Err(Cancelled)` when the run was cancelled; kernels call it between stages.
    pub fn check(&self) -> Result<(), Error> {
        if self.progress.is_cancelled() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
}

pub trait Physics: Send + Sync {
    /// The three-letter code (ADR 12): `dfs`, `sms`, `ems`.
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;

    /// The kinds of study, with their parameters and the dimensions they run in.
    fn studies(&self) -> Vec<StudyDesc>;
    /// The boundary and domain conditions, with their parameters.
    fn bc_kinds(&self) -> Vec<BcDesc>;
    /// The material properties the physics reads (name, unit, form, which study needs it).
    fn properties(&self) -> Vec<PropDesc>;
    /// Physics-level options shared by its studies (`ProblemSpec::params`), for example the plane assumption of a 2D solid.
    fn params(&self) -> Vec<ParamDesc>;
    /// The fields a study produces, with their components in wafer axes.
    fn fields(&self, study: &str) -> Vec<FieldDesc>;
    /// Fields other physics may import (one-way chains only in version 1).
    fn exports(&self) -> Vec<FieldLink>;
    /// Fields this physics takes from other physics (for example `sms` takes a temperature).
    fn imports(&self) -> Vec<FieldLink>;
    /// Exports produced when running a specific `study` (defaults to [`Self::exports`]).
    fn study_exports(&self, _study: &str) -> Vec<FieldLink> {
        self.exports()
    }
    /// Imports consumed when running a specific `study` (defaults to [`Self::imports`]).
    fn study_imports(&self, _study: &str) -> Vec<FieldLink> {
        self.imports()
    }

    /// Checks the physics-specific part of a problem against `mesh`; the structural checks against the descriptors are `validate_spec`.
    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error>;
    /// What the mesher needs to size the elements before any solution exists.
    fn sizing(&self, spec: &ProblemSpec, library: Option<&Arc<femcore::material::Library>>) -> Sizing;
    /// Solves the problem on `mesh`.
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error>;

    /// Evaluates a scalar design objective and its exact discrete adjoint / eigenvalue gradients
    /// on `mesh` with respect to [`SensitivityRequest::params`] (ADR 142).
    fn sensitivity(
        &self,
        _mesh: &Mesh,
        _spec: &ProblemSpec,
        _ctx: &mut Ctx,
        _req: &SensitivityRequest,
    ) -> Result<SensitivityReport, Error> {
        Err(Error::Unsupported(format!("physics '{}' does not implement adjoint sensitivities", self.id())))
    }
}
