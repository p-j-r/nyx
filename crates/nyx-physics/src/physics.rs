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
    /// imported fields by role (`temperature`), values per dof, from physics that ran before this one
    pub inputs: BTreeMap<String, Vec<f64>>,
}

impl<'a> Ctx<'a> {
    pub fn new(resolver: &'a dyn LabelResolver) -> Ctx<'a> {
        Ctx { resolver, progress: Progress::new(), library: None, wavelength: None, inputs: BTreeMap::new() }
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

    /// Checks the physics-specific part of a problem against `mesh`; the structural checks against the descriptors are `validate_spec`.
    fn validate(&self, spec: &ProblemSpec, mesh: &Mesh) -> Result<(), Error>;
    /// What the mesher needs to size the elements before any solution exists.
    fn sizing(&self, spec: &ProblemSpec, library: Option<&Arc<femcore::material::Library>>) -> Sizing;
    /// Solves the problem on `mesh`.
    fn solve(&self, mesh: &Mesh, spec: &ProblemSpec, ctx: &mut Ctx) -> Result<SolveResult, Error>;
}
