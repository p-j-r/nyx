// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! The project: what is saved. A geometry (a `geom` feature list), the mesh settings, one problem (`ProblemSpec` envelope) per physics the
//! project loads, and which physics is active. It is plain serde data; `to_json` writes it with sorted keys so the same project is the same
//! text (and hash) in every build. The project, its feature list and each of its problems carry a `format` key (first value 1, ADR 63); a file
//! of another format, or without the key, is refused with a plain message (there is no migration code).

pub use femcore::meshing::request::PmlModeRequest;
use geom::model::FeatureList;
use nyx_physics::{canonical_json, Error, ProblemSpec, Registry};
use serde::{Deserialize, Serialize};

fn is_auto_pml_mode(m: &PmlModeRequest) -> bool {
    *m == PmlModeRequest::Auto
}
fn default_true() -> bool {
    true
}
fn is_true(b: &bool) -> bool {
    *b
}

/// One layer of a model swept along y into a 3D mesh (the 2D cross-section is extruded with `n` cells across `thickness`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SweepLayer {
    pub thickness: f64,
    pub n: usize,
    /// added to the domain numbers of this layer's cells
    #[serde(default)]
    pub domain_offset: i32,
    /// name of a region of the 2D model: every cell of this layer belongs to that region (a cladding layer above or below a patterned
    /// slab); the offset is then not used
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fill: Option<String>,
    /// optional per-layer region remap (`source_region_name -> target_region_name` in this layer)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub regions: Option<std::collections::BTreeMap<String, String>>,
}

/// How the geometry is meshed.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MeshSettings {
    /// global element size in metres
    pub h: f64,
    /// sweep the 2D cross-section along y through these layers (3D mesh)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sweep: Option<Vec<SweepLayer>>,
    /// adaptive meshing: the factor on the sizes the model asks for (smaller is finer); `None` is the simple mesher
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_scale: Option<f64>,
    /// size the mesh from the physics (points per wavelength by element order, shorter in high-index or slow-wave regions): the engine asks
    /// every physics of the project for its `Sizing` and the mesher resolves the waves (`femcore::meshing::request`, ADR 50, 54)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub physics_sizing: bool,
    /// points per shortest wavelength when `physics_sizing` is on; absent: 12 for order 1, 6 for order 2
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub points_per_wavelength: Option<f64>,
    /// how absorbing layers of the model are meshed (`"auto"` [default], `"structured"`, or `"inline"`)
    #[serde(default, skip_serializing_if = "is_auto_pml_mode")]
    pub pml_mode: PmlModeRequest,
    /// legacy boolean alias (`true` -> `Structured`, `false` -> `Inline` when `pml_mode == Auto`) kept for backward compatibility
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub structured_pml: Option<bool>,
    /// largest distance of a chord from a curved edge or face in metres; absent: h / 50
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_dev: Option<f64>,
    /// 3D: give the cells on a curved surface or fitted sheet arc rim quadratic (isoparametric) geometry (`CurvedTet`, ADR 82, AF-5)
    /// instead of a polyhedron; used by `ems`, `dfs`, and `sheet3d`
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub curved: bool,
    /// 3D sheet stack (`sheet3d`): target minimum number of 2D triangles across narrow coplanar gaps between conductor sheets; `None` uses
    /// only `h` and explicit `add_size` rules
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap_cells: Option<f64>,
    /// 3D sheet stack (`sheet3d`): progressively coarsen the 2D sheet triangulation away from the conductor plane along `z` via crack-free
    /// pure-tetrahedral transition slabs (default `false`)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub z_coarsen: bool,
    /// adaptive / AMR metric gradation slope (`> 0`); `None` uses `0.4`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_grade: Option<f64>,
    /// adaptive / AMR cell cap (`>= 1`); `None` uses `2_000_000`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_max_cells: Option<usize>,
    /// maximum number of cells allowed before aborting meshing; absent: disabled (or `2_000_000` for adaptive meshing)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_cells: Option<usize>,
    /// 3D sheet stack (`sheet3d`): floor on initial vertical z-slab thickness near sheet planes in metres
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dz_min: Option<f64>,
    /// 3D sheet stack (`sheet3d`): maximum total number of horizontal z-slabs across active stack layers
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_z_slabs: Option<usize>,
    /// 3D sheet stack (`sheet3d`): geometric growth ratio for z-slab thickness away from sheet planes (default `1.35`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub z_growth_ratio: Option<f64>,
    /// 2D Ruppert refinement: cap on additional Steiner vertices as a multiple of initial vertices
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_steiner_multiplier: Option<usize>,
    /// fail meshing if 2D Ruppert refinement terminates early or drops unhealed collinear sliver triangles (default `true`)
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub strict_mesh: bool,
    /// 3D sheet stack (`sheet3d`): edge singularity enrichment (`"none"` or `"quarter_point"`); `None` uses study defaults
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub singular_edges: Option<String>,
    /// 3D sheet stack (`sheet3d`): Barsoum quarter-point mid-edge parameter `xi_QP` in `[0.25, 0.50]`; `None` uses `0.25`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub quarter_point_xi: Option<f64>,
    /// 3D sheet stack (`sheet3d`): radical Morrison-Lewis corner grading toward sharp 2D sheet corners (default `false`)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub corner_grading: bool,
    /// target discretization tolerance used to certify error bounds (e.g. Prager–Synge `flux_bounds`) when no `convergence_check` ran
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_tol: Option<f64>,
    /// number of goal-oriented AMR passes for `ems::cavity` or `sms` studies (`None` or `1` disables outer AMR loop)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_passes: Option<usize>,
    /// relative convergence tolerance between consecutive outer AMR passes for `ems::cavity` / `sms`
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_rel_tol: Option<f64>,
    /// cell growth factor per outer AMR pass for `ems::cavity` / `sms` (`> 1.0`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapt_growth: Option<f64>,
    /// retain transient cellwise `ElectroAmrData` indicators on `SolveResult` instead of stripping them after solve
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub retain_amr_data: bool,
    /// mesh only the regions named here (region names or `name*` patterns of the geometry, ADR 135): the cells of the other regions are
    /// dropped after meshing, the labels and names of the model are kept, so a physics can mesh its own subset (the solid regions and their
    /// absorbing layers of a solid-mechanics problem) and still exchange fields by canonical region name. `None` meshes every region
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub active_domains: Option<Vec<String>>,
}

impl Default for MeshSettings {
    fn default() -> Self {
        MeshSettings {
            h: 1e-7,
            sweep: None,
            adapt_scale: None,
            physics_sizing: false,
            points_per_wavelength: None,
            pml_mode: PmlModeRequest::Auto,
            structured_pml: None,
            max_dev: None,
            curved: false,
            gap_cells: None,
            z_coarsen: false,
            adapt_grade: None,
            adapt_max_cells: None,
            max_cells: None,
            dz_min: None,
            max_z_slabs: None,
            z_growth_ratio: None,
            max_steiner_multiplier: None,
            strict_mesh: true,
            singular_edges: None,
            quarter_point_xi: None,
            corner_grading: false,
            target_tol: None,
            adapt_passes: None,
            adapt_rel_tol: None,
            adapt_growth: None,
            retain_amr_data: false,
            active_domains: None,
        }
    }
}

impl MeshSettings {
    /// Resolves effective [`PmlModeRequest`], honoring legacy `structured_pml: Some(true | false)`
    /// when `pml_mode` was left at its default `Auto`.
    pub fn resolved_pml_mode(&self) -> PmlModeRequest {
        match (self.pml_mode, self.structured_pml) {
            (PmlModeRequest::Auto, Some(true)) => PmlModeRequest::Structured,
            (PmlModeRequest::Auto, Some(false)) => PmlModeRequest::Inline,
            (mode, _) => mode,
        }
    }
}

/// The file format of a project (ADR 63).
pub const PROJECT_FORMAT: u32 = 1;

/// The checks of mesh settings that need no mesh, in plain words.
fn check_mesh_settings(s: &MeshSettings) -> Result<(), String> {
    if !(s.h > 0.0 && s.h.is_finite()) {
        return Err("the element size must be a positive finite number".into());
    }
    if s.gap_cells.is_some_and(|gc| !(gc > 0.0 && gc.is_finite())) {
        return Err("`gap_cells` must be a positive finite number".into());
    }
    if s.max_cells.is_some_and(|mc| mc < 1) {
        return Err("`max_cells` must be an integer >= 1".into());
    }
    if s.dz_min.is_some_and(|dz| !(dz > 0.0 && dz.is_finite())) {
        return Err("`dz_min` must be a positive finite number".into());
    }
    if s.max_z_slabs.is_some_and(|ms| ms < 1) {
        return Err("`max_z_slabs` must be an integer >= 1".into());
    }
    if s.z_growth_ratio.is_some_and(|zr| !(zr > 1.0 && zr.is_finite())) {
        return Err("`z_growth_ratio` must be a finite number strictly greater than 1.0".into());
    }
    if let Some(a) = &s.active_domains {
        if a.is_empty() || a.iter().any(|n| n.trim().is_empty()) {
            return Err("`active_domains` must name at least one region, and no empty name".into());
        }
    }
    Ok(())
}

fn current_format() -> u32 {
    PROJECT_FORMAT
}

fn default_coupled_max_iters() -> usize {
    25
}
fn default_coupled_rtol() -> f64 {
    1e-6
}
fn default_anderson_m() -> usize {
    5
}
fn default_source_scale() -> f64 {
    1.0
}
fn is_one_f64(v: &f64) -> bool {
    (*v - 1.0).abs() <= f64::EPSILON
}

/// Explicit opt-in configuration for two-way Anderson-accelerated self-consistent coupled loops
/// (photothermal `ems <-> dfs` or electrothermal `dfs::conduction <-> dfs::heat_stationary`, ADR 139).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoupledLoopConfig {
    /// Maximum outer Anderson/Picard iterations (`>= 1`, default `25`).
    #[serde(default = "default_coupled_max_iters")]
    pub max_iters: usize,
    /// Relative $L_\infty$ temperature convergence tolerance `||T^{k+1} - T^k||_inf / max(1, ||T^k||_inf) <= rtol` (default `1e-6`).
    #[serde(default = "default_coupled_rtol")]
    pub rtol: f64,
    /// Anderson type-II history depth `m` (`0` = plain Picard, default `5`).
    #[serde(default = "default_anderson_m")]
    pub anderson_m: usize,
    /// Optional auxiliary source problem (e.g. `dfs` `conduction` paired with `Project`'s `dfs` `heat_stationary`
    /// for intra-physics electrothermal `V <-> T` coupling within a single `Project`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub auxiliary: Option<ProblemSpec>,
    /// Optional multiplier applied to exported `"heat_source"` (`Q_opt` or `Q_joule`) before thermal assembly (default `1.0`).
    #[serde(default = "default_source_scale", skip_serializing_if = "is_one_f64")]
    pub source_scale: f64,
}

impl Default for CoupledLoopConfig {
    fn default() -> Self {
        Self {
            max_iters: default_coupled_max_iters(),
            rtol: default_coupled_rtol(),
            anderson_m: default_anderson_m(),
            auxiliary: None,
            source_scale: 1.0,
        }
    }
}

impl CoupledLoopConfig {
    pub fn new(max_iters: usize, rtol: f64, anderson_m: usize) -> Self {
        Self { max_iters, rtol, anderson_m, auxiliary: None, source_scale: 1.0 }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Project {
    /// the file format, `1` (ADR 63)
    #[serde(default = "current_format")]
    pub format: u32,
    pub title: String,
    pub geometry: FeatureList,
    pub mesh: MeshSettings,
    /// mesh settings of a physics that is meshed on its own (per-physics meshing, ADR 135), by physics code; a physics not listed uses
    /// `mesh`. Fields between physics on different meshes are transferred by canonical region name (`femcore::transfer`)
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub mesh_by_physics: std::collections::BTreeMap<String, MeshSettings>,
    /// one problem per physics the project loads (one to three)
    #[serde(default)]
    pub physics: Vec<ProblemSpec>,
    /// the code of the physics shown first and solved by default
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub active: Option<String>,
    /// materials the project carries, by name, as TOML text (ADR 29); a problem names these or the library's. They may `extends` a library material
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub materials: std::collections::BTreeMap<String, String>,
    /// when true, out-of-range material temperature evaluations refuse the solve instead of clamping and emitting a warning
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub strict_material_range: bool,
    /// optional explicit opt-in configuration for two-way self-consistent Anderson coupled loops (`ems <-> dfs` or `conduction <-> heat`)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coupled_loop: Option<CoupledLoopConfig>,
}

impl Project {
    /// An empty project: no features, no physics.
    pub fn new(title: &str, dim: u8) -> Project {
        Project {
            format: PROJECT_FORMAT,
            title: title.into(),
            geometry: FeatureList::new(dim),
            mesh: MeshSettings::default(),
            mesh_by_physics: Default::default(),
            physics: vec![],
            active: None,
            materials: Default::default(),
            strict_material_range: false,
            coupled_loop: None,
        }
    }

    /// The mesh settings `physics` is meshed with: its own entry of `mesh_by_physics`, else the project's `mesh`.
    pub fn mesh_for(&self, physics: &str) -> &MeshSettings {
        self.mesh_by_physics.get(physics).unwrap_or(&self.mesh)
    }

    /// Whether `physics` is meshed with settings other than the project's `mesh`.
    pub fn has_own_mesh(&self, physics: &str) -> bool {
        self.mesh_by_physics.get(physics).is_some_and(|s| *s != self.mesh)
    }

    pub fn spec(&self, physics: &str) -> Option<&ProblemSpec> {
        self.physics.iter().find(|s| s.physics == physics)
    }

    pub fn spec_mut(&mut self, physics: &str) -> Option<&mut ProblemSpec> {
        self.physics.iter_mut().find(|s| s.physics == physics)
    }

    /// The spec of a physics the project must have, with the list of those it has in the message when it has not.
    pub fn need_spec(&mut self, physics: &str) -> Result<&mut ProblemSpec, Error> {
        let have = self.physics.iter().map(|s| s.physics.clone()).collect::<Vec<_>>().join(", ");
        self.spec_mut(physics).ok_or_else(|| {
            Error::invalid(format!(
                "the project has no physics '{physics}' (it has: {}); add it first",
                if have.is_empty() { "none" } else { &have }
            ))
        })
    }

    /// The codes to solve when none is named: the active physics, else all of them.
    pub fn default_targets(&self) -> Vec<String> {
        match &self.active {
            Some(a) => vec![a.clone()],
            None => self.physics.iter().map(|s| s.physics.clone()).collect(),
        }
    }

    /// The project as JSON text with sorted keys.
    pub fn to_json(&self) -> String {
        let v = serde_json::to_value(self).expect("a Project serialises");
        serde_json::to_string_pretty(&canonical_json(&v)).expect("a Value prints")
    }

    /// Reads a project file. The `format` key must be present (in the project, the feature list and every problem) and supported: a file of
    /// the layout before the key existed, or of a newer one, is refused (ADR 63), not converted.
    pub fn from_json(text: &str) -> Result<Project, Error> {
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        let v: serde_json::Value =
            serde_json::from_str(text).map_err(|e| Error::invalid(format!("not a nyx project: {e}")))?;
        let need = |what: &str, x: &serde_json::Value| -> Result<(), Error> {
            if x.get("format").is_none() {
                return Err(Error::invalid(format!(
                    "not a nyx project of the current layout: {what} has no \"format\" key (this program reads format {PROJECT_FORMAT}; older files are not converted)"
                )));
            }
            Ok(())
        };
        need("the project", &v)?;
        need("its geometry", &v["geometry"])?;
        for (i, s) in v["physics"].as_array().into_iter().flatten().enumerate() {
            need(&format!("problem {}", i + 1), s)?;
        }
        if let Some(aux) = v.get("coupled_loop").and_then(|c| c.get("auxiliary")).filter(|x| !x.is_null()) {
            need("coupled_loop.auxiliary", aux)?;
        }
        let p: Project = serde_json::from_value(v).map_err(|e| Error::invalid(format!("not a nyx project: {e}")))?;
        p.check_formats()?;
        Ok(p)
    }

    /// Refuses a format this program does not read, in the project, its feature list and its problems.
    pub fn check_formats(&self) -> Result<(), Error> {
        let bad = |what: &str, f: u32, ok: u32| {
            Err(Error::invalid(format!("{what} has format {f}; this program reads format {ok}")))
        };
        if self.format != PROJECT_FORMAT {
            return bad("the project", self.format, PROJECT_FORMAT);
        }
        if self.geometry.format != geom::model::FORMAT {
            return bad("the geometry", self.geometry.format, geom::model::FORMAT);
        }
        for s in &self.physics {
            if s.format != nyx_physics::PROBLEM_FORMAT {
                return bad(&format!("the '{}' problem", s.physics), s.format, nyx_physics::PROBLEM_FORMAT);
            }
        }
        if let Some(aux) = self.coupled_loop.as_ref().and_then(|c| c.auxiliary.as_ref()) {
            if aux.format != nyx_physics::PROBLEM_FORMAT {
                return bad("the coupled_loop auxiliary problem", aux.format, nyx_physics::PROBLEM_FORMAT);
            }
        }
        Ok(())
    }

    /// Everything that can be checked without a mesh: the feature list, the element size, the physics against the registry and each problem
    /// against the descriptors of its physics. The first problem found is returned, in plain words.
    pub fn check(&self, registry: &Registry) -> Result<(), Error> {
        self.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        check_mesh_settings(&self.mesh).map_err(Error::invalid)?;
        for (code, s) in &self.mesh_by_physics {
            if self.spec(code).is_none() {
                return Err(Error::invalid(format!(
                    "`mesh_by_physics` has settings for '{code}', which has no problem in the project"
                )));
            }
            check_mesh_settings(s).map_err(|e| Error::invalid(format!("the mesh of '{code}': {e}")))?;
        }
        for (i, s) in self.physics.iter().enumerate() {
            if self.physics[..i].iter().any(|o| o.physics == s.physics) {
                return Err(Error::invalid(format!("the project has two problems for '{}'", s.physics)));
            }
            let p = crate::physics_for(registry, &s.physics)?;
            nyx_physics::validate_spec(p, s)?;
        }
        if let Some(a) = &self.active {
            if self.spec(a).is_none() {
                return Err(Error::invalid(format!("the active physics '{a}' has no problem in the project")));
            }
        }
        if let Some(cfg) = &self.coupled_loop {
            if cfg.max_iters < 1 {
                return Err(Error::invalid("`coupled_loop.max_iters` must be an integer >= 1"));
            }
            if !(cfg.rtol > 0.0 && cfg.rtol.is_finite()) {
                return Err(Error::invalid("`coupled_loop.rtol` must be a positive finite number"));
            }
            if cfg.anderson_m > 64 {
                return Err(Error::invalid("`coupled_loop.anderson_m` must be <= 64"));
            }
            if !(cfg.source_scale >= 0.0 && cfg.source_scale.is_finite()) {
                return Err(Error::invalid("`coupled_loop.source_scale` must be a non-negative finite number"));
            }
            if let Some(aux) = &cfg.auxiliary {
                let p = crate::physics_for(registry, &aux.physics)?;
                nyx_physics::validate_spec(p, aux)?;
            }
        }
        Ok(())
    }
}
