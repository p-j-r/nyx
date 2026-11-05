//! The project: what is saved. A geometry (a `geom` feature list), the mesh settings, one problem (`ProblemSpec` envelope) per physics the
//! project loads, and which physics is active. It is plain serde data; `to_json` writes it with sorted keys so the same project is the same
//! text (and hash) in every build. The JSON layout is provisional until phase 8 derives the final one (gate G6): there is no `format` key yet.

use geom::model::FeatureList;
use nyx_physics::{canonical_json, Error, ProblemSpec, Registry};
use serde::{Deserialize, Serialize};

/// One layer of a model swept along y into a 3D mesh (the 2D cross-section is extruded with `n` cells across `thickness`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SweepLayer {
    pub thickness: f64,
    pub n: usize,
    /// added to the domain numbers of this layer's cells
    #[serde(default)]
    pub domain_offset: i32,
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
    /// rebuild the absorbing layers of the model as structured layers (quadrilaterals, prisms, hexahedra, cut into simplices)
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub structured_pml: bool,
    /// largest distance of a chord from a curved edge or face in metres; absent: h / 50
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_dev: Option<f64>,
}

impl Default for MeshSettings {
    fn default() -> Self {
        MeshSettings {
            h: 1e-7,
            sweep: None,
            adapt_scale: None,
            physics_sizing: false,
            points_per_wavelength: None,
            structured_pml: false,
            max_dev: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Project {
    pub title: String,
    pub geometry: FeatureList,
    pub mesh: MeshSettings,
    /// one problem per physics the project loads (one to three)
    #[serde(default)]
    pub physics: Vec<ProblemSpec>,
    /// the code of the physics shown first and solved by default
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub active: Option<String>,
    /// materials the project carries, by name, as TOML text (ADR 29); a problem names these or the library's. They may `extends` a library material
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub materials: std::collections::BTreeMap<String, String>,
}

impl Project {
    /// An empty project: no features, no physics.
    pub fn new(title: &str, dim: u8) -> Project {
        Project {
            title: title.into(),
            geometry: FeatureList::new(dim),
            mesh: MeshSettings::default(),
            physics: vec![],
            active: None,
            materials: Default::default(),
        }
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

    pub fn from_json(text: &str) -> Result<Project, Error> {
        serde_json::from_str(text).map_err(|e| Error::invalid(format!("not a nyx project: {e}")))
    }

    /// Everything that can be checked without a mesh: the feature list, the element size, the physics against the registry and each problem
    /// against the descriptors of its physics. The first problem found is returned, in plain words.
    pub fn check(&self, registry: &Registry) -> Result<(), Error> {
        self.geometry.validate().map_err(|e| Error::invalid(format!("geometry: {e}")))?;
        if !(self.mesh.h > 0.0 && self.mesh.h.is_finite()) {
            return Err(Error::invalid("the element size must be a positive finite number"));
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
        Ok(())
    }
}
