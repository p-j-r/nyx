//! The project (what is saved) and the run-time state around it.
//!
//! A `Project` is the geometry (`GeometrySpec`: a `geom` feature list), one `ProblemSpec` per physics, the active physics and the
//! colour map. It is written as XML (`xmlio`) or JSON. `ResultView` wraps a `SolveResult` with what the plots need (selected solution
//! and field, range, deformation scale, cached values).

use crate::colormap::Colormap;
use femcore::fields::Display;
use femcore::problem::*;
use femcore::view::ViewMesh;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// The geometry of a project: one `geom` model (exact regions and edges with names and rules, see `femcore::geomdrive`). A 2D model is
/// the x-z cross-section, meshed with triangles of size `h`; with `extrude` layers it is swept along y (the propagation axis) into
/// tetrahedra, which is how 3D problems are made until `geom` meshes general 3D models.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "serde_json::Value")]
pub struct GeometrySpec {
    pub list: geom::model::FeatureList,
    /// global element size
    pub h: f64,
    pub extrude: Option<Vec<femcore::mesher::Layer>>,
}

impl TryFrom<serde_json::Value> for GeometrySpec {
    type Error = String;

    /// Reads the current form, and the `{"geometry": {"Model": ...}}` form of earlier files. The removed kinds (rectangle, block, shapes,
    /// solids) are named in the error.
    fn try_from(v: serde_json::Value) -> Result<Self, String> {
        #[derive(Deserialize)]
        struct Fields {
            list: geom::model::FeatureList,
            h: f64,
            #[serde(default)]
            extrude: Option<Vec<femcore::mesher::Layer>>,
        }
        let inner = match v.get("geometry") {
            Some(g) => match g.get("Model") {
                Some(m) => m.clone(),
                None => {
                    let kind = g.as_object().and_then(|o| o.keys().next().cloned()).unwrap_or_default();
                    return Err(format!(
                        "this project uses the geometry kind '{kind}', which has been removed: build the geometry as a model (Geometry tab) instead"
                    ));
                }
            },
            None => v,
        };
        let f: Fields = serde_json::from_value(inner).map_err(|e| e.to_string())?;
        Ok(GeometrySpec { list: f.list, h: f.h, extrude: f.extrude })
    }
}

impl GeometrySpec {
    /// A model of the x-z cross-section (a 2D feature list) meshed with element size `h`.
    pub fn model(list: geom::model::FeatureList, h: f64) -> Self {
        GeometrySpec { list, h, extrude: None }
    }

    /// A 2D model swept along y through `layers` (see `femcore::mesher::extrude_y`).
    pub fn model_extruded(list: geom::model::FeatureList, h: f64, layers: Vec<femcore::mesher::Layer>) -> Self {
        GeometrySpec { list, h, extrude: Some(layers) }
    }

    /// A rectangle `lx` wide (x) and `lz` high (z), one region `plate`, with element size `h`; swept `len` along y into a block when `len` is given.
    pub fn plate(lx: f64, lz: f64, h: f64, len: Option<f64>) -> Self {
        let mut m = geom::model::Model::new_2d();
        let r = m.rect("plate", (0.0, 0.0), (lx, lz));
        m.overlay(&[("plate", r)]);
        let layers = len.map(|l| vec![femcore::mesher::Layer { thickness: l, n: ((l / h).round() as usize).clamp(1, 200), domain_offset: 0 }]);
        GeometrySpec { list: m.list, h, extrude: layers }
    }

    /// Evaluates and meshes the model: the mesh with the names of its regions and boundaries.
    pub fn build_named(&self) -> Result<femcore::geomdrive::MeshedGeom, String> {
        if self.list.dim != 2 {
            return Err("a 3D model cannot be meshed yet: model the x-z cross-section in 2D and sweep it along y with layers".into());
        }
        if !(self.h > 0.0 && self.h.is_finite()) {
            return Err("the element size must be a positive finite number".into());
        }
        let mut g = femcore::geomdrive::mesh_geom_2d(femcore::geomdrive::evaluate(&self.list)?, &femcore::geomdrive::GeomMeshOpts::new(self.h))?;
        if let Some(layers) = &self.extrude {
            g.mesh = femcore::mesher::extrude_y(&g.mesh, layers, 0.0)?;
        }
        Ok(g)
    }

    pub fn build(&self) -> Result<femcore::mesh::SimplexMesh, String> {
        Ok(self.build_named()?.mesh)
    }
}


#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Source {
    Generated,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Tab {
    Geometry,
    Materials,
    Physics,
    Study,
    Results,
    /// polarization of the optical modes (Stokes parameters, Poincare sphere)
    Analysis,
}

impl Tab {
    pub const ALL: [Tab; 6] = [Tab::Geometry, Tab::Materials, Tab::Physics, Tab::Study, Tab::Results, Tab::Analysis];
    pub fn name(self) -> &'static str {
        match self {
            Tab::Geometry => "Geometry",
            Tab::Materials => "Materials",
            Tab::Physics => "Physics",
            Tab::Study => "Study",
            Tab::Results => "Results",
            Tab::Analysis => "Analysis",
        }
    }
}

/// Everything that is saved to a `.nyx.json` project file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Project {
    pub version: u32,
    pub geometry: GeometrySpec,
    pub source: Source,
    pub active: Physics,
    pub specs: Vec<ProblemSpec>,
    pub colormap: Colormap,
}

impl Project {
    /// The project as readable XML (`.nyx.xml`), written with quick-xml.
    pub fn to_xml(&self) -> Result<String, String> {
        let v = serde_json::to_value(self).map_err(|e| e.to_string())?;
        femcore::xmlio::value_to_xml("nyx_project", &v)
    }

    pub fn from_xml(xml: &str) -> Result<Project, String> {
        let (root, v) = femcore::xmlio::xml_to_value(xml)?;
        if root != "nyx_project" {
            return Err(format!("not a nyx project: root element is <{root}>"));
        }
        serde_json::from_value(v).map_err(|e| format!("not a nyx project: {e}"))
    }

    /// Reads a `.nyx.xml` or `.nyx.json` project, by file extension.
    pub fn load(path: &std::path::Path) -> Result<Project, String> {
        let s = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut p = if path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("xml")) {
            Project::from_xml(&s)?
        } else {
            serde_json::from_str::<Project>(&s).map_err(|e| format!("not a nyx project: {e}"))?
        };
        let frame = p.specs.first().map(|s| s.materials.frame).unwrap_or_default();
        for ph in Physics::ALL {
            if !p.specs.iter().any(|s| s.physics == ph) {
                let mut ns = ProblemSpec::new(ph);
                ns.materials.frame = frame;
                p.specs.push(ns);
            }
        }
        Ok(p)
    }

    pub fn spec(&self, p: Physics) -> &ProblemSpec {
        self.specs.iter().find(|s| s.physics == p).expect("spec for every physics")
    }
    pub fn spec_mut(&mut self, p: Physics) -> &mut ProblemSpec {
        self.specs.iter_mut().find(|s| s.physics == p).expect("spec for every physics")
    }
    pub fn active_spec(&self) -> &ProblemSpec {
        self.spec(self.active)
    }

    pub fn default_project() -> Self {
        crate::examples::by_index(0)
    }
}

/// Mesh and what we know about its origin.
#[derive(Default)]
pub struct Doc {
    pub mesh: Option<Arc<femcore::mesh::SimplexMesh>>,
    pub view: Option<Arc<ViewMesh>>,
    pub source: String,
    pub mesh_error: Option<String>,
}

/// One result on display, with the plot settings.
pub struct ResultView {
    pub name: String,
    pub physics: Option<Physics>,
    pub result: SolveResult,
    pub display: Display,
    pub sol: usize,
    pub field: usize,
    pub manual_range: Option<(f32, f32)>,
    pub warp: f32,
    pub auto_warp: f32,
    pub values_cache: Option<(usize, usize, Arc<Vec<f32>>)>,
    pub disp_cache: Option<(usize, Arc<Vec<f32>>)>,
}

impl ResultView {
    pub fn new(name: String, physics: Option<Physics>, result: SolveResult, display: Display, dim: usize) -> Self {
        // deformation scale such that the largest displacement is 10 % of the model size
        let mut auto_warp = 0.0;
        if let Some(d) = result.solutions.iter().filter_map(|s| s.disp.as_ref()).next() {
            let maxd = d.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            let (lo, hi) = bbox_of(&display);
            let size = (0..dim).map(|k| (hi[k] - lo[k]) as f64).fold(0.0, f64::max);
            if maxd > 0.0 && size > 0.0 {
                auto_warp = ((if result.solutions.len() > 1 { 0.2 } else { 0.1 }) * size / maxd) as f32;
            }
        }
        let field = result.solutions.first().and_then(|s| s.fields.iter().position(|f| f.name == "|u|" || f.name == "T" || f.name == "V" || f.name == "|E|")).unwrap_or(0);
        let has_disp = result.solutions.iter().any(|s| s.disp.is_some());
        ResultView { name, physics, result, display, sol: 0, field, manual_range: None, warp: if has_disp { auto_warp } else { 0.0 }, auto_warp, values_cache: None, disp_cache: None }
    }

    pub fn solution(&self) -> &Solution {
        &self.result.solutions[self.sol.min(self.result.solutions.len() - 1)]
    }

    pub fn field_values(&mut self) -> Arc<Vec<f32>> {
        let sol = self.sol.min(self.result.solutions.len() - 1);
        let nf = self.result.solutions[sol].fields.len();
        let field = self.field.min(nf.saturating_sub(1));
        if let Some((s, f, v)) = &self.values_cache {
            if *s == sol && *f == field {
                return v.clone();
            }
        }
        let v = Arc::new(self.display.pick(&self.result.solutions[sol].fields[field].values));
        self.values_cache = Some((sol, field, v.clone()));
        v
    }

    /// Displacement of each display point (3 values each) of the current solution.
    pub fn displacement(&mut self) -> Option<Arc<Vec<f32>>> {
        let sol = self.sol.min(self.result.solutions.len() - 1);
        if let Some((s, v)) = &self.disp_cache {
            if *s == sol {
                return Some(v.clone());
            }
        }
        let d = self.result.solutions[sol].disp.as_ref()?;
        let dim = self.result.dim;
        let mut out = vec![0f32; self.display.dofs.len() * 3];
        for (i, &dof) in self.display.dofs.iter().enumerate() {
            for k in 0..dim {
                out[3 * i + k] = d[dof * dim + k] as f32;
            }
        }
        let v = Arc::new(out);
        self.disp_cache = Some((sol, v.clone()));
        Some(v)
    }

    pub fn range(&mut self) -> (f32, f32) {
        if let Some(r) = self.manual_range {
            return r;
        }
        let v = self.field_values();
        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
        for &x in v.iter() {
            lo = lo.min(x);
            hi = hi.max(x);
        }
        if !lo.is_finite() || !hi.is_finite() {
            return (0.0, 1.0);
        }
        if hi - lo < 1e-30 * lo.abs().max(hi.abs()).max(1e-300) || hi == lo {
            let d = lo.abs().max(1e-6) * 0.05;
            return (lo - d, hi + d);
        }
        (lo, hi)
    }
}

fn bbox_of(d: &Display) -> ([f32; 3], [f32; 3]) {
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for p in &d.points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    (lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_example_survives_a_json_round_trip_and_builds() {
        for i in 0..crate::examples::NAMES.len() {
            let p = crate::examples::by_index(i);
            let s = serde_json::to_string_pretty(&p).unwrap();
            let q: Project = serde_json::from_str(&s).unwrap();
            assert_eq!(q.geometry, p.geometry);
            assert_eq!(q.specs, p.specs);
            q.geometry.build().unwrap();
        }
    }
}

#[cfg(test)]
mod xml_tests {
    use super::*;

    #[test]
    fn xml_round_trip_all_examples() {
        for i in 0..crate::examples::NAMES.len() {
            let p = crate::examples::by_index(i);
            let xml = p.to_xml().unwrap();
            let q = Project::from_xml(&xml).unwrap();
            assert_eq!(serde_json::to_value(&p).unwrap(), serde_json::to_value(&q).unwrap(), "example {i}");
        }
    }
}
