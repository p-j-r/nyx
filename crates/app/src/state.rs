//! The project (what is saved) and the run-time state around it.
//!
//! A `Project` is the geometry (`GeometrySpec`: rectangle, block, 2D shapes, 3D solids, or a model file), one `ProblemSpec`
//! per physics, the active physics and the colour map. It is written as XML (`xmlio`) or JSON. `ResultView` wraps a
//! `SolveResult` with what the plots need (selected solution and field, range, deformation scale, cached values).

use crate::colormap::Colormap;
use femcore::fields::Display;
use femcore::problem::*;
use femcore::view::ViewMesh;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Geometry {
    /// structured triangles on the x-z cross-section (width `lx`, height `lz`)
    Rectangle {
        lx: f64,
        #[serde(alias = "ly")]
        lz: f64,
        nx: usize,
        #[serde(alias = "ny")]
        nz: usize,
    },
    Block { lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize },
    /// 2D shapes meshed with triangles, optionally extruded through layers into 3D
    Shapes { spec: femcore::mesher::Spec2d, extrude: Option<Vec<femcore::mesher::Layer>> },
    /// 3D bodies (box, sphere, cylinder, extrusion) meshed with tetrahedra by isosurface stuffing
    Solids { spec: femcore::mesher3d::Spec3d },
    /// A `geom` model: exact regions and edges with names and rules (see `femcore::geomdrive`). A 2D model is the x-z cross-section,
    /// meshed with triangles of size `h`; with `extrude` layers it is swept along y (the propagation axis) into tetrahedra.
    Model { list: geom::model::FeatureList, h: f64, extrude: Option<Vec<femcore::mesher::Layer>> },
}

impl Geometry {
    pub fn dim(&self) -> usize {
        match self {
            Geometry::Rectangle { .. } => 2,
            Geometry::Block { .. } => 3,
            Geometry::Shapes { extrude, .. } => if extrude.is_some() { 3 } else { 2 },
            Geometry::Solids { .. } => 3,
            Geometry::Model { list, extrude, .. } => if extrude.is_some() { 3 } else { list.dim as usize },
        }
    }
}

/// A box (min per axis, then max per axis) that relabels the cells or boundary facets inside it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Region {
    pub bbox: Vec<f64>,
    pub label: i32,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GeometrySpec {
    pub geometry: Geometry,
    pub domains: Vec<Region>,
    pub boundaries: Vec<Region>,
}

impl GeometrySpec {
    pub fn rectangle(lx: f64, lz: f64, nx: usize, nz: usize) -> Self {
        GeometrySpec { geometry: Geometry::Rectangle { lx, lz, nx, nz }, domains: vec![], boundaries: vec![] }
    }
    /// A model of the x-z cross-section (a 2D feature list) meshed with element size `h`.
    pub fn model(list: geom::model::FeatureList, h: f64) -> Self {
        GeometrySpec { geometry: Geometry::Model { list, h, extrude: None }, domains: vec![], boundaries: vec![] }
    }
    /// A 2D model swept along y through `layers` (see `femcore::mesher::extrude_y`).
    pub fn model_extruded(list: geom::model::FeatureList, h: f64, layers: Vec<femcore::mesher::Layer>) -> Self {
        GeometrySpec { geometry: Geometry::Model { list, h, extrude: Some(layers) }, domains: vec![], boundaries: vec![] }
    }
    pub fn shapes() -> Self {
        use femcore::mesher::*;
        let spec = Spec2d { parts: vec![Part { shape: Shape::Rect { x0: 0.0, y0: 0.0, x1: 1.0, y1: 0.5 }, domain: 1, h: None }], h: 0.05, min_angle_deg: 25.0 };
        GeometrySpec { geometry: Geometry::Shapes { spec, extrude: None }, domains: vec![], boundaries: vec![] }
    }
    pub fn solids() -> Self {
        use femcore::mesher3d::*;
        let spec = Spec3d { bodies: vec![Body3d { solid: Solid::Box { min: [0.0, 0.0, 0.0], max: [1.0, 0.5, 0.25] }, domain: 1, h: None }], h: 0.05, alpha: 0.3, snap_features: true, min_layers: 2.0 };
        GeometrySpec { geometry: Geometry::Solids { spec }, domains: vec![], boundaries: vec![] }
    }
    pub fn block(lx: f64, ly: f64, lz: f64, nx: usize, ny: usize, nz: usize) -> Self {
        GeometrySpec { geometry: Geometry::Block { lx, ly, lz, nx, ny, nz }, domains: vec![], boundaries: vec![] }
    }

    pub fn build(&self) -> Result<femcore::mesh::SimplexMesh, String> {
        use femcore::mesh::SimplexMesh;
        let mut mesh = match &self.geometry {
            Geometry::Rectangle { lx, lz, nx, nz } => {
                if !(*lx > 0.0 && *lz > 0.0 && lx.is_finite() && lz.is_finite()) || *nx < 1 || *nz < 1 {
                    return Err("sizes must be positive finite numbers and element counts at least 1".into());
                }
                if nx * nz > 4_000_000 {
                    return Err("that is too many elements (limit 4 million)".into());
                }
                SimplexMesh::rectangle(*lx, *lz, *nx, *nz)
            }
            Geometry::Shapes { spec, extrude } => {
                let m2 = femcore::mesher::mesh_2d(spec)?;
                match extrude {
                    Some(layers) => femcore::mesher::extrude_y(&m2, layers, 0.0)?,
                    None => m2,
                }
            }
            Geometry::Model { list, h, extrude } => {
                if list.dim != 2 {
                    return Err("a 3D model cannot be meshed yet: model the x-z cross-section in 2D and extrude it along y with layers".into());
                }
                let g = femcore::geomdrive::mesh_geom_2d(femcore::geomdrive::evaluate(list)?, &femcore::geomdrive::GeomMeshOpts::new(*h))?;
                match extrude {
                    Some(layers) => femcore::mesher::extrude_y(&g.mesh, layers, 0.0)?,
                    None => g.mesh,
                }
            }
            Geometry::Solids { spec } => femcore::mesher3d::mesh_3d(spec)?.0,
            Geometry::Block { lx, ly, lz, nx, ny, nz } => {
                if !(*lx > 0.0 && *ly > 0.0 && *lz > 0.0 && lx.is_finite() && ly.is_finite() && lz.is_finite()) || *nx < 1 || *ny < 1 || *nz < 1 {
                    return Err("sizes must be positive finite numbers and element counts at least 1".into());
                }
                if nx * ny * nz * 6 > 6_000_000 {
                    return Err("that is too many elements (limit 6 million tetrahedra)".into());
                }
                SimplexMesh::block(*lx, *ly, *lz, *nx, *ny, *nz)
            }
        };
        let dom: Vec<(Vec<f64>, i32)> = self.domains.iter().map(|r| (r.bbox.clone(), r.label)).collect();
        let bnd: Vec<(Vec<f64>, i32)> = self.boundaries.iter().map(|r| (r.bbox.clone(), r.label)).collect();
        mesh.relabel(&dom, &bnd)?;
        Ok(mesh)
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
