//! The run-time state of the window around the engine's project: the mesh on show (`Doc`) and a result on show with its plot settings
//! (`ResultView`). The project itself, with its geometry, mesh settings and one problem per physics, is `nyx_engine::Project` and lives in the
//! engine; the window edits it through engine commands.

use nyx_engine::analysis::ModeAnalysis;
use nyx_engine::facade::physics::SolveResult;
use nyx_engine::facade::{Display, SimplexMesh, ViewMesh};
use nyx_engine::MeshedModel;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
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

/// The mesh on show and what we know about its origin.
#[derive(Default)]
pub struct Doc {
    /// the engine's meshed model: the exchange mesh, the names of regions and boundaries
    pub model: Option<Arc<MeshedModel>>,
    /// the mesh in the numbering the plots read
    pub mesh: Option<Arc<SimplexMesh>>,
    pub view: Option<Arc<ViewMesh>>,
    pub source: String,
    pub mesh_error: Option<String>,
}

impl Doc {
    pub fn dim(&self) -> usize {
        self.mesh.as_ref().map(|m| m.dim).or_else(|| self.view.as_ref().map(|v| v.dim)).unwrap_or(2)
    }

    pub fn set_model(&mut self, model: Arc<MeshedModel>, source: &str) {
        let mesh = model.solver_mesh();
        self.view = Some(Arc::new(ViewMesh::from_simplex(&mesh)));
        self.mesh = Some(mesh);
        self.model = Some(model);
        self.source = source.into();
        self.mesh_error = None;
    }
}

/// A field in the flat form the window reads: scalars keep their name, a vector `E` is `|E|`, `Ex`, `Ez`.
#[derive(Clone, Debug)]
pub struct ViewField {
    pub name: String,
    pub unit: String,
    pub values: Vec<f64>,
}

/// One solution of a result in the form the window reads.
#[derive(Clone, Debug, Default)]
pub struct ViewSolution {
    pub label: String,
    pub fields: Vec<ViewField>,
    /// displacement per dof, `dim` values each (solid mechanics): the deformed shape is drawn from it
    pub disp: Option<Vec<f64>>,
    /// time (s), frequency (Hz) or mode number of this solution
    pub value: Option<f64>,
    /// the effective index (re, im) of a mode study
    pub neff: Option<(f64, f64)>,
}

fn view_solution(s: &nyx_engine::facade::physics::Solution, study: &str, dim: usize) -> ViewSolution {
    let fields = s.flat_fields().into_iter().map(|(name, unit, values)| ViewField { name, unit, values }).collect();
    // the displacement of a real study is the vector field `u`; a complex-frequency mode has the real parts of its components
    let disp = if let Some(u) = s.field("u") {
        Some(u.values.clone())
    } else {
        let axes: &[&str] = if dim == 3 { &["x", "y", "z"] } else { &["x", "z"] };
        let comps: Option<Vec<&Vec<f64>>> = axes.iter().map(|ax| s.field(&format!("Re u{ax}")).map(|f| &f.values)).collect();
        comps.and_then(|comps| {
            let n = comps.first()?.len();
            let mut d = vec![0.0; n * dim];
            for (k, c) in comps.iter().enumerate() {
                for i in 0..n {
                    d[i * dim + k] = c[i];
                }
            }
            Some(d)
        })
    };
    ViewSolution {
        label: s.label.clone(),
        fields,
        disp,
        value: s.value,
        neff: if study == "modes" { s.eigenvalue.map(|e| (e[0], e[1])) } else { None },
    }
}

/// One result on display, with the plot settings.
pub struct ResultView {
    pub name: String,
    /// the code of the physics that made it (`ems`), absent for a file that was opened
    pub physics: Option<String>,
    pub result: Arc<SolveResult>,
    pub sols: Vec<ViewSolution>,
    pub display: Display,
    pub sol: usize,
    pub field: usize,
    pub manual_range: Option<(f32, f32)>,
    pub warp: f32,
    pub auto_warp: f32,
    pub values_cache: Option<(usize, usize, Arc<Vec<f32>>)>,
    pub disp_cache: Option<(usize, Arc<Vec<f32>>)>,
    /// polarization of the modes (Stokes parameters, Poincare sphere); an analysis is computed from a result, it is not part of it
    pub analysis: Option<ModeAnalysis>,
}

impl ResultView {
    pub fn new(name: String, physics: Option<String>, result: Arc<SolveResult>, mesh: &SimplexMesh) -> Self {
        let dim = result.dim as usize;
        let display = Display::new(mesh, result.order as usize);
        let sols: Vec<ViewSolution> = result.solutions.iter().map(|s| view_solution(s, &result.study, dim)).collect();
        // deformation scale such that the largest displacement is 10 % of the model size
        let mut auto_warp = 0.0;
        if let Some(d) = sols.iter().filter_map(|s| s.disp.as_ref()).next() {
            let maxd = d.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            let (lo, hi) = bbox_of(&display);
            let size = (0..dim).map(|k| (hi[k] - lo[k]) as f64).fold(0.0, f64::max);
            if maxd > 0.0 && size > 0.0 {
                auto_warp = ((if sols.len() > 1 { 0.2 } else { 0.1 }) * size / maxd) as f32;
            }
        }
        let field = sols.first().and_then(|s| s.fields.iter().position(|f| f.name == "|u|" || f.name == "T" || f.name == "V" || f.name == "|E|")).unwrap_or(0);
        let has_disp = sols.iter().any(|s| s.disp.is_some());
        ResultView { name, physics, result, sols, display, sol: 0, field, manual_range: None, warp: if has_disp { auto_warp } else { 0.0 }, auto_warp, values_cache: None, disp_cache: None, analysis: None }
    }

    pub fn solution(&self) -> &ViewSolution {
        &self.sols[self.sol.min(self.sols.len() - 1)]
    }

    pub fn field_values(&mut self) -> Arc<Vec<f32>> {
        let sol = self.sol.min(self.sols.len() - 1);
        let nf = self.sols[sol].fields.len();
        let field = self.field.min(nf.saturating_sub(1));
        if let Some((s, f, v)) = &self.values_cache {
            if *s == sol && *f == field {
                return v.clone();
            }
        }
        let v = Arc::new(self.display.pick(&self.sols[sol].fields[field].values));
        self.values_cache = Some((sol, field, v.clone()));
        v
    }

    /// Displacement of each display point (3 values each) of the current solution.
    pub fn displacement(&mut self) -> Option<Arc<Vec<f32>>> {
        let sol = self.sol.min(self.sols.len() - 1);
        if let Some((s, v)) = &self.disp_cache {
            if *s == sol {
                return Some(v.clone());
            }
        }
        let d = self.sols[sol].disp.as_ref()?;
        let dim = self.result.dim as usize;
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
