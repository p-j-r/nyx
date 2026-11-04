//! Error estimation and the solve-estimate-adapt loop for the meshes of `meshing`.
//!
//! The estimator is the standard recovered-Hessian one: the interpolation error of a P1 field is governed by its second derivatives, so
//! the metric that equidistributes it is built from a Hessian recovered from the nodal values (`adapt::metric::lp_metrics`).
//! Two points are specific to nyx:
//!
//! * **Recovery never crosses a material interface.** The normal component of a field jumps there, and a recovery averaged over both sides
//!   smears the jump over the neighbouring cells and asks for refinement everywhere near the interface. The gradient and then the Hessian
//!   are recovered separately for each domain label at a vertex and the resulting metrics are intersected.
//! * **The geometry keeps its say.** The estimator metric is intersected with `meshing::geometry_floor` (`Size` rules, absorbing layers,
//!   curvature), so adaptation to the solution can never drop the resolution the model asked for.
//!
//! `solve_adaptively` runs the loop: solve on the simple mesh, build a metric from the result, adapt, solve again, until the quantity of
//! interest (n_eff, a frequency, a peak value) changes by less than a tolerance between two cycles.

use crate::meshing::geomdrive::{evaluate, mesh_geom_2d, GeomMeshOpts, MeshedGeom};
use crate::mesh::SimplexMesh;
use crate::adapt::metric::{complexity_for_elements, intersect, lp_metrics, map_eig, symmetrize, LpOptions, Sym};
use crate::mesh::interop::simplex_to_adapt;
use crate::adapt::field::BackgroundMetric2;
use crate::meshing::{adapt_to_field, adapt_to_metrics, geometry_floor, Adaptive};
use crate::problem::SolveResult;
use geom::model::FeatureList;
use std::collections::HashMap;
use std::time::Instant;

/// What the metric is built from and how.
#[derive(Clone, Debug)]
pub struct Estimator {
    /// Names of the result fields whose Hessians drive the mesh (`"V"`, `"|E|"`, ...). Empty: the first field of every solution, which by
    /// convention is the primary one. Every solution (every mode) contributes.
    pub fields: Vec<String>,
    /// Norm of the interpolation error that is minimised (2 is the usual choice).
    pub p: f64,
    /// Smallest element size, as a fraction of the global size.
    pub hmin_fraction: f64,
    /// Largest ratio of the long to the short side of an element.
    pub max_aspect: f64,
    /// Treat the named fields of one solution as the components of one vector field: their absolute Hessians are added at each vertex and a
    /// single metric is made from the sum, so the cell target is met once. Without it every field gets its own metric at the full target
    /// and the metrics are intersected, which asks for more cells than the target (each extra component adds its own refinement).
    pub combine: bool,
}

impl Default for Estimator {
    fn default() -> Self {
        Estimator { fields: vec![], p: 2.0, hmin_fraction: 0.005, max_aspect: 20.0, combine: false }
    }
}

impl Estimator {
    /// For the mode fields of electromagnetic studies: the signed Cartesian components of E as one vector field. The magnitude `|E|` is the
    /// default first field, but it has a kink wherever the field changes sign and its recovered Hessian there is noise; the components
    /// are smooth.
    pub fn e_field() -> Self {
        Estimator { fields: vec!["Ex".into(), "Ey".into(), "Ez".into()], combine: true, ..Estimator::default() }
    }
}

/// Recovered Hessian of a P1 field at every vertex of a triangle mesh, intersected over the domains meeting at the vertex.
///
/// For each vertex and each domain label meeting there, a quadratic is fitted by least squares to the nodal values on the patch of
/// vertices of that label within two cells of the vertex (three or four rings when the patch is too small). The fit is exact for any
/// quadratic field on any mesh, which averaging of cell gradients is not on an unstructured mesh. A patch never includes a cell of another
/// label, so a vertex on an interface sees each side separately. The absolute value of each Hessian (eigenvalues made positive) is
/// taken and the sides are intersected.
pub fn recovered_hessians(mesh: &SimplexMesh, u: &[f64]) -> Result<Vec<Sym<2>>, String> {
    if mesh.dim != 2 {
        return Err("Hessian recovery works on 2D meshes".into());
    }
    let n = mesh.n_points();
    if u.len() < n {
        return Err(format!("a field of {} values on a mesh of {n} vertices", u.len()));
    }
    // neighbours of every (vertex, label) through cells of that label
    let mut nbr: HashMap<(u32, i32), Vec<u32>> = HashMap::new();
    for c in 0..mesh.n_cells() {
        let l = mesh.cell_labels[c];
        for &v in mesh.cell(c) {
            let e = nbr.entry((v, l)).or_default();
            for &w in mesh.cell(c) {
                if w != v && !e.contains(&w) {
                    e.push(w);
                }
            }
        }
    }
    let mut keys: Vec<(u32, i32)> = nbr.keys().copied().collect();
    keys.sort_unstable();
    let pt = |i: u32| [mesh.points[2 * i as usize], mesh.points[2 * i as usize + 1]];
    let mut out: Vec<Option<Sym<2>>> = vec![None; n];
    for (v, l) in keys {
        // grow the patch ring by ring until it has enough points for six unknowns
        let mut patch: Vec<u32> = vec![v];
        let mut frontier = vec![v];
        for ring in 0..4 {
            let mut next = vec![];
            for &w in &frontier {
                for &x in &nbr[&(w, l)] {
                    if !patch.contains(&x) {
                        patch.push(x);
                        next.push(x);
                    }
                }
            }
            frontier = next;
            if ring >= 1 && patch.len() >= 9 {
                break;
            }
        }
        let h = fit_hessian(&patch, v, &pt, u).unwrap_or([[0.0; 2]; 2]);
        let abs = map_eig(&symmetrize(&h), f64::abs);
        out[v as usize] = Some(match out[v as usize] {
            None => abs,
            Some(prev) => intersect(&prev, &abs),
        });
    }
    Ok(out.into_iter().map(|m| m.unwrap_or([[0.0; 2]; 2])).collect())
}

/// Least-squares quadratic through the values at `patch` around `centre`; its Hessian, or `None` if the points do not determine a quadratic.
fn fit_hessian(patch: &[u32], centre: u32, pt: &dyn Fn(u32) -> [f64; 2], u: &[f64]) -> Option<[[f64; 2]; 2]> {
    if patch.len() < 6 {
        return None;
    }
    let c = pt(centre);
    let scale = patch.iter().map(|&w| (pt(w)[0] - c[0]).hypot(pt(w)[1] - c[1])).fold(0.0, f64::max);
    if scale <= 0.0 {
        return None;
    }
    let mut a = [[0.0f64; 7]; 6];
    for &w in patch {
        let p = pt(w);
        let (x, z) = ((p[0] - c[0]) / scale, (p[1] - c[1]) / scale);
        let phi = [1.0, x, z, 0.5 * x * x, x * z, 0.5 * z * z];
        for i in 0..6 {
            for j in 0..6 {
                a[i][j] += phi[i] * phi[j];
            }
            a[i][6] += phi[i] * u[w as usize];
        }
    }
    // Gauss elimination with partial pivoting
    for col in 0..6 {
        let piv = (col..6).max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())?;
        if a[piv][col].abs() < 1e-11 {
            return None;
        }
        a.swap(col, piv);
        for r in col + 1..6 {
            let f = a[r][col] / a[col][col];
            for k in col..7 {
                a[r][k] -= f * a[col][k];
            }
        }
    }
    let mut x = [0.0f64; 6];
    for i in (0..6).rev() {
        let mut sum = a[i][6];
        for j in i + 1..6 {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    let s2 = scale * scale;
    Some([[x[3] / s2, x[4] / s2], [x[4] / s2, x[5] / s2]])
}

/// The metric the solution asks for at every vertex of the mesh of `g`, for about `target_cells` triangles, intersected with the floor of the
/// geometry. Fields that are constant or linear (zero Hessian) ask for nothing.
pub fn estimator_metrics(g: &MeshedGeom, result: &SolveResult, est: &Estimator, opts: &GeomMeshOpts, target_cells: f64) -> Result<Vec<Sym<2>>, String> {
    let mesh = &g.mesh;
    let n = mesh.n_points();
    let adapt_mesh = simplex_to_adapt::<2>(mesh)?;
    let vol = adapt_mesh.vertex_volumes();
    let lp = LpOptions {
        p: est.p,
        target_complexity: complexity_for_elements(target_cells, 2),
        hmin: opts.h * est.hmin_fraction,
        hmax: opts.h,
        max_aspect: est.max_aspect,
        ..LpOptions::default()
    };
    let mut asked: Vec<Sym<2>> = vec![[[1.0 / (opts.h * opts.h), 0.0], [0.0, 1.0 / (opts.h * opts.h)]]; n];
    let mut used = 0;
    // `combine`: the sum over the components of |H| at every vertex, per solution
    let mut sum: Vec<Sym<2>> = vec![[[0.0; 2]; 2]; n];
    for sol in &result.solutions {
        let mut any = false;
        let names: Vec<&str> = if est.fields.is_empty() { sol.fields.first().map(|f| f.name.as_str()).into_iter().collect() } else { est.fields.iter().map(String::as_str).collect() };
        for name in names {
            let Some(f) = sol.field(name) else { continue };
            let hess = recovered_hessians(mesh, &f.values)?;
            let norm: f64 = hess.iter().map(|h| h[0][0].abs().max(h[1][1].abs())).fold(0.0, f64::max);
            if !(norm > 1e-300 && norm.is_finite()) {
                continue;
            }
            if est.combine {
                for (acc, h) in sum.iter_mut().zip(&hess) {
                    let a = map_eig(h, f64::abs);
                    for i in 0..2 {
                        for j in 0..2 {
                            acc[i][j] += a[i][j];
                        }
                    }
                }
                any = true;
            } else {
                let m = lp_metrics::<2>(&hess, &vol, &lp);
                for (a, b) in asked.iter_mut().zip(&m.metrics) {
                    *a = intersect(a, b);
                }
                used += 1;
            }
        }
        if est.combine && any {
            used += 1;
        }
    }
    if est.combine && used > 0 {
        let m = lp_metrics::<2>(&sum, &vol, &lp);
        for (a, b) in asked.iter_mut().zip(&m.metrics) {
            *a = intersect(a, b);
        }
    }
    if used == 0 {
        return Err("no result field with a usable Hessian: name the fields of the estimator".into());
    }
    let floor = geometry_floor(g, opts, &Adaptive::default())?;
    Ok(asked.iter().zip(&floor).map(|(a, f)| intersect(a, f)).collect())
}

/// Settings of the loop.
#[derive(Clone, Debug)]
pub struct LoopSettings {
    /// Sizes of the start mesh (always made by the simple mesher at these sizes).
    pub opts: GeomMeshOpts,
    pub estimator: Estimator,
    /// Grading and the cell cap of the adaptation steps (`scale` is not used: the estimator sets the sizes).
    pub adapt: Adaptive,
    /// At most this many adaptation cycles after the start mesh.
    pub max_cycles: usize,
    /// Stop when the quantity changes by less than this fraction between two cycles.
    pub rel_tol: f64,
    /// Target cell count of cycle `k` is `growth^k` times that of the start mesh.
    pub growth: f64,
    /// The final metric of an earlier run on a similar model (a neighbouring point of a sweep). The start mesh is adapted to it at once, so
    /// the loop begins near the resolution the last point needed instead of climbing to it, and usually needs one or two solves.
    pub warm: Option<WarmStart>,
}

/// The mesh and metric an adaptive run ended with.
#[derive(Clone, Debug)]
pub struct WarmStart {
    pub mesh: SimplexMesh,
    pub metrics: Vec<Sym<2>>,
    /// cells obtained per cell requested by the last adaptation (adapting to a graded metric over-delivers); the next run divides its
    /// requests by it
    pub calibration: f64,
}

impl LoopSettings {
    pub fn new(h: f64) -> Self {
        LoopSettings { opts: GeomMeshOpts::new(h), estimator: Estimator::default(), adapt: Adaptive::default(), max_cycles: 5, rel_tol: 1e-5, growth: 2.0, warm: None }
    }
}

/// One cycle of the loop.
#[derive(Clone, Debug)]
pub struct Cycle {
    pub cells: usize,
    pub quantity: f64,
    /// relative change of the quantity from the cycle before
    pub change: Option<f64>,
    pub mesh_seconds: f64,
    pub solve_seconds: f64,
}

pub struct LoopResult {
    /// the last mesh and its solution
    pub geom: MeshedGeom,
    pub result: SolveResult,
    pub cycles: Vec<Cycle>,
    /// true when the change of the quantity fell below the tolerance
    pub converged: bool,
    /// what a neighbouring run can start from (`LoopSettings::warm`)
    pub warm: WarmStart,
}

/// Solve on the start mesh, then adapt to the solution and solve again until the quantity of interest settles. `solve` meshes nothing: it
/// gets a meshed model and returns the result together with the number being tracked (n_eff of the first mode, a peak, an energy).
pub fn solve_adaptively(
    list: &FeatureList,
    s: &LoopSettings,
    solve: &mut dyn FnMut(&MeshedGeom) -> Result<(SolveResult, f64), String>,
) -> Result<LoopResult, String> {
    let t = Instant::now();
    let mut geom = mesh_geom_2d(evaluate(list)?, &s.opts)?;
    let mut warm_used = false;
    if let Some(w) = &s.warm {
        let field = BackgroundMetric2::new_graded(&simplex_to_adapt::<2>(&w.mesh)?, &w.metrics, s.adapt.grade);
        geom = adapt_to_field(geom, &field, &s.adapt)?.0;
        warm_used = true;
    }
    let mut mesh_seconds = t.elapsed().as_secs_f64();
    let cells0 = geom.mesh.n_cells() as f64;
    let mut cycles: Vec<Cycle> = vec![];
    let mut k = 0;
    let mut calib = s.warm.as_ref().map_or(1.0, |w| w.calibration);
    loop {
        let t = Instant::now();
        let (result, q) = solve(&geom)?;
        let change = cycles.last().map(|c: &Cycle| (q - c.quantity).abs() / q.abs().max(1e-300));
        cycles.push(Cycle { cells: geom.mesh.n_cells(), quantity: q, change, mesh_seconds, solve_seconds: t.elapsed().as_secs_f64() });
        let done = change.is_some_and(|c| c < s.rel_tol);
        let last = done || k >= s.max_cycles;
        // cells wanted next; the over-delivery measured at the last adaptation is divided out
        let wanted = if last { geom.mesh.n_cells() as f64 } else if warm_used { geom.mesh.n_cells() as f64 * s.growth.min(1.5) } else { cells0 * s.growth.powi(k as i32 + 1) };
        let requested = wanted / calib;
        let metrics = estimator_metrics(&geom, &result, &s.estimator, &s.opts, requested)?;
        if last {
            let warm = WarmStart { mesh: geom.mesh.clone(), metrics, calibration: calib };
            return Ok(LoopResult { geom, result, cycles, converged: done, warm });
        }
        k += 1;
        let t = Instant::now();
        let (next, report) = adapt_to_metrics(geom, &metrics, &s.adapt)?;
        mesh_seconds = t.elapsed().as_secs_f64();
        geom = next;
        calib = geom.mesh.n_cells() as f64 / requested;
        if report.capped {
            // at the cap: solve once more on this mesh and stop
            let (result, q) = solve(&geom)?;
            let change = cycles.last().map(|c| (q - c.quantity).abs() / q.abs().max(1e-300));
            cycles.push(Cycle { cells: geom.mesh.n_cells(), quantity: q, change, mesh_seconds, solve_seconds: 0.0 });
            let metrics = estimator_metrics(&geom, &result, &s.estimator, &s.opts, geom.mesh.n_cells() as f64)?;
            let warm = WarmStart { mesh: geom.mesh.clone(), metrics, calibration: calib };
            return Ok(LoopResult { geom, result, cycles, converged: false, warm });
        }
    }
}
