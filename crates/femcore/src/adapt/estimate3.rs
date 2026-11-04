//! Error estimation and the solve-estimate-adapt loop for tetrahedral meshes (the 3D counterpart of `estimate`, without the 2D geometry
//! pipeline: the start mesh is any tetrahedral mesh, adaptation is `adapt::drive3`).
//!
//! The estimator is the recovered-Hessian one: for each vertex and each domain label meeting there, a quadratic is fitted by least squares
//! to the nodal values on the patch of vertices of that label (never across a material interface), and the absolute Hessians are
//! combined into a metric for about `target_cells` tetrahedra (`metric::lp_metrics`). `Estimator` is the one of `estimate`
//! (fields by name, `combine` for the components of a vector field).

use crate::adapt::estimate::Estimator;
use crate::mesh::SimplexMesh;
use crate::adapt::drive3::adapt3;
use crate::adapt::field3::BackgroundMetric3;
use crate::adapt::metric::{complexity_for_elements, intersect, isotropic, lp_metrics, map_eig, symmetrize, LpOptions, Sym};
use crate::adapt::ops3::AdaptOptions3;
use crate::adapt::topo3::TopoOptions3;
use crate::adapt::work::{simplex_from_adapt, simplex_to_adapt};
use crate::adapt::PlotFields;
use crate::problem::SolveResult;
use std::collections::HashMap;
use std::time::Instant;

/// Recovered Hessian of a P1 field at every vertex of a tetrahedral mesh, intersected over the domains meeting at the vertex.
pub fn recovered_hessians(mesh: &SimplexMesh, u: &[f64]) -> Result<Vec<Sym<3>>, String> {
    if mesh.dim != 3 {
        return Err("this recovery works on tetrahedral meshes".into());
    }
    let n = mesh.n_points();
    if u.len() < n {
        return Err(format!("a field of {} values on a mesh of {n} vertices", u.len()));
    }
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
    let pt = |i: u32| [mesh.points[3 * i as usize], mesh.points[3 * i as usize + 1], mesh.points[3 * i as usize + 2]];
    let mut out: Vec<Option<Sym<3>>> = vec![None; n];
    for (v, l) in keys {
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
            if ring >= 1 && patch.len() >= 20 {
                break;
            }
        }
        let h = fit_hessian(&patch, v, &pt, u).unwrap_or([[0.0; 3]; 3]);
        let abs = map_eig(&symmetrize(&h), f64::abs);
        out[v as usize] = Some(match out[v as usize] {
            None => abs,
            Some(prev) => intersect(&prev, &abs),
        });
    }
    Ok(out.into_iter().map(|m| m.unwrap_or([[0.0; 3]; 3])).collect())
}

/// Least-squares quadratic through the values at `patch` around `centre`; its Hessian, or `None` if the points do not determine one.
fn fit_hessian(patch: &[u32], centre: u32, pt: &dyn Fn(u32) -> [f64; 3], u: &[f64]) -> Option<[[f64; 3]; 3]> {
    const N: usize = 10;
    if patch.len() < N {
        return None;
    }
    let c = pt(centre);
    let scale = patch.iter().map(|&w| (0..3).map(|k| (pt(w)[k] - c[k]).powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max);
    if scale <= 0.0 {
        return None;
    }
    let mut a = [[0.0f64; N + 1]; N];
    for &w in patch {
        let p = pt(w);
        let (x, y, z) = ((p[0] - c[0]) / scale, (p[1] - c[1]) / scale, (p[2] - c[2]) / scale);
        let phi = [1.0, x, y, z, 0.5 * x * x, 0.5 * y * y, 0.5 * z * z, x * y, x * z, y * z];
        for i in 0..N {
            for j in 0..N {
                a[i][j] += phi[i] * phi[j];
            }
            a[i][N] += phi[i] * u[w as usize];
        }
    }
    for col in 0..N {
        let piv = (col..N).max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap())?;
        if a[piv][col].abs() < 1e-11 {
            return None;
        }
        a.swap(col, piv);
        for r in col + 1..N {
            let f = a[r][col] / a[col][col];
            for k in col..=N {
                a[r][k] -= f * a[col][k];
            }
        }
    }
    let mut x = [0.0f64; N];
    for i in (0..N).rev() {
        let mut sum = a[i][N];
        for j in i + 1..N {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    let s2 = scale * scale;
    Some([[x[4] / s2, x[7] / s2, x[8] / s2], [x[7] / s2, x[5] / s2, x[9] / s2], [x[8] / s2, x[9] / s2, x[6] / s2]])
}

/// The metric the solution asks for at every vertex of `mesh`, for about `target_cells` tetrahedra, with sizes between `hmin_fraction * h`
/// and `h`. Fields that are constant or linear (zero Hessian) ask for nothing. `floor` (a metric per vertex, for example from
/// `physics_size::metrics_d`) is intersected when given.
pub fn estimator_metrics<R: PlotFields>(mesh: &SimplexMesh, result: &R, est: &Estimator, h: f64, target_cells: f64, floor: Option<&[Sym<3>]>) -> Result<Vec<Sym<3>>, String> {
    let n = mesh.n_points();
    let adapt_mesh = simplex_to_adapt::<3>(mesh)?;
    let vol = adapt_mesh.vertex_volumes();
    let lp = LpOptions {
        p: est.p,
        target_complexity: complexity_for_elements(target_cells, 3),
        hmin: h * est.hmin_fraction,
        hmax: h,
        max_aspect: est.max_aspect,
        ..LpOptions::default()
    };
    let mut asked: Vec<Sym<3>> = vec![isotropic::<3>(h); n];
    let mut used = 0;
    for sol in result.plot_fields() {
        let names: Vec<&str> = if est.fields.is_empty() { sol.first().map(|f| f.0.as_str()).into_iter().collect() } else { est.fields.iter().map(String::as_str).collect() };
        let mut sum: Vec<Sym<3>> = vec![[[0.0; 3]; 3]; n];
        let mut any = false;
        for name in names {
            let Some(f) = sol.iter().find(|f| f.0 == name) else { continue };
            let hess = recovered_hessians(mesh, &f.1)?;
            let norm: f64 = hess.iter().map(|m| m[0][0].abs().max(m[1][1].abs()).max(m[2][2].abs())).fold(0.0, f64::max);
            if !(norm > 1e-300 && norm.is_finite()) {
                continue;
            }
            if est.combine {
                for (acc, m) in sum.iter_mut().zip(&hess) {
                    let a = map_eig(m, f64::abs);
                    for i in 0..3 {
                        for j in 0..3 {
                            acc[i][j] += a[i][j];
                        }
                    }
                }
                any = true;
            } else {
                for (a, b) in asked.iter_mut().zip(&lp_metrics::<3>(&hess, &vol, &lp).metrics) {
                    *a = intersect(a, b);
                }
                used += 1;
            }
        }
        if est.combine && any {
            for (a, b) in asked.iter_mut().zip(&lp_metrics::<3>(&sum, &vol, &lp).metrics) {
                *a = intersect(a, b);
            }
            used += 1;
        }
    }
    if used == 0 {
        return Err("no result field with a usable Hessian: name the fields of the estimator".into());
    }
    if let Some(fl) = floor {
        for (a, f) in asked.iter_mut().zip(fl) {
            *a = intersect(a, f);
        }
    }
    Ok(asked)
}

/// Settings of the loop.
#[derive(Clone, Debug)]
pub struct LoopSettings3 {
    /// Global (largest) element size.
    pub h: f64,
    pub estimator: Estimator,
    /// Slope of the size grading of each adaptation step (see `metric::grade`).
    pub grade: f64,
    pub adapt: AdaptOptions3,
    /// At most this many adaptation cycles after the start mesh.
    pub max_cycles: usize,
    /// Stop when the quantity changes by less than this fraction between two cycles.
    pub rel_tol: f64,
    /// Target cell count of cycle `k` is `growth^k` times that of the start mesh.
    pub growth: f64,
}

impl LoopSettings3 {
    pub fn new(h: f64) -> Self {
        LoopSettings3 { h, estimator: Estimator::default(), grade: 0.4, adapt: AdaptOptions3::default(), max_cycles: 4, rel_tol: 1e-4, growth: 2.0 }
    }
}

/// One cycle of the loop.
#[derive(Clone, Debug)]
pub struct Cycle3 {
    pub cells: usize,
    pub quantity: f64,
    pub change: Option<f64>,
    pub mesh_seconds: f64,
    pub solve_seconds: f64,
}

pub struct LoopResult3<R = SolveResult> {
    pub mesh: SimplexMesh,
    pub result: R,
    pub cycles: Vec<Cycle3>,
    pub converged: bool,
}

/// Solve on `start`, adapt to the solution, solve again, until the quantity of interest settles. `solve` gets a mesh and returns the
/// result and the number being tracked. `floor` gives, for a mesh, a metric per vertex that adaptation must always satisfy (physics or
/// geometry sizes); `None` for no floor.
pub fn solve_adaptively3<R: PlotFields>(
    start: &SimplexMesh,
    s: &LoopSettings3,
    floor: Option<&dyn Fn(&SimplexMesh) -> Result<Vec<Sym<3>>, String>>,
    solve: &mut dyn FnMut(&SimplexMesh) -> Result<(R, f64), String>,
) -> Result<LoopResult3<R>, String> {
    let mut mesh = start.clone();
    let cells0 = mesh.n_cells() as f64;
    let mut cycles: Vec<Cycle3> = vec![];
    let mut mesh_seconds = 0.0;
    let mut calib = 1.0;
    let mut k = 0;
    loop {
        let t = Instant::now();
        let (result, q) = solve(&mesh)?;
        let change = cycles.last().map(|c: &Cycle3| (q - c.quantity).abs() / q.abs().max(1e-300));
        cycles.push(Cycle3 { cells: mesh.n_cells(), quantity: q, change, mesh_seconds, solve_seconds: t.elapsed().as_secs_f64() });
        let done = change.is_some_and(|c| c < s.rel_tol);
        if done || k >= s.max_cycles {
            return Ok(LoopResult3 { mesh, result, cycles, converged: done });
        }
        k += 1;
        let requested = cells0 * s.growth.powi(k as i32) / calib;
        let fl = match floor {
            Some(f) => Some(f(&mesh)?),
            None => None,
        };
        let metrics = estimator_metrics(&mesh, &result, &s.estimator, s.h, requested, fl.as_deref())?;
        let t = Instant::now();
        let old = simplex_to_adapt::<3>(&mesh)?;
        let field = BackgroundMetric3::new_graded(&old, &metrics, s.grade);
        let out = adapt3(&old, &field, &TopoOptions3::default(), &s.adapt)?;
        mesh = simplex_from_adapt(&out.mesh);
        mesh_seconds = t.elapsed().as_secs_f64();
        calib = mesh.n_cells() as f64 / requested;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::fixtures::box_tet;

    #[test]
    fn the_recovery_is_exact_for_quadratics_on_an_unstructured_mesh() {
        // adapt a box to a non-uniform size so the mesh has no structure, then recover the Hessian of a quadratic
        let start = box_tet([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]], [3, 3, 3]);
        let field = |x: &[f64; 3]| isotropic::<3>(0.12 + 0.25 * x[0]);
        let out = adapt3(&start, &field, &TopoOptions3::default(), &AdaptOptions3::default()).unwrap();
        let mesh = simplex_from_adapt(&out.mesh);
        let u: Vec<f64> = (0..mesh.n_points())
            .map(|i| {
                let (x, y, z) = (mesh.points[3 * i], mesh.points[3 * i + 1], mesh.points[3 * i + 2]);
                0.5 * (2.0 * x * x + 6.0 * y * y - 4.0 * z * z) + x * y - 3.0 * x * z + 2.0 * y * z + x - y
            })
            .collect();
        let h = recovered_hessians(&mesh, &u).unwrap();
        // the Hessian is [[2,1,-3],[1,6,2],[-3,2,-4]]; the recovery returns |H| (eigenvalues made positive), so compare that
        let exact = map_eig(&[[2.0, 1.0, -3.0], [1.0, 6.0, 2.0], [-3.0, 2.0, -4.0]], f64::abs);
        let worst = h.iter().map(|m| (0..3).flat_map(|i| (0..3).map(move |j| (m[i][j] - exact[i][j]).abs())).fold(0.0, f64::max)).fold(0.0, f64::max);
        assert!(worst < 1e-8, "worst deviation {worst}");
    }
}
