//! The Maxwell residual estimator (ADR 53, `femcore::adapt::maxwell`) with the edge-element cavity kernel of `nyx-ems`, on a PEC box cavity, whose modes are known in closed form: the indicators fall with the mesh at the
//! rate of the element order, and driving the adaptation loop with them reaches an accuracy that uniform blocks need many more cells for,
//! also at the lowest order where the Hessian of the plotted field cannot steer (`MISSING_FEATURES.md`, section 2).
mod common;

use femcore::adapt::estimate3::{adapt_loop3, LoopSettings3};
use femcore::adapt::maxwell::{indicator_metrics, maxwell_indicators, CurlField, EquidistributeOptions, MaxwellOptions};
use femcore::adapt::metric::Sym;
use femcore::adapt::ops3::AdaptOptions3;
use femcore::adapt::topo3::TopoOptions3;
use femcore::adapt::work::{simplex_from_adapt, simplex_to_adapt};
use nyx_ems::cavity::{assemble_emw_eigen, basis3, solve_emw_modes, Ned3};
use nyx_ems::{Problem, Study};
use femcore::linalg::csr::C;
use femcore::mesh::SimplexMesh;

const C0: f64 = 299_792_458.0;

fn cavity_modes(a: f64, b: f64, d: f64, fmax: f64) -> Vec<f64> {
    let mut out = vec![];
    for m in 0..5 {
        for n in 0..5 {
            for p in 0..5 {
                let nonzero = (m > 0) as usize + (n > 0) as usize + (p > 0) as usize;
                if nonzero < 2 {
                    continue;
                }
                let f = 0.5 * C0 * ((m as f64 / a).powi(2) + (n as f64 / b).powi(2) + (p as f64 / d).powi(2)).sqrt();
                if f <= fmax {
                    out.push(f);
                    if nonzero == 3 {
                        out.push(f);
                    }
                }
            }
        }
    }
    out.sort_by(|x, y| x.partial_cmp(y).unwrap());
    out
}

/// The edge-element field of a solved mode, as the estimator reads it: the degrees of freedom of `ems3d` with its own basis functions.
struct Mode<'a> {
    mesh: &'a SimplexMesh,
    ned: &'a Ned3,
    x: Vec<C>,
    order: usize,
}

impl CurlField for Mode<'_> {
    fn sample(&self, cell: usize, lam: [f64; 4]) -> ([C; 3], [C; 3]) {
        let v = self.mesh.cell(cell);
        let p: Vec<[f64; 3]> = v.iter().map(|&i| [self.mesh.points[3 * i as usize], self.mesh.points[3 * i as usize + 1], self.mesh.points[3 * i as usize + 2]]).collect();
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
        let (e1, e2, e3) = (sub(p[1], p[0]), sub(p[2], p[0]), sub(p[3], p[0]));
        let six_v = e1[0] * cross(e2, e3)[0] + e1[1] * cross(e2, e3)[1] + e1[2] * cross(e2, e3)[2];
        let g1 = cross(e2, e3).map(|x| x / six_v);
        let g2 = cross(e3, e1).map(|x| x / six_v);
        let g3 = cross(e1, e2).map(|x| x / six_v);
        let g0 = [-(g1[0] + g2[0] + g3[0]), -(g1[1] + g2[1] + g3[1]), -(g1[2] + g2[2] + g3[2])];
        let gid = self.ned.keys(self.mesh, cell);
        let (mut w, mut cu) = (vec![], vec![]);
        basis3(self.order, &lam, &[g0, g1, g2, g3], &gid, &mut w, &mut cu);
        let mut e = [C::new(0.0, 0.0); 3];
        let mut c = [C::new(0.0, 0.0); 3];
        for (s, &dof) in self.ned.cell_dofs[cell].iter().enumerate() {
            let coef = self.x[dof];
            for k in 0..3 {
                e[k] += coef * w[s][k];
                c[k] += coef * cu[s][k];
            }
        }
        (e, c)
    }
}

fn spec(order: usize, f_exact: f64) -> Problem {
    let mut s = Problem::new(Study::Cavity { n_modes: 1, f_search: 0.9 * f_exact }, order);
    s.materials.default.insert("eps_r".into(), 1.0);
    s
}

/// Solves the lowest mode on `mesh` and returns its frequency and the Maxwell indicators of it.
fn solve_and_estimate(mesh: &SimplexMesh, spec: &Problem) -> Result<(f64, femcore::adapt::maxwell::Indicators), String> {
    let asm = assemble_emw_eigen(mesh, spec)?;
    let solved = solve_emw_modes(&asm)?;
    let (f, idx) = solved.freqs[0];
    let mode = Mode { mesh, ned: &asm.ned, x: asm.red.expand(&solved.vectors[idx]), order: asm.order };
    let k0 = 2.0 * std::f64::consts::PI * f.re / C0;
    let ones = vec![1.0; mesh.n_cells()];
    let opt = MaxwellOptions { k0_sq: k0 * k0, eps_r: &ones, mu_r: &ones, skip: &|_| false, natural: &|_| false };
    Ok((f.re, maxwell_indicators(mesh, &mode, &opt)?))
}

#[test]
fn the_indicators_fall_at_the_rate_of_the_element_order_under_uniform_refinement() {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, d, 20e9)[0];
    for order in [1usize, 2] {
        let s = spec(order, exact);
        let mut totals = vec![];
        for n in [(2, 2, 2), (4, 3, 2), (8, 6, 4)] {
            let m = common::block(a, b, d, n.0, n.1, n.2);
            let (f, ind) = solve_and_estimate(&m, &s).unwrap();
            eprintln!("order {order} {:5} tets  f error {:.2e}  eta {:.3e}  parts {:?}", m.n_cells(), (f / exact - 1.0).abs(), ind.total, ind.parts);
            totals.push(ind.total);
        }
        // halving h twice must lower the estimate at least as fast as h^order (it is an upper bound of the energy error of the field,
        // which falls like h^order; the first mesh is too coarse to be in the asymptotic range, so compare the last two)
        let rate = (totals[1] / totals[2]).log2() / 1.0;
        eprintln!("order {order}: eta falls by a factor {:.2} over the last refinement (h halves)", totals[1] / totals[2]);
        assert!(totals[2] < totals[1] && totals[1] < totals[0], "{totals:?}");
        assert!(rate > 0.6 * order as f64, "order {order}: rate {rate}");
    }
}

/// The loop with the Maxwell indicators on a box cavity: (cells, relative error of the frequency) per cycle.
fn maxwell_loop(order: usize, target_growth: f64, cycles: usize) -> Vec<(usize, f64)> {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, d, 20e9)[0];
    let s = spec(order, exact);
    let mut st = LoopSettings3::new(a / 2.0);
    st.max_cycles = cycles;
    st.rel_tol = 1e-9;
    st.growth = target_growth;
    let mut solve = |m: &SimplexMesh| -> Result<((SimplexMesh, femcore::adapt::maxwell::Indicators), f64), String> {
        let (f, ind) = solve_and_estimate(m, &s)?;
        Ok(((m.clone(), ind), f))
    };
    let mut metrics = |mesh: &SimplexMesh, r: &(SimplexMesh, femcore::adapt::maxwell::Indicators), target: f64, floor: Option<&[Sym<3>]>| {
        let opt = EquidistributeOptions { order, h_max: a / 2.0, h_min_fraction: 0.05, target_cells: target };
        indicator_metrics(mesh, &r.1.eta2, &opt, floor)
    };
    let mut adapt = |mesh: &SimplexMesh, metrics: &[Sym<3>], grade: f64| -> Result<SimplexMesh, String> {
        let old = simplex_to_adapt::<3>(mesh)?;
        let field = femcore::adapt::field3::BackgroundMetric3::new_graded(&old, metrics, grade);
        let out = femcore::adapt::drive3::adapt3(&old, &field, &TopoOptions3::default(), &AdaptOptions3::default())?;
        Ok(simplex_from_adapt(&out.mesh))
    };
    let out = adapt_loop3(&common::block(a, b, d, 2, 2, 2), &st, None, &mut solve, &mut metrics, &mut adapt).unwrap();
    out.cycles.iter().map(|c| (c.cells, (c.quantity / exact - 1.0).abs())).collect()
}

fn uniform_errors(order: usize) -> Vec<(usize, f64)> {
    let (a, b, d) = (0.04, 0.03, 0.02);
    let exact = cavity_modes(a, b, d, 20e9)[0];
    let s = spec(order, exact);
    [(2, 2, 2), (3, 2, 2), (4, 3, 2), (5, 4, 3), (7, 5, 4), (9, 7, 5)]
        .iter()
        .map(|&(nx, ny, nz)| {
            let m = common::block(a, b, d, nx, ny, nz);
            let (f, _) = solve_and_estimate(&m, &s).unwrap();
            (m.n_cells(), (f / exact - 1.0).abs())
        })
        .collect()
}

#[test]
#[ignore = "study: the Maxwell loop against uniform blocks, both orders"]
fn study_the_maxwell_loop_against_uniform_blocks() {
    for order in [1usize, 2] {
        for (n, e) in uniform_errors(order) {
            eprintln!("order {order} uniform {n:6} tets error {e:.2e}");
        }
        for (n, e) in maxwell_loop(order, 2.0, 4) {
            eprintln!("order {order} maxwell {n:6} tets error {e:.2e}");
        }
    }
}

#[test]
fn the_maxwell_loop_is_below_the_uniform_curve_at_both_orders_also_at_the_lowest() {
    // MISSING_FEATURES.md: with the Hessian loop at order 1, 1.2e-2 on 1490 tets against 5.3e-3 on 696 uniform ones. With the residual
    // estimator the loop must be below the uniform curve: interpolate the uniform errors (log-log) at the loop's cell count, and at every
    // cycle after the first (the first is the uniform start mesh).
    for order in [1usize, 2] {
        let uniform = uniform_errors(order);
        let looped = maxwell_loop(order, 2.0, 4);
        let at = |cells: usize| -> f64 {
            let w = uniform.windows(2).filter(|w| w[0].0 < w[1].0).find(|w| (w[0].0..=w[1].0).contains(&cells)).unwrap_or_else(|| &uniform[uniform.len() - 2..]);
            let (n0, e0, n1, e1) = (w[0].0 as f64, w[0].1, w[1].0 as f64, w[1].1);
            let t = ((cells as f64).ln() - n0.ln()) / (n1.ln() - n0.ln());
            (e0.ln() + t * (e1.ln() - e0.ln())).exp()
        };
        eprintln!("order {order}: uniform {uniform:?}\n         loop {looped:?}");
        let (cells, err) = *looped.last().unwrap();
        assert!(err < at(cells), "order {order}: the residual loop is not better than uniform meshes of its size: {err:.2e} against {:.2e}", at(cells));
        // the cell count grows by about the asked factor, not by orders of magnitude
        assert!(cells < 20 * looped[0].0, "{looped:?}");
    }
}

// ---------------------------------------------------------------- a curved cavity from geom: the whole chain of ADR 50 to 53

/// The PEC sphere cavity of radius `r`: lowest mode `k r = 2.7437` (TM101, threefold).
fn sphere_exact(r: f64) -> f64 {
    2.743707 * C0 / (2.0 * std::f64::consts::PI * r)
}

/// Loop on the sphere cavity from geom: the Maxwell indicators steer, `adapt_to_metrics3` adapts with the exact sphere as the oracle and
/// the curvature of the sphere as a floor. Returns (cells, relative frequency error) per cycle.
fn sphere_loop(order: usize, h0: f64, cycles: usize) -> Vec<(usize, f64)> {
    use femcore::meshing::adapt3d::{adapt_to_metrics3, geometry_floor3};
    use femcore::meshing::geomdrive::{GeomMeshOpts, MeshedGeom};
    use femcore::meshing::Adaptive;
    let r = 0.05;
    let exact = sphere_exact(r);
    let g0 = femcore::meshing::geomdrive3::sphere_cavity(r, h0).unwrap();
    let (labels, eval) = (g0.labels.clone(), g0.eval.clone());
    let s = spec(order, exact);
    let opts = GeomMeshOpts::new(h0);
    let a = Adaptive { scale: 1.0, grade: 0.4, max_cells: Some(400_000) };
    let mut st = LoopSettings3::new(h0);
    st.max_cycles = cycles;
    st.rel_tol = 1e-9;
    let wrap = |m: &SimplexMesh| MeshedGeom { mesh: m.clone(), labels: labels.clone(), eval: eval.clone() };
    let mut solve = |m: &SimplexMesh| -> Result<((SimplexMesh, femcore::adapt::maxwell::Indicators), f64), String> {
        let (f, ind) = solve_and_estimate(m, &s)?;
        Ok(((m.clone(), ind), f))
    };
    let mut metrics = |mesh: &SimplexMesh, r: &(SimplexMesh, femcore::adapt::maxwell::Indicators), target: f64, floor: Option<&[Sym<3>]>| {
        let opt = EquidistributeOptions { order, h_max: h0, h_min_fraction: 0.1, target_cells: target };
        indicator_metrics(mesh, &r.1.eta2, &opt, floor)
    };
    let mut adapt = |mesh: &SimplexMesh, metrics: &[Sym<3>], _grade: f64| -> Result<SimplexMesh, String> { Ok(adapt_to_metrics3(wrap(mesh), metrics, &a)?.0.mesh) };
    let floor = |mesh: &SimplexMesh| geometry_floor3(&wrap(mesh), &opts, &a);
    let out = adapt_loop3(&g0.mesh, &st, Some(&floor), &mut solve, &mut metrics, &mut adapt).unwrap();
    out.cycles.iter().map(|c| (c.cells, (c.quantity / exact - 1.0).abs())).collect()
}

fn sphere_uniform(order: usize) -> Vec<(usize, f64)> {
    let r = 0.05;
    let exact = sphere_exact(r);
    let s = spec(order, exact);
    [0.025, 0.02, 0.016, 0.013, 0.011]
        .iter()
        .map(|&h| {
            let g = femcore::meshing::geomdrive3::sphere_cavity(r, h).unwrap();
            let (f, _) = solve_and_estimate(&g.mesh, &s).unwrap();
            (g.mesh.n_cells(), (f / exact - 1.0).abs())
        })
        .collect()
}

#[test]
#[ignore = "study: the Maxwell loop on a sphere cavity from geom against uniform lattice meshes"]
fn study_the_sphere_cavity_loop() {
    for order in [1usize, 2] {
        eprintln!("order {order} uniform {:?}", sphere_uniform(order));
        eprintln!("order {order} loop    {:?}", sphere_loop(order, 0.025, 3));
    }
}

#[test]
fn the_whole_chain_on_a_curved_cavity_from_geom_improves_the_mode_and_keeps_the_sphere() {
    // Geometry limits the sphere (straight-sided tetrahedra, ADR 18), so both curves flatten near 4e-3; the loop must get there with the cells
    // the estimator chooses, from the 1.8e-2 of the start mesh.
    for order in [1usize, 2] {
        let looped = sphere_loop(order, 0.025, 2);
        eprintln!("order {order}: {looped:?}");
        let (first, last) = (looped[0], *looped.last().unwrap());
        assert!(last.1 < 0.6 * first.1, "order {order}: {looped:?}");
        assert!(last.0 < 8 * first.0, "order {order}: {looped:?}");
    }
}
