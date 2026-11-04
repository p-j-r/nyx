//! The adaptation driver: iterate split, collapse, flip and smooth passes until the mesh conforms
//! to the metric field.

use super::field::MetricField;
use super::geom::{Geometry2, PolylineGeometry};
use super::metric::Sym;
use super::ops::{all_edges, collapse_pass, elen, flip_pass, smooth_pass, split_pass, tri_quality, unresolved_frozen};
use super::topo::{TopoOptions, Work2};
use crate::adapt::work::SimplexMesh;
use crate::adapt::timer::Phases;

#[derive(Clone, Debug)]
pub struct AdaptOptions {
    /// Edges longer than this (in the metric) are split. Default `sqrt 2`.
    pub lmax: f64,
    /// Edges shorter than this are collapsed. Default `1 / sqrt 2`; with these two values a split
    /// never produces an edge that the next collapse pass removes.
    pub lmin: f64,
    pub max_iter: usize,
    /// Smoothing sweeps per iteration.
    pub smooth_sweeps: usize,
    /// Relaxation of a smoothing step.
    pub smooth_relax: f64,
    /// A collapse may not create a triangle below `min(collapse_qmin, quality it replaces)`.
    pub collapse_qmin: f64,
    /// During the first `sample_iters` iterations edges are also measured through the field value at their
    /// midpoint, so that features narrower than the starting mesh are found. Later iterations
    /// use the vertex metrics only (cheaper).
    pub sample_iters: usize,
    /// Runaway guard for metric fields that cannot be satisfied: after the first iteration (which
    /// overshoots legitimately) the adaptation stops if an iteration ends with more triangles than
    /// this, and a split cascade is stopped at twenty times it. `None` for no cap.
    pub max_cells: Option<usize>,
    /// Stop when an iteration changes the vertex count by splits plus collapses of less than this
    /// fraction of the triangles (and, with it, any residual ping-pong is ignored).
    pub converge_frac: f64,
    pub topo: TopoOptions,
}

impl Default for AdaptOptions {
    fn default() -> Self {
        AdaptOptions {
            lmax: std::f64::consts::SQRT_2,
            lmin: std::f64::consts::FRAC_1_SQRT_2,
            max_iter: 40,
            smooth_sweeps: 3,
            smooth_relax: 0.7,
            collapse_qmin: 0.2,
            sample_iters: 2,
            max_cells: None,
            converge_frac: 0.002,
            topo: TopoOptions::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct IterStats {
    pub iter: usize,
    pub n_vertices: usize,
    pub n_triangles: usize,
    pub splits: usize,
    pub collapses: usize,
    pub flips: usize,
    pub moves: usize,
    /// Fraction of edges with metric length in `[1/sqrt 2, sqrt 2]`.
    pub in_band: f64,
    pub q_min: f64,
    pub q_mean: f64,
}

#[derive(Clone, Debug, Default)]
pub struct AdaptReport {
    pub iters: Vec<IterStats>,
    /// Frozen edges still longer than `lmax` at the end: they cannot be split, so triangles next
    /// to them stay coarse. Zero when the frozen boundaries were resolved by the input mesh.
    pub unresolved_frozen: usize,
    /// True when the iteration stopped because the changes fell below `converge_frac`; false when
    /// `max_iter` or `max_cells` ended it (the mesh is valid but does not fully conform).
    pub converged: bool,
    /// True if the `max_cells` cap was hit.
    pub capped: bool,
    pub phases: Phases,
}

pub struct Adapted {
    pub mesh: SimplexMesh<2>,
    /// Metric at the vertices of `mesh`.
    pub metrics: Vec<Sym<2>>,
    pub report: AdaptReport,
}

fn stats(w: &Work2, iter: usize, ops: [usize; 4]) -> IterStats {
    let edges = all_edges(w);
    let band = edges
        .iter()
        .filter(|e| {
            let l = elen(w, e[0], e[1]);
            (std::f64::consts::FRAC_1_SQRT_2..=std::f64::consts::SQRT_2).contains(&l)
        })
        .count();
    let (mut qmin, mut qsum, mut nt) = (f64::INFINITY, 0.0, 0);
    let mut verts = std::collections::HashSet::new();
    for (t, tr) in w.tris.iter().enumerate() {
        if !w.talive[t] {
            continue;
        }
        let p = |i: usize| w.pts[tr[i] as usize];
        let q =
            tri_quality([p(0), p(1), p(2)], [&w.met[tr[0] as usize], &w.met[tr[1] as usize], &w.met[tr[2] as usize]]);
        qmin = qmin.min(q);
        qsum += q;
        nt += 1;
        verts.extend(tr.iter().cloned());
    }
    IterStats {
        iter,
        n_vertices: verts.len(),
        n_triangles: nt,
        splits: ops[0],
        collapses: ops[1],
        flips: ops[2],
        moves: ops[3],
        in_band: band as f64 / edges.len().max(1) as f64,
        q_min: qmin,
        q_mean: qsum / nt.max(1) as f64,
    }
}

/// Adapt `mesh` to `field` with the polyline geometry of its own boundary and interfaces.
pub fn adapt2(mesh: &SimplexMesh<2>, field: &dyn MetricField<2>, opt: &AdaptOptions) -> Result<Adapted, String> {
    adapt2_with(mesh, field, opt, PolylineGeometry::from_work)
}

/// Adapt with a geometry built from the classified working mesh (for example a
/// [`PolylineGeometry`] with analytic projectors for curved boundaries).
pub fn adapt2_with<G: Geometry2>(
    mesh: &SimplexMesh<2>,
    field: &dyn MetricField<2>,
    opt: &AdaptOptions,
    make_geom: impl FnOnce(&Work2) -> G,
) -> Result<Adapted, String> {
    let mut report = AdaptReport::default();
    let mut w = report.phases.time("classify", || Work2::from_mesh(mesh, field, &opt.topo))?;
    let geom = make_geom(&w);
    for iter in 0..opt.max_iter {
        let (mut ns, mut nc, mut nf, mut nm) = (0, 0, 0, 0);
        // flips between split passes keep the cascade of bisections from degrading slivers
        for _ in 0..30 {
            // loose safety limit inside the split cascade (the first iteration overshoots legitimately)
            if opt.max_cells.is_some_and(|m| w.n_alive() >= m.saturating_mul(20)) {
                report.capped = true;
                break;
            }
            let n = report.phases.time("split", || split_pass(&mut w, field, &geom, opt, iter < opt.sample_iters));
            ns += n;
            if n == 0 {
                break;
            }
            nf += report.phases.time("flip", || flip_pass(&mut w));
        }
        for _ in 0..10 {
            let n = report.phases.time("collapse", || collapse_pass(&mut w, opt));
            nc += n;
            if n == 0 {
                break;
            }
        }
        for _ in 0..10 {
            let n = report.phases.time("flip", || flip_pass(&mut w));
            nf += n;
            if n == 0 {
                break;
            }
        }
        for _ in 0..opt.smooth_sweeps {
            nm += report.phases.time("smooth", || smooth_pass(&mut w, field, &geom, opt));
        }
        let st = stats(&w, iter, [ns, nc, nf, nm]);
        let done = (ns + nc) as f64 <= opt.converge_frac * st.n_triangles as f64;
        if std::env::var_os("FEMMESH_VERBOSE").is_some() {
            eprintln!(
                "  iter {:2}: {:7} tris +{} -{} flips {} moves {} in-band {:.3} qmin {:.3}",
                st.iter, st.n_triangles, st.splits, st.collapses, st.flips, st.moves, st.in_band, st.q_min
            );
        }
        report.iters.push(st);
        if done {
            report.converged = true;
            break;
        }
        if report.capped
            || (iter >= 1 && opt.max_cells.is_some_and(|m| report.iters.last().is_some_and(|s| s.n_triangles > m)))
        {
            report.capped = true;
            break;
        }
    }
    for _ in 0..opt.smooth_sweeps {
        report.phases.time("flip", || flip_pass(&mut w));
        report.phases.time("smooth", || smooth_pass(&mut w, field, &geom, opt));
    }
    report.unresolved_frozen = unresolved_frozen(&w, opt.lmax);
    w.check().map_err(|e| format!("internal error after adaptation: {e}"))?;
    let (out, metrics) = w.to_mesh();
    out.check().map_err(|e| format!("adapted mesh invalid: {e}"))?;
    Ok(Adapted { mesh: out, metrics, report })
}
