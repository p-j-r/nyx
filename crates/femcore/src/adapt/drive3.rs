//! The 3D adaptation driver: iterate split, collapse, flip and smooth passes until the mesh conforms to the metric field (the 3D
//! counterpart of `drive`, with the operators of `ops3`).

use super::drive::{AdaptReport, IterStats};
use super::field::MetricField;
use super::geom3::{FacetGeometry, Geometry3};
use super::metric::Sym;
use super::ops3::{all_edges_of, collapse_pass, edge_removal_pass, elen, flip_pass, smooth_pass, split_pass, surface_flip_pass, tet_q, unresolved_frozen, AdaptOptions3};
use super::topo3::{TopoOptions3, Work3};
use crate::adapt::work::SimplexMesh;

pub struct Adapted3 {
    pub mesh: SimplexMesh<3>,
    /// Metric at the vertices of `mesh`.
    pub metrics: Vec<Sym<3>>,
    pub report: AdaptReport,
}

fn stats(w: &Work3, iter: usize, ops: [usize; 4]) -> IterStats {
    let edges = all_edges_of(w);
    let band = edges.iter().filter(|e| (std::f64::consts::FRAC_1_SQRT_2..=std::f64::consts::SQRT_2).contains(&elen(w, e[0], e[1]))).count();
    let (mut qmin, mut qsum, mut nt) = (f64::INFINITY, 0.0, 0);
    let mut verts = std::collections::HashSet::new();
    for (t, tt) in w.tets.iter().enumerate() {
        if w.talive[t] {
            let q = tet_q(w, tt);
            qmin = qmin.min(q);
            qsum += q;
            nt += 1;
            verts.extend(tt.iter().copied());
        }
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

/// Adapt `mesh` to `field` with the facet geometry of its own boundary and interfaces. `AdaptReport::iters[..].n_triangles` counts
/// tetrahedra here.
pub fn adapt3(mesh: &SimplexMesh<3>, field: &dyn MetricField<3>, topo: &TopoOptions3, opt: &AdaptOptions3) -> Result<Adapted3, String> {
    adapt3_with(mesh, field, topo, opt, FacetGeometry::from_work)
}

/// Adapt with a geometry built from the classified working mesh (for example a [`FacetGeometry`] with analytic projectors for curved
/// surfaces).
pub fn adapt3_with<G: Geometry3>(
    mesh: &SimplexMesh<3>,
    field: &dyn MetricField<3>,
    topo: &TopoOptions3,
    opt: &AdaptOptions3,
    make_geom: impl FnOnce(&Work3) -> G,
) -> Result<Adapted3, String> {
    let mut report = AdaptReport::default();
    let mut w = report.phases.time("classify", || Work3::from_mesh(mesh, field, topo))?;
    let geom = make_geom(&w);
    // FEMMESH_CHECK=1 validates the working mesh after every pass and names the first one that breaks it
    let check = std::env::var_os("FEMMESH_CHECK").is_some();
    let verify = |w: &Work3, what: &str, iter: usize| -> Result<(), String> {
        if check {
            w.check().map_err(|e| format!("after {what} in iteration {iter}: {e}"))?;
        }
        Ok(())
    };
    for iter in 0..opt.max_iter {
        let (mut ns, mut nc, mut nf, mut nm) = (0, 0, 0, 0);
        for _ in 0..30 {
            if opt.max_cells.is_some_and(|m| w.n_alive() >= m.saturating_mul(20)) {
                report.capped = true;
                break;
            }
            let n = report.phases.time("split", || split_pass(&mut w, field, &geom, opt, iter < opt.sample_iters));
            verify(&w, "split", iter)?;
            ns += n;
            if n == 0 {
                break;
            }
            nf += report.phases.time("flip", || flip_pass(&mut w));
            verify(&w, "flip after split", iter)?;
        }
        for _ in 0..10 {
            let n = report.phases.time("collapse", || collapse_pass(&mut w, opt));
            verify(&w, "collapse", iter)?;
            nc += n;
            if n == 0 {
                break;
            }
        }
        for _ in 0..10 {
            let mut n = report.phases.time("flip", || flip_pass(&mut w));
            verify(&w, "flip", iter)?;
            if opt.edge_removal {
                n += report.phases.time("edge removal", || edge_removal_pass(&mut w));
                verify(&w, "edge removal", iter)?;
                n += report.phases.time("surface flip", || surface_flip_pass(&mut w));
                verify(&w, "surface flip", iter)?;
            }
            nf += n;
            if n == 0 {
                break;
            }
        }
        for _ in 0..opt.smooth_sweeps {
            nm += report.phases.time("smooth", || smooth_pass(&mut w, field, &geom, opt));
            verify(&w, "smooth", iter)?;
        }
        let st = stats(&w, iter, [ns, nc, nf, nm]);
        let done = (ns + nc) as f64 <= opt.converge_frac * st.n_triangles as f64;
        if std::env::var_os("FEMMESH_VERBOSE").is_some() {
            eprintln!("  iter {:2}: {:7} tets +{} -{} flips {} moves {} in-band {:.3} qmin {:.3}", st.iter, st.n_triangles, st.splits, st.collapses, st.flips, st.moves, st.in_band, st.q_min);
        }
        report.iters.push(st);
        if done {
            report.converged = true;
            break;
        }
        if report.capped || (iter >= 1 && opt.max_cells.is_some_and(|m| report.iters.last().is_some_and(|s| s.n_triangles > m))) {
            report.capped = true;
            break;
        }
    }
    for _ in 0..opt.smooth_sweeps {
        report.phases.time("flip", || flip_pass(&mut w));
        report.phases.time("smooth", || smooth_pass(&mut w, field, &geom, opt));
    }
    if std::env::var_os("FEMMESH_VERBOSE").is_some() {
        let (mut long, mut short) = ([0; 2], [0; 2]);
        for e in all_edges_of(&w) {
            let l = elen(&w, e[0], e[1]);
            let surf = (w.vclass[e[0] as usize] != super::topo3::VClass3::Interior) as usize + (w.vclass[e[1] as usize] != super::topo3::VClass3::Interior) as usize;
            let k = (surf > 0) as usize;
            if l > opt.lmax {
                long[k] += 1;
            } else if l < opt.lmin {
                short[k] += 1;
            }
        }
        eprintln!("  left: long [interior, touching surface] {long:?}, short {short:?}");
    }
    report.unresolved_frozen = unresolved_frozen(&w, opt.lmax);
    w.check().map_err(|e| format!("internal error after adaptation: {e}"))?;
    let (out, metrics) = w.to_mesh();
    out.check().map_err(|e| format!("adapted mesh invalid: {e}"))?;
    Ok(Adapted3 { mesh: out, metrics, report })
}
