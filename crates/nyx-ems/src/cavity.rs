// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
//! Electromagnetic eigenmodes of 3D cavities and resonators (the `Eigenfrequency` study of EM waves
//! on a tetrahedral mesh): `curl(1/mu curl E) = (2 pi f / c)^2 eps E` with perfect electric
//! conductor on the exterior boundary and optional perfectly matched layers, kinetic-inductance
//! surface-impedance sheets (`surface_impedance`), interior `pec` sheets, and rank-1 lumped
//! Josephson junctions (`lumped_jj`).
//!
//! Elements are Nedelec elements of the first kind on tetrahedra, tangentially continuous:
//! order 1 has one degree of freedom per edge (the Whitney function `N_ab = l_a grad l_b - l_b grad l_a`),
//! order 2 has 20 per tetrahedron: per edge `N_ab` and the gradient `grad(l_a l_b)`, per face the two
//! functions `l_c N_ab` and `l_a N_bc`. Vertices of an edge or face are always taken in the order of
//! their global numbers, so neighbouring tetrahedra agree on shared edges and faces.
//!
//! With a complex stretch (PML) or surface resistance (`R_s > 0`) the matrices are complex symmetric
//! and the eigenvalues complex; the frequencies then carry a quality factor.

use crate::ibc::{
    assemble_ibc_and_jj_3d, reduce_real_csr, reduce_real_vec, ruiz_equilibrate, ruiz_unscale_and_rayleigh,
    AssembledJjPort, SheetCellKind, MU0,
};
use crate::modes::C0;
use crate::problem::{Problem, Study};
use crate::rom::{build_circuit_rom, build_tls_meta, compute_mode_energy_and_scale};
use femcore::fem::*;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::EigOpts;
use femcore::pml::cell_stretches;
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use rayon::prelude::*;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::time::Instant;

/// Records the time since `t` under the label `what` in the stage log (and on stderr when NYX_PROFILE is set), then restarts `t`.
fn lap(t: &mut Instant, what: &str) {
    femcore::progress::lap(what, t.elapsed().as_secs_f64());
    *t = Instant::now();
}

const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
const FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// Global numbering of the degrees of freedom.
pub struct Ned3 {
    pub order: usize,
    pub n_edges: usize,
    pub n_faces: usize,
    /// sorted vertex triples of all faces
    pub faces: Vec<[u32; 3]>,
    face_id: HashMap<[u32; 3], u32>,
    /// per cell, the global dofs of its basis functions (edge functions first, then face functions)
    pub cell_dofs: Vec<Vec<usize>>,
    /// ordering key of every vertex: edges and faces are oriented by it. Equal to the vertex number,
    /// except that the vertices of periodic boundaries share the key of their counterpart so that the
    /// basis functions of matching edges and faces correspond one to one
    pub vkey: Vec<u32>,
}

fn sorted3(a: u32, b: u32, c: u32) -> [u32; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

impl Ned3 {
    pub fn new(mesh: &Discretisation, order: usize, vkey: Vec<u32>) -> Ned3 {
        let mut face_id: HashMap<[u32; 3], u32> = HashMap::new();
        let mut faces: Vec<[u32; 3]> = vec![];
        let mut cell_faces: Vec<[u32; 4]> = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            let v = mesh.cell(ci);
            let mut cf = [0u32; 4];
            for (k, f) in FACES.iter().enumerate() {
                let key = sorted3(v[f[0]], v[f[1]], v[f[2]]);
                let id = *face_id.entry(key).or_insert_with(|| {
                    faces.push(key);
                    (faces.len() - 1) as u32
                });
                cf[k] = id;
            }
            cell_faces.push(cf);
        }
        let ne = mesh.n_edges();
        let nf = faces.len();
        let nb = if order == 1 { 6 } else { 20 };
        let mut cell_dofs = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            let mut d = Vec::with_capacity(nb);
            let eids = mesh.cell_edge_ids(ci);
            if order == 1 {
                d.extend(eids.iter().map(|&e| e as usize));
            } else {
                for &e in eids {
                    d.push(2 * e as usize);
                    d.push(2 * e as usize + 1);
                }
                for &f in &cell_faces[ci] {
                    d.push(2 * ne + 2 * f as usize);
                    d.push(2 * ne + 2 * f as usize + 1);
                }
            }
            cell_dofs.push(d);
        }
        Ned3 { order, n_edges: ne, n_faces: nf, faces, face_id, cell_dofs, vkey }
    }

    pub fn n_dofs(&self) -> usize {
        if self.order == 1 {
            self.n_edges
        } else {
            2 * self.n_edges + 2 * self.n_faces
        }
    }

    /// Ordering keys of the vertices of cell `ci`.
    pub fn keys(&self, mesh: &Discretisation, ci: usize) -> [u32; 4] {
        let v = mesh.cell(ci);
        [self.vkey[v[0] as usize], self.vkey[v[1] as usize], self.vkey[v[2] as usize], self.vkey[v[3] as usize]]
    }

    /// Dofs of the functions of face `v` (sorted triple of vertex numbers).
    pub fn face_dofs(&self, v: [u32; 3]) -> [usize; 2] {
        let f = self.face_id[&sorted3(v[0], v[1], v[2])] as usize;
        [2 * self.n_edges + 2 * f, 2 * self.n_edges + 2 * f + 1]
    }

    pub fn nb(&self) -> usize {
        if self.order == 1 {
            6
        } else {
            20
        }
    }

    /// Dofs of the tangential trace on a boundary face with vertices `v`: the face functions and the
    /// functions of its three edges.
    pub fn boundary_dofs(&self, mesh: &Discretisation, v: [u32; 3]) -> Vec<usize> {
        let mut out = vec![];
        for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[0], v[2])] {
            let e = mesh.edge_index(a, b);
            if self.order == 1 {
                out.push(e);
            } else {
                out.push(2 * e);
                out.push(2 * e + 1);
            }
        }
        if self.order == 2 {
            let f = self.face_id[&sorted3(v[0], v[1], v[2])] as usize;
            out.push(2 * self.n_edges + 2 * f);
            out.push(2 * self.n_edges + 2 * f + 1);
        }
        out
    }
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn axpy(s: f64, x: [f64; 3], y: [f64; 3]) -> [f64; 3] {
    [s * x[0] + y[0], s * x[1] + y[1], s * x[2] + y[2]]
}
fn scale3(s: f64, x: [f64; 3]) -> [f64; 3] {
    [s * x[0], s * x[1], s * x[2]]
}

/// Basis functions and their curls at the barycentric point `lam` of a tetrahedron with vertex
/// gradients `g` and global vertex numbers `gid`.
pub fn basis3(
    order: usize,
    lam: &[f64; 4],
    g: &[[f64; 3]; 4],
    gid: &[u32; 4],
    w: &mut Vec<[f64; 3]>,
    curl: &mut Vec<[f64; 3]>,
) {
    w.clear();
    curl.clear();
    // Whitney function of the (a, b) pair and its curl
    let whitney = |a: usize, b: usize| -> ([f64; 3], [f64; 3]) {
        (sub3(scale3(lam[a], g[b]), scale3(lam[b], g[a])), scale3(2.0, cross(g[a], g[b])))
    };
    for &(i, j) in EDGES.iter() {
        let (a, b) = if gid[i] < gid[j] { (i, j) } else { (j, i) };
        let (n, cn) = whitney(a, b);
        w.push(n);
        curl.push(cn);
        if order == 2 {
            w.push(axpy(lam[a], g[b], scale3(lam[b], g[a])));
            curl.push([0.0; 3]);
        }
    }
    if order == 2 {
        for f in FACES.iter() {
            let mut t = [f[0], f[1], f[2]];
            t.sort_by_key(|&k| gid[k]);
            let (a, b, cc) = (t[0], t[1], t[2]);
            let (nab, cab) = whitney(a, b);
            let (nbc, cbc) = whitney(b, cc);
            // l_c N_ab and l_a N_bc
            w.push(scale3(lam[cc], nab));
            curl.push(axpy(lam[cc], cab, cross(g[cc], nab)));
            w.push(scale3(lam[a], nbc));
            curl.push(axpy(lam[a], cbc, cross(g[a], nbc)));
        }
    }
}

type C3 = [[C; 3]; 3];

fn assemble_complex<F>(nc: usize, n: usize, f: F) -> CCsr
where
    F: Fn(std::ops::Range<usize>, &mut CCoo) + Sync,
{
    let nt = if nc < 500 { 1 } else { femcore::par::threads_up_to(16) };
    let chunk = nc.div_ceil(nt).max(1);
    let parts: Vec<CCsr> = (0..nt)
        .into_par_iter()
        .map(|t| {
            let mut coo = CCoo::new(n, n);
            f((t * chunk).min(nc)..((t + 1) * chunk).min(nc), &mut coo);
            coo.to_csr()
        })
        .collect();
    let mut it = parts.into_iter();
    let mut acc = it.next().unwrap();
    for p in it {
        acc = acc.add_scaled(&p, c(1.0, 0.0));
    }
    acc
}

/// Degree of the quadrature rule on cells with a curved edge (ADR 82): the integrand is rational in the barycentric coordinates there.
pub const CURVED_QUAD_DEGREE: usize = 6;

/// Adds the contribution of one quadrature point (weight `wq`, basis values `w` and curls `cu`) to the element matrix: the curl-curl
/// matrix for `which == 0`, the mass matrix otherwise. The same arithmetic serves straight and curved cells.
#[allow(clippy::too_many_arguments)]
fn add_point(
    which: usize,
    nb: usize,
    ke: &mut [C],
    w: &[[f64; 3]],
    cu: &[[f64; 3]],
    wq: f64,
    eps: &C3,
    mu_inv: &[C; 3],
) {
    for i in 0..nb {
        if which == 0 {
            let ci_ = cu[i];
            for j in 0..nb {
                let cj = cu[j];
                let mut acc = c(0.0, 0.0);
                for a in 0..3 {
                    acc += mu_inv[a] * (ci_[a] * cj[a]);
                }
                ke[i * nb + j] += acc * wq;
            }
        } else {
            let wi = w[i];
            for j in 0..nb {
                let wj = w[j];
                let mut acc = c(0.0, 0.0);
                for a in 0..3 {
                    for b in 0..3 {
                        acc += eps[a][b] * (wi[a] * wj[b]);
                    }
                }
                ke[i * nb + j] += acc * wq;
            }
        }
    }
}

/// The curl-curl matrix `K` (weights `1/mu` per axis) and the mass matrix `M` (permittivity tensor)
/// over all dofs.
pub fn assemble_maxwell(mesh: &Discretisation, ned: &Ned3, eps: &[C3], mu_inv: &[[C; 3]]) -> (CCsr, CCsr) {
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    // cells with a curved edge (ADR 82) are integrated with the geometry of their ten nodes; every other cell keeps the affine path
    let curved_cell: Vec<bool> = (0..mesh.n_cells()).map(|ci| mesh.cell_is_curved(ci)).collect();
    let cquad = quad_rule_collapsed(CURVED_QUAD_DEGREE);
    let nb = ned.nb();
    let n = ned.n_dofs();
    let nc = mesh.n_cells();
    let run = |which: usize| -> CCsr {
        assemble_complex(nc, n, |range, coo| {
            let (mut w, mut cu) = (Vec::new(), Vec::new());
            let mut ke = vec![c(0.0, 0.0); nb * nb];
            let mut gs = [[0.0; 3]; 4];
            for ci in range {
                let gid = ned.keys(mesh, ci);
                for b in 0..4 {
                    for k in 0..3 {
                        gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                    }
                }
                ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
                if curved_cell[ci] {
                    // quadratic geometry (ADR 82): gradients of the barycentric functions and the weight at every point
                    let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
                    for q in 0..cquad.w.len() {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        let (gq, det) = ct.jacobian(&lam);
                        basis3(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                        let wq = cquad.w[q] * det.abs() / 6.0;
                        add_point(which, nb, &mut ke, &w, &cu, wq, &eps[ci], &mu_inv[ci]);
                    }
                } else {
                    for q in 0..quad.w.len() {
                        let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                        basis3(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                        let wq = quad.w[q] * geo.vol[ci];
                        add_point(which, nb, &mut ke, &w, &cu, wq, &eps[ci], &mu_inv[ci]);
                    }
                }
                coo.push_block(&ned.cell_dofs[ci], &ned.cell_dofs[ci], &ke);
            }
        })
    };
    (run(0), run(1))
}

/// The discrete gradient over all nodes: one column per vertex (the Whitney edge coefficients
/// `phi_b - phi_a`, oriented by the vertex keys, of the gradient of its hat function) and, for order 2,
/// one per edge (the gradient of the edge bubble is the second function of the edge). Its columns lie
/// in the null space of the curl-curl matrix.
pub(crate) fn gradient_full(mesh: &Discretisation, ned: &Ned3) -> Csr {
    let nv = mesh.n_points();
    let ncols = if ned.order == 2 { nv + mesh.n_edges() } else { nv };
    let mut coo = femcore::linalg::csr::Coo::new(ned.n_dofs(), ncols);
    for (e, &(a, b)) in mesh.edges().iter().enumerate() {
        let (lo, hi) = if ned.vkey[a as usize] < ned.vkey[b as usize] { (a, b) } else { (b, a) };
        let dof = if ned.order == 1 { e } else { 2 * e };
        coo.push(dof, hi as usize, 1.0);
        coo.push(dof, lo as usize, -1.0);
        if ned.order == 2 {
            coo.push(2 * e + 1, nv + e, 1.0);
        }
    }
    coo.to_csr()
}

fn to_real(a: &CCsr) -> Csr {
    a.real_part()
}

fn is_real(a: &CCsr) -> bool {
    a.is_real()
}

fn to_complex(a: &Csr) -> CCsr {
    a.to_complex()
}

/// Equivalence classes of the vertices of periodic boundaries: the key of a vertex is the smallest
/// vertex number of its class.
fn vertex_keys(mesh: &Discretisation, matches: &[femcore::periodic::Match]) -> Vec<u32> {
    let n = mesh.n_points();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(p: &mut [u32], mut i: u32) -> u32 {
        while p[i as usize] != i {
            p[i as usize] = p[p[i as usize] as usize];
            i = p[i as usize];
        }
        i
    }
    for m in matches {
        for (&s, &t) in &m.map {
            let (rs, rt) = (find(&mut parent, s), find(&mut parent, t));
            if rs != rt {
                let (lo, hi) = (rs.min(rt), rs.max(rt));
                parent[hi as usize] = lo;
            }
        }
    }
    (0..n as u32).map(|i| find(&mut parent, i)).collect()
}

/// Relations between the dofs of the two sides of a periodic pair (all one to one because the basis
/// functions are oriented by the vertex keys).
fn ned_relations(
    mesh: &Discretisation,
    ned: &Ned3,
    m: &femcore::periodic::Match,
    phase: f64,
) -> Result<Vec<(usize, usize, C)>, String> {
    let f = femcore::periodic::phase_factor(phase);
    let mut out = vec![];
    let mut done_edges = std::collections::BTreeSet::new();
    for &fc in &m.src_facets {
        let v = mesh.facet(fc);
        let vs = [v[0], v[1], v[2]];
        let vt = [m.map[&vs[0]], m.map[&vs[1]], m.map[&vs[2]]];
        for (i, j) in [(0, 1), (1, 2), (0, 2)] {
            let key = (vs[i].min(vs[j]), vs[i].max(vs[j]));
            if !done_edges.insert(key) {
                continue;
            }
            let es = mesh.edge_index(vs[i], vs[j]);
            let et = mesh.edges().binary_search(&(vt[i].min(vt[j]), vt[i].max(vt[j]))).map_err(|_| {
                "an edge of the source boundary has no counterpart on the destination boundary".to_string()
            })?;
            if ned.order == 1 {
                out.push((et, es, f));
            } else {
                out.push((2 * et, 2 * es, f));
                out.push((2 * et + 1, 2 * es + 1, f));
            }
        }
        if ned.order == 2 {
            // the face functions are oriented by the vertex keys, which agree on both sides
            let ds = ned.face_dofs(vs);
            let dt = ned.face_dofs(vt);
            out.push((dt[0], ds[0], f));
            out.push((dt[1], ds[1], f));
        }
    }
    Ok(out)
}

/// Eigenfrequencies of the cavity for a `Cavity` study: the three stages in order. The result carries the complex field of every mode and
/// its frame, which is what the polarization analysis reads (`nyx-analysis`).
pub fn solve_cavity(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let asm = assemble_emw_eigen(mesh, problem)?;
    let solved = solve_emw_modes(&asm)?;
    let mut r = emw_eigen_fields(mesh, problem, &asm, &solved)?;
    r.seconds = t0.elapsed().as_secs_f64();
    Ok(r)
}

/// Stage 1 of the 3D EM eigenfrequency kernel: Nedelec numbering, `K` and `M` (PML stretch, anisotropy, IBC, JJ) reduced by the PEC
/// walls and the periodic pairs, and the gradient matrix that spans the null space of the curl-curl operator.
pub struct EmwAssembly {
    pub order: usize,
    pub n: usize,
    pub nred: usize,
    pub ned: Ned3,
    pub red: femcore::periodic::Reduction,
    pub kc: CCsr,
    pub mc: CCsr,
    pub g_red: CCsr,
    pub g_h: CCsr,
    pub f0: f64,
    pub sigma: f64,
    pub n_modes: usize,
    pub f_search: f64,
    pub has_cqed: bool,
    pub k0_red: Csr,
    pub k_ibc_red: Csr,
    pub c_ibc_red: Csr,
    pub m0_full: Csr,
    pub m0_red: Csr,
    pub jj_ports: Vec<AssembledJjPort>,
    pub cell_sheet: Vec<Option<SheetCellKind>>,
    pub periodic_node_relations: Vec<(usize, usize, C)>,
    pub mu_inv: Vec<[C; 3]>,
}

/// Stage 2 output: complex frequencies (Hz) with the index of their reduced vector, sorted by real part.
pub struct EmwSolve {
    pub freqs: Vec<(C, usize)>,
    pub vectors: Vec<Vec<C>>,
    pub dim: usize,
    pub converged: bool,
}

pub fn assemble_emw_eigen(mesh: &Discretisation, problem: &Problem) -> Result<EmwAssembly, String> {
    if mesh.dim != 3 {
        return Err("the cavity study of EM waves needs a 3D mesh".into());
    }
    let Study::Cavity { n_modes, f_search } = problem.study.clone() else {
        return Err("EM waves on a 3D mesh use the cavity study".into());
    };
    let order = if problem.order == 1 { 1 } else { 2 };
    let tens = problem.materials.cell_complex_permittivity(mesh, Some(1.0))?;
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let has_pml = problem.materials.has_pml();
    let stretch = if has_pml { Some(cell_stretches(mesh, &problem.materials.pml)) } else { None };
    let mut eps: Vec<C3> = Vec::with_capacity(mesh.n_cells());
    let mut mu_inv: Vec<[C; 3]> = Vec::with_capacity(mesh.n_cells());
    for ci in 0..mesh.n_cells() {
        let (sx, sy, sz) =
            stretch.as_ref().map_or((c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1], s[ci][2]));
        // PML: eps' = Lambda eps, mu' = Lambda mu with Lambda = diag(sy sz / sx, sx sz / sy, sx sy / sz)
        let lam = [sy * sz / sx, sx * sz / sy, sx * sy / sz];
        let pml_here = stretch.is_some()
            && (sx.re != 1.0 || sx.im != 0.0 || sy.re != 1.0 || sy.im != 0.0 || sz.re != 1.0 || sz.im != 0.0);
        let (t_re, t_im) = tens[ci];
        // a stretch S = diag(sx, sy, sz) turns eps into det(S) eps_ab / (s_a s_b): the diagonal gets `lam`, off-diagonal terms of an
        // anisotropic crystal are scaled consistently
        let (sv, det) = ([sx, sy, sz], sx * sy * sz);
        let mut e = [[c(0.0, 0.0); 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let ec = c(t_re[a][b], -t_im[a][b]);
                e[a][b] = if pml_here { det * ec / (sv[a] * sv[b]) } else { ec };
            }
        }
        eps.push(e);
        let one = c(1.0, 0.0);
        let mu0 = c(1.0 / mur[ci], 0.0);
        mu_inv.push([
            mu0 / if pml_here { lam[0] } else { one },
            mu0 / if pml_here { lam[1] } else { one },
            mu0 / if pml_here { lam[2] } else { one },
        ]);
    }
    let mut tl = Instant::now();
    // periodic pairs
    let mut matches = vec![];
    for p in &problem.periodic {
        matches.push(femcore::periodic::match_vertices(mesh, p)?);
    }
    let vkey = vertex_keys(mesh, &matches);
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let (k, m) = assemble_maxwell(mesh, &ned, &eps, &mu_inv);
    lap(&mut tl, "emw3d: numbering and assembly");

    // shift: the requested frequency; without one, the frequency of a wave that fits the cavity
    let (lo, hi) = mesh.bbox();
    let size = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(1e-300);
    let f0 = if f_search > 0.0 { f_search } else { C0 / (2.0 * size) };
    let sigma = (2.0 * PI * f0 / C0).powi(2);

    let has_cqed = !problem.surface_impedance.is_empty() || !problem.lumped_jj.is_empty();
    let has_sheets = has_cqed || !problem.pec.is_empty();
    let mut ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d(mesh, &ned, &problem.pec, &problem.surface_impedance, &problem.lumped_jj, f0)?)
    } else {
        None
    };

    // perfect electric conductor on the exterior boundary, except the periodic boundaries, magnetic walls (`Problem::pmc`, natural),
    // `surface_impedance` sheets and `lumped_jj` sheets
    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let ibc_labels: Vec<i32> = problem.surface_impedance.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();

    let mut fixed_set = std::collections::BTreeSet::new();
    let mut ext_pec_vertex = std::collections::BTreeSet::new();
    let mut pec_edge = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        if exterior[f]
            && !periodic_labels.contains(&label)
            && !problem.pmc.contains(&label)
            && !ibc_labels.contains(&label)
            && !jj_labels.contains(&label)
        {
            let fv = mesh.facet(f);
            for d in ned.boundary_dofs(mesh, [fv[0], fv[1], fv[2]]) {
                fixed_set.insert(d);
            }
            for &v in fv {
                ext_pec_vertex.insert(v as usize);
            }
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                pec_edge.insert(mesh.edge_index(a, b));
            }
        }
    }
    let mut cond_dsu_edges: Vec<(usize, usize)> = Vec::new();
    let mut cond_vertices = std::collections::BTreeSet::new();
    if let Some(ref ij) = ibc_jj {
        // Explicit `pec` facets: strong E_|| = 0 elimination in `fixed_set`
        for &fv in &ij.pec_tris {
            for d in ned.boundary_dofs(mesh, fv) {
                fixed_set.insert(d);
            }
            for v in fv {
                cond_vertices.insert(v as usize);
            }
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                pec_edge.insert(mesh.edge_index(a, b));
                cond_dsu_edges.push((a as usize, b as usize));
            }
        }
        // `surface_impedance` facets: keep Nedelec DOFs FREE in `red`
        for &fv in &ij.ibc_tris {
            for v in fv {
                cond_vertices.insert(v as usize);
            }
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                pec_edge.insert(mesh.edge_index(a, b));
                cond_dsu_edges.push((a as usize, b as usize));
            }
        }
        // `lumped_jj` facets: keep Nedelec DOFs FREE in `red`.
        // Include all vertices of `jj_tris` in `cond_vertices` / `cond_dsu_edges` (for floating-component common-mode
        // identification), while for `order == 2` edge bubbles on `rank1` JJ facets only fix edge bubbles that touch
        // `pec`/`ibc` or have `(G^T b_j)_{nv+e} != 0` so interior JJ edge-bubble gradients are projected out by `g_red`.
        let g_full_r = gradient_full(mesh, &ned);
        let gt_r = g_full_r.transpose();
        let nv = mesh.n_points();
        for port in &ij.jj_ports {
            for &fv in &port.triangles {
                for v in fv {
                    cond_vertices.insert(v as usize);
                }
                for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                    cond_dsu_edges.push((a as usize, b as usize));
                    if port.bc.formulation == "distributed_sheet" {
                        pec_edge.insert(mesh.edge_index(a, b));
                    }
                }
            }
            if port.bc.formulation != "distributed_sheet" && order == 2 {
                let gt_b = gt_r.matvec(&port.b_full);
                let b_max = gt_b.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
                if b_max > 0.0 {
                    for e in 0..mesh.n_edges() {
                        if gt_b[nv + e].abs() > 1e-8 * b_max {
                            pec_edge.insert(e);
                        }
                    }
                }
            }
        }
    }

    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    // reduction of the dofs (PEC removed, periodic dofs merged)
    let mut relations = vec![];
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        relations.extend(ned_relations(mesh, &ned, mt, p.phase)?);
    }
    let red = femcore::periodic::Reduction::build(n, &relations, &fixed);
    let pmat = red.prolong();
    let nred = red.n_reduced();
    if n_modes < 1 || n_modes + 1 >= nred {
        return Err("number of modes must be between 1 and the number of free dofs".into());
    }

    let k0_c_red = femcore::periodic::congruence(&pmat, &k);
    let m0_c_red = femcore::periodic::congruence(&pmat, &m);
    let k0_red = to_real(&k0_c_red);
    let m0_red = to_real(&m0_c_red);
    let m0_full = if has_cqed { to_real(&m) } else { Csr::zeros(0, 0) };

    let (kc, mc, k_ibc_red, c_ibc_red, jj_ports, cell_sheet) = if let Some(ref mut ij) = ibc_jj {
        let mut k_tot = k.clone();
        let mut m_tot = m.clone();
        if ij.kc_ibc_full.nnz() > 0 {
            k_tot = k_tot.add_scaled(&ij.kc_ibc_full, c(1.0, 0.0));
        }
        if ij.k_jj_full.nnz() > 0 {
            k_tot = k_tot.add_scaled(&to_complex(&ij.k_jj_full), c(1.0, 0.0));
        }
        if ij.m_jj_full.nnz() > 0 {
            m_tot = m_tot.add_scaled(&to_complex(&ij.m_jj_full), c(1.0, 0.0));
        }
        let kc = if has_cqed { femcore::periodic::congruence(&pmat, &k_tot) } else { k0_c_red };
        let mc = if has_cqed { femcore::periodic::congruence(&pmat, &m_tot) } else { m0_c_red };
        let k_ibc_red = reduce_real_csr(&red, &pmat, &ij.k_ibc_full);
        let c_ibc_red = reduce_real_csr(&red, &pmat, &ij.c_ibc_full);
        for port in &mut ij.jj_ports {
            port.b_red = reduce_real_vec(&red, &port.b_full);
            port.k_jj_red = reduce_real_csr(&red, &pmat, &port.k_jj_red);
            port.m_jj_red = reduce_real_csr(&red, &pmat, &port.m_jj_red);
        }
        (kc, mc, k_ibc_red, c_ibc_red, std::mem::take(&mut ij.jj_ports), std::mem::take(&mut ij.cell_sheet))
    } else {
        (k0_c_red, m0_c_red, Csr::zeros(nred, nred), Csr::zeros(nred, nred), Vec::new(), vec![None; mesh.n_cells()])
    };

    // null space of the curl-curl operator: gradients of the interior (non-PEC, non-IBC, non-JJ) nodal functions with the same periodic structure,
    // plus one common-mode constant-potential vector `1_{G_m}` per completely floating galvanic component `G_m` (not touching exterior PEC).
    let nv = mesh.n_points();
    let n_nodes = if order == 2 { nv + mesh.n_edges() } else { nv };
    let mut periodic_node_relations = vec![];
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        periodic_node_relations.extend(femcore::periodic::nodal_relations(mesh, order, 1, mt, p.phase)?);
    }
    let mut node_rel = periodic_node_relations.clone();
    let mut pec_vertex = ext_pec_vertex.clone();
    if !cond_vertices.is_empty() {
        let mut parent: Vec<usize> = (0..nv).collect();
        let find = |p: &mut [usize], mut i: usize| -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        };
        for &(a, b) in &cond_dsu_edges {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                let (lo, hi) = (ra.min(rb), ra.max(rb));
                parent[hi] = lo;
            }
        }
        let mut comp_touches_ext = std::collections::BTreeSet::new();
        for &v in &cond_vertices {
            if ext_pec_vertex.contains(&v) {
                comp_touches_ext.insert(find(&mut parent, v));
            }
        }
        for &v in &cond_vertices {
            let root = find(&mut parent, v);
            if comp_touches_ext.contains(&root) {
                pec_vertex.insert(v);
            } else if v != root {
                node_rel.push((v, root, c(1.0, 0.0)));
            }
        }
    }
    let mut node_fixed: Vec<usize> = pec_vertex.iter().copied().collect();
    if order == 2 {
        node_fixed.extend(pec_edge.iter().map(|&e| nv + e));
    }
    // nodal_relations numbers nodes as vertices then edge midpoints, as `gradient_full` does
    let nred_nodes = femcore::periodic::Reduction::build(n_nodes, &node_rel, &node_fixed);
    let g_full = to_complex(&gradient_full(mesh, &ned));
    let g_nodes = femcore::periodic::ccsr_mul(&g_full, &nred_nodes.prolong());
    // rows of the reduced dofs: x_full = P x_red has x_red equal to the root dofs
    let g_red = {
        let (mut indptr, mut indices, mut data) = (vec![0usize], vec![], vec![]);
        for &r in &red.roots {
            for q in g_nodes.indptr()[r]..g_nodes.indptr()[r + 1] {
                indices.push(g_nodes.indices()[q] as usize);
                data.push(g_nodes.values()[q]);
            }
            indptr.push(indices.len());
        }
        CCsr::from_raw(nred, nred_nodes.n_reduced(), indptr, indices, data)
    };
    let g_h = femcore::periodic::ccsr_adjoint(&g_red);
    lap(&mut tl, "emw3d: reduction and gradient matrix");

    Ok(EmwAssembly {
        order,
        n,
        nred,
        ned,
        red,
        kc,
        mc,
        g_red,
        g_h,
        f0,
        sigma,
        n_modes,
        f_search,
        has_cqed,
        k0_red,
        k_ibc_red,
        c_ibc_red,
        m0_full,
        m0_red,
        jj_ports,
        cell_sheet,
        periodic_node_relations,
        mu_inv,
    })
}

/// Stage 2: the eigenpairs `K x = k0^2 M x` nearest the shift with the gradients projected out (real symmetric route with
/// shift-invert, or the complex route for absorbing layers, Floquet phases, and lossy `R_s > 0` IBCs).
pub fn solve_emw_modes(asm: &EmwAssembly) -> Result<EmwSolve, String> {
    let (kc, mc, g_red, g_h, sigma, n_modes) = (&asm.kc, &asm.mc, &asm.g_red, &asm.g_h, asm.sigma, asm.n_modes);
    let mut tl = Instant::now();
    // eigenvalue lambda = k0^2: K x = lambda M x
    let mut freqs: Vec<(C, usize)> = vec![];
    let (vectors, dim, converged): (Vec<Vec<C>>, usize, bool);

    if asm.has_cqed {
        let (d_ruiz, kc_eq, mc_eq, g_eq, gh_eq) = ruiz_equilibrate(kc, mc, g_red, sigma);
        if is_real(&kc_eq) && is_real(&mc_eq) && is_real(&g_eq) {
            let (kr, mr, gr) = (to_real(&kc_eq), to_real(&mc_eq), to_real(&g_eq));
            lap(&mut tl, "emw3d: matrices (Ruiz equilibrated)");
            let opts = EigOpts { tol: 1e-12, ..EigOpts::default() };
            let pairs = if gr.ncols() > 0 {
                let gt = gr.transpose();
                let s_mat = femcore::linalg::amg_elastic::spmm(&gt, &femcore::linalg::amg_elastic::spmm(&mr, &gr));
                let s_solver = femcore::linalg::driver::DirectSolver::new(&s_mat, true)?;
                let project = |x: &mut Vec<f64>| {
                    let t = mr.matvec(x);
                    let g = gt.matvec(&t);
                    let y = s_solver.solve(&g);
                    let gy = gr.matvec(&y);
                    for (a, b) in x.iter_mut().zip(&gy) {
                        *a -= b;
                    }
                };
                lap(&mut tl, "emw3d: projection factorisation");
                femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_modes, &opts, Some(&project))?
            } else {
                femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_modes, &opts, None)?
            };
            let mut unscaled = Vec::with_capacity(pairs.vectors.len());
            for (i, v) in pairs.vectors.iter().enumerate() {
                let yc: Vec<C> = v.iter().map(|&x| c(x, 0.0)).collect();
                let (x, f) = ruiz_unscale_and_rayleigh(&d_ruiz, &yc, kc, mc);
                freqs.push((c(f.re, 0.0), i));
                unscaled.push(x);
            }
            vectors = unscaled;
            dim = pairs.dim;
            converged = pairs.converged;
        } else {
            let pairs = if g_eq.ncols() > 0 {
                let s_c = femcore::periodic::ccsr_mul(&gh_eq, &femcore::periodic::ccsr_mul(&mc_eq, &g_eq));
                let s_solver = femcore::linalg::driver::DirectSolver::new(&s_c, false)?;
                let project = |x: &mut Vec<C>| {
                    let t = mc_eq.matvec(x);
                    let g = gh_eq.matvec(&t);
                    let y = s_solver.solve(&g);
                    let gy = g_eq.matvec(&y);
                    for (a, b) in x.iter_mut().zip(&gy) {
                        *a -= *b;
                    }
                };
                lap(&mut tl, "emw3d: projection factorisation");
                femcore::linalg::driver::eigs_complex_projected(
                    &kc_eq,
                    &mc_eq,
                    c(sigma, 0.0),
                    n_modes,
                    1e-11,
                    500,
                    Some(&project),
                )?
            } else {
                femcore::linalg::driver::eigs_complex_projected(
                    &kc_eq,
                    &mc_eq,
                    c(sigma, 0.0),
                    n_modes,
                    1e-11,
                    500,
                    None,
                )?
            };
            let mut unscaled = Vec::with_capacity(pairs.vectors.len());
            for (i, v) in pairs.vectors.iter().enumerate() {
                let (x, f) = ruiz_unscale_and_rayleigh(&d_ruiz, v, kc, mc);
                unscaled.push(x);
                if pairs.residuals.get(i).is_some_and(|r| *r > 1e-5) {
                    continue;
                }
                freqs.push((f, i));
            }
            vectors = unscaled;
            dim = pairs.dim;
            converged = pairs.converged;
        }
    } else if is_real(kc) && is_real(mc) && is_real(g_red) {
        let (kr, mr, gr) = (to_real(kc), to_real(mc), to_real(g_red));
        lap(&mut tl, "emw3d: matrices");
        let pairs = if gr.ncols() > 0 {
            let gt = gr.transpose();
            let s_mat = femcore::linalg::amg_elastic::spmm(&gt, &femcore::linalg::amg_elastic::spmm(&mr, &gr));
            let s_solver = femcore::linalg::driver::DirectSolver::new(&s_mat, true)?;
            let project = |x: &mut Vec<f64>| {
                let t = mr.matvec(x);
                let g = gt.matvec(&t);
                let y = s_solver.solve(&g);
                let gy = gr.matvec(&y);
                for (a, b) in x.iter_mut().zip(&gy) {
                    *a -= b;
                }
            };
            lap(&mut tl, "emw3d: projection factorisation");
            femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_modes, &EigOpts::default(), Some(&project))?
        } else {
            femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_modes, &EigOpts::default(), None)?
        };
        for (i, l) in pairs.values.iter().enumerate() {
            freqs.push((c(l.max(0.0).sqrt() * (C0 / (2.0 * PI)), 0.0), i));
        }
        vectors = pairs.vectors.iter().map(|v| v.iter().map(|&x| c(x, 0.0)).collect()).collect();
        dim = pairs.dim;
        converged = pairs.converged;
    } else {
        // complex: absorbing layers or Floquet phases. Projection x - G S^-1 G^H M x with S = G^H M G
        let pairs = if g_red.ncols() > 0 {
            let s_c = femcore::periodic::ccsr_mul(g_h, &femcore::periodic::ccsr_mul(mc, g_red));
            let s_solver = femcore::linalg::driver::DirectSolver::new(&s_c, false)?;
            let project = |x: &mut Vec<C>| {
                let t = mc.matvec(x);
                let g = g_h.matvec(&t);
                let y = s_solver.solve(&g);
                let gy = g_red.matvec(&y);
                for (a, b) in x.iter_mut().zip(&gy) {
                    *a -= *b;
                }
            };
            lap(&mut tl, "emw3d: projection factorisation");
            femcore::linalg::driver::eigs_complex_projected(kc, mc, c(sigma, 0.0), n_modes, 1e-9, 500, Some(&project))?
        } else {
            femcore::linalg::driver::eigs_complex_projected(kc, mc, c(sigma, 0.0), n_modes, 1e-9, 500, None)?
        };
        for (i, l) in pairs.values.iter().enumerate() {
            if pairs.residuals.get(i).is_some_and(|r| *r > 1e-6) {
                continue;
            }
            let mut f = l.sqrt() * (C0 / (2.0 * PI));
            if f.re < 0.0 {
                f = -f;
            }
            freqs.push((f, i));
        }
        vectors = pairs.vectors.clone();
        dim = pairs.dim;
        converged = pairs.converged;
    }
    lap(&mut tl, "emw3d: eigenvalue iteration");
    freqs.sort_by(|a, b| a.0.re.partial_cmp(&b.0.re).unwrap());
    Ok(EmwSolve { freqs, vectors, dim, converged })
}

/// Stage 3: plot fields and per-cell complex E and H (plus sheet K when sheets are present) of every mode at the quadratic nodes,
/// the frequencies and Q, the frame of the fields for the polarization analysis, and `CircuitRom` + `TlsMeta` when superconducting
/// / cQED conditions are present.
pub fn emw_eigen_fields(
    mesh: &Discretisation,
    problem: &Problem,
    asm: &EmwAssembly,
    solved: &EmwSolve,
) -> Result<SolveResult, String> {
    let (order, n, nred, f0, f_search) = (asm.order, asm.n, asm.nred, asm.f0, asm.f_search);
    let (ned, red) = (&asm.ned, &asm.red);
    let (freqs, vectors, dim, converged) = (&solved.freqs, &solved.vectors, solved.dim, solved.converged);
    // fields at the quadratic nodes
    let geo = geometry(mesh);
    let lamn = local_bary(3, 2);
    let nloc = lamn.len();
    let nc = mesh.n_cells();
    let nb = ned.nb();
    let mut solutions = vec![];
    let mut scalars = vec![];
    let mut mode_energy_rows = Vec::with_capacity(freqs.len());
    let mut max_equip_err = 0.0_f64;
    let has_sheets = asm.has_cqed || !problem.pec.is_empty();
    let (axes, left_handed) = femcore::section::axes_of_3d_periodic(&problem.periodic, mesh);
    for (j, (f, i)) in freqs.iter().enumerate() {
        let x = red.expand(&vectors[*i]);
        let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let mut curl_e = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];
        for ci in 0..nc {
            let gid = ned.keys(mesh, ci);
            for b in 0..4 {
                for kk in 0..3 {
                    gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
                }
            }
            let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
            for (a, l) in lamn.iter().enumerate() {
                let lam = [l[0], l[1], l[2], l[3]];
                match &curved {
                    // gradients at the point for a cell with a curved edge (ADR 82); the plotted nodes stay at the straight positions
                    Some(ct) => basis3(order, &lam, &ct.jacobian(&lam).0, &gid, &mut w, &mut cu),
                    None => basis3(order, &lam, &gs, &gid, &mut w, &mut cu),
                }
                let mut val = [c(0.0, 0.0); 3];
                let mut cval = [c(0.0, 0.0); 3];
                for s in 0..nb {
                    let coef = x[ned.cell_dofs[ci][s]];
                    for kk in 0..3 {
                        val[kk] += coef * w[s][kk];
                        cval[kk] += coef * cu[s][kk];
                    }
                }
                e[ci * nloc + a] = val;
                curl_e[ci * nloc + a] = cval;
            }
        }
        let big =
            e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
        let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
        let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, 2);
        let fields = vec![
            Field::scalar("|E|", "a.u.", avg(mag.clone())),
            Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
            Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
            Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
        ];

        // Physical U_elec = 1 J normalization and exact modal energy partition
        let (s_e, energy_row, equip_err) = compute_mode_energy_and_scale(*f, &vectors[*i], rot, asm);
        max_equip_err = max_equip_err.max(equip_err);

        let mut cells = CellFields::new(mesh, 2);
        let (e_scale, e_unit) = if asm.has_cqed { (rot * s_e, "V/m") } else { (rot * s, "a.u.") };
        for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
            let z: Vec<C> = e.iter().map(|v| v[k] * e_scale).collect();
            cells.push_complex(*name, e_unit, z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
        }

        // Magnetic field H_m = (i / (omega_m mu_0 mu_r)) curl E_m in "A/m" at U_elec = 1 J
        let omega_m = (2.0 * PI * f.re).max(1e-300);
        let h_prefactor = c(0.0, 1.0 / (omega_m * MU0)) * (rot * s_e);
        let mut h_field = vec![[c(0.0, 0.0); 3]; nc * nloc];
        for ci in 0..nc {
            let mui = asm.mu_inv[ci];
            for a in 0..nloc {
                let idx = ci * nloc + a;
                let ce = curl_e[idx];
                h_field[idx] =
                    [ce[0] * mui[0] * h_prefactor, ce[1] * mui[1] * h_prefactor, ce[2] * mui[2] * h_prefactor];
            }
        }
        for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
            let z: Vec<C> = h_field.iter().map(|v| v[k]).collect();
            cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
        }

        // Surface current density Kx, Ky, Kz in "A/m" on cells adjacent to conductive sheets
        if has_sheets {
            let mut k_field = vec![[c(0.0, 0.0); 3]; nc * nloc];
            for ci in 0..nc {
                let Some(sheet_kind) = asm.cell_sheet[ci] else { continue };
                match sheet_kind {
                    SheetCellKind::Ibc { inv_l_eff, g_eff, n_hat } => {
                        let y_s = c(g_eff, -inv_l_eff / omega_m);
                        for a in 0..nloc {
                            let idx = ci * nloc + a;
                            let ep = [e[idx][0] * (rot * s_e), e[idx][1] * (rot * s_e), e[idx][2] * (rot * s_e)];
                            let en = ep[0] * n_hat[0] + ep[1] * n_hat[1] + ep[2] * n_hat[2];
                            let e_par = [ep[0] - en * n_hat[0], ep[1] - en * n_hat[1], ep[2] - en * n_hat[2]];
                            k_field[idx] = [e_par[0] * y_s, e_par[1] * y_s, e_par[2] * y_s];
                        }
                    }
                    SheetCellKind::Jj { port_index, d_hat, width_m } => {
                        let v_p = energy_row.v_mj[port_index];
                        let i_p = energy_row.i_mj[port_index];
                        let cj = asm.jj_ports[port_index].bc.c_j;
                        let i_tot = c(i_p[0], i_p[1]) + c(0.0, omega_m * cj) * c(v_p[0], v_p[1]);
                        let k_amp = i_tot / width_m.max(1e-300);
                        for a in 0..nloc {
                            k_field[ci * nloc + a] = [k_amp * d_hat[0], k_amp * d_hat[1], k_amp * d_hat[2]];
                        }
                    }
                    SheetCellKind::Pec { n_hat } => {
                        for a in 0..nloc {
                            let idx = ci * nloc + a;
                            let hm = h_field[idx];
                            k_field[idx] = [
                                n_hat[1] * hm[2] - n_hat[2] * hm[1],
                                n_hat[2] * hm[0] - n_hat[0] * hm[2],
                                n_hat[0] * hm[1] - n_hat[1] * hm[0],
                            ];
                        }
                    }
                }
            }
            for (k, name) in ["Kx", "Ky", "Kz"].iter().enumerate() {
                let z: Vec<C> = k_field.iter().map(|v| v[k]).collect();
                cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
            }
        }

        mode_energy_rows.push(energy_row);

        if f.im.abs() > 1e-12 * f.re.abs().max(1.0) {
            let q = f.re / (2.0 * f.im.abs());
            solutions.push(Solution {
                label: format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, q),
                fields,
                cell_fields: cell_fields_of(&cells),
                value: Some(f.re),
                eigenvalue: Some([f.re, f.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(q)));
        } else {
            solutions.push(Solution {
                label: format!("mode {}: {:.6} GHz", j + 1, f.re / 1e9),
                fields,
                cell_fields: cell_fields_of(&cells),
                value: Some(f.re),
                eigenvalue: Some([f.re, f.im]),
            });
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(f.re)));
        }
    }
    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let mut log = vec![format!(
        "Nedelec order {order}, {n} dofs ({nred} after constraints), shift {:.4e} Hz, Krylov dimension {dim}{}",
        f0,
        if converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    if !problem.periodic.is_empty() {
        log.push(format!("{} Floquet-periodic pair(s)", problem.periodic.len()));
    }
    if f_search <= 0.0 {
        log.push("no search frequency given: modes nearest the fundamental of a box of this size".into());
    }
    // the frame the polarization is read in: the (lateral, normal, propagation) axes of the problem
    scalars.push(
        nyx_physics::FieldFrame {
            lateral: axes[0],
            normal: axes[1],
            propagation: axes[2],
            left_handed,
            outside: vec![],
        }
        .scalar(),
    );
    if asm.has_cqed {
        let rom = build_circuit_rom(mesh, problem, asm, vectors, freqs, mode_energy_rows)?;
        scalars.push(rom.scalar());
        scalars.push(("equipartition_error".into(), Scalar::Num(max_equip_err)));
    }
    if has_sheets {
        let tls = build_tls_meta(mesh, problem);
        scalars.push(tls.scalar());
    }
    log.push(format!(
        "polarization axes: lateral {}, normal {}, propagation {}",
        ["x", "y", "z"][axes[0]],
        ["x", "y", "z"][axes[1]],
        ["x", "y", "z"][axes[2]]
    ));
    Ok(SolveResult {
        physics: "ems".into(),
        study: "cavity".into(),
        order: 2,
        dim: 3,
        solutions,
        scalars,
        log,
        dofs: n,
        seconds: 0.0,
    })
}
