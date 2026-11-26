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
    assemble_ibc_and_jj_3d_cavity, reduce_real_csr, reduce_real_vec, ruiz_equilibrate, ruiz_unscale_and_rayleigh,
    ruiz_unscale_and_rayleigh_damped, AssembledJjPort, PerturbativeIbcEvaluator, PerturbativeIbcTri, SheetCellKind,
    EPS0, MU0,
};
use crate::mixed_ops::{
    add_maxwell_point_c3, assemble_ibc_and_jj_3d_mixed_cavity, assemble_lumped_ports_3d_mixed, assemble_maxwell_mixed,
    axes_of_3d_periodic_mixed, eigs_mixed_curl_curl_ams, gradient_mixed_full, iso_mu_inv_c3, match_vertices_mixed,
    nodal_relations_mixed_3d, scalar_nullspace_count_mixed, simplex_cell_materials_3d, symmetry_facet_flags_mixed,
    validate_port_contacts_3d_mixed, vertex_keys_mixed, EmPmlTable, MixedCellMaterials,
};
use crate::modes::C0;
use crate::ports::{assemble_lumped_ports_3d, validate_port_contacts_3d, AssembledLumpedPort};
use crate::problem::{symmetry_facet_flags, Problem, Study};
use crate::rom::{build_circuit_rom, build_circuit_rom_mixed, build_tls_meta, compute_mode_energy_and_scale};
use femcore::edge_mixed::{eval_cell_basis_3d, ned_mixed_periodic_relations, nedelec_nloc, NedEvalScratch, NedMixed};
use femcore::fem::*;
use femcore::fem_mixed::{quad_rule_mixed, IsoparametricCell};
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::RCsr as Csr;
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::EigOpts;
use femcore::mesh::mixed::{FacetType, MixedDiscretisation, QUAD_EDGES, TRI_EDGES};
use nyx_physics::{
    cell_fields_of, Discretisation, Field, Mesh, MixedToSimplexProjector, Scalar, Solution, SolveResult,
};
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
    pub(crate) fn empty(order: usize) -> Ned3 {
        Ned3 {
            order,
            n_edges: 0,
            n_faces: 0,
            faces: Vec::new(),
            face_id: HashMap::new(),
            cell_dofs: Vec::new(),
            vkey: Vec::new(),
        }
    }

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

    /// Composite ordering keys `(vkey[v], v)` of cell `ci`, breaking periodic ties deterministically by vertex ID.
    pub fn key_pairs(&self, mesh: &Discretisation, ci: usize) -> [(u32, u32); 4] {
        let v = mesh.cell(ci);
        [
            (self.vkey[v[0] as usize], v[0]),
            (self.vkey[v[1] as usize], v[1]),
            (self.vkey[v[2] as usize], v[2]),
            (self.vkey[v[3] as usize], v[3]),
        ]
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
    let keyed = [(gid[0], 0), (gid[1], 1), (gid[2], 2), (gid[3], 3)];
    basis3_keyed(order, lam, g, &keyed, w, curl);
}

/// Like [`basis3`], taking composite `(vkey, vertex_id)` keys so edges/faces touching periodic
/// counterparts on 1-layer extrusions break `vkey` ties identically across sharing tetrahedra.
pub fn basis3_keyed(
    order: usize,
    lam: &[f64; 4],
    g: &[[f64; 3]; 4],
    gid: &[(u32, u32); 4],
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

pub(crate) type C3 = [[C; 3]; 3];

/// Fixed cell chunk size for thread-count-invariant matrix assembly (ADR 11, FR-29).
const FIXED_CELL_CHUNK: usize = 2048;

/// Merges CSR blocks via a deterministic pairwise binary reduction tree so floating-point summation
/// order is bitwise identical across any `NYX_THREADS >= 1`.
fn reduce_tree_ccsr(mut parts: Vec<CCsr>) -> CCsr {
    let use_par = femcore::par::threads() > 1;
    while parts.len() > 1 {
        let n_pairs = parts.len() / 2;
        let tail = if parts.len() % 2 == 1 { parts.pop() } else { None };
        let mut pairs = Vec::with_capacity(n_pairs);
        let mut it = parts.into_iter();
        for _ in 0..n_pairs {
            let a = it.next().unwrap();
            let b = it.next().unwrap();
            pairs.push((a, b));
        }
        let mut next: Vec<CCsr> = if use_par && pairs.len() > 1 {
            pairs.into_par_iter().map(|(a, b)| a.add_scaled(&b, c(1.0, 0.0))).collect()
        } else {
            pairs.into_iter().map(|(a, b)| a.add_scaled(&b, c(1.0, 0.0))).collect()
        };
        if let Some(t) = tail {
            next.push(t);
        }
        parts = next;
    }
    parts.pop().unwrap()
}

pub(crate) fn assemble_complex<F>(nc: usize, n: usize, f: F) -> CCsr
where
    F: Fn(std::ops::Range<usize>, &mut CCoo) + Sync,
{
    if nc == 0 {
        return CCsr::zeros(n, n);
    }
    let n_chunks = nc.div_ceil(FIXED_CELL_CHUNK);
    let ranges: Vec<std::ops::Range<usize>> =
        (0..n_chunks).map(|k| (k * FIXED_CELL_CHUNK)..((k + 1) * FIXED_CELL_CHUNK).min(nc)).collect();
    let parts: Vec<CCsr> = if n_chunks == 1 || femcore::par::threads() <= 1 {
        ranges
            .into_iter()
            .map(|range| {
                let mut coo = CCoo::new(n, n);
                f(range, &mut coo);
                coo.to_csr()
            })
            .collect()
    } else {
        ranges
            .into_par_iter()
            .map(|range| {
                let mut coo = CCoo::new(n, n);
                f(range, &mut coo);
                coo.to_csr()
            })
            .collect()
    };
    reduce_tree_ccsr(parts)
}

/// Degree of the quadrature rule on cells with a curved edge (ADR 82): the integrand is rational in the barycentric coordinates there.
pub const CURVED_QUAD_DEGREE: usize = 6;

/// Safe physical barycentric gradients on a curved tetrahedron: on Barsoum quarter-point (`xi_QP = 0.25`)
/// singular rim vertices (`det(J) == 0` right on the vertex), pulls `lam` slightly into the cell interior
/// so nodal/surface evaluations remain finite without disturbing smooth curved elements (FR-32).
pub(crate) fn safe_curved_jacobian(ct: &CurvedTet, lam: &[f64; 4], gs: &[[f64; 3]; 4]) -> [[f64; 3]; 4] {
    let (g, det) = ct.jacobian(lam);
    if det.is_finite() && det.abs() > 1e-30 && g.iter().flatten().all(|v| v.is_finite()) {
        return g;
    }
    const EPS: f64 = 1e-3;
    let lam_in = [
        (1.0 - EPS) * lam[0] + EPS * 0.25,
        (1.0 - EPS) * lam[1] + EPS * 0.25,
        (1.0 - EPS) * lam[2] + EPS * 0.25,
        (1.0 - EPS) * lam[3] + EPS * 0.25,
    ];
    let (g_in, det_in) = ct.jacobian(&lam_in);
    if det_in.is_finite() && det_in.abs() > 1e-30 && g_in.iter().flatten().all(|v| v.is_finite()) {
        g_in
    } else {
        *gs
    }
}

/// Adds the contribution of one quadrature point (weight `wq`, basis values `w` and curls `cu`) to the element matrix: the curl-curl
/// matrix for `which == 0`, the mass matrix otherwise. The same arithmetic serves straight and curved cells.
#[allow(clippy::too_many_arguments)]
fn add_point(which: usize, nb: usize, ke: &mut [C], w: &[[f64; 3]], cu: &[[f64; 3]], wq: f64, eps: &C3, mu_inv: &C3) {
    add_maxwell_point_c3(which, nb, ke, w, cu, wq, eps, mu_inv);
}

/// The curl-curl matrix `K` (symmetric `3x3` `mu_inv` tensor) and the mass matrix `M` (`eps` tensor) over all dofs.
pub fn assemble_maxwell(mesh: &Discretisation, ned: &Ned3, eps: &[C3], mu_inv: &[C3]) -> (CCsr, CCsr) {
    assemble_maxwell_with_pml(mesh, ned, eps, mu_inv, None)
}

/// Like [`assemble_maxwell`], evaluating curvilinear PML `pml_jacobian` + `pullback_maxwell` at each Gauss point inside PML cells.
pub fn assemble_maxwell_with_pml(
    mesh: &Discretisation,
    ned: &Ned3,
    eps: &[C3],
    mu_inv: &[C3],
    pml: Option<&EmPmlTable>,
) -> (CCsr, CCsr) {
    let (k, m, _) = assemble_maxwell_with_source(mesh, ned, eps, mu_inv, pml, &femcore::vvuq::UniformSource);
    (k, m)
}

/// Zero-overhead monomorphized 3D Maxwell assembly supporting custom volumetric `SourceField<3>` forcing.
///
/// When `S = femcore::vvuq::UniformSource`, `source.is_uniform()` is a compile-time constant `true`, so no RHS
/// load vector or non-PML physical quadrature coordinates are evaluated.
pub fn assemble_maxwell_with_source<S: femcore::vvuq::SourceField<3>>(
    mesh: &Discretisation,
    ned: &Ned3,
    eps: &[C3],
    mu_inv: &[C3],
    pml: Option<&EmPmlTable>,
    source: &S,
) -> (CCsr, CCsr, Vec<C>) {
    let geo = geometry(mesh);
    let quad = quad_rule(3, 4);
    // cells with a curved edge (ADR 82) are integrated with the geometry of their ten nodes; every other cell keeps the affine path.
    // On Barsoum quarter-point (`has_qp()`, FR-32) meshes, use degree-8 collapsed Gauss-Jacobi quadrature matching `apply_quarter_point_nodes`.
    let curved_cell: Vec<bool> = (0..mesh.n_cells()).map(|ci| mesh.cell_is_curved(ci)).collect();
    let cquad = quad_rule_collapsed(if mesh.has_qp() { 8 } else { CURVED_QUAD_DEGREE });
    let cshape2: Vec<(Vec<f64>, Vec<f64>)> = cquad.lam.iter().map(|l| shape(&[l[0], l[1], l[2], l[3]], 2)).collect();
    let nb = ned.nb();
    let n = ned.n_dofs();
    let nc = mesh.n_cells();
    let run = |which: usize| -> CCsr {
        assemble_complex(nc, n, |range, coo| {
            let (mut w, mut cu) = (Vec::new(), Vec::new());
            let mut ke = vec![c(0.0, 0.0); nb * nb];
            let mut gs = [[0.0; 3]; 4];
            for ci in range {
                let gid = ned.key_pairs(mesh, ci);
                for b in 0..4 {
                    for k in 0..3 {
                        gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                    }
                }
                let is_pml = pml.is_some_and(|p| p.is_pml(ci));
                ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
                if curved_cell[ci] {
                    let nodes10 = mesh.cell_geometry_nodes(ci);
                    let ct = CurvedTet::new(nodes10);
                    for q in 0..cquad.w.len() {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        let (gq, det) = ct.jacobian(&lam);
                        basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                        let wq = cquad.w[q] * det.abs() / 6.0;
                        let (ep, mui) = if let Some(p) = pml.filter(|_| is_pml) {
                            let nv10 = &cshape2[q].0;
                            let mut xq = [0.0; 3];
                            for a in 0..10 {
                                for k in 0..3 {
                                    xq[k] += nv10[a] * nodes10[a][k];
                                }
                            }
                            p.eval_3d(ci, eps[ci], mu_inv[ci], xq)
                        } else {
                            (eps[ci], mu_inv[ci])
                        };
                        add_point(which, nb, &mut ke, &w, &cu, wq, &ep, &mui);
                    }
                } else if let Some(p) = pml.filter(|_| is_pml) {
                    let vs = mesh.cell(ci);
                    let pts = [
                        mesh.point(vs[0] as usize),
                        mesh.point(vs[1] as usize),
                        mesh.point(vs[2] as usize),
                        mesh.point(vs[3] as usize),
                    ];
                    for q in 0..cquad.w.len() {
                        let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                        basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                        let wq = cquad.w[q] * geo.vol[ci];
                        let mut xq = [0.0; 3];
                        for b in 0..4 {
                            for k in 0..3 {
                                xq[k] += lam[b] * pts[b][k];
                            }
                        }
                        let (ep, mui) = p.eval_3d(ci, eps[ci], mu_inv[ci], xq);
                        add_point(which, nb, &mut ke, &w, &cu, wq, &ep, &mui);
                    }
                } else {
                    for q in 0..quad.w.len() {
                        let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                        basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                        let wq = quad.w[q] * geo.vol[ci];
                        add_point(which, nb, &mut ke, &w, &cu, wq, &eps[ci], &mu_inv[ci]);
                    }
                }
                coo.push_block(&ned.cell_dofs[ci], &ned.cell_dofs[ci], &ke);
            }
        })
    };
    let rhs = if source.is_uniform() {
        Vec::new()
    } else {
        let mut f_vol = vec![c(0.0, 0.0); n];
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];
        let quad_vol = quad_rule(3, 6);
        for ci in 0..nc {
            let gid = ned.key_pairs(mesh, ci);
            let dom = mesh.cell_labels[ci] as u32;
            for b in 0..4 {
                for k in 0..3 {
                    gs[b][k] = geo.g[ci * 12 + b * 3 + k];
                }
            }
            if curved_cell[ci] {
                let nodes10 = mesh.cell_geometry_nodes(ci);
                let ct = CurvedTet::new(nodes10);
                for q in 0..cquad.w.len() {
                    let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                    let (gq, det) = ct.jacobian(&lam);
                    basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                    let wq = cquad.w[q] * det.abs() / 6.0;
                    let xq = ct.position(&lam);
                    let fre = source.volume_rhs_re(xq, 0.0, dom);
                    let fim = source.volume_rhs_im(xq, 0.0, dom);
                    for s in 0..nb {
                        let dot_r = w[s][0] * fre[0] + w[s][1] * fre[1] + w[s][2] * fre[2];
                        let dot_i = w[s][0] * fim[0] + w[s][1] * fim[1] + w[s][2] * fim[2];
                        f_vol[ned.cell_dofs[ci][s]] += c(dot_r * wq, dot_i * wq);
                    }
                }
            } else {
                let vs = mesh.cell(ci);
                let pts = [
                    mesh.point(vs[0] as usize),
                    mesh.point(vs[1] as usize),
                    mesh.point(vs[2] as usize),
                    mesh.point(vs[3] as usize),
                ];
                for q in 0..quad_vol.w.len() {
                    let lam = [quad_vol.lam[q][0], quad_vol.lam[q][1], quad_vol.lam[q][2], quad_vol.lam[q][3]];
                    basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                    let wq = quad_vol.w[q] * geo.vol[ci];
                    let mut xq = [0.0; 3];
                    for b in 0..4 {
                        for k in 0..3 {
                            xq[k] += lam[b] * pts[b][k];
                        }
                    }
                    let fre = source.volume_rhs_re(xq, 0.0, dom);
                    let fim = source.volume_rhs_im(xq, 0.0, dom);
                    for s in 0..nb {
                        let dot_r = w[s][0] * fre[0] + w[s][1] * fre[1] + w[s][2] * fre[2];
                        let dot_i = w[s][0] * fim[0] + w[s][1] * fim[1] + w[s][2] * fim[2];
                        f_vol[ned.cell_dofs[ci][s]] += c(dot_r * wq, dot_i * wq);
                    }
                }
            }
        }
        f_vol
    };
    (run(0), run(1), rhs)
}

/// Solution and error norms of the 3D Maxwell curl-curl MMS BVP solver (`solve_maxwell_mms`).
pub struct MaxwellMmsResult {
    pub ned: Ned3,
    pub u: Vec<C>,
    pub l2_rel_error: f64,
    pub curl_rel_error: f64,
    pub hcurl_rel_error: f64,
}

/// Solves `(curl(mu_inv curl E) - k0^2 eps E = f)` on a 3D tetrahedral mesh (straight or `CurvedTet`)
/// using monomorphized `S: femcore::vvuq::SourceField<3>` for volume forcing, non-homogeneous tangential
/// Dirichlet data (`dirichlet_labels`), and optional Robin/Neumann trace data (`robin_labels`).
#[allow(clippy::too_many_arguments)]
pub fn solve_maxwell_mms<S: femcore::vvuq::SourceField<3>>(
    mesh: &Discretisation,
    eps: &[C3],
    mu_inv: &[C3],
    k0: f64,
    order: usize,
    dirichlet_labels: &[i32],
    robin_labels: &[i32],
    robin_alpha: C,
    source: &S,
    exact_e: &(dyn Fn([f64; 3]) -> [C; 3] + Sync),
    exact_curl: &(dyn Fn([f64; 3]) -> [C; 3] + Sync),
) -> Result<MaxwellMmsResult, String> {
    if mesh.dim != 3 {
        return Err("solve_maxwell_mms requires a 3D mesh".into());
    }
    let vkey: Vec<u32> = (0..mesh.n_points() as u32).collect();
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let (k_mat, m_mat, mut rhs) = assemble_maxwell_with_source(mesh, &ned, eps, mu_inv, None, source);
    if rhs.is_empty() {
        rhs.resize(n, c(0.0, 0.0));
    }
    let mut a_mat = k_mat.add_scaled(&m_mat, c(-(k0 * k0), 0.0));

    let face_to_cell = build_face_to_cell_3d(mesh);
    let geo = geometry(mesh);
    let quad2 = quad_rule(2, 6);

    // Helper evaluating facet quadrature points `(xq, wq, n_q, d_dofs, psi_q)` with `CurvedTet` support.
    let eval_facet = |tri_sorted: [u32; 3]| -> Result<(Vec<[f64; 3]>, Vec<f64>, Vec<[f64; 3]>, Vec<usize>, Vec<Vec<[f64; 3]>>), String> {
        let &(ci, lf_raw) = face_to_cell
            .get(&tri_sorted)
            .ok_or_else(|| format!("MMS boundary facet {:?} has no adjacent cell", tri_sorted))?;
        let cv = mesh.cell(ci);
        let mut lf = [0usize; 3];
        for (k, &tv) in tri_sorted.iter().enumerate() {
            lf[k] = cv.iter().position(|&v| v == tv).unwrap();
        }
        let opp = 6 - lf[0] - lf[1] - lf[2];
        let pt = |idx: usize| -> [f64; 3] {
            let p = mesh.point(cv[idx] as usize);
            [p[0], p[1], p[2]]
        };
        let (p0, p1, p2, p_opp) = (pt(lf[0]), pt(lf[1]), pt(lf[2]), pt(opp));
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let cross = |a: [f64; 3], b: [f64; 3]| {
            [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
        };
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let norm = |a: [f64; 3]| dot(a, a).sqrt();
        let mut n_raw = cross(sub(p1, p0), sub(p2, p0));
        if dot(n_raw, sub(p0, p_opp)) < 0.0 {
            n_raw = [-n_raw[0], -n_raw[1], -n_raw[2]];
        }
        let area_flat = 0.5 * norm(n_raw);
        let n_flat = [n_raw[0] / norm(n_raw), n_raw[1] / norm(n_raw), n_raw[2] / norm(n_raw)];

        let bnd_dofs = ned.boundary_dofs(mesh, tri_sorted);
        let cdofs = &ned.cell_dofs[ci];
        let local_slots: Vec<usize> = bnd_dofs.iter().map(|&gd| cdofs.iter().position(|&cd| cd == gd).unwrap()).collect();
        let gid = ned.key_pairs(mesh, ci);
        let mut gs = [[0.0; 3]; 4];
        for b in 0..4 {
            for kk in 0..3 {
                gs[b][kk] = geo.g[ci * 12 + b * 3 + kk];
            }
        }
        let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let (mut xqs, mut wqs, mut nqs, mut psis) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for q in 0..quad2.w.len() {
            let l2 = &quad2.lam[q];
            let mut lam4 = [0.0; 4];
            lam4[lf[0]] = l2[0];
            lam4[lf[1]] = l2[1];
            lam4[lf[2]] = l2[2];
            let (xq, area_q, n_q) = match &curved {
                Some(ct) => {
                    let gq = safe_curved_jacobian(ct, &lam4, &gs);
                    basis3_keyed(order, &lam4, &gq, &gid, &mut w, &mut cu);
                    let d = ct.dx_dlam(&lam4);
                    let nr = cross(sub(d[lf[1]], d[lf[0]]), sub(d[lf[2]], d[lf[0]]));
                    let nl = norm(nr);
                    let mut nq = if nl > 1e-30 { [nr[0] / nl, nr[1] / nl, nr[2] / nl] } else { n_flat };
                    if dot(nq, n_flat) < 0.0 {
                        nq = [-nq[0], -nq[1], -nq[2]];
                    }
                    (ct.position(&lam4), if nl > 1e-30 { 0.5 * nl } else { area_flat }, nq)
                }
                None => {
                    basis3_keyed(order, &lam4, &gs, &gid, &mut w, &mut cu);
                    let xp = [
                        l2[0] * p0[0] + l2[1] * p1[0] + l2[2] * p2[0],
                        l2[0] * p0[1] + l2[1] * p1[1] + l2[2] * p2[1],
                        l2[0] * p0[2] + l2[1] * p1[2] + l2[2] * p2[2],
                    ];
                    (xp, area_flat, n_flat)
                }
            };
            let mut row = Vec::with_capacity(local_slots.len());
            for &s in &local_slots {
                let ws = w[s];
                let wn = dot(ws, n_q);
                row.push([ws[0] - wn * n_q[0], ws[1] - wn * n_q[1], ws[2] - wn * n_q[2]]);
            }
            xqs.push(xq);
            wqs.push(quad2.w[q] * area_q);
            nqs.push(n_q);
            psis.push(row);
        }
        let _ = lf_raw;
        Ok((xqs, wqs, nqs, bnd_dofs, psis))
    };

    // Assemble optional Neumann/Robin boundary terms
    if !robin_labels.is_empty() {
        let mut r_coo = CCoo::new(n, n);
        let mut seen_r = std::collections::BTreeSet::new();
        let coeff_r = c(0.0, k0) * robin_alpha;
        for f in 0..mesh.n_facets() {
            let label = mesh.facet_labels[f];
            if !robin_labels.contains(&label) {
                continue;
            }
            let fv = mesh.facet(f);
            let key = crate::ibc::sorted3(fv[0], fv[1], fv[2]);
            if !seen_r.insert(key) {
                continue;
            }
            let (xqs, wqs, nqs, bnd_dofs, psis) = eval_facet(key)?;
            let nb_f = bnd_dofs.len();
            let mut loc_r = vec![c(0.0, 0.0); nb_f * nb_f];
            for q in 0..wqs.len() {
                let wq = wqs[q];
                let xq = xqs[q];
                let n_q = nqs[q];
                let pq = &psis[q];
                if !source.is_uniform() {
                    let gre = source.neumann_robin_re(xq, n_q, 0.0, label as u32);
                    let gim = source.neumann_robin_im(xq, n_q, 0.0, label as u32);
                    for a in 0..nb_f {
                        let dot_r = pq[a][0] * gre[0] + pq[a][1] * gre[1] + pq[a][2] * gre[2];
                        let dot_i = pq[a][0] * gim[0] + pq[a][1] * gim[1] + pq[a][2] * gim[2];
                        rhs[bnd_dofs[a]] += c(dot_r * wq, dot_i * wq);
                    }
                }
                if coeff_r.norm() > 0.0 {
                    for a in 0..nb_f {
                        for b in 0..nb_f {
                            let dp = pq[a][0] * pq[b][0] + pq[a][1] * pq[b][1] + pq[a][2] * pq[b][2];
                            loc_r[a * nb_f + b] += coeff_r * (dp * wq);
                        }
                    }
                }
            }
            if coeff_r.norm() > 0.0 {
                r_coo.push_block(&bnd_dofs, &bnd_dofs, &loc_r);
            }
        }
        let r_mat = r_coo.to_csr();
        if r_mat.nnz() > 0 {
            a_mat = a_mat.add_scaled(&r_mat, c(1.0, 0.0));
        }
    }

    // Collect Dirichlet boundary DOFs using canonical edge-moment + per-face 2x2 bubble interpolation.
    // Unlike a global L^2(Gamma_D) Gram projection, decoupling edges from faces prevents O(h^{3/2})
    // corner contamination across perpendicular faces while remaining exact for P1/P2 Nedelec on CurvedTet.
    let mut fixed_set = std::collections::BTreeSet::new();
    let mut dir_tris = Vec::new();
    let mut dir_edges = std::collections::BTreeMap::new();
    let mut seen_d = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        if !dirichlet_labels.contains(&label) {
            continue;
        }
        let fv = mesh.facet(f);
        let key = crate::ibc::sorted3(fv[0], fv[1], fv[2]);
        if seen_d.insert(key) {
            dir_tris.push((key, label));
            for (ea, eb) in [(key[0], key[1]), (key[1], key[2]), (key[0], key[2])] {
                let eid = mesh.edge_index(ea, eb);
                dir_edges.entry(eid).or_insert(label);
            }
            for gd in ned.boundary_dofs(mesh, key) {
                fixed_set.insert(gd);
            }
        }
    }
    let fixed: Vec<usize> = fixed_set.into_iter().collect();
    let mut u_dir = vec![c(0.0, 0.0); fixed.len()];
    if !fixed.is_empty() && !source.is_uniform() {
        let mut idx_map = HashMap::with_capacity(fixed.len());
        for (li, &gd) in fixed.iter().enumerate() {
            idx_map.insert(gd, li);
        }
        // 5-point Gauss-Legendre on [0, 1] for exact canonical edge tangential moments
        const GL5_T: [f64; 5] =
            [0.046910077030668004, 0.23076534494715845, 0.5, 0.7692346550528415, 0.9530899229693319];
        const GL5_W: [f64; 5] =
            [0.11846344252809454, 0.23931433524968324, 0.28444444444444444, 0.23931433524968324, 0.11846344252809454];
        for (eid, label) in dir_edges {
            let (v0, v1) = mesh.edges()[eid];
            let (a, b) = if (ned.vkey[v0 as usize], v0) < (ned.vkey[v1 as usize], v1) { (v0, v1) } else { (v1, v0) };
            let pa = mesh.point(a as usize);
            let pb = mesh.point(b as usize);
            let pm =
                mesh.edge_node(eid).unwrap_or([0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])]);
            let (mut m0, mut m1) = (c(0.0, 0.0), c(0.0, 0.0));
            for q in 0..5 {
                let t = GL5_T[q];
                let wq = GL5_W[q];
                let s0 = (1.0 - t) * (1.0 - 2.0 * t);
                let s1 = t * (2.0 * t - 1.0);
                let s2 = 4.0 * t * (1.0 - t);
                let d0 = 4.0 * t - 3.0;
                let d1 = 4.0 * t - 1.0;
                let d2 = 4.0 * (1.0 - 2.0 * t);
                let xq = [
                    s0 * pa[0] + s1 * pb[0] + s2 * pm[0],
                    s0 * pa[1] + s1 * pb[1] + s2 * pm[1],
                    s0 * pa[2] + s1 * pb[2] + s2 * pm[2],
                ];
                let xp = [
                    d0 * pa[0] + d1 * pb[0] + d2 * pm[0],
                    d0 * pa[1] + d1 * pb[1] + d2 * pm[1],
                    d0 * pa[2] + d1 * pb[2] + d2 * pm[2],
                ];
                let ere = source.dirichlet_re(xq, 0.0, label as u32);
                let eim = source.dirichlet_im(xq, 0.0, label as u32);
                let et = c(
                    ere[0] * xp[0] + ere[1] * xp[1] + ere[2] * xp[2],
                    eim[0] * xp[0] + eim[1] * xp[1] + eim[2] * xp[2],
                );
                m0 += et * wq;
                m1 += et * ((1.0 - 2.0 * t) * wq);
            }
            if order == 1 {
                u_dir[idx_map[&eid]] = m0;
            } else {
                u_dir[idx_map[&(2 * eid)]] = m0;
                u_dir[idx_map[&(2 * eid + 1)]] = m1 * 3.0;
            }
        }
        if order == 2 {
            // Canonical De Rham [P_0(f)]^2 face-moment projection in covariant Piola coordinates:
            // int_{\hat{f}} (w - sum u_edge w_edge) . (dx/dxi_a) d\hat{A} = int_{\hat{f}} (E* - sum u_edge w_edge) . (dx/dxi_a) d\hat{A}
            // commutes exactly with curl (curl \Pi_h^{curl} = \Pi_h^{div} curl) on both straight and CurvedTet meshes.
            let mut w_f = Vec::new();
            let mut cu_f = Vec::new();
            let mut gs_f = [[0.0; 3]; 4];
            for (key, label) in dir_tris {
                let &(ci, _) = face_to_cell
                    .get(&key)
                    .ok_or_else(|| format!("MMS boundary facet {:?} has no adjacent cell", key))?;
                let cv = mesh.cell(ci);
                let mut lf = [0usize; 3];
                for (k, &tv) in key.iter().enumerate() {
                    lf[k] = cv.iter().position(|&v| v == tv).unwrap();
                }
                let p0 = mesh.point(cv[lf[0]] as usize);
                let p1 = mesh.point(cv[lf[1]] as usize);
                let p2 = mesh.point(cv[lf[2]] as usize);
                let t0_flat = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                let t1_flat = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                let bnd_dofs = ned.boundary_dofs(mesh, key);
                let cdofs = &ned.cell_dofs[ci];
                let local_slots: Vec<usize> =
                    bnd_dofs.iter().map(|&gd| cdofs.iter().position(|&cd| cd == gd).unwrap()).collect();
                let gid = ned.key_pairs(mesh, ci);
                for b in 0..4 {
                    for kk in 0..3 {
                        gs_f[b][kk] = geo.g[ci * 12 + b * 3 + kk];
                    }
                }
                let curved = mesh.cell_is_curved(ci).then(|| CurvedTet::new(mesh.cell_geometry_nodes(ci)));
                let mut m = [[0.0_f64; 2]; 2];
                let mut r = [c(0.0, 0.0); 2];
                for q in 0..quad2.w.len() {
                    let l2 = &quad2.lam[q];
                    let wq = quad2.w[q];
                    let mut lam4 = [0.0; 4];
                    lam4[lf[0]] = l2[0];
                    lam4[lf[1]] = l2[1];
                    lam4[lf[2]] = l2[2];
                    let (xq, t_vecs) = match &curved {
                        Some(ct) => {
                            let gq = safe_curved_jacobian(ct, &lam4, &gs_f);
                            basis3_keyed(order, &lam4, &gq, &gid, &mut w_f, &mut cu_f);
                            let d = ct.dx_dlam(&lam4);
                            let t0 = [d[lf[1]][0] - d[lf[0]][0], d[lf[1]][1] - d[lf[0]][1], d[lf[1]][2] - d[lf[0]][2]];
                            let t1 = [d[lf[2]][0] - d[lf[0]][0], d[lf[2]][1] - d[lf[0]][1], d[lf[2]][2] - d[lf[0]][2]];
                            (ct.position(&lam4), [t0, t1])
                        }
                        None => {
                            basis3_keyed(order, &lam4, &gs_f, &gid, &mut w_f, &mut cu_f);
                            let xp = [
                                l2[0] * p0[0] + l2[1] * p1[0] + l2[2] * p2[0],
                                l2[0] * p0[1] + l2[1] * p1[1] + l2[2] * p2[1],
                                l2[0] * p0[2] + l2[1] * p1[2] + l2[2] * p2[2],
                            ];
                            (xp, [t0_flat, t1_flat])
                        }
                    };
                    let ere = source.dirichlet_re(xq, 0.0, label as u32);
                    let eim = source.dirichlet_im(xq, 0.0, label as u32);
                    let mut res = [c(ere[0], eim[0]), c(ere[1], eim[1]), c(ere[2], eim[2])];
                    for a in 0..6 {
                        let ua = u_dir[idx_map[&bnd_dofs[a]]];
                        let wa = w_f[local_slots[a]];
                        for k in 0..3 {
                            res[k] -= ua * wa[k];
                        }
                    }
                    for a in 0..2 {
                        let ta = t_vecs[a];
                        let d = c(
                            ta[0] * res[0].re + ta[1] * res[1].re + ta[2] * res[2].re,
                            ta[0] * res[0].im + ta[1] * res[1].im + ta[2] * res[2].im,
                        );
                        r[a] += d * wq;
                        for b in 0..2 {
                            let wb = w_f[local_slots[6 + b]];
                            m[a][b] += wq * (ta[0] * wb[0] + ta[1] * wb[1] + ta[2] * wb[2]);
                        }
                    }
                }
                let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
                let inv_det = 1.0 / det;
                u_dir[idx_map[&bnd_dofs[6]]] = (r[0] * m[1][1] - r[1] * m[0][1]) * inv_det;
                u_dir[idx_map[&bnd_dofs[7]]] = (r[1] * m[0][0] - r[0] * m[1][0]) * inv_det;
            }
        }
    }

    let free = free_dofs(n, &fixed);
    let mut u = vec![c(0.0, 0.0); n];
    for (&d, &val) in fixed.iter().zip(&u_dir) {
        u[d] = val;
    }
    if !free.is_empty() {
        let a_ff = a_mat.submatrix(&free, &free);
        let a_fd = a_mat.submatrix(&free, &fixed);
        let ad = a_fd.matvec(&u_dir);
        let f_free: Vec<C> = free.iter().enumerate().map(|(i, &gd)| rhs[gd] - ad[i]).collect();
        let solver = femcore::linalg::driver::DirectSolver::new(&a_ff, false)?;
        let u_free = solver.solve(&f_free);
        for (i, &gd) in free.iter().enumerate() {
            u[gd] = u_free[i];
        }
    }

    let (l2_rel_error, curl_rel_error, hcurl_rel_error) = maxwell_field_errors(mesh, &ned, &u, exact_e, exact_curl);
    Ok(MaxwellMmsResult { ned, u, l2_rel_error, curl_rel_error, hcurl_rel_error })
}

/// Evaluates continuous relative `L^2`, curl semi-norm, and `H(curl)` errors on curved or straight 3D tets.
pub fn maxwell_field_errors(
    mesh: &Discretisation,
    ned: &Ned3,
    u: &[C],
    exact_e: &(dyn Fn([f64; 3]) -> [C; 3] + Sync),
    exact_curl: &(dyn Fn([f64; 3]) -> [C; 3] + Sync),
) -> (f64, f64, f64) {
    let geo = geometry(mesh);
    let cquad = quad_rule_collapsed(CURVED_QUAD_DEGREE);
    let quad = quad_rule(3, 6);
    let nb = ned.nb();
    let (mut err_e2, mut ref_e2, mut err_c2, mut ref_c2) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];
    for ci in 0..mesh.n_cells() {
        let gid = ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for k in 0..3 {
                gs[b][k] = geo.g[ci * 12 + b * 3 + k];
            }
        }
        let cdofs = &ned.cell_dofs[ci];
        let mut accum_pt = |xq: [f64; 3], wq: f64, w_b: &[[f64; 3]], cu_b: &[[f64; 3]]| {
            let (mut eh, mut ch) = ([c(0.0, 0.0); 3], [c(0.0, 0.0); 3]);
            for s in 0..nb {
                let cf = u[cdofs[s]];
                for k in 0..3 {
                    eh[k] += cf * w_b[s][k];
                    ch[k] += cf * cu_b[s][k];
                }
            }
            let ee = exact_e(xq);
            let ec = exact_curl(xq);
            for k in 0..3 {
                err_e2 += (eh[k] - ee[k]).norm_sqr() * wq;
                ref_e2 += ee[k].norm_sqr() * wq;
                err_c2 += (ch[k] - ec[k]).norm_sqr() * wq;
                ref_c2 += ec[k].norm_sqr() * wq;
            }
        };
        if mesh.cell_is_curved(ci) {
            let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
            for q in 0..cquad.w.len() {
                let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                let (gq, det) = ct.jacobian(&lam);
                basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                accum_pt(ct.position(&lam), cquad.w[q] * det.abs() / 6.0, &w, &cu);
            }
        } else {
            let vs = mesh.cell(ci);
            let pts = [
                mesh.point(vs[0] as usize),
                mesh.point(vs[1] as usize),
                mesh.point(vs[2] as usize),
                mesh.point(vs[3] as usize),
            ];
            for q in 0..quad.w.len() {
                let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                let xq = [
                    (0..4).map(|b| lam[b] * pts[b][0]).sum(),
                    (0..4).map(|b| lam[b] * pts[b][1]).sum(),
                    (0..4).map(|b| lam[b] * pts[b][2]).sum(),
                ];
                accum_pt(xq, quad.w[q] * geo.vol[ci], &w, &cu);
            }
        }
    }
    let l2 = (err_e2 / ref_e2.max(1e-300)).sqrt();
    let curl = (err_c2 / ref_c2.max(1e-300)).sqrt();
    let hcurl = ((err_e2 + err_c2) / (ref_e2 + ref_c2).max(1e-300)).sqrt();
    (l2, curl, hcurl)
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
        let (lo, hi) = if (ned.vkey[a as usize], a) < (ned.vkey[b as usize], b) { (a, b) } else { (b, a) };
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
pub(crate) fn vertex_keys(mesh: &Discretisation, matches: &[femcore::periodic::Match]) -> Vec<u32> {
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
pub(crate) fn ned_relations(
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
    pub kc_no_damp: CCsr,
    pub mc: CCsr,
    pub g_red: CCsr,
    pub g_h: CCsr,
    pub f0: f64,
    pub sigma: f64,
    pub n_modes: usize,
    pub f_search: f64,
    pub has_cqed: bool,
    pub k0_red: Csr,
    pub k_phys_red: Csr,
    pub k_ibc_red: Csr,
    pub c_ibc_red: Csr,
    pub c_port_red: Csr,
    pub c_damp_red: Csr,
    pub m0_full: Csr,
    pub m0_red: Csr,
    pub m_phys_red: Csr,
    pub jj_ports: Vec<AssembledJjPort>,
    pub lumped_ports: Vec<AssembledLumpedPort>,
    pub pert_ibc_tris: Vec<PerturbativeIbcTri>,
    pub ibc_warnings: Vec<String>,
    pub cell_sheet: Vec<Option<SheetCellKind>>,
    pub periodic_node_relations: Vec<(usize, usize, C)>,
    pub mu_inv: Vec<C3>,
    pub pml_cells: Vec<bool>,
    pub eps_phys_re: Vec<[[f64; 3]; 3]>,
    pub d_omega_eps: Vec<[[f64; 3]; 3]>,
    pub mur_phys: Vec<f64>,
    pub y_periodic_ky: Option<f64>,
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
    validate_port_contacts_3d(mesh, problem)?;
    let order = if problem.order == 1 { 1 } else { 2 };
    let tens = problem.materials.cell_complex_permittivity(mesh, Some(1.0))?;
    let eps_phys_re: Vec<[[f64; 3]; 3]> = tens.iter().map(|(re, _)| *re).collect();
    let d_omega_eps = problem.cell_d_omega_eps(mesh)?;
    let mur = problem.materials.cell_values(mesh, "mu_r", Some(1.0))?;
    let (eps, mu_inv) = simplex_cell_materials_3d(mesh, problem)?;
    let (lo, hi) = mesh.bbox();
    let size = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(1e-300);
    let f0 = if f_search > 0.0 { f_search } else { C0 / (2.0 * size) };
    let pml_table = EmPmlTable::from_simplex(mesh, &problem.materials, f0, true, false);
    let pml_cells = pml_table.pml_cells_mask();
    let mut tl = Instant::now();
    // periodic pairs
    let mut matches = vec![];
    let mut y_periodic_ky = None;
    for p in &problem.periodic {
        let mt = femcore::periodic::match_vertices(mesh, p)?;
        let sn = (mt.shift[0] * mt.shift[0] + mt.shift[1] * mt.shift[1] + mt.shift[2] * mt.shift[2]).sqrt();
        if sn > 0.0 && mt.shift[1].abs() > 0.99 * sn && y_periodic_ky.is_none() {
            y_periodic_ky = Some(p.phase / mt.shift[1]);
        }
        matches.push(mt);
    }
    let vkey = vertex_keys(mesh, &matches);
    let ned = Ned3::new(mesh, order, vkey);
    let n = ned.n_dofs();
    let (k, m) = assemble_maxwell_with_pml(mesh, &ned, &eps, &mu_inv, pml_table.has_any().then_some(&pml_table));
    lap(&mut tl, "emw3d: numbering and assembly");

    let sigma = (2.0 * PI * f0 / C0).powi(2);

    let has_cqed =
        !problem.surface_impedance.is_empty() || !problem.lumped_jj.is_empty() || !problem.lumped_ports.is_empty();
    let has_sheets = has_cqed || !problem.pec.is_empty();
    let mut ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d_cavity(
            mesh,
            &ned,
            &problem.pec,
            &problem.surface_impedance,
            &problem.lumped_jj,
            f0,
        )?)
    } else {
        None
    };
    let mut lumped_ports = assemble_lumped_ports_3d(mesh, &ned, &problem.pec, &problem.lumped_ports)?;

    // perfect electric conductor on the exterior boundary, except the periodic boundaries, magnetic walls (`Problem::pmc`, natural),
    // direct `surface_impedance` sheets, `lumped_jj` sheets, and `lumped_port` sheets.
    // Perturbative `surface_impedance` sheets are treated as strong PEC during Stage 1/2 (FR-30).
    let exterior = mesh.exterior_facets();
    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let direct_ibc_labels: Vec<i32> = ibc_jj.as_ref().map_or_else(Vec::new, |ij| ij.direct_ibc_labels.clone());
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();

    let (sym_pec_facets, sym_pmc_facets) = symmetry_facet_flags(mesh, problem.symmetry)?;
    let mut fixed_set = std::collections::BTreeSet::new();
    let mut ext_pec_vertex = std::collections::BTreeSet::new();
    let mut pec_edge = std::collections::BTreeSet::new();
    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        let is_ext_pec = exterior[f]
            && !periodic_labels.contains(&label)
            && (sym_pec_facets[f]
                || (!problem.pmc.contains(&label)
                    && !sym_pmc_facets[f]
                    && !direct_ibc_labels.contains(&label)
                    && !jj_labels.contains(&label)
                    && !lp_labels.contains(&label)));
        if is_ext_pec {
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
    let need_gt_r = !lumped_ports.is_empty()
        || ibc_jj.as_ref().is_some_and(|ij| ij.jj_ports.iter().any(|p| p.bc.formulation != "distributed_sheet"));
    let gt_r_opt = if need_gt_r { Some(gradient_full(mesh, &ned).transpose()) } else { None };
    let nv_mesh = mesh.n_points();
    if let Some(ref ij) = ibc_jj {
        // Explicit `pec` facets and perturbative `surface_impedance` facets (FR-30): strong E_|| = 0 elimination in `fixed_set`
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
        for pt in &ij.pert_ibc_tris {
            let fv = pt.tri;
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
        // Direct `surface_impedance` facets: keep Nedelec DOFs FREE in `red`
        for &fv in &ij.ibc_tris {
            for v in fv {
                cond_vertices.insert(v as usize);
            }
            for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                pec_edge.insert(mesh.edge_index(a, b));
                cond_dsu_edges.push((a as usize, b as usize));
            }
        }
        // Distributed-sheet `lumped_jj` facets carry surface admittance across all sheet vertices/edges.
        for port in &ij.jj_ports {
            if port.bc.formulation == "distributed_sheet" {
                for &fv in &port.triangles {
                    for v in fv {
                        cond_vertices.insert(v as usize);
                    }
                    for (a, b) in [(fv[0], fv[1]), (fv[1], fv[2]), (fv[0], fv[2])] {
                        cond_dsu_edges.push((a as usize, b as usize));
                        pec_edge.insert(mesh.edge_index(a, b));
                    }
                }
            }
        }
    }
    // FR-26: For `rank1` `lumped_jj` and all `lumped_port` sheets, never insert interior or
    // lateral-vacuum vertices into `cond_vertices` / `cond_dsu_edges` (their gradient columns have
    // `(G^T b_j)_v == 0` and must stay in `g_red` so zero-frequency gradient modes are projected out).
    // Bridge only true electrode seam vertices (`|gt_b[v]| > 1e-8 * b_max`, exterior PEC, or conductor).
    let pec_ibc_verts = cond_vertices.clone();
    let register_rank1_port = |triangles: &[[u32; 3]],
                               b_full: &[f64],
                               gt_r: &Csr,
                               cond_vertices: &mut std::collections::BTreeSet<usize>,
                               cond_dsu_edges: &mut Vec<(usize, usize)>,
                               pec_edge: &mut std::collections::BTreeSet<usize>| {
        let gt_b = gt_r.matvec(b_full);
        let b_max = gt_b.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
        let mut port_verts = std::collections::BTreeSet::new();
        for &fv in triangles {
            for v in fv {
                port_verts.insert(v as usize);
            }
        }
        let mut term_verts: Vec<usize> = Vec::new();
        for v in port_verts {
            let is_term = (b_max > 0.0 && gt_b[v].abs() > 1e-8 * b_max)
                || ext_pec_vertex.contains(&v)
                || pec_ibc_verts.contains(&v)
                || cond_vertices.contains(&v);
            if is_term {
                term_verts.push(v);
            }
        }
        for &v in &term_verts {
            cond_vertices.insert(v);
        }
        for win in term_verts.windows(2) {
            cond_dsu_edges.push((win[0], win[1]));
        }
        if order == 2 && b_max > 0.0 {
            for e in 0..mesh.n_edges() {
                if gt_b[nv_mesh + e].abs() > 1e-8 * b_max {
                    pec_edge.insert(e);
                }
            }
        }
    };
    if let Some(ref gt_r) = gt_r_opt {
        if let Some(ref ij) = ibc_jj {
            for port in &ij.jj_ports {
                if port.bc.formulation != "distributed_sheet" {
                    register_rank1_port(
                        &port.triangles,
                        &port.b_full,
                        gt_r,
                        &mut cond_vertices,
                        &mut cond_dsu_edges,
                        &mut pec_edge,
                    );
                }
            }
        }
        for lp in &lumped_ports {
            register_rank1_port(
                &lp.triangles,
                &lp.b_full,
                gt_r,
                &mut cond_vertices,
                &mut cond_dsu_edges,
                &mut pec_edge,
            );
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
    let (k_phys_red, m_phys_red, m0_full) = if pml_cells.iter().any(|&p| p) {
        let mut eps_phys_c: Vec<C3> = Vec::with_capacity(mesh.n_cells());
        let mut mu_phys_c: Vec<C3> = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            if pml_cells[ci] {
                eps_phys_c.push([[c(0.0, 0.0); 3]; 3]);
                mu_phys_c.push([[c(0.0, 0.0); 3]; 3]);
            } else {
                let mut e = [[c(0.0, 0.0); 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        e[a][b] = c(eps_phys_re[ci][a][b], 0.0);
                    }
                }
                eps_phys_c.push(e);
                mu_phys_c.push(iso_mu_inv_c3(mur[ci]));
            }
        }
        let (kp, mp) = assemble_maxwell(mesh, &ned, &eps_phys_c, &mu_phys_c);
        let kp_r = to_real(&femcore::periodic::congruence(&pmat, &kp));
        let mp_r = to_real(&femcore::periodic::congruence(&pmat, &mp));
        let mf = if has_cqed { to_real(&mp) } else { Csr::zeros(0, 0) };
        (kp_r, mp_r, mf)
    } else {
        (k0_red.clone(), m0_red.clone(), if has_cqed { to_real(&m) } else { Csr::zeros(0, 0) })
    };

    let mut c_port_all_full = Csr::zeros(n, n);
    let mut c_port_active_full = Csr::zeros(n, n);
    for lp in &mut lumped_ports {
        lp.b_red = reduce_real_vec(&red, &lp.b_full);
        if lp.c_port_full.nnz() > 0 {
            c_port_all_full = c_port_all_full.add_scaled(&lp.c_port_full, 1.0);
            if lp.bc.complex_damping {
                c_port_active_full = c_port_active_full.add_scaled(&lp.c_port_full, 1.0);
            }
        }
    }
    let c_port_red = reduce_real_csr(&red, &pmat, &c_port_all_full);
    let c_port_active_red = reduce_real_csr(&red, &pmat, &c_port_active_full);

    let k0_f0 = 2.0 * PI * f0 / C0;
    let (kc, kc_no_damp, mc, k_ibc_red, c_ibc_red, c_damp_red, jj_ports, pert_ibc_tris, ibc_warnings, cell_sheet) =
        if let Some(ref mut ij) = ibc_jj {
            let mut k_tot_no_damp = k.clone();
            let mut k_tot = k.clone();
            let mut m_tot = m.clone();
            if ij.k_ibc_full.nnz() > 0 {
                k_tot_no_damp = k_tot_no_damp.add_scaled(&to_complex(&ij.k_ibc_full), c(1.0, 0.0));
            }
            if ij.kc_ibc_full.nnz() > 0 {
                k_tot = k_tot.add_scaled(&ij.kc_ibc_full, c(1.0, 0.0));
            }
            if ij.k_jj_full.nnz() > 0 {
                let kjj_c = to_complex(&ij.k_jj_full);
                k_tot_no_damp = k_tot_no_damp.add_scaled(&kjj_c, c(1.0, 0.0));
                k_tot = k_tot.add_scaled(&kjj_c, c(1.0, 0.0));
            }
            if c_port_active_full.nnz() > 0 {
                k_tot = k_tot.add_scaled(&to_complex(&c_port_active_full), c(0.0, k0_f0));
            }
            if ij.m_jj_full.nnz() > 0 {
                m_tot = m_tot.add_scaled(&to_complex(&ij.m_jj_full), c(1.0, 0.0));
            }
            let kc_no_damp =
                if has_cqed { femcore::periodic::congruence(&pmat, &k_tot_no_damp) } else { k0_c_red.clone() };
            let kc = if has_cqed { femcore::periodic::congruence(&pmat, &k_tot) } else { k0_c_red };
            let mc = if has_cqed { femcore::periodic::congruence(&pmat, &m_tot) } else { m0_c_red };
            let k_ibc_red = reduce_real_csr(&red, &pmat, &ij.k_ibc_full);
            let c_ibc_red = reduce_real_csr(&red, &pmat, &ij.c_ibc_full);
            let c_damp_red = if c_port_active_red.nnz() > 0 {
                c_ibc_red.add_scaled(&c_port_active_red, 1.0)
            } else {
                c_ibc_red.clone()
            };
            for port in &mut ij.jj_ports {
                port.b_red = reduce_real_vec(&red, &port.b_full);
                port.k_jj_red = reduce_real_csr(&red, &pmat, &port.k_jj_red);
                port.m_jj_red = reduce_real_csr(&red, &pmat, &port.m_jj_red);
            }
            (
                kc,
                kc_no_damp,
                mc,
                k_ibc_red,
                c_ibc_red,
                c_damp_red,
                std::mem::take(&mut ij.jj_ports),
                std::mem::take(&mut ij.pert_ibc_tris),
                std::mem::take(&mut ij.warnings),
                std::mem::take(&mut ij.cell_sheet),
            )
        } else {
            let kc_no_damp = k0_c_red.clone();
            let kc = if c_port_active_full.nnz() > 0 {
                let k_tot = k.add_scaled(&to_complex(&c_port_active_full), c(0.0, k0_f0));
                femcore::periodic::congruence(&pmat, &k_tot)
            } else {
                k0_c_red
            };
            (
                kc,
                kc_no_damp,
                m0_c_red,
                Csr::zeros(nred, nred),
                Csr::zeros(nred, nred),
                c_port_active_red,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                vec![None; mesh.n_cells()],
            )
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
    // When there are no PEC vertices and all periodic Floquet phases are 0 mod 2*pi (e.g., all-PMC walls or Gamma-point periodic cell),
    // the constant scalar potential 1_comp lies in ker(G); pin 1 vertex per connected component so G^H M G is non-singular without altering Range(G).
    if pec_vertex.is_empty()
        && nv > 0
        && problem.periodic.iter().all(|p| (femcore::periodic::phase_factor(p.phase) - c(1.0, 0.0)).norm() < 1e-9)
    {
        let mut comp_dsu: Vec<usize> = (0..nv).collect();
        let find_c = |p: &mut [usize], mut i: usize| -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        };
        for &(a, b) in mesh.edges() {
            let (ra, rb) = (find_c(&mut comp_dsu, a as usize), find_c(&mut comp_dsu, b as usize));
            if ra != rb {
                comp_dsu[ra.max(rb)] = ra.min(rb);
            }
        }
        for &(d, s, _) in &node_rel {
            if d < nv && s < nv {
                let (rd, rs) = (find_c(&mut comp_dsu, d), find_c(&mut comp_dsu, s));
                if rd != rs {
                    comp_dsu[rd.max(rs)] = rd.min(rs);
                }
            }
        }
        let mut pinned_roots = std::collections::BTreeSet::new();
        for v in 0..nv {
            let r = find_c(&mut comp_dsu, v);
            if pinned_roots.insert(r) && !node_fixed.contains(&v) {
                node_fixed.push(v);
            }
        }
        node_fixed.sort_unstable();
        node_fixed.dedup();
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
    if g_h.nrows() > 0 {
        let check_orth = |name: &str, b_red: &[f64]| -> Result<(), String> {
            let bnorm = b_red.iter().map(|&v| v * v).sum::<f64>().sqrt();
            if bnorm <= 1e-300 {
                return Ok(());
            }
            let b_c: Vec<C> = b_red.iter().map(|&v| c(v, 0.0)).collect();
            let gtb = g_h.matvec(&b_c);
            let err = gtb.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt() / bnorm;
            if err > 1e-10 {
                return Err(format!(
                    "cavity null-space orthogonality violation on port '{name}': ||G_red^T b_red||_2 / ||b_red||_2 = {err:.3e} > 1e-10"
                ));
            }
            Ok(())
        };
        for port in &jj_ports {
            if port.bc.formulation != "distributed_sheet" {
                check_orth(&port.bc.name, &port.b_red)?;
            }
        }
        for lp in &lumped_ports {
            check_orth(&lp.bc.name, &lp.b_red)?;
        }
    }
    lap(&mut tl, "emw3d: reduction and gradient matrix");

    Ok(EmwAssembly {
        order,
        n,
        nred,
        ned,
        red,
        kc,
        kc_no_damp,
        mc,
        g_red,
        g_h,
        f0,
        sigma,
        n_modes,
        f_search,
        has_cqed,
        k0_red,
        k_phys_red,
        k_ibc_red,
        c_ibc_red,
        c_port_red,
        c_damp_red,
        m0_full,
        m0_red,
        m_phys_red,
        jj_ports,
        lumped_ports,
        pert_ibc_tris,
        ibc_warnings,
        cell_sheet,
        periodic_node_relations,
        mu_inv,
        pml_cells,
        eps_phys_re,
        d_omega_eps,
        mur_phys: mur,
        y_periodic_ky,
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
        let n_extra_lp = asm.lumped_ports.len();
        let n_ask = (n_modes + n_extra_lp).min(asm.nred.saturating_sub(1)).max(n_modes);
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
                femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_ask, &opts, Some(&project))?
            } else {
                femcore::linalg::driver::eigs_sym_indefinite(&kr, &mr, sigma, n_ask, &opts, None)?
            };
            let mut unscaled = Vec::with_capacity(pairs.vectors.len());
            for (i, v) in pairs.vectors.iter().enumerate() {
                let yc: Vec<C> = v.iter().map(|&x| c(x, 0.0)).collect();
                let (x, f) = ruiz_unscale_and_rayleigh(&d_ruiz, &yc, kc, mc);
                unscaled.push(x);
                // When lumped_port(s) are present without complex damping, floating conductors bridged only by
                // lumped_ports produce exact 0 Hz DC electrostatic modes (captured in CircuitRom::Q_static); skip them.
                if n_extra_lp > 0 && f.re <= 1e-3 * asm.f0 {
                    continue;
                }
                freqs.push((c(f.re, 0.0), i));
            }
            if n_extra_lp > 0 && freqs.len() > n_modes {
                freqs.sort_by(|a, b| (a.0.re - asm.f0).abs().total_cmp(&(b.0.re - asm.f0).abs()));
                freqs.truncate(n_modes);
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
                    n_ask,
                    1e-12,
                    500,
                    Some(&project),
                )?
            } else {
                femcore::linalg::driver::eigs_complex_projected(&kc_eq, &mc_eq, c(sigma, 0.0), n_ask, 1e-12, 500, None)?
            };
            let mut unscaled = Vec::with_capacity(pairs.vectors.len());
            for (i, v) in pairs.vectors.iter().enumerate() {
                let (x, f) = if asm.c_damp_red.nnz() > 0 {
                    ruiz_unscale_and_rayleigh_damped(&d_ruiz, v, &asm.kc_no_damp, &asm.c_damp_red, mc)
                } else {
                    ruiz_unscale_and_rayleigh(&d_ruiz, v, kc, mc)
                };
                unscaled.push(x);
                if pairs.residuals.get(i).is_some_and(|r| *r > 1e-5) {
                    continue;
                }
                if n_extra_lp > 0 && f.re <= 1e-3 * asm.f0 {
                    continue;
                }
                freqs.push((f, i));
            }
            if n_extra_lp > 0 && freqs.len() > n_modes {
                freqs.sort_by(|a, b| (a.0.re - asm.f0).abs().total_cmp(&(b.0.re - asm.f0).abs()));
                freqs.truncate(n_modes);
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
    let pert_eval = if !asm.pert_ibc_tris.is_empty() {
        Some(PerturbativeIbcEvaluator::new(mesh, ned, &geo, &asm.mur_phys, &asm.pert_ibc_tris)?)
    } else {
        None
    };
    let mut out_freqs: Vec<(C, usize)> = Vec::with_capacity(freqs.len());
    let (axes, left_handed) = femcore::section::axes_of_3d_periodic(&problem.periodic, mesh);
    for (j, &(f_raw, i)) in freqs.iter().enumerate() {
        let x = red.expand(&vectors[i]);
        let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let mut curl_e = vec![[c(0.0, 0.0); 3]; nc * nloc];
        let (mut w, mut cu) = (Vec::new(), Vec::new());
        let mut gs = [[0.0; 3]; 4];
        for ci in 0..nc {
            let gid = ned.key_pairs(mesh, ci);
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
                    Some(ct) => {
                        let gq = safe_curved_jacobian(ct, &lam, &gs);
                        basis3_keyed(order, &lam, &gq, &gid, &mut w, &mut cu);
                    }
                    None => basis3_keyed(order, &lam, &gs, &gid, &mut w, &mut cu),
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
        let (s_e, mut energy_row, mut equip_err) = compute_mode_energy_and_scale(f_raw, &vectors[i], rot, asm);
        let mut f = f_raw;
        // FR-30: Perturbative Wheeler surface-current kinetic inductance and surface-resistance shifts
        if let Some(ref pe) = pert_eval {
            let scale_phys = rot * s_e;
            let x_phys: Vec<C> = x.iter().map(|&z| z * scale_phys).collect();
            let (sum_lk_ktilde2, sum_rs_ktilde2) = pe.eval_mode(&x_phys);
            let omega_0 = (2.0 * PI * f_raw.re).max(1e-300);
            let k0_0 = omega_0 / C0;
            let du_kin_0 = (EPS0 / (4.0 * k0_0 * k0_0)) * sum_lk_ktilde2;
            let u_mag_0 = (energy_row.u_mag_vol_j + energy_row.u_kin_ibc_j + energy_row.u_ind_jj_j.iter().sum::<f64>())
                .max(1e-300);
            let alpha_kin = (du_kin_0 / u_mag_0).max(0.0);
            let inv_scale = 1.0 / (1.0 + alpha_kin);
            let f_re_new = f_raw.re * inv_scale.sqrt();
            energy_row.u_mag_vol_j *= inv_scale;
            energy_row.u_kin_ibc_j = (energy_row.u_kin_ibc_j + du_kin_0) * inv_scale;
            for u in &mut energy_row.u_ind_jj_j {
                *u *= inv_scale;
            }
            let u_mag_tot = energy_row.u_mag_vol_j + energy_row.u_kin_ibc_j + energy_row.u_ind_jj_j.iter().sum::<f64>();
            equip_err = (energy_row.u_elec_j - u_mag_tot).abs() / energy_row.u_elec_j.max(1e-300);
            let p_cond = (0.5 * inv_scale / (omega_0 * omega_0 * MU0 * MU0)) * sum_rs_ktilde2;
            let u_tot = (energy_row.u_elec_j + u_mag_tot).max(1e-300);
            let df_im = p_cond / (4.0 * PI * u_tot);
            f = c(f_re_new, f_raw.im + df_im);
            energy_row.freq_hz = f.re;
            energy_row.q_factor = (f.im.abs() > 1e-12 * f.re.abs().max(1.0)).then(|| f.re / (2.0 * f.im.abs()));
        }
        out_freqs.push((f, i));
        max_equip_err = max_equip_err.max(equip_err);

        let mut cells = CellFields::new(mesh, 2);
        let (e_scale, e_unit) = if asm.has_cqed { (rot * s_e, "V/m") } else { (rot * s, "a.u.") };
        for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
            let z: Vec<C> = e.iter().map(|v| v[k] * e_scale).collect();
            cells.push_complex(*name, e_unit, z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
        }

        // Magnetic field H_m = (i / (omega_m mu_0 mu_r)) curl E_m scaled consistently with E_m
        // For DC / zero-frequency null modes (e.g. k=0 uniform fields), avoid 1/(1e-300 * MU0)^2 overflow.
        let omega_raw = 2.0 * PI * f.re;
        let omega_m = omega_raw.max(1e-300);
        let h_prefactor =
            if omega_raw * MU0 >= 1e-100 { c(0.0, 1.0 / (omega_raw * MU0)) * e_scale } else { c(0.0, 0.0) };
        let mut h_field = vec![[c(0.0, 0.0); 3]; nc * nloc];
        for ci in 0..nc {
            let mui = c(1.0 / asm.mur_phys[ci], 0.0);
            for a in 0..nloc {
                let idx = ci * nloc + a;
                let ce = curl_e[idx];
                h_field[idx] = [ce[0] * mui * h_prefactor, ce[1] * mui * h_prefactor, ce[2] * mui * h_prefactor];
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
                            let ep = [e[idx][0] * e_scale, e[idx][1] * e_scale, e[idx][2] * e_scale];
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

        let m_metrics = compute_3d_mode_metrics(
            mesh,
            &asm.ned,
            &geo,
            &x,
            c(0.0, 0.0),
            omega_m,
            &asm.pml_cells,
            &asm.eps_phys_re,
            &asm.d_omega_eps,
            &asm.mur_phys,
            asm.y_periodic_ky.is_some(),
        );

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
        scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(f.re)));
        scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(f.im)));
        scalars.push((format!("layer fraction {}", j + 1), Scalar::Num(m_metrics.layer_fraction)));
        scalars.push((format!("Ex share {}", j + 1), Scalar::Num(m_metrics.comp_shares[0])));
        scalars.push((format!("Ey share {}", j + 1), Scalar::Num(m_metrics.comp_shares[1])));
        scalars.push((format!("Ez share {}", j + 1), Scalar::Num(m_metrics.comp_shares[2])));
        for &(label, share) in &m_metrics.domain_shares {
            scalars.push((format!("share in domain {label} {}", j + 1), Scalar::Num(share)));
        }
        if let Some(ky) = asm.y_periodic_ky {
            let k0_m = (2.0 * PI * f.re / C0).max(1e-300);
            scalars.push((format!("n_B {}", j + 1), Scalar::Num(ky / k0_m)));
            scalars.push((format!("v_g {} (m/s)", j + 1), Scalar::Num(m_metrics.vg)));
            if m_metrics.vg.abs() > 1e-9 * C0 {
                scalars.push((format!("n_g {}", j + 1), Scalar::Num(C0 / m_metrics.vg)));
            }
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
    if !asm.pert_ibc_tris.is_empty() {
        log.push(format!("surface_impedance perturbative mode active on {} triangle(s)", asm.pert_ibc_tris.len()));
    }
    log.extend(asm.ibc_warnings.iter().cloned());
    // the frame the polarization is read in: the (lateral, normal, propagation) axes of the problem
    scalars.push(
        nyx_physics::FieldFrame {
            lateral: axes[0],
            normal: axes[1],
            propagation: axes[2],
            left_handed,
            outside: crate::optical::absorbing_domains(problem),
        }
        .scalar(),
    );
    if asm.has_cqed {
        let rom = build_circuit_rom(mesh, problem, asm, vectors, &out_freqs, mode_energy_rows, pert_eval.as_ref())?;
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

/// 3D mode energy shares, PML layer fraction, and (when `compute_vg` is true) single-solve Poynting `v_g`
/// evaluated over non-PML physical cells with unstretched `eps_r`, `d(omega eps_r)/domega`, and `mu_r`.
/// `ky_envelope` is `0` for `cavity` (where `x` already carries the Floquet phase) and `n_B * k_0` for `bloch`
/// (where `x` is the `y`-periodic envelope `e_tilde`).
pub(crate) struct Mode3dMetrics {
    pub vg: f64,
    pub comp_shares: [f64; 3],
    pub domain_shares: Vec<(i32, f64)>,
    pub layer_fraction: f64,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_3d_mode_metrics(
    mesh: &Discretisation,
    ned: &Ned3,
    geo: &Geometry,
    x: &[C],
    ky_envelope: C,
    omega: f64,
    pml_cells: &[bool],
    eps_phys_re: &[[[f64; 3]; 3]],
    d_omega_eps: &[[[f64; 3]; 3]],
    mur_phys: &[f64],
    compute_vg: bool,
) -> Mode3dMetrics {
    let any_pml = pml_cells.iter().any(|&p| p);
    let mut domains: Vec<i32> =
        (0..mesh.n_cells()).filter(|&ci| !pml_cells[ci]).map(|ci| mesh.cell_labels[ci]).collect();
    domains.sort_unstable();
    domains.dedup();
    let compute_vg = compute_vg && (omega * MU0) >= 1e-100;
    let quad = quad_rule(3, 4);
    let cquad = quad_rule_collapsed(if mesh.has_qp() { 8 } else { CURVED_QUAD_DEGREE });
    let nb = ned.nb();
    let j_ky = c(0.0, 1.0) * ky_envelope;

    let mut p_y = 0.0_f64;
    let mut u_em = 0.0_f64;
    let mut comp_e = [0.0_f64; 3];
    let mut dom_e = vec![0.0_f64; domains.len()];
    let mut u_e_phys = 0.0_f64;
    let mut u_e_pml = 0.0_f64;

    let (mut w, mut cu) = (Vec::new(), Vec::new());
    let mut gs = [[0.0; 3]; 4];

    for ci in 0..mesh.n_cells() {
        let gid = ned.key_pairs(mesh, ci);
        for b in 0..4 {
            for k in 0..3 {
                gs[b][k] = geo.g[ci * 12 + b * 3 + k];
            }
        }
        let is_pml = pml_cells[ci];
        let d_idx =
            if is_pml { usize::MAX } else { domains.binary_search(&mesh.cell_labels[ci]).unwrap_or(usize::MAX) };
        let er = &eps_phys_re[ci];
        let de = &d_omega_eps[ci];
        let mur = mur_phys[ci];
        let inv_om_mu = if compute_vg { c(0.0, 1.0 / (omega * MU0 * mur)) } else { c(0.0, 0.0) };

        let mut eval_point = |w_s: &[[f64; 3]], cu_s: &[[f64; 3]], wq: f64| {
            let mut eq = [c(0.0, 0.0); 3];
            let mut cq = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let cf = x[ned.cell_dofs[ci][s]];
                for k in 0..3 {
                    eq[k] += cf * w_s[s][k];
                    cq[k] += cf * cu_s[s][k];
                }
            }
            let mut e_eps_e = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    e_eps_e += er[a][b] * (eq[a].conj() * eq[b]).re;
                }
            }
            if is_pml {
                u_e_pml += e_eps_e * wq;
            } else {
                u_e_phys += e_eps_e * wq;
                for k in 0..3 {
                    comp_e[k] += er[k][k].abs() * eq[k].norm_sqr() * wq;
                }
                if d_idx != usize::MAX {
                    dom_e[d_idx] += e_eps_e * wq;
                }
                if compute_vg {
                    let curl_k = [cq[0] - j_ky * eq[2], cq[1], cq[2] + j_ky * eq[0]];
                    let hq = [curl_k[0] * inv_om_mu, curl_k[1] * inv_om_mu, curl_k[2] * inv_om_mu];
                    p_y += 0.5 * (eq[2].conj() * hq[0] - eq[0].conj() * hq[2]).re * wq;
                    let mut e_de_e = 0.0_f64;
                    for a in 0..3 {
                        for b in 0..3 {
                            e_de_e += de[a][b] * (eq[a].conj() * eq[b]).re;
                        }
                    }
                    let h_mu_h = mur * (hq[0].norm_sqr() + hq[1].norm_sqr() + hq[2].norm_sqr());
                    u_em += 0.25 * (EPS0 * e_de_e + MU0 * h_mu_h) * wq;
                }
            }
        };

        if mesh.cell_is_curved(ci) {
            let ct = CurvedTet::new(mesh.cell_geometry_nodes(ci));
            for q in 0..cquad.w.len() {
                let lam = [cquad.lam[q][0], cquad.lam[q][1], cquad.lam[q][2], cquad.lam[q][3]];
                let (gq, det) = ct.jacobian(&lam);
                basis3_keyed(ned.order, &lam, &gq, &gid, &mut w, &mut cu);
                eval_point(&w, &cu, cquad.w[q] * det.abs() / 6.0);
            }
        } else {
            for q in 0..quad.w.len() {
                let lam = [quad.lam[q][0], quad.lam[q][1], quad.lam[q][2], quad.lam[q][3]];
                basis3_keyed(ned.order, &lam, &gs, &gid, &mut w, &mut cu);
                eval_point(&w, &cu, quad.w[q] * geo.vol[ci]);
            }
        }
    }

    let vg_raw = if compute_vg && u_em > 0.0 { p_y / u_em } else { 0.0 };
    let vg = if vg_raw.abs() <= 1e-9 * C0 { 0.0 } else { vg_raw };
    let comp_tot = comp_e[0] + comp_e[1] + comp_e[2];
    let comp_shares =
        if comp_tot > 0.0 { [comp_e[0] / comp_tot, comp_e[1] / comp_tot, comp_e[2] / comp_tot] } else { [0.0; 3] };
    let dom_tot: f64 = dom_e.iter().sum();
    let domain_shares = domains
        .iter()
        .enumerate()
        .map(|(d, &label)| (label, if dom_tot > 0.0 { dom_e[d] / dom_tot } else { 0.0 }))
        .collect();
    let u_all = u_e_phys + u_e_pml;
    let layer_fraction = if any_pml && u_all > 0.0 { (u_e_pml / u_all).clamp(0.0, 1.0) } else { 0.0 };
    Mode3dMetrics { vg, comp_shares, domain_shares, layer_fraction }
}

// ------------------------------------------------------------------ Native 3D MixedDiscretisation (`Tet + Prism + Hex`) cavity solver

/// Masks of constrained DOFs, edges, vertices, and quad faces on a 3D `MixedDiscretisation` for `Ams` / `AmsP2`.
pub struct MixedCavityMasks {
    pub fixed_dofs: Vec<bool>,
    pub fixed_edges: Vec<bool>,
    pub fixed_verts: Vec<bool>,
    pub fixed_quad_faces: Vec<bool>,
}

/// Solves a 3D `Cavity` study natively on `Tet + Prism + Hex` (`MixedDiscretisation`), projecting Stage-3 display fields onto `SimplexMesh`.
pub fn solve_cavity_mixed(mesh: &Mesh, mixed: &MixedDiscretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let (asm, ned_mixed, _) = assemble_emw_eigen_mixed(mixed, problem, false)?;
    let solved = solve_emw_modes(&asm)?;
    let mut r = emw_eigen_fields_mixed(mesh, mixed, &ned_mixed, problem, &asm, &solved)?;
    r.seconds = t0.elapsed().as_secs_f64();
    Ok(r)
}

/// Solves a 3D `Cavity` study natively on `Tet + Prism + Hex` via iterative Jacobi-Davidson preconditioned by `Ams` (`P1`) or `AmsP2` (`P2`),
/// enforcing the `Rational -> KappaPolynomial { kappa_max: 2.0 }` PML conditioning safeguard.
pub fn solve_cavity_mixed_ams(
    mesh: &Mesh,
    mixed: &MixedDiscretisation,
    problem: &Problem,
) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let (asm, ned_mixed, masks) = assemble_emw_eigen_mixed(mixed, problem, true)?;
    if !problem.periodic.is_empty() {
        return Err(
            "solve_cavity_mixed_ams requires non-periodic boundaries (use solve_cavity_mixed for Floquet cells)".into(),
        );
    }
    let pairs = eigs_mixed_curl_curl_ams(
        mixed,
        &ned_mixed,
        &masks.fixed_dofs,
        &masks.fixed_edges,
        &masks.fixed_verts,
        &masks.fixed_quad_faces,
        &asm.kc,
        &asm.mc,
        c(asm.sigma, 0.0),
        asm.n_modes,
        1e-9,
        120,
    )?;
    let mut freqs = Vec::with_capacity(pairs.values.len());
    for (i, &l) in pairs.values.iter().enumerate() {
        if pairs.residuals.get(i).is_some_and(|r| *r > 1e-4) {
            continue;
        }
        let mut f = l.sqrt() * (C0 / (2.0 * PI));
        if f.re < 0.0 {
            f = -f;
        }
        if f.re <= 1e-3 * asm.f0 {
            continue;
        }
        freqs.push((f, i));
    }
    freqs.sort_by(|a, b| a.0.re.partial_cmp(&b.0.re).unwrap());
    let solved = EmwSolve { freqs, vectors: pairs.vectors, dim: pairs.dim, converged: pairs.converged };
    let mut r = emw_eigen_fields_mixed(mesh, mixed, &ned_mixed, problem, &asm, &solved)?;
    r.seconds = t0.elapsed().as_secs_f64();
    Ok(r)
}

/// Stage 1 on a 3D `MixedDiscretisation` (`Tet + Prism + Hex`).
pub fn assemble_emw_eigen_mixed(
    mesh: &MixedDiscretisation,
    problem: &Problem,
    is_iterative: bool,
) -> Result<(EmwAssembly, NedMixed, MixedCavityMasks), String> {
    if mesh.dim != 3 {
        return Err("the cavity study of EM waves needs a 3D mesh".into());
    }
    let Study::Cavity { n_modes, f_search } = problem.study.clone() else {
        return Err("EM waves on a 3D mesh use the cavity study".into());
    };
    validate_port_contacts_3d_mixed(mesh, problem)?;
    let order = if problem.order == 1 { 1 } else { 2 };
    let mats = MixedCellMaterials::from_problem(mesh, problem)?;
    let (lo, hi) = mesh.bbox();
    let size = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(1e-300);
    let f0 = if f_search > 0.0 { f_search } else { C0 / (2.0 * size) };
    let sigma = (2.0 * PI * f0 / C0).powi(2);

    let pml_table = EmPmlTable::from_mixed(mesh, &problem.materials, f0, true, is_iterative);
    let pml_cells = pml_table.pml_cells_mask();

    let mut tl = Instant::now();
    let mut matches = Vec::with_capacity(problem.periodic.len());
    let mut y_periodic_ky = None;
    for p in &problem.periodic {
        let mt = match_vertices_mixed(mesh, p)?;
        let sn = (mt.shift[0] * mt.shift[0] + mt.shift[1] * mt.shift[1] + mt.shift[2] * mt.shift[2]).sqrt();
        if sn > 0.0 && mt.shift[1].abs() > 0.99 * sn && y_periodic_ky.is_none() {
            y_periodic_ky = Some(p.phase / mt.shift[1]);
        }
        matches.push(mt);
    }
    let vkey = vertex_keys_mixed(mesh, &matches);
    let ned = NedMixed::with_vkey(mesh, order, vkey);
    let n = ned.n_dofs;
    let (k, m) =
        assemble_maxwell_mixed(mesh, &ned, &mats.eps_0, &mats.mu_inv_0, pml_table.has_any().then_some(&pml_table));
    lap(&mut tl, "emw3d(mixed): numbering and assembly");

    let has_cqed =
        !problem.surface_impedance.is_empty() || !problem.lumped_jj.is_empty() || !problem.lumped_ports.is_empty();
    let has_sheets = has_cqed || !problem.pec.is_empty();
    let mut ibc_jj = if has_sheets {
        Some(assemble_ibc_and_jj_3d_mixed_cavity(
            mesh,
            &ned,
            &problem.pec,
            &problem.surface_impedance,
            &problem.lumped_jj,
            f0,
        )?)
    } else {
        None
    };
    let (mut lumped_ports, lp_facets) =
        assemble_lumped_ports_3d_mixed(mesh, &ned, &problem.pec, &problem.lumped_ports)?;

    let periodic_labels: Vec<i32> =
        problem.periodic.iter().flat_map(|p| p.source.iter().chain(p.dest.iter()).copied()).collect();
    let direct_ibc_labels: Vec<i32> =
        ibc_jj.as_ref().map_or_else(Vec::new, |(ij, _, _, _)| ij.direct_ibc_labels.clone());
    let jj_labels: Vec<i32> = problem.lumped_jj.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let lp_labels: Vec<i32> = problem.lumped_ports.iter().flat_map(|b| b.labels.iter().copied()).collect();
    let (sym_pec_facets, sym_pmc_facets) = symmetry_facet_flags_mixed(mesh, problem.symmetry)?;

    let mut fixed_set = std::collections::BTreeSet::new();
    let mut ext_pec_vertex = std::collections::BTreeSet::new();
    let mut pec_edge = std::collections::BTreeSet::new();
    let mut pec_quad_face = std::collections::BTreeSet::new();

    let mark_facet_topology =
        |f: usize, pec_e: &mut std::collections::BTreeSet<usize>, pec_q: &mut std::collections::BTreeSet<usize>| {
            let fv = mesh.facet(f);
            let le: &[(usize, usize)] = match mesh.facet_type(f) {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in le {
                pec_e.insert(mesh.edge_index(fv[i], fv[j]).unwrap());
            }
            if mesh.facet_type(f) == FacetType::Quad4 {
                pec_q.insert(mesh.quad_face_index([fv[0], fv[1], fv[2], fv[3]]).unwrap());
            }
        };

    for f in 0..mesh.n_facets() {
        let label = mesh.facet_labels[f];
        let is_ext_pec = mesh.exterior_facets[f]
            && !periodic_labels.contains(&label)
            && (sym_pec_facets[f]
                || (!problem.pmc.contains(&label)
                    && !sym_pmc_facets[f]
                    && !direct_ibc_labels.contains(&label)
                    && !jj_labels.contains(&label)
                    && !lp_labels.contains(&label)));
        if is_ext_pec {
            for d in ned.facet_dofs(mesh, f) {
                fixed_set.insert(d);
            }
            for &v in mesh.facet(f) {
                ext_pec_vertex.insert(v as usize);
            }
            mark_facet_topology(f, &mut pec_edge, &mut pec_quad_face);
        }
    }

    let mut cond_dsu_edges: Vec<(usize, usize)> = Vec::new();
    let mut cond_vertices = std::collections::BTreeSet::new();
    let need_gt_r = !lumped_ports.is_empty()
        || ibc_jj
            .as_ref()
            .is_some_and(|(ij, _, _, _)| ij.jj_ports.iter().any(|p| p.bc.formulation != "distributed_sheet"));
    let gt_r_opt = if need_gt_r { Some(gradient_mixed_full(mesh, &ned).transpose()) } else { None };
    let nv = mesh.n_points();
    let ne = mesh.n_edges();

    if let Some((ref ij, ref pec_fcs, ref ibc_fcs, ref jj_fcs)) = ibc_jj {
        for &f in pec_fcs {
            for d in ned.facet_dofs(mesh, f) {
                fixed_set.insert(d);
            }
            let fv = mesh.facet(f);
            for &v in fv {
                cond_vertices.insert(v as usize);
            }
            let le: &[(usize, usize)] = match mesh.facet_type(f) {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in le {
                cond_dsu_edges.push((fv[i] as usize, fv[j] as usize));
            }
            mark_facet_topology(f, &mut pec_edge, &mut pec_quad_face);
        }
        for &f in ibc_fcs {
            let fv = mesh.facet(f);
            for &v in fv {
                cond_vertices.insert(v as usize);
            }
            let le: &[(usize, usize)] = match mesh.facet_type(f) {
                FacetType::Edge2 => &[(0, 1)],
                FacetType::Tri3 => &TRI_EDGES,
                FacetType::Quad4 => &QUAD_EDGES,
            };
            for &(i, j) in le {
                cond_dsu_edges.push((fv[i] as usize, fv[j] as usize));
            }
            mark_facet_topology(f, &mut pec_edge, &mut pec_quad_face);
        }
        for (port, p_fcs) in ij.jj_ports.iter().zip(jj_fcs) {
            if port.bc.formulation == "distributed_sheet" {
                for &f in p_fcs {
                    let fv = mesh.facet(f);
                    for &v in fv {
                        cond_vertices.insert(v as usize);
                    }
                    let le: &[(usize, usize)] = match mesh.facet_type(f) {
                        FacetType::Edge2 => &[(0, 1)],
                        FacetType::Tri3 => &TRI_EDGES,
                        FacetType::Quad4 => &QUAD_EDGES,
                    };
                    for &(i, j) in le {
                        cond_dsu_edges.push((fv[i] as usize, fv[j] as usize));
                    }
                    mark_facet_topology(f, &mut pec_edge, &mut pec_quad_face);
                }
            }
        }
    }

    // FR-26: Bridge only true electrode seam vertices on `rank1` `lumped_jj` and `lumped_port` sheets.
    let pec_ibc_verts = cond_vertices.clone();
    let register_rank1_port_mixed =
        |p_fcs: &[usize],
         b_full: &[f64],
         gt_r: &Csr,
         cond_vertices: &mut std::collections::BTreeSet<usize>,
         cond_dsu_edges: &mut Vec<(usize, usize)>,
         pec_edge: &mut std::collections::BTreeSet<usize>,
         pec_quad_face: &mut std::collections::BTreeSet<usize>| {
            let gt_b = gt_r.matvec(b_full);
            let b_max = gt_b.iter().copied().map(f64::abs).fold(0.0_f64, f64::max);
            let mut port_verts = std::collections::BTreeSet::new();
            for &f in p_fcs {
                for &v in mesh.facet(f) {
                    port_verts.insert(v as usize);
                }
            }
            let mut term_verts: Vec<usize> = Vec::new();
            for v in port_verts {
                let is_term = (b_max > 0.0 && gt_b[v].abs() > 1e-8 * b_max)
                    || ext_pec_vertex.contains(&v)
                    || pec_ibc_verts.contains(&v)
                    || cond_vertices.contains(&v);
                if is_term {
                    term_verts.push(v);
                }
            }
            for &v in &term_verts {
                cond_vertices.insert(v);
            }
            for win in term_verts.windows(2) {
                cond_dsu_edges.push((win[0], win[1]));
            }
            if order == 2 && b_max > 0.0 {
                for e in 0..ne {
                    if gt_b[nv + e].abs() > 1e-8 * b_max {
                        pec_edge.insert(e);
                    }
                }
                for qf in 0..mesh.n_quad_faces() {
                    if gt_b[nv + ne + qf].abs() > 1e-8 * b_max {
                        pec_quad_face.insert(qf);
                    }
                }
            }
        };
    if let Some(ref gt_r) = gt_r_opt {
        if let Some((ref ij, _, _, ref jj_fcs)) = ibc_jj {
            for (port, p_fcs) in ij.jj_ports.iter().zip(jj_fcs) {
                if port.bc.formulation != "distributed_sheet" {
                    register_rank1_port_mixed(
                        p_fcs,
                        &port.b_full,
                        gt_r,
                        &mut cond_vertices,
                        &mut cond_dsu_edges,
                        &mut pec_edge,
                        &mut pec_quad_face,
                    );
                }
            }
        }
        for (lp, p_fcs) in lumped_ports.iter().zip(&lp_facets) {
            register_rank1_port_mixed(
                p_fcs,
                &lp.b_full,
                gt_r,
                &mut cond_vertices,
                &mut cond_dsu_edges,
                &mut pec_edge,
                &mut pec_quad_face,
            );
        }
    }

    let fixed: Vec<usize> = fixed_set.iter().copied().collect();
    let mut relations = Vec::new();
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        relations.extend(ned_mixed_periodic_relations(mesh, &ned, &mt.src_facets, &mt.map, p.phase)?);
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
    let (k_phys_red, m_phys_red, m0_full) = if pml_cells.iter().any(|&p| p) {
        let mut eps_phys_c: Vec<C3> = Vec::with_capacity(mesh.n_cells());
        let mut mu_phys_c: Vec<C3> = Vec::with_capacity(mesh.n_cells());
        for ci in 0..mesh.n_cells() {
            if pml_cells[ci] {
                eps_phys_c.push([[c(0.0, 0.0); 3]; 3]);
                mu_phys_c.push([[c(0.0, 0.0); 3]; 3]);
            } else {
                let mut e = [[c(0.0, 0.0); 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        e[a][b] = c(mats.eps_phys_re[ci][a][b], 0.0);
                    }
                }
                eps_phys_c.push(e);
                mu_phys_c.push(iso_mu_inv_c3(mats.mur_phys[ci]));
            }
        }
        let (kp, mp) = assemble_maxwell_mixed(mesh, &ned, &eps_phys_c, &mu_phys_c, None);
        let kp_r = to_real(&femcore::periodic::congruence(&pmat, &kp));
        let mp_r = to_real(&femcore::periodic::congruence(&pmat, &mp));
        let mf = if has_cqed { to_real(&mp) } else { Csr::zeros(0, 0) };
        (kp_r, mp_r, mf)
    } else {
        (k0_red.clone(), m0_red.clone(), if has_cqed { to_real(&m) } else { Csr::zeros(0, 0) })
    };

    let mut c_port_all_full = Csr::zeros(n, n);
    let mut c_port_active_full = Csr::zeros(n, n);
    for lp in &mut lumped_ports {
        lp.b_red = reduce_real_vec(&red, &lp.b_full);
        if lp.c_port_full.nnz() > 0 {
            c_port_all_full = c_port_all_full.add_scaled(&lp.c_port_full, 1.0);
            if lp.bc.complex_damping {
                c_port_active_full = c_port_active_full.add_scaled(&lp.c_port_full, 1.0);
            }
        }
    }
    let c_port_red = reduce_real_csr(&red, &pmat, &c_port_all_full);
    let c_port_active_red = reduce_real_csr(&red, &pmat, &c_port_active_full);

    let k0_f0 = 2.0 * PI * f0 / C0;
    let (kc, kc_no_damp, mc, k_ibc_red, c_ibc_red, c_damp_red, jj_ports, pert_ibc_tris, ibc_warnings, cell_sheet) =
        if let Some((ref mut ij, _, _, _)) = ibc_jj {
            let mut k_tot_no_damp = k.clone();
            let mut k_tot = k.clone();
            let mut m_tot = m.clone();
            if ij.k_ibc_full.nnz() > 0 {
                k_tot_no_damp = k_tot_no_damp.add_scaled(&to_complex(&ij.k_ibc_full), c(1.0, 0.0));
            }
            if ij.kc_ibc_full.nnz() > 0 {
                k_tot = k_tot.add_scaled(&ij.kc_ibc_full, c(1.0, 0.0));
            }
            if ij.k_jj_full.nnz() > 0 {
                let kjj_c = to_complex(&ij.k_jj_full);
                k_tot_no_damp = k_tot_no_damp.add_scaled(&kjj_c, c(1.0, 0.0));
                k_tot = k_tot.add_scaled(&kjj_c, c(1.0, 0.0));
            }
            if c_port_active_full.nnz() > 0 {
                k_tot = k_tot.add_scaled(&to_complex(&c_port_active_full), c(0.0, k0_f0));
            }
            if ij.m_jj_full.nnz() > 0 {
                m_tot = m_tot.add_scaled(&to_complex(&ij.m_jj_full), c(1.0, 0.0));
            }
            let kc_no_damp =
                if has_cqed { femcore::periodic::congruence(&pmat, &k_tot_no_damp) } else { k0_c_red.clone() };
            let kc = if has_cqed { femcore::periodic::congruence(&pmat, &k_tot) } else { k0_c_red };
            let mc = if has_cqed { femcore::periodic::congruence(&pmat, &m_tot) } else { m0_c_red };
            let k_ibc_red = reduce_real_csr(&red, &pmat, &ij.k_ibc_full);
            let c_ibc_red = reduce_real_csr(&red, &pmat, &ij.c_ibc_full);
            let c_damp_red = if c_port_active_red.nnz() > 0 {
                c_ibc_red.add_scaled(&c_port_active_red, 1.0)
            } else {
                c_ibc_red.clone()
            };
            for port in &mut ij.jj_ports {
                port.b_red = reduce_real_vec(&red, &port.b_full);
                port.k_jj_red = reduce_real_csr(&red, &pmat, &port.k_jj_red);
                port.m_jj_red = reduce_real_csr(&red, &pmat, &port.m_jj_red);
            }
            (
                kc,
                kc_no_damp,
                mc,
                k_ibc_red,
                c_ibc_red,
                c_damp_red,
                std::mem::take(&mut ij.jj_ports),
                std::mem::take(&mut ij.pert_ibc_tris),
                std::mem::take(&mut ij.warnings),
                std::mem::take(&mut ij.cell_sheet),
            )
        } else {
            let kc_no_damp = k0_c_red.clone();
            let kc = if c_port_active_full.nnz() > 0 {
                let k_tot = k.add_scaled(&to_complex(&c_port_active_full), c(0.0, k0_f0));
                femcore::periodic::congruence(&pmat, &k_tot)
            } else {
                k0_c_red
            };
            (
                kc,
                kc_no_damp,
                m0_c_red,
                Csr::zeros(nred, nred),
                Csr::zeros(nred, nred),
                c_port_active_red,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                vec![None; mesh.n_cells()],
            )
        };

    // Exact scalar nullspace G_1 (`order == 1`) or G_2 (`order == 2`, `Order2Tensor`: vertices + edges + quad faces + hex bubbles).
    let n_nodes = scalar_nullspace_count_mixed(mesh, order);
    let mut periodic_node_relations = Vec::new();
    for (p, mt) in problem.periodic.iter().zip(&matches) {
        periodic_node_relations.extend(nodal_relations_mixed_3d(mesh, order, mt, p.phase)?);
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
                parent[ra.max(rb)] = ra.min(rb);
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
        node_fixed.extend(pec_quad_face.iter().map(|&qf| nv + ne + qf));
    }
    if pec_vertex.is_empty()
        && nv > 0
        && problem.periodic.iter().all(|p| (femcore::periodic::phase_factor(p.phase) - c(1.0, 0.0)).norm() < 1e-9)
    {
        let mut comp_dsu: Vec<usize> = (0..nv).collect();
        let find_c = |p: &mut [usize], mut i: usize| -> usize {
            while p[i] != i {
                p[i] = p[p[i]];
                i = p[i];
            }
            i
        };
        for &(a, b) in &mesh.edges {
            let (ra, rb) = (find_c(&mut comp_dsu, a as usize), find_c(&mut comp_dsu, b as usize));
            if ra != rb {
                comp_dsu[ra.max(rb)] = ra.min(rb);
            }
        }
        for &(d, s, _) in &node_rel {
            if d < nv && s < nv {
                let (rd, rs) = (find_c(&mut comp_dsu, d), find_c(&mut comp_dsu, s));
                if rd != rs {
                    comp_dsu[rd.max(rs)] = rd.min(rs);
                }
            }
        }
        let mut pinned_roots = std::collections::BTreeSet::new();
        for v in 0..nv {
            let r = find_c(&mut comp_dsu, v);
            if pinned_roots.insert(r) && !node_fixed.contains(&v) {
                node_fixed.push(v);
            }
        }
        node_fixed.sort_unstable();
        node_fixed.dedup();
    }

    let nred_nodes = femcore::periodic::Reduction::build(n_nodes, &node_rel, &node_fixed);
    let g_full = to_complex(&gradient_mixed_full(mesh, &ned));
    let g_nodes = femcore::periodic::ccsr_mul(&g_full, &nred_nodes.prolong());
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
    if g_h.nrows() > 0 {
        let check_orth = |name: &str, b_red: &[f64]| -> Result<(), String> {
            let bnorm = b_red.iter().map(|&v| v * v).sum::<f64>().sqrt();
            if bnorm <= 1e-300 {
                return Ok(());
            }
            let b_c: Vec<C> = b_red.iter().map(|&v| c(v, 0.0)).collect();
            let gtb = g_h.matvec(&b_c);
            let err = gtb.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt() / bnorm;
            if err > 1e-10 {
                return Err(format!(
                    "cavity null-space orthogonality violation on port '{name}': ||G_red^T b_red||_2 / ||b_red||_2 = {err:.3e} > 1e-10"
                ));
            }
            Ok(())
        };
        for port in &jj_ports {
            if port.bc.formulation != "distributed_sheet" {
                check_orth(&port.bc.name, &port.b_red)?;
            }
        }
        for lp in &lumped_ports {
            check_orth(&lp.bc.name, &lp.b_red)?;
        }
    }
    lap(&mut tl, "emw3d(mixed): reduction and gradient matrix");

    let mut fixed_dofs = vec![false; n];
    for &d in &fixed {
        fixed_dofs[d] = true;
    }
    let mut fixed_edges = vec![false; ne];
    for &e in &pec_edge {
        fixed_edges[e] = true;
    }
    let mut fixed_verts = vec![false; nv];
    for &v in &pec_vertex {
        fixed_verts[v] = true;
    }
    let mut fixed_quad_faces = vec![false; mesh.n_quad_faces()];
    for &qf in &pec_quad_face {
        fixed_quad_faces[qf] = true;
    }
    let masks = MixedCavityMasks { fixed_dofs, fixed_edges, fixed_verts, fixed_quad_faces };

    let dummy_ned = Ned3 {
        order,
        n_edges: 0,
        n_faces: 0,
        faces: Vec::new(),
        face_id: HashMap::new(),
        cell_dofs: Vec::new(),
        vkey: Vec::new(),
    };

    Ok((
        EmwAssembly {
            order,
            n,
            nred,
            ned: dummy_ned,
            red,
            kc,
            kc_no_damp,
            mc,
            g_red,
            g_h,
            f0,
            sigma,
            n_modes,
            f_search,
            has_cqed,
            k0_red,
            k_phys_red,
            k_ibc_red,
            c_ibc_red,
            c_port_red,
            c_damp_red,
            m0_full,
            m0_red,
            m_phys_red,
            jj_ports,
            lumped_ports,
            pert_ibc_tris,
            ibc_warnings,
            cell_sheet,
            periodic_node_relations,
            mu_inv: mats.mu_inv_0,
            pml_cells,
            eps_phys_re: mats.eps_phys_re,
            d_omega_eps: mats.d_omega_eps,
            mur_phys: mats.mur_phys,
            y_periodic_ky,
        },
        ned,
        masks,
    ))
}

/// Stage 3 on a 3D `MixedDiscretisation`: evaluates native `NedMixed` fields and projects them onto `SimplexMesh` via `MixedToSimplexProjector`.
pub fn emw_eigen_fields_mixed(
    mesh: &Mesh,
    mixed: &MixedDiscretisation,
    ned: &NedMixed,
    problem: &Problem,
    asm: &EmwAssembly,
    solved: &EmwSolve,
) -> Result<SolveResult, String> {
    let (simplex, proj) = MixedToSimplexProjector::from_mesh(mesh, mixed).map_err(|e| e.to_string())?;
    let (order, n, nred, f0, f_search) = (asm.order, asm.n, asm.nred, asm.f0, asm.f_search);
    let (freqs, vectors, dim, converged) = (&solved.freqs, &solved.vectors, solved.dim, solved.converged);
    let isos: Vec<IsoparametricCell> = (0..mixed.n_cells()).map(|ci| IsoparametricCell::from_mesh(mixed, ci)).collect();
    let nsc = simplex.n_cells();
    let nloc = 10usize;
    let has_sheets = asm.has_cqed || !problem.pec.is_empty();
    let (axes, left_handed) = axes_of_3d_periodic_mixed(&problem.periodic, mixed);

    let mut solutions = Vec::with_capacity(freqs.len());
    let mut scalars = Vec::new();
    let mut mode_energy_rows = Vec::with_capacity(freqs.len());
    let mut max_equip_err = 0.0_f64;
    let pert_eval = if !asm.pert_ibc_tris.is_empty() {
        Some(PerturbativeIbcEvaluator::new_mixed(mixed, ned, &asm.mur_phys, &asm.pert_ibc_tris)?)
    } else {
        None
    };
    let mut out_freqs: Vec<(C, usize)> = Vec::with_capacity(freqs.len());

    for (j, &(f_raw, i)) in freqs.iter().enumerate() {
        let x = asm.red.expand(&vectors[i]);
        let mut scratch = NedEvalScratch::new();
        let mut w = Vec::new();
        let mut cu = Vec::new();
        let sampled = proj.sample_simplex_cell_nodes(mixed, 2, |mc, xi| {
            eval_cell_basis_3d(mixed, ned, mc, &isos[mc], xi, &mut scratch, &mut w, &mut cu);
            let cdofs = ned.cell_dofs(mc);
            let mut val = [c(0.0, 0.0); 3];
            let mut cval = [c(0.0, 0.0); 3];
            for (s, &gd) in cdofs.iter().enumerate() {
                let coef = x[gd];
                for kk in 0..3 {
                    val[kk] += coef * w[s][kk];
                    cval[kk] += coef * cu[s][kk];
                }
            }
            (val, cval)
        });
        let e: Vec<[C; 3]> = sampled.iter().map(|p| p.0).collect();
        let curl_e: Vec<[C; 3]> = sampled.iter().map(|p| p.1).collect();

        let big =
            e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
        let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
        let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
        let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
        let avg = |vals: Vec<f64>| node_average(&simplex, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, 2);
        let fields = vec![
            Field::scalar("|E|", "a.u.", avg(mag.clone())),
            Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
            Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
            Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
        ];

        let (s_e, mut energy_row, mut equip_err) = compute_mode_energy_and_scale(f_raw, &vectors[i], rot, asm);
        let mut f = f_raw;
        if let Some(ref pe) = pert_eval {
            let scale_phys = rot * s_e;
            let x_phys: Vec<C> = x.iter().map(|&z| z * scale_phys).collect();
            let (sum_lk_ktilde2, sum_rs_ktilde2) = pe.eval_mode(&x_phys);
            let omega_0 = (2.0 * PI * f_raw.re).max(1e-300);
            let k0_0 = omega_0 / C0;
            let du_kin_0 = (EPS0 / (4.0 * k0_0 * k0_0)) * sum_lk_ktilde2;
            let u_mag_0 = (energy_row.u_mag_vol_j + energy_row.u_kin_ibc_j + energy_row.u_ind_jj_j.iter().sum::<f64>())
                .max(1e-300);
            let alpha_kin = (du_kin_0 / u_mag_0).max(0.0);
            let inv_scale = 1.0 / (1.0 + alpha_kin);
            let f_re_new = f_raw.re * inv_scale.sqrt();
            energy_row.u_mag_vol_j *= inv_scale;
            energy_row.u_kin_ibc_j = (energy_row.u_kin_ibc_j + du_kin_0) * inv_scale;
            for u in &mut energy_row.u_ind_jj_j {
                *u *= inv_scale;
            }
            let u_mag_tot = energy_row.u_mag_vol_j + energy_row.u_kin_ibc_j + energy_row.u_ind_jj_j.iter().sum::<f64>();
            equip_err = (energy_row.u_elec_j - u_mag_tot).abs() / energy_row.u_elec_j.max(1e-300);
            let p_cond = (0.5 * inv_scale / (omega_0 * omega_0 * MU0 * MU0)) * sum_rs_ktilde2;
            let u_tot = (energy_row.u_elec_j + u_mag_tot).max(1e-300);
            let df_im = p_cond / (4.0 * PI * u_tot);
            f = c(f_re_new, f_raw.im + df_im);
            energy_row.freq_hz = f.re;
            energy_row.q_factor = (f.im.abs() > 1e-12 * f.re.abs().max(1.0)).then(|| f.re / (2.0 * f.im.abs()));
        }
        out_freqs.push((f, i));
        max_equip_err = max_equip_err.max(equip_err);

        let mut cells = CellFields::new(&simplex, 2);
        let (e_scale, e_unit) = if asm.has_cqed { (rot * s_e, "V/m") } else { (rot * s, "a.u.") };
        for (k, name) in ["Ex", "Ey", "Ez"].iter().enumerate() {
            let z: Vec<C> = e.iter().map(|v| v[k] * e_scale).collect();
            cells.push_complex(*name, e_unit, z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
        }

        let omega_raw = 2.0 * PI * f.re;
        let omega_m = omega_raw.max(1e-300);
        let h_prefactor =
            if omega_raw * MU0 >= 1e-100 { c(0.0, 1.0 / (omega_raw * MU0)) * e_scale } else { c(0.0, 0.0) };
        let mut h_field = vec![[c(0.0, 0.0); 3]; nsc * nloc];
        for sc in 0..nsc {
            let mc = proj.parent_cell(sc);
            let mui = c(1.0 / asm.mur_phys[mc], 0.0);
            for a in 0..nloc {
                let idx = sc * nloc + a;
                let ce = curl_e[idx];
                h_field[idx] = [ce[0] * mui * h_prefactor, ce[1] * mui * h_prefactor, ce[2] * mui * h_prefactor];
            }
        }
        for (k, name) in ["Hx", "Hy", "Hz"].iter().enumerate() {
            let z: Vec<C> = h_field.iter().map(|v| v[k]).collect();
            cells.push_complex(*name, "A/m", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
        }

        if has_sheets {
            let mut k_field = vec![[c(0.0, 0.0); 3]; nsc * nloc];
            for sc in 0..nsc {
                let mc = proj.parent_cell(sc);
                let Some(sheet_kind) = asm.cell_sheet[mc] else { continue };
                match sheet_kind {
                    SheetCellKind::Ibc { inv_l_eff, g_eff, n_hat } => {
                        let y_s = c(g_eff, -inv_l_eff / omega_m);
                        for a in 0..nloc {
                            let idx = sc * nloc + a;
                            let ep = [e[idx][0] * e_scale, e[idx][1] * e_scale, e[idx][2] * e_scale];
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
                            k_field[sc * nloc + a] = [k_amp * d_hat[0], k_amp * d_hat[1], k_amp * d_hat[2]];
                        }
                    }
                    SheetCellKind::Pec { n_hat } => {
                        for a in 0..nloc {
                            let idx = sc * nloc + a;
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
        let m_metrics = compute_3d_mode_metrics_mixed(
            mixed,
            ned,
            &isos,
            &x,
            c(0.0, 0.0),
            omega_m,
            &asm.pml_cells,
            &asm.eps_phys_re,
            &asm.d_omega_eps,
            &asm.mur_phys,
            asm.y_periodic_ky.is_some(),
        );

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
        scalars.push((format!("f re {} (Hz)", j + 1), Scalar::Num(f.re)));
        scalars.push((format!("f im {} (Hz)", j + 1), Scalar::Num(f.im)));
        scalars.push((format!("layer fraction {}", j + 1), Scalar::Num(m_metrics.layer_fraction)));
        scalars.push((format!("Ex share {}", j + 1), Scalar::Num(m_metrics.comp_shares[0])));
        scalars.push((format!("Ey share {}", j + 1), Scalar::Num(m_metrics.comp_shares[1])));
        scalars.push((format!("Ez share {}", j + 1), Scalar::Num(m_metrics.comp_shares[2])));
        for &(label, share) in &m_metrics.domain_shares {
            scalars.push((format!("share in domain {label} {}", j + 1), Scalar::Num(share)));
        }
        if let Some(ky) = asm.y_periodic_ky {
            let k0_m = (2.0 * PI * f.re / C0).max(1e-300);
            scalars.push((format!("n_B {}", j + 1), Scalar::Num(ky / k0_m)));
            scalars.push((format!("v_g {} (m/s)", j + 1), Scalar::Num(m_metrics.vg)));
            if m_metrics.vg.abs() > 1e-9 * C0 {
                scalars.push((format!("n_g {}", j + 1), Scalar::Num(C0 / m_metrics.vg)));
            }
        }
    }

    if solutions.is_empty() {
        return Err("no eigenfrequencies converged near the shift".into());
    }
    let mut log = vec![format!(
        "Nedelec mixed order {order}, {n} dofs ({nred} after constraints), shift {:.4e} Hz, Krylov dimension {dim}{}",
        f0,
        if converged { "" } else { " (not every pair converged; those are omitted)" }
    )];
    if !problem.periodic.is_empty() {
        log.push(format!("{} Floquet-periodic pair(s)", problem.periodic.len()));
    }
    if f_search <= 0.0 {
        log.push("no search frequency given: modes nearest the fundamental of a box of this size".into());
    }
    if !asm.pert_ibc_tris.is_empty() {
        log.push(format!("surface_impedance perturbative mode active on {} triangle(s)", asm.pert_ibc_tris.len()));
    }
    log.extend(asm.ibc_warnings.iter().cloned());
    scalars.push(
        nyx_physics::FieldFrame {
            lateral: axes[0],
            normal: axes[1],
            propagation: axes[2],
            left_handed,
            outside: crate::optical::absorbing_domains(problem),
        }
        .scalar(),
    );
    if asm.has_cqed {
        let rom = build_circuit_rom_mixed(
            mixed,
            ned,
            problem,
            asm,
            vectors,
            &out_freqs,
            mode_energy_rows,
            pert_eval.as_ref(),
        )?;
        scalars.push(rom.scalar());
        scalars.push(("equipartition_error".into(), Scalar::Num(max_equip_err)));
    }
    if has_sheets {
        let tls = build_tls_meta(&simplex, problem);
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

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_3d_mode_metrics_mixed(
    mesh: &MixedDiscretisation,
    ned: &NedMixed,
    isos: &[IsoparametricCell],
    x: &[C],
    ky_envelope: C,
    omega: f64,
    pml_cells: &[bool],
    eps_phys_re: &[[[f64; 3]; 3]],
    d_omega_eps: &[[[f64; 3]; 3]],
    mur_phys: &[f64],
    compute_vg: bool,
) -> Mode3dMetrics {
    let any_pml = pml_cells.iter().any(|&p| p);
    let mut domains: Vec<i32> =
        (0..mesh.n_cells()).filter(|&ci| !pml_cells[ci]).map(|ci| mesh.cell_labels[ci]).collect();
    domains.sort_unstable();
    domains.dedup();
    let compute_vg = compute_vg && (omega * MU0) >= 1e-100;
    let j_ky = c(0.0, 1.0) * ky_envelope;

    let mut p_y = 0.0_f64;
    let mut u_em = 0.0_f64;
    let mut comp_e = [0.0_f64; 3];
    let mut dom_e = vec![0.0_f64; domains.len()];
    let mut u_e_phys = 0.0_f64;
    let mut u_e_pml = 0.0_f64;

    let mut scratch = NedEvalScratch::new();
    let mut w = Vec::new();
    let mut cu = Vec::new();

    for ci in 0..mesh.n_cells() {
        let iso = &isos[ci];
        let deg = if iso.curved { CURVED_QUAD_DEGREE } else { 4 };
        let rule = quad_rule_mixed(iso.cell_type, deg);
        let nb = nedelec_nloc(iso.cell_type, ned.order);
        let cdofs = ned.cell_dofs(ci);
        let is_pml = pml_cells[ci];
        let d_idx =
            if is_pml { usize::MAX } else { domains.binary_search(&mesh.cell_labels[ci]).unwrap_or(usize::MAX) };
        let er = &eps_phys_re[ci];
        let de = &d_omega_eps[ci];
        let mur = mur_phys[ci];
        let inv_om_mu = if compute_vg { c(0.0, 1.0 / (omega * MU0 * mur)) } else { c(0.0, 0.0) };

        for (xi, &w_ref) in rule.pts.iter().zip(&rule.w) {
            let jac = eval_cell_basis_3d(mesh, ned, ci, iso, xi, &mut scratch, &mut w, &mut cu);
            let wq = w_ref * jac.det.abs();
            let mut eq = [c(0.0, 0.0); 3];
            let mut cq = [c(0.0, 0.0); 3];
            for s in 0..nb {
                let cf = x[cdofs[s]];
                for k in 0..3 {
                    eq[k] += cf * w[s][k];
                    cq[k] += cf * cu[s][k];
                }
            }
            let mut e_eps_e = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    e_eps_e += er[a][b] * (eq[a].conj() * eq[b]).re;
                }
            }
            if is_pml {
                u_e_pml += e_eps_e * wq;
            } else {
                u_e_phys += e_eps_e * wq;
                for k in 0..3 {
                    comp_e[k] += er[k][k].abs() * eq[k].norm_sqr() * wq;
                }
                if d_idx != usize::MAX {
                    dom_e[d_idx] += e_eps_e * wq;
                }
                if compute_vg {
                    let curl_k = [cq[0] - j_ky * eq[2], cq[1], cq[2] + j_ky * eq[0]];
                    let hq = [curl_k[0] * inv_om_mu, curl_k[1] * inv_om_mu, curl_k[2] * inv_om_mu];
                    p_y += 0.5 * (eq[2].conj() * hq[0] - eq[0].conj() * hq[2]).re * wq;
                    let mut e_de_e = 0.0_f64;
                    for a in 0..3 {
                        for b in 0..3 {
                            e_de_e += de[a][b] * (eq[a].conj() * eq[b]).re;
                        }
                    }
                    let h_mu_h = mur * (hq[0].norm_sqr() + hq[1].norm_sqr() + hq[2].norm_sqr());
                    u_em += 0.25 * (EPS0 * e_de_e + MU0 * h_mu_h) * wq;
                }
            }
        }
    }

    let vg_raw = if compute_vg && u_em > 0.0 { p_y / u_em } else { 0.0 };
    let vg = if vg_raw.abs() <= 1e-9 * C0 { 0.0 } else { vg_raw };
    let comp_tot = comp_e[0] + comp_e[1] + comp_e[2];
    let comp_shares =
        if comp_tot > 0.0 { [comp_e[0] / comp_tot, comp_e[1] / comp_tot, comp_e[2] / comp_tot] } else { [0.0; 3] };
    let dom_tot: f64 = dom_e.iter().sum();
    let domain_shares = domains
        .iter()
        .enumerate()
        .map(|(d, &label)| (label, if dom_tot > 0.0 { dom_e[d] / dom_tot } else { 0.0 }))
        .collect();
    let u_all = u_e_phys + u_e_pml;
    let layer_fraction = if any_pml && u_all > 0.0 { (u_e_pml / u_all).clamp(0.0, 1.0) } else { 0.0 };
    Mode3dMetrics { vg, comp_shares, domain_shares, layer_fraction }
}
