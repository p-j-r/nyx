//! Elastic modes of a body of revolution (a ring, a disk, a wheel) on its half-plane cross-section:
//! the displacement is `u(r, z) exp(-i m theta)` with azimuthal order `m`, and the mesh has `x = r`,
//! `y = z`. The three components are stored as (u_r, u_z, u_theta) so the assembly is that of
//! `sms_kz` with the out-of-plane strain terms replaced by their cylindrical forms,
//!
//! ```text
//!   e_rr = du_r/dr   e_zz = du_z/dz   e_tt = (u_r - i m u_t) / r
//!   g_rz = du_r/dz + du_z/dr   g_zt = du_t/dz - i m u_z / r   g_rt = du_t/dr - u_t / r - i m u_r / r
//! ```
//!
//! and the weight `r dr dz` (the common factor 2 pi dropped). A crystal's stiffness depends on the
//! angle around the axis; the matrices used here are the one given per domain, so only materials
//! with the symmetry of the ring (isotropic, or averaged around the axis with
//! `average_about_axis`) are exact. A mesh must not touch the axis `r = 0`.
//!
//! Absorbing layers are the cylindrical coordinate stretch (ADR 67), the one of `nyx-ems` (ADR 46): inside a layer
//! `r -> r~ = int s_r dr`, so with `s_theta = r~ / r` every `1/r` above becomes `1/(s_theta r)` and the weight `r dr dz`
//! becomes `s_theta r s_r s_z dr dz` (`cylindrical_stretches` builds `[s_r, s_z, s_theta]` per cell). The test function carries
//! the conjugate phase `exp(+i m theta)` only, so the stretch is not conjugated in the test-side strain: the matrices are
//! complex symmetric, not Hermitian, inside a layer. Outside a radial layer `s_theta = 1` and nothing changes.

use femcore::fem::*;
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::material::M6;
use femcore::pml::{cell_stretches, radial_coordinate};
use femcore::props::{Pml, PmlDir};
use nyx_physics::Discretisation;

/// The stiffness (solver Voigt order xx, yy, zz, xy, yz, xz in the frame of the cross-section)
/// averaged over rotations of the crystal about the axis `axis` (0, 1 or 2): what a ring sees when the
/// crystal axes stay fixed while the local frame turns around the ring.
pub fn average_about_axis(c6: &M6, axis: usize, samples: usize) -> M6 {
    // full rank-4 tensor
    let vp = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let mut t = [[[[0.0; 3]; 3]; 3]; 3];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            for (p, q) in [(a, b), (b, a)] {
                for (s, u) in [(cc, d), (d, cc)] {
                    t[p][q][s][u] = c6[i][j];
                }
            }
        }
    }
    let mut avg = [[[[0.0; 3]; 3]; 3]; 3];
    for sidx in 0..samples {
        let th = 2.0 * std::f64::consts::PI * sidx as f64 / samples as f64;
        let (sn, cs) = th.sin_cos();
        let mut r = [[0.0; 3]; 3];
        let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
        r[axis][axis] = 1.0;
        r[a1][a1] = cs;
        r[a1][a2] = -sn;
        r[a2][a1] = sn;
        r[a2][a2] = cs;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let mut s = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                for cc in 0..3 {
                                    for d in 0..3 {
                                        s += r[i][a] * r[j][b] * r[k][cc] * r[l][d] * t[a][b][cc][d];
                                    }
                                }
                            }
                        }
                        avg[i][j][k][l] += s / samples as f64;
                    }
                }
            }
        }
    }
    let mut out = [[0.0; 6]; 6];
    for (i, &(a, b)) in vp.iter().enumerate() {
        for (j, &(cc, d)) in vp.iter().enumerate() {
            out[i][j] = avg[a][b][cc][d];
        }
    }
    out
}

/// Per-cell stretch `[s_r, s_z, s_theta]` of the absorbing layers of an axisymmetric model (mesh x = r, mesh y = z). `s_r` and `s_z` are
/// those of `femcore::pml::cell_stretches`; `s_theta = r~ / r` with the complex radius `r~` of `femcore::pml::radial_coordinate` at the cell
/// centre, as `AxiMaterial::with_pml` does in `nyx-ems` (ADR 46). A cell outside a layer with a radial (`x`) stretch has `s_theta = 1` exactly,
/// so a model without such a layer assembles to the same matrices as before (ADR 67).
pub fn cylindrical_stretches(mesh: &Discretisation, pml: &std::collections::BTreeMap<i32, Pml>) -> Vec<[C; 3]> {
    let mut out = cell_stretches(mesh, pml);
    // radial extent of each layer domain
    let mut boxes: std::collections::BTreeMap<i32, (f64, f64)> = Default::default();
    for ci in 0..mesh.n_cells() {
        let l = mesh.cell_labels[ci];
        if pml.get(&l).is_some_and(|p| p.x != PmlDir::Off) {
            let e = boxes.entry(l).or_insert((f64::INFINITY, f64::NEG_INFINITY));
            for &v in mesh.cell(ci) {
                let r = mesh.point(v as usize)[0];
                e.0 = e.0.min(r);
                e.1 = e.1.max(r);
            }
        }
    }
    for (ci, s) in out.iter_mut().enumerate() {
        let l = mesh.cell_labels[ci];
        let (Some(p), Some(&(lo, hi))) = (pml.get(&l), boxes.get(&l)) else { continue };
        let r = (0..3).map(|k| mesh.point(mesh.cell(ci)[k] as usize)[0]).sum::<f64>() / 3.0;
        s[2] = radial_coordinate(p, r, lo, hi) / r;
    }
    out
}

/// Stiffness matrix over `3 * n_dofs(order)` unknowns (component blocks u_r, u_z, u_theta). `stretch` holds `[s_r, s_z, s_theta]` per cell
/// (`cylindrical_stretches`).
// The B-matrix tables below are written `row * nd + component * nloc + a` for every entry, including the ones with a
// zero factor, so that all lines read alike.
#[allow(clippy::erasing_op)]
pub fn assemble_stiffness_axisym(
    mesh: &Discretisation,
    c6: &[M6],
    order: usize,
    m: f64,
    stretch: Option<&[[C; 3]]>,
) -> Result<CCsr, String> {
    let re = RefElement::new(2, order, 2 * order + 2);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let nd = 3 * nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut grad = vec![0.0; nloc * 2];
    let mut b0 = vec![c(0.0, 0.0); 6 * nd];
    let mut b1 = vec![0.0; 6 * nd];
    let mut b2 = vec![0.0; 6 * nd];
    let mut db = vec![c(0.0, 0.0); 6 * nd];
    let mut ke = vec![c(0.0, 0.0); nd * nd];
    // shape-function values at the quadrature points of P1 give r
    for ci in 0..mesh.n_cells() {
        ke.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let dm = &c6[ci];
        let (sx, sy, st) = stretch.map_or((c(1.0, 0.0), c(1.0, 0.0), c(1.0, 0.0)), |s| (s[ci][0], s[ci][1], s[ci][2]));
        // s_theta = r~ / r: the cylindrical part of the stretch. The `1/r` of the strains becomes `1/(s_theta r)` and the
        // weight `r` becomes `s_theta r`; with s_theta = 1 (every cell outside a radial layer) nothing changes.
        let ist = c(1.0, 0.0) / st;
        let det = sx * sy * st;
        let v = mesh.cell(ci);
        let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
        for q in 0..re.w.len() {
            shape_gradients(&re, &geo, ci, q, &mut grad);
            // r at this point from the vertex shape functions: for P1 the first three are the barycentric coordinates
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            if r <= 0.0 {
                return Err(
                    "the mesh of an axisymmetric problem must lie in r > 0 (it touches or crosses the axis)".into()
                );
            }
            // beta multiplies the trial (u) side, beta_h the test (v) side: the test function carries exp(+i m theta), and
            // only that phase is conjugated, not the complex stretch, so beta_h = -beta (not conj(beta)) inside a layer
            let beta = c(0.0, -m / r) * ist;
            let beta_h = c(0.0, m / r) * ist;
            let inv_r = ist / r;
            b0.iter_mut().for_each(|x| *x = c(0.0, 0.0));
            b1.iter_mut().for_each(|x| *x = 0.0);
            b2.iter_mut().for_each(|x| *x = 0.0);
            for a in 0..nloc {
                let (gx, gy, nv) = (c(grad[a * 2], 0.0) / sx, c(grad[a * 2 + 1], 0.0) / sy, re.n[q][a]);
                // rows: 0 rr, 1 zz, 2 tt, 3 rz, 4 zt, 5 rt ; components 0 u_r, 1 u_z, 2 u_t
                b0[0 * nd + 0 * nloc + a] = gx;
                b0[nd + nloc + a] = gy;
                b2[2 * nd + 0 * nloc + a] = nv;
                b1[2 * nd + 2 * nloc + a] = nv;
                b0[3 * nd + 0 * nloc + a] = gy;
                b0[3 * nd + nloc + a] = gx;
                b0[4 * nd + 2 * nloc + a] = gy;
                b1[4 * nd + nloc + a] = nv;
                b0[5 * nd + 2 * nloc + a] = gx;
                b2[5 * nd + 2 * nloc + a] = -nv;
                b1[5 * nd + 0 * nloc + a] = nv;
            }
            // the tt row has the term u_t * (-i m / r) too: it sits in b1 on the third component
            for s in 0..6 {
                for j in 0..nd {
                    let mut acc = c(0.0, 0.0);
                    for t in 0..6 {
                        acc += (b0[t * nd + j] + beta * b1[t * nd + j] + inv_r * b2[t * nd + j]) * dm[s][t];
                    }
                    db[s * nd + j] = acc;
                }
            }
            for i in 0..nd {
                for j in 0..nd {
                    let mut acc = c(0.0, 0.0);
                    for s in 0..6 {
                        let bh = b0[s * nd + i] + beta_h * b1[s * nd + i] + inv_r * b2[s * nd + i];
                        acc += bh * db[s * nd + j];
                    }
                    ke[i * nd + j] += acc * re.w[q] * r;
                }
            }
        }
        let vol = c(geo.vol[ci], 0.0) * det;
        ke.iter_mut().for_each(|x| *x *= vol);
        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        coo.push_block(&full, &full, &ke);
    }
    Ok(coo.to_csr())
}

/// Barycentric coordinates of quadrature point `q` of the reference element rule.
fn barycentric(re: &RefElement, q: usize, order: usize) -> [f64; 3] {
    // P2 vertex functions are l (2 l - 1); recover l from the vertex shape values by solving
    // N_i = l_i (2 l_i - 1)  ->  not invertible for l < 0.5 on its own, so use the edge functions
    // N_ij = 4 l_i l_j together with sum l = 1 instead. For P1 the shape values are the coordinates.
    if order == 1 {
        return [re.n[q][0], re.n[q][1], re.n[q][2]];
    }
    // P2 ordering: vertices 0..3, then edges (0,1), (0,2), (1,2)
    let nv = &re.n[q];
    // l_i = N_i-vertex + (N_ij + N_ik) / 2 : since N_i = l_i (2 l_i - 1) and N_ij + N_ik = 4 l_i (l_j + l_k) = 4 l_i (1 - l_i),
    // so N_i + (N_ij + N_ik) / 2 = l_i (2 l_i - 1 + 2 - 2 l_i) = l_i
    [nv[0] + (nv[3] + nv[4]) / 2.0, nv[1] + (nv[3] + nv[5]) / 2.0, nv[2] + (nv[4] + nv[5]) / 2.0]
}

/// Mass matrix `rho r` over `3 * n_dofs` unknowns, with the PML stretch determinant `s_r s_z s_theta` as a complex weight.
pub fn assemble_mass_axisym(mesh: &Discretisation, rho: &[f64], order: usize, stretch: Option<&[[C; 3]]>) -> CCsr {
    let re = RefElement::new(2, order, 2 * order + 1);
    let geo = geometry(mesh);
    let nloc = re.nloc;
    let n1 = mesh.n_dofs(order);
    let n = 3 * n1;
    let mut coo = CCoo::new(n, n);
    let (mut scratch, mut full) = (Vec::new(), Vec::new());
    let mut me = vec![c(0.0, 0.0); nloc * nloc];
    let mut block = vec![c(0.0, 0.0); 9 * nloc * nloc];
    for ci in 0..mesh.n_cells() {
        me.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let v = mesh.cell(ci);
        let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
        for q in 0..re.w.len() {
            let lam = barycentric(&re, q, order);
            let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
            for a in 0..nloc {
                for b in 0..nloc {
                    me[a * nloc + b] += c(re.w[q] * re.n[q][a] * re.n[q][b] * r, 0.0);
                }
            }
        }
        // the weight `r` becomes `s_theta r` in a radial layer (third stretch component), as in the stiffness
        let det = stretch.map_or(c(1.0, 0.0), |s| s[ci][0] * s[ci][1] * s[ci][2]);
        let sc = det * (rho[ci] * geo.vol[ci]);
        block.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let nd = 3 * nloc;
        for k in 0..3 {
            for a in 0..nloc {
                for b in 0..nloc {
                    block[(k * nloc + a) * nd + k * nloc + b] = me[a * nloc + b] * sc;
                }
            }
        }
        mesh.cell_dofs(ci, order, &mut scratch);
        full.clear();
        for k in 0..3 {
            full.extend(scratch.iter().map(|&d| d + k * n1));
        }
        coo.push_block(&full, &full, &block[..nd * nd]);
    }
    coo.to_csr()
}
