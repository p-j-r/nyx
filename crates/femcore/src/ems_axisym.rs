//! Electromagnetic eigenmodes of a body of revolution (ring resonators, disks, coaxial cavities) with
//! azimuthal dependence `exp(-i m theta)` on the (r, z) half plane (mesh x = r, mesh z = z; theta is wafer y, so `m > 0` runs along +y): the frequency is
//! sought at a fixed azimuthal order `m`, so `m / (k0 R)` is the effective index of a ring of radius `R`.
//!
//! Unknowns: `(E_r, E_z)` in the same degree-2 Nedelec space as the waveguide mode solver and `E_theta`
//! in P2 Lagrange (the mixed pair is stable for edge elements, so there are no spurious modes except the
//! gradient null space at `k = 0`, far from the shift). The Hermitian pencil is
//! `Int mu^-1 conj(curl F) . curl E r dA = k^2 Int eps conj(F) . E r dA` with
//!
//! ```text
//!   (curl E)_r     = -(i m / r) E_z - dE_theta/dz
//!   (curl E)_theta = dE_r/dz - dE_z/dr
//!   (curl E)_z     = dE_theta/dr + E_theta / r + (i m / r) E_r
//! ```
//!
//! The walls are perfect conductors (tangential `E = 0`); open boundaries are not modelled, which is
//! adequate for well-confined whispering-gallery modes in a large enough domain. The mesh must not
//! touch the axis `r = 0`.

use crate::sbs::fine_rule;
use crate::csparse::{c, eigs_complex, CCoo, CCsr, C};
use crate::ems::{cell_data, edge_basis, tri_g, CellData, C0};
use crate::fem::*;
use crate::mesh::SimplexMesh;
use crate::problem::*;

/// Material data of every cell of an axisymmetric problem (diagonal tensors in the order r, z, theta; complex inside
/// perfectly matched layers).
#[derive(Clone, Debug)]
pub struct AxiMaterial {
    pub eps: Vec<[C; 3]>,
    pub mu: Vec<[C; 3]>,
}

impl AxiMaterial {
    /// Isotropic `eps_r` and `mu_r = 1` per domain (`index`: refractive index per domain label, default 1).
    pub fn from_indices(mesh: &SimplexMesh, index: &std::collections::BTreeMap<i32, f64>) -> AxiMaterial {
        let eps: Vec<[C; 3]> = (0..mesh.n_cells()).map(|ci| { let n = index.get(&mesh.cell_labels[ci]).copied().unwrap_or(1.0); [c(n * n, 0.0); 3] }).collect();
        AxiMaterial { mu: vec![[c(1.0, 0.0); 3]; eps.len()], eps }
    }

    /// Absorbing layers by complex coordinate stretching of `r` and `z` (`Pml::x` is the r direction, `Pml::z` the z direction).
    /// With the stretches `s_r = dr~/dr`, `s_z` and `s_theta = r~/r`, the layer is the material
    /// `eps_i -> eps_i s_1 s_2 s_3 / s_i^2` (and the same for `mu`) in the physical coordinates, `r~` being the integral of
    /// `s_r` from the interface of the layer.
    pub fn with_pml(mut self, mesh: &SimplexMesh, pml: &std::collections::BTreeMap<i32, Pml>) -> AxiMaterial {
        let mut boxes: std::collections::BTreeMap<i32, ([f64; 2], [f64; 2])> = Default::default();
        for ci in 0..mesh.n_cells() {
            let l = mesh.cell_labels[ci];
            if pml.get(&l).is_some_and(|p| p.x != PmlDir::Off || p.z != PmlDir::Off) {
                let e = boxes.entry(l).or_insert(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]));
                for &v in mesh.cell(ci) {
                    let pt = mesh.point(v as usize);
                    for k in 0..2 {
                        e.0[k] = e.0[k].min(pt[k]);
                        e.1[k] = e.1[k].max(pt[k]);
                    }
                }
            }
        }
        for ci in 0..mesh.n_cells() {
            let l = mesh.cell_labels[ci];
            let (Some(p), Some((lo, hi))) = (pml.get(&l), boxes.get(&l)) else { continue };
            let mut cen = [0.0; 2];
            for &v in mesh.cell(ci) {
                for k in 0..2 {
                    cen[k] += mesh.point(v as usize)[k] / 3.0;
                }
            }
            let (sr, sz) = crate::ems_pml::stretch(p, [p.x, p.z], cen, *lo, *hi);
            // r~ at the centre of the cell
            let r = cen[0];
            let len = (hi[0] - lo[0]).max(1e-300);
            let rt = match p.x {
                PmlDir::Off => c(r, 0.0),
                PmlDir::Plus => {
                    let t = ((r - lo[0]) / len).clamp(0.0, 1.0);
                    c(lo[0] + p.re * (r - lo[0]), -p.strength * len * t.powf(p.order + 1.0) / (p.order + 1.0))
                }
                PmlDir::Minus => {
                    let t = ((hi[0] - r) / len).clamp(0.0, 1.0);
                    c(hi[0] - p.re * (hi[0] - r), p.strength * len * t.powf(p.order + 1.0) / (p.order + 1.0))
                }
            };
            let st = rt / r;
            // order r, z, theta
            let lam = [st * sz / sr, sr * st / sz, sr * sz / st];
            for k in 0..3 {
                self.eps[ci][k] *= lam[k];
                self.mu[ci][k] *= lam[k];
            }
        }
        self
    }
}

/// Modes of one azimuthal order.
pub struct AxiModes {
    pub m: f64,
    /// real parts of the frequencies (Hz), ascending in distance from the shift
    pub freqs: Vec<f64>,
    /// complex frequencies (Hz): the imaginary part is the radiation loss (with absorbing layers), `Q = Re / (2 |Im|)`
    pub freqs_c: Vec<C>,
    /// coefficient vectors over all dofs (`net` edge dofs then the P2 dofs of E_theta), zero on the walls
    pub vectors: Vec<Vec<C>>,
    pub order: usize,
    pub converged: bool,
}

/// Fields of one mode at a point.
#[derive(Clone, Copy, Debug, Default)]
pub struct Fields3 {
    pub e: [C; 3],
    pub d: [C; 3],
    pub h: [C; 3],
}

pub(crate) struct AxiCtx<'a> {
    pub mesh: &'a SimplexMesh,
    pub geo: Geometry,
    pub cd: CellData,
    pub order: usize,
}

impl<'a> AxiCtx<'a> {
    pub fn new(mesh: &'a SimplexMesh, order: usize) -> Self {
        AxiCtx { mesh, geo: geometry(mesh), cd: cell_data(mesh, order), order }
    }

    pub fn radius(&self, ci: usize, l: &[f64; 3]) -> f64 {
        let v = self.mesh.cell(ci);
        (0..3).map(|k| l[k] * self.mesh.point(v[k] as usize)[0]).sum()
    }

    /// E (r, z, theta), D and H of a mode (coefficients `x`, azimuthal order `m`, angular frequency `omega`)
    /// at barycentric point `l` of cell `ci`.
    pub fn fields(&self, x: &[C], m: f64, omega: f64, mat: &AxiMaterial, ci: usize, l: &[f64; 3]) -> Fields3 {
        let nloc_e = if self.order == 1 { 3 } else { 8 };
        let mut w = [[0.0; 2]; 8];
        let mut curl = [0.0; 8];
        let g = tri_g(&self.geo, ci);
        edge_basis(l, &g, &self.cd.sign[ci], self.order, &mut w[..nloc_e], &mut curl[..nloc_e]);
        let (mut er, mut ez, mut cth) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        for s in 0..nloc_e {
            let cf = x[self.cd.gdofs[ci][s]];
            er += cf * w[s][0];
            ez += cf * w[s][1];
            cth += cf * curl[s];
        }
        let (vals, dn) = shape(l, self.order);
        let mut pd = Vec::with_capacity(6);
        self.mesh.cell_dofs(ci, self.order, &mut pd);
        let (mut et, mut dtx, mut dty) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
        let gb = ci * 6;
        for a in 0..vals.len() {
            let cf = x[self.cd.net + pd[a]];
            et += cf * vals[a];
            let (mut gx, mut gy) = (0.0, 0.0);
            for b in 0..3 {
                gx += dn[a * 3 + b] * self.geo.g[gb + b * 2];
                gy += dn[a * 3 + b] * self.geo.g[gb + b * 2 + 1];
            }
            dtx += cf * gx;
            dty += cf * gy;
        }
        let r = self.radius(ci, l);
        let im = c(0.0, m);
        // curl components (r, theta, z)
        let cr = -(ez * im) / r - dty;
        let ct = -cth;
        let cz = dtx + et / r + er * im / r;
        let q = c(0.0, 1.0) / (omega * crate::sbs::MU0);
        let mu = mat.mu[ci];
        // H in the order (r, z, theta)
        let h = [cr * q / mu[0], cz * q / mu[1], ct * q / mu[2]];
        let e = [er, ez, et];
        let eps = mat.eps[ci];
        let d = [e[0] * eps[0] * crate::sbs::EPS0, e[1] * eps[1] * crate::sbs::EPS0, e[2] * eps[2] * crate::sbs::EPS0];
        Fields3 { e, d, h }
    }
}

/// Modes with frequency nearest `f_target` at azimuthal order `m` (`nev` of them).
pub fn solve_axisym_modes(mesh: &SimplexMesh, mat: &AxiMaterial, m: f64, f_target: f64, nev: usize, order: usize) -> Result<AxiModes, String> {
    if mesh.dim != 2 {
        return Err("an axisymmetric problem needs a 2D mesh in (r, z)".into());
    }
    if order != 1 && order != 2 {
        return Err("element order must be 1 or 2".into());
    }
    let ctx = AxiCtx::new(mesh, order);
    let (nc, nv, ne) = (mesh.n_cells(), mesh.n_points(), mesh.n_edges());
    let net = ctx.cd.net;
    let nth = mesh.n_dofs(order);
    let n_tot = net + nth;
    let nloc_e = if order == 1 { 3 } else { 8 };
    let nl = mesh.nloc(order);
    let nb = nloc_e + nl;
    let (pts, wts) = fine_rule();
    let mut a_coo = CCoo::new(n_tot, n_tot);
    let mut b_coo = CCoo::new(n_tot, n_tot);
    let mut w = vec![[0.0; 2]; nloc_e];
    let mut curl = vec![0.0; nloc_e];
    let mut pd = Vec::new();
    let mut kl = vec![c(0.0, 0.0); nb * nb];
    let mut ml = vec![c(0.0, 0.0); nb * nb];
    // curl of every basis function at a point: (C_r, C_theta, C_z)
    let mut cb = vec![[c(0.0, 0.0); 3]; nb];
    let mut eb = vec![[0.0; 3]; nb];
    for ci in 0..nc {
        kl.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        ml.iter_mut().for_each(|x| *x = c(0.0, 0.0));
        let g = tri_g(&ctx.geo, ci);
        let area = ctx.geo.vol[ci];
        let (eps, mu) = (mat.eps[ci], mat.mu[ci]);
        for (l, wq) in pts.iter().zip(&wts) {
            let r = ctx.radius(ci, l);
            if r <= 0.0 {
                return Err("the mesh of an axisymmetric problem must lie in r > 0".into());
            }
            edge_basis(l, &g, &ctx.cd.sign[ci], order, &mut w, &mut curl);
            for s in 0..nloc_e {
                cb[s] = [c(0.0, -m / r * w[s][1]), c(-curl[s], 0.0), c(0.0, m / r * w[s][0])];
                eb[s] = [w[s][0], w[s][1], 0.0];
            }
            let (vals, dn) = shape(l, order);
            for a in 0..nl {
                let (mut gx, mut gy) = (0.0, 0.0);
                for b in 0..3 {
                    gx += dn[a * 3 + b] * g[b][0];
                    gy += dn[a * 3 + b] * g[b][1];
                }
                cb[nloc_e + a] = [c(-gy, 0.0), c(0.0, 0.0), c(gx + vals[a] / r, 0.0)];
                eb[nloc_e + a] = [0.0, 0.0, vals[a]];
            }
            let wr = wq * r * area;
            for i in 0..nb {
                for j in 0..nb {
                    // C components are (r, theta, z): mu index (r, z, theta)
                    let kk = cb[i][0].conj() * cb[j][0] / mu[0] + cb[i][1].conj() * cb[j][1] / mu[2] + cb[i][2].conj() * cb[j][2] / mu[1];
                    kl[i * nb + j] += kk * wr;
                    let mm = eps[0] * (eb[i][0] * eb[j][0]) + eps[1] * (eb[i][1] * eb[j][1]) + eps[2] * (eb[i][2] * eb[j][2]);
                    ml[i * nb + j] += mm * wr;
                }
            }
        }
        mesh.cell_dofs(ci, order, &mut pd);
        let mut dofs: Vec<usize> = ctx.cd.gdofs[ci].clone();
        dofs.extend(pd.iter().map(|&p| net + p));
        a_coo.push_block(&dofs, &dofs, &kl);
        b_coo.push_block(&dofs, &dofs, &ml);
    }
    let a_full: CCsr = a_coo.to_csr();
    let b_full: CCsr = b_coo.to_csr();
    // perfectly conducting walls on the exterior facets
    let exterior = mesh.exterior_facets();
    let mut be: Vec<usize> = (0..mesh.n_facets()).filter(|&f| exterior[f]).map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1])).collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets()).filter(|&f| exterior[f]).flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>()).collect();
    bn.sort_unstable();
    bn.dedup();
    let mut fixed: Vec<usize> = be.clone();
    if order == 2 {
        fixed.extend(be.iter().map(|&e| ne + e));
    }
    fixed.extend(bn.iter().map(|&v| net + v));
    if order == 2 {
        fixed.extend(be.iter().map(|&e| net + nv + e));
    }
    let free = free_dofs(n_tot, &fixed);
    let af = a_full.submatrix(&free, &free);
    let bf = b_full.submatrix(&free, &free);
    let k0 = 2.0 * std::f64::consts::PI * f_target / C0;
    let pairs = eigs_complex(&af, &bf, c(k0 * k0, 0.0), nev, 1e-9, 800)?;
    let mut modes: Vec<(C, Vec<C>)> = vec![];
    for (i, lam) in pairs.values.iter().enumerate() {
        if lam.re <= 0.0 || pairs.residuals.get(i).is_some_and(|r| *r > 1e-6) {
            continue;
        }
        let f = lam.sqrt() * (C0 / (2.0 * std::f64::consts::PI));
        let f = if f.re < 0.0 { -f } else { f };
        let mut x = vec![c(0.0, 0.0); n_tot];
        for (a, &dof) in free.iter().enumerate() {
            x[dof] = pairs.vectors[i][a];
        }
        modes.push((f, x));
    }
    modes.sort_by(|a, b| (a.0.re - f_target).abs().partial_cmp(&(b.0.re - f_target).abs()).unwrap());
    if modes.is_empty() {
        return Err("no mode converged near the target frequency".into());
    }
    Ok(AxiModes { m, freqs: modes.iter().map(|x| x.0.re).collect(), freqs_c: modes.iter().map(|x| x.0).collect(), vectors: modes.into_iter().map(|x| x.1).collect(), order, converged: pairs.converged })
}

#[allow(dead_code)]
fn _unused(_: &Physics) {}

/// The electric field of mode `j` at the quadrature points of every cell outside the absorbing layers (`in_pml`), as the samples of a
/// polarization analysis: `ex` = E_r (the in-plane, "TE-like" direction), `ey` = E_z (normal to the plane of the ring), `ez` = E_theta (the
/// direction of propagation). Each point carries the weight `r * area` of the volume of revolution it stands for. The triple (r, z, theta)
/// is left-handed, so the analysis of a ring reads the sense of rotation with `TimeConvention::MinusI`.
pub fn mode_samples(mesh: &SimplexMesh, mat: &AxiMaterial, modes: &AxiModes, j: usize, in_pml: &[bool]) -> Vec<crate::polarization::Sample> {
    let ctx = AxiCtx::new(mesh, modes.order);
    let (pts, wts) = crate::sbs::fine_rule();
    let omega = 2.0 * std::f64::consts::PI * modes.freqs[j];
    let mut out = Vec::with_capacity(mesh.n_cells() * pts.len());
    for ci in 0..mesh.n_cells() {
        if in_pml.get(ci).copied().unwrap_or(false) {
            continue;
        }
        for (l, w) in pts.iter().zip(&wts) {
            let f = ctx.fields(&modes.vectors[j], modes.m, omega, mat, ci, l);
            let cx = |z: C| [z.re, z.im];
            out.push(crate::polarization::Sample { ex: cx(f.e[0]), ey: cx(f.e[1]), ez: cx(f.e[2]), w: w * ctx.geo.vol[ci] * ctx.radius(ci, l) });
        }
    }
    out
}

/// The power `P = 2 Int Re(conj(E_z) H_r - conj(E_r) H_z) dA` of mode `j` along `+theta` (W for the amplitudes of the mode), with the
/// fields of the time convention `e^{j omega t}` (`H = j curl E / (omega mu)`). With the azimuthal dependence `exp(-i m theta)` it is positive
/// for `m > 0` and negative for `m < 0`: the wave runs along `+theta`, wafer `+y`.
pub fn power_along_theta(mesh: &SimplexMesh, mat: &AxiMaterial, modes: &AxiModes, j: usize) -> f64 {
    let ctx = AxiCtx::new(mesh, modes.order);
    let (pts, wts) = fine_rule();
    let omega = 2.0 * std::f64::consts::PI * modes.freqs[j];
    let mut p = 0.0;
    for ci in 0..mesh.n_cells() {
        for (l, w) in pts.iter().zip(&wts) {
            let f = ctx.fields(&modes.vectors[j], modes.m, omega, mat, ci, l);
            p += 2.0 * (f.e[1].conj() * f.h[0] - f.e[0].conj() * f.h[1]).re * w * ctx.geo.vol[ci];
        }
    }
    p
}
