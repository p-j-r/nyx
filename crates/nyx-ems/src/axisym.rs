// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-nyx-Commercial
// Copyright (C) 2026 p-j-r
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

use crate::modes::C0;
use crate::problem::{Problem, Study};
use femcore::edge::{cell_data, edge_basis, tri_g, CellData};
use femcore::fem::*;
use femcore::fields::{local_bary, node_average, CellFields};
use femcore::linalg::csr::{c, CCoo, CCsr, C};
use femcore::linalg::driver::eigs_complex;
use femcore::optical::{fine_rule, EPS0, MU0};
use femcore::pml::{radial_coordinate, stretch};
use femcore::props::{Materials, Pml, PmlDir};
use nyx_physics::{cell_fields_of, Discretisation, Field, Scalar, Solution, SolveResult};
use std::time::Instant;

/// Material data of every cell of an axisymmetric problem (diagonal tensors in the order r, z, theta; complex inside
/// perfectly matched layers).
#[derive(Clone, Debug)]
pub struct AxiMaterial {
    pub eps: Vec<[C; 3]>,
    pub mu: Vec<[C; 3]>,
}

impl AxiMaterial {
    /// Isotropic `eps_r` and `mu_r = 1` per domain (`index`: refractive index per domain label, default 1).
    pub fn from_indices(mesh: &Discretisation, index: &std::collections::BTreeMap<i32, f64>) -> AxiMaterial {
        let eps: Vec<[C; 3]> = (0..mesh.n_cells())
            .map(|ci| {
                let n = index.get(&mesh.cell_labels[ci]).copied().unwrap_or(1.0);
                [c(n * n, 0.0); 3]
            })
            .collect();
        AxiMaterial { mu: vec![[c(1.0, 0.0); 3]; eps.len()], eps }
    }

    /// `eps_r` and `mu_r` per domain from the problem's materials (default 1), absorbing layers (`Materials::pml`) included.
    pub fn from_materials(mesh: &Discretisation, m: &Materials) -> Result<AxiMaterial, String> {
        let eps = m.cell_values(mesh, "eps_r", Some(1.0))?;
        let mu = m.cell_values(mesh, "mu_r", Some(1.0))?;
        let mat = AxiMaterial {
            eps: eps.iter().map(|&e| [c(e, 0.0); 3]).collect(),
            mu: mu.iter().map(|&u| [c(u, 0.0); 3]).collect(),
        };
        Ok(mat.with_pml(mesh, &m.pml))
    }

    /// Absorbing layers by complex coordinate stretching of `r` and `z` (`Pml::x` is the r direction, `Pml::z` the z direction).
    /// With the stretches `s_r = dr~/dr`, `s_z` and `s_theta = r~/r`, the layer is the material
    /// `eps_i -> eps_i s_1 s_2 s_3 / s_i^2` (and the same for `mu`) in the physical coordinates, `r~` being the integral of
    /// `s_r` from the interface of the layer.
    pub fn with_pml(mut self, mesh: &Discretisation, pml: &std::collections::BTreeMap<i32, Pml>) -> AxiMaterial {
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
            let (sr, sz) = stretch(p, [p.x, p.z], cen, *lo, *hi);
            // r~ at the centre of the cell
            let r = cen[0];
            let rt = radial_coordinate(p, r, lo[0], hi[0]);
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
    pub mesh: &'a Discretisation,
    pub geo: Geometry,
    pub cd: CellData,
    pub order: usize,
}

impl<'a> AxiCtx<'a> {
    pub fn new(mesh: &'a Discretisation, order: usize) -> Self {
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
        let q = c(0.0, 1.0) / (omega * MU0);
        let mu = mat.mu[ci];
        // H in the order (r, z, theta)
        let h = [cr * q / mu[0], cz * q / mu[1], ct * q / mu[2]];
        let e = [er, ez, et];
        let eps = mat.eps[ci];
        let d = [e[0] * eps[0] * EPS0, e[1] * eps[1] * EPS0, e[2] * eps[2] * EPS0];
        Fields3 { e, d, h }
    }
}

/// Modes with frequency nearest `f_target` at azimuthal order `m` (`nev` of them).
pub fn solve_axisym_modes(
    mesh: &Discretisation,
    mat: &AxiMaterial,
    m: f64,
    f_target: f64,
    nev: usize,
    order: usize,
) -> Result<AxiModes, String> {
    let asm = assemble_axisym(mesh, mat, m, order)?;
    solve_axisym_pairs(&asm, f_target, nev)
}

/// Stage 1: the curl-curl matrix `A` and the permittivity mass `B` of `A x = (omega/c)^2 B x` on the free dofs (edge
/// functions for E_r, E_z, then nodal E_theta), with perfectly conducting exterior walls.
pub struct AxiAssembly {
    pub m: f64,
    pub order: usize,
    pub a: CCsr,
    pub b: CCsr,
    pub free: Vec<usize>,
    pub n_tot: usize,
}

pub fn assemble_axisym(mesh: &Discretisation, mat: &AxiMaterial, m: f64, order: usize) -> Result<AxiAssembly, String> {
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
                    let kk = cb[i][0].conj() * cb[j][0] / mu[0]
                        + cb[i][1].conj() * cb[j][1] / mu[2]
                        + cb[i][2].conj() * cb[j][2] / mu[1];
                    kl[i * nb + j] += kk * wr;
                    let mm = eps[0] * (eb[i][0] * eb[j][0])
                        + eps[1] * (eb[i][1] * eb[j][1])
                        + eps[2] * (eb[i][2] * eb[j][2]);
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
    let mut be: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| exterior[f])
        .map(|f| mesh.edge_index(mesh.facet(f)[0], mesh.facet(f)[1]))
        .collect();
    be.sort_unstable();
    be.dedup();
    let mut bn: Vec<usize> = (0..mesh.n_facets())
        .filter(|&f| exterior[f])
        .flat_map(|f| mesh.facet(f).iter().map(|&v| v as usize).collect::<Vec<_>>())
        .collect();
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
    Ok(AxiAssembly { m, order, a: a_full.submatrix(&free, &free), b: b_full.submatrix(&free, &free), free, n_tot })
}

/// Stage 2: the `nev` eigenpairs nearest `f_target` (Hz), as complex frequencies and field vectors over all dofs,
/// sorted by closeness to it; unconverged pairs are dropped.
pub fn solve_axisym_pairs(asm: &AxiAssembly, f_target: f64, nev: usize) -> Result<AxiModes, String> {
    let (m, order, free, n_tot) = (asm.m, asm.order, &asm.free, asm.n_tot);
    let k0 = 2.0 * std::f64::consts::PI * f_target / C0;
    let pairs = eigs_complex(&asm.a, &asm.b, c(k0 * k0, 0.0), nev, 1e-9, 800)?;
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
    Ok(AxiModes {
        m,
        freqs: modes.iter().map(|x| x.0.re).collect(),
        freqs_c: modes.iter().map(|x| x.0).collect(),
        vectors: modes.into_iter().map(|x| x.1).collect(),
        order,
        converged: pairs.converged,
    })
}

/// The power `P = 2 Int Re(conj(E_z) H_r - conj(E_r) H_z) dA` of mode `j` along `+theta` (W for the amplitudes of the mode), with the
/// fields of the time convention `e^{j omega t}` (`H = j curl E / (omega mu)`). With the azimuthal dependence `exp(-i m theta)` it is positive
/// for `m > 0` and negative for `m < 0`: the wave runs along `+theta`, wafer `+y`.
pub fn power_along_theta(mesh: &Discretisation, mat: &AxiMaterial, modes: &AxiModes, j: usize) -> f64 {
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

/// Stage 3: the electric field of mode `j` at the local nodes of every cell, as plot fields (`|E|`, `Re Ex`, `Re Ey`, `Re Ez`, wafer names: x is
/// r, z is z, y is theta; the phase is rotated so that the largest component is real; scaled to max |E| = 1) and as the complex components per
/// cell in the same phase and scale.
pub fn axisym_fields_cells(
    mesh: &Discretisation,
    mat: &AxiMaterial,
    modes: &AxiModes,
    j: usize,
) -> (Vec<Field>, CellFields) {
    let order = modes.order;
    let ctx = AxiCtx::new(mesh, order);
    let (nc, nloc) = (mesh.n_cells(), mesh.nloc(order));
    let lam = local_bary(2, order);
    let omega = 2.0 * std::f64::consts::PI * modes.freqs[j];
    let mut e = vec![[c(0.0, 0.0); 3]; nc * nloc];
    for ci in 0..nc {
        for (a, l) in lam.iter().enumerate() {
            e[ci * nloc + a] = ctx.fields(&modes.vectors[j], modes.m, omega, mat, ci, &[l[0], l[1], l[2]]).e;
        }
    }
    let big =
        e.iter().flat_map(|v| v.iter()).cloned().fold(c(0.0, 0.0), |m, z| if z.norm() > m.norm() { z } else { m });
    let rot = if big.norm() > 0.0 { big.conj() / big.norm() } else { c(1.0, 0.0) };
    let mag: Vec<f64> = e.iter().map(|v| (v[0].norm_sqr() + v[1].norm_sqr() + v[2].norm_sqr()).sqrt()).collect();
    let s = 1.0 / mag.iter().cloned().fold(0.0, f64::max).max(1e-300);
    let avg = |vals: Vec<f64>| node_average(mesh, &vals.iter().map(|x| x * s).collect::<Vec<f64>>(), 1, order);
    // the components are stored in the order (r, z, theta) = wafer (x, z, y)
    let fields = vec![
        Field::scalar("|E|", "a.u.", avg(mag)),
        Field::scalar("Re Ex", "a.u.", avg(e.iter().map(|v| (v[0] * rot).re).collect())),
        Field::scalar("Re Ey", "a.u.", avg(e.iter().map(|v| (v[2] * rot).re).collect())),
        Field::scalar("Re Ez", "a.u.", avg(e.iter().map(|v| (v[1] * rot).re).collect())),
    ];
    let mut cells = CellFields::new(mesh, order);
    for (k, name) in [(0, "Ex"), (2, "Ey"), (1, "Ez")] {
        let z: Vec<C> = e.iter().map(|v| v[k] * rot * s).collect();
        cells.push_complex(name, "a.u.", z.iter().map(|x| x.re).collect(), z.iter().map(|x| x.im).collect());
    }
    (fields, cells)
}

/// The modes of a body of revolution for an `AxisymmetricModes` study: the three stages, one solution per mode with its frequency (and Q when
/// absorbing layers make it complex).
pub fn solve_axisymmetric(mesh: &Discretisation, problem: &Problem) -> Result<SolveResult, String> {
    let t0 = Instant::now();
    let Study::AxisymmetricModes { f_search, azimuthal_order, n_modes } = problem.study.clone() else {
        return Err("an axisymmetric problem needs an AxisymmetricModes study".into());
    };
    if f_search <= 0.0 {
        return Err("the search frequency must be positive".into());
    }
    if n_modes < 1 {
        return Err("at least one mode is needed".into());
    }
    let mat = AxiMaterial::from_materials(mesh, &problem.materials)?;
    let asm = assemble_axisym(mesh, &mat, azimuthal_order as f64, problem.order)?;
    let modes = solve_axisym_pairs(&asm, f_search, n_modes)?;
    let mut solutions = vec![];
    let mut scalars = vec![];
    for (j, f) in modes.freqs_c.iter().enumerate() {
        let (fields, cells) = axisym_fields_cells(mesh, &mat, &modes, j);
        let lossy = f.im.abs() > 1e-12 * f.re.abs().max(1.0);
        let label = if lossy {
            format!("mode {}: {:.6e} {:+.3e}i Hz (Q = {:.3e})", j + 1, f.re, f.im, f.re / (2.0 * f.im.abs()))
        } else {
            format!("mode {}: {:.6} GHz", j + 1, f.re / 1e9)
        };
        solutions.push(Solution {
            label,
            fields,
            cell_fields: cell_fields_of(&cells),
            value: Some(f.re),
            ..Default::default()
        });
        if lossy {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Text(format!("{:.8e} {:+.4e}i", f.re, f.im))));
            scalars.push((format!("Q{}", j + 1), Scalar::Num(f.re / (2.0 * f.im.abs()))));
        } else {
            scalars.push((format!("f{} (Hz)", j + 1), Scalar::Num(f.re)));
        }
    }
    let mut log = vec![format!("azimuthal order {azimuthal_order}, {} dofs, shift {:.4e} Hz", asm.n_tot, f_search)];
    if !modes.converged {
        log.push("warning: the eigensolver did not converge every mode to the requested tolerance".into());
    }
    let result = SolveResult {
        physics: "ems".into(),
        study: "axisymmetric_modes".into(),
        order: problem.order as u8,
        dim: 2,
        solutions,
        scalars,
        log,
        seconds: t0.elapsed().as_secs_f64(),
        dofs: asm.n_tot,
    };
    Ok(result)
}
