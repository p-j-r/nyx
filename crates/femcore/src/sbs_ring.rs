//! Brillouin gain of a ring resonator from its axisymmetric cross-section: whispering-gallery pump and
//! Stokes modes (`ems_axisym`, azimuthal orders `m_p`, `m_s`) and the elastic modes of the ring with the
//! phase-matched azimuthal order `M = m_p + m_s` (backward: the Stokes mode runs the other way round) or
//! `M = m_p - m_s` (forward), from `sms_kz::solve_solid_axisym`.
//!
//! The overlap integrals are those of `sbs` with cylindrical fields and strains; every integral over the
//! cross-section carries the weight `r / r_ref` so that the quantities are per unit length of the ring at the
//! reference radius `r_ref` (mean radius of the mode), the optical power `P = 2 Int Re(conj(E_z) H_r - conj(E_r) H_z) dA`
//! is the power crossing the cross-section, and `vg = P / W` is a speed at `r_ref`. The gain, `g0` and the quality
//! factor then follow the conventions of `sbs` with the length `2 pi r_ref`. A crystal's photoelastic
//! tensor and stiffness vary with the angle around the ring; use tensors with the symmetry of the ring or
//! averaged about the axis (`PhotoElastic::average_about_axis`, `sms_axisym::average_about_axis`).

use crate::sbs::*;
use crate::csparse::{c, C};
use crate::ems::C0;
use crate::ems_axisym::*;
use crate::fem::*;
use crate::mesh::SimplexMesh;
use crate::sms_kz::KzModes;
use std::collections::{BTreeMap, BTreeSet};

impl PhotoElastic {
    /// The tensor averaged over rotations of the crystal about the axis `axis` of the (r, z, theta) frame
    /// (index 1 is the ring axis): what a ring sees when the crystal axes stay fixed while the local frame turns.
    pub fn average_about_axis(&self, axis: usize, samples: usize) -> PhotoElastic {
        let mut avg = [[[[0.0; 3]; 3]; 3]; 3];
        let (a1, a2) = ((axis + 1) % 3, (axis + 2) % 3);
        for sidx in 0..samples {
            let th = 2.0 * std::f64::consts::PI * sidx as f64 / samples as f64;
            let (sn, cs) = th.sin_cos();
            let mut r = [[0.0; 3]; 3];
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
                                            s += r[i][a] * r[j][b] * r[k][cc] * r[l][d] * self.p[a][b][cc][d];
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
        PhotoElastic { p: avg }
    }
}

const VP_ORDER: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];

#[derive(Clone, Debug)]
struct PeSample {
    dofs: [u32; 6],
    nv: [f64; 6],
    grad: [[f64; 2]; 6],
    r: f64,
    /// weight: quadrature weight times area times r / r_ref
    w: f64,
    a: [C; 6],
}

#[derive(Clone, Debug)]
struct MbSample {
    dofs: [u32; 6],
    nv: [f64; 6],
    w: f64,
    n: [f64; 2],
    m: C,
}

/// The optical part of the ring coupling.
pub struct RingCoupling {
    pe: Vec<PeSample>,
    mb: Vec<MbSample>,
    pub omega_pump: f64,
    pub omega_stokes: f64,
    pub m_pump: f64,
    pub m_stokes: f64,
    pub p_pump: f64,
    pub p_stokes: f64,
    pub vg_pump: f64,
    pub vg_stokes: f64,
    pub backward: bool,
    pub r_ref: f64,
    norm: f64,
    order: usize,
}

impl RingCoupling {
    /// Azimuthal order of the phase-matched elastic mode.
    pub fn elastic_order(&self) -> f64 {
        if self.backward { self.m_pump + self.m_stokes } else { self.m_pump - self.m_stokes }
    }
    pub fn group_index(&self) -> (f64, f64) {
        (C0 / self.vg_pump.abs(), C0 / self.vg_stokes.abs())
    }
}

pub struct RingInput<'a> {
    pub mesh: &'a SimplexMesh,
    pub material: &'a AxiMaterial,
    pub pump: (&'a AxiModes, usize),
    pub stokes: (&'a AxiModes, usize),
    pub backward: bool,
    pub photoelastic: &'a BTreeMap<i32, PhotoElastic>,
    pub moving: &'a BTreeSet<i32>,
    /// reference radius for the per-length quantities
    pub r_ref: f64,
}

pub fn ring_coupling(inp: &RingInput) -> Result<RingCoupling, String> {
    let mesh = inp.mesh;
    let (mp, jp) = inp.pump;
    let (ms, js) = inp.stokes;
    let order = mp.order;
    if ms.order != order {
        return Err("pump and Stokes modes must use the same element order".into());
    }
    let ctx = AxiCtx::new(mesh, order);
    let mat = inp.material;
    let (xp, xs) = (&mp.vectors[jp], &ms.vectors[js]);
    let (omega_p, omega_s) = (2.0 * std::f64::consts::PI * mp.freqs[jp], 2.0 * std::f64::consts::PI * ms.freqs[js]);
    let sgn = if inp.backward { -1.0 } else { 1.0 };
    // the Stokes mode running the other way round: E_theta, D_theta and H_r, H_z change sign
    let flip = |mut f: Fields3| -> Fields3 {
        f.e[2] *= sgn;
        f.d[2] *= sgn;
        f.h[0] *= sgn;
        f.h[1] *= sgn;
        f
    };
    let (pts, wts) = fine_rule();
    let nloc = mesh.nloc(order);
    let nc = mesh.n_cells();
    let (mut pw_p, mut pw_s, mut en_p, mut en_s) = (0.0, 0.0, 0.0, 0.0);
    for ci in 0..nc {
        let area = ctx.geo.vol[ci];
        for (l, w) in pts.iter().zip(&wts) {
            let r = ctx.radius(ci, l);
            let wr = r / inp.r_ref;
            let ep = ctx.fields(xp, mp.m, omega_p, mat, ci, l);
            let es = flip(ctx.fields(xs, ms.m, omega_s, mat, ci, l));
            // power along theta: (E x H*)_theta = E_z H_r* - E_r H_z*
            let sz = |f: &Fields3| (f.e[1].conj() * f.h[0] - f.e[0].conj() * f.h[1]).re;
            pw_p += 2.0 * sz(&ep) * w * area;
            pw_s += 2.0 * sz(&es) * w * area;
            let eps = mat.eps[ci];
            let en = |f: &Fields3| (f.e[0].norm_sqr() * eps[0].re + f.e[1].norm_sqr() * eps[1].re + f.e[2].norm_sqr() * eps[2].re) * EPS0;
            en_p += 2.0 * en(&ep) * w * area * wr;
            en_s += 2.0 * en(&es) * w * area * wr;
        }
    }
    if pw_p.abs() < 1e-300 || pw_s.abs() < 1e-300 {
        return Err("an optical mode carries no power".into());
    }
    let norm = 1.0 / (pw_p * pw_s).abs().sqrt();
    let mut pe = vec![];
    let mut pd = Vec::with_capacity(6);
    for ci in 0..nc {
        let Some(pt) = inp.photoelastic.get(&mesh.cell_labels[ci]) else { continue };
        mesh.cell_dofs(ci, order, &mut pd);
        let area = ctx.geo.vol[ci];
        let pref = -1.0 / EPS0;
        for (l, w) in pts.iter().zip(&wts) {
            let r = ctx.radius(ci, l);
            let ep = ctx.fields(xp, mp.m, omega_p, mat, ci, l);
            let es = flip(ctx.fields(xs, ms.m, omega_s, mat, ci, l));
            let mut a = [c(0.0, 0.0); 6];
            for (m, &(k, ll)) in VP_ORDER.iter().enumerate() {
                let mut acc = c(0.0, 0.0);
                for aa in 0..3 {
                    for bb in 0..3 {
                        let coef = pt.p[aa][bb][k][ll];
                        if coef != 0.0 {
                            acc += es.d[aa].conj() * ep.d[bb] * coef;
                        }
                    }
                }
                a[m] = acc * if k == ll { pref } else { 2.0 * pref };
            }
            let (vals, dn) = shape(l, order);
            let mut sample = PeSample { dofs: [0; 6], nv: [0.0; 6], grad: [[0.0; 2]; 6], r, w: w * area * r / inp.r_ref, a };
            for aa in 0..nloc {
                sample.dofs[aa] = pd[aa] as u32;
                sample.nv[aa] = vals[aa];
                for k in 0..2 {
                    let mut gk = 0.0;
                    for b in 0..3 {
                        gk += dn[aa * 3 + b] * ctx.geo.g[ci * 6 + b * 2 + k];
                    }
                    sample.grad[aa][k] = gk;
                }
            }
            pe.push(sample);
        }
    }
    // moving boundary on edges between a moving domain and its surroundings
    let mut edge_cells: Vec<Vec<(u32, u8)>> = vec![vec![]; mesh.n_edges()];
    for ci in 0..nc {
        for (loc, &e) in mesh.cell_edge_ids(ci).iter().enumerate() {
            edge_cells[e as usize].push((ci as u32, loc as u8));
        }
    }
    const LOC: [(usize, usize); 3] = [(0, 1), (0, 2), (1, 2)];
    let gx = [0.5 - 0.5 * 0.861_136_311_594_052_6, 0.5 - 0.5 * 0.339_981_043_584_856_3, 0.5 + 0.5 * 0.339_981_043_584_856_3, 0.5 + 0.5 * 0.861_136_311_594_052_6];
    let gw = [0.5 * 0.347_854_845_137_453_9, 0.5 * 0.652_145_154_862_546_1, 0.5 * 0.652_145_154_862_546_1, 0.5 * 0.347_854_845_137_453_9];
    let mut mb = vec![];
    for cells in &edge_cells {
        if cells.len() != 2 {
            continue;
        }
        let (a, b) = {
            let (c0, c1) = (cells[0], cells[1]);
            match (inp.moving.contains(&mesh.cell_labels[c0.0 as usize]), inp.moving.contains(&mesh.cell_labels[c1.0 as usize])) {
                (true, false) => (c0, c1),
                (false, true) => (c1, c0),
                _ => continue,
            }
        };
        let (ca, cb) = (a.0 as usize, b.0 as usize);
        let la = LOC[a.1 as usize];
        let (va, vb) = (mesh.cell(ca), mesh.cell(cb));
        let (v0, v1) = (va[la.0], va[la.1]);
        let (p0, p1) = (mesh.point(v0 as usize), mesh.point(v1 as usize));
        let len = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2)).sqrt();
        let t = [(p1[0] - p0[0]) / len, (p1[1] - p0[1]) / len];
        let pos_b = |v: u32| vb.iter().position(|&x| x == v).unwrap();
        let (ib0, ib1) = (pos_b(v0), pos_b(v1));
        let cen = |cc: usize| {
            let mut s = [0.0; 2];
            for &v in mesh.cell(cc) {
                for k in 0..2 {
                    s[k] += mesh.point(v as usize)[k] / 3.0;
                }
            }
            s
        };
        let (cen_a, cen_b) = (cen(ca), cen(cb));
        let mut n = [t[1], -t[0]];
        if n[0] * (cen_b[0] - cen_a[0]) + n[1] * (cen_b[1] - cen_a[1]) < 0.0 {
            n = [-n[0], -n[1]];
        }
        let diag = |e: [C; 3]| -> [[C; 3]; 3] { [[e[0], c(0.0, 0.0), c(0.0, 0.0)], [c(0.0, 0.0), e[1], c(0.0, 0.0)], [c(0.0, 0.0), c(0.0, 0.0), e[2]]] };
        let blocks = mb_blocks(&diag(mat.eps[ca]), &diag(mat.eps[cb]), t, n);
        for q in 0..4 {
            let s = gx[q];
            let mut l_a = [0.0; 3];
            l_a[la.0] = 1.0 - s;
            l_a[la.1] = s;
            let mut l_b = [0.0; 3];
            l_b[ib0] = 1.0 - s;
            l_b[ib1] = s;
            let avg = |x: &[C], m: f64, om: f64| -> Fields3 {
                let (f1, f2) = (ctx.fields(x, m, om, mat, ca, &l_a), ctx.fields(x, m, om, mat, cb, &l_b));
                let mut r = Fields3::default();
                for k in 0..3 {
                    r.e[k] = (f1.e[k] + f2.e[k]) * 0.5;
                    r.d[k] = (f1.d[k] + f2.d[k]) * 0.5;
                    r.h[k] = (f1.h[k] + f2.h[k]) * 0.5;
                }
                r
            };
            let ep = avg(xp, mp.m, omega_p);
            let es = flip(avg(xs, ms.m, omega_s));
            // in-plane (r, z) components are 0 and 1, theta is 2
            let m = mb_integrand(&blocks, &es.e, &es.d, &ep.e, &ep.d, t, n);
            let (vals, _) = shape(&l_a, order);
            let mut pd_a = Vec::with_capacity(6);
            mesh.cell_dofs(ca, order, &mut pd_a);
            let r = ctx.radius(ca, &l_a);
            let mut sample = MbSample { dofs: [0; 6], nv: [0.0; 6], w: gw[q] * len * r / inp.r_ref, n, m };
            for aa in 0..nloc {
                sample.dofs[aa] = pd_a[aa] as u32;
                sample.nv[aa] = vals[aa];
            }
            mb.push(sample);
        }
    }
    Ok(RingCoupling {
        pe,
        mb,
        omega_pump: omega_p,
        omega_stokes: omega_s,
        m_pump: mp.m,
        m_stokes: ms.m,
        p_pump: pw_p,
        p_stokes: pw_s,
        vg_pump: pw_p / en_p,
        vg_stokes: pw_s / en_s,
        backward: inp.backward,
        r_ref: inp.r_ref,
        norm,
        order,
    })
}

/// Gains of every elastic mode of `mech` (from `sms_kz::solve_solid_axisym` at order `cp.elastic_order()`).
pub fn ring_gains(mesh: &SimplexMesh, cp: &RingCoupling, mech: &KzModes, opts: &GainOptions) -> Vec<GainRow> {
    let order = cp.order;
    let n1 = mesh.n_dofs(order);
    let nloc = mesh.nloc(order);
    let m_az = mech.kz;
    let q = quad_rule(2, 4);
    let geo = geometry(mesh);
    let vals: Vec<Vec<f64>> = q.lam.iter().map(|l| shape(l, order).0).collect();
    let nc = mesh.n_cells();
    let mut pd = Vec::with_capacity(6);
    let mut rows = vec![];
    let length = if opts.length > 0.0 { opts.length } else { 2.0 * std::f64::consts::PI * cp.r_ref };
    for (j, fz) in mech.freqs.iter().enumerate() {
        let u = &mech.vectors[j];
        let (mut mass, mut body, mut pmlm) = (0.0, 0.0, 0.0);
        for ci in 0..nc {
            let rho = mech.rho[ci];
            if rho <= 0.0 {
                continue;
            }
            mesh.cell_dofs(ci, order, &mut pd);
            let v = mesh.cell(ci);
            let xs = [mesh.point(v[0] as usize)[0], mesh.point(v[1] as usize)[0], mesh.point(v[2] as usize)[0]];
            let mut cell_mass = 0.0;
            for (qi, w) in q.w.iter().enumerate() {
                let lam = &q.lam[qi];
                let r = lam[0] * xs[0] + lam[1] * xs[1] + lam[2] * xs[2];
                let mut s = 0.0;
                for k in 0..3 {
                    let mut vv = c(0.0, 0.0);
                    for a in 0..nloc {
                        vv += u[k * n1 + pd[a]] * vals[qi][a];
                    }
                    s += vv.norm_sqr();
                }
                cell_mass += s * w * r / cp.r_ref;
            }
            let m = cell_mass * rho * geo.vol[ci];
            mass += m;
            let lab = mesh.cell_labels[ci];
            if opts.body_domains.contains(&lab) {
                body += m;
            }
            if opts.pml_domains.contains(&lab) {
                pmlm += m;
            }
        }
        let mut f_pe = c(0.0, 0.0);
        for s in &cp.pe {
            let (mut urr, mut urz, mut uzr, mut uzz, mut utr, mut utz) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            let (mut ur, mut uz, mut ut) = (c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nloc {
                let d = s.dofs[a] as usize;
                let (ar, az, at) = (u[d], u[n1 + d], u[2 * n1 + d]);
                urr += ar * s.grad[a][0];
                urz += ar * s.grad[a][1];
                uzr += az * s.grad[a][0];
                uzz += az * s.grad[a][1];
                utr += at * s.grad[a][0];
                utz += at * s.grad[a][1];
                ur += ar * s.nv[a];
                uz += az * s.nv[a];
                ut += at * s.nv[a];
            }
            let im = c(0.0, m_az);
            let inv_r = 1.0 / s.r;
            // tensor strains in the order rr, zz, tt, zt, rt, rz
            let e_rr = urr;
            let e_zz = uzz;
            let e_tt = (ur - im * ut) * inv_r;
            let e_zt = (utz - im * uz * inv_r) * 0.5;
            let e_rt = (utr - ut * inv_r - im * ur * inv_r) * 0.5;
            let e_rz = (urz + uzr) * 0.5;
            let strain = [e_rr, e_zz, e_tt, e_zt, e_rt, e_rz];
            let mut acc = c(0.0, 0.0);
            for m in 0..6 {
                let e = if opts.strain == StrainConvention::Conjugate { strain[m].conj() } else { strain[m] };
                acc += s.a[m] * e;
            }
            f_pe += acc * s.w;
        }
        let mut f_mb = c(0.0, 0.0);
        for s in &cp.mb {
            let (mut ur, mut uz) = (c(0.0, 0.0), c(0.0, 0.0));
            for a in 0..nloc {
                let d = s.dofs[a] as usize;
                ur += u[d] * s.nv[a];
                uz += u[n1 + d] * s.nv[a];
            }
            let un = ur * s.n[0] + uz * s.n[1];
            f_mb += s.m * if opts.mb_conjugate { un.conj() } else { un } * s.w;
        }
        f_pe *= cp.norm;
        f_mb *= cp.norm;
        let omega_m = 2.0 * std::f64::consts::PI * fz.re;
        let inv_rad = if fz.im.abs() > 0.0 { 2.0 * fz.im.abs() / fz.re } else { 0.0 };
        let inv_q = inv_rad + opts.inv_q_extra.get(j).copied().unwrap_or(0.0);
        let q_eig = if inv_q > 0.0 { 1.0 / inv_q } else { f64::INFINITY };
        let q_used = q_eig.min(opts.q_cap);
        let k_eff = opts.mass_factor * mass * omega_m * omega_m;
        let gain = |f: C| 2.0 * cp.omega_pump * (q_used / k_eff) * f.norm_sqr();
        let g0 = |f: C, g: f64| -> f64 {
            let v = (cp.vg_pump * cp.vg_stokes).abs() * (HBAR * cp.omega_pump) * omega_m / (4.0 * length) * (g / q_used);
            v.sqrt() * if f.re > 0.0 { -1.0 } else { 1.0 }
        };
        let (gm, gp, gt) = (gain(f_mb), gain(f_pe), gain(f_pe + f_mb));
        let (g0m, g0p) = (g0(f_mb, gm), g0(f_pe, gp));
        rows.push(GainRow {
            mode: j,
            freq: *fz,
            q_eig,
            q_used,
            f_mb,
            f_pe,
            gain_mb: gm,
            gain_pe: gp,
            gain_total: gt,
            g0_mb: g0m,
            g0_pe: g0p,
            g0_total: (g0m + g0p).abs(),
            frac_body: body / mass.max(1e-300),
            frac_pml: pmlm / mass.max(1e-300),
            m_eff: mass,
            parity: None,
        });
    }
    rows
}
